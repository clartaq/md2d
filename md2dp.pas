{
The MIT License

Copyright 2022 David D. Clark.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
}

program md2dp;

// A molecular dynamics  simulation of a Lennard-Jones fluid.

{$mode objfpc}{$H+}

uses {$IFDEF UNIX}
  Cthreads, {$ENDIF}
  Classes,
  Math,
  SysUtils;

const

  // Simulation parameters
  NUM_PARTICLES = 500;
  NPM1 = NUM_PARTICLES - 1;
  FORCE_CUTOFF = 3.0;
  FORCE_CUTOFF2 = FORCE_CUTOFF * FORCE_CUTOFF;
  // Can't create constants from such complicated expressions in Pascal.
  // P_EAT_CUTOFF          = 4 * (power(FORCE_CUTOFF, -12) - power(FORCE_CUTOFF, -6));
  TARGET_TEMPERATURE = 0.42;
  TIME_STEP = 0.020;
  DT_OVER_2 = 0.5 * TIME_STEP;
  DT_SQUARED_OVER_2 = 0.5 * TIME_STEP * TIME_STEP;
  BOX_WIDTH = 100.0;
  BOX_WIDTH_MINUS_HALF = BOX_WIDTH - 0.5;
  BOX_HEIGHT = 100.0;
  BOX_HEIGHT_MINUS_HALF = BOX_HEIGHT - 0.5;
  WALL_STIFFNESS = 50.0;
  EQUILIBRATION_TIME = 200.0;
  STEPS_BETWEEN_EQUIL_RESCALING = 10;
  PRODUCTION_TIME = 200.0;
  STEPS_BETWEEN_PROD_RESCALING = 100;
  STEPS_PER_PRINTOUT = 50;

type
  TMdArray = array[0..NPM1] of double; //(pred(NUM_PARTICLES))] of double;

var
  // The arrays of particle information.
  x, y, vx, vy, ax, ay: TMdArray;

  StepsAccomplished: integer = 0;
  SampleCount: integer = 0;

  PEatCutoff: double;

  Time: double = 0.0;
  Ke: double = 0.0;
  Pe: double = 0.0;
  CurrentTemp: double = 0.0;
  TotalTemp: double = 0.0;
  AverageTemp: double = 0.0;
  Pressure: double = 0.0;
  TotalPressure: double = 0.0;
  AveragePressure: double = 0.0;

  StartTime, StopTime: TDateTime;

  function Algo647Uniform: double; forward;
  function GaussianDeviateMarsaglia: double; forward;
  procedure ComputeAccelerations; forward;

  // Execute one time step using the Verlet algorithm.
  procedure SingleStep;
  var
    i: integer;
  begin
    for i := 0 to NPM1 do
    begin
      // Update positions.
      x[i] := x[i] + (vx[i] * TIME_STEP) + (ax[i] * DT_SQUARED_OVER_2);
      y[i] := y[i] + (vy[i] * TIME_STEP) + (ay[i] * DT_SQUARED_OVER_2);
      // Update velocities "halfway".
      vx[i] := vx[i] + ax[i] * DT_OVER_2;
      vy[i] := vy[i] + ay[i] * DT_OVER_2;
    end;
    // ComputeAccelerations
    for i := 0 to NPM1 do
    begin
      // Finish updating velocities using the accelerations.
      vx[i] := vx[i] + ax[i] * DT_OVER_2;
      vy[i] := vy[i] + ay[i] * DT_OVER_2;
    end;
  end;

  // Compute accelerations of the molecules from their current positions using
  // the Lennard-Jones potential.
  procedure ComputeAccelerations;
  begin

  end;

  // Reset accumulators and counters used for some measurements.
  procedure ResetMeasurements;
  begin
    TotalTemp := 0.0;
    TotalPressure := 0.0;
    SampleCount := 0;
  end;

  // Calculate and return the instantaneous kinetic energy.
  function CalculateKineticEnergy: double;
  var
    ek: double;
    i: integer;
  begin
    ek := 0.0;
    for i := 0 to NPM1 do
    begin
      ek := ek + 0.5 * (vx[i] * vx[i] + vy[i] * vy[i]);
    end;
    Result := ek;
  end;

  // Compute accumulated property values from sampled values.
  procedure ComputeProperties;
  begin
    SampleCount := SampleCount + 1;
    Ke := CalculateKineticEnergy;
    CurrentTemp := Ke / NUM_PARTICLES;
    TotalTemp := TotalTemp + CurrentTemp;
    AverageTemp := TotalTemp / SampleCount;
    TotalPressure := TotalPressure + Pressure;
    AveragePressure := TotalPressure / SampleCount;
  end;

  // Return a small random number uniformly distributed around 0
  // with a maximum magnitude of epsilon.
  function Nudge(Epsilon: double): double;
  begin
    Result := (Algo647Uniform - 0.5) * Epsilon;
  end;

  // Place the particles in empty locations of a grid on the simulation
  // starting from the lower left. Applies a small amount of "jitter"
  // to break up the regularity of the layout a little. Does not check for
  // overflow of the grid (too many particles to fit).
  procedure PlaceParticles;
  var
    Spacing, HalfSpacing, Bwmhs, Jitter, XPos, YPos: double;
    i: integer;
  begin
    Spacing := 1.3; // Minimum space between particle centers.
    HalfSpacing := Spacing / 2.0;
    Bwmhs := BOX_WIDTH - HalfSpacing; // Box width minus half spacing.
    Jitter := 0.01; // Random amount to break up regularity.
    XPos := HalfSpacing;
    YPos := BOX_HEIGHT - HalfSpacing;

    for i := 0 to NPM1 do
    begin
      x[i] := XPos + Nudge(Jitter);
      y[i] := YPos + Nudge(Jitter);
      XPos := XPos + Spacing;
      if XPos > Bwmhs then
      begin
        XPos := HalfSpacing;
        YPos := YPos - Spacing;
      end;
    end;
  end;

  // Re-scale the velocities to the target temperature.
  procedure RescaleVelocities;
  var
    ScalingFactor, VelocitySquaredSum: double;
    i: integer;
  begin
    VelocitySquaredSum := 0.0;
    for i := 0 to NPM1 do
    begin
      VelocitySquaredSum := VelocitySquaredSum + vx[i] * vx[i] + vy[i] * vy[i];
    end;
    ScalingFactor := 2.0 * NUM_PARTICLES * TARGET_TEMPERATURE / VelocitySquaredSum;
    ScalingFactor := sqrt(ScalingFactor);
    for i := 0 to NPM1 do
    begin
      vx[i] := vx[i] * ScalingFactor;
      vy[i] := vy[i] * ScalingFactor;
    end;
  end;

  // Remove any net momentum from the system of particles.
  procedure RemoveDrift;
  var
    VSumX, VSumY: double;
    i: integer;
  begin
    VSumX := 0.0;
    VSumY := 0.0;
    for i := 0 to NPM1 do
    begin
      VSumX := VSumX + vx[i];
      VSumY := VSumY + vy[i];
    end;

    for i := 0 to NPM1 do
    begin
      vx[i] := vx[i] - VSumX / NUM_PARTICLES;
      vy[i] := vy[i] - VSumY / NUM_PARTICLES;
    end;

  end;

  // Initialize particle positions, scale velocities to target temperature,
  // and remove systematic drift.
  procedure InitializeParticles;
  var
    i: integer;
  begin
    for i := 0 to NPM1 do //pred(NUM_PARTICLES) do
    begin
      x[i] := 0.0;
      y[i] := 0.0;
      vx[i] := GaussianDeviateMarsaglia - 0.5;
      vy[i] := GaussianDeviateMarsaglia - 0.5;
      ax[i] := 0.0;
      ay[i] := 0.0;
    end;

    PlaceParticles;
    RescaleVelocities;
    RemoveDrift;

  end;

  procedure PrintConfig;
  begin
    writeln;
    writeln(' Simulation Configuration:');
    writeln('   Number of particles        : ', NUM_PARTICLES);
    writeln('   Cutoff radius              : ', FORCE_CUTOFF: 0: 2);
    writeln('   Target temperature         : ', TARGET_TEMPERATURE: 0: 3);
    writeln('   Integration step size      : ', TIME_STEP: 0: 3);
    writeln('   Width of simulation area   : ', BOX_WIDTH: 0: 1);
    writeln('   Height of simulation area  : ', BOX_HEIGHT: 0: 1);
    writeln('   Area of simulation         : ', (BOX_WIDTH * BOX_HEIGHT): 0: 1);
    writeln('   Wall stiffness             : ', WALL_STIFFNESS: 0: 1);
    writeln('   Equilibration time         : ', EQUILIBRATION_TIME: 0: 1);
    writeln('      Steps between rescaling : ', STEPS_BETWEEN_EQUIL_RESCALING);
    writeln('   Production time            : ', PRODUCTION_TIME: 0: 1);
    writeln('      Steps between rescaling : ', STEPS_BETWEEN_PROD_RESCALING);
    writeln('   Steps per pringout         : ', STEPS_PER_PRINTOUT);
  end;

  procedure PrintPropertiesHeader;
  begin
    writeln('Time      Temp.   Pressure Tot. E.   Kin. E.   Pot. E.   Steps');
  end;

  procedure PrintProperties;
  begin
    writeln(Time: 7: 3, ',  ', AverageTemp: 5: 3, ',  ', AveragePressure: 6: 4,
      ',  ', (Ke + Pe): 7: 2, ',  ', Ke: 6: 2, ',  ', Pe: 8: 2, ',  ',
      StepsAccomplished);
  end;

  procedure PrintRandomNumbers;
  var
    i: integer;

  begin
    writeln;
    writeln('Random numbers from a uniform distribution.');
    writeln;
    for i := 0 to 9 do
      writeln('Uniform #', i, ': ', Algo647Uniform: 0: 17);

    writeln;
    writeln('Random numbers for a normal distribution.');
    writeln;
    for i := 0 to 9 do
      writeln('Normal #', i, ': ', GaussianDeviateMarsaglia: 0: 17);
  end;


  // Random number stuff.

  // The algo_647_uniform function returns a random uniform deviate in
  // the range [0, 1) calculated using the ACM algorithm 647 as
  // recorded at: http://calgo.acm.org.
  // The algorithm is from the paper "ALGORITHM 647: Implementation and
  // relative Efficiency of Quasirandom Sequence Generators" by Bennet L. Fox,
  // ACM Transactions on Mathematical Software (TOMS), Volume 12, No 4,
  // Dec 1986, pp. 362-376.

  // I have never seen the actual paper - it's behind a paywall - but the
  // FORTRAN version of the algorithm is easy enough to understand.

var
  RngSeed: integer = 12345;

  // Return a random unform deviate as a double.
  function Algo647Uniform: double;
  var
    k: int32;
    partial: int32;

  begin
    k := floor(RngSeed / 127773);
    partial := 16807 * (RngSeed - k * 127773) - k * 2836;
    if (partial < 0) then
      RngSeed := partial + high(int32)
    else
      RngSeed := partial;
    Result := RngSeed * 4.656612875e-10;
  end;

  // The following is an implementation of the Marsaglia polar
  // method for generating random standard normal deviates. See
  // https://en.wikipedia.org/wiki/Marsaglia_polar_method.

var
  has_spare: boolean = False;
  spare: double;

  // Return a random standard normal deviate as a double.
  function GaussianDeviateMarsaglia: double;
  var
    fac, rsq, r1, r2: double;
  begin
    if has_spare then
    begin
      has_spare := False;
      Result := spare;
    end
    else
    begin
      repeat
        r1 := 2.0 * Algo647Uniform - 1.0;
        r2 := 2.0 * Algo647Uniform - 1.0;
        rsq := r1 * r1 + r2 * r2;
      until ((rsq < 1) and not (rsq = 0.0));
      fac := sqrt(-2.0 * ln(rsq) / rsq);
      spare := r1 * fac;
      has_spare := True;
      Result := r2 * fac;
    end;
  end;

begin // main program

  // Do initialization that cannot be done in declarations.
  PEatCutoff := 4 * (power(FORCE_CUTOFF, -12) - power(FORCE_CUTOFF, -6));

  writeln('md2d - MD simulation of 2D argon gas Lennard-Jones system.');
  writeln('This version is written in Pascal and running version ');
  printConfig;

  PrintRandomNumbers;

  writeln;
  writeln('Equilibration');
  PrintPropertiesHeader;
  PrintProperties;
  StartTime := now;
  while Time < EQUILIBRATION_TIME do
  begin
    SingleStep;
    Inc(StepsAccomplished);
    Time := Time + TIME_STEP;
    if (StepsAccomplished mod STEPS_PER_PRINTOUT) = 0 then
    begin
      ComputeProperties;
      PrintProperties;
    end;
    if (StepsAccomplished mod STEPS_BETWEEN_EQUIL_RESCALING) = 0 then
    begin
      RescaleVelocities;
    end;
  end;

  writeln;
  writeln('Production');
  PrintPropertiesHeader;
  while Time < (EQUILIBRATION_TIME + PRODUCTION_TIME) do
  begin
    SingleStep;
    Inc(StepsAccomplished);
    Time := Time + TIME_STEP;
  end;

  // Very crude calculation of steps per second. Includes print time.
  StopTime := now;
  writeln(FormatDateTime('hh.nn.ss.zzz', StopTime - StartTime), ' time taken');
end.
