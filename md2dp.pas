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

uses
{$IFDEF UNIX}
cthreads,
{$ENDIF}
Classes, Math, SysUtils;

const

// Simulation parameters
   NUM_PARTICLES         = 500;
   FORCE_CUTOFF          = 3.0;
   FORCE_CUTOFF2         = FORCE_CUTOFF * FORCE_CUTOFF;
   // Can't create constants from such complicated expressions in Pascal.
   // P_EAT_CUTOFF          = 4 * (power(FORCE_CUTOFF, -12) - power(FORCE_CUTOFF, -6));
   TARGET_TEMPERATURE    = 0.42;
   TIME_STEP             = 0.020;
   DT_OVER_2             = 0.5 * TIME_STEP;
   DT_SQUARED_OVER_2     = 0.5 * TIME_STEP * TIME_STEP;
   BOX_WIDTH             = 100.0;
   BOX_WIDTH_MINUS_HALF  = BOX_WIDTH - 0.5;
   BOX_HEIGHT            = 100.0;
   BOX_HEIGHT_MINUS_HALF = BOX_HEIGHT - 0.5;
   WALL_STIFFNESS        = 50.0;
   EQUILIBRATION_TIME    = 200.0;
   STEPS_BETWEEN_EQUIL_RESCALING = 10;
   PRODUCTION_TIME       = 200.0;
   STEPS_BETWEEN_PROD_RESCALING = 100;
   STEPS_PER_PRINTOUT    = 50;

var
   start, stop  : TDateTime;
   P_EAT_CUTOFF : Double;

function algo_647_uniform: Double; forward;
   function gaussian_deviate_marsaglia: Double; forward;

      procedure printConfig;
      begin
         writeln('printConfig');
         writeln;
         writeln(' Simulation Configuration:');
         writeln('   Number of particles        : ', NUM_PARTICLES);
         writeln('   Cutoff radius              : ', FORCE_CUTOFF:0:2);
         writeln('   Target temperature         : ', TARGET_TEMPERATURE:0:3);
         writeln('   Integration step size      : ', TIME_STEP:0:3);
         writeln('   Width of simulation area   : ', BOX_WIDTH:0:1);
         writeln('   Height of simulation area  : ', BOX_HEIGHT:0:1);
         writeln('   Area of simulation         : ', (BOX_WIDTH * BOX_HEIGHT):0:1);
         writeln('   Wall stiffness             : ', WALL_STIFFNESS:0:1);
         writeln('   Equilibration time         : ', EQUILIBRATION_TIME:0:1);
         writeln('      Steps between rescaling : ', STEPS_BETWEEN_EQUIL_RESCALING);
         writeln('   Production time            : ', PRODUCTION_TIME:0:1);
         writeln('      Steps between rescaling : ', STEPS_BETWEEN_PROD_RESCALING);
         writeln('   Steps per pringout         : ', STEPS_PER_PRINTOUT);
      end;

      procedure printPropertiesHeader;
      begin
         writeln('Time      Temp.   Pressure Tot. E.   Kin. E.   Pot. E.   Steps');
      end;

      procedure printRandomNumbers;
      var
         i : Integer;

      begin
         writeln;
         writeln('Random numbers from a uniform distribution.');
         writeln;
         for i := 0 to 9 do
            writeln('Uniform #', i, ': ', algo_647_uniform:0:17);

         writeln;
         writeln('Random numbers for a normal distribution.');
         writeln;
         for i := 0 to 9 do
            writeln('Normal #', i, ': ', gaussian_deviate_marsaglia:0:17);
      end;


   // Random number stuff.
   //
   // The algo_647_uniform function returns a random uniform deviate in
   // the range [0, 1) calculated using the ACM algorithm 647 as
   // recorded at: http://calgo.acm.org.
   // The algorithm is from the paper "ALGORITHM 647: Implementation and
   // relative Efficiency of Quasirandom Sequence Generators" by Bennet L. Fox,
   // ACM Transactions on Mathematical Software (TOMS), Volume 12, No 4,
   // Dec 1986, pp. 362-376.
   //
   // I have never seen the actual paper - it's behind a paywall - but the
   // FORTRAN version of the algorithm is easy enough to understand.

   var
      rng_seed: Integer = 12345;

   // Return a random unform deviate as a double.
      function algo_647_uniform: Double;
      var
         k       : Int32;
         partial : Int32;

      begin
         k := floor(rng_seed/127773);
         partial := 16807*(rng_seed - k*127773) - k*2836;
         if (partial < 0) then
            rng_seed := partial + high(Int32)
         else
            rng_seed := partial;
         algo_647_uniform := rng_seed*4.656612875e-10;
      end;

   // The following is an implementation of the Marsaglia polar
   // method for generating random standard normal deviates. See
   // https://en.wikipedia.org/wiki/Marsaglia_polar_method.

   var
      has_spare : Boolean = false;
      spare     : Double;

   // Return a random standard normal deviate as a double.
      function gaussian_deviate_marsaglia: Double;
      var
         fac, rsq, r1, r2 : Double;

      begin
         if has_spare then begin
            has_spare := false;
            gaussian_deviate_marsaglia := spare
         end
      else begin
         repeat
            r1 := 2.0*algo_647_uniform - 1.0;
            r2 := 2.0*algo_647_uniform - 1.0;
            rsq := r1*r1 + r2*r2;
         until ( (rsq < 1) and not( rsq = 0.0));
         fac := sqrt(-2.0*ln(rsq)/rsq);
         spare := r1*fac;
         has_spare := true;
         gaussian_deviate_marsaglia := r2*fac;
      end
      end;

      procedure maxvals;
      const width = 20;
      begin
         writeln;
         writeln( 'these change depending on compiler mode:' );
         writeln( '----------------------------------------' );
         writeln( 'maxint:           :', maxint : width );
         writeln( 'high( integer )   :', high( integer ) : width );
         writeln;
         writeln( 'constant, regardless of mode or target: ' );
         writeln( '----------------------------------------' );
         writeln( 'high( int32 )     :', high( int32 ) : width );
         writeln( 'high( int64 )     :', high( int64 ) : width );
         writeln;
         writeln( 'variable, depending on target cpu:' );
         writeln( '----------------------------------------' );
         writeln( 'sizeof( pointer ) :', sizeof( pointer ) : width );
         writeln;
         writeln( 'compile-time definitions:' );
         writeln( '----------------------------------------' );
         {$IFDEF cpu64} writeln( 'cpu64' ); {$ENDIF}
         {$IFDEF cpu32} writeln( 'cpu32' ); {$ENDIF}
         {$IFDEF cpu16} writeln( 'cpu16' ); {$ENDIF}
         writeln;
      end;

   begin // main program
      start := now;

      P_EAT_CUTOFF := 4 * (power(FORCE_CUTOFF, -12) - power(FORCE_CUTOFF, -6));

      writeln('md2d - MD simulation of 2D argon gas Lennard-Jones system.');
      writeln('This version is written in Pascal and running version ');
      printConfig;

      printRandomNumbers;

      writeln;
      writeln('Equilibration');
      writeln;
      printPropertiesHeader;

      writeln;
      writeln('Production');
      writeln;
      printPropertiesHeader;

      //Very crude calculation of steps per second. Includes print time.

      stop := now;
      writeln(FormatDateTime('hh.nn.ss.zzz', stop - start), ' time taken');
   end.
