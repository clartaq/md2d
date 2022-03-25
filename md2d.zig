//
// The MIT License
//
// Copyright 2022 David D. Clark.
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//

//
// A molecular dynamics simulation of a two-dimensional Lennard-Jones fluid.
//
// Author: david
//

const std = @import("std");
const builtin = @import("builtin");
const math = std.math;
const cTime = @cImport(@cInclude("time.h"));

// Simulation parameters.
const num_particles: i32 = 500;
const force_cutoff: f64 = 3.0;
const force_cutoff_2: f64 = force_cutoff * force_cutoff;
const p_eat_cutoff: f64 = 4.0 * (math.pow(f64, force_cutoff, -12.0) - math.pow(f64, force_cutoff, -6.0));
const target_temperature: f64 = 0.42;
const time_step: f64 = 0.020;
const dt_over_2 = 0.5 * time_step;
const dt_squared_over_2 = 0.5 * time_step * time_step;
const box_width: f64 = 100.0;
const box_width_minus_half: f64 = box_width - 0.5;
const box_height: f64 = 100.0;
const box_height_minus_half: f64 = box_height - 0.5;
const wall_stiffness: f64 = 50.0;
const equilibration_time: f64 = 200.0;
const steps_between_equil_rescaling: i32 = 10;
const production_time: f64 = 200.0;
const steps_between_prod_rescaling: i32 = 100;
const steps_per_printout: i32 = 50;

// Other globals.
var x: [num_particles]f64 = std.mem.zeroes([num_particles]f64);
var y: [num_particles]f64 = std.mem.zeroes([num_particles]f64);
var vx: [num_particles]f64 = std.mem.zeroes([num_particles]f64);
var vy: [num_particles]f64 = std.mem.zeroes([num_particles]f64);
var ax: [num_particles]f64 = std.mem.zeroes([num_particles]f64);
var ay: [num_particles]f64 = std.mem.zeroes([num_particles]f64);

var t: f64 = undefined; // Simulation time
var steps_accomplished: i32 = 0;

var v_sum_x: f64 = 0.0;
var v_sum_y: f64 = 0.0;
var k_e: f64 = 0.0;
var p_e: f64 = 0.0;
var current_t: f64 = 0.0;
var total_t: f64 = 0.0;
var average_t: f64 = 0.0;
var pressure: f64 = 0.0;
var total_p: f64 = 0.0;
var average_p: f64 = 0.0;
var sample_count: i32 = 0;

pub fn main() void {
    print("md2d - MD simulation of a 2D argon gas Lennard-Jones system.\n", .{});
    print("This version is written in Zig and running version: {}\n", .{builtin.zig_version});

    printConfig();

    initializeParticles();
    print("\nEquilibration\n", .{});
    printPropertiesHeader();
    printProperties();
    var start: cTime.clock_t = cTime.clock();
    while (t < equilibration_time) {
        singleStep();
        if (steps_accomplished == 0) {
            printProperties();
        }
        steps_accomplished += 1;
        t += time_step;
        if (@mod(steps_accomplished, steps_per_printout) == 0) {
            computeProperties();
            printProperties();
        }
        if (@mod(steps_accomplished, steps_between_equil_rescaling) == 0) {
            rescaleVelocities();
        }
    }

    resetMeasurements();
    rescaleVelocities();
    print("\nProduction\n", .{});
    printPropertiesHeader();
    while (t < (equilibration_time + production_time)) {
        singleStep();
        steps_accomplished += 1;
        t += time_step;
        if (@mod(steps_accomplished, steps_per_printout) == 0) {
            computeProperties();
            printProperties();
        }
        if (@mod(steps_accomplished, steps_between_prod_rescaling) == 0) {
            rescaleVelocities();
        }
    }

    var finish: cTime.clock_t = cTime.clock();
    var elapsed: f64 = @intToFloat(f64, (finish - start)) / @intToFloat(f64, cTime.CLOCKS_PER_SEC);
    print("Elapsed time: {d:.3} seconds, {d:.3} steps per second.\n", .{ elapsed, @intToFloat(f64, steps_accomplished) / elapsed });
    std.os.exit(0);
}

/// Execute one time step using the Verlet algorithm.
fn singleStep() void {
    var i: usize = 0;
    while (i < num_particles) : (i += 1) {
        // Update positions.
        x[i] += (vx[i] * time_step) + (ax[i] * dt_squared_over_2);
        y[i] += (vy[i] * time_step) + (ay[i] * dt_squared_over_2);
        // "Half update" the velocities.
        vx[i] += ax[i] * dt_over_2;
        vy[i] += ay[i] * dt_over_2;
    }
    computeAccelerations();
    // Complete updates to velocities with new accelerations.
    i = 0;
    while (i < num_particles) : (i += 1) {
        vx[i] += ax[i] * dt_over_2;
        vy[i] += ay[i] * dt_over_2;
    }
}

// Compute accelerations of the molecules from their current positions using
// the Lennard-Jones potential.
fn computeAccelerations() void {
    var dx: f64 = undefined;
    var dy: f64 = undefined;
    var dx2: f64 = undefined;
    var dy2: f64 = undefined;
    var r_squared: f64 = undefined;
    var r_squared_inv: f64 = undefined;
    var attract: f64 = undefined;
    var repel: f64 = undefined;
    var f_over_r: f64 = undefined;
    var fx: f64 = undefined;
    var fy: f64 = undefined;
    var wall_force: f64 = 0.0;

    p_e = 0.0;

    // First, check for bouncing against the walls.
    var i: usize = 0;
    while (i < num_particles) : (i += 1) {
        if (x[i] < 0.5) {
            ax[i] = wall_stiffness * (0.5 - x[i]);
            wall_force += ax[i];
            p_e += 0.5 * wall_stiffness * (0.5 - x[i]) * (0.5 - x[i]);
        } else if (x[i] > box_width_minus_half) {
            ax[i] = wall_stiffness * (box_width_minus_half - x[i]);
            wall_force -= ax[i];
            p_e += 0.5 * wall_stiffness * (box_width_minus_half - x[i]) *
                (box_width_minus_half - x[i]);
        } else {
            ax[i] = 0.0;
        }

        if (y[i] < 0.5) {
            ay[i] = wall_stiffness * (0.5 - y[i]);
            wall_force += ay[i];
            p_e += 0.5 * wall_stiffness * (0.5 - y[i]) * (0.5 - y[i]);
        } else if (y[i] > box_height_minus_half) {
            ay[i] = wall_stiffness * (box_height_minus_half - y[i]);
            wall_force -= ay[i];
            p_e += 0.5 * wall_stiffness * (box_height_minus_half - y[i]) *
                (box_height_minus_half - y[i]);
        } else {
            ay[i] = 0.0;
        }
    }

    pressure = wall_force / (4.0 * box_width);

    // Next, compute interactions using the Lennard-Jones potential.
    i = 0;
    var j: usize = 0;
    while (i < num_particles) : (i += 1) {
        j = 0;
        while (j < i) : (j += 1) {
            dx = x[i] - x[j];
            dx2 = dx * dx;
            // Make sure the pair are close enough to bother.
            if (dx2 < force_cutoff_2) {
                dy = y[i] - y[j];
                dy2 = dy * dy;
                if (dy2 < force_cutoff_2) {
                    r_squared = dx2 + dy2;
                    if (r_squared < force_cutoff_2) {
                        r_squared_inv = 1.0 / r_squared;
                        attract = r_squared_inv * r_squared_inv * r_squared_inv;
                        repel = attract * attract;
                        p_e += (4.0 * (repel - attract)) - p_eat_cutoff;
                        f_over_r = 24.0 * ((2.0 * repel) - attract) * r_squared_inv;
                        fx = f_over_r * dx;
                        fy = f_over_r * dy;
                        // Add the force on to i's acceleration.
                        ax[i] += fx;
                        ay[i] += fy;
                        ax[j] -= fx; // Newton's 3rd law.
                        ay[j] -= fy;
                    }
                }
            }
        }
    }
}

/// Reset accumulators and counters used for some measurements.
fn resetMeasurements() void {
    total_t = 0.0;
    total_p = 0.0;
    sample_count = 0;
}

/// Calculate the instantaneous kinetic energy.
fn calculateKineticEnergy() f64 {
    var ek: f64 = 0.0;
    var i: usize = 0;
    while (i < num_particles) : (i += 1) {
        ek += 0.5 * (vx[i] * vx[i] + vy[i] * vy[i]);
    }
    return ek;
}

/// Compute accumulated property values from sampled values.
fn computeProperties() void {
    sample_count += 1;
    k_e = calculateKineticEnergy();
    current_t = k_e / @intToFloat(f64, num_particles);
    total_t += current_t;
    average_t = total_t / @intToFloat(f64, sample_count);
    total_p += pressure;
    average_p = total_p / @intToFloat(f64, sample_count);
}

/// Return a small random number uniformly distributed around 0
/// with a maximum magnitude of epsilon.
fn nudge(epsilon: f64) f64 {
    return (algo647Uniform() - 0.5) * epsilon;
}

/// Place the particles in empty locations of a grid on the simulation
/// box starting from the lower left. Applies a small amount of "jitter"
/// to break up the regularity of the layout a little. Does not check for
/// overflow of the grid (too many particle to fit).
fn placeParticles() void {
    const spacing: f64 = 1.3; // Minimum space between particle centers.
    const half_spacing: f64 = spacing / 2.0;
    const bwmhs: f64 = box_width - half_spacing;
    const jitter: f64 = 0.01; // Random offset to break up regularity.
    var x_pos: f64 = half_spacing;
    var y_pos: f64 = box_height - half_spacing;

    var idx: usize = 0;
    while (idx < num_particles) : (idx += 1) {
        x[idx] = x_pos + nudge(jitter);
        y[idx] = y_pos + nudge(jitter);
        x_pos += spacing;
        if (x_pos > bwmhs) {
            x_pos = half_spacing;
            y_pos -= spacing;
        }
    }
}

/// Rescale the velocities consistent with the target temperature.
fn rescaleVelocities() void {
    var scaling_factor: f64 = undefined;
    var velocity_squared_sum: f64 = 0.0;
    var i: usize = 0;
    while (i < num_particles) : (i += 1) {
        velocity_squared_sum += vx[i] * vx[i] + vy[i] * vy[i];
    }
    scaling_factor = 2.0 * @intToFloat(f64, num_particles) * target_temperature / velocity_squared_sum;
    scaling_factor = math.sqrt(scaling_factor);
    i = 0;
    while (i < num_particles) : (i += 1) {
        vx[i] *= scaling_factor;
        vy[i] *= scaling_factor;
    }
}

/// Remove any net momentum from the system of particles.
fn removeDrift() void {
    var i: usize = 0;
    while (i < num_particles) : (i += 1) {
        v_sum_x += vx[i];
        v_sum_y += vy[i];
    }
    i = 0;
    while (i < num_particles) : (i += 1) {
        vx[i] -= v_sum_x / @intToFloat(f64, num_particles);
        vy[i] -= v_sum_y / @intToFloat(f64, num_particles);
    }
}

/// Initialize particle positions, scale velocities to target temperature,
/// and remove systematic drift.
fn initializeParticles() void {
    var i: usize = 0;
    while (i < num_particles) : (i += 1) {
        x[i] = 0.0;
        y[i] = 0.0;
        vx[i] = gaussianDeviateMarsaglia() - 0.5;
        vy[i] = gaussianDeviateMarsaglia() - 0.5;
        ax[i] = 0.0;
        ay[i] = 0.0;
    }
    placeParticles();
    rescaleVelocities();
    removeDrift();
}

/// Print to stdout, unbuffered, and silently returning on failure.
/// This use of standard output is not really thread safe.
fn print(comptime fmt: []const u8, args: anytype) void {
    const stdout = std.io.getStdOut().writer();
    nosuspend stdout.print(fmt, args) catch return;
}

/// Print out the simulation configuration.
fn printConfig() void {
    print("\nSimulation Configuration:\n", .{});
    print("   Number of particles        : {d}\n", .{num_particles});
    print("   Cutoff radius              : {d}\n", .{force_cutoff});
    print("   Target temperature         : {d}\n", .{target_temperature});
    print("   Integration step size      : {d}\n", .{time_step});
    print("   Width of simulation area   : {d}\n", .{box_width});
    print("   Height of simulation area  : {d}\n", .{box_height});
    print("   Area of simulation         : {d}\n", .{box_width * box_height});
    print("   Wall stiffness             : {d}\n", .{wall_stiffness});
    print("   Equilibration time         : {d}\n", .{equilibration_time});
    print("      Steps between rescaling : {d}\n", .{steps_between_equil_rescaling});
    print("   Production time            : {d}\n", .{production_time});
    print("      Steps between rescaling : {d}\n", .{steps_between_prod_rescaling});
    print("   Steps per print out        : {d}\n", .{steps_per_printout});
}

fn printPropertiesHeader() void {
    print("Time      Temp.   Pressure Tot. E.   Kin. E.   Pot. E.   Steps\n", .{});
}

fn printProperties() void {
    print("{d:>7.3},  {d:>5.3},  {d:>6.4},  {d:>7.2},  {d:>6.2},  {d:>8.2},  {d}\n", .{ t, average_t, average_p, (k_e + p_e), k_e, p_e, steps_accomplished });
}

/// Print some pseudo random numbers. Just for debugging purposes.
fn printRandomNumbers() void {
    print("\nRandom numbers from a uniform distribution:\n", .{});
    var idx: i32 = 0;
    while (idx < 10) : (idx += 1) {
        var rnu: f64 = algo647Uniform();
        print("Uniform #{d} {d}\n", .{ idx, rnu });
    }
    print("\nRandom numbers from a normal distribution:\n", .{});
    idx = 0;
    while (idx < 10) : (idx += 1) {
        var rnn: f64 = gaussianDeviateMarsaglia();
        print("Normal #{d} {d}\n", .{ idx, rnn });
    }
}

// The algo647Uniform function returns a random uniform deviate in
// the range [0, 1) calculated using the ACM algorithm 647 as
// recorded at: http://calgo.acm.org.
// The algorithm is from the paper "ALGORITHM 647: Implementation and
// relative Efficiency of Quasirandom Sequence Generators" by Bennet L. Fox,
// ACM Transactions on Mathematical Software (TOMS), Volume 12, No 4,
// Dec 1986, pp. 362-376.
//
// I have never seen the actual paper - it's behind a paywall - but the
// FORTRAN version of the algorithm is easy enough to understand.

var glbl_seed: i32 = 12345;
const INT32_MAX: i32 = math.maxInt(i32);

/// Return a random uniform deviate in the range 0 <= x < 1.0 as an f64.
fn algo647Uniform() f64 {
    var k: i32 = @divFloor(glbl_seed, 127773);
    // var k: i32 = floor(glbl_seed / 127773);
    var partial: i32 = 16807 * (glbl_seed - k * 127773) - k * 2836;
    if (partial < 0) {
        glbl_seed = partial + INT32_MAX;
    } else {
        glbl_seed = partial;
    }

    return @intToFloat(f64, glbl_seed) * 4.656612875e-10;
}

// The following is an implementation of the Marsaglia polar
// method for generating random standard normal deviates. See
// https://en.wikipedia.org/wiki/Marsaglia_polar_method.
/// Return a random standard normal deviate as an f64.
fn gaussianDeviateMarsaglia() f64 {
    const has_spare = struct {
        var static: bool = false;
    };
    const spare = struct {
        var static: f64 = 0.0;
    };
    var fac: f64 = undefined;
    var rsq: f64 = undefined;
    var r1: f64 = undefined;
    var r2: f64 = undefined;

    if (has_spare.static) {
        has_spare.static = false;
        return spare.static;
    } else {
        while (true) {
            r1 = (2.0 * algo647Uniform()) - 1.0;
            r2 = (2.0 * algo647Uniform()) - 1.0;
            rsq = (r1 * r1) + (r2 * r2);
            if (!((rsq >= 1.0) or (rsq == 0.0))) break;
        }
        fac = math.sqrt(-2.0 * @log(rsq) / rsq);
        spare.static = r1 * fac;
        has_spare.static = true;
        return r2 * fac;
    }
}
