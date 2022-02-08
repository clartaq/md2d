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

// Simulation parameters.
const num_particles: i32 = 500;
const force_cutoff: f64 = 3.0;
const force_cutoff_2: f64 = force_cutoff * force_cutoff;
const p_eat_cutoff: f64 = 4.0 * (math.pow(f64, force_cutoff, -12.0) - math.pow(f64, force_cutoff, -6.0));
const target_temperature: f64 = 0.42;
const time_step: f64 = 0.020;
const dt_over_two = 0.5 * time_step;
const dt_squared_over_2 = 0.5 * time_step * time_step;
const box_width: f64 = 100.0;
const box_widthMinusHalf: f64 = box_width - 0.5;
const box_height: f64 = 100.0;
const box_heightMinusHalf: f64 = box_height - 0.5;
const wall_stiffness: f64 = 50.0;
const equilibration_time: f64 = 100.0;
const steps_between_equil_rescaling: i32 = 10;
const production_time: f64 = 300.0;
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

var measurement_step: i32 = 0;

// Time stuff.
var dt: f64 = undefined; // The time step
var dt_over_2: f64 = undefined;
var dt_squared_over_2: f64 = undefined;

// Energy stuff
var total_energy: f64 = undefined;
var kinetic_energy: f64 = undefined;
var potential_energy: f64 = undefined;
var temp_potential_energy: f64 = undefined;
var potential_energy_sum: f64 = undefined;
var potential_energy_squared_sum: f64 = undefined;
var wall_force: f64 = undefined;

// Properties
var temperature: f64 = undefined;
var temperature_sum: f64 = undefined;
var temperature_squared_sum: f64 = undefined;

var pressure: f64 = undefined;
var instant_pressure: f64 = undefined;
var pressure_sum: f64 = undefined;
var pressure_squared_sum: f64 = undefined;

pub fn main() void {
    var rescaling_counter: i32 = 0;
    print("md2d - MD simulation of a 2D argon gas Lennard-Jones system.\n", .{});
    print("This version is written in Zig version: {}\n", .{builtin.zig_version});

    printConfig();

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

    initialize();
    computeAccelerations();
    initializeMeasurements();
    measureProperties();
    print("\nEquilibration steps, ", .{});
    print("rescaling velocities every {d} steps...\n", .{steps_between_equil_rescaling});
    // var i: usize = 0;
    while (t < equilibration_time) {
        singleStep();
        measureProperties();
        if (@mod(step, steps_per_printout) == 0) {
            printProperties();
        }
        if (rescaling_counter >= steps_between_equil_rescaling) {
            rescaleVelocities();
            rescaling_counter = 0;
        } else {
            rescaling_counter += 1;
        }
    }
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
    var y_pos: f64 = half_spacing;

    var idx: usize = 0;
    while (idx < num_particles) : (idx += 1) {
        x[idx] = x_pos + nudge(jitter);
        y[idx] = y_pos + nudge(jitter);
        // Velocities and accelerations are assumed to be zero at this point.
        // Ignore them for now.
        x_pos += spacing;
        if (x_pos > bwmhs) {
            x_pos = half_spacing;
            y_pos += spacing;
        }
    }
}

/// Initialize the simulation by placing the particles in their initial
/// positions and giving them normally distributed velocities consistent
/// with the target temperature.
fn initialize() void {
    // Initialize some globals.
    t = 0.0;
    dt = step_size;
    dt_over_2 = dt * 0.5;
    dt_squared_over_2 = dt_over_2 * dt;
    steps_accomplished = 0;

    placeParticles();

    // Normally distribute initial velocities.
    var x_velocity_sum: f64 = 0.0;
    var y_velocity_sum: f64 = 0.0;
    var idx: usize = 0;
    while (idx < num_particles) : (idx += 1) {
        vx[idx] = gaussianDeviateMarsaglia();
        x_velocity_sum += vx[idx];
        vy[idx] = gaussianDeviateMarsaglia();
        y_velocity_sum += vy[idx];
    }

    // Zero-out any momentum.
    idx = 0;
    while (idx < num_particles) : (idx += 1) {
        vx[idx] -= x_velocity_sum / @intToFloat(f64, num_particles);
        vy[idx] -= y_velocity_sum / @intToFloat(f64, num_particles);
    }

    rescaleVelocities();
}

fn singleStep() void {
    var i: usize = 0;
    while (i < num_particles) : (i += 1) {
        // Update positions.
        x[i] += (vx[i] * dt) + (ax[i] * dt_squared_over_2);
        y[i] += (vy[i] * dt) + (ay[i] * dt_squared_over_2);
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
    t += dt;
    steps_accomplished += 1;
}

fn checkWalls() void {
    wall_force = 0.0;

    var i: usize = 0;
    while (i < num_particles) : (i += 1) {

        // First compute the x component.
        if (x[i] < 0.5) {
            ax[i] = wall_stiffness * (0.5 - x[i]);
            wall_force += ax[i];
            temp_potential_energy += 0.5 * wall_stiffness * (0.5 - x[i]) * (0.5 - x[i]);
        } else if (x[i] > box_widthMinusHalf) {
            ax[i] = wall_stiffness * (box_widthMinusHalf - x[i]);
            wall_force -= -ax[i];
            temp_potential_energy += wall_stiffness *
                (box_widthMinusHalf - x[i]) * (box_widthMinusHalf - x[i]);
        } else {
            ax[i] = 0.0;
        }

        // Now the y component.
        if (y[i] < 0.5) {
            ay[i] = wall_stiffness * (0.5 - y[i]);
            wall_force += ay[i];
            temp_potential_energy += 0.5 * wall_stiffness * (0.5 - y[i]) * (0.5 - y[i]);
        } else if (y[i] > box_heightMinusHalf) {
            ay[i] = wall_stiffness * (box_heightMinusHalf - y[i]);
            wall_force -= ay[i];
            temp_potential_energy += 0.5 * wall_stiffness * (box_heightMinusHalf - y[i]) * (box_heightMinusHalf - y[i]);
        } else {
            ay[i] = 0.0;
        }

        instant_pressure = wall_force / ((2.0 * box_width) + (2.0 * box_height));
    }
    //    print("temp_potential_energy: {d}\n", .{temp_potential_energy});
    //    print("instant_pressure: {d}\n", .{instant_pressure});
}

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

    temp_potential_energy = 0.0;
    checkWalls();

    // Compute forces of interactions from Lennard-Jones potential.
    var i: usize = 0;
    while (i < num_particles) : (i += 1) {
        var j: usize = 0;
        while (j < i) : (j += 1) {
            dx = x[i] - y[i];
            dx2 = dx * dx;
            if (dx2 < force_cutoff_2) {
                dy = y[i] - y[j];
                dy2 = dy * dy;
                if (dy2 < force_cutoff_2) {
                    r_squared = dx2 + dy2;
                    if (r_squared < force_cutoff_2) {
                        r_squared_inv = 1.0 / r_squared;
                        attract = r_squared_inv * r_squared_inv * r_squared_inv;
                        repel = attract * attract;
                        temp_potential_energy += (4.0 * (repel - attract)) - p_eat_cutoff;
                        f_over_r = 24.0 * ((2.0 * repel) - attract) * r_squared_inv;
                        fx = f_over_r * dx;
                        fy = f_over_r * dy;
                        ax[i] += fx;
                        ay[i] += fy;
                        ax[j] -= fx;
                        ay[j] -= fy;
                    }
                }
            }
        }
    }
}

/// Rescale the velocities consistent with the target temperature.
fn rescaleVelocities() void {
    var scaling_factor: f64 = undefined;
    var velocity_squared_sum: f64 = 0.0;
    var idx: usize = 0;
    while (idx < num_particles) : (idx += 1) {
        velocity_squared_sum += vx[idx] * vx[idx] + vy[idx] * vy[idx];
    }
    scaling_factor = 2.0 * @intToFloat(f64, num_particles) * target_temperature / velocity_squared_sum;
    scaling_factor = math.sqrt(scaling_factor);
    idx = 0;
    while (idx < num_particles) : (idx += 1) {
        vx[idx] *= scaling_factor;
        vy[idx] *= scaling_factor;
    }
}

fn initializeMeasurements() void {
    measurement_step = 0;
    temperature_sum = 0.0;
    temperature_squared_sum = 0.0;
    pressure_sum = 0.0;
    pressure_squared_sum = 0.0;
    potential_energy_sum = 0.0;
    potential_energy_squared_sum = 0.0;
}

fn measureProperties() void {
    kinetic_energy = 0.0;
    var i: usize = 0;
    while (i < num_particles) : (i += 1) {
        kinetic_energy += 0.5 * (vx[i] * vx[i] + vy[i] * vy[i]);
    }
    total_energy = kinetic_energy + potential_energy;
    temperature = kinetic_energy / @intToFloat(f64, num_particles);

    measurement_step += 1;
    temperature_sum += temperature;
    temperature_squared_sum += temperature * temperature;
    potential_energy_sum += potential_energy;
    potential_energy_squared_sum += potential_energy * potential_energy;
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
    print("   Integration step size      : {d}\n", .{step_size});
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

fn printProperties() void {
    print("t: {d:>6.3}, KE: {e:>12.4}, TE: {e:>12.4}, stp: {d}\n", .{ t, kinetic_energy, total_energy, steps_accomplished });
    //   print("t: {d}, TE: {d}, KE: {d}\n", .{ t, total_energy, kinetic_energy });
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
    var partial: i32 = 16807 * (glbl_seed - k * 127773) - k * 2386;
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
