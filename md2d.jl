#
# The MIT License
#
# Copyright 2022 David D. Clark.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#

#
# A molecular dymamics simulation of a two-dimensional Lennard-Jones fluid.
#
# Author: david
#

using Printf
using StaticArrays

# Simulation parameters
const num_particles = 500
const force_cutoff = 3.0
const force_cutoff_2 = force_cutoff * force_cutoff
const p_eat_cutoff = 4.0 * ((force_cutoff^-12) - (force_cutoff^-6))
const target_temperature = 0.42
const time_step = 0.020
const dt_over_2 = 0.5 * time_step
const dt_squared_over_2 = 0.5 * time_step * time_step
const box_width = 100.0 # In units of molecular diameter.
const box_width_minus_half = box_width - 0.5
const box_height = 100.0
const box_height_minus_half = box_height - 0.5
const wall_stiffness = 50.0 # "Spring constant" for walls.
const equilibration_time = 200.0
const steps_between_equil_rescaling = 10
const production_time = 200.0
const steps_between_prod_rescaling = 100
const steps_per_printout = 50.0

const x = @MVector zeros(num_particles)
const y = @MVector zeros(num_particles)
const vx = @MVector zeros(num_particles)
const vy = @MVector zeros(num_particles)
const ax = @MVector zeros(num_particles)
const ay = @MVector zeros(num_particles)

const t = Ref{Float64}(0.0)
const steps_accomplished = Ref{Int64}(0)

const v_sum_x = Ref{Float64}(0.0)
const v_sum_y = Ref{Float64}(0.0)
const k_e = Ref{Float64}(0.0)
const p_e = Ref{Float64}(0.0)
const current_t = Ref{Float64}(0.0)
const total_t = Ref{Float64}(0.0)
const average_t = Ref{Float64}(0.0)
const pressure = Ref{Float64}(0.0)
const total_p = Ref{Float64}(0.0)
const average_p = Ref{Float64}(0.0)
const sample_count = Ref{Int64}(0)

function main()
    @printf "md2d - MD simulation of a 2D argon gas Lennard-Jones system.\n"
    @printf "This version is written in Julia and running version %s.\n" VERSION
    print_config()

    initialize_particles()

    @printf "\nEquilibration\n"
    print_properties_header()
    print_properties()
    start_ns = time_ns()
    while (t[] < equilibration_time)
        single_step()
        steps_accomplished[] += 1
        t[] += time_step
        if ((steps_accomplished[] % steps_per_printout) == 0)
            compute_properties()
            print_properties()
        end
        if ((steps_accomplished[] % steps_between_equil_rescaling) == 0)
            rescale_velocities()
        end
    end

    reset_measurements()
    rescale_velocities()
    @printf "\nProduction\n"
    print_properties_header()
    while (t[] < (equilibration_time + production_time))
        single_step()
        steps_accomplished[] += 1
        t[] += time_step
        if ((steps_accomplished[] % steps_per_printout) == 0)
            compute_properties()
            print_properties()
        end
        if ((steps_accomplished[] % steps_between_prod_rescaling) == 0)
            rescale_velocities()
        end
    end

    # Very crude calculation of steps per second. Includes print time.
    finish_ns = time_ns()
    elapsed_sec = (finish_ns - start_ns) / 1000000000
    @printf "Elapsed time: %.3f seconds, %.3f steps per second\n" elapsed_sec (
        steps_accomplished[] / elapsed_sec
    )
end

"Execute one time step using the Verlet algorithm."
function single_step()
    for i in eachindex(x)
        # Update positions.
        x[i] += (vx[i] * time_step) + (ax[i] * dt_squared_over_2)
        y[i] += (vy[i] * time_step) + (ay[i] * dt_squared_over_2)
        # Update velocities "halfway".
        vx[i] += (ax[i] * dt_over_2)
        vy[i] += (ay[i] * dt_over_2)
    end
    compute_accelerations()
    for i in eachindex(x)
        # Finish updating velocities using the new accelerations.
        vx[i] += (ax[i] * dt_over_2)
        vy[i] += (ay[i] * dt_over_2)
    end
end

"Compute accelerations of the molecules from their current positions using
 the Lennard-Jones potential."
function compute_accelerations()
    p_e[] = 0.0
    wall_force = 0.0

    # First, check for bouncing against the walls.
    for i in eachindex(x)
        if (x[i] < 0.5)
            ax[i] = wall_stiffness * (0.5 - x[i])
            wall_force += ax[i]
            p_e[] += 0.5 * wall_stiffness * (0.5 - x[i]) * (0.5 - x[i])
        elseif (x[i] > box_width_minus_half)
            ax[i] = wall_stiffness * (box_width_minus_half - x[i])
            wall_force -= ax[i]
            p_e[] +=
                0.5 *
                wall_stiffness *
                (box_width_minus_half - x[i]) *
                (box_width_minus_half - x[i])
        else
            ax[i] = 0.0
        end

        if (y[i] < 0.5)
            ay[i] = wall_stiffness * (0.5 - y[i])
            wall_force += ay[i]
            p_e[] = 0.5 * wall_stiffness * (0.5 - y[i]) * (0.5 - y[i])
        elseif (y[i] > box_height_minus_half)
            ay[i] = wall_stiffness * (box_height_minus_half - y[i])
            wall_force -= ay[i]
            p_e[] +=
                0.5 *
                wall_stiffness *
                (box_height_minus_half - y[i]) *
                (box_height_minus_half - y[i])
        else
            ay[i] = 0.0
        end
    end

    pressure[] = wall_force / (4.0 * box_width)

    # Next, compute interactions using the Lennard-Jones potential.
    for i in eachindex(x)
        for j = firstindex(x):(i-1)
            dx = x[i] - x[j]
            dx2 = dx * dx
            # Make sure the pair are close enough to bother.
            if (dx2 < force_cutoff_2)
                dy = y[i] - y[j]
                dy2 = dy * dy
                if (dy2 < force_cutoff_2)
                    r_squared = dx2 + dy2
                    if (r_squared < force_cutoff_2)
                        r_squared_inv = 1.0 / r_squared
                        attract = r_squared_inv * r_squared_inv * r_squared_inv
                        repel = attract * attract
                        p_e[] += (4.0 * (repel - attract)) - p_eat_cutoff
                        f_over_r = 24.0 * ((2.0 * repel) - attract) * r_squared_inv
                        fx = f_over_r * dx
                        fy = f_over_r * dy
                        # Add the force on to i's acceleration.
                        ax[i] += fx
                        ay[i] += fy
                        ax[j] -= fx # Newton's 3rd law
                        ay[j] -= fy
                    end
                end
            end
        end
    end
end

"Reset accumulators and counters used for some measurements."
function reset_measurements()
    total_t[] = 0.0
    total_p[] = 0.0
    sample_count[] = 0
end

"Return the instantaneous kinetic energy."
function calculate_kinetic_energy()
    ek = 0.0
    for i in eachindex(vx)
        ek += 0.5 * (vx[i] * vx[i] + vy[i] * vy[i])
    end
    return ek
end

"Compute accumulated property values from sampled values."
function compute_properties()
    sample_count[] += 1
    k_e[] = calculate_kinetic_energy()
    current_t[] = k_e[] / num_particles
    total_t[] += current_t[]
    average_t[] = total_t[] / sample_count[]
    total_p[] += pressure[]
    average_p[] = total_p[] / sample_count[]
end

const nudge_count = Ref{Int64}(0)

"Return a small random number uniformly distributed around 0
 with a maximum magnitude of epsilon."
function nudge(epsilon)
    return (algo_647_uniform() - 0.5) * epsilon
end

"Place the particles in empty locactions of a grid on the simulation
 box starting from the lower left. Applies a small amount of \"jitter\"
 to break up the regularity of the layout a little. Does not check for
 overflow of the grid (too many particles to fit)."
function place_particles()
    spacing = 1.3 # Minimum space between particle centers.
    half_spacing = spacing / 2.0
    bwmhs = box_width - half_spacing
    jitter = 0.01 # Random amount to break up regularity.
    x_pos = half_spacing
    y_pos = box_height - half_spacing

    for i in eachindex(x)
        x[i] = x_pos + nudge(jitter)
        y[i] = y_pos + nudge(jitter)
        x_pos += spacing
        if (x_pos > bwmhs)
            x_pos = half_spacing
            y_pos -= spacing
        end
    end
end

"Re-scale the velocities to the target temperature."
function rescale_velocities()
    velocity_squared_sum = 0.0
    for i in eachindex(vx)
        velocity_squared_sum += vx[i] * vx[i] + vy[i] * vy[i]
    end
    scaling_factor = 2.0 * num_particles * target_temperature / velocity_squared_sum
    scaling_factor = sqrt(scaling_factor)
    for i in eachindex(vx)
        vx[i] *= scaling_factor
        vy[i] *= scaling_factor
    end
end

"Remove any net momentum from the system of particles."
function remove_drift()
    for i in eachindex(x)
        v_sum_x[] += vx[i]
        v_sum_y[] += vy[i]
    end
    for i in eachindex(x)
        vx[i] -= v_sum_x[] / num_particles
        vy[i] -= v_sum_y[] / num_particles
    end
end

"Initialize particle positions, scale velocities to target temperature,
 and remove systematic drift."
function initialize_particles()
    for i in eachindex(x)
        x[i] = 0.0
        y[i] = 0.0
        vx[i] = gaussian_deviate_marsaglia() - 0.5
        vy[i] = gaussian_deviate_marsaglia() - 0.5
        ax[i] = 0.0
        ay[i] = 0.0
    end
    place_particles()
    rescale_velocities()
    remove_drift()
end

#
# Printing stuff.
#

function print_config()
    @printf "\n Simulation Configuration:\n"
    @printf "   Number of particles        : %d\n" num_particles
    @printf "   Cutoff radius              : %.2f\n" force_cutoff
    @printf "   Target temperature         : %.3f\n" target_temperature
    @printf "   Integration step size      : %.3f\n" time_step
    @printf "   Width of simulation area   : %.1f\n" box_width
    @printf "   Height of simulation area  : %.1f\n" box_height
    @printf "   Area of simulation         : %.1f\n" box_width * box_height
    @printf "   Wall stiffness             : %.1f\n" wall_stiffness
    @printf "   Equilibration time         : %.1f\n" equilibration_time
    @printf "      Steps between rescaling : %d\n" steps_between_equil_rescaling
    @printf "   Production time            : %.1f\n" production_time
    @printf "      Steps between rescaling : %d\n" steps_between_prod_rescaling
    @printf "   Steps per print out        : %d\n" steps_per_printout
end

function print_properties_header()
    @printf "Time      Temp.   Pressure Tot. E.   Kin. E.   Pot. E.   Steps\n"
end

function print_properties()
    total_e = k_e[] + p_e[]
    @printf "%7.3f,  %5.3f,  %6.4f,  %7.2f,  %6.2f,  %8.2f,  %d\n" t[] average_t[] average_p[] total_e k_e[] p_e[] steps_accomplished[]
end

"Print some pseudo random numbers. Just for debugging purposes."
function print_random_numbers()
    @printf "\nRandom numbers from a uniform distribution.\n"
    for i = 0:9
        @printf "Uniform #%i: %.17f\n" i algo_647_uniform()
    end
    @printf "\nRandom numbers from a normal distribution.\n"
    for i = 0:9
        @printf "Normal #%i: %.17f\n" i gaussian_deviate_marsaglia()
    end
end

#
# Random number stuff.
#

# The algo_647_uniform function returns a random uniform deviate in
# the range [0, 1) calculated using the ACM algorithm 647 as
# recorded at: http://calgo.acm.org.
# The algorithm is from the paper "ALGORITHM 647: Implementation and
# relative Efficiency of Quasirandom Sequence Generators" by Bennet L. Fox,
# ACM Transactions on Mathematical Software (TOMS), Volume 12, No 4,
# Dec 1986, pp. 362-376.
#
# I have never seen the actual paper - it's behind a paywall - but the
# FORTRAN version of the algorithm is easy enough to understand.

# Something to hold an mutable seed value for the random number generator.
mutable struct seedType
    seed::Int32
end

# Something to hold an mutable seed value for the random number generator.
glbl_seed = seedType(12345)

"Allow the seed to be set externally."
function algo_647_set_seed(new_seed)
    setproperty!(glbl_seed, :seed, new_seed)
end

"Return a random uniform deviate in the range 0 <= x < 1.0 as a Float64."
function algo_647_uniform()
    k = floor(glbl_seed.seed / 127773)
    partial = 16807 * (glbl_seed.seed - k * 127773) - k * 2836
    if (partial < 0)
        setproperty!(glbl_seed, :seed, partial + typemax(Int32))
    else
        setproperty!(glbl_seed, :seed, partial)
    end
    return glbl_seed.seed * 4.656612875e-10
end

# The following is an implementation of the Marsaglia polar
# method for generating random standard normal deviates. See
# https://en.wikipedia.org/wiki/Marsaglia_polar_method.

mutable struct spareType
    has_spare::Bool
    spare::Float64
end

glbl_spare = spareType(false, 0.0)

"Return a random standard normal deviate as a double."
function gaussian_deviate_marsaglia()
    fac = 0.0
    rsq = 0.0
    r1 = 0.0
    r2 = 0.0
    if (glbl_spare.has_spare)
        setproperty!(glbl_spare, :has_spare, false)
        return glbl_spare.spare
    else
        while true
            r1 = 2.0 * algo_647_uniform() - 1.0
            r2 = 2.0 * algo_647_uniform() - 1.0
            rsq = r1 * r1 + r2 * r2
            (rsq >= 1.0 || rsq == 0.0) || break
        end
        fac = sqrt(-2.0 * log(rsq) / rsq)
        setproperty!(glbl_spare, :spare, r1 * fac)
        setproperty!(glbl_spare, :has_spare, true)
        return r2 * fac
    end
end

main()
