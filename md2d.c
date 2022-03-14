/*
 * The MIT License
 *
 * Copyright 2022 David D. Clark.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

/*
 * A molecular dymamics simulation of a two-dimensional Lennard-Jones fluid.
 *
 * Author: david
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

// Simulation parameters
const int NUM_PARTICLES = 500;
const double FORCE_CUTOFF = 3.0;
const double FORCE_CUTOFF_2 = FORCE_CUTOFF * FORCE_CUTOFF;
double P_EAT_CUTOFF; // = 4 * (pow(FORCE_CUTOFF, -12) - pow(FORCE_CUTOFF, -6));
const double TARGET_TEMPERATURE = 0.42;
const double TIME_STEP = 0.020;
const double DT_OVER_2 = 0.5 * TIME_STEP;
const double DT_SQUARED_OVER_2 = 0.5 * TIME_STEP * TIME_STEP;
const double BOX_WIDTH = 100.0;      // In units of the molecular diameter.
const double BOX_WIDTH_MINUS_HALF = BOX_WIDTH - 0.5;
const double BOX_HEIGHT = 100.0;
const double BOX_HEIGHT_MINUS_HALF = BOX_HEIGHT - 0.5;
const double WALL_STIFFNESS = 50.0;  // "Spring constant" for walls.
const double EQUILIBRATION_TIME = 100.0;
const int STEPS_BETWEEN_EQUIL_RESCALING = 10;
const double PRODUCTION_TIME = 300.0;
const int STEPS_PER_PRINTOUT = 50;

double x[NUM_PARTICLES];
double y[NUM_PARTICLES];
double vx[NUM_PARTICLES];
double vy[NUM_PARTICLES];
double ax[NUM_PARTICLES];
double ay[NUM_PARTICLES];

double t = 0.0;
int steps_accomplished = 0;

double v_sum_x = 0.0;
double v_sum_y = 0.0;
double kE = 0.0;
double pE = 0.0;
double current_t, total_t, average_t = 0.0;
double pressure, total_p, average_p = 0.0;
int sample_count = 0;

static int rng_seed = 12345;

void print_config();
void initialize_particles();
void print_properties_header();
void print_properties();
void single_step();
void compute_accelerations();
void compute_properties();
void rescale_velocities();
void reset_measurements();
double calculate_kinetic_energy();
void print_random_numbers();
void initialize_particles();
void initialize_measurements();
void measure_properties();
double algo_647_uniform();
double gaussian_deviate_marsaglia();

int main() {
  P_EAT_CUTOFF = 4 * (pow(FORCE_CUTOFF, -12) - pow(FORCE_CUTOFF, -6));

  printf("md2d - MD simulation of a 2D argon gas Lennard-Jones system.\n");
  printf("This version is written in C and comiled for Standard %ld\n", __STDC_VERSION__);
  print_config();

  initialize_particles();

  printf("\nEqulibration\n");
  print_properties_header();
  print_properties();
  clock_t start = clock();
  while (t < EQUILIBRATION_TIME) {
    single_step();
    steps_accomplished++;
    t += TIME_STEP;
    if ((steps_accomplished % STEPS_PER_PRINTOUT) == 0) {
      compute_properties();
      print_properties();
    }
    if ((steps_accomplished % STEPS_BETWEEN_EQUIL_RESCALING) == 0) {
      rescale_velocities();
    }
  }

  reset_measurements();
  rescale_velocities();
  printf("\nProduction\n");
  print_properties_header();
  while (t < (EQUILIBRATION_TIME + PRODUCTION_TIME)) {
    single_step();
    steps_accomplished++;
    t += TIME_STEP;
    if ((steps_accomplished % STEPS_PER_PRINTOUT) == 0) {
      compute_properties();
      print_properties();
    }
  }

  // Very crude calculation of steps per second. Includes print time.
  clock_t finish = clock();
  double elapsed = (double) (finish - start)/CLOCKS_PER_SEC;
  printf("Elapsed time: %.3f seconds, %.3f steps per second.\n", elapsed,
         steps_accomplished / elapsed);
  exit(0);
}

/*
 * Execute one time step using the Verlet algorithm.
 */
void single_step() {
  for (int i = 0; i < NUM_PARTICLES; i++) {
    // Update positions.
    x[i] += (vx[i] * TIME_STEP) + (ax[i] * DT_SQUARED_OVER_2);
    y[i] += (vy[i] * TIME_STEP) + (ay[i] * DT_SQUARED_OVER_2);
    // Update velocities "halfway".
    vx[i] += (ax[i] * DT_OVER_2);
    vy[i] += (ay[i] * DT_OVER_2);
  }
  compute_accelerations();
  for (int i = 0; i < NUM_PARTICLES; i++) {
    // Finish updating velocities using the new accellerations.
    vx[i] += (ax[i] * DT_OVER_2);
    vy[i] += (ay[i] * DT_OVER_2);
  }
}

/*
 * Compute accelerations of the molecules from their current positions using
 * the Lennard-Jones potential.
 */
void compute_accelerations() {

  int i, j;
  double dx, dy;
  double dx2, dy2, rSquared, rSquaredInv, attract, repel, fOverR, fx, fy;
  double wall_force = 0.0;

  pE = 0.0;

  // First, check for bouncing against the walls.
  for (i = 0; i < NUM_PARTICLES; i++) {
    if (x[i] < 0.5) {
      ax[i] = WALL_STIFFNESS * (0.5 - x[i]);
      wall_force += ax[i];
      pE += 0.5 * WALL_STIFFNESS * (0.5 - x[i]) * (0.5 - x[i]);
    } else if (x[i] > BOX_WIDTH_MINUS_HALF) {
      ax[i] = WALL_STIFFNESS * (BOX_WIDTH_MINUS_HALF - x[i]);
      wall_force -= ax[i];
      pE += 0.5 * WALL_STIFFNESS * (BOX_WIDTH_MINUS_HALF - x[i])
        * (BOX_WIDTH_MINUS_HALF - x[i]);
    } else {
      ax[i] = 0.0;
    }

    if (y[i] < 0.5) {
      ay[i] = WALL_STIFFNESS * (0.5 - y[i]);
      wall_force += ay[i];
      pE += 0.5 * WALL_STIFFNESS * (0.5 - y[i]) * (0.5 - y[i]);
    } else if (y[i] > BOX_HEIGHT_MINUS_HALF) {
      ay[i] = WALL_STIFFNESS * (BOX_HEIGHT_MINUS_HALF - y[i]);
      wall_force -= ay[i];
      pE += 0.5 * WALL_STIFFNESS * (BOX_HEIGHT_MINUS_HALF - y[i])
        * (BOX_HEIGHT_MINUS_HALF - y[i]);
    } else {
      ay[i] = 0.0;
    }
  }

  pressure = wall_force / (4.0 * BOX_WIDTH);

  // Next, compute interactions using the Lennard-Jones potential.
  for (i = 0; i < NUM_PARTICLES; i++) {
    for (j = 0; j < i; j++) {
      dx = x[i] - x[j];
      dx2 = dx * dx;
      // Make sure the pair are close enough to bother.
      if (dx2 < FORCE_CUTOFF_2) {
        dy = y[i] - y[j];
        dy2 = dy * dy;
        if (dy2 < FORCE_CUTOFF_2) {
          rSquared = dx2 + dy2;
          if (rSquared < FORCE_CUTOFF_2) {
            rSquaredInv = 1.0 / rSquared;
            attract = rSquaredInv * rSquaredInv * rSquaredInv;
            repel = attract * attract;
            pE += (4.0 * (repel - attract)) - P_EAT_CUTOFF;
            fOverR = 24.0 * ((2.0 * repel) - attract) * rSquaredInv;
            fx = fOverR * dx;
            fy = fOverR * dy;
            // Add the force on to i's acceleration.
            ax[i] += fx;
            ay[i] += fy;
            ax[j] -= fx;  // Newton's 3rd law
            ay[j] -= fy;
          }
        }
      }
    }
  }
}

/*
 * Reset accumulators and counters used for some measurements.
 */
void reset_measurements() {
  total_t = 0.0;
  total_p = 0.0;
  sample_count = 0;
}

/*
 * Calculate the instantaneous kinetic energy.
 *
 * @return The instantaneous kinetic energy.
 */
double calculate_kinetic_energy() {
  double ek = 0.0;
  for (int i = 0; i < NUM_PARTICLES; ++i) {
    ek += 0.5 * (vx[i] * vx[i] + vy[i] * vy[i]);
  }
  return ek;
}

/*
 * Compute accumulated property values from sampled values.
 */
void compute_properties() {
  sample_count++;
  kE = calculate_kinetic_energy();
  current_t = kE / NUM_PARTICLES;
  total_t += current_t;
  average_t = total_t / sample_count;
  total_p += pressure;
  average_p = total_p / sample_count;
}

/*
 * Return a small random number uniformly distributed around 0
 * with a maximum magnitude of epsilon.
 */
static int nudge_count = 0;
double nudge(double epsilon) {
  return (algo_647_uniform() - 0.5) * epsilon;
}

/*
 * Place the particles in empty locations of a grid on the simulation
 * box starting from the lower left. Applies a small amount of "jitter"
 * to break up the regularity of the layout a little. Does not check for
 * overflow of the grid (too many particles to fit).
 */
void place_particles() {
  const double spacing = 1.3; // Minimum space between particle centers.
  const double half_spacing = spacing / 2.0;
  const double bwmhs = BOX_WIDTH - half_spacing; // box width minus half spacing.
  const double jitter = 0.01; // Random amount to break up regularity.
  double x_pos = half_spacing;
  double y_pos = BOX_HEIGHT - half_spacing;

  for (int i = 0; i < NUM_PARTICLES; i++) {
    x[i] = x_pos + nudge(jitter);
    y[i] = y_pos + nudge(jitter);
    x_pos += spacing;
    if (x_pos > bwmhs) {
      x_pos = half_spacing;
      y_pos -= spacing;
    }
  }
}

/*
 * Re-scale the velocities to the target temperature.
 */
void rescale_velocities() {
  double scaling_factor;
  double velocity_squared_sum = 0.0;
  for (int i = 0; i < NUM_PARTICLES; i++) {
    velocity_squared_sum += vx[i] * vx[i] + vy[i] * vy[i];
  }
  scaling_factor = 2 * NUM_PARTICLES * TARGET_TEMPERATURE / velocity_squared_sum;
  scaling_factor = sqrt(scaling_factor);
  for (int i = 0; i < NUM_PARTICLES; i++) {
    vx[i] *= scaling_factor;
    vy[i] *= scaling_factor;
  }
}

/*
 * Remove any net momentum from the system of particles.
 */
void remove_drift() {
  for (int i = 0; i < NUM_PARTICLES; i++) {
    v_sum_x += vx[i];
    v_sum_y += vy[i];
  }
  for (int i = 0; i < NUM_PARTICLES; i++) {
    vx[i] -= v_sum_x / NUM_PARTICLES;
    vy[i] -= v_sum_y / NUM_PARTICLES;
  }
}

/*
 * Initialize particle positions, scale velocities to target temperature,
 * and remove systematic drift.
 */
void initialize_particles() {
  for (int i = 0; i < NUM_PARTICLES; i++) {
    x[i] = 0.0;
    y[i] = 0.0;
    vx[i] = gaussian_deviate_marsaglia() - 0.5;
    vy[i] = gaussian_deviate_marsaglia() - 0.5;
    ax[i] = 0.0;
    ay[i] = 0.0;
  }
  place_particles();
  rescale_velocities();
  remove_drift();
}

void print_config() {
  printf("\n Simulation Configuration:\n");
  printf("   Number of particles        : %d\n", NUM_PARTICLES);
  printf("   Cutoff radius              : %.2f\n", FORCE_CUTOFF);
  printf("   Target temperature         : %.3f\n", TARGET_TEMPERATURE);
  printf("   Integration step size      : %.3f\n", TIME_STEP);
  printf("   Width of simulation area   : %.1f\n", BOX_WIDTH);
  printf("   Height of simulation area  : %.1f\n", BOX_HEIGHT);
  printf("   Area of simulation         : %.1f\n", BOX_WIDTH * BOX_HEIGHT);
  printf("   Wall stiffness             : %.1f\n", WALL_STIFFNESS);
  printf("   Equilibration time         : %.1f\n", EQUILIBRATION_TIME);
  printf("      Steps between rescaling : %d\n", STEPS_BETWEEN_EQUIL_RESCALING);
  printf("   Production time            : %.1f\n", PRODUCTION_TIME);
  printf("   Steps per print out        : %d\n", STEPS_PER_PRINTOUT);
}


void print_properties_header() {
  printf("Time      Temp.   Pressure Tot. E.   Kin. E.   Pot. E.   Steps\n");
}

void print_properties() {
  printf("%7.3f,  %5.3f,  %6.4f,  %7.2f,  %6.2f,  %8.2f,  %d\n",
         t, average_t, average_p, (kE + pE), kE, pE, steps_accomplished);
}

void print_random_numbers() {
  printf("\nRandom numbers from a uniform distribution.\n");
  for (int i = 0; i < 10; i++) {
    printf("Uniform #%i: %.17f\n", i, algo_647_uniform());
  }
  printf("\nRandom numbers from a normal distribution.\n");
  for (int i = 0; i < 10; i++) {
    printf("Normal #%i: %.17f\n", i, gaussian_deviate_marsaglia());
  }
}

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

/*
 * Return a random uniform deviate in the range 0 <= x < 1.0 as a double.
 */
double algo_647_uniform() {
  int32_t k = floor(rng_seed/127773);
  int32_t partial = 16807*(rng_seed - k*127773) - k*2836;
  if (partial < 0)
    rng_seed = partial + INT32_MAX;
  else
    rng_seed = partial;

  return rng_seed*4.656612875e-10;
}

// The following is an implementation of the Marsaglia polar
// method for generating random standard normal deviates. See
// https://en.wikipedia.org/wiki/Marsaglia_polar_method.

/*
 * Return a random standard normal deviate as a double.
 */
double gaussian_deviate_marsaglia() {
  static int has_spare = 0;
  static double spare;
  double fac, rsq, r1, r2;
  if (has_spare) {
    has_spare = 0;
    return spare;
  } else {
    do {
      r1 = 2.0*algo_647_uniform() - 1.0;
      r2 = 2.0*algo_647_uniform() - 1.0;
      rsq = r1 * r1 + r2 * r2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac = sqrt(-2.0 * log(rsq) / rsq);
    spare = r1 * fac;
    has_spare = 1;
    return r2 * fac;
  }
}
