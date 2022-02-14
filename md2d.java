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

import java.io.PrintStream;
import java.time.Duration;
import java.time.Instant;

/**
 * A molecular dynamics simulation of a two-dimensional Lennard-Jones fluid.
 *
 * Originally downloaded from
 * http://www.personal.psu.edu/auk183/MolDynamics/LennardJones01.java on 15
 * January, 2015 before being _extensively_ modified.
 *
 * @author david
 */
public class md2d {

    // Simulation parameters
    static final int NUM_PARTICLES = 500;
    static final double FORCE_CUTOFF = 3.0;
    static final double FORCE_CUTOFF_2 = FORCE_CUTOFF * FORCE_CUTOFF;
    static final double P_EAT_CUTOFF = 4 * (Math.pow(FORCE_CUTOFF, -12)
                                            - Math.pow(FORCE_CUTOFF, -6));
    static final double TARGET_TEMPERATURE = 0.42;
    static final double TIME_STEP = 0.020;
    static final double DT_OVER_2 = 0.5 * TIME_STEP;
    static final double DT_SQUARED_OVER_2 = 0.5 * TIME_STEP * TIME_STEP;
    static final double BOX_WIDTH = 100.0;      // In units of the molecular diameter.
    static final double BOX_WIDTH_MINUS_HALF = BOX_WIDTH - 0.5;
    static final double BOX_HEIGHT = 100.0;
    static final double BOX_HEIGHT_MINUS_HALF = BOX_HEIGHT - 0.5;
    static final double WALL_STIFFNESS = 50.0;  // "Spring constant" for walls.
    static final double EQUILIBRATION_TIME = 100.0;
    static final int STEPS_BETWEEN_EQUIL_RESCALING = 10;
    static final double PRODUCTION_TIME = 300.0;
    static final int STEPS_PER_PRINTOUT = 50;

    // Other globals
    static final double[] x = new double[NUM_PARTICLES];
    static final double[] y = new double[NUM_PARTICLES];
    static final double[] vx = new double[NUM_PARTICLES];
    static final double[] vy = new double[NUM_PARTICLES];
    static final double[] ax = new double[NUM_PARTICLES];
    static final double[] ay = new double[NUM_PARTICLES];

    static double t = 0.0;
    static int stepsAccomplished = 0;

    static double vSumX = 0.0;
    static double vSumY = 0.0;
    static double kE = 0.0;
    static double pE = 0.0;
    static double currentT, totalT, averageT = 0.0;
    static double pressure, totalP, averageP = 0.0;
    static int sampleCount = 0;

    public static void main(String[] args) {

        pr.println("Md2dj - MD simulation of a 2D argon gas Lennard-Jones system.");
        pr.println("This version is written in Java and running version "
                   + System.getProperty("java.version") + ".");
        printConfig();

        initializeParticles();

        pr.println("\nEquilibration");
        printPropertiesHeader();
        printProperties();
        Instant start = Instant.now();
        while (t < EQUILIBRATION_TIME) {
            singleStep();
            stepsAccomplished++;
            t += TIME_STEP;
            if ((stepsAccomplished % STEPS_PER_PRINTOUT) == 0) {
                computeProperties();
                printProperties();
            }
            if ((stepsAccomplished % STEPS_BETWEEN_EQUIL_RESCALING) == 0) {
                rescaleVelocities();
            }
        }

        resetMeasurements();
        rescaleVelocities();
        pr.println("\nProduction");
        printPropertiesHeader();
        while (t < (EQUILIBRATION_TIME + PRODUCTION_TIME)) {
            singleStep();
            stepsAccomplished++;
            t += TIME_STEP;
            if ((stepsAccomplished % STEPS_PER_PRINTOUT) == 0) {
                computeProperties();
                printProperties();
            }
        }

        // Very crude calculation of steps per second. Includes print time.
        Instant finish = Instant.now();
        double elapsed = Duration.between(start, finish).toMillis() / 1000.0;
        pr.printf("Elapsed time: %.3f seconds, %.3f steps per second.\n",
                  elapsed, stepsAccomplished / elapsed);
        System.exit(0);
    }

    /**
     * Execute one time step using the Verlet algorithm.
     */
    static void singleStep() {

        for (int i = 0; i < NUM_PARTICLES; i++) {
            // Update positions.
            x[i] += (vx[i] * TIME_STEP) + (ax[i] * DT_SQUARED_OVER_2);
            y[i] += (vy[i] * TIME_STEP) + (ay[i] * DT_SQUARED_OVER_2);
            // Update velocities "halfway".
            vx[i] += (ax[i] * DT_OVER_2);
            vy[i] += (ay[i] * DT_OVER_2);
        }
        computeAccelerations();
        for (int i = 0; i < NUM_PARTICLES; i++) {
            // Finish updating velocities using the new accellerations.
            vx[i] += (ax[i] * DT_OVER_2);
            vy[i] += (ay[i] * DT_OVER_2);
        }
    }

    /**
     * Compute accelerations of the molecules from their current positions using
     * the Lennard-Jones potential.
     */
    static void computeAccelerations() {

        int i, j;
        double dx, dy;
        double dx2, dy2, rSquared, rSquaredInv, attract, repel, fOverR, fx, fy;
        double wallForce = 0.0;

        pE = 0.0;

        // First, check for bouncing against the walls.
        for (i = 0; i < NUM_PARTICLES; i++) {
            if (x[i] < 0.5) {
                ax[i] = WALL_STIFFNESS * (0.5 - x[i]);
                wallForce += ax[i];
                pE += 0.5 * WALL_STIFFNESS * (0.5 - x[i]) * (0.5 - x[i]);
            } else if (x[i] > BOX_WIDTH_MINUS_HALF) {
                ax[i] = WALL_STIFFNESS * (BOX_WIDTH_MINUS_HALF - x[i]);
                wallForce -= ax[i];
                pE += 0.5 * WALL_STIFFNESS * (BOX_WIDTH_MINUS_HALF - x[i])
                    * (BOX_WIDTH_MINUS_HALF - x[i]);
            } else {
                ax[i] = 0.0;
            }

            if (y[i] < 0.5) {
                ay[i] = WALL_STIFFNESS * (0.5 - y[i]);
                wallForce += ay[i];
                pE += 0.5 * WALL_STIFFNESS * (0.5 - y[i]) * (0.5 - y[i]);
            } else if (y[i] > BOX_HEIGHT_MINUS_HALF) {
                ay[i] = (WALL_STIFFNESS * (BOX_HEIGHT_MINUS_HALF - y[i]));
                wallForce -= ay[i];
                pE += 0.5 * WALL_STIFFNESS * (BOX_HEIGHT_MINUS_HALF - y[i])
                    * (BOX_HEIGHT_MINUS_HALF - y[i]);
            } else {
                ay[i] = 0.0;
            }
        }

        pressure = wallForce / (4.0 * BOX_WIDTH);

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
                            fOverR = 24.0 * ((2.0 * repel) - attract)
                                * rSquaredInv;
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

    /**
     * Reset accumulators and counters used for some measurements.
     */
    static void resetMeasurements() {
        totalT = 0.0;
        totalP = 0.0;
        sampleCount = 0;
    }

    /**
     * Calculate the instantaneous kinetic energy.
     *
     * @return The instantaneous kinetic energy.
     */
    static double calculateKineticEnergy() {
        double ek = 0.0;
        for (int i = 0; i < NUM_PARTICLES; ++i) {
            ek += 0.5 * (vx[i] * vx[i] + vy[i] * vy[i]);
        }
        return ek;
    }

    /**
     * Compute accumulated property values from sampled values.
     */
    static void computeProperties() {
        sampleCount++;
        kE = calculateKineticEnergy();
        currentT = kE / NUM_PARTICLES;
        totalT += currentT;
        averageT = totalT / sampleCount;
        totalP += pressure;
        averageP = totalP / sampleCount;
    }

    /**
     * Return a small random number uniformly distributed around 0 with a
     * maximum magnitude of epsilon.
     */
    static double nudge(double epsilon) {
        return (algo647Uniform() - 0.5) * epsilon;
    }

    /**
     * Place the particles in empty locations of a grid on the simulation box
     * starting from the lower left. Applies a small amount of "jitter" to break
     * up the regularity of the layout a little. Does not check for overflow of
     * the grid (too many particles to fit).
     */
    static void placeParticles() {
        double spacing = 1.3; // Minimum space between particle centers.
        double halfSpacing = spacing / 2.0;
        double bwmhs = BOX_WIDTH - halfSpacing; // box width minus half spacing.
        double jitter = 0.01; // Random amount to break up regularity.
        double xPos = halfSpacing;
        double yPos = BOX_HEIGHT - halfSpacing;

        for (int i = 0; i < NUM_PARTICLES; i++) {
            x[i] = xPos + nudge(jitter);
            y[i] = yPos + nudge(jitter);
            xPos += spacing;
            if (xPos > bwmhs) {
                xPos = halfSpacing;
                yPos -= spacing;
            }
        }
    }

    /**
     * Re-scale the velocities to the target temperature.
     */
    static void rescaleVelocities() {
        double scalingFactor;
        double velocitySquaredSum = 0.0;
        for (int i = 0; i < NUM_PARTICLES; i++) {
            velocitySquaredSum += vx[i] * vx[i] + vy[i] * vy[i];
        }
        scalingFactor = 2 * NUM_PARTICLES * TARGET_TEMPERATURE / velocitySquaredSum;
        scalingFactor = Math.sqrt(scalingFactor);
        for (int i = 0; i < NUM_PARTICLES; i++) {
            vx[i] *= scalingFactor;
            vy[i] *= scalingFactor;
        }
    }

    /**
     * Remove any net momentum from the system of particles.
     */
    static void removeDrift() {
        for (int i = 0; i < NUM_PARTICLES; i++) {
            vSumX += vx[i];
            vSumY += vy[i];
        }
        for (int i = 0; i < NUM_PARTICLES; i++) {
            vx[i] -= vSumX / NUM_PARTICLES;
            vy[i] -= vSumY / NUM_PARTICLES;
        }
    }

    /**
     * Initialize particle positions, scale velocities to target temperature,
     * and remove systematic drift.
     */
    static void initializeParticles() {
        for (int i = 0; i < NUM_PARTICLES; i++) {
            x[i] = 0.0;
            y[i] = 0.0;
            vx[i] = gaussianDeviateMarsaglia() - 0.5; //Math.random() - 0.5;
            vy[i] = gaussianDeviateMarsaglia() - 0.5; //Math.random() - 0.5;
            ax[i] = 0.0;
            ay[i] = 0.0;
        }
        placeParticles();
        rescaleVelocities();
        removeDrift();
    }

    //
    // Printing stuff.
    //
    // Just to shorten lines and reduce some typing.
    static PrintStream pr = System.out;

    static void printConfig() {
        pr.println("\n Simulation Configuration:");
        pr.println("   Number of particles        : " + NUM_PARTICLES);
        pr.println("   Cutoff radius              : " + FORCE_CUTOFF);
        pr.println("   Target temperature         : " + TARGET_TEMPERATURE);
        pr.println("   Integration step size      : " + TIME_STEP);
        pr.println("   Width of simulation area   : " + BOX_WIDTH);
        pr.println("   Height of simulation area  : " + BOX_HEIGHT);
        pr.println("   Area of simulation         : " + BOX_WIDTH * BOX_HEIGHT);
        pr.println("   Wall stiffness             : " + WALL_STIFFNESS);
        pr.println("   Equilibration time         : " + EQUILIBRATION_TIME);
        pr.println("      Steps between rescaling : " + STEPS_BETWEEN_EQUIL_RESCALING);
        pr.println("   Production time            : " + PRODUCTION_TIME);
        pr.println("   Steps per print out        : " + STEPS_PER_PRINTOUT);
    }

    static void printPropertiesHeader() {
        pr.printf("Time      Temp.   Pressure Tot. E.   Kin. E.   Pot. E.   Steps\n");
    }

    static void printProperties() {
        pr.printf("%7.3f,  %5.3f,  %6.4f,  %7.2f,  %6.2f,  %8.2f,  %d\n",
                  t, averageT, averageP, (kE + pE), kE, pE, stepsAccomplished);
    }

    /**
     * Print some pseudo random numbers. Just for debugging purposes.
     */
    static void printRandomNumbers() {
        pr.println("\nRandom numbers from a uniform distribution");
        for (int i = 0; i < 10; i++) {
            double rnu = algo647Uniform();
            pr.println("Uniform #" + i + ": " + rnu);
        }

        pr.println("\nRandom numbers from a normal distribution.");
        for (int i = 0; i < 10; i++) {
            double rnn = gaussianDeviateMarsaglia();
            pr.println("Normal #" + i + ": " + rnn);
        }
    }

    //
    // Random number stuff.
    //
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
    static int rngSeed = 12345;

    /**
     * Return a random uniform deviate in the range [0, 1) as a double.
     *
     * @return Random uniform deviate in the range [0, 1) as a double.
     */
    static double algo647Uniform() {
        int k = Math.floorDiv(rngSeed, 127773);
        int partial = 16807 * (rngSeed - k * 127773) - k * 2386;
        if (partial < 0) {
            rngSeed = partial + Integer.MAX_VALUE;
        } else {
            rngSeed = partial;
        }
        return rngSeed * 4.656612875e-10;
    }

    // The following is an implementation of the Marsaglia polar
    // method for generating random standard normal deviates. See
    // https://en.wikipedia.org/wiki/Marsaglia_polar_method.
    private static boolean hasSpare = false;
    private static double spare;

    /**
     * Return a random standard normal deviate as a double.
     *
     * @return A random standard normal deviate as a double.
     */
    static double gaussianDeviateMarsaglia() {
        double fac, rsq, r1, r2;
        if (hasSpare) {
            hasSpare = false;
            return spare;
        } else {
            do {
                r1 = 2.0 * algo647Uniform() - 1.0;
                r2 = 2.0 * algo647Uniform() - 1.0;
                rsq = r1 * r1 + r2 * r2;
            } while (rsq >= 1.0 || rsq == 0.0);
            fac = Math.sqrt(-2.0 * Math.log(rsq) / rsq);
            spare = r1 * fac;
            hasSpare = true;
            return r2 * fac;
        }
    }
}
