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
# A molecular dynamics simulation of a two-dimensional Lennard-Jones fluid.
#
# Author: david
#

# Simulation parameters
(def num-particles 500)
(def force-cutoff 3.0)
(def force-cutoff-2 (* force-cutoff force-cutoff))
(def p-eat-cutoff (* 4.0 (- (math/pow force-cutoff -12.0) (math/pow force-cutoff -6.0))))
(def target-temperature 0.42)
(def time-step 0.020)
(def dt-over-2 (* 0.5 time-step))
(def dt-squared-over-2 (* 0.5 time-step time-step))
(def box-width 100.0)
(def box-width-minus-half (- box-width 0.5))
(def box-height 100.0)
(def box-height-minus-half (- box-height 0.5))
(def wall-stiffness 50.0)
(def equilibration-time 200)
(def steps-between-equil-rescaling 10)
(def production-time 200.0)
(def steps-between-prod-rescaling 100)
(def steps-per-printout 50)

# Other globals
(def x (array/new num-particles))
(def y (array/new num-particles))
(def vx (array/new num-particles))
(def vy (array/new num-particles))
(def ax (array/new num-particles))
(def ay (array/new num-particles))

(var t 0.0)
(var steps-accomplished 0)
(var v-sum-x 0.0)
(var v-sum-y 0.0)
(var k-e 0.0)
(var p-e 0.0)
(var current-T 0.0)
(var total-T 0.0)
(var average-T 0.0)
(var pressure 0.0)
(var total-P 0.0)
(var average-P 0.0)
(var sample-count 0)

(var rng-seed :private 12345)

# "forward" declaration of functions.
(varfn print-config [] nil)
(varfn initialize-properties [] nil)
(varfn print-properties-header [] nil)
(varfn print-properties [] nil)
(varfn single-step [] nil)
(varfn compute-accelerations [] nil)
(varfn compute-properties [] nil)
(varfn rescale-velocities [] nil)
(varfn reset-measurements [] nil)
(varfn calculate-kinetic-energy [] nil)
(varfn print-random-numbers [] nil)
(varfn initialize-particles [] nil)
(varfn initialize-measurements [] nil)
(varfn measure-properties [] nil)
(varfn algo-647-uniform [] nil)
(varfn gaussian-deviate-marsaglia [] nil)
(varfn print-uniform-randoms [] nil)

(defn main []
  (print "md2d - MD simulation of a 2D argon gas Lennard-Jones system.")
  (printf "This version is written in Janet and running version %s" janet/version)
  (print-config)

  (initialize-particles)

  (print "\nEquilibration")
  (print-properties-header)
  (print-properties)
  (def start-time (os/clock))
  (while (< t equilibration-time)
    (single-step)
    (++ steps-accomplished)
    (set t (+ t time-step))
    (if (= 0 (% steps-accomplished steps-per-printout))
      (do (compute-properties)
          (print-properties)))
    (if (= 0 (% steps-accomplished steps-between-equil-rescaling))
      (rescale-velocities)))

  (reset-measurements)
  (rescale-velocities)
  (print "\nProduction")
  (print-properties-header)
  (while (< t (+ equilibration-time production-time))
    (single-step)
    (++ steps-accomplished)
    (set t (+ t time-step))
    (if (= 0 (% steps-accomplished steps-per-printout))
      (do (compute-properties)
          (print-properties)))
    (if (= 0 (% steps-accomplished steps-between-prod-rescaling))
      (rescale-velocities)))

  # Very crude calculation of steps per second. Includes print time.
  (def elapsed-time (- (os/clock) start-time))
  (printf "Elapsed time: %.3f seconds, %.3f steps per second."
          elapsed-time (/ steps-accomplished elapsed-time))
  (os/exit 0))

(varfn single-step
  "Execute one time step using the Verlet algorithm."
  []
  (for i 0 num-particles
    # Update positions.
    (put x i (+ (get x i) (* (get vx i) time-step) (* (get ax i) dt-squared-over-2)))
    (put y i (+ (get y i) (* (get vy i) time-step) (* (get ay i) dt-squared-over-2)))
    # Update velocities "halfway".
    (put vx i (+ (get vx i) (* (get ax i) dt-over-2)))
    (put vy i (+ (get vy i) (* (get ay i) dt-over-2))))
  (compute-accelerations)
  (for i 0 num-particles
    # Finish updating velocities using the new accelerations.
    (put vx i (+ (get vx i) (* (get ax i) dt-over-2)))
    (put vy i (+ (get vy i) (* (get ay i) dt-over-2)))))

(varfn compute-accelerations
  "Compute accelerations of the molecules from their current positions using
the Lennard-Jones potential."
  []
  (var dx 0.0) (var dy 0.0)
  (var dx2 0.0) (var dy2 0.0) (var r-squared 0.0) (var r-squared-inv 0.0)
  (var attract 0.0) (var repel 0.0) (var f-over-r 0.0) (var fx 0.0) (var fy 0.0)
  (var wall-force 0.0)

  (set p-e 0.0)

  # First, check for pressure from particles bouncing off the walls.
  (for i 0 num-particles
    (var xi (get x i))
    (cond
      (< xi 0.5)
      (do
        (put ax i (* wall-stiffness (- 0.5 xi)))
        (+= wall-force (get ax i))
        (+= p-e (* 0.5 wall-stiffness (- 0.5 xi) (- 0.5 xi))))

      (> xi box-width-minus-half)
      (do
        (put ax i (* wall-stiffness (- box-width-minus-half xi)))
        (-= wall-force (get ax i))
        (+= p-e (* 0.5 wall-stiffness (- box-width-minus-half xi)
                   (- box-width-minus-half xi))))

      (put ax i 0.0))

    (var yi (get y i))
    (cond
      (< yi 0.5)
      (do
        (put ay i (* wall-stiffness (- 0.5 yi)))
        (+= wall-force (get ay i))
        (+= p-e (* 0.5 wall-stiffness (- 0.5 yi) (- 0.5 yi))))

      (> yi box-height-minus-half)
      (do
        (put ay i (* wall-stiffness (- box-height-minus-half yi)))
        (-= wall-force (get ay i))
        (+= p-e (* 0.5 wall-stiffness (- box-height-minus-half  yi)
                   (- box-height-minus-half yi))))

      (put ay i 0.0)))

  (set pressure (/ wall-force (* 4.0 box-width)))

  # Next, compute interactions using the Lennard-Jones potential.
  (for i 0 num-particles
    (for j 0 i
      (set dx (- (get x i) (get x j)))
      (set dx2 (* dx dx))
      # Make sure the pair are close enough to bother.
      (if  (< dx2 force-cutoff-2)
        (do
          (set dy (- (get y i) (get y j)))
          (set dy2 (* dy dy))
          (if (< dy2 force-cutoff-2)
            (do
              (set r-squared (+ dx2 dy2))
              (if (< r-squared force-cutoff-2)
                (do
                  (set r-squared-inv (/ 1.0 r-squared))
                  (set attract (* r-squared-inv r-squared-inv r-squared-inv))
                  (set repel (* attract attract))
                  (+= p-e (- (* 4.0 (- repel attract)) p-eat-cutoff))
                  (set f-over-r (* 24.0 (- (* 2.0 repel) attract) r-squared-inv))
                  (set fx (* f-over-r dx))
                  (set fy (* f-over-r dy))
                  # Add the force on to i's accelerations.
                  (put ax i (+ (get ax i) fx))
                  (put ay i (+ (get ay i) fy))
                  (put ax j (- (get ax j) fx)) # Newton's 3rd law
                  (put ay j (- (get ay j) fy)))))))))))

(varfn reset-measurements
  "Reset accumulators and counters for some measurements."
  []
  (set total-T 0.0)
  (set total-P 0.0)
  (set sample-count 0))

(defn calculate-kinetic-energy
  "Calculate the instantaneous kinetic energy."
  []
  (var ek 0.0)
  (for i 0 num-particles
    (+= ek (* 0.5 (+ (* (get vx i) (get vx i)) (* (get vy i) (get vy i))))))
  ek)

(varfn compute-properties
  "Compute accumulated property values from sampled values."
  []
  (++ sample-count)
  (set k-e (calculate-kinetic-energy))
  (set current-T (/ k-e num-particles))
  (+= total-T current-T)
  (set average-T (/ total-T sample-count))
  (+= total-P pressure)
  (set average-P (/ total-P sample-count)))

(defn nudge
  "Return a small random number unformly distributed around 0
with a maximum magniture of epsilon."
  [epsilon]
  (* epsilon (- (algo-647-uniform) 0.5)))

(defn place-particles
  "Place the particles in empty locations of a grid on the simulation
box starting from the lower left. Applies a small amount of 'jitter'
to break up the regularity of the layout a little. Does not check for
overflow of the grid (too many particles to fiel)."
  []
  (def spacing 1.3) # Minimum space between particles.
  (def half-spacing (/ spacing 2.0))
  (def bwmhs (- box-width half-spacing))
  (def jitter 0.01) # Random amount to break up regularity.
  (var x-pos half-spacing)
  (var y-pos (- box-height half-spacing))

  (for i 0 num-particles
    (put x i (+ x-pos (nudge jitter)))
    (put y i (+ y-pos (nudge jitter)))
    (+= x-pos spacing)
    (if (> x-pos bwmhs)
      (do
        (set x-pos half-spacing)
        (-= y-pos spacing)))))

(varfn rescale-velocities
  "Re-scale the velocities to the target temperature."
  []
  (var scaling-factor 0.0)
  (var velocity-squared-sum 0.0)
  (for i 0 num-particles
    (+= velocity-squared-sum (+ (* (get vx i) (get vx i)) (* (get vy i) (get vy i)))))

  (var scaling-factor (/ (* 2.0 num-particles target-temperature) velocity-squared-sum))
  (set scaling-factor (math/sqrt scaling-factor))

  (for i 0 num-particles
    (put vx i (* (get vx i) scaling-factor))
    (put vy i (* (get vy i) scaling-factor))))

(defn remove-drift
  "Remove any net momentum from the system of particles."
  []
  (for i 0 num-particles
    (+= v-sum-x (get vx i))
    (+= v-sum-y (get vy i)))
  (def x-adj (/ v-sum-x num-particles))
  (def y-adj (/ v-sum-y num-particles))
  (for i 0 num-particles
    (put vx i (- (get vx i) x-adj))
    (put vy i (- (get vy i) y-adj))))

(varfn initialize-particles
  "Initialize particle positions, scale velocities to target temperature,
and remove systematic drift."
  []
  (for i 0 num-particles
    (put x i 0.0)
    (put y i 0.0)
    (put vx i (- (gaussian-deviate-marsaglia) 0.5))
    (put vy i (- (gaussian-deviate-marsaglia) 0.5))
    (put ax i 0.0)
    (put ay i 0.0))
  (place-particles)
  (rescale-velocities)
  (remove-drift))

#
# Printing stuff
#

(varfn print-config []
  (print "\n Simulation Configuration:")
  (printf "   Number of particles        : %d" num-particles)
  (printf "   Cutoff radius              : %.2f" force-cutoff)
  (printf "   Target temperature         : %.3f" target-temperature)
  (printf "   Integration step size      : %.3f" time-step)
  (printf "   Width of simulation area   : %.1f" box-width)
  (printf "   Height of sumulation area  : %.1f" box-height)
  (printf "   Area of simulation         : %.1f" (* box-width box-height))
  (printf "   Wall stiffness             : %.1f" wall-stiffness)
  (printf "   Equilibration time         : %.1f" equilibration-time)
  (printf "      Steps between rescaling : %d" steps-between-equil-rescaling)
  (printf "   Production time            : %.1f" production-time)
  (printf "      Steps between rescaling : %d" steps-between-prod-rescaling)
  (printf "   Steps per printout         : %d" steps-per-printout))

(varfn print-properties-header []
  (print "Time      Temp.   Pressure Tot. E.   Kin. E.   Pot. E.   Steps"))

(varfn print-properties []
  (printf "%7.3f,  %5.3f,  %6.4f,  %7.2f,  %6.2f,  %8.2f,  %d"
          t average-T average-P (+ k-e p-e) k-e p-e steps-accomplished))

(varfn print-random-numbers []
  (print "\nRandom numbers from a uniform distribution.")
  (var i 0)
  (while (< i 10)
    (printf "Uniform #%d: %.17f" i (algo-647-uniform))
    (++ i))
  (print "\nRandom numbers from a normal distribution.")
  (set i 0)
  (while (< i 10)
    (printf "Normal #%d: %.17f" i (gaussian-deviate-marsaglia))
    (++ i)))

#
# Random number stuff.
#

# The algo647Uniform function returns a random uniform deviate in
# the range [0, 1) calculated using the ACM algorithm 647 as
# recorded at: http://calgo.acm.org.
# The algorithm is from the paper "ALGORITHM 647: Implementation and
# relative Efficiency of Quasirandom Sequence Generators" by Bennet L. Fox,
# ACM Transactions on Mathematical Software (TOMS), Volume 12, No 4,
# Dec 1986, pp. 362-376.
#
# I have never seen the actual paper - it's behind a paywall - but the
# FORTRAN version of the algorithm is easy enough to understand.
(def int32-max 2147483647)
(varfn algo-647-uniform
  "Return a uniform deviate in the range [0, 1)."
  []
  (let [k (math/floor(/ rng-seed 127773))
        partial (- (* 16807 (- rng-seed (* k 127773))) (* k 2836))]
    (if (< partial 0)
      (set rng-seed (+ partial int32-max))
      (set rng-seed partial))
    (* rng-seed 4.656612875e-10)))

# A more "lispy" version of the uniform random number generator might be:
#(def rng ((fn []
#            (def int32-max 2147483647)
#            (var seed 12345)
#            (fn [] (let [k (math/floor(/ seed 127773))
#                         partial (- (* 16807 (- seed (* k 127773))) (* k 2386))]
#                     (if (< partial 0)
#                       (set seed (+ partial int32-max))
#                       (set seed partial))
#                     (* seed 4.656612875e-10))))))
# Note that the seed variable is captured as a closure.

# The following is an implementation of the Marsaglia polar
# method for generating random standard deviates. Seed
# https://en.wikipedia.org/wiki/Marsaglia_polar_method.

(var has-spare false)
(var spare 0.0)

(varfn gaussian-deviate-marsaglia
  "Return a random standard normal deviate."
  []
  (if has-spare
    (do
      (set has-spare false)
      spare)
    (do
      (var fac 0.0)
      (var rsq 0.0)
      (var r1 0.0)
      (var r2 0.0)
      (while true
        (set r1 (- (* 2.0 (algo-647-uniform)) 1.0))
        (set r2 (- (* 2.0 (algo-647-uniform)) 1.0))
        (set rsq (+ (* r1 r1) (* r2 r2)))
        (if (not (or (>= rsq 1.0) (= rsq 0.0)))
          (break)))
      (var fac (math/sqrt (/ (* -2.0 (math/log rsq)) rsq)))
      (set spare (* r1 fac))
      (set has-spare true)
      (* r2 fac))))

(main)
