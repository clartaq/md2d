;;
;; The MIT License
;;
;; Copyright 2022 David D. Clark.
;;
;; Permission is hereby granted, free of charge, to any person obtaining a copy
;; of this software and associated documentation files (the "Software"), to deal
;; in the Software without restriction, including without limitation the rights
;; to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
;; copies of the Software, and to permit persons to whom the Software is
;; furnished to do so, subject to the following conditions:
;;
;; The above copyright notice and this permission notice shall be included in
;; all copies or substantial portions of the Software.
;;
;; THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
;; IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
;; FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
;; AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
;; LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
;; OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
;; THE SOFTWARE.
;;

;;
;; A molecular dynamics simulation of a two-dimensional Lennard-Jones fluid.
;;
;; Author: david
;;

(import
 (scheme))

;; Simulation parameters
(define num-particles 500)
(define force-cutoff 3.0)
(define force-cutoff-2 (* force-cutoff force-cutoff))
(define p-eat-cutoff (* 4.0 (- (flexpt force-cutoff -12.0) (flexpt force-cutoff -6.0))))
(define target-temperature 0.42)
(define time-step 0.020)
(define dt-over-2 (* 0.5 time-step))
(define dt-squared-over-2 (* 0.5 time-step time-step))
(define box-width 100.0)
(define box-width-minus-half (- box-width 0.5))
(define box-height 100.0)
(define box-height-minus-half (- box-height 0.5))
(define wall-stiffness 50.0)
(define equilibration-time 200.0)
(define steps-between-equil-rescaling 10)
(define production-time 200.0)
(define steps-between-prod-rescaling 100)
(define steps-per-printout 50)

;; Global variables
(define x (make-vector num-particles))
(define y (make-vector num-particles))
(define vx (make-vector num-particles))
(define vy (make-vector num-particles))
(define ax (make-vector num-particles))
(define ay (make-vector num-particles))

(define t 0.0)
(define steps-accomplished 0)
(define v-sum-x 0.0)
(define v-sum-y 0.0)
(define k-e 0.0)
(define p-e 0.0)
(define current-T 0.0)
(define total-T 0.0)
(define average-T 0.0)
(define pressure 0.0)
(define total-P 0.0)
(define average-P 0.0)
(define sample-count 0)

;; Give Scheme some more imperative incrementing, decrementing and looping
;; functionality C.

;; Scalar mutating increment and decrement since Scheme doesn't normally have
;; these for some reason. The two argument versions are somewhat like += and -=.

(define-syntax inc!
  (syntax-rules ()
    ((_ x)   (begin (set! x (+ x 1)) x))
    ((_ x n) (begin (set! x (+ x n)) x))))

(define-syntax dec!
  (syntax-rules ()
    ((_ x)   (inc! x -1))
    ((_ x n) (inc! x (- n)))))

;; Mutating increment/decrement for vectors.
(define-syntax vector-inc!
  (syntax-rules ()
    ((_ v i) (begin (vector-set! v i (+ (vector-ref v i) 1)) (vector-ref v i)))
    ((_ v i n) (begin (vector-set! v i (+ (vector-ref v i) n)) (vector-ref v i)))))

(define-syntax vector-dec!
  (syntax-rules ()
    ((_ v i) (vector-inc! v i -1))
    ((_ v i n) (vector-inc! v i (- n)))))

;; Some looping macros that look like the C equivalents.
(define-syntax while
  (syntax-rules ()
    ((_ pred steps ...)
     (let loop () (when pred steps ... (loop))))))

;; The for loop only increments the loop variable by one on each iteration.
(define-syntax for
  (syntax-rules ()
    ((_ (var from to) steps ...)
     (let loop ([var from])
       (when (< var to)
	 steps ...
	 (loop (+ var 1)))))))

(define (main)
  (let ([start-time 0]
        [elapsed-time 0.0])
    (println "md2d - MD simulation of a 2D argon gas Lennard-Jones system.")
    (println "This version is written in Scheme and running Chez Scheme.")
    (print-config)

    (initialize-particles)

    (println "\nEquilibration")
    (print-properties-header)
    (print-properties)
    (set! start-time (cpu-time))
    (while (< t equilibration-time)
      (single-step)
      (inc! steps-accomplished)
      (inc! t time-step)
      (if (zero? (modulo steps-accomplished steps-per-printout))
          (begin (compute-properties)
                 (print-properties)))
      (if (zero? (modulo steps-accomplished steps-between-equil-rescaling))
          (rescale-velocities)))

    (reset-measurements)
    (rescale-velocities)
    (println "\nProduction")
    (print-properties-header)
    (while (< t (+ equilibration-time production-time))
      (single-step)
      (inc! steps-accomplished)
      (inc! t time-step)
      (if (zero? (modulo steps-accomplished steps-per-printout))
          (begin (compute-properties)
                 (print-properties)))
      (if (zero? (modulo steps-accomplished steps-between-prod-rescaling))
          (rescale-velocities)))

    ;; Very crude calculation of steps per second. Includes print time.
    (set! elapsed-time (/ (- (cpu-time) start-time) 1000.0))
    (println "Elapsed time: " elapsed-time)
    (printf "Elapsed time: ~,3f seconds, ~,3f steps per second.~%"
            elapsed-time (/ steps-accomplished elapsed-time))
    (exit 0)))

(define (single-step)
  "Execute one time step using the Verlet algorithm."
  (for (i 0 num-particles)
       ;; Update positions.
       (vector-inc! x i (+ (* (vector-ref vx i) time-step) (* (vector-ref ax i) dt-squared-over-2)))
       (vector-inc! y i (+ (* (vector-ref vy i) time-step) (* (vector-ref ay i) dt-squared-over-2)))
       ;; Update velocities "halfway".
       (vector-inc! vx i (* (vector-ref ax i) dt-over-2))
       (vector-inc! vy i (* (vector-ref ay i) dt-over-2)))
  (compute-accelerations)
  (for (i 0 num-particles)
       ;; Finish updating velocities using the new accelerations.
       (vector-inc! vx i (* (vector-ref ax i) dt-over-2))
       (vector-inc! vy i (* (vector-ref ay i) dt-over-2))))

(define (compute-accelerations)
  "Compute accelerations of the molecules from their current positions using
the Lennard-Jones potential."
  (set! p-e 0.0)

  ;; First, check for bouncing against the walls.
  (let ([wall-force 0.0])
    (for (i 0 num-particles)
         (let ([xi (vector-ref x i)])
           (cond
            [(< xi 0.5)
             (begin
               (vector-set! ax i (* wall-stiffness (- 0.5 xi)))
               (inc! wall-force (vector-ref ax i))
               (inc! p-e  (* 0.5 wall-stiffness (- 0.5 xi) (- 0.5 xi))))]

            [(> xi box-width-minus-half)
             (begin
               (vector-set! ax i (* wall-stiffness (- box-width-minus-half xi)))
               (dec! wall-force (vector-ref ax i))
               (inc! p-e (* 0.5 wall-stiffness (- box-width-minus-half xi)
                            (- box-width-minus-half xi))))]

            [else (vector-set! ax i 0.0)])))

    (for (i 0 num-particles)
         (let ([yi (vector-ref y i)])
           (cond
            [(< yi 0.5)
             (begin
               (vector-set! ay i (* wall-stiffness (- 0.5 yi)))
               (inc! wall-force (vector-ref ay i))
               (inc! p-e (* 0.5 wall-stiffness (- 0.5 yi) (- 0.5 yi))))]

            [(> yi box-height-minus-half)
             (begin
               (vector-set! ay i (* wall-stiffness (- box-height-minus-half yi)))
               (dec! wall-force (vector-ref ay i))
               (inc! p-e (* 0.5 wall-stiffness (- box-height-minus-half yi)
                            (- box-height-minus-half yi))))]

            [else (vector-set! ay i 0.0)])))

    (set! pressure (/ wall-force (* 4.0 box-width))))

  ;; Next, compute interactions using the Lennard-Jones potential.
  (let ([dx 0.0]
        [dy 0.0]
        [dx2 0.0]
        [dy2 0.0]
        [r-squared 0.0]
        [r-squared-inv 0.0]
        [attract 0.0]
        [repel 0.0]
        [f-over-r 0.0]
        [fx 0.0]
        [fy 0.0])
    (for (i 0 num-particles)
         (for (j 0 i)
              (set! dx (- (vector-ref x i) (vector-ref x j)))
              (set! dx2 (* dx dx))
              ;; Make sure the pair are close enough to bother.
              (if (< dx2 force-cutoff-2)
                  (begin
                    (set! dy (- (vector-ref y i) (vector-ref y j)))
                    (set! dy2 (* dy dy))
                    (if (< dy2 force-cutoff-2)
                        (begin
                          (set! r-squared (+ dx2 dy2))
                          (if (< r-squared force-cutoff-2)
                              (begin
                                (set! r-squared-inv (/ 1.0 r-squared))
                                (set! attract (* r-squared-inv r-squared-inv r-squared-inv))
                                (set! repel (* attract attract))
                                (inc! p-e (- (* 4.0 (- repel attract)) p-eat-cutoff))
                                (set! f-over-r (* 24.0 (- (* 2.0 repel) attract) r-squared-inv))
                                (set! fx (* f-over-r dx))
                                (set! fy (* f-over-r dy))
                                (vector-inc! ax i fx)
                                (vector-inc! ay i fy)
                                (vector-dec! ax j fx)
                                (vector-dec! ay j fy)))))))))))

(define (reset-measurements)
  "Reset accumulators and counters for some measurements."
  (set! total-T 0.0)
  (set! total-P 0.0)
  (set! sample-count 0))

(define (calculate-kinetic-energy)
  "Calculate the instantaneous kinetic energy."
  (let ([ek 0])
    (for (i 0 num-particles)
         (inc! ek  (* 0.5 (+ (* (vector-ref vx i) (vector-ref vx i))
                             (* (vector-ref vy i) (vector-ref vy i))))))
    ek))

(define (compute-properties)
  "Compute accumulated property values from sampled values."
  (inc! sample-count)
  (set! k-e (calculate-kinetic-energy))
  (set! current-T (/ k-e num-particles))
  (inc! total-T current-T)
  (set! average-T (/ total-T sample-count))
  (inc! total-P pressure)
  (set! average-P (/ total-P sample-count)))

(define (nudge epsilon)
  "Return a small random number uniformly distributed around 0
with a maximum magnitude of epsilon."
  (* epsilon (- (algo-647-uniform) 0.5)))

(define (place-particles)
  "Place the particles in empty locations of a grid on the simulation
box starting from the lower left. Applies a small amount of 'jitter'
to break up the regularity of the layout a little. Does not check for
overflow of the grid (too many particles to fiel)."
  (let* ([spacing 1.3] ;; Minimum space between particles.
         [half-spacing (/ spacing 2.0)]
         [bwmhs (- box-width half-spacing)]
         [jitter 0.01] ;; Random amount to break up regularity.
         [x-pos half-spacing]
         [y-pos (- box-height half-spacing)])
    (for (i 0 num-particles)
         (vector-set! x i (+ x-pos (nudge jitter)))
         (vector-set! y i (+ y-pos (nudge jitter)))
         (inc! x-pos spacing)
         (if (> x-pos bwmhs)
             (begin
               (set! x-pos half-spacing)
               (dec! y-pos spacing))))))

(define (rescale-velocities)
  "Re-scale the velocities to the target temperature."
  (let ([scaling-factor 0.0]
        [velocity-squared-sum 0.0])
    (for (i 0 num-particles)
         (inc! velocity-squared-sum
               (+ (* (vector-ref vx i) (vector-ref vx i))
                  (* (vector-ref vy i) (vector-ref vy i)))))
    (set! scaling-factor (/ (* 2.0 num-particles target-temperature)
                            velocity-squared-sum))
    (set! scaling-factor (sqrt scaling-factor))
    (for (i 0 num-particles)
         (vector-set! vx i (* (vector-ref vx i) scaling-factor))
         (vector-set! vy i (* (vector-ref vy i) scaling-factor)))))

(define (remove-drift)
  "Remove any net momentum from the system of particles."
  (for (i 0 num-particles)
       (inc! v-sum-x (vector-ref vx i))
       (inc! v-sum-y (vector-ref vy i)))
  (let ([x-adj (/ v-sum-x num-particles)]
        [y-adj (/ v-sum-y num-particles)])
    (for (i 0 num-particles)
         (vector-dec! vx i x-adj)
         (vector-dec! vy i y-adj))))

(define (initialize-particles)
  "Initialize particle positions, scale velocities to target temperature,
and remove systematic drift."
  (for (i 0 num-particles)
       (vector-set! x i 0.0)
       (vector-set! y i 0.0)
       (vector-set! vx i (- (gaussian-deviate-marsaglia) 0.5))
       (vector-set! vy i (- (gaussian-deviate-marsaglia) 0.5))
       (vector-set! ax i 0.0)
       (vector-set! ay i 0.0))
  (place-particles)
  (rescale-velocities)
  (remove-drift))

;;;
;;; Printing Stuff
;;;

;; Display a series of arguments followed by a newline.
(define (println . args)
  (for-each display args)
  (newline))

(define (print-config)
  (printf "\n Simulation Configuration:~%")
  (printf "    Number of particles        : ~a~%" num-particles)
  (printf "    Cutoff radius              : ~a~%" force-cutoff)
  (printf "    Target temperature         : ~a~%" target-temperature)
  (printf "    Integration step size      : ~a~%" time-step)
  (printf "    Width of simulation area   : ~a~%" box-width)
  (printf "    Height of simulation area  : ~a~%" box-height)
  (printf "    Area of simulation         : ~a~%" (* box-width box-height))
  (printf "    Wall stiffness             : ~a~%" wall-stiffness)
  (printf "    Equilibration time         : ~a~%" equilibration-time)
  (printf "       Steps between rescaling : ~a~%" steps-between-equil-rescaling)
  (printf "    Production time            : ~a~%" production-time)
  (printf "       Steps between rescaling : ~a~%" steps-between-prod-rescaling)
  (printf "    Steps per printout         : ~a~%" steps-per-printout))

(define (print-properties-header)
  (printf "Time      Temp.   Pressure Tot. E.   Kin. E.   Pot. E.   Steps~%"))

(define (print-properties)
  (printf "~7,3f,  ~5,3f,  ~6,4f,  ~7,2f,  ~6,2f,  ~8,2f,  ~a~%"
          t average-T average-P (+ k-e p-e) k-e p-e steps-accomplished))

(define (print-random-numbers)
  (println "\nRandom numbers from a uniform distribution.")
  (for (n 0 10)
       (println "Uniform #" n ": " (algo-647-uniform)))
  (println "\nRandom numbers from a normal distribution.")
  (for (n 0 10)
       (println "Normal #" n ": " (gaussian-deviate-marsaglia))))

;;;
;;; Random number stuff
;;;

;; The algo_647_uniform function returns a random uniform deviate in
;; the range [0, 1) calculated using the ACM algorithm 647 as
;; recorded at: http://calgo.acm.org.
;; The algorithm is from the paper "ALGORITHM 647: Implementation and
;; relative Efficiency of Quasirandom Sequence Generators" by Bennet L. Fox,
;; ACM Transactions on Mathematical Software (TOMS), Volume 12, No 4,
;; Dec 1986, pp. 362-376.
;;
;; I have never seen the actual paper - it's behind a paywall - but the
;; FORTRAN version of the algorithm is easy enough to understand.

;; Return a random uniform deviate in the range 0 <= x < 1.0.
(define algo-647-uniform
  (let ([seed 12345])
    (lambda ()
      (let* ([k (floor (/ seed 127773))]
             [partial (- (* 16807 (- seed (* k 127773))) (* k 2836))])
        (if (< partial 0)
            (set! seed (+ partial 2147483647))
            (set! seed partial))
        (* seed 4.656612875e-10)))))

;; The following is an implementation of the Marsaglia polar
;; method for generating random standard normal deviates. See
;; https://en.wikipedia.org/wiki/Marsaglia_polar_method.

;; Return a random standard normal deviate.
(define gaussian-deviate-marsaglia
  (let* ([has-spare #f]
         [spare 0.0]
         [rsq 0.0] [r1 0.0] [r2 0.0]
         [nxt-rsq (lambda ()
                    (set! r1 (- (* 2.0 (algo-647-uniform)) 1.0))
                    (set! r2 (- (* 2.0 (algo-647-uniform)) 1.0))
                    (set! rsq (+ (* r1 r1) (* r2 r2)))
                    rsq)])
    (lambda ()
      (if has-spare
          (begin
            (set! has-spare #f)
            spare)
          (let loop ([_ (nxt-rsq)])
            (if (or (>= rsq 1.0) (= rsq 0.0))
                (loop (nxt-rsq))
                (let ([fac (sqrt (/ (* -2.0 (log rsq)) rsq))])
                  (set! spare (* r1 fac))
                  (set! has-spare #t)
                  (* r2 fac))))))))

(main)
