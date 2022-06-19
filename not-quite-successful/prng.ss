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
    (lambda new-seed
      (if (pair? new-seed)
          (set! seed (car new-seed))
          (let* ([k (floor (/ seed 127773))]
                 [partial (- (* 16807 (- seed (* k 127773))) (* k 2836))])
            (if (< partial 0)
                (set! seed (+ partial 2147483647))
                (set! seed partial))
            (* seed 4.656612875e-10))))))

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
