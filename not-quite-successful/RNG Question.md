---
author: clartaq
title: Pseudo-Random Numbers in Scheme
date: 2022-10-24T08:57:33.000-04:00
modified: 2022-05-01T16:24:47.417-05:00
tags:
  - Scheme
  - Random Numbers
  - Functional Programming
  - Pseudo-Random Numbers
  - PRNG
---
For reasons that are too dreary to detail, I need to reproduce several pseudo
random number generators in different languages. It mostly involves translating
ancient C and FORTRAN code. No problems with C, Java, Zig, Julia, etc. But for
Lisp-derived languages, I would like to do something more idiomatic.

For example, here is a Scheme interpretation of ACM algorithm 647 for a
portable PRNG to generate samples from a uniform distribution.

```scheme
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
```

The above works and I'm satisfied with the implementation (though eager to hear
suggestions for improvements).

The implementation below of a function to generate normally distributed PRNGs
also works, but looks kind of ugly -- there is so much state to maintain.
Is there a better way to do the following?

```scheme
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
```

In Chez 9.5.6, the first 10 results should be:

```shell
0.4649329739402326
1.4550052699768505
-1.0580380669115383
-0.5790254729247644
0.8434589541668004
-0.6708443571574382
-0.22644041228981196
-0.07818079860601053
-0.7443285279492631
-0.5232388154010481
```


