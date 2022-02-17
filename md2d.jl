"""
# module crng

- Julia version: 1.6.2
- Author: david
- Date: 2021-07-23

# References

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

# Examples

```jldoctest
julia>
```
"""
###module md2d

###module crng

using Printf

function main()
    @printf "md2d - MD simulation of a 2D argon gas Lennard-Jones system.\n"
    @printf "This version is written in Julia. Version: %s.\n" VERSION

    @printf "\nEquilibration\n"

    @printf "\nProduction\n"

    print_random_numbers()
end

#
# Printing stuff.
#

function print_random_numbers()
    @printf "\nRandom numbers from a uniform distribution.\n"
    for i in 0:9
        @printf "Uniform #%i: %.17f\n" i algo_647_uniform()
    end
    @printf "\nRandom numbers from a normal distribution.\n"
    for i in 0:9
        @printf "Normal #%i: %.17f\n" i gaussian_deviate_marsaglia()
    end
end

## Random number stuff.
##
## The algo_647_uniform function returns a random uniform deviate in
## the range [0, 1) calculated using the ACM algorithm 647 as
## recorded at: http://calgo.acm.org.
## The algorithm is from the paper "ALGORITHM 647: Implementation and
## relative Efficiency of Quasirandom Sequence Generators" by Bennet L. Fox,
## ACM Transactions on Mathematical Software (TOMS), Volume 12, No 4,
## Dec 1986, pp. 362-376.
##
## I have never seen the actual paper - it's behind a paywall - but the
## FORTRAN version of the algorithm is easy enough to understand.

# Something to hold an mutable seed value for the random number generator.
mutable struct seedType
    seed::Int32
end

# # Something to hold an mutable seed value for the random number generator.
glbl_seed = seedType(12345)

# Allow the seed to be set externally.
function algo_647_set_seed(new_seed)
    setproperty!(glbl_seed, :seed, new_seed)
end

# Return a random uniform deviate in the range 0 <= x < 1.0 as a Float64.
function algo_647_uniform()
    k = floor(glbl_seed.seed/127773)
    partial = 16807 * (glbl_seed.seed - k * 127773) - k * 2386
    if (partial < 0)
        setproperty!(glbl_seed, :seed, partial + typemax(Int32))
    else
        setproperty!(glbl_seed, :seed, partial)
    end
    return glbl_seed.seed * 4.656612875e-10
end

## The following is an implementation of the Marsaglia polar
## method for generating random standard normal deviates. See
## https://en.wikipedia.org/wiki/Marsaglia_polar_method.

mutable struct spareType
    has_spare::Bool
    spare::Float64
end

glbl_spare = spareType(false, 0.0)

# Return a random standard normal deviate as a double.
function gaussian_deviate_marsaglia()
 # static int has_spare = 0;
 # static double spare;
    fac = 0.0
    rsq = 0.0
    r1 = 0.0
    r2 = 0.0
    if (glbl_spare.has_spare)
        setproperty!(glbl_spare, :has_spare, false)
        return glbl_spare.spare
    else
       while true
           r1 = 2.0*algo_647_uniform() - 1.0
           r2 = 2.0*algo_647_uniform() - 1.0
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
###end of module
