# md2d

**md2d** is a 2D molecular dynamics simulation. The simulation is of argon-like
particles interacting _via_ the Lennard-Jones potential.

The calculations are illustrated
by writing versions of the simulation in six different languages.

## Description

The source files in this directory contain programs to run about the
simplest molecular dynamics simulation possible that generates meaningful
property calculations.
The programs simulate a two dimensional argon gas where the gas atoms have
a Lennard-Jones potential.

The source files illustrate how to do the same simulation using six
different languages: C, Java, Julia, Racket, Scheme, and Zig. (There are 
versions of the program written in the Janet and (Free) Pascal languages,
but they are unfinished and do not produce the same results as the other
programs.)

C and Java are
"traditional" languages that have a history of use for molecular dynamics
programming. Molecular dynamics in C started as an attempt to write
MD code in a more modern language than FORTRAN.

Zig and Julia are newer languages that look like good
candidates for these types of calculation-intensive programs.

Racket and Scheme are members of the Lisp language family. Lisps are not often considered for the types of
calculations used in molecular dynamics, but, with contortion, can be
molded for use.

## What Do These Programs Do?

The programs all calculate the trajectories of 500 disks (2D soft particles) over a period of 400 time steps. The particles interact with each other via a Lennard-Jones potential.

At the beginning of the simulation, the particles are placed in the simulation space in a rectangular array with a little random "jitter" in their locations. The particles are given random velocities for which the entire group has a given mean temperature.

As mentioned, theses programs are about as simple as possible. There
is no fancy graphical output, just print statements to display data to the terminal.

Periodically, the temperature, pressure, total energy, kinetic energy, and potential energy of the system are all calculated and displayed on the terminal. The simulation goes through an "equilibrium" phase where the temperature of the system is reset every few steps. That phase is followed by a "production" phase in which the system evolves with less frequent temperature rescaling.

## Installation

Copy or clone the repository into a convenient directory. Then`cd` into
that directory to begin work with the repository.

## Tools and Development Environment

All development was done on an Intel iMac Pro running macOS Monterey version 12.
It should work similarly on Linux.
Compiling and running on Windows is left as an exercise for the user.

All program editing was done in [Emacs](https://www.gnu.org/software/emacs/)
27.2 using plugins for the major modes for each language.

### C

The C version of the program was written to be compatible with the
[201710 standard](https://en.wikipedia.org/wiki/C17_(C_standard_revision)) for
the language.

This is the version of the program that was used as the prototype when
writing new programs in the other languages.

A small shell script, `xmd2dc.sh` is included to compile and run the program.
The contents
of the shell script are `rm md2dc; gcc -O3 md2d.c -o md2dc; ./md2dc`.
 Obviously, [GCC](https://gcc.gnu.org) is used to compile the program. This is
 the only version of
 the program that is compiled with explicit optimizations.

### Java

The following tools were used during development:

- Java, obviously. I used
[OpenJDK Java version 17.0.2](https://jdk.java.net/17/).

- [Maven](https://maven.apache.org/index.html) 3.8.4 was used to obtain
dependencies and build the program.

- Other development (trivial refactoring) was done using
[IntelliJ IDEA](https://www.jetbrains.com/idea/) 2021.3.

To run the program, simply type `java md2d.java` at the command line.
Alternatively, there is a shell script, `xmd2dja.sh` that can be used to
compile and run the program. The contents of the shell script are
`javac md2d.java; java md2d`.

### Julia

[Julia](https://julialang.org) is a relatively new language intended to
write high performance software.

To run the Julia program, type `julia md2d.jl`.

### Racket

[Racket](https://racket-lang.org) 8.6 was used for development in that language. The Racket program
is a slightly modified version of the Scheme version. DrRacket, a graphical 
IDE, was used for modification, testing, and debugging for the Racket version 
of the program.

The program, `md2dr.rkt` can be loaded into DrRacket and run directly or the 
IDE can be used to create a standalone executable, which is much faster that 
running the program in the IDE.

To run the standalone executable, type `md2dr` at the command line.

Because I like Scheme so much, there is another version written for Racket
that follows the [R7RS](https://small.r7rs.org/attachment/r7rs.pdf) standard for the language.

This is accomplished by using an additional library, [racket-r7rs](https://github.com/lexi-lambda/racket-r7rs), to
add the necessary features to Racket. Follow the instructions at the repository
to install the library.

After installing the library, the R7RS version can be run from the command
line with:

```
racket -I r7rs -f md2d7.scm
```

(I also examined using several different native R7RS implementations but was
disappointed in their performance, level of compliance, or something else.)

### Scheme

[Chez Scheme](https://scheme.com) 9.5.6 was used to develop the Scheme version of the program.
This version follows the R6RS version of the standard.

The Scheme version can be run by typing `chez --program md2d.ss` at
the command line.

### Zig

[Zig](https://ziglang.org) is a young language often described as "A Modern
C" or "A Better C than C". Version 0.9.1 was used
to write the program in the repository.

To run the program, type `zig build-exe md2d.zig --name md2dz; ./md2dz` at
the command line.

## Some Not Entirely Successful Versions

The directory `not-quite-successful` contains implementations in two
additional languages: Janet and Pascal.

All versions of the program run to completion. The Janet and Pascal
versions produce results similar to the other versions, just not
exactly the same. In fact, the results from these two versions are
close to each other. I assume there is some numeric configuration
that I have not figured out that will bring them into agreement
with the other programs.

### Janet

[Janet](https://janet-lang.org) is a young and interesting version of Lisp. It takes many ideas for
a modern Lisp from [Clojure](https://clojure.org).
A version of the program written in this language
is included in the repository, but it is incomplete and does not match
the output of the other programs.

Janet version 1.20.1 was used to write a version of the program in that
language.

The Janet version can be run by typing `janet -l ./md2d` at the command line.

### Pascal

The Pascal version is written in the [Free Pascal](https://www.freepascal.org)
variant of the language, which is similar to
[Delphi](https://www.embarcadero.com/products/Delphi). Free Pascal has a nice
IDE associated with it called [Lazarus](https://www.lazarus-ide.org). I just
can't run the program from the IDE for some reason. It compiles fine but
the output does not show up on any console.

The Pascal version can be compiled to an executable by typing
`fpc md2dp.lpr` at the command line. Then, run it like normal with `./md2dp`.

## Credits


## License

This code is licensed with the
[MIT license](https://opensource.org/licenses/MIT).
