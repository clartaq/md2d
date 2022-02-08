# md2d

## Description

The source files in this directory contain programs to run just about the
simplest molecular dynamics simulation possible that generates meaningful
property calculations.
The programs simulate a two dimensional argon gas where the gas atoms have
a Lennard-Jones potential.

The source files illustrate how to do the same simulation using four
different languages: C, Java, Rust, and Julia.

As mentioned, theses programs are about as simple as possible. There
is no fancy graphical output, just print statements to display data to the
terminal.

### What's the Deal with the Random Number Generators (RNGs)

You may ask yourself why a significant fraction of the code is devoted to
creating RNGs. After all, each of the languages
have standard libraries with facilities to generate random numbers, probably
of better quality than those used here.

The intent is for the programs to produce results that are as close as possible
to the programs in the other languages. To do that, where we use random
numbers, we want them to be identical among the different languages. With
the versions provided here, we can generate the same sequence of random 
numbers in each program, provided we start with the same seed value.

## Tools

- Java and JavaFX, obviously. I used OpenJDK Java version 17.0.2 and 
OpenJFX JavaFX 17.0.2.

- [Maven](https://maven.apache.org/index.html) 3.8.4 was used to obtain dependencies, build the program,
and run it from the command line.

- Other development (trivial refactoring) was done using 
[IntelliJ IDEA](https://www.jetbrains.com/idea/) 2021.3.

## Installing and Running the Project

Copy the repository files to a convenient directory.

To run the project, in the project root directory, execute the following
command line:

    mvn clean javafx:run

The repository contains project files for use with IntelliJ IDEA and can be
run directly from the IDE if so desired.

## Credits


## License

This code is licensed with the [MIT license](https://opensource.org/licenses/MIT).