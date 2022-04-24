# md2d Walkthrough

The following is a short description of how the programs are structured. All
versions have a similar structure.

Because the languages have different conventions for naming variables, like
"kebab" case, "camel" case, and "snake" case, the variables and 
function/method/procedure names are not identical across program versions.
The convetions used here are those of the Java version. You are smart enough
to translate easily to the relevant names in the other languages.

## Style

The program is written in a strongly imperative style. This style derives
almost directly from the first molecular dynamics programs, which were
written in FORTRAN.

There tends to be a lot of global variables and mutation is everywhere.

The "particles" themselves are represented as a combination of array elements.
Not records, not objects, not structures. The "ith" elements of the arrays
represent the x position, y position, x velocity, y velocity, acceleration
along the x axis, and acceleration along the y axis.

"Bookkeeping" variables like the simulation time and step count are held
in global scalars, as are calculated properties like the temperature,
pressure, and energies.

## Structure

### Imports

Each language requires that certain features of its run-time library must
be imported into the program. There are between one and three such imports
following the licensing comment at the top of the file.

### Simulation Parameters

Each version of the program starts with a group of variables that control the 
characteristics of the simulation. These values are held constant during the 
execution of a program.

**NUM_PARTICLES**: The number of particles in the simulation. This variable
determines the size of arrays and length of many program loops.

**FORCE_CUTOFF**: The distance between particles beyond which the interaction
forces are not calculated.

**TARGET_TEMPERATURE**: The temperature, in natural units, that the simulation
attempts to maintain. This particular temperature was chosen because it is
near the transition from a "gas" to a "liquid" for these conditions. At this
temperature, the particles tend to clump into one or two "droplets"
surrounded the other particles in a gas-like state.

**TIME_STEP**: The integration time step. This value was chosen so as to be
big enough for a quick runs without having to worry about runaway
instability.

**BOX_WIDTH** and **BOX_HEIGHT**: The size of the simulation box. For the
moment, these two values should be the same.

**WALL_STIFFNESS**: The "spring constant" representing the springiness of
the walls when particles hit them.

**EQUILIBRATION_TIME**: The simulation  time during
which the conditions are evolving from the starting conditions to an
equilibrium.

**STEP_BETWEEN_EQUIL_RESCALING**: The number of simulation steps that are
executed before the simulation temperature is re-scaled back to the
target temperature.

**PRODUCTION_TIME**: The simulation time after equilibration when generated
values should be stable.

**STEPS_BETWEEN_PROD_RESCSALING**: The number of simulation steps during
the production phase that are executed before conditions are re-scaled
back to the target temperature.

**STEPS_PER_PRINTOUT**: The number of simulation steps executed between
printing out status messages.

### Global Variables

The `x`, `y`, `vx`, `vy`, `ax`, and `ay` vectors contain the state of the
particles -- their positions, velocities, and accelerations.

The `kE` and `pE` (or similarly named) variables contain the values for the
kinetic energy and potential energy, respectively.

The temperature and pressure are calculated from running averages of
instantaneous measurements. `average_t` and `average-p` hold the average
values for pressure and temperature that are displayed as the simulation
progress.

The simulation time is tracked in the variable `t`. The cound of simulation
steps is held in the variable `stepsAccomplished`. The `sampleCount` variable
holds the number of measurements used to calculate the average temperature
and pressure.

### "Forward" Declarations

Some languages have no rules about the ordering of functions in the program
text. Others require that any external functions used must appear before
that use. For those languages the require declaration before use, there
is a section of forward declarations following the global variables.

### Functions/Methods/Procedures

**`main`**: This is the entry point for the program. It begins by printing
out some program information, such as the language that the particular
version was written in and the simulation parameters.

Initializing the particles consists of two main tasks: 1) placing the particles
within the simulation box and 2) giving the particles an initial velocity
based on the desired temperature. This is accomplished by:

**`initializeParticles`**: The method starts by assigning temporary values
to the particle positions and velocities. The accelerations are initialized
to zero at the beginning of the simulation. The velocities are assigned
random values that will be scaled to the correct temperature later by the
**`rescaleVelocities`** method.

As part of the initialization, the particles are placed in a gride near the 
bottom of the simulation box by the **`placeParticles`** method. This method
also applies a small amount of random "jitter" to the particle positions in
order to break up the perfect arrangement of particles.

Once the particles are placed and given initial velocities, any systematic
drift in the particles as a group is removed by a call to the **`removeDrift`**
method.

After initializing the particles, the function runs two loops: one for the
equilibration phase and another for the production phase. The only difference
between these two phases is the frequency of temperature re-scaling.

As the loops run, they periodically display information about the progress
of the simulation and some properties calculated from the state of the
simulation.

Several version of the program keep track of how long it takes to complete
the program and display that information just before the program ends.

(These measurements should not be considered benchmarks. They include the
time taken to display progress messages and so on. They are not an accurate
reflection of the speed of calculations, just sort of general guides.)

**`singleStep`**: As the name implies, this function is responsible for
running a single step of the simulation. It uses the Verlet algorithm to
update the particle positions, velocities, and accelerations.

**`computeAccelerations`**: This is where the heavy work of the simulation
occurs.

## What's the Deal with the Random Number Generators?

You may ask yourself why a significant fraction of the code is devoted to
creating RNGs. After all, each of the languages
have standard libraries with facilities to generate random numbers, probably
of better quality than those used here.

The intent is for the programs to produce results that are as close as possible
to the programs in the other languages. To do that, where we use random
numbers, we want them to be identical among the different languages. With
the versions provided here, we can generate the same sequence of random 
numbers in each program, provided we start with the same seed value.

We need two types of random number generators. The first is to produce a 
sequence of random numbers uniformly distributed between 0 and 1.

The second generator calculates random sample from a standard normal deviate.

### The Scheme Version
If you aren't already familiar with Scheme or some other Lisp-like language, 
the random number generators may look a little odd. The `algo-647-uniform` 
definition implements the same algorithm and produces the same results, but is 
put together differently. In the Scheme version it is not really a function; 
it is a symbol whose value is a function.