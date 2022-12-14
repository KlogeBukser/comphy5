# Computational Physics Project 5: Quantum interference

## DESCRIPTION
- The project simulates a physical system where a quantum state interact with slits in a wall.
- It uses numerical methods to update position and velocity of the quantum state.
- The simulation is written in object oriented C++ code.
- Plotting is written in python3.
- A makefile is used to compile the code, and perform predetermined executions.


## PRODUCE RESULTS
- The code is run from 'src' folder
- make compile: Compiles the c++ program.
- make execute: Makes all relevant results.
- make all: Performs both 'compile' and 'execute'.

- make probability: runs two simulations, one without a middle wall, and one with 2 slits. 
For each time step, the total probability is written to a txt file, along with the time.
These two txt files are read, and the deviation (from 1) of the total probability is plotted as a function of time.
Both graphs are made in the same figure.

- make heat: one simulation is run, with 2 slits. Multiple parameters for time are used.
The model writes the state at each spatial point to file at the times given as parameters. 
These files are then read, and the distribution is plotted as a heatmap in x,y- direction.
The probability density gives the intensity of the color. 
This is also done for the real and imaginary parts of the wave function, seperately.

- make interference: three simulations are run, with 1,2,3 slits in the wall. 
After the simulations is finished, the distribution at x = 0.8 is written to file.
These files are read, and then the density is normalized over the range when x is kept constant.
This is plotted for each of the three simulations separately.

- make animate: This does the same as 'heat', but the state is saved for every timestep.
An animation is made instead of the plots, with the probability density as the color.

## Parameters
- When using the premade run configuration in the makefile, the parameters are chosen for you. 
- The parameters are: instruction, n_slits, h, dt, xc, sig_x, px, yc, sig_y, py, V0, T,
Where instruction is a string that helps the main cpp program decide which type of simulation to run. 
n_slits is the number of slits in the wall. h and dt are the spatial and time-like step sizes respectively.
xc,yc are the x,y coordinates of the center of the state distribution. With px,py as their momentums.
sig_x,sig_y describe the width of the initial distribution. V0 is the strength of the potential in the wall.
T is the total runtime. If more parameter than this are included, they are regarded as addional time parameters.

## CONVENTIONS:
- snake_case (lower case letters separated by an underscore) is used for most variable- and function-names.
- PascalCase (first letter of every word is capitalized, including the very first word) is used for class names.