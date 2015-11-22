# skyward-matlab-rocket-simulator
This is a script developed in MATLAB for the simulation of rocket 6 d.o.f. dynamics.

## How to use it

Just clone the repository. You're going to find two folders, all the *core files* are inside the `sources` folder. 
In the `resources` you're going to find a `config.m` sample file and some aerodynamics dataset of the Skyward Experimental Rocketry's first sounding rocket <a href="http://www.skywarder.eu/blog/rocksanne-i-x/">Rocksanne I-X</a> in the early stages of design (so they can be wrong).

- Fill a `config.m` file and put it inside the `sources` folder.
- Provide the correctly formatted aerodynamics dataset files (see `resources/for006*` files for reference). We used Missile DATCOM '97 software to estimate all the aerodynamics coefficients and <a href="https://github.com/skyward-er/skyward-datcom-parser"> a parser </a> to convert the output text file into MATLAB dataset.
- From the command window run `MAIN`

## `config.m` file
The config file is where you set up all the relevant data about the rocket. So just copy and paste the one inside the `resources`folder and replace the data you need to. The file is fully commented so the variables should be very straightforward.

**PARACHUTES DETAILS**

The script is configured to simulate the descent with a 3 d.o.f. model and two parachutes, a drogue (secondary) parachute and a main parachute. 

**WIND DETAILS**

The wind model is provided by the `windgen.m` function. The function just (buggy) implements a random constant in magnitude wind. It has shown to have few bugs if the starting angle is not zero, so it has to be improved. 

**STOCHASTIC DETAILS**

If the variable `settings.stoch.N` is set to a number greater than one, a *stochastic simulation* is fired, that means running the simulation several times calling the `windgen.m` function each time, so providing different wind magnitudes and orientations. 
The `settings.stoch.parallel` variable implement the `parfor` instead of the classical `for` for the *stochastic simulation* leading to a sensible drop of computation time due to multitasking.

#License
see `sources/LICENSE/LICENSE` file.
