# skyward-matlab-rocket-simulator
This is a code developed in MATLAB for the simulation of 6 d.o.f. rocket  dynamics. <br/>
The simulator predicts 3D trajectory, apogee, forces acting on the rockets, and various other aerodyanmics data. 

## Usage
The simulator requires aerodynamic coefficients computed using Missile DACTOM 97: one with the parameters of
the rocket full of fuel and another with the parameters of the empty rocket. These must be stored in the `\data` folder with the name `rocketname_empty.mat` and `rocketname_full.mat`.
Before starting the simulation, data and parameters to be used must be specified in `config.m` by changing the appropriate field in the `settings` variable, such as
* Rocket name
* Launchpad data (height, inclination)
* Geometry, mass, aerodyanmic properties (length, mass, moments of inertia, reference surface etc..)
* Motor data
* Wind model to be used (constant wind or atmoshwm wind model)
* Integration options
* Type of simulation to be performed (ascend phase + parachute descend, ballistic, stochastic, parachute failure etc..) 
* Plot settings

The appropriate simulation is then started running `start_simulation.m`,based on the settings specified

## Types of simulation
The standard simulation is a 6 d.o.f. dynamics in the ascend phase and 3 d.o.f. descent with a parachute. However it is possible to simulate the descend phase using a 6 d.o.f. dynamics for a ballistic flight or in the case of the parachute failure. <br/>
It also allows the possibility of a stochastic simulation using different values of wind intensity and direction, in order to perform a stastistics analysis of the landing position

## Main output
The simulator records the following data: apogee height, time of flight, maximum acceleration, velocity and Mach reached, run in the launchpad. <br/>
The main plots are: stability margin, altitude, Mach number, forces acting on the rocket, aerodynamic angles, CD, temperature on the nose, body frame velocities, accelerations, eulerian angles and angular rates vs time. <br/> It also plots the 3D trajectory and its projections in the xz,xy,yz planes

More detailed information about the physical model, how to use and how to modify the simulator is available in simulator-documentation.pdf
