# Resources
This folder contains all the input data needed to run the simulator

## missileDATCOM data
This folder must cointain missileDATCOM data for empty and full states.
Saved as "for006_empty.mat" and "for006_full.mat". 
In order to obtain those data from missileDATCOM output file you need to use the parser. 

## "config.m" 
The config file is where you set up all the relevant paramenters about the simulation, such as geometry, engine and wind options.

**STOCHASTIC SIMULATION**
The stochastic simulation only returns a plot of landing points and a histogram of the apogee altitudes distribution.

## "start_simulation.m"
Run this file to start the simulation.