README for public GitHub repository corresponding to melange model in Robel (2017), Nature Communications
Please contact Alex Robel (robel@caltech.edu or alexander.robel@gmail.com) for questions or bugs

The folders in this repository contain code that is used along with MATLAB (for preparing initial
condition) and the LAMMPS-LIGGGHTS discrete element model paired with the DESign sea ice toolbox.
Instructions for proper installation of LIGGGHTS together with the DESIgn toolbox can be found on
Agnieszka Herman's website (http://herman.ocean.ug.edu.pl/LIGGGHTSseaice.html) and the implementation
is detailed in Herman (2016) in Geoscience Model Development.

Using the most recent version of DESIgn on Intel processors, one should be able to successfully run
the scripts in this repository. 

The ConvergeInitialConditions folder includes scrips for making initial
conditions that can then be used as a starting point for all the simulations
in other folders. This is heavily based on utility scripts from DESIgn and
will create unbonded melange with specified size distributions compacted into
a rectangular domain. The scripts used to make initial conditions in the other
folders then start from this initial condition and create melange confined
within a channel bounded by static elements and moving terminus elements.

The SummerCalving and SpringCalving folders include scripts for
running short simulations with prescribed calving events (though which different initial conditions
and different size/duration of calving events). These correspond to Figures 3 and 4 in Robel (2017).
They are short (~hours) and can comfortably be run on a few or just a single processor.

The TerminusAdvance folders includes scripts to run long simulations with
steady terminus advance throughout the runs. These runs are split into 5-day
sub-runs to be able to run on clusters with limitations on job lengths. For
reasonable run times, it is recommended to run each job on 8 to 16 cores (at
least). Potential issues may arise due to bond breaking between cores (causing
the simulation to freeze). This can generally be alleviated by using fewer
cores or shortening the length of sub-runs (though these will require more run
time).

