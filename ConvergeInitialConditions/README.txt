README for ConvergeInitialConditions subfolder

In this folder are example scripts which produce a random initial
condition of unbonded elements and then compresses those elements
into a domain of specified size.

batch.sh is a driver file that runs the MATLAB initialization script
and runs LIGGGHTS with the DESIgn toolbox.

LIGGGHTS_makeinit.m is a MATLAB script which specifies the desired size
distribution of elements and then calls LIGGGHTSinit.m, a modifid version
of a utility from the DESign toolbox that creates an initial condition
file that can be the starting point for a LIGGGHTS simulation. In this
initial condition, the elements have random radii (drawn from the specified
size distribution) and then placed randomly in a very large domain such that
they do not overlap.

in.conv is a run script for LIGGGHTS with the DESIgn toolbox which starts from
an initial condition and then moves the sides of the domain (which is
periodic) inwards to a new specified domain size. This will compress the
sparse initial condition to a more well packed state. After compression, the
simulation continues to run for a period of time to allow for rearrangement
and relaxation of the elements to a low-stress state.
