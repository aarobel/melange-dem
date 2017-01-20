README for TerminusAdvance subfolder

In this folder are scripts which specify and run a LIGGGHTS simulation for
melange under the forcing of a gradually advancing
terminus, with sea ice bonds that dynamically break and  are reformed at the 
beginning of each sub-run (the number and length of which are specified 
in the frs*.sh run scripts). The structure of these scripts (specifically,
breaking the full runs into sub-runs that coincide with the sea ice
reformation time step, rather than using the bond create functionality of the
DESign toolbox) is wholly a consequence of the significant computational
expense of these runs and the runtime limitations of the cluster on which they were
run. (For reference a 5 day run with 1600 elements and time step of 0.002
seconds, take approximately 12 hours on a 16-core node.)

All scripts have four different varieties corresponding to the four lines in
Figure 2 of Robel (2017). NLF scripts suppress and landfast sea ice formation,
GN5 scripts including a narrowing in the fjord geometry.

frs*.sh are the driver scripts, which specify certain parameters for the runs
and then calls makebatch*.sh to create a batch script. When this batch script
is executed, it queues a sub-run on the cluster which each successive run
being dependent on the completion of the previous run. All of the syntax used
in this driver is based on a cluster with SLURM scheduling.

makebatch*.sh scripts control each sub-run. They create an initial condition
file from the final output of the previous sub-run (using
chanbonds_fromstop*.m scripts), during which sea ice bonds are re-frozen based
on proximity between elements. The makeLAMMPS*.m scripts are also called which
create a run script based on parameter specifications in the frs*.sh driver
and then this run script is run with LIGGGHTS and DESign.

makeLAMMPS*.m scripts control the LIGGGHTS run script. Any parameter values
which are not specified in the driver file (frs*.sh) are specified here
(apologies). 

Other scripts, such as assembleresults_fix.m and readdump.m are based on similar
utilities in the DESIgn toolbox (Herman 2016).
