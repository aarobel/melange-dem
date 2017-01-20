README for SpringCalving subfolder

In this folder are scripts which specify and run a LIGGGHTS simulation for
melange initialized from a TerminusAdvance run, including sea ice bonds. The
in.springcalve script specifies a LIGGGHTS simulation with the DESign toolbox
in which a small calving event occurs and then allows for relaxation. This
includes a high temporal resolution model output to capture the rapid wave
propagation in the melange in response to the calving event.

Preferably this script is used in conjunction with TerminusAdvance to produce
an initial condition, however the initial condition used to produce Figure 4
in Robel (2017) is provided as a starting point (it comes from the 2 meter sea
ice GN5 simulation in TerminusAdvance).
