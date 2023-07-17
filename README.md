# Hyades-Planar-Scripts
Here I have some basic scripts that I have used for running planar 1D Hyades simulations of a laser irradiating a sandwhich style target. These are much less sophisticated than my code for the spherical sims (Hyades-ICF-Scripts repo), and generally require the generator scripts to be adapted for each new problem.

File Generation:

_The meshing code here is specifically for a set of simulations I performed, where the target has 4 layers of fixed thickness. It is possible to adapt this code for different thicknesses/layers. Details on how to do this are included in the scripts_

**PlanarMeshing_Git**
This script allows you to mesh the problem. 

**PlanarFileWriter_Git**
Once meshed (with all the variables still open), this script allows you to produce an input deck for the problem.


Analysis:

**PlanarReaderBasic_Git**
A basic file reader for planar hyades sims.

**PlanarReaderWithTracking_Git**
A more complicated file reader which attempts to do shock tracking through the rear two layers. More specialised to my four layer target simulations (but may hopefully come in handy if you're trying to do something similar - even if it needs some adapting!).
