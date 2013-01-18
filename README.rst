bedmachine
==========

Operation IceBridge bedmachine. Or whatever. Just a toy for now.

Uses continuity equation to interpolate ice thickness in a mass conserving way.
Model implementation in FEniCS/dolfin (http://fenicsproject.org/) by Jesse Johnson, University of Montana.

To get the data, run

``$ sh prepare_jakobshavn.sh N ``
``$ sh prepare_79n.sh N ``

and then run the model:

``$ mpirun -np N python scripts/mcb.py ``

where N is the number of cores
