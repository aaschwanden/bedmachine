bedmachine
==========

Operation IceBridge bedmachine. Or whatever.

To get the data, run

$ sh prepare_jakobshavn.sh N

and then run the model:

$ mpirun -np N python scripts/mcb_jakobshavn.py,

where N is the number of cores