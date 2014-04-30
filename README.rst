bedmachine
==========

Operation IceBridge bedmachine. Or whatever. Just a toy for now.

Uses continuity equation to interpolate ice thickness in a mass conserving way.
Model implementation in FEniCS/dolfin (http://fenicsproject.org/) by Jesse Johnson, University of Montana.

To get the data, run

``$ sh prepare_jakobshavn.sh N``

``$ sh prepare_79n.sh N``

and then run the model:

``$ mpirun -np N python scripts/mcb.py``

where N is the number of cores

To connect to the database:

$ ssh -L5433:icebridge.sr.unh.edu:5432 bmachuser@icebridge.sr.unh.edu

If you want to use the database via QGIS: 5432 is the database port which you is being mapped to 5433 locally. So, on your local machine, set up the connection as:

Name: Icebridge
Host: localhost
Port: 5433
Database: icedb
Username: nobody
