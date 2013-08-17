#!/bin/bash

# Copyright (C) 2013 Andy Aschwanden

set -e # exit on error
source /Users/andy/FEniCS/share/dolfin/dolfin.conf 

scale=1.0  # default number of processors
if [ $# -gt 0 ] ; then
  scale="$1"
fi

for gamma in 0.0 0.5 1.0 2.0 5.0 10.0 20.0 50.0 100.0 200.0 500.0 1000.0 2000.0 5000.0 10000.0
do
    for alpha in 0 
    do
	for project in "jakobshavn"
	do
	    python scripts/mcb.py --grid_spacing 250 --velocity_scale $scale --gamma $gamma --alpha $alpha --project $project
	    python scripts/mcb.py --grid_spacing 250 --velocity_scale $scale --gamma $gamma --alpha $alpha --project $project --dhdt
	done
    done
done
