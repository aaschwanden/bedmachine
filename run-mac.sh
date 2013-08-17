#!/bin/bash

# Copyright (C) 2013 Andy Aschwanden

set -e # exit on error
source /Users/andy/FEniCS/share/dolfin/dolfin.conf 

for gamma in 0.0 1.0 2.0 5.0 10.0 20.0 50.0 100.0 200.0 500.0 1000.0 2000.0 5000.0 10000.0 20000.0 50000.0
do
    for scale in 0.5 0.75 1.0 1.25
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
done
