#!/bin/bash

# Copyright (C) 2013 Andy Aschwanden

set -e # exit on error
source /Users/andy/FEniCS/share/dolfin/dolfin.conf 

scale=1.0
if [ $# -gt 0 ] ; then
  scale="$1"
fi

for gamma in 0.0 0.5 1.0 2.0 5.0 10.0 20.0 50.0 100.0
do
    for alpha in 0 1000
    do
	for project in "jakobshavn"
	do
	    python ~/base/bedmachine/scripts/mcb.py --grid_spacing 250 --velocity_scale $scale --gamma $gamma --alpha $alpha --project $project
	    python ~/base/bedmachine/scripts/mcb.py --grid_spacing 250 --velocity_scale $scale --gamma $gamma --alpha $alpha --project $project --dhdt --bmelt
	done
    done
done
