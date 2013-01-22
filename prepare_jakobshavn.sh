#!/bin/bash

# You need to open a port first in a different shell:
# ssh -L5433:icebridge.sr.unh.edu:5432 bmachuser@icebridge.sr.unh.edu

# Kinda non-standard tools needed:
# CDO (https://code.zmaw.de/projects/cdo) e.g. from MacPorts
# Pyresample (https://code.google.com/p/pyresample/)
# nc2cdo.py (part of the PISM distribution)

set -x -e

# run ./prepare.sh 1 if you havent CDO compiled with OpenMP
NN=8  # default number of processors
if [ $# -gt 0 ] ; then
  NN="$1"
fi

# grid spacing
GS=500
if [ $# -gt 1 ] ; then
  GS="$2"
fi


# corners of lon,lat box for query
LON_MIN=-51
LON_MAX=-42
LAT_MIN=68.25
LAT_MAX=70

MOD_VAL=1
MOD_FIELD='gid'

EPSG=3413

PROJECT=jakobshavn
PROJECTC=Jakobshavn
YEARA=2006
YEARE=2012
CRESIS_YEARS=${YEARA}_${YEARE}

# destination area in EPSG:3413 coordinates
X_MIN=-230000.0
X_MAX=80000.0
Y_MIN=-2350000.0
Y_MAX=-2200000.0

source prepare_main.sh