#!/bin/bash

# Copyright (C) 2013 Andy Aschwanden

set -e # exit on error
SPAWNSCRIPT=spawnparam.sh

 SHEBANGLINE="#!/bin/bash"
MPIQUEUELINE="#PBS -q shared"
 MPITIMELINE="#PBS -l walltime=00:20:00"
 MPISIZELINE="#PBS -l nodes=1:ppn=1"
  MPIOUTLINE="#PBS -j oe"
SOURCEFILE="source /center/w/aschwand/FEniCS/share/fenics/fenics.conf"

for gamma in 0.0 1.0 2.0 5.0 10.0 20.0 50.0 100.0 200.0 500.0 1000.0 2000.0 5000.0 10000.0 20000.0 50000.0
do
  for alpha in 0 
  do
      for project in "jakobshavn"
      do
      SCRIPT="do_${alpha}_${gamma}_${project}.sh"
      rm -f $SCRIPT

      # insert preamble
      echo $SHEBANGLINE >> $SCRIPT
      echo >> $SCRIPT # add newline
      echo $MPIQUEUELINE >> $SCRIPT
      echo $MPITIMELINE >> $SCRIPT
      echo $MPISIZELINE >> $SCRIPT
      echo $MPIOUTLINE >> $SCRIPT
      echo >> $SCRIPT # add newline
      echo "cd \$PBS_O_WORKDIR" >> $SCRIPT
      echo $SOURCEFILE >> $SCRIPT

      echo >> $SCRIPT # add newline

      echo "python scripts/mcb.py --grid_spacing 250 --gamma $gamma --alpha $alpha --project $project" >> $SCRIPT

      echo "($SPAWNSCRIPT)  $SCRIPT written"

      done
  done
done



echo
echo "($SPAWNSCRIPT)  use submitparam.sh to submit the scripts"
echo

