#!/bin/bash

# Copyright (C) 2013 Andy Aschwanden

set -e # exit on error
SPAWNSCRIPT=spawnparam.sh

 SHEBANGLINE="#!/bin/bash"
MPIQUEUELINE="#PBS -q shared"
 MPITIMELINE="#PBS -l walltime=06:00:00"
 MPISIZELINE="#PBS -l nodes=1:ppn=1"
  MPIOUTLINE="#PBS -j oe"
SOURCEFILE="source /center/w/aschwand/FEniCS/share/fenics/fenics.conf"

for gamma in 1 5 10 50 100
do

  for alpha in 0
  do
#      for project in "jakobshavn" "79N" "helheim"
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

      echo "python scripts/mcb.py --gamma $gamma --alpha $alpha --project $project" >> $SCRIPT

      echo "($SPAWNSCRIPT)  $SCRIPT written"

      SCRIPT="do_${alpha}_${gamma}_${project}_dhdt.sh"
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

      echo "python scripts/mcb.py --dhdt --gamma $gamma --alpha $alpha --project $project" >> $SCRIPT

      echo "($SPAWNSCRIPT)  $SCRIPT written"

      SCRIPT="do_${alpha}_${gamma}_${project}_bmelt.sh"
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

      echo "python scripts/mcb.py --bmelt --gamma $gamma --alpha $alpha --project $project" >> $SCRIPT

      echo "($SPAWNSCRIPT)  $SCRIPT written"

      SCRIPT="do_${alpha}_${gamma}_${project}_all.sh"
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

      echo "python scripts/mcb.py --bmelt --dhdt --gamma $gamma --alpha $alpha --project $project" >> $SCRIPT

      echo "($SPAWNSCRIPT)  $SCRIPT written"
      done
  done
done



echo
echo "($SPAWNSCRIPT)  use submitparam.sh to submit the scripts"
echo

