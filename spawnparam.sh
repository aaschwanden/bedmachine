#!/bin/bash

# Copyright (C) 2009-2011 Ed Bueler and Andy Aschwanden

#  creates 9 scripts each with NN processors and (potentially) submits them
#    on pacman.arsc.edu

#  needs rawparamscript, trimparam.sh

#  usage: to use NN=64 processors, nodes 16, and duration 16:00:00,
#     $ export PISM_WALLTIME=16:00:00
#     $ export PISM_NODES=16
#     $ ./spawnparam.sh 64
#     (assuming you like the resulting scripts)
#     $ ./submitparam.sh      ### <--- REALLY SUBMITS using qsub

#  see submitparam.sh


set -e # exit on error
SPAWNSCRIPT=spawnparam.sh

#####################################################################
# Default scale used by float functions.

float_scale=3

#####################################################################
# Evaluate a floating point number expression.

function float_eval()
{
    local stat=0
    local result=0.0
    if [[ $# -gt 0 ]]; then
        result=$(echo "scale=$float_scale; $*" | bc -q 2>/dev/null)
        stat=$?
        if [[ $stat -eq 0  &&  -z "$result" ]]; then stat=1; fi
    fi
    echo $result
    return $stat
}

 SHEBANGLINE="#!/bin/bash"
MPIQUEUELINE="#PBS -q standard_4"
 MPITIMELINE="#PBS -l walltime=04:00:00"
 MPISIZELINE="#PBS -l nodes=1:ppn=4"
  MPIOUTLINE="#PBS -j oe"
SOURCEFILE="source /center/w/aschwand/FEniCS/share/fenics/fenics.conf"

for gamma in 1 2 5 10 20 50
do

  for alpha in 0 1 2
  do
#      for project in "jakobshavn" "79N"
      for project in "helheim"
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
      done
  done
done



echo
echo "($SPAWNSCRIPT)  use submitparam.sh to submit the scripts"
echo

