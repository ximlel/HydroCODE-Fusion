#!/bin/bash

export LD_LIBRARY_PATH=lib:$LD_LIBRARY_PATH
#export OMP_STACKSIZE=8192
#export OMP_NUM_THREADS=4

### Run the program
CPath=$(pwd)
#MRun="$HOME/Softwares/MATLAB/R2018a/bin/matlab -nojvm -nodisplay -nosplash -nodesktop"
MRun="octave"
EXE=./hydrocode.out  #EXEcutable program
RSTC=Radial_Symmetry/Two_Component
DI=../../data_in/one-dim/$RSTC

## A3-shell
#:<<!
   cd $DI/A3_shell
   echo "value_start" | $MRun
   cd $CPath
 $EXE $RSTC/A3_shell $RSTC/A3_shell 2_GRP 2
#!
