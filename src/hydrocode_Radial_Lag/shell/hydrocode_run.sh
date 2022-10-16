#!/bin/bash

export LD_LIBRARY_PATH=lib:$LD_LIBRARY_PATH

### Run the program
CPath=$(pwd)
#MRun="$HOME/Softwares/MATLAB/R2018a/bin/matlab -nojvm -nodisplay -nosplash -nodesktop"
MRun="octave"
EXE=./hydrocode.out  #EXEcutable program
RSTC=Radial_Symmetry/Two_Component
DITC=../../data_in/one-dim/$RSTC

## A3-shell
:<<!
   echo "cd $DITC/A3_shell; value_start" | $MRun
 $EXE $RSTC/A3_shell $RSTC/A3_shell 2_GRP 2 42=-2
!
