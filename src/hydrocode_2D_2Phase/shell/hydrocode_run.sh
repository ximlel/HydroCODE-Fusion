#!/bin/bash

shopt -s expand_aliases
export LD_LIBRARY_PATH=lib:$LD_LIBRARY_PATH
export OMP_STACKSIZE=8192
#export OMP_NUM_THREADS=4

### Run the program
CPath=$(pwd)
#alias MRun='~/Softwares/MATLAB/R2018a/bin/matlab -nojvm -nodisplay -nosplash -nodesktop'
alias MRun='octave'
EXE=./hydrocode.out  #EXEcutable program
DI=../../data_in/two-dim

## RP2D_Positive
:<<!
   cd $DI/RP2D_Positive/Config3
   echo "value_start" | MRun
   cd $CPath
 $EXE 
!