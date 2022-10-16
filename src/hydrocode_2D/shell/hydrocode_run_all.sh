#!/bin/bash

export LD_LIBRARY_PATH=lib:$LD_LIBRARY_PATH
#export OMP_STACKSIZE=8192
#export OMP_NUM_THREADS=4

### Run the program
CPath=$(pwd)
#MRun="$HOME/Softwares/MATLAB/R2018a/bin/matlab -nojvm -nodisplay -nosplash -nodesktop"
MRun="octave"
EXE=./hydrocode.out  #EXEcutable program
DI=../../data_in/two-dim

## GRP_Book
#:<<!
 $EXE GRP_Book/6_1_Sod_10_lines   GRP_Book/6_1_Sod_10_lines   1     EUL 33=1
 $EXE GRP_Book/6_1_Sod_10_lines   GRP_Book/6_1_Sod_10_lines   2_GRP EUL 33=1
 $EXE GRP_Book/6_1_Sod_10_columns GRP_Book/6_1_Sod_10_columns 1     EUL 33=1
 $EXE GRP_Book/6_1_Sod_10_columns GRP_Book/6_1_Sod_10_columns 2_GRP EUL 33=1
#!
## RP2D_Positive
#:<<!
   echo "cd $DI/RP2D_Positive/Config3; value_start" | $MRun
 $EXE RP2D_Positive/Config3  RP2D_Positive/Config3  2_GRP EUL
#  echo "cd $DI/RP2D_Positive/Config7; value_start('INPUT')" | $MRun
   echo "cd $DI/RP2D_Positive/Config7; value_start(400)" | $MRun
 $EXE RP2D_Positive/Config7  RP2D_Positive/Config7  2_GRP EUL
   echo "cd $DI/RP2D_Positive/Config12; value_start" | $MRun
 $EXE RP2D_Positive/Config12 RP2D_Positive/Config12 2_GRP EUL
#!
