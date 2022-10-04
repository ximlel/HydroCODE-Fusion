#!/bin/bash

shopt -s expand_aliases
export LD_LIBRARY_PATH=lib:$LD_LIBRARY_PATH
#export OMP_STACKSIZE=8192
#export OMP_NUM_THREADS=4

### Run the program
CPath=$(pwd)
#alias MRun='~/Softwares/MATLAB/R2018a/bin/matlab -nojvm -nodisplay -nosplash -nodesktop'
alias MRun='octave'
EXE=./hydrocode.out  #EXEcutable program
DI=../../data_in/two-dim

## GRP_Book
:<<!
 $EXE GRP_Book/6_1_Sod_10_lines   GRP_Book/6_1_Sod_10_lines   1     EUL 33=1
 $EXE GRP_Book/6_1_Sod_10_lines   GRP_Book/6_1_Sod_10_lines   2_GRP EUL 33=1
 $EXE GRP_Book/6_1_Sod_10_columns GRP_Book/6_1_Sod_10_columns 1     EUL 33=1
 $EXE GRP_Book/6_1_Sod_10_columns GRP_Book/6_1_Sod_10_columns 2_GRP EUL 33=1
!
## RP2D_Positive
:<<!
   cd $DI/RP2D_Positive/Config3
   echo "value_start" | MRun
   cd $CPath
 $EXE RP2D_Positive/Config3  RP2D_Positive/Config3  2_GRP EUL
   cd $DI/RP2D_Positive/Config7
#  echo "value_start('INPUT')" | MRun
   echo "value_start(400)" | MRun
   cd $CPath
 $EXE RP2D_Positive/Config7  RP2D_Positive/Config7  2_GRP EUL
   cd $DI/RP2D_Positive/Config12
   echo "value_start" | MRun
   cd $CPath
 $EXE RP2D_Positive/Config12 RP2D_Positive/Config12 2_GRP EUL
!
