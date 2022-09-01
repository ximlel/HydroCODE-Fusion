#!/bin/bash

ulimit -c unlimited
shopt -s expand_aliases

### Compile the program
# make clean
make


### Run the program
CPath=$(pwd)
MRun=`ls -l`
#alias MRun='~/Softwares/MATLAB/R2018a/bin/matlab -nojvm -nodisplay -nosplash -nodesktop'
alias MRun='octave'
EXEcute=./hydrocode.out

## RP2D_Positive
:<<!
 cd ../../data_in/two-dim/RP2D_Positive/Config3
 echo "value_start" | MRun
 cd $CPath
 $EXEcute RP2D_Positive/Config3  RP2D_Positive/Config3  2 2_GRP EUL
 $EXEcute RP2D_Positive/Config7  RP2D_Positive/Config7  2 2_GRP EUL
 $EXEcute RP2D_Positive/Config12 RP2D_Positive/Config12 2 2_GRP EUL
!


### gprof
gprof -b -A -p -q hydrocode.out gmon.out > pg
gprof -b ./hydrocode.out gmon.out | gprof2dot | dot -Tpng -o pg.png
