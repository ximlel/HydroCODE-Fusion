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

## 2D Riemnnn problem
:<<!
 cd ../../data_in/two-dim/RP2D/RP2D_3_Quad
 echo "value_start" | MRun
 cd $CPath
 $EXEcute RP2D/RP2D_3_Quad RP2D/RP2D_3_Quad 2 2_GRP EUL
!
 $EXEcute GRP_Book/6_1_Sod_10_lines GRP_Book/6_1_Sod_10_lines 2 2_GRP EUL 33=1
 $EXEcute GRP_Book/6_1_Sod_10_columns GRP_Book/6_1_Sod_10_columns 2 2_GRP EUL 33=1
#!
:<<SHARE

SHARE

### gprof
gprof -b -A -p -q hydrocode.out gmon.out > pg
gprof -b ./hydrocode.out gmon.out | gprof2dot | dot -Tpng -o pg.png
