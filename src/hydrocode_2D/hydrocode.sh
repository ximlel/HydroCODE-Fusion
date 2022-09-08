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

## GRP_Book
:<<!
 $EXEcute GRP_Book/6_1_Sod_10_lines   GRP_Book/6_1_Sod_10_lines   2 1     EUL 33=1
 $EXEcute GRP_Book/6_1_Sod_10_lines   GRP_Book/6_1_Sod_10_lines   2 2_GRP EUL 33=1
 $EXEcute GRP_Book/6_1_Sod_10_columns GRP_Book/6_1_Sod_10_columns 2 1     EUL 33=1
 $EXEcute GRP_Book/6_1_Sod_10_columns GRP_Book/6_1_Sod_10_columns 2 2_GRP EUL 33=1
!
## RP2D_Positive
:<<!
 $EXEcute RP2D_Positive/Config3  RP2D_Positive/Config3  2 2_GRP EUL
 cd ../../data_in/two-dim/RP2D_Positive/Config7
#echo "value_start('INPUT')" | MRun
 echo "value_start(400)" | MRun
 cd $CPath
 $EXEcute RP2D_Positive/Config7  RP2D_Positive/Config7  2 2_GRP EUL
 $EXEcute RP2D_Positive/Config12 RP2D_Positive/Config12 2 2_GRP EUL
!


### gprof
# gprof -b -A -p -q hydrocode.out gmon.out > pg
# gprof -b ./hydrocode.out gmon.out | gprof2dot | dot -Tpng -o pg.png

### Valgrind
# valgrind --tool=callgrind --callgrind-out-file=callgrind.out \
$EXEcute RP2D_Positive/Config3  RP2D_Positive/Config3  2 2_GRP EUL
# gprof2dot -f callgrind -s callgrind.out | dot  -Tpng -o callgrind.png

### gcov 
# make html
