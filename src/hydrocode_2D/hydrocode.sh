#!/bin/bash

ulimit -c unlimited
shopt -s expand_aliases
export LD_LIBRARY_PATH=lib:$LD_LIBRARY_PATH

### Compile the program
# make clean
make


### Run the program
CPath=$(pwd)
MRun=`ls -l`
#alias MRun='~/Softwares/MATLAB/R2018a/bin/matlab -nojvm -nodisplay -nosplash -nodesktop'
alias MRun='octave'
EXE=./hydrocode.out  #EXEcutable program

 $EXE GRP_Book/6_1_Sod_10_lines   GRP_Book/6_1_Sod_10_lines   1     EUL 33=1
 $EXE GRP_Book/6_1_Sod_10_lines   GRP_Book/6_1_Sod_10_lines   2_GRP EUL 33=1
 $EXE GRP_Book/6_1_Sod_10_columns GRP_Book/6_1_Sod_10_columns 1     EUL 33=1
 $EXE GRP_Book/6_1_Sod_10_columns GRP_Book/6_1_Sod_10_columns 2_GRP EUL 33=1
 $EXE RP2D_Positive/Config3  RP2D_Positive/Config3  2_GRP EUL
 $EXE RP2D_Positive/Config7  RP2D_Positive/Config7  2_GRP EUL
 $EXE RP2D_Positive/Config12 RP2D_Positive/Config12 2_GRP EUL

## GRP_Book
:<<!
 $EXE GRP_Book/6_1_Sod_10_lines   GRP_Book/6_1_Sod_10_lines   1     EUL 33=1
 $EXE GRP_Book/6_1_Sod_10_lines   GRP_Book/6_1_Sod_10_lines   2_GRP EUL 33=1
 $EXE GRP_Book/6_1_Sod_10_columns GRP_Book/6_1_Sod_10_columns 1     EUL 33=1
 $EXE GRP_Book/6_1_Sod_10_columns GRP_Book/6_1_Sod_10_columns 2_GRP EUL 33=1
!
## RP2D_Positive
:<<!
 $EXE RP2D_Positive/Config3  RP2D_Positive/Config3  2_GRP EUL
 cd ../../data_in/two-dim/RP2D_Positive/Config7
#echo "value_start('INPUT')" | MRun
 echo "value_start(400)" | MRun
 cd $CPath
 $EXE RP2D_Positive/Config7  RP2D_Positive/Config7  2_GRP EUL
 $EXE RP2D_Positive/Config12 RP2D_Positive/Config12 2_GRP EUL
!


### gprof
# gprof -b -A -p -q $EXE gmon.out > pg
# gprof -b $EXE gmon.out | gprof2dot | dot -Tpng -o pg.png

### Valgrind
# valgrind --tool=callgrind --callgrind-out-file=callgrind.out \
# $EXE RP2D_Positive/Config3  RP2D_Positive/Config3 2_GRP EUL
# gprof2dot -f callgrind -s callgrind.out | dot  -Tpng -o callgrind.png

### gcov 
# make html
