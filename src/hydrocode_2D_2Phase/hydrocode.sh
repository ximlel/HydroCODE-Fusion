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
EXE=./hydrocode.out  #EXEcutable program

## RP2D_Positive
:<<!
 cd ../../data_in/two-dim/RP2D_Positive/Config3
 echo "value_start" | MRun
 cd $CPath
 $EXE 
!


### gprof
# gprof -b -A -p -q $EXE gmon.out > pg
# gprof -b $EXE gmon.out | gprof2dot | dot -Tpng -o pg.png

### Valgrind
# valgrind --tool=callgrind --callgrind-out-file=callgrind.out \
$EXE 
# gprof2dot -f callgrind -s callgrind.out | dot  -Tpng -o callgrind.png

### gcov 
# make html
