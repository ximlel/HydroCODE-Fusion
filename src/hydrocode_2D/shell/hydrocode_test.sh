#!/bin/bash

export LD_LIBRARY_PATH=lib:$LD_LIBRARY_PATH
export OMP_STACKSIZE=8192
export OMP_NUM_THREADS=6

### Test the program
EXE=./hydrocode.out  #EXEcutable program

### Valgrind
# valgrind --tool=callgrind --callgrind-out-file=callgrind.out \
$EXE RP2D_Positive/Config3  RP2D_Positive/Config3 2_GRP EUL
# gprof2dot -f callgrind -s callgrind.out | dot  -Tpng -o callgrind.png

### gprof
# gprof -b -A -p -q $EXE gmon.out > pg
# gprof -b $EXE gmon.out | gprof2dot | dot -Tpng -o pg.png

### gcov 
# make html
