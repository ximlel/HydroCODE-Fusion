#!/bin/bash

export LD_LIBRARY_PATH=lib:$LD_LIBRARY_PATH
export OMP_STACKSIZE=8192

### Test the program
EXE=./hydrocode.out  #EXEcutable program
TC=Two_Component

### Valgrind
# valgrind --tool=callgrind --callgrind-out-file=callgrind.out \
$EXE 
# gprof2dot -f callgrind -s callgrind.out | dot  -Tpng -o callgrind.png

### gprof
# gprof -b -A -p -q $EXE gmon.out > pg
# gprof -b $EXE gmon.out | gprof2dot | dot -Tpng -o pg.png

### gcov 
# make html