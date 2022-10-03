#!/bin/bash

export LD_LIBRARY_PATH=lib:$LD_LIBRARY_PATH

### Test the program
EXE=./hydrocode.out  #EXEcutable program

### Valgrind
# valgrind --tool=callgrind --callgrind-out-file=callgrind.out \
$EXE GRP_Book/6_1_LAG GRP_Book/6_1_LAG 1 LAG 1=-1.0
# gprof2dot -f callgrind -s callgrind.out | dot  -Tpng -o callgrind.png

### gprof
# gprof -b -A -p -q $EXE gmon.out > pg
# gprof -b $EXE gmon.out | gprof2dot | dot -Tpng -o pg.png

### gcov 
# make html
