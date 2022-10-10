#!/bin/bash

export LD_LIBRARY_PATH=lib:$LD_LIBRARY_PATH

### Test the program
EXE=./hydrocode.out  #EXEcutable program
RSTC=Radial_Symmetry/Two_Component
TEST="$EXE $RSTC/A3_shell $RSTC/A3_shell 2_GRP 2"

### Perf
# perf record -e cpu-clock -g -F 999 $TEST
# perf script -i perf.data > perf.unfold
# stackcollapse-perf.pl perf.unfold &> perf.floded
# flamegraph.pl perf.floded > perf_flame.svg

### Valgrind
# valgrind --tool=callgrind --callgrind-out-file=callgrind.out $TEST
# gprof2dot -f callgrind -s callgrind.out | dot -Tpng -o callgrind.png

 $TEST

### gprof
# gprof -b -A -p -q $EXE gmon.out > pg
# gprof -b $EXE gmon.out | gprof2dot | dot -Tpng -o pg.png

### gcov 
# make get
# make html
