#!/bin/bash

export LD_LIBRARY_PATH=lib:$LD_LIBRARY_PATH
### Valgrind memcheck --debug
export DEBUGINFOD_URLS="https://debuginfod.archlinux.org"
export G_SLICE=always-malloc
export G_DEBUG=gc-friendly
### OpenMP
#export OMP_STACKSIZE=8192
#export OMP_NUM_THREADS=1

### Test the program
EXE=./hydrocode.out  #EXEcutable program
TEST="$EXE GRP_Book/6_1_Sod_10_lines   GRP_Book/6_1_Sod_10_lines   1     EUL 33=1"

### Perf
# perf record -e cpu-clock -g -F 999 $TEST
# perf script -i perf.data > perf.unfold
# stackcollapse-perf.pl perf.unfold &> perf.floded
# flamegraph.pl perf.floded > perf_flame.svg

### Valgrind
# valgrind --tool=callgrind --callgrind-out-file=callgrind.out $TEST
# gprof2dot -f callgrind -s callgrind.out | dot -Tpng -o callgrind.png
# valgrind --tool=cachegrind --cachegrind-out-file=cachegrind.out $TEST
# valgrind -v --tool=massif --time-unit=B --detailed-freq=1 --keep-debuginfo=yes -s --trace-children=yes --track-fds=yes --massif-out-file=massif.out $TEST
# valgrind -v --tool=memcheck --leak-check=full --show-reachable=yes --track-origins=yes --log-file=memchk.log $TEST

# $TEST

### gprof
# gprof -b -A -p -q $EXE gmon.out > pg
# gprof -b $EXE gmon.out | gprof2dot | dot -Tpng -o pg.png

### gcov 
# make get
# make html
