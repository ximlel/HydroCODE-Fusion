#!/bin/bash

ulimit -c unlimited
export LD_LIBRARY_PATH=lib:$LD_LIBRARY_PATH

### Compile the program
# make clean
make


### Run the program
EXE=./hydrocode.out  #EXEcutable program

 $EXE GRP_Book/6_2_1   GRP_Book/6_2_1   2_GRP EUL 1=-1.0
 $EXE GRP_Book/6_2_3   GRP_Book/6_2_3   2_GRP EUL 1=-1.0
 $EXE GRP_Book/6_1_EUL GRP_direct/9_1_a 1     EUL 1=15
 $EXE GRP_Book/6_1_EUL GRP_direct/9_1_a 2_GRP EUL 1=15
 $EXE GRP_direct/9_1_e GRP_direct/9_1_e 2_GRP EUL

## GRP_Book
:<<!
 $EXE GRP_Book/6_1_LAG GRP_Book/6_1_LAG 1     LAG 1=-1.0
 $EXE GRP_Book/6_1_LAG GRP_Book/6_1_LAG 2_GRP LAG 1=-1.0
 $EXE GRP_Book/6_1_EUL GRP_Book/6_1_EUL 1     EUL 1=-1.0
 $EXE GRP_Book/6_1_EUL GRP_Book/6_1_EUL 2_GRP EUL 1=-1.0
 $EXE GRP_Book/6_2_1   GRP_Book/6_2_1   2_GRP EUL 1=-1.0
 $EXE GRP_Book/6_2_3   GRP_Book/6_2_3   2_GRP EUL 1=-1.0
!
:<<SHARE
 $EXE GRP_Book/6_1_EUL GRP_direct/9_1_a 1     EUL 1=15
 $EXE GRP_Book/6_1_EUL GRP_direct/9_1_a 2_GRP EUL 1=15
SHARE
## GRP_direct
:<<!
 $EXE GRP_direct/9_1_b GRP_direct/9_1_b 2_GRP EUL 41=1.0
 $EXE GRP_direct/9_1_d GRP_direct/9_1_d 2_GRP EUL
 $EXE GRP_direct/9_1_e GRP_direct/9_1_e 2_GRP EUL
!
:<<SHARE
 $EXE GRP_direct/9_1_b Artifi_Heat_Conduct/4_2 1 LAG 1=50 7=0.25
 $EXE GRP_direct/9_1_e Artifi_Heat_Conduct/4_4 1 LAG
SHARE
## Artifi_Heat_Conduct
:<<!
 $EXE Artifi_Heat_Conduct/4_3 Artifi_Heat_Conduct/4_3 1     LAG 7=0.25
 $EXE Artifi_Heat_Conduct/4_3 Artifi_Heat_Conduct/4_3 2_GRP LAG 7=0.25
 $EXE Artifi_Heat_Conduct/4_5 Artifi_Heat_Conduct/4_5 1     LAG 7=0.25
 $EXE Artifi_Heat_Conduct/4_5 Artifi_Heat_Conduct/4_5 2_GRP LAG 7=0.25
!


### gprof
# gprof -b -A -p -q $EXE gmon.out > pg
# gprof -b $EXE gmon.out | gprof2dot | dot -Tpng -o pg.png

### Valgrind
# valgrind --tool=callgrind --callgrind-out-file=callgrind.out \
$EXE GRP_Book/6_1_LAG GRP_Book/6_1_LAG 1 LAG 1=-1.0
# gprof2dot -f callgrind -s callgrind.out | dot  -Tpng -o callgrind.png

### gcov 
# make html
