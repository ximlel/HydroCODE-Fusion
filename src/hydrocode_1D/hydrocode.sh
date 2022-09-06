#!/bin/bash

ulimit -c unlimited

### Compile the program
# make clean
make


### Run the program
EXEcute=./hydrocode.out

## GRP_Book
:<<!
 $EXEcute GRP_Book/6_1_LAG GRP_Book/6_1_LAG 1 1     LAG 1=-1.0
 $EXEcute GRP_Book/6_1_LAG GRP_Book/6_1_LAG 1 2_GRP LAG 1=-1.0
 $EXEcute GRP_Book/6_1_EUL GRP_Book/6_1_EUL 1 1     EUL 1=-1.0
 $EXEcute GRP_Book/6_1_EUL GRP_Book/6_1_EUL 1 2_GRP EUL 1=-1.0
 $EXEcute GRP_Book/6_2_1   GRP_Book/6_2_1   1 2_GRP EUL 1=-1.0
 $EXEcute GRP_Book/6_2_3   GRP_Book/6_2_3   1 2_GRP EUL 1=-1.0
!
:<<SHARE
 $EXEcute GRP_Book/6_1_EUL GRP_direct/9_1_a 1 1     EUL 1=15
 $EXEcute GRP_Book/6_1_EUL GRP_direct/9_1_a 1 2_GRP EUL 1=15
SHARE
## GRP_direct
:<<!
 $EXEcute GRP_direct/9_1_b GRP_direct/9_1_b 1 2_GRP EUL 41=1.0
 $EXEcute GRP_direct/9_1_d GRP_direct/9_1_d 1 2_GRP EUL
 $EXEcute GRP_direct/9_1_e GRP_direct/9_1_e 1 2_GRP EUL
!
:<<SHARE
 $EXEcute GRP_direct/9_1_b Artifi_Heat_Conduct/4_2 1 1 LAG 1=50 7=0.25
 $EXEcute GRP_direct/9_1_e Artifi_Heat_Conduct/4_4 1 1 LAG
SHARE
## Artifi_Heat_Conduct
:<<!
 $EXEcute Artifi_Heat_Conduct/4_3 Artifi_Heat_Conduct/4_3 1 1     LAG 7=0.25
 $EXEcute Artifi_Heat_Conduct/4_3 Artifi_Heat_Conduct/4_3 1 2_GRP LAG 7=0.25
 $EXEcute Artifi_Heat_Conduct/4_5 Artifi_Heat_Conduct/4_5 1 1     LAG 7=0.25
 $EXEcute Artifi_Heat_Conduct/4_5 Artifi_Heat_Conduct/4_5 1 2_GRP LAG 7=0.25
!

### gprof
gprof -b -A -p -q hydrocode.out gmon.out > pg
gprof -b ./hydrocode.out gmon.out | gprof2dot | dot -Tpng -o pg.png

### Valgrind
valgrind --tool=callgrind $EXEcute GRP_Book/6_1_LAG GRP_Book/6_1_LAG 1 1     LAG 1=-1.0
gprof2dot -f callgrind callgrind.out.* | dot  -Tpng -o callgrind.png
rm ./callgrind.out.*
