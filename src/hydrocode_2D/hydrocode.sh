#!/bin/bash

ulimit -c unlimited

### Compile the program
# make clean
make


### Run the program
EXEcute=./hydrocode.out
## 2D Riemnnn problem
#:<<!
 $EXEcute RP2D/RP2D_3_Quad RP2D/RP2D_3_Quad 2 2_GRP EUL
#!
:<<SHARE

SHARE
