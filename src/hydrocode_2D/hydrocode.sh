#!/bin/bash

### Compile the program
# make clean
make


### Run the program
EXEcute=./hydrocode.out
## 2D Riemnnn problem
#:<<!
 $EXEcute RP2D/Riemann_2D3_Quad RP2D/Riemann_2D3_Quad 2 2_GRP EUL
#!
:<<SHARE

SHARE
