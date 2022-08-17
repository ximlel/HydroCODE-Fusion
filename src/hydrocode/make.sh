#!/bin/bash

cd ../..
make
cd ./src/hydrocode

./hydrocode.out Artifi_Heat_Conduct_4_22 1 LAG
