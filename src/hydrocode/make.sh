#!/bin/bash

cd ../..
make
cd ./src/hydrocode

./hydrocode.out GRP_Book_6_1_LAG 1 LAG
