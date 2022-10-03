#!/bin/bash

ulimit -c unlimited

### Compile the program
# make clean
make

### Test the program
sh shell/hydrocode_test.sh

### Run the program
sh shell/hydrocode_run.sh
