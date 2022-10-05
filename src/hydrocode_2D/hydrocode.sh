#!/bin/bash

ulimit -c unlimited

### Compile the program statically
# make clean
make STATIC=1

### Test the program
sh shell/hydrocode_test.sh

### Compile the program dynamically
# make clean
# make

### Run the program
sh shell/hydrocode_run.sh

### Release the program dynamically
make clean
# make RELEASE=1
