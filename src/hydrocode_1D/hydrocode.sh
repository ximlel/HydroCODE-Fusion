#!/bin/bash

ulimit -c unlimited

# make clean
### Compile the program statically
make static

### Test the program
sh shell/hydrocode_test.sh

### Compile the program
make

### Run the program
sh shell/hydrocode_run.sh
