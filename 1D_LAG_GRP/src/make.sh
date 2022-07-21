#!/bin/bash

cd ./file_io/
gcc -c ./comman.c
gcc -c ./_1D_f_io.c -I ../
ar crv file_io.a comman.o _1D_f_io.o

cd ../Riemann_solver/
gcc -c ./Riemann_solver_exact.c -g


cd ../finite_difference_solver/
gcc -c ./GRP_solver_source.c -I ../ -g
gcc -c ./linear_GRP_solver_LAG.c -I ../ -g
ar crv finite_difference_solver.a GRP_solver_source.o linear_GRP_solver_LAG.o
ranlib finite_difference_solver.a

cd ../
gcc -c ./LAG_source.c -g
gcc -o LAG_source ./LAG_source.o ./file_io/file_io.a ./finite_difference_solver/finite_difference_solver.a ./Riemann_solver/Riemann_solver_exact.o -lm

exit 0
