#!/bin/bash

cd ../file_io/
gcc -c ./common.c
gcc -c ./_1D_f_io.c -I ../
ar crvs file_io.a common.o _1D_f_io.o

cd ../Riemann_solver/
gcc -c ./Riemann_solver_exact_Ben.c -g
gcc -c ./Riemann_solver_exact_Toro.c -g
gcc -c ./linear_GRP_solver_LAG.c -I ../ -g
ar crvs Riemann_solver.a Riemann_solver_exact_Ben.o Riemann_solver_exact_Toro.o linear_GRP_solver_LAG.o

cd ../finite_difference_solver/
gcc -c ./Godunov_solver_source.c -I ../ -g
gcc -c ./GRP_solver_source.c -I ../ -g
ar crvs finite_difference_solver.a Godunov_solver_source.o GRP_solver_source.o

cd ../LAG_source
gcc -c ./LAG_source.c -g
gcc -o LAG_source.out ./LAG_source.o ../file_io/file_io.a ../finite_difference_solver/finite_difference_solver.a ../Riemann_solver/Riemann_solver.a -lm

exit 0
