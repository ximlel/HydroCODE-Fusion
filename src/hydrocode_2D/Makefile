CC = gcc
#C compiler
CFLAGS   = -std=c99 -fopenmp -Wall -pedantic -Og -g -pg
#CFLAGS += -fprofile-arcs -ftest-coverage
#C compiler options
CFLAGD = -DHDF5PLOT -DNODATPLOT
#Macro definition
INCLUDE_FOLDER = include
#Inclued folder
LIBS  = -lm -lgcov -lhdf5 -Llib/
#Library files

#Head folder
HEAD = finite_volume flux_calc inter_process riemann_solver file_io tools
#Name of header files or subdirectories
SOURCE = hydrocode
#Name of the main source 

SRC_LIST = sys_pro.c \
	config_handle.c file_out_hdf5.c file_2D_out.c file_2D_in.c io_control.c terminal_io.c \
	hll_2D_solver.c roe_2D_solver.c roe_solver.c roe_HLL_solver.c linear_grp_solver_Edir_G2D.c linear_grp_solver_Edir_Q1D.c riemann_solver_exact_Ben.c riemann_solver_exact_Toro.c \
	bound_cond_slope_limiter_x.c bound_cond_slope_limiter_y.c fluid_var_check.c slope_limiter.c slope_limiter_2D_x.c \
	flux_generator_x.c flux_generator_y.c flux_solver.c \
	grp_solver_2D_EUL_source.c grp_solver_2D_split_EUL_source.c \
#List of source files

include ../hydrocode_1D/hydrocode.mk
