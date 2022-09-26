CC = gcc
#C compiler
CFLAGS   = -std=c99 -Wall -pedantic -Og -g -pg
#CFLAGS += -fprofile-arcs -ftest-coverage
#C compiler options
CFLAGD = -DHDF5PLOT -DNODATPLOT
#Macro definition
INCLUDE_FOLDER = include
#Inclued folder
LIBS  = -lm -lgcov -lhdf5 -Llib/
#Library files

#Head folder
HEAD = finite_volume inter_process riemann_solver file_io tools
#Name of header files or subdirectories
SOURCE = hydrocode
#Name of the main source 

SRC_LIST = sys_pro.c \
	config_handle.c file_out_hdf5.c file_1D_out.c file_1D_in.c io_control.c terminal_io.c file_out_hdf5 \
	linear_grp_solver_LAG.c linear_grp_solver_Edir.c riemann_solver_exact_Ben.c riemann_solver_exact_Toro.c \
	bound_cond_slope_limiter.c slope_limiter.c fluid_var_check.c \
	grp_solver_EUL_source.c grp_solver_LAG_source.c godunov_solver_EUL_source.c godunov_solver_LAG_source.c
#List of source files

include hydrocode.mk
