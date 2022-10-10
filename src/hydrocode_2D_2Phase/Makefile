CC = gcc
#C compiler
CFLAGS = -std=c99 -Wall -pedantic -Og -g -fopenmp #-fopenacc
CFLAGR = -std=c99 -O2 -fopenmp
#C compiler options
CFLAGD = -DMULTIFLUID_BASICS -DMULTIPHASE_BASICS -DNODATPLOT
#Macro definition
INCLUDE_FOLDER = include
#Inclued folder
LDFLAGS = -lm -lhdf5
#Library files

#Head folder
HEAD = finite_volume flux_calc inter_process_BN riemann_solver_BN file_io tools
#Name of header files or subdirectories
SOURCE = hydrocode
#Name of the main source

SRC_LIST = sys_pro.c \
	config_handle.c file_1D_out.c terminal_io.c file_1D_in.c io_control.c \
	sys_pro.c \
	linear_grp_solver_LAG.c linear_grp_solver_Edir.c riemann_solver_exact_Ben.c riemann_solver_exact_Toro.c \
	bound_cond_slope_limiter.c fluid_var_check.c slope_limiter.c \
	grp_solver_EUL_source.c grp_solver_LAG_source.c godunov_solver_EUL_source.c godunov_solver_LAG_source.c
#List of source files

include ../MAKE/hydrocode.mk
