CC = gcc
#C compiler
CFLAGS = -std=c99 -Wall -pedantic -Og -g -fopenmp
CFLAGR = -std=c99 -O2 -fopenmp
#C compiler options
CFLAGD = -DMULTIFLUID_BASICS -DNOVTKPLOT
#Macro definition
INCLUDE_FOLDER = include
#Inclued folder
LDFLAGS = -lm -lhdf5
#Library files

#Head folder
HEAD = finite_volume flux_calc inter_process_unstruct inter_process riemann_solver meshing file_io tools src_cii
#Name of header files or subdirectories
SOURCE = hydrocode
#Name of the main source 

SRC_LIST = except.c mem.c \
	sys_pro.c mat_algo.c \
	config_handle.c file_2D_unstruct_out.c file_2D_in.c io_control.c terminal_io.c \
	ghost_cell.c  mesh_init_free.c  msh_load.c  quad_mesh.c \
	hll_2D_solver.c roe_2D_solver.c roe_solver.c roe_HLL_solver.c linear_grp_solver_Edir_G2D.c linear_grp_solver_Edir_Q1D.c riemann_solver_exact_Ben.c riemann_solver_exact_Toro.c \
	fluid_var_check.c \
	assist_func.c cons_qty_calc.c copy_func.c cell_init_free.c cons_qty_update_P_ave.c slope_limiter_unstruct.c \
	flux_solver.c \
	finite_volume_scheme_unstruct.c
#List of source files

include ../MAKE/hydrocode.mk
