CC = gcc
#C compiler
CFLAGS = -std=c99 -Wall -pedantic -Og -g -fopenmp #-fopenacc
CFLAGR = -std=c99 -O2 -fopenmp
#C compiler options
#CC = /opt/nvidia/hpc_sdk/Linux_x86_64/2022/compilers/bin/nvcc
#CFLAGR = -std=c99 -O2 -acc -mp -ta=multicore -Minfo=accel
#NVIDIA HPC C compiler options
#CC = /opt/intel/oneapi/compiler/latest/linux/bin/icx
#CC = /opt/intel/oneapi/compiler/latest/linux/bin/intel64/icc
#CFLAGR = -std=c99 -O2 -qopenmp -shared-intel
#Intel C compiler options
CFLAGD = -DHDF5PLOT #-DNODATPLOT -DNOTECPLOT
#Macro definition
INCLUDE_FOLDER = include
#Inclued folder
LDFLAGS = -lm -lhdf5
#Library files

#Head folder
HEAD = finite_volume flux_calc inter_process riemann_solver file_io tools
#Name of header files or subdirectories
SOURCE = hydrocode
#Name of the main source

SRC_LIST = sys_pro.c \
	config_handle.c file_out_hdf5.c file_2D_out.c file_2D_in.c io_control.c terminal_io.c \
	hll_2D_solver.c roe_2D_solver.c roe_solver.c roe_hll_solver.c linear_grp_solver_Edir_G2D.c linear_grp_solver_Edir_Q1D.c riemann_solver_exact_Ben.c riemann_solver_exact_Toro.c \
	bound_cond_slope_limiter_x.c bound_cond_slope_limiter_y.c fluid_var_check.c slope_limiter.c slope_limiter_2D_x.c \
	flux_generator_x.c flux_generator_y.c flux_solver.c \
	grp_solver_2D_EUL_source.c grp_solver_2D_split_EUL_source.c
#List of source files

include ../MAKE/hydrocode.mk
