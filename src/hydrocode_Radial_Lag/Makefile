CC = g++
#C compiler
CFLAGS   = -std=c++20 -Wall -pedantic -Og -g -pg
#CFLAGS += -fprofile-arcs -ftest-coverage
#C compiler options
CFLAGD = -DRADIAL_BASICS -DMULTIFLUID_BASICS -D_Bool=bool
#Macro definition
INCLUDE_FOLDER = include
#Inclued folder
LIBS  = -lm -lgcov -lhdf5 -Llib/
#Library files

#Head folder
HEAD = finite_volume inter_process inter_process_cpp riemann_solver meshing file_io tools src_cii
#Name of header files or subdirectories
SOURCE = hydrocode
#Name of the main source 

SRC_LIST = except.c mem.c \
	sys_pro.c \
	config_handle.c file_out_hdf5.c file_radial_out.c file_1D_out.c file_1D_in.c io_control.c terminal_io.c \
	radial_mesh.c \
	linear_grp_solver_radial_LAG.c riemann_solver_starPU.c \
	VIPLimiter.cpp \
	fluid_var_check.c slope_limiter_radial.c slope_VIP_limiter_radial.c \
	grp_solver_radial_LAG_source.c
#List of source files

include ../hydrocode_1D/hydrocode.mk
