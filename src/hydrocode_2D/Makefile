CC = gcc
#C compiler
CFLAGS = -std=c99 -O3 -pg -Wall
#C compiler options
CFLAGS += -D Riemann_solver_exact_single=Riemann_solver_exact_Ben
#Macro definition
SRC = ..
#Directory of sources
INCLUDE = -I $(SRC)/include
#Inclued folder
LIBS = -lm
#Library files

RM = rm -vf
#Delete command

HEAD = finite_volume flux_calc inter_process Riemann_solver file_io tools
#Name of header files or subdirectories
SOURCE = hydrocode
#Name of the main source 

all: modules
	@echo "**********Generate executable file***********"
	$(CC) $(CFLAGS) -o $(SOURCE).out $(SOURCE).c $(addsuffix /*.a,$(addprefix $(SRC)/,$(HEAD))) $(INCLUDE) $(LIBS)
.PHONYP:all

modules:
#Enter each subdirectory
#Call the Makefile in the subdirectory
	@for n in $(HEAD); do \
	( $(MAKE) CFLAGS='$(CFLAGS)' --directory=$(SRC)/$$n )  \
	done;
.PHONYP:modules

clean:
#Clean in the subdirectory
	@$(RM) $(SOURCE).o
	@for n in $(HEAD); do \
	( $(MAKE) --directory=$(SRC)/$$n clean ) \
	done;
.PHONYP:clean
