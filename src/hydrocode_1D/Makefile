CC = gcc
#C compiler
CFLAGS = -std=c99 -Wall -pedantic -O4 -g -pg
#C compiler options
CFLAGD = #-D_Bool=bool
#Macro definition
SRC = ..
#Directory of sources
INCLUDE = -I $(SRC)/include
#Inclued folder
LIBS = -lm
#Library files

RM = rm -vf
#Delete command

HEAD = finite_volume inter_process Riemann_solver file_io tools
#Name of header files or subdirectories
SOURCE = hydrocode
#Name of the main source 

all: modules
	@echo "**********Generate executable file***********"
	$(CC) $(CFLAGS) $(CFLAGD) -o $(SOURCE).out $(SOURCE).c $(addsuffix /*.a,$(addprefix $(SRC)/,$(HEAD))) $(INCLUDE) $(LIBS)
.PHONYP:all

modules:
#Enter each subdirectory
#Call the Makefile in the subdirectory
	@for n in $(HEAD); do \
	( $(MAKE) CC=$(CC) CFLAGS='$(CFLAGS) $(CFLAGD)' --directory=$(SRC)/$$n )  \
	done;
.PHONYP:modules

clean:
#Clean in the subdirectory
	@$(RM) $(SOURCE).o
	@for n in $(HEAD); do \
	( $(MAKE) --directory=$(SRC)/$$n clean ) \
	done;
.PHONYP:clean
