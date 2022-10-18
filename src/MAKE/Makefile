CC = g++
#C compiler
CFLAGS  = -std=c++20 -Wall -pedantic -Og -g -fopenmp #-fopenacc
CFLAGS += --coverage -pg
#C compiler options
CFLAGD = -D_Bool=bool
#Macro definition
SRC = ..
#Directory of sources
COPY_SRC = ../../src~
#Temporary directory of sources
RM = rm -vf
#Delete command

FOLDER = $(filter-out $(shell ls -d $(SRC)/include_*/) $(SRC)/MAKE/, $(shell ls -d $(SRC)/*/))
#Directory of source codes
FOLDER_hydro = $(shell ls -d $(SRC)/hydrocode_*/)
#Directory of hydrocodes

clean:
#Clean in the subdirectory
	@for n in $(FOLDER) ; do \
	( $(MAKE) --directory=$$n clean && \
	$(RM) $$n/*.o $$n/*.slo ) \
	done;
.PHONYP:clean

clean_all: clean 
#Clean all in the subdirectory
	@for n in $(FOLDER_hydro) ; do \
	( $(MAKE) --directory=$$n clean_all ) \
	done;
	@$(RM) -R $(COPY_SRC)
.PHONYP:clean_all

modules:
	@for n in $(filter-out $(FOLDER_hydro), $(FOLDER)) ; do \
	( $(MAKE) CC=$(CC) CFLAGS='$(CFLAGS) $(CFLAGD)' --directory=$$n )  \
	done;
.PHONYP:modules

all:
	@for n in $(FOLDER_hydro); do \
	( $(MAKE) --directory=$$n clean && \
	$(MAKE) --directory=$$n )  \
	done;
.PHONYP:all
