CC = g++
#C compiler
CFLAGS   = -std=c++20 -fopenmp -Wall -pedantic -Og -g -pg
#C compiler options
CFLAGD = -D_Bool=bool
#Macro definition
SRC = .
#Directory of sources
RM = rm -vf
#Delete command

FOLDER = $(filter-out $(shell ls -d include_*/), $(shell ls -d */))
#Directory of source codes
FOLDER_hydro = $(shell ls -d hydrocode_*/)
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
.PHONYP:clean_all

modules:
	@for n in $(filter-out $(FOLDER_hydro), $(FOLDER)) ; do \
	( $(MAKE) CC=$(CC) CFLAGS='$(CFLAGS) $(CFLAGD)' --directory=$(SRC)/$$n )  \
	done;
.PHONYP:modules

all:
	@for n in $(FOLDER_hydro); do \
	( $(MAKE) --directory=$(SRC)/$$n clean && \
	$(MAKE) --directory=$(SRC)/$$n )  \
	done;
.PHONYP:all
