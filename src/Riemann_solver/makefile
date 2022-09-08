CC = gcc
#C compiler
CFLAGS += -std=c99 -Wall
#C compiler options
CFLAGS += -DMULTIFLUID_BASICS
CFLAGD  = -DRiemann_solver_exact_single=Riemann_solver_exact_Ben #-DEXACT_TANGENT_DERIVATIVE
#Macro definition
SRC = ..
#Directory of sources
INCLUDE = -I $(SRC)/include
#Inclued folder

AR = ar
#Archive : creat static libs
ARFLAGS = crsv
#AR options
RM = rm -vf
#Delete command

DIR_NAME = $(shell basename `pwd`)
#Current folder
FILES_C = $(shell ls *.c)
#All .c files

define compile_c_file #Compile all .c files
@for file in $(FILES_C); do \
( echo "$(CC) $(CFLAGS) $(CFLAGD) -c $$file $(INCLUDE)" && $(CC) $(CFLAGS) $(CFLAGD) -c $$file $(INCLUDE) ) \
done;
endef

define rm_o_file #Delete all .o files
@for file in $(FILES_O); do \
( $(RM) $$file ) \
done;
endef

all: compile archive
.PHONYP:all

compile:
	@echo "******************Compiling******************"
	$(call compile_c_file)
.PHONYP:compile

FILES_O = $(shell find ./ -name "*.o" -print)
#All .o files

archive:
	@echo "************Producing static libs************"
	@$(AR) $(ARFLAGS) $(DIR_NAME).a $(FILES_O)
./PHONYP:archive

clean:
	@echo "******************Cleaning*******************"
	$(call rm_o_file)
	@$(RM) $(DIR_NAME).a
	@$(RM) *.gcov *.gcda *.gcno
.PHONYP:clean
