CC = g++
#C compiler
CFLAGS += -std=c++20 -Wall
#C compiler options
CFLAGS += -D_Bool=bool
CFLAGD  = 
#Macro definition
SRC = ..
#Directory of sources
INCLUDE = -I$(SRC)/include -I$(SRC)/include_cpp
#Inclued folder

AR = ar
#Archive : creat static libs
ARFLAGS = crsv
#AR options
RM = rm -vf
#Delete command

DIR_NAME = $(shell basename `pwd`)
#Current folder
FILES_C = $(shell ls *.cpp)
#All .c files

define compile_c_file #Compile all .cpp files
@for file in $(FILES_C); do \
( echo "$(CC) $(CFLAGS) $(CFLAGD) -c $$file $(INCLUDE)" && \
	$(CC) $(CFLAGS) $(CFLAGD) -c $$file $(INCLUDE)  && \
	$(CC) -fPIC $(CFLAGS) $(CFLAGD) -c $$file -o $$file.slo $(INCLUDE) ) \
done;
endef

define rm_o_file #Delete all .o .slo files
@for file in $(FILES_O); do \
( $(RM) $$file ) \
done;
@for file in $(FILES_SLO); do \
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
FILES_SLO = $(shell find ./ -name "*.slo" -print)
#All .slo files

archive:
	@echo "************Producing static libs************"
	@$(AR) $(ARFLAGS) lib$(DIR_NAME).a  $(FILES_O)
	@$(CC) -shared -o lib$(DIR_NAME).so $(FILES_SLO)
./PHONYP:archive

clean:
	@echo "******************Cleaning*******************"
	$(call rm_o_file)
	@$(RM) lib$(DIR_NAME).a lib$(DIR_NAME).so
	@$(RM) *.gcov *.gcda *.gcno
.PHONYP:clean
