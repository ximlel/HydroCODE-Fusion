CC = gcc
#C compiler
CFLAGS += -std=c99 -Wall
#C compiler options
CFLAGS += 
CFLAGD  = 
#Macro definition
SRC = ..
#Directory of sources
INCLUDE_FOLDER = include
INCLUDE = $(addprefix -I, $(addprefix $(SRC)/,$(INCLUDE_FOLDER)))
#Inclued folder

AR = ar
#Archive : creat static libs
ARFLAGS = crs
#AR options
RM = rm -vf
#Delete command

DIR_NAME = $(shell basename `pwd`)
#Current folder
SRC_LIST = $(shell ls *.c)
FILES_C  = $(filter $(SRC_LIST), $(shell ls *.c))
#All .c files
FILES_O  = $(addsuffix .o,$(basename $(FILES_C)))
#All .o files
FILES_SLO = $(addsuffix .slo,$(basename $(FILES_C)))
#All .slo files

define compile_c_file #Compile all .c files
@echo "$(CC) $(CFLAGS) $(CFLAGD)-c {"
@for file in $(FILES_C); do \
( echo "  $$file"  && \
	$(CC) $(CFLAGS) $(CFLAGD) -c $$file $(INCLUDE)  && \
	$(CC) -fPIC $(CFLAGS) $(CFLAGD) -c $$file -o $$(basename $$file .c).slo $(INCLUDE) 2> /dev/null ) \
done;
@echo " } $(INCLUDE)"
endef

define rm_o_file #Delete all .o .slo files
@for file in $(FILES_O) $(FILES_SLO); do \
	( $(RM) $$file ) \
done;
endef

all: lib$(DIR_NAME).a
.PHONYP:all

*.o *.slo: *.c $(addsuffix /*, $(addprefix $(SRC)/,$(INCLUDE_FOLDER)))
	@echo "******************Compiling******************"
	$(call compile_c_file)

lib$(DIR_NAME).a lib$(DIR_NAME).so: *.o
	@echo "******Producing static and dynamic libs******"
	@$(AR) $(ARFLAGS) lib$(DIR_NAME).a  $(FILES_O)
	@echo "$(AR) $(ARFLAGS) lib$(DIR_NAME).a"
	@$(CC) -shared -o lib$(DIR_NAME).so $(FILES_SLO)
	@echo "$(CC) -shared -o lib$(DIR_NAME).so"

clean:
	@echo "******************Cleaning*******************"
	$(call rm_o_file)
	@$(RM) lib$(DIR_NAME).a lib$(DIR_NAME).so
	@$(RM) *.gcov *.gcda *.gcno
.PHONYP:clean
