SRC = ..
#Directory of sources
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
SRC_LIST = $(shell ls *$(C_SUF))
FILES_C  = $(filter $(SRC_LIST), $(shell ls *$(C_SUF)))
#All .c/.cpp files
FILES_O  = $(FILES_C:$(C_SUF)=.o)
#All .o files
FILES_SLO = $(FILES_C:$(C_SUF)=.slo)
#All .slo files

define compile_c_file #Compile all .c/.cpp files
@echo "$(CC) $(CFLAGS) $(CFLAGD)-c {"
@for file in $(FILES_C); do \
( echo "  $$file"  && \
	$(CC) $(CFLAGS) $(CFLAGD) -c $$file $(INCLUDE) && \
	$(CC) -fPIC $(CFLAGS) $(CFLAGD) -c $$file -o $$(basename $$file $(C_SUF)).slo $(INCLUDE) 2> /dev/null \
) \
done;
@echo " } $(INCLUDE)"
endef

define rm_o_file #Delete all .o .slo files
@for file in $(FILES_O) $(FILES_SLO); do \
	( $(RM) $$file ) \
done;
endef

lib$(DIR_NAME).a lib$(DIR_NAME).so all: *.o
	@echo "******Producing static and dynamic libs******"
	@$(AR) $(ARFLAGS) lib$(DIR_NAME).a  $(FILES_O)
	@echo "$(AR) $(ARFLAGS) lib$(DIR_NAME).a"
	@$(CC) -shared -o lib$(DIR_NAME).so $(FILES_SLO)
	@echo "$(CC) -shared -o lib$(DIR_NAME).so"
.PHONYP:all

*.o *.slo: *$(C_SUF) $(addsuffix /*, $(addprefix $(SRC)/,$(INCLUDE_FOLDER)))
	@echo "******************Compiling******************"
	$(call compile_c_file)

clean:
	@echo "******************Cleaning*******************"
	$(call rm_o_file)
	@$(RM) lib$(DIR_NAME).a lib$(DIR_NAME).so
	@$(RM) *.gcov *.gcda *.gcno
.PHONYP:clean
