SRC = ..
#Directory of sources
COPY_SRC = ../../src~
#Temporary directory of sources
INCLUDE = $(addprefix -I, $(addprefix $(SRC)/,$(INCLUDE_FOLDER)))
#Inclued folder

AR = ar
#Archive : creat static libs
ARFLAGS = crs
#AR options
RM = rm -vf
#Delete command
CP = cp -vur
#Copy command
MKDIR = mkdir -pv
#Make directory command

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

define copy_c_file #Copy all .c/.cpp files
@for file in $(FILES_C); do \
	( $(CP) $$file $(COPY_SRC)/$(DIR_NAME) ) \
done;
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

copy:
	@echo "******************Copying********************"
	@$(MKDIR) $(COPY_SRC)/include $(COPY_SRC)/include_cpp
	@$(MKDIR) $(COPY_SRC)/$(DIR_NAME)
	$(call copy_c_file)
	$(call copy_cii_h_file)
	@if [ -f $(SRC)/include/$(DIR_NAME).h ];       then $(CP) $(SRC)/include/$(DIR_NAME).h       $(COPY_SRC)/include;     fi
	@if [ -f $(SRC)/include_cpp/$(DIR_NAME).hpp ]; then $(CP) $(SRC)/include_cpp/$(DIR_NAME).hpp $(COPY_SRC)/include_cpp; fi

clean:
	@echo "******************Cleaning*******************"
	$(call rm_o_file)
	@$(RM) lib$(DIR_NAME).a lib$(DIR_NAME).so
	@$(RM) *.gcov *.gcda *.gcno
.PHONYP:clean
