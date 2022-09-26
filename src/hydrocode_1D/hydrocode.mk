SRC = ..
#Directory of sources
INCLUDE = $(addprefix -I, $(addprefix $(SRC)/,$(INCLUDE_FOLDER)))
#Inclued folder
LIBS += $(addprefix -l, $(HEAD))
#Library files

RM = rm -vf
#Delete command
CP = cp -vu
#Copy command


all: modules libscopy $(SOURCE).out
.PHONYP:all

$(SOURCE).out exe: $(SOURCE).c lib/*.a $(addsuffix /*, $(addprefix $(SRC)/,$(INCLUDE_FOLDER)))
	@echo "**********Generate executable file***********"
	$(CC) $(CFLAGS) $(CFLAGD) -o $(SOURCE).out $(SOURCE).c $(INCLUDE) $(LIBS)
.PHONYP:exe

modules:
#Enter each subdirectory
#Call the Makefile in the subdirectory
	@for n in $(HEAD); do \
	( $(MAKE) CC=$(CC) CFLAGS='$(CFLAGS) $(CFLAGD)' SRC_LIST='$(SRC_LIST)' --directory=$(SRC)/$$n )  \
	done;
.PHONYP:modules

lib/*.a lib/*.so libscopy:
	@mkdir -pv lib
	@$(CP) $(addsuffix /*.a, $(addprefix $(SRC)/,$(HEAD))) lib/
	@$(CP) $(addsuffix /*.so,$(addprefix $(SRC)/,$(HEAD))) lib/
.PHONY: libscopy

get: 
#Generate GCOV code coverage report
	gcov $(SOURCE).c
.PHONY: get

html:
#Generate graphical GCOV code coverage report
	#Create test coverage data file
	lcov -c -d .. -o $(SOURCE).info
	#Information generated by LCOV is generated into HTML
	genhtml -o gcovdir $(SOURCE).info
.PHONYP:html

clean_all: clean
#Clean in the directory
	@$(RM) $(SOURCE).out $(SOURCE).exe
	@$(RM) -R lib gcovdir
	@$(RM) pg.png callgrind.png
.PHONYP:clean_all

clean:
#Clean in the subdirectory
	@for n in $(HEAD) $(HEAD_SO); do \
	( $(MAKE) --directory=$(SRC)/$$n clean ) \
	done;
	@$(RM) lib/*.a
	@$(RM) *.gcov *.gcda *.gcno
	@$(RM) $(SOURCE).info gmon.out pg callgrind.out
.PHONYP:clean
