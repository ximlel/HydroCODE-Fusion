RM = rm -vf
#Delete command

FOLDER = $(shell ls -d */)
#Directory of source codes
FOLDER_hydro = $(shell ls -d hydrocode_*/)
#Directory of hydrocodes

clean_all: clean 
#Clean all in the subdirectory
	@for n in $(FOLDER_hydro) ; do \
	( $(MAKE) --directory=$$n clean_all && \
	$(RM) $$n/*.o $$n/*.slo ) \
	done;
.PHONYP:clean_all

clean:
#Clean in the subdirectory
	@for n in $(FOLDER) ; do \
	( $(MAKE) --directory=$$n clean ) \
	done;
.PHONYP:clean
