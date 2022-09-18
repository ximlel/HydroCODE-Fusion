FOLDER = $(shell ls -d */)
FOLDER_hydrocode = $(shell ls -d hydrocode_*/)

clean_all: clean 
#Clean in the subdirectory
	@for n in $(FOLDER_hydrocode) ; do \
	( $(MAKE) --directory=$$n clean_all ) \
	done;
.PHONYP:clean_all

clean:
#Clean in the subdirectory
	@for n in $(FOLDER) ; do \
	( $(MAKE) --directory=$$n clean ) \
	done;
.PHONYP:clean
