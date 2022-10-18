CC = gcc
#C compiler
CFLAGS += -std=c99 -Wall
#C compiler options
CFLAGS += 
CFLAGD  = 
#Macro definition
INCLUDE_FOLDER = include_cii
#Inclued folder
C_SUF = .c
#Source code file suffix

define copy_cii_h_file #Copy cii .h files
	@$(MKDIR) $(COPY_SRC)/include_cii
@for file in $(FILES_C); do \
	( if [ -f ../include_cii/$$(basename $$file $(C_SUF)).h ]; then $(CP) ../include_cii/$$(basename $$file $(C_SUF)).h $(COPY_SRC)/include_cii; fi ) \
done;
endef

include ../MAKE/module.mk
