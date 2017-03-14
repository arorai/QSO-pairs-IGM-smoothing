SHELL = /bin/sh

prefix=$(QSO_DIR)
src=$(prefix)/src/trans/
include=$(prefix)/include/
lib=$(prefix)/lib


OBJECTS = $(src)/maketau_idl.o $(src)/maketaulight_idl.o $(src)/transtau_idl.o $(src)/lostau_idl.o \
	$(src)/kde_pdf_idl.o $(src)/kde_pdf_weight_idl.o \
	$(src)/qlf_calculator_idl.o

CFLAGS  = -I$(IDL_DIR)/external -I$(include) $(SDSS_CFLAGS) \
	-DCHECK_LEAKS

.c.o :
	$(CC) -c $(CCCHK) $(CFLAGS) $(X_CFLAGS) $*.c

all : $(LIB)/liblya.$(SO_EXT)

$(LIB)/liblya.$(SO_EXT): $(OBJECTS)
	$(LD) $(X_LD_FLAGS) \
	-lm -o $(lib)/liblya.$(SO_EXT) $(OBJECTS)


install : 
	-@ cp *.c $(src)/
	-@ cp *.h $(src)/
	-@ cp Makefile $(src)/

clean :
	- /bin/rm -f *~ core *.o so_locations
