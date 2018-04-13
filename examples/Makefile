#
#  This file is part of MUMPS 4.9, built on Wed Jul 29 10:35:58 UTC 2009
#
topdir = ..
libdir = $(topdir)/lib

default: d

.PHONY: default all s d c z clean
.SECONDEXPANSION:

all:	s d c z

s:	ssimpletest
d:	dsimpletest c_example
c:	csimpletest
z:	zsimpletest


include $(topdir)/Makefile.inc

LIBMUMPS_COMMON = $(libdir)/libmumps_common$(PLAT).a


LIBSMUMPS = $(libdir)/libsmumps$(PLAT).a $(LIBMUMPS_COMMON)

ssimpletest:  $(LIBSMUMPS) $$@.o 
	$(FL) -o $@ $(OPTL) ssimpletest.o $(LIBSMUMPS)  $(LIB) $(LIBBLAS) $(LIBOTHERS) $(LORDERINGS)


LIBDMUMPS = $(libdir)/libdmumps$(PLAT).a $(LIBMUMPS_COMMON)

dsimpletest: $(LIBDMUMPS) $$@.o 
	$(FL) -o $@ $(OPTL) dsimpletest.o $(LIBDMUMPS)  $(LIB) $(LIBBLAS) $(LIBOTHERS) $(LORDERINGS)


LIBCMUMPS = $(libdir)/libcmumps$(PLAT).a $(LIBMUMPS_COMMON)

csimpletest: $(LIBCMUMPS) $$@.o 
	$(FL) -o $@ $(OPTL) csimpletest.o $(LIBCMUMPS)  $(LIB) $(LIBBLAS) $(LIBOTHERS) $(LORDERINGS)


LIBZMUMPS = $(libdir)/libzmumps$(PLAT).a $(LIBMUMPS_COMMON)

zsimpletest: $(LIBZMUMPS) $$@.o 
	$(FL) -o $@ $(OPTL) zsimpletest.o $(LIBZMUMPS)  $(LIB) $(LIBBLAS) $(LIBOTHERS) $(LORDERINGS)


c_example:	$(LIBDMUMPS) $$@.o
	$(FL) -o $@ $(OPTL) $@.o $(LIBDMUMPS) $(LIB) $(LIBBLAS) $(LIBOTHERS) $(LORDERINGS)


.SUFFIXES: .c .F .o
.F.o:
	$(FC) $(OPTF) $(INC) -I. -I$(topdir)/include -c $*.F
.c.o:
	$(CC) $(OPTC) $(INC) -I. -I$(topdir)/include -c $*.c


$(libdir)/libsmumps$(PLAT).a:
	@echo 'Error: you should build the library' $@ 'first'
	exit -1

$(libdir)/libdmumps$(PLAT).a:
	@echo 'Error: you should build the library' $@ 'first'
	exit -1

$(libdir)/libcmumps$(PLAT).a:
	@echo 'Error: you should build the library' $@ 'first'
	exit -1

$(libdir)/libzmumps$(PLAT).a:
	@echo 'Error: you should build the library' $@ 'first'
	exit -1

$(LIBMUMPS_COMMON):
	@echo 'Error: you should build the library' $@ 'first'
	exit -1

clean:
	$(RM) *.o [sdcz]simpletest c_example
