#
#  This file is part of MUMPS 4.8.0, built on Fri Jul 25 14:46:02 2008
#
all: libmpiseq

.PHONY: all libmpiseq clean

include ../Makefile.inc

libmpiseq: libmpiseq$(PLAT).a

libmpiseq$(PLAT).a: mpi.o mpic.o elapse.o
	$(AR) libmpiseq$(PLAT).a mpi.o mpic.o elapse.o
.f.o:
	$(FC) $(OPTF) -c $*.f
.c.o:
	$(CC) $(OPTC) $(CDEFS) -I. -c $*.c

clean:
	$(RM) *.o *.a
