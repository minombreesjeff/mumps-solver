#
#  This file is part of MUMPS 4.8.3, built on Wed Sep 24 17:11:30 UTC 2008
#
all: libmpiseq

.PHONY: all libmpiseq clean

include ../Makefile.inc

libmpiseq: libmpiseq$(PLAT).a

libmpiseq$(PLAT).a: mpi.o mpic.o elapse.o
	$(AR) $@ mpi.o mpic.o elapse.o
	$(RANLIB) $@
.f.o:
	$(FC) $(OPTF) -c $*.f
.c.o:
	$(CC) $(OPTC) $(CDEFS) -I. -c $*.c

clean:
	$(RM) *.o *.a
