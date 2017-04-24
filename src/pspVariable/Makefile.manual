PSPDIR = ..
include $(PSPDIR)/make.inc

.SUFFIXES:
.SUFFIXES: .o .a .f90 .F90

all : pspNode.o psp_dList.o psp_zList.o pspListType.o pspMat.o pspVariable.o 

clean :
	rm -f *.a
	rm -f *.o
	rm -f *.mod

.F90.o :
	$(FORTRAN) $(OPTS) $(FPPFLAGS) -c $< -I$(PSPDIR) -o $@

.f90.o :
	$(FORTRAN) $(OPTS) -c $< -I$(PSPDIR) -o $@
