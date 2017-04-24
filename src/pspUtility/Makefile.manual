PSPDIR = ..
include $(PSPDIR)/make.inc

.SUFFIXES:
.SUFFIXES: .o .a .f90 .F90

all : pspBasicTool.o pspListTool.o psp_spBLAS_Level1.o psp_spBLAS_Level2.o psp_spBLAS_Level3.o psp_spBLAS.o pspUtility.o 

clean :
	rm -f *.a
	rm -f *.o
	rm -f *.mod

.F90.o :
	$(FORTRAN) $(OPTS) $(FPPFLAGS) -c $< -I$(PSPDIR) -o $@

.f90.o :
	$(FORTRAN) $(OPTS) -c $< -I$(PSPDIR) -o $@
