emacmSHELL     = /bin/sh
FORTRAN   = mpiifort
OPTS      = -O3
FPPFLAGS  = -DCONV -DHAVE_PSPBLAS -DHAVE_MPI -DHAVE_LAPACK -DHAVE_SCALAPACK -DHAVE_MKL
ARCH      = ar
ARCHFLAGS = cr
RANLIB    = ranlib

PSPLIBNAME = libPSP.a
PSPBUILDPATH = $(ROOTDIR)/build_PSP

SLAPLIBPATH =

LINALG_INCLUDES = -I${MKLROOT}/include/intel64/lp64
#LINALG_LIBS = -lscalapack -llapack -lblas
#LINALG_LIBS = -llapack -lblas $(SLAPLIBPATH)
LINALG_LIBS = -L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread