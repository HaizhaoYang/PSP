# -*- Automake -*-
#
# Makefile for the pspBLAS package
#
# Copyright (C) 2016 Yann Pouillon
#
# This file is part of the pspBLAS software package. For license information,
# please see the COPYING file in the top-level directory of the source
# distribution.
#

                    # ------------------------------------ #

#
# Main source files
#

# Common source files
psp_srcs = \
  pspBLAS.F90 \
  pspLevel1/pspLevel1.F90 \
  pspLevel2/pspLevel2.F90 \
  pspLevel3/pspLevel3.F90 \
  pspLevel3/pspMspm.F90 \
  pspLevel3/pspMspm_nt.F90 \
  pspLevel3/pspMspm_tt.F90 \
  pspLevel3/pspSpmm_nn.F90 \
  pspLevel3/pspSpmm_tn.F90 \
  pspLevel3/pspSpmSpm.F90 \
  pspLevel3/pspSpmSpm_nt.F90 \
  pspLevel3/pspSpmSpm_tt.F90 \
  pspLevel3/pspGemm.F90 \
  pspLevel3/pspMatSum.F90 \
  pspLevel3/pspMspm_nn.F90 \
  pspLevel3/pspMspm_tn.F90 \
  pspLevel3/pspSpmm.F90 \
  pspLevel3/pspSpmm_nt.F90 \
  pspLevel3/pspSpmm_tt.F90 \
  pspLevel3/pspSpmSpm_nn.F90 \
  pspLevel3/pspSpmSpm_tn.F90 \
  pspMPI/pspMPI.F90 \
  pspUtility/pspBasicTool.F90 \
  pspUtility/pspListTool.F90 \
  pspUtility/psp_spBLAS.F90 \
  pspUtility/psp_spBLAS_Level1.F90 \
  pspUtility/psp_spBLAS_Level2.F90 \
  pspUtility/psp_spBLAS_Level3.F90 \
  pspUtility/pspUtility.F90 \
  pspVariable/psp_dList.F90 \
  pspVariable/pspListType.F90 \
  pspVariable/pspMat.F90 \
  pspVariable/pspNode.F90 \
  pspVariable/pspVariable.F90 \
  pspVariable/psp_zList.F90

# Fortran modules
# TODO: support both lower-case and upper-case file names
psp_f03_mods = \
  psplevel3.$(MODEXT) \
  pspmspm.$(MODEXT) \
  pspmspm_nt.$(MODEXT) \
  pspmspm_tt.$(MODEXT) \
  pspspmm_nn.$(MODEXT) \
  pspspmm_tn.$(MODEXT) \
  pspspmspm.$(MODEXT) \
  pspspmspm_nt.$(MODEXT) \
  pspspmspm_tt.$(MODEXT) \
  pspgemm.$(MODEXT) \
  pspmatsum.$(MODEXT) \
  pspmspm_nn.$(MODEXT) \
  pspmspm_tn.$(MODEXT) \
  pspspmm.$(MODEXT) \
  pspspmm_nt.$(MODEXT) \
  pspspmm_tt.$(MODEXT) \
  pspspmspm_nn.$(MODEXT) \
  pspspmspm_tn.$(MODEXT) \
  pspmpi.$(MODEXT) \
  psplevel2.$(MODEXT) \
  pspblas.$(MODEXT) \
  pspbasictool.$(MODEXT) \
  psplisttool.$(MODEXT) \
  psp_spblas.$(MODEXT) \
  psp_spblas_level1.$(MODEXT) \
  psp_spblas_level2.$(MODEXT) \
  psp_spblas_level3.$(MODEXT) \
  psputility.$(MODEXT) \
  psp_dlist.$(MODEXT) \
  psplisttype.$(MODEXT) \
  pspmat.$(MODEXT) \
  pspnode.$(MODEXT) \
  pspvariable.$(MODEXT) \
  psp_zlist.$(MODEXT) \
  psplevel1.$(MODEXT)

# Libraries to install
lib_LTLIBRARIES = libpspBLAS.la

libpspBLAS_la_SOURCES = $(psp_srcs)
libpspBLAS_la_LDFLAGS = -version-info 0:0:0

# Fortran modules to install
f03moddir = $(includedir)
install-data-local:
	$(INSTALL) -d -m 755 $(DESTDIR)$(f03moddir)
	$(INSTALL) -m 644 $(psp_f03_mods) $(DESTDIR)$(f03moddir)

uninstall-local:
	cd $(DESTDIR)$(f03moddir) && rm -f $(psp_f03_mods)

# Local cleaning
clean-local:
	rm -f *.mod *.MOD *.obj

# Explicit dependencies
include $(srcdir)/dependencies.mf

# Files to mark explicitly for distribution
EXTRA_DIST = dependencies.mf
