# -*- Autoconf -*-
#
# M4 macros for pspBLAS
#
# Copyright (C) 2016 Yann Pouillon
#
# This file is part of the pspBLAS software package. For license information,
# please see the COPYING file in the top-level directory of the pspBLAS source
# distribution.
#

#
# Linear algebra support
#



# PSP_LINALG_DETECT()
# -------------------
#
# Checks that the selected linear algebra libraries properly work.
#
AC_DEFUN([PSP_LINALG_DETECT],[
  dnl Init
  psp_linalg_has_lapack="unknown"
  psp_linalg_has_scalapack="unknown"
  psp_linalg_ok="unknown"

  dnl Prepare environment
  saved_CPPFLAGS="${CPPFLAGS}"
  saved_LIBS="${LIBS}"
  CPPFLAGS="${CPPFLAGS} ${psp_linalg_incs}"
  LIBS="${psp_linalg_libs} ${LIBS}"

  dnl Check BLAS and LAPACK routines
  AC_LANG_PUSH([Fortran])
  AC_MSG_CHECKING([whether linear algebra libraries work])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      call zgemm
      call zhpev
    ]])], [psp_linalg_ok="yes"; psp_linalg_has_lapack="yes"], [psp_linalg_ok="no"])
  AC_MSG_RESULT([${psp_linalg_ok}])
  AC_LANG_POP([Fortran])

  dnl Check ScaLAPACK routines
  AC_LANG_PUSH([Fortran])
  AC_MSG_CHECKING([whether linear algebra libraries have ScaLAPACK])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      call pzheevx
    ]])], [psp_linalg_has_scalapack="yes"], [psp_linalg_has_scalapack="no"])
  AC_MSG_RESULT([${psp_linalg_ok}])
  AC_LANG_POP([Fortran])

  dnl Restore environment
  CPPFLAGS="${saved_CPPFLAGS}"
  LIBS="${saved_LIBS}"
]) # PSP_LINALG_DETECT
