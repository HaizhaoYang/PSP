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
# Fortran compilers support
#



# PSP_FC_EXTENSIONS()
# -------------------
#
# Sets the default extensions of Fortran source files and modules,
# whenever possible.
#
AC_DEFUN([PSP_FC_EXTENSIONS], [
  dnl Set Fortran module extension
  AX_F90_MODULE_EXTENSION
  if test "${ax_cv_f90_modext}" != ""; then
    MODEXT="${ax_cv_f90_modext}"
  else
    MODEXT="mod"
    AC_MSG_NOTICE([setting Fortran module extension to ".${MODEXT}"])
  fi
  AC_SUBST(MODEXT)

  dnl Change the default Fortran extension for tests
  AC_FC_SRCEXT(F90, [psp_fc_src_ok="yes"], [psp_fc_src_ok="no"])
  if test "${psp_fc_src_ok}" != "yes"; then
    AC_MSG_WARN([Fortran file extension could not be changed])
    AC_MSG_WARN([some advanced Fortran tests may fail])
  fi
]) # PSP_FC_EXTENSIONS



# PSP_FC_MOD_CASE()
# -----------------
#
# Checks whether the Fortran compiler creates upper-case or lower-case
# module files.
#
AC_DEFUN([PSP_FC_MOD_CASE],[
  AC_REQUIRE([PSP_FC_EXTENSIONS])

  dnl Init
  fc_mod_lowercase="yes"
  fc_mod_uppercase="no"
  AC_MSG_NOTICE([determining Fortran module case])

  dnl Compile a dummy module
  AC_LANG_PUSH([Fortran])
  AC_COMPILE_IFELSE([[
    module conftest
    end module conftest
  ]], [], [AC_MSG_FAILURE([unable to compile a simple Fortran module])])
  AC_LANG_POP([Fortran])

  dnl Check module file existence
  if test -f "CONFTEST.${MODEXT}"; then
    fc_mod_lowercase="no"
    fc_mod_uppercase="yes"
  elif test ! -f "conftest.${MODEXT}"; then
    AC_MSG_WARN([conftest.${MODEXT} Fortran module could not be found])
  fi

  dnl Clean-up the mess
  rm -f CONFTEST.${MODEXT} conftest.${MODEXT}

  dnl Output final outcome
  AC_MSG_CHECKING([whether Fortran modules are upper-case])
  AC_MSG_RESULT([${fc_mod_uppercase}])
]) # PSP_FC_MOD_CASE
