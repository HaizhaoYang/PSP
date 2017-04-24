# -*- Autoconf -*-
#
# M4 macros for pspBLAS
#
# Copyright (C) 2015 Yann Pouillon
#
# This file is part of the pspBLAS software package. For license information,
# please see the COPYING file in the top-level directory of the pspBLAS source
# distribution.
#

#
# MPI support
#



# PSP_MPI_DETECT()
# ----------------
#
# Checks whether the configured MPI implementation is working.
#
AC_DEFUN([PSP_MPI_DETECT], [
  dnl Init
  psp_mpi_ok="unknown"

  dnl Display current MPI status
  AC_MSG_CHECKING([how MPI parameters have been set])
  AC_MSG_RESULT([${psp_mpi_type}])
  AC_MSG_CHECKING([whether the MPI C compiler is set])
  AC_MSG_RESULT([${psp_mpi_cc_set}])
  AC_MSG_CHECKING([whether the MPI C compiler is wrapped])
  AC_MSG_RESULT([${psp_mpi_cc_wrap}])
  AC_MSG_CHECKING([whether the MPI Fortran compiler is set])
  AC_MSG_RESULT([${psp_mpi_fc_set}])
  AC_MSG_CHECKING([whether the MPI Fortran compiler is wrapped])
  AC_MSG_RESULT([${psp_mpi_fc_wrap}])

  dnl Warn if serial component of wrapped compilers supports MPI
  if test "${psp_mpi_cc_wrap}" = "yes"; then
    AC_MSG_NOTICE([validating that '${psp_sercc}' is indeed serial])
    AC_MSG_NOTICE([please ignore possible warnings about mpi.h not found])
    _PSP_MPI_CHECK_CC([${psp_sercc}])
    if test "${psp_mpi_cc_ok}" = "yes"; then
      AC_MSG_WARN([the serial C compiler is MPI-aware
                    Your current configuration is probably ill-defined.
                    The build will likely fail.])
      sleep 5
    fi
  fi
  if test "${psp_mpi_fc_wrap}" = "yes"; then
    AC_MSG_NOTICE([validating that '${psp_serfc}' is indeed serial])
    _PSP_MPI_CHECK_FC([${psp_serfc}])
    if test "${psp_mpi_fc_ok}" = "yes"; then
      AC_MSG_WARN([the serial Fortran compiler is MPI-aware
                    Your current configuration is probably ill-defined.
                    The build will likely fail.])
      sleep 5
    fi
  fi

  dnl Test MPI compilers
  _PSP_MPI_CHECK_CC([${CC}])
  if test "${psp_mpi_cc_ok}" = "yes"; then
    _PSP_MPI_CHECK_FC([${FC}])
  fi

  dnl Look for mpirun
  dnl FIXME: hard-coded command-line options
  if test "${MPIRUN}" = ""; then
    AC_CHECK_PROGS([MPIRUN], [mpirun mpiexec])
  fi
  if test "${MPIRUN}" != ""; then
    MPIRUN="${MPIRUN} -np 4"
  fi

  dnl Take final decision
  AC_MSG_CHECKING([whether we have a full MPI support])
  if test "${psp_mpi_cc_ok}" = "yes" -a \
          "${psp_mpi_fc_ok}" = "yes"; then
    psp_mpi_ok="yes"
  else
    psp_mpi_ok="no"
  fi
  AC_MSG_RESULT([${psp_mpi_ok}])
]) # PSP_MPI_DETECT



# PSP_MPI_INIT()
# --------------
#
# Initializes MPI parameters.
#
AC_DEFUN([PSP_MPI_INIT], [
  if test "${psp_mpi_enable}" != "no"; then
    AC_MSG_CHECKING([how MPI parameters have been set])
    AC_MSG_RESULT([${psp_mpi_type}])
    if test "${psp_mpi_type}" = "env"; then
      _AC_SRCDIRS(["."])
    fi
    _PSP_MPI_INIT_CC
    _PSP_MPI_INIT_FC
  fi
]) # PSP_MPI_INIT



                    ########################################



# _PSP_MPI_CHECK_CC(CC)
# ---------------------
#
# Check whether the MPI C compiler is working.
#
AC_DEFUN([_PSP_MPI_CHECK_CC], [
  dnl Init
  psp_mpi_cc_ok="unknown"
  psp_mpi_cc_has_funs="unknown"
  psp_mpi_cc_has_hdrs="unknown"

  dnl Prepare environment
  psp_saved_CC="${CC}"
  psp_saved_CC="${CC}"
  CC="$1"
  tmp_mpi_header=mpi.h
  tmp_mpi_cache=AS_TR_SH([ac_cv_header_${tmp_mpi_header}])
  ${as_unset} ${tmp_mpi_cache}

  dnl Look for C includes
  AC_LANG_PUSH([C])
  AC_CHECK_HEADERS([mpi.h],
    [psp_mpi_cc_has_hdrs="yes"], [psp_mpi_cc_has_hdrs="no"])
  AC_LANG_POP([C])

  dnl Look for C functions
  if test "${psp_mpi_cc_has_hdrs}" = "yes"; then
    AC_CHECK_FUNC([MPI_Init], [psp_mpi_cc_has_funs="yes"],
      [psp_mpi_cc_has_funs="no"])
  fi

  dnl Validate C support
  AC_MSG_CHECKING([whether the MPI C compiler works])
  if test "${psp_mpi_cc_has_funs}" = "yes" -a \
          "${psp_mpi_cc_has_hdrs}" = "yes"; then
    psp_mpi_cc_ok="yes"
  else
    psp_mpi_cc_ok="no"
  fi
  AC_MSG_RESULT([${psp_mpi_cc_ok}])

  dnl Restore environment
  CC="${psp_saved_CC}"
  unset tmp_mpi_cache
  unset tmp_mpi_header
]) # _PSP_MPI_CHECK_CC



# _PSP_MPI_CHECK_FC(FC)
# ---------------------
#
# Check whether the MPI Fortran compiler is working.
#
AC_DEFUN([_PSP_MPI_CHECK_FC], [
  dnl Init
  psp_mpi_fc_ok="unknown"
  psp_mpi_fc_has_funs="unknown"
  psp_mpi_fc_has_mods="unknown"

  dnl Prepare environment
  psp_saved_FC="${FC}"
  FC="$1"

  dnl Look for Fortran modules
  AC_LANG_PUSH([Fortran])
  AC_MSG_CHECKING([for a Fortran MPI module])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      use mpi
    ]])], [psp_mpi_fc_has_mods="yes"], [psp_mpi_fc_has_mods="no"])
  AC_MSG_RESULT([${psp_mpi_fc_has_mods}])
  AC_LANG_POP([Fortran])

  dnl Look for Fortran functions
  if test "${psp_mpi_fc_has_mods}" = "yes"; then
    AC_LANG_PUSH([Fortran])
    AC_MSG_CHECKING([for a Fortran MPI_Init])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
        use mpi
        integer :: ierr
        call mpi_init(ierr)
      ]])], [psp_mpi_fc_has_funs="yes"], [psp_mpi_fc_has_funs="no"])
    AC_MSG_RESULT([${psp_mpi_fc_has_funs}])
    AC_LANG_POP([Fortran])
  fi

  dnl Validate Fortran support
  AC_MSG_CHECKING([whether the MPI Fortran compiler works])
  if test "${psp_mpi_fc_has_funs}" = "yes" -a \
          "${psp_mpi_fc_has_mods}" = "yes"; then
    psp_mpi_fc_ok="yes"
  else
    psp_mpi_fc_ok="no"
  fi
  AC_MSG_RESULT([${psp_mpi_fc_ok}])

  dnl Restore environment
  FC="${psp_saved_FC}"
]) # _PSP_MPI_CHECK_FC



# _PSP_MPI_INIT_CC()
# ------------------
#
# Initializes MPI parameters related to the C compiler.
#
AC_DEFUN([_PSP_MPI_INIT_CC], [
  dnl Init
  psp_sercc="${CC}"
  psp_mpicc=""
  psp_mpi_cc_set="no"
  psp_mpi_cc_wrap="unknown"

  dnl Look for a MPI C compiler
  case "${psp_mpi_type}" in

    def)
      psp_mpi_cc_wrap="no"
      ;;

    dir)
      psp_mpicc="${with_mpi}/bin/mpicc"
      if test -x "${psp_mpicc}"; then
        AC_MSG_CHECKING([for an executable MPI C compiler])
        AC_MSG_RESULT([${psp_mpicc}])
        if test "${psp_sercc}" = ""; then
          AC_MSG_NOTICE([setting CC to '${psp_mpicc}'])
          CC="${psp_mpicc}"
          psp_mpi_cc_set="yes"
          psp_mpi_cc_wrap="no"
        else
          psp_mpi_cc_wrap="yes"
        fi
      else
        AC_MSG_ERROR([MPI C compiler not found in ${with_mpi}/bin])
      fi
      ;;

    env|yon)
      if test -n "${MPICC}"; then
        psp_mpicc="${MPICC}"
      else
        AC_CHECK_PROGS([psp_mpicc], [mpicc])
      fi
      if test -n "${psp_sercc}" -a -n "${psp_mpicc}"; then
        psp_mpi_cc_wrap="yes"
      elif test -n "${psp_mpicc}"; then
        AC_MSG_NOTICE([setting CC to '${psp_mpicc}'])
        CC="${psp_mpicc}"
        psp_mpi_cc_set="yes"
        psp_mpi_cc_wrap="no"
      fi
      ;;

  esac

  if test "${psp_mpi_cc_wrap}" = "yes"; then
    _PSP_MPI_CREATE_WRAPPER([CC], [${psp_sercc}], [${psp_mpicc}])
    psp_mpi_cc_set="yes"
  fi
]) # _PSP_MPI_INIT_CC



# _PSP_MPI_INIT_FC()
# ------------------
#
# Initializes MPI parameters related to the Fortran compiler.
#
AC_DEFUN([_PSP_MPI_INIT_FC], [
  dnl Init
  psp_serfc="${FC}"
  psp_mpifc=""
  psp_mpi_fc_set="no"
  psp_mpi_fc_wrap="unknown"

  dnl Look for a MPI Fortran compiler
  case "${psp_mpi_type}" in

    def)
      psp_mpi_fc_wrap="no"
      ;;

    dir)
      psp_mpifc="${with_mpi}/bin/mpif90"
      if test -x "${psp_mpifc}"; then
        AC_MSG_CHECKING([for an executable MPI Fortran compiler])
        AC_MSG_RESULT([${psp_mpifc}])
        if test "${psp_serfc}" = ""; then
          AC_MSG_NOTICE([setting FC to '${psp_mpifc}'])
          FC="${psp_mpifc}"
          psp_mpi_fc_set="yes"
          psp_mpi_fc_wrap="no"
        else
          psp_mpi_fc_wrap="yes"
        fi
      else
        AC_MSG_ERROR([MPI Fortran compiler not found in ${with_mpi}/bin])
      fi
      ;;

    env|yon)
      if test -n "${MPIFC}"; then
        psp_mpifc="${MPIFC}"
      else
        AC_CHECK_PROGS([psp_mpifc], [mpif90 mpif95])
      fi
      if test -n "${psp_serfc}" -a -n "${psp_mpifc}"; then
        psp_mpi_fc_wrap="yes"
      elif test -n "${psp_mpifc}"; then
        AC_MSG_NOTICE([setting FC to '${psp_mpifc}'])
        FC="${psp_mpifc}"
        psp_mpi_fc_set="yes"
        psp_mpi_fc_wrap="no"
      fi
      ;;

  esac

  if test "${psp_mpi_fc_wrap}" = "yes"; then
    _PSP_MPI_CREATE_WRAPPER([FC], [${psp_serfc}], [${psp_mpifc}])
    psp_mpi_fc_set="yes"
  fi
]) # _PSP_MPI_INIT_FC



# _PSP_MPI_CREATE_WRAPPER(COMPILER_TYPE, SERIAL_COMPILER, MPI_COMPILER)
# ---------------------------------------------------------------------
#
# Creates a wrapper for MPI compilers when they can be configured to
# accept different serial compilers.
#
# Note: it is impossible to set two compiler levels with the Autotools,
#       because Automake requires CC, CXX, and FC to be set to
#       the actual compilers.
#
AC_DEFUN([_PSP_MPI_CREATE_WRAPPER], [
  dnl Init
  tmp_comp_name=`echo "$1" | sed -e 's/.*/\L&/'`
  ${MKDIR_P} config/wrappers

  dnl Create file
  cat >config/wrappers/wrap-mpi${tmp_comp_name} <<EOF
#!/bin/sh

$1="$2"
export $1

$3 \[$]{*}
EOF

  dnl Fix permissions
  chmod u+x config/wrappers/wrap-mpi${tmp_comp_name}

  dnl Overwrite compiler setting
  eval tmp_wrapper_path="${ac_abs_top_builddir}/config/wrappers/wrap-mpi${tmp_comp_name}"
  tmp_wrapper_name=`basename "${tmp_wrapper_path}"`
  AC_MSG_NOTICE([wrapping serial and MPI compilers into ${tmp_wrapper_name}])
  $1="${tmp_wrapper_path}"

  dnl Clean-up
  unset tmp_comp_name
  unset tmp_wrapper_name
  unset tmp_wrapper_path
]) # _PSP_MPI_CREATE_WRAPPER
