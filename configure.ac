# Copyright (C) 2016 Yann Pouillon <notifications@materialsevolution.es>
#
# This file is part of the pspBLAS software package. For license information,
# please see the COPYING file in the top-level directory of the source
# distribution.
#

#
# Configure script for the pspBLAS package
#

                    # ------------------------------------ #

#
# Autotools startup
#

# Init Autoconf
AC_PREREQ(2.69)
AC_INIT([pspBLAS], [0.0.0], [http://www.e-cam2020.eu:10080/ESL/omm/issues])
AC_REVISION([Autotools support for pspBLAS])
AC_CONFIG_AUX_DIR(config/gnu)
AC_CONFIG_MACRO_DIR([config/m4])
AC_CONFIG_SRCDIR([src/pspBLAS.F90])
_AC_SRCDIRS([.])

# Init Automake
AC_CANONICAL_TARGET
AM_INIT_AUTOMAKE([1.14 subdir-objects parallel-tests color-tests])
AM_CONFIG_HEADER([config.h])

# Generate timestamp
psp_timestamp=`date '+%Y%m%dT%H%M%S%z'`
AC_SUBST(psp_timestamp)

                    # ------------------------------------ #

#
# System utilities
#

# Check for common programs
AC_PROG_MAKE_SET
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_SED
AC_PROG_AWK
AC_PROG_GREP

                    # ------------------------------------ #

#
# Default settings
#

# Note for developers: you may edit these parameters to tune the behaviour
# of the build system. Please do it with care.

# Lists of options to manage (psp_*_enable_def variables)
#
# Note:
#   * switches (--enable-*) must be set to "yes" or "no"
#   * external dependencies (--with-*) must be set to "yes", "no", or "auto"
#
psp_cfg_switches="debug"
psp_cfg_extdeps="linalg mpi"

# Optional features
psp_debug_enable_def="yes"

# Linear algebra
# Note: default libraries must be empty because some compilers have
#       built-in support.
psp_linalg_enable_def="yes"
psp_linalg_incs_def=""
psp_linalg_libs_def="-lscalapack -lblacs -lblacsCinit -lblacsF77init -llapack -lblas"

# MPI
psp_mpi_enable_def="yes"
psp_mpi_incs_def=""
psp_mpi_libs_def=""

                    # ------------------------------------ #

#
# Check consistency of defaults
#

# Default triggers for optional features must be yes or no, and not empty
for psp_cfg_default in ${psp_cfg_switches}; do
  tmp_default=`eval echo \$\{psp_${psp_cfg_default}_enable_def\}`
  AC_MSG_CHECKING([default value of the ${psp_cfg_default} switch])
  AC_MSG_RESULT([${tmp_default}])
  if test "${tmp_default}" != "no" -a \
          "${tmp_default}" != "yes"; then
    AC_MSG_ERROR([invalid default value: psp_${psp_cfg_default}_enable_def = ${tmp_default}])
  fi
done

# Default triggers for external dependencies must be yes, no, or auto
for psp_cfg_default in ${psp_cfg_extdeps}; do
  tmp_default=`eval echo \$\{psp_${psp_cfg_default}_enable_def\}`
  AC_MSG_CHECKING([default detection behavior for ${psp_cfg_default}])
  AC_MSG_RESULT([${tmp_default}])
  if test "${tmp_default}" != "" -a \
          "${tmp_default}" != "auto" -a \
          "${tmp_default}" != "no" -a \
          "${tmp_default}" != "yes"; then
    AC_MSG_ERROR([invalid default value: psp_${psp_cfg_default}_enable_def = ${tmp_default}])
  fi
done

# Clean-up
unset tmp_default

                    # ------------------------------------ #

#
# Optional features
#

# Conditional compilation: --enable-* options (alphabetical order)
AC_ARG_ENABLE([debug],
  [AS_HELP_STRING([--enable-debug],
    [Enable extra debug information])],
  [psp_debug_enable="${enableval}"; psp_debug_type="yon"],
  [psp_debug_enable="${psp_debug_enable_def}"; psp_debug_type="def"])
AC_SUBST(enable_debug)

                    # ------------------------------------ #

#
# External dependencies
#

# MPI (mandatory for now)
AC_ARG_WITH([mpi],
  [AS_HELP_STRING([--with-mpi],
    [Install prefix of MPI (e.g. /usr/local). The default behaviour is to detect whether the specified compilers properly support MPI and to fall back to serial mode if not. You may use --with-mpi without argument to force MPI detection, in which case detection failures will result in errors, and --without-mpi to disable MPI support completely.])],
  [ if test "${withval}" = "no" -o "${withval}" = "yes"; then
      psp_mpi_enable="${withval}"
      psp_mpi_type="yon"
    else
      psp_mpi_enable="yes"
      psp_mpi_type="dir"
    fi],
  [psp_mpi_enable="${psp_mpi_enable_def}"; psp_mpi_type="def"])
AC_SUBST(with_mpi)

AC_ARG_VAR([MPICC], [MPI-aware C compiler (usually called 'mpicc'). Please note that this variable conflicts with the --with-mpi option.])
if test -n "${MPICC}"; then
  psp_mpi_enable="yes"
  psp_mpi_type="env"
fi
AC_ARG_VAR([MPIFC], [MPI-aware Fortran compiler (usually called 'mpif90'. Please note that this variable conflicts with the --with-mpi option.])
if test -n "${MPIFC}"; then
  psp_mpi_enable="yes"
  psp_mpi_type="env"
fi
AC_ARG_VAR([MPI_INCLUDES],
  [Include flags for MPI. Please note that this variable conflicts with the --with-mpi option.])
if test -n "${MPI_INCLUDES}"; then
  psp_mpi_enable="yes"
  psp_mpi_type="env"
  psp_mpi_incs="${MPI_INCLUDES}"
fi
AC_ARG_VAR([MPI_LIBS],
    [Library flags for MPI. Please note that this variable conflicts with the --with-mpi option.])
if test -n "${MPI_LIBS}"; then
  psp_mpi_enable="yes"
  psp_mpi_type="env"
  psp_mpi_libs="${MPI_LIBS}"
fi

# Linear algebra (mandatory)
AC_ARG_WITH([linalg],
  [AS_HELP_STRING([--with-linalg],
    [Install prefix of the linear algebra libraries (e.g. /usr/local). Since linear algebra is a mandatory dependency, --with-linalg without argument has no effect, while --without-linalg produces an error.])],
  [ if test "${withval}" = "no" -o "${withval}" = "yes"; then
      psp_linalg_enable="${withval}"
      psp_linalg_type="yon"
    else
      psp_linalg_enable="yes"
      psp_linalg_type="dir"
      psp_linalg_incs="-I${withval}/include"
      psp_linalg_libs="-L${withval}/lib"
    fi],
  [psp_linalg_enable="${psp_linalg_enable_def}"; psp_linalg_type="def";
    psp_linalg_incs=""; psp_linalg_libs=""])
AC_SUBST(with_linalg)

AC_ARG_VAR([LINALG_INCLUDES],
  [Include flags for linear algebra. Please note that this variable conflicts with the --with-linalg option.])
if test -n "${LINALG_INCLUDES}"; then
  psp_linalg_enable="yes"
  psp_linalg_type="env"
fi
AC_ARG_VAR([LINALG_LIBS],
    [Library flags for linear algebra. Please note that this variable conflicts with the --with-linalg option.])
if test -n "${LINALG_LIBS}"; then
  psp_linalg_enable="yes"
  psp_linalg_type="env"
fi

                    # ------------------------------------ #

#
# Check option consistency
#

# All --enable-* options must be yes or no
for psp_cfg_option in ${psp_cfg_switches}; do
  tmp_option=`eval echo \$\{enable_${psp_cfg_option}\}`
  if test "${tmp_option}" != "" -a \
          "${tmp_option}" != "no" -a \
          "${tmp_option}" != "yes"; then
    AC_MSG_ERROR([--enable-${psp_cfg_option} must be "yes" or "no"])
  fi
done

# Linear algebra must be enabled
if test "${psp_linalg_enable}" != "yes"; then
  AC_MSG_ERROR([linear algebra support has been disabled
                  while it is mandatory])
fi

# All MPI environment variables must be provided at once
tmp_mpi_vars="${MPICC}${MPIFC}${MPI_INCLUDES}${MPI_LIBS}"
if test -n "${tmp_mpi_vars}" -a \
        \( -z "${MPICC}" -o -z "${MPIFC}" \); then
  AC_MSG_ERROR([MPICC and MPIFC must both be set])
fi

# Environment variables conflict with --with-* options
# Note: for linear algebra, only LINALG_LIBS is relevant
if test -n "${tmp_mpi_vars}" -a -n "${with_mpi}"; then
  AC_MSG_ERROR([conflicting option settings
                  Please use MPICC + MPIFC or --with-mpi,
                  not both.])
fi
if test -n "${LINALG_LIBS}" -a -n "${with_linalg}"; then
  AC_MSG_ERROR([conflicting option settings
                  Please use LINALG_LIBS or --with-linalg,
                  not both.])
fi

# When using environment variables, triggers must be set to yes
if test -n "${tmp_mpi_vars}"; then
  psp_mpi_enable="yes"
  psp_mpi_type="env"
fi

# Clean-up
unset tmp_option
unset tmp_mpi_vars

                    # ------------------------------------ #

#
# Substitute build-system variables
#

# Defaults
AC_SUBST(psp_debug_enable_def)
AC_SUBST(psp_mpi_enable_def)
AC_SUBST(psp_mpi_incs_def)
AC_SUBST(psp_mpi_libs_def)
AC_SUBST(psp_linalg_enable_def)
AC_SUBST(psp_linalg_incs_def)
AC_SUBST(psp_linalg_libs_def)

# Triggers
AC_SUBST(psp_debug_enable)
AC_SUBST(psp_mpi_enable)
AC_SUBST(psp_linalg_enable)

# Initialization types
AC_SUBST(psp_debug_type)
AC_SUBST(psp_mpi_type)
AC_SUBST(psp_linalg_type)

# Internal MPI parameters
AC_SUBST(psp_mpi_cc_set)
AC_SUBST(psp_mpi_cc_wrap)
AC_SUBST(psp_mpi_fc_set)
AC_SUBST(psp_mpi_fc_wrap)

# Build flags for external dependencies
AC_SUBST(psp_mpi_incs)
AC_SUBST(psp_mpi_libs)
AC_SUBST(psp_linalg_incs)
AC_SUBST(psp_linalg_libs)
AC_SUBST(psp_linalg_has_lapack)
AC_SUBST(psp_linalg_has_scalapack)

                    # ------------------------------------ #

#
# Architecture setup
#

# Init MPI parameters, possibly setting and/or wrapping compilers
PSP_MPI_INIT

                    # ------------------------------------ #

#
# C language support
#

# Init C compiler and preprocessor
AC_PROG_CC([mpicc mpiicc icc gcc cc])
AC_PROG_CPP

# C compiler peculiarities (for Libtool)
AM_PROG_CC_C_O

                    # ------------------------------------ #

#
# Fortran language support
#

# Look for the Fortran compiler
if test "${FC}" != "" -a ! -x "${FC}"; then
  psp_fc_probe=`echo "${FC}" | sed -e 's/ .*//'`
  if test ! -x "${psp_fc_probe}"; then
    AC_PATH_PROG([psp_fc_path],[${psp_fc_probe}])
    if test "${psp_fc_path}" = ""; then
      AC_MSG_ERROR([could not run Fortran compiler "${FC}"])
    fi
  fi
fi
AC_PROG_FC([mpif90 mpiifort ifort gfortran g95 xlf95 xlf90 f95 fort ifc efc pgfortran pgf95 lf95 ftn nagfor xlf90 f90 pgf90 pghpf epcf90])

# Fail if no Fortran compiler is available
if test "${FC}" = ""; then
  AC_MSG_ERROR([no Fortran compiler available])
fi

# Set default file extensions
PSP_FC_EXTENSIONS

# Get module file case
PSP_FC_MOD_CASE

                    # ------------------------------------ #

#
# Libtool configuration
#

# Init Libtool (must be done once compilers are fully set)
LT_INIT
LT_PREREQ([2.4.2])
LT_LANG([Fortran])
LTOBJEXT="lo"
AC_SUBST(LTOBJEXT)

                    # ------------------------------------ #

#
# Multicore architecture support
#

# Look for MPI
if test "${psp_mpi_enable}" != "no"; then
  PSP_MPI_DETECT
  if test "${psp_mpi_ok}" = "yes"; then
    AC_DEFINE([HAVE_MPI], 1, [Define to 1 to enable MPI support.])
    psp_mpi_enable="yes"
  else
    if test "${psp_mpi_enable}" = "auto"; then
      AC_MSG_NOTICE([disabling MPI support])
      psp_mpi_enable="no"
    else
      AC_MSG_FAILURE([MPI support is broken - please check your configuration])
    fi
  fi
fi

# From now on, the MPI trigger must be yes or no
if test "${psp_mpi_enable}" != "no" -a "${psp_mpi_enable}" != "yes"; then
  AC_MSG_ERROR([unexpected MPI trigger value: '${psp_mpi_enable}'])
fi

# Let the test suite know about mpirun
TESTS_CMD="${MPIRUN}"
AC_SUBST(TESTS_CMD)

                    # ------------------------------------ #

#
# External dependencies
#

# Set linear algebra includes
if test "${psp_linalg_type}" != "dir"; then
  if test -n "${LINALG_INCLUDES}"; then
    psp_linalg_incs="${LINALG_INCLUDES}"
  else
    psp_linalg_incs="${psp_linalg_incs_def}"
  fi
fi

# Set linear algebra libraries
if test "${psp_linalg_type}" != "dir"; then
  if test -n "${LINALG_LIBS}"; then
    psp_linalg_libs="${LINALG_LIBS}"
  else
    psp_linalg_libs="${psp_linalg_libs_def}"
  fi
fi

# Look for linear algebra
AC_MSG_CHECKING([how linear algebra parameters have been set])
AC_MSG_RESULT([${psp_linalg_type}])
AC_MSG_CHECKING([for linear algebra include flags])
AC_MSG_RESULT([${psp_linalg_incs}])
AC_MSG_CHECKING([for linear algebra library flags])
AC_MSG_RESULT([${psp_linalg_libs}])
PSP_LINALG_DETECT
if test "${psp_linalg_ok}" != "yes" -a "${psp_linalg_type}" = "def"; then
  AC_MSG_NOTICE([retrying with generic linear algebra library flags (Netlib)])
  LIBS="-lscalapack -lblacs -lblacsCinit -lblacsF77init -llapack -lblas ${LIBS}"
  PSP_LINALG_DETECT
fi
if test "${psp_linalg_ok}" = "yes"; then
  CPPFLAGS="${CPPFLAGS} ${psp_linalg_incs}"
  LIBS="${psp_linalg_libs} ${LIBS}"
else
  AC_MSG_FAILURE([linear algebra is missing or incomplete])
fi

                    # ------------------------------------ #

#
# Adjust build parameters according to configure options
#

# Reset option counter
psp_cpp_options=""

# ISO C Bindings are always available nowadays
AC_DEFINE([CBIND], 1, [Define to 1 to use ISO C Bindings])
psp_cpp_options="${psp_cpp_options} CBIND"

# Conditional compilation (alphabetical order)
if test "${psp_debug_enable}" = "yes"; then
  AC_DEFINE([DEBUG_MODE], 1,
    [Define to 1 if you want to enable enhanced debugging information.])
  psp_cpp_options="${psp_cpp_options} DEBUG_MODE"
fi

# External dependencies (alphabetical order)
if test "${psp_linalg_has_lapack}" = "yes"; then
  AC_DEFINE([HAVE_LAPACK], 1, [Define to 1 if you have LAPACK libraries available.])
  psp_cpp_options="${psp_cpp_options} HAVE_LAPACK"
fi
if test "${psp_linalg_has_scalapack}" = "yes"; then
  AC_DEFINE([HAVE_SCALAPACK], 1, [Define to 1 if you have ScaLAPACK libraries available.])
  psp_cpp_options="${psp_cpp_options} HAVE_SCALAPACK"
fi

                    # ------------------------------------ #

# Write down YAML configuration now (allows full report with deferred errors)
AC_OUTPUT([config/data/pspBLAS-config.yaml])

# Report configuration
AC_MSG_NOTICE([])
AC_MSG_NOTICE([Final build parameters])
AC_MSG_NOTICE([----------------------])
AC_MSG_NOTICE([])
AC_MSG_NOTICE([TSTAMP     = ${psp_timestamp}])
AC_MSG_NOTICE([])
AC_MSG_NOTICE([DEBUG    = ${psp_debug_enable} (init: ${psp_debug_type})])
AC_MSG_NOTICE([MPI      = ${psp_mpi_enable} (init: ${psp_mpi_type})])
AC_MSG_NOTICE([LINALG   = ${psp_linalg_enable} (init: ${psp_linalg_type})])
AC_MSG_NOTICE([])
AC_MSG_NOTICE([CPP      = ${CPP}])
AC_MSG_NOTICE([CPPFLAGS = ${CPPFLAGS}])
AC_MSG_NOTICE([FC       = ${CC}])
AC_MSG_NOTICE([MPIFC    = ${MPICC}])
AC_MSG_NOTICE([FC       = ${FC}])
AC_MSG_NOTICE([MPIFC    = ${MPIFC}])
AC_MSG_NOTICE([FCFLAGS  = ${FCFLAGS}])
AC_MSG_NOTICE([LDFLAGS  = ${LDFLAGS}])
AC_MSG_NOTICE([LIBS     = ${LIBS}])
AC_MSG_NOTICE([OPTIONS  = ${psp_cpp_options}])
AC_MSG_NOTICE([])

                    # ------------------------------------ #

#
# Output configuration
#

AC_CONFIG_FILES([
  config/data/pspBLAS.pc
  Makefile
  doc/Makefile
  src/Makefile
  test/Makefile
])

AC_OUTPUT
