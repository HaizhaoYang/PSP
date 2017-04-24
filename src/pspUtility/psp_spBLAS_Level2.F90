!************************************************************************!
!   Copyright (c) 2015-2017, Haizhao Yang                                !
!   All rights reserved.                                                 !
!                                                                        !
!   This file is part of Elemental and is under the BSD 2-Clause License,! 
!   which can be found in the LICENSE file in the root directory, or at  !
!   http://opensource.org/licenses/BSD-2-Clause                          !
!************************************************************************!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

MODULE psp_spBLAS_Level2
  use pspVariable
  use pspBasicTool
  use pspListTool
  use psp_spBLAS_Level1

  ! This module contains sequential sparse BLAS

#ifdef HAVE_MPI
  include 'mpif.h'
#endif

  private

  !**** PARAMS ************************************!

  integer, parameter :: dp=selected_real_kind(15,300)

  complex(dp), parameter :: cmplx_1=(1.0_dp,0.0_dp)
  complex(dp), parameter :: cmplx_i=(0.0_dp,1.0_dp)
  complex(dp), parameter :: cmplx_0=(0.0_dp,0.0_dp)

  !**** INTERFACES ********************************!



contains

END MODULE psp_spBLAS_Level2
