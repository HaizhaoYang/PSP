!************************************************************************!
!   Copyright (c) 2015-2017, Haizhao Yang                                !
!   All rights reserved.                                                 !
!                                                                        !
!   This file is part of PSP and is under the BSD 2-Clause License,! 
!   which can be found in the LICENSE file in the root directory, or at  !
!   http://opensource.org/licenses/BSD-2-Clause                          !
!************************************************************************!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

MODULE pspMspm
  use pspVariable
  use pspUtility
  use pspMPI
  use pspLevel1
  use pspLevel2
use pspMspm_nn, only: psp_gemspm_nn
use pspMspm_nt, only: psp_gemspm_nt
use pspMspm_tn, only: psp_gemspm_tn
use pspMspm_tt, only: psp_gemspm_tt

#ifdef HAVE_MPI
  include 'mpif.h'
#endif

  private

  !**** PARAMS ************************************!

  integer, parameter :: dp=selected_real_kind(15,300)

  complex(dp), parameter :: cmplx_1=(1.0_dp,0.0_dp)
  complex(dp), parameter :: cmplx_i=(0.0_dp,1.0_dp)
  complex(dp), parameter :: cmplx_0=(0.0_dp,0.0_dp)

  !**********************************************!

  integer, external :: numroc ! it is a function to compute local size

  !**** INTERFACES ********************************!

  interface psp_gemspm
     module procedure psp_dgemspm
     module procedure psp_zgemspm
  end interface psp_gemspm

  interface die
     module procedure die
  end interface die

  public :: psp_gemspm

contains

  subroutine psp_dgemspm(M,N,K,A,opA,B,opB,C,alpha,beta)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) ::            M, N, K ! Globally, A of size M by K, B of size K by N
    character(1), intent(in) :: opA ! form of op(A): 'n/N' for A, 't/T/c/C' for A^T
    character(1), intent(in) :: opB ! form of op(B)
    real(dp), intent(in)  ::   alpha, beta ! scalar
    real(dp), intent(in) :: A(:,:)

    !**** INOUT ***********************************!

    type(psp_matrix_spm), intent(in) ::   B ! matrix B
    real(dp), intent(inout) :: C(:,:)

    !**** LOCAL ***********************************!

    !logical :: changeFmtB
    integer :: trA, trB, ot
    integer :: iprow, ipcol, nprow, npcol, mpi_err
    real(dp), allocatable :: tmp(:,:)
    integer :: dims_before(2), dims_after(2)
    integer :: desc_before(9), desc_after(9)

    !**** GLOBAL **********************************!
#ifdef HAVE_MPI
    character(1) :: psp_proc_order

    integer :: psp_mpi_comm_world
    integer :: psp_mpi_size
    integer :: psp_nprow
    integer :: psp_npcol
    integer :: psp_bs_def_row
    integer :: psp_bs_def_col
    integer :: psp_update_rank ! In SUMMA for C=A*B, C is computed by rank psp_update_rank=1 local update
    integer :: psp_bs_num
    integer :: psp_icontxt ! BLACS context handle used by psp
    integer :: psp_mpi_comm_cart, psp_mpi_comm_row, psp_mpi_comm_col


    common /psp_grid2D/ psp_mpi_comm_world, psp_mpi_size, psp_nprow, psp_npcol
    common /psp_grid2D/ psp_bs_def_row, psp_bs_def_col
    common /psp_grid2D/ psp_update_rank ! In SUMMA for C=A*B, C is computed by rank psp_update_rank=1 local update
    common /psp_grid2D/ psp_bs_num
    common /psp_grid2D/ psp_icontxt ! BLACS context handle used by psp
    common /psp_grid2D/ psp_mpi_comm_cart, psp_mpi_comm_row, psp_mpi_comm_col
#endif

    !**********************************************!

    integer, external :: numroc ! it is a function to compute local size

    !**********************************************!
    if (alpha/=0.0_dp) then
       !if (B%str_type=='coo') then
       !   changeFmtB=.true.
       !   call psp_coo2csc(B)
       !else
       !   changeFmtB=.false.
       !end if
       call psp_process_opM(opA,trA)
       call psp_process_opM(opB,trB)
       ! operation table
       if (trA==0 .and. trB==0) then
          ot=1
       else if (trA==0 .and. trB>=1) then
          ot=2
       else if (trA>=1 .and. trB==0) then
          ot=3
       else if (trA>=1 .and. trB>=1) then
          ot=4
       else
          call die('mm_dmultiply: invalid implementation')
       end if

       select case (ot)
       case (1)
          call psp_gemspm_nn(M,N,K,A,opA,B,opB,C,alpha,beta)
       case (2)
          call psp_gemspm_nt(M,N,K,A,opA,B,opB,C,alpha,beta)
       case (3)
          ! TODO: optimize the following code by optimizing the communication according to the 
          ! matrix format and dimension 
          call blacs_gridinfo(psp_icontxt,nprow,npcol,iprow,ipcol)
          if (.true.) then
             ! tmp = transpose(A)
             dims_before(1)=numroc(K,psp_bs_def_row,iprow,0,nprow)
             dims_before(2)=numroc(M,psp_bs_def_col,ipcol,0,npcol)
             call descinit(desc_before,K,M,psp_bs_def_row,psp_bs_def_col,0,0,psp_icontxt,dims_before(1),mpi_err)
             dims_after(1)=numroc(M,psp_bs_def_row,iprow,0,nprow)
             dims_after(2)=numroc(K,psp_bs_def_col,ipcol,0,npcol)
             call descinit(desc_after,M,K,psp_bs_def_row,psp_bs_def_col,0,0,psp_icontxt,dims_after(1),mpi_err)
             allocate(tmp(dims_after(1),dims_after(2)))
             call pdtran(M,K,1.0_dp,A,1,1,desc_before, 0.0_dp,tmp,1,1,desc_after)

             ! compute C = alpha*tmp*B + beta*C
             call psp_gemspm_nn(M,N,K,tmp,'n',B,opB,C,alpha,beta)
          end if
          !call psp_gemspm_tn(M,N,K,A,opA,B,opB,C,alpha,beta)
       case (4)
          call psp_gemspm_tt(M,N,K,A,opA,B,opB,C,alpha,beta)
       end select
       ! change format back to plan
       !if (changeFmtB) then
       !   call psp_csc2coo(B)
       !end if
    else
       if (beta/=0.0_dp) C=beta*C
    end if

  end subroutine psp_dgemspm


  subroutine psp_zgemspm(M,N,K,A,opA,B,opB,C,alpha,beta)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) ::            M, N, K ! Globally, A of size M by K, B of size K by N
    character(1), intent(in) :: opA ! form of op(A): 'n/N' for A, 't/T/c/C' for A^T
    character(1), intent(in) :: opB ! form of op(B)
    complex(dp), intent(in)  ::   alpha, beta ! scalar
    complex(dp), intent(in) :: A(:,:)

    !**** INOUT ***********************************!

    type(psp_matrix_spm), intent(in) ::   B ! matrix B
    complex(dp), intent(inout) :: C(:,:)

    !**** LOCAL ***********************************!
    !logical :: changeFmtB
    integer :: trA, trB, ot
    integer :: iprow, ipcol, nprow, npcol, mpi_err
    complex(dp), allocatable :: tmp(:,:)
    integer :: dims_before(2), dims_after(2)
    integer :: desc_before(9), desc_after(9)

    !**** GLOBAL **********************************!
#ifdef HAVE_MPI
    character(1) :: psp_proc_order

    integer :: psp_mpi_comm_world
    integer :: psp_mpi_size
    integer :: psp_nprow
    integer :: psp_npcol
    integer :: psp_bs_def_row
    integer :: psp_bs_def_col
    integer :: psp_update_rank ! In SUMMA for C=A*B, C is computed by rank psp_update_rank=1 local update
    integer :: psp_bs_num
    integer :: psp_icontxt ! BLACS context handle used by psp
    integer :: psp_mpi_comm_cart, psp_mpi_comm_row, psp_mpi_comm_col


    common /psp_grid2D/ psp_mpi_comm_world, psp_mpi_size, psp_nprow, psp_npcol
    common /psp_grid2D/ psp_bs_def_row, psp_bs_def_col
    common /psp_grid2D/ psp_update_rank ! In SUMMA for C=A*B, C is computed by rank psp_update_rank=1 local update
    common /psp_grid2D/ psp_bs_num
    common /psp_grid2D/ psp_icontxt ! BLACS context handle used by psp
    common /psp_grid2D/ psp_mpi_comm_cart, psp_mpi_comm_row, psp_mpi_comm_col
#endif

    !**********************************************!

    integer, external :: numroc ! it is a function to compute local size

    !**********************************************!
    if (alpha/=cmplx_0) then
       !if (B%str_type=='coo') then
       !   changeFmtB=.true.
       !   call psp_coo2csc(B)
       !else
       !   changeFmtB=.false.
       !end if
       call psp_process_opM(opA,trA)
       call psp_process_opM(opB,trB)
       ! operation table
       if (trA==0 .and. trB==0) then
          ot=1
       else if (trA==0 .and. trB>=1) then
          ot=2
       else if (trA>=1 .and. trB==0) then
          ot=3
       else if (trA>=1 .and. trB>=1) then
          ot=4
       else
          call die('mm_dmultiply: invalid implementation')
       end if

       select case (ot)
       case (1)
          call psp_gemspm_nn(M,N,K,A,opA,B,opB,C,alpha,beta)
       case (2)
          call psp_gemspm_nt(M,N,K,A,opA,B,opB,C,alpha,beta)
       case (3)
          ! TODO: optimize the following code by optimizing the communication according to the 
          ! matrix format and dimension 
          call blacs_gridinfo(psp_icontxt,nprow,npcol,iprow,ipcol)
          if (.true.) then
             ! tmp = transpose(A)
             dims_before(1)=numroc(K,psp_bs_def_row,iprow,0,nprow)
             dims_before(2)=numroc(M,psp_bs_def_col,ipcol,0,npcol)
             call descinit(desc_before,K,M,psp_bs_def_row,psp_bs_def_col,0,0,psp_icontxt,dims_before(1),mpi_err)
             dims_after(1)=numroc(M,psp_bs_def_row,iprow,0,nprow)
             dims_after(2)=numroc(K,psp_bs_def_col,ipcol,0,npcol)
             call descinit(desc_after,M,K,psp_bs_def_row,psp_bs_def_col,0,0,psp_icontxt,dims_after(1),mpi_err)
             allocate(tmp(dims_after(1),dims_after(2)))
             if (trA==1) then
                call pztranc(M,K,cmplx_1,A,1,1,desc_before,cmplx_0,tmp,1,1,desc_after)
             else
                call pztranu(M,K,cmplx_1,A,1,1,desc_before,cmplx_0,tmp,1,1,desc_after)
             end if
 
             ! compute C = alpha*tmp*B + beta*C
             call psp_gemspm_nn(M,N,K,tmp,'n',B,opB,C,alpha,beta)
          end if
          !call psp_gemspm_tn(M,N,K,A,opA,B,opB,C,alpha,beta)
       case (4)
          call psp_gemspm_tt(M,N,K,A,opA,B,opB,C,alpha,beta)
       end select
       ! change format back to plan
       !if (changeFmtB) then
       !   call psp_csc2coo(B)
       !end if
    else
       if (beta/=cmplx_0) C=beta*C
    end if
  end subroutine psp_zgemspm

  subroutine die(message)
    implicit none

    !**** INPUT ***********************************!

    character(*), intent(in), optional :: message

    !**** INTERNAL ********************************!

    logical, save :: log_start=.false.

    integer :: log_unit

    !**********************************************!

    if (log_start) then
       open(newunit=log_unit,file='MatrixSwitch.log',position='append')
    else
       open(newunit=log_unit,file='MatrixSwitch.log',status='replace')
       log_start=.true.
    end if
    write(log_unit,'(a)'), 'FATAL ERROR in matrix_switch!'
    write(log_unit,'(a)'), message
    close(log_unit)
    stop

  end subroutine die


END MODULE pspMspm
