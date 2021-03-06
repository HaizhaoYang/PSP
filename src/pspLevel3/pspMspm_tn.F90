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

MODULE pspMspm_tn
  use pspVariable
  use pspUtility
  use pspMPI
  use pspLevel1
  use pspLevel2

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

  interface psp_gemspm_tn
     module procedure psp_dgemspm_tn
     module procedure psp_zgemspm_tn
  end interface psp_gemspm_tn

  interface die
     module procedure die
  end interface die

  public :: psp_gemspm_tn

contains

  !================================================!
  !    sparse pdgemm: mspm  C=alpha*A'*B+beta*C    !
  !================================================!

  subroutine psp_dgemspm_tn(M,N,K,A,opA,B,opB,C,alpha,beta)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) ::            M, N, K ! Globally, A of size M by K, B of size K by N
    character(1), intent(in) :: opA ! form of op(A): 'n/N' for A, 't/T/c/C' for A^T
    character(1), intent(in) :: opB ! form of op(B)
    real(dp), intent(in)  ::   alpha, beta ! scalar
    real(dp), intent(in) :: A(:,:)
    type(psp_matrix_spm), intent(in) ::   B ! matrix B

    !**** INOUT ***********************************!

    real(dp), intent(inout) :: C(:,:)

    !**** LOCAL ***********************************!

    integer :: kloop, idx_k_row, idx_k_col, idx_pcol, idx_prow, cnt1, cnt2, chunk_sz
    integer :: trA, trB
    integer :: iprow, ipcol, nprow, npcol, mpi_err
    real(dp), allocatable :: A_loc(:,:), C_loc(:,:), CC_loc(:,:)
    real(dp), allocatable :: B_loc_val(:)
    integer, allocatable :: B_loc_idx1(:), B_loc_idx2(:)
    integer :: A_loc_dim(2), B_loc_dim(2), C_loc_dim(2), coords(2)
    integer :: width, glb_st, glb_ed, loc_st, loc_ed, nnz_loc, L

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

    ! obtain grid information, working on the (iprow,ipcol) processor in a grid of size nprow by npcol
    call blacs_gridinfo(psp_icontxt,nprow,npcol,iprow,ipcol)
    ! set up local variables
    ! allocate local matrices
    A_loc_dim(1)=numroc(K,psp_bs_def_row,iprow,0,nprow)
    A_loc_dim(2)=psp_update_rank
    allocate(A_loc(A_loc_dim(1),A_loc_dim(2)))
    A_loc=0.0_dp
    B_loc_dim(1)=numroc(K,psp_bs_def_row,iprow,0,nprow)
    B_loc_dim(2)=numroc(N,psp_bs_def_col,ipcol,0,npcol)
    C_loc_dim(1)=psp_update_rank
    C_loc_dim(2)=numroc(N,psp_bs_def_col,ipcol,0,npcol)
    allocate(C_loc(C_loc_dim(1),C_loc_dim(2)))
    C_loc=0.0_dp
    allocate(CC_loc(C_loc_dim(1),C_loc_dim(2)))
    CC_loc=0.0_dp

    !**** Start the operation ********************!

    ! begin SUMMA
    do kloop = 1, M, psp_update_rank
       idx_k_row = (kloop-1)/psp_bs_def_row+1
       idx_k_col = (kloop-1)/psp_bs_def_col+1
       glb_st=kloop
       glb_ed=min(kloop+psp_update_rank,M+1)-1
       width=glb_ed-glb_st+1
       idx_pcol = mod(idx_k_col-1,npcol)
       if (ipcol==idx_pcol) then
          call psp_idx_glb2loc(glb_st,psp_bs_def_col,npcol,loc_st)
          if (width<psp_update_rank) then
             A_loc=0.0_dp
          endif
          !call psp_copy_m('n',A_loc_dim(1),width,A,1,loc_st,A_loc,1,1,1.0_dp,0.0_dp)
          A_loc(:,1:width)=A(1:A_loc_dim(1),loc_st:loc_st+width-1)
       end if

       ! boardcast in row
       call MPI_Bcast(A_loc, A_loc_dim(1)*A_loc_dim(2), MPI_DOUBLE, idx_pcol, psp_mpi_comm_row,mpi_err)

       ! compute local update of C
       ! C_loc = (A_loc)^t*B
       call psp_sst_gemspm(width,C_loc_dim(2),B_loc_dim(1),opA,opB, &
            1.0_dp,A_loc,1,1,B%row_ind,B%col_ptr,B%dval,C_loc,1,1,0.0_dp)

       idx_prow = mod(idx_k_row-1,nprow ) ! identify the processor owing C(kth block,:), the cart coordinate
       ! boardcast in column
       call MPI_REDUCE(C_loc, CC_loc, C_loc_dim(1)*C_loc_dim(2), MPI_DOUBLE, MPI_SUM, idx_prow, psp_mpi_comm_col, mpi_err)

       if (iprow==idx_prow) then
          call psp_idx_glb2loc(glb_st,psp_bs_def_row,nprow,loc_st)
          !C=beta*C+C_loc
          !call psp_copy_m('n',width,C_loc_dim(2),CC_loc,1,1,C,loc_st,1,alpha,beta)
          C(loc_st:loc_st+width-1,1:C_loc_dim(2))=beta*C(loc_st:loc_st+width-1,1:C_loc_dim(2))+alpha*CC_loc(1:width,1:C_loc_dim(2))
       end if
    enddo

    if (allocated(A_loc)) deallocate(A_loc)
    if (allocated(C_loc)) deallocate(C_loc)
    if (allocated(CC_loc)) deallocate(CC_loc)
    if (allocated(B_loc_val)) deallocate(B_loc_val)
    if (allocated(B_loc_idx1)) deallocate(B_loc_idx1)
    if (allocated(B_loc_idx2)) deallocate(B_loc_idx2)

  end subroutine psp_dgemspm_tn

  !================================================!
  !    sparse pzgemm: mspm  C=alpha*A'*B+beta*C    !
  !================================================!

  subroutine psp_zgemspm_tn(M,N,K,A,opA,B,opB,C,alpha,beta)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) ::            M, N, K ! Globally, A of size M by K, B of size K by N
    character(1), intent(in) :: opA ! form of op(A): 'n/N' for A, 't/T/c/C' for A^T
    character(1), intent(in) :: opB ! form of op(B)
    complex(dp), intent(in)  ::   alpha, beta ! scalar
    complex(dp), intent(in) :: A(:,:)
    type(psp_matrix_spm), intent(in) ::   B ! matrix B

    !**** INOUT ***********************************!

    complex(dp), intent(inout) :: C(:,:)

    !**** LOCAL ***********************************!

    integer :: kloop, idx_k_row, idx_k_col, idx_pcol, idx_prow, cnt1, cnt2, chunk_sz
    integer :: trA, trB
    integer :: iprow, ipcol, nprow, npcol, mpi_err
    complex(dp), allocatable :: A_loc(:,:), C_loc(:,:), CC_loc(:,:)
    complex(dp), allocatable :: B_loc_val(:)
    integer, allocatable :: B_loc_idx1(:), B_loc_idx2(:)
    integer :: A_loc_dim(2), B_loc_dim(2), C_loc_dim(2), coords(2)
    integer :: width, glb_st, glb_ed, loc_st, loc_ed, nnz_loc, L

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

    ! obtain grid information, working on the (iprow,ipcol) processor in a grid of size nprow by npcol
    call blacs_gridinfo(psp_icontxt,nprow,npcol,iprow,ipcol)
    ! set up local variables
    ! allocate local matrices
    A_loc_dim(1)=numroc(K,psp_bs_def_row,iprow,0,nprow)
    A_loc_dim(2)=psp_update_rank
    allocate(A_loc(A_loc_dim(1),A_loc_dim(2)))
    A_loc=cmplx_0
    B_loc_dim(1)=numroc(K,psp_bs_def_row,iprow,0,nprow)
    B_loc_dim(2)=numroc(N,psp_bs_def_col,ipcol,0,npcol)
    C_loc_dim(1)=psp_update_rank
    C_loc_dim(2)=numroc(N,psp_bs_def_col,ipcol,0,npcol)
    allocate(C_loc(C_loc_dim(1),C_loc_dim(2)))
    C_loc=cmplx_0
    allocate(CC_loc(C_loc_dim(1),C_loc_dim(2)))
    CC_loc=cmplx_0

    !**** Start the operation ********************!

    ! begin SUMMA
    do kloop = 1, M, psp_update_rank
       idx_k_row = (kloop-1)/psp_bs_def_row+1
       idx_k_col = (kloop-1)/psp_bs_def_col+1
       glb_st=kloop
       glb_ed=min(kloop+psp_update_rank,M+1)-1
       width=glb_ed-glb_st+1
       idx_pcol = mod(idx_k_col-1,npcol)
       if (ipcol==idx_pcol) then
          call psp_idx_glb2loc(glb_st,psp_bs_def_col,npcol,loc_st)
          if (width<psp_update_rank) then
             A_loc=cmplx_0
          endif
          !call psp_copy_m('n',A_loc_dim(1),width,A,1,loc_st,A_loc,1,1,cmplx_1,cmplx_0)
          A_loc(:,1:width)=A(1:A_loc_dim(1),loc_st:loc_st+width-1)
       end if

       ! boardcast in row
       call MPI_Bcast(A_loc, A_loc_dim(1)*A_loc_dim(2), MPI_DOUBLE_COMPLEX, idx_pcol, psp_mpi_comm_row,mpi_err)

       ! compute local update of C
       ! C_loc = (A_loc)^t*B
       call psp_sst_gemspm(width,C_loc_dim(2),B_loc_dim(1),opA,opB, &
            cmplx_1,A_loc,1,1,B%row_ind,B%col_ptr,B%zval,C_loc,1,1,cmplx_0)

       idx_prow = mod(idx_k_row-1,nprow ) ! identify the processor owing C(kth block,:), the cart coordinate
       ! boardcast in column
       call MPI_REDUCE(C_loc, CC_loc, C_loc_dim(1)*C_loc_dim(2), MPI_DOUBLE_COMPLEX, MPI_SUM, idx_prow, psp_mpi_comm_col, mpi_err)

       if (iprow==idx_prow) then
          call psp_idx_glb2loc(glb_st,psp_bs_def_row,nprow,loc_st)
          !C=beta*C+C_loc
          !call psp_copy_m('n',width,C_loc_dim(2),CC_loc,1,1,C,loc_st,1,alpha,beta)
          C(loc_st:loc_st+width-1,1:C_loc_dim(2))=beta*C(loc_st:loc_st+width-1,1:C_loc_dim(2))+alpha*CC_loc(1:width,1:C_loc_dim(2))
       end if
    enddo

    if (allocated(A_loc)) deallocate(A_loc)
    if (allocated(C_loc)) deallocate(C_loc)
    if (allocated(CC_loc)) deallocate(CC_loc)
    if (allocated(B_loc_val)) deallocate(B_loc_val)
    if (allocated(B_loc_idx1)) deallocate(B_loc_idx1)
    if (allocated(B_loc_idx2)) deallocate(B_loc_idx2)

  end subroutine psp_zgemspm_tn

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


END MODULE pspMspm_tn
