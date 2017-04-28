!************************************************************************!
!   Copyright (c) 2015-2017, Haizhao Yang                                !
!   All rights reserved.                                                 !
!                                                                        !
!   This file is part of PSP and is under the BSD 2-Clause License,! 
!   which can be found in the LICENSE file in the root directory, or at  !
!   http://opensource.org/licenses/BSD-2-Clause                          !
!************************************************************************!

! This code test the pzgemm in scalapack

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

program pzgemmScaling
  use pspBLAS


  implicit none
#ifdef HAVE_MPI
  include 'mpif.h'
#endif

  !**** PARAMS **********************************!

  integer, parameter :: dp=selected_real_kind(15,300)
  complex(dp), parameter :: cmplx_1=(1.0_dp,0.0_dp)
  complex(dp), parameter :: cmplx_i=(0.0_dp,1.0_dp)
  complex(dp), parameter :: cmplx_0=(0.0_dp,0.0_dp)

  !**** VARIABLES *******************************!

  character(1) :: order
  character(21) :: file_name

  logical :: dealloc

  integer :: mpi_err, mpi_size, mpi_rank, m, n, k, iostat, counti, countf, count_rate
  integer :: nprow, npcol, bs_def_row, bs_def_col, icontxt, iprow, ipcol, niter, i, idxr, idxc
  integer :: H_dim(2), S_dim(2), D_dim(2), info, Ht_dim(2)
  integer :: desc_H(9), desc_S(9), desc_D(9), desc_Ht(9)

  complex(dp) :: he, se
  real(dp) :: her,ser, alpha, beta
  real(dp) :: rpt, cpt
  real(dp), allocatable :: H(:,:), S(:,:), D(:,:), Dtrue(:,:), Ht(:,:)
  complex(dp), allocatable :: zH(:,:), zS(:,:), zD(:,:), zDtrue(:,:)
  real(dp) :: dtime, t0, t1, mm_err, err
  integer*4 timeArray(3)    ! Holds the hour, minute, and second
  complex(dp) :: zalpha, zbeta

  type(psp_matrix_spm) :: Hsp, Ssp, Dsp, Dtruesp, Htsp ! sparse matrices in pspBLAS
  type(psp_matrix_spm) :: zHsp, zSsp, zDsp, zDtruesp ! sparse matrices in pspBLAS
  character(3) :: fmtH, fmtS, fmtD ! storage type of the sparse matrix, 'coo' or 'csc'
  real(dp) :: thre ! threshold parameter for converting a dense matrix to a sparse amtrix
  integer, allocatable :: idx1(:), idx2(:)! vectors for sparse matrices
  real(dp), allocatable :: val(:)
  complex(dp), allocatable :: zval(:)
  integer :: desc_Hsp(9), desc_Ssp(9), desc_Dsp(9)
  character(1) :: matdescr(6)
  integer, allocatable :: idx3(:)
  !**********************************************!

  integer, external :: numroc ! it is a function to compute local size

  !**********************************************!
  ! initialization first

  ! initialize information in MPI
  call mpi_init(mpi_err)
  call mpi_comm_size(mpi_comm_world,mpi_size,mpi_err)
  call mpi_comm_rank(mpi_comm_world,mpi_rank,mpi_err)

  ! set up parameters for parallel computing test
  niter=20
  nprow=INT(sqrt(DBLE(mpi_size)))
  npcol=mpi_size/nprow
  order='r' ! important TODO: check how to adapt to different orders
  bs_def_row=10
  bs_def_col=10

  ! initialize information in scalapack
  call blacs_get(-1,0,icontxt)
  call blacs_gridinit(icontxt,order,nprow,npcol)
  ! obtain grid information, working on the (k,l) processor in a grid i by j
  call blacs_gridinfo(icontxt,nprow,npcol,iprow,ipcol)

  ! initialized grid information in pspBLAS
  call psp_gridinit_2D(mpi_comm_world,mpi_size,nprow,order,bs_def_row,bs_def_col,icontxt)

  ! set up parameters for parallel computing test
  niter=1
  matdescr(1)='G'
  matdescr(4)='F'

  !*************************************************************************!
  ! generate test matrices

  if (.true.) then
     ! random matrices
     m=1000 ! global matrix size
     n=700
     k=500

     H_dim(1)=numroc(m,bs_def_row,iprow,0,nprow) ! use 'numroc' to compute the size of a local matrix, row
     H_dim(2)=numroc(k,bs_def_col,ipcol,0,npcol) ! use 'numroc' to compute the size of a local matrix, column
     call descinit(desc_H,m,k,bs_def_row,bs_def_col,0,0,icontxt,H_dim(1),info) ! initialize the descriptor of the global matrix H
     allocate(H(H_dim(1),H_dim(2)))! allocate matrix H
     H=0.0_dp ! assign zero values

     Ht_dim(1)=numroc(k,bs_def_row,iprow,0,nprow) ! use 'numroc' to compute the size of a local matrix, row
     Ht_dim(2)=numroc(m,bs_def_col,ipcol,0,npcol) ! use 'numroc' to compute the size of a local matrix, column
     call descinit(desc_Ht,k,m,bs_def_row,bs_def_col,0,0,icontxt,Ht_dim(1),info) ! initialize the descriptor of the global matrix H
     allocate(Ht(Ht_dim(1),Ht_dim(2)))! allocate matrix H
     Ht=0.0_dp ! assign zero values

     ! initialize and allocate S, D, Dtrue similarly
     S_dim(1)=numroc(k,bs_def_row,iprow,0,nprow)
     S_dim(2)=numroc(n,bs_def_col,ipcol,0,npcol)
     call descinit(desc_S,k,n,bs_def_row,bs_def_col,0,0,icontxt,S_dim(1),info) !
     allocate(S(S_dim(1),S_dim(2)))
     S=0.0_dp


     D_dim(1)=numroc(m,bs_def_row,iprow,0,nprow)
     D_dim(2)=numroc(n,bs_def_col,ipcol,0,npcol)
     call descinit(desc_D,m,n,bs_def_row,bs_def_col,0,0,icontxt,D_dim(1),info) !
     allocate(D(D_dim(1),D_dim(2)))
     D=0.0_dp

     allocate(Dtrue(D_dim(1)+1,D_dim(2)))
     Dtrue=0.0_dp

     ! generate random matrices
     if (.true.) then
        call init_random_seed()
        call RANDOM_NUMBER(H)
        call RANDOM_NUMBER(Ht)
        call RANDOM_NUMBER(S)
     end if
  end if

  !***********************************
  print *, 'test sparse summation in csc format'
  fmtH='csc' ! specify storage format, 'coo' or 'csc'
  fmtS='csc'
  fmtD='csc'
  alpha=1.5_dp
  beta=1.2_dp

  ! first method to generate a sparse matrix: thresholding a dense matrix in MatrixSwitch
  thre = 0.5_dp
  call psp_den2sp_m(H,desc_H,Hsp,fmtH,thre)   ! need to initialize all matrices
  call psp_den2sp_m(Ht,desc_Ht,Htsp,fmtH,thre)   ! need to initialize all matrices
  !  call psp_den2sp_m(S,desc_S,Ssp,fmtS,thre)   ! need to initialize all matrices
  !  call psp_spm_zeros(Dsp,m,n,fmtD,.true.)     ! need to initialize all matrices

  do idxc=1,Ht_dim(2)
     do idxr=1,Ht_dim(1)
        if (abs(Ht(idxr,idxc))<thre) Ht(idxr,idxc) = 0.0_dp
     end do
  end do
  do idxc=1,H_dim(2)
     do idxr=1,H_dim(1)
        if (abs(H(idxr,idxc))<thre) H(idxr,idxc) = 0.0_dp
     end do
  end do
  do idxc=1,S_dim(2)
     do idxr=1,S_dim(1)
        if (abs(S(idxr,idxc))<thre) S(idxr,idxc) = 0.0_dp
     end do
  end do

  print *, 'begin psp'
  t0 = MPI_Wtime()
  do i=1,niter
     call psp_sst_gespmm(m,n,k, 'n','n', &
          1.0_dp,Hsp%row_ind,Hsp%col_ptr,Hsp%dval,S,1,1,D,1,1,0.0_dp)
  enddo
  t1 = MPI_Wtime()
  dtime=(t1-t0)/DBLE(niter)
  if (mpi_rank==0) print *, 'time of PSP = n n ', dtime


  print *, 'set up MKL'

  if (allocated(idx3)) deallocate(idx3)
  allocate(idx3(m))
  idx3(1:m)=Hsp%col_ptr(2:m+1)
  print *, 'begin MKL'
  t0 = MPI_Wtime()
  do i=1,niter
     call mkl_dcscmm('N',m,n,k,1.0_dp,matdescr,Hsp%dval,Hsp%row_ind,&
          Hsp%col_ptr,idx3,S,k,0.0_dp,Dtrue,m+1)
  enddo
  t1 = MPI_Wtime()
  dtime=(t1-t0)/DBLE(niter)
  if (mpi_rank==0) print *, 'time of MKL = n n ', dtime

  print *, 'check results'
  err=MAXVAL(abs(Dtrue(1:m,1:n)-D))
  print *, 'real case              ', err

  print *, 'begin psp'
  t0 = MPI_Wtime()
  do i=1,niter
     call psp_sst_gespmm(m,n,k, 't','n', &
          1.0_dp,Htsp%row_ind,Htsp%col_ptr,Htsp%dval,S,1,1,D,1,1,0.0_dp)
  enddo
  t1 = MPI_Wtime()
  dtime=(t1-t0)/DBLE(niter)
  if (mpi_rank==0) print *, 'time of PSP = t n ', dtime


  print *, 'set up MKL'

  if (allocated(idx3)) deallocate(idx3)
  allocate(idx3(m))
  idx3(1:m)=Htsp%col_ptr(2:m+1)
  print *, 'begin MKL'
  t0 = MPI_Wtime()
  do i=1,niter
     call mkl_dcscmm('T',k,n,m,1.0_dp,matdescr,Htsp%dval,Htsp%row_ind,&
          Htsp%col_ptr,idx3,S,k,0.0_dp,Dtrue,m+1)
  enddo
  t1 = MPI_Wtime()
  dtime=(t1-t0)/DBLE(niter)
  if (mpi_rank==0) print *, 'time of MKL = t n ', dtime

  print *, 'check results'
  err=MAXVAL(abs(Dtrue(1:m,1:n)-D))
  print *, 'real case              ', err



  deallocate(H)
  deallocate(S)
  deallocate(D)
  deallocate(Dtrue)

  call psp_deallocate_spm(Hsp)
  call psp_deallocate_spm(Ssp)
  call psp_deallocate_spm(Dsp)

end program pzgemmScaling
