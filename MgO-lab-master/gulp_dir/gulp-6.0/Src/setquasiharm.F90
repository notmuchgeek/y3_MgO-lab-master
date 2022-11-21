  subroutine setquasiharm(dsqh,mint,mintloc,mcv4)
!
!  Calculates the matrix which maps the atomic free energy
!  derivatives on to strain for the quasiharmonic approxn.
!
!   9/97 Created
!  12/00 Generalised to any number of dimensions
!   2/03 matinv replaced to accelerate speed
!   9/16 Matrix inversion moved to subroutine call
!   2/17 Blocksize added to call to matrix_inversion_library
!   5/17 Corrected to use modified matinv routines
!   5/17 mintloc and mcv4 added to argument list
!   5/17 Parallelisation added
!   8/17 INTEL option added to workaround bugs in the MKL library
!
!  On entry :
!
!    derv2 = internal-internal second derivative matrix
!    dsqh  = external-internal second derivative matrix
!
!  On exit :
!
!    derv2 = inverse of internal-internal second derivatives
!    dsqh  = product of inverse 2nd derivs and external-
!            internal matrix
!
!  Conditions of use:
!
!  GULP is available free of charge to academic institutions
!  and non-commerical establishments only. Copies should be
!  obtained from the author only and should not be distributed
!  in any form by the user to a third party without the express
!  permission of the author. This notice applies to all parts
!  of the program, except any library routines which are
!  distributed with the code for completeness. All rights for
!  such routines remain with the original distributor.
!
!  No claim is made that this program is free from errors and
!  no liability will be accepted for any loss or damage that
!  may result. The user is responsible for checking the validity
!  of their results.
!
!  Copyright Curtin University 2017
!
!  Julian Gale, CIC, Curtin University, August 2017
!
  use current,    only : nstrains
  use derivatives
  use parallel,   only : nblocksize, nprocs
  implicit none
!
!  Passed variables
!
  integer(i4),                   intent(in)    :: mint
  integer(i4),                   intent(in)    :: mintloc
  integer(i4),                   intent(in)    :: mcv4
  real(dp),                      intent(out)   :: dsqh(maxd2,*)
!
!  Local variables
!
#ifdef INTEL
  complex(dpc),  dimension(:,:), allocatable   :: cmat
#endif
  integer(i4)                                  :: i
  integer(i4)                                  :: ifail
  integer(i4)                                  :: j
  integer(i4)                                  :: k
  integer(i4)                                  :: status
  real(dp),      dimension(:,:), allocatable   :: dtmp1
  real(dp),      dimension(:,:), allocatable   :: dtmp2
!************************************
!  Invert second derivative matrix  *
!************************************
  ifail = 0
  if (nprocs.gt.1) then
#ifdef INTEL
!
!  Work around for Intel MKL issues
!
    allocate(cmat(maxd2,mintloc),stat=status)
    if (status/=0) call outofmemory('setquasiharm','cmat')
!
    do i = 1,mintloc
      do j = 1,mint
        cmat(j,i) = dcmplx(derv2(j,i),0.0_dp)
      enddo
    enddo
!
    call cmatrix_inversion_library(mint-3_i4,4_i4,maxd2,cmat,0_i4,ifail)
!
    do i = 1,mintloc
      do j = 1,mint
        derv2(j,i) = dble(cmat(j,i))
      enddo
    enddo
!
    deallocate(cmat,stat=status)
    if (status/=0) call deallocate_error('setquasiharm','cmat')
#else
    call matrix_inversion_library(mint-3_i4,4_i4,maxd2,3_i4*nblocksize,derv2,0_i4,ifail)
#endif
  else
    call matrix_inversion_library(mint-3_i4,4_i4,maxd2,3_i4*nblocksize,derv2,0_i4,ifail)
  endif
  if (ifail.ne.0) then
    call outerror('matrix inversion failed in setquasiharm',0_i4)
    call stopnow('setquasiharm')
  endif
!************************************************************************
!  Multiply inverse second derivatives by mixed strain internal matrix  *
!************************************************************************
  if (nprocs.gt.1) then
    dsqh(1:mint,1:nstrains) = 0.0_dp
    if (mcv4.gt.0) then
!
!  If mcv4 is zero then there is no contribution from this node
!
      do i = 1,nstrains
        do j = 4,mint
          do k = mcv4,mintloc
            dsqh(j,i) = dsqh(j,i) + derv2(j,k)*derv3(k,i)
          enddo
        enddo
      enddo
    endif
!
!  Global sum of dsqh
!
    allocate(dtmp1(mint,nstrains),stat=status)
    if (status/=0) call outofmemory('setquasiharm','dtmp1')
    allocate(dtmp2(mint,nstrains),stat=status)
    if (status/=0) call outofmemory('setquasiharm','dtmp2')
!
    do i = 1,nstrains
      do j = 1,mint
        dtmp1(j,i) = dsqh(j,i)
      enddo
    enddo
!
    call sumall(dtmp1,dtmp2,mint*nstrains,"setquasiharm","dsqh")
!
    do i = 1,nstrains
      do j = 1,mint
        dsqh(j,i) = dtmp2(j,i)
      enddo
    enddo
!
    deallocate(dtmp2,stat=status)
    if (status/=0) call deallocate_error('setquasiharm','dtmp2')
    deallocate(dtmp1,stat=status)
    if (status/=0) call deallocate_error('setquasiharm','dtmp1')
  else
    do i = 1,nstrains
      dsqh(1,i) = 0.0_dp
      dsqh(2,i) = 0.0_dp
      dsqh(3,i) = 0.0_dp
      do j = 4,mint
        dsqh(j,i) = 0.0_dp
        do k = 4,mint
          dsqh(j,i) = dsqh(j,i) + derv2(k,j)*derv3(k,i)
        enddo
      enddo
    enddo
  endif
!
  return
  end
