  subroutine thirdorderfc3
!  subroutine thirdorderfc3(iocptr)
!
!  Analytical calculation of the third order force constants
!
!  NOTE : Needs finishing for BSM 
!
!   4/15 Created from realrecip3d3
!   7/15 Dscg removed as it is not used
!   7/15 Missing terms added to match finite differences
!   8/15 rvghost vectors modified to yield the shortest 
!        equivalent distance
!   8/15 Thresholds added for third order force constants
!   8/15 Fourbody call moved outside i-j loops
!   9/15 Sixbody terms added
!  10/15 Threebody terms moved out of the twobody loop
!  10/15 Real & reciprocal terms moved to a subroutine
!   6/17 Module files renamed to gulp_files
!   6/17 No file written if nfc3 is zero
!   2/18 Trace added
!
!  Assumptions:
!
!  - cores are sorted to come before shells
!  - there are no partial occupancies
!  - this is a single region calculation
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
!  Copyright Curtin University 2018
!
!  Julian Gale, CIC, Curtin University, Feb 2018
!
  use configurations, only : nsuperghost
  use g_constants
  use control
  use current
  use datatypes
  use derivatives
  use feworkspace
  use gulp_files,     only : iout_shengBTE, lshengBTE
  use four
  use general
  use iochannels
  use ksample
  use kspace
  use m_three
  use molecule
  use parallel
  use realvectors
  use shells
  use six
  use symmetry
  use times
  use thresholds
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  use two
  implicit none
!
!  Passed variables
!
!  integer(i4)                                      :: iocptr(*)
!
!  Local variables
!
  integer(i4)                                      :: i
  integer(i4)                                      :: ig
  integer(i4)                                      :: ii
  integer(i4)                                      :: iim
  integer(i4)                                      :: indi
  integer(i4)                                      :: indj
  integer(i4)                                      :: indm
  integer(i4)                                      :: iout2
  integer(i4)                                      :: j
  integer(i4)                                      :: jg
  integer(i4)                                      :: jj
  integer(i4)                                      :: jjm
  integer(i4)                                      :: k
  integer(i4)                                      :: kk
  integer(i4)                                      :: kkm
  integer(i4)                                      :: l
  integer(i4)                                      :: m
  integer(i4)                                      :: maxlim        ! Number of variables including shells
  integer(i4)                                      :: maxlimc       ! Number of variables for cores only
  integer(i4)                                      :: maxlims       ! Number of variables for shells only
  integer(i4)                                      :: maxlimg       ! Number of variables including shells in ghost cell
  integer(i4)                                      :: maxlimcg      ! Number of variables for cores only in ghost cell
  integer(i4)                                      :: mint
  integer(i4)                                      :: mm
  integer(i4)                                      :: nfc3
  integer(i4)                                      :: nghostatom    ! Number of atoms in the ghost cell
  integer(i4)                                      :: nghostcore    ! Number of cores in the ghost cell
  integer(i4)                                      :: nghostcell    ! Number of ghost cells in the full cell
  integer(i4)                                      :: nnonzero
  integer(i4)                                      :: status
  real(dp)                                         :: g_cpu_time
  real(dp),    dimension(:,:),   allocatable, save :: Dsc
  real(dp)                                         :: fcmin
  real(dp),    dimension(:,:,:), allocatable, save :: fc3
  real(dp),    dimension(:,:),   allocatable, save :: rvghost
  real(dp)                                         :: rvtmp(3,3)
  real(dp)                                         :: time1
  real(dp)                                         :: time2
#ifdef TRACE
  call trace_in('thirdorderfc3')
#endif
!
  time1 = g_cpu_time()
!
!  Local variables
!
  fcmin = 1.0d-5
!
!  Find number of ghost cells
!
  nghostcell = nsuperghost(1,ncf)*nsuperghost(2,ncf)*nsuperghost(3,ncf)
!
!  Check the number of atoms is divisable by the number of ghost cells
!
  if (mod(numat,nghostcell).ne.0) then
    call outerror('number of atoms is not compatible with ghost_supercell',0_i4)
    call stopnow('thirdordern')
  endif
!
  nghostatom = numat/nghostcell
  nghostcore = ncore/nghostcell
!
!  Compute second order force constants (uncorrected for shells)
!
  call secondorderfc(.false.)
!********************
!  Allocate memory  *
!********************
  maxlim = 3*numat
  if (nbsmat.gt.0) maxlim = maxlim + numat
!
  maxlimg = 3*nghostatom
  if (nbsmat.gt.0) maxlimg = maxlimg + nghostatom
!
  maxlimc = 3*ncore
  maxlims = 3*nshell
  maxlimcg = 3*nghostcore
!
!  Ensure main arrays are large enough for all species
!
  if (maxlim.gt.maxd2u) then
    maxd2u = maxlim
    call changemaxd2
  endif
  if (maxlim.gt.maxd2) then
    maxd2 = maxlim
    call changemaxd2
  endif
!
!  Dimension third order force constants to hold only the core force constants
!
  allocate(fc3(maxlim,maxlim,maxlimcg),stat=status)
  if (status/=0) call outofmemory('thirdorderfc3','fc3')
!
  nfc3 = 0
  do i = 1,maxlimcg
    do j = 1,maxlim
      do k = 1,maxlim
        fc3(k,j,i) = 0.0_dp
      enddo
    enddo
  enddo
!
  mint = 3*numat
  maxlim = mint
  if (nbsmat.gt.0) maxlim = maxlim + numat
!*******************************************
!  Real and reciprocal space contribution  *
!*******************************************
  call realrecipfc3(maxlim,maxlimcg,fc3)
!****************************
!  Three-body contribution  *
!****************************
  if (nthb.gt.0) then
    call threefc3(maxlim,maxlimcg,fc3)
  endif
!***************************
!  Four-body contribution  *
!***************************
  if (nfor.gt.0) then
    call fourfc3(maxlim,maxlimcg,fc3)
  endif
!**************************
!  Six-body contribution  *
!**************************
  if (nsix.gt.0) then
    call sixfc3(maxlim,maxlimcg,fc3)
  endif
!***************************
!  Shell model correction  *
!***************************
  if (nshell.gt.0) then
!
!  Create array to store shell transformation matrix and workspace
!
    allocate(Dsc(maxlims,maxlimc),stat=status)
    if (status/=0) call outofmemory('thirdorderfc3','Dsc')
!
!  Compute shell transformation matrix
!
    call getshellmatrices(derv2,maxd2,Dsc,maxlims)
!
!  Loop over ghost atom core degrees of freedom and correct core-core-core third derivatives
!
    do m = 1,3*nghostcore
!
!  Term 1 : Dccc' = Dccc - (Dcs*Dss^-1)*Dscc - Dcsc*(Dss^-1*Dsc)
!
      do i = 1,maxlimc
        do j = 1,maxlimc
          do k = 1,maxlims
            fc3(j,i,m) = fc3(j,i,m) - Dsc(k,j)*fc3(maxlimc+k,i,m)
            fc3(j,i,m) = fc3(j,i,m) - Dsc(k,i)*fc3(j,maxlimc+k,m)
          enddo
        enddo
      enddo
!
!  Term 2 : Dccc' = Dccc - (Dcs*Dss^-1)*Dssc*(Dss^-1*Dsc)
!
      do i = 1,maxlimc
        do j = 1,maxlimc
          do k = 1,maxlims
            do l = 1,maxlims
              fc3(j,i,m) = fc3(j,i,m) - Dsc(l,j)*fc3(maxlimc+l,maxlimc+k,m)*Dsc(k,i)
            enddo
          enddo
        enddo
      enddo
    enddo
!
!  Free array that stored the shell transformation matrix and workspace
!
    deallocate(Dsc,stat=status)
    if (status/=0) call deallocate_error('thirdorderfc3','Dsc')
  endif
!
!  Threshold force constants
!
  indm = -3
  do m = 1,nghostcore
    indm = indm + 3
    indi = -3
    do i = 1,ncore
      indi = indi + 3
      indj = -3
      do j = 1,ncore
        indj = indj + 3
        nnonzero = 0
        do mm = 1,3
          do ii = 1,3
            do jj = 1,3
              if (abs(fc3(indj+jj,indi+ii,indm+mm)).lt.thresh_fc3_ind) fc3(indj+jj,indi+ii,indm+mm) = 0.0_dp
              if (abs(fc3(indj+jj,indi+ii,indm+mm)).gt.thresh_fc3_tot) nnonzero = nnonzero + 1
            enddo
          enddo
        enddo
        if (nnonzero.eq.0) then
          do mm = 1,3
            do ii = 1,3
              do jj = 1,3
                fc3(indj+jj,indi+ii,indm+mm) = 0.0_dp
              enddo
            enddo
          enddo
        endif
      enddo
    enddo
  enddo
!
!  Complete on-diagonal blocks by sum rule
!
  indm = -3
  do m = 1,nghostcore
    indm = indm + 3
    i = (m - 1)*nghostcell
    indi = 3*i
    i = i + 1
!
!  Zero block
!
    do mm = 1,3
      do ii = 1,3
        do jj = 1,3
          fc3(indi+jj,indi+ii,indm+mm) = 0.0_dp
        enddo
      enddo
    enddo
!
!  Set diagonal blocks to negative sum of off-diagonals
!
    indj = -3
    do j = 1,ncore
      indj = indj + 3
      if (j.ne.i) then
        do mm = 1,3
          do ii = 1,3
            do jj = 1,3
              fc3(indi+jj,indi+ii,indm+mm) = fc3(indi+jj,indi+ii,indm+mm) - fc3(indj+jj,indi+ii,indm+mm)
            enddo
          enddo
        enddo
      endif
    enddo
  enddo
  if (lshengBTE) then
!************************************
!  Count number of non-zero blocks  *
!************************************
    nfc3 = 0
    indm = -3
    do m = 1,nghostcore
      indm = indm + 3
      indi = -3
      do i = 1,ncore
        indi = indi + 3
        indj = -3
        do j = 1,ncore
          indj = indj + 3
          nnonzero = 0
          do mm = 1,3
            do ii = 1,3
              do jj = 1,3
                if (abs(fc3(indj+jj,indi+ii,indm+mm)).gt.fcmin) nnonzero = nnonzero + 1
              enddo
            enddo
          enddo
          if (nnonzero.gt.0) nfc3 = nfc3 + 1
        enddo
      enddo
    enddo
    allocate(rvghost(3_i4,nghostcell),stat=status)
    if (status/=0) call outofmemory('thirdorderfc3','rvghost')
!
!  Generate original cell prior to supercell
!
    rvtmp(1:3,1) = rv(1:3,1)/dble(nsuperghost(1,ncf))
    rvtmp(1:3,2) = rv(1:3,2)/dble(nsuperghost(2,ncf))
    rvtmp(1:3,3) = rv(1:3,3)/dble(nsuperghost(3,ncf))
!
!  Generate ghost cells
!
    mm = 0
    do ii = 0,nsuperghost(1,ncf)-1
      iim = ii - nsuperghost(1,ncf)
      if (abs(iim).ge.ii) iim = ii
      do jj = 0,nsuperghost(2,ncf)-1
        jjm = jj - nsuperghost(2,ncf)
        if (abs(jjm).ge.jj) jjm = jj
        do kk = 0,nsuperghost(3,ncf)-1
          kkm = kk - nsuperghost(3,ncf)
          if (abs(kkm).ge.kk) kkm = kk
          mm = mm + 1
          rvghost(1,mm) = iim*rvtmp(1,1) + jjm*rvtmp(1,2) + kkm*rvtmp(1,3)
          rvghost(2,mm) = iim*rvtmp(2,1) + jjm*rvtmp(2,2) + kkm*rvtmp(2,3)
          rvghost(3,mm) = iim*rvtmp(3,1) + jjm*rvtmp(3,2) + kkm*rvtmp(3,3)
        enddo
      enddo
    enddo
!************************
!  Output for ShengBTE  *
!************************
    if (nfc3.ne.0) then
!
!  Open file
!
      iout2 = iout_shengBTE + 1
      open(iout2,file='FORCE_CONSTANTS_3RD',status='unknown',form='formatted')
!
!  Write out number of third order force constants 
!
      write(iout2,'(i6)') nfc3
      nfc3 = 0
      indm = -3
      do m = 1,nghostcore
        indm = indm + 3
        indi = -3
        do i = 1,nghostcore
          do ig = 1,nghostcell
            indi = indi + 3
            indj = -3
            do j = 1,nghostcore
              do jg = 1,nghostcell
                indj = indj + 3
                nnonzero = 0
                do mm = 1,3
                  do ii = 1,3
                    do jj = 1,3
                      if (abs(fc3(indj+jj,indi+ii,indm+mm)).gt.fcmin) nnonzero = nnonzero + 1
                    enddo
                  enddo
                enddo
                if (nnonzero.gt.0) then
                  nfc3 = nfc3 + 1
!
!  Write non-zero block
!
                  write(iout2,'(/,i6)') nfc3
                  write(iout2,'(3f12.6)') (rvghost(ii,ig),ii=1,3)
                  write(iout2,'(3f12.6)') (rvghost(jj,jg),jj=1,3)
                  write(iout2,'(3i6)') m,i,j
                  do mm = 1,3
                    do ii = 1,3
                      do jj = 1,3
                        write(iout2,'(3i2,1x,f20.12)') mm,ii,jj,fc3(indj+jj,indi+ii,indm+mm)
                      enddo
                    enddo
                  enddo
                endif
              enddo
            enddo
          enddo
        enddo
      enddo
!
!  Close file
!
      close(iout2)
    endif
!
    deallocate(rvghost,stat=status)
    if (status/=0) call deallocate_error('thirdorderfc3','rvghost')
  endif
!**********************
!  Deallocate memory  *
!**********************
  deallocate(fc3,stat=status)
  if (status/=0) call deallocate_error('thirdorderfc3','fc3')
!
!  Add on CPU time less time spent projecting derivatives
!
  time2 = g_cpu_time()
  tderv3 = tderv3 + time2 - time1 
#ifdef TRACE
  call trace_out('thirdorderfc3')
#endif
!
  return
  end
