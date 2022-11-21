  subroutine thirdordern
!
!  Subroutine for calculating the thirdorder force constants.
!  Note that only the third order force constants for the
!  cores are calculated while the shells, if present, are 
!  relaxed at each point. 
!
!  Numerical version.
!
!   4/15 Created from dynamicn.f90
!   7/15 Modified so that shell optimisation can be avoided
!        at every finite difference step and be replaced by 
!        matrix operations
!   8/15 Thresholds for third order force constants added
!   8/15 rvghost vectors modified to yield the shortest 
!        equivalent distance
!   6/17 No file written if nfc3 is zero
!   2/18 Trace added
!   8/18 rvghost vectors modified to match thirdorderfc3
!
!  Assumptions:
!
!  - cores are sorted to come before shells
!  - there are no partial occupancies
!  - this is a single region calculation
!  - derv3 is already initialised to contain the diagonal 
!    blocks due to a call to phonon 
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
!  Julian Gale, CIC, Curtin University, August 2018
!
  use configurations
  use g_constants
  use control
  use current
  use derivatives
  use gulp_files
  use general,        only : phondiff
  use iochannels
  use parallel
  use shells,         only : ncore, ncsptr, nshell
  use thresholds
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
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
  integer(i4)                                      :: matom
  integer(i4)                                      :: matomg
  integer(i4)                                      :: maxlim        ! Number of variables including shells
  integer(i4)                                      :: maxlima       ! Size of first two dimensions of fc3
  integer(i4)                                      :: maxlimc       ! Number of variables for cores only
  integer(i4)                                      :: maxlims       ! Number of variables for shells only
  integer(i4)                                      :: maxlimg       ! Number of variables including shells in ghost cell
  integer(i4)                                      :: maxlimcg      ! Number of variables for cores only in ghost cell
  integer(i4)                                      :: mcrd
  integer(i4)                                      :: mm
  integer(i4)                                      :: nfc3
  integer(i4)                                      :: nghostatom    ! Number of atoms in the ghost cell
  integer(i4)                                      :: nghostcore    ! Number of cores in the ghost cell
  integer(i4)                                      :: nghostcell    ! Number of ghost cells in the full cell
  integer(i4)                                      :: nnonzero
  integer(i4)                                      :: status
  logical                                          :: lforward
  logical                                          :: lshellopt     ! If true then use shell optimisation 
  real(dp)                                         :: fcmin
  real(dp)                                         :: rstep
  real(dp)                                         :: step
  real(dp),    dimension(:,:),   allocatable, save :: Dsc
  real(dp),    dimension(:,:,:), allocatable, save :: fc3
  real(dp),    dimension(:,:),   allocatable, save :: rvghost
  real(dp)                                         :: rvtmp(3,3)
#ifdef TRACE
  call trace_in('thirdordern')
#endif
!
  fcmin = 1.0d-5
!
!  Do we optimise shells or use matrix operations?
!
  lshellopt = (index(keyword,'shop').ne.0.and.nshell.gt.0)
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
  if (lshellopt) then
    maxlima = maxlimc
  else
    maxlima = maxlim
  endif
  allocate(fc3(maxlima,maxlima,maxlimcg),stat=status)
  if (status/=0) call outofmemory('thirdordern','fc3')
!
  nfc3 = 0
  do i = 1,maxlimcg
    do j = 1,maxlima
      do k = 1,maxlima
        fc3(k,j,i) = 0.0_dp
      enddo
    enddo
  enddo
!
  rstep = 0.5_dp/phondiff
!************************************************
!  Loop over coordinates for finite difference  *
!************************************************
  do m = 1,2*maxlimcg
!
!  Find the degree of freedom and whether this is the forward or backward step
!
    if (m.gt.maxlimcg) then
      lforward = .true.
      mm = m - maxlimcg
    else
      lforward = .false.
      mm = m
    endif
!
!  Here matom and matomg refer to the core number, but particles should have 
!  been sorted so that cores come before shells and so this should be OK. 
!
    matomg = (mm-1)/3 + 1
    mcrd  = mm - 3*(matomg - 1)
!
!  Convert matomg to matom
!
    matom = (matomg - 1)*nghostcell + 1
!
!  Select step
!
    if (lforward) then
      step = phondiff
    else
      step = - phondiff
    endif
!
!  Alter coordinate
!
    if (mcrd.eq.1) then
      xalat(matom) = xalat(matom) + step
      xclat(matom) = xclat(matom) + step
    elseif (mcrd.eq.2) then
      yalat(matom) = yalat(matom) + step
      yclat(matom) = yclat(matom) + step
    elseif (mcrd.eq.3) then
      zalat(matom) = zalat(matom) + step
      zclat(matom) = zclat(matom) + step
    endif
    if (lshellopt) then
!
!  If this core has a shell then displace shell as well to avoid large separations
!
      if (ncsptr(matom).gt.0) then
        l = ncsptr(matom)
        if (mcrd.eq.1) then
          xalat(l) = xalat(l) + step
          xclat(l) = xclat(l) + step
        elseif (mcrd.eq.2) then
          yalat(l) = yalat(l) + step
          yclat(l) = yclat(l) + step
        elseif (mcrd.eq.3) then
          zalat(l) = zalat(l) + step
          zclat(l) = zclat(l) + step
        endif
      endif
!
!  Relax the shells before calculating the force constants
!
      call shellopt
    endif
!
!  Compute second order force constants
!
    call secondorderfc(lshellopt)
!
!  Add gradients to dynamical array
!
    if (lforward) then
      do i = 1,maxlima
        do j = 1,maxlima
          fc3(j,i,mm) = fc3(j,i,mm) + rstep*derv2(j,i)
        enddo
      enddo
    else
      do i = 1,maxlima
        do j = 1,maxlima
          fc3(j,i,mm) = fc3(j,i,mm) - rstep*derv2(j,i)
        enddo
      enddo
    endif
!  
!  Reverse step
!
    if (mcrd.eq.1) then
      xalat(matom) = xalat(matom) - step
      xclat(matom) = xclat(matom) - step
    elseif (mcrd.eq.2) then
      yalat(matom) = yalat(matom) - step
      yclat(matom) = yclat(matom) - step
    elseif (mcrd.eq.3) then
      zalat(matom) = zalat(matom) - step
      zclat(matom) = zclat(matom) - step
    endif
    if (lshellopt) then
!
!  If this core has a shell then displace shell as well to avoid large separations
!
      if (ncsptr(matom).gt.0) then
        l = ncsptr(matom)
        if (mcrd.eq.1) then
          xalat(l) = xalat(l) - step
          xclat(l) = xclat(l) - step
        elseif (mcrd.eq.2) then
          yalat(l) = yalat(l) - step
          yclat(l) = yclat(l) - step
        elseif (mcrd.eq.3) then
          zalat(l) = zalat(l) - step
          zclat(l) = zclat(l) - step
        endif
      endif
    endif
!*************************************
!  End loop over finite differences  *
!*************************************
  enddo
!***************************
!  Shell model correction  *
!***************************
  if (nshell.gt.0.and..not.lshellopt) then
!
!  Compute second order force constants again so that we have the correct values
!
    call secondorderfc(.false.)
!
!  Create array to store shell transformation matrix and workspace
!
    allocate(Dsc(maxlims,maxlimc),stat=status)
    if (status/=0) call outofmemory('thirdordern','Dsc')
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
    if (status/=0) call deallocate_error('thirdordern','Dsc')
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
    if (status/=0) call outofmemory('thirdordern','rvghost')
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
    if (status/=0) call deallocate_error('thirdordern','rvghost')
  endif
!**********************
!  Deallocate memory  *
!**********************
  deallocate(fc3,stat=status)
  if (status/=0) call deallocate_error('thirdordern','fc3')
#ifdef TRACE
  call trace_out('thirdordern')
#endif
!
  return
  end
