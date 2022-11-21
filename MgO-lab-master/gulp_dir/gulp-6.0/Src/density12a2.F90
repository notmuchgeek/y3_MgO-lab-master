  subroutine density12a2(imode,eatom)
!
!  Subroutine for calculating defect MEAM density
!
!  xdefe = region 1 coordinates
!  xclat = region 2 coordinates per unit cell
!
!  imode = 1 => defective region 1 - region 2
!  imode = 2 => perfect region 1 - region 2
!
!   4/09 Created from real12a2
!   8/14 Pair potential for MEAM now evaluated here so that screening function
!        can be used.
!   8/14 Energy passed as argument so that pair potential contribution can be added
!   8/14 neamspec passed psibaskes
!   5/17 Parallel modifications added
!   2/18 Trace added
!   7/20 Separate routine for sumall with 1 argument added
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
!  Copyright Curtin University 2020
!
!  Julian Gale, CIC, Curtin University, July 2020
!
  use control
  use current
  use defects
  use eam
  use general,        only : smallself
  use parallel
  use region2a
  use sutton
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  use two
  implicit none
!
!  Passed variables
!
  integer(i4),                  intent(in)     :: imode
  real(dp),                     intent(inout)  :: eatom
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: j
  integer(i4)                                  :: n
  integer(i4)                                  :: nat1
  integer(i4)                                  :: nat2
  integer(i4)                                  :: nati
  integer(i4)                                  :: natj
  integer(i4)                                  :: neamspeci
  integer(i4)                                  :: neamspecj
  integer(i4)                                  :: nloop
  integer(i4)                                  :: np
  integer(i4)                                  :: npot
  integer(i4), dimension(:), allocatable       :: npotl
  integer(i4)                                  :: npots
  integer(i4)                                  :: ntyp1
  integer(i4)                                  :: ntyp2
  integer(i4)                                  :: ntypi
  integer(i4)                                  :: ntypj
  integer(i4)                                  :: status
  logical                                      :: lfound
  logical                                      :: lmatch
  logical                                      :: lorder12loc
  real(dp)                                     :: cmax
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: cut2
  real(dp)                                     :: ebas
  real(dp)                                     :: eatomloc
  real(dp)                                     :: d1bas
  real(dp)                                     :: d2bas
  real(dp),        dimension(:,:), allocatable :: dscrholoc
  real(dp),        dimension(:,:), allocatable :: dscrhor2dloc
  real(dp),        dimension(:,:), allocatable :: dscrhor2ploc
  real(dp),        dimension(:,:), allocatable :: dtmp
  real(dp)                                     :: oci
  real(dp)                                     :: ocj
  real(dp)                                     :: r
  real(dp)                                     :: r2
  real(dp)                                     :: rcheck
  real(dp)                                     :: rcheck2
  real(dp)                                     :: rdiffc
  real(dp)                                     :: rk
  real(dp)                                     :: rmiddle2
  real(dp)                                     :: rp
  real(dp)                                     :: rsign
  real(dp)                                     :: Sij
  real(dp)                                     :: sctrm1(maxmeamcomponent)
  real(dp)                                     :: sctrm2(maxmeamcomponent)
  real(dp)                                     :: time1
  real(dp)                                     :: time2
  real(dp)                                     :: xal
  real(dp)                                     :: yal
  real(dp)                                     :: zal
  real(dp)                                     :: xcrd
  real(dp)                                     :: ycrd
  real(dp)                                     :: zcrd
  real(dp)                                     :: xdiffc
  real(dp)                                     :: ydiffc
  real(dp)                                     :: zdiffc
#ifdef TRACE
  call trace_in('density12a2')
#endif
!
  time1 = g_cpu_time()
!
!  Allocate local memory
!
  allocate(npotl(npote),stat=status)
  if (status/=0) call outofmemory('density12a2','npotl')
!
!  Decide whether all interactions lie within unit cell and first
!  neighbours - saves time for large systems
!
  cmax = rpmax
!
  if (imode.eq.1) then
    nloop = nreg1
  elseif (imode.eq.2) then
    nloop = nreg1old
  endif
!
  eatomloc = 0.0_dp
!
  if (lMEAM) then
    if (imode.eq.1) then
      allocate(dscrholoc(maxmeamcomponent,nreg1),stat=status)
      if (status/=0) call outofmemory('density12a2','dscrholoc')
      allocate(dscrhor2dloc(maxmeamcomponent,npreg2),stat=status)
      if (status/=0) call outofmemory('density12a2','dscrhor2dloc')
      dscrholoc(1:maxmeamcomponent,1:nreg1) = 0.0_dp
      dscrhor2dloc(1:maxmeamcomponent,1:npreg2) = 0.0_dp
    else
      allocate(dscrhor2ploc(maxmeamcomponent,npreg2),stat=status)
      if (status/=0) call outofmemory('density12a2','dscrhor2ploc')
      dscrhor2ploc(1:maxmeamcomponent,1:npreg2) = 0.0_dp
    endif
  else
    if (imode.eq.1) then
      allocate(dscrholoc(1,nreg1),stat=status)
      if (status/=0) call outofmemory('density12a2','dscrholoc')
      allocate(dscrhor2dloc(1,npreg2),stat=status)
      if (status/=0) call outofmemory('density12a2','dscrhor2dloc')
      dscrholoc(1,1:nreg1) = 0.0_dp
      dscrhor2dloc(1,1:npreg2) = 0.0_dp
    else
      allocate(dscrhor2ploc(1,npreg2),stat=status)
      if (status/=0) call outofmemory('density12a2','dscrhor2ploc')
      dscrhor2ploc(1,1:npreg2) = 0.0_dp
    endif
  endif
!**********************************
!  Region 1 - region 1/2a energy  *
!**********************************
!
!  Outer loop over sites
!
  do i = procid+1,nloop,nprocs
!
!  Set i attributes according to whether it is defective or perfect region 1
!
    if (imode.eq.1) then
      xal = xdefe(i)
      yal = ydefe(i)
      zal = zdefe(i)
      nati = natdefe(i)
      ntypi = ntypdefe(i)
      oci = occdefe(i)
    elseif (imode.eq.2) then
      xal = xperf(i)
      yal = yperf(i)
      zal = zperf(i)
      nati = natp(i)
      ntypi = ntypep(i)
      oci = occp(i)
    endif
!
!  Find distance from the centre
!
    xdiffc = xal - xdc
    ydiffc = yal - ydc
    zdiffc = zal - zdc
    rdiffc = xdiffc*xdiffc + ydiffc*ydiffc + zdiffc*zdiffc
    rdiffc = sqrt(rdiffc)
    rcheck = rdiffc + cmax + 0.5_dp
    rcheck2 = rcheck*rcheck
!
    neamspeci = 0
    lfound = .false.
    do while (neamspeci.lt.neamfnspec.and..not.lfound)
      neamspeci = neamspeci + 1
      if (nati.eq.neamfnnat(neamspeci).and.(ntypi.eq.neamfntyp(neamspeci).or.neamfntyp(neamspeci).eq.0)) then
        lfound = .true.
      endif
    enddo
!***************************
!  Loop over old region 1  *
!***************************
    do j = 1,nreg1old
      natj = natp(j)
      ntypj = ntypep(j)
!
!  Arrange nat1 to be lower than nat2. If they are equal
!  then arrange that ntyp1 is lower than ntyp2
!
      if (nati.eq.natj) then
        nat1 = nati
        nat2 = natj
        if (ntypi.lt.ntypj) then
          lorder12loc = .true.
          ntyp1 = ntypi
          ntyp2 = ntypj
        else
          lorder12loc = .false.
          ntyp1 = ntypj
          ntyp2 = ntypi
        endif
      elseif (nati.lt.natj) then
        lorder12loc = .true.
        nat1 = nati
        nat2 = natj
        ntyp1 = ntypi
        ntyp2 = ntypj
      else
        lorder12loc = .false.
        nat1 = natj
        nat2 = nati
        ntyp1 = ntypj
        ntyp2 = ntypi
      endif
      xcrd = xperf(j) - xal
      ycrd = yperf(j) - yal
      zcrd = zperf(j) - zal
      ocj = occp(j)
!
      neamspecj = 0
      lfound = .false.
      do while (neamspecj.lt.neamfnspec.and..not.lfound)
        neamspecj = neamspecj + 1
        if (natj.eq.neamfnnat(neamspecj).and.(ntypj.eq.neamfntyp(neamspecj).or.neamfntyp(neamspecj).eq.0)) then
          lfound = .true.
        endif
      enddo
!
!  Locate potential number
!  Check whether potential requires specific types
!
      rp = 0.0_dp
      npots = 0
      do n = 1,npote
        if (nptype(n).eq.19.or.nptype(n).eq.45.or.nptype(n).eq.55) then
          if (lmatch(nat1,ntyp1,nspec1(n),nptyp1(n),.true.)) then
            if (lmatch(nat2,ntyp2,nspec2(n),nptyp2(n),.true.)) then
              npots = npots + 1
              npotl(npots)= n
              if (rpot(n).gt.rp) rp = rpot(n)
            endif
          endif
        endif
      enddo
!
!  Generate looping indices
!
      cut2 = rp*rp
!
!  Generate distance squared
!
      r2 = xcrd*xcrd + ycrd*ycrd + zcrd*zcrd
      if (r2.le.cut2.and.r2.gt.smallself) then
!
!  Store vector
!
        r = sqrt(r2)
!
!  Compute unscreened density
!
        call twoden1(1_i4,1_i4,r,xcrd,ycrd,zcrd,npots,npotl,sctrm1,sctrm2,lorder12loc)
        Sij = 1.0_dp
! DEBUG
!  Screening contribution needed here
! DEBUG
!
!  For mode 2 the terms are to be subtracted  = > change sign
!
        if (imode.eq.2) then
          if (lMEAM) then
            sctrm1(1:maxmeamcomponent) = - sctrm1(1:maxmeamcomponent)
            sctrm2(1:maxmeamcomponent) = - sctrm2(1:maxmeamcomponent)
          else
            sctrm1(1) = - sctrm1(1)
            sctrm2(1) = - sctrm2(1)
          endif
        endif
        if (lsuttonc.and.imode.eq.1) then
!
!  Don't need to add terms for perfect region as rho  =  bulk rho
!
          if (lMEAM) then
            if (lorder12loc) then
              dscrholoc(1:maxmeamcomponent,i) = dscrholoc(1:maxmeamcomponent,i) + sctrm1(1:maxmeamcomponent)*ocj
            else
              dscrholoc(1:maxmeamcomponent,i) = dscrholoc(1:maxmeamcomponent,i) + sctrm2(1:maxmeamcomponent)*ocj
            endif
          else
            if (lorder12loc) then
              dscrholoc(1,i) = dscrholoc(1,i) + sctrm1(1)*ocj
            else
              dscrholoc(1,i) = dscrholoc(1,i) + sctrm2(1)*ocj
            endif
          endif
        endif
        if (imode.eq.2) then
          rsign = - 1.0_dp
        else
          rsign =   1.0_dp
        endif
!
!  Compute pairwise contribution to energy
!
        do np = 1,npots
          npot = npotl(np)
          if (nptype(npot).eq.45.or.nptype(npot).eq.55) then
            if (r.gt.rpot2(npot).and.r.le.rpot(npot)) then
              rk = 1.0_dp/r
              call baskes(nati,ntypi,neamspeci,natj,ntypj,neamspecj,npot,r,rk,1.0_dp,ebas,d1bas,d2bas,.false.,.false.)
              eatomloc = eatomloc + rsign*ebas*Sij
            endif
          endif
        enddo
      endif
    enddo
!************************
!  Loop over region 2a  *
!************************
    rmiddle2 = 0.0_dp
    j = 0
    do while (j.lt.npreg2.and.rmiddle2.lt.rcheck2)
      j = j + 1
      xcrd = xr2a(j) - xal
      ycrd = yr2a(j) - yal
      zcrd = zr2a(j) - zal
!
!  If not a molecular calc and component exceeds maximum
!  cut - off then there is nothing to evaluate
!
      if (abs(xcrd).gt.cmax) goto 10
      if (abs(ycrd).gt.cmax) goto 10
      if (abs(zcrd).gt.cmax) goto 10
!
      natj = nr2a(j)
      ntypj = ntr2a(j)
!
!  Arrange nat1 to be lower than nat2. If they are equal
!  then arrange that ntyp1 is lower than ntyp2
!
      if (nati.eq.natj) then
        lorder12loc = .true.
        nat1 = nati
        nat2 = natj
        if (ntypi.lt.ntypj) then
          ntyp1 = ntypi
          ntyp2 = ntypj
        else
          ntyp1 = ntypj
          ntyp2 = ntypi
        endif
      elseif (nati.lt.natj) then
        lorder12loc = .true.
        nat1 = nati
        nat2 = natj
        ntyp1 = ntypi
        ntyp2 = ntypj
      else
        lorder12loc = .false.
        nat1 = natj
        nat2 = nati
        ntyp1 = ntypj
        ntyp2 = ntypi
      endif
      ocj = or2a(j)
!
!  Distance check relative to centre of region 1
!
      xdiffc = xr2a(j) - xdc
      ydiffc = yr2a(j) - ydc
      zdiffc = zr2a(j) - zdc
      rmiddle2 = xdiffc*xdiffc + ydiffc*ydiffc + zdiffc*zdiffc
!
      neamspecj = 0
      lfound = .false.
      do while (neamspecj.lt.neamfnspec.and..not.lfound)
        neamspecj = neamspecj + 1
        if (natj.eq.neamfnnat(neamspecj).and.(ntypj.eq.neamfntyp(neamspecj).or.neamfntyp(neamspecj).eq.0)) then
          lfound = .true.
        endif
      enddo
!
!  Locate potential number
!  Check whether potential requires specific types
!
      rp = 0.0_dp
      npots = 0
      do n = 1,npote
        if (nptype(n).eq.19.or.nptype(n).eq.45.or.nptype(n).eq.55) then
          if (lmatch(nat1,ntyp1,nspec1(n),nptyp1(n),.true.)) then
            if (lmatch(nat2,ntyp2,nspec2(n),nptyp2(n),.true.)) then
              npots = npots + 1
              npotl(npots) = n
              if (rpot(n).gt.rp) rp = rpot(n)
            endif
          endif
        endif
      enddo
!
!  Generate looping indices
!
      cut2 = rp*rp
!
!  Generate distance squared
!
      r2 = xcrd*xcrd + ycrd*ycrd + zcrd*zcrd
      if (r2.le.cut2) then
!
!  Store vector
!
        r = sqrt(r2)
!
!  Compute unscreened density
!
        call twoden1(1_i4,1_i4,r,xcrd,ycrd,zcrd,npots,npotl,sctrm1,sctrm2,lorder12loc)
        Sij = 1.0_dp
! DEBUG
!  Screening contribution needed here
! DEBUG
!
        if (lsuttonc) then
!
!  Don't need to add terms for perfect region as rho = bulk rho
!
          if (lMEAM) then
            if (lorder12loc) then
              if (imode.eq.1) then
                dscrholoc(1:maxmeamcomponent,i) = dscrholoc(1:maxmeamcomponent,i) + sctrm1(1:maxmeamcomponent)*ocj
                dscrhor2dloc(1:maxmeamcomponent,j) = dscrhor2dloc(1:maxmeamcomponent,j) + sctrm2(1:maxmeamcomponent)*oci
              else
                dscrhor2ploc(1:maxmeamcomponent,j) = dscrhor2ploc(1:maxmeamcomponent,j) - sctrm2(1:maxmeamcomponent)*oci
              endif
            else
              if (imode.eq.1) then
                dscrholoc(1:maxmeamcomponent,i) = dscrholoc(1:maxmeamcomponent,i) + sctrm2(1:maxmeamcomponent)*ocj
                dscrhor2dloc(1:maxmeamcomponent,j) = dscrhor2dloc(1:maxmeamcomponent,j) + sctrm1(1:maxmeamcomponent)*oci
              else
                dscrhor2ploc(1:maxmeamcomponent,j) = dscrhor2ploc(1:maxmeamcomponent,j) - sctrm1(1:maxmeamcomponent)*oci
              endif
            endif
          else
            if (lorder12loc) then
              if (imode.eq.1) then
                dscrholoc(1,i) = dscrholoc(1,i) + sctrm1(1)*ocj
                dscrhor2dloc(1,j) = dscrhor2dloc(1,j) + sctrm2(1)*oci
              else
                dscrhor2ploc(1,j) = dscrhor2ploc(1,j) - sctrm2(1)*oci
              endif
            else
              if (imode.eq.1) then
                dscrholoc(1,i) = dscrholoc(1,i) + sctrm2(1)*ocj
                dscrhor2dloc(1,j) = dscrhor2dloc(1,j) + sctrm1(1)*oci
              else
                dscrhor2ploc(1,j) = dscrhor2ploc(1,j) - sctrm1(1)*oci
              endif
            endif
          endif
        endif
!
!  Compute pairwise contribution to energy
!
        do np = 1,npots
          npot = npotl(np)
          if (nptype(npot).eq.45.or.nptype(npot).eq.55) then
            if (r.gt.rpot2(npot).and.r.le.rpot(npot)) then
              rk = 1.0_dp/r
              call baskes(nati,ntypi,neamspeci,natj,ntypj,neamspecj,npot,r,rk,1.0_dp,ebas,d1bas,d2bas,.false.,.false.)
              eatomloc = eatomloc + ebas*Sij
            endif
          endif
        enddo
      endif
10    continue
    enddo
  enddo
  if (nprocs.gt.1) then
!
!  Globalise energy
!
    call sumone(eatomloc,ebas,"density12a2","eatom")
    eatom = eatom + ebas
!
!  Globalise densities
!
    if (lMEAM) then
      allocate(dtmp(maxmeamcomponent,max(nreg1,npreg2)),stat=status)
      if (status/=0) call outofmemory('density12a2','dtmp')
      if (imode.eq.1) then
        call sumall(dscrholoc,dtmp,maxmeamcomponent*nreg1,"density12a2","dscrho")
        dscrho(1:maxmeamcomponent,1:nreg1) = dscrho(1:maxmeamcomponent,1:nreg1) + dtmp(1:maxmeamcomponent,1:nreg1)
        call sumall(dscrhor2dloc,dtmp,maxmeamcomponent*npreg2,"density12a2","dscrhor2d")
        dscrhor2d(1:maxmeamcomponent,1:npreg2) = dscrhor2d(1:maxmeamcomponent,1:npreg2) + dtmp(1:maxmeamcomponent,1:npreg2)
      else
        call sumall(dscrhor2ploc,dtmp,maxmeamcomponent*npreg2,"density12a2","dscrhor2p")
        dscrhor2p(1:maxmeamcomponent,1:npreg2) = dscrhor2p(1:maxmeamcomponent,1:npreg2) + dtmp(1:maxmeamcomponent,1:npreg2)
      endif
    else
      allocate(dtmp(1,max(nreg1,npreg2)),stat=status)
      if (status/=0) call outofmemory('density12a2','dtmp')
      if (imode.eq.1) then
        call sumall(dscrholoc,dtmp,nreg1,"density12a2","dscrho")
        dscrho(1,1:nreg1) = dscrho(1,1:nreg1) + dtmp(1,1:nreg1)
        call sumall(dscrhor2dloc,dtmp,npreg2,"density12a2","dscrhor2d")
        dscrhor2d(1,1:npreg2) = dscrhor2d(1,1:npreg2) + dtmp(1,1:npreg2)
      else
        call sumall(dscrhor2ploc,dtmp,npreg2,"density12a2","dscrhor2p")
        dscrhor2p(1,1:npreg2) = dscrhor2p(1,1:npreg2) + dtmp(1,1:npreg2)
      endif
    endif
    deallocate(dtmp,stat=status)
    if (status/=0) call deallocate_error('density12a2','dtmp')
  else
    eatom = eatom + eatomloc
    if (lMEAM) then
      if (imode.eq.1) then
        dscrho(1:maxmeamcomponent,1:nreg1) = dscrho(1:maxmeamcomponent,1:nreg1) + dscrholoc(1:maxmeamcomponent,1:nreg1)
        dscrhor2d(1:maxmeamcomponent,1:npreg2) = dscrhor2d(1:maxmeamcomponent,1:npreg2) + dscrhor2dloc(1:maxmeamcomponent,1:npreg2)
      else
        dscrhor2p(1:maxmeamcomponent,1:npreg2) = dscrhor2p(1:maxmeamcomponent,1:npreg2) + dscrhor2ploc(1:maxmeamcomponent,1:npreg2)
      endif
    else
      if (imode.eq.1) then
        dscrho(1,1:nreg1) = dscrho(1,1:nreg1) + dscrholoc(1,1:nreg1)
        dscrhor2d(1,1:npreg2) = dscrhor2d(1,1:npreg2) + dscrhor2dloc(1,1:npreg2)
      else
        dscrhor2p(1,1:npreg2) = dscrhor2p(1,1:npreg2) + dscrhor2ploc(1,1:npreg2)
      endif
    endif
  endif
!
!  Free local memory
!
  if (imode.eq.1) then
    deallocate(dscrhor2dloc,stat=status)
    if (status/=0) call deallocate_error('density12a2','dscrhor2dloc')
    deallocate(dscrholoc,stat=status)
    if (status/=0) call deallocate_error('density12a2','dscrholoc')
  else
    deallocate(dscrhor2ploc,stat=status)
    if (status/=0) call deallocate_error('density12a2','dscrhor2ploc')
  endif
  deallocate(npotl,stat=status)
  if (status/=0) call deallocate_error('density12a2','npotl')
!
!  Timing
!
  time2 = g_cpu_time()
  treg1 = treg1 + time2 - time1
#ifdef TRACE
  call trace_out('density12a2')
#endif
!
  return
  end
