  subroutine density12as2(eatom)
!
!  Subroutine for calculating defect MEAM density
!
!  xdefe = region 1 coordinates
!  xclat = region 2 coordinates per unit cell
!
!  imode = 1 => defective region 1 - region 2
!
!   4/09 Created from real12as2
!   6/09 Calls to twobody1 corrected to be twoden
!   8/14 Pair potential for MEAM now evaluated here so that screening function
!        can be used.
!   8/14 Energy passed as argument so that pair potential contribution can be added
!   8/14 Calls to twoden replaced by calls to twoden1
!   8/14 neamspec passed psibaskes
!   2/18 Trace added
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
  use control
  use current
  use defects
  use eam
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
  real(dp),                     intent(inout)  :: eatom
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ii
  integer(i4)                                  :: j
  integer(i4)                                  :: n
  integer(i4)                                  :: nat1
  integer(i4)                                  :: nat2
  integer(i4)                                  :: nati
  integer(i4)                                  :: natj
  integer(i4)                                  :: neamspeci
  integer(i4)                                  :: neamspecj
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
  real(dp)                                     :: d1bas
  real(dp)                                     :: d2bas
  real(dp)                                     :: ocj
  real(dp)                                     :: r
  real(dp)                                     :: r2
  real(dp)                                     :: rcheck
  real(dp)                                     :: rcheck2
  real(dp)                                     :: rdiffc
  real(dp)                                     :: rk
  real(dp)                                     :: rmiddle2
  real(dp)                                     :: rp
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
  call trace_in('density12as2')
#endif
!
  time1 = g_cpu_time()
!
!  Allocate local memory
!
  allocate(npotl(npote),stat=status)
  if (status/=0) call outofmemory('real12as2','npotl')
!
!  Decide whether all interactions lie within unit cell and first
!  neighbours - saves time for large systems
!
  cmax = rpmax
!**********************************
!  Region 1  -  region 1/2a energy  *
!**********************************
!
!  Outer loop over region1 asymmetric unit sites
!
  do i = 1,ndasym
    ii = ndsptr(i)
    xal = xdefe(ii)
    yal = ydefe(ii)
    zal = zdefe(ii)
    nati = natdefe(ii)
    ntypi = ntypdefe(ii)
!
!  Find distance from the centre
!
    xdiffc = xal - xdc
    ydiffc = yal - ydc
    zdiffc = zal - zdc
    rdiffc = xdiffc*xdiffc + ydiffc*ydiffc+zdiffc*zdiffc
    rdiffc = sqrt(rdiffc)
    rcheck = rdiffc + cmax+0.5_dp
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
          if (lMEAM) then
            if (lorder12loc) then
              dscrho(1:maxmeamcomponent,i) = dscrho(1:maxmeamcomponent,i) + sctrm1(1:maxmeamcomponent)*ocj
            else
              dscrho(1:maxmeamcomponent,i) = dscrho(1:maxmeamcomponent,i) + sctrm2(1:maxmeamcomponent)*ocj
            endif
          else
            if (lorder12loc) then
              dscrho(1,i) = dscrho(1,i) + sctrm1(1)*ocj
            else
              dscrho(1,i) = dscrho(1,i) + sctrm2(1)*ocj
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
              eatom = eatom + ebas*Sij
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
      ocj = or2a(j)
!
!  Distance check relative to centre of region 1
!
      xdiffc = xr2a(j) - xdc
      ydiffc = yr2a(j) - ydc
      zdiffc = zr2a(j) - zdc
      rmiddle2 = xdiffc*xdiffc + ydiffc*ydiffc+zdiffc*zdiffc
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
!  Need to divide sctrm by symmetry factor
!
          if (lMEAM) then
            if (lorder12loc) then
              dscrho(1:maxmeamcomponent,i) = dscrho(1:maxmeamcomponent,i) + sctrm1(1:maxmeamcomponent)*ocj
            else
              dscrho(1:maxmeamcomponent,i) = dscrho(1:maxmeamcomponent,i) + sctrm2(1:maxmeamcomponent)*ocj
            endif
          else
            if (lorder12loc) then
              dscrho(1,i) = dscrho(1,i) + sctrm1(1)*ocj
            else
              dscrho(1,i) = dscrho(1,i) + sctrm2(1)*ocj
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
              eatom = eatom + ebas*Sij
            endif
          endif
        enddo
      endif
10    continue
    enddo
  enddo
!
!  End of loop over region 2
!
!
!  Free local memory
!
  deallocate(npotl,stat=status)
  if (status/=0) call deallocate_error('real12as2','npotl')
!
!  Timing
!
  time2 = g_cpu_time()
  treg1 = treg1 + time2 - time1
#ifdef TRACE
  call trace_out('density12as2')
#endif
!
  return
  end
