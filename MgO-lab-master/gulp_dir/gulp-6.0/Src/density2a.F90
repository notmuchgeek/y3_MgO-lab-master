  subroutine density2a(eatom)
!
!  Subroutine for calculating the contribution of the defective
!  region 1 to the region 2a asymmetric unit densities in the MEAM model.
!
!  This routine is only needed when the symmetry adapted
!  algorithm is being used for the defect calculation.
!
!   4/09 Created from real2amany
!   8/14 Pair potential for MEAM now evaluated here so that screening function
!        can be used.
!   8/14 Energy passed as argument so that pair potential contribution can be added
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
  use eam,           only : lMEAM, maxmeamcomponent
  use region2a
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
  integer(i4)                                  :: j
  integer(i4)                                  :: jj
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
  logical                                      :: lmatch
  logical                                      :: lorder12loc
  real(dp)                                     :: cmax
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: cut2
  real(dp)                                     :: ebas
  real(dp)                                     :: d1bas
  real(dp)                                     :: d2bas
  real(dp)                                     :: oci
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
  call trace_in('density2a')
#endif
!
  time1 = g_cpu_time()
!
!  Allocate local memory
!
  allocate(npotl(npote),stat=status)
  if (status/=0) call outofmemory('real2amany','npotl')
!
  cmax = rpmax
!*************************************************
!  Region 1 - region 1/2a density contributions  *
!*************************************************
!
!  Outer loop over sites
!
  do i = 1,nreg1
    xal = xdefe(i)
    yal = ydefe(i)
    zal = zdefe(i)
    nati = natdefe(i)
    ntypi = ntypdefe(i)
    oci = occdefe(i)
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
!************************
!  Loop over region 2a  *
!************************
    rmiddle2 = 0.0_dp
    j = 0
    do while (j.lt.ndpasym2a.and.rmiddle2.lt.rcheck2)
      j = j + 1
      jj = ndsptr2a(j)
      xcrd = xr2a(jj) - xal
      ycrd = yr2a(jj) - yal
      zcrd = zr2a(jj) - zal
!
!  If not a molecular calc and component exceeds maximum
!  cut-off then there is nothing to evaluate
!
      if (abs(xcrd).gt.cmax) cycle
      if (abs(ycrd).gt.cmax) cycle
      if (abs(zcrd).gt.cmax) cycle
      natj = nr2a(jj)
      ntypj = ntr2a(jj)
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
!
!  Distance check relative to centre of region 1
!
      xdiffc = xr2a(jj) - xdc
      ydiffc = yr2a(jj) - ydc
      zdiffc = zr2a(jj) - zdc
      rmiddle2 = xdiffc*xdiffc + ydiffc*ydiffc + zdiffc*zdiffc
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
!*********************************
!  Calculate unscreened density  *
!*********************************
        call twoden1(1_i4,1_i4,r,xcrd,ycrd,zcrd,npots,npotl,sctrm1,sctrm2,lorder12loc)
        Sij = 1.0_dp
! DEBUG
!  Screening contribution needed here
! DEBUG
!
        if (lMEAM) then
          if (lorder12loc) then
            dscrhor2d(1:maxmeamcomponent,j) = dscrhor2d(1:maxmeamcomponent,j) + sctrm2(1:maxmeamcomponent)*oci
          else
            dscrhor2d(1:maxmeamcomponent,j) = dscrhor2d(1:maxmeamcomponent,j) + sctrm1(1:maxmeamcomponent)*oci
          endif
        else
          if (lorder12loc) then
            dscrhor2d(1,j) = dscrhor2d(1,j) + sctrm2(1)*oci
          else
            dscrhor2d(1,j) = dscrhor2d(1,j) + sctrm1(1)*oci
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
  enddo
!
!  End of atom loops
!
!
!  Free local memory
!
  deallocate(npotl,stat=status)
  if (status/=0) call deallocate_error('real2amany','npotl')
!
!  Timing
!
  time2 = g_cpu_time()
  tregm = tregm + time2 - time1
#ifdef TRACE
  call trace_out('density2a')
#endif
!
  return
  end
