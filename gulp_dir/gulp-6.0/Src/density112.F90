  subroutine density112(imode,eatom)
!
!  Subroutine for calculating the region 1 - region 1 MEAM density
!
!  imode = 1 => defective region 1 - defective region 1 
!  imode = 3 => defective region 1 - perfect region 1
!
!  If called with modes 2 & 4 there is nothing to do.
!
!   4/09 Created from real112
!   6/09 Call to twobody1 corrected to be twoden
!   8/14 Pair potential for MEAM now evaluated here so that screening function
!        can be used.
!   8/14 Energy passed as argument so that pair potential contribution can be added
!   8/14 Call to twoden replaced by call to twoden1
!   8/14 neamspec passed psibaskes
!   5/17 Parallelisation added
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
  use parallel
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
  integer(i4)                                  :: ione
  integer(i4)                                  :: j
  integer(i4)                                  :: n
  integer(i4)                                  :: nat1
  integer(i4)                                  :: nat2
  integer(i4)                                  :: nati
  integer(i4)                                  :: natj
  integer(i4)                                  :: neamspeci
  integer(i4)                                  :: neamspecj
  integer(i4)                                  :: nloop1
  integer(i4)                                  :: nloop2
  integer(i4)                                  :: nor
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
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: cut2
  real(dp)                                     :: eatomloc
  real(dp)                                     :: ebas
  real(dp)                                     :: d1bas
  real(dp)                                     :: d2bas
  real(dp),        dimension(:,:), allocatable :: dscrholoc
  real(dp),        dimension(:,:), allocatable :: dtmp
  real(dp)                                     :: oci
  real(dp)                                     :: ocj
  real(dp)                                     :: r
  real(dp)                                     :: r2
  real(dp)                                     :: rk
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
!
  if (imode.ne.1.and.imode.ne.3) return
#ifdef TRACE
  call trace_in('density112')
#endif
!
  time1 = g_cpu_time()
!
!  Allocate local memory
!
  allocate(npotl(npote),stat=status)
  if (status/=0) call outofmemory('density112','npotl')
!
  if (imode.eq.1) then
    ione = 2
    nloop1 = nreg1
  elseif (imode.eq.2) then
    ione = 2
    nloop1 = nreg1old
  elseif (imode.eq.3) then
    ione = 1
    nloop1 = nreg1
  elseif (imode.eq.4) then
    ione = 1
    nloop1 = nreg1old
  endif
!
  eatomloc = 0.0_dp
!
  if (lMEAM) then
    allocate(dscrholoc(maxmeamcomponent,nloop1),stat=status)
    if (status/=0) call outofmemory('density112','dscrholoc')
    dscrholoc(1:maxmeamcomponent,1:nloop1) = 0.0_dp
  else
    allocate(dscrholoc(1,nloop1),stat=status)
    if (status/=0) call outofmemory('density112','dscrholoc')
    dscrholoc(1,1:nloop1) = 0.0_dp
  endif
!
!  Outer loop over sites
!
  do i = procid+ione,nloop1,nprocs
    xal = xdefe(i)
    yal = ydefe(i)
    zal = zdefe(i)
    nati = natdefe(i)
    ntypi = ntypdefe(i)
    oci = occdefe(i)
!
    neamspeci = 0
    lfound = .false.
    do while (neamspeci.lt.neamfnspec.and..not.lfound)
      neamspeci = neamspeci + 1
      if (nati.eq.neamfnnat(neamspeci).and.(ntypi.eq.neamfntyp(neamspeci).or.neamfntyp(neamspeci).eq.0)) then
        lfound = .true.
      endif
    enddo
!
    if (imode.eq.3) then
      nloop2 = nreg1old
    else
      nloop2 = i - 1
    endif
!
!  Inner loop over second site
!
    jloop: do j = 1,nloop2
      if (imode.eq.1) then
        natj = natdefe(j)
        ntypj = ntypdefe(j)
        ocj = occdefe(j)
        xcrd = xdefe(j) - xal
        ycrd = ydefe(j) - yal
        zcrd = zdefe(j) - zal
      else
        natj = natp(j)
        ntypj = ntypep(j)
        ocj = occp(j)
        xcrd = xperf(j) - xal
        ycrd = yperf(j) - yal
        zcrd = zperf(j) - zal
      endif
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
      r2 = xcrd*xcrd + ycrd*ycrd + zcrd*zcrd
      if (r2.lt.cut2) then
!
!  Store vector
!
        nor = 1
        r = sqrt(r2)
      else
        cycle jloop
      endif
!*********************************
!  Calculate unscreened density  *
!*********************************
      call twoden1(nor,1_i4,r,xcrd,ycrd,zcrd,npots,npotl,sctrm1,sctrm2,lorder12loc)
      Sij = 1.0_dp
! DEBUG
!  Screening contribution needed here
! DEBUG
!
!  For mode 3 the terms are to be subtracted  = > change sign
!
      if (imode.eq.3) then
        if (lMEAM) then
          sctrm1(1:maxmeamcomponent) = - sctrm1(1:maxmeamcomponent)
          sctrm2(1:maxmeamcomponent) = - sctrm2(1:maxmeamcomponent)
        else
          sctrm1(1) =  - sctrm1(1)
          sctrm2(1) =  - sctrm2(1)
        endif
      endif
      if (lsuttonc) then
!
!  Don't need to add terms for perfect region as rho  =  bulk rho
!
        if (lMEAM) then
          if (lorder12loc) then
            dscrho(1:maxmeamcomponent,i) = dscrho(1:maxmeamcomponent,i) + sctrm1(1:maxmeamcomponent)*ocj
            if (imode.eq.1) dscrho(1:maxmeamcomponent,j) = dscrho(1:maxmeamcomponent,j) + sctrm2(1:maxmeamcomponent)*oci
          else
            dscrho(1:maxmeamcomponent,i) = dscrho(1:maxmeamcomponent,i) + sctrm2(1:maxmeamcomponent)*ocj
            if (imode.eq.1) dscrho(1:maxmeamcomponent,j) = dscrho(1:maxmeamcomponent,j) + sctrm1(1:maxmeamcomponent)*oci
          endif
        else
          if (lorder12loc) then
            dscrho(1,i) = dscrho(1,i) + sctrm1(1)*ocj
            if (imode.eq.1) dscrho(1,j) = dscrho(1,j) + sctrm2(1)*oci
          else
            dscrho(1,i) = dscrho(1,i) + sctrm2(1)*ocj
            if (imode.eq.1) dscrho(1,j) = dscrho(1,j) + sctrm1(1)*oci
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
    enddo jloop
  enddo
  if (nprocs.gt.1) then
!
!  Globalise energy
!
    call sumone(eatomloc,ebas,"density112","eatom")
    eatom = eatom + ebas
!
!  Globalise densities
!
    if (lMEAM) then
      allocate(dtmp(maxmeamcomponent,nloop1),stat=status)
      if (status/=0) call outofmemory('density112','dtmp')
      call sumall(dscrholoc,dtmp,maxmeamcomponent*nreg1,"density12a2","dscrho")
      dscrho(1:maxmeamcomponent,1:nreg1) = dscrho(1:maxmeamcomponent,1:nreg1) + dtmp(1:maxmeamcomponent,1:nreg1)
    else
      allocate(dtmp(1,nloop1),stat=status)
      if (status/=0) call outofmemory('density112','dtmp')
      call sumall(dscrholoc,dtmp,nloop1,"density112","dscrho")
      dscrho(1,1:nloop1) = dscrho(1,1:nloop1) + dtmp(1,1:nloop1)
    endif
    deallocate(dtmp,stat=status)
    if (status/=0) call deallocate_error('density112','dtmp')
  else
    eatom = eatom + eatomloc
  endif
!
!  Free local memory
!
  deallocate(dscrholoc,stat=status)
  if (status/=0) call deallocate_error('density112','dscrholoc')
  deallocate(npotl,stat=status)
  if (status/=0) call deallocate_error('density112','npotl')
!
!  Timing
!
  time2 = g_cpu_time()
  treg1 = treg1 + time2 - time1
#ifdef TRACE
  call trace_out('density112')
#endif
!
  return
  end
