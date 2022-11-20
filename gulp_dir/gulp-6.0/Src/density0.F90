  subroutine density0(eatom)
!
!  Subroutine for calculating MEAM electron density, including the screening function
!
!   4/09 Created from realmd0
!   4/09 Screening function added
!   8/14 MEAM screening made species specific
!   8/14 Energy passed as argument so that pair potential contribution can be added
!   8/14 All MEAM species (i,j,k) now passed to meanscreen
!   8/14 neamspec passed psibaskes
!   2/18 Trace added
!   5/19 Site energy added
!   5/19 Twobody energy print out added
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
!  Copyright Curtin University 2019
!
!  Julian Gale, CIC, Curtin University, May 2019
!
  use configurations, only : nregionno, nregiontype, QMMMmode
  use control
  use current
  use eam,            only : lMEAM, lMEAMscreen, maxmeamcomponent, meam_Cmax
  use eam,            only : lanyMEAMscreen, neamfnspec
  use element
  use energies,       only : siteenergy
  use general,        only : smallself
  use iochannels,     only : ioout
  use optimisation
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
  real(dp),                     intent(inout)  :: eatom
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: indij
  integer(i4)                                  :: j
  integer(i4)                                  :: k
  integer(i4)                                  :: m
  integer(i4)                                  :: n
  integer(i4)                                  :: nat1
  integer(i4)                                  :: nat2
  integer(i4)                                  :: nati
  integer(i4)                                  :: natj
  integer(i4)                                  :: neamfnspec2
  integer(i4)                                  :: neamspeci
  integer(i4)                                  :: neamspecj
  integer(i4)                                  :: neamspeck
  integer(i4)                                  :: noff
  integer(i4)                                  :: noffm1
  integer(i4)                                  :: noffset
  integer(i4)                                  :: np
  integer(i4), dimension(:), allocatable       :: npotl
  integer(i4)                                  :: npot
  integer(i4)                                  :: npots
  integer(i4)                                  :: nregioni
  integer(i4)                                  :: nregionj
  integer(i4)                                  :: nregiontypi
  integer(i4)                                  :: nregiontypj
  integer(i4)                                  :: ntyp1
  integer(i4)                                  :: ntyp2
  integer(i4)                                  :: ntypi
  integer(i4)                                  :: ntypj
  integer(i4)                                  :: status
  logical                                      :: lmatch
  logical                                      :: lnonzeroSij
  logical                                      :: lopi
  logical                                      :: lopj
  logical                                      :: lorder12loc
  logical                                      :: lpartial
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: cut2
  real(dp)                                     :: dSikjdr(3)      ! Dummy arguments for meamscreen call - not used here
  real(dp)                                     :: ebas
  real(dp)                                     :: ebassum
  real(dp)                                     :: d1bas
  real(dp)                                     :: d2bas
  real(dp)                                     :: oci
  real(dp)                                     :: ocj
  real(dp)                                     :: r2
  real(dp)                                     :: r2ijmid
  real(dp)                                     :: r2ik
  real(dp)                                     :: r2jk
  real(dp)                                     :: rcut2
  real(dp)                                     :: rcutfactor
  real(dp)                                     :: r
  real(dp)                                     :: rk
  real(dp)                                     :: rp
  real(dp)                                     :: sctrm1(maxmeamcomponent)
  real(dp)                                     :: sctrm2(maxmeamcomponent)
  real(dp)                                     :: Sij
  real(dp)                                     :: Sikj
  real(dp),    dimension(:), allocatable       :: sum
  real(dp),    dimension(:), allocatable       :: sum2
  real(dp)                                     :: time1
  real(dp)                                     :: time2
  real(dp)                                     :: tsum0
  real(dp)                                     :: tsuml
  real(dp)                                     :: xal
  real(dp)                                     :: yal
  real(dp)                                     :: zal
  real(dp)                                     :: xcrd
  real(dp)                                     :: ycrd
  real(dp)                                     :: zcrd
  real(dp)                                     :: xij0
  real(dp)                                     :: yij0
  real(dp)                                     :: zij0
  real(dp)                                     :: xik0
  real(dp)                                     :: yik0
  real(dp)                                     :: zik0
  real(dp)                                     :: xjk0
  real(dp)                                     :: yjk0
  real(dp)                                     :: zjk0
!
!  Check that this routine should have been called
!
  if (.not.lsuttonc) return
#ifdef TRACE
  call trace_in('density0')
#endif
!
  time1 = g_cpu_time()
  tsuml = 0.0_dp
!
!  For screened MEAM, set scale factor that determines searching cutoff based on which axis of the ellipse is largest.
!  Note: the cutoff is applied to the mid point of the i-j vector and so this is guaranteed to find all distances.
!
  rcutfactor = 1.0_dp
  if (lanyMEAMscreen) then
    neamfnspec2 = neamfnspec*(neamfnspec+1)/2
    do i = 1,neamfnspec
      do j = 1,neamfnspec2
        if (lMEAMscreen(j,i)) then
          rcutfactor = max(1.0_dp,0.25_dp*(1.0_dp + meam_Cmax(j,i)))
        endif
      enddo
    enddo
  endif
!
!  Opening banner for energy decomposition
!
  if (lPrintTwo) then
    call mpbarrier
    if (ioproc) then
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''  Two : Atom No. 1  Atom No. 2    Baskes short-range energy (eV) '')')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
    endif
  endif
!
!  Allocate local memory
!
  allocate(npotl(npote),stat=status)
  if (status/=0) call outofmemory('realmd0','npotl')
  allocate(sum(numat),stat=status)
  if (status/=0) call outofmemory('realmd0','sum')
  allocate(sum2(numat),stat=status)
  if (status/=0) call outofmemory('realmd0','sum2')
!
!  Outer loop over sites
!
!  Use Brode-Ahlrichs Algorithm
  noff = numat/2 
  noffset = noff
  if (mod(numat,2_i4).eq.0) then
    noffm1 = noff - 1
  else
    noffm1 = noff
  endif
!orig do i=2,numat
  do i = procid+1,numat,nprocs
    if (i.gt.noff) then
      noffset = noffm1
    endif
    xal = xclat(i)
    yal = yclat(i)
    zal = zclat(i)
    nati = nat(i)
    ntypi = nftype(i)
    nregioni = nregionno(nsft+i)
    nregiontypi = nregiontype(nregioni,ncf)
    neamspeci = neamfnspecptr(i)
    oci = occuf(i)
    lopi = (.not.lfreeze.or.lopf(i))
!
!  Inner loop over second site
!
    ebassum = 0.0_dp
    jloop: do m = 1,noffset
      j = mod(i+m-1_i4,numat) + 1
      lopj = (.not.lfreeze.or.lopf(j))
      if (.not.lopi.and..not.lopj) cycle jloop
      nregionj = nregionno(nsft+j)
      nregiontypj = nregiontype(nregionj,ncf)
      neamspecj = neamfnspecptr(j)
!
!  Find index for i-j
!
      if (neamspeci.gt.neamspecj) then
        indij = neamspeci*(neamspeci - 1)/2 + neamspecj
      else
        indij = neamspecj*(neamspecj - 1)/2 + neamspeci
      endif
!  
!  QM/MM handling : i & j are both QM atoms => exclude
!       
      if (QMMMmode(ncf).gt.0) then
        if (nregiontypi.eq.1.and.nregiontypj.eq.1) cycle jloop
      endif
!
      natj = nat(j)
      ntypj = nftype(j)
      ocj = occuf(j)
      xcrd = xclat(j) - xal
      ycrd = yclat(j) - yal
      zcrd = zclat(j) - zal
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
      cut2 = rp*rp
      r2 = xcrd*xcrd + ycrd*ycrd + zcrd*zcrd
!
      if (r2.lt.smallself.or.r2.gt.cut2) then
        cycle jloop
      else
        r = sqrt(r2)
      endif
!*********************************
!  Calculate unscreened density  *
!*********************************
      call twoden1(1_i4,1_i4,r,xcrd,ycrd,zcrd,npots,npotl,sctrm1,sctrm2,lorder12loc)
!*******************************
!  Compute screening function  *
!*******************************
!
!  Initialise screening function to 1
!
      lnonzeroSij = .true.
      Sij = 1.0_dp
!
      if (lanyMEAMscreen) then
!
!  Set cutoffs for possible screening atoms
!
        rcut2 = rcutfactor*r2
!
!  Loop over atoms to search for images that may contribute to the screening
!
        k = 0
        do while (k.lt.numat.and.lnonzeroSij)
          k = k + 1
          neamspeck = neamfnspecptr(k)
          if (lMEAMscreen(indij,neamspeck)) then
!
!  Set vectors between atoms
!
            xik0 = xclat(k) - xal
            yik0 = yclat(k) - yal
            zik0 = zclat(k) - zal
            xjk0 = xik0 - xcrd
            yjk0 = yik0 - ycrd
            zjk0 = zik0 - zcrd
            xij0 = 0.5_dp*(xik0 + xjk0)
            yij0 = 0.5_dp*(yik0 + yjk0)
            zij0 = 0.5_dp*(zik0 + zjk0)
!
!  Compute square of distance to i-j mid point
!
            r2ijmid = xij0*xij0 + yij0*yij0 + zij0*zij0
            if (r2ijmid.lt.rcut2) then
!
!  Complete distances
!
              r2ik = xik0*xik0 + yik0*yik0 + zik0*zik0
              r2jk = xjk0*xjk0 + yjk0*yjk0 + zjk0*zjk0
!
!  Compute screening function
!
              call meamscreen(neamspeci,neamspecj,neamspeck,r2,r2ik,r2jk,Sikj,dSikjdr,lpartial,.false.)
!
!  If screening function contribution is 0, then no need to continue for this pair
!
              if (Sikj.eq.0.0_dp) then
                lnonzeroSij = .false.
                Sij = 0.0_dp
              else
!
!  Multiply total screening product
!
                Sij = Sij*Sikj
              endif
            endif
          endif
!
!  End loop over possible screening atoms
!
        enddo
      endif
!
      if (lnonzeroSij) then
        if (lMEAM) then
          if (lorder12loc) then
            scrho(1:maxmeamcomponent,i) = scrho(1:maxmeamcomponent,i) + ocj*sctrm1(1:maxmeamcomponent)*Sij
            scrho(1:maxmeamcomponent,j) = scrho(1:maxmeamcomponent,j) + oci*sctrm2(1:maxmeamcomponent)*Sij
          else
            scrho(1:maxmeamcomponent,i) = scrho(1:maxmeamcomponent,i) + ocj*sctrm2(1:maxmeamcomponent)*Sij
            scrho(1:maxmeamcomponent,j) = scrho(1:maxmeamcomponent,j) + oci*sctrm1(1:maxmeamcomponent)*Sij
          endif
        else
          if (lorder12loc) then
            scrho(1,i) = scrho(1,i) + ocj*sctrm1(1)*Sij
            scrho(1,j) = scrho(1,j) + oci*sctrm2(1)*Sij
          else
            scrho(1,i) = scrho(1,i) + ocj*sctrm2(1)*Sij
            scrho(1,j) = scrho(1,j) + oci*sctrm1(1)*Sij
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
            call baskes(nati,ntypi,neamspeci,natj,ntypj,neamspecj,npot,r,rk,Sij,ebas,d1bas,d2bas,.false.,.false.)
            eatom = eatom + ebas
!
            siteenergy(i) = siteenergy(i) + 0.5_dp*ebas
            siteenergy(j) = siteenergy(j) + 0.5_dp*ebas
!
            ebassum = ebassum + ebas
          endif
        endif
      enddo
!
      if (lPrintTwo) then
        write(ioout,'(4x,i12,i12,f22.10)') i,j,ebassum
      endif
    enddo jloop
  enddo
!****************
!  Global sums  *
!****************
  tsum0 = g_cpu_time()
  if (nprocs.gt.1) then
    if (lMEAM) then
      do j = 1,maxmeamcomponent
        do i = 1,numat
          sum2(i) = scrho(j,i)
        enddo
        call sumall(sum2,sum,numat,"realmd0","scrho")
        do i = 1,numat
          scrho(j,i) = sum(i)
        enddo
      enddo
    else
      do i = 1,numat
        sum2(i) = scrho(1,i)
      enddo
      call sumall(sum2,sum,numat,"realmd0","scrho")
      do i = 1,numat
        scrho(1,i) = sum(i)
      enddo
    endif
  endif
  tsuml = g_cpu_time() - tsum0
  tsum = tsum + tsuml
!
!  Free local memory
!
  deallocate(sum2,stat=status)
  if (status/=0) call deallocate_error('realmd0','sum2')
  deallocate(sum,stat=status)
  if (status/=0) call deallocate_error('realmd0','sum')
  deallocate(npotl,stat=status)
  if (status/=0) call deallocate_error('realmd0','npotl')
!
!  Closing banner for energy decomposition
!
  if (lPrintTwo) then
    call mpbarrier
    if (ioproc) then
      write(ioout,'(''--------------------------------------------------------------------------------'',/)')
    endif
  endif
!
!  Timing
!
  time2 = g_cpu_time()
  tatom = tatom + time2 - time1 - tsuml
#ifdef TRACE
  call trace_out('density0')
#endif
!
  return
  end
