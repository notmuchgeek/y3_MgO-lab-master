  subroutine density3(eatom,esregion12,esregion2,eattach)
!
!  Subroutine for calculating MEAM electron density, including the screening function
!
!  lfreeze = if .true. only do derivatives for atoms with a degree of freedom
!
!   4/09 Created from realmd3 as a template
!   4/09 Screening function added
!   1/14 vectorpair initialised and saved to avoid memory leak
!   8/14 MEAM screening made species specific
!   8/14 Pair potential for MEAM now evaluated here so that screening function
!        can be used. 
!   8/14 Energy passed as argument so that pair potential contribution can be added
!   8/14 All MEAM species (i,j,k) now passed to meanscreen
!   8/14 neamspec passed psibaskes
!   1/15 Name of subroutine corrected in memory allocation calls
!   2/18 Trace added
!   5/19 Site energy added
!   5/19 Twobody energy print out added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
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
!  Julian Gale, CIC, Curtin University, September 2019
!
  use configurations, only : lsliceatom, nregionno, nregiontype, QMMMmode, nregions
  use control
  use current
  use datatypes
  use eam,            only : lMEAM, lMEAMscreen, maxmeamcomponent, meam_Cmax
  use eam,            only : lanyMEAMscreen, neamfnspec
  use element
  use energies,       only : siteenergy, eregion2region
  use iochannels,     only : ioout
  use optimisation
  use parallel
  use realvectors
  use sutton
  use symmetry
  use times
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  use two
  use vectors,        only : vector_pair
  implicit none
!
!  Passed variables
!
  real(dp),                     intent(inout)  :: eatom
  real(dp),                     intent(inout)  :: eattach
  real(dp),                     intent(inout)  :: esregion12
  real(dp),                     intent(inout)  :: esregion2
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: indij
  integer(i4)                                  :: j
  integer(i4)                                  :: k
  integer(i4)                                  :: kvec
  integer(i4)                                  :: l
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
  integer(i4)                                  :: nmolonly
  integer(i4)                                  :: noff
  integer(i4)                                  :: noffm1
  integer(i4)                                  :: noffset
  integer(i4)                                  :: nor
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
  integer(i4)                                  :: nvec
  integer(i4)                                  :: nvec0
  integer(i4)                                  :: status
  logical                                      :: lattach
  logical                                      :: lcspair
  logical,                                save :: lfirstcall = .true.
  logical                                      :: lmatch
  logical                                      :: lnonzeroSii
  logical                                      :: lnonzeroSij
  logical                                      :: lopi
  logical                                      :: lopj
  logical                                      :: lorder12loc
  logical                                      :: lpartial
  logical                                      :: lreg2one
  logical                                      :: lreg2pair
  logical                                      :: lself
  logical                                      :: lslicei
  logical                                      :: lslicej
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: cut2
  real(dp)                                     :: d1bas
  real(dp)                                     :: d2bas
  real(dp)                                     :: ebas
  real(dp)                                     :: ebassum
  real(dp)                                     :: oci
  real(dp)                                     :: ocj
  real(dp)                                     :: rcut2
  real(dp)                                     :: rcutfactor
  real(dp)                                     :: r
  real(dp)                                     :: rk
  real(dp)                                     :: rp
  real(dp)                                     :: dSikjdr(3)       ! Dummy argument for call to meamscreen - not used here
  real(dp)                                     :: Sii
  real(dp)                                     :: Sij
  real(dp)                                     :: Siki
  real(dp)                                     :: Sikj
  real(dp)                                     :: sctrm1(maxmeamcomponent)
  real(dp)                                     :: sctrm2(maxmeamcomponent)
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
  real(dp)                                     :: xik0
  real(dp)                                     :: yik0
  real(dp)                                     :: zik0
  real(dp)                                     :: xjk0
  real(dp)                                     :: yjk0
  real(dp)                                     :: zjk0
  type(vector_pair),                      save :: vectorpair
!
!  Check that this call is needed
!
  if (.not.lsuttonc) return
#ifdef TRACE
  call trace_in('density3')
#endif
!
!  If this is the first call then initialise vectorpair
!
  if (lfirstcall) then
    lfirstcall = .false.
    call changemaxvectorpair(vectorpair,0_i4)
  endif
!
  time1 = g_cpu_time()
!
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
  if (status/=0) call outofmemory('density3','npotl')
  allocate(sum(numat),stat=status)
  if (status/=0) call outofmemory('density3','sum')
  allocate(sum2(numat),stat=status)
  if (status/=0) call outofmemory('density3','sum2')
!**************************************************************
!  Calculation of densities and screening functions for MEAM  *
!**************************************************************
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
!
  do i = procid+1,numat,nprocs
    if (i.gt.noff) then
      noffset = noffm1
    endif
    xal = xclat(i)
    yal = yclat(i)
    zal = zclat(i)
    nati = nat(i)
    ntypi = nftype(i)
    oci = occuf(i)
    nregioni = nregionno(nsft+nrelf2a(i))
    nregiontypi = nregiontype(nregioni,ncf)
    lopi = (.not.lfreeze.or.lopf(nrelf2a(i)))
    lslicei = lsliceatom(nsft + nrelf2a(i))
    neamspeci = neamfnspecptr(i)
!
!  Start of second atom loop
!
    jloop: do m = 1,noffset
      j = mod(i+m-1_i4,numat) + 1
      natj = nat(j)
      ntypj = nftype(j)
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
        nat2 = nat(j)
        ntyp1 = ntypi
        ntyp2 = nftype(j)
      else
        lorder12loc = .false.
        nat1 = nat(j)
        nat2 = nati
        ntyp1 = nftype(j)
        ntyp2 = ntypi
      endif
      nregionj = nregionno(nsft+nrelf2a(j))
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
!  Freezing flag
!
      lopj = (.not.lfreeze.or.lopf(nrelf2a(j)))
      if (.not.lopi.and..not.lopj) cycle jloop
!
!  Set region 2 pair flag
!
      lreg2one  = .false.
      lreg2pair = .false.
      if (lseok.and.nregions(ncf).ge.2) then
        lreg2pair = (nregioni.eq.2.and.nregionj.eq.2)
        if (.not.lreg2pair) lreg2one = (nregioni.eq.2.or.nregionj.eq.2)
      endif
      lslicej = lsliceatom(nsft + nrelf2a(j))
      lattach = (lslicei.and..not.lslicej.or.lslicej.and..not.lslicei)
!  
!  QM/MM handling : i & j are both QM atoms => exclude
!
      if (QMMMmode(ncf).gt.0) then
        if (nregiontypi.eq.1.and.nregiontypj.eq.1) cycle jloop
      endif
!
      xcrd = xclat(j) - xal
      ycrd = yclat(j) - yal
      zcrd = zclat(j) - zal
      ocj = occuf(j)
!
!  Possible core-shell flag
!
      lcspair = (abs(nat1-nat2).eq.maxele.or.(oci+ocj).lt.1.0001_dp)
!
!  Locate potential number
!  Check whether potential requires specific types
!  Calculate sum of all dispersion terms for pair
!
      rp = 0.0_dp
      npots = 0
      do n = 1,npote
!
!  Screen for density terms and pair potentials
!
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
!  If no valid potentials and charge product is zero
!  then no need to search for distances
!
      if (npots.eq.0) cycle jloop
!
      cut2 = rp*rp
!***********************
!  Find valid vectors  *
!***********************
      if (ndim.eq.3) then
        call rsearch3D(xcrd,ycrd,zcrd,.false.,lcspair,i,j,0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2)
      elseif (ndim.eq.2) then
        call rsearch2D(xcrd,ycrd,zcrd,.false.,lcspair,i,j,0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2)
      elseif (ndim.eq.1) then
        call rsearch1D(xcrd,ycrd,zcrd,.false.,lcspair,i,j,0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2)
      endif
!
!  If there are no valid vectors then cycle
!
      if (nor.eq.0) cycle jloop
!
!  Sqrt distances
!
      do l = 1,nor
        dist(l) = sqrt(dist(l))
      enddo
      ebassum = 0.0_dp
!*******************************
!  Loop over vectors from i-j  *
!*******************************
      norloop: do l = 1,nor
!*******************************
!  Compute unscreened density  *
!*******************************
        call twoden(l,l,npots,npotl,sctrm1,sctrm2,lorder12loc)
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
          rcut2 = rcutfactor*(dist(l))**2
!
!  Loop over atoms to search for images that may contribute to the screening
!
          k = 0
          do while (k.lt.numat.and.lnonzeroSij)
            k = k + 1
            neamspeck = neamfnspecptr(k)
            if (lMEAMscreen(indij,neamspeck)) then
!
!  Set basic vectors between atoms
!
              xik0 = xclat(k) - xal
              yik0 = yclat(k) - yal
              zik0 = zclat(k) - zal
              xjk0 = xik0 - xtmp(l)
              yjk0 = yik0 - ytmp(l)
              zjk0 = zik0 - ztmp(l)
!
!  Find images within cutoffs of both atoms - excluding self images
!
              nvec0 = 0
              call rfindmid(xik0,yik0,zik0,xjk0,yjk0,zjk0,rcut2,.false.,nvec0,nvec,vectorpair)
!
!  Loop over results of search
!
              kvec = 0
              do while (kvec.lt.nvec.and.lnonzeroSij)
                kvec = kvec + 1
!
!  Compute screening function
!
                call meamscreen(neamspeci,neamspecj,neamspeck,dist(l)**2,vectorpair%distance_pair1(kvec), &
                                vectorpair%distance_pair2(kvec),Sikj,dSikjdr,lpartial,.false.)
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
!
!  End loop over images of possible screening atoms
!
              enddo
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
              scrho(1:maxmeamcomponent,i) = scrho(1:maxmeamcomponent,i) + sctrm1(1:maxmeamcomponent)*Sij*ocj
              scrho(1:maxmeamcomponent,j) = scrho(1:maxmeamcomponent,j) + sctrm2(1:maxmeamcomponent)*Sij*oci
              if (lattach) then
                scrho12(1:maxmeamcomponent,i) = scrho12(1:maxmeamcomponent,i) + sctrm1(1:maxmeamcomponent)*Sij*ocj
                scrho12(1:maxmeamcomponent,j) = scrho12(1:maxmeamcomponent,j) + sctrm2(1:maxmeamcomponent)*Sij*oci
              endif
            else
              scrho(1:maxmeamcomponent,i) = scrho(1:maxmeamcomponent,i) + sctrm2(1:maxmeamcomponent)*Sij*ocj
              scrho(1:maxmeamcomponent,j) = scrho(1:maxmeamcomponent,j) + sctrm1(1:maxmeamcomponent)*Sij*oci
              if (lattach) then
                scrho12(1:maxmeamcomponent,i) = scrho12(1:maxmeamcomponent,i) + sctrm2(1:maxmeamcomponent)*Sij*ocj
                scrho12(1:maxmeamcomponent,j) = scrho12(1:maxmeamcomponent,j) + sctrm1(1:maxmeamcomponent)*Sij*oci
              endif
            endif
          else
            if (lorder12loc) then
              scrho(1,i) = scrho(1,i) + sctrm1(1)*Sij*ocj
              scrho(1,j) = scrho(1,j) + sctrm2(1)*Sij*oci
              if (lattach) then
                scrho12(1,i) = scrho12(1,i) + sctrm1(1)*Sij*ocj
                scrho12(1,j) = scrho12(1,j) + sctrm2(1)*Sij*oci
              endif
            else
              scrho(1,i) = scrho(1,i) + sctrm2(1)*Sij*ocj
              scrho(1,j) = scrho(1,j) + sctrm1(1)*Sij*oci
              if (lattach) then
                scrho12(1,i) = scrho12(1,i) + sctrm2(1)*Sij*ocj
                scrho12(1,j) = scrho12(1,j) + sctrm1(1)*Sij*oci
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
            r = dist(l)
            if (r.gt.rpot2(npot).and.r.le.rpot(npot)) then
              rk = 1.0_dp/r
              call baskes(nati,ntypi,neamspeci,natj,ntypj,neamspecj,npot,r,rk,1.0_dp,ebas,d1bas,d2bas,.false.,.false.)
!
              if (lreg2one) then
                esregion12 = esregion12 + ebas*Sij
              elseif (lreg2pair) then
                esregion2 = esregion2 + ebas*Sij
              else
                eatom = eatom + ebas*Sij
              endif
              if (lattach) eattach = eattach + ebas*Sij
!
              eregion2region(nregionj,nregioni) = eregion2region(nregionj,nregioni) + ebas*Sij
!
              siteenergy(i) = siteenergy(i) + 0.5_dp*ebas*Sij
              siteenergy(j) = siteenergy(j) + 0.5_dp*ebas*Sij
!
              ebassum = ebassum + ebas*Sij
            endif
          endif
        enddo
      enddo norloop
!
      if (lPrintTwo) then
        write(ioout,'(4x,i12,i12,f22.10)') i,j,ebassum
      endif
    enddo jloop
  enddo
!*******************
!  Self-term loop  *
!*******************
  iloop: do i = procid+1,numat,nprocs
    lopi = (.not.lfreeze.or.lopf(nrelf2a(i)))
    if (.not.lopi) cycle iloop
    xal = xclat(i)
    yal = yclat(i)
    zal = zclat(i)
    nati = nat(i)
    ntypi = nftype(i)
    oci = occuf(i)
    nregioni = nregionno(nsft+nrelf2a(i))
    nregiontypi = nregiontype(nregioni,ncf)
    neamspeci = neamfnspecptr(i)
!
!  Find index for i-j
!
    indij = neamspeci*(neamspeci + 1)/2
!
!  Set region 2 pair flag
!
    lreg2pair = .false.
    if (lseok.and.nregions(ncf).ge.2) then
      lreg2pair = (nregionno(nsft+nrelf2a(i)).eq.2)
    endif
!
!  QM/MM handling : i is a QM atom => exclude
!
    if (QMMMmode(ncf).gt.0) then
      if (nregiontypi.eq.1) cycle iloop
    endif
!
!  Locate potential number
!  Check whether potential requires specific types
!  Calculate sum of dispersion terms
!
    rp = 0.0_dp
    npots = 0
    do n = 1,npote
      if (nptype(n).eq.19.or.nptype(n).eq.45.or.nptype(n).eq.55) then
        if (lmatch(nati,ntypi,nspec1(n),nptyp1(n),.true.)) then
          if (lmatch(nati,ntypi,nspec2(n),nptyp2(n),.true.)) then
            npots = npots + 1
            npotl(npots) = n
            if (rpot(n).gt.rp) rp = rpot(n)
          endif
        endif
      endif
    enddo
    cut2 = rp*rp
!***********************
!  Find valid vectors  *
!***********************
    if (ndim.eq.3) then
      call rsearch3D(0.0_dp,0.0_dp,0.0_dp,.false.,lcspair,i,i,0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2)
    elseif (ndim.eq.2) then
      call rsearch2D(0.0_dp,0.0_dp,0.0_dp,.false.,lcspair,i,i,0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2)
    elseif (ndim.eq.1) then
      call rsearch1D(0.0_dp,0.0_dp,0.0_dp,.false.,lcspair,i,i,0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2)
    endif
!
!  If there are no vectors then cycle
!
    if (nor.eq.0) cycle iloop
!
!  Sqrt distances
!
    do l = 1,nor
      dist(l) = sqrt(dist(l))
    enddo
    ebassum = 0.0_dp
!*************************************************
!  Loop over vectors from i to images of itself  *
!*************************************************
    norloopself: do l = 1,nor
!*******************************
!  Compute unscreened density  *
!*******************************
      call twoden(l,l,npots,npotl,sctrm1,sctrm2,.true.)
!*******************************
!  Compute screening function  *
!*******************************
!
!  Initialise screening function to 1
!
      lnonzeroSii = .true.
      Sii = 1.0_dp
!
      if (lanyMEAMscreen) then
!
!  Set cutoffs for possible screening atoms
!
        rcut2 = rcutfactor*(dist(l))**2
!
!  Loop over atoms to search for images that may contribute to the screening
!
        k = 0
        do while (k.lt.numat.and.lnonzeroSii)
          k = k + 1
          neamspeck = neamfnspecptr(k)
          if (lMEAMscreen(indij,neamspeck)) then
!
!  Set basic vectors between atoms
!
            xik0 = xclat(k) - xal
            yik0 = yclat(k) - yal
            zik0 = zclat(k) - zal
            xjk0 = xik0 - xtmp(l)
            yjk0 = yik0 - ytmp(l)
            zjk0 = zik0 - ztmp(l)
!
!  Find images within cutoffs of both atoms - excluding self images
!
            nvec0 = 0
            call rfindmid(xik0,yik0,zik0,xjk0,yjk0,zjk0,rcut2,.false.,nvec0,nvec,vectorpair)
!
!  Loop over results of search
!
            kvec = 0
            do while (kvec.lt.nvec.and.lnonzeroSii)
              kvec = kvec + 1
!
!  Compute screening function
!
              call meamscreen(neamspeci,neamspeci,neamspeck,dist(l)**2,vectorpair%distance_pair1(kvec), &
                              vectorpair%distance_pair2(kvec),Siki,dSikjdr,lpartial,.false.)
!
!  If screening function contribution is 0, then no need to continue for this pair
!
              if (Siki.eq.0.0_dp) then
                lnonzeroSii = .false.
                Sii = 0.0_dp
              else
!
!  Multiply total screening product
!
                Sii = Sii*Siki
              endif
!
!  End loop over images of possible screening atoms
!
            enddo
          endif
!
!  End loop over possible screening atoms
!
        enddo
      endif
!
      if (lnonzeroSii) then
        if (lMEAM) then
          scrho(1:maxmeamcomponent,i) = scrho(1:maxmeamcomponent,i) + sctrm1(1:maxmeamcomponent)*Sii*oci
        else
          scrho(1,i) = scrho(1,i) + sctrm1(1)*Sii*oci
        endif
      endif
!
!  Compute pairwise contribution to energy
!
      do np = 1,npots
        npot = npotl(np)
        if (nptype(npot).eq.45.or.nptype(npot).eq.55) then
          r = dist(l)
          if (r.gt.rpot2(npot).and.r.le.rpot(npot)) then
            rk = 1.0_dp/r
            call baskes(nati,ntypi,neamspeci,nati,ntypi,neamspeci,npot,r,rk,1.0_dp,ebas,d1bas,d2bas,.false.,.false.)
            if (lreg2pair) then
              esregion2 = esregion2 + 0.5_dp*ebas*Sii
            else
              eatom = eatom + 0.5_dp*ebas*Sii
            endif
!
            eregion2region(nregioni,nregioni) = eregion2region(nregioni,nregioni) + 0.5_dp*ebas*Sii
!
            siteenergy(i) = siteenergy(i) + 0.5_dp*ebas*Sii
!
            ebassum = ebassum + 0.5_dp*ebas*Sii
          endif
        endif
      enddo
!
    enddo norloopself
!
    if (lPrintTwo) then
      write(ioout,'(4x,i12,i12,f22.10)') i,i,ebassum
    endif
!
!  End of self term loop
!
  enddo iloop
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
        call sumall(sum2,sum,numat,"density3","scrho")
        do i = 1,numat
          scrho(j,i) = sum(i)
        enddo
        do i = 1,numat
          sum2(i) = scrho12(j,i)
        enddo
        call sumall(sum2,sum,numat,"density3","scrho12")
        do i = 1,numat
          scrho12(j,i) = sum(i)
        enddo
      enddo
    else
      do i = 1,numat
        sum2(i) = scrho(1,i)
      enddo
      call sumall(sum2,sum,numat,"density3","scrho")
      do i = 1,numat
        scrho(1,i) = sum(i)
      enddo
      do i = 1,numat
        sum2(i) = scrho12(1,i)
      enddo
      call sumall(sum2,sum,numat,"density3","scrho12")
      do i = 1,numat
        scrho12(1,i) = sum(i)
      enddo
    endif
  endif
  tsuml = g_cpu_time() - tsum0
  tsum = tsum + tsuml
!
!  Free local memory
!
  deallocate(sum2,stat=status)
  if (status/=0) call deallocate_error('density3','sum2')
  deallocate(sum,stat=status)
  if (status/=0) call deallocate_error('density3','sum')
  deallocate(npotl,stat=status)
  if (status/=0) call deallocate_error('density3','npotl')
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
  call trace_out('density3')
#endif
!
  return
  end subroutine density3
