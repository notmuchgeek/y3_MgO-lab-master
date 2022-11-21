  subroutine density1fc(matom,vmatom)
!
!  Subroutine for calculating MEAM electron density, including the screening function
!  Version for unphased force constants by finite differences
!
!   1/15 Created from density3
!   4/17 neamspecj now correctly set
!   2/18 Trace added
!  12/20 Corrections to ensure screening terms use all consistent images
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
  use datatypes
  use eam,            only : lMEAM, lMEAMscreen, maxmeamcomponent, meam_Cmax
  use eam,            only : lanyMEAMscreen, neamfnspec
  use element
  use optimisation
  use parallel
  use realvectors
  use sutton
  use symmetry
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  use two
  use vectors,        only : vector_pair
  implicit none
!
!  Passed variables
!
  integer(i4),                   intent(in)    :: matom
  real(dp),                      intent(in)    :: vmatom(4)
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: indij
  integer(i4)                                  :: j
  integer(i4)                                  :: k
  integer(i4)                                  :: kvec
  integer(i4)                                  :: l
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
  integer(i4)                                  :: nor
  integer(i4), dimension(:), allocatable       :: npotl
  integer(i4)                                  :: npots
  integer(i4)                                  :: ntyp1
  integer(i4)                                  :: ntyp2
  integer(i4)                                  :: ntypi
  integer(i4)                                  :: ntypj
  integer(i4)                                  :: nvec
  integer(i4)                                  :: nvec0
  integer(i4)                                  :: status
  logical                                      :: lchange1
  logical                                      :: lchange2
  logical                                      :: lcspair
  logical,                                save :: lfirstcall = .true.
  logical                                      :: lmatch
  logical                                      :: lnonzeroSii
  logical                                      :: lnonzeroSij
  logical                                      :: lorder12loc
  logical                                      :: lpartial
  logical                                      :: lself
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: cut2
  real(dp)                                     :: oci
  real(dp)                                     :: ocj
  real(dp)                                     :: rcut2
  real(dp)                                     :: rcutfactor
  real(dp)                                     :: rp
  real(dp)                                     :: dSikjdr(3)       ! Dummy argument for call to meamscreen - not used here
  real(dp)                                     :: Sii
  real(dp)                                     :: Sij
  real(dp)                                     :: Siki
  real(dp)                                     :: Sikj
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
  call trace_in('density1fc')
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
!  Allocate local memory
!
  allocate(npotl(npote),stat=status)
  if (status/=0) call outofmemory('density1fc','npotl')
!**************************************************************
!  Calculation of densities and screening functions for MEAM  *
!**************************************************************
!
!  Outer loop over sites
!
  do i = 1,numat
    xal = xclat(i)
    yal = yclat(i)
    zal = zclat(i)
    nati = nat(i)
    ntypi = nftype(i)
    oci = occuf(i)
    neamspeci = neamfnspecptr(i)
!
!  Start of second atom loop
!
    jloop: do j = 1,i-1
      natj = nat(j)
      ntypj = nftype(j)
      neamspecj = neamfnspecptr(i)
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
!
!  Find index for i-j
!
      if (neamspeci.gt.neamspecj) then
        indij = neamspeci*(neamspeci - 1)/2 + neamspecj
      else
        indij = neamspecj*(neamspecj - 1)/2 + neamspeci
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
!  Modify distances if matom equals i or j
!
      if (i.eq.matom) then
        do l = 1,nor
          xtmp(l) = xtmp(l) - (vmatom(1) - xal)
          ytmp(l) = ytmp(l) - (vmatom(2) - yal)
          ztmp(l) = ztmp(l) - (vmatom(3) - zal)
          dist(l) = xtmp(l)**2 + ytmp(l)**2 + ztmp(l)**2
        enddo
      endif
      if (j.eq.matom) then
        do l = 1,nor
          xtmp(l) = xtmp(l) + (vmatom(1) - xclat(j))
          ytmp(l) = ytmp(l) + (vmatom(2) - yclat(j))
          ztmp(l) = ztmp(l) + (vmatom(3) - zclat(j))
          dist(l) = xtmp(l)**2 + ytmp(l)**2 + ztmp(l)**2
        enddo
      endif
!
!  Sqrt distances
!
      do l = 1,nor
        dist(l) = sqrt(dist(l))
      enddo
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
!
              xjk0 = xik0 - xtmp(l)
              yjk0 = yik0 - ytmp(l)
              zjk0 = zik0 - ztmp(l)
!
!  Reverse displacements to xtmp in j-k vector so that cutoffs are consistent for rfindmid call
!
              if (i.eq.matom) then
                xjk0 = xjk0 - (vmatom(1) - xal)
                yjk0 = yjk0 - (vmatom(2) - yal)
                zjk0 = zjk0 - (vmatom(3) - zal)
              endif
              if (j.eq.matom) then
                xjk0 = xjk0 + (vmatom(1) - xclat(j))
                yjk0 = yjk0 + (vmatom(2) - yclat(j))
                zjk0 = zjk0 + (vmatom(3) - zclat(j))
              endif
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
!  Correct distances for matom
!
                lchange1 = .false.
                lchange2 = .false.
                if (i.eq.matom) then
                  vectorpair%xvector_pair1(kvec) = vectorpair%xvector_pair1(kvec) - (vmatom(1) - xclat(i))
                  vectorpair%yvector_pair1(kvec) = vectorpair%yvector_pair1(kvec) - (vmatom(2) - yclat(i))
                  vectorpair%zvector_pair1(kvec) = vectorpair%zvector_pair1(kvec) - (vmatom(3) - zclat(i))
                  lchange1 = .true.
                endif
                if (j.eq.matom) then
                  vectorpair%xvector_pair2(kvec) = vectorpair%xvector_pair2(kvec) - (vmatom(1) - xclat(j))
                  vectorpair%yvector_pair2(kvec) = vectorpair%yvector_pair2(kvec) - (vmatom(2) - yclat(j))
                  vectorpair%zvector_pair2(kvec) = vectorpair%zvector_pair2(kvec) - (vmatom(3) - zclat(j))
                  lchange2 = .true.
                endif
                if (k.eq.matom) then
                  vectorpair%xvector_pair1(kvec) = vectorpair%xvector_pair1(kvec) + (vmatom(1) - xclat(k))
                  vectorpair%yvector_pair1(kvec) = vectorpair%yvector_pair1(kvec) + (vmatom(2) - yclat(k))
                  vectorpair%zvector_pair1(kvec) = vectorpair%zvector_pair1(kvec) + (vmatom(3) - zclat(k))
                  vectorpair%xvector_pair2(kvec) = vectorpair%xvector_pair2(kvec) + (vmatom(1) - xclat(k))
                  vectorpair%yvector_pair2(kvec) = vectorpair%yvector_pair2(kvec) + (vmatom(2) - yclat(k))
                  vectorpair%zvector_pair2(kvec) = vectorpair%zvector_pair2(kvec) + (vmatom(3) - zclat(k))
                  lchange1 = .true.
                  lchange2 = .true.
                endif
!
!  If distances have changed then recompute
!
                if (lchange1) then
                  vectorpair%distance_pair1(kvec) = vectorpair%xvector_pair1(kvec)**2 + &
                                                    vectorpair%yvector_pair1(kvec)**2 + &
                                                    vectorpair%zvector_pair1(kvec)**2
                endif
                if (lchange2) then
                  vectorpair%distance_pair2(kvec) = vectorpair%xvector_pair2(kvec)**2 + &
                                                    vectorpair%yvector_pair2(kvec)**2 + &
                                                    vectorpair%zvector_pair2(kvec)**2
                endif
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
            else
              scrho(1:maxmeamcomponent,i) = scrho(1:maxmeamcomponent,i) + sctrm2(1:maxmeamcomponent)*Sij*ocj
              scrho(1:maxmeamcomponent,j) = scrho(1:maxmeamcomponent,j) + sctrm1(1:maxmeamcomponent)*Sij*oci
            endif
          else
            if (lorder12loc) then
              scrho(1,i) = scrho(1,i) + sctrm1(1)*Sij*ocj
              scrho(1,j) = scrho(1,j) + sctrm2(1)*Sij*oci
            else
              scrho(1,i) = scrho(1,i) + sctrm2(1)*Sij*ocj
              scrho(1,j) = scrho(1,j) + sctrm1(1)*Sij*oci
            endif
          endif
        endif
      enddo norloop
    enddo jloop
  enddo
!*******************
!  Self-term loop  *
!*******************
  iloop: do i = 1,numat
    xal = xclat(i)
    yal = yclat(i)
    zal = zclat(i)
    nati = nat(i)
    ntypi = nftype(i)
    oci = occuf(i)
!
!  Find index for i-j
!
    indij = neamspeci*(neamspeci + 1)/2
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
    enddo norloopself
!
!  End of self term loop
!
  enddo iloop
!
!  Free local memory
!
  deallocate(npotl,stat=status)
  if (status/=0) call deallocate_error('density1fc','npotl')
!
!  Timing
!
  time2 = g_cpu_time()
  tatom = tatom + time2 - time1
#ifdef TRACE
  call trace_out('density1fc')
#endif
!
  return
  end subroutine density1fc
