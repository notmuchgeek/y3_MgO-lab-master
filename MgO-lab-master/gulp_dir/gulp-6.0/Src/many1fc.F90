  subroutine many1fc(maxlhs,d1cell,matom,vmatom)
! 
!  Subroutine for calculating the many-body phonons from the
!  Sutton-Chen potential(s). Requires real space routine to
!  be called first to evaluate the repulsive contribution
!  and the rho parameters. 
!  Unphased finite difference version for fc_supercell option. 
!
!  On entry the array scrho must contain the density at each atomic site.
!
!   1/15 Created from manyfc / manymd3
!   2/15 Debugged for EAM
!   2/15 Cycling of loops for zero density removed for lMEAMden case
!   2/18 Trace added
!  11/18 Cutoff check added for baskes
!  12/20 matom now passed to psiscreenderv instead of i
!  12/20 Incorrect checking of i.ne.j removed for derivative adding
!  12/20 Modified to capture missing terms when matom is k
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
!  Julian Gale, CIC, Curtin University, December 2020
!
  use control
  use current
  use derivatives
  use eam
  use general
  use realvectors,    only : dist, xtmp, ytmp, ztmp
  use realvectors,    only : cellindex
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
  integer(i4),                     intent(in)    :: maxlhs
  integer(i4),                     intent(in)    :: matom
  real(dp),                        intent(out)   :: d1cell(4,maxlhs,*)
  real(dp),                        intent(in)    :: vmatom(4)
!
!  Local variables
!
  integer(i4)                                    :: i
  integer(i4)                                    :: ii
  integer(i4)                                    :: indij
  integer(i4)                                    :: j
  integer(i4)                                    :: k
  integer(i4)                                    :: kvec
  integer(i4)                                    :: m
  integer(i4)                                    :: n  
  integer(i4)                                    :: nat1
  integer(i4)                                    :: nat2
  integer(i4)                                    :: nati
  integer(i4)                                    :: natj
  integer(i4)                                    :: ncindc
  integer(i4)                                    :: ncindm
  integer(i4)                                    :: ncindp
  integer(i4)                                    :: neamspeci
  integer(i4)                                    :: neamspecj
  integer(i4)                                    :: neamspeck
  integer(i4)                                    :: neamfnspec2
  integer(i4)                                    :: nmolonly
  integer(i4)                                    :: nor
  integer(i4)                                    :: np
  integer(i4)                                    :: npartial
  integer(i4)                                    :: npot
  integer(i4)                                    :: npots
  integer(i4)                                    :: ntyp1  
  integer(i4)                                    :: ntyp2
  integer(i4)                                    :: ntypi
  integer(i4)                                    :: ntypj    
  integer(i4), dimension(:),   allocatable       :: npotl
  integer(i4)                                    :: nvec
  integer(i4)                                    :: nvec0
  integer(i4)                                    :: status
  logical,                                  save :: lfirstcall = .true.
  logical                                        :: lchange1
  logical                                        :: lchange2
  logical                                        :: lnonzeroSij
  logical                                        :: lpartial
  logical                                        :: lself 
  logical                                        :: lvalidij
  real(dp)                                       :: g_cpu_time
  real(dp)                                       :: cut2
  real(dp)                                       :: cut2p
  real(dp)                                       :: cut2r
  real(dp)                                       :: deriv(3)
  real(dp)                                       :: drhoij(3,maxmeamcomponent)
  real(dp)                                       :: drhoijs(6,maxmeamcomponent)
  real(dp)                                       :: drhoij2(6,maxmeamcomponent)
  real(dp)                                       :: drhoij2s(21,maxmeamcomponent)
  real(dp)                                       :: drhoij2m(6,3,maxmeamcomponent)
  real(dp)                                       :: drhoji(3,maxmeamcomponent)
  real(dp)                                       :: drhojis(6,maxmeamcomponent)
  real(dp)                                       :: drhoji2(6,maxmeamcomponent)
  real(dp)                                       :: drhoji2s(21,maxmeamcomponent)
  real(dp)                                       :: drhoji2m(6,3,maxmeamcomponent)
  real(dp)                                       :: drhototij(3)
  real(dp)                                       :: drhototijs(6)
  real(dp)                                       :: drhototij2(6)
  real(dp)                                       :: drhototij2s(21)
  real(dp)                                       :: drhototij2m(6,3)
  real(dp)                                       :: drhototij3(10)
  real(dp)                                       :: drhototji(3)
  real(dp)                                       :: drhototjis(6)
  real(dp)                                       :: drhototji2(6)
  real(dp)                                       :: drhototji2s(21)
  real(dp)                                       :: drhototji2m(6,3)
  real(dp)                                       :: drhototji3(10)
  real(dp)                                       :: ebas
  real(dp)                                       :: d1bas
  real(dp)                                       :: d2bas
  real(dp)                                       :: eeam
  real(dp)                                       :: frho(maxmeamcomponent)
  real(dp)                                       :: oci
  real(dp)                                       :: ocj
  real(dp)                                       :: ofct
  real(dp)                                       :: r
  real(dp)                                       :: rk
  real(dp)                                       :: r2
  real(dp)                                       :: rcut2
  real(dp)                                       :: rcutfactor
  real(dp)                                       :: rhoi
  real(dp)                                       :: rhoj
  real(dp)                                       :: rhoij(maxmeamcomponent)
  real(dp)                                       :: rhoji(maxmeamcomponent)
  real(dp)                                       :: rp
  real(dp)                                       :: rscrhoi
  real(dp)                                       :: rscrhoi3
  real(dp)                                       :: rscrhoi5
  real(dp)                                       :: rscrhoj
  real(dp)                                       :: rscrhoj3
  real(dp)                                       :: rscrhoj5
  real(dp)                                       :: scmax
  real(dp)                                       :: Sij
  real(dp)                                       :: Sikj
  real(dp)                                       :: dSikjdr(3)
  real(dp)                                       :: time1
  real(dp)                                       :: time2  
  real(dp)                                       :: xal 
  real(dp)                                       :: yal    
  real(dp)                                       :: zal
  real(dp)                                       :: xcd
  real(dp)                                       :: ycd   
  real(dp)                                       :: zcd
  real(dp)                                       :: xcrd 
  real(dp)                                       :: ycrd    
  real(dp)                                       :: zcrd
  real(dp)                                       :: xik0
  real(dp)                                       :: yik0
  real(dp)                                       :: zik0
  real(dp)                                       :: xjk0
  real(dp)                                       :: yjk0
  real(dp)                                       :: zjk0
  type(screening_atoms)                          :: partial
  type(vector_pair),                        save :: vectorpair
#ifdef TRACE
  call trace_in('many1fc')
#endif
!
  time1 = g_cpu_time()
!***************************
!  Set up local variables  *
!***************************
!
!  Find maximum cut-off radius
!
  scmax = 0.0_dp
  do i = 1,npote
    if (nptype(i).eq.19.or.nptype(i).eq.45.or.nptype(i).eq.55) then
      if (rpot(i).gt.scmax) scmax = rpot(i)
    endif
  enddo
!
!  Set up cutoffs
!
  cut2p = cutp*cutp
  if (lnoreal) goto 1000
!
!  If this is the first call then initialise vectorpair
!
  if (lfirstcall) then
    lfirstcall = .false.
    call changemaxvectorpair(vectorpair,0_i4)
  endif
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
!  Scale density
!
  call eamscalescrho(1_i4)
!
!  Allocate local memory
!
  allocate(npotl(npote),stat=status)
  if (status/=0) call outofmemory('many1fc','npotl')
!
!  Outer loop over i
!
  iloop: do i = 1,numat
!
!  Find EAM species for i
!  
    neamspeci = neamfnspecptr(i)
!
!  Evaluate functional derivatives
!
    if (lMEAMfn) then
      frho(1:maxmeamcomponent) = scrho(1:maxmeamcomponent,i)
      call meamfnderv(neamfn,neamspeci,frho,rhoi,eeam,rscrhoi,rscrhoi3,rscrhoi5,.true.,.false.,.false.)
    else
      call eamfnderv(neamfn,neamspeci,scrho(1,i),eeam,rscrhoi,rscrhoi3,rscrhoi5,.true.,.false.,.false.)
      rhoi = scrho(1,i)
    endif
!
    xal = xclat(i)
    yal = yclat(i)
    zal = zclat(i)
    nati = nat(i)
    ntypi = nftype(i)
    oci = occuf(i)
!
!  Start of second atom loop 
!
    jloop: do j = 1,i
!     
!  Find EAM species for j
!  
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
!  Evaluate functional derivatives
!
      if (lMEAMfn) then
        frho(1:maxmeamcomponent) = scrho(1:maxmeamcomponent,j)
        call meamfnderv(neamfn,neamspecj,frho,rhoj,eeam,rscrhoj,rscrhoj3,rscrhoj5,.true.,.false.,.false.)
      else
        call eamfnderv(neamfn,neamspecj,scrho(1,j),eeam,rscrhoj,rscrhoj3,rscrhoj5,.true.,.false.,.false.)
        rhoj = scrho(1,j)
      endif
!
!  Skip if this is not a MEAM density and the densities are zero
!
      if (.not.lMEAMden.and.rhoi.eq.0.0_dp.and.rhoj.eq.0.0_dp) cycle jloop
!
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
          ntyp1 = ntypi
          ntyp2 = ntypj
        else
          ntyp1 = ntypj
          ntyp2 = ntypi
        endif
      elseif (nati.lt.natj) then
        nat1 = nati
        nat2 = nat(j)
        ntyp1 = ntypi
        ntyp2 = nftype(j)
      else
        nat1 = nat(j)
        nat2 = nati
        ntyp1 = nftype(j)
        ntyp2 = ntypi
      endif
      xcrd = xclat(j) - xal
      ycrd = yclat(j) - yal
      zcrd = zclat(j) - zal
      ocj = occuf(j)
!
      if (i.eq.j) then
        ofct = 0.5_dp*oci*ocj
      else
        ofct = oci*ocj
      endif
!
!  Locate potential number
!  Check whether potential requires specific types
!  Calculate sum of all dispersion terms for pair
!
      npots = 0
      rp = 0.0_dp
      do n = 1,npote
        if (nptype(n).eq.19.or.nptype(n).eq.45.or.nptype(n).eq.55) then
          if (nat1.eq.nspec1(n).and.nat2.eq.nspec2(n)) then
            if ((ntyp1.eq.nptyp1(n).or.nptyp1(n).eq.0).and.(ntyp2.eq.nptyp2(n).or.nptyp2(n).eq.0)) then
              npots = npots + 1
              npotl(npots) = n
              if (rpot(n).gt.rp) rp = rpot(n)
            endif
          endif
        endif
      enddo
!
!  If no valid potentials then there is no need
!  to continue with this pair, unless this is
!  a second derivative calculation, in which
!  case there may be a contribution from triangles
!  of interactions.
!
      if (npots.eq.0) cycle jloop
      cut2r = rp*rp
      if (cut2r.gt.cut2p) cut2r = cut2p
      cut2 = cut2r
      rp = sqrt(cut2)
!*********************************
!  Find valid vectors for i - j  *
!*********************************
      call rfind(xcrd,ycrd,zcrd,cut2,0.0_dp,0.0_dp,.false.,.false.,i,j,0_i4,0_i4,0_i4,0_i4,nmolonly,lself,0_i4,nor)
!
!  Loop over valid vectors 
!   
      do ii = 1,nor
!***************************************************
!  Calculate many-body contribution in real space  *
!***************************************************
        deriv(1:3) = 0.0_dp
        r2 = dist(ii)
        xcd = xtmp(ii)
        ycd = ytmp(ii)
        zcd = ztmp(ii)
!
!  Correct for displacements for matom if i = matom
!
        if (i.eq.matom) then
          xcd = xcd - (vmatom(1) - xclat(matom))
          ycd = ycd - (vmatom(2) - yclat(matom))
          zcd = zcd - (vmatom(3) - zclat(matom))
          r2 = xcd**2 + ycd**2 + zcd**2
        endif
        if (j.eq.matom) then
          xcd = xcd + (vmatom(1) - xclat(matom))
          ycd = ycd + (vmatom(2) - yclat(matom))
          zcd = zcd + (vmatom(3) - zclat(matom))
          r2 = xcd**2 + ycd**2 + zcd**2
        endif
        r = sqrt(r2)
!
!  Find cell index
!
        if (abs(cellindex(1,ii)).gt.nd2cell(1).or. &
            abs(cellindex(2,ii)).gt.nd2cell(2).or. &
            abs(cellindex(3,ii)).gt.nd2cell(3)) then
!
!  Find cell index - if outside user range then assign to central cell
!
          ncindm = nd2central
          ncindp = nd2central
        else
!
!  Compute index
!
          ncindp = nd2cellptr(nd2cell(1)+1+cellindex(1,ii),nd2cell(2)+1+cellindex(2,ii),nd2cell(3)+1+cellindex(3,ii))
          ncindm = nd2cellptr(nd2cell(1)+1-cellindex(1,ii),nd2cell(2)+1-cellindex(2,ii),nd2cell(3)+1-cellindex(3,ii))
        endif
        ncindc = nd2central
!***************************************
!  Valid many-body potentials for i-j  *
!***************************************
        if (lMEAM) then
          rhoij(1:maxmeamcomponent) = 0.0_dp
          rhoji(1:maxmeamcomponent) = 0.0_dp
          drhoij(1:3,1:maxmeamcomponent) = 0.0_dp
          drhoji(1:3,1:maxmeamcomponent) = 0.0_dp
        else
          rhoij(1) = 0.0_dp
          rhoji(1) = 0.0_dp
        endif
        drhototij(1:3) = 0.0_dp
        drhototji(1:3) = 0.0_dp
!
        lvalidij = .false.
        if (npots.gt.0) then
          if (lMEAMden) then
            do m = 1,npots
              npot = npotl(m)
              if (nptype(npot).eq.19) then
                if (r.gt.rpot2(npot).and.r.le.rpot(npot)) then
                  lvalidij = .true.
!**********************************
!  Calculate density derivatives  *
!**********************************
                  call meamrho(nati,ntypi,natj,ntypj,r,rpot(npot),xcd,ycd,zcd,rhoij,rhoji,drhoij,drhoji, &
                               drhoijs,drhojis,drhoij2,drhoji2,drhoij2s,drhoji2s,drhoij2m,drhoji2m, &
                               1.0_dp,1.0_dp,.true.,.false.,.true.,.false.,twopot(1,npot))
                endif
              endif
            enddo
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
              npartial = 0
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
                  xjk0 = xik0 - xtmp(ii)
                  yjk0 = yik0 - ytmp(ii)
                  zjk0 = zik0 - ztmp(ii)
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
                    call meamscreen(neamspeci,neamspecj,neamspeck,r2,vectorpair%distance_pair1(kvec), &
                                    vectorpair%distance_pair2(kvec),Sikj,dSikjdr,lpartial,.true.)
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
                      if (lpartial) then
!
!  If this atom has a screening factor between 0 and 1, we need to keep track of it since it will generate non-zero derivatives
!
                        npartial = npartial + 1
                        if (npartial.gt.partial%sa_maxdim) then
                          call changemaxsa(partial,npartial)
                        endif
                        partial%sa_atom(npartial) = k
                        partial%sa_rij(npartial) = sqrt(r2)
                        partial%sa_rik(npartial) = sqrt(vectorpair%distance_pair1(kvec))
                        partial%sa_rjk(npartial) = sqrt(vectorpair%distance_pair2(kvec))
                        partial%sa_xik(npartial) = vectorpair%xvector_pair1(kvec)
                        partial%sa_yik(npartial) = vectorpair%yvector_pair1(kvec)
                        partial%sa_zik(npartial) = vectorpair%zvector_pair1(kvec)
                        partial%sa_xjk(npartial) = vectorpair%xvector_pair2(kvec)
                        partial%sa_yjk(npartial) = vectorpair%yvector_pair2(kvec)
                        partial%sa_zjk(npartial) = vectorpair%zvector_pair2(kvec)
                        partial%sa_Sikj(npartial) = Sikj
                        partial%sa_dSikjdr(1:3,npartial) = dSikjdr(1:3)
                      endif
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
!
!  End of screening function
!
            endif
!
!  Only do remainder of work if the screening factor is non-zero
!
            if (lnonzeroSij) then
              do m = 1,npots
                npot = npotl(m)
!
!  Pair potential contribution
!
                if (nptype(npot).eq.45.or.nptype(npot).eq.55) then
                  if (r.gt.rpot2(npot).and.r.le.rpot(npot)) then
                    lvalidij = .true.
                    rk = 1.0_dp/r
                    call baskes(nati,ntypi,neamspeci,natj,ntypj,neamspecj,npot,r,rk,1.0_dp,ebas,d1bas,d2bas,.true.,.false.)
                    ebas  = ebas*ofct
                    d1bas = rk*d1bas*ofct
                    deriv(1) = deriv(1) + d1bas*xcd*Sij
                    deriv(2) = deriv(2) + d1bas*ycd*Sij
                    deriv(3) = deriv(3) + d1bas*zcd*Sij
!
                    call psiscreenderv(i,j,npartial,partial,xcd,ycd,zcd,ebas,Sij,.false.)
                  endif
                endif
!
!  Density contribution
!
                if (nptype(npot).eq.19) then
                  if (lanyMEAMscreen) then
                    if (npartial.gt.0) then
!
!  Compute derivative contributions of the screening function
!
                      call meamtotalrhoscreenderv(neamspeci,npartial,partial,xcd,ycd,zcd,scrho(1,matom), &
                                                  rhoi,rscrhoi,rhoij,Sij,.false.)
!
                      do np = 1,npartial
                        k = partial%sa_atom(np)
!
!  i-j contribution
!
                        if (i.ne.matom) then
                          d1cell(1,i,ncindm) = d1cell(1,i,ncindm) - partial%sa_drhototij(1,np)*ofct
                          d1cell(2,i,ncindm) = d1cell(2,i,ncindm) - partial%sa_drhototij(2,np)*ofct
                          d1cell(3,i,ncindm) = d1cell(3,i,ncindm) - partial%sa_drhototij(3,np)*ofct
                        endif
                        if (j.ne.matom) then
                          d1cell(1,j,ncindp) = d1cell(1,j,ncindp) + partial%sa_drhototij(1,np)*ofct
                          d1cell(2,j,ncindp) = d1cell(2,j,ncindp) + partial%sa_drhototij(2,np)*ofct
                          d1cell(3,j,ncindp) = d1cell(3,j,ncindp) + partial%sa_drhototij(3,np)*ofct
                        endif
!
!  i-k contribution
!
                        if (i.ne.matom) then
                          d1cell(1,i,ncindc) = d1cell(1,i,ncindc) - partial%sa_drhototik(1,np)*ofct
                          d1cell(2,i,ncindc) = d1cell(2,i,ncindc) - partial%sa_drhototik(2,np)*ofct
                          d1cell(3,i,ncindc) = d1cell(3,i,ncindc) - partial%sa_drhototik(3,np)*ofct
                        endif
                        if (k.ne.matom) then
                          d1cell(1,k,ncindc) = d1cell(1,k,ncindc) + partial%sa_drhototik(1,np)*ofct
                          d1cell(2,k,ncindc) = d1cell(2,k,ncindc) + partial%sa_drhototik(2,np)*ofct
                          d1cell(3,k,ncindc) = d1cell(3,k,ncindc) + partial%sa_drhototik(3,np)*ofct
                        endif
!
!  j-k contribution
!
                        if (j.ne.matom) then
                          d1cell(1,j,ncindc) = d1cell(1,j,ncindc) - partial%sa_drhototjk(1,np)*ofct
                          d1cell(2,j,ncindc) = d1cell(2,j,ncindc) - partial%sa_drhototjk(2,np)*ofct
                          d1cell(3,j,ncindc) = d1cell(3,j,ncindc) - partial%sa_drhototjk(3,np)*ofct
                        endif
                        if (k.ne.matom) then
                          d1cell(1,k,ncindc) = d1cell(1,k,ncindc) + partial%sa_drhototjk(1,np)*ofct
                          d1cell(2,k,ncindc) = d1cell(2,k,ncindc) + partial%sa_drhototjk(2,np)*ofct
                          d1cell(3,k,ncindc) = d1cell(3,k,ncindc) + partial%sa_drhototjk(3,np)*ofct
                        endif
                      enddo
!
                      call meamtotalrhoscreenderv(neamspecj,npartial,partial,xcd,ycd,zcd,scrho(1,j), &
                                                  rhoj,rscrhoj,rhoji,Sij,.false.)
!
                      do np = 1,npartial
                        k = partial%sa_atom(np)
!
!  i-j contribution
!   
                        if (i.ne.matom) then
                          d1cell(1,i,ncindm) = d1cell(1,i,ncindm) - partial%sa_drhototij(1,np)*ofct
                          d1cell(2,i,ncindm) = d1cell(2,i,ncindm) - partial%sa_drhototij(2,np)*ofct
                          d1cell(3,i,ncindm) = d1cell(3,i,ncindm) - partial%sa_drhototij(3,np)*ofct
                        endif
                        if (j.ne.matom) then
                          d1cell(1,j,ncindp) = d1cell(1,j,ncindp) + partial%sa_drhototij(1,np)*ofct
                          d1cell(2,j,ncindp) = d1cell(2,j,ncindp) + partial%sa_drhototij(2,np)*ofct
                          d1cell(3,j,ncindp) = d1cell(3,j,ncindp) + partial%sa_drhototij(3,np)*ofct
                        endif
!
!  i-k contribution
!
                        if (i.ne.matom) then
                          d1cell(1,i,ncindc) = d1cell(1,i,ncindc) - partial%sa_drhototik(1,np)*ofct
                          d1cell(2,i,ncindc) = d1cell(2,i,ncindc) - partial%sa_drhototik(2,np)*ofct
                          d1cell(3,i,ncindc) = d1cell(3,i,ncindc) - partial%sa_drhototik(3,np)*ofct
                        endif
                        if (k.ne.matom) then
                          d1cell(1,k,ncindc) = d1cell(1,k,ncindc) + partial%sa_drhototik(1,np)*ofct
                          d1cell(2,k,ncindc) = d1cell(2,k,ncindc) + partial%sa_drhototik(2,np)*ofct
                          d1cell(3,k,ncindc) = d1cell(3,k,ncindc) + partial%sa_drhototik(3,np)*ofct
                        endif
!
!  j-k contribution
!
                        if (j.ne.matom) then
                          d1cell(1,j,ncindc) = d1cell(1,j,ncindc) - partial%sa_drhototjk(1,np)*ofct
                          d1cell(2,j,ncindc) = d1cell(2,j,ncindc) - partial%sa_drhototjk(2,np)*ofct
                          d1cell(3,j,ncindc) = d1cell(3,j,ncindc) - partial%sa_drhototjk(3,np)*ofct
                        endif
                        if (k.ne.matom) then
                          d1cell(1,k,ncindc) = d1cell(1,k,ncindc) + partial%sa_drhototjk(1,np)*ofct
                          d1cell(2,k,ncindc) = d1cell(2,k,ncindc) + partial%sa_drhototjk(2,np)*ofct
                          d1cell(3,k,ncindc) = d1cell(3,k,ncindc) + partial%sa_drhototjk(3,np)*ofct
                        endif
                      enddo
                    endif
                  endif
!
!  Scale density and derivatives by screening factor
!
                  drhoij(1:3,1:maxmeamcomponent) = Sij*drhoij(1:3,1:maxmeamcomponent)
                  drhoji(1:3,1:maxmeamcomponent) = Sij*drhoji(1:3,1:maxmeamcomponent)
!
!  Compute total derivatives of MEAM density
!
                  call meamtotalrhoderv(neamspeci,scrho(1,i),rhoi,drhoij,drhototij,drhoijs,drhototijs, &
                                        drhoij2,drhototij2,drhoij2s,drhototij2s,drhoij2m,drhototij2m, &
                                        .false.,.false.)
                  call meamtotalrhoderv(neamspecj,scrho(1,j),rhoj,drhoji,drhototji,drhojis,drhototjis, &
                                        drhoji2,drhototji2,drhoji2s,drhototji2s,drhoji2m,drhototji2m, &
                                        .false.,.false.)
                endif
              enddo
            endif
          else
            do m = 1,npots
              npot = npotl(m)
              if (nptype(npot).eq.19) then
                if (r.gt.rpot2(npot).and.r.le.rpot(npot)) then
                  lvalidij = .true.
                  call eamrho(nati,ntypi,natj,ntypj,r,rpot(npot),xcd,ycd,zcd,rhoij,rhoji,drhototij,drhototji, &
                              drhototijs,drhototjis,drhototij2,drhototji2,drhototij2s,drhototji2s, &
                              drhototij2m,drhototji2m,drhototij3,drhototji3,1.0_dp,1.0_dp,.false.,.true.,.false.,.false., &
                              twopot(1,npot))
                endif
              endif
            enddo
          endif
!
!  Combine derivative terms
!
          deriv(1:3) = deriv(1:3) + (rscrhoi*drhototij(1:3) + rscrhoj*drhototji(1:3))*ofct
        endif
        if (lvalidij) then
!******************************
!  Internal first derivatives *
!******************************
!
!  Only need to do internal derivatives if i not equals j and i = matom
!
          if (i.ne.j) then
            if (i.ne.matom) then
              d1cell(1,i,ncindm) = d1cell(1,i,ncindm) - deriv(1)
              d1cell(2,i,ncindm) = d1cell(2,i,ncindm) - deriv(2)
              d1cell(3,i,ncindm) = d1cell(3,i,ncindm) - deriv(3)
            endif
            if (j.ne.matom) then
              d1cell(1,j,ncindp) = d1cell(1,j,ncindp) + deriv(1)
              d1cell(2,j,ncindp) = d1cell(2,j,ncindp) + deriv(2)
              d1cell(3,j,ncindp) = d1cell(3,j,ncindp) + deriv(3)
            endif
          endif
        endif
!**************************************
!  End of valid distance i-j section  *
!**************************************
      enddo
    enddo jloop
  enddo iloop
!
!  Free local memory
!
  deallocate(npotl,stat=status)
  if (status/=0) call deallocate_error('many1fc','npotl')
!
!  Unscale density
!
  call eamscalescrho(-1_i4)
!
!  Exit point
!
1000 continue
!
!  Timing
!
  time2 = g_cpu_time()
  tmany = tmany + time2 - time1
#ifdef TRACE
  call trace_out('many1fc')
#endif
!
  return
  end
