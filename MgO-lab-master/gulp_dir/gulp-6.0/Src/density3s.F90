  subroutine density3s(eatom,esregion12,esregion2,eattach)
!
!  Subroutine for calculating MEAM electron density, including the screening function, using spatial decomposition
!
!  lfreeze = if .true. only do derivatives for atoms with a degree of freedom
!
!   4/09 Created from realmd3s
!   4/09 Screening function added
!   8/14 MEAM screening made species specific
!   8/14 Energy passed as argument so that pair potential contribution can be added
!   8/14 All MEAM species (i,j,k) now passed to meanscreen
!   8/14 neamspec passed psibaskes
!   9/14 Error in position of endif that was causing use of uninitialised variable
!        has been corrected
!   1/15 Name of subroutine corrected in memory allocation calls
!   2/18 Trace added
!   5/19 Site energy added
!   5/19 Twobody energy print out added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!   8/20 if statement for not cell buffer moved to avoid unnecessary calculation
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
!  Julian Gale, CIC, Curtin University, August 2020
!
  use configurations, only : lsliceatom, nregionno, nregiontype, QMMMmode, nregions
  use control
  use current
  use datatypes
  use eam,            only : lMEAM, lMEAMscreen, maxmeamcomponent, meam_Cmax
  use eam,            only : lanyMEAMscreen, neamfnspec
  use element
  use energies,       only : siteenergy, eregion2region
  use general,        only : smallself
  use iochannels,     only : ioout
  use optimisation
  use parallel
  use spatial
  use sutton
  use symmetry
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
  real(dp),                     intent(inout)  :: eattach
  real(dp),                     intent(inout)  :: esregion12
  real(dp),                     intent(inout)  :: esregion2
!
!  Local variables
!
  integer(i4), dimension(:), allocatable       :: npotl
  integer(i4)                                  :: i
  integer(i4)                                  :: ic
  integer(i4)                                  :: indij
  integer(i4)                                  :: ii
  integer(i4)                                  :: imx
  integer(i4)                                  :: imy
  integer(i4)                                  :: imz
  integer(i4)                                  :: ind
  integer(i4)                                  :: ind2
  integer(i4)                                  :: indn
  integer(i4)                                  :: ixyz
  integer(i4)                                  :: ix
  integer(i4)                                  :: iy
  integer(i4)                                  :: iz
  integer(i4)                                  :: j
  integer(i4)                                  :: jc
  integer(i4)                                  :: jj
  integer(i4)                                  :: jmx
  integer(i4)                                  :: jmy
  integer(i4)                                  :: jmz
  integer(i4)                                  :: jndn
  integer(i4)                                  :: k
  integer(i4)                                  :: kc
  integer(i4)                                  :: kk
  integer(i4)                                  :: maxxy
  integer(i4)                                  :: maxx
  integer(i4)                                  :: n
  integer(i4)                                  :: n1i
  integer(i4)                                  :: n1j
  integer(i4)                                  :: n1k
  integer(i4)                                  :: nat1
  integer(i4)                                  :: nat2
  integer(i4)                                  :: nati
  integer(i4)                                  :: natj
  integer(i4)                                  :: neamfnspec2
  integer(i4)                                  :: neamspeci
  integer(i4)                                  :: neamspecj
  integer(i4)                                  :: neamspeck
  integer(i4)                                  :: ni
  integer(i4)                                  :: nj
  integer(i4)                                  :: nk
  integer(i4)                                  :: np
  integer(i4)                                  :: npot
  integer(i4)                                  :: npots
  integer(i4)                                  :: nregioni
  integer(i4)                                  :: nregionj
  integer(i4)                                  :: nregiontypi
  integer(i4)                                  :: nregiontypj
  integer(i4)                                  :: nsplower(3)
  integer(i4)                                  :: nspupper(3)
  integer(i4)                                  :: ntyp1
  integer(i4)                                  :: ntyp2
  integer(i4)                                  :: ntypi
  integer(i4)                                  :: ntypj
  integer(i4)                                  :: status
  logical                                      :: lattach
  logical                                      :: lmatch
  logical                                      :: lnonzeroSij
  logical                                      :: lopi
  logical                                      :: lopj
  logical                                      :: lorder12loc
  logical                                      :: lpartial
  logical                                      :: lreg2one
  logical                                      :: lreg2pair
  logical                                      :: lslicei
  logical                                      :: lslicej
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: cut2
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
  real(dp)                                     :: dSikjdr(3)       ! Dummy argument for call to meamscreen - not used here
  real(dp)                                     :: Sij
  real(dp)                                     :: Sikj
  real(dp),    dimension(:), allocatable       :: sum
  real(dp),    dimension(:), allocatable       :: sum2
  real(dp)                                     :: time1
  real(dp)                                     :: time2
  real(dp)                                     :: tsum0
  real(dp)                                     :: tsuml
  real(dp)                                     :: xi
  real(dp)                                     :: yi
  real(dp)                                     :: zi
  real(dp)                                     :: xij0
  real(dp)                                     :: yij0
  real(dp)                                     :: zij0
  real(dp)                                     :: xji
  real(dp)                                     :: yji
  real(dp)                                     :: zji
  real(dp)                                     :: xki
  real(dp)                                     :: yki
  real(dp)                                     :: zki
  real(dp)                                     :: xkj
  real(dp)                                     :: ykj
  real(dp)                                     :: zkj
!
!  Check that this call is needed
!
  if (.not.lsuttonc) return
#ifdef TRACE
  call trace_in('density3s')
#endif
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
  if (status/=0) call outofmemory('density3s','npotl')
  allocate(sum(numat),stat=status)
  if (status/=0) call outofmemory('density3s','sum')
  allocate(sum2(numat),stat=status)
  if (status/=0) call outofmemory('density3s','sum2')
!*****************************
!  Compute density for MEAM  *
!*****************************
  maxxy = nspcell(1)*nspcell(2)
  maxx  = nspcell(1)
!
!  Loop over all local spatial cells except buffer regions
!
  do ixyz = 1,ncellpernode
    if (.not.lbuffercell(ixyz)) then
      ind = ncellnodeptr(ixyz)
      ind2 = ind - 1
      iz = ind2/maxxy
      ind2 = ind2 - maxxy*iz
      iy = ind2/maxx
      ix = ind2 - maxx*iy + 1
      iy = iy + 1
      iz = iz + 1
!
!  Set cell search bounds
!  
      nspupper(1) = min(ix+ncellsearch(1),nspcell(1))
      nspupper(2) = min(iy+ncellsearch(2),nspcell(2))
      nspupper(3) = min(iz+ncellsearch(3),nspcell(3))
      nsplower(1) = max(ix-ncellsearch(1),1)
      nsplower(2) = max(iy-ncellsearch(2),1)
      nsplower(3) = max(iz-ncellsearch(3),1)
!
!  Get number of atoms in this cell
!
      ni = nspcellat(ind)
      n1i = nspcellat1ptr(ind)
!
!  Outer loop over atoms within this cell
!
      do ii = 1,ni
        i = nspcellatptr(n1i+ii)
        ic = nspcellatptrcell(n1i+ii)
!
!  Set coordinates of atom i
!
        xi = xinbox(i) + xvec2cell(ic)
        yi = yinbox(i) + yvec2cell(ic)
        zi = zinbox(i) + zvec2cell(ic)
!
!  Set other properties of atom i
!
        nati = nat(i)
        ntypi = nftype(i)
        oci = occuf(i)
        nregioni = nregionno(nsft+nrelf2a(i))
        nregiontypi = nregiontype(nregioni,ncf)
        lopi = (.not.lfreeze.or.lopf(nrelf2a(i)))
        lslicei = lsliceatom(nsft + nrelf2a(i))
        neamspeci = neamfnspecptr(i)
!
!  Loop over neighbouring cells
!
        do imz = nsplower(3),nspupper(3)
          do imy = nsplower(2),nspupper(2)
            do imx = nsplower(1),nspupper(1)
              indn = (imz-1)*maxxy + (imy-1)*maxx + imx
!
!  Loop over atoms within neighbouring cells
!
              nj = nspcellat(indn)
              n1j = nspcellat1ptr(indn)
              ebassum = 0.0_dp
              jloop: do jj = 1,nj
                j = nspcellatptr(n1j+jj)
                jc = nspcellatptrcell(n1j+jj)
!
!  Only calculate lower-half triangular interactions
!
                if (j.le.i) then
!
!  Freezing flag
!
                  lopj = (.not.lfreeze.or.lopf(nrelf2a(j)))
                  if (.not.lopi.and..not.lopj) cycle jloop
!
!  Set coordinate differences and calculate square of distance
!  
                  xji = xvec2cell(jc) + xinbox(j) - xi
                  yji = yvec2cell(jc) + yinbox(j) - yi
                  zji = zvec2cell(jc) + zinbox(j) - zi
                  r2 = xji*xji + yji*yji + zji*zji
!
!  Set species type parameters for atom j
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
!  Set remaining properties for atom type j
!
                  ocj = occuf(j)
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
!  If no valid potentials and charge product is zero
!  then no need to search for distances
!
                  if (npots.eq.0) cycle jloop
                  cut2 = rp*rp
!
!  Check distance against potential cutoff
!
                  if (r2.gt.cut2.or.r2.lt.smallself) cycle jloop
!*********************************
!  Calculate unscreened density  *
!*********************************
                  r = sqrt(r2)
                  call twoden1(1_i4,1_i4,r,xji,yji,zji,npots,npotl,sctrm1,sctrm2,lorder12loc)
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
!  Loop over neighbouring cells to search for atoms that may contribute to the screening
!
                    jmzloop: do jmz = nsplower(3),nspupper(3)
                      do jmy = nsplower(2),nspupper(2)
                        do jmx = nsplower(1),nspupper(1)
                          jndn = (jmz-1)*maxxy + (jmy-1)*maxx + jmx
!
!  Loop over atoms within neighbouring cells
!
                          nk = nspcellat(jndn)
                          n1k = nspcellat1ptr(jndn)
                          kloop: do kk = 1,nk
                            k = nspcellatptr(n1k+kk)
                            kc = nspcellatptrcell(n1k+kk)
                            neamspeck = neamfnspecptr(k)
                            if (lMEAMscreen(indij,neamspeck)) then
!
!  Set vectors between atoms
!
                              xki = xvec2cell(kc) + xinbox(k) - xi
                              yki = yvec2cell(kc) + yinbox(k) - yi
                              zki = zvec2cell(kc) + zinbox(k) - zi
                              xkj = xki - xji
                              ykj = yki - yji
                              zkj = zki - zji
                              xij0 = 0.5_dp*(xki + xkj)
                              yij0 = 0.5_dp*(yki + ykj)
                              zij0 = 0.5_dp*(zki + zkj)
!
!  Compute square of distance to i-j mid point
!
                              r2ijmid = xij0*xij0 + yij0*yij0 + zij0*zij0
                              if (r2ijmid.lt.rcut2) then
!
!  Complete distances
!
                                r2ik = xki*xki + yki*yki + zki*zki
                                r2jk = xkj*xkj + ykj*ykj + zkj*zkj
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
!
!  If screening factor has gone to zero then we can exit all loops over k
!
                                if (.not.lnonzeroSij) exit jmzloop
                              endif
                            endif
                          enddo kloop
!
!  End of loops over cells for k
!
                        enddo 
                      enddo 
                    enddo jmzloop
                  endif
!
                  if (lnonzeroSij) then
                    if (lMEAM) then
                      if (lorder12loc) then
                        scrho(1:maxmeamcomponent,i) = scrho(1:maxmeamcomponent,i) + sctrm1(1:maxmeamcomponent)*ocj*Sij
                        scrho(1:maxmeamcomponent,j) = scrho(1:maxmeamcomponent,j) + sctrm2(1:maxmeamcomponent)*oci*Sij
                        if (lattach) then
                          scrho12(1:maxmeamcomponent,i) = scrho12(1:maxmeamcomponent,i) + sctrm1(1:maxmeamcomponent)*ocj*Sij
                          scrho12(1:maxmeamcomponent,j) = scrho12(1:maxmeamcomponent,j) + sctrm2(1:maxmeamcomponent)*oci*Sij
                        endif
                      else
                        scrho(1:maxmeamcomponent,i) = scrho(1:maxmeamcomponent,i) + sctrm2(1:maxmeamcomponent)*ocj*Sij
                        scrho(1:maxmeamcomponent,j) = scrho(1:maxmeamcomponent,j) + sctrm1(1:maxmeamcomponent)*oci*Sij
                        if (lattach) then
                          scrho12(1:maxmeamcomponent,i) = scrho12(1:maxmeamcomponent,i) + sctrm2(1:maxmeamcomponent)*ocj*Sij
                          scrho12(1:maxmeamcomponent,j) = scrho12(1:maxmeamcomponent,j) + sctrm1(1:maxmeamcomponent)*oci*Sij
                        endif
                      endif
                    else
                      if (lorder12loc) then
                        scrho(1,i) = scrho(1,i) + sctrm1(1)*ocj*Sij
                        scrho(1,j) = scrho(1,j) + sctrm2(1)*oci*Sij
                        if (lattach) then
                          scrho12(1,i) = scrho12(1,i) + sctrm1(1)*ocj*Sij
                          scrho12(1,j) = scrho12(1,j) + sctrm2(1)*oci*Sij
                        endif
                      else
                        scrho(1,i) = scrho(1,i) + sctrm2(1)*ocj*Sij
                        scrho(1,j) = scrho(1,j) + sctrm1(1)*oci*Sij
                        if (lattach) then
                          scrho12(1,i) = scrho12(1,i) + sctrm2(1)*ocj*Sij
                          scrho12(1,j) = scrho12(1,j) + sctrm1(1)*oci*Sij
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
                        call baskes(nati,ntypi,neamspeci,natj,ntypj,neamspecj,npot,r,rk,1.0_dp, &
                                       ebas,d1bas,d2bas,.false.,.false.)
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
                endif
!
                if (lPrintTwo) then
                  write(ioout,'(4x,i12,i12,f22.10)') i,j,ebassum
                endif
!
!  End loop over atom j
!
              enddo jloop
!
!  End loops over neighbouring cells
!
            enddo
          enddo
        enddo
!
!  End loop over atom i
!
      enddo
!
!  End checks on whether cell is required
!
    endif
!
!  End loop over cells on node
!
  enddo
!****************
!  Global sums  *
!****************
  tsum0 = g_cpu_time()
  if (lsuttonc.and.nprocs.gt.1) then
    if (lMEAM) then
      do j = 1,maxmeamcomponent
        do i = 1,numat
          sum2(i) = scrho(j,i) 
        enddo
        call sumall(sum2,sum,numat,"density3s","scrho")
        do i = 1,numat
          scrho(j,i) = sum(i)
        enddo
        do i = 1,numat
          sum2(i) = scrho12(j,i) 
        enddo
        call sumall(sum2,sum,numat,"density3s","scrho12")
        do i = 1,numat
          scrho12(j,i) = sum(i)
        enddo
      enddo
    else
      do i = 1,numat
        sum2(i) = scrho(1,i) 
      enddo
      call sumall(sum2,sum,numat,"density3s","scrho")
      do i = 1,numat
        scrho(1,i) = sum(i)
      enddo
      do i = 1,numat
        sum2(i) = scrho12(1,i) 
      enddo
      call sumall(sum2,sum,numat,"density3s","scrho12")
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
  if (status/=0) call deallocate_error('density3s','sum2')
  deallocate(sum,stat=status)
  if (status/=0) call deallocate_error('density3s','sum')
  deallocate(npotl,stat=status)
  if (status/=0) call deallocate_error('density3s','npotl')
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
  call trace_out('density3s')
#endif
!
  return
  end
