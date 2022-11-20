  subroutine realmd0s(eatom,ereal,eqeq,lgrad1)
!
!  Subroutine for calculating real space energy and gradients using a spatial decomposition
!  algorithm for a finite cluster => 0 D system for MD
!
!  lfreeze = if .true. only do derivatives for atoms with a degree of freedom
!
!   9/15 Created from realmd3s
!   2/18 Trace added
!   5/18 Modified for q ranges
!   8/19 Short range damping of polarisation added
!   8/19 Trap for neemrptr being zero added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!   2/20 Breathing shell contribution now made conditional on lopi
!   3/20 Correction to polarisation handling - polfct check added
!   3/20 Cutoffs introduced for case where two non-charged particles interact
!   8/20 if statement for not cell buffer moved to avoid unnecessary calculation
!  10/20 eqeq now included in site energy and eregion2region
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
!  Julian Gale, CIC, Curtin University, October 2020
!
  use configurations, only : lbsmat, nregionno, nregiontype, QMMMmode
  use g_constants
  use control
  use current
  use datatypes
  use derivatives,    only : xdrv, ydrv, zdrv, raderv
  use derivatives,    only : xregdrv, yregdrv, zregdrv
  use eam,            only : lMEAMden
  use eemdata
  use element
  use energies,       only : siteenergy, eregion2region
  use general,        only : cutw
  use iochannels,     only : ioout
  use kspace
  use mdlogic
  use molecule
  use optimisation
  use parallel
  use polarise
  use potentialxyz
  use realvectors,    only : derivqd
  use shells
  use spatial
  use sutton
  use symmetry
  use thresholds,     only : thresh_q
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  use two
  implicit none
!
!  Passed variables
!
  logical,     intent(in)                      :: lgrad1
  real(dp),    intent(inout)                   :: eatom
  real(dp),    intent(inout)                   :: eqeq
  real(dp),    intent(inout)                   :: ereal
!
!  Local variables
!
  integer(i4), dimension(:), allocatable       :: npotl
  integer(i4)                                  :: i
  integer(i4)                                  :: ic
  integer(i4)                                  :: icx
  integer(i4)                                  :: icy
  integer(i4)                                  :: icz
  integer(i4)                                  :: ii
  integer(i4)                                  :: imx
  integer(i4)                                  :: imy
  integer(i4)                                  :: imz
  integer(i4)                                  :: indmi
  integer(i4)                                  :: indmj
  integer(i4)                                  :: ind
  integer(i4)                                  :: ind2
  integer(i4)                                  :: indn
  integer(i4)                                  :: ix
  integer(i4)                                  :: iy
  integer(i4)                                  :: iz
  integer(i4)                                  :: ixyz
  integer(i4)                                  :: ixx
  integer(i4)                                  :: iyy
  integer(i4)                                  :: izz
  integer(i4)                                  :: j
  integer(i4)                                  :: jxx
  integer(i4)                                  :: jyy
  integer(i4)                                  :: jzz
  integer(i4)                                  :: jc
  integer(i4)                                  :: jcx
  integer(i4)                                  :: jcy
  integer(i4)                                  :: jcz
  integer(i4)                                  :: jj
  integer(i4)                                  :: maxxy
  integer(i4)                                  :: maxx
  integer(i4)                                  :: n
  integer(i4)                                  :: n1i
  integer(i4)                                  :: n1j
  integer(i4)                                  :: nat1
  integer(i4)                                  :: nat2
  integer(i4)                                  :: nati
  integer(i4)                                  :: natj
  integer(i4)                                  :: nbtypeij
  integer(i4)                                  :: nbtypeij2
  integer(i4)                                  :: ni
  integer(i4)                                  :: nj
  integer(i4)                                  :: nmi
  integer(i4)                                  :: nmj
  integer(i4)                                  :: npots
  integer(i4)                                  :: nqri
  integer(i4)                                  :: nqrj
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
  logical                                      :: l2bonds
  logical                                      :: l3bonds
  logical                                      :: lbonded
  logical                                      :: lcspair
  logical                                      :: lgrad1p
  logical                                      :: lmatch
  logical                                      :: lmolok
  logical                                      :: lneedmol
  logical                                      :: lopi
  logical                                      :: lopj
  logical                                      :: lorder12loc
  logical                                      :: lptrmol
  logical                                      :: lQMMMelectro
  real(dp)                                     :: apt
  real(dp)                                     :: bpt
  real(dp)                                     :: c6tot
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: cut2
  real(dp)                                     :: cut2e
  real(dp)                                     :: cut2p
  real(dp)                                     :: cut2r
  real(dp)                                     :: cut2s
  real(dp)                                     :: d0i
  real(dp)                                     :: d0j
  real(dp)                                     :: d1i
  real(dp)                                     :: d1j
  real(dp)                                     :: d2i2
  real(dp)                                     :: d2ij
  real(dp)                                     :: d2j2
  real(dp)                                     :: deriv
  real(dp)                                     :: deriv2
  real(dp)                                     :: deriv3
  real(dp)                                     :: derive0
  real(dp)                                     :: derive
  real(dp)                                     :: derive2
  real(dp)                                     :: derive3
  real(dp)                                     :: dist
  real(dp)                                     :: eatom2
  real(dp)                                     :: ec62
  real(dp)                                     :: eqeq2  
  real(dp)                                     :: ereal2
  real(dp)                                     :: esum 
  real(dp)                                     :: etrm
  real(dp)                                     :: etrm1
  real(dp)                                     :: etrm2
  real(dp)                                     :: factor
  real(dp)                                     :: fct
  real(dp)                                     :: oci
  real(dp)                                     :: ocj
  real(dp)                                     :: ofct
  real(dp)                                     :: ospfct
  real(dp)                                     :: polfct
  real(dp)                                     :: qli
  real(dp)                                     :: qlj
  real(dp)                                     :: r2
  real(dp)                                     :: radi
  real(dp)                                     :: radj
  real(dp)                                     :: radsum
  real(dp)                                     :: rderiv
  real(dp)                                     :: rdiff
  real(dp)                                     :: rp
  real(dp)                                     :: rtrm1
  real(dp)                                     :: rtrm2
  real(dp)                                     :: rtrm3
  real(dp)                                     :: rtrm32
  real(dp)                                     :: sctrm1
  real(dp)                                     :: sctrm2
  real(dp)                                     :: small
  real(dp)                                     :: small2
  real(dp),    dimension(:), allocatable       :: sum
  real(dp),    dimension(:), allocatable       :: sum2
  real(dp)                                     :: time1
  real(dp)                                     :: time2
  real(dp)                                     :: tsum0
  real(dp)                                     :: tsuml
  real(dp)                                     :: xi
  real(dp)                                     :: yi
  real(dp)                                     :: zi
  real(dp)                                     :: xji
  real(dp)                                     :: yji
  real(dp)                                     :: zji
#ifdef TRACE
  call trace_in('realmd0s')
#endif
!
  time1 = g_cpu_time()
  tsuml = 0.0_dp
!
!  Initialise local variables
!
  lgrad1p = (lgrad1.or.lpolar)
  small = 1.0d-12
  small2 = 1.0d-2
!
!  Set up cutoffs
!
  if (lwolf) then
    cut2e = cutw*cutw
  else
    cut2e = 1.0d10
  endif
  cut2p = cutp*cutp
  cut2s = cuts*cuts
!
!  Allocate local memory
!
  allocate(npotl(npote),stat=status)
  if (status/=0) call outofmemory('realmd0s','npotl')
  allocate(sum(numat),stat=status)
  if (status/=0) call outofmemory('realmd0s','sum')
  allocate(sum2(numat),stat=status)
  if (status/=0) call outofmemory('realmd0s','sum2')
!
  if (.not.lnoreal) then
!
!  Opening banner for energy decomposition
!
    if (lPrintTwo) then
      call mpbarrier
      if (ioproc) then
        write(ioout,'(''--------------------------------------------------------------------------------'')')
        write(ioout,'(''  Two : Atom No. 1  Atom No. 2    Short-range energy (eV)   Coulomb energy (eV) '')')
        write(ioout,'(''--------------------------------------------------------------------------------'')')
      endif
    endif
!***************************************************************
!  Atomistic and real space electrostatic component of energy  *
!***************************************************************
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
          call ind2toijk(ic,icx,icy,icz)
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
          qli = qf(i)
          oci = occuf(i)
          nregioni = nregionno(nsft+nrelf2a(i))
          nregiontypi = nregiontype(nregioni,ncf)
          lopi = (.not.lfreeze.or.lopf(nrelf2a(i)))
          if (lbsmat(nrelf2a(i)+nsft)) then
            radi = radf(i)
          else
            radi = 0.0_dp
          endif
          if (leem.and.lmultiqrange) then
            nqri = nqrnow(neemrptr(i))
          else
            nqri = 1
          endif
!
!  Molecule handling for atom i
!
          if (lmol) then
            nmi = natmol(i)
            indmi = nmolind(i)
          endif
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
                jloop: do jj = 1,nj
                  j = nspcellatptr(n1j+jj)
                  jc = nspcellatptrcell(n1j+jj)
                  call ind2toijk(jc,jcx,jcy,jcz)
!
!  Only calculate lower-half triangular interactions
!
                  if (j.lt.i) then
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
!
!  QM/MM handling : i & j are both QM atoms => exclude
!
                    if (QMMMmode(ncf).gt.0) then
                      if (nregiontypi.eq.1.and.nregiontypj.eq.1) cycle jloop
                    endif
!           
!  QM/MM : Set electrostatic embedding flag : If either i or j are QM atoms => exclude electrostatics
!       
                    lQMMMelectro = (QMMMmode(ncf).eq.2.and.(nregiontypi.eq.1.or.nregiontypj.eq.1))
!
!  Set remaining properties for atom type j
!
                    qlj = qf(j)
                    ocj = occuf(j)
                    if (lbsmat(nsft+nrelf2a(j))) then
                      radj = radf(j)
                    else
                      radj = 0.0_dp
                    endif
                    radsum = radi + radj
                    ofct = oci*ocj
                    fct = ofct*angstoev
                    factor = qli*qlj*fct
                    if (lpolar) then
                      polfct = abs(qli*dpolar(j)) + abs(qlj*dpolar(i))
                    else
                      polfct = 0.0_dp
                    endif
                    if (leem.and.lmultiqrange) then
                      nqrj = nqrnow(neemrptr(j))
                    else
                      nqrj = 1
                    endif
!
!  Molecule handling for atom j
!
                    if (lmol) then
                      nmj = natmol(j)
                      indmj = nmolind(j)
                    endif
!
!  Possible core-shell flag
!
                    lcspair = (abs(nat1-nat2).eq.maxele.or.(oci+ocj).lt.1.0001_dp)
                    if (lcspair) then
                      ospfct = oci
                    else
                      ospfct = ofct
                    endif
!
!  Locate potential number
!  Check whether potential requires specific types
!  Calculate sum of all dispersion terms for pair
!
                    rp = 0.0_dp
                    npots = 0
                    c6tot = 0.0_dp
                    lneedmol = (lmol.and..not.lmolq)
                    do n = 1,npote
                      if (lmatch(nat1,ntyp1,nspec1(n),nptyp1(n),.true.)) then
                        if (lmatch(nat2,ntyp2,nspec2(n),nptyp2(n),.true.)) then
                          npots = npots + 1
                          npotl(npots) = n
                          if (mmexc(n).gt.0.or.(lintra(n).and..not.linter(n).or..not.lintra(n).and.linter(n)))  &
                            lneedmol = .true.
                          if (nptype(n).eq.8.or.nptype(n).eq.33) then
                            if (cuts.gt.rp) rp = cuts
                          else
                            if (rpot(n).gt.rp) rp = rpot(n)
                          endif
                        endif
                      endif
                    enddo
!
!  If no valid potentials, charge product is zero, and no polarisation then skip loop
!
                    if (npots.eq.0.and.abs(factor).lt.1.0d-8.and.polfct.lt.1.0d-8) cycle jloop
!
!  Set cutoffs
!
                    cut2r = rp*rp
                    if (cut2r.gt.cut2p) cut2r = cut2p
                    cut2 = cut2r
!
!  If both charges are less than threshold then exclude electrostatics from cutoff
!
                    if (abs(qli*oci)+abs(qlj*ocj).gt.thresh_q) then
                      cut2 = max(cut2,cut2e)
                    endif
!
!  Molecule and bonding checks
!
                    if (lmol) then
                      lmolok = (nmi.eq.nmj.and.nmi.ne.0)
                    else
                      lmolok = .false.
                    endif
!
!  If the Coulomb term and the potentials don't need molecule info for
!  for this pair of atoms, then turn off molecule checking to save time
!
                    if (.not.lneedmol) lmolok = .false.
                    if (lmolok.and.(r2.gt.cut2s.or..not.lcspair)) then
                      ind = indmj - indmi
                      lptrmol = (ind.eq.0)
                      if (.not.lptrmol) then
                        call mindtoijk(indmj,jxx,jyy,jzz)
                        call mindtoijk(indmi,ixx,iyy,izz)
                        jxx = jxx - ixx
                        jyy = jyy - iyy
                        jzz = jzz - izz
                        call samemol(lptrmol,nmi,jxx,jyy,jzz,0_i4,0_i4,0_i4)
                      endif
                      if (lptrmol) then
                        call bonded3(lbonded,l2bonds,l3bonds,nbtypeij,nbtypeij2,i,j,0_i4,0_i4,0_i4)
                      else
                        lbonded   = .false.
                        l2bonds   = .false.
                        l3bonds   = .false.
                        nbtypeij  = 0
                        nbtypeij2 = 0
                      endif
                    else
                      lptrmol   = .false.
                      lbonded   = .false.
                      l2bonds   = .false.
                      l3bonds   = .false.
                      nbtypeij  = 0
                      nbtypeij2 = 0
                    endif
!
                    if (abs(r2-small2).lt.1.0d-12) r2 = small2
                    if (r2.lt.small.or.r2.gt.cut2) then
                      cycle jloop
                    else
!
!  Store vector
!
                      dist = sqrt(r2)
                    endif
!*************************************************
!  Calculate twobody contribution in real space  *
!*************************************************
                    eatom2 = 0.0_dp
                    ereal2 = 0.0_dp
                    ec62 = 0.0_dp
                    dist = sqrt(r2)
                    call twobody1(eatom2,ereal2,ec62,lgrad1p,.false.,.false.,1_i4,1_i4,dist,xji,yji,zji, &
                                  d0i,d0j,deriv,deriv2,deriv3,derive0,derive,derive2,derive3,rderiv,npots, &
                                  npotl,cut2r,cut2e,cut2s,lptrmol,0_i4,factor,ofct,ospfct,radsum,rtrm1, &
                                  rtrm2,rtrm3,rtrm32,sctrm1,sctrm2,qli,qlj,lcspair,.false.,.false., &
                                  lbonded,l2bonds,l3bonds,nbtypeij,nbtypeij2,.false.,lQMMMelectro, &
                                  d1i,d1j,d2i2,d2ij,d2j2)
                    esum = eatom2 + ereal2
                    eatom = eatom + eatom2
                    ereal = ereal + ereal2
!
                    eregion2region(nregionj,nregioni) = eregion2region(nregionj,nregioni) + esum
!
                    siteenergy(i) = siteenergy(i) + 0.5_dp*esum
                    siteenergy(j) = siteenergy(j) + 0.5_dp*esum
!
                    if (lPrintTwo) then
                      write(ioout,'(4x,i12,i12,f22.10,1x,f22.10)') i,j,eatom2,ereal2
                    endif
!
                    if (lDoQDeriv1.and.lgrad1) then
                      d0i = d0i + derive0*qlj
                      d0j = d0j + derive0*qli
                    endif
                    if (leem) then
                      if (lqeq.or.lSandM) then
                        eqeq2 = 0.0_dp
                        if (lqeq) then
                          call qeqbody(eqeq2,lgrad1,.false.,1_i4,1_i4,fct,qli,qlj,nati,natj,nqri,nqrj)
                        elseif (lSandM) then
                          call smbody(eqeq2,lgrad1,.false.,1_i4,1_i4,fct,qli,qlj,nati,natj,nqri,nqrj)
                        endif
                        eqeq = eqeq + eqeq2
!
                        eregion2region(nregionj,nregioni) = eregion2region(nregionj,nregioni) + eqeq2
!
                        siteenergy(i) = siteenergy(i) + 0.5_dp*eqeq2
                        siteenergy(j) = siteenergy(j) + 0.5_dp*eqeq2
                      endif
                    endif
                    if (lsuttonc) then
                      if (.not.lMEAMden) then
                        if (lorder12loc) then
                          scrho(1,i) = scrho(1,i) + sctrm1*ocj
                          scrho(1,j) = scrho(1,j) + sctrm2*oci
                        else
                          scrho(1,i) = scrho(1,i) + sctrm2*ocj
                          scrho(1,j) = scrho(1,j) + sctrm1*oci
                        endif
                      endif
                    endif
!*****************************
!  Charge first derivatives  *
!*****************************
                    if (lgrad1.and.lDoQDeriv1) then
                      call d1charge(i,j,lopi,lopj,1_i4,d0i,d0j)
                    endif
!*************************************
!  Electrostatic potential on-sites  *
!*************************************
                    if (lpolar) then
                      vx(i) = vx(i) - qlj*derivqd(1)*xji
                      vy(i) = vy(i) - qlj*derivqd(1)*yji
                      vz(i) = vz(i) - qlj*derivqd(1)*zji
                      vx(j) = vx(j) + qli*derivqd(1)*xji
                      vy(j) = vy(j) + qli*derivqd(1)*yji
                      vz(j) = vz(j) + qli*derivqd(1)*zji
                    endif
!***********************
!  Radial derivatives  *
!***********************
                    if (lgrad1) then
                      if (radi.gt.0.0_dp.and.lopi) then
                        raderv(i) = raderv(i) + rtrm1
                      endif
                      if (radj.gt.0.0_dp.and.lopj) then
                        raderv(j) = raderv(j) + rtrm1
                      endif
                    endif
!************************
!  Internal Derivatives *
!************************
!
!  First derivatives
!
                    if (lgrad1) then
                      if (lopi) then
                        xdrv(i) = xdrv(i) - deriv*xji
                        ydrv(i) = ydrv(i) - deriv*yji
                        zdrv(i) = zdrv(i) - deriv*zji
                      endif
                      if (lopj) then
                        xdrv(j) = xdrv(j) + deriv*xji
                        ydrv(j) = ydrv(j) + deriv*yji
                        zdrv(j) = zdrv(j) + deriv*zji
                      endif
                      if (nregioni.ne.nregionj) then
                        xregdrv(nregioni) = xregdrv(nregioni) - deriv*xji
                        yregdrv(nregioni) = yregdrv(nregioni) - deriv*yji
                        zregdrv(nregioni) = zregdrv(nregioni) - deriv*zji
                        xregdrv(nregionj) = xregdrv(nregionj) + deriv*xji
                        yregdrv(nregionj) = yregdrv(nregionj) + deriv*yji
                        zregdrv(nregionj) = zregdrv(nregionj) + deriv*zji
                      endif
                    endif
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
!*********************
!  Radial self term  *
!*********************
          if (lbsmat(nsft+nrelf2a(i)).and.lopi) then
!
!  Find self term
!
            eatom2 = 0.0_dp
            do n = 1,npote
              if (nptype(n).eq.14) then
                if (lmatch(nati,ntypi,nspec1(n),nptyp1(n),.true.)) then
                  rdiff = radi - twopot(2,n)
                  apt = twopot(1,n)*oci
                  eatom2 = eatom2 + 0.5_dp*apt*rdiff*rdiff
                  if (lgrad1) then
                    raderv(i) = raderv(i) + apt*rdiff
                  endif
                endif
              elseif (nptype(n).eq.17) then
                if (lmatch(nati,ntypi,nspec1(n),nptyp1(n),.true.)) then
                  rdiff = radi - twopot(3,n)
                  apt = twopot(1,n)*oci
                  bpt = twopot(2,n)
                  etrm1 = exp(bpt*rdiff)
                  etrm2 = 1.0_dp/etrm1
                  etrm = apt*(etrm1 + etrm2)
                  eatom2 = eatom2 + etrm
                  if (lgrad1) then
                    raderv(i) = raderv(i) + apt*bpt*(etrm1 - etrm2)
                  endif
                endif
              elseif (nptype(n).eq.31) then
                if (lmatch(nati,ntypi,nspec1(n),nptyp1(n),.true.)) then
                  rdiff = radi - twopot(3,n)
                  apt = twopot(1,n)*oci
                  bpt = twopot(2,n)
                  etrm1 = exp(bpt*rdiff)
                  etrm = apt*etrm1
                  eatom2 = eatom2 + etrm
                  if (lgrad1) then
                    raderv(i) = raderv(i) + apt*bpt*etrm1
                  endif
                endif
              endif
            enddo
!
            eregion2region(nregioni,nregioni) = eregion2region(nregioni,nregioni) + eatom2
!
            siteenergy(i) = siteenergy(i) + eatom2
            eatom = eatom + eatom2
          endif
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
      if (.not.lMEAMden) then
        do i = 1,numat
          sum2(i) = scrho(1,i) 
        enddo
        call sumall(sum2,sum,numat,"realmd0s","scrho")
        do i = 1,numat
          scrho(1,i) = sum(i)
        enddo
        do i = 1,numat
          sum2(i) = scrho12(1,i) 
        enddo
        call sumall(sum2,sum,numat,"realmd0s","scrho12")
        do i = 1,numat
          scrho12(1,i) = sum(i)
        enddo
      endif
    endif
    tsuml = g_cpu_time() - tsum0
    tsum = tsum + tsuml
  endif
!
!  Closing banner for energy decomposition
!
  if (lPrintTwo) then
    call mpbarrier
    if (ioproc) then
      write(ioout,'(''--------------------------------------------------------------------------------'')')
    endif
  endif
!
!  Free local memory
!
  deallocate(sum2,stat=status)
  if (status/=0) call deallocate_error('realmd0s','sum2')
  deallocate(sum,stat=status)
  if (status/=0) call deallocate_error('realmd0s','sum')
  deallocate(npotl,stat=status)
  if (status/=0) call deallocate_error('realmd0s','npotl')
!
!  Timing
!
  time2 = g_cpu_time()
  tatom = tatom + time2 - time1 - tsuml
#ifdef TRACE
  call trace_out('realmd0s')
#endif
!
  return
  end
