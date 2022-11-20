  subroutine pirealrecips
!
!  Calculates the derivatives of the polarisation energy.
!  New version re-created from pirealrecip after addition of 
!  second derivatives in non-symmetry case.
!
!   7/19 Created from pirealrecip
!   7/19 Parallelisation added
!   8/19 Short range damping of polarisation added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!  10/19 Langevin damping added
!   2/20 Correction to Langevin damped formalism
!   2/20 Modified for rigid molecules
!   4/20 d2r2dsdc added for benefit of rigid molecules
!   4/20 d2xyzdsdc added to cartstrterm arguments
!   4/20 derv3c changes reversed as they are no longer required
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
!  Julian Gale, CIC, Curtin University, April 2020
!
  use g_constants
  use control
  use current
  use datatypes
  use derivatives
  use element,        only : maxele
  use general
  use kspace
  use m_strain,       only : twostrterms, cartstrterm
  use molecule
  use optimisation
  use parallel
  use polarise
  use potentialxyz
  use realvectors
  use shells
  use symmetry
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  use two
  implicit none
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ifail
  integer(i4)                                  :: ii
  integer(i4)                                  :: indm
  integer(i4)                                  :: indmj
  integer(i4)                                  :: ix
  integer(i4)                                  :: iy
  integer(i4)                                  :: iz
  integer(i4)                                  :: ixi
  integer(i4)                                  :: ixj
  integer(i4)                                  :: iyi
  integer(i4)                                  :: iyj
  integer(i4)                                  :: izi
  integer(i4)                                  :: izj
  integer(i4)                                  :: j
  integer(i4)                                  :: jx
  integer(i4)                                  :: jy
  integer(i4)                                  :: jz
  integer(i4)                                  :: k
  integer(i4)                                  :: kk
  integer(i4)                                  :: kl
  integer(i4)                                  :: ks
  integer(i4)                                  :: n
  integer(i4)                                  :: nat1
  integer(i4)                                  :: nat2
  integer(i4)                                  :: nati
  integer(i4)                                  :: natj
  integer(i4)                                  :: nmi
  integer(i4)                                  :: nmj
  integer(i4)                                  :: nmolonly
  integer(i4)                                  :: nor
  integer(i4), dimension(:), allocatable       :: npotl
  integer(i4)                                  :: npots
  integer(i4)                                  :: nr
  integer(i4)                                  :: nreli
  integer(i4)                                  :: nrop
  integer(i4)                                  :: ntyp1
  integer(i4)                                  :: ntyp2
  integer(i4)                                  :: ntypi
  integer(i4)                                  :: ntypj
  integer(i4)                                  :: status
  logical                                      :: lcspair
  logical                                      :: lewaldtype
  logical                                      :: lmolok
  logical                                      :: lneedmol
  logical                                      :: lopi
  logical                                      :: lopj
  logical                                      :: lself
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: cut2
  real(dp)                                     :: cut2e
  real(dp)                                     :: cut2p
  real(dp)                                     :: cut2q
  real(dp)                                     :: cut2r
  real(dp)                                     :: cut2s
  real(dp)                                     :: cut2w
  real(dp)                                     :: d2(3,3)
  real(dp)                                     :: d2k(3,3)
  real(dp)                                     :: d2s(3,6)
  real(dp)                                     :: d3(3,3,3)
  real(dp)                                     :: d3s(3,6,3)
  real(dp)                                     :: d3ss(3,6,6)
  real(dp)                                     :: dxyzds(6,3)
  real(dp)                                     :: d2xyzdsdx(6,3,3)
  real(dp)                                     :: d2xyzds2(6,6,3)
  real(dp)                                     :: eatom
  real(dp)                                     :: ec6
  real(dp)                                     :: ereal
  real(dp)                                     :: Ep
  real(dp)                                     :: dEpdEs
  real(dp)                                     :: d2EpdEs2
  real(dp)                                     :: Es
  real(dp)                                     :: factor
  real(dp)                                     :: fcti
  real(dp)                                     :: oci
  real(dp)                                     :: ocj
  real(dp)                                     :: ofct
  real(dp)                                     :: ospfct
  real(dp)                                     :: poli
  real(dp)                                     :: polj
  real(dp)                                     :: qli
  real(dp)                                     :: qlj
  real(dp)                                     :: qpol
  real(dp)                                     :: rp
  real(dp)                                     :: rstrdloc(6)
  real(dp)                                     :: rtmp(6)
  real(dp)                                     :: rvp(3,3)
  real(dp)                                     :: rvpi(3,3)
  real(dp)                                     :: sctrm1
  real(dp)                                     :: sctrm2
  real(dp)                                     :: time1
  real(dp)                                     :: time2
  real(dp)                                     :: vxi 
  real(dp)                                     :: vyi
  real(dp)                                     :: vzi
  real(dp)                                     :: vxj
  real(dp)                                     :: vyj
  real(dp)                                     :: vzj
  real(dp)                                     :: vxji
  real(dp)                                     :: vyji
  real(dp)                                     :: vzji
  real(dp)                                     :: xal
  real(dp)                                     :: yal
  real(dp)                                     :: zal
  real(dp)                                     :: xcom
  real(dp)                                     :: ycom
  real(dp)                                     :: zcom
  real(dp)                                     :: xcomi
  real(dp)                                     :: ycomi
  real(dp)                                     :: zcomi
  real(dp)                                     :: xcrd
  real(dp)                                     :: ycrd
  real(dp)                                     :: zcrd
  real(dp)                                     :: xdrvij
  real(dp)                                     :: ydrvij
  real(dp)                                     :: zdrvij
#ifdef TRACE
  call trace_in('pirealrecips')
#endif
!
  time1 = g_cpu_time()
!
!  Zero energies although not needed to avoid overflow
!
  eatom = 0.0_dp
  ereal = 0.0_dp
  ec6 = 0.0_dp
!***************************
!  Set up local variables  *
!***************************
!
!  Set the Coulomb term type based on dimensionality :
!
!  1-D         => 1/r
!  2-D and 3-D => erfc(seta*r)/r
!
  lewaldtype = (ndim.ne.1.or.lwolf)
!
!  Set up cutoffs
!
  cut2p = cutp*cutp
  cut2e = rmx2
  cut2s = cuts*cuts
  cut2w = cutw*cutw
  if (lwolf) then        
    cut2q = cut2w        
  else
    cut2q = cut2e
  endif
!
!  Generate non-primitive cell
!
  do i = 1,3
    rvp(1,i) = rv(1,i)
    rvp(2,i) = rv(2,i)
    rvp(3,i) = rv(3,i)
  enddo
  if (ncbl.gt.1) call uncentre(rvp)
!
!  Find inverse of rvp
!
  rvpi(1:3,1:3) = rvp(1:3,1:3)
  call matrix_inversion(rvpi,3_i4,3_i4,rtmp,ifail)
!
!  Allocate local memory
!
  allocate(npotl(npote),stat=status)
  if (status/=0) call outofmemory('pirealrecips','npotl')
!
!  Initialise K vector terms
!
  call setktrmdp
!********************************************************
!  Atomistic and real space electrostatic derivatives   *
!********************************************************
!
!  Outer loop over sites with polarisability on site i
!
  ix = -2
  iy = -1
  iz =  0
!
  iloop: do i = 1,nasym
    lopi = (.not.lfreeze.or.lopf(i))
    if (lopi) then
      ix = ix + 3
      iy = iy + 3
      iz = iz + 3
    endif
!
!  If i has no polarisability or is frozen then skip this atom
!
    if (dpolar(i).eq.0.0_dp.or..not.lopi) cycle iloop
!
!  Inner loop over second site
!
    xal = xalat(i)
    yal = yalat(i)
    zal = zalat(i)
    nati = iatn(i)
    ntypi = natype(i)
    qli = qa(i)
    oci = occua(i)
    fcti = dble(neqv(i))
    nreli = nrela2f(i)
!
!  Molecule handling
!
    if (lmol) then
      nmi = natmol(nreli)
      indm = nmolind(nreli)
      call mindtoijk(indm,ixi,iyi,izi)
!
!  Set COM coordinates
!
      if (lrigid.and.nmi.gt.0) then
        xcomi = molxyz(1,natinmol(nreli),nmi)
        ycomi = molxyz(2,natinmol(nreli),nmi)
        zcomi = molxyz(3,natinmol(nreli),nmi)
      else
        xcomi = 0.0_dp
        ycomi = 0.0_dp
        zcomi = 0.0_dp
      endif
    else
      xcomi = 0.0_dp
      ycomi = 0.0_dp
      zcomi = 0.0_dp
    endif
!
!  Start of second atom loop
!
    jx = -2
    jy = -1
    jz =  0
    jloop: do j = procid+1,numat,nprocs
      lopj = (.not.lfreeze.or.lopf(nrelf2a(j)))
      if (lopj) then
        jx = jx + 3
        jy = jy + 3
        jz = jz + 3
      endif
      natj = nat(j)
      ntypj = nftype(j)
!
!  Arrange nat1 to be lower than nat2. If they are equal then arrange that ntyp1 is lower than ntyp2
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
      qlj = qf(j)
      ocj = occuf(j)
      ofct = oci*ocj
      factor = qli*qlj*ofct*angstoev
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
!  Molecule handling
!
      if (lmol) then
        nmj = natmol(j)
        indmj = nmolind(j)
        call mindtoijk(indmj,ixj,iyj,izj)
        ixj = ixj - ixi
        iyj = iyj - iyi
        izj = izj - izi
        lmolok = (nmi.eq.nmj.and.nmi.gt.0)
!
!  Set COM coordinates
!
        if (lrigid.and.nmj.gt.0) then
          xcom = molxyz(1,natinmol(j),nmj) - xcomi
          ycom = molxyz(2,natinmol(j),nmj) - ycomi
          zcom = molxyz(3,natinmol(j),nmj) - zcomi
        else
          xcom = - xcomi
          ycom = - ycomi
          zcom = - zcomi
        endif
      else
        lmolok = .false.
        xcom = - xcomi
        ycom = - ycomi
        zcom = - zcomi
      endif
!
!  Locate potential number
!  Check whether potential requires specific types
!
      rp = 0.0_dp
      npots = 0
      lneedmol = (lmol.and..not.lmolq)
      do n = 1,npote
        if (nat1.eq.nspec1(n).and.nat2.eq.nspec2(n)) then
          if ((ntyp1.eq.nptyp1(n).or.nptyp1(n).eq.0).and.(ntyp2.eq.nptyp2(n).or.nptyp2(n).eq.0)) then
            npots = npots + 1
            npotl(npots) = n
            if (mmexc(n).gt.0.or.(lintra(n).and..not.linter(n).or..not.lintra(n).and.linter(n))) lneedmol = .true.
            if (nptype(n).eq.8.or.nptype(n).eq.33) then
              if (cuts.gt.rp) rp = cuts
            else
              if (rpot(n).gt.rp) rp = rpot(n)
            endif
          endif
        endif
      enddo
      cut2r = rp*rp
      if (cut2r.gt.cut2p) cut2r = cut2p
      cut2 = cut2r
      if (cut2e.gt.cut2.and.lewald) cut2 = cut2e
      if (cut2w.gt.cut2.and.lwolf) cut2 = cut2w
!
!  If the Coulomb term and the potentials don't need molecule info for
!  for this pair of atoms, then turn off molecule checking to save time
!
      if (.not.lneedmol) lmolok = .false.
!
      nmolonly = 0
      nor = 0
!
!  Zero second derivative arrays
!
      d2(1:3,1:3) = 0.0_dp
      if (lstr) then
        d2s(1:3,1:nstrains) = 0.0_dp
      endif
!
!  If no valid potentials and charge product/polarisability is zero then no need to search for distances
!
      qpol = abs(qlj)*abs(dpolar(i))
      if (npots.eq.0.and.abs(qpol).lt.1.0d-8) cycle jloop
!***********************
!  Find valid vectors  *
!***********************
      if (ndim.eq.3) then
        call rsearch3D(xcrd,ycrd,zcrd,lmolok,lcspair,i,j,nmi,ixj,iyj,izj,nor,nmolonly,lself,cut2)
      elseif (ndim.eq.2) then
        call rsearch2D(xcrd,ycrd,zcrd,lmolok,lcspair,i,j,nmi,ixj,iyj,izj,nor,nmolonly,lself,cut2)
      elseif (ndim.eq.1) then
        call rsearch1D(xcrd,ycrd,zcrd,lmolok,lcspair,i,j,nmi,ixj,iyj,izj,nor,nmolonly,lself,cut2)
      endif
      if (nor.eq.0) cycle jloop
      if (.not.lnoreal) then
!
!  Sqrt distances
!
        do k = 1,nor
          dist(k) = sqrt(dist(k))
        enddo
!
        call twobody(eatom,ereal,ec6,.true.,.true.,.false.,nor,1,npots,npotl,cut2,cut2q,cut2s, &
                     nmolonly,factor,ofct,ospfct,0.0_dp,sctrm1,sctrm2,qli,qlj,lcspair,lewaldtype, &
                     .false.,.false.,.false.)
      endif
!*****************************************************************
!  Calculate reciprocal space contribution to third derivatives  *
!*****************************************************************
      if (lewald.and..not.lnorecip) then
        call reciptrmdp(xcrd,ycrd,zcrd,.false.,lstr,ofct,d2,d2s,d3,d3s,d3ss)
      endif
!****************************
!  Loop over all distances  *
!****************************
      if (lpollangevin) then
        Es = sqrt(vx(i)**2 + vy(i)**2 + vz(i)**2)
        if (Es.gt.1.0d-8) then
          call langevinpol(Es,dpolarmax(i),dpolar(i),Ep,dEpdEs,d2EpdEs2,.true.,.false.)
          poli = qlj*fcti*dEpdEs/angstoev
        else
          poli = qlj*dpolar(i)*fcti/angstoev
        endif
      else
        poli = qlj*dpolar(i)*fcti/angstoev
      endif
      if (.not.lnoreal) then
!
!  Generate products for derivatives
!
        call twostrterms(ndim,maxdis,nor,xtmp,ytmp,ztmp,xcom,ycom,zcom,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.true.)
!
        do kk = 1,nor
          xcrd = xtmp(kk)
          ycrd = ytmp(kk)
          zcrd = ztmp(kk)
!
!  Second derivative terms for first derivatives
!
          d2k(1,1) = derivqd2(kk)*d2r2dx2(kk,1) + derivqd(kk)
          d2k(2,1) = derivqd2(kk)*d2r2dx2(kk,6)
          d2k(3,1) = derivqd2(kk)*d2r2dx2(kk,5)
          d2k(1,2) = derivqd2(kk)*d2r2dx2(kk,6)
          d2k(2,2) = derivqd2(kk)*d2r2dx2(kk,2) + derivqd(kk)
          d2k(3,2) = derivqd2(kk)*d2r2dx2(kk,4)
          d2k(1,3) = derivqd2(kk)*d2r2dx2(kk,5)
          d2k(2,3) = derivqd2(kk)*d2r2dx2(kk,4)
          d2k(3,3) = derivqd2(kk)*d2r2dx2(kk,3) + derivqd(kk)
!
!  Add to cumulative arrays
!
          do ii = 1,3
            d2(1,ii) = d2(1,ii) + d2k(1,ii)
            d2(2,ii) = d2(2,ii) + d2k(2,ii)
            d2(3,ii) = d2(3,ii) + d2k(3,ii)
          enddo
!
!  Strain derivatives
!
          if (lstr) then
            call cartstrterm(ndim,xcrd,ycrd,zcrd,xcom,ycom,zcom,dxyzds,d2xyzdsdx,d2xyzds2,.false.)
            do kl = 1,nstrains
              ks = nstrptr(kl)
!
!  Note we need to use Cartesian strain derivatives rather than d2r2dsdx
!  since the terms have to be calculated consistently in the zero strain limit
!
              d2s(1,kl) = d2s(1,kl) - derivqd(kk)*dxyzds(kl,1) - xcrd*derivqd2(kk)*dr2ds(kk,ks)
              d2s(2,kl) = d2s(2,kl) - derivqd(kk)*dxyzds(kl,2) - ycrd*derivqd2(kk)*dr2ds(kk,ks)
              d2s(3,kl) = d2s(3,kl) - derivqd(kk)*dxyzds(kl,3) - zcrd*derivqd2(kk)*dr2ds(kk,ks)
            enddo
          endif
!****************************
!  End loop over distances  *
!****************************
        enddo
      endif
!
      vxi = vx(i)
      vyi = vy(i)
      vzi = vz(i)
!**************************
!  Coordinate Derivatives *
!**************************
!
!  First derivatives
!
      xdrvij = poli*(vxi*d2(1,1) + vyi*d2(2,1) + vzi*d2(3,1))
      ydrvij = poli*(vxi*d2(1,2) + vyi*d2(2,2) + vzi*d2(3,2))
      zdrvij = poli*(vxi*d2(1,3) + vyi*d2(2,3) + vzi*d2(3,3))
!
      if (j.ne.nreli) then
        xdrv(i) = xdrv(i) - xdrvij
        ydrv(i) = ydrv(i) - ydrvij
        zdrv(i) = zdrv(i) - zdrvij
      endif
!
      if (lstr) then
        do ks = 1,nstrains
          rstrdloc(ks) = - poli*(vxi*d2s(1,ks) + vyi*d2s(2,ks) + vzi*d2s(3,ks))
          rstrd(ks) = rstrd(ks) + rstrdloc(ks)
        enddo
      endif
    enddo jloop
  enddo iloop
!
!  Outer loop over sites with charge acting on site i acting on a polarisable atom
!
  ix = -2
  iy = -1
  iz =  0
!
  iloop2: do i = 1,nasym
    lopi = (.not.lfreeze.or.lopf(i))
    if (lopi) then
      ix = ix + 3
      iy = iy + 3
      iz = iz + 3
    endif
!
!  If i has no charge or is frozen then skip this atom
!
    qli = qa(i)
    if (qli.eq.0.0_dp.or..not.lopi) cycle iloop2
!
!  Inner loop over second site
!
    xal = xalat(i)
    yal = yalat(i)
    zal = zalat(i)
    nati = iatn(i)
    ntypi = natype(i)
    oci = occua(i)
    fcti = dble(neqv(i))
    nreli = nrela2f(i)
!
!  Molecule handling
!
    if (lmol) then
      nmi = natmol(nreli)
      indm = nmolind(nreli)
      call mindtoijk(indm,ixi,iyi,izi)
!
!  Set COM coordinates
!
      if (lrigid.and.nmi.gt.0) then
        xcomi = molxyz(1,natinmol(nreli),nmi)
        ycomi = molxyz(2,natinmol(nreli),nmi)
        zcomi = molxyz(3,natinmol(nreli),nmi)
      else
        xcomi = 0.0_dp
        ycomi = 0.0_dp
        zcomi = 0.0_dp
      endif
    else
      xcomi = 0.0_dp
      ycomi = 0.0_dp
      zcomi = 0.0_dp
    endif
!
!  Start of second atom loop
!
    jx = -2
    jy = -1
    jz =  0
    jloop3: do j = procid+1,numat,nprocs
      lopj = (.not.lfreeze.or.lopf(nrelf2a(j)))
      if (lopj) then
        jx = jx + 3
        jy = jy + 3
        jz = jz + 3
      endif
!
!  If j has no polarisability then skip
!
      if (dpolar(nrelf2a(j)).eq.0.0_dp) cycle jloop3
!
      natj = nat(j)
      ntypj = nftype(j)
!
!  Arrange nat1 to be lower than nat2. If they are equal then arrange that ntyp1 is lower than ntyp2
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
      qlj = qf(j)
      ocj = occuf(j)
      ofct = oci*ocj
      factor = qli*qlj*ofct*angstoev
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
!  Molecule handling
!
      if (lmol) then
        nmj = natmol(j)
        indmj = nmolind(j)
        call mindtoijk(indmj,ixj,iyj,izj)
        ixj = ixj - ixi
        iyj = iyj - iyi
        izj = izj - izi
        lmolok = (nmi.eq.nmj.and.nmi.gt.0)
!
!  Set COM coordinates
!
        if (lrigid.and.nmj.gt.0) then
          xcom = molxyz(1,natinmol(j),nmj) - xcomi
          ycom = molxyz(2,natinmol(j),nmj) - ycomi
          zcom = molxyz(3,natinmol(j),nmj) - zcomi
        else
          xcom = - xcomi
          ycom = - ycomi
          zcom = - zcomi
        endif
      else
        lmolok = .false.
        xcom = - xcomi
        ycom = - ycomi
        zcom = - zcomi
      endif
!
!  Locate potential number
!  Check whether potential requires specific types
!
      rp = 0.0_dp
      npots = 0
      lneedmol = (lmol.and..not.lmolq)
      do n = 1,npote
        if (nat1.eq.nspec1(n).and.nat2.eq.nspec2(n)) then
          if ((ntyp1.eq.nptyp1(n).or.nptyp1(n).eq.0).and.(ntyp2.eq.nptyp2(n).or.nptyp2(n).eq.0)) then
            npots = npots + 1
            npotl(npots) = n
            if (mmexc(n).gt.0.or.(lintra(n).and..not.linter(n).or..not.lintra(n).and.linter(n))) lneedmol = .true.
            if (nptype(n).eq.8.or.nptype(n).eq.33) then
              if (cuts.gt.rp) rp = cuts
            else
              if (rpot(n).gt.rp) rp = rpot(n)
            endif
          endif
        endif
      enddo
      cut2r = rp*rp
      if (cut2r.gt.cut2p) cut2r = cut2p
      cut2 = cut2r
      if (cut2e.gt.cut2.and.lewald) cut2 = cut2e
      if (cut2w.gt.cut2.and.lwolf) cut2 = cut2w
!
!  If the Coulomb term and the potentials don't need molecule info for
!  for this pair of atoms, then turn off molecule checking to save time
!
      if (.not.lneedmol) lmolok = .false.
!
      nmolonly = 0
      nor = 0
!
!  Zero second derivative arrays
!
      d2(1:3,1:3) = 0.0_dp
      if (lstr) then
        d2s(1:3,1:nstrains) = 0.0_dp
      endif
!
!  If no valid potentials and charge product/polarisability is zero then no need to search for distances
!
      qpol = abs(qli)*abs(dpolar(nrelf2a(j)))
      if (npots.eq.0.and.abs(qpol).lt.1.0d-8) cycle jloop3
!***********************
!  Find valid vectors  *
!***********************
      if (ndim.eq.3) then
        call rsearch3D(xcrd,ycrd,zcrd,lmolok,lcspair,i,j,nmi,ixj,iyj,izj,nor,nmolonly,lself,cut2)
      elseif (ndim.eq.2) then
        call rsearch2D(xcrd,ycrd,zcrd,lmolok,lcspair,i,j,nmi,ixj,iyj,izj,nor,nmolonly,lself,cut2)
      elseif (ndim.eq.1) then
        call rsearch1D(xcrd,ycrd,zcrd,lmolok,lcspair,i,j,nmi,ixj,iyj,izj,nor,nmolonly,lself,cut2)
      endif
      if (nor.eq.0) cycle jloop3
      if (.not.lnoreal) then
!
!  Sqrt distances
!
        do k = 1,nor
          dist(k) = sqrt(dist(k))
        enddo
!
        call twobody(eatom,ereal,ec6,.true.,.true.,.false.,nor,1,npots,npotl,cut2,cut2q,cut2s, &
                     nmolonly,factor,ofct,ospfct,0.0_dp,sctrm1,sctrm2,qli,qlj,lcspair,lewaldtype, &
                     .false.,.false.,.false.)
      endif
!*****************************************************************
!  Calculate reciprocal space contribution to third derivatives  *
!*****************************************************************
      if (lewald.and..not.lnorecip) then
        call reciptrmdp(xcrd,ycrd,zcrd,.false.,lstr,ofct,d2,d2s,d3,d3s,d3ss)
      endif
!***************************************************************
!  Rotate potential from equivalent asymmetric unit atom to j  *
!***************************************************************
      nr = nrelf2a(j)
      nrop = nrotop(j)
!
!  Convert potential at i into fractional space
!
      vxi = vx(nr)*rvpi(1,1) + vy(nr)*rvpi(1,2) + vz(nr)*rvpi(1,3)
      vyi = vx(nr)*rvpi(2,1) + vy(nr)*rvpi(2,2) + vz(nr)*rvpi(2,3)
      vzi = vx(nr)*rvpi(3,1) + vy(nr)*rvpi(3,2) + vz(nr)*rvpi(3,3)
!
!  Perform rotation
!
      vxji = rop(1,1,nrop)*vxi + rop(1,2,nrop)*vyi + rop(1,3,nrop)*vzi
      vyji = rop(2,1,nrop)*vxi + rop(2,2,nrop)*vyi + rop(2,3,nrop)*vzi
      vzji = rop(3,1,nrop)*vxi + rop(3,2,nrop)*vyi + rop(3,3,nrop)*vzi
!
!  Convert potential back to Cartesian space
!
      vxj = rvp(1,1)*vxji + rvp(1,2)*vyji + rvp(1,3)*vzji
      vyj = rvp(2,1)*vxji + rvp(2,2)*vyji + rvp(2,3)*vzji
      vzj = rvp(3,1)*vxji + rvp(3,2)*vyji + rvp(3,3)*vzji
!****************************
!  Loop over all distances  *
!****************************
      if (lpollangevin) then
        Es = sqrt(vxj**2 + vyj**2 + vzj**2)
        if (Es.gt.1.0d-8) then
          call langevinpol(Es,dpolarmax(nr),dpolar(nr),Ep,dEpdEs,d2EpdEs2,.true.,.false.)
          polj = qli*dEpdEs*fcti/angstoev
        else
          polj = qli*dpolar(nr)*fcti/angstoev
        endif
      else
        polj = qli*dpolar(nr)*fcti/angstoev
      endif
      if (.not.lnoreal) then
!
!  Generate products for derivatives
!
        call twostrterms(ndim,maxdis,nor,xtmp,ytmp,ztmp,xcom,ycom,zcom,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.true.)
!
        do kk = 1,nor
          xcrd = xtmp(kk)
          ycrd = ytmp(kk)
          zcrd = ztmp(kk)
!
!  Second derivative terms for first derivatives
!
          d2k(1,1) = derivqd2(kk)*d2r2dx2(kk,1) + derivqd(kk)
          d2k(2,1) = derivqd2(kk)*d2r2dx2(kk,6)
          d2k(3,1) = derivqd2(kk)*d2r2dx2(kk,5)
          d2k(1,2) = derivqd2(kk)*d2r2dx2(kk,6)
          d2k(2,2) = derivqd2(kk)*d2r2dx2(kk,2) + derivqd(kk)
          d2k(3,2) = derivqd2(kk)*d2r2dx2(kk,4)
          d2k(1,3) = derivqd2(kk)*d2r2dx2(kk,5)
          d2k(2,3) = derivqd2(kk)*d2r2dx2(kk,4)
          d2k(3,3) = derivqd2(kk)*d2r2dx2(kk,3) + derivqd(kk)
!
!  Add to cumulative arrays
!
          do ii = 1,3
            d2(1,ii) = d2(1,ii) + d2k(1,ii)
            d2(2,ii) = d2(2,ii) + d2k(2,ii)
            d2(3,ii) = d2(3,ii) + d2k(3,ii)
          enddo
!
!  Strain derivatives
!
          if (lstr) then
            call cartstrterm(ndim,xcrd,ycrd,zcrd,xcom,ycom,zcom,dxyzds,d2xyzdsdx,d2xyzds2,.false.)
            do kl = 1,nstrains
              ks = nstrptr(kl)
!
!  Note we need to use Cartesian strain derivatives rather than d2r2dsdx
!  since the terms have to be calculated consistently in the zero strain limit
!
              d2s(1,kl) = d2s(1,kl) - derivqd(kk)*dxyzds(kl,1) - xcrd*derivqd2(kk)*dr2ds(kk,ks)
              d2s(2,kl) = d2s(2,kl) - derivqd(kk)*dxyzds(kl,2) - ycrd*derivqd2(kk)*dr2ds(kk,ks)
              d2s(3,kl) = d2s(3,kl) - derivqd(kk)*dxyzds(kl,3) - zcrd*derivqd2(kk)*dr2ds(kk,ks)
            enddo
          endif
!****************************
!  End loop over distances  *
!****************************
        enddo
      endif
!**************************
!  Coordinate Derivatives *
!**************************
      if (j.ne.nreli) then
!
!  First derivatives
!
        xdrvij = polj*(vxj*d2(1,1) + vyj*d2(2,1) + vzj*d2(3,1))
        ydrvij = polj*(vxj*d2(1,2) + vyj*d2(2,2) + vzj*d2(3,2))
        zdrvij = polj*(vxj*d2(1,3) + vyj*d2(2,3) + vzj*d2(3,3))
!
        xdrv(i) = xdrv(i) + xdrvij
        ydrv(i) = ydrv(i) + ydrvij
        zdrv(i) = zdrv(i) + zdrvij
      endif
    enddo jloop3
  enddo iloop2
!
!  Free local memory
!
  deallocate(npotl,stat=status)
  if (status/=0) call deallocate_error('pirealrecips','npotl')
!
!  Timing
!
  time2 = g_cpu_time()
  tpolar = tpolar + time2 - time1
#ifdef TRACE
  call trace_out('pirealrecips')
#endif
!
  return
  end
