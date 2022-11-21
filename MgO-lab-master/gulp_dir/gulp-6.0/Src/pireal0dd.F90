  subroutine pireal0dd(lgrad2)
!
!  Subroutine for calculating the polarisation contribution
!  to the first derivatives for a finite cluster => 0 D system
!  Distributed memory parallel version.
!
!   7/19 Created from pireal0d
!   8/19 Short range damping of polarisation added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!  10/19 Langevin damping added
!  11/19 Langevin second derivatives added
!   2/20 Correction to Langevin damped formalism
!   2/20 Modified for rigid molecules
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
!  Julian Gale, CIC, Curtin University, February 2020
!
  use configurations, only : nregionno
  use g_constants
  use control
  use current
  use derivatives
  use element,        only : maxele
  use general,        only : cutw
  use molecule
  use optimisation
  use parallel
  use polarise
  use potentialxyz
  use realvectors,    only : derivqd, derivqd2, derivqd3
  use shells
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  use two
  implicit none
!
!  Passed variables
!
  logical,                      intent(in)     :: lgrad2
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ii
  integer(i4)                                  :: ind
  integer(i4)                                  :: indmi
  integer(i4)                                  :: indmj
  integer(i4)                                  :: ix
  integer(i4)                                  :: ixx
  integer(i4)                                  :: iy
  integer(i4)                                  :: iyy
  integer(i4)                                  :: iz
  integer(i4)                                  :: izz
  integer(i4)                                  :: j
  integer(i4)                                  :: jonnode
  integer(i4)                                  :: jx
  integer(i4)                                  :: jxx
  integer(i4)                                  :: jy
  integer(i4)                                  :: jyy
  integer(i4)                                  :: jz
  integer(i4)                                  :: jzz
  integer(i4)                                  :: k
  integer(i4)                                  :: kx
  integer(i4)                                  :: ky
  integer(i4)                                  :: kz
  integer(i4)                                  :: n
  integer(i4)                                  :: nat1
  integer(i4)                                  :: nat2
  integer(i4)                                  :: nati
  integer(i4)                                  :: natj
  integer(i4)                                  :: nbtypeij
  integer(i4)                                  :: nbtypeij2
  integer(i4)                                  :: nmi
  integer(i4)                                  :: nmj
  integer(i4)                                  :: nmolonly
  integer(i4)                                  :: nor
  integer(i4), dimension(:), allocatable       :: npotl
  integer(i4)                                  :: npots
  integer(i4)                                  :: nregioni
  integer(i4)                                  :: nregionj
  integer(i4)                                  :: ntyp1
  integer(i4)                                  :: ntyp2
  integer(i4)                                  :: ntypi
  integer(i4)                                  :: ntypj
  integer(i4)                                  :: status
  logical                                      :: l2bonds
  logical                                      :: l3bonds
  logical                                      :: lbonded
  logical                                      :: lcspair
  logical                                      :: lmolok
  logical                                      :: lneedmol
  logical                                      :: lptrmol
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: cut2p
  real(dp)                                     :: cut2q
  real(dp)                                     :: cut2r
  real(dp)                                     :: cut2s
  real(dp)                                     :: d0i
  real(dp)                                     :: d0j
  real(dp)                                     :: d1iq
  real(dp)                                     :: d1jq
  real(dp)                                     :: d2i2q
  real(dp)                                     :: d2ijq
  real(dp)                                     :: d2j2q
  real(dp)                                     :: d2loc(3,3)
  real(dp)                                     :: d3loc(3,3,3)
  real(dp)                                     :: deriv
  real(dp)                                     :: deriv2
  real(dp)                                     :: deriv3
  real(dp)                                     :: derive0
  real(dp)                                     :: derive
  real(dp)                                     :: derive2
  real(dp)                                     :: derive3
  real(dp)                                     :: dist
  real(dp)                                     :: eatom
  real(dp)                                     :: ec6
  real(dp)                                     :: ereal
  real(dp)                                     :: Ep
  real(dp)                                     :: dEpdEs
  real(dp)                                     :: d2EpdEs2
  real(dp)                                     :: Es
  real(dp)                                     :: factor
  real(dp)                                     :: fct
  real(dp)                                     :: oci
  real(dp)                                     :: ocj
  real(dp)                                     :: ofct
  real(dp)                                     :: ospfct
  real(dp)                                     :: poli
  real(dp)                                     :: pol2i
  real(dp)                                     :: polj
  real(dp)                                     :: pol2j
  real(dp)                                     :: polk
  real(dp)                                     :: pol2k
  real(dp)                                     :: qli
  real(dp)                                     :: qlj
  real(dp)                                     :: qlk
  real(dp)                                     :: r
  real(dp)                                     :: rderiv
  real(dp)                                     :: rp
  real(dp)                                     :: rpd1
  real(dp)                                     :: rpd2
  real(dp)                                     :: rpd3
  real(dp)                                     :: rpd4
  real(dp)                                     :: rpd5
  real(dp)                                     :: rpd6
  real(dp)                                     :: rtrm1
  real(dp)                                     :: rtrm2
  real(dp)                                     :: rtrm3
  real(dp)                                     :: rtrm32
  real(dp)                                     :: sctrm1
  real(dp)                                     :: sctrm2
  real(dp)                                     :: small
  real(dp)                                     :: small2
  real(dp)                                     :: time1
  real(dp)                                     :: time2
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
  real(dp),    dimension(:,:,:), allocatable   :: d2k
  real(dp),    dimension(:,:,:), allocatable   :: d2ksum
#ifdef TRACE
  call trace_in('pireal0dd')
#endif
!
  time1 = g_cpu_time()
!
!  Local variables
!
  small = 1.0d-12
  small2 = 1.0d-2
!
!  Set up cutoffs
!
  cut2p = cutp*cutp
  cut2s = cuts*cuts
  if (lwolf) then        
    cut2q = cutw*cutw
  else
    cut2q = 1.0d12
  endif
!
  if (lnoreal) goto 999
!
!  Allocate local memory
!
  allocate(npotl(npote),stat=status)
  if (status/=0) call outofmemory('pireal0dd','npotl')
  if (lgrad2) then
    allocate(d2k(3,3,numat),stat=status)
    if (status/=0) call outofmemory('pireal0dd','d2k')
    allocate(d2ksum(3,3,numat),stat=status)
    if (status/=0) call outofmemory('pireal0dd','d2ksum')
  endif
!
!  Outer loop over sites with polarisability on site i
!
  ix = -2
  iy = -1
  iz =  0
!
!  Loop over atoms on local node
!
  iloop: do ii = 1,natomsonnode
    i = node2atom(ii)
    ix = ix + 3
    iy = iy + 3
    iz = iz + 3
!
!  If i has no polarisability or charge then skip this atom
!
    qli = qf(i)
    if (dpolar(i).eq.0.0_dp.and.qli.eq.0.0_dp) cycle iloop
!
    xal = xclat(i)
    yal = yclat(i)
    zal = zclat(i)
    nati = nat(i)
    ntypi = nftype(i)
    nregioni = nregionno(nsft+nrelf2a(i))
    oci = occuf(i)
!
!  Molecule handling
!
    if (lmol) then
      nmi = natmol(i)
      indmi = nmolind(i)
!
!  Set COM coordinates
!
      if (lrigid.and.nmi.gt.0) then
        xcomi = molxyz(1,natinmol(i),nmi)
        ycomi = molxyz(2,natinmol(i),nmi)
        zcomi = molxyz(3,natinmol(i),nmi)
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
!  Inner loop over second site
!
    jx = -2
    jy = -1
    jz =  0
    jloop: do j = 1,numat
      jx = jx + 3
      jy = jy + 3
      jz = jz + 3
!
!  Exclude i = j in 0-D
!
      if (j.eq.i) cycle jloop
!
      qlj = qf(j)
!
!  Compute charge x polarisability for pair
!
      if (lpollangevin) then
        Es = sqrt(vx(i)**2 + vy(i)**2 + vz(i)**2)
        if (Es.gt.1.0d-8) then
          call langevinpol(Es,dpolarmax(i),dpolar(i),Ep,dEpdEs,d2EpdEs2,.true.,lgrad2)
          poli = qlj*dEpdEs/angstoev
          if (lgrad2) then
            pol2i = qlj*qlj*d2EpdEs2/angstoev
          endif
        else
          poli = qlj*dpolar(i)/angstoev
          pol2i = 0.0_dp
        endif
        Es = sqrt(vx(j)**2 + vy(j)**2 + vz(j)**2)
        if (Es.gt.1.0d-8) then
          call langevinpol(Es,dpolarmax(j),dpolar(j),Ep,dEpdEs,d2EpdEs2,.true.,lgrad2)
          polj = qli*dEpdEs/angstoev
          if (lgrad2) then
            pol2j = qli*qli*d2EpdEs2/angstoev
          endif
        else
          polj = qli*dpolar(j)/angstoev
          pol2j = 0.0_dp
        endif
      else
        poli = qlj*dpolar(i)/angstoev
        polj = qli*dpolar(j)/angstoev
        pol2i = 0.0_dp
        pol2j = 0.0_dp
      endif
!
!  If there will be no interaction then skip
!
      if (abs(poli).lt.1.0d-12.and.abs(polj).lt.1.0d-12) cycle jloop
!
      natj = nat(j)
      ntypj = nftype(j)
      nregionj = nregionno(nsft+nrelf2a(j))
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
          ntyp1 = ntypi
          ntyp2 = ntypj
        else
          ntyp1 = ntypj
          ntyp2 = ntypi
        endif
      elseif (nati.lt.natj) then
        nat1 = nati
        nat2 = natj
        ntyp1 = ntypi
        ntyp2 = ntypj
      else
        nat1 = natj
        nat2 = nati
        ntyp1 = ntypj
        ntyp2 = ntypi
      endif
      ofct = oci*ocj
      fct = ofct*angstoev
      factor = qli*qlj*fct
!
!  Possible core - shell flag
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
      r = xcrd*xcrd + ycrd*ycrd + zcrd*zcrd
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
!
      if (lmolok.and.(r.gt.cut2s.or..not.lcspair)) then
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
          lbonded = .false.
          l2bonds = .false.
          l3bonds = .false.
        endif
      else
        lptrmol = .false.
        lbonded = .false.
        l2bonds = .false.
        l3bonds = .false.
      endif
      if (abs(r-small2).lt.1.0d-12) r = small2
      if (r.lt.small) then
        cycle jloop
      else
!
!  Store vector
!
        nor = 1
        dist = sqrt(r)
      endif
!*************************************************
!  Calculate twobody contribution in real space  *
!*************************************************
      eatom = 0.0_dp
      ereal = 0.0_dp
      ec6 = 0.0_dp
      call twobody1(eatom,ereal,ec6,.true.,.true.,lgrad2,nor,1_i4,dist,xcrd,ycrd,zcrd,d0i,d0j, &
                    deriv,deriv2,deriv3,derive0,derive,derive2,derive3,rderiv,npots,npotl,0.0_dp, &
                    cut2q,cut2s,lptrmol,nmolonly,factor,ofct,ospfct,0.0_dp,rtrm1,rtrm2,rtrm3,rtrm32, &
                    sctrm1,sctrm2,qli,qlj,lcspair,.false.,.false.,lbonded,l2bonds,l3bonds, &
                    nbtypeij,nbtypeij2,.false.,.false.,d1iq,d1jq,d2i2q,d2ijq,d2j2q)
!**************************
!  Coordinate Derivatives *
!**************************
!
!  First derivatives
!
      rpd1 = xcrd*xcrd
      rpd2 = ycrd*ycrd
      rpd3 = zcrd*zcrd
      rpd4 = ycrd*zcrd
      rpd5 = xcrd*zcrd
      rpd6 = xcrd*ycrd
!
      d2loc(1,1) = derivqd2(1)*rpd1 + derivqd(1)
      d2loc(2,1) = derivqd2(1)*rpd6
      d2loc(3,1) = derivqd2(1)*rpd5
      d2loc(1,2) = derivqd2(1)*rpd6
      d2loc(2,2) = derivqd2(1)*rpd2 + derivqd(1)
      d2loc(3,2) = derivqd2(1)*rpd4
      d2loc(1,3) = derivqd2(1)*rpd5
      d2loc(2,3) = derivqd2(1)*rpd4
      d2loc(3,3) = derivqd2(1)*rpd3 + derivqd(1)
!
      xdrv(i) = xdrv(i) - poli*(vx(i)*d2loc(1,1) + vy(i)*d2loc(2,1) + vz(i)*d2loc(3,1))
      ydrv(i) = ydrv(i) - poli*(vx(i)*d2loc(1,2) + vy(i)*d2loc(2,2) + vz(i)*d2loc(3,2))
      zdrv(i) = zdrv(i) - poli*(vx(i)*d2loc(1,3) + vy(i)*d2loc(2,3) + vz(i)*d2loc(3,3))
!
      xdrv(j) = xdrv(j) + poli*(vx(i)*d2loc(1,1) + vy(i)*d2loc(2,1) + vz(i)*d2loc(3,1))
      ydrv(j) = ydrv(j) + poli*(vx(i)*d2loc(1,2) + vy(i)*d2loc(2,2) + vz(i)*d2loc(3,2))
      zdrv(j) = zdrv(j) + poli*(vx(i)*d2loc(1,3) + vy(i)*d2loc(2,3) + vz(i)*d2loc(3,3))
!
      if (nregioni.ne.nregionj) then
        xregdrv(nregioni) = xregdrv(nregioni) - poli*(vx(i)*d2loc(1,1) + vy(i)*d2loc(2,1) + vz(i)*d2loc(3,1))
        yregdrv(nregioni) = yregdrv(nregioni) - poli*(vx(i)*d2loc(1,2) + vy(i)*d2loc(2,2) + vz(i)*d2loc(3,2))
        zregdrv(nregioni) = zregdrv(nregioni) - poli*(vx(i)*d2loc(1,3) + vy(i)*d2loc(2,3) + vz(i)*d2loc(3,3))
!
        xregdrv(nregionj) = xregdrv(nregionj) + poli*(vx(i)*d2loc(1,1) + vy(i)*d2loc(2,1) + vz(i)*d2loc(3,1))
        yregdrv(nregionj) = yregdrv(nregionj) + poli*(vx(i)*d2loc(1,2) + vy(i)*d2loc(2,2) + vz(i)*d2loc(3,2))
        zregdrv(nregionj) = zregdrv(nregionj) + poli*(vx(i)*d2loc(1,3) + vy(i)*d2loc(2,3) + vz(i)*d2loc(3,3))
      endif
      if (lgrad2) then
!
!  Save d2 term
!
        d2k(1:3,1:3,j) = d2loc(1:3,1:3)
!
!  Second derivatives
!
        d3loc(1,1,1) = derivqd3(1)*rpd1*xcrd + 3.0_dp*xcrd*derivqd2(1)
        d3loc(2,1,1) = derivqd3(1)*rpd1*ycrd + ycrd*derivqd2(1)
        d3loc(3,1,1) = derivqd3(1)*rpd1*zcrd + zcrd*derivqd2(1)
        d3loc(1,2,1) = derivqd3(1)*rpd1*ycrd + ycrd*derivqd2(1)
        d3loc(2,2,1) = derivqd3(1)*rpd2*xcrd + xcrd*derivqd2(1)
        d3loc(3,2,1) = derivqd3(1)*rpd4*xcrd
        d3loc(1,3,1) = derivqd3(1)*rpd1*zcrd + zcrd*derivqd2(1)
        d3loc(2,3,1) = derivqd3(1)*rpd4*xcrd
        d3loc(3,3,1) = derivqd3(1)*rpd3*xcrd + xcrd*derivqd2(1)
!
        d3loc(1,1,2) = derivqd3(1)*rpd1*ycrd + ycrd*derivqd2(1)
        d3loc(2,1,2) = derivqd3(1)*rpd2*xcrd + xcrd*derivqd2(1)
        d3loc(3,1,2) = derivqd3(1)*rpd5*ycrd
        d3loc(1,2,2) = derivqd3(1)*rpd2*xcrd + xcrd*derivqd2(1)
        d3loc(2,2,2) = derivqd3(1)*rpd2*ycrd + 3.0_dp*ycrd*derivqd2(1)
        d3loc(3,2,2) = derivqd3(1)*rpd2*zcrd + zcrd*derivqd2(1)
        d3loc(1,3,2) = derivqd3(1)*rpd5*ycrd
        d3loc(2,3,2) = derivqd3(1)*rpd2*zcrd + zcrd*derivqd2(1)
        d3loc(3,3,2) = derivqd3(1)*rpd3*ycrd + ycrd*derivqd2(1)
!
        d3loc(1,1,3) = derivqd3(1)*rpd1*zcrd + zcrd*derivqd2(1)
        d3loc(2,1,3) = derivqd3(1)*rpd6*zcrd
        d3loc(3,1,3) = derivqd3(1)*rpd3*xcrd + xcrd*derivqd2(1)
        d3loc(1,2,3) = derivqd3(1)*rpd6*zcrd
        d3loc(2,2,3) = derivqd3(1)*rpd2*zcrd + zcrd*derivqd2(1)
        d3loc(3,2,3) = derivqd3(1)*rpd3*ycrd + ycrd*derivqd2(1)
        d3loc(1,3,3) = derivqd3(1)*rpd3*xcrd + xcrd*derivqd2(1)
        d3loc(2,3,3) = derivqd3(1)*rpd3*ycrd + ycrd*derivqd2(1)
        d3loc(3,3,3) = derivqd3(1)*rpd3*zcrd + 3.0_dp*zcrd*derivqd2(1)
!
        if (abs(poli).gt.1.0d-12) then
          derv2(jx,ix) = derv2(jx,ix) - poli*(vx(i)*d3loc(1,1,1) + vy(i)*d3loc(1,1,2) + vz(i)*d3loc(1,1,3) - &
                                              qlj*(d2loc(1,1)*d2loc(1,1) + d2loc(1,2)*d2loc(1,2) + d2loc(1,3)*d2loc(1,3)))
          derv2(jy,ix) = derv2(jy,ix) - poli*(vx(i)*d3loc(2,1,1) + vy(i)*d3loc(2,1,2) + vz(i)*d3loc(2,1,3) - &
                                              qlj*(d2loc(2,1)*d2loc(1,1) + d2loc(2,2)*d2loc(1,2) + d2loc(2,3)*d2loc(1,3)))
          derv2(jz,ix) = derv2(jz,ix) - poli*(vx(i)*d3loc(3,1,1) + vy(i)*d3loc(3,1,2) + vz(i)*d3loc(3,1,3) - &
                                              qlj*(d2loc(3,1)*d2loc(1,1) + d2loc(3,2)*d2loc(1,2) + d2loc(3,3)*d2loc(1,3)))
          derv2(jx,iy) = derv2(jx,iy) - poli*(vx(i)*d3loc(1,2,1) + vy(i)*d3loc(1,2,2) + vz(i)*d3loc(1,2,3) - &
                                              qlj*(d2loc(1,1)*d2loc(2,1) + d2loc(1,2)*d2loc(2,2) + d2loc(1,3)*d2loc(2,3)))
          derv2(jy,iy) = derv2(jy,iy) - poli*(vx(i)*d3loc(2,2,1) + vy(i)*d3loc(2,2,2) + vz(i)*d3loc(2,2,3) - &
                                              qlj*(d2loc(2,1)*d2loc(2,1) + d2loc(2,2)*d2loc(2,2) + d2loc(2,3)*d2loc(2,3)))
          derv2(jz,iy) = derv2(jz,iy) - poli*(vx(i)*d3loc(3,2,1) + vy(i)*d3loc(3,2,2) + vz(i)*d3loc(3,2,3) - &
                                              qlj*(d2loc(3,1)*d2loc(2,1) + d2loc(3,2)*d2loc(2,2) + d2loc(3,3)*d2loc(2,3)))
          derv2(jx,iz) = derv2(jx,iz) - poli*(vx(i)*d3loc(1,3,1) + vy(i)*d3loc(1,3,2) + vz(i)*d3loc(1,3,3) - &
                                              qlj*(d2loc(1,1)*d2loc(3,1) + d2loc(1,2)*d2loc(3,2) + d2loc(1,3)*d2loc(3,3)))
          derv2(jy,iz) = derv2(jy,iz) - poli*(vx(i)*d3loc(2,3,1) + vy(i)*d3loc(2,3,2) + vz(i)*d3loc(2,3,3) - &
                                              qlj*(d2loc(2,1)*d2loc(3,1) + d2loc(2,2)*d2loc(3,2) + d2loc(2,3)*d2loc(3,3)))
          derv2(jz,iz) = derv2(jz,iz) - poli*(vx(i)*d3loc(3,3,1) + vy(i)*d3loc(3,3,2) + vz(i)*d3loc(3,3,3) - &
                                              qlj*(d2loc(3,1)*d2loc(3,1) + d2loc(3,2)*d2loc(3,2) + d2loc(3,3)*d2loc(3,3)))
!
          if (lpollangevin) then
            derv2(jx,ix) = derv2(jx,ix) + pol2i*(vx(i)*d2loc(1,1) + vy(i)*d2loc(2,1) + vz(i)*d2loc(3,1))* &
                                                (vx(i)*d2loc(1,1) + vy(i)*d2loc(2,1) + vz(i)*d2loc(3,1))
            derv2(jy,ix) = derv2(jy,ix) + pol2i*(vx(i)*d2loc(1,1) + vy(i)*d2loc(2,1) + vz(i)*d2loc(3,1))* &
                                                (vx(i)*d2loc(1,2) + vy(i)*d2loc(2,2) + vz(i)*d2loc(3,2))
            derv2(jz,ix) = derv2(jz,ix) + pol2i*(vx(i)*d2loc(1,1) + vy(i)*d2loc(2,1) + vz(i)*d2loc(3,1))* &
                                                (vx(i)*d2loc(1,3) + vy(i)*d2loc(2,3) + vz(i)*d2loc(3,3))
            derv2(jx,iy) = derv2(jx,iy) + pol2i*(vx(i)*d2loc(1,2) + vy(i)*d2loc(2,2) + vz(i)*d2loc(3,2))* &
                                                (vx(i)*d2loc(1,1) + vy(i)*d2loc(2,1) + vz(i)*d2loc(3,1))
            derv2(jy,iy) = derv2(jy,iy) + pol2i*(vx(i)*d2loc(1,2) + vy(i)*d2loc(2,2) + vz(i)*d2loc(3,2))* &
                                                (vx(i)*d2loc(1,2) + vy(i)*d2loc(2,2) + vz(i)*d2loc(3,2))
            derv2(jz,iy) = derv2(jz,iy) + pol2i*(vx(i)*d2loc(1,2) + vy(i)*d2loc(2,2) + vz(i)*d2loc(3,2))* &
                                                (vx(i)*d2loc(1,3) + vy(i)*d2loc(2,3) + vz(i)*d2loc(3,3))
            derv2(jx,iz) = derv2(jx,iz) + pol2i*(vx(i)*d2loc(1,3) + vy(i)*d2loc(2,3) + vz(i)*d2loc(3,3))* &
                                                (vx(i)*d2loc(1,1) + vy(i)*d2loc(2,1) + vz(i)*d2loc(3,1))
            derv2(jy,iz) = derv2(jy,iz) + pol2i*(vx(i)*d2loc(1,3) + vy(i)*d2loc(2,3) + vz(i)*d2loc(3,3))* &
                                                (vx(i)*d2loc(1,2) + vy(i)*d2loc(2,2) + vz(i)*d2loc(3,2))
            derv2(jz,iz) = derv2(jz,iz) + pol2i*(vx(i)*d2loc(1,3) + vy(i)*d2loc(2,3) + vz(i)*d2loc(3,3))* &
                                                (vx(i)*d2loc(1,3) + vy(i)*d2loc(2,3) + vz(i)*d2loc(3,3))
          endif
        endif
!
        if (abs(polj).gt.1.0d-12) then
          derv2(jx,ix) = derv2(jx,ix) + polj*(vx(j)*d3loc(1,1,1) + vy(j)*d3loc(1,1,2) + vz(j)*d3loc(1,1,3) + &
                                              qli*(d2loc(1,1)*d2loc(1,1) + d2loc(1,2)*d2loc(1,2) + d2loc(1,3)*d2loc(1,3)))
          derv2(jy,ix) = derv2(jy,ix) + polj*(vx(j)*d3loc(2,1,1) + vy(j)*d3loc(2,1,2) + vz(j)*d3loc(2,1,3) + &
                                              qli*(d2loc(2,1)*d2loc(1,1) + d2loc(2,2)*d2loc(1,2) + d2loc(2,3)*d2loc(1,3)))
          derv2(jz,ix) = derv2(jz,ix) + polj*(vx(j)*d3loc(3,1,1) + vy(j)*d3loc(3,1,2) + vz(j)*d3loc(3,1,3) + &
                                              qli*(d2loc(3,1)*d2loc(1,1) + d2loc(3,2)*d2loc(1,2) + d2loc(3,3)*d2loc(1,3)))
          derv2(jx,iy) = derv2(jx,iy) + polj*(vx(j)*d3loc(1,2,1) + vy(j)*d3loc(1,2,2) + vz(j)*d3loc(1,2,3) + &
                                              qli*(d2loc(1,1)*d2loc(2,1) + d2loc(1,2)*d2loc(2,2) + d2loc(1,3)*d2loc(2,3)))
          derv2(jy,iy) = derv2(jy,iy) + polj*(vx(j)*d3loc(2,2,1) + vy(j)*d3loc(2,2,2) + vz(j)*d3loc(2,2,3) + &
                                              qli*(d2loc(2,1)*d2loc(2,1) + d2loc(2,2)*d2loc(2,2) + d2loc(2,3)*d2loc(2,3)))
          derv2(jz,iy) = derv2(jz,iy) + polj*(vx(j)*d3loc(3,2,1) + vy(j)*d3loc(3,2,2) + vz(j)*d3loc(3,2,3) + &
                                              qli*(d2loc(3,1)*d2loc(2,1) + d2loc(3,2)*d2loc(2,2) + d2loc(3,3)*d2loc(2,3)))
          derv2(jx,iz) = derv2(jx,iz) + polj*(vx(j)*d3loc(1,3,1) + vy(j)*d3loc(1,3,2) + vz(j)*d3loc(1,3,3) + &
                                              qli*(d2loc(1,1)*d2loc(3,1) + d2loc(1,2)*d2loc(3,2) + d2loc(1,3)*d2loc(3,3)))
          derv2(jy,iz) = derv2(jy,iz) + polj*(vx(j)*d3loc(2,3,1) + vy(j)*d3loc(2,3,2) + vz(j)*d3loc(2,3,3) + &
                                              qli*(d2loc(2,1)*d2loc(3,1) + d2loc(2,2)*d2loc(3,2) + d2loc(2,3)*d2loc(3,3)))
          derv2(jz,iz) = derv2(jz,iz) + polj*(vx(j)*d3loc(3,3,1) + vy(j)*d3loc(3,3,2) + vz(j)*d3loc(3,3,3) + &
                                              qli*(d2loc(3,1)*d2loc(3,1) + d2loc(3,2)*d2loc(3,2) + d2loc(3,3)*d2loc(3,3)))
!
          if (lpollangevin) then
            derv2(jx,ix) = derv2(jx,ix) + pol2j*(vx(j)*d2loc(1,1) + vy(j)*d2loc(2,1) + vz(j)*d2loc(3,1))* & 
                                                (vx(j)*d2loc(1,1) + vy(j)*d2loc(2,1) + vz(j)*d2loc(3,1))
            derv2(jy,ix) = derv2(jy,ix) + pol2j*(vx(j)*d2loc(1,1) + vy(j)*d2loc(2,1) + vz(j)*d2loc(3,1))* & 
                                                (vx(j)*d2loc(1,2) + vy(j)*d2loc(2,2) + vz(j)*d2loc(3,2))
            derv2(jz,ix) = derv2(jz,ix) + pol2j*(vx(j)*d2loc(1,1) + vy(j)*d2loc(2,1) + vz(j)*d2loc(3,1))* & 
                                                (vx(j)*d2loc(1,3) + vy(j)*d2loc(2,3) + vz(j)*d2loc(3,3))
            derv2(jx,iy) = derv2(jx,iy) + pol2j*(vx(j)*d2loc(1,2) + vy(j)*d2loc(2,2) + vz(j)*d2loc(3,2))* & 
                                                (vx(j)*d2loc(1,1) + vy(j)*d2loc(2,1) + vz(j)*d2loc(3,1))
            derv2(jy,iy) = derv2(jy,iy) + pol2j*(vx(j)*d2loc(1,2) + vy(j)*d2loc(2,2) + vz(j)*d2loc(3,2))* & 
                                                (vx(j)*d2loc(1,2) + vy(j)*d2loc(2,2) + vz(j)*d2loc(3,2))
            derv2(jz,iy) = derv2(jz,iy) + pol2j*(vx(j)*d2loc(1,2) + vy(j)*d2loc(2,2) + vz(j)*d2loc(3,2))* & 
                                                (vx(j)*d2loc(1,3) + vy(j)*d2loc(2,3) + vz(j)*d2loc(3,3))
            derv2(jx,iz) = derv2(jx,iz) + pol2j*(vx(j)*d2loc(1,3) + vy(j)*d2loc(2,3) + vz(j)*d2loc(3,3))* & 
                                                (vx(j)*d2loc(1,1) + vy(j)*d2loc(2,1) + vz(j)*d2loc(3,1))
            derv2(jy,iz) = derv2(jy,iz) + pol2j*(vx(j)*d2loc(1,3) + vy(j)*d2loc(2,3) + vz(j)*d2loc(3,3))* & 
                                                (vx(j)*d2loc(1,2) + vy(j)*d2loc(2,2) + vz(j)*d2loc(3,2))
            derv2(jz,iz) = derv2(jz,iz) + pol2j*(vx(j)*d2loc(1,3) + vy(j)*d2loc(2,3) + vz(j)*d2loc(3,3))* & 
                                                (vx(j)*d2loc(1,3) + vy(j)*d2loc(2,3) + vz(j)*d2loc(3,3))
          endif
        endif
      endif
    enddo jloop
!
    if (lgrad2) then
!***********************************************************************************
!  Threebody contribution to second derivatives from j/k with polarisability of i  *
!***********************************************************************************
!
!  Loop over j
!
      jx = -2
      jy = -1
      jz =  0
      jloop2: do j = 1,numat
        jx = jx + 3
        jy = jy + 3
        jz = jz + 3
!
        if (j.eq.i) cycle jloop2
!
        qlj = qf(j)
!
!  Loop over k
!
        kx = -2
        ky = -1
        kz =  0
        kloop: do k = 1,j-1
          kx = kx + 3
          ky = ky + 3
          kz = kz + 3
!
          if (k.eq.i) cycle kloop
!
          qlk = qf(k)
          polk = qlj*qlk*dpolar(i)/angstoev
          if (lpollangevin) then
            Es = sqrt(vx(i)**2 + vy(i)**2 + vz(i)**2)
            if (Es.gt.1.0d-8) then
              call langevinpol(Es,dpolarmax(i),dpolar(i),Ep,dEpdEs,d2EpdEs2,.true.,lgrad2)
              polk = qlj*qlk*dEpdEs/angstoev
              if (lgrad2) then
                pol2k = qlj*qlk*d2EpdEs2/angstoev
              endif
            else
              polk = qlj*qlk*dpolar(i)/angstoev
              pol2k = 0.0_dp
            endif
          else
            polk = qlj*qlk*dpolar(i)/angstoev
            pol2k = 0.0_dp
          endif
!
          if (abs(polk).gt.1.0d-12) then
            derv2(jx,ix) = derv2(jx,ix) + polk*(d2k(1,1,j)*d2k(1,1,k) + d2k(1,2,j)*d2k(1,2,k) + d2k(1,3,j)*d2k(1,3,k))
            derv2(jy,ix) = derv2(jy,ix) + polk*(d2k(2,1,j)*d2k(1,1,k) + d2k(2,2,j)*d2k(1,2,k) + d2k(2,3,j)*d2k(1,3,k))
            derv2(jz,ix) = derv2(jz,ix) + polk*(d2k(3,1,j)*d2k(1,1,k) + d2k(3,2,j)*d2k(1,2,k) + d2k(3,3,j)*d2k(1,3,k))
            derv2(jx,iy) = derv2(jx,iy) + polk*(d2k(1,1,j)*d2k(2,1,k) + d2k(1,2,j)*d2k(2,2,k) + d2k(1,3,j)*d2k(2,3,k))
            derv2(jy,iy) = derv2(jy,iy) + polk*(d2k(2,1,j)*d2k(2,1,k) + d2k(2,2,j)*d2k(2,2,k) + d2k(2,3,j)*d2k(2,3,k))
            derv2(jz,iy) = derv2(jz,iy) + polk*(d2k(3,1,j)*d2k(2,1,k) + d2k(3,2,j)*d2k(2,2,k) + d2k(3,3,j)*d2k(2,3,k))
            derv2(jx,iz) = derv2(jx,iz) + polk*(d2k(1,1,j)*d2k(3,1,k) + d2k(1,2,j)*d2k(3,2,k) + d2k(1,3,j)*d2k(3,3,k))
            derv2(jy,iz) = derv2(jy,iz) + polk*(d2k(2,1,j)*d2k(3,1,k) + d2k(2,2,j)*d2k(3,2,k) + d2k(2,3,j)*d2k(3,3,k))
            derv2(jz,iz) = derv2(jz,iz) + polk*(d2k(3,1,j)*d2k(3,1,k) + d2k(3,2,j)*d2k(3,2,k) + d2k(3,3,j)*d2k(3,3,k))
!
            derv2(kx,ix) = derv2(kx,ix) + polk*(d2k(1,1,k)*d2k(1,1,j) + d2k(1,2,k)*d2k(1,2,j) + d2k(1,3,k)*d2k(1,3,j))
            derv2(ky,ix) = derv2(ky,ix) + polk*(d2k(2,1,k)*d2k(1,1,j) + d2k(2,2,k)*d2k(1,2,j) + d2k(2,3,k)*d2k(1,3,j))
            derv2(kz,ix) = derv2(kz,ix) + polk*(d2k(3,1,k)*d2k(1,1,j) + d2k(3,2,k)*d2k(1,2,j) + d2k(3,3,k)*d2k(1,3,j))
            derv2(kx,iy) = derv2(kx,iy) + polk*(d2k(1,1,k)*d2k(2,1,j) + d2k(1,2,k)*d2k(2,2,j) + d2k(1,3,k)*d2k(2,3,j))
            derv2(ky,iy) = derv2(ky,iy) + polk*(d2k(2,1,k)*d2k(2,1,j) + d2k(2,2,k)*d2k(2,2,j) + d2k(2,3,k)*d2k(2,3,j))
            derv2(kz,iy) = derv2(kz,iy) + polk*(d2k(3,1,k)*d2k(2,1,j) + d2k(3,2,k)*d2k(2,2,j) + d2k(3,3,k)*d2k(2,3,j))
            derv2(kx,iz) = derv2(kx,iz) + polk*(d2k(1,1,k)*d2k(3,1,j) + d2k(1,2,k)*d2k(3,2,j) + d2k(1,3,k)*d2k(3,3,j))
            derv2(ky,iz) = derv2(ky,iz) + polk*(d2k(2,1,k)*d2k(3,1,j) + d2k(2,2,k)*d2k(3,2,j) + d2k(2,3,k)*d2k(3,3,j))
            derv2(kz,iz) = derv2(kz,iz) + polk*(d2k(3,1,k)*d2k(3,1,j) + d2k(3,2,k)*d2k(3,2,j) + d2k(3,3,k)*d2k(3,3,j))
!
            if (lpollangevin) then
              derv2(jx,ix) = derv2(jx,ix) + pol2k*(vx(i)*d2k(1,1,j) + vy(i)*d2k(2,1,j) + vz(i)*d2k(3,1,j))* &
                                                  (vx(i)*d2k(1,1,k) + vy(i)*d2k(2,1,k) + vz(i)*d2k(3,1,k))
              derv2(jx,iy) = derv2(jx,iy) + pol2k*(vx(i)*d2k(1,1,j) + vy(i)*d2k(2,1,j) + vz(i)*d2k(3,1,j))* &
                                                  (vx(i)*d2k(1,2,k) + vy(i)*d2k(2,2,k) + vz(i)*d2k(3,2,k))
              derv2(jx,iz) = derv2(jx,iz) + pol2k*(vx(i)*d2k(1,1,j) + vy(i)*d2k(2,1,j) + vz(i)*d2k(3,1,j))* &
                                                  (vx(i)*d2k(1,3,k) + vy(i)*d2k(2,3,k) + vz(i)*d2k(3,3,k))
              derv2(jy,ix) = derv2(jy,ix) + pol2k*(vx(i)*d2k(1,2,j) + vy(i)*d2k(2,2,j) + vz(i)*d2k(3,2,j))* &
                                                  (vx(i)*d2k(1,1,k) + vy(i)*d2k(2,1,k) + vz(i)*d2k(3,1,k))
              derv2(jy,iy) = derv2(jy,iy) + pol2k*(vx(i)*d2k(1,2,j) + vy(i)*d2k(2,2,j) + vz(i)*d2k(3,2,j))* &
                                                  (vx(i)*d2k(1,2,k) + vy(i)*d2k(2,2,k) + vz(i)*d2k(3,2,k))
              derv2(jy,iz) = derv2(jy,iz) + pol2k*(vx(i)*d2k(1,2,j) + vy(i)*d2k(2,2,j) + vz(i)*d2k(3,2,j))* &
                                                  (vx(i)*d2k(1,3,k) + vy(i)*d2k(2,3,k) + vz(i)*d2k(3,3,k))
              derv2(jz,ix) = derv2(jz,ix) + pol2k*(vx(i)*d2k(1,3,j) + vy(i)*d2k(2,3,j) + vz(i)*d2k(3,3,j))* &
                                                  (vx(i)*d2k(1,1,k) + vy(i)*d2k(2,1,k) + vz(i)*d2k(3,1,k))
              derv2(jz,iy) = derv2(jz,iy) + pol2k*(vx(i)*d2k(1,3,j) + vy(i)*d2k(2,3,j) + vz(i)*d2k(3,3,j))* &
                                                  (vx(i)*d2k(1,2,k) + vy(i)*d2k(2,2,k) + vz(i)*d2k(3,2,k))
              derv2(jz,iz) = derv2(jz,iz) + pol2k*(vx(i)*d2k(1,3,j) + vy(i)*d2k(2,3,j) + vz(i)*d2k(3,3,j))* &
                                                  (vx(i)*d2k(1,3,k) + vy(i)*d2k(2,3,k) + vz(i)*d2k(3,3,k))
!
              derv2(kx,ix) = derv2(kx,ix) + pol2k*(vx(i)*d2k(1,1,k) + vy(i)*d2k(2,1,k) + vz(i)*d2k(3,1,k))* &
                                                  (vx(i)*d2k(1,1,j) + vy(i)*d2k(2,1,j) + vz(i)*d2k(3,1,j))
              derv2(kx,iy) = derv2(kx,iy) + pol2k*(vx(i)*d2k(1,1,k) + vy(i)*d2k(2,1,k) + vz(i)*d2k(3,1,k))* &
                                                  (vx(i)*d2k(1,2,j) + vy(i)*d2k(2,2,j) + vz(i)*d2k(3,2,j))
              derv2(kx,iz) = derv2(kx,iz) + pol2k*(vx(i)*d2k(1,1,k) + vy(i)*d2k(2,1,k) + vz(i)*d2k(3,1,k))* &
                                                  (vx(i)*d2k(1,3,j) + vy(i)*d2k(2,3,j) + vz(i)*d2k(3,3,j))
              derv2(ky,ix) = derv2(ky,ix) + pol2k*(vx(i)*d2k(1,2,k) + vy(i)*d2k(2,2,k) + vz(i)*d2k(3,2,k))* &
                                                  (vx(i)*d2k(1,1,j) + vy(i)*d2k(2,1,j) + vz(i)*d2k(3,1,j))
              derv2(ky,iy) = derv2(ky,iy) + pol2k*(vx(i)*d2k(1,2,k) + vy(i)*d2k(2,2,k) + vz(i)*d2k(3,2,k))* &
                                                  (vx(i)*d2k(1,2,j) + vy(i)*d2k(2,2,j) + vz(i)*d2k(3,2,j))
              derv2(ky,iz) = derv2(ky,iz) + pol2k*(vx(i)*d2k(1,2,k) + vy(i)*d2k(2,2,k) + vz(i)*d2k(3,2,k))* &
                                                  (vx(i)*d2k(1,3,j) + vy(i)*d2k(2,3,j) + vz(i)*d2k(3,3,j))
              derv2(kz,ix) = derv2(kz,ix) + pol2k*(vx(i)*d2k(1,3,k) + vy(i)*d2k(2,3,k) + vz(i)*d2k(3,3,k))* &
                                                  (vx(i)*d2k(1,1,j) + vy(i)*d2k(2,1,j) + vz(i)*d2k(3,1,j))
              derv2(kz,iy) = derv2(kz,iy) + pol2k*(vx(i)*d2k(1,3,k) + vy(i)*d2k(2,3,k) + vz(i)*d2k(3,3,k))* &
                                                  (vx(i)*d2k(1,2,j) + vy(i)*d2k(2,2,j) + vz(i)*d2k(3,2,j))
              derv2(kz,iz) = derv2(kz,iz) + pol2k*(vx(i)*d2k(1,3,k) + vy(i)*d2k(2,3,k) + vz(i)*d2k(3,3,k))* &
                                                  (vx(i)*d2k(1,3,j) + vy(i)*d2k(2,3,j) + vz(i)*d2k(3,3,j))
            endif
          endif
        enddo kloop
      enddo jloop2
    endif
  enddo iloop
!
  if (lgrad2) then
!************************************************************************************
!  Outer loop over sites with polarisability on site i for second derivative terms  *
!************************************************************************************
    ix = -2
    iy = -1
    iz =  0
!
    iloop3: do i = 1,numat
      ix = ix + 3
      iy = iy + 3
      iz = iz + 3
!
!  If i has no polarisability then skip this atom
!
      if (dpolar(i).eq.0.0_dp) cycle iloop3
!
!  Initialise d2k array
!
      d2k(1:3,1:3,1:numat) = 0.0_dp
!
      xal = xclat(i)
      yal = yclat(i)
      zal = zclat(i)
      nati = nat(i)
      ntypi = nftype(i)
      nregioni = nregionno(nsft+nrelf2a(i))
      qli = qf(i)
      oci = occuf(i)
!
!  Molecule handling
!
      if (lmol) then
        nmi = natmol(i)
        indmi = nmolind(i)
!
!  Set COM coordinates
!
        if (lrigid.and.nmi.gt.0) then
          xcomi = molxyz(1,natinmol(i),nmi)
          ycomi = molxyz(2,natinmol(i),nmi)
          zcomi = molxyz(3,natinmol(i),nmi)
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
!  Inner loop over second site to generate j terms for atoms on this node
!
      jloop3: do jonnode = 1,natomsonnode
        j = node2atom(jonnode)
!
!  Exclude i = j in 0-D
!
        if (j.eq.i) cycle jloop3
!
        natj = nat(j)
        ntypj = nftype(j)
        nregionj = nregionno(nsft+nrelf2a(j))
        qlj = qf(j)
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
            ntyp1 = ntypi
            ntyp2 = ntypj
          else
            ntyp1 = ntypj
            ntyp2 = ntypi
          endif
        elseif (nati.lt.natj) then
          nat1 = nati
          nat2 = natj
          ntyp1 = ntypi
          ntyp2 = ntypj
        else
          nat1 = natj
          nat2 = nati
          ntyp1 = ntypj
          ntyp2 = ntypi
        endif
        ofct = oci*ocj
        fct = ofct*angstoev
        factor = qli*qlj*fct
!
!  Possible core - shell flag
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
        r = xcrd*xcrd + ycrd*ycrd + zcrd*zcrd
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
!
        if (lmolok.and.(r.gt.cut2s.or..not.lcspair)) then
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
            lbonded = .false.
            l2bonds = .false.
            l3bonds = .false.
          endif
        else
          lptrmol = .false.
          lbonded = .false.
          l2bonds = .false.
          l3bonds = .false.
        endif
        if (abs(r-small2).lt.1.0d-12) r = small2
        if (r.lt.small) then
          cycle jloop3
        else
!
!  Store vector
!
          nor = 1
          dist = sqrt(r)
        endif
!*************************************************
!  Calculate twobody contribution in real space  *
!*************************************************
        eatom = 0.0_dp
        ereal = 0.0_dp
        ec6 = 0.0_dp
        call twobody1(eatom,ereal,ec6,.true.,.true.,lgrad2,nor,1_i4,dist,xcrd,ycrd,zcrd,d0i,d0j, &
                      deriv,deriv2,deriv3,derive0,derive,derive2,derive3,rderiv,npots,npotl,0.0_dp, &
                      cut2q,cut2s,lptrmol,nmolonly,factor,ofct,ospfct,0.0_dp,rtrm1,rtrm2,rtrm3,rtrm32, &
                      sctrm1,sctrm2,qli,qlj,lcspair,.false.,.false.,lbonded,l2bonds,l3bonds, &
                      nbtypeij,nbtypeij2,.false.,.false.,d1iq,d1jq,d2i2q,d2ijq,d2j2q)
!**************************
!  Coordinate Derivatives *
!**************************
!
!  First derivatives
!
        rpd1 = xcrd*xcrd
        rpd2 = ycrd*ycrd
        rpd3 = zcrd*zcrd
        rpd4 = ycrd*zcrd
        rpd5 = xcrd*zcrd
        rpd6 = xcrd*ycrd
!
        d2loc(1,1) = derivqd2(1)*rpd1 + derivqd(1)
        d2loc(2,1) = derivqd2(1)*rpd6
        d2loc(3,1) = derivqd2(1)*rpd5
        d2loc(1,2) = derivqd2(1)*rpd6
        d2loc(2,2) = derivqd2(1)*rpd2 + derivqd(1)
        d2loc(3,2) = derivqd2(1)*rpd4
        d2loc(1,3) = derivqd2(1)*rpd5
        d2loc(2,3) = derivqd2(1)*rpd4
        d2loc(3,3) = derivqd2(1)*rpd3 + derivqd(1)
!
!  Save d2 term
!
        d2k(1:3,1:3,j) = d2loc(1:3,1:3)
!
!  Skip to here if frozen pair
!
      enddo jloop3
!***************************************
!  Now globalise d2k across all nodes  *
!***************************************
      call sumall(d2k,d2ksum,9_i4*numat,"preal0dd","d2k")
      d2k(1:3,1:3,1:numat) = d2ksum(1:3,1:3,1:numat)
!***********************************************************************************
!  Threebody contribution to second derivatives from j/k with polarisability of i  *
!***********************************************************************************
!
!  Loop over j
!
      jx = -2
      jy = -1
      jz =  0
      jloop4: do jonnode = 1,natomsonnode
        j = node2atom(jonnode)
!
        jx = jx + 3
        jy = jy + 3
        jz = jz + 3
!
        if (j.eq.i) cycle jloop4
!
        qlj = qf(j)
!
!  Loop over k
!
        kx = -2
        ky = -1
        kz =  0
        kloop3: do k = 1,numat
          kx = kx + 3
          ky = ky + 3
          kz = kz + 3
!
          if (k.eq.i.or.k.eq.j) cycle kloop3
!
          qlk = qf(k)
          if (lpollangevin) then
            Es = sqrt(vx(i)**2 + vy(i)**2 + vz(i)**2)
            if (Es.gt.1.0d-8) then
              call langevinpol(Es,dpolarmax(i),dpolar(i),Ep,dEpdEs,d2EpdEs2,.true.,lgrad2)
              polk = qlj*qlk*dEpdEs/angstoev
              if (lgrad2) then
                pol2k = qlj*qlk*d2EpdEs2/angstoev
              endif
            else
              polk = qlj*qlk*dpolar(i)/angstoev
              pol2k = 0.0_dp
            endif
          else
            polk = qlj*qlk*dpolar(i)/angstoev
            pol2k = 0.0_dp
          endif
!
          if (abs(polk).gt.1.0d-12) then
            derv2(ix,jx) = derv2(ix,jx) + polk*(d2k(1,1,k)*d2k(1,1,j) + d2k(1,2,k)*d2k(1,2,j) + d2k(1,3,k)*d2k(1,3,j))
            derv2(iy,jx) = derv2(iy,jx) + polk*(d2k(2,1,k)*d2k(1,1,j) + d2k(2,2,k)*d2k(1,2,j) + d2k(2,3,k)*d2k(1,3,j))
            derv2(iz,jx) = derv2(iz,jx) + polk*(d2k(3,1,k)*d2k(1,1,j) + d2k(3,2,k)*d2k(1,2,j) + d2k(3,3,k)*d2k(1,3,j))
            derv2(ix,jy) = derv2(ix,jy) + polk*(d2k(1,1,k)*d2k(2,1,j) + d2k(1,2,k)*d2k(2,2,j) + d2k(1,3,k)*d2k(2,3,j))
            derv2(iy,jy) = derv2(iy,jy) + polk*(d2k(2,1,k)*d2k(2,1,j) + d2k(2,2,k)*d2k(2,2,j) + d2k(2,3,k)*d2k(2,3,j))
            derv2(iz,jy) = derv2(iz,jy) + polk*(d2k(3,1,k)*d2k(2,1,j) + d2k(3,2,k)*d2k(2,2,j) + d2k(3,3,k)*d2k(2,3,j))
            derv2(ix,jz) = derv2(ix,jz) + polk*(d2k(1,1,k)*d2k(3,1,j) + d2k(1,2,k)*d2k(3,2,j) + d2k(1,3,k)*d2k(3,3,j))
            derv2(iy,jz) = derv2(iy,jz) + polk*(d2k(2,1,k)*d2k(3,1,j) + d2k(2,2,k)*d2k(3,2,j) + d2k(2,3,k)*d2k(3,3,j))
            derv2(iz,jz) = derv2(iz,jz) + polk*(d2k(3,1,k)*d2k(3,1,j) + d2k(3,2,k)*d2k(3,2,j) + d2k(3,3,k)*d2k(3,3,j))
!
            if (lpollangevin) then
              derv2(ix,jx) = derv2(ix,jx) + pol2k*(vx(i)*d2k(1,1,k) + vy(i)*d2k(2,1,k) + vz(i)*d2k(3,1,k))* &
                                                  (vx(i)*d2k(1,1,j) + vy(i)*d2k(2,1,j) + vz(i)*d2k(3,1,j))
              derv2(ix,jy) = derv2(ix,jy) + pol2k*(vx(i)*d2k(1,1,k) + vy(i)*d2k(2,1,k) + vz(i)*d2k(3,1,k))* &
                                                  (vx(i)*d2k(1,2,j) + vy(i)*d2k(2,2,j) + vz(i)*d2k(3,2,j))
              derv2(ix,jz) = derv2(ix,jz) + pol2k*(vx(i)*d2k(1,1,k) + vy(i)*d2k(2,1,k) + vz(i)*d2k(3,1,k))* &
                                                  (vx(i)*d2k(1,3,j) + vy(i)*d2k(2,3,j) + vz(i)*d2k(3,3,j))
              derv2(iy,jx) = derv2(iy,jx) + pol2k*(vx(i)*d2k(1,2,k) + vy(i)*d2k(2,2,k) + vz(i)*d2k(3,2,k))* &
                                                  (vx(i)*d2k(1,1,j) + vy(i)*d2k(2,1,j) + vz(i)*d2k(3,1,j))
              derv2(iy,jy) = derv2(iy,jy) + pol2k*(vx(i)*d2k(1,2,k) + vy(i)*d2k(2,2,k) + vz(i)*d2k(3,2,k))* &
                                                  (vx(i)*d2k(1,2,j) + vy(i)*d2k(2,2,j) + vz(i)*d2k(3,2,j))
              derv2(iy,jz) = derv2(iy,jz) + pol2k*(vx(i)*d2k(1,2,k) + vy(i)*d2k(2,2,k) + vz(i)*d2k(3,2,k))* &
                                                  (vx(i)*d2k(1,3,j) + vy(i)*d2k(2,3,j) + vz(i)*d2k(3,3,j))
              derv2(iz,jx) = derv2(iz,jx) + pol2k*(vx(i)*d2k(1,3,k) + vy(i)*d2k(2,3,k) + vz(i)*d2k(3,3,k))* &
                                                  (vx(i)*d2k(1,1,j) + vy(i)*d2k(2,1,j) + vz(i)*d2k(3,1,j))
              derv2(iz,jy) = derv2(iz,jy) + pol2k*(vx(i)*d2k(1,3,k) + vy(i)*d2k(2,3,k) + vz(i)*d2k(3,3,k))* &
                                                  (vx(i)*d2k(1,2,j) + vy(i)*d2k(2,2,j) + vz(i)*d2k(3,2,j))
              derv2(iz,jz) = derv2(iz,jz) + pol2k*(vx(i)*d2k(1,3,k) + vy(i)*d2k(2,3,k) + vz(i)*d2k(3,3,k))* &
                                                  (vx(i)*d2k(1,3,j) + vy(i)*d2k(2,3,j) + vz(i)*d2k(3,3,j))
            endif
          endif
        enddo kloop3
!
        kx = -2
        ky = -1
        kz =  0
        kloop4: do k = 1,numat
          kx = kx + 3
          ky = ky + 3
          kz = kz + 3
!
          if (k.eq.i.or.k.eq.j) cycle kloop4
!
          qlk = qf(k)
          if (lpollangevin) then
            Es = sqrt(vx(i)**2 + vy(i)**2 + vz(i)**2)
            if (Es.gt.1.0d-8) then
              call langevinpol(Es,dpolarmax(i),dpolar(i),Ep,dEpdEs,d2EpdEs2,.true.,lgrad2)
              polk = qlj*qlk*dEpdEs/angstoev
              if (lgrad2) then
                pol2k = qlj*qlk*d2EpdEs2/angstoev
              endif
            else
              polk = qlj*qlk*dpolar(i)/angstoev
              pol2k = 0.0_dp
            endif
          else
            polk = qlj*qlk*dpolar(i)/angstoev
            pol2k = 0.0_dp
          endif
!
          if (abs(polk).gt.1.0d-12) then
            derv2(kx,jx) = derv2(kx,jx) - polk*(d2k(1,1,k)*d2k(1,1,j) + d2k(1,2,k)*d2k(1,2,j) + d2k(1,3,k)*d2k(1,3,j))
            derv2(ky,jx) = derv2(ky,jx) - polk*(d2k(2,1,k)*d2k(1,1,j) + d2k(2,2,k)*d2k(1,2,j) + d2k(2,3,k)*d2k(1,3,j))
            derv2(kz,jx) = derv2(kz,jx) - polk*(d2k(3,1,k)*d2k(1,1,j) + d2k(3,2,k)*d2k(1,2,j) + d2k(3,3,k)*d2k(1,3,j))
            derv2(kx,jy) = derv2(kx,jy) - polk*(d2k(1,1,k)*d2k(2,1,j) + d2k(1,2,k)*d2k(2,2,j) + d2k(1,3,k)*d2k(2,3,j))
            derv2(ky,jy) = derv2(ky,jy) - polk*(d2k(2,1,k)*d2k(2,1,j) + d2k(2,2,k)*d2k(2,2,j) + d2k(2,3,k)*d2k(2,3,j))
            derv2(kz,jy) = derv2(kz,jy) - polk*(d2k(3,1,k)*d2k(2,1,j) + d2k(3,2,k)*d2k(2,2,j) + d2k(3,3,k)*d2k(2,3,j))
            derv2(kx,jz) = derv2(kx,jz) - polk*(d2k(1,1,k)*d2k(3,1,j) + d2k(1,2,k)*d2k(3,2,j) + d2k(1,3,k)*d2k(3,3,j))
            derv2(ky,jz) = derv2(ky,jz) - polk*(d2k(2,1,k)*d2k(3,1,j) + d2k(2,2,k)*d2k(3,2,j) + d2k(2,3,k)*d2k(3,3,j))
            derv2(kz,jz) = derv2(kz,jz) - polk*(d2k(3,1,k)*d2k(3,1,j) + d2k(3,2,k)*d2k(3,2,j) + d2k(3,3,k)*d2k(3,3,j))
!
            if (lpollangevin) then
              derv2(kx,jx) = derv2(kx,jx) - pol2k*(vx(i)*d2k(1,1,k) + vy(i)*d2k(2,1,k) + vz(i)*d2k(3,1,k))* &
                                                  (vx(i)*d2k(1,1,j) + vy(i)*d2k(2,1,j) + vz(i)*d2k(3,1,j))
              derv2(kx,jy) = derv2(kx,jy) - pol2k*(vx(i)*d2k(1,1,k) + vy(i)*d2k(2,1,k) + vz(i)*d2k(3,1,k))* &
                                                  (vx(i)*d2k(1,2,j) + vy(i)*d2k(2,2,j) + vz(i)*d2k(3,2,j))
              derv2(kx,jz) = derv2(kx,jz) - pol2k*(vx(i)*d2k(1,1,k) + vy(i)*d2k(2,1,k) + vz(i)*d2k(3,1,k))* &
                                                  (vx(i)*d2k(1,3,j) + vy(i)*d2k(2,3,j) + vz(i)*d2k(3,3,j))
              derv2(ky,jx) = derv2(ky,jx) - pol2k*(vx(i)*d2k(1,2,k) + vy(i)*d2k(2,2,k) + vz(i)*d2k(3,2,k))* &
                                                  (vx(i)*d2k(1,1,j) + vy(i)*d2k(2,1,j) + vz(i)*d2k(3,1,j))
              derv2(ky,jy) = derv2(ky,jy) - pol2k*(vx(i)*d2k(1,2,k) + vy(i)*d2k(2,2,k) + vz(i)*d2k(3,2,k))* &
                                                  (vx(i)*d2k(1,2,j) + vy(i)*d2k(2,2,j) + vz(i)*d2k(3,2,j))
              derv2(ky,jz) = derv2(ky,jz) - pol2k*(vx(i)*d2k(1,2,k) + vy(i)*d2k(2,2,k) + vz(i)*d2k(3,2,k))* &
                                                  (vx(i)*d2k(1,3,j) + vy(i)*d2k(2,3,j) + vz(i)*d2k(3,3,j))
              derv2(kz,jx) = derv2(kz,jx) - pol2k*(vx(i)*d2k(1,3,k) + vy(i)*d2k(2,3,k) + vz(i)*d2k(3,3,k))* &
                                                  (vx(i)*d2k(1,1,j) + vy(i)*d2k(2,1,j) + vz(i)*d2k(3,1,j))
              derv2(kz,jy) = derv2(kz,jy) - pol2k*(vx(i)*d2k(1,3,k) + vy(i)*d2k(2,3,k) + vz(i)*d2k(3,3,k))* &
                                                  (vx(i)*d2k(1,2,j) + vy(i)*d2k(2,2,j) + vz(i)*d2k(3,2,j))
              derv2(kz,jz) = derv2(kz,jz) - pol2k*(vx(i)*d2k(1,3,k) + vy(i)*d2k(2,3,k) + vz(i)*d2k(3,3,k))* &
                                                  (vx(i)*d2k(1,3,j) + vy(i)*d2k(2,3,j) + vz(i)*d2k(3,3,j))
            endif
          endif
        enddo kloop4
      enddo jloop4
    enddo iloop3
  endif
!
!  Free local memory
!
  if (lgrad2) then
    deallocate(d2ksum,stat=status)
    if (status/=0) call deallocate_error('pireal0dd','d2ksum')
    deallocate(d2k,stat=status)
    if (status/=0) call deallocate_error('pireal0dd','d2k')
  endif
  deallocate(npotl,stat=status)
  if (status/=0) call deallocate_error('pireal0dd','npotl')
!
!  End of real space part - perform general tasks
!
999 continue
!
!  Timing
!
  time2 = g_cpu_time()
  tpolar = tpolar + time2 - time1
#ifdef TRACE
  call trace_out('pireal0dd')
#endif
!
  return
  end
