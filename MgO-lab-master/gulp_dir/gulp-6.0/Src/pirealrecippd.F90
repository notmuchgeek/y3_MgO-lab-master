  subroutine pirealrecippd
!
!  Calculates the second derivatives of the polarisation energy.
!  Distributed memory parallel version for phonons at gamma point.
!
!   7/19 Created from pirealrecipd
!   8/19 Trapping use of uninitialised d2a and d2as added
!   8/19 Short range damping of polarisation added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!  10/19 Langevin damping added
!  11/19 Langevin second derivatives added
!   2/20 Correction to Langevin damped formalism
!   2/20 Modified for rigid molecules
!   4/20 derv3c and d2r2dsdc added for benefit of rigid molecules
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
  use m_strain,       only : twostrterms
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
  integer(i4)                                  :: ii
  integer(i4)                                  :: ionnode
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
  integer(i4)                                  :: jonnode
  integer(i4)                                  :: jx
  integer(i4)                                  :: jy
  integer(i4)                                  :: jz
  integer(i4)                                  :: k
  integer(i4)                                  :: kk
  integer(i4)                                  :: kx
  integer(i4)                                  :: ky
  integer(i4)                                  :: kz
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
  integer(i4)                                  :: ntyp1
  integer(i4)                                  :: ntyp2
  integer(i4)                                  :: ntypi
  integer(i4)                                  :: ntypj
  integer(i4)                                  :: status
  logical                                      :: lcspair
  logical                                      :: lewaldtype
  logical                                      :: lmolok
  logical                                      :: lneedmol
  logical                                      :: lself
  logical,     dimension(:),     allocatable   :: ld2ok
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
  real(dp)                                     :: eatom
  real(dp)                                     :: ec6
  real(dp)                                     :: ereal
  real(dp)                                     :: Ep
  real(dp)                                     :: dEpdEs
  real(dp)                                     :: d2EpdEs2
  real(dp)                                     :: Es
  real(dp)                                     :: factor
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
  real(dp)                                     :: rp
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
  real(dp)                                     :: xal
  real(dp)                                     :: yal
  real(dp)                                     :: zal
  real(dp)                                     :: xcom
  real(dp)                                     :: ycom
  real(dp)                                     :: zcom
  real(dp)                                     :: xcrd
  real(dp)                                     :: xcomi
  real(dp)                                     :: ycomi
  real(dp)                                     :: zcomi
  real(dp)                                     :: ycrd
  real(dp)                                     :: zcrd
  real(dp),    dimension(:,:,:), allocatable   :: d2a
  real(dp),    dimension(:,:,:), allocatable   :: d2asum
#ifdef TRACE
  call trace_in('pirealrecippd')
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
!  Allocate local memory
!
  allocate(ld2ok(numat),stat=status)
  if (status/=0) call outofmemory('pirealrecippd','ld2ok')
  allocate(npotl(npote),stat=status)
  if (status/=0) call outofmemory('pirealrecippd','npotl')
  allocate(d2a(3,3,numat),stat=status)
  if (status/=0) call outofmemory('pirealrecippd','d2a')
  allocate(d2asum(3,3,numat),stat=status)
  if (status/=0) call outofmemory('pirealrecippd','d2asum')
!
!  Initialise K vector terms
!
  call setktrmdp
!***************************************************************
!  Atomistic and real space electrostatic second derivatives   *
!***************************************************************
!
!  Outer loop over sites with polarisability on site i
!
  ix = -2
  iy = -1
  iz =  0
!
!  Loop over atoms on local node
!
  iloop: do ionnode = 1,natomsonnode
    i = node2atom(ionnode)
!
    ix = ix + 3
    iy = iy + 3
    iz = iz + 3
!
!  If i has no polarisability or charge then skip this atom
!
    qli = qf(i)
    if (dpolar(nrelf2a(i)).eq.0.0_dp.and.qli.eq.0.0_dp) cycle iloop
!
!  Inner loop over second site
!
    xal = xclat(i)
    yal = yclat(i)
    zal = zclat(i)
    nati = nat(i)
    ntypi = nftype(i)
    oci = occuf(i)
!
!  Molecule handling
!
    if (lmol) then
      nmi = natmol(i)
      indm = nmolind(i)
      call mindtoijk(indm,ixi,iyi,izi)
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
!  Initialise flag as to whether d2a(s) arrays are set
!
    ld2ok(1:numat) = .false.
!
!  Start of second atom loop
!
    jx = -2
    jy = -1
    jz =  0
    jloop: do j = 1,numat
      jx = jx + 3
      jy = jy + 3
      jz = jz + 3
!
      qlj = qf(j)
!
!  Compute charge x polarisability for pair
!
      if (lpollangevin) then
        Es = sqrt(vx(i)**2 + vy(i)**2 + vz(i)**2)
        if (Es.gt.1.0d-8) then
          call langevinpol(Es,dpolarmax(nrelf2a(i)),dpolar(nrelf2a(i)),Ep,dEpdEs,d2EpdEs2,.true.,.true.)
          poli = qlj*dEpdEs/angstoev
          pol2i = qlj*qlj*d2EpdEs2/angstoev
        else
          poli = qlj*dpolar(nrelf2a(i))/angstoev
          pol2i = 0.0_dp
        endif
!
        Es = sqrt(vx(j)**2 + vy(j)**2 + vz(j)**2)
        if (Es.gt.1.0d-8) then
          call langevinpol(Es,dpolarmax(nrelf2a(j)),dpolar(nrelf2a(j)),Ep,dEpdEs,d2EpdEs2,.true.,.true.)
          polj = qli*dEpdEs/angstoev
          pol2j = qli*qli*d2EpdEs2/angstoev
        else
          polj = qli*dpolar(nrelf2a(j))/angstoev
          pol2j = 0.0_dp
        endif
      else
        poli = qlj*dpolar(nrelf2a(i))/angstoev
        polj = qli*dpolar(nrelf2a(j))/angstoev
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
!
!  Zero third derivative arrays
!
      d3(1:3,1:3,1:3) = 0.0_dp
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
!
!  If no valid distances then skip rest of this loop
!
      if (nor.eq.0) cycle jloop
!
      if (.not.lnoreal) then
!
!  Sqrt distances
!
        do k = 1,nor
          dist(k) = sqrt(dist(k))
        enddo
!
        call twobody(eatom,ereal,ec6,.true.,.true.,.true.,nor,1,npots,npotl,cut2,cut2q,cut2s, &
                     nmolonly,factor,ofct,ospfct,0.0_dp,sctrm1,sctrm2,qli,qlj,lcspair,lewaldtype, &
                     .false.,.false.,.false.)
      endif
!*****************************************************************
!  Calculate reciprocal space contribution to third derivatives  *
!*****************************************************************
      if (lewald.and..not.lnorecip) then
        call reciptrmdp(xcrd,ycrd,zcrd,.true.,.false.,ofct,d2,d2s,d3,d3s,d3ss)
      endif
!****************************
!  Loop over all distances  *
!****************************
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
!  Third derivative terms for second derivatives
!
          d3(1,1,1) = d3(1,1,1) + derivqd3(kk)*d2r2dx2(kk,1)*xcrd + 3.0_dp*xcrd*derivqd2(kk)
          d3(2,1,1) = d3(2,1,1) + derivqd3(kk)*d2r2dx2(kk,1)*ycrd + ycrd*derivqd2(kk)
          d3(3,1,1) = d3(3,1,1) + derivqd3(kk)*d2r2dx2(kk,1)*zcrd + zcrd*derivqd2(kk)
          d3(1,2,1) = d3(1,2,1) + derivqd3(kk)*d2r2dx2(kk,1)*ycrd + ycrd*derivqd2(kk)
          d3(2,2,1) = d3(2,2,1) + derivqd3(kk)*d2r2dx2(kk,2)*xcrd + xcrd*derivqd2(kk)
          d3(3,2,1) = d3(3,2,1) + derivqd3(kk)*d2r2dx2(kk,4)*xcrd
          d3(1,3,1) = d3(1,3,1) + derivqd3(kk)*d2r2dx2(kk,1)*zcrd + zcrd*derivqd2(kk)
          d3(2,3,1) = d3(2,3,1) + derivqd3(kk)*d2r2dx2(kk,4)*xcrd
          d3(3,3,1) = d3(3,3,1) + derivqd3(kk)*d2r2dx2(kk,3)*xcrd + xcrd*derivqd2(kk)
!
          d3(1,1,2) = d3(1,1,2) + derivqd3(kk)*d2r2dx2(kk,1)*ycrd + ycrd*derivqd2(kk)
          d3(2,1,2) = d3(2,1,2) + derivqd3(kk)*d2r2dx2(kk,2)*xcrd + xcrd*derivqd2(kk)
          d3(3,1,2) = d3(3,1,2) + derivqd3(kk)*d2r2dx2(kk,5)*ycrd
          d3(1,2,2) = d3(1,2,2) + derivqd3(kk)*d2r2dx2(kk,2)*xcrd + xcrd*derivqd2(kk)
          d3(2,2,2) = d3(2,2,2) + derivqd3(kk)*d2r2dx2(kk,2)*ycrd + 3.0_dp*ycrd*derivqd2(kk)
          d3(3,2,2) = d3(3,2,2) + derivqd3(kk)*d2r2dx2(kk,2)*zcrd + zcrd*derivqd2(kk)
          d3(1,3,2) = d3(1,3,2) + derivqd3(kk)*d2r2dx2(kk,5)*ycrd
          d3(2,3,2) = d3(2,3,2) + derivqd3(kk)*d2r2dx2(kk,2)*zcrd + zcrd*derivqd2(kk)
          d3(3,3,2) = d3(3,3,2) + derivqd3(kk)*d2r2dx2(kk,3)*ycrd + ycrd*derivqd2(kk)
!
          d3(1,1,3) = d3(1,1,3) + derivqd3(kk)*d2r2dx2(kk,1)*zcrd + zcrd*derivqd2(kk)
          d3(2,1,3) = d3(2,1,3) + derivqd3(kk)*d2r2dx2(kk,6)*zcrd
          d3(3,1,3) = d3(3,1,3) + derivqd3(kk)*d2r2dx2(kk,3)*xcrd + xcrd*derivqd2(kk)
          d3(1,2,3) = d3(1,2,3) + derivqd3(kk)*d2r2dx2(kk,6)*zcrd
          d3(2,2,3) = d3(2,2,3) + derivqd3(kk)*d2r2dx2(kk,2)*zcrd + zcrd*derivqd2(kk)
          d3(3,2,3) = d3(3,2,3) + derivqd3(kk)*d2r2dx2(kk,3)*ycrd + ycrd*derivqd2(kk)
          d3(1,3,3) = d3(1,3,3) + derivqd3(kk)*d2r2dx2(kk,3)*xcrd + xcrd*derivqd2(kk)
          d3(2,3,3) = d3(2,3,3) + derivqd3(kk)*d2r2dx2(kk,3)*ycrd + ycrd*derivqd2(kk)
          d3(3,3,3) = d3(3,3,3) + derivqd3(kk)*d2r2dx2(kk,3)*zcrd + 3.0_dp*zcrd*derivqd2(kk)
!****************************
!  End loop over distances  *
!****************************
        enddo
      endif
!**************************
!  Coordinate Derivatives *
!**************************
      vxi = vx(i)
      vyi = vy(i)
      vzi = vz(i)
!
      vxj = vx(j)
      vyj = vy(j)
      vzj = vz(j)
!
!  Save d2 term
!
      ld2ok(j) = .true.
      d2a(1:3,1:3,j) = d2(1:3,1:3)
!
!  Second derivatives - excluding diagonal block
!
      if (j.ne.i) then
        derv2(jx,ix) = derv2(jx,ix) - poli*(vxi*d3(1,1,1) + vyi*d3(1,1,2) + vzi*d3(1,1,3) - &
                                            qlj*(d2(1,1)*d2(1,1) + d2(1,2)*d2(1,2) + d2(1,3)*d2(1,3)))
        derv2(jy,ix) = derv2(jy,ix) - poli*(vxi*d3(2,1,1) + vyi*d3(2,1,2) + vzi*d3(2,1,3) - &
                                            qlj*(d2(2,1)*d2(1,1) + d2(2,2)*d2(1,2) + d2(2,3)*d2(1,3)))
        derv2(jz,ix) = derv2(jz,ix) - poli*(vxi*d3(3,1,1) + vyi*d3(3,1,2) + vzi*d3(3,1,3) - &
                                            qlj*(d2(3,1)*d2(1,1) + d2(3,2)*d2(1,2) + d2(3,3)*d2(1,3)))
        derv2(jx,iy) = derv2(jx,iy) - poli*(vxi*d3(1,2,1) + vyi*d3(1,2,2) + vzi*d3(1,2,3) - &
                                            qlj*(d2(1,1)*d2(2,1) + d2(1,2)*d2(2,2) + d2(1,3)*d2(2,3)))
        derv2(jy,iy) = derv2(jy,iy) - poli*(vxi*d3(2,2,1) + vyi*d3(2,2,2) + vzi*d3(2,2,3) - &
                                            qlj*(d2(2,1)*d2(2,1) + d2(2,2)*d2(2,2) + d2(2,3)*d2(2,3)))
        derv2(jz,iy) = derv2(jz,iy) - poli*(vxi*d3(3,2,1) + vyi*d3(3,2,2) + vzi*d3(3,2,3) - &
                                            qlj*(d2(3,1)*d2(2,1) + d2(3,2)*d2(2,2) + d2(3,3)*d2(2,3)))
        derv2(jx,iz) = derv2(jx,iz) - poli*(vxi*d3(1,3,1) + vyi*d3(1,3,2) + vzi*d3(1,3,3) - &
                                            qlj*(d2(1,1)*d2(3,1) + d2(1,2)*d2(3,2) + d2(1,3)*d2(3,3)))
        derv2(jy,iz) = derv2(jy,iz) - poli*(vxi*d3(2,3,1) + vyi*d3(2,3,2) + vzi*d3(2,3,3) - &
                                            qlj*(d2(2,1)*d2(3,1) + d2(2,2)*d2(3,2) + d2(2,3)*d2(3,3)))
        derv2(jz,iz) = derv2(jz,iz) - poli*(vxi*d3(3,3,1) + vyi*d3(3,3,2) + vzi*d3(3,3,3) - &
                                            qlj*(d2(3,1)*d2(3,1) + d2(3,2)*d2(3,2) + d2(3,3)*d2(3,3)))
!
        if (lpollangevin) then
          derv2(jx,ix) = derv2(jx,ix) + pol2i*(vxi*d2(1,1) + vyi*d2(2,1) + vzi*d2(3,1))* &
                                              (vxi*d2(1,1) + vyi*d2(2,1) + vzi*d2(3,1))
          derv2(jy,ix) = derv2(jy,ix) + pol2i*(vxi*d2(1,1) + vyi*d2(2,1) + vzi*d2(3,1))* &
                                              (vxi*d2(1,2) + vyi*d2(2,2) + vzi*d2(3,2))
          derv2(jz,ix) = derv2(jz,ix) + pol2i*(vxi*d2(1,1) + vyi*d2(2,1) + vzi*d2(3,1))* &
                                              (vxi*d2(1,3) + vyi*d2(2,3) + vzi*d2(3,3))
          derv2(jx,iy) = derv2(jx,iy) + pol2i*(vxi*d2(1,2) + vyi*d2(2,2) + vzi*d2(3,2))* &
                                              (vxi*d2(1,1) + vyi*d2(2,1) + vzi*d2(3,1))
          derv2(jy,iy) = derv2(jy,iy) + pol2i*(vxi*d2(1,2) + vyi*d2(2,2) + vzi*d2(3,2))* &
                                              (vxi*d2(1,2) + vyi*d2(2,2) + vzi*d2(3,2))
          derv2(jz,iy) = derv2(jz,iy) + pol2i*(vxi*d2(1,2) + vyi*d2(2,2) + vzi*d2(3,2))* &
                                              (vxi*d2(1,3) + vyi*d2(2,3) + vzi*d2(3,3))
          derv2(jx,iz) = derv2(jx,iz) + pol2i*(vxi*d2(1,3) + vyi*d2(2,3) + vzi*d2(3,3))* &
                                              (vxi*d2(1,1) + vyi*d2(2,1) + vzi*d2(3,1))
          derv2(jy,iz) = derv2(jy,iz) + pol2i*(vxi*d2(1,3) + vyi*d2(2,3) + vzi*d2(3,3))* &
                                              (vxi*d2(1,2) + vyi*d2(2,2) + vzi*d2(3,2))
          derv2(jz,iz) = derv2(jz,iz) + pol2i*(vxi*d2(1,3) + vyi*d2(2,3) + vzi*d2(3,3))* &
                                              (vxi*d2(1,3) + vyi*d2(2,3) + vzi*d2(3,3))
        endif
!
        derv2(jx,ix) = derv2(jx,ix) + polj*(vxj*d3(1,1,1) + vyj*d3(1,1,2) + vzj*d3(1,1,3) + &
                                            qli*(d2(1,1)*d2(1,1) + d2(1,2)*d2(1,2) + d2(1,3)*d2(1,3)))
        derv2(jy,ix) = derv2(jy,ix) + polj*(vxj*d3(2,1,1) + vyj*d3(2,1,2) + vzj*d3(2,1,3) + &
                                            qli*(d2(2,1)*d2(1,1) + d2(2,2)*d2(1,2) + d2(2,3)*d2(1,3)))
        derv2(jz,ix) = derv2(jz,ix) + polj*(vxj*d3(3,1,1) + vyj*d3(3,1,2) + vzj*d3(3,1,3) + &
                                            qli*(d2(3,1)*d2(1,1) + d2(3,2)*d2(1,2) + d2(3,3)*d2(1,3)))
        derv2(jx,iy) = derv2(jx,iy) + polj*(vxj*d3(1,2,1) + vyj*d3(1,2,2) + vzj*d3(1,2,3) + &
                                            qli*(d2(1,1)*d2(2,1) + d2(1,2)*d2(2,2) + d2(1,3)*d2(2,3)))
        derv2(jy,iy) = derv2(jy,iy) + polj*(vxj*d3(2,2,1) + vyj*d3(2,2,2) + vzj*d3(2,2,3) + &
                                            qli*(d2(2,1)*d2(2,1) + d2(2,2)*d2(2,2) + d2(2,3)*d2(2,3)))
        derv2(jz,iy) = derv2(jz,iy) + polj*(vxj*d3(3,2,1) + vyj*d3(3,2,2) + vzj*d3(3,2,3) + &
                                            qli*(d2(3,1)*d2(2,1) + d2(3,2)*d2(2,2) + d2(3,3)*d2(2,3)))
        derv2(jx,iz) = derv2(jx,iz) + polj*(vxj*d3(1,3,1) + vyj*d3(1,3,2) + vzj*d3(1,3,3) + &
                                            qli*(d2(1,1)*d2(3,1) + d2(1,2)*d2(3,2) + d2(1,3)*d2(3,3)))
        derv2(jy,iz) = derv2(jy,iz) + polj*(vxj*d3(2,3,1) + vyj*d3(2,3,2) + vzj*d3(2,3,3) + &
                                            qli*(d2(2,1)*d2(3,1) + d2(2,2)*d2(3,2) + d2(2,3)*d2(3,3)))
        derv2(jz,iz) = derv2(jz,iz) + polj*(vxj*d3(3,3,1) + vyj*d3(3,3,2) + vzj*d3(3,3,3) + &
                                            qli*(d2(3,1)*d2(3,1) + d2(3,2)*d2(3,2) + d2(3,3)*d2(3,3)))
!
        if (lpollangevin) then
          derv2(jx,ix) = derv2(jx,ix) + pol2j*(vxj*d2(1,1) + vyj*d2(2,1) + vzj*d2(3,1))* &
                                              (vxj*d2(1,1) + vyj*d2(2,1) + vzj*d2(3,1))
          derv2(jy,ix) = derv2(jy,ix) + pol2j*(vxj*d2(1,1) + vyj*d2(2,1) + vzj*d2(3,1))* &
                                              (vxj*d2(1,2) + vyj*d2(2,2) + vzj*d2(3,2))
          derv2(jz,ix) = derv2(jz,ix) + pol2j*(vxj*d2(1,1) + vyj*d2(2,1) + vzj*d2(3,1))* &
                                              (vxj*d2(1,3) + vyj*d2(2,3) + vzj*d2(3,3))
          derv2(jx,iy) = derv2(jx,iy) + pol2j*(vxj*d2(1,2) + vyj*d2(2,2) + vzj*d2(3,2))* &
                                              (vxj*d2(1,1) + vyj*d2(2,1) + vzj*d2(3,1))
          derv2(jy,iy) = derv2(jy,iy) + pol2j*(vxj*d2(1,2) + vyj*d2(2,2) + vzj*d2(3,2))* &
                                              (vxj*d2(1,2) + vyj*d2(2,2) + vzj*d2(3,2))
          derv2(jz,iy) = derv2(jz,iy) + pol2j*(vxj*d2(1,2) + vyj*d2(2,2) + vzj*d2(3,2))* &
                                              (vxj*d2(1,3) + vyj*d2(2,3) + vzj*d2(3,3))
          derv2(jx,iz) = derv2(jx,iz) + pol2j*(vxj*d2(1,3) + vyj*d2(2,3) + vzj*d2(3,3))* &
                                              (vxj*d2(1,1) + vyj*d2(2,1) + vzj*d2(3,1))
          derv2(jy,iz) = derv2(jy,iz) + pol2j*(vxj*d2(1,3) + vyj*d2(2,3) + vzj*d2(3,3))* &
                                              (vxj*d2(1,2) + vyj*d2(2,2) + vzj*d2(3,2))
          derv2(jz,iz) = derv2(jz,iz) + pol2j*(vxj*d2(1,3) + vyj*d2(2,3) + vzj*d2(3,3))* &
                                              (vxj*d2(1,3) + vyj*d2(2,3) + vzj*d2(3,3))
        endif
      endif
    enddo jloop
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
!  If not ld2ok then skip
!
      if (.not.ld2ok(j)) cycle jloop2
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
!  If not ld2ok then skip
!
        if (.not.ld2ok(k)) cycle kloop
!
        qlk = qf(k)
        if (lpollangevin) then
          Es = sqrt(vx(i)**2 + vy(i)**2 + vz(i)**2)
          if (Es.gt.1.0d-8) then
            call langevinpol(Es,dpolarmax(nrelf2a(i)),dpolar(nrelf2a(i)),Ep,dEpdEs,d2EpdEs2,.true.,.true.)
            polk = qlj*qlk*dEpdEs/angstoev
            pol2k = qlj*qlk*d2EpdEs2/angstoev
          else
            polk = qlj*qlk*dpolar(nrelf2a(i))/angstoev
            pol2k = 0.0_dp
          endif
        else
          polk = qlj*qlk*dpolar(nrelf2a(i))/angstoev
          pol2k = 0.0_dp
        endif
!
        derv2(jx,ix) = derv2(jx,ix) + polk*(d2a(1,1,j)*d2a(1,1,k) + d2a(1,2,j)*d2a(1,2,k) + d2a(1,3,j)*d2a(1,3,k))
        derv2(jy,ix) = derv2(jy,ix) + polk*(d2a(2,1,j)*d2a(1,1,k) + d2a(2,2,j)*d2a(1,2,k) + d2a(2,3,j)*d2a(1,3,k))
        derv2(jz,ix) = derv2(jz,ix) + polk*(d2a(3,1,j)*d2a(1,1,k) + d2a(3,2,j)*d2a(1,2,k) + d2a(3,3,j)*d2a(1,3,k))
        derv2(jx,iy) = derv2(jx,iy) + polk*(d2a(1,1,j)*d2a(2,1,k) + d2a(1,2,j)*d2a(2,2,k) + d2a(1,3,j)*d2a(2,3,k))
        derv2(jy,iy) = derv2(jy,iy) + polk*(d2a(2,1,j)*d2a(2,1,k) + d2a(2,2,j)*d2a(2,2,k) + d2a(2,3,j)*d2a(2,3,k))
        derv2(jz,iy) = derv2(jz,iy) + polk*(d2a(3,1,j)*d2a(2,1,k) + d2a(3,2,j)*d2a(2,2,k) + d2a(3,3,j)*d2a(2,3,k))
        derv2(jx,iz) = derv2(jx,iz) + polk*(d2a(1,1,j)*d2a(3,1,k) + d2a(1,2,j)*d2a(3,2,k) + d2a(1,3,j)*d2a(3,3,k))
        derv2(jy,iz) = derv2(jy,iz) + polk*(d2a(2,1,j)*d2a(3,1,k) + d2a(2,2,j)*d2a(3,2,k) + d2a(2,3,j)*d2a(3,3,k))
        derv2(jz,iz) = derv2(jz,iz) + polk*(d2a(3,1,j)*d2a(3,1,k) + d2a(3,2,j)*d2a(3,2,k) + d2a(3,3,j)*d2a(3,3,k))

!
        derv2(kx,ix) = derv2(kx,ix) + polk*(d2a(1,1,k)*d2a(1,1,j) + d2a(1,2,k)*d2a(1,2,j) + d2a(1,3,k)*d2a(1,3,j))
        derv2(ky,ix) = derv2(ky,ix) + polk*(d2a(2,1,k)*d2a(1,1,j) + d2a(2,2,k)*d2a(1,2,j) + d2a(2,3,k)*d2a(1,3,j))
        derv2(kz,ix) = derv2(kz,ix) + polk*(d2a(3,1,k)*d2a(1,1,j) + d2a(3,2,k)*d2a(1,2,j) + d2a(3,3,k)*d2a(1,3,j))
        derv2(kx,iy) = derv2(kx,iy) + polk*(d2a(1,1,k)*d2a(2,1,j) + d2a(1,2,k)*d2a(2,2,j) + d2a(1,3,k)*d2a(2,3,j))
        derv2(ky,iy) = derv2(ky,iy) + polk*(d2a(2,1,k)*d2a(2,1,j) + d2a(2,2,k)*d2a(2,2,j) + d2a(2,3,k)*d2a(2,3,j))
        derv2(kz,iy) = derv2(kz,iy) + polk*(d2a(3,1,k)*d2a(2,1,j) + d2a(3,2,k)*d2a(2,2,j) + d2a(3,3,k)*d2a(2,3,j))
        derv2(kx,iz) = derv2(kx,iz) + polk*(d2a(1,1,k)*d2a(3,1,j) + d2a(1,2,k)*d2a(3,2,j) + d2a(1,3,k)*d2a(3,3,j))
        derv2(ky,iz) = derv2(ky,iz) + polk*(d2a(2,1,k)*d2a(3,1,j) + d2a(2,2,k)*d2a(3,2,j) + d2a(2,3,k)*d2a(3,3,j))
        derv2(kz,iz) = derv2(kz,iz) + polk*(d2a(3,1,k)*d2a(3,1,j) + d2a(3,2,k)*d2a(3,2,j) + d2a(3,3,k)*d2a(3,3,j))
!
        if (lpollangevin) then
          derv2(jx,ix) = derv2(jx,ix) + pol2k*(vx(i)*d2a(1,1,j) + vy(i)*d2a(2,1,j) + vz(i)*d2a(3,1,j))* &
                                              (vx(i)*d2a(1,1,k) + vy(i)*d2a(2,1,k) + vz(i)*d2a(3,1,k))
          derv2(jx,iy) = derv2(jx,iy) + pol2k*(vx(i)*d2a(1,1,j) + vy(i)*d2a(2,1,j) + vz(i)*d2a(3,1,j))* &
                                              (vx(i)*d2a(1,2,k) + vy(i)*d2a(2,2,k) + vz(i)*d2a(3,2,k))
          derv2(jx,iz) = derv2(jx,iz) + pol2k*(vx(i)*d2a(1,1,j) + vy(i)*d2a(2,1,j) + vz(i)*d2a(3,1,j))* &
                                              (vx(i)*d2a(1,3,k) + vy(i)*d2a(2,3,k) + vz(i)*d2a(3,3,k))
          derv2(jy,ix) = derv2(jy,ix) + pol2k*(vx(i)*d2a(1,2,j) + vy(i)*d2a(2,2,j) + vz(i)*d2a(3,2,j))* &
                                              (vx(i)*d2a(1,1,k) + vy(i)*d2a(2,1,k) + vz(i)*d2a(3,1,k))
          derv2(jy,iy) = derv2(jy,iy) + pol2k*(vx(i)*d2a(1,2,j) + vy(i)*d2a(2,2,j) + vz(i)*d2a(3,2,j))* &
                                              (vx(i)*d2a(1,2,k) + vy(i)*d2a(2,2,k) + vz(i)*d2a(3,2,k))
          derv2(jy,iz) = derv2(jy,iz) + pol2k*(vx(i)*d2a(1,2,j) + vy(i)*d2a(2,2,j) + vz(i)*d2a(3,2,j))* &
                                              (vx(i)*d2a(1,3,k) + vy(i)*d2a(2,3,k) + vz(i)*d2a(3,3,k))
          derv2(jz,ix) = derv2(jz,ix) + pol2k*(vx(i)*d2a(1,3,j) + vy(i)*d2a(2,3,j) + vz(i)*d2a(3,3,j))* &
                                              (vx(i)*d2a(1,1,k) + vy(i)*d2a(2,1,k) + vz(i)*d2a(3,1,k))
          derv2(jz,iy) = derv2(jz,iy) + pol2k*(vx(i)*d2a(1,3,j) + vy(i)*d2a(2,3,j) + vz(i)*d2a(3,3,j))* &
                                              (vx(i)*d2a(1,2,k) + vy(i)*d2a(2,2,k) + vz(i)*d2a(3,2,k))
          derv2(jz,iz) = derv2(jz,iz) + pol2k*(vx(i)*d2a(1,3,j) + vy(i)*d2a(2,3,j) + vz(i)*d2a(3,3,j))* &
                                              (vx(i)*d2a(1,3,k) + vy(i)*d2a(2,3,k) + vz(i)*d2a(3,3,k))
!
          derv2(kx,ix) = derv2(kx,ix) + pol2k*(vx(i)*d2a(1,1,k) + vy(i)*d2a(2,1,k) + vz(i)*d2a(3,1,k))* &
                                              (vx(i)*d2a(1,1,j) + vy(i)*d2a(2,1,j) + vz(i)*d2a(3,1,j))
          derv2(kx,iy) = derv2(kx,iy) + pol2k*(vx(i)*d2a(1,1,k) + vy(i)*d2a(2,1,k) + vz(i)*d2a(3,1,k))* &
                                              (vx(i)*d2a(1,2,j) + vy(i)*d2a(2,2,j) + vz(i)*d2a(3,2,j))
          derv2(kx,iz) = derv2(kx,iz) + pol2k*(vx(i)*d2a(1,1,k) + vy(i)*d2a(2,1,k) + vz(i)*d2a(3,1,k))* &
                                              (vx(i)*d2a(1,3,j) + vy(i)*d2a(2,3,j) + vz(i)*d2a(3,3,j))
          derv2(ky,ix) = derv2(ky,ix) + pol2k*(vx(i)*d2a(1,2,k) + vy(i)*d2a(2,2,k) + vz(i)*d2a(3,2,k))* &
                                              (vx(i)*d2a(1,1,j) + vy(i)*d2a(2,1,j) + vz(i)*d2a(3,1,j))
          derv2(ky,iy) = derv2(ky,iy) + pol2k*(vx(i)*d2a(1,2,k) + vy(i)*d2a(2,2,k) + vz(i)*d2a(3,2,k))* &
                                              (vx(i)*d2a(1,2,j) + vy(i)*d2a(2,2,j) + vz(i)*d2a(3,2,j))
          derv2(ky,iz) = derv2(ky,iz) + pol2k*(vx(i)*d2a(1,2,k) + vy(i)*d2a(2,2,k) + vz(i)*d2a(3,2,k))* &
                                              (vx(i)*d2a(1,3,j) + vy(i)*d2a(2,3,j) + vz(i)*d2a(3,3,j))
          derv2(kz,ix) = derv2(kz,ix) + pol2k*(vx(i)*d2a(1,3,k) + vy(i)*d2a(2,3,k) + vz(i)*d2a(3,3,k))* &
                                              (vx(i)*d2a(1,1,j) + vy(i)*d2a(2,1,j) + vz(i)*d2a(3,1,j))
          derv2(kz,iy) = derv2(kz,iy) + pol2k*(vx(i)*d2a(1,3,k) + vy(i)*d2a(2,3,k) + vz(i)*d2a(3,3,k))* &
                                              (vx(i)*d2a(1,2,j) + vy(i)*d2a(2,2,j) + vz(i)*d2a(3,2,j))
          derv2(kz,iz) = derv2(kz,iz) + pol2k*(vx(i)*d2a(1,3,k) + vy(i)*d2a(2,3,k) + vz(i)*d2a(3,3,k))* &
                                              (vx(i)*d2a(1,3,j) + vy(i)*d2a(2,3,j) + vz(i)*d2a(3,3,j))
        endif
      enddo kloop
    enddo jloop2
  enddo iloop
!************************************************************************************
!  Outer loop over sites with polarisability on site i for second derivative terms  *
!************************************************************************************
  ix = -2
  iy = -1
  iz =  0
!
  iloop2: do i = 1,numat
    ix = ix + 3
    iy = iy + 3
    iz = iz + 3
!
!  If i has no polarisability then skip this atom
!
    if (dpolar(nrelf2a(i)).eq.0.0_dp) cycle iloop2
!
!  Initialise d2a arrays
!
    d2a(1:3,1:3,1:numat) = 0.0_dp
!
!  Inner loop over second site
!
    xal = xclat(i)
    yal = yclat(i)
    zal = zclat(i)
    nati = nat(i)
    ntypi = nftype(i)
    qli = qf(i)
    oci = occuf(i)
!
!  Molecule handling
!
    if (lmol) then
      nmi = natmol(i)
      indm = nmolind(i)
      call mindtoijk(indm,ixi,iyi,izi)
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
!  Start of second atom loop
!
    jx = -2
    jy = -1
    jz =  0
!
!  Inner loop over second site to generate j terms for atoms on this node
!
    jloop3: do jonnode = 1,natomsonnode
      j = node2atom(jonnode)
!
      jx = jx + 3
      jy = jy + 3
      jz = jz + 3
!
      qlj = qf(j)
!
!  Compute charge x polarisability for pair
!
      if (lpollangevin) then
        Es = sqrt(vx(i)**2 + vy(i)**2 + vz(i)**2)
        if (Es.gt.1.0d-8) then
          call langevinpol(Es,dpolarmax(nrelf2a(i)),dpolar(nrelf2a(i)),Ep,dEpdEs,d2EpdEs2,.true.,.true.)
          poli = qlj*dEpdEs/angstoev
          pol2i = qlj*qlj*d2EpdEs2/angstoev
        else
          poli = qlj*dpolar(nrelf2a(i))/angstoev
          pol2i = 0.0_dp
        endif
      else
        poli = qlj*dpolar(nrelf2a(i))/angstoev
        pol2i = 0.0_dp
      endif
!
!  If there will be no interaction then skip
!
      if (abs(poli).lt.1.0d-12) cycle jloop3
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
!
!  Zero third derivative arrays
!
      d3(1:3,1:3,1:3) = 0.0_dp
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
        call twobody(eatom,ereal,ec6,.true.,.true.,.true.,nor,1,npots,npotl,cut2,cut2q,cut2s, &
                     nmolonly,factor,ofct,ospfct,0.0_dp,sctrm1,sctrm2,qli,qlj,lcspair,lewaldtype, &
                     .false.,.false.,.false.)
      endif
!*****************************************************************
!  Calculate reciprocal space contribution to third derivatives  *
!*****************************************************************
      if (lewald.and..not.lnorecip) then
        call reciptrmdp(xcrd,ycrd,zcrd,.true.,.false.,ofct,d2,d2s,d3,d3s,d3ss)
      endif
!****************************
!  Loop over all distances  *
!****************************
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
!****************************
!  End loop over distances  *
!****************************
        enddo
      endif
!
!  Save d2 term
!  
      d2a(1:3,1:3,j) = d2(1:3,1:3)
    enddo jloop3
!***************************************
!  Now globalise d2k across all nodes  *
!***************************************
    call sumall(d2a,d2asum,9_i4*numat,"pirealrecippd","d2a")
    d2a(1:3,1:3,1:numat) = d2asum(1:3,1:3,1:numat)
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
        if (k.eq.i) cycle kloop3
        if (k.eq.j) cycle kloop3
!
        qlk = qf(k)
        if (lpollangevin) then
          Es = sqrt(vx(i)**2 + vy(i)**2 + vz(i)**2)
          if (Es.gt.1.0d-8) then
            call langevinpol(Es,dpolarmax(nrelf2a(i)),dpolar(nrelf2a(i)),Ep,dEpdEs,d2EpdEs2,.true.,.true.)
            polk = qlj*qlk*dEpdEs/angstoev
            pol2k = qlj*qlk*d2EpdEs2/angstoev
          else
            polk = qlj*qlk*dpolar(nrelf2a(i))/angstoev
            pol2k = 0.0_dp
          endif
        else
          polk = qlj*qlk*dpolar(nrelf2a(i))/angstoev
          pol2k = 0.0_dp
        endif
!
        derv2(ix,jx) = derv2(ix,jx) + polk*(d2a(1,1,k)*d2a(1,1,j) + d2a(1,2,k)*d2a(1,2,j) + d2a(1,3,k)*d2a(1,3,j))
        derv2(iy,jx) = derv2(iy,jx) + polk*(d2a(2,1,k)*d2a(1,1,j) + d2a(2,2,k)*d2a(1,2,j) + d2a(2,3,k)*d2a(1,3,j))
        derv2(iz,jx) = derv2(iz,jx) + polk*(d2a(3,1,k)*d2a(1,1,j) + d2a(3,2,k)*d2a(1,2,j) + d2a(3,3,k)*d2a(1,3,j))
        derv2(ix,jy) = derv2(ix,jy) + polk*(d2a(1,1,k)*d2a(2,1,j) + d2a(1,2,k)*d2a(2,2,j) + d2a(1,3,k)*d2a(2,3,j))
        derv2(iy,jy) = derv2(iy,jy) + polk*(d2a(2,1,k)*d2a(2,1,j) + d2a(2,2,k)*d2a(2,2,j) + d2a(2,3,k)*d2a(2,3,j))
        derv2(iz,jy) = derv2(iz,jy) + polk*(d2a(3,1,k)*d2a(2,1,j) + d2a(3,2,k)*d2a(2,2,j) + d2a(3,3,k)*d2a(2,3,j))
        derv2(ix,jz) = derv2(ix,jz) + polk*(d2a(1,1,k)*d2a(3,1,j) + d2a(1,2,k)*d2a(3,2,j) + d2a(1,3,k)*d2a(3,3,j))
        derv2(iy,jz) = derv2(iy,jz) + polk*(d2a(2,1,k)*d2a(3,1,j) + d2a(2,2,k)*d2a(3,2,j) + d2a(2,3,k)*d2a(3,3,j))
        derv2(iz,jz) = derv2(iz,jz) + polk*(d2a(3,1,k)*d2a(3,1,j) + d2a(3,2,k)*d2a(3,2,j) + d2a(3,3,k)*d2a(3,3,j))
!
        if (lpollangevin) then
          derv2(ix,jx) = derv2(ix,jx) + pol2k*(vx(i)*d2a(1,1,k) + vy(i)*d2a(2,1,k) + vz(i)*d2a(3,1,k))* &
                                              (vx(i)*d2a(1,1,j) + vy(i)*d2a(2,1,j) + vz(i)*d2a(3,1,j))
          derv2(ix,jy) = derv2(ix,jy) + pol2k*(vx(i)*d2a(1,1,k) + vy(i)*d2a(2,1,k) + vz(i)*d2a(3,1,k))* &
                                              (vx(i)*d2a(1,2,j) + vy(i)*d2a(2,2,j) + vz(i)*d2a(3,2,j))
          derv2(ix,jz) = derv2(ix,jz) + pol2k*(vx(i)*d2a(1,1,k) + vy(i)*d2a(2,1,k) + vz(i)*d2a(3,1,k))* &
                                              (vx(i)*d2a(1,3,j) + vy(i)*d2a(2,3,j) + vz(i)*d2a(3,3,j))
          derv2(iy,jx) = derv2(iy,jx) + pol2k*(vx(i)*d2a(1,2,k) + vy(i)*d2a(2,2,k) + vz(i)*d2a(3,2,k))* &
                                              (vx(i)*d2a(1,1,j) + vy(i)*d2a(2,1,j) + vz(i)*d2a(3,1,j))
          derv2(iy,jy) = derv2(iy,jy) + pol2k*(vx(i)*d2a(1,2,k) + vy(i)*d2a(2,2,k) + vz(i)*d2a(3,2,k))* &
                                              (vx(i)*d2a(1,2,j) + vy(i)*d2a(2,2,j) + vz(i)*d2a(3,2,j))
          derv2(iy,jz) = derv2(iy,jz) + pol2k*(vx(i)*d2a(1,2,k) + vy(i)*d2a(2,2,k) + vz(i)*d2a(3,2,k))* &
                                              (vx(i)*d2a(1,3,j) + vy(i)*d2a(2,3,j) + vz(i)*d2a(3,3,j))
          derv2(iz,jx) = derv2(iz,jx) + pol2k*(vx(i)*d2a(1,3,k) + vy(i)*d2a(2,3,k) + vz(i)*d2a(3,3,k))* &
                                              (vx(i)*d2a(1,1,j) + vy(i)*d2a(2,1,j) + vz(i)*d2a(3,1,j))
          derv2(iz,jy) = derv2(iz,jy) + pol2k*(vx(i)*d2a(1,3,k) + vy(i)*d2a(2,3,k) + vz(i)*d2a(3,3,k))* &
                                              (vx(i)*d2a(1,2,j) + vy(i)*d2a(2,2,j) + vz(i)*d2a(3,2,j))
          derv2(iz,jz) = derv2(iz,jz) + pol2k*(vx(i)*d2a(1,3,k) + vy(i)*d2a(2,3,k) + vz(i)*d2a(3,3,k))* &
                                              (vx(i)*d2a(1,3,j) + vy(i)*d2a(2,3,j) + vz(i)*d2a(3,3,j))
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
            call langevinpol(Es,dpolarmax(nrelf2a(i)),dpolar(nrelf2a(i)),Ep,dEpdEs,d2EpdEs2,.true.,.true.)
            polk = qlj*qlk*dEpdEs/angstoev
            pol2k = qlj*qlk*d2EpdEs2/angstoev
          else
            polk = qlj*qlk*dpolar(nrelf2a(i))/angstoev
            pol2k = 0.0_dp
          endif
        else
          polk = qlj*qlk*dpolar(nrelf2a(i))/angstoev
          pol2k = 0.0_dp
        endif
!
        derv2(kx,jx) = derv2(kx,jx) - polk*(d2a(1,1,k)*d2a(1,1,j) + d2a(1,2,k)*d2a(1,2,j) + d2a(1,3,k)*d2a(1,3,j))
        derv2(ky,jx) = derv2(ky,jx) - polk*(d2a(2,1,k)*d2a(1,1,j) + d2a(2,2,k)*d2a(1,2,j) + d2a(2,3,k)*d2a(1,3,j))
        derv2(kz,jx) = derv2(kz,jx) - polk*(d2a(3,1,k)*d2a(1,1,j) + d2a(3,2,k)*d2a(1,2,j) + d2a(3,3,k)*d2a(1,3,j))
        derv2(kx,jy) = derv2(kx,jy) - polk*(d2a(1,1,k)*d2a(2,1,j) + d2a(1,2,k)*d2a(2,2,j) + d2a(1,3,k)*d2a(2,3,j))
        derv2(ky,jy) = derv2(ky,jy) - polk*(d2a(2,1,k)*d2a(2,1,j) + d2a(2,2,k)*d2a(2,2,j) + d2a(2,3,k)*d2a(2,3,j))
        derv2(kz,jy) = derv2(kz,jy) - polk*(d2a(3,1,k)*d2a(2,1,j) + d2a(3,2,k)*d2a(2,2,j) + d2a(3,3,k)*d2a(2,3,j))
        derv2(kx,jz) = derv2(kx,jz) - polk*(d2a(1,1,k)*d2a(3,1,j) + d2a(1,2,k)*d2a(3,2,j) + d2a(1,3,k)*d2a(3,3,j))
        derv2(ky,jz) = derv2(ky,jz) - polk*(d2a(2,1,k)*d2a(3,1,j) + d2a(2,2,k)*d2a(3,2,j) + d2a(2,3,k)*d2a(3,3,j))
        derv2(kz,jz) = derv2(kz,jz) - polk*(d2a(3,1,k)*d2a(3,1,j) + d2a(3,2,k)*d2a(3,2,j) + d2a(3,3,k)*d2a(3,3,j))
!
        if (lpollangevin) then
          derv2(kx,jx) = derv2(kx,jx) - pol2k*(vx(i)*d2a(1,1,k) + vy(i)*d2a(2,1,k) + vz(i)*d2a(3,1,k))* &
                                              (vx(i)*d2a(1,1,j) + vy(i)*d2a(2,1,j) + vz(i)*d2a(3,1,j))
          derv2(kx,jy) = derv2(kx,jy) - pol2k*(vx(i)*d2a(1,1,k) + vy(i)*d2a(2,1,k) + vz(i)*d2a(3,1,k))* &
                                              (vx(i)*d2a(1,2,j) + vy(i)*d2a(2,2,j) + vz(i)*d2a(3,2,j))
          derv2(kx,jz) = derv2(kx,jz) - pol2k*(vx(i)*d2a(1,1,k) + vy(i)*d2a(2,1,k) + vz(i)*d2a(3,1,k))* &
                                              (vx(i)*d2a(1,3,j) + vy(i)*d2a(2,3,j) + vz(i)*d2a(3,3,j))
          derv2(ky,jx) = derv2(ky,jx) - pol2k*(vx(i)*d2a(1,2,k) + vy(i)*d2a(2,2,k) + vz(i)*d2a(3,2,k))* &
                                              (vx(i)*d2a(1,1,j) + vy(i)*d2a(2,1,j) + vz(i)*d2a(3,1,j))
          derv2(ky,jy) = derv2(ky,jy) - pol2k*(vx(i)*d2a(1,2,k) + vy(i)*d2a(2,2,k) + vz(i)*d2a(3,2,k))* &
                                              (vx(i)*d2a(1,2,j) + vy(i)*d2a(2,2,j) + vz(i)*d2a(3,2,j))
          derv2(ky,jz) = derv2(ky,jz) - pol2k*(vx(i)*d2a(1,2,k) + vy(i)*d2a(2,2,k) + vz(i)*d2a(3,2,k))* &
                                              (vx(i)*d2a(1,3,j) + vy(i)*d2a(2,3,j) + vz(i)*d2a(3,3,j))
          derv2(kz,jx) = derv2(kz,jx) - pol2k*(vx(i)*d2a(1,3,k) + vy(i)*d2a(2,3,k) + vz(i)*d2a(3,3,k))* &
                                              (vx(i)*d2a(1,1,j) + vy(i)*d2a(2,1,j) + vz(i)*d2a(3,1,j))
          derv2(kz,jy) = derv2(kz,jy) - pol2k*(vx(i)*d2a(1,3,k) + vy(i)*d2a(2,3,k) + vz(i)*d2a(3,3,k))* &
                                              (vx(i)*d2a(1,2,j) + vy(i)*d2a(2,2,j) + vz(i)*d2a(3,2,j))
          derv2(kz,jz) = derv2(kz,jz) - pol2k*(vx(i)*d2a(1,3,k) + vy(i)*d2a(2,3,k) + vz(i)*d2a(3,3,k))* &
                                              (vx(i)*d2a(1,3,j) + vy(i)*d2a(2,3,j) + vz(i)*d2a(3,3,j))
        endif
      enddo kloop4
    enddo jloop4
  enddo iloop2
!
!  Free local memory
!
  deallocate(d2asum,stat=status)
  if (status/=0) call deallocate_error('pirealrecippd','d2asum')
  deallocate(d2a,stat=status)
  if (status/=0) call deallocate_error('pirealrecippd','d2a')
  deallocate(npotl,stat=status)
  if (status/=0) call deallocate_error('pirealrecippd','npotl')
  deallocate(ld2ok,stat=status)
  if (status/=0) call deallocate_error('pirealrecippd','ld2ok')
!
!  Timing
!
  time2 = g_cpu_time()
  tpolar = tpolar + time2 - time1
#ifdef TRACE
  call trace_out('pirealrecippd')
#endif
!
  return
  end
