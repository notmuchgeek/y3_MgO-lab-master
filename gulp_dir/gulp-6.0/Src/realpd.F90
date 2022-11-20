  subroutine realpd(xkv,ykv,zkv)
!
!  Routine for calculating the phased second derivatives for use in a phonon calculation.
!  Distributed memory version.
!
!   7/12 Created from realp
!  12/14 rtrm1 changed from scalar to array
!   9/16 cputime renamed to g_cpu_time
!   9/16 constants renamed to g_constants
!   9/16 lorder12loc removed as it is not needed
!  12/16 derv2dk added back
!   2/18 Trace added
!   5/18 Modified for q ranges
!   8/19 Modifications for Intel compiler
!   8/19 Trap for neemrptr being zero added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!   3/20 Electrostatic cutoff only included where the either charge exceeds the threshold
!   4/20 Rigid molecule modifications added
!   5/20 Phasing based on centre of mass added
!   5/20 Phasing of group velocities modified for rigid molecules
!   6/20 Shells excluded from centre of mass phasing
!   6/20 Call to selfterm replaced with selftermp
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
!  Julian Gale, CIC, Curtin University, June 2020
!
  use configurations, only : lbsmat
  use g_constants
  use control
  use current
  use datatypes
  use derivatives
  use eemdata
  use element
  use general,        only : cutw
  use kspace
  use m_strain,       only : twostrterms
  use molecule
  use parallel
  use realvectors
  use shells
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
  real(dp),    intent(in)                      :: xkv
  real(dp),    intent(in)                      :: ykv
  real(dp),    intent(in)                      :: zkv
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ifs
  integer(i4)                                  :: iis
  integer(i4)                                  :: iloc
  integer(i4)                                  :: indi
  integer(i4)                                  :: indj
  integer(i4)                                  :: indm
  integer(i4)                                  :: indmj
  integer(i4)                                  :: indri
  integer(i4)                                  :: indrif
  integer(i4)                                  :: ix
  integer(i4)                                  :: iy
  integer(i4)                                  :: iz
  integer(i4)                                  :: ixf
  integer(i4)                                  :: iyf
  integer(i4)                                  :: izf
  integer(i4)                                  :: ixi
  integer(i4)                                  :: iyi
  integer(i4)                                  :: izi
  integer(i4)                                  :: ixj
  integer(i4)                                  :: iyj
  integer(i4)                                  :: izj
  integer(i4)                                  :: j
  integer(i4)                                  :: jx
  integer(i4)                                  :: jy
  integer(i4)                                  :: jz
  integer(i4)                                  :: k
  integer(i4)                                  :: maxlim
  integer(i4)                                  :: mint
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
  integer(i4)                                  :: nqri
  integer(i4)                                  :: nqrj
  integer(i4)                                  :: ntyp1
  integer(i4)                                  :: ntyp2
  integer(i4)                                  :: ntypi
  integer(i4)                                  :: ntypj
  integer(i4)                                  :: status
  logical                                      :: lmatch
  logical                                      :: lc6loc
  logical                                      :: lcspair
  logical                                      :: lewaldtype
  logical                                      :: lmolok
  logical                                      :: lneedmol
  logical                                      :: lself
  complex(dpc)                                 :: cdk(3)
  real(dp)                                     :: c6tot
  real(dp)                                     :: cosk
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: cut2
  real(dp)                                     :: cut2e
  real(dp)                                     :: cut2p
  real(dp)                                     :: cut2q
  real(dp)                                     :: cut2r
  real(dp)                                     :: cut2s
  real(dp)                                     :: cut2w
  real(dp)                                     :: d1ix
  real(dp)                                     :: d1iy
  real(dp)                                     :: d1iz
  real(dp)                                     :: d1jx
  real(dp)                                     :: d1jy
  real(dp)                                     :: d1jz
  real(dp)                                     :: d2k
  real(dp)                                     :: d2ks
  real(dp)                                     :: d2self
  real(dp)                                     :: derive0self
  real(dp)                                     :: derive0selfi
  real(dp)                                     :: derive0selfj
  real(dp)                                     :: dk
  real(dp)                                     :: dks
  real(dp)                                     :: eatom
  real(dp)                                     :: ec6
  real(dp)                                     :: eqeq  
  real(dp)                                     :: ereal
  real(dp)                                     :: erecip
  real(dp)                                     :: factor
  real(dp)                                     :: fct
  real(dp)                                     :: hfactor
  real(dp)                                     :: oci
  real(dp)                                     :: ocj
  real(dp)                                     :: ofct
  real(dp)                                     :: ospfct
  real(dp)                                     :: oneij
  real(dp)                                     :: qli
  real(dp)                                     :: qlj
  real(dp)                                     :: radi
  real(dp)                                     :: radj
  real(dp)                                     :: radsum
  real(dp)                                     :: ritrm
  real(dp)                                     :: rp
  real(dp)                                     :: rqeq2
  real(dp)                                     :: rrtrm
  real(dp)                                     :: sctrm1
  real(dp)                                     :: sctrm2
  real(dp)                                     :: sink
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
#ifdef TRACE
  call trace_in('realpd')
#endif
!
  time1 = g_cpu_time()
!
!  Local variables
!
  lc6loc = (lc6.and.ndim.eq.3)
  rqeq2 = rqeq*rqeq
  mint = 3*numat
  maxlim = mint
  if (nbsmat.gt.0) maxlim = maxlim + numat
  eatom = 0.0_dp
  ereal = 0.0_dp
  erecip = 0.0_dp
  ec6 = 0.0_dp
  eqeq = 0.0_dp
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
!  Set the Coulomb term type based on dimensionality :
!
!  1-D         => 1/r
!  2-D and 3-D => erfc(seta*r)/r
!
  lewaldtype = (ndim.ne.1.or.lwolf)
!
!  Allocate local memory
!
  allocate(npotl(npote),stat=status)
  if (status/=0) call outofmemory('realp','npotl')
!
  if (.not.lnoreal) then
!***************************************************************
!  Atomistic and real space electrostatic second derivatives   *
!***************************************************************
!
!  Outer loop over sites
!
    do iloc = 1,natomsonnode
      i = node2atom(iloc)
      xal = xclat(i)
      yal = yclat(i)
      zal = zclat(i)
      nati = nat(i)
      ntypi = nftype(i)
      qli = qf(i)
      oci = occuf(i)
      if (lbsmat(nsft+nrelf2a(i))) then
        radi = radf(i)
        indri  = 3*natomsonnode + iloc
        indrif = 3*numat + i
        indi = 3*(i - 1)
        ixf = indi + 1
        iyf = indi + 2
        izf = indi + 3
      else
        radi = 0.0_dp
      endif
      if (leem.and.lmultiqrange.and.neemrptr(i).ne.0) then
        nqri = nqrnow(neemrptr(i))
      else
        nqri = 1
      endif
      indi = 3*(iloc - 1)
      ix = indi + 1
      iy = indi + 2
      iz = indi + 3
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
        if (lrigid.and.nmi.gt.0.and.nati.le.maxele) then
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
      do j = 1,numat
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
        qlj = qf(j)
        ocj = occuf(j)
        if (lbsmat(nsft+nrelf2a(j))) then
          radj = radf(j)
        else
          radj = 0.0_dp
        endif
        radsum = radi + radj
        ofct = oci*ocj
        indj = 3*(j - 1)
        jx = indj + 1
        jy = indj + 2
        jz = indj + 3
        fct = ofct*angstoev
        factor = qli*qlj*fct
        if (leem.and.lmultiqrange.and.neemrptr(j).ne.0) then
          nqrj = nqrnow(neemrptr(j))
        else
          nqrj = 1
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
          if (lrigid.and.nmj.gt.0.and.natj.le.maxele) then
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
        c6tot = 0.0_dp
        lneedmol = (lmol.and..not.lmolq)
        do n = 1,npote
          if (lmatch(nat1,ntyp1,nspec1(n),nptyp1(n),.true.)) then
            if (lmatch(nat2,ntyp2,nspec2(n),nptyp2(n),.true.)) then
              npots = npots + 1
              npotl(npots) = n
              if (mmexc(n).gt.0.or.(lintra(n).and..not.linter(n).or..not.lintra(n).and.linter(n))) lneedmol = .true.
              if (nptype(n).eq.8.or.nptype(n).eq.33) then
                if (cuts.gt.rp) rp = cuts
              elseif (lc6loc) then
                if (nptype(n).eq.1.or.nptype(n).eq.7) then
                  c6tot = c6tot + twopot(3,n)
                  if (repcut(n).gt.rp) rp = repcut(n)
                elseif (nptype(n).eq.2.or.nptype(n).eq.21) then
                  c6tot = c6tot + twopot(2,n)
                  if (repcut(n).gt.rp) rp = repcut(n)
                else
                  if (rpot(n).gt.rp) rp = rpot(n)
                endif
              else
                if (rpot(n).gt.rp) rp = rpot(n)
              endif
            endif
          endif
        enddo
!
!  If no valid potentials and charge product is zero
!  then no need to search for distances
!
        if (npots.eq.0.and.abs(factor).lt.1.0d-8) goto 1120
        cut2r = rp*rp
        if (cut2r.gt.cut2p) cut2r = cut2p
        cut2 = cut2r
!
!  If both charges are less than threshold then exclude electrostatics from cutoff
!
        if (abs(qli*oci)+abs(qlj*ocj).gt.thresh_q) then
          if (cut2e.gt.cut2.and.lewald) cut2 = cut2e
          if (cut2w.gt.cut2.and.lwolf) cut2 = cut2w
        endif
        if (lqeq.or.lSandM) cut2 = max(cut2,rqeq2)
!  
!  If the Coulomb term and the potentials don't need molecule info for
!  for this pair of atoms, then turn off molecule checking to save time
!
        if (.not.lneedmol) lmolok = .false.
!
        d2self = 0.0_dp
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
        derive0self = 0.0_dp
        if (lself) call selftermp(xkv,ykv,zkv,xcom,ycom,zcom,derive0self,factor,fct,ofct,ospfct,1.0_dp,npotl,npots,c6tot, &
                                 d2self,i,j,ix,jx,1.0_dp,0.5_dp,lewaldtype,qli,qlj)
!
        if (nor.eq.0) goto 1110
!
        do k = 1,nor
          deriv2(k) = 0.0_dp
        enddo
!
!  Sqrt distances
!
        do k = 1,nor
          dist(k) = sqrt(dist(k))
        enddo
!
        call twobody(eatom,ereal,ec6,.true.,.true.,.false.,nor,1_i4,npots,npotl,cut2r,cut2q,cut2s, &
                     nmolonly,factor,ofct,ospfct,radsum,sctrm1,sctrm2,qli,qlj,lcspair,lewaldtype, &
                     .false.,.false.,.false.)
        if (leem.or.lDoQDeriv2) then
          do k = 1,nor
            d0i(k) = d0i(k) + derive0(k)*qlj
            d0j(k) = d0j(k) + derive0(k)*qli
          enddo
        endif
        if (leem) then
          if (lqeq) then
            call qeqbody(eqeq,.true.,.true.,nor,1_i4,fct,qli,qlj,nati,natj,nqri,nqrj)
          elseif (lSandM) then
            call smbody(eqeq,.true.,.true.,nor,1_i4,fct,qli,qlj,nati,natj,nqri,nqrj)
          endif
        endif
!
!  Generate products for derivatives
!
        call twostrterms(ndim,maxdis,nor,xtmp,ytmp,ztmp,xcom,ycom,zcom,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.true.)
!***************************************
!  Generate phased second derivatives  *
!***************************************
!
!  Set diagonal blocks to be difference from gamma point
!
        if (i.eq.j) then
          oneij = 1.0_dp
        else
          oneij = 0.0_dp
        endif
        do k = 1,nor
          if (lphasecom) then
            cosk = xkv*(xtmp(k) - xcom) + ykv*(ytmp(k) - ycom) + zkv*(ztmp(k) - zcom)
          else
            cosk = xkv*xtmp(k) + ykv*ytmp(k) + zkv*ztmp(k)
          endif
          sink = sin(cosk)
          cosk = cos(cosk) - oneij
          dk = deriv(k)*cosk
          d2k = deriv2(k)*cosk
          dks = deriv(k)*sink
          d2ks = deriv2(k)*sink
!
          derv2(jx,ix) = derv2(jx,ix) - d2k*d2r2dx2(k,1)
          derv2(jy,ix) = derv2(jy,ix) - d2k*d2r2dx2(k,6)
          derv2(jz,ix) = derv2(jz,ix) - d2k*d2r2dx2(k,5)
          derv2(jx,iy) = derv2(jx,iy) - d2k*d2r2dx2(k,6)
          derv2(jy,iy) = derv2(jy,iy) - d2k*d2r2dx2(k,2)
          derv2(jz,iy) = derv2(jz,iy) - d2k*d2r2dx2(k,4)
          derv2(jx,iz) = derv2(jx,iz) - d2k*d2r2dx2(k,5)
          derv2(jy,iz) = derv2(jy,iz) - d2k*d2r2dx2(k,4)
          derv2(jz,iz) = derv2(jz,iz) - d2k*d2r2dx2(k,3)
          derv2(jx,ix) = derv2(jx,ix) - dk
          derv2(jy,iy) = derv2(jy,iy) - dk
          derv2(jz,iz) = derv2(jz,iz) - dk
!
          dervi(jx,ix) = dervi(jx,ix) - d2ks*d2r2dx2(k,1)
          dervi(jy,ix) = dervi(jy,ix) - d2ks*d2r2dx2(k,6)
          dervi(jz,ix) = dervi(jz,ix) - d2ks*d2r2dx2(k,5)
          dervi(jx,iy) = dervi(jx,iy) - d2ks*d2r2dx2(k,6)
          dervi(jy,iy) = dervi(jy,iy) - d2ks*d2r2dx2(k,2)
          dervi(jz,iy) = dervi(jz,iy) - d2ks*d2r2dx2(k,4)
          dervi(jx,iz) = dervi(jx,iz) - d2ks*d2r2dx2(k,5)
          dervi(jy,iz) = dervi(jy,iz) - d2ks*d2r2dx2(k,4)
          dervi(jz,iz) = dervi(jz,iz) - d2ks*d2r2dx2(k,3)
          dervi(jx,ix) = dervi(jx,ix) - dks
          dervi(jy,iy) = dervi(jy,iy) - dks
          dervi(jz,iz) = dervi(jz,iz) - dks
!
          if (lgroupvelocity) then
!
!  Group velocities
!
            if (lphasecom) then
              cdk(1) = dcmplx(d2k*(xtmp(k)-xcom),d2ks*(xtmp(k)-xcom))*dcmplx(0.0_dp,1.0_dp)
              cdk(2) = dcmplx(d2k*(ytmp(k)-ycom),d2ks*(ytmp(k)-ycom))*dcmplx(0.0_dp,1.0_dp)
              cdk(3) = dcmplx(d2k*(ztmp(k)-zcom),d2ks*(ztmp(k)-zcom))*dcmplx(0.0_dp,1.0_dp)
            else
              cdk(1) = dcmplx(d2k*xtmp(k),d2ks*xtmp(k))*dcmplx(0.0_dp,1.0_dp)
              cdk(2) = dcmplx(d2k*ytmp(k),d2ks*ytmp(k))*dcmplx(0.0_dp,1.0_dp)
              cdk(3) = dcmplx(d2k*ztmp(k),d2ks*ztmp(k))*dcmplx(0.0_dp,1.0_dp)
            endif
!
            derv2dk(1:3,jx,ix) = derv2dk(1:3,jx,ix) - cdk(1:3)*dcmplx(d2r2dx2(k,1),0.0_dp)
            derv2dk(1:3,jy,ix) = derv2dk(1:3,jy,ix) - cdk(1:3)*dcmplx(d2r2dx2(k,6),0.0_dp)
            derv2dk(1:3,jz,ix) = derv2dk(1:3,jz,ix) - cdk(1:3)*dcmplx(d2r2dx2(k,5),0.0_dp)
            derv2dk(1:3,jx,iy) = derv2dk(1:3,jx,iy) - cdk(1:3)*dcmplx(d2r2dx2(k,6),0.0_dp)
            derv2dk(1:3,jy,iy) = derv2dk(1:3,jy,iy) - cdk(1:3)*dcmplx(d2r2dx2(k,2),0.0_dp)
            derv2dk(1:3,jz,iy) = derv2dk(1:3,jz,iy) - cdk(1:3)*dcmplx(d2r2dx2(k,4),0.0_dp)
            derv2dk(1:3,jx,iz) = derv2dk(1:3,jx,iz) - cdk(1:3)*dcmplx(d2r2dx2(k,5),0.0_dp)
            derv2dk(1:3,jy,iz) = derv2dk(1:3,jy,iz) - cdk(1:3)*dcmplx(d2r2dx2(k,4),0.0_dp)
            derv2dk(1:3,jz,iz) = derv2dk(1:3,jz,iz) - cdk(1:3)*dcmplx(d2r2dx2(k,3),0.0_dp)
!
            if (lphasecom) then
              cdk(1) = dcmplx(dk*(xtmp(k)-xcom),dks*(xtmp(k)-xcom))*dcmplx(0.0_dp,1.0_dp)
              cdk(2) = dcmplx(dk*(ytmp(k)-ycom),dks*(ytmp(k)-ycom))*dcmplx(0.0_dp,1.0_dp)
              cdk(3) = dcmplx(dk*(ztmp(k)-zcom),dks*(ztmp(k)-zcom))*dcmplx(0.0_dp,1.0_dp)
            else
              cdk(1) = dcmplx(dk*xtmp(k),dks*xtmp(k))*dcmplx(0.0_dp,1.0_dp)
              cdk(2) = dcmplx(dk*ytmp(k),dks*ytmp(k))*dcmplx(0.0_dp,1.0_dp)
              cdk(3) = dcmplx(dk*ztmp(k),dks*ztmp(k))*dcmplx(0.0_dp,1.0_dp)
            endif
!
            derv2dk(1:3,jx,ix) = derv2dk(1:3,jx,ix) - cdk(1:3)
            derv2dk(1:3,jy,iy) = derv2dk(1:3,jy,iy) - cdk(1:3)
            derv2dk(1:3,jz,iz) = derv2dk(1:3,jz,iz) - cdk(1:3)
          endif
!
!  rpd arrays no longer needed - use to save cos/sin for use in EEM/QEq
!
          rpd(k,1) = cosk
          rpd(k,2) = sink
        enddo
        if (radsum.gt.0.0_dp) then
!
!  Radial components
!
          do k = 1,nor
            if (lphasecom) then
              cosk = xkv*(xtmp(k) - xcom) + ykv*(ytmp(k) - ycom) + zkv*(ztmp(k) - zcom)
            else
              cosk = xkv*xtmp(k) + ykv*ytmp(k) + zkv*ztmp(k)
            endif
            sink = sin(cosk)
            cosk = cos(cosk) - oneij
            rrtrm = rderiv(k)*cosk
            ritrm = rderiv(k)*sink
            if (radi.gt.0.0_dp) then
              iis = indri
              ifs = indrif
              derv2(ixf,iis) = derv2(ixf,iis) - rrtrm*xtmp(k)
              derv2(iyf,iis) = derv2(iyf,iis) - rrtrm*ytmp(k)
              derv2(izf,iis) = derv2(izf,iis) - rrtrm*ztmp(k)
              derv2(jx,iis)  = derv2(jx,iis)  + rrtrm*xtmp(k)
              derv2(jy,iis)  = derv2(jy,iis)  + rrtrm*ytmp(k)
              derv2(jz,iis)  = derv2(jz,iis)  + rrtrm*ztmp(k)
              dervi(ixf,iis) = dervi(ixf,iis) - ritrm*xtmp(k)
              dervi(iyf,iis) = dervi(iyf,iis) - ritrm*ytmp(k)
              dervi(izf,iis) = dervi(izf,iis) - ritrm*ztmp(k)
              dervi(jx,iis)  = dervi(jx,iis)  + ritrm*xtmp(k)
              dervi(jy,iis)  = dervi(jy,iis)  + ritrm*ytmp(k)
              dervi(jz,iis)  = dervi(jz,iis)  + ritrm*ztmp(k)
              derv2(ifs,iis) = derv2(ifs,iis) + rtrm2(k)*cosk
              dervi(ifs,iis) = dervi(ifs,iis) + rtrm2(k)*sink
!
              if (lgroupvelocity) then
!
!  Group velocities
!
                if (lphasecom) then
                  cdk(1) = dcmplx(rrtrm*(xtmp(k)-xcom),ritrm*(xtmp(k)-xcom))*dcmplx(0.0_dp,1.0_dp)
                  cdk(2) = dcmplx(rrtrm*(ytmp(k)-ycom),ritrm*(ytmp(k)-ycom))*dcmplx(0.0_dp,1.0_dp)
                  cdk(3) = dcmplx(rrtrm*(ztmp(k)-zcom),ritrm*(ztmp(k)-zcom))*dcmplx(0.0_dp,1.0_dp)
                else
                  cdk(1) = dcmplx(rrtrm*xtmp(k),ritrm*xtmp(k))*dcmplx(0.0_dp,1.0_dp)
                  cdk(2) = dcmplx(rrtrm*ytmp(k),ritrm*ytmp(k))*dcmplx(0.0_dp,1.0_dp)
                  cdk(3) = dcmplx(rrtrm*ztmp(k),ritrm*ztmp(k))*dcmplx(0.0_dp,1.0_dp)
                endif
!
                derv2dk(1:3,jx,iis) = derv2dk(1:3,jx,iis) + cdk(1:3)*dcmplx(xtmp(k),0.0_dp)
                derv2dk(1:3,jy,iis) = derv2dk(1:3,jy,iis) + cdk(1:3)*dcmplx(ytmp(k),0.0_dp)
                derv2dk(1:3,jz,iis) = derv2dk(1:3,jz,iis) + cdk(1:3)*dcmplx(ztmp(k),0.0_dp)
              endif
            endif
          enddo
        endif
!
!  If nor = 0, then rejoin here since there can still be a contribution to 
!  the second derivatives from the variable charges and the self term.
!
1110    continue
!********************************************
!  Variable charge contribution to phonons  *
!********************************************
        if (lDoQDeriv2) then
!
!  Calculate phased terms
!
          hfactor = 1.0_dp
          if (i.eq.j) hfactor = 0.5_dp
          d2self = hfactor*d2self
          derive0selfi = hfactor*derive0self*qlj
          derive0selfj = hfactor*derive0self*qli
          d1ix = 0.0_dp
          d1iy = 0.0_dp
          d1iz = 0.0_dp
          d1jx = 0.0_dp
          d1jy = 0.0_dp
          d1jz = 0.0_dp
          do k = 1,nor
            d0i(k) = d0i(k)*hfactor
            d0j(k) = d0j(k)*hfactor
            d1i(k) = d1i(k)*hfactor
            d1j(k) = d1j(k)*hfactor
            d2i2(k) = d2i2(k)*hfactor
            d2ij(k) = d2ij(k)*hfactor
            d2j2(k) = d2j2(k)*hfactor
            d1ix = d1ix + d1i(k)*xtmp(k)
            d1iy = d1iy + d1i(k)*ytmp(k)
            d1iz = d1iz + d1i(k)*ztmp(k)
            d1jx = d1jx + d1j(k)*xtmp(k)
            d1jy = d1jy + d1j(k)*ytmp(k)
            d1jz = d1jz + d1j(k)*ztmp(k)
          enddo
!
!  Apply variable charge correction to second derivatives
!
          call d2chargep(i,j,nor,ix,iy,iz,jx,jy,jz,xkv,ykv,zkv,d0i,d0j,d1ix,d1iy,d1iz, &
                         d1jx,d1jy,d1jz,d2i2,d2ij,d2j2,d2self,derive0selfi,derive0selfj,.true.)
        endif
1120    continue
      enddo
    enddo
!
!  End of real space part - perform general tasks
!
  endif
!
!  Free local memory
!
  deallocate(npotl,stat=status)
  if (status/=0) call deallocate_error('realp','npotl')
!
!  Timing
!
  time2 = g_cpu_time()
  tatom = tatom + time2 - time1
#ifdef TRACE
  call trace_out('realpd')
#endif
!
  return
  end
