  subroutine real0dd(eatom,ereal,eqeq,lgrad1,lgrad2)
!
!  Subroutine for calculating the real space energy for a finite cluster => 0 D system
!  Distributed second derivatives version. Freezing not allowed.
!
!   5/13 Created from real0d
!   9/16 cputime renamed to g_cpu_time
!   9/16 constants renamed to g_constants
!   9/16 lorder12 renamed to lorder12loc
!   1/17 Globalisation of scrho added
!   7/17 Calls to d1charge changed to d1charged
!   7/17 Site energies added
!   2/18 Trace added
!   5/18 Modified for q ranges
!   8/19 Short range damping of polarisation added
!   8/19 Trap for neemrptr being zero added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!   3/20 Cutoffs introduced for case where two non-charged particles interact
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
  use derivatives
  use eam,            only : lMEAMden
  use eemdata
  use element
  use energies,       only : eregion2region, siteenergy
  use general,        only : cutw
  use iochannels,     only : ioout
  use molecule
  use optimisation
  use parallel
  use polarise
  use potentialxyz
  use realvectors,    only : derivqd
  use shells
  use splinedata
  use sutton
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
  logical,     intent(in)                      :: lgrad2
  real(dp),    intent(inout)                   :: eatom
  real(dp),    intent(inout)                   :: ereal
  real(dp),    intent(inout)                   :: eqeq
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ii
  integer(i4)                                  :: iar
  integer(i4)                                  :: ind
  integer(i4)                                  :: indmi
  integer(i4)                                  :: indmj
  integer(i4)                                  :: indri
  integer(i4)                                  :: indrif
  integer(i4)                                  :: indrj
  integer(i4)                                  :: ix
  integer(i4)                                  :: iy
  integer(i4)                                  :: iz
  integer(i4)                                  :: ixx
  integer(i4)                                  :: iyy
  integer(i4)                                  :: izz
  integer(i4)                                  :: j
  integer(i4)                                  :: jx
  integer(i4)                                  :: jy
  integer(i4)                                  :: jz
  integer(i4)                                  :: jxx
  integer(i4)                                  :: jyy
  integer(i4)                                  :: jzz
  integer(i4)                                  :: k
  integer(i4)                                  :: m
  integer(i4)                                  :: n
  integer(i4)                                  :: nat1
  integer(i4)                                  :: nat2
  integer(i4)                                  :: nati
  integer(i4)                                  :: natj
  integer(i4)                                  :: nbtypeij
  integer(i4)                                  :: nbtypeij2
  integer(i4)                                  :: nff
  integer(i4)                                  :: nffl
  integer(i4)                                  :: nfi
  integer(i4)                                  :: nfj
  integer(i4)                                  :: nmi
  integer(i4)                                  :: nmj
  integer(i4)                                  :: nor
  integer(i4)                                  :: npot
  integer(i4), dimension(:), allocatable       :: npotl
  integer(i4)                                  :: npots
  integer(i4)                                  :: npt
  integer(i4)                                  :: nqri
  integer(i4)                                  :: nqrj
  integer(i4)                                  :: nregioni
  integer(i4)                                  :: nregionj
  integer(i4)                                  :: nregiontypi
  integer(i4)                                  :: nregiontypj
  integer(i4)                                  :: ntyp1
  integer(i4)                                  :: ntyp2
  integer(i4)                                  :: ntypi
  integer(i4)                                  :: ntypj
  integer(i4)                                  :: status
  logical                                      :: l2bonds
  logical                                      :: l3bonds
  logical                                      :: lbonded
  logical                                      :: lbreathe
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
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: cut2
  real(dp)                                     :: cut2p
  real(dp)                                     :: cut2q
  real(dp)                                     :: cut2r
  real(dp)                                     :: cut2s
  real(dp)                                     :: d0i
  real(dp)                                     :: d0j
  real(dp)                                     :: d1i
  real(dp)                                     :: d1j
  real(dp)                                     :: d1ix
  real(dp)                                     :: d1iy
  real(dp)                                     :: d1iz
  real(dp)                                     :: d1jx
  real(dp)                                     :: d1jy
  real(dp)                                     :: d1jz
  real(dp)                                     :: d2i2
  real(dp)                                     :: d2ij
  real(dp)                                     :: d2j2
  real(dp)                                     :: d2self
  real(dp)                                     :: deriv
  real(dp)                                     :: deriv2
  real(dp)                                     :: deriv3
  real(dp)                                     :: derive0
  real(dp)                                     :: derive
  real(dp)                                     :: derive2
  real(dp)                                     :: derive3
  real(dp)                                     :: dist
  real(dp)                                     :: eatoml
  real(dp)                                     :: ec6
  real(dp)                                     :: ec6l
  real(dp)                                     :: eqeql
  real(dp)                                     :: ereall
  real(dp)                                     :: esum
  real(dp)                                     :: etrm
  real(dp)                                     :: etrm1
  real(dp)                                     :: etrm2
  real(dp)                                     :: factor
  real(dp)                                     :: fct
  real(dp)                                     :: oci
  real(dp)                                     :: ocj
  real(dp)                                     :: ofct
  real(dp)                                     :: ofctij
  real(dp)                                     :: ospfct
  real(dp)                                     :: qli
  real(dp)                                     :: qlj
  real(dp)                                     :: r
  real(dp)                                     :: radi
  real(dp)                                     :: radj
  real(dp)                                     :: radsum
  real(dp)                                     :: rderiv
  real(dp)                                     :: rdiff
  real(dp)                                     :: rp
  real(dp)                                     :: rpdl(6)
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
  real(dp)                                     :: xal
  real(dp)                                     :: yal
  real(dp)                                     :: zal
  real(dp)                                     :: xcrd
  real(dp)                                     :: ycrd
  real(dp)                                     :: zcrd
#ifdef TRACE
  call trace_in('real0dd')
#endif
!
  time1 = g_cpu_time()
!
!  Local variables
!
  small = 1.0d-12
  small2 = 1.0d-2
  tsuml = 0.0_dp
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
!  Zero energy
!
  ec6 = 0.0_dp
!
!  Allocate local memory
!
  allocate(npotl(npote),stat=status)
  if (status/=0) call outofmemory('real0dd','npotl')
!
!  Set flag as to whether first derivatives must be calculated in
!  twobody, either for energy derivatives or the site potential
!
  lgrad1p = (lgrad1.or.lpolar)
  if (lnoreal) goto 999
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
!
!  Find total number of variable atoms
!
  nff = numat
  nffl = natomsonnode
!
!  Outer loop over atoms on node
!
  ix = - 2
  iy = - 1
  iz =   0
  nfi =   0
  do ii = 1,natomsonnode
    i = node2atom(ii)
    xal = xclat(i)
    yal = yclat(i)
    zal = zclat(i)
    nati = nat(i)
    ntypi = nftype(i)
    nregioni = nregionno(nsft+i)
    nregiontypi = nregiontype(nregioni,ncf)
    qli = qf(i)
    oci = occuf(i)
    lopi = .true.
    ix = ix + 3
    iy = iy + 3
    iz = iz + 3
    nfi = nfi + 1
    if (lbsmat(i+nsft)) then
      radi = radf(i)
      indri = 3*nffl + nfi
      indrif = 3*nff + i
    else
      radi = 0.0_dp
    endif
    if (leem.and.lmultiqrange.and.neemrptr(i).ne.0) then
      nqri = nqrnow(neemrptr(i))
    else
      nqri = 1
    endif
!
!  Molecule handling
!
    if (lmol) then
      nmi = natmol(i)
      indmi = nmolind(i)
    endif
!
!  Inner loop over second site
!
    jx = - 2
    jy = - 1
    jz =   0
    nfj =   0
    jloop: do j = 1,numat
      lopj = .true.
      jx = jx + 3
      jy = jy + 3
      jz = jz + 3
      nfj = nfj + 1
      nregionj = nregionno(nsft+j)
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
      natj = nat(j)
      ntypj = nftype(j)
      qlj = qf(j)
      ocj = occuf(j)
      if (lbsmat(nsft+j)) then
        radj = radf(j)
        indrj = 3*nff + nfj
      else
        radj = 0.0_dp
      endif
      xcrd = xclat(j) - xal
      ycrd = yclat(j) - yal
      zcrd = zclat(j) - zal
      radsum = radi + radj
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
      if (i.eq.j) then
        ofct = 0.5_dp*oci*ocj
        ofctij = 1.0_dp
      else
        ofct = oci*ocj
        ofctij = 0.5_dp
      endif
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
      endif
!
!  Locate potential number
!  Check whether potential requires specific types
!
      rp = 0.0_dp
      npots = 0
      lneedmol = (lmol.and..not.lmolq)
      do n = 1,npote
        if (lmatch(nat1,ntyp1,nspec1(n),nptyp1(n),.true.)) then
          if (lmatch(nat2,ntyp2,nspec2(n),nptyp2(n),.true.)) then
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
        cut2 = max(cut2,cut2q)
      endif
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
      if (abs(r-small2).lt.1.0d-12) r = small2
      if (r.lt.small) then
!
!  Core-shell spring constant at zero distance - correct second derivative matrix
!
        if (lgrad2) then
          do k = 1,npots
            npot = npotl(k)
            npt = nptype(npot)
            if (npt.eq.5.or.npt.eq.8.or.npt.eq.33) then
              apt = twopot(1,npot)*ospfct
              derv2(jx,ix) = derv2(jx,ix) - apt
              derv2(jy,iy) = derv2(jy,iy) - apt
              derv2(jz,iz) = derv2(jz,iz) - apt
            endif
          enddo
        endif
        cycle jloop
      elseif (r.gt.cut2) then
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
      eatoml = 0.0_dp
      ereall = 0.0_dp
      ec6l = 0.0_dp
      call twobody1(eatoml,ereall,ec6l,lgrad1p,lgrad2,.false.,nor,1_i4,dist,xcrd,ycrd,zcrd,d0i,d0j, &
                    deriv,deriv2,deriv3,derive0,derive,derive2,derive3,rderiv,npots,npotl,cut2r, &
                    cut2q,cut2s,lptrmol,0_i4,factor,ofct,ospfct,radsum,rtrm1,rtrm2,rtrm3,rtrm32,sctrm1, &
                    sctrm2,qli,qlj,lcspair,.false.,.false.,lbonded,l2bonds,l3bonds,nbtypeij, &
                    nbtypeij2,.false.,lQMMMelectro,d1i,d1j,d2i2,d2ij,d2j2)
!
      eatoml  = ofctij*eatoml
      ereall  = ofctij*ereall
      ec6l    = ofctij*ec6l
!
      eatom = eatom + eatoml
      ereal = ereal + ereall
      ec6 = ec6 + ec6l
!
      esum = eatoml + ec6l + ereall
      eregion2region(nregionj,nregioni) = eregion2region(nregionj,nregioni) + esum
!
      siteenergy(i) = siteenergy(i) + esum
!
      if (lPrintTwo) then
        write(ioout,'(4x,i12,i12,f22.10,1x,f22.10)') i,j,eatoml+ec6l,ereall
      endif
!
      if ((lDoQDeriv1.or.lDoQDeriv2).and.lgrad1) then
        d0i = d0i + derive0*qlj
        d0j = d0j + derive0*qli
      endif
      if (leem) then
        eqeql = 0.0_dp
        if (lqeq) then
          call qeqbody1(eqeql,lgrad1,lgrad2,nor,1_i4,dist,deriv,deriv2, &
                        fct,qli,qlj,nati,natj,nqri,nqrj,d1i,d1j,d2i2,d2ij,d2j2)
        elseif (lSandM) then
          call smbody1(eqeql,lgrad1,lgrad2,nor,1_i4,dist,deriv,deriv2, &
                        fct,qli,qlj,nati,natj,nqri,nqrj,d1i,d1j,d2i2,d2ij,d2j2)
        endif
!
        eqeql  = ofctij*eqeql
        eqeq = eqeq + eqeql
!
        eregion2region(nregionj,nregioni) = eregion2region(nregionj,nregioni) + eqeql
!
        siteenergy(i) = siteenergy(i) + eqeql
      endif
      if (lsuttonc) then
        if (.not.lMEAMden) then
          if (lorder12loc) then
            scrho(1,i) = scrho(1,i) + ocj*sctrm1
          else
            scrho(1,i) = scrho(1,i) + ocj*sctrm2
          endif
        endif
      endif
!***********************
!  Radial derivatives  *
!***********************
      if (lgrad1) then
        if (radi.gt.0.0_dp) then
          raderv(i) = raderv(i) + rtrm1
          if (lgrad2) then
            derv2(indrif,indri) = derv2(indrif,indri) + rtrm2
          endif
        endif
        if (radj.gt.0.0_dp) then
          if (lgrad2) then
            if (radi.gt.0.0_dp) then
              derv2(indrj,indri) = derv2(indrj,indri) + rtrm2
            endif
          endif
        endif
      endif
!*****************************
!  Charge first derivatives  *
!*****************************
      if (lgrad1.and.lDoQDeriv1) then
        call d1charged(ii,i,lopi,1_i4,d0i)
      endif
!*************************************
!  Electrostatic potential on-sites  *
!*************************************
      if (lpolar) then
        vx(i) = vx(i) - qlj*derivqd(1)*xcrd
        vy(i) = vy(i) - qlj*derivqd(1)*ycrd
        vz(i) = vz(i) - qlj*derivqd(1)*zcrd
      endif
!**************************
!  Coordinate Derivatives *
!**************************
!
!  First derivatives
!
      if (lgrad1) then
        xdrv(i) = xdrv(i) - deriv*xcrd
        ydrv(i) = ydrv(i) - deriv*ycrd
        zdrv(i) = zdrv(i) - deriv*zcrd
        if (lgrad2.and.radi.gt.0.0_dp) then
          derv2(indrif,ix) = derv2(indrif,ix) - rderiv*xcrd
          derv2(indrif,iy) = derv2(indrif,iy) - rderiv*ycrd
          derv2(indrif,iz) = derv2(indrif,iz) - rderiv*zcrd
        endif
        if (nregioni.ne.nregionj) then
          xregdrv(nregioni) = xregdrv(nregioni) - deriv*xcrd
          yregdrv(nregioni) = yregdrv(nregioni) - deriv*ycrd
          zregdrv(nregioni) = zregdrv(nregioni) - deriv*zcrd
        endif
      endif
!
!  Second derivatives
!
      if (lgrad2) then
        rpdl(1) = xcrd*xcrd
        rpdl(2) = ycrd*ycrd
        rpdl(3) = zcrd*zcrd
        rpdl(4) = ycrd*zcrd
        rpdl(5) = xcrd*zcrd
        rpdl(6) = xcrd*ycrd
!
        derv2(jx,ix) = derv2(jx,ix) - deriv2*rpdl(1)
        derv2(jy,ix) = derv2(jy,ix) - deriv2*rpdl(6)
        derv2(jz,ix) = derv2(jz,ix) - deriv2*rpdl(5)
        derv2(jx,iy) = derv2(jx,iy) - deriv2*rpdl(6)
        derv2(jy,iy) = derv2(jy,iy) - deriv2*rpdl(2)
        derv2(jz,iy) = derv2(jz,iy) - deriv2*rpdl(4)
        derv2(jx,iz) = derv2(jx,iz) - deriv2*rpdl(5)
        derv2(jy,iz) = derv2(jy,iz) - deriv2*rpdl(4)
        derv2(jz,iz) = derv2(jz,iz) - deriv2*rpdl(3)
        derv2(jx,ix) = derv2(jx,ix) - deriv
        derv2(jy,iy) = derv2(jy,iy) - deriv
        derv2(jz,iz) = derv2(jz,iz) - deriv
!
!  Coordinate - radius mixed
!
        if (radi.gt.0.0_dp) then
          derv2(jx,indri) = derv2(jx,indri) + rderiv*xcrd
          derv2(jy,indri) = derv2(jy,indri) + rderiv*ycrd
          derv2(jz,indri) = derv2(jz,indri) + rderiv*zcrd
        endif
!*******************************************************
!  Variable charge contribution to second derivatives  *
!*******************************************************
        if (lDoQDeriv2) then
          d2self = 0.0_dp
          d1ix = d1i*xcrd
          d1iy = d1i*ycrd
          d1iz = d1i*zcrd
          d1jx = d1j*xcrd
          d1jy = d1j*ycrd
          d1jz = d1j*zcrd
          call d2charge(i,j,1_i4,ix,iy,iz,jx,jy,jz,lopi,lopj,d0i,d0j,d1ix,d1iy,d1iz,d1jx,d1jy,d1jz, &
                        d1i,d1j,d2i2,d2ij,d2j2,d2self,0.0_dp,0.0_dp,.true.,.true.)
        endif
      endif
!
!  Skip to here if frozen pair
!
    enddo jloop
  enddo
!
!  Breathing shell self terms
!
  nfi = 0
  iloop: do ii = 1,natomsonnode
    i = node2atom(ii)
    lopi = .true.
    nregioni = nregionno(nsft+i)
    nregiontypi = nregiontype(nregioni,ncf)
!
!  QM/MM handling : i is a QM atom => exclude
!
    if (QMMMmode(ncf).gt.0) then
      if (nregiontypi.eq.1) cycle iloop
    endif
!
    nfi = nfi + 1
    iar = nsft + nrelf2a(i)
    lbreathe = lbsmat(iar)
    if (lbreathe) then
      nati = nat(i)
      ntypi = nftype(i)
      oci = occuf(i)
      radi = radf(i)
      indri = 3*nff + nfi
      indrif = 3*nffl + i
      if (nati.gt.maxele) nati = nati - maxele
!******************************
!  Breathing shell self term  *
!******************************
      eatoml = 0.0_dp
      do m = 1,npote
        if (nptype(m).eq.14) then
          if (lmatch(nati,ntypi,nspec1(m),nptyp1(m),.true.)) then
            apt = twopot(1,m)*oci
            rdiff = radi - twopot(2,m)
            eatoml = eatoml + 0.5_dp*apt*rdiff*rdiff
            if (lgrad1) then
              raderv(i) = raderv(i) + apt*rdiff
              if (lgrad2) then
                derv2(indrif,indri) = derv2(indrif,indri) + apt
              endif
            endif
          endif
        elseif (nptype(m).eq.17) then
          if (lmatch(nati,ntypi,nspec1(m),nptyp1(m),.true.)) then
            apt = twopot(1,m)*oci
            bpt = twopot(2,m)
            rdiff = radi - twopot(3,m)
            etrm1 = exp(bpt*rdiff)
            etrm2 = 1.0_dp/etrm1
            etrm = apt*(etrm1 + etrm2)
            eatoml = eatoml + etrm
            if (lgrad1) then
              raderv(i) = raderv(i) + apt*bpt*(etrm1 - etrm2)
              if (lgrad2) then
                derv2(indrif,indri) = derv2(indrif,indri) + bpt*bpt*etrm
              endif
            endif
          endif
        elseif (nptype(m).eq.31) then
          if (lmatch(nati,ntypi,nspec1(m),nptyp1(m),.true.)) then
            apt = twopot(1,m)*oci
            bpt = twopot(2,m)
            rdiff = radi - twopot(3,m)
            etrm1 = exp(bpt*rdiff)
            etrm = apt*etrm1
            eatoml = eatoml + etrm
            if (lgrad1) then
              raderv(i) = raderv(i) + apt*bpt*etrm1
              if (lgrad2) then
                derv2(indrif,indri) = derv2(indrif,indri) + bpt*bpt*etrm
              endif
            endif
          endif
        endif
!
        eregion2region(nregioni,nregioni) = eregion2region(nregioni,nregioni) + eatoml
!
        siteenergy(i) = siteenergy(i) + eatoml
!
      enddo
!
!  Add energy to total
!     
      eatom = eatom + eatoml*ofctij
    endif
  enddo iloop
  if (nprocs.gt.1) then
!****************
!  Global sums  *
!****************
    if (lsuttonc) then
      tsum0 = g_cpu_time()
!
      allocate(sum(numat),stat=status)
      if (status/=0) call outofmemory('real0dd','sum')
      allocate(sum2(numat),stat=status)
      if (status/=0) call outofmemory('real0dd','sum2')
!
      if (.not.lMEAMden) then
        do i = 1,numat
          sum2(i) = scrho(1,i)
        enddo
        call sumall(sum2,sum,numat,"real0dd","scrho")
        do i = 1,numat
          scrho(1,i) = sum(i)
        enddo
        do i = 1,numat
          sum2(i) = scrho12(1,i)
        enddo
        call sumall(sum2,sum,numat,"real0dd","scrho12")
        do i = 1,numat
          scrho12(1,i) = sum(i)
        enddo
      endif
!
      deallocate(sum2,stat=status)
      if (status/=0) call deallocate_error('real0dd','sum2')
      deallocate(sum,stat=status)
      if (status/=0) call deallocate_error('real0dd','sum')
!
      tsuml = g_cpu_time() - tsum0
      tsum = tsum + tsuml
    endif
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
!  End of real space part - perform general tasks
!
999 continue
!
!  Add any C6 terms to ereal
!
  ereal = ereal + ec6
!
!  Free local memory
!
  deallocate(npotl,stat=status)
  if (status/=0) call deallocate_error('real0dd','npotl')
!
!  Timing
!
  time2 = g_cpu_time()
  tatom = tatom + time2 - time1 - tsuml
#ifdef TRACE   
  call trace_out('real0dd')
#endif
!
  return
  end
