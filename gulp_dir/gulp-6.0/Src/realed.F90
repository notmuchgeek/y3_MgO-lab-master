  subroutine realed(eatom,ereal,erecip,ec6,eqeq,esregion12,esregion2,eattach,lgrad1,lgrad2)
!
!  Subroutine for calculating real space energy and up to second derivatives.
!  Distributed second derivatives version. Freezing not allowed. 
!
!   5/13 Created from reale.
!  12/14 rtrm1 changed from scalar to array
!   9/16 cputime renamed to g_cpu_time
!   9/16 constants renamed to g_constants
!   9/16 lorder12 renamed to lorder12loc
!   9/16 ofctij included for rstrd/sderv2 due to no longer being
!        upper half triangular pairwise sum
!  12/16 eatom from breathing shell self term added to total energy
!   1/17 Globalisation of scrho added
!   7/17 Calls to d1charge changed to d1charged
!   7/17 Debugged in parallel
!   7/17 Site energies added
!   8/17 Reference to eatoml in self energy corrected to ereall
!   2/18 Trace added
!   5/18 Modified for q ranges
!   8/18 Call to twostrterms introduced for setting of rpd
!   9/18 Strain module introduced
!   9/18 Second derivatives completed here rather than in strfin
!   5/19 Skipping of ij loop now checks for polarisation
!   8/19 Short range damping of polarisation added
!   8/19 Trap for neemrptr being zero added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!  11/19 Rigid molecule modifications added
!   3/20 Electrostatic cutoff only included where the either charge exceeds the threshold
!   4/20 derv3c and d2r2dsdc added for benefit of rigid molecules
!   4/20 derv3c changes reversed as they are no longer required
!   7/20 Modifications for gfortran v10
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
  use configurations, only : lbsmat, lsliceatom, nregions, nregionno, nregiontype, QMMMmode
  use g_constants
  use control
  use current
  use datatypes
  use derivatives
  use eam,            only : lMEAMden
  use eemdata
  use element
  use energies,       only : eregion2region, siteenergy
  use general,        only : cutw
  use iochannels,     only : ioout
  use kspace
  use m_strain,       only : twostrterms
  use mdlogic
  use molecule
  use optimisation
  use parallel
  use polarise
  use potentialxyz
  use realvectors
  use shells
  use sutton
  use symmetry
  use thresholds,    only : thresh_q
  use times
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  use two
  implicit none
!
!  Passed variables
!
  logical,     intent(in)                      :: lgrad1
  logical,     intent(in)                      :: lgrad2
  real(dp),    intent(inout)                   :: eatom
  real(dp),    intent(inout)                   :: eattach
  real(dp),    intent(inout)                   :: ec6
  real(dp),    intent(inout)                   :: eqeq
  real(dp),    intent(inout)                   :: ereal
  real(dp),    intent(inout)                   :: erecip
  real(dp),    intent(inout)                   :: esregion12
  real(dp),    intent(inout)                   :: esregion2
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ii
  integer(i4)                                  :: indm
  integer(i4)                                  :: indmj
  integer(i4)                                  :: indri
  integer(i4)                                  :: indrif
  integer(i4)                                  :: indrj
  integer(i4)                                  :: ix
  integer(i4)                                  :: iy
  integer(i4)                                  :: iz
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
  integer(i4)                                  :: kk
  integer(i4)                                  :: kl
  integer(i4)                                  :: ks
  integer(i4)                                  :: kt
  integer(i4)                                  :: m
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
  integer(i4)                                  :: nff
  integer(i4)                                  :: nffl
  integer(i4)                                  :: nfi
  integer(i4)                                  :: nfj
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
  logical                                      :: lattach
  logical                                      :: lc6loc
  logical                                      :: lcspair
  logical                                      :: ldog1
  logical                                      :: ldog2
  logical                                      :: lewaldtype
  logical                                      :: lgrad1p
  logical                                      :: lmatch
  logical                                      :: lmdl
  logical                                      :: lmolok
  logical                                      :: lneedmol
  logical                                      :: lopi
  logical                                      :: lopj
  logical                                      :: lorder12loc
  logical                                      :: lQMMMelectro
  logical                                      :: lreg2one
  logical                                      :: lreg2pair
  logical                                      :: lself
  logical                                      :: lsg1
  logical                                      :: lsg2
  logical                                      :: lslicei
  logical                                      :: lslicej
  real(dp)                                     :: apt
  real(dp)                                     :: bpt
  real(dp)                                     :: c6tot
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: cut2
  real(dp)                                     :: cut2e
  real(dp)                                     :: cut2p
  real(dp)                                     :: cut2q
  real(dp)                                     :: cut2r
  real(dp)                                     :: cut2s
  real(dp)                                     :: cut2w
  real(dp)                                     :: d1is(6)
  real(dp)                                     :: d1ix
  real(dp)                                     :: d1iy
  real(dp)                                     :: d1iz
  real(dp)                                     :: d1js(6)
  real(dp)                                     :: d1jx
  real(dp)                                     :: d1jy
  real(dp)                                     :: d1jz
  real(dp)                                     :: d2self
  real(dp)                                     :: derive0self
  real(dp)                                     :: derive0selfi
  real(dp)                                     :: derive0selfj
  real(dp)                                     :: dqj(1)
  real(dp)                                     :: eatoml
  real(dp)                                     :: ec6l
  real(dp)                                     :: eqeql
  real(dp)                                     :: ereall
  real(dp)                                     :: erecipl
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
  real(dp)                                     :: polfct
  real(dp)                                     :: qli
  real(dp)                                     :: qlj
  real(dp)                                     :: radi
  real(dp)                                     :: radj
  real(dp)                                     :: radsum
  real(dp)                                     :: rdiff
  real(dp)                                     :: rp
  real(dp)                                     :: rqeq2
  real(dp)                                     :: rstrdloc(6)
  real(dp)                                     :: rt2
  real(dp)                                     :: sctrm1
  real(dp)                                     :: sctrm2
  real(dp),    dimension(:), allocatable       :: sum
  real(dp),    dimension(:), allocatable       :: sum2
  real(dp)                                     :: time1
  real(dp)                                     :: time2
  real(dp)                                     :: tsum0
  real(dp)                                     :: tsuml
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
  call trace_in('realed')
#endif
!
  time1 = g_cpu_time()
!
!  Set local variables
!
  lgrad1p = (lgrad1.or.lpolar)
  lc6loc = (lc6.and.ndim.eq.3)
  lmdl = lmd
  if (.not.lgrad1) lmdl = .false.
  rqeq2 = rqeq*rqeq
  tsuml = 0.0_dp
!
!  Set the Coulomb term type based on dimensionality :
!
!  1-D         => 1/r
!  2-D and 3-D => erfc(seta*r)/r
!
  lewaldtype = (ndim.ne.1.or.lwolf)
!
  lsg1 = (lgrad1.and.lstr)
  lsg2 = (lgrad2.and.lstr)
!
!  Find number of unfrozen atoms
!
  nff = numat
  nffl = natomsonnode
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
  allocate(npotl(npote),stat=status)
  if (status/=0) call outofmemory('realed','npotl')
!
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
!***************************************************************
!  Atomistic and real space electrostatic component of energy  *
!***************************************************************
  ix = - 2
  iy = - 1
  iz = 0
  nfi = 0
!
!  Loop over atoms on local node
!
  do ii = 1,natomsonnode
    i = node2atom(ii)
    xal = xclat(i)
    yal = yclat(i)
    zal = zclat(i)
    nati = nat(i)
    ntypi = nftype(i)
    qli = qf(i)
    oci = occuf(i)
    nregioni = nregionno(nsft+nrelf2a(i))
    nregiontypi = nregiontype(nregioni,ncf)
    lopi = .true.
    lslicei = lsliceatom(nsft + nrelf2a(i))
    nfi = nfi + 1
    ix = ix + 3
    iy = iy + 3
    iz = iz + 3
    if (lbsmat(nrelf2a(i)+nsft)) then
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
    jx = - 2
    jy = - 1
    jz = 0
    nfj = 0
    jloop: do j = 1,numat
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
      lopj = .true.
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
!  QM/MM : Set electrostatic embedding flag : If either i or j are QM atoms => exclude electrostatics
!       
      lQMMMelectro = (QMMMmode(ncf).eq.2.and.(nregiontypi.eq.1.or.nregiontypj.eq.1))
!
      nfj = nfj + 1
      jx = jx + 3
      jy = jy + 3
      jz = jz + 3
!
      xcrd = xclat(j) - xal
      ycrd = yclat(j) - yal
      zcrd = zclat(j) - zal
      qlj = qf(j)
      ocj = occuf(j)
      if (lbsmat(nsft+nrelf2a(j))) then
        radj = radf(j)
        indrj = 3*nff + nfj
      else
        radj = 0.0_dp
      endif
      radsum = radi + radj
      if (i.eq.j) then
        ofct = 0.5_dp*oci*ocj
        ofctij = 1.0_dp
      else
        ofct = oci*ocj
        ofctij = 0.5_dp
      endif
      fct = ofct*angstoev
      factor = qli*qlj*fct
      if (lpolar) then
        polfct = abs(qli*dpolar(nrelf2a(j))) + abs(qlj*dpolar(nrelf2a(i)))
      else
        polfct = 0.0_dp
      endif
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
!  If no valid potentials, charge product is zero, and no polarisation then skip loop
!
      if (npots.eq.0.and.abs(factor).lt.1.0d-8.and.polfct.lt.1.0d-8) then
        nor = 0
        goto 1110
      endif
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
      ereall = 0.0_dp
      erecipl = 0.0_dp
      ec6l = 0.0_dp
      derive0self = 0.0_dp
      if (lself) call selfterm(ereall,erecipl,ec6l,derive0self,factor,fct,ofct,ospfct,1.0_dp,npotl,npots, &
                               c6tot,d2self,lgrad1,lgrad2,i,j,ix,jx,1.0_dp,0.5_dp,lewaldtype,qli,qlj)
!
      ereall  = ofctij*ereall
      erecipl = ofctij*erecipl
      ec6l    = ofctij*ec6l
!
      esum = ereall + erecipl + ec6l
      if (lreg2one) then
        esregion12 = esregion12 + esum
      elseif (lreg2pair) then
        esregion2 = esregion2 + esum
      else
        ereal = ereal + ereall
        erecip = erecip + erecipl
        ec6 = ec6 + ec6l
      endif
!
      eregion2region(nregionj,nregioni) = eregion2region(nregionj,nregioni) + esum
!
      siteenergy(i) = siteenergy(i) + esum
!
      if (lattach) eattach = eattach + esum
!*******************************************
!  Charge first derivatives for self-term  *
!*******************************************
      if (lself.and.lgrad1.and.lDoQDeriv1) then
        dqj(1) = 2.0_dp*derive0self*qlj
        call d1charged(ii,i,lopi,1_i4,dqj)
      endif
!
      ldog1 = (lmdl.or.lsg1.or.lbsmat(nsft+nrelf2a(i)).or.leem.or.lDoQDeriv1)
      ldog2 = (lsg2.or.lbsmat(nsft+nrelf2a(i)).or.leem.or.lDoQDeriv2)
      if (nor.eq.0) goto 1110
!
!  Sqrt distances
!
      do k = 1,nor
        dist(k) = sqrt(dist(k))
      enddo
!*************************************************
!  Calculate twobody contribution in real space  *
!*************************************************
      eatoml = 0.0_dp
      ereall = 0.0_dp
      ec6l = 0.0_dp
      call twobody(eatoml,ereall,ec6l,lgrad1p,lgrad2,.false.,nor,1_i4,npots,npotl,cut2r, &
                   cut2q,cut2s,nmolonly,factor,ofct,ospfct,radsum,sctrm1,sctrm2,qli,qlj, &
                   lcspair,lewaldtype,.false.,.false.,lQMMMelectro)
!
      eatoml = ofctij*eatoml
      ereall = ofctij*ereall
      ec6l   = ofctij*ec6l
!
      esum = eatoml + ereall + ec6l
      if (lreg2one) then
        esregion12 = esregion12 + esum
      elseif (lreg2pair) then
        esregion2 = esregion2 + esum
      else
        eatom = eatom + eatoml
        ereal = ereal + ereall
        ec6 = ec6 + ec6l
      endif
      if (lattach) eattach = eattach + esum
!
      eregion2region(nregionj,nregioni) = eregion2region(nregionj,nregioni) + esum
!
      siteenergy(i) = siteenergy(i) + esum
!
      if (lPrintTwo) then
        write(ioout,'(4x,i12,i12,f22.10,1x,f22.10)') i,j,eatoml+ec6l,ereall
      endif
!
      if ((lDoQDeriv1.or.lDoQDeriv2).and.lgrad1) then
        do k = 1,nor
          d0i(k) = d0i(k) + derive0(k)*qlj
          d0j(k) = d0j(k) + derive0(k)*qli
        enddo
      endif
      if (leem) then
        if (lqeq.or.lSandM) then
          eqeql = 0.0_dp
          if (lqeq) then
            call qeqbody(eqeql,lgrad1,lgrad2,nor,1_i4,fct,qli,qlj,nati,natj,nqri,nqrj)
          elseif (lSandM) then
            call smbody(eqeql,lgrad1,lgrad2,nor,1_i4,fct,qli,qlj,nati,natj,nqri,nqrj)
          endif
          esum = eqeql
          if (lreg2one) then
            esregion12 = esregion12 + esum
          elseif (lreg2pair) then
            esregion2 = esregion2 + esum 
          else
            eqeq = eqeq + eqeql
          endif
          if (lattach) eattach = eattach + esum
!
          eregion2region(nregionj,nregioni) = eregion2region(nregionj,nregioni) + esum
!
          siteenergy(i) = siteenergy(i) + esum
        endif
      endif
      if (lsuttonc) then
        if (.not.lMEAMden) then
          if (lorder12loc) then
            scrho(1,i) = scrho(1,i) + sctrm1*ocj
            if (lattach) then
              scrho12(1,i) = scrho12(1,i) + sctrm1*ocj
            endif
          else
            scrho(1,i) = scrho(1,i) + sctrm2*ocj
            if (lattach) then
              scrho12(1,i) = scrho12(1,i) + sctrm2*ocj
            endif
          endif
        endif
      endif
!
!  Generate products for derivatives
!
      if (lmdl.or.lsg1.or.lgrad2) then
        call twostrterms(ndim,maxdis,nor,xtmp,ytmp,ztmp,xcom,ycom,zcom,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,lgrad2)
      endif
!*****************************
!  Charge first derivatives  *
!*****************************
      if (lgrad1.and.lDoQDeriv1) then
        call d1charged(ii,i,lopi,nor,d0i)
!
!  For self term this needs to be applied twice
!
        if (i.eq.j) call d1charged(ii,i,lopi,nor,d0i)
      endif
!*************************************
!  Electrostatic potential on-sites  *
!*************************************
      if (lpolar) then
        do k = 1,nor
          vx(i) = vx(i) - qlj*derivqd(k)*xtmp(k)
          vy(i) = vy(i) - qlj*derivqd(k)*ytmp(k)
          vz(i) = vz(i) - qlj*derivqd(k)*ztmp(k)
        enddo
        if (lattach) then
          do k = 1,nor
            vx12(i) = vx12(i) - qlj*derivqd(k)*xtmp(k)
            vy12(i) = vy12(i) - qlj*derivqd(k)*ytmp(k)
            vz12(i) = vz12(i) - qlj*derivqd(k)*ztmp(k)
          enddo
        endif
      endif
!***********************
!  Radial derivatives  *
!***********************
      if (lgrad1) then
        if (lgrad2) then
          rt2 = 0.0_dp
          do k = 1,nor
            rt2 = rt2 + rtrm2(k)
          enddo
        endif
        if (radi.gt.0.0_dp) then
          do k = 1,nor
            raderv(i) = raderv(i) + rtrm1(k)
          enddo
          if (lgrad2) then
            derv2(indrif,indri) = derv2(indrif,indri) + rt2
          endif
        endif
        if (lgrad2) then
          if (radi.gt.0.0_dp.and.radj.gt.0.0_dp) then
            derv2(indrj,indri) = derv2(indrj,indri) + rt2
          endif
        endif
      endif
!************************
!  Internal Derivatives *
!************************
!
!  First derivatives
!
      if (lgrad1) then
        do k = 1,nor
          xdrv(i) = xdrv(i) - deriv(k)*xtmp(k)
          ydrv(i) = ydrv(i) - deriv(k)*ytmp(k)
          zdrv(i) = zdrv(i) - deriv(k)*ztmp(k)
        enddo
        if (lgrad2.and.radi.gt.0.0_dp) then
          do k = 1,nor
            derv2(indrif,ix) = derv2(indrif,ix) - rderiv(k)*xtmp(k)
            derv2(indrif,iy) = derv2(indrif,iy) - rderiv(k)*ytmp(k)
            derv2(indrif,iz) = derv2(indrif,iz) - rderiv(k)*ztmp(k)
          enddo
        endif
        if (nregioni.ne.nregionj) then
          do k = 1,nor
            xregdrv(nregioni) = xregdrv(nregioni) - deriv(k)*xtmp(k)
            yregdrv(nregioni) = yregdrv(nregioni) - deriv(k)*ytmp(k)
            zregdrv(nregioni) = zregdrv(nregioni) - deriv(k)*ztmp(k)
          enddo
        endif
      endif
!***********************
!  Second derivatives  *
!***********************
      if (lgrad2) then
!
!  Coordinate only
!
        do k = 1,nor
          derv2(jx,ix) = derv2(jx,ix) - deriv2(k)*d2r2dx2(k,1)
          derv2(jy,ix) = derv2(jy,ix) - deriv2(k)*d2r2dx2(k,6)
          derv2(jz,ix) = derv2(jz,ix) - deriv2(k)*d2r2dx2(k,5)
          derv2(jx,iy) = derv2(jx,iy) - deriv2(k)*d2r2dx2(k,6)
          derv2(jy,iy) = derv2(jy,iy) - deriv2(k)*d2r2dx2(k,2)
          derv2(jz,iy) = derv2(jz,iy) - deriv2(k)*d2r2dx2(k,4)
          derv2(jx,iz) = derv2(jx,iz) - deriv2(k)*d2r2dx2(k,5)
          derv2(jy,iz) = derv2(jy,iz) - deriv2(k)*d2r2dx2(k,4)
          derv2(jz,iz) = derv2(jz,iz) - deriv2(k)*d2r2dx2(k,3)
          derv2(jx,ix) = derv2(jx,ix) - deriv(k)
          derv2(jy,iy) = derv2(jy,iy) - deriv(k)
          derv2(jz,iz) = derv2(jz,iz) - deriv(k)
        enddo
!
!  Coordinate - radius mixed
!
        if (radi.gt.0.0_dp) then
          do k = 1,nor
            derv2(jx,indri) = derv2(jx,indri) + rderiv(k)*xtmp(k)
            derv2(jy,indri) = derv2(jy,indri) + rderiv(k)*ytmp(k)
            derv2(jz,indri) = derv2(jz,indri) + rderiv(k)*ztmp(k)
          enddo
        endif
      endif
!***********************
!  Strain derivatives  *
!***********************
      if (lsg2) then
!
!  Mixed derivatives
!
        do kl = 1,nstrains
          ks = nstrptr(kl)
          do k = 1,nor
            derv3(ix,kl) = derv3(ix,kl) - deriv(k)*d2r2dsdx(k,ks,1) - xtmp(k)*deriv2(k)*dr2ds(k,ks)
            derv3(iy,kl) = derv3(iy,kl) - deriv(k)*d2r2dsdx(k,ks,2) - ytmp(k)*deriv2(k)*dr2ds(k,ks)
            derv3(iz,kl) = derv3(iz,kl) - deriv(k)*d2r2dsdx(k,ks,3) - ztmp(k)*deriv2(k)*dr2ds(k,ks)
          enddo
        enddo
        if (radi.gt.0.0_dp) then
          do kl = 1,nstrains
            ks = nstrptr(kl)
            do k = 1,nor
              derv3(indri,kl) = derv3(indri,kl) + rderiv(k)*dr2ds(k,ks)
            enddo
          enddo
        endif
      endif
!
!  Strain only terms 
!
!  First derivatives 
!
      if (lsg1.or.lgrad2) then
        rstrdloc(1:nstrains) = 0.0_dp
        do kl = 1,nstrains
          ks = nstrptr(kl)
          do k = 1,nor
            rstrdloc(kl) = rstrdloc(kl) + deriv(k)*dr2ds(k,ks)
          enddo
          rstrd(kl) = rstrd(kl) + ofctij*rstrdloc(kl)
        enddo
        if (latomicstress) then
          if (i.eq.j) then
            do kl = 1,nstrains
              atomicstress(kl,i) = atomicstress(kl,i) + rstrdloc(kl)
            enddo
          else
            do kl = 1,nstrains
              atomicstress(kl,i) = atomicstress(kl,i) + 0.5_dp*rstrdloc(kl)
            enddo
          endif
        endif
!
!  Second derivatives
!
        if (lsg2) then
          do kk = 1,nstrains
            ks = nstrptr(kk)
            do kl = 1,nstrains
              kt = nstrptr(kl)
              do k = 1,nor
                sderv2(kl,kk) = sderv2(kl,kk) + ofctij*deriv2(k)*dr2ds(k,kt)*dr2ds(k,ks)
                sderv2(kl,kk) = sderv2(kl,kk) + ofctij*deriv(k)*d2r2ds2(k,kt,ks)
              enddo
            enddo
          enddo
        endif
      endif
1110  continue
!*******************************************************
!  Variable charge contribution to second derivatives  *
!*******************************************************
      if (lgrad2.and.lDoQDeriv2) then
        derive0selfi = derive0self*qlj
        derive0selfj = derive0self*qli
        d1ix = 0.0_dp
        d1iy = 0.0_dp
        d1iz = 0.0_dp
        do k = 1,nor
          d1ix = d1ix + d1i(k)*xtmp(k)
          d1iy = d1iy + d1i(k)*ytmp(k)
          d1iz = d1iz + d1i(k)*ztmp(k)
        enddo
        if (lstr) then
          do kl = 1,nstrains
            ks = nstrptr(kl)
            d1is(kl) = 0.0_dp
            do k = 1,nor
              d1is(kl) = d1is(kl) + d1i(k)*dr2ds(k,ks)
            enddo
          enddo
        endif
        call d2charge(i,j,nor,ix,iy,iz,jx,jy,jz,lopi,lopj,d0i,d0j,d1ix,d1iy,d1iz, &
                      d1jx,d1jy,d1jz,d1is,d1js,d2i2,d2ij,d2j2,d2self,derive0selfi, &
                      derive0selfj,.true.,.true.)
      endif
      if (lbsmat(nsft+nrelf2a(i))) then
!*********************
!  Radial self term  *
!*********************
        eatoml = 0.0_dp
        do m = 1,npote
          if (nptype(m).eq.14) then
            if (lmatch(nati,ntypi,nspec1(m),nptyp1(m),.true.)) then
              rdiff = radi - twopot(2,m)
              apt = twopot(1,m)*oci
              eatoml = eatoml + 0.5_dp*apt*rdiff*rdiff
              if (ldog1) then
                raderv(i) = raderv(i) + apt*rdiff
                if (ldog2) then
                  derv2(indrif,indri) = derv2(indrif,indri) + apt
                endif
              endif
            endif
          elseif (nptype(m).eq.17) then
            if (lmatch(nati,ntypi,nspec1(m),nptyp1(m),.true.)) then
              rdiff = radi - twopot(3,m)
              apt = twopot(1,m)*oci
              bpt = twopot(2,m)
              etrm1 = exp(bpt*rdiff)
              etrm2 = 1.0_dp/etrm1
              etrm = apt*(etrm1 + etrm2)
              eatoml = eatoml + etrm
              if (ldog1) then
                raderv(i) = raderv(i) + apt*bpt*(etrm1 - etrm2)
                if (ldog2) then
                  derv2(indrif,indri) = derv2(indrif,indri) + bpt*bpt*etrm
                endif
              endif
            endif
          elseif (nptype(m).eq.31) then
            if (lmatch(nati,ntypi,nspec1(m),nptyp1(m),.true.)) then
              rdiff = radi - twopot(3,m)
              apt = twopot(1,m)*oci
              bpt = twopot(2,m)
              etrm1 = exp(bpt*rdiff)
              etrm = apt*etrm1
              eatoml = eatoml + etrm
              if (ldog1) then
                raderv(i) = raderv(i) + apt*bpt*etrm1
                if (ldog2) then
                  derv2(indrif,indri) = derv2(indrif,indri) + bpt*bpt*etrm
                endif
              endif
            endif
          endif
        enddo
!
!  Add energy to total
!
        eatom = eatom + eatoml*ofctij
      endif
!
    enddo jloop
  enddo
  if (nprocs.gt.1) then
!****************
!  Global sums  *
!****************
    if (lsuttonc) then
      tsum0 = g_cpu_time()
!
      allocate(sum(numat),stat=status)
      if (status/=0) call outofmemory('realed','sum')
      allocate(sum2(numat),stat=status)
      if (status/=0) call outofmemory('realed','sum2')
!
      if (.not.lMEAMden) then
        do i = 1,numat
          sum2(i) = scrho(1,i)
        enddo
        call sumall(sum2,sum,numat,"realed","scrho")
        do i = 1,numat
          scrho(1,i) = sum(i)
        enddo
        do i = 1,numat
          sum2(i) = scrho12(1,i)
        enddo
        call sumall(sum2,sum,numat,"realed","scrho12")
        do i = 1,numat
          scrho12(1,i) = sum(i)
        enddo
      endif
!
      deallocate(sum2,stat=status)
      if (status/=0) call deallocate_error('realed','sum2')
      deallocate(sum,stat=status)
      if (status/=0) call deallocate_error('realed','sum')
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
      write(ioout,'(''--------------------------------------------------------------------------------'',/)')
    endif
  endif
!
!  End of real space part - perform general tasks
!
999 continue
!
!  Free local memory
!
  deallocate(npotl,stat=status)
  if (status/=0) call deallocate_error('realed','npotl')
!
!  Timing
!
  time2 = g_cpu_time()
  tatom = tatom + time2 - time1 - tsuml
#ifdef TRACE
  call trace_out('realed')
#endif
!
  return
  end
