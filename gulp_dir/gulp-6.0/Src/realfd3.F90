  subroutine realfd3(matom,eatom,ereal,erecip,ec6,eqeq,lgrad1)
!
!  Subroutine for calculating real space energy and gradients.
!  Finite difference version that focuses on derivatives of matom.
!
!  12/17 Created from realmd3
!   2/18 Trace added
!   5/18 Modified for q ranges
!   8/18 Call to twostrterms introduced for setting of rpd
!   9/18 Strain module introduced
!   8/19 Corrections to arguments in call to twobody
!   8/19 Trap for neemrptr being zero added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!  11/19 Rigid molecule modifications added
!   3/20 Electrostatic cutoff only included where the either charge exceeds the threshold
!   4/20 derv3c and d2r2dsdc added for benefit of rigid molecules
!   4/20 derv3c changes reversed as they are no longer required
!   7/20 Modifications for gfortran v10
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
!  Julian Gale, CIC, Curtin University, July 2020
!
  use configurations, only : lbsmat,nregionno,nregiontype,QMMMmode
  use g_constants
  use control
  use current
  use datatypes
  use derivatives
  use eam,            only : lMEAMden
  use eemdata
  use element
  use general,        only : cutw
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
  integer(i4), intent(in)                      :: matom
  logical,     intent(in)                      :: lgrad1
  real(dp),    intent(inout)                   :: eatom
  real(dp),    intent(inout)                   :: ec6
  real(dp),    intent(inout)                   :: eqeq
  real(dp),    intent(inout)                   :: ereal
  real(dp),    intent(inout)                   :: erecip
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: indm
  integer(i4)                                  :: indmj
  integer(i4)                                  :: ix
  integer(i4)                                  :: ixi
  integer(i4)                                  :: iyi
  integer(i4)                                  :: izi
  integer(i4)                                  :: ixj
  integer(i4)                                  :: iyj
  integer(i4)                                  :: izj
  integer(i4)                                  :: j
  integer(i4)                                  :: jx
  integer(i4)                                  :: k
  integer(i4)                                  :: kl
  integer(i4)                                  :: ks
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
  logical                                      :: lc6loc
  logical                                      :: lcspair
  logical                                      :: ldog1
  logical                                      :: ldog2
  logical                                      :: lewaldtype
  logical                                      :: lmatch
  logical                                      :: lmdl
  logical                                      :: lmolok
  logical                                      :: lneedmol
  logical                                      :: lorder12loc
  logical                                      :: lQMMMelectro
  logical                                      :: lself
  logical                                      :: lsg1  
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
  real(dp)                                     :: d2self
  real(dp)                                     :: derive0self
  real(dp)                                     :: dqi(1)
  real(dp)                                     :: dqj(1)
  real(dp)                                     :: eatom2
  real(dp)                                     :: ec62 
  real(dp)                                     :: eqeq2  
  real(dp)                                     :: ereal2
  real(dp)                                     :: erecip2
  real(dp)                                     :: etrm 
  real(dp)                                     :: etrm1
  real(dp)                                     :: etrm2
  real(dp)                                     :: factor
  real(dp)                                     :: fct
  real(dp)                                     :: oci
  real(dp)                                     :: ocj
  real(dp)                                     :: ofct
  real(dp)                                     :: ospfct
  real(dp)                                     :: qli
  real(dp)                                     :: qlj
  real(dp)                                     :: radi
  real(dp)                                     :: radj
  real(dp)                                     :: radsum
  real(dp)                                     :: rdiff
  real(dp)                                     :: rp
  real(dp)                                     :: rqeq2
  real(dp)                                     :: rstrdloc(6)
  real(dp)                                     :: sctrm1
  real(dp)                                     :: sctrm2
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
  call trace_in('realfd3')
#endif
!
  time1 = g_cpu_time()
!
!  Local variables
!
  lmdl = lmd
  if (.not.lgrad1) lmdl = .false.
  lc6loc = (lc6.and.ndim.eq.3)
  rqeq2 = rqeq*rqeq
!
!  Set the Coulomb term type based on dimensionality :
!
!  1-D         => 1/r
!  2-D and 3-D => erfc(seta*r)/r
!
  lewaldtype = (ndim.ne.1.or.lwolf)
  lsg1 = (lgrad1.and.lstr)
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
  allocate(npotl(npote),stat=status)
  if (status/=0) call outofmemory('realfd3','npotl')
!
  if (.not.lnoreal) then
!***************************************************************
!  Atomistic and real space electrostatic component of energy  *
!***************************************************************
!
!  Only compute interactions with matom
!
    i = matom
    xal = xclat(i)
    yal = yclat(i)
    zal = zclat(i)
    nati = nat(i)
    ntypi = nftype(i)
    qli = qf(i)
    oci = occuf(i)
    nregioni = nregionno(nsft+nrelf2a(i))
    nregiontypi = nregiontype(nregioni,ncf)
    if (lbsmat(nrelf2a(i)+nsft)) then
      radi = radf(i)
    else
      radi = 0.0_dp
    endif
!
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
    jloop: do j = 1,numat
!
!  Exclude self term
!
      if (i.eq.j) cycle jloop
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
              elseif (nptype(n).eq.57) then
                c6tot = c6tot + twopot(3,n)*twopot(4,n)*twopot(5,n)**6
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
      if (npots.eq.0.and.abs(factor).lt.1.0d-8) cycle jloop
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
      nmolonly = 0
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
      ereal2 = 0.0_dp
      erecip2 = 0.0_dp
      ec62 = 0.0_dp
      derive0self = 0.0_dp
      if (lself) call selfterm(ereal2,erecip2,ec62,derive0self,factor,fct,ofct,ospfct,1.0_dp,npotl,npots, &
                               c6tot,d2self,lgrad1,.false.,i,j,ix,jx,1.0_dp,0.5_dp,lewaldtype,qli,qlj)
      ereal = ereal + ereal2
      erecip = erecip + erecip2
      ec6 = ec6 + ec62
!*******************************************
!  Charge first derivatives for self-term  *
!*******************************************
      if (lself.and.lgrad1.and.lDoQDeriv1) then
        dqi(1) = derive0self*qli
        dqj(1) = derive0self*qlj
        call d1charge(i,j,.true.,.true.,1_i4,dqj,dqi)
      endif
!
      if (nor.eq.0) cycle jloop
!
!  Sqrt distances
!
      do k = 1,nor
        dist(k) = sqrt(dist(k))
      enddo
!*************************************************
!  Calculate twobody contribution in real space  *
!*************************************************
      eatom2 = 0.0_dp
      ereal2 = 0.0_dp
      ec62 = 0.0_dp
      call twobody(eatom2,ereal2,ec62,lgrad1,.false.,.false.,nor,1_i4,npots,npotl,cut2r, &
                   cut2q,cut2s,nmolonly,factor,ofct,ospfct,radsum,sctrm1,sctrm2,qli,qlj, &
                   lcspair,lewaldtype,.false.,.false.,lQMMMelectro)
      eatom = eatom + eatom2
      ereal = ereal + ereal2
      ec6 = ec6 + ec62
!
      if ((lDoQDeriv1.or.lDoQDeriv2).and.lgrad1) then
        do k = 1,nor
          d0i(k) = d0i(k) + derive0(k)*qlj
          d0j(k) = d0j(k) + derive0(k)*qli
        enddo
      endif
      if (leem) then
        if (lqeq.or.lSandM) then
          eqeq2 = 0.0_dp
          if (lqeq) then
            call qeqbody(eqeq2,lgrad1,.false.,nor,1_i4,fct,qli,qlj,nati,natj,nqri,nqrj)
          elseif (lSandM) then
            call smbody(eqeq2,lgrad1,.false.,nor,1_i4,fct,qli,qlj,nati,natj,nqri,nqrj)
          endif
          eqeq = eqeq + eqeq2
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
!
!  Generate products for derivatives
!
      if (lmdl.or.lsg1) then
        call twostrterms(ndim,maxdis,nor,xtmp,ytmp,ztmp,xcom,ycom,zcom,dr2ds,rpd,d2r2dsdx,d2r2ds2,.false.)
      endif
!*****************************
!  Charge first derivatives  *
!*****************************
      if (lgrad1.and.lDoQDeriv1) then
        call d1charge(i,j,.true.,.true.,nor,d0i,d0j)
      endif
!***********************
!  Radial derivatives  *
!***********************
      if (lgrad1) then
        if (radi.gt.0.0_dp) then
          do k = 1,nor
            raderv(i) = raderv(i) + rtrm1(k)
          enddo
        endif
        if (radj.gt.0.0_dp) then
          do k = 1,nor
            raderv(j) = raderv(j) + rtrm1(k)
          enddo
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
          xdrv(j) = xdrv(j) + deriv(k)*xtmp(k)
          ydrv(j) = ydrv(j) + deriv(k)*ytmp(k)
          zdrv(j) = zdrv(j) + deriv(k)*ztmp(k)
        enddo
      endif
!***********************
!  Strain derivatives  *
!***********************
!
!  Strain only terms 
!
!  First derivatives 
!
      if (lsg1) then
        rstrdloc(1:nstrains) = 0.0_dp
        do kl = 1,nstrains
          ks = nstrptr(kl)
          do k = 1,nor
            rstrdloc(kl) = rstrdloc(kl) + deriv(k)*dr2ds(k,ks)
          enddo
          rstrd(kl) = rstrd(kl) + rstrdloc(kl)
        enddo
      endif
    enddo jloop
!**************
!  Self-term  *
!**************
    if (lbsmat(nrelf2a(i)+nsft)) then
      radi = radf(i)
      radsum = 2.0_dp*radi
    else
      radi = 0.0_dp
      radsum = 0.0_dp
    endif
!
!  Factor of half for self term
!
    ofct = 0.5_dp*oci*oci
    ospfct = ofct
    fct = ofct*angstoev
    factor = qli*qli*fct
!
!  QM/MM handling : i is a QM atom => exclude
!
    if (QMMMmode(ncf).eq.0.or.nregiontypi.ne.1) then
!           
!  QM/MM : Set electrostatic embedding flag : If i is a QM atom => exclude electrostatics
!       
      lQMMMelectro = (QMMMmode(ncf).eq.2.and.nregiontypi.eq.1)
!
!  Molecule handling
!
      if (lmol) then
        ixj = 0
        iyj = 0
        izj = 0
        nmi = natmol(i)
        lmolok = (nmi.gt.0)
      else
        lmolok = .false.
      endif
!
!  COM - zero for self term
!
      xcom = 0.0_dp
      ycom = 0.0_dp
      zcom = 0.0_dp
!
!  Locate potential number
!  Check whether potential requires specific types
!  Calculate sum of dispersion terms
!
      rp = 0.0_dp
      npots = 0
      c6tot = 0.0_dp
      lneedmol = (lmol.and..not.lmolq)
      do n = 1,npote
        if (lmatch(nati,ntypi,nspec1(n),nptyp1(n),.true.)) then
          if (lmatch(nati,ntypi,nspec2(n),nptyp2(n),.true.)) then
            npots = npots + 1
            npotl(npots) = n
            if (mmexc(n).gt.0.or.(lintra(n).and..not.linter(n).or..not.lintra(n).and.linter(n))) lneedmol = .true.
            if (rpot(n).gt.rp) rp = rpot(n)
            if (nptype(n).eq.1.or.nptype(n).eq.7) then
              c6tot = c6tot + twopot(3,n)
            elseif (nptype(n).eq.2.or.nptype(n).eq.21) then
              c6tot = c6tot + twopot(2,n)
            elseif (nptype(n).eq.57) then
              c6tot = c6tot + twopot(3,n)*twopot(4,n)*twopot(5,n)**6
            endif
          endif
        endif
      enddo
      cut2r = rp*rp
      if (cut2r.gt.cut2p) cut2r = cut2p
      cut2 = cut2r
!
!  If both charges are less than threshold then exclude electrostatics from cutoff
!
      if (abs(qli*oci).gt.thresh_q) then
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
        call rsearch3D(0.0_dp,0.0_dp,0.0_dp,lmolok,lcspair,i,i,nmi,ixj,iyj,izj,nor,nmolonly,lself,cut2)
      elseif (ndim.eq.2) then
        call rsearch2D(0.0_dp,0.0_dp,0.0_dp,lmolok,lcspair,i,i,nmi,ixj,iyj,izj,nor,nmolonly,lself,cut2)
      elseif (ndim.eq.1) then
        call rsearch1D(0.0_dp,0.0_dp,0.0_dp,lmolok,lcspair,i,i,nmi,ixj,iyj,izj,nor,nmolonly,lself,cut2)
      endif
!
      if (lself) then
        ereal2 = 0.0_dp
        erecip2 = 0.0_dp
        ec62 = 0.0_dp
        derive0self = 0.0_dp
        call selfterm(ereal2,erecip2,ec62,derive0self,factor,fct,ofct,ospfct,1.0_dp,npotl,npots, &
                      c6tot,d2self,lgrad1,.false.,i,i,ix,ix,1.0_dp,0.5_dp,lewaldtype,qli,qlj)
        ereal = ereal + ereal2
        erecip = erecip + erecip2
        ec6 = ec6 + ec62
      endif
!*******************************************
!  Charge first derivatives for self-term  *
!*******************************************
      if (lself.and.lgrad1.and.lDoQDeriv1) then
        dqi(1) = derive0self*qli
        call d1charge(i,i,.true.,.true.,1_i4,dqi,dqi)
      endif
!
      if (nor.gt.0) then
!
        ldog1 = (lmdl.or.lsg1.or.lbsmat(nsft+nrelf2a(i)).or.leem)
        ldog2 = (lbsmat(nsft+nrelf2a(i)).or.leem)
!
!  Sqrt distances
!
        do k = 1,nor
          dist(k) = sqrt(dist(k))
        enddo
!*************************************************
!  Calculate twobody contribution in real space  *
!*************************************************
        eatom2 = 0.0_dp
        ereal2 = 0.0_dp
        ec62 = 0.0_dp
        call twobody(eatom2,ereal2,ec62,ldog1,ldog2,.false.,nor,1_i4,npots,npotl,cut2r,cut2q,cut2s,nmolonly, &
                     factor,ofct,ospfct,radsum,sctrm1,sctrm2,qli,qli,.false.,lewaldtype,.false.,.false., &
                     lQMMMelectro)
        eatom = eatom + eatom2
        ereal = ereal + ereal2
        ec6 = ec6 + ec62
!
        if (lDoQDeriv1.and.lgrad1) then
          do k = 1,nor
            d0i(k) = d0i(k) + derive0(k)*qli
            d0j(k) = d0j(k) + derive0(k)*qli
          enddo
        endif
        if (leem.and.ldog1) then
          if (lqeq.or.lSandM) then
            eqeq2 = 0.0_dp
            if (lqeq) then
              call qeqbody(eqeq2,ldog1,ldog2,nor,1_i4,fct,qli,qli,nati,nati,nqri,nqri)
            elseif (lSandM) then
              call smbody(eqeq2,ldog1,ldog2,nor,1_i4,fct,qli,qli,nati,nati,nqri,nqri)
            endif
            eqeq = eqeq + eqeq2
          endif
        endif
!
!  Many-body contribution - multiply by 2 to correct for factor of 1/2 in ofct which isn't needed here
!
        if (lsuttonc) then
          if (.not.lMEAMden) then
            scrho(1,i) = scrho(1,i) + sctrm1*oci
          endif
        endif
!
!  Generate products for derivatives
!
        if (lmdl.or.lsg1) then
          call twostrterms(ndim,maxdis,nor,xtmp,ytmp,ztmp,xcom,ycom,zcom,dr2ds,rpd,d2r2dsdx,d2r2ds2,.false.)
        endif
!*****************************
!  Charge first derivatives  *
!*****************************
        if (lgrad1.and.lDoQDeriv1) then
          call d1charge(i,i,.true.,.true.,nor,d0i,d0i)
        endif
      endif
      if (lbsmat(nsft+nrelf2a(i))) then
!***********************
!  Radial derivatives  *
!***********************
!
!  Find self term
!
        eatom2 = 0.0_dp
        do m = 1,npote
          if (nptype(m).eq.14) then
            if (lmatch(nati,ntypi,nspec1(m),nptyp1(m),.true.)) then
              rdiff = radi - twopot(2,m)
              apt = twopot(1,m)*oci
              etrm = 0.5_dp*apt*rdiff*rdiff
              eatom2 = eatom2 + etrm
              if (lgrad1) then
                raderv(i) = raderv(i) + apt*rdiff
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
              eatom2 = eatom2 + etrm
              if (lgrad1) then
                raderv(i) = raderv(i) + apt*bpt*(etrm1 - etrm2)
              endif
            endif
          elseif (nptype(m).eq.31) then
            if (lmatch(nati,ntypi,nspec1(m),nptyp1(m),.true.)) then
              rdiff = radi - twopot(3,m)
              apt = twopot(1,m)*oci
              bpt = twopot(2,m)
              etrm1 = exp(bpt*rdiff)
              etrm = apt*etrm1
              eatom2 = eatom2 + etrm
              if (lgrad1) then
                raderv(i) = raderv(i) + apt*bpt*etrm1
              endif
            endif
          endif
        enddo
        eatom = eatom + eatom2
        if (lgrad1) then
          do k = 1,nor
            raderv(i) = raderv(i) + 2.0_dp*rtrm1(k)
          enddo
        endif
      endif
!***********************
!  Strain derivatives  *
!***********************
      if (lsg1) then
        rstrdloc(1:nstrains) = 0.0_dp
        do kl = 1,nstrains
          ks = nstrptr(kl)
          do k = 1,nor
            rstrdloc(kl) = rstrdloc(kl) + deriv(k)*dr2ds(k,ks)
          enddo
          rstrd(kl) = rstrd(kl) + rstrdloc(kl)
        enddo
      endif
    endif
!
!  Exit point
!
  endif
!
!  Free local memory
!
  deallocate(npotl,stat=status)
  if (status/=0) call deallocate_error('realfd3','npotl')
!
!  Timing
!
  time2 = g_cpu_time()
  tatom = tatom + time2 - time1
#ifdef TRACE
  call trace_out('realfd3')
#endif
!
  return
  end
