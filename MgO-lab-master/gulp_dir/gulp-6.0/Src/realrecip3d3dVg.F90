  subroutine realrecip3d3dVg(nkp,mcv,eigr,maxd2l,iocptr,grueneisen)
!
!  Supplies the first derivatives of the phonon frequencies for a solid with respect to volume
!  Gamma point version.
!
!  NOTE : Needs finishing for BSM 
!
!   1/18 Created from realrecip3d3dV
!   2/18 Trace added
!   8/19 Modifications for Intel compiler
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!   3/20 Electrostatic cutoff only included where the either charge exceeds the threshold
!   4/20 Mass arrays now use fmass and rfmass
!
!  nmanyk = number of three-/four-body term K atoms whose derivatives
!          depend on the i-j block
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
  use configurations, only : lbsmat
  use g_constants
  use control
  use current
  use datatypes
  use derivatives
  use element,        only : maxele
  use feworkspace
  use four
  use frequencies,    only : rfmass
  use general
  use iochannels
  use ksample
  use kspace
  use m_three
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
  integer(i4),                    intent(in)   :: iocptr(*)
  integer(i4),                    intent(in)   :: maxd2l
  integer(i4),                    intent(in)   :: mcv
  integer(i4),                    intent(in)   :: nkp
  real(dp),                       intent(in)   :: eigr(maxd2l,*)
  real(dp),                       intent(out)  :: grueneisen(*)
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ii
  integer(i4)                                  :: iimn
  integer(i4)                                  :: iimx
  integer(i4)                                  :: indi
  integer(i4)                                  :: indj
  integer(i4)                                  :: indm
  integer(i4)                                  :: indmj
  integer(i4)                                  :: ioc
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
  integer(i4)                                  :: jj
  integer(i4)                                  :: jjmn
  integer(i4)                                  :: jjmx
  integer(i4)                                  :: joc
  integer(i4)                                  :: jx
  integer(i4)                                  :: jy
  integer(i4)                                  :: jz
  integer(i4)                                  :: kk
  integer(i4)                                  :: kkmn
  integer(i4)                                  :: kkmx
  integer(i4)                                  :: max1l
  integer(i4)                                  :: max1l1
  integer(i4)                                  :: max1u
  integer(i4)                                  :: max2l
  integer(i4)                                  :: max2l1
  integer(i4)                                  :: max2u
  integer(i4)                                  :: max3l
  integer(i4)                                  :: max3l1
  integer(i4)                                  :: max3u
  integer(i4)                                  :: maxlim
  integer(i4)                                  :: mint
  integer(i4)                                  :: n
  integer(i4)                                  :: naddx
  integer(i4)                                  :: naddy
  integer(i4)                                  :: naddz
  integer(i4)                                  :: nat1
  integer(i4)                                  :: nat2
  integer(i4)                                  :: nati
  integer(i4)                                  :: natj
  integer(i4)                                  :: nforkl
  integer(i4)                                  :: nkvec0
  integer(i4)                                  :: nmanyk
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
  integer(i4)                                  :: nwords
  integer(i4)                                  :: status
  logical                                      :: l111
  logical                                      :: lcspair
  logical                                      :: lc6loc
  logical                                      :: lcorei
  logical                                      :: lcorej
  logical                                      :: lewaldtype
  logical                                      :: lfor
  logical                                      :: lthb
  logical                                      :: lmany
  logical                                      :: lmatch
  logical                                      :: lmolok
  logical                                      :: lnadd
  logical                                      :: lneedmol 
  logical                                      :: lorder12loc
  logical                                      :: lsamemol
  real(dp)                                     :: c6tot
  real(dp)                                     :: cmax
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: cut2
  real(dp)                                     :: cut2e
  real(dp)                                     :: cut2p
  real(dp)                                     :: cut2q
  real(dp)                                     :: cut2r
  real(dp)                                     :: cut2s
  real(dp)                                     :: cut2w
  real(dp)                                     :: d3k(3,3,3)
  real(dp)                                     :: d3s(3,3,3)
  real(dp)                                     :: d3ks(3,3,3)
  real(dp)                                     :: eatom
  real(dp)                                     :: ec6
  real(dp)                                     :: emax
  real(dp)                                     :: ereal
  real(dp)                                     :: factor
  real(dp),    dimension(:), allocatable       :: ktrm6
  real(dp),    dimension(:), allocatable       :: ktrm0
  real(dp),    dimension(:), allocatable       :: ktrm20
  real(dp),    dimension(:), allocatable       :: ktrm60
  real(dp),    dimension(:), allocatable       :: ktrm62
  real(dp),    dimension(:), allocatable       :: ktrm620
  real(dp)                                     :: oci
  real(dp)                                     :: ocj
  real(dp)                                     :: ofct
  real(dp)                                     :: ospfct
  real(dp)                                     :: qli
  real(dp)                                     :: qlj
  real(dp)                                     :: perpr
  real(dp)                                     :: projr
  real(dp)                                     :: r2
  real(dp)                                     :: ra
  real(dp)                                     :: radi
  real(dp)                                     :: radj
  real(dp)                                     :: radsum
  real(dp)                                     :: rb
  real(dp)                                     :: rc
  real(dp)                                     :: rdl1
  real(dp)                                     :: rdl2
  real(dp)                                     :: rdl3
  real(dp)                                     :: rdl4
  real(dp)                                     :: rdl5
  real(dp)                                     :: rdl6
  real(dp)                                     :: rmassi
  real(dp)                                     :: rmassj
  real(dp)                                     :: rp
  real(dp)                                     :: rpres
  real(dp)                                     :: rpres2
  real(dp)                                     :: rpres3
  real(dp)                                     :: ru1x
  real(dp)                                     :: ru1y
  real(dp)                                     :: ru1z
  real(dp)                                     :: ru2x
  real(dp)                                     :: ru2y
  real(dp)                                     :: ru2z
  real(dp)                                     :: ru3x
  real(dp)                                     :: ru3y
  real(dp)                                     :: ru3z
  real(dp)                                     :: sctrm1
  real(dp)                                     :: sctrm2
  real(dp)                                     :: small
  real(dp)                                     :: t1
  real(dp)                                     :: t2
  real(dp)                                     :: time1
  real(dp)                                     :: time2
  real(dp)                                     :: tproj0
  real(dp)                                     :: xal
  real(dp)                                     :: yal
  real(dp)                                     :: zal
  real(dp)                                     :: xcd
  real(dp)                                     :: ycd
  real(dp)                                     :: zcd
  real(dp)                                     :: xcd2
  real(dp)                                     :: ycd2
  real(dp)                                     :: xcdi
  real(dp)                                     :: ycdi
  real(dp)                                     :: zcdi
  real(dp)                                     :: xcdj
  real(dp)                                     :: ycdj
  real(dp)                                     :: zcdj
  real(dp)                                     :: xcrd
  real(dp)                                     :: ycrd
  real(dp)                                     :: zcrd
  real(dp)                                     :: xd2
  real(dp)                                     :: yd2
  real(dp)                                     :: zd2
  real(dp)                                     :: xd3
  real(dp)                                     :: yd3
  real(dp)                                     :: zd3
#ifdef TRACE
  call trace_in('realrecip3d3dVg')
#endif
!
  time1 = g_cpu_time()
  tproj0 = tproj
  lc6loc = (lc6.and.ndim.eq.3)
!
!  Initialise volume derivatives of frequencies
!
  do i = 1,mcv
    grueneisen(i) = 0.0_dp
  enddo
!
!  Local variables
!
  small = 1.0d-12
  lthb = (nthb.gt.0)
  lfor = (nfor.gt.0)
  lmany = (lthb.or.lfor)
!  
!  Set the Coulomb term type based on dimensionality :
!     
!  1-D         => 1/r
!  2-D and 3-D => erfc(seta*r)/r 
!
  lewaldtype = (ndim.ne.1.or.lwolf)
!
  mint = 3*numat
  maxlim = mint
  if (nbsmat.gt.0) maxlim = maxlim + numat
!
!  Zero energies although not needed to avoid overflow
!
  eatom = 0.0_dp
  ereal = 0.0_dp
  ec6 = 0.0_dp
!
  lnadd = .true.
  if (nadd.eq.0) then
    lnadd = .false.
    if (lra) then
      nadd = 1
    else
      if (ndim.eq.3) then
        if (alpha.lt.30.0_dp.or.beta.lt.30.0_dp.or.gamma.lt.30.0_dp) then
          nadd = 5   
        elseif (alpha.gt.150.0_dp.or.beta.gt.150.0_dp.or.gamma.gt.150.0_dp) then
          nadd = 5
        elseif (alpha.lt.50.0_dp.or.beta.lt.50.0_dp.or.gamma.lt.50.0_dp) then
          nadd = 4
        elseif (alpha.gt.130.0_dp.or.beta.gt.130.0_dp.or.gamma.gt.130.0_dp) then
          nadd = 4
        elseif (alpha.lt.70.0_dp.or.beta.lt.70.0_dp.or.gamma.lt.70.0_dp) then
          nadd = 3
        elseif (alpha.gt.110.0_dp.or.beta.gt.110.0_dp.or.gamma.gt.110.0_dp) then
          nadd = 3
        else
          nadd = 2
        endif
      elseif (ndim.eq.2) then
        if (alpha.lt.30.0_dp.or.alpha.gt.150.0_dp) then
          nadd = 4
        elseif (alpha.lt.50.0_dp.or.alpha.gt.130.0_dp) then
          nadd = 3
        elseif (alpha.lt.70.0_dp.or.alpha.gt.110.0_dp) then
          nadd = 2
        else
          nadd = 1
        endif
      elseif (ndim.eq.1) then
        nadd = 2
      endif
    endif
  endif
!***************************
!  Set up local variables  *
!***************************
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
!  Decide whether all interactions lie within unit cell and first
!  neighbours - saves time for large systems
!
  emax = sqrt(cut2q)
  cmax = max(rpmax,emax)
  l111 = .true.
  if (ndim.eq.3) then
    ra = 1.0_dp/a
    rb = 1.0_dp/b
    rc = 1.0_dp/c
    if (a.lt.cmax.or.b.lt.cmax.or.c.lt.cmax) l111 = .false.
    if (alpha.lt.80.0_dp.or.beta.lt.80.0_dp.or.gamma.lt.80.0_dp) l111 = .false.
    if (alpha.gt.100.0_dp.or.beta.gt.100.0_dp.or.gamma.gt.100.0_dp) l111 = .false.
    if (l111) then
      iimn = 1
      iimx = 1
      jjmn = 1
      jjmx = 1
      kkmn = 1
      kkmx = 1
    else
      naddx = nadd
      naddy = nadd
      naddz = 1
    endif
  elseif (ndim.eq.2) then
    ra = 1.0_dp/a
    rb = 1.0_dp/b
    rc = 0.0_dp
    if (a.lt.cmax.or.b.lt.cmax) l111 = .false.
    if (alpha.lt.80.0_dp.or.alpha.gt.100.0_dp) l111 = .false.
    if (l111) then
      iimn = 1
      iimx = 1
      jjmn = 1
      jjmx = 1
      kkmn = 0
      kkmx = 0
    else
      naddx = nadd
      naddy = nadd
      naddz = 0
    endif
  elseif (ndim.eq.1) then
    ra = 1.0_dp/a
    rb = 0.0_dp
    rc = 0.0_dp
    if (a.lt.cmax) l111 = .false.
    if (l111) then
      iimn = 1
      iimx = 1
      jjmn = 0
      jjmx = 0
      kkmn = 0
      kkmx = 0
    else
      naddx = nadd
      naddy = 0
      naddz = 0
    endif
  endif
!
!  Create unit vectors
!
  ru1x = r1x*ra
  ru1y = r1y*ra
  ru1z = r1z*ra
  ru2x = r2x*rb
  ru2y = r2y*rb
  ru2z = r2z*rb
  ru3x = r3x*rc
  ru3y = r3y*rc
  ru3z = r3z*rc
!
!  Allocate local memory
!
  allocate(npotl(npote),stat=status)
  if (status/=0) call outofmemory('realrecip3d3dVg','npotl')
!**********************************************************
!  Calculate K vectors for reciprocal space contribution  *
!**********************************************************
!
!  Allocate memory
!
  nwords = maxkvec + 1
  allocate(ktrm0(nwords),stat=status)
  if (status/=0) call outofmemory('realrecip3d3dVg','ktrm0')
  allocate(ktrm20(nwords),stat=status)
  if (status/=0) call outofmemory('realrecip3d3dVg','ktrm20')
  if (lc6loc) then
    allocate(ktrm6(nwords),stat=status)
    if (status/=0) call outofmemory('realrecip3d3dVg','ktrm6')
    allocate(ktrm60(nwords),stat=status)
    if (status/=0) call outofmemory('realrecip3d3dVg','ktrm60')
    allocate(ktrm62(nwords),stat=status)
    if (status/=0) call outofmemory('realrecip3d3dVg','ktrm62')
    allocate(ktrm620(nwords),stat=status)
    if (status/=0) call outofmemory('realrecip3d3dVg','ktrm620')
  else
    allocate(ktrm6(1),stat=status)
    if (status/=0) call outofmemory('realrecip3d3dVg','ktrm6')
    allocate(ktrm60(1),stat=status)
    if (status/=0) call outofmemory('realrecip3d3dVg','ktrm60')
    allocate(ktrm62(1),stat=status)
    if (status/=0) call outofmemory('realrecip3d3dVg','ktrm62')
    allocate(ktrm620(1),stat=status)
    if (status/=0) call outofmemory('realrecip3d3dVg','ktrm620')
  endif
!
!  Calculate terms
!
  if (ndim.gt.1) then
    call setktrmd3(nkp,nkvec0,ktrm6,ktrm0,ktrm20,ktrm60,ktrm62,ktrm620,.false.)
  endif
!***************************************************************
!  Atomistic and real space electrostatic second derivatives   *
!***************************************************************
!
!  Outer loop over sites
!
  do i = 1,numat
!
!  Inner loop over second site
!
    xal = xclat(i)
    yal = yclat(i)
    zal = zclat(i)
    nati = nat(i)
    ntypi = nftype(i)
    lcorei = (nati.le.maxele)
    qli = qf(i)
    oci = occuf(i)
    ioc = iocptr(i)
    if (lcorei) rmassi = rfmass(3*(ioc-1)+1)
    if (lbsmat(nsft+nrelf2a(i))) then
      radi = radf(i)
    else
      radi = 0.0_dp
    endif
    indi = 3*(ioc-1)
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
    endif
!
!  Start of second atom loop
!
    do j = 1,i
      natj = nat(j)
      ntypj = nftype(j)
      lcorej = (natj.le.maxele)
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
      xcrd = xclat(j) - xal
      ycrd = yclat(j) - yal
      zcrd = zclat(j) - zal
      qlj = qf(j)
      ocj = occuf(j)
      joc = iocptr(j)
      if (lcorej) rmassj = rfmass(3*(joc-1)+1)
      if (lbsmat(nsft+nrelf2a(j))) then
        radj = radf(j)
      else
        radj = 0.0_dp
      endif
      radsum = radi + radj
      ofct = oci*ocj
      indj = 3*(joc - 1)
      jx = indj + 1
      jy = indj + 2
      jz = indj + 3
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
      else
        lmolok = .false.
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
              elseif (nptype(n).eq.2) then
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
      rp = sqrt(cut2)
!  
!  If the Coulomb term and the potentials don't need molecule info for
!  for this pair of atoms, then turn off molecule checking to save time
!
      if (.not.lneedmol) lmolok = .false.
!
      nmolonly = 0
      nor = 0
!
!  Zero third derivative arrays
!
      d3s(1:3,1:3,1:3) = 0.0_dp
!
!  If no valid potentials and charge product is zero
!  then no need to search for distances
!
      if (npots.eq.0.and.abs(factor).lt.1.0d-8) then
        if (lmany) then
          goto 1100
        else
          goto 1110
        endif
      endif
!*******************************************************
!  Loop over unit cells to find interatomic distances  *
!*******************************************************
!
!  Set up bonding
!
      if (lmolok) call getbonds(i,j)
!
      if (lra) then
!
!  Right angled cell 
!
        if (l111) then
          xcd = xcrd - (iimn+1)*r1x
          do ii = -iimn,iimx
            xcd = xcd + r1x
            xcd2 = xcd*xcd
            ycd = ycrd - (jjmn+1)*r2y
            do jj = - jjmn,jjmx
              ycd = ycd + r2y
              ycd2 = ycd*ycd
              zcd = zcrd - (kkmn+1)*r3z
              do kk = - kkmn,kkmx
                zcd = zcd + r3z
                r2 = xcd2 + ycd2 + zcd*zcd
                if (nor.ge.maxdis) then
                  maxdis = nor + 100
                  call changemaxdis
                endif
!
!  Molecule - check index
!
                if (lmolok.and.(r2.gt.cut2s.or..not.lcspair)) then
                  call bonded3c(lbonded(nor+1),l2bonds(nor+1),l3bonds(nor+1),nbotype(nor+1),nbotype2(nor+1),ii,jj,kk)
                  lsamemol = (lbonded(nor+1).or.l2bonds(nor+1).or.l3bonds(nor+1))
                  if (.not.lsamemol) then
                    call samemol(lsamemol,nmi,ii,jj,kk,ixj,iyj,izj)
                  endif
                  lptrmol(nor+1) = lsamemol
                  if (lsamemol) then
                    if (r2.gt.cut2q) nmolonly = nor + 1
                  else
                    lbonded(nor+1)  = .false.
                    l2bonds(nor+1)  = .false.
                    l3bonds(nor+1)  = .false.
                    nbotype(nor+1)  = 0
                    nbotype2(nor+1) = 0
                  endif
                else
                  lptrmol(nor+1)  = .false.
                  lbonded(nor+1)  = .false.
                  l2bonds(nor+1)  = .false.
                  l3bonds(nor+1)  = .false.
                  nbotype(nor+1)  = 0
                  nbotype2(nor+1) = 0
                endif
                if (r2.gt.small.and.(r2.le.cut2.or.lptrmol(nor+1))) then
!
!  Store vector - use in all terms
!
                  nor = nor + 1
                  dist(nor) = r2
                  xtmp(nor) = xcd
                  ytmp(nor) = ycd
                  ztmp(nor) = zcd
                endif
              enddo
            enddo
          enddo
        else
!
!  Non-l111
!
          max1u = (rp-xcrd)*ra + naddx
          max1l = (rp+xcrd)*ra + naddx
          max1l1 = max1l + 1
          xcd = xcrd - max1l1*r1x
          do ii = -max1l,max1u
            xcd = xcd + r1x
            if (abs(xcd).lt.rp) then
              xcd2 = xcd*xcd
              rpres2 = rp*rp - xcd2
              rpres = sqrt(rpres2)
              max2u = (rpres-ycrd)*rb + naddy
              max2l = (rpres+ycrd)*rb + naddy
              max2l1 = max2l + 1
              ycd = ycrd - max2l1*r2y
              do jj = -max2l,max2u
                ycd = ycd + r2y
                ycd2 = ycd*ycd
                rpres3 = rpres2 - ycd2
                if (rpres3.gt.0.0_dp) then
                  rpres3 = sqrt(rpres3)
                  max3u = (rpres3-zcrd)*rc + naddz
                  max3l = (rpres3+zcrd)*rc + naddz
                  max3l1 = max3l + 1
                  zcd = zcrd - max3l1*r3z
                  do kk = - max3l,max3u
                    zcd = zcd + r3z
                    r2 = xcd2 + ycd2 + zcd*zcd
                    if (nor.ge.maxdis) then
                      maxdis = nor + 100
                      call changemaxdis
                    endif
!
!  Molecule - check index
!
                    if (lmolok.and.(r2.gt.cut2s.or..not.lcspair)) then
                      call bonded3c(lbonded(nor+1),l2bonds(nor+1),l3bonds(nor+1),nbotype(nor+1),nbotype2(nor+1),ii,jj,kk)
                      lsamemol = (lbonded(nor+1).or.l2bonds(nor+1).or.l3bonds(nor+1))
                      if (.not.lsamemol) then
                        call samemol(lsamemol,nmi,ii,jj,kk,ixj,iyj,izj)
                      endif
                      lptrmol(nor+1) = lsamemol
                      if (lsamemol) then
                        if (r2.gt.cut2q) nmolonly = nor + 1
                      else
                        lbonded(nor+1) = .false.
                        l2bonds(nor+1) = .false.
                        l3bonds(nor+1) = .false.
                      endif
                    else
                      lptrmol(nor+1) = .false.
                      lbonded(nor+1) = .false.
                      l2bonds(nor+1) = .false.
                      l3bonds(nor+1) = .false.
                    endif
                    if (r2.gt.small.and.(r2.le.cut2.or.lptrmol(nor+1))) then
!
!  Store vector - use in all terms
!
                      nor = nor + 1
                      dist(nor) = r2
                      xtmp(nor) = xcd
                      ytmp(nor) = ycd
                      ztmp(nor) = zcd
                    endif
                  enddo
                endif
              enddo
            endif
          enddo
        endif
      else
!
!  General cell
!
        if (l111) then
          xcdi = xcrd - (iimn + 1)*r1x
          ycdi = ycrd - (iimn + 1)*r1y
          zcdi = zcrd - (iimn + 1)*r1z
          do ii = - iimn,iimx
            xcdi = xcdi + r1x
            ycdi = ycdi + r1y
            zcdi = zcdi + r1z
            xcdj = xcdi - (jjmn + 1)*r2x
            ycdj = ycdi - (jjmn + 1)*r2y
            zcdj = zcdi - (jjmn + 1)*r2z
            do jj = - jjmn,jjmx
              xcdj = xcdj + r2x
              ycdj = ycdj + r2y
              zcdj = zcdj + r2z
              xcd = xcdj - (kkmn + 1)*r3x
              ycd = ycdj - (kkmn + 1)*r3y
              zcd = zcdj - (kkmn + 1)*r3z
              do kk = - kkmn,kkmx
                xcd = xcd + r3x
                ycd = ycd + r3y
                zcd = zcd + r3z
                r2 = xcd*xcd + ycd*ycd + zcd*zcd
                if (nor.ge.maxdis) then
                  maxdis = nor + 100
                  call changemaxdis
                endif
!
!  Molecule - check index
!
                if (lmolok.and.(r2.gt.cut2s.or..not.lcspair)) then
                  call bonded3c(lbonded(nor+1),l2bonds(nor+1),l3bonds(nor+1),nbotype(nor+1),nbotype2(nor+1),ii,jj,kk)
                  lsamemol = (lbonded(nor+1).or.l2bonds(nor+1).or.l3bonds(nor+1))
                  if (.not.lsamemol) then
                    call samemol(lsamemol,nmi,ii,jj,kk,ixj,iyj,izj)
                  endif
                  lptrmol(nor+1) = lsamemol
                  if (lsamemol) then
                    if (r2.gt.cut2q) nmolonly = nor + 1
                  else
                    lbonded(nor+1)  = .false.
                    l2bonds(nor+1)  = .false.
                    l3bonds(nor+1)  = .false.
                    nbotype(nor+1)  = 0
                    nbotype2(nor+1) = 0
                  endif
                else
                  lptrmol(nor+1)  = .false.
                  lbonded(nor+1)  = .false.
                  l2bonds(nor+1)  = .false.
                  l3bonds(nor+1)  = .false.
                  nbotype(nor+1)  = 0
                  nbotype2(nor+1) = 0
                endif
                if (r2.gt.small.and.(r2.le.cut2.or.lptrmol(nor+1))) then
!
!  Store vector - use in all terms
!
                  nor = nor + 1
                  dist(nor) = r2
                  xtmp(nor) = xcd
                  ytmp(nor) = ycd
                  ztmp(nor) = zcd
                endif
              enddo
            enddo
          enddo
        else
!
!  Non - l111
!
          projr = xcrd*ru1x + ycrd*ru1y + zcrd*ru1z
          max1u = (rp - projr)*ra + naddx
          max1l = (rp + projr)*ra + naddx
          max1l1 = max1l + 1
          xcdi = xcrd - max1l1*r1x
          ycdi = ycrd - max1l1*r1y
          zcdi = zcrd - max1l1*r1z
          do ii = - max1l,max1u
            xcdi = xcdi + r1x
            ycdi = ycdi + r1y
            zcdi = zcdi + r1z
!
            projr = xcdi*ru2x + ycdi*ru2y + zcdi*ru2z
            max2u = (rp - projr)*rb + naddy
            max2l = (rp + projr)*rb + naddy
            max2l1 = max2l + 1
            xcdj = xcdi - max2l1*r2x
            ycdj = ycdi - max2l1*r2y
            zcdj = zcdi - max2l1*r2z
            do jj = - max2l,max2u
              xcdj = xcdj + r2x
              ycdj = ycdj + r2y
              zcdj = zcdj + r2z
!
              projr = xcdj*ru3x + ycdj*ru3y + zcdj*ru3z
              perpr = xcdj*xcdj + ycdj*ycdj + zcdj*zcdj
              perpr = perpr - projr*projr
              perpr = cut2 - perpr
              perpr = sqrt(abs(perpr))
              max3u = (perpr - projr)*rc + naddz
              max3l = (perpr + projr)*rc + naddz
              max3l1 = max3l + 1
!
              xcd = xcdj - max3l1*r3x
              ycd = ycdj - max3l1*r3y
              zcd = zcdj - max3l1*r3z
              do kk = - max3l,max3u
                xcd = xcd + r3x
                ycd = ycd + r3y
                zcd = zcd + r3z
                r2 = xcd*xcd + ycd*ycd + zcd*zcd
                if (nor.ge.maxdis) then
                  maxdis = nor + 100
                  call changemaxdis
                endif
!
!  Molecule - check index
!
                if (lmolok.and.(r2.gt.cut2s.or..not.lcspair)) then
                  call bonded3c(lbonded(nor+1),l2bonds(nor+1),l3bonds(nor+1),nbotype(nor+1),nbotype2(nor+1),ii,jj,kk)
                  lsamemol = (lbonded(nor+1).or.l2bonds(nor+1).or.l3bonds(nor+1))
                  if (.not.lsamemol) then
                    call samemol(lsamemol,nmi,ii,jj,kk,ixj,iyj,izj)
                  endif
                  lptrmol(nor+1) = lsamemol
                  if (lsamemol) then
                    if (r2.gt.cut2q) nmolonly = nor + 1
                  else
                    lbonded(nor+1)  = .false.
                    l2bonds(nor+1)  = .false.
                    l3bonds(nor+1)  = .false.
                    nbotype(nor+1)  = 0
                    nbotype2(nor+1) = 0
                  endif
                else
                  lptrmol(nor+1)  = .false.
                  lbonded(nor+1)  = .false.
                  l2bonds(nor+1)  = .false.
                  l3bonds(nor+1)  = .false.
                  nbotype(nor+1)  = 0
                  nbotype2(nor+1) = 0
                endif
                if (r2.gt.small.and.(r2.le.cut2.or.lptrmol(nor+1))) then
!
!  Store vector - use in all terms
!
                  nor = nor + 1
                  dist(nor) = r2
                  xtmp(nor) = xcd
                  ytmp(nor) = ycd
                  ztmp(nor) = zcd
                endif
              enddo
            enddo
          enddo
        endif
      endif
!
!  Sqrt distances
!
      do kk = 1,nor
        dist(kk) = sqrt(dist(kk))
      enddo
!
!  Twobody call
!
      call twobody(eatom,ereal,ec6,.true.,.true.,.true.,nor,1_i4,npots,npotl,cut2,cut2q, &
                   cut2s,nmolonly,factor,ofct,ospfct,radsum,sctrm1,sctrm2,qli,qlj, &
                   lcspair,lewaldtype,.false.,.false.,.false.)
!*****************************************************************
!  Calculate reciprocal space contribution to third derivatives  *
!*****************************************************************
      if (ndim.gt.1) then
        call reciptrmd3dVg(i,j,c6tot,ktrm6,ktrm62,d3s)
      endif
      if (index(keyword,'real').eq.0) then
!****************************
!  Loop over all distances  *
!****************************
        do kk = 1,nor
!
!  Generate products for derivatives
!
          xcrd = xtmp(kk)
          ycrd = ytmp(kk)
          zcrd = ztmp(kk)
          rdl1 = xcrd*xcrd
          rdl2 = ycrd*ycrd
          rdl3 = zcrd*zcrd
          rdl4 = ycrd*zcrd
          rdl5 = xcrd*zcrd
          rdl6 = xcrd*ycrd
          xd2 = xcrd*deriv2(kk)
          yd2 = ycrd*deriv2(kk)
          zd2 = zcrd*deriv2(kk)
          xd3 = xcrd*deriv3(kk)
          yd3 = ycrd*deriv3(kk)
          zd3 = zcrd*deriv3(kk)
!
!  Calculate third derivative matrix - first term
!
          d3k(1,1,1) = rdl1*xd3
          d3k(2,1,1) = rdl6*xd3
          d3k(3,1,1) = rdl5*xd3
          d3k(2,2,1) = rdl2*xd3
          d3k(3,2,1) = rdl4*xd3
          d3k(3,3,1) = rdl3*xd3
          d3k(2,2,2) = rdl2*yd3
          d3k(3,2,2) = rdl4*yd3
          d3k(3,3,2) = rdl3*yd3
          d3k(3,3,3) = rdl3*zd3
!
!  Add second term to third derivative matrix - now that the first term has been used for strain terms
!
          d3k(1,1,1) = d3k(1,1,1) + 3.0_dp*xd2
          d3k(2,1,1) = d3k(2,1,1) + yd2
          d3k(3,1,1) = d3k(3,1,1) + zd2
          d3k(2,2,1) = d3k(2,2,1) + xd2
          d3k(3,3,1) = d3k(3,3,1) + xd2
          d3k(2,2,2) = d3k(2,2,2) + 3.0_dp*yd2
          d3k(3,2,2) = d3k(3,2,2) + zd2
          d3k(3,3,2) = d3k(3,3,2) + yd2
          d3k(3,3,3) = d3k(3,3,3) + 3.0_dp*zd2
!
!  Strain derivatives
!
          d3ks(1,1,1) = d3k(1,1,1)*xcrd
          d3ks(2,1,1) = d3k(2,1,1)*xcrd
          d3ks(3,1,1) = d3k(3,1,1)*xcrd
          d3ks(2,2,1) = d3k(2,2,1)*xcrd
          d3ks(3,2,1) = d3k(3,2,1)*xcrd
          d3ks(3,3,1) = d3k(3,3,1)*xcrd
          d3ks(1,1,2) = d3k(2,1,1)*ycrd
          d3ks(2,1,2) = d3k(2,2,1)*ycrd
          d3ks(3,1,2) = d3k(3,2,1)*ycrd
          d3ks(2,2,2) = d3k(2,2,2)*ycrd
          d3ks(3,2,2) = d3k(3,2,2)*ycrd
          d3ks(3,3,2) = d3k(3,3,2)*ycrd
          d3ks(1,1,3) = d3k(3,1,1)*zcrd
          d3ks(2,1,3) = d3k(3,2,1)*zcrd
          d3ks(3,1,3) = d3k(3,3,1)*zcrd
          d3ks(2,2,3) = d3k(3,2,2)*zcrd
          d3ks(3,2,3) = d3k(3,3,2)*zcrd
          d3ks(3,3,3) = d3k(3,3,3)*zcrd
!
!  Calculate real component of third derivative matrix summed with unphased component for diagonal blocks
!
          do ii = 1,3
            d3s(1,1,ii) = d3s(1,1,ii) + d3ks(1,1,ii)
            d3s(2,1,ii) = d3s(2,1,ii) + d3ks(2,1,ii)
            d3s(3,1,ii) = d3s(3,1,ii) + d3ks(3,1,ii)
            d3s(2,2,ii) = d3s(2,2,ii) + d3ks(2,2,ii)
            d3s(3,2,ii) = d3s(3,2,ii) + d3ks(3,2,ii)
            d3s(3,3,ii) = d3s(3,3,ii) + d3ks(3,3,ii)
          enddo
!****************************
!  End loop over distances  *
!****************************
        enddo
      endif
1100  continue
      xcrd = xclat(j) - xal
      ycrd = yclat(j) - yal
      zcrd = zclat(j) - zal
!****************************
!  1-D electrostatic terms  *
!****************************
      if (ndim.eq.1) then
        call real1d3dVg(i,j,nati,natj,d3s,xcrd,ycrd,zcrd)
      endif
!***************************
!  Symmetrise d3 matrices  *
!***************************
      do ii = 1,3
        d3s(1,2,ii) = d3s(2,1,ii)
        d3s(1,3,ii) = d3s(3,1,ii)
        d3s(2,3,ii) = d3s(3,2,ii)
      enddo
!****************************
!  Three-body contribution  *
!****************************
      nmanyk = 0
      nforkl = 0
      if (lthb) then
        call three3d3dVg(i,j,nati,ntypi,natj,ntypj,d3s, &
          xal,yal,zal,xcrd,ycrd,zcrd,nmanyk)
      endif
!***************************
!  Four-body contribution  *
!***************************
      if (lfor) then
        call four3d3dVg(i,j,d3s,nmanyk,nforkl)
      endif
!******************************************************************************
!  Project contributions on to phonon modes to compute Grueneisen parameters  *
!******************************************************************************
      t1 = g_cpu_time()
      call projd3dVg(mcv,d3s,i,ix,iy,iz,j,jx,jy,jz,lcorei,lcorej, &
                     rmassi,rmassj,derv2,eigr,maxd2l,grueneisen,d33s, &
                     lmany,nmanyk,d34s,nforkl)
      t2 = g_cpu_time()
      tproj = tproj + t2 - t1
1110  continue
    enddo
  enddo
!
  if (.not.lnadd) nadd = 0
!************************
!  Free dynamic memory  *
!************************
  deallocate(ktrm620,stat=status)
  if (status/=0) call deallocate_error('realrecip3d3dVg','ktrm620')
  deallocate(ktrm62,stat=status)
  if (status/=0) call deallocate_error('realrecip3d3dVg','ktrm62')
  deallocate(ktrm60,stat=status)
  if (status/=0) call deallocate_error('realrecip3d3dVg','ktrm60')
  deallocate(ktrm6,stat=status)
  if (status/=0) call deallocate_error('realrecip3d3dVg','ktrm6')
  deallocate(ktrm20,stat=status)
  if (status/=0) call deallocate_error('realrecip3d3dVg','ktrm20')
  deallocate(ktrm0,stat=status)
  if (status/=0) call deallocate_error('realrecip3d3dVg','ktrm0')
!
!  Free local memory
!
  deallocate(npotl,stat=status)
  if (status/=0) call deallocate_error('realrecip3d3dVg','npotl')
!
!  Add on CPU time less time spent projecting derivatives
!
  time2 = g_cpu_time()
  tderv3 = tderv3 + time2 - time1 - tproj + tproj0
#ifdef TRACE
  call trace_out('realrecip3d3dVg')
#endif
!
  return
  end
