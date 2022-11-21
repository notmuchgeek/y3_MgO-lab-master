  subroutine fourimpsd2(eimp,lgrad1,lgrad2)
!
!  Subroutine for four-body improper torsion energy. 
!  Symmetry adapted version
!
!   7/13 Created from foursd2
!  12/14 Corrected so that ll vector is added to 42 and not 43
!   3/15 Error in one term in internal second derivatives fixed
!   1/17 ioptptr replaced by noptatrptr
!   2/18 Trace added
!   9/18 Strain module introduced
!   9/18 Call to fourstrterms replaced with more general realstrterms
!   4/19 Call to lmatchpair modified to include wildcard argument
!   4/19 Missing extra terms for strain derivatives added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!   4/20 Rigid molecule modifications added
!   5/20 Centre of mass handling added
!   7/20 lneedmol set to true for rigid molecule case
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
  use g_constants
  use control,        only : lrigid
  use current
  use derivatives
  use four
  use iochannels,     only : ioout
  use m_strain,       only : realstrterms
  use mdlogic
  use molecule
  use optimisation
  use parallel,       only : ioproc
  use symmetry
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  real(dp), intent(inout)                      :: eimp
  logical,  intent(in)                         :: lgrad1
  logical,  intent(in)                         :: lgrad2
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ii
  integer(i4)                                  :: iloop
  integer(i4)                                  :: isgn
  integer(i4)                                  :: indmi
  integer(i4)                                  :: indmj
  integer(i4)                                  :: indmk
  integer(i4)                                  :: indml
  integer(i4), dimension(:), allocatable       :: ioptaptr
  integer(i4)                                  :: ix
  integer(i4)                                  :: iy
  integer(i4)                                  :: iz
  integer(i4)                                  :: ixf
  integer(i4)                                  :: ixi
  integer(i4)                                  :: iyi
  integer(i4)                                  :: izi
  integer(i4)                                  :: ixj
  integer(i4)                                  :: iyj
  integer(i4)                                  :: izj
  integer(i4)                                  :: ixk
  integer(i4)                                  :: iyk
  integer(i4)                                  :: izk
  integer(i4)                                  :: ixl
  integer(i4)                                  :: iyl
  integer(i4)                                  :: izl
  integer(i4)                                  :: ixx
  integer(i4)                                  :: iyy
  integer(i4)                                  :: izz
  integer(i4)                                  :: j
  integer(i4)                                  :: jj
  integer(i4)                                  :: jx
  integer(i4)                                  :: jy
  integer(i4)                                  :: jz
  integer(i4)                                  :: jxf
  integer(i4)                                  :: jxx
  integer(i4)                                  :: jyy
  integer(i4)                                  :: jzz
  integer(i4)                                  :: k
  integer(i4)                                  :: kb(6,6)
  integer(i4)                                  :: ki
  integer(i4)                                  :: kj
  integer(i4)                                  :: kk
  integer(i4)                                  :: kl
  integer(i4)                                  :: kloop
  integer(i4)                                  :: ks
  integer(i4)                                  :: kt
  integer(i4)                                  :: kx
  integer(i4)                                  :: ky
  integer(i4)                                  :: kz
  integer(i4)                                  :: kxf
  integer(i4)                                  :: kxx
  integer(i4)                                  :: kyy
  integer(i4)                                  :: kzz
  integer(i4)                                  :: l
  integer(i4)                                  :: l4
  integer(i4)                                  :: li
  integer(i4)                                  :: lk
  integer(i4)                                  :: ll
  integer(i4)                                  :: lloop
  integer(i4)                                  :: lu
  integer(i4)                                  :: lx
  integer(i4)                                  :: lxf
  integer(i4)                                  :: ly
  integer(i4)                                  :: lz
  integer(i4)                                  :: n
  integer(i4)                                  :: n11
  integer(i4)                                  :: n11a
  integer(i4)                                  :: n21
  integer(i4)                                  :: n21a
  integer(i4)                                  :: n31
  integer(i4)                                  :: n31a
  integer(i4)                                  :: n41
  integer(i4)                                  :: n41a
  integer(i4)                                  :: n12
  integer(i4)                                  :: n12a
  integer(i4)                                  :: n22
  integer(i4)                                  :: n22a
  integer(i4)                                  :: n32
  integer(i4)                                  :: n32a
  integer(i4)                                  :: n42
  integer(i4)                                  :: n42a
  integer(i4)                                  :: n1x
  integer(i4)                                  :: n1xa
  integer(i4)                                  :: n1ya
  integer(i4)                                  :: n1za
  integer(i4)                                  :: n2x
  integer(i4)                                  :: n2xa
  integer(i4)                                  :: n2ya
  integer(i4)                                  :: n2za
  integer(i4)                                  :: n3x
  integer(i4)                                  :: n3xa
  integer(i4)                                  :: n3ya
  integer(i4)                                  :: n3za
  integer(i4)                                  :: n4x
  integer(i4)                                  :: n4xa
  integer(i4)                                  :: n4ya
  integer(i4)                                  :: n4za
  integer(i4)                                  :: n3vec(3,4)
  integer(i4), dimension(:), allocatable       :: natmiddle
  integer(i4), dimension(:), allocatable       :: ntypmiddle
  integer(i4)                                  :: nbtypeji
  integer(i4)                                  :: nbtypejk
  integer(i4)                                  :: nbtypekl
  integer(i4)                                  :: nbtypeji2
  integer(i4)                                  :: nbtypejk2
  integer(i4)                                  :: nbtypekl2
  integer(i4)                                  :: neqi
  integer(i4)                                  :: neqj
  integer(i4)                                  :: neqk
  integer(i4)                                  :: neql
  integer(i4)                                  :: nfortype
  integer(i4)                                  :: nforimp
  integer(i4)                                  :: ni
  integer(i4)                                  :: nil
  integer(i4)                                  :: niltor
  integer(i4)                                  :: nj
  integer(i4)                                  :: nk
  integer(i4)                                  :: nl 
  integer(i4)                                  :: nmi 
  integer(i4)                                  :: nmiddle
  integer(i4)                                  :: nmj 
  integer(i4)                                  :: nmk 
  integer(i4)                                  :: nml 
  integer(i4)                                  :: nn
  integer(i4)                                  :: npha
  integer(i4), dimension(:), allocatable       :: nptrnforimp
  integer(i4)                                  :: nreli
  integer(i4)                                  :: nrelj
  integer(i4)                                  :: nrelk
  integer(i4)                                  :: nrell
  integer(i4)                                  :: nt1
  integer(i4)                                  :: nt2
  integer(i4)                                  :: nt3
  integer(i4)                                  :: nt4
  integer(i4)                                  :: ntyp1
  integer(i4)                                  :: ntyp2
  integer(i4)                                  :: ntyp3
  integer(i4)                                  :: ntyp4
  integer(i4)                                  :: ntypi
  integer(i4)                                  :: ntypj
  integer(i4)                                  :: ntypk
  integer(i4)                                  :: ntypl
  integer(i4)                                  :: nuniquei
  integer(i4)                                  :: nuniquek
  integer(i4)                                  :: nuniquel
  integer(i4), dimension(:), allocatable       :: nuniqueiptr
  integer(i4), dimension(:), allocatable       :: nuniquekptr
  integer(i4), dimension(:), allocatable       :: nuniquelptr
  integer(i4)                                  :: status
  logical                                      :: l2bondsij
  logical                                      :: l2bondsjk
  logical                                      :: l2bondskl
  logical                                      :: lasu1
  logical                                      :: lasu2
  logical                                      :: lasu3
  logical                                      :: lasu4
  logical                                      :: lallbtyp
  logical                                      :: lanybtyp
  logical                                      :: lanyneedmol
  logical                                      :: lbondedij
  logical                                      :: lbondedjk
  logical                                      :: lbondedkl
  logical                                      :: lbtyp
  logical                                      :: libond
  logical                                      :: liok
  logical                                      :: limatch1
  logical                                      :: lintra_only
  logical                                      :: linter_only
  logical                                      :: ljmatch2
  logical                                      :: ljkmatch
  logical                                      :: lkbond
  logical                                      :: lkmatch3
  logical                                      :: llbond
  logical                                      :: lmatch
  logical                                      :: lmatchany
  logical                                      :: lmatchpair
  logical                                      :: lmolloc
  logical                                      :: lmolok
  logical                                      :: lmolokjk
  logical                                      :: lneedmol
  logical                                      :: lopi
  logical                                      :: lopj
  logical                                      :: lopk
  logical                                      :: lopl
  logical                                      :: lsamemolij
  logical                                      :: lsamemoljk
  logical                                      :: lsamemolkl
  logical                                      :: lsg1
  logical                                      :: lunique
  real(dp)                                     :: cut
  real(dp)                                     :: cutmax
  real(dp)                                     :: e1d(6)
  real(dp)                                     :: e2d(21)
  real(dp)                                     :: e3d(1)
  real(dp)                                     :: eterm
  real(dp)                                     :: fpoly(5)
  real(dp)                                     :: oci
  real(dp)                                     :: ocj
  real(dp)                                     :: ock
  real(dp)                                     :: ocl 
  real(dp)                                     :: ofct 
  real(dp)                                     :: phi0 
  real(dp)                                     :: phi0o
  real(dp)                                     :: r21
  real(dp)                                     :: r212
  real(dp)                                     :: r31
  real(dp)                                     :: r312
  real(dp)                                     :: r32
  real(dp)                                     :: r322
  real(dp)                                     :: r41
  real(dp)                                     :: r412
  real(dp)                                     :: r42
  real(dp)                                     :: r422
  real(dp)                                     :: r43
  real(dp)                                     :: r432
  real(dp)                                     :: rkforloc
  real(dp)                                     :: rko
  real(dp)                                     :: rn
  real(dp)                                     :: dr2ds(6,6)
  real(dp)                                     :: d2r2dx2(3,3,6)
  real(dp)                                     :: d2r2ds2(6,6,6)
  real(dp)                                     :: d2r2dsdx(6,3,6)
  real(dp)                                     :: rtfct
  real(dp)                                     :: t12
  real(dp)                                     :: t13
  real(dp)                                     :: t14
  real(dp)                                     :: t23
  real(dp)                                     :: t24
  real(dp)                                     :: t34
  real(dp)                                     :: temp(6)
  real(dp)                                     :: tr1
  real(dp)                                     :: tr2
  real(dp)                                     :: tr2max
  real(dp)                                     :: tr3
  real(dp)                                     :: tr4
  real(dp)                                     :: vec(3,3,4)
  real(dp)                                     :: x21
  real(dp)                                     :: y21
  real(dp)                                     :: z21
  real(dp)                                     :: x21t
  real(dp)                                     :: y21t
  real(dp)                                     :: z21t
  real(dp)                                     :: x31
  real(dp)                                     :: y31
  real(dp)                                     :: z31
  real(dp)                                     :: x41
  real(dp)                                     :: y41
  real(dp)                                     :: z41
  real(dp)                                     :: x32
  real(dp)                                     :: y32
  real(dp)                                     :: z32
  real(dp)                                     :: x32t
  real(dp)                                     :: y32t
  real(dp)                                     :: z32t
  real(dp)                                     :: x42
  real(dp)                                     :: y42
  real(dp)                                     :: z42
  real(dp)                                     :: x43
  real(dp)                                     :: y43
  real(dp)                                     :: z43
  real(dp)                                     :: x42t
  real(dp)                                     :: y42t
  real(dp)                                     :: z42t
  real(dp)                                     :: xc1t
  real(dp)                                     :: yc1t
  real(dp)                                     :: zc1t
  real(dp)                                     :: xc2
  real(dp)                                     :: yc2
  real(dp)                                     :: zc2
  real(dp)                                     :: xc3
  real(dp)                                     :: yc3
  real(dp)                                     :: zc3
  real(dp)                                     :: xc3t
  real(dp)                                     :: yc3t
  real(dp)                                     :: zc3t
  real(dp)                                     :: xc4t
  real(dp)                                     :: yc4t
  real(dp)                                     :: zc4t
  real(dp)                                     :: xcom(6)
  real(dp)                                     :: ycom(6)
  real(dp)                                     :: zcom(6)
  real(dp)                                     :: xcomj
  real(dp)                                     :: ycomj
  real(dp)                                     :: zcomj
  real(dp)                                     :: xv4(6)
  real(dp)                                     :: yv4(6)
  real(dp)                                     :: zv4(6)
!
  data kb/1,2,3,4,5,6,2,7,8,9,10,11,3,8,12,13,14,15,4,9,13,16,17,18,5,10,14,17,19,20,6,11,15,18,20,21/
  data n3vec/1,2,3,1,4,5,2,4,6,3,5,6/
#ifdef TRACE
  call trace_in('fourimpsd2')
#endif
!
  lsg1 = (lstr.and.lgrad1)
  lmolloc = (nmol.gt.0)
!
!  Allocate local memory
!
  allocate(natmiddle(2*nfor),stat=status)
  if (status/=0) call outofmemory('fourimpsd2','natmiddle')
  allocate(ntypmiddle(2*nfor),stat=status)
  if (status/=0) call outofmemory('fourimpsd2','ntypmiddle')
  allocate(nptrnforimp(nfor),stat=status)
  if (status/=0) call outofmemory('fourimpsd2','nptrnforimp')
  allocate(nuniqueiptr(numat),stat=status)
  if (status/=0) call outofmemory('fourimpsd2','nuniqueiptr')
  allocate(nuniquekptr(numat),stat=status)
  if (status/=0) call outofmemory('fourimpsd2','nuniquekptr')
  allocate(nuniquelptr(numat),stat=status)
  if (status/=0) call outofmemory('fourimpsd2','nuniquelptr')
  allocate(ioptaptr(nasym),stat=status)
  if (status/=0) call outofmemory('fourimpsd2','ioptaptr')
!
!  Initialise centre of mass arrays to zero
!
  xcom(1:6) = 0.0_dp
  ycom(1:6) = 0.0_dp
  zcom(1:6) = 0.0_dp
!
  lu = 0
  do i = 1,nasym
    if (lopf(i).or..not.lfreeze) then
      lu = lu + 1
      ioptaptr(i) = lu
    else
      ioptaptr(i) = 0
    endif
  enddo
!
!  Check how many four-body potentials are of out of plane type and how many aren't
!
  nforimp = 0
  do n = 1,nfor
    if (mmfexc(n).eq.2) then
      nforimp = nforimp + 1
      nptrnforimp(nforimp) = n
    endif
  enddo
!
!  Find out if any require molecule information and whether any potential is of bonded type
!
  lallbtyp = .false.
  lanybtyp = .false.
  lanyneedmol = .false.
  do n = 1,nfor
    if (mmfexc(n).eq.2) then
      lbtyp = (mmfexc(n).ge.1)
      lintra_only = (lfintra(n).and..not.lfinter(n))
      linter_only = (lfinter(n).and..not.lfintra(n))
      lneedmol = (lintra_only.or.linter_only.or.lbtyp.or.lrigid)
      if (lneedmol) lanyneedmol = .true.
      if (lbtyp) lanybtyp = .true.
      if (.not.lbtyp) lallbtyp = .false.
    endif
  enddo
!
!  Build a list of middle atom species types for potentials
!
  nmiddle = 0
  do n = 1,nfor
    if (mmfexc(n).eq.2) then
      if (.not.lmatchany(nfspec2(n),nfptyp2(n),nmiddle,natmiddle,ntypmiddle,.true.)) then
        nmiddle = nmiddle + 1
        natmiddle(nmiddle) = nfspec2(n)
        ntypmiddle(nmiddle) = nfptyp2(n)
      endif
      if (.not.lmatchany(nfspec3(n),nfptyp3(n),nmiddle,natmiddle,ntypmiddle,.true.)) then
        nmiddle = nmiddle + 1
        natmiddle(nmiddle) = nfspec3(n)
        ntypmiddle(nmiddle) = nfptyp3(n)
      endif
    endif
  enddo
!
!  Find maximum cutoff distance for ends and middle atoms
!
  cutmax = 0.0_dp
  tr2max = 0.0_dp
  do n = 1,nfor
    if (mmfexc(n).eq.2) then
      cut = max(for1(n)+for2(n),for2(n)+for3(n),for1(n)+for3(n),for4(n))
      cutmax = max(cut,cutmax)
      tr2 = for2(n)**2
      tr2max = max(tr2,tr2max)
    endif
  enddo
!
!  Opening banner for energy decomposition
!
  if (lPrintFour) then
    call mpbarrier
    if (ioproc) then
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''  IMP  : Atom No. 1  Atom No. 2  Atom No. 3  Atom No. 4    Torsion energy (eV)  '')')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
    endif
  endif
!********************************
!  Loop over middle site 2 / j  *
!********************************
  ljloop: do j = 1,numat
    nj = nat(j)
    ntypj = nftype(j)
    nrelj = nrelf2a(j)
    lasu2 = (nrela2f(nrelj).eq.j)
!
!  For energy only call only need terms where lasu2 is true
!
    if (.not.lgrad1.and..not.lasu2) cycle ljloop
!
    lopj = (lopf(nrelj).or..not.lfreeze)
    if (lopj) then
      if (lasu2) then
        jx = 3*ioptaptr(nrelj) - 2
        jy = jx + 1
        jz = jx + 2
      endif
      jxf = 3*noptatrptr(j) - 2
    endif
!
!  Check whether species may be valid
!
    if (.not.lmatchany(nj,ntypj,nmiddle,natmiddle,ntypmiddle,.true.)) cycle ljloop
!  
!  Set remaining properties for atom j
! 
    ocj = occuf(j)
    neqj = neqv(nrelj)
!
!  Set loop range for k
!
    if (lallbtyp) then
      lkbond = .true.
      if (nbonds(j).gt.0) then
        nuniquek = 1
        nuniquekptr(1) = nbonded(1,j)
        do kloop = 2,nbonds(j)
          lunique = .true.
          do lu = 1,nuniquek
            if (nbonded(kloop,j).eq.nuniquekptr(lu)) lunique = .false.
          enddo
          if (lunique) then
            nuniquek = nuniquek + 1
            nuniquekptr(nuniquek) = nbonded(kloop,j)
          endif
        enddo
        kloop = nuniquek
      else
        kloop = 0
      endif
    else
      lkbond = .false.
      kloop = numat
    endif
!
!  Skip if kloop is zero
!
    if (kloop.eq.0) cycle ljloop
!
!  Molecule handling
!
    if (lmolloc.and.lanyneedmol) then
      nmj = natmol(j)
      if (ndim.gt.0) then
        indmj = nmolind(j)
        call mindtoijk(indmj,ixj,iyj,izj)
      endif
!
!  Set COM coordinates
!
      if (lrigid.and.nmj.gt.0) then
        xcomj = molxyz(1,natinmol(j),nmj)
        ycomj = molxyz(2,natinmol(j),nmj)
        zcomj = molxyz(3,natinmol(j),nmj)
      else
        xcomj = 0.0_dp
        ycomj = 0.0_dp
        zcomj = 0.0_dp
      endif
    else
      xcomj = 0.0_dp
      ycomj = 0.0_dp
      zcomj = 0.0_dp
    endif
    xc2 = xclat(j)
    yc2 = yclat(j)
    zc2 = zclat(j)
!***************************************
!  Loop over second middle site 3 / k  *
!***************************************
    lkloop: do lk = 1,kloop
      if (lkbond) then
        k = nuniquekptr(lk)
      else
        k = lk
      endif
      nk = nat(k)
      ntypk = nftype(k)
!
!  Check whether species may be valid
!
      if (.not.lmatchany(nk,ntypk,nmiddle,natmiddle,ntypmiddle,.true.)) cycle lkloop
      if (.not.lmatchpair(nj,ntypj,nk,ntypk,nfor,nfspec2,nfptyp2,nfspec3,nfptyp3,.true.)) cycle lkloop
!
      lopk = (lopf(nrelf2a(k)).or..not.lfreeze)
      nrelk = nrelf2a(k)
      lasu3 = (nrela2f(nrelk).eq.k)
      if (lopk) then
        if (lasu3) then
          kx = 3*ioptaptr(nrelk) - 2
          ky = kx + 1
          kz = kx + 2
        endif
        kxf = 3*noptatrptr(k) - 2
      endif
!
      ock = occuf(k)
      neqk = neqv(nrelk)
!
      if (lmolloc.and.lanyneedmol) then
!
!  Molecule handling
!
        nmk = natmol(k)
        if (ndim.gt.0) then
          indmk = nmolind(k)
          call mindtoijk(indmk,ixk,iyk,izk)
          ixk = ixk - ixj
          iyk = iyk - iyj
          izk = izk - izj
        endif
        lmolokjk = (nmj.eq.nmk.and.nmj.gt.0)
!
!  Set COM coordinates
!
        if (lrigid.and.nmk.gt.0.and.lmolokjk) then
          xcom(4) = molxyz(1,natinmol(k),nmk) - xcomj
          ycom(4) = molxyz(2,natinmol(k),nmk) - ycomj
          zcom(4) = molxyz(3,natinmol(k),nmk) - zcomj
        else
          xcom(4) = - xcomj
          ycom(4) = - ycomj
          zcom(4) = - zcomj
        endif
      else
        lmolokjk = .false.
        xcom(4) = - xcomj
        ycom(4) = - ycomj
        zcom(4) = - zcomj
      endif
      xc3t = xclat(k)
      yc3t = yclat(k)
      zc3t = zclat(k)
      x32t = xc3t - xc2
      y32t = yc3t - yc2
      z32t = zc3t - zc2
!
!  Check r32 is OK
!  Loop over cell vectors
!
      do 120 jj = 1,iimax
        r322 = (xvec1cell(jj)+x32t)**2 + (yvec1cell(jj)+y32t)**2 + (zvec1cell(jj)+z32t)**2
        if (r322.lt.1d-12) goto 120
        if (k.eq.j.and.jj.eq.iimid) goto 120
!
!  Molecule checking
!
        lbondedjk = .false.
        if (lmolokjk) then
          if (ndim.eq.0) then
            if (lanybtyp) then
              call bonded(lbondedjk,l2bondsjk,nbtypejk,nbtypejk2,j,k,0_i4,0_i4,0_i4)
            endif
          else
            call lintoijk(jxx,jyy,jzz,jj,imaxl,jmaxl,kmaxl)
            if (lanybtyp) then
              call bonded(lbondedjk,l2bondsjk,nbtypejk,nbtypejk2,j,k,jxx,jyy,jzz)
              lsamemoljk = (lbondedjk.or.l2bondsjk)
            else
              lsamemoljk = .false.
            endif
            if (.not.lsamemoljk) then
              call samemol(lsamemoljk,nmj,jxx,jyy,jzz,ixk,iyk,izk)
            endif
          endif
        endif
!
!  Distance checking
!
        if (r322.gt.tr2max.and.(.not.lanybtyp.or..not.lbondedjk)) goto 120
!
        xc3 = xc3t + xvec1cell(jj)
        yc3 = yc3t + yvec1cell(jj)
        zc3 = zc3t + zvec1cell(jj)
        x32 = x32t + xvec1cell(jj)
        y32 = y32t + yvec1cell(jj)
        z32 = z32t + zvec1cell(jj)
!
!  Set counter for number of valid i/l end atom combinations
!
        niltor = 0
!***********************************
!  Loop over four-body potentials  *
!***********************************
        pots: do nn = 1,nforimp
          n = nptrnforimp(nn)
          nfortype = nforty(n)
          nt1 = nfspec1(n)
          nt2 = nfspec2(n)
          nt3 = nfspec3(n)
          nt4 = nfspec4(n)
          ntyp1 = nfptyp1(n)
          ntyp2 = nfptyp2(n)
          ntyp3 = nfptyp3(n)
          ntyp4 = nfptyp4(n)
          tr1 = for1(n)**2
          tr2 = for2(n)**2
          tr3 = for3(n)**2
          tr4 = for4(n)**2
          lbtyp = (mmfexc(n).ge.1)
          lintra_only = (lfintra(n).and..not.lfinter(n))
          linter_only = (lfinter(n).and..not.lfintra(n))
          lneedmol = (lintra_only.or.linter_only.or.lbtyp.or.lrigid)
!
!  Reset lmolok to initial state for j-k pair for each potential
!
          lmolok = lmolokjk
!************************************
!  Validate potential for j-k pair  *
!************************************
!
!  Check whether j and k are allowed for n
!
          ljmatch2 = lmatch(nj,ntypj,nt2,ntyp2,.true.)
          lkmatch3 = lmatch(nk,ntypk,nt3,ntyp3,.true.)
!
!  Check whether j-k or k-j orders are OK for 2-3
!
          ljkmatch = (ljmatch2.and.lkmatch3)
!
!  If no pair of matches can be found then cycle
!
          if (.not.ljkmatch) cycle pots
!
!  Distance checking for j-k
!
          if (lbtyp) then
            if (.not.lbondedjk) cycle pots
          else
            if (r322.gt.tr2) cycle pots
          endif
!       
!  Check for intra and but not in same molecule
!       
          if (lintra_only.and..not.lmolok) cycle pots
          if (lbtyp.and..not.lmolok) cycle pots
!                 
!  Molecule checking
!
          if (lmolok) then
            if (ndim.eq.0) then
              if (linter_only) cycle pots
              if (lbtyp) then
                if (.not.lbondedjk) cycle pots
!               
!  Check central bond type for correct order
!
                if (n4botype(1,n).gt.0) then
                  if (n4botype(1,n).ne.nbtypejk) cycle pots
                endif
                if (n4botype(2,n).ne.nbtypejk2) cycle pots
              endif
            else
              if (lbtyp) then
                if (.not.lbondedjk) cycle pots
!               
!  Check central bond type for correct order
!
                if (n4botype(1,n).gt.0) then
                  if (n4botype(1,n).ne.nbtypejk) cycle pots
                endif
                if (n4botype(2,n).ne.nbtypejk2) cycle pots
              endif
              if (lintra_only.and..not.lsamemoljk) cycle pots
              if (linter_only.and.lsamemoljk) cycle pots
            endif
          endif
!
!  Set loop range for i
!
          if (lbtyp) then
            libond = .true.
            if (nbonds(j).gt.0) then
              nuniquei = 1
              nuniqueiptr(1) = nbonded(1,j)
              do iloop = 2,nbonds(j)
                lunique = .true.
                do lu = 1,nuniquei
                  if (nbonded(iloop,j).eq.nuniqueiptr(lu)) lunique = .false.
                enddo
                if (lunique) then
                  nuniquei = nuniquei + 1
                  nuniqueiptr(nuniquei) = nbonded(iloop,j)
                endif
              enddo
              iloop = nuniquei
            else
              iloop = 0
            endif
          else
            libond = .false.
            iloop = numat
          endif
!
!  Skip if iloop is zero
!
          if (iloop.eq.0) cycle pots
!*****************************
!  Loop over end site 1 / i  *
!*****************************
          liloop: do li = 1,iloop
            if (libond) then
              i = nuniqueiptr(li)
            else
              i = li
            endif
            ni = nat(i)
            ntypi = nftype(i)
!
!  Check whether i matches either of types 1 and 4
!
            limatch1 = lmatch(ni,ntypi,nt1,ntyp1,.true.)
!
!  Is i allowed for type 1?
!
            liok = limatch1
            if (.not.liok) cycle liloop
!
            nreli = nrelf2a(i)
            lasu1 = (nrela2f(nreli).eq.i)
!
!  Set properties for atom i 
!
            oci = occuf(i)
            neqi = neqv(nreli)
            lopi = (lopf(nrelf2a(i)).or..not.lfreeze)
            if (lopi) then
              if (lasu1) then
                ix = 3*ioptaptr(nreli) - 2
              endif
              ixf = 3*noptatrptr(i) - 2
            endif
!
!  Molecule handling
!
            if (lmolloc.and.lneedmol) then
              nmi = natmol(i)
              if (ndim.gt.0) then
                indmi = nmolind(i)
                call mindtoijk(indmi,ixi,iyi,izi)
                ixi = ixi - ixj
                iyi = iyi - iyj
                izi = izi - izj
              endif
              lmolok = (nmj.eq.nmi.and.nmj.gt.0)
!
!  Set COM coordinates
!
              if (lrigid.and.nmi.gt.0.and.lmolok) then
                xcom(1) = xcomj - molxyz(1,natinmol(i),nmi)
                ycom(1) = ycomj - molxyz(2,natinmol(i),nmi)
                zcom(1) = zcomj - molxyz(3,natinmol(i),nmi)
              else
                xcom(1) = xcomj
                ycom(1) = ycomj
                zcom(1) = zcomj
              endif
              xcom(2) = xcom(4) + xcom(1)
              ycom(2) = ycom(4) + ycom(1)
              zcom(2) = zcom(4) + zcom(1)
            else
              lmolok = .false.
              xcom(1) = xcomj
              ycom(1) = ycomj
              zcom(1) = zcomj
              xcom(2) = xcom(4) + xcom(1)
              ycom(2) = ycom(4) + ycom(1)
              zcom(2) = zcom(4) + zcom(1)
            endif
!
!  Check for intra and but not in same molecule
!
            if (lintra_only.and..not.lmolok) cycle liloop
            if (lbtyp.and..not.lmolok) cycle liloop
!
            xc1t = xclat(i)
            yc1t = yclat(i)
            zc1t = zclat(i)
            x21t = xc2 - xc1t
            y21t = yc2 - yc1t
            z21t = zc2 - zc1t
!
!  Check r21 is OK
!  Loop over cell vectors
!
            iiloop: do ii = 1,iimax
              r212 = (-xvec1cell(ii) + x21t)**2 + (-yvec1cell(ii) + y21t)**2 + (-zvec1cell(ii) + z21t)**2
              if (r212.lt.1d-12) cycle iiloop
!
!  Prevent atoms i and k being the same atom
!
              if (k.eq.i.and.ii.eq.jj) cycle iiloop
!
!  Molecule checking
!
              lbondedij = .false.
              if (lmolok) then
                if (ndim.eq.0) then
                  if (linter_only) cycle iiloop
                  if (lbtyp) then
                    call bonded(lbondedij,l2bondsij,nbtypeji,nbtypeji2,j,i,0_i4,0_i4,0_i4)
                    if (.not.lbondedij) cycle iiloop
                  endif
                else
                  call lintoijk(ixx,iyy,izz,ii,imaxl,jmaxl,kmaxl)
                  if (lbtyp) then
                    call bonded(lbondedij,l2bondsij,nbtypeji,nbtypeji2,j,i,ixx,iyy,izz)
                    if (.not.lbondedij) cycle iiloop
                    lsamemolij = (lbondedij.or.l2bondsij)
                  else
                    lsamemolij = .false.
                  endif
                  if (.not.lsamemolij) then
                    call samemol(lsamemolij,nmj,ixx,iyy,izz,ixi,iyi,izi)
                  endif
                  if (lintra_only.and..not.lsamemolij) cycle iiloop
                  if (linter_only.and.lsamemolij) cycle iiloop
                endif
              endif
!
!  Distance checking
!
              if (lbtyp) then
                if (.not.lbondedij) cycle iiloop
              else
                if (r212.gt.tr1) cycle iiloop
              endif
!
              x21 = x21t - xvec1cell(ii)
              y21 = y21t - yvec1cell(ii)
              z21 = z21t - zvec1cell(ii)
!
!  Check r31 is OK
!
              x31 = x32 + x21
              y31 = y32 + y21
              z31 = z32 + z21
              r312 = x31*x31 + y31*y31 + z31*z31
              if (r312.lt.1.0d-12) cycle iiloop
!
!  Set loop range for l
!
              if (lbtyp) then
                llbond = .true.
                if (nbonds(j).gt.0) then
                  nuniquel = 1
                  nuniquelptr(1) = nbonded(1,j)
                  do lloop = 2,nbonds(j)
                    lunique = .true.
                    do lu = 1,nuniquel
                      if (nbonded(lloop,j).eq.nuniquelptr(lu)) lunique = .false.
                    enddo
                    if (lunique) then
                      nuniquel = nuniquel + 1
                      nuniquelptr(nuniquel) = nbonded(lloop,j)
                    endif
                  enddo
                  lloop = nuniquel
                else
                  lloop = 0
                endif
              else
                llbond = .false.
                lloop = numat
              endif
!
!  Skip if lloop is zero
!
              if (lloop.eq.0) cycle iiloop
!**********************************
!  Loop over last end site 4 / l  *
!**********************************
              l4loop: do l4 = 1,lloop
                if (llbond) then
                  l = nuniquelptr(l4)
                else
                  l = l4
                endif
                nrell = nrelf2a(l)
                lasu4 = (nrela2f(nrell).eq.l)
!
!  If none of the atoms are in the asymmetric unit then cycle
!
                if (.not.lasu1.and..not.lasu2.and..not.lasu3.and..not.lasu4) cycle l4loop
!
                nl = nat(l)
                ntypl = nftype(l)
!
!  Check l is allowed for n
!             
                if (.not.lmatch(nl,ntypl,nt4,ntyp4,.true.)) cycle l4loop
!
!  Set properties for atom l
!
                ocl = occuf(l)
                lopl = (lopf(nrelf2a(l)).or..not.lfreeze)
                if (lopl) then
                  if (lasu4) then
                    lx = 3*ioptaptr(nrell) - 2
                  endif
                  lxf = 3*noptatrptr(l) - 2
                endif
!
!  If lfreeze=.true. and no atoms have any variables
!  then skip this four body term
!
                if (.not.lopi.and..not.lopj.and..not.lopk.and..not.lopl) cycle l4loop
!
                if (lmolloc.and.lanyneedmol) then
!
!  Molecule handling
!
                  nml = natmol(l)
                  if (ndim.gt.0) then
                    indml = nmolind(l)
                    call mindtoijk(indml,ixl,iyl,izl)
                    ixl = ixl - ixj
                    iyl = iyl - iyj
                    izl = izl - izj
                  endif
                  lmolok = (nmj.eq.nml.and.nmj.gt.0)
!
!  Set COM coordinates
!
                  if (lrigid.and.nml.gt.0.and.lmolok) then
                    xcom(5) = molxyz(1,natinmol(l),nml) - xcomj
                    ycom(5) = molxyz(2,natinmol(l),nml) - ycomj
                    zcom(5) = molxyz(3,natinmol(l),nml) - zcomj
                  else
                    xcom(5) = - xcomj
                    ycom(5) = - ycomj
                    zcom(5) = - zcomj
                  endif
                  xcom(3) = xcom(5) + xcom(1)
                  ycom(3) = ycom(5) + ycom(1)
                  zcom(3) = zcom(5) + zcom(1)
                  xcom(6) = xcom(5) - xcom(4)
                  ycom(6) = ycom(5) - ycom(4)
                  zcom(6) = zcom(5) - zcom(4)
                else
                  lmolok = .false.
                  xcom(5) = - xcomj
                  ycom(5) = - ycomj
                  zcom(5) = - zcomj
                  xcom(3) = xcom(5) + xcom(1)
                  ycom(3) = ycom(5) + ycom(1)
                  zcom(3) = zcom(5) + zcom(1)
                  xcom(6) = xcom(5) - xcom(4)
                  ycom(6) = ycom(5) - ycom(4)
                  zcom(6) = zcom(5) - zcom(4)
                endif
!
!  Check for intra and but not in same molecule
!
                if (lintra_only.and..not.lmolok) cycle l4loop
                if (lbtyp.and..not.lmolok) cycle l4loop
!
                xc4t = xclat(l)
                yc4t = yclat(l)
                zc4t = zclat(l)
                x42t = xc4t - xc2
                y42t = yc4t - yc2
                z42t = zc4t - zc2
!
!  Check r42 is OK
!  Loop over cell vectors
!
                llloop: do ll = 1,iimax
                  r422 = (xvec1cell(ll)+x42t)**2 + (yvec1cell(ll)+y42t)**2 + (zvec1cell(ll)+z42t)**2
                  if (r422.lt.1d-12) cycle llloop
!
!  Molecule checking
!
                  lbondedkl = .false.
                  if (lmolok) then
                    if (ndim.eq.0) then
                      if (linter_only) cycle llloop
                      if (lbtyp) then
                        call bonded(lbondedkl,l2bondskl,nbtypekl,nbtypekl2,j,l,0_i4,0_i4,0_i4)
                        if (.not.lbondedkl) cycle llloop
                      endif
                    else
                      call lintoijk(kxx,kyy,kzz,ll,imaxl,jmaxl,kmaxl)
                      if (lbtyp) then
                        call bonded(lbondedkl,l2bondskl,nbtypekl,nbtypekl2,j,l,kxx,kyy,kzz)
                        if (.not.lbondedkl) cycle llloop
                        lsamemolkl = (lbondedkl.or.l2bondskl)
                      else
                        lsamemolkl = .false.
                      endif
                      if (.not.lsamemolkl) then
                        call samemol(lsamemolkl,nmj,kxx,kyy,kzz,ixl,iyl,izl)
                      endif
                      if (lintra_only.and..not.lsamemolkl) cycle llloop
                      if (linter_only.and.lsamemolkl) cycle llloop
                    endif
                  endif
!
!  Distance checking
!
                  if (lbtyp) then
                    if (.not.lbondedkl) cycle llloop
                  else
                    if (r422.gt.tr3) cycle llloop
                  endif
!
                  x42 = x42t + xvec1cell(ll)
                  y42 = y42t + yvec1cell(ll)
                  z42 = z42t + zvec1cell(ll)
!
!  Check r41 is OK
!
                  x41 = x42 + x21
                  y41 = y42 + y21
                  z41 = z42 + z21
                  r412 = x41*x41 + y41*y41 + z41*z41
                  if (r412.lt.1.0d-12) cycle llloop
!
!  Check r43 is OK
!
                  x43 = x42 - x32
                  y43 = y42 - y32
                  z43 = z42 - z32
                  r432 = x43*x43 + y43*y43 + z43*z43
                  if (r432.lt.1.0d-12) cycle llloop
!*********************************
!  Valid four-body term located  *
!*********************************
!
!  Finish calculating distances
!
                  r21 = sqrt(r212)
                  r31 = sqrt(r312)
                  r32 = sqrt(r322)
                  r41 = sqrt(r412)
                  r42 = sqrt(r422)
                  r43 = sqrt(r432)
!
!  Store information into iltor arrays
!
                  niltor = niltor + 1
                  if (niltor.gt.maxiltor) then
                    maxiltor = niltor + 3
                    call changemaxiltor
                  endif
!
                  nfortor(niltor) = n
!
                  iltor(1,niltor) = i
                  iltor(2,niltor) = l
                  ilftor(1,niltor) = ixf
                  ilftor(2,niltor) = lxf
                  ilxtor(1,niltor) = ix
                  ilxtor(2,niltor) = lx
!
                  riltor(1,niltor) = r21
                  riltor(2,niltor) = r31
                  riltor(3,niltor) = r41
                  riltor(4,niltor) = r42
                  riltor(5,niltor) = r43
!
                  xiltor(1,niltor) = x21
                  yiltor(1,niltor) = y21
                  ziltor(1,niltor) = z21
                  xiltor(2,niltor) = x31
                  yiltor(2,niltor) = y31
                  ziltor(2,niltor) = z31
                  xiltor(3,niltor) = x41
                  yiltor(3,niltor) = y41
                  ziltor(3,niltor) = z41
                  xiltor(4,niltor) = x42
                  yiltor(4,niltor) = y42
                  ziltor(4,niltor) = z42
                  xiltor(5,niltor) = x43
                  yiltor(5,niltor) = y43
                  ziltor(5,niltor) = z43
!
                  if (lrigid) then
                    xcomtor(1,niltor) = xcom(1)
                    ycomtor(1,niltor) = ycom(1)
                    zcomtor(1,niltor) = zcom(1)
                    xcomtor(2,niltor) = xcom(2)
                    ycomtor(2,niltor) = ycom(2)
                    zcomtor(2,niltor) = zcom(2)
                    xcomtor(3,niltor) = xcom(3)
                    ycomtor(3,niltor) = ycom(3)
                    zcomtor(3,niltor) = zcom(3)
                    xcomtor(4,niltor) = xcom(5)
                    ycomtor(4,niltor) = ycom(5)
                    zcomtor(4,niltor) = zcom(5)
                    xcomtor(5,niltor) = xcom(6)
                    ycomtor(5,niltor) = ycom(6)
                    zcomtor(5,niltor) = zcom(6)
                  endif
!
                  neqiltor(niltor) = neqi
                  oiltor(niltor) = oci*ocj*ock*ocl
                  lopiltor(1,niltor) = lopi
                  lopiltor(2,niltor) = lopl
!
!  End of inner loops over atoms and cell vectors
!
                enddo llloop
              enddo l4loop
            enddo iiloop
          enddo liloop
!
!  End loop over potentials
!
        enddo pots
!*******************************
!  Loop over i/l combinations  *
!*******************************
        do nil = 1,niltor
!
!  Return values to local variables
!
          n = nfortor(nil)
!
          i = iltor(1,nil)
          l = iltor(2,nil)
          nreli = nrelf2a(i)
          nrell = nrelf2a(l)
          neqi = neqv(nreli)
          neql = neqv(nrell)
!
          ixf = ilftor(1,nil)
          lxf = ilftor(2,nil)
!
          ix = ilxtor(1,nil)
          iy = ix + 1
          iz = ix + 2
          lx = ilxtor(2,nil)
          ly = lx + 1
          lz = lx + 2
!
          r21 = riltor(1,nil)
          r31 = riltor(2,nil)
          r41 = riltor(3,nil)
          r42 = riltor(4,nil)
          r43 = riltor(5,nil)
!
          x21 = xiltor(1,nil)
          y21 = yiltor(1,nil)
          z21 = ziltor(1,nil)
          x31 = xiltor(2,nil)
          y31 = yiltor(2,nil)
          z31 = ziltor(2,nil)
          x41 = xiltor(3,nil)
          y41 = yiltor(3,nil)
          z41 = ziltor(3,nil)
          x42 = xiltor(4,nil)
          y42 = yiltor(4,nil)
          z42 = ziltor(4,nil)
          x43 = xiltor(5,nil)
          y43 = yiltor(5,nil)
          z43 = ziltor(5,nil)
!
          if (lrigid) then
            xcom(1) = xcomtor(1,nil)
            ycom(1) = ycomtor(1,nil)
            zcom(1) = zcomtor(1,nil)
            xcom(2) = xcomtor(2,nil)
            ycom(2) = ycomtor(2,nil)
            zcom(2) = zcomtor(2,nil)
            xcom(3) = xcomtor(3,nil)
            ycom(3) = ycomtor(3,nil)
            zcom(3) = zcomtor(3,nil)
            xcom(5) = xcomtor(4,nil)
            ycom(5) = ycomtor(4,nil)
            zcom(5) = zcomtor(4,nil)
            xcom(6) = xcomtor(5,nil)
            ycom(6) = ycomtor(5,nil)
            zcom(6) = zcomtor(5,nil)
          endif
!
          ofct = oiltor(nil)
          lopi = lopiltor(1,nil)
          lopl = lopiltor(2,nil)
          lasu1 = (nrela2f(nreli).eq.i)
          lasu4 = (nrela2f(nrell).eq.l)
!
!  Set terms for potentials
!
          rkforloc = fork(n)
!
!  If this is Dreiding mode then divide force constant by number of torsions
!
          if (lfdreiding(n)) then
            rkforloc = rkforloc/dble(niltor)
          endif
          npha = 0
          nfortype = nforty(n)
          if (nfortype.eq.1) then
            npha = npfor(n)
            if (npha.gt.0) then
              isgn = 1
            else
              isgn = - 1
            endif
            npha = abs(npha)
            phi0 = forpoly(1,n)*degtorad
          elseif (nfortype.eq.4.or.nfortype.eq.6.or.nfortype.eq.7) then
            npha = npfor(n)
            if (npha.gt.0) then
              isgn = 1
            else
              isgn = - 1
            endif
            npha = abs(npha)
            if (nfortype.eq.6) then
              phi0 = forpoly(1,n)*degtorad
            else
              phi0 = forpoly(1,n)
            endif
            if (nfortype.eq.6.or.nfortype.eq.7) then
              fpoly(2:4) = forpoly(2:4,n)
            endif
          elseif (nfortype.eq.8.or.nfortype.eq.9) then
            npha = npfor(n)
            if (npha.gt.0) then
              isgn = 1
            else
              isgn = - 1
            endif
            npha = abs(npha)
            if (nfortype.eq.8) then
              phi0 = forpoly(1,n)*degtorad
            else
              phi0 = forpoly(1,n)
            endif
            fpoly(2) = forpoly(2,n)
            fpoly(3) = for1(n)
            fpoly(4) = for2(n)
            fpoly(5) = for3(n)
          elseif (nfortype.eq.2) then
            npha = npfor(n)
          elseif (nfortype.eq.5) then
            phi0 = forpoly(1,n)*degtorad
          elseif (nfortype.eq.10.or.nfortype.eq.17) then
            fpoly(1) = forpoly(1,n)*degtorad
            fpoly(2) = forpoly(2,n)*degtorad
          elseif (nfortype.eq.13) then
            npha = abs(npfor(n))
            phi0 = forpoly(1,n)*degtorad
          endif
          rn = dble(npha)
!
!  Scaling of terms
!
          rko = rkforloc*ofct
          phi0o = phi0
          if (nfortype.eq.2) then
            do kk = 1,npha
              fpoly(kk) = forpoly(kk,n)*ofct
            enddo
          elseif (nfortype.eq.4.or.nfortype.eq.7.or.nfortype.eq.9) then
            phi0o = phi0*ofct
          endif
!
!  Call subroutine to calculate energy and derivatives
!
          call fourbody(n,nfortype,r21,r31,r41,r32,r42,r43,eterm,e1d,e2d,e3d, &
                        rko,rn,phi0o,isgn,fpoly,lgrad1,lgrad2,.false.)
          rtfct = dble(neqj)
          if (lasu2) then
            eimp = eimp + eterm*rtfct
!
!  Output energy contribution
!
            if (lPrintFour) then
              write(ioout,'(4x,4i12,1x,f22.10)') i,j,k,l,eterm*rtfct
            endif
          endif
!**************************
!  Torsional derivatives  *
!**************************
          if (lgrad1) then
!
!  Set up strain products
!
            if (lsg1) then
              xv4(1) = x21
              yv4(1) = y21
              zv4(1) = z21
              xv4(2) = x31
              yv4(2) = y31
              zv4(2) = z31
              xv4(3) = x41
              yv4(3) = y41
              zv4(3) = z41
              xv4(4) = x32
              yv4(4) = y32
              zv4(4) = z32
              xv4(5) = x42
              yv4(5) = y42
              zv4(5) = z42
              xv4(6) = x43
              yv4(6) = y43
              zv4(6) = z43
              call realstrterms(ndim,6_i4,6_i4,xv4,yv4,zv4,xcom,ycom,zcom,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,lgrad2)
            endif
          endif
          if (lgrad2) then
            n1x = ixf
            n2x = jxf
            n3x = kxf
            n4x = lxf
            n1xa = ix
            n1ya = iy
            n1za = iz
            n2xa = jx
            n2ya = jy
            n2za = jz
            n3xa = kx
            n3ya = ky
            n3za = kz
            n4xa = lx
            n4ya = ly
            n4za = lz
          endif
!***********************
!  Strain derivatives  *
!***********************
          if (lsg1) then
            if (lasu2) then
!
!  First strain derivatives
!
              do kl = 1,nstrains
                ks = nstrptr(kl)
                rstrd(kl) = rstrd(kl) + e1d(1)*dr2ds(ks,1)*rtfct
                rstrd(kl) = rstrd(kl) + e1d(2)*dr2ds(ks,2)*rtfct
                rstrd(kl) = rstrd(kl) + e1d(3)*dr2ds(ks,3)*rtfct
                rstrd(kl) = rstrd(kl) + e1d(4)*dr2ds(ks,4)*rtfct
                rstrd(kl) = rstrd(kl) + e1d(5)*dr2ds(ks,5)*rtfct
                rstrd(kl) = rstrd(kl) + e1d(6)*dr2ds(ks,6)*rtfct
              enddo
            endif
!
!  Second strain derivatives
!
!  Strain only
!
            if (lgrad2.and.lstr) then
              do kk = 1,nstrains
                ks = nstrptr(kk)
                do kl = 1,6
                  temp(kl) = e2d(kb(1,kl))*dr2ds(ks,1)
                  temp(kl) = temp(kl) + e2d(kb(2,kl))*dr2ds(ks,2)
                  temp(kl) = temp(kl) + e2d(kb(3,kl))*dr2ds(ks,3)
                  temp(kl) = temp(kl) + e2d(kb(4,kl))*dr2ds(ks,4)
                  temp(kl) = temp(kl) + e2d(kb(5,kl))*dr2ds(ks,5)
                  temp(kl) = temp(kl) + e2d(kb(6,kl))*dr2ds(ks,6)
                enddo
                if (lasu2) then
!
!  Strain-strain second derivatives
!
                  do kl = 1,nstrains
                    kt = nstrptr(kl)
                    sderv2(kl,kk) = sderv2(kl,kk) + dr2ds(kt,1)*temp(1)*rtfct
                    sderv2(kl,kk) = sderv2(kl,kk) + dr2ds(kt,2)*temp(2)*rtfct
                    sderv2(kl,kk) = sderv2(kl,kk) + dr2ds(kt,3)*temp(3)*rtfct
                    sderv2(kl,kk) = sderv2(kl,kk) + dr2ds(kt,4)*temp(4)*rtfct
                    sderv2(kl,kk) = sderv2(kl,kk) + dr2ds(kt,5)*temp(5)*rtfct
                    sderv2(kl,kk) = sderv2(kl,kk) + dr2ds(kt,6)*temp(6)*rtfct
!
                    sderv2(kl,kk) = sderv2(kl,kk) + d2r2ds2(kt,ks,1)*e1d(1)*rtfct
                    sderv2(kl,kk) = sderv2(kl,kk) + d2r2ds2(kt,ks,2)*e1d(2)*rtfct
                    sderv2(kl,kk) = sderv2(kl,kk) + d2r2ds2(kt,ks,3)*e1d(3)*rtfct
                    sderv2(kl,kk) = sderv2(kl,kk) + d2r2ds2(kt,ks,4)*e1d(4)*rtfct
                    sderv2(kl,kk) = sderv2(kl,kk) + d2r2ds2(kt,ks,5)*e1d(5)*rtfct
                    sderv2(kl,kk) = sderv2(kl,kk) + d2r2ds2(kt,ks,6)*e1d(6)*rtfct
                  enddo
                endif
!  
!  Mixed derivatives
!                   
                if (lasu1) then
                  derv3(n1xa,kk) = derv3(n1xa,kk) - (x21*temp(1) + x31*temp(2) + x41*temp(3))*dble(neqi)
                  derv3(n1ya,kk) = derv3(n1ya,kk) - (y21*temp(1) + y31*temp(2) + y41*temp(3))*dble(neqi)
                  derv3(n1za,kk) = derv3(n1za,kk) - (z21*temp(1) + z31*temp(2) + z41*temp(3))*dble(neqi)
!
                  derv3(n1xa,kk) = derv3(n1xa,kk) - (e1d(1)*d2r2dsdx(ks,1,1) + &
                    e1d(2)*d2r2dsdx(ks,1,2) + e1d(3)*d2r2dsdx(ks,1,3))*dble(neqi)
                  derv3(n1ya,kk) = derv3(n1ya,kk) - (e1d(1)*d2r2dsdx(ks,2,1) + &
                    e1d(2)*d2r2dsdx(ks,2,2) + e1d(3)*d2r2dsdx(ks,2,3))*dble(neqi)
                  derv3(n1za,kk) = derv3(n1za,kk) - (e1d(1)*d2r2dsdx(ks,3,1) + &
                    e1d(2)*d2r2dsdx(ks,3,2) + e1d(3)*d2r2dsdx(ks,3,3))*dble(neqi)
                endif
                if (lasu2) then
                  derv3(n2xa,kk) = derv3(n2xa,kk) - (x32*temp(4) - x21*temp(1) + x42*temp(5))*dble(neqj)
                  derv3(n2ya,kk) = derv3(n2ya,kk) - (y32*temp(4) - y21*temp(1) + y42*temp(5))*dble(neqj)
                  derv3(n2za,kk) = derv3(n2za,kk) - (z32*temp(4) - z21*temp(1) + z42*temp(5))*dble(neqj)
!
                  derv3(n2xa,kk) = derv3(n2xa,kk) - (e1d(4)*d2r2dsdx(ks,1,4) - &
                    e1d(1)*d2r2dsdx(ks,1,1) + e1d(5)*d2r2dsdx(ks,1,5))*dble(neqj)
                  derv3(n2ya,kk) = derv3(n2ya,kk) - (e1d(4)*d2r2dsdx(ks,2,4) - &
                    e1d(1)*d2r2dsdx(ks,2,1) + e1d(5)*d2r2dsdx(ks,2,5))*dble(neqj)
                  derv3(n2za,kk) = derv3(n2za,kk) - (e1d(4)*d2r2dsdx(ks,3,4) - &
                    e1d(1)*d2r2dsdx(ks,3,1) + e1d(5)*d2r2dsdx(ks,3,5))*dble(neqj)
                endif
                if (lasu3) then
                  derv3(n3xa,kk) = derv3(n3xa,kk) + (x32*temp(4) - x43*temp(6) + x31*temp(2))*dble(neqk)
                  derv3(n3ya,kk) = derv3(n3ya,kk) + (y32*temp(4) - y43*temp(6) + y31*temp(2))*dble(neqk)
                  derv3(n3za,kk) = derv3(n3za,kk) + (z32*temp(4) - z43*temp(6) + z31*temp(2))*dble(neqk)
!
                  derv3(n3xa,kk) = derv3(n3xa,kk) + (e1d(4)*d2r2dsdx(ks,1,4) - &
                    e1d(6)*d2r2dsdx(ks,1,6) + e1d(2)*d2r2dsdx(ks,1,2))*dble(neqk)
                  derv3(n3ya,kk) = derv3(n3ya,kk) + (e1d(4)*d2r2dsdx(ks,2,4) - &
                    e1d(6)*d2r2dsdx(ks,2,6) + e1d(2)*d2r2dsdx(ks,2,2))*dble(neqk)
                  derv3(n3za,kk) = derv3(n3za,kk) + (e1d(4)*d2r2dsdx(ks,3,4) - &
                    e1d(6)*d2r2dsdx(ks,3,6) + e1d(2)*d2r2dsdx(ks,3,2))*dble(neqk)
                endif
                if (lasu4) then
                  derv3(n4xa,kk) = derv3(n4xa,kk) + (x43*temp(6) + x42*temp(5) + x41*temp(3))*dble(neql)
                  derv3(n4ya,kk) = derv3(n4ya,kk) + (y43*temp(6) + y42*temp(5) + y41*temp(3))*dble(neql)
                  derv3(n4za,kk) = derv3(n4za,kk) + (z43*temp(6) + z42*temp(5) + z41*temp(3))*dble(neql)
!
                  derv3(n4xa,kk) = derv3(n4xa,kk) + (e1d(6)*d2r2dsdx(ks,1,6) + &
                    e1d(5)*d2r2dsdx(ks,1,5) + e1d(3)*d2r2dsdx(ks,1,3))*dble(neql)
                  derv3(n4ya,kk) = derv3(n4ya,kk) + (e1d(6)*d2r2dsdx(ks,2,6) + &
                    e1d(5)*d2r2dsdx(ks,2,5) + e1d(3)*d2r2dsdx(ks,2,3))*dble(neql)
                  derv3(n4za,kk) = derv3(n4za,kk) + (e1d(6)*d2r2dsdx(ks,3,6) + &
                    e1d(5)*d2r2dsdx(ks,3,5) + e1d(3)*d2r2dsdx(ks,3,3))*dble(neql)
                endif
              enddo
            endif
          endif
!*************************
!  Internal derivatives  *
!*************************
          if (lgrad1) then
            if (lasu1) then
              xdrv(nreli) = xdrv(nreli) - (x21*e1d(1) + x31*e1d(2) + x41*e1d(3))*dble(neqi)
              ydrv(nreli) = ydrv(nreli) - (y21*e1d(1) + y31*e1d(2) + y41*e1d(3))*dble(neqi)
              zdrv(nreli) = zdrv(nreli) - (z21*e1d(1) + z31*e1d(2) + z41*e1d(3))*dble(neqi)
            endif
!
            if (lasu2) then
              xdrv(nrelj) = xdrv(nrelj) - (x32*e1d(4) - x21*e1d(1) + x42*e1d(5))*dble(neqj)
              ydrv(nrelj) = ydrv(nrelj) - (y32*e1d(4) - y21*e1d(1) + y42*e1d(5))*dble(neqj)
              zdrv(nrelj) = zdrv(nrelj) - (z32*e1d(4) - z21*e1d(1) + z42*e1d(5))*dble(neqj)
            endif
!
            if (lasu3) then
              xdrv(nrelk) = xdrv(nrelk) + (x32*e1d(4) - x43*e1d(6) + x31*e1d(2))*dble(neqk)
              ydrv(nrelk) = ydrv(nrelk) + (y32*e1d(4) - y43*e1d(6) + y31*e1d(2))*dble(neqk)
              zdrv(nrelk) = zdrv(nrelk) + (z32*e1d(4) - z43*e1d(6) + z31*e1d(2))*dble(neqk)
            endif
!
            if (lasu4) then
              xdrv(nrell) = xdrv(nrell) + (x43*e1d(6) + x42*e1d(5) + x41*e1d(3))*dble(neql)
              ydrv(nrell) = ydrv(nrell) + (y43*e1d(6) + y42*e1d(5) + y41*e1d(3))*dble(neql)
              zdrv(nrell) = zdrv(nrell) + (z43*e1d(6) + z42*e1d(5) + z41*e1d(3))*dble(neql)
            endif
          endif
          if (lgrad2) then
!
!  New vector array between atoms to handle sign 
!
!  Atom 1
!
            vec(1,1,1) = -x21
            vec(2,1,1) = -y21
            vec(3,1,1) = -z21
            vec(1,2,1) = -x31
            vec(2,2,1) = -y31
            vec(3,2,1) = -z31
            vec(1,3,1) = -x41
            vec(2,3,1) = -y41
            vec(3,3,1) = -z41
!
!  Atom 2
!
            vec(1,1,2) = x21
            vec(2,1,2) = y21
            vec(3,1,2) = z21
            vec(1,2,2) = -x32
            vec(2,2,2) = -y32
            vec(3,2,2) = -z32
            vec(1,3,2) = -x42
            vec(2,3,2) = -y42
            vec(3,3,2) = -z42
!
!  Atom 3
!
            vec(1,1,3) = x31
            vec(2,1,3) = y31
            vec(3,1,3) = z31
            vec(1,2,3) = x32
            vec(2,2,3) = y32
            vec(3,2,3) = z32
            vec(1,3,3) = -x43
            vec(2,3,3) = -y43
            vec(3,3,3) = -z43
!
!  Atom 4
!
            vec(1,1,4) = x41
            vec(2,1,4) = y41
            vec(3,1,4) = z41
            vec(1,2,4) = x42
            vec(2,2,4) = y42
            vec(3,2,4) = z42
            vec(1,3,4) = x43
            vec(2,3,4) = y43
            vec(3,3,4) = z43
!
!  Loop over first coordinate
!
            do kk = 1,3
              n11 = n1x - 1 + kk
              n21 = n2x - 1 + kk
              n31 = n3x - 1 + kk
              n41 = n4x - 1 + kk
              n11a = n1xa - 1 + kk
              n21a = n2xa - 1 + kk
              n31a = n3xa - 1 + kk
              n41a = n4xa - 1 + kk
!
!  First term
!
              if (lasu1.and.lopi) then
                if (lopj) then
                  derv2(n21,n11a) = derv2(n21,n11a) - e1d(1)*dble(neqi)
                else
                  derv2(n11,n11a) = derv2(n11,n11a) - e1d(1)*dble(neqi)
                endif
                if (lopk) then
                  derv2(n31,n11a) = derv2(n31,n11a) - e1d(2)*dble(neqi)
                else
                  derv2(n11,n11a) = derv2(n11,n11a) - e1d(2)*dble(neqi)
                endif
                if (lopl) then
                  derv2(n41,n11a) = derv2(n41,n11a) - e1d(3)*dble(neqi)
                else
                  derv2(n11,n11a) = derv2(n11,n11a) - e1d(3)*dble(neqi)
                endif
              endif
!
              if (lasu2.and.lopj) then
                if (lopi) then
                  derv2(n11,n21a) = derv2(n11,n21a) - e1d(1)*dble(neqj)
                else
                  derv2(n21,n21a) = derv2(n21,n21a) - e1d(1)*dble(neqj)
                endif
                if (lopk) then
                  derv2(n31,n21a) = derv2(n31,n21a) - e1d(4)*dble(neqj)
                else
                  derv2(n21,n21a) = derv2(n21,n21a) - e1d(4)*dble(neqj)
                endif
                if (lopl) then
                  derv2(n41,n21a) = derv2(n41,n21a) - e1d(5)*dble(neqj)
                else
                  derv2(n21,n21a) = derv2(n21,n21a) - e1d(5)*dble(neqj)
                endif
              endif
!
              if (lasu3.and.lopk) then
                if (lopi) then
                  derv2(n11,n31a) = derv2(n11,n31a) - e1d(2)*dble(neqk)
                else
                  derv2(n31,n31a) = derv2(n31,n31a) - e1d(2)*dble(neqk)
                endif
                if (lopj) then
                  derv2(n21,n31a) = derv2(n21,n31a) - e1d(4)*dble(neqk)
                else
                  derv2(n31,n31a) = derv2(n31,n31a) - e1d(4)*dble(neqk)
                endif
                if (lopl) then
                  derv2(n41,n31a) = derv2(n41,n31a) - e1d(6)*dble(neqk)
                else
                  derv2(n31,n31a) = derv2(n31,n31a) - e1d(6)*dble(neqk)
                endif
              endif
              if (lasu4.and.lopl) then
                if (lopi) then
                  derv2(n11,n41a) = derv2(n11,n41a) - e1d(3)*dble(neql)
                else
                  derv2(n41,n41a) = derv2(n41,n41a) - e1d(3)*dble(neql)
                endif
                if (lopj) then
                  derv2(n21,n41a) = derv2(n21,n41a) - e1d(5)*dble(neql)
                else
                  derv2(n41,n41a) = derv2(n41,n41a) - e1d(5)*dble(neql)
                endif
                if (lopk) then
                  derv2(n31,n41a) = derv2(n31,n41a) - e1d(6)*dble(neql)
                else
                  derv2(n41,n41a) = derv2(n41,n41a) - e1d(6)*dble(neql)
                endif
              endif
!
!  Loop over second coordinate
!
              do kl = 1,3
                n12 = n1x - 1 + kl
                n22 = n2x - 1 + kl
                n32 = n3x - 1 + kl
                n42 = n4x - 1 + kl
                n12a = n1xa - 1 + kl
                n22a = n2xa - 1 + kl
                n32a = n3xa - 1 + kl
                n42a = n4xa - 1 + kl
!
!  Sum over vectors atom-atom second derivatives
!
                t12 = 0.0_dp
                t13 = 0.0_dp
                t14 = 0.0_dp
                t23 = 0.0_dp
                t24 = 0.0_dp
                t34 = 0.0_dp
                do ki = 1,3
                  do kj = 1,3
                    t12 = t12 + vec(kk,ki,1)*vec(kl,kj,2)*e2d(kb(n3vec(ki,1),n3vec(kj,2)))
                    t13 = t13 + vec(kk,ki,1)*vec(kl,kj,3)*e2d(kb(n3vec(ki,1),n3vec(kj,3)))
                    t14 = t14 + vec(kk,ki,1)*vec(kl,kj,4)*e2d(kb(n3vec(ki,1),n3vec(kj,4)))
                    t23 = t23 + vec(kk,ki,2)*vec(kl,kj,3)*e2d(kb(n3vec(ki,2),n3vec(kj,3)))
                    t24 = t24 + vec(kk,ki,2)*vec(kl,kj,4)*e2d(kb(n3vec(ki,2),n3vec(kj,4)))
                    t34 = t34 + vec(kk,ki,3)*vec(kl,kj,4)*e2d(kb(n3vec(ki,3),n3vec(kj,4)))
                  enddo
                enddo
!
!  Add terms to total second derivative matrix
!
                if (lasu1.and.lopi) then
                  if (lopj) then
                    derv2(n22,n11a) = derv2(n22,n11a) + t12*dble(neqi)
                  else
                    derv2(n12,n11a) = derv2(n12,n11a) + t12*dble(neqi)
                  endif
                  if (lopk) then
                    derv2(n32,n11a) = derv2(n32,n11a) + t13*dble(neqi)
                  else
                    derv2(n12,n11a) = derv2(n12,n11a) + t13*dble(neqi)
                  endif
                  if (lopl) then
                    derv2(n42,n11a) = derv2(n42,n11a) + t14*dble(neqi)
                  else
                    derv2(n12,n11a) = derv2(n12,n11a) + t14*dble(neqi)
                  endif
                endif
                if (lasu2.and.lopj) then
                  if (lopi) then
                    derv2(n11,n22a) = derv2(n11,n22a) + t12*dble(neqj)
                  else
                    derv2(n21,n22a) = derv2(n21,n22a) + t12*dble(neqj)
                  endif
                  if (lopk) then
                    derv2(n32,n21a) = derv2(n32,n21a) + t23*dble(neqj)
                  else
                    derv2(n22,n21a) = derv2(n22,n21a) + t23*dble(neqj)
                  endif
                  if (lopl) then
                    derv2(n42,n21a) = derv2(n42,n21a) + t24*dble(neqj)
                  else
                    derv2(n22,n21a) = derv2(n22,n21a) + t24*dble(neqj)
                  endif
                endif
                if (lasu3.and.lopk) then
                  if (lopi) then
                    derv2(n11,n32a) = derv2(n11,n32a) + t13*dble(neqk)
                  else
                    derv2(n31,n32a) = derv2(n31,n32a) + t13*dble(neqk)
                  endif
                  if (lopj) then
                    derv2(n21,n32a) = derv2(n21,n32a) + t23*dble(neqk)
                  else
                    derv2(n31,n32a) = derv2(n31,n32a) + t23*dble(neqk)
                  endif
                  if (lopl) then
                    derv2(n42,n31a) = derv2(n42,n31a) + t34*dble(neqk)
                  else
                    derv2(n32,n31a) = derv2(n32,n31a) + t34*dble(neqk)
                  endif
                endif
                if (lasu4.and.lopl) then
                  if (lopi) then
                    derv2(n11,n42a) = derv2(n11,n42a) + t14*dble(neql)
                  else
                    derv2(n41,n42a) = derv2(n41,n42a) + t14*dble(neql)
                  endif
                  if (lopj) then
                    derv2(n21,n42a) = derv2(n21,n42a) + t24*dble(neql)
                  else
                    derv2(n41,n42a) = derv2(n41,n42a) + t24*dble(neql)
                  endif
                  if (lopk) then
                    derv2(n31,n42a) = derv2(n31,n42a) + t34*dble(neql)
                  else
                    derv2(n41,n42a) = derv2(n41,n42a) + t34*dble(neql)
                  endif
                endif
              enddo
            enddo
          endif
        enddo
120   continue
    enddo lkloop
  enddo ljloop
!
!  Closing banner for energy decomposition
!
  if (lPrintFour) then
    call mpbarrier
    if (ioproc) then
      write(ioout,'(''--------------------------------------------------------------------------------'',/)')
    endif
  endif
!
!  Free local memory
!
  deallocate(ioptaptr,stat=status)
  if (status/=0) call deallocate_error('fourimpsd2','ioptaptr')
  deallocate(nuniquelptr,stat=status)
  if (status/=0) call deallocate_error('fourimpsd2','nuniquelptr')
  deallocate(nuniquekptr,stat=status)
  if (status/=0) call deallocate_error('fourimpsd2','nuniquekptr')
  deallocate(nuniqueiptr,stat=status)
  if (status/=0) call deallocate_error('fourimpsd2','nuniqueiptr')
  deallocate(nptrnforimp,stat=status)
  if (status/=0) call deallocate_error('fourimpsd2','nptrnforimp')
  deallocate(ntypmiddle,stat=status)
  if (status/=0) call deallocate_error('fourimpsd2','ntypmiddle')
  deallocate(natmiddle,stat=status)
  if (status/=0) call deallocate_error('fourimpsd2','natmiddle')
#ifdef TRACE
  call trace_out('fourimpsd2')
#endif
!
  return
  end
