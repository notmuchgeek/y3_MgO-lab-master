  subroutine fourfc3(maxlim,maxlimcg,fc3)
!
!  Subroutine for third derivatives of four-body energy.
!  Called from thirdorderfc3.
!
!   5/15 Created from four3d3 / threefc3
!   2/18 Trace added
!   4/19 Call to lmatchpair modified to include wildcard argument
!   9/19 lnotorXduplicates flag added
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
!  Copyright Curtin University 2019
!
!  Julian Gale, CIC, Curtin University, September 2019
!
  use g_constants
  use configurations, only : nsuperghost
  use control,        only : lnotorXduplicates
  use current
  use derivatives
  use element,        only : maxele
  use four
  use mdlogic
  use molecule
  use optimisation
  use symmetry
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4),                   intent(in)    :: maxlim
  integer(i4),                   intent(in)    :: maxlimcg
  real(dp),                      intent(inout) :: fc3(maxlim,maxlim,maxlimcg)
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ii
  integer(i4)                                  :: iloop
  integer(i4)                                  :: iloopmin
  integer(i4)                                  :: imax
  integer(i4)                                  :: isgn
  integer(i4)                                  :: indmi
  integer(i4)                                  :: indmj
  integer(i4)                                  :: indmk
  integer(i4)                                  :: indml
  integer(i4),                            save :: ind2ptr(6,6)
  integer(i4),                            save :: ind3ptr(6,6,6)
  integer(i4)                                  :: ia
  integer(i4)                                  :: ic
  integer(i4)                                  :: ig
  integer(i4)                                  :: igx
  integer(i4)                                  :: igxyz
  integer(i4)                                  :: ii1
  integer(i4)                                  :: ii2
  integer(i4)                                  :: ii3
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
  integer(i4)                                  :: ix
  integer(i4)                                  :: ixx
  integer(i4)                                  :: iyy
  integer(i4)                                  :: izz
  integer(i4)                                  :: j
  integer(i4)                                  :: ja
  integer(i4)                                  :: jc
  integer(i4)                                  :: jg
  integer(i4)                                  :: jgx
  integer(i4)                                  :: jj
  integer(i4)                                  :: jmax
  integer(i4)                                  :: jx
  integer(i4)                                  :: jxx
  integer(i4)                                  :: jxyz
  integer(i4)                                  :: jyy
  integer(i4)                                  :: jzz
  integer(i4)                                  :: k
  integer(i4)                                  :: ka
  integer(i4)                                  :: kc
  integer(i4)                                  :: kg
  integer(i4)                                  :: kgx
  integer(i4)                                  :: kk
  integer(i4)                                  :: kloop
  integer(i4)                                  :: kloopmin
  integer(i4)                                  :: kmax
  integer(i4)                                  :: kx
  integer(i4)                                  :: kxx
  integer(i4)                                  :: kxyz
  integer(i4)                                  :: kyy
  integer(i4)                                  :: kzz
  integer(i4)                                  :: l
  integer(i4)                                  :: l4
  integer(i4)                                  :: lg
  integer(i4)                                  :: lgx
  integer(i4)                                  :: li
  integer(i4)                                  :: lk
  integer(i4)                                  :: ll
  integer(i4)                                  :: lloop
  integer(i4)                                  :: lloopmin
  integer(i4)                                  :: lx
  integer(i4)                                  :: lu
  integer(i4),                            save :: maxvector = 27
  integer(i4)                                  :: n
  integer(i4)                                  :: n2
  integer(i4), dimension(:), allocatable       :: natmiddle
  integer(i4), dimension(:), allocatable       :: ntypmiddle
  integer(i4)                                  :: nbtypeji
  integer(i4)                                  :: nbtypejk
  integer(i4)                                  :: nbtypekl
  integer(i4)                                  :: nbtypeji2
  integer(i4)                                  :: nbtypejk2
  integer(i4)                                  :: nbtypekl2
  integer(i4)                                  :: nfortype
  integer(i4)                                  :: nfornonoop
  integer(i4)                                  :: nghostcell    ! Number of ghost cells in the full cell
  integer(i4)                                  :: ni
  integer(i4)                                  :: nil
  integer(i4)                                  :: nil2
  integer(i4)                                  :: niltor
  integer(i4)                                  :: nimproper
  integer(i4)                                  :: nj
  integer(i4)                                  :: nk
  integer(i4)                                  :: nl 
  integer(i4)                                  :: nm
  integer(i4)                                  :: nmatch
  integer(i4)                                  :: nmi
  integer(i4)                                  :: nmiddle
  integer(i4)                                  :: nmid
  integer(i4)                                  :: nmj 
  integer(i4)                                  :: nmk 
  integer(i4)                                  :: nml 
  integer(i4)                                  :: nn
  integer(i4)                                  :: noofp
  integer(i4)                                  :: npha
  integer(i4), dimension(:), allocatable       :: nptrnfornonoop
  integer(i4)                                  :: nt1
  integer(i4)                                  :: nt2
  integer(i4)                                  :: nt3
  integer(i4)                                  :: nt4
  integer(i4)                                  :: ntmp
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
  integer(i4)                                  :: nvector
  integer(i4)                                  :: nwild1
  integer(i4)                                  :: nwild2
  integer(i4)                                  :: status
  logical                                      :: l2bondsij
  logical                                      :: l2bondsjk
  logical                                      :: l2bondskl
  logical                                      :: lallbtyp
  logical                                      :: lanybtyp
  logical                                      :: lanyneedmol
  logical                                      :: lbondedij
  logical                                      :: lbondedjk
  logical                                      :: lbondedkl
  logical                                      :: lbtyp
  logical                                      :: lexactmatch
  logical,                                save :: lfirst = .true.
  logical                                      :: lghosti
  logical                                      :: lghostj
  logical                                      :: lghostk
  logical                                      :: lghostl
  logical                                      :: libond
  logical                                      :: liok
  logical                                      :: limatch1
  logical                                      :: limatch4
  logical                                      :: lintra_only
  logical                                      :: linter_only
  logical                                      :: ljmatch2
  logical                                      :: ljmatch3
  logical                                      :: ljkmatch
  logical                                      :: lkbond
  logical                                      :: lkjmatch
  logical                                      :: lkmatch2
  logical                                      :: lkmatch3
  logical                                      :: llbond
  logical                                      :: lmatch
  logical                                      :: lmatchany
  logical                                      :: lmatchpair
  logical                                      :: lmeither
  logical                                      :: lmolloc
  logical                                      :: lmolok
  logical                                      :: lmolokjk
  logical                                      :: lneedmol
  logical                                      :: lsamemolij
  logical                                      :: lsamemoljk
  logical                                      :: lsamemolkl
  logical                                      :: lswitchil
  logical                                      :: lswitchjk
  logical                                      :: ltsyme_exact
  logical                                      :: lunique
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: cut
  real(dp)                                     :: cutmax
  real(dp),                               save :: delta(4,4,6)
  real(dp)                                     :: e1d(6)
  real(dp)                                     :: e2d(21)
  real(dp)                                     :: e3d(56)
  real(dp)                                     :: eterm
  real(dp)                                     :: fc33(3,4,3,4,3,4)    ! Full third derivative matrix
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
  real(dp)                                     :: rc(3,4,6)
  real(dp)                                     :: rkforloc
  real(dp)                                     :: rko
  real(dp)                                     :: rn
  real(dp)                                     :: rtmp
  real(dp)                                     :: time1
  real(dp)                                     :: time2
  real(dp)                                     :: tr1
  real(dp)                                     :: tr2
  real(dp)                                     :: tr2max
  real(dp)                                     :: tr3
  real(dp)                                     :: tr4
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
  real(dp)                                     :: x43t
  real(dp)                                     :: y43t
  real(dp)                                     :: z43t
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
  real(dp), dimension(:), allocatable,    save :: xvec
  real(dp), dimension(:), allocatable,    save :: yvec
  real(dp), dimension(:), allocatable,    save :: zvec
!
  data ind2ptr/1_i4,2_i4,3_i4,4_i4,5_i4,6_i4,2_i4,7_i4,8_i4,9_i4,10_i4,11_i4, &
               3_i4,8_i4,12_i4,13_i4,14_i4,15_i4,4_i4,9_i4,13_i4,16_i4,17_i4, &
               18_i4,5_i4,10_i4,14_i4,17_i4,19_i4,20_i4,6_i4,11_i4,15_i4,18_i4,20_i4,21_i4/
#ifdef TRACE
  call trace_in('fourfc3')
#endif
  if (lfirst) then
    lfirst = .false.
!
!  Set up ind3ptr
!
    ii = 0
    do i = 1,6
      do j = i,6
        do k = j,6
          ii = ii + 1
          ind3ptr(k,j,i) = ii
          ind3ptr(j,k,i) = ii
          ind3ptr(k,i,j) = ii
          ind3ptr(i,k,j) = ii
          ind3ptr(j,i,k) = ii
          ind3ptr(i,j,k) = ii
        enddo
      enddo
    enddo
!
!  Set up delta term array
!
    delta(1:4,1:4,1:6) = 0.0_dp
!
    delta(1,1,1) = 1.0_dp
    delta(2,1,1) = -1.0_dp
    delta(1,2,1) = -1.0_dp
    delta(2,2,1) = 1.0_dp
!
    delta(1,1,2) = 1.0_dp
    delta(3,1,2) = -1.0_dp
    delta(1,3,2) = -1.0_dp
    delta(3,3,2) = 1.0_dp
!
    delta(1,1,3) = 1.0_dp
    delta(4,1,3) = -1.0_dp
    delta(1,4,3) = -1.0_dp
    delta(4,4,3) = 1.0_dp
!
    delta(2,2,4) = 1.0_dp
    delta(3,2,4) = -1.0_dp
    delta(2,3,4) = -1.0_dp
    delta(3,3,4) = 1.0_dp
!
    delta(2,2,5) = 1.0_dp
    delta(4,2,5) = -1.0_dp
    delta(2,4,5) = -1.0_dp
    delta(4,4,5) = 1.0_dp
!
    delta(3,3,6) = 1.0_dp
    delta(3,4,6) = -1.0_dp
    delta(4,3,6) = -1.0_dp
    delta(4,4,6) = 1.0_dp
  endif
!
  time1 = g_cpu_time()
  lmolloc = (nmol.gt.0)
!
!  Find number of ghost cells
!
  nghostcell = nsuperghost(1,ncf)*nsuperghost(2,ncf)*nsuperghost(3,ncf)
!
!  Allocate local memory
!
  allocate(natmiddle(2*nfor),stat=status)
  if (status/=0) call outofmemory('fourfc3','natmiddle')
  allocate(ntypmiddle(2*nfor),stat=status)
  if (status/=0) call outofmemory('fourfc3','ntypmiddle')
  allocate(nptrnfornonoop(nfor),stat=status)
  if (status/=0) call outofmemory('fourfc3','nptrnfornonoop')
  allocate(nuniqueiptr(numat),stat=status)
  if (status/=0) call outofmemory('fourfc3','nuniqueiptr')
  allocate(nuniquekptr(numat),stat=status)
  if (status/=0) call outofmemory('fourfc3','nuniquekptr')
  allocate(nuniquelptr(numat),stat=status)
  if (status/=0) call outofmemory('fourfc3','nuniquelptr')
  if (ndim.gt.0) then
    allocate(xvec(maxvector),stat=status)
    if (status/=0) call outofmemory('fourfc3','xvec')
    allocate(yvec(maxvector),stat=status)
    if (status/=0) call outofmemory('fourfc3','yvec')
    allocate(zvec(maxvector),stat=status)
    if (status/=0) call outofmemory('fourfc3','zvec')
  else
    allocate(xvec(1),stat=status)
    if (status/=0) call outofmemory('fourfc3','xvec')
    allocate(yvec(1),stat=status)
    if (status/=0) call outofmemory('fourfc3','yvec')
    allocate(zvec(1),stat=status)
    if (status/=0) call outofmemory('fourfc3','zvec')
  endif
!
!  Check how many four-body potentials are of out of plane type and how many aren't
!
  noofp = 0
  nimproper = 0
  nfornonoop = 0
  do n = 1,nfor
    if (loutofplane(n)) then
      noofp = noofp + 1
    elseif (mmfexc(n).eq.2) then
      nimproper = nimproper + 1
    else
      nfornonoop = nfornonoop + 1
      nptrnfornonoop(nfornonoop) = n
    endif
  enddo
!
!  Find out if any require molecule information and whether any potential is of bonded type
!
  lallbtyp = .true.
  lanybtyp = .false.
  lanyneedmol = .false.
  do n = 1,nfor
    if (.not.loutofplane(n)) then
      lbtyp = (mmfexc(n).ge.1)
      lintra_only = (lfintra(n).and..not.lfinter(n))
      linter_only = (lfinter(n).and..not.lfintra(n))
      lneedmol = (lintra_only.or.linter_only.or.lbtyp)
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
    if (.not.loutofplane(n)) then
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
    if (.not.loutofplane(n)) then
      cut = for1(n) + for2(n) + for3(n)
      if (for4(n).gt.0.0_dp) cut = for4(n)
      cutmax = max(cut,cutmax)
      tr2 = for2(n)**2
      tr2max = max(tr2,tr2max)
    endif
  enddo
!
!  Create lattice vectors
!
  if (ndim.gt.0) then
    call rtlist(nvector,cutmax,xvec,yvec,zvec,imax,jmax,kmax,nmid,maxvector)
    if (nvector.gt.maxvector) then
!
!  Too many vectors
!
      deallocate(zvec,stat=status)
      if (status/=0) call deallocate_error('fourfc3','zvec')
      deallocate(yvec,stat=status)
      if (status/=0) call deallocate_error('fourfc3','yvec')
      deallocate(xvec,stat=status)
      if (status/=0) call deallocate_error('fourfc3','xvec')
      maxvector = nint(1.1*nvector)
      allocate(xvec(maxvector),stat=status)
      if (status/=0) call outofmemory('fourfc3','xvec')
      allocate(yvec(maxvector),stat=status)
      if (status/=0) call outofmemory('fourfc3','yvec')
      allocate(zvec(maxvector),stat=status)
      if (status/=0) call outofmemory('fourfc3','zvec')
      call rtlist(nvector,cutmax,xvec,yvec,zvec,imax,jmax,kmax,nmid,maxvector)
    endif
  else
    nvector = 1
    nmid = 1
    xvec(1) = 0.0_dp
    yvec(1) = 0.0_dp
    zvec(1) = 0.0_dp
  endif
!****************************************************************************
!  If there are no non out of plane potentials then we can skip everything  *
!****************************************************************************
  if (nfornonoop.eq.0) goto 5
!********************************
!  Loop over middle site 2 / j  *
!********************************
  ljloop:  do j = 1,numat
    nj = nat(j)
    ntypj = nftype(j)
    ocj = occuf(j)
!
!  Check whether species may be valid
!
    if (.not.lmatchany(nj,ntypj,nmiddle,natmiddle,ntypmiddle,.true.)) cycle ljloop
!
!  Is j a ghost atom?
!
    lghostj = (mod(j-1,nghostcell).eq.0.and.nj.le.maxele)
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
      kloopmin = 1
    else
      lkbond = .false.
      kloopmin = j
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
    endif
    xc2 = xclat(j)
    yc2 = yclat(j)
    zc2 = zclat(j)
!***************************************
!  Loop over second middle site 3 / k  *
!***************************************
    lkloop: do lk = kloopmin,kloop
      if (lkbond) then
        k = nuniquekptr(lk)
      else
        k = lk
      endif
!
!  Only do upper triangular set of atoms
!
      if (k.lt.j) cycle lkloop
!
      nk = nat(k)
      ntypk = nftype(k)
!
!  Is k a ghost atom?
!
      lghostk = (mod(k-1,nghostcell).eq.0.and.nk.le.maxele)
!
!  Check whether species may be valid
!
      if (.not.lmatchany(nk,ntypk,nmiddle,natmiddle,ntypmiddle,.true.)) cycle lkloop
      if (.not.lmatchpair(nj,ntypj,nk,ntypk,nfor,nfspec2,nfptyp2,nfspec3,nfptyp3,.true.)) cycle lkloop
!
      ock = occuf(k)
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
      else
        lmolokjk = .false.
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
      do 120 jj = 1,nvector
        r322 = (xvec(jj)+x32t)**2 + (yvec(jj)+y32t)**2 + (zvec(jj)+z32t)**2
        if (r322.lt.1d-12) goto 120
        if (k.eq.j.and.jj.eq.nmid) goto 120
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
            call lintoijk(jxx,jyy,jzz,jj,imax,jmax,kmax)
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
        xc3 = xc3t + xvec(jj)
        yc3 = yc3t + yvec(jj)
        zc3 = zc3t + zvec(jj)
        x32 = x32t + xvec(jj)
        y32 = y32t + yvec(jj)
        z32 = z32t + zvec(jj)
!
!  Set counter for number of valid i/l end atom combinations
!
        niltor = 0
!***********************************
!  Loop over four-body potentials  *
!***********************************
        pots: do nn = 1,nfornonoop
          n = nptrnfornonoop(nn)
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
          ltsyme_exact = lexactmatch(nt1,ntyp1,nt4,ntyp4)
          lbtyp = (mmfexc(n).ge.1)
          lintra_only = (lfintra(n).and..not.lfinter(n))
          linter_only = (lfinter(n).and..not.lfintra(n))
          lneedmol = (lintra_only.or.linter_only.or.lbtyp)
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
          ljmatch3 = lmatch(nj,ntypj,nt3,ntyp3,.true.)
          lkmatch2 = lmatch(nk,ntypk,nt2,ntyp2,.true.)
          lkmatch3 = lmatch(nk,ntypk,nt3,ntyp3,.true.)
!
!  Check whether j-k or k-j orders are OK for 2-3
!
          ljkmatch = (ljmatch2.and.lkmatch3)
          lkjmatch = (ljmatch3.and.lkmatch2)
!
!  If no pair of matches can be found then cycle
!
          if (.not.ljkmatch.and..not.lkjmatch) cycle pots
          lswitchil = .false.
          if (.not.ljkmatch) then
!
!  If j-k doesn't match, but k-j does then swap terms
!
            ntmp = nt2
            nt2 = nt3
            nt3 = ntmp
            ntmp = ntyp2
            ntyp2 = ntyp3
            ntyp3 = ntmp
            rtmp = tr1
            tr1 = tr3
            tr3 = rtmp
            if (.not.ltsyme_exact) then
              ntmp = nt1
              nt1 = nt4
              nt4 = ntmp
              ntmp = ntyp1
              ntyp1 = ntyp4
              ntyp4 = ntmp
            endif
            lswitchil = .true.
          endif
!
!  Set flag indicating whether middle atoms could be matched either way round
!
          lmeither = (ljkmatch.and.lkjmatch)
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
            iloopmin = 1
          else
            libond = .false.
            iloopmin = 1
            iloop = numat
          endif
!
!  Skip if iloop is zero
!
          if (iloop.eq.0) cycle pots
!*****************************
!  Loop over end site 1 / i  *
!*****************************
          liloop: do li = iloopmin,iloop
            if (libond) then
              i = nuniqueiptr(li)
            else
              i = li
            endif
            ni = nat(i)
            ntypi = nftype(i)
            oci = occuf(i)
!
!  Check whether i matches either of types 1 and 4
!
            limatch1 = lmatch(ni,ntypi,nt1,ntyp1,.true.)
            limatch4 = lmatch(ni,ntypi,nt4,ntyp4,.true.)
!
!  Is i allowed for type 1, or type 4 if the middle atoms can be switched?
!
            liok = (limatch1.or.(limatch4.and.lmeither))
            if (.not.liok) cycle liloop
!
!  Is i a ghost atom?
!
            lghosti = (mod(i-1,nghostcell).eq.0.and.ni.le.maxele)
!
            lswitchjk = .false.
            if (.not.limatch1.and.(limatch4.and.lmeither)) then
!
!  Switch round order of torsional atoms
!
              ntmp = nt1
              nt1 = nt4
              nt4 = ntmp
              ntmp = ntyp1
              ntyp1 = ntyp4
              ntyp4 = ntmp
              rtmp = tr1
              tr1 = tr3
              tr3 = rtmp
              lswitchjk = .true.
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
            else
              lmolok = .false.
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
            iiloop: do ii = 1,nvector
              r212 = (-xvec(ii) + x21t)**2 + (-yvec(ii) + y21t)**2 + (-zvec(ii) + z21t)**2
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
                  call lintoijk(ixx,iyy,izz,ii,imax,jmax,kmax)
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
              x21 = x21t - xvec(ii)
              y21 = y21t - yvec(ii)
              z21 = z21t - zvec(ii)
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
                if (nbonds(k).gt.0) then
                  nuniquel = 1
                  nuniquelptr(1) = nbonded(1,k)
                  do lloop = 2,nbonds(k)
                    lunique = .true.
                    do lu = 1,nuniquel
                      if (nbonded(lloop,k).eq.nuniquelptr(lu)) lunique = .false.
                    enddo
                    if (lunique) then
                      nuniquel = nuniquel + 1
                      nuniquelptr(nuniquel) = nbonded(lloop,k)
                    endif
                  enddo
                  lloop = nuniquel
                else
                  lloop = 0
                endif
                lloopmin = 1
              else
                llbond = .false.
                lloopmin = 1
                lloop = numat
              endif
!
!  Skip if lloop is zero
!
              if (lloop.eq.0) cycle iiloop
!**********************************
!  Loop over last end site 4 / l  *
!**********************************
              l4loop: do l4 = lloopmin,lloop
                if (llbond) then
                  l = nuniquelptr(l4)
                else
                  l = l4
                endif
                nl = nat(l)
                ntypl = nftype(l)
                ocl = occuf(l)
!
!  Check l is allowed for n
!             
                if (.not.lmatch(nl,ntypl,nt4,ntyp4,.true.)) cycle l4loop
!
!  Is l a ghost atom?
!
                lghostl = (mod(l-1,nghostcell).eq.0.and.nl.le.maxele)
!
!  If none of the atoms are ghosts then cycle
!
                if (.not.lghosti.and..not.lghostj.and..not.lghostk.and..not.lghostl) cycle l4loop
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
                else
                  lmolok = .false.
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
                x43t = xc4t - xc3
                y43t = yc4t - yc3
                z43t = zc4t - zc3
!
!  Check r43 is OK
!  Loop over cell vectors
!
                llloop: do ll = 1,nvector
                  r432 = (xvec(ll)+x43t)**2 + (yvec(ll)+y43t)**2 + (zvec(ll)+z43t)**2
                  if (r432.lt.1d-12) cycle llloop
!
!  Molecule checking
!
                  lbondedkl = .false.
                  if (lmolok) then
                    if (ndim.eq.0) then
                      if (linter_only) cycle llloop
                      if (lbtyp) then
                        call bonded(lbondedkl,l2bondskl,nbtypekl,nbtypekl2,k,l,0_i4,0_i4,0_i4)
                        if (.not.lbondedkl) cycle llloop
                      endif
                    else
                      call lintoijk(kxx,kyy,kzz,ll,imax,jmax,kmax)
                      if (lbtyp) then
                        call bonded(lbondedkl,l2bondskl,nbtypekl,nbtypekl2,k,l,kxx-jxx,kyy-jyy,kzz-jzz)
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
                    if (r432.gt.tr3) cycle llloop
                  endif
!
                  x43 = x43t + xvec(ll)
                  y43 = y43t + yvec(ll)
                  z43 = z43t + zvec(ll)
!
!  Check r41 is OK
!
                  x41 = x43 + x32 + x21
                  y41 = y43 + y32 + y21
                  z41 = z43 + z32 + z21
                  r412 = x41*x41 + y41*y41 + z41*z41
                  if (r412.gt.tr4.and.tr4.gt.0.0_dp.and..not.lbtyp) cycle llloop
                  if (r412.lt.1.0d-12) cycle llloop
!
!  Check r42 is OK
!
                  x42 = x32 + x43
                  y42 = y32 + y43
                  z42 = z32 + z43
                  r422 = x42*x42 + y42*y42 + z42*z42
                  if (r422.lt.1.0d-12) cycle llloop
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
                  liltorswitch(niltor) = lswitchil
                  ljktorswitch(niltor) = lswitchjk
!
                  iltor(1,niltor) = i
                  iltor(2,niltor) = l
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
                  oiltor(niltor) = oci*ocj*ock*ocl
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
!  Remove wildcard duplicates  *
!*******************************
        if (lnotorXduplicates) then
          lkeeptor(1:niltor) = .true.
          do nil = 1,niltor-1
            if (lkeeptor(nil)) then
              i = iltor(1,nil)
              l = iltor(2,nil)
              n = nfortor(nil)
              nwild1 = 0
              if (nfspec1(n).eq.maxele) nwild1 = nwild1 + 1
              if (nfspec4(n).eq.maxele) nwild1 = nwild1 + 1
!
!  Loop over remaining torsions looking for matches
!
              nmatch = 0
              do nil2 = nil+1,niltor
                if (iltor(1,nil2).eq.i.and.iltor(2,nil2).eq.l) then
!
!  Add duplicate potential to the list
!
                  nmatch = nmatch + 1
                  n2 = nfortor(nil2)
                  nwild2 = 0
                  if (nfspec1(n2).eq.maxele) nwild2 = nwild2 + 1
                  if (nfspec4(n2).eq.maxele) nwild2 = nwild2 + 1
                  nduptor(nmatch) = nil2
                  nwildduptor(nmatch) = nwild2
                endif
              enddo
            endif
!
!  If there are matches then set the flag to which one will be kept
!
            nkeepfor = nil
            if (nmatch.gt.0) then
              do nm = 1,nmatch
                if (nwildduptor(nm).lt.nwild1) then
                  nwild1 = nwildduptor(nm)
                  nkeepfor = nduptor(nm)
                endif
              enddo
!
!  Change flags for potentials not to be kept
!
              if (nkeepfor.eq.nil) then
                do nm = 1,nmatch
                  lkeeptor(nduptor(nm)) = .false.
                enddo
              else
                lkeeptor(nil) = .false.
                do nm = 1,nmatch
                  if (nkeepfor.ne.nduptor(nm)) then
                    lkeeptor(nduptor(nm)) = .false.
                  endif
                enddo
              endif
            endif
          enddo
!
!  Count the number of actual torsions to be computed
!
          nkeepfor = 0
          do nil = 1,niltor
            if (lkeeptor(nil)) then
              nkeepfor = nkeepfor + 1
              nkeeptor(nkeepfor) = nil
            endif
          enddo
        else
          nkeepfor = 0
          do nil = 1,niltor
            nkeepfor = nkeepfor + 1
            nkeeptor(nkeepfor) = nil
          enddo
        endif
!*******************************
!  Loop over i/l combinations  *
!*******************************
        do nil2 = 1,nkeepfor
          nil = nkeeptor(nil2)
!
!  Return values to local variables
!
          n = nfortor(nil)
          lswitchil = liltorswitch(nil)
          lswitchjk = ljktorswitch(nil)
!
          i = iltor(1,nil)
          l = iltor(2,nil)
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
          ofct = oiltor(nil)
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
!  Switch terms if necessary
!
          if (lswitchil) then
            if (nfortype.eq.8.or.nfortype.eq.9) then
              rtmp = fpoly(3)
              fpoly(3) = fpoly(5)
              fpoly(5) = rtmp
            elseif (nfortype.eq.6.or.nfortype.eq.7) then
              rtmp = fpoly(2)
              fpoly(2) = fpoly(4)
              fpoly(4) = rtmp
            elseif (nfortype.eq.10.or.nfortype.eq.17) then
              rtmp = fpoly(2)
              fpoly(2) = fpoly(1)
              fpoly(1) = rtmp
            endif
          endif
          if (lswitchjk) then
            if (nfortype.eq.8.or.nfortype.eq.9) then
              rtmp = fpoly(3)
              fpoly(3) = fpoly(5)
              fpoly(5) = rtmp
            elseif (nfortype.eq.6.or.nfortype.eq.7) then
              rtmp = fpoly(2)
              fpoly(2) = fpoly(4)
              fpoly(4) = rtmp
            elseif (nfortype.eq.10.or.nfortype.eq.17) then
              rtmp = fpoly(2)
              fpoly(2) = fpoly(1)
              fpoly(1) = rtmp
            endif
          endif
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
                        rko,rn,phi0o,isgn,fpoly,.true.,.true.,.true.)
!         
!  Define pointers to elements of second derivative matrices
!           
          ix = 3*(i-1) + 1
          jx = 3*(j-1) + 1
          kx = 3*(k-1) + 1
          lx = 3*(l-1) + 1
!
!  Is i a ghost atom?
!
          lghosti = (mod(i-1,nghostcell).eq.0.and.nat(i).le.maxele)
          if (lghosti) then
            ig = (i - 1)/nghostcell + 1
            igx = 3*(ig-1) + 1
          endif
!
!  Is j a ghost atom?
!
          lghostj = (mod(j-1,nghostcell).eq.0.and.nat(j).le.maxele)
          if (lghostj) then
            jg = (j - 1)/nghostcell + 1
            jgx = 3*(jg-1) + 1
          endif
!
!  Is k a ghost atom?
!
          lghostk = (mod(k-1,nghostcell).eq.0.and.nat(k).le.maxele)
          if (lghostk) then
            kg = (k - 1)/nghostcell + 1
            kgx = 3*(kg-1) + 1
          endif
!
!  Is l a ghost atom?
!
          lghostl = (mod(l-1,nghostcell).eq.0.and.nat(l).le.maxele)
          if (lghostl) then
            lg = (l - 1)/nghostcell + 1
            lgx = 3*(lg-1) + 1
          endif
!********************************
!  Calculate third derivatives  *
!********************************
!
!  Initialise full third derivative matrix
!
          fc33(1:3,1:4,1:3,1:4,1:3,1:4) = 0.0_dp
!
!  Set up array of Cartesian components for vectors
!
          rc(1:3,1:4,1:6) = 0.0_dp
!
          rc(1,1,1) = - x21
          rc(2,1,1) = - y21
          rc(3,1,1) = - z21
          rc(1,2,1) = x21
          rc(2,2,1) = y21
          rc(3,2,1) = z21
!
          rc(1,1,2) = - x31
          rc(2,1,2) = - y31
          rc(3,1,2) = - z31
          rc(1,3,2) = x31
          rc(2,3,2) = y31
          rc(3,3,2) = z31
!
          rc(1,1,3) = - x41
          rc(2,1,3) = - y41
          rc(3,1,3) = - z41
          rc(1,4,3) = x41
          rc(2,4,3) = y41
          rc(3,4,3) = z41
!
          rc(1,2,4) = - x32
          rc(2,2,4) = - y32
          rc(3,2,4) = - z32
          rc(1,3,4) = x32
          rc(2,3,4) = y32
          rc(3,3,4) = z32
!
          rc(1,2,5) = - x42
          rc(2,2,5) = - y42
          rc(3,2,5) = - z42
          rc(1,4,5) = x42
          rc(2,4,5) = y42
          rc(3,4,5) = z42
!
          rc(1,3,6) = - x43
          rc(2,3,6) = - y43
          rc(3,3,6) = - z43
          rc(1,4,6) = x43
          rc(2,4,6) = y43
          rc(3,4,6) = z43
!
!  Loop over atoms and Cartesian components to compute third derivatives
!
          do ia = 1,4   ! First atom
            do ic = 1,3   ! First atom Cartesian component
              do ja = 1,4   ! Second atom
                do jc = 1,3   ! Second atom Cartesian component
                  do ka = 1,4   ! Third atom
                    do kc = 1,3   ! Third atom Cartesian component
                      do ii1 = 1,6  ! Vectors between atoms
!
                        do ii2 = 1,6  ! Vectors between atoms
                          do ii3 = 1,6  ! Vectors between atoms
                            fc33(kc,ka,jc,ja,ic,ia) = fc33(kc,ka,jc,ja,ic,ia) + &
                              e3d(ind3ptr(ii3,ii2,ii1))*rc(kc,ka,ii3)*rc(jc,ja,ii2)*rc(ic,ia,ii1)
                          enddo
!
                          if (ic.eq.jc) then
                            fc33(kc,ka,jc,ja,ic,ia) = fc33(kc,ka,jc,ja,ic,ia) + &
                              e2d(ind2ptr(ii2,ii1))*delta(ja,ia,ii1)*rc(kc,ka,ii2)
                          endif
                          if (ic.eq.kc) then
                            fc33(kc,ka,jc,ja,ic,ia) = fc33(kc,ka,jc,ja,ic,ia) + &
                              e2d(ind2ptr(ii2,ii1))*delta(ka,ia,ii1)*rc(jc,ja,ii2)
                          endif
                          if (jc.eq.kc) then
                            fc33(kc,ka,jc,ja,ic,ia) = fc33(kc,ka,jc,ja,ic,ia) + &
                              e2d(ind2ptr(ii2,ii1))*delta(ka,ja,ii1)*rc(ic,ia,ii2)
                          endif
                        enddo
                      enddo
!
!  End of loops over atoms and Cartesian components
!
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
!
!  Copy relevant parts of full third derivative matrix to fc3
!
          if (lghosti) then
            igxyz = igx - 1
            do ic = 1,3
              igxyz = igxyz + 1
              do ja = 1,4
                if (ja.eq.1) then
                  jxyz = ix - 1
                elseif (ja.eq.2) then
                  jxyz = jx - 1
                elseif (ja.eq.3) then
                  jxyz = kx - 1
                else
                  jxyz = lx - 1
                endif
                do jc = 1,3
                  jxyz = jxyz + 1
                  do ka = 1,4
                    if (ka.eq.1) then
                      kxyz = ix - 1
                    elseif (ka.eq.2) then
                      kxyz = jx - 1
                    elseif (ka.eq.3) then
                      kxyz = kx - 1
                    else
                      kxyz = lx - 1
                    endif
                    do kc = 1,3
                      kxyz = kxyz + 1
                      fc3(kxyz,jxyz,igxyz) = fc3(kxyz,jxyz,igxyz) + fc33(kc,ka,jc,ja,ic,1)
                    enddo
                  enddo
                enddo
              enddo
            enddo
          endif
!
          if (lghostj) then
            igxyz = jgx - 1
            do ic = 1,3
              igxyz = igxyz + 1
              do ja = 1,4
                if (ja.eq.1) then
                  jxyz = ix - 1
                elseif (ja.eq.2) then
                  jxyz = jx - 1
                elseif (ja.eq.3) then
                  jxyz = kx - 1
                else
                  jxyz = lx - 1
                endif
                do jc = 1,3
                  jxyz = jxyz + 1
                  do ka = 1,4
                    if (ka.eq.1) then
                      kxyz = ix - 1
                    elseif (ka.eq.2) then
                      kxyz = jx - 1
                    elseif (ka.eq.3) then
                      kxyz = kx - 1
                    else
                      kxyz = lx - 1
                    endif
                    do kc = 1,3
                      kxyz = kxyz + 1
                      fc3(kxyz,jxyz,igxyz) = fc3(kxyz,jxyz,igxyz) + fc33(kc,ka,jc,ja,ic,2)
                    enddo
                  enddo
                enddo
              enddo
            enddo
          endif
!
          if (lghostk) then
            igxyz = kgx - 1 
            do ic = 1,3 
              igxyz = igxyz + 1 
              do ja = 1,4 
                if (ja.eq.1) then
                  jxyz = ix - 1 
                elseif (ja.eq.2) then
                  jxyz = jx - 1 
                elseif (ja.eq.3) then
                  jxyz = kx - 1
                else
                  jxyz = lx - 1 
                endif
                do jc = 1,3 
                  jxyz = jxyz + 1 
                  do ka = 1,4 
                    if (ka.eq.1) then
                      kxyz = ix - 1 
                    elseif (ka.eq.2) then
                      kxyz = jx - 1 
                    elseif (ka.eq.3) then
                      kxyz = kx - 1
                    else
                      kxyz = lx - 1 
                    endif
                    do kc = 1,3 
                      kxyz = kxyz + 1 
                      fc3(kxyz,jxyz,igxyz) = fc3(kxyz,jxyz,igxyz) + fc33(kc,ka,jc,ja,ic,3)
                    enddo
                  enddo
                enddo
              enddo
            enddo
          endif
!
          if (lghostl) then
            igxyz = lgx - 1 
            do ic = 1,3 
              igxyz = igxyz + 1 
              do ja = 1,4 
                if (ja.eq.1) then
                  jxyz = ix - 1 
                elseif (ja.eq.2) then
                  jxyz = jx - 1 
                elseif (ja.eq.3) then
                  jxyz = kx - 1
                else
                  jxyz = lx - 1 
                endif
                do jc = 1,3 
                  jxyz = jxyz + 1 
                  do ka = 1,4 
                    if (ka.eq.1) then
                      kxyz = ix - 1 
                    elseif (ka.eq.2) then
                      kxyz = jx - 1 
                    elseif (ka.eq.3) then
                      kxyz = kx - 1
                    else
                      kxyz = lx - 1 
                    endif
                    do kc = 1,3 
                      kxyz = kxyz + 1 
                      fc3(kxyz,jxyz,igxyz) = fc3(kxyz,jxyz,igxyz) + fc33(kc,ka,jc,ja,ic,4)
                    enddo
                  enddo
                enddo
              enddo
            enddo
          endif
!
        enddo
120   continue
    enddo lkloop
  enddo ljloop
!
!  End of outer loops
!
5 continue
!****************************
!  Out of plane potentials  *
!****************************
  if (noofp.gt.0) call fouroopfc3(maxlim,maxlimcg,fc3,ind2ptr,ind3ptr,delta)
!************************
!  Improper potentials  *
!************************
  if (nimproper.gt.0) call fourimpfc3(maxlim,maxlimcg,fc3,ind2ptr,ind3ptr,delta)
!
!  Free local memory
!
  deallocate(zvec,stat=status)
  if (status/=0) call deallocate_error('fourfc3','zvec')
  deallocate(yvec,stat=status)
  if (status/=0) call deallocate_error('fourfc3','yvec')
  deallocate(xvec,stat=status)
  if (status/=0) call deallocate_error('fourfc3','xvec')
  deallocate(nuniquelptr,stat=status)
  if (status/=0) call deallocate_error('fourfc3','nuniquelptr')
  deallocate(nuniquekptr,stat=status)
  if (status/=0) call deallocate_error('fourfc3','nuniquekptr')
  deallocate(nuniqueiptr,stat=status)
  if (status/=0) call deallocate_error('fourfc3','nuniqueiptr')
  deallocate(nptrnfornonoop,stat=status)
  if (status/=0) call deallocate_error('fourfc3','nptrnfornonoop')
  deallocate(ntypmiddle,stat=status)
  if (status/=0) call deallocate_error('fourfc3','ntypmiddle')
  deallocate(natmiddle,stat=status)
  if (status/=0) call deallocate_error('fourfc3','natmiddle')
!
!  Timing
!
  time2 = g_cpu_time()
  tfour = tfour + time2 - time1
#ifdef TRACE
  call trace_out('fourfc3')
#endif
!
  return
  end
