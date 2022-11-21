  subroutine four1fc(maxlhs,d1cell,matom,vmatom)
!
!  Subroutine for four-body first derivatives saved by supercell
!  for finite differencing to get the unphased force constants
!
!  12/14 Created from fourfc
!   2/18 Trace added
!   4/19 Call to lmatchpair modified to include wildcard argument
!   8/19 Corrections to arguments passed to routines called
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
  use control,       only : lnotorXduplicates
  use current
  use derivatives
  use element,       only : maxele
  use four
  use molecule
  use times
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4),                     intent(in)    :: maxlhs
  integer(i4),                     intent(in)    :: matom
  real(dp),                        intent(out)   :: d1cell(4,maxlhs,*)
  real(dp),                        intent(in)    :: vmatom(4)
!
!  Local variables
!
  integer(i4)                                    :: i
  integer(i4)                                    :: ii
  integer(i4)                                    :: iloop
  integer(i4)                                    :: iloopmin
  integer(i4)                                    :: imax
  integer(i4)                                    :: isgn
  integer(i4)                                    :: indmi
  integer(i4)                                    :: indmj
  integer(i4)                                    :: indmk
  integer(i4)                                    :: indml
  integer(i4), dimension(:,:), allocatable, save :: ivec
  integer(i4)                                    :: ixi
  integer(i4)                                    :: iyi
  integer(i4)                                    :: izi
  integer(i4)                                    :: ixj
  integer(i4)                                    :: iyj
  integer(i4)                                    :: izj
  integer(i4)                                    :: ixk
  integer(i4)                                    :: iyk
  integer(i4)                                    :: izk
  integer(i4)                                    :: ixl
  integer(i4)                                    :: iyl
  integer(i4)                                    :: izl
  integer(i4)                                    :: ixx
  integer(i4)                                    :: iyy
  integer(i4)                                    :: izz
  integer(i4)                                    :: j
  integer(i4)                                    :: jj
  integer(i4)                                    :: jmax
  integer(i4)                                    :: jxx
  integer(i4)                                    :: jyy
  integer(i4)                                    :: jzz
  integer(i4)                                    :: k
  integer(i4)                                    :: kk
  integer(i4)                                    :: kloop
  integer(i4)                                    :: kloopmin
  integer(i4)                                    :: kmax
  integer(i4)                                    :: kxx
  integer(i4)                                    :: kyy
  integer(i4)                                    :: kzz
  integer(i4)                                    :: l
  integer(i4)                                    :: l4
  integer(i4)                                    :: li
  integer(i4)                                    :: lk
  integer(i4)                                    :: ll
  integer(i4)                                    :: lloop
  integer(i4)                                    :: lloopmin
  integer(i4)                                    :: lu
  integer(i4),                              save :: maxvector = 27
  integer(i4)                                    :: n
  integer(i4)                                    :: n2
  integer(i4), dimension(:),   allocatable, save :: natmiddle
  integer(i4), dimension(:),   allocatable, save :: ntypmiddle
  integer(i4)                                    :: nbtypeji
  integer(i4)                                    :: nbtypejk
  integer(i4)                                    :: nbtypekl
  integer(i4)                                    :: nbtypeji2
  integer(i4)                                    :: nbtypejk2
  integer(i4)                                    :: nbtypekl2
  integer(i4)                                    :: ncind12m
  integer(i4)                                    :: ncind12p
  integer(i4)                                    :: ncind13m
  integer(i4)                                    :: ncind13p
  integer(i4)                                    :: ncind14m
  integer(i4)                                    :: ncind14p
  integer(i4)                                    :: ncind23m
  integer(i4)                                    :: ncind23p
  integer(i4)                                    :: ncind24m
  integer(i4)                                    :: ncind24p
  integer(i4)                                    :: ncind34m
  integer(i4)                                    :: ncind34p
  integer(i4)                                    :: nfortype
  integer(i4)                                    :: nfornonoop
  integer(i4)                                    :: ni
  integer(i4)                                    :: nil
  integer(i4)                                    :: nil2
  integer(i4)                                    :: niltor
  integer(i4)                                    :: nimproper
  integer(i4)                                    :: nj
  integer(i4)                                    :: nk
  integer(i4)                                    :: nl 
  integer(i4)                                    :: nm
  integer(i4)                                    :: nmatch
  integer(i4)                                    :: nmi
  integer(i4)                                    :: nmiddle
  integer(i4)                                    :: nmid
  integer(i4)                                    :: nmj 
  integer(i4)                                    :: nmk 
  integer(i4)                                    :: nml 
  integer(i4)                                    :: nn
  integer(i4)                                    :: noofp
  integer(i4)                                    :: npha
  integer(i4), dimension(:),   allocatable, save :: nptrnfornonoop
  integer(i4)                                    :: nt1
  integer(i4)                                    :: nt2
  integer(i4)                                    :: nt3
  integer(i4)                                    :: nt4
  integer(i4)                                    :: ntmp
  integer(i4)                                    :: ntyp1
  integer(i4)                                    :: ntyp2
  integer(i4)                                    :: ntyp3
  integer(i4)                                    :: ntyp4
  integer(i4)                                    :: ntypi
  integer(i4)                                    :: ntypj
  integer(i4)                                    :: ntypk
  integer(i4)                                    :: ntypl
  integer(i4)                                    :: nuniquei
  integer(i4)                                    :: nuniquek
  integer(i4)                                    :: nuniquel
  integer(i4), dimension(:),   allocatable, save :: nuniqueiptr
  integer(i4), dimension(:),   allocatable, save :: nuniquekptr
  integer(i4), dimension(:),   allocatable, save :: nuniquelptr
  integer(i4)                                    :: nvector
  integer(i4)                                    :: nwild1
  integer(i4)                                    :: nwild2
  integer(i4)                                    :: status
  logical                                        :: l2bondsij
  logical                                        :: l2bondsjk
  logical                                        :: l2bondskl
  logical                                        :: lallbtyp
  logical                                        :: lanybtyp
  logical                                        :: lanyneedmol
  logical                                        :: lbondedij
  logical                                        :: lbondedjk
  logical                                        :: lbondedkl
  logical                                        :: lbtyp
  logical                                        :: lexactmatch
  logical                                        :: libond
  logical                                        :: liok
  logical                                        :: limatch1
  logical                                        :: limatch4
  logical                                        :: lintra_only
  logical                                        :: linter_only
  logical                                        :: ljmatch2
  logical                                        :: ljmatch3
  logical                                        :: ljkmatch
  logical                                        :: lkbond
  logical                                        :: lkjmatch
  logical                                        :: lkmatch2
  logical                                        :: lkmatch3
  logical                                        :: llbond
  logical                                        :: lmatch
  logical                                        :: lmatchany
  logical                                        :: lmatchpair
  logical                                        :: lmeither
  logical                                        :: lmolloc
  logical                                        :: lmolok
  logical                                        :: lmolokjk
  logical                                        :: lneedmol
  logical                                        :: lsamemolij
  logical                                        :: lsamemoljk
  logical                                        :: lsamemolkl
  logical                                        :: lswitchil
  logical                                        :: lswitchjk
  logical                                        :: ltsyme_exact
  logical                                        :: lunique
  real(dp)                                       :: g_cpu_time
  real(dp)                                       :: cut
  real(dp)                                       :: cutmax
  real(dp)                                       :: e1d(6)
  real(dp)                                       :: e2d(21)
  real(dp)                                       :: e3d(56)
  real(dp)                                       :: eterm
  real(dp)                                       :: fpoly(5)
  real(dp)                                       :: oci
  real(dp)                                       :: ocj
  real(dp)                                       :: ock
  real(dp)                                       :: ocl 
  real(dp)                                       :: ofct 
  real(dp)                                       :: phi0 
  real(dp)                                       :: phi0o
  real(dp)                                       :: r21
  real(dp)                                       :: r21d
  real(dp)                                       :: r212
  real(dp)                                       :: r212d
  real(dp)                                       :: r31
  real(dp)                                       :: r31d
  real(dp)                                       :: r312
  real(dp)                                       :: r312d
  real(dp)                                       :: r32d
  real(dp)                                       :: r322
  real(dp)                                       :: r322d
  real(dp)                                       :: r41
  real(dp)                                       :: r41d
  real(dp)                                       :: r412
  real(dp)                                       :: r412d
  real(dp)                                       :: r42
  real(dp)                                       :: r42d
  real(dp)                                       :: r422
  real(dp)                                       :: r422d
  real(dp)                                       :: r43
  real(dp)                                       :: r43d
  real(dp)                                       :: r432
  real(dp)                                       :: r432d
  real(dp)                                       :: rkforloc
  real(dp)                                       :: rko
  real(dp)                                       :: rn
  real(dp)                                       :: rtmp
  real(dp)                                       :: time1
  real(dp)                                       :: time2
  real(dp)                                       :: tr1
  real(dp)                                       :: tr2
  real(dp)                                       :: tr2max
  real(dp)                                       :: tr3
  real(dp)                                       :: tr4
  real(dp)                                       :: x21
  real(dp)                                       :: y21
  real(dp)                                       :: z21
  real(dp)                                       :: x21d
  real(dp)                                       :: y21d
  real(dp)                                       :: z21d
  real(dp)                                       :: x21t
  real(dp)                                       :: y21t
  real(dp)                                       :: z21t
  real(dp)                                       :: x31
  real(dp)                                       :: y31
  real(dp)                                       :: z31
  real(dp)                                       :: x31d
  real(dp)                                       :: y31d
  real(dp)                                       :: z31d
  real(dp)                                       :: x41
  real(dp)                                       :: y41
  real(dp)                                       :: z41
  real(dp)                                       :: x41d
  real(dp)                                       :: y41d
  real(dp)                                       :: z41d
  real(dp)                                       :: x32
  real(dp)                                       :: y32
  real(dp)                                       :: z32
  real(dp)                                       :: x32d
  real(dp)                                       :: y32d
  real(dp)                                       :: z32d
  real(dp)                                       :: x32t
  real(dp)                                       :: y32t
  real(dp)                                       :: z32t
  real(dp)                                       :: x42
  real(dp)                                       :: y42
  real(dp)                                       :: z42
  real(dp)                                       :: x42d
  real(dp)                                       :: y42d
  real(dp)                                       :: z42d
  real(dp)                                       :: x43
  real(dp)                                       :: y43
  real(dp)                                       :: z43
  real(dp)                                       :: x43d
  real(dp)                                       :: y43d
  real(dp)                                       :: z43d
  real(dp)                                       :: x43t
  real(dp)                                       :: y43t
  real(dp)                                       :: z43t
  real(dp)                                       :: xc1t
  real(dp)                                       :: yc1t
  real(dp)                                       :: zc1t
  real(dp)                                       :: xc2
  real(dp)                                       :: yc2
  real(dp)                                       :: zc2
  real(dp)                                       :: xc3
  real(dp)                                       :: yc3
  real(dp)                                       :: zc3
  real(dp)                                       :: xc3t
  real(dp)                                       :: yc3t
  real(dp)                                       :: zc3t
  real(dp)                                       :: xc4t
  real(dp)                                       :: yc4t
  real(dp)                                       :: zc4t
  real(dp),    dimension(:),   allocatable, save :: xvec
  real(dp),    dimension(:),   allocatable, save :: yvec
  real(dp),    dimension(:),   allocatable, save :: zvec
#ifdef TRACE
  call trace_in('four1fc')
#endif
!
  time1 = g_cpu_time()
  lmolloc = (nmol.gt.0)
!
!  Allocate local memory
!
  allocate(natmiddle(2*nfor),stat=status)
  if (status/=0) call outofmemory('four1fc','natmiddle')
  allocate(ntypmiddle(2*nfor),stat=status)
  if (status/=0) call outofmemory('four1fc','ntypmiddle')
  allocate(nptrnfornonoop(nfor),stat=status)
  if (status/=0) call outofmemory('four1fc','nptrnfornonoop')
  allocate(nuniqueiptr(numat),stat=status)
  if (status/=0) call outofmemory('four1fc','nuniqueiptr')
  allocate(nuniquekptr(numat),stat=status)
  if (status/=0) call outofmemory('four1fc','nuniquekptr')
  allocate(nuniquelptr(numat),stat=status)
  if (status/=0) call outofmemory('four1fc','nuniquelptr')
  if (ndim.gt.0) then
    allocate(ivec(3,maxvector),stat=status)
    if (status/=0) call outofmemory('four1fc','ivec')
    allocate(xvec(maxvector),stat=status)
    if (status/=0) call outofmemory('four1fc','xvec')
    allocate(yvec(maxvector),stat=status)
    if (status/=0) call outofmemory('four1fc','yvec')
    allocate(zvec(maxvector),stat=status)
    if (status/=0) call outofmemory('four1fc','zvec')
  else
    allocate(ivec(3,1),stat=status)
    if (status/=0) call outofmemory('four1fc','ivec')
    allocate(xvec(1),stat=status)
    if (status/=0) call outofmemory('four1fc','xvec')
    allocate(yvec(1),stat=status)
    if (status/=0) call outofmemory('four1fc','yvec')
    allocate(zvec(1),stat=status)
    if (status/=0) call outofmemory('four1fc','zvec')
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
    call rtclist(nvector,cutmax,xvec,yvec,zvec,ivec,imax,jmax,kmax,nmid,maxvector)
    if (nvector.gt.maxvector) then
!
!  Too many vectors
!
      deallocate(zvec,stat=status)
      if (status/=0) call deallocate_error('four1fc','zvec')
      deallocate(yvec,stat=status)
      if (status/=0) call deallocate_error('four1fc','yvec')
      deallocate(xvec,stat=status)
      if (status/=0) call deallocate_error('four1fc','xvec')
      deallocate(ivec,stat=status)
      if (status/=0) call deallocate_error('four1fc','ivec')
      maxvector = nint(1.1*nvector)
      allocate(ivec(3,maxvector),stat=status)
      if (status/=0) call outofmemory('four1fc','ivec')
      allocate(xvec(maxvector),stat=status)
      if (status/=0) call outofmemory('four1fc','xvec')
      allocate(yvec(maxvector),stat=status)
      if (status/=0) call outofmemory('four1fc','yvec')
      allocate(zvec(maxvector),stat=status)
      if (status/=0) call outofmemory('four1fc','zvec')
      call rtclist(nvector,cutmax,xvec,yvec,zvec,ivec,imax,jmax,kmax,nmid,maxvector)
    endif
  else
    nvector = 1
    nmid = 1
    ivec(1,1) = 0
    ivec(2,1) = 0
    ivec(3,1) = 0
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
!  Find cell index for 2-3
!
        if (abs(ivec(1,jj)).gt.nd2cell(1).or. &
            abs(ivec(2,jj)).gt.nd2cell(2).or. &
            abs(ivec(3,jj)).gt.nd2cell(3)) then
!
!  Find cell index - if outside user range then assign to central cell
!
          ncind23m = nd2central
          ncind23p = nd2central
        else
!
!  Compute index
!
          ncind23p = nd2cellptr(nd2cell(1)+1-ivec(1,jj),nd2cell(2)+1-ivec(2,jj),nd2cell(3)+1-ivec(3,jj))
          ncind23m = nd2cellptr(nd2cell(1)+1+ivec(1,jj),nd2cell(2)+1+ivec(2,jj),nd2cell(3)+1+ivec(3,jj))
        endif
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
!
!  Find cell index for 1-2
!
              if (abs(ivec(1,ii)).gt.nd2cell(1).or. &
                  abs(ivec(2,ii)).gt.nd2cell(2).or. &
                  abs(ivec(3,ii)).gt.nd2cell(3)) then
!
!  Find cell index - if outside user range then assign to central cell
!
                ncind12m = nd2central
                ncind12p = nd2central
              else
!
!  Compute index
!
                ncind12p = nd2cellptr(nd2cell(1)+1-ivec(1,ii),nd2cell(2)+1-ivec(2,ii),nd2cell(3)+1-ivec(3,ii))
                ncind12m = nd2cellptr(nd2cell(1)+1+ivec(1,ii),nd2cell(2)+1+ivec(2,ii),nd2cell(3)+1+ivec(3,ii))
              endif
!
!  Find cell index for 1-3
!
              if (abs(ivec(1,jj)-ivec(1,ii)).gt.nd2cell(1).or. &
                  abs(ivec(2,jj)-ivec(2,ii)).gt.nd2cell(2).or. &
                  abs(ivec(3,jj)-ivec(3,ii)).gt.nd2cell(3)) then
!
!  Find cell index - if outside user range then assign to central cell
!
                ncind13m = nd2central
                ncind13p = nd2central
              else
!
!  Compute index
!
                ncind13p = nd2cellptr(nd2cell(1)+1+ivec(1,jj)-ivec(1,ii), &
                                      nd2cell(2)+1+ivec(2,jj)-ivec(2,ii), &
                                      nd2cell(3)+1+ivec(3,jj)-ivec(3,ii))
                ncind13m = nd2cellptr(nd2cell(1)+1-ivec(1,jj)+ivec(1,ii), &
                                      nd2cell(2)+1-ivec(2,jj)+ivec(2,ii), &
                                      nd2cell(3)+1-ivec(3,jj)+ivec(3,ii))
              endif
!**********************************
!  Loop over last end site 4 / l  *
!**********************************
              l4loop: do l4 = lloopmin,lloop
                if (llbond) then
                  l = nuniquelptr(l4)
                else
                  l = l4
                endif
!
!  Check that one of the atoms is the matom
!
                if (i.ne.matom.and.j.ne.matom.and.k.ne.matom.and.l.ne.matom) cycle l4loop
!
                nl = nat(l)
                ntypl = nftype(l)
                ocl = occuf(l)
!
!  Check l is allowed for n
!             
                if (.not.lmatch(nl,ntypl,nt4,ntyp4,.true.)) cycle l4loop
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
!
!  Find cell index for 1-4
!
                  if (abs(ivec(1,ll)-ivec(1,ii)).gt.nd2cell(1).or. &
                      abs(ivec(2,ll)-ivec(2,ii)).gt.nd2cell(2).or. &
                      abs(ivec(3,ll)-ivec(3,ii)).gt.nd2cell(3)) then
!
!  Find cell index - if outside user range then assign to central cell
!
                    ncind14m = nd2central
                    ncind14p = nd2central
                  else
!
!  Compute index
!
                    ncind14p = nd2cellptr(nd2cell(1)+1-ivec(1,ll)+ivec(1,ii), &
                                          nd2cell(2)+1-ivec(2,ll)+ivec(2,ii), &
                                          nd2cell(3)+1-ivec(3,ll)+ivec(3,ii))
                    ncind14m = nd2cellptr(nd2cell(1)+1+ivec(1,ll)-ivec(1,ii), &
                                          nd2cell(2)+1+ivec(2,ll)-ivec(2,ii), &
                                          nd2cell(3)+1+ivec(3,ll)-ivec(3,ii))
                  endif
!
!  Find cell index for 2-4
!
                  if (abs(ivec(1,ll)).gt.nd2cell(1).or. &
                      abs(ivec(2,ll)).gt.nd2cell(2).or. &
                      abs(ivec(3,ll)).gt.nd2cell(3)) then
!
!  Find cell index - if outside user range then assign to central cell
!
                    ncind24m = nd2central
                    ncind24p = nd2central
                  else
!
!  Compute index
!
                    ncind24p = nd2cellptr(nd2cell(1)+1+ivec(1,ll), &
                                          nd2cell(2)+1+ivec(2,ll), &
                                          nd2cell(3)+1+ivec(3,ll))
                    ncind24m = nd2cellptr(nd2cell(1)+1-ivec(1,ll), &
                                          nd2cell(2)+1-ivec(2,ll), &
                                          nd2cell(3)+1-ivec(3,ll))
                  endif
!
!  Find cell index for 3-4
!
                  if (abs(ivec(1,ll)-ivec(1,jj)).gt.nd2cell(1).or. &
                      abs(ivec(2,ll)-ivec(2,jj)).gt.nd2cell(2).or. &
                      abs(ivec(3,ll)-ivec(3,jj)).gt.nd2cell(3)) then
!
!  Find cell index - if outside user range then assign to central cell
!
                    ncind34m = nd2central
                    ncind34p = nd2central
                  else
!
!  Compute index
!
                    ncind34p = nd2cellptr(nd2cell(1)+1+ivec(1,ll)-ivec(1,jj), &
                                          nd2cell(2)+1+ivec(2,ll)-ivec(2,jj), &
                                          nd2cell(3)+1+ivec(3,ll)-ivec(3,jj))
                    ncind34m = nd2cellptr(nd2cell(1)+1-ivec(1,ll)+ivec(1,jj), &
                                          nd2cell(2)+1-ivec(2,ll)+ivec(2,jj), &
                                          nd2cell(3)+1-ivec(3,ll)+ivec(3,jj))
                  endif
!*********************************
!  Valid four-body term located  *
!*********************************
!
!  Make a local copy of the vectors
!
                  x21d = x21
                  y21d = y21
                  z21d = z21
                  x31d = x31
                  y31d = y31
                  z31d = z31
                  x41d = x41
                  y41d = y41
                  z41d = z41
                  x42d = x42
                  y42d = y42
                  z42d = z42
                  x43d = x43
                  y43d = y43
                  z43d = z43
                  r212d = r212
                  r312d = r312
                  r412d = r412
                  r422d = r422
                  r432d = r432
!
!  Correct vectors for matom displacement
!
                  if (i.eq.matom) then
                    x21d = x21 - (vmatom(1) - xc1t)
                    y21d = y21 - (vmatom(2) - yc1t)
                    z21d = z21 - (vmatom(3) - zc1t)
                    x31d = x31 - (vmatom(1) - xc1t)
                    y31d = y31 - (vmatom(2) - yc1t)
                    z31d = z31 - (vmatom(3) - zc1t)
                    x41d = x41 - (vmatom(1) - xc1t)
                    y41d = y41 - (vmatom(2) - yc1t)
                    z41d = z41 - (vmatom(3) - zc1t)
                    r212d = x21d**2 + y21d**2 + z21d**2
                    r312d = x31d**2 + y31d**2 + z31d**2
                    r412d = x41d**2 + y41d**2 + z41d**2
                  endif
                  if (j.eq.matom) then
                    x21d = x21 + (vmatom(1) - xc2)
                    y21d = y21 + (vmatom(2) - yc2)
                    z21d = z21 + (vmatom(3) - zc2)
                    x42d = x42 - (vmatom(1) - xc2)
                    y42d = y42 - (vmatom(2) - yc2)
                    z42d = z42 - (vmatom(3) - zc2)
                    r212d = x21d**2 + y21d**2 + z21d**2
                    r422d = x42d**2 + y42d**2 + z42d**2
                  endif
                  if (k.eq.matom) then
                    x31d = x31 + (vmatom(1) - xc3t)
                    y31d = y31 + (vmatom(2) - yc3t)
                    z31d = z31 + (vmatom(3) - zc3t)
                    x43d = x43 - (vmatom(1) - xc3t)
                    y43d = y43 - (vmatom(2) - yc3t)
                    z43d = z43 - (vmatom(3) - zc3t)
                    r312d = x31d**2 + y31d**2 + z31d**2
                    r432d = x43d**2 + y43d**2 + z43d**2
                  endif
                  if (l.eq.matom) then
                    x41d = x41 + (vmatom(1) - xc4t)
                    y41d = y41 + (vmatom(2) - yc4t)
                    z41d = z41 + (vmatom(3) - zc4t)
                    x42d = x42 + (vmatom(1) - xc4t)
                    y42d = y42 + (vmatom(2) - yc4t)
                    z42d = z42 + (vmatom(3) - zc4t)
                    x43d = x43 + (vmatom(1) - xc4t)
                    y43d = y43 + (vmatom(2) - yc4t)
                    z43d = z43 + (vmatom(3) - zc4t)
                    r412d = x41d**2 + y41d**2 + z41d**2
                    r422d = x42d**2 + y42d**2 + z42d**2
                    r432d = x43d**2 + y43d**2 + z43d**2
                  endif
!
!  Finish calculating distances
!
                  r21d = sqrt(r212d)
                  r31d = sqrt(r312d)
                  r41d = sqrt(r412d)
                  r42d = sqrt(r422d)
                  r43d = sqrt(r432d)
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
                  nctor(1,1,niltor) = ncind12p
                  nctor(2,1,niltor) = ncind12m
                  nctor(1,2,niltor) = ncind13p
                  nctor(2,2,niltor) = ncind13m
                  nctor(1,3,niltor) = ncind14p
                  nctor(2,3,niltor) = ncind14m
                  nctor(1,4,niltor) = ncind23p
                  nctor(2,4,niltor) = ncind23m
                  nctor(1,5,niltor) = ncind24p
                  nctor(2,5,niltor) = ncind24m
                  nctor(1,6,niltor) = ncind34p
                  nctor(2,6,niltor) = ncind34m
!
                  riltor(1,niltor) = r21d
                  riltor(2,niltor) = r31d
                  riltor(3,niltor) = r41d
                  riltor(4,niltor) = r42d
                  riltor(5,niltor) = r43d
!
                  xiltor(1,niltor) = x21d
                  yiltor(1,niltor) = y21d
                  ziltor(1,niltor) = z21d
                  xiltor(2,niltor) = x31d
                  yiltor(2,niltor) = y31d
                  ziltor(2,niltor) = z31d
                  xiltor(3,niltor) = x41d
                  yiltor(3,niltor) = y41d
                  ziltor(3,niltor) = z41d
                  xiltor(4,niltor) = x42d
                  yiltor(4,niltor) = y42d
                  ziltor(4,niltor) = z42d
                  xiltor(5,niltor) = x43d
                  yiltor(5,niltor) = y43d
                  ziltor(5,niltor) = z43d
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
!
!  Make a local copy of the 3-2 vector
!
        x32d = x32
        y32d = y32
        z32d = z32
        r322d = r322
!
!  Correct 3-2 vector for matom displacement
!
        if (j.eq.matom) then
          x32d = x32 - (vmatom(1) - xc2)
          y32d = y32 - (vmatom(2) - yc2)
          z32d = z32 - (vmatom(3) - zc2)
          r322d = x32d**2 + y32d**2 + z32d**2
        endif
        if (k.eq.matom) then
          x32d = x32 + (vmatom(1) - xc3t)
          y32d = y32 + (vmatom(2) - yc3t)
          z32d = z32 + (vmatom(3) - zc3t)
          r322d = x32d**2 + y32d**2 + z32d**2
        endif
!
!  Finish calculating 3-2 distance
!
        r32d = sqrt(r322d)
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
          ncind12p = nctor(1,1,nil)
          ncind12m = nctor(2,1,nil)
          ncind13p = nctor(1,2,nil)
          ncind13m = nctor(2,2,nil)
          ncind14p = nctor(1,3,nil)
          ncind14m = nctor(2,3,nil)
          ncind23p = nctor(1,4,nil)
          ncind23m = nctor(2,4,nil)
          ncind24p = nctor(1,5,nil)
          ncind24m = nctor(2,5,nil)
          ncind34p = nctor(1,6,nil)
          ncind34m = nctor(2,6,nil)
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
          call fourbody(n,nfortype,r21,r31,r41,r32d,r42,r43,eterm,e1d,e2d,e3d, &
                        rko,rn,phi0o,isgn,fpoly,.true.,.false.,.false.)
!
          if (i.eq.matom) then
            d1cell(1,j,ncind12p) = d1cell(1,j,ncind12p) - x32d*e1d(4) + x21*e1d(1) - x42*e1d(5)
            d1cell(2,j,ncind12p) = d1cell(2,j,ncind12p) - y32d*e1d(4) + y21*e1d(1) - y42*e1d(5)
            d1cell(3,j,ncind12p) = d1cell(3,j,ncind12p) - z32d*e1d(4) + z21*e1d(1) - z42*e1d(5)
            d1cell(1,k,ncind13p) = d1cell(1,k,ncind13p) + x32d*e1d(4) - x43*e1d(6) + x31*e1d(2)
            d1cell(2,k,ncind13p) = d1cell(2,k,ncind13p) + y32d*e1d(4) - y43*e1d(6) + y31*e1d(2)
            d1cell(3,k,ncind13p) = d1cell(3,k,ncind13p) + z32d*e1d(4) - z43*e1d(6) + z31*e1d(2)
            d1cell(1,l,ncind14m) = d1cell(1,l,ncind14m) + x43*e1d(6) + x42*e1d(5) + x41*e1d(3)
            d1cell(2,l,ncind14m) = d1cell(2,l,ncind14m) + y43*e1d(6) + y42*e1d(5) + y41*e1d(3)
            d1cell(3,l,ncind14m) = d1cell(3,l,ncind14m) + z43*e1d(6) + z42*e1d(5) + z41*e1d(3)
          endif
          if (j.eq.matom) then
            d1cell(1,i,ncind12m) = d1cell(1,i,ncind12m) - x21*e1d(1) - x31*e1d(2) - x41*e1d(3)
            d1cell(2,i,ncind12m) = d1cell(2,i,ncind12m) - y21*e1d(1) - y31*e1d(2) - y41*e1d(3)
            d1cell(3,i,ncind12m) = d1cell(3,i,ncind12m) - z21*e1d(1) - z31*e1d(2) - z41*e1d(3)
            d1cell(1,k,ncind23m) = d1cell(1,k,ncind23m) + x32d*e1d(4) - x43*e1d(6) + x31*e1d(2)
            d1cell(2,k,ncind23m) = d1cell(2,k,ncind23m) + y32d*e1d(4) - y43*e1d(6) + y31*e1d(2)
            d1cell(3,k,ncind23m) = d1cell(3,k,ncind23m) + z32d*e1d(4) - z43*e1d(6) + z31*e1d(2)
            d1cell(1,l,ncind24p) = d1cell(1,l,ncind24p) + x43*e1d(6) + x42*e1d(5) + x41*e1d(3)
            d1cell(2,l,ncind24p) = d1cell(2,l,ncind24p) + y43*e1d(6) + y42*e1d(5) + y41*e1d(3)
            d1cell(3,l,ncind24p) = d1cell(3,l,ncind24p) + z43*e1d(6) + z42*e1d(5) + z41*e1d(3)
          endif
          if (k.eq.matom) then
            d1cell(1,i,ncind13m) = d1cell(1,i,ncind13m) - x21*e1d(1) - x31*e1d(2) - x41*e1d(3)
            d1cell(2,i,ncind13m) = d1cell(2,i,ncind13m) - y21*e1d(1) - y31*e1d(2) - y41*e1d(3)
            d1cell(3,i,ncind13m) = d1cell(3,i,ncind13m) - z21*e1d(1) - z31*e1d(2) - z41*e1d(3)
            d1cell(1,j,ncind23p) = d1cell(1,j,ncind23p) - x32d*e1d(4) + x21*e1d(1) - x42*e1d(5)
            d1cell(2,j,ncind23p) = d1cell(2,j,ncind23p) - y32d*e1d(4) + y21*e1d(1) - y42*e1d(5)
            d1cell(3,j,ncind23p) = d1cell(3,j,ncind23p) - z32d*e1d(4) + z21*e1d(1) - z42*e1d(5)
            d1cell(1,l,ncind34p) = d1cell(1,l,ncind34p) + x43*e1d(6) + x42*e1d(5) + x41*e1d(3)
            d1cell(2,l,ncind34p) = d1cell(2,l,ncind34p) + y43*e1d(6) + y42*e1d(5) + y41*e1d(3)
            d1cell(3,l,ncind34p) = d1cell(3,l,ncind34p) + z43*e1d(6) + z42*e1d(5) + z41*e1d(3)
          endif
          if (l.eq.matom) then
            d1cell(1,i,ncind14p) = d1cell(1,i,ncind14p) - x21*e1d(1) - x31*e1d(2) - x41*e1d(3)
            d1cell(2,i,ncind14p) = d1cell(2,i,ncind14p) - y21*e1d(1) - y31*e1d(2) - y41*e1d(3)
            d1cell(3,i,ncind14p) = d1cell(3,i,ncind14p) - z21*e1d(1) - z31*e1d(2) - z41*e1d(3)
            d1cell(1,j,ncind24m) = d1cell(1,j,ncind24m) - x32d*e1d(4) + x21*e1d(1) - x42*e1d(5)
            d1cell(2,j,ncind24m) = d1cell(2,j,ncind24m) - y32d*e1d(4) + y21*e1d(1) - y42*e1d(5)
            d1cell(3,j,ncind24m) = d1cell(3,j,ncind24m) - z32d*e1d(4) + z21*e1d(1) - z42*e1d(5)
            d1cell(1,k,ncind34m) = d1cell(1,k,ncind34m) + x32d*e1d(4) - x43*e1d(6) + x31*e1d(2)
            d1cell(2,k,ncind34m) = d1cell(2,k,ncind34m) + y32d*e1d(4) - y43*e1d(6) + y31*e1d(2)
            d1cell(3,k,ncind34m) = d1cell(3,k,ncind34m) + z32d*e1d(4) - z43*e1d(6) + z31*e1d(2)
          endif
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
  if (noofp.gt.0) call fouroop1fc(maxlhs,d1cell,matom,vmatom)
!************************
!  Improper potentials  *
!************************
  if (nimproper.gt.0) call fourimp1fc(maxlhs,d1cell,matom,vmatom)
!
!  Free local memory
!
  deallocate(zvec,stat=status)
  if (status/=0) call deallocate_error('four1fc','zvec')
  deallocate(yvec,stat=status)
  if (status/=0) call deallocate_error('four1fc','yvec')
  deallocate(xvec,stat=status)
  if (status/=0) call deallocate_error('four1fc','xvec')
  deallocate(ivec,stat=status)
  if (status/=0) call deallocate_error('four1fc','ivec')
  deallocate(nuniquelptr,stat=status)
  if (status/=0) call deallocate_error('four1fc','nuniquelptr')
  deallocate(nuniquekptr,stat=status)
  if (status/=0) call deallocate_error('four1fc','nuniquekptr')
  deallocate(nuniqueiptr,stat=status)
  if (status/=0) call deallocate_error('four1fc','nuniqueiptr')
  deallocate(nptrnfornonoop,stat=status)
  if (status/=0) call deallocate_error('four1fc','nptrnfornonoop')
  deallocate(ntypmiddle,stat=status)
  if (status/=0) call deallocate_error('four1fc','ntypmiddle')
  deallocate(natmiddle,stat=status)
  if (status/=0) call deallocate_error('four1fc','natmiddle')
!
!  Timing
!
  time2 = g_cpu_time()
  tfour = tfour + time2 - time1
#ifdef TRACE
  call trace_out('four1fc')
#endif
!
  return
  end
