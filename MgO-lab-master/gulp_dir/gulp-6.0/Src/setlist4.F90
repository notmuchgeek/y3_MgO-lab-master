  subroutine setlist4
!
!  Generate list of valid four body terms for MD
!
!  nlist4md = number of valid four body terms
!  nforptr  = pointer to potential number
!  ilind    = pointer to terminal atoms
!  jkind    = pointer to middle atoms
!  icell41  = pointer to cell vector
!  icell42  = pointer to cell vector
!  icell43  = pointer to cell vector
!
!  Strategy - sift by potential first, then cutoffs
!
!   1/95 Intra/inter-molecular specification added
!   2/95 Bonded specification added for four-body terms
!   3/95 Periodic molecule corrections added
!   4/97 Modifications added to allow for out of plane pots
!   2/01 lintoijk calls now use imaxl,jmaxl and kmaxl
!   4/01 Checking altered so that bonding takes precedence
!        over distances
!   6/01 Setting of lsamemol altered
!   9/01 lmolq calculations accelerated using lneedmol
!  11/02 Wildcard atoms added
!  12/02 Cell indices corrected for 1-D and 2-D cases
!   5/04 Search over lattice vectors generalised
!  10/05 Inversion potential added
!   6/06 Inversion squared potential added
!   9/06 Order of atom search changed to allow for Dreiding
!   9/06 Dreiding scheme for force constant added as an option
!  10/06 imaxl, jmaxl & kmaxl renamed to avoid conflict with global symbols
!  10/06 Setting of icell41/2/3 now uses imax/jmax/kmax rather than fixed values
!   1/07 Wildcard handling in lmatch calls corrected
!   2/07 Bonding types and test added
!   4/07 Code reordered so that atom loops are on the outside and potentials on the inside
!   4/07 Screening of j & k for potential middle species added
!   4/07 Bond type checking extended to ndim > 0
!   5/07 QM/MM schemes added
!   6/07 Bug due to niltor not being set fixed
!   6/07 lmolok reset to initial value for each potential loop
!  10/07 Angle-angle cross potential added
!  10/07 Error in checking of exocyclic attribute for bonds corrected
!  10/07 Handling of n4botype(2,n) modified
!  12/07 Unused variables removed
!   5/08 UFFoop potential added
!  11/08 loutofplane array introduced to indicate whether a potential is out-of-plane 
!        type or not to replace potential number checking
!  11/08 Setting of maximum cutoffs now only looks at non-out of plane potentials
!  11/08 New logic for matching species introduced to handle case of the same element
!        being the middle atom, but one with a specific type and the other without.
!  11/08 ixl etc defined relative to ixj etc rather than ixk to correct error
!   5/10 Code modified to increase speed for bonded potentials by only searching over
!        bonded atoms
!   7/13 Distance check logic changed for bonded potentials
!   7/13 Improper torsion type added
!   9/13 Separate looping added for all bonded potential case to get speed up
!   4/19 Call to lmatchpair modified to include wildcard argument
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!   9/19 lnotorXduplicates added
!  12/19 nlist4md0 zeroed to avoid uninitialised use depending on algorithm
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
!  Julian Gale, CIC, Curtin University, December 2019
!
  use configurations, only : nregionno, nregiontype, QMMMmode
  use control,        only : lnotorXduplicates
  use current
  use element,        only : maxele
  use four
  use molecule
  use parallel
  use times
  implicit none
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ia
  integer(i4)                                  :: ii
  integer(i4)                                  :: iloop
  integer(i4)                                  :: iloopmin
  integer(i4)                                  :: imax
  integer(i4)                                  :: indmi
  integer(i4)                                  :: indmj
  integer(i4)                                  :: indmk
  integer(i4)                                  :: indml
  integer(i4)                                  :: ind
  integer(i4)                                  :: ix
  integer(i4)                                  :: iy
  integer(i4)                                  :: iz
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
  integer(i4)                                  :: jmax
  integer(i4)                                  :: jxx
  integer(i4)                                  :: jyy
  integer(i4)                                  :: jzz
  integer(i4)                                  :: k
  integer(i4)                                  :: kloop
  integer(i4)                                  :: kmax
  integer(i4)                                  :: kxx
  integer(i4)                                  :: kyy
  integer(i4)                                  :: kzz
  integer(i4)                                  :: l
  integer(i4)                                  :: l4
  integer(i4)                                  :: la
  integer(i4)                                  :: li
  integer(i4)                                  :: lk
  integer(i4)                                  :: ll
  integer(i4)                                  :: lloop
  integer(i4)                                  :: lloopmin
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
  integer(i4)                                  :: nfornonoop
  integer(i4)                                  :: ni
  integer(i4)                                  :: nil
  integer(i4)                                  :: nil2
  integer(i4)                                  :: niltor
  integer(i4)                                  :: nimproper
  integer(i4)                                  :: nj
  integer(i4)                                  :: nk
  integer(i4)                                  :: nl 
  integer(i4)                                  :: nlist4md0
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
  integer(i4), dimension(:), allocatable       :: nptrnfornonoop
  integer(i4)                                  :: nregioni
  integer(i4)                                  :: nregionj
  integer(i4)                                  :: nregionk
  integer(i4)                                  :: nregionl
  integer(i4)                                  :: nregiontypi
  integer(i4)                                  :: nregiontypj
  integer(i4)                                  :: nregiontypk
  integer(i4)                                  :: nregiontypl
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
  integer(i4), dimension(:),   allocatable     :: nuniqueiind
  integer(i4), dimension(:),   allocatable     :: nuniquekind
  integer(i4), dimension(:),   allocatable     :: nuniquelind
  integer(i4), dimension(:),   allocatable     :: nuniqueiptr
  integer(i4), dimension(:),   allocatable     :: nuniquekptr
  integer(i4), dimension(:),   allocatable     :: nuniquelptr
  integer(i4), dimension(:,:), allocatable     :: nuniquektyp
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
  logical                                      :: libond
  logical                                      :: liok
  logical                                      :: limatch1
  logical                                      :: limatch4
  logical                                      :: limproper
  logical                                      :: lintra_only
  logical                                      :: linter_only
  logical                                      :: ljmatch2
  logical                                      :: ljmatch3
  logical                                      :: ljkmatch
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
  logical                                      :: ltsyme_exact
  logical                                      :: lunique
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: cut
  real(dp)                                     :: cutmax
  real(dp)                                     :: r212
  real(dp)                                     :: r312
  real(dp)                                     :: r322
  real(dp)                                     :: r412
  real(dp)                                     :: r422
  real(dp)                                     :: r432
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
  real(dp), dimension(:), allocatable          :: xvec
  real(dp), dimension(:), allocatable          :: yvec
  real(dp), dimension(:), allocatable          :: zvec
!
  time1 = g_cpu_time()
  nlist4md = 0
  nlist4md0 = 0
  lmolloc = (nmol.gt.0)
!
!  Allocate local memory
!
  allocate(natmiddle(2*nfor),stat=status)
  if (status/=0) call outofmemory('setlist4','natmiddle')
  allocate(ntypmiddle(2*nfor),stat=status)
  if (status/=0) call outofmemory('setlist4','ntypmiddle')
  allocate(nptrnfornonoop(nfor),stat=status)
  if (status/=0) call outofmemory('setlist4','nptrnfornonoop')
  allocate(nuniqueiind(numat),stat=status)
  if (status/=0) call outofmemory('setlist4','nuniqueiind')
  allocate(nuniquekind(numat),stat=status)
  if (status/=0) call outofmemory('setlist4','nuniquekind')
  allocate(nuniquelind(numat),stat=status)
  if (status/=0) call outofmemory('setlist4','nuniquelind')
  allocate(nuniqueiptr(numat),stat=status)
  if (status/=0) call outofmemory('setlist4','nuniqueiptr')
  allocate(nuniquekptr(numat),stat=status)
  if (status/=0) call outofmemory('setlist4','nuniquekptr')
  allocate(nuniquelptr(numat),stat=status)
  if (status/=0) call outofmemory('setlist4','nuniquelptr')
  allocate(nuniquektyp(2,numat),stat=status)
  if (status/=0) call outofmemory('setlist4','nuniquektyp')
  if (ndim.gt.0) then
    allocate(xvec(maxvector),stat=status)
    if (status/=0) call outofmemory('setlist4','xvec')
    allocate(yvec(maxvector),stat=status)
    if (status/=0) call outofmemory('setlist4','yvec')
    allocate(zvec(maxvector),stat=status)
    if (status/=0) call outofmemory('setlist4','zvec')
  else
    allocate(xvec(1),stat=status)
    if (status/=0) call outofmemory('setlist4','xvec')
    allocate(yvec(1),stat=status)
    if (status/=0) call outofmemory('setlist4','yvec')
    allocate(zvec(1),stat=status)
    if (status/=0) call outofmemory('setlist4','zvec')
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
    if (.not.loutofplane(n).and.mmfexc(n).ne.2) then
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
    if (.not.loutofplane(n).and.mmfexc(n).ne.2) then
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
    if (.not.loutofplane(n).and.mmfexc(n).ne.2) then
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
      if (status/=0) call deallocate_error('setlist4','zvec')
      deallocate(yvec,stat=status)
      if (status/=0) call deallocate_error('setlist4','yvec')
      deallocate(xvec,stat=status)
      if (status/=0) call deallocate_error('setlist4','xvec')
      maxvector = nint(1.1*nvector)
      allocate(xvec(maxvector),stat=status)
      if (status/=0) call outofmemory('setlist4','xvec')
      allocate(yvec(maxvector),stat=status)
      if (status/=0) call outofmemory('setlist4','yvec')
      allocate(zvec(maxvector),stat=status)
      if (status/=0) call outofmemory('setlist4','zvec')
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
  if (lallbtyp) then
!******************************************************************************************
!  If all torsions are of bonded type then we can use pure connectivity to find torsions  *
!******************************************************************************************
!
!  Loop over middle site 2 / j
!
    ljbloop: do j = 1,numat
!
!  Skip if there are no bonds for this atom
!
      if (nbonds(j).eq.0) cycle ljbloop
!
      nj = nat(j)
      ntypj = nftype(j)
      nregionj = nregionno(nsft+nrelf2a(j))
      nregiontypj = nregiontype(nregionj,ncf)
!
!  Check whether species may be valid
!
      if (.not.lmatchany(nj,ntypj,nmiddle,natmiddle,ntypmiddle,.true.)) cycle ljbloop
!
!  Set loop range for k
!
      nuniquek = 1
      nuniquekptr(1) = nbonded(1,j)
      nuniquekind(1) = nbondind(1,j)
      nuniquektyp(1,1) = nbondedtype(1,1,j)
      nuniquektyp(2,1) = nbondedtype(2,1,j)
      do kloop = 2,nbonds(j)
        lunique = .true.
        do lu = 1,nuniquek
          if (nbonded(kloop,j).eq.nuniquekptr(lu)) lunique = .false.
        enddo
        if (lunique) then
          nuniquek = nuniquek + 1
          nuniquekptr(nuniquek) = nbonded(kloop,j)
          nuniquekind(nuniquek) = nbondind(kloop,j)
          nuniquektyp(1,nuniquek) = nbondedtype(1,kloop,j)
          nuniquektyp(2,nuniquek) = nbondedtype(2,kloop,j)
        endif
      enddo
!
!  Loop over second middle site 3 / k
!
      lkbloop: do lk = 1,nuniquek
        k = nuniquekptr(lk)
!
!  Only do upper triangular set of atoms
!
        if (k.lt.j) cycle lkbloop
!
        nk = nat(k)
        ntypk = nftype(k)
        nregionk = nregionno(nsft+nrelf2a(k))
        nregiontypk = nregiontype(nregionk,ncf)
!
!  Check whether species may be valid
!
        if (.not.lmatchany(nk,ntypk,nmiddle,natmiddle,ntypmiddle,.true.)) cycle lkbloop
        if (.not.lmatchpair(nj,ntypj,nk,ntypk,nfor,nfspec2,nfptyp2,nfspec3,nfptyp3,.true.)) cycle lkbloop
!
!  Set counter for number of valid i/l end atom combinations
!
        niltor = 0
        nlist4md0 = nlist4md
!
!  Loop over four-body potentials
!
        bpots: do nn = 1,nfornonoop
          n = nptrnfornonoop(nn)
          nt1 = nfspec1(n)
          nt2 = nfspec2(n)
          nt3 = nfspec3(n)
          nt4 = nfspec4(n)
          ntyp1 = nfptyp1(n)
          ntyp2 = nfptyp2(n)
          ntyp3 = nfptyp3(n)
          ntyp4 = nfptyp4(n)
          ltsyme_exact = lexactmatch(nt1,ntyp1,nt4,ntyp4)
          limproper = (mmfexc(n).eq.2)
!
!  QM/MM handling : j & k are both QM atoms and potential is of bonded type => exclude
!
          if (QMMMmode(ncf).gt.0) then
            if (nregiontypj.eq.1.and.nregiontypk.eq.1) cycle bpots
          endif
!
!  Validate potential for j-k pair
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
          if (.not.ljkmatch.and..not.lkjmatch) cycle bpots
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
            if (.not.ltsyme_exact) then
              ntmp = nt1
              nt1 = nt4
              nt4 = ntmp
              ntmp = ntyp1
              ntyp1 = ntyp4
              ntyp4 = ntmp
            endif
          endif
!
!  Set flag indicating whether middle atoms could be matched either way round
!
          lmeither = (ljkmatch.and.lkjmatch)
!               
!  Check central bond type for correct order
!
          if (n4botype(1,n).gt.0) then
            if (n4botype(1,n).ne.nuniquektyp(1,lk)) cycle bpots
          endif
          if (n4botype(2,n).ne.nuniquektyp(2,lk)) cycle bpots
!
!  Set loop range for i
!
          if (limproper.and.lkjmatch) then
            if (nbonds(k).gt.0) then
              nuniquei = 1
              nuniqueiptr(1) = nbonded(1,k)
              nuniqueiind(1) = nbondind(1,k)
              do iloop = 2,nbonds(k)
                lunique = .true.
                do lu = 1,nuniquei
                  if (nbonded(iloop,k).eq.nuniqueiptr(lu)) lunique = .false.
                enddo
                if (lunique) then
                  nuniquei = nuniquei + 1
                  nuniqueiptr(nuniquei) = nbonded(iloop,k)
                  nuniqueiind(nuniquei) = nbondind(iloop,k)
                endif
              enddo
              iloop = nuniquei
            else
              iloop = 0
            endif
          else
            if (nbonds(j).gt.0) then
              nuniquei = 1
              nuniqueiptr(1) = nbonded(1,j)
              nuniqueiind(1) = nbondind(1,j)
              do iloop = 2,nbonds(j)
                lunique = .true.
                do lu = 1,nuniquei
                  if (nbonded(iloop,j).eq.nuniqueiptr(lu)) lunique = .false.
                enddo
                if (lunique) then
                  nuniquei = nuniquei + 1
                  nuniqueiptr(nuniquei) = nbonded(iloop,j)
                  nuniqueiind(nuniquei) = nbondind(iloop,j)
                endif
              enddo
              iloop = nuniquei
            else
              iloop = 0
            endif
          endif
!
!  Skip if iloop is zero
!
          if (iloop.eq.0) cycle bpots
!
!  Loop over end site 1 / i
!
          libloop: do li = 1,iloop
            i = nuniqueiptr(li)
            ni = nat(i)
            ntypi = nftype(i)
            nregioni = nregionno(nsft+nrelf2a(i))
            nregiontypi = nregiontype(nregioni,ncf)
!
!  Check whether i matches either of types 1 and 4
!
            limatch1 = lmatch(ni,ntypi,nt1,ntyp1,.true.)
            limatch4 = lmatch(ni,ntypi,nt4,ntyp4,.true.)
!
!  Is i allowed for type 1, or type 4 if the middle atoms can be switched?
!
            liok = (limatch1.or.(limatch4.and.lmeither))
            if (.not.liok) cycle libloop
!
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
            endif
!
!  Prevent atoms i and k being the same atom
!
            if (i.eq.k) cycle libloop
!
!  Set loop range for l
!
            if (limproper.and.ljkmatch) then
              if (nbonds(j).gt.0) then
                nuniquel = 1
                nuniquelptr(1) = nbonded(1,j)
                nuniquelind(1) = nbondind(1,j)
                do lloop = 2,nbonds(j)
                  lunique = .true.
                  do lu = 1,nuniquel
                    if (nbonded(lloop,j).eq.nuniquelptr(lu)) lunique = .false.
                  enddo
                  if (lunique) then
                    nuniquel = nuniquel + 1
                    nuniquelptr(nuniquel) = nbonded(lloop,j)
                    nuniquelind(nuniquel) = nbondind(lloop,j)
                  endif
                enddo
                lloop = nuniquel
              else
                lloop = 0
              endif
            else
              if (nbonds(k).gt.0) then
                nuniquel = 1
                nuniquelptr(1) = nbonded(1,k)
                nuniquelind(1) = nbondind(1,k)
                do lloop = 2,nbonds(k)
                  lunique = .true.
                  do lu = 1,nuniquel
                    if (nbonded(lloop,k).eq.nuniquelptr(lu)) lunique = .false.
                  enddo
                  if (lunique) then
                    nuniquel = nuniquel + 1
                    nuniquelptr(nuniquel) = nbonded(lloop,k)
                    nuniquelind(nuniquel) = nbondind(lloop,k)
                  endif
                enddo
                lloop = nuniquel
              else
                lloop = 0
              endif
            endif
!
!  Skip if lloop is zero
!
            if (lloop.eq.0) cycle libloop
!
!  Loop over last end site 4 / l
!
            l4bloop: do l4 = 1,lloop
              l = nuniquelptr(l4)
              nl = nat(l)
              ntypl = nftype(l)
              nregionl = nregionno(nsft+nrelf2a(l))
              nregiontypl = nregiontype(nregionl,ncf)
!
!  Prevent atoms j and l being the same atom
!
              if (j.eq.l) cycle l4bloop
!
!  Check l is allowed for n
!             
              if (.not.lmatch(nl,ntypl,nt4,ntyp4,.true.)) cycle l4bloop
!
!  QM/MM handling : i, j, k & l are all QM atoms => exclude
!               
              if (QMMMmode(ncf).gt.0) then
                if (nregiontypi.eq.1.and.nregiontypj.eq.1.and.nregiontypk.eq.1.and.nregiontypl.eq.1) cycle l4bloop
              endif
!
!  Valid four-body term located
!
              nlist4md = nlist4md + 1
              niltor   = niltor   + 1
              if (nlist4md.gt.maxlist4) then
                maxlist4 = nlist4md + 100
                call changemaxlist4
              endif
              ilind(nlist4md) = i + l*(numat+1)
              jkind(nlist4md) = j + k*(numat+1)
              nforptr(nlist4md) = n
              if (ndim.gt.0) then
!
!  Shift icell43 using icell42
!
                ii = nuniquekind(lk)
                jzz = (ii/100) - 5
                ii = ii - 100*(jzz+5)
                jyy = (ii/10) - 5
                ii = ii - 10*(jyy+5)
                jxx = ii - 5
!
                ii = nuniquelind(l4)
                kzz = (ii/100) - 5
                ii = ii - 100*(kzz+5)
                kyy = (ii/10) - 5
                ii = ii - 10*(kyy+5)
                kxx = ii - 5
!
                icell41(nlist4md) = nuniqueiind(li)
                icell42(nlist4md) = nuniquekind(lk)
                icell43(nlist4md) = (kxx+jxx) + 5 + 10*(kyy+jyy+5) + 100*(kzz+jzz+5)
              else
                icell41(nlist4md) = 0
                icell42(nlist4md) = 0 
                icell43(nlist4md) = 0
              endif
!
!  End of inner loops over atoms and cell vectors
!
            enddo l4bloop
          enddo libloop
!
!  End loop over potentials
!
        enddo bpots
!*******************************
!  Remove wildcard duplicates  *
!*******************************
        if (lnotorXduplicates) then
          lkeeptor(1:niltor) = .true.
          do nil = 1,niltor-1
            if (lkeeptor(nil)) then
              ind = ilind(nlist4md0+nil)
              l = ind/(numat+1)
              i = ind - l*(numat+1)
              n = nforptr(nlist4md0+nil)
              nwild1 = 0
              if (nfspec1(n).eq.maxele) nwild1 = nwild1 + 1
              if (nfspec4(n).eq.maxele) nwild1 = nwild1 + 1
!
!  Loop over remaining torsions looking for matches
!
              nmatch = 0
              do nil2 = nil+1,niltor
                ind = ilind(nlist4md0+nil2)
                la = ind/(numat+1)
                ia = ind - la*(numat+1)
                if (ia.eq.i.and.la.eq.l) then
!
!  Add duplicate potential to the list
!
                  nmatch = nmatch + 1
                  n2 = nfortor(nlist4md0+nil2)
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
!*********************************************************
!  Condense torsion list for those that need to be kept  *
!*********************************************************
          if (nkeepfor.lt.niltor) then
            do nil = 1,nkeepfor
              nil2 = nkeeptor(nil)
              ilind(nlist4md0+nil) = ilind(nlist4md0+nil2)
              jkind(nlist4md0+nil) = jkind(nlist4md0+nil2)
              nforptr(nlist4md0+nil) = nforptr(nlist4md0+nil2)
              icell41(nlist4md0+nil) = icell41(nlist4md0+nil2)
              icell42(nlist4md0+nil) = icell42(nlist4md0+nil2)
              icell43(nlist4md0+nil) = icell43(nlist4md0+nil2)
            enddo
          endif
!
!  Assign number of i/l combinations to preceeding potentials
!
          do l = 1,nkeepfor
            ilnum(nlist4md+1-l) = nkeepfor
          enddo
        else
!
!  Assign number of i/l combinations to preceeding potentials
!
          do l = 1,niltor
            ilnum(nlist4md+1-l) = niltor
          enddo
        endif
      enddo lkbloop
    enddo ljbloop
  else
!****************************************************
!  List of torsions contains non-bonded potentials  *
!****************************************************
!
!  Loop over middle site 2 / j
!
    ljloop:  do j = 1,numat
      nj = nat(j)
      ntypj = nftype(j)
      nregionj = nregionno(nsft+nrelf2a(j))
      nregiontypj = nregiontype(nregionj,ncf)
!
!  Check whether species may be valid
!
      if (.not.lmatchany(nj,ntypj,nmiddle,natmiddle,ntypmiddle,.true.)) cycle ljloop
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
!
!  Loop over second middle site 3 / k
!
      lkloop: do k = j,numat
        nk = nat(k)
        ntypk = nftype(k)
        nregionk = nregionno(nsft+nrelf2a(k))
        nregiontypk = nregiontype(nregionk,ncf)
!
!  Check whether species may be valid
!
        if (.not.lmatchany(nk,ntypk,nmiddle,natmiddle,ntypmiddle,.true.)) cycle lkloop
        if (.not.lmatchpair(nj,ntypj,nk,ntypk,nfor,nfspec2,nfptyp2,nfspec3,nfptyp3,.true.)) cycle lkloop
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
        jjloop: do jj = 1,nvector
          r322 = (xvec(jj)+x32t)**2 + (yvec(jj)+y32t)**2 + (zvec(jj)+z32t)**2
          if (r322.lt.1d-12) cycle jjloop
          if (k.eq.j.and.jj.eq.nmid) cycle jjloop
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
          if (r322.gt.tr2max.and.(.not.lanybtyp.or..not.lbondedjk)) cycle jjloop
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
!
!  Loop over four-body potentials
!
          pots: do nn = 1,nfornonoop
            n = nptrnfornonoop(nn)
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
            limproper = (mmfexc(n).eq.2)
            lintra_only = (lfintra(n).and..not.lfinter(n))
            linter_only = (lfinter(n).and..not.lfintra(n))
            lneedmol = (lintra_only.or.linter_only.or.lbtyp)
!
!  Reset lmolok to initial state for j-k pair for each potential
!
            lmolok = lmolokjk
!
!  QM/MM handling : j & k are both QM atoms and potential is of bonded type => exclude
!
            if (QMMMmode(ncf).gt.0) then
              if (nregiontypj.eq.1.and.nregiontypk.eq.1.and.lbtyp) cycle pots
            endif
!
!  Validate potential for j-k pair
!
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
            if (limproper.and.lkjmatch) then
              libond = .true.
              if (nbonds(k).gt.0) then
                nuniquei = 1
                nuniqueiptr(1) = nbonded(1,k)
                do iloop = 2,nbonds(k)
                  lunique = .true.
                  do lu = 1,nuniquei
                    if (nbonded(iloop,k).eq.nuniqueiptr(lu)) lunique = .false.
                  enddo
                  if (lunique) then
                    nuniquei = nuniquei + 1
                    nuniqueiptr(nuniquei) = nbonded(iloop,k)
                  endif
                enddo
                iloop = nuniquei
              else
                iloop = 0
              endif
              iloopmin = 1
            elseif (lbtyp) then
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
!
!  Loop over end site 1 / i
!
            liloop: do li = iloopmin,iloop
              if (libond) then
                i = nuniqueiptr(li)
              else
                i = li
              endif
              ni = nat(i)
              ntypi = nftype(i)
              nregioni = nregionno(nsft+nrelf2a(i))
              nregiontypi = nregiontype(nregioni,ncf)
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
                    if (limproper.and.lkjmatch) then
                      call bonded(lbondedij,l2bondsij,nbtypeji,nbtypeji2,k,i,0_i4,0_i4,0_i4)
                      if (.not.lbondedij) cycle iiloop
                    elseif (lbtyp) then
                      call bonded(lbondedij,l2bondsij,nbtypeji,nbtypeji2,j,i,0_i4,0_i4,0_i4)
                      if (.not.lbondedij) cycle iiloop
                    endif
                  else
                    call lintoijk(ixx,iyy,izz,ii,imax,jmax,kmax)
                    if (limproper.and.lkjmatch) then
                      call bonded(lbondedij,l2bondsij,nbtypeji,nbtypeji2,k,i,ixx-jxx,iyy-jzz,izz-jzz)
                      if (.not.lbondedij) cycle iiloop
                      lsamemolij = (lbondedij.or.l2bondsij)
                    elseif (lbtyp) then
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
                if (limproper.and.ljkmatch) then
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
                  lloopmin = 1
                elseif (lbtyp) then
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
!  Loop over last end site 4 / l
!
                l4loop: do l4 = lloopmin,lloop
                  if (llbond) then
                    l = nuniquelptr(l4)
                  else
                    l = l4
                  endif
                  nl = nat(l)
                  ntypl = nftype(l)
                  nregionl = nregionno(nsft+nrelf2a(l))
                  nregiontypl = nregiontype(nregionl,ncf)
!
!  Check l is allowed for n
!             
                  if (.not.lmatch(nl,ntypl,nt4,ntyp4,.true.)) cycle l4loop
!
!  QM/MM handling : i, j, k & l are all QM atoms => exclude
!               
                  if (QMMMmode(ncf).gt.0) then
                    if (nregiontypi.eq.1.and.nregiontypj.eq.1.and.nregiontypk.eq.1.and.nregiontypl.eq.1) cycle l4loop
                  endif
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
                        if (limproper.and.ljkmatch) then
                          call bonded(lbondedkl,l2bondskl,nbtypekl,nbtypekl2,j,l,0_i4,0_i4,0_i4)
                          if (.not.lbondedkl) cycle llloop
                        elseif (lbtyp) then
                          call bonded(lbondedkl,l2bondskl,nbtypekl,nbtypekl2,k,l,0_i4,0_i4,0_i4)
                          if (.not.lbondedkl) cycle llloop
                        endif
                      else
                        call lintoijk(kxx,kyy,kzz,ll,imax,jmax,kmax)
                        if (limproper.and.ljkmatch) then
                          call bonded(lbondedkl,l2bondskl,nbtypekl,nbtypekl2,j,l,kxx,kyy,kzz)
                          if (.not.lbondedkl) cycle llloop
                          lsamemolkl = (lbondedkl.or.l2bondskl)
                        elseif (lbtyp) then
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
!  Valid four-body term located
!
                    nlist4md = nlist4md + 1
                    niltor   = niltor   + 1
                    if (nlist4md.gt.maxlist4) then
                      maxlist4 = nlist4md + 100
                      call changemaxlist4
                    endif
                    ilind(nlist4md) = i + l*(numat+1)
                    jkind(nlist4md) = j + k*(numat+1)
                    nforptr(nlist4md) = n
                    if (ndim.gt.0) then
                      ind = ii - 1
                      ix = (ind/((2*jmax + 1)*(2*kmax + 1)))
                      ind = ind - ix*(2*jmax + 1)*(2*kmax + 1)
                      iy = (ind/(2*kmax + 1))
                      ind = ind - iy*(2*kmax + 1)
                      iz = ind - kmax
                      iy = iy - jmax
                      ix = ix - imax
                      icell41(nlist4md) = ix + 5 + 10*(iy+5) + 100*(iz+5)
                      ind = jj - 1
                      ix = (ind/((2*jmax + 1)*(2*kmax + 1)))
                      ind = ind - ix*(2*jmax + 1)*(2*kmax + 1)
                      iy = (ind/(2*kmax + 1))
                      ind = ind - iy*(2*kmax + 1)
                      iz = ind - kmax
                      iy = iy - jmax
                      ix = ix - imax
                      icell42(nlist4md) = ix + 5 + 10*(iy+5) + 100*(iz+5)
                      ind = ll - 1
                      ix = (ind/((2*jmax + 1)*(2*kmax + 1)))
                      ind = ind - ix*(2*jmax + 1)*(2*kmax + 1)
                      iy = (ind/(2*kmax + 1))
                      ind = ind - iy*(2*kmax + 1)
                      iz = ind - kmax
                      iy = iy - jmax
                      ix = ix - imax
                      icell43(nlist4md) = ix + 5 + 10*(iy+5) + 100*(iz+5)
                    else
                      icell41(nlist4md) = 0
                      icell42(nlist4md) = 0 
                      icell43(nlist4md) = 0
                    endif
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
                ind = ilind(nlist4md0+nil)
                l = ind/(numat+1)
                i = ind - l*(numat+1)
                n = nforptr(nlist4md0+nil)
                nwild1 = 0
                if (nfspec1(n).eq.maxele) nwild1 = nwild1 + 1
                if (nfspec4(n).eq.maxele) nwild1 = nwild1 + 1
!
!  Loop over remaining torsions looking for matches
!
                nmatch = 0
                do nil2 = nil+1,niltor
                  ind = ilind(nlist4md0+nil2)
                  la = ind/(numat+1)
                  ia = ind - la*(numat+1)
                  if (ia.eq.i.and.la.eq.l) then
!
!  Add duplicate potential to the list
!
                    nmatch = nmatch + 1
                    n2 = nfortor(nlist4md0+nil2)
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
!*********************************************************
!  Condense torsion list for those that need to be kept  *
!*********************************************************
            if (nkeepfor.lt.niltor) then
              do nil = 1,nkeepfor
                nil2 = nkeeptor(nil)
                ilind(nlist4md0+nil) = ilind(nlist4md0+nil2)
                jkind(nlist4md0+nil) = jkind(nlist4md0+nil2)
                nforptr(nlist4md0+nil) = nforptr(nlist4md0+nil2)
                icell41(nlist4md0+nil) = icell41(nlist4md0+nil2)
                icell42(nlist4md0+nil) = icell42(nlist4md0+nil2)
                icell43(nlist4md0+nil) = icell43(nlist4md0+nil2)
              enddo
            endif
!
!  Assign number of i/l combinations to preceeding potentials
!
            do l = 1,nkeepfor
              ilnum(nlist4md+1-l) = nkeepfor
            enddo
          else
!
!  Assign number of i/l combinations to preceeding potentials
!
            do l = 1,niltor
              ilnum(nlist4md+1-l) = niltor
            enddo
          endif
        enddo jjloop
      enddo lkloop
    enddo ljloop
  endif
!
!  End of outer loops
!
5 continue
!
!  Free local memory
!
  deallocate(zvec,stat=status)
  if (status/=0) call deallocate_error('setlist4','zvec')
  deallocate(yvec,stat=status)
  if (status/=0) call deallocate_error('setlist4','yvec')
  deallocate(xvec,stat=status)
  if (status/=0) call deallocate_error('setlist4','xvec')
  deallocate(nuniquektyp,stat=status)
  if (status/=0) call deallocate_error('setlist4','nuniquektyp')
  deallocate(nuniquelptr,stat=status)
  if (status/=0) call deallocate_error('setlist4','nuniquelptr')
  deallocate(nuniquekptr,stat=status)
  if (status/=0) call deallocate_error('setlist4','nuniquekptr')
  deallocate(nuniqueiptr,stat=status)
  if (status/=0) call deallocate_error('setlist4','nuniqueiptr')
  deallocate(nuniquelind,stat=status)
  if (status/=0) call deallocate_error('setlist4','nuniquelind')
  deallocate(nuniquekind,stat=status)
  if (status/=0) call deallocate_error('setlist4','nuniquekind')
  deallocate(nuniqueiind,stat=status)
  if (status/=0) call deallocate_error('setlist4','nuniqueiind')
  deallocate(nptrnfornonoop,stat=status)
  if (status/=0) call deallocate_error('setlist4','nptrnfornonoop')
  deallocate(ntypmiddle,stat=status)
  if (status/=0) call deallocate_error('setlist4','ntypmiddle')
  deallocate(natmiddle,stat=status)
  if (status/=0) call deallocate_error('setlist4','natmiddle')
!****************************
!  Out of plane potentials  *
!****************************
  if (noofp.gt.0) call setlistoop
!**********************
!  Improper torsions  *
!**********************
  if (nimproper.gt.0) call setlistimp
!
!  Timing
!
  time2 = g_cpu_time()
  tfour = tfour + time2 - time1
!
  return
  end
