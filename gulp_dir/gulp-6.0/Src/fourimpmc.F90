  subroutine fourimpmc(eimp,ltrialatom)
!
!  Subroutine for four-body improper torsion energy - subset version 
!
!   7/13 Created from fourmc
!  12/14 Corrected so that ll vector is added to 42 and not 43
!   2/18 Trace added
!   4/19 Call to lmatchpair modified to include wildcard argument
!   8/19 ntrialatom removed from arguments
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
!  Julian Gale, CIC, Curtin University, August 2019
!
  use g_constants
  use current
  use four
  use iochannels,     only : ioout
  use molecule
  use parallel
  use symmetry
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  logical,     intent(in)                      :: ltrialatom(numat)
  real(dp), intent(inout)                      :: eimp
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
  integer(i4)                                  :: jxx
  integer(i4)                                  :: jyy
  integer(i4)                                  :: jzz
  integer(i4)                                  :: k
  integer(i4)                                  :: kk
  integer(i4)                                  :: kloop
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
  integer(i4)                                  :: n
  integer(i4), dimension(:), allocatable       :: natmiddle
  integer(i4), dimension(:), allocatable       :: ntypmiddle
  integer(i4)                                  :: nbtypeji
  integer(i4)                                  :: nbtypejk
  integer(i4)                                  :: nbtypekl
  integer(i4)                                  :: nbtypeji2
  integer(i4)                                  :: nbtypejk2
  integer(i4)                                  :: nbtypekl2
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
  logical                                      :: lsamemolij
  logical                                      :: lsamemoljk
  logical                                      :: lsamemolkl
  logical                                      :: lunique
  real(dp)                                     :: cut
  real(dp)                                     :: cutmax
  real(dp)                                     :: e1d(6)
  real(dp)                                     :: e2d(21)
  real(dp)                                     :: e3d(56)
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
#ifdef TRACE
  call trace_in('fourimpmc')
#endif
!
  lmolloc = (nmol.gt.0)
!
!  Allocate local memory
!
  allocate(natmiddle(2*nfor),stat=status)
  if (status/=0) call outofmemory('fourimpmc','natmiddle')
  allocate(ntypmiddle(2*nfor),stat=status)
  if (status/=0) call outofmemory('fourimpmc','ntypmiddle')
  allocate(nptrnforimp(nfor),stat=status)
  if (status/=0) call outofmemory('fourimpmc','nptrnforimp')
  allocate(nuniqueiptr(numat),stat=status)
  if (status/=0) call outofmemory('fourimpmc','nuniqueiptr')
  allocate(nuniquekptr(numat),stat=status)
  if (status/=0) call outofmemory('fourimpmc','nuniquekptr')
  allocate(nuniquelptr(numat),stat=status)
  if (status/=0) call outofmemory('fourimpmc','nuniquelptr')
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
  lallbtyp = .true.
  lanybtyp = .false.
  lanyneedmol = .false.
  do n = 1,nfor
    if (mmfexc(n).eq.2) then
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
  ljloop: do j = procid+1,numat,nprocs
    nj = nat(j)
    ntypj = nftype(j)
!
!  Check whether species may be valid
!
    if (.not.lmatchany(nj,ntypj,nmiddle,natmiddle,ntypmiddle,.true.)) cycle ljloop
    ocj = occuf(j)
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
!      if (ndim.gt.0) then
        indmj = nmolind(j)
        call mindtoijk(indmj,ixj,iyj,izj)
!      endif
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
      jjloop: do jj = 1,iimax
        r322 = (xvec1cell(jj)+x32t)**2 + (yvec1cell(jj)+y32t)**2 + (zvec1cell(jj)+z32t)**2
        if (r322.lt.1d-12) cycle jjloop
        if (k.eq.j.and.jj.eq.iimid) cycle jjloop
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
        if (r322.gt.tr2max.and.(.not.lanybtyp.or..not.lbondedjk)) cycle jjloop
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
            oci = occuf(i)
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
                nl = nat(l)
                ntypl = nftype(l)
                ocl = occuf(l)
!
!  Check l is allowed for n
!             
                if (.not.lmatch(nl,ntypl,nt4,ntyp4,.true.)) cycle l4loop
!
!  If no atoms are in trial set then skip
!
                if (.not.ltrialatom(i).and..not.ltrialatom(j).and..not.ltrialatom(k).and..not.ltrialatom(l)) cycle l4loop
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
                x42t = xc4t - xc2
                y42t = yc4t - yc2
                z42t = zc4t - zc2
!
!  Check r43 is OK
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
                        rko,rn,phi0o,isgn,fpoly,.false.,.false.,.false.)
          eimp = eimp + eterm
!
!  Output energy contribution
!
          if (lPrintFour) then
            write(ioout,'(4x,4i12,1x,f22.10)') i,j,k,l,eterm
          endif
        enddo
      enddo jjloop
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
  deallocate(nuniquelptr,stat=status)
  if (status/=0) call deallocate_error('fourimpmc','nuniquelptr')
  deallocate(nuniquekptr,stat=status)
  if (status/=0) call deallocate_error('fourimpmc','nuniquekptr')
  deallocate(nuniqueiptr,stat=status)
  if (status/=0) call deallocate_error('fourimpmc','nuniqueiptr')
  deallocate(nptrnforimp,stat=status)
  if (status/=0) call deallocate_error('fourimpmc','nptrnforimp')
  deallocate(ntypmiddle,stat=status)
  if (status/=0) call deallocate_error('fourimpmc','ntypmiddle')
  deallocate(natmiddle,stat=status)
  if (status/=0) call deallocate_error('fourimpmc','natmiddle')
#ifdef TRACE
  call trace_out('fourimpmc')
#endif
!
  return
  end
