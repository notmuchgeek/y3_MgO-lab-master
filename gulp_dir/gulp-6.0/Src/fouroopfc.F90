  subroutine fouroopfc
!
!  Subroutine for four-body phonons from out of plane potentials.
!  Unphased version.
!
!  Strategy - sift by potential first, then cutoffs
!
!  10/14 Created from fouroopp
!   8/17 Angle/force constant permutations corrected for n4ty 14 & 16
!   8/17 Correction for degtorad factors added for potentials with angles
!   2/18 Trace added
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
!  Copyright Curtin University 2018
!
!  Julian Gale, CIC, Curtin University, Feb 2018
!
  use current
  use derivatives
  use four
  use g_constants,    only : degtorad
  use molecule
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ii
  integer(i4)                                  :: isgn
  integer(i4)                                  :: indm
  integer(i4)                                  :: indmj
  integer(i4)                                  :: indmk
  integer(i4)                                  :: indml
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
  integer(i4)                                  :: jjmax
  integer(i4)                                  :: jloop
  integer(i4)                                  :: jx
  integer(i4)                                  :: jxx
  integer(i4)                                  :: jyy
  integer(i4)                                  :: jzz
  integer(i4)                                  :: k
  integer(i4)                                  :: kb(6,6)
  integer(i4)                                  :: ki
  integer(i4)                                  :: kj
  integer(i4)                                  :: kk
  integer(i4)                                  :: kl
  integer(i4)                                  :: kmax
  integer(i4)                                  :: kx
  integer(i4)                                  :: kxx
  integer(i4)                                  :: kyy
  integer(i4)                                  :: kzz
  integer(i4)                                  :: l
  integer(i4)                                  :: l4
  integer(i4)                                  :: lj
  integer(i4)                                  :: lk
  integer(i4)                                  :: ll
  integer(i4)                                  :: lu
  integer(i4)                                  :: llmax
  integer(i4)                                  :: lmax
  integer(i4)                                  :: lx
  integer(i4)                                  :: n
  integer(i4)                                  :: n11
  integer(i4)                                  :: n21
  integer(i4)                                  :: n31
  integer(i4)                                  :: n41
  integer(i4)                                  :: n22
  integer(i4)                                  :: n32
  integer(i4)                                  :: n42
  integer(i4)                                  :: n1x
  integer(i4)                                  :: n2x
  integer(i4)                                  :: n3x
  integer(i4)                                  :: n4x
  integer(i4)                                  :: n3vec(3,4)
  integer(i4)                                  :: nbtypeij
  integer(i4)                                  :: nbtypeik
  integer(i4)                                  :: nbtypeil
  integer(i4)                                  :: nbtypeij2
  integer(i4)                                  :: nbtypeik2
  integer(i4)                                  :: nbtypeil2
  integer(i4)                                  :: ncind12m
  integer(i4)                                  :: ncind12p
  integer(i4)                                  :: ncind13m
  integer(i4)                                  :: ncind13p
  integer(i4)                                  :: ncind14m
  integer(i4)                                  :: ncind14p
  integer(i4)                                  :: ncind23m
  integer(i4)                                  :: ncind23p
  integer(i4)                                  :: ncind24m
  integer(i4)                                  :: ncind24p
  integer(i4)                                  :: ncind34m
  integer(i4)                                  :: ncind34p
  integer(i4)                                  :: nfortype
  integer(i4)                                  :: ni
  integer(i4)                                  :: nj
  integer(i4)                                  :: nk
  integer(i4)                                  :: nl 
  integer(i4)                                  :: nmi 
  integer(i4)                                  :: nmj 
  integer(i4)                                  :: nmk 
  integer(i4)                                  :: nml 
  integer(i4)                                  :: nt1
  integer(i4)                                  :: nt2
  integer(i4)                                  :: nt3
  integer(i4)                                  :: nt4
  integer(i4)                                  :: ntp2
  integer(i4)                                  :: ntp3
  integer(i4)                                  :: ntp4
  integer(i4)                                  :: ntyp1
  integer(i4)                                  :: ntyp2
  integer(i4)                                  :: ntyp3
  integer(i4)                                  :: ntyp4
  integer(i4)                                  :: ntypi
  integer(i4)                                  :: ntypj
  integer(i4)                                  :: ntypk
  integer(i4)                                  :: ntypl
  integer(i4)                                  :: nunique
  integer(i4), dimension(:), allocatable       :: nuniqueptr
  integer(i4)                                  :: status
  logical                                      :: l2bonds
  logical                                      :: lbonded
  logical                                      :: lbtyp
  logical                                      :: lintra_only
  logical                                      :: linter_only
  logical                                      :: ljbond
  logical                                      :: lmatch
  logical                                      :: lmatch2
  logical                                      :: lmatch3
  logical                                      :: lmatchanyof2
  logical                                      :: lmatchanyof3
  logical                                      :: lmolok
  logical                                      :: lneedmol
  logical                                      :: lsamemol
  logical                                      :: lunique
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
  real(dp)                                     :: rkfor
  real(dp)                                     :: rkfor4
  real(dp)                                     :: rko
  real(dp)                                     :: rn
  real(dp)                                     :: t12
  real(dp)                                     :: t13
  real(dp)                                     :: t14
  real(dp)                                     :: t23
  real(dp)                                     :: t24
  real(dp)                                     :: t34
  real(dp)                                     :: tr1
  real(dp)                                     :: tr2
  real(dp)                                     :: tr3
  real(dp)                                     :: tr1min
  real(dp)                                     :: tr2min
  real(dp)                                     :: tr3min
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
  real(dp)                                     :: x31t
  real(dp)                                     :: y31t
  real(dp)                                     :: z31t
  real(dp)                                     :: x41
  real(dp)                                     :: y41
  real(dp)                                     :: z41
  real(dp)                                     :: x41t
  real(dp)                                     :: y41t
  real(dp)                                     :: z41t
  real(dp)                                     :: x32
  real(dp)                                     :: y32
  real(dp)                                     :: z32
  real(dp)                                     :: x42
  real(dp)                                     :: y42
  real(dp)                                     :: z42
  real(dp)                                     :: x43
  real(dp)                                     :: y43
  real(dp)                                     :: z43
  real(dp)                                     :: xc1
  real(dp)                                     :: yc1
  real(dp)                                     :: zc1
  real(dp)                                     :: xc2t
  real(dp)                                     :: yc2t
  real(dp)                                     :: zc2t
  real(dp)                                     :: xc3t
  real(dp)                                     :: yc3t
  real(dp)                                     :: zc3t
  real(dp)                                     :: xc4t
  real(dp)                                     :: yc4t
  real(dp)                                     :: zc4t
!
  data kb/1,2,3,4,5,6,2,7,8,9,10,11,3,8,12,13,14,15,4,9,13,16,17, &
          18,5,10,14,17,19,20,6,11,15,18,20,21/
  data n3vec/1,2,3,1,4,5,2,4,6,3,5,6/
#ifdef TRACE
  call trace_in('fouroopfc')
#endif
!
!  Allocate local memory
!
  allocate(nuniqueptr(numat),stat=status)
  if (status/=0) call outofmemory('fouroopfc','nuniqueptr')
!*************************
!  Loop over potentials  *
!*************************
  pots: do n = 1,nfor
    nfortype = nforty(n)
    if (.not.loutofplane(n)) cycle pots
    ntp2 = 2
    ntp3 = 3
    ntp4 = 4
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
    tr1min = for1min(n)**2
    tr2min = for2min(n)**2
    tr3min = for3min(n)**2
    lbtyp = (mmfexc(n).eq.1)
    rkfor = fork(n)
    rkfor4 = forpoly(1,n)
    lintra_only = (lfintra(n).and..not.lfinter(n))
    linter_only = (lfinter(n).and..not.lfintra(n))
    lneedmol = (lintra_only.or.linter_only.or.lbtyp)
!********************************
!  Loop over middle site 1 / i  *
!********************************
    ix = -2
    iy = -1
    iz = 0
    liloop: do i = 1,numat
      ni = nat(i)
      ntypi = nftype(i)
      oci = occuf(i)
      ix = ix + 3
      iy = iy + 3
      iz = iz + 3
!
!  Check i is allowed for n
!
      if (.not.lmatch(ni,ntypi,nt1,ntyp1,.true.)) cycle liloop
!
!  Only 3 bonds check
!
      if (lbtyp.and.lonly3oop(n).and.nbonds(i).ne.3) cycle liloop
!
!  Set loop range for other atoms
!
      if (lbtyp) then
        ljbond = .true.
        if (nbonds(i).gt.0) then
          nunique = 1
          nuniqueptr(1) = nbonded(1,i)
          do jloop = 2,nbonds(i)
            lunique = .true.
            do lu = 1,nunique
              if (nbonded(jloop,i).eq.nuniqueptr(lu)) lunique = .false.
            enddo
            if (lunique) then
              nunique = nunique + 1
              nuniqueptr(nunique) = nbonded(jloop,i)
            endif
          enddo
          jloop = nunique
        else
          jloop = 0
        endif
      else
        ljbond = .false.
        jloop = numat
      endif
!
!  Skip if jloop is zero
!
      if (jloop.eq.0) cycle liloop
!
!  i has been accepted
!
      xc1 = xclat(i)
      yc1 = yclat(i)
      zc1 = zclat(i)
!
!  Molecule handling
!
      if (lmol.and.lneedmol) then
        nmi = natmol(i)
        if (ndim.gt.0) then
          indm = nmolind(i)
          call mindtoijk(indm,ixi,iyi,izi)
        endif
      endif
!***********************************
!  Loop over first end site 2 / j  *
!***********************************
      ljloop: do lj = 1,jloop
        if (ljbond) then
          j = nuniqueptr(lj)
        else
          j = lj
        endif
        nj = nat(j)
        ntypj = nftype(j)
!
!  Check j is allowed for n
!
        lmatch3 = lmatchanyof3(nj,ntypj,ntp2,nt2,ntyp2,tr1,tr1min,ntp3,nt3,ntyp3,tr2,tr2min,ntp4,nt4,ntyp4,tr3,tr3min)
        if (.not.lmatch3) cycle ljloop
!
!  Set properties for atom j
!
        ocj = occuf(j)
        jx = 3*(j - 1) + 1
!
        if (lmol.and.lneedmol) then
!
!  Molecule handling
!
          nmj = natmol(j)
          if (ndim.gt.0) then
            indmj = nmolind(j)
            call mindtoijk(indmj,ixj,iyj,izj)
            ixj = ixj - ixi
            iyj = iyj - iyi
            izj = izj - izi
          endif
          lmolok = (nmi.eq.nmj.and.nmi.gt.0)
        else
          lmolok = .false.
        endif
!
!  Check for intra and but not in same molecule
!
        if (lintra_only.and..not.lmolok) cycle ljloop
        if (lbtyp.and..not.lmolok) cycle ljloop
!
        xc2t = xclat(j)
        yc2t = yclat(j)
        zc2t = zclat(j)
        x21t = xc2t - xc1
        y21t = yc2t - yc1
        z21t = zc2t - zc1
!
!  Check r21 is OK
!  Loop over cell vectors
!
        iiloop: do ii = 1,iimax
          r212 = (xvec1cell(ii)+x21t)**2 + (yvec1cell(ii)+y21t)**2 + (zvec1cell(ii)+z21t)**2
          if (r212.lt.1d-10) cycle iiloop
!
!  Molecule checking
!
          lbonded = .false.
          if (lmolok) then
            if (ndim.eq.0) then
              if (linter_only) cycle iiloop
              if (lbtyp) then
                call bonded(lbonded,l2bonds,nbtypeij,nbtypeij2,i,j,0_i4,0_i4,0_i4)
                if (.not.lbonded) cycle iiloop
              endif
            else
              call lintoijk(ixx,iyy,izz,ii,imaxl,jmaxl,kmaxl)
              if (lbtyp) then
                call bonded(lbonded,l2bonds,nbtypeij,nbtypeij2,i,j,ixx,iyy,izz)
                if (.not.lbonded) cycle iiloop
                lsamemol = (lbonded.or.l2bonds)
              else
                lsamemol = .false.
              endif
              if (.not.lsamemol) then
                call samemol(lsamemol,nmi,ixx,iyy,izz,ixj,iyj,izj)
              endif
              if (lintra_only.and..not.lsamemol) cycle iiloop
              if (linter_only.and.lsamemol) cycle iiloop
            endif
          endif
!
!  Distance checking
!
          if ((r212.gt.tr1.or.r212.lt.tr1min).and.(.not.lbtyp.or..not.lbonded)) cycle iiloop
          r21 = sqrt(r212)
          x21 = x21t + xvec1cell(ii)
          y21 = y21t + yvec1cell(ii)
          z21 = z21t + zvec1cell(ii)
!
!  Find cell index for 1-2
!
          if (abs(ivec1cell(1,ii)).gt.nd2cell(1).or. &
              abs(ivec1cell(2,ii)).gt.nd2cell(2).or. &
              abs(ivec1cell(3,ii)).gt.nd2cell(3)) then
!
!  Find cell index - if outside user range then assign to central cell
!
            ncind12m = nd2central
            ncind12p = nd2central
          else
!
!  Compute index
!
            ncind12p = nd2cellptr(nd2cell(1)+1+ivec1cell(1,ii), &
                                  nd2cell(2)+1+ivec1cell(2,ii), &
                                  nd2cell(3)+1+ivec1cell(3,ii))
            ncind12m = nd2cellptr(nd2cell(1)+1-ivec1cell(1,ii), &
                                  nd2cell(2)+1-ivec1cell(2,ii), &
                                  nd2cell(3)+1-ivec1cell(3,ii))
          endif
!
          if (ndim.eq.0) then
            kmax = lj - 1
          else
            kmax = lj
          endif
!************************************
!  Loop over second end site 3 / k  *
!************************************
          lkloop: do lk = 1,kmax
            if (ljbond) then
              k = nuniqueptr(lk)
            else
              k = lk
            endif
            nk = nat(k)
            ntypk = nftype(k)
!
!  Check k is allowed for n
!
            lmatch2 = lmatchanyof2(nk,ntypk,ntp3,nt3,ntyp3,tr2,tr2min,ntp4,nt4,ntyp4,tr3,tr3min)
            if (.not.lmatch2) cycle lkloop
!
!  Set properties of atom k
!
            ock = occuf(k)
            kx = 3*(k - 1) + 1
!
            if (lmol.and.lneedmol) then
!
!  Molecule handling
!
              nmk = natmol(k)
              if (ndim.gt.0) then
                indmk = nmolind(k)
                call mindtoijk(indmk,ixk,iyk,izk)
                ixk = ixk - ixi
                iyk = iyk - iyi
                izk = izk - izi
              endif
              lmolok = (nmi.eq.nmk.and.nmi.gt.0)
            else
              lmolok = .false.
            endif
!
!  Check for intra and but not in same molecule
!
            if (lintra_only.and..not.lmolok) cycle lkloop
            if (lbtyp.and..not.lmolok) cycle lkloop
            xc3t = xclat(k)
            yc3t = yclat(k)
            zc3t = zclat(k)
            x31t = xc3t - xc1
            y31t = yc3t - yc1
            z31t = zc3t - zc1
!
            if (j.eq.k) then
              jjmax = ii - 1
            else
              jjmax = iimax
            endif
!
!  Check r31 is OK
!  Loop over cell vectors
!
            jjloop: do jj = 1,jjmax
              r312 = (xvec1cell(jj)+x31t)**2 + (yvec1cell(jj)+y31t)**2 + (zvec1cell(jj)+z31t)**2
              if (r312.lt.1d-10) cycle jjloop
!  
!  Prevent atoms i and k being the same atom
!             
              if (k.eq.i.and.(jj.eq.iimid)) cycle jjloop
!  
!  Prevent atoms j and k being the same atom
!             
              if (k.eq.j.and.jj.eq.ii) cycle jjloop
!
!  Molecule checking
!
              lbonded = .false.
              if (lmolok) then
                if (ndim.eq.0) then
                  if (linter_only) cycle jjloop
                  if (lbtyp) then
                    call bonded(lbonded,l2bonds,nbtypeik,nbtypeik2,i,k,0_i4,0_i4,0_i4)
                    if (.not.lbonded) cycle jjloop
                  endif
                else
                  call lintoijk(jxx,jyy,jzz,jj,imaxl,jmaxl,kmaxl)
                  if (lbtyp) then
                    call bonded(lbonded,l2bonds,nbtypeik,nbtypeik2,i,k,jxx,jyy,jzz)
                    if (.not.lbonded) cycle jjloop
                    lsamemol = (lbonded.or.l2bonds)
                  else
                    lsamemol = .false.
                  endif
                  if (.not.lsamemol) then
                    call samemol(lsamemol,nmi,jxx,jyy,jzz,ixk,iyk,izk)
                  endif
                  if (lintra_only.and..not.lsamemol) cycle jjloop
                  if (linter_only.and.lsamemol) cycle jjloop
                endif
              endif
!
!  Distance checking
!
              if ((r312.gt.tr2.or.r312.lt.tr2min).and.(.not.lbtyp.or..not.lbonded)) cycle jjloop
              r31 = sqrt(r312)
              x31 = x31t + xvec1cell(jj)
              y31 = y31t + yvec1cell(jj)
              z31 = z31t + zvec1cell(jj)
!
!  Find cell index for 1-3
!
              if (abs(ivec1cell(1,jj)).gt.nd2cell(1).or. &
                  abs(ivec1cell(2,jj)).gt.nd2cell(2).or. &
                  abs(ivec1cell(3,jj)).gt.nd2cell(3)) then
!
!  Find cell index - if outside user range then assign to central cell
!
                ncind13m = nd2central
                ncind13p = nd2central
              else
!
!  Compute index
!
                ncind13p = nd2cellptr(nd2cell(1)+1+ivec1cell(1,jj), &
                                      nd2cell(2)+1+ivec1cell(2,jj), &
                                      nd2cell(3)+1+ivec1cell(3,jj))
                ncind13m = nd2cellptr(nd2cell(1)+1-ivec1cell(1,jj), &
                                      nd2cell(2)+1-ivec1cell(2,jj), &
                                      nd2cell(3)+1-ivec1cell(3,jj))
              endif
!
!  Find cell index for 2-3
!             
              if (abs(ivec1cell(1,jj)-ivec1cell(1,ii)).gt.nd2cell(1).or. &
                  abs(ivec1cell(2,jj)-ivec1cell(2,ii)).gt.nd2cell(2).or. &
                  abs(ivec1cell(3,jj)-ivec1cell(3,ii)).gt.nd2cell(3)) then
!             
!  Find cell index - if outside user range then assign to central cell
!  
                ncind23m = nd2central
                ncind23p = nd2central
              else
!                 
!  Compute index
!                 
                ncind23p = nd2cellptr(nd2cell(1)+1+ivec1cell(1,jj)-ivec1cell(1,ii), &
                                      nd2cell(2)+1+ivec1cell(2,jj)-ivec1cell(2,ii), &
                                      nd2cell(3)+1+ivec1cell(3,jj)-ivec1cell(3,ii))
                ncind23m = nd2cellptr(nd2cell(1)+1-ivec1cell(1,jj)+ivec1cell(1,ii), &
                                      nd2cell(2)+1-ivec1cell(2,jj)+ivec1cell(2,ii), &
                                      nd2cell(3)+1-ivec1cell(3,jj)+ivec1cell(3,ii))
              endif
!
              if (ndim.eq.0) then
                lmax = lk - 1
              else
                lmax = lk
              endif
!**********************************
!  Loop over last end site 4 / l  *
!**********************************
              l4loop: do l4 = 1,lmax
                if (ljbond) then
                  l = nuniqueptr(l4)
                else
                  l = l4
                endif
                nl = nat(l)
                ntypl = nftype(l)
!
!  Check l is allowed for n
!
                if (.not.lmatch(nl,ntypl,nt4,ntyp4,.true.)) cycle l4loop
!
!  Set properties of atom l
!
                ocl = occuf(l)
                lx = 3*(l - 1) + 1
!
                if (lmol.and.lneedmol) then
!
!  Molecule handling
!
                  nml = natmol(l)
                  if (ndim.gt.0) then
                    indml = nmolind(l)
                    call mindtoijk(indml,ixl,iyl,izl)
                    ixl = ixl - ixi
                    iyl = iyl - iyi
                    izl = izl - izi
                  endif
                  lmolok = (nmi.eq.nml.and.nmi.gt.0)
                else
                  lmolok = .false.
                endif
!
!  Check for intra and but not in same molecule
!
                if (lintra_only.and..not.lmolok) cycle l4loop
                if (lbtyp.and..not.lmolok) cycle l4loop
                xc4t = xclat(l)
                yc4t = yclat(l)
                zc4t = zclat(l)
                x41t = xc4t - xc1
                y41t = yc4t - yc1
                z41t = zc4t - zc1
!
                if (k.eq.l) then
                  llmax = jj - 1
                else
                  llmax = iimax
                endif
!
!  Check r41 is OK
!  Loop over cell vectors
!
                llloop: do ll = 1,llmax
                  r412 = (xvec1cell(ll)+x41t)**2 + (yvec1cell(ll)+y41t)**2 + (zvec1cell(ll)+z41t)**2
                  if (r412.lt.1d-10) cycle llloop
!
!  Prevent atoms i and l being the same atom
!
                  if (l.eq.i.and.ll.eq.iimid) cycle llloop
!
!  Prevent atoms j and l being the same atom
!
                  if (l.eq.j.and.ll.eq.ii) cycle llloop
!
!  Prevent atoms k and l being the same atom
!
                  if (l.eq.k.and.ll.eq.jj) cycle llloop
!
!  Molecule checking
!
                  lbonded = .false.
                  if (lmolok) then
                    if (ndim.eq.0) then
                      if (linter_only) cycle llloop
                      if (lbtyp) then
                        call bonded(lbonded,l2bonds,nbtypeil,nbtypeil2,i,l,0_i4,0_i4,0_i4)
                        if (.not.lbonded) cycle llloop
                      endif
                    else
                      call lintoijk(kxx,kyy,kzz,ll,imaxl,jmaxl,kmaxl)
                      if (lbtyp) then
                        call bonded(lbonded,l2bonds,nbtypeil,nbtypeil2,i,l,kxx,kyy,kzz)
                        if (.not.lbonded) cycle llloop
                        lsamemol = (lbonded.or.l2bonds)
                      else
                        lsamemol = .false.
                      endif
                      if (.not.lsamemol) then
                        call samemol(lsamemol,nmi,kxx,kyy,kzz,ixl,iyl,izl)
                      endif
                      if (lintra_only.and..not.lsamemol) cycle llloop
                      if (linter_only.and.lsamemol) cycle llloop
                    endif
                  endif
!
!  Distance checking
!
                  if ((r412.gt.tr3.or.r412.lt.tr3min).and.(.not.lbtyp.or..not.lbonded)) cycle llloop
!*********************************
!  Valid four-body term located  *
!*********************************
!
!  Define remaining vectors between atoms
!
                  x41 = x41t + xvec1cell(ll)
                  y41 = y41t + yvec1cell(ll)
                  z41 = z41t + zvec1cell(ll)
                  r412 = x41*x41 + y41*y41 + z41*z41
                  r41 = sqrt(r412)
                  x43 = x41 - x31
                  y43 = y41 - y31
                  z43 = z41 - z31
                  r432 = x43*x43 + y43*y43 + z43*z43
                  r43 = sqrt(r432)
                  x32 = x31 - x21
                  y32 = y31 - y21
                  z32 = z31 - z21
                  r322 = x32*x32 + y32*y32 + z32*z32
                  r32 = sqrt(r322)
                  x42 = x43 + x32
                  y42 = y43 + y32
                  z42 = z43 + z32
                  r422 = x42*x42 + y42*y42 + z42*z42
                  r42 = sqrt(r422)
!
!  Find cell index for 1-4
!
                  if (abs(ivec1cell(1,ll)).gt.nd2cell(1).or. &
                      abs(ivec1cell(2,ll)).gt.nd2cell(2).or. &
                      abs(ivec1cell(3,ll)).gt.nd2cell(3)) then
!
!  Find cell index - if outside user range then assign to central cell
!
                    ncind14m = nd2central
                    ncind14p = nd2central
                  else
!
!  Compute index
!
                    ncind14p = nd2cellptr(nd2cell(1)+1+ivec1cell(1,ll), &
                                          nd2cell(2)+1+ivec1cell(2,ll), &
                                          nd2cell(3)+1+ivec1cell(3,ll))
                    ncind14m = nd2cellptr(nd2cell(1)+1-ivec1cell(1,ll), &
                                          nd2cell(2)+1-ivec1cell(2,ll), &
                                          nd2cell(3)+1-ivec1cell(3,ll))
                  endif
!
!  Find cell index for 2-4
!
                  if (abs(ivec1cell(1,ll)-ivec1cell(1,ii)).gt.nd2cell(1).or. &
                      abs(ivec1cell(2,ll)-ivec1cell(2,ii)).gt.nd2cell(2).or. &
                      abs(ivec1cell(3,ll)-ivec1cell(3,ii)).gt.nd2cell(3)) then
!
!  Find cell index - if outside user range then assign to central cell
!
                    ncind24m = nd2central
                    ncind24p = nd2central
                  else
!
!  Compute index
!
                    ncind24p = nd2cellptr(nd2cell(1)+1+ivec1cell(1,ll)-ivec1cell(1,ii), &
                                          nd2cell(2)+1+ivec1cell(2,ll)-ivec1cell(2,ii), &
                                          nd2cell(3)+1+ivec1cell(3,ll)-ivec1cell(3,ii))
                    ncind24m = nd2cellptr(nd2cell(1)+1-ivec1cell(1,ll)+ivec1cell(1,ii), &
                                          nd2cell(2)+1-ivec1cell(2,ll)+ivec1cell(2,ii), &
                                          nd2cell(3)+1-ivec1cell(3,ll)+ivec1cell(3,ii))
                  endif
!
!  Find cell index for 3-4
!
                  if (abs(ivec1cell(1,ll)-ivec1cell(1,jj)).gt.nd2cell(1).or. &
                      abs(ivec1cell(2,ll)-ivec1cell(2,jj)).gt.nd2cell(2).or. &
                      abs(ivec1cell(3,ll)-ivec1cell(3,jj)).gt.nd2cell(3)) then
!
!  Find cell index - if outside user range then assign to central cell
!
                    ncind34m = nd2central
                    ncind34p = nd2central
                  else
!
!  Compute index
!
                    ncind34p = nd2cellptr(nd2cell(1)+1+ivec1cell(1,ll)-ivec1cell(1,jj), &
                                          nd2cell(2)+1+ivec1cell(2,ll)-ivec1cell(2,jj), &
                                          nd2cell(3)+1+ivec1cell(3,ll)-ivec1cell(3,jj))
                    ncind34m = nd2cellptr(nd2cell(1)+1-ivec1cell(1,ll)+ivec1cell(1,jj), &
                                          nd2cell(2)+1-ivec1cell(2,ll)+ivec1cell(2,jj), &
                                          nd2cell(3)+1-ivec1cell(3,ll)+ivec1cell(3,jj))
                  endif
!*****************************************************
!  Calculate derivative terms for first derivatives  *
!*****************************************************
                  ofct = oci*ocj*ock*ocl
                  rko = rkfor*ofct
                  if (nfortype.eq.14.or.nfortype.eq.16) then
                    if (ntp2.eq.2.and.ntp3.eq.3) then
                      fpoly(1) = forpoly(1,n)*ofct
                      fpoly(2) = forpoly(2,n)*ofct
                      fpoly(3) = forpoly(3,n)*degtorad
                      fpoly(4) = forpoly(4,n)*degtorad
                      fpoly(5) = forpoly(5,n)*degtorad
                    elseif (ntp2.eq.2.and.ntp3.eq.4) then
                      fpoly(1) = forpoly(2,n)*ofct
                      fpoly(2) = forpoly(1,n)*ofct
                      fpoly(3) = forpoly(4,n)*degtorad
                      fpoly(4) = forpoly(3,n)*degtorad
                      fpoly(5) = forpoly(5,n)*degtorad
                    elseif (ntp2.eq.3.and.ntp3.eq.2) then
                      fpoly(1) = forpoly(2,n)*ofct
                      fpoly(2) = forpoly(1,n)*ofct
                      fpoly(3) = forpoly(3,n)*degtorad
                      fpoly(4) = forpoly(5,n)*degtorad
                      fpoly(5) = forpoly(4,n)*degtorad
                    elseif (ntp2.eq.3.and.ntp3.eq.4) then
                      fpoly(1) = forpoly(2,n)*ofct
                      fpoly(2) = rko
                      rko = forpoly(1,n)*ofct
                      fpoly(3) = forpoly(5,n)*degtorad
                      fpoly(4) = forpoly(3,n)*degtorad
                      fpoly(5) = forpoly(4,n)*degtorad
                    elseif (ntp2.eq.4.and.ntp3.eq.2) then
                      fpoly(1) = rko
                      rko = forpoly(2,n)*ofct
                      fpoly(2) = forpoly(1,n)*ofct
                      fpoly(3) = forpoly(4,n)*degtorad
                      fpoly(4) = forpoly(5,n)*degtorad
                      fpoly(5) = forpoly(3,n)*degtorad
                    elseif (ntp2.eq.4.and.ntp3.eq.3) then
                      fpoly(2) = rko
                      rko = forpoly(2,n)*ofct
                      fpoly(1) = forpoly(1,n)*ofct
                      fpoly(3) = forpoly(5,n)*degtorad
                      fpoly(4) = forpoly(4,n)*degtorad
                      fpoly(5) = forpoly(3,n)*degtorad
                    endif
                  elseif (nfortype.eq.12) then
                    fpoly(1) = forpoly(1,n)*degtorad
                  elseif (nfortype.eq.15) then
                    fpoly(1) = forpoly(1,n)
                    fpoly(2) = forpoly(2,n)
                    fpoly(3) = forpoly(3,n)
                  else
                    fpoly(1) = rkfor4*ofct
                    fpoly(2) = forpoly(2,n)
                  endif
                  call fourbody(n,nfortype,r21,r31,r41,r32,r42,r43,eterm,e1d,e2d,e3d,rko,rn,phi0, &
                                isgn,fpoly,.true.,.true.,.false.)
!
!  Define pointers to elements of second derivative matrices
!
                  n1x = ix
                  n2x = jx
                  n3x = kx
                  n4x = lx
!********************************
!  Generate second derivatives  *
!********************************
!
!  New vector array between atoms to handle sign
!
!  Atom 1
!
                  vec(1,1,1) = - x21
                  vec(2,1,1) = - y21
                  vec(3,1,1) = - z21
                  vec(1,2,1) = - x31
                  vec(2,2,1) = - y31
                  vec(3,2,1) = - z31
                  vec(1,3,1) = - x41
                  vec(2,3,1) = - y41
                  vec(3,3,1) = - z41
!
!  Atom 2
!
                  vec(1,1,2) = x21
                  vec(2,1,2) = y21
                  vec(3,1,2) = z21
                  vec(1,2,2) = - x32
                  vec(2,2,2) = - y32
                  vec(3,2,2) = - z32
                  vec(1,3,2) = - x42
                  vec(2,3,2) = - y42
                  vec(3,3,2) = - z42
!
!  Atom 3
!
                  vec(1,1,3) = x31
                  vec(2,1,3) = y31
                  vec(3,1,3) = z31
                  vec(1,2,3) = x32
                  vec(2,2,3) = y32
                  vec(3,2,3) = z32
                  vec(1,3,3) = - x43
                  vec(2,3,3) = - y43
                  vec(3,3,3) = - z43
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
!
!  First term
!
                    d2cell(n21,n11,ncind12p) = d2cell(n21,n11,ncind12p) - e1d(1)
                    d2cell(n31,n11,ncind13p) = d2cell(n31,n11,ncind13p) - e1d(2)
                    d2cell(n41,n11,ncind14p) = d2cell(n41,n11,ncind14p) - e1d(3)
                    d2cell(n31,n21,ncind23p) = d2cell(n31,n21,ncind23p) - e1d(4)
                    d2cell(n41,n21,ncind24p) = d2cell(n41,n21,ncind24p) - e1d(5)
                    d2cell(n41,n31,ncind34p) = d2cell(n41,n31,ncind34p) - e1d(6)
!
                    d2cell(n11,n21,ncind12m) = d2cell(n11,n21,ncind12m) - e1d(1)
                    d2cell(n11,n31,ncind13m) = d2cell(n11,n31,ncind13m) - e1d(2)
                    d2cell(n11,n41,ncind14m) = d2cell(n11,n41,ncind14m) - e1d(3)
                    d2cell(n21,n31,ncind23m) = d2cell(n21,n31,ncind23m) - e1d(4)
                    d2cell(n21,n41,ncind24m) = d2cell(n21,n41,ncind24m) - e1d(5)
                    d2cell(n31,n41,ncind34m) = d2cell(n31,n41,ncind34m) - e1d(6)
!
!  Loop over second coordinate
!
                    do kl = 1,3
                      n22 = n2x - 1 + kl
                      n32 = n3x - 1 + kl
                      n42 = n4x - 1 + kl
!
!  Sum over vectors atom - atom second derivatives
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
!  Add phased terms to second derivative matrix
!
                      d2cell(n22,n11,ncind12p) = d2cell(n22,n11,ncind12p) + t12
                      d2cell(n32,n11,ncind13p) = d2cell(n32,n11,ncind13p) + t13
                      d2cell(n42,n11,ncind14p) = d2cell(n42,n11,ncind14p) + t14
                      d2cell(n32,n21,ncind23p) = d2cell(n32,n21,ncind23p) + t23
                      d2cell(n42,n21,ncind24p) = d2cell(n42,n21,ncind24p) + t24
                      d2cell(n42,n31,ncind34p) = d2cell(n42,n31,ncind34p) + t34
!
                      d2cell(n11,n22,ncind12m) = d2cell(n11,n22,ncind12m) + t12
                      d2cell(n11,n32,ncind13m) = d2cell(n11,n32,ncind13m) + t13
                      d2cell(n11,n42,ncind14m) = d2cell(n11,n42,ncind14m) + t14
                      d2cell(n21,n32,ncind23m) = d2cell(n21,n32,ncind23m) + t23
                      d2cell(n21,n42,ncind24m) = d2cell(n21,n42,ncind24m) + t24
                      d2cell(n31,n42,ncind34m) = d2cell(n31,n42,ncind34m) + t34
                    enddo
                  enddo
!
!  End of inner loops over atoms and cell vectors
!
                enddo llloop
              enddo l4loop
            enddo jjloop
          enddo lkloop
        enddo iiloop
      enddo ljloop
    enddo liloop
!
!  End of outer loops
!
  enddo pots
!
!  Free local memory
!
  deallocate(nuniqueptr,stat=status)
  if (status/=0) call deallocate_error('fouroopfc','nuniqueptr')
!
!  All tidying up of derivatives is handled by four so we can just return here
!
#ifdef TRACE
  call trace_out('fouroopfc')
#endif
  return
  end
