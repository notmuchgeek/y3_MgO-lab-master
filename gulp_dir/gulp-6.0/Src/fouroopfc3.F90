  subroutine fouroopfc3(maxlim,maxlimcg,fc3,ind2ptr,ind3ptr,delta)
!
!  Subroutine for third derivatives of four-body energy from out of plane potentials.
!  Called from fourfc3/thirdorderfc3.
!
!   8/15 Created from fouroopp and fourfc3
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
  use configurations, only : nsuperghost
  use current
  use derivatives
  use element,        only : maxele
  use four
  use g_constants,    only : degtorad
  use molecule
  use optimisation
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4),                   intent(in)    :: maxlim
  integer(i4),                   intent(in)    :: maxlimcg
  integer(i4),                   intent(in)    :: ind2ptr(6,6)
  integer(i4),                   intent(in)    :: ind3ptr(6,6,6)
  real(dp),                      intent(in)    :: delta(4,4,6)
  real(dp),                      intent(inout) :: fc3(maxlim,maxlim,maxlimcg)
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ia
  integer(i4)                                  :: ic
  integer(i4)                                  :: ii
  integer(i4)                                  :: ii1
  integer(i4)                                  :: ii2
  integer(i4)                                  :: ii3
  integer(i4)                                  :: isgn
  integer(i4)                                  :: indm
  integer(i4)                                  :: indmj
  integer(i4)                                  :: indmk
  integer(i4)                                  :: indml
  integer(i4)                                  :: ig
  integer(i4)                                  :: igx
  integer(i4)                                  :: igxyz
  integer(i4)                                  :: ix
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
  integer(i4)                                  :: ja
  integer(i4)                                  :: jc
  integer(i4)                                  :: jg
  integer(i4)                                  :: jgx
  integer(i4)                                  :: jj
  integer(i4)                                  :: jjmax
  integer(i4)                                  :: jloop
  integer(i4)                                  :: jx
  integer(i4)                                  :: jxx
  integer(i4)                                  :: jyy
  integer(i4)                                  :: jzz
  integer(i4)                                  :: jxyz
  integer(i4)                                  :: k
  integer(i4)                                  :: ka
  integer(i4)                                  :: kc
  integer(i4)                                  :: kg
  integer(i4)                                  :: kgx
  integer(i4)                                  :: kmax
  integer(i4)                                  :: kx
  integer(i4)                                  :: kxx
  integer(i4)                                  :: kyy
  integer(i4)                                  :: kzz
  integer(i4)                                  :: kxyz
  integer(i4)                                  :: l
  integer(i4)                                  :: l4
  integer(i4)                                  :: lg
  integer(i4)                                  :: lgx
  integer(i4)                                  :: lj
  integer(i4)                                  :: lk
  integer(i4)                                  :: ll
  integer(i4)                                  :: lu
  integer(i4)                                  :: llmax
  integer(i4)                                  :: lmax
  integer(i4)                                  :: lx
  integer(i4)                                  :: n
  integer(i4)                                  :: nbtypeij
  integer(i4)                                  :: nbtypeik
  integer(i4)                                  :: nbtypeil
  integer(i4)                                  :: nbtypeij2
  integer(i4)                                  :: nbtypeik2
  integer(i4)                                  :: nbtypeil2
  integer(i4)                                  :: nfortype
  integer(i4)                                  :: nghostcell    ! Number of ghost cells in the full cell
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
  logical                                      :: lghosti
  logical                                      :: lghostj
  logical                                      :: lghostk
  logical                                      :: lghostl
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
  real(dp)                                     :: rkfor
  real(dp)                                     :: rkfor4
  real(dp)                                     :: rko
  real(dp)                                     :: rn
  real(dp)                                     :: tr1
  real(dp)                                     :: tr2
  real(dp)                                     :: tr3
  real(dp)                                     :: tr1min
  real(dp)                                     :: tr2min
  real(dp)                                     :: tr3min
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
#ifdef TRACE
  call trace_in('fouroopfc3')
#endif
!
!  Find number of ghost cells
!
  nghostcell = nsuperghost(1,ncf)*nsuperghost(2,ncf)*nsuperghost(3,ncf)
!
!  Allocate local memory
!
  allocate(nuniqueptr(numat),stat=status)
  if (status/=0) call outofmemory('fouroopfc3','nuniqueptr')
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
    liloop: do i = 1,numat
      ni = nat(i)
      ntypi = nftype(i)
      oci = occuf(i)
!
!  Check i is allowed for n
!
      if (.not.lmatch(ni,ntypi,nt1,ntyp1,.true.)) cycle liloop
!
!  Only 3 bonds check
!
      if (lbtyp.and.lonly3oop(n).and.nbonds(i).ne.3) cycle liloop
!
!  Is i a ghost atom?
!
      lghosti = (mod(i-1,nghostcell).eq.0.and.ni.le.maxele)
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
      ix = 3*(i-1) + 1
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
!  Is j a ghost atom?
!
        lghostj = (mod(j-1,nghostcell).eq.0.and.nj.le.maxele)
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
!  Is k a ghost atom?
!
            lghostk = (mod(k-1,nghostcell).eq.0.and.nk.le.maxele)
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
!  Is l a ghost atom?
!
                lghostl = (mod(l-1,nghostcell).eq.0.and.nl.le.maxele)
!
!  If none of the atoms are ghosts then cycle
!
                if (.not.lghosti.and..not.lghostj.and..not.lghostk.and..not.lghostl) cycle l4loop
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
                                isgn,fpoly,.true.,.true.,.true.)
!
!  Define pointers to elements of second derivative matrices
!
                  if (lghosti) then
                    ig = (i - 1)/nghostcell + 1
                    igx = 3*(ig-1) + 1
                  endif
!
                  if (lghostj) then
                    jg = (j - 1)/nghostcell + 1
                    jgx = 3*(jg-1) + 1
                  endif
!
                  if (lghostk) then
                    kg = (k - 1)/nghostcell + 1
                    kgx = 3*(kg-1) + 1
                  endif
!
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
  if (status/=0) call deallocate_error('fouroopfc3','nuniqueptr')
#ifdef TRACE
  call trace_out('fouroopfc3')
#endif
!
  return
  end
