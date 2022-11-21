  subroutine sixpd(xkv,ykv,zkv)
!
!  Subroutine to calculate six-body potential contribution to phonons
!  Distributed memory parallel version.
!
!   1/17 Created from sixp.f
!   1/17 Parallelisation added
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
  use control,       only : lgroupvelocity
  use current
  use derivatives
  use molecule
  use parallel
  use six
  use times,         only : tsix
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  real(dp), intent(in)                         :: xkv
  real(dp), intent(in)                         :: ykv
  real(dp), intent(in)                         :: zkv
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ia
  integer(i4)                                  :: ialp
  integer(i4)                                  :: ib
  integer(i4)                                  :: ibet
  integer(i4)                                  :: ii
  integer(i4)                                  :: iloc
  integer(i4)                                  :: indm
  integer(i4)                                  :: indmj
  integer(i4)                                  :: indmk
  integer(i4)                                  :: indml
  integer(i4)                                  :: indmm
  integer(i4)                                  :: indmn
  integer(i4)                                  :: indvec
  integer(i4)                                  :: isatom(6)
  integer(i4)                                  :: isxyz(3,6)
  integer(i4)                                  :: isxyzl(3)
  integer(i4)                                  :: ivec
  integer(i4)                                  :: iveca(2,15)
  integer(i4)                                  :: ixm
  integer(i4)                                  :: ixn
  integer(i4)                                  :: iym
  integer(i4)                                  :: iyn
  integer(i4)                                  :: izm
  integer(i4)                                  :: izn
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
  integer(i4)                                  :: jb
  integer(i4)                                  :: jj
  integer(i4)                                  :: jmax
  integer(i4)                                  :: jmin
  integer(i4)                                  :: jvec
  integer(i4)                                  :: jxx
  integer(i4)                                  :: jyy
  integer(i4)                                  :: jzz
  integer(i4)                                  :: k
  integer(i4)                                  :: kb2(6,6)
  integer(i4)                                  :: kloc
  integer(i4)                                  :: kxx
  integer(i4)                                  :: kyy
  integer(i4)                                  :: kzz
  integer(i4)                                  :: l
  integer(i4)                                  :: liimax
  integer(i4)                                  :: ll
  integer(i4)                                  :: lmax
  integer(i4)                                  :: lmin
  integer(i4)                                  :: m
  integer(i4)                                  :: mm
  integer(i4)                                  :: mn
  integer(i4)                                  :: mxx
  integer(i4)                                  :: myy
  integer(i4)                                  :: mzz
  integer(i4)                                  :: n
  integer(i4)                                  :: nbtypeij
  integer(i4)                                  :: nbtypeik
  integer(i4)                                  :: nbtypeil
  integer(i4)                                  :: nbtypejm
  integer(i4)                                  :: nbtypejn
  integer(i4)                                  :: nbtypeij2
  integer(i4)                                  :: nbtypeik2
  integer(i4)                                  :: nbtypeil2
  integer(i4)                                  :: nbtypejm2
  integer(i4)                                  :: nbtypejn2
  integer(i4)                                  :: nsixtype
  integer(i4)                                  :: ni
  integer(i4)                                  :: niimax
  integer(i4)                                  :: nj
  integer(i4)                                  :: nk
  integer(i4)                                  :: nl 
  integer(i4)                                  :: nm 
  integer(i4)                                  :: nn 
  integer(i4)                                  :: nmi 
  integer(i4)                                  :: nmj 
  integer(i4)                                  :: nmk 
  integer(i4)                                  :: nml 
  integer(i4)                                  :: nmm 
  integer(i4)                                  :: nmn 
  integer(i4)                                  :: nmax
  integer(i4)                                  :: nmin
  integer(i4)                                  :: np
  integer(i4)                                  :: nt2
  integer(i4)                                  :: nt3
  integer(i4)                                  :: nt4
  integer(i4)                                  :: nt5
  integer(i4)                                  :: nt6
  integer(i4)                                  :: ntyp2
  integer(i4)                                  :: ntyp3
  integer(i4)                                  :: ntyp4
  integer(i4)                                  :: ntyp5
  integer(i4)                                  :: ntyp6
  integer(i4)                                  :: ntypi
  integer(i4)                                  :: ntypj
  integer(i4)                                  :: ntypk
  integer(i4)                                  :: ntypl
  integer(i4)                                  :: ntypm
  integer(i4)                                  :: ntypn
  integer(i4)                                  :: nxx
  integer(i4)                                  :: nyy
  integer(i4)                                  :: nzz
  logical                                      :: l2bonds
  logical                                      :: lbonded
  logical                                      :: lbtyp
  logical                                      :: linter_only
  logical                                      :: lintra_only
  logical                                      :: lmatch
  logical                                      :: lmolok
  logical                                      :: lmolokij
  logical                                      :: lmolokik
  logical                                      :: lmolokil
  logical                                      :: lmolokjm
  logical                                      :: lmolokjn
  logical                                      :: lneedmol
  logical                                      :: lsamemol
  logical                                      :: ltsym12
  logical                                      :: ltsym34
  logical                                      :: ltsym56
  complex(dpc)                                 :: cdk(3)
  real(dp)                                     :: cosvec(15)
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: d2k
  real(dp)                                     :: d2ks
  real(dp)                                     :: d2trm
  real(dp)                                     :: e1d(15)
  real(dp)                                     :: e2d(120)
  real(dp)                                     :: e3d(1)
  real(dp)                                     :: eterm
  real(dp)                                     :: oci
  real(dp)                                     :: ocj
  real(dp)                                     :: ock
  real(dp)                                     :: ocl 
  real(dp)                                     :: ocm 
  real(dp)                                     :: ocn 
  real(dp)                                     :: ofct 
  real(dp)                                     :: phase
  real(dp)                                     :: r212
  real(dp)                                     :: r312
  real(dp)                                     :: r412
  real(dp)                                     :: r522
  real(dp)                                     :: r622
  real(dp)                                     :: rksix
  real(dp)                                     :: rko
  real(dp)                                     :: sinvec(15)
  real(dp)                                     :: time1
  real(dp)                                     :: time2
  real(dp)                                     :: tr1
  real(dp)                                     :: tr2
  real(dp)                                     :: tr3
  real(dp)                                     :: tr4
  real(dp)                                     :: tr5
  real(dp)                                     :: ttr2
  real(dp)                                     :: ttr3
  real(dp)                                     :: ttr4
  real(dp)                                     :: ttr5
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
  real(dp)                                     :: x52
  real(dp)                                     :: y52
  real(dp)                                     :: z52
  real(dp)                                     :: x52t
  real(dp)                                     :: y52t
  real(dp)                                     :: z52t
  real(dp)                                     :: x62
  real(dp)                                     :: y62
  real(dp)                                     :: z62
  real(dp)                                     :: x62t
  real(dp)                                     :: y62t
  real(dp)                                     :: z62t
  real(dp)                                     :: xc2t
  real(dp)                                     :: yc2t
  real(dp)                                     :: zc2t
  real(dp)                                     :: xc3t
  real(dp)                                     :: yc3t
  real(dp)                                     :: zc3t
  real(dp)                                     :: xc4t
  real(dp)                                     :: yc4t
  real(dp)                                     :: zc4t
  real(dp)                                     :: xc5t
  real(dp)                                     :: yc5t
  real(dp)                                     :: zc5t
  real(dp)                                     :: xc6t
  real(dp)                                     :: yc6t
  real(dp)                                     :: zc6t
  real(dp)                                     :: xvec(15)
  real(dp)                                     :: yvec(15)
  real(dp)                                     :: zvec(15)
  real(dp)                                     :: sdist(15)
  real(dp)                                     :: svec(3,6,6)
  real(dp)                                     :: sxyz(3,6)
#ifdef TRACE
  call trace_in('sixpd')
#endif
!
  time1 = g_cpu_time()
!
!  Set up vector to atom mapping
!
  indvec = 0
  do ivec = 1,5
    do jvec = ivec+1,6
      indvec = indvec + 1
      iveca(1,indvec) = jvec
      iveca(2,indvec) = ivec
    enddo
  enddo
  indvec = 0
  do ivec = 1,5
    do jvec = ivec+1,6
      indvec = indvec + 1
      kb2(jvec,ivec) = indvec
      kb2(ivec,jvec) = indvec
    enddo
  enddo
!*************************
!  Loop over potentials  *
!*************************
  do np = 1,nsix
    nsixtype = nsixty(np)
    tr1 = six1(np)**2
    ttr2 = six2(np)**2
    ttr3 = six3(np)**2
    ttr4 = six4(np)**2
    ttr5 = six5(np)**2
    ltsym12 = (lmatch(nsspec1(np),nsptyp1(np),nsspec2(np),nsptyp2(np),.true.).or. &
               lmatch(nsspec2(np),nsptyp2(np),nsspec1(np),nsptyp1(np),.true.))
    lbtyp = (mmsexc(np).eq.1)
    rksix = sixk(np)
    lintra_only = (lsintra(np).and..not.lsinter(np))
    linter_only = (lsinter(np).and..not.lsintra(np))
    lneedmol = (lintra_only.or.linter_only.or.lbtyp)
!********************************
!  Loop over middle site 1 / i  *
!********************************
    iloop: do iloc = 1,natomsonnode
      i = node2atom(iloc)
      ni = nat(i)
      ntypi = nftype(i)
      oci = occuf(i)
!
!  Check i is allowed for n
!
      if (lmatch(ni,ntypi,nsspec1(np),nsptyp1(np),.true.)) then
        nt2 = nsspec2(np)
        ntyp2 = nsptyp2(np)
        nt3 = nsspec3(np)
        nt4 = nsspec4(np)
        nt5 = nsspec5(np)
        nt6 = nsspec6(np)
        ntyp3 = nsptyp3(np)
        ntyp4 = nsptyp4(np)
        ntyp5 = nsptyp5(np)
        ntyp6 = nsptyp6(np)
        tr2 = ttr2
        tr3 = ttr3
        tr4 = ttr4
        tr5 = ttr5
      elseif (lmatch(ni,ntypi,nsspec2(np),nsptyp2(np),.true.)) then
        nt2 = nsspec1(np)
        ntyp2 = nsptyp1(np)
        nt3 = nsspec5(np)
        nt4 = nsspec6(np)
        nt5 = nsspec3(np)
        nt6 = nsspec4(np)
        ntyp3 = nsptyp5(np)
        ntyp4 = nsptyp6(np)
        ntyp5 = nsptyp3(np)
        ntyp6 = nsptyp4(np)
        tr2 = ttr4
        tr3 = ttr5
        tr4 = ttr2
        tr5 = ttr3
      else
        cycle iloop
      endif
!
!  Set flags as to whether potentials are symmetric at the two ends
!
      ltsym34 = (lmatch(nt3,ntyp3,nt4,ntyp4,.true.).or.lmatch(nt4,ntyp4,nt3,ntyp3,.true.))
      ltsym56 = (lmatch(nt5,ntyp5,nt6,ntyp6,.true.).or.lmatch(nt6,ntyp6,nt5,ntyp5,.true.))
!
!  i has been accepted
!
      isatom(1) = i
      isxyz(1,1) = 3*(i - 1) + 1
      isxyz(2,1) = 3*(i - 1) + 2
      isxyz(3,1) = 3*(i - 1) + 3
!
      isxyzl(1) = 3*(iloc - 1) + 1
      isxyzl(2) = 3*(iloc - 1) + 2
      isxyzl(3) = 3*(iloc - 1) + 3
!
      sxyz(1,1) = xclat(i)
      sxyz(2,1) = yclat(i)
      sxyz(3,1) = zclat(i)
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
!
!  Loop over all j in this algorithm
!
      jmin = 1
      jmax = numat
!********************************
!  Loop over middle site 2 / j  *
!********************************
      jloop: do j = jmin,jmax
        nj = nat(j)
        ntypj = nftype(j)
        ocj = occuf(j)
!
!  Check j is allowed for n
!
        if (.not.lmatch(nj,ntypj,nt2,ntyp2,.true.)) cycle jloop
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
        if (lintra_only.and..not.lmolok) cycle jloop
        if (lbtyp.and..not.lmolok) cycle jloop
        xc2t = xclat(j)
        yc2t = yclat(j)
        zc2t = zclat(j)
        x21t = xc2t - sxyz(1,1)
        y21t = yc2t - sxyz(2,1)
        z21t = zc2t - sxyz(3,1)
!
        lmolokij = lmolok
!
!  Check r21 is OK
!  Loop over cell vectors
!
        iiloop: do ii = 1,iimax
          x21 = x21t + xvec1cell(ii)
          y21 = y21t + yvec1cell(ii)
          z21 = z21t + zvec1cell(ii)
          r212 = x21*x21 + y21*y21 + z21*z21
          if (r212.lt.1d-10) cycle iiloop
!
!  Reset lmolok to preloop value
!
          lmolok = lmolokij
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
!
!  Check central bond type for correct order
!
                if (n6botype(1,np).gt.0) then
                  if (n6botype(1,np).ne.nbtypeij) cycle iiloop
                endif
                if (n6botype(2,np).gt.0) then
                  if (n6botype(2,np).ne.nbtypeij2) cycle iiloop
                endif
              endif
            else
              call lintoijk(ixx,iyy,izz,ii,imaxl,jmaxl,kmaxl)
              if (lbtyp) then
                call bonded(lbonded,l2bonds,nbtypeij,nbtypeij2,i,j,ixx,iyy,izz)
                if (.not.lbonded) cycle iiloop
!
!  Check central bond type for correct order
!
                if (n6botype(1,np).gt.0) then
                  if (n6botype(1,np).ne.nbtypeij) cycle iiloop
                endif
                if (n6botype(2,np).gt.0) then
                  if (n6botype(2,np).ne.nbtypeij2) cycle iiloop
                endif
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
          if (r212.gt.tr1.and.(.not.lbtyp.or..not.lbonded)) cycle iiloop
!
!  j has been accepted
!
          isatom(2) = j
          isxyz(1,2) = 3*(j - 1) + 1
          isxyz(2,2) = 3*(j - 1) + 2
          isxyz(3,2) = 3*(j - 1) + 3
          sxyz(1,2) = x21 + sxyz(1,1)
          sxyz(2,2) = y21 + sxyz(2,1)
          sxyz(3,2) = z21 + sxyz(3,1)
!*****************************************
!  Loop over first end site bonded to i  *
!*****************************************
          kloop: do k = 1,numat
            nk = nat(k)
            ntypk = nftype(k)
            ock = occuf(k)
!
!  Check k is allowed for n
!
            if (.not.lmatch(nk,ntypk,nt3,ntyp3,.true.)) cycle kloop
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
            if (lintra_only.and..not.lmolok) cycle kloop
            if (lbtyp.and..not.lmolok) cycle kloop
            xc3t = xclat(k)
            yc3t = yclat(k)
            zc3t = zclat(k)
            x31t = xc3t - sxyz(1,1)
            y31t = yc3t - sxyz(2,1)
            z31t = zc3t - sxyz(3,1)
!
            lmolokik = lmolok
!
!  Check r31 is OK
!  Loop over cell vectors
!
            jjloop: do jj = 1,iimax
              x31 = x31t + xvec1cell(jj)
              y31 = y31t + yvec1cell(jj)
              z31 = z31t + zvec1cell(jj)
              r312 = x31*x31 + y31*y31 + z31*z31
              if (r312.lt.1d-10) cycle jjloop
!
!  Prevent atoms i and k being the same atom
!
              if (k.eq.i.and.jj.eq.iimid) cycle jjloop
!
!  Prevent atoms j and k being the same atom
!
              if (k.eq.j.and.jj.eq.ii) cycle jjloop
!
!  Reset lmolok to preloop value
!
              lmolok = lmolokik
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
              if (r312.gt.tr2.and.(.not.lbtyp.or..not.lbonded)) cycle jjloop
!
!  k has been accepted
!
              isatom(3) = k
              isxyz(1,3) = 3*(k - 1) + 1
              isxyz(2,3) = 3*(k - 1) + 2
              isxyz(3,3) = 3*(k - 1) + 3
              sxyz(1,3) = x31 + sxyz(1,1)
              sxyz(2,3) = y31 + sxyz(2,1)
              sxyz(3,3) = z31 + sxyz(3,1)
!************************************
!  Loop over second end site for i  *
!************************************
!
!  Set l looping indices
!
              if (ltsym34) then
                lmin = 1
                lmax = k
              else
                lmin = 1
                lmax = numat
              endif
              lloop: do l = lmin,lmax
                nl = nat(l)
                ntypl = nftype(l)
                ocl = occuf(l)
!
!  Check l is allowed for n
!
                if (.not.lmatch(nl,ntypl,nt4,ntyp4,.true.)) cycle lloop
!
!  Molecule handling
!
                if (lmol.and.lneedmol) then
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
                if (lintra_only.and..not.lmolok) cycle lloop
                if (lbtyp.and..not.lmolok) cycle lloop
                xc4t = xclat(l)
                yc4t = yclat(l)
                zc4t = zclat(l)
                x41t = xc4t - sxyz(1,1)
                y41t = yc4t - sxyz(2,1)
                z41t = zc4t - sxyz(3,1)
                if (ltsym34.and.k.eq.l) then
                  liimax = jj - 1
                else
                  liimax = iimax
                endif
!
                lmolokil = lmolok
!
!  Check r41 is OK
!  Loop over cell vectors
!
                llloop: do ll = 1,liimax
                  x41 = x41t + xvec1cell(ll)
                  y41 = y41t + yvec1cell(ll)
                  z41 = z41t + zvec1cell(ll)
                  r412 = x41*x41 + y41*y41 + z41*z41
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
!  Reset lmolok to preloop value
!
                  lmolok = lmolokil
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
                  if (r412.gt.tr3.and.(.not.lbtyp.or..not.lbonded)) cycle llloop
!
!  l has been accepted
!
                  isatom(4) = l
                  isxyz(1,4) = 3*(l - 1) + 1
                  isxyz(2,4) = 3*(l - 1) + 2
                  isxyz(3,4) = 3*(l - 1) + 3
                  sxyz(1,4) = x41 + sxyz(1,1)
                  sxyz(2,4) = y41 + sxyz(2,1)
                  sxyz(3,4) = z41 + sxyz(3,1)
!*****************************************
!  Loop over first end site bonded to j  *
!*****************************************
                  mloop: do m = 1,numat
                    nm = nat(m)
                    ntypm = nftype(m)
                    ocm = occuf(m)
!
!  Check m is allowed for n
!
                    if (.not.lmatch(nm,ntypm,nt5,ntyp5,.true.)) cycle mloop
                    if (lmol.and.lneedmol) then
!
!  Molecule handling
!
                      nmm = natmol(m)
                      if (ndim.gt.0) then
                        indmm = nmolind(m)
                        call mindtoijk(indmm,ixm,iym,izm)
                        ixm = ixm - ixj
                        iym = iym - iyj
                        izm = izm - izj
                      endif
                      lmolok = (nmi.eq.nmm.and.nmi.gt.0)
                    else
                      lmolok = .false.
                    endif
!
!  Check for intra and but not in same molecule
!
                    if (lintra_only.and..not.lmolok) cycle mloop
                    if (lbtyp.and..not.lmolok) cycle mloop
                    xc5t = xclat(m)
                    yc5t = yclat(m)
                    zc5t = zclat(m)
                    x52t = xc5t - sxyz(1,2)
                    y52t = yc5t - sxyz(2,2)
                    z52t = zc5t - sxyz(3,2)
!
                    lmolokjm = lmolok
!
!  Check r52 is OK
!  Loop over cell vectors
!
                    mmloop: do mm = 1,iimax
                      x52 = x52t + xvec1cell(mm)
                      y52 = y52t + yvec1cell(mm)
                      z52 = z52t + zvec1cell(mm)
                      r522 = x52*x52 + y52*y52 + z52*z52
                      if (r522.lt.1d-10) cycle mmloop
!
!  Prevent atoms i and m being the same atom
!
                      if (m.eq.i.and.mm.eq.iimid) cycle mmloop
!
!  Prevent atoms j and m being the same atom
!
                      if (m.eq.j.and.mm.eq.ii) cycle mmloop
!
!  Reset lmolok to preloop value
!
                      lmolok = lmolokjm
!
!  Molecule checking
!
                      lbonded = .false.
                      if (lmolok) then
                        if (ndim.eq.0) then
                          if (linter_only) cycle mmloop
                          if (lbtyp) then
                            call bonded(lbonded,l2bonds,nbtypejm,nbtypejm2,j,m,0_i4,0_i4,0_i4)
                            if (.not.lbonded) cycle mmloop
                          endif
                        else
                          call lintoijk(mxx,myy,mzz,mm,imaxl,jmaxl,kmaxl)
                          if (lbtyp) then
                            call bonded(lbonded,l2bonds,nbtypejm,nbtypejm2,j,m,mxx-ixx,myy-iyy,mzz-izz)
                            if (.not.lbonded) cycle mmloop
                            lsamemol = (lbonded.or.l2bonds)
                          else
                            lsamemol = .false.
                          endif
                          if (.not.lsamemol) then
                            call samemol(lsamemol,nmi,jxx,jyy,jzz,ixm,iym,izm)
                          endif
                          if (lintra_only.and..not.lsamemol) cycle mmloop
                          if (linter_only.and.lsamemol) cycle mmloop
                        endif
                      endif
!
!  Distance checking
!
                      if (r522.gt.tr4.and.(.not.lbtyp.or..not.lbonded)) cycle mmloop
!
!  m has been accepted
!
                      isatom(5) = m
                      isxyz(1,5) = 3*(m - 1) + 1
                      isxyz(2,5) = 3*(m - 1) + 2
                      isxyz(3,5) = 3*(m - 1) + 3
                      sxyz(1,5) = x52 + sxyz(1,2)
                      sxyz(2,5) = y52 + sxyz(2,2)
                      sxyz(3,5) = z52 + sxyz(3,2)
!************************************
!  Loop over second end site for j  *
!************************************
!
!  Set n looping indices
!
                      if (ltsym56) then
                        nmin = 1
                        nmax = m
                      else
                        nmin = 1
                        nmax = numat
                      endif
                      nloop: do n = nmin,nmax
                        nn = nat(n)
                        ntypn = nftype(n)
                        ocn = occuf(n)
!
!  Check n is allowed for n
!
                        if (.not.lmatch(nn,ntypn,nt6,ntyp6,.true.)) cycle nloop
!
!  Molecule handling
!
                        if (lmol.and.lneedmol) then
                          nmn = natmol(n)
                          if (ndim.gt.0) then
                            indmn = nmolind(n)
                            call mindtoijk(indmn,ixn,iyn,izn)
                            ixn = ixn - ixj
                            iyn = iyn - iyj
                            izn = izn - izj
                          endif
                          lmolok = (nmi.eq.nmn.and.nmi.gt.0)
                        else
                          lmolok = .false.
                        endif
!
!  Check for intra and but not in same molecule
!
                        if (lintra_only.and..not.lmolok) cycle nloop
                        if (lbtyp.and..not.lmolok) cycle nloop
                        xc6t = xclat(n)
                        yc6t = yclat(n)
                        zc6t = zclat(n)
                        x62t = xc6t - sxyz(1,2)
                        y62t = yc6t - sxyz(2,2)
                        z62t = zc6t - sxyz(3,2)
                        if (ltsym56.and.m.eq.n) then
                          niimax = mm - 1
                        else
                          niimax = iimax
                        endif
!
                        lmolokjn = lmolok
!
!  Check r62 is OK
!  Loop over cell vectors
!
                        mnloop: do mn = 1,niimax
                          x62 = x62t + xvec1cell(mn)
                          y62 = y62t + yvec1cell(mn)
                          z62 = z62t + zvec1cell(mn)
                          r622 = x62*x62 + y62*y62 + z62*z62
                          if (r622.lt.1d-10) cycle mnloop
!
!  Prevent atoms i and n being the same atom
!
                          if (n.eq.i.and.mn.eq.iimid) cycle mnloop
!
!  Prevent atoms j and n being the same atom
!
                          if (n.eq.j.and.mn.eq.ii) cycle mnloop
!
!  Prevent atoms m and n being the same atom
!
                          if (n.eq.m.and.mn.eq.mm) cycle mnloop
!
!  Reset lmolok to preloop value
!
                          lmolok = lmolokjn
!
!  Molecule checking
!
                          lbonded = .false.
                          if (lmolok) then
                            if (ndim.eq.0) then
                              if (linter_only) cycle mnloop
                              if (lbtyp) then
                                call bonded(lbonded,l2bonds,nbtypejn,nbtypejn2,j,n,0_i4,0_i4,0_i4)
                                if (.not.lbonded) cycle mnloop
                              endif
                            else
                              call lintoijk(nxx,nyy,nzz,mn,imaxl,jmaxl,kmaxl)
                              if (lbtyp) then
                                call bonded(lbonded,l2bonds,nbtypejn,nbtypejn2,j,n,nxx-ixx,nyy-iyy,nzz-izz)
                                if (.not.lbonded) cycle mnloop
                                lsamemol = (lbonded.or.l2bonds)
                              else
                                lsamemol = .false.
                              endif
                              if (.not.lsamemol) then
                                call samemol(lsamemol,nmi,jxx,jyy,jzz,ixn,iyn,izn)
                              endif
                              if (lintra_only.and..not.lsamemol) cycle mnloop
                              if (linter_only.and.lsamemol) cycle mnloop
                            endif
                          endif
!
!  Distance checking
!
                          if (r622.gt.tr5.and.(.not.lbtyp.or..not.lbonded)) cycle mnloop
!
!  n has been accepted
!
                          isatom(6) = n
                          isxyz(1,6) = 3*(n - 1) + 1
                          isxyz(2,6) = 3*(n - 1) + 2
                          isxyz(3,6) = 3*(n - 1) + 3
                          sxyz(1,6) = x62 + sxyz(1,2)
                          sxyz(2,6) = y62 + sxyz(2,2)
                          sxyz(3,6) = z62 + sxyz(3,2)
!********************************
!  Valid six-body term located  *
!********************************
!
!  Calculate vectors and remaining distances
!
                          do ivec = 1,6
                            do jvec = 1,6
                              svec(1,jvec,ivec) = sxyz(1,ivec) - sxyz(1,jvec)
                              svec(2,jvec,ivec) = sxyz(2,ivec) - sxyz(2,jvec)
                              svec(3,jvec,ivec) = sxyz(3,ivec) - sxyz(3,jvec)
                            enddo
                          enddo
                          indvec = 0
                          do ivec = 1,5
                            do jvec = ivec+1,6
                              indvec = indvec + 1
                              sdist(indvec) = svec(1,jvec,ivec)**2 + svec(2,jvec,ivec)**2 + svec(3,jvec,ivec)**2
                              sdist(indvec) = sqrt(sdist(indvec))
                            enddo
                          enddo
!
                          ofct = oci*ocj*ock*ocl*ocm*ocn
                          rko = rksix*ofct
                          call sixbody(np,nsixtype,sdist,eterm,e1d,e2d,e3d,rko,.true.,.true.,.false.)
!*************************
!  Internal derivatives  *
!*************************
!
!  Compute phase factors
!
                          indvec = 0
                          do ivec = 1,5
                            do jvec = ivec+1,6
                              indvec = indvec + 1
                              phase = svec(1,jvec,ivec)*xkv + svec(2,jvec,ivec)*ykv + svec(3,jvec,ivec)*zkv
                              sinvec(indvec) = sin(phase)
                              if (isatom(ivec).eq.isatom(jvec)) then
                                cosvec(indvec) = cos(phase) - 1.0_dp
                              else
                                cosvec(indvec) = cos(phase) 
                              endif
                              xvec(indvec) = svec(1,jvec,ivec)
                              yvec(indvec) = svec(2,jvec,ivec)
                              zvec(indvec) = svec(3,jvec,ivec)
                            enddo
                          enddo
!
                          indvec = 0
                          do ivec = 1,15
                            ja = iveca(1,ivec)
                            ia = iveca(2,ivec)
!
                            if (ia.eq.1) then
                              derv2(isxyz(1,ja),isxyzl(1)) = derv2(isxyz(1,ja),isxyzl(1)) - e1d(ivec)*cosvec(ivec)
                              derv2(isxyz(2,ja),isxyzl(2)) = derv2(isxyz(2,ja),isxyzl(2)) - e1d(ivec)*cosvec(ivec)
                              derv2(isxyz(3,ja),isxyzl(3)) = derv2(isxyz(3,ja),isxyzl(3)) - e1d(ivec)*cosvec(ivec)
                              if (ia.gt.ja) then
                                dervi(isxyz(1,ja),isxyzl(1)) = dervi(isxyz(1,ja),isxyzl(1)) - e1d(ivec)*sinvec(ivec)
                                dervi(isxyz(2,ja),isxyzl(2)) = dervi(isxyz(2,ja),isxyzl(2)) - e1d(ivec)*sinvec(ivec)
                                dervi(isxyz(3,ja),isxyzl(3)) = dervi(isxyz(3,ja),isxyzl(3)) - e1d(ivec)*sinvec(ivec)
                              else
                                dervi(isxyz(1,ja),isxyzl(1)) = dervi(isxyz(1,ja),isxyzl(1)) + e1d(ivec)*sinvec(ivec)
                                dervi(isxyz(2,ja),isxyzl(2)) = dervi(isxyz(2,ja),isxyzl(2)) + e1d(ivec)*sinvec(ivec)
                                dervi(isxyz(3,ja),isxyzl(3)) = dervi(isxyz(3,ja),isxyzl(3)) + e1d(ivec)*sinvec(ivec)
                              endif
                              if (lgroupvelocity) then
!
!  Group velocities
!
                                d2k  = e1d(ivec)*cosvec(ivec)
                                d2ks = e1d(ivec)*sinvec(ivec)
                                cdk(1) = dcmplx(d2k*xvec(ivec),d2ks*xvec(ivec))*dcmplx(0.0_dp,1.0_dp)
                                cdk(2) = dcmplx(d2k*yvec(ivec),d2ks*yvec(ivec))*dcmplx(0.0_dp,1.0_dp)
                                cdk(3) = dcmplx(d2k*zvec(ivec),d2ks*zvec(ivec))*dcmplx(0.0_dp,1.0_dp)
!
                                if (ia.gt.ja) then
                                  derv2dk(1:3,isxyz(1,ja),isxyzl(1)) = derv2dk(1:3,isxyz(1,ja),isxyzl(1)) - cdk(1:3)
                                  derv2dk(1:3,isxyz(2,ja),isxyzl(2)) = derv2dk(1:3,isxyz(2,ja),isxyzl(2)) - cdk(1:3)
                                  derv2dk(1:3,isxyz(3,ja),isxyzl(3)) = derv2dk(1:3,isxyz(3,ja),isxyzl(3)) - cdk(1:3)
                                else
                                  derv2dk(1:3,isxyz(1,ja),isxyzl(1)) = derv2dk(1:3,isxyz(1,ja),isxyzl(1)) - conjg(cdk(1:3))
                                  derv2dk(1:3,isxyz(2,ja),isxyzl(2)) = derv2dk(1:3,isxyz(2,ja),isxyzl(2)) - conjg(cdk(1:3))
                                  derv2dk(1:3,isxyz(3,ja),isxyzl(3)) = derv2dk(1:3,isxyz(3,ja),isxyzl(3)) - conjg(cdk(1:3))
                                endif
                              endif
                            endif
!
                            if (ja.eq.1) then
                              derv2(isxyz(1,ia),isxyzl(1)) = derv2(isxyz(1,ia),isxyzl(1)) - e1d(ivec)*cosvec(ivec)
                              derv2(isxyz(2,ia),isxyzl(2)) = derv2(isxyz(2,ia),isxyzl(2)) - e1d(ivec)*cosvec(ivec)
                              derv2(isxyz(3,ia),isxyzl(3)) = derv2(isxyz(3,ia),isxyzl(3)) - e1d(ivec)*cosvec(ivec)
                              if (ia.gt.ja) then
                                dervi(isxyz(1,ia),isxyzl(1)) = dervi(isxyz(1,ia),isxyzl(1)) + e1d(ivec)*sinvec(ivec)
                                dervi(isxyz(2,ia),isxyzl(2)) = dervi(isxyz(2,ia),isxyzl(2)) + e1d(ivec)*sinvec(ivec)
                                dervi(isxyz(3,ia),isxyzl(3)) = dervi(isxyz(3,ia),isxyzl(3)) + e1d(ivec)*sinvec(ivec)
                              else
                                dervi(isxyz(1,ia),isxyzl(1)) = dervi(isxyz(1,ia),isxyzl(1)) - e1d(ivec)*sinvec(ivec)
                                dervi(isxyz(2,ia),isxyzl(2)) = dervi(isxyz(2,ia),isxyzl(2)) - e1d(ivec)*sinvec(ivec)
                                dervi(isxyz(3,ia),isxyzl(3)) = dervi(isxyz(3,ia),isxyzl(3)) - e1d(ivec)*sinvec(ivec)
                              endif
                              if (lgroupvelocity) then
!
!  Group velocities
!
                                d2k  = e1d(ivec)*cosvec(ivec)
                                d2ks = e1d(ivec)*sinvec(ivec)
                                cdk(1) = dcmplx(d2k*xvec(ivec),d2ks*xvec(ivec))*dcmplx(0.0_dp,1.0_dp)
                                cdk(2) = dcmplx(d2k*yvec(ivec),d2ks*yvec(ivec))*dcmplx(0.0_dp,1.0_dp)
                                cdk(3) = dcmplx(d2k*zvec(ivec),d2ks*zvec(ivec))*dcmplx(0.0_dp,1.0_dp)
!
                                if (ia.gt.ja) then
                                  derv2dk(1:3,isxyz(1,ia),isxyzl(1)) = derv2dk(1:3,isxyz(1,ia),isxyzl(1)) - conjg(cdk(1:3))
                                  derv2dk(1:3,isxyz(2,ia),isxyzl(2)) = derv2dk(1:3,isxyz(2,ia),isxyzl(2)) - conjg(cdk(1:3))
                                  derv2dk(1:3,isxyz(3,ia),isxyzl(3)) = derv2dk(1:3,isxyz(3,ia),isxyzl(3)) - conjg(cdk(1:3))
                                else
                                  derv2dk(1:3,isxyz(1,ia),isxyzl(1)) = derv2dk(1:3,isxyz(1,ia),isxyzl(1)) - cdk(1:3)
                                  derv2dk(1:3,isxyz(2,ia),isxyzl(2)) = derv2dk(1:3,isxyz(2,ia),isxyzl(2)) - cdk(1:3)
                                  derv2dk(1:3,isxyz(3,ia),isxyzl(3)) = derv2dk(1:3,isxyz(3,ia),isxyzl(3)) - cdk(1:3)
                                endif
                              endif
                            endif
!
                            do jvec = ivec,15
                              indvec = indvec + 1 
                              jb = iveca(1,jvec)
                              ib = iveca(2,jvec)
                              do ialp = 1,3
                                do ibet = 1,3
                                  d2trm = svec(ialp,ja,ia)*svec(ibet,jb,ib)*e2d(indvec)
                                  if (ia.eq.1) then
                                    if (ia.gt.ib) then
                                      derv2(isxyz(ibet,ib),isxyzl(ialp)) = derv2(isxyz(ibet,ib),isxyzl(ialp)) +  &
                                          d2trm*cosvec(kb2(ib,ia))
                                      dervi(isxyz(ibet,ib),isxyzl(ialp)) = dervi(isxyz(ibet,ib),isxyzl(ialp)) +  &
                                          d2trm*sinvec(kb2(ib,ia))
                                    elseif (ia.lt.ib) then
                                      derv2(isxyz(ibet,ib),isxyzl(ialp)) = derv2(isxyz(ibet,ib),isxyzl(ialp)) +  &
                                          d2trm*cosvec(kb2(ib,ia))
                                      dervi(isxyz(ibet,ib),isxyzl(ialp)) = dervi(isxyz(ibet,ib),isxyzl(ialp)) -  &
                                          d2trm*sinvec(kb2(ib,ia))
                                    endif
!
                                    if (lgroupvelocity.and.ia.ne.ib) then
!
!  Group velocities
!
                                      d2k  = d2trm*cosvec(kb2(ib,ia))
                                      d2ks = d2trm*sinvec(kb2(ib,ia))
                                      cdk(1) = dcmplx(d2k*xvec(kb2(ib,ia)),d2ks*xvec(kb2(ib,ia)))*dcmplx(0.0_dp,1.0_dp)
                                      cdk(2) = dcmplx(d2k*yvec(kb2(ib,ia)),d2ks*yvec(kb2(ib,ia)))*dcmplx(0.0_dp,1.0_dp)
                                      cdk(3) = dcmplx(d2k*zvec(kb2(ib,ia)),d2ks*zvec(kb2(ib,ia)))*dcmplx(0.0_dp,1.0_dp)
!
                                      if (ia.gt.ib) then
                                        derv2dk(1:3,isxyz(ibet,ib),isxyzl(ialp)) = derv2dk(1:3,isxyz(ibet,ib),isxyzl(ialp)) +  &
                                            cdk(1:3)
                                      elseif (ia.lt.ib) then
                                        derv2dk(1:3,isxyz(ibet,ib),isxyzl(ialp)) = derv2dk(1:3,isxyz(ibet,ib),isxyzl(ialp)) +  &
                                            conjg(cdk(1:3))
                                      endif
                                    endif
!
                                    if (ia.gt.jb) then
                                      derv2(isxyz(ibet,jb),isxyzl(ialp)) = derv2(isxyz(ibet,jb),isxyzl(ialp)) -  &
                                          d2trm*cosvec(kb2(jb,ia))
                                      dervi(isxyz(ibet,jb),isxyzl(ialp)) = dervi(isxyz(ibet,jb),isxyzl(ialp)) -  &
                                          d2trm*sinvec(kb2(jb,ia))
                                    elseif (ia.lt.jb) then
                                      derv2(isxyz(ibet,jb),isxyzl(ialp)) = derv2(isxyz(ibet,jb),isxyzl(ialp)) -  &
                                          d2trm*cosvec(kb2(jb,ia))
                                      dervi(isxyz(ibet,jb),isxyzl(ialp)) = dervi(isxyz(ibet,jb),isxyzl(ialp)) +  &
                                          d2trm*sinvec(kb2(jb,ia))
                                    endif
!
                                    if (lgroupvelocity.and.ia.ne.jb) then
!
!  Group velocities
!
                                      d2k  = d2trm*cosvec(kb2(jb,ia))
                                      d2ks = d2trm*sinvec(kb2(jb,ia))
                                      cdk(1) = dcmplx(d2k*xvec(kb2(jb,ia)),d2ks*xvec(kb2(jb,ia)))*dcmplx(0.0_dp,1.0_dp)
                                      cdk(2) = dcmplx(d2k*yvec(kb2(jb,ia)),d2ks*yvec(kb2(jb,ia)))*dcmplx(0.0_dp,1.0_dp)
                                      cdk(3) = dcmplx(d2k*zvec(kb2(jb,ia)),d2ks*zvec(kb2(jb,ia)))*dcmplx(0.0_dp,1.0_dp)
!
                                      if (ia.gt.jb) then
                                        derv2dk(1:3,isxyz(ibet,jb),isxyzl(ialp)) = derv2dk(1:3,isxyz(ibet,jb),isxyzl(ialp)) -  &
                                            cdk(1:3)
                                      elseif (ia.lt.jb) then
                                        derv2dk(1:3,isxyz(ibet,jb),isxyzl(ialp)) = derv2dk(1:3,isxyz(ibet,jb),isxyzl(ialp)) -  &
                                            conjg(cdk(1:3))
                                      endif
                                    endif   
                                  endif   
!
                                  if (ja.eq.1) then
                                    if (ja.gt.jb) then
                                      derv2(isxyz(ibet,jb),isxyzl(ialp)) = derv2(isxyz(ibet,jb),isxyzl(ialp)) +  &
                                          d2trm*cosvec(kb2(jb,ja))
                                      dervi(isxyz(ibet,jb),isxyzl(ialp)) = dervi(isxyz(ibet,jb),isxyzl(ialp)) +  &
                                          d2trm*sinvec(kb2(jb,ja))
                                    elseif (ja.lt.jb) then
                                      derv2(isxyz(ibet,jb),isxyzl(ialp)) = derv2(isxyz(ibet,jb),isxyzl(ialp)) +  &
                                          d2trm*cosvec(kb2(jb,ja))
                                      dervi(isxyz(ibet,jb),isxyzl(ialp)) = dervi(isxyz(ibet,jb),isxyzl(ialp)) -  &
                                          d2trm*sinvec(kb2(jb,ja))
                                    endif
!
                                    if (lgroupvelocity.and.ja.ne.jb) then
!                                   
!  Group velocities                     
!                                   
                                      d2k  = d2trm*cosvec(kb2(jb,ja))
                                      d2ks = d2trm*sinvec(kb2(jb,ja))
                                      cdk(1) = dcmplx(d2k*xvec(kb2(jb,ja)),d2ks*xvec(kb2(jb,ja)))*dcmplx(0.0_dp,1.0_dp)
                                      cdk(2) = dcmplx(d2k*yvec(kb2(jb,ja)),d2ks*yvec(kb2(jb,ja)))*dcmplx(0.0_dp,1.0_dp)
                                      cdk(3) = dcmplx(d2k*zvec(kb2(jb,ja)),d2ks*zvec(kb2(jb,ja)))*dcmplx(0.0_dp,1.0_dp)
!                                       
                                      if (ja.gt.jb) then
                                        derv2dk(1:3,isxyz(ibet,jb),isxyzl(ialp)) = derv2dk(1:3,isxyz(ibet,jb),isxyzl(ialp)) +  &
                                            cdk(1:3)
                                      elseif (ja.lt.jb) then
                                        derv2dk(1:3,isxyz(ibet,jb),isxyzl(ialp)) = derv2dk(1:3,isxyz(ibet,jb),isxyzl(ialp)) +  &
                                            conjg(cdk(1:3))
                                      endif
                                    endif
!
                                    if (ja.gt.ib) then
                                      derv2(isxyz(ibet,ib),isxyzl(ialp)) = derv2(isxyz(ibet,ib),isxyzl(ialp)) -  &
                                          d2trm*cosvec(kb2(ib,ja))
                                      dervi(isxyz(ibet,ib),isxyzl(ialp)) = dervi(isxyz(ibet,ib),isxyzl(ialp)) -  &
                                          d2trm*sinvec(kb2(ib,ja))
                                    elseif (ja.lt.ib) then
                                      derv2(isxyz(ibet,ib),isxyzl(ialp)) = derv2(isxyz(ibet,ib),isxyzl(ialp)) -  &
                                          d2trm*cosvec(kb2(ib,ja))
                                      dervi(isxyz(ibet,ib),isxyzl(ialp)) = dervi(isxyz(ibet,ib),isxyzl(ialp)) +  &
                                          d2trm*sinvec(kb2(ib,ja))
                                    endif
!
                                    if (lgroupvelocity.and.ja.ne.ib) then
!                                   
!  Group velocities                     
!                                   
                                      d2k  = d2trm*cosvec(kb2(ib,ja))
                                      d2ks = d2trm*sinvec(kb2(ib,ja))
                                      cdk(1) = dcmplx(d2k*xvec(kb2(ib,ja)),d2ks*xvec(kb2(ib,ja)))*dcmplx(0.0_dp,1.0_dp)
                                      cdk(2) = dcmplx(d2k*yvec(kb2(ib,ja)),d2ks*yvec(kb2(ib,ja)))*dcmplx(0.0_dp,1.0_dp)
                                      cdk(3) = dcmplx(d2k*zvec(kb2(ib,ja)),d2ks*zvec(kb2(ib,ja)))*dcmplx(0.0_dp,1.0_dp)
!                                       
                                      if (ja.gt.ib) then
                                        derv2dk(1:3,isxyz(ibet,ib),isxyzl(ialp)) = derv2dk(1:3,isxyz(ibet,ib),isxyzl(ialp)) -  &
                                            cdk(1:3)
                                      elseif (ja.lt.ib) then
                                        derv2dk(1:3,isxyz(ibet,ib),isxyzl(ialp)) = derv2dk(1:3,isxyz(ibet,ib),isxyzl(ialp)) -  &
                                            conjg(cdk(1:3))
                                      endif
                                    endif
                                  endif
!
                                  if (ivec.ne.jvec) then
                                    if (ib.eq.1) then
                                      if (ia.gt.ib) then
                                        derv2(isxyz(ialp,ia),isxyzl(ibet)) = derv2(isxyz(ialp,ia),isxyzl(ibet)) +  &
                                          d2trm*cosvec(kb2(ia,ib))
                                        dervi(isxyz(ialp,ia),isxyzl(ibet)) = dervi(isxyz(ialp,ia),isxyzl(ibet)) -  &
                                          d2trm*sinvec(kb2(ia,ib))
                                      elseif (ia.lt.ib) then
                                        derv2(isxyz(ialp,ia),isxyzl(ibet)) = derv2(isxyz(ialp,ia),isxyzl(ibet)) +  &
                                          d2trm*cosvec(kb2(ia,ib))
                                        dervi(isxyz(ialp,ia),isxyzl(ibet)) = dervi(isxyz(ialp,ia),isxyzl(ibet)) +  &
                                          d2trm*sinvec(kb2(ia,ib))
                                      endif
!
                                      if (lgroupvelocity.and.ia.ne.ib) then
!                                   
!  Group velocities                     
!                                   
                                        d2k  = d2trm*cosvec(kb2(ia,ib))
                                        d2ks = d2trm*sinvec(kb2(ia,ib))
                                        cdk(1) = dcmplx(d2k*xvec(kb2(ia,ib)),d2ks*xvec(kb2(ia,ib)))*dcmplx(0.0_dp,1.0_dp)
                                        cdk(2) = dcmplx(d2k*yvec(kb2(ia,ib)),d2ks*yvec(kb2(ia,ib)))*dcmplx(0.0_dp,1.0_dp)
                                        cdk(3) = dcmplx(d2k*zvec(kb2(ia,ib)),d2ks*zvec(kb2(ia,ib)))*dcmplx(0.0_dp,1.0_dp)
!                                       
                                        if (ia.gt.ib) then
                                          derv2dk(1:3,isxyz(ialp,ia),isxyzl(ibet)) = derv2dk(1:3,isxyz(ialp,ia),isxyzl(ibet)) +  &
                                            conjg(cdk(1:3))
                                        elseif (ia.lt.ib) then
                                          derv2dk(1:3,isxyz(ialp,ia),isxyzl(ibet)) = derv2dk(1:3,isxyz(ialp,ia),isxyzl(ibet)) +  &
                                            cdk(1:3)
                                        endif
                                      endif
                                    endif
!
                                    if (jb.eq.1) then
                                      if (ia.gt.jb) then
                                        derv2(isxyz(ialp,ia),isxyzl(ibet)) = derv2(isxyz(ialp,ia),isxyzl(ibet)) -  &
                                          d2trm*cosvec(kb2(ia,jb))
                                        dervi(isxyz(ialp,ia),isxyzl(ibet)) = dervi(isxyz(ialp,ia),isxyzl(ibet)) +  &
                                          d2trm*sinvec(kb2(ia,jb))
                                      elseif (ia.lt.jb) then
                                        derv2(isxyz(ialp,ia),isxyzl(ibet)) = derv2(isxyz(ialp,ia),isxyzl(ibet)) -  &
                                          d2trm*cosvec(kb2(ia,jb))
                                        dervi(isxyz(ialp,ia),isxyzl(ibet)) = dervi(isxyz(ialp,ia),isxyzl(ibet)) -  &
                                          d2trm*sinvec(kb2(ia,jb))
                                      endif
!
                                      if (lgroupvelocity.and.ia.ne.jb) then
!
!  Group velocities
!
                                        d2k  = d2trm*cosvec(kb2(ia,jb))
                                        d2ks = d2trm*sinvec(kb2(ia,jb))
                                        cdk(1) = dcmplx(d2k*xvec(kb2(ia,jb)),d2ks*xvec(kb2(ia,jb)))*dcmplx(0.0_dp,1.0_dp)
                                        cdk(2) = dcmplx(d2k*yvec(kb2(ia,jb)),d2ks*yvec(kb2(ia,jb)))*dcmplx(0.0_dp,1.0_dp)
                                        cdk(3) = dcmplx(d2k*zvec(kb2(ia,jb)),d2ks*zvec(kb2(ia,jb)))*dcmplx(0.0_dp,1.0_dp)
!
                                        if (ia.gt.jb) then
                                          derv2dk(1:3,isxyz(ialp,ia),isxyzl(ibet)) = derv2dk(1:3,isxyz(ialp,ia),isxyzl(ibet)) -  &
                                            conjg(cdk(1:3))
                                        elseif (ia.lt.jb) then
                                          derv2dk(1:3,isxyz(ialp,ia),isxyzl(ibet)) = derv2dk(1:3,isxyz(ialp,ia),isxyzl(ibet)) -  &
                                            cdk(1:3)
                                        endif
                                      endif
!
                                      if (ja.gt.jb) then
                                        derv2(isxyz(ialp,ja),isxyzl(ibet)) = derv2(isxyz(ialp,ja),isxyzl(ibet)) +  &
                                          d2trm*cosvec(kb2(ja,jb))
                                        dervi(isxyz(ialp,ja),isxyzl(ibet)) = dervi(isxyz(ialp,ja),isxyzl(ibet)) -  &
                                          d2trm*sinvec(kb2(ja,jb))
                                      elseif (ja.lt.jb) then
                                        derv2(isxyz(ialp,ja),isxyzl(ibet)) = derv2(isxyz(ialp,ja),isxyzl(ibet)) +  &
                                          d2trm*cosvec(kb2(ja,jb))
                                        dervi(isxyz(ialp,ja),isxyzl(ibet)) = dervi(isxyz(ialp,ja),isxyzl(ibet)) +  &
                                          d2trm*sinvec(kb2(ja,jb))
                                      endif
!
                                      if (lgroupvelocity.and.ja.ne.jb) then
!
!  Group velocities
!
                                        d2k  = d2trm*cosvec(kb2(ja,jb))
                                        d2ks = d2trm*sinvec(kb2(ja,jb))
                                        cdk(1) = dcmplx(d2k*xvec(kb2(ja,jb)),d2ks*xvec(kb2(ja,jb)))*dcmplx(0.0_dp,1.0_dp)
                                        cdk(2) = dcmplx(d2k*yvec(kb2(ja,jb)),d2ks*yvec(kb2(ja,jb)))*dcmplx(0.0_dp,1.0_dp)
                                        cdk(3) = dcmplx(d2k*zvec(kb2(ja,jb)),d2ks*zvec(kb2(ja,jb)))*dcmplx(0.0_dp,1.0_dp)
!
                                        if (ja.gt.jb) then
                                          derv2dk(1:3,isxyz(ialp,ja),isxyzl(ibet)) = derv2dk(1:3,isxyz(ialp,ja),isxyzl(ibet)) +  &
                                            conjg(cdk(1:3))
                                        elseif (ja.lt.jb) then
                                          derv2dk(1:3,isxyz(ialp,ja),isxyzl(ibet)) = derv2dk(1:3,isxyz(ialp,ja),isxyzl(ibet)) +  &
                                            cdk(1:3)
                                        endif
                                      endif
                                    endif
!
                                    if (ib.eq.1) then
                                      if (ja.gt.ib) then
                                        derv2(isxyz(ialp,ja),isxyzl(ibet)) = derv2(isxyz(ialp,ja),isxyzl(ibet)) -  &
                                          d2trm*cosvec(kb2(ja,ib))
                                        dervi(isxyz(ialp,ja),isxyzl(ibet)) = dervi(isxyz(ialp,ja),isxyzl(ibet)) +  &
                                          d2trm*sinvec(kb2(ja,ib))
                                      elseif (ja.lt.ib) then
                                        derv2(isxyz(ialp,ja),isxyzl(ibet)) = derv2(isxyz(ialp,ja),isxyzl(ibet)) -  &
                                          d2trm*cosvec(kb2(ja,ib))
                                        dervi(isxyz(ialp,ja),isxyzl(ibet)) = dervi(isxyz(ialp,ja),isxyzl(ibet)) -  &
                                          d2trm*sinvec(kb2(ja,ib))
                                      endif
!
                                      if (lgroupvelocity.and.ja.ne.ib) then
!
!  Group velocities
!
                                        d2k  = d2trm*cosvec(kb2(ja,ib))
                                        d2ks = d2trm*sinvec(kb2(ja,ib))
                                        cdk(1) = dcmplx(d2k*xvec(kb2(ja,ib)),d2ks*xvec(kb2(ja,ib)))*dcmplx(0.0_dp,1.0_dp)
                                        cdk(2) = dcmplx(d2k*yvec(kb2(ja,ib)),d2ks*yvec(kb2(ja,ib)))*dcmplx(0.0_dp,1.0_dp)
                                        cdk(3) = dcmplx(d2k*zvec(kb2(ja,ib)),d2ks*zvec(kb2(ja,ib)))*dcmplx(0.0_dp,1.0_dp)
!
                                        if (ja.gt.ib) then
                                          derv2dk(1:3,isxyz(ialp,ja),isxyzl(ibet)) = derv2dk(1:3,isxyz(ialp,ja),isxyzl(ibet)) -  &
                                            conjg(cdk(1:3))
                                        elseif (ja.lt.ib) then
                                          derv2dk(1:3,isxyz(ialp,ja),isxyzl(ibet)) = derv2dk(1:3,isxyz(ialp,ja),isxyzl(ibet)) -  &
                                            cdk(1:3)
                                        endif
                                      endif
                                    endif
                                  endif
                                enddo
                              enddo
                            enddo
                          enddo
!
!  End of inner loops over atoms and cell vectors
!
                        enddo mnloop
                      enddo nloop
                    enddo mmloop
                  enddo mloop
                enddo llloop
              enddo lloop
            enddo jjloop
          enddo kloop
        enddo iiloop
      enddo jloop
    enddo iloop
!***********************************************
!  Loop over middle site 1 / i - end atom case *
!***********************************************
    iloop2: do i = 1,numat
      ni = nat(i)
      ntypi = nftype(i)
      oci = occuf(i)
!
!  Check i is allowed for n
!
      if (lmatch(ni,ntypi,nsspec1(np),nsptyp1(np),.true.)) then
        nt2 = nsspec2(np)
        ntyp2 = nsptyp2(np)
        nt3 = nsspec3(np)
        nt4 = nsspec4(np)
        nt5 = nsspec5(np)
        nt6 = nsspec6(np)
        ntyp3 = nsptyp3(np)
        ntyp4 = nsptyp4(np)
        ntyp5 = nsptyp5(np)
        ntyp6 = nsptyp6(np)
        tr2 = ttr2
        tr3 = ttr3
        tr4 = ttr4
        tr5 = ttr5
      elseif (lmatch(ni,ntypi,nsspec2(np),nsptyp2(np),.true.)) then
        nt2 = nsspec1(np)
        ntyp2 = nsptyp1(np)
        nt3 = nsspec5(np)
        nt4 = nsspec6(np)
        nt5 = nsspec3(np)
        nt6 = nsspec4(np)
        ntyp3 = nsptyp5(np)
        ntyp4 = nsptyp6(np)
        ntyp5 = nsptyp3(np)
        ntyp6 = nsptyp4(np)
        tr2 = ttr4
        tr3 = ttr5
        tr4 = ttr2
        tr5 = ttr3
      else
        cycle iloop2
      endif
!
!  Set flags as to whether potentials are symmetric at the two ends
!
      ltsym34 = (lmatch(nt3,ntyp3,nt4,ntyp4,.true.).or.lmatch(nt4,ntyp4,nt3,ntyp3,.true.))
      ltsym56 = (lmatch(nt5,ntyp5,nt6,ntyp6,.true.).or.lmatch(nt6,ntyp6,nt5,ntyp5,.true.))
!
!  i has been accepted
!
      isatom(1) = i
      isxyz(1,1) = 3*(i - 1) + 1
      isxyz(2,1) = 3*(i - 1) + 2
      isxyz(3,1) = 3*(i - 1) + 3
!
      sxyz(1,1) = xclat(i)
      sxyz(2,1) = yclat(i)
      sxyz(3,1) = zclat(i)
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
!
!  Loop over all j in this algorithm
!
      jmin = 1
      jmax = numat
!********************************
!  Loop over middle site 2 / j  *
!********************************
      jloop2: do j = jmin,jmax
        nj = nat(j)
        ntypj = nftype(j)
        ocj = occuf(j)
!
!  Check j is allowed for n
!
        if (.not.lmatch(nj,ntypj,nt2,ntyp2,.true.)) cycle jloop2
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
        if (lintra_only.and..not.lmolok) cycle jloop2
        if (lbtyp.and..not.lmolok) cycle jloop2
        xc2t = xclat(j)
        yc2t = yclat(j)
        zc2t = zclat(j)
        x21t = xc2t - sxyz(1,1)
        y21t = yc2t - sxyz(2,1)
        z21t = zc2t - sxyz(3,1)
!
        lmolokij = lmolok
!
!  Check r21 is OK
!  Loop over cell vectors
!
        iiloop2: do ii = 1,iimax
          x21 = x21t + xvec1cell(ii)
          y21 = y21t + yvec1cell(ii)
          z21 = z21t + zvec1cell(ii)
          r212 = x21*x21 + y21*y21 + z21*z21
          if (r212.lt.1d-10) cycle iiloop2
!
!  Reset lmolok to preloop2 value
!
          lmolok = lmolokij
!
!  Molecule checking
!
          lbonded = .false.
          if (lmolok) then
            if (ndim.eq.0) then
              if (linter_only) cycle iiloop2
              if (lbtyp) then
                call bonded(lbonded,l2bonds,nbtypeij,nbtypeij2,i,j,0_i4,0_i4,0_i4)
                if (.not.lbonded) cycle iiloop2
!
!  Check central bond type for correct order
!
                if (n6botype(1,np).gt.0) then
                  if (n6botype(1,np).ne.nbtypeij) cycle iiloop2
                endif
                if (n6botype(2,np).gt.0) then
                  if (n6botype(2,np).ne.nbtypeij2) cycle iiloop2
                endif
              endif
            else
              call lintoijk(ixx,iyy,izz,ii,imaxl,jmaxl,kmaxl)
              if (lbtyp) then
                call bonded(lbonded,l2bonds,nbtypeij,nbtypeij2,i,j,ixx,iyy,izz)
                if (.not.lbonded) cycle iiloop2
!
!  Check central bond type for correct order
!
                if (n6botype(1,np).gt.0) then
                  if (n6botype(1,np).ne.nbtypeij) cycle iiloop2
                endif
                if (n6botype(2,np).gt.0) then
                  if (n6botype(2,np).ne.nbtypeij2) cycle iiloop2
                endif
                lsamemol = (lbonded.or.l2bonds)
              else
                lsamemol = .false.
              endif
              if (.not.lsamemol) then
                call samemol(lsamemol,nmi,ixx,iyy,izz,ixj,iyj,izj)
              endif
              if (lintra_only.and..not.lsamemol) cycle iiloop2
              if (linter_only.and.lsamemol) cycle iiloop2
            endif
          endif
!
!  Distance checking
!
          if (r212.gt.tr1.and.(.not.lbtyp.or..not.lbonded)) cycle iiloop2
!
!  j has been accepted
!
          isatom(2) = j
          isxyz(1,2) = 3*(j - 1) + 1
          isxyz(2,2) = 3*(j - 1) + 2
          isxyz(3,2) = 3*(j - 1) + 3
          sxyz(1,2) = x21 + sxyz(1,1)
          sxyz(2,2) = y21 + sxyz(2,1)
          sxyz(3,2) = z21 + sxyz(3,1)
!*****************************************
!  Loop over first end site bonded to i  *
!*****************************************
          kloop2: do kloc = 1,natomsonnode
            k = node2atom(kloc)
            nk = nat(k)
            ntypk = nftype(k)
            ock = occuf(k)
!
!  Check k is allowed for n
!
            if (.not.lmatch(nk,ntypk,nt3,ntyp3,.true.)) cycle kloop2
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
            if (lintra_only.and..not.lmolok) cycle kloop2
            if (lbtyp.and..not.lmolok) cycle kloop2
            xc3t = xclat(k)
            yc3t = yclat(k)
            zc3t = zclat(k)
            x31t = xc3t - sxyz(1,1)
            y31t = yc3t - sxyz(2,1)
            z31t = zc3t - sxyz(3,1)
!
            lmolokik = lmolok
!
!  Check r31 is OK
!  Loop over cell vectors
!
            jjloop2: do jj = 1,iimax
              x31 = x31t + xvec1cell(jj)
              y31 = y31t + yvec1cell(jj)
              z31 = z31t + zvec1cell(jj)
              r312 = x31*x31 + y31*y31 + z31*z31
              if (r312.lt.1d-10) cycle jjloop2
!
!  Prevent atoms i and k being the same atom
!
              if (k.eq.i.and.jj.eq.iimid) cycle jjloop2
!
!  Prevent atoms j and k being the same atom
!
              if (k.eq.j.and.jj.eq.ii) cycle jjloop2
!
!  Reset lmolok to preloop2 value
!
              lmolok = lmolokik
!
!  Molecule checking
!
              lbonded = .false.
              if (lmolok) then
                if (ndim.eq.0) then
                  if (linter_only) cycle jjloop2
                  if (lbtyp) then
                    call bonded(lbonded,l2bonds,nbtypeik,nbtypeik2,i,k,0_i4,0_i4,0_i4)
                    if (.not.lbonded) cycle jjloop2
                  endif
                else
                  call lintoijk(jxx,jyy,jzz,jj,imaxl,jmaxl,kmaxl)
                  if (lbtyp) then
                    call bonded(lbonded,l2bonds,nbtypeik,nbtypeik2,i,k,jxx,jyy,jzz)
                    if (.not.lbonded) cycle jjloop2
                    lsamemol = (lbonded.or.l2bonds)
                  else
                    lsamemol = .false.
                  endif
                  if (.not.lsamemol) then
                    call samemol(lsamemol,nmi,jxx,jyy,jzz,ixk,iyk,izk)
                  endif
                  if (lintra_only.and..not.lsamemol) cycle jjloop2
                  if (linter_only.and.lsamemol) cycle jjloop2
                endif
              endif
!
!  Distance checking
!
              if (r312.gt.tr2.and.(.not.lbtyp.or..not.lbonded)) cycle jjloop2
!
!  k has been accepted
!
              isatom(3) = k
              isxyz(1,3) = 3*(k - 1) + 1
              isxyz(2,3) = 3*(k - 1) + 2
              isxyz(3,3) = 3*(k - 1) + 3
!
              isxyzl(1) = 3*(kloc - 1) + 1
              isxyzl(2) = 3*(kloc - 1) + 2
              isxyzl(3) = 3*(kloc - 1) + 3
!
              sxyz(1,3) = x31 + sxyz(1,1)
              sxyz(2,3) = y31 + sxyz(2,1)
              sxyz(3,3) = z31 + sxyz(3,1)
!************************************
!  Loop over second end site for i  *
!************************************
!
!  Set l looping indices
!
              lmin = 1
              lmax = numat
!
              lloop2: do l = lmin,lmax
                nl = nat(l)
                ntypl = nftype(l)
                ocl = occuf(l)
!
!  Check l is allowed for n
!
                if (.not.lmatch(nl,ntypl,nt4,ntyp4,.true.)) cycle lloop2
!
!  Molecule handling
!
                if (lmol.and.lneedmol) then
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
                if (lintra_only.and..not.lmolok) cycle lloop2
                if (lbtyp.and..not.lmolok) cycle lloop2
                xc4t = xclat(l)
                yc4t = yclat(l)
                zc4t = zclat(l)
                x41t = xc4t - sxyz(1,1)
                y41t = yc4t - sxyz(2,1)
                z41t = zc4t - sxyz(3,1)
                if (ltsym34.and.k.eq.l) then
                  liimax = jj - 1
                else
                  liimax = iimax
                endif
!
                lmolokil = lmolok
!
!  Check r41 is OK
!  Loop over cell vectors
!
                llloop2: do ll = 1,liimax
                  x41 = x41t + xvec1cell(ll)
                  y41 = y41t + yvec1cell(ll)
                  z41 = z41t + zvec1cell(ll)
                  r412 = x41*x41 + y41*y41 + z41*z41
                  if (r412.lt.1d-10) cycle llloop2
!
!  Prevent atoms i and l being the same atom
!
                  if (l.eq.i.and.ll.eq.iimid) cycle llloop2
!
!  Prevent atoms j and l being the same atom
!
                  if (l.eq.j.and.ll.eq.ii) cycle llloop2
!
!  Prevent atoms k and l being the same atom
!
                  if (l.eq.k.and.ll.eq.jj) cycle llloop2
!
!  Reset lmolok to preloop value
!
                  lmolok = lmolokil
!
!  Molecule checking
!
                  lbonded = .false.
                  if (lmolok) then
                    if (ndim.eq.0) then
                      if (linter_only) cycle llloop2
                      if (lbtyp) then
                        call bonded(lbonded,l2bonds,nbtypeil,nbtypeil2,i,l,0_i4,0_i4,0_i4)
                        if (.not.lbonded) cycle llloop2
                      endif
                    else
                      call lintoijk(kxx,kyy,kzz,ll,imaxl,jmaxl,kmaxl)
                      if (lbtyp) then
                        call bonded(lbonded,l2bonds,nbtypeil,nbtypeil2,i,l,kxx,kyy,kzz)
                        if (.not.lbonded) cycle llloop2
                        lsamemol = (lbonded.or.l2bonds)
                      else
                        lsamemol = .false.
                      endif
                      if (.not.lsamemol) then
                        call samemol(lsamemol,nmi,kxx,kyy,kzz,ixl,iyl,izl)
                      endif
                      if (lintra_only.and..not.lsamemol) cycle llloop2
                      if (linter_only.and.lsamemol) cycle llloop2
                    endif
                  endif
!
!  Distance checking
!
                  if (r412.gt.tr3.and.(.not.lbtyp.or..not.lbonded)) cycle llloop2
!
!  l has been accepted
!
                  isatom(4) = l
                  isxyz(1,4) = 3*(l - 1) + 1
                  isxyz(2,4) = 3*(l - 1) + 2
                  isxyz(3,4) = 3*(l - 1) + 3
                  sxyz(1,4) = x41 + sxyz(1,1)
                  sxyz(2,4) = y41 + sxyz(2,1)
                  sxyz(3,4) = z41 + sxyz(3,1)
!*****************************************
!  Loop over first end site bonded to j  *
!*****************************************
                  mloop2: do m = 1,numat
                    nm = nat(m)
                    ntypm = nftype(m)
                    ocm = occuf(m)
!
!  Check m is allowed for n
!
                    if (.not.lmatch(nm,ntypm,nt5,ntyp5,.true.)) cycle mloop2
                    if (lmol.and.lneedmol) then
!
!  Molecule handling
!
                      nmm = natmol(m)
                      if (ndim.gt.0) then
                        indmm = nmolind(m)
                        call mindtoijk(indmm,ixm,iym,izm)
                        ixm = ixm - ixj
                        iym = iym - iyj
                        izm = izm - izj
                      endif
                      lmolok = (nmi.eq.nmm.and.nmi.gt.0)
                    else
                      lmolok = .false.
                    endif
!
!  Check for intra and but not in same molecule
!
                    if (lintra_only.and..not.lmolok) cycle mloop2
                    if (lbtyp.and..not.lmolok) cycle mloop2
                    xc5t = xclat(m)
                    yc5t = yclat(m)
                    zc5t = zclat(m)
                    x52t = xc5t - sxyz(1,2)
                    y52t = yc5t - sxyz(2,2)
                    z52t = zc5t - sxyz(3,2)
!
                    lmolokjm = lmolok
!
!  Check r52 is OK
!  Loop over cell vectors
!
                    mmloop2: do mm = 1,iimax
                      x52 = x52t + xvec1cell(mm)
                      y52 = y52t + yvec1cell(mm)
                      z52 = z52t + zvec1cell(mm)
                      r522 = x52*x52 + y52*y52 + z52*z52
                      if (r522.lt.1d-10) cycle mmloop2
!
!  Prevent atoms i and m being the same atom
!
                      if (m.eq.i.and.mm.eq.iimid) cycle mmloop2
!
!  Prevent atoms j and m being the same atom
!
                      if (m.eq.j.and.mm.eq.ii) cycle mmloop2
!
!  Reset lmolok to preloop value
!
                      lmolok = lmolokjm
!
!  Molecule checking
!
                      lbonded = .false.
                      if (lmolok) then
                        if (ndim.eq.0) then
                          if (linter_only) cycle mmloop2
                          if (lbtyp) then
                            call bonded(lbonded,l2bonds,nbtypejm,nbtypejm2,j,m,0_i4,0_i4,0_i4)
                            if (.not.lbonded) cycle mmloop2
                          endif
                        else
                          call lintoijk(mxx,myy,mzz,mm,imaxl,jmaxl,kmaxl)
                          if (lbtyp) then
                            call bonded(lbonded,l2bonds,nbtypejm,nbtypejm2,j,m,mxx-ixx,myy-iyy,mzz-izz)
                            if (.not.lbonded) cycle mmloop2
                            lsamemol = (lbonded.or.l2bonds)
                          else
                            lsamemol = .false.
                          endif
                          if (.not.lsamemol) then
                            call samemol(lsamemol,nmi,jxx,jyy,jzz,ixm,iym,izm)
                          endif
                          if (lintra_only.and..not.lsamemol) cycle mmloop2
                          if (linter_only.and.lsamemol) cycle mmloop2
                        endif
                      endif
!
!  Distance checking
!
                      if (r522.gt.tr4.and.(.not.lbtyp.or..not.lbonded)) cycle mmloop2
!
!  m has been accepted
!
                      isatom(5) = m
                      isxyz(1,5) = 3*(m - 1) + 1
                      isxyz(2,5) = 3*(m - 1) + 2
                      isxyz(3,5) = 3*(m - 1) + 3
                      sxyz(1,5) = x52 + sxyz(1,2)
                      sxyz(2,5) = y52 + sxyz(2,2)
                      sxyz(3,5) = z52 + sxyz(3,2)
!************************************
!  Loop over second end site for j  *
!************************************
!
!  Set n looping indices
!
                      if (ltsym56) then
                        nmin = 1
                        nmax = m
                      else
                        nmin = 1
                        nmax = numat
                      endif
                      nloop2: do n = nmin,nmax
                        nn = nat(n)
                        ntypn = nftype(n)
                        ocn = occuf(n)
!
!  Check n is allowed for n
!
                        if (.not.lmatch(nn,ntypn,nt6,ntyp6,.true.)) cycle nloop2
!
!  Molecule handling
!
                        if (lmol.and.lneedmol) then
                          nmn = natmol(n)
                          if (ndim.gt.0) then
                            indmn = nmolind(n)
                            call mindtoijk(indmn,ixn,iyn,izn)
                            ixn = ixn - ixj
                            iyn = iyn - iyj
                            izn = izn - izj
                          endif
                          lmolok = (nmi.eq.nmn.and.nmi.gt.0)
                        else
                          lmolok = .false.
                        endif
!
!  Check for intra and but not in same molecule
!
                        if (lintra_only.and..not.lmolok) cycle nloop2
                        if (lbtyp.and..not.lmolok) cycle nloop2
                        xc6t = xclat(n)
                        yc6t = yclat(n)
                        zc6t = zclat(n)
                        x62t = xc6t - sxyz(1,2)
                        y62t = yc6t - sxyz(2,2)
                        z62t = zc6t - sxyz(3,2)
                        if (ltsym56.and.m.eq.n) then
                          niimax = mm - 1
                        else
                          niimax = iimax
                        endif
!
                        lmolokjn = lmolok
!
!  Check r62 is OK
!  Loop over cell vectors
!
                        mnloop2: do mn = 1,niimax
                          x62 = x62t + xvec1cell(mn)
                          y62 = y62t + yvec1cell(mn)
                          z62 = z62t + zvec1cell(mn)
                          r622 = x62*x62 + y62*y62 + z62*z62
                          if (r622.lt.1d-10) cycle mnloop2
!
!  Prevent atoms i and n being the same atom
!
                          if (n.eq.i.and.mn.eq.iimid) cycle mnloop2
!
!  Prevent atoms j and n being the same atom
!
                          if (n.eq.j.and.mn.eq.ii) cycle mnloop2
!
!  Prevent atoms m and n being the same atom
!
                          if (n.eq.m.and.mn.eq.mm) cycle mnloop2
!
!  Reset lmolok to preloop value
!
                          lmolok = lmolokjn
!
!  Molecule checking
!
                          lbonded = .false.
                          if (lmolok) then
                            if (ndim.eq.0) then
                              if (linter_only) cycle mnloop2
                              if (lbtyp) then
                                call bonded(lbonded,l2bonds,nbtypejn,nbtypejn2,j,n,0_i4,0_i4,0_i4)
                                if (.not.lbonded) cycle mnloop2
                              endif
                            else
                              call lintoijk(nxx,nyy,nzz,mn,imaxl,jmaxl,kmaxl)
                              if (lbtyp) then
                                call bonded(lbonded,l2bonds,nbtypejn,nbtypejn2,j,n,nxx-ixx,nyy-iyy,nzz-izz)
                                if (.not.lbonded) cycle mnloop2
                                lsamemol = (lbonded.or.l2bonds)
                              else
                                lsamemol = .false.
                              endif
                              if (.not.lsamemol) then
                                call samemol(lsamemol,nmi,jxx,jyy,jzz,ixn,iyn,izn)
                              endif
                              if (lintra_only.and..not.lsamemol) cycle mnloop2
                              if (linter_only.and.lsamemol) cycle mnloop2
                            endif
                          endif
!
!  Distance checking
!
                          if (r622.gt.tr5.and.(.not.lbtyp.or..not.lbonded)) cycle mnloop2
!
!  n has been accepted
!
                          isatom(6) = n
                          isxyz(1,6) = 3*(n - 1) + 1
                          isxyz(2,6) = 3*(n - 1) + 2
                          isxyz(3,6) = 3*(n - 1) + 3
                          sxyz(1,6) = x62 + sxyz(1,2)
                          sxyz(2,6) = y62 + sxyz(2,2)
                          sxyz(3,6) = z62 + sxyz(3,2)
!********************************
!  Valid six-body term located  *
!********************************
!
!  Calculate vectors and remaining distances
!
                          do ivec = 1,6
                            do jvec = 1,6
                              svec(1,jvec,ivec) = sxyz(1,ivec) - sxyz(1,jvec)
                              svec(2,jvec,ivec) = sxyz(2,ivec) - sxyz(2,jvec)
                              svec(3,jvec,ivec) = sxyz(3,ivec) - sxyz(3,jvec)
                            enddo
                          enddo
                          indvec = 0
                          do ivec = 1,5
                            do jvec = ivec+1,6
                              indvec = indvec + 1
                              sdist(indvec) = svec(1,jvec,ivec)**2 + svec(2,jvec,ivec)**2 + svec(3,jvec,ivec)**2
                              sdist(indvec) = sqrt(sdist(indvec))
                            enddo
                          enddo
!
                          ofct = oci*ocj*ock*ocl*ocm*ocn
                          rko = rksix*ofct
                          call sixbody(np,nsixtype,sdist,eterm,e1d,e2d,e3d,rko,.true.,.true.,.false.)
!*************************
!  Internal derivatives  *
!*************************
!
!  Compute phase factors
!
                          indvec = 0
                          do ivec = 1,5
                            do jvec = ivec+1,6
                              indvec = indvec + 1
                              phase = svec(1,jvec,ivec)*xkv + svec(2,jvec,ivec)*ykv + svec(3,jvec,ivec)*zkv
                              sinvec(indvec) = sin(phase)
                              if (isatom(ivec).eq.isatom(jvec)) then
                                cosvec(indvec) = cos(phase) - 1.0_dp
                              else
                                cosvec(indvec) = cos(phase) 
                              endif
                              xvec(indvec) = svec(1,jvec,ivec)
                              yvec(indvec) = svec(2,jvec,ivec)
                              zvec(indvec) = svec(3,jvec,ivec)
                            enddo
                          enddo
!
                          indvec = 0
                          do ivec = 1,15
                            ja = iveca(1,ivec)
                            ia = iveca(2,ivec)
!
                            if (ia.eq.3) then
                              derv2(isxyz(1,ja),isxyzl(1)) = derv2(isxyz(1,ja),isxyzl(1)) - e1d(ivec)*cosvec(ivec)
                              derv2(isxyz(2,ja),isxyzl(2)) = derv2(isxyz(2,ja),isxyzl(2)) - e1d(ivec)*cosvec(ivec)
                              derv2(isxyz(3,ja),isxyzl(3)) = derv2(isxyz(3,ja),isxyzl(3)) - e1d(ivec)*cosvec(ivec)
                              if (ia.gt.ja) then
                                dervi(isxyz(1,ja),isxyzl(1)) = dervi(isxyz(1,ja),isxyzl(1)) - e1d(ivec)*sinvec(ivec)
                                dervi(isxyz(2,ja),isxyzl(2)) = dervi(isxyz(2,ja),isxyzl(2)) - e1d(ivec)*sinvec(ivec)
                                dervi(isxyz(3,ja),isxyzl(3)) = dervi(isxyz(3,ja),isxyzl(3)) - e1d(ivec)*sinvec(ivec)
                              else
                                dervi(isxyz(1,ja),isxyzl(1)) = dervi(isxyz(1,ja),isxyzl(1)) + e1d(ivec)*sinvec(ivec)
                                dervi(isxyz(2,ja),isxyzl(2)) = dervi(isxyz(2,ja),isxyzl(2)) + e1d(ivec)*sinvec(ivec)
                                dervi(isxyz(3,ja),isxyzl(3)) = dervi(isxyz(3,ja),isxyzl(3)) + e1d(ivec)*sinvec(ivec)
                              endif
                              if (lgroupvelocity) then
!
!  Group velocities
!
                                d2k  = e1d(ivec)*cosvec(ivec)
                                d2ks = e1d(ivec)*sinvec(ivec)
                                cdk(1) = dcmplx(d2k*xvec(ivec),d2ks*xvec(ivec))*dcmplx(0.0_dp,1.0_dp)
                                cdk(2) = dcmplx(d2k*yvec(ivec),d2ks*yvec(ivec))*dcmplx(0.0_dp,1.0_dp)
                                cdk(3) = dcmplx(d2k*zvec(ivec),d2ks*zvec(ivec))*dcmplx(0.0_dp,1.0_dp)
!
                                if (ia.gt.ja) then
                                  derv2dk(1:3,isxyz(1,ja),isxyzl(1)) = derv2dk(1:3,isxyz(1,ja),isxyzl(1)) - cdk(1:3)
                                  derv2dk(1:3,isxyz(2,ja),isxyzl(2)) = derv2dk(1:3,isxyz(2,ja),isxyzl(2)) - cdk(1:3)
                                  derv2dk(1:3,isxyz(3,ja),isxyzl(3)) = derv2dk(1:3,isxyz(3,ja),isxyzl(3)) - cdk(1:3)
                                else
                                  derv2dk(1:3,isxyz(1,ja),isxyzl(1)) = derv2dk(1:3,isxyz(1,ja),isxyzl(1)) - conjg(cdk(1:3))
                                  derv2dk(1:3,isxyz(2,ja),isxyzl(2)) = derv2dk(1:3,isxyz(2,ja),isxyzl(2)) - conjg(cdk(1:3))
                                  derv2dk(1:3,isxyz(3,ja),isxyzl(3)) = derv2dk(1:3,isxyz(3,ja),isxyzl(3)) - conjg(cdk(1:3))
                                endif
                              endif
                            endif
!
                            if (ja.eq.3) then
                              derv2(isxyz(1,ia),isxyzl(1)) = derv2(isxyz(1,ia),isxyzl(1)) - e1d(ivec)*cosvec(ivec)
                              derv2(isxyz(2,ia),isxyzl(2)) = derv2(isxyz(2,ia),isxyzl(2)) - e1d(ivec)*cosvec(ivec)
                              derv2(isxyz(3,ia),isxyzl(3)) = derv2(isxyz(3,ia),isxyzl(3)) - e1d(ivec)*cosvec(ivec)
                              if (ia.gt.ja) then
                                dervi(isxyz(1,ia),isxyzl(1)) = dervi(isxyz(1,ia),isxyzl(1)) + e1d(ivec)*sinvec(ivec)
                                dervi(isxyz(2,ia),isxyzl(2)) = dervi(isxyz(2,ia),isxyzl(2)) + e1d(ivec)*sinvec(ivec)
                                dervi(isxyz(3,ia),isxyzl(3)) = dervi(isxyz(3,ia),isxyzl(3)) + e1d(ivec)*sinvec(ivec)
                              else
                                dervi(isxyz(1,ia),isxyzl(1)) = dervi(isxyz(1,ia),isxyzl(1)) - e1d(ivec)*sinvec(ivec)
                                dervi(isxyz(2,ia),isxyzl(2)) = dervi(isxyz(2,ia),isxyzl(2)) - e1d(ivec)*sinvec(ivec)
                                dervi(isxyz(3,ia),isxyzl(3)) = dervi(isxyz(3,ia),isxyzl(3)) - e1d(ivec)*sinvec(ivec)
                              endif
                              if (lgroupvelocity) then
!
!  Group velocities
!
                                d2k  = e1d(ivec)*cosvec(ivec)
                                d2ks = e1d(ivec)*sinvec(ivec)
                                cdk(1) = dcmplx(d2k*xvec(ivec),d2ks*xvec(ivec))*dcmplx(0.0_dp,1.0_dp)
                                cdk(2) = dcmplx(d2k*yvec(ivec),d2ks*yvec(ivec))*dcmplx(0.0_dp,1.0_dp)
                                cdk(3) = dcmplx(d2k*zvec(ivec),d2ks*zvec(ivec))*dcmplx(0.0_dp,1.0_dp)
!
                                if (ia.gt.ja) then
                                  derv2dk(1:3,isxyz(1,ia),isxyzl(1)) = derv2dk(1:3,isxyz(1,ia),isxyzl(1)) - conjg(cdk(1:3))
                                  derv2dk(1:3,isxyz(2,ia),isxyzl(2)) = derv2dk(1:3,isxyz(2,ia),isxyzl(2)) - conjg(cdk(1:3))
                                  derv2dk(1:3,isxyz(3,ia),isxyzl(3)) = derv2dk(1:3,isxyz(3,ia),isxyzl(3)) - conjg(cdk(1:3))
                                else
                                  derv2dk(1:3,isxyz(1,ia),isxyzl(1)) = derv2dk(1:3,isxyz(1,ia),isxyzl(1)) - cdk(1:3)
                                  derv2dk(1:3,isxyz(2,ia),isxyzl(2)) = derv2dk(1:3,isxyz(2,ia),isxyzl(2)) - cdk(1:3)
                                  derv2dk(1:3,isxyz(3,ia),isxyzl(3)) = derv2dk(1:3,isxyz(3,ia),isxyzl(3)) - cdk(1:3)
                                endif
                              endif
                            endif
!
                            do jvec = ivec,15
                              indvec = indvec + 1 
                              jb = iveca(1,jvec)
                              ib = iveca(2,jvec)
                              do ialp = 1,3
                                do ibet = 1,3
                                  d2trm = svec(ialp,ja,ia)*svec(ibet,jb,ib)*e2d(indvec)
                                  if (ia.eq.3) then
                                    if (ia.gt.ib) then
                                      derv2(isxyz(ibet,ib),isxyzl(ialp)) = derv2(isxyz(ibet,ib),isxyzl(ialp)) +  &
                                          d2trm*cosvec(kb2(ib,ia))
                                      dervi(isxyz(ibet,ib),isxyzl(ialp)) = dervi(isxyz(ibet,ib),isxyzl(ialp)) +  &
                                          d2trm*sinvec(kb2(ib,ia))
                                    elseif (ia.lt.ib) then
                                      derv2(isxyz(ibet,ib),isxyzl(ialp)) = derv2(isxyz(ibet,ib),isxyzl(ialp)) +  &
                                          d2trm*cosvec(kb2(ib,ia))
                                      dervi(isxyz(ibet,ib),isxyzl(ialp)) = dervi(isxyz(ibet,ib),isxyzl(ialp)) -  &
                                          d2trm*sinvec(kb2(ib,ia))
                                    endif
!
                                    if (lgroupvelocity.and.ia.ne.ib) then
!
!  Group velocities
!
                                      d2k  = d2trm*cosvec(kb2(ib,ia))
                                      d2ks = d2trm*sinvec(kb2(ib,ia))
                                      cdk(1) = dcmplx(d2k*xvec(kb2(ib,ia)),d2ks*xvec(kb2(ib,ia)))*dcmplx(0.0_dp,1.0_dp)
                                      cdk(2) = dcmplx(d2k*yvec(kb2(ib,ia)),d2ks*yvec(kb2(ib,ia)))*dcmplx(0.0_dp,1.0_dp)
                                      cdk(3) = dcmplx(d2k*zvec(kb2(ib,ia)),d2ks*zvec(kb2(ib,ia)))*dcmplx(0.0_dp,1.0_dp)
!
                                      if (ia.gt.ib) then
                                        derv2dk(1:3,isxyz(ibet,ib),isxyzl(ialp)) = derv2dk(1:3,isxyz(ibet,ib),isxyzl(ialp)) +  &
                                            cdk(1:3)
                                      elseif (ia.lt.ib) then
                                        derv2dk(1:3,isxyz(ibet,ib),isxyzl(ialp)) = derv2dk(1:3,isxyz(ibet,ib),isxyzl(ialp)) +  &
                                            conjg(cdk(1:3))
                                      endif
                                    endif
!
                                    if (ia.gt.jb) then
                                      derv2(isxyz(ibet,jb),isxyzl(ialp)) = derv2(isxyz(ibet,jb),isxyzl(ialp)) -  &
                                          d2trm*cosvec(kb2(jb,ia))
                                      dervi(isxyz(ibet,jb),isxyzl(ialp)) = dervi(isxyz(ibet,jb),isxyzl(ialp)) -  &
                                          d2trm*sinvec(kb2(jb,ia))
                                    elseif (ia.lt.jb) then
                                      derv2(isxyz(ibet,jb),isxyzl(ialp)) = derv2(isxyz(ibet,jb),isxyzl(ialp)) -  &
                                          d2trm*cosvec(kb2(jb,ia))
                                      dervi(isxyz(ibet,jb),isxyzl(ialp)) = dervi(isxyz(ibet,jb),isxyzl(ialp)) +  &
                                          d2trm*sinvec(kb2(jb,ia))
                                    endif
!
                                    if (lgroupvelocity.and.ia.ne.jb) then
!
!  Group velocities
!
                                      d2k  = d2trm*cosvec(kb2(jb,ia))
                                      d2ks = d2trm*sinvec(kb2(jb,ia))
                                      cdk(1) = dcmplx(d2k*xvec(kb2(jb,ia)),d2ks*xvec(kb2(jb,ia)))*dcmplx(0.0_dp,1.0_dp)
                                      cdk(2) = dcmplx(d2k*yvec(kb2(jb,ia)),d2ks*yvec(kb2(jb,ia)))*dcmplx(0.0_dp,1.0_dp)
                                      cdk(3) = dcmplx(d2k*zvec(kb2(jb,ia)),d2ks*zvec(kb2(jb,ia)))*dcmplx(0.0_dp,1.0_dp)
!
                                      if (ia.gt.jb) then
                                        derv2dk(1:3,isxyz(ibet,jb),isxyzl(ialp)) = derv2dk(1:3,isxyz(ibet,jb),isxyzl(ialp)) -  &
                                            cdk(1:3)
                                      elseif (ia.lt.jb) then
                                        derv2dk(1:3,isxyz(ibet,jb),isxyzl(ialp)) = derv2dk(1:3,isxyz(ibet,jb),isxyzl(ialp)) -  &
                                            conjg(cdk(1:3))
                                      endif
                                    endif   
                                  endif   
!
                                  if (ja.eq.3) then
                                    if (ja.gt.jb) then
                                      derv2(isxyz(ibet,jb),isxyzl(ialp)) = derv2(isxyz(ibet,jb),isxyzl(ialp)) +  &
                                          d2trm*cosvec(kb2(jb,ja))
                                      dervi(isxyz(ibet,jb),isxyzl(ialp)) = dervi(isxyz(ibet,jb),isxyzl(ialp)) +  &
                                          d2trm*sinvec(kb2(jb,ja))
                                    elseif (ja.lt.jb) then
                                      derv2(isxyz(ibet,jb),isxyzl(ialp)) = derv2(isxyz(ibet,jb),isxyzl(ialp)) +  &
                                          d2trm*cosvec(kb2(jb,ja))
                                      dervi(isxyz(ibet,jb),isxyzl(ialp)) = dervi(isxyz(ibet,jb),isxyzl(ialp)) -  &
                                          d2trm*sinvec(kb2(jb,ja))
                                    endif
!
                                    if (lgroupvelocity.and.ja.ne.jb) then
!                                   
!  Group velocities                     
!                                   
                                      d2k  = d2trm*cosvec(kb2(jb,ja))
                                      d2ks = d2trm*sinvec(kb2(jb,ja))
                                      cdk(1) = dcmplx(d2k*xvec(kb2(jb,ja)),d2ks*xvec(kb2(jb,ja)))*dcmplx(0.0_dp,1.0_dp)
                                      cdk(2) = dcmplx(d2k*yvec(kb2(jb,ja)),d2ks*yvec(kb2(jb,ja)))*dcmplx(0.0_dp,1.0_dp)
                                      cdk(3) = dcmplx(d2k*zvec(kb2(jb,ja)),d2ks*zvec(kb2(jb,ja)))*dcmplx(0.0_dp,1.0_dp)
!                                       
                                      if (ja.gt.jb) then
                                        derv2dk(1:3,isxyz(ibet,jb),isxyzl(ialp)) = derv2dk(1:3,isxyz(ibet,jb),isxyzl(ialp)) +  &
                                            cdk(1:3)
                                      elseif (ja.lt.jb) then
                                        derv2dk(1:3,isxyz(ibet,jb),isxyzl(ialp)) = derv2dk(1:3,isxyz(ibet,jb),isxyzl(ialp)) +  &
                                            conjg(cdk(1:3))
                                      endif
                                    endif
!
                                    if (ja.gt.ib) then
                                      derv2(isxyz(ibet,ib),isxyzl(ialp)) = derv2(isxyz(ibet,ib),isxyzl(ialp)) -  &
                                          d2trm*cosvec(kb2(ib,ja))
                                      dervi(isxyz(ibet,ib),isxyzl(ialp)) = dervi(isxyz(ibet,ib),isxyzl(ialp)) -  &
                                          d2trm*sinvec(kb2(ib,ja))
                                    elseif (ja.lt.ib) then
                                      derv2(isxyz(ibet,ib),isxyzl(ialp)) = derv2(isxyz(ibet,ib),isxyzl(ialp)) -  &
                                          d2trm*cosvec(kb2(ib,ja))
                                      dervi(isxyz(ibet,ib),isxyzl(ialp)) = dervi(isxyz(ibet,ib),isxyzl(ialp)) +  &
                                          d2trm*sinvec(kb2(ib,ja))
                                    endif
!
                                    if (lgroupvelocity.and.ja.ne.ib) then
!                                   
!  Group velocities                     
!                                   
                                      d2k  = d2trm*cosvec(kb2(ib,ja))
                                      d2ks = d2trm*sinvec(kb2(ib,ja))
                                      cdk(1) = dcmplx(d2k*xvec(kb2(ib,ja)),d2ks*xvec(kb2(ib,ja)))*dcmplx(0.0_dp,1.0_dp)
                                      cdk(2) = dcmplx(d2k*yvec(kb2(ib,ja)),d2ks*yvec(kb2(ib,ja)))*dcmplx(0.0_dp,1.0_dp)
                                      cdk(3) = dcmplx(d2k*zvec(kb2(ib,ja)),d2ks*zvec(kb2(ib,ja)))*dcmplx(0.0_dp,1.0_dp)
!                                       
                                      if (ja.gt.ib) then
                                        derv2dk(1:3,isxyz(ibet,ib),isxyzl(ialp)) = derv2dk(1:3,isxyz(ibet,ib),isxyzl(ialp)) -  &
                                            cdk(1:3)
                                      elseif (ja.lt.ib) then
                                        derv2dk(1:3,isxyz(ibet,ib),isxyzl(ialp)) = derv2dk(1:3,isxyz(ibet,ib),isxyzl(ialp)) -  &
                                            conjg(cdk(1:3))
                                      endif
                                    endif
                                  endif
!
                                  if (ivec.ne.jvec) then
                                    if (ib.eq.3) then
                                      if (ia.gt.ib) then
                                        derv2(isxyz(ialp,ia),isxyzl(ibet)) = derv2(isxyz(ialp,ia),isxyzl(ibet)) +  &
                                          d2trm*cosvec(kb2(ia,ib))
                                        dervi(isxyz(ialp,ia),isxyzl(ibet)) = dervi(isxyz(ialp,ia),isxyzl(ibet)) -  &
                                          d2trm*sinvec(kb2(ia,ib))
                                      elseif (ia.lt.ib) then
                                        derv2(isxyz(ialp,ia),isxyzl(ibet)) = derv2(isxyz(ialp,ia),isxyzl(ibet)) +  &
                                          d2trm*cosvec(kb2(ia,ib))
                                        dervi(isxyz(ialp,ia),isxyzl(ibet)) = dervi(isxyz(ialp,ia),isxyzl(ibet)) +  &
                                          d2trm*sinvec(kb2(ia,ib))
                                      endif
!
                                      if (lgroupvelocity.and.ia.ne.ib) then
!                                   
!  Group velocities                     
!                                   
                                        d2k  = d2trm*cosvec(kb2(ia,ib))
                                        d2ks = d2trm*sinvec(kb2(ia,ib))
                                        cdk(1) = dcmplx(d2k*xvec(kb2(ia,ib)),d2ks*xvec(kb2(ia,ib)))*dcmplx(0.0_dp,1.0_dp)
                                        cdk(2) = dcmplx(d2k*yvec(kb2(ia,ib)),d2ks*yvec(kb2(ia,ib)))*dcmplx(0.0_dp,1.0_dp)
                                        cdk(3) = dcmplx(d2k*zvec(kb2(ia,ib)),d2ks*zvec(kb2(ia,ib)))*dcmplx(0.0_dp,1.0_dp)
!                                       
                                        if (ia.gt.ib) then
                                          derv2dk(1:3,isxyz(ialp,ia),isxyzl(ibet)) = derv2dk(1:3,isxyz(ialp,ia),isxyzl(ibet)) +  &
                                            conjg(cdk(1:3))
                                        elseif (ia.lt.ib) then
                                          derv2dk(1:3,isxyz(ialp,ia),isxyzl(ibet)) = derv2dk(1:3,isxyz(ialp,ia),isxyzl(ibet)) +  &
                                            cdk(1:3)
                                        endif
                                      endif
                                    endif
!
                                    if (jb.eq.3) then
                                      if (ia.gt.jb) then
                                        derv2(isxyz(ialp,ia),isxyzl(ibet)) = derv2(isxyz(ialp,ia),isxyzl(ibet)) -  &
                                          d2trm*cosvec(kb2(ia,jb))
                                        dervi(isxyz(ialp,ia),isxyzl(ibet)) = dervi(isxyz(ialp,ia),isxyzl(ibet)) +  &
                                          d2trm*sinvec(kb2(ia,jb))
                                      elseif (ia.lt.jb) then
                                        derv2(isxyz(ialp,ia),isxyzl(ibet)) = derv2(isxyz(ialp,ia),isxyzl(ibet)) -  &
                                          d2trm*cosvec(kb2(ia,jb))
                                        dervi(isxyz(ialp,ia),isxyzl(ibet)) = dervi(isxyz(ialp,ia),isxyzl(ibet)) -  &
                                          d2trm*sinvec(kb2(ia,jb))
                                      endif
!
                                      if (lgroupvelocity.and.ia.ne.jb) then
!
!  Group velocities
!
                                        d2k  = d2trm*cosvec(kb2(ia,jb))
                                        d2ks = d2trm*sinvec(kb2(ia,jb))
                                        cdk(1) = dcmplx(d2k*xvec(kb2(ia,jb)),d2ks*xvec(kb2(ia,jb)))*dcmplx(0.0_dp,1.0_dp)
                                        cdk(2) = dcmplx(d2k*yvec(kb2(ia,jb)),d2ks*yvec(kb2(ia,jb)))*dcmplx(0.0_dp,1.0_dp)
                                        cdk(3) = dcmplx(d2k*zvec(kb2(ia,jb)),d2ks*zvec(kb2(ia,jb)))*dcmplx(0.0_dp,1.0_dp)
!
                                        if (ia.gt.jb) then
                                          derv2dk(1:3,isxyz(ialp,ia),isxyzl(ibet)) = derv2dk(1:3,isxyz(ialp,ia),isxyzl(ibet)) -  &
                                            conjg(cdk(1:3))
                                        elseif (ia.lt.jb) then
                                          derv2dk(1:3,isxyz(ialp,ia),isxyzl(ibet)) = derv2dk(1:3,isxyz(ialp,ia),isxyzl(ibet)) -  &
                                            cdk(1:3)
                                        endif
                                      endif
!
                                      if (ja.gt.jb) then
                                        derv2(isxyz(ialp,ja),isxyzl(ibet)) = derv2(isxyz(ialp,ja),isxyzl(ibet)) +  &
                                          d2trm*cosvec(kb2(ja,jb))
                                        dervi(isxyz(ialp,ja),isxyzl(ibet)) = dervi(isxyz(ialp,ja),isxyzl(ibet)) -  &
                                          d2trm*sinvec(kb2(ja,jb))
                                      elseif (ja.lt.jb) then
                                        derv2(isxyz(ialp,ja),isxyzl(ibet)) = derv2(isxyz(ialp,ja),isxyzl(ibet)) +  &
                                          d2trm*cosvec(kb2(ja,jb))
                                        dervi(isxyz(ialp,ja),isxyzl(ibet)) = dervi(isxyz(ialp,ja),isxyzl(ibet)) +  &
                                          d2trm*sinvec(kb2(ja,jb))
                                      endif
!
                                      if (lgroupvelocity.and.ja.ne.jb) then
!
!  Group velocities
!
                                        d2k  = d2trm*cosvec(kb2(ja,jb))
                                        d2ks = d2trm*sinvec(kb2(ja,jb))
                                        cdk(1) = dcmplx(d2k*xvec(kb2(ja,jb)),d2ks*xvec(kb2(ja,jb)))*dcmplx(0.0_dp,1.0_dp)
                                        cdk(2) = dcmplx(d2k*yvec(kb2(ja,jb)),d2ks*yvec(kb2(ja,jb)))*dcmplx(0.0_dp,1.0_dp)
                                        cdk(3) = dcmplx(d2k*zvec(kb2(ja,jb)),d2ks*zvec(kb2(ja,jb)))*dcmplx(0.0_dp,1.0_dp)
!
                                        if (ja.gt.jb) then
                                          derv2dk(1:3,isxyz(ialp,ja),isxyzl(ibet)) = derv2dk(1:3,isxyz(ialp,ja),isxyzl(ibet)) +  &
                                            conjg(cdk(1:3))
                                        elseif (ja.lt.jb) then
                                          derv2dk(1:3,isxyz(ialp,ja),isxyzl(ibet)) = derv2dk(1:3,isxyz(ialp,ja),isxyzl(ibet)) +  &
                                            cdk(1:3)
                                        endif
                                      endif
                                    endif
!
                                    if (ib.eq.3) then
                                      if (ja.gt.ib) then
                                        derv2(isxyz(ialp,ja),isxyzl(ibet)) = derv2(isxyz(ialp,ja),isxyzl(ibet)) -  &
                                          d2trm*cosvec(kb2(ja,ib))
                                        dervi(isxyz(ialp,ja),isxyzl(ibet)) = dervi(isxyz(ialp,ja),isxyzl(ibet)) +  &
                                          d2trm*sinvec(kb2(ja,ib))
                                      elseif (ja.lt.ib) then
                                        derv2(isxyz(ialp,ja),isxyzl(ibet)) = derv2(isxyz(ialp,ja),isxyzl(ibet)) -  &
                                          d2trm*cosvec(kb2(ja,ib))
                                        dervi(isxyz(ialp,ja),isxyzl(ibet)) = dervi(isxyz(ialp,ja),isxyzl(ibet)) -  &
                                          d2trm*sinvec(kb2(ja,ib))
                                      endif
!
                                      if (lgroupvelocity.and.ja.ne.ib) then
!
!  Group velocities
!
                                        d2k  = d2trm*cosvec(kb2(ja,ib))
                                        d2ks = d2trm*sinvec(kb2(ja,ib))
                                        cdk(1) = dcmplx(d2k*xvec(kb2(ja,ib)),d2ks*xvec(kb2(ja,ib)))*dcmplx(0.0_dp,1.0_dp)
                                        cdk(2) = dcmplx(d2k*yvec(kb2(ja,ib)),d2ks*yvec(kb2(ja,ib)))*dcmplx(0.0_dp,1.0_dp)
                                        cdk(3) = dcmplx(d2k*zvec(kb2(ja,ib)),d2ks*zvec(kb2(ja,ib)))*dcmplx(0.0_dp,1.0_dp)
!
                                        if (ja.gt.ib) then
                                          derv2dk(1:3,isxyz(ialp,ja),isxyzl(ibet)) = derv2dk(1:3,isxyz(ialp,ja),isxyzl(ibet)) -  &
                                            conjg(cdk(1:3))
                                        elseif (ja.lt.ib) then
                                          derv2dk(1:3,isxyz(ialp,ja),isxyzl(ibet)) = derv2dk(1:3,isxyz(ialp,ja),isxyzl(ibet)) -  &
                                            cdk(1:3)
                                        endif
                                      endif
                                    endif
                                  endif
                                enddo
                              enddo
                            enddo
                          enddo
!
!  End of inner loops over atoms and cell vectors
!
                        enddo mnloop2
                      enddo nloop2
                    enddo mmloop2
                  enddo mloop2
                enddo llloop2
              enddo lloop2
            enddo jjloop2
          enddo kloop2
        enddo iiloop2
      enddo jloop2
    enddo iloop2
!
!  End of outer loop over potentials
!
  enddo
!
!  Timing
!
  time2 = g_cpu_time()
  tsix = tsix + time2 - time1
#ifdef TRACE
  call trace_out('sixpd')
#endif
!
  return
  end
