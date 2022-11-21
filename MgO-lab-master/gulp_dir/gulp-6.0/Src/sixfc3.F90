  subroutine sixfc3(maxlim,maxlimcg,fc3)
!
!  Subroutine to calculate third derivatives of six-body potentials.
!  Called from thirdorderfc3.
!
!   9/15 Created from sixp/fourfc3
!  10/15 Initialisation of ind3ptr & ind2ptr corrected
!   1/17 Looping indices for j now triangular for all cases
!   1/17 Cell image indices changed to correct a bug
!   1/17 lmolok reset to pre-loop value for each loop over cell images
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
  use molecule
  use six
  use times,          only : tsix
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
  integer(i4)                                  :: ia
  integer(i4)                                  :: ic
  integer(i4)                                  :: ig
  integer(i4)                                  :: igx
  integer(i4)                                  :: igxyz
  integer(i4)                                  :: ii
  integer(i4)                                  :: ii1
  integer(i4)                                  :: ii2
  integer(i4)                                  :: ii3
  integer(i4),                            save :: ind2ptr(15,15)
  integer(i4),                            save :: ind3ptr(15,15,15)
  integer(i4)                                  :: indm
  integer(i4)                                  :: indmj
  integer(i4)                                  :: indmk
  integer(i4)                                  :: indml
  integer(i4)                                  :: indmm
  integer(i4)                                  :: indmn
  integer(i4)                                  :: indvec
  integer(i4)                                  :: isatom(6)
  integer(i4)                                  :: ivec
  integer(i4)                                  :: iveca(2,15)
  integer(i4)                                  :: ix
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
  integer(i4)                                  :: jc
  integer(i4)                                  :: jg
  integer(i4)                                  :: jgx
  integer(i4)                                  :: jj
  integer(i4)                                  :: jmax
  integer(i4)                                  :: jmin
  integer(i4)                                  :: jvec
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
  integer(i4)                                  :: kvec
  integer(i4)                                  :: kx
  integer(i4)                                  :: kxx
  integer(i4)                                  :: kxyz
  integer(i4)                                  :: kyy
  integer(i4)                                  :: kzz
  integer(i4)                                  :: l
  integer(i4)                                  :: lg
  integer(i4)                                  :: lgx
  integer(i4)                                  :: liimax
  integer(i4)                                  :: ll
  integer(i4)                                  :: lmax
  integer(i4)                                  :: lmin
  integer(i4)                                  :: lx
  integer(i4)                                  :: m
  integer(i4)                                  :: mg
  integer(i4)                                  :: mgx
  integer(i4)                                  :: mm
  integer(i4)                                  :: mn
  integer(i4)                                  :: mx
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
  integer(i4)                                  :: ng
  integer(i4)                                  :: nghostcell    ! Number of ghost cells in the full cell
  integer(i4)                                  :: ngx
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
  integer(i4)                                  :: nx
  integer(i4)                                  :: nxx
  integer(i4)                                  :: nyy
  integer(i4)                                  :: nzz
  logical                                      :: l2bonds
  logical                                      :: lbonded
  logical                                      :: lbtyp
  logical,                                save :: lfirst = .true.
  logical                                      :: lghosti
  logical                                      :: lghostj
  logical                                      :: lghostk
  logical                                      :: lghostl
  logical                                      :: lghostm
  logical                                      :: lghostn
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
  real(dp)                                     :: g_cpu_time
  real(dp),                               save :: delta(6,6,15)
  real(dp)                                     :: e1d(15)
  real(dp)                                     :: e2d(120)
  real(dp)                                     :: e3d(680)
  real(dp)                                     :: fc33(3,6,3,6,3,6)    ! Full third derivative matrix
  real(dp)                                     :: eterm
  real(dp)                                     :: oci
  real(dp)                                     :: ocj
  real(dp)                                     :: ock
  real(dp)                                     :: ocl 
  real(dp)                                     :: ocm 
  real(dp)                                     :: ocn 
  real(dp)                                     :: ofct 
  real(dp)                                     :: r212
  real(dp)                                     :: r312
  real(dp)                                     :: r412
  real(dp)                                     :: r522
  real(dp)                                     :: r622
  real(dp)                                     :: rc(3,6,15)
  real(dp)                                     :: rksix
  real(dp)                                     :: rko
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
  call trace_in('sixfc3')
#endif
!
  time1 = g_cpu_time()
!
!  Find number of ghost cells
!
  nghostcell = nsuperghost(1,ncf)*nsuperghost(2,ncf)*nsuperghost(3,ncf)
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
  if (lfirst) then
    lfirst = .false.
!
!  ind2ptr
!
    indvec = 0
    do ivec = 1,15
      do jvec = ivec,15
        indvec = indvec + 1
        ind2ptr(jvec,ivec) = indvec
        ind2ptr(ivec,jvec) = indvec
      enddo
    enddo
!
!  ind3ptr
!
    indvec = 0
    do ivec = 1,15
      do jvec = ivec,15
        do kvec = jvec,15
          indvec = indvec + 1
          ind3ptr(kvec,jvec,ivec) = indvec
          ind3ptr(jvec,kvec,ivec) = indvec
          ind3ptr(kvec,ivec,jvec) = indvec
          ind3ptr(ivec,kvec,jvec) = indvec
          ind3ptr(jvec,ivec,kvec) = indvec
          ind3ptr(ivec,jvec,kvec) = indvec
        enddo
      enddo
    enddo
!
!  Set up delta term array
!
    delta(1:6,1:6,1:15) = 0.0_dp
!
    indvec = 0
    do ivec = 1,5
      do jvec = ivec+1,6
        indvec = indvec + 1
        delta(ivec,ivec,indvec) = 1.0_dp
        delta(ivec,jvec,indvec) = -1.0_dp
        delta(jvec,ivec,indvec) = -1.0_dp
        delta(jvec,jvec,indvec) = 1.0_dp
      enddo
    enddo
  endif
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
    do 10 i = 1,numat
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
        goto 10
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
      sxyz(1,1) = xclat(i)
      sxyz(2,1) = yclat(i)
      sxyz(3,1) = zclat(i)
!
!  Is i a ghost atom?
!
      ix = 3*(i-1) + 1
      lghosti = (mod(i-1,nghostcell).eq.0.and.ni.le.maxele)
      if (lghosti) then
        ig = (i - 1)/nghostcell + 1
        igx = 3*(ig-1) + 1
      endif
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
!  Set j looping indices
!
      jmin = 1
      jmax = i
!********************************
!  Loop over middle site 2 / j  *
!********************************
      do 20 j = jmin,jmax
        nj = nat(j)
        ntypj = nftype(j)
        ocj = occuf(j)
!
!  Check j is allowed for n
!
        if (.not.lmatch(nj,ntypj,nt2,ntyp2,.true.)) goto 20
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
        if (lintra_only.and..not.lmolok) goto 20
        if (lbtyp.and..not.lmolok) goto 20
        xc2t = xclat(j)
        yc2t = yclat(j)
        zc2t = zclat(j)
        x21t = xc2t - sxyz(1,1)
        y21t = yc2t - sxyz(2,1)
        z21t = zc2t - sxyz(3,1)
!
!  Is j a ghost atom?
!
        jx = 3*(j-1) + 1
        lghostj = (mod(j-1,nghostcell).eq.0.and.nj.le.maxele)
        if (lghostj) then
          jg = (j - 1)/nghostcell + 1 
          jgx = 3*(jg-1) + 1 
        endif
!
        lmolokij = lmolok
!
!  Check r21 is OK
!  Loop over cell vectors
!
        do 120 ii = 1,iimax
          x21 = x21t + xvec1cell(ii)
          y21 = y21t + yvec1cell(ii)
          z21 = z21t + zvec1cell(ii)
          r212 = x21*x21 + y21*y21 + z21*z21
          if (r212.lt.1d-10) goto 120
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
              if (linter_only) goto 120
              if (lbtyp) then
                call bonded(lbonded,l2bonds,nbtypeij,nbtypeij2,i,j,0_i4,0_i4,0_i4)
                if (.not.lbonded) goto 120
!
!  Check central bond type for correct order
!
                if (n6botype(1,np).gt.0) then
                  if (n6botype(1,np).ne.nbtypeij) goto 120
                endif
                if (n6botype(2,np).gt.0) then
                  if (n6botype(2,np).ne.nbtypeij2) goto 120
                endif
              endif
            else
              call lintoijk(ixx,iyy,izz,ii,imaxl,jmaxl,kmaxl)
              if (lbtyp) then
                call bonded(lbonded,l2bonds,nbtypeij,nbtypeij2,i,j,ixx,iyy,izz)
                if (.not.lbonded) goto 120
!
!  Check central bond type for correct order
!
                if (n6botype(1,np).gt.0) then
                  if (n6botype(1,np).ne.nbtypeij) goto 120
                endif
                if (n6botype(2,np).gt.0) then
                  if (n6botype(2,np).ne.nbtypeij2) goto 120
                endif
                lsamemol = (lbonded.or.l2bonds)
              else
                lsamemol = .false.
              endif
              if (.not.lsamemol) then
                call samemol(lsamemol,nmi,ixx,iyy,izz,ixj,iyj,izj)
              endif
              if (lintra_only.and..not.lsamemol) goto 120
              if (linter_only.and.lsamemol) goto 120
            endif
          endif
!
!  Distance checking
!
          if (r212.gt.tr1.and.(.not.lbtyp.or..not.lbonded)) goto 120
!
!  j has been accepted
!
          isatom(2) = j
          sxyz(1,2) = x21 + sxyz(1,1)
          sxyz(2,2) = y21 + sxyz(2,1)
          sxyz(3,2) = z21 + sxyz(3,1)
!*****************************************
!  Loop over first end site bonded to i  *
!*****************************************
          do 30 k = 1,numat
            nk = nat(k)
            ntypk = nftype(k)
            ock = occuf(k)
!
!  Check k is allowed for n
!
            if (.not.lmatch(nk,ntypk,nt3,ntyp3,.true.)) goto 30
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
            if (lintra_only.and..not.lmolok) goto 30
            if (lbtyp.and..not.lmolok) goto 30
            xc3t = xclat(k)
            yc3t = yclat(k)
            zc3t = zclat(k)
            x31t = xc3t - sxyz(1,1)
            y31t = yc3t - sxyz(2,1)
            z31t = zc3t - sxyz(3,1)
!
!  Is k a ghost atom?
!
            kx = 3*(k-1) + 1
            lghostk = (mod(k-1,nghostcell).eq.0.and.nk.le.maxele)
            if (lghostk) then
              kg = (k - 1)/nghostcell + 1
              kgx = 3*(kg-1) + 1
            endif
!
            lmolokik = lmolok
!
!  Check r31 is OK
!  Loop over cell vectors
!
            do 130 jj = 1,iimax
              x31 = x31t + xvec1cell(jj)
              y31 = y31t + yvec1cell(jj)
              z31 = z31t + zvec1cell(jj)
              r312 = x31*x31 + y31*y31 + z31*z31
              if (r312.lt.1d-10) goto 130
!
!  Prevent atoms i and k being the same atom
!
              if (k.eq.i.and.jj.eq.iimid) goto 130
!
!  Prevent atoms j and k being the same atom
!
              if (k.eq.j.and.jj.eq.ii) goto 130
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
                  if (linter_only) goto 130
                  if (lbtyp) then
                    call bonded(lbonded,l2bonds,nbtypeik,nbtypeik2,i,k,0_i4,0_i4,0_i4)
                    if (.not.lbonded) goto 130
                  endif
                else
                  call lintoijk(jxx,jyy,jzz,jj,imaxl,jmaxl,kmaxl)
                  if (lbtyp) then
                    call bonded(lbonded,l2bonds,nbtypeik,nbtypeik2,i,k,jxx,jyy,jzz)
                    if (.not.lbonded) goto 130
                    lsamemol = (lbonded.or.l2bonds)
                  else
                    lsamemol = .false.
                  endif
                  if (.not.lsamemol) then
                    call samemol(lsamemol,nmi,jxx,jyy,jzz,ixk,iyk,izk)
                  endif
                  if (lintra_only.and..not.lsamemol) goto 130
                  if (linter_only.and.lsamemol) goto 130
                endif
              endif
!
!  Distance checking
!
              if (r312.gt.tr2.and.(.not.lbtyp.or..not.lbonded)) goto 130
!
!  k has been accepted
!
              isatom(3) = k
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
              do 40 l = lmin,lmax
                nl = nat(l)
                ntypl = nftype(l)
                ocl = occuf(l)
!
!  Check l is allowed for n
!
                if (.not.lmatch(nl,ntypl,nt4,ntyp4,.true.)) goto 40
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
                if (lintra_only.and..not.lmolok) goto 40
                if (lbtyp.and..not.lmolok) goto 40
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
!  Is l a ghost atom?
!  
                lx = 3*(l-1) + 1
                lghostl = (mod(l-1,nghostcell).eq.0.and.nl.le.maxele)
                if (lghostl) then
                  lg = (l - 1)/nghostcell + 1
                  lgx = 3*(lg-1) + 1
                endif
!
                lmolokil = lmolok
!
!  Check r41 is OK
!  Loop over cell vectors
!
                do 140 ll = 1,liimax
                  x41 = x41t + xvec1cell(ll)
                  y41 = y41t + yvec1cell(ll)
                  z41 = z41t + zvec1cell(ll)
                  r412 = x41*x41 + y41*y41 + z41*z41
                  if (r412.lt.1d-10) goto 140
!
!  Prevent atoms i and l being the same atom
!
                  if (l.eq.i.and.ll.eq.iimid) goto 140
!
!  Prevent atoms j and l being the same atom
!
                  if (l.eq.j.and.ll.eq.ii) goto 140
!
!  Prevent atoms k and l being the same atom
!
                  if (l.eq.k.and.ll.eq.jj) goto 140
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
                      if (linter_only) goto 140
                      if (lbtyp) then
                        call bonded(lbonded,l2bonds,nbtypeil,nbtypeil2,i,l,0_i4,0_i4,0_i4)
                        if (.not.lbonded) goto 140
                      endif
                    else
                      call lintoijk(kxx,kyy,kzz,ll,imaxl,jmaxl,kmaxl)
                      if (lbtyp) then
                        call bonded(lbonded,l2bonds,nbtypeil,nbtypeil2,i,l,kxx,kyy,kzz)
                        if (.not.lbonded) goto 140
                        lsamemol = (lbonded.or.l2bonds)
                      else
                        lsamemol = .false.
                      endif
                      if (.not.lsamemol) then
                        call samemol(lsamemol,nmi,kxx,kyy,kzz,ixl,iyl,izl)
                      endif
                      if (lintra_only.and..not.lsamemol) goto 140
                      if (linter_only.and.lsamemol) goto 140
                    endif
                  endif
!
!  Distance checking
!
                  if (r412.gt.tr3.and.(.not.lbtyp.or..not.lbonded)) goto 140
!
!  l has been accepted
!
                  isatom(4) = l
                  sxyz(1,4) = x41 + sxyz(1,1)
                  sxyz(2,4) = y41 + sxyz(2,1)
                  sxyz(3,4) = z41 + sxyz(3,1)
!*****************************************
!  Loop over first end site bonded to j  *
!*****************************************
                  do 50 m = 1,numat
                    nm = nat(m)
                    ntypm = nftype(m)
                    ocm = occuf(m)
!
!  Check m is allowed for n
!
                    if (.not.lmatch(nm,ntypm,nt5,ntyp5,.true.)) goto 50
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
                    if (lintra_only.and..not.lmolok) goto 50
                    if (lbtyp.and..not.lmolok) goto 50
                    xc5t = xclat(m)
                    yc5t = yclat(m)
                    zc5t = zclat(m)
                    x52t = xc5t - sxyz(1,2)
                    y52t = yc5t - sxyz(2,2)
                    z52t = zc5t - sxyz(3,2)
!
!  Is m a ghost atom?
! 
                    mx = 3*(m-1) + 1
                    lghostm = (mod(m-1,nghostcell).eq.0.and.nm.le.maxele)
                    if (lghostm) then
                      mg = (m - 1)/nghostcell + 1
                      mgx = 3*(mg-1) + 1
                    endif
!
                    lmolokjm = lmolok
!
!  Check r52 is OK
!  Loop over cell vectors
!
                    do 150 mm = 1,iimax
                      x52 = x52t + xvec1cell(mm)
                      y52 = y52t + yvec1cell(mm)
                      z52 = z52t + zvec1cell(mm)
                      r522 = x52*x52 + y52*y52 + z52*z52
                      if (r522.lt.1d-10) goto 150
!
!  Prevent atoms i and m being the same atom
!
                      if (m.eq.i.and.mm.eq.iimid) goto 150
!
!  Prevent atoms j and m being the same atom
!
                      if (m.eq.j.and.mm.eq.ii) goto 150
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
                          if (linter_only) goto 150
                          if (lbtyp) then
                            call bonded(lbonded,l2bonds,nbtypejm,nbtypejm2,j,m,0_i4,0_i4,0_i4)
                            if (.not.lbonded) goto 150
                          endif
                        else
                          call lintoijk(mxx,myy,mzz,mm,imaxl,jmaxl,kmaxl)
                          if (lbtyp) then
                            call bonded(lbonded,l2bonds,nbtypejm,nbtypejm2,j,m,mxx-ixx,myy-iyy,mzz-izz)
                            if (.not.lbonded) goto 150
                            lsamemol = (lbonded.or.l2bonds)
                          else
                            lsamemol = .false.
                          endif
                          if (.not.lsamemol) then
                            call samemol(lsamemol,nmi,jxx,jyy,jzz,ixm,iym,izm)
                          endif
                          if (lintra_only.and..not.lsamemol) goto 150
                          if (linter_only.and.lsamemol) goto 150
                        endif
                      endif
!
!  Distance checking
!
                      if (r522.gt.tr4.and.(.not.lbtyp.or..not.lbonded)) goto 150
!
!  m has been accepted
!
                      isatom(5) = m
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
                      do 60 n = nmin,nmax
                        nn = nat(n)
                        ntypn = nftype(n)
                        ocn = occuf(n)
!
!  Check n is allowed for n
!
                        if (.not.lmatch(nn,ntypn,nt6,ntyp6,.true.)) goto 60
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
                        if (lintra_only.and..not.lmolok) goto 60
                        if (lbtyp.and..not.lmolok) goto 60
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
!  Is n a ghost atom?
!
                        nx = 3*(n-1) + 1
                        lghostn = (mod(n-1,nghostcell).eq.0.and.nn.le.maxele)
                        if (lghostn) then
                          ng = (n - 1)/nghostcell + 1
                          ngx = 3*(ng-1) + 1
                        endif
!
                        lmolokjn = lmolok
!
!  Check r62 is OK
!  Loop over cell vectors
!
                        do 160 mn = 1,niimax
                          x62 = x62t + xvec1cell(mn)
                          y62 = y62t + yvec1cell(mn)
                          z62 = z62t + zvec1cell(mn)
                          r622 = x62*x62 + y62*y62 + z62*z62
                          if (r622.lt.1d-10) goto 160
!
!  Prevent atoms i and n being the same atom
!
                          if (n.eq.i.and.mn.eq.iimid) goto 160
!
!  Prevent atoms j and n being the same atom
!
                          if (n.eq.j.and.mn.eq.ii) goto 160
!
!  Prevent atoms m and n being the same atom
!
                          if (n.eq.m.and.mn.eq.mm) goto 160
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
                              if (linter_only) goto 160
                              if (lbtyp) then
                                call bonded(lbonded,l2bonds,nbtypejn,nbtypejn2,j,n,0_i4,0_i4,0_i4)
                                if (.not.lbonded) goto 160
                              endif
                            else
                              call lintoijk(nxx,nyy,nzz,mn,imaxl,jmaxl,kmaxl)
                              if (lbtyp) then
                                call bonded(lbonded,l2bonds,nbtypejn,nbtypejn2,j,n,nxx-ixx,nyy-iyy,nzz-izz)
                                if (.not.lbonded) goto 160
                                lsamemol = (lbonded.or.l2bonds)
                              else
                                lsamemol = .false.
                              endif
                              if (.not.lsamemol) then
                                call samemol(lsamemol,nmi,jxx,jyy,jzz,ixn,iyn,izn)
                              endif
                              if (lintra_only.and..not.lsamemol) goto 160
                              if (linter_only.and.lsamemol) goto 160
                            endif
                          endif
!
!  Distance checking
!
                          if (r622.gt.tr5.and.(.not.lbtyp.or..not.lbonded)) goto 160
!
!  n has been accepted
!
                          isatom(6) = n
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
                          call sixbody(np,nsixtype,sdist,eterm,e1d,e2d,e3d,rko,.true.,.true.,.true.)
!
                          indvec = 0
                          do ivec = 1,5
                            do jvec = ivec+1,6
                              indvec = indvec + 1
                              xvec(indvec) = svec(1,jvec,ivec)
                              yvec(indvec) = svec(2,jvec,ivec)
                              zvec(indvec) = svec(3,jvec,ivec)
                            enddo
                          enddo
!********************************
!  Calculate third derivatives  *
!********************************
!
!  Initialise full third derivative matrix
!
                          fc33(1:3,1:6,1:3,1:6,1:3,1:6) = 0.0_dp
!
!  Set up array of Cartesian components for vectors
!
                          rc(1:3,1:6,1:15) = 0.0_dp
!
                          do ivec = 1,15
                            ja = iveca(1,ivec)
                            ia = iveca(2,ivec)
                            rc(1,ja,ivec) = - xvec(indvec)
                            rc(2,ja,ivec) = - yvec(indvec)
                            rc(3,ja,ivec) = - zvec(indvec)
                            rc(1,ia,ivec) = xvec(indvec)
                            rc(2,ia,ivec) = yvec(indvec)
                            rc(3,ia,ivec) = zvec(indvec)
                          enddo
!
!  Loop over atoms and Cartesian components to compute third derivatives
!
                          do ia = 1,6   ! First atom
                            do ic = 1,3   ! First atom Cartesian component
                              do ja = 1,6   ! Second atom
                                do jc = 1,3   ! Second atom Cartesian component
                                  do ka = 1,6   ! Third atom
                                    do kc = 1,3   ! Third atom Cartesian component
                                      do ii1 = 1,15  ! Vectors between atoms
!
                                        do ii2 = 1,15  ! Vectors between atoms
                                          do ii3 = 1,15  ! Vectors between atoms
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
                              do ja = 1,6
                                if (ja.eq.1) then
                                  jxyz = ix - 1
                                elseif (ja.eq.2) then
                                  jxyz = jx - 1
                                elseif (ja.eq.3) then
                                  jxyz = kx - 1
                                elseif (ja.eq.4) then
                                  jxyz = lx - 1
                                elseif (ja.eq.5) then
                                  jxyz = mx - 1
                                else
                                  jxyz = nx - 1
                                endif
                                do jc = 1,3
                                  jxyz = jxyz + 1
                                  do ka = 1,6
                                    if (ka.eq.1) then
                                      kxyz = ix - 1
                                    elseif (ka.eq.2) then
                                      kxyz = jx - 1
                                    elseif (ka.eq.3) then
                                      kxyz = kx - 1
                                    elseif (ka.eq.4) then
                                      kxyz = lx - 1
                                    elseif (ka.eq.5) then
                                      kxyz = mx - 1
                                    else
                                      kxyz = nx - 1
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
                              do ja = 1,6
                                if (ja.eq.1) then
                                  jxyz = ix - 1
                                elseif (ja.eq.2) then
                                  jxyz = jx - 1
                                elseif (ja.eq.3) then
                                  jxyz = kx - 1
                                elseif (ja.eq.4) then
                                  jxyz = lx - 1
                                elseif (ja.eq.5) then
                                  jxyz = mx - 1
                                else
                                  jxyz = nx - 1
                                endif
                                do jc = 1,3
                                  jxyz = jxyz + 1
                                  do ka = 1,6
                                    if (ka.eq.1) then
                                      kxyz = ix - 1
                                    elseif (ka.eq.2) then
                                      kxyz = jx - 1
                                    elseif (ka.eq.3) then
                                      kxyz = kx - 1
                                    elseif (ka.eq.4) then
                                      kxyz = lx - 1
                                    elseif (ka.eq.5) then
                                      kxyz = mx - 1
                                    else
                                      kxyz = nx - 1
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
                              do ja = 1,6 
                                if (ja.eq.1) then
                                  jxyz = ix - 1 
                                elseif (ja.eq.2) then
                                  jxyz = jx - 1 
                                elseif (ja.eq.3) then
                                  jxyz = kx - 1
                                elseif (ja.eq.4) then
                                  jxyz = lx - 1
                                elseif (ja.eq.5) then
                                  jxyz = mx - 1
                                else
                                  jxyz = nx - 1 
                                endif
                                do jc = 1,3 
                                  jxyz = jxyz + 1 
                                  do ka = 1,6 
                                    if (ka.eq.1) then
                                      kxyz = ix - 1 
                                    elseif (ka.eq.2) then
                                      kxyz = jx - 1 
                                    elseif (ka.eq.3) then
                                      kxyz = kx - 1
                                    elseif (ka.eq.4) then
                                      kxyz = lx - 1
                                    elseif (ka.eq.5) then
                                      kxyz = mx - 1
                                    else
                                      kxyz = nx - 1 
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
                              do ja = 1,6 
                                if (ja.eq.1) then
                                  jxyz = ix - 1 
                                elseif (ja.eq.2) then
                                  jxyz = jx - 1 
                                elseif (ja.eq.3) then
                                  jxyz = kx - 1
                                elseif (ja.eq.4) then
                                  jxyz = lx - 1
                                elseif (ja.eq.5) then
                                  jxyz = mx - 1
                                else
                                  jxyz = nx - 1 
                                endif
                                do jc = 1,3 
                                  jxyz = jxyz + 1 
                                  do ka = 1,6 
                                    if (ka.eq.1) then
                                      kxyz = ix - 1 
                                    elseif (ka.eq.2) then
                                      kxyz = jx - 1 
                                    elseif (ka.eq.3) then
                                      kxyz = kx - 1
                                    elseif (ka.eq.4) then
                                      kxyz = lx - 1
                                    elseif (ka.eq.5) then
                                      kxyz = mx - 1
                                    else
                                      kxyz = nx - 1 
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
                          if (lghostm) then
                            igxyz = mgx - 1 
                            do ic = 1,3 
                              igxyz = igxyz + 1 
                              do ja = 1,6 
                                if (ja.eq.1) then
                                  jxyz = ix - 1 
                                elseif (ja.eq.2) then
                                  jxyz = jx - 1 
                                elseif (ja.eq.3) then
                                  jxyz = kx - 1
                                elseif (ja.eq.4) then
                                  jxyz = lx - 1
                                elseif (ja.eq.5) then
                                  jxyz = mx - 1
                                else
                                  jxyz = nx - 1 
                                endif
                                do jc = 1,3 
                                  jxyz = jxyz + 1 
                                  do ka = 1,6 
                                    if (ka.eq.1) then
                                      kxyz = ix - 1 
                                    elseif (ka.eq.2) then
                                      kxyz = jx - 1 
                                    elseif (ka.eq.3) then
                                      kxyz = kx - 1
                                    elseif (ka.eq.4) then
                                      kxyz = lx - 1
                                    elseif (ka.eq.5) then
                                      kxyz = mx - 1
                                    else
                                      kxyz = nx - 1 
                                    endif
                                    do kc = 1,3 
                                      kxyz = kxyz + 1 
                                      fc3(kxyz,jxyz,igxyz) = fc3(kxyz,jxyz,igxyz) + fc33(kc,ka,jc,ja,ic,5)
                                    enddo
                                  enddo
                                enddo
                              enddo
                            enddo
                          endif
!
                          if (lghostn) then
                            igxyz = ngx - 1 
                            do ic = 1,3 
                              igxyz = igxyz + 1 
                              do ja = 1,6 
                                if (ja.eq.1) then
                                  jxyz = ix - 1 
                                elseif (ja.eq.2) then
                                  jxyz = jx - 1 
                                elseif (ja.eq.3) then
                                  jxyz = kx - 1
                                elseif (ja.eq.4) then
                                  jxyz = lx - 1
                                elseif (ja.eq.5) then
                                  jxyz = mx - 1
                                else
                                  jxyz = nx - 1 
                                endif
                                do jc = 1,3 
                                  jxyz = jxyz + 1 
                                  do ka = 1,6 
                                    if (ka.eq.1) then
                                      kxyz = ix - 1 
                                    elseif (ka.eq.2) then
                                      kxyz = jx - 1 
                                    elseif (ka.eq.3) then
                                      kxyz = kx - 1
                                    elseif (ka.eq.4) then
                                      kxyz = lx - 1
                                    elseif (ka.eq.5) then
                                      kxyz = mx - 1
                                    else
                                      kxyz = nx - 1 
                                    endif
                                    do kc = 1,3 
                                      kxyz = kxyz + 1 
                                      fc3(kxyz,jxyz,igxyz) = fc3(kxyz,jxyz,igxyz) + fc33(kc,ka,jc,ja,ic,6)
                                    enddo
                                  enddo
                                enddo
                              enddo
                            enddo
                          endif
!
!  End of inner loops over atoms and cell vectors
!
160                       continue
60                     continue
150                   continue
50                 continue
140               continue
40             continue
130           continue
30         continue
120       continue
20     continue
10   continue
!
!  End of outer loops
!
  enddo
!
!  Timing
!
  time2 = g_cpu_time()
  tsix = tsix + time2 - time1
#ifdef TRACE
  call trace_out('sixfc3')
#endif
!
  return
  end
