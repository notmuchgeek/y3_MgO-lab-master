  subroutine sixsd(esix,lgrad1)
!
!  Subroutine for six-body potentials with symmetry - first derivatives only
!
!   7/06 Created based on sixsd2
!   2/07 Bonding types and test added
!   5/07 Position of goto for freezing corrected
!  12/07 Unused variables removed
!   5/12 Atomic stresses added
!   5/12 Atomic stresses removed for routines involving symmetry
!   1/17 lPrintSix added
!   1/17 Cell image indices changed to correct a bug
!   1/17 lmolok reset to pre-loop value for each loop over cell images
!   1/17 Typo in adding to rstrd in loops to rstrdl corrected
!   1/17 Corrections made to end atom code
!   2/18 Trace added
!   9/18 Strain module introduced
!   9/18 Call to sixstrterms replaced with more general realstrterms
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
  use control,        only : lrigid
  use current
  use derivatives
  use iochannels,     only : ioout
  use m_strain,       only : realstrterms
  use mdlogic
  use molecule
  use optimisation
  use parallel,       only : ioproc
  use six
  use symmetry
  use times,          only : tsix
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  real(dp), intent(inout)                      :: esix
  logical,  intent(in)                         :: lgrad1
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: iai
  integer(i4)                                  :: ii
  integer(i4)                                  :: indm
  integer(i4)                                  :: indmj
  integer(i4)                                  :: indmk
  integer(i4)                                  :: indml
  integer(i4)                                  :: indmm
  integer(i4)                                  :: indmn
  integer(i4)                                  :: indvec
  integer(i4)                                  :: ivec
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
  integer(i4)                                  :: jj
  integer(i4)                                  :: jvec
  integer(i4)                                  :: jxx
  integer(i4)                                  :: jyy
  integer(i4)                                  :: jzz
  integer(i4)                                  :: k
  integer(i4)                                  :: kl
  integer(i4)                                  :: ks
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
  integer(i4)                                  :: nbtypejk
  integer(i4)                                  :: nbtypejl
  integer(i4)                                  :: nbtypejm
  integer(i4)                                  :: nbtypejn
  integer(i4)                                  :: nbtypekm
  integer(i4)                                  :: nbtypekn
  integer(i4)                                  :: nbtypeij2
  integer(i4)                                  :: nbtypeik2
  integer(i4)                                  :: nbtypeil2
  integer(i4)                                  :: nbtypejk2
  integer(i4)                                  :: nbtypejl2
  integer(i4)                                  :: nbtypejm2
  integer(i4)                                  :: nbtypejn2
  integer(i4)                                  :: nbtypekm2
  integer(i4)                                  :: nbtypekn2
  integer(i4)                                  :: neqi
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
  integer(i4)                                  :: nt1
  integer(i4)                                  :: nt2
  integer(i4)                                  :: nt3
  integer(i4)                                  :: nt4
  integer(i4)                                  :: nt5
  integer(i4)                                  :: nt6
  integer(i4)                                  :: ntyp1
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
  logical                                      :: lmolokjk
  logical                                      :: lmolokjl
  logical                                      :: lmolokjm
  logical                                      :: lmolokjn
  logical                                      :: lmolokkm
  logical                                      :: lmolokkn
  logical                                      :: lneedmol
  logical                                      :: lopi
  logical                                      :: lopj
  logical                                      :: lopk
  logical                                      :: lopl
  logical                                      :: lopm
  logical                                      :: lopn
  logical                                      :: lsamemol
  logical                                      :: lsg1
  logical                                      :: ltsym34
  logical                                      :: ltsym56
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: e1d(15)
  real(dp)                                     :: e2d(1)
  real(dp)                                     :: e3d(1)
  real(dp)                                     :: esixl
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
  real(dp)                                     :: r512
  real(dp)                                     :: r522
  real(dp)                                     :: r612
  real(dp)                                     :: r622
  real(dp)                                     :: rksix
  real(dp)                                     :: rko
  real(dp)                                     :: dr2ds(6,15)
  real(dp)                                     :: d2r2dx2(3,3,15)
  real(dp)                                     :: d2r2ds2(6,6,15)
  real(dp)                                     :: d2r2dsdx(6,3,15)
  real(dp)                                     :: rstrdl(6)
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
  real(dp)                                     :: x51
  real(dp)                                     :: y51
  real(dp)                                     :: z51
  real(dp)                                     :: x52
  real(dp)                                     :: y52
  real(dp)                                     :: z52
  real(dp)                                     :: x52t
  real(dp)                                     :: y52t
  real(dp)                                     :: z52t
  real(dp)                                     :: x61
  real(dp)                                     :: y61
  real(dp)                                     :: z61
  real(dp)                                     :: x62
  real(dp)                                     :: y62
  real(dp)                                     :: z62
  real(dp)                                     :: x62t
  real(dp)                                     :: y62t
  real(dp)                                     :: z62t
  real(dp)                                     :: sdist(15)
  real(dp)                                     :: svec(3,6,6)
  real(dp)                                     :: sxyz(3,6)
  real(dp)                                     :: xcom(15)
  real(dp)                                     :: ycom(15)
  real(dp)                                     :: zcom(15)
  real(dp)                                     :: xcomi
  real(dp)                                     :: ycomi
  real(dp)                                     :: zcomi
  real(dp)                                     :: xv6(15)
  real(dp)                                     :: yv6(15)
  real(dp)                                     :: zv6(15)
#ifdef TRACE
  call trace_in('sixsd')
#endif
!
  time1 = g_cpu_time()
  lsg1 = (lstr.and.lgrad1)
  esixl = 0.0_dp
  if (lstr) then
    do kl = 1,nstrains
      rstrdl(kl) = 0.0_dp
    enddo
  endif
!
!  Initialise centre of mass arrays to zero since there should be no intermolecular terms here
!
  xcom(1:15) = 0.0_dp
  ycom(1:15) = 0.0_dp
  zcom(1:15) = 0.0_dp
!
!  Opening banner for energy decomposition
!
  if (lPrintSix) then
    call mpbarrier
    if (ioproc) then
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''  Six  :  Mid Atom 1 + End 3/4    Mid Atom 2 + End 5/6     Sixbody energy (eV)  '')')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
    endif
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
    lbtyp = (mmsexc(np).eq.1)
    rksix = sixk(np)
    lintra_only = (lsintra(np).and..not.lsinter(np))
    linter_only = (lsinter(np).and..not.lsintra(np))
    lneedmol = (lintra_only.or.linter_only.or.lbtyp.or.lrigid)
!********************************
!  Loop over middle site 1 / i  *
!********************************
    do 10 iai = 1,nasym
      i = nrela2f(iai)
      ni = iatn(iai)
      ntypi = natype(iai)
      neqi = neqv(iai)
      oci = occua(iai)*dble(neqi)
      lopi = (lopf(iai).or..not.lfreeze)
!
!  Check i is allowed for n
!
      if (lmatch(ni,ntypi,nsspec1(np),nsptyp1(np),.true.)) then
        nt1 = nsspec1(np)
        nt2 = nsspec2(np)
        ntyp1 = nsptyp1(np)
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
        nt1 = nsspec2(np)
        nt2 = nsspec1(np)
        ntyp1 = nsptyp2(np)
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
!
!  Set COM coordinates
!
        if (lrigid.and.nmi.gt.0) then
          xcomi = molxyz(1,natinmol(i),nmi)
          ycomi = molxyz(2,natinmol(i),nmi)
          zcomi = molxyz(3,natinmol(i),nmi)
        else
          xcomi = 0.0_dp
          ycomi = 0.0_dp
          zcomi = 0.0_dp
        endif
      else
        xcomi = 0.0_dp
        ycomi = 0.0_dp
        zcomi = 0.0_dp
      endif
!********************************
!  Loop over middle site 2 / j  *
!********************************
      do 20 j = 1,numat
        nj = nat(j)
        ntypj = nftype(j)
        ocj = occuf(j)
        lopj = (lopf(nrelf2a(j)).or..not.lfreeze)
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
!
!  Set COM coordinates
!
          if (lrigid.and.nmj.gt.0.and.lmolok) then
            xcom(1) = molxyz(1,natinmol(j),nmj) - xcomi
            ycom(1) = molxyz(2,natinmol(j),nmj) - ycomi
            zcom(1) = molxyz(3,natinmol(j),nmj) - zcomi
          else
            xcom(1) = - xcomi
            ycom(1) = - ycomi
            zcom(1) = - zcomi
          endif
        else
          lmolok = .false.
          xcom(1) = - xcomi
          ycom(1) = - ycomi
          zcom(1) = - zcomi
        endif
!
!  Check for intra and but not in same molecule
!
        if (lintra_only.and..not.lmolok) goto 20
        if (lbtyp.and..not.lmolok) goto 20
        x21t = xclat(j) - sxyz(1,1)
        y21t = yclat(j) - sxyz(2,1)
        z21t = zclat(j) - sxyz(3,1)
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
            lopk = (lopf(nrelf2a(k)).or..not.lfreeze)
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
!
!  Set COM coordinates
!
              if (lrigid.and.nmk.gt.0.and.lmolok) then
                xcom(2) = molxyz(1,natinmol(k),nmk) - xcomi
                ycom(2) = molxyz(2,natinmol(k),nmk) - ycomi
                zcom(2) = molxyz(3,natinmol(k),nmk) - zcomi
              else
                xcom(2) = - xcomi
                ycom(2) = - ycomi
                zcom(2) = - zcomi
              endif
              xcom(6) = xcom(2) - xcom(1)
              ycom(6) = ycom(2) - ycom(1)
              zcom(6) = zcom(2) - zcom(1)
            else
              lmolok = .false.
              xcom(2) = - xcomi
              ycom(2) = - ycomi
              zcom(2) = - zcomi
              xcom(6) = xcom(2) - xcom(1)
              ycom(6) = ycom(2) - ycom(1)
              zcom(6) = zcom(2) - zcom(1)
            endif
!
!  Check for intra and but not in same molecule
!
            if (lintra_only.and..not.lmolok) goto 30
            if (lbtyp.and..not.lmolok) goto 30
            x31t = xclat(k) - sxyz(1,1)
            y31t = yclat(k) - sxyz(2,1)
            z31t = zclat(k) - sxyz(3,1)
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
                lopl = (lopf(nrelf2a(l)).or..not.lfreeze)
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
!
!  Set COM coordinates
!
                  if (lrigid.and.nml.gt.0.and.lmolok) then
                    xcom(3) = molxyz(1,natinmol(l),nml) - xcomi
                    ycom(3) = molxyz(2,natinmol(l),nml) - ycomi
                    zcom(3) = molxyz(3,natinmol(l),nml) - zcomi
                  else
                    xcom(3) = - xcomi
                    ycom(3) = - ycomi
                    zcom(3) = - zcomi
                  endif
                  xcom(7) = xcom(3) - xcom(1)
                  ycom(7) = ycom(3) - ycom(1)
                  zcom(7) = zcom(3) - zcom(1)
                  xcom(10) = xcom(3) - xcom(2)
                  ycom(10) = ycom(3) - ycom(2)
                  zcom(10) = zcom(3) - zcom(2)
                else
                  lmolok = .false.
                  xcom(3) = - xcomi
                  ycom(3) = - ycomi
                  zcom(3) = - zcomi
                  xcom(7) = xcom(3) - xcom(1)
                  ycom(7) = ycom(3) - ycom(1)
                  zcom(7) = zcom(3) - zcom(1)
                  xcom(10) = xcom(3) - xcom(2)
                  ycom(10) = ycom(3) - ycom(2)
                  zcom(10) = zcom(3) - zcom(2)
                endif
!
!  Check for intra and but not in same molecule
!
                if (lintra_only.and..not.lmolok) goto 40
                if (lbtyp.and..not.lmolok) goto 40
                x41t = xclat(l) - sxyz(1,1)
                y41t = yclat(l) - sxyz(2,1)
                z41t = zclat(l) - sxyz(3,1)
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
                    lopm = (lopf(nrelf2a(m)).or..not.lfreeze)
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
!
!  Set COM coordinates
!
                      if (lrigid.and.nmm.gt.0.and.lmolok) then
                        xcom(4) = molxyz(1,natinmol(m),nmm) - xcomi
                        ycom(4) = molxyz(2,natinmol(m),nmm) - ycomi
                        zcom(4) = molxyz(3,natinmol(m),nmm) - zcomi
                      else
                        xcom(4) = - xcomi
                        ycom(4) = - ycomi
                        zcom(4) = - zcomi
                      endif
                      xcom(8) = xcom(4) - xcom(1)
                      ycom(8) = ycom(4) - ycom(1)
                      zcom(8) = zcom(4) - zcom(1)
                      xcom(11) = xcom(4) - xcom(2)
                      ycom(11) = ycom(4) - ycom(2)
                      zcom(11) = zcom(4) - zcom(2)
                      xcom(13) = xcom(4) - xcom(3)
                      ycom(13) = ycom(4) - ycom(3)
                      zcom(13) = zcom(4) - zcom(3)
                    else
                      lmolok = .false.
                      xcom(4) = - xcomi
                      ycom(4) = - ycomi
                      zcom(4) = - zcomi
                      xcom(8) = xcom(4) - xcom(1)
                      ycom(8) = ycom(4) - ycom(1)
                      zcom(8) = zcom(4) - zcom(1)
                      xcom(11) = xcom(4) - xcom(2)
                      ycom(11) = ycom(4) - ycom(2)
                      zcom(11) = zcom(4) - zcom(2)
                      xcom(13) = xcom(4) - xcom(3)
                      ycom(13) = ycom(4) - ycom(3)
                      zcom(13) = zcom(4) - zcom(3)
                    endif
!
!  Check for intra and but not in same molecule
!
                    if (lintra_only.and..not.lmolok) goto 50
                    if (lbtyp.and..not.lmolok) goto 50
                    x52t = xclat(m) - sxyz(1,2)
                    y52t = yclat(m) - sxyz(2,2)
                    z52t = zclat(m) - sxyz(3,2)
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
                        lopn = (lopf(nrelf2a(n)).or..not.lfreeze)
!
!  If lfreeze=.true. and no atoms have any variables then skip this four body term
!
                        if (.not.lopi.and..not.lopj.and..not.lopk.and..not.lopl.and..not.lopm.and..not.lopn) goto 60
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
!
!  Set COM coordinates
!
                          if (lrigid.and.nmn.gt.0.and.lmolok) then
                            xcom(5) = molxyz(1,natinmol(n),nmn) - xcomi
                            ycom(5) = molxyz(2,natinmol(n),nmn) - ycomi
                            zcom(5) = molxyz(3,natinmol(n),nmn) - zcomi
                          else
                            xcom(5) = - xcomi
                            ycom(5) = - ycomi
                            zcom(5) = - zcomi
                          endif
                          xcom(9) = xcom(5) - xcom(1)
                          ycom(9) = ycom(5) - ycom(1)
                          zcom(9) = zcom(5) - zcom(1)
                          xcom(12) = xcom(5) - xcom(2)
                          ycom(12) = ycom(5) - ycom(2)
                          zcom(12) = zcom(5) - zcom(2)
                          xcom(14) = xcom(5) - xcom(3)
                          ycom(14) = ycom(5) - ycom(3)
                          zcom(14) = zcom(5) - zcom(3)
                          xcom(15) = xcom(5) - xcom(4)
                          ycom(15) = ycom(5) - ycom(4)
                          zcom(15) = zcom(5) - zcom(4)
                        else
                          lmolok = .false.
                          xcom(5) = - xcomi
                          ycom(5) = - ycomi
                          zcom(5) = - zcomi
                          xcom(9) = xcom(5) - xcom(1)
                          ycom(9) = ycom(5) - ycom(1)
                          zcom(9) = zcom(5) - zcom(1)
                          xcom(12) = xcom(5) - xcom(2)
                          ycom(12) = ycom(5) - ycom(2)
                          zcom(12) = zcom(5) - zcom(2)
                          xcom(14) = xcom(5) - xcom(3)
                          ycom(14) = ycom(5) - ycom(3)
                          zcom(14) = zcom(5) - zcom(3)
                          xcom(15) = xcom(5) - xcom(4)
                          ycom(15) = ycom(5) - ycom(4)
                          zcom(15) = zcom(5) - zcom(4)
                        endif
!
!  Check for intra and but not in same molecule
!
                        if (lintra_only.and..not.lmolok) goto 60
                        if (lbtyp.and..not.lmolok) goto 60
                        x62t = xclat(n) - sxyz(1,2)
                        y62t = yclat(n) - sxyz(2,2)
                        z62t = zclat(n) - sxyz(3,2)
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
                              xv6(indvec) = svec(1,jvec,ivec)
                              yv6(indvec) = svec(2,jvec,ivec)
                              zv6(indvec) = svec(3,jvec,ivec)
                              sdist(indvec) = svec(1,jvec,ivec)**2 + svec(2,jvec,ivec)**2 + svec(3,jvec,ivec)**2
                              sdist(indvec) = sqrt(sdist(indvec))
                            enddo
                          enddo
!
                          ofct = oci*ocj*ock*ocl*ocm*ocn
                          rko = rksix*ofct
                          call sixbody(np,nsixtype,sdist,eterm,e1d,e2d,e3d,rko,lgrad1,.false.,.false.)
                          esixl = esixl + eterm
!
!  Output energy contribution - note the factor of half is included here as per the end of this routine
!
                          if (lPrintSix) then
                            write(ioout,'(4x,6i8,1x,f22.10)') i,k,l,j,m,n,0.5_dp*eterm
                          endif
!***********************************
!  Cross out of plane derivatives  *
!***********************************
!
!  Set up strain products
!
                          if (lsg1) then
                            call realstrterms(ndim,15_i4,15_i4,xv6,yv6,zv6,xcom,ycom,zcom, &
                              dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
                          endif
!***********************
!  Strain derivatives  *
!***********************
                          if (lsg1) then
!
!  First strain derivatives
!
                            do ivec = 1,15
                              do kl = 1,nstrains
                                ks = nstrptr(kl)
                                rstrdl(kl) = rstrdl(kl) + e1d(ivec)*dr2ds(ks,ivec)
                              enddo
                            enddo
                          endif
!*************************
!  Internal derivatives  *
!*************************
                          if (lgrad1.and.lopi) then
                            indvec = 0
                            do jvec = 2,6
                              indvec = indvec + 1
                              xdrv(iai) = xdrv(iai) + svec(1,jvec,1)*e1d(indvec)
                              ydrv(iai) = ydrv(iai) + svec(2,jvec,1)*e1d(indvec)
                              zdrv(iai) = zdrv(iai) + svec(3,jvec,1)*e1d(indvec)
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
!*********************************************************************************************
!  Second search for potentials - only required to compute the derivatives of the end atoms  *
!*********************************************************************************************
    if (lgrad1) then
!************************************************
!  Loop over end site 1 in asymmetric unit / i  *
!************************************************
      do 15 iai = 1,nasym
        i = nrela2f(iai)
        ni = iatn(iai)
        ntypi = natype(iai)
        neqi = neqv(iai)
        oci = occua(iai)*dble(neqi)
        lopi = (lopf(iai).or..not.lfreeze)
!
!  Check i is allowed for n
!
        if (lmatch(ni,ntypi,nsspec3(np),nsptyp3(np),.true.)) then
          nt1 = nsspec1(np)
          nt2 = nsspec2(np)
          ntyp1 = nsptyp1(np)
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
        elseif (lmatch(ni,ntypi,nsspec4(np),nsptyp4(np),.true.)) then
          nt1 = nsspec1(np)
          nt2 = nsspec2(np)
          ntyp1 = nsptyp1(np)
          ntyp2 = nsptyp2(np)
          nt3 = nsspec4(np)
          nt4 = nsspec3(np)
          nt5 = nsspec5(np)
          nt6 = nsspec6(np)
          ntyp3 = nsptyp4(np)
          ntyp4 = nsptyp3(np)
          ntyp5 = nsptyp5(np)
          ntyp6 = nsptyp6(np)
          tr2 = ttr3
          tr3 = ttr2
          tr4 = ttr4
          tr5 = ttr5
        elseif (lmatch(ni,ntypi,nsspec5(np),nsptyp5(np),.true.)) then
          nt1 = nsspec2(np)
          nt2 = nsspec1(np)
          ntyp1 = nsptyp2(np)
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
        elseif (lmatch(ni,ntypi,nsspec6(np),nsptyp6(np),.true.)) then
          nt1 = nsspec2(np)
          nt2 = nsspec1(np)
          ntyp1 = nsptyp2(np)
          ntyp2 = nsptyp1(np)
          nt3 = nsspec6(np)
          nt4 = nsspec5(np)
          nt5 = nsspec3(np)
          nt6 = nsspec4(np)
          ntyp3 = nsptyp6(np)
          ntyp4 = nsptyp5(np)
          ntyp5 = nsptyp3(np)
          ntyp6 = nsptyp4(np)
          tr2 = ttr5
          tr3 = ttr4
          tr4 = ttr2
          tr5 = ttr3
        else
          goto 15
        endif
!
!  Set flags as to whether potentials are symmetric at the two ends
!
        ltsym56 = (lmatch(nt5,ntyp5,nt6,ntyp6,.true.).or.lmatch(nt6,ntyp6,nt5,ntyp5,.true.))
!
!  Molecule handling
!
        if (lmol.and.lneedmol) then
          nmi = natmol(i)
          if (ndim.gt.0) then
            indm = nmolind(i)
            call mindtoijk(indm,ixi,iyi,izi)
          endif
!
!  Set COM coordinates
!
          if (lrigid.and.nmi.gt.0) then
            xcomi = molxyz(1,natinmol(i),nmi)
            ycomi = molxyz(2,natinmol(i),nmi)
            zcomi = molxyz(3,natinmol(i),nmi)
          else
            xcomi = 0.0_dp
            ycomi = 0.0_dp
            zcomi = 0.0_dp
          endif
        else
          xcomi = 0.0_dp
          ycomi = 0.0_dp
          zcomi = 0.0_dp
        endif
!********************************
!  Loop over middle site 2 / j  *
!********************************
        do 25 j = 1,numat
          nj = nat(j)
          ntypj = nftype(j)
          ocj = occuf(j)
          lopj = (lopf(nrelf2a(j)).or..not.lfreeze)
!
!  Check j is allowed for n
!
          if (.not.lmatch(nj,ntypj,nt1,ntyp1,.true.)) goto 25
          if (lmol.and.lneedmol) then
!
!  Molecule handling
!
            nmj = natmol(j)
            if (ndim.gt.0) then
              indmj = nmolind(j)
              call mindtoijk(indmj,ixj,iyj,izj)
              ixi = ixi - ixj
              iyi = iyi - iyj
              izi = izi - izj
            endif
            lmolok = (nmi.eq.nmj.and.nmi.gt.0)
!
!  Set COM coordinates
!
            if (lrigid.and.nmj.gt.0.and.lmolok) then
              xcom(1) = molxyz(1,natinmol(j),nmj) - xcomi
              ycom(1) = molxyz(2,natinmol(j),nmj) - ycomi
              zcom(1) = molxyz(3,natinmol(j),nmj) - zcomi
            else
              xcom(1) = - xcomi
              ycom(1) = - ycomi
              zcom(1) = - zcomi
            endif
          else
            lmolok = .false.
            xcom(1) = - xcomi
            ycom(1) = - ycomi
            zcom(1) = - zcomi
          endif
!
!  Check for intra and but not in same molecule
!
          if (lintra_only.and..not.lmolok) goto 25
          if (lbtyp.and..not.lmolok) goto 25
          x31t = xclat(i) - xclat(j)
          y31t = yclat(i) - yclat(j)
          z31t = zclat(i) - zclat(j)
!
          lmolokij = lmolok
!
!  Check r21 is OK
!  Loop over cell vectors
!
          do 125 ii = 1,iimax
            x31 = x31t + xvec1cell(ii)
            y31 = y31t + yvec1cell(ii)
            z31 = z31t + zvec1cell(ii)
            r312 = x31*x31 + y31*y31 + z31*z31
            if (r312.lt.1d-10) goto 125
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
                if (linter_only) goto 125
                if (lbtyp) then
                  call bonded(lbonded,l2bonds,nbtypeij,nbtypeij2,j,i,0_i4,0_i4,0_i4)
                  if (.not.lbonded) goto 125
                endif
              else
                call lintoijk(ixx,iyy,izz,ii,imaxl,jmaxl,kmaxl)
                if (lbtyp) then
                  call bonded(lbonded,l2bonds,nbtypeij,nbtypeij2,j,i,ixx,iyy,izz)
                  if (.not.lbonded) goto 125
                  lsamemol = (lbonded.or.l2bonds)
                else
                  lsamemol = .false.
                endif
                if (.not.lsamemol) then
                  call samemol(lsamemol,nmi,ixx,iyy,izz,ixi,iyi,izi)
                endif
                if (lintra_only.and..not.lsamemol) goto 125
                if (linter_only.and.lsamemol) goto 125
              endif
            endif
!
!  Distance checking
!
            if (r312.gt.tr2.and.(.not.lbtyp.or..not.lbonded)) goto 125
!
!  j has been accepted
!
            sxyz(1,1) = xclat(j)
            sxyz(2,1) = yclat(j)
            sxyz(3,1) = zclat(j)
!
            sxyz(1,3) = sxyz(1,1) + x31
            sxyz(2,3) = sxyz(2,1) + y31
            sxyz(3,3) = sxyz(3,1) + z31
!**************************************
!  Loop over middle site bonded to j  *
!**************************************
            do 35 k = 1,numat
              nk = nat(k)
              ntypk = nftype(k)
              ock = occuf(k)
              lopk = (lopf(nrelf2a(k)).or..not.lfreeze)
!
!  Check k is allowed for n
!
              if (.not.lmatch(nk,ntypk,nt2,ntyp2,.true.)) goto 35
              if (lmol.and.lneedmol) then
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
                lmolok = (nmj.eq.nmk.and.nmj.gt.0)
!
!  Set COM coordinates
!
                if (lrigid.and.nmk.gt.0.and.lmolok) then
                  xcom(2) = molxyz(1,natinmol(k),nmk) - xcomi
                  ycom(2) = molxyz(2,natinmol(k),nmk) - ycomi
                  zcom(2) = molxyz(3,natinmol(k),nmk) - zcomi
                else
                  xcom(2) = - xcomi
                  ycom(2) = - ycomi
                  zcom(2) = - zcomi
                endif
                xcom(6) = xcom(2) - xcom(1)
                ycom(6) = ycom(2) - ycom(1)
                zcom(6) = zcom(2) - zcom(1)
              else
                lmolok = .false.
                xcom(2) = - xcomi
                ycom(2) = - ycomi
                zcom(2) = - zcomi
                xcom(6) = xcom(2) - xcom(1)
                ycom(6) = ycom(2) - ycom(1)
                zcom(6) = zcom(2) - zcom(1)
              endif
!
!  Check for intra and but not in same molecule
!
              if (lintra_only.and..not.lmolok) goto 35
              if (lbtyp.and..not.lmolok) goto 35
              x21t = xclat(k) - sxyz(1,1)
              y21t = yclat(k) - sxyz(2,1)
              z21t = zclat(k) - sxyz(3,1)
!
              lmolokjk = lmolok
!
!  Check r31 is OK
!  Loop over cell vectors
!
              do 135 jj = 1,iimax
                x21 = x21t + xvec1cell(jj)
                y21 = y21t + yvec1cell(jj)
                z21 = z21t + zvec1cell(jj)
                r212 = x21*x21 + y21*y21 + z21*z21
                if (r212.lt.1d-10) goto 135
!
!  Prevent atoms i and k being the same atom
!
                if (k.eq.i.and.jj.eq.ii) goto 135
!
!  Prevent atoms j and k being the same atom
!
                if (k.eq.j.and.jj.eq.iimid) goto 135
!
!  Reset lmolok to preloop value
!
                lmolok = lmolokjk
!
!  Molecule checking
!
                lbonded = .false.
                if (lmolok) then
                  if (ndim.eq.0) then
                    if (linter_only) goto 135
                    if (lbtyp) then
                      call bonded(lbonded,l2bonds,nbtypejk,nbtypejk2,j,k,0_i4,0_i4,0_i4)
                      if (.not.lbonded) goto 135
                    endif
                  else
                    call lintoijk(jxx,jyy,jzz,jj,imaxl,jmaxl,kmaxl)
                    if (lbtyp) then
                      call bonded(lbonded,l2bonds,nbtypejk,nbtypejk2,j,k,jxx,jyy,jzz)
                      if (.not.lbonded) goto 135
                      lsamemol = (lbonded.or.l2bonds)
                    else
                      lsamemol = .false.
                    endif
                    if (.not.lsamemol) then
                      call samemol(lsamemol,nmj,jxx,jyy,jzz,ixk,iyk,izk)
                    endif
                    if (lintra_only.and..not.lsamemol) goto 135
                    if (linter_only.and.lsamemol) goto 135
                  endif
                endif
!
!  Distance checking
!
                if (r212.gt.tr1.and.(.not.lbtyp.or..not.lbonded)) goto 135
!
!  k has been accepted
!
                sxyz(1,2) = x21 + sxyz(1,1)
                sxyz(2,2) = y21 + sxyz(2,1)
                sxyz(3,2) = z21 + sxyz(3,1)
!************************************
!  Loop over second end site for j  *
!************************************
!
!  Set l looping indices
!
                do 45 l = 1,numat
                  nl = nat(l)
                  ntypl = nftype(l)
                  ocl = occuf(l)
                  lopl = (lopf(nrelf2a(l)).or..not.lfreeze)
!
!  Check l is allowed for n
!
                  if (.not.lmatch(nl,ntypl,nt4,ntyp4,.true.)) goto 45
!
!  Molecule handling
!
                  if (lmol.and.lneedmol) then
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
                      xcom(3) = molxyz(1,natinmol(l),nml) - xcomi
                      ycom(3) = molxyz(2,natinmol(l),nml) - ycomi
                      zcom(3) = molxyz(3,natinmol(l),nml) - zcomi
                    else
                      xcom(3) = - xcomi
                      ycom(3) = - ycomi
                      zcom(3) = - zcomi
                    endif
                    xcom(7) = xcom(3) - xcom(1)
                    ycom(7) = ycom(3) - ycom(1)
                    zcom(7) = zcom(3) - zcom(1)
                    xcom(10) = xcom(3) - xcom(2)
                    ycom(10) = ycom(3) - ycom(2)
                    zcom(10) = zcom(3) - zcom(2)
                  else
                    lmolok = .false.
                    xcom(3) = - xcomi
                    ycom(3) = - ycomi
                    zcom(3) = - zcomi
                    xcom(7) = xcom(3) - xcom(1)
                    ycom(7) = ycom(3) - ycom(1)
                    zcom(7) = zcom(3) - zcom(1)
                    xcom(10) = xcom(3) - xcom(2)
                    ycom(10) = ycom(3) - ycom(2)
                    zcom(10) = zcom(3) - zcom(2)
                  endif
!
!  Check for intra and but not in same molecule
!
                  if (lintra_only.and..not.lmolok) goto 45
                  if (lbtyp.and..not.lmolok) goto 45
                  x41t = xclat(l) - sxyz(1,1)
                  y41t = yclat(l) - sxyz(2,1)
                  z41t = zclat(l) - sxyz(3,1)
!
                  lmolokjl = lmolok
!
!  Check r41 is OK
!  Loop over cell vectors
!
                  do 145 ll = 1,iimax
                    x41 = x41t + xvec1cell(ll)
                    y41 = y41t + yvec1cell(ll)
                    z41 = z41t + zvec1cell(ll)
                    r412 = x41*x41 + y41*y41 + z41*z41
                    if (r412.lt.1d-10) goto 145
!
!  Prevent atoms i and l being the same atom
!
                    if (l.eq.i.and.ll.eq.ii) goto 145
!
!  Prevent atoms j and l being the same atom
!
                    if (l.eq.j.and.ll.eq.iimid) goto 145
!
!  Prevent atoms k and l being the same atom
!
                    if (l.eq.k.and.ll.eq.jj) goto 145
!
!  Reset lmolok to preloop value
!
                    lmolok = lmolokjl
!
!  Molecule checking
!
                    lbonded = .false.
                    if (lmolok) then
                      if (ndim.eq.0) then
                        if (linter_only) goto 145
                        if (lbtyp) then
                          call bonded(lbonded,l2bonds,nbtypejl,nbtypejl2,j,l,0_i4,0_i4,0_i4)
                          if (.not.lbonded) goto 145
                        endif
                      else
                        call lintoijk(kxx,kyy,kzz,ll,imaxl,jmaxl,kmaxl)
                        if (lbtyp) then
                          call bonded(lbonded,l2bonds,nbtypejl,nbtypejl2,j,l,kxx,kyy,kzz)
                          if (.not.lbonded) goto 145
                          lsamemol = (lbonded.or.l2bonds)
                        else
                          lsamemol = .false.
                        endif
                        if (.not.lsamemol) then
                          call samemol(lsamemol,nmj,kxx,kyy,kzz,ixl,iyl,izl)
                        endif
                        if (lintra_only.and..not.lsamemol) goto 145
                        if (linter_only.and.lsamemol) goto 145
                      endif
                    endif
!
!  Distance checking
!
                    if (r412.gt.tr3.and.(.not.lbtyp.or..not.lbonded)) goto 145
!
!  l has been accepted
!
                    sxyz(1,4) = x41 + sxyz(1,1)
                    sxyz(2,4) = y41 + sxyz(2,1)
                    sxyz(3,4) = z41 + sxyz(3,1)
!*****************************************
!  Loop over first end site bonded to k  *
!*****************************************
                    do 55 m = 1,numat
                      nm = nat(m)
                      ntypm = nftype(m)
                      ocm = occuf(m)
                      lopm = (lopf(nrelf2a(m)).or..not.lfreeze)
!
!  Check m is allowed for n
!
                      if (.not.lmatch(nm,ntypm,nt5,ntyp5,.true.)) goto 55
                      if (lmol.and.lneedmol) then
!
!  Molecule handling
!
                        nmm = natmol(m)
                        if (ndim.gt.0) then
                          indmm = nmolind(m)
                          call mindtoijk(indmm,ixm,iym,izm)
                          ixm = ixm - ixk
                          iym = iym - iyk
                          izm = izm - izk
                        endif
                        lmolok = (nmj.eq.nmm.and.nmj.gt.0)
!
!  Set COM coordinates
!
                        if (lrigid.and.nmm.gt.0.and.lmolok) then
                          xcom(4) = molxyz(1,natinmol(m),nmm) - xcomi
                          ycom(4) = molxyz(2,natinmol(m),nmm) - ycomi
                          zcom(4) = molxyz(3,natinmol(m),nmm) - zcomi
                        else
                          xcom(4) = - xcomi
                          ycom(4) = - ycomi
                          zcom(4) = - zcomi
                        endif
                        xcom(8) = xcom(4) - xcom(1)
                        ycom(8) = ycom(4) - ycom(1)
                        zcom(8) = zcom(4) - zcom(1)
                        xcom(11) = xcom(4) - xcom(2)
                        ycom(11) = ycom(4) - ycom(2)
                        zcom(11) = zcom(4) - zcom(2)
                        xcom(13) = xcom(4) - xcom(3)
                        ycom(13) = ycom(4) - ycom(3)
                        zcom(13) = zcom(4) - zcom(3)
                      else
                        lmolok = .false.
                        xcom(4) = - xcomi
                        ycom(4) = - ycomi
                        zcom(4) = - zcomi
                        xcom(8) = xcom(4) - xcom(1)
                        ycom(8) = ycom(4) - ycom(1)
                        zcom(8) = zcom(4) - zcom(1)
                        xcom(11) = xcom(4) - xcom(2)
                        ycom(11) = ycom(4) - ycom(2)
                        zcom(11) = zcom(4) - zcom(2)
                        xcom(13) = xcom(4) - xcom(3)
                        ycom(13) = ycom(4) - ycom(3)
                        zcom(13) = zcom(4) - zcom(3)
                      endif
!
!  Check for intra and but not in same molecule
!
                      if (lintra_only.and..not.lmolok) goto 55
                      if (lbtyp.and..not.lmolok) goto 55
                      x52t = xclat(m) - sxyz(1,2)
                      y52t = yclat(m) - sxyz(2,2)
                      z52t = zclat(m) - sxyz(3,2)
!
                      lmolokkm = lmolok
!
!  Check r52 is OK
!  Loop over cell vectors
!
                      do 155 mm = 1,iimax
                        x52 = x52t + xvec1cell(mm)
                        y52 = y52t + yvec1cell(mm)
                        z52 = z52t + zvec1cell(mm)
                        r522 = x52*x52 + y52*y52 + z52*z52
                        if (r522.lt.1.0d-10) goto 155
!
!  Prevent atoms j and m being the same atom
!
                        if (m.eq.j) then
                          x51 = x52 + x21
                          y51 = y52 + y21
                          z51 = z52 + z21
                          r512 = x51**2 + y51**2 + z51**2
                          if (r512.lt.1.0d-10) goto 155
                        endif
!
!  Prevent atoms k and m being the same atom
!
                        if (m.eq.k.and.mm.eq.iimid) goto 155
!
!  Reset lmolok to preloop value
!
                        lmolok = lmolokkm
!
!  Molecule checking
!
                        lbonded = .false.
                        if (lmolok) then
                          if (ndim.eq.0) then
                            if (linter_only) goto 155
                            if (lbtyp) then
                              call bonded(lbonded,l2bonds,nbtypekm,nbtypekm2,k,m,0_i4,0_i4,0_i4)
                              if (.not.lbonded) goto 155
                            endif
                          else
                            call lintoijk(mxx,myy,mzz,mm,imaxl,jmaxl,kmaxl)
                            if (lbtyp) then
                              call bonded(lbonded,l2bonds,nbtypekm,nbtypekm2,k,m,mxx-jxx,myy-jyy,mzz-jzz)
                              if (.not.lbonded) goto 155
                              lsamemol = (lbonded.or.l2bonds)
                            else
                              lsamemol = .false.
                            endif
                            if (.not.lsamemol) then
                              call samemol(lsamemol,nmk,jxx,jyy,jzz,ixm,iym,izm)
                            endif
                            if (lintra_only.and..not.lsamemol) goto 155
                            if (linter_only.and.lsamemol) goto 155
                          endif
                        endif
!
!  Distance checking
!
                        if (r522.gt.tr4.and.(.not.lbtyp.or..not.lbonded)) goto 155
!
!  m has been accepted
!
                        sxyz(1,5) = x52 + sxyz(1,2)
                        sxyz(2,5) = y52 + sxyz(2,2)
                        sxyz(3,5) = z52 + sxyz(3,2)
!************************************
!  Loop over second end site for k  *
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
                        do 65 n = nmin,nmax
                          nn = nat(n)
                          ntypn = nftype(n)
                          ocn = occuf(n)
                          lopn = (lopf(nrelf2a(n)).or..not.lfreeze)
!
!  Check n is allowed for n
!
                          if (.not.lmatch(nn,ntypn,nt6,ntyp6,.true.)) goto 65
!
!  Molecule handling
!
                          if (lmol.and.lneedmol) then
                            nmn = natmol(n)
                            if (ndim.gt.0) then
                              indmn = nmolind(n)
                              call mindtoijk(indmn,ixn,iyn,izn)
                              ixn = ixn - ixk
                              iyn = iyn - iyk
                              izn = izn - izk
                            endif
                            lmolok = (nmj.eq.nmn.and.nmj.gt.0)
!
!  Set COM coordinates
!
                            if (lrigid.and.nmn.gt.0.and.lmolok) then
                              xcom(5) = molxyz(1,natinmol(n),nmn) - xcomi
                              ycom(5) = molxyz(2,natinmol(n),nmn) - ycomi
                              zcom(5) = molxyz(3,natinmol(n),nmn) - zcomi
                            else
                              xcom(5) = - xcomi
                              ycom(5) = - ycomi
                              zcom(5) = - zcomi
                            endif
                            xcom(9) = xcom(5) - xcom(1)
                            ycom(9) = ycom(5) - ycom(1)
                            zcom(9) = zcom(5) - zcom(1)
                            xcom(12) = xcom(5) - xcom(2)
                            ycom(12) = ycom(5) - ycom(2)
                            zcom(12) = zcom(5) - zcom(2)
                            xcom(14) = xcom(5) - xcom(3)
                            ycom(14) = ycom(5) - ycom(3)
                            zcom(14) = zcom(5) - zcom(3)
                            xcom(15) = xcom(5) - xcom(4)
                            ycom(15) = ycom(5) - ycom(4)
                            zcom(15) = zcom(5) - zcom(4)
                          else
                            lmolok = .false.
                            xcom(5) = - xcomi
                            ycom(5) = - ycomi
                            zcom(5) = - zcomi
                            xcom(9) = xcom(5) - xcom(1)
                            ycom(9) = ycom(5) - ycom(1)
                            zcom(9) = zcom(5) - zcom(1)
                            xcom(12) = xcom(5) - xcom(2)
                            ycom(12) = ycom(5) - ycom(2)
                            zcom(12) = zcom(5) - zcom(2)
                            xcom(14) = xcom(5) - xcom(3)
                            ycom(14) = ycom(5) - ycom(3)
                            zcom(14) = zcom(5) - zcom(3)
                            xcom(15) = xcom(5) - xcom(4)
                            ycom(15) = ycom(5) - ycom(4)
                            zcom(15) = zcom(5) - zcom(4)
                          endif
!
!  Check for intra and but not in same molecule
!
                          if (lintra_only.and..not.lmolok) goto 65
                          if (lbtyp.and..not.lmolok) goto 65
                          x62t = xclat(n) - sxyz(1,2)
                          y62t = yclat(n) - sxyz(2,2)
                          z62t = zclat(n) - sxyz(3,2)
                          if (ltsym56.and.m.eq.n) then
                            niimax = mm - 1
                          else
                            niimax = iimax
                          endif
!
                          lmolokkn = lmolok
!
!  Check r62 is OK
!  Loop over cell vectors
!
                          do 165 mn = 1,niimax
                            x62 = x62t + xvec1cell(mn)
                            y62 = y62t + yvec1cell(mn)
                            z62 = z62t + zvec1cell(mn)
                            r622 = x62*x62 + y62*y62 + z62*z62
                            if (r622.lt.1d-10) goto 165
!
!  Prevent atoms j and n being the same atom
!
                            if (n.eq.j) then
                              x61 = x62 + x21
                              y61 = y62 + y21
                              z61 = z62 + z21
                              r612 = x61**2 + y61**2 + z61**2
                              if (r612.lt.1.0d-10) goto 165
                            endif
!
!  Prevent atoms k and n being the same atom
!
                            if (n.eq.k.and.mn.eq.iimid) goto 165
!
!  Prevent atoms m and n being the same atom
!
                            if (n.eq.m.and.mn.eq.mm) goto 165
!
!  Reset lmolok to preloop value
!
                            lmolok = lmolokkn
!
!  Molecule checking
!
                            lbonded = .false.
                            if (lmolok) then
                              if (ndim.eq.0) then
                                if (linter_only) goto 165
                                if (lbtyp) then
                                  call bonded(lbonded,l2bonds,nbtypekn,nbtypekn2,k,n,0_i4,0_i4,0_i4)
                                  if (.not.lbonded) goto 165
                                endif
                              else
                                call lintoijk(nxx,nyy,nzz,mn,imaxl,jmaxl,kmaxl)
                                if (lbtyp) then
                                  call bonded(lbonded,l2bonds,nbtypekn,nbtypekn2,k,n,nxx-jxx,nyy-jyy,nzz-jzz)
                                  if (.not.lbonded) goto 165
                                  lsamemol = (lbonded.or.l2bonds)
                                else
                                  lsamemol = .false.
                                endif
                                if (.not.lsamemol) then
                                  call samemol(lsamemol,nmi,jxx,jyy,jzz,ixn,iyn,izn)
                                endif
                                if (lintra_only.and..not.lsamemol) goto 165
                                if (linter_only.and.lsamemol) goto 165
                              endif
                            endif
!
!  Distance checking
!
                            if (r622.gt.tr5.and.(.not.lbtyp.or..not.lbonded)) goto 165
!
!  n has been accepted
!
                            sxyz(1,6) = x62 + sxyz(1,2)
                            sxyz(2,6) = y62 + sxyz(2,2)
                            sxyz(3,6) = z62 + sxyz(3,2)
!
!  If lfreeze=.true. and no atoms have any variables then skip this four body term
!
                            if (.not.lopi.and..not.lopj.and..not.lopk.and..not.lopl.and..not.lopm.and..not.lopn) goto 65
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
                                xv6(indvec) = svec(1,jvec,ivec)
                                yv6(indvec) = svec(2,jvec,ivec)
                                zv6(indvec) = svec(3,jvec,ivec)
                                sdist(indvec) = svec(1,jvec,ivec)**2 + svec(2,jvec,ivec)**2 + svec(3,jvec,ivec)**2
                                sdist(indvec) = sqrt(sdist(indvec))
                              enddo
                            enddo
!
                            ofct = oci*ocj*ock*ocl*ocm*ocn
                            rko = rksix*ofct
                            call sixbody(np,nsixtype,sdist,eterm,e1d,e2d,e3d,rko,lgrad1,.false.,.false.)
!********************************************
!  Cross out of plane internal derivatives  *
!********************************************
                            xdrv(iai) = xdrv(iai) + svec(1,1,3)*e1d(2)
                            ydrv(iai) = ydrv(iai) + svec(2,1,3)*e1d(2)
                            zdrv(iai) = zdrv(iai) + svec(3,1,3)*e1d(2)
                            xdrv(iai) = xdrv(iai) + svec(1,2,3)*e1d(6)
                            ydrv(iai) = ydrv(iai) + svec(2,2,3)*e1d(6)
                            zdrv(iai) = zdrv(iai) + svec(3,2,3)*e1d(6)
                            xdrv(iai) = xdrv(iai) + svec(1,4,3)*e1d(10)
                            ydrv(iai) = ydrv(iai) + svec(2,4,3)*e1d(10)
                            zdrv(iai) = zdrv(iai) + svec(3,4,3)*e1d(10)
                            xdrv(iai) = xdrv(iai) + svec(1,5,3)*e1d(11)
                            ydrv(iai) = ydrv(iai) + svec(2,5,3)*e1d(11)
                            zdrv(iai) = zdrv(iai) + svec(3,5,3)*e1d(11)
                            xdrv(iai) = xdrv(iai) + svec(1,6,3)*e1d(12)
                            ydrv(iai) = ydrv(iai) + svec(2,6,3)*e1d(12)
                            zdrv(iai) = zdrv(iai) + svec(3,6,3)*e1d(12)
!
!  End of inner loops over atoms and cell vectors
!
165                         continue
65                       continue
155                     continue
55                   continue
145                 continue
45               continue
135             continue
35           continue
125         continue
25       continue
15     continue
    endif
!
!  End of outer loops
!
  enddo
!
!  Closing banner for energy decomposition
!
  if (lPrintSix) then
    call mpbarrier
    if (ioproc) then
      write(ioout,'(''--------------------------------------------------------------------------------'',/)')
    endif
  endif
!
!  Multiply energy by a factor of a half to allow for double counting in loops
!
  esix = esix + 0.5_dp*esixl
  if (lstr) then
    do kl = 1,nstrains
      rstrd(kl) = rstrd(kl) + 0.5_dp*rstrdl(kl)
    enddo
  endif
!
!  Timing
!
  time2 = g_cpu_time()
  tsix = tsix + time2 - time1
#ifdef TRACE
  call trace_out('sixsd')
#endif
!
  return
  end
