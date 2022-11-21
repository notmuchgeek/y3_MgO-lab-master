  subroutine improper(iout,imode,nphitot)
!
!  Calculates improper torsion angles and prints them out for valid four-body terms
!
!  On entry:
!
!  iout      = channel number for I/O - must be open already
!  imode     = mode number : 0 => standard GULP output
!                            1 => LAMMPS input file output
!
!  On exit:
!
!  nphitot   = total number of torsion terms over all potentials
!
!   7/13 Created from torsion
!   7/13 Number of improper torsions corrected by subtracting number of other torsions
!   2/18 Trace added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
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
  use configurations, only : nregionno, nregiontype, QMMMmode
  use g_constants
  use current
  use element
  use four
  use iochannels
  use molecule
  use spatial,        only : lspatialok, xinbox, yinbox, zinbox
  use symmetry
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4),                     intent(in)    :: iout
  integer(i4),                     intent(in)    :: imode
  integer(i4),                     intent(inout) :: nphitot
!
!  Local variables
!
  character(len=5)                               :: lab1
  character(len=5)                               :: lab2
  character(len=5)                               :: lab3
  character(len=5)                               :: lab4
  character(len=80)                              :: string
  integer(i4)                                    :: i
  integer(i4)                                    :: ii
  integer(i4)                                    :: indmi
  integer(i4)                                    :: indmj
  integer(i4)                                    :: indmk
  integer(i4)                                    :: indml
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
  integer(i4)                                    :: jxx
  integer(i4)                                    :: jyy
  integer(i4)                                    :: jzz
  integer(i4)                                    :: k
  integer(i4)                                    :: kxx
  integer(i4)                                    :: kyy
  integer(i4)                                    :: kzz
  integer(i4)                                    :: l
  integer(i4)                                    :: ll
  integer(i4)                                    :: n
  integer(i4)                                    :: nbtypeji
  integer(i4)                                    :: nbtypejk
  integer(i4)                                    :: nbtypekl
  integer(i4)                                    :: nbtypeji2
  integer(i4)                                    :: nbtypejk2
  integer(i4)                                    :: nbtypekl2
  integer(i4)                                    :: ni
  integer(i4)                                    :: nj
  integer(i4)                                    :: nk
  integer(i4)                                    :: nl 
  integer(i4)                                    :: nmi 
  integer(i4)                                    :: nmj 
  integer(i4)                                    :: nmk 
  integer(i4)                                    :: nml 
  integer(i4)                                    :: nimproper
  integer(i4)                                    :: nphi
  integer(i4)                                    :: nphitotin
  integer(i4)                                    :: nregioni
  integer(i4)                                    :: nregionj
  integer(i4)                                    :: nregionk
  integer(i4)                                    :: nregionl
  integer(i4)                                    :: nregiontypi
  integer(i4)                                    :: nregiontypj
  integer(i4)                                    :: nregiontypk
  integer(i4)                                    :: nregiontypl
  integer(i4)                                    :: nt1
  integer(i4)                                    :: nt2
  integer(i4)                                    :: nt3
  integer(i4)                                    :: nt4
  integer(i4)                                    :: ntyp1
  integer(i4)                                    :: ntyp2
  integer(i4)                                    :: ntyp3
  integer(i4)                                    :: ntyp4
  integer(i4)                                    :: ntypi
  integer(i4)                                    :: ntypj
  integer(i4)                                    :: ntypk
  integer(i4)                                    :: ntypl
  logical                                        :: l2bonds
  logical                                        :: lbonded
  logical                                        :: lbtyp
  logical                                        :: liok
  logical                                        :: limatch1
  logical                                        :: lintra_only
  logical                                        :: linter_only
  logical                                        :: ljmatch2
  logical                                        :: ljkmatch
  logical                                        :: lkmatch3
  logical                                        :: lmatch
  logical                                        :: lmolloc
  logical                                        :: lmolok
  logical                                        :: lneedmol
  logical                                        :: lsamemol
  real(dp)                                       :: c12
  real(dp)                                       :: c34
  real(dp)                                       :: cos2
  real(dp)                                       :: cos3
  real(dp)                                       :: cosphi
  real(dp)                                       :: phi
  real(dp)                                       :: r21
  real(dp)                                       :: r31
  real(dp)                                       :: r32
  real(dp)                                       :: r41
  real(dp)                                       :: r42
  real(dp)                                       :: r43
  real(dp)                                       :: rvp21
  real(dp)                                       :: rvp34
  real(dp)                                       :: tr1
  real(dp)                                       :: tr2
  real(dp)                                       :: tr3
  real(dp)                                       :: tr4
  real(dp)                                       :: vp21x
  real(dp)                                       :: vp21y
  real(dp)                                       :: vp21z
  real(dp)                                       :: vp34x
  real(dp)                                       :: vp34y
  real(dp)                                       :: vp34z
  real(dp)                                       :: x21
  real(dp)                                       :: y21
  real(dp)                                       :: z21
  real(dp)                                       :: x21t
  real(dp)                                       :: y21t
  real(dp)                                       :: z21t
  real(dp)                                       :: x31
  real(dp)                                       :: y31
  real(dp)                                       :: z31
  real(dp)                                       :: x41
  real(dp)                                       :: y41
  real(dp)                                       :: z41
  real(dp)                                       :: x32
  real(dp)                                       :: y32
  real(dp)                                       :: z32
  real(dp)                                       :: x32t
  real(dp)                                       :: y32t
  real(dp)                                       :: z32t
  real(dp)                                       :: x42
  real(dp)                                       :: y42
  real(dp)                                       :: z42
  real(dp)                                       :: x43
  real(dp)                                       :: y43
  real(dp)                                       :: z43
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
#ifdef TRACE
  call trace_in('improper')
#endif
!
  lmolloc = (nmol.gt.0)
!
!  Check mode that has been input
!
  if (imode.lt.0.or.imode.gt.1) then
    call outerror('invalid mode passed to improper subroutine',0_i4)
    call stopnow('torsion')
  endif
!
!  Define error string
!
  string = 'angle in improper torsion has become 180 degrees'
!
!  Check how many four-body potentials are of improper type
!
  nimproper = 0
  do n = 1,nfor
    if (mmfexc(n).eq.2) then
      nimproper = nimproper + 1
    endif
  enddo
!
!  If there are no impropers then return
!
  if (nimproper.eq.0) then
#ifdef TRACE
    call trace_out('improper')
#endif
    return
  endif
!
  nphitotin = nphitot
!
  if (imode.eq.0) then
    write(iout,'(/,''  Analysis of improper four-body terms for primitive cell:'',/)')
    write(iout,'(''--------------------------------------------------------------------------------'')')
    write(iout,'(''  Potential   Atom 1     Atom 2     Atom 3     Atom 4         Phi (degrees)'')')
    write(iout,'(''--------------------------------------------------------------------------------'')')
  endif 
!*************************
!  Loop over potentials  *
!*************************
  pots: do n = 1,nfor
    nphi = 0
    if (mmfexc(n).ne.2) cycle pots
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
!********************************
!  Loop over middle site 2 / j  *
!********************************
    jloop: do j = 1,numat
      nj = nat(j)
      ntypj = nftype(j)
      nregionj = nregionno(nsft+nrelf2a(j))       
      nregiontypj = nregiontype(nregionj,ncf)
!
!  Check whether j is allowed for n
!
      ljmatch2 = lmatch(nj,ntypj,nt2,ntyp2,.true.)
!
!  If no matches is found then cycle
!
      if (.not.ljmatch2) cycle jloop
!
!  j has been accepted
!
      if (lspatialok) then
        xc2 = xinbox(j)
        yc2 = yinbox(j)
        zc2 = zinbox(j)
      else
        xc2 = xclat(j)
        yc2 = yclat(j)
        zc2 = zclat(j)
      endif
!
!  Molecule handling
!
      if (lmolloc.and.lneedmol) then
        nmj = natmol(j)
        if (ndim.gt.0) then
          indmj = nmolind(j)
          call mindtoijk(indmj,ixj,iyj,izj)
        endif
      endif
      call label(nj,ntypj,lab2)
!***************************************
!  Loop over second middle site 3 / k  *
!***************************************
      kloop: do k = 1,numat
        nk = nat(k)
        ntypk = nftype(k)
        nregionk = nregionno(nsft+nrelf2a(k))
        nregiontypk = nregiontype(nregionk,ncf)
!
!  Check whether j and k are allowed for n
!
        lkmatch3 = lmatch(nk,ntypk,nt3,ntyp3,.true.)
!
!  Check whether j-k or k-j orders are OK for 2-3
!
        ljkmatch = (ljmatch2.and.lkmatch3)
!
!  If no pair of matches can be found then cycle
!
        if (.not.ljkmatch) cycle kloop
!
!  QM/MM handling : j & k are both QM atoms and potential is of bonded type => exclude
!
        if (QMMMmode(ncf).gt.0) then
          if (nregiontypj.eq.1.and.nregiontypk.eq.1.and.lbtyp) cycle kloop
        endif
        if (lmolloc.and.lneedmol) then
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
        else
          lmolok = .false.
        endif
!
!  Check for intra and but not in same molecule
!
        if (lintra_only.and..not.lmolok) cycle kloop
        if (lbtyp.and..not.lmolok) cycle kloop
        if (lspatialok) then
          xc3t = xinbox(k)
          yc3t = yinbox(k)
          zc3t = zinbox(k)
        else
          xc3t = xclat(k)
          yc3t = yclat(k)
          zc3t = zclat(k)
        endif
        x32t = xc3t - xc2
        y32t = yc3t - yc2
        z32t = zc3t - zc2
        call label(nk,ntypk,lab3)
!
!  Check r32 is OK
!  Loop over cell vectors
!
        jjloop: do jj = 1,iimax
          r32 = (xvec1cell(jj)+x32t)**2 + (yvec1cell(jj)+y32t)**2 + (zvec1cell(jj)+z32t)**2
          if (r32.lt.1d-12) cycle jjloop
          if (k.eq.j.and.jj.eq.iimid) cycle jjloop
!
!  Molecule checking
!
          lbonded = .false.
          if (lmolok) then
            if (ndim.eq.0) then
              if (linter_only) cycle jjloop
              if (lbtyp) then
                call bonded(lbonded,l2bonds,nbtypejk,nbtypejk2,j,k,0_i4,0_i4,0_i4)
                if (.not.lbonded) cycle jjloop
!
!  Check central bond type for correct order
!
                if (n4botype(1,n).gt.0) then
                  if (n4botype(1,n).ne.nbtypejk) cycle jjloop
                endif
                if (n4botype(2,n).ne.nbtypejk2) cycle jjloop
              endif
            else
              call lintoijk(jxx,jyy,jzz,jj,imaxl,jmaxl,kmaxl)
              if (lbtyp) then
                call bonded(lbonded,l2bonds,nbtypejk,nbtypejk2,j,k,jxx,jyy,jzz)
                if (.not.lbonded) cycle jjloop
!
!  Check central bond type for correct order
!
                if (n4botype(1,n).gt.0) then
                  if (n4botype(1,n).ne.nbtypejk) cycle jjloop
                endif
                if (n4botype(2,n).ne.nbtypejk2) cycle jjloop
                lsamemol = (lbonded.or.l2bonds)
              else
                lsamemol = .false.
              endif
              if (.not.lsamemol) then
                call samemol(lsamemol,nmj,jxx,jyy,jzz,ixk,iyk,izk)
              endif
              if (lintra_only.and..not.lsamemol) cycle jjloop
              if (linter_only.and.lsamemol) cycle jjloop
            endif
          endif
!
!  Distance checking
!
          if (lbtyp) then
            if (.not.lbonded) cycle jjloop
          else
            if (r32.gt.tr2) cycle jjloop
          endif
!
          r32 = sqrt(r32)
          xc3 = xc3t + xvec1cell(jj)
          yc3 = yc3t + yvec1cell(jj)
          zc3 = zc3t + zvec1cell(jj)
          x32 = x32t + xvec1cell(jj)
          y32 = y32t + yvec1cell(jj)
          z32 = z32t + zvec1cell(jj)
!*****************************
!  Loop over end site 1 / i  *
!*****************************
          iloop: do i = 1,numat
            ni = nat(i)
            ntypi = nftype(i)
            nregioni = nregionno(nsft+nrelf2a(i))
            nregiontypi = nregiontype(nregioni,ncf)
!
!  Check whether i matches either of types 1 and 4
!
            limatch1 = lmatch(ni,ntypi,nt1,ntyp1,.true.)
!
!  Is i allowed for type 1?
!
            liok = limatch1
            if (.not.liok) cycle iloop
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
            if (lintra_only.and..not.lmolok) cycle iloop
            if (lbtyp.and..not.lmolok) cycle iloop
!
            if (lspatialok) then
              xc1t = xinbox(i)
              yc1t = yinbox(i)
              zc1t = zinbox(i)
            else
              xc1t = xclat(i)
              yc1t = yclat(i)
              zc1t = zclat(i)
            endif
            x21t = xc2 - xc1t
            y21t = yc2 - yc1t
            z21t = zc2 - zc1t
            call label(ni,ntypi,lab1)
!
!  Check r21 is OK
!  Loop over cell vectors
!
            iiloop: do ii = 1,iimax
              r21 = (-xvec1cell(ii) + x21t)**2 + (-yvec1cell(ii) + y21t)**2 + (-zvec1cell(ii) + z21t)**2
              if (r21.lt.1d-12) cycle iiloop
!
!  Prevent atoms i and k being the same atom
!
              if (k.eq.i.and.ii.eq.jj) cycle iiloop
!
!  Molecule checking
!
              lbonded = .false.
              if (lmolok) then
                if (ndim.eq.0) then
                  if (linter_only) cycle iiloop
                  if (lbtyp) then
                    call bonded(lbonded,l2bonds,nbtypeji,nbtypeji2,j,i,0_i4,0_i4,0_i4)
                    if (.not.lbonded) cycle iiloop
                  endif
                else
                  call lintoijk(ixx,iyy,izz,ii,imaxl,jmaxl,kmaxl)
                  if (lbtyp) then
                    call bonded(lbonded,l2bonds,nbtypeji,nbtypeji2,j,i,ixx,iyy,izz)
                    if (.not.lbonded) cycle iiloop
                    lsamemol = (lbonded.or.l2bonds)
                  else
                    lsamemol = .false.
                  endif
                  if (.not.lsamemol) then
                    call samemol(lsamemol,nmj,ixx,iyy,izz,ixi,iyi,izi)
                  endif
                  if (lintra_only.and..not.lsamemol) cycle iiloop
                  if (linter_only.and.lsamemol) cycle iiloop
                endif
              endif
!
!  Distance checking
!
              if (lbtyp) then
                if (.not.lbonded) cycle iiloop
              else
                if (r21.gt.tr1) cycle iiloop
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
              r31 = x31*x31 + y31*y31 + z31*z31
              if (r31.lt.1.0d-8) cycle iiloop
!
              r21 = sqrt(r21)
              r31 = sqrt(r31)
!**********************************
!  Loop over last end site 4 / l  *
!**********************************
              lloop: do l = 1,numat
                nl = nat(l)
                ntypl = nftype(l)
                nregionl = nregionno(nsft+nrelf2a(l))
                nregiontypl = nregiontype(nregionl,ncf)
!
!  Check l is allowed for n
!
                if (.not.lmatch(nl,ntypl,nt4,ntyp4,.true.)) cycle lloop
!
!  QM/MM handling : i, j, k & l are all QM atoms => exclude
!               
                if (QMMMmode(ncf).gt.0) then
                  if (nregiontypi.eq.1.and.nregiontypj.eq.1.and.nregiontypk.eq.1.and.nregiontypl.eq.1) cycle lloop
                endif
!
                if (lmolloc.and.lneedmol) then
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
                if (lintra_only.and..not.lmolok) cycle lloop
                if (lbtyp.and..not.lmolok) cycle lloop
                if (lspatialok) then
                  xc4t = xinbox(l)
                  yc4t = yinbox(l)
                  zc4t = zinbox(l)
                else
                  xc4t = xclat(l)
                  yc4t = yclat(l)
                  zc4t = zclat(l)
                endif
                x43t = xc4t - xc3
                y43t = yc4t - yc3
                z43t = zc4t - zc3
                call label(nl,ntypl,lab4)
!
!  Check r43 is OK
!  Loop over cell vectors
!
                llloop: do ll = 1,iimax
                  r43 = (xvec1cell(ll)+x43t)**2 + (yvec1cell(ll)+y43t)**2 + (zvec1cell(ll)+z43t)**2
                  if (r43.lt.1d-12) cycle llloop
!
!  Molecule checking
!
                  lbonded = .false.
                  if (lmolok) then
                    if (ndim.eq.0) then
                      if (linter_only) cycle llloop
                      if (lbtyp) then
                        call bonded(lbonded,l2bonds,nbtypekl,nbtypekl2,j,l,0_i4,0_i4,0_i4)
                        if (.not.lbonded) cycle llloop
                      endif
                    else
                      call lintoijk(kxx,kyy,kzz,ll,imaxl,jmaxl,kmaxl)
                      if (lbtyp) then
                        call bonded(lbonded,l2bonds,nbtypekl,nbtypekl2,j,l,kxx,kyy,kzz)
                        if (.not.lbonded) cycle llloop
                        lsamemol = (lbonded.or.l2bonds)
                      else
                        lsamemol = .false.
                      endif
                      if (.not.lsamemol) then
                        call samemol(lsamemol,nmj,kxx,kyy,kzz,ixl,iyl,izl)
                      endif
                      if (lintra_only.and..not.lsamemol) cycle llloop
                      if (linter_only.and.lsamemol) cycle llloop
                    endif
                  endif
!
!  Distance checking
!
                  if (lbtyp) then
                    if (.not.lbonded) cycle llloop
                  else
                    if (r43.gt.tr3) cycle llloop
                  endif
!
                  r43 = sqrt(r43)
                  x43 = x43t + xvec1cell(ll)
                  y43 = y43t + yvec1cell(ll)
                  z43 = z43t + zvec1cell(ll)
!
!  Check r41 is OK
!
                  x41 = x43 + x32 + x21
                  y41 = y43 + y32 + y21
                  z41 = z43 + z32 + z21
                  r41 = x41*x41 + y41*y41 + z41*z41
                  if (r41.gt.tr4.and.tr4.gt.0.0_dp.and..not.lbtyp) cycle llloop
                  if (r41.lt.1.0d-8) cycle llloop
!
!  Check r42 is OK
!
                  x42 = x32 + x43
                  y42 = y32 + y43
                  z42 = z32 + z43
                  r42 = x42*x42 + y42*y42 + z42*z42
                  if (r42.lt.1.0d-8) cycle llloop
!*********************************
!  Valid four-body term located  *
!*********************************
!                   
!  Generate perpendicular end vectors
!                   
                  c12 = (x21*x32+y21*y32+z21*z32)
                  cos2 = - c12/(r21*r32)
                  c12 = c12/(r32*r32)
                  if (cos2.le.-0.99999999) then
                    call outerror(string,0_i4)
                    call stopnow('torsion')
                  endif
                  if (cos2.ge.0.99999999) then
                    call outerror(string,0_i4)
                    call stopnow('torsion')
                  endif
                  vp21x = c12*x32 - x21
                  vp21y = c12*y32 - y21
                  vp21z = c12*z32 - z21
                  rvp21 = (vp21x*vp21x) + (vp21y*vp21y) + (vp21z*vp21z)
                  rvp21 = sqrt(rvp21)
!
                  c34 = (x43*x32+y43*y32+z43*z32)
                  cos3 = - c34/(r43*r32)
                  c34 = c34/(r32*r32)
                  if (cos3.le.-0.99999999) then 
                    call outerror(string,0_i4)
                    call stopnow('torsion')
                  endif
                  if (cos3.ge.0.99999999) then
                    call outerror(string,0_i4)
                    call stopnow('torsion')
                  endif
                  vp34x = x43 - c34*x32
                  vp34y = y43 - c34*y32
                  vp34z = z43 - c34*z32
                  rvp34 = (vp34x*vp34x) + (vp34y*vp34y) + (vp34z*vp34z)
                  rvp34 = sqrt(rvp34)
!
                  cosphi = (vp21x*vp34x+vp21y*vp34y+vp21z*vp34z)
                  cosphi = cosphi/(rvp21*rvp34)
                  if (cosphi.gt.1.0_dp) cosphi = 1.0_dp
                  if (cosphi.lt.-1.0_dp) cosphi = - 1.0_dp
                  phi = acos(cosphi)*radtodeg
!
                  nphi = nphi + 1
                  if (imode.eq.0) then
!
!  GULP standard I/O
!
                    if (nphi.eq.1) then
                      write(iout,'(5x,i3,4x,4(1x,a5,1x,i4),7x,f10.6)') n,lab1,i,lab2,j,lab3,k,lab4,l,phi
                    else
                      write(iout,'(12x,4(1x,a5,1x,i4),7x,f10.6)') lab1,i,lab2,j,lab3,k,lab4,l,phi
                    endif
                  elseif (imode.eq.1) then
!
!  LAMMPS file I/O
!
                    write(iout,'(6i8)') nphi,n,i,j,k,l
                  endif
!
!  End of inner loops over atoms and cell vectors
!
                enddo llloop
              enddo lloop
            enddo iiloop
          enddo iloop
        enddo jjloop
      enddo kloop
    enddo jloop
    if (imode.eq.0) then
      if (nphi.eq.0) then
        write(iout,'(5x,i3,9x,''No angles found'')') n
      endif
    endif
    nphitot = nphitot + nphi
    if (imode.eq.0) then
      write(iout,'(''--------------------------------------------------------------------------------'')')
    endif
!
!  End of outer loop over potentials
!
  enddo pots
  if (imode.eq.0) then
    if (nphitot.gt.nphitotin) then
      write(iout,'(''  Total number of improper torsions = '',i8)') nphitot - nphitotin
      write(iout,'(''--------------------------------------------------------------------------------'')')
    endif
    write(iout,'(/)')
  endif
#ifdef TRACE
  call trace_out('improper')
#endif
!
  return
  end
