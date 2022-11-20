  subroutine molind
!
!  Determines cell indices for molecules
!
!   3/95 Determination of molecule dimensionality added
!   2/01 lintoijk calls now use imaxl,jmaxl and kmaxl
!   6/01 Call for checking for wrapping of coordinates added
!   5/04 Modified to use inbox coordinates if lspatialok is true
!   8/06 Pointers nullified on first call
!   1/07 Modified for noautobond option
!   4/07 Searching for bonded atoms now uses nbonds number rather
!        than checking for a non-zero type
!   4/07 Handling of case where connectivity is input for several 
!        images of same atom by checking when a connection has
!        already been used.
!   5/07 Call to connectwrap excluded if .not.lmodco
!   6/07 Algorithm to find atom nearest to origin changed to be 
!        more robust for large distances.
!   6/09 Modified to use new atom in molecule pointers
!   5/10 Timing added
!   3/14 Pointers nullified on declaration line
!   2/18 Trace added
!   2/20 Use of inbox coordinates for spatial removed
!  12/20 Trap for case where no searching for valid cell indices is
!        required added to speed up code for large systems
!  12/20 Use of molecule specific pointer added to connection lists
!   1/21 If condition on molecule number not being zero added before
!        testing molecule connectivity
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
!  Copyright Curtin University 2021
!
!  Julian Gale, CIC, Curtin University, January 2021
!
  use control,        only : lmodco
  use current
  use element
  use iochannels
  use molecule
  use parallel
  use reallocate
  use species
  use times,          only : tmol
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ic
  integer(i4)                                  :: icm
  integer(i4)                                  :: ierror
  integer(i4)                                  :: ii
  integer(i4)                                  :: iimin
  integer(i4)                                  :: iii
  integer(i4), dimension(:), pointer,     save :: iimptr => null()
  integer(i4)                                  :: iio
  integer(i4)                                  :: ind
  integer(i4)                                  :: ix
  integer(i4)                                  :: iy
  integer(i4)                                  :: iz
  integer(i4)                                  :: ixx
  integer(i4)                                  :: iyy
  integer(i4)                                  :: izz
  integer(i4)                                  :: j
  integer(i4)                                  :: jj
  integer(i4)                                  :: jjo
  integer(i4)                                  :: k
  integer(i4)                                  :: kk
  integer(i4)                                  :: kko
  integer(i4)                                  :: na
  integer(i4)                                  :: nact
  integer(i4)                                  :: nassign
  integer(i4)                                  :: nasslast
  integer(i4)                                  :: nbon
  integer(i4)                                  :: nbond
  integer(i4)                                  :: nc1
  integer(i4)                                  :: nc2
  integer(i4)                                  :: ni
  integer(i4)                                  :: nj
  integer(i4)                                  :: nj1
  integer(i4)                                  :: nk
  integer(i4)                                  :: nmi
  integer(i4)                                  :: nmlind
  integer(i4)                                  :: nmlindo
  integer(i4)                                  :: nno
  integer(i4)                                  :: npn
  integer(i4)                                  :: nptr
  integer(i4)                                  :: npy
  integer(i4), dimension(:), allocatable       :: npno
  integer(i4), dimension(:), allocatable       :: npyes
  integer(i4)                                  :: nyes
  integer(i4)                                  :: status
  logical                                      :: lassigned
  logical                                      :: lfound
  logical                                      :: lskip
  logical                                      :: lxcross
  logical                                      :: lycross
  logical                                      :: lzcross
  logical,     dimension(:), allocatable       :: lAlreadyUsed
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: r2
  real(dp)                                     :: r2min
  real(dp)                                     :: rcut
  real(dp)                                     :: ri
  real(dp)                                     :: rij
  real(dp)                                     :: rj
  real(dp)                                     :: rjk
  real(dp)                                     :: rk
  real(dp)                                     :: rmin
  real(dp)                                     :: t1
  real(dp)                                     :: t2
  real(dp)                                     :: xal
  real(dp)                                     :: yal
  real(dp)                                     :: zal
  real(dp)                                     :: xcd
  real(dp)                                     :: ycd
  real(dp)                                     :: zcd
  real(dp)                                     :: xcdi
  real(dp)                                     :: ycdi
  real(dp)                                     :: zcdi
  real(dp)                                     :: xcrd
  real(dp)                                     :: ycrd
  real(dp)                                     :: zcrd
#ifdef TRACE
  call trace_in('molind')
#endif
!
  t1 = g_cpu_time()
!
!  Allocate local memory
!
  call realloc(iimptr,maxbond,ierror)
  if (ierror.ne.0) call outofmemory('molind','iimptr')
  allocate(npno(numat),stat=status)
  if (status/=0) call outofmemory('molind','npno')
  allocate(npyes(numat),stat=status)
  if (status/=0) call outofmemory('molind','npyes')
  if (nconnect.gt.0) then
    allocate(lAlreadyUsed(nconnect),stat=status)
    if (status/=0) call outofmemory('molind','lAlreadyUsed')
  endif
!
!  Handle periodic wrapping of connection indices
!
  if (nconnect.gt.0.and.lmodco) call connectwrap
!
!  Note - higher tolerance can be used here than in initial search
!  for bonds, as the atoms in the molecule have been identified
!  already and this allows for greater distortions during energy
!  minimisation.
!
  do i = 1,numat
    nmolind(i) = 0
  enddo
!**************************************
!  Obtain cell indices for each atom  *
!**************************************
  do i = 1,nmol
!
!  lxcross etc are logicals indicating whether a boundary
!  crossing has occur so far for this molecule.
!
    lxcross = .false.
    lycross = .false.
    lzcross = .false.
    if (nmolatom(i).gt.0) then
!
!  Check whether connections have been specified across cell boundaries between
!  the same atom - this implies periodicity.
!
      if (nmolconnect(i).gt.0) then
        do icm = 1,nmolconnect(i)
          ic = nmolconnectptr(icm,i)
          do j = 1,nmolatom(i)
            if (n1connect(ic).eq.j) then
              if (n2connect(ic).eq.j) then
                call mindtoijk(nconnectind(ic),ii,jj,kk)
                if (ii.ne.0) lxcross = .true.
                if (jj.ne.0) lycross = .true.
                if (kk.ne.0) lzcross = .true.
              endif
            endif
          enddo
        enddo
      endif
!
!  For no auto-bond case with mod of coordinates then assume user has provided the correct input
!
      if (lnoautobond.and..not.lmodco.and.nconnect.gt.0) then
        nc1 = 0
        nc2 = 0
        do ic = 1,nconnect
          if (nconnectcfg(ic).eq.ncf) then
            nc1 = nc1 + 1
            if (nconnectind(ic).eq.555) then
              nc2 = nc2 + 1
            endif
          endif
        enddo
        lskip = (nc1.eq.nc2)
        if (lskip) then
          do j = 1,nmolatom(i)
            nmolind(nmollist(nmolptr(i)+j)) = 555
          enddo 
        endif
      else
        lskip = .false.
      endif
!
      if (.not.lskip) then
!
!  Locate most central atom of molecule
!
        do j = 1,nmolatom(i)
          na = nmollist(nmolptr(i)+j)
          if (ndim.eq.3) then
            rij = (xfrac(na)-0.5)**2 + (yfrac(na)-0.5)**2 + (zfrac(na)-0.5)**2
          elseif (ndim.eq.2) then
            rij = (xfrac(na)-0.5)**2 + (yfrac(na)-0.5)**2
          elseif (ndim.eq.1) then
            rij = (xfrac(na)-0.5)**2
          endif
          if (j.eq.1) then
            nptr = 1
            rmin = rij
          elseif (rij.lt.rmin) then
            nptr = j
            rmin = rij
          endif
        enddo
!
!  Recursive search to locate cell index for each atom of molecule
!
        nact = nmollist(nmolptr(i)+nptr)
        nmolind(nact) = 555
        nyes = 1
        nassign = 1
        npyes(1) = nptr
        nno = 0
        do j = 1,nmolatom(i)
          if (j.ne.nptr) then
            nno = nno + 1
            npno(nno) = j
          endif
        enddo
 10     nasslast = nassign
        do j = 1,nno
          npn = nmollist(nmolptr(i)+npno(j))
          xal = xclat(npn)
          yal = yclat(npn)
          zal = zclat(npn)
          nj = nat(npn)
          if (nj.gt.maxele) nj = nj - maxele
          rj = rcov(nj)
          lassigned = .false.
          nmlindo = 0
          do k = 1,nyes
            npy = nmollist(nmolptr(i)+npyes(k))
            nk = nat(npy)
            if (nk.gt.maxele) nk = nk - maxele
            rk = rcov(nk)
            rcut = rtol*(rj+rk)
            rcut = rcut*rcut
            ind = nmolind(npy)
            iz = (ind/100) - 5
            ind = ind - 100*(iz+5)
            iy = (ind/10) - 5
            ind = ind - 10*(iy+5)
            ix = ind - 5
!
!  Check for input-driven connectivity first
!
            if (nmolconnect(i).gt.0) then
              do icm = 1,nmolconnect(i)
                ic = nmolconnectptr(icm,i)
                if (n1connect(ic).eq.npy.and.n2connect(ic).eq.npn) then
                  if (nconnectind(ic).gt.0) then
                    call mindtoijk(nconnectind(ic),ii,jj,kk)
                  else
!
!  Find nearest image
!
                    xcdi = xal - xclat(npy)
                    ycdi = yal - yclat(npy)
                    zcdi = zal - zclat(npy)
                    r2min = 1000000.0_dp
                    iimin = 0
                    do iii = 1,iimax
                      xcrd = xcdi + xvec1cell(iii)
                      ycrd = ycdi + yvec1cell(iii)
                      zcrd = zcdi + zvec1cell(iii)
                      r2 = xcrd*xcrd + ycrd*ycrd + zcrd*zcrd
                      if (r2.lt.r2min) then
                        r2min = r2
                        iimin = iii
                      endif
                    enddo
                    call lintoijk(ii,jj,kk,iimin,imaxl,jmaxl,kmaxl)
                  endif
                  nmlind = (ii+ix+5)+10*(jj+iy+5)+100*(kk+iz+5)
                  if (lassigned.and.nmlind.ne.nmlindo) then
                    if (iio.ne.(ii+ix)) lxcross = .true.
                    if (jjo.ne.(jj+iy)) lycross = .true.
                    if (kko.ne.(kk+iz)) lzcross = .true.
                  elseif (.not.lassigned) then
                    nmolind(npn) = nmlind
                    nassign = nassign + 1
                    lassigned = .true.
                    nmlindo = nmlind
                    iio = ii + ix
                    jjo = jj + iy
                    kko = kk + iz
                  endif
                elseif (n1connect(ic).eq.npn.and.n2connect(ic).eq.npy) then
                  if (nconnectind(ic).gt.0) then
                    call mindtoijk(1110_i4-nconnectind(ic),ii,jj,kk)
                  else
!
!  Find nearest image
!
                    xcdi = xal - xclat(npy)
                    ycdi = yal - yclat(npy)
                    zcdi = zal - zclat(npy)
                    r2min = 1000000.0_dp
                    iimin = 0
                    do iii = 1,iimax
                      xcrd = xcdi + xvec1cell(iii)
                      ycrd = ycdi + yvec1cell(iii)
                      zcrd = zcdi + zvec1cell(iii)
                      r2 = xcrd*xcrd + ycrd*ycrd + zcrd*zcrd
                      if (r2.lt.r2min) then
                        r2min = r2
                        iimin = iii
                      endif
                    enddo
                    call lintoijk(ii,jj,kk,iimin,imaxl,jmaxl,kmaxl)
                  endif
                  nmlind = (ii+ix+5) + 10*(jj+iy+5) + 100*(kk+iz+5)
                  if (lassigned.and.nmlind.ne.nmlindo) then
                    if (iio.ne.(ii+ix)) lxcross = .true.
                    if (jjo.ne.(jj+iy)) lycross = .true.
                    if (kko.ne.(kk+iz)) lzcross = .true.
                  elseif (.not.lassigned) then
                    nmolind(npn) = nmlind
                    nassign = nassign + 1
                    lassigned = .true.
                    nmlindo = nmlind
                    iio = ii + ix
                    jjo = jj + iy
                    kko = kk + iz
                  endif
                endif
              enddo
            endif
            if (.not.lnoautobond) then
!
!  Now use interatomic distance checks
!
              xcdi = xal - xclat(npy)
              ycdi = yal - yclat(npy)
              zcdi = zal - zclat(npy)
              do iii = 1,iimax
                xcrd = xcdi + xvec1cell(iii)
                ycrd = ycdi + yvec1cell(iii)
                zcrd = zcdi + zvec1cell(iii)
                rjk = xcrd*xcrd + ycrd*ycrd + zcrd*zcrd
                if (rjk.lt.rcut.and.rjk.gt.0.0_dp) then
                  call lintoijk(ii,jj,kk,iii,imaxl,jmaxl,kmaxl)
                  nmlind = (ii+ix+5) + 10*(jj+iy+5) + 100*(kk+iz+5)
                  if (lassigned.and.nmlind.ne.nmlindo) then
                    if (iio.ne.(ii+ix)) lxcross = .true.
                    if (jjo.ne.(jj+iy)) lycross = .true.
                    if (kko.ne.(kk+iz)) lzcross = .true.
                  elseif (.not.lassigned) then
                    nmolind(npn) = nmlind
                    nassign = nassign + 1
                    lassigned = .true.
                    nmlindo = nmlind
                    iio = ii + ix
                    jjo = jj + iy
                    kko = kk + iz
                  endif
                  if (lassigned.and.lxcross.and.lycross.and.lzcross) goto 20
                endif
              enddo
            endif
          enddo
 20       continue
!
!  End loop over unindexed sites
!
        enddo
        if (nmolatom(i).eq.1.and..not.lnoautobond) then
!
!  Check to see whether unassigned atoms are bonded to periodic replications of themselves
!
          npn = nmollist(nmolptr(i)+1)
          xal = xclat(npn)
          yal = yclat(npn)
          zal = zclat(npn)
          nj = nat(npn)
          if (nj.gt.maxele) nj = nj - maxele
          rj = rcov(nj)
          lassigned = .false.
          rcut = 2.0*rtol*rj
          rcut = rcut*rcut
          xcdi = 0.0_dp
          ycdi = 0.0_dp
          zcdi = 0.0_dp
          do iii = 1,iimax
            xcrd = xcdi + xvec1cell(iii)
            ycrd = ycdi + yvec1cell(iii)
            zcrd = zcdi + zvec1cell(iii)
            rjk = xcrd*xcrd + ycrd*ycrd + zcrd*zcrd
            if (rjk.lt.rcut.and.rjk.gt.0.0_dp.and.iii.ne.iimid) then
              call lintoijk(ii,jj,kk,iii,imaxl,jmaxl,kmaxl)
              if (ii.ne.0) lxcross = .true.
              if (jj.ne.0) lycross = .true.
              if (kk.ne.0) lzcross = .true.
              if (lxcross.and.lycross.and.lzcross) goto 30
            endif
          enddo
 30       continue
!
!  Check to see that more assignments have been made - if not
!  this means that the a bond has increased beyond the tolerance.
!
        else
!
!  Set up
!
          nno = 0
          nyes = 0
          do j = 1,nmolatom(i)
            if (nmolind(nmollist(nmolptr(i)+j)).eq.0) then
              nno = nno + 1
              npno(nno) = j
            else
              nyes = nyes + 1
              npyes(nyes) = j
            endif
          enddo
          if (nassign.eq.nasslast) then
            call outerror('Bond length in molecule has exceeded tolerance',0_i4)
            if (ioproc) then
              write(ioout,'(''!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'')')
              write(ioout,'(''!!  Configuration number = '',i4)')ncf
              write(ioout,'(''!!  Molecule      number = '',i4)')i
              write(ioout,'(''!!  Atoms in molecule with missing bonds:'')')
              write(ioout,'(''!!  '',10(2x,i4))')(nmollist(nmolptr(i)+npno(j)),j=1,nno)
              write(ioout,'(''!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'',/)')
              call outed
            endif
            call stopnow('molind')
          endif
        endif
!
!  Latest round of assignments made
!  If there are any unassigned sites repeat loop
!
        if (nassign.ne.nmolatom(i).and.nassign.ne.nasslast) goto 10
      endif
!
!  End of skip for special case
!
    endif
!
!  Assign dimensionality
!
    if (lxcross) then
      if (lycross) then
        if (lzcross) then
          moldim(i) = 3
        else
          moldim(i) = 2
          moldimi(i) = 1
        endif
      elseif (lzcross) then
        moldim(i) = 2
        moldimi(i) = 2
      else
        moldim(i) = 1
        moldimi(i) = 1
      endif
    elseif (lycross) then
      if (lzcross) then
        moldim(i) = 2
        moldimi(i) = 3
      else
        moldim(i) = 1
        moldimi(i) = 2
      endif
    elseif (lzcross) then
      moldim(i) = 1
      moldimi(i) = 3
    endif
!
!  End of loop over molecules
!
  enddo
!******************************
!  Find cell index for bonds  *
!******************************
  do i = 1,numat
    xal = xclat(i)
    yal = yclat(i)
    zal = zclat(i)
    ni = nat(i)
    if (ni.gt.maxele) ni = ni - maxele
    ri = rcov(ni)
    nmi = natmol(i)
    nbond = 1
    if (nmi.gt.0) then
      if (nmolconnect(nmi).gt.0) then
        j = nbonded(nbond,i)
        lAlreadyUsed(1:nmolconnect(nmi)) = .false.
        do while (j.ne.0.and.nbond.le.nbonds(i))
          lfound = .false.
          icm = 0
          do while (icm.lt.nmolconnect(nmi).and..not.lfound)
            icm = icm + 1
            ic = nmolconnectptr(icm,nmi)
            if (nconnectcfg(ic).eq.ncf.and..not.lAlreadyUsed(icm)) then
              if (n1connect(ic).eq.i.and.n2connect(ic).eq.j) then
                lfound = .true.
                lAlreadyUsed(icm) = .true.
                if (nconnectind(ic).gt.0) then
                  nbondind(nbond,i) = nconnectind(ic)
                else
!
!  Find nearest image
!
                  xcdi = xclat(j) - xal
                  ycdi = yclat(j) - yal
                  zcdi = zclat(j) - zal
                  r2min = 1000000.0_dp
                  iimin = 0
                  do iii = 1,iimax
                    xcrd = xcdi + xvec1cell(iii)
                    ycrd = ycdi + yvec1cell(iii)
                    zcrd = zcdi + zvec1cell(iii)
                    r2 = xcrd*xcrd + ycrd*ycrd + zcrd*zcrd
                    if (r2.lt.r2min) then
                      r2min = r2
                      iimin = iii
                    endif
                  enddo
                  call lintoijk(ii,jj,kk,iimin,imaxl,jmaxl,kmaxl)
                  nbondind(nbond,i) = (ii+5) + 10*(jj+5) + 100*(kk+5)
                endif
                iimptr(nbond) = 0
              elseif (n2connect(ic).eq.i.and.n1connect(ic).eq.j) then
                lfound = .true.
                lAlreadyUsed(icm) = .true.
                if (nconnectind(ic).gt.0) then
                  nbondind(nbond,i) = 1110 - nconnectind(ic)
                else
!
!  Find nearest image
!
                  xcdi = xclat(j) - xal
                  ycdi = yclat(j) - yal
                  zcdi = zclat(j) - zal
                  r2min = 1000000.0_dp
                  iimin = 0
                  do iii = 1,iimax
                    xcrd = xcdi + xvec1cell(iii)
                    ycrd = ycdi + yvec1cell(iii)
                    zcrd = zcdi + zvec1cell(iii)
                    r2 = xcrd*xcrd + ycrd*ycrd + zcrd*zcrd
                    if (r2.lt.r2min) then
                      r2min = r2
                      iimin = iii
                    endif
                  enddo
                  call lintoijk(ii,jj,kk,iimin,imaxl,jmaxl,kmaxl)
                  nbondind(nbond,i) = (ii+5)+10*(jj+5)+100*(kk+5)
                endif
                iimptr(nbond) = 0
              endif
            endif
          enddo
          nbond = nbond + 1
          if (nbond.gt.maxbond) then
            maxbond = nbond + 2
            call changemaxbond
            call realloc(iimptr,maxbond,ierror)
            if (ierror.ne.0) call outofmemory('molind','iimptr')
          endif
          j = nbonded(nbond,i)
        enddo
      endif
    endif
    if (.not.lnoautobond.and.ri.ne.0.0_dp) then
!
!  Loop over bonds
!
      j = nbonded(nbond,i)
      do while (j.gt.0.and.nbond.le.nbonds(i))
        xcd = xclat(j)
        ycd = yclat(j)
        zcd = zclat(j)
        nj1 = nat(j)
        nj = nj1
        if (nj.gt.maxele) nj = nj - maxele
        rj = rcov(nj)
        rcut = rtol*(ri+rj)
        rcut = rcut*rcut
        if (rj.ne.0.0_dp) then
          xcdi = xcd - xal
          ycdi = ycd - yal
          zcdi = zcd - zal
!
!  Set starting point in cell vector loop to avoid same index
!  as previous occurance of the same atom
!
          iimin = 1
          lfound = .false.
          nbon = nbond - 1
          do while (.not.lfound.and.nbon.gt.0) 
            if (nbonded(nbon,i).eq.j) then
              iimin = iimptr(nbon) + 1
              lfound = .true.
            endif
            nbon = nbon - 1
          enddo
          lfound = .false.
          ii = iimin
          do while (ii.le.iimax.and..not.lfound)
            xcrd = xcdi + xvec1cell(ii)
            ycrd = ycdi + yvec1cell(ii)
            zcrd = zcdi + zvec1cell(ii)
            rij = xcrd*xcrd + ycrd*ycrd + zcrd*zcrd
            if (rij.le.rcut.and.rij.gt.0.0_dp.and.(i.ne.j.or.ii.ne.iimid)) then
!
!  Valid bond
!
              lfound = .true.
              iimptr(nbond) = ii
              if (ii.eq.iimid) then
                nbondind(nbond,i) = 555
              else
                call lintoijk(ixx,iyy,izz,ii,imaxl,jmaxl,kmaxl)
                ixx = ixx + 5
                iyy = iyy + 5
                izz = izz + 5
                nbondind(nbond,i) = ixx + 10*iyy + 100*izz
              endif
            endif
            ii = ii + 1
          enddo
        endif
        nbond = nbond + 1
        j = nbonded(nbond,i)
      enddo
    endif
  enddo
!
!  Free local memory
!
  if (nconnect.gt.0) then
    deallocate(lAlreadyUsed,stat=status)
    if (status/=0) call deallocate_error('molind','lAlreadyUsed')
  endif
  deallocate(npyes,stat=status)
  if (status/=0) call deallocate_error('molind','npyes')
  deallocate(npno,stat=status)
  if (status/=0) call deallocate_error('molind','npno')
  call realloc(iimptr,0_i4,ierror)
!
  t2 = g_cpu_time()
  tmol = tmol + t2 - t1
#ifdef TRACE
  call trace_out('molind')
#endif
!
  return
  end
