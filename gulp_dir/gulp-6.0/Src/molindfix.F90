  subroutine molindfix
!
!  Determines cell indices for molecules
!
!   3/95 Determination of molecule dimensionality added
!   8/95 This version of molind created which preserves
!        molecules and bonds regardless of bond cutoff
!   2/01 lintoijk calls now use imaxl,jmaxl and kmaxl
!   4/01 Modified to handle fixed connectivity
!   5/04 Modified to use inbox coordinates if lspatialok is true
!   4/07 Searching for bonded atoms now uses nbonds number rather
!        than checking for a non-zero type
!   4/07 Handling of case where connectivity is input for several
!        images of same atom by checking when a connection has
!        already been used.
!   5/07 Call to connectwrap excluded if .not.lmodco
!   6/07 Algorithm to find atom nearest to origin changed to be 
!        more robust for large distances.
!   5/10 Timing added
!   5/10 Modified to accelerate assignment of nmolind
!   2/18 Trace added
!   2/20 Use of inbox coordinates for spatial removed
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
  use molecule
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
  integer(i4)                                  :: ixx
  integer(i4)                                  :: iyy
  integer(i4)                                  :: izz
  integer(i4)                                  :: j
  integer(i4)                                  :: jj
  integer(i4)                                  :: k
  integer(i4)                                  :: kk
  integer(i4)                                  :: na
  integer(i4)                                  :: nact
  integer(i4)                                  :: nassign
  integer(i4)                                  :: nasslast
  integer(i4)                                  :: nbon
  integer(i4)                                  :: nbond
  integer(i4)                                  :: nbondedatom
  integer(i4)                                  :: ni
  integer(i4)                                  :: ninc
  integer(i4)                                  :: nj
  integer(i4)                                  :: nj1
  integer(i4)                                  :: nmi
  integer(i4)                                  :: nmlind
  integer(i4)                                  :: nnew
  integer(i4)                                  :: nnewnext
  integer(i4)                                  :: nno
  integer(i4)                                  :: npn
  integer(i4)                                  :: nptr
  integer(i4)                                  :: npy
  integer(i4), dimension(:), allocatable       :: npnew
  integer(i4), dimension(:), allocatable       :: npnewnext
  integer(i4), dimension(:), allocatable       :: npno
  integer(i4), dimension(:), allocatable       :: npno2numat
  integer(i4), dimension(:), allocatable       :: npyes
  integer(i4)                                  :: nyes
  integer(i4)                                  :: status
  logical,     dimension(:), allocatable       :: lAlreadyUsed
  logical                                      :: lfindconnect
  logical                                      :: lfound
  logical                                      :: lmolconn
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: r2
  real(dp)                                     :: r2min
  real(dp)                                     :: ri
  real(dp)                                     :: rij
  real(dp)                                     :: rj
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
  call trace_in('molindfix')
#endif
!
  t1 = g_cpu_time()
!
!  Allocate local memory
!
  call realloc(iimptr,maxbond,ierror)
  if (ierror.ne.0) call outofmemory('molindfix','iimptr')
  allocate(npnew(numat),stat=status)
  if (status/=0) call outofmemory('molindfix','npnew')
  allocate(npnewnext(numat),stat=status)
  if (status/=0) call outofmemory('molindfix','npnewnext')
  allocate(npno(numat),stat=status)
  if (status/=0) call outofmemory('molindfix','npno')
  allocate(npno2numat(numat),stat=status)
  if (status/=0) call outofmemory('molindfix','npno2numat')
  allocate(npyes(numat),stat=status)
  if (status/=0) call outofmemory('molindfix','npyes')
  if (nconnect.gt.0) then
    allocate(lAlreadyUsed(nconnect),stat=status)
    if (status/=0) call outofmemory('molind','lAlreadyUsed')
  endif
!
!  Handle periodic wrapping of connection indices
!
  if (nconnect.gt.0.and.lmodco) call connectwrap
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
    j = nbonded(1,i)
    lfindconnect = (ri.ne.0.0_dp)
    lmolconn = .false.
    if (nmi.ne.0) then
      if (nmolconnect(nmi).ne.0) then
        lfindconnect = .true.
        lmolconn = .true.
      endif
    endif
    if (lfindconnect) then
!
!  Loop over bonds
!
      do while (j.gt.0.and.nbond.le.nbonds(i))
        lfound = .false.
        if (lmolconn) then
          lAlreadyUsed(1:nmolconnect(nmi)) = .false.
          icm = 0
          do while (icm.lt.nmolconnect(nmi).and..not.lfound)
            icm = icm + 1
            ic = nmolconnectptr(icm,nmi)
            if (nconnectcfg(ic).eq.ncf.and..not.lAlreadyUsed(ic)) then
              if (n1connect(ic).eq.i.and.n2connect(ic).eq.j) then
                nbonded(nbond,i) = n2connect(ic)
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
                lfound = .true.
                lAlreadyUsed(icm) = .true.
                iimptr(nbond) = 0
              elseif (n2connect(ic).eq.i.and.n1connect(ic).eq.j) then
                nbonded(nbond,i) = n1connect(ic)
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
                lfound = .true.
                lAlreadyUsed(icm) = .true.
              endif
            endif
          enddo
        endif
        if (.not.lfound) then
          xcd = xclat(j)
          ycd = yclat(j)
          zcd = zclat(j)
          nj1 = nat(j)
          nj = nj1
          if (nj.gt.maxele) nj = nj - maxele
          rj = rcov(nj)
          if (rj.ne.0.0_dp) then
            xcdi = xcd - xal
            ycdi = ycd - yal
            zcdi = zcd - zal
!
!  Find shortest distance to bonded atom image
!  excluding those already located
!
            rmin = 1.0d10
            do ii = 1,iimax
              xcrd = xcdi + xvec1cell(ii)
              ycrd = ycdi + yvec1cell(ii)
              zcrd = zcdi + zvec1cell(ii)
              rij = xcrd*xcrd + ycrd*ycrd + zcrd*zcrd
              if (rij.le.rmin.and.rij.gt.0.0_dp.and.(i.ne.j.or.ii.ne.iimid)) then
!
!  Check to see whether this image has already been found
!
                lfound = .false.
                nbon = nbond - 1
                do while (.not.lfound.and.nbon.gt.0) 
                  if (nbonded(nbon,i).eq.j.and.iimptr(nbon).eq.ii) then
                    lfound = .true.
                  endif
                  nbon = nbon - 1
                enddo
!
!  Valid bond
!
                if (.not.lfound) then
                  iimptr(nbond) = ii
                  rmin = rij
                endif
              endif
            enddo
            ii = iimptr(nbond)
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
        endif
        nbond = nbond + 1
        if (nbond.gt.maxbond) then
          maxbond = nbond + 2
          call changemaxbond
          call realloc(iimptr,maxbond,ierror)
          if (ierror.ne.0) call outofmemory('molindfix','iimptr')
        endif
        j = nbonded(nbond,i)
      enddo
    endif
  enddo
  do i = 1,numat
    nmolind(i) = 0
  enddo
!**************************************
!  Obtain cell indices for each atom  *
!**************************************
!
!  Initialisation of npno2numat
!
  npno2numat(1:numat) = 0
!
!  Loop over molecules
!
  do i = 1,nmol
    if (nmolatom(i).gt.0) then
!
!  Locate most central atom of molecule
!
      ninc = nmolptr(i)
      do j = 1,nmolatom(i)
        na = nmollist(ninc+j)
        if (ndim.eq.3) then
          rij = (xfrac(na) - 0.5)**2 + (yfrac(na)-0.5)**2 + (zfrac(na)-0.5)**2
        elseif (ndim.eq.2) then
          rij = (xfrac(na) - 0.5)**2 + (yfrac(na)-0.5)**2
        elseif (ndim.eq.1) then
          rij = (xfrac(na) - 0.5)**2
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
!  npno       = list of atoms yet to be found
!  npyes      = list of atoms that have been found
!  npnew      = list of atoms found in the last cycle
!  npnewnext  = list of atoms found in the current cycle that will become npnew in the next cycle
!
      nact = nmollist(ninc+nptr)
      nmolind(nact) = 555
      nassign = 1
      nyes = 1
      npyes(1) = nptr
      nnew = 1
      npnew(1) = nptr
      nno = 0
      do j = 1,nmolatom(i)
        if (j.ne.nptr) then
          nno = nno + 1
          npno(nno) = j
          npno2numat(nmollist(ninc+j)) = nno
        endif
      enddo
10    nasslast = nassign
!
!  Reset list of new atoms for next cycle
!
      nnewnext = 0
!
!  Loop over new atoms from last cycle to look for bonds
!
      do j = 1,nnew
        npy = nmollist(ninc+npnew(j))
        do k = 1,nbonds(npy)
          nbondedatom = nbonded(k,npy)
!
!  Is this atom in the no list?
!
          if (npno2numat(nbondedatom).ne.0) then
            npn = npno2numat(nbondedatom)
!
!  Is index still zero? (A previous atom in this shell might have bonded it)
!
            if (nmolind(nbondedatom).eq.0) then
!
!  Find index number
!
              nmlind = nmolind(npy) + nbondind(k,npy) - 555
              nmolind(nbondedatom) = nmlind
              nassign = nassign + 1
!
!  Add to list of new atoms for next cycle
!
              nnewnext = nnewnext + 1
              npnewnext(nnewnext) = npno(npn)
            endif
          endif
        enddo
!
!  End loop over unindexed sites
!
      enddo
!
!  Reinitialise npno2numat using old npno data
!
      do j = 1,nno
        npno2numat(nmollist(ninc+npno(nno))) = 0
      enddo
!
!  Copy nnewnext lists to nnew lists
!
      do j = 1,nnewnext
        npnew(j) = npnewnext(j)
      enddo
      nnew = nnewnext
!
!  Rebuild list of no and yes atoms
!
      nno = 0
      nyes = 0
      do j = 1,nmolatom(i)
        if (nmolind(nmollist(ninc+j)).eq.0) then
          nno = nno + 1
          npno(nno) = j
          npno2numat(nmollist(ninc+j)) = nno
        else
          nyes = nyes + 1
          npyes(nyes) = j
        endif
      enddo
!
!  Have we assigned all the atoms yet? If not, then loop again
!
      if (nassign.ne.nmolatom(i).and.nassign.ne.nasslast) goto 10
    endif
!
!  End of loop over molecules
!
  enddo
!
!  Free local memory
!
  if (nconnect.gt.0) then
    deallocate(lAlreadyUsed,stat=status)
    if (status/=0) call deallocate_error('molind','lAlreadyUsed')
  endif
  deallocate(npyes,stat=status)
  if (status/=0) call deallocate_error('molindfix','npyes')
  deallocate(npno2numat,stat=status)
  if (status/=0) call deallocate_error('molindfix','npno2numat')
  deallocate(npno,stat=status)
  if (status/=0) call deallocate_error('molindfix','npno')
  deallocate(npnewnext,stat=status)
  if (status/=0) call deallocate_error('molindfix','npnewnext')
  deallocate(npnew,stat=status)
  if (status/=0) call deallocate_error('molindfix','npnew')
  call realloc(iimptr,0_i4,ierror)
!
  t2 = g_cpu_time()
  tmol = tmol + t2 - t1
#ifdef TRACE
  call trace_out('molindfix')
#endif
!
  return
  end
