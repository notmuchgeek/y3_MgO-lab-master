  subroutine destroy(natom)
!
!  Destroys an atom and removes all data relating to it
!
!  On entry :
!
!  natom = atom in asymmetric unit whose symmetry related
!          images and itself will be removed
!
!   6/01 Pointer for moving atoms corrected
!  10/02 Correction to condensing of nbonded added
!   5/06 Handling of nspecptr added
!   4/07 Searching for bonded atoms now uses nbonds number rather
!        than checking for a non-zero type
!  12/08 scrho, and defect versions, converted to 2-D array for MEAM
!   5/11 Algorithm changed for handling of nbonds/nbonded/nbondind
!        to correct a bug.
!   3/17 Modifications to allow for new ordering of variables
!   2/18 Trace added
!   5/18 nbondqb added
!   2/19 x0 removed
!   3/19 iopt changed to ioptindex and iopttype
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!  10/19 Langevin damping of dipoles added
!   7/20 Arguments passed to collect have (1) removed
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
  use cellmultipole,  only : nboxat
  use configurations, only : lopfi
  use current
  use eam,            only : lMEAM, maxmeamcomponent
  use molecule,       only : natmol, nmolind
  use optimisation
  use polarise,       only : dpolar, qpolar, dpolarmax
  use potentialxyz
  use shells,         only : ncsptr, nshptr
  use sutton,         only : scrho, scrho12
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  use velocities 
  implicit none
!
!  Passed variables
!
  integer(i4)                            :: natom
!
!  Local variables
!
  integer(i4)                            :: neqvatom
  integer(i4)                            :: i
  integer(i4)                            :: ii
  integer(i4)                            :: ind
  integer(i4)                            :: j
  integer(i4)                            :: k
  integer(i4), dimension(:), allocatable :: iptr,itmp,itmp2
  integer(i4)                            :: status
  logical                                :: lreduce_nbonds
  logical,     dimension(:), allocatable :: ltmp
  logical,     dimension(:), allocatable :: ltmp2
  real(dp),    dimension(:), allocatable :: rtmp, rtmp2
#ifdef TRACE
  call trace_in('destroy')
#endif
!
!  Allocate local memory
!
  allocate(iptr(numat),stat=status)
  if (status/=0) call outofmemory('destroy','iptr')
  allocate(itmp(numat),stat=status)
  if (status/=0) call outofmemory('destroy','itmp')
  allocate(itmp2(numat),stat=status)
  if (status/=0) call outofmemory('destroy','itmp2')
  allocate(ltmp(numat),stat=status)
  if (status/=0) call outofmemory('destroy','ltmp')
  allocate(ltmp2(numat),stat=status)
  if (status/=0) call outofmemory('destroy','ltmp2')
  allocate(rtmp(numat),stat=status)
  if (status/=0) call outofmemory('destroy','rtmp')
  allocate(rtmp2(numat),stat=status)
  if (status/=0) call outofmemory('destroy','rtmp2')
!
!  Find all images in full cell and sort to the end
!
  neqvatom = neqv(natom)
  do i = 1,natom-1
    iptr(i) = i
  enddo
  do i = natom,numat-1
    iptr(i) = i + 1
  enddo
  iptr(numat) = natom
!
!  Sort bonding data before compression
!
!  Remove atoms no longer present from bonding list
!
  do i = 1,numat
    lreduce_nbonds = .false.
    do j = 1,nbonds(i)
      if (nbonded(j,i).eq.natom) then
        do k = j+1,nbonds(i)
          nbonded(k-1,i) = nbonded(k,i)
          nbondqb(k-1,i) = nbondqb(k,i)
          nbondind(k-1,i) = nbondind(k,i)
        enddo
        lreduce_nbonds = .true.
      endif
    enddo
    if (lreduce_nbonds) then
      nbonds(i) = nbonds(i) - 1
    endif
  enddo
!
!  Move information in nbonded before nbonds is sorted
!
  do i = natom+1,numat
    do j = 1,nbonds(i)
      nbonded(j,i-1) = nbonded(j,i)
      nbondqb(j,i-1) = nbondqb(j,i)
      nbondind(j,i-1) = nbondind(j,i)
    enddo
  enddo
!
!  Sort data
!
  call collect(numat,c6f,rtmp,iptr)
  call collect(numat,cnf,rtmp,iptr)
  call icollect(numat,icosx,itmp,iptr)
  call icollect(numat,icosy,itmp,iptr)
  call icollect(numat,icosz,itmp,iptr)
  call collect(numat,mass,rtmp,iptr)
  call icollect(numat,nat,itmp,iptr)
  call icollect(numat,natmol,itmp,iptr)
  call icollect(numat,nbonds,itmp,iptr)
  call icollect(numat,nboxat,itmp,iptr)
  call icollect(numat,ncsptr,itmp,iptr)
  call icollect(numat,nftype,itmp,iptr)
  call icollect(numat,nmolind,itmp,iptr)
  call icollect(numat,nrelf2a,itmp,iptr)
  call icollect(numat,nrotop,itmp,iptr)
  call icollect(numat,nshptr,itmp,iptr)
  call collect(numat,occuf,rtmp,iptr)
  call collect(numat,oxf,rtmp,iptr)
  call collect(numat,qf,rtmp,iptr)
  call collect(numat,radf,rtmp,iptr)
  call collect(numat,rmass,rtmp,iptr)
  call collect(numat,velx,rtmp,iptr)
  call collect(numat,vely,rtmp,iptr)
  call collect(numat,velz,rtmp,iptr)
  call collect(numat,xclat,rtmp,iptr)
  call collect(numat,yclat,rtmp,iptr)
  call collect(numat,zclat,rtmp,iptr)
  call collect(numat,xfrac,rtmp,iptr)
  call collect(numat,yfrac,rtmp,iptr)
  call collect(numat,zfrac,rtmp,iptr)
  call collect(numat,x2,rtmp,iptr)
  call collect(numat,y2,rtmp,iptr)
  call collect(numat,z2,rtmp,iptr)
  call collect(numat,x3,rtmp,iptr)
  call collect(numat,y3,rtmp,iptr)
  call collect(numat,z3,rtmp,iptr)
  call collect(numat,x4,rtmp,iptr)
  call collect(numat,y4,rtmp,iptr)
  call collect(numat,z4,rtmp,iptr)
  call collect(numat,x5,rtmp,iptr)
  call collect(numat,y5,rtmp,iptr)
  call collect(numat,z5,rtmp,iptr)
!
!  If atom in bonded list is greater than natom, then shift number down by 1
!
  do i = 1,numat-1
    do j = 1,nbonds(i)
      if (nbonded(j,i).gt.natom) nbonded(j,i) = nbonded(j,i) - 1
    enddo
  enddo
!
!  Move data to overwrite atom in asymmetric unit
!
  do i = 1,nasym-1
    iptr(i) = i
  enddo
  do i = natom,nasym-1
    iptr(i) = i + 1
  enddo
  iptr(nasym) = natom
!
!  Sort data
!
  call collect(nasym,c6a,rtmp,iptr)
  call collect(nasym,cna,rtmp,iptr)
  call icollect(nasym,iatn,itmp,iptr)
  call lcollect(nasym,lopf,ltmp,iptr)
  call icollect(nasym,natype,itmp,iptr)
  call icollect(nasym,nspecptr,itmp,iptr)
  call icollect(nasym,neqv,itmp,iptr)
  call icollect(nasym,nrela2f,itmp,iptr)
  call collect(nasym,occua,rtmp,iptr)
  call collect(nasym,oxa,rtmp,iptr)
  call collect(nasym,dpolar,rtmp,iptr)
  call collect(nasym,dpolarmax,rtmp,iptr)
  call collect(nasym,qpolar,rtmp,iptr)
  call collect(nasym,qa,rtmp,iptr)
  call collect(nasym,rada,rtmp,iptr)
  call collect(nasym,vx,rtmp,iptr)
  call collect(nasym,vy,rtmp,iptr)
  call collect(nasym,vz,rtmp,iptr)
  call collect(nasym,vx12,rtmp,iptr)
  call collect(nasym,vy12,rtmp,iptr)
  call collect(nasym,vz12,rtmp,iptr)
!
  if (lMEAM) then
    do j = 1,maxmeamcomponent
      do i = 1,nasym
        rtmp2(i) = scrho(j,i)
      enddo
      call collect(nasym,rtmp2,rtmp,iptr)
      do i = 1,nasym
        scrho(j,i) = rtmp(i)
        rtmp2(i) = scrho12(j,i)
      enddo
      call collect(nasym,rtmp2,rtmp,iptr)
      do i = 1,nasym
        scrho12(j,i) = rtmp(i)
      enddo
    enddo
  else
    do i = 1,nasym
      rtmp2(i) = scrho(1,i)
    enddo
    call collect(nasym,rtmp2,rtmp,iptr)
    do i = 1,nasym
      scrho(1,i) = rtmp(i)
      rtmp2(i) = scrho12(1,i)
    enddo
    call collect(nasym,rtmp2,rtmp,iptr)
    do i = 1,nasym
      scrho12(1,i) = rtmp(i)
    enddo
  endif
!
  do i = 1,6
    do j = 1,nasym
      rtmp2(j) = v2xyz(i,j)
    enddo
    call collect(nasym,rtmp2,rtmp,iptr)
    do j = 1,nasym
      v2xyz(i,j) = rtmp(j)
    enddo
    do j = 1,nasym
      rtmp2(j) = v2xyz12(i,j)
    enddo
    call collect(nasym,rtmp2,rtmp,iptr)
    do j = 1,nasym
      v2xyz12(i,j) = rtmp(j)
    enddo
  enddo
  call collect(nasym,xalat,rtmp,iptr)
  call collect(nasym,yalat,rtmp,iptr)
  call collect(nasym,zalat,rtmp,iptr)
  call collect(nasym,xstore,rtmp,iptr)
  call collect(nasym,ystore,rtmp,iptr)
  call collect(nasym,zstore,rtmp,iptr)
  call collect(nasym,rstore,rtmp,iptr)
!
!  Optimisation variables - sort lopfi and rebuild iopt
!
  lopfi(1:3*nasym) = .false.
  do i = 1,nvar
    ii = ioptindex(i)
    if (iopttype(i).eq.iopt_xf) then
      lopfi(3*ii-2) = .true.
    elseif (iopttype(i).eq.iopt_yf) then
      lopfi(3*ii-1) = .true.
    elseif (iopttype(i).eq.iopt_zf) then
      lopfi(3*ii) = .true.
    endif
  enddo
  do i = 1,nasym
    ltmp2(i) = lopfi(3*i-2)
  enddo
  call lcollect(nasym,ltmp2(1),ltmp,iptr)
  do i = 1,nasym
    lopfi(3*i-2) = ltmp2(i)
  enddo
  do i = 1,nasym
    ltmp2(i) = lopfi(3*i-1)
  enddo
  call lcollect(nasym,ltmp2,ltmp,iptr)
  do i = 1,nasym
    lopfi(3*i-1) = ltmp2(i)
  enddo
  do i = 1,nasym
    ltmp2(i) = lopfi(3*i)
  enddo
  call lcollect(nasym,ltmp2,ltmp,iptr)
  do i = 1,nasym
    lopfi(3*i) = ltmp2(i)
  enddo
  if (loldvarorder) then
!
!  Old variable order
!
    nvar = ncell
  else
!
!  New variable order
!
    nvar = 0
  endif
  ind = 0
  do i = 1,nasym
    do j = 1,3
      ind = ind + 1
      if (lopfi(ind)) then
        nvar = nvar + 1
        ioptindex(nvar) = i
        if (j.eq.1) then
          iopttype(nvar) = iopt_xf
        elseif (j.eq.2) then
          iopttype(nvar) = iopt_yf
        else
          iopttype(nvar) = iopt_zf
        endif
      endif
    enddo
  enddo
!
!  Coordinate variables
!
  call collect(nasym,xafrac,rtmp,iptr)
  call collect(nasym,yafrac,rtmp,iptr)
  call collect(nasym,zafrac,rtmp,iptr)
!
  if (nbsm.gt.0) then
!
!  Reorder radii
!
    call collect(nasym,rada,rtmp,iptr)
!
!  Shuffle radial optimisation variables
!
    lopfi(1:nasym) = .false.
    do i = 1,nvar
      ii = ioptindex(i)
      if (iopttype(i).eq.iopt_radius) then
        lopfi(ii) = .true.
      endif
    enddo
    call lcollect(nasym,lopfi,ltmp,iptr)
!
    do i = 1,nasym
      if (lopfi(i)) then
        nvar = nvar + 1
        ioptindex(nvar) = i
        iopttype(nvar) = iopt_radius
      endif
    enddo
  endif
!
!  For new variable order add the cell parameters to the end of the variable list
!  
  if (.not.loldvarorder) then
    ninternalmax = nvar
    do i = ncellmin,ncellmax
      nvar = nvar + 1
      ioptindex(nvar) = ioptindex(i)
      iopttype(nvar)  = iopttype(i)
    enddo
    ncellmin = ninternalmax + 1
    ncellmax = ncellmin + ncell - 1
  endif
!
!  Correct totals
!
  nasym = nasym - 1
  numat = numat - neqvatom
!
!  Free local memory
!
  deallocate(rtmp2,stat=status)
  if (status/=0) call deallocate_error('destroy','rtmp2')
  deallocate(rtmp,stat=status)
  if (status/=0) call deallocate_error('destroy','rtmp')
  deallocate(ltmp2,stat=status)
  if (status/=0) call deallocate_error('destroy','ltmp2')
  deallocate(ltmp,stat=status)
  if (status/=0) call deallocate_error('destroy','ltmp')
  deallocate(itmp2,stat=status)
  if (status/=0) call deallocate_error('destroy','itmp2')
  deallocate(itmp,stat=status)
  if (status/=0) call deallocate_error('destroy','itmp')
  deallocate(iptr,stat=status)
  if (status/=0) call deallocate_error('destroy','iptr')
#ifdef TRACE
  call trace_out('destroy')
#endif
!
  return
  end
