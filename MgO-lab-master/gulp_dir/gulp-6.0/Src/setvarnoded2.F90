  subroutine setvarnoded2(ldefect)
!
!  Sets up mapping between variables and nodes for use in distributed second derivative calculations
!  Also creates an alternative mapping for atoms that focuses on variables
!
!   2/17 Created from setatomnoded2
!   3/17 node2atom corrected to node2atomv
!   3/17 ncellonnode added
!   3/17 ninternalonnode added
!   5/17 Input flag added to indicate whether this is a defect call or not
!   5/17 Parallel defect modifications added to connect variables to atoms
!   5/17 Most of routine removed for defects as it is not needed
!   3/19 iopt replaced by ioptindex and iopttype
!   5/20 Rigid molecule modifications added
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
!  Julian Gale, CIC, Curtin University, May 2020
!
  use current
  use optimisation
  use parallel
  use reallocate
  implicit none
!
!  Passed variables
!
  logical,                         intent(in)       :: ldefect
!
!  Local variables
!
  integer(i4)                                       :: i
  integer(i4)                                       :: icount
  integer(i4)                                       :: ityp
  integer(i4)                                       :: ivar
  integer(i4)                                       :: node
  integer(i4)                                       :: status
  logical,          dimension(:), allocatable       :: latomneeded
!
!  Determine parallel distribution of variables for block cyclic form
!
  if (nprocs.gt.1) then
    nvaronnode = 0
    node = 0
    icount = 0
!
!  If block size hasn't been input then choose a value based on the number of atoms versus processors
!
    if (nblocksizevar.eq.0) then
      nblocksizevar = 3
    endif
!
    do i = 1,nvar
      icount = icount + 1
      nvar2node(i) = node
      if (node.eq.procid) then
        nvaronnode = nvaronnode + 1
        node2var(nvaronnode) = i
        nvar2local(i) = nvaronnode
      else
        nvar2local(i) = 0
      endif
      if (icount.eq.nblocksizevar) then
        icount = 0
        node = node + 1
        if (node.eq.nprocs) node = 0
      endif
    enddo
  else
    nvaronnode = nvar
    do i = 1,nvar
      nvar2node(i) = 0
      node2var(i) = i
      nvar2local(i) = i
    enddo
  endif
  if (ldefect) then
!********************************************
!  Defect mapping of variables to region 1  *
!********************************************
  else
!******************************************
!  Bulk mapping of variables to region 1  *
!******************************************
    allocate(latomneeded(numat),stat=status)
    if (status/=0) call outofmemory('setvarnoded2','latomneeded')
!
!  Having determined the distribution of variables work out the atoms that 
!  need to be local to the current processor.
!
    if (nprocs.gt.1) then
      latomneeded(1:numat) = .false.
      ncellonnode = 0
      ncellmaxonnode = 0
      ncellminonnode = 0
      ninternalonnode = 0
      ninternalatmonnode = 0
      ninternalmaxonnode = 0
      ninternalminonnode = 0
      do i = 1,nvaronnode
        ityp = iopttype(node2var(i))
        ivar = ioptindex(node2var(i))
        if (ityp.eq.iopt_radius) then
!
!  Breathing core or shell radius
!
          ninternalonnode = ninternalonnode + 1
          latomneeded(ivar) = .true.
          if (ninternalminonnode.eq.0) ninternalminonnode = i
          ninternalmaxonnode = i
        elseif (ityp.eq.iopt_xf.or.ityp.eq.iopt_yf.or.ityp.eq.iopt_zf) then
!
!  Core or shell coordinate
!
          ninternalonnode = ninternalonnode + 1
          ninternalatmonnode = ninternalatmonnode + 1
          latomneeded(ivar) = .true.
          if (ninternalminonnode.eq.0) ninternalminonnode = i
          ninternalmaxonnode = i
       elseif (ityp.eq.iopt_xcom.or.ityp.eq.iopt_ycom.or.ityp.eq.iopt_zcom) then
!
!  Centre of mass coordinate
!
          ninternalonnode = ninternalonnode + 1
          latomneeded(ivar) = .true.
          if (ninternalminonnode.eq.0) ninternalminonnode = i
          ninternalmaxonnode = i
       elseif (ityp.eq.iopt_xqtn.or.ityp.eq.iopt_yqtn.or.ityp.eq.iopt_zqtn) then
!
!  Quaternion coordinate
!
          ninternalonnode = ninternalonnode + 1
          latomneeded(ivar) = .true.
          if (ninternalminonnode.eq.0) ninternalminonnode = i
          ninternalmaxonnode = i
        elseif (ityp.eq.iopt_strain.or.ityp.eq.iopt_cell) then
!
!  Strain or cell
!
          ncellonnode = ncellonnode + 1
          if (ncellminonnode.eq.0) ncellminonnode = i
          ncellmaxonnode = i
        endif
      enddo
    else
      ncellonnode = ncell
      ncellmaxonnode = ncellmax
      ncellminonnode = ncellmin
      ninternalonnode = ninternal
      ninternalatmonnode = ninternalatm
      ninternalmaxonnode = ninternalmax
      ninternalminonnode = ninternalmin
      latomneeded(1:numat) = .true.
    endif
!
!  Build atom lists for variables
!
    natomsonnodev = 0
    do i = 1,numat
      if (latomneeded(i)) then
        natomsonnodev = natomsonnodev + 1
        atom2nodev(i) = procid             ! This won't be unique since several nodes may share an atom
        atom2localv(i) = natomsonnodev
        node2atomv(natomsonnodev) = i
      endif
    enddo
!
    deallocate(latomneeded,stat=status)
    if (status/=0) call deallocate_error('setvarnoded2','latomneeded')
!
!  Update maxatloc
!
    if (natomsonnodev.gt.maxatloc) then
      maxatloc = natomsonnodev
      call changemaxatloc
    endif
  endif
!
  return
  end
