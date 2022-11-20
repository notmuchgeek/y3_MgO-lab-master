  subroutine setatomnoded2
!
!  Sets up mapping between atoms and nodes for use in distributed second derivative calculations
!
!   7/12 Created from code in phonon
!  11/12 Re-merged and separated from phonon atom pointer distribution
!   9/13 Name changed to distinguish this routine from setatomnodes
!   9/13 Calculation of cores and shells on node added
!   9/16 Definition of ncoonnodeptr & nshonnodeptr changed
!   2/17 maxatloc added
!   2/17 Atom centric arrays substituted for overall distribution ones
!   7/17 Check that blocksize is a factor of the number of particles
!        added otherwise scalapack may object
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
!  Copyright Curtin University 2017
!
!  Julian Gale, CIC, Curtin University, July 2017
!
  use current
  use element,     only : maxele
  use parallel
  use shells
  implicit none
!
!  Local variables
!
  integer(i4)                                       :: i
  integer(i4)                                       :: icount
  integer(i4)                                       :: ii
  integer(i4)                                       :: node
!
!  Having determined the global information now construct parallel decomposition if needed
!
!  1) Divide cores that will be needed for second derivative calculation in a block cyclic decomposition
!  2) Assign shells of these cores to the same node
!
  if (nprocs.gt.1) then
    natomsonnodea = 0
    node = 0
    icount = 0
!
    if (nblocksize.eq.0) then
!
!  If block size hasn't been input then choose a value based on the number of atoms versus processors
!
      nblocksize = 1
    else
!
!  Check that blocksize is a factor of the number of atoms
!
      if (mod(numat,nblocksize).ne.0) then
        call outerror('Blocksize is not a factor of the number of atoms',0_i4)
        call stopnow('setatomnoded2')
      endif
!
!  Check that blocksize is a factor of the number of cores
!
      if (mod(ncore,nblocksize).ne.0) then
        call outerror('Blocksize is not a factor of the number of cores',0_i4)
        call stopnow('setatomnoded2')
      endif
!
!  Check that blocksize is a factor of the number of shells
!
      if (nshell.gt.0) then
        if (mod(nshell,nblocksize).ne.0) then
          call outerror('Blocksize is not a factor of the number of shells',0_i4)
          call stopnow('setatomnoded2')
        endif
      endif
    endif
!
    do i = 1,numat
      icount = icount + 1
      atom2nodea(i) = node
      if (node.eq.procid) then
        natomsonnodea = natomsonnodea + 1
        node2atoma(natomsonnodea) = i
        atom2locala(i) = natomsonnodea
      else
        atom2locala(i) = 0
      endif
      if (icount.eq.nblocksize) then
        icount = 0
        node = node + 1
        if (node.eq.nprocs) node = 0
      endif
    enddo
!
!  Set up shell pointer array for properties section
!
    ncoreonnode  = 0
    nshellonnode = 0
    do ii = 1,natomsonnodea
      i = node2atoma(ii)
      if (nat(i).gt.maxele) then
        nshellonnode = nshellonnode + 1
        nshonnodeptr(nshellonnode) = ii
      else
        ncoreonnode = ncoreonnode + 1
        ncoonnodeptr(ncoreonnode) = ii
      endif
    enddo
  else
    natomsonnodea = numat
    do i = 1,numat
      atom2nodea(i) = 0
      node2atoma(i) = i
      atom2locala(i) = i
    enddo
    ncoreonnode = ncore
    nshellonnode = nshell
    ncoonnodeptr(1:ncore)  = ncoptr(1:ncore)
    nshonnodeptr(1:nshell) = nshptr(1:nshell)
  endif
!
!  Update maxatloc
!
  if (natomsonnodea.gt.maxatloc) then
    maxatloc = natomsonnodea
    call changemaxatloc
  endif
!
  return
  end
