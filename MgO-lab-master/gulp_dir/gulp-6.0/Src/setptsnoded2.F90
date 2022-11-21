  subroutine setptsnoded2
!
!  Sets up mapping between SAS points and nodes for use in distributed second derivative calculations
!
!   4/17 Created from setvarnoded2
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
!  Julian Gale, CIC, Curtin University, April 2017
!
  use cosmic
  use current
  use parallel
  use reallocate
  implicit none
!
!  Local variables
!
  integer(i4)                                       :: i
  integer(i4)                                       :: icount
  integer(i4)                                       :: node
!
!  Determine parallel distribution of SAS points for block cyclic form
!
  if (nprocs.gt.1) then
    nptsonnode = 0
    node = 0
    icount = 0
!
!  If block size hasn't been input then choose a value based on the number of atoms versus processors
!
    if (nblocksizesas.eq.0) then
      nblocksizesas = 12
    endif
!
    do i = 1,npts
      icount = icount + 1
      npts2node(i) = node
      if (node.eq.procid) then
        nptsonnode = nptsonnode + 1
        node2pts(nptsonnode) = i
        npts2local(i) = nptsonnode
      else
        npts2local(i) = 0
      endif
      if (icount.eq.nblocksizesas) then
        icount = 0
        node = node + 1
        if (node.eq.nprocs) node = 0
      endif
    enddo
  else
    nptsonnode = npts
    do i = 1,npts
      npts2node(i) = 0
      node2pts(i) = i
      npts2local(i) = i
    enddo
  endif
!
  return
  end
