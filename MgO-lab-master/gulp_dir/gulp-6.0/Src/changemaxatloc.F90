  subroutine changemaxatloc
!
!  Alters the size of the arrays associated with maxatloc
!
!   2/17 Created from changemaxat
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
!  Julian Gale, CIC, Curtin University, February 2017
!
  use current
  use derivatives,    only : nqatoms, nqatomcell, nqatomptr, d2qdxyz2, d2qdxyzs, d2qds2, maxqatoms
  use derivatives,    only : qatomxyz, dqds
  use reallocate
  implicit none
!
!  Local variables
!
  integer(i4)       :: ierror, i, maxqatoms2
  integer(i4), save :: oldmaxatloc = 0
!
  call realloc(dqds,6_i4,maxatloc,ierror)
  if (ierror.ne.0) call outofmemory('changemaxatloc','dqds')
  call realloc(d2qdxyzs,6_i4,3_i4*maxqatoms,maxatloc,ierror)
  if (ierror.ne.0) call outofmemory('changemaxatloc','d2qdxyzs')
  maxqatoms2 = (3_i4*maxqatoms + 3_i4)*(3_i4*maxqatoms + 6_i4)/2_i4
  call realloc(d2qdxyz2,maxqatoms2,maxatloc,ierror)
  if (ierror.ne.0) call outofmemory('changemaxatloc','d2qdxyz2')
  call realloc(d2qds2,21_i4,maxatloc,ierror)
  if (ierror.ne.0) call outofmemory('changemaxatloc','d2qds2')
  call realloc(nqatomcell,3_i4,maxqatoms,maxatloc,ierror)
  if (ierror.ne.0) call outofmemory('changemaxatloc','nqatomcell')
  call realloc(nqatomptr,maxqatoms,maxatloc,ierror)
  if (ierror.ne.0) call outofmemory('changemaxatloc','nqatomptr')
  call realloc(nqatoms,maxatloc,ierror)
  if (ierror.ne.0) call outofmemory('changemaxatloc','nqatoms')
  call realloc(qatomxyz,3_i4,maxqatoms,maxatloc,ierror)
  if (ierror.ne.0) call outofmemory('changemaxatloc','qatomxyz')
!
!  Initialise new parts of data arrays
!
  if (maxatloc.gt.oldmaxatloc) then
    do i = oldmaxatloc+1,maxatloc
      nqatoms(i) = 0
    enddo
  endif
!
!  Save current value of maxatloc for next call
!
  oldmaxatloc = maxatloc
!
  return
  end
