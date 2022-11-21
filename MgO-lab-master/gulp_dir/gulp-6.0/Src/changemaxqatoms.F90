  subroutine changemaxqatoms
!
!  Alters the size of the arrays associated with maxqatoms
!
!   1/09 Integer datatypes all explicitly declared
!   2/15 nqatomcell and qatomxyz added
!   2/17 maxat changed to maxatloc 
!   7/17 maxat no longer explicitly imported
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
  use current,        only : maxatloc
  use derivatives,    only : maxqatoms, nqatomptr, d2qdxyzs, nqatomcell, qatomxyz
  use reallocate
  implicit none
!
!  Local variables
!
  integer(i4)       :: ierror
!
  call realloc(nqatomcell,3_i4,maxqatoms,maxatloc,ierror)
  if (ierror.ne.0) call outofmemory('changemaxqatoms','nqatomcell')
  call realloc(nqatomptr,maxqatoms,maxatloc,ierror)
  if (ierror.ne.0) call outofmemory('changemaxqatoms','nqatomptr')
  call realloc(d2qdxyzs,6_i4,3_i4*maxqatoms,maxatloc,ierror)
  if (ierror.ne.0) call outofmemory('changemaxqatoms','d2qdxyzs')
  call realloc(qatomxyz,3_i4,maxqatoms,maxatloc,ierror)
  if (ierror.ne.0) call outofmemory('changemaxqatoms','qatomxyz')
!
  return
  end
