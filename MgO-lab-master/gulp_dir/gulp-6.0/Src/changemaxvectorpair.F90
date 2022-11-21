  subroutine changemaxvectorpair(vec,maxdim)
!
!  Alters the size of the arrays within the type vector_pair
!
!   4/09 Created
!   1/14 Modified to include option to nullify pointers
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
!  Copyright Curtin University 2014
!
!  Julian Gale, CIC, Curtin University, January 2014
!
  use vectors,      only : vector_pair
  use reallocate
  implicit none
!
!  Passed variables
!
  type(vector_pair), intent(inout) :: vec
  integer(i4),       intent(in)    :: maxdim
!
!  Local variables
!
  integer(i4)       :: ierror
!
!  if maxdim =< 0 then this is an initialisation call
!
  if (maxdim.le.0) then
    vec%maxdim_pair = 0
    nullify(vec%distance_pair1)
    nullify(vec%distance_pair2)
    nullify(vec%xvector_pair1)
    nullify(vec%yvector_pair1)
    nullify(vec%zvector_pair1)
    nullify(vec%xvector_pair2)
    nullify(vec%yvector_pair2)
    nullify(vec%zvector_pair2)
  else
!
    vec%maxdim_pair = maxdim
!
    call realloc(vec%distance_pair1,maxdim,ierror)
    if (ierror.ne.0) call outofmemory('changemaxvectorpair','distance_pair1')
    call realloc(vec%distance_pair2,maxdim,ierror)
    if (ierror.ne.0) call outofmemory('changemaxvectorpair','distance_pair2')
    call realloc(vec%xvector_pair1,maxdim,ierror)
    if (ierror.ne.0) call outofmemory('changemaxvectorpair','xvector_pair1')
    call realloc(vec%yvector_pair1,maxdim,ierror)
    if (ierror.ne.0) call outofmemory('changemaxvectorpair','yvector_pair1')
    call realloc(vec%zvector_pair1,maxdim,ierror)
    if (ierror.ne.0) call outofmemory('changemaxvectorpair','zvector_pair1')
    call realloc(vec%xvector_pair2,maxdim,ierror)
    if (ierror.ne.0) call outofmemory('changemaxvectorpair','xvector_pair2')
    call realloc(vec%yvector_pair2,maxdim,ierror)
    if (ierror.ne.0) call outofmemory('changemaxvectorpair','yvector_pair2')
    call realloc(vec%zvector_pair2,maxdim,ierror)
    if (ierror.ne.0) call outofmemory('changemaxvectorpair','zvector_pair2')
  endif
!
  return
  end
