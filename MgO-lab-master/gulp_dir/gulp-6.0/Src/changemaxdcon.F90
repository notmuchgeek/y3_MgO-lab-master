  subroutine changemaxdcon
!
!  Alters the size of the arrays associated with maxdcon
!
!   3/19 Change of defect constraints to have index and type
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
!  Julian Gale, CIC, Curtin University, March 2019
!
  use defects
  use fitting
  use reallocate
  implicit none
!
  integer(i4)       :: ierror
!
!  Configuration data
!
  call realloc(dconco,maxdcon,ierror)
  if (ierror.ne.0) call outofmemory('changemaxdcon','dconco')
  call realloc(ncdfixind,maxdcon,ierror)
  if (ierror.ne.0) call outofmemory('changemaxdcon','ncdfixind')
  call realloc(ncdfixtyp,maxdcon,ierror)
  if (ierror.ne.0) call outofmemory('changemaxdcon','ncdfixtyp')
  call realloc(ncdvarind,maxdcon,ierror)
  if (ierror.ne.0) call outofmemory('changemaxdcon','ncdvarind')
  call realloc(ncdvartyp,maxdcon,ierror)
  if (ierror.ne.0) call outofmemory('changemaxdcon','ncdvartyp')
!
  return
  end
