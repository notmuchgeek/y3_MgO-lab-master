  subroutine changemaxcon
!
!  Alters the size of the arrays associated with maxcon
!
!   3/19 Constraint arrays changed to have index and type
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
  use current
  use reallocate
  implicit none
!
  integer(i4)       :: ierror
!
!  Configuration data
!
  call realloc(conadd,maxcon,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcon','conadd')
  call realloc(conco,maxcon,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcon','conco')
  call realloc(ncfixind,maxcon,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcon','ncfixind')
  call realloc(ncfixtyp,maxcon,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcon','ncfixtyp')
  call realloc(ncvarind,maxcon,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcon','ncvarind')
  call realloc(ncvartyp,maxcon,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcon','ncvartyp')
!
  return
  end
