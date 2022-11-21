  subroutine changemaxd1
!
!  Alters the size of the arrays associated with maxd1 or maxd1u
!
!   3/09 Non-radial arrays added
!   6/09 Virial arrays added
!   4/12 xvir, yvir and zvir removed
!  11/18 Non-radial arrays removed since they are no longer needed
!   5/20 xfdrv, yfdrv, zfdrv added
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
  use derivatives
  use reallocate
  implicit none
!
  integer(i4) :: ierror
!
  call realloc(raderv,maxd1,ierror)
  if (ierror.ne.0) call outofmemory('changemaxd1','raderv')
  call realloc(xdrv,maxd1,ierror)
  if (ierror.ne.0) call outofmemory('changemaxd1','xdrv')
  call realloc(ydrv,maxd1,ierror)
  if (ierror.ne.0) call outofmemory('changemaxd1','ydrv')
  call realloc(zdrv,maxd1,ierror)
  if (ierror.ne.0) call outofmemory('changemaxd1','zdrv')
  call realloc(xfdrv,maxd1,ierror)
  if (ierror.ne.0) call outofmemory('changemaxd1','xfdrv')
  call realloc(yfdrv,maxd1,ierror)
  if (ierror.ne.0) call outofmemory('changemaxd1','yfdrv')
  call realloc(zfdrv,maxd1,ierror)
  if (ierror.ne.0) call outofmemory('changemaxd1','zfdrv')
!
  return
  end
