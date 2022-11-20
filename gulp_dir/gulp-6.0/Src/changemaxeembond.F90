  subroutine changemaxeembond
!
!  Alters the size of the arrays associated with maxeembond
!
!   4/18 Created from changemaxbondQ
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
!  Copyright Curtin University 2018
!
!  Julian Gale, CIC, Curtin University, April 2018
!
  use eembonds
  use reallocate
  implicit none
!
!  Local variables
!
  integer(i4)       :: ierror
!
!  Split bond charge data
!
  call realloc(neembonded,2_i4,maxeembond,ierror)
  if (ierror.ne.0) call outofmemory('changemaxeembond','neembonded')
  call realloc(qbond,maxeembond,ierror)
  if (ierror.ne.0) call outofmemory('changemaxeembond','qbond')
!
  return
  end
