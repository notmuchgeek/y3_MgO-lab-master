  subroutine changemaxnboZ
!
!  Alters the size of the arrays associated with maxnboZ
!
!   3/16 Created from changemaxnboA
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
!  Copyright Curtin University 2016
!
!  Julian Gale, CIC, Curtin University, March 2016
!
  use bondorderdata
  use reallocate
  implicit none
!
!  Local variables
!
  integer(i4)       :: ierror, i
  integer(i4), save :: oldmaxnboZ = 0
!
!  Bond order potential data
!
  call realloc(BOccoeffZ,2_i4,maxnboZ,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnboZ','BOccoeffZ')
  call realloc(BOecoeffZ,maxnboZ,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnboZ','BOecoeffZ')
  call realloc(BOzcoeffZ,maxnboZ,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnboZ','BOzcoeffZ')
  call realloc(nBOspecZ,maxnboZ,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnboZ','nBOspecZ')
  call realloc(nBOtypZ,maxnboZ,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnboZ','nBOtypZ')
!
!  Initialise defaults for new part of array
!
  if (maxnboZ.gt.oldmaxnboZ) then
    do i = oldmaxnboZ+1,maxnboZ
      call initmaxnbozdefaults(i)
    enddo
  endif
!
!  Save current value of maxnboA for next call
!
  oldmaxnboZ = maxnboZ
!
  return
  end
