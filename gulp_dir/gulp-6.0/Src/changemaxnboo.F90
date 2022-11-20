  subroutine changemaxnboo
!
!  Alters the size of the arrays associated with maxnboO
!
!  11/20 Created from changemaxnbopot
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
!  Julian Gale, CIC, Curtin University, November 2020
!
  use bondorderdata
  use reallocate
  implicit none
!
!  Local variables
!
  integer(i4)       :: ierror, i
  integer(i4), save :: oldmaxnboO = 0
!
!  Bond order potential data
!
  call realloc(BOecoeffA,maxnboO,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnboO','BOecoeffA')
  call realloc(BOecoeffR,maxnboO,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnboO','BOecoeffR')
  call realloc(BOncoeffA,maxnboO,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnboO','BOncoeffA')
  call realloc(BOncoeffR,maxnboO,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnboO','BOncoeffR')
  call realloc(nBOspec0,maxnboO,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnboO','nBOspec0')
  call realloc(nBOtyp0,maxnboO,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnboO','nBOtyp0')
!
!  Initialise defaults for new part of array
!
  if (maxnboO.gt.oldmaxnboO) then
    do i = oldmaxnboO+1,maxnboO
      call initmaxnboOdefaults(i)
    enddo
  endif
!
!  Save current value of maxnbopot for next call
!
  oldmaxnboO = maxnboO
!
  return
  end
