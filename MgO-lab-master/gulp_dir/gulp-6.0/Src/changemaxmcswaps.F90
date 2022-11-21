  subroutine changemaxmcswaps
!
!  Alters the size of the arrays associated with maxmcswaps
!
!   5/16 Created from changemaxmcswapspec
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
!  Julian Gale, CIC, Curtin University, May 2016
!
  use montecarlo
  use reallocate
  implicit none
!
  integer(i4)       :: ierror, i
  integer(i4), save :: oldmaxmcswaps   = 0
!
  call realloc(nmcswapnat,maxmcswapspec,maxmcswaps,ierror)
  if (ierror.ne.0) call outofmemory('changemaxmcswaps','nmcswapnat')
  call realloc(nmcswaptype,maxmcswapspec,maxmcswaps,ierror)
  if (ierror.ne.0) call outofmemory('changemaxmcswaps','nmcswaptype')
  call realloc(nmcswapspec,maxmcswaps,ierror)
  if (ierror.ne.0) call outofmemory('changemaxmcswaps','nmcswapspec')
  call realloc(nmcswappair,maxmcswaps,ierror)
  if (ierror.ne.0) call outofmemory('changemaxmcswaps','nmcswappair')
  call realloc(nswapable,maxmcswaps,ierror)
  if (ierror.ne.0) call outofmemory('changemaxmcswaps','nswapable')
  call realloc(pswap,maxmcswaps,ierror)
  if (ierror.ne.0) call outofmemory('changemaxmcswaps','pswap')
  call realloc(lmcswapany,maxmcswaps,ierror)
  if (ierror.ne.0) call outofmemory('changemaxmcswaps','lmcswapany')
!
!  Initialise new parts of data arrays
!
  if (maxmcswaps.gt.oldmaxmcswaps) then
    do i = oldmaxmcswaps+1,maxmcswaps
      call initmaxmcswapsdefaults(i)
    enddo
  endif
!
!  Save current value of maxmcswaps for next call
!
  oldmaxmcswaps = maxmcswaps
!
  return
  end
