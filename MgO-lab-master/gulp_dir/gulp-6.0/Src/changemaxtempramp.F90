  subroutine changemaxtempramp
!
!  Alters the size of the arrays associated with maxtempramp
!
!   3/19 Created from changemaxcfg
!   3/19 Multiple temperature ramps added
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
  use configurations
  use current
  use reallocate
  implicit none
!
!  Local variables
!
  integer(i4)       :: ierror, i
  integer(i4), save :: oldmaxtempramp = 0
!
!  Configuration data
!
  call realloc(ntemperaturestep,maxtempramp,ierror)
  if (ierror.ne.0) call outofmemory('changemaxtempstp','ntemperaturestep')
  call realloc(ntemperaturestepstart,maxtempramp,ierror)
  if (ierror.ne.0) call outofmemory('changemaxtempstp','ntemperaturestepstart')
  call realloc(ntemperaturestepstop,maxtempramp,ierror)
  if (ierror.ne.0) call outofmemory('changemaxtempstp','ntemperaturestepstop')
  call realloc(ntempstp,maxtempramp,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxtempstp','ntempstp')
  call realloc(ntempstpstart,maxtempramp,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxtempstp','ntempstpstart')
  call realloc(tempstp,maxtempramp,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxtempstp','tempstp')
  call realloc(temperaturestart,maxtempramp,ierror)
  if (ierror.ne.0) call outofmemory('changemaxtempstp','temperaturestart')
  call realloc(temperaturestep,maxtempramp,ierror)
  if (ierror.ne.0) call outofmemory('changemaxtempstp','temperaturestep')
!
!  Initialise defaults for new part of array
!
  if (maxtempramp.gt.oldmaxtempramp) then
    do i = oldmaxtempramp+1,maxtempramp
      ntempstp(i,1:maxcfg) = 0
      ntempstpstart(i,1:maxcfg) = 0
      tempstp(i,1:maxcfg) = 0.0_dp
    enddo
  endif
!
!  Save current value of maxtempstp for next call
!
  oldmaxtempramp = maxtempramp
!
  return
  end
