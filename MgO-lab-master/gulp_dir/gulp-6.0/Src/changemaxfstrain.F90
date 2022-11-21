  subroutine changemaxfstrain
!
!  Alters the size of the arrays associated with maxfstrain
!
!  12/13 Created from changemaxfstress
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
!  Copyright Curtin University 2013
!
!  Julian Gale, CIC, Curtin University, December 2013
!
  use configurations
  use observables
  use reallocate
  implicit none
!
  integer(i4)       :: ierror, i
  integer(i4), save :: oldmaxfstrain = 0
!
!  Configuration data
!
  call realloc(nfstraincfg,maxfstrain,ierror)
  if (ierror.ne.0) call outofmemory('changemaxfstrain','nfstraincfg')
  call realloc(nfstraint,maxfstrain,ierror)
  if (ierror.ne.0) call outofmemory('changemaxfstrain','nfstraint')
  call realloc(fstrain,maxfstrain,ierror)
  if (ierror.ne.0) call outofmemory('changemaxfstrain','fstrain')
  call realloc(fstrainweight,maxfstrain,ierror)
  if (ierror.ne.0) call outofmemory('changemaxfstrain','fstrainweight')
!
!  Initialise new parts of data arrays
!
  if (maxfstrain.gt.oldmaxfstrain) then
    do i = oldmaxfstrain+1,maxfstrain
      call initmaxfstraindefaults(i)
    enddo
  endif
!
!  Save current value of maxfstrain for next call
!
  oldmaxfstrain = maxfstrain
!
  return
  end
