  subroutine changemaxnboA
!
!  Alters the size of the arrays associated with maxnboA
!
!   6/10 Option to for second atom specifier added to attractive terms.
!   9/10 Initialisations now performed in a subroutine
!   1/14 lBOzrlA added
!   9/15 BOccoeffA made into a 2-D array
!  11/20 Tersoff reorganised
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
  integer(i4), save :: oldmaxnboA = 0
!
!  Bond order potential data
!
  call realloc(BOccoeffA,5_i4,maxnboA,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnboA','BOccoeffA')
  call realloc(BOhcoeffA,maxnboA,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnboA','BOhcoeffA')
  call realloc(BOlcoeffA,maxnboA,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnboA','BOlcoeffA')
  call realloc(BOmcoeffA,maxnboA,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnboA','BOmcoeffA')
  call realloc(BOocoeffA,maxnboA,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnboA','BOocoeffA')
  call realloc(nBOspecA1,maxnboA,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnboA','nBOspecA1')
  call realloc(nBOspecA2,maxnboA,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnboA','nBOspecA2')
  call realloc(nBOtypA1,maxnboA,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnboA','nBOtypA1')
  call realloc(nBOtypA2,maxnboA,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnboA','nBOtypA2')
  call realloc(nBOtypeA,maxnboA,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnboA','nBOtypeA')
  call realloc(lBOzrlA,maxnboA,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnboA','lBOzrlA')
!
!  Initialise defaults for new part of array
!
  if (maxnboA.gt.oldmaxnboA) then
    do i = oldmaxnboA+1,maxnboA
      call initmaxnboadefaults(i)
    enddo
  endif
!
!  Save current value of maxnboA for next call
!
  oldmaxnboA = maxnboA
!
  return
  end
