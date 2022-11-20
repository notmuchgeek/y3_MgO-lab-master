  subroutine changemaxqrange
!
!  Alters the size of the arrays associated with maxqrange
!
!   5/18 Created
!   6/18 e0range added
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
!  Julian Gale, CIC, Curtin University, June 2018
!
  use element,      only : maxele
  use eemdata
  use reallocate
  implicit none
!
  integer(i4)       :: ierror, i
  integer(i4), save :: oldmaxqrange = 0
!
!  Electronegativity data for charge ranges
!
  call realloc(nqrange,maxele,maxeemtype,ierror)
  if (ierror.ne.0) call outofmemory('changemaxqrange','nqrange')
  call realloc(nqrangetype,maxqrange,maxele,maxeemtype,ierror)
  if (ierror.ne.0) call outofmemory('changemaxqrange','nqrangetype')
  call realloc(qrangemax,maxqrange,maxele,maxeemtype,ierror)
  if (ierror.ne.0) call outofmemory('changemaxqrange','qrangemax')
  call realloc(qrangemin,maxqrange,maxele,maxeemtype,ierror)
  if (ierror.ne.0) call outofmemory('changemaxqrange','qrangemin')
  call realloc(e0range,maxqrange,maxele,maxeemtype,ierror)
  if (ierror.ne.0) call outofmemory('changemaxqrange','e0range')
  call realloc(q0range,maxqrange,maxele,maxeemtype,ierror)
  if (ierror.ne.0) call outofmemory('changemaxqrange','q0range')
  call realloc(chirange,maxqrange,maxele,maxeemtype,ierror)
  if (ierror.ne.0) call outofmemory('changemaxqrange','chirange')
  call realloc(murange,maxqrange,maxele,maxeemtype,ierror)
  if (ierror.ne.0) call outofmemory('changemaxqrange','murange')
  call realloc(radrange,maxqrange,maxele,maxeemtype,ierror)
  if (ierror.ne.0) call outofmemory('changemaxqrange','radrange')
  call realloc(znucrange,maxqrange,maxele,maxeemtype,ierror)
  if (ierror.ne.0) call outofmemory('changemaxqrange','znucrange')
  call realloc(zetarange,maxqrange,maxele,maxeemtype,ierror)
  if (ierror.ne.0) call outofmemory('changemaxqrange','zetarange')
!
!  Initialise defaults for new part of array
!
  if (maxqrange.gt.oldmaxqrange) then
    do i = oldmaxqrange+1,maxqrange
      call initmaxqrangedefaults(i)
    enddo
  endif
!
!  Save current value of maxqrange for next call
!
  oldmaxqrange = maxqrange
!
  return
  end
