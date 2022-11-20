  subroutine mindtoijk(indin,ii,jj,kk)
!
!  Converts molecule index number to x,y and z component indices
!
!   2/18 Trace added
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
!  Julian Gale, CIC, Curtin University, Feb 2018
!
  use datatypes
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)  :: indin
  integer(i4), intent(out) :: ii
  integer(i4), intent(out) :: jj
  integer(i4), intent(out) :: kk
!
!  Local variables
!
  integer(i4)              :: ind
#ifdef TRACE
  call trace_in('mindtoijk')
#endif
!
  ind = indin
  kk = (ind/100)
  ind = ind - 100*kk
  jj = (ind/10)
  ii = ind - 10*jj - 5
  jj = jj - 5
  kk = kk - 5
#ifdef TRACE
  call trace_out('mindtoijk')
#endif
!
  return
  end
