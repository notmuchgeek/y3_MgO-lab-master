  subroutine collect(n,x,y,iptr)
!
!  Shuffle location of floating point variables
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
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)    :: n
  real(dp),    intent(inout) :: x(*)
  real(dp),    intent(out)   :: y(*)
  integer(i4), intent(in)    :: iptr(*)
!
!  Local variables
!
  integer(i4)                :: i
#ifdef TRACE
  call trace_in('collect')
#endif
!
  do i = 1,n
    y(i) = x(i)
  enddo
  do i = 1,n
    x(i) = y(iptr(i))
  enddo
#ifdef TRACE
  call trace_out('collect')
#endif
  return
  end
!
  subroutine icollect(n,ix,iy,iptr)
!
!  Shuffle location of integer variables
!
  use datatypes
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)    :: n
  integer(i4), intent(inout) :: ix(*)
  integer(i4), intent(out)   :: iy(*)
  integer(i4), intent(in)    :: iptr(*)
!
!  Local variables
!
  integer(i4) :: i
#ifdef TRACE
  call trace_in('icollect')
#endif
!
  do i = 1,n
    iy(i) = ix(i)
  enddo
  do i = 1,n
    ix(i) = iy(iptr(i))
  enddo
#ifdef TRACE
  call trace_out('icollect')
#endif
  return
  end
!
  subroutine lcollect(n,lx,ly,iptr)
!
!  Shuffle location of logical variables
!
  use datatypes
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)    :: n
  logical,     intent(inout) :: lx(*)
  logical,     intent(out)   :: ly(*)
  integer(i4), intent(in)    :: iptr(*)
!
!  Local variables
!
  integer(i4) :: i
#ifdef TRACE
  call trace_in('lcollect')
#endif
!
  do i = 1,n
    ly(i) = lx(i)
  enddo
  do i = 1,n
    lx(i) = ly(iptr(i))
  enddo
#ifdef TRACE
  call trace_out('lcollect')
#endif
  return
  end
