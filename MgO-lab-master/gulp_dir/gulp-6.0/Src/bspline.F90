  recursive subroutine bspline_M(n,u,mnu)
!
!  Compute the b-spline function M_n(u)
!
!   5/16 Created
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
!  Julian Gale, Curtin University, Feb 2018
!
  use datatypes
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)  :: n            ! Order of B-spline
  real(dp),    intent(in)  :: u            ! Real number for B-spline
  real(dp),    intent(out) :: mnu          ! Return value of M_n(u)
!
!  Local variables
!
  integer(i4)              :: n1           ! n - 1 for passing to next level
  real(dp)                 :: mnu1         ! M_n-1(u)
  real(dp)                 :: mnu2         ! M_n-1(u-1)
  real(dp)                 :: u1           ! u - 1
!
!  Check bounds of u
!
  if (u.le.0.0_dp.or.u.ge.dble(n)) then
    mnu = 0.0_dp
    return
  endif
#ifdef TRACE
  call trace_in('bspline')
#endif
!
!  Set u - 1
!
  u1 = u - 1.0_dp
!
  if (n.eq.2) then
!
!  If n = 2 then this is the lowest level, so compute
!
    mnu = 1.0_dp - abs(u1)
  else
!
!  If n > 2 then use recursion
!
    n1 = n - 1
    call bspline_M(n1,u,mnu1)
    call bspline_M(n1,u1,mnu2)
!
    mnu = (u*mnu1 + (dble(n) - u)*mnu2)/dble(n1)
  endif
#ifdef TRACE
  call trace_out('bspline')
#endif
!
  return
  end
