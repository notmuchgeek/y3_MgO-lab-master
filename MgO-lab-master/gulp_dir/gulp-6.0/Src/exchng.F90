  subroutine linmin_exchng(a,b,x,y,t,q,n)
!
!  The contents of a, c, t, and x are stored in b, d, q, and y!
!  This is a dedicated routine, it is called by linmin only.
!
  use datatypes
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)  :: n
  real(dp),    intent(in)  :: a
  real(dp),    intent(out) :: b
  real(dp),    intent(in)  :: t
  real(dp),    intent(out) :: q
  real(dp),    intent(in)  :: x(*)
  real(dp),    intent(out) :: y(*)
!
!  Local variables
!
  integer(i4)    :: i
#ifdef TRACE
  call trace_in('exchng')
#endif
!
  b = a
  q = t
  do i = 1,n
    y(i) = x(i)
  enddo
#ifdef TRACE
  call trace_out('exchng')
#endif
!
  return
  end
