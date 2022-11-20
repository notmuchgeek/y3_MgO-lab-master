  subroutine fefunctn(iflag,n,xc,fc,gc,hessian,nhwords,lhess2D)
!
!  Supplies the function and first derivatives of the free energy
!  using numerical derivatives. Main purposes is for checking
!  analytical derivatives.
!
!   8/97 Created from funct to replace old routine from 
!        numerical version.
!   9/97 If iflag=2 then need to keep second derivatives for
!        hessian calculation
!  11/07 Unused variables removed
!   3/17 nhwords and lhess2D added as arguments
!   1/18 Trace added
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
!  Julian Gale, CIC, Curtin University, January 2018
!
  use control
  use general
#ifdef TRACE
  use trace,       only : trace_in, trace_out
#endif
  implicit none
!
!  Passed arguments
!
  integer(i4), intent(inout) :: iflag
  integer(i4), intent(in)    :: n
  integer(i4), intent(in)    :: nhwords
  logical,     intent(in)    :: lhess2D
  real(dp),    intent(out)   :: fc
  real(dp),    intent(out)   :: xc(*)
  real(dp),    intent(out)   :: gc(*)
  real(dp),    intent(out)   :: hessian(nhwords,*)
!
!  Local variables
!
  integer(i4)                :: i
  real(dp)                   :: g_cpu_time
  real(dp)                   :: fcb
  real(dp)                   :: fcf
  real(dp),             save :: tdmax = 0.0_dp
  real(dp)                   :: t1
  real(dp)                   :: t2
  real(dp)                   :: xci
#ifdef TRACE
  call trace_in('fefunctn')
#endif
!
  t1 = g_cpu_time()
!
!  Derivative flags
!
  lfirst = .true.
!**************************
!  Calculate free energy  *
!**************************
  call fefunct(0_i4,n,xc,fc,gc,hessian,nhwords,lhess2D)
!**********************************
!  Calculate numerical gradients  *
!**********************************
  do i = 1,n
    xci = xc(i)
!
!  Forward
!
    xc(i) = xci + findiff
    call fefunct(0_i4,n,xc,fcf,gc,hessian,nhwords,lhess2D)
!
!  Backward
!
    xc(i) = xci - findiff
    call fefunct(0_i4,n,xc,fcb,gc,hessian,nhwords,lhess2D)
!
!  Calculate gradient
!
    gc(i) = (fcf-fcb)/(2.0_dp*findiff)
!
!  Restore xc
!
    xc(i) = xci
  enddo
!************************************
!  Calculate free energy for return *
!************************************
  call fefunct(0_i4,n,xc,fc,gc,hessian,nhwords,lhess2D)
!*******************
!  CPU time check  *
!*******************
  t2 = g_cpu_time()
  if (t2-t1.gt.tdmax) tdmax = t2 - t1
  if (timmax.gt.0.0_dp) then
    if ((timmax-t2).lt.tdmax) iflag = -1
  endif
#ifdef TRACE
  call trace_out('fefunctn')
#endif
!
  return
  end
