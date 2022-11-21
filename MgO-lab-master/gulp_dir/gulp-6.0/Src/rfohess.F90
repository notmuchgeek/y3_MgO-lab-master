  subroutine rfohess(hessian,maxhess,lhess2D,nvar)
!
!  Form hessian from symmetry reduced second derivative matrix
!
!   6/05 Intent added
!   2/17 nmin removed from arguments
!   2/17 maxhess & lhess2D added as arguments
!   3/17 Parallelisation added
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
  use derivatives
  use parallel
#ifdef TRACE
  use trace,       only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)  :: nvar
  integer(i4), intent(in)  :: maxhess
  logical,     intent(in)  :: lhess2D
  real(dp),    intent(out) :: hessian(maxhess,*)
!
!  Local variables
!
  integer(i4)              :: i
  integer(i4)              :: ind
  integer(i4)              :: j
#ifdef TRACE
  call trace_in('rfohess')
#endif
!
!  Transfer hessian to correct matrix
!
  if (lhess2D) then
    do i = 1,nvaronnode
      do j = 1,nvar
        hessian(j,i) = derv2(j,i)
      enddo
    enddo
  else
    ind = 0
    do i = 1,nvar
      do j = 1,i
        ind = ind + 1
        hessian(ind,1) = derv2(j,i)
      enddo
    enddo
  endif
#ifdef TRACE
  call trace_out('rfohess')
#endif
!
  return
  end
