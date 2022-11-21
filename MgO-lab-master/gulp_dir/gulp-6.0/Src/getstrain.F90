!****************
!  3-D systems  *
!****************
  subroutine getstrain3D(rv,xstr)
!
!  Computes the strain required to map original unit cell to current unit cell
!
!   6/18 Created
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
  use configurations, only : rvcfg
  use current,        only : ncf
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  real(dp),   intent(in)    :: rv(3,3)
  real(dp),   intent(out)   :: xstr(6)
!
!  Local variables
!
  integer(i4)               :: i
  integer(i4)               :: ifail
  integer(i4)               :: j
  integer(i4)               :: k
  real(dp)                  :: rvinv(3,3)
  real(dp)                  :: strmat(3,3)
  real(dp)                  :: wrk(6)
#ifdef TRACE
  call trace_in('getstrain3D')
#endif
!
!  Compute the inverse of the lattice vectors for the reference configuration
!
  rvinv(1:3,1:3) = rvcfg(1:3,1:3,ncf)
  call matrix_inversion(rvinv,3_i4,3_i4,wrk,ifail)
!
!  Multiply the matrices together to find the strain matrix
!
  do i = 1,3
    do j = 1,3
      strmat(j,i) = 0.0_dp
      do k = 1,3
       strmat(j,i) = strmat(j,i) + rv(j,k)*rvinv(k,i)
      enddo
    enddo
  enddo
!
!  Assign strains 
!
  xstr(1) = strmat(1,1) - 1.0_dp
  xstr(2) = strmat(2,2) - 1.0_dp
  xstr(3) = strmat(3,3) - 1.0_dp
  xstr(4) = strmat(2,3) + strmat(3,2)
  xstr(5) = strmat(1,3) + strmat(3,1)
  xstr(6) = strmat(1,2) + strmat(2,1)
#ifdef TRACE
  call trace_out('getstrain3D')
#endif
!
  return
  end
!****************
!  2-D systems  *
!****************
  subroutine getstrain2D(rv,xstr)
!
!  Computes the strain required to map original surface unit cell to current surface unit cell
!
!   6/18 Created
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
  use configurations, only : rvcfg
  use current,        only : ncf
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  real(dp),   intent(in)    :: rv(3,2)
  real(dp),   intent(out)   :: xstr(3)
!
!  Local variables
!
  integer(i4)               :: i
  integer(i4)               :: ifail
  integer(i4)               :: j
  integer(i4)               :: k
  real(dp)                  :: rvinv(3,2)
  real(dp)                  :: strmat(2,2)
  real(dp)                  :: wrk(4)
#ifdef TRACE
  call trace_in('getstrain2D')
#endif
!
!  Compute the inverse of the lattice vectors for the reference configuration
!
  rvinv(1:3,1:2) = rvcfg(1:3,1:2,ncf)
  call matrix_inversion(rvinv,3_i4,2_i4,wrk,ifail)
!
!  Multiply the matrices together to find the strain matrix
!
  do i = 1,2
    do j = 1,2
      strmat(j,i) = 0.0_dp
      do k = 1,2
       strmat(j,i) = strmat(j,i) + rv(j,k)*rvinv(k,i)
      enddo
    enddo
  enddo
!
!  Assign strains 
!
  xstr(1) = strmat(1,1) - 1.0_dp
  xstr(2) = strmat(2,2) - 1.0_dp
  xstr(3) = strmat(1,2) + strmat(2,1)
#ifdef TRACE
  call trace_out('getstrain2D')
#endif
!
  return
  end
!****************
!  1-D systems  *
!****************
  subroutine getstrain1D(rv,xstr)
!
!  Computes the strain required to map original unit cell to current unit cell
!
!   6/18 Created
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
!  Julian Gale, Curtin University, June 2018
!
  use datatypes
  use configurations, only : rvcfg
  use current,        only : ncf
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
! 
  real(dp),   intent(in)    :: rv(3,1)
  real(dp),   intent(out)   :: xstr(1)
#ifdef TRACE
  call trace_in('getstrain1D')
#endif
!
!  Apply strain
!
  xstr(1) = rv(1,1)/rvcfg(1,1,ncf) - 1.0_dp
#ifdef TRACE
  call trace_out('getstrain1D')
#endif
!
  return
  end
