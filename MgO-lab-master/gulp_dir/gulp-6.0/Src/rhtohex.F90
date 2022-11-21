  subroutine rhtohex
!
!  Convert rhombohedral cell to hexagonal form
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
  use configurations
  use g_constants
  use current
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Local variables
!
  integer(i4)     :: i
  integer(i4)     :: j
  real(dp)        :: ang
  real(dp)        :: bcosphi
  real(dp)        :: r
  real(dp)        :: rroot3
#ifdef TRACE
  call trace_in('rhtohex')
#endif
!
!  Work out hexagonal cell parameters
!
!  a(hex) = 2*a(rh)*sin(alpha/2)
!
  ang = 0.5_dp*alpha/radtodeg
  r = a*sin(ang)
  a = 2.0_dp*r
!
!  c(hex)=3*a(rh)*cos(phi) , phi = angle to cell diagonal
!
  rroot3 = 1.0_dp/sqrt(3.0_dp)
  r = 2.0_dp*r*rroot3
  bcosphi = b*b - r*r
  bcosphi = sqrt(bcosphi)
  c = 3.0_dp*bcosphi
  b = a
  alpha = 90.0_dp
  beta = 90.0_dp
  gamma = 120.0_dp
!
!  Convert cell parameters to cell
!
  call cell3D(rv,a,b,c,alpha,beta,gamma)
  do i = 1,3
    do j = 1,3
      rvcfg(j,i,ncf) = rv(j,i)
    enddo
  enddo
#ifdef TRACE
  call trace_out('rhtohex')
#endif
!
  return
  end
