  subroutine closestimage(x,y,z)
!
!  Find nearest image of x, y, z
!
!  On entry :
!
!    x = X component of vector as supplied
!    y = Y component of vector as supplied
!    z = Z component of vector as supplied
!
!  On exit :
!
!    x = X component of nearest image to origin
!    y = Y component of nearest image to origin
!    z = Z component of nearest image to origin
!
!   8/01 Created
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
  use current
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  real(dp)    :: x
  real(dp)    :: y
  real(dp)    :: z
!
!  Local variables
!
  integer(i4) :: ii
  real(dp)    :: r
  real(dp)    :: rmin
  real(dp)    :: xmin
  real(dp)    :: ymin
  real(dp)    :: zmin
  real(dp)    :: xd
  real(dp)    :: yd
  real(dp)    :: zd
#ifdef TRACE
  call trace_in('closestimage')
#endif
!
  rmin = 1.0d10
  xmin = x
  ymin = y
  zmin = z
!
!  Loop over unit cells
!
  do ii = 1,iimax
    xd = x + xvec1cell(ii)
    yd = y + yvec1cell(ii)
    zd = z + zvec1cell(ii)
    r = xd*xd + yd*yd + zd*zd
    if (r .lt. rmin) then
      rmin = r
      xmin = xd
      ymin = yd
      zmin = zd
    endif
  enddo
  x = xmin
  y = ymin
  z = zmin
#ifdef TRACE
  call trace_out('closestimage')
#endif
!
  return
  end
