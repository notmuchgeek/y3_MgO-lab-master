  subroutine setrotation(axes,drR)
!
!  Subroutine to setup derivatives of the rotation matrix about an axis 
!  NB: These are all evaluated for theta = 0
!
!   4/20 Created
!   5/20 Second derivatives removed as they are not needed at present
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
!  Copyright Curtin University 2020
!
!  Julian Gale, CIC, Curtin University, May 2020
!
  use datatypes
#ifdef TRACE
  use trace,            only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  real(dp),      intent(in)   :: axes(3,3)       ! Set of 3 axes about which to compute the rotation
  real(dp),      intent(out)  :: drR(3,3,3)      ! First derivatives of the rotation matrix
!
!  Local variables
!
  integer(i4)                 :: i
#ifdef TRACE
  call trace_in('setrotation')
#endif
!
!  Loop over axes
!
  do i = 1,3
!
!  First derivative with respect to an axis
!
    drR(1,1,i) = 0.0_dp
    drR(2,1,i) = - axes(3,i)
    drR(3,1,i) = axes(2,i)
    drR(1,2,i) = axes(3,i)
    drR(2,2,i) = 0.0_dp
    drR(3,2,i) = - axes(1,i)
    drR(1,3,i) = - axes(2,i)
    drR(2,3,i) = axes(1,i)
    drR(3,3,i) = 0.0_dp
  enddo
!
#ifdef TRACE
  call trace_out('setrotation')
#endif
!
  return
  end
!
  subroutine getrotationmatrix(axis,theta,rotm)
!
!  Subroutine to setup the rotation matrix about an axis 
!
!   5/20 Created
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
!  Copyright Curtin University 2020
!
!  Julian Gale, CIC, Curtin University, May 2020
!
  use datatypes
#ifdef TRACE
  use trace,            only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  real(dp),      intent(in)   :: axis(3)      ! Rotation axis
  real(dp),      intent(in)   :: theta        ! Angle of rotation in radians
  real(dp),      intent(out)  :: rotm(3,3)    ! Rotation matrix
!
!  Local variables
!
  real(dp)                    :: costh
  real(dp)                    :: sinth
#ifdef TRACE
  call trace_in('getrotationmatrix')
#endif
  costh = cos(theta)
  sinth = sin(theta)
!
  rotm(1,1) = costh + (axis(1)**2)*(1.0_dp - costh)
  rotm(2,1) = axis(1)*axis(2)*(1.0_dp - costh) - axis(3)*sinth
  rotm(3,1) = axis(1)*axis(3)*(1.0_dp - costh) + axis(2)*sinth
  rotm(1,2) = axis(2)*axis(1)*(1.0_dp - costh) + axis(3)*sinth
  rotm(2,2) = costh + (axis(2)**2)*(1.0_dp - costh)
  rotm(3,2) = axis(2)*axis(3)*(1.0_dp - costh) - axis(1)*sinth
  rotm(1,3) = axis(3)*axis(1)*(1.0_dp - costh) - axis(2)*sinth
  rotm(2,3) = axis(3)*axis(2)*(1.0_dp - costh) + axis(1)*sinth
  rotm(3,3) = costh + (axis(3)**2)*(1.0_dp - costh)
!
#ifdef TRACE
  call trace_out('getrotationmatrix')
#endif
!
  return
  end
