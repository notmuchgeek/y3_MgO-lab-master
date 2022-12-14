  subroutine getdipole2D(dipolez)
!
!  Calculate the dipole moment in the z direction for a surface.
!
!   2/01 Created
!   2/08 Corrected for partial occupancy
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
!  Julian Gale, CIC, Curtin University, February 2018
!
  use current
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  real(dp)         :: dipolez
!
!  Local variables
!
  integer(i4)      :: i
#ifdef TRACE
  call trace_in('getdipole2D')
#endif
!
!  Loop over atoms and find dipole
!
  dipolez = 0.0_dp
  do i = 1,numat
    dipolez = dipolez + occuf(i)*qf(i)*zclat(i)
  enddo
#ifdef TRACE
  call trace_out('getdipole2D')
#endif
  return
  end
!
  subroutine getdipole1D(dipoley,dipolez)
!
!  Calculate the dipole moment in the y and z directions for a polymer.
!
!   2/01 Created
!   2/08 Corrected for partial occupancy
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
!  Julian Gale, CIC, Curtin University, February 2018
!
  use current
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
!
  implicit none
!
!  Passed variables
!
  real(dp)         :: dipoley
  real(dp)         :: dipolez
!
!  Local variables
!
  integer(i4)      :: i
#ifdef TRACE
  call trace_in('getdipole1D')
#endif
!
!  Loop over atoms and find dipole
!
  dipoley = 0.0_dp
  dipolez = 0.0_dp
  do i = 1,numat
    dipoley = dipoley + occuf(i)*qf(i)*yclat(i)
    dipolez = dipolez + occuf(i)*qf(i)*zclat(i)
  enddo
#ifdef TRACE
  call trace_out('getdipole1D')
#endif
  return
  end
!
  subroutine getdipole0D(dipolex,dipoley,dipolez)
!
!  Calculate the dipole moment in the x, y and z directions 
!  for a cluster.
!
!   8/01 Created
!   2/08 Corrected for partial occupancy
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
!  Julian Gale, CIC, Curtin University, February 2018
!
  use current
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
!
  implicit none
!
!  Passed variables
!
  real(dp)         :: dipolex
  real(dp)         :: dipoley
  real(dp)         :: dipolez
!
!  Local variables
!
  integer(i4)      :: i
#ifdef TRACE
  call trace_in('getdipole0D')
#endif
!
!  Loop over atoms and find dipole
!
  dipolex = 0.0_dp
  dipoley = 0.0_dp
  dipolez = 0.0_dp
  do i = 1,numat
    dipolex = dipolex + occuf(i)*qf(i)*xclat(i)
    dipoley = dipoley + occuf(i)*qf(i)*yclat(i)
    dipolez = dipolez + occuf(i)*qf(i)*zclat(i)
  enddo
#ifdef TRACE
  call trace_out('getdipole0D')
#endif
!
  return
  end
