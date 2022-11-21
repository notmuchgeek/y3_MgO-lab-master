  subroutine projthbk0r(i,j,ia,ib,nmanyk,nptrmanyk,d33,dt11,dt21,dt31,dt12,dt22,dt32,dt13,dt23,dt33)
!
!  Calculates the Raman susceptibility contribution from the
!  third atom k for the i-j dynamical matrix block.
!
!   9/13 Created from projthbk0
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
  use properties,     only : ramanasus
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4),  intent(in)    :: i
  integer(i4),  intent(in)    :: ia
  integer(i4),  intent(in)    :: ib
  integer(i4),  intent(in)    :: j
  integer(i4),  intent(in)    :: nmanyk
  integer(i4),  intent(in)    :: nptrmanyk(*)
  real(dp),     intent(in)    :: d33(54,*)
  real(dp),     intent(in)    :: dt11
  real(dp),     intent(in)    :: dt21
  real(dp),     intent(in)    :: dt31
  real(dp),     intent(in)    :: dt12
  real(dp),     intent(in)    :: dt22
  real(dp),     intent(in)    :: dt32
  real(dp),     intent(in)    :: dt13
  real(dp),     intent(in)    :: dt23
  real(dp),     intent(in)    :: dt33
!
!  Local variables
!
  integer(i4)                 :: ii
  integer(i4)                 :: k
  real(dp)                    :: d2dx
  real(dp)                    :: d2dy
  real(dp)                    :: d2dz
#ifdef TRACE
  call trace_in('projthbk0r')
#endif
!
!  Loop over three body k atoms
!
  do ii = 1,nmanyk
    k = nptrmanyk(ii)
!
!  Calculate i - k component - no assumption that d33 is symmetric
!
    d2dx = dt11*d33(1,ii)  + dt21*d33(2,ii)  + dt31*d33(3,ii) +  &
           dt12*d33(4,ii)  + dt22*d33(5,ii)  + dt32*d33(6,ii) +  &
           dt13*d33(7,ii)  + dt23*d33(8,ii)  + dt33*d33(9,ii)
    d2dy = dt11*d33(10,ii) + dt21*d33(11,ii) + dt31*d33(12,ii) +  &
           dt12*d33(13,ii) + dt22*d33(14,ii) + dt32*d33(15,ii) +  &
           dt13*d33(16,ii) + dt23*d33(17,ii) + dt33*d33(18,ii)
    d2dz = dt11*d33(19,ii) + dt21*d33(20,ii) + dt31*d33(21,ii) +  &
           dt12*d33(22,ii) + dt22*d33(23,ii) + dt32*d33(24,ii) +  &
           dt13*d33(25,ii) + dt23*d33(26,ii) + dt33*d33(27,ii)
!   
!  Halve term so that it can be applied symmetrically
!
    d2dx = 0.5_dp*d2dx
    d2dy = 0.5_dp*d2dy
    d2dz = 0.5_dp*d2dz
!
!  Add on derivative contributions for i - k
!
    ramanasus(ib,ia,1,i) = ramanasus(ib,ia,1,i) + d2dx
    ramanasus(ib,ia,2,i) = ramanasus(ib,ia,2,i) + d2dy
    ramanasus(ib,ia,3,i) = ramanasus(ib,ia,3,i) + d2dz
    ramanasus(ia,ib,1,i) = ramanasus(ia,ib,1,i) + d2dx
    ramanasus(ia,ib,2,i) = ramanasus(ia,ib,2,i) + d2dy
    ramanasus(ia,ib,3,i) = ramanasus(ia,ib,3,i) + d2dz
!
    ramanasus(ib,ia,1,k) = ramanasus(ib,ia,1,k) - d2dx
    ramanasus(ib,ia,2,k) = ramanasus(ib,ia,2,k) - d2dy
    ramanasus(ib,ia,3,k) = ramanasus(ib,ia,3,k) - d2dz
    ramanasus(ia,ib,1,k) = ramanasus(ia,ib,1,k) - d2dx
    ramanasus(ia,ib,2,k) = ramanasus(ia,ib,2,k) - d2dy
    ramanasus(ia,ib,3,k) = ramanasus(ia,ib,3,k) - d2dz
!
!  Calculate j - k component  -  no assumption that d33 is symmetric
!
    d2dx = dt11*d33(28,ii) + dt21*d33(29,ii) + dt31*d33(30,ii) +  &
           dt12*d33(31,ii) + dt22*d33(32,ii) + dt32*d33(33,ii) +  &
           dt13*d33(34,ii) + dt23*d33(35,ii) + dt33*d33(36,ii)
    d2dy = dt11*d33(37,ii) + dt21*d33(38,ii) + dt31*d33(39,ii) +  &
           dt12*d33(40,ii) + dt22*d33(41,ii) + dt32*d33(42,ii) +  &
           dt13*d33(43,ii) + dt23*d33(44,ii) + dt33*d33(45,ii)
    d2dz = dt11*d33(46,ii) + dt21*d33(47,ii) + dt31*d33(48,ii) +  &
           dt12*d33(49,ii) + dt22*d33(50,ii) + dt32*d33(51,ii) +  &
           dt13*d33(52,ii) + dt23*d33(53,ii) + dt33*d33(54,ii)
!   
!  Halve term so that it can be applied symmetrically
!
    d2dx = 0.5_dp*d2dx
    d2dy = 0.5_dp*d2dy
    d2dz = 0.5_dp*d2dz
!
!  Add on derivative contributions for j - k
!
    ramanasus(ib,ia,1,j) = ramanasus(ib,ia,1,j) + d2dx
    ramanasus(ib,ia,2,j) = ramanasus(ib,ia,2,j) + d2dy
    ramanasus(ib,ia,3,j) = ramanasus(ib,ia,3,j) + d2dz
    ramanasus(ia,ib,1,j) = ramanasus(ia,ib,1,j) + d2dx
    ramanasus(ia,ib,2,j) = ramanasus(ia,ib,2,j) + d2dy
    ramanasus(ia,ib,3,j) = ramanasus(ia,ib,3,j) + d2dz
!
    ramanasus(ib,ia,1,k) = ramanasus(ib,ia,1,k) - d2dx
    ramanasus(ib,ia,2,k) = ramanasus(ib,ia,2,k) - d2dy
    ramanasus(ib,ia,3,k) = ramanasus(ib,ia,3,k) - d2dz
    ramanasus(ia,ib,1,k) = ramanasus(ia,ib,1,k) - d2dx
    ramanasus(ia,ib,2,k) = ramanasus(ia,ib,2,k) - d2dy
    ramanasus(ia,ib,3,k) = ramanasus(ia,ib,3,k) - d2dz
  enddo
#ifdef TRACE
  call trace_out('projthbk0r')
#endif
!
  return
  end
