  subroutine projfork0r(ia,ib,nforkl,nptrfork,nptrforl,d34,dt11,dt21,dt31,dt12,dt22,dt32,dt13,dt23,dt33)
!
!  Calculates the Raman susceptibility contribution from the
!  third atom k to fourth atom l distance for the i-j dynamical 
!  matrix block in a four-body interaction.
!
!   9/13 Created from projfork0
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
  integer(i4),  intent(in)    :: ia
  integer(i4),  intent(in)    :: ib
  integer(i4),  intent(in)    :: nforkl
  integer(i4),  intent(in)    :: nptrfork(*)
  integer(i4),  intent(in)    :: nptrforl(*)
  real(dp),     intent(in)    :: d34(27,*)
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
  integer(i4)                 :: l
  real(dp)                    :: d2dx
  real(dp)                    :: d2dy
  real(dp)                    :: d2dz
#ifdef TRACE
  call trace_in('projfork0r')
#endif
!
!  Loop over three body k atoms
!
  do ii = 1,nforkl
    k = nptrfork(ii)
    l = nptrforl(ii)
!
!  Calculate k-l component - no assumption that d34 is symmetric
!
    d2dx = dt11*d34(1,ii) + dt21*d34(2,ii) + dt31*d34(3,ii) +  &
           dt12*d34(4,ii) + dt22*d34(5,ii) + dt32*d34(6,ii) +  &
           dt13*d34(7,ii) + dt23*d34(8,ii) + dt33*d34(9,ii)
    d2dy = dt11*d34(10,ii) + dt21*d34(11,ii) + dt31*d34(12,ii) +  &
           dt12*d34(13,ii) + dt22*d34(14,ii) + dt32*d34(15,ii) +  &
           dt13*d34(16,ii) + dt23*d34(17,ii) + dt33*d34(18,ii)
    d2dz = dt11*d34(19,ii) + dt21*d34(20,ii) + dt31*d34(21,ii) +  &
           dt12*d34(22,ii) + dt22*d34(23,ii) + dt32*d34(24,ii) +  &
           dt13*d34(25,ii) + dt23*d34(26,ii) + dt33*d34(27,ii)
!
!  Halve term so that it can be applied symmetrically
!
    d2dx = 0.5_dp*d2dx
    d2dy = 0.5_dp*d2dy
    d2dz = 0.5_dp*d2dz
!
!  Add on derivative contributions for k - l
!
    ramanasus(ib,ia,1,k) = ramanasus(ib,ia,1,k) + d2dx
    ramanasus(ib,ia,2,k) = ramanasus(ib,ia,2,k) + d2dy
    ramanasus(ib,ia,3,k) = ramanasus(ib,ia,3,k) + d2dz
    ramanasus(ia,ib,1,k) = ramanasus(ia,ib,1,k) + d2dx
    ramanasus(ia,ib,2,k) = ramanasus(ia,ib,2,k) + d2dy
    ramanasus(ia,ib,3,k) = ramanasus(ia,ib,3,k) + d2dz
!
    ramanasus(ib,ia,1,l) = ramanasus(ib,ia,1,l) - d2dx
    ramanasus(ib,ia,2,l) = ramanasus(ib,ia,2,l) - d2dy
    ramanasus(ib,ia,3,l) = ramanasus(ib,ia,3,l) - d2dz
    ramanasus(ia,ib,1,l) = ramanasus(ia,ib,1,l) - d2dx
    ramanasus(ia,ib,2,l) = ramanasus(ia,ib,2,l) - d2dy
    ramanasus(ia,ib,3,l) = ramanasus(ia,ib,3,l) - d2dz
  enddo
#ifdef TRACE
  call trace_out('projfork0r')
#endif
!
  return
  end
