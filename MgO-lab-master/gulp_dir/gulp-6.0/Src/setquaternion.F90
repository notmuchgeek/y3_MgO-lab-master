  subroutine setquaternion(nm,mQ,rQ,drQ,drQ2,lgrad1,lgrad2)
!
!  Subroutine to setup quaternion rotation matrix and optionally its derivatives
!
!  10/19 Created
!  12/19 Second derivatives added
!   3/20 Molecule number added to setquaternion arguments
!   5/20 Standard frame removed
!   7/20 New symmetry adapted algorithm added for quaternions
!   7/20 Check for unphysical e0 added
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
!  Julian Gale, CIC, Curtin University, July 2020
!
  use datatypes
  use molecule,           only : nmolf2a, molQsym
  use symmetry,           only : lsymopt
#ifdef TRACE
  use trace,              only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  logical,       intent(in)   :: lgrad1        ! If true then compute first derivatives of rQ
  logical,       intent(in)   :: lgrad2        ! If true then compute second derivatives of rQ
  integer(i4),   intent(in)   :: nm            ! Molecule number in the full unit cell
  real(dp),      intent(in)   :: mQ(3)         ! Three components of quaternion vectors
  real(dp),      intent(out)  :: rQ(3,3)       ! Rotation matrix
  real(dp),      intent(out)  :: drQ(3,3,3)    ! First derivatives of the rotation matrix if lgrad1 is true
  real(dp),      intent(out)  :: drQ2(3,3,3,3) ! Second derivatives of the rotation matrix if lgrad2 is true
!
!  Local variables
!
  integer(i4)                 :: i
  integer(i4)                 :: ix
  integer(i4)                 :: j
  integer(i4)                 :: jx
  integer(i4)                 :: k
  integer(i4)                 :: nma
  real(dp)                    :: e0
  real(dp)                    :: e1
  real(dp)                    :: e2
  real(dp)                    :: e3
  real(dp)                    :: re0
  real(dp)                    :: re03
  real(dp)                    :: rQtmp(3,3)
#ifdef TRACE
  call trace_in('setquaternion')
#endif
!
!  Setup 4 components of quaternion based on norm = 1
!
  e1 = mQ(1)
  e2 = mQ(2)
  e3 = mQ(3)
  e0 = sqrt(1.0_dp - e1**2 - e2**2 - e3**2)
!
!  Check that e0 is sensible
!
  if (isnan(e0)) then
    call outerror('quaternions have become unphysical - q0 is complex',0_i4)
    call stopnow('setquaternion')
  endif
!
!  Build rotation matrix
!
  rQ(1,1) = e0**2 + e1**2 - e2**2 - e3**2
  rQ(2,1) = 2.0_dp*(e1*e2 - e0*e3)
  rQ(3,1) = 2.0_dp*(e1*e3 + e0*e2)
  rQ(1,2) = 2.0_dp*(e1*e2 + e0*e3)
  rQ(2,2) = e0**2 - e1**2 + e2**2 - e3**2
  rQ(3,2) = 2.0_dp*(e2*e3 - e0*e1)
  rQ(1,3) = 2.0_dp*(e1*e3 - e0*e2)
  rQ(2,3) = 2.0_dp*(e2*e3 + e0*e1)
  rQ(3,3) = e0**2 - e1**2 - e2**2 + e3**2
!
  if (lsymopt) then
!
!  Find molecule in asymmetric unit 
!
    nma = nmolf2a(nm)
!
!  Multiply molQsym by rotation matrix for asymmetric unit to generate rotation matrix for symmetry related molecule
!
    rQtmp(1:3,1:3) = rQ(1:3,1:3)
    do i = 1,3
      do j = 1,3
        rQ(j,i) = 0.0_dp
        do k = 1,3
          rQ(j,i) = rQ(j,i) + molQsym(j,k,nm)*rQtmp(k,i)
        enddo
      enddo
    enddo
  endif
!
  drQ(1:3,1:3,1:3) = 0.0_dp
  drQ2(1:3,1:3,1:3,1:3) = 0.0_dp
!
  if (lgrad1) then
    if (abs(e0).gt.1.0d-10) then
      re0 = 1.0_dp/e0
    else
      call outerror('Quaternion component e0 is too close to zero',0_i4)
      call stopnow('setquaternion')
    endif
!*****************************************
!  First derivatives of rotation matrix  *
!*****************************************
!
!  Multiply by middle index of original vector
!
!  drx/de1
!
    drQ(1,2,1) = 2.0_dp*(e2 + e1*e3*re0)
    drQ(1,3,1) = 2.0_dp*(e3 - e1*e2*re0)
!
!  dry/de1
!
    drQ(1,1,2) = 2.0_dp*(e2 - e1*e3*re0)
    drQ(1,2,2) = - 4.0_dp*e1
    drQ(1,3,2) = 2.0_dp*(-e0 + e1*e1*re0)
!
!  drz/de1
!
    drQ(1,1,3) = 2.0_dp*(e3 + e1*e2*re0)
    drQ(1,2,3) = 2.0_dp*(e0 - e1*e1*re0)
    drQ(1,3,3) = - 4.0_dp*e1
!
!  drx/de2
!
    drQ(2,1,1) = - 4.0_dp*e2
    drQ(2,2,1) = 2.0_dp*(e1 + e2*e3*re0)
    drQ(2,3,1) = 2.0_dp*(e0 - e2*e2*re0)
!
!  dry/de2
!
    drQ(2,1,2) = 2.0_dp*(e1 - e2*e3*re0)
    drQ(2,3,2) = 2.0_dp*(e3 + e1*e2*re0)
!
!  drz/de2
!
    drQ(2,1,3) = 2.0_dp*(-e0 + e2*e2*re0)
    drQ(2,2,3) = 2.0_dp*(e3 - e1*e2*re0)
    drQ(2,3,3) = - 4.0_dp*e2
!
!  drx/de3
!
    drQ(3,1,1) = - 4.0_dp*e3
    drQ(3,2,1) = 2.0_dp*(-e0 + e3*e3*re0)
    drQ(3,3,1) = 2.0_dp*(e1 - e2*e3*re0)
!
!  dry/de3
!
    drQ(3,1,2) = 2.0_dp*(e0 - e3*e3*re0)
    drQ(3,2,2) = - 4.0_dp*e3
    drQ(3,3,2) = 2.0_dp*(e2 + e1*e3*re0)
!
!  drz/de3
!
    drQ(3,1,3) = 2.0_dp*(e1 + e2*e3*re0)
    drQ(3,2,3) = 2.0_dp*(e2 - e1*e3*re0)
!******************************************
!  Second derivatives of rotation matrix  *
!******************************************
    if (lgrad2) then
      re03 = re0**3
!
!  drx/de1
!
!     drQ(1,2,1) = 2.0_dp*(e2 + e1*e3*re0)
      drQ2(1,1,2,1) = 2.0_dp*(e3*re0 + e1*e1*e3*re03)
      drQ2(2,1,2,1) = 2.0_dp*(1.0_dp + e1*e2*e3*re03)
      drQ2(3,1,2,1) = 2.0_dp*(e1*re0 + e1*e3*e3*re03)

!     drQ(1,3,1) = 2.0_dp*(e3 - e1*e2*re0)
      drQ2(1,1,3,1) = - 2.0_dp*(e2*re0 + e1*e1*e2*re03)
      drQ2(2,1,3,1) = - 2.0_dp*(e1*re0 + e1*e2*e2*re03)
      drQ2(3,1,3,1) = 2.0_dp*(1.0_dp - e1*e2*e3*re03)
!
!  dry/de1
!
!     drQ(1,1,2) = 2.0_dp*(e2 - e1*e3*re0)
      drQ2(1,1,1,2) = - 2.0_dp*(e3*re0 + e1*e1*e3*re03)
      drQ2(2,1,1,2) = 2.0_dp*(1.0_dp - e1*e2*e3*re03)
      drQ2(3,1,1,2) = - 2.0_dp*(e1*re0 + e1*e3*e3*re03)

!     drQ(1,2,2) = - 4.0_dp*e1
      drQ2(1,1,2,2) = - 4.0_dp

!     drQ(1,3,2) = 2.0_dp*(-e0 + e1*e1*re0)
      drQ2(1,1,3,2) = 2.0_dp*(e1*re0 + 2.0_dp*e1*re0 + e1*e1*e1*re03)
      drQ2(2,1,3,2) = 2.0_dp*(e2*re0 + e1*e1*e2*re03)
      drQ2(3,1,3,2) = 2.0_dp*(e3*re0 + e1*e1*e3*re03)
!
!  drz/de1
!
!     drQ(1,1,3) = 2.0_dp*(e3 + e1*e2*re0)
      drQ2(1,1,1,3) = 2.0_dp*(e2*re0 + e1*e1*e2*re03)
      drQ2(2,1,1,3) = 2.0_dp*(e1*re0 + e1*e2*e2*re03)
      drQ2(3,1,1,3) = 2.0_dp*(1.0_dp + e1*e2*e3*re03)

!     drQ(1,2,3) = 2.0_dp*(e0 - e1*e1*re0)
      drQ2(1,1,2,3) = - 2.0_dp*(e1*re0 + 2.0_dp*e1*re0 + e1*e1*e1*re03)
      drQ2(2,1,2,3) = - 2.0_dp*(e2*re0 + e1*e1*e2*re03)
      drQ2(3,1,2,3) = - 2.0_dp*(e3*re0 + e1*e1*e3*re03)

!     drQ(1,3,3) = - 4.0_dp*e1
      drQ2(1,1,3,3) = - 4.0_dp
!
!  drx/de2
!
!     drQ(2,1,1) = - 4.0_dp*e2
      drQ2(2,2,1,1) = - 4.0_dp

!     drQ(2,2,1) = 2.0_dp*(e1 + e2*e3*re0)
      drQ2(1,2,2,1) = 2.0_dp*(1.0_dp + e1*e2*e3*re03)
      drQ2(2,2,2,1) = 2.0_dp*(e3*re0 + e2*e2*e3*re03)
      drQ2(3,2,2,1) = 2.0_dp*(e2*re0 + e2*e3*e3*re03)

!     drQ(2,3,1) = 2.0_dp*(e0 - e2*e2*re0)
      drQ2(1,2,3,1) = - 2.0_dp*(e1*re0 + e2*e2*e1*re03)
      drQ2(2,2,3,1) = - 2.0_dp*(e2*re0 + 2.0_dp*e2*re0 + e2*e2*e2*re03)
      drQ2(3,2,3,1) = - 2.0_dp*(e3*re0 + e2*e2*e3*re03)
!
!  dry/de2
!
!     drQ(2,1,2) = 2.0_dp*(e1 - e2*e3*re0)
      drQ2(1,2,1,2) = 2.0_dp*(1.0_dp - e1*e2*e3*re03)
      drQ2(2,2,1,2) = - 2.0_dp*(e3*re0 + e2*e2*e3*re03)
      drQ2(3,2,1,2) = - 2.0_dp*(e2*re0 + e2*e3*e3*re03)

!     drQ(2,3,2) = 2.0_dp*(e3 + e1*e2*re0)
      drQ2(1,2,3,2) = 2.0_dp*(e2*re0 + e1*e1*e2*re03)
      drQ2(2,2,3,2) = 2.0_dp*(e1*re0 + e1*e2*e2*re03)
      drQ2(3,2,3,2) = 2.0_dp*(1.0_dp + e1*e2*e3*re03)
!
!  drz/de2
!
!     drQ(2,1,3) = 2.0_dp*(-e0 + e2*e2*re0)
      drQ2(1,2,1,3) = 2.0_dp*(e1*re0 + e1*e2*e2*re03)
      drQ2(2,2,1,3) = 2.0_dp*(e2*re0 + 2.0_dp*e2*re0 + e2*e2*e2*re03)
      drQ2(3,2,1,3) = 2.0_dp*(e3*re0 + e2*e2*e3*re03)

!     drQ(2,2,3) = 2.0_dp*(e3 - e1*e2*re0)
      drQ2(1,2,2,3) = - 2.0_dp*(e2*re0 + e1*e1*e2*re03)
      drQ2(2,2,2,3) = - 2.0_dp*(e1*re0 + e1*e2*e2*re03)
      drQ2(3,2,2,3) = 2.0_dp*(1.0_dp - e1*e2*e3*re03)

!     drQ(2,3,3) = - 4.0_dp*e2
      drQ2(2,2,3,3) = - 4.0_dp
!
!  drx/de3
!
!     drQ(3,1,1) = - 4.0_dp*e3
      drQ2(3,3,1,1) = - 4.0_dp

!     drQ(3,2,1) = 2.0_dp*(-e0 + e3*e3*re0)
      drQ2(1,3,2,1) = 2.0_dp*(e1*re0 + e1*e3*e3*re03)
      drQ2(2,3,2,1) = 2.0_dp*(e2*re0 + e2*e3*e3*re03)
      drQ2(3,3,2,1) = 2.0_dp*(e3*re0 + 2.0_dp*e3*re0 + e3*e3*e3*re03)

!     drQ(3,3,1) = 2.0_dp*(e1 - e2*e3*re0)
      drQ2(1,3,3,1) = 2.0_dp*(1.0_dp - e1*e2*e3*re03)
      drQ2(2,3,3,1) = - 2.0_dp*(e3*re0 + e2*e2*e3*re03)
      drQ2(3,3,3,1) = - 2.0_dp*(e2*re0 + e2*e3*e3*re03)
!
!  dry/de3
!
!     drQ(3,1,2) = 2.0_dp*(e0 - e3*e3*re0)
      drQ2(1,3,1,2) = - 2.0_dp*(e1*re0 + e1*e3*e3*re03)
      drQ2(2,3,1,2) = - 2.0_dp*(e2*re0 + e2*e3*e3*re03)
      drQ2(3,3,1,2) = - 2.0_dp*(e3*re0 + 2.0_dp*e3*re0 + e3*e3*e3*re03)

!     drQ(3,2,2) = - 4.0_dp*e3
      drQ2(3,3,2,2) = - 4.0_dp

!     drQ(3,3,2) = 2.0_dp*(e2 + e1*e3*re0)
      drQ2(1,3,3,2) = 2.0_dp*(e3*re0 + e1*e1*e3*re03)
      drQ2(2,3,3,2) = 2.0_dp*(1.0_dp + e1*e2*e3*re03)
      drQ2(3,3,3,2) = 2.0_dp*(e1*re0 + e1*e3*e3*re03)
!
!  drz/de3
!
!     drQ(3,1,3) = 2.0_dp*(e1 + e2*e3*re0)
      drQ2(1,3,1,3) = 2.0_dp*(1.0_dp + e1*e2*e3*re03)
      drQ2(2,3,1,3) = 2.0_dp*(e3*re0 + e2*e2*e3*re03)
      drQ2(3,3,1,3) = 2.0_dp*(e2*re0 + e2*e3*e3*re03)

!     drQ(3,2,3) = 2.0_dp*(e2 - e1*e3*re0)
      drQ2(1,3,2,3) = - 2.0_dp*(e3*re0 + e1*e1*e3*re03)
      drQ2(2,3,2,3) = 2.0_dp*(1.0_dp - e1*e2*e3*re03)
      drQ2(3,3,2,3) = - 2.0_dp*(e1*re0 + e1*e3*e3*re03)
    endif
!
    if (lsymopt) then
!
!  Multiply drQ by rotation matrix for asymmetric unit to generate rotation matrix for symmetry related molecule
!
      do ix = 1,3
        rQtmp(1:3,1:3) = drQ(ix,1:3,1:3)
        do i = 1,3
          do j = 1,3
            drQ(ix,j,i) = 0.0_dp
            do k = 1,3
              drQ(ix,j,i) = drQ(ix,j,i) + molQsym(j,k,nm)*rQtmp(k,i)
            enddo
          enddo
        enddo
      enddo
!
      if (lgrad2) then
        do ix = 1,3
          do jx = 1,3
            rQtmp(1:3,1:3) = drQ2(jx,ix,1:3,1:3)
            do i = 1,3
              do j = 1,3
                drQ2(jx,ix,j,i) = 0.0_dp
                do k = 1,3
                  drQ2(jx,ix,j,i) = drQ2(jx,ix,j,i) + molQsym(j,k,nm)*rQtmp(k,i)
                enddo
              enddo
            enddo
          enddo
        enddo
      endif
    endif
  endif
!
#ifdef TRACE
  call trace_out('setquaternion')
#endif
!
  return
  end
