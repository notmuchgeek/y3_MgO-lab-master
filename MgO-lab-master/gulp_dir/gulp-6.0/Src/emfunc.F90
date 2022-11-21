  subroutine emfunc(lm,u,x,signx,y,z,a,emf,d1emf,d2emf,d3emf,lgrad1,lgrad2,lgrad3)
!
!  Calculates the Euler-Maclaurin integrals required by the 1-D Coulomb sum.
!  NB: This version doesn't have explicit strain derivatives and so doesn't
!      need to support rigid molecules.
!
!   9/01 Created
!  10/01 Derivatives added
!  10/01 Derivative algorithm simplified using recursive relationship of
!        W terms w.r.t. differentation
!  12/01 First strain derivative added
!  12/01 Second derivatives added
!   5/02 Second strain derivatives added
!   5/02 Third derivatives added
!   7/05 Style updated
!   2/18 Trace added
!  12/18 Finite strain derivatives added
!   3/20 Modified to remove strains following creation of emfuncs
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
!  Julian Gale, CIC, Curtin University, March 2020
!
  use datatypes
#ifdef TRACE
  use trace,       only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)  :: lm
  real(dp),    intent(in)  :: a
  real(dp),    intent(out) :: d1emf(3)
  real(dp),    intent(out) :: d2emf(6)
  real(dp),    intent(out) :: d3emf(10)
  real(dp),    intent(out) :: emf
  real(dp),    intent(in)  :: signx
  real(dp),    intent(in)  :: u
  real(dp),    intent(in)  :: x
  real(dp),    intent(in)  :: y
  real(dp),    intent(in)  :: z
  logical,     intent(in)  :: lgrad1
  logical,     intent(in)  :: lgrad2
  logical,     intent(in)  :: lgrad3
!
!  Local variables
!
  integer(i4)              :: i
  integer(i4)              :: ii
  integer(i4)              :: maxlmorder
  integer(i4)              :: maxlmorderd1
  integer(i4)              :: maxlmorderd2
  real(dp)                 :: alpha
  real(dp)                 :: ap
  real(dp)                 :: g
  real(dp)                 :: d1g(3)
  real(dp)                 :: d2g(6)
  real(dp)                 :: ecoeff(5)
  real(dp)                 :: gtrm1
  real(dp)                 :: gtrm2
  real(dp)                 :: raptrm
  real(dp)                 :: ux
  real(dp)                 :: w(0:11)
  real(dp)                 :: dw(3,0:11)
  real(dp)                 :: d2w(6,0:11)
  real(dp)                 :: d3w(10,0:11)
#ifdef TRACE
  call trace_in('emfunc')
#endif
!
  if (lm.gt.5) then
    call outerror(' Order of Euler-MacLaurin expansion exceeds 5!',0_i4)
    call stopnow('emfunc')
  endif
!
  ecoeff(1) = -0.5_dp/12.0_dp
  ecoeff(2) =  0.875_dp/720.0_dp
  ecoeff(3) = -0.96875_dp/30240.0_dp
  ecoeff(4) =  0.9921875_dp/1209600.0_dp
  ecoeff(5) = -0.998046875_dp/47900160.0_dp
!
  alpha = y*y + z*z
  ux = u + x
!
!  Set order of W to calculate to - increase by 1 for each level of derivatives
!
  maxlmorder = 2*lm - 1
  maxlmorderd1 = 2*lm - 1
  maxlmorderd2 = 2*lm - 1
  if (lgrad1) then
    maxlmorder = maxlmorder + 1
    if (lgrad2) then
      maxlmorder = maxlmorder + 1
      maxlmorderd1 = maxlmorderd1 + 1
      if (lgrad3) then
        maxlmorder = maxlmorder + 1
        maxlmorderd1 = maxlmorderd1 + 1
        maxlmorderd2 = maxlmorderd2 + 1
      endif
    endif
  endif
!
  g = (ux*ux+alpha)
  raptrm = 1.0_dp/g
  w(0) = sqrt(raptrm)
  w(1) = -ux*w(0)*raptrm
  do i = 2,maxlmorder
    w(i) = (-dble(2*i-1)*ux*w(i-1) - dble((i-1)**2)*w(i-2))*raptrm
  enddo
!
!  Function
!
  emf = 0.0_dp
  ap = a
  do i = 1,lm
    emf = emf - ecoeff(i)*ap*w(2*i-1)
    ap = ap*a*a
  enddo
!
!  Derivatives
!
  if (lgrad1) then
!
!  First derivatives of g
!
    d1g(1) = 2.0_dp*ux*signx
    d1g(2) = 2.0_dp*y
    d1g(3) = 2.0_dp*z
!
!  Cartesian derivatives of W
!
    do i = 0,maxlmorderd1
      dw(1,i) = w(i+1)*signx
    enddo
    dw(2,0) = - w(0)*raptrm*y
    dw(3,0) = - w(0)*raptrm*z
    dw(2,1) = -(2.0_dp*w(1)*y + ux*dw(2,0))*raptrm
    dw(3,1) = -(2.0_dp*w(1)*z + ux*dw(3,0))*raptrm
    do i = 2,maxlmorderd1
      dw(2,i) = - (2.0_dp*w(i)*y + (dble(2*i-1)*ux*dw(2,i-1) + dble((i-1)**2)*dw(2,i-2)))*raptrm
      dw(3,i) = - (2.0_dp*w(i)*z + (dble(2*i-1)*ux*dw(3,i-1) + dble((i-1)**2)*dw(3,i-2)))*raptrm
    enddo
!
    if (lgrad2) then
!
!  Second derivatives of g
!
      d2g(1) = 2.0_dp
      d2g(2) = 2.0_dp
      d2g(3) = 2.0_dp
      d2g(4) = 0.0_dp
      d2g(5) = 0.0_dp
      d2g(6) = 0.0_dp
!
!  Cartesian - cartesian derivatives of W
!
      d2w(1,0) = w(2)
      d2w(2,0) = w(0)*raptrm*(3.0_dp*y*y*raptrm - 1.0_dp)
      d2w(3,0) = w(0)*raptrm*(3.0_dp*z*z*raptrm - 1.0_dp)
      d2w(4,0) = w(0)*raptrm*3.0_dp*y*z*raptrm
      d2w(5,0) = - w(1)*signx*raptrm*z
      d2w(6,0) = - w(1)*signx*raptrm*y
      d2w(1,1) = w(3)
      d2w(2,1) = - raptrm*(4.0_dp*raptrm*y*dw(2,0) + ux*d2w(2,0) + 2.0_dp*w(1))
      d2w(3,1) = - raptrm*(4.0_dp*raptrm*z*dw(3,0) + ux*d2w(3,0) + 2.0_dp*w(1))
      d2w(4,1) = - raptrm*(4.0_dp*raptrm*z*dw(2,0) + ux*d2w(4,0))
      d2w(5,1) = -(2.0_dp*w(2)*z*signx + signx*dw(3,0) + ux*d2w(5,0))*raptrm
      d2w(6,1) = -(2.0_dp*w(2)*y*signx + signx*dw(2,0) + ux*d2w(6,0))*raptrm
      do i = 2,maxlmorderd2
        d2w(1,i) = w(i+2)
        d2w(2,i) = - 4.0_dp*raptrm*y*dw(2,i) - (dble(2*i-1)*ux*d2w(2,i-1) + dble((i-1)**2)*d2w(2,i-2) + 2.0*w(i))*raptrm
        d2w(3,i) = - 4.0_dp*raptrm*z*dw(3,i) - (dble(2*i-1)*ux*d2w(3,i-1) + dble((i-1)**2)*d2w(3,i-2) + 2.0*w(i))*raptrm
        d2w(4,i) = - 4.0_dp*raptrm*z*dw(2,i) - (dble(2*i-1)*ux*d2w(4,i-1) + dble((i-1)**2)*d2w(4,i-2))*raptrm
        d2w(5,i) = - (2.0_dp*w(i+1)*z*signx + (dble(2*i-1)*ux*d2w(5,i-1) +  &
          signx*dble(2*i-1)*dw(3,i-1) + dble((i-1)**2)*d2w(5,i-2)))*raptrm
        d2w(6,i) = - (2.0_dp*w(i+1)*y*signx + (dble(2*i-1)*ux*d2w(6,i-1) +  &
          signx*dble(2*i-1)*dw(2,i-1) + dble((i-1)**2)*d2w(6,i-2)))*raptrm
      enddo
      if (lgrad3) then
        gtrm1 = 0.75_dp*raptrm*raptrm*w(0)
        gtrm2 = - 2.5_dp*gtrm1*raptrm
!
!  Cartesian - cartesian - cartesian derivatives of W
!
        d3w(1,0) = gtrm2*d1g(1)*d1g(1)*d1g(1) + gtrm1*(3.0_dp*d2g(1)*d1g(1))
        d3w(2,0) = gtrm2*d1g(2)*d1g(1)*d1g(1) + gtrm1*(2.0_dp*d2g(6)*d1g(1) + d2g(1)*d1g(2))
        d3w(3,0) = gtrm2*d1g(3)*d1g(1)*d1g(1) + gtrm1*(2.0_dp*d2g(5)*d1g(1) + d2g(1)*d1g(3))
        d3w(4,0) = gtrm2*d1g(2)*d1g(2)*d1g(1) + gtrm1*(2.0_dp*d2g(6)*d1g(2) + d2g(2)*d1g(1))
        d3w(5,0) = gtrm2*d1g(3)*d1g(2)*d1g(1) + gtrm1*(d2g(4)*d1g(1) + d2g(5)*d1g(2) + d2g(6)*d1g(3))
        d3w(6,0) = gtrm2*d1g(3)*d1g(3)*d1g(1) + gtrm1*(2.0_dp*d2g(5)*d1g(1) + d2g(3)*d1g(1))
        d3w(7,0) = gtrm2*d1g(2)*d1g(2)*d1g(2) + gtrm1*(3.0_dp*d2g(2)*d1g(2))
        d3w(8,0) = gtrm2*d1g(3)*d1g(2)*d1g(2) + gtrm1*(2.0_dp*d2g(4)*d1g(2) + d2g(2)*d1g(3))
        d3w(9,0) = gtrm2*d1g(3)*d1g(3)*d1g(2) + gtrm1*(2.0_dp*d2g(4)*d1g(3) + d2g(3)*d1g(2))
        d3w(10,0)= gtrm2*d1g(3)*d1g(3)*d1g(3) + gtrm1*(3.0_dp*d2g(3)*d1g(3))
!
        d3w(1,1) = -6.0_dp*ux*(dw(1,0)*dw(1,0)*dw(1,0) + w(0)*(3.0_dp*d2w(1,0)*dw(1,0) +  &
          0.5_dp*w(0)*d3w(1,0))) - 3.0_dp*signx*w(0)*(6.0_dp*dw(1,0)*dw(1,0) + w(0)*(3.0_dp*d2w(1,0)))
        d3w(2,1) = -6.0_dp*ux*(dw(2,0)*dw(1,0)*dw(1,0) + w(0)*(2.0_dp*d2w(6,0)*dw(1,0) + d2w(1,0)*dw(2,0) + &
          0.5_dp*w(0)*d3w(2,0))) - 6.0_dp*signx*w(0)*(2.0_dp*dw(2,0)*dw(1,0) + w(0)*d2w(6,0))
        d3w(3,1) = -6.0_dp*ux*(dw(3,0)*dw(1,0)*dw(1,0) + w(0)*(2.0_dp*d2w(5,0)*dw(1,0) + d2w(1,0)*dw(3,0) + &
          0.5_dp*w(0)*d3w(3,0))) - 6.0_dp*signx*w(0)*(2.0_dp*dw(3,0)*dw(1,0) + w(0)*d2w(5,0))
        d3w(4,1) = -6.0_dp*ux*(dw(2,0)*dw(2,0)*dw(1,0) + w(0)*(2.0_dp*d2w(6,0)*dw(2,0) + d2w(2,0)*dw(1,0) + &
          0.5_dp*w(0)*d3w(4,0))) - 6.0_dp*signx*w(0)*(dw(2,0)*dw(2,0) + 0.5_dp*w(0)*d2w(2,0))
        d3w(5,1) = -6.0_dp*ux*(dw(3,0)*dw(2,0)*dw(1,0) + w(0)*(d2w(4,0)*dw(1,0) + d2w(5,0)*dw(2,0) + &
          d2w(6,0)*dw(3,0) + 0.5_dp*w(0)*d3w(5,0))) - 6.0_dp*signx*w(0)*(dw(3,0)*dw(2,0) + 0.5_dp*w(0)*d2w(4,0))
        d3w(6,1) = -6.0_dp*ux*(dw(3,0)*dw(3,0)*dw(1,0) + w(0)*(2.0_dp*d2w(5,0)*dw(3,0) + d2w(3,0)*dw(1,0) + &
          0.5_dp*w(0)*d3w(6,0))) - 6.0_dp*signx*w(0)*(dw(3,0)*dw(3,0) + 0.5_dp*w(0)*d2w(3,0))
        d3w(7,1) = -6.0_dp*ux*(dw(2,0)*dw(2,0)*dw(2,0) + w(0)*(3.0_dp*d2w(2,0)*dw(2,0) + 0.5_dp*w(0)*d3w(7,0)))
        d3w(8,1) = -6.0_dp*ux*(dw(3,0)*dw(2,0)*dw(2,0) + w(0)*(2.0_dp*d2w(5,0)*dw(2,0) + d2w(2,0)*dw(3,0) + &
          0.5_dp*w(0)*d3w(8,0)))
        d3w(9,1) = -6.0_dp*ux*(dw(3,0)*dw(3,0)*dw(2,0) + w(0)*(2.0_dp*d2w(5,0)*dw(3,0) + d2w(3,0)*dw(2,0) + &
          0.5_dp*w(0)*d3w(9,0)))
        d3w(10,1) = -6.0_dp*ux*(dw(3,0)*dw(3,0)*dw(3,0) + w(0)*(3.0_dp*d2w(3,0)*dw(3,0) + 0.5_dp*w(0)*d3w(10,0)))
!
        do i = 2,2*lm-1
          d3w(1,i) = - raptrm*(3.0_dp*(d2w(1,i)*d1g(1) + dw(1,i)*d2g(1)) + &
            dble(2*i-1)*(signx*(3.0_dp*d2w(1,i-1)) + ux*d3w(1,i-1)) + dble((i-1)**2)*d3w(1,i-2))
          d3w(2,i) = - raptrm*( &
            2.0_dp*d2w(6,i)*d1g(1) + d2w(1,i)*d1g(2) + 2.0_dp*dw(1,i)*d2g(6) + dw(2,i)*d2g(1) + &
            dble(2*i-1)*(signx*(2.0_dp*d2w(6,i-1)) + ux*d3w(2,i-1)) + dble((i-1)**2)*d3w(2,i-2))
          d3w(3,i) = - raptrm*(2.0_dp*d2w(5,i)*d1g(1) + d2w(1,i)*d1g(3) + &
            2.0_dp*dw(1,i)*d2g(5) + dw(3,i)*d2g(1) + dble(2*i-1)*(signx*(2.0_dp*d2w(5,i-1)) +  &
            ux*d3w(3,i-1)) + dble((i-1)**2)*d3w(3,i-2))
          d3w(4,i) = - raptrm*(2.0_dp*d2w(6,i)*d1g(2) + d2w(2,i)*d1g(1) + 2.0_dp*dw(2,i)*d2g(6) + dw(1,i)*d2g(2) + &
            dble(2*i-1)*(signx*(d2w(2,i-1)) + ux*d3w(4,i-1)) + dble((i-1)**2)*d3w(4,i-2))
          d3w(5,i) = - raptrm*(d2w(4,i)*d1g(1) + d2w(5,i)*d1g(2) + d2w(6,i)*d1g(3) + dw(1,i)*d2g(4) +  &
            dw(2,i)*d2g(5) + dw(3,i)*d2g(6) + dble(2*i-1)*(signx*(d2w(4,i-1)) + ux*d3w(5,i-1)) + dble((i-1)**2)*d3w(5,i-2))
          d3w(6,i) = - raptrm*(2.0_dp*d2w(5,i)*d1g(3) + d2w(3,i)*d1g(1) + 2.0_dp*dw(3,i)*d2g(5) + dw(1,i)*d2g(3) + &
            dble(2*i-1)*(signx*(d2w(3,i-1)) + ux*d3w(6,i-1)) + dble((i-1)**2)*d3w(6,i-2))
          d3w(7,i) = - raptrm*(3.0_dp*(d2w(2,i)*d1g(2) + dw(2,i)*d2g(2)) + &
            dble(2*i-1)*(ux*d3w(7,i-1)) + dble((i-1)**2)*d3w(7,i-2))
          d3w(8,i) = - raptrm*(2.0_dp*d2w(4,i)*d1g(2) + d2w(2,i)*d1g(3) + &
            2.0_dp*dw(2,i)*d2g(4) + dw(3,i)*d2g(2) + dble(2*i-1)*(ux*d3w(8,i-1)) + dble((i-1)**2)*d3w(8,i-2))
          d3w(9,i) = - raptrm*(2.0_dp*d2w(4,i)*d1g(3) + d2w(3,i)*d1g(2) + 2.0_dp*dw(3,i)*d2g(4) + dw(2,i)*d2g(3) + &
            dble(2*i-1)*(ux*d3w(9,i-1)) + dble((i-1)**2)*d3w(9,i-2))
          d3w(10,i) = - raptrm*(3.0_dp*(d2w(3,i)*d1g(3) + dw(3,i)*d2g(3)) + &
            dble(2*i-1)*(ux*d3w(10,i-1)) + dble((i-1)**2)*d3w(10,i-2))
        enddo
      endif
    endif
!
!  Initialise total derivatives to zero
!
    d1emf(1:3) = 0.0_dp
    if (lgrad2) then
      d2emf(1:6) = 0.0_dp
      if (lgrad3) then
        d3emf(1:10) = 0.0_dp
      endif
    endif
    ap = a
    do i = 1,lm
!
!  First derivatives : Cartesian
!
      d1emf(1) = d1emf(1) - ecoeff(i)*ap*dw(1,2*i-1)
      d1emf(2) = d1emf(2) - ecoeff(i)*ap*dw(2,2*i-1)
      d1emf(3) = d1emf(3) - ecoeff(i)*ap*dw(3,2*i-1)
      if (lgrad2) then
!
!  Second derivatives : cartesian - cartesian
!
        d2emf(1) = d2emf(1) - ecoeff(i)*ap*d2w(1,2*i-1)
        d2emf(2) = d2emf(2) - ecoeff(i)*ap*d2w(2,2*i-1)
        d2emf(3) = d2emf(3) - ecoeff(i)*ap*d2w(3,2*i-1)
        d2emf(4) = d2emf(4) - ecoeff(i)*ap*d2w(4,2*i-1)
        d2emf(5) = d2emf(5) - ecoeff(i)*ap*d2w(5,2*i-1)
        d2emf(6) = d2emf(6) - ecoeff(i)*ap*d2w(6,2*i-1)
        if (lgrad3) then
!
!  Third derivatives : cartesian - cartesian - cartesian
!
          do ii = 1,10
            d3emf(ii) = d3emf(ii) - ecoeff(i)*ap*d3w(ii,2*i-1)
          enddo
        endif
      endif
      ap = ap*a*a
    enddo
  endif
#ifdef TRACE
  call trace_out('emfunc')
#endif
!
  return
  end
!
  subroutine emfuncs(lm,u,x,signx,y,z,xcom,a,emf,d1emf,d2emf,d1emfs, &
                     d2emfs,d2emfm,d3emf,d3emfm,lgrad1,lgrad2,lgrad3)
!
!  Calculates the Euler-Maclaurin integrals required by the 1-D Coulomb sum.
!  NB: This version has explicit strain derivatives
!
!   3/20 Created from emfunc
!   3/20 Centre of mass corrections added for rigid molecules
!   4/20 Corrections to handle finite strain
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
!  Julian Gale, CIC, Curtin University, April 2020
!
  use datatypes
  use current,     only : strain
  use derivatives, only : lfinitestrain
#ifdef TRACE
  use trace,       only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)  :: lm
  real(dp),    intent(in)  :: a
  real(dp),    intent(out) :: d1emf(3)
  real(dp),    intent(out) :: d2emf(6)
  real(dp),    intent(out) :: d2emfm(3)
  real(dp),    intent(out) :: d3emf(10)
  real(dp),    intent(out) :: d3emfm(6)
  real(dp),    intent(out) :: d1emfs
  real(dp),    intent(out) :: d2emfs
  real(dp),    intent(out) :: emf
  real(dp),    intent(in)  :: signx
  real(dp),    intent(in)  :: u
  real(dp),    intent(in)  :: x
  real(dp),    intent(in)  :: y
  real(dp),    intent(in)  :: z
  real(dp),    intent(in)  :: xcom
  logical,     intent(in)  :: lgrad1
  logical,     intent(in)  :: lgrad2
  logical,     intent(in)  :: lgrad3
!
!  Local variables
!
  integer(i4)              :: i
  integer(i4)              :: ii
  integer(i4)              :: maxlmorder
  integer(i4)              :: maxlmorderd1
  integer(i4)              :: maxlmorderd2
  real(dp)                 :: alpha
  real(dp)                 :: ap
  real(dp)                 :: dwde
  real(dp)                 :: g
  real(dp)                 :: d1g(3)
  real(dp)                 :: d1gs
  real(dp)                 :: d2g(6)
  real(dp)                 :: d2gm(3)
  real(dp)                 :: d2gs
  real(dp)                 :: ecoeff(5)
  real(dp)                 :: gtrm1
  real(dp)                 :: gtrm2
  real(dp)                 :: raptrm
  real(dp)                 :: trm1s
  real(dp)                 :: dads
  real(dp)                 :: d2ads2
  real(dp)                 :: ux
  real(dp)                 :: duxds
  real(dp)                 :: d2uxds2
  real(dp)                 :: d2uxdsdx
  real(dp)                 :: w(0:11)
  real(dp)                 :: dw(3,0:11)
  real(dp)                 :: dws(0:11)
  real(dp)                 :: d2w(6,0:11)
  real(dp)                 :: d2wm(3,0:11)
  real(dp)                 :: d2ws(0:11)
  real(dp)                 :: d3w(10,0:11)
  real(dp)                 :: d3ws(6,0:11)
#ifdef TRACE
  call trace_in('emfuncs')
#endif
!
  if (lm.gt.5) then
    call outerror(' Order of Euler-MacLaurin expansion exceeds 5!',0_i4)
    call stopnow('emfuncs')
  endif
!
  ecoeff(1) = -0.5_dp/12.0_dp
  ecoeff(2) =  0.875_dp/720.0_dp
  ecoeff(3) = -0.96875_dp/30240.0_dp
  ecoeff(4) =  0.9921875_dp/1209600.0_dp
  ecoeff(5) = -0.998046875_dp/47900160.0_dp
!
  alpha = y*y + z*z
  ux = u + x
!
!  Set order of W to calculate to - increase by 1 for each level of derivatives
!
  maxlmorder = 2*lm - 1
  maxlmorderd1 = 2*lm - 1
  maxlmorderd2 = 2*lm - 1
  if (lgrad1) then
    maxlmorder = maxlmorder + 1
    if (lgrad2) then
      maxlmorder = maxlmorder + 1
      maxlmorderd1 = maxlmorderd1 + 1
      if (lgrad3) then
        maxlmorder = maxlmorder + 1
        maxlmorderd1 = maxlmorderd1 + 1
        maxlmorderd2 = maxlmorderd2 + 1
      endif
    endif
  endif
!
  g = (ux*ux + alpha)
  raptrm = 1.0_dp/g
  w(0) = sqrt(raptrm)
  w(1) = -ux*w(0)*raptrm
  do i = 2,maxlmorder
    w(i) = (-dble(2*i-1)*ux*w(i-1) - dble((i-1)**2)*w(i-2))*raptrm
  enddo
!
!  Function
!
  emf = 0.0_dp
  ap = a
  do i = 1,lm
    emf = emf - ecoeff(i)*ap*w(2*i-1)
    ap = ap*a*a
  enddo
!
!  Derivatives
!
  if (lgrad1) then
    if (lfinitestrain) then
      dads  = 1.0_dp/(1.0_dp + strain(1))  ! dads is the derivative of a with respect to strain, divided by a
      duxds = (ux - xcom)/(1.0_dp + strain(1))
    else
      dads  = 1.0_dp
      duxds = ux - xcom
    endif
!
!  First derivatives of g
!
    d1g(1) = 2.0_dp*ux*signx
    d1g(2) = 2.0_dp*y
    d1g(3) = 2.0_dp*z
    d1gs   = 2.0_dp*ux*duxds
!
!  Cartesian derivatives of W
!
    do i = 0,maxlmorderd1
      dw(1,i) = w(i+1)*signx
    enddo
    dw(2,0) = - w(0)*raptrm*y
    dw(3,0) = - w(0)*raptrm*z
    dw(2,1) = -(2.0_dp*w(1)*y + ux*dw(2,0))*raptrm
    dw(3,1) = -(2.0_dp*w(1)*z + ux*dw(3,0))*raptrm
    do i = 2,maxlmorderd1
      dw(2,i) = - (2.0_dp*w(i)*y + (dble(2*i-1)*ux*dw(2,i-1) + dble((i-1)**2)*dw(2,i-2)))*raptrm
      dw(3,i) = - (2.0_dp*w(i)*z + (dble(2*i-1)*ux*dw(3,i-1) + dble((i-1)**2)*dw(3,i-2)))*raptrm
    enddo
!
!  Strain derivatives of W
!
    do i = 0,maxlmorderd1
      dws(i)  = w(i+1)*duxds
    enddo
!
    if (lgrad2) then
      if (lfinitestrain) then
        d2ads2  = 0.0_dp
        d2uxds2 = 0.0_dp
        d2uxdsdx = signx/(1.0_dp + strain(1))
      else
        d2ads2  = 1.0_dp     ! Second derivative of a w.r.t. strain divided by a
        d2uxds2 = duxds
        d2uxdsdx = signx
      endif
!
!  Second derivatives of g
!
      d2g(1) = 2.0_dp
      d2g(2) = 2.0_dp
      d2g(3) = 2.0_dp
      d2g(4) = 0.0_dp
      d2g(5) = 0.0_dp
      d2g(6) = 0.0_dp
!
      d2gs    = 2.0_dp*(duxds*duxds + ux*d2uxds2)
      d2gm(1) = 2.0_dp*duxds*signx + 2.0_dp*ux*d2uxdsdx
      d2gm(2) = 0.0_dp
      d2gm(3) = 0.0_dp
!
!  Strain - strain derivatives of W
!
      do i = 0,maxlmorderd2
        d2ws(i) = dws(i+1)*duxds + w(i+1)*d2uxds2
      enddo
!
!  Cartesian - strain derivatives of W
!
      do i = 0,maxlmorderd2
        d2wm(1,i) = dw(1,i+1)*duxds + w(i+1)*d2uxdsdx
        d2wm(2,i) = dw(2,i+1)*duxds
        d2wm(3,i) = dw(3,i+1)*duxds
      enddo
!
!  Cartesian - cartesian derivatives of W
!
      d2w(1,0) = w(2)
      d2w(2,0) = w(0)*raptrm*(3.0_dp*y*y*raptrm - 1.0_dp)
      d2w(3,0) = w(0)*raptrm*(3.0_dp*z*z*raptrm - 1.0_dp)
      d2w(4,0) = w(0)*raptrm*3.0_dp*y*z*raptrm
      d2w(5,0) = - w(1)*signx*raptrm*z
      d2w(6,0) = - w(1)*signx*raptrm*y
      d2w(1,1) = w(3)
      d2w(2,1) = - raptrm*(4.0_dp*raptrm*y*dw(2,0) + ux*d2w(2,0) + 2.0_dp*w(1))
      d2w(3,1) = - raptrm*(4.0_dp*raptrm*z*dw(3,0) + ux*d2w(3,0) + 2.0_dp*w(1))
      d2w(4,1) = - raptrm*(4.0_dp*raptrm*z*dw(2,0) + ux*d2w(4,0))
      d2w(5,1) = -(2.0_dp*w(2)*z*signx + signx*dw(3,0) + ux*d2w(5,0))*raptrm
      d2w(6,1) = -(2.0_dp*w(2)*y*signx + signx*dw(2,0) + ux*d2w(6,0))*raptrm
      do i = 2,maxlmorderd2
        d2w(1,i) = w(i+2)
        d2w(2,i) = - 4.0_dp*raptrm*y*dw(2,i) - (dble(2*i-1)*ux*d2w(2,i-1) + dble((i-1)**2)*d2w(2,i-2) + 2.0*w(i))*raptrm
        d2w(3,i) = - 4.0_dp*raptrm*z*dw(3,i) - (dble(2*i-1)*ux*d2w(3,i-1) + dble((i-1)**2)*d2w(3,i-2) + 2.0*w(i))*raptrm
        d2w(4,i) = - 4.0_dp*raptrm*z*dw(2,i) - (dble(2*i-1)*ux*d2w(4,i-1) + dble((i-1)**2)*d2w(4,i-2))*raptrm
        d2w(5,i) = - (2.0_dp*w(i+1)*z*signx + (dble(2*i-1)*ux*d2w(5,i-1) +  &
          signx*dble(2*i-1)*dw(3,i-1) + dble((i-1)**2)*d2w(5,i-2)))*raptrm
        d2w(6,i) = - (2.0_dp*w(i+1)*y*signx + (dble(2*i-1)*ux*d2w(6,i-1) +  &
          signx*dble(2*i-1)*dw(2,i-1) + dble((i-1)**2)*d2w(6,i-2)))*raptrm
      enddo
      if (lgrad3) then
        gtrm1 = 0.75_dp*raptrm*raptrm*w(0)
        gtrm2 = - 2.5_dp*gtrm1*raptrm
!
!  Cartesian - cartesian - cartesian derivatives of W
!
        d3w(1,0) = gtrm2*d1g(1)*d1g(1)*d1g(1) + gtrm1*(3.0_dp*d2g(1)*d1g(1))
        d3w(2,0) = gtrm2*d1g(2)*d1g(1)*d1g(1) + gtrm1*(2.0_dp*d2g(6)*d1g(1) + d2g(1)*d1g(2))
        d3w(3,0) = gtrm2*d1g(3)*d1g(1)*d1g(1) + gtrm1*(2.0_dp*d2g(5)*d1g(1) + d2g(1)*d1g(3))
        d3w(4,0) = gtrm2*d1g(2)*d1g(2)*d1g(1) + gtrm1*(2.0_dp*d2g(6)*d1g(2) + d2g(2)*d1g(1))
        d3w(5,0) = gtrm2*d1g(3)*d1g(2)*d1g(1) + gtrm1*(d2g(4)*d1g(1) + d2g(5)*d1g(2) + d2g(6)*d1g(3))
        d3w(6,0) = gtrm2*d1g(3)*d1g(3)*d1g(1) + gtrm1*(2.0_dp*d2g(5)*d1g(1) + d2g(3)*d1g(1))
        d3w(7,0) = gtrm2*d1g(2)*d1g(2)*d1g(2) + gtrm1*(3.0_dp*d2g(2)*d1g(2))
        d3w(8,0) = gtrm2*d1g(3)*d1g(2)*d1g(2) + gtrm1*(2.0_dp*d2g(4)*d1g(2) + d2g(2)*d1g(3))
        d3w(9,0) = gtrm2*d1g(3)*d1g(3)*d1g(2) + gtrm1*(2.0_dp*d2g(4)*d1g(3) + d2g(3)*d1g(2))
        d3w(10,0)= gtrm2*d1g(3)*d1g(3)*d1g(3) + gtrm1*(3.0_dp*d2g(3)*d1g(3))
!
        d3w(1,1) = -6.0_dp*ux*(dw(1,0)*dw(1,0)*dw(1,0) + w(0)*(3.0_dp*d2w(1,0)*dw(1,0) +  &
          0.5_dp*w(0)*d3w(1,0))) - 3.0_dp*signx*w(0)*(6.0_dp*dw(1,0)*dw(1,0) + w(0)*(3.0_dp*d2w(1,0)))
        d3w(2,1) = -6.0_dp*ux*(dw(2,0)*dw(1,0)*dw(1,0) + w(0)*(2.0_dp*d2w(6,0)*dw(1,0) + d2w(1,0)*dw(2,0) + &
          0.5_dp*w(0)*d3w(2,0))) - 6.0_dp*signx*w(0)*(2.0_dp*dw(2,0)*dw(1,0) + w(0)*d2w(6,0))
        d3w(3,1) = -6.0_dp*ux*(dw(3,0)*dw(1,0)*dw(1,0) + w(0)*(2.0_dp*d2w(5,0)*dw(1,0) + d2w(1,0)*dw(3,0) + &
          0.5_dp*w(0)*d3w(3,0))) - 6.0_dp*signx*w(0)*(2.0_dp*dw(3,0)*dw(1,0) + w(0)*d2w(5,0))
        d3w(4,1) = -6.0_dp*ux*(dw(2,0)*dw(2,0)*dw(1,0) + w(0)*(2.0_dp*d2w(6,0)*dw(2,0) + d2w(2,0)*dw(1,0) + &
          0.5_dp*w(0)*d3w(4,0))) - 6.0_dp*signx*w(0)*(dw(2,0)*dw(2,0) + 0.5_dp*w(0)*d2w(2,0))
        d3w(5,1) = -6.0_dp*ux*(dw(3,0)*dw(2,0)*dw(1,0) + w(0)*(d2w(4,0)*dw(1,0) + d2w(5,0)*dw(2,0) + &
          d2w(6,0)*dw(3,0) + 0.5_dp*w(0)*d3w(5,0))) - 6.0_dp*signx*w(0)*(dw(3,0)*dw(2,0) + 0.5_dp*w(0)*d2w(4,0))
        d3w(6,1) = -6.0_dp*ux*(dw(3,0)*dw(3,0)*dw(1,0) + w(0)*(2.0_dp*d2w(5,0)*dw(3,0) + d2w(3,0)*dw(1,0) + &
          0.5_dp*w(0)*d3w(6,0))) - 6.0_dp*signx*w(0)*(dw(3,0)*dw(3,0) + 0.5_dp*w(0)*d2w(3,0))
        d3w(7,1) = -6.0_dp*ux*(dw(2,0)*dw(2,0)*dw(2,0) + w(0)*(3.0_dp*d2w(2,0)*dw(2,0) + 0.5_dp*w(0)*d3w(7,0)))
        d3w(8,1) = -6.0_dp*ux*(dw(3,0)*dw(2,0)*dw(2,0) + w(0)*(2.0_dp*d2w(5,0)*dw(2,0) + d2w(2,0)*dw(3,0) + &
          0.5_dp*w(0)*d3w(8,0)))
        d3w(9,1) = -6.0_dp*ux*(dw(3,0)*dw(3,0)*dw(2,0) + w(0)*(2.0_dp*d2w(5,0)*dw(3,0) + d2w(3,0)*dw(2,0) + &
          0.5_dp*w(0)*d3w(9,0)))
        d3w(10,1) = -6.0_dp*ux*(dw(3,0)*dw(3,0)*dw(3,0) + w(0)*(3.0_dp*d2w(3,0)*dw(3,0) + 0.5_dp*w(0)*d3w(10,0)))
!
        do i = 2,2*lm-1
          d3w(1,i) = - raptrm*(3.0_dp*(d2w(1,i)*d1g(1) + dw(1,i)*d2g(1)) + &
            dble(2*i-1)*(signx*(3.0_dp*d2w(1,i-1)) + ux*d3w(1,i-1)) + dble((i-1)**2)*d3w(1,i-2))
          d3w(2,i) = - raptrm*( &
            2.0_dp*d2w(6,i)*d1g(1) + d2w(1,i)*d1g(2) + 2.0_dp*dw(1,i)*d2g(6) + dw(2,i)*d2g(1) + &
            dble(2*i-1)*(signx*(2.0_dp*d2w(6,i-1)) + ux*d3w(2,i-1)) + dble((i-1)**2)*d3w(2,i-2))
          d3w(3,i) = - raptrm*(2.0_dp*d2w(5,i)*d1g(1) + d2w(1,i)*d1g(3) + &
            2.0_dp*dw(1,i)*d2g(5) + dw(3,i)*d2g(1) + dble(2*i-1)*(signx*(2.0_dp*d2w(5,i-1)) +  &
            ux*d3w(3,i-1)) + dble((i-1)**2)*d3w(3,i-2))
          d3w(4,i) = - raptrm*(2.0_dp*d2w(6,i)*d1g(2) + d2w(2,i)*d1g(1) + 2.0_dp*dw(2,i)*d2g(6) + dw(1,i)*d2g(2) + &
            dble(2*i-1)*(signx*(d2w(2,i-1)) + ux*d3w(4,i-1)) + dble((i-1)**2)*d3w(4,i-2))
          d3w(5,i) = - raptrm*(d2w(4,i)*d1g(1) + d2w(5,i)*d1g(2) + d2w(6,i)*d1g(3) + dw(1,i)*d2g(4) +  &
            dw(2,i)*d2g(5) + dw(3,i)*d2g(6) + dble(2*i-1)*(signx*(d2w(4,i-1)) + ux*d3w(5,i-1)) + dble((i-1)**2)*d3w(5,i-2))
          d3w(6,i) = - raptrm*(2.0_dp*d2w(5,i)*d1g(3) + d2w(3,i)*d1g(1) + 2.0_dp*dw(3,i)*d2g(5) + dw(1,i)*d2g(3) + &
            dble(2*i-1)*(signx*(d2w(3,i-1)) + ux*d3w(6,i-1)) + dble((i-1)**2)*d3w(6,i-2))
          d3w(7,i) = - raptrm*(3.0_dp*(d2w(2,i)*d1g(2) + dw(2,i)*d2g(2)) + &
            dble(2*i-1)*(ux*d3w(7,i-1)) + dble((i-1)**2)*d3w(7,i-2))
          d3w(8,i) = - raptrm*(2.0_dp*d2w(4,i)*d1g(2) + d2w(2,i)*d1g(3) + &
            2.0_dp*dw(2,i)*d2g(4) + dw(3,i)*d2g(2) + dble(2*i-1)*(ux*d3w(8,i-1)) + dble((i-1)**2)*d3w(8,i-2))
          d3w(9,i) = - raptrm*(2.0_dp*d2w(4,i)*d1g(3) + d2w(3,i)*d1g(2) + 2.0_dp*dw(3,i)*d2g(4) + dw(2,i)*d2g(3) + &
            dble(2*i-1)*(ux*d3w(9,i-1)) + dble((i-1)**2)*d3w(9,i-2))
          d3w(10,i) = - raptrm*(3.0_dp*(d2w(3,i)*d1g(3) + dw(3,i)*d2g(3)) + &
            dble(2*i-1)*(ux*d3w(10,i-1)) + dble((i-1)**2)*d3w(10,i-2))
        enddo
!
!  Strain - cartesian - cartesian derivatives of W
!
        do i = 0,2*lm-1
          d3ws(1,i) = d2w(1,i+1)*u + d3w(1,i)*x*signx + 2.0_dp*d2w(1,i)*signx
          d3ws(2,i) = d2w(2,i+1)*u + d3w(4,i)*x*signx
          d3ws(3,i) = d2w(3,i+1)*u + d3w(6,i)*x*signx
          d3ws(4,i) = d2w(4,i+1)*u + d3w(5,i)*x*signx
          d3ws(5,i) = d2w(5,i+1)*u + d3w(3,i)*x*signx + d2w(5,i)*signx
          d3ws(6,i) = d2w(6,i+1)*u + d3w(2,i)*x*signx + d2w(6,i)*signx
        enddo
      endif
    endif
!
!  Initialise total derivatives to zero
!
    d1emf(1:3) = 0.0_dp
    d1emfs = 0.0_dp
    if (lgrad2) then
      d2emfs = 0.0_dp
      d2emf(1:6) = 0.0_dp
      d2emfm(1:3) = 0.0_dp
      if (lgrad3) then
        d3emf(1:10) = 0.0_dp
        d3emfm(1:6) = 0.0_dp
      endif
    endif
    ap = a
    do i = 1,lm
!
!  First derivatives : Cartesian
!
      d1emf(1) = d1emf(1) - ecoeff(i)*ap*dw(1,2*i-1)
      d1emf(2) = d1emf(2) - ecoeff(i)*ap*dw(2,2*i-1)
      d1emf(3) = d1emf(3) - ecoeff(i)*ap*dw(3,2*i-1)
!
!  First derivatives : strain
!
      dwde = dw(1,2*i-1)*x*signx + w(2*i)*u
      trm1s = ecoeff(i)*ap*(dble(2*i-1)*w(2*i-1)*dads + dws(2*i-1))
      d1emfs = d1emfs - trm1s
      if (lgrad2) then
!
!  Second derivatives : cartesian - cartesian
!
        d2emf(1) = d2emf(1) - ecoeff(i)*ap*d2w(1,2*i-1)
        d2emf(2) = d2emf(2) - ecoeff(i)*ap*d2w(2,2*i-1)
        d2emf(3) = d2emf(3) - ecoeff(i)*ap*d2w(3,2*i-1)
        d2emf(4) = d2emf(4) - ecoeff(i)*ap*d2w(4,2*i-1)
        d2emf(5) = d2emf(5) - ecoeff(i)*ap*d2w(5,2*i-1)
        d2emf(6) = d2emf(6) - ecoeff(i)*ap*d2w(6,2*i-1)
!
!  Second derivatives : strain - strain
!
        d2emfs = d2emfs - ecoeff(i)*ap*(dble(2*i-1)*dws(2*i-1)*dads + d2ws(2*i-1)) &
                        - ecoeff(i)*ap*dble(2*i-1)*dads*(dble(2*i-2)*w(2*i-1)*dads + dws(2*i-1)) &
                        - ecoeff(i)*ap*dble(2*i-1)*w(2*i-1)*d2ads2
!
!  Second derivatives : cartesian - strain
!
        d2emfm(1) = d2emfm(1) - ecoeff(i)*ap*(dads*dw(1,2*i-1) + d2wm(1,2*i-1))
        d2emfm(2) = d2emfm(2) - ecoeff(i)*ap*(dads*dw(2,2*i-1) + d2wm(2,2*i-1))
        d2emfm(3) = d2emfm(3) - ecoeff(i)*ap*(dads*dw(3,2*i-1) + d2wm(3,2*i-1))
!
        if (lgrad3) then
!
!  Third derivatives : cartesian - cartesian - cartesian
!
          do ii = 1,10
            d3emf(ii) = d3emf(ii) - ecoeff(i)*ap*d3w(ii,2*i-1)
          enddo
!
!  Third derivatives : cartesian - cartesian - strain
!
          do ii = 1,6
            d3emfm(ii) = d3emfm(ii) - ecoeff(i)*ap*(d3ws(ii,2*i-1) + dble(2*i-1)*d2w(ii,2*i-1))
          enddo
        endif
      endif
      ap = ap*a*a
    enddo
  endif
#ifdef TRACE
  call trace_out('emfuncs')
#endif
!
  return
  end
