  subroutine setktrm4(esum4)
!
!  Sets k vector related terms needed for reciprocal space sums.
!
!  11/04 Intent added and sqrt pi taken from module
!  12/07 Unused variables removed
!  10/13 Hardwired maximum for k vector index replaced with maxindk
!   3/14 derfc changed to g_derfc for benefit of ChemShell
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
!  Copyright Curtin University 2014
!
!  Julian Gale, CIC, Curtin University, March 2014
!
  use g_constants
  use current
  use kspace
  use symmetry
  implicit none
!
!  Passed variables
!
  real(dp),   intent(inout) :: esum4
!
!  Local variables
!
  integer(i4)             :: e
  integer(i4)             :: f
  integer(i4)             :: g
  integer(i4)             :: i
  integer(i4)             :: idk
  real(dp)                :: arge
  real(dp)                :: cerf
  real(dp)                :: const
  real(dp)                :: ctrm
  real(dp)                :: g_derfc
  real(dp)                :: factor
  real(dp)                :: kvv(3)
  real(dp)                :: rk
  real(dp)                :: rk2
  real(dp)                :: rkt
  real(dp)                :: rrk2
  real(dp)                :: vol
  real(dp)                :: volume
  real(dp)                :: xpon
!
  eta4 = 0.25_dp/eta
  vol = volume(rv)
  const = 2.0_dp*pi*sqrtpi/vol
!
!  Add zero wavevector term
!
  esum4 = esum4 + const*seta
!
!  Calculate and store k - vectors
!
  if (lra) then
    kvv(1) = kv(1,1)
    kvv(2) = kv(2,2)
    kvv(3) = kv(3,3)
    do i = 1,nkvec
      idk = indk(i)
      e = (idk/maxindk3) - maxindk
      if (e.eq.0) then
        factor = 1.0_dp
      else
        factor = 2.0_dp
      endif
      idk = idk - (e + maxindk)*maxindk3
      f = (idk/maxindk2) - maxindk
      g = idk - (f + maxindk)*maxindk2 - maxindk
      xrk(i) = e*kvv(1)
      yrk(i) = f*kvv(2)
      zrk(i) = g*kvv(3)
      rk2 = xrk(i)*xrk(i)
      rk2 = rk2 + yrk(i)*yrk(i)
      rk2 = rk2 + zrk(i)*zrk(i)
      rrk2 = sqrt(rk2)
      rk = 0.5_dp*rrk2
      arge =  - rk2*eta4
      xpon = exp(arge)
      rkt = rk/seta
      cerf = g_derfc(rkt)
      ktrm(i) = factor*(seta*xpon - rk*sqrtpi*cerf)
      sine(i) = 2.0_dp*factor*cerf/rrk2
    enddo
  else
    do i = 1,nkvec
      idk = indk(i)
      e = (idk/maxindk3) - maxindk
      idk = idk - (e + maxindk)*maxindk3
      f = (idk/maxindk2) - maxindk
      g = idk - (f + maxindk)*maxindk2 - maxindk
      factor = 2.0_dp
      if (e.eq.0.and.nkangle.eq.1) then
        factor = 1.0_dp
      elseif (f.eq.0.and.nkangle.eq.2) then
        factor = 1.0_dp
      elseif (g.eq.0.and.nkangle.eq.3) then
        factor = 1.0_dp
      elseif (nkangle.eq.0) then
        factor = 1.0_dp
      endif
      xrk(i) = e*kv(1,1) + f*kv(1,2) + g*kv(1,3)
      yrk(i) = e*kv(2,1) + f*kv(2,2) + g*kv(2,3)
      zrk(i) = e*kv(3,1) + f*kv(3,2) + g*kv(3,3)
      rk2 = xrk(i)*xrk(i)
      rk2 = rk2 + yrk(i)*yrk(i)
      rk2 = rk2 + zrk(i)*zrk(i)
      rrk2 = sqrt(rk2)
      rk = 0.5_dp*rrk2
      arge =  - rk2*eta4
      xpon = exp(arge)
      rkt = rk/seta
      cerf = g_derfc(rkt)
      ktrm(i) = factor*(seta*xpon - rk*sqrtpi*cerf)
      sine(i) = 2.0_dp*factor*cerf/rrk2
    enddo
  endif
!
!  Store terms for 1/r**4 sum in ktrm()
!
  do i = 1,nkvec
    ktrm(i) = ktrm(i)*const
  enddo
!
!  Store terms for ra*rb/r**6 in sine()
!
  ctrm =  - pi*pi/vol
  do i = 1,nkvec
    sine(i) = sine(i)*ctrm
  enddo
!
  return
  end
