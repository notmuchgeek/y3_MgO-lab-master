  subroutine reciptrmd3dV(i,j,nkp,c6tot,nkvec0,ktrm0,ktrm20, &
    ktrm6,ktrm60,ktrm62,ktrm620,d3s,d3rs,d3is)
!
!  Calculates the reciprocal space contribution to the third
!  derivatives of frequencies respect to volume when called
!  from realrecip3d3dV.
!
!  nkp    = k point to be calculated
!
!  The use of right angled saving must be disabled in kindex when
!  calling this routine as the k point breaks the mirror plane
!  symmetry.
!
!   1/18 Created from reciptrmd3
!   2/18 Trace added
!   3/20 Tolerance for ldoc6 made global
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
  use g_constants
  use current
  use ksample
  use kspace
  use thresholds,    only : thresh_c6
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4)        :: i
  integer(i4)        :: j
  integer(i4)        :: nkp
  integer(i4)        :: nkvec0
  real(dp)           :: c6tot
  real(dp)           :: ktrm6(*)
  real(dp)           :: ktrm62(*)
  real(dp)           :: ktrm0(*)
  real(dp)           :: ktrm20(*)
  real(dp)           :: ktrm60(*)
  real(dp)           :: ktrm620(*)
  real(dp)           :: d3is(3,3,3)
  real(dp)           :: d3rs(3,3,3)
  real(dp)           :: d3s(3,3,3)
!
!  Local variables
!
  integer(i4)        :: kk
  integer(i4)        :: iv
  logical            :: ldoc6
  logical            :: lgamma
  real(dp)           :: arg
  real(dp)           :: cosa
  real(dp)           :: cosq
  real(dp)           :: coss
  real(dp)           :: dt1
  real(dp)           :: oci
  real(dp)           :: ocj
  real(dp)           :: qfct
  real(dp)           :: qli
  real(dp)           :: qlj
  real(dp)           :: rk1x
  real(dp)           :: rk1y
  real(dp)           :: rk1z
  real(dp)           :: rk2x
  real(dp)           :: rk2y
  real(dp)           :: rk2z
  real(dp)           :: rk3x
  real(dp)           :: rk3y
  real(dp)           :: rk3z
  real(dp)           :: sina
  real(dp)           :: sinq
  real(dp)           :: sins
  real(dp)           :: sktrm1
  real(dp)           :: sktrm2
  real(dp)           :: tmp1
  real(dp)           :: tmp2
  real(dp)           :: tmp3
  real(dp)           :: tmp4
  real(dp)           :: tmp5
  real(dp)           :: tmp6
  real(dp)           :: xd
  real(dp)           :: yd
  real(dp)           :: zd
  real(dp)           :: xk
  real(dp)           :: yk
  real(dp)           :: zk
  real(dp)           :: xkv
  real(dp)           :: ykv
  real(dp)           :: zkv
  real(dp)           :: xrkk
  real(dp)           :: yrkk
  real(dp)           :: zrkk
#ifdef TRACE
  call trace_in('reciptrmd3dV')
#endif
!
!  Assign local variables
!
  rk1x = kv(1,1)
  rk1y = kv(2,1)
  rk1z = kv(3,1)
  rk2x = kv(1,2)
  rk2y = kv(2,2)
  rk2z = kv(3,2)
  rk3x = kv(1,3)
  rk3y = kv(2,3)
  rk3z = kv(3,3)
!
!  Calculate basic k vector
!
  xk = xkpt(nkp)
  yk = ykpt(nkp)
  zk = zkpt(nkp)
  lgamma = ((abs(xk)+abs(yk)+abs(zk)).lt.1.0d-8)
  xkv = xk*rk1x + yk*rk2x + zk*rk3x
  ykv = xk*rk1y + yk*rk2y + zk*rk3y
  zkv = xk*rk1z + yk*rk2z + zk*rk3z
!************************************************
!  Algorithm for cases where dispersion cannot  *
!  be factorised into one centre terms          *
!************************************************
  oci = occuf(i)
  qli = qf(i)*oci
  ocj = occuf(j)
  qlj = qf(j)*ocj
  qfct = qli*qlj
  c6tot = c6tot*oci*ocj
  ldoc6 = (abs(c6tot).gt.thresh_c6)
!
!  Find relative vector between atoms
!
  xd = xclat(j) - xclat(i)
  yd = yclat(j) - yclat(i)
  zd = zclat(j) - zclat(i)
!********************************
!  Calculate third derivatives  *
!********************************
  do iv = 1,nkvec
    xrkk = xrk(iv)
    yrkk = yrk(iv)
    zrkk = zrk(iv)
!
!  3-D
!
    arg = (xrkk-xkv)*xd + (yrkk-ykv)*yd + (zrkk-zkv)*zd
    if (ldoc6) then
      cosa = cos(arg)
      sina = sin(arg)
      sktrm1 = (ktrm(iv)*qfct - ktrm6(iv)*c6tot)
      sktrm2 = (ktrms(iv)*qfct - ktrm62(iv)*c6tot)
      cosq = cosa*sktrm1
      sinq = sina*sktrm1
      coss = cosa*sktrm2
      sins = sina*sktrm2
    else
      cosa = cos(arg)*qfct
      sina = sin(arg)*qfct
      cosq = cosa*ktrm(iv)
      sinq = sina*ktrm(iv)
      coss = cosa*ktrms(iv)
      sins = sina*ktrms(iv)
    endif
!
    tmp1 = xrkk*xrkk
    tmp2 = yrkk*yrkk
    tmp3 = zrkk*zrkk
    tmp4 = yrkk*zrkk
    tmp5 = xrkk*zrkk
    tmp6 = xrkk*yrkk
!
! Calculate real and imaginary third derivatives
!
! no z component
!
    d3rs(1,1,1) = d3rs(1,1,1) + coss*tmp1*tmp1 + 3.0_dp*cosq*tmp1
    d3rs(2,1,1) = d3rs(2,1,1) + coss*tmp6*tmp1 + 2.0_dp*cosq*tmp6
    d3rs(2,2,1) = d3rs(2,2,1) + coss*tmp2*tmp1 + cosq*tmp2
    d3rs(1,1,2) = d3rs(1,1,2) + coss*tmp1*tmp2 + cosq*tmp1
    d3rs(2,1,2) = d3rs(2,1,2) + coss*tmp6*tmp2 + 2.0_dp*cosq*tmp6
    d3rs(2,2,2) = d3rs(2,2,2) + coss*tmp2*tmp2 + 3.0_dp*cosq*tmp2
!
    d3is(1,1,1) = d3is(1,1,1) - sins*tmp1*tmp1 - 3.0_dp*sinq*tmp1
    d3is(2,1,1) = d3is(2,1,1) - sins*tmp6*tmp1 - 2.0_dp*sinq*tmp6
    d3is(2,2,1) = d3is(2,2,1) - sins*tmp2*tmp1 - sinq*tmp2
    d3is(1,1,2) = d3is(1,1,2) - sins*tmp1*tmp2 - sinq*tmp1
    d3is(2,1,2) = d3is(2,1,2) - sins*tmp6*tmp2 - 2.0_dp*sinq*tmp6
    d3is(2,2,2) = d3is(2,2,2) - sins*tmp2*tmp2 - 3.0_dp*sinq*tmp2
!
    d3rs(3,1,1) = d3rs(3,1,1) + coss*tmp5*tmp1 + 2.0_dp*cosq*tmp5
    d3rs(3,2,1) = d3rs(3,2,1) + coss*tmp4*tmp1 + cosq*tmp4
    d3rs(3,3,1) = d3rs(3,3,1) + coss*tmp3*tmp1 + cosq*tmp3
    d3rs(3,1,2) = d3rs(3,1,2) + coss*tmp5*tmp2 + cosq*tmp5
    d3rs(3,2,2) = d3rs(3,2,2) + coss*tmp4*tmp2 + 2.0_dp*cosq*tmp4
    d3rs(3,3,2) = d3rs(3,3,2) + coss*tmp3*tmp2 + cosq*tmp3
    d3rs(1,1,3) = d3rs(1,1,3) + coss*tmp1*tmp3 + cosq*tmp1
    d3rs(2,1,3) = d3rs(2,1,3) + coss*tmp6*tmp3 + cosq*tmp6
    d3rs(2,2,3) = d3rs(2,2,3) + coss*tmp2*tmp3 + cosq*tmp2
    d3rs(3,1,3) = d3rs(3,1,3) + coss*tmp5*tmp3 + 2.0_dp*cosq*tmp5
    d3rs(3,2,3) = d3rs(3,2,3) + coss*tmp4*tmp3 + 2.0_dp*cosq*tmp4
    d3rs(3,3,3) = d3rs(3,3,3) + coss*tmp3*tmp3 + 3.0_dp*cosq*tmp3
!
    d3is(3,1,1) = d3is(3,1,1) - sins*tmp5*tmp1 - 2.0_dp*sinq*tmp5
    d3is(3,2,1) = d3is(3,2,1) - sins*tmp4*tmp1 - sinq*tmp4
    d3is(3,3,1) = d3is(3,3,1) - sins*tmp3*tmp1 - sinq*tmp3
    d3is(3,1,2) = d3is(3,1,2) - sins*tmp5*tmp2 - sinq*tmp5
    d3is(3,2,2) = d3is(3,2,2) - sins*tmp4*tmp2 - 2.0_dp*sinq*tmp4
    d3is(3,3,2) = d3is(3,3,2) - sins*tmp3*tmp2 - sinq*tmp3
    d3is(1,1,3) = d3is(1,1,3) - sins*tmp1*tmp3 - sinq*tmp1
    d3is(2,1,3) = d3is(2,1,3) - sins*tmp6*tmp3 - sinq*tmp6
    d3is(2,2,3) = d3is(2,2,3) - sins*tmp2*tmp3 - sinq*tmp2
    d3is(3,1,3) = d3is(3,1,3) - sins*tmp5*tmp3 - 2.0_dp*sinq*tmp5
    d3is(3,2,3) = d3is(3,2,3) - sins*tmp4*tmp3 - 2.0_dp*sinq*tmp4
    d3is(3,3,3) = d3is(3,3,3) - sins*tmp3*tmp3 - 3.0_dp*sinq*tmp3
  enddo
  if (.not.lgamma) then
    do iv = 1,nkvec0
      xrkk = xrk0(iv)
      yrkk = yrk0(iv)
      zrkk = zrk0(iv)
      arg = xrkk*xd + yrkk*yd + zrkk*zd
      if (ldoc6) then
        cosa = cos(arg)
        sina = sin(arg)
        sktrm1 = (ktrm0(iv)*qfct - ktrm60(iv)*c6tot)
        sktrm2 = (ktrm20(iv)*qfct - ktrm620(iv)*c6tot)
        cosq = cosa*sktrm1
        sinq = sina*sktrm1
        coss = cosa*sktrm2
        sins = sina*sktrm2
      else
        cosa = cos(arg)*qfct
        sina = sin(arg)*qfct
        cosq = cosa*ktrm0(iv)
        sinq = sina*ktrm0(iv)
        coss = cosa*ktrm20(iv)
        sins = sina*ktrm20(iv)
      endif
      tmp1 = xrkk*xrkk
      tmp2 = yrkk*yrkk
      tmp3 = zrkk*zrkk
      tmp4 = yrkk*zrkk
      tmp5 = xrkk*zrkk
      tmp6 = xrkk*yrkk
!
!  Calculate real third derivatives
!
      dt1 = coss*tmp1
      d3s(1,1,1) = d3s(1,1,1) + dt1*tmp1
      d3s(2,1,1) = d3s(2,1,1) + dt1*tmp6
      d3s(3,1,1) = d3s(3,1,1) + dt1*tmp5
      d3s(2,2,1) = d3s(2,2,1) + dt1*tmp2
      d3s(3,2,1) = d3s(3,2,1) + dt1*tmp4
      d3s(3,3,1) = d3s(3,3,1) + dt1*tmp3
      d3s(1,1,1) = d3s(1,1,1) + 3.0_dp*cosq*tmp1
      d3s(2,1,1) = d3s(2,1,1) + 2.0_dp*cosq*tmp6
      d3s(3,1,1) = d3s(3,1,1) + 2.0_dp*cosq*tmp5
      d3s(2,2,1) = d3s(2,2,1) + cosq*tmp2
      d3s(3,2,1) = d3s(3,2,1) + cosq*tmp4
      d3s(3,3,1) = d3s(3,3,1) + cosq*tmp3
      dt1 = coss*tmp2
      d3s(1,1,2) = d3s(1,1,2) + dt1*tmp1
      d3s(2,1,2) = d3s(2,1,2) + dt1*tmp6
      d3s(3,1,2) = d3s(3,1,2) + dt1*tmp5
      d3s(2,2,2) = d3s(2,2,2) + dt1*tmp2
      d3s(3,2,2) = d3s(3,2,2) + dt1*tmp4
      d3s(3,3,2) = d3s(3,3,2) + dt1*tmp3
      d3s(1,1,2) = d3s(1,1,2) + cosq*tmp1
      d3s(2,1,2) = d3s(2,1,2) + 2.0_dp*cosq*tmp6
      d3s(3,1,2) = d3s(3,1,2) + cosq*tmp5
      d3s(2,2,2) = d3s(2,2,2) + 3.0_dp*cosq*tmp2
      d3s(3,2,2) = d3s(3,2,2) + 2.0_dp*cosq*tmp4
      d3s(3,3,2) = d3s(3,3,2) + cosq*tmp3
      dt1 = coss*tmp3
      d3s(1,1,3) = d3s(1,1,3) + dt1*tmp1
      d3s(2,1,3) = d3s(2,1,3) + dt1*tmp6
      d3s(3,1,3) = d3s(3,1,3) + dt1*tmp5
      d3s(2,2,3) = d3s(2,2,3) + dt1*tmp2
      d3s(3,2,3) = d3s(3,2,3) + dt1*tmp4
      d3s(3,3,3) = d3s(3,3,3) + dt1*tmp3
      d3s(1,1,3) = d3s(1,1,3) + cosq*tmp1
      d3s(2,1,3) = d3s(2,1,3) + cosq*tmp6
      d3s(3,1,3) = d3s(3,1,3) + 2.0_dp*cosq*tmp5
      d3s(2,2,3) = d3s(2,2,3) + cosq*tmp2
      d3s(3,2,3) = d3s(3,2,3) + 2.0_dp*cosq*tmp4
      d3s(3,3,3) = d3s(3,3,3) + 3.0_dp*cosq*tmp3
    enddo
  else
!
!  Gamma point - just copy values
!
    do kk = 1,3
      d3s(1,1,kk) = d3rs(1,1,kk)
      d3s(2,1,kk) = d3rs(2,1,kk)
      d3s(3,1,kk) = d3rs(3,1,kk)
      d3s(2,2,kk) = d3rs(2,2,kk)
      d3s(3,2,kk) = d3rs(3,2,kk)
      d3s(3,3,kk) = d3rs(3,3,kk)
    enddo
  endif
#ifdef TRACE
  call trace_out('reciptrmd3dV')
#endif
!
  return
  end
