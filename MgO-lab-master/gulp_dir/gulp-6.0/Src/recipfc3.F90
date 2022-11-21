  subroutine recipfc3(i,j,c6tot,ktrm6,d3)
!
!  Calculates the reciprocal space contribution to the third
!  derivatives when called from thirdorderfc3
!
!  The use of right angled saving must be disabled in kindex when
!  calling this routine as the k point breaks the mirror plane
!  symmetry.
!
!   4/15 Created from reciptrmd3
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
  real(dp)           :: c6tot
  real(dp)           :: ktrm6(*)
  real(dp)           :: d3(3,3,3)
!
!  Local variables
!
  integer(i4)        :: iv
  logical            :: ldoc6
  real(dp)           :: arg
  real(dp)           :: cosa
  real(dp)           :: cosr
  real(dp)           :: d11m22
  real(dp)           :: d13m24
  real(dp)           :: d13p24
  real(dp)           :: d3t1dz3
  real(dp)           :: darg1
  real(dp)           :: darg2
  real(dp)           :: g_derfc
  real(dp)           :: derfc1
  real(dp)           :: derfc2
  real(dp)           :: dexp1
  real(dp)           :: dexp2
  real(dp)           :: dexp3
  real(dp)           :: dexp4
  real(dp)           :: dexpz
  real(dp)           :: dtrm2
  real(dp)           :: etaz
  real(dp)           :: etaz2
  real(dp)           :: kexperfc
  real(dp)           :: kvec
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
  real(dp)           :: sinr
  real(dp)           :: sktrm1
  real(dp)           :: tmp1
  real(dp)           :: tmp2
  real(dp)           :: tmp3
  real(dp)           :: tmp4
  real(dp)           :: tmp5
  real(dp)           :: tmp6
  real(dp)           :: trmzz
  real(dp)           :: trmzzz
  real(dp)           :: xd
  real(dp)           :: yd
  real(dp)           :: zd
  real(dp)           :: xrkk
  real(dp)           :: yrkk
  real(dp)           :: zrkk
  real(dp)           :: ztrm1
#ifdef TRACE
  call trace_in('recipfc3')
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
!
!  Define K vector independent constants and T1 surface derivatives
!
  if (ndim.eq.2) then
!
!  Define constants
!
    rpieta = 1.0_dp / sqrt(pi * eta)
    rhseta = 0.5_dp / seta
!
!  Evaluate K -> 0 limit of 2-D sum
!
    etaz = seta*zd
    etaz2 = etaz * etaz
    dexpz = exp(-etaz2)*angstoev
    dtrm2 = - qfct*vol4pi*tweatpi*dexpz
    d3t1dz3 = 2.0_dp*qfct*vol4pi*tweatpi*dexpz*eta*zd
    d3(3,3,3)  = d3(3,3,3)  + d3t1dz3
  endif
!********************************
!  Calculate third derivatives  *
!********************************
  do iv = 1,nkvec
    xrkk = xrk(iv)
    yrkk = yrk(iv)
    zrkk = zrk(iv)
    if (ndim.eq.3) then
!
!  3-D
!
      arg = xrkk*xd + yrkk*yd + zrkk*zd
      if (ldoc6) then
        cosa = cos(arg)
        sina = sin(arg)
        sktrm1 = (ktrm(iv)*qfct - ktrm6(iv)*c6tot)
        sinq = sina*sktrm1
      else
        cosa = cos(arg)*qfct
        sina = sin(arg)*qfct
        sinq = sina*ktrm(iv)
      endif
    elseif (ndim.eq.2) then
!
!  2-D
!
      arg = xrkk*xd + yrkk*yd
      cosa = cos(arg)*qfct
      sina = sin(arg)*qfct
!           
      kvec = kmod(iv)
      dexp1 = exp(kvec*zd)
      dexp2 = 1.0_dp/dexp1
      darg1 = kvec*rhseta + etaz
      darg2 = kvec*rhseta - etaz
      dexp3 = exp(-(darg1)**2)
      dexp4 = exp(-(darg2)**2)
      derfc1 = g_derfc(darg1)
      derfc2 = g_derfc(darg2)
      d11m22 = dexp1*derfc1 - dexp2*derfc2
      d13m24 = dexp1*dexp3 - dexp2*dexp4
      d13p24 = dexp1*dexp3 + dexp2*dexp4
      kexperfc = dexp1*derfc1 + dexp2*derfc2
      ztrm1 = kvec*d11m22 - tweatpi*d13m24
      trmzz = (kvec*kvec*kexperfc - 2.0_dp*tweatpi*(kvec*d13p24 -  &
               seta*(darg1*dexp1*dexp3 + darg2*dexp2*dexp4)))
      trmzzz = kvec*kvec*kvec*d11m22 - tweatpi*(3.0_dp*kvec*kvec*d13m24  &
      - 6.0_dp*kvec*seta*(dexp1*darg1*dexp3-dexp2*darg2*dexp4) &
      + 4.0*eta*(dexp1*darg1*darg1*dexp3-dexp2*darg2*darg2*dexp4) &
      - 2.0_dp*eta*d13m24)
!            
      cosr = cosa*ktrm(iv)
      sinr = sina*ktrm(iv)
      sinq = sinr*kexperfc
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
    d3(1,1,1) = d3(1,1,1) + sinq*tmp1*xrkk
    d3(2,1,1) = d3(2,1,1) + sinq*tmp6*xrkk
    d3(2,2,1) = d3(2,2,1) + sinq*tmp2*xrkk
    d3(2,2,2) = d3(2,2,2) + sinq*tmp2*yrkk
    if (ndim.eq.3) then
      d3(3,1,1) = d3(3,1,1) + sinq*tmp5*xrkk
      d3(3,2,1) = d3(3,2,1) + sinq*tmp4*xrkk
      d3(3,3,1) = d3(3,3,1) + sinq*tmp3*xrkk
      d3(3,2,2) = d3(3,2,2) + sinq*tmp4*yrkk
      d3(3,3,2) = d3(3,3,2) + sinq*tmp3*yrkk
      d3(3,3,3) = d3(3,3,3) + sinq*tmp3*zrkk
    elseif (ndim.eq.2) then
! 1 z component
      d3(3,1,1) = d3(3,1,1) - cosr*tmp1*ztrm1
      d3(3,2,1) = d3(3,2,1) - cosr*tmp6*ztrm1
      d3(3,2,2) = d3(3,2,2) - cosr*tmp2*ztrm1
! 2 z components
      d3(3,3,1) = d3(3,3,1) - sinr*trmzz*xrkk
      d3(3,3,2) = d3(3,3,2) - sinr*trmzz*yrkk
! 3 z components
      d3(3,3,3) = d3(3,3,3) + cosr*trmzzz
    endif
  enddo
#ifdef TRACE
  call trace_out('recipfc3')
#endif
!
  return
  end
