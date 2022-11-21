  subroutine calcf2D(kvec,zd,seta,f,f1z,f1g,f2zz,f2zg,f2gg,f3zzz,f3zzg,f3zgg, &
                     lgrad1,lgrad2,lgrad3)
!
!  Computes function f needed for 2D Parry sum
!  Here f(G,z) = exp(|G|.r).erfc((|G|/(2.seta) + seta*z) + 
!                exp(-|G|.r).erfc((|G|/(2.seta) - seta*z)
!
!  Input:
!
!    kvec = modulus of reciprocal lattice vector
!    zd   = z coordinate
!    seta = square root of eta
!
!  Output:
!
!    Function and derivatives as required. Note, that the third derivative
!    of f with respect to 3 reciprocal lattice vectors is not required at
!    present and so not computed here.
!
!   6/19 Created
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
!  Copyright Curtin University 2019
!
!  Julian Gale, Curtin University, June  2019
!
  use g_constants
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  logical,     intent(in)   :: lgrad1       ! If true then compute first derivatives
  logical,     intent(in)   :: lgrad2       ! If true then compute second derivatives
  logical,     intent(in)   :: lgrad3       ! If true then compute third derivatives
  real(dp),    intent(in)   :: kvec         ! Modulus of reciprocal lattice vector
  real(dp),    intent(in)   :: zd           ! z coordinate
  real(dp),    intent(in)   :: seta         ! Square root of eta
  real(dp),    intent(out)  :: f            ! Function
  real(dp),    intent(out)  :: f1z          ! First derivative of f w.r.t. z
  real(dp),    intent(out)  :: f1g          ! First derivative of f w.r.t. |G| divided by |G|
  real(dp),    intent(out)  :: f2zz         ! Second derivative of f w.r.t. z (x 2)
  real(dp),    intent(out)  :: f2zg         ! Second derivative of f w.r.t. z and |G| divided by |G|
  real(dp),    intent(out)  :: f2gg         ! Second derivative of f w.r.t. |G| divided by |G| (x 2)
  real(dp),    intent(out)  :: f3zzz        ! Third derivative of f w.r.t. z (x 3)
  real(dp),    intent(out)  :: f3zzg        ! Third derivative of f w.r.t. z (x 2) and |G| divided by |G|
  real(dp),    intent(out)  :: f3zgg        ! Third derivative of f w.r.t. z and |G| divided by |G| (x 2)
!
!  Local variables
!
  real(dp)                  :: darg1
  real(dp)                  :: darg2
  real(dp)                  :: derfc1
  real(dp)                  :: derfc2
  real(dp)                  :: d1erfc1g
  real(dp)                  :: d1erfc2g
  real(dp)                  :: d1erfc1z
  real(dp)                  :: d1erfc2z
  real(dp)                  :: d2erfc1gg
  real(dp)                  :: d2erfc2gg
  real(dp)                  :: d2erfc1zg
  real(dp)                  :: d2erfc2zg
  real(dp)                  :: d2erfc1zz
  real(dp)                  :: d2erfc2zz
  real(dp)                  :: d3erfc1zgg
  real(dp)                  :: d3erfc2zgg
  real(dp)                  :: d3erfc1zzg
  real(dp)                  :: d3erfc2zzg
  real(dp)                  :: d3erfc1zzz
  real(dp)                  :: d3erfc2zzz
  real(dp)                  :: dexp1
  real(dp)                  :: dexp2
  real(dp)                  :: dexp3
  real(dp)                  :: dexp4
  real(dp)                  :: d1exp1g
  real(dp)                  :: d1exp2g
  real(dp)                  :: d1exp1z
  real(dp)                  :: d1exp2z
  real(dp)                  :: d2exp1gg
  real(dp)                  :: d2exp2gg
  real(dp)                  :: d2exp1zg
  real(dp)                  :: d2exp2zg
  real(dp)                  :: d2exp1zz
  real(dp)                  :: d2exp2zz
  real(dp)                  :: d3exp1zgg
  real(dp)                  :: d3exp2zgg
  real(dp)                  :: d3exp1zzg
  real(dp)                  :: d3exp2zzg
  real(dp)                  :: d3exp1zzz
  real(dp)                  :: d3exp2zzz
  real(dp)                  :: g_derfc
  real(dp)                  :: etaz             ! seta*zd
  real(dp)                  :: rhseta           ! 1/(2*seta)
  real(dp)                  :: rkvec            ! 1/|G|
  real(dp)                  :: tworpi           ! 2/sqrt(pi)
#ifdef TRACE
  call trace_in('getf2D')
#endif
!
!  Set up terms needed for f
!
  rhseta = 0.5_dp/seta
  etaz   = seta*zd
  dexp1 = exp(kvec*zd)
  dexp2 = exp(-kvec*zd)
  darg1 = kvec*rhseta + etaz
  darg2 = kvec*rhseta - etaz
  derfc1 = g_derfc(darg1)
  derfc2 = g_derfc(darg2)
!
!  Function
!
  f = dexp1*derfc1 + dexp2*derfc2
!
!  First derivatives
!
  if (lgrad1) then
    rkvec = 1.0_dp/kvec
    tworpi = 2.0_dp/sqrt(pi)
    dexp3 = tworpi*exp(-(darg1)**2)
    dexp4 = tworpi*exp(-(darg2)**2)
!
    d1exp1z = kvec*dexp1
    d1exp2z = - kvec*dexp2
    d1exp1g = zd*dexp1*rkvec
    d1exp2g = - zd*dexp2*rkvec
    d1erfc1z = - seta*dexp3
    d1erfc2z = + seta*dexp4
    d1erfc1g = - rhseta*dexp3*rkvec
    d1erfc2g = - rhseta*dexp4*rkvec
!
    f1z = dexp1*d1erfc1z + dexp2*d1erfc2z + d1exp1z*derfc1 + d1exp2z*derfc2
    f1g = dexp1*d1erfc1g + dexp2*d1erfc2g + d1exp1g*derfc1 + d1exp2g*derfc2
!
!  Second derivatives
!
    if (lgrad2) then
      d2exp1zz = dexp1*kvec*kvec
      d2exp2zz = dexp2*kvec*kvec
      d2exp1zg = dexp1*(zd + rkvec)
      d2exp2zg = dexp2*(zd - rkvec)
      d2exp1gg = dexp1*zd*rkvec*rkvec*(zd - rkvec)
      d2exp2gg = dexp2*zd*rkvec*rkvec*(zd + rkvec)
      d2erfc1zz = - 2.0_dp*d1erfc1z*seta*darg1
      d2erfc2zz = + 2.0_dp*d1erfc2z*seta*darg2
      d2erfc1zg = - 2.0_dp*d1erfc1z*rhseta*darg1*rkvec
      d2erfc2zg = - 2.0_dp*d1erfc2z*rhseta*darg2*rkvec
      d2erfc1gg = - d1erfc1g*rkvec*(rkvec + 2.0_dp*rhseta*darg1)
      d2erfc2gg = - d1erfc2g*rkvec*(rkvec + 2.0_dp*rhseta*darg2)
!
      f2zz = dexp1*d2erfc1zz + dexp2*d2erfc2zz + 2.0_dp*(d1exp1z*d1erfc1z + d1exp2z*d1erfc2z) + &
             d2exp1zz*derfc1 + d2exp2zz*derfc2
      f2zg = dexp1*d2erfc1zg + dexp2*d2erfc2zg + d1exp1g*d1erfc1z + d1exp2g*d1erfc2z + &
             d2exp1zg*derfc1 + d2exp2zg*derfc2 + d1exp1z*d1erfc1g + d1exp2z*d1erfc2g
      f2gg = dexp1*d2erfc1gg + dexp2*d2erfc2gg + 2.0_dp*(d1exp1g*d1erfc1g + d1exp2g*d1erfc2g) + &
             d2exp1gg*derfc1 + d2exp2gg*derfc2
!
!  Third derivatives
!
      if (lgrad3) then
        d3exp1zzz = dexp1*kvec*kvec*kvec
        d3exp2zzz = - dexp2*kvec*kvec*kvec
        d3exp1zzg = dexp1*(2.0_dp + zd*kvec)
        d3exp2zzg = dexp2*(2.0_dp - zd*kvec)
        d3exp1zgg = dexp1*rkvec*(zd*rkvec - rkvec*rkvec + zd*zd)
        d3exp2zgg = dexp2*rkvec*(zd*rkvec + rkvec*rkvec - zd*zd)
        d3erfc1zzz = + seta*seta*d1erfc1z*(4.0_dp*darg1*darg1 - 2.0_dp)
        d3erfc2zzz = + seta*seta*d1erfc2z*(4.0_dp*darg2*darg2 - 2.0_dp)
        d3erfc1zzg = + seta*seta*d1erfc1g*(4.0_dp*darg1*darg1 - 2.0_dp)
        d3erfc2zzg = + seta*seta*d1erfc2g*(4.0_dp*darg2*darg2 - 2.0_dp)
        d3erfc1zgg = + seta*rkvec*d1erfc1g*(rhseta*(4.0_dp*darg1*darg1 - 2.0_dp) + 2.0_dp*rkvec*darg1)
        d3erfc2zgg = - seta*rkvec*d1erfc2g*(rhseta*(4.0_dp*darg2*darg2 - 2.0_dp) + 2.0_dp*rkvec*darg2)
!
        f3zzz = dexp1*d3erfc1zzz + dexp2*d3erfc2zzz + 3.0_dp*(d1exp1z*d2erfc1zz + d1exp2z*d2erfc2zz) + &
                d3exp1zzz*derfc1 + d3exp2zzz*derfc2 + 3.0_dp*(d2exp1zz*d1erfc1z + d2exp2zz*d1erfc2z)
        f3zzg = dexp1*d3erfc1zzg + dexp2*d3erfc2zzg + 2.0_dp*(d1exp1z*d2erfc1zg + d1exp2z*d2erfc2zg) + &
                d3exp1zzg*derfc1 + d3exp2zzg*derfc2 + 2.0_dp*(d2exp1zg*d1erfc1z + d2exp2zg*d1erfc2z) + &
                d1exp1g*d2erfc1zz + d1exp2g*d2erfc2zz + d2exp1zz*d1erfc1g + d2exp2zz*d1erfc2g
        f3zgg = dexp1*d3erfc1zgg + dexp2*d3erfc2zgg + 2.0_dp*(d1exp1g*d2erfc1zg + d1exp2g*d2erfc2zg) + &
                d3exp1zgg*derfc1 + d3exp2zgg*derfc2 + 2.0_dp*(d2exp1zg*d1erfc1g + d2exp2zg*d1erfc2g) + &
                d1exp1z*d2erfc1gg + d1exp2z*d2erfc2gg + d2exp1gg*d1erfc1z + d2exp2gg*d1erfc2z
      endif
    endif
  endif
#ifdef TRACE
  call trace_out('getf2D')
#endif
!
  return
  end
