  subroutine psibaskes(nati,ntypi,neamspeci,natj,ntypj,neamspecj,npot,r,rk,psfct, &
                       eatm,deriv,deriv2,lgrad1,lgrad2)
!
!  Subroutine for calculating psi(R) in MEAM twobody potential
!
!  On return deriv and deriv2 contain the complete first and
!  second derivative terms.
!
!   8/14 Created from twobody
!   8/14 Modified to use sum or exponential density expressions
!   8/14 Contribution of second neighbours in reference density
!        now allowed for based on Z2S term
!   8/14 Trap added for small densities
!   8/14 Atom types for pair potential added
!   8/14 Density parameters now taken from density arrays
!   8/14 neamspec for i and j now passed in 
!   8/14 rho0 replaced by eamfnpar(3,)
!   8/14 rho now calculated using call to psirho
!   8/14 Weighting of different MEAM functional contributions added
!   2/18 Trace added
!   1/21 Trapping of small densities corrected for second derivatives
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
!  Copyright Curtin University 2021
!
!  Julian Gale, CIC, Curtin University, January 2021
!
  use eam
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  use two
  implicit none
!
!  Passed variables
!
  integer(i4),    intent(in)      :: nati
  integer(i4),    intent(in)      :: natj
  integer(i4),    intent(in)      :: neamspeci
  integer(i4),    intent(in)      :: neamspecj
  integer(i4),    intent(in)      :: ntypi
  integer(i4),    intent(in)      :: ntypj
  integer(i4),    intent(in)      :: npot
  logical,        intent(in)      :: lgrad1
  logical,        intent(in)      :: lgrad2
  real(dp),       intent(in)      :: r      ! Distance
  real(dp),       intent(in)      :: rk     ! Inverse distance
  real(dp),       intent(in)      :: psfct  ! Scale factor
  real(dp),       intent(out)     :: eatm   ! Energy of potential
  real(dp),       intent(out)     :: deriv  ! First derivative of potential divided by r
  real(dp),       intent(out)     :: deriv2 ! Second derivative of potential 
!
!  Local variables
!
  logical                         :: lmatch
  logical                         :: lswap
  real(dp)                        :: bpt
  real(dp)                        :: d
  real(dp)                        :: Eci
  real(dp)                        :: Ecj
  real(dp)                        :: Ecij
  real(dp)                        :: explam
  real(dp)                        :: exptrm
  real(dp)                        :: drhoidr
  real(dp)                        :: drhojdr
  real(dp)                        :: d2rhoidr2
  real(dp)                        :: d2rhojdr2
  real(dp)                        :: gam
  real(dp)                        :: lam
  real(dp)                        :: logtrmi
  real(dp)                        :: logtrmj
  real(dp)                        :: r0
  real(dp)                        :: rho_threshold
  real(dp)                        :: rhoi
  real(dp)                        :: rhoj
  real(dp)                        :: rho0i
  real(dp)                        :: rho0j
  real(dp)                        :: rrhoi
  real(dp)                        :: rrhoj
  real(dp)                        :: trm0
  real(dp)                        :: twoZ
  real(dp)                        :: dtrm0dr
  real(dp)                        :: strm
  real(dp)                        :: dstrmdr
  real(dp)                        :: d2strmdr2
  real(dp)                        :: xi
  real(dp)                        :: xj
  real(dp)                        :: dxidr
  real(dp)                        :: dxjdr
  real(dp)                        :: d2xidr2
  real(dp)                        :: d2xjdr2
  real(dp)                        :: sca
  real(dp)                        :: wi
  real(dp)                        :: wj
  real(dp)                        :: Z2S
#ifdef TRACE
  call trace_in('psibaskes')
#endif
!
!  Check potential types
!
  if (nptype(npot).ne.45.and.nptype(npot).ne.55) then
    eatm = 0.0_dp
    if (lgrad1) then
      deriv = 0.0_dp
      if (lgrad2) then
        deriv2 = 0.0_dp
      endif
    endif
#ifdef TRACE
    call trace_out('psibaskes')
#endif
    return
  endif
!
!  Set threshold for neglecting density contributions
!
  rho_threshold = 1.0d-32
!
!  1NN terms
!
  twoZ = - 2.0_dp*psfct/twopot(5,npot) ! - (2/Z)
  Ecij = twopot(1,npot)  ! Ec for i-j pair
  Eci  = eamfnpar(1,neamspeci)*eamfnpar(2,neamspeci)  ! Ec*A for i
  Ecj  = eamfnpar(1,neamspecj)*eamfnpar(2,neamspecj)  ! Ec*A for j
  bpt  = twopot(3,npot)  ! alpha
  r0   = twopot(4,npot)  ! r0
  d    = twopot(7,npot)  ! d
  gam  = twopot(8,npot)  ! gamma
  lam  = twopot(9,npot)  ! lambda
!
!  Extra 2NN terms
!
  Z2S = twopot(10,npot)  ! Z2 x S => number of second neighbours multiplied by screening factor
  sca = twopot(11,npot)  ! Scale factor to convert first neighbour distance to second neighbour one
!
!  Find weight of functionals taking care with regard to the order
!
  if (lmatch(nati,ntypi,nspec1(npot),nptyp1(npot),.true.).and. &
      lmatch(natj,ntypj,nspec2(npot),nptyp2(npot),.true.)) then
    wi = twopot(6,npot)
    wj = 1.0_dp - wi
    lswap = .false.
  else
    wj = twopot(6,npot)
    wi = 1.0_dp - wj
    lswap = .true.
  endif
!************************************
!  Compute the reference densities  *
!************************************
  call psirho(nati,ntypi,neamspeci,natj,ntypj,neamspecj,npot,r,rhoi,rhoj,drhoidr,drhojdr, &
              d2rhoidr2,d2rhojdr2,lgrad1,lgrad2)
  if (rhoi.gt.rho_threshold) then
    rho0i = eamfnpar(3,neamspeci)
    xi = (rhoi/rho0i)
    logtrmi = log(xi)
    if (lgrad1) then
      rrhoi = 1.0_dp/rhoi
      dxidr = 0.5_dp*rrhoi*drhoidr/rho0i
      if (lgrad2) then
        d2xidr2 = 0.5_dp*rrhoi*(d2rhoidr2 - 0.5_dp*rrhoi*rrhoi*drhoidr*drhoidr)/rho0i
      endif
    endif
  else
    xi = 0.0_dp
    logtrmi = 0.0_dp
    if (lgrad1) then
      dxidr = 0.0_dp
      if (lgrad2) then
        d2xidr2 = 0.0_dp
      endif
    endif
  endif
  if (rhoj.gt.rho_threshold) then
    rho0j = eamfnpar(3,neamspecj)
    xj = (rhoj/rho0j)
    logtrmj = log(xj)
    if (lgrad1) then
      rrhoj = 1.0_dp/rhoj
      dxjdr = 0.5_dp*rrhoj*drhojdr/rho0j
      if (lgrad2) then
        d2xjdr2 = 0.5_dp*rrhoj*(d2rhojdr2 - 0.5_dp*rrhoj*rrhoj*drhojdr*drhojdr)/rho0j
      endif
    endif
  else
    xj = 0.0_dp
    logtrmj = 0.0_dp
    if (lgrad1) then
      dxjdr = 0.0_dp
      if (lgrad2) then
        d2xjdr2 = 0.0_dp
      endif
    endif
  endif
!
!  Trap small densities that might cause numerical instability
!
  if ((rhoi+rhoj).lt.rho_threshold) then
#ifdef TRACE
    call trace_out('psibaskes')
#endif
    return
  endif
!
!  Common terms
!
  dtrm0dr = bpt/r0
  trm0 = dtrm0dr*(r - r0)
  exptrm = exp(-trm0)
!
  if (nptype(npot).eq.45) then
!***************
!  Baskes pot  *
!***************
    if (abs(d).gt.1.0d-12) then
      strm = (1.0_dp + trm0*(1.0_dp + (trm0**2)*d))
      eatm = twoZ*Ecij*strm*exptrm
      if (rhoi.gt.rho_threshold) then
        eatm = eatm + twoZ*Eci*wi*xi*logtrmi
      endif
      if (rhoj.gt.rho_threshold) then
        eatm = eatm + twoZ*Ecj*wj*xj*logtrmj
      endif
      if (lgrad1) then
        dstrmdr = (1.0_dp + (trm0**2)*(3.0_dp*d))*dtrm0dr
        deriv = twoZ*Ecij*(dstrmdr - strm*dtrm0dr)*exptrm
        if (rhoi.gt.rho_threshold) then
          deriv = deriv + twoZ*Eci*wi*dxidr*(1.0_dp + logtrmi)
        endif
        if (rhoj.gt.rho_threshold) then
          deriv = deriv + twoZ*Ecj*wj*dxjdr*(1.0_dp + logtrmj)
        endif
        if (lgrad2) then
          d2strmdr2 = (6.0_dp*trm0*d)*dtrm0dr*dtrm0dr
          deriv2 = twoZ*Ecij*exptrm*(d2strmdr2 - 2.0_dp*dstrmdr*dtrm0dr + strm*dtrm0dr**2)
          if (rhoi.gt.rho_threshold) then
            deriv2 = deriv2 + twoZ*Eci*wi*(d2xidr2*(1.0_dp + logtrmi) + dxidr*dxidr/xi)
          endif
          if (rhoj.gt.rho_threshold) then
            deriv2 = deriv2 + twoZ*Ecj*wj*(d2xjdr2*(1.0_dp + logtrmj) + dxjdr*dxjdr/xj)
          endif
        endif
      endif
    else
      strm = (1.0_dp + trm0)
      eatm = twoZ*Ecij*strm*exptrm
      if (rhoi.gt.rho_threshold) then
        eatm = eatm + twoZ*Eci*wi*xi*logtrmi
      endif
      if (rhoj.gt.rho_threshold) then
        eatm = eatm + twoZ*Ecj*wj*xj*logtrmj
      endif
      if (lgrad1) then
        dstrmdr = dtrm0dr
        deriv = twoZ*Ecij*(dstrmdr - strm*dtrm0dr)*exptrm
        if (rhoi.gt.rho_threshold) then
          deriv = deriv + twoZ*Eci*wi*dxidr*(1.0_dp + logtrmi)
        endif
        if (rhoj.gt.rho_threshold) then
          deriv = deriv + twoZ*Ecj*wj*dxjdr*(1.0_dp + logtrmj)
        endif
        if (lgrad2) then
          d2strmdr2 = 0.0_dp
          deriv2 = twoZ*Ecij*exptrm*(d2strmdr2 - 2.0_dp*dstrmdr*dtrm0dr + strm*dtrm0dr**2)
          if (rhoi.gt.rho_threshold) then
            deriv2 = deriv2 + twoZ*Eci*wi*(d2xidr2*(1.0_dp + logtrmi) + dxidr*dxidr/xi)
          endif
          if (rhoj.gt.rho_threshold) then
            deriv2 = deriv2 + twoZ*Ecj*wj*(d2xjdr2*(1.0_dp + logtrmj) + dxjdr*dxjdr/xj)
          endif
        endif
      endif
    endif
  elseif (nptype(npot).eq.55) then
!******************
!  Baskes a4 pot  *
!******************
    explam = exp(-lam*trm0**2)
    if (abs(d).gt.1.0d-12) then
      strm = (1.0_dp + trm0*(1.0_dp + (trm0**2)*(d + gam*trm0*explam*rk)))
      eatm = twoZ*Ecij*strm*exptrm
      if (rhoi.gt.rho_threshold) then
        eatm = eatm + twoZ*Eci*wi*xi*logtrmi
      endif
      if (rhoj.gt.rho_threshold) then
        eatm = eatm + twoZ*Ecj*wj*xj*logtrmj
      endif
      if (lgrad1) then
        dstrmdr = (1.0_dp + (trm0**2)*(3.0_dp*d + gam*trm0*explam*rk*(4.0_dp - 2.0_dp*lam*trm0**2)))*dtrm0dr - &
                  gam*(trm0**4)*explam*rk**2
        deriv = twoZ*Ecij*(dstrmdr - strm*dtrm0dr)*exptrm
        if (rhoi.gt.rho_threshold) then
          deriv = deriv + twoZ*Eci*wi*dxidr*(1.0_dp + logtrmi)
        endif
        if (rhoj.gt.rho_threshold) then
          deriv = deriv + twoZ*Ecj*wj*dxjdr*(1.0_dp + logtrmj)
        endif
        if (lgrad2) then
          d2strmdr2 = (dtrm0dr**2)*(6.0_dp*d*trm0 + gam*explam*(trm0**2)*rk*(12.0_dp - 18.0_dp*lam*trm0**2 + &
                                    4.0_dp*lam*lam*trm0**4)) + &
                      dtrm0dr*gam*rk*rk*explam*(trm0**3)*(6.0_dp*lam*trm0**2 - 4.0_dp) + &
                      2.0_dp*explam*gam*rk*rk*rk*trm0**4
          deriv2 = twoZ*Ecij*exptrm*(d2strmdr2 - 2.0_dp*dstrmdr*dtrm0dr + strm*dtrm0dr**2)
          if (rhoi.gt.rho_threshold) then
            deriv2 = deriv2 + twoZ*Eci*wi*(d2xidr2*(1.0_dp + logtrmi) + dxidr*dxidr/xi)
          endif
          if (rhoj.gt.rho_threshold) then
            deriv2 = deriv2 + twoZ*Ecj*wj*(d2xjdr2*(1.0_dp + logtrmj) + dxjdr*dxjdr/xj)
          endif
        endif
      endif
    else
      strm = (1.0_dp + trm0*(1.0_dp + (trm0**3)*gam*explam*rk))
      eatm = twoZ*Ecij*strm*exptrm
      if (rhoi.gt.rho_threshold) then
        eatm = eatm + twoZ*Eci*wi*xi*logtrmi
      endif
      if (rhoj.gt.rho_threshold) then
        eatm = eatm + twoZ*Ecj*wj*xj*logtrmj
      endif
      if (lgrad1) then
        dstrmdr = (1.0_dp + (trm0**2)*(gam*trm0*explam*rk*(4.0_dp - 2.0_dp*lam*trm0**2)))*dtrm0dr - &
                  gam*(trm0**4)*explam*rk**2
        deriv = twoZ*Ecij*(dstrmdr - strm*dtrm0dr)*exptrm
        if (rhoi.gt.rho_threshold) then
          deriv = deriv + twoZ*Eci*wi*dxidr*(1.0_dp + logtrmi)
        endif
        if (rhoj.gt.rho_threshold) then
          deriv = deriv + twoZ*Ecj*wj*dxjdr*(1.0_dp + logtrmj)
        endif
        if (lgrad2) then
          d2strmdr2 = (dtrm0dr**2)*(gam*explam*(trm0**2)*rk*(12.0_dp - 18.0_dp*lam*trm0**2 + &
                                    4.0_dp*lam*lam*trm0**4)) + &
                      dtrm0dr*gam*rk*rk*explam*(trm0**3)*(6.0_dp*lam*trm0**2 - 4.0_dp) + &
                      2.0_dp*explam*gam*rk*rk*rk*trm0**4
          deriv2 = twoZ*Ecij*exptrm*(d2strmdr2 - 2.0_dp*dstrmdr*dtrm0dr + strm*dtrm0dr**2)
          if (rhoi.gt.rho_threshold) then
            deriv2 = deriv2 + twoZ*Eci*wi*(d2xidr2*(1.0_dp + logtrmi) + dxidr*dxidr/xi)
          endif
          if (rhoj.gt.rho_threshold) then
            deriv2 = deriv2 + twoZ*Ecj*wj*(d2xjdr2*(1.0_dp + logtrmj) + dxjdr*dxjdr/xj)
          endif
        endif
      endif
    endif
  endif
#ifdef TRACE
  call trace_out('psibaskes')
#endif
!
  return
  end
