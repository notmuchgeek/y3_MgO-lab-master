  subroutine EDIP_threebody(nspeci,nspecj,nspeck,xij,yij,zij,rij,xik,yik,zik,rik,Zi,e3i, &
                            de3idr,de3idZi,d2e3idr2,d2e3idrdZi,d2e3idZi2,lgrad1,lgrad2)
!
!  Subroutine to calculate the threebody term for the EDIP potential
!
!  On entry : 
!
!  nspeci          = species number of atom type i 
!  nspecj          = species number of atom type j 
!  nspeck          = species number of atom type k 
!  xij             = x component of i->j vector
!  yij             = y component of i->j vector
!  zij             = z component of i->j vector
!  rij             = length of i->j vector
!  xik             = x component of i->k vector
!  yik             = y component of i->k vector
!  zik             = z component of i->k vector
!  rik             = length of i->k vector
!  Zi              = coordination number of pivot atom i
!  lgrad1          = if .true. calculate the first derivatives
!  lgrad2          = if .true. calculate the second derivatives
!
!  On exit :
!
!  e3i             = the value of the threebody energy
!  de3idr(3)       = the first derivatives of e3i w.r.t. to the three different interatomic vectors
!  de3idZi         = the first derivative of e3i w.r.t. to Zi
!  d2e3idr2(6)     = the second derivatives of e3i w.r.t. to the three different interatomic vectors
!  d2e3idrdZi(3)   = the second derivatives of e3i w.r.t. to the three different interatomic vectors and Zi
!  d2e3idZi2       = the second derivative of e3i w.r.t. to Zi twice
!
!  10/10 Created from gtheta
!  10/10 Cutoffs added for exponential terms
!  10/10 EDIP linear threebody modifications added
!   7/13 Call to EDIP_expfn corrected so that Zi is removed as an argument
!   1/14 Modified to allow for original form of 3-body interaction as well
!        as Nigel's newer one. 
!   7/14 Modified so that kq term becomes an option to allow backwards
!        compatibility
!  10/14 Correction of muZ to - mu in original expression for dhargdZi
!  10/14 Correction to derivatives of lambda for original formulation
!   1/18 Second derivatives added
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
  use edipdata
  use g_constants,    only : pi
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)             :: nspeci
  integer(i4), intent(in)             :: nspecj
  integer(i4), intent(in)             :: nspeck
  real(dp),    intent(in)             :: Zi
  real(dp),    intent(in)             :: rij
  real(dp),    intent(in)             :: rik
  real(dp),    intent(in)             :: xij
  real(dp),    intent(in)             :: yij
  real(dp),    intent(in)             :: zij
  real(dp),    intent(in)             :: xik
  real(dp),    intent(in)             :: yik
  real(dp),    intent(in)             :: zik
  real(dp),    intent(out)            :: e3i
  real(dp),    intent(out)            :: de3idr(3)
  real(dp),    intent(out)            :: de3idZi
  real(dp),    intent(out)            :: d2e3idr2(6)
  real(dp),    intent(out)            :: d2e3idrdZi(3)
  real(dp),    intent(out)            :: d2e3idZi2
  logical,     intent(in)             :: lgrad1
  logical,     intent(in)             :: lgrad2
!
!  Local variables
!
  integer(i4)                         :: indij
  integer(i4)                         :: indik
  integer(i4)                         :: indjk
  real(dp)                            :: cZ2
  real(dp)                            :: dcZ2dZ
  real(dp)                            :: d2cZ2dZ2
  real(dp)                            :: costheta
  real(dp)                            :: cos1d(3)
  real(dp)                            :: cos2d(6)
  real(dp)                            :: cutoffij
  real(dp)                            :: cutoffij2
  real(dp)                            :: cutoffik
  real(dp)                            :: cutoffik2
  real(dp)                            :: dgijdrij
  real(dp)                            :: dgikdrik
  real(dp)                            :: dgijdZi
  real(dp)                            :: dgikdZi
  real(dp)                            :: d2gijdrij2
  real(dp)                            :: d2gijdrijdZi
  real(dp)                            :: d2gijdZi2
  real(dp)                            :: d2gikdrik2
  real(dp)                            :: d2gikdrikdZi
  real(dp)                            :: d2gikdZi2
  real(dp)                            :: dharg
  real(dp)                            :: d2harg
  real(dp)                            :: dhargdZi
  real(dp)                            :: dhijkdcostheta
  real(dp)                            :: dhijkdcosthetaA
  real(dp)                            :: dhijkdcosthetaB
  real(dp)                            :: dhijkdZi
  real(dp)                            :: dhijkdZiA
  real(dp)                            :: dhijkdZiB
  real(dp)                            :: d2hijkdcostheta2
  real(dp)                            :: d2hijkdcosthetadZi
  real(dp)                            :: d2hijkdZi2
  real(dp)                            :: d2hijkdcostheta2B
  real(dp)                            :: d2hijkdcosthetadZiB
  real(dp)                            :: d2hijkdZi2B
  real(dp)                            :: dlambdadZi
  real(dp)                            :: d2lambdadZi2
  real(dp)                            :: dtauZdZ
  real(dp)                            :: d2tauZdZ2
  real(dp)                            :: eta
  real(dp)                            :: exph
  real(dp)                            :: expharg
  real(dp)                            :: exphargeta
  real(dp)                            :: gij
  real(dp)                            :: gik
  real(dp)                            :: harg
  real(dp)                            :: hijk
  real(dp)                            :: hijkA
  real(dp)                            :: hijkB
  real(dp)                            :: lambda
  real(dp)                            :: mu
  real(dp)                            :: muZ
  real(dp)                            :: Qh
  real(dp)                            :: kQh2
  real(dp)                            :: rQh
  real(dp)                            :: rrij
  real(dp)                            :: rrik
  real(dp)                            :: rrij2
  real(dp)                            :: rrik2
  real(dp)                            :: rrij4
  real(dp)                            :: rrik4
  real(dp)                            :: r2ij
  real(dp)                            :: r2ik
  real(dp)                            :: r2jk
  real(dp)                            :: tauZ
  real(dp)                            :: x
  real(dp)                            :: xjk
  real(dp)                            :: yjk
  real(dp)                            :: zjk
  real(dp)                            :: Z0
#ifdef TRACE
  call trace_in('edip_threebody')
#endif
!
!  Compute index of i-j pair
!
  if (nspeci.ge.nspecj) then
    indij = nspeci*(nspeci - 1) + nspecj
  else
    indij = nspecj*(nspecj - 1) + nspeci
  endif
!
!  Compute index of i-k pair
!
  if (nspeci.ge.nspeck) then
    indik = nspeci*(nspeci - 1) + nspeck
  else
    indik = nspeck*(nspeck - 1) + nspeci
  endif
!
!  Compute index of j-k pair
!
  if (nspecj.ge.nspeck) then
    indjk = nspecj*(nspecj - 1) + nspeck
  else
    indjk = nspeck*(nspeck - 1) + nspecj
  endif
!*******************************
!  Compute and check cut-offs  *
!*******************************
  cutoffij = EDIP2a(indij) + EDIP2aprime(indij)*Zi
  cutoffik = EDIP2a(indik) + EDIP2aprime(indik)*Zi
  cutoffij2 = cutoffij + EDIP3gamma0(indjk,nspeci)*EDIPaccuracy2drmax
  cutoffik2 = cutoffik + EDIP3gamma0(indjk,nspeci)*EDIPaccuracy2drmax
!
  if (rij.lt.cutoffij2.and.rik.lt.cutoffik2) then
!
!  Calculate j->k vector
!
    xjk = xik - xij
    yjk = yik - yij
    zjk = zik - zij
!
!  Calculate interatomic distances
!
    r2ij = xij*xij + yij*yij + zij*zij
    r2ik = xik*xik + yik*yik + zik*zik
    r2jk = xjk*xjk + yjk*yjk + zjk*zjk
!
!  Calculate reciprocal distances
!
    rrij = 1.0_dp/rij
    rrik = 1.0_dp/rik
!
!  Calculate cos(theta) using cosine rule
!
    costheta = 0.5_dp*(r2ij + r2ik - r2jk)*rrij*rrik
!**********************
!  Compute lambda(Z)  *
!**********************
    if (lEDIP3orig(indjk,nspeci)) then
!
!  Original
!
      lambda = EDIP3lambda0(indjk,nspeci)
    else
!
!  New
!
      Z0 = EDIP3Z0(indjk,nspeci)
      lambda = EDIP3lambda0(indjk,nspeci)*exp(-EDIP3lambdap(indjk,nspeci)*(Zi-Z0)**2)
    endif
!*******************
!  Compute g(r,Z)  *
!*******************
    call EDIP_expfn(rij,cutoffij,EDIP2aprime(indij),EDIP3gamma0(indjk,nspeci),gij,dgijdrij,dgijdZi, &
                    d2gijdrij2,d2gijdrijdZi,d2gijdZi2,lgrad1,lgrad2)
    call EDIP_expfn(rik,cutoffik,EDIP2aprime(indik),EDIP3gamma0(indjk,nspeci),gik,dgikdrik,dgikdZi, &
                    d2gikdrik2,d2gikdrikdZi,d2gikdZi2,lgrad1,lgrad2)
    gij = gij*EDIP3gammap(indjk,nspeci)
    gik = gik*EDIP3gammap(indjk,nspeci)
!*******************
!  Compute tau(Z)  *
!*******************
    call EDIP_tau(Zi,tauZ,dtauZdZ,d2tauZdZ2,lgrad1,lgrad2)
!***********************
!  Compute h(theta,Z)  *
!***********************
    if (lEDIP3orig(indjk,nspeci)) then
!
!  Original
!
      eta = EDIP3Z0(indjk,nspeci)
      mu = EDIP3lambdap(indjk,nspeci)
      muZ = - mu*Zi
      Qh = EDIP3q(indjk,nspeci)*exp(muZ)
      harg = Qh*(costheta + tauZ)**2
      expharg = exp(-harg)
      hijk = 1.0_dp - expharg + eta*harg
      if (lgrad1) then
        dharg = 2.0_dp*Qh*(costheta + tauZ)
        exphargeta = expharg + eta
        dhijkdcostheta = dharg*exphargeta
!
!  Following line has been corrected muZ -> - mu
!
        dhargdZi = - mu*harg + dharg*dtauZdZ
        dhijkdZi = dhargdZi*exphargeta
        if (lgrad2) then
          d2harg = 2.0_dp*Qh
          d2hijkdcostheta2   = exphargeta*d2harg
          d2hijkdcosthetadZi = exphargeta*(d2harg*dtauZdZ - mu*dharg) - &
                               expharg*dhargdZi*dharg
          d2hijkdZi2         = exphargeta*(dharg*d2tauZdZ2 + d2harg*dtauZdZ**2 - &
                                           2.0_dp*mu*dharg*dtauZdZ + mu*mu*harg) - &
                               expharg*dhargdZi*dhargdZi
        endif
      endif
    elseif (lEDIP3mod(indjk,nspeci)) then
!
!  Modified new form
!
      if (Zi.ge.3) then
        Qh = EDIP3q(indjk,nspeci)
        rQh = 1.0_dp/Qh
        exph = exp(-Qh*(costheta + tauZ)**2)
        hijk = rQh*(1.0_dp - exph)
        if (lgrad1) then
          dhijkdcostheta = 2.0_dp*exph*(costheta + tauZ)
          dhijkdZi = dhijkdcostheta*dtauZdZ
          if (lgrad2) then
            d2hijkdcostheta2   = exph*(2.0_dp - 4.0_dp*Qh*(costheta + tauZ)**2)
            d2hijkdcosthetadZi = d2hijkdcostheta2*dtauZdZ
            d2hijkdZi2         = d2hijkdcosthetadZi*dtauZdZ + dhijkdcostheta*d2tauZdZ2
          endif
        endif
      elseif (Zi.le.2) then
        kQh2 = EDIP3kq2(indjk,nspeci)
        hijk = kQh2*(1.0_dp + costheta)
        if (lgrad1) then
          dhijkdcostheta = kQh2
          dhijkdZi = 0.0_dp
          if (lgrad2) then
            d2hijkdcostheta2   = 0.0_dp
            d2hijkdcosthetadZi = 0.0_dp
            d2hijkdZi2         = 0.0_dp
          endif
        endif
      else
        Qh = EDIP3q(indjk,nspeci)
        rQh = 1.0_dp/Qh
        exph = exp(-Qh*(costheta + tauZ)**2)
!
        x = (Zi - 2.0_dp)
        cZ2 = 0.5_dp*(1.0_dp + cos(pi*x))
!
        kQh2 = EDIP3kq2(indjk,nspeci)
        hijkA = kQh2*(1.0_dp + costheta)
        hijkB = rQh*(1.0_dp - exph)
        hijk = cZ2*hijkA + (1.0_dp - cZ2)*hijkB
        if (lgrad1) then
          dcZ2dZ = - 0.5_dp*pi*sin(pi*x)
!
          dhijkdcosthetaA = kQh2
          dhijkdZiA = 0.0_dp
!
          dhijkdcosthetaB = 2.0_dp*exph*(costheta + tauZ)
          dhijkdZiB = dhijkdcosthetaB*dtauZdZ
!
          dhijkdcostheta = cZ2*dhijkdcosthetaA + (1.0_dp - cZ2)*dhijkdcosthetaB
          dhijkdZi = cZ2*dhijkdZiA + (1.0_dp - cZ2)*dhijkdZiB + &
                     dcZ2dZ*(hijkA - hijkB)
          if (lgrad2) then
            d2cZ2dZ2 = - 0.5_dp*pi*pi*cos(pi*x)
!
            d2hijkdcostheta2B   = exph*(2.0_dp - 4.0_dp*Qh*(costheta + tauZ)**2)
            d2hijkdcosthetadZiB = d2hijkdcostheta2B*dtauZdZ
            d2hijkdZi2B         = d2hijkdcosthetadZiB*dtauZdZ + dhijkdcosthetaB*d2tauZdZ2
!
            d2hijkdcostheta2   = (1.0_dp - cZ2)*d2hijkdcostheta2B
            d2hijkdcosthetadZi = (1.0_dp - cZ2)*d2hijkdcosthetadZiB + dcZ2dZ*(dhijkdcosthetaA - dhijkdcosthetaB)
            d2hijkdZi2         = (1.0_dp - cZ2)*d2hijkdZi2B + 2.0_dp*dcZ2dZ*(dhijkdcosthetaA - dhijkdcosthetaB) + &
                                 d2cZ2dZ2*(hijkA - hijkB)
          endif
        endif
      endif
    else
!
!  New without kq
!
      Qh = EDIP3q(indjk,nspeci)
      rQh = 1.0_dp/Qh
      exph = exp(-Qh*(costheta + tauZ)**2)
      hijk = rQh*(1.0_dp - exph)
      if (lgrad1) then
        dhijkdcostheta = 2.0_dp*exph*(costheta + tauZ)
        dhijkdZi = dhijkdcostheta*dtauZdZ
        if (lgrad2) then
          d2hijkdcostheta2   = exph*(2.0_dp - 4.0_dp*Qh*(costheta + tauZ)**2)
          d2hijkdcosthetadZi = d2hijkdcostheta2*dtauZdZ
          d2hijkdZi2         = d2hijkdcosthetadZi*dtauZdZ + dhijkdcostheta*d2tauZdZ2
        endif
      endif
    endif
!**************************************
!  Compute total energy contribution  *
!**************************************
    e3i = lambda*gij*gik*hijk
!**************************
!  Calculate derivatives  *
!**************************
    if (lgrad1) then
!
!  Initialise first derivatives of e3i
!
      de3idr(1:3) = 0.0_dp
      de3idZi = 0.0_dp
!
!  First derivatives of cos(theta)
!
!  1 = ij
!  2 = ik
!  3 = jk
!
      rrij2 = rrij*rrij
      rrik2 = rrik*rrik
      cos1d(1) = rrij*rrik - costheta*rrij2
      cos1d(2) = rrij*rrik - costheta*rrik2
      cos1d(3) = - rrij*rrik
      if (lgrad2) then
!
!  Second derivatives of cos(theta)
!
!  1 = 11
!  2 = 21
!  3 = 31
!  4 = 22
!  5 = 32
!  6 = 33
!
        rrij4 = rrij2*rrij2
        rrik4 = rrik2*rrik2
        cos2d(1) = - 2.0_dp*rrij2*rrij*rrik + 3.0_dp*costheta*rrij4
        cos2d(2) = costheta*rrij2*rrik2 - rrij*rrik*(rrij2+rrik2)
        cos2d(3) = rrij2*rrij*rrik
        cos2d(4) = - 2.0_dp*rrik2*rrik*rrij + 3.0_dp*costheta*rrik4
        cos2d(5) = rrik2*rrij*rrik
        cos2d(6) = 0.0_dp
      endif
!
!  Derivatives of g(r,Z)
!
      dgijdrij = dgijdrij*EDIP3gammap(indjk,nspeci)
      dgikdrik = dgikdrik*EDIP3gammap(indjk,nspeci)
      dgijdZi = dgijdZi*EDIP3gammap(indjk,nspeci)
      dgikdZi = dgikdZi*EDIP3gammap(indjk,nspeci)
!
      de3idZi = de3idZi + lambda*gik*hijk*dgijdZi
      de3idZi = de3idZi + lambda*gij*hijk*dgikdZi
      de3idr(1) = de3idr(1) + lambda*gik*hijk*dgijdrij*rrij
      de3idr(2) = de3idr(2) + lambda*gij*hijk*dgikdrik*rrik
!
!  Derivatives of h(theta,Z)
!
      de3idr(1:3) = de3idr(1:3) + lambda*gij*gik*dhijkdcostheta*cos1d(1:3)
      de3idZi = de3idZi + lambda*gij*gik*dhijkdZi
!
      if (.not.lEDIP3orig(indjk,nspeci)) then
!
!  Derivative of lambda - original method has zero derivative
!
        dlambdadZi = - 2.0_dp*EDIP3lambdap(indjk,nspeci)*(Zi - Z0)*lambda
        de3idZi = de3idZi + gij*gik*hijk*dlambdadZi
      endif
      if (lgrad2) then
!
!  Second derivatives
!
!  Initialise derivatives
!
        d2e3idr2(1:6) = 0.0_dp
        d2e3idrdZi(1:3) = 0.0_dp
        d2e3idZi2 = 0.0_dp
!
!  Multiple g terms by constants
!
        d2gijdrij2   = d2gijdrij2*EDIP3gammap(indjk,nspeci)
        d2gikdrik2   = d2gikdrik2*EDIP3gammap(indjk,nspeci)
        d2gijdrijdZi = d2gijdrijdZi*EDIP3gammap(indjk,nspeci)
        d2gikdrikdZi = d2gikdrikdZi*EDIP3gammap(indjk,nspeci)
        d2gijdZi2    = d2gijdZi2*EDIP3gammap(indjk,nspeci)
        d2gikdZi2    = d2gikdZi2*EDIP3gammap(indjk,nspeci)
!*****************************************
!  Distance-distance second derivatives  *
!*****************************************
!
!  Angle derivatives
!
        d2e3idr2(1:6) = d2e3idr2(1:6) + lambda*gij*gik*dhijkdcostheta*cos2d(1:6)
!
        d2e3idr2(1) = d2e3idr2(1) + lambda*gij*gik*d2hijkdcostheta2*cos1d(1)*cos1d(1)
        d2e3idr2(2) = d2e3idr2(2) + lambda*gij*gik*d2hijkdcostheta2*cos1d(2)*cos1d(1)
        d2e3idr2(3) = d2e3idr2(3) + lambda*gij*gik*d2hijkdcostheta2*cos1d(3)*cos1d(1)
        d2e3idr2(4) = d2e3idr2(4) + lambda*gij*gik*d2hijkdcostheta2*cos1d(2)*cos1d(2)
        d2e3idr2(5) = d2e3idr2(5) + lambda*gij*gik*d2hijkdcostheta2*cos1d(3)*cos1d(2)
        d2e3idr2(6) = d2e3idr2(6) + lambda*gij*gik*d2hijkdcostheta2*cos1d(3)*cos1d(3)
!
        d2e3idr2(1) = d2e3idr2(1) + 2.0_dp*lambda*dhijkdcostheta*cos1d(1)*dgijdrij*gik
        d2e3idr2(2) = d2e3idr2(2) + lambda*dhijkdcostheta*(cos1d(2)*dgijdrij*gik + cos1d(1)*gij*dgikdrik)
        d2e3idr2(3) = d2e3idr2(3) + lambda*dhijkdcostheta*cos1d(3)*dgijdrij*gik
        d2e3idr2(4) = d2e3idr2(4) + 2.0_dp*lambda*dhijkdcostheta*cos1d(2)*gij*dgikdrik
        d2e3idr2(5) = d2e3idr2(5) + lambda*dhijkdcostheta*cos1d(3)*gij*dgikdrik
!
!  IJ case
!
        d2e3idr2(1) = d2e3idr2(1) + lambda*d2gijdrij2*gik*hijk
!
!  IK case
!
        d2e3idr2(4) = d2e3idr2(4) + lambda*gij*d2gikdrik2*hijk
!
!  IJ-IK mixed case
!
        d2e3idr2(2) = d2e3idr2(2) + lambda*dgijdrij*dgikdrik*hijk
!*********************************************
!  Distance-coordination second derivatives  *
!*********************************************
!
!  Angle derivatives
!
        d2e3idrdZi(1:3) = d2e3idrdZi(1:3) + lambda*d2hijkdcosthetadZi*cos1d(1:3)*gij*gik
!
        d2e3idrdZi(1) = d2e3idrdZi(1) + lambda*dhijkdcostheta*cos1d(1)*(dgijdZi*gik + gij*dgikdZi)
        d2e3idrdZi(2) = d2e3idrdZi(2) + lambda*dhijkdcostheta*cos1d(2)*(dgijdZi*gik + gij*dgikdZi)
        d2e3idrdZi(3) = d2e3idrdZi(3) + lambda*dhijkdcostheta*cos1d(3)*(dgijdZi*gik + gij*dgikdZi)
!
!  IJ case
!
        d2e3idrdZi(1) = d2e3idrdZi(1) + lambda*dhijkdZi*dgijdrij*gik
        d2e3idrdZi(1) = d2e3idrdZi(1) + lambda*hijk*dgijdrij*dgikdZi
!
!  IK case
!
        d2e3idrdZi(2) = d2e3idrdZi(2) + lambda*dhijkdZi*gij*dgikdrik
        d2e3idrdZi(2) = d2e3idrdZi(2) + lambda*hijk*dgikdrik*dgijdZi
!*************************************************
!  Coordination-coordination second derivatives  *
!*************************************************
!
!  Angle derivatives
!
        d2e3idZi2 = d2e3idZi2 + lambda*d2hijkdZi2*gij*gik
        d2e3idZi2 = d2e3idZi2 + 2.0_dp*lambda*dhijkdZi*dgijdZi*gik
        d2e3idZi2 = d2e3idZi2 + 2.0_dp*lambda*dhijkdZi*gij*dgikdZi
!
        d2e3idZi2 = d2e3idZi2 + lambda*hijk*d2gijdZi2*gik
        d2e3idZi2 = d2e3idZi2 + 2.0_dp*lambda*hijk*dgijdZi*dgikdZi
        d2e3idZi2 = d2e3idZi2 + lambda*hijk*gij*d2gikdZi2
!
        if (.not.lEDIP3orig(indjk,nspeci)) then
!
!  Second derivative of lambda - original method has zero derivative
!
          d2lambdadZi2 = lambda*(4.0_dp*(EDIP3lambdap(indjk,nspeci)*(Zi - Z0))**2 - &
                                 2.0_dp*EDIP3lambdap(indjk,nspeci))
!
          d2e3idrdZi(1) = d2e3idrdZi(1) + gij*gik*dhijkdcostheta*cos1d(1)*dlambdadZi
          d2e3idrdZi(2) = d2e3idrdZi(2) + gij*gik*dhijkdcostheta*cos1d(2)*dlambdadZi
          d2e3idrdZi(3) = d2e3idrdZi(3) + gij*gik*dhijkdcostheta*cos1d(3)*dlambdadZi
!
          d2e3idrdZi(1) = d2e3idrdZi(1) + dgijdrij*gik*dlambdadZi
          d2e3idrdZi(2) = d2e3idrdZi(2) + gij*dgikdrik*dlambdadZi
!
          d2e3idZi2 = d2e3idZi2 + gij*gik*hijk*d2lambdadZi2
          d2e3idZi2 = d2e3idZi2 + gij*gik*dhijkdZi*dlambdadZi
          d2e3idZi2 = d2e3idZi2 + dgijdZi*gik*hijk*dlambdadZi
          d2e3idZi2 = d2e3idZi2 + gij*dgikdZi*hijk*dlambdadZi
        endif
      endif
    endif
  else
    e3i = 0.0_dp
    if (lgrad1) then
      de3idr(1:3) = 0.0_dp
      de3idZi = 0.0_dp
      if (lgrad2) then
        d2e3idr2(1:6) = 0.0_dp
        d2e3idrdZi(1:3) = 0.0_dp
        d2e3idZi2 = 0.0_dp
      endif
    endif
  endif
#ifdef TRACE
  call trace_out('edip_threebody')
#endif
!
  return
  end
!
  subroutine EDIP_tau(Z,tauZ,dtauZdZ,d2tauZdZ2,lgrad1,lgrad2)
!
!  Subroutine to calculate the tau term for the EDIP potential
!
!  On entry : 
!
!  Zi              = coordination number of pivot atom i
!  lgrad1          = if .true. calculate the first derivative
!  lgrad2          = if .true. calculate the second derivative
!
!  On exit :
!
!  tauZ            = the value of tau
!  dtauZdZ         = the first derivatives of tauZ w.r.t. Zi
!  d2tauZdZ2       = the second derivatives of tauZ w.r.t. Zi
!
!  10/10 Created 
!   1/18 Second derivatives added
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
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  real(dp),    intent(in)             :: Z
  real(dp),    intent(out)            :: tauZ
  real(dp),    intent(out)            :: dtauZdZ
  real(dp),    intent(out)            :: d2tauZdZ2
  logical,     intent(in)             :: lgrad1
  logical,     intent(in)             :: lgrad2
!
!  Local variables
!
  real(dp)                            :: dtanhxdx
  real(dp)                            :: d2tanhxdx2
  real(dp)                            :: exp2x
  real(dp)                            :: tanhx
  real(dp)                            :: twelth
  real(dp)                            :: x
#ifdef TRACE
  call trace_in('edip_tau')
#endif
!
!  Find range
!
  if (Z.le.2.0_dp) then
    tauZ = 1.0_dp
    if (lgrad1) then
      dtauZdZ = 0.0_dp
      if (lgrad2) then
        d2tauZdZ2 = 0.0_dp
      endif
    endif
  elseif (Z.ge.6.0_dp) then
    tauZ = 0.0_dp
    if (lgrad1) then
      dtauZdZ = 0.0_dp
      if (lgrad2) then
        d2tauZdZ2 = 0.0_dp
      endif
    endif
  else
    twelth = 1.0_dp/12.0_dp
    x = 6.0_dp*(Z - 2.5_dp)
    exp2x = exp(2.0_dp*x)
    tanhx = (exp2x - 1.0_dp)/(exp2x + 1.0_dp)
    tauZ = 1.0_dp - twelth*Z*(1.0_dp + tanhx)
    if (lgrad1) then
      dtanhxdx = 6.0_dp*(1.0_dp - tanhx*tanhx)
      dtauZdZ = - twelth*((1.0_dp + tanhx) + Z*dtanhxdx)
      if (lgrad2) then
        d2tanhxdx2 = -12.0_dp*tanhx*dtanhxdx
        d2tauZdZ2 = - twelth*(2.0_dp*dtanhxdx + Z*d2tanhxdx2)
      endif
    endif
  endif
#ifdef TRACE
  call trace_out('edip_tau')
#endif
!
  return
  end
