  subroutine GthetaBO(xij,yij,zij,xik,yik,zik,nBOtype,c,h,Gijk,dGijkdr, &
    d2Gijkdr2,d3Gijkdr3,lgrad1,lgrad2,lgrad3)
!
!  Subroutine to calculate the G(theta) function for Bond Order potentials
!
!  On entry : 
!
!  xij             = x component of i->j vector
!  yij             = y component of i->j vector
!  zij             = z component of i->j vector
!  xik             = x component of i->k vector
!  yik             = y component of i->k vector
!  zik             = z component of i->k vector
!  nBOtype         = integer pointer to type of bond order expression
!  c               = c coefficients / d coefficient
!  h               = h coefficient
!  lgrad1          = if .true. calculate the first derivative
!  lgrad2          = if .true. calculate the second derivative
!  lgrad3          = if .true. calculate the third derivative
!
!  On exit :
!
!  Gijk            = the value of the function G(theta)
!  dGijkdr(3)      = the first derivatives of G(theta) w.r.t. the
!                    three different interatomic vectors
!  d2Gijkdr2(6)    = the second derivatives of G(theta) w.r.t. the
!                    six combinations of 2 different vectors
!  d3Gijkdr3(10)   = the third derivatives of G(theta) w.r.t. the
!                    ten combinations of 3 different vectors
!
!  11/03 Created from gtheta
!  11/04 Pi accessed from module
!  12/07 Unused variables removed
!   9/15 d removed from argument list
!   9/15 nBOtype added to argument list
!   3/16 Modified to handle case where c2 = 0
!   1/18 Modified for MMP potential
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
  use bondorderdata
  use iochannels
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)             :: nBOtype
  real(dp),    intent(in)             :: xij
  real(dp),    intent(in)             :: yij
  real(dp),    intent(in)             :: zij
  real(dp),    intent(in)             :: xik
  real(dp),    intent(in)             :: yik
  real(dp),    intent(in)             :: zik
  real(dp),    intent(in)             :: c(5)
  real(dp),    intent(in)             :: h
  real(dp),    intent(out)            :: Gijk
  real(dp),    intent(out)            :: dGijkdr(3)
  real(dp),    intent(out)            :: d2Gijkdr2(6)
  real(dp),    intent(out)            :: d3Gijkdr3(10)
  logical,     intent(in)             :: lgrad1
  logical,     intent(in)             :: lgrad2
  logical,     intent(in)             :: lgrad3
!
!  Local variables
!
  real(dp)                            :: bGo
  real(dp)                            :: c0
  real(dp)                            :: c1
  real(dp)                            :: c2
  real(dp)                            :: c3
  real(dp)                            :: c4
  real(dp)                            :: c5
  real(dp)                            :: costheta
  real(dp)                            :: cos2theta
  real(dp)                            :: cos3theta
  real(dp)                            :: cos4theta
  real(dp)                            :: cos1d(3)
  real(dp)                            :: cos2d(6)
  real(dp)                            :: cos3d(10)
  real(dp)                            :: d2
  real(dp)                            :: dcos2theta
  real(dp)                            :: dcos3theta
  real(dp)                            :: dcos4theta
  real(dp)                            :: d2cos2theta
  real(dp)                            :: d2cos3theta
  real(dp)                            :: d2cos4theta
  real(dp)                            :: d3cos3theta
  real(dp)                            :: d3cos4theta
  real(dp)                            :: dGijkdcostheta
  real(dp)                            :: d2Gijkdcostheta2
  real(dp)                            :: d3Gijkdcostheta3
  real(dp)                            :: Ga
  real(dp)                            :: dGadcostheta
  real(dp)                            :: d2Gadcostheta2
  real(dp)                            :: d3Gadcostheta3
  real(dp)                            :: Go
  real(dp)                            :: dGodcostheta
  real(dp)                            :: d2Godcostheta2
  real(dp)                            :: d3Godcostheta3
  real(dp)                            :: gsub
  real(dp)                            :: eGa
  real(dp)                            :: rij
  real(dp)                            :: rik
  real(dp)                            :: rrij
  real(dp)                            :: rrik
  real(dp)                            :: rrij2
  real(dp)                            :: rrik2
  real(dp)                            :: rrij4
  real(dp)                            :: rrik4
  real(dp)                            :: r2ij
  real(dp)                            :: r2ik
  real(dp)                            :: r2jk
  real(dp)                            :: rtrm
  real(dp)                            :: trm
  real(dp)                            :: xjk
  real(dp)                            :: yjk
  real(dp)                            :: zjk
#ifdef TRACE
  call trace_in('gthetabo')
#endif
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
  rij = sqrt(r2ij)
  rik = sqrt(r2ik)
!
!  Calculate reciprocal distances
!
  rrij = 1.0_dp/rij
  rrik = 1.0_dp/rik
!
!  Calculate cos(theta) using cosine rule
!
  costheta = 0.5_dp*(r2ij + r2ik - r2jk)*rrij*rrik
!
!  Set remaining terms for convenience
!
  c2 = c(1)**2
  d2 = c(2)**2
!***********************
!  Calculate G(theta)  *
!***********************
  if (nBOtype.eq.2) then
!*********************
!  Standard Tersoff  *
!*********************
!
!  Calculate G for all other ranges and pivots
!
    if (c2.gt.0.0_dp) then
      trm = d2 + (h - costheta)**2
      rtrm = 1.0_dp/trm
      Gijk = 1.0_dp + c2/d2 - c2*rtrm
      if (lgrad1) then
        dGijkdcostheta = - 2.0_dp*c2*rtrm*rtrm*(h - costheta)
        if (lgrad2) then
          d2Gijkdcostheta2 = - 8.0_dp*c2*rtrm*rtrm*rtrm*(h - costheta)*(h - costheta) + &
            2.0_dp*c2*rtrm*rtrm
          if (lgrad3) then
            d3Gijkdcostheta3 = - 24.0_dp*c2*rtrm*rtrm*rtrm*rtrm*(h - costheta)**3 + &
              16.0_dp*c2*rtrm*rtrm*rtrm*(h - costheta) + 4.0_dp*c2*rtrm*rtrm*rtrm*(h - costheta)
          endif
        endif
      endif
    else
      Gijk = 1.0_dp
      if (lgrad1) then
        dGijkdcostheta = 0.0_dp
        if (lgrad2) then
          d2Gijkdcostheta2 = 0.0_dp
          if (lgrad3) then
            d3Gijkdcostheta3 = 0.0_dp
          endif
        endif
      endif
    endif
  elseif (nBOtype.eq.3) then
!********************
!  Kumagai variant  *
!********************
!
!  Calculate G for all other ranges and pivots
!
    c2 = c(2)
    c3 = c(3)
    c4 = c(4)
    c5 = c(5)
!
    bGo = 1.0_dp/(c3 + (h - costheta)**2)
    Go = bGo*c2*(h - costheta)**2
    eGa = exp(-c5*(h - costheta)**2)
    Ga = 1.0_dp + c4*eGa
!
    Gijk = c(1) + Go*Ga
    if (lgrad1) then
      dGodcostheta = 2.0_dp*bGo*c2*(- (h - costheta) + bGo*(h - costheta)**3)
      dGadcostheta = 2.0_dp*c4*eGa*c5*(h - costheta)
!
      dGijkdcostheta = dGodcostheta*Ga + Go*dGadcostheta
      if (lgrad2) then
        d2Godcostheta2 = 2.0_dp*bGo*c2*(1.0_dp - 5.0_dp*bGo*(h - costheta)**2 + &
                         6.0_dp*bGo*bGo*(h - costheta)**4)
        d2Gadcostheta2 = 4.0_dp*c4*eGa*c5*c5*(h - costheta)**2 - 2.0_dp*c4*eGa*c5
!
        d2Gijkdcostheta2 = d2Godcostheta2*Ga + 2.0_dp*dGodcostheta*dGadcostheta + Go*d2Gadcostheta2
        if (lgrad3) then
          d3Godcostheta3 = 2.0_dp*bGo*bGo*c2*(12.0_dp*(h - costheta) - 44.0_dp*bGo*(h - costheta)**3 + &
                           36.0_dp*bGo*bGo*(h - costheta)**5)
          d3Gadcostheta3 = 4.0_dp*c4*eGa*c5*c5*(h - costheta)*(2.0_dp*c5*(h - costheta)**2 - 3.0_dp)
!
          d3Gijkdcostheta3 = d3Godcostheta3*Ga + 3.0_dp*d2Godcostheta2*dGadcostheta + &
            3.0_dp*dGodcostheta*d2Gadcostheta2 + Go*d3Gadcostheta3
        endif
      endif
    endif
  elseif (nBOtype.eq.4) then
!****************
!  MMP variant  *
!****************
!
!  Use multiple angle formulae
!
    cos2theta = 2.0_dp*costheta**2 - 1.0_dp
    cos3theta = costheta*(4.0_dp*costheta**2 - 3.0_dp)
    cos4theta = 8.0_dp*(costheta**2)*(costheta**2 - 1.0_dp) + 1.0_dp
    c0 = c(1)
    c1 = c(2)
    c2 = c(3)
    c3 = c(4)
    c4 = c(5)
    gsub = c0 + c1*costheta + c2*cos2theta + c3*cos3theta + c4*cos4theta
    Gijk = gsub**2
    if (lgrad1) then
      dcos2theta = 4.0_dp*costheta
      dcos3theta = 12.0_dp*costheta**2 - 3.0_dp
      dcos4theta = 8.0_dp*costheta*(4.0_dp*costheta**2 - 2.0_dp)
      dGijkdcostheta = 2.0_dp*gsub*(c1 + c2*dcos2theta + c3*dcos3theta + c4*dcos4theta)
      if (lgrad2) then
        d2cos2theta = 4.0_dp
        d2cos3theta = 24.0_dp*costheta
        d2cos4theta = 8.0_dp*(12.0_dp*costheta**2 - 2.0_dp)
        d2Gijkdcostheta2 = 2.0_dp*gsub*(c2*d2cos2theta + c3*d2cos3theta + c4*d2cos4theta) + &
                           2.0_dp*(c1 + c2*dcos2theta + c3*dcos3theta + c4*dcos4theta)**2
        if (lgrad3) then
          d3cos3theta = 24.0_dp
          d3cos4theta = 192.0_dp*costheta
          d3Gijkdcostheta3 = 2.0_dp*gsub*(c3*d3cos3theta + c4*d3cos4theta) + &
                             4.0_dp*(c1 + c2*dcos2theta + c3*dcos3theta + c4*dcos4theta)* &
                                    (c2*d2cos2theta + c3*d2cos3theta + c4*d2cos4theta)
        endif
      endif
    endif
  else
    Gijk = 0.0_dp
    if (lgrad1) then
      dGijkdcostheta = 0.0_dp
      if (lgrad2) then
        d2Gijkdcostheta2 = 0.0_dp
        if (lgrad3) then
          d3Gijkdcostheta3 = 0.0_dp
        endif
      endif
    endif
  endif
!**************************
!  Calculate derivatives  *
!**************************
  if (lgrad1) then
    rrij2 = rrij*rrij
    rrik2 = rrik*rrik
    rrij4 = rrij2*rrij2
    rrik4 = rrik2*rrik2
!
!  First derivatives of cos(theta)
!
!  1 = ij
!  2 = ik
!  3 = jk
!
    cos1d(1) = rrij*rrik - costheta*rrij2
    cos1d(2) = rrij*rrik - costheta*rrik2
    cos1d(3) = - rrij*rrik
!
!  First derivatives of Gijk
!
    dGijkdr(1:3) = dGijkdcostheta*cos1d(1:3)
!
    if (lgrad2) then
!
!  Second derivatives of cos(theta)
!
!  1 = ij/ij
!  2 = ik/ij
!  3 = jk/ij
!  4 = ik/ik
!  5 = jk/ik
!  6 = jk/jk
!
      cos2d(1) = - 2.0_dp*rrij2*rrij*rrik + 3.0_dp*costheta*rrij4
      cos2d(2) = costheta*rrij2*rrik2 - rrij*rrik*(rrij2 + rrik2)
      cos2d(3) = rrij2*rrij*rrik
      cos2d(4) = - 2.0_dp*rrik2*rrik*rrij + 3.0_dp*costheta*rrik4
      cos2d(5) = rrik2*rrij*rrik
      cos2d(6) = 0.0_dp
!
!  Second derivatives of Gijk
!
      d2Gijkdr2(1) = dGijkdcostheta*cos2d(1) + d2Gijkdcostheta2*cos1d(1)*cos1d(1)
      d2Gijkdr2(2) = dGijkdcostheta*cos2d(2) + d2Gijkdcostheta2*cos1d(2)*cos1d(1)
      d2Gijkdr2(3) = dGijkdcostheta*cos2d(3) + d2Gijkdcostheta2*cos1d(3)*cos1d(1)
      d2Gijkdr2(4) = dGijkdcostheta*cos2d(4) + d2Gijkdcostheta2*cos1d(2)*cos1d(2)
      d2Gijkdr2(5) = dGijkdcostheta*cos2d(5) + d2Gijkdcostheta2*cos1d(3)*cos1d(2)
      d2Gijkdr2(6) = dGijkdcostheta*cos2d(6) + d2Gijkdcostheta2*cos1d(3)*cos1d(3)
      if (lgrad3) then
!
!  Third derivatives of cos(theta)
!
!  1 = 111
!  2 = 211
!  3 = 311
!  4 = 221
!  5 = 321
!  6 = 331
!  7 = 222
!  8 = 322
!  9 = 332
! 10 = 333
!
        cos3d(1) = rrij4*(9.0_dp*rrij*rrik - 15.0_dp*costheta*rrij2)
        cos3d(2) = rrij2*(2.0_dp*rrij*rrik2*rrik + 3.0_dp*rrij2*rrij*rrik -  &
                   3.0_dp*costheta*rrij2*rrik2)
        cos3d(3) = - 3.0_dp*rrij4*rrij*rrik
        cos3d(4) = rrik2*(2.0_dp*rrik*rrij2*rrij + 3.0_dp*rrik2*rrik*rrij -  &
                   3.0_dp*costheta*rrik2*rrij2)
        cos3d(5) = - rrij2*rrij*rrik2*rrik
        cos3d(6) = 0.0_dp
        cos3d(7) = rrik4*(9.0_dp*rrik*rrij - 15.0_dp*costheta*rrik2)
        cos3d(8) = - 3.0_dp*rrik4*rrik*rrij
        cos3d(9) = 0.0_dp
        cos3d(10) = 0.0_dp
!
!  Third derivatives of Gijk
!
        d3Gijkdr3(1) = d3Gijkdcostheta3*cos1d(1)*cos1d(1)*cos1d(1) &
          + d2Gijkdcostheta2*(3.0_dp*cos2d(1)*cos1d(1)) &
          + dGijkdcostheta*cos3d(1)
        d3Gijkdr3(2) = d3Gijkdcostheta3*cos1d(2)*cos1d(1)*cos1d(1) &
          + d2Gijkdcostheta2*(2.0_dp*cos2d(2)*cos1d(1) &
          + cos2d(1)*cos1d(2)) + dGijkdcostheta*cos3d(2)
        d3Gijkdr3(3) = d3Gijkdcostheta3*cos1d(3)*cos1d(1)*cos1d(1) &
          + d2Gijkdcostheta2*(2.0_dp*cos2d(3)*cos1d(1) &
          + cos2d(1)*cos1d(3)) + dGijkdcostheta*cos3d(3)
        d3Gijkdr3(4) = d3Gijkdcostheta3*cos1d(2)*cos1d(2)*cos1d(1) &
          + d2Gijkdcostheta2*(2.0_dp*cos2d(2)*cos1d(2) &
          + cos2d(4)*cos1d(1)) + dGijkdcostheta*cos3d(4)
        d3Gijkdr3(5) = d3Gijkdcostheta3*cos1d(3)*cos1d(2)*cos1d(1) &
          + d2Gijkdcostheta2*(cos2d(5)*cos1d(1) + cos2d(3)*cos1d(2) &
          + cos2d(2)*cos1d(3)) + dGijkdcostheta*cos3d(5)
        d3Gijkdr3(6) = d3Gijkdcostheta3*cos1d(3)*cos1d(3)*cos1d(1) &
          + d2Gijkdcostheta2*(2.0_dp*cos2d(3)*cos1d(3) &
          + cos2d(6)*cos1d(1)) + dGijkdcostheta*cos3d(6)
        d3Gijkdr3(7) = d3Gijkdcostheta3*cos1d(2)*cos1d(2)*cos1d(2) &
          + d2Gijkdcostheta2*(3.0_dp*cos2d(4)*cos1d(2)) &
          + dGijkdcostheta*cos3d(7)
        d3Gijkdr3(8) = d3Gijkdcostheta3*cos1d(3)*cos1d(2)*cos1d(2) &
          + d2Gijkdcostheta2*(2.0_dp*cos2d(5)*cos1d(2) &
          + cos2d(4)*cos1d(3)) + dGijkdcostheta*cos3d(8)
        d3Gijkdr3(9) = d3Gijkdcostheta3*cos1d(3)*cos1d(3)*cos1d(2) &
          + d2Gijkdcostheta2*(2.0_dp*cos2d(5)*cos1d(3) &
          + cos2d(6)*cos1d(2)) + dGijkdcostheta*cos3d(9)
        d3Gijkdr3(10) = d3Gijkdcostheta3*cos1d(3)*cos1d(3)*cos1d(3) &
          + d2Gijkdcostheta2*(3.0_dp*cos2d(6)*cos1d(3)) &
          + dGijkdcostheta*cos3d(10)
      endif
    endif
  endif
#ifdef TRACE
  call trace_out('gthetabo')
#endif
!
  return
  end
