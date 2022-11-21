  subroutine EDIP_twobody(nzij,rij,Zi,e2i,de2idrij,de2idZi,d2e2idrij2,d2e2idrijdZi,d2e2idZi2,lgrad1,lgrad2)
!
!  Subroutine to calculate the twobody energy for EDIP. 
!
!  On entry : 
!
!  nzij            = index for pair of EDIP species
!  rij             = distance between i and j
!  Zi              = coordination number for atom i
!  lgrad1          = if .true. compute the first derivative 
!  lgrad2          = if .true. compute the first derivative 
!
!  On exit :
!
!  e2i             = two-body energy for i
!  de2idrij        = first derivative of e2i w.r.t. rij x 1/rij if lgrad1
!  de2idZi         = first derivative of e2i w.r.t. Zi if lgrad1
!  d2e2idrij2      = second derivative of e2i w.r.t. rij x 1/rij if lgrad2
!  d2e2idrijdZi    = second derivative of e2i w.r.t. rij x 1/rij and Zi (mixed) if lgrad2
!  d2e2idZi2       = second derivative of e2i w.r.t. Zi if lgrad2
!
!   9/10 Created
!  10/10 Modified so that only i contribution (rather than i & j) returned
!  10/10 Cutoff applied to 2-body interaction based on a + aprime*Z
!   1/14 EDIP2p introduced a generalised exponent
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
  use EDIPdata
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)             :: nzij
  real(dp),    intent(in)             :: rij
  real(dp),    intent(in)             :: Zi
  real(dp),    intent(out)            :: e2i
  real(dp),    intent(out)            :: de2idrij
  real(dp),    intent(out)            :: de2idZi
  real(dp),    intent(out)            :: d2e2idrij2
  real(dp),    intent(out)            :: d2e2idrijdZi
  real(dp),    intent(out)            :: d2e2idZi2
  logical,     intent(in)             :: lgrad1
  logical,     intent(in)             :: lgrad2
!
!  Local variables
!
  real(dp)                            :: btrm
  real(dp)                            :: cutoff
  real(dp)                            :: cutoff2
  real(dp)                            :: dexp1idrij
  real(dp)                            :: dexp1idZi
  real(dp)                            :: d2exp1idrij2
  real(dp)                            :: d2exp1idrijdZi
  real(dp)                            :: d2exp1idZi2
  real(dp)                            :: exp1i
  real(dp)                            :: exp2i
  real(dp)                            :: rrij
#ifdef TRACE
  call trace_in('edip_twobody')
#endif
!
!  Initialise value and derivatives 
!
  e2i = 0.0_dp
  if (lgrad1) then
    de2idrij = 0.0_dp
    de2idZi  = 0.0_dp
    if (lgrad2) then
      d2e2idrij2   = 0.0_dp
      d2e2idrijdZi = 0.0_dp
      d2e2idZi2    = 0.0_dp
    endif
  endif
!
!  If this is a valid pair then calculate twobody energy
!
  if (lEDIPpairOK(nzij)) then
    cutoff  = EDIP2a(nzij) + EDIP2aprime(nzij)*Zi
    cutoff2 = cutoff + EDIP2sigma(nzij)*EDIPaccuracy2drmax
    if (rij.lt.cutoff2) then
      rrij = 1.0_dp/rij
      call EDIP_expfn(rij,cutoff,EDIP2aprime(nzij),EDIP2sigma(nzij),exp1i,dexp1idrij,dexp1idZi, &
                      d2exp1idrij2,d2exp1idrijdZi,d2exp1idZi2,lgrad1,lgrad2)
      exp1i = exp1i*EDIP2epsilon(nzij)
      exp2i = exp(-EDIP2beta(nzij)*Zi*Zi)
      btrm  = (EDIP2B(nzij)*rrij)**EDIP2p(nzij)
      e2i = (btrm - exp2i)*exp1i
      if (lgrad1) then
        dexp1idrij = dexp1idrij*EDIP2epsilon(nzij)
        dexp1idZi  = dexp1idZi*EDIP2epsilon(nzij)
!
        de2idrij = ((btrm - exp2i)*dexp1idrij - exp1i*(EDIP2p(nzij)*btrm*rrij))*rrij
        de2idZi  = (btrm - exp2i)*dexp1idZi + 2.0_dp*exp1i*EDIP2beta(nzij)*Zi*exp2i
!
        if (lgrad2) then
          d2exp1idrij2   = d2exp1idrij2*EDIP2epsilon(nzij)
          d2exp1idrijdZi = d2exp1idrijdZi*EDIP2epsilon(nzij)
          d2exp1idZi2    = d2exp1idZi2*EDIP2epsilon(nzij)
!
          d2e2idrij2   = (btrm - exp2i)*d2exp1idrij2 - 2.0_dp*dexp1idrij*(EDIP2p(nzij)*btrm*rrij) &
                         + exp1i*(EDIP2p(nzij)*(EDIP2p(nzij) + 1.0_dp)*btrm*rrij*rrij)
          d2e2idrijdZi = (btrm - exp2i)*d2exp1idrijdZi - dexp1idZi*(EDIP2p(nzij)*btrm*rrij) &
                         + dexp1idrij*2.0_dp*EDIP2beta(nzij)*Zi*exp2i
          d2e2idZi2    = (btrm - exp2i)*d2exp1idZi2 + 4.0_dp*dexp1idZi*EDIP2beta(nzij)*Zi*exp2i &
                         + exp1i*2.0_dp*EDIP2beta(nzij)*exp2i*(1.0_dp - 2.0_dp*EDIP2beta(nzij)*Zi*Zi)
        endif
      endif
    endif
  endif
#ifdef TRACE
  call trace_out('edip_twobody')
#endif
!
  return
  end
!
  subroutine EDIP_expfn(rij,cutoff,aprime,sigma,exp2i,dexp2idrij,dexp2idZi, &
                        d2exp2idrij2,d2exp2idrijdZi,d2exp2idZi2,lgrad1,lgrad2)
!
!  Subroutine to calculate the exponential cutoff function for EDIP. 
!
!  On entry : 
!
!  rij             = distance between i and j
!  cutoff          = basic cut-off between pair within accuracy tolerance
!  aprime          = multiplier for Zi in cutoff
!  sigma           = numerator in exponential
!  lgrad1          = if .true. compute the first derivatives
!  lgrad2          = if .true. compute the second derivatives
!
!  On exit :
!
!  exp2i           = exponential function for i-j pair
!  dexp2idrij      = first derivative of exp2i w.r.t. rij if lgrad1
!  dexp2idZi       = first derivative of exp2i w.r.t. Zi if lgrad1
!  d2exp2idrij2    = second derivative of exp2i w.r.t. rij if lgrad2
!  d2exp2idrijdZi  = second derivative of exp2i w.r.t. rij and Zi (mixed) if lgrad2
!  d2exp2idZi2     = second derivative of exp2i w.r.t. Zi if lgrad2
!
!  10/10 Created from EDIP_twobody
!  11/12 Zi removed from argument list as it is not used
!  10/14 Tapered derivatives corrected
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
  use g_constants,    only : pi
  use EDIPdata
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  real(dp),    intent(in)             :: rij
  real(dp),    intent(in)             :: cutoff
  real(dp),    intent(in)             :: aprime
  real(dp),    intent(in)             :: sigma
  real(dp),    intent(out)            :: exp2i
  real(dp),    intent(out)            :: dexp2idrij
  real(dp),    intent(out)            :: dexp2idZi
  real(dp),    intent(out)            :: d2exp2idrij2
  real(dp),    intent(out)            :: d2exp2idrijdZi
  real(dp),    intent(out)            :: d2exp2idZi2
  logical,     intent(in)             :: lgrad1
  logical,     intent(in)             :: lgrad2
!
!  Local variables
!
  real(dp)                            :: cosx
  real(dp)                            :: sinx
  real(dp)                            :: cutoff1
  real(dp)                            :: cutoff2
  real(dp)                            :: exp1i
  real(dp)                            :: rd
  real(dp)                            :: rrzi
  real(dp)                            :: t
  real(dp)                            :: dtdr
  real(dp)                            :: d2tdr2
  real(dp)                            :: x
#ifdef TRACE
  call trace_in('edip_expfn')
#endif
!
!  Initialise value and derivatives 
!
  exp2i = 0.0_dp
  if (lgrad1) then
    dexp2idrij = 0.0_dp
    dexp2idZi  = 0.0_dp
    if (lgrad2) then
      d2exp2idrij2   = 0.0_dp
      d2exp2idrijdZi = 0.0_dp
      d2exp2idZi2    = 0.0_dp
    endif
  endif
!
!  If this is a valid pair then calculate twobody energy
!
  cutoff2 = cutoff + sigma*EDIPaccuracy2drmax
  if (rij.lt.cutoff2) then
    cutoff1 = cutoff + sigma*EDIPaccuracy1drmax
    rrzi = 1.0_dp/(rij - cutoff)
    exp1i = exp(sigma*rrzi)
    if (rij.lt.cutoff1) then
!
!  Function only
!
      exp2i = exp1i
      if (lgrad1) then
        dexp2idrij = - exp2i*sigma*rrzi*rrzi
        dexp2idZi  = - dexp2idrij*aprime
        if (lgrad2) then
          d2exp2idrij2   = exp2i*sigma*rrzi*rrzi*rrzi*(2.0_dp + sigma*rrzi)
          d2exp2idrijdZi = - d2exp2idrij2*aprime
          d2exp2idZi2    = d2exp2idrij2*aprime*aprime
        endif
      endif
    else
!
!  Multiply by taper function 
!
      rd = pi/(cutoff2 - cutoff1)
      x = rd*(rij - cutoff1)
      cosx = cos(x)
      t = 0.5_dp*(1.0_dp + cosx)
      exp2i = exp1i*t
      if (lgrad1) then
        sinx = sin(x)
        dtdr = - 0.5_dp*sinx*rd
        dexp2idrij = - exp1i*sigma*rrzi*rrzi*t + exp1i*dtdr
        dexp2idZi  = - dexp2idrij*aprime
        if (lgrad2) then
          d2tdr2 = - 0.5_dp*cosx*rd*rd
          d2exp2idrij2   = exp2i*sigma*rrzi*rrzi*rrzi*(2.0_dp + sigma*rrzi)*t + &
                           2.0_dp*dexp2idrij*dtdr + exp1i*d2tdr2
          d2exp2idrijdZi = - d2exp2idrij2*aprime
          d2exp2idZi2    = d2exp2idrij2*aprime*aprime
        endif
      endif
    endif
  endif
#ifdef TRACE
  call trace_out('edip_expfn')
#endif
!
  return
  end
