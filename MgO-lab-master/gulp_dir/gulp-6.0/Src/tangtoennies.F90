  subroutine tangtoennies(kmax,b,r,f,dfdr,d2fdr2,d3fdr3,lgrad1,lgrad2,lgrad3)
!
!  Subroutine to calculate the Tang-Toennies damping function.
!
!  On entry : 
!
!  kmax            = maximum power exponent
!  b               = constant that multiplies the distance
!  r               = distance
!  lgrad1          = if .true. compute the first derivative 
!  lgrad2          = if .true. compute the second derivative 
!  lgrad3          = if .true. compute the third derivative 
!
!  On exit :
!
!  f               = damping factor
!  dfdr            = first derivative of f w.r.t. r if lgrad1
!  d2fdr1          = second derivative of f w.r.t. r if lgrad2
!  d3fdr3          = third derivative of f w.r.t. r if lgrad3
!
!   8/19 Created
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
!  Julian Gale, CIC, Curtin University, August 2019
!
  use datatypes
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)             :: kmax
  real(dp),    intent(in)             :: b
  real(dp),    intent(in)             :: r
  real(dp),    intent(out)            :: f
  real(dp),    intent(out)            :: dfdr
  real(dp),    intent(out)            :: d2fdr2
  real(dp),    intent(out)            :: d3fdr3
  logical,     intent(in)             :: lgrad1
  logical,     intent(in)             :: lgrad2
  logical,     intent(in)             :: lgrad3
!
!  Local variables
!
  integer(i4)                         :: k
  real(dp)                            :: br
  real(dp)                            :: brk
  real(dp)                            :: brkm1
  real(dp)                            :: brkm2
  real(dp)                            :: brkm3
  real(dp)                            :: expbr
  real(dp)                            :: kfactorial
  real(dp)                            :: sum
  real(dp)                            :: dsumdr
  real(dp)                            :: d2sumdr2
  real(dp)                            :: d3sumdr3
#ifdef TRACE
  call trace_in('tangtoennies')
#endif
!
!  Calculate function and derivatives 
!
  br = b*r
  expbr = exp(-br)
  sum = 1.0_dp
  kfactorial = 1.0_dp
  brk = 1.0_dp
  if (lgrad1) then
    brkm1 = 1.0_dp
    dsumdr = 0.0_dp
    if (lgrad2) then
      brkm2 = 1.0_dp
      d2sumdr2 = 0.0_dp
      if (lgrad3) then
        brkm3 = 1.0_dp
        d3sumdr3 = 0.0_dp
      endif
    endif
  endif
  do k = 1,kmax
    kfactorial = kfactorial*dble(k)
    brk = brk*br
    sum = sum + brk/kfactorial
    if (lgrad1) then
      dsumdr = dsumdr + dble(k)*brkm1/kfactorial
      brkm1 = brkm1*br
      if (lgrad2.and.k.gt.1) then
        d2sumdr2 = d2sumdr2 + dble(k*(k-1))*brkm2/kfactorial
        brkm2 = brkm2*br
        if (lgrad3.and.k.gt.2) then
          d3sumdr3 = d3sumdr3 + dble(k*(k-1)*(k-2))*brkm3/kfactorial
          brkm3 = brkm3*br
        endif
      endif
    endif
  enddo
  f = 1.0_dp - expbr*sum
  if (lgrad1) then
    dfdr = expbr*b*(sum - dsumdr)
    if (lgrad2) then
      d2fdr2 = - expbr*b*b*(sum - 2.0_dp*dsumdr + d2sumdr2)
      if (lgrad3) then
        d3fdr3 = expbr*b*b*b*(sum - 3.0_dp*(dsumdr - d2sumdr2) - d3sumdr3)
      endif
    endif
  endif
#ifdef TRACE
  call trace_out('tangtoennies')
#endif
!
  return
  end
