  subroutine lambda(nmode,nlower,nupper,nvptr,nactual,rlam,tmp,eig)
!
!  Finds the optimum value of lambda in the RFO method
!
!   5/13 Trap for nonsensical negative rlam added
!   9/15 Modified to ensure that lambda corresponding to lowest function
!        is returned even when convergence fails
!   9/15 Limit on lowest value removed since this can prevent the use
!        of the correct value if the first eigenvalue is very negative
!   2/18 Trace added
!  10/18 Parallel I/O handling added
!
!  nmode  = index of mode with eigenvalue greater than that required
!  rlam   = lambda value
!  nlower = lowest mode to be included
!  nupper = upper mode to be included
!  tmp    = transformed gradient vectors dot products with themselves
!  eig    = eigenvalues of Hessian
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
!  Julian Gale, CIC, Curtin University, October 2018
!
  use control
  use iochannels
  use parallel,      only : ioproc
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4)     :: j
  integer(i4)     :: jj
  integer(i4)     :: maxcyc
  integer(i4)     :: nactual
  integer(i4)     :: ncyc
  integer(i4)     :: nlower
  integer(i4)     :: nmode
  integer(i4)     :: nupper
  integer(i4)     :: nvptr(*)
  logical         :: lconverged
  real(dp)        :: d1
  real(dp)        :: d11
  real(dp)        :: d2
  real(dp)        :: d22
  real(dp)        :: diff
  real(dp)        :: eig(*)
  real(dp)        :: fn
  real(dp)        :: fnlowest
  real(dp)        :: fnlast
  real(dp)        :: rlam
  real(dp)        :: rlamlowest
  real(dp)        :: rtrm
  real(dp)        :: step
  real(dp)        :: stepmx
  real(dp)        :: tmp(*)
  real(dp)        :: trm
#ifdef TRACE
  call trace_in('lambda')
#endif
!
!  Find initial lambda value
!
  if (nmode.eq.1) then
    rlam = eig(nvptr(1)) - 1.0_dp
  elseif (nmode.gt.nactual) then
    rlam = eig(nvptr(nactual)) + 1.0_dp
  else
    rlam = eig(nvptr(nmode)) + eig(nvptr(nmode-1))
    rlam = 0.5_dp*rlam - 0.001_dp
  endif
!
!  Iterate to find optimum lambda
!
  maxcyc = 100
  fnlast = 0.0_dp
  fnlowest = 1.0d6
  rlamlowest = 1.0_dp
  fn = 0.0_dp
  lconverged = .false.
  ncyc = 0
  do while (.not.lconverged.and.ncyc.lt.maxcyc)
    ncyc = ncyc + 1
    fnlast = fn
!
!  Find function, first and second derivatives
!
    fn = rlam
    d1 = 1.0_dp
    d2 = 0.0_dp
    do j = nlower,nupper
      jj = nvptr(j)
      trm = tmp(jj)
      diff = rlam - eig(jj)
      if (abs(diff).lt.1.0d-8) diff = sign(1.0d-8,diff)
      rtrm = 1.0_dp/diff
      trm = trm*rtrm
      fn = fn - trm
      trm = trm*rtrm
      d1 = d1 + trm
      trm = trm*rtrm
      d2 = d2 - 2.0_dp*trm
    enddo
    d11 = 2.0_dp*d1*fn
    d22 = 2.0_dp*(d1*d1+fn*d2)
    fn = fn*fn
    step = - d11/abs(d22)
!
!  Derive stepmx
!
    if (step.lt.0.0_dp) then
      if (nmode.gt.nupper) then
        stepmx = 0.5_dp*(rlam-eig(nvptr(nmode-1)))
      else
        stepmx = 100.0_dp
      endif
    else
      if (nmode.gt.nupper) then
        stepmx = 100.0_dp
      else
        stepmx = 0.5_dp*(eig(nvptr(nmode))-rlam)
      endif
    endif
    if (abs(step).gt.stepmx) then
      step = sign(stepmx,step)
    endif
!
!  Apply step
!
    rlam = rlam + step
    if (ldebug.and.ioproc) write(ioout,'('' Current lambda and fn = '',2f12.6)') rlam,fn
    lconverged = (abs(fn).lt.1.0d-4)
!
!  Save best result so far
!
    if (fn.lt.fnlowest) then
      fnlowest = fn
      rlamlowest = rlam
    endif
  enddo
!
  if (ldebug.and.ioproc) then
    if (lconverged) then
      write(ioout,'(''  ** Converged in '',i2,'' cycles: '')') ncyc
    else
      write(ioout,'(''  ** No lower point found after '',i2,'' cycles: '')') ncyc
    endif
    write(ioout,'(''  ** Final difference = '',f12.6)') fn
    write(ioout,'(''  ** Lambda value = '',f12.6)') rlam
  endif
!
!  Check that best value is being returned
!
  if (fn.gt.fnlowest) then
    rlam = rlamlowest
    fn = fnlowest
  endif
#ifdef TRACE
  call trace_out('lambda')
#endif
!
  return
  end
