  subroutine numhess(nvar,xc,fc,gc,hesinv,maxhess,lhess2D,xvar,gvar)
!
!  Numerical estimation of hessian diagonal elements
!
!   2/17 nmin removed from arguments
!   2/17 maxhess & lhess2D added as arguments
!   3/17 Parallel handling of hesinv added
!   1/18 Trace added
!   7/20 Trapping of geometry problem added on return from funct
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
  use iochannels
  use parallel
#ifdef TRACE
  use trace,       only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)    :: nvar
  integer(i4), intent(in)    :: maxhess
  logical,     intent(in)    :: lhess2D
  real(dp),    intent(inout) :: gc(*)
  real(dp),    intent(out)   :: gvar(*)
  real(dp),    intent(out)   :: hesinv(maxhess,*)
  real(dp),    intent(inout) :: xc(*)
  real(dp),    intent(out)   :: xvar(*)
!
!  Local variables
!
  integer(i4)                :: i
  integer(i4)                :: ii
  integer(i4)                :: iflag
  integer(i4)                :: il
  real(dp)                   :: del
  real(dp)                   :: deltag
  real(dp)                   :: deltax
  real(dp)                   :: fc
  real(dp)                   :: funct2
  real(dp)                   :: ggd
  real(dp)                   :: pmstep
  real(dp)                   :: tdel
#ifdef TRACE
  call trace_in('numhess')
#endif
!
  del = 0.01_dp
  tdel = 0.06_dp
!
!  Generate new position for differencing
!
  do i = 1,nvar
    xvar(i) = xc(i) - sign(del,gc(i))
  enddo
  iflag = 1
!
!  Calculate function and gradients
!
  call funct(iflag,nvar,xvar,funct2,gvar)
!
!  Check return flag for problem with geometry
!
  if (iflag.eq.-2) then
    call outerror('geometry has become unphysical during numerical hessian',0_i4)
    call stopnow('numhess')
  endif
!
  ii = 0
!
!  Calculate diagonal elements
!
  do il = 1,nvaronnode
    i = node2var(il)
    deltag = gc(i) - gvar(i)
    deltax = xc(i) - xvar(i)
    if (abs(deltag).lt.1.d-12) goto 70
    ggd = abs(gc(i))
    if (funct2.lt.fc) ggd = abs(gvar(i))
    if (lhess2D) then
      hesinv(i,il) = deltax/deltag
      if (hesinv(i,il).lt.0.0_dp.and.ggd.lt.1.d-12) goto 70
      if (hesinv(i,il).lt.0.0_dp) hesinv(i,il) = tdel/ggd
    else
      ii = ii + i 
      hesinv(ii,1) = deltax/deltag
      if (hesinv(ii,1).lt.0.0_dp.and.ggd.lt.1.d-12) goto 70
      if (hesinv(ii,1).lt.0.0_dp) hesinv(ii,1) = tdel/ggd
    endif
    goto 80
70  if (lhess2D) then
      hesinv(i,il) = 0.01_dp
    else
      hesinv(ii,1) = 0.01_dp
    endif
80  continue
    if (ggd.lt.1.d-12) ggd = 1.d-12
    pmstep = abs(0.1_dp/ggd)
    if (lhess2D) then
      if (hesinv(i,il).gt.pmstep) hesinv(i,il) = pmstep
    else
      if (hesinv(ii,1).gt.pmstep) hesinv(ii,1) = pmstep
    endif
  enddo
!
!  If energy has been lowered by move, accept new position as current one
!
  if (funct2.lt.fc) then
    fc = funct2
    do i = 1,nvar
      xc(i) = xvar(i)
      gc(i) = gvar(i)
    enddo
  endif
  if (ioproc) then
    write(ioout,'(''  ** Hessian recalculated - numerical/diagonal only **'')')
  endif
#ifdef TRACE
  call trace_out('numhess')
#endif
!
  return
  end
