  subroutine optshe(istep,etot,iter,spring,bspring)
!
!  Subroutine for optimizing the positions of shells during an MD run.
!
!  Assumes that all shell positions are variable.
!
!   2/97 Modified for breathing shells (JDG)
!   1/05 Logical to compute interatomic vector table added in call
!        to mdfunct. This table only really needs setting once. (JDG)
!   4/05 Extrapolation of previous shell positions added to improve optimisation (JDG)
!   7/05 Use of ratiom replaced by ladiabatic logical
!   5/07 nbsmptr replaced by nbsptr
!  10/08 New default algorithm added in which the core-shell vector is used to initialise
!        the shell positions for the new step
!  10/08 Extrapolation is now made on core-shell vector
!   8/11 Call to mdfunct modified to pass step number
!   8/11 Arguments into routine now include step number 
!  10/13 Modifed to call shmin which uses a conjugate gradient minimiser
!  10/13 Extrapolation modified so that there is a check for excessive displacements
!        If this occurs then the last displacement is used instead.
!   2/18 Trace added
!
!  Initially written by J.R. Hill, January 1996
!
!  Last modified to by Julian Gale, Feb 2018
!
  use control
  use current
  use derivatives
  use iochannels
  use mdlogic,            only : ladiabatic
  use moldyn
  use optimisation
  use parallel
  use shells
  use shellextrapolation
#ifdef TRACE
  use trace,              only : trace_in, trace_out
#endif
  use velocities
  implicit none
!
!  Passed variables
!
  integer(i4),     intent(in)    :: istep           ! Step number of MD run
  integer(i4),     intent(out)   :: iter            ! Number of iterations 
  real(dp),        intent(in)    :: bspring(*)      ! Spring constants for breathing shells
  real(dp),        intent(out)   :: etot            ! Final energy
  real(dp),        intent(in)    :: spring(*)       ! Spring constants for shells
!
!  Local variables
!
  integer(i4)                    :: i
  integer(i4)                    :: ifail
  integer(i4)                    :: icp
  integer(i4)                    :: isp
  integer(i4)                    :: m
  integer(i4),              save :: nextrapol = 0
  integer(i4)                    :: status
  logical                        :: loutput
  real(dp)                       :: disp_lim2
  real(dp)                       :: rcs
  real(dp)                       :: rrmi
  real(dp),    allocatable       :: tmp1(:)
  real(dp),    allocatable       :: tmp2(:)
#ifdef TRACE
  call trace_in('optshe')
#endif
!
  loutput = (lsopt.and.ioproc)
!
  if (lextrapolateshells.and.nextrapol.gt.3) then
!
!  Extrapolate shell positions using rational functions 
!  only start once 3 points are already known though
!
    m = min(nextrapol,maxextrapol)
    allocate(tmp1(m),stat=status)
    if (status/=0) call outofmemory('optshe','tmp1')
    allocate(tmp2(m),stat=status)
    if (status/=0) call outofmemory('optshe','tmp2')
!
!  Set displacement limit based on 80% of cuts squared
!
    disp_lim2 = 0.8_dp*cuts**2
!
!  Extrapolate forward core-shell vector
!
    do i = 1,nshell
      isp = nshptr(i)
      icp = ncsptr(isp)
      call ratfn(m,xshellsave(1,i),xalat(isp),tmp1,tmp2)
      call ratfn(m,yshellsave(1,i),yalat(isp),tmp1,tmp2)
      call ratfn(m,zshellsave(1,i),zalat(isp),tmp1,tmp2)
      rcs = xalat(isp)**2 + yalat(isp)**2 + zalat(isp)**2
      if (rcs.gt.disp_lim2) then
!
!  Reset shell displacement to last point
!
        xalat(isp) = xshellsave(1,i)
        yalat(isp) = yshellsave(1,i)
        zalat(isp) = zshellsave(1,i)
      endif
!
!  Add on core position
!
      xalat(isp) = xalat(isp) + xalat(icp)
      yalat(isp) = yalat(isp) + yalat(icp)
      zalat(isp) = zalat(isp) + zalat(icp)
    enddo
!
    deallocate(tmp2,stat=status)
    if (status/=0) call deallocate_error('optshe','tmp2')
    deallocate(tmp1,stat=status)
    if (status/=0) call deallocate_error('optshe','tmp1')
  elseif (nextrapol.eq.1) then
!
!  Position shell based on previous core-shell vector
!
    do i = 1,nshell
      isp = nshptr(i)
      icp = ncsptr(isp)
      xalat(isp) = xalat(icp) + xshellsave(1,i)
      yalat(isp) = yalat(icp) + yshellsave(1,i)
      zalat(isp) = zalat(icp) + zshellsave(1,i)
    enddo
  endif
!
!  Call conjugate gradient minimiser
!
  call shmin(etot,spring,bspring,sgtol,loutput,iter,ifail)
!
  if (.not.ladiabatic) then
    do i = 1,numat
      if (lopf(i).and..not.lfix(i)) then
        rrmi = rmass(i)
        x2(i) = - rrmi*xdrv(i)*stpsqh
        y2(i) = - rrmi*ydrv(i)*stpsqh
        z2(i) = - rrmi*zdrv(i)*stpsqh
      endif
    enddo
#ifdef TRACE
    call trace_out('optshe')
#endif
    return
  endif
  if (lextrapolateshells) then
!
!  Store shell positions for future extrapolation
!
    nextrapol = nextrapol + 1
    if (nextrapol.gt.maxextrapol) then
      do m = 2,maxextrapol
        do i = 1,nshell
          xshellsave(m-1,i) = xshellsave(m,i)
          yshellsave(m-1,i) = yshellsave(m,i)
          zshellsave(m-1,i) = zshellsave(m,i)
        enddo
      enddo
    endif
    m = min(nextrapol,maxextrapol)
    do i = 1,nshell
      isp = nshptr(i)
      icp = ncsptr(isp)
      xshellsave(m,i) = xalat(isp) - xalat(icp)
      yshellsave(m,i) = yalat(isp) - yalat(icp)
      zshellsave(m,i) = zalat(isp) - zalat(icp)
    enddo
  else
!
!  Store position of shell relative to its core
!
    nextrapol = 1
    do i = 1,nshell
      isp = nshptr(i)
      icp = ncsptr(isp)
      xshellsave(1,i) = xalat(isp) - xalat(icp)
      yshellsave(1,i) = yalat(isp) - yalat(icp)
      zshellsave(1,i) = zalat(isp) - zalat(icp)
    enddo
  endif
#ifdef TRACE
  call trace_out('optshe')
#endif
!
  return
  end
