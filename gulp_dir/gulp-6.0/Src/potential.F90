  subroutine electrostatic_potential
!
!  Subroutine for evaluating electrostatic potential at a
!  series of general points. Used in electrostatic potential
!  surface fitting.
!
!   7/00 referencing of xsite,ysite,zsite shift back by npotpt0
!   7/00 arrays filled for 0-D case as well
!   8/01 Call to epot0/3 replaced with generic call to epot
!   3/14 Name changed from potential to electrostatic_potential
!        for the benefit of ChemShell
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
  use control
  use current
  use potentialpoints
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Local variables
!
  integer(i4)                               :: i
  integer(i4)                               :: nsite
  integer(i4)                               :: status
  real(dp)                                  :: g_cpu_time
  real(dp)                                  :: efg
  real(dp)                                  :: time1
  real(dp)                                  :: time2
  real(dp)                                  :: vx
  real(dp)                                  :: vy
  real(dp)                                  :: vz
  real(dp), dimension(:), allocatable       :: xsite
  real(dp), dimension(:), allocatable       :: ysite
  real(dp), dimension(:), allocatable       :: zsite
#ifdef TRACE
  call trace_in('potential')
#endif
!
  time1 = g_cpu_time()
!
!  Set pointer to points
!
  npotpt0 = 0
  do i = 1,ncf - 1
    npotpt0 = npotpt0 + npotptcfg(i)
  enddo
  nsite = npotptcfg(ncf)
!
!  Allocate local memory
!
  allocate(xsite(nsite),stat=status)
  if (status/=0) call outofmemory('electrostatic_potential','xsite')
  allocate(ysite(nsite),stat=status)
  if (status/=0) call outofmemory('electrostatic_potential','ysite')
  allocate(zsite(nsite),stat=status)
  if (status/=0) call outofmemory('electrostatic_potential','zsite')
!
!  Generate cartesian coordinates
!
  if (ndim.eq.3) then
!**************
!  3D system  *
!**************
    do i = npotpt0 + 1,npotpt0 + nsite
      xsite(i-npotpt0) = xpotpt(i)*r1x + ypotpt(i)*r2x + zpotpt(i)*r3x
      ysite(i-npotpt0) = xpotpt(i)*r1y + ypotpt(i)*r2y + zpotpt(i)*r3y
      zsite(i-npotpt0) = xpotpt(i)*r1z + ypotpt(i)*r2z + zpotpt(i)*r3z
    enddo
  elseif (ndim.eq.2) then
!**************
!  2D system  *
!**************
    do i = npotpt0 + 1,npotpt0 + nsite
      xsite(i-npotpt0) = xpotpt(i)*r1x + ypotpt(i)*r2x
      ysite(i-npotpt0) = xpotpt(i)*r1y + ypotpt(i)*r2y
      zsite(i-npotpt0) = zpotpt(i)
    enddo
  elseif (ndim.eq.1) then
!**************
!  1D system  *
!**************
    do i = npotpt0 + 1,npotpt0 + nsite
      xsite(i-npotpt0) = xpotpt(i)*r1x
      ysite(i-npotpt0) = ypotpt(i)
      zsite(i-npotpt0) = zpotpt(i)
    enddo
  elseif (ndim.eq.0) then
!**************
!  0D system  *
!**************
    do i = npotpt0 + 1,npotpt0 + nsite
      xsite(i-npotpt0) = xpotpt(i)
      ysite(i-npotpt0) = ypotpt(i)
      zsite(i-npotpt0) = zpotpt(i)
    enddo
  endif
!************************
!  Calculate potential  *
!************************
  call epot(.true.,nsite,vpotpt(npotpt0+1),xsite,ysite,zsite,.false.,vx,vy,vz,.false.,efg,.false.)
!
!  Free local memory
!
  deallocate(zsite,stat=status)
  if (status/=0) call deallocate_error('electrostatic_potential','zsite')
  deallocate(ysite,stat=status)
  if (status/=0) call deallocate_error('electrostatic_potential','ysite')
  deallocate(xsite,stat=status)
  if (status/=0) call deallocate_error('electrostatic_potential','xsite')
!
!  Timing
!
  time2 = g_cpu_time()
  tion = tion + time2 - time1
#ifdef TRACE
  call trace_out('potential')
#endif
!
  return
  end
