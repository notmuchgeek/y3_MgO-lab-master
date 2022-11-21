  subroutine mdvvcorrect
!
!  Performs the corrector step of the Velocity Verlet algorithm
!
!   7/97 Created from mdcorrect
!  10/02 Constraint force modifications added
!   7/05 lfirststep argument added to mdke call
!   7/08 Constraint code moved to subroutine
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
  use current
  use derivatives
  use moldyn
  use optimisation
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  use velocities
  implicit none
!
!  Local variables
!
  integer(i4) :: i
  real(dp)    :: rrmi
  real(dp)    :: sfac2
  real(dp)    :: velrsq
  real(dp)    :: xacc
  real(dp)    :: yacc
  real(dp)    :: zacc
#ifdef TRACE
  call trace_in('mdvvcorrect')
#endif
!
  if (nensemble(ncf).eq.2) then
!********
!  NVT  *
!********
    sfac0 = sfac
!
!  Estimate change in thermostat variable
!
    call mdke(velrsq,.false.)
    svel = qtemp(ncf)*((ekin/smdfctt)-1.0_dp)
    sfac = sfac0 + svel
    sfac2 = 0.5_dp*(sfac0 + sfac)
    do i = 1,numat
      if (lopf(i).and..not.lfix(i)) then
        rrmi = -rmass(i)*stpsqh
        xacc = rrmi*xdrv(i)
        yacc = rrmi*ydrv(i)
        zacc = rrmi*zdrv(i)
        velx(i) = velx(i) + x2(i) + xacc - sfac2*velx(i)
        vely(i) = vely(i) + y2(i) + yacc - sfac2*vely(i)
        velz(i) = velz(i) + z2(i) + zacc - sfac2*velz(i)
        x2(i) = xacc
        y2(i) = yacc
        z2(i) = zacc
      endif
    enddo
  elseif (nensemble(ncf).eq.1) then
!********
!  NVE  *
!********
    do i = 1,numat
      if (lopf(i).and..not.lfix(i)) then
        rrmi = -rmass(i)*stpsqh
        xacc = rrmi*xdrv(i)
        yacc = rrmi*ydrv(i)
        zacc = rrmi*zdrv(i)
        velx(i) = velx(i) + x2(i) + xacc
        vely(i) = vely(i) + y2(i) + yacc
        velz(i) = velz(i) + z2(i) + zacc
        x2(i) = xacc
        y2(i) = yacc
        z2(i) = zacc
      endif
    enddo
  endif
!
!  Distance constraint
!
  if (lmdconstrain(ncf)) call mddistconpost
#ifdef TRACE
  call trace_out('mdvvcorrect')
#endif
!
  return
  end
