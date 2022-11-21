  subroutine mdtemperature(nstep,targettemperature)
!
!  Sets temperature for current MD step
!
!  12/07 Target temperature now made a return argument
!   3/08 Modified for new MD integrator
!   5/09 Bug in temperature setting for case when nstep = ntemperaturestep+ntemperaturestepstart
!  11/11 Modifications made for variable temperature with stochastic thermostat
!   2/18 Trace added
!   3/19 Multiple temperature ramps added
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
!  Copyright CIC, Curtin University
!
!  Julian Gale, CIC Curtin University, March 2019
!
  use configurations,     only : tempcfg
  use current
  use moldyn
  use m_pr
  use pr_random_module,   only : gauss_rand
#ifdef TRACE
  use trace,              only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)  :: nstep
  real(dp),    intent(out) :: targettemperature
!
!  Local variables
!
  integer(i4)              :: ntr
  logical                  :: lfound
#ifdef TRACE
  call trace_in('mdtemperature')
#endif
!
!  Determine current temperature
!
  targettemperature = tempcfg(ncf)
  if (ntemperatureramp.gt.0) then
    if (nstep.gt.ntemperaturestepstop(ntemperatureramp)) then
!
!  Past last ramp step
!
      ntr = ntemperatureramp
      targettemperature = temperaturestart(ntr) + dble(ntemperaturestep(ntr))*temperaturestep(ntr)
    else
!
!  During ramp steps
!
      lfound = .false.
      do ntr = 1,ntemperatureramp
        if (nstep.gt.ntemperaturestepstart(ntr).and.nstep.le.ntemperaturestepstop(ntr)) then
          lfound = .true.
          targettemperature = temperaturestart(ntr) + dble(nstep-ntemperaturestepstart(ntr))*temperaturestep(ntr)
        endif
      enddo
      if (.not.lfound) then
!
!  Between ramp steps
!
        do ntr = 2,ntemperatureramp
          if (nstep.gt.ntemperaturestepstop(ntr-1).and.nstep.le.ntemperaturestepstart(ntr)) then
            targettemperature = temperaturestart(ntr-1) + dble(ntemperaturestep(ntr-1))*temperaturestep(ntr-1)
          endif
        enddo
      endif
    endif
  endif
!
!  Set dependent variables
!
  smdfctt = smdfct*targettemperature
  velfctt = velfct*targettemperature
  if (targettemperature.ne.0.0_dp) psfctt = psfct/targettemperature
  if (nmdintegrator.eq.4) then
    pr_target_temp = targettemperature
!
!  Temperature related quantities
!
    pr_ekin = 0.5_dp*dble(ndof_atm)*pr_boltz*targettemperature
    pr_temp = 2.0_dp*pr_ekin/pr_boltz/dble(ndof_atm)
!
!  Target kinetic energies
!
    pr_ekintarget_atm  = 0.5_dp*pr_boltz*pr_target_temp
    pr_ekintarget_baro = 0.5_dp*pr_boltz*pr_target_temp
!
!  Barostat
!
    if (lpr_baro_iso) then
      bmass_iso = 2.0_dp*pr_ekintarget_atm*dble(ndof_atm)*taub**2
      if (.not.lpisoinput) then
        p_iso = sqrt(pr_ekintarget_baro*dble(ndof_baro)*bmass_iso)*gauss_rand()
      endif
      pr_ekinbaro = 0.5_dp/bmass_iso*p_iso**2
    elseif (lpr_baro_flx) then
      bmass_flx = 2.0_dp*pr_ekintarget_atm*dble(ndof_atm)*taub**2/3.0_dp
      pr_ekinbaro = 0.0_dp
    endif
  endif
#ifdef TRACE
  call trace_out('mdtemperature')
#endif
!
  return
  end
