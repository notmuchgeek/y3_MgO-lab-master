  subroutine vibenergy(mcv,nllkpt,wkpt,freq,maxfd,rnokpt,lprint,fc)
!
!  Calculates the vibrational energetics from the phonon modes and outputs them
!
!   9/13 Created from phonon
!  11/13 Trap for inverting rkt added when T < 1.0d-6
!  12/16 Modified to handle cluster case
!   7/17 lprinloc moved to only wrap output and not calculation
!   8/17 linear now set for cases other than 0-D
!  10/17 fhenergy moved to energies module
!   2/18 Trace added
!   3/19 Multiple temperature ramps added
!   8/19 Correction to handling of zero point energy in equipartition free energy
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
  use configurations
  use g_constants
  use control
  use current
  use element
  use energies,          only : fhenergy
  use general
  use iochannels
  use parallel
  use properties,        only : cv, entropy
#ifdef TRACE
  use trace,             only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4),  intent(in)                          :: mcv                 ! Number of modes
  integer(i4),  intent(in)                          :: nllkpt              ! Number of kpoints for this configuration
  integer(i4),  intent(in)                          :: maxfd               ! Maximum first dimension of the frequency array
  logical,      intent(in)                          :: lprint              ! If true then output results
  real(dp),     intent(in)                          :: fc                  ! Internal energy
  real(dp),     intent(in)                          :: wkpt(nllkpt)        ! Weights for k points
  real(dp),     intent(in)                          :: freq(maxfd,nllkpt)  ! Vibrational frequencies for the k points
  real(dp),     intent(in)                          :: rnokpt              ! Inverse sum of k point weights
!
!  Local variables
!
  integer(i4)                                       :: i
  integer(i4)                                       :: mcvmax
  integer(i4)                                       :: mcvmin
  integer(i4)                                       :: nimag
  integer(i4)                                       :: nk
  integer(i4)                                       :: nt
  integer(i4)                                       :: nt0
  integer(i4)                                       :: ntr
  integer(i4)                                       :: status
  logical                                           :: linear
  logical                                           :: lnozero
  logical                                           :: lprinloc
  real(dp)                                          :: cmfact
  real(dp)                                          :: cv2
  real(dp)                                          :: ent2
  real(dp)                                          :: factor
  real(dp)                                          :: fe_equipartition
  real(dp)                                          :: freqmin
  real(dp)                                          :: kinenergy
  real(dp)                                          :: rkt
  real(dp)                                          :: rmode
  real(dp)                                          :: rtlnz
  real(dp),     dimension(:),     allocatable       :: rtmp2
  real(dp)                                          :: tem
  real(dp)                                          :: trm
  real(dp)                                          :: trm1
  real(dp)                                          :: trmcv
  real(dp)                                          :: trmen
  real(dp)                                          :: trmfe
  real(dp)                                          :: trmfe_eq
  real(dp)                                          :: trmzp
  real(dp)                                          :: trmzp_eq
  real(dp),     dimension(:),     allocatable       :: w1
  real(dp),     dimension(:),     allocatable       :: w2
  real(dp),     dimension(:),     allocatable       :: w3
  real(dp)                                          :: wk
  real(dp)                                          :: zpe
#ifdef TRACE
  call trace_in('vibenergy')
#endif
!
!  Set logicals
!
  lprinloc = (lprint.and.ioproc)
  lnozero  = (index(keyword,'noze').ne.0)
!
!  Allocate local pointer arrays
!
  allocate(rtmp2(mcv),stat=status)
  if (status/=0) call outofmemory('vibenergy','rtmp2')
  allocate(w1(mcv),stat=status)
  if (status/=0) call outofmemory('vibenergy','w1')
  allocate(w2(mcv),stat=status)
  if (status/=0) call outofmemory('vibenergy','w2')
  allocate(w3(mcv),stat=status)
  if (status/=0) call outofmemory('vibenergy','w3')
!
!  Set bounds on vibrational modes to consider
!
  if (ndim.eq.0) then
!
!  Cluster
!
    call moltype(linear)
    if (linear) then
      mcvmin = 6
    else
      mcvmin = 7
    endif
  else
!
!  Periodic system
!
    linear = .false.
    mcvmin = 1
    mcvmax = mcv
  endif
!
!  Check for imaginary modes and exclude them
!
  nimag = 0
  do i = 1,mcv
    if (freq(i,1).lt.-0.5_dp) nimag = nimag + 1
  enddo
  if (linear) then
    mcvmin = mcvmin + max(0,nimag-2)
  else
    mcvmin = mcvmin + max(0,nimag-3)
  endif
  mcvmax = mcv
  if (minmode.ne.1) mcvmin = minmode
  if (maxmode.ne.0) mcvmax = maxmode
!*************************************
!  Output phonon related properties  *
!*************************************
  if (ntemperatureramp.gt.0) then
    do ntr = 1,ntemperatureramp
      tem = temperaturestart(ntr)
!
!  For first ramp need to start from initial temperature
!  For subsequent ramps this would be the same as the end of the previous ramp
!
      if (ntr.eq.1) then
        nt0 = 0
      else
        nt0 = 1
      endif
      do nt = nt0,ntemperaturestep(ntr)
!
!  Set temperature 
!
        tem = temperaturestart(ntr) + dble(nt)*temperaturestep(ntr)
!***************************************
!  Evaluate phonon related properties  *
!***************************************
!
!  Initialise thermodynamic properties
!
        zpe = 0.0_dp
        entropy = 0.0_dp
        fhenergy = 0.0_dp
        kinenergy = 0.0_dp
        rtlnz = 0.0_dp
        cv = 0.0_dp
        fe_equipartition = 0.0_dp
        rkt = boltz*tem
        if (tem.gt.1.0d-6) then
          cmfact = planck*speedl/rkt
        else
          cmfact = 0.0_dp
        endif
!
!  Loop over K points
!
        do nk = 1,nllkpt
          wk = wkpt(nk)*rnokpt
          if (tem.gt.1.0d-6) then
!
!  Scale frequencies to hw/kT
!
            do i = 1,mcv
              rtmp2(i) = cmfact*freq(i,nk)
            enddo
!
!  Store exp(x) in w1, exp(x)-1 in w2 and 1/(exp(x)-1) in w3
!
            do i = 1,mcv
              if (rtmp2(i).lt.12.0_dp) then
                w1(i) = exp(rtmp2(i))
                w2(i) = w1(i) - 1.0_dp
                if (abs(w2(i)).gt.0.0_dp) w3(i) = 1.0_dp/w2(i)
              else
                w1(i) = exp(-rtmp2(i))
                w3(i) = exp(-rtmp2(i))
              endif
            enddo
!
!  Zero point energy
!
            factor = 0.5_dp*wk*rkt/evtoj
            trmzp = 0.0_dp
            do i = mcvmin,mcvmax
              if (rtmp2(i).gt.cmfact) trmzp = trmzp + rtmp2(i)
            enddo
            trmzp_eq = 0.5_dp*wk*trmzp
            trmzp = factor*trmzp
!
!  Kinetic energy
!
            do i = mcvmin,mcvmax
              kinenergy = kinenergy + factor*rtmp2(i)*(0.5_dp + w3(i))
            enddo
!
!  Entropy and free energy
!
            trmfe = 0.0_dp
            trmen = 0.0_dp
            freqmin = cmfact
            do i = mcvmin,mcvmax
              trm1 = rtmp2(i)
              if (trm1.gt.freqmin) then
                trm1 = 1.0_dp - exp(-trm1)
                trmfe = trmfe + log(trm1)
                trmen = trmen + rtmp2(i)*w3(i)
              endif
            enddo
!
!  Equipartition free energy
!
            trmfe_eq = 0.0_dp
            rmode = 0.0_dp
            do i = mcvmin,mcvmax
              trm1 = rtmp2(i)
              if (trm1.gt.freqmin) then
                trmfe_eq = trmfe_eq + log(trm1) - 1.0_dp
                rmode = rmode + 1.0_dp
              endif
            enddo
!
            factor = 2.0_dp*factor
            trm = factor*trmfe
            if (lnozero) then
              fhenergy = fhenergy + trm 
            else
              fhenergy = fhenergy + trm + trmzp
              zpe = zpe + trmzp
            endif
            rtlnz = rtlnz - trm
            factor = factor/tem
            entropy = entropy + factor*trmen
            if (lnozero) then
              fe_equipartition = fe_equipartition + wk*trmfe_eq - trmzp_eq
            else
              fe_equipartition = fe_equipartition + wk*trmfe_eq
            endif
!
!  Heat capacity - constant volume
!
            trmcv = 0.0_dp
            do i = mcvmin,mcvmax
              if (rtmp2(i).gt.freqmin) then
                if (rtmp2(i).lt.12.0_dp) then
                  trmcv = trmcv + rtmp2(i)*rtmp2(i)*w1(i)*w3(i)*w3(i)
                else
                  trmcv = trmcv + rtmp2(i)*rtmp2(i)*w3(i)
                endif
              endif
            enddo
            cv = cv + factor*trmcv
          elseif (.not.lnozero) then
!
!  Zero point energy
!
            factor = 0.5_dp*wk*planck*speedl/evtoj
            trmzp = 0.0_dp
            do i = mcvmin,mcvmax
              if (freq(i,nk).gt.1.0_dp) trmzp = trmzp + freq(i,nk)
            enddo
            trmzp = factor*trmzp
            zpe = zpe + trmzp
          endif
!
!  End of loop over K points
!
        enddo
        if (lprinloc) then
          if (ndim.eq.0) then
            write(ioout,'(''  Vibrational properties (for cluster):  Temperature  =  '',f10.3,'' K'')') tem
          else
            write(ioout,'(''  Phonon properties (per mole of unit cells): Temperature = '',f10.3,'' K'')') tem
          endif
          write(ioout,'(''--------------------------------------------------------------------------------'')')
          write(ioout,'(''  Zero point energy            = '',f15.6,'' eV'')') zpe
          if (tem.gt.1.0d-06) then
            trmen = fhenergy - zpe
            entropy = entropy - trmen/tem
            ent2 = entropy*evtoj*avogadro
            cv2 = cv*evtoj*avogadro
            write(ioout,'(''  Entropy                      = '',f15.6,'' eV/K'')') entropy
            write(ioout,'(''                               = '',f15.6,'' J/(mol.K)'')') ent2
            write(ioout,'(''  Helmholtz free-energy        = '',f15.6,'' eV'')') fhenergy+fc
            write(ioout,'(''                               = '',f15.6,'' kJmol-1'')') (fhenergy+fc)*evtoj*avogadro*0.001_dp
            write(ioout,'(''  Free energy (equipartition)  = '',f15.6,'' eV'')') fc + (fe_equipartition + rmode)*rkt/evtoj 
            write(ioout,'(''  - T*S       (equipartition)  = '',f15.6,'' eV'')') fe_equipartition*rkt/evtoj
            write(ioout,'(''  Uvib        (equipartition)  = '',f15.6,'' eV'')') rmode*rkt/evtoj
            write(ioout,'(''  Mean kinetic energy          = '',f15.6,'' eV'')') kinenergy
            write(ioout,'(''  Heat capacity - const volume = '',f15.6,'' eV/K'')') cv
            write(ioout,'(''                               = '',f15.6,'' J/(mol.K)'')') cv2
          endif
          write(ioout,'(''--------------------------------------------------------------------------------'')')
        endif
      enddo
    enddo
  else
!
!  Single temperature
!
    tem = temperature
!***************************************
!  Evaluate phonon related properties  *
!***************************************
!
!  Initialise thermodynamic properties
!
    zpe = 0.0_dp
    entropy = 0.0_dp
    fhenergy = 0.0_dp
    kinenergy = 0.0_dp
    rtlnz = 0.0_dp
    cv = 0.0_dp
    fe_equipartition = 0.0_dp
    rkt = boltz*tem
    if (tem.gt.1.0d-6) then
      cmfact = planck*speedl/rkt
    else
      cmfact = 0.0_dp
    endif
!
!  Loop over K points
!
    do nk = 1,nllkpt
      wk = wkpt(nk)*rnokpt
      if (tem.gt.1.0d-6) then
!
!  Scale frequencies to hw/kT
!
        do i = 1,mcv
          rtmp2(i) = cmfact*freq(i,nk)
        enddo
!
!  Store exp(x) in w1, exp(x)-1 in w2 and 1/(exp(x)-1) in w3
!
        do i = 1,mcv
          if (rtmp2(i).lt.12.0_dp) then
            w1(i) = exp(rtmp2(i))
            w2(i) = w1(i) - 1.0_dp
            if (abs(w2(i)).gt.0.0_dp) w3(i) = 1.0_dp/w2(i)
          else
            w1(i) = exp(-rtmp2(i))
            w3(i) = exp(-rtmp2(i))
          endif
        enddo
!
!  Zero point energy
!
        factor = 0.5_dp*wk*rkt/evtoj
        trmzp = 0.0_dp
        do i = mcvmin,mcvmax
          if (rtmp2(i).gt.cmfact) trmzp = trmzp + rtmp2(i)
        enddo
        trmzp_eq = 0.5_dp*wk*trmzp
        trmzp = factor*trmzp
!
!  Kinetic energy
!
        do i = mcvmin,mcvmax
          kinenergy = kinenergy + factor*rtmp2(i)*(0.5_dp + w3(i))
        enddo
!
!  Entropy and free energy
!
        trmfe = 0.0_dp
        trmen = 0.0_dp
        freqmin = cmfact
        do i = mcvmin,mcvmax
          trm1 = rtmp2(i)
          if (trm1.gt.freqmin) then
            trm1 = 1.0_dp - exp(-trm1)
            trmfe = trmfe + log(trm1)
            trmen = trmen + rtmp2(i)*w3(i)
          endif
        enddo
!
!  Equipartition free energy
!
        trmfe_eq = 0.0_dp
        rmode = 0.0_dp
        do i = mcvmin,mcvmax
          trm1 = rtmp2(i)
          if (trm1.gt.freqmin) then
            trmfe_eq = trmfe_eq + log(trm1) - 1.0_dp
            rmode = rmode + 1.0_dp
          endif
        enddo
!
        factor = 2.0_dp*factor
        trm = factor*trmfe
        if (lnozero) then
          fhenergy = fhenergy + trm 
        else
          fhenergy = fhenergy + trm + trmzp
          zpe = zpe + trmzp
        endif
        rtlnz = rtlnz - trm
        factor = factor/tem
        entropy = entropy + factor*trmen
        if (lnozero) then
          fe_equipartition = fe_equipartition + wk*trmfe_eq - trmzp_eq
        else
          fe_equipartition = fe_equipartition + wk*trmfe_eq
        endif
!
!  Heat capacity - constant volume
!
        trmcv = 0.0_dp
        do i = mcvmin,mcvmax
          if (rtmp2(i).gt.freqmin) then
            if (rtmp2(i).lt.12.0_dp) then
              trmcv = trmcv + rtmp2(i)*rtmp2(i)*w1(i)*w3(i)*w3(i)
            else
              trmcv = trmcv + rtmp2(i)*rtmp2(i)*w3(i)
            endif
          endif
        enddo
        cv = cv + factor*trmcv
      elseif (.not.lnozero) then
!
!  Zero point energy
!
        factor = 0.5_dp*wk*planck*speedl/evtoj
        trmzp = 0.0_dp
        do i = mcvmin,mcvmax
          if (freq(i,nk).gt.1.0_dp) trmzp = trmzp + freq(i,nk)
        enddo
        trmzp = factor*trmzp
        zpe = zpe + trmzp
      endif
!
!  End of loop over K points
!
    enddo
    if (lprinloc) then
      if (ndim.eq.0) then
        write(ioout,'(''  Vibrational properties (for cluster):  Temperature  =  '',f10.3,'' K'')') tem
      else
        write(ioout,'(''  Phonon properties (per mole of unit cells): Temperature = '',f10.3,'' K'')') tem
      endif
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''  Zero point energy            = '',f15.6,'' eV'')') zpe
      if (tem.gt.1.0d-06) then
        trmen = fhenergy - zpe
        entropy = entropy - trmen/tem
        ent2 = entropy*evtoj*avogadro
        cv2 = cv*evtoj*avogadro
        write(ioout,'(''  Entropy                      = '',f15.6,'' eV/K'')') entropy
        write(ioout,'(''                               = '',f15.6,'' J/(mol.K)'')') ent2
        write(ioout,'(''  Helmholtz free-energy        = '',f15.6,'' eV'')') fhenergy + fc
        write(ioout,'(''                               = '',f15.6,'' kJmol-1'')') (fhenergy+fc)*evtoj*avogadro*0.001_dp
        write(ioout,'(''  Free energy (equipartition)  = '',f15.6,'' eV'')') fc + (fe_equipartition + rmode)*rkt/evtoj 
        write(ioout,'(''  - T*S       (equipartition)  = '',f15.6,'' eV'')') fe_equipartition*rkt/evtoj
        write(ioout,'(''  Uvib        (equipartition)  = '',f15.6,'' eV'')') rmode*rkt/evtoj
        write(ioout,'(''  Mean kinetic energy          = '',f15.6,'' eV'')') kinenergy
        write(ioout,'(''  Heat capacity - const volume = '',f15.6,'' eV/K'')') cv
        write(ioout,'(''                               = '',f15.6,'' J/(mol.K)'')') cv2
      endif
      write(ioout,'(''--------------------------------------------------------------------------------'')')
    endif
  endif
!
!  Free local memory
!
  deallocate(w3,stat=status)
  if (status/=0) call deallocate_error('vibenergy','w3')
  deallocate(w2,stat=status)
  if (status/=0) call deallocate_error('vibenergy','w2')
  deallocate(w1,stat=status)
  if (status/=0) call deallocate_error('vibenergy','w1')
  deallocate(rtmp2,stat=status)
  if (status/=0) call deallocate_error('vibenergy','rtmp2')
#ifdef TRACE
  call trace_out('vibenergy')
#endif
!
  return
  end
