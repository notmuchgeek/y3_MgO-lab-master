  subroutine gulp_options
!
!  Main subroutine from which options are called
!
!   4/01 Created from main routine gulp 
!  11/02 Scan point number passed to optim
!   7/05 Call to eem added after optimisation to output final charges
!  11/06 Modifications for NEB added
!   1/07 Gasteiger charges added
!   5/07 Bond increment charges added
!   7/07 Metadynamics added
!   8/08 Metadynamics removed as option and introduced as a bias
!        to a standard MD run.
!   2/09 Modified to accommodate new version of FoX and gulp_cml
!   6/09 Renamed to gulp_options to avoid conflict in chemshell
!   4/10 Flag added to defect call to signal whether the defect
!        calculation is OK. 
!   4/10 QTPie charge transfer option added
!   9/11 lgrad1 argument added to eem call
!   8/13 Outfile only called by IO processor
!   1/15 Table of energies from translate steps added
!   8/15 Call to reaxffcharges added to output charges from reaxFF
!   8/15 Repeat calls to optim added for optlower case
!   9/15 translate noise added
!   7/17 Calls to eemdm added for parallel case
!   1/18 Trace added
!   2/18 Alamode interface added
!   2/18 Thermal noise added to translate option
!   3/18 Sign option added to translate
!   4/18 Scan_cell option added
!   4/18 Call to eemsplit added
!   5/18 Strain sub-option added to scan_cell
!   6/18 lscansuccess added
!   6/18 Strain scan now uses strain cell option
!   6/18 Format of output banner corrected
!  10/18 Handling of units for energy added
!   1/19 Size of confignum string adjusted
!   7/20 lkcal and lkjmol now accessed from module
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
!  Copyright Curtin University 2020
!
!  Julian Gale, CIC, Curtin University, July 2020
!
  use configurations
  use control
  use current
  use defects
  use element,         only : lgasteiger
  use g_constants,     only : kjmtoev, kcaltoev
  use general
  use genetic
  use gulp_cml,        only : lcml, gulp_cml_StartModule, gulp_cml_EndModule, gulp_cml_structure, &
                              gulp_cml_startstep, gulp_cml_endstep
  use iochannels
  use mdlogic
  use m_alamode,       only : gen_fc_alamode, anharmonic_alamode
  use optimisation
  use parallel
  use potentialgrid
  use potentialsites,  only : npotsitescfg
  use scan
#ifdef TRACE
  use trace,           only : trace_in, trace_out
#endif

  implicit none
!
!  Local variables
!
  character(len=8)                             :: confignum
  character(len=8)                             :: unitword
  integer(i4)                                  :: i
  integer(i4)                                  :: ic
  integer(i4)                                  :: ii
  integer(i4)                                  :: mc
  integer(i4)                                  :: mcfg
  integer(i4)                                  :: n
  integer(i4)                                  :: nsc
  integer(i4)                                  :: ntat
  integer(i4)                                  :: ntr
  integer(i4)                                  :: status
  logical                                      :: ldefectOK
  logical,     allocatable, dimension(:), save :: lscansuccess
  logical                                      :: lselect
  logical                                      :: lnofp
  real(dp)                                     :: rn
  real(dp),    allocatable, dimension(:), save :: ecscan
  real(dp),    allocatable, dimension(:), save :: etranslate
  real(dp)                                     :: GULP_grandom
  real(dp)                                     :: GULP_random
  real(dp)                                     :: cstep(6)
  real(dp)                                     :: tsign
  real(dp)                                     :: units
  real(dp)                                     :: xnoise
  real(dp)                                     :: ynoise
  real(dp)                                     :: znoise
  real(dp)                                     :: xtr
  real(dp)                                     :: ytr
  real(dp)                                     :: ztr
#ifdef TRACE
  call trace_in('gulp_options')
#endif
!
  if (lkcal) then
    units = 1.0_dp/kcaltoev
    unitword = '(kcal)  '
  elseif (lkjmol) then
    units = 1.0_dp/kjmtoev
    unitword = '(kJ/mol)'
  else
    units = 1.0_dp
    unitword = '(eV)    '
  endif
!***********************
!  Configuration loop  *
!***********************
  do ic = 1,ncfg
    ncf = ic
!
!  Put nreg1 to zero so that dumpdur knows not to dump
!  explicit region 1 until defect calc has started.
!
    nreg1 = 0
    if (ioproc) then
      write(ioout,'(/,''********************************************************************************'')')
      if (names(ncf)(1:1).ne.' ') then
        write(ioout,'(''*  Output for configuration '',i3,'' : '',a45,''*'')') ncf,names(ncf)(1:45)
      else
        write(ioout,'(''*  Output for configuration '',i3,48x,''*'')') ncf
      endif
      write(ioout,'(''********************************************************************************'',/)')
    endif
!
!  New CML
!  Start a module for this configuration
!
    if (lcml) then
      write(confignum,'(i8)') ncf
      call gulp_cml_StartModule(title='Gulp configuration')
    endif

!*************************
!  Structure prediction  *
!*************************
    mcfg = 1
    if (lpredict.and.ndimen(ic).eq.3) then
      call setup(.true.)
      call predict(mcfg)
    endif
!**************************************************************
!  Loop over the number of configurations found from predict  *
!**************************************************************
    do mc = 1,mcfg
      if (lpredict.and.ndimen(ic).eq.3) then
        ncf = ic
        call setup(.false.)
        call gabcfg(mc,mcfg)
      endif
      if (ncscan(ncf).gt.0) then
!
!  Allocate space to store final energies
!
        allocate(ecscan(0:ncscan(ncf)),stat=status)
        if (status/=0) call outofmemory('options','ecscan')
        allocate(lscansuccess(0:ncscan(ncf)),stat=status)
        if (status/=0) call outofmemory('options','lscansuccess')
!*******************
!  Scan over cell  *
!*******************
        call setup(.true.)
        nsc = ncscan(ncf)
        lnofp = (index(keyword,'nofi').ne.0)
!
!  Set the vector steps
!
        rn = 1.0_dp/dble(nsc)
        if (ndimen(ncf).eq.3) then
          cstep(1:6) = cscan(1:6,ncf)*rn
        elseif (ndimen(ncf).eq.2) then
          cstep(1:3) = cscan(1:3,ncf)*rn
        elseif (ndimen(ncf).eq.1) then
          cstep(1) = cscan(1,ncf)*rn
        endif
        if (ioproc) then
          write(ioout,'(/,''  Repeat cell scan for '',i6,'' points'')') nsc+1
          if (lcscanstrain(ncf)) then
            write(ioout,'(/,''  Scan vector for cell strains:'',/)')
            if (ndimen(ncf).eq.3) then
              write(ioout,'(''  1(xx) = '',f12.6,''  2(yy) = '',f12.6,''  3(zz) = '',f12.6)') &
                cscan(1,ncf),cscan(2,ncf),cscan(3,ncf)
              write(ioout,'(''  4(yz) = '',f12.6,''  5(xz) = '',f12.6,''  6(xy) = '',f12.6,/)') &
                cscan(4,ncf),cscan(5,ncf),cscan(6,ncf)
            elseif (ndimen(ncf).eq.2) then
              write(ioout,'(''  1(xx) = '',f12.6,''  2(yy) = '',f12.6,''  3(xy) = '',f12.6,/)') &
                cscan(1,ncf),cscan(2,ncf),cscan(3,ncf)
            elseif (ndimen(ncf).eq.1) then
              write(ioout,'(''  1(xx) = '',f12.6,/)') cscan(1,ncf)
            endif
          else
            write(ioout,'(/,''  Scan vector for cell parameters:'',/)')
            if (ndimen(ncf).eq.3) then
              write(ioout,'(''  a     = '',f12.6,''  b    = '',f12.6,''  c     = '',f12.6)') &
                cscan(1,ncf),cscan(2,ncf),cscan(3,ncf)
              write(ioout,'(''  alpha = '',f12.6,''  beta = '',f12.6,''  gamma = '',f12.6,/)') &
                cscan(4,ncf),cscan(5,ncf),cscan(6,ncf)
            elseif (ndimen(ncf).eq.2) then
              write(ioout,'(''  a     = '',f12.6,''  b    = '',f12.6,''  alpha = '',f12.6,/)') &
                cscan(1,ncf),cscan(2,ncf),cscan(3,ncf)
            elseif (ndimen(ncf).eq.1) then
              write(ioout,'(''  a     = '',f12.6,/)') cscan(1,ncf)
            endif
          endif
        endif
!**************************
!  Loop over cell points  *
!**************************
        do n = 0,nsc
!
! New CML
! Start a step for the point of the scan
! and add an Initialisation module
!
          if (lcml) then
             Call gulp_cml_startstep(cycle=n, type="gulp:energyscan")
             Call gulp_cml_StartModule(title='Initialisation', dictRef='eminerals:init')
             Call gulp_cml_structure(ncf)
             Call gulp_cml_endModule
          endif
          if (.not.lnofp.or.n.ne.0) then
            if (ioproc) then
              write(ioout,'(/,''################################################################################'')')
              write(ioout,'(''#  Output for cell point '',i5,49x,''#'')') n
              write(ioout,'(''################################################################################'',/)')
            endif
            lselect = .false.
            if (lrelax.and..not.lopt) then
              if (ioproc) then
                write(ioout,'(/,''  Properties for fitted structure :'',/)')
              endif
            endif
            if (lqbond) then
              lselect = .true.
              call setup(.true.)
              call bondq(.true.)
            endif
            if (lgasteiger) then
              lselect = .true.
              call setup(.true.)
              call gasteiger(.true.)
            endif
            if (leem) then
              lselect = .true.
              call setup(.true.)
              if (lqtpie) then
                call qtpie(.true.)
              else
                if (nprocs.gt.1) then
                  if (leembond) then
                    call eemsplitd(.true.,.false.,.false.)
                  else
                    call eemdm(.true.,.false.,.false.)
                  endif
                else
                  if (leembond) then
                    call eemsplit(.true.,.false.,.false.)
                  else
                    call eem(.true.,.false.,.false.)
                  endif
                endif
              endif
            endif
            if (.not.lnoenergy) then
              lselect = .true.
              call setup(.true.)
              call optim(.true.,.true.,n)
!
!  Lower repeat optimisation if needed
!
              if (loptlower.and.llowered(ncf)) then
                if (ioproc) then
                  write(ioout,'(/,''--------------------------------------------------------------------------------'')')
                  write(ioout,'(''|  Repeat optimisation after lower for point '',i5,29x,''|'')') n
                  write(ioout,'(''--------------------------------------------------------------------------------'',/)')
                endif
                call optim(.true.,.true.,n)
!
!  Reset configuration lower flag for the benefit of translate
!
                llowered(ncf) = .false.
              endif
!
!  Save energy for configuration
!
              ecscan(n) = energycfg(ncf)
              lscansuccess(n) = loptsuccess
              if (leem) then
                call setup(.true.)
                if (lqtpie) then
                  call qtpie(.true.)
                else
                  if (nprocs.gt.1) then
                    if (leembond) then
                      call eemsplitd(.true.,.false.,.false.)
                    else
                      call eemdm(.true.,.false.,.false.)
                    endif
                  else
                    if (leembond) then
                      call eemsplit(.true.,.false.,.false.)
                    else
                      call eem(.true.,.false.,.false.)
                    endif
                  endif
                endif
              elseif (lreaxFF.and.lreaxFFQ) then
                call reaxffcharges
              endif
            endif
            if (lalamode) then
              lselect = .true.
              call setup(.true.)
              call gen_fc_alamode
              call anharmonic_alamode
            endif
            if (lpot.or.lefg) then
              lselect = .true.
              call setup(.true.)
              call latpot
            endif
            if (npotsitescfg(ic).gt.0) then
              lselect = .true.
              call setup(.true.)
              call sitepot
            endif
            if (lneb) then 
              lselect = .true.
              call setup(.true.)
              if (index(keyword,'sync').ne.0) then
                call runsync
              else
                call runneb
              endif
            endif
            if (ldefect.and.ndimen(ic).eq.3) then
!$$$$$$$$$$$$$$$$$$$$$$$
!  Defect calculation  $
!$$$$$$$$$$$$$$$$$$$$$$$
              lselect = .true.
              call setup(.true.)
              if (.not.lnoenergy) call setdr2mat
!
!  Setup defect regions
!
              call defect(ldefectOK)
              if (ldefectOK) then
!
!  Function call
!
                if (.not.lnoenergy) then
                  call defopt(.true.)
                endif
!
!  Defect site potentials
!
                if (lpot) call defpot
!
!  Move region 2a ions into region with final displacements
!
                if (reg2a1(ncf).gt.0.0_dp) call move2a1
              endif
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!  End of defect calculation  $
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            endif
            if (.not.lselect) then
              if (ioproc) then
                write(ioout,'(''  **** No options selected for configuration '',i3,''  ****'')')ncf
              endif
            endif
          endif
!
!  Shift cell for next step - if vector has no component specified then preserve value from last optimisation step.
!
          if (n.ne.nsc) then
            if (lcscanstrain(ncf)) then
!################
!  Cell strain  #
!################
              straincfg(1:nstrains,ncf) = straincfg(1:nstrains,ncf) + cstep(1:nstrains)
            else
!####################
!  Cell parameters  #
!####################
              if (ndimen(ncf).eq.3) then
                a = a + cstep(1)
                b = b + cstep(2)
                c = c + cstep(3)
                if (a.lt.0.0_dp.or.b.lt.0.0_dp.or.c.lt.0.0_dp) then
                  call outerror('cell length has gone below zero',0_i4)
                  call stopnow('options')
                endif
                alpha = alpha + cstep(4)
                beta  = beta  + cstep(5)
                gamma = gamma + cstep(6)
                if (alpha.lt.0.0_dp.or.beta.lt.0.0_dp.or.gamma.lt.0.0_dp) then
                  call outerror('cell angle has gone below zero',0_i4)
                  call stopnow('options')
                endif
                call cell3D(rv,a,b,c,alpha,beta,gamma)
                do i = 1,3
                  rvcfg(1:3,i,ncf) = rv(1:3,i)
                enddo
              elseif (ndimen(ncf).eq.2) then
                a = a + cstep(1)
                b = b + cstep(2)
                if (a.lt.0.0_dp.or.b.lt.0.0_dp) then
                  call outerror('cell length has gone below zero',0_i4)
                  call stopnow('options')
                endif
                alpha = alpha + cstep(3)
                if (alpha.lt.0.0_dp) then
                  call outerror('cell angle has gone below zero',0_i4)
                  call stopnow('options')
                endif
                call cell2D(rv,a,b,alpha)
                do i = 1,2
                  rvcfg(1:3,i,ncf) = rv(1:3,i)
                enddo
              elseif (ndimen(ncf).eq.1) then
                a = a + cstep(1)
                if (a.lt.0.0_dp) then
                  call outerror('cell length has gone below zero',0_i4)
                  call stopnow('options')
                endif 
                call cell1D(rv,a)
                rvcfg(1:3,1,ncf) = rv(1:3,1)
              endif
            endif
          endif
! 
!  New CML
!  Having done that step close it
!
          if (lcml) call gulp_cml_endstep
        enddo
!
!  Output table of final energies
!
        if (ioproc.and..not.lnoenergy) then
          write(ioout,'(/,''  Final energies for configurations from scan_cell option:'',/)')
          write(ioout,'(''--------------------------------------------------------------------------------'')')
          if (lopt) then
            write(ioout,'(''    Cell scan step number                           Energy '',a8,''  Opt OK?    '')') unitword
          else
            write(ioout,'(''    Cell scan step number                           Energy '',a8,''             '')') unitword
          endif
          write(ioout,'(''--------------------------------------------------------------------------------'')')
          if (lopt) then
            do n = 0,nsc
              write(ioout,'(2x,i12,30x,f20.8,6x,l5)') n,ecscan(n)*units,lscansuccess(n)
            enddo
          else
            do n = 0,nsc
              write(ioout,'(2x,i12,30x,f20.8)') n,ecscan(n)*units
            enddo
          endif
          write(ioout,'(''--------------------------------------------------------------------------------'',/)')
        endif
!
!  Deallocate space to store final energies
!
        deallocate(lscansuccess,stat=status)
        if (status/=0) call deallocate_error('options','lscansuccess')
        deallocate(ecscan,stat=status)
        if (status/=0) call deallocate_error('options','ecscan')
      elseif (ntran(ncf).gt.0) then
!
!  Allocate space to store final energies
!
        allocate(etranslate(0:ntran(ncf)),stat=status)
        if (status/=0) call outofmemory('options','etranslate')
!********************************
!  Scan over atom translations  *
!********************************
        call setup(.true.)
        ntr = ntran(ncf)
        lnofp = (index(keyword,'nofi').ne.0)
!
!  Set the vector steps
!
        rn = 1.0_dp/dble(ntr)
        xtr = xtran(ncf)*rn
        ytr = ytran(ncf)*rn
        ztr = ztran(ncf)*rn
!
!  Count number of atoms being translated
!
        ntat = 0
        do i = 1,nasym
          if (ltranat(i+nsft)) then
            ntat = ntat + 1
          endif
        enddo
!
!  Check that there some atoms to be translated
!
        if (ntat.eq.0) then
          call outerror('no atoms specified for translation',0_i4)
          call stopnow('options')
        endif
        if (ioproc) then
          write(ioout,'(/,''  Repeat calculation for '',i3,'' points'')') ntr+1
          write(ioout,'(/,''  Translation vector :'',/)')
          write(ioout,'(''  x = '',f12.6,''  y = '',f12.6,''  z = '',f12.6,/)') xtran(ncf),ytran(ncf),ztran(ncf)
          if (ltrantherm(ncf)) then
            write(ioout,'(''  Thermal noise will be added during translation up to '',f8.4,'' Ang'')') trantherm(ncf)
          elseif (ltrannoise(ncf)) then
            write(ioout,'(''  Random noise will be added to displacements scaled by '',f8.6)') trannoise(ncf)
          endif
        endif
!*********************
!  Loop over points  *
!*********************
        do n = 0,ntr
!
! New CML
! Start a step for the point of the scan
! and add an Initialisation module
!
          if (lcml) then
             Call gulp_cml_startstep(cycle=n, type="gulp:energyscan")
             Call gulp_cml_StartModule(title='Initialisation', dictRef='eminerals:init')
             Call gulp_cml_structure(ncf)
             Call gulp_cml_endModule
          endif
          if (.not.lnofp.or.n.ne.0) then
            if (ioproc) then
              write(ioout,'(/,''################################################################################'')')
              write(ioout,'(''#  Output for point '',i5,54x,''#'')')n
              write(ioout,'(''################################################################################'',/)')
            endif
            lselect = .false.
            if (lrelax.and..not.lopt) then
              if (ioproc) then
                write(ioout,'(/,''  Properties for fitted structure :'',/)')
              endif
            endif
            if (lqbond) then
              lselect = .true.
              call setup(.true.)
              call bondq(.true.)
            endif
            if (lgasteiger) then
              lselect = .true.
              call setup(.true.)
              call gasteiger(.true.)
            endif
            if (leem) then
              lselect = .true.
              call setup(.true.)
              if (lqtpie) then
                call qtpie(.true.)
              else
                if (nprocs.gt.1) then
                  if (leembond) then
                    call eemsplitd(.true.,.false.,.false.)
                  else
                    call eemdm(.true.,.false.,.false.)
                  endif
                else
                  if (leembond) then
                    call eemsplit(.true.,.false.,.false.)
                  else
                    call eem(.true.,.false.,.false.)
                  endif
                endif
              endif
            endif
            if (.not.lnoenergy) then
              lselect = .true.
              call setup(.true.)
              call optim(.true.,.true.,n)
!
!  Lower repeat optimisation if needed
!
              if (loptlower.and.llowered(ncf)) then
                if (ioproc) then
                  write(ioout,'(/,''--------------------------------------------------------------------------------'')')
                  write(ioout,'(''|  Repeat optimisation after lower for point '',i5,29x,''|'')') n
                  write(ioout,'(''--------------------------------------------------------------------------------'',/)')
                endif
                call optim(.true.,.true.,n)
!
!  Reset configuration lower flag for the benefit of translate
!
                llowered(ncf) = .false.
              endif
!
!  Save energy for configuration
!
              etranslate(n) = energycfg(ncf)
              if (leem) then
                call setup(.true.)
                if (lqtpie) then
                  call qtpie(.true.)
                else
                  if (nprocs.gt.1) then
                    if (leembond) then
                      call eemsplitd(.true.,.false.,.false.)
                    else
                      call eemdm(.true.,.false.,.false.)
                    endif
                  else
                    if (leembond) then
                      call eemsplit(.true.,.false.,.false.)
                    else
                      call eem(.true.,.false.,.false.)
                    endif
                  endif
                endif
              elseif (lreaxFF.and.lreaxFFQ) then
                call reaxffcharges
              endif
            endif
            if (lalamode) then
              lselect = .true.
              call setup(.true.)
              call gen_fc_alamode
              call anharmonic_alamode
            endif
            if (lpot.or.lefg) then
              lselect = .true.
              call setup(.true.)
              call latpot
            endif
            if (npotsitescfg(ic).gt.0) then
              lselect = .true.
              call setup(.true.)
              call sitepot
            endif
            if (lneb) then 
              lselect = .true.
              call setup(.true.)
              if (index(keyword,'sync').ne.0) then
                call runsync
              else
                call runneb
              endif
            endif
            if (ldefect.and.ndimen(ic).eq.3) then
!$$$$$$$$$$$$$$$$$$$$$$$
!  Defect calculation  $
!$$$$$$$$$$$$$$$$$$$$$$$
              lselect = .true.
              call setup(.true.)
              if (.not.lnoenergy) call setdr2mat
!
!  Setup defect regions
!
              call defect(ldefectOK)
              if (ldefectOK) then
!
!  Function call
!
                if (.not.lnoenergy) then
                  call defopt(.true.)
                endif
!
!  Defect site potentials
!
                if (lpot) call defpot
!
!  Move region 2a ions into region with final displacements
!
                if (reg2a1(ncf).gt.0.0_dp) call move2a1
              endif
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!  End of defect calculation  $
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            endif
            if (.not.lselect) then
              if (ioproc) then
                write(ioout,'(''  **** No options selected for configuration '',i3,''  ****'')')ncf
              endif
            endif
          endif
!
!  Shift coordinates for next point - if vector has
!  no component in a particular cartesian direction
!  then preserve value from last optimisation step.
!
          if (n.ne.ntr) then
            ii = 0
            do i = 1,nasym
              if (ltranat(i+nsft)) then
                if (ltranatminus(i+nsft)) then
                  tsign = -1.0_dp
                else
                  tsign = 1.0_dp
                endif
                ii = ii + 1
                if (ltrannoise(ncf).and..not.lregionrigid(nregionno(nsft+i),ncf)) then
!
!  Only apply random shifts to non-rigid region atoms
!
                  xnoise = GULP_random(iseed,2_i4)
                  ynoise = GULP_random(iseed,2_i4)
                  znoise = GULP_random(iseed,2_i4)
                  xcfg(i+nsft) = xcfg(i+nsft) + tsign*xtr*(1.0_dp + xnoise*trannoise(ncf))
                  ycfg(i+nsft) = ycfg(i+nsft) + tsign*ytr*(1.0_dp + ynoise*trannoise(ncf))
                  zcfg(i+nsft) = zcfg(i+nsft) + tsign*ztr*(1.0_dp + znoise*trannoise(ncf))
                else
                  xcfg(i+nsft) = xcfg(i+nsft) + tsign*xtr
                  ycfg(i+nsft) = ycfg(i+nsft) + tsign*ytr
                  zcfg(i+nsft) = zcfg(i+nsft) + tsign*ztr
                endif
              endif
              if (ltrantherm(ncf).and..not.lregionrigid(nregionno(nsft+i),ncf)) then
!
!  Apply thermal noise to all atoms that are not in a rigid region
!
                xnoise = GULP_grandom(iseed,trantherm(ncf))
                ynoise = GULP_grandom(iseed,trantherm(ncf))
                znoise = GULP_grandom(iseed,trantherm(ncf))
                if (ndimen(ncf).eq.3) then
                  xcfg(i+nsft) = xcfg(i+nsft) + xnoise/a
                  ycfg(i+nsft) = ycfg(i+nsft) + ynoise/b
                  zcfg(i+nsft) = zcfg(i+nsft) + znoise/c
                elseif (ndimen(ncf).eq.2) then
                  xcfg(i+nsft) = xcfg(i+nsft) + xnoise/a
                  ycfg(i+nsft) = ycfg(i+nsft) + ynoise/b
                  zcfg(i+nsft) = zcfg(i+nsft) + znoise
                elseif (ndimen(ncf).eq.1) then
                  xcfg(i+nsft) = xcfg(i+nsft) + xnoise/a
                  ycfg(i+nsft) = ycfg(i+nsft) + ynoise
                  zcfg(i+nsft) = zcfg(i+nsft) + znoise
                else
                  xcfg(i+nsft) = xcfg(i+nsft) + xnoise
                  ycfg(i+nsft) = ycfg(i+nsft) + ynoise
                  zcfg(i+nsft) = zcfg(i+nsft) + znoise
                endif
              endif
            enddo
          endif
! 
!  New CML
!  Having done that step close it
!
          if (lcml) call gulp_cml_endstep
        enddo
!
!  Output table of final energies
!
        if (ioproc.and..not.lnoenergy) then
          write(ioout,'(/,''  Final energies for configurations from translate option:'',/)')
          write(ioout,'(''--------------------------------------------------------------------------------'')')
          write(ioout,'(''    Translate step number                           Energy '',a8,''             '')') unitword
          write(ioout,'(''--------------------------------------------------------------------------------'')')
          do n = 0,ntr
            write(ioout,'(2x,i12,30x,f20.8)') n,etranslate(n)*units
          enddo
          write(ioout,'(''--------------------------------------------------------------------------------'',/)')
        endif
!
!  Deallocate space to store final energies
!
        deallocate(etranslate,stat=status)
        if (status/=0) call deallocate_error('options','etranslate')
      else
!***********************************
!  Single pass over configuration  *
!***********************************
! 
!  New CML 
!  For this single pass write an Initialisation module
!
        if (lcml) then
          write(confignum,'(i8)') ncf
          Call gulp_cml_StartModule(title='Initialisation', dictRef='eminerals:init')
          Call gulp_cml_structure(ncf)
          Call gulp_cml_endModule
        endif
        lselect = .false.
        if (lrelax.and..not.lopt) then
          if (ioproc) then
            write(ioout,'(/,''  Properties for fitted structure :'',/)')
          endif
        endif
        if (lqbond) then
          lselect = .true.
          call setup(.true.)
          call bondq(.true.)
        endif
        if (lgasteiger) then
          lselect = .true.
          call setup(.true.)
          call gasteiger(.true.)
        endif
        if (leem) then
          lselect = .true.
          call setup(.true.)
          if (lqtpie) then
            call qtpie(.true.)
          else
            if (nprocs.gt.1) then
              if (leembond) then
                call eemsplitd(.true.,.false.,.false.)
              else
                call eemdm(.true.,.false.,.false.)
              endif
            else
              if (leembond) then
                call eemsplit(.true.,.false.,.false.)
              else
                call eem(.true.,.false.,.false.)
              endif
            endif
          endif
        endif
        if (.not.lnoenergy) then
          lselect = .true.
          call setup(.true.)
          call optim(.true.,.false.,0_i4)
!
!  Lower repeat optimisation if needed
!
          if (loptlower.and.llowered(ncf)) then
            if (ioproc) then
              write(ioout,'(/,''--------------------------------------------------------------------------------'')')
              write(ioout,'(''|  Repeat optimisation after lower for configuration '',i5,21x,''|'')') ncf
              write(ioout,'(''--------------------------------------------------------------------------------'',/)')
            endif
            call optim(.true.,.false.,0_i4)
          endif
          if (leem) then
            call setup(.true.)
            if (lqtpie) then
              call qtpie(.true.)
            else
              if (nprocs.gt.1) then
                if (leembond) then
                  call eemsplitd(.true.,.false.,.false.)
                else
                  call eemdm(.true.,.false.,.false.)
                endif
              else
                if (leembond) then
                  call eemsplit(.true.,.false.,.false.)
                else
                  call eem(.true.,.false.,.false.)
                endif
              endif
            endif
          elseif (lreaxFF.and.lreaxFFQ) then
            call reaxffcharges
          endif
        endif
        if (lalamode) then
          lselect = .true.
          call setup(.true.)
          call gen_fc_alamode
          call anharmonic_alamode
        endif
        if (lpot.or.lefg) then
          lselect = .true.
          call setup(.true.)
          call latpot
        endif
        if (npotsitescfg(ic).gt.0) then
          lselect = .true.
          call setup(.true.)
          call sitepot
        endif
        if (nxpg(ic).ne.0) then
          lselect = .true.
          call setup(.true.)
          call potgrid
        endif
        if (lneb) then
          lselect = .true.
          call setup(.true.)
          if (index(keyword,'sync').ne.0) then
            call runsync
          else
            call runneb
          endif
        endif
        if (lmc) then
          lselect = .true.
          call runmc
        endif
        if (lmd) then
          lselect = .true.
          call runmd
        endif
        if (ldefect.and.ndimen(ic).eq.3) then
!$$$$$$$$$$$$$$$$$$$$$$$
!  Defect calculation  $
!$$$$$$$$$$$$$$$$$$$$$$$
          lselect = .true.
          call setup(.true.)
          if (.not.lnoenergy) call setdr2mat
!
!  Setup defect regions
!
          call defect(ldefectOK) 
          if (ldefectOK) then
!
!  Function call
!
            if (.not.lnoenergy) then
              call defopt(.true.)
            endif
!
!  Defect site potentials
!
            if (lpot) call defpot
!
!  Move region 2a ions into region with final displacements
!
            if (reg2a1(ncf).gt.0.0_dp) call move2a1
          endif
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!  End of defect calculation  $
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        endif
        if (.not.lselect) then
          if (ioproc) then
            write(ioout,'(''  **** No options selected for configuration '',i3,'' ****'')') ncf
          endif
        endif
      endif
    enddo
    if (ioproc.and.lpredict) call outfile
!
! New CML
! Close the module for the configuration.
! (Finalisation etc written by the called 
! subs I'm afraid - I can not find a way a
! around this).
    if (lcml) then
      Call gulp_cml_EndModule
    endif
  enddo
#ifdef TRACE
  call trace_out('gulp_options')
#endif
!
  return
  end
