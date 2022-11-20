  module m_plumed
!
!  Module that contains subroutines that communicate with PLUMED-2
!
!   1/19 maxwordlength changes added
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
!  Julian Gale, Curtin University, January 2019
!

  use datatypes
  use general,     only : time0
  use gulp_lengths
  use interupt,    only : controlC_fit, sigint
  use iochannels
  use parallel,    only : ioproc

  implicit none

  character(len=maxwordlength), save :: plumedinput       ! Name of plumed input file
  character(len=maxwordlength), save :: plumedlog         ! Name of plumed log file
  logical,                      save :: lplumed           ! If true then call plumed
  logical,                      save :: lplumed_available ! If true then plumed can be called
  real(dp),                     save :: tplumed = 0.0_dp  ! Time used by plumed

  CONTAINS

  subroutine init_plumed
!
!  Initialises PLUMED at the start of a run
!
!   3/16 Created from runmd
!
!  For now Plumed is called from a single node in serial mode
!
#ifdef PLUMED
  use configurations, only : rvcfg, xcfg, ycfg, zcfg, radcfg, ndimen
#endif
  use control
  use current
#ifdef PLUMED
  use g_constants,    only : kjmtoev
#endif
  use general
  use moldyn
  use parallel
  implicit none
#ifdef PLUMED
!
!  Local variables
!
  integer(i4)                               :: nplumedavailable
  real*8                                    :: evtokj
!
  if (lplumed) then
!**************************
!  Is PLUMED available ?  *
!**************************
    call plumed_f_installed(nplumedavailable)
    if (nplumedavailable.le.0) then
      call outerror('PLUMED is requested but not available',0_i4)
      call stopnow('init_plumed')
    endif
!**************************
!  PLUMED initialisation  *
!**************************
!
!  Set up values for passing to Plumed
!
    evtokj = 1.0_dp/kjmtoev
    if (ioproc) then
!
!  Create the Plumed object
!
      call plumed_f_gcreate()
!
!  Pass important data to Plumed
!
      call plumed_f_gcmd("setRealPrecision"//char(0),8)                      ! Pass real precision number
      call plumed_f_gcmd("setMDEnergyUnits"//char(0),evtokj)                 ! Pass energy conversion factor to kJ/mol
      call plumed_f_gcmd("setMDLengthUnits"//char(0),0.1d0)                  ! Pass length conversion factor to nm
      call plumed_f_gcmd("setMDTimeUnits"//char(0),1.0d0)                    ! Pass time conversion factor to ps
      call plumed_f_gcmd("setPlumedDat"//char(0),trim(plumedinput)//char(0)) ! Plumed input file name
      call plumed_f_gcmd("setLogFile"//char(0),trim(plumedlog)//char(0))     ! Plumed log file name
!      if (nprocs.gt.1) then
!        call plumed_f_gcmd("setMPIFComm"//char(0),MPI_comm_GULP)             ! MPI communicator
!      endif
      call plumed_f_gcmd("setNatoms"//char(0),numat)                         ! Number of atoms
      call plumed_f_gcmd("setMDEngine"//char(0),"gulp")                      ! Pass name of MD engine as a label
      call plumed_f_gcmd("setTimestep"//char(0),tstep)                       ! Timestep
!
!  Initialise Plumed
!
      call plumed_f_gcmd("init"//char(0),0)
    endif
  endif
#endif

  return
  end subroutine init_plumed

  subroutine finalise_plumed
!
!  Stops PLUMED at the end of a run
!
!   3/16 Created from init_plumed
!
!  For now Plumed is called from a single node in serial mode
!
  implicit none
!
#ifdef PLUMED
!**************************
!  PLUMED initialisation  *
!**************************
  if (lplumed) then
!
!  Destroy the Plumed object
!
    if (ioproc) then
      call plumed_f_gfinalize()
    endif
  endif
#endif

  return
  end subroutine finalise_plumed

  subroutine plumed_energy_force(istep,etot,xdrv,ydrv,zdrv)
!
!  Computes the metadynamics contribution to the energy & forces
!
!   3/16 Created 
!   4/16 Sign of forces corrected
!
!  For now Plumed is called from a single node in serial mode
!
#ifdef PLUMED
  use current,       only : ndim, rv, numat
  use current,       only : xclat, yclat, zclat, qf, mass
  use derivatives,   only : virial
  use m_pr,          only : virial_m
  use symmetry,      only : lstr
#endif
  implicit none
!
!  Passed variables
!
  integer(i4),  intent(in)    :: istep
  real(dp),     intent(inout) :: etot
  real(dp),     intent(inout) :: xdrv(*)
  real(dp),     intent(inout) :: ydrv(*)
  real(dp),     intent(inout) :: zdrv(*)
#ifdef PLUMED
!
!  Local variables
!
  integer(i4)                 :: i
  integer(i4)                 :: j
  real*8,                save :: pvirial(3,3)
  real*8, allocatable,   save :: pforcex(:)
  real*8, allocatable,   save :: pforcey(:)
  real*8, allocatable,   save :: pforcez(:)
  real(dp)                    :: t1
  real(dp)                    :: t2
  real(dp)                    :: g_cpu_time
!************************************
!  PLUMED metadynamics calculation  *
!************************************
  if (lplumed) then
    pvirial(1:3,1:3) = 0.0d0
    if (ioproc) then
      t1 = g_cpu_time()
!
!  Initialise force arrays
!
      allocate(pforcex(numat))
      allocate(pforcey(numat))
      allocate(pforcez(numat))
      pforcex(1:numat) = 0.0d0
      pforcey(1:numat) = 0.0d0
      pforcez(1:numat) = 0.0d0
!
!  Pass current data
!
      call plumed_f_gcmd("setStep"//char(0),istep)                     ! Time step number
      call plumed_f_gcmd("setPositionsX"//char(0),xclat)               ! X coordinates of atoms
      call plumed_f_gcmd("setPositionsY"//char(0),yclat)               ! Y coordinates of atoms
      call plumed_f_gcmd("setPositionsZ"//char(0),zclat)               ! Z coordinates of atoms
      call plumed_f_gcmd("setMasses"//char(0),mass)                    ! Masses of atoms
      call plumed_f_gcmd("setCharges"//char(0),qf)                     ! Charges of atoms
      if (ndim.eq.3_i4) then
        call plumed_f_gcmd("setBox"//char(0),rv)                       ! Unit cell
        call plumed_f_gcmd("setVirial"//char(0),pvirial)               ! Virial
      endif
!
!  Pass quantities that Plumed needs to modify
!
      call plumed_f_gcmd("setEnergy"//char(0),etot)
      call plumed_f_gcmd("setForcesX"//char(0),pforcex)
      call plumed_f_gcmd("setForcesY"//char(0),pforcey)
      call plumed_f_gcmd("setForcesZ"//char(0),pforcez)
!
!  Setup calculation in Plumed
!
      call plumed_f_gcmd("prepareCalc"//char(0),0)
      call plumed_f_gcmd("prepareDependencies"//char(0),0)
      call plumed_f_gcmd("shareData"//char(0),0)
!
!  Compute contribution from Plumed
!
      call plumed_f_gcmd("performCalc"//char(0),0)
!
!  Add derivatives (negative forces) from Plumed to those from GULP
!
      do i = 1,numat
        xdrv(i) = xdrv(i) - pforcex(i)
        ydrv(i) = ydrv(i) - pforcey(i)
        zdrv(i) = zdrv(i) - pforcez(i)
      enddo
      if (lstr) then
!
!  Add virial from Plumed to those from GULP
!
        do i = 1,3
          do j = 1,3
            virial_m(j,i) = virial_m(j,i) + pvirial(j,i)
          enddo
        enddo
        select case(ndim)
          case(1)
            virial = pvirial(1,1)
          case(2)
            virial = pvirial(1,1) + pvirial(2,2)
          case(3)
            virial = pvirial(1,1) + pvirial(2,2) + pvirial(3,3)
        end select
      endif
!
      deallocate(pforcez)
      deallocate(pforcey)
      deallocate(pforcex)
      t2 = g_cpu_time()
      tplumed = tplumed + t2 - t1
    endif
  endif
#endif

  return
  end subroutine plumed_energy_force

  end module m_plumed
