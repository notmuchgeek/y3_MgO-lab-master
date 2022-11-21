  subroutine gulp_setup
!
!  Performs setup tasks for GULP
!
!   9/06 Created from gulp.F
!   3/07 initcomms renamed to GULP_initcomms
!   8/07 Call to GULP_initcomms removed to avoid MPI error
!   1/08 Declaration of ierror wrapped with ifdef
!   2/09 Modified to accommodate new version of FoX and gulp_cml
!   6/09 Call to banner changed to gulp_banner
!   6/09 Renamed to gulp_setup for consistency with other names
!   6/09 Output of extra species info added for PDF/CML
!   9/09 Number for buffer channel now accessed from module
!   7/11 Version incremented to 4.0
!   8/11 Output of hostname added to setup info
!  12/12 Duplicate use of iochannels removed
!   9/15 Opening of defect channels 41, 42, 48 removed since these
!        are no longer used
!   7/16 Modified to allow for one off allocation of pkim_model
!   4/17 ChemShell interaction modified
!   6/17 Old ChemShell restored as a compile option
!   1/18 Memory allocations checked after keywords have been set
!        through a second call to initmemory
!   1/18 Trace added
!   5/18 EEM setup added
!   8/18 Version number updated
!   8/18 KIM handling removed for version 2.
!   7/20 Version number updated
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
!  Julian Gale, CIC, Curtin University, August 2020
!
  use configurations
  use control
  use current
  use general
  use gulp_cml,        only : lcml, gulp_cml_init, gulp_cml_outkey
  use gulp_cml_phonon, only : gulp_cml_outspec
  use gulpchemsh
  use iochannels,      only : iotmp, ioout
  use m_pdfneutron,    only : lpdfout, outpdfin
  use parallel
#ifdef TRACE
  use trace,           only : trace_in, trace_out
#endif
#ifdef ACCELRYS
  use license
#endif
  implicit none
!
  character(len=40) :: hostname
  integer(i4)       :: iline
#ifdef ACCELRYS
  integer(i4)       :: ierror
#endif
  integer(i4)       :: hlength
  integer(i4)       :: status
  logical           :: lopened
  ! ChemShell additions
  logical           :: opc_root
  external opc_root
#ifndef OLDCS
  integer           :: ipunch
  integer           :: i
#endif
#ifdef MPI
  include 'mpif.h'
#endif
#ifdef TRACE
  call trace_in('gulp_setup')
#endif

!*************************
!  Nullify all pointers  *
!*************************
  call nullpointer
!*************************
!  Initialise memory     *
!*************************
  call initmemory
!*************************
!  Initialise variables  *
!*************************
  iline = 0
  call initial
!******************
!  Output header  *
!******************
  version = 6.0_dp
#ifdef ACCELRYS
  call license_checkout(nprocs,ierror)
  if (ierror.ne.0) call gulpfinish
! Set traps for signals
  
  call setup_traps()
#endif
  call gulp_banner
!**************************
!  Get local information  *
!**************************
  call local
!****************************
!  Set element information  *
!****************************
  call setele(lopened)
!*******************************************************
!  Pre-processor passes before reading input properly  *
!*******************************************************
  call channels(iotmp,.true.)
  call firstpass
  rewind(iotmp)
  call secondpass
  rewind(iotmp)
!*********************
!  Read in keywords  *
!*********************
  call getkeyword(iline)
  call setkeyword
!***********************
!  Process main input  *
!***********************
  call inword(iline)
!*******************************************************
!  Set keywords again in case any were in the library  *
!*******************************************************
  call setkeyword
  close(iotmp,status='delete')
!****************************************
!  Set charge equilibration parameters  *
!****************************************
  call seteem
!***************************
!  Output keyword details  *
!***************************
  call outkey
!*********************************************************
!  Re-initialise memory as keywords may change settings  *
!*********************************************************
  call initmemory
!**************
!  Site name  *
!**************
  if (ioproc) then
    if (site.ne.' ') then
      write(ioout,'(''* '',a76,'' *'')') site(1:76)
      write(ioout,'(''********************************************************************************'')')
    endif
    call datetime(1_i4)
    write(ioout,'(''  Number of CPUs = '',i5,/)') nprocs
    hostname = ' '
    call get_environment_variable('HOSTNAME',hostname,hlength,status)
    if (status.eq.0) then
      write(ioout,'(''  Host name      = '',a40,/)') hostname
    elseif (status.eq.1) then
      call get_environment_variable('HOST',hostname,hlength,status)
      if (status.eq.0) then
        write(ioout,'(''  Host name      = '',a40,/)') hostname
      endif
    endif
  endif
!***********************
!  CML initialisation  *
!***********************
  if (lcml) then
    call gulp_cml_init
    call gulp_cml_outkey
  endif
!**************************
!  One off initial set up *
!**************************
  call setcfg
!*******************
!  Species output  *
!*******************
  if (ioproc) then
    call outspec
    if (lcml) call gulp_cml_outspec
    if (lpdfout) call outspec_pdf
  endif
!*****************************
!  Electronegativity output  *
!*****************************
  if (ioproc) then
    call outeem
  endif
!***********************************************************************
!  Check PDF related settings, convert frequencies and output PDF info *
!***********************************************************************
  if (lpdfout) call outpdfin
!************************
!  Output general info  *
!************************
  if (ioproc) call outgen
!*******************************
!  Output polarisability info  *
!*******************************
  if (ioproc) call outpolar
!*********************
!  Check potentials  *
!*********************
  call checkpot
!**********************
!  Output potentials  *
!**********************
  call outpot
!*****************************
!  Check options for method  *
!*****************************
  call methodok
!*****************************
!  ChemShell charge-only run *
!*****************************
#ifdef OLDCS
  if (ichemsh_qm .eq. 99 .and. opc_root()) then
    call ExportGulpCharges(numat, qlcfg)
  endif
#else
  if (ichemsh_output .eq. 1 .and. opc_root()) then

    ipunch = 7

    open(ipunch,file='gulp.charges',form='formatted')

    write(ipunch,*) "block=matrix records=0"
    write(ipunch,*) "block=matrix_title records=1"
    write(ipunch,*) "Gulp charges"

    write(ipunch,101) numat, 1, numat

    do i = 1,numat
        write(ipunch,100) qlcfg(i)
    enddo

    close(unit=ipunch)

  endif
100  format(f28.14)
101  format('block=dense_real_matrix records=',i6,' dimensions=',2i6)
#endif
#ifdef TRACE
  call trace_out('gulp_setup')
#endif
!
  return
  end
