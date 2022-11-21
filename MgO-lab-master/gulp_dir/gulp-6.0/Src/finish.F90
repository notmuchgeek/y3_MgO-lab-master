  subroutine finish
!
!  Performs tidying up and general output tasks at the end of a run.
!
!   6/03 XML modifications added
!  11/04 Intent added
!  11/06 Unit 4 scratch file now deleted at finish
!   3/07 Chemshell modifications added
!   2/09 Old xml calls removed
!   9/15 Closing of defect channels removed as they are no
!        longer used
!   9/16 lclose dummy argument removed
!   4/17 ChemShell interaction modified
!   6/17 Old ChemShell restored as a compile option
!   1/18 Trace added
!   8/18 Modified for version 2 of OpenKIM
!   9/18 Modified to allow for multiple models in OpenKIM
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
!  Julian Gale, CIC, Curtin University, September 2018
!
  use, intrinsic :: iso_c_binding
  use control
  use general,     only : nwarn
  use gulpchemsh
  use iochannels
#ifdef KIM
  use kim_models
#endif
  use parallel
  use reallocate
#ifdef TRACE
  use trace,       only : trace_in, trace_out
#endif
  implicit none
!
!  Local variables
!
#ifdef KIM
  integer(i4)           :: nm
#endif
#ifdef TRACE
  call trace_in('finish')
#endif
#ifdef KIM
!
!  OpenKIM finishing tasks
!
  if (lkim_model) then
!
!  Finalise KIM model
!  
    do nm = 1,nkimmodel
      call kim_model_destroy(kim_model(nm))
    enddo
  endif
#endif
!
  if (ioproc) then
!*************************
!  Remove scratch files  *
!*************************
    close(4,status='delete')
!***********************
!  Output peak memory  *
!***********************
    call printmemory
!***************************
!  Output timing analysis  *
!***************************
    call outtime
!*****************
!  Output files  *
!*****************
    call outfile
    if (nwarn.gt.1) then
      write(ioout,'(/,''  **** GULP has completed with '',i3,'' warnings - beware! ****'',/)') nwarn
    elseif (nwarn.eq.1) then
      write(ioout,'(/,''  **** GULP has completed with 1 warning - beware! ****'',/)')
    endif
  endif
  call datetime(3_i4)
!*********************************************************
!  Close down message-passing (if appropriate) and stop  *
!*********************************************************
!
! Chemshell - no mpfinish call
!
#ifdef TRACE
  call trace_out('finish')
#endif
#ifdef OLDCS
  if (ichemsh_qm < 0) call mpfinish
#else
  if (ichemsh_link.eq.0) call mpfinish
#endif
!
  return
  end
