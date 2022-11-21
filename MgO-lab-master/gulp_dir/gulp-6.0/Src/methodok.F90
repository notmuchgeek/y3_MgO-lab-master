  subroutine methodok
!
!  Check whether requested options are compatible with the
!  derivatives that are available for the present potential
!  model.
!
!   4/01 Created 
!   9/04 Check to make sure that EEM / nboQ are not used in
!        the same run
!   7/05 Streitz and Mintmire modifications added
!   5/06 Check for switch to second derivative optimisation
!        method in parallel and disable
!   7/06 Check for whether Mott-Littleton method is OK added
!   4/08 Phonon/frequency calculations allowed to proceed for
!        case when lnoanald2 is true if finite differences are
!        specified
!   5/08 Finite differences enabled when eigen keyword is specified
!   5/09 Free energy minimisation disabled for EAM/MEAM
!   6/09 Modified so that properties can be calculated with finite
!        differences where analytic second derivatives are unavailable
!   3/10 Incorrect flagging of error for lsuttonc case removed
!   6/10 Check for incompatibility between solvation model and other
!        options added
!   9/12 Pacha added
!   9/13 Checks on third derivatives added for Raman option
!   2/14 Trap for variable charges in parallel removed
!   8/16 Check that FFTW3 is available for SPME
!   4/17 Second derivatives enabled in parallel for rfo and frequencies
!   5/17 Free energy calculations enabled in parallel
!   5/17 Trap on parallel with minimisation change removed
!   7/17 Trap on charge second derivatives added
!   7/17 Trap on breathing shell second derivatives added
!   7/17 Trap on partial occupancy second derivatives added
!   7/17 Trap on COSMO/COSMIC second derivatives added
!   7/17 Trap on makeeigenarrays keyword for PDF in parallel added
!   7/17 Trap on use of force constant supercell with variable charges
!  12/17 Trap for models that are not currently enabled for fastfd added
!   1/18 Intensities now allowed for finite difference case
!   1/18 Trap on ghostcell phonons in parallel added
!   1/18 Trap on grueneisen for models without analytic third derivatives added
!   1/18 Trace added
!   6/18 Restriction on charge derivatives in parallel for EEM removed as long
!        as iterative algorithm is not being used.
!   7/18 Restriction reintroduced on charge second derivatives in parallel 
!   5/19 Finite differences enabled for RFO
!   5/19 Finite difference flag split for first and second derivatives
!   5/20 Check for compatibility with rigid molecules
!   6/20 Error message updated to mention rigid as a possible issue
!   7/20 Further checks added for compatability with rigid molecules
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
  use bondorderdata,  only : nboQ, nboQ0, nbopot
  use cellmultipole,  only : icmm
  use configurations, only : lbsmat, nasum, occucfg
  use control
  use cosmic,         only : lcosmic
  use derivatives,    only : lfcsupercell
  use element,        only : lqeq, lSandM, lpacha
  use general,        only : phondiff, lfinitediff2
  use iochannels
  use kim_models,     only : lkim_model
  use m_fft3d,        only : lfftw3
  use m_pdfneutron,   only : lmakeeigarray
  use mdlogic,        only : lmd
  use optimisation,   only : lminch, mintype
  use parallel
  use polarise,       only : lpolar
  use spme,           only : lspme
  use sutton,         only : lsuttonc
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Local variables
!
  logical            :: lerror
  logical            :: lneed2nd
  real(dp)           :: totocc
#ifdef TRACE
  call trace_in('methodok')
#endif
!
  lerror = .false. 
!
  if (lnoanald2) then
!**********************************************
!  Analytical second derivatives unavailable  *
!**********************************************
    if (ldefect) lerror = .true.
    if (leigen.and.phondiff.lt.1.0d-12) lerror = .true.
    if (lfreq.and.phondiff.lt.1.0d-12) lerror = .true.
    if (lphon.and.phondiff.lt.1.0d-12) lerror = .true.
    if (linten.and.phondiff.lt.1.0d-12) lerror = .true.
    if (lrfo.and..not.lfinitediff2) lerror = .true.
    if (lerror) then
      if (ioproc) then
        call outerror('analytical second derivatives unavailable',0_i4)
        write(ioout,'('' The current potential model is incompatible with keywords:'',/)')
        write(ioout,'('' defect eigenvectors free frequency intensity phonon property rfo'',/)')
        write(ioout,'('' Try using finite differences'',/)')
      endif
      call stopnow('methodok')
    endif
  endif
  if (lnoanald3) then
!*********************************************
!  Analytical third derivatives unavailable  *
!*********************************************
    if (lfree) lerror = .true.
    if (lraman) lerror = .true.
    if (lgrueneisen) lerror = .true.
    if (lerror) then
      if (ioproc) then
        call outerror('analytical third derivatives unavailable',0_i4)
        if (lrigid) then
          write(ioout,'('' The use of rigid and/or current potential model is incompatible with keywords:'',/)')
          write(ioout,'('' free raman grueneisen '',/)')
        else
          write(ioout,'('' The current potential model is incompatible with keywords:'',/)')
          write(ioout,'('' free raman grueneisen '',/)')
        endif
      endif
      call stopnow('methodok')
    endif
  endif
  if (lnomottlittleton.and.ldefect) then  
!**************************************
!  Mott-Littleton method unavailable  *
!**************************************
    if (ioproc) then
      call outerror('Mott-Littleton method unavailable',0_i4)
      write(ioout,'('' The current potential model is incompatible with keywords:'',/)')
      write(ioout,'('' defect'',/)')
    endif
    call stopnow('methodok')
  endif
!*************************
!  Parallel unavailable  *
!*************************
  if (nprocs.gt.1) then
    lneed2nd = .false.
!
!  Do we need second derivatives?
!
    if (lopt.and.(.not.lconj.and..not.llbfgs)) lneed2nd = .true.
    if (leigen) lneed2nd = .true.
    if (lfreq) lneed2nd = .true.
    if (linten) lneed2nd = .true.
    if (lphon) lneed2nd = .true.
    if (lprop) lneed2nd = .true.
    if (lrfo) lneed2nd = .true.
    if (lfree) lneed2nd = .true.
!
!  Ensure that any minimizer change won't require second derivatives
!
    if (lminch.and.mintype.ne.5) lneed2nd = .true.
!
!  Variable charges
!
    if (lneed2nd.and.lDoQDeriv2) then
      if (ioproc) then
        call outerror('variable charge second derivatives unavailable in parallel',0_i4)
      endif
      call stopnow('methodok')
    endif
!
!  Breathing shells
!
    if (lneed2nd.and.any(lbsmat)) then
      if (ioproc) then
        call outerror('breathing shell second derivatives unavailable in parallel',0_i4)
      endif
      call stopnow('methodok')
    endif
!
!  Partial occupancy
!
    totocc = sum(occucfg(1:nasum))
    if (lneed2nd.and.(abs(totocc-dble(nasum)).gt.1.0d-3)) then
      if (ioproc) then
        call outerror('partial occupancy second derivatives unavailable in parallel',0_i4)
      endif
      call stopnow('methodok')
    endif
!
!  COSMO/COSMIC solvation
!
    if (lneed2nd.and.lcosmic) then
      if (ioproc) then
        call outerror('COSMIC second derivatives unavailable in parallel',0_i4)
      endif
      call stopnow('methodok')
    endif
    if (lneed2nd.and.lcosmo) then
      if (ioproc) then
        call outerror('COSMO second derivatives unavailable in parallel',0_i4)
      endif
      call stopnow('methodok')
    endif
!
!  PDF makeeigenarrays keyword
!
    if (lmakeeigarray) then
      if (ioproc) then
        call outerror('PDF unavailable in parallel with makeEigenArrays keyword',0_i4)
      endif
      call stopnow('methodok')
    endif
  endif
!*****************************************
!  Mutually exclusive potential options  *
!*****************************************
  if (nboQ.gt.0.and.leem) then
!
!  EEM & Bond order charge setting
!
    if (ioproc) then
      call outerror('EEM and bond order charges are mutually exclusive',0_i4)
    endif
    call stopnow('methodok')
  endif
!
  if (nboQ.gt.0.and.icmm.gt.0) then
!
!  Bond order charges and cell multipole method
!
    if (ioproc) then
      call outerror('Bond order charges not available with cell multipoles',0_i4)
    endif
    call stopnow('methodok')
  endif
  if (nboQ.gt.0.and.ldefect) then
!
!  Bond order charges and defect calculations
!
    if (ioproc) then
      call outerror('Bond order charges not available with defect calculations',0_i4)
    endif
    call stopnow('methodok')
  endif
!
!  Variable charges with force constant supercell
!
  if (lfcsupercell.and.lDoQDeriv2) then
    if (ioproc) then
      call outerror('charge second derivatives unavailable in FC supercell',0_i4)
    endif
    call stopnow('methodok')
  endif
!
!  Ghostcell algorithm is only available in serial at present
!
  if (lfcsupercell.and.lghost.and.nprocs.gt.1) then
    if (ioproc) then
      call outerror('ghostcell phonons unavailable in parallel',0_i4)
    endif
    call stopnow('methodok')
  elseif (lfcsupercell.and.lghost.and.lrigid) then
    if (ioproc) then
      call outerror('ghostcell phonons unavailable for rigid molecules',0_i4)
    endif
    call stopnow('methodok')
  endif
!*****************************************************************************
!  Solvation and variable charge models - no SCF contribution yet available  *
!*****************************************************************************
  if (lcosmo) then
    if (leem.or.lSandM.or.lqeq.or.lpacha) then
      if (ioproc) then
        call outerror('Solvation not yet compatible with variable charges',0_i4)
      endif
      call stopnow('methodok')
    elseif (lreaxFF) then
      if (ioproc) then
        call outerror('Solvation not yet compatible with ReaxFF variable charges',0_i4)
      endif
      call stopnow('methodok')
    endif
  endif
!*********************
!  SPME needs FFTW3  *
!*********************
  if (lspme.and..not.lfftw3) then
    if (ioproc) then
      call outerror('SPME cannot be used as GULP was not compiled with FFTW3',0_i4)
    endif
    call stopnow('methodok')
  endif
!********************************************
!  FastFD keyword - models not implemented  *
!********************************************
  if (lfastfd) then
    if (lcosmo) then
      if (ioproc) then
        call outerror('solvation models not yet implemented for fastfd',0_i4)
      endif
      call stopnow('methodok')
    endif
    if (lpolar) then
      if (ioproc) then
        call outerror('polarisation not yet implemented for fastfd',0_i4)
      endif
      call stopnow('methodok')
    endif
    if ((nboQ+nboQ0).gt.0) then
      if (ioproc) then
        call outerror('bondorder charges not yet implemented for fastfd',0_i4)
      endif
      call stopnow('methodok')
    endif
  endif
!********************
!  Rigid molecules  *
!********************
  if (lrigid) then
    if (lmd) then
      if (ioproc) then
        call outerror('molecular dynamics not yet implemented for rigid molecules',0_i4)
      endif
      call stopnow('methodok')
    endif
    if (ldefect) then
      if (ioproc) then
        call outerror('defect calculations not yet implemented for rigid molecules',0_i4)
      endif
      call stopnow('methodok')
    endif
    if (leem.or.lSandM.or.lqeq.or.lpacha) then
      if (ioproc) then
        call outerror('variable charges not yet implemented for rigid molecules',0_i4)
      endif
      call stopnow('methodok')
    endif
    if (lreaxFF) then
      if (ioproc) then
        call outerror('ReaxFF not yet implemented for rigid molecules',0_i4)
      endif
      call stopnow('methodok')
    endif
    if (lreaxFF) then
      if (ioproc) then
        call outerror('ReaxFF not yet implemented for rigid molecules',0_i4)
      endif
      call stopnow('methodok')
    endif
    if (lbrenner) then
      if (ioproc) then
        call outerror('REBO not yet implemented for rigid molecules',0_i4)
      endif
      call stopnow('methodok')
    endif
    if (lEDIP) then
      if (ioproc) then
        call outerror('EDIP not yet implemented for rigid molecules',0_i4)
      endif
      call stopnow('methodok')
    endif
    if (lsuttonc) then
      if (ioproc) then
        call outerror('EAM and MEAM not yet implemented for rigid molecules',0_i4)
      endif
      call stopnow('methodok')
    endif
    if (lkim_model) then
      if (ioproc) then
        call outerror('Open KIM not yet implemented for rigid molecules',0_i4)
      endif
      call stopnow('methodok')
    endif
    if (nbopot.gt.0) then
      if (ioproc) then
        call outerror('Bond order potentials not yet implemented for rigid molecules',0_i4)
      endif
      call stopnow('methodok')
    endif
  endif
#ifdef TRACE
  call trace_out('methodok')
#endif
!
  return
  end
