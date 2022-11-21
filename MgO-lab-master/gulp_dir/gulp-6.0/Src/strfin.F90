  subroutine strfin(lgrad2)
!
!  Completes strain derivatives - should only be called if lgrad1 is true
!
!   5/95 Modifications added to handle symmetrised second derivatives
!   6/95 Symmetry adaption of strain moved here from energy
!   6/97 Corrections for strain-strain second derivatives at non-zero
!        strain added.
!   8/98 Symmetrisation of strains moved into subroutine shared with
!        fefunct
!   4/02 Referencing of derv3 corrected to allow for frozen atoms
!   5/02 Printing for debugging moved to separate routine
!  10/04 Intent added
!  10/04 oldel option removed as no longer used
!  11/07 Modifications for reaxFF Coulomb term added
!   1/08 lreaxFFqreal removed
!   3/09 Non-radial contributions now subtracted from
!        derv3 corrections
!   1/10 Addition of rstrd to strderv corrected for case where lwolf = true
!   3/11 Copying of strderv to stresses added
!   9/13 Modified for parallel execution
!  12/13 Stresses now fully computed here rather than just copying strderv
!   7/17 Stresses computed only for 3D since otherwise the volume is undefined
!   2/18 Trace added
!   9/18 Handling of second derivatives changed. Contributions are added 
!        directly in the computing subroutine to allow for finite strain.
!   9/18 Correction to sderv2 added for non lstraincell case to maintain
!        consistency with old values
!  11/18 Use of reference cell for volume in stress calculation added
!  11/18 Non-radial arrays removed since they are no longer needed
!  11/18 Finite strain flag introduced instead of lstraincell
!   5/19 Real space second derivatives should now be complete and not need addition here
!   5/19 Reciprocal space second derivatives should now be complete and not need addition here
!   7/20 Modifications for gfortran v10
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
  use configurations,   only : lanisotropicpresscfg, rvcfg
  use g_constants,      only : evtoj
  use control
  use current
  use derivatives
  use iochannels
  use numbers
  use parallel
  use symmetry
#ifdef TRACE
  use trace,            only : trace_in, trace_out
#endif
  implicit none
!
!  Passed arguments
!
  logical, intent(in)     :: lgrad2
!
!  Local variables
!
  integer(i4)             :: i
  integer(i4)             :: j
  real(dp)                :: rvol
  real(dp)                :: rvtmp(3,3)
  real(dp)                :: vol
  real(dp)                :: volume
!
  if (.not.lstr) return
#ifdef TRACE
  call trace_in('strfin')
#endif
!**************************************************
!  Symmetrisation of strain vectors and matrices  *
!**************************************************
  if (lsymderv) then
    call strsym
  endif
!****************************************************************
!  Add atomistic first strain derivatives to total derivatives  *
!****************************************************************
  do i = 1,nstrains
    strderv(i) = strderv(i) + rstrd(i)
    stresses(i) = strderv(i)
  enddo
  if (ndim.eq.3) then
!*******************************************************************
!  Convert strain derivatives to stresses and adjust units to GPa  *
!*******************************************************************
    if (lfinitestrain) then
!
!  Use reference cell for this algorithm
!
      rvtmp(1:3,1:3) = rvcfg(1:3,1:3,ncf)
      vol = volume(rvtmp)
    else
      vol = volume(rv)
    endif
    rvol = 10.0_dp*evtoj*1.0d20/vol
    stresses(1:6) = rvol*stresses(1:6)
!
!  Correct stresses for any applied pressure
!
    stresses(1:3) = stresses(1:3) - press
    if (lanisotropicpresscfg(ncf)) then
      stresses(1:6) = stresses(1:6) - anisotropicpress(1:6)
    endif
  endif
  if (lgrad2) then
!*********************************
!  Symmetrise elastic constants  *
!*********************************
    do i = 1,nstrains
      do j = 1,i-1
        sderv2(j,i) = sderv2(i,j)
      enddo
    enddo
!
!  Equivalence sdrv2 and sderv2
!
    do i = 1,nstrains
      do j = 1,nstrains
        sdrv2(j,i) = sderv2(j,i)
      enddo
    enddo
  endif
#ifdef TRACE
  call trace_out('strfin')
#endif
!
  return
  end
