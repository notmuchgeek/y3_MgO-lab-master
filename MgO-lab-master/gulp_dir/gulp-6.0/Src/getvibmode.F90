  subroutine getvibmode(nobsmodeptr,maxd2,eigr,nfreqptr,overlap_max)
!
!  Find the vibrational mode with maximal overlap with a 
!  given set of eigenvectors stored in fobsmode
!
!   1/12 Created
!   7/17 Atom number check moved outside loop over modes
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
  use observables,    only : fobsmode, nobsmodeat
  use partial,        only : ncfoc
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)  :: nobsmodeptr    ! Pointer to mode in fobsmode
  integer(i4), intent(in)  :: maxd2          ! Left-hand dimension of eigr array
  integer(i4), intent(out) :: nfreqptr       ! Pointer to frequency with maximal overlap
  real(dp),    intent(in)  :: eigr(maxd2,*)  ! Eigenvector array
  real(dp),    intent(out) :: overlap_max    ! Maximum overlap for mode
!
!  Local variables
!
  integer(i4)              :: i
  integer(i4)              :: ind
  integer(i4)              :: mode
  real(dp)                 :: overlap
#ifdef TRACE
  call trace_in('getvibmode')
#endif
!
!  Check that number of atoms matches
!
  if (nobsmodeat(ncf).ne.ncfoc) then
    call outerror('mismatch between number of atoms for mode in fit',0_i4)
    call stopnow('getvibmode')
  endif
!
!  Loop over modes to find best overlap
!
  nfreqptr = 0
  overlap_max = 0.0_dp
  do mode = 1,3*ncfoc
    overlap = 0.0_dp
    ind = 0
    do i = 1,ncfoc
      overlap = overlap + fobsmode(1,i,nobsmodeptr)*eigr(ind+1,mode) + &
                          fobsmode(2,i,nobsmodeptr)*eigr(ind+2,mode) + &
                          fobsmode(3,i,nobsmodeptr)*eigr(ind+3,mode)
      ind = ind + 3
    enddo
!
!  Take the absolute of the overlap since we don't care about the sign
!
    overlap = abs(overlap)
    if (overlap.gt.overlap_max) then
      overlap_max = overlap
      nfreqptr = mode
    endif
  enddo
!
!  Check that mode was found
!
  if (nfreqptr.eq.0) then
    call outerror('no mode with greater than zero overlap found',0_i4)
    call stopnow('getvibmode')
  endif
#ifdef TRACE
  call trace_out('getvibmode')
#endif
!
  return
  end
