  subroutine getvibmoded(nobsmodeptr,maxd2,eigr,mcvloc,mcvptr,nfreqptr,overlap_max)
!
!  Find the vibrational mode with maximal overlap with a 
!  given set of eigenvectors stored in fobsmode
!
!  Distributed memory parallel version
!
!   7/17 Created from getvibmode
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
  integer(i4),                     intent(in)  :: nobsmodeptr    ! Pointer to mode in fobsmode
  integer(i4),                     intent(in)  :: maxd2          ! Left-hand dimension of eigr array
  integer(i4),                     intent(in)  :: mcvloc         ! Local number of modes
  integer(i4),                     intent(in)  :: mcvptr(mcvloc) ! Pointer from local to global modes
  integer(i4),                     intent(out) :: nfreqptr       ! Pointer to frequency with maximal overlap
  real(dp),                        intent(in)  :: eigr(maxd2,*)  ! Eigenvector array
  real(dp),                        intent(out) :: overlap_max    ! Maximum overlap for mode
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ind
  integer(i4)                                  :: m
  integer(i4)                                  :: mode
  integer(i4)                                  :: status
  real(dp),    dimension(:), allocatable       :: otmp
  real(dp),    dimension(:), allocatable       :: overlap
#ifdef TRACE
  call trace_in('getvibmoded')
#endif
!
!  Check that number of atoms matches
!
  if (nobsmodeat(ncf).ne.ncfoc) then
    call outerror('mismatch between number of atoms for mode in fit',0_i4)
    call stopnow('getvibmode')
  endif
!
!  Allocate workspace
!
  allocate(overlap(3*ncfoc),stat=status)
  if (status/=0) call outofmemory('getvibmoded','overlap')
  allocate(otmp(3*ncfoc),stat=status)
  if (status/=0) call outofmemory('getvibmoded','otmp')
!
  overlap(1:3*ncfoc) = 0.0_dp
!
!  Compute overlaps with modes
!
  do m = 1,mcvloc
    mode = mcvptr(m)
    ind = 0
    do i = 1,ncfoc
      overlap(mode) = overlap(mode) + fobsmode(1,i,nobsmodeptr)*eigr(ind+1,m) + &
                                      fobsmode(2,i,nobsmodeptr)*eigr(ind+2,m) + &
                                      fobsmode(3,i,nobsmodeptr)*eigr(ind+3,m)
      ind = ind + 3
    enddo
!
!  Take the absolute of the overlap since we don't care about the sign
!
    overlap(mode) = abs(overlap(mode))
  enddo
!
!  Global sum over overlaps
!
  call sumall(overlap,otmp,3_i4*ncfoc,"getvibmoded","overlap")
!
!  Loop over modes to find best overlap
!
  nfreqptr = 0
  overlap_max = 0.0_dp
  do mode = 1,3*ncfoc
    if (otmp(mode).gt.overlap_max) then
      overlap_max = otmp(mode)
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
!
!  Free workspace
!
  deallocate(otmp,stat=status)
  if (status/=0) call deallocate_error('getvibmoded','otmp')
  deallocate(overlap,stat=status)
  if (status/=0) call deallocate_error('getvibmoded','overlap')
#ifdef TRACE
  call trace_out('getvibmoded')
#endif
!
  return
  end
