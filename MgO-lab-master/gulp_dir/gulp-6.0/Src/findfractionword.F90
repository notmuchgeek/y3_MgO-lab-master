  subroutine findfractionword(fnumber,lfraction,fractionword)
!
!  Tests to see if a number equates to a fraction and if so returns the string
!  Only searches up to 9ths
!
!   5/15 Created
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
  use datatypes
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  character(len=3),                    intent(out)  :: fractionword ! String with fraction
  logical,                             intent(out)  :: lfraction    ! If true then fraction has been found
  real(dp),                            intent(in)   :: fnumber      ! Floating point number between 0 and 1
!
!  Local variables
!
  integer(i4)                                       :: i
  integer(i4)                                       :: j
  real(dp)                                          :: diff
  real(dp)                                          :: thresh
#ifdef TRACE
  call trace_in('findfractionword')
#endif
!
!  Initialise return variables
!
  lfraction = .false.
  fractionword = ' '
!
  thresh = 1.0d-8
!
!  Loop over fractions
!
  do i = 2,9
    do j = 1,i-1
      diff = fnumber - dble(j)/dble(i)
      if (abs(diff).lt.thresh) then
        lfraction = .true.
        write(fractionword(1:1),'(i1)') j
        fractionword(2:2) = '/'
        write(fractionword(3:3),'(i1)') i
      endif
    enddo
  enddo
#ifdef TRACE
  call trace_out('findfractionword')
#endif
  return
  end
