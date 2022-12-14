  subroutine inttobin(i,ibinary,imax)
!
!  Convert an integer number to a binary integer representation
!
!  i       = integer number
!  ibinary = integer array of zeros and ones giving binary value of i
!  imax    = integer length of binary representation (maximum power of two)
!
!   1/18 Trace added
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
!  Julian Gale, CIC, Curtin University, January 2018
!  Scott Woodley, October 1998.
!
  use parallel
#ifdef TRACE
  use trace,      only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4) :: i
  integer(i4) :: ibinary(*)
  integer(i4) :: imax
!
!  Local variables
!
  integer(i4) :: idiv
  integer(i4) :: iloc
  integer(i4) :: itwo
  integer(i4) :: j
#ifdef TRACE
  call trace_in('inttobin')
#endif
!
  iloc = i
  itwo = 2**imax
  if (i.gt.itwo.or.i.lt.0) then
    call outerror('Internal error in binary - 2',0_i4)
    call stopnow('binary')
  else
    do j = imax,1,-1
      itwo = itwo/2
      idiv = iloc/itwo
      iloc = iloc - idiv*itwo
      ibinary(j) = idiv
    enddo
  endif
#ifdef TRACE
  call trace_out('inttobin')
#endif
!
  return
  end

  subroutine bintoint(i,ibinary,imax)
!
!  Convert a binary integer representation to an integer number
!
!  i       = integer number
!  ibinary = integer array of zeros and ones giving binary value of i
!  imax    = integer length of binary representation (maximum power
!            of two)
!
!  Julian Gale, CIC, Curtin University, July 2005.
!
  use datatypes
#ifdef TRACE
  use trace,      only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4) :: i
  integer(i4) :: ibinary(*)
  integer(i4) :: imax
!
!  Local variables
!
  integer(i4) :: itwo
  integer(i4) :: j
#ifdef TRACE
  call trace_in('bintoint')
#endif
!
  itwo = 1
  i = 0
  do j = 1,imax
    i = i + ibinary(j)*itwo
    itwo = itwo*2
  enddo
#ifdef TRACE
  call trace_out('bintoint')
#endif
!
  return
  end
