  subroutine ind1toijk(indin,ii,jj,kk)
!
!  Converts linear index to x,y and z components for
!  unit cell in 111 case.
!
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
  integer(i4), intent(in)  :: indin
  integer(i4), intent(out) :: ii
  integer(i4), intent(out) :: jj
  integer(i4), intent(out) :: kk
!
!  Local variables
!
  integer(i4)              :: ind
#ifdef TRACE
  call trace_in('ind1toijk')
#endif
!
  ind = indin - 1
  ii = ind/9
  ind = ind - 9*ii
  jj = ind/3
  kk = ind - 3*jj + 1
  jj = jj + 1
  ii = ii + 1
#ifdef TRACE
  call trace_out('ind1toijk')
#endif
!
  return
  end
!
  subroutine ind2toijk(indin,ii,jj,kk)
!
!  Converts linear index to x,y and z components for
!  unit cell in 222 case.
!
!   4/08 Created from ind1toijk
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
  integer(i4), intent(in)  :: indin
  integer(i4), intent(out) :: ii
  integer(i4), intent(out) :: jj
  integer(i4), intent(out) :: kk
!
!  Local variables
!
  integer(i4)              :: ind
#ifdef TRACE
  call trace_in('ind2toijk')
#endif
!
  ind = indin - 1
  ii = ind/25
  ind = ind - 25*ii
  jj = ind/5
  kk = ind - 5*jj + 1
  jj = jj + 1
  ii = ii + 1
#ifdef TRACE
  call trace_out('ind2toijk')
#endif
!
  return
  end
