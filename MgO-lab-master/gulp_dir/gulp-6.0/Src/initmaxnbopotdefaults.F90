  subroutine initmaxnbopotdefaults(i)
!
!  Initialises the arrays associated with maxnbopot
!
!   9/10 Created from changemax routine
!  12/15 nBOtapertype added
!  11/20 Tersoff reorganised
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
!  Julian Gale, CIC, Curtin University, November 2020
!
  use bondorderdata
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in) :: i
!
!  Initialise defaults for new part of array
!
  if (i.ge.1.and.i.le.maxnbopot) then
    BOacoeff(i) = 0.0_dp
    BObcoeff(i) = 0.0_dp
    BOchiA(i) = 1.0_dp
    BOchiR(i) = 1.0_dp
    BOzacoeff(i) = 0.0_dp
    BOzbcoeff(i) = 0.0_dp
    BOcombi(i) = .false.
    rBOmax(i) = 0.0_dp
    rBOmin(i) = 0.0_dp
    nBOtapertype(i) = 1
    nBOtyp1(i) = 0
    nBOtyp2(i) = 0
    nBOtypeT(i) = 1
  endif
!
  return
  end
