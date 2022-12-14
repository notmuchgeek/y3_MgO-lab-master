  subroutine initmaxeamspecdefaults(i)
!
!  Initialises the arrays associated with maxeamspec
!
!   9/10 Created from changemax routine
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
!  Copyright Curtin University 2014
!
!  Julian Gale, CIC, Curtin University, August 2014
!
  use eam
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in) :: i
!
  integer(i4)             :: j
!
!  Initialise defaults for new part of array
!
  if (i.ge.1.and.i.le.maxeamspec) then
    symboleamspec(1:2,i) = ' '
    eamdenfile(i) = ' '
    ndenfncomp(i) = 0
    do j = 1,maxeamden
      ndenfn(j,i) = 0
      neammeamorder(j,i) = 1
      denpar(1:16,1:maxmeamorder,j,i) = 0.0_dp
    enddo
    eamalloy(1,i) = 1.0_dp
    eamalloy(2,i) = 0.0_dp
    eamtaperdrho(i) = 0.0_dp
    eamtaperrho(i) = 0.0_dp
    lmeamspec(i) = .false.
  endif
!
  return
  end
