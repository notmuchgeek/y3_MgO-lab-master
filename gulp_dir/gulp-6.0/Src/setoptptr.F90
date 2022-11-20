  subroutine setoptptr(lfreezeok)
!
!  Sets up pointer mapping atoms to be optimised
!
!   1/17 Created 
!   7/17 lfreezeok now passed in
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
!  Copyright Curtin University 2017
!
!  Julian Gale, CIC, Curtin University, July 2017
!
  use current
  use optimisation
  use parallel
  implicit none
!
!  Passed variables
!
  logical,                  intent(in)              :: lfreezeok
!
!  Local variables
!
  integer(i4)                                       :: i
  logical                                           :: lopi
!
!  Initialise values
!
  noptat = 0
  noptatloc = 0
  noptatrptr(1:numat) = 0
  noptatlocrptr(1:numat) = 0
!
!  Loop over atoms 
!
  do i = 1,numat
    lopi = (.not.lfreezeok.or.lopf(i))
    if (lopi) then
      noptat = noptat + 1
      noptatptr(noptat) = i
      noptatrptr(i) = noptat
      if (atom2local(i).gt.0) then
        noptatloc = noptatloc + 1
        noptatlocptr(noptatloc) = i
        noptatlocrptr(i) = noptatloc
      endif
    endif
  enddo
!
  return
  end
