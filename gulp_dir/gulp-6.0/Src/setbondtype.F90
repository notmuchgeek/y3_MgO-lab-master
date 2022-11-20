  subroutine setbondtype(nati,ntypi,natj,ntypj,bondtype1,bondtype2)
!
!  Determines default bond type for species pair
!
!   3/15 Created
!   8/18 Default bond type set to single as per help text
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
!  Julian Gale, CIC, Curtin University, August 2018
!
  use current
  use element
  use molecule
  implicit none
!
!  Passed variables
!
  integer(i4),                    intent(in)   :: nati
  integer(i4),                    intent(in)   :: natj
  integer(i4),                    intent(in)   :: ntypi
  integer(i4),                    intent(in)   :: ntypj
  integer(i4),                    intent(out)  :: bondtype1
  integer(i4),                    intent(out)  :: bondtype2
!
!  Local variables
!
  integer(i4)                                  :: ii
  integer(i4)                                  :: nti
  integer(i4)                                  :: ntj
  integer(i4)                                  :: ntyi
  integer(i4)                                  :: ntyj
  logical                                      :: lfound
!
!  Initialise bondtype values 
!
  bondtype1 = 1
  bondtype2 = 1
!
!  Check whether there is a default bond type for this pair
!
  if (nbondtype.gt.0) then
    if (nati.eq.natj) then
      nti = nati
      ntj = natj
      if (ntypi.lt.ntypj) then
        ntyi = ntypi
        ntyj = ntypj
      else
        ntyi = ntypj
        ntyj = ntypi
      endif
    elseif (nati.lt.natj) then
      nti = natj
      ntj = nati
      ntyi = ntypi
      ntyj = ntypj
    else
      nti = nati
      ntj = natj
      ntyi = ntypj
      ntyj = ntypi
    endif
    lfound = .false.
    ii = 0
    do while (.not.lfound.and.(ii.lt.nbondtype))
      ii = ii + 1
      if (natbondtype(1,ii).eq.nti.and.natbondtype(2,ii).eq.ntj) then
        lfound = (ntypbondtype(1,ii).eq.ntyi.and.ntypbondtype(2,ii).eq.ntyj)
      endif
    enddo
    if (lfound) then
      bondtype1 = nbondtypeptr(1,ii)
      bondtype2 = nbondtypeptr(2,ii)
    endif
  endif
!
  return
  end
