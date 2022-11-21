  subroutine findspacegroup(ngroup)
!
!  Finds the high symmetry space group for a given structure that is initially in P 1
!  Note that this is a very simplistic trial and error algorithm for now!
!
!   8/13 Created from testspacegroup
!   8/13 Pre-screening for cell type added
!   8/13 Search over origin settings added
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
!  On exit:
!
!  ngroup = Best guess at space group for structure
!
  use current
  use symmetry
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(out)  :: ngroup
!
!  Local variables
!
  integer(i4)               :: ifsosave
  integer(i4)               :: ifhrsave
  integer(i4)               :: ngroupmax
  integer(i4)               :: ngroupmin
  integer(i4)               :: norigin
  integer(i4)               :: nspcgsave
  logical                   :: lvalid
#ifdef TRACE
  call trace_in('findspacegroup')
#endif
!
!  Initialise return flag
!
  lvalid = .false.
!
!  Save input space group
!
  nspcgsave = nspcg(ncf)
  ifhrsave  = ifhr(ncf)
  ifsosave  = ifso(ncf)
!
!  Find cell type
!
  call findcelltype(a,b,c,alpha,beta,gamma,ictype,icfhr)
!
!  Set range of space groups based on cell type
!
  if (ictype.eq.1) then
    ngroupmin = 1
    ngroupmax = 2
  elseif (ictype.eq.2) then
    ngroupmin = 3
    ngroupmax = 15
  elseif (ictype.eq.3) then
    ngroupmin = 16
    ngroupmax = 74
  elseif (ictype.eq.4) then
    ngroupmin = 75
    ngroupmax = 142
  elseif (ictype.eq.5) then
    ngroupmin = 143
    ngroupmax = 194
  elseif (ictype.eq.6) then
    ngroupmin = 195
    ngroupmax = 230
  endif
!
!  Try to find valid space group based on cell type
!
  ngroup = ngroupmax + 1
  do while (ngroup.gt.ngroupmin.and..not.lvalid)
    ngroup = ngroup - 1
!
!  Loop over standard origin settings
!
    norigin = 0
    do while (norigin.lt.2.and..not.lvalid)
      norigin = norigin + 1
!
!  Set space group number for now
!
      nspcg(ncf) = ngroup
      iflags(ncf) = 0
      ifso(ncf) = norigin
      ifhr(ncf) = icfhr
!
!  Set up symmetry operators for this space group
!
      call symmet
!
!  Call subroutine to test space group
!
      call testspacegroup(ngroup,lvalid)
!
!  End of loop over origin settings
!
    enddo
!
!  End of loop over symmetry operators
!
  enddo
!
  if (.not.lvalid) then
!
!  Cell type has failed and so search over all space groups
!
!  Loop over space groups from 230 to 2 looking for valid space groups
!
!  Take the first one that is found since this should be the highest symmetry space group most of the time
!
    lvalid = .false.
    ngroup = 231
    do while (ngroup.gt.1.and..not.lvalid)
      ngroup = ngroup - 1
!
!  Loop over standard origin settings
!
      norigin = 0
      do while (norigin.lt.2.and..not.lvalid)
        norigin = norigin + 1
!
!  Set space group number for now
!
        nspcg(ncf) = ngroup
        iflags(ncf) = 0
        ifso(ncf) = norigin
        ifhr(ncf) = 0
!
!  Set up symmetry operators for this space group
!
        call symmet
!
!  Call subroutine to test space group
!
        call testspacegroup(ngroup,lvalid)
!
!  End of loop over origin settings
!
      enddo
!
!  End of loop over symmetry operators
!
    enddo 
  endif  
!
!  If no space group found then set to P 1
!
  if (.not.lvalid) then
    ngroup = 1
  endif  
!
!  Reset space group
!
  nspcg(ncf) = nspcgsave
  ifhr(ncf)  = ifhrsave
  ifso(ncf)  = ifsosave
#ifdef TRACE
  call trace_out('findspacegroup')
#endif
!
  return
  end 
