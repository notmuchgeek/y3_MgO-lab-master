  subroutine changemaxnpts
!
!  Alters the size of the arrays associated with maxnpts
!  for one-dimensional arrays
!
!  12/08 Migrated to version 3.5 and converted to f90 format
!   9/10 Initialisations now performed in a subroutine
!   4/17 Parallel arrays added
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
!  Julian Gale, CIC, Curtin University, April 2017
!
  use cosmic
  use parallel,      only : npts2local, npts2node, node2pts
  use reallocate
  implicit none
!
  integer(i4)       :: ierror, i
  integer(i4), save :: oldmaxnpts = 0
!
  call realloc(sphere1,3_i4,maxnpts,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnpts','sphere1')
  call realloc(cosmoatomptr,maxnpts,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnpts','cosmoatomptr')
  call realloc(cosmoBq,maxnpts,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnpts','cosmoBq')
  call realloc(sas,3_i4,maxnpts,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnpts','sas')
  call realloc(segweight,maxnpts,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnpts','segweight')
  call realloc(nnearseg,maxnpts,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnpts','nnearseg')
  call realloc(nnearsegptr,maxnearseg,maxnpts,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnpts','nnearsegptr')
  call realloc(nnearsegptrcell,maxnearseg,maxnpts,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnpts','nnearsegptrcell')
  call realloc(npwt,maxnpts,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnpts','npwt')
  call realloc(npwtptr,maxnpwt,maxnpts,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnpts','npwtptr')
!
  call realloc(node2pts,maxnpts,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnpts','node2pts')
  call realloc(npts2local,maxnpts,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnpts','npts2local')
  call realloc(npts2node,maxnpts,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnpts','npts2node')
!     
!  Initialise new parts of data arrays
!     
  if (maxnpts.gt.oldmaxnpts) then
    do i = oldmaxnpts+1,maxnpts
      call initmaxnptsdefaults(i)
    enddo
  endif
!
!  Save current value of maxnpts for next call
!
  oldmaxnpts = maxnpts
!
  return
  end
!
  subroutine changemaxcosmoA
!
!  Alters the size of the arrays associated with maxnpts2
!  for two-dimensional arrays
!
!   1/05 cosmoB removed
!   4/16 Trap for integer precision being exceed added
!   4/17 Change for possible 2D structure of cosmoA made
!   4/17 Name changed from changemaxnpts2 to changemaxcosmoA
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
!  Julian Gale, Curtin University, April 2017
!
  use cosmic
  use reallocate
  implicit none
!
  integer(i4)       :: ierror
!
  call realloc(cosmoA,maxcosmoA,maxcosmoAu,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcosmoA','cosmoA')
!
  return
  end
