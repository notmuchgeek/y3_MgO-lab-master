  subroutine changemaxd2
!
!  Alters the size of the arrays associated with maxd2 or maxd2u
!
!  NB: The array derv2dk, which is complex, only needs to be 
!      allocated if group velocities are to be computed.
!
!  10/14 d2cell added
!  12/19 Rigid molecule arrays added
!   4/20 derv3c added
!   4/20 derv3c changes reversed as they are no longer required
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
!  Julian Gale, CIC, Curtin University, April 2020
!
  use control,     only : lgroupvelocity, lrigid
  use derivatives
  use molecule,    only : maxmol
  use reallocate
  implicit none
!
  integer(i4) :: ierror
!
  call realloc(diagblock,maxd2,3_i4,ierror)
  if (ierror.ne.0) call outofmemory('changemaxd2','diagblock')
  call realloc(derv2,maxd2,maxd2u,ierror)
  if (ierror.ne.0) call outofmemory('changemaxd2','derv2')
  call realloc(derv2d,maxd2u,ierror)
  if (ierror.ne.0) call outofmemory('changemaxd2','derv2d')
  call realloc(derv3,maxd2,6_i4,ierror)
  if (ierror.ne.0) call outofmemory('changemaxd2','derv3')
  call realloc(dervi,maxd2,maxd2u,ierror)
  if (ierror.ne.0) call outofmemory('changemaxd2','dervi')
  if (lgroupvelocity) then
    call realloc(derv2dk,3_i4,maxd2,maxd2u,ierror)
    if (ierror.ne.0) call outofmemory('changemaxd2','derv2dk')
  endif
  if (lfcsupercell) then
    call realloc(d2cell,maxd2,maxd2u,maxd2cells,ierror)
    if (ierror.ne.0) call outofmemory('changemaxd2','d2cell')
  endif
  if (lrigid) then
    call realloc(molQCdrv,maxd2,3_i4,maxmol,ierror)
    if (ierror.ne.0) call outofmemory('changemaxd2','molQCdrv')
    call realloc(molQCdri,maxd2,3_i4,maxmol,ierror)
    if (ierror.ne.0) call outofmemory('changemaxd2','molQCdri')
    call realloc(molTCdrv,maxd2,3_i4,maxmol,ierror)
    if (ierror.ne.0) call outofmemory('changemaxd2','molTCdrv')
    call realloc(molTCdri,maxd2,3_i4,maxmol,ierror)
    if (ierror.ne.0) call outofmemory('changemaxd2','molTCdri')
    if (lgroupvelocity) then
      call realloc(molQCdk,3_i4,maxd2,3_i4,maxmol,ierror)
      if (ierror.ne.0) call outofmemory('changemaxd2','molQCdk')
      call realloc(molTCdk,3_i4,maxd2,3_i4,maxmol,ierror)
      if (ierror.ne.0) call outofmemory('changemaxd2','molTCdk')
    endif
  endif
!
  return
  end
