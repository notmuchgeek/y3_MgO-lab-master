  subroutine changemaxmol
!
!  Alters the size of the arrays associated with maxmol
!
!   7/07 GCMC molecule flag added
!   6/09 New molecule indexing arrays added
!   9/10 Initialisations now performed in a subroutine
!   2/19 Rigid molecules added
!   8/19 molatom option added
!  10/19 nmola2f added
!  12/19 Second derivatives added
!  12/19 nmolasymptr added
!   2/20 nmoleqv added
!   3/20 molQeig and molQsym added
!   4/20 molI added
!   4/20 Restart arrays added for rigid molecules
!   5/20 molaxes added
!   5/20 Group velocity arrays added
!   6/20 nmolcore arrays added
!  12/20 nmolconnectptr added
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
!  Julian Gale, CIC, Curtin University, December 2020
!
  use configurations, only : maxcfg
  use control,        only : lrigid, lgroupvelocity
  use molecule
  use derivatives,    only : molQdrv, molTdrv, molQTdrv, molQQdrv
  use derivatives,    only : molQCdrv, molTCdrv, molTTdrv, maxd2
  use derivatives,    only : molQSdrv, molTSdrv
  use derivatives,    only : molQTdri, molQQdri, molQCdri, molTCdri, molTTdri
  use derivatives,    only : molQTdk, molQQdk, molQCdk, molTCdk, molTTdk
  use reallocate
  implicit none
!
  integer(i4)       :: ierror, i
  integer(i4), save :: oldmaxmol = 0
!
  call realloc(molcom,3_i4,maxmol,ierror)
  if (ierror.ne.0) call outofmemory('changemaxmol','molcom')
  call realloc(molcoma,3_i4,maxmol,ierror)
  if (ierror.ne.0) call outofmemory('changemaxmol','molcoma')
  call realloc(moldim,maxmol,ierror)
  if (ierror.ne.0) call outofmemory('changemaxmol','moldim')
  call realloc(moldimi,maxmol,ierror)
  if (ierror.ne.0) call outofmemory('changemaxmol','moldimi')
  call realloc(molgcmc,maxmol,ierror)
  if (ierror.ne.0) call outofmemory('changemaxmol','molgcmc')
  call realloc(molI,3_i4,maxmol,ierror)
  if (ierror.ne.0) call outofmemory('changemaxmol','molI')
  call realloc(molQ,3_i4,maxmol,ierror)
  if (ierror.ne.0) call outofmemory('changemaxmol','molQ')
  call realloc(molQa,3_i4,maxmol,ierror)
  if (ierror.ne.0) call outofmemory('changemaxmol','molQa')
  call realloc(molQd1,3_i4,3_i4,3_i4,maxmol,ierror)
  if (ierror.ne.0) call outofmemory('changemaxmol','molQd1')
  call realloc(molQxyz,3_i4,maxmolat,maxmol,ierror)
  if (ierror.ne.0) call outofmemory('changemaxmol','molQxyz')
  call realloc(molxyz,3_i4,maxmolat,maxmol,ierror)
  if (ierror.ne.0) call outofmemory('changemaxmol','molxyz')
  call realloc(molcomcfg,3_i4,maxmol,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxmol','molcomcfg')
  call realloc(molQcfg,3_i4,maxmol,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxmol','molQcfg')
  call realloc(nmolasymno,maxmol,ierror)
  if (ierror.ne.0) call outofmemory('changemaxmol','nmolasymno')
  call realloc(nmolasymptr,maxmolat,maxmol,ierror)
  if (ierror.ne.0) call outofmemory('changemaxmol','nmolasymptr')
  call realloc(nmolasymeqvptr,maxmoleqv,maxmol,ierror)
  if (ierror.ne.0) call outofmemory('changemaxmol','nmolasymeqvptr')
  call realloc(nmola2f,maxmol,ierror)
  if (ierror.ne.0) call outofmemory('changemaxmol','nmola2f')
  call realloc(nmolf2a,maxmol,ierror)
  if (ierror.ne.0) call outofmemory('changemaxmol','nmolf2a')
  call realloc(nmoleqv,maxmol,ierror)
  if (ierror.ne.0) call outofmemory('changemaxmol','nmoleqv')
  call realloc(nmolatom,maxmol,ierror)
  if (ierror.ne.0) call outofmemory('changemaxmol','nmolatom')
  call realloc(nmolatomcfg,maxmol,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxmol','nmolatomcfg')
  call realloc(nmolcore,maxmol,ierror)
  if (ierror.ne.0) call outofmemory('changemaxmol','nmolcore')
  call realloc(nmolcorecfg,maxmol,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxmol','nmolcorecfg')
  call realloc(nmolptr,maxmol,ierror)
  if (ierror.ne.0) call outofmemory('changemaxmol','nmolptr')
  call realloc(nmolconnect,maxmol,ierror)
  if (ierror.ne.0) call outofmemory('changemaxmol','nmolconnect')
  call realloc(nmolconnectptr,maxconnectpermol,maxmol,ierror)
  if (ierror.ne.0) call outofmemory('changemaxmol','nmolconnectptr')
  call realloc(lgcmcmol,maxmol,ierror)
  if (ierror.ne.0) call outofmemory('changemaxmol','lgcmcmol')
!
!  Rigid molecule derivatives
!
  if (lrigid) then
    call realloc(molaxes,3_i4,3_i4,maxmol,ierror)
    if (ierror.ne.0) call outofmemory('changemaxmol','molaxes')
    call realloc(molQdrv,3_i4,maxmol,ierror)
    if (ierror.ne.0) call outofmemory('changemaxmol','molQdrv')
    call realloc(molTdrv,3_i4,maxmol,ierror)
    if (ierror.ne.0) call outofmemory('changemaxmol','molTdrv')
    call realloc(molQeig,3_i4,3_i4,maxmol,ierror)
    if (ierror.ne.0) call outofmemory('changemaxmol','molQeig')
    call realloc(molQsym,3_i4,3_i4,maxmol,ierror)
    if (ierror.ne.0) call outofmemory('changemaxmol','molQsym')
!
    call realloc(molQCdrv,maxd2,3_i4,maxmol,ierror)
    if (ierror.ne.0) call outofmemory('changemaxmol','molQCdrv')
    call realloc(molQCdri,maxd2,3_i4,maxmol,ierror)
    if (ierror.ne.0) call outofmemory('changemaxmol','molQCdri')
    call realloc(molQSdrv,6_i4,3_i4,maxmol,ierror)
    if (ierror.ne.0) call outofmemory('changemaxmol','molQSdrv')
    call realloc(molQQdrv,3_i4,3_i4,maxmol,maxmol,ierror)
    if (ierror.ne.0) call outofmemory('changemaxmol','molQQdrv')
    call realloc(molQQdri,3_i4,3_i4,maxmol,maxmol,ierror)
    if (ierror.ne.0) call outofmemory('changemaxmol','molQQdri')
    call realloc(molQTdrv,3_i4,3_i4,maxmol,maxmol,ierror)
    if (ierror.ne.0) call outofmemory('changemaxmol','molQTdrv')
    call realloc(molQTdri,3_i4,3_i4,maxmol,maxmol,ierror)
    if (ierror.ne.0) call outofmemory('changemaxmol','molQTdri')
    call realloc(molTCdrv,maxd2,3_i4,maxmol,ierror)
    if (ierror.ne.0) call outofmemory('changemaxmol','molTCdrv')
    call realloc(molTCdri,maxd2,3_i4,maxmol,ierror)
    if (ierror.ne.0) call outofmemory('changemaxmol','molTCdri')
    call realloc(molTSdrv,6_i4,3_i4,maxmol,ierror)
    if (ierror.ne.0) call outofmemory('changemaxmol','molTSdrv')
    call realloc(molTTdrv,3_i4,3_i4,maxmol,maxmol,ierror)
    if (ierror.ne.0) call outofmemory('changemaxmol','molTTdrv')
    call realloc(molTTdri,3_i4,3_i4,maxmol,maxmol,ierror)
    if (ierror.ne.0) call outofmemory('changemaxmol','molTTdri')
!
    if (lgroupvelocity) then
      call realloc(molQCdk,3_i4,maxd2,3_i4,maxmol,ierror)
      if (ierror.ne.0) call outofmemory('changemaxmol','molQCdk')
      call realloc(molQQdk,3_i4,3_i4,3_i4,maxmol,maxmol,ierror)
      if (ierror.ne.0) call outofmemory('changemaxmol','molQQdk')
      call realloc(molQTdk,3_i4,3_i4,3_i4,maxmol,maxmol,ierror)
      if (ierror.ne.0) call outofmemory('changemaxmol','molQTdk')
      call realloc(molTCdk,3_i4,maxd2,3_i4,maxmol,ierror)
      if (ierror.ne.0) call outofmemory('changemaxmol','molTCdk')
      call realloc(molTTdk,3_i4,3_i4,3_i4,maxmol,maxmol,ierror)
      if (ierror.ne.0) call outofmemory('changemaxmol','molTTdk')
    endif
  endif
!
!  Initialise new parts of data arrays
!
  if (maxmol.gt.oldmaxmol) then
    do i = oldmaxmol+1,maxmol
      call initmaxmoldefaults(i)
    enddo
  endif
!
!  Save current value of maxmol for next call
!
  oldmaxmol = maxmol
!
  return
  end
