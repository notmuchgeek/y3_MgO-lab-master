  subroutine changemaxEDIPspec
!
!  Alters the size of the arrays associated with maxEDIPspec
!
!   9/10 Created from changemaxreaxFFspec
!   1/14 EDIP2p, lEDIP3orig and lEDIPpi added
!   7/14 lEDIP3mod added
!   1/19 Use of general string reallocate added
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
!  Copyright Curtin University 2019
!
!  Julian Gale, CIC, Curtin University, January 2019
!
  use library
  use reallocate
  use EDIPdata
  implicit none
!
  integer(i4)       :: ierror, i
  integer(i4)       :: maxEDIPspec2
  integer(i4), save :: oldmaxEDIPspec = 0
!
!  Some things depend on pairs of species
!
  maxEDIPspec2 = maxEDIPspec*(maxEDIPspec + 1)/2
!
!  Species data
!
  call realloc_ch(5_i4,symbolEDIPspec,maxEDIPspec,ierror)
  if (ierror.ne.0) call outofmemory('changemaxEDIPspec','symbolEDIPspec')
  call realloc(natEDIPspec,maxEDIPspec,ierror)
  if (ierror.ne.0) call outofmemory('changemaxEDIPspec','natEDIPspec')
  call realloc(ntypEDIPspec,maxEDIPspec,ierror)
  if (ierror.ne.0) call outofmemory('changemaxEDIPspec','ntypEDIPspec')
  call realloc(lEDIPpairOK,maxEDIPspec2,ierror)
  if (ierror.ne.0) call outofmemory('changemaxEDIPspec','lEDIPpairOK')
  call realloc(lEDIPtriadOK,maxEDIPspec2,maxEDIPspec,ierror)
  if (ierror.ne.0) call outofmemory('changemaxEDIPspec','lEDIPtriadOK')
  call realloc(lEDIPpi,maxEDIPspec2,ierror)
  if (ierror.ne.0) call outofmemory('changemaxEDIPspec','lEDIPpi')
  call realloc(EDIPrmaxpair,maxEDIPspec2,ierror)
  if (ierror.ne.0) call outofmemory('changemaxEDIPspec','EDIPrmaxpair')
  call realloc(EDIPrmax,maxEDIPspec,ierror)
  if (ierror.ne.0) call outofmemory('changemaxEDIPspec','EDIPrmax')
  call realloc(EDIPfhigh,maxEDIPspec2,ierror)
  if (ierror.ne.0) call outofmemory('changemaxEDIPspec','EDIPfhigh')
  call realloc(EDIPflow,maxEDIPspec2,ierror)
  if (ierror.ne.0) call outofmemory('changemaxEDIPspec','EDIPflow')
  call realloc(EDIPphigh,maxEDIPspec2,ierror)
  if (ierror.ne.0) call outofmemory('changemaxEDIPspec','EDIPphigh')
  call realloc(EDIPplow,maxEDIPspec2,ierror)
  if (ierror.ne.0) call outofmemory('changemaxEDIPspec','EDIPplow')
  call realloc(EDIPalpha,maxEDIPspec2,ierror)
  if (ierror.ne.0) call outofmemory('changemaxEDIPspec','EDIPalpha')
  call realloc(EDIPZdih,maxEDIPspec2,ierror)
  if (ierror.ne.0) call outofmemory('changemaxEDIPspec','EDIPZdih')
  call realloc(EDIPZrep,maxEDIPspec2,ierror)
  if (ierror.ne.0) call outofmemory('changemaxEDIPspec','EDIPZrep')
  call realloc(EDIPc0,maxEDIPspec2,ierror)
  if (ierror.ne.0) call outofmemory('changemaxEDIPspec','EDIPc0')
  call realloc(EDIP2epsilon,maxEDIPspec2,ierror)
  if (ierror.ne.0) call outofmemory('changemaxEDIPspec','EDIP2epsilon')
  call realloc(EDIP2a,maxEDIPspec2,ierror)
  if (ierror.ne.0) call outofmemory('changemaxEDIPspec','EDIP2a')
  call realloc(EDIP2aprime,maxEDIPspec2,ierror)
  if (ierror.ne.0) call outofmemory('changemaxEDIPspec','EDIP2aprime')
  call realloc(EDIP2B,maxEDIPspec2,ierror)
  if (ierror.ne.0) call outofmemory('changemaxEDIPspec','EDIP2B')
  call realloc(EDIP2p,maxEDIPspec2,ierror)
  if (ierror.ne.0) call outofmemory('changemaxEDIPspec','EDIP2p')
  call realloc(EDIP2beta,maxEDIPspec2,ierror)
  if (ierror.ne.0) call outofmemory('changemaxEDIPspec','EDIP2beta')
  call realloc(EDIP2sigma,maxEDIPspec2,ierror)
  if (ierror.ne.0) call outofmemory('changemaxEDIPspec','EDIP2sigma')
  call realloc(lEDIP3mod,maxEDIPspec2,maxEDIPspec,ierror)
  if (ierror.ne.0) call outofmemory('changemaxEDIPspec','lEDIP3mod')
  call realloc(lEDIP3orig,maxEDIPspec2,maxEDIPspec,ierror)
  if (ierror.ne.0) call outofmemory('changemaxEDIPspec','lEDIP3orig')
  call realloc(EDIP3lambda0,maxEDIPspec2,maxEDIPspec,ierror)
  if (ierror.ne.0) call outofmemory('changemaxEDIPspec','EDIP3lambda0')
  call realloc(EDIP3lambdap,maxEDIPspec2,maxEDIPspec,ierror)
  if (ierror.ne.0) call outofmemory('changemaxEDIPspec','EDIP3lambdap')
  call realloc(EDIP3gamma0,maxEDIPspec2,maxEDIPspec,ierror)
  if (ierror.ne.0) call outofmemory('changemaxEDIPspec','EDIP3gamma0')
  call realloc(EDIP3gammap,maxEDIPspec2,maxEDIPspec,ierror)
  if (ierror.ne.0) call outofmemory('changemaxEDIPspec','EDIP3gammap')
  call realloc(EDIP3Z0,maxEDIPspec2,maxEDIPspec,ierror)
  if (ierror.ne.0) call outofmemory('changemaxEDIPspec','EDIP3Z0')
  call realloc(EDIP3q,maxEDIPspec2,maxEDIPspec,ierror)
  if (ierror.ne.0) call outofmemory('changemaxEDIPspec','EDIP3q')
  call realloc(EDIP3kq2,maxEDIPspec2,maxEDIPspec,ierror)
  if (ierror.ne.0) call outofmemory('changemaxEDIPspec','EDIP3kq2')
!
!  Initialise new parts of data arrays
!
  if (maxEDIPspec.gt.oldmaxEDIPspec) then
    do i = oldmaxEDIPspec+1,maxEDIPspec
      call initmaxEDIPspecdefaults(i)
    enddo
  endif
!
!  Save current value of maxEDIPspec for next call
!
  oldmaxEDIPspec = maxEDIPspec
!
  return
  end
