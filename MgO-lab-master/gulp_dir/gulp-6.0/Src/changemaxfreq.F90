  subroutine changemaxfreq(nfq)
!
!  Alters the size of the arrays associated with maxfreq
!
!   4/11 Created from changemaxat
!   8/14 eigv and groupvelocity added
!   8/16 IR intensity array added
!   1/18 Grueneisen parameters added
!   4/20 maxfqat changed to maxfreq since it is no longer just atom based
!   4/20 fmass added
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
  use datatypes
  use control,        only : lgroupvelocity, lgrueneisen
  use frequencies,    only : freq, eigv, maxfkpt, maxfreq, lStoreEig, IRintensity
  use frequencies,    only : groupvelocity, grueneisen, fmass, rfmass
  use reallocate
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in) :: nfq
!
!  Local variables
!
  integer(i4)             :: ierror
!
  if (nfq.gt.maxfreq) then
    maxfreq = nfq
!
    call realloc(fmass,maxfreq,ierror)
    if (ierror.ne.0) call outofmemory('changemaxfreq','fmass')
    call realloc(rfmass,maxfreq,ierror)
    if (ierror.ne.0) call outofmemory('changemaxfreq','rfmass')
    call realloc(freq,maxfreq,maxfkpt,ierror)
    if (ierror.ne.0) call outofmemory('changemaxfreq','freq')
    call realloc(IRintensity,maxfreq,maxfkpt,ierror)
    if (ierror.ne.0) call outofmemory('changemaxfreq','IRintensity')
!
    if (lStoreEig) then
      call realloc(eigv,maxfreq,maxfreq,maxfkpt,ierror)
      if (ierror.ne.0) call outofmemory('changemaxfreq','eigv')
    endif
    if (lgroupvelocity) then
      call realloc(groupvelocity,3_i4,maxfreq,maxfkpt,ierror)
      if (ierror.ne.0) call outofmemory('changemaxfreq','groupvelocity')
    endif
    if (lgrueneisen) then
      call realloc(grueneisen,maxfreq,maxfkpt,ierror)
      if (ierror.ne.0) call outofmemory('changemaxfreq','grueneisen')
    endif
  endif
!
  return
  end
