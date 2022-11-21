  subroutine changemaxr1at
!
!  Alters the size of the arrays associated with maxr1at
!
!   4/04 Maximum dimension of idopt increased to allow for breathing shells
!   5/08 Defect bonding arrays changed to match structure of perfect ones
!  12/08 scrho, and defect versions, converted to 2-D array for MEAM
!   9/10 Initialisations now performed in a subroutine
!   3/17 Parallel variable arrays added
!   5/17 Parallel arrays added for nreg1 distribution over nodes
!   8/17 Modified so that changemaxr1at and changemaxat don't interfere
!        with the size of common arrays
!   3/19 Change of idopt to idoptindex and idopttype
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
!  Julian Gale, CIC, Curtin University, March 2019
!
  use current,    only : maxbond
  use defects
  use eam,        only : maxmeamcomponent
  use parallel,   only : node2var, nvar2node, nvar2local
  use parallel,   only : node2reg1, reg12node, reg12local
  use reallocate
  implicit none
!
  integer(i4)       :: ierror, i
  integer(i4), save :: oldmaxr1at = 0
!
  call realloc(idoptindex,4_i4*maxr1at,ierror)
  if (ierror.ne.0) call outofmemory('changemaxr1at','idoptindex')
  call realloc(idopttype,4_i4*maxr1at,ierror)
  if (ierror.ne.0) call outofmemory('changemaxr1at','idopttype')
  call realloc(inddfix,maxr1at,ierror)
  if (ierror.ne.0) call outofmemory('changemaxr1at','inddfix')
  call realloc(ldefbsmat,maxr1at,ierror)
  if (ierror.ne.0) call outofmemory('changemaxr1at','ldefbsmat')
  call realloc(ldfix,maxr1at,ierror)
  if (ierror.ne.0) call outofmemory('changemaxr1at','ldfix')
  call realloc(lr1created,maxr1at,ierror)
  if (ierror.ne.0) call outofmemory('changemaxr1at','lr1created')
  call realloc(ndefmol,maxr1at,ierror)
  if (ierror.ne.0) call outofmemory('changemaxr1at','ndefmol')
  call realloc(ndefmolp,maxr1at,ierror)
  if (ierror.ne.0) call outofmemory('changemaxr1at','ndefmolp')
  call realloc(ldqmatom,maxr1at,ierror)
  if (ierror.ne.0) call outofmemory('changemaxr1at','ldqmatom')
  call realloc(natdefe,maxr1at,ierror)
  if (ierror.ne.0) call outofmemory('changemaxr1at','natdefe')
  call realloc(natp,maxr1at,ierror)
  if (ierror.ne.0) call outofmemory('changemaxr1at','natp')
  call realloc(nbondsdef,maxr1at,ierror)
  if (ierror.ne.0) call outofmemory('changemaxr1at','nbondsdef')
  call realloc(nbondeddef,maxbond,maxr1at,ierror)
  if (ierror.ne.0) call outofmemory('changemaxr1at','nbondeddef')
  call realloc(ndefind,maxr1at,ierror)
  if (ierror.ne.0) call outofmemory('changemaxr1at','ndefind')
  call realloc(ndefindp,maxr1at,ierror)
  if (ierror.ne.0) call outofmemory('changemaxr1at','ndefindp')
  call realloc(ndeqv,maxr1at,ierror)
  if (ierror.ne.0) call outofmemory('changemaxr1at','ndeqv')
  call realloc(ndrel,maxr1at,ierror)
  if (ierror.ne.0) call outofmemory('changemaxr1at','ndrel')
  call realloc(ndrelop,maxr1at,ierror)
  if (ierror.ne.0) call outofmemory('changemaxr1at','ndrelop')
  call realloc(ndsptr,maxr1at,ierror)
  if (ierror.ne.0) call outofmemory('changemaxr1at','ndsptr')
  call realloc(npsite,maxr1at,ierror)
  if (ierror.ne.0) call outofmemory('changemaxr1at','npsite')
  call realloc(nptrr1,maxr1at,ierror)
  if (ierror.ne.0) call outofmemory('changemaxr1at','nptrr1')
  call realloc(nreldef,maxr1at,ierror)
  if (ierror.ne.0) call outofmemory('changemaxr1at','nreldef')
  call realloc(ntypdefe,maxr1at,ierror)
  if (ierror.ne.0) call outofmemory('changemaxr1at','ntypdefe')
  call realloc(ntypep,maxr1at,ierror)
  if (ierror.ne.0) call outofmemory('changemaxr1at','ntypep')
  call realloc(dscrho,maxmeamcomponent,maxr1at,ierror)
  if (ierror.ne.0) call outofmemory('changemaxr1at','dscrho')
  call realloc(occdefe,maxr1at,ierror)
  if (ierror.ne.0) call outofmemory('changemaxr1at','occdefe')
  call realloc(occp,maxr1at,ierror)
  if (ierror.ne.0) call outofmemory('changemaxr1at','occp')
  call realloc(qdefe,maxr1at,ierror)
  if (ierror.ne.0) call outofmemory('changemaxr1at','qdefe')
  call realloc(qp,maxr1at,ierror)
  if (ierror.ne.0) call outofmemory('changemaxr1at','qp')
  call realloc(radefe,maxr1at,ierror)
  if (ierror.ne.0) call outofmemory('changemaxr1at','radefe')
  call realloc(xdefe,maxr1at,ierror)
  if (ierror.ne.0) call outofmemory('changemaxr1at','xdefe')
  call realloc(ydefe,maxr1at,ierror)
  if (ierror.ne.0) call outofmemory('changemaxr1at','ydefe')
  call realloc(zdefe,maxr1at,ierror)
  if (ierror.ne.0) call outofmemory('changemaxr1at','zdefe')
  call realloc(xperf,maxr1at,ierror)
  if (ierror.ne.0) call outofmemory('changemaxr1at','xperf')
  call realloc(yperf,maxr1at,ierror)
  if (ierror.ne.0) call outofmemory('changemaxr1at','yperf')
  call realloc(zperf,maxr1at,ierror)
  if (ierror.ne.0) call outofmemory('changemaxr1at','zperf')
!
  if (4_i4*maxr1at.gt.size(node2var)) then
    call realloc(node2var,4_i4*maxr1at,ierror)
    if (ierror.ne.0) call outofmemory('changemaxr1at','node2var')
    call realloc(nvar2local,4_i4*maxr1at,ierror)
    if (ierror.ne.0) call outofmemory('changemaxr1at','nvar2local')
    call realloc(nvar2node,4_i4*maxr1at,ierror)
    if (ierror.ne.0) call outofmemory('changemaxr1at','nvar2node')
  endif
!
  call realloc(node2reg1,maxr1at,ierror)
  if (ierror.ne.0) call outofmemory('changemaxr1at','node2reg1')
  call realloc(reg12local,maxr1at,ierror)
  if (ierror.ne.0) call outofmemory('changemaxr1at','reg12local')
  call realloc(reg12node,maxr1at,ierror)
  if (ierror.ne.0) call outofmemory('changemaxr1at','reg12node')
!
!  Initialise new parts of data arrays
!
  if (maxr1at.gt.oldmaxr1at) then
    do i = oldmaxr1at+1,maxr1at
      call initmaxr1atdefaults(i)
    enddo
  endif
!
!  Save current value of maxat for next call
!
  oldmaxr1at = maxr1at
!
  return
  end
