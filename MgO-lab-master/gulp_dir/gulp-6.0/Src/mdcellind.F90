  subroutine mdcellind(ni,icos1,icos2,icos3)
!
!  Corrects cell indices for MD based on
!  changes in directional cosines. Applies
!  to molecule indices and list terms.
!
!   3/95 Bonding list indices now corrected
!   5/96 Bug in bonded index changes fixed
!  10/06 Fourlist handling modified in line with changes
!   4/07 Searching for bonded atoms now uses nbonds number rather
!        than checking for a non-zero type
!   6/09 Module name changed from three to m_three
!  10/11 Corrections made to handling of cell indices for fourbody potentials
!   9/13 Modified to speed up index update for atoms bonded to ni
!   9/13 Only update the 3 and 4 body list terms that will be computed on the local node
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
  use control
  use current
  use four
  use m_three
  use molecule
  use parallel,       only : nprocs, procid
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in) :: ni         ! Atom whose index has changed
  integer(i4), intent(in) :: icos1
  integer(i4), intent(in) :: icos2
  integer(i4), intent(in) :: icos3
!
!  Local variables
!
  integer(i4)             :: i
  integer(i4)             :: icm
  integer(i4)             :: ii
  integer(i4)             :: ind
  integer(i4)             :: ixd
  integer(i4)             :: iyd
  integer(i4)             :: izd
  integer(i4)             :: j
  integer(i4)             :: jcm
  integer(i4)             :: k
  integer(i4)             :: l
  integer(i4)             :: n
  integer(i4)             :: numat1
  logical                 :: lfound
  logical                 :: llist
#ifdef TRACE
  call trace_in('mdcellind')
#endif
!
!  See if cosines have changed
!
  ixd = icos1 - icosx(ni)
  iyd = icos2 - icosy(ni)
  izd = icos3 - icosz(ni)
!
!  If no change then return
!
  if ((abs(ixd)+abs(iyd)+abs(izd)).eq.0) then
#ifdef TRACE
    call trace_out('mdcellind')
#endif
    return
  endif
!
  llist = (index(keyword,'noli').eq.0)
  if (natmol(ni).gt.0) then
!*************
!  Molecule  *
!*************
    ind = nmolind(ni)
    call indsft(ind,ixd,iyd,izd,1_i4)
    nmolind(ni) = ind
!************
!  Bonding  *
!************
!
!  Change in bonding index for other atoms
!
    do icm = 1,nbonds(ni)
      i = nbonded(icm,ni)
      lfound = .false.
      jcm = 0
      do while (.not.lfound.and.jcm.lt.nbonds(i))
        jcm = jcm + 1
        if (nbonded(jcm,i).eq.ni) then
          lfound = .true.
          ind = nbondind(jcm,i)
          call indsft(ind,ixd,iyd,izd,1_i4)
          nbondind(jcm,i) = ind
        endif
      enddo
    enddo
!
!  Change in bonding indices for ni
!
    icm = nbonded(1,ni)
    jcm = 1
    do while (icm.gt.0.and.jcm.le.nbonds(ni))
      ind = nbondind(jcm,ni)
      call indsft(ind,ixd,iyd,izd,-1_i4)
      nbondind(jcm,ni) = ind
      jcm = jcm + 1
      icm = nbonded(jcm,ni)
    enddo
  endif
!
  if (llist.and.nlist3md.gt.0) then
!********************
!  Three-body list  *
!********************
    do i = procid+1,nlist3md,nprocs
      ii = i3ind(i)
      j = j3ind(i)
      k = k3ind(i)
      if (ni.eq.j.or.ni.eq.k) then
        if (ni.eq.j) then
          ind = icell31(i)
        else
          ind = icell32(i)
        endif
        call indsft(ind,ixd,iyd,izd,1_i4)
        if (ni.eq.j) then
          icell31(i) = ind
        else
          icell32(i) = ind
        endif
      elseif (ni.eq.ii) then
        ind = icell31(i)
        call indsft(ind,ixd,iyd,izd,-1_i4)
        icell31(i) = ind
        ind = icell32(i)
        call indsft(ind,ixd,iyd,izd,-1_i4)
        icell32(i) = ind
      endif
    enddo
  endif
!
  if (llist.and.nlist4md.gt.0) then
    numat1 = (numat + 1)
!*******************
!  Four-body list  *
!*******************
    do i = procid+1,nlist4md,nprocs
      n = nforptr(i)
      ind = ilind(i)
      l = ind/numat1
      ii = ind - l*numat1
      ind = jkind(i)
      k = ind/numat1
      j = ind - k*numat1
      if (loutofplane(n)) then
!
!  Out of plane
!
        if (ni.eq.ii) then
          ind = icell41(i)
          call indsft(ind,ixd,iyd,izd,-1_i4)
          icell41(i) = ind
          ind = icell42(i)
          call indsft(ind,ixd,iyd,izd,-1_i4)
          icell42(i) = ind
          ind = icell43(i)
          call indsft(ind,ixd,iyd,izd,-1_i4)
          icell43(i) = ind
        elseif (ni.eq.j) then
          ind = icell41(i)
          call indsft(ind,ixd,iyd,izd,1_i4)
          icell41(i) = ind
        elseif (ni.eq.k) then
          ind = icell42(i)
          call indsft(ind,ixd,iyd,izd,1_i4)
          icell42(i) = ind
        elseif (ni.eq.l) then
          ind = icell43(i)
          call indsft(ind,ixd,iyd,izd,1_i4)
          icell43(i) = ind
        endif
      else
!
!  Standard torsional potential
!
        if (ni.eq.ii) then
          ind = icell41(i)
          call indsft(ind,ixd,iyd,izd,1_i4)
          icell41(i) = ind
        elseif (ni.eq.j) then
          ind = icell41(i)
          call indsft(ind,ixd,iyd,izd,-1_i4)
          icell41(i) = ind
          ind = icell42(i)
          call indsft(ind,ixd,iyd,izd,-1_i4)
          icell42(i) = ind
          ind = icell43(i)
          call indsft(ind,ixd,iyd,izd,-1_i4)
          icell43(i) = ind
        elseif (ni.eq.k) then
          ind = icell42(i)
          call indsft(ind,ixd,iyd,izd,1_i4)
          icell42(i) = ind
        elseif (ni.eq.l) then
          ind = icell43(i)
          call indsft(ind,ixd,iyd,izd,1_i4)
          icell43(i) = ind
        endif
      endif
    enddo
  endif
!
!********************************
!  Correct stored cell indices  *
!********************************
  icosx(ni) = icos1
  icosy(ni) = icos2
  icosz(ni) = icos3
#ifdef TRACE
  call trace_out('mdcellind')
#endif
!
  return
  end
