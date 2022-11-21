  subroutine reorder(nc,iptr)
!
!  Reorders atoms according to a pointer which is passed as an argument
!
!  On entry :
!
!  nc   = configuration to reorder
!  iptr = pointer indicating the new atom number for each atom
!
!  11/01 Created from sort.f
!  11/01 Slice atom logical added to data to reorder
!   8/02 Reordering of forcecfg added
!   2/04 Reordering of time dependent force data added
!   4/04 Reordering of ltranat added
!   5/04 Reordering of NEB final coordinates added
!   8/04 Trap for zero atom case added
!   5/06 Handling of nspecptr added
!  11/06 Reordering of NEB final coordinates added
!  11/06 nebfinalradius added
!   2/07 Connectivity information reordering added
!   5/07 Out of bounds check added for connectivity atom numbers
!   7/15 External potential added
!   2/18 Trace added
!   3/18 Sign option added to translate
!   3/19 Constraint arrays changed to have index and type
!   8/19 molatom option added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!   7/20 Handling of nmollistcfg corrected for symmetry case
!   7/20 Modifications for gfortran v10
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
!  Julian Gale, CIC, Curtin University, July 2020
!
  use control
  use configurations
  use current
  use defects
  use moldyn,          only : lfix, lmdconstrain, nmdconstrainatom
  use molecule
  use g_neb,           only : nebfinalxyz, nebfinalradius
  use observables
  use optimisation
  use projectdos
  use scan,            only : ltranat, ltranatminus
#ifdef TRACE
  use trace,           only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4)                                  :: nc
  integer(i4)                                  :: iptr(*)
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ii
  integer(i4)                                  :: ij
  integer(i4)                                  :: ind
  integer(i4), dimension(:), allocatable       :: itmp
  integer(i4), dimension(:), allocatable       :: itmp2
  integer(i4)                                  :: j
  integer(i4)                                  :: k
  integer(i4)                                  :: n
  integer(i4)                                  :: naddi
  integer(i4)                                  :: naddj
  integer(i4)                                  :: naddk
  integer(i4)                                  :: ni
  integer(i4)                                  :: ni2
  integer(i4)                                  :: nj
  integer(i4)                                  :: nj2
  integer(i4)                                  :: nk
  integer(i4)                                  :: nk2
  integer(i4)                                  :: noffset
  integer(i4)                                  :: np
  integer(i4)                                  :: npc
  integer(i4)                                  :: npi
  integer(i4)                                  :: npfirst
  integer(i4)                                  :: nplast
  integer(i4)                                  :: npifirst
  integer(i4)                                  :: npilast
  integer(i4)                                  :: nsum
  integer(i4), dimension(:), allocatable       :: nrela2fsorted
  integer(i4)                                  :: status
  logical,     dimension(:), allocatable       :: ltmp
  logical,     dimension(:), allocatable       :: ltmp2
  real(dp),    dimension(:), allocatable       :: rtmp
  real(dp),    dimension(:), allocatable       :: tmp
!
!  Get number of atoms for this configuration
!
  nasym = nascfg(nc)
!
!  If number of atoms is zero then return as there is nothing to do
!
  if (nasym.eq.0) return
#ifdef TRACE
  call trace_in('reorder')
#endif
!
!  Allocate local memory
!
  allocate(itmp(nasym),stat=status)
  if (status/=0) call outofmemory('reorder','itmp')
  allocate(itmp2(nasym),stat=status)
  if (status/=0) call outofmemory('reorder','itmp2')
  allocate(rtmp(nasym),stat=status)
  if (status/=0) call outofmemory('reorder','rtmp')
  allocate(tmp(nasym),stat=status)
  if (status/=0) call outofmemory('reorder','tmp')
  allocate(ltmp(3*nasym),stat=status)
  if (status/=0) call outofmemory('reorder','ltmp')
  allocate(ltmp2(nasym),stat=status)
  if (status/=0) call outofmemory('reorder','ltmp2')
!
!  Swap simple atom referenced data around
!
  itmp(1:nasym) = natcfg(nsft+1:nsft+nasym)
  call icollect(nasym,itmp,itmp2,iptr)
  natcfg(nsft+1:nsft+nasym) = itmp(1:nasym)
!
  itmp(1:nasym) = ntypcfg(nsft+1:nsft+nasym)
  call icollect(nasym,itmp,itmp2,iptr)
  ntypcfg(nsft+1:nsft+nasym) = itmp(1:nasym)
!
  itmp(1:nasym) = nregionno(nsft+1:nsft+nasym)
  call icollect(nasym,itmp,itmp2,iptr)
  nregionno(nsft+1:nsft+nasym) = itmp(1:nasym)
!
  itmp(1:nasym) = nspecptrcfg(nsft+1:nsft+nasym)
  call icollect(nasym,itmp,itmp2,iptr)
  nspecptrcfg(nsft+1:nsft+nasym) = itmp(1:nasym)
!
  rtmp(1:nasym) = cncfg(nsft+1:nsft+nasym)
  call collect(nasym,rtmp,tmp,iptr)
  cncfg(nsft+1:nsft+nasym) = rtmp(1:nasym)
!
  rtmp(1:nasym) = xcfg(nsft+1:nsft+nasym)
  call collect(nasym,rtmp,tmp,iptr)
  xcfg(nsft+1:nsft+nasym) = rtmp(1:nasym)
!
  rtmp(1:nasym) = ycfg(nsft+1:nsft+nasym)
  call collect(nasym,rtmp,tmp,iptr)
  ycfg(nsft+1:nsft+nasym) = rtmp(1:nasym)
!
  rtmp(1:nasym) = zcfg(nsft+1:nsft+nasym)
  call collect(nasym,rtmp,tmp,iptr)
  zcfg(nsft+1:nsft+nasym) = rtmp(1:nasym)
!
  rtmp(1:nasym) = qlcfg(nsft+1:nsft+nasym)
  call collect(nasym,rtmp,tmp,iptr)
  qlcfg(nsft+1:nsft+nasym) = rtmp(1:nasym)
!
  rtmp(1:nasym) = occucfg(nsft+1:nsft+nasym)
  call collect(nasym,rtmp,tmp,iptr)
  occucfg(nsft+1:nsft+nasym) = rtmp(1:nasym)
!
  rtmp(1:nasym) = oxcfg(nsft+1:nsft+nasym)
  call collect(nasym,rtmp,tmp,iptr)
  oxcfg(nsft+1:nsft+nasym) = rtmp(1:nasym)
!
  rtmp(1:nasym) = radcfg(nsft+1:nsft+nasym)
  call collect(nasym,rtmp,tmp,iptr)
  radcfg(nsft+1:nsft+nasym) = rtmp(1:nasym)
!
  ltmp(1:nasym) = lbsmat(nsft+1:nsft+nasym)
  call lcollect(nasym,ltmp,ltmp2,iptr)
  lbsmat(nsft+1:nsft+nasym) = ltmp(1:nasym)
!
  ltmp(1:nasym) = lqmatom(nsft+1:nsft+nasym)
  call lcollect(nasym,ltmp,ltmp2,iptr)
  lqmatom(nsft+1:nsft+nasym) = ltmp(1:nasym)
!
  ltmp(1:nasym) = lsliceatom(nsft+1:nsft+nasym)
  call lcollect(nasym,ltmp,ltmp2,iptr)
  lsliceatom(nsft+1:nsft+nasym) = ltmp(1:nasym)
!
  ltmp(1:nasym) = ltranat(nsft+1:nsft+nasym)
  call lcollect(nasym,ltmp,ltmp2,iptr)
  ltranat(nsft+1:nsft+nasym) = ltmp(1:nasym)
!
  ltmp(1:nasym) = ltranatminus(nsft+1:nsft+nasym)
  call lcollect(nasym,ltmp,ltmp2,iptr)
  ltranatminus(nsft+1:nsft+nasym) = ltmp(1:nasym)
!
  ltmp(1:nasym) = lfix(nsft+1:nsft+nasym)
  call lcollect(nasym,ltmp,ltmp2,iptr)
  lfix(nsft+1:nsft+nasym) = ltmp(1:nasym)
!
!  Swap data with more complex dependendancy on the atom numbers
!
!  Change MD constraint if present
!
  if (lmdconstrain(nc)) then
    nmdconstrainatom(1,nc) = iptr(nmdconstrainatom(1,nc))
    nmdconstrainatom(2,nc) = iptr(nmdconstrainatom(2,nc))
  endif
!
!  Change around lopfi
!
  noffset = 3*nsft
  do i = 1,3*nasym
    ltmp(i) = lopfi(noffset+i)
  enddo
  do i = 1,nasym
    lopfi(noffset+3*(i-1)+1) = ltmp(3*(iptr(i)-1)+1)
    lopfi(noffset+3*(i-1)+2) = ltmp(3*(iptr(i)-1)+2)
    lopfi(noffset+3*(i-1)+3) = ltmp(3*(iptr(i)-1)+3)
  enddo
!
!  Change NEB final coordinates and radii
!
  do i = 1,3
    do j = 1,nasym
      rtmp(j) = nebfinalxyz(i,j,nc)
    enddo
    call collect(nasym,rtmp,tmp,iptr)
    do j = 1,nasym
      nebfinalxyz(i,j,nc) = rtmp(j)
    enddo
  enddo
  do j = 1,nasym
    rtmp(j) = nebfinalradius(j,nc)
  enddo
  call collect(nasym,rtmp,tmp,iptr)
  do j = 1,nasym
    nebfinalradius(j,nc) = rtmp(j)
  enddo
!
!  Change external potentials
!
  do j = 1,nasym
    rtmp(j) = extpotcfg(nsft+j)
  enddo
  call collect(nasym,rtmp,tmp,iptr)
  do j = 1,nasym
    extpotcfg(nsft+j) = rtmp(j)
  enddo
!
!  Change external forces
!
  do i = 1,3
    do j = 1,nasym
      rtmp(j) = forcecfg(i,nsft+j)
    enddo
    call collect(nasym,rtmp,tmp,iptr)
    do j = 1,nasym
      forcecfg(i,nsft+j) = rtmp(j)
    enddo
  enddo
  do i = 1,3
    do j = 1,nasym
      ltmp(j) = ltdforcecfg(i,nsft+j)
    enddo
    call lcollect(nasym,ltmp,ltmp(nasym+1),iptr)
    do j = 1,nasym
      ltdforcecfg(i,nsft+j) = ltmp(j)
    enddo
  enddo
  do i = 1,3
    do k = 1,3
      do j = 1,nasym
        rtmp(j) = tdforcecfg(k,i,nsft+j)
      enddo
      call collect(nasym,rtmp,tmp,iptr)
      do j = 1,nasym
        tdforcecfg(k,i,nsft+j) = rtmp(j)
      enddo
    enddo
  enddo
!
!  Change defect atom pointers
!
  do i = 1,ndef
    if (ndefcfg(i).eq.nc) then
      if (ndeftyp(i).eq.1.or.ndeftyp(i).eq.11.or.ndeftyp(i).eq.12) then
        ni = nint(xdef(i))
        do k = 1,nasym
          if (iptr(k).eq.ni) then
            xdef(i) = k
          endif
        enddo
      endif
    endif
  enddo
!
!  Change defect centre pointer
!
  if (ndcentyp(nc).eq.1.or.ndcentyp(nc).eq.2) then
    ni = nint(xdcent(nc))
    do k = 1,nasym
      if (iptr(k).eq.ni) then
        xdcent(nc) = k
      endif
    enddo
  endif
!
!  Change fitted gradient pointers
!
  do i = 1,nfgrad
    if (nfgracfg(i).eq.nc) then
      ni = nfgrat(i)
      do k = 1,nasym
        if (iptr(k).eq.ni) then
          nfgrat(i) = k
        endif
      enddo
    endif
  enddo
!
!  Change any constraints
!
  if (ncontot.gt.0) then
    do i = 1,ncontot
      if (nconcfg(i).eq.nc) then
        if (ncfixtypcfg(i).eq.iopt_radius) then
!
!  Breathing shell constraint
!
          ii = ncfixindcfg(i)
          ij = ncvarindcfg(i)
          do j = 1,nasym
            if (ii.eq.iptr(j)) then
              ncfixindcfg(i) = ncfixindcfg(i) + j - ii
            elseif (ij.eq.iptr(j)) then
              ncvarindcfg(i) = ncvarindcfg(i) +j - ij
            endif
          enddo
        elseif (ncfixtypcfg(i).eq.iopt_xf.or.ncfixtypcfg(i).eq.iopt_yf.or.ncfixtypcfg(i).eq.iopt_zf) then
!
!  Coordinate constraint
!
          ii = ncfixindcfg(i)
          ij = ncvarindcfg(i)
          do j = 1,nasym
            if (ii.eq.iptr(j)) then
              ncfixindcfg(i) = ncfixindcfg(i) + j - ii
            elseif (ij.eq.iptr(j)) then
              ncvarindcfg(i) = ncvarindcfg(i) + j - ij
            endif
          enddo
        endif
      endif
    enddo
  endif
!
!  Change connectivity/molecule pointers
!
  if (nconnect.gt.0.or.nmolcfg(nc).gt.0) then
    if (nasym.ne.numat) then
!
!  Count position of asymmetric unit atoms in full cell for symmetry case
!
      allocate(nrela2fsorted(nasym),stat=status)
      if (status/=0) call outofmemory('reorder','nrela2fsorted')
      do i = 1,nasym
        nrela2fsorted(iptr(i)) = neqv(i)
      enddo
      nsum = 0
      do i = 1,nasym
        ni2 = nrela2fsorted(i)
        nrela2fsorted(i) = nsum + 1
        nsum = nsum + ni2
      enddo
    endif
    if (nconnect.gt.0) then
      do n = 1,nconnect
        if (nconnectcfg(n).eq.nc) then
          i = n1connect(n)
          j = n2connect(n)
!
!  Check that i & j are in bounds
!
          if (i.gt.numat.or.j.gt.numat) then
            call outerror('Atom number in connectivity is greater than number of atoms',0_i4)
            call stopnow('reorder')
          endif
          if (nasym.eq.numat) then
!
!  Simple no symmetry sort
!
            n1connect(n) = iptr(i)
            n2connect(n) = iptr(j)
          else
!
!  More complex symmetry adapted sort
!
            ni = nrelf2a(i)
            nj = nrelf2a(j)
            ni2 = iptr(ni)
            nj2 = iptr(nj)
            naddi = i - nrela2f(ni)
            naddj = j - nrela2f(nj)
            n1connect(n) = nrela2fsorted(ni2) + naddi
            n2connect(n) = nrela2fsorted(nj2) + naddj
          endif
        endif
      enddo
    endif
!
!  Change molatom pointers
!
    if (nmolcfg(nc).gt.0) then
      ind = 0
      do i = 1,nmolcfg(nc)
        do j = 1,nmolatomcfg(i,nc)
          k = nmollistcfg(ind+j,nc)
!
!  Check that atom is within bounds
!
          if (k.gt.numat) then
            call outerror('Atom number in nmollistcfg is greater than number of atoms',0_i4)
            call stopnow('reorder')
          endif
          if (nasym.eq.numat) then
!
!  Simple no symmetry sort
!
            nmollistcfg(ind+j,nc) = iptr(k)
          else
!
!  More complex symmetry adapted sort
!
            nk = nrelf2a(k)
            nk2 = iptr(nk)
            naddk = k - nrela2f(nk)
            nmollistcfg(ind+j,nc) = nrela2fsorted(nk2) + naddk
          endif
        enddo
        ind = ind + nmolatomcfg(i,nc)
      enddo
    endif
!
    if (nasym.ne.numat) then
      deallocate(nrela2fsorted,stat=status)
      if (status/=0) call deallocate_error('reorder','nrela2fsorted')
    endif
  endif
!
!  Change weights - code not correct because of cell constraints
!  also may not be worth including as people work out which
!  observable to weight by observable table
!
!      if (lfit) then
!        nlower = noffset + nobs + 6*nc
!        nupper = nlower + 3*nasym
!        do i = nlower+1,nupper
!          tmp(i-nlower) = weight(i)
!        enddo
!        do i = 1,nasym
!          weight(nlower+3*(i-1)+1) = tmp(3*(iptr(i)-1)+1)
!          weight(nlower+3*(i-1)+2) = tmp(3*(iptr(i)-1)+2)
!          weight(nlower+3*(i-1)+3) = tmp(3*(iptr(i)-1)+3)
!        enddo
!      endif
!
!  Change projections
!
  if ((nprojcfg(nc)-nprojdef(nc)).gt.0) then
    npfirst = 1
    npifirst = 1
    ii = 0
    do i = 1,nc-1
      npc = nprojcfg(i)
      npfirst = npfirst + npc
      do j = 1,npc
        npifirst = npifirst + nprojit(ii+j)
      enddo
      ii = ii + npc
    enddo
    nplast = npfirst + nprojcfg(nc) - 1
    npilast = npifirst
    do i = 1,nprojcfg(nc)
      npilast = npilast + nprojit(ii+i)
    enddo
    npilast = npilast - 1
    do np = npfirst,nplast
      if (nprojdb(np).eq.1) then
        do npi = npifirst,npilast
          if (nprojptr(npi).eq.np) then
            if (nprojtyp(npi).gt.99) then
              ind = nprojnat(npi)
              do k = 1,nasym
                if (iptr(k).eq.ind) nprojnat(npi) = k
              enddo
            endif
          endif
        enddo
      endif
    enddo
  endif
!
!  Free local memory
!
  deallocate(ltmp2,stat=status)
  if (status/=0) call deallocate_error('reorder','ltmp2')
  deallocate(ltmp,stat=status)
  if (status/=0) call deallocate_error('reorder','ltmp')
  deallocate(tmp,stat=status)
  if (status/=0) call deallocate_error('reorder','tmp')
  deallocate(rtmp,stat=status)
  if (status/=0) call deallocate_error('reorder','rtmp')
  deallocate(itmp2,stat=status)
  if (status/=0) call deallocate_error('reorder','itmp2')
  deallocate(itmp,stat=status)
  if (status/=0) call deallocate_error('reorder','itmp')
#ifdef TRACE
  call trace_out('reorder')
#endif
!
  return
  end
