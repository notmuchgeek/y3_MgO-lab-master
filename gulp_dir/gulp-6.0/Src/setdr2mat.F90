  subroutine setdr2mat
!
!  Generate second derivative elements needed for region 2 energy
!  Now called from main routine so that the second derivative
!  matrices can be handled dynamically.
!
!  Channel 44 is used to save defect matrices
!
!   6/95 Modified to allow for additive defect constraints
!  11/95 Defect flags modified for correct restart
!   5/07 Handling of partial occupancy reduction of dervi added
!   5/17 Parallelisation added
!   6/17 Corrections to parallelisation made
!   1/18 Trace added
!   7/20 sendall2D added to handled communication of derv3
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
!  Julian Gale, Curtin University, July 2020
!
  use control
  use current
  use derivatives
  use iochannels
  use parallel
  use partial
  use properties
  use shells,      only : nshell, nshptr
#ifdef TRACE
  use trace,       only : trace_in, trace_out
#endif
  implicit none
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ii
  integer(i4)                                  :: is
  integer(i4)                                  :: ind
  integer(i4)                                  :: indf
  integer(i4)                                  :: indi
  integer(i4)                                  :: indj
  integer(i4)                                  :: indx
  integer(i4)                                  :: indy
  integer(i4)                                  :: indz
  integer(i4)                                  :: inds
  integer(i4)                                  :: indss
  integer(i4)                                  :: j
  integer(i4)                                  :: jj
  integer(i4)                                  :: js
  integer(i4)                                  :: k
  integer(i4)                                  :: naprob
  integer(i4)                                  :: ns
  integer(i4)                                  :: nss
  integer(i4)                                  :: status
  real(dp)                                     :: qj
  real(dp),   dimension(:), allocatable        :: dtmp
#ifdef TRACE
  call trace_in('setdr2mat')
#endif
!****************************************
!  Use Bulk Lattice second derivatives  *
!****************************************
!
!  Restore from disk
!
  if (lrest) then
    if (ioproc) then
      open(44,form='unformatted',status='old',err=10)
      goto 20
10    call outerror('restart file fort.44 not found',0_i4)
      call stopnow('setdr2mat')
20    rewind(44)
      j = 1
      do i = 1,3
        read(44,err=30)(derv3(j,i),j=1,3*numat)
      enddo
      do i = 1,6
        read(44,err=30)(sderv2(j,i),j=1,6)
      enddo
      do i = 1,3
        read(44,err=30)(diconh(j,i),j=1,3)
      enddo
      do i = 1,3
        read(44,err=30)(dicons(j,i),j=1,3)
      enddo
    endif
    if (nprocs.gt.1) then
      call sendall2D(derv3,maxd2,3_i4*numat,3_i4,ioproc,"setdr2mat","derv3")
      call sendall(sderv2,36_i4,ioproc,"setdr2mat","sderv2")
      call sendall(diconh,9_i4,ioproc,"setdr2mat","sderv2")
      call sendall(dicons,9_i4,ioproc,"setdr2mat","sderv2")
    endif
  else
!
!  Second derivative matrix will have already been inverted during call to property.
!
    if (lshello) then
      ns = 0
      derv3(1:3*numat,1:3) = 0.0_dp
      do is = 1,nshellonnode
        ind = 3*(node2atom(nshonnodeptr(is))-1)
        ns = ns + 1
        inds = 3*(ns-1)
        nss = 0
        do js = 1,nshell
          j = nshptr(js)
          nss = nss + 1
          qj = qf(j)*occuf(j)
          indss = 3*(nss - 1)
          do ii = 1,3
            do jj = 1,3
              derv3(ind+jj,ii) = derv3(ind+jj,ii) + qj*dervi(indss+jj,inds+ii)
            enddo
          enddo
        enddo
      enddo
      if (nprocs.gt.1) then
        allocate(dtmp(3*numat),stat=status)
        if (status/=0) call outofmemory('setdr2mat','dtmp')
!
!  Globalise derv3
!
        call sumall(derv3(1,1),dtmp,3_i4*numat,"setdr2mat","derv3")
        derv3(1:3*numat,1) = dtmp(1:3*numat)
        call sumall(derv3(1,2),dtmp,3_i4*numat,"setdr2mat","derv3")
        derv3(1:3*numat,2) = dtmp(1:3*numat)
        call sumall(derv3(1,3),dtmp,3_i4*numat,"setdr2mat","derv3")
        derv3(1:3*numat,3) = dtmp(1:3*numat)
!
        deallocate(dtmp,stat=status)
        if (status/=0) call deallocate_error('setdr2mat','dtmp')
      endif
    else
!
!  Need to zero elements of dervi surrounding matrix to account
!  for atom that was excluded during inversion, otherwise radial
!  terms can get mixed in.
!
      if (nprocs.gt.1) then
!
!  In parallel, partial occupancy is not currently allowed
!
        if (atom2local(numat).gt.0) then
          indx = 3*(natomsonnode-1) + 1
          indy = indx + 1
          indz = indy + 1
          do i = 1,3*numat
            dervi(i,indx) = 0.0_dp
            dervi(i,indy) = 0.0_dp
            dervi(i,indz) = 0.0_dp
          enddo
        endif
!
        indx = 3*(numat-1) + 1
        indy = indx + 1
        indz = indy + 1
        if (atom2local(numat).gt.0) then
          do i = 1,3*(natomsonnode-1)
            dervi(indx,i) = 0.0_dp
            dervi(indy,i) = 0.0_dp
            dervi(indz,i) = 0.0_dp
          enddo
        else
          do i = 1,3*natomsonnode
            dervi(indx,i) = 0.0_dp
            dervi(indy,i) = 0.0_dp
            dervi(indz,i) = 0.0_dp
          enddo
        endif
!
        derv3(1:3*numat,1:3) = 0.0_dp
!
        do i = 1,natomsonnode
          ind = 3*(i - 1)
          indf = 3*(node2atom(i) - 1)
          do j = 1,numat
            qj = qf(j)*occuf(j)
            indj = 3*(j - 1)
            do ii = 1,3
              do jj = 1,3
                derv3(indf+jj,ii) = derv3(indf+jj,ii) + qj*dervi(indj+jj,ind+ii)
              enddo
            enddo
          enddo
        enddo
!
        allocate(dtmp(3*numat),stat=status)
        if (status/=0) call outofmemory('setdr2mat','dtmp')
!
!  Globalise derv3
!
        call sumall(derv3(1,1),dtmp,3_i4*numat,"setdr2mat","derv3")
        derv3(1:3*numat,1) = dtmp(1:3*numat)
        call sumall(derv3(1,2),dtmp,3_i4*numat,"setdr2mat","derv3")
        derv3(1:3*numat,2) = dtmp(1:3*numat)
        call sumall(derv3(1,3),dtmp,3_i4*numat,"setdr2mat","derv3")
        derv3(1:3*numat,3) = dtmp(1:3*numat)
!
        deallocate(dtmp,stat=status)
        if (status/=0) call deallocate_error('setdr2mat','dtmp')
      else
        indx = 3*(ncsfoc-1) + 1
        indy = indx + 1
        indz = indy + 1
        do i = 1,3*ncsfoc
          dervi(i,indx) = 0.0_dp
          dervi(i,indy) = 0.0_dp
          dervi(i,indz) = 0.0_dp
        enddo
        do i = 1,3*(ncsfoc-1)
          dervi(indx,i) = 0.0_dp
          dervi(indy,i) = 0.0_dp
          dervi(indz,i) = 0.0_dp
        enddo
        do i = 1,numat
          ind = 3*(i - 1)
          indi = 3*(iocptr(i) - 1)
          do j = 1,3
            do k = 1,3
              derv3(ind+k,j) = 0.0_dp
            enddo
          enddo
          do j = 1,numat
            qj = qf(j)*occuf(j)
            indj = 3*(iocptr(j) - 1)
            do ii = 1,3
              do jj = 1,3
                derv3(ind+jj,ii) = derv3(ind+jj,ii) + qj*dervi(indj+jj,indi+ii)
              enddo
            enddo
          enddo
        enddo
      endif
    endif
    if (lsave.and.ioproc) then
      if (.not.lrest) open(44,form='unformatted',status='unknown')
      rewind(44)
      do i = 1,3
        write(44)(derv3(j,i),j=1,3*numat)
      enddo
      do i = 1,6
        write(44)(sderv2(j,i),j=1,6)
      enddo
      do i = 1,3
        write(44)(diconh(j,i),j=1,3)
      enddo
      do i = 1,3
        write(44)(dicons(j,i),j=1,3)
      enddo
    endif
  endif
#ifdef TRACE
  call trace_out('setdr2mat')
#endif
  return
30 call outerror('restart file fort.44 is of an incorrect size',0_i4)
  naprob = j - 1 + (i - 1)*3*numat
  if (ioproc) then
    write(ioout,'(''**** Probable system size of fort.44 = '',i5,'' atoms        ****'',/)') naprob
  endif
#ifdef TRACE
  call trace_out('setdr2mat')
#endif
  call stopnow('setdr2mat')
!
  end
