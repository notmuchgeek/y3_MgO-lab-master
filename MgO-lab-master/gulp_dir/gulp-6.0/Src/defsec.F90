  subroutine defsec
!
!  Generates variable reduced second derivative matrix
!
!  Can work in one of two modes:
!
!  (a) Reduce symmetry generated derv2(3*nreg1,3*ndasym) 
!  (b) Full nosymmetry generated derv2(3*nreg1,3*nreg1)
!
!  Both reduced to derv2(nvar,nvar)
!
!  11/07 Unused variables removed
!   5/17 Parallel second derivative modifications made
!   2/18 Trace added
!   3/18 Parallel I/O corrected
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
  use control
  use current
  use defects
  use derivatives
  use general
  use iochannels
  use optimisation
  use parallel
  use times
#ifdef TRACE
  use trace,       only : trace_in, trace_out
#endif
  use transform
  implicit none
#ifdef MPI
  include 'mpif.h'
#endif
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ia
  integer(i4)                                  :: ic
  integer(i4)                                  :: id
  integer(i4)                                  :: idloc
#ifdef MPI
  integer(i4)                                  :: ii
#endif
  integer(i4)                                  :: iloc
  integer(i4)                                  :: ind
  integer(i4)                                  :: ix
  integer(i4)                                  :: iy
  integer(i4)                                  :: iz
  integer(i4)                                  :: j
  integer(i4)                                  :: jj
  integer(i4)                                  :: jjx
  integer(i4)                                  :: jk
  integer(i4)                                  :: jx
  integer(i4)                                  :: jy
  integer(i4)                                  :: jz
  integer(i4)                                  :: jd
  integer(i4)                                  :: k
  integer(i4)                                  :: maxlim
  integer(i4)                                  :: maxlimloc
  integer(i4)                                  :: maxlims
  integer(i4)                                  :: neqj
  integer(i4)                                  :: status
  logical                                      :: lredlocal
  real(dp),    dimension(:), allocatable       :: tmpx
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: d2(3,3)
  real(dp)                                     :: t1
  real(dp)                                     :: t2
#ifdef MPI
  integer                                      :: MPIerror
  integer(i4)                                  :: inodec
  integer(i4)                                  :: inodev
  integer(i4)                                  :: nlocalcomm
  integer(i4)                                  :: ncopy
  integer(i4)                                  :: nrecv
  integer(i4)                                  :: nsend
  integer                                      :: ntag
  integer                                      :: nnode
  integer                                      :: ntmp
  integer(i4), dimension(:),   allocatable     :: nrecvnode     ! Pointer to node to receive from
  integer(i4), dimension(:),   allocatable     :: nsendnode     ! Pointer to node to send to
  integer(i4), dimension(:,:), allocatable     :: ncopyptr      ! Pointer information for data to be copied
  integer(i4), dimension(:,:), allocatable     :: nrecvptr      ! Pointer information for data to be received
  integer(i4), dimension(:,:), allocatable     :: nsendptr      ! Pointer information for data to be sent
  integer(i4), dimension(:),   allocatable     :: Request       ! Array for requests to MPI
  integer(i4), dimension(:,:), allocatable     :: StatMPI       ! Array for status from MPI
  real(dp),    dimension(:,:), allocatable     :: datarecv
  real(dp),    dimension(:,:), allocatable     :: datasend
#endif
#ifdef TRACE
  call trace_in('defsec')
#endif
!
  t1 = g_cpu_time()
  maxlim = 3*nreg1
  maxlimloc = 3*nreg1onnode
  if (ldbsm) then
    maxlim = maxlim + nreg1
    maxlimloc = maxlimloc + nreg1onnode
  endif
  if (index(keyword,'derv2').ne.0) then
    if (ioproc) then
      write(ioout,'(/,''  Second derivative matrix before reduction : (eV/Angstrom**2)'',/)')
    endif
    call mpbarrier
    if (ld2sym) then
      maxlims = 3*ndasym
      if (ldbsm) maxlims = maxlims + ndasym
      do i = 1,maxlim + 3
        write(ioout,'(9f10.4)')(derv2(i,j),j=1,maxlims)
      enddo
    else
#ifdef MPI
      if (lioproconly) then
!
!  Allocate temporary workspace for communication
!
        ntmp = 3*(maxlim+3)
        ntag = 1
        allocate(datasend(maxlim+3,3_i4),stat=status)
        if (status/=0) call outofmemory('defsec','dtmp')
        allocate(Request(1),stat=status)
        if (status/=0) call outofmemory('defsec','Request')
        allocate(StatMPI(MPI_Status_Size,1),stat=status)
        if (status/=0) call outofmemory('defsec','StatMPI')
!
        do i = 1,nreg1
          iloc = reg12local(i)
          if (reg12node(i).gt.0) then
!
!  Post receive
!
            if (ioproc) then
              nnode = reg12node(i)
              call MPI_IRecv(datasend,ntmp,MPI_double_precision,nnode, &
                             ntag,MPI_Comm_World,Request,MPIerror)
            endif
!
!  Pass data to ioproc for writing
!
            if (iloc.gt.0) then
              ind = 3*(iloc-1)
              do ii = 1,3
                do j = 1,maxlim+3
                  datasend(j,ii) = derv2(j,ind+ii)
                enddo
              enddo
!
!  Post send
!
              call MPI_ISend(datasend,ntmp,MPI_double_precision,0, &
                             ntag,MPI_Comm_World,Request,MPIerror)
            endif
            if (ioproc.or.iloc.gt.0) then
              call MPI_WaitAll(1,Request,StatMPI,MPIerror)
            endif
            if (ioproc) then
!
!  Write on I/O node
!
              do ii = 1,3
                write(ioout,'(9f10.4)')(datasend(j,ii),j=1,maxlim+3)
              enddo
            endif
          else
            if (iloc.gt.0) then
              ind = 3*(iloc-1) 
              write(ioout,'(9f10.4)')(derv2(j,ind+1),j=1,maxlim+3)
              write(ioout,'(9f10.4)')(derv2(j,ind+2),j=1,maxlim+3)
              write(ioout,'(9f10.4)')(derv2(j,ind+3),j=1,maxlim+3)
            endif
          endif
        enddo
!
        deallocate(StatMPI,stat=status)
        if (status/=0) call deallocate_error('defsec','StatMPI')
        deallocate(Request,stat=status)
        if (status/=0) call deallocate_error('defsec','Request')
        deallocate(datasend,stat=status)
        if (status/=0) call deallocate_error('defsec','datasend')
      else
#endif
        do i = 1,nreg1
          iloc = reg12local(i)
          if (iloc.gt.0) then
            ind = 3*(iloc-1) 
            write(ioout,'(9f10.4)')(derv2(j,ind+1),j=1,maxlim+3)
            write(ioout,'(9f10.4)')(derv2(j,ind+2),j=1,maxlim+3)
            write(ioout,'(9f10.4)')(derv2(j,ind+3),j=1,maxlim+3)
          endif
          call mpbarrier
        enddo
#ifdef MPI
      endif
#endif
    endif
    if (ioproc) then
      write(ioout,'(/)')
    endif
    call mpbarrier
  endif
!
!  If no reduction is needed then return
!
  if (nvar.eq.maxlim) goto 10
!
!  Check size of second derivative arrays
!
  if (maxlimloc+nvaronnode.gt.maxd2u) then
    maxd2u = maxlimloc + nvaronnode
    call changemaxd2
  endif
  if (maxlim.gt.maxd2) then
    maxd2 = maxlim
    call changemaxd2
  endif
!
  if (.not.ldsym) then
!************************************
!  Reduce second derivative matrix  *
!************************************
!
!  Is reduction purely local?
!
    lredlocal = .true.
    if (nprocs.gt.1) then
      i = 0
      do while (lredlocal.and.i.lt.nvar)
        i = i + 1
        ia = idoptindex(i)
        if (nvar2local(i).ne.reg12local(ia)) lredlocal = .false.
      enddo
    endif
!
    if (lredlocal) then
!
!  Local reduction only
!
      if (nprocs.gt.1) then
        do iloc = 1,nvaronnode
          ia = idoptindex(node2var(iloc))
          if (idopttype(node2var(iloc)).eq.idopt_dx) then
            ic = 1
          elseif (idopttype(node2var(iloc)).eq.idopt_dy) then
            ic = 2
          elseif (idopttype(node2var(iloc)).eq.idopt_dz) then
            ic = 3
          endif
          idloc = 3*(reg12local(ia)-1) + ic
          do j = 1,nvar
            if (idopttype(j).eq.idopt_dx) then
              jd = 3*idoptindex(j) - 2
            elseif (idopttype(j).eq.idopt_dy) then
              jd = 3*idoptindex(j) - 1
            elseif (idopttype(j).eq.idopt_dz) then
              jd = 3*idoptindex(j)
            elseif (idopttype(j).eq.idopt_dradius) then
              jd = 3*nreg1 + idoptindex(j)
            endif
            derv2(j,iloc) = derv2(jd,idloc)
          enddo
        enddo
      else
        do i = 1,nvar
          if (idopttype(i).eq.idopt_dx) then
            id = 3*idoptindex(i) - 2
          elseif (idopttype(i).eq.idopt_dy) then
            id = 3*idoptindex(i) - 1
          elseif (idopttype(i).eq.idopt_dz) then
            id = 3*idoptindex(i)
          elseif (idopttype(i).eq.idopt_dradius) then
            id = 3*nreg1 + idoptindex(i)
          endif
          do j = 1,nvar
            if (idopttype(j).eq.idopt_dx) then
              jd = 3*idoptindex(j) - 2
            elseif (idopttype(j).eq.idopt_dy) then
              jd = 3*idoptindex(j) - 1
            elseif (idopttype(j).eq.idopt_dz) then
              jd = 3*idoptindex(j)
            elseif (idopttype(j).eq.idopt_dradius) then
              jd = 3*nreg1 + idoptindex(j)
            endif
            derv2(j,i) = derv2(jd,id)
          enddo
        enddo
      endif
#ifdef MPI
    else
!
!  Parallel reduction
!
!  Compress derivatives within columns to reduce amount of data that has to be transfered
!
      do i = 1,3*nreg1onnode
        do j = 1,nvar
          jd = idoptindex(j)
          if (idopttype(j).eq.idopt_dx) then
            jd = 3*idoptindex(j) - 2
          elseif (idopttype(j).eq.idopt_dy) then
            jd = 3*idoptindex(j) - 1
          elseif (idopttype(j).eq.idopt_dz) then
            jd = 3*idoptindex(j)
          elseif (idopttype(j).eq.idopt_dradius) then
            jd = 3*nreg1 + idoptindex(j)
          endif
          derv2(j,i) = derv2(jd,i)
        enddo
      enddo
!
!  Allocate arrays for keeping track of data to be moved
!
      allocate(ncopyptr(2_i4,nvaronnode),stat=status)
      if (status/=0) call outofmemory('defsec','ncopyptr')
      allocate(nrecvnode(nvaronnode),stat=status)
      if (status/=0) call outofmemory('defsec','nrecvnode')
      allocate(nrecvptr(2_i4,nvaronnode),stat=status)
      if (status/=0) call outofmemory('defsec','nrecvptr')
      allocate(nsendnode(3_i4*nreg1onnode),stat=status)
      if (status/=0) call outofmemory('defsec','nsendnode')
      allocate(nsendptr(2_i4,3_i4*nreg1onnode),stat=status)
      if (status/=0) call outofmemory('defsec','nsendptr')
!
!  Find out how many columns have to be sent and received and where they are going to or coming from
!
      ncopy = 0
      nrecv = 0
      nsend = 0
      do i = 1,nvar
!
!  inodev is node for variable
!
        inodev = nvar2node(i)
!
!  inodec is the node for the internal derivatives in derv2
!
        ia = idoptindex(i)
        if (idopttype(i).eq.idopt_dx) then
          ic = 1
        elseif (idopttype(i).eq.idopt_dy) then
          ic = 2
        elseif (idopttype(i).eq.idopt_dz) then
          ic = 3
        endif
        inodec = reg12node(ia)
!
!  convert ii to local version
!
        ii = 3*(reg12local(ia) - 1) + ic
!
        if (inodev.eq.procid.and.inodec.eq.procid) then
!
!  Initial and final location of column are on the same node => copy
!
          ncopy = ncopy + 1
          ncopyptr(1,ncopy) = nvar2local(i)
          ncopyptr(2,ncopy) = ii
        elseif (inodec.eq.procid) then
!
!  Initial location of column is on this node => send
!
          nsend = nsend + 1
          nsendnode(nsend) = inodev
          nsendptr(1,nsend) = ii
          nsendptr(2,nsend) = i
        elseif (inodev.eq.procid) then
!
!  Final location of column is on this node => receive
!
          nrecv = nrecv + 1
          nrecvnode(nrecv) = inodec
          nrecvptr(1,nrecv) = nvar2local(i)
          nrecvptr(2,nrecv) = i
        endif
      enddo
!
!  Allocate storage for sending / receiving
!
      allocate(datarecv(nvar,nrecv),stat=status)
      if (status/=0) call outofmemory('defsec','datarecv')
      allocate(datasend(nvar,nsend),stat=status)
      if (status/=0) call outofmemory('defsec','datasend')
      allocate(Request(nrecv+nsend),stat=status)
      if (status/=0) call outofmemory('defsec','Request')
      allocate(StatMPI(MPI_Status_Size,nrecv+nsend),stat=status)
      if (status/=0) call outofmemory('defsec','StatMPI')
!
!  Fill arrays for data sending
!
      do i = 1,nsend
        ii = nsendptr(1,i)
        datasend(1:nvar,i) = derv2(1:nvar,ii)
      enddo
!
!  Send data and post receives
!
      nlocalcomm = 0
      do i = 1,nrecv
        nlocalcomm = nlocalcomm + 1
        call MPI_IRecv(datarecv(1,i),nvar,MPI_double_precision,nrecvnode(i), &
                       nrecvptr(2,i),MPI_Comm_World,Request(nlocalcomm),MPIerror)
      enddo
      do i = 1,nsend
        nlocalcomm = nlocalcomm + 1
        call MPI_ISend(datasend(1,i),nvar,MPI_double_precision,nsendnode(i), &
                       nsendptr(2,i),MPI_Comm_World,Request(nlocalcomm),MPIerror)
      enddo
!
!  Copy data that is local into the correct location
!
      do i = 1,ncopy
        do j = 1,nvar
          derv2(j,ncopyptr(1,i)) = derv2(j,ncopyptr(2,i))
        enddo
      enddo
!
!  Wait for data transfer to finish
!
      call MPI_WaitAll(nlocalcomm,Request,StatMPI,MPIerror)
!
!  Place received data into the right locations
!
      do i = 1,nrecv
        ii = nrecvptr(1,i)
        derv2(1:nvar,ii) = datarecv(1:nvar,i)
      enddo
!
!  Deallocate storage for sending
!
      deallocate(StatMPI,stat=status)
      if (status/=0) call deallocate_error('defsec','StatMPI')
      deallocate(Request,stat=status)
      if (status/=0) call deallocate_error('defsec','Request')
      deallocate(datasend,stat=status)
      if (status/=0) call deallocate_error('defsec','datasend')
      deallocate(datarecv,stat=status)
      if (status/=0) call deallocate_error('defsec','datarecv')
      deallocate(nsendptr,stat=status)
      if (status/=0) call deallocate_error('defsec','nsendptr')
      deallocate(nsendnode,stat=status)
      if (status/=0) call deallocate_error('defsec','nsendnode')
      deallocate(nrecvptr,stat=status)
      if (status/=0) call deallocate_error('defsec','nrecvptr')
      deallocate(nrecvnode,stat=status)
      if (status/=0) call deallocate_error('defsec','nrecvnode')
      deallocate(ncopyptr,stat=status)
      if (status/=0) call deallocate_error('defsec','ncopyptr')
#endif
    endif
  else
    allocate(tmpx(maxlim),stat=status)
    if (status/=0) call outofmemory('defsec','tmpx')
    if (ld2sym) then
!****************************
!  Symmetry adapted method  *
!****************************
!
!  Reduce matrix to symmetry reduced full second derivative form
!
      do i = 1,ndasym
        ix = 3*(i - 1) + 1
        iy = ix + 1
        iz = ix + 2
        do j = 1,ndasym
          jj = ndsptr(j)
          jx = 3*(j - 1) + 1
          jy = jx + 1
          jz = jx + 2
          neqj = ndeqv(j)
          jjx = 3*(jj - 1)
          d2(1,1) = 0.0_dp
          d2(2,1) = 0.0_dp
          d2(3,1) = 0.0_dp
          d2(1,2) = 0.0_dp
          d2(2,2) = 0.0_dp
          d2(3,2) = 0.0_dp
          d2(1,3) = 0.0_dp
          d2(2,3) = 0.0_dp
          d2(3,3) = 0.0_dp
          do jk = 1,3*neqj
            d2(1,1) = d2(1,1) + derv2(jjx + jk,ix)*tmat(jjx + jk,ix)
            d2(2,1) = d2(2,1) + derv2(jjx + jk,iy)*tmat(jjx + jk,ix)
            d2(3,1) = d2(3,1) + derv2(jjx + jk,iz)*tmat(jjx + jk,ix)
            d2(1,2) = d2(1,2) + derv2(jjx + jk,ix)*tmat(jjx + jk,iy)
            d2(2,2) = d2(2,2) + derv2(jjx + jk,iy)*tmat(jjx + jk,iy)
            d2(3,2) = d2(3,2) + derv2(jjx + jk,iz)*tmat(jjx + jk,iy)
            d2(1,3) = d2(1,3) + derv2(jjx + jk,ix)*tmat(jjx + jk,iz)
            d2(2,3) = d2(2,3) + derv2(jjx + jk,iy)*tmat(jjx + jk,iz)
            d2(3,3) = d2(3,3) + derv2(jjx + jk,iz)*tmat(jjx + jk,iz)
          enddo
          derv2(jx,ix) = d2(1,1)
          derv2(jy,ix) = d2(1,2)
          derv2(jz,ix) = d2(1,3)
          derv2(jx,iy) = d2(2,1)
          derv2(jy,iy) = d2(2,2)
          derv2(jz,iy) = d2(2,3)
          derv2(jx,iz) = d2(3,1)
          derv2(jy,iz) = d2(3,2)
          derv2(jz,iz) = d2(3,3)
        enddo
      enddo
!
!  Symmetrise matrix
!
      maxlim = 3*ndasym
!
!  Reduce to nvar x nvar form
!
      if (ndcon.eq.0) then
        do i = 1,nvar
          id = idoptindex(i)
          if (idopttype(i).eq.idopt_dx) then
            id = 3*idoptindex(i) - 2
          elseif (idopttype(i).eq.idopt_dy) then
            id = 3*idoptindex(i) - 1
          elseif (idopttype(i).eq.idopt_dz) then
            id = 3*idoptindex(i)
          elseif (idopttype(i).eq.idopt_dradius) then
            id = 3*ndasym + idoptindex(i)
          endif
          do j = 1,nvar
            jd = idoptindex(j)
            if (idopttype(j).eq.idopt_dx) then
              jd = 3*idoptindex(j) - 2
            elseif (idopttype(j).eq.idopt_dy) then
              jd = 3*idoptindex(j) - 1
            elseif (idopttype(j).eq.idopt_dz) then
              jd = 3*idoptindex(j)
            elseif (idopttype(j).eq.idopt_dradius) then
              jd = 3*ndasym + idoptindex(j)
            endif
            derv2(j,i) = derv2(jd,id)
          enddo
        enddo
      else
        if ((nvar + maxlim).le.maxd2) then
          do i = 1,nvar
            do j = 1,maxlim
              derv2(j,maxlim + i) = 0.0_dp
              do k = 1,maxlim
                derv2(j,maxlim + i) = derv2(j,maxlim + i) + tmat(k,maxlim + i)*derv2(k,j)
              enddo
            enddo
          enddo
          do i = 1,nvar
            do j = 1,nvar
              derv2(j,i) = 0.0_dp
              do k = 1,maxlim
                derv2(j,i) = derv2(j,i) + tmat(k,maxlim + j)*derv2(k,maxlim + i)
              enddo
            enddo
          enddo
        else
          do i = 1,nvar
            do j = 1,maxlim
              tmpx(j) = 0.0_dp
              do k = 1,maxlim
                tmpx(j) = tmpx(j) + tmat(k,maxlim + i)*derv2(j,k)
              enddo
            enddo
            do j = 1,nvar
              derv2(j,i) = 0.0_dp
              do k = 1,maxlim
                derv2(j,i) = derv2(j,i) + tmpx(k)*tmat(k,maxlim + j)
              enddo
            enddo
          enddo
        endif
      endif
    else
!***********************
!  No symmetry method  *
!***********************
      if ((2*nvar).le.maxd2) then
        do i = 1,nvar
          do j = 1,maxlim
            tmat(j,nvar + i) = 0.0_dp
            do k = 1,maxlim
              tmat(j,nvar + i) = tmat(j,nvar + i) + tmat(k,i)*derv2(k,j)
            enddo
          enddo
        enddo
        do i = 1,nvar
          do j = 1,nvar
            derv2(j,i) = 0.0_dp
            do k = 1,maxlim
              derv2(j,i) = derv2(j,i) + tmat(k,j)*tmat(k,nvar + i)
            enddo
          enddo
        enddo
      else
        do i = 1,nvar
          do j = 1,maxlim
            tmpx(j) = 0.0_dp
            do k = 1,maxlim
              tmpx(j) = tmpx(j) + tmat(k,i)*derv2(j,k)
            enddo
          enddo
          do j = 1,nvar
            derv2(j,i) = 0.0_dp
            do k = 1,maxlim
              derv2(j,i) = derv2(j,i) + tmpx(k)*tmat(k,j)
            enddo
          enddo
        enddo
      endif
    endif
    deallocate(tmpx,stat=status)
    if (status/=0) call deallocate_error('defsec','tmpx')
  endif
10 if (index(keyword,'derv2').ne.0) then
    if (ioproc) then
      write(ioout,'(/,''  Reduced second derivative matrix : (eV/Angstrom**2)'',/)')
    endif
#ifdef MPI
    if (lioproconly) then
!
!  Allocate temporary workspace for communication
!
      ntmp = nvar
      ntag = 1
      allocate(tmpx(nvar),stat=status)
      if (status/=0) call outofmemory('defsec','temp')
      allocate(Request(1),stat=status)
      if (status/=0) call outofmemory('defsec','Request')
      allocate(StatMPI(MPI_Status_Size,1),stat=status)
      if (status/=0) call outofmemory('defsec','StatMPI')
!
      do i = 1,nvar
        iloc = nvar2local(i)
        if (nvar2node(i).ne.0_i4) then
!
!  Post receive
!
          if (ioproc) then
            nnode = nvar2node(i)
            call MPI_IRecv(tmpx,ntmp,MPI_double_precision,nnode, &
                           ntag,MPI_Comm_World,Request,MPIerror)
          endif
!
!  Pass data to ioproc for writing
!
          if (iloc.gt.0) then
            tmpx(1:nvar) = derv2(1:nvar,iloc)
!
!  Post send
!
            call MPI_ISend(tmpx,ntmp,MPI_double_precision,0, &
                           ntag,MPI_Comm_World,Request,MPIerror)
          endif
          if (ioproc.or.iloc.gt.0) then
            call MPI_WaitAll(1,Request,StatMPI,MPIerror)
          endif
          if (ioproc) then
!
!  Write on I/O node
!
            write(ioout,'(9f10.4)')(tmpx(j),j=1,nvar)
          endif
        else
          if (iloc.gt.0) then
            write(ioout,'(9f10.4)')(derv2(j,iloc),j = 1,nvar)
          endif
        endif
      enddo
!
      deallocate(StatMPI,stat=status)
      if (status/=0) call deallocate_error('defsec','StatMPI')
      deallocate(Request,stat=status)
      if (status/=0) call deallocate_error('defsec','Request')
      deallocate(tmpx,stat=status)
      if (status/=0) call deallocate_error('defsec','tmpx')
    else
#endif
      call mpbarrier
      do i = 1,nvar
        iloc = nvar2local(i)
        if (iloc.gt.0) then
          write(ioout,'(9f10.4)')(derv2(j,iloc),j = 1,nvar)
        endif
        call mpbarrier
      enddo
#ifdef MPI
    endif
#endif
    if (ioproc) then
      write(ioout,'(/)')
    endif
    call mpbarrier
  endif
!
  t2 = g_cpu_time()
  tdel = t2 - t1
  thes = thes + tdel
#ifdef TRACE
  call trace_out('defsec')
#endif
!
  return
  end
