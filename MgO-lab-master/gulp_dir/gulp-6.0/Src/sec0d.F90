  subroutine sec0d
!
!  Generate symmetry adapted cluster second derivative matrix.
!
!  Distributed memory version.
!
!  NB: Algorithm needs nblocksizevar to be 3 x nblocksize
!
!   1/17 Created from sec0
!   3/17 fix_atom option added
!   4/17 Parallelisation now working without tmat for general fixed atom
!   4/17 Numat-numat tmat algorithm now implemented in parallel
!   7/17 ld for idest corrected to be consistent with tmat
!   8/17 Routine wrapped in ifdef since it won't work or be called without MPI
!   2/18 Trace added
!   3/18 Parallel I/O corrected
!   3/19 iopt replaced by ioptindex and iopttype
!   4/19 Format for symmetrised second derivatives updated
!   5/20 Rigid molecule modifications added
!   6/20 Corrections to non-tmat algorithm 
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
!  Julian Gale, CIC, Curtin University, June 2020
!
  use control
  use current
  use derivatives
  use iochannels
  use molecule
#ifdef MPI
  use optimisation
#endif
  use parallel
  use times
#ifdef TRACE
  use trace,        only : trace_in, trace_out
#endif
  use transform
  implicit none
#ifdef MPI
  include 'mpif.h'
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ii
  integer(i4)                                  :: iloc
  integer(i4)                                  :: indj
  integer(i4)                                  :: indk
  integer(i4)                                  :: ino
  integer(i4)                                  :: iproc
  integer(i4)                                  :: ix
  integer(i4)                                  :: j
  integer(i4)                                  :: jj
  integer(i4)                                  :: jx
  integer(i4)                                  :: k
  integer(i4)                                  :: kk
  integer(i4)                                  :: n3f
  integer(i4)                                  :: n3floc
  integer(i4)                                  :: ni
  integer(i4)                                  :: nj
  integer(i4)                                  :: nmf
  integer(i4)                                  :: n3mf
  integer(i4)                                  :: status
!
  integer                                      :: ides2(9)
  integer                                      :: idesd(9)
  integer                                      :: idest(9)
  integer                                      :: ifails
  integer                                      :: ld
  integer                                      :: nb
  integer                                      :: ndiml
  integer                                      :: ndimr
  integer                                      :: ntag
  integer                                      :: nnode
  integer                                      :: ntmp
!
  integer                                      :: MPIerror
  integer(i4)                                  :: ia
  integer(i4)                                  :: ic
  integer(i4)                                  :: inodec
  integer(i4)                                  :: inodev
  integer(i4)                                  :: nlocalcomm
  integer(i4)                                  :: ncopy
  integer(i4)                                  :: nrecv
  integer(i4)                                  :: nsend
  integer(i4), dimension(:),   allocatable     :: nrecvnode     ! Pointer to node to receive from
  integer(i4), dimension(:),   allocatable     :: nsendnode     ! Pointer to node to send to
  integer(i4), dimension(:,:), allocatable     :: ncopyptr      ! Pointer information for data to be copied
  integer(i4), dimension(:,:), allocatable     :: nrecvptr      ! Pointer information for data to be received
  integer(i4), dimension(:,:), allocatable     :: nsendptr      ! Pointer information for data to be sent
  integer(i4), dimension(:),   allocatable     :: Request       ! Array for requests to MPI
  integer(i4), dimension(:,:), allocatable     :: StatMPI       ! Array for status from MPI
!
  logical                                      :: lmolQi
  logical                                      :: lmolQj
  logical                                      :: lmolTi
  logical                                      :: lmolTj
  logical                                      :: lsdebug
  logical                                      :: ltmat
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: t1
  real(dp)                                     :: t2
  real(dp), dimension(:),   allocatable        :: tmp2
  real(dp), dimension(:,:), allocatable        :: tmp2D
  real(dp), dimension(:,:), allocatable        :: datarecv
  real(dp), dimension(:,:), allocatable        :: datasend
#ifdef TRACE
  call trace_in('sec0d')
#endif
!
  t1 = g_cpu_time()
  lsdebug = (index(keyword,'derv2').ne.0)
  n3f = 3*numat
  n3floc = 3*natomsonnode
  if (nbsm.gt.0) then
    n3f = n3f + numat
    n3floc = n3floc + natomsonnode
  endif
!
  if (lrigid) then
    nmf = nmol
    n3mf = 3*nmf
  endif
!
!  Work out whether full tmat multiplication is needed - use faster method to handle P1 fully relaxed case
!
  if (ncon.eq.0) then
    ltmat = .false.
  else
    ltmat = .true.
  endif
!
!  Print out rigid molecule second derivatives
!
  if (lsdebug.and.ioproc) then
    write(ioout,'(/)')
    if (lrigid) then
      write(ioout,'(/,''  Second Derivative Matrix : Rigid molecule translations :'',/)')
      do i = 1,nmol
        do j = 1,nmol
          write(ioout,'(2x,2i5,'' x'',1x,3(f9.4))') i,j,(molTTdrv(1,jj,i,j),jj=1,3)
          write(ioout,'(2x,2i5,'' y'',1x,3(f9.4))') i,j,(molTTdrv(2,jj,i,j),jj=1,3)
          write(ioout,'(2x,2i5,'' z'',1x,3(f9.4))') i,j,(molTTdrv(3,jj,i,j),jj=1,3)
        enddo
      enddo
      write(ioout,'(/,''  Second Derivative Matrix : Rigid molecule quaternions :'',/)')
      do i = 1,nmol
        do j = 1,nmol
          write(ioout,'(2x,2i5,'' x'',1x,3(f9.4))') i,j,(molQQdrv(1,jj,i,j),jj=1,3)
          write(ioout,'(2x,2i5,'' y'',1x,3(f9.4))') i,j,(molQQdrv(2,jj,i,j),jj=1,3)
          write(ioout,'(2x,2i5,'' z'',1x,3(f9.4))') i,j,(molQQdrv(3,jj,i,j),jj=1,3)
        enddo
      enddo
      write(ioout,'(/,''  Second Derivative Matrix : Rigid molecule quaternion-translation :'',/)')
      do i = 1,nmol
        do j = 1,nmol
          write(ioout,'(2x,2i5,'' x'',1x,3(f9.4))') i,j,(molQTdrv(1,jj,i,j),jj=1,3)
          write(ioout,'(2x,2i5,'' y'',1x,3(f9.4))') i,j,(molQTdrv(2,jj,i,j),jj=1,3)
          write(ioout,'(2x,2i5,'' z'',1x,3(f9.4))') i,j,(molQTdrv(3,jj,i,j),jj=1,3)
        enddo
      enddo
      write(ioout,'(/)')
    endif
  endif
!********************************************
!  Rigid molecules - check for linear case  *
!********************************************
  if (lrigid) then
    do i = 1,nmol
      do ix = 1,3
!
!  If molecule is linear or there are no rotational forces that set the diagonal element to be a dummy value
!  so that the Hessian can be inverted
!
        if (abs(molQQdrv(ix,ix,i,i)).lt.0.001_dp) then
          do jx = 1,3
            molQQdrv(jx,ix,i,i) = 0.0_dp
            molQQdrv(ix,jx,i,i) = 0.0_dp
          enddo
          molQQdrv(ix,ix,i,i) = 1.0_dp
        endif
      enddo
    enddo
  endif
!******************************************
!  Internal derivatives :                 *
!  Symmetry transform second derivatives  *
!******************************************
  if (ltmat) then
    allocate(tmp2D(maxn3f,ninternalatmonnode),stat=status)
    if (status/=0) call outofmemory('sec0d','tmp2D')
!
!  Set up Blacs descriptors for matrices
!
    nb = nblocksizevar
    ifails = 0
    ndiml = n3f
    ndimr = ninternalatm
    ld = n3f
    call descinit( ides2, ndiml, ndimr, nb, nb, 0, 0, iBlacsContext, ld, ifails )
!
    if (ifails.ne.0) then
      call outerror('initialisation in descinit failed',0_i4)
      call stopnow('sec0d')
    endif
!
    nb = nblocksize
    ndiml = n3f
    ndimr = n3f
    ld = maxd2
    call descinit( idesd, ndiml, ndimr, 3*nb, 3*nb, 0, 0, iBlacsContext, ld, ifails )
!
    if (ifails.ne.0) then
      call outerror('initialisation in descinit failed',0_i4)
      call stopnow('sec0d')
    endif
!
    nb = nblocksizevar
    ndiml = n3f
    ndimr = ninternalatm
    ld = maxn3f
    call descinit( idest, ndiml, ndimr, nb, nb, 0, 0, iBlacsContext, ld, ifails )
!
    if (ifails.ne.0) then
      call outerror('initialisation in descinit failed',0_i4)
      call stopnow('sec0d')
    endif
!
!  Perform parallel matrix-matrix multiply : D = T^t*D*T
!
    call pdgemm('n','n',n3f,ninternalatm,n3f,1.0d0,derv2,1,1,idesd,tmat,1,1,idest,0.0d0,tmp2D,1,1,ides2)
    call pdgemm('t','n',ninternalatm,ninternalatm,n3f,1.0d0,tmat,1,1,idest,tmp2D,1,1,ides2,0.0d0,derv2,1,1,idesd)
!
    deallocate(tmp2D,stat=status)
    if (status/=0) call deallocate_error('sec0d','tmp2D')
!
    if (lrigid) then
!------------------------------------------
!  Rigid molecule - rigid molecule terms  |
!------------------------------------------
      allocate(tmp2(2*n3mf),stat=status)
      if (status/=0) call outofmemory('sec0d','tmp2')
!
!  Translation - translation
!
      ii = ninternalatmonnode
      do i = 1,ninternalmolT
        if (nvar2node(ninternalatm+i).eq.procid) then
          ii = ii + 1
          tmp2(1:n3mf) = 0.0_dp
          do j = 1,nmf
            indj = 3*(j-1)
            do jj = 1,3
              do k = 1,nmf
                indk = 3*(k-1)
                do kk = 1,3
                  tmp2(indj+jj) = tmp2(indj+jj) + molTTdrv(jj,kk,j,k)*tmatT(indk+kk,i)
                enddo
              enddo
            enddo
          enddo
!
          do j = 1,ninternalmolT
            jj = ninternalatm + j
            derv2(jj,ii) = 0.0_dp
            do k = 1,nmf
              indk = 3*(k-1)
              do kk = 1,3
                derv2(jj,ii) = derv2(jj,ii) + tmp2(indk+kk)*tmatT(indk+kk,j)
              enddo
            enddo
          enddo
        endif
      enddo
!
!  Translation - quaternion
!
      ii = ninternalatmonnode
      do i = 1,ninternalmolT
        if (nvar2node(ninternalatm+i).eq.procid) then
          ii = ii + 1
          tmp2(1:n3mf) = 0.0_dp
          do j = 1,nmf
            indj = 3*(j-1)
            do jj = 1,3
              do k = 1,nmf
                indk = 3*(k-1)
                do kk = 1,3
                  tmp2(indj+jj) = tmp2(indj+jj) + molQTdrv(jj,kk,j,k)*tmatT(indk+kk,i)
                enddo
              enddo
            enddo
          enddo
!
          do j = 1,ninternalmolQ
            jj = ninternalatm + ninternalmolT + j
            derv2(jj,ii) = 0.0_dp
            do k = 1,nmf
              indk = 3*(k-1)
              do kk = 1,3
                derv2(jj,ii) = derv2(jj,ii) + tmp2(indk+kk)*tmatQ(indk+kk,j)
              enddo
            enddo
          enddo
        endif
      enddo
!
!  Quaternion - translation
!
      do i = 1,ninternalmolQ
        if (nvar2node(ninternalatm+ninternalmolT+i).eq.procid) then
          ii = nvar2local(ninternalatm+ninternalmolT+i)
          tmp2(1:n3mf) = 0.0_dp
          do j = 1,nmf
            indj = 3*(j-1)
            do jj = 1,3
              do k = 1,nmf
                indk = 3*(k-1)
                do kk = 1,3
                  tmp2(indj+jj) = tmp2(indj+jj) + molQTdrv(kk,jj,k,j)*tmatQ(indk+kk,i)
                enddo
              enddo
            enddo
          enddo
!
          do j = 1,ninternalmolT
            jj = ninternalatm + j
            derv2(jj,ii) = 0.0_dp
            do k = 1,nmf
              indk = 3*(k-1)
              do kk = 1,3
                derv2(jj,ii) = derv2(jj,ii) + tmp2(indk+kk)*tmatT(indk+kk,j)
              enddo
            enddo
          enddo
        endif
      enddo
!
!  Quaternion - quaternion
!
      do i = 1,ninternalmolQ
        if (nvar2node(ninternalatm+ninternalmolT+i).eq.procid) then
          ii = nvar2local(ninternalatm+ninternalmolT+i)
          tmp2(1:n3mf) = 0.0_dp
          do j = 1,nmf
            indj = 3*(j-1)
            do jj = 1,3
              do k = 1,nmf
                indk = 3*(k-1)
                do kk = 1,3
                  tmp2(indj+jj) = tmp2(indj+jj) + molQQdrv(jj,kk,j,k)*tmatQ(indk+kk,i)
                enddo
              enddo
            enddo
          enddo
!
          do j = 1,ninternalmolQ
            jj = ninternalatm + ninternalmolT + j
            derv2(jj,ii) = 0.0_dp
            do k = 1,nmf
              indk = 3*(k-1)
              do kk = 1,3
                derv2(jj,ii) = derv2(jj,ii) + tmp2(indk+kk)*tmatQ(indk+kk,j)
              enddo
            enddo
          enddo
        endif
      enddo
!--------------------------------
!  Rigid molecule - atom terms  |
!--------------------------------
      do i = 1,ninternalatm
        if (nvar2node(i).eq.procid) then
          ii = nvar2local(i)
          tmp2(1:2*n3mf) = 0.0_dp
          do j = 1,nmf
            indj = 3*(j-1)
            do jj = 1,3
              do k = 1,n3f
                tmp2(indj+jj)      = tmp2(indj+jj)      + molTCdrv(k,jj,j)*tmat(k,ii)
                tmp2(n3mf+indj+jj) = tmp2(n3mf+indj+jj) + molQCdrv(k,jj,j)*tmat(k,ii)
              enddo
            enddo
          enddo
        endif
!
!  Send to all nodes
!
        call sendall(tmp2,2_i4*n3mf,nvar2node(i),"sec0d","tmp2")
!
        if (nvar2node(i).eq.procid) then
          ii = nvar2local(i)
          do j = 1,ninternalmolT
            jj = ninternalatm + j
            derv2(jj,ii) = 0.0_dp
            do k = 1,nmf
              indk = 3*(k-1)
              do kk = 1,3
                derv2(jj,ii) = derv2(jj,ii) + tmp2(indk+kk)*tmatT(indk+kk,j)
              enddo
            enddo
          enddo
          do j = 1,ninternalmolQ
            jj = ninternalatm + ninternalmolT + j
            derv2(jj,ii) = 0.0_dp
            do k = 1,nmf
              indk = 3*(k-1)
              do kk = 1,3
                derv2(jj,ii) = derv2(jj,ii) + tmp2(n3mf+indk+kk)*tmatQ(indk+kk,j)
              enddo
            enddo
          enddo
        endif
!
        jj = ninternalatmonnode
        do j = 1,ninternalmolT
          if (nvar2node(ninternalatm+j).eq.procid) then
            jj = jj + 1
            derv2(i,jj) = 0.0_dp
            do k = 1,nmf
              indk = 3*(k-1)
              do kk = 1,3
                derv2(i,jj) = derv2(i,jj) + tmp2(indk+kk)*tmatT(indk+kk,j)
              enddo
            enddo
          endif
        enddo
        do j = 1,ninternalmolQ
          if (nvar2node(ninternalatm+ninternalmolT+j).eq.procid) then
            jj = nvar2local(ninternalatm+ninternalmolT+j)
            derv2(i,jj) = 0.0_dp
            do k = 1,nmf
              indk = 3*(k-1)
              do kk = 1,3
                derv2(i,jj) = derv2(i,jj) + tmp2(n3mf+indk+kk)*tmatQ(indk+kk,j)
              enddo
            enddo
          endif
        enddo
      enddo
!
      deallocate(tmp2,stat=status)
      if (status/=0) call deallocate_error('sec0d','tmp2')
    endif
  elseif (ninternal.eq.(n3f-3).and.nfixatom.eq.numat.and..not.lrigid) then
!
!  All atoms are variable apart from the last one and so there is nothing to do
!
    continue
  else
!
!  Compress derivatives within columns to reduce amount of data that has to be transfered
!
    do i = 1,3*natomsonnode
      do j = 1,ninternal
        nj = ioptindex(j)
        if (iopttype(j).eq.iopt_xf) then
          jj = 3*nasymnomolptr(nj) - 2
          derv2(j,i) = derv2(jj,i)
        elseif (iopttype(j).eq.iopt_yf) then
          jj = 3*nasymnomolptr(nj) - 1
          derv2(j,i) = derv2(jj,i)
        elseif (iopttype(j).eq.iopt_zf) then
          jj = 3*nasymnomolptr(nj)
          derv2(j,i) = derv2(jj,i)
        elseif (iopttype(j).eq.iopt_radius) then
          jj = 3*numat + nasymnomolptr(nj)
          derv2(j,i) = derv2(jj,i)
        endif
      enddo
    enddo
!
!  Allocate arrays for keeping track of data to be moved
!
    allocate(ncopyptr(2_i4,ninternalatmonnode),stat=status)
    if (status/=0) call outofmemory('sec0d','ncopyptr')
    allocate(nrecvnode(ninternalatmonnode),stat=status)
    if (status/=0) call outofmemory('sec0d','nrecvnode')
    allocate(nrecvptr(2_i4,ninternalatmonnode),stat=status)
    if (status/=0) call outofmemory('sec0d','nrecvptr')
    allocate(nsendnode(3_i4*natomsonnode),stat=status)
    if (status/=0) call outofmemory('sec0d','nsendnode')
    allocate(nsendptr(2_i4,3_i4*natomsonnode),stat=status)
    if (status/=0) call outofmemory('sec0d','nsendptr')
!
!  Find out how many columns have to be sent and received and where they are going to or coming from
!
    ncopy = 0
    nrecv = 0
    nsend = 0
    do i = 1,ninternalatm
!
!  inodev is node for variable
!
      inodev = nvar2node(i)
!
!  inodec is the node for the internal derivatives in derv2
!
      ia = nasymnomolptr(ioptindex(i))
      if (iopttype(i).eq.iopt_xf) then
        ic = 1
      elseif (iopttype(i).eq.iopt_yf) then
        ic = 2
      elseif (iopttype(i).eq.iopt_zf) then
        ic = 3
      endif
      inodec = atom2node(ia)
!
!  convert ii to local version
!
      ii = 3*(atom2local(ia) - 1) + ic
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
    allocate(datarecv(ninternalatm,nrecv),stat=status)
    if (status/=0) call outofmemory('sec0d','datarecv')
    allocate(datasend(ninternalatm,nsend),stat=status)
    if (status/=0) call outofmemory('sec0d','datasend')
    allocate(Request(nrecv+nsend),stat=status)
    if (status/=0) call outofmemory('sec0d','Request')
    allocate(StatMPI(MPI_Status_Size,nrecv+nsend),stat=status)
    if (status/=0) call outofmemory('sec0d','StatMPI')
!
!  Fill arrays for data sending
!
    do i = 1,nsend
      ii = nsendptr(1,i)
      datasend(1:ninternalatm,i) = derv2(1:ninternalatm,ii)
    enddo
!
!  Send data and post receives
!
    nlocalcomm = 0
    do i = 1,nrecv
      nlocalcomm = nlocalcomm + 1
      call MPI_IRecv(datarecv(1,i),ninternalatm,MPI_double_precision,nrecvnode(i), &
                     nrecvptr(2,i),MPI_Comm_World,Request(nlocalcomm),MPIerror)
    enddo
    do i = 1,nsend
      nlocalcomm = nlocalcomm + 1
      call MPI_ISend(datasend(1,i),ninternalatm,MPI_double_precision,nsendnode(i), &
                     nsendptr(2,i),MPI_Comm_World,Request(nlocalcomm),MPIerror)
    enddo
!
!  Copy data that is local into the correct location
!
    do i = 1,ncopy
      do j = 1,ninternalatm
        derv2(j,ncopyptr(1,i)) = derv2(j,ncopyptr(2,i))
      enddo
    enddo
!
!  While waiting for data to transfer, fill in molecule terms
!
    if (lrigid) then
      do ino = 1,ninternalonnode
        i = node2var(ino)
        lmolQi = .false.
        lmolTi = .false.
        ni = ioptindex(i)
        if (iopttype(i).eq.iopt_xf) then
          ii = 3*nasymnomolptr(ni) - 2
        elseif (iopttype(i).eq.iopt_yf) then
          ii = 3*nasymnomolptr(ni) - 1
        elseif (iopttype(i).eq.iopt_zf) then
          ii = 3*nasymnomolptr(ni)
        elseif (iopttype(i).eq.iopt_radius) then
          ii = 3*numat + nasymnomolptr(ni)
        elseif (iopttype(i).eq.iopt_xcom) then
          lmolTi = .true.
          ii = 1
        elseif (iopttype(i).eq.iopt_ycom) then
          lmolTi = .true.
          ii = 2
        elseif (iopttype(i).eq.iopt_zcom) then
          lmolTi = .true.
          ii = 3
        elseif (iopttype(i).eq.iopt_xqtn) then
          lmolQi = .true.
          ii = 1
        elseif (iopttype(i).eq.iopt_yqtn) then
          lmolQi = .true.
          ii = 2
        elseif (iopttype(i).eq.iopt_zqtn) then
          lmolQi = .true.
          ii = 3
        endif
        do j = 1,ninternal
          lmolQj = .false.
          lmolTj = .false.
          nj = ioptindex(j)
          if (iopttype(j).eq.iopt_xf) then
            jj = 3*nasymnomolptr(nj) - 2
          elseif (iopttype(j).eq.iopt_yf) then
            jj = 3*nasymnomolptr(nj) - 1
          elseif (iopttype(j).eq.iopt_zf) then
            jj = 3*nasymnomolptr(nj)
          elseif (iopttype(j).eq.iopt_radius) then
            jj = 3*numat + nasymnomolptr(nj)
          elseif (iopttype(j).eq.iopt_xcom) then
            lmolTj = .true.
            jj = 1
          elseif (iopttype(j).eq.iopt_ycom) then
            lmolTj = .true.
            jj = 2
          elseif (iopttype(j).eq.iopt_zcom) then
            lmolTj = .true.
            jj = 3
          elseif (iopttype(j).eq.iopt_xqtn) then
            lmolQj = .true.
            jj = 1
          elseif (iopttype(j).eq.iopt_yqtn) then
            lmolQj = .true.
            jj = 2
          elseif (iopttype(j).eq.iopt_zqtn) then
            lmolQj = .true.
            jj = 3
          endif
          if (lmolQi) then
            if (lmolQj) then
              derv2(j,ino) = molQQdrv(jj,ii,nj,ni)
            elseif (lmolTj) then
              derv2(j,ino) = molQTdrv(ii,jj,ni,nj)
            else
              derv2(j,ino) = molQCdrv(jj,ii,ni)
            endif
          elseif (lmolTi) then
            if (lmolQj) then
              derv2(j,ino) = molQTdrv(jj,ii,nj,ni)
            elseif (lmolTj) then
              derv2(j,ino) = molTTdrv(jj,ii,nj,ni)
            else
              derv2(j,ino) = molTCdrv(jj,ii,ni)
            endif
          else
            if (lmolQj) then
              derv2(j,ino) = molQCdrv(ii,jj,nj)
            elseif (lmolTj) then
              derv2(j,ino) = molTCdrv(ii,jj,nj)
            endif
          endif
        enddo
      enddo
    endif
!
!  Wait for data transfer to finish
!
    call MPI_WaitAll(nlocalcomm,Request,StatMPI,MPIerror)
!
!  Place received data into the right locations
!
    do i = 1,nrecv
      ii = nrecvptr(1,i)
      derv2(1:ninternalatm,ii) = datarecv(1:ninternalatm,i)
    enddo
!
!  Deallocate storage for sending
!
    deallocate(StatMPI,stat=status)
    if (status/=0) call deallocate_error('sec0d','StatMPI')
    deallocate(Request,stat=status)
    if (status/=0) call deallocate_error('sec0d','Request')
    deallocate(datasend,stat=status)
    if (status/=0) call deallocate_error('sec0d','datasend')
    deallocate(datarecv,stat=status)
    if (status/=0) call deallocate_error('sec0d','datarecv')
    deallocate(nsendptr,stat=status)
    if (status/=0) call deallocate_error('sec0d','nsendptr')
    deallocate(nsendnode,stat=status)
    if (status/=0) call deallocate_error('sec0d','nsendnode')
    deallocate(nrecvptr,stat=status)
    if (status/=0) call deallocate_error('sec0d','nrecvptr')
    deallocate(nrecvnode,stat=status)
    if (status/=0) call deallocate_error('sec0d','nrecvnode')
    deallocate(ncopyptr,stat=status)
    if (status/=0) call deallocate_error('sec0d','ncopyptr')
  endif
  if (lsdebug) then
    if (ioproc) then
      write(ioout,'(/,''  Symmetrised Second Derivative Matrix  :'',/)')
    endif
    if (lioproconly) then
!
!  Allocate temporary workspace for communication
!
      ntmp = nvar
      ntag = 1
      allocate(tmp2(nvar),stat=status)
      if (status/=0) call outofmemory('sec0d','tmp2')
      allocate(Request(1),stat=status)
      if (status/=0) call outofmemory('sec0d','Request')
      allocate(StatMPI(MPI_Status_Size,1),stat=status)
      if (status/=0) call outofmemory('sec0d','StatMPI')
!
      do i = 1,nvar
        iloc = nvar2local(i)
        if (nvar2node(i).ne.0_i4) then
!
!  Post receive
!
          if (ioproc) then
            nnode = nvar2node(i)
            call MPI_IRecv(tmp2,ntmp,MPI_double_precision,nnode, &
                           ntag,MPI_Comm_World,Request,MPIerror)
          endif
!
!  Pass data to ioproc for writing
!
          if (iloc.gt.0) then
            tmp2(1:nvar) = derv2(1:nvar,iloc)
!
!  Post send
!
            call MPI_ISend(tmp2,ntmp,MPI_double_precision,0, &
                           ntag,MPI_Comm_World,Request,MPIerror)
          endif
          if (ioproc.or.iloc.gt.0) then
            call MPI_WaitAll(1,Request,StatMPI,MPIerror)
          endif
          if (ioproc) then
!
!  Write on I/O node
!
            write(ioout,'(2x,9(f9.4))')(tmp2(j),j=1,nvar)
          endif
        else
          if (ioproc) then
            write(ioout,'(2x,9(f9.4))')(derv2(j,iloc),j=1,nvar)
          endif
        endif
      enddo
!
      deallocate(StatMPI,stat=status)
      if (status/=0) call deallocate_error('sec0d','StatMPI')
      deallocate(Request,stat=status)
      if (status/=0) call deallocate_error('sec0d','Request')
      deallocate(tmp2,stat=status)
      if (status/=0) call deallocate_error('sec0d','tmp2')
    else
      do iproc = 0,nprocs-1
        call mpbarrier
        if (procid.eq.iproc) then
          write(ioout,'(2x,9(i6,3x))')(node2var(j),j=1,nvaronnode)
          do i = 1,nvar
            write(ioout,'(2x,9(f9.4))')(derv2(i,j),j=1,nvaronnode)
          enddo
          write(ioout,'(/)')
        endif
      enddo
    endif
    call mpbarrier
    if (ioproc) then
      write(ioout,'(/)')
    endif
  endif
!
  t2 = g_cpu_time()
  thes = thes + t2 - t1
#ifdef TRACE
  call trace_out('sec0d')
#endif
#else
  call outerror('sec0d called when not compiled with MPI',0_i4)
  call stopnow('sec0d')
#endif
!
  return
  end
