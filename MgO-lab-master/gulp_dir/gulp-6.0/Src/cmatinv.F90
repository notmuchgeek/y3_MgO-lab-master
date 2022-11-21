!*********************************************************************************
!  Serial or parallel complex matrix inversion of serial or block cyclic matrix  *
!*********************************************************************************
  subroutine cmatrix_inversion_library(n,nmin,ldm,matrix,nproc0,ifail)
!
!  Calculates inverse of the complex matrix input using Scalapack/Blacs
!  in parallel or Lapack in serial.
!
!  11/16 Created from matrix_inversion_library
!   5/17 Modified to allow for nmin being used
!   7/17 nproc0 added
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
  use iochannels
  use parallel
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4),  intent(in)                      :: n             ! Size of matrix
  integer(i4),  intent(in)                      :: nmin          ! Lower starting point for n
  integer(i4),  intent(in)                      :: ldm           ! First dimension of matrix
  complex(dpc), intent(inout)                   :: matrix(ldm,*) ! Matrix to be inverted
  integer(i4),  intent(in)                      :: nproc0        ! Processor that has the first element of matrix
  integer(i4),  intent(out)                     :: ifail         ! Flag indicating status of call
!
!  Local variables in Scalapack/Blacs integer precision
!
#ifdef MPI
  integer                                       :: ifails
  integer                                       :: idesc(9)
  integer                                       :: ldms
  integer,      dimension(:), allocatable       :: iwrk
  integer                                       :: liwrk
  integer                                       :: lcwrk
  integer                                       :: nb
  integer                                       :: np
  integer                                       :: np0
  integer                                       :: ns
#endif
  integer,      dimension(:), allocatable       :: ipivot
!
!  Local variables in GULP precision
!
  integer(i4)                                   :: i
  integer(i4)                                   :: j
  integer(i4)                                   :: k
  integer(i4)                                   :: status
  complex(dpc), dimension(:), allocatable       :: cpacked
  real(dp)                                      :: g_cpu_time
  real(dp)                                      :: t1p
  real(dp)                                      :: t2p
  complex(dpc), dimension(:), allocatable       :: cwrk
#ifdef TRACE
  call trace_in('cmatrix_inversion_library')
#endif
!
  t1p = g_cpu_time()
!
!  Initialise status flag
!
  ifail = 0
!
  if (nprocs.gt.1) then
#ifdef MPI
!***************************************
!  Parallel inversion using scalapack  *
!***************************************
    allocate(ipivot(4*n),stat=status)
    if (status/=0) call outofmemory('cmatrix_inversion','ipivot')
!
!  Set up Blacs descriptor
!
    nb = nblocksize
    np = nprocs
    np0 = nproc0
    ns = n + nmin - 1
    ifails = 0
    ldms = ldm
    call descinit( idesc, ns, ns, 3*nb, 3*nb, 0, np0, iBlacsContext, ldms, ifails )
!
    if (ifails.ne.0) then
      call outerror('initialisation in descinit failed',0_i4)
      call stopnow('cmatrix_inversion')
    endif
    ifail = ifails
!
!  Check for failure
!
    if (ifail.ne.0) then
      call outerror('Blacs has failed to initialise a descriptor',0_i4)
      call stopnow('cmatrix_inversion')
    endif
!     
!  Factorise matrix using Scalapack
!     
    call pzgetrf(n,n,matrix,nmin,nmin,idesc,ipivot,ifails)  
    ifail = ifails
    if (ifail.eq.0) then
!
!  Initial dummy allocation of workspace
!
      lcwrk = 1
      liwrk = 1
      allocate(iwrk(liwrk),stat=status)
      if (status/=0) call outofmemory('cmatrix_inversion','iwrk')
      allocate(cwrk(lcwrk),stat=status)
      if (status/=0) call outofmemory('cmatrix_inversion','cwrk')
!
!  Query to find workspace needed
!
      lcwrk = -1
      liwrk = -1
      ifails = 0
      call pzgetri(n,matrix,nmin,nmin,idesc,ipivot,cwrk,lcwrk,iwrk,liwrk,ifails)  
!
!  Use double the amount suggested by the query to avoid out of bounds issues
!
      lcwrk = 2*nint(real(cwrk(1)))
      liwrk = 2*iwrk(1)
!
!  Reallocate workspace to size needed
!
      deallocate(cwrk,stat=status)
      if (status/=0) call deallocate_error('cmatrix_inversion','cwrk')
      deallocate(iwrk,stat=status)
      if (status/=0) call deallocate_error('cmatrix_inversion','iwrk')
      allocate(iwrk(liwrk),stat=status)
      if (status/=0) call outofmemory('cmatrix_inversion','iwrk')
      allocate(cwrk(lcwrk),stat=status)
      if (status/=0) call outofmemory('cmatrix_inversion','cwrk')
!
!  Form inverse
!
      call pzgetri(n,matrix,nmin,nmin,idesc,ipivot,cwrk,lcwrk,iwrk,liwrk,ifails)  
      ifail = ifails
!
!  Deallocate workspace
!
      deallocate(cwrk,stat=status)
      if (status/=0) call deallocate_error('cmatrix_inversion','cwrk')
      deallocate(iwrk,stat=status)
      if (status/=0) call deallocate_error('cmatrix_inversion','iwrk')
    endif
    deallocate(ipivot,stat=status)
    if (status/=0) call deallocate_error('cmatrix_inversion','ipivot')
#else
    call outerror('Parallel matrix inversion called without MPI',0_i4)
    call stopnow('cmatrix_inversion')
#endif
  else
!**********************************
!  Serial inversion using lapack  *
!**********************************
!
!  Allocate workspace for inversion
!
    allocate(cpacked(n*(n+1)/2),stat=status)
    if (status/=0) call outofmemory('cmatrix_inversion','cpacked')
    allocate(ipivot(n),stat=status)
    if (status/=0) call outofmemory('cmatrix_inversion','ipivot')
    allocate(cwrk(2*n),stat=status)
    if (status/=0) call outofmemory('cmatrix_inversion','cwrk')
!
!  Transfer data to packed storage
!     
    k = 0                                    
    do i = nmin,nmin+n-1
      do j = nmin,i
        k = k + 1                            
        cpacked(k) = matrix(j,i)
      enddo    
    enddo  
!     
!  Factorise matrix
!     
    call zsptrf('U',n-nmin+1,cpacked,ipivot,ifail)  
    if (ifail.eq.0) then
!
!  Form inverse
!
      call zsptri('U',n-nmin+1,cpacked,ipivot,cwrk,ifail)
!
!  Transfer data back
!
      k = 0
      do i = nmin,nmin+n-1
        do j = nmin,i
          k = k + 1
          matrix(j,i) = cpacked(k)            
          matrix(i,j) = conjg(cpacked(k))
        enddo   
      enddo 
    endif
!           
!  Free workspace
!     
    deallocate(cwrk,stat=status)
    if (status/=0) call deallocate_error('cmatrix_inversion','cwrk')
    deallocate(ipivot,stat=status)
    if (status/=0) call deallocate_error('cmatrix_inversion','ipivot')
    deallocate(cpacked,stat=status)
    if (status/=0) call deallocate_error('cmatrix_inversion','cpacked')
  endif
!
  t2p = g_cpu_time()
#ifdef TRACE
  call trace_out('cmatrix_inversion_library')
#endif
!
  return
  end
!**********************************************************************
!  Serial or parallel complex matrix inversion of shell-shell matrix  *
!**********************************************************************
  subroutine cmatrix_inversion_shells(n,nmin,ldm,matrix,nsh,nshloc,ifail)
!
!  Calculates inverse of the shell-shell complex matrix using Scalapack/Blacs.
!  in parallel or Lapack in serial.
!
!  For the serial case this is just a matter of calling the 
!  standard matrix inversion routine. 
!  For the parallel case the matrix has to be reorganised to 
!  conform to the block cyclic distribution required by 
!  Scalapack. This is because the data is block cyclic in the
!  total number of particles, but not necessarily in the shells.
!
!  11/16 Created from matrix_inversion_shells
!   5/17 Use of nmin now made to handle case where this changes
!        the data distribution between nodes.
!   7/17 Algorithm changed to using first shell as offset in
!        processor grid.
!   7/17 nmin > 1 now assumes that data is as per block cyclic
!        distribution and hasn't been moved
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
  use iochannels
  use parallel
#ifdef MPI
  use shells,       only : ncore, nshptr
#endif
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
#ifdef MPI
  include 'mpif.h'
#endif
!
!  Passed variables
!
  integer(i4),  intent(in)                      :: n             ! Size of matrix
  integer(i4),  intent(in)                      :: nmin          ! Lower starting point for n
  integer(i4),  intent(in)                      :: ldm           ! First dimension of matrix
  integer(i4),  intent(in)                      :: nsh           ! Number of shells (global)
  integer(i4),  intent(in)                      :: nshloc        ! Number of shells (local)
  complex(dpc), intent(inout)                   :: matrix(ldm,*) ! Matrix to be inverted
  integer(i4),  intent(out)                     :: ifail         ! Flag indicating status of call
!
!  Local variables 
!
#ifdef MPI
  integer(i4)                                   :: i
  integer(i4)                                   :: icount
  integer                                       :: MPIerror
  integer(i4)                                   :: ir
  integer(i4)                                   :: is
  integer(i4),  dimension(:),   allocatable     :: itmp
  integer(i4)                                   :: j
  integer(i4)                                   :: nblocks
  integer(i4)                                   :: nlocalcomm
  integer(i4)                                   :: ncores
  integer(i4)                                   :: nploop
  integer(i4)                                   :: nmaxt
  integer(i4),  dimension(:),   allocatable     :: nshnode       ! Number of shells per node for assessing pattern
  integer(i4),  dimension(:),   allocatable     :: nshbcyc       ! Number of shells per node according to block cyclic pattern
  integer(i4),  dimension(:),   allocatable     :: nsh2loc
  integer(i4),  dimension(:,:), allocatable     :: nsh2node
  integer(i4),  dimension(:,:), allocatable     :: nsh2nodeptr
  integer(i4),  dimension(:),   allocatable     :: nnoderecv     ! Pointer to node that is receiving
  integer(i4),  dimension(:),   allocatable     :: nnoderecvptr  ! Pointer to columns on node that is receiving
  integer(i4),  dimension(:),   allocatable     :: nnodesend     ! Pointer to node that is sending
  integer(i4),  dimension(:),   allocatable     :: nnodesendptr  ! Pointer to columns on node that is sending
  integer(i4),  dimension(:),   allocatable     :: Request       ! Array for requests to MPI
  integer(i4),  dimension(:,:), allocatable     :: StatMPI       ! Array for status from MPI
  integer(i4)                                   :: ntransfer
  integer(i4)                                   :: node
  integer(i4)                                   :: nproc0
  integer(i4)                                   :: status
  logical                                       :: lblockcyclic
  complex(dpc), dimension(:,:), allocatable     :: matrix2       ! Copy of matrix for block cyclic redistribution
#endif
  real(dp)                                      :: g_cpu_time
  real(dp)                                      :: t1p
  real(dp)                                      :: t2p
#ifdef TRACE
  call trace_in('cmatrix_inversion_shells')
#endif
!
  t1p = g_cpu_time()
!
!  Initialise status flag
!
  ifail = 0
!
  if (nprocs.gt.1) then
#ifdef MPI
!
!  Check whether we can use existing data distribution
!
    if (nmin.eq.1) then
      ncores = ncore
      nblocks = ncores/(nblocksize*nprocs)
      ncores = ncores - nblocks*nblocksize*nprocs
      nproc0 = ncores/nblocksize
      ncores = ncores - nproc0*nblocksize + 1
    else
!
!  If nmin isn't one then we presume that data is positioned in array as per block cyclic arrangement
!
      nproc0 = 0
      ncores = 1
    endif
!
!  If ncores is 1 then the distribution is properly aligned with a
!  block cyclic distribution beginning on process nproc0
!
    lblockcyclic = (ncores.eq.1)
!
    if (.not.lblockcyclic) then
!
!  Data is not block cyclic and so needs redistributing
!
!  Allocate workspace
!
      allocate(nsh2node(2,nsh),stat=status)
      if (status/=0) call outofmemory('cmatinv_shells','nsh2node')
      allocate(nsh2nodeptr(2,nsh),stat=status)
      if (status/=0) call outofmemory('cmatinv_shells','nsh2nodeptr')
      allocate(nsh2loc(0:nprocs-1),stat=status)
      if (status/=0) call outofmemory('cmatinv_shells','nsh2loc')
      allocate(nshnode(0:nprocs-1),stat=status)
      if (status/=0) call outofmemory('cmatinv_shells','nshnode')
      allocate(nshbcyc(0:nprocs-1),stat=status)
      if (status/=0) call outofmemory('cmatinv_shells','nshnode')
      allocate(itmp(0:nprocs-1),stat=status)
      if (status/=0) call outofmemory('cmatinv_shells','itmp')
!*********************************************************************
!  Check whether current data conforms to block cyclic distribution  *
!*********************************************************************
      nshnode(0:nprocs-1) = 0
      nshnode(procid) = nshloc
!
!  Globalise nshnode
!
      call isumall(nshnode(0),itmp(0),nprocs,"cmatinv_shells","nshnode")
      nshnode(0:nprocs-1) = itmp(0:nprocs-1)
!
!  Compute block cyclic distribution of shells
!
      nshbcyc(0:nprocs-1) = 0
      if (nmin.eq.1) then
!
!  For nmin start from beginning
!
        node = 0
        icount = 0
      else
!
!  For nmin > 1 update node and icount to reflect starting point in distribution
!
        ncore = nmin - 1
        nblocks = ncore/(3_i4*nblocksize)
        nploop = nblocks/nprocs
        node = nblocks - nploop*nprocs
        icount = ncore - nblocks*(3_i4*nblocksize)
      endif
!
      nsh2loc(0:nprocs-1) = 0
      do i = 1,nsh
        nshbcyc(node) = nshbcyc(node) + 1
        nsh2node(1,i) = atom2node(nshptr(i))
        nsh2node(2,i) = node
        nsh2nodeptr(1,i) = atom2local(nshptr(i))
!
!  Find processor that has this shell
!
        nsh2loc(nsh2node(1,i)) = nsh2loc(nsh2node(1,i)) + 1
        nsh2nodeptr(1,i) = nsh2loc(nsh2node(1,i))
        nsh2nodeptr(2,i) = nshbcyc(node)
        icount = icount + 1
        if (icount.eq.nblocksize) then
          icount = 0
          node = node + 1
          if (node.eq.nprocs) node = 0
        endif
      enddo
!
!  Find maximum number of transfers
!
      nmaxt = 0
      do i = 1,nsh
        if (nsh2node(1,i).ne.nsh2node(2,i)) then
          nmaxt = nmaxt + 1
        endif
      enddo
!***********************************************************************
!  If data is not block cyclic then create a new matrix and rearrange  *
!***********************************************************************
      allocate(nnoderecv(nmaxt),stat=status)
      if (status/=0) call outofmemory('cmatinv_shells','nnoderecv')
      allocate(nnoderecvptr(nmaxt),stat=status)
      if (status/=0) call outofmemory('cmatinv_shells','nnoderecvptr')
      allocate(nnodesend(nmaxt),stat=status)
      if (status/=0) call outofmemory('cmatinv_shells','nnodesend')
      allocate(nnodesendptr(nmaxt),stat=status)
      if (status/=0) call outofmemory('cmatinv_shells','nnodesendptr')
      allocate(Request(2*nmaxt),stat=status)
      if (status/=0) call outofmemory('cmatinv_shells','Request')
      allocate(StatMPI(MPI_Status_Size,2*nmaxt),stat=status)
      if (status/=0) call outofmemory('cmatinv_shells','StatMPI')
!
!  Create matrix copy to hold new distribution of 3 x number of shells on node after redistribution
!
      allocate(matrix2(ldm,3*nshbcyc(procid)),stat=status)
      if (status/=0) call outofmemory('cmatinv_shells','matrix2')
!
!  Work out how to redistribute columns
!
      ntransfer = 0
      do i = 1,nsh
        if (nsh2node(1,i).ne.nsh2node(2,i)) then
          ntransfer = ntransfer + 1
          nnoderecv(ntransfer) = nsh2node(2,i)
          nnodesend(ntransfer) = nsh2node(1,i)
          nnoderecvptr(ntransfer) = nsh2nodeptr(2,i)
          nnodesendptr(ntransfer) = nsh2nodeptr(1,i)
        endif
      enddo
!
!  Set up data redistribution using MPI
!
      nlocalcomm = 0
      do i = 1,ntransfer
        if (procid.eq.nnoderecv(i)) then
          nlocalcomm = nlocalcomm + 1
          call MPI_IRecv(matrix2(1,3*(nnoderecvptr(i)-1)+1),3*ldm,MPI_double_complex,nnodesend(i), &
                         i,MPI_Comm_World,Request(nlocalcomm),MPIerror)
        elseif (procid.eq.nnodesend(i)) then
          nlocalcomm = nlocalcomm + 1
          call MPI_ISend(matrix(1,3*(nnodesendptr(i)-1)+1),3*ldm,MPI_double_complex,nnoderecv(i), &
                         i,MPI_Comm_World,Request(nlocalcomm),MPIerror)
        endif
      enddo
!
!  Copy the local elements of matrix while waiting for the data to transfer
!
      do i = 1,nsh
        if (nsh2node(1,i).eq.procid) then
          if (nsh2node(1,i).eq.nsh2node(2,i)) then
            do j = 1,3
              matrix2(1:ldm,3*(nsh2nodeptr(2,i)-1)+j) = matrix(1:ldm,3*(nsh2nodeptr(1,i)-1)+j)
            enddo
          endif
        endif
      enddo
!
!  Wait for data transfer to complete
!
      call MPI_WaitAll(nlocalcomm,Request,StatMPI,MPIerror) 
!
!  Invert in block cyclic form
!
      call cmatrix_inversion_library(n,1_i4,ldm,matrix2,0_i4,ifail)
!
!  Set up redistribution back to original distribution
!
      nlocalcomm = 0
      do i = 1,ntransfer
        if (procid.eq.nnoderecv(i)) then
          nlocalcomm = nlocalcomm + 1
          call MPI_ISend(matrix2(1,3*(nnoderecvptr(i)-1)+1),3*ldm,MPI_double_complex,nnodesend(i), &
                         i,MPI_Comm_World,Request(nlocalcomm),MPIerror)
        elseif (procid.eq.nnodesend(i)) then
          nlocalcomm = nlocalcomm + 1
          call MPI_IRecv(matrix(1,3*(nnodesendptr(i)-1)+1),3*ldm,MPI_double_complex,nnoderecv(i), &
                         i,MPI_Comm_World,Request(nlocalcomm),MPIerror)
        endif
      enddo
!
!  Copy the local elements of matrix2 back to matrix while waiting for the data to transfer
!
      do i = 1,nsh
        if (nsh2node(1,i).eq.procid) then
          if (nsh2node(1,i).eq.nsh2node(2,i)) then
            do j = 1,3
              matrix(1:ldm,3*(nsh2nodeptr(1,i)-1)+j) = matrix2(1:ldm,3*(nsh2nodeptr(2,i)-1)+j)
            enddo
          endif
        endif
      enddo
!
!  Wait for data transfer to complete
!
      call MPI_WaitAll(nlocalcomm,Request,StatMPI,MPIerror) 
!
!  Deallocate workspace
!
      deallocate(matrix2,stat=status)
      if (status/=0) call deallocate_error('cmatinv_shells','matrix2')
      deallocate(StatMPI,stat=status)
      if (status/=0) call deallocate_error('cmatinv_shells','StatMPI')
      deallocate(Request,stat=status)
      if (status/=0) call deallocate_error('cmatinv_shells','Request')
      deallocate(nnodesendptr,stat=status)
      if (status/=0) call deallocate_error('cmatinv_shells','nnodesendptr')
      deallocate(nnodesend,stat=status)
      if (status/=0) call deallocate_error('cmatinv_shells','nnodesend')
      deallocate(nnoderecvptr,stat=status)
      if (status/=0) call deallocate_error('cmatinv_shells','nnoderecvptr')
      deallocate(nnoderecv,stat=status)
      if (status/=0) call deallocate_error('cmatinv_shells','nnoderecv')
      deallocate(itmp,stat=status)
      if (status/=0) call deallocate_error('cmatinv_shells','itmp')
      deallocate(nshbcyc,stat=status)
      if (status/=0) call deallocate_error('cmatinv_shells','nshbcyc')
      deallocate(nshnode,stat=status)
      if (status/=0) call deallocate_error('cmatinv_shells','nshnode')
      deallocate(nsh2loc,stat=status)
      if (status/=0) call deallocate_error('cmatinv_shells','nsh2loc')
      deallocate(nsh2nodeptr,stat=status)
      if (status/=0) call deallocate_error('cmatinv_shells','nsh2nodeptr')
      deallocate(nsh2node,stat=status)
      if (status/=0) call deallocate_error('cmatinv_shells','nsh2node')
    else
!***************************************
!  Parallel inversion using scalapack  *
!***************************************
      call cmatrix_inversion_library(n,nmin,ldm,matrix,nproc0,ifail)
    endif
#else
    call outerror('Parallel matrix inversion called without MPI',0_i4)
    call stopnow('cmatinv_shells')
#endif
  else
!****************************************
!  Serial inversion using main routine  *
!****************************************
    call cmatrix_inversion_library(n,nmin,ldm,matrix,0_i4,ifail)
  endif
!
  t2p = g_cpu_time()
#ifdef TRACE
  call trace_out('cmatrix_inversion_shells')
#endif
!
  return
  end
