!*****************************************************
!  Serial matrix inversion using standalone routine  *
!*****************************************************
  subroutine matrix_inversion(a,ia,n,wrk,ifail)
!
!  Matrix inverter
!
!  On entry :
!
!  a     = matrix to be inverted
!  ia    = lower dimension of a
!  n     = actual size of matrix to be inverted
!  wrk   = workspace array of length 2*n
!
!  On exit :
!
!  a     = inverse matrix
!  ifail = 0, if OK
!
!   3/14 Renamed from matinv for benefit of ChemShell
!   2/18 Trace added
!
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4)    :: ia
  integer(i4)    :: ifail
  integer(i4)    :: n
  real(dp)       :: a(ia,*)
  real(dp)       :: wrk(*)
!
!  Local variables
!
  integer(i4)    :: i
  integer(i4)    :: i1
  integer(i4)    :: ins
  integer(i4)    :: j
  integer(i4)    :: jj1
  integer(i4)    :: k
  integer(i4)    :: kk
  integer(i4)    :: l
  integer(i4)    :: l1
  integer(i4)    :: m
  integer(i4)    :: nupper
  real(dp)       :: acc
  real(dp)       :: aloc
  real(dp)       :: g_cpu_time
  real(dp)       :: t
  real(dp)       :: t1
  real(dp)       :: t2
#ifdef TRACE
  call trace_in('matrix_inversion')
#endif
!
  t1 = g_cpu_time()
  acc = 1.0d-8
  ifail = 0
  do j = 1,n
    if (j.ne.1) then
      call mxm1(a(j,1),ia,a(1,j),ia,wrk,ia,n-j+1_i4,j-1_i4)
      nupper = n - j + 1
      do l = 1,nupper
        a(j+l-1,j) = a(j+l-1,j) + wrk(l)
      enddo
    endif
    t = abs(a(j,j))
    k = j
    if (j.ne.n) then
      do i = j+1,n
        if (abs(a(i,j)).gt.t) then
          t = abs(a(i,j))
          k = i
        endif
      enddo
    endif
    wrk(j+n) = dble(k)
    if (t.le.acc) then
      ifail = 1
      t2 = g_cpu_time()
      tmati = tmati + t2 - t1
#ifdef TRACE
      call trace_out('matrix_inversion')
#endif
      return
    endif
    if (k.ne.j) then
      do m = 1,n
        t = a(j,m)
        a(j,m) = a(k,m)
        a(k,m) = t
      enddo
    endif
    a(j,j) = 1.0_dp/a(j,j)
    if (j.ne.n) then
      if (j.ne.1) then
        call mxm2(a(1,j+1),ia,a(j,1),ia,wrk,n-j,j-1_i4)
        nupper = n - j
        do l1 = 1,nupper
          a(j,j+l1) = a(j,j+l1) + wrk(l1)
        enddo
      endif
      t = - a(j,j)
      nupper = n - (j+1)
      do i1 = j+1,n
        a(j,i1) = t*a(j,i1)
      enddo
    endif
  enddo
!
!  Use cminv method to solve for a**-1
!
  do k = 2,n
    nupper = k - 1
    do m = 1,nupper
      wrk(m) = 0.0_dp
    enddo
    do j = 1,k-1
      aloc = a(k,j)
      do m = 1,j
        wrk(m) = wrk(m)-aloc*a(j,m)
      enddo
    enddo
    aloc = a(k,k)
    nupper = k - 1
    do m = 1,nupper
      a(k,m) = wrk(m)*aloc
    enddo
  enddo
!
!  Now back substitution
!
  k = n
  do kk = 2,n
    k = k - 1
    jj1 = kk - 1
    call mxm2(a(k+1,1),ia,a(k,k+1),ia,wrk,n,jj1)
    do l = 1,k
      wrk(l) = wrk(l)+a(k,l)
    enddo
    do j = 1,n
      a(k,j) = wrk(j)
    enddo
  enddo
!
!  Multiply solution by inverse of permutation matrix
!
  k = n + 1
  do i = 1,n
    k = k-1
    ins = int(wrk(k+n))
    if (ins.ne.k) then
      do j = 1,n
        wrk(j) = a(j,ins)
        a(j,ins) = a(j,k)
        a(j,k) = wrk(j)
      enddo
    endif
  enddo
  t2 = g_cpu_time()
  tmati = tmati + t2 - t1
#ifdef TRACE
  call trace_out('matrix_inversion')
#endif
  return
  end
!******************************
!  Ancillary matrix routines  *
!******************************
  subroutine mxm1(rarr,ir,sarr,jr,tarr,kr,id1,id2)
!
!  Matrix multiplier
!
  use datatypes
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4)    :: id1
  integer(i4)    :: id2
  integer(i4)    :: ir
  integer(i4)    :: jr
  integer(i4)    :: kr
  real(dp)       :: rarr(*)
  real(dp)       :: sarr(*)
  real(dp)       :: tarr(*)
!
!  Local variables
!
  integer(i4)    :: i
  integer(i4)    :: ia
  integer(i4)    :: j
  real(dp)       :: sum
#ifdef TRACE
  call trace_in('mxm1')
#endif
!
  do i = 1,id1
    ia = i - ir
    sum = 0.0_dp
    do j = 1,id2
      ia = ia + ir
      sum = sum + rarr(ia)*sarr(j)
    enddo
    tarr(i) = sum
  enddo
#ifdef TRACE
  call trace_out('mxm1')
#endif
  return
  end
  subroutine mxm2(rarr,ic,sarr,jc,tarr,id1,id2)
!
!  Matrix multiplier
!
  use datatypes
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4)    :: id1
  integer(i4)    :: id2
  integer(i4)    :: ic
  integer(i4)    :: jc
  real(dp)       :: rarr(*)
  real(dp)       :: sarr(*)
  real(dp)       :: tarr(*)
!
!  Local variables
!
  integer(i4)    :: i
  integer(i4)    :: ia
  integer(i4)    :: ira
  integer(i4)    :: j
  integer(i4)    :: ja
  real(dp)       :: sum
#ifdef TRACE
  call trace_in('mxm2')
#endif
!
  ira = 1 - ic
  do i = 1,id1
    ira = ira + ic
    ia = ira - 1
    ja = 1 - jc
    sum = 0.0_dp
    do j = 1,id2
      ja = ja + jc
      sum = sum + rarr(ia+j)*sarr(ja)
    enddo
    tarr(i) = sum
  enddo
#ifdef TRACE
  call trace_out('mxm2')
#endif
  return
  end
!*************************************************************************
!  Serial or parallel matrix inversion of serial or block cyclic matrix  *
!*************************************************************************
  subroutine matrix_inversion_library(n,nmin,ldm,nblock,matrix,nproc0,ifail)
!
!  Calculates inverse of the matrix input using Scalapack/Blacs
!  in parallel or Lapack in serial.
!
!   9/16 Created
!   2/17 Blocksize argument added
!   2/17 nblock now used instead of nblocksize
!   5/17 Handling of nmin corrected 
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
  integer(i4), intent(in)                      :: n             ! Size of matrix
  integer(i4), intent(in)                      :: nmin          ! Lower starting point for n
  integer(i4), intent(in)                      :: ldm           ! First dimension of matrix
  integer(i4), intent(in)                      :: nblock        ! Block size if in parallel
  real(dp),    intent(inout)                   :: matrix(ldm,*) ! Matrix to be inverted
  integer(i4), intent(in)                      :: nproc0        ! Processor that has the first element of matrix
  integer(i4), intent(out)                     :: ifail         ! Flag indicating status of call
!
!  Local variables in Scalapack/Blacs integer precision
!
#ifdef MPI
  integer                                      :: ifails
  integer                                      :: idesc(9)
  integer                                      :: ldms
  integer,     dimension(:), allocatable       :: iwrk
  integer                                      :: liwrk
  integer                                      :: lwrk
  integer                                      :: nb
  integer                                      :: nm
  integer                                      :: np
  integer                                      :: np0
  integer                                      :: ns
#endif
  integer,     dimension(:), allocatable       :: ipivot
!
!  Local variables in GULP precision
!
  integer(i4)                                  :: i
  integer(i4)                                  :: j
  integer(i4)                                  :: k
  integer(i4)                                  :: status
  real(dp),    dimension(:), allocatable       :: dpacked
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: t1p
  real(dp)                                     :: t2p
  real(dp),    dimension(:), allocatable       :: wrk
#ifdef TRACE
  call trace_in('matrix_inversion_library')
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
    if (status/=0) call outofmemory('matrix_inversion','ipivot')
!
!  Set local block size
!
    nb = nblock
    np = nprocs
    np0 = nproc0
    ns = n + nmin - 1
    nm = nmin
    ifails = 0
    ldms = ldm
!
!  Set up Blacs descriptor
!
    call descinit( idesc, ns, ns, nb, nb, 0, np0, iBlacsContext, ldms, ifails )
!
    if (ifails.ne.0) then
      call outerror('initialisation in descinit failed',0_i4)
      call stopnow('matrix_inversion')
    endif
    ifail = ifails
!
!  Check for failure
!
    if (ifail.ne.0) then
      call outerror('Blacs has failed to initialise a descriptor',0_i4)
      call stopnow('matrix_inversion')
    endif
!     
!  Factorise matrix using Scalapack
!     
    call pdgetrf(n,n,matrix,nmin,nmin,idesc,ipivot,ifails)  
    ifail = ifails
    if (ifail.eq.0) then
!
!  Initial dummy allocation of workspace
!
      lwrk = 1
      liwrk = 1
      allocate(iwrk(liwrk),stat=status)
      if (status/=0) call outofmemory('matrix_inversion','iwrk')
      allocate(wrk(lwrk),stat=status)
      if (status/=0) call outofmemory('matrix_inversion','wrk')
!
!  Query to find workspace needed
!
      lwrk = -1
      liwrk = -1
      ifails = 0
      call pdgetri(n,matrix,nmin,nmin,idesc,ipivot,wrk,lwrk,iwrk,liwrk,ifails)  
!
!  Use double the amount suggested by the query to avoid out of bounds issues
!
      lwrk = 2*nint(wrk(1))
      liwrk = 2*iwrk(1)
!
!  Reallocate workspace to size needed
!
      deallocate(wrk,stat=status)
      if (status/=0) call deallocate_error('matrix_inversion','wrk')
      deallocate(iwrk,stat=status)
      if (status/=0) call deallocate_error('matrix_inversion','iwrk')
      allocate(iwrk(liwrk),stat=status)
      if (status/=0) call outofmemory('matrix_inversion','iwrk')
      allocate(wrk(lwrk),stat=status)
      if (status/=0) call outofmemory('matrix_inversion','wrk')
!
!  Form inverse
!
      call pdgetri(n,matrix,nmin,nmin,idesc,ipivot,wrk,lwrk,iwrk,liwrk,ifails)  
      ifail = ifails
!
!  Deallocate workspace
!
      deallocate(wrk,stat=status)
      if (status/=0) call deallocate_error('matrix_inversion','wrk')
      deallocate(iwrk,stat=status)
      if (status/=0) call deallocate_error('matrix_inversion','iwrk')
    endif
    deallocate(ipivot,stat=status)
    if (status/=0) call deallocate_error('matrix_inversion','ipivot')
#else
    call outerror('Parallel matrix inversion called without MPI',0_i4)
    call stopnow('matrix_inversion')
#endif
  else
!**********************************
!  Serial inversion using lapack  *
!**********************************
!
!  Allocate workspace for inversion
!
    allocate(dpacked(n*(n+1)/2),stat=status)
    if (status/=0) call outofmemory('matrix_inversion','dpacked')
    allocate(ipivot(n),stat=status)
    if (status/=0) call outofmemory('matrix_inversion','ipivot')
    allocate(wrk(3*n),stat=status)
    if (status/=0) call outofmemory('matrix_inversion','wrk')
!
!  Transfer data to packed storage
!     
    k = 0                                    
    do i = nmin,nmin+n-1
      do j = nmin,i
        k = k + 1                            
        dpacked(k) = matrix(j,i)
      enddo    
    enddo  
!     
!  Factorise matrix
!     
    call dsptrf('U',n,dpacked,ipivot,ifail)  
    if (ifail.eq.0) then
!
!  Form inverse
!
      call dsptri('U',n,dpacked,ipivot,wrk,ifail)
!
!  Transfer data back
!
      k = 0
      do i = nmin,nmin+n-1
        do j = nmin,i
          k = k + 1
          matrix(j,i) = dpacked(k)            
          matrix(i,j) = dpacked(k)
        enddo   
      enddo 
    endif
!
!  Free workspace
!
    deallocate(wrk,stat=status)
    if (status/=0) call deallocate_error('matrix_inversion','wrk')
    deallocate(ipivot,stat=status)
    if (status/=0) call deallocate_error('matrix_inversion','ipivot')
    deallocate(dpacked,stat=status)
    if (status/=0) call deallocate_error('matrix_inversion','dpacked')
  endif
!
  t2p = g_cpu_time()
!
#ifdef TRACE
  call trace_out('matrix_inversion_library')
#endif
  return
  end
!**************************************************************
!  Serial or parallel matrix inversion of shell-shell matrix  *
!**************************************************************
  subroutine matrix_inversion_shells(n,nmin,ldm,matrix,nsh,nshloc,ifail)
!
!  Calculates inverse of the shell-shell matrix using Scalapack/Blacs.
!  in parallel or Lapack in serial.
!
!  For the serial case this is just a matter of calling the 
!  standard matrix inversion routine. 
!  For the parallel case the matrix has to be reorganised to 
!  conform to the block cyclic distribution required by 
!  Scalapack. This is because the data is block cyclic in the
!  total number of particles, but not necessarily in the shells.
!
!   9/16 Created 
!   2/17 Blocksize added to call to matrix_inversion_library
!   5/17 Use of nmin now made to handle case where this changes
!        the data distribution between nodes. 
!   5/17 Arguments added to handle bulk vs defect call
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
  use trace,        only : trace_in, trace_out
#endif
  implicit none
#ifdef MPI
  include 'mpif.h'
#endif
!
!  Passed variables
!
  integer(i4), intent(in)                      :: n             ! Size of matrix
  integer(i4), intent(in)                      :: nmin          ! Lower starting point for n
  integer(i4), intent(in)                      :: ldm           ! First dimension of matrix
  integer(i4), intent(in)                      :: nsh           ! Number of shells (global)
  integer(i4), intent(in)                      :: nshloc        ! Number of shells (local)
  real(dp),    intent(inout)                   :: matrix(ldm,*) ! Matrix to be inverted
  integer(i4), intent(out)                     :: ifail         ! Flag indicating status of call
!
!  Local variables 
!
#ifdef MPI
  integer(i4)                                  :: i
  integer(i4)                                  :: icount
  integer(i4), dimension(:), allocatable       :: itmp
  integer                                      :: MPIerror
  integer(i4)                                  :: j
  integer(i4)                                  :: nlocalcomm
  integer(i4)                                  :: nblocks
  integer(i4)                                  :: ncores
  integer(i4)                                  :: nploop
  integer(i4)                                  :: nmaxt
  integer(i4), dimension(:),   allocatable     :: nshnode       ! Number of shells per node for assessing pattern
  integer(i4), dimension(:),   allocatable     :: nshbcyc       ! Number of shells per node according to block cyclic pattern
  integer(i4), dimension(:),   allocatable     :: nsh2loc
  integer(i4), dimension(:,:), allocatable     :: nsh2node
  integer(i4), dimension(:,:), allocatable     :: nsh2nodeptr
  integer(i4), dimension(:),   allocatable     :: nnoderecv     ! Pointer to node that is receiving
  integer(i4), dimension(:),   allocatable     :: nnoderecvptr  ! Pointer to columns on node that is receiving
  integer(i4), dimension(:),   allocatable     :: nnodesend     ! Pointer to node that is sending
  integer(i4), dimension(:),   allocatable     :: nnodesendptr  ! Pointer to columns on node that is sending
  integer(i4), dimension(:),   allocatable     :: Request       ! Array for requests to MPI
  integer(i4), dimension(:,:), allocatable     :: StatMPI       ! Array for status from MPI
  integer(i4)                                  :: node
  integer(i4)                                  :: nproc0
  integer(i4)                                  :: ntransfer
  integer(i4)                                  :: status
  logical                                      :: lblockcyclic
  real(dp),    dimension(:,:), allocatable     :: matrix2       ! Copy of matrix for block cyclic redistribution
#endif
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: t1p
  real(dp)                                     :: t2p
#ifdef TRACE
  call trace_in('matrix_inversion_shells')
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
      if (status/=0) call outofmemory('matinv_shells','nsh2node')
      allocate(nsh2nodeptr(2,nsh),stat=status)
      if (status/=0) call outofmemory('matinv_shells','nsh2nodeptr')
      allocate(nsh2loc(0:nprocs-1),stat=status)
      if (status/=0) call outofmemory('matinv_shells','nsh2loc')
      allocate(nshnode(0:nprocs-1),stat=status)
      if (status/=0) call outofmemory('matinv_shells','nshnode')
      allocate(nshbcyc(0:nprocs-1),stat=status)
      if (status/=0) call outofmemory('matinv_shells','nshnode')
      allocate(itmp(0:nprocs-1),stat=status)
      if (status/=0) call outofmemory('matinv_shells','itmp')
!*********************************************************************
!  Check whether current data conforms to block cyclic distribution  *
!*********************************************************************
      nshnode(0:nprocs-1) = 0
      nshnode(procid) = nshloc
!
!  Globalise nshnode
!
      call isumall(nshnode(0),itmp(0),nprocs,"matinv_shells","nshnode")
      nshnode(0:nprocs-1) = itmp(0:nprocs-1)
!
!  Compute block cyclic distribution of shells
!
      nshbcyc(0:nprocs-1) = 0
!
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
!
      allocate(nnoderecv(nmaxt),stat=status)
      if (status/=0) call outofmemory('matinv_shells','nnoderecv')
      allocate(nnoderecvptr(nmaxt),stat=status)
      if (status/=0) call outofmemory('matinv_shells','nnoderecvptr')
      allocate(nnodesend(nmaxt),stat=status)
      if (status/=0) call outofmemory('matinv_shells','nnodesend')
      allocate(nnodesendptr(nmaxt),stat=status)
      if (status/=0) call outofmemory('matinv_shells','nnodesendptr')
      allocate(Request(2*nmaxt),stat=status)
      if (status/=0) call outofmemory('matinv_shells','Request')
      allocate(StatMPI(MPI_Status_Size,2*nmaxt),stat=status)
      if (status/=0) call outofmemory('matinv_shells','StatMPI')
!
!  Create matrix copy to hold new distribution of 3 x number of shells on node after redistribution
!
      allocate(matrix2(ldm,3*nshbcyc(procid)),stat=status)
      if (status/=0) call outofmemory('matinv_shells','matrix2')
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
          call MPI_IRecv(matrix2(1,3*(nnoderecvptr(i)-1)+1),3*ldm,MPI_double_precision,nnodesend(i), &
                         i,MPI_Comm_World,Request(nlocalcomm),MPIerror)
        elseif (procid.eq.nnodesend(i)) then
          nlocalcomm = nlocalcomm + 1
          call MPI_ISend(matrix(1,3*(nnodesendptr(i)-1)+1),3*ldm,MPI_double_precision,nnoderecv(i), &
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
      call matrix_inversion_library(n,1_i4,ldm,3_i4*nblocksize,matrix2,0_i4,ifail)
!
!  Set up redistribution back to original distribution
!
      nlocalcomm = 0
      do i = 1,ntransfer
        if (procid.eq.nnoderecv(i)) then
          nlocalcomm = nlocalcomm + 1
          call MPI_ISend(matrix2(1,3*(nnoderecvptr(i)-1)+1),3*ldm,MPI_double_precision,nnodesend(i), &
                         i,MPI_Comm_World,Request(nlocalcomm),MPIerror)
        elseif (procid.eq.nnodesend(i)) then
          nlocalcomm = nlocalcomm + 1
          call MPI_IRecv(matrix(1,3*(nnodesendptr(i)-1)+1),3*ldm,MPI_double_precision,nnoderecv(i), &
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
      if (status/=0) call deallocate_error('matinv_shells','matrix2')
      deallocate(StatMPI,stat=status)
      if (status/=0) call deallocate_error('matinv_shells','StatMPI')
      deallocate(Request,stat=status)
      if (status/=0) call deallocate_error('matinv_shells','Request')
      deallocate(nnodesendptr,stat=status)
      if (status/=0) call deallocate_error('matinv_shells','nnodesendptr')
      deallocate(nnodesend,stat=status)
      if (status/=0) call deallocate_error('matinv_shells','nnodesend')
      deallocate(nnoderecvptr,stat=status)
      if (status/=0) call deallocate_error('matinv_shells','nnoderecvptr')
      deallocate(nnoderecv,stat=status)
      if (status/=0) call deallocate_error('matinv_shells','nnoderecv')
      deallocate(itmp,stat=status)
      if (status/=0) call deallocate_error('matinv_shells','itmp')
      deallocate(nshbcyc,stat=status)
      if (status/=0) call deallocate_error('matinv_shells','nshbcyc')
      deallocate(nshnode,stat=status)
      if (status/=0) call deallocate_error('matinv_shells','nshnode')
      deallocate(nsh2loc,stat=status)
      if (status/=0) call deallocate_error('matinv_shells','nsh2loc')
      deallocate(nsh2nodeptr,stat=status)
      if (status/=0) call deallocate_error('matinv_shells','nsh2nodeptr')
      deallocate(nsh2node,stat=status)
      if (status/=0) call deallocate_error('matinv_shells','nsh2node')
    else
!***************************************
!  Parallel inversion using scalapack  *
!***************************************
      call matrix_inversion_library(n,nmin,ldm,3_i4*nblocksize,matrix,nproc0,ifail)
    endif
#else
    call outerror('Parallel matrix inversion called without MPI',0_i4)
    call stopnow('matinv_shells')
#endif
  else
!****************************************
!  Serial inversion using main routine  *
!****************************************
    call matrix_inversion_library(n,nmin,ldm,3_i4*nblocksize,matrix,0_i4,ifail)
  endif
!
  t2p = g_cpu_time()
#ifdef TRACE
  call trace_out('matrix_inversion_shells')
#endif
!
  return
  end
