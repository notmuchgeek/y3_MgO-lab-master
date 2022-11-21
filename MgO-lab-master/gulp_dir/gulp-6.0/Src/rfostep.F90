  subroutine rfostep(pvect,hessian,maxhess,lhess2D,xc,fc,gc,eig, &
                     cosm,pnlast,nvar,jnrst,lfrst,loffridge,imode)
!
!  Determine the optimum minimisation/maximisation step according
!  to the RFO approach of Simons et al.
!
!  (1) Diagonalise the hessian, storing the eigenvectors in derv2
!  (2) Find lambda values for two P-RFO problems
!  (3) Form pvect, containing optimum displacement vector.
!
!  cosm = cosm of angle between search direction and gradients
!  eig    = eigenvalue array
!  gc     = gradient vector
!  hessian= hessian matrix in lower-half triangular form
!  jnrst  = cycle counter since hessian was last reset
!  lfrst  = flag to indicate whether to set mode for following
!  nupdate= no. of cycles of updating before hessian recalculation
!  nvar   = number of variables
!  pvect  = search direction
!
!  10/04 Style updated & rsp replaced by lapack call
!   6/05 Deallocation order reversed
!   3/07 linmin renamed to olinmin
!   9/15 Checking of gradients now uses a larger and more consistent
!        value for numerical safety
!   9/15 Trap added for step > stepmax that led to infinite loop problem
!   9/15 Separate stepmax value added for RFO
!   5/16 Trap for small values of rlam vs eig
!   5/16 Check that sign of pvect is consistent with that of the gradient
!   2/17 nmin removed from arguments
!   2/17 maxhess & lhess2D added as arguments
!   3/17 Parallelisation added
!   4/17 Hessian copied to derv2 for parallel diagonalisation
!   4/17 Output of eigenvectors of Hessian now handled in parallel
!   4/17 Reordering of modes during mode following corrected
!   3/18 Parallel I/O corrected
!   9/18 Call to sendall removed
!  10/18 RFO control variables added and stepmaxrfo renamed 
!  12/19 rfotoleig made negative as per the help text
!   4/20 Check on size of derv2 added
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
  use control
  use derivatives
  use general
  use iochannels
  use optimisation
  use parallel
  use times
  implicit none
#ifdef MPI
  include 'mpif.h'
#endif
!
!  Passed variables
!
  integer(i4),      intent(in)                 :: imode
  integer(i4),      intent(out)                :: jnrst
  integer(i4),      intent(in)                 :: nvar
  integer(i4),      intent(in)                 :: maxhess
  logical,          intent(inout)              :: lfrst
  logical,          intent(in)                 :: lhess2D
  logical,          intent(inout)              :: loffridge
  real(dp),         intent(out)                :: cosm
  real(dp),         intent(out)                :: eig(*)
  real(dp)                                     :: fc
  real(dp)                                     :: gc(*)
  real(dp)                                     :: hessian(maxhess,*)
  real(dp)                                     :: pnlast
  real(dp)                                     :: pvect(*)
  real(dp)                                     :: xc(*)
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ierr
  integer(i4)                                  :: igroup
  integer(i4)                                  :: ihdim
  integer(i4)                                  :: ii
  integer(i4)                                  :: il
  integer(i4)                                  :: indi
#ifdef MPI
  integer(i4)                                  :: indig
  integer(i4)                                  :: iproc
#endif
  integer(i4)                                  :: ip
  integer(i4)                                  :: iresid
  integer(i4)                                  :: j
  integer(i4)                                  :: k
  integer(i4)                                  :: nactual
  integer(i4)                                  :: nb
  integer(i4)                                  :: nlower
  integer(i4)                                  :: nmode
  integer(i4)                                  :: nneg
  integer(i4)                                  :: nnegloc
  integer(i4)                                  :: nupper
  integer(i4), dimension(:), allocatable       :: nvptr
  integer(i4), dimension(:), allocatable       :: nvd2ptr
  integer(i4)                                  :: status
  logical                                      :: lokf
  logical                                      :: lsavehess
  real(dp)                                     :: alpha
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: ddot
#ifdef MPI
  real(dp),    dimension(:,:), allocatable     :: dtmp
#endif
  real(dp)                                     :: fcold
  real(dp)                                     :: gl
  real(dp)                                     :: gnrm
  real(dp),    dimension(:,:), allocatable     :: hessave
  real(dp)                                     :: olap
  real(dp)                                     :: omax
  real(dp)                                     :: pdotg
  real(dp)                                     :: pnorm
  real(dp)                                     :: pscal
  real(dp)                                     :: rlam
  real(dp)                                     :: step
  real(dp)                                     :: t1
  real(dp)                                     :: t2
  real(dp),    dimension(:), allocatable       :: overlaps
  real(dp),    dimension(:), allocatable       :: sum1
  real(dp),    dimension(:), allocatable       :: sum2
  real(dp),    dimension(:), allocatable       :: tmp1
  real(dp),    dimension(:), allocatable       :: tmp2
  real(dp),    dimension(:), allocatable       :: wrk
  real(dp)                                     :: trm1
  real(dp)                                     :: trm2
!
  integer                                      :: lwrk
#ifdef MPI
!
!  Local variables in Scalapack/Blacs/MPI integer precision
!
  integer                                      :: ifails
  integer                                      :: idesd(9)
  integer                                      :: ldm
  integer                                      :: nv
  integer                                      :: ntag
  integer                                      :: nnode
  integer                                      :: ntmp
  integer                                      :: Request
  integer                                      :: MPIerror
  integer,     dimension(:), allocatable       :: StatMPI       ! Array for status from MPI
#endif
!
!  Functions
!
  integer(i4)                                  :: ilaenv
!
  t1 = g_cpu_time()
  lsavehess = (nupdate.gt.1.and..not.lhess2D)
  if (lsavehess) then
    ihdim = nvar*(nvar+1)/2
    allocate(hessave(ihdim,1),stat=status)
    if (status/=0) call outofmemory('rfostep','hessave')
  endif
!
!  Allocate local memory
!
  allocate(nvptr(nvar),stat=status)
  if (status/=0) call outofmemory('rfostep','nvptr')
  allocate(nvd2ptr(nvar),stat=status)
  if (status/=0) call outofmemory('rfostep','nvd2ptr')
  allocate(tmp1(3*nvar),stat=status)
  if (status/=0) call outofmemory('rfostep','tmp1')
  allocate(tmp2(nvar),stat=status)
  if (status/=0) call outofmemory('rfostep','tmp2')
  if (nprocs.gt.1) then
    allocate(sum1(2*nvar),stat=status)
    if (status/=0) call outofmemory('rfostep','sum1')
    allocate(sum2(2*nvar),stat=status)
    if (status/=0) call outofmemory('rfostep','sum2')
  endif
!
  cosm = 0.0_dp
  pnorm = 0.0_dp
  gnrm = 0.0_dp
!*************************************
!  Save hessian to avoid corruption  *
!*************************************
  if (lsavehess) then
    do i = 1,ihdim
      hessave(i,1) = hessian(i,1)
    enddo
  endif
#ifdef MPI
  if (nprocs.gt.1) then
!
!  Set local block size
!
    nb = nblocksizevar
    nv = nvar
    ifails = 0
!
!  Set up Blacs descriptors
!
    ldm = maxd2
    call descinit( idesd, nv, nv, nb, nb, 0, 0, iBlacsContext, ldm, ifails )
!
    if (ifails.ne.0) then
      call outerror('initialisation in descinit failed',0_i4)
      call stopnow('rfostep')
    endif
!
!  Check dimensions of derv2
!
    if (nvaronnode.gt.maxd2u.or.nvar.gt.maxd2) then
      maxd2u = max(nvaronnode,maxd2u)
      maxd2  = max(nvar,maxd2)
      call changemaxd2
    endif
!
!  Copy hessian to derv2 
!
    derv2(1:nvar,1:nvaronnode) = hessian(1:nvar,1:nvaronnode)
!
!  Allocate workspace for vectors
!
    allocate(hessave(maxd2,nvaronnode),stat=status)
    if (status/=0) call outofmemory('rfostep','hessave')
    allocate(wrk(1),stat=status)
    if (status/=0) call outofmemory('rfostep','wrk')
!
!  Workspace query
!
    lwrk = -1
    call pdsyev('V','U',nv,derv2,1,1,idesd,eig,hessave,1,1,idesd,wrk,lwrk,ierr)
    lwrk = 2*nint(wrk(1))
    deallocate(wrk,stat=status)
    if (status/=0) call deallocate_error('rfostep','wrk')
    allocate(wrk(lwrk),stat=status)
    if (status/=0) call outofmemory('rfostep','wrk')
!
    call pdsyev('V','U',nv,derv2,1,1,idesd,eig,hessave,1,1,idesd,wrk,lwrk,ierr)
!
    if (ierr.gt.0) then
      call outerror('diagonalisation of Hessian failed in RFO',0_i4)
      call stopnow('rfostep')
    endif
!
!  Copy eigenvectors back to derv2 to match serial code
!
    derv2(1:nvar,1:nvaronnode) = hessave(1:nvar,1:nvaronnode)
!
!  Free workspace for diagonalisation
!
    deallocate(wrk,stat=status)
    if (status/=0) call deallocate_error('rfostep','wrk')
    deallocate(hessave,stat=status)
    if (status/=0) call deallocate_error('rfostep','hessave')
  else
#endif
!************************
!  Diagonalise hessian  *
!************************
    if (lhess2D) then
!
!  Copy Hessian to derv2 so that eigenvectors are returned in the same array as for 1-D case
!
      derv2(1:nvar,1:nvar) = hessian(1:nvar,1:nvar)
!
!  Compute the optimal workspace for diagonalisation
!
      nb = ilaenv(1_i4,'DSYTRD','U',nvar,-1_i4,-1_i4,-1_i4)
      lwrk = max(1_i4,(nb+2_i4)*nvar)
      allocate(wrk(lwrk),stat=status)
      if (status/=0) call outofmemory('rfostep','wrk')
!
      call dsyev('V','U',nvar,derv2,maxd2,eig,wrk,lwrk,ierr)
!
      deallocate(wrk,stat=status)
      if (status/=0) call deallocate_error('rfostep','wrk')
    else
      call dspev('V','U',nvar,hessian(1,1),eig,derv2,maxd2,tmp1,ierr)
    endif
    if (ierr.gt.0) then
      call outerror('diagonalisation of Hessian failed in RFO',0_i4)
      call stopnow('rfostep')
    endif
#ifdef MPI
  endif
#endif
!****************************
!  Check hessian structure  *
!****************************
  nneg = 0
  do i = 1,nvar
    if (eig(i).lt.-rfotoleig) nneg = nneg + 1
  enddo
  if (lopprt.and.ioproc) then
    if (nneg.eq.morder) then
      write(ioout,'(''  ** Hessian has required structure'')')
    else
      write(ioout,'(''  ** Hessian has wrong structure'')')
      write(ioout,'(''  ** Imaginary eigenvectors = '',i3)') nneg
    endif
  endif
  if (nneg.gt.morder.and.loffridge) then
    loffridge = .false.
!
!  Too many imaginary modes - perform off ridge line minimisation
!
    alpha = 0.1_dp
    if (nprocs.gt.1) then
!
!  Find node with eigenvector and place in pvect
!
      if (nvar2node(nneg).eq.procid) then
        nnegloc = nvar2local(nneg)
        do i = 1,nvar
          pvect(i) = derv2(i,nnegloc)
        enddo
      endif
      call sendall(pvect,nvar,nvar2node(nneg),"rfostep","pvect")
    else
      do i = 1,nvar
        pvect(i) = derv2(i,nneg)
      enddo
    endif
    pnorm = sqrt(ddot(nvar,pvect,1_i4,pvect,1_i4))
    pdotg = ddot(nvar,pvect,1_i4,gc,1_i4)
    if (pnorm.gt.1.5_dp*pnlast) then
      pscal = 1.5_dp*pnlast/pnorm
      pscal = sign(pscal,pdotg)
      do i = 1,nvar
        pvect(i) = pvect(i)*pscal
      enddo
    elseif (pdotg.lt.0.0_dp) then
      pscal = sign(1.0_dp,pdotg)
      do i = 1,nvar
        pvect(i) = pvect(i)*pscal
      enddo
    endif
    fcold = fc
    call olinmin(xc,alpha,pvect,nvar,fc,lokf,gc,imode)
    if (lokf) then
      jnrst = nupdate + 1
    else
      fc = fcold
      call daxpy(nvar,-alpha,pvect,1_i4,xc,1_i4)
    endif
  endif
!********************
!  Restore hessian  *
!********************
  if (lsavehess) then
    do i = 1,ihdim
      hessian(i,1) = hessave(i,1)
    enddo
  endif
!****************************************************
!  Transform gradient into the local hessian modes  *
!****************************************************
  tmp1(1:nvar) = 0.0_dp
  tmp2(1:nvar) = 0.0_dp
  do il = 1,nvaronnode
    i = node2var(il)
    do j = 1,nvar
      tmp1(i) = tmp1(i) + derv2(j,il)*gc(j)
    enddo
    tmp2(i) = tmp1(i)*tmp1(i)
  enddo
  if (nprocs.gt.1) then
!
!  Sum tmp1 & tmp2
!
    sum1(1:nvar) = tmp1(1:nvar)
    sum1(nvar+1:2*nvar) = tmp2(1:nvar)
    call sumall(sum1,sum2,2_i4*nvar,"tmp1","rfostep")
    tmp1(1:nvar) = sum2(1:nvar)
    tmp2(1:nvar) = sum2(nvar+1:2*nvar)
  endif
!*******************************
!  Output Hessian if required  *
!*******************************
  if (index(keyword,'hess').ne.0) then
#ifdef MPI
    if (nprocs.eq.1) then
#endif
!
!  Serial
!
      write(ioout,'(/,''  Diagonalised Hessian analysis : '',/)')
      igroup = nvar/3
      iresid = nvar - igroup*3
      indi = 0
      if (igroup.gt.0) then
        do i = 1,igroup
          write(ioout,'(''  Eigenvalues'',3f20.6)') (eig(indi+j),j=1,3)
          write(ioout,'(''  Gradient   '',3f20.6)') (tmp1(indi+j),j=1,3)
          write(ioout,'(/,''  Eigenvectors:'')')
          do j = 1,nvar
            write(ioout,'(3x,i6,4x,3f20.6)') j,(derv2(j,indi+k),k=1,3)
          enddo
          indi = indi + 3
          write(ioout,'(/)')
        enddo
      endif
      if (iresid.gt.0) then
        write(ioout,'(''  Eigenvalues'',3f20.6)') (eig(indi+j),j=1,iresid)
        write(ioout,'(''  Gradient   '',3f20.6)') (tmp1(indi+j),j=1,iresid)
        write(ioout,'(/,''  Eigenvectors:'')')
        do j = 1,nvar
          write(ioout,'(3x,i6,4x,3f20.6)') j,(derv2(j,indi+k),k=1,iresid)
        enddo
      endif
      write(ioout,'(/)')
#ifdef MPI
    else
!
!  Parallel
!
      call mpbarrier
      if (ioproc) then
        write(ioout,'(/,''  Diagonalised Hessian analysis : '',/)')
      endif
      if (lioproconly) then
!
!  Allocate temporary workspace for communication
!
        ntmp = 3_i4*nvar
        ntag = 1
        allocate(dtmp(nvar,3),stat=status)
        if (status/=0) call outofmemory('rfostep','dtmp')
        allocate(StatMPI(MPI_Status_Size),stat=status)
        if (status/=0) call outofmemory('rfostep','StatMPI')
      endif
      igroup = nvar/nblocksizevar
      iresid = nvar - igroup*nblocksizevar
      indi = 0
      indig = 0
      if (igroup.gt.0) then
        do i = 1,igroup
          ii = (i-1)/nprocs
          iproc = i - nprocs*ii - 1
          if (lioproconly.and.(iproc.ne.0_i4)) then
!
!  Post receive
!
            if (ioproc) then
              nnode = iproc
              call MPI_IRecv(dtmp,ntmp,MPI_double_precision,nnode, &
                             ntag,MPI_Comm_World,Request,MPIerror)
            endif
!
!  Pass data to ioproc for writing
!
            if (procid.eq.iproc) then
              do k = 1,3
                dtmp(1:nvar,k) = derv2(1:nvar,indi+k)
              enddo
              indi = indi + 3
!
!  Post send
!
              call MPI_ISend(dtmp,ntmp,MPI_double_precision,0, &
                             ntag,MPI_Comm_World,Request,MPIerror)
            endif
            if (ioproc.or.procid.eq.iproc) then
              call MPI_WaitAll(1,Request,StatMPI,MPIerror)
            endif
            if (ioproc) then
!
!  Write on I/O node
!
              write(ioout,'(''  Eigenvalues'',3f20.6)') (eig(indig+j),j=1,3)
              write(ioout,'(''  Gradient   '',3f20.6)') (tmp1(indig+j),j=1,3)
              write(ioout,'(/,''  Eigenvectors:'')')
              do j = 1,nvar
                write(ioout,'(3x,i6,4x,3f20.6)') j,(dtmp(j,k),k=1,3)
              enddo
              write(ioout,'(/)')
            endif
          else
            call mpbarrier
            if (procid.eq.iproc) then
              write(ioout,'(''  Eigenvalues'',3f20.6)') (eig(indig+j),j=1,3)
              write(ioout,'(''  Gradient   '',3f20.6)') (tmp1(indig+j),j=1,3)
              write(ioout,'(/,''  Eigenvectors:'')')
              do j = 1,nvar
                write(ioout,'(3x,i6,4x,3f20.6)') j,(derv2(j,indi+k),k=1,3)
              enddo
              indi = indi + 3
              write(ioout,'(/)')
            endif
          endif
          indig = indig + 3
        enddo
      endif
      if (iresid.gt.0) then
        ii = igroup/nprocs
        iproc = igroup - nprocs*ii
        if (lioproconly.and.(iproc.ne.0_i4)) then
!
!  Post receive
!
          if (ioproc) then
            ntmp = nvar*iresid
            nnode = iproc
            call MPI_IRecv(dtmp,ntmp,MPI_double_precision,nnode, &
                           ntag,MPI_Comm_World,Request,MPIerror)
          endif
!
!  Pass data to ioproc for writing
!
          if (procid.eq.iproc) then
            do k = 1,iresid
              dtmp(1:nvar,k) = derv2(1:nvar,indi+k)
            enddo
!
!  Post send
!
            ntmp = nvar*iresid
            call MPI_ISend(dtmp,ntmp,MPI_double_precision,0, &
                           ntag,MPI_Comm_World,Request,MPIerror)
          endif
          if (ioproc.or.procid.eq.iproc) then
            call MPI_WaitAll(1,Request,StatMPI,MPIerror)
          endif
          if (ioproc) then
!
!  Write on I/O node
!
            write(ioout,'(''  Eigenvalues'',3f20.6)') (eig(indig+j),j=1,iresid)
            write(ioout,'(''  Gradient   '',3f20.6)') (tmp1(indig+j),j=1,iresid)
            write(ioout,'(/,''  Eigenvectors:'')')
            do j = 1,nvar
              write(ioout,'(3x,i6,4x,3f20.6)') j,(dtmp(j,k),k=1,iresid)
            enddo
            write(ioout,'(/)')
          endif
        else
          call mpbarrier
          if (procid.eq.iproc) then
            write(ioout,'(''  Eigenvalues'',3f20.6)') (eig(indig+j),j=1,iresid)
            write(ioout,'(''  Gradient   '',3f20.6)') (tmp1(indig+j),j=1,iresid)
            write(ioout,'(/,''  Eigenvectors:'')')
            do j = 1,nvar
              write(ioout,'(3x,i6,4x,3f20.6)') j,(derv2(j,indi+k),k=1,iresid)
            enddo
          endif
        endif
      endif
      call mpbarrier
      if (lioproconly) then
        deallocate(StatMPI,stat=status)
        if (status/=0) call deallocate_error('rfostep','StatMPI')
        deallocate(dtmp,stat=status)
        if (status/=0) call deallocate_error('rfostep','dtmp')
      endif
      if (ioproc) then
        write(ioout,'(/)')
      endif
    endif
#endif
  endif
!
!  Analyse local modes to exclude those with zero gradients
!
  nactual = 0
  do i = 1,nvar
    if (abs(tmp1(i)).gt.rfotolgrad) then
      nactual = nactual + 1
      nvptr(nactual) = i
    endif
    nvd2ptr(i) = i
  enddo
!*******************
!  Mode following  *
!*******************
  if (mode.gt.0) then
    allocate(overlaps(nvar),stat=status)
    if (status/=0) call outofmemory('rfostep','overlaps')
!
    if (lfrst) then
!
!  First time vector is given by mode number
!
      ii = nvptr(mode)
      lfrst = .false.
    else
!
!  Compute overlaps with mode
!
      overlaps(1:nactual) = 0.0_dp
      do i = 1,nactual
        ip = nvptr(i)
        if (nvar2node(ip).eq.procid) then
          il = nvar2local(ip)
          olap = 0.0_dp
          do j = 1,nvar
            olap = olap + rmode(j)*derv2(j,il)
          enddo
          overlaps(i) = abs(olap)
        endif
      enddo
!
!  If parallel then globalise overlaps
!
      if (nprocs.gt.1) then
        call sumall(overlaps,sum1,nvar,"overlaps","rfostep")
        overlaps(1:nvar) = sum1(1:nvar)
      endif
!
!  Subsequently find maximum overlap with previous mode
!
      omax = 0.0_dp
      ii = 1
      do i = 1,nactual
        ip = nvptr(i)
        olap = overlaps(i)
        if (olap.gt.omax) then
          omax = olap
          ii = ip
        endif
      enddo
!
      deallocate(overlaps,stat=status)
      if (status/=0) call deallocate_error('rfostep','overlaps')
    endif
!
!  Swap mode to be followed to be mode 1
!
    trm1 = tmp1(ii)
    trm2 = tmp2(ii)
!
!  Instead of moving columns of derv2 just keep track of where they are using a pointer
!
    if (nvar2node(nvd2ptr(ii)).eq.procid) then
      il = nvar2local(nvd2ptr(ii))
      do i = 1,nvar
        rmode(i) = derv2(i,il)
      enddo
    endif
!
!  In parallel communicate rmode
!
    if (nprocs.gt.1) then
      call sendall(rmode,nvar,nvar2node(nvd2ptr(ii)),"rfostep","rmode")
    endif
!
    do i = ii,2,-1
      tmp1(i) = tmp1(i-1)
      tmp2(i) = tmp2(i-1)
      nvd2ptr(i) = i - 1
    enddo
    tmp1(1) = trm1
    tmp2(1) = trm2
    nvd2ptr(1) = ii
!
!  Reset pointer to valid modes
!
    nactual = 0
    do i = 1,nvar
      if (abs(tmp1(i)).gt.rfotolgrad) then
        nactual = nactual + 1
        nvptr(nactual) = i
      endif
    enddo
  endif
!***************
!  Zero pvect  *
!***************
  do i = 1,nvar
    pvect(i) = 0.0_dp
  enddo
!*********************************
!  Generate displacement vector  *
!*********************************
!
!  Minimisation modes
!
  nmode = morder + 1
  nlower = morder + 1
  nupper = nactual
  if (nlower.gt.nupper+1) then
    call outerror('no. of actual local modes is less than order',0_i4)
    call stopnow('rfostep')
  endif
  if (nlower.le.nupper) then
    call lambda(nmode,nlower,nupper,nvptr,nactual,rlam,tmp2,eig)
    do i = nlower,nupper
      ii = nvptr(i)
      gl = tmp1(ii)
!
!  Trap small differences between rlam and eig
!
      if (abs(rlam-eig(ii)).lt.1.0d-8) then
        step = sign(1.0d-3,gl)
      else
        step = gl/(rlam-eig(ii))
      endif
! 
!  Trap steps that exceed maximum size
!
      if (abs(step).gt.rfostepmax) step = sign(rfostepmax,step)
      cosm = cosm - step*gl
      pnorm = pnorm + step*step
      gnrm = gnrm + gl*gl
      if (nvar2node(nvd2ptr(ii)).eq.procid) then
!
!  Mode is local to this node so add contribution to pvect
!
        il = nvar2local(nvd2ptr(ii))
        do j = 1,nvar
          pvect(j) = pvect(j) + step*derv2(j,il)
        enddo
      endif
    enddo
    if (ldebug.and.ioproc) then
      write(ioout,'(''  ** Lambda for minimisation = '',f12.6)') rlam
    endif
  endif
!
!  Maximisation modes
!
  if (morder.gt.0) then
    nmode = morder + 1
    nlower = 1
    nupper = morder
    call lambda(nmode,nlower,nupper,nvptr,nactual,rlam,tmp2,eig)
    do i = nlower,nupper
      ii = nvptr(i)
      gl = tmp1(ii)
!
!  Trap small differences between rlam and eig
!
      if (abs(rlam-eig(ii)).lt.1.0d-8) then
        step = sign(1.0d-3,gl)
      else
        step = gl/(rlam-eig(ii))
      endif
! 
!  Trap steps that exceed maximum size
!
      if (abs(step).gt.rfostepmax) step = sign(rfostepmax,step)
      cosm = cosm + step*gl
      pnorm = pnorm + step*step
      gnrm = gnrm + gl*gl
      if (nvar2node(nvd2ptr(ii)).eq.procid) then
!
!  Mode is local to this node so add contribution to pvect
!
        il = nvar2local(nvd2ptr(ii))
        do j = 1,nvar
          pvect(j) = pvect(j) + step*derv2(j,il)
        enddo
      endif
    enddo
    if (ldebug.and.ioproc) then
      write(ioout,'(''  ** Lambda for maximisation = '',f12.6)') rlam
    endif
  endif
  gnrm = sqrt(gnrm)
  pnorm = sqrt(pnorm)
  cosm = cosm/(pnorm*gnrm)
!
!  Globalise pvect if in parallel
!
  if (nprocs.gt.1) then
    call sumall(pvect,sum1,nvar,"pvect","rfostep")
    pvect(1:nvar) = sum1(1:nvar)
  endif
!
!  Free local memory
!
  if (nprocs.gt.1) then
    deallocate(sum2,stat=status)
    if (status/=0) call deallocate_error('rfostep','sum2')
    deallocate(sum1,stat=status)
    if (status/=0) call deallocate_error('rfostep','sum1')
  endif
  deallocate(tmp2,stat=status)
  if (status/=0) call deallocate_error('rfostep','tmp2')
  deallocate(tmp1,stat=status)
  if (status/=0) call deallocate_error('rfostep','tmp1')
  deallocate(nvd2ptr,stat=status)
  if (status/=0) call deallocate_error('rfostep','nvd2ptr')
  deallocate(nvptr,stat=status)
  if (status/=0) call deallocate_error('rfostep','nvptr')
  if (lsavehess) then
    deallocate(hessave,stat=status)
    if (status/=0) call deallocate_error('rfostep','hessave')
  endif
!
!  Timing
!
  t2 = g_cpu_time()
  tdel = t2 - t1
  thes = thes + tdel
!
  return
  end
