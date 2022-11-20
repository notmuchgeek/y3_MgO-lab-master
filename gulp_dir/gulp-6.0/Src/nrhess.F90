  subroutine nrhess(hesinv,maxhess,lhess2D,diag,nvar,ldiag)
!
!  Calculate hessian matrix by one of two methods:
!
!    (1) Exact inverse of second derivative matrix
!    (2) Diagonal elements only
!
!  The first corresponds to Newton-Raphson, while the second is only
!  employed if the Hessian is singular.
!
!  12/07 Unused variables removed
!   2/17 nmin removed from arguments
!   2/17 maxhess & lhess2D added as arguments
!   3/17 Parallelisation added
!   4/17 Parallel handling of output added
!   4/17 Globalise of hessian diagonal elements added
!   6/17 Default on-diagonal element when matrix is ill-condition and a value
!        is zero now set to one.
!   1/18 Trace added
!   3/18 Parallel I/O corrected
!  12/20 Changes to accommodate gfortran v10
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
!  Julian Gale, CIC, Curtin University, December 2020
!
  use control
  use derivatives
  use general
  use iochannels
  use parallel
  use times
#ifdef TRACE
  use trace,       only : trace_in, trace_out
#endif
  implicit none
#ifdef MPI
  include 'mpif.h'
#endif
!
!  Passed variables
!
  integer(i4)                                  :: nvar
  integer(i4)                                  :: maxhess
  logical                                      :: ldiag
  logical                                      :: lhess2D
  real(dp)                                     :: diag(*)
  real(dp)                                     :: hesinv(maxhess,*)
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ihdim
  integer(i4)                                  :: il
  integer(i4)                                  :: ind
  integer(i4)                                  :: info
  integer(i4)                                  :: j
  integer(i4), dimension(:), allocatable       :: kpvt
  integer(i4)                                  :: status
  logical                                      :: lhdebug
  real(dp)                                     :: g_cpu_time
#ifdef MPI
  real(dp)                                     :: htmp(1)
#endif
  real(dp)                                     :: t1
  real(dp)                                     :: t2
  real(dp)                                     :: t1i
  real(dp)                                     :: t2i
  real(dp),    dimension(:), allocatable       :: temp
  real(dp),    dimension(:), allocatable       :: wrk
!
  integer                                      :: nb
  integer                                      :: lwrk
#ifdef MPI
!
!  Local variables in Scalapack/Blacs/MPI integer precision
!
  integer                                      :: ifails
  integer                                      :: idesh(9)
  integer,     dimension(:), allocatable       :: ipivot
  integer,     dimension(:), allocatable       :: iwrk
  integer                                      :: liwrk
  integer                                      :: ldm
  integer                                      :: np
  integer                                      :: nv
  integer                                      :: MPIerror
  integer                                      :: ntag
  integer                                      :: nnode
  integer                                      :: ntmp
  integer                                      :: Request
  integer,     dimension(:), allocatable       :: StatMPI       ! Array for status from MPI
#endif
!
!  Functions
!
  integer(i4)                                  :: ilaenv
#ifdef TRACE
  call trace_in('nrhess')
#endif
!
  t1 = g_cpu_time()
  lhdebug = (index(keyword,'hess').ne.0)
  ldiag = .false.
!
!  Pack hessian and save diagonal elements
!
  if (lhess2D) then
    diag(1:nvar) = 0.0_dp
    do il = 1,nvaronnode
      i = node2var(il)
      do j = 1,nvar
        hesinv(j,il) = derv2(j,il)
      enddo
      diag(i) = derv2(i,il)
    enddo
    if (nprocs.gt.1) then
      allocate(temp(nvar),stat=status)
      if (status/=0) call outofmemory('nrhess','temp')
!
      call sumall(diag,temp,nvar,"diag","nrhess")
      diag(1:nvar) = temp(1:nvar)
!
      deallocate(temp,stat=status)
      if (status/=0) call deallocate_error('nrhess','temp')
    endif
  else
    ind = 0
    do i = 1,nvar
      do j = 1,i
        ind = ind + 1
        hesinv(ind,1) = derv2(j,i)
      enddo
      diag(i) = derv2(i,i)
    enddo
  endif
  if (lhdebug) then
    if (ioproc) then
      write(ioout,'(/,'' Hessian Matrix :'',/)')
    endif
#ifdef MPI
    if (lioproconly.and.lhess2D) then
!
!  Allocate temporary workspace for communication
!
      ntmp = nvar
      ntag = 1
      allocate(temp(nvar),stat=status)
      if (status/=0) call outofmemory('nrhess','temp')
      allocate(StatMPI(MPI_Status_Size),stat=status)
      if (status/=0) call outofmemory('nrhess','StatMPI')
    endif
    call mpbarrier
#endif
    if (lhess2D) then
#ifdef MPI
      if (nprocs.gt.1) then
        do i = 1,nvar
          il = nvar2local(i)
          if (lioproconly.and.nvar2node(i).ne.0_i4) then
!
!  Post receive
!
            if (ioproc) then
              nnode = nvar2node(i)
              call MPI_IRecv(temp,ntmp,MPI_double_precision,nnode, &
                             ntag,MPI_Comm_World,Request,MPIerror)
            endif
!
!  Pass data to ioproc for writing
!
            if (il.gt.0) then
              temp(1:nvar) = hesinv(1:nvar,il)
!
!  Post send
!
              call MPI_ISend(temp,ntmp,MPI_double_precision,0, &
                             ntag,MPI_Comm_World,Request,MPIerror)
            endif
            if (ioproc.or.il.gt.0) then
              call MPI_WaitAll(1,Request,StatMPI,MPIerror)
            endif
            if (ioproc) then
!
!  Write on I/O node
!
              write(ioout,'(2x,10(f12.4))')(temp(j),j=1,nvar)
            endif
          else
            if (il.gt.0) then
              write(ioout,'(2x,10(f12.4))')(hesinv(j,il),j=1,nvar)
            endif
          endif
          call mpbarrier
        enddo
      else
#endif
        do i = 1,nvar
          write(ioout,'(2x,10(f12.4))')(hesinv(j,i),j=1,nvar)
        enddo
#ifdef MPI
      endif
#endif
    else
      do i = 1,nvar
        ind = i*(i-1)/2
        write(ioout,'(2x,10(f12.4))')(hesinv(ind+j,1),j=1,i)
      enddo
    endif
#ifdef MPI
    if (lioproconly.and.lhess2D) then
      deallocate(StatMPI,stat=status)
      if (status/=0) call deallocate_error('nrhess','StatMPI')
      deallocate(temp,stat=status)
      if (status/=0) call deallocate_error('nrhess','temp')
    endif
#endif
  endif
  t1i = g_cpu_time()
#ifdef MPI
  if (nprocs.gt.1) then
!***************************************
!  Parallel inversion using scalapack  *
!***************************************
    allocate(ipivot(4*nvar),stat=status)
    if (status/=0) call outofmemory('nrhess','ipivot')
!
!  Set local block size
!
    nb = nblocksizevar
    np = nprocs
    nv = nvar
    ifails = 0
!
!  Set up Blacs descriptors
!
    ldm = maxhess
    call descinit( idesh, nv, nv, nb, nb, 0, 0, iBlacsContext, ldm, ifails )
!
    if (ifails.ne.0) then
      call outerror('initialisation in descinit failed',0_i4)
      call stopnow('nrhess')
    endif
!
!  Factorise matrix using Scalapack
!
    call pdgetrf(nv,nv,hesinv,1,1,idesh,ipivot,ifails)
    info = ifails
!
!  Check for singularities
!
    if (info.gt.0) then
      call outwarning('ill conditioned Hessian - using diagonal elements only',0_i4)
      ldiag = .true.
      hesinv(1:nvar,1:nvaronnode) = 0.0_dp
      do il = 1,nvaronnode
        i = node2var(il)
        if (abs(diag(i)).lt.1.0d-6) diag(i) = 1.0d-6
        hesinv(i,il) = 1.0_dp/diag(i)
      enddo
    else
!***********************
!  Complete inversion  *
!***********************
!
!  Initial dummy allocation of workspace
!
      lwrk = 1
      liwrk = 1
      allocate(iwrk(liwrk),stat=status)
      if (status/=0) call outofmemory('nrhess','iwrk')
      allocate(wrk(lwrk),stat=status)
      if (status/=0) call outofmemory('nrhess','wrk')
!
!  Query to find workspace needed
!
      lwrk = -1
      liwrk = -1
      ifails = 0
      call pdgetri(nv,hesinv,1,1,idesh,ipivot,wrk,lwrk,iwrk,liwrk,ifails)
!
!  Use double the amount suggested by the query to avoid out of bounds issues
!
      lwrk = 2*nint(wrk(1))
      liwrk = 2*iwrk(1)
!
!  Reallocate workspace to size needed
!
      deallocate(wrk,stat=status)
      if (status/=0) call deallocate_error('nrhess','wrk')
      deallocate(iwrk,stat=status)
      if (status/=0) call deallocate_error('nrhess','iwrk')
      allocate(iwrk(liwrk),stat=status)
      if (status/=0) call outofmemory('nrhess','iwrk')
      allocate(wrk(lwrk),stat=status)
      if (status/=0) call outofmemory('nrhess','wrk')
!
!  Form inverse
!
      call pdgetri(nv,hesinv,1,1,idesh,ipivot,wrk,lwrk,iwrk,liwrk,ifails)
!
!  Deallocate workspace
!
      deallocate(wrk,stat=status)
      if (status/=0) call deallocate_error('nrhess','wrk')
      deallocate(iwrk,stat=status)
      if (status/=0) call deallocate_error('nrhess','iwrk')
    endif
    deallocate(ipivot,stat=status)
    if (status/=0) call deallocate_error('nrhess','ipivot')
  else
#endif
!
!  Factorise matrix
!
    allocate(kpvt(nvar),stat=status)
    if (status/=0) call outofmemory('nrhess','kpvt')
    if (lhess2D) then
      call dgetrf(nvar,nvar,hesinv,maxhess,kpvt,info)
    else
      call dsptrf('U',nvar,hesinv,kpvt,info)
    endif
!
!  Check for singularities
!
    if (info.gt.0) then
      call outwarning('ill conditioned Hessian - using diagonal elements only',0_i4)
      ldiag = .true.
      if (lhess2D) then
        hesinv(1:nvar,1:nvar) = 0.0_dp
        do i = 1,nvar
          if (abs(diag(i)).lt.1.0d-6) diag(i) = 1.0d-6
          hesinv(i,i) = 1.0_dp/diag(i)
        enddo
      else
        ihdim = nvar*(nvar+1)/2
        do i = 1,ihdim
          hesinv(i,1) = 0.0_dp
        enddo
        do i = 1,nvar
          ind = i*(i+1)/2
          if (abs(diag(i)).lt.1.0d-6) diag(i) = 1.0d-6
          hesinv(ind,1) = 1.0_dp/diag(i)
        enddo
      endif
    else
!
!  Complete inversion
!
      if (lhess2D) then
!
!  Find block size
!
        nb = ilaenv( 1, 'DGETRI', ' ', nvar, -1, -1, -1 )
!
        lwrk = nb*nvar
        allocate(wrk(lwrk),stat=status)
        if (status/=0) call outofmemory('nrhess','wrk')
!
        call dgetri(nvar,hesinv,maxhess,kpvt,wrk,lwrk,info)
!
        deallocate(wrk,stat=status)
        if (status/=0) call deallocate_error('nrhess','wrk')
      else
        allocate(temp(3*nvar),stat=status)
        if (status/=0) call outofmemory('nrhess','temp')
        call dsptri('U',nvar,hesinv,kpvt,temp,info)
        deallocate(temp,stat=status)
        if (status/=0) call deallocate_error('nrhess','temp')
      endif
    endif
    deallocate(kpvt,stat=status)
    if (status/=0) call deallocate_error('nrhess','kpvt')
#ifdef MPI
  endif
#endif
  t2i = g_cpu_time()
  tmati = tmati + t2i - t1i
  if (lhdebug) then
    if (index(keyword,'hessd').ne.0) then
      if (ioproc) then
        write(ioout,'(/,'' Inverse Hessian Diagonal Elements :'',/)')
      endif
#ifdef MPI
      if (lioproconly.and.lhess2D) then
!
!  Allocate temporary workspace for communication
!
        ntmp = 1
        ntag = 1
        allocate(StatMPI(MPI_Status_Size),stat=status)
        if (status/=0) call outofmemory('nrhess','StatMPI')
      endif
      call mpbarrier
#endif
      if (lhess2D) then
#ifdef MPI
        if (nprocs.gt.1) then
          do i = 1,nvar
            il = nvar2local(i)
            if (lioproconly.and.nvar2node(i).ne.0_i4) then
!
!  Post receive
!
              if (ioproc) then
                nnode = nvar2node(i)
                call MPI_IRecv(htmp,ntmp,MPI_double_precision,nnode, &
                               ntag,MPI_Comm_World,Request,MPIerror)
              endif
!
!  Pass data to ioproc for writing
!
              if (il.gt.0) then
                htmp(1) = hesinv(i,il)
!
!  Post send
!
                call MPI_ISend(htmp,ntmp,MPI_double_precision,0, &
                               ntag,MPI_Comm_World,Request,MPIerror)
              endif
              if (ioproc.or.il.gt.0) then
                call MPI_WaitAll(1,Request,StatMPI,MPIerror)
              endif
              if (ioproc) then
!
!  Write on I/O node
!
                write(ioout,'(2x,f12.8)') htmp(1)
              endif
            else
              if (il.gt.0) then
                write(ioout,'(2x,f12.8)') hesinv(i,il)
              endif
            endif
            call mpbarrier
          enddo
        else
#endif
          do i = 1,nvar
            write(ioout,'(2x,f12.8)') hesinv(i,i)
          enddo
#ifdef MPI
        endif
#endif
      else
        do i = 1,nvar
          ind = i*(i+1)/2
          write(ioout,'(2x,f12.8)') hesinv(ind,1)
        enddo
      endif
    else
      if (ioproc) then
        write(ioout,'(/,'' Inverse Hessian Matrix :'',/)')
      endif
#ifdef MPI
      if (lioproconly.and.lhess2D) then
!
!  Allocate temporary workspace for communication
!
        ntag = 1
        allocate(temp(nvar),stat=status)
        if (status/=0) call outofmemory('nrhess','temp')
        allocate(StatMPI(MPI_Status_Size),stat=status)
        if (status/=0) call outofmemory('nrhess','StatMPI')
      endif
      call mpbarrier
#endif
      if (lhess2D) then
#ifdef MPI
        if (nprocs.gt.1) then
          do i = 1,nvar
            il = nvar2local(i)
            if (lioproconly.and.nvar2node(i).ne.0_i4) then
!
!  Post receive
!
              if (ioproc) then
                ntmp = i
                nnode = nvar2node(i)
                call MPI_IRecv(temp,ntmp,MPI_double_precision,nnode, &
                               ntag,MPI_Comm_World,Request,MPIerror)
              endif
!
!  Pass data to ioproc for writing
!
              if (il.gt.0) then
                ntmp = i
                temp(1:i) = hesinv(1:i,il)
!
!  Post send
!
                call MPI_ISend(temp,ntmp,MPI_double_precision,0, &
                               ntag,MPI_Comm_World,Request,MPIerror)
              endif
              if (ioproc.or.il.gt.0) then
                call MPI_WaitAll(1,Request,StatMPI,MPIerror)
              endif
              if (ioproc) then
!
!  Write on I/O node
!
                write(ioout,'(2x,10(f12.8))')(temp(j),j=1,i)
              endif
            else
              if (il.gt.0) then
                write(ioout,'(2x,10(f12.8))')(hesinv(j,il),j=1,i)
              endif
            endif
            call mpbarrier
          enddo
        else
#endif
          do i = 1,nvar
            write(ioout,'(2x,10(f12.8))')(hesinv(j,i),j=1,i)
          enddo
#ifdef MPI
        endif
#endif
      else
        do i = 1,nvar
          ind = i*(i-1)/2
          write(ioout,'(2x,10(f12.8))')(hesinv(ind+j,1),j=1,i)
        enddo
      endif
    endif
#ifdef MPI
    if (lioproconly.and.lhess2D) then
      deallocate(StatMPI,stat=status)
      if (status/=0) call deallocate_error('nrhess','StatMPI')
      deallocate(temp,stat=status)
      if (status/=0) call deallocate_error('nrhess','temp')
    endif
#endif
    if (ioproc) then
      write(ioout,'(/)')
    endif
  endif
!
  t2 = g_cpu_time()
  tdel = t2 - t1
  thes = thes + tdel - t2i + t1i
#ifdef TRACE
  call trace_out('nrhess')
#endif
!
  return
  end
