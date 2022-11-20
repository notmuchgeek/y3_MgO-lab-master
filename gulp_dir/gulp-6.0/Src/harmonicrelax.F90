  subroutine harmonicrelax(xc,fc,gc,hessian,maxhess,lhess2D,imode)
!
!  Estimate the energy change on relaxation using the harmonic approximation
!
!  imode   = 1 => bulk calculation
!  imode   = 2 => defect calculation
!
!   3/14 Created from minimise
!   1/17 Calls to sec0/sec3 modified to handle distributed memory
!   2/17 nmin removed from arguments
!   2/17 nvar removed from arguments
!   3/17 maxhess & lhess2D added as arguments
!   4/17 Parallelisation of second derivatives added
!   2/18 Trace added
!   7/20 Trapping of geometry problem flag on return from funct added
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
!  Julian Gale, CIC, Curtin University, July 2020
!
#ifdef MPI
  use configurations, only : maxvar
#endif
  use control
  use current
  use defects
  use derivatives,    only : derv2
  use energies,       only : erelax, efreeze, fcsave
  use iochannels
  use optimisation
  use parallel
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif

  implicit none
!
!  Passed variables
!
  integer(i4),          intent(in)            :: imode
  integer(i4),          intent(in)            :: maxhess
  logical,              intent(in)            :: lhess2D
  real(dp),             intent(inout)         :: fc
  real(dp),             intent(inout)         :: gc(nvar)
  real(dp),             intent(inout)         :: hessian(maxhess,*)
  real(dp),             intent(inout)         :: xc(nvar)
!
!  Local variables
!
  integer(i4)                                 :: i
  integer(i4)                                 :: ind
  integer(i4)                                 :: info
  integer(i4)                                 :: iflag
  integer(i4)                                 :: j
  integer(i4)                                 :: status
  integer(i4),              allocatable, save :: kpvt(:)
  real(dp)                                    :: ddot
  real(dp)                                    :: fcin
  real(dp)                                    :: funct1
  real(dp),                 allocatable, save :: pvect(:)
  real(dp),                 allocatable, save :: wrk(:)
!
  integer                                     :: nb
  integer                                     :: lwrk
#ifdef MPI
!
!  Local variables in Scalapack/Blacs integer precision
!
  integer                                     :: ifails
  integer                                     :: idesg(9)
  integer                                     :: idesh(9)
  integer,     dimension(:), allocatable      :: ipivot
  integer,     dimension(:), allocatable      :: iwrk
  integer                                     :: liwrk
  integer                                     :: ldm
  integer                                     :: np
  integer                                     :: nv
#endif
!
!  Functions
!
  integer(i4)                                 :: ilaenv
!
!  If number of variables is zero then relaxation energy is zero
!
  if (nvar.eq.0) then
    erelax = 0.0_dp
    return
  endif
#ifdef TRACE
  call trace_in('harmonicrelax')
#endif
!
!  Allocate local memory
!
  allocate(pvect(nvar),stat=status)
  if (status/=0) call outofmemory('harmonicrelax','pvect')
  if (.not.lhess2D.or.nprocs.eq.1) then
    allocate(kpvt(nvar),stat=status)
    if (status/=0) call outofmemory('harmonicrelax','kpvt')
  endif
!
!  Calculate initial function, gradients and necessary second derivatives
!
  iflag = 2
  fcin = fc
  if (imode.eq.1) then
    call funct(iflag,nvar,xc,fc,gc)
!
!  Check return flag for problem with geometry
!
    if (iflag.eq.-2) then
      call outerror('geometry has become unphysical during harmonic relaxation',0_i4)
      call stopnow('harmonicrelax')
    endif
  else
    call deffun(iflag,nvar,xc,fc,gc)
  endif
  fcsave = fc
  if (lfreeze) then
    efreeze = fcin - fc
    fc = fc + efreeze
  else
    efreeze = 0.0_dp
  endif
!***************************
!  Generate exact hessian  *
!***************************
  iflag = 2
  if (imode.eq.1) then
    call funct(iflag,nvar,xc,funct1,gc)
!
!  Check return flag for problem with geometry
!
    if (iflag.eq.-2) then
      call outerror('geometry has become unphysical during harmonic relaxation',0_i4)
      call stopnow('harmonicrelax')
    endif
  else
    call deffun(iflag,nvar,xc,funct1,gc)
  endif
!
!  Set Hessian
!
  if (imode.eq.1) then
    if (ndim.ge.1) then
      if (lfreeze) then
        call sec3f
      else
        if (nprocs.gt.1) then
          call sec3d
        else
          call sec3
        endif
      endif
    else
      if (lfreeze) then
        call sec0f
      else
        if (nprocs.gt.1) then
          call sec0d
        else
          call sec0
        endif
      endif
    endif
  else
    call defsec
  endif
!
!  Pack or copy hessian 
!
  if (lhess2D) then
    hessian(1:nvar,1:nvaronnode) = derv2(1:nvar,1:nvaronnode)
  else
    ind = 0
    do i = 1,nvar
      do j = 1,i
        ind = ind + 1
        hessian(ind,1) = derv2(j,i)
      enddo
    enddo
  endif
!
!  Factorise matrix
!
  if (lhess2D) then
#ifdef MPI
    if (nprocs.gt.1) then
!***************************************
!  Parallel inversion using scalapack  *
!***************************************
      allocate(ipivot(4*nvar),stat=status)
      if (status/=0) call outofmemory('harmonicrelax','ipivot')
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
        call stopnow('harmonicrelax')
      endif
!
!  Factorise matrix using Scalapack
!
      call pdgetrf(nv,nv,hessian,1,1,idesh,ipivot,ifails)
      info = ifails
    else
#endif
      call dgetrf(nvar,nvar,hessian,maxhess,kpvt,info)
#ifdef MPI
    endif
#endif
  else
    call dsptrf('U',nvar,hessian,kpvt,info)
  endif
!
!  Check for singularities
!
  if (info.gt.0) then
    call outerror('hessian inversion failed during estimation of relaxation energy',0_i4)
    call stopnow('harmonicrelax')
  else
!***********************
!  Complete inversion  *
!***********************
    if (lhess2D) then
#ifdef MPI
      if (nprocs.gt.1) then
!
!  Initial dummy allocation of workspace
!
        lwrk = 1
        liwrk = 1
        allocate(iwrk(liwrk),stat=status)
        if (status/=0) call outofmemory('harmonicrelax','iwrk')
        allocate(wrk(lwrk),stat=status)
        if (status/=0) call outofmemory('harmonicrelax','wrk')
!
!  Query to find workspace needed
!
        lwrk = -1
        liwrk = -1
        ifails = 0
        call pdgetri(nv,hessian,1,1,idesh,ipivot,wrk,lwrk,iwrk,liwrk,ifails)
!
!  Use double the amount suggested by the query to avoid out of bounds issues
!
        lwrk = 2*nint(wrk(1))
        liwrk = 2*iwrk(1)
!
!  Reallocate workspace to size needed
!
        deallocate(wrk,stat=status)
        if (status/=0) call deallocate_error('harmonicrelax','wrk')
        deallocate(iwrk,stat=status)
        if (status/=0) call deallocate_error('harmonicrelax','iwrk')
        allocate(iwrk(liwrk),stat=status)
        if (status/=0) call outofmemory('harmonicrelax','iwrk')
        allocate(wrk(lwrk),stat=status)
        if (status/=0) call outofmemory('harmonicrelax','wrk')
!
!  Form inverse
!
        call pdgetri(nv,hessian,1,1,idesh,ipivot,wrk,lwrk,iwrk,liwrk,ifails)
!
!  Deallocate workspace
!
        deallocate(wrk,stat=status)
        if (status/=0) call deallocate_error('harmonicrelax','wrk')
        deallocate(iwrk,stat=status)
        if (status/=0) call deallocate_error('harmonicrelax','iwrk')
      else
#endif
!
!  Find block size
!
        nb = ilaenv( 1, 'DGETRI', ' ', nvar, -1, -1, -1 )
        lwrk = nb*nvar
!
        allocate(wrk(lwrk),stat=status)
        if (status/=0) call outofmemory('harmonicrelax','wrk')
!
        call dgetri(nvar,hessian,maxhess,kpvt,wrk,lwrk,info)
!
        deallocate(wrk,stat=status)
        if (status/=0) call deallocate_error('harmonicrelax','wrk')
#ifdef MPI
      endif
#endif
    else
      allocate(wrk(3*nvar),stat=status)
      if (status/=0) call outofmemory('nrhess','wrk')
!
      call dsptri('U',nvar,hessian,kpvt,wrk,info)
!
      deallocate(wrk,stat=status)
      if (status/=0) call deallocate_error('harmonicrelax','wrk')
    endif
  endif
#ifdef MPI
  if (allocated(ipivot)) then
    deallocate(ipivot,stat=status)
    if (status/=0) call deallocate_error('harmonicrelax','ipivot')
  endif
#endif
!
!  Multiply inverse Hessian by gradient to get displacement vector
!
  if (lhess2D) then
#ifdef MPI
    if (nprocs.gt.1) then
!**************************************************
!  Call pblas routine for matrix-vector multiply  *
!**************************************************
!
!  Set local block size
!
      nb = nblocksizevar
      nv = nvar
      ifails = 0
!
      ldm = maxvar
      call descinit( idesg, nv, 1, nb, nb, 0, 0, iBlacsContext, ldm, ifails )
!
      if (ifails.ne.0) then
        call outerror('initialisation in descinit failed',0_i4)
        call stopnow('harmonicrelax')
      endif
!
!  Call pblas
!
      call pdgemv('N',nvar,nvar,1.0d0,hessian,1,1,idesh,gc,1,1,idesg,1,0.0d0,pvect,1,1,idesg,1)
!
!  Broadcase the result from node 0
!
      call sendall(pvect,nvar,0_i4,"harmonicrelax","pvect")
    else
#endif
      call dgemv('N',nvar,nvar,1.0_dp,hessian,maxhess,gc,1_i4,0.0_dp,pvect,1_i4)
#ifdef MPI
    endif
#endif
  else
    call dspmv('U',nvar,1.0_dp,hessian,gc,1_i4,0.0_dp,pvect,1_i4)
  endif
!
!  Compute predicted energy charge within the harmonic approximation
!
  erelax = - ddot(nvar,pvect,1_i4,gc,1_i4)
  erelax = 0.5_dp*erelax
!
!  Free local memory
!
  if (.not.lhess2D.or.nprocs.eq.1) then
    deallocate(kpvt,stat=status)
    if (status/=0) call deallocate_error('harmonicrelax','kpvt')
  endif
  deallocate(pvect,stat=status)
  if (status/=0) call deallocate_error('harmonicrelax','pvect')
#ifdef TRACE
  call trace_out('harmonicrelax')
#endif
!
  return
  end
