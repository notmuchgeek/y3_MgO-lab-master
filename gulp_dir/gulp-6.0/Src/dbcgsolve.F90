  subroutine dbcgsolve(mode,n,nloc,nlocptr,A,maxA,b,x,tol,itmax,iter,err)
!
!  Dense version of bcgsolve
!
!  If mode = 0 => general matrix
!  If mode > 0 => matrix is of special implicit form for EEM
!
  use datatypes
  use iochannels
  use control,       only : keyword
  use parallel,      only : ioproc
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif

  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)    :: mode               ! Flag that specifies the matrix form to be used
  integer(i4), intent(in)    :: n                  ! Total dimension of problem
  integer(i4), intent(in)    :: nloc               ! Number of columns of A local to this node
  integer(i4), intent(in)    :: nlocptr(nloc)      ! Pointer from local columns to global number
  integer(i4), intent(in)    :: itmax              ! Maximum number of iterations
  integer(i4), intent(out)   :: iter               ! Number of iterations taken
  integer(i4), intent(in)    :: maxA               ! LHS dimension of listA and A
  real(dp),    intent(in)    :: A(maxA,nloc)       ! Values of non-zero elements of column of A
  real(dp),    intent(in)    :: tol                ! Tolerance
  real(dp),    intent(in)    :: b(n)               ! Left-hand vector
  real(dp),    intent(inout) :: x(n)               ! Solution vector : guess on input, actual on output
  real(dp),    intent(out)   :: err                ! Error flag
!
!  Local variables
!
  integer(i4)                :: j
  logical                    :: converged
  real(dp)                   :: ak
  real(dp)                   :: akden
  real(dp)                   :: bk
  real(dp)                   :: bkden
  real(dp)                   :: bknum
  real(dp)                   :: bnrm
  real(dp)                   :: ddot
  real(dp),    allocatable   :: p(:)
  real(dp),    allocatable   :: pp(:)
  real(dp),    allocatable   :: r(:)
  real(dp),    allocatable   :: rr(:)
  real(dp),    allocatable   :: z(:)
  real(dp),    allocatable   :: zz(:)
#ifdef TRACE
  call trace_in('dbcgsolve')
#endif
!
!  Allocate workspace
!
  allocate(p(n))
  allocate(pp(n))
  allocate(r(n))
  allocate(rr(n))
  allocate(z(n))
  allocate(zz(n))

  call denseAxV(mode,n,nloc,nlocptr,A,maxA,x,r,.false.)

  do j = 1,n
    r(j)  = b(j) - r(j)
    rr(j) = r(j)
  enddo

  bnrm = ddot(n,b,1_i4,b,1_i4)
  call denseAdiagprecon(mode,n,nloc,nlocptr,A,maxA,r,z)
!
!  Main loop
!
  iter = 0
  converged = .false.
  do while (iter.le.itmax.and..not.converged)
    iter = iter + 1
    call denseAdiagprecon(mode,n,nloc,nlocptr,A,maxA,rr,zz)
    bknum = ddot(n,z,1_i4,rr,1_i4)
    if (iter.eq.1) then
      do j = 1,n
        p(j) = z(j)
        pp(j) = zz(j)
      enddo
    else
      bk = bknum/bkden
      do j = 1,n
        p(j) = bk*p(j) + z(j)
        pp(j) = bk*pp(j) + zz(j)
      enddo
    endif
    bkden = bknum
    call denseAxV(mode,n,nloc,nlocptr,A,maxA,p,z,.false.)
    akden = ddot(n,z,1_i4,pp,1_i4)
    ak = bknum/akden
    call denseAxV(mode,n,nloc,nlocptr,A,maxA,pp,zz,.true.)
!
    call daxpy(n,ak,p,1_i4,x,1_i4)
    call daxpy(n,-ak,z,1_i4,r,1_i4)
    call daxpy(n,-ak,zz,1_i4,rr,1_i4)
!
    call denseAdiagprecon(mode,n,nloc,nlocptr,A,maxA,r,z)
!
    err = ddot(n,r,1_i4,r,1_i4)/bnrm
    err = sqrt(err)
    if (ioproc.and.index(keyword,'verb').ne.0) then
      write(ioout,'('' BCG iteration = '',i4,'' Error = '',f16.12)') iter,err
    endif
    converged = (err.lt.tol)
  enddo
!
!  Free workspace
!
  deallocate(zz)
  deallocate(z)
  deallocate(rr)
  deallocate(r)
  deallocate(pp)
  deallocate(p)
#ifdef TRACE
  call trace_out('dbcgsolve')
#endif

  end subroutine dbcgsolve

  subroutine denseAdiagprecon(mode,n,nloc,nlocptr,A,maxA,Vin,Vout)
!
!  Perform the precondition of a vector by the diagonal of a dense matrix
!
!  NB: If mode > 0 then code is specific to the purpose here in that 
!  the last row and column are assumed to be 1, except the diagonal 
!  element which is zero. This avoids having to explicitly store 
!  the constraint terms.
!
  use datatypes
  use parallel,      only : nprocs, procid
  use times,         only : tsum
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
  integer(i4),  intent(in)    :: mode
  integer(i4),  intent(in)    :: n
  integer(i4),  intent(in)    :: nloc
  integer(i4),  intent(in)    :: nlocptr(nloc)
  integer(i4),  intent(in)    :: maxA
  real(dp),     intent(in)    :: A(maxA,nloc)
  real(dp),     intent(in)    :: Vin(n)
  real(dp),     intent(out)   :: Vout(n)
!
  integer(i4)                 :: i
  integer(i4)                 :: io
  real(dp),     allocatable   :: Vloc(:)
  real(dp)                    :: g_cpu_time
  real(dp)                    :: tsuml
#ifdef TRACE
  call trace_in('denseAdiagprecon')
#endif
!
!  This code is specific to the form of the matrix being used here!!
!
  Vout(1:n) = 0.0_dp
  do i = 1,nloc
    io = nlocptr(i)
    Vout(io) = Vin(io)/A(io,i)
  enddo
  if (mode.gt.0.and.procid.eq.0) then
    Vout(n) = Vin(n)
  endif
  if (nprocs.gt.1) then
!
!  Globalization of Vout
!
    tsuml = g_cpu_time()
    allocate(Vloc(n))
    call sumall(Vout,Vloc,n,"denseAdiagprecon","Vout") 
    Vout(1:n) = Vloc(1:n)
    deallocate(Vloc)
    tsum = tsum + g_cpu_time() - tsuml
  endif
#ifdef TRACE
  call trace_out('denseAdiagprecon')
#endif

  end subroutine denseAdiagprecon

  subroutine denseAxV(mode,n,nloc,nlocptr,A,maxA,Vin,Vout,ltranspose)
!
!  Perform the multiplication of a dense matrix by a vector
!
!  NB: If mode > 0 then code is specific to the purpose here in that 
!  the last row and column are assumed to be 1, except the diagonal 
!  element which is zero. This avoids having to explicitly store 
!  the constraint terms.
! 
  use datatypes
  use parallel,      only : nprocs, procid
  use times,         only : tsum
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
  integer(i4),  intent(in)    :: mode
  integer(i4),  intent(in)    :: n
  integer(i4),  intent(in)    :: nloc
  integer(i4),  intent(in)    :: nlocptr(nloc)
  integer(i4),  intent(in)    :: maxA
  logical,      intent(in)    :: ltranspose
  real(dp),     intent(in)    :: A(maxA,*)
  real(dp),     intent(in)    :: Vin(n)
  real(dp),     intent(out)   :: Vout(n)

  integer(i4)                 :: i
  integer(i4)                 :: il
  integer(i4)                 :: j
  real(dp),     allocatable   :: Vloc(:)
!
  real(dp)                    :: g_cpu_time
  real(dp)                    :: ddot
  real(dp)                    :: tsuml
#ifdef TRACE
  call trace_in('denseAxV')
#endif
!
  if (ltranspose) then
    if (nprocs.eq.1) then
      do i = 1,n
        Vout(i) = 0.0_dp
        Vout(i) = ddot(n,A(i,1),maxA,Vin,1_i4)
      enddo
    else
      Vout(1:n) = 0.0_dp
      do il = 1,nloc
        i = nlocptr(il)
        do j = 1,n
          Vout(j) = Vout(j) + A(j,il)*Vin(i)
        enddo
      enddo
      if (mode.gt.0.and.procid.eq.0) then
        do j = 1,n
          Vout(j) = Vout(j) + A(j,nloc+1)*Vin(n)
        enddo
      endif
    endif
  else
    if (nprocs.eq.1) then
      do i = 1,n
        Vout(i) = ddot(n,A(1,i),1_i4,Vin,1_i4)
      enddo
    else
      Vout(1:n) = 0.0_dp
      do il = 1,nloc
        i = nlocptr(il)
        Vout(i) = ddot(n,A(1,il),1_i4,Vin,1_i4)
      enddo
      if (mode.gt.0.and.procid.eq.0) then
        Vout(n) = ddot(n,A(1,nloc+1),1_i4,Vin,1_i4)
      endif
    endif
  endif
  if (nprocs.gt.1) then
    allocate(Vloc(n))
!
!  Globalization of Vout
!
    tsuml = g_cpu_time()
    call sumall(Vout,Vloc,n,"denseAxV","Vout") 
    Vout(1:n) = Vloc(1:n)
    tsum = tsum + g_cpu_time() - tsuml
    deallocate(Vloc)
  endif
#ifdef TRACE
  call trace_out('denseAxV')
#endif

  end subroutine denseAxV
