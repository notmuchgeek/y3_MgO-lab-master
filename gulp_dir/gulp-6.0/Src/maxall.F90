  subroutine maxall(partial,gmax,count,caller,desc)
!
!  Uses MPI_allreduce to perform a global max for "count" real*8 items
!

!
!  Modules
!
  use parallel
  implicit none
!
!  Passed variables
!
  integer(i4)      :: count
  real(dp)         :: gmax(count)
  real(dp)         :: partial(count)
  character(len=*) :: caller, desc
!
!  Other variables
!
#ifdef MPI
  include 'mpif.h'
  integer*4           :: ierror,icode,lenstring,ic,i
  real*8, allocatable :: partial_mpi(:)
  real*8, allocatable :: gmax_mpi(:)
  character error_string(MPI_max_error_string)
!
  if (nprocs.gt.1) then
!
!  Copy values to MPI compatible variables
!
    allocate(partial_mpi(count))
    allocate(gmax_mpi(count))
    do i = 1,count
      partial_mpi(i) = partial(i)
    enddo
!
    call MPI_allreduce(partial_mpi,gmax_mpi,count,MPI_double_precision,MPI_max,MPI_comm_GULP,ierror)
!
!  Copy values back from MPI compatible variables
!
    do i = 1,count
      gmax(i) = gmax_mpi(i)
    enddo
    deallocate(gmax_mpi)
    deallocate(partial_mpi)

!  Check for errors
    if (ierror.ne.MPI_success) then
      write(*,"(5a)")  &
       "Call to MPI_allreduce failed in subroutine maxall,", &
       " called from subroutine ", caller, " for ", desc
      icode = ierror
      call MPI_error_string(icode,error_string,lenstring,ierror)
      write(*,"(a)") error_string 
      call MPI_abort(MPI_comm_GULP,1,ierror)
    endif
  elseif (nprocs.eq.1) then
    do ic = 1,count
      gmax(ic) = partial(ic)
    enddo
  endif
#else
  integer ic
  do ic = 1,count
    gmax(ic) = partial(ic)
  enddo
#endif
  return
  end
