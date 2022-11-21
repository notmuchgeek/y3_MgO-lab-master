  subroutine imaxall(ipartial,imax,count,caller,desc)
!
!  Uses MPI_allreduce to perform a global max for "count" integer*4 items
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
  integer(i4)      :: imax(count)
  integer(i4)      :: ipartial(count)
  character(len=*) :: caller, desc
!
!  Other variables
!
#ifdef MPI
  include 'mpif.h'
  integer*4              :: ierror,icode,lenstring,ic,i
  integer*4, allocatable :: ipartial_mpi(:)
  integer*4, allocatable :: imax_mpi(:)
  character error_string(MPI_max_error_string)
!
  if (nprocs.gt.1) then
!
!  Copy values to MPI compatible variables
!
    allocate(ipartial_mpi(count))
    allocate(imax_mpi(count))
    do i = 1,count
      ipartial_mpi(i) = ipartial(i)
    enddo
!
    call MPI_allreduce(ipartial_mpi,imax_mpi,count,MPI_integer,MPI_max,MPI_comm_GULP,ierror)
!
!  Copy values back from MPI compatible variables
!
    do i = 1,count
      imax(i) = imax_mpi(i)
    enddo
    deallocate(imax_mpi)
    deallocate(ipartial_mpi)

!  Check for errors
    if (ierror.ne.MPI_success) then
      write(*,"(5a)")  &
       "Call to MPI_allreduce failed in subroutine imaxall,", &
       " called from subroutine ", caller, " for ", desc
      icode = ierror
      call MPI_error_string(icode,error_string,lenstring,ierror)
      write(*,"(a)") error_string 
      call MPI_abort(MPI_comm_GULP,1,ierror)
    endif
  elseif (nprocs.eq.1) then
    do ic = 1,count
      imax(ic) = ipartial(ic)
    enddo
  endif
#else
  integer ic
  do ic = 1,count
    imax(ic) = ipartial(ic)
  enddo
#endif
  return
  end
