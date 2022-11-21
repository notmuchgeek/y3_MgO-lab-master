  subroutine sumall(partial,sum,count,caller,desc)
!
!  Uses MPI_allreduce to perform a global sum for "count" real*8 items
!

!
!  Modules
!
  use parallel
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4)      :: count
  real(dp)         :: sum(count)
  real(dp)         :: partial(count)
  character(len=*) :: caller, desc
!
!  Other variables
!
#ifdef MPI
  include 'mpif.h'
  integer*4           :: ierror,icode,lenstring,ic,i
  real*8, allocatable :: partial_mpi(:)
  real*8, allocatable :: sum_mpi(:)
  character error_string(MPI_max_error_string)
#else
  integer ic
#endif
!
#ifdef TRACE
  call trace_in('sumall')
#endif
#ifdef MPI
  if (nprocs.gt.1) then
!
!  Copy values to MPI compatible variables
!
    allocate(partial_mpi(count))
    allocate(sum_mpi(count))
    do i = 1,count
      partial_mpi(i) = partial(i)
    enddo
!
    call MPI_allreduce(partial_mpi,sum_mpi,count,MPI_double_precision,MPI_sum,MPI_comm_GULP,ierror)
!
!  Copy values back from MPI compatible variables
!
    do i = 1,count
      sum(i) = sum_mpi(i)
    enddo
    deallocate(sum_mpi)
    deallocate(partial_mpi)

!  Check for errors
    if (ierror.ne.MPI_success) then
      write(*,"(5a)")  &
       "Call to MPI_allreduce failed in subroutine sumall,", &
       " called from subroutine ", caller, " for ", desc
      icode = ierror
      call MPI_error_string(icode,error_string,lenstring,ierror)
      write(*,"(a)") error_string 
      call MPI_abort(MPI_comm_GULP,1,ierror)
    endif
  elseif (nprocs.eq.1) then
    do ic = 1,count
      sum(ic) = partial(ic)
    enddo
  endif
#else
  do ic = 1,count
    sum(ic) = partial(ic)
  enddo
#endif
#ifdef TRACE
  call trace_out('sumall')
#endif
  return
  end

  subroutine isumall(ipartial,isum,count,caller,desc)
!
!  Uses MPI_allreduce to perform a global sum for "count" integer*4 items
!

!
!  Modules
!
  use parallel
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4)      :: count
  integer(i4)      :: isum(count)
  integer(i4)      :: ipartial(count)
  character(len=*) :: caller, desc
!
!  Other variables
!
#ifdef MPI
  include 'mpif.h'
  integer*4              :: ierror,icode,lenstring,ic,i
  integer*4, allocatable :: ipartial_mpi(:)
  integer*4, allocatable :: isum_mpi(:)
  character error_string(MPI_max_error_string)
#else
  integer ic
#endif
!
#ifdef TRACE
  call trace_in('isumall')
#endif
#ifdef MPI
  if (nprocs.gt.1) then
!
!  Copy values to MPI compatible variables
!
    allocate(ipartial_mpi(count))
    allocate(isum_mpi(count))
    do i = 1,count
      ipartial_mpi(i) = ipartial(i)
    enddo
!
    call MPI_allreduce(ipartial_mpi,isum_mpi,count,MPI_integer,MPI_sum,MPI_comm_GULP,ierror)
!
!  Copy values back from MPI compatible variables
!
    do i = 1,count
      isum(i) = isum_mpi(i)
    enddo
    deallocate(isum_mpi)
    deallocate(ipartial_mpi)

!  Check for errors
    if (ierror.ne.MPI_success) then
      write(*,"(5a)")  &
       "Call to MPI_allreduce failed in subroutine isumall,", &
       " called from subroutine ", caller, " for ", desc
      icode = ierror
      call MPI_error_string(icode,error_string,lenstring,ierror)
      write(*,"(a)") error_string 
      call MPI_abort(MPI_comm_GULP,1,ierror)
    endif
  elseif (nprocs.eq.1) then
    do ic = 1,count
      isum(ic) = ipartial(ic)
    enddo
  endif
#else
  do ic = 1,count
    isum(ic) = ipartial(ic)
  enddo
#endif
#ifdef TRACE
  call trace_out('isumall')
#endif
  return
  end

  subroutine landall(lpartial,land,count,caller,desc)
!
!  Uses MPI_allreduce to perform a global logical "and" for "count" logical items
!

!
!  Modules
!
  use parallel
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4)      :: count
  logical          :: land(count)
  logical          :: lpartial(count)
  character(len=*) :: caller, desc
!
!  Other variables
!
#ifdef MPI
  include 'mpif.h'
  integer*4              :: ierror,icode,lenstring,ic,i
  character error_string(MPI_max_error_string)
#else
  integer ic
#endif
#ifdef TRACE
  call trace_in('landall')
#endif
!
#ifdef MPI
  if (nprocs.gt.1) then
!
    call MPI_allreduce(lpartial,land,count,MPI_logical,MPI_land,MPI_comm_GULP,ierror)

!  Check for errors
    if (ierror.ne.MPI_success) then
      write(*,"(5a)")  &
       "Call to MPI_allreduce failed in subroutine landall,", &
       " called from subroutine ", caller, " for ", desc
      icode = ierror
      call MPI_error_string(icode,error_string,lenstring,ierror)
      write(*,"(a)") error_string 
      call MPI_abort(MPI_comm_GULP,1,ierror)
    endif
  elseif (nprocs.eq.1) then
    do ic = 1,count
      land(ic) = lpartial(ic)
    enddo
  endif
#else
  do ic = 1,count
    land(ic) = lpartial(ic)
  enddo
#endif
#ifdef TRACE
  call trace_out('landall')
#endif
  return
  end

  subroutine csumall(partial,sum,count,caller,desc)
!
!  Uses MPI_allreduce to perform a global sum for "count" complex*16 items
!

!
!  Modules
!
  use parallel
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4)      :: count
  complex(dpc)     :: sum(count)
  complex(dpc)     :: partial(count)
  character(len=*) :: caller, desc
!
!  Other variables
!
#ifdef MPI
  include 'mpif.h'
  integer*4           :: ierror,icode,lenstring,ic,i
  complex*16, allocatable :: partial_mpi(:)
  complex*16, allocatable :: sum_mpi(:)
  character error_string(MPI_max_error_string)
#else
  integer ic
#endif
!
#ifdef TRACE
  call trace_in('csumall')
#endif
#ifdef MPI
  if (nprocs.gt.1) then
!
!  Copy values to MPI compatible variables
!
    allocate(partial_mpi(count))
    allocate(sum_mpi(count))
    do i = 1,count
      partial_mpi(i) = partial(i)
    enddo
!
    call MPI_allreduce(partial_mpi,sum_mpi,count,MPI_double_complex,MPI_sum,MPI_comm_GULP,ierror)
!
!  Copy values back from MPI compatible variables
!
    do i = 1,count
      sum(i) = sum_mpi(i)
    enddo
    deallocate(sum_mpi)
    deallocate(partial_mpi)

!  Check for errors
    if (ierror.ne.MPI_success) then
      write(*,"(5a)")  &
       "Call to MPI_allreduce failed in subroutine csumall,", &
       " called from subroutine ", caller, " for ", desc
      icode = ierror
      call MPI_error_string(icode,error_string,lenstring,ierror)
      write(*,"(a)") error_string 
      call MPI_abort(MPI_comm_GULP,1,ierror)
    endif
  elseif (nprocs.eq.1) then
    do ic = 1,count
      sum(ic) = partial(ic)
    enddo
  endif
#else
  do ic = 1,count
    sum(ic) = partial(ic)
  enddo
#endif
#ifdef TRACE
  call trace_out('csumall')
#endif
  return
  end
!
  subroutine sumone(partial,sum,caller,desc)
!
!  Uses MPI_allreduce to perform a global sum for "count" real*8 items
!

!
!  Modules
!
  use parallel
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  real(dp)         :: sum
  real(dp)         :: partial
  character(len=*) :: caller, desc
!
!  Other variables
!
#ifdef MPI
  include 'mpif.h'
  integer*4           :: ierror,icode,lenstring
  real*8              :: partial_mpi
  real*8              :: sum_mpi
  character error_string(MPI_max_error_string)
#endif
!
#ifdef TRACE
  call trace_in('sumall')
#endif
#ifdef MPI
  if (nprocs.gt.1) then
!
!  Copy values to MPI compatible variables
!
    partial_mpi = partial
!
    call MPI_allreduce(partial_mpi,sum_mpi,1,MPI_double_precision,MPI_sum,MPI_comm_GULP,ierror)
!
!  Copy values back from MPI compatible variables
!
    sum = sum_mpi

!  Check for errors
    if (ierror.ne.MPI_success) then
      write(*,"(5a)")  &
       "Call to MPI_allreduce failed in subroutine sumall,", &
       " called from subroutine ", caller, " for ", desc
      icode = ierror
      call MPI_error_string(icode,error_string,lenstring,ierror)
      write(*,"(a)") error_string 
      call MPI_abort(MPI_comm_GULP,1,ierror)
    endif
  elseif (nprocs.eq.1) then
    sum = partial
  endif
#else
  sum = partial
#endif
#ifdef TRACE
  call trace_out('sumall')
#endif
  return
  end
