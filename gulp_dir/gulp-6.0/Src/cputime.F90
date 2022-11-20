  function g_cpu_time()
  use datatypes
  implicit none
!
!  Passed variables
!
  real(dp)           :: g_cpu_time
#ifdef MPI
!***************
!  MPI timing  *
!***************
  include 'mpif.h'
  integer(i4),  save :: first = 1
  real(dp),     save :: start

  if (first.eq.1) then
    first = 0
    start = MPI_wtime()
    g_cpu_time = 0.0_dp
  else
    g_cpu_time = MPI_wtime() - start
  endif
#else
!***************
!  F90 timing  *
!***************
  call cpu_time(g_cpu_time)
#endif
  return
  end
!**************************************************************
!  Dummy routines for flush - on VMS, HP and IBM call fflush  *
!**************************************************************
  subroutine gflush(iout)
#ifdef NAG
!For NAG compiler only:  
use f90_unix_io,only:flush
#endif
  use datatypes
  implicit none
  integer(i4) :: iout
#ifdef FFLUSH
  call fflush(iout)
#elif defined FLUSHWITHERRORFLAG
  integer(i4) :: ierr
  call flush(iout,ierr)
#elif defined FLUSH
  call flush(iout)
#endif
  return
  end
