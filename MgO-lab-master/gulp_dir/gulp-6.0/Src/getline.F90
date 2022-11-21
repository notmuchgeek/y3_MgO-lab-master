  subroutine gulp_getline(ioch,line,lend)
!
!  Reads in one line of input file and broadcasts it to all processors as necessary
!
!   2/17 I/O channel is now an argument so that the routine can be
!        used for multiple purposes
!
!  Julian Gale, CIC, Curtin University, February 2017
!

!
!  Modules
!
  use gulpinput
  use iochannels
  use parallel
  implicit none
!
!  Passed variables
!
  character(len=maxlinelength), intent(out) :: line
  integer(i4),                  intent(in)  :: ioch
  logical,                      intent(out) :: lend
!
!  Local variables
!
#ifdef MPI
  include 'mpif.h'
  integer                                   :: ierr
  integer                                   :: maxlinelen

  lend = .false.

!
!  Copy to MPI data type compatible variable
!
  maxlinelen = int(maxlinelength)

  line = ' '
  if (procid.eq.0) then
    read(ioch,'(a)',end=100) line
    call MPI_bcast(lend,1,MPI_logical,0,MPI_comm_GULP,ierr)
#ifdef USE_CVF
    call MPI_bcast_str(line,maxlinelen,MPI_character,0,MPI_comm_GULP,ierr)
#else
    call MPI_bcast(line,maxlinelen,MPI_character,0,MPI_comm_GULP,ierr)
#endif
    return
100   continue
!  End of input, need to broadcast this to all processors
    lend = .true.
    call MPI_bcast(lend,1,MPI_logical,0,MPI_comm_GULP,ierr)
  else
    call MPI_bcast(lend,1,MPI_logical,0,MPI_comm_GULP,ierr)
    if(lend) return
#ifdef USE_CVF
    call MPI_bcast_str(line,maxlinelen,MPI_character,0,MPI_comm_GULP,ierr)
#else
    call MPI_bcast(line,maxlinelen,MPI_character,0,MPI_comm_GULP,ierr)
#endif
  endif
#else
  lend = .false.
  line = ' '
  read(ioch,'(a)',end=100) line
  return
100   lend = .true.
#endif
  end
