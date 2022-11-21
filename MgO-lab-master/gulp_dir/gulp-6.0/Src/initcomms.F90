  subroutine GULP_initcomms(MPI_comm_in)
!
!  Initialises MPI if necessary, finds own taskid and number of tasks
!
!  Modified to handle possible precision issues - pass local scalars
!  to get procid/nprocs for benefit of Cray
!
!  12/03 Silent option added in which ioproc is set to false
!   3/07 Chemshell modifications added and renamed 
!   3/09 MPI communicator changed to MPI_comm_GULP based on passed argument
!   6/09 MPI barrier changed to MPI_comm_world instead of MPI_comm_GULP
!  11/16 Blacs initialisation added
!   4/17 ChemShell interaction modified
!   6/17 Old ChemShell restored as a compile option
!   9/20 Blacs context now set to be MPI_Comm_GULP for ChemShell
!
!  Julian Gale, CIC, Curtin University, September 2020
!

!
!  Modules
!
  use iochannels
  use gulpchemsh
  use parallel

  implicit none
!
!  Passed variables
!
  integer*4,   intent(in)   :: MPI_comm_in
!
#ifdef MPI
  include 'mpif.h'
  integer ierr,lprocid,lnprocs
  logical lmpiinit

#ifdef OLDCS
  if (ichemsh_qm .lt. 0) then
#else
  if (ichemsh_link .eq. 0) then
#endif
!
!  Non-ChemShell case - initialise MPI
!
    call MPI_Initialized(lmpiinit, ierr)
    if (.not. lmpiinit) then
       call MPI_init(ierr)
!
!  Set communicator for MPI based on MPI_comm_world
!
       MPI_comm_GULP = MPI_comm_world
    end if
  else
!
!  ChemShell case - MPI is assumed to be already running and communicator set using argument
!
    MPI_comm_GULP = MPI_comm_in
  endif

  call MPI_comm_rank(MPI_comm_GULP,lprocid,ierr)
  call MPI_comm_size(MPI_comm_GULP,lnprocs,ierr)

  procid  = lprocid
  nprocs  = lnprocs
!
!  Initialise Blacs for use by pblas/scalapack
!
  iBlacsContext = MPI_comm_GULP
  call blacs_gridinit( iBlacsContext, 'C', 1, lnprocs)
#else
  procid  = 0
  nprocs  = 1
#endif
  if (lsilent) then
    ioproc = .false.
  else
    ioproc = (procid.eq.0)
  endif
#ifdef MPI
  call MPI_barrier(MPI_comm_GULP,ierr)
#endif

  return
  end
