  subroutine pdiaggd(mcv,mcvloc,maxd2,derv2,maxeigr,eigr,freq,fscale,lvectors,lprint,ifail)
!
!  Calculate eigenvalues and eigenvectors of the dynamical matrix
!  and subsequent calculation of properties of eigenvectors.
!
!  Gamma point only version with distributed memory
!
!  11/16 Created from pdiaggd
!  12/16 maxeigr added as a separate argument to maxd2
!   2/17 Corrected so that idesce is not set up when eigenvectors are not needed
!        otherwise there will be an error due to the dimensions of the array
!   8/17 Routine wrapped in ifdef since it won't work or be called without MPI
!   2/18 Trace added
!
!  On entry :
!
!  mcv      = number of modes (global)
!  mcvloc   = number of modes (local)
!  maxd2    = left-hand dimension of derv2
!  maxeigr  = left-hand dimension of eigr
!  derv2    = mass-weighted dynamical matrix
!  fscale   = scale factor for frequencies to convert to wavenumbers
!  lvectors = if .true. then calculate eigenvectors
!  lprint   = if .true. print warnings
!
!  On exit :
!
!  eigr     = eigenvectors of dynamical matrix (if lvectors is true)
!  freq     = frequencies of vibration in wavenumbers
!  ifail    = flag indicating success or failure
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
!  Copyright Curtin University 2018
!
!  Julian Gale, CIC, Curtin University, Feb 2018
!
  use g_constants
  use current
  use element
  use parallel
  use species
  use times
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4),  intent(out)                    :: ifail
  integer(i4),  intent(in)                     :: maxd2
  integer(i4),  intent(in)                     :: maxeigr
  integer(i4),  intent(in)                     :: mcv
  integer(i4),  intent(in)                     :: mcvloc
  logical,      intent(in)                     :: lprint
  logical,      intent(in)                     :: lvectors
  real(dp),     intent(in)                     :: derv2(maxd2,*)
  real(dp),     intent(out)                    :: eigr(maxeigr,*)
  real(dp),     intent(in)                     :: fscale
  real(dp),     intent(out)                    :: freq(*)
#ifdef MPI
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: j
  integer(i4)                                  :: status
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: root
  real(dp)                                     :: t1
  real(dp)                                     :: t2
  real(dp),    dimension(:,:), allocatable     :: wrk2
  real(dp),    dimension(:),   allocatable     :: wrk
!
!  Blacs / Scalapack integers
!
  integer                                      :: idescd(9)
  integer                                      :: idesce(9)
  integer                                      :: ifails
  integer                                      :: ld
  integer                                      :: lwrk
  integer                                      :: n
  integer                                      :: nb
  integer                                      :: nloc
  integer                                      :: nsize
#ifdef TRACE
  call trace_in('pdiaggd')
#endif
!
  t1 = g_cpu_time()
!
  n = mcv
  nloc = mcvloc
  nb = 3*nblocksize
!
!  Set up Blacs descriptors for arrays
!
  ld = maxd2
  call descinit(idescd, n, n, nb, nb, 0, 0, iBlacsContext, ld, ifails)
!
  if (ifails.ne.0) then
    call outerror('initialisation in descinit failed - idescd ',0_i4)
    call stopnow('pdiaggd')
  endif
!
!  Only need eigr BLACS descriptor if vectors are to be computed
!
  if (lvectors) then
    ld = maxeigr
    call descinit(idesce, n, n, nb, nb, 0, 0, iBlacsContext, ld, ifails)
!
    if (ifails.ne.0) then
      call outerror('initialisation in descinit failed - idesce ',0_i4)
      call stopnow('pdiaggd')
    endif
  endif
!
!  Find size of work space array
!
  if (lvectors) then
    nsize = max((nb*(nb-1))/2,(n+nloc)*nb) + nb*nb
    lwrk = 5*n + n*max(1,nloc) + max(2*n,nsize) + 1
  else
    lwrk = 5*n + max(nb*(n+1),3*nb) + 1
  endif
!
!  Allocate work space
!
  allocate(wrk(lwrk),stat=status)
  if (status/=0) call outofmemory('pdiaggd','wrk')
  allocate(wrk2(maxd2,mcvloc),stat=status)
  if (status/=0) call outofmemory('pdiaggd','wrk2')
!
!  Make copy of derv2 to avoid overwriting
!
  do i = 1,mcvloc
    do j = 1,mcv
      wrk2(j,i) = derv2(j,i)
    enddo
  enddo
!
!  Call Scalapack eigensolver
!
  if (lvectors) then
    call pdsyev( 'V', 'U', n, wrk2, 1, 1, idescd, freq, eigr, 1, 1, idesce, wrk, lwrk, ifails )
  else
    call pdsyev( 'N', 'U', n, wrk2, 1, 1, idescd, freq, wrk2, 1, 1, idescd, wrk, lwrk, ifails )
  endif
  ifail = ifails
!
!  Free workspace
!
  deallocate(wrk2,stat=status)
  if (status/=0) call deallocate_error('pdiaggd','wrk2')
  deallocate(wrk,stat=status)
  if (status/=0) call deallocate_error('pdiaggd','wrk')
!
  t2 = g_cpu_time()
  tdiag = tdiag + t2 - t1
!
!  Trap diagonalisation failure
!
  if (ifail.gt.0) then
    call outerror('diagonalisation of dynamical matrix failed',0_i4)
    call stopnow('pdiaggd')
  endif
!
!  Convert frequency units - imaginary freqs denoted by negative no.
!
  do i = 1,mcv
    root = freq(i)
    if (root.ge.0.0_dp) then
      freq(i) = sqrt(root)*fscale
    else
      root = abs(root)
      freq(i) = - sqrt(root)*fscale
    endif
  enddo
#else
  call outerror('pdiaggd called when not compiled with MPI',0_i4)
  call stopnow('pdiaggd')
#endif
#ifdef TRACE
  call trace_out('pdiaggd')
#endif
!
  return
  end
