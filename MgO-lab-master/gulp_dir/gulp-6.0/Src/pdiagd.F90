  subroutine pdiagd(mcv,mcvloc,maxd2,derv2,dervi,maxeigc,eigc,freq,fscale,lvectors,lprint,ifail)
!
!  Calculate eigenvalues and eigenvectors of the dynamical matrix
!  and subsequent calculation of properties of eigenvectors.
!
!  General k point version version with distributed memory
!
!  12/16 Created from pdiaggd
!  12/16 maxeigc added as a separate argument to maxd2
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
!  maxeigc  = left-hand dimension of derv2 
!  derv2    = mass-weighted dynamical matrix (real part)
!  dervi    = mass-weighted dynamical matrix (complex part)
!  fscale   = scale factor for frequencies to convert to wavenumbers
!  lvectors = if .true. then calculate eigenvectors
!  lprint   = if .true. print warnings
!
!  On exit :
!
!  eigc     = eigenvectors of dynamical matrix (if lvectors is true)
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
  integer(i4),  intent(in)                     :: maxeigc
  integer(i4),  intent(in)                     :: mcv
  integer(i4),  intent(in)                     :: mcvloc
  logical,      intent(in)                     :: lprint
  logical,      intent(in)                     :: lvectors
  real(dp),     intent(in)                     :: derv2(maxd2,*)
  real(dp),     intent(in)                     :: dervi(maxd2,*)
  complex(dpc), intent(out)                    :: eigc(maxeigc,*)
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
  complex(dpc), dimension(:,:), allocatable    :: cwrk2
  complex(dpc), dimension(:),   allocatable    :: rwrk
  complex(dpc), dimension(:),   allocatable    :: wrk
!
!  Blacs / Scalapack integers
!
  integer                                      :: idescd(9)
  integer                                      :: idesce(9)
  integer                                      :: ifails
  integer                                      :: ld
  integer                                      :: lrwrk
  integer                                      :: lwrk
  integer                                      :: n
  integer                                      :: nb
  integer                                      :: np0
  integer                                      :: nq0
  integer                                      :: nloc
  integer                                      :: numroc
#ifdef TRACE
  call trace_in('pdiagd')
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
    call stopnow('pdiagd')
  endif
!
  if (lvectors) then
    ld = maxeigc
    call descinit(idesce, n, n, nb, nb, 0, 0, iBlacsContext, ld, ifails)
!
    if (ifails.ne.0) then
      call outerror('initialisation in descinit failed - idesce ',0_i4)
      call stopnow('pdiagd')
    endif
  endif
!
!  Find size of work space array
!
  np0 = numroc( n, nb, 0, 0, 1)
  nq0 = numroc( max( n, nb, 2 ), nb, 0, 0, nprocs )
  if (lvectors) then
    lwrk = 3*n + n*n + nb*(np0 + nq0 + nb)
    lrwrk = 4*n
  else
    lwrk = 3*n + max(nb*(np0+1),3) + 1
    lrwrk = 2*n
  endif
!
!  Allocate work space
!
  allocate(wrk(lwrk),stat=status)
  if (status/=0) call outofmemory('pdiagd','wrk')
  allocate(rwrk(lrwrk),stat=status)
  if (status/=0) call outofmemory('pdiagd','rwrk')
  allocate(cwrk2(maxd2,mcvloc),stat=status)
  if (status/=0) call outofmemory('pdiagd','cwrk2')
!
!  Make copy of derv2 to avoid overwriting
!
  do i = 1,mcvloc
    do j = 1,mcv
      cwrk2(j,i) = dcmplx(derv2(j,i),dervi(j,i))
    enddo
  enddo
!
!  Call Scalapack eigensolver
!
  if (lvectors) then
    call pzheev( 'V', 'U', n, cwrk2, 1, 1, idescd, freq, eigc, 1, 1, idesce, wrk, lwrk, rwrk, lrwrk, ifails )
  else
    call pzheev( 'N', 'U', n, cwrk2, 1, 1, idescd, freq, cwrk2, 1, 1, idescd, wrk, lwrk, rwrk, lrwrk, ifails )
  endif
  ifail = ifails
!
!  Free workspace
!
  deallocate(cwrk2,stat=status)
  if (status/=0) call deallocate_error('pdiagd','cwrk2')
  deallocate(rwrk,stat=status)
  if (status/=0) call deallocate_error('pdiagd','rwrk')
  deallocate(wrk,stat=status)
  if (status/=0) call deallocate_error('pdiagd','wrk')
!
  t2 = g_cpu_time()
  tdiag = tdiag + t2 - t1
!
!  Trap diagonalisation failure
!
  if (ifail.gt.0) then
    call outerror('diagonalisation of dynamical matrix failed',0_i4)
    call stopnow('pdiagd')
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
  call outerror('pdiagd called when not compiled with MPI',0_i4)
  call stopnow('pdiagd')
#endif
#ifdef TRACE
  call trace_out('pdiagd')
#endif
!
  return
  end
