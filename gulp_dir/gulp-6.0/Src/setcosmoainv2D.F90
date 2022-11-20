  subroutine setcosmoainv2D
!
!  Calculates the inverse of the COSMO A matrix between points on the SAS
!  2D matrix format for cosmoA version.
!
!   4/17 Created from setcosmoainv
!   4/17 Parallelisation added
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
!  Copyright Curtin University 2017
!
!  Julian Gale, CIC, Curtin University, April 2017
!
  use cosmic
  use parallel
  use times,   only : tmati
  implicit none
!
!  Local variables
!
  integer(i4)                                    :: i
  integer(i4)                                    :: info
  integer(i4), dimension(:),   allocatable, save :: ipivot
  integer(i4)                                    :: j
  integer(i4)                                    :: lwrk
  integer(i4)                                    :: nblock
  integer(i4)                                    :: status
  real(dp),    dimension(:),   allocatable, save :: wrk
  real(dp)                                       :: g_cpu_time
  real(dp)                                       :: t1
  real(dp)                                       :: t2
!
!  Functions
!
  integer(i4)                                    :: ilaenv
#ifdef MPI
  integer(i4), dimension(:),   allocatable, save :: iwrk
!
!  Local variables in Scalapack/Blacs integer precision
!
  integer                                        :: ifails
  integer                                        :: idesA(9)
  integer                                        :: ldm
  integer                                        :: liwrk
  integer                                        :: nb
  integer                                        :: np
  integer                                        :: nv
#endif
!
!  If npts is zero then return
!
  if (npts.eq.0) return
!
  t1 = g_cpu_time()
#ifdef MPI
  if (nprocs.gt.1) then
    allocate(ipivot(2*npts),stat=status)
    if (status/=0) call outofmemory('setcosmoainv2D','ipivot')
!
!  Set local block size
!
    nb = nblocksizesas
    np = nprocs
    nv = npts
    ifails = 0
!
!  Set up Blacs descriptors
!
    ldm = maxcosmoA
    call descinit( idesA, nv, nv, nb, nb, 0, 0, iBlacsContext, ldm, ifails )
!
    if (ifails.ne.0) then
      call outerror('initialisation in descinit failed',0_i4)
      call stopnow('setcosmoainv2D')
    endif
!
!  Factorise matrix using Scalapack
!
    call pdgetrf(nv,nv,cosmoA,1,1,idesA,ipivot,ifails)
!
    if (ifails.ne.0) then
      call outerror('inversion of cosmoA has failed in pdgetrf',0_i4)
      call stopnow('setcosmoainv2D')
    endif
!***********************
!  Complete inversion  *
!***********************
!
!  Initial dummy allocation of workspace
!
    lwrk = 1
    liwrk = 1
    allocate(iwrk(liwrk),stat=status)
    if (status/=0) call outofmemory('setcosmoainv2D','iwrk')
    allocate(wrk(lwrk),stat=status)
    if (status/=0) call outofmemory('setcosmoainv2D','wrk')
!
!  Query to find workspace needed
!
    lwrk = -1
    liwrk = -1
    ifails = 0
    call pdgetri(nv,cosmoA,1,1,idesA,ipivot,wrk,lwrk,iwrk,liwrk,ifails)
!
    lwrk = nint(wrk(1))
    liwrk = iwrk(1)
!
!  Reallocate workspace to size needed
!
    deallocate(wrk,stat=status)
    if (status/=0) call deallocate_error('setcosmoainv2D','wrk')
    deallocate(iwrk,stat=status)
    if (status/=0) call deallocate_error('setcosmoainv2D','iwrk')
    allocate(iwrk(liwrk),stat=status)
    if (status/=0) call outofmemory('setcosmoainv2D','iwrk')
    allocate(wrk(lwrk),stat=status)
    if (status/=0) call outofmemory('setcosmoainv2D','wrk')
!
!  Form inverse
!
    call pdgetri(nv,cosmoA,1,1,idesA,ipivot,wrk,lwrk,iwrk,liwrk,ifails)
!
!
!  Deallocate workspace
!
    deallocate(wrk,stat=status)
    if (status/=0) call deallocate_error('setcosmoainv2D','wrk')
    deallocate(iwrk,stat=status)
    if (status/=0) call deallocate_error('setcosmoainv2D','iwrk')
    deallocate(ipivot,stat=status)
    if (status/=0) call deallocate_error('setcosmoainv2D','ipivot')
  else
#endif
!
!  Compute optimal workspace
!
    nblock = ilaenv(1,'DSYTRF','U',npts,-1,-1,-1)
    lwrk = nblock*npts
!
!  Allocate memory only needed for matrix operations
!
    allocate(ipivot(npts),stat=status)
    if (status/=0) call outofmemory('setcosmoainv2D','ipivot')
    allocate(wrk(lwrk),stat=status)
    if (status/=0) call outofmemory('setcosmoainv2D','wrk')
!
!  Invert COSMO A-matrix
!
    call dsytrf('U',npts,cosmoA,maxcosmoA,ipivot,wrk,lwrk,info)
    if (info.ne.0) then
      call outerror('inversion of cosmoA has failed in dsytrf',0_i4)
      call stopnow('setcosmoainv2D')
    endif
    call dsytri('U',npts,cosmoA,maxcosmoA,ipivot,wrk,info)
!
!  Copy upper triangular part to lower triangle
!
    do i = 2,npts
      do j = 1,i-1
        cosmoA(i,j) = cosmoA(j,i)
      enddo
    enddo
!
!  Free local memory
!
    deallocate(wrk,stat=status)
    if (status/=0) call deallocate_error('setcosmoainv2D','wrk')
    deallocate(ipivot,stat=status)
    if (status/=0) call deallocate_error('setcosmoainv2D','ipivot')
#ifdef MPI
  endif
#endif
!
  t2 = g_cpu_time()
  tmati = tmati + t2 - t1
!
  return
  end
