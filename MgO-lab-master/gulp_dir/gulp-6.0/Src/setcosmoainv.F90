  subroutine setcosmoainv
!
!  Calculates the inverse of the COSMO A matrix between points on the SAS
!
!   5/03 Separated from setcosmoamat
!  11/04 Timing of inversions added to global inversion total
!  12/08 Migrated to version 3.5 and converted to f90 format
!  10/14 Trap for npts = 0 added
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
!  Copyright Curtin University 2014
!
!  Julian Gale, CIC, Curtin University, October 2014
!
  use cosmic
  use times,   only : tmati
  implicit none
!
!  Local variables
!
  integer(i4)                                    :: info
  integer(i4), dimension(:),   allocatable, save :: ipivot
  integer(i4)                                    :: status
  real(dp),    dimension(:),   allocatable, save :: wrk
  real(dp)                                       :: g_cpu_time
  real(dp)                                       :: t1
  real(dp)                                       :: t2
!
!  If npts is zero then return
!
  if (npts.eq.0) return
!
  t1 = g_cpu_time()
!
!  Allocate memory only needed for matrix operations
!
  allocate(ipivot(npts),stat=status)
  if (status/=0) call outofmemory('setcosmoainv','ipivot')
  allocate(wrk(3*maxnpts),stat=status)
  if (status/=0) call outofmemory('setcosmoainv','wrk')
!
!  Invert COSMO A-matrix
!
  call dsptrf('U',npts,cosmoA,ipivot,info)
  call dsptri('U',npts,cosmoA,ipivot,wrk,info)
!
!  Free local memory
!
  deallocate(wrk,stat=status)
  if (status/=0) call deallocate_error('setcosmoainv','wrk')
  deallocate(ipivot,stat=status)
  if (status/=0) call deallocate_error('setcosmoainv','ipivot')
!
  t2 = g_cpu_time()
  tmati = tmati + t2 - t1
!
  return
  end
