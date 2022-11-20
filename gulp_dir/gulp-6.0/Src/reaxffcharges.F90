  subroutine reaxFFcharges
!
!  Calculates the charges for ReaxFF and outputs them
!
!   8/15 Created from reaxffmd.f90
!   6/17 nboatomptr passed to setreaxffQiter
!   1/18 Trace added
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
!  Julian Gale, CIC, Curtin University, January 2018
!
  use datatypes
  use control,        only : lreaxFFQ, literativeQ
  use current
  use reaxFFdata
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Local variables
!
  integer(i4)                                    :: i
  integer(i4)                                    :: j
  integer(i4)                                    :: nati
  integer(i4), dimension(:),   allocatable, save :: nbos
  integer(i4), dimension(:),   allocatable, save :: nbosptr
  integer(i4)                                    :: nboatom
  integer(i4), dimension(:),   allocatable, save :: nboatomptr
  integer(i4), dimension(:),   allocatable, save :: nboatomRptr
  integer(i4)                                    :: ntypi
  integer(i4)                                    :: status
  real(dp)                                       :: edummy
  real(dp)                                       :: g_cpu_time
  real(dp)                                       :: t1
  real(dp)                                       :: t2
  real(dp)                                       :: t3
  real(dp)                                       :: t4
#ifdef TRACE
  call trace_in('reaxffcharges')
#endif
!
  t1 = g_cpu_time()
!
  allocate(nboatomptr(numat),stat=status)
  if (status/=0) call outofmemory('reaxFFcharges','nboatomptr')
  allocate(nboatomRptr(numat),stat=status)
  if (status/=0) call outofmemory('reaxFFcharges','nboatomRptr')
!
!  Allocate memory 
!
  allocate(nbos(numat),stat=status)
  if (status/=0) call outofmemory('reaxFFcharges','nbos')
  allocate(nbosptr(numat),stat=status)
  if (status/=0) call outofmemory('reaxFFcharges','nbosptr')
!*************************************
!  Find cut-off radii for all atoms  *
!*************************************
  nbosptr(1:numat) = 0
  nboatom = 0
  do i = 1,numat
    nati = nat(i)
    ntypi = nftype(i)
!
!  Check repulsive terms and build atom pointers to species
!
    nbos(i) = 0
    do j = 1,nreaxFFspec
      if (nati.eq.natreaxFFspec(j).and.(ntypi.eq.ntypreaxFFspec(j).or.ntypreaxFFspec(j).eq.0)) then
        nbos(i) = nbos(i) + 1
        nbosptr(i) = j
      endif
    enddo
!
!  Check number of species for now
!
    if (nbos(i).gt.1) then
      call outerror('Multiple species per atom not yet allowed for in reaxFF',0_i4)
      call stopnow('reaxFFcharges')
    elseif (nbos(i).eq.1) then
      nboatom = nboatom + 1
      nboatomptr(i) = nboatom
      nboatomRptr(nboatom) = i
    endif
  enddo
!***************************
!  Compute ReaxFF charges  *
!***************************
  t3 = g_cpu_time()
  if (lreaxFFQ) then
    edummy = 0.0_dp
    if (literativeQ) then
      call setreaxffQiter(nboatom,nboatomptr,nboatomRptr,nbos,nbosptr,qreaxFF,edummy,.true.)
    else
      call setreaxffQ(nboatom,nboatomRptr,nbos,nbosptr,qreaxFF,edummy,.true.)
    endif
  endif
  t4 = g_cpu_time()
!
!  Free local memory
!
  deallocate(nbosptr,stat=status)
  if (status/=0) call deallocate_error('reaxFFcharges','nbosptr')
  deallocate(nbos,stat=status)
  if (status/=0) call deallocate_error('reaxFFcharges','nbos')
  deallocate(nboatomRptr,stat=status)
  if (status/=0) call deallocate_error('reaxFFcharges','nboatomRptr')
  deallocate(nboatomptr,stat=status)
  if (status/=0) call deallocate_error('reaxFFcharges','nboatomptr')
!
  t2 = g_cpu_time()
  treaxFF = treaxFF + t2 - t1 - (t4 - t3)
  teem    = teem + t4 - t3
#ifdef TRACE
  call trace_out('reaxffcharges')
#endif
!
  return
  end
