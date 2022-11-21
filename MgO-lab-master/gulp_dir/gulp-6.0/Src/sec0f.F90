  subroutine sec0f
!
!  Generate symmetry adapted cluster second derivative matrix.
!
!  Freezing now added.
!
!   8/95 Modifications added to avoid using transformation
!        matrix unless constraints are present as this is
!        much slower.
!   7/02 Referencing of iopt corrected
!   8/02 Error in iopt for breathing shell corrected
!   7/05 Deallocation cleaned up
!   3/17 fix_atom option added
!   2/18 Trace added
!   3/19 iopt replaced by ioptindex and iopttype
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
!  Copyright Curtin University 2019
!
!  Julian Gale, CIC, Curtin University, March 2019
!
  use control
  use current
  use derivatives
  use iochannels
  use optimisation
  use parallel
  use times
#ifdef TRACE
  use trace,       only : trace_in, trace_out
#endif
  use transform
  implicit none
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ii
  integer(i4), dimension(:), allocatable       :: ioptfindex
  integer(i4)                                  :: j
  integer(i4)                                  :: jj
  integer(i4)                                  :: k
  integer(i4)                                  :: n3f
  integer(i4)                                  :: nff
  integer(i4), dimension(:), allocatable       :: nptr
  integer(i4)                                  :: status
  logical                                      :: lsdebug
  logical                                      :: ltmat
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: t1
  real(dp)                                     :: t2
  real(dp),    dimension(:), allocatable       :: tmp2
  real(dp)                                     :: tr1
#ifdef TRACE
  call trace_in('sec0f')
#endif
!
  t1 = g_cpu_time()
  lsdebug = (index(keyword,'derv2').ne.0)
!
!  Free local memory
!
  allocate(ioptfindex(nvar),stat=status)
  if (status/=0) call outofmemory('sec0f','ioptfindex')
  allocate(nptr(numat),stat=status)
  if (status/=0) call outofmemory('sec0f','nptr')
!
  nff = 0
  do i = 1,numat
    if (lopf(i)) then
      nff = nff + 1
      nptr(i) = nff
    else
      nptr(i) = 0
    endif
  enddo
  n3f = 3*nff
  if (nbsm.gt.0) n3f = n3f + nff
!
!  Convert iopt values to new reference system of atoms
!
  do i = 1,nvar
    ii = ioptindex(i)
    if (iopttype(i).eq.iopt_radius) then
      ioptfindex(i) = nptr(ii)
    elseif (iopttype(i).eq.iopt_xf) then
      ioptfindex(i) = nptr(ii)
    elseif (iopttype(i).eq.iopt_yf) then
      ioptfindex(i) = nptr(ii)
    elseif (iopttype(i).eq.iopt_zf) then
      ioptfindex(i) = nptr(ii)
    endif
  enddo
!
!  Work out whether full tmat multiplication is needed - use
!  faster method to handle P1 fully relaxed case
!
  if (ncon.eq.0) then
    ltmat = .false.
  else
    ltmat = .true.
  endif
!******************************************
!                                         *
!  Internal derivatives :                 *
!                                         *
!  Symmetry transform second derivatives  *
!******************************************
  if (ltmat) then
    if ((2*nvar).le.maxd2) then
      do i = 1,nvar
        do j = 1,n3f
          tr1 = 0.0_dp
          do k = 1,n3f
            tr1 = tr1 + tmat(k,i)*derv2(k,j)
          enddo
          tmat(j,nvar+i) = tr1
        enddo
      enddo
      do i = 1,nvar
        do j = 1,nvar
          derv2(j,i) = 0.0_dp
          do k = 1,n3f
            derv2(j,i) = derv2(j,i) + tmat(k,j)*tmat(k,nvar+i)
          enddo
        enddo
      enddo
    else
      allocate(tmp2(n3f),stat=status)
      if (status/=0) call outofmemory('sec0f','tmp2')
      do i = 1,nvar
        do j = 1,n3f
          tmp2(j) = 0.0_dp
          do k = 1,n3f
            tmp2(j) = tmp2(j) + tmat(k,i)*derv2(j,k)
          enddo
        enddo
        do j = 1,nvar
          derv2(j,i) = 0.0_dp
          do k = 1,n3f
            derv2(j,i) = derv2(j,i) + tmp2(k)*tmat(k,j)
          enddo
        enddo
      enddo
      deallocate(tmp2,stat=status)
      if (status/=0) call deallocate_error('sec0f','tmp2')
    endif
  elseif (nvar.ne.n3f) then
    do i = 1,nvar
      if (iopttype(i).eq.iopt_xf) then
        ii = 3*ioptfindex(i) - 2
      elseif (iopttype(i).eq.iopt_yf) then
        ii = 3*ioptfindex(i) - 1
      elseif (iopttype(i).eq.iopt_zf) then
        ii = 3*ioptfindex(i)
      elseif (iopttype(i).eq.iopt_radius) then
        ii = 3*nff + ioptfindex(i)
      endif
      do j = 1,nvar
        if (iopttype(j).eq.iopt_xf) then
          jj = 3*ioptfindex(j) - 2
        elseif (iopttype(j).eq.iopt_yf) then
          jj = 3*ioptfindex(j) - 1
        elseif (iopttype(j).eq.iopt_zf) then
          jj = 3*ioptfindex(j)
        elseif (iopttype(j).eq.iopt_radius) then
          jj = 3*nff + ioptfindex(j)
        endif
        derv2(j,i) = derv2(jj,ii)
      enddo
    enddo
  endif
  if (lsdebug.and.ioproc) then
    write(ioout,'(/,''  Symmetrised Second Derivative Matrix  :'',/)')
    do i = 1,nvar
      write(ioout,'(2x,10(f14.4))')(derv2(i,j),j=1,i)
    enddo
    write(ioout,'(/)')
  endif
!
!  Free local memory
!
  deallocate(nptr,stat=status)
  if (status/=0) call deallocate_error('sec0f','nptr')
  deallocate(ioptfindex,stat=status)
  if (status/=0) call deallocate_error('sec0f','ioptfindex')
!
!  Timing
!
  t2 = g_cpu_time()
  thes = thes + t2 - t1
#ifdef TRACE
  call trace_out('sec0f')
#endif
  return
  end
