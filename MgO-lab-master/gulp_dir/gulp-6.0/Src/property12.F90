  subroutine property12(lprint)
!
!  Calculate the properties of 1-D and 2-D systems.
!
!  10/13 Created from property0
!   3/14 n3f defined
!  12/14 nbs, nbss and nbsptr moved to shells modul
!   2/17 Blocksize added to call to matrix_inversion_library
!   1/17 Trace added
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
  use control,        only : lraman
  use current
  use derivatives,    only : derv2, dervi, maxd2
  use element
  use general,        only : nwarn
  use iochannels
  use parallel
  use partial
  use shells,         only : nshell, nshptr, nbs, nbss, nbsptr
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  logical,                            intent(in) :: lprint
!
!  Local variables
!
  integer(i4)                                    :: i
  integer(i4)                                    :: ifail
  integer(i4)                                    :: indk
  integer(i4)                                    :: indl
  integer(i4)                                    :: iptr
  integer(i4)                                    :: ix
  integer(i4)                                    :: iy
  integer(i4)                                    :: iz
  integer(i4)                                    :: ixs
  integer(i4)                                    :: iys
  integer(i4)                                    :: izs
  integer(i4)                                    :: j
  integer(i4)                                    :: jptr
  integer(i4)                                    :: jx
  integer(i4)                                    :: jy
  integer(i4)                                    :: jz
  integer(i4)                                    :: jxs
  integer(i4)                                    :: jys
  integer(i4)                                    :: jzs
  integer(i4)                                    :: k
  integer(i4)                                    :: msvar
  integer(i4)                                    :: n
  integer(i4)                                    :: n3f
  integer(i4)                                    :: nbsc
  integer(i4)                                    :: ni
  integer(i4)                                    :: nj
  integer(i4)                                    :: status
  real(dp)                                       :: g_cpu_time
  real(dp)                                       :: qak
  real(dp),    dimension(:,:), allocatable       :: qD
  real(dp),    dimension(:),   allocatable       :: qshell
  real(dp)                                       :: t1p
  real(dp)                                       :: t2p
#ifdef TRACE
  call trace_in('property12')
#endif
!
  t1p = g_cpu_time()
!
  n3f = 3*numat
!**************************
!  Raman susceptibilites  *
!**************************
  if (lraman.and.nshell.gt.0) then
!
!  Allocate local memory
!
    allocate(qshell(numat),stat=status)
    if (status/=0) call outofmemory('property12','qshell')
    allocate(qD(3*nsfoc,3),stat=status)
    if (status/=0) call outofmemory('property12','qD')
!
!  Collect shell second derivative terms
!
    msvar = 3*nshell
    do i = 1,nshell
      ni = nshptr(i)
      qshell(i) = qf(ni)*occuf(ni)
      ix = 3*(ni-1) + 1
      iy = ix + 1
      iz = ix + 2
      ixs = 3*(i-1) + 1
      iys = ixs + 1
      izs = ixs + 2
      do j = 1,nshell
        nj = nshptr(j)
        jx = 3*(nj-1) + 1
        jy = jx + 1
        jz = jx + 2
        jxs = 3*(j-1) + 1
        jys = jxs + 1
        jzs = jxs + 2
        dervi(jxs,ixs) = derv2(jx,ix)
        dervi(jys,ixs) = derv2(jy,ix)
        dervi(jzs,ixs) = derv2(jz,ix)
        dervi(jxs,iys) = derv2(jx,iy)
        dervi(jys,iys) = derv2(jy,iy)
        dervi(jzs,iys) = derv2(jz,iy)
        dervi(jxs,izs) = derv2(jx,iz)
        dervi(jys,izs) = derv2(jy,iz)
        dervi(jzs,izs) = derv2(jz,iz)
      enddo
    enddo
    if (nbss.gt.0) then
!
!  Collect radial second derivatives
!
      if (lpocc) then
        do i = 1,nshell
          ni = nshptr(i)
          do j = 1,nshell
            nj = nshptr(j)
            dervi(msvar+j,msvar+i) = derv2(n3f+nj,n3f+ni)
          enddo
          do j = 1,nshell
            nj = nshptr(j)
            jx = 3*(nj-1) + 1
            jy = jx + 1
            jz = jx + 2
            jxs = 3*(j-1) + 1
            jys = jxs + 1
            jzs = jxs + 2
            dervi(jxs,msvar+i) = derv2(jx,n3f+ni)
            dervi(jys,msvar+i) = derv2(jy,n3f+ni)
            dervi(jzs,msvar+i) = derv2(jz,n3f+ni)
            dervi(msvar+i,jxs) = derv2(n3f+ni,jx)
            dervi(msvar+i,jys) = derv2(n3f+ni,jy)
            dervi(msvar+i,jzs) = derv2(n3f+ni,jz)
          enddo
        enddo
      else
        nbsc = nbs - nbss
        do i = 1,nbss
          iptr = n3f + nbsptr(nbsc+i)
          do j = 1,nbss
            jptr = n3f + nbsptr(nbsc+j)
            dervi(msvar+j,msvar+i) = derv2(jptr,iptr)
          enddo
          do j = 1,nshell
            nj = nshptr(j)
            jx = 3*(nj-1) + 1
            jy = jx + 1
            jz = jx + 2
            jxs = 3*(j-1) + 1
            jys = jxs + 1
            jzs = jxs + 2
            dervi(jxs,msvar+i) = derv2(jx,iptr)
            dervi(jys,msvar+i) = derv2(jy,iptr)
            dervi(jzs,msvar+i) = derv2(jz,iptr)
          enddo
        enddo
      endif
    endif
!
!  Compress derivatives for partial occupancies
!
    if (lpocc) then
      call compressd1(qshell,0_i4,nsfoc,nshell,iocshptr)
      call compressd2(dervi,maxd2,0_i4,nsfoc,nbsfoc,nshell,iocshptr,ibocshptr)
      n = 3*nsfoc + nbsfoc
    else
      n = msvar + nbss
    endif
!
!  Invert second derivative matrix
!
    ifail = 1
    call matrix_inversion_library(n,1_i4,maxd2,3_i4*nblocksize,dervi,0_i4,ifail)
!
    if (ifail.ne.0) then
      nwarn = nwarn + 1
      call outwarning('Properties cannot be calculated - matrix is singular',0_i4)
      goto 999
    endif
!
!  Multiply inverse second derivatives by charge vectors:
!
!  First multiply by one charge vector and store in dervi for use in Raman calculation, if needed
!
    do k = 1,nsfoc
      qak = qshell(k)
      indk = 3*(k-1)
      do indl = 1,3*nsfoc
        dervi(indl,indk+1) = dervi(indl,indk+1)*qak
      enddo
      do indl = 1,3*nsfoc
        dervi(indl,indk+2) = dervi(indl,indk+2)*qak
      enddo
      do indl = 1,3*nsfoc
        dervi(indl,indk+3) = dervi(indl,indk+3)*qak
      enddo
    enddo
!
!  Compack dervi down into qD
!
    qD(1:3*nsfoc,1:3) = 0.0_dp
    do k = 1,nsfoc
      indk = 3*(k-1)
      do i = 1,3
        do indl = 1,3*nsfoc
          qD(indl,i) = qD(indl,i) + dervi(indl,indk+i)
        enddo
      enddo
    enddo
!
!  Call routine to compute the third derivatives and generate the susceptibility tensors
!
    call raman3(qD,3*nsfoc)
!
!  Deallocate memory
!
    deallocate(qD,stat=status)
    if (status/=0) call deallocate_error('property12','qD')
    deallocate(qshell,stat=status)
    if (status/=0) call deallocate_error('property12','qshell')
  endif
999 continue
!***************
!  Exit tasks  *
!***************
!
!  Timings
!
  t2p = g_cpu_time()
  tprop = t2p - t1p + tprop
#ifdef TRACE
  call trace_out('property12')
#endif
!
  return
  end
