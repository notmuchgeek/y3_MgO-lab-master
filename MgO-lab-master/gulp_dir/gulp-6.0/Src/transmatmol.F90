  subroutine transmatmolT
!
!  Create transformation matrix which maps full set of internal
!  variables onto asymmetric unit variables for rigid molecules
!  - translation variables
!
!  N.B. derv2 is used as a scratch array here to save space
!  as transmatmol is called prior to calculating the second
!  derivative matrix
!
!  N.B. rigid molecules are incompatible with freezing and so this
!  is not allowed for here
!
!  Matrix are stored as : tmatT(3*nmol,3*nmolasym) - translation of molecules
!
!   2/20 Created from transmat
!   3/20 Made translation specific
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
!  Copyright Curtin University 2020
!
!  Julian Gale, CIC, Curtin University, March 2020
!
  use control
  use current
  use derivatives
  use iochannels
  use molecule
  use optimisation
  use parallel
  use symmetry
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  use transform
  implicit none
!
!  Local variables
!
  integer(i4)                                          :: i
  integer(i4)                                          :: iint
  integer(i4)                                          :: ind
  integer(i4)                                          :: io
  integer(i4)                                          :: ix
  integer(i4)                                          :: iy
  integer(i4)                                          :: iz
  integer(i4)                                          :: j
  integer(i4)                                          :: jx
  integer(i4)                                          :: jy
  integer(i4)                                          :: jz
  integer(i4)                                          :: k
  integer(i4)                                          :: mv
  integer(i4)                                          :: n3a
  integer(i4)                                          :: n3f
  integer(i4)                                          :: neqi
  integer(i4)                                          :: nfa
  integer(i4)                                          :: nff
  integer(i4)                                          :: ni
  integer(i4)                                          :: status
  logical                                              :: ltmat
  real(dp),            dimension(:), allocatable       :: cmat
  real(dp)                                             :: g_cpu_time
  real(dp)                                             :: t1
  real(dp)                                             :: t2
#ifdef TRACE
  call trace_in('transmatmolT')
#endif
!
  t1  =  g_cpu_time()
!
!  Work out whether full tmat multiplication is needed - use
!  faster method to handle P1 fully relaxed case
!
  if (ncon.eq.0.and.(nspcg(ncf).le.1.and.ngocfg(ncf).eq.1)) then
    ltmat = .false.
  elseif (ninternalmol.eq.0) then
    ltmat = .false.
  else
    ltmat = .true.
  endif
!
!  If tmat is not needed and this is a cluster return now
!
  if (.not.ltmat.and.ndim.eq.0) goto 999
!
!  Allocate local memory
!
  allocate(cmat(3*nmolasym),stat=status)
  if (status/=0) call outofmemory('transmatmolT','cmat')
!
  nfa = nmolasym
  nff = nmol
!
  n3a = 3*nfa
  n3f = 3*nff
!
!  If tmat is not needed for internals then exit now
!
  if (.not.ltmat) goto 999
!
!  Check size of tmat
!
  if (n3a+ninternalmol.gt.maxn3ma) then
    maxn3ma = n3a + ninternal
    call changemaxn3ma
  endif
  if (n3f.gt.maxn3mf) then
    maxn3mf = n3f
    call changemaxn3mf
  endif
  if (ninternalmolT.gt.0) then
    if (.not.lsymderv2) then
      if (n3a.gt.maxd2u) then
        maxd2u = n3a
        call changemaxd2
      endif
      if (n3f.gt.maxd2) then
        maxd2 = n3f
        call changemaxd2
      endif
!**************************
!  Numat-numat algorithm  *
!**************************
!
!  Zero array
!
      derv2(1:n3f,1:n3a) = 0.0_dp
!************************************
!  Collect transformation elements  *
!************************************
!
!  Coordinates
!
      if (nspcg(ncf).gt.1.or.ngocfg(ncf).gt.1) then
        ix = - 2
        iy = - 1
        iz =   0
        jx = - 2
        jy = - 1
        jz =   0
        do i = 1,nfa
          ix = ix + 3
          iy = iy + 3
          iz = iz + 3
          neqi = nmoleqv(i)
          ind  = nmola2f(i) - 1
          do ni = 1,neqi
            j = nmolasymeqvptr(ni,i)
            mv = nrotop(nmollist(nmolptr(j)+1))
            jx = jx + 3
            jy = jy + 3
            jz = jz + 3
            derv2(jx,ix) = rop(1,1,mv)
            derv2(jy,ix) = rop(2,1,mv)
            derv2(jz,ix) = rop(3,1,mv)
            derv2(jx,iy) = rop(1,2,mv)
            derv2(jy,iy) = rop(2,2,mv)
            derv2(jz,iy) = rop(3,2,mv)
            derv2(jx,iz) = rop(1,3,mv)
            derv2(jy,iz) = rop(2,3,mv)
            derv2(jz,iz) = rop(3,3,mv)
          enddo
        enddo
      else
        do i = 1,n3a
          derv2(i,i) = 1.0_dp
        enddo
      endif
    else
!**************************
!  Nasym-numat algorithm  *
!**************************
!
!  Trap this algorithm since it's not yet implemented
!
      call outerror('symmetry adapted second derivative algorithm for rigid molecules not implemented',0_i4)
      call stopnow('transmatmolT')
    endif
!*********************************************
!  Loop over molecule translation variables  *
!*********************************************
    do iint = 1,ninternalmolT
      i = ninternalmin + ninternalatm + iint - 1
      do j = 1,3*nmolasym
        cmat(j) = 0.0_dp
      enddo
!
!  Variable mapping matrix
!
      io = ioptindex(i)
      if (iopttype(i).eq.iopt_xcom) then
        cmat(3*io-2) = 1.0_dp
      elseif (iopttype(i).eq.iopt_ycom) then
        cmat(3*io-1) = 1.0_dp
      elseif (iopttype(i).eq.iopt_zcom) then
        cmat(3*io) = 1.0_dp
      endif
!
!  Apply constraints
!
      if (ncon.gt.0) then
        do k = 1,ncon
          if (ncvartyp(k).eq.iopttype(i).and.ncvarind(k).eq.ioptindex(i)) then
            io = ncfixind(k)
            if (ncfixtyp(k).eq.iopt_xcom) then
              cmat(3*io-2) = conco(k)
            elseif (ncfixtyp(k).eq.iopt_ycom) then
              cmat(3*io-1) = conco(k)
            elseif (ncfixtyp(k).eq.iopt_zcom) then
              cmat(3*io) = conco(k)
            endif
          endif
        enddo
      endif
!
      if (.not.lsymderv2) then
!
!  Combine transformation matrices
!
        do j = 1,n3f
          tmatT(j,iint) = 0.0_dp
          do k = 1,n3a
            tmatT(j,iint) = tmatT(j,iint) + derv2(j,k)*cmat(k)
          enddo
        enddo
      else
        do j = 1,n3a
          tmatT(j,iint+n3a) = cmat(j)
        enddo
      endif
    enddo
  endif
!***********************
!  Debugging printing  *
!***********************
  if (index(keyword,'tmat').ne.0.and.ioproc) then
    if (lsymderv2) then
!
!  Nasym-numat method
!
      write(ioout,'(/,''  Transformation matrix for rigid molecule translations : '',/)')
      do i = 1,n3a
        write(ioout,'(2x,12f6.3)')(tmatT(j,i),j=1,n3f)
      enddo
      write(ioout,'(/,''  Reduction matrix for rigid molecule translations : '',/)')
      do i = 1,n3a
        write(ioout,'(2x,12f6.3)')(tmatT(i,n3a+j),j=1,ninternalmolT)
      enddo
      write(ioout,'(/)')
    else
!
!  Numat-numat method
!
      write(ioout,'(/,''  Transformation matrix for rigid molecule translations : '',/)')
      do i = 1,ninternalmolT
        write(ioout,'(2x,12f6.3)')(tmatT(j,i),j=1,n3f)
      enddo
      write(ioout,'(/)')
    endif
  endif
!
!  Exit point
!
999 continue
!
!  Free local memory
!
  if (allocated(cmat)) then
    deallocate(cmat,stat=status)
    if (status/=0) call deallocate_error('transmatmolT','cmat')
  end if
!
!  Timing
!
  t2 = g_cpu_time()
  ttmat = ttmat + t2 - t1
#ifdef TRACE
  call trace_out('transmatmolT')
#endif
!
  return
  end
!
  subroutine transmatmolQ
!
!  Create transformation matrix which maps full set of internal
!  variables onto asymmetric unit variables for rigid molecules
!  - quaternion variables
!
!  N.B. derv2 is used as a scratch array here to save space
!  as transmatmol is called prior to calculating the second
!  derivative matrix
!
!  N.B. rigid molecules are incompatible with freezing and so this
!  is not allowed for here
!
!  Matrix are stored as : tmatQ(3*nmol,3*nmolasym) - translation of molecules
!
!   3/20 Created from transmatmolT
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
!  Copyright Curtin University 2020
!
!  Julian Gale, CIC, Curtin University, March 2020
!
  use control
  use current
  use derivatives
  use iochannels
  use molecule
  use optimisation
  use parallel
  use symmetry
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  use transform
  implicit none
!
!  Local variables
!
  integer(i4)                                          :: i
  integer(i4)                                          :: iint
  integer(i4)                                          :: ind
  integer(i4)                                          :: io
  integer(i4)                                          :: ix
  integer(i4)                                          :: iy
  integer(i4)                                          :: iz
  integer(i4)                                          :: j
  integer(i4)                                          :: jx
  integer(i4)                                          :: jy
  integer(i4)                                          :: jz
  integer(i4)                                          :: k
  integer(i4)                                          :: n3a
  integer(i4)                                          :: n3f
  integer(i4)                                          :: neqi
  integer(i4)                                          :: nfa
  integer(i4)                                          :: nff
  integer(i4)                                          :: ni
  integer(i4)                                          :: status
  logical                                              :: ltmat
  real(dp),            dimension(:), allocatable       :: cmat
  real(dp)                                             :: g_cpu_time
  real(dp)                                             :: t1
  real(dp)                                             :: t2
#ifdef TRACE
  call trace_in('transmatmolQ')
#endif
!
  t1  =  g_cpu_time()
!
!  Work out whether full tmat multiplication is needed - use
!  faster method to handle P1 fully relaxed case
!
  if (ncon.eq.0.and.(nspcg(ncf).le.1.and.ngocfg(ncf).eq.1)) then
    ltmat = .false.
  elseif (ninternalmol.eq.0) then
    ltmat = .false.
  else
    ltmat = .true.
  endif
!
!  If tmat is not needed and this is a cluster return now
!
  if (.not.ltmat.and.ndim.eq.0) goto 999
!
!  Allocate local memory
!
  allocate(cmat(3*nmolasym),stat=status)
  if (status/=0) call outofmemory('transmatmolQ','cmat')
!
  nfa = nmolasym
  nff = nmol
!
  n3a = 3*nfa
  n3f = 3*nff
!
!  If tmat is not needed for internals then exit now
!
  if (.not.ltmat) goto 999
!
!  Check size of tmat
!
  if (n3a+ninternalmol.gt.maxn3ma) then
    maxn3ma = n3a + ninternal
    call changemaxn3ma
  endif
  if (n3f.gt.maxn3mf) then
    maxn3mf = n3f
    call changemaxn3mf
  endif
  if (ninternalmolQ.gt.0) then
    if (.not.lsymderv2) then
      if (n3a.gt.maxd2u) then
        maxd2u = n3a
        call changemaxd2
      endif
      if (n3f.gt.maxd2) then
        maxd2 = n3f
        call changemaxd2
      endif
!**************************
!  Numat-numat algorithm  *
!**************************
!
!  Zero array
!
      derv2(1:n3f,1:n3a) = 0.0_dp
!************************************
!  Collect transformation elements  *
!************************************
      if (nspcg(ncf).gt.1.or.ngocfg(ncf).gt.1) then
        ix = - 2
        iy = - 1
        iz =   0
        jx = - 2
        jy = - 1
        jz =   0
        do i = 1,nmolasym
          ix = ix + 3
          iy = iy + 3
          iz = iz + 3
          neqi = nmoleqv(i)
          ind  = nmola2f(i) - 1
          do ni = 1,neqi
            jx = jx + 3
            jy = jy + 3
            jz = jz + 3
!
!  Symmetry-adapted quaternions are already referenced to the values in the asymmetric unit
!  => no symmetry operator needed
!
            derv2(jx,ix) = 1.0_dp
            derv2(jy,ix) = 0.0_dp
            derv2(jz,ix) = 0.0_dp
            derv2(jx,iy) = 0.0_dp
            derv2(jy,iy) = 1.0_dp
            derv2(jz,iy) = 0.0_dp
            derv2(jx,iz) = 0.0_dp
            derv2(jy,iz) = 0.0_dp
            derv2(jz,iz) = 1.0_dp
          enddo
        enddo
      else
        do i = 1,n3a
          derv2(i,i) = 1.0_dp
        enddo
      endif
    else
!**************************
!  Nasym-numat algorithm  *
!**************************
!
!  Trap this algorithm since it's not yet implemented
!
      call outerror('symmetry adapted second derivative algorithm for rigid molecules not implemented',0_i4)
      call stopnow('transmatmolQ')
    endif
!*********************************************
!  Loop over molecule quaternion variables  *
!*********************************************
    do iint = 1,ninternalmolQ
      i = ninternalmin + ninternalatm + ninternalmolT + iint - 1
      do j = 1,3*nmolasym
        cmat(j) = 0.0_dp
      enddo
!
!  Variable mapping matrix
!
      io = ioptindex(i)
      if (iopttype(i).eq.iopt_xqtn) then
        cmat(3*io-2) = 1.0_dp
      elseif (iopttype(i).eq.iopt_yqtn) then
        cmat(3*io-1) = 1.0_dp
      elseif (iopttype(i).eq.iopt_zqtn) then
        cmat(3*io) = 1.0_dp
      endif
!
!  Apply constraints
!
      if (ncon.gt.0) then
        do k = 1,ncon
          if (ncvartyp(k).eq.iopttype(i).and.ncvarind(k).eq.ioptindex(i)) then
            io = ncfixind(k)
            if (ncfixtyp(k).eq.iopt_xqtn) then
              cmat(3*io-2) = conco(k)
            elseif (ncfixtyp(k).eq.iopt_yqtn) then
              cmat(3*io-1) = conco(k)
            elseif (ncfixtyp(k).eq.iopt_zqtn) then
              cmat(3*io) = conco(k)
            endif
          endif
        enddo
      endif
!
      if (.not.lsymderv2) then
!
!  Combine transformation matrices
!
        do j = 1,n3f
          tmatQ(j,iint) = 0.0_dp
          do k = 1,n3a
            tmatQ(j,iint) = tmatQ(j,iint) + derv2(j,k)*cmat(k)
          enddo
        enddo
      else
        do j = 1,n3a
          tmatQ(j,iint+n3a) = cmat(j)
        enddo
      endif
    enddo
  endif
!***********************
!  Debugging printing  *
!***********************
  if (index(keyword,'tmat').ne.0.and.ioproc) then
    if (lsymderv2) then
!
!  Nasym-numat method
!
      write(ioout,'(/,''  Transformation matrix for rigid molecule quaternions : '',/)')
      do i = 1,n3a
        write(ioout,'(2x,12f6.3)')(tmatQ(j,i),j=1,n3f)
      enddo
      write(ioout,'(/,''  Reduction matrix for rigid molecule quaternions : '',/)')
      do i = 1,n3a
        write(ioout,'(2x,12f6.3)')(tmatQ(i,n3a+j),j=1,ninternalmolQ)
      enddo
      write(ioout,'(/)')
    else
!
!  Numat-numat method
!
      write(ioout,'(/,''  Transformation matrix for rigid molecule quaternions : '',/)')
      do i = 1,ninternalmolQ
        write(ioout,'(2x,12f6.3)')(tmatQ(j,i),j=1,n3f)
      enddo
      write(ioout,'(/)')
    endif
  endif
!
!  Exit point
!
999 continue
!
!  Free local memory
!
  if (allocated(cmat)) then
    deallocate(cmat,stat=status)
    if (status/=0) call deallocate_error('transmatmolQ','cmat')
  end if
!
!  Timing
!
  t2 = g_cpu_time()
  ttmat = ttmat + t2 - t1
#ifdef TRACE
  call trace_out('transmatmolQ')
#endif
!
  return
  end
