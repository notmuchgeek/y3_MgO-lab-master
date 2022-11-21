  subroutine transmatd
!
!  Create transformation matrix which maps full set of internal
!  variables onto asymmetric unit variables
!
!  Distributed memory version
!
!  N.B. derv2 is used as a scratch array here to save space
!  as transmat is called prior to calculating the second
!  derivative matrix
!  NB: Not enabled for freezing or BSM at present 
!
!  Matrix is stored as : tmat(3*numat,3*nasym)
!  Strain matrix as    : stmat(6,ncell)
!
!   1/17 Created from transmat
!   3/17 Modifications made to allow for new variable order in iopt
!   4/17 Freezing removed since this is turned off for parallel second derivatives
!   5/17 Symmetry enabled for lsymderv2 false case
!   7/17 Trap added in case number of variables on node is zero
!   2/18 Trace added
!   3/19 iopt replaced by ioptindex and iopttype
!   4/19 Correction to breathing shell cmat
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!   5/20 Constraint handling corrected
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
!  Julian Gale, CIC, Curtin University, May 2020
!
  use control
  use current
  use derivatives
  use iochannels
  use optimisation
  use parallel
  use symmetry
  use times
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  use transform
  implicit none
!
!  Local variables
!
  character(len = 80)                                  :: string
  integer(i4)                                          :: i
  integer(i4)                                          :: ifail
  integer(i4)                                          :: ii
  integer(i4)                                          :: iio
  integer(i4)                                          :: iint
  integer(i4)                                          :: iloc
  integer(i4)                                          :: ind
  integer(i4)                                          :: io
  integer(i4)                                          :: iv
  integer(i4)                                          :: ix
  integer(i4)                                          :: iy
  integer(i4)                                          :: iz
  integer(i4)                                          :: j
  integer(i4)                                          :: jj
  integer(i4)                                          :: jk
  integer(i4)                                          :: jo
  integer(i4)                                          :: jx
  integer(i4)                                          :: jy
  integer(i4)                                          :: jz
  integer(i4)                                          :: jkx
  integer(i4)                                          :: jky
  integer(i4)                                          :: jkz
  integer(i4)                                          :: k
  integer(i4)                                          :: ki
  integer(i4)                                          :: kj
  integer(i4)                                          :: l
  integer(i4)                                          :: lo
  integer(i4)                                          :: mv
  integer(i4)                                          :: n3a
  integer(i4)                                          :: n3f
  integer(i4)                                          :: neqi
  integer(i4)                                          :: neqj
  integer(i4)                                          :: nreli
  integer(i4)                                          :: status
  logical                                              :: lfound
  logical                                              :: ltmat
  real(dp),            dimension(:), allocatable       :: cmat
  real(dp)                                             :: cop(3,3)
  real(dp)                                             :: copi(3,3)
  real(dp)                                             :: g_cpu_time
  real(dp)                                             :: diff
  real(dp)                                             :: diffx
  real(dp)                                             :: diffy
  real(dp)                                             :: diffz
  real(dp)                                             :: rri
  real(dp)                                             :: rtmp(6)
  real(dp)                                             :: t1
  real(dp)                                             :: t2
  real(dp)                                             :: xdiff
  real(dp)                                             :: ydiff
  real(dp)                                             :: zdiff
  real(dp)                                             :: xj1
  real(dp)                                             :: yj1
  real(dp)                                             :: zj1
  real(dp)                                             :: xj
  real(dp)                                             :: yj
  real(dp)                                             :: zj
  real(dp)                                             :: xjj
  real(dp)                                             :: yjj
  real(dp)                                             :: zjj
#ifdef TRACE
  call trace_in('transmatd')
#endif
!
  t1  =  g_cpu_time()
!
!  Work out whether full tmat multiplication is needed - use
!  faster method to handle P1 fully relaxed case
!
  if (ncon.eq.0.and.(nspcg(ncf).le.1.and.ngocfg(ncf).eq.1)) then
    ltmat = .false.
  elseif (ninternal.eq.0) then
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
  allocate(cmat(5*nasym),stat=status)
  if (status/=0) call outofmemory('transmatd','cmat')
!
  n3a = 3*nasym
  n3f = 3*numat
!************************************
!  Loop over cell strain variables  *
!************************************
  if (ncell.gt.0) then
    do i = 1,ncell
      do j = 1,nstrains
        stmat(j,i) = 0.0_dp
      enddo
      iv = ncellmin - 1 + i
      io = ioptindex(iv)
      stmat(io,i) = 1.0_dp
!
!  Apply constraints
!
      if (ncon.gt.0) then
        do k = 1,ncon
          ii = ncvarind(k)
          if (ncvartyp(k).eq.iopt_strain.or.ncvartyp(k).eq.iopt_cell) then
            io = ncfixind(k)
            if (ioptindex(iv).eq.ii) stmat(io,i) = conco(k)
          endif
        enddo
      endif
    enddo
  endif
!
!  If tmat is not needed for internals then exit now
!
  if (.not.ltmat) goto 999
!
!  Check size of tmat
!
  if (lsymderv2) then
    if (n3a.gt.maxn3a) then
      maxn3a = n3a
      call changemaxn3a
    endif
    if (n3f.gt.maxn3f) then
      maxn3f = n3f
      call changemaxn3f
    endif
  else
    if (ninternalonnode.gt.maxn3a) then
      maxn3a = ninternalonnode
      call changemaxn3a
    endif
    if (n3f.gt.maxn3f) then
      maxn3f = n3f
      call changemaxn3f
    endif
  endif
  if (ninternal.gt.0) then
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
        do i = 1,nasym
          ix = ix + 3
          iy = iy + 3
          iz = iz + 3
          neqi = neqv(i)
          ind  = nrela2f(i) - 1
          do j = 1,neqi
            jx = jx + 3
            jy = jy + 3
            jz = jz + 3
            mv = nrotop(ind+j)
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
      endif
    else
!**************************
!  Nasym-numat algorithm  *
!**************************
!
!  Set up symmetry operator tables
!
      call symprod
!
!  Centring handling - generate fractional coordinates
!  for non-centred cell for rotation matching
!
      if (ncbl.gt.1) then
        cop(1,1) = w1(ncbl,1,1)
        cop(2,1) = w1(ncbl,2,1)
        cop(3,1) = w1(ncbl,3,1)
        cop(1,2) = w1(ncbl,1,2)
        cop(2,2) = w1(ncbl,2,2)
        cop(3,2) = w1(ncbl,3,2)
        cop(1,3) = w1(ncbl,1,3)
        cop(2,3) = w1(ncbl,2,3)
        cop(3,3) = w1(ncbl,3,3)
        copi(1,1) = cop(1,1)
        copi(2,1) = cop(2,1)
        copi(3,1) = cop(3,1)
        copi(1,2) = cop(1,2)
        copi(2,2) = cop(2,2)
        copi(3,2) = cop(3,2)
        copi(1,3) = cop(1,3)
        copi(2,3) = cop(2,3)
        copi(3,3) = cop(3,3)
        call matrix_inversion(copi,3_i4,3_i4,rtmp,ifail)
      endif
!****************************
!  Symmetry adapted method  *
!****************************
!
!  Zero tmat
!
      tmat(1:n3f,1:n3a) = 0.0_dp
!
!  Loop over asymmetric unit site
!
      ix = -2
      iy = -1
      iz = 0
      do i = 1,nasym
        ix = ix + 3
        iy = iy + 3
        iz = iz + 3
        neqi = neqv(i)
        rri = 1.0_dp/neqi
        nreli = nrela2f(i)
        ki = nreli - 1
        jx = - 2
        jy = - 1
        jz =   0
        do j = 1,nasym
          neqj = neqv(j)
          kj = nrela2f(j) - 1
          ii = ki
          do k = 1,neqi
            ii = ii + 1
            io = nrotop(ii)
            iio = inverse(io)
            jj = kj
            do l = 1,neqj
              jj = jj + 1
              jo = nrotop(jj)
              xdiff = xfrac(jj) - xfrac(ii)
              ydiff = yfrac(jj) - yfrac(ii)
              zdiff = zfrac(jj) - zfrac(ii)
              if (ncbl.gt.1) then
                xj = xdiff*copi(1,1) + ydiff*copi(1,2) + zdiff*copi(1,3)
                yj = xdiff*copi(2,1) + ydiff*copi(2,2) + zdiff*copi(2,3)
                zj = xdiff*copi(3,1) + ydiff*copi(3,2) + zdiff*copi(3,3)
              else
                xj = xdiff
                yj = ydiff
                zj = zdiff
              endif
!
!  Transformation element
!  Inverse of io * jo
!
              lo = iptab(iio,jo)
!
!  Find equivalent derv2 matrix element in asymmetric unit
!
              lfound = .false.
              jk = 0
              do while (jk.lt.neqj.and..not.lfound)
                jk = jk + 1
                lfound = (nrotop(jk+kj).eq.lo)
              enddo
              if (.not.lfound) then
                jk = 0
                xjj = xj*rop(1,1,iio) + yj*rop(1,2,iio) + zj*rop(1,3,iio)
                yjj = xj*rop(2,1,iio) + yj*rop(2,2,iio) + zj*rop(2,3,iio)
                zjj = xj*rop(3,1,iio) + yj*rop(3,2,iio) + zj*rop(3,3,iio)
                if (ncbl.gt.1) then
                  xj1 = xjj
                  yj1 = yjj
                  zj1 = zjj
                  xjj = xj1*cop(1,1) + yj1*cop(1,2) + zj1*cop(1,3)
                  yjj = xj1*cop(2,1) + yj1*cop(2,2) + zj1*cop(2,3)
                  zjj = xj1*cop(3,1) + yj1*cop(3,2) + zj1*cop(3,3)
                endif
!
!  Place vector relative to atom i in asymmetric unit
!
                xjj = xjj + xfrac(nreli)
                yjj = yjj + yfrac(nreli)
                zjj = zjj + zfrac(nreli)
!
!  Place within range 0-1
!
                xjj = mod(xjj+3.0_dp,1.0_dp)
                yjj = mod(yjj+3.0_dp,1.0_dp)
                zjj = mod(zjj+3.0_dp,1.0_dp)
                do while (jk.lt.neqj.and..not.lfound)
                  jk = jk+1
                  diffx = abs(xjj-xfrac(jk+kj))
                  diffy = abs(yjj-yfrac(jk+kj))
                  diffz = abs(zjj-zfrac(jk+kj))
                  diffx = mod(diffx,1.0_dp)
                  diffy = mod(diffy,1.0_dp)
                  diffz = mod(diffz,1.0_dp)
                  if (diffx.gt.0.5_dp) diffx = 1.0_dp - diffx
                  if (diffy.gt.0.5_dp) diffy = 1.0_dp - diffy
                  if (diffz.gt.0.5_dp) diffz = 1.0_dp - diffz
                  diff = diffx + diffy + diffz
                  lfound = (diff.lt.1.0d-4)
                enddo
              endif
              if (.not.lfound) then
                if (ioproc) then
                  string  =  'symmetry atom not found in transmatd'
                  call outerror(string,0_i4)
                endif
                call stopnow('transmatd')
              endif
              jkx = jx + 3*jk
              jky = jy + 3*jk
              jkz = jz + 3*jk
!
!  Add term to tmat
!
              tmat(jkx,ix) = tmat(jkx,ix) + rop(1,1,lo)*rri
              tmat(jky,ix) = tmat(jky,ix) + rop(2,1,lo)*rri
              tmat(jkz,ix) = tmat(jkz,ix) + rop(3,1,lo)*rri
              tmat(jkx,iy) = tmat(jkx,iy) + rop(1,2,lo)*rri
              tmat(jky,iy) = tmat(jky,iy) + rop(2,2,lo)*rri
              tmat(jkz,iy) = tmat(jkz,iy) + rop(3,2,lo)*rri
              tmat(jkx,iz) = tmat(jkx,iz) + rop(1,3,lo)*rri
              tmat(jky,iz) = tmat(jky,iz) + rop(2,3,lo)*rri
              tmat(jkz,iz) = tmat(jkz,iz) + rop(3,3,lo)*rri
            enddo
          enddo
          jx = jx + 3*neqj
          jy = jy + 3*neqj
          jz = jz + 3*neqj
        enddo
      enddo
    endif
!*********************************
!  Loop over internal variables  *
!*********************************
    if (ninternalatmonnode.gt.0) then
      ii = 0
      do iint = 1,ninternalatm
        i = iint + ninternalmin - 1
        if (nvar2local(i).gt.0) then
          ii = ii + 1
          cmat(1:5*nasym) = 0.0_dp
!
!  Variable mapping matrix
!
          io = ioptindex(i)
          if (iopttype(i).eq.iopt_xf) then
            ind = 3*nasymnomolptr(io) - 2
            cmat(ind) = 1.0_dp
          elseif (iopttype(i).eq.iopt_yf) then
            ind = 3*nasymnomolptr(io) - 1
            cmat(ind) = 1.0_dp
          elseif (iopttype(i).eq.iopt_zf) then
            ind = 3*nasymnomolptr(io)
            cmat(ind) = 1.0_dp
          elseif (iopttype(i).eq.iopt_radius) then
            cmat(n3a+io) = 1.0_dp
          endif
!
!  Apply constraints
!
          if (ncon.gt.0) then
            do k = 1,ncon
              if (ncvartyp(k).eq.iopttype(i).and.ncvarind(k).eq.ioptindex(i)) then
                io = ncfixind(k)
                if (ncfixtyp(k).eq.iopt_xf) then
                  ind = 3*nasymnomolptr(io) - 2
                  cmat(ind) = conco(k)
                elseif (ncfixtyp(k).eq.iopt_yf) then
                  ind = 3*nasymnomolptr(io) - 1
                  cmat(ind) = conco(k)
                elseif (ncfixtyp(k).eq.iopt_zf) then
                  ind = 3*nasymnomolptr(io)
                  cmat(ind) = conco(k)
                elseif (ncfixtyp(k).eq.iopt_radius) then
                  cmat(3*nasym+io) = conco(k)
                endif
              endif
            enddo
          endif
!
!  Reduce cmat from 3*nasym to nasym
!
          if (lsymderv2) then
            do j = 1,n3a
              tmat(j,n3a+ii) = cmat(j)
            enddo
          else
!
!  Combine transformation matrices: In case without symmetry this can be
!  simplified since the matrix in derv2 is diagonal.
!
            if (nspcg(ncf).gt.1.or.ngocfg(ncf).gt.1) then
              do j = 1,n3f
                tmat(j,ii) = 0.0_dp
                do k = 1,n3a
                  tmat(j,ii) = tmat(j,ii) + derv2(j,k)*cmat(k)
                enddo
              enddo
            else
              do j = 1,n3f
                tmat(j,ii) = cmat(j)
              enddo
            endif
          endif
        endif
      enddo
    endif
  endif
!***********************
!  Debugging printing  *
!***********************
  if (index(keyword,'tmat').ne.0) then
    call mpbarrier
    if (lsymderv2) then
!
!  Nasym-numat method
!
      if (ioproc) then
        write(ioout,'(/,''  Transformation matrix : '',/)')
        do i = 1,n3a
          write(ioout,'(2x,12f6.3)')(tmat(j,i),j=1,n3f)
        enddo
        write(ioout,'(/,''  Reduction matrix : '',/)')
      endif
      call mpbarrier
      do i = 1,ninternal
        iloc = nvar2local(i)
        if (iloc.gt.0) then
          write(ioout,'(2x,12f6.3)')(tmat(j,n3a+i),j=1,n3a)
        endif
        call mpbarrier
      enddo
      if (ioproc) then
        write(ioout,'(/)')
      endif
    else
!
!  Numat-numat method
!
      if (ioproc) then
        write(ioout,'(/,''  Transformation matrix : '',/)')
      endif
      call mpbarrier
      do i = ninternalmin,ninternalmax
        ii = nvar2local(i)
        if (ii.gt.0) then
          ii = ii - ninternalminonnode + 1
          write(ioout,'(2x,12f6.3)')(tmat(j,ii),j=1,n3f)
        endif
        call mpbarrier
      enddo
      call mpbarrier
      if (ioproc) then
        write(ioout,'(/)')
      endif
      if (ncell.gt.0.and.ioproc) then
        write(ioout,'(/,''  Strain transformation matrix : '',/)')
        do i = 1,ncell
          write(ioout,'(2x,6f6.3)')(stmat(j,i),j=1,nstrains)
        enddo
        write(ioout,'(/)')
      endif
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
    if (status/=0) call deallocate_error('transmatd','cmat')
  endif
!
!  Timing
!
  t2 = g_cpu_time()
  ttmat = ttmat + t2 - t1
#ifdef TRACE
  call trace_out('transmatd')
#endif
!
  return
  end
