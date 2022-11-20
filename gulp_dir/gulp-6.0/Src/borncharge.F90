  subroutine borncharge(lprint)
!
!  Calculates the Born effective charges for a system
!
!   4/02 Created from property
!  11/02 Basic code for calculation of Born charges for EEM/QEq added
!        but needs phase issue to be addressed
!   1/03 Full charge tensor now output
!   2/03 Matrix inversion accelerated through packed storage
!   3/03 Partial occupancy correction added
!   7/03 Handling of symmetry for output added
!   6/05 Style updated for continuations
!   5/07 Shell pointers moved to module
!   6/07 Modified to handle 2-D case where charges are not to be computed
!        for region 2.
!   2/09 Modified to accommodate new version of FoX and gulp_cml
!   9/16 Modified to allow for parallel second derivatives
!   9/16 Shell pointers modified
!   5/17 nshell and nshellonnode added to subroutine call
!   7/17 Parallelisation added for dqdxyz
!   7/17 Multiple region modifications added
!   7/17 Fix for core only case
!   2/18 Trace added
!   4/18 Output format changed to allow for more atoms
!   5/18 Output format corrected to add extra space on lines without x,y,z
!   7/18 Global sum now only performed if lbornloc is true
!   8/19 Reference to nshell corrected to nshellr1
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!   6/20 Handling of shells in rigid molecules added
!   6/20 Inclusion of rigid molecules shells added
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
!  Julian Gale, CIC, Curtin University, June 2020
!
  use configurations, only : nregionno
  use control,        only : leem, lrigid
  use current
  use derivatives
  use element
  use gulp_cml,       only : lcml
  use gulp_cml_props, only : gulp_cml_print_born_charges
  use iochannels
  use molecule
  use parallel
  use partial
  use properties
  use shells
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  logical,     intent(in)                      :: lprint
!
!  Local variables
!
  character(len=5)                             :: lab
  integer(i4)                                  :: i
  integer(i4)                                  :: ifail
  integer(i4)                                  :: ii
  integer(i4)                                  :: inat
  integer(i4)                                  :: is
  integer(i4)                                  :: itype
  integer(i4)                                  :: ix
  integer(i4)                                  :: iy
  integer(i4)                                  :: iz
  integer(i4)                                  :: ixs
  integer(i4)                                  :: iys
  integer(i4)                                  :: izs
  integer(i4)                                  :: j
  integer(i4)                                  :: jfoc
  integer(i4)                                  :: jj
  integer(i4)                                  :: jx
  integer(i4)                                  :: jy
  integer(i4)                                  :: jz
  integer(i4)                                  :: jxs
  integer(i4)                                  :: jys
  integer(i4)                                  :: jzs
  integer(i4)                                  :: k
  integer(i4)                                  :: kx
  integer(i4)                                  :: ky
  integer(i4)                                  :: kz
  integer(i4)                                  :: msvar
  integer(i4)                                  :: n
  integer(i4)                                  :: ni
  integer(i4)                                  :: nj
  integer(i4)                                  :: n3f
  integer(i4)                                  :: n3fu
  integer(i4)                                  :: nshellr1onnode
  integer(i4), dimension(:),     allocatable   :: nshellr1ptr
  integer(i4)                                  :: numatr1
  integer(i4)                                  :: status
  logical                                      :: lbornloc
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: t1p
  real(dp)                                     :: t2p
  real(dp),    dimension(:),     allocatable   :: qshell
  real(dp),    dimension(:,:,:), allocatable   :: bqtmp
  real(dp)                                     :: xi
  real(dp)                                     :: yi
  real(dp)                                     :: zi
#ifdef TRACE
  call trace_in('borncharge')
#endif
!
  t1p = g_cpu_time()
  lbornloc = (.not.leem)
!
!  Set numatr1 to be sum of relevant cores and shells
!
  numatr1 = ncorer1 + nshellr1
!
!  Check second derivative memory
!
  n3f = 3*numat
  if (nprocs.gt.1) then
    n3fu = 3_i4*natomsonnode
  else
    n3fu = n3f
  endif
  if (n3f.gt.maxd2.or.n3fu.gt.maxd2u) then
    maxd2u = max(n3fu,maxd2u)
    maxd2  = max(n3f,maxd2)
    call changemaxd2
  endif
!
!  Allocate local memory
!
  allocate(qshell(numatr1),stat=status)
  if (status/=0) call outofmemory('borncharge','qshell')
  if (nshellr1.gt.0) then
    allocate(nshellr1ptr(nshellr1),stat=status)
    if (status/=0) call outofmemory('borncharge','nshellr1ptr')
  endif
!***************************
!  Born effective charges  *
!***************************
  if (lbornloc) then
!
!  Initialise Born charges using core contribution
!
    do i = 1,numat
      if (nregionno(nsft+nrelf2a(i)).eq.1.and.nat(i).le.maxele) then
        bornq(1,1,i) = qf(i)*occuf(i)/dble(nprocs)
        bornq(2,1,i) = 0.0_dp
        bornq(3,1,i) = 0.0_dp
        bornq(1,2,i) = 0.0_dp
        bornq(2,2,i) = qf(i)*occuf(i)/dble(nprocs)
        bornq(3,2,i) = 0.0_dp
        bornq(1,3,i) = 0.0_dp
        bornq(2,3,i) = 0.0_dp
        bornq(3,3,i) = qf(i)*occuf(i)/dble(nprocs)
      else
        bornq(1:3,1:3,i) = 0.0_dp
      endif
    enddo
    if (nshell.gt.0) then
!
!  Born charges equal normal charges if there are no shells present
!
!  Collect shell second derivative terms
!
      do i = 1,nshellr1
        ni = nshptr(i)
        qshell(i) = qf(ni)*occuf(ni)
      enddo
      msvar = 3*nshellr1
      nshellr1onnode = 0
!
      do i = 1,nshellonnode
        ii = nshptr(i)
        if (nregionno(nsft+nrelf2a(ii)).eq.1) then
          nshellr1onnode = nshellr1onnode + 1
          nshellr1ptr(nshellr1onnode) = i
          is = nshonnodeptr(i)
          ix = 3*(is-1) + 1
          iy = ix + 1
          iz = ix + 2
          ixs = 3*(nshellr1onnode-1) + 1
          iys = ixs + 1
          izs = ixs + 2
          do j = 1,nshellr1
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
        endif
      enddo
!
      if (lpocc) then
        call compressd1(qshell,0_i4,nsfoc,nshellr1,iocshptr)
        if (nprocs.gt.1) then
          call outerror('compressd2 not yet parallelised',0_i4)
          call stopnow('property3')
        endif
        call compressd2(dervi,maxd2,0_i4,nsfoc,0_i4,nshellr1,iocshptr,ibocshptr)
        n = 3*nsfoc
      else
        n = msvar
      endif
!************************************
!  Invert second derivative matrix  *
!************************************
      ifail = 0
      call matrix_inversion_shells(n,1_i4,maxd2,dervi,nshellr1,nshellr1onnode,ifail)
!
      if (ifail.gt.0) then
        call outwarning('Born charges cannot be calculated - matrix is singular',0_i4)
        goto 999
      endif
!
!  Multiply inverse shell second derivatives by core-shell second derivatives
!
      if (nprocs.gt.1) then
        ix = - 2
        iy = - 1
        iz =   0
        do i = 1,numat
          ix = ix + 3
          iy = iy + 3
          iz = iz + 3
          if (nregionno(nsft+nrelf2a(i)).eq.1.and.nat(i).le.maxele) then
            do j = 1,nshellr1onnode
              jj = nshellr1ptr(j)
              nj = nshonnodeptr(jj)
              jx = 3*(jj - 1) + 1
              jy = jx + 1
              jz = jy + 1
              jxs = 3*(nj - 1) + 1
              jys = jxs + 1
              jzs = jxs + 2
              kx = - 2
              ky = - 1
              kz =   0
              do k = 1,nshellr1
                kx = kx + 3
                ky = ky + 3
                kz = kz + 3
                bornq(1,1,i) = bornq(1,1,i) -  &
                  qshell(k)*(dervi(kx,jx)*derv2(ix,jxs) + &
                             dervi(kx,jy)*derv2(ix,jys) + &
                             dervi(kx,jz)*derv2(ix,jzs))
                bornq(1,2,i) = bornq(1,2,i) -  &
                  qshell(k)*(dervi(ky,jx)*derv2(ix,jxs) + &
                             dervi(ky,jy)*derv2(ix,jys) + &
                             dervi(ky,jz)*derv2(ix,jzs))
                bornq(1,3,i) = bornq(1,3,i) -  &
                  qshell(k)*(dervi(kz,jx)*derv2(ix,jxs) + &
                             dervi(kz,jy)*derv2(ix,jys) + &
                             dervi(kz,jz)*derv2(ix,jzs))
                bornq(2,1,i) = bornq(2,1,i) -  &
                  qshell(k)*(dervi(kx,jx)*derv2(iy,jxs) + &
                             dervi(kx,jy)*derv2(iy,jys) + &
                             dervi(kx,jz)*derv2(iy,jzs))
                bornq(2,2,i) = bornq(2,2,i) -  &
                  qshell(k)*(dervi(ky,jx)*derv2(iy,jxs) + &
                             dervi(ky,jy)*derv2(iy,jys) + &
                             dervi(ky,jz)*derv2(iy,jzs))
                bornq(2,3,i) = bornq(2,3,i) -  &
                  qshell(k)*(dervi(kz,jx)*derv2(iy,jxs) + &
                             dervi(kz,jy)*derv2(iy,jys) + &
                             dervi(kz,jz)*derv2(iy,jzs))
                bornq(3,1,i) = bornq(3,1,i) -  &
                  qshell(k)*(dervi(kx,jx)*derv2(iz,jxs) + &
                             dervi(kx,jy)*derv2(iz,jys) + &
                             dervi(kx,jz)*derv2(iz,jzs))
                bornq(3,2,i) = bornq(3,2,i) -  &
                  qshell(k)*(dervi(ky,jx)*derv2(iz,jxs) + &
                             dervi(ky,jy)*derv2(iz,jys) + &
                             dervi(ky,jz)*derv2(iz,jzs))
                bornq(3,3,i) = bornq(3,3,i) -  &
                  qshell(k)*(dervi(kz,jx)*derv2(iz,jxs) + &
                             dervi(kz,jy)*derv2(iz,jys) + &
                             dervi(kz,jz)*derv2(iz,jzs))
              enddo
            enddo
          endif
        enddo
      else
        if (lrigid) then
          ix = - 2
          iy = - 1
          iz =   0
          do i = 1,numat
            ix = ix + 3
            iy = iy + 3
            iz = iz + 3
            if (nregionno(nsft+nrelf2a(i)).eq.1.and.nat(i).le.maxele) then
              do j = 1,nshellr1
                nj = nshptr(j)
                jx = 3*(j - 1) + 1
                jy = jx + 1
                jz = jx + 2
                jxs = 3*(nj - 1) + 1
                jys = jxs + 1
                jzs = jxs + 2
                kx = - 2
                ky = - 1
                kz =   0
                do k = 1,nshellr1
                  kx = kx + 3
                  ky = ky + 3
                  kz = kz + 3
                  bornq(1,1,i) = bornq(1,1,i) -  &
                    qshell(k)*(dervi(kx,jx)*derv2(jxs,ix) + &
                               dervi(kx,jy)*derv2(jys,ix) + &
                               dervi(kx,jz)*derv2(jzs,ix))
                  bornq(1,2,i) = bornq(1,2,i) -  &
                    qshell(k)*(dervi(ky,jx)*derv2(jxs,ix) + &
                               dervi(ky,jy)*derv2(jys,ix) + &
                               dervi(ky,jz)*derv2(jzs,ix))
                  bornq(1,3,i) = bornq(1,3,i) -  &
                    qshell(k)*(dervi(kz,jx)*derv2(jxs,ix) + &
                               dervi(kz,jy)*derv2(jys,ix) + &
                               dervi(kz,jz)*derv2(jzs,ix))
                  bornq(2,1,i) = bornq(2,1,i) -  &
                    qshell(k)*(dervi(kx,jx)*derv2(jxs,iy) + &
                               dervi(kx,jy)*derv2(jys,iy) + &
                               dervi(kx,jz)*derv2(jzs,iy))
                  bornq(2,2,i) = bornq(2,2,i) -  &
                    qshell(k)*(dervi(ky,jx)*derv2(jxs,iy) + &
                               dervi(ky,jy)*derv2(jys,iy) + &
                               dervi(ky,jz)*derv2(jzs,iy))
                  bornq(2,3,i) = bornq(2,3,i) -  &
                    qshell(k)*(dervi(kz,jx)*derv2(jxs,iy) + &
                               dervi(kz,jy)*derv2(jys,iy) + &
                               dervi(kz,jz)*derv2(jzs,iy))
                  bornq(3,1,i) = bornq(3,1,i) -  &
                    qshell(k)*(dervi(kx,jx)*derv2(jxs,iz) + &
                               dervi(kx,jy)*derv2(jys,iz) + &
                               dervi(kx,jz)*derv2(jzs,iz))
                  bornq(3,2,i) = bornq(3,2,i) -  &
                    qshell(k)*(dervi(ky,jx)*derv2(jxs,iz) + &
                               dervi(ky,jy)*derv2(jys,iz) + &
                               dervi(ky,jz)*derv2(jzs,iz))
                  bornq(3,3,i) = bornq(3,3,i) -  &
                    qshell(k)*(dervi(kz,jx)*derv2(jxs,iz) + &
                               dervi(kz,jy)*derv2(jys,iz) + &
                               dervi(kz,jz)*derv2(jzs,iz))
                enddo
              enddo
            endif
          enddo
        else
          ix = - 2
          iy = - 1
          iz =   0
          do i = 1,numat
            ix = ix + 3
            iy = iy + 3
            iz = iz + 3
            if (nregionno(nsft+nrelf2a(i)).eq.1.and.nat(i).le.maxele) then
              do j = 1,nshellr1
                nj = nshptr(j)
                jfoc = iocshptr(j)
                jx = 3*(jfoc - 1) + 1
                jy = jx + 1
                jz = jy + 1
                jxs = 3*(nj - 1) + 1
                jys = jxs + 1
                jzs = jxs + 2
                kx = - 2
                ky = - 1
                kz =   0
                do k = 1,nsfoc
                  kx = kx + 3
                  ky = ky + 3
                  kz = kz + 3
                  bornq(1,1,i) = bornq(1,1,i) -  &
                    qshell(k)*(dervi(kx,jx)*derv2(jxs,ix) + &
                               dervi(kx,jy)*derv2(jys,ix) + &
                               dervi(kx,jz)*derv2(jzs,ix))
                  bornq(1,2,i) = bornq(1,2,i) -  &
                    qshell(k)*(dervi(ky,jx)*derv2(jxs,ix) + &
                               dervi(ky,jy)*derv2(jys,ix) + &
                               dervi(ky,jz)*derv2(jzs,ix))
                  bornq(1,3,i) = bornq(1,3,i) -  &
                    qshell(k)*(dervi(kz,jx)*derv2(jxs,ix) + &
                               dervi(kz,jy)*derv2(jys,ix) + &
                               dervi(kz,jz)*derv2(jzs,ix))
                  bornq(2,1,i) = bornq(2,1,i) -  &
                    qshell(k)*(dervi(kx,jx)*derv2(jxs,iy) + &
                               dervi(kx,jy)*derv2(jys,iy) + &
                               dervi(kx,jz)*derv2(jzs,iy))
                  bornq(2,2,i) = bornq(2,2,i) -  &
                    qshell(k)*(dervi(ky,jx)*derv2(jxs,iy) + &
                               dervi(ky,jy)*derv2(jys,iy) + &
                               dervi(ky,jz)*derv2(jzs,iy))
                  bornq(2,3,i) = bornq(2,3,i) -  &
                    qshell(k)*(dervi(kz,jx)*derv2(jxs,iy) + &
                               dervi(kz,jy)*derv2(jys,iy) + &
                               dervi(kz,jz)*derv2(jzs,iy))
                  bornq(3,1,i) = bornq(3,1,i) -  &
                    qshell(k)*(dervi(kx,jx)*derv2(jxs,iz) + &
                               dervi(kx,jy)*derv2(jys,iz) + &
                               dervi(kx,jz)*derv2(jzs,iz))
                  bornq(3,2,i) = bornq(3,2,i) -  &
                    qshell(k)*(dervi(ky,jx)*derv2(jxs,iz) + &
                               dervi(ky,jy)*derv2(jys,iz) + &
                               dervi(ky,jz)*derv2(jzs,iz))
                  bornq(3,3,i) = bornq(3,3,i) -  &
                    qshell(k)*(dervi(kz,jx)*derv2(jxs,iz) + &
                               dervi(kz,jy)*derv2(jys,iz) + &
                               dervi(kz,jz)*derv2(jzs,iz))
                enddo
              enddo
            endif
          enddo
        endif
      endif
    endif
!**************************************************************
!  Charge equilibration scheme contributions to Born charges  *
!**************************************************************
    if (leem) then
      do ii = 1,natomsonnode
        i = node2atom(ii)
        if (nregionno(nsft+nrelf2a(i)).eq.1) then
          kx = - 2
          ky = - 1
          kz =   0
          xi = xclat(i)
          yi = yclat(i)
          zi = zclat(i)
          do k = 1,numatr1
            kx = kx + 3
            ky = ky + 3
            kz = kz + 3
            bornq(1,1,k) = bornq(1,1,k) + dqdxyz(kx,ii)*xi
            bornq(2,1,k) = bornq(2,1,k) + dqdxyz(kx,ii)*yi
            bornq(3,1,k) = bornq(3,1,k) + dqdxyz(kx,ii)*zi
            bornq(1,2,k) = bornq(1,2,k) + dqdxyz(ky,ii)*xi
            bornq(2,2,k) = bornq(2,2,k) + dqdxyz(ky,ii)*yi
            bornq(3,2,k) = bornq(3,2,k) + dqdxyz(ky,ii)*zi
            bornq(1,3,k) = bornq(1,3,k) + dqdxyz(kz,ii)*xi
            bornq(2,3,k) = bornq(2,3,k) + dqdxyz(kz,ii)*yi
            bornq(3,3,k) = bornq(3,3,k) + dqdxyz(kz,ii)*zi
          enddo
        endif
      enddo
    endif
!**********************************
!  Global sum of bornq if needed  *
!**********************************
    if (nprocs.gt.1) then
      allocate(bqtmp(3,3,numat),stat=status)
      if (status/=0) call outofmemory('borncharge','bqtmp')
!
      call sumall(bornq,bqtmp,9_i4*numat,"borncharge","bornq")
      bornq(1:3,1:3,1:numat) = bqtmp(1:3,1:3,1:numat)
!
      deallocate(bqtmp,stat=status)
      if (status/=0) call deallocate_error('borncharge','bqtmp')
    endif
  endif
!**********************
!  Output properties  *
!**********************
  if (lprint.and.ioproc) then
    if (lbornloc) then
      if (lcml) call gulp_cml_print_born_charges(bornq)
      write(ioout,'(/,''  Born effective charge tensors : '',/)')
      write(ioout,'(''-------------------------------------------------------------------------------'')')
      write(ioout,'(''  Atom             x           y             z'')')
      write(ioout,'(''-------------------------------------------------------------------------------'')')
      n = 0
      if (ndim.eq.2) then
        do i = 1,numat
          inat = nat(i)  
          if (nregionno(nsft+nrelf2a(i)).eq.1.and.inat.le.maxele) then
            n = n + 1
            itype = nftype(i)
            call label(inat,itype,lab)
            write(ioout,'(i5,1x,a5,1x,''x '',3(2x,f10.4))') n,lab,bornq(1,1,i),bornq(2,1,i),bornq(3,1,i)
            write(ioout,'(12x,''y '',3(2x,f10.4))') bornq(1,2,i),bornq(2,2,i),bornq(3,2,i)
            write(ioout,'(12x,''z '',3(2x,f10.4))') bornq(1,3,i),bornq(2,3,i),bornq(3,3,i)
            write(ioout,'(''-------------------------------------------------------------------------------'')')
          endif
        enddo
      else
        do i = 1,nasym
          ii = nrela2f(i)
          inat = nat(ii)  
          if (inat.le.maxele) then
            n = n + 1
            itype = nftype(ii)
            call label(inat,itype,lab)
            write(ioout,'(i5,1x,a5,1x,''x '',3(2x,f10.4))') n,lab,bornq(1,1,ii),bornq(2,1,ii),bornq(3,1,ii)
            write(ioout,'(12x,''y '',3(2x,f10.4))') bornq(1,2,ii),bornq(2,2,ii),bornq(3,2,ii)
            write(ioout,'(12x,''z '',3(2x,f10.4))') bornq(1,3,ii),bornq(2,3,ii),bornq(3,3,ii)
            write(ioout,'(''-------------------------------------------------------------------------------'')')
          endif
        enddo
      endif
      write(ioout,'(/)')
    endif
    call gflush(ioout)
  endif
!***************
!  Exit tasks  *
!***************
999 continue
!
!  Timings
!
  t2p = g_cpu_time()
  tprop = t2p - t1p + tprop
!
!  Deallocate memory
!
  if (nshellr1.gt.0) then
    deallocate(nshellr1ptr,stat=status)
    if (status/=0) call deallocate_error('borncharge','nshellr1ptr')
  endif
  deallocate(qshell,stat=status)
  if (status/=0) call deallocate_error('borncharge','qshell')
#ifdef TRACE
  call trace_out('borncharge')
#endif
!
  return
  end
