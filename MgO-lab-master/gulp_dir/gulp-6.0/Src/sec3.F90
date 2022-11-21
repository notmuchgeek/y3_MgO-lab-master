  subroutine sec3
!
!  Generate symmetry adapted second derivative matrix in
!  fractional units and strain from cartesian matrix for
!  the full unit cell.
!
!   5/95 Modifications now added for symmetry adapted second
!        derivative matrix.
!   8/95 Modification added to avoid use of transformation
!        unless necessary as this slows jobs down a lot.
!  11/96 Freezing removed as sec3f is a separate routine now
!   2/01 Modifications for general dimensionality added
!  10/03 Modified to handle cell parameter variables
!   1/04 Typos in rotation of breathing shells corrected
!   1/04 Modification of nuclear-radius cross terms
!        corrected in derv2 & radius-radius terms
!   1/04 Rotation of coordinate - radius block introduced
!        for symmetry adapted case
!  11/06 Cell parameter optimisation added
!  12/07 Arguments to setcellderv1D corrected
!   3/17 fix_atom option added
!   3/17 Modifications made to allow for new variable order in iopt
!   2/18 Trace added
!   4/18 Bug corrected for case where nfixatom is zero due to flags
!   3/19 iopt replaced by ioptindex and iopttype
!   4/19 Format for symmetrised second derivatives updated
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!  12/19 Rigid molecule modifications added
!   3/20 Changes to ninternal to separate rigid molecule variables
!   4/20 Non-orthogonal cell modifications added for rigid molecules
!   5/20 Corrections to molQTdrv handling
!   5/20 Trapping of molecules with no on-diagonal molQQdrv added
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
  use configurations, only : lbsmat
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
  integer(i4)                                  :: i
  integer(i4)                                  :: ibs
  integer(i4)                                  :: id
  integer(i4)                                  :: ii
  integer(i4)                                  :: indi
  integer(i4)                                  :: indii
  integer(i4)                                  :: indj
  integer(i4)                                  :: indk
  integer(i4)                                  :: io
  integer(i4)                                  :: ix
  integer(i4)                                  :: iy
  integer(i4)                                  :: iz
  integer(i4)                                  :: j
  integer(i4)                                  :: jd
  integer(i4)                                  :: jj
  integer(i4)                                  :: jjx
  integer(i4)                                  :: jk
  integer(i4)                                  :: jo
  integer(i4)                                  :: jx
  integer(i4)                                  :: jy
  integer(i4)                                  :: jz
  integer(i4)                                  :: k
  integer(i4)                                  :: kk
  integer(i4)                                  :: l
  integer(i4)                                  :: n
  integer(i4)                                  :: ni
  integer(i4)                                  :: nj
  integer(i4)                                  :: nm
  integer(i4)                                  :: nma
  integer(i4)                                  :: nmf
  integer(i4)                                  :: nn
  integer(i4)                                  :: n3a
  integer(i4)                                  :: n3f
  integer(i4)                                  :: n3ma
  integer(i4)                                  :: n3mf
  integer(i4)                                  :: ncellfull
  integer(i4)                                  :: neqj
  integer(i4)                                  :: nintmin
  integer(i4)                                  :: nloop
  integer(i4)                                  :: nofa
  integer(i4)                                  :: noff
  integer(i4)                                  :: status
  logical                                      :: lfullra
  logical                                      :: lmolQi
  logical                                      :: lmolTi
  logical                                      :: lmolQj
  logical                                      :: lmolTj
  logical                                      :: lsdebug
  logical                                      :: ltmat
  real(dp)                                     :: celldrv2(6,6)
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: d2(3,3)
  real(dp)                                     :: d2x
  real(dp)                                     :: d2y
  real(dp)                                     :: d2z
  real(dp)                                     :: d3tmp(6)
  real(dp)                                     :: dum
  real(dp)                                     :: dx
  real(dp)                                     :: dy
  real(dp)                                     :: dz
  real(dp)                                     :: r11
  real(dp)                                     :: r12
  real(dp)                                     :: r13
  real(dp)                                     :: r22
  real(dp)                                     :: r23
  real(dp)                                     :: r33
  real(dp)                                     :: r1xl
  real(dp)                                     :: r2yl
  real(dp)                                     :: r3zl
  real(dp)                                     :: rmat(3,3)
  real(dp)                                     :: sum
  real(dp)                                     :: t1
  real(dp)                                     :: t2
  real(dp)                                     :: tmp1(6,6)
  real(dp), dimension(:), allocatable          :: tmp2
#ifdef TRACE
  call trace_in('sec3')
#endif
!
  t1 = g_cpu_time()
  lsdebug = (index(keyword,'derv2').ne.0)
  n3f = 3*numat
  n3a = 3*nasym
  noff = n3f
  nofa = n3a
  if (nbsm.gt.0) then
    n3f = n3f + numat
    n3a = n3a + nasym
  endif
!
  if (lrigid) then
    nma = nmolasym
    nmf = nmol
    n3ma = 3*nma
    n3mf = 3*nmf
  endif
!
!  Work out whether full tmat multiplication is needed - use faster method to handle P1 fully relaxed case
!
  if (ncon.eq.0.and.(nspcg(ncf).le.1.and.ngocfg(ncf).eq.1)) then
    ltmat = .false.
  elseif (ninternal.eq.0) then
    ltmat = .false.
  else
    ltmat = .true.
  endif
  if (lsymderv2) then
    nloop = nasym
  else
    nloop = numat
  endif
!
!  Set up cell scaling factors
!
  if (ndim.eq.3) then
    do i = 1,3
      rmat(i,1) = rv(i,1)
      rmat(i,2) = rv(i,2)
      rmat(i,3) = rv(i,3)
    enddo
  elseif (ndim.eq.2) then
    do i = 1,2
      rmat(i,1) = rv(i,1)
      rmat(i,2) = rv(i,2)
      rmat(i,3) = 0.0_dp
      rmat(3,i) = 0.0_dp
    enddo
    rmat(3,3) = 1.0_dp
  elseif (ndim.eq.1) then
    do i = 1,3
      rmat(i,1) = 0.0_dp
      rmat(i,2) = 0.0_dp
      rmat(i,3) = 0.0_dp
    enddo
    rmat(1,1) = rv(1,1)
    rmat(2,2) = 1.0_dp
    rmat(3,3) = 1.0_dp
  endif
!
!  Generate full cell vectors. Also if full cell is right angled then use 
!  speed up tricks for this special case (lfullra=.true.)
!
  if (ncbl.gt.1) then
    call uncentre(rmat)
    sum = abs(rv(1,2)) + abs(rv(2,1)) + abs(rv(1,3)) + abs(rv(3,1)) + abs(rv(2,3)) + abs(rv(3,2))
    lfullra = (sum.lt.1.0d-6)
  else
    lfullra = lra
  endif
!
!  Print out full second derivatives
!
  if (lsdebug.and.ioproc) then
    write(ioout,'(/,''  Second Derivative Matrix  :'',/)')
    do i = 1,n3f
      write(ioout,'(2x,9(f9.4))')(derv2(i,j),j=1,(3*nloop+nbsmat))
    enddo
    write(ioout,'(/)')
    if (lrigid) then
      write(ioout,'(/,''  Second Derivative Matrix : Rigid molecule translation-translation :'',/)')
      do i = 1,nmf
        do j = 1,nmf
          write(ioout,'(2x,2i5,'' x'',1x,3(f12.4))') i,j,(molTTdrv(1,jj,i,j),jj=1,3)
          write(ioout,'(2x,2i5,'' y'',1x,3(f12.4))') i,j,(molTTdrv(2,jj,i,j),jj=1,3)
          write(ioout,'(2x,2i5,'' z'',1x,3(f12.4))') i,j,(molTTdrv(3,jj,i,j),jj=1,3)
        enddo
      enddo
      write(ioout,'(/,''  Second Derivative Matrix : Rigid molecule quaternion-quaternion :'',/)')
      do i = 1,nmf
        do j = 1,nmf
          write(ioout,'(2x,2i5,'' x'',1x,3(f12.4))') i,j,(molQQdrv(1,jj,i,j),jj=1,3)
          write(ioout,'(2x,2i5,'' y'',1x,3(f12.4))') i,j,(molQQdrv(2,jj,i,j),jj=1,3)
          write(ioout,'(2x,2i5,'' z'',1x,3(f12.4))') i,j,(molQQdrv(3,jj,i,j),jj=1,3)
        enddo
      enddo
      write(ioout,'(/,''  Second Derivative Matrix : Rigid molecule quaternion-translation :'',/)')
      do i = 1,nmf
        do j = 1,nmf
          write(ioout,'(2x,2i5,'' x'',1x,3(f12.4))') i,j,(molQTdrv(1,jj,i,j),jj=1,3)
          write(ioout,'(2x,2i5,'' y'',1x,3(f12.4))') i,j,(molQTdrv(2,jj,i,j),jj=1,3)
          write(ioout,'(2x,2i5,'' z'',1x,3(f12.4))') i,j,(molQTdrv(3,jj,i,j),jj=1,3)
        enddo
      enddo
      if (ndim.gt.0) then
        write(ioout,'(/,''  Second Derivative Matrix : Rigid molecule translation-strain :'',/)')
        do i = 1,nmf
          do j = 1,nstrains
            write(ioout,'(2x,2i5,1x,3(f12.4))') i,j,(molTSdrv(j,jj,i),jj=1,3)
          enddo
        enddo
        write(ioout,'(/,''  Second Derivative Matrix : Rigid molecule quaternion-strain :'',/)')
        do i = 1,nmf
          do j = 1,nstrains
            write(ioout,'(2x,2i5,1x,3(f12.4))') i,j,(molQSdrv(j,jj,i),jj=1,3)
          enddo
        enddo
      endif
      write(ioout,'(/)')
    endif
  endif
!********************************************
!  Rigid molecules - check for linear case  *
!********************************************
  if (lrigid) then
    do i = 1,nmf
      do ix = 1,3
!
!  If molecule is linear or there are no rotational forces that set the diagonal element to be a dummy value 
!  so that the Hessian can be inverted
!
        if (abs(molQQdrv(ix,ix,i,i)).lt.0.001_dp) then
          do jx = 1,3
            molQQdrv(jx,ix,i,i) = 0.0_dp
            molQQdrv(ix,jx,i,i) = 0.0_dp
          enddo
          molQQdrv(ix,ix,i,i) = 1.0_dp
        endif
      enddo
    enddo
  endif
!***********************************************************************
!  Transform internal cartesian second derivatives to internal system  *
!  Needed in all algorithms                                            *
!***********************************************************************
  if (lfullra) then
    r1xl = rmat(1,1)
    r2yl = rmat(2,2)
    r3zl = rmat(3,3)
    r11 = r1xl*r1xl
    r12 = r1xl*r2yl
    r13 = r1xl*r3zl
    r22 = r2yl*r2yl
    r23 = r2yl*r3zl
    r33 = r3zl*r3zl
    do i = 1,nloop
      indi = 3*(i-1)
!
!  Coordinate
!
      do j = 1,numat
        indj = 3*(j-1)
        derv2(indj+1,indi+1) = r11*derv2(indj+1,indi+1)
        derv2(indj+1,indi+2) = r12*derv2(indj+1,indi+2)
        derv2(indj+1,indi+3) = r13*derv2(indj+1,indi+3)
        derv2(indj+2,indi+1) = r12*derv2(indj+2,indi+1)
        derv2(indj+2,indi+2) = r22*derv2(indj+2,indi+2)
        derv2(indj+2,indi+3) = r23*derv2(indj+2,indi+3)
        derv2(indj+3,indi+1) = r13*derv2(indj+3,indi+1)
        derv2(indj+3,indi+2) = r23*derv2(indj+3,indi+2)
        derv2(indj+3,indi+3) = r33*derv2(indj+3,indi+3)
      enddo
!
!  Radial
!
      if (lsymderv2) then
        ii = i
      else
        ii = nrelf2a(i)
      endif
      if (lbsmat(nsft+ii)) then
        do j = 1,numat
          derv2(noff+j,indi+1) = r1xl*derv2(noff+j,indi+1)
          derv2(noff+j,indi+2) = r2yl*derv2(noff+j,indi+2)
          derv2(noff+j,indi+3) = r3zl*derv2(noff+j,indi+3)
        enddo
        indii = 3*nloop + i
        do j = 1,numat
          jj = 3*(j-1)
          jx = jj + 1
          jy = jj + 2
          jz = jj + 3
          derv2(jx,indii) = r1xl*derv2(jx,indii)
          derv2(jy,indii) = r2yl*derv2(jy,indii)
          derv2(jz,indii) = r3zl*derv2(jz,indii)
        enddo
      endif
    enddo
!
!  Rigid molecule derivatives
!
    if (lrigid) then
      do nm = 1,nmol
        do nn = 1,nmol
!
!  Molecule - molecule
!
          molTTdrv(1,1,nn,nm) = molTTdrv(1,1,nn,nm)*r11
          molTTdrv(2,1,nn,nm) = molTTdrv(2,1,nn,nm)*r12
          molTTdrv(3,1,nn,nm) = molTTdrv(3,1,nn,nm)*r13
          molTTdrv(1,2,nn,nm) = molTTdrv(1,2,nn,nm)*r12
          molTTdrv(2,2,nn,nm) = molTTdrv(2,2,nn,nm)*r22
          molTTdrv(3,2,nn,nm) = molTTdrv(3,2,nn,nm)*r23
          molTTdrv(1,3,nn,nm) = molTTdrv(1,3,nn,nm)*r13
          molTTdrv(2,3,nn,nm) = molTTdrv(2,3,nn,nm)*r23
          molTTdrv(3,3,nn,nm) = molTTdrv(3,3,nn,nm)*r33
!
          molQTdrv(1,1,nn,nm) = molQTdrv(1,1,nn,nm)*r1xl
          molQTdrv(2,1,nn,nm) = molQTdrv(2,1,nn,nm)*r1xl
          molQTdrv(3,1,nn,nm) = molQTdrv(3,1,nn,nm)*r1xl
          molQTdrv(1,2,nn,nm) = molQTdrv(1,2,nn,nm)*r2yl
          molQTdrv(2,2,nn,nm) = molQTdrv(2,2,nn,nm)*r2yl
          molQTdrv(3,2,nn,nm) = molQTdrv(3,2,nn,nm)*r2yl
          molQTdrv(1,3,nn,nm) = molQTdrv(1,3,nn,nm)*r3zl
          molQTdrv(2,3,nn,nm) = molQTdrv(2,3,nn,nm)*r3zl
          molQTdrv(3,3,nn,nm) = molQTdrv(3,3,nn,nm)*r3zl
        enddo
!
!  Coordinate - molecule mixed
!
        do i = 1,nloop
          indi = 3*(i-1)
          molTCdrv(indi+1,1,nm) = molTCdrv(indi+1,1,nm)*r11
          molTCdrv(indi+2,1,nm) = molTCdrv(indi+2,1,nm)*r12
          molTCdrv(indi+3,1,nm) = molTCdrv(indi+3,1,nm)*r13
          molTCdrv(indi+1,2,nm) = molTCdrv(indi+1,2,nm)*r12
          molTCdrv(indi+2,2,nm) = molTCdrv(indi+2,2,nm)*r22
          molTCdrv(indi+3,2,nm) = molTCdrv(indi+3,2,nm)*r23
          molTCdrv(indi+1,3,nm) = molTCdrv(indi+1,3,nm)*r13
          molTCdrv(indi+2,3,nm) = molTCdrv(indi+2,3,nm)*r23
          molTCdrv(indi+3,3,nm) = molTCdrv(indi+3,3,nm)*r33
!
          molQCdrv(indi+1,1,nm) = molQCdrv(indi+1,1,nm)*r1xl
          molQCdrv(indi+2,1,nm) = molQCdrv(indi+2,1,nm)*r1xl
          molQCdrv(indi+3,1,nm) = molQCdrv(indi+3,1,nm)*r1xl
          molQCdrv(indi+1,2,nm) = molQCdrv(indi+1,2,nm)*r2yl
          molQCdrv(indi+2,2,nm) = molQCdrv(indi+2,2,nm)*r2yl
          molQCdrv(indi+3,2,nm) = molQCdrv(indi+3,2,nm)*r2yl
          molQCdrv(indi+1,3,nm) = molQCdrv(indi+1,3,nm)*r3zl
          molQCdrv(indi+2,3,nm) = molQCdrv(indi+2,3,nm)*r3zl
          molQCdrv(indi+3,3,nm) = molQCdrv(indi+3,3,nm)*r3zl
        enddo
      enddo
    endif
  else
    if (lsymderv2) then
      do i = 1,nasym
        indi = 3*(i - 1)
!
!  Coordinates
!
        do j = 1,numat
          indj = 3*(j - 1)
          do ii = 1,3
            do jj = 1,3
              tmp1(jj,ii) = derv2(indj+1,indi+ii)*rmat(1,jj) + &
                            derv2(indj+2,indi+ii)*rmat(2,jj) + &
                            derv2(indj+3,indi+ii)*rmat(3,jj)
            enddo
          enddo
          do ii = 1,3
            do jj = 1,3
              derv2(indj+jj,indi+ii) = rmat(1,ii)*tmp1(jj,1) + &
                                       rmat(2,ii)*tmp1(jj,2) + &
                                       rmat(3,ii)*tmp1(jj,3)
            enddo
          enddo
        enddo
!
!  Radial
!
        if (lbsmat(nsft+i)) then
          do j = 1,numat
            dx = derv2(noff+j,indi+1)
            dy = derv2(noff+j,indi+2)
            dz = derv2(noff+j,indi+3)
            derv2(noff+j,indi+1) = dx*rmat(1,1) + dy*rmat(2,1) + dz*rmat(3,1)
            derv2(noff+j,indi+2) = dx*rmat(1,2) + dy*rmat(2,2) + dz*rmat(3,2)
            derv2(noff+j,indi+3) = dx*rmat(1,3) + dy*rmat(2,3) + dz*rmat(3,3)
          enddo
          indii = 3*nasym + i
          do j = 1,numat
            jj = 3*(j - 1)
            jx = jj + 1
            jy = jj + 2
            jz = jj + 3
            dx = derv2(jx,indii)
            dy = derv2(jy,indii)
            dz = derv2(jz,indii)
            derv2(jx,indii) = dx*rmat(1,1) + dy*rmat(2,1) + dz*rmat(3,1)
            derv2(jy,indii) = dx*rmat(1,2) + dy*rmat(2,2) + dz*rmat(3,2)
            derv2(jz,indii) = dx*rmat(1,3) + dy*rmat(2,3) + dz*rmat(3,3)
          enddo
        endif
      enddo
!
!  Rigid molecule derivatives
!
      if (lrigid) then
        do nm = 1,nmolasym
          do nn = 1,nmol
!
!  Molecule - molecule
!
            do ii = 1,3
              do jj = 1,3
                tmp1(jj,ii) = molTTdrv(1,ii,nn,nm)*rmat(1,jj) + &
                              molTTdrv(2,ii,nn,nm)*rmat(2,jj) + &
                              molTTdrv(3,ii,nn,nm)*rmat(3,jj)
              enddo
            enddo
            do ii = 1,3
              do jj = 1,3
                molTTdrv(jj,ii,nn,nm) = rmat(1,ii)*tmp1(jj,1) + &
                                        rmat(2,ii)*tmp1(jj,2) + &
                                        rmat(3,ii)*tmp1(jj,3)
              enddo
            enddo
          enddo
!
          do ii = 1,3
            dx = molQTdrv(1,ii,nn,nm)
            dy = molQTdrv(2,ii,nn,nm)
            dz = molQTdrv(3,ii,nn,nm)
            molQTdrv(1,ii,nn,nm) = dx*rmat(1,1) + dy*rmat(2,1) + dz*rmat(3,1)
            molQTdrv(2,ii,nn,nm) = dx*rmat(1,2) + dy*rmat(2,2) + dz*rmat(3,2)
            molQTdrv(3,ii,nn,nm) = dx*rmat(1,3) + dy*rmat(2,3) + dz*rmat(3,3)
          enddo
!
!  Coordinate - molecule mixed
!
          do i = 1,nloop
            indi = 3*(i-1)
            do ii = 1,3
              do jj = 1,3
                tmp1(jj,ii) = molTCdrv(indi+1,ii,nm)*rmat(1,jj) + &
                              molTCdrv(indi+2,ii,nm)*rmat(2,jj) + &
                              molTCdrv(indi+3,ii,nm)*rmat(3,jj)
              enddo
            enddo
            do ii = 1,3
              do jj = 1,3
                molTCdrv(indi+jj,ii,nm) = rmat(1,ii)*tmp1(jj,1) + &
                                          rmat(2,ii)*tmp1(jj,2) + &
                                          rmat(3,ii)*tmp1(jj,3)
              enddo
            enddo
!
            do ii = 1,3
              dx = molQCdrv(indi+1,ii,nm)
              dy = molQCdrv(indi+2,ii,nm)
              dz = molQCdrv(indi+3,ii,nm)
              molQCdrv(indi+1,ii,nm) = dx*rmat(1,1) + dy*rmat(2,1) + dz*rmat(3,1)
              molQCdrv(indi+2,ii,nm) = dx*rmat(1,2) + dy*rmat(2,2) + dz*rmat(3,2)
              molQCdrv(indi+3,ii,nm) = dx*rmat(1,3) + dy*rmat(2,3) + dz*rmat(3,3)
            enddo
          enddo
        enddo
      endif
    else
      do i = 1,numat
        indi = 3*(i-1)
!
!  Coordinates
!
        do j = 1,i
          indj = 3*(j-1)
          do ii = 1,3
            do jj = 1,3
              tmp1(ii,jj) = derv2(indi+ii,indj+1)*rmat(1,jj) + &
                            derv2(indi+ii,indj+2)*rmat(2,jj) + &
                            derv2(indi+ii,indj+3)*rmat(3,jj)
            enddo
          enddo
          do ii = 1,3
            do jj = 1,3
              derv2(indi+ii,indj+jj) = rmat(1,ii)*tmp1(1,jj) + &
                                       rmat(2,ii)*tmp1(2,jj) + &
                                       rmat(3,ii)*tmp1(3,jj)
              derv2(indj+jj,indi+ii) = derv2(indi+ii,indj+jj)
            enddo
          enddo
        enddo
!
!  Radial
!
        if (lbsmat(nsft+nrelf2a(i))) then
          do j = 1,numat
            dx = derv2(noff+j,indi+1)
            dy = derv2(noff+j,indi+2)
            dz = derv2(noff+j,indi+3)
            derv2(noff+j,indi+1) = dx*rmat(1,1) + dy*rmat(2,1) + dz*rmat(3,1)
            derv2(noff+j,indi+2) = dx*rmat(1,2) + dy*rmat(2,2) + dz*rmat(3,2)
            derv2(noff+j,indi+3) = dx*rmat(1,3) + dy*rmat(2,3) + dz*rmat(3,3)
          enddo
          indii = 3*numat + i
          do j = 1,numat
            jj = 3*(j - 1)
            jx = jj + 1
            jy = jj + 2
            jz = jj + 3
            dx = derv2(jx,indii)
            dy = derv2(jy,indii)
            dz = derv2(jz,indii)
            derv2(jx,indii) = dx*rmat(1,1) + dy*rmat(2,1) + dz*rmat(3,1)
            derv2(jy,indii) = dx*rmat(1,2) + dy*rmat(2,2) + dz*rmat(3,2)
            derv2(jz,indii) = dx*rmat(1,3) + dy*rmat(2,3) + dz*rmat(3,3)
          enddo
        endif
      enddo
! 
!  Rigid molecule derivatives
! 
      if (lrigid) then
        do nm = 1,nmol
!
!  Molecule - molecule
!
          do nn = 1,nmol
            do ii = 1,3
              do jj = 1,3
                tmp1(jj,ii) = molTTdrv(1,ii,nn,nm)*rmat(1,jj) + &
                              molTTdrv(2,ii,nn,nm)*rmat(2,jj) + &
                              molTTdrv(3,ii,nn,nm)*rmat(3,jj)
              enddo
            enddo
            do ii = 1,3
              do jj = 1,3
                molTTdrv(jj,ii,nn,nm) = rmat(1,ii)*tmp1(jj,1) + &
                                        rmat(2,ii)*tmp1(jj,2) + &
                                        rmat(3,ii)*tmp1(jj,3)
              enddo
            enddo
!
            do ii = 1,3
              do jj = 1,3
                tmp1(jj,ii) = molQTdrv(jj,1,nn,nm)*rmat(1,ii) + &
                              molQTdrv(jj,2,nn,nm)*rmat(2,ii) + &
                              molQTdrv(jj,3,nn,nm)*rmat(3,ii)
              enddo
            enddo
            do ii = 1,3
              do jj = 1,3
                molQTdrv(jj,ii,nn,nm) = tmp1(jj,ii)
              enddo
            enddo
          enddo
!
!  Coordinate - molecule mixed
!
          do i = 1,nloop
            indi = 3*(i-1)
            do ii = 1,3
              do jj = 1,3
                tmp1(jj,ii) = molTCdrv(indi+1,ii,nm)*rmat(1,jj) + &
                              molTCdrv(indi+2,ii,nm)*rmat(2,jj) + &
                              molTCdrv(indi+3,ii,nm)*rmat(3,jj)
              enddo
            enddo
            do ii = 1,3
              do jj = 1,3
                molTCdrv(indi+jj,ii,nm) = rmat(1,ii)*tmp1(jj,1) + &
                                          rmat(2,ii)*tmp1(jj,2) + &
                                          rmat(3,ii)*tmp1(jj,3)
              enddo
            enddo
!
            do ii = 1,3
              do jj = 1,3
                tmp1(ii,jj) = molQCdrv(indi+1,jj,nm)*rmat(1,ii) + &
                              molQCdrv(indi+2,jj,nm)*rmat(2,ii) + &
                              molQCdrv(indi+3,jj,nm)*rmat(3,ii)
              enddo
            enddo
            do ii = 1,3
              do jj = 1,3
                molQCdrv(indi+jj,ii,nm) = tmp1(jj,ii)
              enddo
            enddo
          enddo
        enddo
      endif
    endif
  endif
!******************************************
!                                         *
!  Internal derivatives :                 *
!                                         *
!  Symmetry transform second derivatives  *
!******************************************
  if (ltmat) then
    if (.not.lsymderv2) then
!**************************
!  Numat-numat algorithm  *
!**************************
      allocate(tmp2(n3f),stat=status)
      if (status/=0) call outofmemory('sec3','tmp2')
!
!  Atom coordinate - atom coordinate
!
      do i = 1,ninternalatm
        do j = 1,n3f
          tmp2(j) = 0.0_dp
          do k = 1,n3f
            tmp2(j) = tmp2(j) + tmat(k,i)*derv2(j,k)
          enddo
        enddo
        do j = 1,ninternalatm
          derv2(j,i) = 0.0_dp
          do k = 1,n3f
            derv2(j,i) = derv2(j,i) + tmp2(k)*tmat(k,j)
          enddo
        enddo
      enddo
!
      if (lrigid) then
!------------------------------------------
!  Rigid molecule - rigid molecule terms  |
!------------------------------------------
!
!  Translation - translation
!
        do i = 1,ninternalmolT
          ii = ninternalatm + i
          tmp2(1:n3mf) = 0.0_dp
          do j = 1,nmf
            indj = 3*(j-1)
            do jj = 1,3
              do k = 1,nmf
                indk = 3*(k-1)
                do kk = 1,3
                  tmp2(indj+jj) = tmp2(indj+jj) + molTTdrv(jj,kk,j,k)*tmatT(indk+kk,i)
                enddo
              enddo
            enddo
          enddo
!
          do j = 1,ninternalmolT
            jj = ninternalatm + j
            derv2(jj,ii) = 0.0_dp
            do k = 1,nmf
              indk = 3*(k-1)
              do kk = 1,3
                derv2(jj,ii) = derv2(jj,ii) + tmp2(indk+kk)*tmatT(indk+kk,j)
              enddo
            enddo
          enddo
        enddo
!
!  Quaternion - translation
!
        do i = 1,ninternalmolT
          ii = ninternalatm + i
          tmp2(1:n3mf) = 0.0_dp
          do j = 1,nmf
            indj = 3*(j-1)
            do jj = 1,3
              do k = 1,nmf
                indk = 3*(k-1)
                do kk = 1,3
                  tmp2(indj+jj) = tmp2(indj+jj) + molQTdrv(jj,kk,j,k)*tmatT(indk+kk,i)
                enddo
              enddo
            enddo
          enddo
!
          do j = 1,ninternalmolQ
            jj = ninternalatm + ninternalmolT + j
            derv2(ii,jj) = 0.0_dp
            derv2(jj,ii) = 0.0_dp
            do k = 1,nmf
              indk = 3*(k-1)
              do kk = 1,3
                derv2(ii,jj) = derv2(ii,jj) + tmp2(indk+kk)*tmatQ(indk+kk,j)
                derv2(jj,ii) = derv2(jj,ii) + tmp2(indk+kk)*tmatQ(indk+kk,j)
              enddo
            enddo
          enddo
        enddo
!
!  Quaternion - quaternion
!
        do i = 1,ninternalmolQ
          ii = ninternalatm + ninternalmolT + i
          tmp2(1:n3mf) = 0.0_dp
          do j = 1,nmf
            indj = 3*(j-1)
            do jj = 1,3
              do k = 1,nmf
                indk = 3*(k-1)
                do kk = 1,3
                  tmp2(indj+jj) = tmp2(indj+jj) + molQQdrv(jj,kk,j,k)*tmatQ(indk+kk,i)
                enddo
              enddo
            enddo
          enddo
!
          do j = 1,ninternalmolQ
            jj = ninternalatm + ninternalmolT + j
            derv2(jj,ii) = 0.0_dp
            do k = 1,nmf
              indk = 3*(k-1)
              do kk = 1,3
                derv2(jj,ii) = derv2(jj,ii) + tmp2(indk+kk)*tmatQ(indk+kk,j)
              enddo
            enddo
          enddo
        enddo
!--------------------------------
!  Rigid molecule - atom terms  |
!--------------------------------
        do i = 1,ninternalatm
          tmp2(1:n3mf) = 0.0_dp
          do j = 1,nmf
            indj = 3*(j-1)
            do jj = 1,3
              do k = 1,n3f
                tmp2(indj+jj) = tmp2(indj+jj) + molTCdrv(k,jj,j)*tmat(k,i)
              enddo
            enddo
          enddo
          do j = 1,ninternalmolT
            jj = ninternalatm + j
            derv2(jj,i) = 0.0_dp
            derv2(i,jj) = 0.0_dp
            do k = 1,nmf
              indk = 3*(k-1)
              do kk = 1,3
                derv2(jj,i) = derv2(jj,i) + tmp2(indk+kk)*tmatT(indk+kk,j)
                derv2(i,jj) = derv2(i,jj) + tmp2(indk+kk)*tmatT(indk+kk,j)
              enddo
            enddo
          enddo
        enddo
!
        do i = 1,ninternalatm
          tmp2(1:n3mf) = 0.0_dp
          do j = 1,nmf
            indj = 3*(j-1)
            do jj = 1,3
              do k = 1,n3f
                tmp2(indj+jj) = tmp2(indj+jj) + molQCdrv(k,jj,j)*tmat(k,i)
              enddo
            enddo
          enddo
          do j = 1,ninternalmolQ
            jj = ninternalatm + ninternalmolT + j
            derv2(jj,i) = 0.0_dp
            derv2(i,jj) = 0.0_dp
            do k = 1,nmf
              indk = 3*(k-1)
              do kk = 1,3
                derv2(jj,i) = derv2(jj,i) + tmp2(indk+kk)*tmatQ(indk+kk,j)
                derv2(i,jj) = derv2(i,jj) + tmp2(indk+kk)*tmatQ(indk+kk,j)
              enddo
            enddo
          enddo
        enddo
      endif
!
      deallocate(tmp2,stat=status)
      if (status/=0) call deallocate_error('sec3','tmp2')
    else
!**************************
!  Nasym-numat algorithm  *
!**************************
!
!  Reduce matrix to symmetry reduced full second derivative form
!
      do i = 1,nasym
        ix = 3*(i-1) + 1
        iy = ix + 1
        iz = ix + 2
        do j = 1,i
          jj = nrela2f(j)
          jx = 3*(j-1) + 1
          jy = jx + 1
          jz = jx + 2
          neqj = neqv(j)
          jjx = 3*(jj-1)
          d2(1,1) = 0.0_dp
          d2(2,1) = 0.0_dp
          d2(3,1) = 0.0_dp
          d2(1,2) = 0.0_dp
          d2(2,2) = 0.0_dp
          d2(3,2) = 0.0_dp
          d2(1,3) = 0.0_dp
          d2(2,3) = 0.0_dp
          d2(3,3) = 0.0_dp
          do jk = 1,3*neqj
            d2(1,1) = d2(1,1) + derv2(jjx+jk,ix)*tmat(jjx+jk,ix)
            d2(2,1) = d2(2,1) + derv2(jjx+jk,iy)*tmat(jjx+jk,ix)
            d2(3,1) = d2(3,1) + derv2(jjx+jk,iz)*tmat(jjx+jk,ix)
            d2(1,2) = d2(1,2) + derv2(jjx+jk,ix)*tmat(jjx+jk,iy)
            d2(2,2) = d2(2,2) + derv2(jjx+jk,iy)*tmat(jjx+jk,iy)
            d2(3,2) = d2(3,2) + derv2(jjx+jk,iz)*tmat(jjx+jk,iy)
            d2(1,3) = d2(1,3) + derv2(jjx+jk,ix)*tmat(jjx+jk,iz)
            d2(2,3) = d2(2,3) + derv2(jjx+jk,iy)*tmat(jjx+jk,iz)
            d2(3,3) = d2(3,3) + derv2(jjx+jk,iz)*tmat(jjx+jk,iz)
          enddo
          derv2(jx,ix) = d2(1,1)
          derv2(jy,ix) = d2(1,2)
          derv2(jz,ix) = d2(1,3)
          derv2(jx,iy) = d2(2,1)
          derv2(jy,iy) = d2(2,2)
          derv2(jz,iy) = d2(2,3)
          derv2(jx,iz) = d2(3,1)
          derv2(jy,iz) = d2(3,2)
          derv2(jz,iz) = d2(3,3)
        enddo
!
!  Breathing shells - coordinate-radius
!
        if (nbsm.gt.0.and.lbsmat(nsft+i)) then
          ibs = nofa + i
          do j = 1,nasym
            jj = nrela2f(j)
            neqj = neqv(j)
            d2x = 0.0_dp
            d2y = 0.0_dp
            d2z = 0.0_dp
            kk = 3*(jj-1)
            do k = 1,neqj
              io = nrotop(jj+k-1)
              d2x = d2x + derv2(1+kk,ibs)*rop(1,1,io) + derv2(2+kk,ibs)*rop(2,1,io) + derv2(3+kk,ibs)*rop(3,1,io)
              d2y = d2y + derv2(1+kk,ibs)*rop(1,2,io) + derv2(2+kk,ibs)*rop(2,2,io) + derv2(3+kk,ibs)*rop(3,2,io)
              d2z = d2z + derv2(1+kk,ibs)*rop(1,3,io) + derv2(2+kk,ibs)*rop(2,3,io) + derv2(3+kk,ibs)*rop(3,3,io)
              kk = kk + 3
            enddo
            indj = 3*(j-1)
            derv2(indj+1,ibs) = d2x
            derv2(indj+2,ibs) = d2y
            derv2(indj+3,ibs) = d2z
          enddo
        endif
      enddo
!
!  Breathing shells - radius-radius
!
      if (nbsm.gt.0) then
        do i = 1,nasym
          do j = 1,numat
            jj = nrelf2a(j)
            if (lbsmat(nsft+i).and.lbsmat(nsft+jj)) then
              if (j.eq.nrela2f(jj)) then
                derv2(nofa+jj,nofa+i) = derv2(noff+j,nofa+i)
              else
                derv2(nofa+jj,nofa+i) = derv2(nofa+jj,nofa+i) + derv2(noff+j,nofa+i)
              endif
            else
              derv2(nofa+jj,nofa+i) = 0.0_dp
            endif
          enddo
        enddo
      endif
!
!  Symmetrise matrix
!
      do i = 2,n3a
        do j = 1,i-1
          derv2(i,j) = derv2(j,i)
        enddo
      enddo
!
!  Reduce to nvar x nvar form
!
      if (ncon.eq.0) then
        do i = ninternalmin,ninternalmax
          if (iopttype(i).eq.iopt_xf) then
            id = 3*ioptindex(i) - 2
          elseif (iopttype(i).eq.iopt_yf) then
            id = 3*ioptindex(i) - 1
          elseif (iopttype(i).eq.iopt_zf) then
            id = 3*ioptindex(i)
          elseif (iopttype(i).eq.iopt_radius) then
            id = 3*nasym + ioptindex(i)
          endif
          do j = ninternalmin,ninternalmax
            if (iopttype(j).eq.iopt_xf) then
              jd = 3*ioptindex(j) - 2
            elseif (iopttype(j).eq.iopt_yf) then
              jd = 3*ioptindex(j) - 1
            elseif (iopttype(j).eq.iopt_zf) then
              jd = 3*ioptindex(j)
            elseif (iopttype(j).eq.iopt_radius) then
              jd = 3*nasym + ioptindex(j)
            endif
            derv2(j,i) = derv2(jd,id)
          enddo
        enddo
      else
        allocate(tmp2(ninternalatm),stat=status)
        if (status/=0) call outofmemory('sec3','tmp2')
        do i = 1,n3a
          do j = 1,ninternalatm
            tmp2(j) = 0.0_dp
            do k = 1,n3a
              tmp2(j) = tmp2(j) + tmat(k,n3a+j)*derv2(k,i)
            enddo
          enddo
          do j = 1,ninternalatm
            derv2(j,i) = tmp2(j)
          enddo
        enddo
        do i = 1,ninternalatm
          do j = 1,ninternalatm
            tmp2(j) = 0.0_dp
            do k = 1,n3a
              tmp2(j) = tmp2(j) + tmat(k,n3a+j)*derv2(i,k)
            enddo
          enddo
          do j = 1,ninternalatm
            derv2(i,j) = tmp2(j)
          enddo
        enddo
        deallocate(tmp2,stat=status)
        if (status/=0) call deallocate_error('sec3','tmp2')
      endif
    endif
  elseif (nfixatom.gt.0.and.(ninternal.eq.(n3f-3).and.ninternal.gt.0).and..not.lrigid) then
!
!  All atoms are variable bar one & so only need to copy from fixed atom onwards
!
    nintmin = 3*(nfixatom-1) 
    do i = nintmin+1,ninternal
      do j = nintmin+1,ninternal
        derv2(j,i) = derv2(j+3,i+3)
      enddo
    enddo
    do i = nintmin+1,ninternal
      do j = 1,nintmin
        derv2(j,i) = derv2(j,i+3)
      enddo
    enddo
    do i = 1,nintmin
      do j = nintmin+1,ninternal
        derv2(j,i) = derv2(j+3,i)
      enddo
    enddo
  else
    do i = ninternalmin,ninternalmax
      lmolQi = .false.
      lmolTi = .false.
      ni = ioptindex(i)
      if (iopttype(i).eq.iopt_xf) then
        ii = 3*nasymnomolptr(ni) - 2
      elseif (iopttype(i).eq.iopt_yf) then
        ii = 3*nasymnomolptr(ni) - 1
      elseif (iopttype(i).eq.iopt_zf) then
        ii = 3*nasymnomolptr(ni)
      elseif (iopttype(i).eq.iopt_radius) then
        ii = 3*nasym + nasymnomolptr(ni)
      elseif (iopttype(i).eq.iopt_xcom) then
        lmolTi = .true.
        ii = 1
      elseif (iopttype(i).eq.iopt_ycom) then
        lmolTi = .true.
        ii = 2
      elseif (iopttype(i).eq.iopt_zcom) then
        lmolTi = .true.
        ii = 3
      elseif (iopttype(i).eq.iopt_xqtn) then
        lmolQi = .true.
        ii = 1
      elseif (iopttype(i).eq.iopt_yqtn) then
        lmolQi = .true.
        ii = 2
      elseif (iopttype(i).eq.iopt_zqtn) then
        lmolQi = .true.
        ii = 3
      endif
      do j = ninternalmin,ninternalmax
        lmolQj = .false.
        lmolTj = .false.
        nj = ioptindex(j)
        if (iopttype(j).eq.iopt_xf) then
          jj = 3*nasymnomolptr(nj) - 2
        elseif (iopttype(j).eq.iopt_yf) then
          jj = 3*nasymnomolptr(nj) - 1
        elseif (iopttype(j).eq.iopt_zf) then
          jj = 3*nasymnomolptr(nj)
        elseif (iopttype(j).eq.iopt_radius) then
          jj = 3*nasym + nasymnomolptr(nj)
        elseif (iopttype(j).eq.iopt_xcom) then
          lmolTj = .true.
          jj = 1
        elseif (iopttype(j).eq.iopt_ycom) then
          lmolTj = .true.
          jj = 2
        elseif (iopttype(j).eq.iopt_zcom) then
          lmolTj = .true.
          jj = 3
        elseif (iopttype(j).eq.iopt_xqtn) then
          lmolQj = .true.
          jj = 1
        elseif (iopttype(j).eq.iopt_yqtn) then
          lmolQj = .true.
          jj = 2
        elseif (iopttype(j).eq.iopt_zqtn) then
          lmolQj = .true.
          jj = 3
        endif
        if (lmolQi) then
          if (lmolQj) then
            derv2(j,i) = molQQdrv(jj,ii,nj,ni)
          elseif (lmolTj) then
            derv2(j,i) = molQTdrv(ii,jj,ni,nj)
          else
            derv2(j,i) = molQCdrv(jj,ii,ni)
          endif
        elseif (lmolTi) then
          if (lmolQj) then
            derv2(j,i) = molQTdrv(jj,ii,nj,ni)
          elseif (lmolTj) then
            derv2(j,i) = molTTdrv(jj,ii,nj,ni)
          else
            derv2(j,i) = molTCdrv(jj,ii,ni)
          endif
        else
          if (lmolQj) then
            derv2(j,i) = molQCdrv(ii,jj,nj)
          elseif (lmolTj) then
            derv2(j,i) = molTCdrv(ii,jj,nj)
          else
            derv2(j,i) = derv2(jj,ii)
          endif
        endif
      enddo
    enddo
  endif
!
!  Check size of derv2
!
  if (ninternal+ncell.gt.maxd2u) then
    maxd2u = ninternal + ncell
    call changemaxd2
  endif
  if (ninternal+ncell.gt.maxd2) then
    maxd2 = ninternal + ncell
    call changemaxd2
  endif
  if (loldvarorder.and.ncell.gt.0) then
!
!  Shift position of second derivatives 
!
    do i = ninternal,1,-1
      do j = ninternal,1,-1
        derv2(j+ncell,i+ncell) = derv2(j,i)
      enddo
    enddo
  endif
!**********************
!  Mixed derivatives  *
!**********************
  if (lstr) then
!
!  Transform to fractional derivatives
!
    if (lfullra) then
      do jj = 1,nstrains
        do i = 1,nloop
          indi = 3*(i-1)
          derv3(indi+1,jj) = derv3(indi+1,jj)*r1xl
          derv3(indi+2,jj) = derv3(indi+2,jj)*r2yl
          derv3(indi+3,jj) = derv3(indi+3,jj)*r3zl
        enddo
      enddo
      if (lrigid) then
        do nm = 1,nmol
!
!  Molecule - strain
!
          do jj = 1,nstrains
            molTSdrv(jj,1,nm) = molTSdrv(jj,1,nm)*r1xl
            molTSdrv(jj,2,nm) = molTSdrv(jj,2,nm)*r2yl
            molTSdrv(jj,3,nm) = molTSdrv(jj,3,nm)*r3zl
          enddo
        enddo
      endif
    else
      do i = 1,nloop
        indi = 3*(i-1)
        do ii = 1,3
          do jj = 1,nstrains
            tmp1(ii,jj) = derv3(indi+1,jj)*rmat(1,ii) + &
                          derv3(indi+2,jj)*rmat(2,ii) + &
                          derv3(indi+3,jj)*rmat(3,ii)
          enddo
        enddo
        do ii = 1,3
          do jj = 1,nstrains
            derv3(indi+ii,jj) = tmp1(ii,jj)
          enddo
        enddo
      enddo
      if (lrigid) then
        do nm = 1,nmol
!
!  Molecule - strain
!
          do ii = 1,3
            do jj = 1,nstrains
              tmp1(ii,jj) = molTSdrv(jj,1,nm)*rmat(1,ii) + &
                            molTSdrv(jj,2,nm)*rmat(2,ii) + &
                            molTSdrv(jj,3,nm)*rmat(3,ii)
            enddo
          enddo
          do ii = 1,3
            do jj = 1,nstrains
              molTSdrv(jj,ii,nm) = tmp1(ii,jj)
            enddo
          enddo
        enddo
      endif
    endif
!
    if (loptcellpar) then
!
!  Transform to cell parameter derivatives
!
      if (ndim.eq.3) then
        call setcellderv3D(.true.,.false.,.true.)
      elseif (ndim.eq.2) then
        call setcellderv2D(.true.,.false.,.true.)
      elseif (ndim.eq.1) then
        call setcellderv1D(.true.,.true.)
      endif
      do n = 1,ninternal
        do i = 1,nstrains
          d3tmp(i) = derv3(n,i)
          derv3(n,i) = 0.0_dp
        enddo
        do i = 1,nstrains
          do j = 1,nstrains
            derv3(n,j) = derv3(n,j) + d3tmp(i)*cderv(j,i)
          enddo    
        enddo
      enddo
    endif
!*****************************************
!  Symmetry transform mixed derivatives  *
!*****************************************
    if (ltmat) then
      if (.not.lsymderv2) then
        allocate(tmp2(n3f),stat=status)
        if (status/=0) call outofmemory('sec3','tmp2')
        do i = 1,nstrains
          do j = 1,n3f
            tmp2(j) = derv3(j,i)
          enddo
          do j = 1,ninternalatm
            derv3(j,i) = 0.0_dp
            do k = 1,n3f
              derv3(j,i) = derv3(j,i) + tmat(k,j)*tmp2(k)
            enddo
          enddo
        enddo
!
        if (lrigid) then
          if (ninternalmolT.gt.0) then
            do i = 1,nstrains
              do j = 1,nmf
                indj = 3*(j-1)
                do jj = 1,3
                  tmp2(indj+jj) = molTSdrv(i,jj,j)
                enddo
              enddo
              do j = 1,ninternalmolT
                jj = ninternalatm + j
                derv3(jj,i) = 0.0_dp
                do k = 1,n3mf
                  derv3(jj,i) = derv3(jj,i) + tmp2(k)*tmatT(k,j)
                enddo
              enddo
            enddo
          endif
!
          if (ninternalmolQ.gt.0) then
            do i = 1,nstrains
              do j = 1,nmf
                indj = 3*(j-1)
                do jj = 1,3
                  tmp2(indj+jj) = molQSdrv(i,jj,j)
                enddo
              enddo
              do j = 1,ninternalmolQ
                jj = ninternalatm + ninternalmolT + j
                derv3(jj,i) = 0.0_dp
                do k = 1,n3mf
                  derv3(jj,i) = derv3(jj,i) + tmp2(k)*tmatQ(k,j)
                enddo
              enddo
            enddo
          endif
        endif
!
        deallocate(tmp2,stat=status)
        if (status/=0) call deallocate_error('sec3','tmp2')
      else
        if (ncon.eq.0) then
          do i = 1,nstrains
            do j = ninternalmin,ninternalmax
              if (iopttype(j).eq.iopt_xf) then
                jd = 3*ioptindex(j) - 2
              elseif (iopttype(j).eq.iopt_yf) then
                jd = 3*ioptindex(j) - 1
              elseif (iopttype(j).eq.iopt_zf) then
                jd = 3*ioptindex(j)
              elseif (iopttype(j).eq.iopt_radius) then
                jd = 3*nasym + ioptindex(j)
              endif
              derv3(j,i) = derv3(jd,i)
            enddo
          enddo
        else
          allocate(tmp2(ninternalatm),stat=status)
          if (status/=0) call outofmemory('sec3','tmp2')
          do i = 1,nstrains
            do j = 1,ninternalatm
              tmp2(j) = 0.0_dp
              do k = 1,n3a
                tmp2(j) = tmp2(j) + tmat(k,n3a+j)*derv3(k,i)
              enddo
            enddo
            do j = 1,ninternalatm
              derv3(j,i) = tmp2(j)
            enddo
          enddo
          deallocate(tmp2,stat=status)
          if (status/=0) call deallocate_error('sec3','tmp2')
        endif
      endif
    elseif (nfixatom.gt.0.and.(ninternal.eq.(n3f-3).and.ninternal.gt.0).and..not.lrigid) then
!
!  All atoms are variable bar one & so only need to copy from fixed atom onwards
!
      nintmin = 3*(nfixatom-1)
      do i = 1,nstrains
        do j = nintmin+1,ninternal
          derv3(j,i) = derv3(j+3,i)
        enddo
      enddo
    else
      do i = 1,nstrains
        do j = ninternalmin,ninternalmax
          nj = ioptindex(j)
          if (iopttype(j).eq.iopt_xf) then
            jd = 3*nasymnomolptr(nj) - 2
            derv3(j,i) = derv3(jd,i)
          elseif (iopttype(j).eq.iopt_yf) then
            jd = 3*nasymnomolptr(nj) - 1
            derv3(j,i) = derv3(jd,i)
          elseif (iopttype(j).eq.iopt_zf) then
            jd = 3*nasymnomolptr(nj)
            derv3(j,i) = derv3(jd,i)
          elseif (iopttype(j).eq.iopt_radius) then
            jd = 3*nasym + nasymnomolptr(nj)
            derv3(j,i) = derv3(jd,i)
          elseif (iopttype(j).eq.iopt_xcom) then
            jd = 1
            derv3(j,i) = molTSdrv(i,jd,nj)
          elseif (iopttype(j).eq.iopt_ycom) then
            jd = 2
            derv3(j,i) = molTSdrv(i,jd,nj)
          elseif (iopttype(j).eq.iopt_zcom) then
            jd = 3
            derv3(j,i) = molTSdrv(i,jd,nj)
          elseif (iopttype(j).eq.iopt_xqtn) then
            jd = 1
            derv3(j,i) = molQSdrv(i,jd,nj)
          elseif (iopttype(j).eq.iopt_yqtn) then
            jd = 2
            derv3(j,i) = molQSdrv(i,jd,nj)
          elseif (iopttype(j).eq.iopt_zqtn) then
            jd = 3
            derv3(j,i) = molQSdrv(i,jd,nj)
          endif
        enddo
      enddo
    endif
!
!  Reduce derv3 by strain reduction matrix
!
    allocate(tmp2(ninternal),stat=status)
    if (status/=0) call outofmemory('sec3','tmp2')
    do i = 1,ncell
      do j = 1,ninternal
        tmp2(j) = 0.0_dp
        do k = 1,nstrains
          tmp2(j) = tmp2(j) + stmat(k,i)*derv3(j,k)
        enddo
      enddo
      do j = 1,ninternal
        derv3(j,i) = tmp2(j)
      enddo
    enddo
    deallocate(tmp2,stat=status)
    if (status/=0) call deallocate_error('sec3','tmp2')
    if (loldvarorder) then
!
!  Move into place
!
      do i = 1,ncell
        do j = 1,ninternal
          derv2(ncell+j,i) = derv3(j,i)
          derv2(i,j+ncell) = derv3(j,i)
        enddo
      enddo
    else
      do i = 1,ncell
        ii = ncellmin - 1 + i
        do j = 1,ninternal
          derv2(j,ii) = derv3(j,i)
          derv2(ii,j) = derv3(j,i)
        enddo
      enddo
    endif
!******************************
!  Strain second derivatives  *
!******************************
    if (loptcellpar) then
      if (ndim.eq.3) then
        ncellfull = 6
      elseif (ndim.eq.2) then
        ncellfull = 3
      else
        ncellfull = 1
      endif
!
!  Convert strain second derivatives to cell parameter second derivatives
!
     celldrv2(1:ncellfull,1:ncellfull) = 0.0_dp
      do i = 1,ncellfull
        do j = 1,ncellfull
          do k = 1,nstrains
            do l = 1,nstrains
              celldrv2(j,i) = celldrv2(j,i) + sdrv2(l,k)*cderv(j,l)*cderv(i,k)
            enddo
          enddo
        enddo
      enddo
      sdrv2(1:ncellfull,1:ncellfull) = celldrv2(1:ncellfull,1:ncellfull)
    endif
    if (ncon.eq.0) then
      do i = ncellmin,ncellmax
        io = ioptindex(i)
        do j = ncellmin,ncellmax
          jo = ioptindex(j)
          dum = sdrv2(jo,io)
          derv2(j,i) = dum
          derv2(i,j) = dum
        enddo
      enddo
    else
      do i = 1,ncell
        do j = 1,nstrains
          tmp1(j,i) = 0.0_dp
          do k = 1,nstrains
            tmp1(j,i) = tmp1(j,i) + sdrv2(j,k)*stmat(k,i)
          enddo
        enddo
      enddo
      do i = 1,ncell
        ii = ncellmin - 1 + i
        do j = 1,ncell
          jj = ncellmin - 1 + j
          derv2(jj,ii) = 0.0_dp
          do k = 1,nstrains
            derv2(jj,ii) = derv2(jj,ii) + tmp1(k,j)*stmat(k,i)
          enddo
        enddo
      enddo
    endif
  endif
  if (lsdebug.and.ioproc) then
    write(ioout,'(/,''  Symmetrised Second Derivative Matrix  :'',/)')
    do i = 1,nvar
      write(ioout,'(2x,9(f9.4))')(derv2(j,i),j=1,i)
    enddo
    write(ioout,'(/)')
  endif
!
  t2 = g_cpu_time()
  thes = thes + t2 - t1
#ifdef TRACE
  call trace_out('sec3')
#endif
!
  return
  end
