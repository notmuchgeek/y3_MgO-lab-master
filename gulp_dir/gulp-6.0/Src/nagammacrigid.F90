  subroutine nagammacrigid(mtvrptr,rfmass,qfrac,maxd2,derv2,dervi,xk,yk,zk)
!
!  Calculates the non-analytic correction to the second
!  derivatives based on the Born effective charges.
!  Rigid molecule/complex version.
!
!  On entry :
!
!  mtvrptr  = pointer from dynamical matrix row global index to local index
!  rfmass   = inverse square root of atomic masses
!  qfrac    = fractional approach direction to gamma point
!  maxd2    = left-hand dimension of derv2
!  derv2      = uncorrected dynamical matrix (real)
!  dervi      = uncorrected dynamical matrix (complex)
!  xk         = x component of k point in reciprocal space
!  yk         = y component of k point in reciprocal space
!  zk         = z component of k point in reciprocal space
!
!  On exit :
!
!  derv2  = corrected dynamical matrix (real)
!  dervi  = corrected dynamical matrix (complex)
!
!   5/20 Created from nagammac and nagammarigid
!   6/20 Parallel modifications made
!   6/20 nmolcore changes added
!   6/20 Changes made for shell model
!   6/20 Correction to indices for atom - atom/molecule terms
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
  use g_constants
  use configurations, only : nsuperghost
  use control,        only : lkfull
  use current
  use element,        only : maxele
  use molecule
  use parallel
  use phononatoms
  use properties
  use shells
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4),                    intent(in)    :: mtvrptr(*)
  integer(i4),                    intent(in)    :: maxd2
  real(dp),                       intent(in)    :: qfrac(3)
  real(dp),                       intent(in)    :: rfmass(*)
  real(dp),                       intent(inout) :: derv2(maxd2,*)
  real(dp),                       intent(inout) :: dervi(maxd2,*)
  real(dp),                       intent(in)    :: xk
  real(dp),                       intent(in)    :: yk
  real(dp),                       intent(in)    :: zk
!
!  Local variables
!
  integer(i4)                                   :: i
  integer(i4)                                   :: ia
  integer(i4)                                   :: ii
  integer(i4)                                   :: iloc
  integer(i4)                                   :: indi
  integer(i4)                                   :: indm
  integer(i4)                                   :: ix
  integer(i4)                                   :: iy
  integer(i4)                                   :: iz
  integer(i4)                                   :: ixm
  integer(i4)                                   :: iym
  integer(i4)                                   :: izm
  integer(i4)                                   :: j
  integer(i4)                                   :: j1
  integer(i4)                                   :: jj
  integer(i4)                                   :: jndm
  integer(i4)                                   :: jx
  integer(i4)                                   :: jy
  integer(i4)                                   :: jz
  integer(i4)                                   :: jxm
  integer(i4)                                   :: jym
  integer(i4)                                   :: jzm
  integer(i4)                                   :: k
  integer(i4)                                   :: k1
  integer(i4)                                   :: l
  integer(i4)                                   :: m
  integer(i4)                                   :: n
  integer(i4)                                   :: nm
  integer(i4)                                   :: nn
  integer(i4)                                   :: status
  real(dp)                                      :: cosk
  real(dp)                                      :: dQ(3)
  real(dp)                                      :: drR(3,3,3)
  real(dp)                                      :: drS(3,3,3)
  real(dp)                                      :: drRi(3,3)
  real(dp)                                      :: drRj(3,3)
  real(dp)                                      :: d2(3,3)
  real(dp),   dimension(:,:), allocatable       :: qz
  real(dp)                                      :: factor
  real(dp)                                      :: inveps
  real(dp)                                      :: kvf(3,3)
  real(dp)                                      :: qcart(3)
  real(dp)                                      :: sink
  real(dp)                                      :: vol
  real(dp)                                      :: volume
  real(dp)                                      :: xkv
  real(dp)                                      :: ykv
  real(dp)                                      :: zkv
#ifdef TRACE
  call trace_in('nagammacrigid')
#endif
!
!  Calculate constant = 4*pi/V
!
  vol = volume(rv)
  factor = 4.0_dp*pi/vol
!
!  Select appropriate K vectors
!
  if (lkfull) then
    call kvector3Df(kvf)
  else
    kvf(1:3,1:3) = kv(1:3,1:3)
  endif
!
!  Correct K vectors for ghost supercell
!
  kvf(1:3,1) = kvf(1:3,1)*dble(nsuperghost(1,ncf))
  kvf(1:3,2) = kvf(1:3,2)*dble(nsuperghost(2,ncf))
  kvf(1:3,3) = kvf(1:3,3)*dble(nsuperghost(3,ncf))
! 
!  Calculate Cartesian space direction
! 
  qcart(1) = qfrac(1)*kvf(1,1) + qfrac(2)*kvf(1,2) + qfrac(3)*kvf(1,3)
  qcart(2) = qfrac(1)*kvf(2,1) + qfrac(2)*kvf(2,2) + qfrac(3)*kvf(2,3)
  qcart(3) = qfrac(1)*kvf(3,1) + qfrac(2)*kvf(3,2) + qfrac(3)*kvf(3,3)
!
!  Build arrays of Q.Z
!
  allocate(qz(3,ncore),stat=status)
  if (status/=0) call outofmemory('nagammacrigid','qz')
  qz(1:3,1:ncore) = 0.0_dp
  do i = 1,ncore
    qz(1,i) = qz(1,i) + qcart(1)*bornq(1,1,i) + qcart(2)*bornq(2,1,i) + qcart(3)*bornq(3,1,i)
    qz(2,i) = qz(2,i) + qcart(1)*bornq(1,2,i) + qcart(2)*bornq(2,2,i) + qcart(3)*bornq(3,2,i)
    qz(3,i) = qz(3,i) + qcart(1)*bornq(1,3,i) + qcart(2)*bornq(2,3,i) + qcart(3)*bornq(3,3,i)
  enddo
!
!  Calculate Q.epsilon(static).Q and invert
!
  inveps = qcart(1)*diconh(1,1)*qcart(1) + qcart(2)*diconh(2,1)*qcart(1) + qcart(3)*diconh(3,1)*qcart(1) + &
           qcart(1)*diconh(1,2)*qcart(2) + qcart(2)*diconh(2,2)*qcart(2) + qcart(3)*diconh(3,2)*qcart(2) + &
           qcart(1)*diconh(1,3)*qcart(3) + qcart(2)*diconh(2,3)*qcart(3) + qcart(3)*diconh(3,3)*qcart(3)
  inveps = 1.0_dp/inveps
!
!  Loop over second derivative matrix elements adding terms
!
  factor = factor*inveps*angstoev
!
!  Compute point in BZ in Cartesian form
!
  xkv = xk*kvf(1,1) + yk*kvf(1,2) + zk*kvf(1,3)
  ykv = xk*kvf(2,1) + yk*kvf(2,2) + zk*kvf(2,3)
  zkv = xk*kvf(3,1) + yk*kvf(3,2) + zk*kvf(3,3)
!
!  Atom contribution excluding molecules
!
  ix = - 2
  iy = - 1
  iz =   0
  do iloc = 1,numatnomolonnode
    ii = numatnomolonnodeptr(iloc)
    if (nat(ii).le.maxele) then
      indm = 3*(iloc-1)
      ix = ix + 3
      iy = iy + 3
      iz = iz + 3
      jx = - 2
      jy = - 1
      jz =   0
      do j = 1,ncorenomol
        jj = ncorenomolptr(j)
        jx = jx + 3
        jy = jy + 3
        jz = jz + 3
!
!  Compute phase factor for this K point
!
        cosk = xkv*(xclat(jj)-xclat(ii)) + ykv*(yclat(jj)-yclat(ii)) + zkv*(zclat(jj)-zclat(ii))
        sink = sin(cosk)
        cosk = cos(cosk)
!
        derv2(jx,ix) = derv2(jx,ix) + factor*rfmass(indm+1)*rfmass(jx)*qz(1,jj)*qz(1,ii)*cosk
        derv2(jy,ix) = derv2(jy,ix) + factor*rfmass(indm+1)*rfmass(jy)*qz(2,jj)*qz(1,ii)*cosk
        derv2(jz,ix) = derv2(jz,ix) + factor*rfmass(indm+1)*rfmass(jz)*qz(3,jj)*qz(1,ii)*cosk
        derv2(jx,iy) = derv2(jx,iy) + factor*rfmass(indm+2)*rfmass(jx)*qz(1,jj)*qz(2,ii)*cosk
        derv2(jy,iy) = derv2(jy,iy) + factor*rfmass(indm+2)*rfmass(jy)*qz(2,jj)*qz(2,ii)*cosk
        derv2(jz,iy) = derv2(jz,iy) + factor*rfmass(indm+2)*rfmass(jz)*qz(3,jj)*qz(2,ii)*cosk
        derv2(jx,iz) = derv2(jx,iz) + factor*rfmass(indm+3)*rfmass(jx)*qz(1,jj)*qz(3,ii)*cosk
        derv2(jy,iz) = derv2(jy,iz) + factor*rfmass(indm+3)*rfmass(jy)*qz(2,jj)*qz(3,ii)*cosk
        derv2(jz,iz) = derv2(jz,iz) + factor*rfmass(indm+3)*rfmass(jz)*qz(3,jj)*qz(3,ii)*cosk
!
        dervi(jx,ix) = dervi(jx,ix) + factor*rfmass(indm+1)*rfmass(jx)*qz(1,jj)*qz(1,ii)*sink
        dervi(jy,ix) = dervi(jy,ix) + factor*rfmass(indm+1)*rfmass(jy)*qz(2,jj)*qz(1,ii)*sink
        dervi(jz,ix) = dervi(jz,ix) + factor*rfmass(indm+1)*rfmass(jz)*qz(3,jj)*qz(1,ii)*sink
        dervi(jx,iy) = dervi(jx,iy) + factor*rfmass(indm+2)*rfmass(jx)*qz(1,jj)*qz(2,ii)*sink
        dervi(jy,iy) = dervi(jy,iy) + factor*rfmass(indm+2)*rfmass(jy)*qz(2,jj)*qz(2,ii)*sink
        dervi(jz,iy) = dervi(jz,iy) + factor*rfmass(indm+2)*rfmass(jz)*qz(3,jj)*qz(2,ii)*sink
        dervi(jx,iz) = dervi(jx,iz) + factor*rfmass(indm+3)*rfmass(jx)*qz(1,jj)*qz(3,ii)*sink
        dervi(jy,iz) = dervi(jy,iz) + factor*rfmass(indm+3)*rfmass(jy)*qz(2,jj)*qz(3,ii)*sink
        dervi(jz,iz) = dervi(jz,iz) + factor*rfmass(indm+3)*rfmass(jz)*qz(3,jj)*qz(3,ii)*sink
      enddo
!
!  Molecule - non-molecule : translation
!
      do nn = 1,nmol
        jx = 3*ncorenomol + 6*(nn - 1) + 1
        jy = jx + 1
        jz = jx + 2
        do n = 1,nmolcore(nn)
          j = nmollist(nmolptr(nn)+n)
!
!  Compute phase factor for this K point
!
          if (lphasecom) then
            cosk = xkv*(xclat(j)-xclat(ii)-molxyz(1,n,nn)) + &
                   ykv*(yclat(j)-yclat(ii)-molxyz(2,n,nn)) + &
                   zkv*(zclat(j)-zclat(ii)-molxyz(3,n,nn))
          else
            cosk = xkv*(xclat(j)-xclat(ii)) + ykv*(yclat(j)-yclat(ii)) + zkv*(zclat(j)-zclat(ii))
          endif
          sink = sin(cosk)
          cosk = cos(cosk)
!
!  molTCdrv
!
          derv2(jx,ix) = derv2(jx,ix) + factor*rfmass(indm+1)*rfmass(jx)*qz(1,j)*qz(1,ii)*cosk
          derv2(jy,ix) = derv2(jy,ix) + factor*rfmass(indm+1)*rfmass(jy)*qz(2,j)*qz(1,ii)*cosk
          derv2(jz,ix) = derv2(jz,ix) + factor*rfmass(indm+1)*rfmass(jz)*qz(3,j)*qz(1,ii)*cosk
          derv2(jx,iy) = derv2(jx,iy) + factor*rfmass(indm+2)*rfmass(jx)*qz(1,j)*qz(2,ii)*cosk
          derv2(jy,iy) = derv2(jy,iy) + factor*rfmass(indm+2)*rfmass(jy)*qz(2,j)*qz(2,ii)*cosk
          derv2(jz,iy) = derv2(jz,iy) + factor*rfmass(indm+2)*rfmass(jz)*qz(3,j)*qz(2,ii)*cosk
          derv2(jx,iz) = derv2(jx,iz) + factor*rfmass(indm+3)*rfmass(jx)*qz(1,j)*qz(3,ii)*cosk
          derv2(jy,iz) = derv2(jy,iz) + factor*rfmass(indm+3)*rfmass(jy)*qz(2,j)*qz(3,ii)*cosk
          derv2(jz,iz) = derv2(jz,iz) + factor*rfmass(indm+3)*rfmass(jz)*qz(3,j)*qz(3,ii)*cosk
!
          dervi(jx,ix) = dervi(jx,ix) + factor*rfmass(indm+1)*rfmass(jx)*qz(1,j)*qz(1,ii)*sink
          dervi(jy,ix) = dervi(jy,ix) + factor*rfmass(indm+1)*rfmass(jy)*qz(2,j)*qz(1,ii)*sink
          dervi(jz,ix) = dervi(jz,ix) + factor*rfmass(indm+1)*rfmass(jz)*qz(3,j)*qz(1,ii)*sink
          dervi(jx,iy) = dervi(jx,iy) + factor*rfmass(indm+2)*rfmass(jx)*qz(1,j)*qz(2,ii)*sink
          dervi(jy,iy) = dervi(jy,iy) + factor*rfmass(indm+2)*rfmass(jy)*qz(2,j)*qz(2,ii)*sink
          dervi(jz,iy) = dervi(jz,iy) + factor*rfmass(indm+2)*rfmass(jz)*qz(3,j)*qz(2,ii)*sink
          dervi(jx,iz) = dervi(jx,iz) + factor*rfmass(indm+3)*rfmass(jx)*qz(1,j)*qz(3,ii)*sink
          dervi(jy,iz) = dervi(jy,iz) + factor*rfmass(indm+3)*rfmass(jy)*qz(2,j)*qz(3,ii)*sink
          dervi(jz,iz) = dervi(jz,iz) + factor*rfmass(indm+3)*rfmass(jz)*qz(3,j)*qz(3,ii)*sink
        enddo
      enddo
!
!  Molecule - non-molecule : rotation
!
      do nn = 1,nmol
        jndm = 3*ncorenomol + 6*(nn - 1) + 3
        jx = jndm + 1
        jy = jndm + 2
        jz = jndm + 3
!
!  Set up rotation matrix derivatives
!
        call setrotation(molaxes(1,1,nn),drS)
!
        do n = 1,nmolcore(nn)
          j = nmollist(nmolptr(nn)+n)
!
!  Multiply rotation derivatives about axes by vector to atom
!
          drRj(1:3,1:3) = 0.0_dp
          do l = 1,3
            do k = 1,3
              drRj(1,l) = drRj(1,l) + drS(1,k,l)*molxyz(k,n,nn)
              drRj(2,l) = drRj(2,l) + drS(2,k,l)*molxyz(k,n,nn)
              drRj(3,l) = drRj(3,l) + drS(3,k,l)*molxyz(k,n,nn)
            enddo
          enddo
!
          d2(1,1) = factor*qz(1,j)*qz(1,ii)
          d2(2,1) = factor*qz(2,j)*qz(1,ii)
          d2(3,1) = factor*qz(3,j)*qz(1,ii)
          d2(1,2) = factor*qz(1,j)*qz(2,ii)
          d2(2,2) = factor*qz(2,j)*qz(2,ii)
          d2(3,2) = factor*qz(3,j)*qz(2,ii)
          d2(1,3) = factor*qz(1,j)*qz(3,ii)
          d2(2,3) = factor*qz(2,j)*qz(3,ii)
          d2(3,3) = factor*qz(3,j)*qz(3,ii)
!
!  Compute phase factor for this K point
!
          if (lphasecom) then
            cosk = xkv*(xclat(j)-xclat(ii)-molxyz(1,n,nn)) + &
                   ykv*(yclat(j)-yclat(ii)-molxyz(2,n,nn)) + &
                   zkv*(zclat(j)-zclat(ii)-molxyz(3,n,nn))
          else
            cosk = xkv*(xclat(j)-xclat(ii)) + ykv*(yclat(j)-yclat(ii)) + zkv*(zclat(j)-zclat(ii))
          endif
          sink = sin(cosk)
          cosk = cos(cosk)
!
!  molQCdrv
!
          do k = 1,3
            derv2(jndm+k,ix) = derv2(jndm+k,ix) + d2(1,1)*drRj(1,k)*rfmass(jndm+k)*rfmass(indm+1)*cosk
            derv2(jndm+k,iy) = derv2(jndm+k,iy) + d2(2,1)*drRj(1,k)*rfmass(jndm+k)*rfmass(indm+2)*cosk
            derv2(jndm+k,iz) = derv2(jndm+k,iz) + d2(3,1)*drRj(1,k)*rfmass(jndm+k)*rfmass(indm+3)*cosk
            derv2(jndm+k,ix) = derv2(jndm+k,ix) + d2(1,2)*drRj(2,k)*rfmass(jndm+k)*rfmass(indm+1)*cosk
            derv2(jndm+k,iy) = derv2(jndm+k,iy) + d2(2,2)*drRj(2,k)*rfmass(jndm+k)*rfmass(indm+2)*cosk
            derv2(jndm+k,iz) = derv2(jndm+k,iz) + d2(3,2)*drRj(2,k)*rfmass(jndm+k)*rfmass(indm+3)*cosk
            derv2(jndm+k,ix) = derv2(jndm+k,ix) + d2(1,3)*drRj(3,k)*rfmass(jndm+k)*rfmass(indm+1)*cosk
            derv2(jndm+k,iy) = derv2(jndm+k,iy) + d2(2,3)*drRj(3,k)*rfmass(jndm+k)*rfmass(indm+2)*cosk
            derv2(jndm+k,iz) = derv2(jndm+k,iz) + d2(3,3)*drRj(3,k)*rfmass(jndm+k)*rfmass(indm+3)*cosk
!
            dervi(jndm+k,ix) = dervi(jndm+k,ix) + d2(1,1)*drRj(1,k)*rfmass(jndm+k)*rfmass(indm+1)*sink
            dervi(jndm+k,iy) = dervi(jndm+k,iy) + d2(2,1)*drRj(1,k)*rfmass(jndm+k)*rfmass(indm+2)*sink
            dervi(jndm+k,iz) = dervi(jndm+k,iz) + d2(3,1)*drRj(1,k)*rfmass(jndm+k)*rfmass(indm+3)*sink
            dervi(jndm+k,ix) = dervi(jndm+k,ix) + d2(1,2)*drRj(2,k)*rfmass(jndm+k)*rfmass(indm+1)*sink
            dervi(jndm+k,iy) = dervi(jndm+k,iy) + d2(2,2)*drRj(2,k)*rfmass(jndm+k)*rfmass(indm+2)*sink
            dervi(jndm+k,iz) = dervi(jndm+k,iz) + d2(3,2)*drRj(2,k)*rfmass(jndm+k)*rfmass(indm+3)*sink
            dervi(jndm+k,ix) = dervi(jndm+k,ix) + d2(1,3)*drRj(3,k)*rfmass(jndm+k)*rfmass(indm+1)*sink
            dervi(jndm+k,iy) = dervi(jndm+k,iy) + d2(2,3)*drRj(3,k)*rfmass(jndm+k)*rfmass(indm+2)*sink
            dervi(jndm+k,iz) = dervi(jndm+k,iz) + d2(3,3)*drRj(3,k)*rfmass(jndm+k)*rfmass(indm+3)*sink
          enddo
        enddo
      enddo
    endif
  enddo
!----------------------------------------------
!  Molecule translational second derivatives  |
!----------------------------------------------
  do nm = 1,nmol
    iloc = ncorenomol + 2*(nm-1) + 1
    indi = 3*(iloc-1)
    if (mtvrptr(indi+1).ne.0) then
      indm = mtvrptr(indi+1) - 1
      ixm = indm + 1
      iym = indm + 2
      izm = indm + 3
      do m = 1,nmolcore(nm)
        i = nmollist(nmolptr(nm)+m)
!
!  Molecule - molecule : translation
!
        do nn = 1,nmol
          jxm = 3*ncorenomol + 6*(nn - 1) + 1
          jym = jxm + 1
          jzm = jxm + 2
          do n = 1,nmolcore(nn)
            j = nmollist(nmolptr(nn)+n)
!
!  Compute phase factor for this K point
!
            if (lphasecom) then
              cosk = xkv*(xclat(j)-xclat(i)-molxyz(1,n,nn)+molxyz(1,m,nm)) + &
                     ykv*(yclat(j)-yclat(i)-molxyz(2,n,nn)+molxyz(2,m,nm)) + &
                     zkv*(zclat(j)-zclat(i)-molxyz(3,n,nn)+molxyz(3,m,nm))
            else
              cosk = xkv*(xclat(j)-xclat(i)) + ykv*(yclat(j)-yclat(i)) + zkv*(zclat(j)-zclat(i))
            endif
            sink = sin(cosk)
            cosk = cos(cosk)
!
!  molTTdrv
!
            derv2(jxm,ixm) = derv2(jxm,ixm) + factor*rfmass(indi+1)*rfmass(jxm)*qz(1,j)*qz(1,i)*cosk
            derv2(jym,ixm) = derv2(jym,ixm) + factor*rfmass(indi+1)*rfmass(jym)*qz(2,j)*qz(1,i)*cosk
            derv2(jzm,ixm) = derv2(jzm,ixm) + factor*rfmass(indi+1)*rfmass(jzm)*qz(3,j)*qz(1,i)*cosk
            derv2(jxm,iym) = derv2(jxm,iym) + factor*rfmass(indi+2)*rfmass(jxm)*qz(1,j)*qz(2,i)*cosk
            derv2(jym,iym) = derv2(jym,iym) + factor*rfmass(indi+2)*rfmass(jym)*qz(2,j)*qz(2,i)*cosk
            derv2(jzm,iym) = derv2(jzm,iym) + factor*rfmass(indi+2)*rfmass(jzm)*qz(3,j)*qz(2,i)*cosk
            derv2(jxm,izm) = derv2(jxm,izm) + factor*rfmass(indi+3)*rfmass(jxm)*qz(1,j)*qz(3,i)*cosk
            derv2(jym,izm) = derv2(jym,izm) + factor*rfmass(indi+3)*rfmass(jym)*qz(2,j)*qz(3,i)*cosk
            derv2(jzm,izm) = derv2(jzm,izm) + factor*rfmass(indi+3)*rfmass(jzm)*qz(3,j)*qz(3,i)*cosk
!
            dervi(jxm,ixm) = dervi(jxm,ixm) + factor*rfmass(indi+1)*rfmass(jxm)*qz(1,j)*qz(1,i)*sink
            dervi(jym,ixm) = dervi(jym,ixm) + factor*rfmass(indi+1)*rfmass(jym)*qz(2,j)*qz(1,i)*sink
            dervi(jzm,ixm) = dervi(jzm,ixm) + factor*rfmass(indi+1)*rfmass(jzm)*qz(3,j)*qz(1,i)*sink
            dervi(jxm,iym) = dervi(jxm,iym) + factor*rfmass(indi+2)*rfmass(jxm)*qz(1,j)*qz(2,i)*sink
            dervi(jym,iym) = dervi(jym,iym) + factor*rfmass(indi+2)*rfmass(jym)*qz(2,j)*qz(2,i)*sink
            dervi(jzm,iym) = dervi(jzm,iym) + factor*rfmass(indi+2)*rfmass(jzm)*qz(3,j)*qz(2,i)*sink
            dervi(jxm,izm) = dervi(jxm,izm) + factor*rfmass(indi+3)*rfmass(jxm)*qz(1,j)*qz(3,i)*sink
            dervi(jym,izm) = dervi(jym,izm) + factor*rfmass(indi+3)*rfmass(jym)*qz(2,j)*qz(3,i)*sink
            dervi(jzm,izm) = dervi(jzm,izm) + factor*rfmass(indi+3)*rfmass(jzm)*qz(3,j)*qz(3,i)*sink
          enddo
        enddo
!
!  Molecule - molecule : translation - rotation
!
        do nn = 1,nmol
          jndm = 3*ncorenomol + 6*(nn - 1) + 3
!
!  Set up rotation matrix derivatives
!
          call setrotation(molaxes(1,1,nn),drS)
!
          do n = 1,nmolcore(nn)
            j = nmollist(nmolptr(nn)+n)
!
!  Multiply rotation derivatives about axes by vector to atom
!
            drRj(1:3,1:3) = 0.0_dp
            do l = 1,3
              do k = 1,3
                drRj(1,l) = drRj(1,l) + drS(1,k,l)*molxyz(k,n,nn)
                drRj(2,l) = drRj(2,l) + drS(2,k,l)*molxyz(k,n,nn)
                drRj(3,l) = drRj(3,l) + drS(3,k,l)*molxyz(k,n,nn)
              enddo
            enddo
!
            d2(1,1) = factor*qz(1,j)*qz(1,i)
            d2(2,1) = factor*qz(2,j)*qz(1,i)
            d2(3,1) = factor*qz(3,j)*qz(1,i)
            d2(1,2) = factor*qz(1,j)*qz(2,i)
            d2(2,2) = factor*qz(2,j)*qz(2,i)
            d2(3,2) = factor*qz(3,j)*qz(2,i)
            d2(1,3) = factor*qz(1,j)*qz(3,i)
            d2(2,3) = factor*qz(2,j)*qz(3,i)
            d2(3,3) = factor*qz(3,j)*qz(3,i)
!
!  Compute phase factor for this K point
!
            if (lphasecom) then
              cosk = xkv*(xclat(j)-xclat(i)-molxyz(1,n,nn)+molxyz(1,m,nm)) + &
                     ykv*(yclat(j)-yclat(i)-molxyz(2,n,nn)+molxyz(2,m,nm)) + &
                     zkv*(zclat(j)-zclat(i)-molxyz(3,n,nn)+molxyz(3,m,nm))
            else
              cosk = xkv*(xclat(j)-xclat(i)) + ykv*(yclat(j)-yclat(i)) + zkv*(zclat(j)-zclat(i))
            endif
            sink = sin(cosk)
            cosk = cos(cosk)
!
!  molQTdrv
!
            do ia = 1,3
              dQ(1) = d2(1,ia)
              dQ(2) = d2(2,ia)
              dQ(3) = d2(3,ia)
              do k = 1,3
                derv2(jndm+k,indm+ia) = derv2(jndm+k,indm+ia) + (drRj(1,k)*dQ(1) &
                                                              + drRj(2,k)*dQ(2) &
                                                              + drRj(3,k)*dQ(3)) &
                                                              *rfmass(indi+ia)*rfmass(jndm+k)*cosk
 
                dervi(jndm+k,indm+ia) = dervi(jndm+k,indm+ia) + (drRj(1,k)*dQ(1) &
                                                              + drRj(2,k)*dQ(2) &
                                                              + drRj(3,k)*dQ(3)) &
                                                              *rfmass(indi+ia)*rfmass(jndm+k)*sink
              enddo
            enddo
          enddo
        enddo
!
!  Molecule - non-molecule
!
        do nn = 1,ncorenomol
          j = ncorenomolptr(nn)
          jx = 3*(nn - 1) + 1
          jy = 3*(nn - 1) + 2
          jz = 3*(nn - 1) + 3
!
!  Compute phase factor for this K point
!
          if (lphasecom) then
            cosk = xkv*(xclat(j)-xclat(i)+molxyz(1,m,nm)) + &
                   ykv*(yclat(j)-yclat(i)+molxyz(2,m,nm)) + &
                   zkv*(zclat(j)-zclat(i)+molxyz(3,m,nm))
          else
            cosk = xkv*(xclat(j)-xclat(i)) + ykv*(yclat(j)-yclat(i)) + zkv*(zclat(j)-zclat(i))
          endif
          sink = sin(cosk)
          cosk = cos(cosk)
!
!  molTCdrv
!
          derv2(jx,ixm) = derv2(jx,ixm) + factor*rfmass(indi+1)*rfmass(jx)*qz(1,j)*qz(1,i)*cosk
          derv2(jy,ixm) = derv2(jy,ixm) + factor*rfmass(indi+1)*rfmass(jy)*qz(2,j)*qz(1,i)*cosk
          derv2(jz,ixm) = derv2(jz,ixm) + factor*rfmass(indi+1)*rfmass(jz)*qz(3,j)*qz(1,i)*cosk
          derv2(jx,iym) = derv2(jx,iym) + factor*rfmass(indi+2)*rfmass(jx)*qz(1,j)*qz(2,i)*cosk
          derv2(jy,iym) = derv2(jy,iym) + factor*rfmass(indi+2)*rfmass(jy)*qz(2,j)*qz(2,i)*cosk
          derv2(jz,iym) = derv2(jz,iym) + factor*rfmass(indi+2)*rfmass(jz)*qz(3,j)*qz(2,i)*cosk
          derv2(jx,izm) = derv2(jx,izm) + factor*rfmass(indi+3)*rfmass(jx)*qz(1,j)*qz(3,i)*cosk
          derv2(jy,izm) = derv2(jy,izm) + factor*rfmass(indi+3)*rfmass(jy)*qz(2,j)*qz(3,i)*cosk
          derv2(jz,izm) = derv2(jz,izm) + factor*rfmass(indi+3)*rfmass(jz)*qz(3,j)*qz(3,i)*cosk
!
          dervi(jx,ixm) = dervi(jx,ixm) + factor*rfmass(indi+1)*rfmass(jx)*qz(1,j)*qz(1,i)*sink
          dervi(jy,ixm) = dervi(jy,ixm) + factor*rfmass(indi+1)*rfmass(jy)*qz(2,j)*qz(1,i)*sink
          dervi(jz,ixm) = dervi(jz,ixm) + factor*rfmass(indi+1)*rfmass(jz)*qz(3,j)*qz(1,i)*sink
          dervi(jx,iym) = dervi(jx,iym) + factor*rfmass(indi+2)*rfmass(jx)*qz(1,j)*qz(2,i)*sink
          dervi(jy,iym) = dervi(jy,iym) + factor*rfmass(indi+2)*rfmass(jy)*qz(2,j)*qz(2,i)*sink
          dervi(jz,iym) = dervi(jz,iym) + factor*rfmass(indi+2)*rfmass(jz)*qz(3,j)*qz(2,i)*sink
          dervi(jx,izm) = dervi(jx,izm) + factor*rfmass(indi+3)*rfmass(jx)*qz(1,j)*qz(3,i)*sink
          dervi(jy,izm) = dervi(jy,izm) + factor*rfmass(indi+3)*rfmass(jy)*qz(2,j)*qz(3,i)*sink
          dervi(jz,izm) = dervi(jz,izm) + factor*rfmass(indi+3)*rfmass(jz)*qz(3,j)*qz(3,i)*sink
        enddo
      enddo
    endif
  enddo
!-------------------------------------------
!  Molecule rotational second derivatives  |
!-------------------------------------------
  do nm = 1,nmol
    iloc = ncorenomol + 2*(nm-1) + 2
    indi = 3*(iloc-1)
    if (mtvrptr(indi+1).ne.0) then
      indm = mtvrptr(indi+1) - 1
!
!  Set up rotation matrix derivatives
!
      call setrotation(molaxes(1,1,nm),drR)
      do m = 1,nmolcore(nm)
        i = nmollist(nmolptr(nm)+m)
!
!  Multiply rotation derivatives about axes by vector to atom
!
        drRi(1:3,1:3) = 0.0_dp
        do l = 1,3
          do k = 1,3
            drRi(1,l) = drRi(1,l) + drR(1,k,l)*molxyz(k,m,nm)
            drRi(2,l) = drRi(2,l) + drR(2,k,l)*molxyz(k,m,nm)
            drRi(3,l) = drRi(3,l) + drR(3,k,l)*molxyz(k,m,nm)
          enddo
        enddo
!
!  Loop over non-molecule atoms
!
        do nn = 1,ncorenomol
          j = ncorenomolptr(nn)
          jx = 3*(nn - 1) + 1
          jy = 3*(nn - 1) + 2
          jz = 3*(nn - 1) + 3
!
          d2(1,1) = factor*qz(1,j)*qz(1,i)
          d2(2,1) = factor*qz(2,j)*qz(1,i)
          d2(3,1) = factor*qz(3,j)*qz(1,i)
          d2(1,2) = factor*qz(1,j)*qz(2,i)
          d2(2,2) = factor*qz(2,j)*qz(2,i)
          d2(3,2) = factor*qz(3,j)*qz(2,i)
          d2(1,3) = factor*qz(1,j)*qz(3,i)
          d2(2,3) = factor*qz(2,j)*qz(3,i)
          d2(3,3) = factor*qz(3,j)*qz(3,i)
!
!  Compute phase factor for this K point
!
          if (lphasecom) then
            cosk = xkv*(xclat(j)-xclat(i)+molxyz(1,m,nm)) + &
                   ykv*(yclat(j)-yclat(i)+molxyz(2,m,nm)) + &
                   zkv*(zclat(j)-zclat(i)+molxyz(3,m,nm))
          else
            cosk = xkv*(xclat(j)-xclat(i)) + ykv*(yclat(j)-yclat(i)) + zkv*(zclat(j)-zclat(i))
          endif
          sink = sin(cosk)
          cosk = cos(cosk)
!
          do k = 1,3
!
!  molQCdrv
!
            derv2(jx,indm+k) = derv2(jx,indm+k) + d2(1,1)*drRi(1,k)*rfmass(jx)*rfmass(indi+k)*cosk
            derv2(jy,indm+k) = derv2(jy,indm+k) + d2(2,1)*drRi(1,k)*rfmass(jy)*rfmass(indi+k)*cosk
            derv2(jz,indm+k) = derv2(jz,indm+k) + d2(3,1)*drRi(1,k)*rfmass(jz)*rfmass(indi+k)*cosk
            derv2(jx,indm+k) = derv2(jx,indm+k) + d2(1,2)*drRi(2,k)*rfmass(jx)*rfmass(indi+k)*cosk
            derv2(jy,indm+k) = derv2(jy,indm+k) + d2(2,2)*drRi(2,k)*rfmass(jy)*rfmass(indi+k)*cosk
            derv2(jz,indm+k) = derv2(jz,indm+k) + d2(3,2)*drRi(2,k)*rfmass(jz)*rfmass(indi+k)*cosk
            derv2(jx,indm+k) = derv2(jx,indm+k) + d2(1,3)*drRi(3,k)*rfmass(jx)*rfmass(indi+k)*cosk
            derv2(jy,indm+k) = derv2(jy,indm+k) + d2(2,3)*drRi(3,k)*rfmass(jy)*rfmass(indi+k)*cosk
            derv2(jz,indm+k) = derv2(jz,indm+k) + d2(3,3)*drRi(3,k)*rfmass(jz)*rfmass(indi+k)*cosk
!
            dervi(jx,indm+k) = dervi(jx,indm+k) + d2(1,1)*drRi(1,k)*rfmass(jx)*rfmass(indi+k)*sink
            dervi(jy,indm+k) = dervi(jy,indm+k) + d2(2,1)*drRi(1,k)*rfmass(jy)*rfmass(indi+k)*sink
            dervi(jz,indm+k) = dervi(jz,indm+k) + d2(3,1)*drRi(1,k)*rfmass(jz)*rfmass(indi+k)*sink
            dervi(jx,indm+k) = dervi(jx,indm+k) + d2(1,2)*drRi(2,k)*rfmass(jx)*rfmass(indi+k)*sink
            dervi(jy,indm+k) = dervi(jy,indm+k) + d2(2,2)*drRi(2,k)*rfmass(jy)*rfmass(indi+k)*sink
            dervi(jz,indm+k) = dervi(jz,indm+k) + d2(3,2)*drRi(2,k)*rfmass(jz)*rfmass(indi+k)*sink
            dervi(jx,indm+k) = dervi(jx,indm+k) + d2(1,3)*drRi(3,k)*rfmass(jx)*rfmass(indi+k)*sink
            dervi(jy,indm+k) = dervi(jy,indm+k) + d2(2,3)*drRi(3,k)*rfmass(jy)*rfmass(indi+k)*sink
            dervi(jz,indm+k) = dervi(jz,indm+k) + d2(3,3)*drRi(3,k)*rfmass(jz)*rfmass(indi+k)*sink
          enddo
        enddo
!
!  Loop over second molecule
!
        do nn = 1,nmol
          jndm = 3*ncorenomol + 6*(nn - 1) + 3
!
!  Set up rotation matrix derivatives
!
          call setrotation(molaxes(1,1,nn),drS)
!
          do n = 1,nmolcore(nn)
            j = nmollist(nmolptr(nn)+n)
!
!  Multiply rotation derivatives about axes by vector to atom
!
            drRj(1:3,1:3) = 0.0_dp
            do l = 1,3
              do k = 1,3
                drRj(1,l) = drRj(1,l) + drS(1,k,l)*molxyz(k,n,nn)
                drRj(2,l) = drRj(2,l) + drS(2,k,l)*molxyz(k,n,nn)
                drRj(3,l) = drRj(3,l) + drS(3,k,l)*molxyz(k,n,nn)
              enddo
            enddo
!*********************************************
!  Mixed rotation - translation derivatives  *
!*********************************************
            jx = 3*ncorenomol + 6*(nn - 1)
!
            d2(1,1) = factor*qz(1,j)*qz(1,i)
            d2(2,1) = factor*qz(2,j)*qz(1,i)
            d2(3,1) = factor*qz(3,j)*qz(1,i)
            d2(1,2) = factor*qz(1,j)*qz(2,i)
            d2(2,2) = factor*qz(2,j)*qz(2,i)
            d2(3,2) = factor*qz(3,j)*qz(2,i)
            d2(1,3) = factor*qz(1,j)*qz(3,i)
            d2(2,3) = factor*qz(2,j)*qz(3,i)
            d2(3,3) = factor*qz(3,j)*qz(3,i)
!
!  Compute phase factor for this K point
!
            if (lphasecom) then
              cosk = xkv*(xclat(j)-xclat(i)-molxyz(1,n,nn)+molxyz(1,m,nm)) + &
                     ykv*(yclat(j)-yclat(i)-molxyz(2,n,nn)+molxyz(2,m,nm)) + &
                     zkv*(zclat(j)-zclat(i)-molxyz(3,n,nn)+molxyz(3,m,nm))
            else
              cosk = xkv*(xclat(j)-xclat(i)) + ykv*(yclat(j)-yclat(i)) + zkv*(zclat(j)-zclat(i))
            endif
            sink = sin(cosk)
            cosk = cos(cosk)
!
!  molQTdrv
!
            do ia = 1,3
              dQ(1) = d2(ia,1)
              dQ(2) = d2(ia,2)
              dQ(3) = d2(ia,3)
              do k = 1,3
                derv2(jx+ia,indm+k) = derv2(jx+ia,indm+k) + (drRi(1,k)*dQ(1) &
                                                            + drRi(2,k)*dQ(2) &
                                                            + drRi(3,k)*dQ(3)) &
                                                            *rfmass(indi+k)*rfmass(jx+ia)*cosk
 
                dervi(jx+ia,indm+k) = dervi(jx+ia,indm+k) + (drRi(1,k)*dQ(1) &
                                                            + drRi(2,k)*dQ(2) &
                                                            + drRi(3,k)*dQ(3)) &
                                                            *rfmass(indi+k)*rfmass(jx+ia)*sink
              enddo
            enddo
!************************************
!  Rotation - rotation derivatives  *
!************************************
!
!  molQQdrv
!
            do j1 = 1,3
              do k1 = 1,3
                derv2(jndm+k1,indm+j1) = derv2(jndm+k1,indm+j1) + (d2(1,1)*drRj(1,k1)*drRi(1,j1) &
                                                                + d2(2,1)*drRj(2,k1)*drRi(1,j1) &
                                                                + d2(3,1)*drRj(3,k1)*drRi(1,j1) &
                                                                + d2(1,2)*drRj(1,k1)*drRi(2,j1) &
                                                                + d2(2,2)*drRj(2,k1)*drRi(2,j1) &
                                                                + d2(3,2)*drRj(3,k1)*drRi(2,j1) &
                                                                + d2(1,3)*drRj(1,k1)*drRi(3,j1) &
                                                                + d2(2,3)*drRj(2,k1)*drRi(3,j1) &
                                                                + d2(3,3)*drRj(3,k1)*drRi(3,j1)) &
                                                                  *rfmass(jndm+k1)*rfmass(indi+j1)*cosk
!
                dervi(jndm+k1,indm+j1) = dervi(jndm+k1,indm+j1) + (d2(1,1)*drRj(1,k1)*drRi(1,j1) &
                                                                + d2(2,1)*drRj(2,k1)*drRi(1,j1) &
                                                                + d2(3,1)*drRj(3,k1)*drRi(1,j1) &
                                                                + d2(1,2)*drRj(1,k1)*drRi(2,j1) &
                                                                + d2(2,2)*drRj(2,k1)*drRi(2,j1) &
                                                                + d2(3,2)*drRj(3,k1)*drRi(2,j1) &
                                                                + d2(1,3)*drRj(1,k1)*drRi(3,j1) &
                                                                + d2(2,3)*drRj(2,k1)*drRi(3,j1) &
                                                                + d2(3,3)*drRj(3,k1)*drRi(3,j1)) &
                                                                  *rfmass(jndm+k1)*rfmass(indi+j1)*sink
              enddo
            enddo
          enddo
        enddo
      enddo
    endif
  enddo
!
!  Free local memory
!
  deallocate(qz,stat=status)
  if (status/=0) call deallocate_error('nagammacrigid','qz')
#ifdef TRACE
  call trace_out('nagammacrigid')
#endif
!
  return
  end
