  subroutine nagammac(ncfoc,nphonatc,nphonatptr,iocptr,rfmass,qfrac,maxd2,derv2,dervi,xk,yk,zk)
!
!  Calculates the non-analytic correction to the second
!  derivatives based on the Born effective charges.
!  Complex version.
!
!   7/13 Created from nagamma
!  12/16 Modified for parallel second derivatives
!   1/18 Modified to handle ghostcell keyword
!   2/18 Trace added
!   3/20 Location of angstoev changed to current
!   4/20 Mass arrays now use fmass and rfmass
!
!  On entry :
!
!  ncfoc      = no. of reduced core sites 
!  nphonatc   = no. of cores included in phonon calculation
!  nphonatptr = pointer from nphonatc to atom number 
!  iocptr     = pointer to core sites from reduced list
!  rfmass     = inverse square root of atomic masses
!  qfrac      = fractional approach direction to gamma point
!  maxd2      = left-hand dimension of derv2
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
!  Julian Gale, CIC, Curtin University, April 2020
!
  use configurations, only : nsuperghost
  use g_constants
  use control,        only : lkfull, lghost
  use current,        only : bornq, kv, rv, xclat, yclat, zclat, occuf, ncf, angstoev
  use parallel
  use phononatoms,    only : nphonatonnodec, nphonatonnodecptr
  use properties
  use shells
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4),                    intent(in)    :: iocptr(*)
  integer(i4),                    intent(in)    :: maxd2
  integer(i4),                    intent(in)    :: ncfoc
  integer(i4),                    intent(in)    :: nphonatc
  integer(i4),                    intent(in)    :: nphonatptr(nphonatc)
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
  integer(i4)                                   :: ii
  integer(i4)                                   :: indi
  integer(i4)                                   :: ip
  integer(i4)                                   :: ix
  integer(i4)                                   :: iy
  integer(i4)                                   :: iz
  integer(i4)                                   :: j
  integer(i4)                                   :: jj
  integer(i4)                                   :: jx
  integer(i4)                                   :: jy
  integer(i4)                                   :: jz
  integer(i4)                                   :: nghostcell
  integer(i4)                                   :: status
  real(dp),   dimension(:,:), allocatable       :: qz
  real(dp)                                      :: cosk
  real(dp)                                      :: factor
  real(dp)                                      :: inveps
  real(dp)                                      :: kvf(3,3)
  real(dp)                                      :: mfactor
  real(dp)                                      :: oci
  real(dp)                                      :: ocj
  real(dp)                                      :: qcart(3)
  real(dp)                                      :: sink
  real(dp)                                      :: vol
  real(dp)                                      :: volume
  real(dp)                                      :: xkv
  real(dp)                                      :: ykv
  real(dp)                                      :: zkv
#ifdef TRACE
  call trace_in('nagammac')
#endif
!
!  Calculate constant = 4*pi/V
!
  vol = volume(rv)
  factor = 4.0_dp*pi/vol
!
!  Find number of ghost cells
!
  if (lghost) then
    nghostcell = nsuperghost(1,ncf)*nsuperghost(2,ncf)*nsuperghost(3,ncf)
  else
    nghostcell = 1
  endif
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
  allocate(qz(3,ncfoc),stat=status)
  if (status/=0) call outofmemory('nagamma','qz')
  qz(1:3,1:ncfoc) = 0.0_dp
  do i = 1,nphonatc,nghostcell
    if (lghost) then
      ii = (iocptr(i) - 1)/nghostcell + 1
    else
      ii = iocptr(i)
    endif
    qz(1,ii) = qz(1,ii) + qcart(1)*bornq(1,1,i) + qcart(2)*bornq(2,1,i) + qcart(3)*bornq(3,1,i)
    qz(2,ii) = qz(2,ii) + qcart(1)*bornq(1,2,i) + qcart(2)*bornq(2,2,i) + qcart(3)*bornq(3,2,i)
    qz(3,ii) = qz(3,ii) + qcart(1)*bornq(1,3,i) + qcart(2)*bornq(2,3,i) + qcart(3)*bornq(3,3,i)
  enddo
!
!  Calculate Q.epsilon(static).Q and invert
!
  inveps = qcart(1)*diconh(1,1)*qcart(1) + qcart(2)*diconh(2,1)*qcart(1) + qcart(3)*diconh(3,1)*qcart(1) + &
           qcart(1)*diconh(1,2)*qcart(2) + qcart(2)*diconh(2,2)*qcart(2) + qcart(3)*diconh(3,2)*qcart(2) + &
           qcart(1)*diconh(1,3)*qcart(3) + qcart(2)*diconh(2,3)*qcart(3) + qcart(3)*diconh(3,3)*qcart(3)
  inveps = 1.0_dp/inveps
!
!  Compute point in BZ in Cartesian form
!
  xkv = xk*kvf(1,1) + yk*kvf(1,2) + zk*kvf(1,3)
  ykv = xk*kvf(2,1) + yk*kvf(2,2) + zk*kvf(2,3)
  zkv = xk*kvf(3,1) + yk*kvf(3,2) + zk*kvf(3,3)
!
!  Loop over second derivative matrix elements adding terms
!
  factor = factor*inveps*angstoev
  if (nprocs.gt.1) then
    ix = -2
    iy = -1
    iz =  0
    do ip = 1,nphonatonnodec,nghostcell
      i  = nphonatonnodecptr(ip)
      ii = iocptr(i)
      indi = 3*(i-1)
! DEBUG - the following needs handling for partial occupancy case
! DEBUG - ghost case needs handling in parallel
      ix = ix + 3
      iy = iy + 3
      iz = iz + 3
      oci = occuf(i)
      do j = 1,ncfoc
        jj = iocptr(j)
        jx = 3*(jj-1) + 1
        jy = 3*(jj-1) + 2
        jz = 3*(jj-1) + 3
        ocj = occuf(j)
!
!  Compute phase factor for this K point
!
        cosk = xkv*(xclat(jj)-xclat(i)) + ykv*(yclat(jj)-yclat(i)) + zkv*(zclat(jj)-zclat(i))
        sink = sin(cosk)
        cosk = cos(cosk)
!
        mfactor = factor*oci*ocj
!
        derv2(jx,ix) = derv2(jx,ix) + mfactor*rfmass(indi+1)*rfmass(jx)*qz(1,j)*qz(1,i)*cosk
        derv2(jy,ix) = derv2(jy,ix) + mfactor*rfmass(indi+1)*rfmass(jy)*qz(2,j)*qz(1,i)*cosk
        derv2(jz,ix) = derv2(jz,ix) + mfactor*rfmass(indi+1)*rfmass(jz)*qz(3,j)*qz(1,i)*cosk
        derv2(jx,iy) = derv2(jx,iy) + mfactor*rfmass(indi+2)*rfmass(jx)*qz(1,j)*qz(2,i)*cosk
        derv2(jy,iy) = derv2(jy,iy) + mfactor*rfmass(indi+2)*rfmass(jy)*qz(2,j)*qz(2,i)*cosk
        derv2(jz,iy) = derv2(jz,iy) + mfactor*rfmass(indi+2)*rfmass(jz)*qz(3,j)*qz(2,i)*cosk
        derv2(jx,iz) = derv2(jx,iz) + mfactor*rfmass(indi+3)*rfmass(jx)*qz(1,j)*qz(3,i)*cosk
        derv2(jy,iz) = derv2(jy,iz) + mfactor*rfmass(indi+3)*rfmass(jy)*qz(2,j)*qz(3,i)*cosk
        derv2(jz,iz) = derv2(jz,iz) + mfactor*rfmass(indi+3)*rfmass(jz)*qz(3,j)*qz(3,i)*cosk
!
        dervi(jx,ix) = dervi(jx,ix) + mfactor*rfmass(indi+1)*rfmass(jx)*qz(1,j)*qz(1,i)*sink
        dervi(jy,ix) = dervi(jy,ix) + mfactor*rfmass(indi+1)*rfmass(jy)*qz(2,j)*qz(1,i)*sink
        dervi(jz,ix) = dervi(jz,ix) + mfactor*rfmass(indi+1)*rfmass(jz)*qz(3,j)*qz(1,i)*sink
        dervi(jx,iy) = dervi(jx,iy) + mfactor*rfmass(indi+2)*rfmass(jx)*qz(1,j)*qz(2,i)*sink
        dervi(jy,iy) = dervi(jy,iy) + mfactor*rfmass(indi+2)*rfmass(jy)*qz(2,j)*qz(2,i)*sink
        dervi(jz,iy) = dervi(jz,iy) + mfactor*rfmass(indi+2)*rfmass(jz)*qz(3,j)*qz(2,i)*sink
        dervi(jx,iz) = dervi(jx,iz) + mfactor*rfmass(indi+3)*rfmass(jx)*qz(1,j)*qz(3,i)*sink
        dervi(jy,iz) = dervi(jy,iz) + mfactor*rfmass(indi+3)*rfmass(jy)*qz(2,j)*qz(3,i)*sink
        dervi(jz,iz) = dervi(jz,iz) + mfactor*rfmass(indi+3)*rfmass(jz)*qz(3,j)*qz(3,i)*sink
      enddo
    enddo
  else
    ix = - 2
    iy = - 1
    iz =   0
    do i = 1,ncfoc
      ii = iocptr(i)
      ix = ix + 3
      iy = iy + 3
      iz = iz + 3
      jx = - 2
      jy = - 1
      jz =   0
      oci = occuf(ii)
      do j = 1,ncfoc
        jj = iocptr(j)
        jx = jx + 3
        jy = jy + 3
        jz = jz + 3
        ocj = occuf(jj)
!
!  Compute phase factor for this K point
!
        cosk = xkv*(xclat(jj)-xclat(ii)) + ykv*(yclat(jj)-yclat(ii)) + zkv*(zclat(jj)-zclat(ii))
        sink = sin(cosk)
        cosk = cos(cosk)
!
        mfactor = factor*oci*ocj
!
        derv2(jx,ix) = derv2(jx,ix) + mfactor*rfmass(ix)*rfmass(jx)*qz(1,j)*qz(1,i)*cosk
        derv2(jy,ix) = derv2(jy,ix) + mfactor*rfmass(ix)*rfmass(jy)*qz(2,j)*qz(1,i)*cosk
        derv2(jz,ix) = derv2(jz,ix) + mfactor*rfmass(ix)*rfmass(jz)*qz(3,j)*qz(1,i)*cosk
        derv2(jx,iy) = derv2(jx,iy) + mfactor*rfmass(iy)*rfmass(jx)*qz(1,j)*qz(2,i)*cosk
        derv2(jy,iy) = derv2(jy,iy) + mfactor*rfmass(iy)*rfmass(jy)*qz(2,j)*qz(2,i)*cosk
        derv2(jz,iy) = derv2(jz,iy) + mfactor*rfmass(iy)*rfmass(jz)*qz(3,j)*qz(2,i)*cosk
        derv2(jx,iz) = derv2(jx,iz) + mfactor*rfmass(iz)*rfmass(jx)*qz(1,j)*qz(3,i)*cosk
        derv2(jy,iz) = derv2(jy,iz) + mfactor*rfmass(iz)*rfmass(jy)*qz(2,j)*qz(3,i)*cosk
        derv2(jz,iz) = derv2(jz,iz) + mfactor*rfmass(iz)*rfmass(jz)*qz(3,j)*qz(3,i)*cosk
!
        dervi(jx,ix) = dervi(jx,ix) + mfactor*rfmass(ix)*rfmass(jx)*qz(1,j)*qz(1,i)*sink
        dervi(jy,ix) = dervi(jy,ix) + mfactor*rfmass(ix)*rfmass(jy)*qz(2,j)*qz(1,i)*sink
        dervi(jz,ix) = dervi(jz,ix) + mfactor*rfmass(ix)*rfmass(jz)*qz(3,j)*qz(1,i)*sink
        dervi(jx,iy) = dervi(jx,iy) + mfactor*rfmass(iy)*rfmass(jx)*qz(1,j)*qz(2,i)*sink
        dervi(jy,iy) = dervi(jy,iy) + mfactor*rfmass(iy)*rfmass(jy)*qz(2,j)*qz(2,i)*sink
        dervi(jz,iy) = dervi(jz,iy) + mfactor*rfmass(iy)*rfmass(jz)*qz(3,j)*qz(2,i)*sink
        dervi(jx,iz) = dervi(jx,iz) + mfactor*rfmass(iz)*rfmass(jx)*qz(1,j)*qz(3,i)*sink
        dervi(jy,iz) = dervi(jy,iz) + mfactor*rfmass(iz)*rfmass(jy)*qz(2,j)*qz(3,i)*sink
        dervi(jz,iz) = dervi(jz,iz) + mfactor*rfmass(iz)*rfmass(jz)*qz(3,j)*qz(3,i)*sink
      enddo
    enddo
  endif
!
!  Free local memory
!
  deallocate(qz,stat=status)
  if (status/=0) call deallocate_error('nagamma','qz')
#ifdef TRACE
  call trace_out('nagammac')
#endif
!
  return
  end
