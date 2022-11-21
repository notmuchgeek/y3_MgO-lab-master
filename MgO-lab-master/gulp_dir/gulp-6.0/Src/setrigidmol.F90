  subroutine setrigidmol
!
!  Sets up quantities for rigid molecules including:
!
!  1) Centre of molecule - stored in molcom
!  2) Reference orientation - stored in molxyz
!  3) Quaternions
!
!  NB: For centred cells with symmetry: molcom has centre of mass as fractions of primitive cell
!                                       molcoma has centre of mass as fractions of full cell
!
!   2/19 Created
!  10/19 Symmetry modifications added
!  10/19 molcom now expressed in fractional where appropriate
!   3/20 lreorient replaced by lmolstdframe using keyword
!   3/20 Variable atom order handling added
!   3/20 Separate array introduce for molxyz in quaternion frame of reference
!   3/20 Modified to handle finite strains
!   4/20 lmolcom_mass flag added to control whether centre is mass weighted or not
!   4/20 Restart handling added
!   4/20 Centre of mass for restart now in fractional
!   5/20 Standard frame removed as no longer required
!   6/20 nmolcore changes added
!   7/20 molQsym setting removed as this is now handled in symupdatemol
!   7/20 Handling of centred cells with symmetry added for centre of mass
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
!  Julian Gale, CIC, Curtin University, July 2020
!
  use configurations,     only : rvcfg
  use control
  use current
  use element
  use iochannels
  use molecule
  use parallel
  use reallocate
  use species,            only : massspec
  use symmetry
  implicit none
!
!  Local variables
!
  integer(i4)                                    :: i
  integer(i4)                                    :: icosxi
  integer(i4)                                    :: icosyi
  integer(i4)                                    :: icoszi
  integer(i4)                                    :: ifail
  integer(i4)                                    :: indm
  integer(i4)                                    :: ix
  integer(i4)                                    :: iy
  integer(i4)                                    :: iz
  integer(i4)                                    :: j
  integer(i4)                                    :: mv
  integer(i4)                                    :: n
  integer(i4)                                    :: n0
  integer(i4)                                    :: nm
  real(dp)                                       :: inertia(3,3)
  real(dp)                                       :: dx
  real(dp)                                       :: dy
  real(dp)                                       :: dz
  real(dp)                                       :: fp(3)
  real(dp)                                       :: mi
  real(dp)                                       :: totalmass
  real(dp)                                       :: w1inv(3,3)
  real(dp)                                       :: wrk(9)
  real(dp)                                       :: xci
  real(dp)                                       :: yci
  real(dp)                                       :: zci
  real(dp)                                       :: xcom
  real(dp)                                       :: ycom
  real(dp)                                       :: zcom
!
!  Set up inverse of centring matrix if required
!
  if (lsymopt) then
    if (ngocfg(ncf).le.1) then
      w1inv(1:3,1:3) = w1(ncbl,1:3,1:3)
      call matrix_inversion(w1inv,3_i4,3_i4,wrk,ifail)
    endif
  endif
!************************
!  Loop over molecules  *
!************************
  n0 = 0
  do nm = 1,nmol
    if (lmolrigid) then
!############################
!  Use restart information  #
!############################
!
!  Set quaternion components and centre of mass
!
      molQ(1:3,nm) = molQcfg(1:3,nm,ncf)
      molcom(1:3,nm) = molcomcfg(1:3,nm,ncf)
!
!  Set reference frame for coordinates
!
      do n = 1,nmolcore(nm)
        do ix = 1,3
          molQxyz(ix,n,nm) = molQxyzcfg(ix,n0+n,ncf)
          molxyz(ix,n,nm) = molQxyz(ix,n,nm)
        enddo
      enddo
!
!  Increment total number of atoms across all molecules so far
!
      n0 = n0 + nmolatom(nm)
!
!  Compute moment of inertia tensor
!
      inertia(1:3,1:3) = 0.0_dp
      do n = 1,nmolcore(nm)
        i = nmollist(nmolptr(nm)+n)
        mi = occuf(i)*massspec(nspecptr(nrelf2a(i)))
        dx = molxyz(1,n,nm)
        dy = molxyz(2,n,nm)
        dz = molxyz(3,n,nm)
        inertia(1,1) = inertia(1,1) + mi*dy*dy + mi*dz*dz
        inertia(2,1) = inertia(2,1) - mi*dx*dy
        inertia(3,1) = inertia(3,1) - mi*dx*dz
        inertia(1,2) = inertia(1,2) - mi*dy*dx
        inertia(2,2) = inertia(2,2) + mi*dx*dx + mi*dz*dz
        inertia(3,2) = inertia(3,2) - mi*dy*dz
        inertia(1,3) = inertia(1,3) - mi*dz*dx
        inertia(2,3) = inertia(2,3) - mi*dz*dy
        inertia(3,3) = inertia(3,3) + mi*dx*dx + mi*dy*dy
      enddo
!
!  Diagonalise tensor to get principal component moments of inertia
!
      call dsyev('V','U',3_i4,inertia,3_i4,molI(1,nm),wrk,9_i4,ifail)
    else
!#####################################################
!  Compute reference orientation and centre of mass  #
!#####################################################
!
!  Set up coordinates of atoms within the molecule
!
      xcom = 0.0_dp
      ycom = 0.0_dp
      zcom = 0.0_dp
      totalmass = 0.0_dp
!
      do n = 1,nmolcore(nm)
        i = nmollist(nmolptr(nm)+n)
!
!  To handle finite strain case generate coordinates in unstrained cell
!
        if (ndim.eq.3) then
          molxyz(1,n,nm) = xfrac(i)*rvcfg(1,1,ncf) + yfrac(i)*rvcfg(1,2,ncf) + zfrac(i)*rvcfg(1,3,ncf)
          molxyz(2,n,nm) = xfrac(i)*rvcfg(2,1,ncf) + yfrac(i)*rvcfg(2,2,ncf) + zfrac(i)*rvcfg(2,3,ncf)
          molxyz(3,n,nm) = xfrac(i)*rvcfg(3,1,ncf) + yfrac(i)*rvcfg(3,2,ncf) + zfrac(i)*rvcfg(3,3,ncf)
        elseif (ndim.eq.2) then
          molxyz(1,n,nm) = xfrac(i)*rvcfg(1,1,ncf) + yfrac(i)*rvcfg(1,2,ncf)
          molxyz(2,n,nm) = xfrac(i)*rvcfg(2,1,ncf) + yfrac(i)*rvcfg(2,2,ncf)
          molxyz(3,n,nm) = zfrac(i)
        elseif (ndim.eq.1) then
          molxyz(1,n,nm) = xfrac(i)*rvcfg(1,1,ncf)
          molxyz(2,n,nm) = yfrac(i)
          molxyz(3,n,nm) = zfrac(i)
        else
          molxyz(1,n,nm) = xfrac(i)
          molxyz(2,n,nm) = yfrac(i)
          molxyz(3,n,nm) = zfrac(i)
        endif
!
!  Find cell image for atom
!
        indm = nmolind(i)
        call mindtoijk(indm,ix,iy,iz)
!
        molxyz(1,n,nm) = molxyz(1,n,nm) + dble(ix)*rv(1,1) &
                                        + dble(iy)*rv(1,2) &
                                        + dble(iz)*rv(1,3)
        molxyz(2,n,nm) = molxyz(2,n,nm) + dble(ix)*rv(2,1) &
                                        + dble(iy)*rv(2,2) &
                                        + dble(iz)*rv(2,3)
        molxyz(3,n,nm) = molxyz(3,n,nm) + dble(ix)*rv(3,1) &
                                        + dble(iy)*rv(3,2) &
                                        + dble(iz)*rv(3,3)
!
        xci = xclat(i) + dble(ix)*rv(1,1) &
                       + dble(iy)*rv(1,2) &
                       + dble(iz)*rv(1,3)
        yci = yclat(i) + dble(ix)*rv(2,1) &
                       + dble(iy)*rv(2,2) &
                       + dble(iz)*rv(2,3)
        zci = zclat(i) + dble(ix)*rv(3,1) &
                       + dble(iy)*rv(3,2) &
                       + dble(iz)*rv(3,3)
!
!  Centre of molecule for real position
!
        if (lmolcom_mass) then
          mi = occuf(i)*massspec(nspecptr(nrelf2a(i)))
          xcom = xcom + xci*mi
          ycom = ycom + yci*mi
          zcom = zcom + zci*mi
          totalmass = totalmass + mi
        else
          xcom = xcom + xci
          ycom = ycom + yci
          zcom = zcom + zci
        endif
      enddo
!
!  Find centre of mass for reference position
!
      molcom(1:3,nm) = 0.0_dp
      if (lmolcom_mass) then
        do n = 1,nmolcore(nm)
          i = nmollist(nmolptr(nm)+n)
          mi = occuf(i)*massspec(nspecptr(nrelf2a(i)))
          molcom(1,nm) = molcom(1,nm) + molxyz(1,n,nm)*mi
          molcom(2,nm) = molcom(2,nm) + molxyz(2,n,nm)*mi
          molcom(3,nm) = molcom(3,nm) + molxyz(3,n,nm)*mi
        enddo
        molcom(1:3,nm) = molcom(1:3,nm)/totalmass
      else
        do n = 1,nmolcore(nm)
          molcom(1,nm) = molcom(1,nm) + molxyz(1,n,nm)
          molcom(2,nm) = molcom(2,nm) + molxyz(2,n,nm)
          molcom(3,nm) = molcom(3,nm) + molxyz(3,n,nm)
        enddo
        molcom(1:3,nm) = molcom(1:3,nm)/dble(nmolcore(nm))
      endif
!
!  Find centre of mass for actual position
!
      if (lmolcom_mass) then
        xcom = xcom/totalmass
        ycom = ycom/totalmass
        zcom = zcom/totalmass
      else
        xcom = xcom/dble(nmolcore(nm))
        ycom = ycom/dble(nmolcore(nm))
        zcom = zcom/dble(nmolcore(nm))
      endif
!
!  Subtract COM position from atom coordinates based on reference frame
!
      do n = 1,nmolcore(nm)
        molxyz(1,n,nm) = molxyz(1,n,nm) - molcom(1,nm) 
        molxyz(2,n,nm) = molxyz(2,n,nm) - molcom(2,nm)
        molxyz(3,n,nm) = molxyz(3,n,nm) - molcom(3,nm)
      enddo
!
!  Compute moment of inertia tensor
!
      inertia(1:3,1:3) = 0.0_dp
      do n = 1,nmolcore(nm)
        mi = occuf(i)*massspec(nspecptr(nrelf2a(i)))
        dx = molxyz(1,n,nm)
        dy = molxyz(2,n,nm)
        dz = molxyz(3,n,nm)
        inertia(1,1) = inertia(1,1) + mi*dy*dy + mi*dz*dz
        inertia(2,1) = inertia(2,1) - mi*dx*dy
        inertia(3,1) = inertia(3,1) - mi*dx*dz
        inertia(1,2) = inertia(1,2) - mi*dy*dx
        inertia(2,2) = inertia(2,2) + mi*dx*dx + mi*dz*dz
        inertia(3,2) = inertia(3,2) - mi*dy*dz
        inertia(1,3) = inertia(1,3) - mi*dz*dx
        inertia(2,3) = inertia(2,3) - mi*dz*dy
        inertia(3,3) = inertia(3,3) + mi*dx*dx + mi*dy*dy
      enddo
!
!  Diagonalise tensor to get principal component moments of inertia
!
      call dsyev('V','U',3_i4,inertia,3_i4,molI(1,nm),wrk,9_i4,ifail)
!
!  Save eigenvectors of inertia tensor to set reference frame for moment of inertia
!
      molQeig(1:3,1:3,nm) = inertia(1:3,1:3)
!
!  Reset COM based on actual strained position
!
      molcom(1,nm) = xcom
      molcom(2,nm) = ycom
      molcom(3,nm) = zcom
!
!  Initialise restart info for the molecule present
!
      molcomcfg(1:3,nm,ncf) = molcom(1:3,nm)
!
!  Convert COM from Cartesian to fractional
!
      call cart2frac(ndim,molcom(1,nm),molcom(2,nm),molcom(3,nm),rv,molcom(1,nm),molcom(2,nm),molcom(3,nm),icosxi,icosyi,icoszi)
!
!  Initialise quaternion components as zero for original orientation
!
      molQ(1:3,nm) = 0.0_dp
!
!  Copy initial orientation into the frame for quaternions
!
      do n = 1,nmolcore(nm)
        do ix = 1,3
          molQxyz(ix,n,nm) = molxyz(ix,n,nm)
        enddo
      enddo
!
!  Set up molecule info in configuration arrays if necessary
!
      if (nmolcfg(ncf).lt.nmol) then
        nmolcfg(ncf) = nmolcfg(ncf) + 1
        nmolatomcfg(nmolcfg(ncf),ncf) = nmolatom(nm)
        nmolcorecfg(nmolcfg(ncf),ncf) = nmolcore(nm)
        do i = 1,nmolatom(nm)
          nmollistcfg(n0+i,ncf) = nmollist(nmolptr(nm)+i)
        enddo
        nmolatomtotcfg(ncf) = nmolatomtotcfg(ncf) + nmolatom(nm)
      endif
!
!  Place reference frame for coordinates and initial quaternions into restart arrays
!
      molQcfg(1:3,nm,ncf) = molQ(1:3,nm)
      do n = 1,nmolcore(nm)
        do ix = 1,3
          molQxyzcfg(ix,n0+n,ncf) = molQxyz(ix,n,nm)
        enddo
      enddo
      n0 = n0 + nmolatom(nm)
    endif
  enddo
!
!  Set values for asymmetric unit
!
  do n = 1,nmolasym
    nm = nmola2f(n)
!
!  Handle centring for symmetry case if required
!
    if (lsymopt) then
      if (ngocfg(ncf).le.1) then
        fp(1:3) = molcom(1:3,nm)
!
!  The use of the shift is unlikely to be necessary as the first atom tends not to have one
!
        mv = nrotop(nmollist(nmolptr(nm)+1))
        molcoma(1:3,n) = - vit(1:3,mv)
!
        do i = 1,3
          do j = 1,3
            molcoma(i,n) = molcoma(i,n) + w1inv(i,j)*fp(j)
          enddo
        enddo
      endif
    else
      molcoma(1:3,n) = molcom(1:3,nm)
    endif
    molQa(1:3,n) = molQ(1:3,nm)
  enddo
!
  return
  end
