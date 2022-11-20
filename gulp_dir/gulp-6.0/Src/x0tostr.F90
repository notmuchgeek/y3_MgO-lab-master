  subroutine x0tostr
!
!  Subroutine to convert linear structure array to main structure arrays
!
!   8/97 Created from part of energy.f
!  12/00 2-D modifications added
!  11/04 Inverse of cell parameters set here
!  11/06 lfirst argument added to equpos call
!   6/09 Module name changed from three to m_three
!  10/17 Modified so that absolute coordinates are not overwritten for MD
!   1/18 Trace added
!   2/19 x0 removed
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!  10/19 Rigid molecules added
!  12/19 Call to setquaternion modified for second derivatives
!   3/20 Molecule number added to setquaternion arguments
!   3/20 Fractional coordinates for dump added
!   4/20 molxyz now set instead of using local array ra
!   5/20 Set up of molaxes for current frame added
!   5/20 cart2frac changed to cart2frac_nomod
!   6/20 nmolcore changes added
!   7/20 Symmetry adaption of molecular quaternions changed
!        to allow for non-orthogonal cell case
!   7/20 Modifications for gfortran v10
!   7/20 Call to setquaternion replaced in non-symmetry case with setting molQsym to a unit matrix
!   7/20 Modified so that only atoms that relate to rigid molecules have their coordinates overwritten
!        after rigid molecule generation
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
  use configurations
  use control,          only : lrigid, lstraincell
  use current
  use four
  use m_three
  use mdlogic,          only : lmd
  use moldyn,           only : labscoany, labsco
  use molecule
  use species,          only : massspec
  use symmetry
#ifdef TRACE
  use trace,            only : trace_in, trace_out
#endif
  implicit none
!
!  Local variables
!
  integer(i4)            :: i
  integer(i4)            :: ia
  integer(i4)            :: ifail
  integer(i4)            :: indm
  integer(i4)            :: ix
  integer(i4)            :: iy
  integer(i4)            :: iz
  integer(i4)            :: j
  integer(i4)            :: n
  integer(i4)            :: nm
  integer(i4)            :: nma
  integer(i4)            :: nr
  integer(i4)            :: nsi
  real(dp)               :: rQ(3,3)
  real(dp)               :: drQ(3,3,3)
  real(dp)               :: drQ2(3,3,3,3)
  real(dp)               :: dx
  real(dp)               :: dy
  real(dp)               :: dz
  real(dp)               :: mI(3)
  real(dp)               :: rmassi
  real(dp)               :: rvtmp(3,3)
  real(dp)               :: xcom
  real(dp)               :: ycom
  real(dp)               :: zcom
  real(dp)               :: xci
  real(dp)               :: yci
  real(dp)               :: zci
  real(dp)               :: xfi
  real(dp)               :: yfi
  real(dp)               :: zfi
  real(dp)               :: wrk(9)
#ifdef TRACE
  call trace_in('x0tostr')
#endif
!
  if (ndim.eq.3) then
    call uncell3D(rv,a,b,c,alpha,beta,gamma)
    if (a.gt.1.0d-12) then
      recipa = 1.0_dp/a
    else
      recipa = 0.0_dp
    endif
    if (b.gt.1.0d-12) then
      recipb = 1.0_dp/b
    else
      recipb = 0.0_dp
    endif
    if (c.gt.1.0d-12) then
      recipc = 1.0_dp/c
    else
      recipc = 0.0_dp
    endif
    r1x = rv(1,1)
    r1y = rv(2,1)
    r1z = rv(3,1)
    r2x = rv(1,2)
    r2y = rv(2,2)
    r2z = rv(3,2)
    r3x = rv(1,3)
    r3y = rv(2,3)
    r3z = rv(3,3)
    call rlist
  elseif (ndim.eq.2) then
    call uncell2D(rv,a,b,alpha)
    if (a.gt.1.0d-12) then
      recipa = 1.0_dp/a
    else
      recipa = 0.0_dp
    endif
    if (b.gt.1.0d-12) then
      recipb = 1.0_dp/b
    else
      recipb = 0.0_dp
    endif
    r1x = rv(1,1)
    r1y = rv(2,1)
    r2x = rv(1,2)
    r2y = rv(2,2)
    call rlist
  elseif (ndim.eq.1) then
    call uncell1D(rv,a)
    if (a.gt.1.0d-12) then
      recipa = 1.0_dp/a
    else
      recipa = 0.0_dp
    endif
    r1x = rv(1,1)
    call rlist
  endif
!*************************************
!  Substitute parameters into place  *
!*************************************
!
!  Generate full coordinate set
!
  if (lsymopt) then
    call equpos(.true.,.false.)
  else
    do i = 1,numat
      xfrac(i) = xafrac(i)
      yfrac(i) = yafrac(i)
      zfrac(i) = zafrac(i)
    enddo
    if (nbsmat.gt.0) then
      do i = 1,numat
        radf(i) = rada(i)
      enddo
    endif
  endif
!
!  Convert cell parameters and internal coordinates into cartesian coordinates
!
  if (ndim.eq.3) then
    do i = 1,numat
      xclat(i) = xfrac(i)*r1x + yfrac(i)*r2x + zfrac(i)*r3x
      yclat(i) = xfrac(i)*r1y + yfrac(i)*r2y + zfrac(i)*r3y
      zclat(i) = xfrac(i)*r1z + yfrac(i)*r2z + zfrac(i)*r3z
    enddo
  elseif (ndim.eq.2) then
    do i = 1,numat
      xclat(i) = xfrac(i)*r1x + yfrac(i)*r2x
      yclat(i) = xfrac(i)*r1y + yfrac(i)*r2y
      zclat(i) = zfrac(i)
    enddo
  elseif (ndim.eq.1) then
    do i = 1,numat
      xclat(i) = xfrac(i)*r1x
      yclat(i) = yfrac(i)
      zclat(i) = zfrac(i)
    enddo
  else
    do i = 1,numat
      xclat(i) = xfrac(i)
      yclat(i) = yfrac(i)
      zclat(i) = zfrac(i)
    enddo
  endif
  if (lsymopt) then
    do i = 1,nasym
      nr = nrela2f(i)
      xalat(i) = xclat(nr)
      yalat(i) = yclat(nr)
      zalat(i) = zclat(nr)
    enddo
  else
!
!  Avoid overwriting absolute coordinates for MD
!
    if (lmd.and.labscoany) then
      do i = 1,nasym
        if (.not.labsco(i)) then
          xalat(i) = xclat(i)
          yalat(i) = yclat(i)
          zalat(i) = zclat(i)
        endif
      enddo
    else
      do i = 1,nasym
        xalat(i) = xclat(i)
        yalat(i) = yclat(i)
        zalat(i) = zclat(i)
      enddo
    endif
  endif
  if (lrigid) then
!
!  Rigid molecules
!
    if (lsymopt) then
!
!  Symmetry generation of centres of mass and quaternions
!
      call symupdatemol
    else
      do nm = 1,nmol
        molcom(1:3,nm) = molcoma(1:3,nm)
        molQ(1:3,nm) = molQa(1:3,nm)
!
!  Set rotation matrix for molecule
!
        molQsym(1:3,1:3,nm) = 0.0_dp
        molQsym(1,1,nm) = 1.0_dp
        molQsym(2,2,nm) = 1.0_dp
        molQsym(3,3,nm) = 1.0_dp
      enddo
    endif
!
!  Copy initial fractional coordinates for dumping
!
    do i = 1,numat
      xfdmp(i) = xfrac(i)
      yfdmp(i) = yfrac(i)
      zfdmp(i) = zfrac(i)
    enddo
!
!  Build full molecules
!
    do nm = 1,nmol
!
!  Set symmetry equivalent molecule in the asymmetric unit
!
      nma = nmolf2a(nm)
!
!  Convert fractional coordinates for COM to Cartesian
!
      if (ndim.eq.3) then
        xcom = molcom(1,nm)*r1x + molcom(2,nm)*r2x + molcom(3,nm)*r3x
        ycom = molcom(1,nm)*r1y + molcom(2,nm)*r2y + molcom(3,nm)*r3y
        zcom = molcom(1,nm)*r1z + molcom(2,nm)*r2z + molcom(3,nm)*r3z
      elseif (ndim.eq.2) then
        xcom = molcom(1,nm)*r1x + molcom(2,nm)*r2x
        ycom = molcom(1,nm)*r1y + molcom(2,nm)*r2y
        zcom = molcom(3,nm)
      elseif (ndim.eq.1) then
        xcom = molcom(1,nm)*r1x
        ycom = molcom(2,nm)
        zcom = molcom(3,nm)
      else
        xcom = molcom(1,nm)
        ycom = molcom(2,nm)
        zcom = molcom(3,nm)
      endif
!
!  Set rotation matrix for molecule
!
      call setquaternion(nm,molQa(1,nma),rQ,drQ,drQ2,.false.,.false.)
!
      do n = 1,nmolcore(nm)
        i = nmollist(nmolptr(nm)+n)
!
!  Find cell image for atom
!
        indm = nmolind(i)
        call mindtoijk(indm,ix,iy,iz)
!
!  Apply quaternions to molecule orientation - use asymmetric unit reference to maintain symmetry
!
        do j = 1,3
          molxyz(j,n,nm) = rQ(j,1)*molQxyz(1,n,nma) + &
                           rQ(j,2)*molQxyz(2,n,nma) + &
                           rQ(j,3)*molQxyz(3,n,nma)
        enddo
!
        xclat(i) = xcom + molxyz(1,n,nm) - dble(ix)*rv(1,1) &
                                         - dble(iy)*rv(1,2) &
                                         - dble(iz)*rv(1,3)
        yclat(i) = ycom + molxyz(2,n,nm) - dble(ix)*rv(2,1) &
                                         - dble(iy)*rv(2,2) &
                                         - dble(iz)*rv(2,3)
        zclat(i) = zcom + molxyz(3,n,nm) - dble(ix)*rv(3,1) &
                                         - dble(iy)*rv(3,2) &
                                         - dble(iz)*rv(3,3)
!
!  Convert molecule atom coordinates back to fractional
!
        xci = xclat(i)
        yci = yclat(i)
        zci = zclat(i)
        call cart2frac_nomod(ndim,xci,yci,zci,rv,xfi,yfi,zfi)
        xfrac(i) = xfi
        yfrac(i) = yfi
        zfrac(i) = zfi
!
!  Handle finite strains
!
        if (lstraincell) then
          xci = xclat(i) - xcom
          yci = yclat(i) - ycom
          zci = zclat(i) - zcom
          rvtmp(1:3,1:3) = rvcfg(1:3,1:3,ncf)
          call cart2frac_nomod(ndim,xci,yci,zci,rvtmp,xfi,yfi,zfi)
          xfdmp(i) = molcom(1,nm) + xfi
          yfdmp(i) = molcom(2,nm) + yfi
          zfdmp(i) = molcom(3,nm) + zfi
        else
          xfdmp(i) = xfrac(i)
          yfdmp(i) = yfrac(i)
          zfdmp(i) = zfrac(i)
        endif
!
!  Set asymmetric unit coordinates
!
        ia = nrelf2a(i)
        if (nrela2f(ia).eq.i) then
          xafrac(ia) = xfrac(i)
          yafrac(ia) = yfrac(i)
          zafrac(ia) = zfrac(i)
!
          xalat(ia) = xclat(i)
          yalat(ia) = yclat(i)
          zalat(ia) = zclat(i)
        endif
      enddo
    enddo
!
!  Set current frame for molecular axes - needed for properties and phonons
!
    do nm = 1,nmol
      molaxes(1:3,1:3,nm) = 0.0_dp
      do n = 1,nmolcore(nm)
        i = nmollist(nmolptr(nm)+n)
        nsi = nspecptr(nrelf2a(i))
        rmassi = massspec(nsi)*occuf(i)
!
        dx = molxyz(1,n,nm)
        dy = molxyz(2,n,nm)
        dz = molxyz(3,n,nm)
!
        molaxes(1,1,nm) = molaxes(1,1,nm) + rmassi*dy*dy + rmassi*dz*dz
        molaxes(2,1,nm) = molaxes(2,1,nm) - rmassi*dx*dy
        molaxes(3,1,nm) = molaxes(3,1,nm) - rmassi*dx*dz
        molaxes(1,2,nm) = molaxes(1,2,nm) - rmassi*dy*dx
        molaxes(2,2,nm) = molaxes(2,2,nm) + rmassi*dx*dx + rmassi*dz*dz
        molaxes(3,2,nm) = molaxes(3,2,nm) - rmassi*dy*dz
        molaxes(1,3,nm) = molaxes(1,3,nm) - rmassi*dz*dx
        molaxes(2,3,nm) = molaxes(2,3,nm) - rmassi*dz*dy
        molaxes(3,3,nm) = molaxes(3,3,nm) + rmassi*dx*dx + rmassi*dy*dy
      enddo
!
!  Diagonalise tensor to get axes for moment of inertia tensor
!
      call dsyev('V','U',3_i4,molaxes(1,1,nm),3_i4,mI,wrk,9_i4,ifail)
    enddo
  endif
#ifdef TRACE
  call trace_out('x0tostr')
#endif
!
  return
  end
