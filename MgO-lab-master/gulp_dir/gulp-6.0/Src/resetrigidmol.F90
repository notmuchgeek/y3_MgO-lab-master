  subroutine resetrigidmol(nm)
!
!  Resets orientation for rigid molecules based on last successfully generated set of coordinates
!  This routine is called when the original quaternions have become unphysical
!  i.e. q0 is complex
!
!   7/20 Created from setrigidmol
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
  use current
  use element
  use molecule
  use parallel
  use species,            only : massspec
  implicit none
!
!  Passed variables
!
  integer(i4),                    intent(in)     :: nm   ! Molecule whose quaternions need to be reset
!
!  Local variables
!
  integer(i4)                                    :: i
  integer(i4)                                    :: ifail
  integer(i4)                                    :: ix
  integer(i4)                                    :: n
  integer(i4)                                    :: n0
  integer(i4)                                    :: nma
  real(dp)                                       :: dx
  real(dp)                                       :: dy
  real(dp)                                       :: dz
  real(dp)                                       :: inertia(3,3)
  real(dp)                                       :: mi
  real(dp)                                       :: wrk(9)
!############################################################################
!  Set current orientation as the reference orientation and centre of mass  #
!############################################################################
!
!  Copy current orientation to be the initial frame for quaternions
!
  do n = 1,nmolcore(nm)
    molQxyz(1:3,n,nm) = molxyz(1:3,n,nm)
  enddo
!
!  Recompute moment of inertia tensor for current orientation
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
!
!  Save eigenvectors of inertia tensor to set reference frame for moment of inertia
!
  molQeig(1:3,1:3,nm) = inertia(1:3,1:3)
!
!  Reinitialise quaternion components as zero for current orientation
!
  molQ(1:3,nm) = 0.0_dp
!
!  Place new reference frame for coordinates and initial quaternions into restart arrays
!
  molQcfg(1:3,nm,ncf) = molQ(1:3,nm)
  n0 = 0
  do n = 1,nm-1
    n0 = n0 + nmolatom(n)
  enddo
  do n = 1,nmolcore(nm)
    do ix = 1,3
      molQxyzcfg(ix,n0+n,ncf) = molQxyz(ix,n,nm)
    enddo
  enddo
!
!  Set values for asymmetric unit if this is the first image
!
  nma = nmolf2a(nm)
  if (nmola2f(nma).eq.nm) then
    molQa(1:3,nma) = molQ(1:3,nm)
  endif
!
  return
  end
