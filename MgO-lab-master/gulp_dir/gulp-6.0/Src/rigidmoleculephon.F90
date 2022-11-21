  subroutine rigidmoleculephon(lcomplex)
!
!  Set second derivatives for rigid molecules during a phonon calculation
!  NB: Uses rotation angles rather than quaternions.
!
!   4/20 Created from rigidmoleculedrv
!   5/20 Corrections made to molQTdrv and updating local axes
!   5/20 molQQself removed as well as second derivatives in setrotation
!   5/20 axes moved to global array molaxes
!   5/20 Sign of complex term adjusted to ensure correct phasing
!   5/20 Parallel modifications added
!   6/20 nmolcore changes added
!   6/20 Phasing of molQQ terms corrected back
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
  use current
  use derivatives
  use molecule
  use parallel
  implicit none
!
!  Passed variables
!
  logical,                        intent(in)    :: lcomplex    ! If true then compute the imaginary part of second derivatives
!
!  Local variables
!
  integer(i4)                                   :: i
  integer(i4)                                   :: ii
#ifdef MPI
  integer(i4)                                   :: ind
#endif
  integer(i4)                                   :: ia
  integer(i4)                                   :: ix
  integer(i4)                                   :: iy
  integer(i4)                                   :: iz
  integer(i4)                                   :: j
  integer(i4)                                   :: j1
  integer(i4)                                   :: jx
  integer(i4)                                   :: jy
  integer(i4)                                   :: jz
  integer(i4)                                   :: k
  integer(i4)                                   :: k1
  integer(i4)                                   :: l
  integer(i4)                                   :: m
  integer(i4)                                   :: n
  integer(i4)                                   :: n3f
  integer(i4)                                   :: nm
  integer(i4)                                   :: nn
#ifdef MPI
  integer(i4)                                   :: nsize
  integer(i4)                                   :: status
#endif
  real(dp)                                      :: drR(3,3,3)
  real(dp)                                      :: drS(3,3,3)
  real(dp)                                      :: drRi(3,3)
  real(dp)                                      :: drRj(3,3)
#ifdef MPI
  real(dp),   dimension(:),   allocatable       :: sumin
  real(dp),   dimension(:),   allocatable       :: sumout
#endif
!
!  Shells are alredy process for phonons and so only handle cores
!
  if (ncorenomol.gt.0) then
    n3f = 3*ncorenomolptr(ncorenomol)
!
!  Initialise derivatives of molecules : Note that Q here will be used for R
!
    molQCdrv(1:n3f,1:3,1:nmol) = 0.0_dp
    molTCdrv(1:n3f,1:3,1:nmol) = 0.0_dp
  endif
  molQQdrv(1:3,1:3,1:nmol,1:nmol) = 0.0_dp
  molQTdrv(1:3,1:3,1:nmol,1:nmol) = 0.0_dp
  molTTdrv(1:3,1:3,1:nmol,1:nmol) = 0.0_dp
!
  if (lcomplex) then
    if (ncorenomol.gt.0) then
      molQCdri(1:n3f,1:3,1:nmol) = 0.0_dp
      molTCdri(1:n3f,1:3,1:nmol) = 0.0_dp
    endif
    molQQdri(1:3,1:3,1:nmol,1:nmol) = 0.0_dp
    molQTdri(1:3,1:3,1:nmol,1:nmol) = 0.0_dp
    molTTdri(1:3,1:3,1:nmol,1:nmol) = 0.0_dp
  endif
!
!  Sum of translational second derivatives
!
  if (nprocs.gt.1) then
!-----------------------
!  Parallel algorithm  |
!-----------------------
    do ii = 1,natomsonnode
      i = node2atom(ii)
      ix = 3*(ii - 1) + 1
      iy = 3*(ii - 1) + 2
      iz = 3*(ii - 1) + 3
!
!  Map atom back to molecule
!
      nm = natmol(i)
      if (nm.gt.0) then
        m = natinmol(i)
        if (m.le.nmolcore(nm)) then
!
!  Molecule - molecule
!
          do nn = 1,nmol
            do n = 1,nmolcore(nn)
              j = nmollist(nmolptr(nn)+n)
              jx = 3*(j - 1) + 1
              jy = 3*(j - 1) + 2
              jz = 3*(j - 1) + 3
              molTTdrv(1,1,nn,nm) = molTTdrv(1,1,nn,nm) + derv2(jx,ix)
              molTTdrv(2,1,nn,nm) = molTTdrv(2,1,nn,nm) + derv2(jy,ix)
              molTTdrv(3,1,nn,nm) = molTTdrv(3,1,nn,nm) + derv2(jz,ix)
              molTTdrv(1,2,nn,nm) = molTTdrv(1,2,nn,nm) + derv2(jx,iy)
              molTTdrv(2,2,nn,nm) = molTTdrv(2,2,nn,nm) + derv2(jy,iy)
              molTTdrv(3,2,nn,nm) = molTTdrv(3,2,nn,nm) + derv2(jz,iy)
              molTTdrv(1,3,nn,nm) = molTTdrv(1,3,nn,nm) + derv2(jx,iz)
              molTTdrv(2,3,nn,nm) = molTTdrv(2,3,nn,nm) + derv2(jy,iz)
              molTTdrv(3,3,nn,nm) = molTTdrv(3,3,nn,nm) + derv2(jz,iz)
            enddo
          enddo
!
!  Molecule - non-molecule
!
          do nn = 1,ncorenomol
            j = ncorenomolptr(nn)
            jx = 3*(j - 1) + 1
            jy = 3*(j - 1) + 2
            jz = 3*(j - 1) + 3
            molTCdrv(jx,1,nm) = molTCdrv(jx,1,nm) + derv2(jx,ix)
            molTCdrv(jy,1,nm) = molTCdrv(jy,1,nm) + derv2(jy,ix)
            molTCdrv(jz,1,nm) = molTCdrv(jz,1,nm) + derv2(jz,ix)
            molTCdrv(jx,2,nm) = molTCdrv(jx,2,nm) + derv2(jx,iy)
            molTCdrv(jy,2,nm) = molTCdrv(jy,2,nm) + derv2(jy,iy)
            molTCdrv(jz,2,nm) = molTCdrv(jz,2,nm) + derv2(jz,iy)
            molTCdrv(jx,3,nm) = molTCdrv(jx,3,nm) + derv2(jx,iz)
            molTCdrv(jy,3,nm) = molTCdrv(jy,3,nm) + derv2(jy,iz)
            molTCdrv(jz,3,nm) = molTCdrv(jz,3,nm) + derv2(jz,iz)
          enddo
        endif
      endif
    enddo
!
    if (lcomplex) then
      do ii = 1,natomsonnode
        i = node2atom(ii)
        ix = 3*(ii - 1) + 1 
        iy = 3*(ii - 1) + 2 
        iz = 3*(ii - 1) + 3 
!
!  Map atom back to molecule
!
        nm = natmol(i)
        if (nm.gt.0) then
          m = natinmol(i)
          if (m.le.nmolcore(nm)) then
!
!  Molecule - molecule
!
            do nn = 1,nmol
              do n = 1,nmolcore(nn)
                j = nmollist(nmolptr(nn)+n)
                jx = 3*(j - 1) + 1
                jy = 3*(j - 1) + 2
                jz = 3*(j - 1) + 3
!
                molTTdri(1,1,nn,nm) = molTTdri(1,1,nn,nm) + dervi(jx,ix)
                molTTdri(2,1,nn,nm) = molTTdri(2,1,nn,nm) + dervi(jy,ix)
                molTTdri(3,1,nn,nm) = molTTdri(3,1,nn,nm) + dervi(jz,ix)
                molTTdri(1,2,nn,nm) = molTTdri(1,2,nn,nm) + dervi(jx,iy)
                molTTdri(2,2,nn,nm) = molTTdri(2,2,nn,nm) + dervi(jy,iy)
                molTTdri(3,2,nn,nm) = molTTdri(3,2,nn,nm) + dervi(jz,iy)
                molTTdri(1,3,nn,nm) = molTTdri(1,3,nn,nm) + dervi(jx,iz)
                molTTdri(2,3,nn,nm) = molTTdri(2,3,nn,nm) + dervi(jy,iz)
                molTTdri(3,3,nn,nm) = molTTdri(3,3,nn,nm) + dervi(jz,iz)
              enddo
            enddo
!
!  Molecule - non-molecule
!
            do nn = 1,ncorenomol
              j = ncorenomolptr(nn)
              jx = 3*(j - 1) + 1
              jy = 3*(j - 1) + 2
              jz = 3*(j - 1) + 3
!
              molTCdri(jx,1,nm) = molTCdri(jx,1,nm) + dervi(jx,ix)
              molTCdri(jy,1,nm) = molTCdri(jy,1,nm) + dervi(jy,ix)
              molTCdri(jz,1,nm) = molTCdri(jz,1,nm) + dervi(jz,ix)
              molTCdri(jx,2,nm) = molTCdri(jx,2,nm) + dervi(jx,iy)
              molTCdri(jy,2,nm) = molTCdri(jy,2,nm) + dervi(jy,iy)
              molTCdri(jz,2,nm) = molTCdri(jz,2,nm) + dervi(jz,iy)
              molTCdri(jx,3,nm) = molTCdri(jx,3,nm) + dervi(jx,iz)
              molTCdri(jy,3,nm) = molTCdri(jy,3,nm) + dervi(jy,iz)
              molTCdri(jz,3,nm) = molTCdri(jz,3,nm) + dervi(jz,iz)
            enddo
          endif
        endif
      enddo
    endif
  else
!---------------------
!  Serial algorithm  |
!---------------------
    do nm = 1,nmol
      do m = 1,nmolcore(nm)
        i = nmollist(nmolptr(nm)+m)
        ix = 3*(i - 1) + 1
        iy = 3*(i - 1) + 2
        iz = 3*(i - 1) + 3
!
!  Molecule - molecule
!
        do nn = 1,nmol
          do n = 1,nmolcore(nn)
            j = nmollist(nmolptr(nn)+n)
            jx = 3*(j - 1) + 1
            jy = 3*(j - 1) + 2
            jz = 3*(j - 1) + 3
            molTTdrv(1,1,nn,nm) = molTTdrv(1,1,nn,nm) + derv2(jx,ix)
            molTTdrv(2,1,nn,nm) = molTTdrv(2,1,nn,nm) + derv2(jy,ix)
            molTTdrv(3,1,nn,nm) = molTTdrv(3,1,nn,nm) + derv2(jz,ix)
            molTTdrv(1,2,nn,nm) = molTTdrv(1,2,nn,nm) + derv2(jx,iy)
            molTTdrv(2,2,nn,nm) = molTTdrv(2,2,nn,nm) + derv2(jy,iy)
            molTTdrv(3,2,nn,nm) = molTTdrv(3,2,nn,nm) + derv2(jz,iy)
            molTTdrv(1,3,nn,nm) = molTTdrv(1,3,nn,nm) + derv2(jx,iz)
            molTTdrv(2,3,nn,nm) = molTTdrv(2,3,nn,nm) + derv2(jy,iz)
            molTTdrv(3,3,nn,nm) = molTTdrv(3,3,nn,nm) + derv2(jz,iz)
          enddo
        enddo
!
!  Molecule - non-molecule
!
        do nn = 1,ncorenomol
          j = ncorenomolptr(nn)
          jx = 3*(j - 1) + 1
          jy = 3*(j - 1) + 2
          jz = 3*(j - 1) + 3
          molTCdrv(jx,1,nm) = molTCdrv(jx,1,nm) + derv2(jx,ix)
          molTCdrv(jy,1,nm) = molTCdrv(jy,1,nm) + derv2(jy,ix)
          molTCdrv(jz,1,nm) = molTCdrv(jz,1,nm) + derv2(jz,ix)
          molTCdrv(jx,2,nm) = molTCdrv(jx,2,nm) + derv2(jx,iy)
          molTCdrv(jy,2,nm) = molTCdrv(jy,2,nm) + derv2(jy,iy)
          molTCdrv(jz,2,nm) = molTCdrv(jz,2,nm) + derv2(jz,iy)
          molTCdrv(jx,3,nm) = molTCdrv(jx,3,nm) + derv2(jx,iz)
          molTCdrv(jy,3,nm) = molTCdrv(jy,3,nm) + derv2(jy,iz)
          molTCdrv(jz,3,nm) = molTCdrv(jz,3,nm) + derv2(jz,iz)
        enddo
      enddo
    enddo
!
    if (lcomplex) then
      do nm = 1,nmol
        do m = 1,nmolcore(nm)
          i = nmollist(nmolptr(nm)+m)
          ix = 3*(i - 1) + 1
          iy = 3*(i - 1) + 2
          iz = 3*(i - 1) + 3
!
!  Molecule - molecule
!
          do nn = 1,nmol
            do n = 1,nmolcore(nn)
              j = nmollist(nmolptr(nn)+n)
              jx = 3*(j - 1) + 1
              jy = 3*(j - 1) + 2
              jz = 3*(j - 1) + 3
!
              molTTdri(1,1,nn,nm) = molTTdri(1,1,nn,nm) + dervi(jx,ix)
              molTTdri(2,1,nn,nm) = molTTdri(2,1,nn,nm) + dervi(jy,ix)
              molTTdri(3,1,nn,nm) = molTTdri(3,1,nn,nm) + dervi(jz,ix)
              molTTdri(1,2,nn,nm) = molTTdri(1,2,nn,nm) + dervi(jx,iy)
              molTTdri(2,2,nn,nm) = molTTdri(2,2,nn,nm) + dervi(jy,iy)
              molTTdri(3,2,nn,nm) = molTTdri(3,2,nn,nm) + dervi(jz,iy)
              molTTdri(1,3,nn,nm) = molTTdri(1,3,nn,nm) + dervi(jx,iz)
              molTTdri(2,3,nn,nm) = molTTdri(2,3,nn,nm) + dervi(jy,iz)
              molTTdri(3,3,nn,nm) = molTTdri(3,3,nn,nm) + dervi(jz,iz)
            enddo
          enddo
!
!  Molecule - non-molecule
!
          do nn = 1,ncorenomol
            j = ncorenomolptr(nn)
            jx = 3*(j - 1) + 1
            jy = 3*(j - 1) + 2
            jz = 3*(j - 1) + 3
!
            molTCdri(jx,1,nm) = molTCdri(jx,1,nm) + dervi(jx,ix)
            molTCdri(jy,1,nm) = molTCdri(jy,1,nm) + dervi(jy,ix)
            molTCdri(jz,1,nm) = molTCdri(jz,1,nm) + dervi(jz,ix)
            molTCdri(jx,2,nm) = molTCdri(jx,2,nm) + dervi(jx,iy)
            molTCdri(jy,2,nm) = molTCdri(jy,2,nm) + dervi(jy,iy)
            molTCdri(jz,2,nm) = molTCdri(jz,2,nm) + dervi(jz,iy)
            molTCdri(jx,3,nm) = molTCdri(jx,3,nm) + dervi(jx,iz)
            molTCdri(jy,3,nm) = molTCdri(jy,3,nm) + dervi(jy,iz)
            molTCdri(jz,3,nm) = molTCdri(jz,3,nm) + dervi(jz,iz)
          enddo
        enddo
      enddo
    endif
  endif
!
!  Loop over molecules to handle rotation derivatives
!
  if (nprocs.gt.1) then
!-----------------------
!  Parallel algorithm  |
!-----------------------
    do ii = 1,natomsonnode
      i = node2atom(ii)
      ix = 3*(ii - 1) + 1
      iy = 3*(ii - 1) + 2
      iz = 3*(ii - 1) + 3
!
!  Map atom back to molecule
!
      nm = natmol(i)
      if (nm.gt.0) then
        m = natinmol(i)
        if (m.le.nmolcore(nm)) then
!
!  Set up rotation matrix derivatives
!
          call setrotation(molaxes(1,1,nm),drR)
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
            jx = 3*(j - 1) + 1
            jy = 3*(j - 1) + 2
            jz = 3*(j - 1) + 3
            do k = 1,3
              molQCdrv(jx,k,nm) = molQCdrv(jx,k,nm) + derv2(jx,ix)*drRi(1,k)
              molQCdrv(jy,k,nm) = molQCdrv(jy,k,nm) + derv2(jy,ix)*drRi(1,k)
              molQCdrv(jz,k,nm) = molQCdrv(jz,k,nm) + derv2(jz,ix)*drRi(1,k)
              molQCdrv(jx,k,nm) = molQCdrv(jx,k,nm) + derv2(jx,iy)*drRi(2,k)
              molQCdrv(jy,k,nm) = molQCdrv(jy,k,nm) + derv2(jy,iy)*drRi(2,k)
              molQCdrv(jz,k,nm) = molQCdrv(jz,k,nm) + derv2(jz,iy)*drRi(2,k)
              molQCdrv(jx,k,nm) = molQCdrv(jx,k,nm) + derv2(jx,iz)*drRi(3,k)
              molQCdrv(jy,k,nm) = molQCdrv(jy,k,nm) + derv2(jy,iz)*drRi(3,k)
              molQCdrv(jz,k,nm) = molQCdrv(jz,k,nm) + derv2(jz,iz)*drRi(3,k)
            enddo
          enddo
!
          if (lcomplex) then
            do nn = 1,ncorenomol
              j = ncorenomolptr(nn)
              jx = 3*(j - 1) + 1
              jy = 3*(j - 1) + 2
              jz = 3*(j - 1) + 3
              do k = 1,3
                molQCdri(jx,k,nm) = molQCdri(jx,k,nm) + dervi(jx,ix)*drRi(1,k)
                molQCdri(jy,k,nm) = molQCdri(jy,k,nm) + dervi(jy,ix)*drRi(1,k)
                molQCdri(jz,k,nm) = molQCdri(jz,k,nm) + dervi(jz,ix)*drRi(1,k)
                molQCdri(jx,k,nm) = molQCdri(jx,k,nm) + dervi(jx,iy)*drRi(2,k)
                molQCdri(jy,k,nm) = molQCdri(jy,k,nm) + dervi(jy,iy)*drRi(2,k)
                molQCdri(jz,k,nm) = molQCdri(jz,k,nm) + dervi(jz,iy)*drRi(2,k)
                molQCdri(jx,k,nm) = molQCdri(jx,k,nm) + dervi(jx,iz)*drRi(3,k)
                molQCdri(jy,k,nm) = molQCdri(jy,k,nm) + dervi(jy,iz)*drRi(3,k)
                molQCdri(jz,k,nm) = molQCdri(jz,k,nm) + dervi(jz,iz)*drRi(3,k)
              enddo
            enddo
          endif
!
!  Loop over second molecule
!
          do nn = 1,nmol
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
              jx = 3*(j - 1)
              do ia = 1,3
                do k = 1,3
                  molQTdrv(k,ia,nm,nn) = molQTdrv(k,ia,nm,nn) + drRi(1,k)*derv2(jx+ia,ix) + &
                                                                drRi(2,k)*derv2(jx+ia,iy) + &
                                                                drRi(3,k)*derv2(jx+ia,iz)
                enddo
              enddo
!
              if (lcomplex) then
                do ia = 1,3
                  do k = 1,3
                    molQTdri(k,ia,nm,nn) = molQTdri(k,ia,nm,nn) + drRi(1,k)*dervi(jx+ia,ix) + &
                                                                  drRi(2,k)*dervi(jx+ia,iy) + &
                                                                  drRi(3,k)*dervi(jx+ia,iz)
                  enddo
                enddo
              endif
!************************************
!  Rotation - rotation derivatives  *
!************************************
              jx = 3*(j - 1) + 1
              jy = 3*(j - 1) + 2
              jz = 3*(j - 1) + 3
!
!  Term 1
!
              do j1 = 1,3
                do k1 = 1,3
                  molQQdrv(k1,j1,nn,nm) = molQQdrv(k1,j1,nn,nm) + derv2(jx,ix)*drRj(1,k1)*drRi(1,j1) &
                                                                + derv2(jy,ix)*drRj(2,k1)*drRi(1,j1) &
                                                                + derv2(jz,ix)*drRj(3,k1)*drRi(1,j1) &
                                                                + derv2(jx,iy)*drRj(1,k1)*drRi(2,j1) &
                                                                + derv2(jy,iy)*drRj(2,k1)*drRi(2,j1) &
                                                                + derv2(jz,iy)*drRj(3,k1)*drRi(2,j1) &
                                                                + derv2(jx,iz)*drRj(1,k1)*drRi(3,j1) &
                                                                + derv2(jy,iz)*drRj(2,k1)*drRi(3,j1) &
                                                                + derv2(jz,iz)*drRj(3,k1)*drRi(3,j1)
                enddo
              enddo
!
              if (lcomplex) then
                do j1 = 1,3
                  do k1 = 1,3
                    molQQdri(k1,j1,nn,nm) = molQQdri(k1,j1,nn,nm) + dervi(jx,ix)*drRj(1,k1)*drRi(1,j1) &
                                                                  + dervi(jy,ix)*drRj(2,k1)*drRi(1,j1) &
                                                                  + dervi(jz,ix)*drRj(3,k1)*drRi(1,j1) &
                                                                  + dervi(jx,iy)*drRj(1,k1)*drRi(2,j1) &
                                                                  + dervi(jy,iy)*drRj(2,k1)*drRi(2,j1) &
                                                                  + dervi(jz,iy)*drRj(3,k1)*drRi(2,j1) &
                                                                  + dervi(jx,iz)*drRj(1,k1)*drRi(3,j1) &
                                                                  + dervi(jy,iz)*drRj(2,k1)*drRi(3,j1) &
                                                                  + dervi(jz,iz)*drRj(3,k1)*drRi(3,j1)
                  enddo
                enddo
              endif
            enddo
          enddo
        endif
      endif
    enddo
  else
!---------------------
!  Serial algorithm  |
!---------------------
    do nm = 1,nmol
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
        ix = 3*(i - 1) + 1
        iy = 3*(i - 1) + 2
        iz = 3*(i - 1) + 3
!
!  Loop over non-molecule atoms
!
        do nn = 1,ncorenomol
          j = ncorenomolptr(nn)
          jx = 3*(j - 1) + 1
          jy = 3*(j - 1) + 2
          jz = 3*(j - 1) + 3
          do k = 1,3
            molQCdrv(jx,k,nm) = molQCdrv(jx,k,nm) + derv2(jx,ix)*drRi(1,k)
            molQCdrv(jy,k,nm) = molQCdrv(jy,k,nm) + derv2(jy,ix)*drRi(1,k)
            molQCdrv(jz,k,nm) = molQCdrv(jz,k,nm) + derv2(jz,ix)*drRi(1,k)
            molQCdrv(jx,k,nm) = molQCdrv(jx,k,nm) + derv2(jx,iy)*drRi(2,k)
            molQCdrv(jy,k,nm) = molQCdrv(jy,k,nm) + derv2(jy,iy)*drRi(2,k)
            molQCdrv(jz,k,nm) = molQCdrv(jz,k,nm) + derv2(jz,iy)*drRi(2,k)
            molQCdrv(jx,k,nm) = molQCdrv(jx,k,nm) + derv2(jx,iz)*drRi(3,k)
            molQCdrv(jy,k,nm) = molQCdrv(jy,k,nm) + derv2(jy,iz)*drRi(3,k)
            molQCdrv(jz,k,nm) = molQCdrv(jz,k,nm) + derv2(jz,iz)*drRi(3,k)
          enddo
        enddo
!
        if (lcomplex) then
          do nn = 1,ncorenomol
            j = ncorenomolptr(nn)
            jx = 3*(j - 1) + 1
            jy = 3*(j - 1) + 2
            jz = 3*(j - 1) + 3
            do k = 1,3
              molQCdri(jx,k,nm) = molQCdri(jx,k,nm) + dervi(jx,ix)*drRi(1,k)
              molQCdri(jy,k,nm) = molQCdri(jy,k,nm) + dervi(jy,ix)*drRi(1,k)
              molQCdri(jz,k,nm) = molQCdri(jz,k,nm) + dervi(jz,ix)*drRi(1,k)
              molQCdri(jx,k,nm) = molQCdri(jx,k,nm) + dervi(jx,iy)*drRi(2,k)
              molQCdri(jy,k,nm) = molQCdri(jy,k,nm) + dervi(jy,iy)*drRi(2,k)
              molQCdri(jz,k,nm) = molQCdri(jz,k,nm) + dervi(jz,iy)*drRi(2,k)
              molQCdri(jx,k,nm) = molQCdri(jx,k,nm) + dervi(jx,iz)*drRi(3,k)
              molQCdri(jy,k,nm) = molQCdri(jy,k,nm) + dervi(jy,iz)*drRi(3,k)
              molQCdri(jz,k,nm) = molQCdri(jz,k,nm) + dervi(jz,iz)*drRi(3,k)
            enddo
          enddo
        endif
!
!  Loop over second molecule
!
        do nn = 1,nmol
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
            jx = 3*(j - 1)
            do ia = 1,3
              do k = 1,3
                molQTdrv(k,ia,nm,nn) = molQTdrv(k,ia,nm,nn) + (drRi(1,k)*derv2(jx+ia,ix) + &
                                                              drRi(2,k)*derv2(jx+ia,iy) + &
                                                              drRi(3,k)*derv2(jx+ia,iz))
              enddo
            enddo
!
            if (lcomplex) then
              do ia = 1,3
                do k = 1,3
                  molQTdri(k,ia,nm,nn) = molQTdri(k,ia,nm,nn) + (drRi(1,k)*dervi(jx+ia,ix) + &
                                                                drRi(2,k)*dervi(jx+ia,iy) + &
                                                                drRi(3,k)*dervi(jx+ia,iz))
                enddo
              enddo
            endif
!************************************
!  Rotation - rotation derivatives  *
!************************************
            jx = 3*(j - 1) + 1
            jy = 3*(j - 1) + 2
            jz = 3*(j - 1) + 3
!
!  Term 1
!
            do j1 = 1,3
              do k1 = 1,3
                molQQdrv(k1,j1,nn,nm) = molQQdrv(k1,j1,nn,nm) + derv2(jx,ix)*drRj(1,k1)*drRi(1,j1) &
                                                              + derv2(jy,ix)*drRj(2,k1)*drRi(1,j1) &
                                                              + derv2(jz,ix)*drRj(3,k1)*drRi(1,j1) &
                                                              + derv2(jx,iy)*drRj(1,k1)*drRi(2,j1) &
                                                              + derv2(jy,iy)*drRj(2,k1)*drRi(2,j1) &
                                                              + derv2(jz,iy)*drRj(3,k1)*drRi(2,j1) &
                                                              + derv2(jx,iz)*drRj(1,k1)*drRi(3,j1) &
                                                              + derv2(jy,iz)*drRj(2,k1)*drRi(3,j1) &
                                                              + derv2(jz,iz)*drRj(3,k1)*drRi(3,j1)
              enddo
            enddo
!
            if (lcomplex) then
              do j1 = 1,3
                do k1 = 1,3
                  molQQdri(k1,j1,nn,nm) = molQQdri(k1,j1,nn,nm) + dervi(jx,ix)*drRj(1,k1)*drRi(1,j1) &
                                                                + dervi(jy,ix)*drRj(2,k1)*drRi(1,j1) &
                                                                + dervi(jz,ix)*drRj(3,k1)*drRi(1,j1) &
                                                                + dervi(jx,iy)*drRj(1,k1)*drRi(2,j1) &
                                                                + dervi(jy,iy)*drRj(2,k1)*drRi(2,j1) &
                                                                + dervi(jz,iy)*drRj(3,k1)*drRi(2,j1) &
                                                                + dervi(jx,iz)*drRj(1,k1)*drRi(3,j1) &
                                                                + dervi(jy,iz)*drRj(2,k1)*drRi(3,j1) &
                                                                + dervi(jz,iz)*drRj(3,k1)*drRi(3,j1)
                enddo
              enddo
            endif
          enddo
        enddo
      enddo
    enddo
  endif
#ifdef MPI
!
!  Globalise molecule derivatives
!
  if (nprocs.gt.1) then
    if (ncorenomol.gt.0) then
      nsize = max(27*nmol*nmol,3*nmol*n3f)
    else
      nsize = 27*nmol*nmol
    endif
    allocate(sumin(nsize),stat=status)
    if (status/=0) call outofmemory('rigidmoleculephon','sumin')
    allocate(sumout(nsize),stat=status)
    if (status/=0) call outofmemory('rigidmoleculephon','sumout')
    if (ncorenomol.gt.0) then
!
!  molQCdrv
!
      ind = 0
      do nm = 1,nmol
        do ix = 1,3
          sumin(ind+1:ind+n3f) = molQCdrv(1:n3f,ix,nm)
          ind = ind + n3f
        enddo
      enddo
      call sumall(sumin,sumout,ind,"rigidmoleculephon","molQCdrv")
      ind = 0
      do nm = 1,nmol
        do ix = 1,3
          molQCdrv(1:n3f,ix,nm) = sumout(ind+1:ind+n3f)
          ind = ind + n3f
        enddo
      enddo
!
!  molTCdrv
!
      ind = 0
      do nm = 1,nmol
        do ix = 1,3
          sumin(ind+1:ind+n3f) = molTCdrv(1:n3f,ix,nm)
          ind = ind + n3f
        enddo
      enddo
      call sumall(sumin,sumout,ind,"rigidmoleculephon","molTCdrv")
      ind = 0
      do nm = 1,nmol
        do ix = 1,3
          molTCdrv(1:n3f,ix,nm) = sumout(ind+1:ind+n3f)
          ind = ind + n3f
        enddo
      enddo
    endif
!
!  molQQdrv, molQTdrv, molTTdrv
!
    ind = 0
    do nm = 1,nmol
      do nn = 1,nmol
        do ix = 1,3
          do jx = 1,3
            sumin(ind+1) = molQQdrv(jx,ix,nn,nm)
            sumin(ind+2) = molQTdrv(jx,ix,nn,nm)
            sumin(ind+3) = molTTdrv(jx,ix,nn,nm)
            ind = ind + 3
          enddo
        enddo
      enddo
    enddo
    call sumall(sumin,sumout,ind,"rigidmoleculephon","molQQdrv")
    ind = 0
    do nm = 1,nmol
      do nn = 1,nmol
        do ix = 1,3
          do jx = 1,3
            molQQdrv(jx,ix,nn,nm) = sumout(ind+1)
            molQTdrv(jx,ix,nn,nm) = sumout(ind+2)
            molTTdrv(jx,ix,nn,nm) = sumout(ind+3)
            ind = ind + 3
          enddo
        enddo
      enddo
    enddo
!
    if (lcomplex) then
      if (ncorenomol.gt.0) then
!
!  molQCdri
!
        ind = 0
        do nm = 1,nmol
          do ix = 1,3
            sumin(ind+1:ind+n3f) = molQCdri(1:n3f,ix,nm)
            ind = ind + n3f
          enddo
        enddo
        call sumall(sumin,sumout,ind,"rigidmoleculephon","molQCdri")
        ind = 0
        do nm = 1,nmol
          do ix = 1,3
            molQCdri(1:n3f,ix,nm) = sumout(ind+1:ind+n3f)
            ind = ind + n3f
          enddo
        enddo
!
!  molTCdri
!
        ind = 0
        do nm = 1,nmol
          do ix = 1,3
            sumin(ind+1:ind+n3f) = molTCdri(1:n3f,ix,nm)
            ind = ind + n3f
          enddo
        enddo
        call sumall(sumin,sumout,ind,"rigidmoleculephon","molTCdri")
        ind = 0
        do nm = 1,nmol
          do ix = 1,3
            molTCdri(1:n3f,ix,nm) = sumout(ind+1:ind+n3f)
            ind = ind + n3f
          enddo
        enddo
      endif
!
!  molQQdri, molQTdri, molTTdri
!
      ind = 0
      do nm = 1,nmol
        do nn = 1,nmol
          do ix = 1,3
            do jx = 1,3
              sumin(ind+1) = molQQdri(jx,ix,nn,nm)
              sumin(ind+2) = molQTdri(jx,ix,nn,nm)
              sumin(ind+3) = molTTdri(jx,ix,nn,nm)
              ind = ind + 3
            enddo
          enddo
        enddo
      enddo
      call sumall(sumin,sumout,ind,"rigidmoleculephon","molQQdri")
      ind = 0
      do nm = 1,nmol
        do nn = 1,nmol
          do ix = 1,3
            do jx = 1,3
              molQQdri(jx,ix,nn,nm) = sumout(ind+1)
              molQTdri(jx,ix,nn,nm) = sumout(ind+2)
              molTTdri(jx,ix,nn,nm) = sumout(ind+3)
              ind = ind + 3
            enddo
          enddo
        enddo
      enddo
    endif
!
    deallocate(sumout,stat=status)
    if (status/=0) call deallocate_error('rigidmoleculephon','sumout')
    deallocate(sumin,stat=status)
    if (status/=0) call deallocate_error('rigidmoleculephon','sumin')
  endif
#endif
!
  return
  end
