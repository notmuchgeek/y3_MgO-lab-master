  subroutine rigidmoleculedrv(lgrad2)
!
!  Set derivatives for rigid molecule translations
!
!   2/19 Created
!  10/19 Quaternion derivatives added
!  12/19 Second derivatives added
!   3/20 Molecule number added to setquaternion arguments
!   4/20 Modified so that derv3c is used in the finite strain case
!   4/20 Transformation of derv3 to correct strain-Cartesian mixed second derivatives
!        added for finite strain case
!   4/20 Correction to molQQdrv
!   4/20 Storing of molQQdrv term 2 added for use in phonons
!   5/20 molQQself removed
!   5/20 Correction to molQTdrv
!   5/20 References to nasym removed since rigid molecule derivatives don't currently
!        use symmetry for second derivatives
!   5/20 Parallel modifications added
!   6/20 nmolcore changes added
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
  use derivatives
  use m_strain,       only : strainfull, straininverse
  use molecule
  use parallel
  use symmetry,       only : lstr, lsymopt
  implicit none
!
!  Passed variables
!
  logical,                      intent(in)    :: lgrad2
!
!  Local variables
!
  integer(i4)                                 :: i
  integer(i4)                                 :: ia
#ifdef MPI
  integer(i4)                                 :: ind
#endif
  integer(i4)                                 :: ii
  integer(i4)                                 :: ix
  integer(i4)                                 :: iy
  integer(i4)                                 :: iz
  integer(i4)                                 :: j
  integer(i4)                                 :: j1
  integer(i4)                                 :: jx
  integer(i4)                                 :: jy
  integer(i4)                                 :: jz
  integer(i4)                                 :: k
  integer(i4)                                 :: k1
  integer(i4)                                 :: l
  integer(i4)                                 :: m
  integer(i4)                                 :: n
  integer(i4)                                 :: n3f
  integer(i4)                                 :: nm
  integer(i4)                                 :: nma
  integer(i4)                                 :: nn
  integer(i4)                                 :: nna
  integer(i4)                                 :: ns
#ifdef MPI
  integer(i4)                                 :: nsize
  integer(i4)                                 :: status
#endif
  real(dp)                                    :: derv3l(3,6)
  real(dp)                                    :: rQ(3,3)
  real(dp)                                    :: rQm2(3,3)
  real(dp)                                    :: drQ(3,3,3)
  real(dp)                                    :: drQm2(3,3,3)
  real(dp)                                    :: drQ2(3,3,3,3)
  real(dp)                                    :: drQ2m2(3,3,3,3)
  real(dp)                                    :: drQi(3,3)
  real(dp)                                    :: drQj(3,3)
  real(dp)                                    :: drQ2i(3,3,3)
  real(dp)                                    :: ricom(3)
  real(dp)                                    :: rjcom(3)
#ifdef MPI
  real(dp), dimension(:),   allocatable       :: sumin
  real(dp), dimension(:),   allocatable       :: sumout
#endif
  real(dp)                                    :: tmp(3)
!
  n3f = 3*numat
!
!  Initialise derivatives of molecules
!
  molTdrv(1:3,1:nmol) = 0.0_dp
  molQdrv(1:3,1:nmol) = 0.0_dp
  if (lgrad2) then
    molQCdrv(1:n3f,1:3,1:nmol) = 0.0_dp
    molTCdrv(1:n3f,1:3,1:nmol) = 0.0_dp
    molQQdrv(1:3,1:3,1:nmol,1:nmol) = 0.0_dp
    molQTdrv(1:3,1:3,1:nmol,1:nmol) = 0.0_dp
    molTTdrv(1:3,1:3,1:nmol,1:nmol) = 0.0_dp
    if (lstr) then
      molQSdrv(1:nstrains,1:3,1:nmol) = 0.0_dp
      molTSdrv(1:nstrains,1:3,1:nmol) = 0.0_dp
    endif
  endif
!
!  Compute translational derivatives of molecules by summing those of component atoms
!
  do nm = 1,nmol
    do m = 1,nmolcore(nm)
      i = nmollist(nmolptr(nm)+m)
      molTdrv(1,nm) = molTdrv(1,nm) + xdrv(i)
      molTdrv(2,nm) = molTdrv(2,nm) + ydrv(i)
      molTdrv(3,nm) = molTdrv(3,nm) + zdrv(i)
    enddo
  enddo
!
!  Sum of translational second derivatives
!
  if (lgrad2) then
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
!
!  Molecule - molecule
!
          if (m.le.nmolcore(nm)) then
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
            do nn = 1,numatnomol
              j = numatnomolptr(nn)
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
!
!  Molecule strain
!
            if (lstr) then
              do n = 1,nstrains
                molTSdrv(n,1,nm) = molTSdrv(n,1,nm) + derv3(ix,n)
                molTSdrv(n,2,nm) = molTSdrv(n,2,nm) + derv3(iy,n)
                molTSdrv(n,3,nm) = molTSdrv(n,3,nm) + derv3(iz,n)
              enddo
            endif
          endif
        endif
      enddo
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
          do nn = 1,numatnomol
            j = numatnomolptr(nn)
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
!
!  Molecule strain
!
          if (lstr) then
            do n = 1,nstrains
              molTSdrv(n,1,nm) = molTSdrv(n,1,nm) + derv3(ix,n)
              molTSdrv(n,2,nm) = molTSdrv(n,2,nm) + derv3(iy,n)
              molTSdrv(n,3,nm) = molTSdrv(n,3,nm) + derv3(iz,n)
            enddo
          endif
        enddo
      enddo
    endif
  endif
  if (lgrad2.and.nprocs.gt.1) then
!-----------------------
!  Parallel algorithm  |
!-----------------------
!
!  Loop over molecules to handle quaternion first derivatives on all nodes
!
    do nm = 1,nmol
!
!  Set up quaternion matrix
!
      call setquaternion(nm,molQ(1,nm),rQ,drQ,drQ2,.true.,lgrad2)
!
      do m = 1,nmolcore(nm)
        i = nmollist(nmolptr(nm)+m)
!
!  Set original intra molecular vector
!
        if (lsymopt) then
          nma = nmolf2a(nm)
          ricom(1:3) = molQxyz(1:3,m,nma)
        else
          ricom(1:3) = molQxyz(1:3,m,nm)
        endif
!
!  Transform by derivatives of quaternion rotation matrix
!
        drQi(1:3,1:3) = 0.0_dp
        do l = 1,3
          do k = 1,3
            drQi(k,1) = drQi(k,1) + drQ(k,1,l)*ricom(l)
            drQi(k,2) = drQi(k,2) + drQ(k,2,l)*ricom(l)
            drQi(k,3) = drQi(k,3) + drQ(k,3,l)*ricom(l)
          enddo
        enddo
!
        do k = 1,3
          molQdrv(k,nm) = molQdrv(k,nm) + drQi(k,1)*xdrv(i) + drQi(k,2)*ydrv(i) + drQi(k,3)*zdrv(i)
        enddo
      enddo
    enddo
!
!  Now for the second derivatives...
!
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
!  Set up quaternion matrix
!
          call setquaternion(nm,molQ(1,nm),rQ,drQ,drQ2,.true.,lgrad2)
!
!  Set original intra molecular vector
!
          if (lsymopt) then
            nma = nmolf2a(nm)
            ricom(1:3) = molQxyz(1:3,m,nma)
          else
            ricom(1:3) = molQxyz(1:3,m,nm)
          endif
!
!  Transform by derivatives of quaternion rotation matrix
!
          drQi(1:3,1:3) = 0.0_dp
          do l = 1,3
            do k = 1,3
              drQi(k,1) = drQi(k,1) + drQ(k,1,l)*ricom(l)
              drQi(k,2) = drQi(k,2) + drQ(k,2,l)*ricom(l)
              drQi(k,3) = drQi(k,3) + drQ(k,3,l)*ricom(l)
            enddo
          enddo
!
!  Loop over non-molecule atoms
!
          do nn = 1,numatnomol
            j = numatnomolptr(nn)
            jx = 3*(j - 1) + 1
            jy = 3*(j - 1) + 2
            jz = 3*(j - 1) + 3
            do k = 1,3
              molQCdrv(jx,k,nm) = molQCdrv(jx,k,nm) + derv2(jx,ix)*drQi(k,1)
              molQCdrv(jy,k,nm) = molQCdrv(jy,k,nm) + derv2(jy,ix)*drQi(k,1)
              molQCdrv(jz,k,nm) = molQCdrv(jz,k,nm) + derv2(jz,ix)*drQi(k,1)
              molQCdrv(jx,k,nm) = molQCdrv(jx,k,nm) + derv2(jx,iy)*drQi(k,2)
              molQCdrv(jy,k,nm) = molQCdrv(jy,k,nm) + derv2(jy,iy)*drQi(k,2)
              molQCdrv(jz,k,nm) = molQCdrv(jz,k,nm) + derv2(jz,iy)*drQi(k,2)
              molQCdrv(jx,k,nm) = molQCdrv(jx,k,nm) + derv2(jx,iz)*drQi(k,3)
              molQCdrv(jy,k,nm) = molQCdrv(jy,k,nm) + derv2(jy,iz)*drQi(k,3)
              molQCdrv(jz,k,nm) = molQCdrv(jz,k,nm) + derv2(jz,iz)*drQi(k,3)
            enddo
          enddo
!
!  Loop over second molecule
!
          do nn = 1,nmol
!
!  Set up quaternion matrix for second molecule
!
            call setquaternion(nn,molQ(1,nn),rQm2,drQm2,drQ2m2,.true.,.false.)
!
            do n = 1,nmolcore(nn)
              j = nmollist(nmolptr(nn)+n)
!**********************************************
!  Mixed quaternion-translation derivatives  *
!**********************************************
              jx = 3*(j - 1)
              do ia = 1,3
                do k = 1,3
                  molQTdrv(k,ia,nm,nn) = molQTdrv(k,ia,nm,nn) + drQi(k,1)*derv2(jx+ia,ix) + &
                                                                drQi(k,2)*derv2(jx+ia,iy) + &
                                                                drQi(k,3)*derv2(jx+ia,iz)
                enddo
              enddo
!**************************************
!  Quaternion-quaternion derivatives  *
!**************************************
              jx = 3*(j - 1) + 1
              jy = 3*(j - 1) + 2
              jz = 3*(j - 1) + 3
!
!  Set original intra molecular vector
!
              if (lsymopt) then
                nna = nmolf2a(nn)
                rjcom(1:3) = molQxyz(1:3,n,nna)
              else
                rjcom(1:3) = molQxyz(1:3,n,nn)
              endif
!
!  Transform by derivatives of quaternion rotation matrix
!
              drQj(1:3,1:3) = 0.0_dp
              do l = 1,3
                do k = 1,3
                  drQj(k,1) = drQj(k,1) + drQm2(k,1,l)*rjcom(l)
                  drQj(k,2) = drQj(k,2) + drQm2(k,2,l)*rjcom(l)
                  drQj(k,3) = drQj(k,3) + drQm2(k,3,l)*rjcom(l)
                enddo
              enddo
!
              do j1 = 1,3
                do k1 = 1,3
!
!  Term 1
!
                  molQQdrv(k1,j1,nn,nm) = molQQdrv(k1,j1,nn,nm) + derv2(jx,ix)*drQj(k1,1)*drQi(j1,1) &
                                                                + derv2(jy,ix)*drQj(k1,2)*drQi(j1,1) &
                                                                + derv2(jz,ix)*drQj(k1,3)*drQi(j1,1) &
                                                                + derv2(jx,iy)*drQj(k1,1)*drQi(j1,2) &
                                                                + derv2(jy,iy)*drQj(k1,2)*drQi(j1,2) &
                                                                + derv2(jz,iy)*drQj(k1,3)*drQi(j1,2) &
                                                                + derv2(jx,iz)*drQj(k1,1)*drQi(j1,3) &
                                                                + derv2(jy,iz)*drQj(k1,2)*drQi(j1,3) &
                                                                + derv2(jz,iz)*drQj(k1,3)*drQi(j1,3)
                enddo
              enddo
            enddo
          enddo
!
!  Transform by second derivatives of quaternion rotation matrix
!
          drQ2i(1:3,1:3,1:3) = 0.0_dp
          do l = 1,3
            do j = 1,3
              do k = 1,3
                drQ2i(k,j,1) = drQ2i(k,j,1) + drQ2(k,j,1,l)*ricom(l)
                drQ2i(k,j,2) = drQ2i(k,j,2) + drQ2(k,j,2,l)*ricom(l)
                drQ2i(k,j,3) = drQ2i(k,j,3) + drQ2(k,j,3,l)*ricom(l)
              enddo
            enddo
          enddo
!
!  Term 2 - self only
!
          do j = 1,3
            do k = 1,3
              molQQdrv(k,j,nm,nm) = molQQdrv(k,j,nm,nm) + drQ2i(k,j,1)*xdrv(i) &
                                                        + drQ2i(k,j,2)*ydrv(i) &
                                                        + drQ2i(k,j,3)*zdrv(i)
            enddo
          enddo
!
!  Strain mixed derivatives
!
          if (lstr) then
!
!  Remove cell volume factor from the mixed strain derivatives
!
            if (lfinitestrain) then
              derv3l(1,1:nstrains) = derv3(ix,1:nstrains)
              derv3l(2,1:nstrains) = derv3(iy,1:nstrains)
              derv3l(3,1:nstrains) = derv3(iz,1:nstrains)
!
!  Transform mixed derivatives by strain matrix
!
              do ns = 1,nstrains
                tmp(1:ndim) = 0.0_dp
                do j = 1,ndim
                  do l = 1,ndim
                    tmp(j) = tmp(j) + derv3l(l,ns)*strainfull(l,j)
                  enddo
                enddo
                derv3l(1:ndim,ns) = tmp(1:ndim)
              enddo
!
!  Remove strain derivatives
!
              if (ndim.eq.3) then
                derv3l(1,1) =  derv3l(1,1) - xdrv(i)
                derv3l(2,2) =  derv3l(2,2) - ydrv(i)
                derv3l(3,3) =  derv3l(3,3) - zdrv(i)
                derv3l(1,5) =  derv3l(1,5) - 0.5_dp*zdrv(i)
                derv3l(1,6) =  derv3l(1,6) - 0.5_dp*ydrv(i)
                derv3l(2,4) =  derv3l(2,4) - 0.5_dp*zdrv(i)
                derv3l(2,6) =  derv3l(2,6) - 0.5_dp*xdrv(i)
                derv3l(3,4) =  derv3l(3,4) - 0.5_dp*ydrv(i)
                derv3l(3,5) =  derv3l(3,5) - 0.5_dp*xdrv(i)
              elseif (ndim.eq.2) then
                derv3l(1,1) =  derv3l(1,1) - xdrv(i)
                derv3l(2,2) =  derv3l(2,2) - ydrv(i)
                derv3l(1,3) =  derv3l(1,3) - 0.5_dp*ydrv(i)
                derv3l(2,3) =  derv3l(2,3) - 0.5_dp*xdrv(i)
              elseif (ndim.eq.1) then
                derv3l(1,1) =  derv3l(1,1) - xdrv(i)
              endif
!
!  Transform mixed derivatives by inverse strain matrix
!
              do ns = 1,nstrains
                tmp(1:ndim) = 0.0_dp
                do j = 1,ndim
                  do l = 1,ndim
                    tmp(j) = tmp(j) + derv3l(l,ns)*straininverse(l,j)
                  enddo
                enddo
                derv3l(1:ndim,ns) = tmp(1:ndim)
              enddo
            else
              derv3l(1,1:nstrains) = derv3(ix,1:nstrains)
              derv3l(2,1:nstrains) = derv3(iy,1:nstrains)
              derv3l(3,1:nstrains) = derv3(iz,1:nstrains)
!
              if (ndim.eq.3) then
                derv3l(1,1) =  derv3(ix,1) - xdrv(i)
                derv3l(2,2) =  derv3(iy,2) - ydrv(i)
                derv3l(3,3) =  derv3(iz,3) - zdrv(i)
                derv3l(1,5) =  derv3(ix,5) - 0.5_dp*zdrv(i)
                derv3l(1,6) =  derv3(ix,6) - 0.5_dp*ydrv(i)
                derv3l(2,4) =  derv3(iy,4) - 0.5_dp*zdrv(i)
                derv3l(2,6) =  derv3(iy,6) - 0.5_dp*xdrv(i)
                derv3l(3,4) =  derv3(iz,4) - 0.5_dp*ydrv(i)
                derv3l(3,5) =  derv3(iz,5) - 0.5_dp*xdrv(i)
              elseif (ndim.eq.2) then
                derv3l(1,1) =  derv3(ix,1) - xdrv(i)
                derv3l(2,2) =  derv3(iy,2) - ydrv(i)
                derv3l(1,3) =  derv3(ix,3) - 0.5_dp*ydrv(i)
                derv3l(2,3) =  derv3(iy,3) - 0.5_dp*xdrv(i)
              elseif (ndim.eq.1) then
                derv3l(1,1) =  derv3(ix,1) - xdrv(i)
              endif
            endif
!
            do n = 1,nstrains
              do k = 1,3
                molQSdrv(n,k,nm) = molQSdrv(n,k,nm) + derv3l(1,n)*drQi(k,1)
                molQSdrv(n,k,nm) = molQSdrv(n,k,nm) + derv3l(2,n)*drQi(k,2)
                molQSdrv(n,k,nm) = molQSdrv(n,k,nm) + derv3l(3,n)*drQi(k,3)
              enddo
            enddo
          endif
        endif
      endif
    enddo
  else
!---------------------
!  Serial algorithm  |
!---------------------
!
!  Loop over molecules to handle quaternion derivatives
!
    do nm = 1,nmol
!
!  Set up quaternion matrix
!
      call setquaternion(nm,molQ(1,nm),rQ,drQ,drQ2,.true.,lgrad2)
!
      do m = 1,nmolcore(nm)
        i = nmollist(nmolptr(nm)+m)
!
!  Set original intra molecular vector
!
        if (lsymopt) then
          nma = nmolf2a(nm)
          ricom(1:3) = molQxyz(1:3,m,nma)
        else
          ricom(1:3) = molQxyz(1:3,m,nm)
        endif
!
!  Transform by derivatives of quaternion rotation matrix
!
        drQi(1:3,1:3) = 0.0_dp
        do l = 1,3
          do k = 1,3
            drQi(k,1) = drQi(k,1) + drQ(k,1,l)*ricom(l)
            drQi(k,2) = drQi(k,2) + drQ(k,2,l)*ricom(l)
            drQi(k,3) = drQi(k,3) + drQ(k,3,l)*ricom(l)
          enddo
        enddo
!
        do k = 1,3
          molQdrv(k,nm) = molQdrv(k,nm) + drQi(k,1)*xdrv(i) + drQi(k,2)*ydrv(i) + drQi(k,3)*zdrv(i)
        enddo
!
        if (lgrad2) then
          ix = 3*(i - 1) + 1
          iy = 3*(i - 1) + 2
          iz = 3*(i - 1) + 3
!
!  Loop over non-molecule atoms
!
          do nn = 1,numatnomol
            j = numatnomolptr(nn)
            jx = 3*(j - 1) + 1
            jy = 3*(j - 1) + 2
            jz = 3*(j - 1) + 3
            do k = 1,3
              molQCdrv(jx,k,nm) = molQCdrv(jx,k,nm) + derv2(jx,ix)*drQi(k,1)
              molQCdrv(jy,k,nm) = molQCdrv(jy,k,nm) + derv2(jy,ix)*drQi(k,1)
              molQCdrv(jz,k,nm) = molQCdrv(jz,k,nm) + derv2(jz,ix)*drQi(k,1)
              molQCdrv(jx,k,nm) = molQCdrv(jx,k,nm) + derv2(jx,iy)*drQi(k,2)
              molQCdrv(jy,k,nm) = molQCdrv(jy,k,nm) + derv2(jy,iy)*drQi(k,2)
              molQCdrv(jz,k,nm) = molQCdrv(jz,k,nm) + derv2(jz,iy)*drQi(k,2)
              molQCdrv(jx,k,nm) = molQCdrv(jx,k,nm) + derv2(jx,iz)*drQi(k,3)
              molQCdrv(jy,k,nm) = molQCdrv(jy,k,nm) + derv2(jy,iz)*drQi(k,3)
              molQCdrv(jz,k,nm) = molQCdrv(jz,k,nm) + derv2(jz,iz)*drQi(k,3)
            enddo
          enddo
!
!  Loop over second molecule
!
          do nn = 1,nmol
!
!  Set up quaternion matrix for second molecule
!
            call setquaternion(nn,molQ(1,nn),rQm2,drQm2,drQ2m2,.true.,.false.)
!
            do n = 1,nmolcore(nn)
              j = nmollist(nmolptr(nn)+n)
!**********************************************
!  Mixed quaternion-translation derivatives  *
!**********************************************
              jx = 3*(j - 1)
              do ia = 1,3
                do k = 1,3
                  molQTdrv(k,ia,nm,nn) = molQTdrv(k,ia,nm,nn) + drQi(k,1)*derv2(jx+ia,ix) + &
                                                                drQi(k,2)*derv2(jx+ia,iy) + &
                                                                drQi(k,3)*derv2(jx+ia,iz)
                enddo
              enddo
!**************************************
!  Quaternion-quaternion derivatives  *
!**************************************
              jx = 3*(j - 1) + 1
              jy = 3*(j - 1) + 2
              jz = 3*(j - 1) + 3
!
!  Set original intra molecular vector
!
              if (lsymopt) then
                nna = nmolf2a(nn)
                rjcom(1:3) = molQxyz(1:3,n,nna)
              else
                rjcom(1:3) = molQxyz(1:3,n,nn)
              endif
!
!  Transform by derivatives of quaternion rotation matrix
!
              drQj(1:3,1:3) = 0.0_dp
              do l = 1,3
                do k = 1,3
                  drQj(k,1) = drQj(k,1) + drQm2(k,1,l)*rjcom(l)
                  drQj(k,2) = drQj(k,2) + drQm2(k,2,l)*rjcom(l)
                  drQj(k,3) = drQj(k,3) + drQm2(k,3,l)*rjcom(l)
                enddo
              enddo
!
              do j1 = 1,3
                do k1 = 1,3
!
!  Term 1
!
                  molQQdrv(k1,j1,nn,nm) = molQQdrv(k1,j1,nn,nm) + derv2(jx,ix)*drQj(k1,1)*drQi(j1,1) &
                                                                + derv2(jy,ix)*drQj(k1,2)*drQi(j1,1) &
                                                                + derv2(jz,ix)*drQj(k1,3)*drQi(j1,1) &
                                                                + derv2(jx,iy)*drQj(k1,1)*drQi(j1,2) &
                                                                + derv2(jy,iy)*drQj(k1,2)*drQi(j1,2) &
                                                                + derv2(jz,iy)*drQj(k1,3)*drQi(j1,2) &
                                                                + derv2(jx,iz)*drQj(k1,1)*drQi(j1,3) &
                                                                + derv2(jy,iz)*drQj(k1,2)*drQi(j1,3) &
                                                                + derv2(jz,iz)*drQj(k1,3)*drQi(j1,3)
                enddo
              enddo
            enddo
          enddo
!
!  Transform by second derivatives of quaternion rotation matrix
!
          drQ2i(1:3,1:3,1:3) = 0.0_dp
          do l = 1,3
            do j = 1,3
              do k = 1,3
                drQ2i(k,j,1) = drQ2i(k,j,1) + drQ2(k,j,1,l)*ricom(l)
                drQ2i(k,j,2) = drQ2i(k,j,2) + drQ2(k,j,2,l)*ricom(l)
                drQ2i(k,j,3) = drQ2i(k,j,3) + drQ2(k,j,3,l)*ricom(l)
              enddo
            enddo
          enddo
!
!  Term 2 - self only
!
          do j = 1,3
            do k = 1,3
              molQQdrv(k,j,nm,nm) = molQQdrv(k,j,nm,nm) + drQ2i(k,j,1)*xdrv(i) &
                                                        + drQ2i(k,j,2)*ydrv(i) &
                                                        + drQ2i(k,j,3)*zdrv(i)
            enddo
          enddo
!
!  Strain mixed derivatives
!
          if (lstr) then
!
!  Remove cell volume factor from the mixed strain derivatives
!
            if (lfinitestrain) then
              derv3l(1,1:nstrains) = derv3(ix,1:nstrains)
              derv3l(2,1:nstrains) = derv3(iy,1:nstrains)
              derv3l(3,1:nstrains) = derv3(iz,1:nstrains)
!
!  Transform mixed derivatives by strain matrix
!
              do ns = 1,nstrains
                tmp(1:ndim) = 0.0_dp
                do j = 1,ndim
                  do l = 1,ndim
                    tmp(j) = tmp(j) + derv3l(l,ns)*strainfull(l,j)
                  enddo
                enddo
                derv3l(1:ndim,ns) = tmp(1:ndim)
              enddo
!
!  Remove strain derivatives
!
              if (ndim.eq.3) then
                derv3l(1,1) =  derv3l(1,1) - xdrv(i)
                derv3l(2,2) =  derv3l(2,2) - ydrv(i)
                derv3l(3,3) =  derv3l(3,3) - zdrv(i)
                derv3l(1,5) =  derv3l(1,5) - 0.5_dp*zdrv(i)
                derv3l(1,6) =  derv3l(1,6) - 0.5_dp*ydrv(i)
                derv3l(2,4) =  derv3l(2,4) - 0.5_dp*zdrv(i)
                derv3l(2,6) =  derv3l(2,6) - 0.5_dp*xdrv(i)
                derv3l(3,4) =  derv3l(3,4) - 0.5_dp*ydrv(i)
                derv3l(3,5) =  derv3l(3,5) - 0.5_dp*xdrv(i)
              elseif (ndim.eq.2) then
                derv3l(1,1) =  derv3l(1,1) - xdrv(i)
                derv3l(2,2) =  derv3l(2,2) - ydrv(i)
                derv3l(1,3) =  derv3l(1,3) - 0.5_dp*ydrv(i)
                derv3l(2,3) =  derv3l(2,3) - 0.5_dp*xdrv(i)
              elseif (ndim.eq.1) then
                derv3l(1,1) =  derv3l(1,1) - xdrv(i)
              endif
!
!  Transform mixed derivatives by inverse strain matrix
!
              do ns = 1,nstrains
                tmp(1:ndim) = 0.0_dp
                do j = 1,ndim
                  do l = 1,ndim
                    tmp(j) = tmp(j) + derv3l(l,ns)*straininverse(l,j)
                  enddo
                enddo
                derv3l(1:ndim,ns) = tmp(1:ndim)
              enddo
            else
              derv3l(1,1:nstrains) = derv3(ix,1:nstrains)
              derv3l(2,1:nstrains) = derv3(iy,1:nstrains)
              derv3l(3,1:nstrains) = derv3(iz,1:nstrains)
!
              if (ndim.eq.3) then
                derv3l(1,1) =  derv3(ix,1) - xdrv(i)
                derv3l(2,2) =  derv3(iy,2) - ydrv(i)
                derv3l(3,3) =  derv3(iz,3) - zdrv(i)
                derv3l(1,5) =  derv3(ix,5) - 0.5_dp*zdrv(i)
                derv3l(1,6) =  derv3(ix,6) - 0.5_dp*ydrv(i)
                derv3l(2,4) =  derv3(iy,4) - 0.5_dp*zdrv(i)
                derv3l(2,6) =  derv3(iy,6) - 0.5_dp*xdrv(i)
                derv3l(3,4) =  derv3(iz,4) - 0.5_dp*ydrv(i)
                derv3l(3,5) =  derv3(iz,5) - 0.5_dp*xdrv(i)
              elseif (ndim.eq.2) then
                derv3l(1,1) =  derv3(ix,1) - xdrv(i)
                derv3l(2,2) =  derv3(iy,2) - ydrv(i)
                derv3l(1,3) =  derv3(ix,3) - 0.5_dp*ydrv(i)
                derv3l(2,3) =  derv3(iy,3) - 0.5_dp*xdrv(i)
              elseif (ndim.eq.1) then
                derv3l(1,1) =  derv3(ix,1) - xdrv(i)
              endif
            endif
!
            do n = 1,nstrains
              do k = 1,3
                molQSdrv(n,k,nm) = molQSdrv(n,k,nm) + derv3l(1,n)*drQi(k,1)
                molQSdrv(n,k,nm) = molQSdrv(n,k,nm) + derv3l(2,n)*drQi(k,2)
                molQSdrv(n,k,nm) = molQSdrv(n,k,nm) + derv3l(3,n)*drQi(k,3)
              enddo
            enddo
          endif
        endif
      enddo
    enddo
  endif
#ifdef MPI
!
!  Globalise molecule derivatives
!
  if (lgrad2.and.nprocs.gt.1) then
    nsize = max(27*nmol*nmol,3*nmol*n3f,6*nmol*nstrains)
    allocate(sumin(nsize),stat=status)
    if (status/=0) call outofmemory('rigidmoleculedrv','sumin')
    allocate(sumout(nsize),stat=status)
    if (status/=0) call outofmemory('rigidmoleculedrv','sumout')
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
    call sumall(sumin,sumout,ind,"rigidmoleculedrv","molQCdrv")
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
    call sumall(sumin,sumout,ind,"rigidmoleculedrv","molTCdrv")
    ind = 0
    do nm = 1,nmol
      do ix = 1,3
        molTCdrv(1:n3f,ix,nm) = sumout(ind+1:ind+n3f)
        ind = ind + n3f
      enddo
    enddo
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
    call sumall(sumin,sumout,ind,"rigidmoleculedrv","molQQdrv")
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
    if (lstr) then
!   
!  molQSdrv, molTSdrv
!
      ind = 0
      do nm = 1,nmol
        do ix = 1,3
          do jx = 1,nstrains
            sumin(ind+1) = molQSdrv(jx,ix,nm)
            sumin(ind+2) = molTSdrv(jx,ix,nm)
            ind = ind + 2
          enddo
        enddo
      enddo
      call sumall(sumin,sumout,ind,"rigidmoleculedrv","molQSdrv")
      ind = 0
      do nm = 1,nmol
        do ix = 1,3
          do jx = 1,nstrains
            molQSdrv(jx,ix,nm) = sumout(ind+1)
            molTSdrv(jx,ix,nm) = sumout(ind+2)
            ind = ind + 2
          enddo
        enddo
      enddo
    endif
!
    deallocate(sumout,stat=status)
    if (status/=0) call deallocate_error('rigidmoleculedrv','sumout')
    deallocate(sumin,stat=status)
    if (status/=0) call deallocate_error('rigidmoleculedrv','sumin')
  endif
#endif
!
  return
  end
