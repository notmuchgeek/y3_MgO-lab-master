  subroutine rigidmoleculeprop
!
!  Set second derivatives for rigid molecules during a property calculation
!  NB: Uses rotation angles rather than quaternions.
!
!   5/20 Created from rigidmoleculephon
!   5/20 Correction for systems that use symmetry by changing xdrv for xfdrv etc
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
!  Julian Gale, CIC, Curtin University, June 2020
!
  use current
  use derivatives
  use m_strain,       only : strainfull, straininverse
  use molecule
  use parallel
  use symmetry,       only : lstr
  implicit none
!
!  Local variables
!
  integer(i4)                                   :: i
  integer(i4)                                   :: ia
  integer(i4)                                   :: ii
#ifdef MPI
  integer(i4)                                   :: ind
#endif
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
  integer(i4)                                   :: ns
#ifdef MPI
  integer(i4)                                   :: nsize
  integer(i4)                                   :: status
#endif
  real(dp)                                      :: derv3l(3,6)
  real(dp)                                      :: drR(3,3,3)
  real(dp)                                      :: drS(3,3,3)
  real(dp)                                      :: drRi(3,3)
  real(dp)                                      :: drRj(3,3)
#ifdef MPI
  real(dp),   dimension(:),   allocatable       :: sumin
  real(dp),   dimension(:),   allocatable       :: sumout
#endif
  real(dp)                                      :: tmp(3)
!
  n3f = 3*numat
!
!  Initialise derivatives of molecules : Note that Q here will be used for R
!
  molQCdrv(1:n3f,1:3,1:nmol) = 0.0_dp
  molTCdrv(1:n3f,1:3,1:nmol) = 0.0_dp
  molQQdrv(1:3,1:3,1:nmol,1:nmol) = 0.0_dp
  molQTdrv(1:3,1:3,1:nmol,1:nmol) = 0.0_dp
  molTTdrv(1:3,1:3,1:nmol,1:nmol) = 0.0_dp
  if (lstr) then
    molQSdrv(1:nstrains,1:3,1:nmol) = 0.0_dp
    molTSdrv(1:nstrains,1:3,1:nmol) = 0.0_dp
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
!
!  Sum of translational second derivatives
!
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
  if (nprocs.gt.1) then
!-----------------------
!  Parallel algorithm  |
!-----------------------
!
!  Loop over atoms to handle rotation derivatives
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
          do nn = 1,numatnomol
            j = numatnomolptr(nn)
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
                derv3l(1,1) =  derv3l(1,1) - xfdrv(i)
                derv3l(2,2) =  derv3l(2,2) - yfdrv(i)
                derv3l(3,3) =  derv3l(3,3) - zfdrv(i)
                derv3l(1,5) =  derv3l(1,5) - 0.5_dp*zfdrv(i)
                derv3l(1,6) =  derv3l(1,6) - 0.5_dp*yfdrv(i)
                derv3l(2,4) =  derv3l(2,4) - 0.5_dp*zfdrv(i)
                derv3l(2,6) =  derv3l(2,6) - 0.5_dp*xfdrv(i)
                derv3l(3,4) =  derv3l(3,4) - 0.5_dp*yfdrv(i)
                derv3l(3,5) =  derv3l(3,5) - 0.5_dp*xfdrv(i)
              elseif (ndim.eq.2) then
                derv3l(1,1) =  derv3l(1,1) - xfdrv(i)
                derv3l(2,2) =  derv3l(2,2) - yfdrv(i)
                derv3l(1,3) =  derv3l(1,3) - 0.5_dp*yfdrv(i)
                derv3l(2,3) =  derv3l(2,3) - 0.5_dp*xfdrv(i)
              elseif (ndim.eq.1) then
                derv3l(1,1) =  derv3l(1,1) - xfdrv(i)
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
                derv3l(1,1) =  derv3(ix,1) - xfdrv(i)
                derv3l(2,2) =  derv3(iy,2) - yfdrv(i)
                derv3l(3,3) =  derv3(iz,3) - zfdrv(i)
                derv3l(1,5) =  derv3(ix,5) - 0.5_dp*zfdrv(i)
                derv3l(1,6) =  derv3(ix,6) - 0.5_dp*yfdrv(i)
                derv3l(2,4) =  derv3(iy,4) - 0.5_dp*zfdrv(i)
                derv3l(2,6) =  derv3(iy,6) - 0.5_dp*xfdrv(i)
                derv3l(3,4) =  derv3(iz,4) - 0.5_dp*yfdrv(i)
                derv3l(3,5) =  derv3(iz,5) - 0.5_dp*xfdrv(i)
              elseif (ndim.eq.2) then
                derv3l(1,1) =  derv3(ix,1) - xfdrv(i)
                derv3l(2,2) =  derv3(iy,2) - yfdrv(i)
                derv3l(1,3) =  derv3(ix,3) - 0.5_dp*yfdrv(i)
                derv3l(2,3) =  derv3(iy,3) - 0.5_dp*xfdrv(i)
              elseif (ndim.eq.1) then
                derv3l(1,1) =  derv3(ix,1) - xfdrv(i)
              endif
            endif
!
            do n = 1,nstrains
              do k = 1,3
                molQSdrv(n,1,nm) = molQSdrv(n,1,nm) + derv3l(k,n)*drRi(k,1)
                molQSdrv(n,2,nm) = molQSdrv(n,2,nm) + derv3l(k,n)*drRi(k,2)
                molQSdrv(n,3,nm) = molQSdrv(n,3,nm) + derv3l(k,n)*drRi(k,3)
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
!  Loop over molecules to handle rotation derivatives
!
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
        do nn = 1,numatnomol
          j = numatnomolptr(nn)
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
              derv3l(1,1) =  derv3l(1,1) - xfdrv(i)
              derv3l(2,2) =  derv3l(2,2) - yfdrv(i)
              derv3l(3,3) =  derv3l(3,3) - zfdrv(i)
              derv3l(1,5) =  derv3l(1,5) - 0.5_dp*zfdrv(i)
              derv3l(1,6) =  derv3l(1,6) - 0.5_dp*yfdrv(i)
              derv3l(2,4) =  derv3l(2,4) - 0.5_dp*zfdrv(i)
              derv3l(2,6) =  derv3l(2,6) - 0.5_dp*xfdrv(i)
              derv3l(3,4) =  derv3l(3,4) - 0.5_dp*yfdrv(i)
              derv3l(3,5) =  derv3l(3,5) - 0.5_dp*xfdrv(i)
            elseif (ndim.eq.2) then
              derv3l(1,1) =  derv3l(1,1) - xfdrv(i)
              derv3l(2,2) =  derv3l(2,2) - yfdrv(i)
              derv3l(1,3) =  derv3l(1,3) - 0.5_dp*yfdrv(i)
              derv3l(2,3) =  derv3l(2,3) - 0.5_dp*xfdrv(i)
            elseif (ndim.eq.1) then
              derv3l(1,1) =  derv3l(1,1) - xfdrv(i)
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
              derv3l(1,1) =  derv3(ix,1) - xfdrv(i)
              derv3l(2,2) =  derv3(iy,2) - yfdrv(i)
              derv3l(3,3) =  derv3(iz,3) - zfdrv(i)
              derv3l(1,5) =  derv3(ix,5) - 0.5_dp*zfdrv(i)
              derv3l(1,6) =  derv3(ix,6) - 0.5_dp*yfdrv(i)
              derv3l(2,4) =  derv3(iy,4) - 0.5_dp*zfdrv(i)
              derv3l(2,6) =  derv3(iy,6) - 0.5_dp*xfdrv(i)
              derv3l(3,4) =  derv3(iz,4) - 0.5_dp*yfdrv(i)
              derv3l(3,5) =  derv3(iz,5) - 0.5_dp*xfdrv(i)
            elseif (ndim.eq.2) then
              derv3l(1,1) =  derv3(ix,1) - xfdrv(i)
              derv3l(2,2) =  derv3(iy,2) - yfdrv(i)
              derv3l(1,3) =  derv3(ix,3) - 0.5_dp*yfdrv(i)
              derv3l(2,3) =  derv3(iy,3) - 0.5_dp*xfdrv(i)
            elseif (ndim.eq.1) then
              derv3l(1,1) =  derv3(ix,1) - xfdrv(i)
            endif
          endif
!
          do n = 1,nstrains
            do k = 1,3
              molQSdrv(n,1,nm) = molQSdrv(n,1,nm) + derv3l(k,n)*drRi(k,1)
              molQSdrv(n,2,nm) = molQSdrv(n,2,nm) + derv3l(k,n)*drRi(k,2)
              molQSdrv(n,3,nm) = molQSdrv(n,3,nm) + derv3l(k,n)*drRi(k,3)
            enddo
          enddo
        endif
      enddo
    enddo
  endif
#ifdef MPI
!
!  Globalise molecule derivatives
!
  if (nprocs.gt.1) then
    nsize = max(27*nmol*nmol,3*nmol*n3f,6*nmol*nstrains)
    allocate(sumin(nsize),stat=status)
    if (status/=0) call outofmemory('rigidmoleculeprop','sumin')
    allocate(sumout(nsize),stat=status)
    if (status/=0) call outofmemory('rigidmoleculeprop','sumout')
    if (numatnomol.gt.0) then
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
      call sumall(sumin,sumout,ind,"rigidmoleculeprop","molQCdrv")
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
      call sumall(sumin,sumout,ind,"rigidmoleculeprop","molTCdrv")
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
    call sumall(sumin,sumout,ind,"rigidmoleculeprop","molQQdrv")
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
      call sumall(sumin,sumout,ind,"rigidmoleculeprop","molQSdrv")
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
    if (status/=0) call deallocate_error('rigidmoleculeprop','sumout')
    deallocate(sumin,stat=status)
    if (status/=0) call deallocate_error('rigidmoleculeprop','sumin')
  endif
#endif
!
  return
  end
