  subroutine rigidmoleculegv
!
!  Set second derivatives for rigid molecules during calculation of group velocities
!  NB: Uses rotation angles rather than quaternions.
!
!   5/20 Created from rigidmoleculephon
!   6/20 Parallelisation added
!   6/20 nmolcore changes added
!   6/20 Phasing of molQQdk corrected back
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
#ifdef MPI
  integer(i4)                                   :: kx
#endif
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
  complex(dpc)                                  :: dQ(3,3)
#ifdef MPI
  complex(dpc), dimension(:),   allocatable     :: sumin
  complex(dpc), dimension(:),   allocatable     :: sumout
#endif
  real(dp)                                      :: drR(3,3,3)
  real(dp)                                      :: drS(3,3,3)
  real(dp)                                      :: drRi(3,3)
  real(dp)                                      :: drRj(3,3)
!
  n3f = 3*numat
!
!  Initialise derivatives of molecules : Note that Q here will be used for R
!
  molQCdk(1:3,1:n3f,1:3,1:nmol) = 0.0_dpc
  molTCdk(1:3,1:n3f,1:3,1:nmol) = 0.0_dpc
  molQQdk(1:3,1:3,1:3,1:nmol,1:nmol) = 0.0_dpc
  molQTdk(1:3,1:3,1:3,1:nmol,1:nmol) = 0.0_dpc
  molTTdk(1:3,1:3,1:3,1:nmol,1:nmol) = 0.0_dpc
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
!
              molTTdk(1:3,1,1,nn,nm) = molTTdk(1:3,1,1,nn,nm) + derv2dk(1:3,jx,ix)
              molTTdk(1:3,2,1,nn,nm) = molTTdk(1:3,2,1,nn,nm) + derv2dk(1:3,jy,ix)
              molTTdk(1:3,3,1,nn,nm) = molTTdk(1:3,3,1,nn,nm) + derv2dk(1:3,jz,ix)
              molTTdk(1:3,1,2,nn,nm) = molTTdk(1:3,1,2,nn,nm) + derv2dk(1:3,jx,iy)
              molTTdk(1:3,2,2,nn,nm) = molTTdk(1:3,2,2,nn,nm) + derv2dk(1:3,jy,iy)
              molTTdk(1:3,3,2,nn,nm) = molTTdk(1:3,3,2,nn,nm) + derv2dk(1:3,jz,iy)
              molTTdk(1:3,1,3,nn,nm) = molTTdk(1:3,1,3,nn,nm) + derv2dk(1:3,jx,iz)
              molTTdk(1:3,2,3,nn,nm) = molTTdk(1:3,2,3,nn,nm) + derv2dk(1:3,jy,iz)
              molTTdk(1:3,3,3,nn,nm) = molTTdk(1:3,3,3,nn,nm) + derv2dk(1:3,jz,iz)
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
            molTCdk(1:3,jx,1,nm) = molTCdk(1:3,jx,1,nm) + derv2dk(1:3,jx,ix)
            molTCdk(1:3,jy,1,nm) = molTCdk(1:3,jy,1,nm) + derv2dk(1:3,jy,ix)
            molTCdk(1:3,jz,1,nm) = molTCdk(1:3,jz,1,nm) + derv2dk(1:3,jz,ix)
            molTCdk(1:3,jx,2,nm) = molTCdk(1:3,jx,2,nm) + derv2dk(1:3,jx,iy)
            molTCdk(1:3,jy,2,nm) = molTCdk(1:3,jy,2,nm) + derv2dk(1:3,jy,iy)
            molTCdk(1:3,jz,2,nm) = molTCdk(1:3,jz,2,nm) + derv2dk(1:3,jz,iy)
            molTCdk(1:3,jx,3,nm) = molTCdk(1:3,jx,3,nm) + derv2dk(1:3,jx,iz)
            molTCdk(1:3,jy,3,nm) = molTCdk(1:3,jy,3,nm) + derv2dk(1:3,jy,iz)
            molTCdk(1:3,jz,3,nm) = molTCdk(1:3,jz,3,nm) + derv2dk(1:3,jz,iz)
          enddo
        endif
      endif
    enddo
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
          do nn = 1,ncorenomol
            j = ncorenomolptr(nn)
            jx = 3*(j - 1) + 1
            jy = 3*(j - 1) + 2
            jz = 3*(j - 1) + 3
            do k = 1,3
              molQCdk(1:3,jx,k,nm) = molQCdk(1:3,jx,k,nm) + derv2dk(1:3,jx,ix)*drRi(1,k)
              molQCdk(1:3,jy,k,nm) = molQCdk(1:3,jy,k,nm) + derv2dk(1:3,jy,ix)*drRi(1,k)
              molQCdk(1:3,jz,k,nm) = molQCdk(1:3,jz,k,nm) + derv2dk(1:3,jz,ix)*drRi(1,k)
              molQCdk(1:3,jx,k,nm) = molQCdk(1:3,jx,k,nm) + derv2dk(1:3,jx,iy)*drRi(2,k)
              molQCdk(1:3,jy,k,nm) = molQCdk(1:3,jy,k,nm) + derv2dk(1:3,jy,iy)*drRi(2,k)
              molQCdk(1:3,jz,k,nm) = molQCdk(1:3,jz,k,nm) + derv2dk(1:3,jz,iy)*drRi(2,k)
              molQCdk(1:3,jx,k,nm) = molQCdk(1:3,jx,k,nm) + derv2dk(1:3,jx,iz)*drRi(3,k)
              molQCdk(1:3,jy,k,nm) = molQCdk(1:3,jy,k,nm) + derv2dk(1:3,jy,iz)*drRi(3,k)
              molQCdk(1:3,jz,k,nm) = molQCdk(1:3,jz,k,nm) + derv2dk(1:3,jz,iz)*drRi(3,k)
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
!
              do ia = 1,3
                do k = 1,3
                  molQTdk(1:3,k,ia,nm,nn) = molQTdk(1:3,k,ia,nm,nn) + drRi(1,k)*derv2dk(1:3,jx+ia,ix) + &
                                                                      drRi(2,k)*derv2dk(1:3,jx+ia,iy) + &
                                                                      drRi(3,k)*derv2dk(1:3,jx+ia,iz)
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
                  molQQdk(1:3,k1,j1,nn,nm) = molQQdk(1:3,k1,j1,nn,nm) + derv2dk(1:3,jx,ix)*drRj(1,k1)*drRi(1,j1) &
                                                                      + derv2dk(1:3,jy,ix)*drRj(2,k1)*drRi(1,j1) &
                                                                      + derv2dk(1:3,jz,ix)*drRj(3,k1)*drRi(1,j1) &
                                                                      + derv2dk(1:3,jx,iy)*drRj(1,k1)*drRi(2,j1) &
                                                                      + derv2dk(1:3,jy,iy)*drRj(2,k1)*drRi(2,j1) &
                                                                      + derv2dk(1:3,jz,iy)*drRj(3,k1)*drRi(2,j1) &
                                                                      + derv2dk(1:3,jx,iz)*drRj(1,k1)*drRi(3,j1) &
                                                                      + derv2dk(1:3,jy,iz)*drRj(2,k1)*drRi(3,j1) &
                                                                      + derv2dk(1:3,jz,iz)*drRj(3,k1)*drRi(3,j1)
                enddo
              enddo
            enddo
          enddo
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
            molTTdk(1:3,1,1,nn,nm) = molTTdk(1:3,1,1,nn,nm) + derv2dk(1:3,jx,ix)
            molTTdk(1:3,2,1,nn,nm) = molTTdk(1:3,2,1,nn,nm) + derv2dk(1:3,jy,ix)
            molTTdk(1:3,3,1,nn,nm) = molTTdk(1:3,3,1,nn,nm) + derv2dk(1:3,jz,ix)
            molTTdk(1:3,1,2,nn,nm) = molTTdk(1:3,1,2,nn,nm) + derv2dk(1:3,jx,iy)
            molTTdk(1:3,2,2,nn,nm) = molTTdk(1:3,2,2,nn,nm) + derv2dk(1:3,jy,iy)
            molTTdk(1:3,3,2,nn,nm) = molTTdk(1:3,3,2,nn,nm) + derv2dk(1:3,jz,iy)
            molTTdk(1:3,1,3,nn,nm) = molTTdk(1:3,1,3,nn,nm) + derv2dk(1:3,jx,iz)
            molTTdk(1:3,2,3,nn,nm) = molTTdk(1:3,2,3,nn,nm) + derv2dk(1:3,jy,iz)
            molTTdk(1:3,3,3,nn,nm) = molTTdk(1:3,3,3,nn,nm) + derv2dk(1:3,jz,iz)
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
          molTCdk(1:3,jx,1,nm) = molTCdk(1:3,jx,1,nm) + derv2dk(1:3,jx,ix)
          molTCdk(1:3,jy,1,nm) = molTCdk(1:3,jy,1,nm) + derv2dk(1:3,jy,ix)
          molTCdk(1:3,jz,1,nm) = molTCdk(1:3,jz,1,nm) + derv2dk(1:3,jz,ix)
          molTCdk(1:3,jx,2,nm) = molTCdk(1:3,jx,2,nm) + derv2dk(1:3,jx,iy)
          molTCdk(1:3,jy,2,nm) = molTCdk(1:3,jy,2,nm) + derv2dk(1:3,jy,iy)
          molTCdk(1:3,jz,2,nm) = molTCdk(1:3,jz,2,nm) + derv2dk(1:3,jz,iy)
          molTCdk(1:3,jx,3,nm) = molTCdk(1:3,jx,3,nm) + derv2dk(1:3,jx,iz)
          molTCdk(1:3,jy,3,nm) = molTCdk(1:3,jy,3,nm) + derv2dk(1:3,jy,iz)
          molTCdk(1:3,jz,3,nm) = molTCdk(1:3,jz,3,nm) + derv2dk(1:3,jz,iz)
        enddo
      enddo
    enddo
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
        do nn = 1,ncorenomol
          j = ncorenomolptr(nn)
          jx = 3*(j - 1) + 1
          jy = 3*(j - 1) + 2
          jz = 3*(j - 1) + 3
          do k = 1,3
            molQCdk(1:3,jx,k,nm) = molQCdk(1:3,jx,k,nm) + derv2dk(1:3,jx,ix)*drRi(1,k)
            molQCdk(1:3,jy,k,nm) = molQCdk(1:3,jy,k,nm) + derv2dk(1:3,jy,ix)*drRi(1,k)
            molQCdk(1:3,jz,k,nm) = molQCdk(1:3,jz,k,nm) + derv2dk(1:3,jz,ix)*drRi(1,k)
            molQCdk(1:3,jx,k,nm) = molQCdk(1:3,jx,k,nm) + derv2dk(1:3,jx,iy)*drRi(2,k)
            molQCdk(1:3,jy,k,nm) = molQCdk(1:3,jy,k,nm) + derv2dk(1:3,jy,iy)*drRi(2,k)
            molQCdk(1:3,jz,k,nm) = molQCdk(1:3,jz,k,nm) + derv2dk(1:3,jz,iy)*drRi(2,k)
            molQCdk(1:3,jx,k,nm) = molQCdk(1:3,jx,k,nm) + derv2dk(1:3,jx,iz)*drRi(3,k)
            molQCdk(1:3,jy,k,nm) = molQCdk(1:3,jy,k,nm) + derv2dk(1:3,jy,iz)*drRi(3,k)
            molQCdk(1:3,jz,k,nm) = molQCdk(1:3,jz,k,nm) + derv2dk(1:3,jz,iz)*drRi(3,k)
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
              dQ(1:3,1) = derv2dk(1:3,ix,jx+ia)
              dQ(1:3,2) = derv2dk(1:3,iy,jx+ia)
              dQ(1:3,3) = derv2dk(1:3,iz,jx+ia)
              do k = 1,3
                molQTdk(1:3,k,ia,nm,nn) = molQTdk(1:3,k,ia,nm,nn) + drRi(1,k)*dQ(1:3,1) &
                                                                  + drRi(2,k)*dQ(1:3,2) &
                                                                  + drRi(3,k)*dQ(1:3,3)
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
                molQQdk(1:3,k1,j1,nn,nm) = molQQdk(1:3,k1,j1,nn,nm) + derv2dk(1:3,jx,ix)*drRj(1,k1)*drRi(1,j1) &
                                                                    + derv2dk(1:3,jy,ix)*drRj(2,k1)*drRi(1,j1) &
                                                                    + derv2dk(1:3,jz,ix)*drRj(3,k1)*drRi(1,j1) &
                                                                    + derv2dk(1:3,jx,iy)*drRj(1,k1)*drRi(2,j1) &
                                                                    + derv2dk(1:3,jy,iy)*drRj(2,k1)*drRi(2,j1) &
                                                                    + derv2dk(1:3,jz,iy)*drRj(3,k1)*drRi(2,j1) &
                                                                    + derv2dk(1:3,jx,iz)*drRj(1,k1)*drRi(3,j1) &
                                                                    + derv2dk(1:3,jy,iz)*drRj(2,k1)*drRi(3,j1) &
                                                                    + derv2dk(1:3,jz,iz)*drRj(3,k1)*drRi(3,j1)
              enddo
            enddo
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
    nsize = max(81*nmol*nmol,9*nmol*n3f)
    allocate(sumin(nsize),stat=status)
    if (status/=0) call outofmemory('rigidmoleculegv','sumin')
    allocate(sumout(nsize),stat=status)
    if (status/=0) call outofmemory('rigidmoleculegv','sumout')
    if (numatnomol.gt.0) then
!
!  molQCdk
!
      ind = 0
      do nm = 1,nmol
        do ix = 1,3
          do jx = 1,3
            sumin(ind+1:ind+n3f) = molQCdk(jx,1:n3f,ix,nm)
            ind = ind + n3f
          enddo
        enddo
      enddo
      call csumall(sumin,sumout,ind,"rigidmoleculegv","molQCdk")
      ind = 0
      do nm = 1,nmol
        do ix = 1,3
          do jx = 1,3
            molQCdk(jx,1:n3f,ix,nm) = sumout(ind+1:ind+n3f)
            ind = ind + n3f
          enddo
        enddo
      enddo
!
!  molTCdk
!
      ind = 0
      do nm = 1,nmol
        do ix = 1,3
          do jx = 1,3
            sumin(ind+1:ind+n3f) = molTCdk(jx,1:n3f,ix,nm)
            ind = ind + n3f
          enddo
        enddo
      enddo
      call csumall(sumin,sumout,ind,"rigidmoleculegv","molTCdk")
      ind = 0
      do nm = 1,nmol
        do ix = 1,3
          do jx = 1,3
            molTCdk(jx,1:n3f,ix,nm) = sumout(ind+1:ind+n3f)
            ind = ind + n3f
          enddo
        enddo
      enddo
    endif
!
!  molQQdk, molQTdk, molTTdk
!
    ind = 0
    do nm = 1,nmol
      do nn = 1,nmol
        do ix = 1,3
          do jx = 1,3
            do kx = 1,3
              sumin(ind+1) = molQQdk(kx,jx,ix,nn,nm)
              sumin(ind+2) = molQTdk(kx,jx,ix,nn,nm)
              sumin(ind+3) = molTTdk(kx,jx,ix,nn,nm)
              ind = ind + 3
            enddo
          enddo
        enddo
      enddo
    enddo
    call csumall(sumin,sumout,ind,"rigidmoleculegv","molQQdk")
    ind = 0
    do nm = 1,nmol
      do nn = 1,nmol
        do ix = 1,3
          do jx = 1,3
            do kx = 1,3
              molQQdk(kx,jx,ix,nn,nm) = sumout(ind+1)
              molQTdk(kx,jx,ix,nn,nm) = sumout(ind+2)
              molTTdk(kx,jx,ix,nn,nm) = sumout(ind+3)
              ind = ind + 3
            enddo
          enddo
        enddo
      enddo
    enddo
!
    deallocate(sumout,stat=status)
    if (status/=0) call deallocate_error('rigidmoleculegv','sumout')
    deallocate(sumin,stat=status)
    if (status/=0) call deallocate_error('rigidmoleculegv','sumin')
  endif
#endif
!
  return
  end
