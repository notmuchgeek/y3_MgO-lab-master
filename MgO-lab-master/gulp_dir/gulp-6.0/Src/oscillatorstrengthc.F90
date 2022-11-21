  subroutine oscillatorstrengthc(mtv,mtvloc,mtvptr,nphonatc,nphonatptr,ncfoc,iocptr,eigc,maxd2,oscstrength)
!
!  Calculates the oscillator strengths for the modes.
!  Complex version of oscillatorstrength
!
!  12/16 Created from oscillatorstrength
!   2/18 Trace added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!   5/20 Rigid molecule modifications added
!   6/20 nmolcore changes added
!
!  On entry :
!
!  mtv         = no. of modes
!  mtvloc      = no. of modes on local node
!  mtvptr      = pointer from local to global mode
!  nphonatc    = total number of cores in relevant regions
!  nphonatptr  = pointer from reduce to full atom set
!  ncfoc       = number of condensed core sites
!  iocptr      = pointer that connects sites to condensed sites
!  eigc        = eigenvectors of dynamical matrix - complex
!  maxd2       = left-hand dimension of eigc
!
!  On exit : 
!
!  oscstrength = oscillator strengths for each mode as a 3 x 3 tensor per mode
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
  use control,        only : lrigid
  use current
  use element
  use frequencies
  use iochannels
  use molecule
  use parallel,       only : nprocs
  use species,        only : massspec
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4),                   intent(in)    :: iocptr(*)
  integer(i4),                   intent(in)    :: maxd2
  integer(i4),                   intent(in)    :: mtv                   ! Total number of modes = mcv + mmv
  integer(i4),                   intent(in)    :: mtvloc
  integer(i4),                   intent(in)    :: mtvptr(mtvloc)
  integer(i4),                   intent(in)    :: nphonatc
  integer(i4),                   intent(in)    :: nphonatptr(*)
  integer(i4),                   intent(in)    :: ncfoc
  complex(dpc),                  intent(in)    :: eigc(maxd2,*)
  real(dp),                      intent(out)   :: oscstrength(3,3,mtv)
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ind
  integer(i4)                                  :: ix
  integer(i4)                                  :: iy
  integer(i4)                                  :: iz
  integer(i4)                                  :: j
  integer(i4)                                  :: k
  integer(i4)                                  :: l
  integer(i4)                                  :: m
  integer(i4)                                  :: mm
  integer(i4)                                  :: status
  real(dp)                                     :: drQ(3,3)
  complex(dpc)                                 :: osx
  complex(dpc)                                 :: osy
  complex(dpc)                                 :: osz
  complex(dpc)                                 :: reig(3)
  real(dp)                                     :: trmj
  real(dp),     dimension(:,:,:), allocatable  :: sum3
#ifdef TRACE
  call trace_in('oscillatorstrengthc')
#endif
!**********************************
!  Oscillator strengths per mode  *
!**********************************
!
!  Zero arrays
!
  oscstrength(1:3,1:3,1:mtv) = 0.0_dp
  if (lrigid) then
!
!  Loop over modes
!
    do m = 1,mtvloc
      mm = mtvptr(m)
      osx = 0.0_dpc
      osy = 0.0_dpc
      osz = 0.0_dpc
!
!  Loop over atoms
!
      do i = 1,ncorenomol
        ix = 3*(i-1) + 1
        iy = ix + 1
        iz = ix + 2
        j = ncorenomolptr(i)
!
!  Multiple inverse mass weighted eigenvectors by Born charges
!
        osx = osx + (eigc(ix,m)*bornq(1,1,j)*rfmass(ix) + &
                     eigc(iy,m)*bornq(2,1,j)*rfmass(iy) + &
                     eigc(iz,m)*bornq(3,1,j)*rfmass(iz))*occuf(j)
        osy = osy + (eigc(ix,m)*bornq(1,2,j)*rfmass(ix) + &
                     eigc(iy,m)*bornq(2,2,j)*rfmass(iy) + &
                     eigc(iz,m)*bornq(3,2,j)*rfmass(iz))*occuf(j)
        osz = osz + (eigc(ix,m)*bornq(1,3,j)*rfmass(ix) + &
                     eigc(iy,m)*bornq(2,3,j)*rfmass(iy) + &
                     eigc(iz,m)*bornq(3,3,j)*rfmass(iz))*occuf(j)
      enddo
!
!  Loop over atoms in molecules
!
      do i = 1,nmol
        ind = 6*(i-1) + 3*ncorenomol
        do j = 1,nmolcore(i)
          k = nmollist(nmolptr(i)+j)
!
!  Multiple inverse mass weighted eigenvectors for translation by Born charges
!
          osx = osx + (eigc(ind+1,m)*bornq(1,1,k)*rfmass(ind+1) + &
                       eigc(ind+2,m)*bornq(2,1,k)*rfmass(ind+2) + &
                       eigc(ind+3,m)*bornq(3,1,k)*rfmass(ind+3))*occuf(k)
          osy = osy + (eigc(ind+1,m)*bornq(1,2,k)*rfmass(ind+1) + &
                       eigc(ind+2,m)*bornq(2,2,k)*rfmass(ind+2) + &
                       eigc(ind+3,m)*bornq(3,2,k)*rfmass(ind+3))*occuf(k)
          osz = osz + (eigc(ind+1,m)*bornq(1,3,k)*rfmass(ind+1) + &
                       eigc(ind+2,m)*bornq(2,3,k)*rfmass(ind+2) + &
                       eigc(ind+3,m)*bornq(3,3,k)*rfmass(ind+3))*occuf(k)
!
!  Take cross product of axes with atom vector
!
          do l = 1,3
            drQ(1,l) = molaxes(2,l,i)*molxyz(3,j,i) - molaxes(2,l,i)*molxyz(3,j,i)
            drQ(2,l) = molaxes(3,l,i)*molxyz(1,j,i) - molaxes(1,l,i)*molxyz(3,j,i)
            drQ(3,l) = molaxes(1,l,i)*molxyz(2,j,i) - molaxes(2,l,i)*molxyz(1,j,i)
          enddo
!
!  Multiply rotation contribution by eigenvector
!
          reig(1:3) = 0.0_dpc
          do l = 1,3
            reig(l) = reig(l) + drQ(l,1)*eigc(ind+4,m)*rfmass(ind+4)
            reig(l) = reig(l) + drQ(l,2)*eigc(ind+5,m)*rfmass(ind+5)
            reig(l) = reig(l) + drQ(l,3)*eigc(ind+6,m)*rfmass(ind+6)
          enddo
!
          osx = osx + (reig(1)*bornq(1,1,k) + &
                       reig(2)*bornq(2,1,k) + &
                       reig(3)*bornq(3,1,k))*occuf(k)
          osy = osy + (reig(1)*bornq(1,2,k) + &
                       reig(2)*bornq(2,2,k) + &
                       reig(3)*bornq(3,2,k))*occuf(k)
          osz = osz + (reig(1)*bornq(1,3,k) + &
                       reig(2)*bornq(2,3,k) + &
                       reig(3)*bornq(3,3,k))*occuf(k)
        enddo
      enddo
!
!  Rigid molecule version of algorithm
!
      oscstrength(1,1,mm) = dble(conjg(osx)*osx)
      oscstrength(2,1,mm) = dble(conjg(osy)*osx)
      oscstrength(3,1,mm) = dble(conjg(osz)*osx)
      oscstrength(1,2,mm) = dble(conjg(osx)*osy)
      oscstrength(2,2,mm) = dble(conjg(osy)*osy)
      oscstrength(3,2,mm) = dble(conjg(osz)*osy)
      oscstrength(1,3,mm) = dble(conjg(osx)*osz)
      oscstrength(2,3,mm) = dble(conjg(osy)*osz)
      oscstrength(3,3,mm) = dble(conjg(osz)*osz)
    enddo
  else
!
!  Loop over modes
!
    do m = 1,mtvloc
      mm = mtvptr(m)
      osx = 0.0_dpc
      osy = 0.0_dpc
      osz = 0.0_dpc
!
!  Loop over full sites
!
      do i = 1,ncfoc
        ix = 3*(i-1) + 1
        iy = ix + 1
        iz = ix + 2
!
!  Find all cores associated with full site
!
        do j = 1,nphonatc
          if (iocptr(j).eq.i) then
            trmj = occuf(j)/sqrt(massspec(nspecptr(nrelf2a(nphonatptr(j)))))
!
!  Multiple inverse mass weighted eigenvectors by Born charges
!
            osx = osx + (eigc(ix,m)*bornq(1,1,j) + &
                        eigc(iy,m)*bornq(2,1,j) + &
                        eigc(iz,m)*bornq(3,1,j))*trmj
            osy = osy + (eigc(ix,m)*bornq(1,2,j) + &
                        eigc(iy,m)*bornq(2,2,j) + &
                        eigc(iz,m)*bornq(3,2,j))*trmj
            osz = osz + (eigc(ix,m)*bornq(1,3,j) + &
                        eigc(iy,m)*bornq(2,3,j) + &
                        eigc(iz,m)*bornq(3,3,j))*trmj
          endif
        enddo
      enddo
      oscstrength(1,1,mm) = dble(conjg(osx)*osx)
      oscstrength(2,1,mm) = dble(conjg(osy)*osx)
      oscstrength(3,1,mm) = dble(conjg(osz)*osx)
      oscstrength(1,2,mm) = dble(conjg(osx)*osy)
      oscstrength(2,2,mm) = dble(conjg(osy)*osy)
      oscstrength(3,2,mm) = dble(conjg(osz)*osy)
      oscstrength(1,3,mm) = dble(conjg(osx)*osz)
      oscstrength(2,3,mm) = dble(conjg(osy)*osz)
      oscstrength(3,3,mm) = dble(conjg(osz)*osz)
    enddo
  endif
  if (nprocs.gt.1) then
!
!  Globalise oscillator strengths
!
    allocate(sum3(3,3,mtv),stat=status)
    if (status/=0) call outofmemory('oscillatorstrengthc','sum3')
!
    call sumall(oscstrength,sum3,9_i4*mtv,"oscillatorstrengthc","oscstrength")
    oscstrength(1:3,1:3,1:mtv) = sum3(1:3,1:3,1:mtv)
!
    deallocate(sum3,stat=status)
    if (status/=0) call deallocate_error('oscillatorstrengthc','sum3')
  endif
#ifdef TRACE
  call trace_out('oscillatorstrengthc')
#endif
!
  return
  end
