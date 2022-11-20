  subroutine oscillatorstrength(mtv,nphonatc,nphonatptr,ncfoc,iocptr,eigr,eigi,maxd2,oscstrength)
!
!  Calculates the oscillator strengths for the modes
!
!   3/02 Created from peigeng
!   7/02 Region 1 only modifications added
!   5/06 Mass now uses species values
!   8/12 Site occupancy introduced to scale oscillator strength
!   8/12 Created by adding imaginary eigenvectors to oscillatorstrengthg
!   1/18 Modified for ghostcell case
!   2/18 Trace added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!   5/20 Rigid molecule modifications added
!   6/20 nmolcore changes added
!
!  On entry :
!
!  mtv         = no. of modes
!  nphonatc    = total number of cores in relevant regions
!  nphonatptr  = pointer from reduce to full atom set
!  ncfoc       = number of condensed core sites
!  iocptr      = pointer that connects sites to condensed sites
!  eigr        = eigenvectors of dynamical matrix - real part
!  eigi        = eigenvectors of dynamical matrix - imaginary part
!  maxd2       = left-hand dimension of eigr
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
  use configurations, only : nsuperghost
  use control,        only : lghost, lrigid
  use current
  use element
  use frequencies
  use iochannels
  use molecule
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
  integer(i4),                   intent(in)    :: nphonatc
  integer(i4),                   intent(in)    :: nphonatptr(*)
  integer(i4),                   intent(in)    :: ncfoc
  real(dp),                      intent(in)    :: eigr(maxd2,mtv)
  real(dp),                      intent(in)    :: eigi(maxd2,mtv)
  real(dp),                      intent(out)   :: oscstrength(3,3,mtv)
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ind
  integer(i4)                                  :: iocj
  integer(i4)                                  :: ix
  integer(i4)                                  :: iy
  integer(i4)                                  :: iz
  integer(i4)                                  :: j
  integer(i4)                                  :: k
  integer(i4)                                  :: l
  integer(i4)                                  :: m
  integer(i4)                                  :: nghostcell
  real(dp)                                     :: drQ(3,3)
  real(dp)                                     :: reigi(3)
  real(dp)                                     :: reigr(3)
  real(dp)                                     :: osxi
  real(dp)                                     :: osyi
  real(dp)                                     :: oszi
  real(dp)                                     :: osxr
  real(dp)                                     :: osyr
  real(dp)                                     :: oszr
  real(dp)                                     :: trmj
#ifdef TRACE
  call trace_in('oscillatorstrength')
#endif
!**********************************
!  Oscillator strengths per mode  *
!**********************************
!
!  Find number of ghost cells
!
  if (lghost) then
    nghostcell = nsuperghost(1,ncf)*nsuperghost(2,ncf)*nsuperghost(3,ncf)
  else
    nghostcell = 1
  endif
  if (lrigid) then
!
!  Rigid molecule version of algorithm
!
!  Loop over modes
!
    do m = 1,mtv
      osxr = 0.0_dp
      osyr = 0.0_dp
      oszr = 0.0_dp
      osxi = 0.0_dp
      osyi = 0.0_dp
      oszi = 0.0_dp
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
        osxr = osxr + (eigr(ix,m)*bornq(1,1,j)*rfmass(ix) + &
                       eigr(iy,m)*bornq(2,1,j)*rfmass(iy) + &
                       eigr(iz,m)*bornq(3,1,j)*rfmass(iz))*occuf(j)
        osyr = osyr + (eigr(ix,m)*bornq(1,2,j)*rfmass(ix) + &
                       eigr(iy,m)*bornq(2,2,j)*rfmass(iy) + &
                       eigr(iz,m)*bornq(3,2,j)*rfmass(iz))*occuf(j)
        oszr = oszr + (eigr(ix,m)*bornq(1,3,j)*rfmass(ix) + &
                       eigr(iy,m)*bornq(2,3,j)*rfmass(iy) + &
                       eigr(iz,m)*bornq(3,3,j)*rfmass(iz))*occuf(j)
!
        osxi = osxi + (eigi(ix,m)*bornq(1,1,j)*rfmass(ix) + &
                       eigi(iy,m)*bornq(2,1,j)*rfmass(iy) + &
                       eigi(iz,m)*bornq(3,1,j)*rfmass(iz))*occuf(j)
        osyi = osyi + (eigi(ix,m)*bornq(1,2,j)*rfmass(ix) + &
                       eigi(iy,m)*bornq(2,2,j)*rfmass(iy) + &
                       eigi(iz,m)*bornq(3,2,j)*rfmass(iz))*occuf(j)
        oszi = oszi + (eigi(ix,m)*bornq(1,3,j)*rfmass(ix) + &
                       eigi(iy,m)*bornq(2,3,j)*rfmass(iy) + &
                       eigi(iz,m)*bornq(3,3,j)*rfmass(iz))*occuf(j)
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
          osxr = osxr + (eigr(ind+1,m)*bornq(1,1,k)*rfmass(ind+1) + &
                         eigr(ind+2,m)*bornq(2,1,k)*rfmass(ind+2) + &
                         eigr(ind+3,m)*bornq(3,1,k)*rfmass(ind+3))*occuf(k)
          osyr = osyr + (eigr(ind+1,m)*bornq(1,2,k)*rfmass(ind+1) + &
                         eigr(ind+2,m)*bornq(2,2,k)*rfmass(ind+2) + &
                         eigr(ind+3,m)*bornq(3,2,k)*rfmass(ind+3))*occuf(k)
          oszr = oszr + (eigr(ind+1,m)*bornq(1,3,k)*rfmass(ind+1) + &
                         eigr(ind+2,m)*bornq(2,3,k)*rfmass(ind+2) + &
                         eigr(ind+3,m)*bornq(3,3,k)*rfmass(ind+3))*occuf(k)
!
          osxi = osxi + (eigi(ind+1,m)*bornq(1,1,k)*rfmass(ind+1) + &
                         eigi(ind+2,m)*bornq(2,1,k)*rfmass(ind+2) + &
                         eigi(ind+3,m)*bornq(3,1,k)*rfmass(ind+3))*occuf(k)
          osyi = osyi + (eigi(ind+1,m)*bornq(1,2,k)*rfmass(ind+1) + &
                         eigi(ind+2,m)*bornq(2,2,k)*rfmass(ind+2) + &
                         eigi(ind+3,m)*bornq(3,2,k)*rfmass(ind+3))*occuf(k)
          oszi = oszi + (eigi(ind+1,m)*bornq(1,3,k)*rfmass(ind+1) + &
                         eigi(ind+2,m)*bornq(2,3,k)*rfmass(ind+2) + &
                         eigi(ind+3,m)*bornq(3,3,k)*rfmass(ind+3))*occuf(k)
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
          reigr(1:3) = 0.0_dp
          reigi(1:3) = 0.0_dp
          do l = 1,3
            reigr(l) = reigr(l) + drQ(l,1)*eigr(ind+4,m)*rfmass(ind+4)
            reigr(l) = reigr(l) + drQ(l,2)*eigr(ind+5,m)*rfmass(ind+5)
            reigr(l) = reigr(l) + drQ(l,3)*eigr(ind+6,m)*rfmass(ind+6)
!
            reigi(l) = reigi(l) + drQ(l,1)*eigi(ind+4,m)*rfmass(ind+4)
            reigi(l) = reigi(l) + drQ(l,2)*eigi(ind+5,m)*rfmass(ind+5)
            reigi(l) = reigi(l) + drQ(l,3)*eigi(ind+6,m)*rfmass(ind+6)
          enddo
!
          osxr = osxr + (reigr(1)*bornq(1,1,k) + &
                         reigr(2)*bornq(2,1,k) + &
                         reigr(3)*bornq(3,1,k))*occuf(k)
          osyr = osyr + (reigr(1)*bornq(1,2,k) + &
                         reigr(2)*bornq(2,2,k) + &
                         reigr(3)*bornq(3,2,k))*occuf(k)
          oszr = oszr + (reigr(1)*bornq(1,3,k) + &
                         reigr(2)*bornq(2,3,k) + &
                         reigr(3)*bornq(3,3,k))*occuf(k)
!
          osxi = osxi + (reigi(1)*bornq(1,1,k) + &
                         reigi(2)*bornq(2,1,k) + &
                         reigi(3)*bornq(3,1,k))*occuf(k)
          osyi = osyi + (reigi(1)*bornq(1,2,k) + &
                         reigi(2)*bornq(2,2,k) + &
                         reigi(3)*bornq(3,2,k))*occuf(k)
          oszi = oszi + (reigi(1)*bornq(1,3,k) + &
                         reigi(2)*bornq(2,3,k) + &
                         reigi(3)*bornq(3,3,k))*occuf(k)
        enddo
      enddo
!
      oscstrength(1,1,m) = osxr*osxr + osxi*osxi
      oscstrength(2,1,m) = osyr*osxr + osyi*osxi
      oscstrength(3,1,m) = oszr*osxr + oszi*osxi
      oscstrength(1,2,m) = osxr*osyr + osxi*osyi
      oscstrength(2,2,m) = osyr*osyr + osyi*osyi
      oscstrength(3,2,m) = oszr*osyr + oszi*osyi
      oscstrength(1,3,m) = osxr*oszr + osxi*oszi
      oscstrength(2,3,m) = osyr*oszr + osyi*oszi
      oscstrength(3,3,m) = oszr*oszr + oszi*oszi
    enddo
  else
!
!  Loop over modes
!
    do m = 1,mtv
      osxr = 0.0_dp
      osyr = 0.0_dp
      oszr = 0.0_dp
      osxi = 0.0_dp
      osyi = 0.0_dp
      oszi = 0.0_dp
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
        do j = 1,nphonatc,nghostcell
          iocj = (iocptr(j) - 1)/nghostcell + 1
          if (iocj.eq.i) then
            trmj = occuf(j)/sqrt(massspec(nspecptr(nrelf2a(nphonatptr(j)))))
!
!  Multiple inverse mass weighted eigenvectors by Born charges
!
            osxr = osxr + (eigr(ix,m)*bornq(1,1,j) + &
                         eigr(iy,m)*bornq(2,1,j) + &
                         eigr(iz,m)*bornq(3,1,j))*trmj
            osyr = osyr + (eigr(ix,m)*bornq(1,2,j) + &
                         eigr(iy,m)*bornq(2,2,j) + &
                         eigr(iz,m)*bornq(3,2,j))*trmj
            oszr = oszr + (eigr(ix,m)*bornq(1,3,j) + &
                         eigr(iy,m)*bornq(2,3,j) + &
                         eigr(iz,m)*bornq(3,3,j))*trmj
!
            osxi = osxi + (eigi(ix,m)*bornq(1,1,j) + &
                         eigi(iy,m)*bornq(2,1,j) + &
                         eigi(iz,m)*bornq(3,1,j))*trmj
            osyi = osyi + (eigi(ix,m)*bornq(1,2,j) + &
                         eigi(iy,m)*bornq(2,2,j) + &
                         eigi(iz,m)*bornq(3,2,j))*trmj
            oszi = oszi + (eigi(ix,m)*bornq(1,3,j) + &
                         eigi(iy,m)*bornq(2,3,j) + &
                         eigi(iz,m)*bornq(3,3,j))*trmj
          endif
        enddo
      enddo
      oscstrength(1,1,m) = osxr*osxr + osxi*osxi
      oscstrength(2,1,m) = osyr*osxr + osyi*osxi
      oscstrength(3,1,m) = oszr*osxr + oszi*osxi
      oscstrength(1,2,m) = osxr*osyr + osxi*osyi
      oscstrength(2,2,m) = osyr*osyr + osyi*osyi
      oscstrength(3,2,m) = oszr*osyr + oszi*osyi
      oscstrength(1,3,m) = osxr*oszr + osxi*oszi
      oscstrength(2,3,m) = osyr*oszr + osyi*oszi
      oscstrength(3,3,m) = oszr*oszr + oszi*oszi
    enddo
  endif
#ifdef TRACE
  call trace_out('oscillatorstrength')
#endif
!
  return
  end
