  subroutine d1charged(iloc,i,lopi,nor,d0i)
!
!  Calculates the contribution to the first derivatives from the bond
!  order charge derivatives. 
!  Distibuted memory parallel version
!
!   7/17 Created from d1charge
!   2/18 Trace added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
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
!  Copyright Curtin University 2019
!
!  Julian Gale, CIC, Curtin University, September 2019
!
  use configurations, only : nregionno
  use control,        only : latomicstress
  use current,        only : nstrains, nrelf2a, nsft
  use derivatives
  use optimisation,   only : lfreeze, lopf
  use symmetry,       only : lstr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in) :: i
  integer(i4), intent(in) :: iloc
  integer(i4), intent(in) :: nor
  logical,     intent(in) :: lopi
  real(dp),    intent(in) :: d0i(*)
!
!  Local variables
!
  integer(i4)             :: ix
  integer(i4)             :: iy
  integer(i4)             :: iz
  integer(i4)             :: k
  integer(i4)             :: kl
  integer(i4)             :: kx
  integer(i4)             :: ky
  integer(i4)             :: kz
  integer(i4)             :: n
  integer(i4)             :: nregioni
  integer(i4)             :: nregionk
  logical                 :: lopk
  real(dp)                :: qpotsumi
#ifdef TRACE
  call trace_in('d1charged')
#endif
!
!  Loop over distances collecting total Coulomb potential
!
  qpotsumi = 0.0_dp
  do n = 1,nor
    qpotsumi = qpotsumi + d0i(n)
  enddo
!
!  Derivatives for atom i
!
  if (lopi) then
    ix = 3*(i - 1) + 1
    iy = ix + 1
    iz = iy + 1
    xdrv(i) = xdrv(i) + qpotsumi*dqdxyz(ix,iloc)
    ydrv(i) = ydrv(i) + qpotsumi*dqdxyz(iy,iloc)
    zdrv(i) = zdrv(i) + qpotsumi*dqdxyz(iz,iloc)
  endif
!
  nregioni = nregionno(nsft+nrelf2a(i))
  xregdrv(nregioni) = xregdrv(nregioni) + qpotsumi*dqdxyz(ix,iloc)
  yregdrv(nregioni) = yregdrv(nregioni) + qpotsumi*dqdxyz(iy,iloc)
  zregdrv(nregioni) = zregdrv(nregioni) + qpotsumi*dqdxyz(iz,iloc)
!
  do n = 1,nqatoms(iloc)
    k = nqatomptr(n,iloc)
    kx = 3*(k - 1) + 1
    ky = kx + 1
    kz = ky + 1
    lopk = (.not.lfreeze.or.lopf(nrelf2a(k)))
    if (lopk) then
      xdrv(k) = xdrv(k) + qpotsumi*dqdxyz(kx,iloc)
      ydrv(k) = ydrv(k) + qpotsumi*dqdxyz(ky,iloc)
      zdrv(k) = zdrv(k) + qpotsumi*dqdxyz(kz,iloc)
    endif
    nregionk = nregionno(nsft+nrelf2a(k))
    xregdrv(nregionk) = xregdrv(nregionk) + qpotsumi*dqdxyz(kx,iloc)
    yregdrv(nregionk) = yregdrv(nregionk) + qpotsumi*dqdxyz(ky,iloc)
    zregdrv(nregionk) = zregdrv(nregionk) + qpotsumi*dqdxyz(kz,iloc)
  enddo
!
!  Strain derivatives
!
  if (lstr) then
    do kl = 1,nstrains
      rstrd(kl) = rstrd(kl) + qpotsumi*dqds(kl,iloc)
    enddo
    if (latomicstress) then
      do kl = 1,nstrains
        atomicstress(kl,i) = atomicstress(kl,i) + qpotsumi*dqds(kl,iloc)
      enddo
    endif
  endif
#ifdef TRACE
  call trace_out('d1charged')
#endif
!
  return
  end
