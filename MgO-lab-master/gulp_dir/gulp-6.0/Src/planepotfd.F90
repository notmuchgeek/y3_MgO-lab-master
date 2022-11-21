  subroutine planepotfd(matom,eplane,lgrad1)
!
!  Subroutine for calculating the energy from the plane potentials - 1st deriv only
!  Finite difference version that focuses on derivatives of matom.
!
!  12/17 Created from planepotmd
!   2/18 Trace added
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
!  Copyright Curtin University 2018
!
!  Julian Gale, CIC, Curtin University, Feb 2018
!
  use configurations, only : nregionno
  use control
  use current
  use derivatives
  use mdlogic
  use optimisation
  use parallel
  use plane
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                      :: matom
  real(dp),    intent(inout)                   :: eplane
  logical,     intent(in)                      :: lgrad1
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: m
  integer(i4)                                  :: n  
  integer(i4)                                  :: nati
  integer(i4)                                  :: npp
  integer(i4)                                  :: nregioni
  integer(i4)                                  :: ntypi
  logical                                      :: lmatch
  real(dp)                                     :: Atrm
  real(dp)                                     :: Btrm
  real(dp)                                     :: d1trm
  real(dp)                                     :: deltaz
  real(dp)                                     :: etrm
  real(dp)                                     :: oci      
  real(dp)                                     :: rdeltaz
  real(dp)                                     :: zp
#ifdef TRACE
  call trace_in('planepotfd')
#endif
!
!  Initialise the energy
!
  eplane = 0.0_dp
! 
!  Only do case of i = matom
! 
  i = matom
!
!  Set attributes of i
!
  nati = nat(i)
  ntypi = nftype(i)
  oci = occuf(i)
  nregioni = nregionno(nsft+i)
!
!  Loop over plane potentials
!
  etrm = 0.0_dp
  do npp = 1,nplanepot
!
!  Does this potential apply to this atom
!
    if (lmatch(nati,ntypi,natplanepot(npp),ntypplanepot(npp),.true.)) then
!
!  Select potential type
!
      if (nplanepottype(npp).eq.1) then
!
!  Lennard-Jones
!
        zp = planepot(1,npp)
        deltaz = zclat(i) - zp
!
!  Check distance against cutoffs
!
        if (abs(deltaz).gt.planepotrmin(npp).and.abs(deltaz).le.planepotrmax(npp)) then
          m = nplanepotpower(1,npp)
          n = nplanepotpower(2,npp)
          rdeltaz = 1.0_dp/deltaz
          Atrm = planepot(2,npp)*rdeltaz**m
          Btrm = planepot(3,npp)*rdeltaz**n
          etrm = etrm + oci*(Atrm - Btrm)
          if (lgrad1) then
            d1trm = - oci*rdeltaz*(dble(m)*Atrm - dble(n)*Btrm)
            zdrv(i) = zdrv(i) + d1trm
          endif
        endif
      endif
    endif
  enddo
!
!  Region handling
!
  eplane = eplane + etrm
#ifdef TRACE
  call trace_out('planepotfd')
#endif
!
  return
  end
