  subroutine setpolar
!
!  Set up polarisability for each atom in structure
!
!   5/00 created
!   4/19 Analytic second derivatives now allowed
!   8/19 Set cutoff for polarisation damping
!  10/19 Langevin damping of dipoles added
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
!  Julian Gale, CIC, Curtin University, October 2019
!
  use control
  use current
  use element, only : maxele
  use parallel
  use polarise
  implicit none
!
!  Local variables
!
  integer(i4)      :: i
  integer(i4)      :: j
  integer(i4)      :: natp
  integer(i4)      :: ntypp
  real(dp)         :: dpoli
  real(dp)         :: dpolmaxi
  real(dp)         :: qpoli
!
!  Set flag indicating that third derivatives are not available with polarisability at the moment.
!
  lnoanald3 = .true.
!
!  Initialise polarisabilities to zero 
!
  do i = 1,nasym
    dpolar(i) = 0.0_dp
    qpolar(i) = 0.0_dp
    dpolarmax(i) = 0.0_dp
  enddo
!
!  Loop over polarisability species assigning to atoms
!
  do i = 1,npolspec
    natp = natpolspec(i)
    ntypp = ntyppolspec(i)
    dpoli = dpolspec(i)
    qpoli = qpolspec(i)
    dpolmaxi = dpolmaxspec(i)
    do j = 1,nasym
      if (iatn(j).le.maxele) then
!
!  Cores only
!
        if (iatn(j).eq.natp) then
          if (ntypp.eq.0.or.ntypp.eq.natype(j)) then
            dpolar(j) = dpolar(j) + dpoli
            qpolar(j) = qpolar(j) + qpoli
            dpolarmax(j) = dpolarmax(j) + dpolmaxi
          endif
        endif
      endif
    enddo
  enddo
!
  if (lpoldamp) then
    rpdamp = 20.0_dp
  endif
!
  return
  end
