  subroutine initmaxspecdefaults(i)
!
!  Initialises the arrays associated with maxspec
!
!   9/10 Created from changemax routine
!   9/12 NMR properties added
!  12/12 NMR f value added
!   8/13 symspec added
!   2/15 lmm3se added
!   3/15 Gasteiger parameters added
!   1/16 Species specific VDW radius added for cosmo
!   6/19 Spin added
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
  use library
  use polarise
  use montecarlo
  use shells
  use species
  use two
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in) :: i
!
!  Initialise defaults for new part of array
!
  if (i.ge.1.and.i.le.maxspec) then
    symspec(i) = ' '
    libspec(i) = ' '
    lbrspec(i) = .false.
    ldefshspec(i) = .false.
    lgastinlibspec(i) = .false.
    lgastinspec(i) = .false.
    linspec(i) = .false.
    lmassinspec(i) = .false.
    lnmrinspec(i) = .false.
    lqinspec(i) = .false.
    lmask(i) = .false.
    lmm3se(i) = .false.
    lspininspec(i) = .false.
    lvdwinspec(i) = .false.
    gastspec(1:4,i) = 0.0_dp
    massspec(i) = 0.0_dp
    nmrspec(1,i) = 0.0_dp
    nmrspec(2,i) = 0.0_dp
    nmrspec(3,i) = 1.0_dp
    radspec(i) = 0.0_dp
    spinspec(i) = 0.0_dp
    vdwspec(i) = 0.0_dp
    natpolspec(i) = 0
    natratiomspec(i) = 0
    dpolmaxspec(i) = 0.0_dp
    dpolspec(i) = 0.0_dp
    qpolspec(i) = 0.0_dp
    ratiomspec(i) = 0.0_dp
  endif
!
  return
  end
