  subroutine mcstrain(mode)
!
!  MC routine for cell strain. 
!
!  mode = if mode = 1, then create new trial strain
!         if mode = 2, then undo previous strain
!
!   5/07 Created
!   6/07 Call to x0tostrcentroid added for undo since x0 array is
!        modified to correct fractional coordinates in here.
!        Not needed for trial step since this is done in 
!        call to energy.
!  12/07 Unused variables removed
!   1/08 random -> GULP_random
!   5/14 Trap on strain number being 0 added
!   2/18 Trace added
!   8/18 Modified for changes to lstraincell algorithm
!   8/18 Adding 1 to strains 1-3 removed
!   2/19 x0 removed
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
!  Julian Gale, CIC, Curtin University, February 2019
!
  use control,       only : lstraincell
  use current
  use general
  use genetic,       only : iseed
  use molecule
  use montecarlo
  use parallel
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
!
!  Passed variables
!
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                  :: mode
!
!  Local variables
!
  integer(i4),                        save :: ntostrain
  real(dp),                           save :: deltastrain = 0.0_dp
!
  logical                                  :: lgeometryOK
  real(dp)                                 :: randnum
  real(dp)                                 :: GULP_random
#ifdef TRACE
  call trace_in('mcstrain')
#endif
!
  if (mode.eq.2) then
!*************************************
!  Mode 2 : Undo last displacements  *
!*************************************
    if (.not.lstraincell) strain(1:nstrains) = 0.0_dp
    strain(ntostrain) = - deltastrain/(1.0_dp + deltastrain)
    call x0strain(lgeometryOK)
    call x0tostrcentroid
  else
!*******************************
!  Mode 1 : New displacements  *
!*******************************
    if (.not.lstraincell) strain(1:nstrains) = 0.0_dp
!
!  Choose cell parameter to strain
!
    randnum = GULP_random(iseed,1_i4)
    ntostrain = nstrainable*randnum + 1_i4
    if (ntostrain.gt.nstrainable) ntostrain = nstrainable
    if (ntostrain.lt.1) ntostrain = 1
    ntostrain = nptrstrainable(ntostrain)
!
!  Find displacement to apply
!
    deltastrain = smaxmc*GULP_random(iseed,2_i4)
!
!  Apply strain to configuration array
!
    strain(ntostrain) = strain(ntostrain) + deltastrain
    call x0strain(lgeometryOK)
  endif
#ifdef TRACE
  call trace_out('mcstrain')
#endif
!
  return
  end
