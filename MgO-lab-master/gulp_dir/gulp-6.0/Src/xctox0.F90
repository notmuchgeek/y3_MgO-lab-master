  subroutine xctox0(n,xc,lgeometryOK)
!
!  Called by funct / fefunct to set linear structure array
!  from optimisation variables array xc.
!
!   8/97 Created from funct
!  12/00 Generalised for 0 to 3-D
!   6/01 Order of cell operations corrected and lra test added
!   5/02 Check for small cell added
!  10/03 Cell parameter variables option added
!   5/04 lmodco option introduced
!  11/04 Inverse cell parameters set
!   5/07 Application of cell strain moved to subroutine
!  12/07 Unused variables removed
!   5/08 Geometry check flag added as argument
!   9/15 Added number of cells for lmodco increased from 10 to 1000
!   1/18 Trace added
!   6/18 Strain cell option added
!   8/18 Modified due to changes in lstraincell algorithm
!   8/18 Adding 1 to strains 1-3 removed
!   3/19 x0 removed
!   3/19 Constraint arrays changed to have index and type
!  10/19 Rigid molecule constraints added
!   5/20 Application of constraints now moved to subroutine
!   7/20 Checking of quaternions added for setting of lgeometryOK
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
  use configurations
  use control,        only : lmodco, lstraincell, lrigid
  use current
  use molecule
  use optimisation
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4),  intent(in)                     :: n
  logical,      intent(out)                    :: lgeometryOK
  real(dp),     intent(in)                     :: xc(*)
!
!  Local variables
!
  integer(i4)                                  :: i
  logical                                      :: lgok
#ifdef TRACE
  call trace_in('xctox0')
#endif
!
  lgeometryOK = .true.
!
!  First substitute parameters into place
!
  if (ndim.gt.0) then
    do i = 1,3
      rv(1,i) = rvcfg(1,i,ncf)
      rv(2,i) = rvcfg(2,i,ncf)
      rv(3,i) = rvcfg(3,i,ncf)
    enddo
    if (.not.loptcellpar.and..not.lstraincell) then
      strain(1:nstrains) = 0.0_dp
    endif
  endif
!
!  Transfer variables to configuration
!
  call vartocfg(n,xc)
!
!  For rigid molecules check quaternions
!
  if (lrigid) then
    call checkquaternions(lgok,.false.)
    if (.not.lgok) lgeometryOK = .false.
  endif
!
!  Make sure all fractional coords are between 0 and 1
!
  if (lmodco) then
    if (ndim.eq.3) then
      do i = 1,nasym
        xafrac(i) = mod(xafrac(i)+1000.0_dp,1.0_dp)
        yafrac(i) = mod(yafrac(i)+1000.0_dp,1.0_dp)
        zafrac(i) = mod(zafrac(i)+1000.0_dp,1.0_dp)
      enddo
    elseif (ndim.eq.2) then
      do i = 1,nasym
        xafrac(i) = mod(xafrac(i)+1000.0_dp,1.0_dp)
        yafrac(i) = mod(yafrac(i)+1000.0_dp,1.0_dp)
      enddo
    elseif (ndim.eq.1) then
      do i = 1,nasym
        xafrac(i) = mod(xafrac(i)+1000.0_dp,1.0_dp)
      enddo
    endif
  endif
!**********************
!  Apply constraints  *
!**********************
  if (ncon.gt.0) then
    call applyconstraints
  endif
  if (ndim.gt.0) then
!*************************
!  Apply strain to cell  *
!*************************
    call x0strain(lgok)
    if (.not.lgok) lgeometryOK = .false.
  endif
#ifdef TRACE
  call trace_out('xctox0')
#endif
!
  return
  end
