  subroutine findcelltype(a,b,c,alpha,beta,gamma,ictype,icfhr)
!
!   Finds the cell type based on the cell parameters
!
!   8/13 Created from cellcheck
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
  use datatypes
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4),                         intent(out)  :: icfhr        ! 0 => hexagonal, 1 => rhombohedral, if ictype = 5
  integer(i4),                         intent(out)  :: ictype       ! Cell type indicator
  real(dp),                            intent(in)   :: a
  real(dp),                            intent(in)   :: b
  real(dp),                            intent(in)   :: c
  real(dp),                            intent(in)   :: alpha
  real(dp),                            intent(in)   :: beta
  real(dp),                            intent(in)   :: gamma
!
!  Local variables
!
  logical                                           :: la90
  logical                                           :: lb90
  logical                                           :: lc90
  logical                                           :: lc120
  logical                                           :: laeqb
  logical                                           :: laeqc
  logical                                           :: laleqbe
  logical                                           :: laleqga
  logical                                           :: lcellok
#ifdef TRACE
  call trace_in('findcelltype')
#endif
!
!  Check cell parameters are consistent
!
  lcellok = .false.
  la90 = (abs(90.0_dp-alpha).lt.1.0d-6)
  lb90 = (abs(90.0_dp-beta).lt.1.0d-6)
  lc90 = (abs(90.0_dp-gamma).lt.1.0d-6)
  lc120 = (abs(120.0_dp-gamma).lt.1.0d-6)
  laeqb = (abs(a-b).lt.1.0d-6)
  laeqc = (abs(a-c).lt.1.0d-6)
  laleqbe = (abs(alpha-beta).lt.1.0d-6)
  laleqga = (abs(alpha-gamma).lt.1.0d-6)
!
!  Set default cell as triclinic
!
  ictype = 1
  icfhr  = 0
!
!  Right angled cells
!
  if (la90.and.lb90.and.lc90) then
    if (laeqb.and.laeqc) then
!
!  Cubic
!
      ictype = 6
    elseif (laeqb) then
!
!  Tetragonal
!
      ictype = 4
    else
!
!  Orthorhombic
!
      ictype = 3
    endif
  elseif ((la90.and.lb90.and.lc120).and.laeqb) then
!
!  Hexagonal
!
    ictype = 5
    icfhr = 0
  elseif ((la90.and.lb90).or.(la90.and.lc90).or.(lb90.and.lc90)) then
!
!  Monoclinic
!
    ictype = 2
  elseif (laleqbe.and.laleqga.and.laeqb.and.laeqc) then
!
!  Rhombohedral
!
    ictype = 5
    icfhr = 1
  endif
#ifdef TRACE
  call trace_out('findcelltype')
#endif
!
  return
  end
