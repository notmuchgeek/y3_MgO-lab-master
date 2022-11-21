  subroutine defvartocfg(n,xc)
!
!  Transfers values from xc to current configuration based on idoptindex/idopttype
!  Defect version
!
!   3/19 Created 
!   4/19 Symmetry modifications made
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
!  Julian Gale, CIC, Curtin University, April 2019
!
  use control,        only : ldsym
  use current
  use defects
  use optimisation,   only : idopt_dx, idopt_dy, idopt_dz, idopt_dradius
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4),  intent(in)                     :: n
  real(dp),     intent(in)                     :: xc(*)
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ind
  real(dp)                                     :: var
#ifdef TRACE
  call trace_in('defvartocfg')
#endif
!************************
!  Loop over variables  *
!************************
  do i = 1,n
    ind = idoptindex(i)
!
!  Symmetry adapt index
!
    if (ldsym) ind = ndsptr(ind)
!
    var = xc(i)
    if (idopttype(i).eq.idopt_dx) then
!
!  Defect x coordinate
!
      xdefe(ind) = var
    elseif (idopttype(i).eq.idopt_dy) then
!
!  Defect y coordinate
!
      ydefe(ind) = var
    elseif (idopttype(i).eq.idopt_dz) then
!
!  Defect z coordinate
!
      zdefe(ind) = var
    elseif (idopttype(i).eq.idopt_dradius) then
!
!  Defect radius
!
      radefe(ind) = var
    else
      call outerror('unknown optimisation variable type',0_i4)
      call stopnow('defvartocfg')
    endif
  enddo
#ifdef TRACE
  call trace_out('defvartocfg')
#endif
!
  return
  end
!
  subroutine defcfgtovar(n,xc)
!
!  Transfers values from current configuration to var based on idoptindex/idopttype
!  Defect version
!
!   3/19 Created 
!   4/19 Symmetry modifications made
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
!  Julian Gale, CIC, Curtin University, April 2019
!
  use control,        only : ldsym
  use current
  use defects
  use optimisation,   only : idopt_dx, idopt_dy, idopt_dz, idopt_dradius
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4),  intent(in)                     :: n
  real(dp),     intent(out)                    :: xc(*)
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ind
  real(dp)                                     :: var
#ifdef TRACE
  call trace_in('defcfgtovar')
#endif
!************************
!  Loop over variables  *
!************************
  do i = 1,n
    ind = idoptindex(i)
!
!  Symmetry adapt index
!
    if (ldsym) ind = ndsptr(ind)
!
    if (idopttype(i).eq.idopt_dx) then
!
!  Defect x coordinate
!
      var = xdefe(ind)
    elseif (idopttype(i).eq.idopt_dy) then
!
!  Defect y coordinate
!
      var = ydefe(ind)
    elseif (idopttype(i).eq.idopt_dz) then
!
!  Defect z coordinate
!
      var = zdefe(ind)
    elseif (idopttype(i).eq.idopt_dradius) then
!
!  Defect radius
!
      var = radefe(ind)
    else
      call outerror('unknown optimisation variable type',0_i4)
      call stopnow('defcfgtovar')
    endif
!
!  Set xc value
!
    xc(i) = var
  enddo
#ifdef TRACE
  call trace_out('defcfgtovar')
#endif
!
  return
  end
