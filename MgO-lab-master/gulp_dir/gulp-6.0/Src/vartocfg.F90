  subroutine vartocfg(n,xc)
!
!  Transfers values from xc to current configuration based on ioptindex/iopttype
!
!   3/19 Created 
!  10/19 Rigid molecules added
!  12/19 Correction to referencing of internal derivatives for rigid molecule case
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
!  Julian Gale, CIC, Curtin University, December 2019
!
  use current
  use molecule,       only : molcoma, molQa
  use optimisation,   only : iopt_cell, iopt_strain
  use optimisation,   only : iopt_xf, iopt_yf, iopt_zf, iopt_radius
  use optimisation,   only : iopt_xcom, iopt_ycom, iopt_zcom
  use optimisation,   only : iopt_xqtn, iopt_yqtn, iopt_zqtn
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
  call trace_in('vartocfg')
#endif
!************************
!  Loop over variables  *
!************************
  do i = 1,n
    ind = ioptindex(i)
    var = xc(i)
    if (iopttype(i).eq.iopt_cell) then
!
!  Cell parameter
!
      if (ind.eq.1) then
        a = var
      elseif (ind.eq.2) then
        b = var
      elseif (ind.eq.3) then
        c = var
      elseif (ind.eq.4) then
        alpha = var
      elseif (ind.eq.5) then
        beta = var
      elseif (ind.eq.6) then
        gamma = var
      endif
    elseif (iopttype(i).eq.iopt_strain) then
!
!  Strain
!
      strain(ind) = var
    elseif (iopttype(i).eq.iopt_xf) then
!
!  x fractional coordinate
!
      xafrac(nasymnomolptr(ind)) = var
    elseif (iopttype(i).eq.iopt_yf) then
!
!  y fractional coordinate
!
      yafrac(nasymnomolptr(ind)) = var
    elseif (iopttype(i).eq.iopt_zf) then
!
!  z fractional coordinate
!
      zafrac(nasymnomolptr(ind)) = var
    elseif (iopttype(i).eq.iopt_radius) then
!
!  Radius
!
      rada(nasymnomolptr(ind)) = var
    elseif (iopttype(i).eq.iopt_xcom) then
!
!  x molecule centre of mass
!
      molcoma(1,ind) = var
    elseif (iopttype(i).eq.iopt_ycom) then
!
!  y molecule centre of mass
!
      molcoma(2,ind) = var
    elseif (iopttype(i).eq.iopt_zcom) then
!
!  z molecule centre of mass
!
      molcoma(3,ind) = var
    elseif (iopttype(i).eq.iopt_xqtn) then
!
!  x molecule quaternion
!
      molQa(1,ind) = var
    elseif (iopttype(i).eq.iopt_yqtn) then
!
!  y molecule quaternion
!
      molQa(2,ind) = var
    elseif (iopttype(i).eq.iopt_zqtn) then
!
!  z molecule quaternion
!
      molQa(3,ind) = var
    else
      call outerror('unknown optimisation variable type',0_i4)
      call stopnow('vartocfg')
    endif
  enddo
#ifdef TRACE
  call trace_out('vartocfg')
#endif
!
  return
  end
!
  subroutine cfgtovar(n,xc)
!
!  Transfers values from current configuration to var based on ioptindex/iopttype
!
!  iopttype indices:
!
!    1 = cell parameter
!    2 = strain component
!    3 = x fractional coordinate
!    4 = y fractional coordinate
!    5 = z fractional coordinate
!    6 = radius
!    7 = x molecule centre of mass
!    8 = y molecule centre of mass
!    9 = z molecule centre of mass
!   10 = x molecule quaternion
!   11 = y molecule quaternion
!   12 = z molecule quaternion
!
!   3/19 Created 
!  10/19 Rigid molecules added
!  12/19 Correction to referencing of internal derivatives for rigid molecule case
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
!  Julian Gale, CIC, Curtin University, December 2019
!
  use current
  use molecule,       only : molcoma, molQa
  use optimisation,   only : iopt_cell, iopt_strain
  use optimisation,   only : iopt_xf, iopt_yf, iopt_zf, iopt_radius
  use optimisation,   only : iopt_xcom, iopt_ycom, iopt_zcom
  use optimisation,   only : iopt_xqtn, iopt_yqtn, iopt_zqtn
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
  call trace_in('cfgtovar')
#endif
!************************
!  Loop over variables  *
!************************
  do i = 1,n
    ind = ioptindex(i)
    if (iopttype(i).eq.iopt_cell) then
!
!  Cell parameter
!
      if (ind.eq.1) then
        var = a
      elseif (ind.eq.2) then
        var = b
      elseif (ind.eq.3) then
        var = c
      elseif (ind.eq.4) then
        var = alpha
      elseif (ind.eq.5) then
        var = beta
      elseif (ind.eq.6) then
        var = gamma
      endif
    elseif (iopttype(i).eq.iopt_strain) then
!
!  Strain
!
      var = strain(ind)
    elseif (iopttype(i).eq.iopt_xf) then
!
!  x fractional coordinate
!
      var = xafrac(nasymnomolptr(ind))
    elseif (iopttype(i).eq.iopt_yf) then
!
!  y fractional coordinate
!
      var = yafrac(nasymnomolptr(ind))
    elseif (iopttype(i).eq.iopt_zf) then
!
!  z fractional coordinate
!
      var = zafrac(nasymnomolptr(ind))
    elseif (iopttype(i).eq.iopt_radius) then
!
!  Radius
!
      var = rada(nasymnomolptr(ind))
    elseif (iopttype(i).eq.iopt_xcom) then
!
!  x molecule centre of mass
!
      var = molcoma(1,ind)
    elseif (iopttype(i).eq.iopt_ycom) then
!
!  y molecule centre of mass
!
      var = molcoma(2,ind)
    elseif (iopttype(i).eq.iopt_zcom) then
!
!  z molecule centre of mass
!
      var = molcoma(3,ind)
    elseif (iopttype(i).eq.iopt_xqtn) then
!
!  x molecule quaternion
!
      var = molQa(1,ind)
    elseif (iopttype(i).eq.iopt_yqtn) then
!
!  y molecule quaternion
!
      var = molQa(2,ind)
    elseif (iopttype(i).eq.iopt_zqtn) then
!
!  z molecule quaternion
!
      var = molQa(3,ind)
    else
      call outerror('unknown optimisation variable type',0_i4)
      call stopnow('cfgtovar')
    endif
!
!  Set xc value
!
    xc(i) = var
  enddo
#ifdef TRACE
  call trace_out('cfgtovar')
#endif
!
  return
  end
