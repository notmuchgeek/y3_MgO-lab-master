  subroutine x0strain(lgeometryOK)
!
!  Applies a strain to the cell based on contents of x0.
!
!   5/07 Created from part of xctox0
!   5/08 Flag added to indicate whether the geometry is OK
!        to replace internal error trap on cell parameters
!   1/18 Trace added
!   8/18 Adding 1 to strains 1-3 removed
!   9/18 Strain module added
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
  use configurations
  use control,        only : lstraincell
  use current
  use general,        only : cellmin
  use m_strain,       only : getinversestrain
  use optimisation,   only : loptcellpar
  use symmetry,       only : lra
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  logical,                        intent(out)  :: lgeometryOK
!
!  Local variables
!
  integer(i4)                                  :: i
  real(dp)                                     :: sum
#ifdef TRACE
  call trace_in('x0strain')
#endif
!
  if (ndim.gt.0) then
    if (loptcellpar) then
!**************************
!  Generate cell vectors  *
!**************************
      if (ndim.eq.3) then
        call cell3D(rv,a,b,c,alpha,beta,gamma)
        call celltype(ictype,icfhr)
      elseif (ndim.eq.2) then
        call cell2D(rv,a,b,alpha)
      elseif (ndim.eq.1) then
        call cell1D(rv,a)
      endif
    else
!******************
!  Apply strains  *
!******************
      if (lstraincell) then
        do i = 1,ndim
          rv(1,i) = rvcfg(1,i,ncf)
          rv(2,i) = rvcfg(2,i,ncf)
          rv(3,i) = rvcfg(3,i,ncf)
        enddo
      endif
      if (ndim.eq.3) then
        call strain3D(strain,rv)
        call celltype(ictype,icfhr)
      elseif (ndim.eq.2) then
        call strain2D(strain,rv)
      elseif (ndim.eq.1) then
        call strain1D(strain,rv)
      endif
    endif
    if (lstraincell) then
!
!  Generate inverse strain matrix for lstraincell algorithm
!
      call getinversestrain(strain)
    endif
  endif
!
!  Set parameters related to cell
!
  if (ndim.eq.3) then
    r1x = rv(1,1)
    r1y = rv(2,1)
    r1z = rv(3,1)
    r2x = rv(1,2)
    r2y = rv(2,2)
    r2z = rv(3,2)
    r3x = rv(1,3)
    r3y = rv(2,3)
    r3z = rv(3,3)
    call uncell3D(rv,a,b,c,alpha,beta,gamma)
    if (a.gt.1.0d-12) then
      recipa = 1.0_dp/a
    else
      recipa = 0.0_dp
    endif
    if (b.gt.1.0d-12) then
      recipb = 1.0_dp/b
    else
      recipb = 0.0_dp
    endif
    if (c.gt.1.0d-12) then
      recipc = 1.0_dp/c
    else
      recipc = 0.0_dp
    endif
    sum = abs(r2x) + abs(r1y) + abs(r3x) + abs(r1z) + abs(r3y) + abs(r2z)
    lra = (sum.lt.1.0d-6)
    call rlist
  elseif (ndim.eq.2) then
    r1x = rv(1,1)
    r1y = rv(2,1)
    r2x = rv(1,2)
    r2y = rv(2,2)
    call uncell2D(rv,a,b,alpha)
    if (a.gt.1.0d-12) then
      recipa = 1.0_dp/a
    else
      recipa = 0.0_dp
    endif
    if (b.gt.1.0d-12) then
      recipb = 1.0_dp/b
    else
      recipb = 0.0_dp
    endif
    lra = (abs(alpha-90.0_dp).lt.1.0d-6)
    call rlist
  elseif (ndim.eq.1) then
    r1x = rv(1,1)
    call uncell1D(rv,a)
    if (a.gt.1.0d-12) then
      recipa = 1.0_dp/a
    else
      recipa = 0.0_dp
    endif
    lra = .true.
    call rlist
  endif
!
!  Check for extreme cell parameters
!
  lgeometryOK = .true.
  if (ndim.eq.3) then
    if (abs(a).lt.cellmin.or.abs(b).lt.cellmin.or.abs(c).lt.cellmin) then
      call outerror('Cell parameter has fallen below allowed limit',0_i4)
      lgeometryOK = .false.
    endif
  elseif (ndim.eq.2) then
    if (abs(a).lt.cellmin.or.abs(b).lt.cellmin) then
      call outerror('Cell parameter has fallen below allowed limit',0_i4)
      lgeometryOK = .false.
    endif
  elseif (ndim.eq.1) then
    if (abs(a).lt.cellmin) then
      call outerror('Cell parameter has fallen below allowed limit',0_i4)
      lgeometryOK = .false.
    endif
  endif
#ifdef TRACE
  call trace_out('x0strain')
#endif
!
  return
  end
!
  subroutine x0strainf(x0str,lgeometryOK)
!
!  Applies a strain to the cell based on contents of x0str.
!
!   9/17 Created from x0strain
!   8/18 Modified due to changes in lstraincell algorithm
!   8/18 Adding 1 to strains 1-3 removed
!   7/19 Resetting of rv for straincell removed since this now
!        done in the calling routine
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
!  Julian Gale, CIC, Curtin University, July 2019
!
  use configurations
  use control,        only : lstraincell
  use current
  use general,        only : cellmin
  use m_strain,       only : getinversestrain
  use symmetry,       only : lra
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  logical,                        intent(out)  :: lgeometryOK
  real(dp),                       intent(in)   :: x0str(6)
!
!  Local variables
!
  real(dp)                                     :: sum
#ifdef TRACE
  call trace_in('x0strainf')
#endif
!
  if (ndim.gt.0) then
!******************
!  Apply strains  *
!******************
    if (ndim.eq.3) then
      call strain3D(x0str,rv)
      call celltype(ictype,icfhr)
    elseif (ndim.eq.2) then
      call strain2D(x0str,rv)
    elseif (ndim.eq.1) then
      call strain1D(x0str,rv)
    endif
    if (lstraincell) then
!
!  Generate inverse strain matrix for lstraincell algorithm
!
      call getinversestrain(x0str)
    endif
  endif
!
!  Set parameters related to cell
!
  if (ndim.eq.3) then
    r1x = rv(1,1)
    r1y = rv(2,1)
    r1z = rv(3,1)
    r2x = rv(1,2)
    r2y = rv(2,2)
    r2z = rv(3,2)
    r3x = rv(1,3)
    r3y = rv(2,3)
    r3z = rv(3,3)
    call uncell3D(rv,a,b,c,alpha,beta,gamma)
    if (a.gt.1.0d-12) then
      recipa = 1.0_dp/a
    else
      recipa = 0.0_dp
    endif
    if (b.gt.1.0d-12) then
      recipb = 1.0_dp/b
    else
      recipb = 0.0_dp
    endif
    if (c.gt.1.0d-12) then
      recipc = 1.0_dp/c
    else
      recipc = 0.0_dp
    endif
    sum = abs(r2x) + abs(r1y) + abs(r3x) + abs(r1z) + abs(r3y) + abs(r2z)
    lra = (sum.lt.1.0d-6)
    call rlist
  elseif (ndim.eq.2) then
    r1x = rv(1,1)
    r1y = rv(2,1)
    r2x = rv(1,2)
    r2y = rv(2,2)
    call uncell2D(rv,a,b,alpha)
    if (a.gt.1.0d-12) then
      recipa = 1.0_dp/a
    else
      recipa = 0.0_dp
    endif
    if (b.gt.1.0d-12) then
      recipb = 1.0_dp/b
    else
      recipb = 0.0_dp
    endif
    lra = (abs(alpha-90.0_dp).lt.1.0d-6)
    call rlist
  elseif (ndim.eq.1) then
    r1x = rv(1,1)
    call uncell1D(rv,a)
    if (a.gt.1.0d-12) then
      recipa = 1.0_dp/a
    else
      recipa = 0.0_dp
    endif
    lra = .true.
    call rlist
  endif
!
!  Check for extreme cell parameters
!
  lgeometryOK = .true.
  if (ndim.eq.3) then
    if (abs(a).lt.cellmin.or.abs(b).lt.cellmin.or.abs(c).lt.cellmin) then
      call outerror('Cell parameter has fallen below allowed limit',0_i4)
      lgeometryOK = .false.
    endif
  elseif (ndim.eq.2) then
    if (abs(a).lt.cellmin.or.abs(b).lt.cellmin) then
      call outerror('Cell parameter has fallen below allowed limit',0_i4)
      lgeometryOK = .false.
    endif
  elseif (ndim.eq.1) then
    if (abs(a).lt.cellmin) then
      call outerror('Cell parameter has fallen below allowed limit',0_i4)
      lgeometryOK = .false.
    endif
  endif
#ifdef TRACE
  call trace_out('x0strainf')
#endif
!
  return
  end
