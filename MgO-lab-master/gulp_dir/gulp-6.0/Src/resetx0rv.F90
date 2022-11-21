  subroutine resetx0rv(n,xc)
!
!  Applies a strain to the cell based on contents of x0
!  stores in rvcfg and then resets the strains to 0.
!  NB Not to be used with lstraincell
!
!  12/18 Created from part of x0strain
!   3/19 Bug in re-initialisation of x0 fixed
!   3/19 x0 removed
!   3/19 iopt replaced by ioptindex and iopttype
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
!  Julian Gale, CIC, Curtin University, March 2019
!
  use configurations
  use control,        only : lstraincell
  use current
  use optimisation,   only : loptcellpar, iopt_strain
  use symmetry,       only : lra
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4),                   intent(in)    :: n
  real(dp),                      intent(inout) :: xc(n)
!
!  Local variables
!
  integer(i4)                                  :: i
  real(dp)                                     :: sum
!
!  Trap straincell or cell parameter algorithms
!
  if (lstraincell.or.loptcellpar) return
#ifdef TRACE
  call trace_in('resetx0rv')
#endif
!
  if (ndim.gt.0) then
!******************
!  Apply strains  *
!******************
    rv(1:3,1:ndim) = rvcfg(1:3,1:ndim,ncf)
    if (ndim.eq.3) then
      call strain3D(strain,rv)
      call celltype(ictype,icfhr)
    elseif (ndim.eq.2) then
      call strain2D(strain,rv)
    elseif (ndim.eq.1) then
      call strain1D(strain,rv)
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
!  Store cell
!
  rvcfg(1:3,1:ndim,ncf) = rv(1:3,1:ndim)
!
!  Zero strains
!
  strain(1:nstrains) = 0.0_dp
  do i = 1,nvar
    if (iopttype(i).eq.iopt_strain) then
      xc(i) = 0.0_dp
    endif
  enddo
#ifdef TRACE
  call trace_out('resetx0rv')
#endif
!
  return
  end
