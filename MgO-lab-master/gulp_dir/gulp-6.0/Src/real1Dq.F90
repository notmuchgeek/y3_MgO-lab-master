  subroutine real1Dq(ereal,qtot,lgrad1,lgrad2)
!
!  This subroutine calculates the electrostatic self-interaction energy of a 1-D system.
!
!   5/02 Created from real1D
!   5/02 hfunc call modified for third derivatives
!   3/09 Explicit 1.0d-15 replaced by global value smallself from general module
!   5/12 Atomic stresses added
!   2/18 Trace added
!   5/19 sderv2 contributions corrected for change in strfin
!   3/20 Location of angstoev changed to current
!   3/20 Call to emfunc change due to creation of separate emfuncs for strains
!   4/20 Finite strain corrections added
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
!  No claim is made that this program is free from errors and
!  no liability will be accepted for any loss or damage that
!  may result. The user is responsible for checking the validity
!  of their results.
!
!  Copyright Curtin University 2020
!
!  Julian Gale, CIC, Curtin University, April 2020
!
  use control,       only : lnoreal, latomicstress
  use current
  use derivatives
  use general,       only : nemorder, smallself
  use parallel
  use qmedata,       only : maxloop
  use symmetry,      only : lstr
  use times
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed arguments
!
  real(dp) :: ereal
  real(dp) :: qtot
  logical  :: lgrad1
  logical  :: lgrad2
!
!  Local variables
!
  integer   :: i
  integer   :: m
  real(dp)  :: acell
  real(dp)  :: assum
  real(dp)  :: g_cpu_time
  real(dp)  :: d0
  real(dp)  :: d1
  real(dp)  :: d1s
  real(dp)  :: d2s
  real(dp)  :: d2
  real(dp)  :: dads
  real(dp)  :: d2ads2
  real(dp)  :: dh1(3)
  real(dp)  :: dh1s
  real(dp)  :: d2h1(6)
  real(dp)  :: d2h1m(3)
  real(dp)  :: d2h1s
  real(dp)  :: d3h1(10)
  real(dp)  :: d3h1m(6)
  real(dp)  :: esum
  real(dp)  :: esumem
  real(dp)  :: esumh
  real(dp)  :: e1
  real(dp)  :: h1
  real(dp)  :: lna
  real(dp)  :: qi
  real(dp)  :: r
  real(dp)  :: rr
  real(dp)  :: t1, t2
  real(dp)  :: u
!
!  If noreal specified, return
!
  if (lnoreal) then
    ereal = 0.0_dp
    return
  endif
#ifdef TRACE
  call trace_in('real1Dq')
#endif
!
  t1 = g_cpu_time()
  if (lgrad1) then
!
!  Set up strain derivatives of cell parameter 
!
    if (lfinitestrain) then
      dads = a/(1.0_dp + strain(1))
      if (lgrad2) d2ads2 = 0.0_dp
    else
      dads = a
      if (lgrad2) d2ads2 = a
    endif
  endif
!********************************************************
!  Calculate Coulomb sum converged to desired accuracy  *
!********************************************************
!
!  Loop over number of cells in sum
!
  esum = 0.0_dp
  assum = 0.0_dp
  do m = 1,maxloop(1)
!
!  Direct sum component over neutral cells
!
    acell = dble(m)*a
    r = acell*acell
    if (r.gt.smallself) then
      r = sqrt(r)
      rr = 1.0_dp/r
      d0 = qtot*rr
      esum = esum - d0
      if (lgrad1.and.lstr) then
        d1 = qtot*angstoev*rr*rr*dble(m)*dads
        rstrd(1) = rstrd(1) + d1
        assum = assum + rr
        if (lgrad2) then
          d2 = qtot*angstoev*rr*rr*(-(2.0_dp*rr)*(dble(m)*dads)**2 + dble(m)*d2ads2)
          sderv2(1,1) = sderv2(1,1) + d2
        endif
      endif
    endif
  enddo
!
!  Neutralising terms
!
!  Background
!
!  and
!
!  Euler-MacLaurin component
!
  lna = log(a)
  u = (dble(maxloop(1))+0.5_dp)*a
  if (maxloop(1).gt.0) then
    call hfunc(u,0.0_dp,1.0_dp,0.0_dp,0.0_dp,h1,dh1,d2h1,d3h1,.false.,.false.,.false.)
    esumh = qtot*(h1 - lna)/a
    call emfunc(nemorder,u,0.0_dp,1.0_dp,0.0_dp,0.0_dp,a,e1,dh1,d2h1,d3h1,.false.,.false.,.false.)
    esumem = - qtot*e1
  else
    esumh = 0.0_dp
    esumem = 0.0_dp
  endif
!
!  Sum up terms
!
  ereal = esum + esumh + esumem
  ereal = ereal*angstoev
!
!  Derivatives of Euler-MacLaurin terms since these are not cumulative
!
  if (lgrad1.and.lstr.and.maxloop(1).gt.0) then
    call hfuncs(u,0.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp,h1,dh1,dh1s,d2h1,d2h1s,d2h1m, &
                    d3h1,lgrad1,lgrad2,.false.)
    d1 = qtot*angstoev
    d1s = (h1 - lna)*dads/a**2 - (dh1s - dads/a)/a
    assum = assum - d1s
    d1s = d1s*d1
    rstrd(1) = rstrd(1) - d1s
    if (lgrad2) then
      d2s = - 2.0_dp*(h1 - lna)*dads*dads/a**3 + 2.0_dp*(dh1s - dads/a)*dads/a**2 - &
            (d2h1s + dads*dads/a**2)/a + (h1 - lna)*d2ads2/a**2 + d2ads2/a**2
      sderv2(1,1) = sderv2(1,1) - d1*d2s
    endif
    call emfuncs(nemorder,u,0.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp,a,e1,dh1,d2h1,dh1s,d2h1s,d2h1m, &
                 d3h1,d3h1m,lgrad1,lgrad2,.false.)
    rstrd(1) = rstrd(1) - qtot*dh1s*angstoev
    assum = assum - dh1s
    if (lgrad2) then
      sderv2(1,1) = sderv2(1,1) - qtot*d2h1s*angstoev
    endif
  endif
  if (latomicstress) then
    assum = angstoev*assum
    do i = procid+1,numat,nprocs
      qi = qf(i)*occuf(i)
      atomicstress(1,i) = atomicstress(1,i) + assum*qi
    enddo
  endif
!
  t2 = g_cpu_time()
  tatom = tatom + t2 - t1
#ifdef TRACE
  call trace_out('real1Dq')
#endif
!
  return
  end
