  subroutine real1fc3(i,j,nati,natj,d3,xji,yji,zji)
!
!  Subroutine for third derivatives of 1-D electrostatic energy
!  for i-j pair. Called from thirdorderfc3. 
!
!   4/15 Created from real1d3
!   2/18 Trace added
!   3/20 Location of angstoev changed to current
!   3/20 Call to emfunc change due to creation of separate emfuncs for strains
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
!  Julian Gale, CIC, Curtin University, March 2020
!
  use current
  use element,       only : maxele
  use general,       only : nemorder, smallself
  use qmedata,       only : maxloop
  use shells,        only : cuts
  use times
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed arguments
!
  integer(i4), intent(in)    :: i
  integer(i4), intent(in)    :: j
  integer(i4), intent(in)    :: nati
  integer(i4), intent(in)    :: natj
  real(dp),    intent(inout) :: d3(3,3,3)
  real(dp),    intent(in)    :: xji
  real(dp),    intent(in)    :: yji
  real(dp),    intent(in)    :: zji
!
!  Local variables
!
  integer                    :: m
  logical                    :: lcspair
  real(dp)                   :: acell
  real(dp)                   :: g_cpu_time
  real(dp)                   :: cut2s
  real(dp)                   :: d0
  real(dp)                   :: d1
  real(dp)                   :: d2
  real(dp)                   :: dh1(3)
  real(dp)                   :: dh2(3)
  real(dp)                   :: d2h1(6)
  real(dp)                   :: d2h2(6)
  real(dp)                   :: d3h1(10)
  real(dp)                   :: d3h2(10)
  real(dp)                   :: d3l
  real(dp)                   :: e1
  real(dp)                   :: e2
  real(dp)                   :: h1
  real(dp)                   :: h2
  real(dp)                   :: oci     
  real(dp)                   :: ocj 
  real(dp)                   :: qi  
  real(dp)                   :: qj
  real(dp)                   :: qij
  real(dp)                   :: r
  real(dp)                   :: rcut
  real(dp)                   :: rr
  real(dp)                   :: t1, t2
  real(dp)                   :: u
  real(dp)                   :: x
  real(dp)                   :: y
  real(dp)                   :: z
#ifdef TRACE
  call trace_in('real1fc3')
#endif
!
  t1 = g_cpu_time()
!
!  Set up local variables
!
  cut2s = cuts*cuts
  oci = occuf(i)
  qi = qf(i)*oci
  ocj = occuf(j)
  qj = qf(j)*ocj
  qij = qi*qj
  lcspair = (abs(nati-natj).eq.maxele.or.(oci+ocj).lt.1.0001d0)
  if (lcspair) then
    rcut = cut2s
  else
    rcut = smallself
  endif
  y = yji
  z = zji
!
!  Loop over number of cells in sum
!
  do m = -maxloop(1),maxloop(1)
!
!  Direct sum component over neutral cells
!
    acell = dble(m)*a
    x = acell + xji
    r = x*x + y*y + z*z
    if (r.gt.rcut) then
      r = sqrt(r)
      rr = 1.0_dp/r
      d0 = qij*angstoev*rr
      d1 = d0*rr*rr
      d2 = 3.0_dp*d1*rr*rr
      d3l = - 5.0_dp*d2*rr*rr
!
!  Calculate components of third derivative matrix
!
      d3(1,1,1) = d3(1,1,1) + x*x*x*d3l + 3.0_dp*x*d2
      d3(2,1,1) = d3(2,1,1) + x*y*x*d3l + 2.0_dp*y*d2
      d3(3,1,1) = d3(3,1,1) + x*z*x*d3l + 2.0_dp*z*d2
      d3(2,2,1) = d3(2,2,1) + y*y*x*d3l + 2.0_dp*x*d2
      d3(3,2,1) = d3(3,2,1) + y*z*x*d3l
      d3(3,3,1) = d3(3,3,1) + z*z*x*d3l + 2.0_dp*x*d2
      d3(2,2,2) = d3(2,2,2) + y*y*y*d3l + 3.0_dp*y*d2
      d3(3,2,2) = d3(3,2,2) + y*z*y*d3l + 2.0_dp*z*d2
      d3(3,3,2) = d3(3,3,2) + z*z*y*d3l + 2.0_dp*y*d2
      d3(3,3,3) = d3(3,3,3) + z*z*z*d3l + 3.0_dp*z*d2
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
  if (maxloop(1).gt.0) then
    u = (dble(maxloop(1))+0.5_dp)*a
    x = xji
!
!  H term
!
    call hfunc(u,+x,1.0_dp,y,z,h1,dh1,d2h1,d3h1,.true.,.true.,.true.)
    call hfunc(u,-x,-1.0_dp,y,z,h2,dh2,d2h2,d3h2,.true.,.true.,.true.)
    d2 = qij*angstoev/a
!
    d3(1,1,1) = d3(1,1,1) - d2*(d3h1(1) + d3h2(1))
    d3(2,1,1) = d3(2,1,1) - d2*(d3h1(2) + d3h2(2))
    d3(3,1,1) = d3(3,1,1) - d2*(d3h1(3) + d3h2(3))
    d3(2,2,1) = d3(2,2,1) - d2*(d3h1(4) + d3h2(4))
    d3(3,2,1) = d3(3,2,1) - d2*(d3h1(5) + d3h2(5))
    d3(3,3,1) = d3(3,3,1) - d2*(d3h1(6) + d3h2(6))
    d3(2,2,2) = d3(2,2,2) - d2*(d3h1(7) + d3h2(7))
    d3(3,2,2) = d3(3,2,2) - d2*(d3h1(8) + d3h2(8))
    d3(3,3,2) = d3(3,3,2) - d2*(d3h1(9) + d3h2(9))
    d3(3,3,3) = d3(3,3,3) - d2*(d3h1(10)+ d3h2(10))
!
!  E-M term
!
    call emfunc(nemorder,u,+x,1.0_dp,y,z,a,e1,dh1,d2h1,d3h1,.true.,.true.,.true.)
    call emfunc(nemorder,u,-x,-1.0_dp,y,z,a,e2,dh2,d2h2,d3h2,.true.,.true.,.true.)
    d2 = qij*angstoev
    d3(1,1,1) = d3(1,1,1) + d2*(d3h1(1) + d3h2(1))
    d3(2,1,1) = d3(2,1,1) + d2*(d3h1(2) + d3h2(2))
    d3(3,1,1) = d3(3,1,1) + d2*(d3h1(3) + d3h2(3))
    d3(2,2,1) = d3(2,2,1) + d2*(d3h1(4) + d3h2(4))
    d3(3,2,1) = d3(3,2,1) + d2*(d3h1(5) + d3h2(5))
    d3(3,3,1) = d3(3,3,1) + d2*(d3h1(6) + d3h2(6))
    d3(2,2,2) = d3(2,2,2) + d2*(d3h1(7) + d3h2(7))
    d3(3,2,2) = d3(3,2,2) + d2*(d3h1(8) + d3h2(8))
    d3(3,3,2) = d3(3,3,2) + d2*(d3h1(9) + d3h2(9))
    d3(3,3,3) = d3(3,3,3) + d2*(d3h1(10)+ d3h2(10))
  endif
!
  t2 = g_cpu_time()
  tatom = tatom + t2 - t1
#ifdef TRACE
  call trace_out('real1fc3')
#endif
!
  return
  end
