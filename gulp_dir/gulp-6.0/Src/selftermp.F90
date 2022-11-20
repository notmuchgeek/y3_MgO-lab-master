  subroutine selftermp(xkv,ykv,zkv,xcom,ycom,zcom,derive0self,factor,fct,ofct,ospfct,dfct,npotl,npots, &
                      c6tot,d2self,i,j,ix,jx,escale,c6scale,lewaldtype,qli,qlj)
!
!  Calculates self term contributions to second derivatives for phonons.
!
!   6/20 Created from selfterm
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
!  Julian Gale, CIC, Curtin University, June 2020
!
  use control
  use current
  use derivatives, only : derv2, dervi
  use general,     only : etaw
  use kspace
  use molecule,    only : lphasecom
  use numbers,     only : third
  use shells,      only : ncsptr
#ifdef TRACE
  use trace,       only : trace_in, trace_out
#endif
  use two
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)    :: i            ! Atom number for first atom i
  integer(i4), intent(in)    :: ix           ! Coordinate index for i / x
  integer(i4), intent(in)    :: j            ! Atom number for second atom j
  integer(i4), intent(in)    :: jx           ! Coordinate index for j / x
  integer(i4), intent(in)    :: npotl(*)     ! Array containing a list of valid potentials
  integer(i4), intent(in)    :: npots        ! Number of valid potentials
  logical,     intent(in)    :: lewaldtype   ! Flag to indicate whether Ewald sum is being performed
  real(dp),    intent(in)    :: c6scale
  real(dp),    intent(in)    :: c6tot
  real(dp),    intent(inout) :: d2self
  real(dp),    intent(inout) :: derive0self
  real(dp),    intent(in)    :: dfct
  real(dp),    intent(in)    :: escale
  real(dp),    intent(in)    :: factor       ! Charges of i and j times fct
  real(dp),    intent(in)    :: fct          ! ofct times 1/r -> eV conversion factor
  real(dp),    intent(in)    :: ofct         ! Product of occupancies of i and j
  real(dp),    intent(in)    :: ospfct       ! Occupancy of i
  real(dp),    intent(in)    :: qli          ! Charge of i
  real(dp),    intent(in)    :: qlj          ! Charge of j
  real(dp),    intent(in)    :: xcom
  real(dp),    intent(in)    :: ycom
  real(dp),    intent(in)    :: zcom
  real(dp),    intent(in)    :: xkv
  real(dp),    intent(in)    :: ykv
  real(dp),    intent(in)    :: zkv
!
!  Local variables
!
  integer(i4)                :: iy
  integer(i4)                :: iz
  integer(i4)                :: jy
  integer(i4)                :: jz
  integer(i4)                :: k
  integer(i4)                :: npot
  integer(i4)                :: npt
  logical                    :: lc6loc
  real(dp)                   :: apt
  real(dp)                   :: c6prod
  real(dp)                   :: c6self1
  real(dp)                   :: c6self1d2
  real(dp)                   :: cosk
  real(dp)                   :: sink
  real(dp)                   :: eta3
  real(dp)                   :: setrm
  real(dp)                   :: twoeta3
#ifdef TRACE
  call trace_in('selftermp')
#endif
!
  lc6loc = (lc6.and.ndim.eq.3)
!
!  Set up second derivative position pointers
!
  if (i.ne.j) then
    iy = ix + 1
    iz = ix + 2
    jy = jx + 1
    jz = jx + 2
!
!  Compute phase factor for second derivatives
!
    if (lrigid.and.lphasecom) then
      cosk = xkv*(-xcom) + ykv*(-ycom) + zkv*(-zcom)
      sink = sin(cosk)
      cosk = cos(cosk)
    else
      sink = 0.0_dp
      cosk = 1.0_dp
    endif
  endif
!
!  Core-shell spring constant at zero distant correct second derivative matrix
!
  if (i.ne.j.and.ncsptr(i).eq.j) then
    do k = 1,npots
      npot = npotl(k)
      npt = nptype(npot)
      if (npt.eq.5.or.npt.eq.8.or.npt.eq.33) then
        apt = twopot(1,npot)*ospfct*dfct
        derv2(jx,ix) = derv2(jx,ix) - apt*cosk
        derv2(jy,iy) = derv2(jy,iy) - apt*cosk
        derv2(jz,iz) = derv2(jz,iz) - apt*cosk
!
        dervi(jx,ix) = dervi(jx,ix) - apt*sink
        dervi(jy,iy) = dervi(jy,iy) - apt*sink
        dervi(jz,iz) = dervi(jz,iz) - apt*sink
      endif
    enddo
  endif
  if (lewaldtype) then
!
!  Remove self second derivative from Wolf sum 
!
    if (lwolf) then
      twoeta3 = 2.0_dp*etaw*etaw*third
      setrm = factor*tweatpi
      d2self = d2self - fct*tweatpi
      derive0self = derive0self - tweatpi*fct*escale
    else
!
!  Remove self energy/second derivative from Ewald sums (charge and C6) 
!
      twoeta3 = 2.0_dp*eta*third
      setrm = factor*tweatpi
      d2self = d2self - fct*tweatpi
      derive0self = derive0self - tweatpi*fct*escale
    endif
    if (lc6loc) then
      eta3 = eta*eta*eta
      c6self1 = eta3*third*c6scale
      if (lc6loc.and.lc6one) then
        c6prod = ofct*c6f(i)*c6f(j)
      else
        c6prod = ofct*c6tot
      endif
    endif
    if (i.ne.j) then
      setrm = setrm*twoeta3
      if (lc6loc) then
        c6self1d2 = 0.25_dp*eta3*eta
        setrm = setrm - c6prod*c6self1d2
      endif
      setrm = dfct*setrm
!
      derv2(jx,ix) = derv2(jx,ix) - setrm*cosk
      derv2(jy,iy) = derv2(jy,iy) - setrm*cosk
      derv2(jz,iz) = derv2(jz,iz) - setrm*cosk
!
      dervi(jx,ix) = dervi(jx,ix) - setrm*sink
      dervi(jy,iy) = dervi(jy,iy) - setrm*sink
      dervi(jz,iz) = dervi(jz,iz) - setrm*sink
    endif
  endif
#ifdef TRACE
  call trace_out('selftermp')
#endif
!
  return
  end
