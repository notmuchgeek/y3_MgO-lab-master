  subroutine ramantrm(j,rmnx,rmny,rmnz)
!
!  Subroutine for calculating Raman intensity terms
!  Simple approximation based on upon sum of unit
!  vectors over bonded atoms. Only include cores in
!  the sum.
!
!   7/05 Intent added
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
  use current
  use element
  use molecule
  use species
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)  :: j
  real(dp),    intent(out) :: rmnx
  real(dp),    intent(out) :: rmny
  real(dp),    intent(out) :: rmnz
!
!  Local variables
!
  integer(i4)            :: ii
  integer(i4)            :: k
  integer(i4)            :: nj
  integer(i4)            :: nk
  real(dp)               :: cut2
  real(dp)               :: r
  real(dp)               :: rcut
  real(dp)               :: rj
  real(dp)               :: rk
  real(dp)               :: xcd
  real(dp)               :: ycd
  real(dp)               :: zcd
  real(dp)               :: xcdi
  real(dp)               :: ycdi
  real(dp)               :: zcdi
  real(dp)               :: xcrd
  real(dp)               :: ycrd
  real(dp)               :: zcrd
#ifdef TRACE
  call trace_in('ramantrm')
#endif
!
!  Initialise Raman terms
!
  rmnx = 0.0_dp
  rmny = 0.0_dp
  rmnz = 0.0_dp
!
!  Get attributes of j, the target atom - must be a core
!
  xcd = xclat(j)
  ycd = yclat(j)
  zcd = zclat(j)
  nj = nat(j)
  rj = rcov(nj)
!
!  If radius of j is zero then there can be no bonds so we are finished here
!
  if (rj.eq.0.0_dp) then
#ifdef TRACE
    call trace_out('ramantrm')
#endif
    return
  endif
!
!  Loop over second site
!
  do k = 1,numat
    nk = nat(k)
    if (nk.le.maxele) then
!
!  k is a core
!
      rk = rcov(nk)
      rcut = rtol*(rj + rk)
      cut2 = rcut*rcut
      xcdi = xclat(k) - xcd
      ycdi = yclat(k) - ycd
      zcdi = zclat(k) - zcd
!
!  Loop over unit cells
!
      do ii = 1,iimax
        xcrd = xcdi + xvec1cell(ii)
        ycrd = ycdi + yvec1cell(ii)
        zcrd = zcdi + zvec1cell(ii)
        r = xcrd*xcrd + ycrd*ycrd + zcrd*zcrd
        if (r.le.cut2.and.r.gt.0.0001_dp) then
          r = sqrt(r)
!
!  Add Raman contribution
!
          rmnx = rmnx + xcrd
          rmny = rmny + ycrd
          rmnz = rmnz + zcrd
        endif
!
!  End of loop over cell vectors
!
      enddo
    endif
!
!  End loop over k
!
  enddo
#ifdef TRACE
  call trace_out('ramantrm')
#endif
!
  return
  end
