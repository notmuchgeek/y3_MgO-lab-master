  subroutine sumderv2p(n,n2i,n2,ldbsm)
!
!  Completes second derivative matrix : diagonal blocks = sum of the off diagonal blocks.
!  Called by both energy and defener. Parallel version.
!
!  ldbsm = flag for whether breathing shells are present in defect calc
!
!  11/02 Modified to add on on-diagonal elements
!   5/03 Modified for region 3
!   7/07 Diagonal blocks for fixed sites added
!   9/13 Second argument added to handle parallel case
!   5/17 Modified so that outer loop doesn't require use of ldbsm
!        and that indj for ldbsm is get from n2
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
  use derivatives
  use parallel
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in) :: n         ! Number of local atoms for outer loop
  integer(i4), intent(in) :: n2i(*)    ! Array that points from local atom to global index
  integer(i4), intent(in) :: n2        ! Number of atoms for inner loop
  logical,     intent(in) :: ldbsm
!
!  Local variables
!
  integer(i4)             :: i
  integer(i4)             :: ii
  integer(i4)             :: indi
  integer(i4)             :: indig
  integer(i4)             :: indj
  integer(i4)             :: ix
  integer(i4)             :: iy
  integer(i4)             :: iz
  integer(i4)             :: ixg
  integer(i4)             :: iyg
  integer(i4)             :: izg
  integer(i4)             :: j
  integer(i4)             :: jx
  integer(i4)             :: jy
  integer(i4)             :: jz
#ifdef TRACE
  call trace_in('sumderv2p')
#endif
!*************************
!  Sum over coordinates  *
!*************************
  do ii = 1,n
    i = n2i(ii)
    indi = 3*(ii-1)
    ix = indi + 1
    iy = indi + 2
    iz = indi + 3
    indig = 3*(i-1)
    ixg = indig + 1
    iyg = indig + 2
    izg = indig + 3
    derv2(ixg,ix) = derv2d(ix)
    derv2(ixg,iy) = 0.0_dp
    derv2(ixg,iz) = 0.0_dp
    derv2(iyg,ix) = 0.0_dp
    derv2(iyg,iy) = derv2d(iy)
    derv2(iyg,iz) = 0.0_dp
    derv2(izg,ix) = 0.0_dp
    derv2(izg,iy) = 0.0_dp
    derv2(izg,iz) = derv2d(iz)
    do j = 1,n2
      if (i.ne.j) then
        if (ldbsm.and.j.eq.n2) then
          indj = 4*(n2-1)
        else
          indj = 3*(j-1)
        endif
        jx = indj + 1
        jy = indj + 2
        jz = indj + 3
        derv2(ixg,ix) = derv2(ixg,ix) - derv2(jx,ix)
        derv2(ixg,iy) = derv2(ixg,iy) - derv2(jy,ix)
        derv2(ixg,iz) = derv2(ixg,iz) - derv2(jz,ix)
        derv2(iyg,ix) = derv2(iyg,ix) - derv2(jx,iy)
        derv2(iyg,iy) = derv2(iyg,iy) - derv2(jy,iy)
        derv2(iyg,iz) = derv2(iyg,iz) - derv2(jz,iy)
        derv2(izg,ix) = derv2(izg,ix) - derv2(jx,iz)
        derv2(izg,iy) = derv2(izg,iy) - derv2(jy,iz)
        derv2(izg,iz) = derv2(izg,iz) - derv2(jz,iz)
      endif
    enddo
  enddo
#ifdef TRACE
  call trace_out('sumderv2p')
#endif
!
  return
  end
