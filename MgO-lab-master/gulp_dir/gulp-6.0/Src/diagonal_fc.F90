  subroutine diagonal_fc(lused1)
!
!  Subroutine for calculating the on-diagonal blocks of the
!  dynamical matrix from the force constants arranged by 
!  cell
!
!  11/14 Created from dynamicn_fc
!   4/17 Order of loop & if changed for efficiency
!   8/17 Handling of 1-D case corrected by calling real1Dp
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
  use derivatives
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  logical,                         intent(in)    :: lused1       ! If true then this the first derivative algorithm
!
!  Local variables
!
  integer(i4)                                    :: i
  integer(i4)                                    :: ix
  integer(i4)                                    :: iy
  integer(i4)                                    :: iz
  integer(i4)                                    :: j
  integer(i4)                                    :: jx
  integer(i4)                                    :: jy
  integer(i4)                                    :: jz
  integer(i4)                                    :: maxlim
  integer(i4)                                    :: mint
  integer(i4)                                    :: n
#ifdef TRACE
  call trace_in('diagonal_fc')
#endif
!
!  Check dimensions of derv2
!
  mint = 3*numat
  maxlim = mint
  if (nbsmat.gt.0) maxlim = maxlim + numat
  if (maxlim.gt.maxd2u) then
    maxd2u = maxlim
    call changemaxd2
  endif
  if (maxlim.gt.maxd2) then
    maxd2 = maxlim
    call changemaxd2
  endif
!
!  Zero second derivatives
!
  derv2(1:maxlim,1:maxlim) = 0.0_dp
!
!  If this is the first derivative algorithm then the self-terms need to be added
!
  if (lused1) then
    call realselffc
  endif
!
!  Calculate reciprocal space contribution to the diagonal elements
!
  if (lewald.and.ndim.gt.0) then
    if (ndim.gt.1) then
      call kindex
    endif
    if (ndim.eq.3) then
      call recip3Dp(0.0_dp,0.0_dp,0.0_dp)
    elseif (ndim.eq.2) then
      call recip2Dp(0.0_dp,0.0_dp)
    elseif (ndim.eq.1) then
      call real1Dp(0.0_dp)
    endif
    ix = -2
    iy = -1
    iz =  0
    do i = 1,numat
      ix = ix + 3
      iy = iy + 3
      iz = iz + 3
      jx = -2
      jy = -1
      jz =  0
      do j = 1,numat
        jx = jx + 3
        jy = jy + 3
        jz = jz + 3
        if (j.ne.i) then
          derv2(ix,ix) = derv2(ix,ix) - derv2(jx,ix)
          derv2(iy,ix) = derv2(iy,ix) - derv2(jy,ix)
          derv2(iz,ix) = derv2(iz,ix) - derv2(jz,ix)
          derv2(ix,iy) = derv2(ix,iy) - derv2(jx,iy)
          derv2(iy,iy) = derv2(iy,iy) - derv2(jy,iy)
          derv2(iz,iy) = derv2(iz,iy) - derv2(jz,iy)
          derv2(ix,iz) = derv2(ix,iz) - derv2(jx,iz)
          derv2(iy,iz) = derv2(iy,iz) - derv2(jy,iz)
          derv2(iz,iz) = derv2(iz,iz) - derv2(jz,iz)
        endif
      enddo
    enddo
  endif
!
!  Set diagonal blocks from force constant contributions
!
  ix = -2
  iy = -1
  iz =  0
  do i = 1,numat
    ix = ix + 3
    iy = iy + 3
    iz = iz + 3
    jx = -2
    jy = -1
    jz =  0
    do j = 1,numat
      jx = jx + 3
      jy = jy + 3
      jz = jz + 3
      if (j.ne.i) then
        do n = 1,nd2cells
          derv2(ix,ix) = derv2(ix,ix) - d2cell(jx,ix,n)
          derv2(iy,ix) = derv2(iy,ix) - d2cell(jy,ix,n)
          derv2(iz,ix) = derv2(iz,ix) - d2cell(jz,ix,n)
          derv2(ix,iy) = derv2(ix,iy) - d2cell(jx,iy,n)
          derv2(iy,iy) = derv2(iy,iy) - d2cell(jy,iy,n)
          derv2(iz,iy) = derv2(iz,iy) - d2cell(jz,iy,n)
          derv2(ix,iz) = derv2(ix,iz) - d2cell(jx,iz,n)
          derv2(iy,iz) = derv2(iy,iz) - d2cell(jy,iz,n)
          derv2(iz,iz) = derv2(iz,iz) - d2cell(jz,iz,n)
        enddo
      endif
    enddo
  enddo
#ifdef TRACE
  call trace_out('diagonal_fc')
#endif
!
  return
  end
