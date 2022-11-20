  subroutine diagonal_fcd(lused1)
!
!  Subroutine for calculating the on-diagonal blocks of the
!  dynamical matrix from the force constants arranged by 
!  cell
!  Distributed memory second derivative parallel version.
!
!   4/17 Created from diagonal_fc
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
  use parallel
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
  integer(i4)                                    :: ixg
  integer(i4)                                    :: iyg
  integer(i4)                                    :: izg
  integer(i4)                                    :: j
  integer(i4)                                    :: jx
  integer(i4)                                    :: jy
  integer(i4)                                    :: jz
  integer(i4)                                    :: maxlim
  integer(i4)                                    :: maxlimloc
  integer(i4)                                    :: mint
  integer(i4)                                    :: mintloc
  integer(i4)                                    :: n
#ifdef TRACE
  call trace_in('diagonal_fcd')
#endif
!
!  Check dimensions of derv2
!
  mint = 3*numat
  mintloc = 3*natomsonnode
  maxlim = mint
  maxlimloc = mintloc
  if (nbsmat.gt.0) then
    maxlim = maxlim + numat
    maxlimloc = maxlimloc + natomsonnode
  endif
  if (maxlimloc.gt.maxd2u) then
    maxd2u = maxlimloc
    call changemaxd2
  endif
  if (maxlim.gt.maxd2) then
    maxd2 = maxlim
    call changemaxd2
  endif
!
!  Zero second derivatives
!
  derv2(1:maxlim,1:maxlimloc) = 0.0_dp
!
!  If this is the first derivative algorithm then the self-terms need to be added
!
  if (lused1) then
    call realselffcd
  endif
!
!  Calculate reciprocal space contribution to the diagonal elements
!
  if (lewald.and.ndim.gt.0) then
    if (ndim.gt.1) then
      call kindex
    endif
    if (ndim.eq.3) then
      call recip3Dpd(0.0_dp,0.0_dp,0.0_dp)
    elseif (ndim.eq.2) then
      call recip2Dpd(0.0_dp,0.0_dp)
    elseif (ndim.eq.1) then
      call real1Dpd(0.0_dp)
    endif
    ix = -2
    iy = -1
    iz =  0
    do i = 1,natomsonnode
      ix = ix + 3
      iy = iy + 3
      iz = iz + 3
      ixg = 3*(node2atom(i) - 1) + 1
      iyg = ixg + 1
      izg = ixg + 2
      jx = -2
      jy = -1
      jz =  0
      do j = 1,numat
        jx = jx + 3
        jy = jy + 3
        jz = jz + 3
        if (j.ne.node2atom(i)) then
          derv2(ixg,ix) = derv2(ixg,ix) - derv2(jx,ix)
          derv2(iyg,ix) = derv2(iyg,ix) - derv2(jy,ix)
          derv2(izg,ix) = derv2(izg,ix) - derv2(jz,ix)
          derv2(ixg,iy) = derv2(ixg,iy) - derv2(jx,iy)
          derv2(iyg,iy) = derv2(iyg,iy) - derv2(jy,iy)
          derv2(izg,iy) = derv2(izg,iy) - derv2(jz,iy)
          derv2(ixg,iz) = derv2(ixg,iz) - derv2(jx,iz)
          derv2(iyg,iz) = derv2(iyg,iz) - derv2(jy,iz)
          derv2(izg,iz) = derv2(izg,iz) - derv2(jz,iz)
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
  do i = 1,natomsonnode
    ix = ix + 3
    iy = iy + 3
    iz = iz + 3
    ixg = 3*(node2atom(i) - 1) + 1
    iyg = ixg + 1
    izg = ixg + 2
    jx = -2
    jy = -1
    jz =  0
    do j = 1,numat
      jx = jx + 3
      jy = jy + 3
      jz = jz + 3
      if (j.ne.node2atom(i)) then
        do n = 1,nd2cells
          derv2(ixg,ix) = derv2(ixg,ix) - d2cell(jx,ix,n)
          derv2(iyg,ix) = derv2(iyg,ix) - d2cell(jy,ix,n)
          derv2(izg,ix) = derv2(izg,ix) - d2cell(jz,ix,n)
          derv2(ixg,iy) = derv2(ixg,iy) - d2cell(jx,iy,n)
          derv2(iyg,iy) = derv2(iyg,iy) - d2cell(jy,iy,n)
          derv2(izg,iy) = derv2(izg,iy) - d2cell(jz,iy,n)
          derv2(ixg,iz) = derv2(ixg,iz) - d2cell(jx,iz,n)
          derv2(iyg,iz) = derv2(iyg,iz) - d2cell(jy,iz,n)
          derv2(izg,iz) = derv2(izg,iz) - d2cell(jz,iz,n)
        enddo
      endif
    enddo
  enddo
#ifdef TRACE
  call trace_out('diagonal_fcd')
#endif
!
  return
  end
