  subroutine rpair(i,j,indri,ix,iy,iz,jx,jy,jz,lgrad2,xcrd,ycrd, &
    zcrd,deriv,deriv2,rderiv,rtrm1,rtrm2,radi,radj)
!
!  Subroutine performs derivative operations - called only from
!  real12a to save code duplication. Only called if lgrad1=.true.
!
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
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4)    :: i
  integer(i4)    :: indri
  integer(i4)    :: ix
  integer(i4)    :: iy
  integer(i4)    :: iz
  integer(i4)    :: j
  integer(i4)    :: jx
  integer(i4)    :: jy
  integer(i4)    :: jz
  logical        :: lgrad2
  real(dp)       :: deriv
  real(dp)       :: deriv2
  real(dp)       :: radi
  real(dp)       :: radj
  real(dp)       :: rderiv
  real(dp)       :: rtrm1
  real(dp)       :: rtrm2
  real(dp)       :: xcrd
  real(dp)       :: ycrd
  real(dp)       :: zcrd
!
!  Local variables
!
  real(dp)       :: rpd1
  real(dp)       :: rpd2
  real(dp)       :: rpd3
  real(dp)       :: rpd4
  real(dp)       :: rpd5
  real(dp)       :: rpd6
#ifdef TRACE
  call trace_in('rpair')
#endif
!***********************
!  Radial derivatives  *
!***********************
  if (radi.gt.0.0_dp) then
    raderv(i) = raderv(i) + rtrm1
    if (lgrad2) then
      derv2(indri,indri) = derv2(indri,indri) + rtrm2
    endif
  endif
!************************
!  Internal Derivatives *
!************************
!
!  First derivatives
!
  xdrv(i) = xdrv(i) - deriv*xcrd
  ydrv(i) = ydrv(i) - deriv*ycrd
  zdrv(i) = zdrv(i) - deriv*zcrd
!
!  Second derivatives
!
  if (lgrad2) then
    if (radi.gt.0.0_dp) then
      derv2(indri,ix) = derv2(indri,ix) - rderiv*xcrd
      derv2(indri,iy) = derv2(indri,iy) - rderiv*ycrd
      derv2(indri,iz) = derv2(indri,iz) - rderiv*zcrd
    endif
    rpd1 = xcrd*xcrd
    rpd2 = ycrd*ycrd
    rpd3 = zcrd*zcrd
    rpd4 = ycrd*zcrd
    rpd5 = xcrd*zcrd
    rpd6 = xcrd*ycrd
    derv2(ix,jx) = derv2(ix,jx) - deriv2*rpd1
    derv2(ix,jy) = derv2(ix,jy) - deriv2*rpd6
    derv2(ix,jz) = derv2(ix,jz) - deriv2*rpd5
    derv2(iy,jy) = derv2(iy,jy) - deriv2*rpd2
    derv2(iy,jz) = derv2(iy,jz) - deriv2*rpd4
    derv2(iz,jz) = derv2(iz,jz) - deriv2*rpd3
    derv2(ix,jx) = derv2(ix,jx) - deriv
    derv2(iy,jy) = derv2(iy,jy) - deriv
    derv2(iz,jz) = derv2(iz,jz) - deriv
  endif
#ifdef TRACE
  call trace_out('rpair')
#endif
!
  return
  end
!
  subroutine rpaird(i,j,indri,indrif,ix,iy,iz,jx,jy,jz,lgrad2, &
    xcrd,ycrd,zcrd,deriv,deriv2,rderiv,rtrm1,rtrm2,radi,radj)
!
!  Subroutine performs derivative operations - called only from
!  real12ad to save code duplication. Only called if lgrad1=.true.
!
!   5/17 Created from rpair for parallel case
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
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4)    :: i
  integer(i4)    :: indri
  integer(i4)    :: indrif
  integer(i4)    :: ix
  integer(i4)    :: iy
  integer(i4)    :: iz
  integer(i4)    :: j
  integer(i4)    :: jx
  integer(i4)    :: jy
  integer(i4)    :: jz
  logical        :: lgrad2
  real(dp)       :: deriv
  real(dp)       :: deriv2
  real(dp)       :: radi
  real(dp)       :: radj
  real(dp)       :: rderiv
  real(dp)       :: rtrm1
  real(dp)       :: rtrm2
  real(dp)       :: xcrd
  real(dp)       :: ycrd
  real(dp)       :: zcrd
!
!  Local variables
!
  real(dp)       :: rpd1
  real(dp)       :: rpd2
  real(dp)       :: rpd3
  real(dp)       :: rpd4
  real(dp)       :: rpd5
  real(dp)       :: rpd6
#ifdef TRACE
  call trace_in('rpaird')
#endif
!***********************
!  Radial derivatives  *
!***********************
  if (radi.gt.0.0_dp) then
    raderv(i) = raderv(i) + rtrm1
    if (lgrad2) then
      derv2(indrif,indri) = derv2(indrif,indri) + rtrm2
    endif
  endif
!************************
!  Internal Derivatives *
!************************
!
!  First derivatives
!
  xdrv(i) = xdrv(i) - deriv*xcrd
  ydrv(i) = ydrv(i) - deriv*ycrd
  zdrv(i) = zdrv(i) - deriv*zcrd
!
!  Second derivatives
!
  if (lgrad2) then
    if (radi.gt.0.0_dp) then
      derv2(indrif,ix) = derv2(indrif,ix) - rderiv*xcrd
      derv2(indrif,iy) = derv2(indrif,iy) - rderiv*ycrd
      derv2(indrif,iz) = derv2(indrif,iz) - rderiv*zcrd
    endif
    rpd1 = xcrd*xcrd
    rpd2 = ycrd*ycrd
    rpd3 = zcrd*zcrd
    rpd4 = ycrd*zcrd
    rpd5 = xcrd*zcrd
    rpd6 = xcrd*ycrd
    derv2(jx,ix) = derv2(jx,ix) - deriv2*rpd1
    derv2(jx,iy) = derv2(jx,iy) - deriv2*rpd6
    derv2(jx,iz) = derv2(jx,iz) - deriv2*rpd5
    derv2(jy,iy) = derv2(jy,iy) - deriv2*rpd2
    derv2(jy,iz) = derv2(jy,iz) - deriv2*rpd4
    derv2(jz,iz) = derv2(jz,iz) - deriv2*rpd3
    derv2(jx,ix) = derv2(jx,ix) - deriv
    derv2(jy,iy) = derv2(jy,iy) - deriv
    derv2(jz,iz) = derv2(jz,iz) - deriv
  endif
#ifdef TRACE
  call trace_out('rpaird')
#endif
!
  return
  end
