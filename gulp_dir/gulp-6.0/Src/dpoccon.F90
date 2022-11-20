  subroutine dpoccon(iop)
!
!  This routine creates constraints that are needed to ensure
!  that species with partial occupancies on the same site
!  follow each other's motions.
!
!  iop = array of integer flags according to whether a parameter 
!        can be varied or not. Passed from setcfg.
!
!   6/95 Initially created.
!   2/18 Trace added
!   3/19 Change of defect constraints to have index and type
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
  use control
  use current
  use defects
  use element,        only : maxele
  use optimisation
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4)        :: iop(*)
!
!  Local variables
!
  integer(i4)        :: i
  integer(i4)        :: ii
  integer(i4)        :: indjx
  integer(i4)        :: indjy
  integer(i4)        :: indjz
  integer(i4)        :: indx
  integer(i4)        :: indy
  integer(i4)        :: indz
  integer(i4)        :: j
  integer(i4)        :: jj
  integer(i4)        :: k
  integer(i4)        :: nati
  logical            :: lcore
  real(dp)           :: r
  real(dp)           :: xal
  real(dp)           :: yal
  real(dp)           :: zal
  real(dp)           :: xd
  real(dp)           :: yd
  real(dp)           :: zd
!
  if (index(keyword,'noco').ne.0) return
#ifdef TRACE
  call trace_in('dpoccon')
#endif
!*******************************************************************
!  Check for partial occupancy sites which need to be constrained  *
!*******************************************************************
  do i = 2,ndasym
    ii = ndsptr(i)
    nati = natdefe(ii)
    lcore = (nati.le.maxele)
    xal = xdefe(ii)
    yal = ydefe(ii)
    zal = zdefe(ii)
    indx = 3*(i-1) + 1
    indy = indx + 1
    indz = indy + 1
!
!  Loop over previous atoms to find any on the same site
!
    do j = 1,i-1
      jj = ndsptr(j)
!
!  Check that both species are of the same type
!
      if (natdefe(jj).le.maxele.and.lcore.or.natdefe(jj).gt.maxele.and..not.lcore) then
!
!  Check distance 
!
        xd = xdefe(jj) - xal
        yd = ydefe(jj) - yal
        zd = zdefe(jj) - zal
        r = xd*xd + yd*yd + zd*zd
        if (r.lt.1.0d-10) then
          indjx = 3*(j-1) + 1
          indjy = indjx + 1
          indjz = indjy + 1
          if (iop(indjx).eq.1) then
!
!  Add new constraint putting x coordinate of i equal to that of x
!
            ndcon = ndcon + 1
            if (ndcon.ge.maxdcon) then
              maxdcon  =  ndcon + 10
              call changemaxdcon
            endif
            iop(indx) = 0
            ncdvarind(ndcon) = j
            ncdvartyp(ndcon) = idopt_dx
            ncdfixind(ndcon) = i
            ncdfixtyp(ndcon) = idopt_dx
            dconco(ndcon) = 1.0_dp
!
!  Check to see whether any other coordinates of i are
!  dependent on the x coordinate of i
!
            do k = 1,ndcon-1
              if (ncdvarind(k).eq.indx.and.ncdvartyp(k).eq.idopt_dx) then
                ncdvarind(k) = j
              endif
            enddo
          endif
          if (iop(indjy).eq.1) then
!
!  Add new constraint putting y coordinate of i equal to that of y
!
            ndcon = ndcon + 1
            if (ndcon.ge.maxdcon) then
              maxdcon  =  ndcon + 10
              call changemaxdcon
            endif
            iop(indy) = 0
            ncdvarind(ndcon) = j
            ncdvartyp(ndcon) = idopt_dy
            ncdfixind(ndcon) = i
            ncdfixtyp(ndcon) = idopt_dy
            dconco(ndcon) = 1.0_dp
!
!  Check to see whether any other coordinates of i are
!  dependent on the y coordinate of i
!
            do k = 1,ndcon-1
              if (ncdvarind(k).eq.indy.and.ncdvartyp(k).eq.idopt_dy) then
                ncdvarind(k) = j
              endif
            enddo
          endif
          if (iop(indjz).eq.1) then
!
!  Add new constraint putting z coordinate of i equal to that of z
!
            ndcon = ndcon + 1
            if (ndcon.ge.maxdcon) then
              maxdcon  =  ndcon + 10
              call changemaxdcon
            endif
            iop(indz) = 0
            ncdvarind(ndcon) = j
            ncdvartyp(ndcon) = idopt_dz
            ncdfixind(ndcon) = i
            ncdfixtyp(ndcon) = idopt_dz
            dconco(ndcon) = 1.0_dp
!
!  Check to see whether any other coordinates of i are
!  dependent on the z coordinate of i
!
            do k = 1,ndcon-1
              if (ncdvarind(k).eq.indz.and.ncdvartyp(k).eq.idopt_dz) then
                ncdvarind(k) = j
              endif
            enddo
          endif
        endif
      endif
    enddo
  enddo
#ifdef TRACE
  call trace_out('dpoccon')
#endif
!
  return
  end
