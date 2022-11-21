  subroutine poccon(loptx,lopty,loptz,loptr)
!
!  This routine creates constraints that are needed to ensure
!  that species with partial occupancies on the same site
!  follow each other's motions.
!
!  lopt = array of logical flags according to whether a
!         parameter can be varied or not. Passed from setcfg.
!
!   6/95 Initially created.
!   5/03 Style updated
!   4/05 Intent added
!   2/06 Bug in referencing of arrays for X case fixed
!   7/13 Handling of breathing shells added
!   2/18 Trace added
!   3/19 Constraint arrays changed to have index and type
!   3/19 ltmp changed for individual arrays for x, y and z
!  10/19 Modified for rigid molecules
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
!  Julian Gale, CIC, Curtin University, October 2019
!
  use configurations
  use control
  use current
  use element,        only : maxele
  use optimisation
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  logical,     intent(inout) :: loptr(*)
  logical,     intent(inout) :: loptx(*)
  logical,     intent(inout) :: lopty(*)
  logical,     intent(inout) :: loptz(*)
!
!  Local variables
!
  integer(i4)                :: i
  integer(i4)                :: ii
  integer(i4)                :: j
  integer(i4)                :: jj
  integer(i4)                :: k
  integer(i4)                :: nati
  logical                    :: lcore
  real(dp)                   :: r
  real(dp)                   :: xal
  real(dp)                   :: yal
  real(dp)                   :: zal
  real(dp)                   :: xd
  real(dp)                   :: yd
  real(dp)                   :: zd
!
  if (index(keyword,'noco').ne.0) return
#ifdef TRACE
  call trace_in('poccon')
#endif
!*******************************************************************
!  Check for partial occupancy sites which need to be constrained  *
!*******************************************************************
  do ii = 2,nasymnomol
    i = nasymnomolptr(ii)
    nati = iatn(i)
    lcore = (nati.le.maxele)
    xal = xalat(i)
    yal = yalat(i)
    zal = zalat(i)
!
!  Loop over previous atoms to find any on the same site
!
    do jj = 1,ii-1
      j = nasymnomolptr(jj)
!
!  Check that both species are of the same type
!
      if (iatn(j).le.maxele.and.lcore.or.iatn(j).gt.maxele.and..not.lcore) then
!
!  Check distance 
!
        xd = xalat(j) - xal
        yd = yalat(j) - yal
        zd = zalat(j) - zal
        r = xd*xd + yd*yd + zd*zd
        if (r.lt.1.0d-10) then
          if (loptx(jj)) then
            if (ncontot.ge.maxcontot) then
              maxcontot = ncontot + 10
              call changemaxcontot
            endif
!
!  Add new constraint putting x coordinate of i equal to that of x
!
            if (ncf.lt.ncfg) then
              do k = ncontot,n1con(ncf+1),-1
                ncvarindcfg(k+1) = ncvarindcfg(k)
                ncvartypcfg(k+1) = ncvartypcfg(k)
                ncfixindcfg(k+1) = ncfixindcfg(k)
                ncfixtypcfg(k+1) = ncfixtypcfg(k)
                concocfg(k+1) = concocfg(k)
                nconcfg(k+1) = nconcfg(k)
                conaddcfg(k+1) = conaddcfg(k)
              enddo
              do k = ncf+1,ncfg
                n1con(k) = n1con(k) + 1
              enddo
            endif
            ncontot = ncontot + 1
            loptx(ii) = .false.
            ncon = ncon + 1
            ncvarindcfg(ncfst+ncon) = jj
            ncvartypcfg(ncfst+ncon) = iopt_xf
            ncfixindcfg(ncfst+ncon) = ii
            ncfixtypcfg(ncfst+ncon) = iopt_xf
            concocfg(ncfst+ncon) = 1.0_dp
            conaddcfg(ncfst+ncon) = 0.0_dp
            nconcfg(ncfst+ncon) = ncf
!
!  Check to see whether any other coordinates of i are dependent on the x coordinate of i
!
            do k = 1,ncon-1
              if (ncvarindcfg(ncfst+k).eq.ii.and.ncvartypcfg(ncfst+k).eq.iopt_xf) then
                ncvarindcfg(ncfst+k) = jj
              endif
            enddo
          endif
          if (lopty(jj)) then
            if (ncontot.ge.maxcontot) then
              maxcontot = ncontot + 10
              call changemaxcontot
            endif
!
!  Add new constraint putting y coordinate of i equal to that of y
!
            if (ncf.lt.ncfg) then
              do k = ncontot,n1con(ncf+1),-1
                ncvarindcfg(k+1) = ncvarindcfg(k)
                ncvartypcfg(k+1) = ncvartypcfg(k)
                ncfixindcfg(k+1) = ncfixindcfg(k)
                ncfixtypcfg(k+1) = ncfixtypcfg(k)
                concocfg(k+1) = concocfg(k)
                nconcfg(k+1) = nconcfg(k)
                conaddcfg(k+1) = conaddcfg(k)
              enddo
              do k = ncf+1,ncfg
                n1con(k) = n1con(k) + 1
              enddo
            endif
            ncontot = ncontot + 1
            lopty(ii) = .false.
            ncon = ncon + 1
            ncvarindcfg(ncfst+ncon) = jj
            ncvartypcfg(ncfst+ncon) = iopt_yf
            ncfixindcfg(ncfst+ncon) = ii
            ncfixtypcfg(ncfst+ncon) = iopt_yf
            concocfg(ncfst+ncon) = 1.0_dp
            conaddcfg(ncfst+ncon) = 0.0_dp
            nconcfg(ncfst+ncon) = ncf
!
!  Check to see whether any other coordinates of i are dependent on the y coordinate of i
!
            do k = 1,ncon-1
              if (ncvarindcfg(ncfst+k).eq.ii.and.ncvartypcfg(ncfst+k).eq.iopt_yf) then
                ncvarindcfg(ncfst+k) = jj
              endif
            enddo
          endif
          if (loptz(jj)) then
            if (ncontot.ge.maxcontot) then
              maxcontot = ncontot + 10
              call changemaxcontot
            endif
!
!  Add new constraint putting z coordinate of i equal to that of z
!
            if (ncf.lt.ncfg) then
              do k = ncontot,n1con(ncf+1),-1
                ncvarindcfg(k+1) = ncvarindcfg(k)
                ncvartypcfg(k+1) = ncvartypcfg(k)
                ncfixindcfg(k+1) = ncfixindcfg(k)
                ncfixtypcfg(k+1) = ncfixtypcfg(k)
                concocfg(k+1) = concocfg(k)
                nconcfg(k+1) = nconcfg(k)
                conaddcfg(k+1) = conaddcfg(k)
              enddo
              do k = ncf+1,ncfg
                n1con(k) = n1con(k) + 1
              enddo
            endif
            ncontot = ncontot + 1
            loptz(ii) = .false.
            ncon = ncon + 1
            ncvarindcfg(ncfst+ncon) = jj
            ncvartypcfg(ncfst+ncon) = iopt_zf
            ncfixindcfg(ncfst+ncon) = ii
            ncfixtypcfg(ncfst+ncon) = iopt_zf
            concocfg(ncfst+ncon) = 1.0_dp
            conaddcfg(ncfst+ncon) = 0.0_dp
            nconcfg(ncfst+ncon) = ncf
!
!  Check to see whether any other coordinates of i are dependent on the z coordinate of i
!
            do k = 1,ncon-1
              if (ncvarindcfg(ncfst+k).eq.ii.and.ncvartypcfg(ncfst+k).eq.iopt_zf) then
                ncvarindcfg(ncfst+k) = jj
              endif
            enddo
          endif
!
!  Breathing shell constraints
!
          if (lbsmat(nsft+i).and.lbsmat(nsft+j)) then
            if (loptr(jj)) then
              if (ncontot.ge.maxcontot) then
                maxcontot = ncontot + 10
                call changemaxcontot
              endif
!
!  Add new constraint putting radius of i equal to that of j
!
              if (ncf.lt.ncfg) then
                do k = ncontot,n1con(ncf+1),-1
                  ncvarindcfg(k+1) = ncvarindcfg(k)
                  ncvartypcfg(k+1) = ncvartypcfg(k)
                  ncfixindcfg(k+1) = ncfixindcfg(k)
                  ncfixtypcfg(k+1) = ncfixtypcfg(k)
                  concocfg(k+1) = concocfg(k)
                  nconcfg(k+1) = nconcfg(k)
                  conaddcfg(k+1) = conaddcfg(k)
                enddo
                do k = ncf+1,ncfg
                  n1con(k) = n1con(k) + 1
                enddo
              endif
              ncontot = ncontot + 1
              loptr(ii) = .false.
              ncon = ncon + 1
              ncvarindcfg(ncfst+ncon) = jj
              ncvartypcfg(ncfst+ncon) = iopt_radius
              ncfixindcfg(ncfst+ncon) = ii
              ncfixtypcfg(ncfst+ncon) = iopt_radius
              concocfg(ncfst+ncon) = 1.0_dp
              conaddcfg(ncfst+ncon) = 0.0_dp
              nconcfg(ncfst+ncon) = ncf
!
!  Check to see whether any other radii are dependent on the radius of i
!
              do k = 1,ncon-1
                if (ncvarindcfg(ncfst+k).eq.ii.and.ncvartypcfg(ncfst+k).eq.iopt_radius) then
                  ncvarindcfg(ncfst+k) = jj
                endif
              enddo
            endif
          endif
        endif
      endif
    enddo
  enddo
#ifdef TRACE
  call trace_out('poccon')
#endif
!
  return
  end
