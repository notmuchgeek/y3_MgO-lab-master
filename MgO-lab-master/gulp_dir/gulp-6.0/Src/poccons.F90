  subroutine poccons(loptx,lopty,loptz,loptr)
!
!  This routine creates constraints that are needed to ensure
!  that species with partial occupancies on the same site
!  follow each other's motions.
!  Spatial decomposition version.
!
!  lopt = array of logical flags according to whether a
!         parameter can be varied or not. Passed from
!         setcfg.
!
!   6/95 Initially created.
!   5/03 Style updated
!   4/05 Intent added
!   2/06 Bug in referencing of arrays for X case fixed
!   8/11 Algorithm completely changed to match that in setoccptrs
!        to fix a bug for parallel case.
!  11/11 lfound removed from search as this variable was not set
!   7/13 Handling of breathing shells added
!   2/18 Trace added
!   3/19 Constraint arrays changed to have index and type
!   3/19 ltmp changed for individual arrays for x, y and z
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
  use configurations
  use control
  use current
  use element,        only : maxele
  use optimisation
  use spatial,        only : nspcellat,nspcell2atptr,nspcellat1ptr,nspcellatptr
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
  integer(i4)                :: ind
  integer(i4)                :: indjr
  integer(i4)                :: indjx
  integer(i4)                :: indjy
  integer(i4)                :: indjz
  integer(i4)                :: ixyz
  integer(i4)                :: j
  integer(i4)                :: k
  integer(i4)                :: nati
  integer(i4)                :: ni
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
  call trace_in('poccons')
#endif
!*******************************************************************
!  Check for partial occupancy sites which need to be constrained  *
!*******************************************************************
  do i = 1,numat
    nati = nat(i)
    lcore = (nati.le.maxele)
    xal = xclat(i)
    yal = yclat(i)
    zal = zclat(i)
!
!  Translate atom to spatial decomposition cell
!
    ixyz = nspcell2atptr(i)
    ind = nspcellat1ptr(ixyz)
!
!  Loop over previous atoms to find any on the same site
!
    ni = 0
    do while (ni.lt.nspcellat(ixyz))
      ni = ni + 1
      j = nspcellatptr(ind+ni)
      if (j.lt.i) then
!
!  Check that both species are of the same type
!
        if (iatn(j).le.maxele.and.lcore.or.iatn(j).gt.maxele.and..not.lcore) then
!
!  Check distance 
!
          xd = xclat(j) - xal
          yd = yclat(j) - yal
          zd = zclat(j) - zal
          r = xd*xd + yd*yd + zd*zd
          if (r.lt.1.0d-10) then
            indjx = 3*(j-1) + nstrains + 1
            indjy = indjx + 1
            indjz = indjy + 1
            indjr = 3*numat + nstrains + j
            if (loptx(j)) then
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
              loptx(i) = .false.
              ncon = ncon + 1
              ncvarindcfg(ncfst+ncon) = j
              ncvartypcfg(ncfst+ncon) = iopt_xf
              ncfixindcfg(ncfst+ncon) = i
              ncfixtypcfg(ncfst+ncon) = iopt_xf
              concocfg(ncfst+ncon) = 1.0_dp
              conaddcfg(ncfst+ncon) = 0.0_dp
              nconcfg(ncfst+ncon) = ncf
!
!  Check to see whether any other coordinates of i are dependent on the x coordinate of i
!
              do k = 1,ncon-1
                if (ncvarindcfg(ncfst+k).eq.i.and.ncvartypcfg(ncfst+k).eq.iopt_xf) then
                  ncvarindcfg(ncfst+k) = j
                endif
              enddo
            endif
            if (lopty(j)) then
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
              lopty(i) = .false.
              ncon = ncon + 1
              ncvarindcfg(ncfst+ncon) = j
              ncvartypcfg(ncfst+ncon) = iopt_yf
              ncfixindcfg(ncfst+ncon) = i
              ncfixtypcfg(ncfst+ncon) = iopt_yf
              concocfg(ncfst+ncon) = 1.0_dp
              conaddcfg(ncfst+ncon) = 0.0_dp
              nconcfg(ncfst+ncon) = ncf
!
!  Check to see whether any other coordinates of i are dependent on the y coordinate of i
!
              do k = 1,ncon-1
                if (ncvarindcfg(ncfst+k).eq.i.and.ncvartypcfg(ncfst+k).eq.iopt_yf) then
                  ncvarindcfg(ncfst+k) = j
                endif
              enddo
            endif
            if (loptz(j)) then
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
              loptz(i) = .false.
              ncon = ncon + 1
              ncvarindcfg(ncfst+ncon) = j
              ncvartypcfg(ncfst+ncon) = iopt_zf
              ncfixindcfg(ncfst+ncon) = i
              ncfixtypcfg(ncfst+ncon) = iopt_zf
              concocfg(ncfst+ncon) = 1.0_dp
              conaddcfg(ncfst+ncon) = 0.0_dp
              nconcfg(ncfst+ncon) = ncf
!
!  Check to see whether any other coordinates of i are dependent on the z coordinate of i
!
              do k = 1,ncon-1
                if (ncvarindcfg(ncfst+k).eq.i.and.ncvartypcfg(ncfst+k).eq.iopt_zf) then
                  ncvarindcfg(ncfst+k) = j
                endif
              enddo
            endif
!
!  Breathing shells
!
            if (lbsmat(nsft+i).and.lbsmat(nsft+j)) then
              if (loptr(j)) then
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
                loptr(i) = .false.
                ncon = ncon + 1
                ncvarindcfg(ncfst+ncon) = j
                ncvartypcfg(ncfst+ncon) = iopt_radius
                ncfixindcfg(ncfst+ncon) = i
                ncfixtypcfg(ncfst+ncon) = iopt_radius
                concocfg(ncfst+ncon) = 1.0_dp
                conaddcfg(ncfst+ncon) = 0.0_dp
                nconcfg(ncfst+ncon) = ncf
!
!  Check to see whether any other radii are dependent on the radius of i
!
                do k = 1,ncon-1
                  if (ncvarindcfg(ncfst+k).eq.i.and.ncvartypcfg(ncfst+k).eq.iopt_radius) then
                    ncvarindcfg(ncfst+k) = j
                  endif
                enddo
              endif
            endif
          endif
        endif
      endif
    enddo
  enddo
#ifdef TRACE
  call trace_out('poccons')
#endif
!
  return
  end
