  subroutine gamaxmin(xc,xmax,xmin,imode)
!
!  Create maximum and minimum limits for gafit or gaopt
!  imode = indicates whether call is from gafit or gaopt
!          1 => gafit
!          2 => gaopt
!
!  10/98 Codes for fitting variables simplified
!  11/06 Error in addressing of xmaxcfg and xmincfg fixed
!   1/09 Use of nfvar for two-body potentials modified
!   2/18 Trace added
!   3/19 iopt replaced by ioptindex and iopttype
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
  use current
  use fitting
  use genetic,       only : xmaxcfg, xmincfg
  use optimisation
  use parallel
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  use two
  implicit none
!
!  Passed variables
!
  integer(i4)        :: imode
  real(dp)           :: xc(*)
  real(dp)           :: xmax(*)
  real(dp)           :: xmin(*)
!
!  Local variables
!
  integer(i4)        :: i
  integer(i4)        :: ind
  integer(i4)        :: nf
  integer(i4)        :: np
  integer(i4)        :: nt
  integer(i4)        :: nv
#ifdef TRACE
  call trace_in('gamaxmin')
#endif
!
  if (nbsm.gt.0) then
    call outerror('breathing shells not allowed with GA fitting',0_i4)
    call stopnow('gamaxmin')
  endif
  if (imode.eq.1) then
    do i = 1,nfit
      nt = nftyp(i)
      nf = nfpot(i)
      nv = nfvar(i)
      if (nt.eq.2.and.nv.le.5) then
        np = nptype(nf)
        ind = (np-1)*4 + nf
        if (ind.eq.2.or.ind.eq.22.or.ind.eq.38) then
          if (xmin(i).eq.0.0_dp) then
            xmin(i) = 0.5_dp*xc(i)
          else
            xmin(i) = xmin(i)/scale(i)
          endif
        else
          if (xmin(i).ne.0.0_dp) then
            xmin(i) = xmin(i)/scale(i)
          endif
        endif
        if (xmax(i).eq.0.0_dp) then
          xmax(i) = 2.0_dp*xc(i)
        else
          xmax(i) = xmax(i)/scale(i)
        endif
      else
        if (xmax(i).eq.0.0_dp) then
          xmax(i) = 2.0_dp*xc(i)
        else
          xmax(i) = xmax(i)/scale(i)
        endif
        if (xmin(i).ne.0.0_dp) then
          xmin(i) = xmin(i)/scale(i)
        endif
      endif
    enddo
  else
    do i = 1,nvar
      if (iopttype(i).eq.iopt_xf) then
        if (xmax(i).eq.0.0_dp) xmax(i) = xmaxcfg(1,ncf)
        if (xmin(i).eq.0.0_dp) xmin(i) = xmincfg(1,ncf)
      elseif (iopttype(i).eq.iopt_yf) then
        if (xmax(i).eq.0.0_dp) xmax(i) = xmaxcfg(2,ncf)
        if (xmin(i).eq.0.0_dp) xmin(i) = xmincfg(2,ncf)
      elseif (iopttype(i).eq.iopt_zf) then
        if (xmax(i).eq.0.0_dp) xmax(i) = xmaxcfg(3,ncf)
        if (xmin(i).eq.0.0_dp) xmin(i) = xmincfg(3,ncf)
      else
        if (xmax(i).eq.0.0_dp) xmax(i) = 1.5_dp
      endif
    enddo
  endif
#ifdef TRACE
  call trace_out('gamaxmin')
#endif
!
  return
  end
