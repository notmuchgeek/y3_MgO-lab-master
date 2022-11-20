  subroutine reciptrmdp(xd,yd,zd,lgrad3,lstr,fct,d2,d2s,d3,d3s,d3ss)
!
!  Calculates the reciprocal space contribution to the derivatives of the polarisation energy
!
!   6/00 Created from recipktrmd3
!   2/01 Modifications for 2-D added
!   3/01 Derivatives for 2-D completed
!   3/03 Modified to allow for symmetry
!  11/07 Unused variables cleaned up
!   3/14 derfc changed to g_derfc for benefit of ChemShell
!   2/18 Trace added
!   4/19 lnorecip added
!   4/19 Third derivatives added as an option
!   5/19 Sign of d2s changed
!   5/19 Extra terms for d2s corrected
!   5/19 Final corrections made to 3D terms
!   5/19 2D second derivatives added for z term in reciprocal space
!   6/19 Call to calcf2d introduced for reciprocal space function and derivatives
!   7/19 Finite strain modifications added
!   8/19 Correction to 2D z-only term with finite strain 
!   8/19 Correction to d2s changes for non-finite strain case
!   3/20 Location of angstoev changed to current
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
!  Julian Gale, CIC, Curtin University, March 2020
!
  use g_constants,   only : pi
  use control,       only : lnorecip
  use current
  use derivatives,   only : lfinitestrain
  use kspace
  use m_strain,      only : gstrterms, strainddetds, straindet, gxyzterms
  use m_strain,      only : straind2detds2
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  logical,  intent(in)     :: lgrad3          ! If true then compute third derivatives
  logical,  intent(in)     :: lstr            ! If true then compute strain derivatives
  real(dp), intent(inout)  :: d2(3,3)         ! Cartesian second derivative matrix
  real(dp), intent(inout)  :: d2s(3,6)        ! Mixed Cartesian-strain second derivative matrix
  real(dp), intent(inout)  :: d3(3,3,3)       ! Cartesian third derivative matrix (if lgrad3 is true)
  real(dp), intent(inout)  :: d3s(3,6,3)      ! Cartesian-strain third derivative matrix (if lgrad3 is true)
  real(dp), intent(inout)  :: d3ss(3,6,6)     ! Cartesian-strain third derivative matrix (if lgrad3 is true)
  real(dp), intent(in)     :: fct
  real(dp), intent(in)     :: xd
  real(dp), intent(in)     :: yd
  real(dp), intent(in)     :: zd
!
!  Local variables
!
  integer(i4)              :: iv
  integer(i4)              :: kk
  integer(i4)              :: kl
  real(dp)                 :: arg
  real(dp)                 :: cosa
  real(dp)                 :: cosq
  real(dp)                 :: coss
  real(dp)                 :: g_derf
  real(dp)                 :: derfez
  real(dp)                 :: dexpz
  real(dp)                 :: dtrm1
  real(dp)                 :: dtrm2
  real(dp)                 :: etaz
  real(dp)                 :: etaz2
  real(dp)                 :: f
  real(dp)                 :: f1z
  real(dp)                 :: f1g
  real(dp)                 :: f2zz
  real(dp)                 :: f2zg
  real(dp)                 :: f2gg
  real(dp)                 :: f3zzz
  real(dp)                 :: f3zzg
  real(dp)                 :: f3zgg
  real(dp)                 :: kvec
  real(dp)                 :: rkvec
  real(dp)                 :: sina
  real(dp)                 :: sinq
  real(dp)                 :: sins
  real(dp)                 :: sins2
  real(dp)                 :: tmp1
  real(dp)                 :: tmp2
  real(dp)                 :: tmp3
  real(dp)                 :: tmp4
  real(dp)                 :: tmp5
  real(dp)                 :: tmp6
  real(dp)                 :: xrkk
  real(dp)                 :: yrkk
  real(dp)                 :: zrkk
#ifdef TRACE
  call trace_in('reciptrmdp')
#endif
!
!  If keyword is forcing reciprocal space to be skipped then goto exit point
!
  if (lnorecip) goto 999
!
!  Define constants
!
  rpieta = 1.0_dp/sqrt(pi*eta)
  rhseta = 0.5_dp/seta
!
  if (ndim.eq.2) then
!
!  First term - K vector independent
!
    etaz = seta*zd
    etaz2 = etaz*etaz
    dexpz  = exp(-etaz2)
    derfez = g_derf(etaz)
    dtrm2 = vol4pi*tweatpi*dexpz*angstoev*fct
    d2(3,3) = d2(3,3) - dtrm2
    if (lgrad3) then
      d3(3,3,3) = d3(3,3,3) + 2.0_dp*dtrm2*eta*zd
    endif
!
    if (lstr) then
      dtrm1 = vol4pi*derfez*angstoev*fct
      if (lfinitestrain) then
        do kk = 1,nstrains
          d2s(3,kk) = d2s(3,kk) - dtrm1*strainddetds(kk)*straindet
        enddo
      else
        d2s(3,1) = d2s(3,1) - dtrm1
        d2s(3,2) = d2s(3,2) - dtrm1
      endif
      if (lgrad3) then
        if (lfinitestrain) then
          do kk = 1,nstrains
            d3s(3,kk,3) = d3s(3,kk,3) - dtrm2*strainddetds(kk)*straindet
          enddo
        else
          d3s(3,1,3) = d3s(3,1,3) - dtrm2
          d3s(3,2,3) = d3s(3,2,3) - dtrm2
        endif
!
        if (lfinitestrain) then
          do kk = 1,nstrains
            do kl = 1,nstrains
              d3ss(3,kl,kk) = d3ss(3,kl,kk) - dtrm1*straindet*(2.0_dp*strainddetds(kl)*strainddetds(kk) &
                                              *straindet - straind2detds2(kl,kk))
            enddo
          enddo
        else
          d3ss(3,1,1) = d3ss(3,1,1) - dtrm1
          d3ss(3,2,1) = d3ss(3,2,1) - dtrm1
          d3ss(3,1,2) = d3ss(3,1,2) - dtrm1
          d3ss(3,2,2) = d3ss(3,2,2) - dtrm1
        endif
      endif
    endif
  endif
!
!  Set up strain derivatives
!
  if (lstr) then
    call gstrterms(ndim,maxkvec,nkvec,xrk,yrk,zrk,dg2ds,d2g2dx2,d2g2ds2,lgrad3)
    call gxyzterms(ndim,maxkvec,nkvec,xrk,yrk,zrk,dgds,d2gds2,lgrad3)
  endif
!**************************
!  Calculate derivatives  *
!**************************
  do iv = 1,nkvec
    xrkk = xrk(iv)
    yrkk = yrk(iv)
    zrkk = zrk(iv)
    if (ndim.eq.3) then
      arg = xrkk*xd + yrkk*yd + zrkk*zd
      cosa = cos(arg)
      sina = sin(arg)
      cosq = cosa*ktrm(iv)
      sinq = sina*ktrm(iv)
      coss = cosa*ktrms(iv)
      sins = sina*ktrms(iv)
      if (lgrad3) then
        sins2 = sina*ktrms2(iv)
      endif
    elseif (ndim.eq.2) then
      kvec = kmod(iv)
      rkvec = 1.0_dp/kvec
      arg = xrkk*xd + yrkk*yd
      cosa = cos(arg)
      sina = sin(arg)
      cosq = cosa*ktrm(iv)
      sinq = sina*ktrm(iv)
!
      call calcf2D(kvec,zd,seta,f,f1z,f1g,f2zz,f2zg,f2gg,f3zzz,f3zzg,f3zgg, &
                   .true.,.true.,lgrad3)
    endif
!
    tmp1 = xrkk*xrkk
    tmp2 = yrkk*yrkk
    tmp3 = zrkk*zrkk
    tmp4 = yrkk*zrkk
    tmp5 = xrkk*zrkk
    tmp6 = xrkk*yrkk
!
!  Calculate Cartesian second derivatives
!
    if (ndim.eq.3) then
      d2(1,1) = d2(1,1) - cosq*tmp1*fct
      d2(2,1) = d2(2,1) - cosq*tmp6*fct
      d2(2,2) = d2(2,2) - cosq*tmp2*fct
      d2(3,1) = d2(3,1) - cosq*tmp5*fct
      d2(3,2) = d2(3,2) - cosq*tmp4*fct
      d2(3,3) = d2(3,3) - cosq*tmp3*fct
    elseif (ndim.eq.2) then
      d2(1,1) = d2(1,1) - cosq*tmp1*f*fct
      d2(2,1) = d2(2,1) - cosq*tmp6*f*fct
      d2(2,2) = d2(2,2) - cosq*tmp2*f*fct
      d2(3,1) = d2(3,1) - sinq*xrkk*f1z*fct
      d2(3,2) = d2(3,2) - sinq*yrkk*f1z*fct
      d2(3,3) = d2(3,3) + cosq*f2zz*fct
    endif
    if (lgrad3) then
!
!  Calculate Cartesian third derivatives
!
      if (ndim.eq.3) then
        d3(1,1,1) = d3(1,1,1) + sinq*tmp1*xrkk*fct
        d3(2,1,1) = d3(2,1,1) + sinq*tmp6*xrkk*fct
        d3(2,2,1) = d3(2,2,1) + sinq*tmp2*xrkk*fct
        d3(2,2,2) = d3(2,2,2) + sinq*tmp2*yrkk*fct
        d3(3,1,1) = d3(3,1,1) + sinq*tmp5*xrkk*fct
        d3(3,2,1) = d3(3,2,1) + sinq*tmp4*xrkk*fct
        d3(3,3,1) = d3(3,3,1) + sinq*tmp3*xrkk*fct
        d3(3,2,2) = d3(3,2,2) + sinq*tmp4*yrkk*fct
        d3(3,3,2) = d3(3,3,2) + sinq*tmp3*yrkk*fct
        d3(3,3,3) = d3(3,3,3) + sinq*tmp3*zrkk*fct
      elseif (ndim.eq.2) then
! 0 z component
        d3(1,1,1) = d3(1,1,1) + sinq*tmp1*xrkk*f*fct
        d3(2,1,1) = d3(2,1,1) + sinq*tmp6*xrkk*f*fct
        d3(2,2,1) = d3(2,2,1) + sinq*tmp2*xrkk*f*fct
        d3(2,2,2) = d3(2,2,2) + sinq*tmp2*yrkk*f*fct
! 1 z component
        d3(3,1,1) = d3(3,1,1) - cosq*tmp1*f1z*fct
        d3(3,2,1) = d3(3,2,1) - cosq*tmp6*f1z*fct
        d3(3,2,2) = d3(3,2,2) - cosq*tmp2*f1z*fct
! 2 z components
        d3(3,3,1) = d3(3,3,1) - sinq*xrkk*f2zz*fct
        d3(3,3,2) = d3(3,3,2) - sinq*yrkk*f2zz*fct
! 3 z components
        d3(3,3,3) = d3(3,3,3) + cosq*f3zzz*fct
      endif
    endif
    if (lstr) then
      if (ndim.eq.3) then
        do kk = 1,nstrains
          d2s(1,kk) = d2s(1,kk) + xrkk*dg2ds(iv,kk)*sins*fct
          d2s(2,kk) = d2s(2,kk) + yrkk*dg2ds(iv,kk)*sins*fct
          d2s(3,kk) = d2s(3,kk) + zrkk*dg2ds(iv,kk)*sins*fct
        enddo
      elseif (ndim.eq.2) then
        do kk = 1,nstrains
          d2s(1,kk) = d2s(1,kk) - xrkk*sinq*dg2ds(iv,kk)*fct*(f*rkvec**2 - f1g)
          d2s(2,kk) = d2s(2,kk) - yrkk*sinq*dg2ds(iv,kk)*fct*(f*rkvec**2 - f1g)
          d2s(3,kk) = d2s(3,kk) - cosq*dg2ds(iv,kk)*fct*(f2zg - f1z*rkvec**2)
        enddo
      endif
      if (lfinitestrain) then
!
!  Subtract volume/area derivative terms
!
        if (ndim.eq.3) then
          do kk = 1,nstrains
            d2s(1,kk) = d2s(1,kk) - xrkk*sinq*fct*strainddetds(kk)*straindet
            d2s(2,kk) = d2s(2,kk) - yrkk*sinq*fct*strainddetds(kk)*straindet
            d2s(3,kk) = d2s(3,kk) - zrkk*sinq*fct*strainddetds(kk)*straindet
          enddo
        elseif (ndim.eq.2) then
          do kk = 1,nstrains
            d2s(1,kk) = d2s(1,kk) - xrkk*sinq*f*fct*strainddetds(kk)*straindet
            d2s(2,kk) = d2s(2,kk) - yrkk*sinq*f*fct*strainddetds(kk)*straindet
            d2s(3,kk) = d2s(3,kk) + cosq*f1z*fct*strainddetds(kk)*straindet
          enddo
        endif
      else
        if (ndim.eq.3) then
          do kk = 1,3
            d2s(1,kk) = d2s(1,kk) - xrkk*sinq*fct
            d2s(2,kk) = d2s(2,kk) - yrkk*sinq*fct
            d2s(3,kk) = d2s(3,kk) - zrkk*sinq*fct
          enddo
        elseif (ndim.eq.2) then
          do kk = 1,2
            d2s(1,kk) = d2s(1,kk) - xrkk*sinq*f*fct
            d2s(2,kk) = d2s(2,kk) - yrkk*sinq*f*fct
            d2s(3,kk) = d2s(3,kk) + cosq*f1z*fct
          enddo
        endif
      endif
!
!  Derivatives of reciprocal lattice component with respect to strain
!
      if (ndim.eq.3) then
        d2s(1,1:nstrains) = d2s(1,1:nstrains) + dgds(iv,1,1:nstrains)*sinq*fct
        d2s(2,1:nstrains) = d2s(2,1:nstrains) + dgds(iv,2,1:nstrains)*sinq*fct
        d2s(3,1:nstrains) = d2s(3,1:nstrains) + dgds(iv,3,1:nstrains)*sinq*fct
      elseif (ndim.eq.2) then
        d2s(1,1:nstrains) = d2s(1,1:nstrains) + dgds(iv,1,1:nstrains)*sinq*f*fct
        d2s(2,1:nstrains) = d2s(2,1:nstrains) + dgds(iv,2,1:nstrains)*sinq*f*fct
      endif
      if (lgrad3) then
!
!  Calculate Cartesian/strain/Cartesian third derivatives
!
        if (ndim.eq.3) then
          do kk = 1,nstrains
            d3s(1,kk,1) = d3s(1,kk,1) + tmp1*dg2ds(iv,kk)*coss*fct
            d3s(2,kk,1) = d3s(2,kk,1) + tmp6*dg2ds(iv,kk)*coss*fct
            d3s(3,kk,1) = d3s(3,kk,1) + tmp5*dg2ds(iv,kk)*coss*fct
            d3s(1,kk,2) = d3s(1,kk,2) + tmp6*dg2ds(iv,kk)*coss*fct
            d3s(2,kk,2) = d3s(2,kk,2) + tmp2*dg2ds(iv,kk)*coss*fct
            d3s(3,kk,2) = d3s(3,kk,2) + tmp4*dg2ds(iv,kk)*coss*fct
            d3s(1,kk,3) = d3s(1,kk,3) + tmp5*dg2ds(iv,kk)*coss*fct
            d3s(2,kk,3) = d3s(2,kk,3) + tmp4*dg2ds(iv,kk)*coss*fct
            d3s(3,kk,3) = d3s(3,kk,3) + tmp3*dg2ds(iv,kk)*coss*fct
          enddo
!
          if (lfinitestrain) then
!
!  Subtract volume derivative terms
!
            do kk = 1,nstrains
              d3s(1,kk,1) = d3s(1,kk,1) - tmp1*cosq*fct*strainddetds(kk)*straindet
              d3s(2,kk,1) = d3s(2,kk,1) - tmp6*cosq*fct*strainddetds(kk)*straindet
              d3s(3,kk,1) = d3s(3,kk,1) - tmp5*cosq*fct*strainddetds(kk)*straindet
              d3s(1,kk,2) = d3s(1,kk,2) - tmp6*cosq*fct*strainddetds(kk)*straindet
              d3s(2,kk,2) = d3s(2,kk,2) - tmp2*cosq*fct*strainddetds(kk)*straindet
              d3s(3,kk,2) = d3s(3,kk,2) - tmp4*cosq*fct*strainddetds(kk)*straindet
              d3s(1,kk,3) = d3s(1,kk,3) - tmp5*cosq*fct*strainddetds(kk)*straindet
              d3s(2,kk,3) = d3s(2,kk,3) - tmp4*cosq*fct*strainddetds(kk)*straindet
              d3s(3,kk,3) = d3s(3,kk,3) - tmp3*cosq*fct*strainddetds(kk)*straindet
            enddo
          else
            do kk = 1,3
              d3s(1,kk,1) = d3s(1,kk,1) - tmp1*cosq*fct
              d3s(2,kk,1) = d3s(2,kk,1) - tmp6*cosq*fct
              d3s(3,kk,1) = d3s(3,kk,1) - tmp5*cosq*fct
              d3s(1,kk,2) = d3s(1,kk,2) - tmp6*cosq*fct
              d3s(2,kk,2) = d3s(2,kk,2) - tmp2*cosq*fct
              d3s(3,kk,2) = d3s(3,kk,2) - tmp4*cosq*fct
              d3s(1,kk,3) = d3s(1,kk,3) - tmp5*cosq*fct
              d3s(2,kk,3) = d3s(2,kk,3) - tmp4*cosq*fct
              d3s(3,kk,3) = d3s(3,kk,3) - tmp3*cosq*fct
            enddo
          endif
!
          d3s(1,1:nstrains,1) = d3s(1,1:nstrains,1) + xrkk*dgds(iv,1,1:nstrains)*cosq*fct
          d3s(2,1:nstrains,1) = d3s(2,1:nstrains,1) + xrkk*dgds(iv,2,1:nstrains)*cosq*fct
          d3s(3,1:nstrains,1) = d3s(3,1:nstrains,1) + xrkk*dgds(iv,3,1:nstrains)*cosq*fct
          d3s(1,1:nstrains,2) = d3s(1,1:nstrains,2) + yrkk*dgds(iv,1,1:nstrains)*cosq*fct
          d3s(2,1:nstrains,2) = d3s(2,1:nstrains,2) + yrkk*dgds(iv,2,1:nstrains)*cosq*fct
          d3s(3,1:nstrains,2) = d3s(3,1:nstrains,2) + yrkk*dgds(iv,3,1:nstrains)*cosq*fct
          d3s(1,1:nstrains,3) = d3s(1,1:nstrains,3) + zrkk*dgds(iv,1,1:nstrains)*cosq*fct
          d3s(2,1:nstrains,3) = d3s(2,1:nstrains,3) + zrkk*dgds(iv,2,1:nstrains)*cosq*fct
          d3s(3,1:nstrains,3) = d3s(3,1:nstrains,3) + zrkk*dgds(iv,3,1:nstrains)*cosq*fct
        elseif (ndim.eq.2) then
!
!  Derivatives for mod G
!
          do kk = 1,nstrains
            d3s(1,kk,1) = d3s(1,kk,1) + tmp1*cosq*dg2ds(iv,kk)*fct*(f1g - f*rkvec**2)
            d3s(2,kk,1) = d3s(2,kk,1) + tmp6*cosq*dg2ds(iv,kk)*fct*(f1g - f*rkvec**2)
            d3s(3,kk,1) = d3s(3,kk,1) + xrkk*sinq*dg2ds(iv,kk)*fct*(f2zg - f1z*rkvec**2)
            d3s(1,kk,2) = d3s(1,kk,2) + tmp6*cosq*dg2ds(iv,kk)*fct*(f1g - f*rkvec**2)
            d3s(2,kk,2) = d3s(2,kk,2) + tmp2*cosq*dg2ds(iv,kk)*fct*(f1g - f*rkvec**2)
            d3s(3,kk,2) = d3s(3,kk,2) + yrkk*sinq*dg2ds(iv,kk)*fct*(f2zg - f1z*rkvec**2)
            d3s(1,kk,3) = d3s(1,kk,3) + xrkk*sinq*dg2ds(iv,kk)*fct*(f2zg - f1z*rkvec**2)
            d3s(2,kk,3) = d3s(2,kk,3) + yrkk*sinq*dg2ds(iv,kk)*fct*(f2zg - f1z*rkvec**2)
            d3s(3,kk,3) = d3s(3,kk,3) - cosq*dg2ds(iv,kk)*fct*(f3zzg - f2zz*rkvec**2)
          enddo
!
!  Derivatives from components of G
!
          do kk = 1,nstrains
            d3s(1,kk,1) = d3s(1,kk,1) + cosq*f*fct*xrkk*dgds(iv,1,kk)
            d3s(2,kk,1) = d3s(2,kk,1) + cosq*f*fct*xrkk*dgds(iv,2,kk)
            d3s(1,kk,2) = d3s(1,kk,2) + cosq*f*fct*yrkk*dgds(iv,1,kk)
            d3s(2,kk,2) = d3s(2,kk,2) + cosq*f*fct*yrkk*dgds(iv,2,kk)
            d3s(1,kk,3) = d3s(1,kk,3) + sinq*f1z*fct*dgds(iv,1,kk)
            d3s(2,kk,3) = d3s(2,kk,3) + sinq*f1z*fct*dgds(iv,2,kk)
          enddo
!
!  Area derivatives
!
          if (lfinitestrain) then
            do kk = 1,nstrains
              d3s(1,kk,1) = d3s(1,kk,1) - tmp1*cosq*f*fct*strainddetds(kk)*straindet
              d3s(2,kk,1) = d3s(2,kk,1) - tmp6*cosq*f*fct*strainddetds(kk)*straindet
              d3s(3,kk,1) = d3s(3,kk,1) - xrkk*sinq*f1z*fct*strainddetds(kk)*straindet
              d3s(1,kk,2) = d3s(1,kk,2) - tmp6*cosq*f*fct*strainddetds(kk)*straindet
              d3s(2,kk,2) = d3s(2,kk,2) - tmp2*cosq*f*fct*strainddetds(kk)*straindet
              d3s(3,kk,2) = d3s(3,kk,2) - yrkk*sinq*f1z*fct*strainddetds(kk)*straindet
              d3s(1,kk,3) = d3s(1,kk,3) - xrkk*sinq*f1z*fct*strainddetds(kk)*straindet
              d3s(2,kk,3) = d3s(2,kk,3) - yrkk*sinq*f1z*fct*strainddetds(kk)*straindet
              d3s(3,kk,3) = d3s(3,kk,3) + cosq*f2zz*fct*strainddetds(kk)*straindet
            enddo
          else
            do kk = 1,2
              d3s(1,kk,1) = d3s(1,kk,1) - tmp1*cosq*f*fct
              d3s(2,kk,1) = d3s(2,kk,1) - tmp6*cosq*f*fct
              d3s(3,kk,1) = d3s(3,kk,1) - xrkk*sinq*f1z*fct
              d3s(1,kk,2) = d3s(1,kk,2) - tmp6*cosq*f*fct
              d3s(2,kk,2) = d3s(2,kk,2) - tmp2*cosq*f*fct
              d3s(3,kk,2) = d3s(3,kk,2) - yrkk*sinq*f1z*fct
              d3s(1,kk,3) = d3s(1,kk,3) - xrkk*sinq*f1z*fct
              d3s(2,kk,3) = d3s(2,kk,3) - yrkk*sinq*f1z*fct
              d3s(3,kk,3) = d3s(3,kk,3) + cosq*f2zz*fct
            enddo
          endif
        endif
!
!  Calculate Cartesian/strain/strain third derivatives
!
        if (ndim.eq.3) then
          do kk = 1,nstrains
            do kl = 1,nstrains
              d3ss(1,kl,kk) = d3ss(1,kl,kk) - xrkk*d2g2ds2(iv,kl,kk)*sins*fct
              d3ss(1,kl,kk) = d3ss(1,kl,kk) - xrkk*dg2ds(iv,kl)*dg2ds(iv,kk)*sins2*fct
              d3ss(1,kl,kk) = d3ss(1,kl,kk) - (dg2ds(iv,kl)*dgds(iv,1,kk) + dg2ds(iv,kk)*dgds(iv,1,kl))*sins*fct
              d3ss(1,kl,kk) = d3ss(1,kl,kk) - d2gds2(iv,1,kl,kk)*sinq*fct
!
              d3ss(2,kl,kk) = d3ss(2,kl,kk) - yrkk*d2g2ds2(iv,kl,kk)*sins*fct
              d3ss(2,kl,kk) = d3ss(2,kl,kk) - yrkk*dg2ds(iv,kl)*dg2ds(iv,kk)*sins2*fct
              d3ss(2,kl,kk) = d3ss(2,kl,kk) - (dg2ds(iv,kl)*dgds(iv,2,kk) + dg2ds(iv,kk)*dgds(iv,2,kl))*sins*fct
              d3ss(2,kl,kk) = d3ss(2,kl,kk) - d2gds2(iv,2,kl,kk)*sinq*fct
!
              d3ss(3,kl,kk) = d3ss(3,kl,kk) - zrkk*d2g2ds2(iv,kl,kk)*sins*fct
              d3ss(3,kl,kk) = d3ss(3,kl,kk) - zrkk*dg2ds(iv,kl)*dg2ds(iv,kk)*sins2*fct
              d3ss(3,kl,kk) = d3ss(3,kl,kk) - (dg2ds(iv,kl)*dgds(iv,3,kk) + dg2ds(iv,kk)*dgds(iv,3,kl))*sins*fct
              d3ss(3,kl,kk) = d3ss(3,kl,kk) - d2gds2(iv,3,kl,kk)*sinq*fct
            enddo
          enddo
!
!  Volume and volume mixed derivatives
!
          if (lfinitestrain) then
            do kk = 1,nstrains
              do kl = 1,nstrains
                d3ss(1,kl,kk) = d3ss(1,kl,kk) + xrkk*sins*fct*straindet*(dg2ds(iv,kl)*strainddetds(kk) + &
                                                                         dg2ds(iv,kk)*strainddetds(kl))
                d3ss(1,kl,kk) = d3ss(1,kl,kk) + sinq*fct*straindet*(dgds(iv,1,kl)*strainddetds(kk) + &
                                                                    dgds(iv,1,kk)*strainddetds(kl))
                d3ss(1,kl,kk) = d3ss(1,kl,kk) - xrkk*sinq*fct*(2.0_dp*strainddetds(kk)*strainddetds(kl)*straindet**2 - &
                                                               straind2detds2(kl,kk)*straindet)
!
                d3ss(2,kl,kk) = d3ss(2,kl,kk) + yrkk*sins*fct*straindet*(dg2ds(iv,kl)*strainddetds(kk) + &
                                                                         dg2ds(iv,kk)*strainddetds(kl))
                d3ss(2,kl,kk) = d3ss(2,kl,kk) + sinq*fct*straindet*(dgds(iv,2,kl)*strainddetds(kk) + &
                                                                    dgds(iv,2,kk)*strainddetds(kl))
                d3ss(2,kl,kk) = d3ss(2,kl,kk) - yrkk*sinq*fct*(2.0_dp*strainddetds(kk)*strainddetds(kl)*straindet**2 - &
                                                               straind2detds2(kl,kk)*straindet)
!
                d3ss(3,kl,kk) = d3ss(3,kl,kk) + zrkk*sins*fct*straindet*(dg2ds(iv,kl)*strainddetds(kk) + &
                                                                         dg2ds(iv,kk)*strainddetds(kl))
                d3ss(3,kl,kk) = d3ss(3,kl,kk) + sinq*fct*straindet*(dgds(iv,3,kl)*strainddetds(kk) + &
                                                                    dgds(iv,3,kk)*strainddetds(kl))
                d3ss(3,kl,kk) = d3ss(3,kl,kk) - zrkk*sinq*fct*(2.0_dp*strainddetds(kk)*strainddetds(kl)*straindet**2 - &
                                                               straind2detds2(kl,kk)*straindet)
              enddo
            enddo
          else
            do kk = 1,nstrains
              do kl = 1,nstrains
                if (kk.le.3) then
                  d3ss(1,kl,kk) = d3ss(1,kl,kk) + xrkk*dg2ds(iv,kl)*sins*fct
                  d3ss(1,kl,kk) = d3ss(1,kl,kk) + dgds(iv,1,kl)*sinq*fct
!
                  d3ss(2,kl,kk) = d3ss(2,kl,kk) + yrkk*dg2ds(iv,kl)*sins*fct
                  d3ss(2,kl,kk) = d3ss(2,kl,kk) + dgds(iv,2,kl)*sinq*fct
!
                  d3ss(3,kl,kk) = d3ss(3,kl,kk) + zrkk*dg2ds(iv,kl)*sins*fct
                  d3ss(3,kl,kk) = d3ss(3,kl,kk) + dgds(iv,3,kl)*sinq*fct
!
                  if (kl.le.3) then
                    d3ss(1,kl,kk) = d3ss(1,kl,kk) - xrkk*sinq*fct
                    d3ss(2,kl,kk) = d3ss(2,kl,kk) - yrkk*sinq*fct
                    d3ss(3,kl,kk) = d3ss(3,kl,kk) - zrkk*sinq*fct
                  endif
                endif
                if (kl.le.3) then
                  d3ss(1,kl,kk) = d3ss(1,kl,kk) + xrkk*dg2ds(iv,kk)*sins*fct
                  d3ss(1,kl,kk) = d3ss(1,kl,kk) + dgds(iv,1,kk)*sinq*fct
!
                  d3ss(2,kl,kk) = d3ss(2,kl,kk) + yrkk*dg2ds(iv,kk)*sins*fct
                  d3ss(2,kl,kk) = d3ss(2,kl,kk) + dgds(iv,2,kk)*sinq*fct
!
                  d3ss(3,kl,kk) = d3ss(3,kl,kk) + zrkk*dg2ds(iv,kk)*sins*fct
                  d3ss(3,kl,kk) = d3ss(3,kl,kk) + dgds(iv,3,kk)*sinq*fct
                endif
              enddo
            enddo
          endif
        elseif (ndim.eq.2) then
          do kk = 1,nstrains
            do kl = 1,nstrains
              d3ss(1,kl,kk) = d3ss(1,kl,kk) - sinq*fct*(f1g - f*rkvec**2)*(dg2ds(iv,kl)*dgds(iv,1,kk) + &
                                                                           dg2ds(iv,kk)*dgds(iv,1,kl) + &
                                                                           xrkk*d2g2ds2(iv,kl,kk))
              d3ss(1,kl,kk) = d3ss(1,kl,kk) - sinq*fct*(f2gg - 2.0_dp*f1g*rkvec**2 + 3.0_dp*f*rkvec**4)*&
                                                                           dg2ds(iv,kl)*dg2ds(iv,kk)*xrkk
              d3ss(1,kl,kk) = d3ss(1,kl,kk) - sinq*fct*f*d2gds2(iv,1,kl,kk)
!
              d3ss(2,kl,kk) = d3ss(2,kl,kk) - sinq*fct*(f1g - f*rkvec**2)*(dg2ds(iv,kl)*dgds(iv,2,kk) + &
                                                                           dg2ds(iv,kk)*dgds(iv,2,kl) + &
                                                                           yrkk*d2g2ds2(iv,kl,kk))
              d3ss(2,kl,kk) = d3ss(2,kl,kk) - sinq*fct*(f2gg - 2.0_dp*f1g*rkvec**2 + 3.0_dp*f*rkvec**4)*&
                                                                           dg2ds(iv,kl)*dg2ds(iv,kk)*yrkk
              d3ss(2,kl,kk) = d3ss(2,kl,kk) - sinq*fct*f*d2gds2(iv,2,kl,kk)
!
              d3ss(3,kl,kk) = d3ss(3,kl,kk) + cosq*fct*(f2zg - f1z*rkvec**2)*d2g2ds2(iv,kl,kk)
              d3ss(3,kl,kk) = d3ss(3,kl,kk) + cosq*fct*(f3zgg - 2.0_dp*f2zg*rkvec**2 + &
                                                        3.0_dp*f1z*rkvec**4)*dg2ds(iv,kl)*dg2ds(iv,kk)
            enddo
          enddo
!
!  Area and area mixed derivatives
!
          if (lfinitestrain) then
            do kk = 1,nstrains
              do kl = 1,nstrains
                d3ss(1,kl,kk) = d3ss(1,kl,kk) + sinq*fct*f*straindet*(dgds(iv,1,kl)*strainddetds(kk) + &
                                                                      dgds(iv,1,kk)*strainddetds(kl))
                d3ss(1,kl,kk) = d3ss(1,kl,kk) + xrkk*sinq*fct*(f1g - f*rkvec**2)*straindet*( &
                                                dg2ds(iv,kl)*strainddetds(kk) + dg2ds(iv,kk)*strainddetds(kl))
                d3ss(1,kl,kk) = d3ss(1,kl,kk) - sinq*fct*xrkk*f*(2.0_dp*strainddetds(kk)*strainddetds(kl)*straindet**2 - &
                                                               straind2detds2(kl,kk)*straindet)
!
                d3ss(2,kl,kk) = d3ss(2,kl,kk) + sinq*fct*f*straindet*(dgds(iv,2,kl)*strainddetds(kk) + &
                                                                      dgds(iv,2,kk)*strainddetds(kl))
                d3ss(2,kl,kk) = d3ss(2,kl,kk) + yrkk*sinq*fct*(f1g - f*rkvec**2)*straindet*( &
                                                dg2ds(iv,kl)*strainddetds(kk) + dg2ds(iv,kk)*strainddetds(kl))
                d3ss(2,kl,kk) = d3ss(2,kl,kk) - sinq*fct*yrkk*f*(2.0_dp*strainddetds(kk)*strainddetds(kl)*straindet**2 - &
                                                               straind2detds2(kl,kk)*straindet)
!
                d3ss(3,kl,kk) = d3ss(3,kl,kk) - cosq*fct*(f2zg - f1z*rkvec**2)*straindet*( &
                                                dg2ds(iv,kl)*strainddetds(kk) + dg2ds(iv,kk)*strainddetds(kl))
                d3ss(3,kl,kk) = d3ss(3,kl,kk) + cosq*fct*f1z*(2.0_dp*strainddetds(kk)*strainddetds(kl)*straindet**2 - &
                                                               straind2detds2(kl,kk)*straindet)
              enddo
            enddo
          else
            do kk = 1,nstrains
              do kl = 1,nstrains
                if (kk.le.2) then
                  d3ss(1,kl,kk) = d3ss(1,kl,kk) + sinq*fct*(f*dgds(iv,1,kl) + xrkk*(f1g - f*rkvec**2)*dg2ds(iv,kl))
                  d3ss(2,kl,kk) = d3ss(2,kl,kk) + sinq*fct*(f*dgds(iv,2,kl) + yrkk*(f1g - f*rkvec**2)*dg2ds(iv,kl))
                  d3ss(3,kl,kk) = d3ss(3,kl,kk) - cosq*fct*(f2zg - f1z*rkvec**2)*dg2ds(iv,kl)
                  if (kl.le.2) then
                    d3ss(1,kl,kk) = d3ss(1,kl,kk) - sinq*fct*xrkk*f
                    d3ss(2,kl,kk) = d3ss(2,kl,kk) - sinq*fct*yrkk*f
                    d3ss(3,kl,kk) = d3ss(3,kl,kk) + cosq*fct*f1z
                  endif
                endif
                if (kl.le.2) then
                  d3ss(1,kl,kk) = d3ss(1,kl,kk) + sinq*fct*(f*dgds(iv,1,kk) + xrkk*(f1g - f*rkvec**2)*dg2ds(iv,kk))
                  d3ss(2,kl,kk) = d3ss(2,kl,kk) + sinq*fct*(f*dgds(iv,2,kk) + yrkk*(f1g - f*rkvec**2)*dg2ds(iv,kk))
                  d3ss(3,kl,kk) = d3ss(3,kl,kk) - cosq*fct*(f2zg - f1z*rkvec**2)*dg2ds(iv,kk)
                endif
              enddo
            enddo
          endif
        endif
      endif
    endif
  enddo
!
!  Symmetrise d2
!
  d2(1,2) = d2(2,1)
  d2(1,3) = d2(3,1)
  d2(2,3) = d2(3,2)
!
  if (lgrad3) then
!
!  Symmetrise d3
!
    d3(1,1,2) = d3(2,1,1)
    d3(1,2,1) = d3(2,1,1)
    d3(1,1,3) = d3(3,1,1)
    d3(1,3,1) = d3(3,1,1)
    d3(1,2,2) = d3(2,2,1)
    d3(2,1,2) = d3(2,2,1)
    d3(2,3,1) = d3(3,2,1)
    d3(1,3,2) = d3(3,2,1)
    d3(3,1,2) = d3(3,2,1)
    d3(1,2,3) = d3(3,2,1)
    d3(2,1,3) = d3(3,2,1)
    d3(1,3,3) = d3(3,3,1)
    d3(3,1,3) = d3(3,3,1)
    d3(2,2,3) = d3(3,2,2)
    d3(2,3,2) = d3(3,2,2)
    d3(2,3,3) = d3(3,3,2)
    d3(3,2,3) = d3(3,3,2)
  endif
!
!  Exit point
!
999 continue
#ifdef TRACE
  call trace_out('reciptrmdp')
#endif
!
  return
  end
