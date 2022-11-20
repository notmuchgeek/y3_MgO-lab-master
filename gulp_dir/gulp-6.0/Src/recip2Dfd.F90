  subroutine recip2Dfd(matom,erecip,lgrad1)
!
!  Calculation of electrostatic potential and derivatives,
!  including strain derivatives. Reciprocal space part.
!  Two-dimensional version using the Parry sum.
!  Finite difference version that focuses on derivatives of matom.
!
!  12/17 Created from recip2D
!   2/18 Trace added
!   9/18 Handling of lstraincell algorithm added
!   9/18 Strain module introduced
!  11/18 Finite strain flag introduced instead of lstraincell
!  12/19 Rigid molecules added
!   1/20 Correction to com coordinate setting
!   3/20 Use of charge pointer added
!   3/20 Rigid molecule correction added
!   4/20 d2xyzdsdc added to cartstrterm arguments
!   4/20 derv3c changes reversed as they are no longer required
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
!  Julian Gale, CIC, Curtin University, April 2020
!
  use g_constants
  use control
  use current
  use derivatives
  use kspace
  use m_strain,       only : gstrterms, strainddetds, straindet, cartstrterm, gxyzterms
  use molecule
  use optimisation
  use parallel
  use polarise
  use potentialxyz
  use symmetry
  use thresholds,     only : thresh_q
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  use two
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                        :: matom
  real(dp),    intent(inout)                     :: erecip
  logical,     intent(in)                        :: lgrad1
!
!  Local variables
!
  integer(i4)                                    :: i
  integer(i4)                                    :: idk
  integer(i4)                                    :: ii
  integer(i4)                                    :: imin
  integer(i4)                                    :: iv
  integer(i4)                                    :: j
  integer(i4)                                    :: jj
  integer(i4)                                    :: kk
  integer(i4)                                    :: nmi
  integer(i4)                                    :: nmj
  logical                                        :: lsg1
  real(dp)                                       :: accf2
  real(dp)                                       :: arg
  real(dp)                                       :: argtest
  real(dp)                                       :: cosa
  real(dp)                                       :: cosq
  real(dp)                                       :: g_cpu_time
  real(dp)                                       :: csinq
  real(dp)                                       :: d0term
  real(dp)                                       :: darg1
  real(dp)                                       :: darg2
  real(dp)                                       :: g_derf
  real(dp)                                       :: g_derfc
  real(dp)                                       :: derfc1
  real(dp)                                       :: derfc2
  real(dp)                                       :: derfez
  real(dp)                                       :: dexp1
  real(dp)                                       :: dexp2
  real(dp)                                       :: dexp3
  real(dp)                                       :: dexp4
  real(dp)                                       :: dexpz
  real(dp)                                       :: dtrm1
  real(dp)                                       :: dGrds
  real(dp)                                       :: drxyzds(6,3)
  real(dp)                                       :: d2rxyzdsdx(6,3,3)
  real(dp)                                       :: d2rxyzds2(6,6,3)
  real(dp)                                       :: etaz
  real(dp)                                       :: etaz2
  real(dp)                                       :: eztrm
  real(dp)                                       :: factor
  real(dp)                                       :: fct
  real(dp)                                       :: Gmax
  real(dp)                                       :: Gmin
  real(dp)                                       :: gseta
  real(dp)                                       :: kexperfc
  real(dp)                                       :: kvec
  real(dp)                                       :: kvv(3)
  real(dp)                                       :: oci
  real(dp)                                       :: ocj
  real(dp)                                       :: qfct
  real(dp)                                       :: qli
  real(dp)                                       :: qlj
  real(dp)                                       :: rk2
  real(dp)                                       :: rktemp
  real(dp)                                       :: rkvec
  real(dp)                                       :: sina
  real(dp)                                       :: sineq
  real(dp)                                       :: sinqx
  real(dp)                                       :: sinqy
  real(dp)                                       :: sinqz
  real(dp)                                       :: smallestG
  real(dp)                                       :: strm1
  real(dp)                                       :: time0
  real(dp)                                       :: time1
  real(dp)                                       :: twoqv
  real(dp)                                       :: xci
  real(dp)                                       :: yci
  real(dp)                                       :: zci
  real(dp)                                       :: xcomi
  real(dp)                                       :: ycomi
  real(dp)                                       :: zcomi
  real(dp)                                       :: xcom
  real(dp)                                       :: ycom
  real(dp)                                       :: zcom
  real(dp)                                       :: xd
  real(dp)                                       :: yd
  real(dp)                                       :: zd
  real(dp)                                       :: xrkk
  real(dp)                                       :: yrkk
  real(dp)                                       :: ztrm1
#ifdef TRACE
  call trace_in('recip2Dfd')
#endif
!
  time0 = g_cpu_time()
  lsg1 = (lstr.and.lgrad1)
  accf2 = accf*accf
  argtest = sqrt(3.0_dp+0.5_dp*accf2) - sqrt(3.0_dp)
  smallestG = min(abs(kv(1,1)),abs(kv(2,2)))
!**********
!  Setup  *
!**********
  if (lra) then
    kvv(1) = kv(1,1)
    kvv(2) = kv(2,2)
    do iv = 1,nkvec
      idk = indk(iv)
      ii = (idk/maxindk3) - maxindk
      if (ii.eq.0) then
        factor = 1.0_dp
      else
        factor = 2.0_dp
      endif
      idk = idk - (ii+maxindk)*maxindk3
      jj = (idk/maxindk2) - maxindk
      xrk(iv) = ii*kvv(1)
      yrk(iv) = jj*kvv(2)
      rk2 = xrk(iv)*xrk(iv) + yrk(iv)*yrk(iv)
      kmod(iv) = sqrt(rk2)
      ktrm(iv) = 0.5_dp*vol4pi*factor/kmod(iv)
    enddo
  else
    do iv = 1,nkvec
      idk = indk(iv)
      ii = (idk/maxindk3) - maxindk
      idk = idk - (ii+maxindk)*maxindk3
      jj = (idk/maxindk2) - maxindk
      factor = 2.0_dp
      if (ii.eq.0.and.nkangle.eq.1) then
        factor = 1.0_dp
      elseif (nkangle.eq.0) then
        factor = 1.0_dp
      endif
      xrk(iv) = ii*kv(1,1) + jj*kv(1,2)
      yrk(iv) = ii*kv(2,1) + jj*kv(2,2)
      rk2 = xrk(iv)*xrk(iv) + yrk(iv)*yrk(iv)
      kmod(iv) = sqrt(rk2)
      ktrm(iv) = 0.5_dp*vol4pi*factor/kmod(iv)
    enddo
  endif
!*******************************************
!  Sort K vectors by increasing magnitude  *
!*******************************************
  do i = 1,nkvec
    Gmin = 2.0_dp*rradmax
    do j = i,nkvec
      if (kmod(j).lt.Gmin) then
        imin = j
        Gmin = kmod(j)
      endif
    enddo
    rktemp = kmod(i)
    kmod(i) = kmod(imin)
    kmod(imin) = rktemp
    rktemp = ktrm(i)
    ktrm(i) = ktrm(imin)
    ktrm(imin) = rktemp
    rktemp = xrk(i)
    xrk(i) = xrk(imin)
    xrk(imin) = rktemp
    rktemp = yrk(i)
    yrk(i) = yrk(imin)
    yrk(imin) = rktemp
  enddo
!**************************
!  End of set-up section  *
!**************************
  if (lnorecip) goto 999
!
!  Define constants
!
  rpieta = 1.0_dp/sqrt(pi * eta)
  rhseta = 0.5_dp/seta
!
!  Build products of K vector components
!
  if (lsg1) then
    call gstrterms(ndim,maxkvec,nkvec,xrk,yrk,zrk,dg2ds,d2g2dx2,d2g2ds2,.false.)
    if (lrigid) then
      call gxyzterms(ndim,maxkvec,nkvec,xrk,yrk,zrk,dgds,d2gds2,.false.)
    endif
  endif
!
!  Only do i = matom case
!
  i = matom
!
  oci = occuf(i)*angstoev
  qli = qf(i)*oci
!
!  Check charge magnitude is worth computing
!
  if (abs(qli).gt.thresh_q) then
    xci = xclat(i)
    yci = yclat(i)
    zci = zclat(i)
!
    nmi = natmol(i)
    if (lrigid.and.nmi.gt.0) then
      xcomi = molxyz(1,natinmol(i),nmi)
      ycomi = molxyz(2,natinmol(i),nmi)
      zcomi = molxyz(3,natinmol(i),nmi)
    else
      xcomi = 0.0_dp
      ycomi = 0.0_dp
      zcomi = 0.0_dp
    endif
!
    jloop: do jj = 1,ncharge
      j = nchargeptr(jj)
      ocj = occuf(j)
      if (i.eq.j) then
        fct = 0.5_dp
      else 
        fct = 1.0_dp
      endif
      qlj = qf(j)*ocj*fct
      if (i.ne.j) then
        nmj = natmol(j)
!
        if (lrigid) then
          if (nmj.gt.0) then
            xcom = molxyz(1,natinmol(j),nmj) - xcomi
            ycom = molxyz(2,natinmol(j),nmj) - ycomi
            zcom = molxyz(3,natinmol(j),nmj) - zcomi
          else
            xcom = - xcomi
            ycom = - ycomi
            zcom = - zcomi
          endif
        else
          xcom = - xcomi
          ycom = - ycomi
          zcom = - zcomi
        endif
!
!  Find relative vector between atoms
!
        xd = xclat(j) - xci
        yd = yclat(j) - yci
        zd = zclat(j) - zci
!
        qfct = qli*qlj
        etaz = seta*zd
        etaz2 = etaz * etaz
!
        if (lrigid.and.lsg1) then
          call cartstrterm(ndim,xd,yd,zd,xcom,ycom,zcom,drxyzds,d2rxyzdsdx,d2rxyzds2,.false.)
        endif
!
!  First term - K vector independent
!
        derfez = g_derf(etaz)
        dexpz  = exp(-etaz2)
        twoqv  = qfct*vol4pi
        eztrm  = twoqv*(zd*derfez + dexpz*rpieta)
        erecip = erecip - eztrm
        if (lgrad1) then
          dtrm1 = - twoqv*derfez
        endif
!
!  Add in charge derivatives
!
        if (lgrad1.and.lDoQDeriv1) then
          d0term = - vol4pi*(zd*derfez + dexpz*rpieta)*fct*angstoev
          call d1charge(i,j,.true.,.true.,1_i4,d0term*qf(j),d0term*qf(i))
        endif
!
!  Find local kvector cut-off
!
        if (abs(etaz).gt.argtest) then
          Gmax = abs(accf2/zd)
        else
          Gmax = sqrt(4.0_dp*eta*(accf2-etaz2))
        endif
!
!  Second term - K vector dependent
!
        csinq = 0.0_dp
        if (lgrad1) then
          sinqx = 0.0_dp
          sinqy = 0.0_dp
          sinqz = 0.0_dp
          if (Gmax.ge.smallestG) then
            iv = 1
            kvec = kmod(iv)
            do while (iv.le.nkvec.and.kvec.le.Gmax)
              xrkk = xrk(iv)
              yrkk = yrk(iv)
              arg = xrkk*xd + yrkk*yd
              cosa = cos(arg)
              sina = sin(arg)*qfct
              cosq = cosa*ktrm(iv)
              sineq = sina*ktrm(iv)
              dexp1 = exp(kvec*zd)
              dexp2 = 1.0_dp/dexp1
              darg1 = kvec*rhseta + etaz
              darg2 = kvec*rhseta - etaz
              dexp3 = exp(-(darg1)**2)
              dexp4 = exp(-(darg2)**2)
              derfc1 = g_derfc(darg1)
              derfc2 = g_derfc(darg2)
              kexperfc = dexp1*derfc1 + dexp2*derfc2
!
!  Energy
!
              csinq = csinq + cosq*kexperfc
!
!  First derivatives with respect to atoms
!
              cosq = cosq*qfct
              sinqx = sinqx - sineq*xrkk*kexperfc
              sinqy = sinqy - sineq*yrkk*kexperfc
              ztrm1 = (kvec*(dexp1*derfc1-dexp2*derfc2) - tweatpi*(dexp1*dexp3-dexp2*dexp4))
              sinqz = sinqz + cosq*ztrm1
              if (lsg1) then
!
!  Strain first derivatives
!
                rkvec = 1.0_dp/kvec
                strm1 = rkvec*(-rkvec*kexperfc + zd*(dexp1*derfc1-dexp2*derfc2) - rpieta*(dexp1*dexp3+dexp2*dexp4))
                do kk = 1,nstrains
                  strderv(kk) = strderv(kk) + strm1*cosq*dg2ds(iv,kk)
                enddo
              endif
              iv = iv + 1
              kvec = kmod(iv)
            enddo
          endif
        else
          if (Gmax.ge.smallestG) then
            iv = 1
            kvec = kmod(iv)
            do while (iv.le.nkvec.and.kvec.le.Gmax)
              arg = xrk(iv)*xd + yrk(iv)*yd
              cosa = cos(arg)
              dexp1 = exp(kvec*zd)
              dexp2 = 1.0_dp/dexp1
              gseta = kvec*rhseta
              kexperfc = dexp1*g_derfc(gseta+etaz) + dexp2*g_derfc(gseta-etaz)
              csinq = csinq + cosa*ktrm(iv)*kexperfc
              iv = iv + 1
              kvec = kmod(iv)
            enddo
          endif
        endif
      else
        qfct = qli*qlj
!
!  First term - K vector independent
!
        twoqv = qfct*vol4pi
        erecip = erecip - twoqv*rpieta
!
!  Add in charge derivatives
!
        if (lgrad1.and.lDoQDeriv1) then
          d0term = - vol4pi*rpieta*fct*angstoev
          call d1charge(i,j,.true.,.true.,1_i4,d0term*qf(j),d0term*qf(i))
        endif
!
!  Second term - K vector dependent
!
        csinq = 0.0_dp
        if (lgrad1) then
          sinqx = 0.0_dp
          sinqy = 0.0_dp
          sinqz = 0.0_dp
          do iv = 1,nkvec
            xrkk = xrk(iv)
            yrkk = yrk(iv)
            cosq = ktrm(iv)
            kvec = kmod(iv)
            darg1 = kvec*rhseta
            dexp1 = exp(-(darg1)**2)
            derfc1 = g_derfc(darg1)
            kexperfc = 2.0_dp*derfc1
!
!  Energy
!
            csinq = csinq + cosq*kexperfc
!
            cosq = cosq*qfct
            if (lsg1) then
!
!  Strain first derivatives
!
              rkvec = 1.0_dp/kvec
              strm1 = rkvec*(-rkvec*kexperfc - 2.0_dp*rpieta*dexp1)
              do kk = 1,nstrains
                strderv(kk) = strderv(kk) + strm1*cosq*dg2ds(iv,kk)
              enddo
!
              if (lrigid) then
                do kk = 1,nstrains
                  dGrds = xrkk*drxyzds(kk,1) + yrkk*drxyzds(kk,2) + &
                          xd*dgds(iv,1,kk) + yd*dgds(iv,2,kk)
                  strderv(kk) = strderv(kk) - sineq*dGrds*kexperfc
                enddo
              endif
            endif
          enddo
        else
          do iv = 1,nkvec
            kvec = kmod(iv)
            gseta = kvec*rhseta
            kexperfc = 2.0_dp*g_derfc(gseta)
            csinq = csinq + ktrm(iv)*kexperfc
          enddo
        endif
      endif
!
!  Lattice energy
!
      erecip = erecip + csinq*qfct
!
!  Add in charge derivatives
!
      if (lgrad1.and.lDoQDeriv1) then
        d0term = csinq*fct*angstoev
        call d1charge(i,j,.true.,.true.,1_i4,d0term*qf(j),d0term*qf(i))
      endif
!
!  Internal derivatives
!
      if (lgrad1.and.i.ne.j) then
        xdrv(i) = xdrv(i) - sinqx
        xdrv(j) = xdrv(j) + sinqx
        ydrv(i) = ydrv(i) - sinqy
        ydrv(j) = ydrv(j) + sinqy
!
!  Because dtrm1 doesn't depend on kvector only do on 1 processor
!
        if (procid.eq.0) then
          zdrv(i) = zdrv(i) - sinqz - dtrm1
          zdrv(j) = zdrv(j) + sinqz + dtrm1
        else
          zdrv(i) = zdrv(i) - sinqz
          zdrv(j) = zdrv(j) + sinqz
        endif
      endif
    enddo jloop
!****************************************************
!  Complete strain derivatives in reciprocal space  *
!****************************************************
    if (lsg1) then
      if (lfinitestrain) then
        strderv(1) = strderv(1) - strainddetds(1)*straindet*erecip
        strderv(2) = strderv(2) - strainddetds(2)*straindet*erecip
        strderv(3) = strderv(3) - strainddetds(3)*straindet*erecip
      else
        strderv(1) = strderv(1) - erecip
        strderv(2) = strderv(2) - erecip
      endif
    endif
  endif
!***************
!  Exit point  *
!***************
999 continue
!
!  Timing
!
  time1 = g_cpu_time()
  tion = tion + time1 - time0
#ifdef TRACE
  call trace_out('recip2Dfd')
#endif
!
  return
  end
