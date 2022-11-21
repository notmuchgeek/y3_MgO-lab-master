  subroutine recip3Dp(xkv,ykv,zkv)
!
!  Calculate phased second derivatives in reciprocal space
!  for phonon calculations.
!
!  The use of right angled saving must be disabled in kindex when
!  calling this routine as the k point breaks the mirror plane
!  symmetry.
!
!  NOTE : outer loop over atom algorithm has been removed in
!         this version and needs restoring with gamma point
!         subtraction added
!
!   8/95 Ewald sum for dispersion terms added
!   8/97 Set up of K vectors removed as should be done in advance
!        by kindex
!   9/97 Diagonal blocks of dynamical matrix set to be difference
!        between k point and gamma point
!   8/98 Derv2/dervi storage switched to normal triangle and
!        call to d2chargep added for EEM/QEq
!   8/98 ESFF Lennard-Jones form now allowed for
!   8/98 lfirst parameter removed as call to recipp should
!        always be proceeded by call to recip
!   7/00 Dynamic memory allocation added
!   4/01 Passed argument changed from K point number to K point 
!        coordinate
!  11/02 K vector now passed d2chargep
!  12/02 nkaddx/nkaddy introduced to handle extreme angles
!   9/04 Modifications for generalisation of variable charges added
!  11/04 Sqrt pi taken from module
!  12/07 Unused variables removed
!   3/14 derfc changed to g_derfc for benefit of ChemShell
!   8/14 Group velocities added as an option
!   9/14 Derivatives of ktrm added for group velocities
!  10/14 Symmetrisation of reciprocal space second derivatives added
!   2/15 MM3buck added
!   2/18 Trace added
!   3/20 Use of charge pointer added
!   3/20 Tolerance for ldoc6 made global
!   5/20 Rigid molecule changes added
!   5/20 Phasing based on centre of mass added
!   5/20 Bugs in group velocities fixed
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
!  Julian Gale, CIC, Curtin University, July 2020
!
  use g_constants
  use control
  use current
  use derivatives
  use ksample
  use kspace
  use molecule
  use parallel
  use shells
  use symmetry
  use thresholds,     only : thresh_c6, thresh_q
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  use two
  implicit none
!
!  Passed variables
!
  real(dp),     intent(in)                    :: xkv
  real(dp),     intent(in)                    :: ykv
  real(dp),     intent(in)                    :: zkv
!
!  Local variables
!
  integer(i4)                                 :: i
  integer(i4)                                 :: ii
  integer(i4)                                 :: iv
  integer(i4)                                 :: ix
  integer(i4)                                 :: iy
  integer(i4)                                 :: iz
  integer(i4)                                 :: j
  integer(i4)                                 :: jj
  integer(i4)                                 :: jx
  integer(i4)                                 :: jy
  integer(i4)                                 :: jz
  integer(i4)                                 :: kk
  integer(i4)                                 :: max1l
  integer(i4)                                 :: max1l1
  integer(i4)                                 :: max1u
  integer(i4)                                 :: max2l
  integer(i4)                                 :: max2l1
  integer(i4)                                 :: max2u
  integer(i4)                                 :: max3l
  integer(i4)                                 :: max3l1
  integer(i4)                                 :: max3u
  integer(i4)                                 :: n
  integer(i4)                                 :: nati
  integer(i4)                                 :: natj
  integer(i4)                                 :: ni
  integer(i4)                                 :: nj
  integer(i4)                                 :: nkaddx
  integer(i4)                                 :: nkaddy
  integer(i4)                                 :: nkvec0
  integer(i4)                                 :: nmi
  integer(i4)                                 :: nmj
  integer(i4)                                 :: ns
  integer(i4)                                 :: ntypj
  integer(i4)                                 :: ntypi
  integer(i4)                                 :: status
  logical                                     :: lc6loc
  logical                                     :: ldoc6
  logical                                     :: lgamma
  complex(dpc)                                :: cdk(3,6)
  real(dp)                                    :: anglemax
  real(dp)                                    :: anglemin
  real(dp)                                    :: arg
  real(dp)                                    :: arge
  real(dp)                                    :: c6t1
  real(dp)                                    :: c6t2
  real(dp)                                    :: c6t3
  real(dp)                                    :: c6t4
  real(dp)                                    :: c6t2dk
  real(dp)                                    :: c6t3dk
  real(dp)                                    :: c6t4dk
  real(dp)                                    :: c6tot
  real(dp)                                    :: cos6
  real(dp)                                    :: cosa
  real(dp)                                    :: cosq
  real(dp)                                    :: g_cpu_time
  real(dp)                                    :: d1ix
  real(dp)                                    :: d1iy
  real(dp)                                    :: d1iz
  real(dp)                                    :: d1jx
  real(dp)                                    :: d1jy
  real(dp)                                    :: d1jz
  real(dp)                                    :: d2trmr
  real(dp)                                    :: d2trmi
  real(dp)                                    :: d2trmx
  real(dp)                                    :: d2trmy
  real(dp)                                    :: d2trmz
  real(dp)                                    :: d2qi(6)
  real(dp)                                    :: d2qr(6)
  real(dp)                                    :: d2self
  real(dp)                                    :: dkdk
  real(dp)                                    :: g_derfc
  real(dp)                                    :: fct
  real(dp), dimension(:,:), allocatable       :: ktrmdk
  real(dp), dimension(:,:), allocatable       :: ktrm6dk
  real(dp), dimension(:),   allocatable       :: ktrm0
  real(dp), dimension(:),   allocatable       :: ktrm3
  real(dp), dimension(:),   allocatable       :: ktrm6
  real(dp), dimension(:),   allocatable       :: ktrm60
  real(dp)                                    :: oci
  real(dp)                                    :: ocj
  real(dp)                                    :: perpk
  real(dp)                                    :: projk
  real(dp)                                    :: qfct
  real(dp)                                    :: qli
  real(dp)                                    :: qlj
  real(dp)                                    :: rk
  real(dp)                                    :: rk2
  real(dp)                                    :: rk20
  real(dp)                                    :: rketa2
  real(dp)                                    :: rkv
  real(dp)                                    :: rkv2
  real(dp)                                    :: rk1x
  real(dp)                                    :: rk1y
  real(dp)                                    :: rk1z
  real(dp)                                    :: rk2x
  real(dp)                                    :: rk2y
  real(dp)                                    :: rk2z
  real(dp)                                    :: rk3x
  real(dp)                                    :: rk3y
  real(dp)                                    :: rk3z
  real(dp)                                    :: rkk1
  real(dp)                                    :: rkk2
  real(dp)                                    :: rkk3
  real(dp)                                    :: rlmx
  real(dp)                                    :: rlmx2
  real(dp)                                    :: rlmxy
  real(dp)                                    :: rrk2
  real(dp)                                    :: ruk1x
  real(dp)                                    :: ruk1y
  real(dp)                                    :: ruk1z
  real(dp)                                    :: ruk2x
  real(dp)                                    :: ruk2y
  real(dp)                                    :: ruk2z
  real(dp)                                    :: ruk3x
  real(dp)                                    :: ruk3y
  real(dp)                                    :: ruk3z
  real(dp)                                    :: sin6
  real(dp)                                    :: sina
  real(dp)                                    :: sinq
  real(dp)                                    :: time0
  real(dp)                                    :: time1
  real(dp), dimension(:,:), allocatable       :: tmp6
  real(dp)                                    :: trmc
  real(dp)                                    :: xci
  real(dp)                                    :: yci
  real(dp)                                    :: zci
  real(dp)                                    :: xcomi
  real(dp)                                    :: ycomi
  real(dp)                                    :: zcomi
  real(dp)                                    :: xcom
  real(dp)                                    :: ycom
  real(dp)                                    :: zcom
  real(dp)                                    :: xd
  real(dp)                                    :: yd
  real(dp)                                    :: zd
  real(dp)                                    :: xkd
  real(dp)                                    :: ykd
  real(dp)                                    :: zkd
  real(dp)                                    :: xke
  real(dp)                                    :: yke
  real(dp)                                    :: zke
  real(dp)                                    :: xkf
  real(dp)                                    :: ykf
  real(dp)                                    :: zkf
  real(dp)                                    :: xrkk
  real(dp)                                    :: yrkk
  real(dp)                                    :: zrkk
  real(dp)                                    :: xpon
  real(dp)                                    :: xxk
  real(dp)                                    :: xxk2
  real(dp)                                    :: yyk
  real(dp)                                    :: yyk2
  real(dp)                                    :: zzk
  real(dp)                                    :: zero
!
  if (lnorecip) return
!
#ifdef TRACE
  call trace_in('recip3Dp')
#endif
  time0 = g_cpu_time()
  lc6loc = (lc6.and.ndim.eq.3)
!
!  Initialise local variables
!
  zero = 1.0d-10
!
!  Assign local variables
!
  rk1x = kv(1,1)
  rk1y = kv(2,1)
  rk1z = kv(3,1)
  rk2x = kv(1,2)
  rk2y = kv(2,2)
  rk2z = kv(3,2)
  rk3x = kv(1,3)
  rk3y = kv(2,3)
  rk3z = kv(3,3)
  rkk1 = rk1x*rk1x + rk1y*rk1y + rk1z*rk1z
  rkk2 = rk2x*rk2x + rk2y*rk2y + rk2z*rk2z
  rkk3 = rk3x*rk3x + rk3y*rk3y + rk3z*rk3z
  rkk1 = sqrt(rkk1)
  rkk1 = 1.0_dp/rkk1
  rkk2 = sqrt(rkk2)
  rkk2 = 1.0_dp/rkk2
  rkk3 = sqrt(rkk3)
  rkk3 = 1.0_dp/rkk3
  ruk1x = rkk1*rk1x
  ruk1y = rkk1*rk1y
  ruk1z = rkk1*rk1z
  ruk2x = rkk2*rk2x
  ruk2y = rkk2*rk2y
  ruk2z = rkk2*rk2z
  ruk3x = rkk3*rk3x
  ruk3y = rkk3*rk3y
  ruk3z = rkk3*rk3z
!
!  Set gamma point flag
!
  lgamma = ((abs(xkv)+abs(ykv)+abs(zkv)).lt.1.0d-8)
!
  rkv2 = xkv*xkv + ykv*ykv + zkv*zkv
  rkv = sqrt(rkv2)
  rlmx = rkv + rradmax
  rlmx2 = rlmx*rlmx
  nkvec = 0
  nkvec0 = 0
!*************************
!  Find valid k vectors  *
!*************************
!
!  Determine maximum looping indices, allowing for k point
!
  if (lra) then
    max1u = (rlmx-xkv)*rkk1 + 1
    max1l = (rlmx+xkv)*rkk1 + 1
    max1l1 = max1l + 1
    xxk = xkv - max1l1*rk1x
    do ii = - max1l,max1u
      xxk = xxk + rk1x
      if (abs(xxk).lt.rlmx) then
        xxk2 = xxk*xxk
        rlmxy = rlmx2 - xxk2
        rlmxy = sqrt(rlmxy)
        max2u = (rlmxy-ykv)*rkk2 + 1
        max2l = (rlmxy+ykv)*rkk2 + 1
        max2l1 = max2l + 1
!
        yyk = ykv - max2l1*rk2y
        do jj = - max2l,max2u
          yyk = yyk + rk2y
!
          yyk2 = yyk*yyk
          perpk = rlmx2 - xxk2 - yyk2
          if (perpk.ge.0.0_dp) then
            perpk = sqrt(perpk)
            max3u = (perpk-zkv)*rkk3 + 1
            max3l = (perpk+zkv)*rkk3 + 1
            max3l1 = max3l + 1
!
            zzk = zkv - max3l1*rk3z
            do kk = - max3l,max3u
              zzk = zzk + rk3z
              rk2 = xxk2 + yyk2 + zzk*zzk
!
!  Test for zero wavevector and exclude vectors outside the search
!  region, given by rradmx.
!
              if (rk2.gt.zero.and.rk2.le.rlmx2) then
                nkvec = nkvec + 1
                if (nkvec.gt.maxkvec) then
                  maxkvec = nkvec + 100
                  call changemaxkvec
                endif
                xrk(nkvec) = xxk
                yrk(nkvec) = yyk
                zrk(nkvec) = zzk
              endif
              xkd = xxk - xkv
              ykd = yyk - ykv
              zkd = zzk - zkv
              rk20 = xkd*xkd + ykd*ykd + zkd*zkd
              if (rk20.gt.zero.and.rk20.le.rlmx2) then
                nkvec0 = nkvec0 + 1
                if (nkvec0.gt.maxkvec) then
                  maxkvec = nkvec0 + 100
                  call changemaxkvec
                endif
                xrk0(nkvec0) = xkd
                yrk0(nkvec0) = ykd
                zrk0(nkvec0) = zkd
              endif
!
!  End of ii,jj,kk loops over reciprocal lattice vectors
!
            enddo
          endif
        enddo
      endif
    enddo
  else
!
!  Set amounts to add to looping indices
!
    anglemax = max(alpha,beta)
    anglemax = max(anglemax,gamma)
    anglemin = min(alpha,beta)
    anglemin = min(anglemin,gamma)
    if (anglemax.gt.150.0_dp.or.anglemin.lt.30.0_dp) then
      nkaddx = 8
      nkaddy = 8
    elseif (anglemax.gt.140.0_dp.or.anglemin.lt.40.0_dp) then
      nkaddx = 7
      nkaddy = 7
    elseif (anglemax.gt.130.0_dp.or.anglemin.lt.50.0_dp) then
      nkaddx = 6
      nkaddy = 6
    elseif (anglemax.gt.120.0_dp.or.anglemin.lt.60.0_dp) then
      nkaddx = 5
      nkaddy = 5
    elseif (anglemax.gt.110.0_dp.or.anglemin.lt.70.0_dp) then
      nkaddx = 4
      nkaddy = 4
    elseif (anglemax.gt.100.0_dp.or.anglemin.lt.80.0_dp) then
      nkaddx = 3
      nkaddy = 3
    else
      nkaddx = 2
      nkaddy = 2
    endif
    projk = xkv*ruk1x + ykv*ruk1y + zkv*ruk1z
    max1u = (rlmx-projk)*rkk1 + nkaddx
    max1l = (rlmx+projk)*rkk1 + nkaddx
    max1l1 = max1l + 1
!
    xke = xkv - max1l1*rk1x
    yke = ykv - max1l1*rk1y
    zke = zkv - max1l1*rk1z
    do ii = -max1l,max1u
      xke = xke + rk1x
      yke = yke + rk1y
      zke = zke + rk1z
!
      projk = xke*ruk2x + yke*ruk2y + zke*ruk2z
      max2u = (rlmx-projk)*rkk2 + nkaddy
      max2l = (rlmx+projk)*rkk2 + nkaddy
      max2l1 = max2l + 1
!
      xkf = xke - max2l1*rk2x
      ykf = yke - max2l1*rk2y
      zkf = zke - max2l1*rk2z
      do jj = -max2l,max2u
        xkf = xkf + rk2x
        ykf = ykf + rk2y
        zkf = zkf + rk2z
!
        projk = xkf*ruk3x + ykf*ruk3y + zkf*ruk3z
        perpk = xkf*xkf + ykf*ykf + zkf*zkf
        perpk = perpk - projk*projk
        perpk = rlmx2 - perpk
        if (perpk.ge.0.0_dp) then
          perpk = sqrt(perpk)
          max3u = (perpk-projk)*rkk3 + 1
          max3l = (perpk+projk)*rkk3 + 1
          max3l1 = max3l + 1
!
          xxk = xkf - max3l1*rk3x
          yyk = ykf - max3l1*rk3y
          zzk = zkf - max3l1*rk3z
          do kk = -max3l,max3u
            xxk = xxk + rk3x
            yyk = yyk + rk3y
            zzk = zzk + rk3z
            rk2 = xxk*xxk + yyk*yyk + zzk*zzk
!
!  Test for zero wavevector and exclude vectors outside the search
!  region, given by rradmx.
!
            if (rk2.gt.zero.and.rk2.le.rlmx2) then
              nkvec = nkvec + 1
              if (nkvec.gt.maxkvec) then
                maxkvec = nkvec + 100
                call changemaxkvec
              endif
              xrk(nkvec) = xxk
              yrk(nkvec) = yyk
              zrk(nkvec) = zzk
            endif
            xkd = xxk - xkv
            ykd = yyk - ykv
            zkd = zzk - zkv
            rk20 = xkd*xkd + ykd*ykd + zkd*zkd
            if (rk20.gt.zero.and.rk20.le.rlmx2) then
              nkvec0 = nkvec0 + 1
              if (nkvec0.gt.maxkvec) then
                maxkvec = nkvec0 + 100
                call changemaxkvec
              endif
              xrk0(nkvec0) = xkd
              yrk0(nkvec0) = ykd
              zrk0(nkvec0) = zkd
            endif
!
!  End of ii,jj,kk loops over reciprocal lattice vectors
!
          enddo
        endif
      enddo
    enddo
  endif
!
!  Allocate local memory now that maxkvec is known for sure
!
  allocate(ktrm0(maxkvec),stat=status)
  if (status/=0) call outofmemory('recip3Dp','ktrm0')
  allocate(ktrm3(maxkvec),stat=status)
  if (status/=0) call outofmemory('recip3Dp','ktrm3')
  allocate(ktrm6(maxkvec),stat=status)
  if (status/=0) call outofmemory('recip3Dp','ktrm6')
  allocate(ktrm60(maxkvec),stat=status)
  if (status/=0) call outofmemory('recip3Dp','ktrm60')
  allocate(tmp6(6,maxkvec),stat=status)
  if (status/=0) call outofmemory('recip3Dp','tmp6')
  if (lgroupvelocity) then
    allocate(ktrmdk(3_i4,maxkvec),stat=status)
    if (status/=0) call outofmemory('recip3Dp','ktrmdk')
    allocate(ktrm6dk(3_i4,maxkvec),stat=status)
    if (status/=0) call outofmemory('recip3Dp','ktrm6dk')
  endif
!
!  Generate reciprocal space term
!
  if (lc6loc) then
    c6t1 = vol4pi*sqrtpi/48.0_dp
!
!  Reciprocal space self term
!
    do i = 1,nkvec
      rk2 = xrk(i)*xrk(i) + yrk(i)*yrk(i) + zrk(i)*zrk(i)
      arge = - rk2*eta4
      xpon = exp(arge)
      rrk2 = 1.0_dp/rk2
      ktrm(i) = xpon*vol4pi*rrk2*angstoev
      rk = sqrt(rk2)
      rketa2 = 0.5_dp*rk/seta
      c6t2 = sqrtpi*g_derfc(rketa2)
      rketa2 = 1.0_dp/rketa2
      c6t3 = 0.5_dp*rketa2**3-rketa2
      c6t3 = c6t3*xpon
      c6t4 = c6t1*rk2*rk
      ktrm6(i) = c6t4*(c6t2 + c6t3)
      if (lgroupvelocity) then
        dkdk = - 2.0_dp*ktrm(i)*(rrk2 + eta4)
        ktrmdk(1,i) = dkdk*xrk(i)
        ktrmdk(2,i) = dkdk*yrk(i)
        ktrmdk(3,i) = dkdk*zrk(i)
!
        c6t2dk = - xpon/seta
        c6t3dk = - (1.5_dp*rketa2**2 - 1.0_dp)*2.0_dp*seta*rrk2*xpon - 2.0_dp*rk*eta4*c6t3
        c6t4dk = 3.0_dp*c6t1*rk2
        ktrm6dk(1,i) = (c6t4*(c6t2dk + c6t3dk) + c6t4dk*(c6t2 + c6t3))*xrk(i)*sqrt(rrk2)
        ktrm6dk(2,i) = (c6t4*(c6t2dk + c6t3dk) + c6t4dk*(c6t2 + c6t3))*yrk(i)*sqrt(rrk2)
        ktrm6dk(3,i) = (c6t4*(c6t2dk + c6t3dk) + c6t4dk*(c6t2 + c6t3))*zrk(i)*sqrt(rrk2)
      endif
    enddo
    if (.not.lgamma) then
      do i = 1,nkvec0
        rk2 = xrk0(i)*xrk0(i) + yrk0(i)*yrk0(i) + zrk0(i)*zrk0(i)
        arge = - rk2*eta4
        xpon = exp(arge)
        rrk2 = 1.0_dp/rk2
        ktrm0(i) = xpon*vol4pi*rrk2*angstoev
        rk = sqrt(rk2)
        rketa2 = 0.5_dp*rk/seta
        c6t2 = sqrtpi*g_derfc(rketa2)
        rketa2 = 1.0_dp/rketa2
        c6t3 = 0.5_dp*rketa2**3 - rketa2
        c6t3 = c6t3*xpon
        c6t4 = c6t1*rk2*rk
        ktrm60(i) = c6t4*(c6t2 + c6t3)
      enddo
    endif
  else
    do i = 1,nkvec
      rk2 = xrk(i)*xrk(i) + yrk(i)*yrk(i) + zrk(i)*zrk(i)
      arge = - rk2*eta4
      xpon = exp(arge)
      rrk2 = 1.0_dp/rk2
      ktrm(i) = xpon*vol4pi*rrk2*angstoev
      if (lgroupvelocity) then
        dkdk = - 2.0_dp*ktrm(i)*(rrk2 + eta4)
        ktrmdk(1,i) = dkdk*xrk(i)
        ktrmdk(2,i) = dkdk*yrk(i)
        ktrmdk(3,i) = dkdk*zrk(i)
      endif
    enddo
    if (.not.lgamma) then
      do i = 1,nkvec0
        rk2 = xrk0(i)*xrk0(i) + yrk0(i)*yrk0(i) + zrk0(i)*zrk0(i)
        arge = - rk2*eta4
        xpon = exp(arge)
        rrk2 = 1.0_dp/rk2
        ktrm0(i) = xpon*vol4pi*rrk2*angstoev
      enddo
    endif
  endif
!************************************************************************************
!  Algorithm for cases where dispersion cannot be factorised into one centre terms  *
!************************************************************************************
  do iv = 1,nkvec
    tmp6(1,iv) = xrk(iv)*xrk(iv)
    tmp6(2,iv) = yrk(iv)*yrk(iv)
    tmp6(3,iv) = zrk(iv)*zrk(iv)
    tmp6(4,iv) = yrk(iv)*zrk(iv)
    tmp6(5,iv) = xrk(iv)*zrk(iv)
    tmp6(6,iv) = xrk(iv)*yrk(iv)
  enddo
  if (lc6loc) then
!------------------
!  Charge and C6  |
!------------------
    do i = 1,numat
      oci = occuf(i)
      qli = qf(i)*oci
      nati = nat(i)
      ntypi = nftype(i)
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
      ni = 3*(i-1)
      ix = ni + 1
      iy = ix + 1
      iz = iy + 1
      jloop: do j = 1,i
        ocj = occuf(j)
        qlj = qf(j)*ocj
        natj = nat(j)
        ntypj = nftype(j)
!
!  Find C6 term for pair
!
        c6tot = 0.0_dp
        do n = 1,npote
          if (nati.eq.nspec1(n).and.natj.eq.nspec2(n)) then
            if ((ntypi.eq.nptyp1(n).or.nptyp1(n).eq.0).and.(ntypj.eq.nptyp2(n).or.nptyp2(n).eq.0)) then
              if (nptype(n).eq.1.or.nptype(n).eq.7) then
                c6tot = c6tot + twopot(3,n)
              elseif (nptype(n).eq.2) then
                c6tot = c6tot + twopot(2,n)
              elseif (nptype(n).eq.57) then
                c6tot = c6tot + twopot(3,n)*twopot(4,n)*twopot(5,n)**6
              endif
            endif
          elseif (natj.eq.nspec1(n).and.nati.eq.nspec2(n)) then
            if ((ntypj.eq.nptyp1(n).or.nptyp1(n).eq.0).and.(ntypi.eq.nptyp2(n).or.nptyp2(n).eq.0)) then
              if (nptype(n).eq.1.or.nptype(n).eq.7) then
                c6tot = c6tot + twopot(3,n)
              elseif (nptype(n).eq.2) then
                c6tot = c6tot + twopot(2,n)
              elseif (nptype(n).eq.57) then
                c6tot = c6tot + twopot(3,n)*twopot(4,n)*twopot(5,n)**6
              endif
            endif
          endif
        enddo
        ldoc6 = (lc6loc.and.abs(c6tot).gt.thresh_c6)
!
!  Check whether this pair is worth doing
!
        if (abs(qli*qlj).lt.thresh_q.and..not.ldoc6) cycle jloop
!
        nmj = natmol(j)
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
        qfct = qli*qlj
        d2qr(1:6) = 0.0_dp
        d2qi(1:6) = 0.0_dp
        if (lgroupvelocity) then
          cdk(1:3,1:6) = 0.0_dpc
        endif
        if (ldoc6) then
          c6tot = c6tot*oci*ocj
          do iv = 1,nkvec
            xrkk = xrk(iv)
            yrkk = yrk(iv)
            zrkk = zrk(iv)
            if (lphasecom) then
              arg = xrkk*xd - xkv*(xd-xcom) + yrkk*yd - ykv*(yd-ycom) + zrkk*zd - zkv*(zd-zcom)
            else
              arg = (xrkk-xkv)*xd + (yrkk-ykv)*yd + (zrkk-zkv)*zd
            endif
            cosa = cos(arg)
            sina = sin(arg)
            cosq = cosa*ktrm(iv)*qfct
            cos6 = cosa*ktrm6(iv)*c6tot
            sinq = sina*ktrm(iv)*qfct
            sin6 = sina*ktrm6(iv)*c6tot
            d2trmr = cosq - cos6
            d2trmi = sinq - sin6
            do ns = 1,6
              d2qr(ns) = d2qr(ns) + d2trmr*tmp6(ns,iv)
              d2qi(ns) = d2qi(ns) - d2trmi*tmp6(ns,iv)
            enddo
            if (lgroupvelocity) then
              d2trmx = qfct*ktrmdk(1,iv) - c6tot*ktrm6dk(1,iv)
              d2trmy = qfct*ktrmdk(2,iv) - c6tot*ktrm6dk(2,iv)
              d2trmz = qfct*ktrmdk(3,iv) - c6tot*ktrm6dk(3,iv)
              do ns = 1,6
!
!  Terms from derivatives of G terms
!
                cdk(1,ns) = cdk(1,ns) + dcmplx(cosa*d2trmx*tmp6(ns,iv),-sina*d2trmx*tmp6(ns,iv))
                cdk(2,ns) = cdk(2,ns) + dcmplx(cosa*d2trmy*tmp6(ns,iv),-sina*d2trmy*tmp6(ns,iv))
                cdk(3,ns) = cdk(3,ns) + dcmplx(cosa*d2trmz*tmp6(ns,iv),-sina*d2trmz*tmp6(ns,iv))
              enddo
!
!  Terms from derivatives of G_alpha.G_beta terms
!
              cdk(1,1) = cdk(1,1) + dcmplx(2.0_dp*d2trmr*xrk(iv),-2.0_dp*d2trmi*xrk(iv))
              cdk(2,2) = cdk(2,2) + dcmplx(2.0_dp*d2trmr*yrk(iv),-2.0_dp*d2trmi*yrk(iv))
              cdk(3,3) = cdk(3,3) + dcmplx(2.0_dp*d2trmr*zrk(iv),-2.0_dp*d2trmi*zrk(iv))
              cdk(2,4) = cdk(2,4) + dcmplx(d2trmr*zrk(iv),-d2trmi*zrk(iv))
              cdk(3,4) = cdk(3,4) + dcmplx(d2trmr*yrk(iv),-d2trmi*yrk(iv))
              cdk(1,5) = cdk(1,5) + dcmplx(d2trmr*zrk(iv),-d2trmi*zrk(iv))
              cdk(3,5) = cdk(3,5) + dcmplx(d2trmr*xrk(iv),-d2trmi*xrk(iv))
              cdk(1,6) = cdk(1,6) + dcmplx(d2trmr*yrk(iv),-d2trmi*yrk(iv))
              cdk(2,6) = cdk(2,6) + dcmplx(d2trmr*xrk(iv),-d2trmi*xrk(iv))
!
!  Terms from derivative with respect to arg
!
              if (lrigid) then
                trmc = qfct*ktrm(iv) - c6tot*ktrm6(iv)
                do ns = 1,6
                  cdk(1,ns) = cdk(1,ns) + dcmplx(-sina*xcom*trmc*tmp6(ns,iv),-cosa*xcom*trmc*tmp6(ns,iv))
                  cdk(2,ns) = cdk(2,ns) + dcmplx(-sina*ycom*trmc*tmp6(ns,iv),-cosa*ycom*trmc*tmp6(ns,iv))
                  cdk(3,ns) = cdk(3,ns) + dcmplx(-sina*zcom*trmc*tmp6(ns,iv),-cosa*zcom*trmc*tmp6(ns,iv))
                enddo
              endif
            endif
          enddo
          if (i.eq.j) then
!
!  Subtract gamma point terms
!
            if (.not.lgamma) then
              do iv = 1,nkvec0
                xrkk = xrk0(iv)
                yrkk = yrk0(iv)
                zrkk = zrk0(iv)
                arg = xrkk*xd + yrkk*yd + zrkk*zd
                cosa = cos(arg)
                cosq = cosa*ktrm0(iv)*qfct
                cos6 = cosa*ktrm60(iv)*c6tot
                d2trmr = cosq - cos6
                d2qr(1) = d2qr(1) - d2trmr*xrkk*xrkk
                d2qr(2) = d2qr(2) - d2trmr*yrkk*yrkk
                d2qr(3) = d2qr(3) - d2trmr*zrkk*zrkk
                d2qr(4) = d2qr(4) - d2trmr*yrkk*zrkk
                d2qr(5) = d2qr(5) - d2trmr*xrkk*zrkk
                d2qr(6) = d2qr(6) - d2trmr*xrkk*yrkk
              enddo
            else
              d2qr(1:6) = 0.0_dp
            endif
          endif
        else
          do iv = 1,nkvec
            xrkk = xrk(iv)
            yrkk = yrk(iv)
            zrkk = zrk(iv)
            if (lphasecom) then
              arg = xrkk*xd - xkv*(xd-xcom) + yrkk*yd - ykv*(yd-ycom) + zrkk*zd - zkv*(zd-zcom)
            else
              arg = (xrkk-xkv)*xd + (yrkk-ykv)*yd + (zrkk-zkv)*zd
            endif
            cosa = cos(arg)*qfct
            sina = sin(arg)*qfct
            cosq = cosa*ktrm(iv)
            sinq = sina*ktrm(iv)
            do ns = 1,6
              d2qr(ns) = d2qr(ns) + cosq*tmp6(ns,iv)
              d2qi(ns) = d2qi(ns) - sinq*tmp6(ns,iv)
            enddo
            if (lgroupvelocity) then
              do ns = 1,6
!
!  Terms from derivatives of G terms
!
                cdk(1,ns) = cdk(1,ns) + dcmplx(cosa*ktrmdk(1,iv)*tmp6(ns,iv),-sina*ktrmdk(1,iv)*tmp6(ns,iv))
                cdk(2,ns) = cdk(2,ns) + dcmplx(cosa*ktrmdk(2,iv)*tmp6(ns,iv),-sina*ktrmdk(2,iv)*tmp6(ns,iv))
                cdk(3,ns) = cdk(3,ns) + dcmplx(cosa*ktrmdk(3,iv)*tmp6(ns,iv),-sina*ktrmdk(3,iv)*tmp6(ns,iv))
              enddo
!
!  Terms from derivatives of G_alpha.G_beta terms
!
              cdk(1,1) = cdk(1,1) + dcmplx(2.0_dp*cosq*xrk(iv),-2.0_dp*sinq*xrk(iv))
              cdk(2,2) = cdk(2,2) + dcmplx(2.0_dp*cosq*yrk(iv),-2.0_dp*sinq*yrk(iv))
              cdk(3,3) = cdk(3,3) + dcmplx(2.0_dp*cosq*zrk(iv),-2.0_dp*sinq*zrk(iv))
              cdk(2,4) = cdk(2,4) + dcmplx(cosq*zrk(iv),-sinq*zrk(iv))
              cdk(3,4) = cdk(3,4) + dcmplx(cosq*yrk(iv),-sinq*yrk(iv))
              cdk(1,5) = cdk(1,5) + dcmplx(cosq*zrk(iv),-sinq*zrk(iv))
              cdk(3,5) = cdk(3,5) + dcmplx(cosq*xrk(iv),-sinq*xrk(iv))
              cdk(1,6) = cdk(1,6) + dcmplx(cosq*yrk(iv),-sinq*yrk(iv))
              cdk(2,6) = cdk(2,6) + dcmplx(cosq*xrk(iv),-sinq*xrk(iv))
!
!  Terms from derivative with respect to arg
!
              if (lrigid) then
                trmc = ktrm(iv)
                do ns = 1,6
                  cdk(1,ns) = cdk(1,ns) + dcmplx(-sina*xcom*trmc*tmp6(ns,iv),-cosa*xcom*trmc*tmp6(ns,iv))
                  cdk(2,ns) = cdk(2,ns) + dcmplx(-sina*ycom*trmc*tmp6(ns,iv),-cosa*ycom*trmc*tmp6(ns,iv))
                  cdk(3,ns) = cdk(3,ns) + dcmplx(-sina*zcom*trmc*tmp6(ns,iv),-cosa*zcom*trmc*tmp6(ns,iv))
                enddo
              endif
            endif
          enddo
          if (i.eq.j) then
!
!  Subtract gamma point terms
!
            if (.not.lgamma) then
              do iv = 1,nkvec0
                xrkk = xrk0(iv)
                yrkk = yrk0(iv)
                zrkk = zrk0(iv)
                arg = xrkk*xd + yrkk*yd + zrkk*zd
                cosa = cos(arg)*qfct
                cosq = cosa*ktrm0(iv)
                d2qr(1) = d2qr(1) - cosq*xrkk*xrkk
                d2qr(2) = d2qr(2) - cosq*yrkk*yrkk
                d2qr(3) = d2qr(3) - cosq*zrkk*zrkk
                d2qr(4) = d2qr(4) - cosq*yrkk*zrkk
                d2qr(5) = d2qr(5) - cosq*xrkk*zrkk
                d2qr(6) = d2qr(6) - cosq*xrkk*yrkk
              enddo
            else
              d2qr(1:6) = 0.0_dp
            endif
          endif
        endif
        nj = 3*(j-1)
!
        derv2(nj+1,ix) = derv2(nj+1,ix) + d2qr(1)
        derv2(nj+2,ix) = derv2(nj+2,ix) + d2qr(6)
        derv2(nj+3,ix) = derv2(nj+3,ix) + d2qr(5)
        derv2(nj+1,iy) = derv2(nj+1,iy) + d2qr(6)
        derv2(nj+2,iy) = derv2(nj+2,iy) + d2qr(2)
        derv2(nj+3,iy) = derv2(nj+3,iy) + d2qr(4)
        derv2(nj+1,iz) = derv2(nj+1,iz) + d2qr(5)
        derv2(nj+2,iz) = derv2(nj+2,iz) + d2qr(4)
        derv2(nj+3,iz) = derv2(nj+3,iz) + d2qr(3)
!
        dervi(nj+1,ix) = dervi(nj+1,ix) + d2qi(1)
        dervi(nj+2,ix) = dervi(nj+2,ix) + d2qi(6)
        dervi(nj+3,ix) = dervi(nj+3,ix) + d2qi(5)
        dervi(nj+1,iy) = dervi(nj+1,iy) + d2qi(6)
        dervi(nj+2,iy) = dervi(nj+2,iy) + d2qi(2)
        dervi(nj+3,iy) = dervi(nj+3,iy) + d2qi(4)
        dervi(nj+1,iz) = dervi(nj+1,iz) + d2qi(5)
        dervi(nj+2,iz) = dervi(nj+2,iz) + d2qi(4)
        dervi(nj+3,iz) = dervi(nj+3,iz) + d2qi(3)
!
        if (lgroupvelocity) then
!
!  Group velocities
!
          derv2dk(1,nj+1,ix) = derv2dk(1,nj+1,ix) + cdk(1,1)
          derv2dk(2,nj+1,ix) = derv2dk(2,nj+1,ix) + cdk(2,1)
          derv2dk(3,nj+1,ix) = derv2dk(3,nj+1,ix) + cdk(3,1)
!
          derv2dk(1,nj+2,ix) = derv2dk(1,nj+2,ix) + cdk(1,6)
          derv2dk(2,nj+2,ix) = derv2dk(2,nj+2,ix) + cdk(2,6)
          derv2dk(3,nj+2,ix) = derv2dk(3,nj+2,ix) + cdk(3,6)
!
          derv2dk(1,nj+3,ix) = derv2dk(1,nj+3,ix) + cdk(1,5)
          derv2dk(2,nj+3,ix) = derv2dk(2,nj+3,ix) + cdk(2,5)
          derv2dk(3,nj+3,ix) = derv2dk(3,nj+3,ix) + cdk(3,5)
!
          derv2dk(1,nj+1,iy) = derv2dk(1,nj+1,iy) + cdk(1,6)
          derv2dk(2,nj+1,iy) = derv2dk(2,nj+1,iy) + cdk(2,6)
          derv2dk(3,nj+1,iy) = derv2dk(3,nj+1,iy) + cdk(3,6)
!
          derv2dk(1,nj+2,iy) = derv2dk(1,nj+2,iy) + cdk(1,2)
          derv2dk(2,nj+2,iy) = derv2dk(2,nj+2,iy) + cdk(2,2)
          derv2dk(3,nj+2,iy) = derv2dk(3,nj+2,iy) + cdk(3,2)
!
          derv2dk(1,nj+3,iy) = derv2dk(1,nj+3,iy) + cdk(1,4)
          derv2dk(2,nj+3,iy) = derv2dk(2,nj+3,iy) + cdk(2,4)
          derv2dk(3,nj+3,iy) = derv2dk(3,nj+3,iy) + cdk(3,4)
!
          derv2dk(1,nj+1,iz) = derv2dk(1,nj+1,iz) + cdk(1,5)
          derv2dk(2,nj+1,iz) = derv2dk(2,nj+1,iz) + cdk(2,5)
          derv2dk(3,nj+1,iz) = derv2dk(3,nj+1,iz) + cdk(3,5)
!
          derv2dk(1,nj+2,iz) = derv2dk(1,nj+2,iz) + cdk(1,4)
          derv2dk(2,nj+2,iz) = derv2dk(2,nj+2,iz) + cdk(2,4)
          derv2dk(3,nj+2,iz) = derv2dk(3,nj+2,iz) + cdk(3,4)
!
          derv2dk(1,nj+3,iz) = derv2dk(1,nj+3,iz) + cdk(1,3)
          derv2dk(2,nj+3,iz) = derv2dk(2,nj+3,iz) + cdk(2,3)
          derv2dk(3,nj+3,iz) = derv2dk(3,nj+3,iz) + cdk(3,3)
        endif
      enddo jloop
    enddo
  else
!----------------
!  Charge only  |
!----------------
    do ii = 1,ncharge
      i = nchargeptr(ii)
      oci = occuf(i)
      qli = qf(i)*oci
      nati = nat(i)
      ntypi = nftype(i)
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
      ni = 3*(i-1)
      ix = ni + 1
      iy = ix + 1
      iz = iy + 1
      do jj = 1,ii
        j = nchargeptr(jj)
        ocj = occuf(j)
        qlj = qf(j)*ocj
        natj = nat(j)
        ntypj = nftype(j)
!
        nmj = natmol(j)
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
        qfct = qli*qlj
        d2qr(1:6) = 0.0_dp
        d2qi(1:6) = 0.0_dp
        if (lgroupvelocity) then
          cdk(1:3,1:6) = 0.0_dpc
        endif
        do iv = 1,nkvec
          xrkk = xrk(iv)
          yrkk = yrk(iv)
          zrkk = zrk(iv)
          if (lphasecom) then
            arg = xrkk*xd - xkv*(xd-xcom) + yrkk*yd - ykv*(yd-ycom) + zrkk*zd - zkv*(zd-zcom)
          else
            arg = (xrkk-xkv)*xd + (yrkk-ykv)*yd + (zrkk-zkv)*zd
            endif
          cosa = cos(arg)*qfct
          sina = sin(arg)*qfct
          cosq = cosa*ktrm(iv)
          sinq = sina*ktrm(iv)
          do ns = 1,6
            d2qr(ns) = d2qr(ns) + cosq*tmp6(ns,iv)
            d2qi(ns) = d2qi(ns) - sinq*tmp6(ns,iv)
          enddo
          if (lgroupvelocity) then
            do ns = 1,6
!
!  Terms from derivatives of G terms
!
              cdk(1,ns) = cdk(1,ns) + dcmplx(cosa*ktrmdk(1,iv)*tmp6(ns,iv),-sina*ktrmdk(1,iv)*tmp6(ns,iv))
              cdk(2,ns) = cdk(2,ns) + dcmplx(cosa*ktrmdk(2,iv)*tmp6(ns,iv),-sina*ktrmdk(2,iv)*tmp6(ns,iv))
              cdk(3,ns) = cdk(3,ns) + dcmplx(cosa*ktrmdk(3,iv)*tmp6(ns,iv),-sina*ktrmdk(3,iv)*tmp6(ns,iv))
            enddo
!
!  Terms from derivatives of G_alpha.G_beta terms
!
            cdk(1,1) = cdk(1,1) + dcmplx(2.0_dp*cosq*xrk(iv),-2.0_dp*sinq*xrk(iv))
            cdk(2,2) = cdk(2,2) + dcmplx(2.0_dp*cosq*yrk(iv),-2.0_dp*sinq*yrk(iv))
            cdk(3,3) = cdk(3,3) + dcmplx(2.0_dp*cosq*zrk(iv),-2.0_dp*sinq*zrk(iv))
            cdk(2,4) = cdk(2,4) + dcmplx(cosq*zrk(iv),-sinq*zrk(iv))
            cdk(3,4) = cdk(3,4) + dcmplx(cosq*yrk(iv),-sinq*yrk(iv))
            cdk(1,5) = cdk(1,5) + dcmplx(cosq*zrk(iv),-sinq*zrk(iv))
            cdk(3,5) = cdk(3,5) + dcmplx(cosq*xrk(iv),-sinq*xrk(iv))
            cdk(1,6) = cdk(1,6) + dcmplx(cosq*yrk(iv),-sinq*yrk(iv))
            cdk(2,6) = cdk(2,6) + dcmplx(cosq*xrk(iv),-sinq*xrk(iv))
!
!  Terms from derivative with respect to arg
!
            if (lrigid) then
              trmc = ktrm(iv)
              do ns = 1,6
                cdk(1,ns) = cdk(1,ns) + dcmplx(-sina*xcom*trmc*tmp6(ns,iv),-cosa*xcom*trmc*tmp6(ns,iv))
                cdk(2,ns) = cdk(2,ns) + dcmplx(-sina*ycom*trmc*tmp6(ns,iv),-cosa*ycom*trmc*tmp6(ns,iv))
                cdk(3,ns) = cdk(3,ns) + dcmplx(-sina*zcom*trmc*tmp6(ns,iv),-cosa*zcom*trmc*tmp6(ns,iv))
              enddo
            endif
          endif
        enddo
        if (i.eq.j) then
!
!  Subtract gamma point terms
!
          if (.not.lgamma) then
            do iv = 1,nkvec0
              xrkk = xrk0(iv)
              yrkk = yrk0(iv)
              zrkk = zrk0(iv)
              arg = xrkk*xd + yrkk*yd + zrkk*zd
              cosa = cos(arg)*qfct
              cosq = cosa*ktrm0(iv)
              d2qr(1) = d2qr(1) - cosq*xrkk*xrkk
              d2qr(2) = d2qr(2) - cosq*yrkk*yrkk
              d2qr(3) = d2qr(3) - cosq*zrkk*zrkk
              d2qr(4) = d2qr(4) - cosq*yrkk*zrkk
              d2qr(5) = d2qr(5) - cosq*xrkk*zrkk
              d2qr(6) = d2qr(6) - cosq*xrkk*yrkk
            enddo
          else
            d2qr(1:6) = 0.0_dp
          endif
        endif
        nj = 3*(j-1)
!
        derv2(nj+1,ix) = derv2(nj+1,ix) + d2qr(1)
        derv2(nj+2,ix) = derv2(nj+2,ix) + d2qr(6)
        derv2(nj+3,ix) = derv2(nj+3,ix) + d2qr(5)
        derv2(nj+1,iy) = derv2(nj+1,iy) + d2qr(6)
        derv2(nj+2,iy) = derv2(nj+2,iy) + d2qr(2)
        derv2(nj+3,iy) = derv2(nj+3,iy) + d2qr(4)
        derv2(nj+1,iz) = derv2(nj+1,iz) + d2qr(5)
        derv2(nj+2,iz) = derv2(nj+2,iz) + d2qr(4)
        derv2(nj+3,iz) = derv2(nj+3,iz) + d2qr(3)
!
        dervi(nj+1,ix) = dervi(nj+1,ix) + d2qi(1)
        dervi(nj+2,ix) = dervi(nj+2,ix) + d2qi(6)
        dervi(nj+3,ix) = dervi(nj+3,ix) + d2qi(5)
        dervi(nj+1,iy) = dervi(nj+1,iy) + d2qi(6)
        dervi(nj+2,iy) = dervi(nj+2,iy) + d2qi(2)
        dervi(nj+3,iy) = dervi(nj+3,iy) + d2qi(4)
        dervi(nj+1,iz) = dervi(nj+1,iz) + d2qi(5)
        dervi(nj+2,iz) = dervi(nj+2,iz) + d2qi(4)
        dervi(nj+3,iz) = dervi(nj+3,iz) + d2qi(3)
!
        if (lgroupvelocity) then
!
!  Group velocities
!
          derv2dk(1,nj+1,ix) = derv2dk(1,nj+1,ix) + cdk(1,1)
          derv2dk(2,nj+1,ix) = derv2dk(2,nj+1,ix) + cdk(2,1)
          derv2dk(3,nj+1,ix) = derv2dk(3,nj+1,ix) + cdk(3,1)
!
          derv2dk(1,nj+2,ix) = derv2dk(1,nj+2,ix) + cdk(1,6)
          derv2dk(2,nj+2,ix) = derv2dk(2,nj+2,ix) + cdk(2,6)
          derv2dk(3,nj+2,ix) = derv2dk(3,nj+2,ix) + cdk(3,6)
!
          derv2dk(1,nj+3,ix) = derv2dk(1,nj+3,ix) + cdk(1,5)
          derv2dk(2,nj+3,ix) = derv2dk(2,nj+3,ix) + cdk(2,5)
          derv2dk(3,nj+3,ix) = derv2dk(3,nj+3,ix) + cdk(3,5)
!
          derv2dk(1,nj+1,iy) = derv2dk(1,nj+1,iy) + cdk(1,6)
          derv2dk(2,nj+1,iy) = derv2dk(2,nj+1,iy) + cdk(2,6)
          derv2dk(3,nj+1,iy) = derv2dk(3,nj+1,iy) + cdk(3,6)
!
          derv2dk(1,nj+2,iy) = derv2dk(1,nj+2,iy) + cdk(1,2)
          derv2dk(2,nj+2,iy) = derv2dk(2,nj+2,iy) + cdk(2,2)
          derv2dk(3,nj+2,iy) = derv2dk(3,nj+2,iy) + cdk(3,2)
!
          derv2dk(1,nj+3,iy) = derv2dk(1,nj+3,iy) + cdk(1,4)
          derv2dk(2,nj+3,iy) = derv2dk(2,nj+3,iy) + cdk(2,4)
          derv2dk(3,nj+3,iy) = derv2dk(3,nj+3,iy) + cdk(3,4)
!
          derv2dk(1,nj+1,iz) = derv2dk(1,nj+1,iz) + cdk(1,5)
          derv2dk(2,nj+1,iz) = derv2dk(2,nj+1,iz) + cdk(2,5)
          derv2dk(3,nj+1,iz) = derv2dk(3,nj+1,iz) + cdk(3,5)
!
          derv2dk(1,nj+2,iz) = derv2dk(1,nj+2,iz) + cdk(1,4)
          derv2dk(2,nj+2,iz) = derv2dk(2,nj+2,iz) + cdk(2,4)
          derv2dk(3,nj+2,iz) = derv2dk(3,nj+2,iz) + cdk(3,4)
!
          derv2dk(1,nj+3,iz) = derv2dk(1,nj+3,iz) + cdk(1,3)
          derv2dk(2,nj+3,iz) = derv2dk(2,nj+3,iz) + cdk(2,3)
          derv2dk(3,nj+3,iz) = derv2dk(3,nj+3,iz) + cdk(3,3)
        endif
      enddo
    enddo
  endif

  if (lDoQDeriv2) then
!**************************************************************
!  Calculation of charge derivative contribution for EEM/QEq  *
!**************************************************************
!
!  To save space :
!  d1i  is stored in trm3
!  d2i2 is stored in csin
!  d2ij is stored in argc
!  d2j2 is stored in sine
!
    do iv = 1,nkvec
      csin(iv) = 0.0_dp
      sine(iv) = 0.0_dp
    enddo
    ix = - 2
    iy = - 1
    iz =   0
    do i = 1,numat
      oci = occuf(i)
      qli = qf(i)*oci
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
      ix = ix + 3
      iy = iy + 3
      iz = iz + 3
      jx = - 2
      jy = - 1
      jz =   0
      do j = 1,i
        ocj = occuf(j)
        if (i.eq.j) then
          fct = 0.5_dp
        else
          fct = 1.0_dp
        endif
        qlj = qf(j)*ocj
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
        jx = jx + 3
        jy = jy + 3
        jz = jz + 3
!
!  Find relative vector between atoms
!
        xd = xclat(j) - xci
        yd = yclat(j) - yci
        zd = zclat(j) - zci
        do iv = 1,nkvec
          xrkk = xrk(iv)
          yrkk = yrk(iv)
          zrkk = zrk(iv)
          if (lphasecom) then
            arg = xrkk*(xd - xcom) + yrkk*(yd - ycom) + zrkk*(zd - zcom)
          else
            arg = xrkk*xd + yrkk*yd + zrkk*zd
          endif
          cosa = cos(arg)*fct
          sina = sin(arg)*fct
! d2E/dqi.dqj
          argc(iv) = cosa*oci*ocj*ktrm(iv)
! d2E/dqi.dr
          ktrm3(iv) = - ktrm(iv)*sina
! dE/dqi  
          ktrm6(iv) = cosa*oci*qlj*ktrm(iv)
! dE/dqj  
          ktrm60(iv) = cosa*qli*ocj*ktrm(iv)
        enddo
!
!  Call d2charge
!
        d2self = 0.0_dp
        d1ix = 0.0_dp
        d1iy = 0.0_dp
        d1iz = 0.0_dp
        do iv = 1,nkvec
          d1ix = d1ix + ktrm3(iv)*xrk(iv)
          d1iy = d1iy + ktrm3(iv)*yrk(iv)
          d1iz = d1iz + ktrm3(iv)*zrk(iv)
        enddo
        d1jx = d1ix*qli*ocj
        d1jy = d1iy*qli*ocj
        d1jz = d1iz*qli*ocj
        d1ix = d1ix*qlj*oci
        d1iy = d1iy*qlj*oci
        d1iz = d1iz*qlj*oci
        call d2chargep(i,j,nkvec,ix,iy,iz,jx,jy,jz,xkv,ykv,zkv,ktrm6,ktrm60, &
                       d1ix,d1iy,d1iz,d1jx,d1jy,d1jz,csin,argc,sine,d2self, &
                       0.0_dp,0.0_dp,.false.)
      enddo
    enddo
  endif
!
!  Free local memory
!
  if (lgroupvelocity) then
    deallocate(ktrm6dk,stat=status)
    if (status/=0) call deallocate_error('recip3Dp','ktrm6dk')
    deallocate(ktrmdk,stat=status)
    if (status/=0) call deallocate_error('recip3Dp','ktrmdk')
  endif
  deallocate(tmp6,stat=status)
  if (status/=0) call deallocate_error('recip3Dp','tmp6')
  deallocate(ktrm60,stat=status)
  if (status/=0) call deallocate_error('recip3Dp','ktrm60')
  deallocate(ktrm6,stat=status)
  if (status/=0) call deallocate_error('recip3Dp','ktrm6')
  deallocate(ktrm3,stat=status)
  if (status/=0) call deallocate_error('recip3Dp','ktrm3')
  deallocate(ktrm0,stat=status)
  if (status/=0) call deallocate_error('recip3Dp','ktrm0')
!
!  Symmetrise second derivative matrix - used to be done in realp, but this is
!  unsafe if the noreal keyword is specified or if using the fc_supercell option.
!
  do i = 2,3*numat
    do j = 1,i-1
      derv2(i,j) = derv2(j,i)
      dervi(i,j) = - dervi(j,i)
    enddo
  enddo
!
!  Timing
!
  time1 = g_cpu_time()
  tion = tion + time1 - time0
#ifdef TRACE
  call trace_out('recip3Dp')
#endif
!
  return
  end
