  subroutine recip3Dfd(matom,erecip,ec6,lgrad1)
!
!  Calculation of electrostatic potential and derivatives,
!  including strain derivatives. Reciprocal space part.
!  Finite difference version that focuses on derivatives of matom.
!
!  12/17 Created from recip3Dmd
!   2/18 Trace added
!   9/18 Handling of lstraincell algorithm added
!   9/18 Strain module introduced
!  11/18 Finite strain flag introduced instead of lstraincell
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!  12/19 Rigid molecules added
!   1/20 Correction to com coordinate setting
!   3/20 Tolerance for ldoc6 made global
!   3/20 Use of charge pointer added
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
!  Julian Gale, CIC, Curtin University, July 2020
!
  use configurations, only : nregionno
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
  integer(i4), intent(in)                     :: matom
  logical,     intent(in)                     :: lgrad1
  real(dp),    intent(out)                    :: ec6
  real(dp),    intent(out)                    :: erecip
!
!  Local variables
!
  integer(i4)                                 :: i
  integer(i4)                                 :: idk
  integer(i4)                                 :: ii
  integer(i4)                                 :: iv
  integer(i4)                                 :: j
  integer(i4)                                 :: jj
  integer(i4)                                 :: k
  integer(i4)                                 :: kk
  integer(i4)                                 :: kl
  integer(i4)                                 :: n
  integer(i4)                                 :: nati
  integer(i4)                                 :: natj
  integer(i4)                                 :: nmi
  integer(i4)                                 :: nmj
  integer(i4)                                 :: nregioni
  integer(i4)                                 :: nregionj
  integer(i4)                                 :: ntypj
  integer(i4)                                 :: ntypi
  integer(i4)                                 :: status
  logical                                     :: lc6loc
  logical                                     :: ldoc6
  logical                                     :: lsg1
  real(dp)                                    :: arg
  real(dp)                                    :: arge
  real(dp)                                    :: c6self2
  real(dp)                                    :: c6t1
  real(dp)                                    :: c6t2
  real(dp)                                    :: c6t3
  real(dp)                                    :: c6t4
  real(dp)                                    :: c6tot
  real(dp)                                    :: cos6
  real(dp)                                    :: cosa
  real(dp)                                    :: cosq
  real(dp)                                    :: g_cpu_time
  real(dp)                                    :: csin6
  real(dp), dimension(:),   allocatable       :: csink
  real(dp), dimension(:),   allocatable       :: csink6
  real(dp)                                    :: csinq
  real(dp)                                    :: d1trm
  real(dp)                                    :: d3trm
  real(dp)                                    :: dGrds
  real(dp)                                    :: drxyzds(6,3)
  real(dp)                                    :: d2rxyzdsdx(6,3,3)
  real(dp)                                    :: d2rxyzds2(6,6,3)
  real(dp)                                    :: g_derfc
  real(dp)                                    :: esum
  real(dp)                                    :: factor
  real(dp)                                    :: fct
  real(dp), dimension(:),   allocatable       :: ktrm4
  real(dp), dimension(:),   allocatable       :: ktrm6
  real(dp), dimension(:),   allocatable       :: ktrm62
  real(dp)                                    :: kvv(3)
  real(dp)                                    :: oci
  real(dp)                                    :: ocj
  real(dp), dimension(:),   allocatable       :: phsq
  real(dp)                                    :: qfct
  real(dp)                                    :: qli
  real(dp)                                    :: qlj
  real(dp)                                    :: rangstoev
  real(dp)                                    :: reta
  real(dp)                                    :: rk
  real(dp)                                    :: rk2
  real(dp)                                    :: rketa2
  real(dp)                                    :: rrk2
  real(dp)                                    :: sina
  real(dp), dimension(:),   allocatable       :: sinek
  real(dp), dimension(:),   allocatable       :: sinek6
  real(dp)                                    :: sineq
  real(dp)                                    :: sinqx
  real(dp)                                    :: sinqy
  real(dp)                                    :: sinqz
  real(dp)                                    :: strdervloc(6)
  real(dp)                                    :: strm1
  real(dp)                                    :: time0
  real(dp)                                    :: time1
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
  real(dp)                                    :: xrkk
  real(dp)                                    :: yrkk
  real(dp)                                    :: zrkk
  real(dp)                                    :: xpon
#ifdef TRACE
  call trace_in('recip3Dfd')
#endif
!
  time0 = g_cpu_time()
!
!  Initialise energies
!
  ec6 = 0.0_dp
  erecip = 0.0_dp
  lsg1 = (lstr.and.lgrad1)
  lc6loc = (lc6.and.ndim.eq.3)
!
!  Allocate local memory
!
  allocate(csink(nkvec),stat=status)
  if (status/=0) call outofmemory('recip3Dfd','csink')
  allocate(sinek(nkvec),stat=status)
  if (status/=0) call outofmemory('recip3Dfd','sinek')
  allocate(csink6(nkvec),stat=status)
  if (status/=0) call outofmemory('recip3Dfd','csink6')
  allocate(sinek6(nkvec),stat=status)
  if (status/=0) call outofmemory('recip3Dfd','sinek6')
  allocate(phsq(nkvec),stat=status)
  if (status/=0) call outofmemory('recip3Dfd','phsq')
  allocate(ktrm4(nkvec),stat=status)
  if (status/=0) call outofmemory('recip3Dfd','ktrm4')
  allocate(ktrm6(nkvec),stat=status)
  if (status/=0) call outofmemory('recip3Dfd','ktrm6')
  allocate(ktrm62(nkvec),stat=status)
  if (status/=0) call outofmemory('recip3Dfd','ktrm62')
!
!  Setup
!
  eta4 = 0.25_dp/eta
  reta = eta4/eta
  rangstoev = 1.0_dp/angstoev
!*******************************
!  Sum 1/r**6 + coulomb terms  *
!*******************************
  if (lc6loc) then
    c6t1 = vol4pi*rangstoev*sqrtpi/48.0_dp
!
!  Reciprocal space self term
!
    c6self2 = 4.0_dp*c6t1*eta*seta
    if (lra) then
      kvv(1) = kv(1,1)
      kvv(2) = kv(2,2)
      kvv(3) = kv(3,3)
      do iv = 1,nkvec
        idk = indk(iv)
        ii = (idk/maxindk3) - maxindk
        if (ii.eq.0) then
          factor = 1.0_dp
        else
          factor =  2.0_dp
        endif
        idk = idk - (ii+maxindk)*maxindk3
        jj = (idk/maxindk2) - maxindk
        kk = idk - (jj+maxindk)*maxindk2 - maxindk
        xrk(iv) = ii*kvv(1)
        yrk(iv) = jj*kvv(2)
        zrk(iv) = kk*kvv(3)
        rk2 = xrk(iv)*xrk(iv) + yrk(iv)*yrk(iv) + zrk(iv)*zrk(iv)
        arge = - rk2*eta4
        xpon = exp(arge)
        rrk2 = 1.0_dp/rk2
        ktrm(iv) = xpon*vol4pi*factor*rrk2
        rk = sqrt(rk2)
        rketa2 = 0.5_dp*rk/seta
        c6t2 = sqrtpi*g_derfc(rketa2)
        rketa2 = 1.0_dp/rketa2
        c6t3 = 0.5_dp*rketa2**3- rketa2
        c6t3 = c6t3*xpon
        c6t4 = c6t1*rk2*rk*factor
        ktrm6(iv) = c6t4*(c6t2+c6t3)
        if (lsg1) then
          ktrms(iv) = - 2.0_dp*ktrm(iv)*(4.0_dp*eta+rk2)*eta4*rrk2
          ktrm62(iv) = 3.0*ktrm6(iv)*rrk2
          ktrm62(iv) = ktrm62(iv) - c6t4*xpon*(12.0*eta*seta*rrk2*rrk2)/rk
        endif
      enddo
    else
      do iv = 1,nkvec
        idk = indk(iv)
        ii = (idk/maxindk3) - maxindk
        idk = idk - (ii+maxindk)*maxindk3
        jj = (idk/maxindk2) - maxindk
        kk = idk - (jj+maxindk)*maxindk2 - maxindk
        factor = 2.0_dp
        if (ii.eq.0.and.nkangle.eq.1) then
          factor = 1.0_dp
        elseif (jj.eq.0.and.nkangle.eq.2) then
          factor = 1.0_dp
        elseif (kk.eq.0.and.nkangle.eq.3) then
          factor = 1.0_dp
        elseif (nkangle.eq.0) then
          factor = 1.0_dp
        endif
        xrk(iv) = ii*kv(1,1) + jj*kv(1,2) + kk*kv(1,3)
        yrk(iv) = ii*kv(2,1) + jj*kv(2,2) + kk*kv(2,3)
        zrk(iv) = ii*kv(3,1) + jj*kv(3,2) + kk*kv(3,3)
        rk2 = xrk(iv)*xrk(iv) + yrk(iv)*yrk(iv) + zrk(iv)*zrk(iv)
        arge = - rk2*eta4
        xpon = exp(arge)
        rrk2 = 1.0_dp/rk2
        ktrm(iv) = xpon*vol4pi*factor*rrk2
        rk = sqrt(rk2)
        rketa2 = 0.5_dp*rk/seta
        c6t2 = sqrtpi*g_derfc(rketa2)
        rketa2 = 1.0_dp/rketa2
        c6t3 = 0.5_dp*rketa2**3 - rketa2
        c6t3 = c6t3*xpon
        c6t4 = c6t1*rk2*rk*factor
        ktrm6(iv) = c6t4*(c6t2 + c6t3)
        if (lsg1) then
          ktrms(iv) = - 2.0_dp*ktrm(iv)*(4.0_dp*eta+rk2)*eta4*rrk2
          ktrm62(iv) = 3.0_dp*ktrm6(iv)*rrk2
          ktrm62(iv) = ktrm62(iv) - c6t4*xpon*(12.0*eta*seta*rrk2*rrk2)/rk
        endif
      enddo
    endif
  else
!*****************
!  Coulomb only  *
!*****************
    if (lra) then
      kvv(1) = kv(1,1)
      kvv(2) = kv(2,2)
      kvv(3) = kv(3,3)
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
        kk = idk - (jj+maxindk)*maxindk2 - maxindk
        xrk(iv) = ii*kvv(1)
        yrk(iv) = jj*kvv(2)
        zrk(iv) = kk*kvv(3)
        rk2 = xrk(iv)*xrk(iv) + yrk(iv)*yrk(iv) + zrk(iv)*zrk(iv)
        arge = - rk2*eta4
        xpon = exp(arge)
        rrk2 = 1.0_dp/rk2
        ktrm(iv) = xpon*vol4pi*factor*rrk2
        if (lsg1) then
          ktrms(iv) = - 2.0_dp*ktrm(iv)*(4.0_dp*eta+rk2)*eta4*rrk2
        endif
      enddo
    else
      do iv = 1,nkvec
        idk = indk(iv)
        ii = (idk/maxindk3) - maxindk
        idk = idk - (ii+maxindk)*maxindk3
        jj = (idk/maxindk2) - maxindk
        kk = idk - (jj+maxindk)*maxindk2 - maxindk
        factor = 2.0_dp
        if (ii.eq.0.and.nkangle.eq.1) then
          factor = 1.0_dp
        elseif (jj.eq.0.and.nkangle.eq.2) then
          factor = 1.0_dp
        elseif (kk.eq.0.and.nkangle.eq.3) then
          factor = 1.0_dp
        elseif (nkangle.eq.0) then
          factor = 1.0_dp
        endif
        xrk(iv) = ii*kv(1,1) + jj*kv(1,2) + kk*kv(1,3)
        yrk(iv) = ii*kv(2,1) + jj*kv(2,2) + kk*kv(2,3)
        zrk(iv) = ii*kv(3,1) + jj*kv(3,2) + kk*kv(3,3)
        rk2 = xrk(iv)*xrk(iv) + yrk(iv)*yrk(iv) + zrk(iv)*zrk(iv)
        arge = - rk2*eta4
        xpon = exp(arge)
        rrk2 = 1.0_dp/rk2
        ktrm(iv) = xpon*vol4pi*factor*rrk2
        if (lsg1) then
          ktrms(iv) = - 2.0_dp*ktrm(iv)*(4.0_dp*eta + rk2)*eta4*rrk2
        endif
      enddo
    endif
  endif
!
!  End of set-up section
!
  if (lnorecip) goto 999
!
  if (lsg1) then
    call gstrterms(ndim,maxkvec,nkvec,xrk,yrk,zrk,dg2ds,d2g2dx2,d2g2ds2,.false.)
    if (lrigid) then
      call gxyzterms(ndim,maxkvec,nkvec,xrk,yrk,zrk,dgds,d2gds2,.false.)
    endif
  endif
!
!  Only do case where i = matom
!
  i = matom
  oci = occuf(i)*angstoev
  qli = qf(i)*oci
!
!  Check whether we can skip this atom based on the charge
!
  if (abs(qli).lt.thresh_q.and..not.lc6loc) goto 999
!
  nati = nat(i)
  ntypi = nftype(i)
  nmi = natmol(i)
  xci = xclat(i)
  yci = yclat(i)
  zci = zclat(i)
  nregioni = nregionno(nsft+nrelf2a(i))
!
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
  do j = 1,numat
    ocj = occuf(j)
    if (i.eq.j) then
      ocj = 0.5_dp*ocj
    endif
    qlj = qf(j)*ocj
    natj = nat(j)
    ntypj = nftype(j)
    nmj = natmol(j)
    nregionj = nregionno(nsft+nrelf2a(j))
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
    if (lc6loc) then
!
!  Find C6 term for pair
!
      c6tot = 0.0_dp
      do n = 1,npote
        if (nati.eq.nspec1(n).and.natj.eq.nspec2(n)) then
          if ((ntypi.eq.nptyp1(n).or.nptyp1(n).eq.0).and.(ntypj.eq.nptyp2(n).or.nptyp2(n).eq.0)) then
            if (nptype(n).eq.1.or.nptype(n).eq.7) then
              c6tot = c6tot + twopot(3,n)
            elseif (nptype(n).eq.2.or.nptype(n).eq.21) then
              c6tot = c6tot + twopot(2,n)
            elseif (nptype(n).eq.57) then
              c6tot = c6tot + twopot(3,n)*twopot(4,n)*twopot(5,n)**6
            endif
          endif
        elseif (natj.eq.nspec1(n).and.nati.eq.nspec2(n)) then
          if ((ntypj.eq.nptyp1(n).or.nptyp1(n).eq.0).and.(ntypi.eq.nptyp2(n).or.nptyp2(n).eq.0)) then
            if (nptype(n).eq.1.or.nptype(n).eq.7) then
              c6tot = c6tot + twopot(3,n)
            elseif (nptype(n).eq.2.or.nptype(n).eq.21) then
              c6tot = c6tot + twopot(2,n)
            elseif (nptype(n).eq.57) then
              c6tot = c6tot + twopot(3,n)*twopot(4,n)*twopot(5,n)**6
            endif
          endif
        endif
      enddo
      ldoc6 = (abs(c6tot).gt.thresh_c6)
    else
      ldoc6 = .false.
    endif
!
!  Find relative vector between atoms
!
    xd = xclat(j) - xci
    yd = yclat(j) - yci
    zd = zclat(j) - zci
!
    if (lrigid.and.lsg1) then
      call cartstrterm(ndim,xd,yd,zd,xcom,ycom,zcom,drxyzds,d2rxyzdsdx,d2rxyzds2,.false.)
    endif
!
    qfct = qli*qlj
    csinq = 0.0_dp
    if (lgrad1) then
      sinqx = 0.0_dp
      sinqy = 0.0_dp
      sinqz = 0.0_dp
      if (lsg1) then
        strdervloc(1:6) = 0.0_dp
      endif
    endif
    if (ldoc6) then
      c6tot = c6tot*oci*ocj
      if (lgrad1) then
        csin6 = 0.0_dp
        do iv = 1,nkvec
          xrkk = xrk(iv)
          yrkk = yrk(iv)
          zrkk = zrk(iv)
          arg = xrkk*xd + yrkk*yd + zrkk*zd
          cosa = cos(arg)
          sina = sin(arg)
          cosq = cosa*ktrm(iv)*qfct
          cos6 = cosa*ktrm6(iv)*c6tot
          csinq = csinq + cosq
          csin6 = csin6 + cos6
          d1trm = (ktrm(iv)*qfct - ktrm6(iv)*c6tot)*sina
          sinqx = sinqx + d1trm*xrkk
          sinqy = sinqy + d1trm*yrkk
          sinqz = sinqz + d1trm*zrkk
          if (lsg1) then
            strm1 = (ktrms(iv)*qfct - ktrm62(iv)*c6tot)
            d3trm = strm1*sina
            strm1 = strm1*cosa
            if (lrigid) then
              do k = 1,6
                dGrds = xrkk*drxyzds(k,1) + yrkk*drxyzds(k,2) + zrkk*drxyzds(k,3) + &
                        xd*dgds(iv,1,k) + yd*dgds(iv,2,k) + zd*dgds(iv,3,k)
                strdervloc(k) = strdervloc(k) - d1trm*dGrds
              enddo
            endif
            do k = 1,6
              strdervloc(k) = strdervloc(k) + strm1*dg2ds(iv,k)
            enddo
          endif
        enddo
      else
        csin6 = 0.0_dp
        do iv = 1,nkvec
          arg = xrk(iv)*xd + yrk(iv)*yd + zrk(iv)*zd
          cosa = cos(arg)
          csinq = csinq + cosa*ktrm(iv)
          csin6 = csin6 + cosa*ktrm6(iv)
        enddo
        csinq = csinq*qfct
        csin6 = csin6*c6tot
      endif
!
!  Lattice energy
!
      erecip = erecip + csinq
      ec6 = ec6 - (csin6 + c6self2*c6tot)
    else
      if (lgrad1) then
        do iv = 1,nkvec
          xrkk = xrk(iv)
          yrkk = yrk(iv)
          zrkk = zrk(iv)
          arg = xrkk*xd + yrkk*yd + zrkk*zd
          cosa = cos(arg)*qfct
          sina = sin(arg)*qfct
          cosq = cosa*ktrm(iv)
          sineq = sina*ktrm(iv)
          csinq = csinq + cosq
          sinqx = sinqx + sineq*xrkk
          sinqy = sinqy + sineq*yrkk
          sinqz = sinqz + sineq*zrkk
          if (lsg1) then
            strm1 = ktrms(iv)*cosa
            d3trm = ktrms(iv)*sina
            if (lrigid) then
              do k = 1,6
                dGrds = xrkk*drxyzds(k,1) + yrkk*drxyzds(k,2) + zrkk*drxyzds(k,3) + &
                        xd*dgds(iv,1,k) + yd*dgds(iv,2,k) + zd*dgds(iv,3,k)
                strdervloc(k) = strdervloc(k) - sineq*dGrds
              enddo
            endif
            do k = 1,6
              strdervloc(k) = strdervloc(k) + strm1*dg2ds(iv,k)
            enddo
          endif
        enddo
      else
        csinq = 0.0_dp
        do iv = 1,nkvec
          arg = xrk(iv)*xd + yrk(iv)*yd + zrk(iv)*zd
          cosa = cos(arg)
          csinq = csinq + cosa*ktrm(iv)
        enddo
        csinq = csinq*qfct
      endif
!
!  Lattice energy
!
      erecip = erecip + csinq
    endif
    if (lsg1) then
!
!  Strain terms
!
      do kl = 1,nstrains
        strderv(kl) = strderv(kl) + strdervloc(kl)
      enddo
    endif
!
!  Internal derivatives
!
    if (lgrad1.and.i.ne.j) then
      xdrv(i) = xdrv(i) + sinqx
      ydrv(i) = ydrv(i) + sinqy
      zdrv(i) = zdrv(i) + sinqz
      xdrv(j) = xdrv(j) - sinqx
      ydrv(j) = ydrv(j) - sinqy
      zdrv(j) = zdrv(j) - sinqz
    endif
  enddo
!**********************************
!  Bond order charge derivatives  *
!**********************************
  if (lgrad1.and.lDoQDeriv1) then
    oci = occuf(i)
    qli = qf(i)*oci
    fct = angstoev
    xci = xclat(i)
    yci = yclat(i)
    zci = zclat(i)
    do j = 1,numat
      ocj = occuf(j)
      if (i.eq.j) then
        fct = 0.5_dp*fct
      endif
      qlj = qf(j)
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
        arg = xrkk*xd + yrkk*yd + zrkk*zd
        cosa = cos(arg)*fct
        argc(iv) = cosa*oci*ocj*ktrm(iv)*qlj
        phsq(iv) = cosa*oci*ocj*ktrm(iv)*qli
      enddo
      call d1charge(i,j,.true.,.true.,nkvec,argc,phsq)
    enddo
  endif
!****************************************************
!  Complete strain derivatives in reciprocal space  *
!****************************************************
  if (lsg1) then
    if (lc6loc) then
      esum = erecip + ec6
    else
      esum = erecip
    endif
    if (lfinitestrain) then
      do i = 1,6
        strderv(i) = strderv(i) - strainddetds(i)*straindet*esum
      enddo
    else
      strderv(1) = strderv(1) - esum
      strderv(2) = strderv(2) - esum
      strderv(3) = strderv(3) - esum
    endif
  endif
!***************
!  Exit point  *
!***************
999 continue
!
!  Free local memory
!
  deallocate(ktrm62,stat=status)
  if (status/=0) call deallocate_error('recip3Dfd','ktrm62')
  deallocate(ktrm6,stat=status)
  if (status/=0) call deallocate_error('recip3Dfd','ktrm6')
  deallocate(ktrm4,stat=status)
  if (status/=0) call deallocate_error('recip3Dfd','ktrm4')
  deallocate(phsq,stat=status)
  if (status/=0) call deallocate_error('recip3Dfd','phsq')
  deallocate(sinek6,stat=status)
  if (status/=0) call deallocate_error('recip3Dfd','sinek6')
  deallocate(csink6,stat=status)
  if (status/=0) call deallocate_error('recip3Dfd','csink6')
  deallocate(sinek,stat=status)
  if (status/=0) call deallocate_error('recip3Dfd','sinek')
  deallocate(csink,stat=status)
  if (status/=0) call deallocate_error('recip3Dfd','csink')
!
!  Timing
!
  time1 = g_cpu_time()
  tion = tion + time1 - time0
#ifdef TRACE
  call trace_out('recip3Dfd')
#endif
!
  return
  end
