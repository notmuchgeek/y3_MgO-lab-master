  subroutine recip3Dmd(erecip,ec6,lgrad1)
!
!  Calculation of electrostatic potential and derivatives,
!  including strain derivatives. Reciprocal space part.
!  MD version
!
!  Freezing now added.
!
!   9/13 Created from recip3D
!   9/13 veck algorithm made an option only
!  10/13 Hardwired maximum for k vector index replaced with maxindk
!   3/14 derfc changed to g_derfc for benefit of ChemShell
!   2/15 MM3buck added
!   4/15 Name of routine updated in memory calls
!   2/18 Trace added
!   9/18 Handling of lstraincell algorithm added
!   9/18 Strain module introduced
!  11/18 Finite strain flag introduced instead of lstraincell
!   7/19 Global sum of vx, vy, vz moved to polarisation
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!  11/19 Rigid molecules added
!   1/20 Correction to com coordinate setting
!   3/20 Use of charge pointer added
!   3/20 Tolerance for ldoc6 made global
!   4/20 d2xyzdsdc added to cartstrterm arguments
!   4/20 derv3c changes reversed as they are no longer required
!   7/20 Separate routine for sumall with 1 argument added
!   7/20 Corrected for site energies when energy only
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
  use energies,       only : siteenergy, eregion2region
  use kspace
  use m_strain,       only : gstrterms, strainddetds, straindet, cartstrterm, gxyzterms
  use molecule
  use optimisation
  use parallel
  use polarise
  use potentialxyz
  use symmetry
  use thresholds,     only : thresh_c6
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  use two
  implicit none
!
!  Passed variables
!
  logical,  intent(in)                        :: lgrad1
  real(dp), intent(out)                       :: ec6
  real(dp), intent(out)                       :: erecip
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
  integer(i4)                                 :: kvec0
  integer(i4)                                 :: n
  integer(i4)                                 :: nati
  integer(i4)                                 :: natj
  integer(i4)                                 :: nlocalkvec
  integer(i4)                                 :: nloopi
  integer(i4)                                 :: nmi
  integer(i4)                                 :: nmj
  integer(i4)                                 :: nprock
  integer(i4)                                 :: nregioni
  integer(i4)                                 :: nregionj
  integer(i4)                                 :: nremainder
  integer(i4)                                 :: ntypj
  integer(i4)                                 :: ntypi
  integer(i4)                                 :: status
  logical                                     :: lc6loc
  logical                                     :: ldoc6
  logical                                     :: lopi
  logical                                     :: lopj
  logical                                     :: lsg1
  logical                                     :: lveck
  real(dp)                                    :: arg
  real(dp)                                    :: argck
  real(dp)                                    :: arge
  real(dp)                                    :: c6i
  real(dp)                                    :: c6self2
  real(dp)                                    :: c6t1
  real(dp)                                    :: c6t2
  real(dp)                                    :: c6t3
  real(dp)                                    :: c6t4
  real(dp)                                    :: c6tot
  real(dp)                                    :: cos6
  real(dp)                                    :: cosa
  real(dp)                                    :: cosi
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
  real(dp)                                    :: esum6
  real(dp)                                    :: factor
  real(dp)                                    :: fct
  real(dp), dimension(:),   allocatable       :: ktrm4
  real(dp), dimension(:),   allocatable       :: ktrm6
  real(dp), dimension(:),   allocatable       :: ktrm62
  real(dp)                                    :: kvv(3)
  real(dp)                                    :: oci
  real(dp)                                    :: ocj
  real(dp), dimension(:),   allocatable       :: phsq
  real(dp)                                    :: phsqk
  real(dp)                                    :: phsqk6
  real(dp)                                    :: phsqksum
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
  real(dp)                                    :: sini
  real(dp), dimension(:),   allocatable       :: sinek
  real(dp), dimension(:),   allocatable       :: sinek6
  real(dp)                                    :: sineq
  real(dp)                                    :: sinqx
  real(dp)                                    :: sinqy
  real(dp)                                    :: sinqz
  real(dp)                                    :: ssum(6)
  real(dp)                                    :: strdervloc(6)
  real(dp)                                    :: strm1
  real(dp), dimension(:),   allocatable       :: sum
  real(dp)                                    :: time0
  real(dp)                                    :: time1
  real(dp)                                    :: trmk
  real(dp)                                    :: trmks
  real(dp)                                    :: trmk6
  real(dp)                                    :: trmkr(6)
  real(dp)                                    :: trmkr6(6)
  real(dp)                                    :: tsum0
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
  call trace_in('recip3Dmd')
#endif
!
  time0 = g_cpu_time()
!
!  Initialise energies
!
  ec6 = 0.0_dp
  erecip = 0.0_dp
  lveck = .false.
  lsg1 = (lstr.and.lgrad1)
  lc6loc = (lc6.and.ndim.eq.3)
!
!  Keyword to turn on veck 
!
  if (index(keyword,'veck').gt.0) lveck = .true.
!
!  If C6 terms cannot use one centre decomposition then avoid lveck algorithm
!
  if (lc6loc.and..not.lc6one) lveck = .false.
!
!  If site energies are required then don't use lveck
!
  if (lsiteenergy) lveck = .false.
!
  if (lveck) then
!
!  Distribute - use a mix of k parallelism and atom parallelism
!
    nlocalkvec = nkvec
    nprock = nprocs
    kvec0 = procid + 1
  else
!
!  Distribute kvec loops
!
    kvec0 = procid + 1
    nprock = nprocs
    nlocalkvec = (nkvec/nprocs)
    nremainder = nkvec - nlocalkvec*nprocs
    if (procid.lt.nremainder) nlocalkvec = nlocalkvec + 1
  endif
!
!  Allocate local memory
!
  allocate(csink(nlocalkvec),stat=status)
  if (status/=0) call outofmemory('recip3Dmd','csink')
  allocate(sinek(nlocalkvec),stat=status)
  if (status/=0) call outofmemory('recip3Dmd','sinek')
  allocate(csink6(nlocalkvec),stat=status)
  if (status/=0) call outofmemory('recip3Dmd','csink6')
  allocate(sinek6(nlocalkvec),stat=status)
  if (status/=0) call outofmemory('recip3Dmd','sinek6')
  allocate(phsq(nlocalkvec),stat=status)
  if (status/=0) call outofmemory('recip3Dmd','phsq')
  allocate(ktrm4(nlocalkvec),stat=status)
  if (status/=0) call outofmemory('recip3Dmd','ktrm4')
  allocate(ktrm6(nlocalkvec),stat=status)
  if (status/=0) call outofmemory('recip3Dmd','ktrm6')
  allocate(ktrm62(nlocalkvec),stat=status)
  if (status/=0) call outofmemory('recip3Dmd','ktrm62')
  allocate(sum(max(nlocalkvec,nstrains)),stat=status)
  if (status/=0) call outofmemory('recip3Dmd','sum')
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
    if (lveck) then
      c6self2 = 4.0_dp*c6t1*eta*seta
    else
      c6self2 = 4.0_dp*c6t1*eta*seta/dble(nprock)
    endif
    iv = 0
!
!  For veck algorithm we need to zero terms so that sum works
!
    if (lveck) then
      xrk(1:nlocalkvec) = 0.0_dp
      yrk(1:nlocalkvec) = 0.0_dp
      zrk(1:nlocalkvec) = 0.0_dp
      ktrm(1:nlocalkvec) = 0.0_dp
      ktrm6(1:nlocalkvec) = 0.0_dp
    endif
    if (lra) then
      kvv(1) = kv(1,1)
      kvv(2) = kv(2,2)
      kvv(3) = kv(3,3)
      do i = kvec0,nkvec,nprock
        if (lveck) then
          iv = i
        else
          iv = iv + 1
        endif
        idk = indk(i)
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
      do i = kvec0,nkvec,nprock
        if (lveck) then
          iv = i
        else
          iv = iv + 1
        endif
        idk = indk(i)
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
    if (lveck) then
!
!  For veck algorithm we need to globalise some K terms
!
      call sumall(xrk,sum,nlocalkvec,"recip3Dmd","xrk")
      xrk(1:nlocalkvec) = sum(1:nlocalkvec)
      call sumall(yrk,sum,nlocalkvec,"recip3Dmd","yrk")
      yrk(1:nlocalkvec) = sum(1:nlocalkvec)
      call sumall(zrk,sum,nlocalkvec,"recip3Dmd","zrk")
      zrk(1:nlocalkvec) = sum(1:nlocalkvec)
      call sumall(ktrm,sum,nlocalkvec,"recip3Dmd","ktrm")
      ktrm(1:nlocalkvec) = sum(1:nlocalkvec)
      call sumall(ktrm6,sum,nlocalkvec,"recip3Dmd","ktrm6")
      ktrm6(1:nlocalkvec) = sum(1:nlocalkvec)
    endif
  else
!*****************
!  Coulomb only  *
!*****************
    iv = 0
!
!  For veck algorithm we need to zero terms so that sum works
!
    if (lveck) then
      xrk(1:nlocalkvec) = 0.0_dp
      yrk(1:nlocalkvec) = 0.0_dp
      zrk(1:nlocalkvec) = 0.0_dp
      ktrm(1:nlocalkvec) = 0.0_dp
    endif
    if (lra) then
      kvv(1) = kv(1,1)
      kvv(2) = kv(2,2)
      kvv(3) = kv(3,3)
      do i = kvec0,nkvec,nprock
        if (lveck) then
          iv = i
        else
          iv = iv + 1
        endif
        idk = indk(i)
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
      do i = kvec0,nkvec,nprock
        if (lveck) then
          iv = i
        else
          iv = iv + 1
        endif
        idk = indk(i)
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
    if (lveck) then
!
!  For veck algorithm we need to globalise some K terms
!
      call sumall(xrk,sum,nlocalkvec,"recip3Dmd","xrk")
      xrk(1:nlocalkvec) = sum(1:nlocalkvec)
      call sumall(yrk,sum,nlocalkvec,"recip3Dmd","yrk")
      yrk(1:nlocalkvec) = sum(1:nlocalkvec)
      call sumall(zrk,sum,nlocalkvec,"recip3Dmd","zrk")
      zrk(1:nlocalkvec) = sum(1:nlocalkvec)
      call sumall(ktrm,sum,nlocalkvec,"recip3Dmd","ktrm")
      ktrm(1:nlocalkvec) = sum(1:nlocalkvec)
    endif
  endif
!
!  End of set-up section
!
  if (lnorecip) goto 999
!*****************************
!  Vectorise over k vectors  *
!*****************************
  if (lveck) then
    if (lc6loc) then
!------------
!  C6 case  |
!------------
      do iv = 1,nlocalkvec
        csink(iv)  = 0.0_dp
        sinek(iv)  = 0.0_dp
        csink6(iv) = 0.0_dp
        sinek6(iv) = 0.0_dp
      enddo
      c6tot = 0.0_dp
      do ii = procid+1,nchargec6,nprocs
        i = nchargec6ptr(ii)
        qli = qf(i)*occuf(i)
        c6i = c6f(i)*occuf(i)
        c6tot = c6tot + c6i
        do iv = 1,nlocalkvec
          argck = xrk(iv)*xclat(i) + yrk(iv)*yclat(i) + zrk(iv)*zclat(i)
          cosi  = cos(argck)
          sini  = sin(argck)
          csink(iv)  = csink(iv)  + qli*cosi
          sinek(iv)  = sinek(iv)  + qli*sini
          csink6(iv) = csink6(iv) + c6i*cosi
          sinek6(iv) = sinek6(iv) + c6i*sini
        enddo
      enddo
!
!  Global sum over kvector contributions
!
      call sumall(csink,sum,nlocalkvec,"recip3Dmd","csink")
      csink(1:nlocalkvec) = sum(1:nlocalkvec)
      call sumall(sinek,sum,nlocalkvec,"recip3Dmd","sinek")
      sinek(1:nlocalkvec) = sum(1:nlocalkvec)
      call sumall(csink6,sum,nlocalkvec,"recip3Dmd","csink6")
      csink6(1:nlocalkvec) = sum(1:nlocalkvec)
      call sumall(sinek6,sum,nlocalkvec,"recip3Dmd","sinek6")
      sinek6(1:nlocalkvec) = sum(1:nlocalkvec)
      call sumone(c6tot,esum6,"recip3Dmd","c6tot")
      c6tot = esum6
!
      esum  = 0.0_dp
      esum6 = 0.0_dp
!***********************
!  Strain derivatives  *
!***********************
      if (lsg1) then
        ssum(1:6) = 0.0_dp
!
        call gstrterms(ndim,maxkvec,nlocalkvec,xrk,yrk,zrk,dg2ds,d2g2dx2,d2g2ds2,.false.)
        if (lrigid) then
          call gxyzterms(ndim,maxkvec,nlocalkvec,xrk,yrk,zrk,dgds,d2gds2,.false.)
        endif
!
        do iv = procid+1,nlocalkvec,nprocs
          phsqk  = (csink(iv)**2 + sinek(iv)**2)
          phsqk6 = (csink6(iv)**2 + sinek6(iv)**2)
          trmk  = ktrm(iv)
          trmk6 = ktrm6(iv)
          esum  = esum  + trmk*phsqk
          esum6 = esum6 + trmk6*phsqk6
          trmks = ktrms(iv)*phsqk - ktrm62(iv)*phsqk6
!
          if (lrigid) then
!
!  Rigid molecules
!
            xrkk = xrk(iv)
            yrkk = yrk(iv)
            zrkk = zrk(iv)
!
            trmkr(1:6) = 0.0_dp
            trmkr6(1:6) = 0.0_dp
!
            do i = 1,numat
              nmi = natmol(i)
              if (nmi.gt.0) then
                argck = xrkk*xclat(i) + yrkk*yclat(i) + zrkk*zclat(i)
                cosi  = cos(argck)
                sini  = sin(argck)
                qli = qf(i)*occuf(i)
                c6i = c6f(i)*occuf(i)
                xcomi = molxyz(1,natinmol(i),nmi)
                ycomi = molxyz(2,natinmol(i),nmi)
                zcomi = molxyz(3,natinmol(i),nmi)
                call cartstrterm(ndim,xclat(i),yclat(i),zclat(i),xcomi,ycomi,zcomi,drxyzds,d2rxyzdsdx,d2rxyzds2,.false.)
!
                do j = 1,6
                  dGrds = xrkk*drxyzds(j,1) + yrkk*drxyzds(j,2) + zrkk*drxyzds(j,3) + &
                          xclat(i)*dgds(iv,1,j) + yclat(i)*dgds(iv,2,j) + zclat(i)*dgds(iv,3,j)
                  trmkr(j) = trmkr(j) - (csink(iv)*sini - sinek(iv)*cosi)*qli*dGrds
                  trmkr6(j) = trmkr6(j) - (csink6(iv)*sini - sinek6(iv)*cosi)*c6i*dGrds
                enddo
              endif
            enddo
            do j = 1,6
              ssum(j) = ssum(j) + trmks*dg2ds(iv,j) + 2.0_dp*(trmk*trmkr(j) - trmk6*trmkr6(j))
            enddo
          else
            do j = 1,6
              ssum(j) = ssum(j) + trmks*dg2ds(iv,j)
            enddo
          endif
        enddo
        do j = 1,6
          strderv(j) = strderv(j) + 0.5_dp*angstoev*ssum(j)
        enddo
      else
        do iv = procid+1,nlocalkvec,nprocs
          phsqk  = (csink(iv)**2 + sinek(iv)**2)
          phsqk6 = (csink6(iv)**2 + sinek6(iv)**2)
          esum  = esum  + ktrm(iv)*phsqk
          esum6 = esum6 + ktrm6(iv)*phsqk6
        enddo
      endif
!*******************
!  Lattice energy  *
!*******************
      erecip = erecip + 0.5_dp*angstoev*esum
      ec6 = ec6 - 0.5_dp*angstoev*esum6
!
!  Only add self term on one node
!
      if (procid.eq.0) ec6 = ec6 - 0.5_dp*angstoev*c6self2*c6tot*c6tot
!************************
!  Internal derivatives *
!************************
      if (lgrad1) then
!
!  First derivatives
!
        do ii = procid+1,nchargec6,nprocs
          i = nchargec6ptr(ii)
          oci = occuf(i)
          qli = qf(i)*oci
          c6i = c6f(i)*oci
          nregioni = nregionno(nsft+nrelf2a(i))
          do iv = 1,nlocalkvec
            argck = xrk(iv)*xclat(i) + yrk(iv)*yclat(i) + zrk(iv)*zclat(i)
            cosi = cos(argck)
            sini = sin(argck)
            phsqk  = qli*(cosi*sinek(iv)  - sini*csink(iv))
            phsqk6 = c6i*(cosi*sinek6(iv) - sini*csink6(iv))
            phsqksum = angstoev*(ktrm(iv)*phsqk - ktrm6(iv)*phsqk6)
!
            esum = angstoev*(qli*ktrm(iv)*(cosi*csink(iv) + sini*sinek(iv)) - &
                             c6i*ktrm6(iv)*(cosi*csink6(iv) + sini*sinek6(iv)))
            eregion2region(nregioni,nregioni) = eregion2region(nregioni,nregioni) + esum
!
            siteenergy(i) = siteenergy(i) + 0.5_dp*esum
!
            xdrv(i) = xdrv(i) + phsqksum*xrk(iv)
            ydrv(i) = ydrv(i) + phsqksum*yrk(iv)
            zdrv(i) = zdrv(i) + phsqksum*zrk(iv)
!
            xregdrv(nregioni) = xregdrv(nregioni) + phsqksum*xrk(iv)
            yregdrv(nregioni) = yregdrv(nregioni) + phsqksum*yrk(iv)
            zregdrv(nregioni) = zregdrv(nregioni) + phsqksum*zrk(iv)
          enddo
!
!  Self term
!
          esum = angstoev*c6i*c6tot*c6self2
          eregion2region(nregioni,nregioni) = eregion2region(nregioni,nregioni) - esum
          siteenergy(i) = siteenergy(i) - 0.5_dp*esum
        enddo
      elseif (lsiteenergy) then
        do ii = procid+1,nchargec6,nprocs
          i = nchargec6ptr(ii)
          oci = occuf(i)
          qli = qf(i)*oci
          c6i = c6f(i)*oci
          nregioni = nregionno(nsft+nrelf2a(i))
          do iv = 1,nlocalkvec
            argck = xrk(iv)*xclat(i) + yrk(iv)*yclat(i) + zrk(iv)*zclat(i)
            cosi = cos(argck)
            sini = sin(argck)
!
            esum = angstoev*(qli*ktrm(iv)*(cosi*csink(iv) + sini*sinek(iv)) - &
                             c6i*ktrm6(iv)*(cosi*csink6(iv) + sini*sinek6(iv)))
            eregion2region(nregioni,nregioni) = eregion2region(nregioni,nregioni) + esum
!
            siteenergy(i) = siteenergy(i) + 0.5_dp*esum
          enddo
!
!  Self term
!
          esum = angstoev*c6i*c6tot*c6self2
          eregion2region(nregioni,nregioni) = eregion2region(nregioni,nregioni) - esum
          siteenergy(i) = siteenergy(i) - 0.5_dp*esum
        enddo
      endif
    else
!---------------
!  Non-C6 case |
!---------------
      do iv = 1,nlocalkvec
        csink(iv) = 0.0_dp
        sinek(iv) = 0.0_dp
      enddo
      do ii = procid+1,ncharge,nprocs
        i = nchargeptr(ii)
        qli = qf(i)*occuf(i)
        do iv = 1,nlocalkvec
          argck = xrk(iv)*xclat(i) + yrk(iv)*yclat(i) + zrk(iv)*zclat(i)
          csink(iv) = csink(iv) + qli*cos(argck)
          sinek(iv) = sinek(iv) + qli*sin(argck)
        enddo
      enddo
!
!  Global sum over kvector contributions
!
      call sumall(csink,sum,nlocalkvec,"recip3Dmd","csink")
      csink(1:nlocalkvec) = sum(1:nlocalkvec)
      call sumall(sinek,sum,nlocalkvec,"recip3Dmd","sinek")
      sinek(1:nlocalkvec) = sum(1:nlocalkvec)
!
      esum = 0.0_dp
!***********************
!  Strain derivatives  *
!***********************
      if (lsg1) then
        ssum(1:6) = 0.0_dp
!
        call gstrterms(ndim,maxkvec,nlocalkvec,xrk,yrk,zrk,dg2ds,d2g2dx2,d2g2ds2,.false.)
        if (lrigid) then
          call gxyzterms(ndim,maxkvec,nlocalkvec,xrk,yrk,zrk,dgds,d2gds2,.false.)
        endif
!
        do iv = procid+1,nlocalkvec,nprocs
          phsqk = (csink(iv)**2 + sinek(iv)**2)
          trmk  = ktrm(iv)
          esum = esum + trmk*phsqk
          trmks = ktrms(iv)*phsqk
          if (lrigid) then
!
!  Rigid molecules
!
            xrkk = xrk(iv)
            yrkk = yrk(iv)
            zrkk = zrk(iv)
!
            trmkr(1:6) = 0.0_dp
!
            do i = 1,numat
              nmi = natmol(i)
              if (nmi.gt.0) then
                argck = xrkk*xclat(i) + yrkk*yclat(i) + zrkk*zclat(i)
                cosi  = cos(argck)
                sini  = sin(argck)
                qli = qf(i)*occuf(i)
                xcomi = molxyz(1,natinmol(i),nmi)
                ycomi = molxyz(2,natinmol(i),nmi)
                zcomi = molxyz(3,natinmol(i),nmi)
                call cartstrterm(ndim,xclat(i),yclat(i),zclat(i),xcomi,ycomi,zcomi,drxyzds,d2rxyzdsdx,d2rxyzds2,.false.)
!
                do j = 1,6
                  dGrds = xrkk*drxyzds(j,1) + yrkk*drxyzds(j,2) + zrkk*drxyzds(j,3) + &
                          xclat(i)*dgds(iv,1,j) + yclat(i)*dgds(iv,2,j) + zclat(i)*dgds(iv,3,j)
                  trmkr(j) = trmkr(j) - (csink(iv)*sini - sinek(iv)*cosi)*qli*dGrds
                enddo
              endif
            enddo
            do j = 1,6
              ssum(j) = ssum(j) + trmks*dg2ds(iv,j) + 2.0_dp*trmk*trmkr(j)
            enddo
          else
            do j = 1,6
              ssum(j) = ssum(j) + trmks*dg2ds(iv,j)
            enddo
          endif
        enddo
        do j = 1,6
          strderv(j) = strderv(j) + 0.5_dp*angstoev*ssum(j)
        enddo
      else
        do iv = procid+1,nlocalkvec,nprocs
          phsqk = (csink(iv)**2 + sinek(iv)**2)
          esum = esum + ktrm(iv)*phsqk
        enddo
      endif
!*******************
!  Lattice energy  *
!*******************
      erecip = erecip + 0.5_dp*angstoev*esum
!************************
!  Internal derivatives *
!************************
      if (lgrad1) then
!
!  First derivatives
!
        do ii = procid+1,ncharge,nprocs
          i = nchargeptr(ii)
          oci = occuf(i)
          qli = qf(i)*oci
          nregioni = nregionno(nsft+nrelf2a(i))
          do iv = 1,nlocalkvec
            argck = xrk(iv)*xclat(i) + yrk(iv)*yclat(i) + zrk(iv)*zclat(i)
            cosi = cos(argck)
            sini = sin(argck)
            phsqk = qli*(cosi*sinek(iv) - sini*csink(iv))
            phsqksum = angstoev*ktrm(iv)*phsqk
!
            esum = qli*angstoev*ktrm(iv)*(cosi*csink(iv) + sini*sinek(iv))
            eregion2region(nregioni,nregioni) = eregion2region(nregioni,nregioni) + esum
!
            siteenergy(i) = siteenergy(i) + 0.5_dp*esum
!
            xdrv(i) = xdrv(i) + phsqksum*xrk(iv)
            ydrv(i) = ydrv(i) + phsqksum*yrk(iv)
            zdrv(i) = zdrv(i) + phsqksum*zrk(iv)
!
            xregdrv(nregioni) = xregdrv(nregioni) + phsqksum*xrk(iv)
            yregdrv(nregioni) = yregdrv(nregioni) + phsqksum*yrk(iv)
            zregdrv(nregioni) = zregdrv(nregioni) + phsqksum*zrk(iv)
          enddo
        enddo
      elseif (lsiteenergy) then
        do ii = procid+1,ncharge,nprocs
          i = nchargeptr(ii)
          oci = occuf(i)
          qli = qf(i)*oci
          nregioni = nregionno(nsft+nrelf2a(i))
          do iv = 1,nlocalkvec
            argck = xrk(iv)*xclat(i) + yrk(iv)*yclat(i) + zrk(iv)*zclat(i)
            cosi = cos(argck)
            sini = sin(argck)
!
            esum = qli*angstoev*ktrm(iv)*(cosi*csink(iv) + sini*sinek(iv))
            eregion2region(nregioni,nregioni) = eregion2region(nregioni,nregioni) + esum
!
            siteenergy(i) = siteenergy(i) + 0.5_dp*esum
          enddo
        enddo
      endif
    endif
  elseif ((lc6loc.and..not.lc6one).or.(lstr.and.latomicstress)) then
!************************************************************************************
!  Algorithm for cases where dispersion cannot be factorised into one centre terms  *
!  and/or atomic stresses are to be calculated.                                     *
!************************************************************************************
    if (lsg1) then
      call gstrterms(ndim,maxkvec,nlocalkvec,xrk,yrk,zrk,dg2ds,d2g2dx2,d2g2ds2,.false.)
      if (lrigid) then
        call gxyzterms(ndim,maxkvec,nlocalkvec,xrk,yrk,zrk,dgds,d2gds2,.false.)
      endif
    endif
    if (lc6loc) then
      nloopi = numat
    else
      nloopi = ncharge
    endif
    do ii = 1,nloopi
      if (lc6loc) then
        i = ii
      else
        i = nchargeptr(ii)
      endif
      oci = occuf(i)*angstoev
      qli = qf(i)*oci
      nati = nat(i)
      ntypi = nftype(i)
      nmi = natmol(i)
      xci = xclat(i)
      yci = yclat(i)
      zci = zclat(i)
      nregioni = nregionno(nsft+nrelf2a(i))
      lopi = (.not.lfreeze.or.lopf(nrelf2a(i)))
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
      do jj = 1,ii
        if (lc6loc) then
          j = jj
        else
          j = nchargeptr(jj)
        endif
        ocj = occuf(j)
        if (i.eq.j) then
          ocj = 0.5_dp*ocj
        endif
        qlj = qf(j)*ocj
        natj = nat(j)
        ntypj = nftype(j)
        nmj = natmol(j)
        nregionj = nregionno(nsft+nrelf2a(j))
        lopj = (.not.lfreeze.or.lopf(nrelf2a(j)))
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
            do iv = 1,nlocalkvec
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
                    strdervloc(k) = strdervloc(k) - sina*ktrm(iv)*qfct*dGrds + sina*ktrm6(iv)*c6tot*dGrds
                  enddo
                endif
                do k = 1,6
                  strdervloc(k) = strdervloc(k) + strm1*dg2ds(iv,k)
                enddo
              endif
            enddo
          else
            csin6 = 0.0_dp
            do iv = 1,nlocalkvec
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
!
          esum = csinq - csin6 - c6self2*c6tot
          eregion2region(nregionj,nregioni) = eregion2region(nregionj,nregioni) + esum
!
          siteenergy(i) = siteenergy(i) + 0.5_dp*esum
          siteenergy(j) = siteenergy(j) + 0.5_dp*esum
        else
          if (lgrad1) then
            do iv = 1,nlocalkvec
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
            do iv = 1,nlocalkvec
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
!
          eregion2region(nregionj,nregioni) = eregion2region(nregionj,nregioni) + csinq
!
          esum = csinq
          siteenergy(i) = siteenergy(i) + 0.5_dp*esum
          siteenergy(j) = siteenergy(j) + 0.5_dp*esum
        endif
        if (lsg1) then
!
!  Strain terms
!
          do kl = 1,nstrains
            strderv(kl) = strderv(kl) + strdervloc(kl)
          enddo
          if (latomicstress) then
            do kl = 1,nstrains
              atomicstress(kl,i) = atomicstress(kl,i) + 0.5_dp*strdervloc(kl)
              atomicstress(kl,j) = atomicstress(kl,j) + 0.5_dp*strdervloc(kl)
            enddo
            do kl = 1,3
              atomicstress(kl,i) = atomicstress(kl,i) - 0.5_dp*esum
              atomicstress(kl,j) = atomicstress(kl,j) - 0.5_dp*esum
            enddo
          endif
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
          if (nregioni.ne.nregionj) then
            xregdrv(nregioni) = xregdrv(nregioni) + sinqx
            yregdrv(nregioni) = yregdrv(nregioni) + sinqy
            zregdrv(nregioni) = zregdrv(nregioni) + sinqz
            xregdrv(nregionj) = xregdrv(nregionj) - sinqx
            yregdrv(nregionj) = yregdrv(nregionj) - sinqy
            zregdrv(nregionj) = zregdrv(nregionj) - sinqz
          endif
        endif
      enddo
    enddo
  else
!***********************************
!  Vectorise over number of atoms  *
!***********************************
    if (lsg1) then
      call gstrterms(ndim,maxkvec,nlocalkvec,xrk,yrk,zrk,dg2ds,d2g2dx2,d2g2ds2,.false.)
      if (lrigid) then
        call gxyzterms(ndim,maxkvec,nlocalkvec,xrk,yrk,zrk,dgds,d2gds2,.false.)
      endif
    endif
    if (lc6loc.and.lc6one) then
!------------
!  C6 case  |
!------------
      do iv = 1,nlocalkvec
        xrkk = xrk(iv)
        yrkk = yrk(iv)
        zrkk = zrk(iv)
        trmk = ktrm(iv)*angstoev
        trmk6 = ktrm6(iv)*angstoev
        csink(iv) = 0.0_dp
        sinek(iv) = 0.0_dp
        csink6(iv) = 0.0_dp
        sinek6(iv) = 0.0_dp
        do ii = 1,nchargec6
          i = nchargec6ptr(ii)
          oci = occuf(i)
          qli = qf(i)*oci
          c6i = c6f(i)*oci
          argc(i) = xrkk*xclat(i) + yrkk*yclat(i) + zrkk*zclat(i)
          csin(i) = cos(argc(i))
          sine(i) = sin(argc(i))
          csink6(iv) = csink6(iv) + csin(i)*c6i
          sinek6(iv) = sinek6(iv) + sine(i)*c6i
          csink(iv) = csink(iv) + csin(i)*qli
          sinek(iv) = sinek(iv) + sine(i)*qli
        enddo
        phsqk = (csink(iv)**2 + sinek(iv)**2)
        phsqk6 = (csink6(iv)**2 + sinek6(iv)**2)
!**********************
!  Lattice energy     *
!**********************
!  erecip has a factor of a half but this multiplied later
        erecip = erecip + trmk*phsqk
        ec6 = ec6 - 0.5_dp*trmk6*phsqk6
!**********************
!  Strain derivatives *
!**********************
        if (lsg1) then
          trmks = 0.5_dp*angstoev*(ktrms(iv)*phsqk - ktrm62(iv)*phsqk6)
!
          if (lrigid) then
!
!  Rigid molecules
!
            trmkr(1:6) = 0.0_dp
            trmkr6(1:6) = 0.0_dp
            do i = 1,numat
              nmi = natmol(i)
              if (nmi.gt.0) then
                qli = qf(i)*occuf(i)
                c6i = c6f(i)*occuf(i)
                xcomi = molxyz(1,natinmol(i),nmi)
                ycomi = molxyz(2,natinmol(i),nmi)
                zcomi = molxyz(3,natinmol(i),nmi)
                call cartstrterm(ndim,xclat(i),yclat(i),zclat(i),xcomi,ycomi,zcomi,drxyzds,d2rxyzdsdx,d2rxyzds2,.false.)
!
                do j = 1,6
                  dGrds = xrkk*drxyzds(j,1) + yrkk*drxyzds(j,2) + zrkk*drxyzds(j,3) + &
                          xclat(i)*dgds(iv,1,j) + yclat(i)*dgds(iv,2,j) + zclat(i)*dgds(iv,3,j)
                  trmkr(j) = trmkr(j) - (csink(iv)*sine(i) - sinek(iv)*csin(i))*qli*dGrds
                  trmkr6(j) = trmkr6(j) - (csink6(iv)*sine(i) - sinek6(iv)*csin(i))*c6i*dGrds
                enddo
              endif
            enddo
!
            do i = 1,6
              strderv(i) = strderv(i) + trmks*dg2ds(iv,i) + trmk*trmkr(i) - trmk6*trmkr6(i)
            enddo
          else
            do i = 1,6
              strderv(i) = strderv(i) + trmks*dg2ds(iv,i)
            enddo
          endif
        endif
!************************
!  Internal derivatives *
!************************
        if (lgrad1) then
          do ii = 1,nchargec6
            i = nchargec6ptr(ii)
            oci = occuf(i)
            qli = qf(i)*oci
            nregioni = nregionno(nsft+nrelf2a(i))
            c6i = c6f(i)*oci
!
            phsqk = qli*(csin(i)*sinek(iv) - sine(i)*csink(iv))
            phsqk6 = c6i*(csin(i)*sinek6(iv) - sine(i)*csink6(iv))
            phsqksum = phsqk*trmk - phsqk6*trmk6
!
            esum = (trmk*qli*(csin(i)*csink(iv) + sine(i)*sinek(iv)) - &
                    trmk6*c6i*(csin(i)*csink6(iv) + sine(i)*sinek6(iv)))
            eregion2region(nregioni,nregioni) = eregion2region(nregioni,nregioni) + esum
!
            siteenergy(i) = siteenergy(i) + 0.5_dp*esum
            xdrv(i) = xdrv(i) + phsqksum*xrkk
            ydrv(i) = ydrv(i) + phsqksum*yrkk
            zdrv(i) = zdrv(i) + phsqksum*zrkk
!
            xregdrv(nregioni) = xregdrv(nregioni) + phsqksum*xrkk
            yregdrv(nregioni) = yregdrv(nregioni) + phsqksum*yrkk
            zregdrv(nregioni) = zregdrv(nregioni) + phsqksum*zrkk
          enddo
        elseif (lsiteenergy) then
          do ii = 1,nchargec6
            i = nchargec6ptr(ii)
            oci = occuf(i)
            nregioni = nregionno(nsft+nrelf2a(i))
            qli = qf(i)*oci
            c6i = c6f(i)*oci
!
            esum = (trmk*qli*(csin(i)*csink(iv) + sine(i)*sinek(iv)) - &
                    trmk6*c6i*(csin(i)*csink6(iv) + sine(i)*sinek6(iv)))
            eregion2region(nregioni,nregioni) = eregion2region(nregioni,nregioni) + esum
            siteenergy(i) = siteenergy(i) + 0.5_dp*esum
          enddo
        endif
!
!  End of loop over k vectors
!
      enddo
!
!  Self term
!
      c6tot = 0.0_dp
      do ii = 1,nchargec6
        i = nchargec6ptr(ii)
        c6tot = c6tot + c6f(i)*occuf(i)
      enddo
!
!  Only add self term on one node
!
      ec6 = ec6 - 0.5_dp*c6self2*c6tot*c6tot*angstoev
      if (lgrad1) then
        do ii = 1,nchargec6
          i = nchargec6ptr(ii)
          oci = occuf(i)
          nregioni = nregionno(nsft+nrelf2a(i))
          c6i = c6f(i)*oci
          esum = angstoev*c6i*c6tot*c6self2
          eregion2region(nregioni,nregioni) = eregion2region(nregioni,nregioni) - esum
          siteenergy(i) = siteenergy(i) - 0.5_dp*esum
        enddo
      endif
    else
!----------
!  No C6  |
!----------
      do iv = 1,nlocalkvec
        csink(iv) = 0.0_dp
        sinek(iv) = 0.0_dp
        xrkk = xrk(iv)
        yrkk = yrk(iv)
        zrkk = zrk(iv)
        trmk = ktrm(iv)*angstoev
        do ii = 1,ncharge
          i = nchargeptr(ii)
          qli = qf(i)*occuf(i)
          argc(i) = xrkk*xclat(i) + yrkk*yclat(i) + zrkk*zclat(i)
          csin(i) = qli*cos(argc(i))
          sine(i) = qli*sin(argc(i))
          csink(iv) = csink(iv) + csin(i)
          sinek(iv) = sinek(iv) + sine(i)
        enddo
        phsqk = (csink(iv)**2 + sinek(iv)**2)
!**********************
!  Lattice energy     *
!**********************
!  erecip has a factor of a half but this multiplied later
        erecip = erecip + trmk*phsqk
!**********************
!  Strain derivatives *
!**********************
        if (lsg1) then
          trmks = 0.5_dp*angstoev*ktrms(iv)*phsqk
          if (lrigid) then
!
!  Rigid molecules
!
            trmkr(1:6) = 0.0_dp
            do i = 1,numat
              nmi = natmol(i)
              if (nmi.gt.0) then
                xcomi = molxyz(1,natinmol(i),nmi)
                ycomi = molxyz(2,natinmol(i),nmi)
                zcomi = molxyz(3,natinmol(i),nmi)
                call cartstrterm(ndim,xclat(i),yclat(i),zclat(i),xcomi,ycomi,zcomi,drxyzds,d2rxyzdsdx,d2rxyzds2,.false.)
!
                do j = 1,nstrains
                  dGrds = xrkk*drxyzds(j,1) + yrkk*drxyzds(j,2) + zrkk*drxyzds(j,3) + &
                          xclat(i)*dgds(iv,1,j) + yclat(i)*dgds(iv,2,j) + zclat(i)*dgds(iv,3,j)
                  trmkr(j) = trmkr(j) - csink(iv)*sine(i)*dGrds + sinek(iv)*csin(i)*dGrds
                enddo
              endif
            enddo
!
            do i = 1,nstrains
              strderv(i) = strderv(i) + trmks*dg2ds(iv,i) + trmk*trmkr(i)
            enddo
          else
            do i = 1,nstrains
              strderv(i) = strderv(i) + trmks*dg2ds(iv,i)
            enddo
          endif
        endif
!************************
!  Internal derivatives *
!************************
        if (lgrad1) then
          do ii = 1,ncharge
            i = nchargeptr(ii)
            oci = occuf(i)
            nregioni = nregionno(nsft+nrelf2a(i))
!
            phsqk = (csin(i)*sinek(iv) - sine(i)*csink(iv))
            phsqksum = trmk*phsqk
!
            esum = trmk*(csin(i)*csink(iv) + sine(i)*sinek(iv))
            eregion2region(nregioni,nregioni) = eregion2region(nregioni,nregioni) + esum
!
            siteenergy(i) = siteenergy(i) + 0.5_dp*esum
!
            xdrv(i) = xdrv(i) + phsqksum*xrkk
            ydrv(i) = ydrv(i) + phsqksum*yrkk
            zdrv(i) = zdrv(i) + phsqksum*zrkk
!
            xregdrv(nregioni) = xregdrv(nregioni) + phsqksum*xrkk
            yregdrv(nregioni) = yregdrv(nregioni) + phsqksum*yrkk
            zregdrv(nregioni) = zregdrv(nregioni) + phsqksum*zrkk
          enddo
        elseif (lsiteenergy) then
          do ii = 1,ncharge
            i = nchargeptr(ii)
            oci = occuf(i)
            nregioni = nregionno(nsft+nrelf2a(i))
!
            esum = trmk*(csin(i)*csink(iv) + sine(i)*sinek(iv))
            eregion2region(nregioni,nregioni) = eregion2region(nregioni,nregioni) + esum
!
            siteenergy(i) = siteenergy(i) + 0.5_dp*esum
          enddo
        endif
!
!  End of loop over k vectors
!
      enddo
    endif
!
!  Multiply factor of half into erecip
!
    erecip = 0.5_dp*erecip
  endif
!**********************************
!  Bond order charge derivatives  *
!**********************************
  if (lgrad1.and.lDoQDeriv1) then
    do i = 1,numat
      oci = occuf(i)
      qli = qf(i)*oci
      fct = angstoev
      xci = xclat(i)
      yci = yclat(i)
      zci = zclat(i)
      lopi = (.not.lfreeze.or.lopf(nrelf2a(i)))
      do j = 1,i
        ocj = occuf(j)
        if (i.eq.j) then
          fct = 0.5_dp*fct
        endif
        qlj = qf(j)
        lopj = (.not.lfreeze.or.lopf(nrelf2a(j)))
!
!  Find relative vector between atoms
!
        xd = xclat(j) - xci
        yd = yclat(j) - yci
        zd = zclat(j) - zci
        do iv = 1,nlocalkvec
          xrkk = xrk(iv)
          yrkk = yrk(iv)
          zrkk = zrk(iv)
          arg = xrkk*xd + yrkk*yd + zrkk*zd
          cosa = cos(arg)*fct
          argc(iv) = cosa*oci*ocj*ktrm(iv)*qlj
          phsq(iv) = cosa*oci*ocj*ktrm(iv)*qli
        enddo
        call d1charge(i,j,lopi,lopj,nlocalkvec,argc,phsq)
      enddo
    enddo
  endif
!****************
!  Global sums  *
!****************
  tsum0 = g_cpu_time()
  if (lc6loc) then
    call sumone(erecip+ec6,esum,"recip3Dmd","ec6")
  else
    call sumone(erecip,esum,"recip3Dmd","erecip")
  endif       
  if (lsg1) then
    call sumall(strderv,sum,6_i4,"recip3Dmd","strderv")
    do i = 1,6
      strderv(i) = sum(i)
    enddo
  endif
  tsum = tsum + g_cpu_time() - tsum0
!****************************************************
!  Complete strain derivatives in reciprocal space  *
!****************************************************
  if (lsg1) then
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
  if (lpolar) then
!*******************************
!  Electric field calculation  *
!*******************************
!
!  Start loop over k vectors
!
    do iv = 1,nlocalkvec
      csink(iv) = 0.0_dp
      sinek(iv) = 0.0_dp
      xrkk = xrk(iv)
      yrkk = yrk(iv)
      zrkk = zrk(iv)
      trmk = ktrm(iv)*angstoev
      do i = 1,numat
        xci = xclat(i)
        yci = yclat(i)
        zci = zclat(i)
        argc(i) = xrkk*xci + yrkk*yci + zrkk*zci
        csin(i) = cos(argc(i))
        sine(i) = sin(argc(i))
      enddo
!
      do ii = 1,ncharge
        i = nchargeptr(ii)
        qli = qf(i)*occuf(i)
        csink(iv) = csink(iv) + csin(i)*qli
        sinek(iv) = sinek(iv) + sine(i)*qli
      enddo
!
      do i = 1,numat
!
!  Electric field
!
        phsqk = (csin(i)*sinek(iv) - sine(i)*csink(iv))
        phsqksum = trmk*phsqk
        vx(i) = vx(i) + phsqksum*xrkk
        vy(i) = vy(i) + phsqksum*yrkk
        vz(i) = vz(i) + phsqksum*zrkk
      enddo
!
!  End of loop over k vectors
!
    enddo
  endif
!***************
!  Exit point  *
!***************
999 continue
!
!  Free local memory
!
  deallocate(sum,stat=status)
  if (status/=0) call deallocate_error('recip3Dmd','sum')
  deallocate(ktrm62,stat=status)
  if (status/=0) call deallocate_error('recip3Dmd','ktrm62')
  deallocate(ktrm6,stat=status)
  if (status/=0) call deallocate_error('recip3Dmd','ktrm6')
  deallocate(ktrm4,stat=status)
  if (status/=0) call deallocate_error('recip3Dmd','ktrm4')
  deallocate(phsq,stat=status)
  if (status/=0) call deallocate_error('recip3Dmd','phsq')
  deallocate(sinek6,stat=status)
  if (status/=0) call deallocate_error('recip3Dmd','sinek6')
  deallocate(csink6,stat=status)
  if (status/=0) call deallocate_error('recip3Dmd','csink6')
  deallocate(sinek,stat=status)
  if (status/=0) call deallocate_error('recip3Dmd','sinek')
  deallocate(csink,stat=status)
  if (status/=0) call deallocate_error('recip3Dmd','csink')
!
!  Timing
!
  time1 = g_cpu_time()
  tion = tion + time1 - time0
#ifdef TRACE
  call trace_out('recip3Dmd')
#endif
!
  return
  end
