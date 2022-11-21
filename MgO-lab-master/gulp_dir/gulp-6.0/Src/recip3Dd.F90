  subroutine recip3Dd(erecip,ec6,lgrad1,lgrad2)
!
!  Calculation of electrostatic potential and derivatives,
!  including strain derivatives. Reciprocal space part.
!  Distributed memory second derivatives. No freezing.
!
!   5/13 Created from recip3D
!   9/16 cputime renamed to g_cpu_time
!   9/16 constants renamed to g_constants
!   9/16 Hardwired maximum for k vector index replaced with maxindk
!  12/16 derfc changed to g_derfc for benefit of ChemShell
!   7/17 Calls to d1charge changed to d1charged
!   7/17 Further debugging and changes for c6 algorithm tests
!   2/18 Trace added
!   9/18 Handling of lstraincell algorithm added
!   9/18 Strain module introduced
!  11/18 Finite strain flag introduced instead of lstraincell
!   5/19 Modified so that strfin additions to sderv2 are no longer needed
!   7/19 Global sum of vx, vy, vz moved to polarisation
!   8/19 Correction to sderv2 for c6 case
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!  12/19 Rigid molecules added
!   1/20 Correction to com coordinate setting
!   1/20 Correction to strain first derivatives for one algorithm
!   2/20 Second derivatives updated for rigid molecules
!   3/20 Tolerance for ldoc6 made global
!   3/20 Use of charge pointer added
!   4/20 d2xyzdsdc added to cartstrterm arguments
!   4/20 derv3c added for benefit of rigid molecules
!   4/20 derv3c changes reversed as they are no longer required
!   6/20 Bug in second derivatives for one algorithm fixed
!   6/20 Rigid molecule changes added to third algorithm to speed up calculation
!   6/20 Correction to derv3 for second algorithm with C6 terms
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
  use m_strain,       only : gstrterms, strainddetds, straindet, straind2detds2, cartstrterm, gxyzterms
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
  logical,  intent(in)                        :: lgrad1
  logical,  intent(in)                        :: lgrad2
  real(dp), intent(out)                       :: ec6
  real(dp), intent(out)                       :: erecip
!
!  Local variables
!
  integer(i4)                                 :: i
  integer(i4)                                 :: idk
  integer(i4)                                 :: ii
  integer(i4)                                 :: iloc
  integer(i4)                                 :: iv
  integer(i4)                                 :: ix
  integer(i4)                                 :: iy
  integer(i4)                                 :: iz
  integer(i4)                                 :: j
  integer(i4)                                 :: jj
  integer(i4)                                 :: jx
  integer(i4)                                 :: jy
  integer(i4)                                 :: jz
  integer(i4)                                 :: k
  integer(i4)                                 :: kk
  integer(i4)                                 :: kl
  integer(i4)                                 :: ks
  integer(i4)                                 :: kvec0
  integer(i4)                                 :: ll
  integer(i4)                                 :: n
  integer(i4)                                 :: nati
  integer(i4)                                 :: natj
  integer(i4)                                 :: nlocalkvec
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
  real(dp)                                    :: c6j
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
  real(dp)                                    :: csink
  real(dp)                                    :: csink6
  real(dp)                                    :: csinq
  real(dp)                                    :: csprod
  real(dp)                                    :: d1c
  real(dp)                                    :: d1c6
  real(dp)                                    :: d1s
  real(dp)                                    :: d1s6
  real(dp)                                    :: d1is(6)
  real(dp)                                    :: d1ix
  real(dp)                                    :: d1iy
  real(dp)                                    :: d1iz
  real(dp)                                    :: d1js(6)
  real(dp)                                    :: d1jx
  real(dp)                                    :: d1jy
  real(dp)                                    :: d1jz
  real(dp)                                    :: d1trm
  real(dp)                                    :: d2trm
  real(dp)                                    :: d2trmc
  real(dp)                                    :: d2trmc6
  real(dp)                                    :: d2trms
  real(dp)                                    :: d2trms6
  real(dp)                                    :: d21q
  real(dp)                                    :: d22q
  real(dp)                                    :: d23q
  real(dp)                                    :: d24q
  real(dp)                                    :: d25q
  real(dp)                                    :: d26q
  real(dp)                                    :: d2self
  real(dp)                                    :: d3kk
  real(dp)                                    :: d3trm
  real(dp)                                    :: dGrds
  real(dp)                                    :: dGrds1
  real(dp)                                    :: dGrds2
  real(dp)                                    :: drxyzds(6,3)
  real(dp)                                    :: d2rxyzdsdx(6,3,3)
  real(dp)                                    :: d2rxyzds2(6,6,3)
  real(dp)                                    :: g_derfc
  real(dp)                                    :: eltrm
  real(dp)                                    :: esum
  real(dp)                                    :: factor
  real(dp)                                    :: fct
  real(dp), dimension(:),   allocatable       :: ktrm3
  real(dp), dimension(:),   allocatable       :: ktrm4
  real(dp), dimension(:),   allocatable       :: ktrm6
  real(dp), dimension(:),   allocatable       :: ktrm62
  real(dp), dimension(:),   allocatable       :: ktrm63
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
  real(dp)                                    :: sdrvc(6)
  real(dp)                                    :: sdrvs(6)
  real(dp)                                    :: sdrvc6(6)
  real(dp)                                    :: sdrvs6(6)
  real(dp)                                    :: sina
  real(dp)                                    :: sini
  real(dp)                                    :: sinek
  real(dp)                                    :: sinek6
  real(dp)                                    :: sineq
  real(dp)                                    :: sinqx
  real(dp)                                    :: sinqy
  real(dp)                                    :: sinqz
  real(dp)                                    :: strdervloc(6)
  real(dp)                                    :: strm1
  real(dp), dimension(:),   allocatable       :: sum
  real(dp)                                    :: time0
  real(dp)                                    :: time1
  real(dp)                                    :: trmk
  real(dp)                                    :: trmks
  real(dp)                                    :: trmk2
  real(dp)                                    :: trmk3
  real(dp)                                    :: trmk6
  real(dp)                                    :: trmk62
  real(dp)                                    :: trmk63
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
  call trace_in('recip3Dd')
#endif
!
  time0 = g_cpu_time()
!
!  Initialise energies
!
  ec6 = 0.0_dp
  erecip = 0.0_dp
!
  lveck = (nkvec.ge.numat)
  lsg1 = (lstr.and.lgrad1)
  lc6loc = (lc6.and.ndim.eq.3)
!
!  Modification added to increase speed - outer loop over k vectors is faster for derivatives as
!  recalculation of cos and sin is avoided.
!
  if (lgrad1.or.lgrad2) lveck = .false.
!
!  If Ewald sum for dispersion then don't use lveck as this requires more vector storage
!
  if (lc6loc) lveck = .false.
!
!  If site energies are required then don't use lveck
!
  if (lsiteenergy) lveck = .false.
!
!  Distribute kvec loops for lveck algorithm
!
  if (lveck) then
    kvec0 = procid + 1
    nprock = nprocs
    nlocalkvec = (nkvec/nprocs)
    nremainder = nkvec - nlocalkvec*nprocs
    if (procid.lt.nremainder) nlocalkvec = nlocalkvec + 1
  else
    kvec0 = 1
    nprock = 1
    nlocalkvec = nkvec
  endif
!
!  Allocate local memory
!
  allocate(phsq(nlocalkvec),stat=status)
  if (status/=0) call outofmemory('recip3Dd','phsq')
  allocate(ktrm3(nlocalkvec),stat=status)
  if (status/=0) call outofmemory('recip3Dd','ktrm3')
  allocate(ktrm4(nlocalkvec),stat=status)
  if (status/=0) call outofmemory('recip3Dd','ktrm4')
  allocate(ktrm6(nlocalkvec),stat=status)
  if (status/=0) call outofmemory('recip3Dd','ktrm6')
  allocate(ktrm62(nlocalkvec),stat=status)
  if (status/=0) call outofmemory('recip3Dd','ktrm62')
  allocate(ktrm63(nlocalkvec),stat=status)
  if (status/=0) call outofmemory('recip3Dd','ktrm63')
  allocate(sum(nstrains),stat=status)
  if (status/=0) call outofmemory('recip3Dd','sum')
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
    c6self2 = 4.0_dp*c6t1*eta*seta/dble(nprock)
    iv = 0
    if (lra) then
      kvv(1) = kv(1,1)
      kvv(2) = kv(2,2)
      kvv(3) = kv(3,3)
      do i = kvec0,nkvec,nprock
        iv = iv + 1
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
          if (lgrad2) then
            ktrm3(iv) = (ktrm(iv)*reta - 4.0_dp*ktrms(iv)*rrk2)
            ktrm63(iv) = 3.0_dp*c6t4*c6t2*rrk2*rrk2
          endif
        endif
      enddo
    else
      do i = kvec0,nkvec,nprock
        iv = iv + 1
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
          if (lgrad2) then
            ktrm3(iv) = (ktrm(iv)*reta - 4.0_dp*ktrms(iv)*rrk2)
            ktrm63(iv) = 3.0_dp*c6t4*c6t2*rrk2*rrk2
          endif
        endif
      enddo
    endif
  else
!*****************
!  Coulomb only  *
!*****************
    iv = 0
    if (lra) then
      kvv(1) = kv(1,1)
      kvv(2) = kv(2,2)
      kvv(3) = kv(3,3)
      do i = kvec0,nkvec,nprock
        iv = iv + 1
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
          if (lgrad2) then
            ktrm3(iv) = (ktrm(iv)*reta - 4.0_dp*ktrms(iv)*rrk2)
          endif
        endif
      enddo
    else
      do i = kvec0,nkvec,nprock
        iv = iv + 1
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
          if (lgrad2) then
            ktrm3(iv) = (ktrm(iv)*reta - 4.0_dp*ktrms(iv)*rrk2)
          endif
        endif
      enddo
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
!
!  Start loop over cluster atom - unit cell atom pairs
!
    do iv = 1,nlocalkvec
      csin(iv) = 0.0_dp
      sine(iv) = 0.0_dp
    enddo
    do i = 1,numat
      qli = qf(i)*occuf(i)
      xci = xclat(i)
      yci = yclat(i)
      zci = zclat(i)
      do iv = 1,nlocalkvec
        argc(iv) = xrk(iv)*xci + yrk(iv)*yci + zrk(iv)*zci
        csin(iv) = csin(iv) + qli*cos(argc(iv))
        sine(iv) = sine(iv) + qli*sin(argc(iv))
      enddo
    enddo
    do iv = 1,nlocalkvec
      phsq(iv) = (csin(iv)*csin(iv) + sine(iv)*sine(iv))
    enddo
!**********************
!  Lattice energy     *
!**********************
    erecip = 0.0_dp
    do iv = 1,nlocalkvec
      erecip = erecip + ktrm(iv)*phsq(iv)
    enddo
    erecip = 0.5_dp*erecip*angstoev
  elseif ((lc6loc.and..not.lc6one).or.(lstr.and.latomicstress)) then
!************************************************************************************
!  Algorithm for cases where dispersion cannot be factorised into one centre terms  *
!  and/or atomic stresses or rigid molecule strains are to be calculated            *
!************************************************************************************
    if (lsg1.or.lgrad2) then
      call gstrterms(ndim,maxkvec,nkvec,xrk,yrk,zrk,dg2ds,d2g2dx2,d2g2ds2,lgrad2)
      if (lrigid) then
        call gxyzterms(ndim,maxkvec,nlocalkvec,xrk,yrk,zrk,dgds,d2gds2,lgrad2)
      endif
    endif
    if (lc6loc) then
!------------------
!  Charge and C6  |
!------------------
      ix = - 2
      iy = - 1
      iz =   0
      do ii = 1,natomsonnode
        i = node2atom(ii)
        oci = occuf(i)*angstoev
        qli = qf(i)*oci
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
        ix = ix + 3
        iy = iy + 3
        iz = iz + 3
        jx = - 2
        jy = - 1
        jz =   0
        do j = 1,numat
          ocj = occuf(j)
          qlj = qf(j)*ocj
          natj = nat(j)
          ntypj = nftype(j)
          nmj = natmol(j)
          nregionj = nregionno(nsft+nrelf2a(j))
          lopj = .true.
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
                endif
              endif
            elseif (natj.eq.nspec1(n).and.nati.eq.nspec2(n)) then
              if ((ntypj.eq.nptyp1(n).or.nptyp1(n).eq.0).and.(ntypi.eq.nptyp2(n).or.nptyp2(n).eq.0)) then
                if (nptype(n).eq.1.or.nptype(n).eq.7) then
                  c6tot = c6tot + twopot(3,n)
                elseif (nptype(n).eq.2.or.nptype(n).eq.21) then
                  c6tot = c6tot + twopot(2,n)
                endif
              endif
            endif
          enddo
          ldoc6 = (abs(c6tot).gt.thresh_c6)
!
!  Find relative vector between atoms
!
          xd = xclat(j) - xci
          yd = yclat(j) - yci
          zd = zclat(j) - zci
!
          if (lrigid.and.lsg1) then
            call cartstrterm(ndim,xd,yd,zd,xcom,ycom,zcom,drxyzds,d2rxyzdsdx,d2rxyzds2,lgrad2)
          endif
!
          qfct = qli*qlj
          csinq = 0.0_dp
          if (lgrad1) then
            sinqx = 0.0_dp
            sinqy = 0.0_dp
            sinqz = 0.0_dp
            if (lgrad2) then
              d21q = 0.0_dp
              d22q = 0.0_dp
              d23q = 0.0_dp
              d24q = 0.0_dp
              d25q = 0.0_dp
              d26q = 0.0_dp
            endif
            if (lsg1) then
              strdervloc(1:6) = 0.0_dp
            endif
          endif
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
              if (lgrad2) then
                d2trm = cosq - cos6
                d21q = d21q + d2trm*d2g2dx2(iv,1)
                d22q = d22q + d2trm*d2g2dx2(iv,2)
                d23q = d23q + d2trm*d2g2dx2(iv,3)
                d24q = d24q + d2trm*d2g2dx2(iv,4)
                d25q = d25q + d2trm*d2g2dx2(iv,5)
                d26q = d26q + d2trm*d2g2dx2(iv,6)
              endif
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
                if (lgrad2) then
                  eltrm = 0.5_dp*(ktrm3(iv)*qfct - ktrm63(iv)*c6tot)*cosa
                  do kk = 1,6
                    if (lrigid) then
                      dGrds1 = xrkk*drxyzds(kk,1) + yrkk*drxyzds(kk,2) + zrkk*drxyzds(kk,3) + &
                               xd*dgds(iv,1,kk) + yd*dgds(iv,2,kk) + zd*dgds(iv,3,kk)
                    endif
                    do ll = kk,6
                      if (lrigid) then
                        dGrds2 = xrkk*drxyzds(ll,1) + yrkk*drxyzds(ll,2) + zrkk*drxyzds(ll,3) + &
                                 xd*dgds(iv,1,ll) + yd*dgds(iv,2,ll) + zd*dgds(iv,3,ll)
                        sderv2(ll,kk) = sderv2(ll,kk) - 0.5_dp*(ktrm(iv)*qfct - ktrm6(iv)*c6tot)*(cosa*dGrds1*dGrds2 + sina*( &
                                                        xrkk*d2rxyzds2(ll,kk,1) + xd*d2gds2(iv,1,ll,kk) + &
                                                        yrkk*d2rxyzds2(ll,kk,2) + yd*d2gds2(iv,2,ll,kk) + &
                                                        zrkk*d2rxyzds2(ll,kk,3) + zd*d2gds2(iv,3,ll,kk) + &
                                                        dgds(iv,1,ll)*drxyzds(kk,1) + dgds(iv,1,kk)*drxyzds(ll,1) + &
                                                        dgds(iv,2,ll)*drxyzds(kk,2) + dgds(iv,2,kk)*drxyzds(ll,2) + &
                                                        dgds(iv,3,ll)*drxyzds(kk,3) + dgds(iv,3,kk)*drxyzds(ll,3)))
                        sderv2(ll,kk) = sderv2(ll,kk) - 0.5_dp*d3trm*(dGrds1*dg2ds(iv,ll) + dGrds2*dg2ds(iv,kk))
                      endif
                      sderv2(ll,kk) = sderv2(ll,kk) + eltrm*dg2ds(iv,kk)*dg2ds(iv,ll) + 0.5_dp*strm1*d2g2ds2(iv,ll,kk)
                    enddo
                    d3kk = d3trm*dg2ds(iv,kk)
                    if (lrigid) then
                      derv3(ix,kk) = derv3(ix,kk) + xrkk*d3kk + (cosq-cos6)*xrkk*dGrds1 + sineq*( &
                        xrkk*d2rxyzdsdx(kk,1,1) + dgds(iv,1,kk) + yrkk*d2rxyzdsdx(kk,2,1) + zrkk*d2rxyzdsdx(kk,3,1))
                      derv3(iy,kk) = derv3(iy,kk) + yrkk*d3kk + (cosq-cos6)*yrkk*dGrds1 + sineq*( &
                        xrkk*d2rxyzdsdx(kk,1,2) + yrkk*d2rxyzdsdx(kk,2,2) + dgds(iv,2,kk) + zrkk*d2rxyzdsdx(kk,3,2))
                      derv3(iz,kk) = derv3(iz,kk) + zrkk*d3kk + (cosq-cos6)*zrkk*dGrds1 + sineq*( &
                        xrkk*d2rxyzdsdx(kk,1,3) + yrkk*d2rxyzdsdx(kk,2,3) + zrkk*d2rxyzdsdx(kk,3,3) + dgds(iv,3,kk))
                    else
                      derv3(ix,kk) = derv3(ix,kk) + xrkk*d3kk
                      derv3(iy,kk) = derv3(iy,kk) + yrkk*d3kk
                      derv3(iz,kk) = derv3(iz,kk) + zrkk*d3kk
                    endif
                  enddo
                endif
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
          erecip = erecip + 0.5_dp*csinq
          ec6 = ec6 - 0.5_dp*(csin6 + c6self2*c6tot)
!
          esum = csinq - csin6 - c6self2*c6tot
          eregion2region(nregionj,nregioni) = eregion2region(nregionj,nregioni) + 0.5_dp*esum
!
          siteenergy(i) = siteenergy(i) + 0.5_dp*esum
!
          if (lsg1) then
!
!  Strain terms
!
            do kl = 1,nstrains
              strderv(kl) = strderv(kl) + 0.5_dp*strdervloc(kl)
            enddo
            if (latomicstress) then
              do kl = 1,nstrains
                atomicstress(kl,i) = atomicstress(kl,i) + 0.5_dp*strdervloc(kl)
              enddo
              do kl = 1,3
                atomicstress(kl,i) = atomicstress(kl,i) - 0.5_dp*esum
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
            if (nregioni.ne.nregionj) then
              xregdrv(nregioni) = xregdrv(nregioni) + sinqx
              yregdrv(nregioni) = yregdrv(nregioni) + sinqy
              zregdrv(nregioni) = zregdrv(nregioni) + sinqz
            endif
!
            if (lgrad2) then
              derv2(jx,ix) = derv2(jx,ix) + d21q
              derv2(jy,ix) = derv2(jy,ix) + d26q
              derv2(jz,ix) = derv2(jz,ix) + d25q
              derv2(jx,iy) = derv2(jx,iy) + d26q
              derv2(jy,iy) = derv2(jy,iy) + d22q
              derv2(jz,iy) = derv2(jz,iy) + d24q
              derv2(jx,iz) = derv2(jx,iz) + d25q
              derv2(jy,iz) = derv2(jy,iz) + d24q
              derv2(jz,iz) = derv2(jz,iz) + d23q
            endif
          endif
        enddo
      enddo
    else
!----------------
!  Charge only  |
!----------------
      ix = - 2
      iy = - 1
      iz =   0
      iiloop: do ii = 1,natomsonnode
        i = node2atom(ii)
        oci = occuf(i)*angstoev
        qli = qf(i)*oci
        ix = ix + 3
        iy = iy + 3
        iz = iz + 3
!
!  Check charge magnitude
!
        if (abs(qli).lt.thresh_q) cycle iiloop
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
        jx = - 2
        jy = - 1
        jz =   0
        do jj = 1,ncharge
          j = nchargeptr(jj)
          ocj = occuf(j)
          qlj = qf(j)*ocj
          natj = nat(j)
          ntypj = nftype(j)
          nmj = natmol(j)
          nregionj = nregionno(nsft+nrelf2a(j))
          lopj = .true.
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
          jx = 3*(j - 1) + 1
          jy = jx + 1
          jz = jx + 2
!
!  Find relative vector between atoms
!
          xd = xclat(j) - xci
          yd = yclat(j) - yci
          zd = zclat(j) - zci
!
          if (lrigid.and.lsg1) then
            call cartstrterm(ndim,xd,yd,zd,xcom,ycom,zcom,drxyzds,d2rxyzdsdx,d2rxyzds2,lgrad2)
          endif
!
          qfct = qli*qlj
          csinq = 0.0_dp
          if (lgrad1) then
            sinqx = 0.0_dp
            sinqy = 0.0_dp
            sinqz = 0.0_dp
            if (lgrad2) then
              d21q = 0.0_dp
              d22q = 0.0_dp
              d23q = 0.0_dp
              d24q = 0.0_dp
              d25q = 0.0_dp
              d26q = 0.0_dp
            endif
            if (lsg1) then
              strdervloc(1:6) = 0.0_dp
            endif
          endif
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
              if (lgrad2) then
                d21q = d21q + cosq*d2g2dx2(iv,1)
                d22q = d22q + cosq*d2g2dx2(iv,2)
                d23q = d23q + cosq*d2g2dx2(iv,3)
                d24q = d24q + cosq*d2g2dx2(iv,4)
                d25q = d25q + cosq*d2g2dx2(iv,5)
                d26q = d26q + cosq*d2g2dx2(iv,6)
              endif
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
                if (lgrad2) then
                  eltrm = 0.5_dp*ktrm3(iv)*cosa
                  do kk = 1,6
                    if (lrigid) then
                      dGrds1 = xrkk*drxyzds(kk,1) + yrkk*drxyzds(kk,2) + zrkk*drxyzds(kk,3) + &
                               xd*dgds(iv,1,kk) + yd*dgds(iv,2,kk) + zd*dgds(iv,3,kk)
                    endif
                    do ll = kk,6
                      if (lrigid) then
                        dGrds2 = xrkk*drxyzds(ll,1) + yrkk*drxyzds(ll,2) + zrkk*drxyzds(ll,3) + &
                                 xd*dgds(iv,1,ll) + yd*dgds(iv,2,ll) + zd*dgds(iv,3,ll)
                        sderv2(ll,kk) = sderv2(ll,kk) - 0.5_dp*(cosq*dGrds1*dGrds2 + sineq*( &
                                                        xrkk*d2rxyzds2(ll,kk,1) + xd*d2gds2(iv,1,ll,kk) + &
                                                        yrkk*d2rxyzds2(ll,kk,2) + yd*d2gds2(iv,2,ll,kk) + &
                                                        zrkk*d2rxyzds2(ll,kk,3) + zd*d2gds2(iv,3,ll,kk) + &
                                                        dgds(iv,1,ll)*drxyzds(kk,1) + dgds(iv,1,kk)*drxyzds(ll,1) + &
                                                        dgds(iv,2,ll)*drxyzds(kk,2) + dgds(iv,2,kk)*drxyzds(ll,2) + &
                                                        dgds(iv,3,ll)*drxyzds(kk,3) + dgds(iv,3,kk)*drxyzds(ll,3)))
                        sderv2(ll,kk) = sderv2(ll,kk) - 0.5_dp*d3trm*(dGrds1*dg2ds(iv,ll) + dGrds2*dg2ds(iv,kk))
                      endif
                      sderv2(ll,kk) = sderv2(ll,kk) + eltrm*dg2ds(iv,kk)*dg2ds(iv,ll) + 0.5_dp*strm1*d2g2ds2(iv,ll,kk)
                    enddo
                    d3kk = dg2ds(iv,kk)*d3trm
                    if (lrigid) then
                      derv3(ix,kk) = derv3(ix,kk) + xrkk*d3kk + cosq*xrkk*dGrds1 + sineq*( &
                        xrkk*d2rxyzdsdx(kk,1,1) + dgds(iv,1,kk) + yrkk*d2rxyzdsdx(kk,2,1) + zrkk*d2rxyzdsdx(kk,3,1))
                      derv3(iy,kk) = derv3(iy,kk) + yrkk*d3kk + cosq*yrkk*dGrds1 + sineq*( &
                        xrkk*d2rxyzdsdx(kk,1,2) + yrkk*d2rxyzdsdx(kk,2,2) + dgds(iv,2,kk) + zrkk*d2rxyzdsdx(kk,3,2))
                      derv3(iz,kk) = derv3(iz,kk) + zrkk*d3kk + cosq*zrkk*dGrds1 + sineq*( &
                        xrkk*d2rxyzdsdx(kk,1,3) + yrkk*d2rxyzdsdx(kk,2,3) + zrkk*d2rxyzdsdx(kk,3,3) + dgds(iv,3,kk))
                    else
                      derv3(ix,kk) = derv3(ix,kk) + xrkk*d3kk
                      derv3(iy,kk) = derv3(iy,kk) + yrkk*d3kk
                      derv3(iz,kk) = derv3(iz,kk) + zrkk*d3kk
                    endif
                  enddo
                endif
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
          erecip = erecip + 0.5_dp*csinq
!
          eregion2region(nregionj,nregioni) = eregion2region(nregionj,nregioni) + 0.5_dp*csinq
!
          esum = csinq
          siteenergy(i) = siteenergy(i) + 0.5_dp*esum
!
          if (lsg1) then
!
!  Strain terms
!
            do kl = 1,nstrains
              strderv(kl) = strderv(kl) + 0.5_dp*strdervloc(kl)
            enddo
            if (latomicstress) then
              do kl = 1,nstrains
                atomicstress(kl,i) = atomicstress(kl,i) + 0.5_dp*strdervloc(kl)
              enddo
              do kl = 1,3
                atomicstress(kl,i) = atomicstress(kl,i) - 0.5_dp*esum
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
            if (nregioni.ne.nregionj) then
              xregdrv(nregioni) = xregdrv(nregioni) + sinqx
              yregdrv(nregioni) = yregdrv(nregioni) + sinqy
              zregdrv(nregioni) = zregdrv(nregioni) + sinqz
            endif
!
            if (lgrad2) then
              derv2(jx,ix) = derv2(jx,ix) + d21q
              derv2(jy,ix) = derv2(jy,ix) + d26q
              derv2(jz,ix) = derv2(jz,ix) + d25q
              derv2(jx,iy) = derv2(jx,iy) + d26q
              derv2(jy,iy) = derv2(jy,iy) + d22q
              derv2(jz,iy) = derv2(jz,iy) + d24q
              derv2(jx,iz) = derv2(jx,iz) + d25q
              derv2(jy,iz) = derv2(jy,iz) + d24q
              derv2(jz,iz) = derv2(jz,iz) + d23q
            endif
          endif
        enddo
      enddo iiloop
    endif
  else
!***********************************
!  Vectorise over number of atoms  *
!***********************************
    if (lsg1.or.lgrad2) then
      call gstrterms(ndim,maxkvec,nkvec,xrk,yrk,zrk,dg2ds,d2g2dx2,d2g2ds2,lgrad2)
      if (lrigid) then
        call gxyzterms(ndim,maxkvec,nkvec,xrk,yrk,zrk,dgds,d2gds2,lgrad2)
      endif
    endif
!
!  Start loop over k vectors
!
    do iv = 1,nkvec
      csink = 0.0_dp
      sinek = 0.0_dp
      xrkk = xrk(iv)
      yrkk = yrk(iv)
      zrkk = zrk(iv)
      trmk = ktrm(iv)*angstoev
      if (lsg1) then
        trmk2 = ktrms(iv)*angstoev
        if (lgrad2) trmk3 = ktrm3(iv)*angstoev
      endif
      if (lc6loc.and.lc6one) then
        trmk6 = ktrm6(iv)*angstoev
        if (lsg1) trmk62 = ktrm62(iv)*angstoev
        if (lgrad2) trmk63 = ktrm63(iv)*angstoev
        csink6 = 0.0_dp
        sinek6 = 0.0_dp
        do ii = 1,nchargec6
          i = nchargec6ptr(ii)
          oci = occuf(i)
          qli = qf(i)*oci
          c6i = c6f(i)*oci
          xci = xclat(i)
          yci = yclat(i)
          zci = zclat(i)
          argc(i) = xrkk*xci + yrkk*yci + zrkk*zci
          csin(i) = cos(argc(i))
          sine(i) = sin(argc(i))
          csink6 = csink6 + csin(i)*c6i
          sinek6 = sinek6 + sine(i)*c6i
          csink = csink + csin(i)*qli
          sinek = sinek + sine(i)*qli
        enddo
        phsqk6 = (csink6*csink6 + sinek6*sinek6)
      else
        do ii = 1,ncharge
          i = nchargeptr(ii)
          qli = qf(i)*occuf(i)
          xci = xclat(i)
          yci = yclat(i)
          zci = zclat(i)
          argc(i) = xrkk*xci + yrkk*yci + zrkk*zci
          csin(i) = cos(argc(i))
          sine(i) = sin(argc(i))
          csink = csink + csin(i)*qli
          sinek = sinek + sine(i)*qli
        enddo
      endif
      phsqk = (csink*csink + sinek*sinek)
!**********************
!  Lattice energy     *
!**********************
!  erecip has a factor of a half but this multiplied later
      erecip = erecip + trmk*phsqk
      if (lc6loc.and.lc6one) ec6 = ec6 - 0.5_dp*trmk6*phsqk6
!**********************
!  Strain derivatives *
!**********************
      if (lsg1) then
        if (lc6loc.and.lc6one) then
          trmks = 0.5_dp*(trmk2*phsqk - trmk62*phsqk6)
          if (lgrad2) eltrm = 0.5_dp*(trmk3*phsqk - trmk63*phsqk6)
        else
          trmks = 0.5_dp*trmk2*phsqk
          if (lgrad2) eltrm = 0.5_dp*trmk3*phsqk
        endif
!
        if (lrigid) then
!
!  First loop to build extra term for strain derivatives
!
          sdrvc(1:6) = 0.0_dp
          sdrvs(1:6) = 0.0_dp
!
          if (lc6loc.and.lc6one) then
            sdrvc6(1:6) = 0.0_dp
            sdrvs6(1:6) = 0.0_dp
!
            do ii = 1,nchargec6
              i = nchargec6ptr(ii)
              nmi = natmol(i)
              if (nmi.gt.0) then
                qli = qf(i)*occuf(i)
                c6i = c6f(i)*oci
                d1c = qli*csin(i)
                d1s = qli*sine(i)
                d1c6 = c6i*csin(i)
                d1s6 = c6i*sine(i)
!
                xd = xclat(i)
                yd = yclat(i)
                zd = zclat(i)
!
                xcom = molxyz(1,natinmol(i),nmi)
                ycom = molxyz(2,natinmol(i),nmi)
                zcom = molxyz(3,natinmol(i),nmi)
!
                call cartstrterm(ndim,xd,yd,zd,xcom,ycom,zcom,drxyzds,d2rxyzdsdx,d2rxyzds2,lgrad2)
!
                do k = 1,6
                  dGrds = xrkk*drxyzds(k,1) + yrkk*drxyzds(k,2) + zrkk*drxyzds(k,3) + &
                          xd*dgds(iv,1,k) + yd*dgds(iv,2,k) + zd*dgds(iv,3,k)
                  sdrvc(k) = sdrvc(k) + d1c*dGrds
                  sdrvs(k) = sdrvs(k) + d1s*dGrds
                  sdrvc6(k) = sdrvc6(k) + d1c6*dGrds
                  sdrvs6(k) = sdrvs6(k) + d1s6*dGrds
                enddo
              endif
            enddo
          else
            do ii = 1,ncharge
              i = nchargeptr(ii)
              nmi = natmol(i)
              if (nmi.gt.0) then
                qli = qf(i)*occuf(i)
                d1c = qli*csin(i)
                d1s = qli*sine(i)
!
                xd = xclat(i)
                yd = yclat(i)
                zd = zclat(i)
!
                xcom = molxyz(1,natinmol(i),nmi)
                ycom = molxyz(2,natinmol(i),nmi)
                zcom = molxyz(3,natinmol(i),nmi)
!
                call cartstrterm(ndim,xd,yd,zd,xcom,ycom,zcom,drxyzds,d2rxyzdsdx,d2rxyzds2,lgrad2)
!
                do k = 1,6
                  dGrds = xrkk*drxyzds(k,1) + yrkk*drxyzds(k,2) + zrkk*drxyzds(k,3) + &
                          xd*dgds(iv,1,k) + yd*dgds(iv,2,k) + zd*dgds(iv,3,k)
                  sdrvc(k) = sdrvc(k) + d1c*dGrds
                  sdrvs(k) = sdrvs(k) + d1s*dGrds
                enddo
              endif
            enddo
          endif
!
!  Compute extra terms for strain derivatives from rigid molecules
!
          do ii = 1,nchargec6
            i = nchargec6ptr(ii)
            nmi = natmol(i)
            if (nmi.gt.0) then
              qli = qf(i)*occuf(i)
              d1trm = trmk*qli*(sinek*csin(i) - csink*sine(i))
              if (lc6loc.and.lc6one) then
                c6i = c6f(i)*oci
                d1trm = d1trm - trmk6*c6i*(sinek6*csin(i) - csink6*sine(i))
              endif
!
              xd = xclat(i)
              yd = yclat(i)
              zd = zclat(i)
!
              xcom = molxyz(1,natinmol(i),nmi)
              ycom = molxyz(2,natinmol(i),nmi)
              zcom = molxyz(3,natinmol(i),nmi)
!
              call cartstrterm(ndim,xd,yd,zd,xcom,ycom,zcom,drxyzds,d2rxyzdsdx,d2rxyzds2,lgrad2)
!
              do k = 1,6
                dGrds = xrkk*drxyzds(k,1) + yrkk*drxyzds(k,2) + zrkk*drxyzds(k,3) + &
                        xd*dgds(iv,1,k) + yd*dgds(iv,2,k) + zd*dgds(iv,3,k)
                strderv(k) = strderv(k) + d1trm*dGrds
              enddo
!
              if (lgrad2) then
                d2trm = trmk*qli*(- sinek*sine(i) - csink*csin(i))
                d2trmc = 0.5_dp*trmk*qli*csin(i)
                d2trms = 0.5_dp*trmk*qli*sine(i)
                d3trm = trmk2*qli*(sinek*csin(i) - csink*sine(i))
                if (lc6loc.and.lc6one) then
                  d2trm = d2trm - trmk6*c6i*(- sinek6*sine(i) - csink6*csin(i))
                  d2trmc6 = 0.5_dp*trmk6*c6i*csin(i)
                  d2trms6 = 0.5_dp*trmk6*c6i*sine(i)
                  d3trm = d3trm - trmk62*c6i*(sinek6*csin(i) - csink6*sine(i))
                endif
!
                do kk = 1,6
                  dGrds1 = xrkk*drxyzds(kk,1) + yrkk*drxyzds(kk,2) + zrkk*drxyzds(kk,3) + &
                           xd*dgds(iv,1,kk) + yd*dgds(iv,2,kk) + zd*dgds(iv,3,kk)
                  do ll = kk,6
                    dGrds2 = xrkk*drxyzds(ll,1) + yrkk*drxyzds(ll,2) + zrkk*drxyzds(ll,3) + &
                             xd*dgds(iv,1,ll) + yd*dgds(iv,2,ll) + zd*dgds(iv,3,ll)
                    sderv2(ll,kk) = sderv2(ll,kk) + d2trm*dGrds1*dGrds2
                    sderv2(ll,kk) = sderv2(ll,kk) + (d2trmc*sdrvc(kk) + d2trms*sdrvs(kk))*dGrds2
                    sderv2(ll,kk) = sderv2(ll,kk) + (d2trmc*sdrvc(ll) + d2trms*sdrvs(ll))*dGrds1
                    sderv2(ll,kk) = sderv2(ll,kk) + d1trm*( &
                                                    xrkk*d2rxyzds2(ll,kk,1) + xd*d2gds2(iv,1,ll,kk) + &
                                                    yrkk*d2rxyzds2(ll,kk,2) + yd*d2gds2(iv,2,ll,kk) + &
                                                    zrkk*d2rxyzds2(ll,kk,3) + zd*d2gds2(iv,3,ll,kk) + &
                                                    dgds(iv,1,ll)*drxyzds(kk,1) + dgds(iv,1,kk)*drxyzds(ll,1) + &
                                                    dgds(iv,2,ll)*drxyzds(kk,2) + dgds(iv,2,kk)*drxyzds(ll,2) + &
                                                    dgds(iv,3,ll)*drxyzds(kk,3) + dgds(iv,3,kk)*drxyzds(ll,3))
                    sderv2(ll,kk) = sderv2(ll,kk) + d3trm*(dGrds1*dg2ds(iv,ll) + dGrds2*dg2ds(iv,kk))
!
                    if (lc6loc.and.lc6one) then
                      sderv2(ll,kk) = sderv2(ll,kk) - (d2trmc6*sdrvc6(kk) + d2trms6*sdrvs6(kk))*dGrds2
                      sderv2(ll,kk) = sderv2(ll,kk) - (d2trmc6*sdrvc6(ll) + d2trms6*sdrvs6(ll))*dGrds1
                    endif
                  enddo
                enddo
!
!  Mixed derivatives
!
                iloc = atom2local(i)
                if (iloc.gt.0) then
                  ix = 3*(iloc-1) + 1
                  iy = ix + 1
                  iz = ix + 2
!
                  do kk = 1,6
                    dGrds = xrkk*drxyzds(kk,1) + yrkk*drxyzds(kk,2) + zrkk*drxyzds(kk,3) + &
                            xd*dgds(iv,1,kk) + yd*dgds(iv,2,kk) + zd*dgds(iv,3,kk)
!
                    derv3(ix,kk) = derv3(ix,kk) + 2.0_dp*xrkk*(d2trms*sdrvs(kk) + d2trmc*sdrvc(kk))
                    derv3(iy,kk) = derv3(iy,kk) + 2.0_dp*yrkk*(d2trms*sdrvs(kk) + d2trmc*sdrvc(kk))
                    derv3(iz,kk) = derv3(iz,kk) + 2.0_dp*zrkk*(d2trms*sdrvs(kk) + d2trmc*sdrvc(kk))
!
                    derv3(ix,kk) = derv3(ix,kk) + d2trm*xrkk*dGrds
                    derv3(iy,kk) = derv3(iy,kk) + d2trm*yrkk*dGrds
                    derv3(iz,kk) = derv3(iz,kk) + d2trm*zrkk*dGrds
!
                    derv3(ix,kk) = derv3(ix,kk) + d1trm*( &
                      xrkk*d2rxyzdsdx(kk,1,1) + dgds(iv,1,kk) + yrkk*d2rxyzdsdx(kk,2,1) + zrkk*d2rxyzdsdx(kk,3,1))
                    derv3(iy,kk) = derv3(iy,kk) + d1trm*( &
                      xrkk*d2rxyzdsdx(kk,1,2) + yrkk*d2rxyzdsdx(kk,2,2) + dgds(iv,2,kk) + zrkk*d2rxyzdsdx(kk,3,2))
                    derv3(iz,kk) = derv3(iz,kk) + d1trm*( &
                      xrkk*d2rxyzdsdx(kk,1,3) + yrkk*d2rxyzdsdx(kk,2,3) + zrkk*d2rxyzdsdx(kk,3,3) + dgds(iv,3,kk))
!
                    if (lc6loc.and.lc6one) then
                      derv3(ix,kk) = derv3(ix,kk) - 2.0_dp*xrkk*(d2trms6*sdrvs6(kk) + d2trmc6*sdrvc6(kk))
                      derv3(iy,kk) = derv3(iy,kk) - 2.0_dp*yrkk*(d2trms6*sdrvs6(kk) + d2trmc6*sdrvc6(kk))
                      derv3(iz,kk) = derv3(iz,kk) - 2.0_dp*zrkk*(d2trms6*sdrvs6(kk) + d2trmc6*sdrvc6(kk))
                    endif
                  enddo
                endif
              endif
            endif
          enddo
        endif
!
        do i = 1,6
          strderv(i) = strderv(i) + trmks*dg2ds(iv,i)
        enddo
        if (lgrad2) then
!
!  Second derivatives
!
!  General terms
!
          do i = 1,6
            do j = i,6
              sderv2(j,i) = sderv2(j,i) + eltrm*dg2ds(iv,j)*dg2ds(iv,i) + trmks*d2g2ds2(iv,j,i)
            enddo
          enddo
        endif
      endif
!************************
!  Internal derivatives *
!************************
!
!  First and second derivatives
!
      if (lgrad1) then
        ix = - 2
        iy = - 1
        iz =   0
        do ii = 1,natomsonnode
          i = node2atom(ii)
          oci = occuf(i)
          qli = qf(i)*oci
          nregioni = nregionno(nsft+nrelf2a(i))
          if (lc6loc.and.lc6one) c6i = c6f(i)*oci
!
!  Excursion into second derivatives
!
          if (lgrad2) then
            ix = ix + 3
            iy = iy + 3
            iz = iz + 3
            do jj = 1,nchargec6
              j = nchargec6ptr(jj)
              jx = 3*(j - 1) + 1
              jy = jx + 1
              jz = jx + 2
              ocj = occuf(j)
              qlj = qf(j)*ocj
              csprod = (csin(i)*csin(j) + sine(i)*sine(j))
              if (lc6loc.and.lc6one) then
                c6j = c6f(j)*ocj
                argck = csprod*(trmk*qli*qlj - trmk6*c6i*c6j)
              else
                argck = trmk*csprod*qli*qlj
              endif
              derv2(jx,ix) = derv2(jx,ix) + argck*d2g2dx2(iv,1)
              derv2(jy,iy) = derv2(jy,iy) + argck*d2g2dx2(iv,2)
              derv2(jz,iz) = derv2(jz,iz) + argck*d2g2dx2(iv,3)
              derv2(jy,ix) = derv2(jy,ix) + argck*d2g2dx2(iv,6)
              derv2(jz,ix) = derv2(jz,ix) + argck*d2g2dx2(iv,5)
              derv2(jz,iy) = derv2(jz,iy) + argck*d2g2dx2(iv,4)
            enddo
          endif
!
!  Return to first derivatives
!
          if (lc6loc.and.lc6one) then
            phsqk = qli*(csin(i)*sinek - sine(i)*csink)
            phsqk6 = c6i*(csin(i)*sinek6 - sine(i)*csink6)
            phsqksum = phsqk*trmk - phsqk6*trmk6
!
            esum = (trmk*qli*(csin(i)*csink + sine(i)*sinek) - &
                    trmk6*c6i*(csin(i)*csink6 + sine(i)*sinek6))
            eregion2region(nregioni,nregioni) = eregion2region(nregioni,nregioni) + esum
!
            siteenergy(i) = siteenergy(i) + 0.5_dp*esum
          else
            phsqk = qli*(csin(i)*sinek - sine(i)*csink)
            phsqksum = trmk*phsqk
!
            esum = trmk*(csin(i)*csink + sine(i)*sinek)
            eregion2region(nregioni,nregioni) = eregion2region(nregioni,nregioni) + esum
!
            siteenergy(i) = siteenergy(i) + 0.5_dp*esum
          endif
          xdrv(i) = xdrv(i) + phsqksum*xrkk
          ydrv(i) = ydrv(i) + phsqksum*yrkk
          zdrv(i) = zdrv(i) + phsqksum*zrkk
!
          xregdrv(nregioni) = xregdrv(nregioni) + phsqksum*xrkk
          yregdrv(nregioni) = yregdrv(nregioni) + phsqksum*yrkk
          zregdrv(nregioni) = zregdrv(nregioni) + phsqksum*zrkk
!
!  Mixed strain - internal derivatives
!
          if (lgrad2.and.lsg1) then
            if (lc6loc.and.lc6one) then
              argck = phsqk*trmk2 - phsqk6*trmk62
            else
              argck = phsqk*trmk2
            endif
            do kk = 1,6
              derv3(ix,kk) = derv3(ix,kk) + dg2ds(iv,kk)*xrkk*argck
              derv3(iy,kk) = derv3(iy,kk) + dg2ds(iv,kk)*yrkk*argck
              derv3(iz,kk) = derv3(iz,kk) + dg2ds(iv,kk)*zrkk*argck
            enddo
          endif
        enddo
      endif
!*******************************
!  End of loop over k vectors  *
!*******************************
    enddo
    if (lc6loc.and.lc6one) then
!
!  Self term
!
      c6tot = 0.0_dp
      do i = 1,numat
        c6tot = c6tot + c6f(i)*occuf(i)
      enddo
      ec6 = ec6 - 0.5_dp*c6self2*c6tot*c6tot*angstoev
    endif
    if (lgrad2) then
      ix = - 2
      iy = - 1
      iz =   0
      do ii = 1,natomsonnode
        ix = ix + 3
        iy = iy + 3
        iz = iz + 3
        jx = - 2
        jy = - 1
        jz =   0
        do j = 1,numat
          jx = jx + 3
          jy = jy + 3
          jz = jz + 3
          derv2(jx,iy) = derv2(jy,ix)
          derv2(jx,iz) = derv2(jz,ix)
          derv2(jy,iz) = derv2(jz,iy)
        enddo
      enddo
    endif
!
!  Multiply factor of half into erecip
!
    erecip = 0.5_dp*erecip
!
!  If nprocs > 1 then divide erecip & sderv2 by number processors as they are already the global sum
!
    if (nprocs.gt.1) then
      erecip = erecip/dble(nprocs)
      esum = esum/dble(nprocs)
      if (lc6loc) ec6 = ec6/dble(nprocs)
      if (lstr) then
        sderv2 = sderv2/dble(nprocs)
        strderv = strderv/dble(nprocs)
      endif
    endif
  endif
!**********************************
!  Bond order charge derivatives  *
!**********************************
  if (lgrad1.and.lDoQDeriv1) then
    do ii = 1,natomsonnode
      i = node2atom(ii)
      oci = occuf(i)
      qli = qf(i)*oci
      fct = angstoev
      xci = xclat(i)
      yci = yclat(i)
      zci = zclat(i)
      lopi = .true.
      do j = 1,numat
        ocj = occuf(j)
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
        enddo
        call d1charged(ii,i,lopi,nkvec,argc)
      enddo
    enddo
  endif
  if (lDoQDeriv2.and.lgrad2) then
!***********************************************************************
!  Calculation of charge derivative contribution for variable charges  *
!***********************************************************************
!
!  To save space :
!  d1i  is stored in ktrm3
!  d1j  is stored in ktrm4
!  d2i2 is stored in csin
!  d2ij is stored in argc 
!  d2j2 is stored in sine
!
    call gstrterms(ndim,maxkvec,nkvec,xrk,yrk,zrk,dg2ds,d2g2dx2,d2g2ds2,lgrad2)
    do iv = 1,nkvec
      csin(iv) = 0.0_dp
      sine(iv) = 0.0_dp
    enddo
    ix = - 2
    iy = - 1
    iz =   0
    do ii = 1,natomsonnode
      oci = occuf(i)
      qli = qf(i)*oci
      fct = angstoev
      xci = xclat(i)
      yci = yclat(i)
      zci = zclat(i)
      lopi = .true.
      ix = ix + 3
      iy = iy + 3
      iz = iz + 3
      jx = - 2
      jy = - 1
      jz =   0
      do j = 1,numat
        ocj = occuf(j)
        if (i.eq.j) then
          fct = 0.5_dp*fct
        endif
        qlj = qf(j)*ocj
        lopj = .true.
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
          arg = xrkk*xd + yrkk*yd + zrkk*zd
          cosa = cos(arg)*fct
          sina = sin(arg)*fct
! d2E/dqi.dqj
          argc(iv) = cosa*oci*ocj*ktrm(iv)
! d2E/dq.dr
          ktrm3(iv) = - ktrm(iv)*sina
! d2E/dqi.de
          ktrm6(iv) = - ktrms(iv)*oci*qlj*cosa
! d2E/dqj.de
          ktrm62(iv) = - ktrms(iv)*ocj*qli*cosa
! dE/dqi
          phsq(iv) = cosa*oci*qlj*ktrm(iv)
! dE/dqj
          ktrm4(iv) = cosa*qli*ocj*ktrm(iv)
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
        if (lstr) then
          do kl = 1,nstrains
            ks = nstrptr(kl)
            d1is(kl) = 0.0_dp
            d1js(kl) = 0.0_dp
            do iv = 1,nkvec
              d1is(kl) = d1is(kl) - ktrm6(iv)*dg2ds(iv,ks)
              d1js(kl) = d1js(kl) - ktrm62(iv)*dg2ds(iv,ks)
            enddo
          enddo
        endif
        call d2charge(i,j,nkvec,ix,iy,iz,jx,jy,jz,lopi,lopj,phsq,ktrm4,d1ix,d1iy,d1iz, &
                      d1jx,d1jy,d1jz,d1is,d1js,csin,argc,sine,d2self,0.0_dp,0.0_dp,.false., &
                      .false.)
      enddo
    enddo
  endif
!***********************************************************
!  Complete strain second derivatives in reciprocal space  *
!***********************************************************
  if (lc6loc) then
    esum = erecip + ec6
  else
    esum = erecip
  endif
  if (lsg1) then
    if (lgrad2) then
!
!  Volume corrections to mixed second derivatives
!
      ix = - 2
      iy = - 1
      iz =   0
      if (lfinitestrain) then
        do ii = 1,natomsonnode
          i = node2atom(ii)
          ix = ix + 3
          iy = iy + 3
          iz = iz + 3
          do j = 1,6
            derv3(ix,j) = derv3(ix,j) - xdrv(i)*strainddetds(j)*straindet
            derv3(iy,j) = derv3(iy,j) - ydrv(i)*strainddetds(j)*straindet
            derv3(iz,j) = derv3(iz,j) - zdrv(i)*strainddetds(j)*straindet
          enddo
        enddo
      else
        do ii = 1,natomsonnode
          i = node2atom(ii)
          ix = ix + 3
          iy = iy + 3
          iz = iz + 3
          derv3(ix,1) = derv3(ix,1) - xdrv(i)
          derv3(iy,1) = derv3(iy,1) - ydrv(i)
          derv3(iz,1) = derv3(iz,1) - zdrv(i)
          derv3(ix,2) = derv3(ix,2) - xdrv(i)
          derv3(iy,2) = derv3(iy,2) - ydrv(i)
          derv3(iz,2) = derv3(iz,2) - zdrv(i)
          derv3(ix,3) = derv3(ix,3) - xdrv(i)
          derv3(iy,3) = derv3(iy,3) - ydrv(i)
          derv3(iz,3) = derv3(iz,3) - zdrv(i)
        enddo
      endif
!
!  Volume corrections to strain second derivatives
!
      if (lfinitestrain) then
        do i = 1,6
          do j = 1,6
            sderv2(j,i) = sderv2(j,i) - strderv(j)*strainddetds(i)*straindet &
                                      - strderv(i)*strainddetds(j)*straindet &
                                      + 2.0_dp*esum*strainddetds(j)*strainddetds(i)*straindet**2 &
                                      - esum*straind2detds2(j,i)*straindet
          enddo
        enddo
      else
        do i = 1,3
          do j = 1,3
            sderv2(j,i) = sderv2(j,i) - strderv(j) - strderv(i)
            sderv2(j,i) = sderv2(j,i) + esum
          enddo
          do j = 4,6
            sderv2(j,i) = sderv2(j,i) - strderv(j)
            sderv2(i,j) = sderv2(i,j) - strderv(j)
          enddo
        enddo
      endif
    endif
    if (lfinitestrain) then
!
!  Volume first derivatives
!
      do i = 1,6
        strderv(i) = strderv(i) - strainddetds(i)*straindet*esum
      enddo
    else
      strderv(1) = strderv(1) - esum
      strderv(2) = strderv(2) - esum
      strderv(3) = strderv(3) - esum
    endif
  endif
!****************
!  Global sums  *
!****************
  if (.not.lveck) then
    tsum0 = g_cpu_time()
    if (lc6loc) then
      call sumone(erecip+ec6,esum,"recip3Dd","ec6")
    else
      call sumone(erecip,esum,"recip3Dd","erecip")
    endif
    if (lsg1) then
      call sumall(strderv,sum,6_i4,"recip3Dd","strderv")
      do i = 1,6
        strderv(i) = sum(i)
      enddo
    endif
    tsum = tsum + g_cpu_time() - tsum0
  elseif (lsg1) then
    if (lc6loc.and.lc6loc) then
      esum = erecip + ec6
    else
      esum = erecip
    endif
  endif
  if (lpolar) then
!*******************************
!  Electric field calculation  *
!*******************************
!
!  Start loop over k vectors
!
    do iv = 1,nkvec
      csink = 0.0_dp
      sinek = 0.0_dp
      xrkk = xrk(iv)
      yrkk = yrk(iv)
      zrkk = zrk(iv)
      trmk = ktrm(iv)*angstoev
      do i = 1,numat
        qli = qf(i)*occuf(i)
        xci = xclat(i)
        yci = yclat(i)
        zci = zclat(i)
        argc(i) = xrkk*xci + yrkk*yci + zrkk*zci
        cosi = cos(argc(i))
        sini = sin(argc(i))
        csin(i) = cosi
        sine(i) = sini
        csink = csink + cosi*qli
        sinek = sinek + sini*qli
      enddo
      do ii = 1,natomsonnode
        i = node2atom(ii)
!
!  Electric field
!
        phsqk = (csin(i)*sinek - sine(i)*csink)
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
  if (status/=0) call deallocate_error('recip3Dd','sum')
  deallocate(ktrm63,stat=status)
  if (status/=0) call deallocate_error('recip3Dd','ktrm63')
  deallocate(ktrm62,stat=status)
  if (status/=0) call deallocate_error('recip3Dd','ktrm62')
  deallocate(ktrm6,stat=status)
  if (status/=0) call deallocate_error('recip3Dd','ktrm6')
  deallocate(ktrm4,stat=status)
  if (status/=0) call deallocate_error('recip3Dd','ktrm4')
  deallocate(ktrm3,stat=status)
  if (status/=0) call deallocate_error('recip3Dd','ktrm3')
  deallocate(phsq,stat=status)
  if (status/=0) call deallocate_error('recip3Dd','phsq')
!
!  Timing
!
  time1 = g_cpu_time()
  tion = tion + time1 - time0
#ifdef TRACE
  call trace_out('recip3Dd')
#endif
!
  return
  end
