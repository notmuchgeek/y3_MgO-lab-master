  subroutine recipsd2(erecip,ec6,lgrad1,lgrad2)
!
!  Calculation of electrostatic potential and derivatives,
!  including strain derivatives. Reciprocal space part.
!  Vector version
!
!  Freezing now added.
!
!   4/95 Speed up added to avoid recalculation of cos and sine
!   5/95 Algorithm changed to aid the above - lveck is .false.
!        if either first or second derivatives are being calculated
!        as this appears to give large speed ups.
!   8/95 Ewald sum for 1/r**6 added
!   8/95 Derivatives for lveck removed as this is only used for the
!        energy - simplifies routine greatly!
!  11/96 Compression of second derivatives added when lfreeze=.true.
!        Because i(opt)-j(frozen) d2 blocks are stored in i(asym)-
!        i(full) block, it is necessary to exclude self-terms in the
!        the second derivatives.
!  12/97 Modification of energy/derivatives made purely local
!   4/98 ESFF Lennard-Jones form now allowed for
!   7/00 Dynamic memory allocation added
!   7/00 lfirst removed for safety
!   2/01 Electric field calculation added for polarisation
!   4/02 derv3 indexing for strain correction fixed for lfreeze case
!   9/04 Charge first derivatives added
!  10/04 oldel option removed as no longer used
!  11/04 Sqrt pi taken from module
!  12/07 Unused variables removed
!   1/09 Arrays apot, bpot, cpot & dpot replaced by new 2-D array twopot
!   5/09 Adding of contribution to sderv2 regardless of lgrad2 trapped
!   5/12 Atomic stresses added
!   5/12 Atomic stresses removed for routines involving symmetry
!  11/12 Unused variable removed
!  10/13 Hardwired maximum for k vector index replaced with maxindk
!   3/14 derfc changed to g_derfc for benefit of ChemShell
!   2/15 MM3buck added
!   2/18 Trace added
!   9/18 Handling of lstraincell algorithm added
!   9/18 Strain module introduced
!  11/18 Finite strain flag introduced instead of lstraincell
!   5/19 Modified so that strfin additions to sderv2 are no longer needed
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!  12/19 Rigid molecules added
!   1/20 Correction to com coordinate setting
!   3/20 Tolerance for ldoc6 made global
!   3/20 Use of charge pointer added
!   4/20 d2xyzdsdc added to cartstrterm arguments
!   4/20 derv3c added for benefit of rigid molecules
!   4/20 derv3c changes reversed as they are no longer required
!   7/20 Corrections to derv3 for C6 case where lc6one is false
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
  use kspace
  use m_strain,       only : gstrterms, strainddetds, straindet, straind2detds2, cartstrterm, gxyzterms
  use molecule
  use optimisation
  use polarise
  use potentialxyz
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
  logical,     intent(in)                     :: lgrad1
  logical,     intent(in)                     :: lgrad2
  real(dp),    intent(out)                    :: ec6
  real(dp),    intent(out)                    :: erecip
!
!  Local variables
!
  integer(i4)                                 :: e
  integer(i4)                                 :: f
  integer(i4)                                 :: g
  integer(i4)                                 :: i
  integer(i4)                                 :: idk
  integer(i4)                                 :: ii
  integer(i4)                                 :: iv
  integer(i4)                                 :: ix
  integer(i4)                                 :: iy
  integer(i4)                                 :: iz
  integer(i4)                                 :: ixf
  integer(i4)                                 :: iyf
  integer(i4)                                 :: izf
  integer(i4)                                 :: ixfo
  integer(i4)                                 :: iyfo
  integer(i4)                                 :: izfo
  integer(i4)                                 :: j
  integer(i4)                                 :: jx
  integer(i4)                                 :: jy
  integer(i4)                                 :: jz
  integer(i4)                                 :: jxc
  integer(i4)                                 :: jyc
  integer(i4)                                 :: jzc
  integer(i4)                                 :: k
  integer(i4)                                 :: kk
  integer(i4)                                 :: ll
  integer(i4)                                 :: n
  integer(i4)                                 :: nati
  integer(i4)                                 :: natj
  integer(i4)                                 :: nmi
  integer(i4)                                 :: nmj
  integer(i4)                                 :: nreli
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
  real(dp)                                    :: cosq
  real(dp)                                    :: g_cpu_time
  real(dp)                                    :: csin6
  real(dp)                                    :: csink
  real(dp)                                    :: csink6
  real(dp)                                    :: csinq
  real(dp)                                    :: csprod
  real(dp)                                    :: d1trm
  real(dp)                                    :: d2trm
  real(dp)                                    :: d3trm
  real(dp)                                    :: d21q
  real(dp)                                    :: d22q
  real(dp)                                    :: d23q
  real(dp)                                    :: d24q
  real(dp)                                    :: d25q
  real(dp)                                    :: d26q
  real(dp)                                    :: d3kk
  real(dp)                                    :: dGrds
  real(dp)                                    :: dGrds1
  real(dp)                                    :: dGrds2
  real(dp)                                    :: drxyzds(6,3)
  real(dp)                                    :: d2rxyzdsdx(6,3,3)
  real(dp)                                    :: d2rxyzds2(6,6,3)
  real(dp)                                    :: g_derfc
  real(dp)                                    :: eltrm
  real(dp)                                    :: erecipr
  real(dp)                                    :: esum
  real(dp)                                    :: factor
  real(dp), dimension(:),   allocatable       :: ktrm3
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
  real(dp)                                    :: rneqi
  real(dp)                                    :: rrk2
  real(dp)                                    :: sina
  real(dp)                                    :: sinek
  real(dp)                                    :: sinek6
  real(dp)                                    :: sineq
  real(dp)                                    :: sinqx
  real(dp)                                    :: sinqy
  real(dp)                                    :: sinqz
  real(dp)                                    :: strm1
  real(dp)                                    :: time0
  real(dp)                                    :: time1
  real(dp)                                    :: trmk
  real(dp)                                    :: trmk2
  real(dp)                                    :: trmk3
  real(dp)                                    :: trmk6
  real(dp)                                    :: trmk62
  real(dp)                                    :: trmk63
  real(dp)                                    :: trmks
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
  call trace_in('recipsd2')
#endif
!
  time0 = g_cpu_time()
!
!  Zero energy
!
  erecip = 0.0_dp
  ec6 = 0.0_dp
  lveck = (nkvec.ge.numat)
  lsg1 = (lstr.and.lgrad1)
  lc6loc = (lc6.and.ndim.eq.3)
!
!  Modification added to increase speed - outer loop
!  over k vectors is faster for derivatives as
!  recalculation of cos and sin is avoided.
!
  if (lgrad1.or.lgrad2) lveck = .false.
!
!  If Ewald sum for dispersion then don't use lveck
!  as this requires more vector storage
!
  if (lc6loc) lveck = .false.
!
!  Allocate local memory
!
  allocate(phsq(nkvec),stat=status)
  if (status/=0) call outofmemory('recipsd2','phsq')
  allocate(ktrm3(nkvec),stat=status)
  if (status/=0) call outofmemory('recipsd2','ktrm3')
  allocate(ktrm6(nkvec),stat=status)
  if (status/=0) call outofmemory('recipsd2','ktrm6')
  allocate(ktrm62(nkvec),stat=status)
  if (status/=0) call outofmemory('recipsd2','ktrm62')
  allocate(ktrm63(nkvec),stat=status)
  if (status/=0) call outofmemory('recipsd2','ktrm63')
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
      do i = 1,nkvec
        idk = indk(i)
        e = (idk/maxindk3) - maxindk
        if (e.eq.0) then
          factor = 1.0_dp
        else
          factor = 2.0_dp
        endif
        idk = idk - (e + maxindk)*maxindk3
        f = (idk/maxindk2) - maxindk
        g = idk - (f + maxindk)*maxindk2 - maxindk
        xrk(i) = e*kvv(1)
        yrk(i) = f*kvv(2)
        zrk(i) = g*kvv(3)
        rk2 = xrk(i)*xrk(i) + yrk(i)*yrk(i) + zrk(i)*zrk(i)
        arge = - rk2*eta4
        xpon = exp(arge)
        rrk2 = 1.0_dp/rk2
        ktrm(i) = xpon*vol4pi*factor*rrk2
        rk = sqrt(rk2)
        rketa2 = 0.5_dp*rk/seta
        c6t2 = sqrtpi*g_derfc(rketa2)
        rketa2 = 1.0_dp/rketa2
        c6t3 = 0.5_dp*rketa2**3 - rketa2
        c6t3 = c6t3*xpon
        c6t4 = c6t1*rk2*rk*factor
        ktrm6(i) = c6t4*(c6t2 + c6t3)
        if (lgrad1) then
          ktrms(i) = - 2.0_dp*ktrm(i)*(4.0_dp*eta + rk2)*eta4*rrk2
          ktrm62(i) = 3.0_dp*ktrm6(i)*rrk2
          ktrm62(i) = ktrm62(i) - c6t4*xpon*(12.0_dp*eta*seta*rrk2*rrk2)/rk
          if (lgrad2) then
            ktrm3(i) = (ktrm(i)*reta-4.0_dp*ktrms(i)*rrk2)
            ktrm63(i) = 3.0_dp*c6t4*c6t2*rrk2*rrk2
          endif
        endif
      enddo
    else
      do i = 1,nkvec
        idk = indk(i)
        e = (idk/maxindk3) - maxindk
        idk = idk - (e + maxindk)*maxindk3
        f = (idk/maxindk2) - maxindk
        g = idk - (f + maxindk)*maxindk2 - maxindk
        factor = 2.0_dp
        if (e.eq.0.and.nkangle.eq.1) then
          factor = 1.0_dp
        elseif (f.eq.0.and.nkangle.eq.2) then
          factor = 1.0_dp
        elseif (g.eq.0.and.nkangle.eq.3) then
          factor = 1.0_dp
        elseif (nkangle.eq.0) then
          factor = 1.0_dp
        endif
        xrk(i) = e*kv(1,1) + f*kv(1,2) + g*kv(1,3)
        yrk(i) = e*kv(2,1) + f*kv(2,2) + g*kv(2,3)
        zrk(i) = e*kv(3,1) + f*kv(3,2) + g*kv(3,3)
        rk2 = xrk(i)*xrk(i) + yrk(i)*yrk(i) + zrk(i)*zrk(i)
        arge = - rk2*eta4
        xpon = exp(arge)
        rrk2 = 1.0_dp/rk2
        ktrm(i) = xpon*vol4pi*factor*rrk2
        rk = sqrt(rk2)
        rketa2 = 0.5_dp*rk/seta
        c6t2 = sqrtpi*g_derfc(rketa2)
        rketa2 = 1.0_dp/rketa2
        c6t3 = 0.5_dp*rketa2**3 - rketa2
        c6t3 = c6t3*xpon
        c6t4 = c6t1*rk2*rk*factor
        ktrm6(i) = c6t4*(c6t2 + c6t3)
        if (lgrad1) then
          ktrms(i) = - 2.0_dp*ktrm(i)*(4.0_dp*eta + rk2)*eta4*rrk2
          ktrm62(i) = 3.0_dp*ktrm6(i)*rrk2
          ktrm62(i) = ktrm62(i) - c6t4*xpon*(12.0_dp*eta*seta*rrk2*rrk2)/rk
          if (lgrad2) then
            ktrm3(i) = (ktrm(i)*reta - 4.0_dp*ktrms(i)*rrk2)
            ktrm63(i) = 3.0_dp*c6t4*c6t2*rrk2*rrk2
          endif
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
      do i = 1,nkvec
        idk = indk(i)
        e = (idk/maxindk3) - maxindk
        if (e.eq.0) then
          factor = 1.0_dp
        else
          factor = 2.0_dp
        endif
        idk = idk - (e + maxindk)*maxindk3
        f = (idk/maxindk2) - maxindk
        g = idk - (f + maxindk)*maxindk2 - maxindk
        xrk(i) = e*kvv(1)
        yrk(i) = f*kvv(2)
        zrk(i) = g*kvv(3)
        rk2 = xrk(i)*xrk(i) + yrk(i)*yrk(i) + zrk(i)*zrk(i)
        arge = - rk2*eta4
        xpon = exp(arge)
        rrk2 = 1.0_dp/rk2
        ktrm(i) = xpon*vol4pi*factor*rrk2
        if (lgrad1) then
          ktrms(i) = - 2.0_dp*ktrm(i)*(4.0_dp*eta + rk2)*eta4*rrk2
          if (lgrad2) then
            ktrm3(i) = (ktrm(i)*reta - 4.0_dp*ktrms(i)*rrk2)
          endif
        endif
      enddo
    else
      do i = 1,nkvec
        idk = indk(i)
        e = (idk/maxindk3) - maxindk
        idk = idk - (e + maxindk)*maxindk3
        f = (idk/maxindk2) - maxindk
        g = idk - (f + maxindk)*maxindk2 - maxindk
        factor = 2.0_dp
        if (e.eq.0.and.nkangle.eq.1) then
          factor = 1.0_dp
        elseif (f.eq.0.and.nkangle.eq.2) then
          factor = 1.0_dp
        elseif (g.eq.0.and.nkangle.eq.3) then
          factor = 1.0_dp
        elseif (nkangle.eq.0) then
          factor = 1.0_dp
        endif
        xrk(i) = e*kv(1,1) + f*kv(1,2) + g*kv(1,3)
        yrk(i) = e*kv(2,1) + f*kv(2,2) + g*kv(2,3)
        zrk(i) = e*kv(3,1) + f*kv(3,2) + g*kv(3,3)
        rk2 = xrk(i)*xrk(i) + yrk(i)*yrk(i) + zrk(i)*zrk(i)
        arge = - rk2*eta4
        xpon = exp(arge)
        rrk2 = 1.0_dp/rk2
        ktrm(i) = xpon*vol4pi*factor*rrk2
        if (lgrad1) then
          ktrms(i) = - 2.0_dp*ktrm(i)*(4.0_dp*eta+rk2)*eta4*rrk2
          if (lgrad2) then
            ktrm3(i) = (ktrm(i)*reta - 4.0_dp*ktrms(i)*rrk2)
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
    do iv = 1,nkvec
      csin(iv) = 0.0_dp
      sine(iv) = 0.0_dp
    enddo
    do ii = 1,ncharge
      i = nchargeptr(ii)
      qli = qf(i)*occuf(i)
      xci = xclat(i)
      yci = yclat(i)
      zci = zclat(i)
      do iv = 1,nkvec
        argc(iv) = xrk(iv)*xci + yrk(iv)*yci+argc(iv) + zrk(iv)*zci+argc(iv)
      enddo
      do iv = 1,nkvec
        csin(iv) = csin(iv) + qli*cos(argc(iv))
        sine(iv) = sine(iv) + qli*sin(argc(iv))
      enddo
    enddo
    do iv = 1,nkvec
      phsq(iv) = (csin(iv)*csin(iv) + sine(iv)*sine(iv))
    enddo
!**********************
!  Lattice energy     *
!**********************
    erecipr = 0.0_dp
    do iv = 1,nkvec
      erecipr = erecipr + ktrm(iv)*phsq(iv)
    enddo
    erecip = erecip + 0.5_dp*angstoev*erecipr
  elseif ((lc6loc.and..not.lc6one).or.(lstr.and.lrigid).or.lDoQDeriv1) then
!************************************************************************************
!  Algorithm for cases where dispersion cannot be factorised into one centre terms  *
!  rigid molecule strains are to be calculated                                      *
!************************************************************************************
    if (lsg1.or.lgrad2) then
      call gstrterms(ndim,maxkvec,nkvec,xrk,yrk,zrk,dg2ds,d2g2dx2,d2g2ds2,lgrad2)
      if (lrigid) then
        call gxyzterms(ndim,maxkvec,nkvec,xrk,yrk,zrk,dgds,d2gds2,lgrad2)
      endif
    endif
    ix = - 2
    iy = - 1
    iz =   0
    ixf = 1
    iyf = 2
    izf = 3
    ixfo = 1
    iyfo = 2
    izfo = 3
    do i = 1,nasym
      oci = occua(i)*dble(neqv(i))*angstoev
      qli = qa(i)*oci
      nati = iatn(i)
      ntypi = natype(i)
      nmi = natmol(nrela2f(i))
      xci = xalat(i)
      yci = yalat(i)
      zci = zalat(i)
      nreli = nrela2f(i)
      lopi = lopf(i)
!
      if (lrigid.and.nmi.gt.0) then
        xcomi = molxyz(1,natinmol(nrela2f(i)),nmi)
        ycomi = molxyz(2,natinmol(nrela2f(i)),nmi)
        zcomi = molxyz(3,natinmol(nrela2f(i)),nmi)
      else
        xcomi = 0.0_dp
        ycomi = 0.0_dp
        zcomi = 0.0_dp
      endif
!
      if (.not.lfreeze.or.lopi) then
        ix = ix + 3
        iy = iy + 3
        iz = iz + 3
        ixf = ixfo
        iyf = iyfo
        izf = izfo
        ixfo = ixfo + 3*neqv(i)
        iyfo = iyfo + 3*neqv(i)
        izfo = izfo + 3*neqv(i)
      endif
      jx = - 2
      jy = - 1
      jz =   0
      jloop: do j = 1,numat
        ocj = occuf(j)
        qlj = qf(j)*ocj
        natj = nat(j)
        ntypj = nftype(j)
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
        lopj = lopf(nrelf2a(j))
        if (.not.lfreeze.or.lopj) then
          jx = jx + 3
          jy = jy + 3
          jz = jz + 3
        endif
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
        endif
        ldoc6 = (lc6loc.and.abs(c6tot).gt.thresh_c6)
!
!  Check whether this pair is worth doing
!
        qfct = qli*qlj
        if (abs(qfct).lt.thresh_q.and..not.ldoc6) cycle jloop
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
              cosq = cosa*ktrm(iv)
              cos6 = cosa*ktrm6(iv)*c6tot
              csinq = csinq + cosq
              cosq = cosq*qfct
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
                strm1 = 0.5_dp*strm1*cosa
                if (lrigid) then
                  do k = 1,6
                    dGrds = xrkk*drxyzds(k,1) + yrkk*drxyzds(k,2) + zrkk*drxyzds(k,3) + &
                            xd*dgds(iv,1,k) + yd*dgds(iv,2,k) + zd*dgds(iv,3,k)
                    strderv(k) = strderv(k) - 0.5_dp*d1trm*dGrds
                  enddo
                endif
                do k = 1,6
                  strderv(k) = strderv(k) + strm1*dg2ds(iv,k)
                enddo
                if (lgrad2) then
                  eltrm = 0.5_dp*(ktrm3(iv)*qfct-ktrm63(iv)*c6tot)*cosa
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
                      sderv2(ll,kk) = sderv2(ll,kk) + eltrm*dg2ds(iv,kk)*dg2ds(iv,ll) + strm1*d2g2ds2(iv,ll,kk)
                    enddo
                    if (.not.lfreeze.or.lopi) then
                      d3kk = d3trm*dg2ds(iv,kk)
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
            csin6 = csin6*c6tot
          endif
!*******************
!  Lattice energy  *
!*******************
          erecip = erecip + 0.5_dp*csinq*qfct
          ec6 = ec6 - 0.5_dp*(csin6 + c6self2*c6tot)
          esum = 0.5_dp*csinq*qfct - 0.5_dp*(csin6 + c6self2*c6tot)
!*****************************
!  Charge first derivatives  *
!*****************************
!              if (lgrad1.and.lDoQDeriv1) then
!                call d1charges(i,lopi,1_i4,0.5_dp*csinq*qlj)
!              endif
        else
          if (lgrad1) then
            do iv = 1,nkvec
              xrkk = xrk(iv)
              yrkk = yrk(iv)
              zrkk = zrk(iv)
              arg = xrkk*xd + yrkk*yd + zrkk*zd
              cosa = cos(arg)
              sina = sin(arg)*qfct
              cosq = cosa*ktrm(iv)
              cosa = cosa*qfct
              csinq = csinq + cosq
              cosq = cosq*qfct
              sineq = sina*ktrm(iv)
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
                strm1 = 0.5_dp*ktrms(iv)*cosa
                d3trm = ktrms(iv)*sina
                if (lrigid) then
                  do k = 1,6
                    dGrds = xrkk*drxyzds(k,1) + yrkk*drxyzds(k,2) + zrkk*drxyzds(k,3) + &
                            xd*dgds(iv,1,k) + yd*dgds(iv,2,k) + zd*dgds(iv,3,k)
                    strderv(k) = strderv(k) - 0.5_dp*sineq*dGrds
                  enddo
                endif
                do k = 1,6
                  strderv(k) = strderv(k) + strm1*dg2ds(iv,k)
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
                      sderv2(ll,kk) = sderv2(ll,kk) + eltrm*dg2ds(iv,kk)*dg2ds(iv,ll) + strm1*d2g2ds2(iv,ll,kk)
                    enddo
                    if (.not.lfreeze.or.lopi) then
                      d3kk = d3trm*dg2ds(iv,kk)
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
          endif
!*******************
!  Lattice energy  *
!*******************
          erecip = erecip + 0.5_dp*csinq*qfct
          esum = 0.5_dp*csinq*qfct
!*****************************
!  Charge first derivatives  *
!*****************************
!              if (lgrad1.and.lDoQDeriv1) then
!                call d1charges(i,lopi,1_i4,0.5_dp*csinq*qlj)
!              endif
        endif
!
!  Internal derivatives
!
        if (lgrad1.and.(.not.lfreeze.or.lopi)) then
          xdrv(i) = xdrv(i) + sinqx
          ydrv(i) = ydrv(i) + sinqy
          zdrv(i) = zdrv(i) + sinqz
          if (lgrad2.and.j.ne.nreli) then
            if (.not.lfreeze.or.lopj) then
              jxc = jx
              jyc = jy
              jzc = jz
            else
              jxc = ixf
              jyc = iyf
              jzc = izf
            endif
            derv2(jxc,ix) = derv2(jxc,ix) + d21q
            derv2(jyc,ix) = derv2(jyc,ix) + d26q
            derv2(jzc,ix) = derv2(jzc,ix) + d25q
            derv2(jxc,iy) = derv2(jxc,iy) + d26q
            derv2(jyc,iy) = derv2(jyc,iy) + d22q
            derv2(jzc,iy) = derv2(jzc,iy) + d24q
            derv2(jxc,iz) = derv2(jxc,iz) + d25q
            derv2(jyc,iz) = derv2(jyc,iz) + d24q
            derv2(jzc,iz) = derv2(jzc,iz) + d23q
          endif
        endif
      enddo jloop
    enddo
  else
!***********************************
!  Vectorise over number of atoms  *
!***********************************
    if (lsg1.or.lgrad2) then
      call gstrterms(ndim,maxkvec,nkvec,xrk,yrk,zrk,dg2ds,d2g2dx2,d2g2ds2,lgrad2)
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
      if (lgrad1) then
        trmk2 = ktrms(iv)*angstoev
        if (lgrad2) then
          trmk3 = ktrm3(iv)*angstoev
        endif
      endif
      if (lc6loc.and.lc6one) then
        trmk6 = ktrm6(iv)*angstoev
        if (lgrad1) then
          trmk62 = ktrm62(iv)*angstoev
          if (lgrad2) then
            trmk63 = ktrm63(iv)*angstoev
          endif
        endif
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
          csin(i) = cos(argc(i))*qli
          sine(i) = sin(argc(i))*qli
          csink = csink + csin(i)
          sinek = sinek + sine(i)
        enddo
      endif
      phsqk = (csink*csink + sinek*sinek)
!**********************
!  Lattice energy     *
!**********************
      erecip = erecip + 0.5_dp*trmk*phsqk
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
              sderv2(j,i) = sderv2(j,i) + eltrm*dg2ds(iv,i)*dg2ds(iv,j) + trmks*d2g2ds2(iv,j,i)
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
        ixf = 1
        iyf = 2
        izf = 3
        ixfo = 1
        iyfo = 2
        izfo = 3
        do i = 1,nasym
          lopi = lopf(i)
          if (.not.lfreeze.or.lopi) then
            rneqi = dble(neqv(i))
            oci = occua(i)*rneqi
            qli = qa(i)*oci
            nreli = nrela2f(i)
            if (lc6loc.and.lc6one) c6i = c6a(i)*oci
!
!  Excursion into second derivatives
!
            if (lgrad2) then
              ix = ix + 3
              iy = iy + 3
              iz = iz + 3
              ixf = ixfo
              iyf = iyfo
              izf = izfo
              ixfo = ixfo + 3*neqv(i)
              iyfo = iyfo + 3*neqv(i)
              izfo = izfo + 3*neqv(i)
              jx = - 2
              jy = - 1
              jz =   0
              do j = 1,numat
                lopj = lopf(nrelf2a(j))
                if (.not.lfreeze.or.lopj) then
                  jx = jx + 3
                  jy = jy + 3
                  jz = jz + 3
                  jxc = jx
                  jyc = jy
                  jzc = jz
                else
                  jxc = ixf
                  jyc = iyf
                  jzc = izf
                endif
                if (j.ne.nreli) then
                  ocj = occuf(j)
                  qlj = qf(j)*ocj
                  if (lc6loc.and.lc6one) c6j = c6f(j)*ocj
                  csprod = (csin(nreli)*csin(j) + sine(nreli)*sine(j))
                  if (lc6loc.and.lc6one) then
                    argck = csprod*(trmk*qli*qlj - trmk6*c6i*c6j)
                  else
                    argck = rneqi*trmk*csprod
                  endif
                  derv2(jxc,ix) = derv2(jxc,ix) + argck*d2g2dx2(iv,1)
                  derv2(jyc,ix) = derv2(jyc,ix) + argck*d2g2dx2(iv,6)
                  derv2(jzc,ix) = derv2(jzc,ix) + argck*d2g2dx2(iv,5)
                  derv2(jxc,iy) = derv2(jxc,iy) + argck*d2g2dx2(iv,6)
                  derv2(jyc,iy) = derv2(jyc,iy) + argck*d2g2dx2(iv,2)
                  derv2(jzc,iy) = derv2(jzc,iy) + argck*d2g2dx2(iv,4)
                  derv2(jxc,iz) = derv2(jxc,iz) + argck*d2g2dx2(iv,5)
                  derv2(jyc,iz) = derv2(jyc,iz) + argck*d2g2dx2(iv,4)
                  derv2(jzc,iz) = derv2(jzc,iz) + argck*d2g2dx2(iv,3)
                endif
              enddo
            endif
!
!  Return to first derivatives
!
            if (lc6loc.and.lc6one) then
              phsqk = qli*(csin(nreli)*sinek - sine(nreli)*csink)
              phsqk6 = c6i*(csin(nreli)*sinek6 - sine(nreli)*csink6)
              phsqksum = phsqk*trmk - phsqk6*trmk6
            else
              phsqk = rneqi*(csin(nreli)*sinek - sine(nreli)*csink)
              phsqksum = trmk*phsqk
            endif
            xdrv(i) = xdrv(i) + phsqksum*xrkk
            ydrv(i) = ydrv(i) + phsqksum*yrkk
            zdrv(i) = zdrv(i) + phsqksum*zrkk
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
    if (lgrad2) then
!
!  Volume corrections to mixed second derivatives
!
      ix = - 2
      iy = - 1
      iz =   0
      if (lfinitestrain) then
        do i = 1,numat
          lopi = lopf(nrelf2a(i))
          if (.not.lfreeze.or.lopi) then
            ix = ix + 3
            iy = iy + 3
            iz = iz + 3
            do j = 1,6
              derv3(ix,j) = derv3(ix,j) - xdrv(i)*strainddetds(j)*straindet
              derv3(iy,j) = derv3(iy,j) - ydrv(i)*strainddetds(j)*straindet
              derv3(iz,j) = derv3(iz,j) - zdrv(i)*strainddetds(j)*straindet
            enddo
          endif
        enddo
      else
        do i = 1,numat
          lopi = lopf(nrelf2a(i))
          if (.not.lfreeze.or.lopi) then
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
          endif
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
        xci = xclat(i)
        yci = yclat(i)
        zci = zclat(i)
        argc(i) = xrkk*xci + yrkk*yci + zrkk*zci
        csin(i) = cos(argc(i))
        sine(i) = sin(argc(i))
      enddo
      do ii = 1,ncharge
        i = nchargeptr(ii)
        qli = qf(i)*occuf(i)
        csink = csink + csin(i)*qli
        sinek = sinek + sine(i)*qli
      enddo
      do i = 1,nasym
!
!  Electric field
!
        nreli = nrela2f(i)
        phsqk = (csin(nreli)*sinek - sine(nreli)*csink)
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
  deallocate(ktrm63,stat=status)
  if (status/=0) call deallocate_error('recipsd2','ktrm63')
  deallocate(ktrm62,stat=status)
  if (status/=0) call deallocate_error('recipsd2','ktrm62')
  deallocate(ktrm6,stat=status)
  if (status/=0) call deallocate_error('recipsd2','ktrm6')
  deallocate(ktrm3,stat=status)
  if (status/=0) call deallocate_error('recipsd2','ktrm3')
  deallocate(phsq,stat=status)
  if (status/=0) call deallocate_error('recipsd2','phsq')
!
!  Timing
!
  time1 = g_cpu_time()
  tres = tres + time1 - time0
#ifdef TRACE
  call trace_out('recipsd2')
#endif
!
  return
  end
