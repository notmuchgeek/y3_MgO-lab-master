  subroutine recipsd(erecip,ec6,lgrad1)
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
!  12/97 Modification of energy/derivatives made purely local
!   4/98 ESFF Lennard-Jones form now allowed for
!   3/99 Parallel modifications added
!   7/00 Dynamic memory allocation added
!   7/00 lfirst removed for safety
!   2/01 Electric field calculation added for polarisation
!   4/01 Handle of K vector arrays in parallel corrected
!  11/02 Parallel changes made
!   9/04 Charge derivatives added
!  11/04 Sqrt pi taken from module
!  12/07 Unused variables removed
!   1/09 Arrays apot, bpot, cpot & dpot replaced by new 2-D array twopot
!   1/12 Syntax for passing energies to sumall changed for NAG compatibility. (EB)
!   5/12 Atomic stresses added
!   5/12 Atomic stresses removed for routines involving symmetry
!  11/12 Unused varible removed
!  10/13 Hardwired maximum for k vector index replaced with maxindk
!   3/14 derfc changed to g_derfc for benefit of ChemShell
!   2/15 MM3buck added
!   2/18 Trace added
!   9/18 Handling of lstraincell algorithm added
!   9/18 Strain module introduced
!  11/18 Finite strain flag introduced instead of lstraincell
!   7/19 Global sum of vx, vy, vz moved to polarisation
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!  12/19 Rigid molecules added
!   1/20 Correction to com coordinate setting
!   3/20 Tolerance for ldoc6 made global
!   3/20 Use of charge pointer added
!   4/20 d2xyzdsdc added to cartstrterm arguments
!   4/20 derv3c changes reversed as they are no longer required
!   7/20 Separate routine for sumall with 1 argument added
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
  use m_strain,       only : gstrterms, strainddetds, straindet, cartstrterm, gxyzterms
  use molecule
  use optimisation
  use parallel
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
  logical,    intent(in)                      :: lgrad1
  real(dp)                                    :: ec6
  real(dp)                                    :: erecip
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
  integer(i4)                                 :: j
  integer(i4)                                 :: k
  integer(i4)                                 :: kvec0
  integer(i4)                                 :: n
  integer(i4)                                 :: nati
  integer(i4)                                 :: natj
  integer(i4)                                 :: nlocalkvec
  integer(i4)                                 :: nmi
  integer(i4)                                 :: nmj
  integer(i4)                                 :: nprock
  integer(i4)                                 :: nreli
  integer(i4)                                 :: nremainder
  integer(i4)                                 :: ntypj
  integer(i4)                                 :: ntypi
  integer(i4)                                 :: status
  logical                                     :: lc6loc
  logical                                     :: ldoc6
  logical                                     :: lopi
  logical                                     :: lsg1
  logical                                     :: lveck
  real(dp)                                    :: arg
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
  real(dp)                                    :: cosq
  real(dp)                                    :: g_cpu_time
  real(dp)                                    :: csin6
  real(dp)                                    :: csink
  real(dp)                                    :: csink6
  real(dp)                                    :: csinq
  real(dp)                                    :: d1trm
  real(dp)                                    :: dGrds
  real(dp)                                    :: drxyzds(6,3)
  real(dp)                                    :: d2rxyzdsdx(6,3,3)
  real(dp)                                    :: d2rxyzds2(6,6,3)
  real(dp)                                    :: g_derfc
  real(dp)                                    :: esum
  real(dp)                                    :: factor
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
  real(dp), dimension(:),   allocatable       :: sum
  real(dp)                                    :: time0
  real(dp)                                    :: time1
  real(dp)                                    :: trmk
  real(dp)                                    :: trmk2
  real(dp)                                    :: trmk6
  real(dp)                                    :: trmk62
  real(dp)                                    :: trmks
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
  call trace_in('recipsd')
#endif
!
  time0 = g_cpu_time()
  lveck = (nkvec.ge.numat)
  lsg1 = (lstr.and.lgrad1)
  lc6loc = (lc6.and.ndim.eq.3)
!
!  Modification added to increase speed - outer loop
!  over k vectors is faster for derivatives as
!  recalculation of cos and sin is avoided.
!
  if (lgrad1) lveck = .false.
!
!  If Ewald sum for dispersion then don't use lveck
!  as this requires more vector storage
!
  if (lc6loc) lveck = .false.
!
!  Distribute kvec loops
!
  kvec0 = procid + 1
  nprock = nprocs
  nlocalkvec = (nkvec/nprocs)
  nremainder = nkvec - nlocalkvec*nprocs
  if (procid.lt.nremainder) nlocalkvec = nlocalkvec + 1
!
!  Allocate local memory
!
  allocate(phsq(nlocalkvec),stat=status)
  if (status/=0) call outofmemory('recipsd','phsq')
  allocate(ktrm6(nlocalkvec),stat=status)
  if (status/=0) call outofmemory('recipsd','ktrm6')
  allocate(ktrm62(nlocalkvec),stat=status)
  if (status/=0) call outofmemory('recipsd','ktrm62')
  allocate(sum(6_i4),stat=status)
  if (status/=0) call outofmemory('recipsd','sum')
!
!  Setup
!
  eta4 = 0.25_dp/eta
  rangstoev = 1.0_dp/angstoev
!*******************************
!  Sum 1/r**6 + coulomb terms  *
!*******************************
  if (lc6loc) then
    c6t1 = vol4pi*rangstoev*sqrtpi/48.0_dp
!
!  Reciprocal space self term
!
    c6self2 = 4.0_dp*c6t1*eta*seta/dble(nprocs)
    iv = 0
    if (lra) then
      kvv(1) = kv(1,1)
      kvv(2) = kv(2,2)
      kvv(3) = kv(3,3)
      do i = kvec0,nkvec,nprock
        iv = iv + 1
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
        xrk(iv) = e*kvv(1)
        yrk(iv) = f*kvv(2)
        zrk(iv) = g*kvv(3)
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
        if (lgrad1) then
          ktrms(iv) = - 2.0_dp*ktrm(iv)*(4.0_dp*eta + rk2)*eta4*rrk2
          ktrm62(iv) = 3.0_dp*ktrm6(iv)*rrk2
          ktrm62(iv) = ktrm62(iv) - c6t4*xpon*(12.0_dp*eta*seta*rrk2*rrk2)/rk
        endif
      enddo
    else
      do i = kvec0,nkvec,nprock
        iv = iv + 1
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
        xrk(iv) = e*kv(1,1) + f*kv(1,2) + g*kv(1,3)
        yrk(iv) = e*kv(2,1) + f*kv(2,2) + g*kv(2,3)
        zrk(iv) = e*kv(3,1) + f*kv(3,2) + g*kv(3,3)
        rk2 = xrk(iv)*xrk(iv) + yrk(iv)*yrk(iv) + zrk(iv)*zrk(iv)
        arge = - rk2*eta4
        xpon = exp(arge)
        rrk2 = 1.0_dp/rk2
        ktrm(iv) = xpon*vol4pi*factor*rrk2
        rk = sqrt(rk2)
        rketa2 = 0.5_dp*rk/seta
        c6t2 = sqrtpi*g_derfc(rketa2)
        rketa2 = 1.0_dp/rketa2
        c6t3 = 0.5_dp*rketa2**3-rketa2
        c6t3 = c6t3*xpon
        c6t4 = c6t1*rk2*rk*factor
        ktrm6(iv) = c6t4*(c6t2 + c6t3)
        if (lgrad1) then
          ktrms(iv) = - 2.0_dp*ktrm(iv)*(4.0_dp*eta+rk2)*eta4*rrk2
          ktrm62(iv) = 3.0_dp*ktrm6(iv)*rrk2
          ktrm62(iv) = ktrm62(iv) - c6t4*xpon*(12.0_dp*eta*seta*rrk2*rrk2)/rk
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
        e = (idk/maxindk3) - maxindk
        if (e.eq.0) then
          factor = 1.0_dp
        else
          factor = 2.0_dp
        endif
        idk = idk - (e + maxindk)*maxindk3
        f = (idk/maxindk2) - maxindk
        g = idk - (f + maxindk)*maxindk2 - maxindk
        xrk(iv) = e*kvv(1)
        yrk(iv) = f*kvv(2)
        zrk(iv) = g*kvv(3)
        rk2 = xrk(iv)*xrk(iv) + yrk(iv)*yrk(iv) + zrk(iv)*zrk(iv)
        arge = - rk2*eta4
        xpon = exp(arge)
        rrk2 = 1.0_dp/rk2
        ktrm(iv) = xpon*vol4pi*factor*rrk2
        if (lgrad1) then
          ktrms(iv) = - 2.0_dp*ktrm(iv)*(4.0_dp*eta+rk2)*eta4*rrk2
        endif
      enddo
    else
      do i = kvec0,nkvec,nprock
        iv = iv + 1
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
        xrk(iv) = e*kv(1,1) + f*kv(1,2) + g*kv(1,3)
        yrk(iv) = e*kv(2,1) + f*kv(2,2) + g*kv(2,3)
        zrk(iv) = e*kv(3,1) + f*kv(3,2) + g*kv(3,3)
        rk2 = xrk(iv)*xrk(iv) + yrk(iv)*yrk(iv) + zrk(iv)*zrk(iv)
        arge = - rk2*eta4
        xpon = exp(arge)
        rrk2 = 1.0_dp/rk2
        ktrm(iv) = xpon*vol4pi*factor*rrk2
        if (lgrad1) then
          ktrms(iv) = - 2.0_dp*ktrm(iv)*(4.0_dp*eta+rk2)*eta4*rrk2
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
    do ii = 1,ncharge
      i = nchargeptr(ii)
      qli = qf(i)*occuf(i)
      xci = xclat(i)
      yci = yclat(i)
      zci = zclat(i)
      do iv = 1,nlocalkvec
        argc(iv) = xrk(iv)*xci + yrk(iv)*yci + zrk(iv)*zci
      enddo
      do iv = 1,nlocalkvec
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
  elseif ((lc6loc.and..not.lc6one).or.(lstr.and.lrigid).or.lDoQDeriv1) then
!************************************************************************************
!  Algorithm for cases where dispersion cannot be factorised into one centre terms  *
!  rigid molecule strains are to be calculated                                      *
!************************************************************************************
    if (lsg1) then
      call gstrterms(ndim,maxkvec,nlocalkvec,xrk,yrk,zrk,dg2ds,d2g2dx2,d2g2ds2,.false.)
      if (lrigid) then
        call gxyzterms(ndim,maxkvec,nlocalkvec,xrk,yrk,zrk,dgds,d2gds2,.false.)
      endif
    endif
    do i = 1,nasym
      oci = occua(i)*angstoev*dble(neqv(i))
      qli = qa(i)*oci
      nati = iatn(i)
      ntypi = natype(i)
      nmi = natmol(nrela2f(i))
      xci = xalat(i)
      yci = yalat(i)
      zci = zalat(i)
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
        c6tot = 0.0_dp
        if (lc6loc) then
!
!  Find C6 term for pair
!
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
!  Check whether term is worth computing
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
          call cartstrterm(ndim,xd,yd,zd,xcom,ycom,zcom,drxyzds,d2rxyzdsdx,d2rxyzds2,.false.)
        endif
!
        qfct = qli*qlj
        csinq = 0.0_dp
        if (lgrad1) then
          sinqx = 0.0_dp
          sinqy = 0.0_dp
          sinqz = 0.0_dp
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
              cosq = cosa*ktrm(iv)
              cos6 = cosa*ktrm6(iv)
              csinq = csinq + cosq
              csin6 = csin6 + cos6
              d1trm = (ktrm(iv)*qfct - ktrm6(iv)*c6tot)*sina
              sinqx = sinqx + d1trm*xrkk
              sinqy = sinqy + d1trm*yrkk
              sinqz = sinqz + d1trm*zrkk
              if (lsg1) then
                strm1 = 0.5_dp*(ktrms(iv)*qfct - ktrm62(iv)*c6tot)*cosa
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
          endif
!*******************
!  Lattice energy  *
!*******************
          erecip = erecip + 0.5_dp*csinq*qfct
          ec6 = ec6 - 0.5_dp*(csin6 + c6self2)*c6tot
          esum = 0.5_dp*csinq*qfct - 0.5_dp*(csin6 + c6self2)*c6tot
!*****************************
!  Charge first derivatives  *
!*****************************
!              if (lgrad1.and.lDoQDeriv1) then
!                call d1charges(i,lopi,1_i4,0.5_dp*csinq*qlj*dble(neqv(i)))
!              endif
        else
          if (lgrad1) then
            do iv = 1,nlocalkvec
              xrkk = xrk(iv)
              yrkk = yrk(iv)
              zrkk = zrk(iv)
              arg = xrkk*xd + yrkk*yd + zrkk*zd
              cosa = cos(arg)
              sina = sin(arg)*qfct
              cosq = cosa*ktrm(iv)
              csinq = csinq + cosq
              sineq = sina*ktrm(iv)
              sinqx = sinqx + sineq*xrkk
              sinqy = sinqy + sineq*yrkk
              sinqz = sinqz + sineq*zrkk
              if (lsg1) then
                strm1 = 0.5_dp*ktrms(iv)*cosa*qfct
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
              endif
            enddo
          else
            csinq = 0.0_dp
            do iv = 1,nlocalkvec
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
!                call d1charges(i,lopi,1_i4,0.5_dp*csinq*qlj*dble(neqv(i)))
!              endif
        endif
!
!  Internal derivatives
!
        if (lgrad1.and.(.not.lfreeze.or.lopi)) then
          xdrv(i) = xdrv(i) + sinqx
          ydrv(i) = ydrv(i) + sinqy
          zdrv(i) = zdrv(i) + sinqz
        endif
      enddo jloop
    enddo
  else
!***********************************
!  Vectorise over number of atoms  *
!***********************************
    if (lsg1) then
      call gstrterms(ndim,maxkvec,nlocalkvec,xrk,yrk,zrk,dg2ds,d2g2dx2,d2g2ds2,.false.)
    endif
!
!  Start loop over k vectors
!
    do iv = 1,nlocalkvec
      csink = 0.0_dp
      sinek = 0.0_dp
      xrkk = xrk(iv)
      yrkk = yrk(iv)
      zrkk = zrk(iv)
      trmk = ktrm(iv)*angstoev
      if (lgrad1) then
        trmk2 = ktrms(iv)*angstoev
      endif
      if (lc6loc.and.lc6one) then
        trmk6 = ktrm6(iv)*angstoev
        if (lgrad1) then
          trmk62 = ktrm62(iv)*angstoev
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
        else
          trmks = 0.5_dp*trmk2*phsqk
        endif
        do i = 1,6
          strderv(i) = strderv(i) + trmks*dg2ds(iv,i)
        enddo
      endif
!************************
!  Internal derivatives *
!************************
      if (lgrad1) then
        do i = 1,nasym
          lopi = lopf(i)
          nreli = nrela2f(i)
          if (.not.lfreeze.or.lopi) then
            rneqi = dble(neqv(i))
            oci = occua(i)*rneqi
            qli = qa(i)*oci
            if (lc6loc.and.lc6one) then
              c6i = c6a(i)*oci
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
      do ii = 1,nchargec6
        i = nchargec6ptr(ii)
        c6tot = c6tot + c6f(i)*occuf(i)
      enddo
      ec6 = ec6 - 0.5_dp*c6self2*c6tot*c6tot*angstoev
    endif
  endif
!
!  Complete strain derivatives in reciprocal space
!
  if (lsg1) then
    tsum0 = g_cpu_time()
    if (lsg1) then
      call sumall(strderv,sum,6_i4,"recipsd","strderv")
      do i = 1,6
        strderv(i) = sum(i)
      enddo
    endif
    if (lc6loc) then
      call sumone(erecip+ec6,esum,"recipsd","erecip")
    else
      call sumone(erecip,esum,"recipsd","erecip")
    endif
    tsum = tsum + g_cpu_time() - tsum0
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
  deallocate(sum,stat=status)
  if (status/=0) call deallocate_error('recipsd','sum')
  deallocate(ktrm62,stat=status)
  if (status/=0) call deallocate_error('recipsd','ktrm62')
  deallocate(ktrm6,stat=status)
  if (status/=0) call deallocate_error('recipsd','ktrm6')
  deallocate(phsq,stat=status)
  if (status/=0) call deallocate_error('recipsd','phsq')
!
!  Timing
!
  time1 = g_cpu_time()
  tres = tres + time1 - time0
#ifdef TRACE
  call trace_out('recipsd')
#endif
!
  return
  end
