  subroutine recip2D(erecip,esregion12,esregion2,eattach,lgrad1,lgrad2)
!
!  Calculation of electrostatic potential and derivatives,
!  including strain derivatives. Reciprocal space part.
!  Two-dimensional version using the Parry sum.
!
!   2/01 Created from recip
!   2/01 Proper cut-offs added and i=j term handled
!   2/01 Sorting of kvectors for efficient screening added
!   4/01 Calculation of region2 self energy term added
!  11/01 Attachment energy added
!   4/02 derv3 indexing for strain correction fixed for lfreeze case
!   8/02 Surface energy calculation algorithm changed
!  11/02 Parallel modifications made
!   9/04 Charge first derivatives added
!   9/04 Setting of iz fixed for derivative subtraction from derv3
!   6/07 Parallel bug in derivatives fixed due to dtrm1 being added on
!        all nodes and then summed.
!   6/09 Site energy added
!  11/09 Region derivatives added
!  11/11 Region-region energy contributions stored
!   4/12 xvir, yvir and zvir removed
!   5/12 Atomic stresses added
!  10/13 Hardwired maximum for k vector index replaced with maxindk
!   3/14 derfc changed to g_derfc for benefit of ChemShell
!   2/18 Trace added
!   9/18 Handling of lstraincell algorithm added
!   9/18 Strain module introduced
!  11/18 Finite strain flag introduced instead of lstraincell
!   1/19 Charge second derivatives modified for lstraincell
!   5/19 Modified so that strfin additions to sderv2 are no longer needed
!   7/19 Global sum of vx, vy, vz moved to polarisation
!   7/19 esum replaces erecip in strain derivatives as this requires
!        the global sum in parallel
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!  12/19 Rigid molecules added
!   1/20 Correction to com coordinate setting
!   3/20 Second derivatives added for rigid molecules
!   3/20 Use of charge pointer added
!   4/20 d2xyzdsdc added to cartstrterm arguments
!   4/20 derv3c added for benefit of rigid molecules
!   4/20 derv3c changes reversed as they are no longer required
!   6/20 Shells excluded from centre of mass phasing
!   6/20 Wrapping of derv3 with lopi/lopj
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
  use configurations, only : lsliceatom, nregions, nregionno
  use g_constants
  use control
  use current
  use derivatives
  use element
  use energies,       only : siteenergy, eregion2region
  use kspace
  use m_strain,       only : gstrterms, strainddetds, straindet, straind2detds2, cartstrterm, gxyzterms
  use molecule
  use optimisation
  use parallel
  use polarise
  use potentialxyz
  use symmetry
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  use two
  implicit none
!
!  Passed variables
!
  real(dp),    intent(inout)                     :: erecip
  real(dp),    intent(inout)                     :: esregion12
  real(dp),    intent(inout)                     :: esregion2
  real(dp),    intent(inout)                     :: eattach
  logical,     intent(in)                        :: lgrad1
  logical,     intent(in)                        :: lgrad2
!
!  Local variables
!
  integer(i4)                                    :: i
  integer(i4)                                    :: idk
  integer(i4)                                    :: ii
  integer(i4)                                    :: imin
  integer(i4)                                    :: iv
  integer(i4)                                    :: ix
  integer(i4)                                    :: iy
  integer(i4)                                    :: iz
  integer(i4)                                    :: ixc
  integer(i4)                                    :: iyc
  integer(i4)                                    :: izc
  integer(i4)                                    :: j
  integer(i4)                                    :: jj
  integer(i4)                                    :: jx
  integer(i4)                                    :: jy
  integer(i4)                                    :: jz
  integer(i4)                                    :: jxc
  integer(i4)                                    :: jyc
  integer(i4)                                    :: jzc
  integer(i4)                                    :: kk
  integer(i4)                                    :: kvec0
  integer(i4)                                    :: ll
  integer(i4)                                    :: nlocalkvec
  integer(i4)                                    :: nmi
  integer(i4)                                    :: nmj
  integer(i4)                                    :: nprock
  integer(i4)                                    :: nremainder
  integer(i4)                                    :: nregioni
  integer(i4)                                    :: nregionj
  integer(i4)                                    :: status
  logical                                        :: lattach
  logical                                        :: lopi
  logical                                        :: lopj
  logical                                        :: lreg2one
  logical                                        :: lreg2pair
  logical                                        :: lsg1
  logical                                        :: lslicei
  logical                                        :: lslicej
  real(dp)                                       :: accf2
  real(dp)                                       :: arg
  real(dp)                                       :: argtest
  real(dp)                                       :: cosa
  real(dp)                                       :: cosq
  real(dp)                                       :: g_cpu_time
  real(dp)                                       :: csinq
  real(dp)                                       :: d0term
  real(dp)                                       :: d21q
  real(dp)                                       :: d22q
  real(dp)                                       :: d23q
  real(dp)                                       :: d24q
  real(dp)                                       :: d25q
  real(dp)                                       :: d26q
  real(dp)                                       :: d2self
  real(dp)                                       :: d2trm1dij
  real(dp)                                       :: d2trm1diz
  real(dp)                                       :: d2trm1djz
  real(dp)                                       :: d3kxy
  real(dp)                                       :: d3kz
  real(dp)                                       :: darg1
  real(dp)                                       :: darg2
  real(dp)                                       :: g_derf
  real(dp)                                       :: g_derfc
  real(dp)                                       :: d1is(3)
  real(dp)                                       :: d1ix
  real(dp)                                       :: d1iy
  real(dp)                                       :: d1iz
  real(dp)                                       :: d1js(3)
  real(dp)                                       :: d1jx
  real(dp)                                       :: d1jy
  real(dp)                                       :: d1jz
  real(dp)                                       :: derfc1
  real(dp)                                       :: derfc2
  real(dp)                                       :: derfez
  real(dp)                                       :: dexp1
  real(dp)                                       :: dexp2
  real(dp)                                       :: dexp3
  real(dp)                                       :: dexp4
  real(dp)                                       :: dexpz
  real(dp)                                       :: dtrm1
  real(dp)                                       :: dtrm1di
  real(dp)                                       :: dtrm1dj
  real(dp)                                       :: dtrm2
  real(dp)                                       :: dGrds
  real(dp)                                       :: dGrds1
  real(dp)                                       :: dGrds2
  real(dp)                                       :: d2Grds12
  real(dp)                                       :: drxyzds(6,3)
  real(dp)                                       :: d2rxyzdsdx(6,3,3)
  real(dp)                                       :: d2rxyzds2(6,6,3)
  real(dp)                                       :: eltrm
  real(dp)                                       :: eltrm1
  real(dp)                                       :: eltrm2
  real(dp)                                       :: esum
  real(dp)                                       :: etaz
  real(dp)                                       :: etaz2
  real(dp)                                       :: etrm
  real(dp)                                       :: eztrm
  real(dp)                                       :: factor
  real(dp)                                       :: fct
  real(dp)                                       :: Gmax
  real(dp)                                       :: Gmin
  real(dp)                                       :: gseta
  real(dp)                                       :: kexperfc
  real(dp)                                       :: kvec
  real(dp)                                       :: kvv(3)
  real(dp),    dimension(:),   allocatable       :: ktrm3
  real(dp),    dimension(:),   allocatable       :: ktrm3z
  real(dp),    dimension(:),   allocatable       :: ktrm6
  real(dp),    dimension(:),   allocatable       :: ktrm62
  real(dp),    dimension(:),   allocatable       :: ktrm63
  real(dp)                                       :: oci
  real(dp)                                       :: ocj
  real(dp),    dimension(:),   allocatable       :: phsq
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
  real(dp),    dimension(:),   allocatable       :: sum
  real(dp)                                       :: sztrm2
  real(dp)                                       :: time0
  real(dp)                                       :: time1
  real(dp)                                       :: tsum0
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
  call trace_in('recip2D')
#endif
!
  time0 = g_cpu_time()
  lsg1 = (lstr.and.lgrad1)
  accf2 = accf*accf
  argtest = sqrt(3.0_dp+0.5_dp*accf2) - sqrt(3.0_dp)
  smallestG = min(abs(kv(1,1)),abs(kv(2,2)))
!
!  Distribute kvec loops
!
  if (lgrad2) then
    kvec0 = 1
    nprock = 1
    nlocalkvec = nkvec
  else
    kvec0 = procid + 1
    nprock = nprocs
    nlocalkvec = (nkvec/nprocs)
    nremainder = nkvec - nlocalkvec*nprocs
    if (procid.lt.nremainder) nlocalkvec = nlocalkvec + 1
  endif
!
!  Allocate local memory
!
  allocate(sum(nstrains),stat=status)
  if (status/=0) call outofmemory('recip2D','sum')
!**********
!  Setup  *
!**********
  iv = 0
  if (lra) then
    kvv(1) = kv(1,1)
    kvv(2) = kv(2,2)
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
      xrk(iv) = ii*kvv(1)
      yrk(iv) = jj*kvv(2)
      rk2 = xrk(iv)*xrk(iv) + yrk(iv)*yrk(iv)
      kmod(iv) = sqrt(rk2)
      ktrm(iv) = 0.5_dp*vol4pi*factor/kmod(iv)
    enddo
  else
    do i = kvec0,nkvec,nprock
      iv = iv + 1
      idk = indk(i)
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
  do i = 1,nlocalkvec
    Gmin = 2.0_dp*rradmax
    do j = i,nlocalkvec
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
  if (lsg1.or.lgrad2) then
    call gstrterms(ndim,maxkvec,nlocalkvec,xrk,yrk,zrk,dg2ds,d2g2dx2,d2g2ds2,lgrad2)
    if (lrigid) then
      call gxyzterms(ndim,maxkvec,nlocalkvec,xrk,yrk,zrk,dgds,d2gds2,lgrad2)
    endif
  endif
  do ii = 1,ncharge
    i = nchargeptr(ii)
    oci = occuf(i)*angstoev
    qli = qf(i)*oci
    xci = xclat(i)
    yci = yclat(i)
    zci = zclat(i)
    nmi = natmol(i)
    nregioni = nregionno(nsft+i)
    lopi = (.not.lfreeze.or.lopf(i))
    lslicei = lsliceatom(nsft + i)
!
    if (lrigid.and.nmi.gt.0.and.nat(i).le.maxele) then
      xcomi = molxyz(1,natinmol(i),nmi)
      ycomi = molxyz(2,natinmol(i),nmi)
      zcomi = molxyz(3,natinmol(i),nmi)
    else
      xcomi = 0.0_dp
      ycomi = 0.0_dp
      zcomi = 0.0_dp
    endif
!
    if (lopi) then
      ix = 3*(noptatrptr(i) - 1) + 1
      iy = ix + 1
      iz = ix + 2
    endif
    jloop: do jj = 1,ii
      j = nchargeptr(jj)
      ocj = occuf(j)
      if (i.eq.j) then
        fct = 0.5_dp
      else 
        fct = 1.0_dp
      endif
      qlj = qf(j)*ocj*fct
      lopj = (.not.lfreeze.or.lopf(j))
      nmj = natmol(j)
      nregionj = nregionno(nsft+j)
!
!  If freezing is turned on then skip i and j pairs that don't change
!
      if (.not.lopi.and..not.lopj) cycle jloop
!
!  Set region 2 pair flag
!
      nregionj = nregionno(nsft+j)
      lreg2one  = .false.
      lreg2pair = .false.
      if (lseok.and.nregions(ncf).ge.2) then
        lreg2pair = (nregioni.eq.2.and.nregionj.eq.2)
        if (.not.lreg2pair) lreg2one = (nregioni.eq.2.or.nregionj.eq.2)
      endif
      lslicej = lsliceatom(nsft + j)
      lattach = (lslicei.and..not.lslicej.or.lslicej.and..not.lslicei)
!
      if (lrigid) then
        if (nmj.gt.0.and.nat(j).le.maxele) then
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
      if (lopj) then
        jx = 3*(noptatrptr(j) - 1) + 1
        jy = jx + 1
        jz = jx + 2
      endif
      if (i.ne.j) then
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
        etaz = seta*zd
        etaz2 = etaz*etaz
!
!  First term - K vector independent
!
        derfez = g_derf(etaz)
        dexpz  = exp(-etaz2)
        twoqv  = qfct*vol4pi
        eztrm  = twoqv*(zd*derfez + dexpz*rpieta)
        if (procid.eq.0) then
          if (lreg2one) then
            esregion12 = esregion12 - eztrm
          elseif (lreg2pair) then
            esregion2 = esregion2 - eztrm
          else
            erecip = erecip - eztrm
          endif
          if (lattach) eattach = eattach - eztrm
!
          eregion2region(nregionj,nregioni) = eregion2region(nregionj,nregioni) - eztrm
!
          siteenergy(i) = siteenergy(i) - 0.5_dp*eztrm
          siteenergy(j) = siteenergy(j) - 0.5_dp*eztrm
        endif
        if (lgrad1) then
          dtrm1 = - twoqv*derfez
          if (lgrad2) then
            dtrm2 = - twoqv*tweatpi*dexpz
          endif
        endif
!
!  Add in charge derivatives
!
        if (lgrad1.and.lDoQDeriv1) then
          d0term = - vol4pi*(zd*derfez + dexpz*rpieta)*fct*angstoev
          call d1charge(i,j,lopi,lopj,1_i4,d0term*qf(j),d0term*qf(i))
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
          if (lgrad2) then
            d21q = 0.0_dp
            d22q = 0.0_dp
            d23q = 0.0_dp
            d24q = 0.0_dp
            d25q = 0.0_dp
            d26q = 0.0_dp
          endif
          if (Gmax.ge.smallestG) then
            iv = 1
            kvec = kmod(iv)
            do while (iv.le.nlocalkvec.and.kvec.le.Gmax)
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
              if (lgrad2) then
!
!  Second derivatives with respect to atoms
!
                d21q = d21q + cosq*d2g2dx2(iv,1)*kexperfc
                d22q = d22q + cosq*d2g2dx2(iv,2)*kexperfc
                d23q = d23q - cosq*(kvec*kvec*kexperfc - 2.0_dp*tweatpi*(kvec*(dexp1*dexp3+dexp2*dexp4) &
                  - seta*(darg1*dexp1*dexp3+darg2*dexp2*dexp4)))
                d24q = d24q + yrkk*sineq*ztrm1
                d25q = d25q + xrkk*sineq*ztrm1
                d26q = d26q + cosq*d2g2dx2(iv,3)*kexperfc
              endif
              if (lsg1) then
!
!  Strain first derivatives
!
                rkvec = 1.0_dp/kvec
                strm1 = rkvec*(-rkvec*kexperfc + zd*(dexp1*derfc1-dexp2*derfc2) - rpieta*(dexp1*dexp3+dexp2*dexp4))
                do kk = 1,nstrains
                  strderv(kk) = strderv(kk) + strm1*cosq*dg2ds(iv,kk)
                enddo
!
                if (lrigid) then
                  if (latomicstress) then
                    do kk = 1,nstrains
                      dGrds = xrkk*drxyzds(kk,1) + yrkk*drxyzds(kk,2) + &
                              xd*dgds(iv,1,kk) + yd*dgds(iv,2,kk) 
                      strderv(kk) = strderv(kk) - sineq*dGrds*kexperfc
                      atomicstress(kk,i) = atomicstress(kk,i) - 0.5_dp*sineq*dGrds*kexperfc
                      atomicstress(kk,j) = atomicstress(kk,j) - 0.5_dp*sineq*dGrds*kexperfc
                    enddo
                  else
                    do kk = 1,nstrains
                      dGrds = xrkk*drxyzds(kk,1) + yrkk*drxyzds(kk,2) + &
                              xd*dgds(iv,1,kk) + yd*dgds(iv,2,kk)
                      strderv(kk) = strderv(kk) - sineq*dGrds*kexperfc
                    enddo
                  endif
                endif
!
                if (latomicstress) then
                  do kk = 1,nstrains
                    atomicstress(kk,i) = atomicstress(kk,i) + 0.5_dp*strm1*cosq*dg2ds(iv,kk)
                    atomicstress(kk,j) = atomicstress(kk,j) + 0.5_dp*strm1*cosq*dg2ds(iv,kk)
                  enddo
                  if (lfinitestrain) then
                    atomicstress(1,i) = atomicstress(1,i) - 0.5_dp*strainddetds(1)*straindet*eztrm
                    atomicstress(2,i) = atomicstress(2,i) - 0.5_dp*strainddetds(2)*straindet*eztrm
                    atomicstress(3,i) = atomicstress(3,i) - 0.5_dp*strainddetds(3)*straindet*eztrm
                    atomicstress(1,j) = atomicstress(1,j) - 0.5_dp*strainddetds(1)*straindet*eztrm
                    atomicstress(2,j) = atomicstress(2,j) - 0.5_dp*strainddetds(2)*straindet*eztrm
                    atomicstress(3,j) = atomicstress(3,j) - 0.5_dp*strainddetds(3)*straindet*eztrm
                  else
                    atomicstress(1,i) = atomicstress(1,i) - 0.5_dp*eztrm
                    atomicstress(2,i) = atomicstress(2,i) - 0.5_dp*eztrm
                    atomicstress(1,j) = atomicstress(1,j) - 0.5_dp*eztrm
                    atomicstress(2,j) = atomicstress(2,j) - 0.5_dp*eztrm
                  endif
                endif
!
                if (lgrad2) then
                  eltrm1 = - 2.0_dp*strm1*cosq*rkvec*rkvec
                  eltrm2 = cosq*rkvec*rkvec*(dexp1*derfc1*(rkvec*(rkvec-zd)+zd*zd) + dexp2*derfc2*(rkvec*(rkvec+zd)+zd*zd) + &
                    rpieta*(dexp1*dexp3*(rkvec-zd+0.5*kvec/eta) + dexp2*dexp4*(rkvec+zd+0.5*kvec/eta)))
                  eltrm = (eltrm1 + eltrm2)
                  sztrm2 = cosq*rkvec*(zd*kvec*kexperfc + tweatpi*rkvec*(dexp1*dexp3 - dexp2*dexp4))
                  do kk = 1,nstrains
                    if (lrigid) then
                      dGrds1 = xrkk*drxyzds(kk,1) + yrkk*drxyzds(kk,2) + xd*dgds(iv,1,kk) + yd*dgds(iv,2,kk)
                    endif
                    do ll = 1,nstrains
                      if (lrigid) then
                        dGrds2 = xrkk*drxyzds(ll,1) + yrkk*drxyzds(ll,2) + xd*dgds(iv,1,ll) + yd*dgds(iv,2,ll)
                        d2Grds12 = dgds(iv,1,ll)*drxyzds(kk,1) + dgds(iv,2,ll)*drxyzds(kk,2) + &
                                   dgds(iv,1,kk)*drxyzds(ll,1) + dgds(iv,2,kk)*drxyzds(ll,2) + &
                                   xrkk*d2rxyzds2(ll,kk,1) + yrkk*d2rxyzds2(ll,kk,2) + &
                                   xd*d2gds2(iv,1,kk,ll) + yd*d2gds2(iv,2,kk,ll)
!
                        sderv2(kk,ll) = sderv2(kk,ll) - sineq*strm1*(dg2ds(iv,kk)*dGrds2 + dg2ds(iv,ll)*dGrds1)
                        sderv2(kk,ll) = sderv2(kk,ll) - cosq*kexperfc*dGrds1*dGrds2
                        sderv2(kk,ll) = sderv2(kk,ll) - sineq*kexperfc*d2Grds12
                      endif
                      sderv2(kk,ll) = sderv2(kk,ll) + eltrm*dg2ds(iv,kk)*dg2ds(iv,ll) + strm1*cosq*d2g2ds2(iv,ll,kk)
                    enddo
!
                    d3kxy = dg2ds(iv,kk)*strm1*sineq
                    d3kz  = dg2ds(iv,kk)*sztrm2
!
                    if (lopi) then
                      derv3(ix,kk) = derv3(ix,kk) + xrkk*d3kxy
                      derv3(iy,kk) = derv3(iy,kk) + yrkk*d3kxy
                      derv3(iz,kk) = derv3(iz,kk) - d3kz
                    endif
                    if (lopj) then
                      derv3(jx,kk) = derv3(jx,kk) - xrkk*d3kxy
                      derv3(jy,kk) = derv3(jy,kk) - yrkk*d3kxy
                      derv3(jz,kk) = derv3(jz,kk) + d3kz
                    endif
!
                    if (lrigid) then
                      derv3(ix,kk) = derv3(ix,kk) + xrkk*cosq*kexperfc*dGrds1
                      derv3(iy,kk) = derv3(iy,kk) + yrkk*cosq*kexperfc*dGrds1
                      derv3(iz,kk) = derv3(iz,kk) + sineq*ztrm1*dGrds1
                      derv3(jx,kk) = derv3(jx,kk) - xrkk*cosq*kexperfc*dGrds1
                      derv3(jy,kk) = derv3(jy,kk) - yrkk*cosq*kexperfc*dGrds1
                      derv3(jz,kk) = derv3(jz,kk) - sineq*ztrm1*dGrds1
                    endif
                  enddo
                endif
              endif
              iv = iv + 1
              kvec = kmod(iv)
            enddo
          endif
        else
          if (Gmax.ge.smallestG) then
            iv = 1
            kvec = kmod(iv)
            do while (iv.le.nlocalkvec.and.kvec.le.Gmax)
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
!***************
!  i equals j  *
!***************
        qfct = qli*qlj
!
!  First term - K vector independent
!
        twoqv = qfct*vol4pi
        if (procid.eq.0) then
          if (lreg2one) then
            esregion12 = esregion12 - twoqv*rpieta
          elseif (lreg2pair) then
            esregion2 = esregion2 - twoqv*rpieta
          else
            erecip = erecip - twoqv*rpieta
          endif
          if (lattach) eattach = eattach - twoqv*rpieta
!
          eregion2region(nregionj,nregioni) = eregion2region(nregionj,nregioni) - twoqv*rpieta
!
          siteenergy(i) = siteenergy(i) - 0.5_dp*twoqv*rpieta
          siteenergy(j) = siteenergy(j) - 0.5_dp*twoqv*rpieta
        endif
        if (lgrad2) then
          dtrm2 = - twoqv*tweatpi
        endif
!
!  Add in charge derivatives
!
        if (lgrad1.and.lDoQDeriv1) then
          d0term = - vol4pi*rpieta*fct*angstoev
          call d1charge(i,j,lopi,lopj,1_i4,d0term*qf(j),d0term*qf(i))
        endif
!
!  Second term - K vector dependent
!
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
          do iv = 1,nlocalkvec
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
            if (lgrad2) then
!
!  Second derivatives with respect to atoms
!
              d21q = d21q + cosq*d2g2dx2(iv,1)*kexperfc
              d22q = d22q + cosq*d2g2dx2(iv,2)*kexperfc
              d23q = d23q - cosq*(kvec*kvec*kexperfc - 4.0_dp*tweatpi*(kvec*dexp1 - seta*darg1*dexp1))
              d26q = d26q + cosq*d2g2dx2(iv,3)*kexperfc
            endif
            if (lsg1) then
!
!  Strain first derivatives
!
              rkvec = 1.0_dp/kvec
              strm1 = rkvec*(-rkvec*kexperfc - 2.0_dp*rpieta*dexp1)
              do kk = 1,nstrains
                strderv(kk) = strderv(kk) + strm1*cosq*dg2ds(iv,kk)
              enddo
              if (latomicstress) then
                do kk = 1,nstrains
                  atomicstress(kk,i) = atomicstress(kk,i) + 0.5_dp*strm1*cosq*dg2ds(iv,kk)
                  atomicstress(kk,j) = atomicstress(kk,j) + 0.5_dp*strm1*cosq*dg2ds(iv,kk)
                enddo
                if (lfinitestrain) then
                  atomicstress(1,i) = atomicstress(1,i) - 0.5_dp*strainddetds(1)*straindet*twoqv*rpieta
                  atomicstress(2,i) = atomicstress(2,i) - 0.5_dp*strainddetds(2)*straindet*twoqv*rpieta
                  atomicstress(3,i) = atomicstress(3,i) - 0.5_dp*strainddetds(3)*straindet*twoqv*rpieta
                  atomicstress(1,j) = atomicstress(1,j) - 0.5_dp*strainddetds(1)*straindet*twoqv*rpieta
                  atomicstress(2,j) = atomicstress(2,j) - 0.5_dp*strainddetds(2)*straindet*twoqv*rpieta
                  atomicstress(3,j) = atomicstress(3,j) - 0.5_dp*strainddetds(3)*straindet*twoqv*rpieta
                else
                  atomicstress(1,i) = atomicstress(1,i) - 0.5_dp*twoqv*rpieta
                  atomicstress(2,i) = atomicstress(2,i) - 0.5_dp*twoqv*rpieta
                  atomicstress(1,j) = atomicstress(1,j) - 0.5_dp*twoqv*rpieta
                  atomicstress(2,j) = atomicstress(2,j) - 0.5_dp*twoqv*rpieta
                endif
              endif
              if (lgrad2) then
                eltrm1 = - 2.0_dp*strm1*cosq*rkvec*rkvec
                eltrm2 = cosq*rkvec*rkvec*(2.0_dp*derfc1*(rkvec*rkvec) + rpieta*(2.0_dp*dexp1*(rkvec+0.5_dp*kvec/eta)))
                eltrm = (eltrm1 + eltrm2)
                do kk = 1,nstrains
                  do ll = 1,nstrains
                    sderv2(kk,ll) = sderv2(kk,ll) + eltrm*dg2ds(iv,kk)*dg2ds(iv,ll) + strm1*cosq*d2g2ds2(iv,ll,kk)
                  enddo
                enddo
              endif
            endif
          enddo
        else
          do iv = 1,nlocalkvec
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
      if (lreg2one) then
        esregion12 = esregion12 + csinq*qfct
      elseif (lreg2pair) then
        esregion2 = esregion2 + csinq*qfct
      else
        erecip = erecip + csinq*qfct
      endif
      if (lattach) eattach = eattach + csinq*qfct
!
      eregion2region(nregionj,nregioni) = eregion2region(nregionj,nregioni) + csinq*qfct
!
      siteenergy(i) = siteenergy(i) + 0.5_dp*csinq*qfct
      siteenergy(j) = siteenergy(j) + 0.5_dp*csinq*qfct
!
!  Add in charge derivatives
!
      if (lgrad1.and.lDoQDeriv1) then
        d0term = csinq*fct*angstoev
        call d1charge(i,j,lopi,lopj,1_i4,d0term*qf(j),d0term*qf(i))
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
        if (nregioni.ne.nregionj) then
          xregdrv(nregioni) = xregdrv(nregioni) - sinqx
          xregdrv(nregionj) = xregdrv(nregionj) + sinqx
          yregdrv(nregioni) = yregdrv(nregioni) - sinqy
          yregdrv(nregionj) = yregdrv(nregionj) + sinqy
        endif
!
!  Because dtrm1 doesn't depend on kvector only do on 1 processor
!
        if (procid.eq.0) then
          zdrv(i) = zdrv(i) - sinqz - dtrm1
          zdrv(j) = zdrv(j) + sinqz + dtrm1
          if (nregioni.ne.nregionj) then
            zregdrv(nregioni) = zregdrv(nregioni) - sinqz - dtrm1
            zregdrv(nregionj) = zregdrv(nregionj) + sinqz + dtrm1
          endif
        else
          zdrv(i) = zdrv(i) - sinqz
          zdrv(j) = zdrv(j) + sinqz
          if (nregioni.ne.nregionj) then
            zregdrv(nregioni) = zregdrv(nregioni) - sinqz
            zregdrv(nregionj) = zregdrv(nregionj) + sinqz
          endif
        endif
        if (lgrad2) then
          if (lopi.and.lopj) then
            ixc = ix
            iyc = iy
            izc = iz
            jxc = jx
            jyc = jy
            jzc = jz
          elseif (lopi) then
            ixc = ix
            iyc = iy
            izc = iz
            jxc = ix
            jyc = iy
            jzc = iz
          elseif (lopj) then
            ixc = jx
            iyc = jy         
            izc = jz
            jxc = jx
            jyc = jy
            jzc = jz
          endif
          derv2(jxc,ixc) = derv2(jxc,ixc) + d21q
          derv2(jyc,ixc) = derv2(jyc,ixc) + d26q
          derv2(jzc,ixc) = derv2(jzc,ixc) + d25q
          derv2(jxc,iyc) = derv2(jxc,iyc) + d26q
          derv2(jyc,iyc) = derv2(jyc,iyc) + d22q
          derv2(jzc,iyc) = derv2(jzc,iyc) + d24q
          derv2(jxc,izc) = derv2(jxc,izc) + d25q
          derv2(jyc,izc) = derv2(jyc,izc) + d24q
          derv2(jzc,izc) = derv2(jzc,izc) + d23q - dtrm2
        endif
      endif
    enddo jloop
  enddo
!****************
!  Global sums  *
!****************
  if (.not.lgrad2) then
    tsum0 = g_cpu_time()
    call sumone(erecip,sum,"recip2D","erecip")
    esum = sum(1)
    if (lsg1) then
      call sumall(strderv,sum,nstrains,"recip2D","strderv")
      do i = 1,nstrains
        strderv(i) = sum(i)
      enddo
    endif
    tsum = tsum + g_cpu_time() - tsum0
  elseif (lsg1) then
    esum = erecip
  endif
  if (lDoQDeriv2.and.lgrad2) then
    allocate(phsq(nlocalkvec),stat=status)
    if (status/=0) call outofmemory('recip2D','phsq')
    allocate(ktrm3(nlocalkvec),stat=status)
    if (status/=0) call outofmemory('recip2D','ktrm3')
    allocate(ktrm3z(nlocalkvec),stat=status)
    if (status/=0) call outofmemory('recip2D','ktrm3z')
    allocate(ktrm6(nlocalkvec),stat=status)
    if (status/=0) call outofmemory('recip2D','ktrm6')
    allocate(ktrm62(nlocalkvec),stat=status)
    if (status/=0) call outofmemory('recip2D','ktrm62')
    allocate(ktrm63(nlocalkvec),stat=status)
    if (status/=0) call outofmemory('recip2D','ktrm63')
!*********************************************************************************
!  Calculation of variable charge derivative contribution to second derivatives  *
!*********************************************************************************
!
!  To save space :
!  d2i2 is stored in csin
!  d2ij is stored in argc 
!  d2j2 is stored in sine
!
    call gstrterms(ndim,maxkvec,nlocalkvec,xrk,yrk,zrk,dg2ds,d2g2dx2,d2g2ds2,lgrad2)
    do iv = 1,nlocalkvec
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
      lopi = (.not.lfreeze.or.lopf(i))
      if (lopi) then
        ix = ix + 3
        iy = iy + 3
        iz = iz + 3
      endif
      jx = - 2
      jy = - 1
      jz =   0
      do j = 1,i
        ocj = occuf(j)
        if (i.eq.j) then
          fct = 0.5_dp*angstoev
        else
          fct = angstoev
        endif
        qlj = qf(j)*ocj
        lopj = (.not.lfreeze.or.lopf(j))
        if (lopj) then
          jx = jx + 3
          jy = jy + 3
          jz = jz + 3
        endif
!
!  Find relative vector between atoms
!
        xd = xclat(j) - xci
        yd = yclat(j) - yci
        zd = zclat(j) - zci
        etaz = seta*zd
        etaz2 = etaz*etaz
!
!  First term - K vector independent
!
        derfez = g_derf(etaz)
        dexpz  = exp(-etaz2)
        etrm   = - vol4pi*(zd*derfez + dexpz*rpieta)*fct
        dtrm1  = - vol4pi*derfez*fct
! d2E/dqi.dqj
        d2trm1dij = etrm*oci*ocj
! d2E/dqi.dz
        d2trm1diz = dtrm1*oci*qlj
! d2E/dqj.dz
        d2trm1djz = dtrm1*ocj*qli
! dE/dqi
        dtrm1di = etrm*oci*qlj
! dE/dqj
        dtrm1dj = etrm*qli*ocj
!
!  Find local kvector cut-off
!
        if (abs(etaz).gt.argtest) then
          Gmax = abs(accf2/zd)
        else
          Gmax = sqrt(4.0_dp*eta*(accf2-etaz2))
        endif
        if (Gmax.ge.smallestG) then
          do iv = 1,nlocalkvec
            kvec = kmod(iv)
            if (kvec.le.Gmax) then
              xrkk = xrk(iv)
              yrkk = yrk(iv)
              arg = xrkk*xd + yrkk*yd
              cosa = cos(arg)*fct*ktrm(iv)
              sina = sin(arg)*fct*ktrm(iv)
              dexp1 = exp(kvec*zd)
              dexp2 = 1.0_dp/dexp1
              darg1 = kvec*rhseta + etaz
              darg2 = kvec*rhseta - etaz
              dexp3 = exp(-(darg1)**2)
              dexp4 = exp(-(darg2)**2)
              derfc1 = g_derfc(darg1)
              derfc2 = g_derfc(darg2)
              kexperfc = dexp1*derfc1 + dexp2*derfc2
              ztrm1 = (kvec*(dexp1*derfc1-dexp2*derfc2) - tweatpi*(dexp1*dexp3-dexp2*dexp4))
              rkvec = 1.0_dp/kvec
              strm1 = rkvec*(-rkvec*kexperfc + zd*(dexp1*derfc1-dexp2*derfc2) - rpieta*(dexp1*dexp3+dexp2*dexp4))
! d2E/dqi.dqj
              argc(iv) = cosa*oci*ocj*kexperfc
! d2E/dqi.dr
              ktrm3(iv) = - sina*kexperfc
! d2E/dqi.dz
              ktrm3z(iv) = cosa*ztrm1
! d2E/dqi.de
              ktrm6(iv) = - oci*qlj*cosa*strm1
! d2E/dqj.de
              ktrm62(iv) = - ocj*qli*cosa*strm1
! dE/dqi
              phsq(iv) = cosa*oci*qlj*kexperfc
! dE/dqj
              ktrm63(iv) = cosa*qli*ocj*kexperfc
            else
              argc(iv) = 0.0_dp
              ktrm3(iv) = 0.0_dp
              ktrm3z(iv) = 0.0_dp
              ktrm6(iv) = 0.0_dp
              ktrm62(iv) = 0.0_dp
              phsq(iv) = 0.0_dp
              ktrm63(iv) = 0.0_dp
            endif
          enddo
        else
          do iv = 1,nlocalkvec
            argc(iv) = 0.0_dp
            ktrm3(iv) = 0.0_dp
            ktrm3z(iv) = 0.0_dp
            ktrm6(iv) = 0.0_dp
            ktrm62(iv) = 0.0_dp
            phsq(iv) = 0.0_dp
            ktrm63(iv) = 0.0_dp
          enddo
        endif
!
!  Call d2charge
!
        d2self = 0.0_dp
        d1ix = 0.0_dp
        d1iy = 0.0_dp
        d1iz = 0.0_dp
        do iv = 1,nlocalkvec
          d1ix = d1ix + ktrm3(iv)*xrk(iv)
          d1iy = d1iy + ktrm3(iv)*yrk(iv)
          d1iz = d1iz + ktrm3z(iv)
        enddo
!
        d1jx = d1ix*qli*ocj
        d1jy = d1iy*qli*ocj
        d1jz = d1iz*qli*ocj
        d1ix = d1ix*qlj*oci
        d1iy = d1iy*qlj*oci
        d1iz = d1iz*qlj*oci
        if (lstr) then
          do kk = 1,nstrains
            d1is(kk) = 0.0_dp
            d1js(kk) = 0.0_dp
            do iv = 1,nlocalkvec
              d1is(kk) = d1is(kk) - ktrm6(iv)*dg2ds(iv,kk)
              d1js(kk) = d1js(kk) - ktrm62(iv)*dg2ds(iv,kk)
            enddo
          enddo
        endif
        call d2charge2D(i,j,nlocalkvec,ix,iy,iz,jx,jy,jz,lopi,lopj,phsq,ktrm63, &
                        d1ix,d1iy,d1iz,d1jx,d1jy,d1jz,d1is,d1js,csin,argc,sine, &
                        d2self,d2trm1dij,d2trm1diz,d2trm1djz,dtrm1di,dtrm1dj)
      enddo
    enddo
    deallocate(ktrm63,stat=status)
    if (status/=0) call deallocate_error('recip2D','ktrm63')
    deallocate(ktrm62,stat=status)
    if (status/=0) call deallocate_error('recip2D','ktrm62')
    deallocate(ktrm6,stat=status)
    if (status/=0) call deallocate_error('recip2D','ktrm6')
    deallocate(ktrm3z,stat=status)
    if (status/=0) call deallocate_error('recip2D','ktrm3z')
    deallocate(ktrm3,stat=status)
    if (status/=0) call deallocate_error('recip2D','ktrm3')
    deallocate(phsq,stat=status)
    if (status/=0) call deallocate_error('recip2D','phsq')
  endif
!****************************************************
!  Complete strain derivatives in reciprocal space  *
!****************************************************
  if (lsg1) then
    if (lgrad2) then
!
!  Area corrections to mixed second derivatives
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
            do j = 1,3
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
          endif
        enddo
      endif
!
!  Area corrections to strain second derivatives
!
      if (lfinitestrain) then
        do i = 1,3
          do j = 1,3
            sderv2(j,i) = sderv2(j,i) - strderv(j)*strainddetds(i)*straindet &
                                      - strderv(i)*strainddetds(j)*straindet &
                                      + 2.0_dp*esum*strainddetds(j)*strainddetds(i)*straindet**2 &
                                      - esum*straind2detds2(j,i)*straindet
          enddo
        enddo
      else
        do i = 1,2
          do j = 1,2
            sderv2(j,i) = sderv2(j,i) - strderv(j) - strderv(i)
            sderv2(j,i) = sderv2(j,i) + esum
          enddo
          sderv2(3,i) = sderv2(3,i) - strderv(3)
          sderv2(i,3) = sderv2(i,3) - strderv(3)
        enddo
      endif
    endif
    if (lfinitestrain) then
      strderv(1) = strderv(1) - strainddetds(1)*straindet*esum
      strderv(2) = strderv(2) - strainddetds(2)*straindet*esum
      strderv(3) = strderv(3) - strainddetds(3)*straindet*esum
    else
      strderv(1) = strderv(1) - esum
      strderv(2) = strderv(2) - esum
    endif
  endif
  if (lpolar) then
!*******************************
!  Electric field calculation  *
!*******************************
    do i = 1,numat
      do jj = 1,ncharge
        j = nchargeptr(jj)
        xd = xclat(j) - xclat(i)
        yd = yclat(j) - yclat(i)
        zd = zclat(j) - zclat(i)
        qlj = qf(j)*occuf(j)*angstoev
!
        etaz = seta*zd
        etaz2 = etaz * etaz
        derfez = g_derf(etaz)
        dexpz  = exp(-etaz2)
        twoqv  = vol4pi
        dtrm1 = - twoqv*derfez
        vz(i) = vz(i) - dtrm1*qlj
!
        if (abs(etaz).gt.argtest) then
          Gmax = abs(accf2/zd)
        else
          Gmax = sqrt(4.0_dp*eta*(accf2-etaz2))
        endif
!
        if (Gmax.ge.smallestG) then
          do iv = 1,nlocalkvec
            kvec = kmod(iv)
            if (kvec.le.Gmax) then
              arg = xrk(iv)*xd + yrk(iv)*yd
              cosa = cos(arg)*ktrm(iv)*qlj
              sina = sin(arg)*ktrm(iv)*qlj
              dexp1 = exp(kvec*zd)
              dexp2 = 1.0_dp/dexp1
              darg1 = kvec*rhseta + etaz
              darg2 = kvec*rhseta - etaz
              derfc1 = g_derfc(darg1)
              derfc2 = g_derfc(darg2)
              kexperfc = dexp1*derfc1 + dexp2*derfc2
              dexp3 = exp(-(darg1)**2)
              dexp4 = exp(-(darg2)**2)
              ztrm1 = (kvec*(dexp1*derfc1-dexp2*derfc2) - tweatpi*(dexp1*dexp3-dexp2*dexp4))
              vx(i) = vx(i) + sina*kexperfc*xrk(iv)
              vy(i) = vy(i) + sina*kexperfc*yrk(iv)
              vz(i) = vz(i) - cosa*ztrm1
              if (lattach) then
                vx12(i) = vx12(i) + sina*kexperfc*xrk(iv)
                vy12(i) = vy12(i) + sina*kexperfc*yrk(iv)
                vz12(i) = vz12(i) - cosa*ztrm1
              endif
            endif
          enddo
        endif
      enddo
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
  if (status/=0) call deallocate_error('recip2D','sum')
!
!  Timing
!
  time1 = g_cpu_time()
  tion = tion + time1 - time0
#ifdef TRACE
  call trace_out('recip2D')
#endif
!
  return
  end
