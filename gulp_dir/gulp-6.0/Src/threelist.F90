  subroutine threelist(ethb,lgrad1)
!
!  Subroutine for three-body energy for MD using list.
!
!  Not available for Exponential and Axilrod-Teller as these
!  are normally intermolecular three-body potentials.
!
!   1/95 Intra/intermolecular specification added
!   1/95 K3 added for three-body potential
!   2/95 SW three-body potentials added
!   3/95 Bcross potential added
!   3/96 Urey-Bradley potential added
!   4/96 Exponentially decaying Vessal form added
!   6/96 General theta0 added to SW3
!   9/97 Bug in sw3 derivatives fixed
!   3/98 Cosine-harmonic form added
!   4/98 Small constant added to rtrm1/rtrm2 in SW3 to avoid overflow
!   4/98 Error in derivatives for distance dependent potentials when
!        theta = 180 corrected
!   5/98 Potential dependent parts placed in separated subroutine
!   6/98 Murrell-Mottram potential added
!  10/98 BAcross potential added
!  10/98 Conversion of theta to rad now done in here
!   3/99 Parallel modifications introduced
!   8/99 Linear-three potential added
!   5/01 Modifications added for rhos in sw3
!   6/01 Passing of cut-offs to threebody fixed
!  10/02 Bcoscross potential added
!  11/02 Parallel changes made
!   6/04 Sign of virial corrected
!   9/04 New arguments added to threebody
!  10/05 Hydrogen-bond potential added
!  12/05 ESFF equatorial potential added
!   9/06 Theta tapering added
!   1/07 UFF3 potential added
!  12/07 Unused variables removed & handling of rho3 corrected
!  11/08 BAcoscross form added
!  11/08 Option to output energy terms added
!   3/08 3coulomb potential added
!   6/09 Site energy and virials added
!   6/09 Module name changed from three to m_three
!   7/09 Modifications for exp2 potential added
!   5/10 g3coulomb potential added
!  10/11 Strain derivatives added
!  11/11 Region-region energy contributions stored
!   4/12 Explicit virial calculation removed as no longer needed
!   4/12 xvir, yvir and zvir removed
!   5/12 Atomic stresses added
!   2/13 BAGcross potential added
!   6/13 BALcross potential added
!   1/15 k3 and k4 terms added for cosine form
!   1/15 Modified Stillinger-Weber 3-body added
!   3/15 MM3angle added
!   8/15 Garofalini form of sw3 added
!   2/18 Trace added
!   9/18 Strain module added
!   9/18 Call to threestrterms modified to allow for second derivatives
!   9/18 Call to threestrterms replaced with more general realstrterms
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!  11/19 ppp3body added
!   4/20 Rigid molecule modifications added
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
  use configurations, only : nregionno
  use g_constants
  use control,        only : lDoQDeriv1, latomicstress, lrigid
  use current
  use derivatives
  use energies,       only : siteenergy, eregion2region
  use iochannels,     only : ioout
  use m_strain,       only : realstrterms
  use m_three
  use molecule
  use numbers,        only : third
  use parallel
  use species,        only : spinspec
  use symmetry,       only : lstr
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  real(dp),    intent(inout)                :: ethb
  logical,     intent(in)                   :: lgrad1
!
!  Local variables
!
  integer(i4)                               :: i
  integer(i4)                               :: ii
  integer(i4)                               :: ix
  integer(i4)                               :: iy
  integer(i4)                               :: iz
  integer(i4)                               :: j
  integer(i4)                               :: jj
  integer(i4)                               :: jx
  integer(i4)                               :: jy
  integer(i4)                               :: jz
  integer(i4)                               :: k
  integer(i4)                               :: kl
  integer(i4)                               :: ks
  integer(i4)                               :: m
  integer(i4)                               :: n
  integer(i4)                               :: n3ty
  integer(i4)                               :: nj
  integer(i4)                               :: nlast
  integer(i4)                               :: nmi
  integer(i4)                               :: nmj
  integer(i4)                               :: nmk
  integer(i4)                               :: nregioni
  integer(i4)                               :: nregionj
  integer(i4)                               :: nregionk
  integer(i4)                               :: nspeci
  integer(i4)                               :: nspecj
  integer(i4)                               :: nspeck
  integer(i4)                               :: nt2
  integer(i4)                               :: ntyp2
  integer(i4)                               :: ntypj
  logical                                   :: lsg1
  real(dp)                                  :: ang
  real(dp)                                  :: g_cpu_time
  real(dp)                                  :: d0i
  real(dp)                                  :: d0j
  real(dp)                                  :: d0k
  real(dp)                                  :: d1q(3,3)
  real(dp)                                  :: d2q(6)
  real(dp)                                  :: dot
  real(dp)                                  :: e2d(1)
  real(dp)                                  :: e3d(1)
  real(dp)                                  :: ed11
  real(dp)                                  :: ed12
  real(dp)                                  :: ed13
  real(dp)                                  :: ethb1
  real(dp)                                  :: oci
  real(dp)                                  :: ocj
  real(dp)                                  :: ock
  real(dp)                                  :: ofct
  real(dp)                                  :: one
  real(dp)                                  :: qli
  real(dp)                                  :: qlj
  real(dp)                                  :: qlk
  real(dp)                                  :: r12
  real(dp)                                  :: r122
  real(dp)                                  :: r13
  real(dp)                                  :: r132
  real(dp)                                  :: r23
  real(dp)                                  :: r232
  real(dp)                                  :: rho1
  real(dp)                                  :: rho2
  real(dp)                                  :: rho3
  real(dp)                                  :: rho4
  real(dp)                                  :: rho5
  real(dp)                                  :: rk32
  real(dp)                                  :: rk33
  real(dp)                                  :: rk34
  real(dp)                                  :: rkthb
  real(dp)                                  :: rkthb3
  real(dp)                                  :: rkthb4
  real(dp)                                  :: rktmp
  real(dp)                                  :: rstrdloc(6)
  real(dp)                                  :: ro1
  real(dp)                                  :: ro2
  real(dp)                                  :: ro3
  real(dp)                                  :: ro4
  real(dp)                                  :: ro5
  real(dp)                                  :: dr2ds(6,3)
  real(dp)                                  :: d2r2dx2(3,3,3)
  real(dp)                                  :: d2r2ds2(6,6,3)
  real(dp)                                  :: d2r2dsdx(6,3,3)
  real(dp)                                  :: the0
  real(dp)                                  :: time1
  real(dp)                                  :: time2
  real(dp)                                  :: tr11
  real(dp)                                  :: tr21
  real(dp)                                  :: tr31
  real(dp)                                  :: ttr11
  real(dp)                                  :: ttr21
  real(dp)                                  :: x21
  real(dp)                                  :: y21
  real(dp)                                  :: z21
  real(dp)                                  :: x23
  real(dp)                                  :: y23
  real(dp)                                  :: z23
  real(dp)                                  :: x31
  real(dp)                                  :: y31
  real(dp)                                  :: z31
  real(dp)                                  :: xc1
  real(dp)                                  :: yc1
  real(dp)                                  :: zc1
  real(dp)                                  :: xcom(3)
  real(dp)                                  :: ycom(3)
  real(dp)                                  :: zcom(3)
  real(dp)                                  :: xcomi
  real(dp)                                  :: ycomi
  real(dp)                                  :: zcomi
  real(dp)                                  :: xv3(3)
  real(dp)                                  :: yv3(3)
  real(dp)                                  :: zv3(3)
!
  if (nlist3md.eq.0) return
#ifdef TRACE
  call trace_in('threelist')
#endif
  time1 = g_cpu_time()
  lsg1 = (lstr.and.lgrad1)
!
!  Opening banner for energy decomposition
!
  if (lPrintThree) then
    call mpbarrier
    if (ioproc) then
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''  Three: Atom No. 1  Atom No. 2  Atom No. 3              Threebody energy (eV)  '')')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
    endif
  endif
!
!  Initialisation
!
  ethb = 0.0_dp
  one = 1.0_dp
  nlast = 0
!
!  Loop over valid terms
!
  do m = procid+1,nlist3md,nprocs
    n = nthbptr(m)
    i = i3ind(m)
    j = j3ind(m)
    k = k3ind(m)
    ii = icell31(m)
    if (ii.eq.0) then
      ix = 0
      iy = 0
      iz = 0
    else
      iz = (ii/100) - 5
      ii = ii - 100*(iz+5)
      iy = (ii/10) - 5
      ii = ii - 10*(iy+5)
      ix = ii - 5
    endif
    jj = icell32(m)
    if (jj.eq.0) then
      jx = 0
      jy = 0
      jz = 0
    else
      jz = (jj/100) - 5
      jj = jj-100*(jz+5)
      jy = (jj/10) - 5
      jj = jj - 10*(jy+5)
      jx = jj - 5
    endif
    if (n.ne.nlast) then
      n3ty = nthrty(n)
      nt2 = ntspec2(n)
      ntyp2 = ntptyp2(n)
      tr11 = thr1(n)
      tr21 = thr2(n)
      tr31 = thr3(n)
      rkthb = thbk(n)
      rkthb3 = 0.0_dp
      rkthb4 = 0.0_dp
      ro1  =  0.0_dp
      ro2  =  0.0_dp
      ro3  =  0.0_dp
      ro4  =  0.0_dp
      ro5  =  0.0_dp
      if (n3ty.eq.2) then
        the0 = theta(n)*degtorad
        ro1 = thrho1(n)
        ro2 = thrho2(n)
        if (ro1.ne.0.0_dp) ro1 = 1.0_dp/ro1
        if (ro2.ne.0.0_dp) ro2 = 1.0_dp/ro2
      elseif (n3ty.eq.1) then
        the0 = theta(n)*degtorad
        rkthb3 = thrho2(n)/6.0_dp
        rkthb4 = thrho1(n)/24.0_dp
      elseif (n3ty.eq.5.or.n3ty.eq.25) then
        the0 = cos(theta(n)*degtorad)
        ro1 = thrho1(n)
        ro2 = thrho2(n)
        ro3 = thrho3(n)
      elseif (n3ty.eq.6) then
        the0 = 0.0_dp
        ro1 = thrho1(n)
        ro2 = thrho2(n)
      elseif (n3ty.eq.7) then
        the0 = theta(n)
      elseif (n3ty.eq.8) then
        the0 = theta(n)*degtorad
        ro1 = thrho1(n)
        ro2 = thrho2(n)
        if (ro1.ne.0.0_dp) ro1 = 1.0_dp/ro1
        if (ro2.ne.0.0_dp) ro2 = 1.0_dp/ro2
        the0 = the0 - pi
        the0 = the0*the0
        rkthb = 0.25_dp*rkthb/the0
      elseif (n3ty.eq.9) then
        the0 = cos(theta(n)*degtorad)
        rkthb3 = thrho2(n)/6.0_dp
        rkthb4 = thrho1(n)/24.0_dp
      elseif (n3ty.eq.10) then
        the0 = theta(n)
        ro1 = thrho1(n)
        ro2 = thrho2(n)
        ro3 = thrho3(n)
      elseif (n3ty.eq.11) then
        rkthb3 = thrho3(n)
        ro1 = thrho1(n)
        ro2 = thrho2(n)
        the0 = theta(n)*degtorad
      elseif (n3ty.eq.12) then
        the0 = theta(n)
        rkthb3 = nint(thrho1(n))
      elseif (n3ty.eq.13) then
        the0 = theta(n)
        rkthb3 = thrho3(n)
        ro1 = thrho1(n)
        ro2 = thrho2(n)  
      elseif (n3ty.eq.14) then
        the0 = cos(theta(n)*degtorad)
        ro1 = thrho1(n)
        ro2 = thrho2(n)
        ro3 = thrho3(n)
      elseif (n3ty.eq.15) then
        rkthb3 = theta(n)
        ro1 = thrho1(n)
        ro2 = thrho2(n)
        ro3 = thrho3(n)
      elseif (n3ty.eq.16) then
        the0 = theta(n)
        ro1 = thrho1(n)
        ro2 = thrho2(n)
      elseif (n3ty.eq.17) then
        the0 = theta(n)*degtorad
        rkthb4 = 1.0_dp/(2.0_dp*sin(the0))**2
        rkthb3 = - 4.0_dp*rkthb4*cos(the0)
        the0 = rkthb4*(2.0_dp*cos(the0)**2 + 1.0_dp)
      elseif (n3ty.eq.18) then
        rkthb3 = thrho3(n)
        ro1 = thrho1(n)
        ro2 = thrho2(n)
        the0 = cos(theta(n)*degtorad)
      elseif (n3ty.eq.20) then
        the0 = 0.0_dp
        ro1 = theta(n)
        ro2 = thrho2(n)
        ro4 = thrho1(n)
        ro5 = thrho3(n)
      elseif (n3ty.eq.21) then
        the0 = theta(n)
      elseif (n3ty.eq.22) then
        ro1 = thrho1(n)
        ro2 = thrho2(n)
        the0 = theta(n)*degtorad
      elseif (n3ty.eq.23) then
        ro1 = thrho1(n)
        ro2 = thrho2(n)
        the0 = theta(n)
      elseif (n3ty.eq.24) then
        the0 = theta(n)*degtorad
      elseif (n3ty.eq.27) then
        the0 = 0.0_dp
        ro1 = thrho1(n)
        ro2 = thrho2(n)
        ro3 = thrho3(n)
      endif
      nlast = n
    endif
!
!  First atom
!
    nregioni = nregionno(nsft+i)
    nspeci = nspecptr(nrelf2a(i))
    oci = occuf(i)
    qli = qf(i)
    xc1 = xclat(i)
    yc1 = yclat(i)
    zc1 = zclat(i)
!
!  Second atom
!
    nj = nat(j)
    ntypj = nftype(j)
    nregionj = nregionno(nsft+j)
    nspecj = nspecptr(nrelf2a(j))
    ocj = occuf(j)
    qlj = qf(j)
    if (nj.eq.nt2.and.(ntyp2.eq.0.or.ntyp2.eq.ntypj)) then
      rho1 = ro1
      rho2 = ro2
      rho3 = ro3
      rho4 = ro4
      rho5 = ro5
      ttr11 = tr11
      ttr21 = tr21
    else
      rho1 = ro2
      rho2 = ro1
      rho3 = ro3
      rho4 = ro5
      rho5 = ro4
      ttr11 = tr21
      ttr21 = tr11
      if (n3ty.eq.11.or.n3ty.eq.18) then
        rktmp = rkthb
        rkthb = rkthb3
        rkthb3 = rktmp
      endif
    endif
!
!  Rigid molecule handling
!
    if (lrigid) then
      nmi = natmol(i)
!
!  Set COM coordinates
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
      nmj = natmol(j)
      if (nmj.gt.0) then
        xcom(1) = molxyz(1,natinmol(j),nmj) - xcomi
        ycom(1) = molxyz(2,natinmol(j),nmj) - ycomi
        zcom(1) = molxyz(3,natinmol(j),nmj) - zcomi
      else
        xcom(1) = - xcomi
        ycom(1) = - ycomi
        zcom(1) = - zcomi
      endif
!
      nmk = natmol(k)
      if (nmk.gt.0) then
        xcom(2) = molxyz(1,natinmol(k),nmk) - xcomi
        ycom(2) = molxyz(2,natinmol(k),nmk) - ycomi
        zcom(2) = molxyz(3,natinmol(k),nmk) - zcomi
      else
        xcom(2) = - xcomi
        ycom(2) = - ycomi
        zcom(2) = - zcomi
      endif
      xcom(3) = xcom(2) - xcom(1)
      ycom(3) = ycom(2) - ycom(1)
      zcom(3) = zcom(2) - zcom(1)
    else
      xcom(1:3) = 0.0_dp
      ycom(1:3) = 0.0_dp
      zcom(1:3) = 0.0_dp
    endif

    x21 = xclat(j) - xc1 + ix*r1x + iy*r2x + iz*r3x
    y21 = yclat(j) - yc1 + ix*r1y + iy*r2y + iz*r3y
    z21 = zclat(j) - zc1 + ix*r1z + iy*r2z + iz*r3z
    r122 = x21*x21 + y21*y21 + z21*z21
    r12 = sqrt(r122)
!
!  Third atom
!
    nregionk = nregionno(nsft+k)
    nspeck = nspecptr(nrelf2a(k))
    ock = occuf(k)
    qlk = qf(k)
    x31 = xclat(k) - xc1 + jx*r1x + jy*r2x + jz*r3x
    y31 = yclat(k) - yc1 + jx*r1y + jy*r2y + jz*r3y
    z31 = zclat(k) - zc1 + jx*r1z + jy*r2z + jz*r3z
    r132 = x31*x31 + y31*y31 + z31*z31
    r13 = sqrt(r132)
    x23 = x31 - x21
    y23 = y31 - y21
    z23 = z31 - z21
!
!  No need to sqrt r23 as only needed as r23**2
!
    r232 = x23**2 + y23**2 + z23**2
    if (n3ty.ne.6.and.n3ty.ne.7) then
      dot = x21*x31 + y21*y31 + z21*z31
      dot = dot/(r12*r13)
      if (abs(dot).gt.0.999999999999_dp) dot = sign(one,dot)
      if (n3ty.eq.9) then
        ang = dot
      else
        ang = acos(dot)
      endif
    else
      dot = 0.0_dp
      ang = 0.0_dp
    endif
    r23 = sqrt(r232)
    ofct = oci*ocj*ock
    if (n3ty.eq.19) then
      rk32 = rkthb*ofct*qf(j)*qf(k)
    else
      rk32 = rkthb*ofct
    endif
    if (n3ty.eq.12.or.n3ty.eq.17) then
      rk33 = rkthb3
    elseif (n3ty.eq.26) then
      rk33 = spinspec(nspeci)*spinspec(nspecj)*spinspec(nspeck)
    else
      rk33 = rkthb3*ofct
    endif
    if (n3ty.eq.17) then
      rk34 = rkthb4
    else
      rk34 = rkthb4*ofct
    endif
    if (n3ty.eq.15) then
      rho1 = thrho1(n)
      rho2 = thrho2(n)
      rho3 = thrho3(n)
    elseif (n3ty.eq.16) then
      rho1 = thrho1(n)
      rho2 = thrho2(n)
    endif
!*****************************************************
!  Calculate derivatives with respect to potentials  *
!*****************************************************
    call threebody(1_i4,n3ty,r12,r13,r23,ed11,ed12,ed13,ethb1,e2d,e3d,ttr11,ttr21,tr31, &
                   rho1,rho2,rho3,rho4,rho5,rk32,rk33,rk34,the0,ang,dot,.true.,.false.,.false., &
                   n,qli,qlj,qlk,d0i,d0j,d0k,d1q,d2q,lthetataper(n),thetatapermin(n), &
                   thetatapermax(n))
    ethb = ethb + ethb1
!
    eregion2region(nregionj,nregioni) = eregion2region(nregionj,nregioni) + third*ethb1
    eregion2region(nregionk,nregioni) = eregion2region(nregionk,nregioni) + third*ethb1
    eregion2region(nregionk,nregionj) = eregion2region(nregionk,nregionj) + third*ethb1
!
    siteenergy(i) = siteenergy(i) + third*ethb1
    siteenergy(j) = siteenergy(j) + third*ethb1
    siteenergy(k) = siteenergy(k) + third*ethb1
!
!  Output energy contribution
!
    if (lPrintThree) then
      write(ioout,'(4x,3i12,8x,f27.10)') i,j,k,ethb1
    endif
!***********************
!  Strain derivatives  *
!***********************
    if (lsg1) then
!
!  Set up strain products
!
      xv3(1) = x21
      xv3(2) = x31
      xv3(3) = x23
      yv3(1) = y21
      yv3(2) = y31
      yv3(3) = y23
      zv3(1) = z21
      zv3(2) = z31
      zv3(3) = z23
      call realstrterms(ndim,3_i4,3_i4,xv3,yv3,zv3,xcom,ycom,zcom,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
!
!  First strain derivatives
!
      rstrdloc(1:nstrains) = 0.0_dp
      do kl = 1,nstrains
        ks = nstrptr(kl)
        rstrdloc(kl) = rstrdloc(kl) + ed11*dr2ds(ks,1)
        rstrdloc(kl) = rstrdloc(kl) + ed12*dr2ds(ks,2)
        rstrdloc(kl) = rstrdloc(kl) + ed13*dr2ds(ks,3)
        rstrd(kl) = rstrd(kl) + rstrdloc(kl)
      enddo
      if (latomicstress) then
        do kl = 1,nstrains
          atomicstress(kl,i) = atomicstress(kl,i) + third*rstrdloc(kl)
          atomicstress(kl,j) = atomicstress(kl,j) + third*rstrdloc(kl)
          atomicstress(kl,k) = atomicstress(kl,k) + third*rstrdloc(kl)
        enddo
      endif
    endif
    if (lgrad1) then
!*************************
!  Internal derivatives  *
!*************************
      xdrv(i) = xdrv(i) - x21*ed11 - x31*ed12
      ydrv(i) = ydrv(i) - y21*ed11 - y31*ed12
      zdrv(i) = zdrv(i) - z21*ed11 - z31*ed12
      xdrv(j) = xdrv(j) + x21*ed11 - x23*ed13
      ydrv(j) = ydrv(j) + y21*ed11 - y23*ed13
      zdrv(j) = zdrv(j) + z21*ed11 - z23*ed13
      xdrv(k) = xdrv(k) + x31*ed12 + x23*ed13
      ydrv(k) = ydrv(k) + y31*ed12 + y23*ed13
      zdrv(k) = zdrv(k) + z31*ed12 + z23*ed13
!   
!  Charge derivatives
!   
      if (lDoQDeriv1) then
        call d1charge3(i,j,k,.true.,.true.,.true.,1_i4,d0i,d0j,d0k)
      endif
    endif
  enddo
!
!  Closing banner for energy decomposition
!
  if (lPrintThree) then
    call mpbarrier
    if (ioproc) then
      write(ioout,'(''--------------------------------------------------------------------------------'',/)')
    endif
  endif
!
!  Timing
!
  time2 = g_cpu_time()
  tthree = tthree + time2 - time1
#ifdef TRACE
  call trace_out('threelist')
#endif
!
  return
  end
