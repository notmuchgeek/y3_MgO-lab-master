  subroutine real1Dfd(matom,erealin,lgrad1)
!
!  This subroutine calculates the electrostatic energy of a 1-D system 
!  in real space based on the algorithm implemented in CRYSTAL. A 
!  neutralising uniform background charge density is applied and then 
!  subtracted again. Because a sum over neutral unit cells is required, 
!  this code is kept separate from the other real space routines. The 
!  conventional real space routines are used to handle the Coulomb
!  subtraction issues in order to keep this routine simple.
!  Finite difference version that focuses on derivatives of matom.
!
!  12/17 Created from real1Dmd
!   2/18 Trace added
!   9/18 Partially modified for lstraincell algorithm
!  11/18 Finite strain flag introduced instead of lstraincell
!  12/18 Sign of first derivatives changed to match usual convention
!  12/18 real1strterm used in place of twostrterms
!  12/19 Rigid molecule modifications added
!   3/20 Location of angstoev changed to current
!   3/20 Centre of mass positions passed to hfunc/emfunc 
!   3/20 Call to emfunc change due to creation of separate emfuncs for strains
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
!  No claim is made that this program is free from errors and
!  no liability will be accepted for any loss or damage that
!  may result. The user is responsible for checking the validity
!  of their results.
!
!  Copyright Curtin University 2020
!
!  Julian Gale, CIC, Curtin University, March 2020
!
  use configurations, only : nregionno, nregiontype, QMMMmode
  use control,        only : lnoreal, lDoQDeriv1, lrigid
  use current
  use derivatives
  use element,        only : maxele
  use general,        only : nmaxcells, nemorder, smallself
  use kspace,         only : accf1D
  use m_strain,       only : real1strterm
  use molecule
  use qmedata,        only : maxloop
  use shells,         only : cuts
  use symmetry,       only : lstr
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed arguments
!
  integer(i4),                 intent(in)    :: matom
  real(dp),                    intent(inout) :: erealin
  logical,                     intent(in)    :: lgrad1
!
!  Local variables
!
  integer(i4)                                :: i
  integer(i4)                                :: j 
  integer(i4)                                :: m
  integer(i4)                                :: nati
  integer(i4)                                :: natj
  integer(i4)                                :: nmi
  integer(i4)                                :: nmj
  integer(i4)                                :: nregioni
  integer(i4)                                :: nregionj
  integer(i4)                                :: nregiontypi
  integer(i4)                                :: nregiontypj
  logical                                    :: lconverged
  logical                                    :: lcspair
  real(dp)                                   :: accf
  real(dp)                                   :: acell
  real(dp)                                   :: g_cpu_time
  real(dp)                                   :: cut2s
  real(dp)                                   :: d0
  real(dp)                                   :: d0i
  real(dp)                                   :: d0j
  real(dp)                                   :: d0term
  real(dp)                                   :: d1
  real(dp)                                   :: d1s
  real(dp)                                   :: dh1(3)
  real(dp)                                   :: dh2(3)
  real(dp)                                   :: dh1s
  real(dp)                                   :: dh2s
  real(dp)                                   :: d2h1(6)
  real(dp)                                   :: d2h2(6)
  real(dp)                                   :: d2h1m(3)
  real(dp)                                   :: d2h2m(3)
  real(dp)                                   :: d2h1s
  real(dp)                                   :: d2h2s
  real(dp)                                   :: d3h1(10)
  real(dp)                                   :: d3h1m(6)
  real(dp)                                   :: dr2ds(6)
  real(dp)                                   :: d2r2dx2(3,3)
  real(dp)                                   :: d2r2ds2(6,6)
  real(dp)                                   :: d2r2dsdx(6,3)
  real(dp)                                   :: dads
  real(dp)                                   :: ediff
  real(dp)                                   :: elast
  real(dp)                                   :: ereal
  real(dp)                                   :: esum
  real(dp)                                   :: esumem
  real(dp)                                   :: esumh
  real(dp)                                   :: esum12
  real(dp)                                   :: esumem12
  real(dp)                                   :: esumh12
  real(dp)                                   :: esum2
  real(dp)                                   :: esumem2
  real(dp)                                   :: esumh2
  real(dp)                                   :: e1
  real(dp)                                   :: e2
  real(dp)                                   :: h1
  real(dp)                                   :: h2
  real(dp)                                   :: lna
  real(dp)                                   :: oci
  real(dp)                                   :: ocj
  real(dp)                                   :: qi
  real(dp)                                   :: qj
  real(dp)                                   :: qii
  real(dp)                                   :: qij
  real(dp)                                   :: r
  real(dp)                                   :: rcut
  real(dp)                                   :: rr
  real(dp)                                   :: t1, t2
  real(dp)                                   :: u
  real(dp)                                   :: x
  real(dp)                                   :: y
  real(dp)                                   :: z
  real(dp)                                   :: xcom
  real(dp)                                   :: ycom
  real(dp)                                   :: zcom
  real(dp)                                   :: xcomi
  real(dp)                                   :: ycomi
  real(dp)                                   :: zcomi
!
!  If noreal specified, return
!
  if (lnoreal) then
    erealin = 0.0_dp
    return
  endif
#ifdef TRACE
  call trace_in('real1Dfd')
#endif
!
  t1 = g_cpu_time()
!********************************************************
!  Calculate Coulomb sum converged to desired accuracy  *
!********************************************************
!
!  Loop over number of cells in sum
!
  accf = 10.0**(-accf1D)
  lna = log(a)
  cut2s = cuts*cuts
  m = - 1
  lconverged = .false.
  elast = 0.0_dp
  esum = 0.0_dp
  esum12 = 0.0_dp
  esum2 = 0.0_dp
!
  do while (m.lt.nmaxcells.and..not.lconverged) 
    m = m + 1
!********************************************
!  Direct sum component over neutral cells  *
!********************************************
    acell = dble(m)*a
!
!  Only compute i = matom case
!
    i = matom
    nati = nat(i)
    oci = occuf(i)
    qi = qf(i)*oci
    nmi = natmol(i)
    nregioni = nregionno(nsft+i)
    nregiontypi = nregiontype(nregioni,ncf)
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
    jloop: do j = 1,numat
!
!  Exclude self term
!
      if (i.eq.j) cycle jloop
!
      natj = nat(j)
      ocj = occuf(j)
      qj = qf(j)*ocj
      nmj = natmol(j)
      nregionj = nregionno(nsft+j)
      nregiontypj = nregiontype(nregionj,ncf)
      lcspair = (abs(nati-natj).eq.maxele.or.(oci+ocj).lt.1.0001_dp)
      if (lcspair) then
        rcut = cut2s
      else
        rcut = smallself
      endif
!  
!  QM/MM handling : i & j are both QM atoms => exclude
!
      if (QMMMmode(ncf).gt.0) then
        if (nregiontypi.eq.1.and.nregiontypj.eq.1) cycle jloop
!
!  QM/MM : Electrostatic embedding : If either i or j are QM atoms => exclude 
!
        if (QMMMmode(ncf).eq.2.and.(nregiontypi.eq.1.or.nregiontypj.eq.1)) cycle jloop
      endif
!
!  Set COM coordinates
!
      if (lrigid.and.nmj.gt.0) then
        xcom = molxyz(1,natinmol(j),nmj) - xcomi
        ycom = molxyz(2,natinmol(j),nmj) - ycomi
        zcom = molxyz(3,natinmol(j),nmj) - zcomi
      else
        xcom = - xcomi
        ycom = - ycomi
        zcom = - zcomi
      endif
!
      x = acell + xclat(j) - xclat(i)
      y = yclat(j) - yclat(i)
      z = zclat(j) - zclat(i)
      r = x*x + y*y + z*z
      if (r.gt.rcut) then
        r = sqrt(r)
        rr = 1.0_dp/r
        d0 = qi*qj*rr*angstoev
        esum = esum + d0
        if (lgrad1) then
          d1 = - d0*rr*rr
          xdrv(i) = xdrv(i) - d1*x
          ydrv(i) = ydrv(i) - d1*y
          zdrv(i) = zdrv(i) - d1*z
          xdrv(j) = xdrv(j) + d1*x
          ydrv(j) = ydrv(j) + d1*y
          zdrv(j) = zdrv(j) + d1*z
          if (lstr) then
            call real1strterm(ndim,x,y,z,xcom,ycom,zcom,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
            rstrd(1) = rstrd(1) + d1*dr2ds(1)
          endif
          if (lDoQDeriv1) then
            d0i = qj*rr*angstoev
            d0j = qi*rr*angstoev
            call d1charge(i,j,.true.,.true.,1_i4,d0i,d0j)
          endif
        endif
      endif
      if (m.gt.0) then
        x = - acell + xclat(j) - xclat(i)
        r = x*x + y*y + z*z
        if (r.gt.rcut) then
          r = sqrt(r)
          rr = 1.0_dp/r
          d0 = qi*qj*rr*angstoev
          esum = esum + d0
          if (lgrad1) then
            d1 = - d0*rr*rr
            xdrv(i) = xdrv(i) - d1*x
            ydrv(i) = ydrv(i) - d1*y
            zdrv(i) = zdrv(i) - d1*z
            xdrv(j) = xdrv(j) + d1*x
            ydrv(j) = ydrv(j) + d1*y
            zdrv(j) = zdrv(j) + d1*z
            if (lstr) then
              call real1strterm(ndim,x,y,z,xcom,ycom,zcom,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
              rstrd(1) = rstrd(1) + d1*dr2ds(1)
            endif
            if (lDoQDeriv1) then
              d0i = qj*rr*angstoev
              d0j = qi*rr*angstoev
              call d1charge(i,j,.true.,.true.,1_i4,d0i,d0j)
            endif
          endif
        endif
      endif
    enddo jloop
!**********************
!  Self interactions  *
!**********************
    oci = occuf(i)
    qi = qf(i)*oci
    nregioni = nregionno(nsft+i)
    nregiontypi = nregiontype(nregioni,ncf)
!  
!  QM/MM handling : i is a QM atom => exclude
!
    if (QMMMmode(ncf).eq.0.or.nregiontypi.ne.1) then
      r = abs(acell)
      if (r.gt.smallself) then
        rr = 1.0_dp/r
        d0 = qi*qi*rr*angstoev
!
!  Set region 2 pair flags
!        
        esum = esum + d0
        if (lgrad1.and.lstr) then
          d1 = - d0*rr*rr
          call real1strterm(ndim,acell,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
          rstrd(1) = rstrd(1) + d1*dr2ds(1)
        endif
        if (lgrad1.and.lDoQDeriv1) then
          d0i = qi*rr*angstoev
          call d1charge(i,i,.true.,.true.,1_i4,d0i,d0i)
        endif
      endif
    endif
!
!  Neutralising terms
!
!  Background
!
!  and
!
!  Euler-MacLaurin component
!
    esumh = 0.0_dp
    esumh12 = 0.0_dp
    esumh2 = 0.0_dp
    esumem = 0.0_dp
    esumem12 = 0.0_dp
    esumem2 = 0.0_dp
!
    if (m.gt.0) then
      u = (dble(m)+0.5_dp)*a
      jloop2: do j = 1,numat
        ocj = occuf(j)
        qj = qf(j)*ocj
        nregionj = nregionno(nsft+j)
        nregiontypj = nregiontype(nregionj,ncf)
!  
!  QM/MM handling : i & j are both QM atoms => exclude
!
        if (QMMMmode(ncf).gt.0) then
          if (nregiontypi.eq.1.and.nregiontypj.eq.1) cycle jloop2
!
!  QM/MM : Electrostatic embedding : If either i or j are QM atoms => exclude 
!
          if (QMMMmode(ncf).eq.2.and.(nregiontypi.eq.1.or.nregiontypj.eq.1)) cycle jloop2
        endif
!
        qij = qi*qj*angstoev
        x = xclat(j) - xclat(i)
        y = yclat(j) - yclat(i)
        z = zclat(j) - zclat(i)
        call hfunc(u,+x,1.0_dp,y,z,h1,dh1,d2h1,d3h1,.false.,.false.,.false.)
        call hfunc(u,-x,-1.0_dp,y,z,h2,dh2,d2h2,d3h1,.false.,.false.,.false.)
        esumh = esumh - qij*(h1 + h2 - 2.0_dp*lna)/a
!
        call emfunc(nemorder,u,+x,1.0_dp,y,z,a,e1,dh1,d2h1,d3h1,.false.,.false.,.false.)
        call emfunc(nemorder,u,-x,-1.0_dp,y,z,a,e2,dh2,d2h2,d3h1,.false.,.false.,.false.)
        esumem = esumem + qij*(e1 + e2)
      enddo jloop2
!
!  QM/MM handling : i is a QM atom => exclude
!
      if (QMMMmode(ncf).eq.0.or.nregiontypi.ne.1) then
        qii = qi*qi*angstoev
        call hfunc(u,0.0_dp,1.0_dp,0.0_dp,0.0_dp,h1,dh1,d2h1,d3h1,.false.,.false.,.false.)
        esumh = esumh - qii*(h1 - lna)/a
!
        call emfunc(nemorder,u,0.0_dp,1.0_dp,0.0_dp,0.0_dp,a,e1,dh1,d2h1,d3h1,.false.,.false.,.false.)
        esumem = esumem + qii*e1
      endif
    endif
!
!  Sum up terms
!
    ereal = esum + esumh + esumem
!
!  Compare to energy with previous number of cells
!  and check for convergence to the required
!  accuracy.
!
    if (abs(ereal).lt.accf) lconverged = .true.
    if (.not.lconverged) then
      ediff = abs((ereal - elast)/ereal)
      lconverged = (ediff.lt.accf)
    endif
    elast = ereal
  enddo
!
!  Divide ereal by number of processors to compensate for the fact that it has already been summed
!
  erealin = erealin + ereal
!
!  Save number of cells needed
!
  maxloop(1) = m
!
!  Derivatives of Euler-MacLaurin terms since these are not cumulative
!
  if (lgrad1.and.m.gt.0) then
!
!  Set up strain derivatives if lfinitestrain
!
    if (lfinitestrain) then
      dads = a/(1.0_dp + strain(1))
    else
      dads = a
    endif
!
    u = (dble(m)+0.5_dp)*a
    jloop3: do j = 1,numat
!
!  Exclude self term 
!
      if (i.eq.j) cycle jloop3
!
      ocj = occuf(j)
      qj = qf(j)*ocj
      nregionj = nregionno(nsft+j)
      nregiontypj = nregiontype(nregionj,ncf)
!  
!  QM/MM handling : i & j are both QM atoms => exclude
!
      if (QMMMmode(ncf).gt.0) then
        if (nregiontypi.eq.1.and.nregiontypj.eq.1) cycle jloop3
!
!  QM/MM : Electrostatic embedding : If either i or j are QM atoms => exclude 
!
        if (QMMMmode(ncf).eq.2.and.(nregiontypi.eq.1.or.nregiontypj.eq.1)) cycle jloop3
      endif
      qij = qi*qj*angstoev
      x = xclat(j) - xclat(i)
      y = yclat(j) - yclat(i)
      z = zclat(j) - zclat(i)
      r = sqrt(x*x + y*y + z*z)
      call hfuncs(u,+x,1.0_dp,y,z,0.0_dp,h1,dh1,dh1s,d2h1,d2h1s,d2h1m,d3h1,lgrad1,.false.,.false.)
      call hfuncs(u,-x,-1.0_dp,y,z,0.0_dp,h2,dh2,dh2s,d2h2,d2h2s,d2h2m,d3h1,lgrad1,.false.,.false.)
!
      d1 = qij/a
      xdrv(i) = xdrv(i) + d1*(dh1(1) + dh2(1))
      ydrv(i) = ydrv(i) + d1*(dh1(2) + dh2(2))
      zdrv(i) = zdrv(i) + d1*(dh1(3) + dh2(3))
      xdrv(j) = xdrv(j) - d1*(dh1(1) + dh2(1))
      ydrv(j) = ydrv(j) - d1*(dh1(2) + dh2(2))
      zdrv(j) = zdrv(j) - d1*(dh1(3) + dh2(3))
      if (lstr) then
        d1s = (h1 + h2 - 2.0_dp*lna)*dads/a**2 - (dh1s + dh2s - 2.0_dp*dads/a)/a
        d1s = d1s*qij
        rstrd(1) = rstrd(1) + d1s
      endif
!
!  Compute variable charge term = E without charges
!
      d0term = - angstoev*(h1 + h2 - 2.0_dp*lna)/a
!
      call emfuncs(nemorder,u,+x,1.0_dp,y,z,0.0_dp,a,e1,dh1,d2h1,dh1s, &
                   d2h1s,d2h1m,d3h1,d3h1m,lgrad1,.false.,.false.)
      call emfuncs(nemorder,u,-x,-1.0_dp,y,z,0.0_dp,a,e2,dh2,d2h2,dh2s, &
                   d2h2s,d2h2m,d3h1,d3h1m,lgrad1,.false.,.false.)
!
      d1 = qij
      xdrv(i) = xdrv(i) - d1*(dh1(1) + dh2(1))
      ydrv(i) = ydrv(i) - d1*(dh1(2) + dh2(2))
      zdrv(i) = zdrv(i) - d1*(dh1(3) + dh2(3))
      xdrv(j) = xdrv(j) + d1*(dh1(1) + dh2(1))
      ydrv(j) = ydrv(j) + d1*(dh1(2) + dh2(2))
      zdrv(j) = zdrv(j) + d1*(dh1(3) + dh2(3))
      if (lstr) then
        rstrd(1) = rstrd(1) + d1*(dh1s + dh2s)
      endif
      d0term = d0term + (e1 + e2)*angstoev 
      if (lDoQDeriv1) then
        d0i = qj*d0term
        d0j = qi*d0term
        call d1charge(i,j,.true.,.true.,1_i4,d0i,d0j)
      endif
    enddo jloop3
    if (lstr.or.lDoQDeriv1) then
!***************************
!  Self-interaction terms  *
!***************************
!
!  QM/MM handling : i is a QM atom => exclude
!
      if (QMMMmode(ncf).eq.0.or.nregiontypi.ne.1) then
        qii = qi*qi*angstoev
        call hfuncs(u,0.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp,h1,dh1,dh1s,d2h1,d2h1s,d2h1m, &
                    d3h1,lgrad1,.false.,.false.)
!
        if (lstr) then
          d1s = (h1 - lna)*dads/a**2 - (dh1s - dads/a)/a
          d1s = d1s*qii
          rstrd(1) = rstrd(1) + d1s
        endif
!
!  Compute variable charge term = E without charges
!
        d0term = - angstoev*(h1 - lna)/a
        call emfuncs(nemorder,u,0.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp,a,e1,dh1,d2h1, &
                     dh1s,d2h1s,d2h1m,d3h1,d3h1m,lgrad1,.false.,.false.)
!
        if (lstr) then
          rstrd(1) = rstrd(1) + qii*dh1s*angstoev
        endif
        d0term = d0term + e1*angstoev
        if (lDoQDeriv1) then
          d0i = qi*d0term
          call d1charge(i,i,.true.,.true.,1_i4,d0i,d0i)
        endif
      endif
    endif
  endif
!
  t2 = g_cpu_time()
  tatom = tatom + t2 - t1
#ifdef TRACE
  call trace_out('real1Dfd')
#endif
!
  return
  end
