  subroutine recip3Dspme(erecip,lgrad1)
!
!  Calculation of electrostatic potential and derivatives,
!  including strain derivatives. Reciprocal space part.
!  Smooth particle mesh Ewald (SPME) version.
!
!   6/16 Created from recip3Dmd
!   8/16 First derivatives added and debugged
!   2/18 Trace added
!   9/18 Handling of lstraincell algorithm added
!   9/18 Strain module introduced
!  11/18 Finite strain flag introduced instead of lstraincell
!   1/19 Unused used variables removed
!   3/20 Use of charge pointer added
!   7/20 Separate routine for sumall with 1 argument added
!
!  NB: This algorithm doesn't support multiple regions or
!      the calculation of site energies
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
  use kspace,         only : eta, rradmax, dg2ds
  use m_fft3D,        only : cq, initfft3D, fft3Dforward, fft3Dbackward
  use m_strain,       only : gstrterms, strainddetds, straindet
  use optimisation
  use parallel
  use polarise
  use spme
  use symmetry
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  logical,        intent(in)                        :: lgrad1
  real(dp),       intent(out)                       :: erecip
!
!  Local variables
!
  integer(i4),                                   save :: halfK1
  integer(i4),                                   save :: halfK2
  integer(i4),                                   save :: halfK3
  integer(i4)                                         :: i
  integer(i4)                                         :: idim
  integer(i4)                                         :: ifail
  integer(i4)                                         :: ii
  integer(i4)                                         :: iu
  integer(i4),                                   save :: K1
  integer(i4),                                   save :: K2
  integer(i4),                                   save :: K3
  integer(i4),                                   save :: K3loc
  integer(i4),                                   save :: K3start
  integer(i4)                                         :: k
  integer(i4)                                         :: knot1
  integer(i4)                                         :: knot2
  integer(i4)                                         :: knot3
  integer(i4)                                         :: littlek1
  integer(i4)                                         :: littlek2
  integer(i4)                                         :: littlek3
  integer(i4)                                         :: m1
  integer(i4)                                         :: m1lower
  integer(i4)                                         :: m1p
  integer(i4)                                         :: m2
  integer(i4)                                         :: m2p
  integer(i4)                                         :: m3
  integer(i4)                                         :: m3l
  integer(i4)                                         :: m3p
  integer(i4)                                         :: nm3
  integer(i4)                                         :: nm3loc
  integer(i4),  dimension(:),       allocatable, save :: nm3ptr
  integer(i4),                                   save :: ncflast = 0
  integer(i4)                                         :: status
  logical,                                       save :: lcelldone = .false.
  logical                                             :: lincreasesize
  logical                                             :: lsg1
  logical                                             :: lneedsetup
  complex(dpc), dimension(:),       allocatable, save :: cexp1
  complex(dpc), dimension(:),       allocatable, save :: cexp2
  complex(dpc), dimension(:),       allocatable, save :: cexp3
  complex(dpc)                                        :: crsum
  complex(dpc)                                        :: bofm1
  complex(dpc)                                        :: bofm2
  complex(dpc)                                        :: bofm3
  real(dp),     dimension(:),       allocatable, save :: b2mod1
  real(dp),     dimension(:),       allocatable, save :: b2mod2
  real(dp),     dimension(:),       allocatable, save :: b2mod3
  real(dp)                                            :: arge
  real(dp)                                            :: cq2
  real(dp),     dimension(:,:,:),   allocatable, save :: Cspme
  real(dp),     dimension(:,:,:,:), allocatable, save :: Vspme    ! Stores virial components 
  real(dp),     dimension(:,:,:),   allocatable, save :: rq
  real(dp),     dimension(:),       allocatable       :: dmnu1
  real(dp),     dimension(:),       allocatable       :: dmnu2
  real(dp),     dimension(:),       allocatable       :: dmnu3
  real(dp)                                            :: g_cpu_time
  real(dp)                                            :: dtmp
  real(dp)                                            :: duconv
  real(dp)                                            :: eta4
  real(dp)                                            :: expphasec
  real(dp)                                            :: expphases
  real(dp)                                            :: esum
  real(dp)                                            :: kspmetrm
  real(dp)                                            :: dkspmetrm
  real(dp),     dimension(:),       allocatable       :: mnknot
  real(dp),     dimension(:,:),     allocatable       :: mnu1
  real(dp),     dimension(:,:),     allocatable       :: mnu2
  real(dp),     dimension(:,:),     allocatable       :: mnu3
  real(dp)                                            :: phase
  real(dp)                                            :: qi
  real(dp)                                            :: r11
  real(dp)                                            :: r12
  real(dp)                                            :: r13
  real(dp)                                            :: r21
  real(dp)                                            :: r22
  real(dp)                                            :: r23
  real(dp)                                            :: r31
  real(dp)                                            :: r32
  real(dp)                                            :: r33
  real(dp)                                            :: reta
  real(dp)                                            :: rk2
  real(dp)                                            :: rrk2
  real(dp)                                            :: rmat(3,3)
  real(dp),                                      save :: rpivol
  real(dp)                                            :: rqx
  real(dp)                                            :: rqy
  real(dp)                                            :: rqz
  real(dp)                                            :: rrmx2
  real(dp)                                            :: strdervloc(6)
  real(dp),     dimension(:),       allocatable       :: sum
  real(dp)                                            :: time0
  real(dp)                                            :: time1
  real(dp)                                            :: tsum0
  real(dp)                                            :: twopioverK
  real(dp)                                            :: u
  real(dp)                                            :: u1i
  real(dp)                                            :: u2i
  real(dp)                                            :: u3i
  real(dp)                                            :: wrk(6)
  real(dp)                                            :: xrk
  real(dp)                                            :: yrk
  real(dp)                                            :: zrk
  real(dp)                                            :: dg2dsl(1,6)
  real(dp)                                            :: d2g2ds2l(1,6,6)
  real(dp)                                            :: d2g2dx2l(1,6)
  real(dp)                                            :: x1(1,1)
  real(dp)                                            :: y1(1,1)
  real(dp)                                            :: z1(1,1)
  real(dp)                                            :: xpon
  real(dp)                                            :: vol
  real(dp)                                            :: volume
#ifdef TRACE
  call trace_in('recip3Dspme')
#endif
!
  time0 = g_cpu_time()
!
!  Initialise energies
!
  erecip = 0.0_dp
  lsg1 = (lstr.and.lgrad1)
  if (lsg1) then
!
!  Initialise strain derivatives
!
    strdervloc(1:6) = 0.0_dp
  endif
!
!  Allocate local memory
!
  allocate(mnknot(nBsplineorder),stat=status)
  if (status/=0) call outofmemory('recip3Dmdspme','mnknot')
  allocate(mnu1(nBsplineorder,numat),stat=status)
  if (status/=0) call outofmemory('recip3Dmdspme','mnu1')
  allocate(mnu2(nBsplineorder,numat),stat=status)
  if (status/=0) call outofmemory('recip3Dmdspme','mnu2')
  allocate(mnu3(nBsplineorder,numat),stat=status)
  if (status/=0) call outofmemory('recip3Dmdspme','mnu3')
  allocate(nm3ptr(nBsplineorder),stat=status)
  if (status/=0) call outofmemory('recip3Dmdspme','nm3ptr')
  if (lgrad1) then
    allocate(dmnu1(nBsplineorder),stat=status)
    if (status/=0) call outofmemory('recip3Dmdspme','dmnu1')
    allocate(dmnu2(nBsplineorder),stat=status)
    if (status/=0) call outofmemory('recip3Dmdspme','dmnu2')
    allocate(dmnu3(nBsplineorder),stat=status)
    if (status/=0) call outofmemory('recip3Dmdspme','dmnu3')
  endif
  allocate(sum(max(numat,nstrains)),stat=status)
  if (status/=0) call outofmemory('recip3Dmdspme','sum')
!
!  Setup
!
  eta4 = 0.25_dp/eta
  reta = eta4/eta
  rrmx2 = rradmax*rradmax
!
!  End of set-up section
!
  if (lnorecip) goto 999
!
!  Only do general set up for grid if the cell changes
!
  if (ncf.ne.ncflast) lcelldone = .false.
  lneedsetup = (ncf.ne.ncflast.or.ncell.gt.0.or.(lsg1.and..not.lcelldone))
  ncflast = ncf
  if (lneedsetup) then
!*********************************************************
!  Set up that is independent of the atomic coordinates  *
!*********************************************************
!
!  Determine grid size
!
! NB for now grid size is read in and stored in nqkgrid
    K1 = nqkgrid(1,ncf)
    K2 = nqkgrid(2,ncf)
    K3 = nqkgrid(3,ncf)
!
    halfK1 = K1/2_i4
    halfK2 = K2/2_i4
    halfK3 = K3/2_i4
!
!  Check array dimensions
!
    lincreasesize = .false.
    do idim = 1,3
      if (nqkgrid(idim,ncf).gt.maxqkgrid(idim)) then
        lincreasesize = .true.
        maxqkgrid(idim) = nqkgrid(idim,ncf)
      endif
    enddo
!
!  Initialise FFT
!
    call initfft3D(K1,K2,K3,K3loc,K3start)
!
    if (lincreasesize.or..not.allocated(Cspme)) then
      allocate(Cspme(K1,K2,K3loc),stat=status)
      if (status/=0) call outofmemory('recip3Dmdspme','Cspme')
      allocate(rq(K1,K2,K3loc),stat=status)
      if (status/=0) call outofmemory('recip3Dmdspme','rq')
      allocate(b2mod1(K1),stat=status)
      if (status/=0) call outofmemory('recip3Dmdspme','b2mod1')
      allocate(b2mod2(K2),stat=status)
      if (status/=0) call outofmemory('recip3Dmdspme','b2mod2')
      allocate(b2mod3(K3),stat=status)
      if (status/=0) call outofmemory('recip3Dmdspme','b2mod3')
      allocate(cexp1(K1),stat=status)
      if (status/=0) call outofmemory('recip3Dmdspme','cexp1')
      allocate(cexp2(K2),stat=status)
      if (status/=0) call outofmemory('recip3Dmdspme','cexp2')
      allocate(cexp3(K3),stat=status)
      if (status/=0) call outofmemory('recip3Dmdspme','cexp3')
    endif
    if (lsg1) then
      if (lincreasesize.or..not.allocated(Vspme)) then
        lcelldone = .true.
        allocate(Vspme(6,K1,K2,K3loc),stat=status)
        if (status/=0) call outofmemory('recip3Dmdspme','Vspme')
      endif
    endif
!
!  Compute B spline coefficients for 1 to n - 1 at knots
!
    mnknot(1) = 0.0_dp
    mnknot(nBsplineorder) = 0.0_dp
    do iu = 1,nBsplineorder-1
      u = dble(iu)
      call bspline_M(nBsplineorder,u,mnknot(iu))
    enddo
!
!  Set up complex exponential phase factors for each cell vector
!
    twopioverK = 2.0_dp*pi/dble(K1)
    do m1 = 0,K1-1
      phase = twopioverK*dble(m1)
      expphasec = cos(phase)
      expphases = sin(phase)
      cexp1(m1+1) = dcmplx(expphasec,expphases)
    enddo
!
    twopioverK = 2.0_dp*pi/dble(K2)
    do m2 = 0,K2-1
      phase = twopioverK*dble(m2)
      expphasec = cos(phase)
      expphases = sin(phase)
      cexp2(m2+1) = dcmplx(expphasec,expphases)
    enddo
!
    twopioverK = 2.0_dp*pi/dble(K3)
    do m3 = 0,K3-1
      phase = twopioverK*dble(m3)
      expphasec = cos(phase)
      expphases = sin(phase)
      cexp3(m3+1) = dcmplx(expphasec,expphases)
    enddo
!
!  Compute b_i(m_mi) coefficients
!
    do m1 = 0,K1-1
      crsum = 0.0_dpc
      do iu = 0,nBsplineorder-2
        crsum = crsum + dcmplx(mnknot(iu+1),0.0_dp)*cexp1(mod(m1*iu,K1)+1)
      enddo
      bofm1 = cexp1(mod(m1*(nBsplineorder-1),K1)+1)/crsum
      b2mod1(m1+1) = dble(bofm1*dconjg(bofm1))
    enddo
!
    do m2 = 0,K2-1
      crsum = 0.0_dpc
      do iu = 0,nBsplineorder-2
        crsum = crsum + dcmplx(mnknot(iu+1),0.0_dp)*cexp2(mod(m2*iu,K2)+1)
      enddo
      bofm2 = cexp2(mod(m2*(nBsplineorder-1),K2)+1)/crsum
      b2mod2(m2+1) = dble(bofm2*dconjg(bofm2))
    enddo
!
    do m3 = 0,K3-1
      crsum = 0.0_dpc
      do iu = 0,nBsplineorder-2
        crsum = crsum + dcmplx(mnknot(iu+1),0.0_dp)*cexp3(mod(m3*iu,K3)+1)
      enddo
      bofm3 = cexp3(mod(m3*(nBsplineorder-1),K3)+1)/crsum
      b2mod3(m3+1) = dble(bofm3*dconjg(bofm3))
    enddo
!
!  Compute Cspme(m1,m2,m3)
!
    vol = volume(rv)
    rpivol = 2.0_dp*pi/vol
    do m3l = 0,K3loc-1
      m3 = K3start + m3l
      if (m3.le.halfK3) then
        m3p = m3
      else
        m3p = m3 - K3
      endif
      do m2 = 0,K2-1
        if (m2.le.halfK2) then
          m2p = m2
        else
          m2p = m2 - K2
        endif
!
!  Set lower bound of m1 to avoid m = 0 term
!
        if (m1.eq.0.and.m2.eq.0) then
          m1lower = 1
        else
          m1lower = 0
        endif
        do m1 = m1lower,K1-1
          if (m1.le.halfK1) then
            m1p = m1
          else
            m1p = m1 - K1
          endif
!
          xrk = m1p*kv(1,1) + m2p*kv(1,2) + m3p*kv(1,3)
          yrk = m1p*kv(2,1) + m2p*kv(2,2) + m3p*kv(2,3)
          zrk = m1p*kv(3,1) + m2p*kv(3,2) + m3p*kv(3,3)
!
          rk2 = xrk*xrk + yrk*yrk + zrk*zrk
          arge = - rk2*eta4
          xpon = exp(arge)
          rrk2 = 1.0_dp/rk2
          kspmetrm = xpon*rrk2*b2mod1(m1+1)*b2mod2(m2+1)*b2mod3(m3+1)
          Cspme(m1+1,m2+1,m3l+1) = kspmetrm
          if (lsg1) then
            dkspmetrm = - 2.0_dp*(eta4 + rrk2)*kspmetrm
!
            x1(1,1) = xrk 
            y1(1,1) = yrk 
            z1(1,1) = zrk 
            call gstrterms(ndim,1_i4,1_i4,x1,y1,z1,dg2dsl,d2g2dx2l,d2g2ds2l,.false.)
!
            Vspme(1,m1+1,m2+1,m3l+1) = dg2ds(1,1)*dkspmetrm
            Vspme(2,m1+1,m2+1,m3l+1) = dg2ds(1,2)*dkspmetrm
            Vspme(3,m1+1,m2+1,m3l+1) = dg2ds(1,3)*dkspmetrm
            Vspme(4,m1+1,m2+1,m3l+1) = dg2ds(1,4)*dkspmetrm
            Vspme(5,m1+1,m2+1,m3l+1) = dg2ds(1,5)*dkspmetrm
            Vspme(6,m1+1,m2+1,m3l+1) = dg2ds(1,6)*dkspmetrm
          endif
        enddo
      enddo
    enddo
!
!  Set Cspme for m = 0 case
!
    if (K3start.eq.0.and.K3loc.ge.1) then
      Cspme(1,1,1) = 0.0_dp
      if (lsg1) then
        Vspme(1:6,1,1,1) = 0.0_dp
      endif
    endif
  endif

!*************************************************
!  Calculate the charge matrix in real space, Q  *
!*************************************************
  rq = 0.0_dp
!
!  Loop over atoms
!
  do ii = 1,ncharge
    i = nchargeptr(ii)
    qi = qf(i)*occuf(i)
    u1i = xfrac(i)*dble(K1)
    u2i = yfrac(i)*dble(K2)
    u3i = zfrac(i)*dble(K3)
!
!  Find lowest knot point relevant to fractional coordinate
!
    knot1 = int(u1i)
    knot2 = int(u2i)
    knot3 = int(u3i)
!
    if (mod(nBsplineorder,2).eq.0.0) then
!
!  Even order -> odd number of knots
!
      if ((u1i-dble(knot1)).gt.0.5_dp) knot1 = knot1 + 1
      if ((u2i-dble(knot2)).gt.0.5_dp) knot2 = knot2 + 1
      if ((u3i-dble(knot3)).gt.0.5_dp) knot3 = knot3 + 1
!
      knot1 = knot1 - nBsplineorder - 1
      knot2 = knot2 - nBsplineorder - 1
      knot3 = knot3 - nBsplineorder - 1
    else
!
!  Odd order -> even number of knots
!
      knot1 = knot1 - nBsplineorder - 1
      knot2 = knot2 - nBsplineorder - 1
      knot3 = knot3 - nBsplineorder - 1
    endif
    if (nprocs.gt.1) then
!
!  Find if any of the values of m3 are on this node
!
      nm3loc = 0
      do m3 = 1,nBsplineorder
        littlek3 = knot3 + m3
        if (littlek3.ge.K3) littlek3 = littlek3 - K3
        if (littlek3.lt.0)  littlek3 = littlek3 + K3
!
!  Check that element is local in parallel
!
        m3l = littleK3 + 1 - K3start
        if (m3l.gt.0.and.m3l.le.K3loc) then
          nm3loc = nm3loc + 1
          nm3ptr(nm3loc) = m3
        endif
      enddo
    else
      nm3loc = 0
      do m3 = 1,nBsplineorder
        nm3loc = nm3loc + 1
        nm3ptr(nm3loc) = m3
      enddo
    endif
!
!  Only go further if there are local elements
!
    if (nm3loc.gt.0) then
!
!  Compute B-splines for atom
!
      do m1 = 1,nBsplineorder
        littlek1 = knot1 + m1
        u = u1i - dble(littlek1)
        call bspline_M(nBsplineorder,u,mnu1(m1,i))
      enddo
      do m2 = 1,nBsplineorder
        littlek2 = knot2 + m2
        u = u2i - dble(littlek2)
        call bspline_M(nBsplineorder,u,mnu2(m2,i))
      enddo
      do m3 = 1,nBsplineorder
        littlek3 = knot3 + m3
        u = u3i - dble(littlek3)
        call bspline_M(nBsplineorder,u,mnu3(m3,i))
      enddo
!
!  Loop over cells / k values
!
      do nm3 = 1,nm3loc
        m3 = nm3ptr(nm3)
        littlek3 = knot3 + m3
        if (littlek3.ge.K3) littlek3 = littlek3 - K3
        if (littlek3.lt.0)  littlek3 = littlek3 + K3
!
!  Check that element is local in parallel
!
        m3l = littleK3 + 1 - K3start
!
        do m2 = 1,nBsplineorder
          littlek2 = knot2 + m2
          if (littlek2.ge.K2) littlek2 = littlek2 - K2
          if (littlek2.lt.0)  littlek2 = littlek2 + K2
!
          do m1 = 1,nBsplineorder
            littlek1 = knot1 + m1
            if (littlek1.ge.K1) littlek1 = littlek1 - K1
            if (littlek1.lt.0)  littlek1 = littlek1 + K1
!
!  Combine values for Q
!
            rq(littlek1+1,littlek2+1,m3l) = rq(littlek1+1,littlek2+1,m3l) + &
              qi*mnu1(m1,i)*mnu2(m2,i)*mnu3(m3,i)
          enddo
!
        enddo
!
      enddo
!
!  End of check on nm3loc
!
    endif
!
  enddo
!
!  Generate complex version of rq for FFT
!
  do m3 = 1,K3loc
    do m2 = 1,K2
      do m1 = 1,K1
        cq(m1,m2,m3) = dcmplx(rq(m1,m2,m3),0.0_dp)
      enddo
    enddo
  enddo
!************************************************
!  Use FFTs to transform Q to reciprocal space  *
!************************************************
  call fft3Dforward
!****************************************
!  Compute the reciprocal space energy  *
!****************************************
  do m3 = 1,K3loc
    do m2 = 1,K2
      do m1 = 1,K1
        erecip = erecip + Cspme(m1,m2,m3)*dble(cq(m1,m2,m3)*conjg(cq(m1,m2,m3)))
      enddo
    enddo
  enddo
!
!  Multiply by unit factors and constants to get to eV
!
  erecip = erecip*rpivol*angstoev
  if (lgrad1) then
!************************
!  Compute derivatives  *
!************************
!
!  Strain derivatives
!
    if (lsg1) then
!
!  Use cq to compute the strain derivative components first
!
      do m3 = 1,K3loc
        do m2 = 1,K2
          do m1 = 1,K1
            cq2 = dble(cq(m1,m2,m3)*conjg(cq(m1,m2,m3)))
            do k = 1,6
              strdervloc(k) = strdervloc(k) + Vspme(k,m1,m2,m3)*cq2
            enddo
          enddo
        enddo
      enddo
!
!  Convert units of strain derivatives
!
      do k = 1,6
        strdervloc(k) = strdervloc(k)*rpivol*angstoev
      enddo
    endif
!
!  Create complex copy of Cspme in cq
!
    do m3 = 1,K3loc
      do m2 = 1,K2
        do m1 = 1,K1
          cq(m1,m2,m3) = Cspme(m1,m2,m3)*cq(m1,m2,m3)
        enddo
      enddo
    enddo
!
!  Transform back to real space
!
    call fft3Dbackward
!
!  Generate real version of cq after FFT
!
    do m3 = 1,K3loc
      do m2 = 1,K2
        do m1 = 1,K1
          rq(m1,m2,m3) = real(cq(m1,m2,m3))
        enddo
      enddo
    enddo
!
!  Unit conversion factor for derivatives
!
    duconv = rpivol*angstoev*2.0d0
!
!  Invert cell matrix for use in derivatives
!
    rmat(1,1) = r1x
    rmat(2,1) = r1y
    rmat(3,1) = r1z
    rmat(1,2) = r2x
    rmat(2,2) = r2y
    rmat(3,2) = r2z
    rmat(1,3) = r3x
    rmat(2,3) = r3y
    rmat(3,3) = r3z
    call matrix_inversion(rmat,3_i4,3_i4,wrk,ifail)
!
!  Loop over atoms
!
    do ii = 1,ncharge
      i = nchargeptr(ii)
      qi = qf(i)*occuf(i)*duconv
      u1i = xfrac(i)*dble(K1)
      u2i = yfrac(i)*dble(K2)
      u3i = zfrac(i)*dble(K3)
!
      r11 = rmat(1,1)*qi*dble(K1)
      r12 = rmat(1,2)*qi*dble(K1)
      r13 = rmat(1,3)*qi*dble(K1)
      r21 = rmat(2,1)*qi*dble(K2)
      r22 = rmat(2,2)*qi*dble(K2)
      r23 = rmat(2,3)*qi*dble(K2)
      r31 = rmat(3,1)*qi*dble(K3)
      r32 = rmat(3,2)*qi*dble(K3)
      r33 = rmat(3,3)*qi*dble(K3)
!
!  Find lowest knot point relevant to fractional coordinate
!
      knot1 = int(u1i)
      knot2 = int(u2i)
      knot3 = int(u3i)
!
      if (mod(nBsplineorder,2).eq.0.0) then
!
!  Even order -> odd number of knots
!
        if ((u1i-dble(knot1)).gt.0.5_dp) knot1 = knot1 + 1
        if ((u2i-dble(knot2)).gt.0.5_dp) knot2 = knot2 + 1
        if ((u3i-dble(knot3)).gt.0.5_dp) knot3 = knot3 + 1
!
        knot1 = knot1 - nBsplineorder - 1
        knot2 = knot2 - nBsplineorder - 1
        knot3 = knot3 - nBsplineorder - 1
      else
!
!  Odd order -> even number of knots
!
        knot1 = knot1 - nBsplineorder - 1
        knot2 = knot2 - nBsplineorder - 1
        knot3 = knot3 - nBsplineorder - 1
      endif
      if (nprocs.gt.1) then
!
!  Find if any of the values of m3 are on this node
!
        nm3loc = 0
        do m3 = 1,nBsplineorder
          littlek3 = knot3 + m3
          if (littlek3.ge.K3) littlek3 = littlek3 - K3
          if (littlek3.lt.0)  littlek3 = littlek3 + K3
!
!  Check that element is local in parallel
!
          m3l = littleK3 + 1 - K3start
          if (m3l.gt.0.and.m3l.le.K3loc) then
            nm3loc = nm3loc + 1
            nm3ptr(nm3loc) = m3
          endif
        enddo
      else
        nm3loc = 0
        do m3 = 1,nBsplineorder
          nm3loc = nm3loc + 1
          nm3ptr(nm3loc) = m3
        enddo
      endif
!
!  Only go further if there are local elements
!
      if (nm3loc.gt.0) then
!
!  Compute first derivatives of B-splines for atom
!
!  Compute M_(n-1)(u) and M_(n-1)(u-1) & then take the difference
!
        do m1 = 1,nBsplineorder
          littlek1 = knot1 + m1
          u = u1i - dble(littlek1)
          call bspline_M(nBsplineorder-1,u,dmnu1(m1))
          u = u - 1.0_dp
          call bspline_M(nBsplineorder-1,u,dtmp)
          dmnu1(m1) = dmnu1(m1) - dtmp
        enddo
        do m2 = 1,nBsplineorder
          littlek2 = knot2 + m2
          u = u2i - dble(littlek2)
          call bspline_M(nBsplineorder-1,u,dmnu2(m2))
          u = u - 1.0_dp
          call bspline_M(nBsplineorder-1,u,dtmp)
          dmnu2(m2) = dmnu2(m2) - dtmp
        enddo
        do m3 = 1,nBsplineorder
          littlek3 = knot3 + m3
          u = u3i - dble(littlek3)
          call bspline_M(nBsplineorder-1,u,dmnu3(m3))
          u = u - 1.0_dp
          call bspline_M(nBsplineorder-1,u,dtmp)
          dmnu3(m3) = dmnu3(m3) - dtmp
        enddo
!
!  Loop over cells / k values
!
        do nm3 = 1,nm3loc
          m3 = nm3ptr(nm3)
          littlek3 = knot3 + m3
          if (littlek3.ge.K3) littlek3 = littlek3 - K3
          if (littlek3.lt.0)  littlek3 = littlek3 + K3
!
!  Check that element is local in parallel
!
          m3l = littleK3 + 1 - K3start
!
          do m2 = 1,nBsplineorder
            littlek2 = knot2 + m2
            if (littlek2.ge.K2) littlek2 = littlek2 - K2
            if (littlek2.lt.0)  littlek2 = littlek2 + K2
!
            do m1 = 1,nBsplineorder
              littlek1 = knot1 + m1
              if (littlek1.ge.K1) littlek1 = littlek1 - K1
              if (littlek1.lt.0)  littlek1 = littlek1 + K1
!
              rqx = rq(littlek1+1,littlek2+1,m3l)*dmnu1(m1)*mnu2(m2,i)*mnu3(m3,i)
              rqy = rq(littlek1+1,littlek2+1,m3l)*mnu1(m1,i)*dmnu2(m2)*mnu3(m3,i)
              rqz = rq(littlek1+1,littlek2+1,m3l)*mnu1(m1,i)*mnu2(m2,i)*dmnu3(m3)
!
!  Combine values for Q
!
              xdrv(i) = xdrv(i) + rqx*r11 + rqy*r21 + rqz*r31
              ydrv(i) = ydrv(i) + rqx*r12 + rqy*r22 + rqz*r32
              zdrv(i) = zdrv(i) + rqx*r13 + rqy*r23 + rqz*r33
            enddo
          enddo
        enddo
!
!  End of check on nm3loc
!
      endif
    enddo
  endif
!****************
!  Global sums  *
!****************
  if (nprocs.gt.1) then
    tsum0 = g_cpu_time()
    call sumone(erecip,esum,"recip3Dspme","erecip")
    if (lsg1) then
      call sumall(strdervloc,sum,6_i4,"recip3Dspme","strderv")
      do i = 1,6
        strderv(i) = strderv(i) + sum(i)
      enddo
    endif
    tsum = tsum + g_cpu_time() - tsum0
  else
    esum = erecip
    do i = 1,6
      strderv(i) = strderv(i) + strdervloc(i)
    enddo
  endif
!************************************************************************************
!  Subtract average sum of all gradients to reduce any rippling effect of the grid  *
!************************************************************************************
  if (lgrad1) then
    sum(1) = 0.0_dp
    sum(2) = 0.0_dp
    sum(3) = 0.0_dp
    do i = 1,numat
      sum(1) = sum(1) + xdrv(i)
      sum(2) = sum(2) + ydrv(i)
      sum(3) = sum(3) + zdrv(i)
    enddo
    sum(1:3) = sum(1:3)/dble(numat)
    do i = 1,numat
      xdrv(i) = xdrv(i) - sum(1)
      ydrv(i) = ydrv(i) - sum(2)
      zdrv(i) = zdrv(i) - sum(3)
    enddo
  endif
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
!***************
!  Exit point  *
!***************
999 continue
!
!  Free local memory
!
  deallocate(sum,stat=status)
  if (status/=0) call deallocate_error('recip3Dmdspme','sum')
  if (lgrad1) then
    deallocate(dmnu3,stat=status)
    if (status/=0) call deallocate_error('recip3Dmdspme','dmnu3')
    deallocate(dmnu2,stat=status)
    if (status/=0) call deallocate_error('recip3Dmdspme','dmnu2')
    deallocate(dmnu1,stat=status)
    if (status/=0) call deallocate_error('recip3Dmdspme','dmnu1')
  endif
  deallocate(nm3ptr,stat=status)
  if (status/=0) call deallocate_error('recip3Dmdspme','nm3ptr')
  deallocate(mnu3,stat=status)
  if (status/=0) call deallocate_error('recip3Dmdspme','mnu3')
  deallocate(mnu2,stat=status)
  if (status/=0) call deallocate_error('recip3Dmdspme','mnu2')
  deallocate(mnu1,stat=status)
  if (status/=0) call deallocate_error('recip3Dmdspme','mnu1')
  deallocate(mnknot,stat=status)
  if (status/=0) call deallocate_error('recip3Dmdspme','mnknot')
!
!  Timing
!
  time1 = g_cpu_time()
  tion = tion + time1 - time0
#ifdef TRACE
  call trace_out('recip3Dspme')
#endif
!
  return
  end
