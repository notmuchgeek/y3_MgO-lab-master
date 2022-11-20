  subroutine dchargesplitd(lprint,lrecalcA,lstrainin,nbloc,nblocptr,potl,nmax)
!
!  Calculates the first derivative of the charge with respect
!  to the coordinates of the atoms. This version is specifically
!  for EEM with split bond charges.
!  Distributed memory parallel version.
!
!   5/18 Created from dchargesplit and dcharged
!   5/18 Multiple qranges added
!   7/18 dqs now made into a local temporary array
!   8/19 Strain module added
!   9/18 Dimension of dqs corrected
!  11/18 Finite strain flag introduced instead of lstraincell
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!  12/19 Rigid molecule modifications added
!
!  On entry:
!
!    potl      = inverse matrix calculated during EEM/QEq
!    nmax      = left-hand dimension of potl
!    lprint    = logical indicating whether printed output is wanted
!    lrecalcA  = if .true. this forces the recalculation of A
!    lstrainin = if .true. then strain derivatives are calculated
!                and derv3 is overwritten
!    nbloc     = number of local columns of potl
!    nblocptr  = pointer from local columns to global
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
!  Copyright Curtin University 2019
!
!  Julian Gale, CIC, Curtin University, December 2019
!
#ifdef MPI
  use configurations, only : nregionno
#endif
  use g_constants
  use control
  use current
  use derivatives
  use eembonds
  use eemdata
  use element
#ifdef MPI
  use field,          only : lfieldcfg, ntdfieldcfg
  use general,        only : cutw, etaw
#endif
  use iochannels
  use kspace
#ifdef MPI
  use m_strain,       only : strainddetds, straindet
#endif
  use m_strain,       only : gstrterms, real1strterm
  use parallel
  use qmedata
  use shells
  use symmetry
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
#ifdef MPI
  include 'mpif.h'
#endif
!
!  Passed variables
!
  integer(i4), intent(in)                      :: nmax
  integer(i4), intent(in)                      :: nbloc
  integer(i4), intent(in)                      :: nblocptr(*)
  logical,     intent(in)                      :: lprint
  logical,     intent(in)                      :: lrecalcA
  logical,     intent(in)                      :: lstrainin
  real(dp),    intent(inout)                   :: potl(nmax,nbloc)
#ifdef MPI
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ib1
  integer(i4)                                  :: ib2
  integer(i4)                                  :: ieem
  integer(i4)                                  :: ifail
  integer(i4)                                  :: ii
  integer(i4)                                  :: il
  integer(i4)                                  :: iloc
  integer(i4)                                  :: ind
  integer(i4)                                  :: indj
  integer(i4)                                  :: inode
  integer(i4)                                  :: is
  integer(i4)                                  :: iv
  integer(i4)                                  :: j
  integer(i4)                                  :: jb1
  integer(i4)                                  :: jb2
  integer(i4)                                  :: jeem
  integer(i4)                                  :: jj
  integer(i4)                                  :: jl
  integer(i4)                                  :: jloc
  integer(i4)                                  :: jx
  integer(i4)                                  :: jy
  integer(i4)                                  :: jz
  integer(i4)                                  :: k
  integer(i4)                                  :: kk
  integer(i4)                                  :: kl
  integer(i4)                                  :: ml1
  integer(i4)                                  :: ml2
  integer(i4)                                  :: ml3
  integer(i4)                                  :: n
  integer(i4)                                  :: nbi
  integer(i4)                                  :: nbinode
  integer(i4)                                  :: nbj
  integer(i4)                                  :: nbq
  integer(i4)                                  :: nbqj
  integer(i4)                                  :: nbqo
  integer(i4)                                  :: nbqjo
  integer(i4)                                  :: ndqmax
  integer(i4)                                  :: neemfull
  integer(i4), dimension(:),   allocatable     :: neemfullptr
  integer(i4), dimension(:),   allocatable     :: neemfullrptr
  integer(i4)                                  :: ni
  integer(i4)                                  :: nj
  integer(i4)                                  :: npqni
  integer(i4)                                  :: npqnj
  integer(i4)                                  :: nqr
  integer(i4)                                  :: nqrj
  integer(i4)                                  :: nqrib1
  integer(i4)                                  :: nqrib2
  integer(i4)                                  :: ns
  integer(i4)                                  :: status
!
  integer                                      :: MPIerror
  integer                                      :: ntag
  integer                                      :: nnode
  integer                                      :: ntmp
  integer(i4), dimension(:),   allocatable     :: nbtmp
  integer                                      :: Request
  integer(i4), dimension(:),   allocatable     :: StatMPI       ! Array for status from MPI
!
  logical                                      :: lfound
  logical                                      :: lstrain
  real(dp)                                     :: accf2
  real(dp)                                     :: arg
  real(dp)                                     :: argtest
  real(dp)                                     :: cosa
  real(dp)                                     :: costrm
  real(dp)                                     :: cuts2
  real(dp)                                     :: d2zetaii
  real(dp)                                     :: d2zetaij
  real(dp)                                     :: d2zetajj
  real(dp)                                     :: d2zetari
  real(dp)                                     :: d2zetarj
  real(dp)                                     :: d2gamr2
  real(dp)                                     :: d2gamifj
  real(dp)                                     :: d2gamjfi
  real(dp)                                     :: darg1
  real(dp)                                     :: darg2
  real(dp)                                     :: dchis(6)
  real(dp),    dimension(:),   allocatable     :: dchix
  real(dp),    dimension(:),   allocatable     :: dchiy
  real(dp),    dimension(:),   allocatable     :: dchiz
  real(dp)                                     :: g_derf
  real(dp)                                     :: g_derfc
  real(dp)                                     :: derfc1
  real(dp)                                     :: derfc2
  real(dp)                                     :: derfez
  real(dp)                                     :: dexp1
  real(dp)                                     :: dexp2
  real(dp)                                     :: dexp3
  real(dp)                                     :: dexp4
  real(dp)                                     :: dexpz
  real(dp)                                     :: dqme(3)
  real(dp)                                     :: d2qme(6)
  real(dp)                                     :: dgam
  real(dp)                                     :: dgamifj
  real(dp)                                     :: dgamjfi
  real(dp),    dimension(:,:), allocatable     :: dqbondds
  real(dp),    dimension(:,:), allocatable     :: dqbonddxyz
  real(dp),    dimension(:,:), allocatable     :: dqs
  real(dp),    dimension(:),   allocatable     :: dqx
  real(dp),    dimension(:),   allocatable     :: dqy
  real(dp),    dimension(:),   allocatable     :: dqz
  real(dp),    dimension(:),   allocatable     :: dqtmp1
  real(dp),    dimension(:),   allocatable     :: dqtmp2
  real(dp),    dimension(:,:), allocatable     :: dqbloc
  real(dp),    dimension(:,:), allocatable     :: dqbsloc
  real(dp)                                     :: dr2ds(6)
  real(dp)                                     :: d2r2dx2(3,3)
  real(dp)                                     :: d2r2ds2(6,6)
  real(dp)                                     :: d2r2dsdx(6,3)
  real(dp)                                     :: dtrm1
  real(dp)                                     :: dtrm1i
  real(dp)                                     :: dtrm1zn
  real(dp)                                     :: dtrm1j
  real(dp)                                     :: dzetai
  real(dp)                                     :: dzetaj
  real(dp)                                     :: errfcn
  real(dp)                                     :: etaloc
  real(dp)                                     :: etaz
  real(dp)                                     :: etaz2
  real(dp)                                     :: etrm
  real(dp)                                     :: fieldx
  real(dp)                                     :: fieldy
  real(dp)                                     :: fieldz
  real(dp)                                     :: gam
  real(dp)                                     :: gamifj
  real(dp)                                     :: gamjfi
  real(dp)                                     :: Gmax
  real(dp)                                     :: kexperfc
  real(dp)                                     :: kvec
  real(dp)                                     :: qme
  real(dp)                                     :: qi
  real(dp)                                     :: qli
  real(dp)                                     :: qlii
  real(dp)                                     :: qlj
  real(dp)                                     :: qljj
  real(dp)                                     :: rconv
  real(dp)                                     :: rexp
  real(dp)                                     :: rkvec
  real(dp)                                     :: rl
  real(dp)                                     :: rmax
  real(dp)                                     :: rqeq2
  real(dp)                                     :: rr2
  real(dp)                                     :: rrr
  real(dp)                                     :: rv2
  real(dp)                                     :: rx
  real(dp)                                     :: rxi
  real(dp)                                     :: rxj
  real(dp)                                     :: rxk
  real(dp)                                     :: ry
  real(dp)                                     :: ry2
  real(dp)                                     :: ryi
  real(dp)                                     :: ryj
  real(dp)                                     :: ryk
  real(dp)                                     :: rz
  real(dp)                                     :: rz2
  real(dp)                                     :: rzi
  real(dp)                                     :: rzj
  real(dp)                                     :: rzk
  real(dp)                                     :: setaloc
  real(dp)                                     :: sina
  real(dp)                                     :: sineq
  real(dp)                                     :: smallestG
  real(dp)                                     :: strm
  real(dp)                                     :: strm1
  real(dp)                                     :: strm2
  real(dp)                                     :: sum
  real(dp)                                     :: sumx
  real(dp)                                     :: sumy
  real(dp)                                     :: sumz
  real(dp)                                     :: trmi
  real(dp)                                     :: xci
  real(dp)                                     :: yci
  real(dp)                                     :: zci
  real(dp)                                     :: xd
  real(dp)                                     :: yd
  real(dp)                                     :: zd
  real(dp)                                     :: zetai
  real(dp)                                     :: zetaj
  real(dp)                                     :: znucj
  real(dp)                                     :: ztrm1
  real(dp),    dimension(:),   allocatable     :: zb
#ifdef TRACE
  call trace_in('dchargesplitd')
#endif
!
!  Set local strain flag
!
  lstrain = (lstrainin.and.nstrains.gt.0)
!
!  Check the memory for the square arrays
!
  if (numat.gt.maxd2) then
    maxd2 = numat
    call changemaxd2
  endif
!
!  Check the memory for the charge derivatives
!
  if (natomsonnode.gt.maxd2qu) then
    maxd2qu = natomsonnode
    call changemaxd2q
  endif
  if (3*numat.gt.maxd2q) then
    maxd2q = 3*numat
    call changemaxd2q
  endif
!
!  Zero dq/dxyz and dq/ds
!
  do i = 1,natomsonnode
    do j = 1,3*numat
      dqdxyz(j,i) = 0.0_dp
    enddo
  enddo
  if (lstrain) then
    do i = 1,natomsonnode
      do j = 1,nstrains
        dqds(j,i) = 0.0_dp
      enddo
    enddo
  endif
!
!  If electrostatics have been turned off then there is no point in proceeding.
!
  if (.not.lDoElectrostatics) then
#ifdef TRACE
    call trace_out('dchargesplitd')
#endif
    return
  endif
!
!  Create local array for bond charge derivatives
!
  allocate(dqbonddxyz(3*numat,nbloc),stat=status)
  if (status/=0) call outofmemory('dchargesplitd','dqbonddxyz')
  if (lstrain) then
    allocate(dqbondds(nstrains,nbloc),stat=status)
    if (status/=0) call outofmemory('dchargesplitd','dqbondds')
  endif
!
!  Zero dqbond/dxyz and dqbond/ds
!
  do i = 1,nbloc
    do j = 1,3*numat
      dqbonddxyz(j,i) = 0.0_dp
    enddo
  enddo
  if (lstrain) then
    do i = 1,nbloc
      do j = 1,nstrains
        dqbondds(j,i) = 0.0_dp
      enddo
    enddo
  endif
!
  cuts2 = cuts*cuts
  rconv = 1.0_dp/autoangs
!
!  Find number of EEM active atoms
!
  allocate(neemfullptr(numat),stat=status)
  if (status/=0) call outofmemory('dchargesplitd','neemfullptr')
  allocate(neemfullrptr(numat),stat=status)
  if (status/=0) call outofmemory('dchargesplitd','neemfullrptr')
!
  neemfull = 0
  do i = 1,numat
    ni = nat(i)
    qi = qf(i)
    if (lelementOK(nat(i)).and.nregionno(nsft+nrelf2a(i)).eq.1) then
      neemfull = neemfull + 1
      neemfullptr(neemfull) = i
      neemfullrptr(i) = neemfull
      if (lmultiqrange) then
        if (nqrange(ni,neemtype).gt.1) then
          lfound = .false.
          nqr = 0 
          do while (.not.lfound.and.nqr.lt.nqrange(ni,neemtype))
            nqr = nqr + 1 
            if (nqrangetype(nqr,ni,neemtype).eq.3) then
              lfound = (qi.ge.qrangemin(nqr,ni,neemtype).and.qi.lt.qrangemax(nqr,ni,neemtype))
            elseif (nqrangetype(nqr,ni,neemtype).eq.2) then
              lfound = (qi.le.qrangemax(nqr,ni,neemtype))
            elseif (nqrangetype(nqr,ni,neemtype).eq.1) then
              lfound = (qi.ge.qrangemin(nqr,ni,neemtype))
            endif
          enddo
          if (.not.lfound) nqr = 1 
          nqrnow(neemfull) = nqr
        else
          nqrnow(neemfull) = 1
        endif
      endif
    endif
  enddo
!
  if (lfieldcfg(ncf).or.(ntdfieldcfg(ncf).gt.0)) then
!
!  Get electric field if needed
!
    call electricfieldparts(fieldx,fieldy,fieldz)
  endif
  if (lrecalcA) then
!****************
!  Symopt case  *
!****************
    allocate(zb(numat),stat=status)
    if (status/=0) call outofmemory('dchargesplitd','zb')
!
!  First generate 2 centre terms
!
    zb(1:numat) = 0.0_dp
    if (lnoqeem) then
      potl(1:neembond,1:nbloc) = 0.0_dp
    else
      call genpotsplitd(nbloc,nblocptr,potl,nmax,zb)
!
!  Now create matrix for split bonds by reducing left-hand dimension
!  Right-hand dimension is already based on local split bonds
!
      do il = 1,nbloc
        zb(1:numat) = potl(1:numat,il)
        potl(1:neembond,il) = 0.0_dp
        do j = 1,neembond
          jb1 = neembonded(1,j)
          jb2 = neembonded(2,j)
          potl(j,il) = potl(j,il) + zb(jb1) - zb(jb2)
        enddo
      enddo
    endif
!
!  Add one centre terms
!
    if (nbloc.gt.0) then
      do il = 1,nbloc
        i = nblocptr(il)
        ib1 = neembonded(1,i)
        ib2 = neembonded(2,i)
        if (lmultiqrange) then
          nqrib1 = nqrnow(neemrptr(ib1))
          nqrib2 = nqrnow(neemrptr(ib2))
        else
          nqrib1 = 1
          nqrib2 = 1
        endif
        do j = 1,neembond
          jb1 = neembonded(1,j)
          jb2 = neembonded(2,j)
          if (ib1.eq.jb1) then
            potl(j,il) = potl(j,il) + 2.0_dp*murange(nqrib1,iatn(ib1),neemtype)
          endif
          if (ib1.eq.jb2) then
            potl(j,il) = potl(j,il) - 2.0_dp*murange(nqrib1,iatn(ib1),neemtype)
          endif
          if (ib2.eq.jb1) then
            potl(j,il) = potl(j,il) - 2.0_dp*murange(nqrib2,iatn(ib2),neemtype)
          endif
          if (ib2.eq.jb2) then
            potl(j,il) = potl(j,il) + 2.0_dp*murange(nqrib2,iatn(ib2),neemtype)
          endif
        enddo
      enddo
    endif
!     
!  Matrix inversion 
!     
    ifail = 0
    n = neembond
    call matrix_inversion_library(n,1_i4,nmax,nblocksize,potl,0_i4,ifail)
!
!  Was inversion successful?
!
    if (ifail.ne.0) then
      call outerror('matrix inversion failed in dchargesplitd',0_i4)
      call stopnow('dchargesplitd')
    endif
!
    deallocate(zb,stat=status)
    if (status/=0) call deallocate_error('dchargesplitd','zb')
  endif
  rqeq2 = rqeq*rqeq
!
!  Allocate local memory
!
  ndqmax = max(neembond,numat)
  allocate(dqx(ndqmax),stat=status)
  if (status/=0) call outofmemory('dchargesplitd','dqx')
  allocate(dqy(ndqmax),stat=status)
  if (status/=0) call outofmemory('dchargesplitd','dqy')
  allocate(dqz(ndqmax),stat=status)
  if (status/=0) call outofmemory('dchargesplitd','dqz')
!
  allocate(dqtmp1(max(3_i4,nstrains)*ndqmax+nstrains),stat=status)
  if (status/=0) call outofmemory('dchargesplitd','dqtmp1')
  allocate(dqtmp2(max(3_i4,nstrains)*ndqmax+nstrains),stat=status)
  if (status/=0) call outofmemory('dchargesplitd','dqtmp2')
!
  if (lSandM) then
    allocate(dchix(numat),stat=status)
    if (status/=0) call outofmemory('dchargesplitd','chix')
    allocate(dchiy(numat),stat=status)
    if (status/=0) call outofmemory('dchargesplitd','chiy')
    allocate(dchiz(numat),stat=status)
    if (status/=0) call outofmemory('dchargesplitd','chiz')
  endif
  if (lstrain) then
    allocate(dqs(ndqmax,nstrains),stat=status)
    if (status/=0) call outofmemory('dchargesplitd','dqs')
!
!  Zero temporary strain derivative arrays
!
    do i = 1,nstrains
      do j = 1,ndqmax
        dqs(j,i) = 0.0_dp
      enddo
    enddo
  endif
  if (ndim.gt.1.or.lwolf) then
!******************
!  Periodic case  *
!******************
!
!  Define constants
!
    if (lwolf) then
      radmax = cutw
      etaloc = etaw*etaw
      setaloc = etaw
    elseif (lewald) then
      if (ndim.eq.2) then
        rpieta = 1.0_dp / sqrt(pi * eta)
        rhseta = 0.5_dp / seta
        accf2 = accf*accf
        argtest = sqrt(3.0+0.5*accf2) - sqrt(3.0)
        smallestG = min(kv(1,1),kv(2,2))
      endif
      radmax = accf/seta
      eta4 = 0.25_dp/eta
      etaloc = eta
      setaloc = seta
    else
      radmax = 0.0_dp
      eta4 = 0.25_dp
      etaloc = eta
      setaloc = seta
    endif
!
    if (lqeq.or.lSandM) then
      rmax = max(radmax,rqeq)
    else
      rmax = radmax
    endif
    rmax2 = rmax*rmax
!
!  Estimate upper limits for looping
!
    if (ndim.eq.3) then
      rv2 = rv(1,1)**2 + rv(2,1)**2 + rv(3,1)**2
      rv2 = sqrt(rv2)
      ml1 = rmax/rv2 + 1
      rv2 = rv(1,2)**2 + rv(2,2)**2 + rv(3,2)**2
      rv2 = sqrt(rv2)
      ml2 = rmax/rv2 + 1
      rv2 = rv(1,3)**2 + rv(2,3)**2 + rv(3,3)**2
      rv2 = sqrt(rv2)
      ml3 = rmax/rv2 + 1
    elseif (ndim.eq.2) then
      rv2 = rv(1,1)**2 + rv(2,1)**2
      rv2 = sqrt(rv2)
      ml1 = rmax/rv2 + 1
      rv2 = rv(1,2)**2 + rv(2,2)**2
      rv2 = sqrt(rv2)
      ml2 = rmax/rv2 + 1
      ml3 = 0
    elseif (ndim.eq.1) then
      rv2 = rv(1,1)
      ml1 = rmax/rv2 + 1
      ml2 = 0
      ml3 = 0
    elseif (ndim.eq.0) then
      ml1 = 0
      ml2 = 0
      ml3 = 0
    endif
!**********************************
!  Reciprocal space contribution  *
!**********************************
    if (lnorecip.or.lwolf) goto 5
!
!  Set up strain terms for reciprocal space
!
    if (lstrain) then
      call gstrterms(ndim,maxkvec,nkvec,xrk,yrk,zrk,dg2ds,d2g2dx2,d2g2ds2,.false.)
    endif
!
!  Start loop over atoms for coordinate differentiation
!
    ieem = 0
    do i = 1,numat
      xci = xclat(i)
      yci = yclat(i)
      zci = zclat(i)
      qli = qf(i)
      ind = 3*(i - 1)
      if (lelementOK(nat(i)).and.nregionno(nsft+nrelf2a(i)).eq.1) ieem = ieem + 1
!
!  Zero temporary derivative arrays
!
      dqx(1:neembond) = 0.0_dp
      dqy(1:neembond) = 0.0_dp
      dqz(1:neembond) = 0.0_dp
!
      do jloc = 1,natomsonnode
        j = node2atom(jloc)
        jeem = neemfullrptr(j)
        xd = xclat(j) - xci
        yd = yclat(j) - yci
        zd = zclat(j) - zci
        qlj = qf(j)
        if (lelementOK(nat(j)).and.nregionno(nsft+nrelf2a(j)).eq.1) jeem = jeem + 1
!
!  Evaluate derivative of reciprocal space elements of A
!
        if (ndim.eq.3) then
          do iv = 1,nkvec
            argc(iv) = xrk(iv)*xd + yrk(iv)*yd + zrk(iv)*zd
            sine(iv) = sin(argc(iv))*ktrm(iv)
            do nbi = 1,nbonds(i)
              nbqo = nbondqb(nbi,i)
              nbq = abs(nbqo)
              if (nbqo.gt.0) then
                dqx(nbq) = dqx(nbq) + sine(iv)*xrk(iv)*qlj
                dqy(nbq) = dqy(nbq) + sine(iv)*yrk(iv)*qlj
                dqz(nbq) = dqz(nbq) + sine(iv)*zrk(iv)*qlj
                if (lstrain) then
                  costrm = cos(argc(iv))*angstoev*qlj
                  strm1 = costrm*ktrms(iv)
                  strm2 = costrm*ktrm(iv)
                  if (lfinitestrain) then
                    do is = 1,nstrains
                      dqs(nbq,is) = dqs(nbq,is) + strm1*dg2ds(iv,is) - strm2*strainddetds(is)*straindet
                    enddo
                  else
                    do is = 1,nstrains
                      dqs(nbq,is) = dqs(nbq,is) + strm1*dg2ds(iv,is)
                    enddo
                    dqs(nbq,1) = dqs(nbq,1) - strm2
                    dqs(nbq,2) = dqs(nbq,2) - strm2
                    dqs(nbq,3) = dqs(nbq,3) - strm2
                  endif
                endif
              elseif (nbqo.lt.0) then
                dqx(nbq) = dqx(nbq) - sine(iv)*xrk(iv)*qlj
                dqy(nbq) = dqy(nbq) - sine(iv)*yrk(iv)*qlj
                dqz(nbq) = dqz(nbq) - sine(iv)*zrk(iv)*qlj
                if (lstrain) then
                  costrm = cos(argc(iv))*angstoev*qlj
                  strm1 = costrm*ktrms(iv)
                  strm2 = costrm*ktrm(iv)
                  if (lfinitestrain) then
                    do is = 1,nstrains
                      dqs(nbq,is) = dqs(nbq,is) - strm1*dg2ds(iv,is) + strm2*strainddetds(is)*straindet
                    enddo
                  else
                    do is = 1,nstrains
                      dqs(nbq,is) = dqs(nbq,is) - strm1*dg2ds(iv,is)
                    enddo
                    dqs(nbq,1) = dqs(nbq,1) + strm2
                    dqs(nbq,2) = dqs(nbq,2) + strm2
                    dqs(nbq,3) = dqs(nbq,3) + strm2
                  endif
                endif
              endif
            enddo
            do nbj = 1,nbonds(j)
              nbqo = nbondqb(nbj,j)
              nbq = abs(nbqo)
              if (nbqo.gt.0) then
                dqx(nbq) = dqx(nbq) + sine(iv)*xrk(iv)*qli
                dqy(nbq) = dqy(nbq) + sine(iv)*yrk(iv)*qli
                dqz(nbq) = dqz(nbq) + sine(iv)*zrk(iv)*qli
              elseif (nbqo.lt.0) then
                dqx(nbq) = dqx(nbq) - sine(iv)*xrk(iv)*qli
                dqy(nbq) = dqy(nbq) - sine(iv)*yrk(iv)*qli
                dqz(nbq) = dqz(nbq) - sine(iv)*zrk(iv)*qli
              endif
            enddo
          enddo
        elseif (ndim.eq.2) then
!
!  First term - K vector independent
!
          etaz = seta*zd
          etaz2 = etaz*etaz
          derfez = g_derf(etaz)
          dexpz  = exp(-etaz2)
          etrm   = - vol4pi*(zd*derfez + dexpz*rpieta)*angstoev
          dtrm1  = - vol4pi*derfez
          do nbi = 1,nbonds(i)
            nbqo = nbondqb(nbi,i)
            nbq = abs(nbqo)
            if (nbqo.gt.0) then
              dqz(nbq) = dqz(nbq) - dtrm1*qlj
              if (lstrain) then
                dqs(nbq,1) = dqs(nbq,1) - etrm*qlj
                dqs(nbq,2) = dqs(nbq,2) - etrm*qlj
              endif
            elseif (nbqo.lt.0) then
              dqz(nbq) = dqz(nbq) + dtrm1*qlj
              if (lstrain) then
                dqs(nbq,1) = dqs(nbq,1) + etrm*qlj
                dqs(nbq,2) = dqs(nbq,2) + etrm*qlj
              endif
            endif
          enddo
          do nbj = 1,nbonds(j)
            nbqo = nbondqb(nbj,j)
            nbq = abs(nbqo)
            if (nbqo.gt.0) then
              dqz(nbq) = dqz(nbq) - dtrm1*qli
            elseif (nbqo.lt.0) then
              dqz(nbq) = dqz(nbq) + dtrm1*qli
            endif
          enddo
!
!  Find local kvector cut-off
!
          if (abs(etaz).gt.argtest) then
            Gmax = abs(accf2/zd)
          else
            Gmax = sqrt(4.0*eta*(accf2-etaz2))
          endif
          if (Gmax.ge.smallestG) then
            do iv = 1,nkvec
              kvec = kmod(iv)
              if (kvec.le.Gmax) then
                arg = xrk(iv)*xd + yrk(iv)*yd
                sina = sin(arg)*ktrm(iv)
                cosa = cos(arg)*ktrm(iv)
                dexp1 = exp(kvec*zd)
                dexp2 = 1.0_dp/dexp1
                darg1 = kvec*rhseta + etaz
                darg2 = kvec*rhseta - etaz
                dexp3 = exp(-(darg1)**2)
                dexp4 = exp(-(darg2)**2)
                derfc1 = g_derfc(darg1)
                derfc2 = g_derfc(darg2)
                kexperfc = dexp1*derfc1 + dexp2*derfc2
                sineq = sina*kexperfc
                ztrm1 = (kvec*(dexp1*derfc1-dexp2*derfc2) - tweatpi*(dexp1*dexp3-dexp2*dexp4))
                do nbi = 1,nbonds(i)
                  nbqo = nbondqb(nbi,i)
                  nbq = abs(nbqo)
                  if (nbqo.gt.0) then
                    dqx(nbq) = dqx(nbq) + sineq*xrk(iv)*qlj
                    dqy(nbq) = dqy(nbq) + sineq*yrk(iv)*qlj
                    dqz(nbq) = dqz(nbq) - cosa*ztrm1*qlj
                    if (lstrain) then
                      costrm = cosa*angstoev*qlj
                      rkvec = 1.0_dp/kvec
                      strm = rkvec*(-rkvec*kexperfc + zd*(dexp1*derfc1-dexp2*derfc2) - &
                             rpieta*(dexp1*dexp3+dexp2*dexp4))
                      strm1 = costrm*strm
                      strm2 = costrm*kexperfc
                      if (lfinitestrain) then
                        do is = 1,nstrains
                          dqs(nbq,is) = dqs(nbq,is) + strm1*dg2ds(iv,is) - strm2*strainddetds(is)*straindet
                        enddo
                      else
                        do is = 1,nstrains
                          dqs(nbq,is) = dqs(nbq,is) + strm1*dg2ds(iv,is)
                        enddo
                        dqs(nbq,1) = dqs(nbq,1) - strm2
                        dqs(nbq,2) = dqs(nbq,2) - strm2
                      endif
                    endif
                  elseif (nbqo.lt.0) then
                    dqx(nbq) = dqx(nbq) - sineq*xrk(iv)*qlj
                    dqy(nbq) = dqy(nbq) - sineq*yrk(iv)*qlj
                    dqz(nbq) = dqz(nbq) + cosa*ztrm1*qlj
                    if (lstrain) then
                      costrm = cosa*angstoev*qlj
                      rkvec = 1.0_dp/kvec
                      strm = rkvec*(-rkvec*kexperfc + zd*(dexp1*derfc1-dexp2*derfc2) - &
                             rpieta*(dexp1*dexp3+dexp2*dexp4))
                      strm1 = costrm*strm
                      strm2 = costrm*kexperfc
                      if (lfinitestrain) then
                        do is = 1,nstrains
                          dqs(nbq,is) = dqs(nbq,is) - strm1*dg2ds(iv,is) + strm2*strainddetds(is)*straindet
                        enddo
                      else
                        do is = 1,nstrains
                          dqs(nbq,is) = dqs(nbq,is) - strm1*dg2ds(iv,is)
                        enddo
                        dqs(nbq,1) = dqs(nbq,1) + strm2
                        dqs(nbq,2) = dqs(nbq,2) + strm2
                      endif
                    endif
                  endif
                enddo
                do nbj = 1,nbonds(j)
                  nbqo = nbondqb(nbj,j)
                  nbq = abs(nbqo)
                  if (nbqo.gt.0) then
                    dqx(nbq) = dqx(nbq) + sineq*xrk(iv)*qli
                    dqy(nbq) = dqy(nbq) + sineq*yrk(iv)*qli
                    dqz(nbq) = dqz(nbq) - cosa*ztrm1*qli
                  elseif (nbqo.lt.0) then
                    dqx(nbq) = dqx(nbq) - sineq*xrk(iv)*qli
                    dqy(nbq) = dqy(nbq) - sineq*yrk(iv)*qli
                    dqz(nbq) = dqz(nbq) + cosa*ztrm1*qli
                  endif
                enddo
              endif
            enddo
          endif
        endif
!
!  End inner atom loop
!
      enddo
!
!  Global sum of dqx etc
!
      dqtmp1(1:neembond) = dqx(1:neembond)
      dqtmp1(neembond+1:2*neembond) = dqy(1:neembond)
      dqtmp1(2*neembond+1:3*neembond) = dqz(1:neembond)
!
      call sumall(dqtmp1,dqtmp2,3_i4*neembond,"dchargesplitd","dqtmp")
!
      dqx(1:neembond) = dqtmp2(1:neembond)
      dqy(1:neembond) = dqtmp2(neembond+1:2*neembond)
      dqz(1:neembond) = dqtmp2(2*neembond+1:3*neembond)
!
!  Multiply dqx/dqy/dqz by inverse matrix
!
      do jl = 1,nbloc
        sumx = 0.0_dp
        sumy = 0.0_dp
        sumz = 0.0_dp
        do k = 1,neembond
          sumx = sumx + dqx(k)*potl(k,jl)
          sumy = sumy + dqy(k)*potl(k,jl)
          sumz = sumz + dqz(k)*potl(k,jl)
        enddo
        sumx = sumx*angstoev
        sumy = sumy*angstoev
        sumz = sumz*angstoev
        dqbonddxyz(ind+1,jl) = dqbonddxyz(ind+1,jl) - sumx
        dqbonddxyz(ind+2,jl) = dqbonddxyz(ind+2,jl) - sumy
        dqbonddxyz(ind+3,jl) = dqbonddxyz(ind+3,jl) - sumz
      enddo
    enddo
5   continue
!*************************
!  Real space summation  *
!*************************
    if (lnoreal.and.lstrain) then
!
!  Global sum of dqs
!
      do kl = 1,nstrains
        do k = 1,neembond
          dqtmp1((kl-1)*nstrains+k) = dqs(k,kl)
        enddo
      enddo
!
      call sumall(dqtmp1,dqtmp2,nstrains*neembond,"dchargesplitd","dqstmp")
!
      do kl = 1,nstrains
        do k = 1,neembond
          dqs(k,kl) = dqtmp2((kl-1)*nstrains+k)
        enddo
      enddo
!
!  Strain terms
!
      do jl = 1,nbloc
        do kl = 1,nstrains
          sum = 0.0_dp
          do k = 1,neembond
            sum = sum + dqs(k,kl)*potl(k,jl)
          enddo
          dqbondds(kl,jl) = dqbondds(kl,jl) - sum
        enddo
      enddo
      goto 135
    endif
!
!  Start loop over atoms for coordinate differentiation
!
    ieem = 0
    do i = 1,numat
      qlii = qf(i)
      xci = xclat(i)
      yci = yclat(i)
      zci = zclat(i)
      ind = 3*(i-1)
      ni = nat(i)
      if (lelementOK(nat(i)).and.nregionno(nsft+nrelf2a(i)).eq.1) then
        ieem = ieem + 1
!
!  If QEq work out principal quantum number
!
        if (lqeq) then
          if (lmultiqrange) then
            nqr = nqrnow(ieem)
          else
            nqr = 1
          endif
          if (ni.le.2) then
            npqni = 1
          elseif (ni.le.10) then
            npqni = 2
          elseif (ni.le.18) then
            npqni = 3
          elseif (ni.le.36) then
            npqni = 4
          elseif (ni.le.54) then
            npqni = 5
          elseif (ni.le.86) then
            npqni = 6
          else
            npqni = 7
          endif
          zetai = 0.5_dp*qeqlambda*(2*npqni+1)/radrange(nqr,ni,neemtype)
        elseif (lSandM) then   
          if (lmultiqrange) then
            nqr = nqrnow(ieem)
          else
            nqr = 1
          endif
          zetai = zetarange(nqr,ni,neemtype)
        endif
      endif
!
!  Zero temporary derivative arrays
!
      dqx(1:neembond) = 0.0_dp
      dqy(1:neembond) = 0.0_dp
      dqz(1:neembond) = 0.0_dp
      if (lSandM) then
        dchix(1:numat) = 0.0_dp
        dchiy(1:numat) = 0.0_dp
        dchiz(1:numat) = 0.0_dp
        dchis(1:6) = 0.0_dp
      endif
!
      do jloc = 1,natomsonnode
        j = node2atom(jloc)
        jeem = neemfullrptr(j)
        qljj = qf(j)
        qlj = qljj*occuf(j)
        rx = xclat(j) - xci
        ry = yclat(j) - yci
        rz = zclat(j) - zci
        nj = nat(j)
!
!  If QEq work out principal quantum number
!
        if (lqeq) then
          if (lmultiqrange) then
            nqrj = nqrnow(jeem)
          else
            nqrj = 1
          endif
          if (nj.le.2) then
            npqnj = 1
          elseif (nj.le.10) then
            npqnj = 2
          elseif (nj.le.18) then
            npqnj = 3
          elseif (nj.le.36) then
            npqnj = 4
          elseif (nj.le.54) then
            npqnj = 5
          elseif (nj.le.86) then
            npqnj = 6
          else
            npqnj = 7
          endif
          zetaj = 0.5_dp*qeqlambda*(2*npqnj+1)/radrange(nqrj,nj,neemtype)
        elseif (lSandM) then
          if (lmultiqrange) then
            nqrj = nqrnow(jeem)
          else
            nqrj = 1
          endif
          zetaj = zetarange(nqrj,nj,neemtype)
          znucj = znucrange(nqrj,nj,neemtype)
        endif
!
!  Loop over cell vectors
!
        rxi = rx - (ml1 + 1)*r1x
        ryi = ry - (ml1 + 1)*r1y
        rzi = rz - (ml1 + 1)*r1z
        do ii = -ml1,ml1
          rxi = rxi + r1x
          ryi = ryi + r1y
          rzi = rzi + r1z
          rxj = rxi - (ml2+1)*r2x
          ryj = ryi - (ml2+1)*r2y
          rzj = rzi - (ml2+1)*r2z
          do jj = -ml2,ml2
            rxj = rxj + r2x
            ryj = ryj + r2y
            rzj = rzj + r2z
            rxk = rxj - (ml3+1)*r3x
            ryk = ryj - (ml3+1)*r3y
            rzk = rzj - (ml3+1)*r3z
            do 120 kk = -ml3,ml3
              rxk = rxk + r3x
              ryk = ryk + r3y
              rzk = rzk + r3z
!
!  Calculate distance squared
!
              rr2 = rxk*rxk + ryk*ryk + rzk*rzk
!
!  Exclude distances outside maximum cutoff
!
              if (rr2.gt.rmax2) goto 120
!
!  Trap self term
!
              if (rr2.lt.1.0d-15) goto 120
              rl = sqrt(rr2)
              rrr = 1.0_dp/rl
              if (rr2.lt.cuts2) then
!
!  Core-shell interaction
!
                dtrm1 = rrr*rrr*rrr
                dtrm1i = dtrm1
                dtrm1j = dtrm1
                dtrm1zn = 0.0_dp
              elseif (rr2.lt.rqeq2.and.lqeq) then
!
!  Calculate Coulomb interaction according to QEq scheme and
!  subtract 1/r term from Ewald sum.
!
                call gammas(npqni,npqnj,zetai,zetaj,rl,gam,dgam,dzetai,dzetaj, &
                            d2zetaii,d2zetaij,d2zetajj,d2zetari,d2zetarj,d2gamr2)
                dtrm1 = (rrr*rrr+dgam)*rrr
                dtrm1i = dtrm1
                dtrm1j = dtrm1
              elseif (rr2.lt.rqeq2.and.lSandM) then
                call gammasm(zetai,zetaj,rl,gam,dgam,d2gamr2,gamifj,gamjfi,dgamifj,dgamjfi,d2gamifj,d2gamjfi)
                dtrm1 = dgam*rrr
                dtrm1i = dtrm1
                dtrm1j = dtrm1
                dtrm1zn = znucj*(dgamjfi - dgam)*rrr
              else
                dtrm1i = 0.0_dp
                dtrm1j = 0.0_dp
                dtrm1zn = 0.0_dp
              endif
!
!  Complementary error function
!
              errfcn = g_derfc(setaloc*rl)
              trmi = errfcn/rl
              rexp = tweatpi*exp(-etaloc*rr2)
              dtrm1 = (trmi+rexp)*rrr*rrr
              dtrm1i = dtrm1i - dtrm1
              dtrm1j = dtrm1j - dtrm1
!
!  First derivatives of matrix elements
!
              dtrm1i = dtrm1i*angstoev
              dtrm1j = dtrm1j*angstoev
              if (lstrain) then
                call real1strterm(ndim,rxk,ryk,rzk,0.0_dp,0.0_dp,0.0_dp,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
              endif
              if (lSandM) then
                dtrm1zn = dtrm1zn*angstoev
                if (lelementOK(nat(i)).and.nregionno(nsft+nrelf2a(i)).eq.1) then
                  dchix(ieem) = dchix(ieem) - dtrm1zn*rxk
                  dchiy(ieem) = dchiy(ieem) - dtrm1zn*ryk
                  dchiz(ieem) = dchiz(ieem) - dtrm1zn*rzk
                  if (lstrain) then
                    do is = 1,nstrains
                      ns = nstrptr(is)
                      dchis(is) = dchis(is) + dtrm1zn*dr2ds(ns)
                    enddo
                  endif
                endif
                if (lelementOK(nat(j)).and.nregionno(nsft+nrelf2a(j)).eq.1) then
                  dchix(jeem) = dchix(jeem) + dtrm1zn*rxk
                  dchiy(jeem) = dchiy(jeem) + dtrm1zn*ryk
                  dchiz(jeem) = dchiz(jeem) + dtrm1zn*rzk
                endif
              endif
!
              do nbi = 1,nbonds(i)
                nbqo = nbondqb(nbi,i)
                nbq = abs(nbqo)
                if (nbqo.gt.0) then
                  dqx(nbq) = dqx(nbq) - qljj*dtrm1i*rxk
                  dqy(nbq) = dqy(nbq) - qljj*dtrm1i*ryk
                  dqz(nbq) = dqz(nbq) - qljj*dtrm1i*rzk
                  if (lstrain) then
                    do is = 1,nstrains
                      ns = nstrptr(is)
                      dqs(nbq,is) = dqs(nbq,is) + qljj*dtrm1i*dr2ds(ns)
                    enddo
                  endif
                elseif (nbqo.lt.0) then
                  dqx(nbq) = dqx(nbq) + qljj*dtrm1i*rxk
                  dqy(nbq) = dqy(nbq) + qljj*dtrm1i*ryk
                  dqz(nbq) = dqz(nbq) + qljj*dtrm1i*rzk
                  if (lstrain) then
                    do is = 1,nstrains
                      ns = nstrptr(is)
                      dqs(nbq,is) = dqs(nbq,is) - qljj*dtrm1i*dr2ds(ns)
                    enddo
                  endif
                endif
              enddo
              do nbj = 1,nbonds(j)
                nbqo = nbondqb(nbj,j)
                nbq = abs(nbqo)
                if (nbqo.gt.0) then
                  dqx(nbq) = dqx(nbq) - qlii*dtrm1j*rxk
                  dqy(nbq) = dqy(nbq) - qlii*dtrm1j*ryk
                  dqz(nbq) = dqz(nbq) - qlii*dtrm1j*rzk
                elseif (nbqo.lt.0) then
                  dqx(nbq) = dqx(nbq) + qlii*dtrm1j*rxk
                  dqy(nbq) = dqy(nbq) + qlii*dtrm1j*ryk
                  dqz(nbq) = dqz(nbq) + qlii*dtrm1j*rzk
                endif
              enddo
!
!  End of loops over lattice vectors
!
120         continue
          enddo
        enddo
!
!  End of loop over inner atom 
!
      enddo
!
!  Global sum of dqx, dqy, dqz
!
      dqtmp1(1:neembond) = dqx(1:neembond)
      dqtmp1(neembond+1:2*neembond) = dqy(1:neembond)
      dqtmp1(2*neembond+1:3*neembond) = dqz(1:neembond)
!
      call sumall(dqtmp1,dqtmp2,3_i4*neembond,"dchargesplitd","dqtmp")
!
      dqx(1:neembond) = dqtmp2(1:neembond)
      dqy(1:neembond) = dqtmp2(neembond+1:2*neembond)
      dqz(1:neembond) = dqtmp2(2*neembond+1:3*neembond)
!
!  Multiply dqx/dqy/dqz by inverse matrix
!
      do jl = 1,nbloc
        sumx = 0.0_dp
        sumy = 0.0_dp
        sumz = 0.0_dp
        do k = 1,neembond
          sumx = sumx + dqx(k)*potl(k,jl)
          sumy = sumy + dqy(k)*potl(k,jl)
          sumz = sumz + dqz(k)*potl(k,jl)
        enddo
        dqbonddxyz(ind+1,jl) = dqbonddxyz(ind+1,jl) - sumx
        dqbonddxyz(ind+2,jl) = dqbonddxyz(ind+2,jl) - sumy
        dqbonddxyz(ind+3,jl) = dqbonddxyz(ind+3,jl) - sumz
      enddo
      if (lSandM) then
!
!  Global sum of dchix, dchiy, dchiz, dchis
!
        dqtmp1(1:numat) = dchix(1:numat)
        dqtmp1(numat+1:2*numat) = dchiy(1:numat)
        dqtmp1(2*numat+1:3*numat) = dchiz(1:numat)
!
        if (lstrain) then
          dqtmp1(3*numat+1:3*numat+nstrains) = dchis(1:nstrains)
          call sumall(dqtmp1,dqtmp2,3_i4*numat+nstrains,"dchargesplitd","dchitmp")
        else
          call sumall(dqtmp1,dqtmp2,3_i4*numat,"dchargesplitd","dchitmp")
        endif
!
        dchix(1:numat) = dqtmp2(1:numat)
        dchiy(1:numat) = dqtmp2(numat+1:2*numat)
        dchiz(1:numat) = dqtmp2(2*numat+1:3*numat)
!
        if (lstrain) then
          dchis(1:nstrains) = dqtmp2(3*numat+1:3*numat+nstrains)
        endif
!
!  For S & M, multiply dchix/dchiy/dchiz by inverse matrix
!
        do nbi = 1,nbonds(i)
          nbqo = nbondqb(nbi,i)
          nbq = abs(nbqo)
          if (nbqo.gt.0) then
            do jl = 1,nbloc
              do k = 1,neemfull
                indj = 3*(neemfullptr(k) - 1)
                dqbonddxyz(indj+1,jl) = dqbonddxyz(indj+1,jl) - dchix(k)*potl(nbq,jl)
                dqbonddxyz(indj+2,jl) = dqbonddxyz(indj+2,jl) - dchiy(k)*potl(nbq,jl)
                dqbonddxyz(indj+3,jl) = dqbonddxyz(indj+3,jl) - dchiz(k)*potl(nbq,jl)
              enddo
            enddo
          elseif (nbqo.lt.0) then
            do jl = 1,nbloc
              do k = 1,neemfull
                indj = 3*(neemfullptr(k) - 1)
                dqbonddxyz(indj+1,jl) = dqbonddxyz(indj+1,jl) + dchix(k)*potl(nbq,jl)
                dqbonddxyz(indj+2,jl) = dqbonddxyz(indj+2,jl) + dchiy(k)*potl(nbq,jl)
                dqbonddxyz(indj+3,jl) = dqbonddxyz(indj+3,jl) + dchiz(k)*potl(nbq,jl)
              enddo
            enddo
          endif
        enddo
      endif
!
!  Strain terms
!
      if (lSandM.and.lstrain) then
!
!  For S & M, multiply dchis by inverse matrix
!                 
        do nbi = 1,nbonds(i)
          nbqo = nbondqb(nbi,i)
          nbq = abs(nbqo)
          if (nbqo.gt.0) then
            do jl = 1,nbloc
              do kl = 1,nstrains
                dqbondds(kl,jl) = dqbondds(kl,jl) - dchis(kl)*potl(nbq,jl)
              enddo
            enddo
          elseif (nbqo.lt.0) then
            do jl = 1,nbloc
              do kl = 1,nstrains
                dqbondds(kl,jl) = dqbondds(kl,jl) + dchis(kl)*potl(nbq,jl)
              enddo
            enddo
          endif
        enddo
      endif
!
!  End loop over i
!
    enddo
    if (lstrain) then
!
!  Global sum of dqs
!
      do kl = 1,nstrains
        do k = 1,neembond
          dqtmp1((kl-1)*neembond+k) = dqs(k,kl)
        enddo
      enddo
!
      call sumall(dqtmp1,dqtmp2,nstrains*neembond,"dchargesplitd","dqstmp")
!
      do kl = 1,nstrains
        do k = 1,neembond
          dqs(k,kl) = dqtmp2((kl-1)*neembond+k)
        enddo
      enddo
!
!  Strain terms
!
      do jl = 1,nbloc
        do kl = 1,nstrains
          sum = 0.0_dp
          do k = 1,neembond
            sum = sum + dqs(k,kl)*potl(k,jl)
          enddo
          dqbondds(kl,jl) = dqbondds(kl,jl) - sum
        enddo
      enddo
    endif
!**********************
!  End periodic case  *
!**********************
  else
!***********************
!  Cluster case / 1-D  *
!***********************
    if (lnoreal) goto 135
    if (ndim.eq.1) then
      call setmaxcell1D(maxloop(1))
      ml1 = maxloop(1)
    else
      ml1 = 0
      r1x = 0.0_dp
    endif
!
!  Start loop over cluster atoms
!
    ieem = 0
    do i = 1,numat
      qlii = qf(i)
      xci = xclat(i)
      yci = yclat(i)
      zci = zclat(i)
      ind = 3*(i-1)
      ni = nat(i)
      if (lelementOK(nat(i)).and.nregionno(nsft+nrelf2a(i)).eq.1) then
        ieem = ieem + 1
!
!  If QEq work out principal quantum number
!
        if (lqeq) then
          if (lmultiqrange) then
            nqr = nqrnow(ieem)
          else
            nqr = 1
          endif
          if (ni.le.2) then
            npqni = 1
          elseif (ni.le.10) then
            npqni = 2
          elseif (ni.le.18) then
            npqni = 3
          elseif (ni.le.36) then
            npqni = 4
          elseif (ni.le.54) then
            npqni = 5
          elseif (ni.le.86) then
            npqni = 6
          else
            npqni = 7
          endif
          zetai = 0.5_dp*qeqlambda*(2*npqni+1)/radrange(nqr,ni,neemtype)
        elseif (lSandM) then
          if (lmultiqrange) then
            nqr = nqrnow(ieem)
          else
            nqr = 1
          endif
          zetai = zetarange(nqr,ni,neemtype)
        endif
      endif
!
!  Zero temporary derivative arrays
!
      dqx(1:neembond) = 0.0_dp
      dqy(1:neembond) = 0.0_dp
      dqz(1:neembond) = 0.0_dp
      if (lSandM) then
        dchix(1:numat) = 0.0_dp
        dchiy(1:numat) = 0.0_dp
        dchiz(1:numat) = 0.0_dp
        dchis(1:6) = 0.0_dp
      endif
!
!  Loop over other atoms and build daij/d(alpha)
!
      do jloc = 1,natomsonnode
        j = node2atom(jloc)
        jeem = neemfullrptr(j)
        if (i.ne.j.or.lstrain) then
          qljj = qf(j)
          qlj = qljj*occuf(j)
          nj = nat(j)
!
!  If QEq work out principal quantum number
!
          if (lqeq) then
            if (lmultiqrange) then
              nqrj = nqrnow(jeem)
            else
              nqrj = 1
            endif
            if (nj.le.2) then
              npqnj = 1
            elseif (nj.le.10) then
              npqnj = 2
            elseif (nj.le.18) then
              npqnj = 3
            elseif (nj.le.36) then
              npqnj = 4
            elseif (nj.le.54) then
              npqnj = 5
            elseif (nj.le.86) then
              npqnj = 6
            else
              npqnj = 7
            endif
            zetaj = 0.5_dp*qeqlambda*(2*npqnj+1)/radrange(nqrj,nj,neemtype)
          elseif (lSandM) then
            if (lmultiqrange) then
              nqrj = nqrnow(jeem)
            else
              nqrj = 1
            endif
            zetaj = zetarange(nqrj,nj,neemtype)
            znucj = znucrange(nqrj,nj,neemtype)
          endif
!
!  Find relative vector between atoms
!
          rx = xclat(j) - xci
          ry = yclat(j) - yci
          rz = zclat(j) - zci
!
!  Calculate Euler-Maclaurin correction to 1-D sum
!
          if (ndim.eq.1) then
            qme = 0.0_dp
            dqme(1) = 0.0_dp
            dqme(2) = 0.0_dp
            dqme(3) = 0.0_dp
            call qmatrix1D(rx,ry,rz,.true.,.false.,qme,dqme,d2qme)
            dqme(1) = dqme(1)*angstoev
            dqme(2) = dqme(2)*angstoev
            dqme(3) = dqme(3)*angstoev
            do nbi = 1,nbonds(i)
              nbqo = nbondqb(nbi,i)
              nbq = abs(nbqo)
              if (nbqo.gt.0) then
                dqx(nbq) = dqx(nbq) - qljj*dqme(1)
                dqy(nbq) = dqy(nbq) - qljj*dqme(2)
                dqz(nbq) = dqz(nbq) - qljj*dqme(3)
              elseif (nbqo.lt.0) then
                dqx(nbq) = dqx(nbq) + qljj*dqme(1)
                dqy(nbq) = dqy(nbq) + qljj*dqme(2)
                dqz(nbq) = dqz(nbq) + qljj*dqme(3)
              endif
            enddo
            do nbj = 1,nbonds(j)
              nbqo = nbondqb(nbj,j)
              nbq = abs(nbqo)
              if (nbqo.gt.0) then
                dqx(nbq) = dqx(nbq) - qlii*dqme(1)
                dqy(nbq) = dqy(nbq) - qlii*dqme(2)
                dqz(nbq) = dqz(nbq) - qlii*dqme(3)
              elseif (nbqo.lt.0) then
                dqx(nbq) = dqx(nbq) + qlii*dqme(1)
                dqy(nbq) = dqy(nbq) + qlii*dqme(2)
                dqz(nbq) = dqz(nbq) + qlii*dqme(3)
              endif
            enddo
          endif
!
!  Calculate distances for search
!
          rx = rx - (ml1+1)*r1x
          ry2 = ry*ry
          rz2 = rz*rz
!
!  Loop over lattice vectors
!
          do ii = -ml1,ml1
            rx = rx + r1x
            rr2 = rx*rx + ry2 + rz2
            rl = sqrt(rr2)
            if (rr2.gt.cuts2) then
              rrr = 1.0_dp/rl
              if (lqeq.and.rl.lt.rqeq) then
!***************
!  QEq scheme  *
!***************
                call gammas(npqni,npqnj,zetai,zetaj,rl,gam,dgam,dzetai,dzetaj, &
                            d2zetaii,d2zetaij,d2zetajj,d2zetari,d2zetarj,d2gamr2)
                dtrm1 = dgam*rrr
                dtrm1i = dtrm1
                dtrm1j = dtrm1
              elseif (lSandM.and.rl.lt.rqeq) then
!****************************
!  Streitz-Mintmire scheme  *
!****************************
                call gammasm(zetai,zetaj,rl,gam,dgam,d2gamr2,gamifj,gamjfi,dgamifj,dgamjfi,d2gamifj,d2gamjfi)
                dtrm1 = (dgam - rrr*rrr)*rrr
                dtrm1i = dtrm1
                dtrm1j = dtrm1
                dtrm1zn = znucj*(dgamjfi - dgam)*rrr
              else
!***************
!  EEM scheme  *
!***************
                dtrm1 = - rrr*rrr*rrr
                dtrm1i = dtrm1
                dtrm1j = dtrm1
                dtrm1zn = 0.0_dp
              endif
!
!  First derivatives of matrix elements
!
              dtrm1i = dtrm1i*angstoev
              dtrm1j = dtrm1j*angstoev
              if (lstrain) then
                call real1strterm(ndim,rx,ry,rz,0.0_dp,0.0_dp,0.0_dp,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
              endif
              if (lSandM) then
                dtrm1zn = dtrm1zn*angstoev
                if (lelementOK(nat(i)).and.nregionno(nsft+nrelf2a(i)).eq.1) then
                  dchix(ieem) = dchix(ieem) - dtrm1zn*rx
                  dchiy(ieem) = dchiy(ieem) - dtrm1zn*ry
                  dchiz(ieem) = dchiz(ieem) - dtrm1zn*rz
                endif
                if (lelementOK(nat(j)).and.nregionno(nsft+nrelf2a(j)).eq.1) then
                  dchix(jeem) = dchix(jeem) + dtrm1zn*rx
                  dchiy(jeem) = dchiy(jeem) + dtrm1zn*ry
                  dchiz(jeem) = dchiz(jeem) + dtrm1zn*rz
                endif
                if (lstrain) then
                  dchis(1) = dchis(1) + dtrm1zn*dr2ds(1)
                endif
              endif
              do nbi = 1,nbonds(i)
                nbqo = nbondqb(nbi,i)
                nbq = abs(nbqo)
                if (nbqo.gt.0) then
                  dqx(nbq) = dqx(nbq) - qljj*dtrm1i*rx
                  dqy(nbq) = dqy(nbq) - qljj*dtrm1i*ry
                  dqz(nbq) = dqz(nbq) - qljj*dtrm1i*rz
                  if (lstrain) then
                    dqs(nbq,1) = dqs(nbq,1) + qljj*dtrm1i*dr2ds(1)
                  endif
                elseif (nbqo.lt.0) then
                  dqx(nbq) = dqx(nbq) + qljj*dtrm1i*rx
                  dqy(nbq) = dqy(nbq) + qljj*dtrm1i*ry
                  dqz(nbq) = dqz(nbq) + qljj*dtrm1i*rz
                  if (lstrain) then
                    dqs(nbq,1) = dqs(nbq,1) - qljj*dtrm1i*dr2ds(1)
                  endif
                endif
              enddo
              do nbj = 1,nbonds(j)
                nbqo = nbondqb(nbj,j)
                nbq = abs(nbqo)
                if (nbqo.gt.0) then
                  dqx(nbq) = dqx(nbq) - qlii*dtrm1j*rx
                  dqy(nbq) = dqy(nbq) - qlii*dtrm1j*ry
                  dqz(nbq) = dqz(nbq) - qlii*dtrm1j*rz
                elseif (nbqo.lt.0) then
                  dqx(nbq) = dqx(nbq) + qlii*dtrm1j*rx
                  dqy(nbq) = dqy(nbq) + qlii*dtrm1j*ry
                  dqz(nbq) = dqz(nbq) + qlii*dtrm1j*rz
                endif
              enddo
            endif
          enddo
!       
!  End of loop over lattice vectors
!
        endif
      enddo
!
!  Global sum of dqx, dqy, dqz
!
      dqtmp1(1:neembond) = dqx(1:neembond)
      dqtmp1(neembond+1:2*neembond) = dqy(1:neembond)
      dqtmp1(2*neembond+1:3*neembond) = dqz(1:neembond)
!
      call sumall(dqtmp1,dqtmp2,3_i4*neembond,"dchargesplitd","dqtmp")
!
      dqx(1:neembond) = dqtmp2(1:neembond)
      dqy(1:neembond) = dqtmp2(neembond+1:2*neembond)
      dqz(1:neembond) = dqtmp2(2*neembond+1:3*neembond)
!
!  Multiply dqx/dqy/dqz by inverse matrix 
!
      do jl = 1,nbloc
        sumx = 0.0_dp
        sumy = 0.0_dp
        sumz = 0.0_dp
        do k = 1,neembond
          sumx = sumx + dqx(k)*potl(k,jl)
          sumy = sumy + dqy(k)*potl(k,jl)
          sumz = sumz + dqz(k)*potl(k,jl)
        enddo
        dqbonddxyz(ind+1,jl) = dqbonddxyz(ind+1,jl) - sumx
        dqbonddxyz(ind+2,jl) = dqbonddxyz(ind+2,jl) - sumy
        dqbonddxyz(ind+3,jl) = dqbonddxyz(ind+3,jl) - sumz
      enddo
      if (lSandM) then
!
!  Global sum of dchix, dchiy, dchiz
!
        dqtmp1(1:numat) = dchix(1:numat)
        dqtmp1(numat+1:2*numat) = dchiy(1:numat)
        dqtmp1(2*numat+1:3*numat) = dchiz(1:numat)
!
        if (lstrain) then
          dqtmp1(3*numat+1:3*numat+nstrains) = dchis(1:nstrains)
          call sumall(dqtmp1,dqtmp2,3_i4*numat+nstrains,"dchargesplitd","dchitmp")
        else
          call sumall(dqtmp1,dqtmp2,3_i4*numat,"dchargesplitd","dchitmp")
        endif
!
        dchix(1:numat) = dqtmp2(1:numat)
        dchiy(1:numat) = dqtmp2(numat+1:2*numat)
        dchiz(1:numat) = dqtmp2(2*numat+1:3*numat)
!
        if (lstrain) then
          dchis(1:nstrains) = dqtmp2(3*numat+1:3*numat+nstrains)
        endif
!
!  For S & M, multiply dchix/dchiy/dchiz by inverse matrix
!
        do nbi = 1,nbonds(i)
          nbqo = nbondqb(nbi,i)
          nbq = abs(nbqo)
          if (nbqo.gt.0) then
            do jl = 1,nbloc
              do k = 1,neemfull
                indj = 3*(neemfullptr(k) - 1)
                dqbonddxyz(indj+1,jl) = dqbonddxyz(indj+1,jl) - dchix(k)*potl(nbq,jl)
                dqbonddxyz(indj+2,jl) = dqbonddxyz(indj+2,jl) - dchiy(k)*potl(nbq,jl)
                dqbonddxyz(indj+3,jl) = dqbonddxyz(indj+3,jl) - dchiz(k)*potl(nbq,jl)
              enddo
            enddo
          elseif (nbqo.lt.0) then
            do jl = 1,nbloc
              do k = 1,neemfull
                indj = 3*(neemfullptr(k) - 1)
                dqbonddxyz(indj+1,jl) = dqbonddxyz(indj+1,jl) + dchix(k)*potl(nbq,jl)
                dqbonddxyz(indj+2,jl) = dqbonddxyz(indj+2,jl) + dchiy(k)*potl(nbq,jl)
                dqbonddxyz(indj+3,jl) = dqbonddxyz(indj+3,jl) + dchiz(k)*potl(nbq,jl)
              enddo
            enddo
          endif
        enddo
      endif
      if (lSandM.and.lstrain) then
!
!  For S & M, multiply dchis by inverse matrix
!
        do nbi = 1,nbonds(i)
          nbqo = nbondqb(nbi,i)
          nbq = abs(nbqo)
          if (nbqo.gt.0) then
            do jl = 1,nbloc
              do kl = 1,nstrains
                dqbondds(kl,jl) = dqbondds(kl,jl) - dchis(kl)*potl(nbq,jl)
              enddo
            enddo
          elseif (nbqo.lt.0) then
            do jl = 1,nbloc
              do kl = 1,nstrains
                dqbondds(kl,jl) = dqbondds(kl,jl) + dchis(kl)*potl(nbq,jl)
              enddo
            enddo
          endif
        enddo
      endif
!
!  Outer loop over atoms
!
    enddo
    if (lstrain) then
!
!  Global sum of dqs
!
      do kl = 1,nstrains
        do k = 1,neembond
          dqtmp1((kl-1)*neembond+k) = dqs(k,kl)
        enddo
      enddo
!
      call sumall(dqtmp1,dqtmp2,nstrains*neembond,"dchargesplitd","dqstmp")
!
      do kl = 1,nstrains
        do k = 1,neembond
          dqs(k,kl) = dqtmp2((kl-1)*neembond+k)
        enddo
      enddo
!         
!  Strain terms
!           
      do jl = 1,nbloc
        do kl = 1,nstrains
          sum = 0.0_dp
          do k = 1,neembond
            sum = sum + dqs(k,kl)*potl(k,jl)
          enddo
          dqbondds(kl,jl) = dqbondds(kl,jl) - sum
        enddo
      enddo
    endif
!******************************
!  End of cluster / 1-D case  *
!******************************
  endif
135 continue
!******************************************************
!  Electric field contribution to charge derivatives  *
!******************************************************
  if (lfieldcfg(ncf).or.(ntdfieldcfg(ncf).gt.0)) then
!
!  For electric field, multiply dchix/dchiy/dchiz by inverse matrix
!
    do il = 1,nbloc
      do j = 1,numat
        if (lelementOK(nat(j)).and.nregionno(nsft+nrelf2a(j)).eq.1) then
          indj = 3*(neemfullptr(j) - 1)
          do nbj = 1,nbonds(j)
            nbqjo = nbondqb(nbj,j)
            nbqj = abs(nbqjo)
            if (nbqjo.gt.0) then
              dqbonddxyz(indj+1,il) = dqbonddxyz(indj+1,il) - fieldx*potl(nbqj,il)
              dqbonddxyz(indj+2,il) = dqbonddxyz(indj+2,il) - fieldy*potl(nbqj,il)
              dqbonddxyz(indj+3,il) = dqbonddxyz(indj+3,il) - fieldz*potl(nbqj,il)
            elseif (nbqjo.lt.0) then
              dqbonddxyz(indj+1,il) = dqbonddxyz(indj+1,il) + fieldx*potl(nbqj,il)
              dqbonddxyz(indj+2,il) = dqbonddxyz(indj+2,il) + fieldy*potl(nbqj,il)
              dqbonddxyz(indj+3,il) = dqbonddxyz(indj+3,il) + fieldz*potl(nbqj,il)
            endif
          enddo
        endif
      enddo
    enddo
  endif
!************************************************************
!  Convert bond charge derivatives into charge derivatives  *
!************************************************************
  allocate(nbtmp(neembond+1),stat=status)
  if (status/=0) call outofmemory('dchargesplitd','nbtmp')
  allocate(dqbloc(3*numat,nbloc+nblocksize),stat=status)
  if (status/=0) call outofmemory('dchargesplitd','dqbloc')
  if (lstrain) then
    allocate(dqbsloc(nstrains,nbloc+nblocksize),stat=status)
    if (status/=0) call outofmemory('dchargesplitd','dqbsloc')
  endif
!
  nbtmp(1:neembond+1) = 0
!
!  For now the sharing of split bond charge derivatives is handled in the simplest way
!  by looping over nodes and broadcasting.
!
  do inode = 0,nprocs-1
!
!  Broadcast number of elements on this node and pointer
!
    if (procid.eq.inode) then
      nbtmp(1) = nbloc
      nbtmp(2:nbloc+1) = nblocptr(1:nbloc)
    endif
    call isendall(nbtmp,neembond+1_i4,inode,"dchargesplitd","nbtmp")
!
!  Set number of elements for current node
!
    nbinode = nbtmp(1)
!
!  Send local part of split bond charge derivatives
!
    if (procid.eq.inode) then
      do il = 1,nbloc
        do j = 1,3*numat
          dqbloc(j,il) = dqbonddxyz(j,il)
        enddo
      enddo
    endif
    call sendall(dqbloc,3_i4*numat*nbinode,inode,"dchargesplitd","dqbloc")
!
    if (lstrain) then
      if (procid.eq.inode) then
        do il = 1,nbloc
          do kl = 1,nstrains
            dqbsloc(kl,il) = dqbondds(kl,il)
          enddo
        enddo
      endif
      call sendall(dqbsloc,nstrains*nbinode,inode,"dchargesplitd","dqbsloc")
    endif
!
!  Populate local elements of dqdxyz
!
    do il = 1,nbinode
      i = nbtmp(1+il)
      ib1 = atom2local(neembonded(1,i))
      ib2 = atom2local(neembonded(2,i))
      if (ib1.gt.0) then
        do j = 1,3*numat
          dqdxyz(j,ib1) = dqdxyz(j,ib1) + dqbloc(j,il)
        enddo
      endif
      if (ib2.gt.0) then
        do j = 1,3*numat
          dqdxyz(j,ib2) = dqdxyz(j,ib2) - dqbloc(j,il)
        enddo
      endif
    enddo
    if (lstrain) then
!
!  Populate local elements of dqds
!
      do il = 1,nbinode
        i = nbtmp(1+il)
        ib1 = atom2local(neembonded(1,i))
        ib2 = atom2local(neembonded(2,i))
        if (ib1.gt.0) then
          do kl = 1,nstrains
            dqds(kl,ib1) = dqds(kl,ib1) + dqbsloc(kl,il)
          enddo
        endif
        if (ib2.gt.0) then
          do kl = 1,nstrains
            dqds(kl,ib2) = dqds(kl,ib2) - dqbsloc(kl,il)
          enddo
        endif
      enddo
    endif
!
!  End of loop over nodes
!
  enddo
!
  if (lstrain) then
    deallocate(dqbsloc,stat=status)
    if (status/=0) call deallocate_error('dchargesplitd','dqbsloc')
  endif
  deallocate(dqbloc,stat=status)
  if (status/=0) call deallocate_error('dchargesplitd','dqbloc')
  deallocate(nbtmp,stat=status)
  if (status/=0) call deallocate_error('dchargesplitd','nbtmp')
!***********************************************************************************
!  Enforce sum rules that total charge derivative equals zero for each coordinate  *
!***********************************************************************************
!
!  Sum elements and count number of active atoms
!
  dqx(1:numat) = 0.0_dp
  dqy(1:numat) = 0.0_dp
  dqz(1:numat) = 0.0_dp
!
  do iloc = 1,natomsonnode
    jx = - 2
    jy = - 1
    jz =   0
    do j = 1,numat
      jx = jx + 3
      jy = jy + 3
      jz = jz + 3
      dqx(j) = dqx(j) + dqdxyz(jx,iloc)
      dqy(j) = dqy(j) + dqdxyz(jy,iloc)
      dqz(j) = dqz(j) + dqdxyz(jz,iloc)
    enddo
  enddo
!
!  Globalise dqx / dqy / dqz 
!
  call sumall(dqx,dqtmp1,numat,"dcharged","dqx")
  dqx(1:numat) = dqtmp1(1:numat)
  call sumall(dqy,dqtmp1,numat,"dcharged","dqx")
  dqy(1:numat) = dqtmp1(1:numat)
  call sumall(dqz,dqtmp1,numat,"dcharged","dqx")
  dqz(1:numat) = dqtmp1(1:numat)
!
!  Average error and subtract from active elements
!
  if (neemfull.gt.0) then
    do i = 1,numat
      dqx(i) = - dqx(i)/dble(neemfull)
      dqy(i) = - dqy(i)/dble(neemfull)
      dqz(i) = - dqz(i)/dble(neemfull)
    enddo
    do iloc = 1,natomsonnode
      i = node2atom(iloc)
      if (lelementOK(nat(i)).and.nregionno(nsft+nrelf2a(i)).eq.1) then
        jx = - 2
        jy = - 1
        jz =   0
        do j = 1,numat
          jx = jx + 3
          jy = jy + 3
          jz = jz + 3
          dqdxyz(jx,iloc) = dqdxyz(jx,iloc) + dqx(j)
          dqdxyz(jy,iloc) = dqdxyz(jy,iloc) + dqy(j)
          dqdxyz(jz,iloc) = dqdxyz(jz,iloc) + dqz(j)
        enddo
      endif
    enddo
  endif
!
!  Free local memory
!
  if (lstrain) then
    deallocate(dqs,stat=status)
    if (status/=0) call deallocate_error('dchargesplitd','dqs')
  endif
  if (lSandM) then
    deallocate(dchiz,stat=status)
    if (status/=0) call deallocate_error('dchargesplitd','dchiz')
    deallocate(dchiy,stat=status)
    if (status/=0) call deallocate_error('dchargesplitd','dchiy')
    deallocate(dchix,stat=status)
    if (status/=0) call deallocate_error('dchargesplitd','dchix')
  endif
  deallocate(dqtmp2,stat=status)
  if (status/=0) call deallocate_error('dchargesplitd','dqtmp2')
  deallocate(dqtmp1,stat=status)
  if (status/=0) call deallocate_error('dchargesplitd','dqtmp1')
  deallocate(dqz,stat=status)
  if (status/=0) call deallocate_error('dchargesplitd','dqz')
  deallocate(dqy,stat=status)
  if (status/=0) call deallocate_error('dchargesplitd','dqy')
  deallocate(dqx,stat=status)
  if (status/=0) call deallocate_error('dchargesplitd','dqx')
  deallocate(neemfullrptr,stat=status)
  if (status/=0) call deallocate_error('dchargesplitd','neemfullrptr')
  deallocate(neemfullptr,stat=status)
  if (status/=0) call deallocate_error('dchargesplitd','neemfullptr')
!
  if (lstrain) then
    deallocate(dqbondds,stat=status)
    if (status/=0) call deallocate_error('dchargesplitd','dqbondds')
  endif
  deallocate(dqbonddxyz,stat=status)
  if (status/=0) call deallocate_error('dchargesplitd','dqbonddxyz')
!
  if (lprint) then
    call mpbarrier
    if (ndim.eq.0) then
      if (ioproc) then
        write(ioout,'(/,''  First derivatives of charge distribution : (Angstroms**-1)'',/)')
      endif
    else
      if (ioproc) then
        write(ioout,'(/,''  First derivatives of charge distribution : '',/)')
        write(ioout,'(''  Strain :'',/)')
        write(ioout,'(''  Atom   '',6(i10))') (j,j=1,nstrains)
      endif
      if (lioproconly) then
!
!  Allocate temporary workspace for communication
!
        ntmp = nstrains
        ntag = 1
        allocate(dqtmp2(ntmp),stat=status)
        if (status/=0) call outofmemory('dcharged','dqtmp2')
        allocate(StatMPI(MPI_Status_Size),stat=status)
        if (status/=0) call outofmemory('dcharged','StatMPI')
      endif
      call mpbarrier
      do i = 1,numat
        iloc = atom2local(i)
        if (lioproconly.and.atom2node(i).ne.0_i4) then
!
!  Post receive
!
          if (ioproc) then
            nnode = atom2node(i)
            call MPI_IRecv(dqtmp2,ntmp,MPI_double_precision,nnode, &
                           ntag,MPI_Comm_World,Request,MPIerror)
          endif
!
!  Pass data to ioproc for writing
!
          if (iloc.gt.0) then
            dqtmp2(1:nstrains) = dqds(1:nstrains,iloc)
!
!  Post send
!
            call MPI_ISend(dqtmp2,ntmp,MPI_double_precision,0, &
                           ntag,MPI_Comm_World,Request,MPIerror)
          endif
          if (ioproc.or.iloc.gt.0) then
            call MPI_WaitAll(1,Request,StatMPI,MPIerror)
          endif
          if (ioproc) then
!
!  Write on I/O node
!
            write(ioout,'(i6,4x,6f10.6)') i,(dqtmp2(j),j=1,nstrains)
          endif
        else
          if (iloc.gt.0) then
            write(ioout,'(i6,4x,6f10.6)') i,(dqds(j,atom2local(i)),j=1,nstrains)
          endif
        endif
        call mpbarrier
      enddo
      if (lioproconly) then
        deallocate(StatMPI,stat=status)
        if (status/=0) call deallocate_error('dcharged','StatMPI')
        deallocate(dqtmp2,stat=status)
        if (status/=0) call deallocate_error('dcharged','dqtmp2')
      endif
      if (ioproc) then
        write(ioout,'(/,''  Coordinate (Angstroms**-1) :'',/)')
      endif
    endif
!
    if (lioproconly) then
!
!  Allocate temporary workspace for communication
!
      ntmp = 3*numat
      ntag = 1
      allocate(dqtmp2(ntmp),stat=status)
      if (status/=0) call outofmemory('dcharged','dqtmp2')
      allocate(StatMPI(MPI_Status_Size),stat=status)
      if (status/=0) call outofmemory('dcharged','StatMPI')
    endif
    call mpbarrier
    do i = 1,numat
      iloc = atom2local(i)
      if (lioproconly.and.atom2node(i).ne.0_i4) then
!
!  Post receive
!
        if (ioproc) then
          nnode = atom2node(i)
          call MPI_IRecv(dqtmp2,ntmp,MPI_double_precision,nnode, &
                         ntag,MPI_Comm_World,Request,MPIerror)
        endif
!
!  Pass data to ioproc for writing
!
        if (iloc.gt.0) then
          dqtmp2(1:3*numat) = dqdxyz(1:3*numat,iloc)
!
!  Post send
!
          call MPI_ISend(dqtmp2,ntmp,MPI_double_precision,0, &
                         ntag,MPI_Comm_World,Request,MPIerror)
        endif
        if (ioproc.or.iloc.gt.0) then
          call MPI_WaitAll(1,Request,StatMPI,MPIerror)
        endif
        if (ioproc) then
!
!  Write on I/O node
!
          write(ioout,'(''  Atom   '',i10)') i
          indj = 0
          do j = 1,numat
            write(ioout,'(i6,'' x '',4x,f10.6)') j,dqtmp2(indj+1)
            write(ioout,'(i6,'' y '',4x,f10.6)') j,dqtmp2(indj+2)
            write(ioout,'(i6,'' z '',4x,f10.6)') j,dqtmp2(indj+3)
            indj = indj + 3
          enddo
          write(ioout,'(/)')
        endif
      else
        if (iloc.gt.0) then
          write(ioout,'(''  Atom   '',i10)') i
          indj = 0
          do j = 1,numat
            write(ioout,'(i6,'' x '',4x,f10.6)') j,dqdxyz(indj+1,iloc)
            write(ioout,'(i6,'' y '',4x,f10.6)') j,dqdxyz(indj+2,iloc)
            write(ioout,'(i6,'' z '',4x,f10.6)') j,dqdxyz(indj+3,iloc)
            indj = indj + 3
          enddo
          write(ioout,'(/)')
        endif
      endif
      call mpbarrier
    enddo
    if (lioproconly) then
      deallocate(StatMPI,stat=status)
      if (status/=0) call deallocate_error('dcharged','StatMPI')
      deallocate(dqtmp2,stat=status)
      if (status/=0) call deallocate_error('dcharged','dqtmp2')
    endif
  endif
#ifdef TRACE
  call trace_out('dchargesplitd')
#endif
#else
  call outerror('dchargesplitd called when not compiled with MPI',0_i4)
  call stopnow('dchargesplitd')
#endif
!
  return
  end
