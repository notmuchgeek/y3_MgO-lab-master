  subroutine dcharged(lprint,emat,maxe,lrecalcA,lstrainin)
!
!  Calculates the first derivative of the charge with respect
!  to the coordinates of the atoms.
!
!  NB: Might be problems for fixed charge case!!
!
!  Distributed memory parallel version.
!
!   2/17 Created from dcharge
!   7/17 Further work on parallelisation 
!   8/17 Routine wrapped in ifdef since it won't work or be called without MPI
!   2/18 Trace added
!   4/18 Parallel I/O corrected
!   5/18 Local strain flag set based on whether nstrains is non-zero
!   5/18 Check on maxd2u corrected to allow for parallel distribution
!   5/18 Multiple qranges added
!   7/18 Array for equations now passed in
!   7/18 dqs now made into a local temporary array
!   9/18 Strain module added
!  11/18 Finite strain flag introduced instead of lstraincell
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!  12/19 Rigid molecule modifications added
!   2/21 Correction to cluster case for setting of jeem
!
!  On entry:
!
!    emat      = inverse matrix calculated during EEM/QEq
!    lprint    = logical indicating whether printed output is wanted
!    lrecalcA  = if .true. this forces the recalculation of A
!    lstrainin = if .true. then strain derivatives are calculated
!                if dimensionality is greater than 0
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
!  Copyright Curtin University 2021
!
!  Julian Gale, CIC, Curtin University, February 2021
!
#ifdef MPI
  use configurations, only : nregionno
#endif
  use g_constants
  use control
  use current
  use derivatives
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
  integer(i4), intent(in)                      :: maxe
  real(dp),    intent(inout)                   :: emat(maxe,*)
  logical,     intent(in)                      :: lprint
  logical,     intent(in)                      :: lrecalcA
  logical,     intent(in)                      :: lstrainin
#ifdef MPI
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ieem
  integer(i4)                                  :: ieemf
  integer(i4)                                  :: ifail
  integer(i4)                                  :: ii
  integer(i4)                                  :: iloc
  integer(i4)                                  :: ind
  integer(i4)                                  :: indj
  integer(i4)                                  :: is
  integer(i4)                                  :: iv
  integer(i4)                                  :: j
  integer(i4)                                  :: jeem
  integer(i4)                                  :: jj
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
  integer(i4)                                  :: neemfull
  integer(i4)                                  :: neemfullloc
  integer(i4), dimension(:), allocatable       :: neemfullptr
  integer(i4), dimension(:), allocatable       :: neemfullrptr
  integer(i4)                                  :: ni
  integer(i4)                                  :: nj
  integer(i4)                                  :: npqni
  integer(i4)                                  :: npqnj
  integer(i4)                                  :: nqr
  integer(i4)                                  :: nqrj
  integer(i4)                                  :: ns
  integer(i4)                                  :: status
!
  integer                                      :: MPIerror
  integer                                      :: idesd(9)
  integer                                      :: idesq(9)
  integer                                      :: idess(9)
  integer                                      :: idest(9)
  integer                                      :: ifails
  integer                                      :: ld
  integer                                      :: nb
  integer                                      :: ncs
  integer                                      :: ntag
  integer                                      :: nnode
  integer                                      :: ntmp
  integer                                      :: Request
  integer(i4), dimension(:), allocatable       :: StatMPI       ! Array for status from MPI
!
  logical                                      :: lfound
  logical                                      :: lhi
  logical                                      :: lhj
  logical                                      :: lhpresent
  logical                                      :: lstrain
  real(dp)                                     :: accf2
  real(dp)                                     :: arg
  real(dp)                                     :: argtest
  real(dp),    dimension(:), allocatable       :: chitmp
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
  real(dp),    dimension(:,:), allocatable     :: dchis
  real(dp),    dimension(:),   allocatable     :: dchix
  real(dp),    dimension(:),   allocatable     :: dchiy
  real(dp),    dimension(:),   allocatable     :: dchiz
  real(dp),    dimension(:,:), allocatable     :: dsmtrm
  real(dp),    dimension(:),   allocatable     :: dtmp
  real(dp),    dimension(:),   allocatable     :: dtmp2
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
  real(dp),    dimension(:,:), allocatable     :: dqtmp
  real(dp),    dimension(:,:), allocatable     :: dqs
  real(dp),    dimension(:,:), allocatable     :: dqx
  real(dp),    dimension(:,:), allocatable     :: dqy
  real(dp),    dimension(:,:), allocatable     :: dqz
  real(dp),    dimension(:,:), allocatable     :: dstmp
  real(dp)                                     :: dtrm1
  real(dp)                                     :: dtrm1i
  real(dp)                                     :: dtrm1zn
  real(dp)                                     :: dtrm1j
  real(dp)                                     :: dzetai
  real(dp)                                     :: dzetaj
  real(dp)                                     :: dr2ds(6)
  real(dp)                                     :: d2r2dx2(3,3)
  real(dp)                                     :: d2r2ds2(6,6)
  real(dp)                                     :: d2r2dsdx(6,3)
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
  real(dp)                                     :: rjfac
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
  real(dp)                                     :: trmi
  real(dp)                                     :: xci
  real(dp)                                     :: yci
  real(dp)                                     :: zci
  real(dp)                                     :: xd
  real(dp)                                     :: yd
  real(dp)                                     :: zd
  real(dp)                                     :: zetah0
  real(dp)                                     :: zetai
  real(dp)                                     :: zetaj
  real(dp)                                     :: znucj
  real(dp)                                     :: ztrm1
#ifdef TRACE
  call trace_in('dcharged')
#endif
!
!  Set local strain flag
!
  lstrain = (lstrainin.and.nstrains.gt.0)
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
  if (ndim.gt.0) then
    do i = 1,natomsonnode
      do j = 1,nstrains
        dqds(j,i) = 0.0_dp
      enddo
    enddo
  endif
!
!  For S and M there is a need to store a large array for terms
!
  if (lSandM) then
    allocate(dsmtrm(3*numat,numat),stat=status)
    if (status/=0) call outofmemory('dcharged','dsmtrm')
    allocate(dtmp(3*numat),stat=status)
    if (status/=0) call outofmemory('dcharged','dtmp')
!
    dsmtrm(1:3*numat,1:numat) = 0.0_dp
  endif
!
!  If electrostatics have been turned off then there is no point in proceeding.
!
  if (.not.lDoElectrostatics) then
#ifdef TRACE
    call trace_out('dcharged')
#endif
    return
  endif
!
  cuts2 = cuts*cuts
  rconv = 1.0_dp/autoangs
!
!  Is this QEq with H present?
!
  lhpresent = .false.
  if (lqeq) then
    i = 0
    do while (.not.lhpresent.and.i.lt.numat)
      i = i + 1
      if (nat(i).eq.1) lhpresent = .true.
    enddo
  endif
!
!  Find number of EEM active atoms
!
  allocate(neemfullptr(numat),stat=status)
  if (status/=0) call outofmemory('dcharged','neemfullptr')
  allocate(neemfullrptr(numat),stat=status)
  if (status/=0) call outofmemory('dcharged','neemfullrptr')
!
  neemfullrptr(1:numat) = 0
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
  neemfullloc = 0
  do iloc = 1,natomsonnode
    i = node2atom(iloc)
    if (lelementOK(nat(i)).and.nregionno(nsft+nrelf2a(i)).eq.1) then
      neemfullloc = neemfullloc + 1
    endif
  enddo
!
  if (lfieldcfg(ncf).or.(ntdfieldcfg(ncf).gt.0)) then
!
!  Get electric field if needed
!
    call electricfieldparts(fieldx,fieldy,fieldz)
  endif
  if ((lqeq.and.lhpresent).or.lrecalcA) then
!*************************************
!  Hydrogen/QEq case or Symopt case  *
!*************************************
    allocate(chitmp(numat),stat=status)
    if (status/=0) call outofmemory('dcharged','chitmp')
!
!  If hydrogen is present then we need to calculate A + (dA/dq).q
!  invert this matrix to use in place of A**-1 below.
!
!  First generate 2 centre terms
!
    chitmp(1:numat) = 0.0_dp
    if (lnoqeem) then
      emat(1:numat,1:natomsonnode) = 0.0_dp
    else
      if (lqeq.and.lhpresent) then
        call genpotdm(emat,maxe,chitmp,2_i4)
      else
        call genpotdm(emat,maxe,chitmp,1_i4)
      endif
    endif
!
!  Reduce potential contributions to neemfull x neemfull
!
    ieem = 0
    do iloc = 1,natomsonnode
      i = node2atom(iloc)
      if (lelementOK(nat(i)).and.nregionno(nsft+nrelf2a(i)).eq.1) then
        ieem = ieem + 1
        jeem = 0
        do j = 1,numat
          if (lelementOK(nat(j)).and.nregionno(nsft+nrelf2a(j)).eq.1) then
            jeem = jeem + 1
            emat(jeem,ieem) = emat(j,iloc)
          endif
        enddo
      endif
    enddo
!
!  Add one centre terms
!
    if (lqeq) then
      ieem = 0
      do iloc = 1,natomsonnode
        i = node2atom(iloc)
        ni = nat(i)
        if (lelementOK(ni).and.nregionno(nsft+nrelf2a(i)).eq.1) then
          ieem = ieem + 1
          ieemf = neemfullrptr(i)
          if (lmultiqrange) then
            nqr = nqrnow(ieemf)
          else
            nqr = 1
          endif
          if (nat(i).ne.1) then
            emat(ieemf,ieem) = emat(ieemf,ieem) + 2.0_dp*murange(nqr,ni,neemtype)*occuf(i)
          else
!
!  For hydrogen charge dependant factor must be introduced
!  which includes both Aii and dAii/dqi
!
            zetah0 = 0.529177_dp*0.75_dp/radrange(nqr,1,neemtype)
            rjfac = 1.0_dp + (2.0_dp*qf(i)/zetah0)
            emat(ieemf,ieem) = emat(ieemf,ieem) + 2.0_dp*murange(nqr,1,neemtype)*occuf(i)*rjfac
          endif 
        endif 
      enddo
    else
      ieem = 0
      do iloc = 1,natomsonnode
        i = node2atom(iloc)
        ni = nat(i)
        if (lelementOK(ni).and.nregionno(nsft+nrelf2a(i)).eq.1) then
          ieem = ieem + 1
          ieemf = neemfullrptr(i)
          if (lmultiqrange) then
            nqr = nqrnow(ieemf)
          else
            nqr = 1
          endif
          emat(ieemf,ieem) = emat(ieemf,ieem) + 2.0_dp*murange(nqr,ni,neemtype)*occuf(i)
        endif
      enddo
    endif
!
!  Complete constraint terms
!
    ieem = 0
    do iloc = 1,natomsonnode
      i = node2atom(iloc)
      if (lelementOK(nat(i)).and.nregionno(nsft+nrelf2a(i)).eq.1) then
        ieem = ieem + 1
        ieemf = neemfullrptr(i)
        emat(ieemf,neemfullloc+1) = occuf(i)
        emat(neemfull+1,ieem) = 1.0_dp
      endif
    enddo
    emat(neemfull+1,neemfullloc+1) = 0.0_dp
!     
!  Matrix inversion 
!     
    n = neemfull + 1
    call matrix_inversion_library(n,1_i4,maxe,nblocksize,emat,0_i4,ifail)
!
!  Was inversion successful?
!
    if (ifail.ne.0) then
      call outerror('matrix inversion failed in dcharged',0_i4)
      call stopnow('dcharged')
    endif
!
    deallocate(chitmp,stat=status)
    if (status/=0) call deallocate_error('dcharged','chitmp')
!**************************
!  End hydrogen/QEq case  *
!**************************
  endif
  rqeq2 = rqeq*rqeq
!
!  Allocate local memory
!
  allocate(dqtmp(max(6_i4,numat),max(2_i4,natomsonnode)),stat=status)
  if (status/=0) call outofmemory('dcharged','dqtmp')
  allocate(dqx(numat,natomsonnode),stat=status)
  if (status/=0) call outofmemory('dcharged','dqx')
  allocate(dqy(numat,natomsonnode),stat=status)
  if (status/=0) call outofmemory('dcharged','dqy')
  allocate(dqz(numat,natomsonnode),stat=status)
  if (status/=0) call outofmemory('dcharged','dqz')
  if (lSandM) then
    allocate(dchix(numat),stat=status)
    if (status/=0) call outofmemory('dcharged','dchix')
    allocate(dchiy(numat),stat=status)
    if (status/=0) call outofmemory('dcharged','dchiy')
    allocate(dchiz(numat),stat=status)
    if (status/=0) call outofmemory('dcharged','dchiz')
    if (lstrain) then
      allocate(dchis(6,natomsonnode),stat=status)
      if (status/=0) call outofmemory('dcharged','dchis')
      allocate(dstmp(6,natomsonnode),stat=status)
      if (status/=0) call outofmemory('dcharged','dstmp')
    endif
  endif
!
!  Set up Blacs descriptors for matrices
!
  nb = nblocksize
  ifails = 0
  ncs = numat
  ld = maxe
  call descinit( idesd, ncs, ncs, nb, nb, 0, 0, iBlacsContext, ld, ifails )
!
  if (ifails.ne.0) then
    call outerror('initialisation in descinit failed',0_i4)
    call stopnow('dcharged')
  endif
!
  ifails = 0
  ncs = numat
  ld = numat
  call descinit( idesq, ncs, ncs, nb, nb, 0, 0, iBlacsContext, ld, ifails )
!
  if (ifails.ne.0) then
    call outerror('initialisation in descinit failed',0_i4)
    call stopnow('dcharged')
  endif
!
  ifails = 0
  ncs = numat
  ld = max(6,numat)
  call descinit( idest, ncs, ncs, nb, nb, 0, 0, iBlacsContext, ld, ifails )
!
  if (ifails.ne.0) then
    call outerror('initialisation in descinit failed',0_i4)
    call stopnow('dcharged')
  endif
!
  if (lSandM.and.lstrain) then
    ifails = 0
    ncs = numat
    ld = 6
    call descinit( idess, 6, ncs, nb, nb, 0, 0, iBlacsContext, ld, ifails )
  endif
!
  if (ifails.ne.0) then
    call outerror('initialisation in descinit failed',0_i4)
    call stopnow('dcharged')
  endif
!
  if (lstrain) then
    allocate(dqs(numat,nstrains),stat=status)
    if (status/=0) call outofmemory('dcharged','dqs')
!
!  Zero temporary strain derivative arrays
!
    do i = 1,nstrains
      do j = 1,numat
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
    if (lstrain) then
      call gstrterms(ndim,maxkvec,nkvec,xrk,yrk,zrk,dg2ds,d2g2dx2,d2g2ds2,.false.)
    endif
!
!  Start loop over atoms for coordinate differentiation
!
    ieem = 0
    do iloc = 1,natomsonnode
      i = node2atom(iloc)
      xci = xclat(i)
      yci = yclat(i)
      zci = zclat(i)
      qli = qf(i)
      ind = 3*(i - 1)
      if (lelementOK(nat(i)).and.nregionno(nsft+nrelf2a(i)).eq.1) then
        ieem = ieem + 1
        ieemf = neemfullrptr(i)
      endif
!
!  Zero temporary derivative arrays
!
      do j = 1,numat
        dqx(j,iloc) = 0.0_dp
        dqy(j,iloc) = 0.0_dp
        dqz(j,iloc) = 0.0_dp
      enddo
      jeem = 0
      do j = 1,numat
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
            if (lelementOK(nat(i)).and.nregionno(nsft+nrelf2a(i)).eq.1) then
              dqx(ieemf,iloc) = dqx(ieemf,iloc) + sine(iv)*xrk(iv)*qlj
              dqy(ieemf,iloc) = dqy(ieemf,iloc) + sine(iv)*yrk(iv)*qlj
              dqz(ieemf,iloc) = dqz(ieemf,iloc) + sine(iv)*zrk(iv)*qlj
            endif
            if (lelementOK(nat(j)).and.nregionno(nsft+nrelf2a(j)).eq.1) then
              dqx(jeem,iloc) = dqx(jeem,iloc) + sine(iv)*xrk(iv)*qli
              dqy(jeem,iloc) = dqy(jeem,iloc) + sine(iv)*yrk(iv)*qli
              dqz(jeem,iloc) = dqz(jeem,iloc) + sine(iv)*zrk(iv)*qli
            endif
            if (lstrain) then
              if (lelementOK(nat(i)).and.nregionno(nsft+nrelf2a(i)).eq.1) then
                costrm = cos(argc(iv))*angstoev*qlj
                strm1 = costrm*ktrms(iv)
                strm2 = costrm*ktrm(iv)
                if (lfinitestrain) then
                  do is = 1,nstrains
                    dqs(ieemf,is) = dqs(ieemf,is) + strm1*dg2ds(iv,is) - strm2*strainddetds(is)*straindet
                  enddo
                else
                  do is = 1,nstrains
                    dqs(ieemf,is) = dqs(ieemf,is) + strm1*dg2ds(iv,is)
                  enddo
                  dqs(ieemf,1) = dqs(ieemf,1) - strm2
                  dqs(ieemf,2) = dqs(ieemf,2) - strm2
                  dqs(ieemf,3) = dqs(ieemf,3) - strm2
                endif
              endif
            endif
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
          if (lelementOK(nat(i)).and.nregionno(nsft+nrelf2a(i)).eq.1) then
            dqz(ieemf,iloc) = dqz(ieemf,iloc) - dtrm1*qlj
          endif
          if (lelementOK(nat(j)).and.nregionno(nsft+nrelf2a(j)).eq.1) then
            dqz(jeem,iloc) = dqz(jeem,iloc) - dtrm1*qli
          endif
          if (lstrain) then
            if (lelementOK(nat(i)).and.nregionno(nsft+nrelf2a(i)).eq.1) then
              dqs(ieemf,1) = dqs(ieemf,1) - etrm*qlj
              dqs(ieemf,2) = dqs(ieemf,2) - etrm*qlj
            endif
          endif
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
                if (lelementOK(nat(i)).and.nregionno(nsft+nrelf2a(i)).eq.1) then
                  dqx(ieemf,iloc) = dqx(ieemf,iloc) + sineq*xrk(iv)*qlj
                  dqy(ieemf,iloc) = dqy(ieemf,iloc) + sineq*yrk(iv)*qlj
                  dqz(ieemf,iloc) = dqz(ieemf,iloc) - cosa*ztrm1*qlj
                endif
                if (lelementOK(nat(j)).and.nregionno(nsft+nrelf2a(j)).eq.1) then
                  dqx(jeem,iloc) = dqx(jeem,iloc) + sineq*xrk(iv)*qli
                  dqy(jeem,iloc) = dqy(jeem,iloc) + sineq*yrk(iv)*qli
                  dqz(jeem,iloc) = dqz(jeem,iloc) - cosa*ztrm1*qli
                endif
                if (lstrain) then
                  if (lelementOK(nat(i)).and.nregionno(nsft+nrelf2a(i)).eq.1) then
                    costrm = cosa*angstoev*qlj
                    rkvec = 1.0_dp/kvec
                    strm = rkvec*(-rkvec*kexperfc + zd*(dexp1*derfc1-dexp2*derfc2) - &
                           rpieta*(dexp1*dexp3+dexp2*dexp4))
                    strm1 = costrm*strm
                    strm2 = costrm*kexperfc
                    if (lfinitestrain) then
                      do is = 1,nstrains
                        dqs(ieemf,is) = dqs(ieemf,is) + strm1*dg2ds(iv,is) - strm2*strainddetds(is)*straindet
                      enddo
                    else
                      do is = 1,nstrains
                        dqs(ieemf,is) = dqs(ieemf,is) + strm1*dg2ds(iv,is)
                      enddo
                      dqs(ieemf,1) = dqs(ieemf,1) - strm2
                      dqs(ieemf,2) = dqs(ieemf,2) - strm2
                    endif
                  endif
                endif
              endif
            enddo
          endif
        endif
!
!  End inner atom loop
!
      enddo
!
!  End outer atom loop
!
    enddo
!
!  Global matrix multiplies
!
    call pdgemm('t','n',numat,numat,numat,1.0d0,dqx,1,1,idesq,emat,1,1,idesd,0.0d0,dqtmp,1,1,idest)
    dqx(1:numat,1:natomsonnode) = dqtmp(1:numat,1:natomsonnode)
    call pdgemm('t','n',numat,numat,numat,1.0d0,dqy,1,1,idesq,emat,1,1,idesd,0.0d0,dqtmp,1,1,idest)
    dqy(1:numat,1:natomsonnode) = dqtmp(1:numat,1:natomsonnode)
    call pdgemm('t','n',numat,numat,numat,1.0d0,dqz,1,1,idesq,emat,1,1,idesd,0.0d0,dqtmp,1,1,idest)
    dqz(1:numat,1:natomsonnode) = dqtmp(1:numat,1:natomsonnode)
!
!  Multiply dqx/dqy/dqz by inverse matrix
!
    do iloc = 1,natomsonnode
      i = node2atom(iloc)
      if (lelementOK(nat(i)).and.nregionno(nsft+nrelf2a(i)).eq.1) then
        do j = 1,numat
          ind = 3*(j - 1)
          dqdxyz(ind+1,iloc) = dqdxyz(ind+1,iloc) - dqx(j,iloc)*angstoev
          dqdxyz(ind+2,iloc) = dqdxyz(ind+2,iloc) - dqy(j,iloc)*angstoev
          dqdxyz(ind+3,iloc) = dqdxyz(ind+3,iloc) - dqz(j,iloc)*angstoev
        enddo
      endif
    enddo
5   continue
!*************************
!  Real space summation  *
!*************************
    if (lnoreal.and.lstrain) then
!
!  Globalise dqs
!
      do kl = 1,nstrains
        do j = 1,numat
          dqtmp(j,1) = dqs(j,kl)
        enddo
        call sumall(dqtmp(1,1),dqtmp(1,2),numat,"dcharged","dqs")
        do j = 1,numat
          dqs(j,kl) = dqtmp(j,2)
        enddo
      enddo
!
!  Strain terms
!
      jeem = 0
      do jloc = 1,natomsonnode
        j = node2atom(jloc)
        if (lelementOK(nat(j)).and.nregionno(nsft+nrelf2a(j)).eq.1) then
          jeem = jeem + 1
          do kl = 1,nstrains
            sum = 0.0_dp
            do k = 1,neemfull
              sum = sum + dqs(k,kl)*emat(k,jeem)
            enddo
            dqds(kl,jloc) = dqds(kl,jloc) - sum
          enddo
        endif
      enddo
      goto 135
    endif
!
!  Start loop over atoms for coordinate differentiation
!
    ieem = 0
    do iloc = 1,natomsonnode
      i = node2atom(iloc)
      qlii = qf(i)
      xci = xclat(i)
      yci = yclat(i)
      zci = zclat(i)
      ind = 3*(i-1)
      ni = nat(i)
      if (lelementOK(nat(i)).and.nregionno(nsft+nrelf2a(i)).eq.1) then
        ieem = ieem + 1
        ieemf = neemfullrptr(i)
!
!  If QEq work out principal quantum number
!
        if (lqeq) then
          if (lmultiqrange) then
            nqr = nqrnow(ieemf)
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
          lhi = (ni.eq.1)
          if (lhi) then
!
!  Special case for hydrogen
!
            zetai = zetai + qlii*rconv
          endif
        elseif (lSandM) then   
          if (lmultiqrange) then
            nqr = nqrnow(ieemf)
          else
            nqr = 1
          endif
          zetai = zetarange(nqr,ni,neemtype)
        endif
      endif
!
!  Zero temporary derivative arrays
!
      do j = 1,numat
        dqx(j,iloc) = 0.0_dp
        dqy(j,iloc) = 0.0_dp
        dqz(j,iloc) = 0.0_dp
      enddo
      if (lSandM) then
        dchix(1:numat) = 0.0_dp
        dchiy(1:numat) = 0.0_dp
        dchiz(1:numat) = 0.0_dp
        if (lstrain) then
          dchis(1:6,iloc) = 0.0_dp
        endif
      endif
      jeem = 0
      do j = 1,numat
        qljj = qf(j)
        qlj = qljj*occuf(j)
        rx = xclat(j) - xci
        ry = yclat(j) - yci
        rz = zclat(j) - zci
        nj = nat(j)
        if (lelementOK(nat(j)).and.nregionno(nsft+nrelf2a(j)).eq.1) then
          jeem = jeem + 1
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
            lhj = (nj.eq.1)
            if (lhj) then
!
!  Special case for hydrogen
!
              zetaj = zetaj + qljj*rconv
            endif
          elseif (lSandM) then
            if (lmultiqrange) then
              nqrj = nqrnow(jeem)
            else
              nqrj = 1
            endif
            zetaj = zetarange(nqrj,nj,neemtype)
            znucj = znucrange(nqrj,nj,neemtype)
          endif
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
                if (lhi) then
                  dtrm1i = dtrm1i + qlii*d2zetari*rrr*rconv
                endif
                if (lhj) then
                  dtrm1j = dtrm1j + qljj*d2zetarj*rrr*rconv
                endif
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
              if (lelementOK(nat(i)).and.nregionno(nsft+nrelf2a(i)).eq.1) then
                dqx(ieemf,iloc) = dqx(ieemf,iloc) - qljj*dtrm1i*rxk
                dqy(ieemf,iloc) = dqy(ieemf,iloc) - qljj*dtrm1i*ryk
                dqz(ieemf,iloc) = dqz(ieemf,iloc) - qljj*dtrm1i*rzk
              endif
              if (lelementOK(nat(j)).and.nregionno(nsft+nrelf2a(j)).eq.1) then
                dqx(jeem,iloc) = dqx(jeem,iloc) - qlii*dtrm1j*rxk
                dqy(jeem,iloc) = dqy(jeem,iloc) - qlii*dtrm1j*ryk
                dqz(jeem,iloc) = dqz(jeem,iloc) - qlii*dtrm1j*rzk
              endif
              if (lSandM) then
                dtrm1zn = dtrm1zn*angstoev
                if (lelementOK(nat(i)).and.nregionno(nsft+nrelf2a(i)).eq.1) then
                  dchix(ieemf) = dchix(ieemf) - dtrm1zn*rxk
                  dchiy(ieemf) = dchiy(ieemf) - dtrm1zn*ryk
                  dchiz(ieemf) = dchiz(ieemf) - dtrm1zn*rzk
                endif
                if (lelementOK(nat(j)).and.nregionno(nsft+nrelf2a(j)).eq.1) then
                  dchix(jeem) = dchix(jeem) + dtrm1zn*rxk
                  dchiy(jeem) = dchiy(jeem) + dtrm1zn*ryk
                  dchiz(jeem) = dchiz(jeem) + dtrm1zn*rzk
                endif
              endif
              if (lstrain) then
                if (lelementOK(nat(i)).and.nregionno(nsft+nrelf2a(i)).eq.1) then
!
!  Strain
!
                  call real1strterm(ndim,rxk,ryk,rzk,0.0_dp,0.0_dp,0.0_dp,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
                  do is = 1,nstrains
                    ns = nstrptr(is)
                    dqs(ieemf,is) = dqs(ieemf,is) + qljj*dtrm1i*dr2ds(ns)
                  enddo
                  if (lSandM) then
                    do is = 1,nstrains
                      ns = nstrptr(is)
                      dchis(is,iloc) = dchis(is,iloc) + dtrm1zn*dr2ds(ns)
                    enddo
                  endif
                endif
              endif
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
      if (lSandM) then
!
!  For S & M, multiply dchix/dchiy/dchiz by inverse matrix
!
        jeem = 0
        do j = 1,numat
          if (lelementOK(nat(j)).and.nregionno(nsft+nrelf2a(j)).eq.1) then
            jeem = jeem + 1
            do k = 1,neemfull
              indj = 3*(neemfullptr(k) - 1)
              dsmtrm(indj+1,j) = dsmtrm(indj+1,j) - dchix(k)*emat(jeem,ieem)
              dsmtrm(indj+2,j) = dsmtrm(indj+2,j) - dchiy(k)*emat(jeem,ieem)
              dsmtrm(indj+3,j) = dsmtrm(indj+3,j) - dchiz(k)*emat(jeem,ieem)
            enddo
          endif
        enddo
      endif
!
      if (lfieldcfg(ncf).or.(ntdfieldcfg(ncf).gt.0)) then
!
!  For electric field, multiply dchix/dchiy/dchiz by inverse matrix
!
        jeem = 0
        do j = 1,numat
          if (lelementOK(nat(j)).and.nregionno(nsft+nrelf2a(j)).eq.1) then
            jeem = jeem + 1
            indj = 3*(j - 1)
            dqdxyz(indj+1,iloc) = dqdxyz(indj+1,iloc) - fieldx*emat(jeem,ieem)
            dqdxyz(indj+2,iloc) = dqdxyz(indj+2,iloc) - fieldy*emat(jeem,ieem)
            dqdxyz(indj+3,iloc) = dqdxyz(indj+3,iloc) - fieldz*emat(jeem,ieem)
          endif
        enddo
      endif
!
!  End of loop over outer atom
!
    enddo
!
!  Global matrix multiplies
!
    call pdgemm('t','n',numat,numat,numat,1.0d0,dqx,1,1,idesq,emat,1,1,idesd,0.0d0,dqtmp,1,1,idest)
    dqx(1:numat,1:natomsonnode) = dqtmp(1:numat,1:natomsonnode)
    call pdgemm('t','n',numat,numat,numat,1.0d0,dqy,1,1,idesq,emat,1,1,idesd,0.0d0,dqtmp,1,1,idest)
    dqy(1:numat,1:natomsonnode) = dqtmp(1:numat,1:natomsonnode)
    call pdgemm('t','n',numat,numat,numat,1.0d0,dqz,1,1,idesq,emat,1,1,idesd,0.0d0,dqtmp,1,1,idest)
    dqz(1:numat,1:natomsonnode) = dqtmp(1:numat,1:natomsonnode)
!
!  Multiply dqx/dqy/dqz by inverse matrix
!
    do iloc = 1,natomsonnode
      i = node2atom(iloc)
      if (lelementOK(nat(i)).and.nregionno(nsft+nrelf2a(i)).eq.1) then
        do j = 1,numat
          ind = 3*(j - 1)
          dqdxyz(ind+1,iloc) = dqdxyz(ind+1,iloc) - dqx(j,iloc)
          dqdxyz(ind+2,iloc) = dqdxyz(ind+2,iloc) - dqy(j,iloc)
          dqdxyz(ind+3,iloc) = dqdxyz(ind+3,iloc) - dqz(j,iloc)
        enddo
      endif
    enddo
!
    if (lSandM) then
!
!  Global sum of dsmtrm
!
      do i = 1,numat
        call sumall(dsmtrm(1,i),dtmp,3_i4*numat,"dcharged","dsmtrm")
        dsmtrm(1:3*numat,i) = dtmp(1:3*numat)
      enddo
!
!  For S & M, multiply dchix/dchiy/dchiz by inverse matrix
!
      ieem = 0
      do iloc = 1,natomsonnode
        i = node2atom(iloc)
        if (lelementOK(nat(i)).and.nregionno(nsft+nrelf2a(i)).eq.1) then
          ieem = ieem + 1
          do k = 1,neemfull
            ind = 3*(neemfullptr(k) - 1)
            dqdxyz(ind+1,iloc) = dqdxyz(ind+1,iloc) + dsmtrm(ind+1,i)
            dqdxyz(ind+2,iloc) = dqdxyz(ind+2,iloc) + dsmtrm(ind+2,i)
            dqdxyz(ind+3,iloc) = dqdxyz(ind+3,iloc) + dsmtrm(ind+3,i)
          enddo
        endif
      enddo
    endif
!
    if (lstrain) then
      if (lSandM) then
        call pdgemm('n','n',6,numat,numat,1.0d0,dchis,1,1,idess,emat,1,1,idesd,0.0d0,dstmp,1,1,idess) 
        dchis(1:6,1:natomsonnode) = dstmp(1:6,1:natomsonnode)
!
!  Strain terms
!
!  For S & M, multiply dchis by inverse matrix
!                 
        do iloc = 1,natomsonnode
          i = node2atom(iloc)
          if (lelementOK(nat(i)).and.nregionno(nsft+nrelf2a(i)).eq.1) then
            do kl = 1,nstrains
              dqds(kl,iloc) = dqds(kl,iloc) - dchis(kl,iloc)
            enddo
          endif
        enddo
      endif
!
!  Globalise dqs
!
      do kl = 1,nstrains
        do j = 1,numat
          dqtmp(j,1) = dqs(j,kl)
        enddo
        call sumall(dqtmp(1,1),dqtmp(1,2),numat,"dcharged","dqs")
        do j = 1,numat
          dqs(j,kl) = dqtmp(j,2) 
        enddo
      enddo
!
!  Strain terms
!
      jeem = 0
      do jloc = 1,natomsonnode
        j = node2atom(jloc)
        if (lelementOK(nat(j)).and.nregionno(nsft+nrelf2a(j)).eq.1) then
          jeem = jeem + 1
          do kl = 1,nstrains
            sum = 0.0_dp
            do k = 1,neemfull
              sum = sum + dqs(k,kl)*emat(k,jeem)
            enddo
            dqds(kl,jloc) = dqds(kl,jloc) - sum
          enddo
        endif
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
    do iloc = 1,natomsonnode
      i = node2atom(iloc)
      qlii = qf(i)
      xci = xclat(i)
      yci = yclat(i)
      zci = zclat(i)
      ind = 3*(i-1)
      ni = nat(i)
      if (lelementOK(nat(i)).and.nregionno(nsft+nrelf2a(i)).eq.1) then
        ieem = ieem + 1
        ieemf = neemfullrptr(i)
!
!  If QEq work out principal quantum number
!
        if (lqeq) then
          if (lmultiqrange) then
            nqr = nqrnow(ieemf)
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
          lhi = (ni.eq.1)
          if (lhi) then
!
!  Special case for hydrogen
!
            zetai = zetai + qlii*rconv
          endif
        elseif (lSandM) then
          if (lmultiqrange) then
            nqr = nqrnow(ieemf)
          else
            nqr = 1
          endif
          zetai = zetarange(nqr,ni,neemtype)
        endif
      endif
!
!  Zero temporary derivative arrays
!
      dqx(1:numat,iloc) = 0.0_dp
      dqy(1:numat,iloc) = 0.0_dp
      dqz(1:numat,iloc) = 0.0_dp
      if (lSandM) then
        dchix(1:numat) = 0.0_dp
        dchiy(1:numat) = 0.0_dp
        dchiz(1:numat) = 0.0_dp
        if (lstrain) then
          dchis(1,iloc) = 0.0_dp
        endif
      endif
      jeem = 0
!
!  Loop over other atoms and build daij/d(alpha)
!
      do j = 1,numat
        if (lelementOK(nat(j)).and.nregionno(nsft+nrelf2a(j)).eq.1) then
          jeem = jeem + 1
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
              lhj = (nj.eq.1)
              if (lhj) then
!
!  Special case for hydrogen
!
                zetaj = zetaj + qljj*rconv
              endif
            elseif (lSandM) then
              if (lmultiqrange) then
                nqrj = nqrnow(jeem)
              else
                nqrj = 1
              endif
              zetaj = zetarange(nqrj,nj,neemtype)
              znucj = znucrange(nqrj,nj,neemtype)
            endif
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
            if (lelementOK(nat(i)).and.nregionno(nsft+nrelf2a(i)).eq.1) then
              dqx(ieemf,iloc) = dqx(ieemf,iloc) - qljj*dqme(1)
              dqy(ieemf,iloc) = dqy(ieemf,iloc) - qljj*dqme(2)
              dqz(ieemf,iloc) = dqz(ieemf,iloc) - qljj*dqme(3)
            endif
            if (lelementOK(nat(j)).and.nregionno(nsft+nrelf2a(j)).eq.1) then
              dqx(jeem,iloc) = dqx(jeem,iloc) - qlii*dqme(1)
              dqy(jeem,iloc) = dqy(jeem,iloc) - qlii*dqme(2)
              dqz(jeem,iloc) = dqz(jeem,iloc) - qlii*dqme(3)
            endif
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
                if (lhi) then
                  dtrm1i = dtrm1i + qlii*d2zetari*rrr*rconv
                endif
                if (lhj) then
                  dtrm1j = dtrm1j + qljj*d2zetarj*rrr*rconv
                endif
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
              if (lelementOK(nat(i)).and.nregionno(nsft+nrelf2a(i)).eq.1) then
                dqx(ieemf,iloc) = dqx(ieemf,iloc) - qljj*dtrm1i*rx
                dqy(ieemf,iloc) = dqy(ieemf,iloc) - qljj*dtrm1i*ry
                dqz(ieemf,iloc) = dqz(ieemf,iloc) - qljj*dtrm1i*rz
              endif
              if (lelementOK(nat(j)).and.nregionno(nsft+nrelf2a(j)).eq.1) then
                dqx(jeem,iloc) = dqx(jeem,iloc) - qlii*dtrm1j*rx
                dqy(jeem,iloc) = dqy(jeem,iloc) - qlii*dtrm1j*ry
                dqz(jeem,iloc) = dqz(jeem,iloc) - qlii*dtrm1j*rz
              endif
              if (lSandM) then
                dtrm1zn = dtrm1zn*angstoev
                if (lelementOK(nat(i)).and.nregionno(nsft+nrelf2a(i)).eq.1) then
                  dchix(ieemf) = dchix(ieemf) - dtrm1zn*rx
                  dchiy(ieemf) = dchiy(ieemf) - dtrm1zn*ry
                  dchiz(ieemf) = dchiz(ieemf) - dtrm1zn*rz
                endif
                if (lelementOK(nat(j)).and.nregionno(nsft+nrelf2a(j)).eq.1) then
                  dchix(jeem) = dchix(jeem) + dtrm1zn*rx
                  dchiy(jeem) = dchiy(jeem) + dtrm1zn*ry
                  dchiz(jeem) = dchiz(jeem) + dtrm1zn*rz
                endif
              endif
            endif
            if (lstrain) then
!             
!  Strain 
!             
              if (lelementOK(nat(i)).and.nregionno(nsft+nrelf2a(i)).eq.1) then
                call real1strterm(ndim,rx,ry,rz,0.0_dp,0.0_dp,0.0_dp,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
                dqs(ieemf,1) = dqs(ieemf,1) + qljj*dtrm1i*dr2ds(1)
                if (lSandM) then
                  dchis(1,iloc) = dchis(1,iloc) + dtrm1zn*dr2ds(1)
                endif
              endif
            endif
          enddo
!       
!  End of loop over lattice vectors
!
        endif
      enddo
!
      if (lSandM) then
!
!  For S & M, multiply dchix/dchiy/dchiz by inverse matrix
!
        jeem = 0
        do j = 1,numat
          if (lelementOK(nat(j)).and.nregionno(nsft+nrelf2a(j)).eq.1) then
            jeem = jeem + 1
            do k = 1,neemfull
              indj = 3*(neemfullptr(k) - 1)
              dsmtrm(indj+1,j) = dsmtrm(indj+1,j) - dchix(k)*emat(jeem,ieem)
              dsmtrm(indj+2,j) = dsmtrm(indj+2,j) - dchiy(k)*emat(jeem,ieem)
              dsmtrm(indj+3,j) = dsmtrm(indj+3,j) - dchiz(k)*emat(jeem,ieem)
            enddo
          endif
        enddo
      endif
!
      if (lfieldcfg(ncf).or.(ntdfieldcfg(ncf).gt.0)) then
!
!  For electric field, multiply dchix/dchiy/dchiz by inverse matrix
!
        jeem = 0
        do j = 1,numat
          if (lelementOK(nat(j)).and.nregionno(nsft+nrelf2a(j)).eq.1) then
            jeem = jeem + 1
            indj = 3*(j - 1)
            dqdxyz(indj+1,iloc) = dqdxyz(indj+1,iloc) - fieldx*emat(jeem,ieem)
            dqdxyz(indj+2,iloc) = dqdxyz(indj+2,iloc) - fieldy*emat(jeem,ieem)
            dqdxyz(indj+3,iloc) = dqdxyz(indj+3,iloc) - fieldz*emat(jeem,ieem)
          endif
        enddo
      endif
!
!  End of outer loop over atoms
!
    enddo
!
!  Global matrix multiplies
!
    call pdgemm('t','n',numat,numat,numat,1.0d0,dqx,1,1,idesq,emat,1,1,idesd,0.0d0,dqtmp,1,1,idest)
    dqx(1:numat,1:natomsonnode) = dqtmp(1:numat,1:natomsonnode)
    call pdgemm('t','n',numat,numat,numat,1.0d0,dqy,1,1,idesq,emat,1,1,idesd,0.0d0,dqtmp,1,1,idest)
    dqy(1:numat,1:natomsonnode) = dqtmp(1:numat,1:natomsonnode)
    call pdgemm('t','n',numat,numat,numat,1.0d0,dqz,1,1,idesq,emat,1,1,idesd,0.0d0,dqtmp,1,1,idest)
    dqz(1:numat,1:natomsonnode) = dqtmp(1:numat,1:natomsonnode)
!
!  Multiply dqx/dqy/dqz by inverse matrix
!
    do iloc = 1,natomsonnode
      i = node2atom(iloc)
      if (lelementOK(nat(i)).and.nregionno(nsft+nrelf2a(i)).eq.1) then
        do j = 1,numat
          ind  = 3*(j - 1)
          dqdxyz(ind+1,iloc) = dqdxyz(ind+1,iloc) - dqx(j,iloc)
          dqdxyz(ind+2,iloc) = dqdxyz(ind+2,iloc) - dqy(j,iloc)
          dqdxyz(ind+3,iloc) = dqdxyz(ind+3,iloc) - dqz(j,iloc)
        enddo
      endif
    enddo
!
    if (lSandM) then
!
!  Global sum of dsmtrm
!
      do i = 1,numat
        call sumall(dsmtrm(1,i),dtmp,3_i4*numat,"dcharged","dsmtrm")
        dsmtrm(1:3*numat,i) = dtmp(1:3*numat)
      enddo
!
!  For S & M, multiply dchix/dchiy/dchiz by inverse matrix
!
      ieem = 0
      do iloc = 1,natomsonnode
        i = node2atom(iloc)
        if (lelementOK(nat(i)).and.nregionno(nsft+nrelf2a(i)).eq.1) then
          ieem = ieem + 1
          do k = 1,neemfull
            ind = 3*(neemfullptr(k) - 1)
            dqdxyz(ind+1,iloc) = dqdxyz(ind+1,iloc) + dsmtrm(ind+1,i)
            dqdxyz(ind+2,iloc) = dqdxyz(ind+2,iloc) + dsmtrm(ind+2,i)
            dqdxyz(ind+3,iloc) = dqdxyz(ind+3,iloc) + dsmtrm(ind+3,i)
          enddo
        endif
      enddo
    endif
!
    if (lstrain) then
      if (lSandM) then
        call pdgemm('n','n',1,numat,numat,1.0d0,dchis,1,1,idess,emat,1,1,idesd,0.0d0,dstmp,1,1,idess)
        dchis(1,1:natomsonnode) = dstmp(1,1:natomsonnode)
!
!  Strain terms
!
!  For S & M, multiply dchis by inverse matrix
!
        do iloc = 1,natomsonnode
          i = node2atom(iloc)
          if (lelementOK(nat(i)).and.nregionno(nsft+nrelf2a(i)).eq.1) then
            dqds(1,iloc) = dqds(1,iloc) - dchis(1,iloc)
          endif
        enddo
      endif
!
!  Globalise dqs
!
      do j = 1,numat
        dqtmp(j,1) = dqs(j,1)
      enddo
      call sumall(dqtmp(1,1),dqtmp(1,2),numat,"dcharged","dqs")
      do j = 1,numat
        dqs(j,1) = dqtmp(j,2) 
      enddo
!
!  Strain terms
!
      jeem = 0
      do jloc = 1,natomsonnode
        j = node2atom(jloc)
        if (lelementOK(nat(j)).and.nregionno(nsft+nrelf2a(j)).eq.1) then
          jeem = jeem + 1
          sum = 0.0_dp
          do k = 1,neemfull
            sum = sum + dqs(k,1)*emat(k,jeem)
          enddo
          dqds(1,jloc) = dqds(1,jloc) - sum
        endif
      enddo
    endif
!******************************
!  End of cluster / 1-D case  *
!******************************
  endif
135 continue
!***********************************************************************************
!  Enforce sum rules that total charge derivative equals zero for each coordinate  *
!***********************************************************************************
!
!  Sum elements and count number of active atoms
!
  do i = 1,numat
    dqx(i,1) = 0.0_dp
    dqy(i,1) = 0.0_dp
    dqz(i,1) = 0.0_dp
  enddo
  do iloc = 1,natomsonnode
    jx = - 2
    jy = - 1
    jz =   0
    do j = 1,numat
      jx = jx + 3
      jy = jy + 3
      jz = jz + 3
      dqx(j,1) = dqx(j,1) + dqdxyz(jx,iloc)
      dqy(j,1) = dqy(j,1) + dqdxyz(jy,iloc)
      dqz(j,1) = dqz(j,1) + dqdxyz(jz,iloc)
    enddo
  enddo
!
!  Globalise dqx / dqy / dqz 
!
  call sumall(dqx,dqtmp(1,1),numat,"dcharged","dqx")
  dqx(1:numat,1) = dqtmp(1:numat,1)
  call sumall(dqy,dqtmp(1,1),numat,"dcharged","dqx")
  dqy(1:numat,1) = dqtmp(1:numat,1)
  call sumall(dqz,dqtmp(1,1),numat,"dcharged","dqx")
  dqz(1:numat,1) = dqtmp(1:numat,1)
!
!  Average error and subtract from active elements
!
  if (neemfull.gt.0) then
    do i = 1,numat
      dqx(i,1) = - dqx(i,1)/dble(neemfull)
      dqy(i,1) = - dqy(i,1)/dble(neemfull)
      dqz(i,1) = - dqz(i,1)/dble(neemfull)
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
          dqdxyz(jx,iloc) = dqdxyz(jx,iloc) + dqx(j,1)
          dqdxyz(jy,iloc) = dqdxyz(jy,iloc) + dqy(j,1)
          dqdxyz(jz,iloc) = dqdxyz(jz,iloc) + dqz(j,1)
        enddo
      endif
    enddo
  endif
!
!  Free local memory
!
  if (lstrain) then
    deallocate(dqs,stat=status)
    if (status/=0) call deallocate_error('dcharged','dqs')
  endif
  if (lSandM) then
    if (lstrain) then
      deallocate(dstmp,stat=status)
      if (status/=0) call deallocate_error('dcharged','dstmp')
      deallocate(dchis,stat=status)
      if (status/=0) call deallocate_error('dcharged','dchis')
    endif
    deallocate(dchiz,stat=status)
    if (status/=0) call deallocate_error('dcharged','dchiz')
    deallocate(dchiy,stat=status)
    if (status/=0) call deallocate_error('dcharged','dchiy')
    deallocate(dchix,stat=status)
    if (status/=0) call deallocate_error('dcharged','dchix')
  endif
  deallocate(dqz,stat=status)
  if (status/=0) call deallocate_error('dcharged','dqz')
  deallocate(dqy,stat=status)
  if (status/=0) call deallocate_error('dcharged','dqy')
  deallocate(dqx,stat=status)
  if (status/=0) call deallocate_error('dcharged','dqx')
  deallocate(dqtmp,stat=status)
  if (status/=0) call deallocate_error('dcharged','dqtmp')
  deallocate(neemfullrptr,stat=status)
  if (status/=0) call deallocate_error('dcharged','neemfullrptr')
  deallocate(neemfullptr,stat=status)
  if (status/=0) call deallocate_error('dcharged','neemfullptr')
!
  if (lSandM) then
    deallocate(dtmp,stat=status)
    if (status/=0) call deallocate_error('dcharged','dtmp')
    deallocate(dsmtrm,stat=status)
    if (status/=0) call deallocate_error('dcharged','dsmtrm')
  endif
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
        allocate(dtmp2(ntmp),stat=status)
        if (status/=0) call outofmemory('dcharged','dtmp2')
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
            call MPI_IRecv(dtmp2,ntmp,MPI_double_precision,nnode, &
                           ntag,MPI_Comm_World,Request,MPIerror)
          endif
!
!  Pass data to ioproc for writing
!
          if (iloc.gt.0) then
            dtmp2(1:nstrains) = dqds(1:nstrains,iloc)
!
!  Post send
!
            call MPI_ISend(dtmp2,ntmp,MPI_double_precision,0, &
                           ntag,MPI_Comm_World,Request,MPIerror)
          endif
          if (ioproc.or.iloc.gt.0) then
            call MPI_WaitAll(1,Request,StatMPI,MPIerror)
          endif
          if (ioproc) then
!
!  Write on I/O node
!
            write(ioout,'(i6,4x,6f10.6)') i,(dtmp2(j),j=1,nstrains)
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
        deallocate(dtmp2,stat=status)
        if (status/=0) call deallocate_error('dcharged','dtmp2')
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
      allocate(dtmp2(ntmp),stat=status)
      if (status/=0) call outofmemory('dcharged','dtmp2')
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
          call MPI_IRecv(dtmp2,ntmp,MPI_double_precision,nnode, &
                         ntag,MPI_Comm_World,Request,MPIerror)
        endif
!
!  Pass data to ioproc for writing
!
        if (iloc.gt.0) then
          dtmp2(1:3*numat) = dqdxyz(1:3*numat,iloc)
!
!  Post send
!
          call MPI_ISend(dtmp2,ntmp,MPI_double_precision,0, &
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
            write(ioout,'(i6,'' x '',4x,f10.6)') j,dtmp2(indj+1)
            write(ioout,'(i6,'' y '',4x,f10.6)') j,dtmp2(indj+2)
            write(ioout,'(i6,'' z '',4x,f10.6)') j,dtmp2(indj+3)
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
      deallocate(dtmp2,stat=status)
      if (status/=0) call deallocate_error('dcharged','dtmp2')
    endif
  endif
#ifdef TRACE
  call trace_out('dcharged')
#endif
#else
  call outerror('dcharged called when not compiled with MPI',0_i4)
  call stopnow('dcharged')
#endif
!
  return
  end
