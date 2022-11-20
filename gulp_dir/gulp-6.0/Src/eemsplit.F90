  subroutine eemsplit(lmain,lgrad1,lgrad2)
!
!  Subroutine for performing electronegativity equilisation calcns
!  using the split charge approach.
!
!  NB: Does not currently allow iterative QEq for hydrogen or
!      partial occupancy
!
!  lmain = .true.  => call from gulp, one off calc with output
!        = .false. => call from energy, needed for charges
!
!  If second derivatives are being used then tmat must be stored
!  on channel 54 while charges are calculated.
!
!   4/18 Created from eem
!   5/18 Call to dcharge replaced with equivalent for split bonds
!   5/18 Call to eemsplitd added for parallel runs
!   5/18 oldeem and lelementOK handling moved to seteem
!   5/18 Multiple qranges added
!   6/18 Check that number of bond charges does not equal atoms 
!        in molecule added
!   6/18 Extra flag passed to dbcgsolve
!   6/18 e0range added
!   6/18 Setting of bonds now calls seteembond
!   6/18 EEM matrix no longer overwrites derv2/dervi to avoid issues
!        when using numerical derivatives
!   7/18 lDoChargeDerv now replaced by ldcharge
!   7/18 emat and ematb passed to dchargesplit routine
!   7/18 Changed so that no symmetry is used in calculation of charges
!   7/18 Trap for zero molecules added
!   7/18 Coordinate handling for lsymopt case fixed
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
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
!  Julian Gale, CIC, Curtin University, September 2019
!
  use control
  use configurations
  use current
  use derivatives
  use element
  use eembonds
  use eemdata
  use energies
  use field,           only : lfieldcfg, ntdfieldcfg
  use iochannels
  use mdlogic,         only : lmd
  use moldyn,          only : labscoany, labsco
  use molecule,        only : nmol
  use parallel
  use partial
  use symmetry
  use times
#ifdef TRACE
  use trace,           only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  logical, intent(in)                          :: lgrad1
  logical, intent(in)                          :: lgrad2
  logical, intent(in)                          :: lmain
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ib1
  integer(i4)                                  :: ib2
  integer(i4)                                  :: ifail
  integer(i4)                                  :: ii
  integer(i4)                                  :: iter
  integer(i4)                                  :: j
  integer(i4)                                  :: jb1
  integer(i4)                                  :: jb2
  integer(i4)                                  :: jj
  integer(i4)                                  :: n
  integer(i4),                            save :: ncfold = 0
  integer(i4)                                  :: neemfoc
  integer(i4)                                  :: ni
  integer(i4)                                  :: niter
  integer(i4)                                  :: nitereem
  integer(i4)                                  :: nj
  integer(i4)                                  :: nloc
  integer(i4)                                  :: nqr
  integer(i4)                                  :: nqrib1
  integer(i4)                                  :: nqrib2
  integer(i4), dimension(:), allocatable       :: nlocptr
  integer(i4), dimension(:), allocatable       :: nqrlast
  integer(i4)                                  :: nmax
  integer(i4)                                  :: nr
  integer(i4)                                  :: nregioni
  integer(i4)                                  :: status
  logical,     dimension(:), allocatable       :: leemfoc
  logical                                      :: lconverged
  logical                                      :: lfound
  logical                                      :: literate
  logical                                      :: lqchange
  logical                                      :: lqiter
  real(dp)                                     :: chii
  real(dp),    dimension(:,:), allocatable     :: emat
  real(dp),    dimension(:,:), allocatable     :: ematb
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: err
  real(dp)                                     :: eself_before
  real(dp)                                     :: q0i
  real(dp)                                     :: qguesstot
  real(dp)                                     :: qd
  real(dp)                                     :: qdiff
  real(dp)                                     :: qi
  real(dp),    dimension(:), allocatable       :: qnmr
  real(dp)                                     :: qsum
  real(dp)                                     :: qtot
  real(dp)                                     :: rmui
  real(dp)                                     :: rnguess
  real(dp)                                     :: time1
  real(dp)                                     :: time2
  real(dp),    dimension(:), allocatable       :: oldqf
  real(dp),    dimension(:), allocatable       :: vfield
  real(dp),    dimension(:), allocatable       :: z
  real(dp),    dimension(:), allocatable       :: z2
#ifdef TRACE
  call trace_in('eemsplit')
#endif
!
!  If this is a parallel run then call distributed memory version
!
  if (nprocs.gt.1) then
    call eemsplitd(lmain,lgrad1,lgrad2)
#ifdef TRACE
    call trace_out('eemsplit')
#endif
    return
  endif
!
  time1 = g_cpu_time()
!
!  Allocate local memory
!
  allocate(leemfoc(numat),stat=status)
  if (status/=0) call outofmemory('eemsplit','leemfoc')
  if (lmultiqrange) then
    allocate(nqrlast(numat),stat=status)
    if (status/=0) call outofmemory('eemsplit','nqrlast')
  endif
!
!  Set flag for iterative charges - can't be used for all algorithms
!
  lqiter = literativeQ
  if (lgrad2.or.(lgrad1.and.(lfieldcfg(ncf).or.(ntdfieldcfg(ncf).gt.0))).or.ldcharge) lqiter = .false.
  if (lallbonds) lqiter = .true.
!
  qsum = 0.0_dp
  qtot = 0.0_dp
  neem = 0
  neemrptr(1:numat) = 0
!
!  Check elements
!
  do i = 1,numat
    ni = nat(i)
    qi = qf(i)
    if (lelementOK(ni).and.nregionno(nsft+nrelf2a(i)).eq.1) then
      neem = neem + 1
      neemptr(neem) = i
      neemrptr(i) = neem
      qsum = qsum - qf(i)*occuf(i)
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
          nqrlast(neem) = nqr
          nqrnow(neem) = nqr
        else
          nqrlast(neem) = 1
          nqrnow(neem) = 1
        endif
      endif
    elseif (ni.gt.maxele) then
      call outerror('cannot use EEM with shells present',0_i4)
      call stopnow('eemsplit')
    else
      qtot = qtot + qf(i)*occuf(i)
    endif
  enddo
!
!  Check that there are molecules otherwise the charges will be zero
!
  if (nmol.eq.0) goto 100
!
!  Set up split bond variables
!
  call seteembond(lallbonds)
!
!  Now find the number of fully occupied sites for EEM/QEq
!
  leemfoc(1:numat) = .false.
  do i = 1,neem
    ii = iocptr(neemptr(i))
    leemfoc(ii) = .true.
  enddo
  neemfoc = 0
  do i = 1,ncfoc
    if (leemfoc(i)) neemfoc = neemfoc + 1
  enddo
!
!  For now split bonds is not enabled with partial occupancy
!
  if (neemfoc.ne.neem) then
    call outerror('partial occupancy not currently enabled in eemsplit',0_i4)
    call stopnow('eemsplit')
  endif
!
!  Check the memory for the linear arrays
!
  nmax = max(neembond,numat)
  if (nmax.gt.maxat) then
    maxat = nmax
    call changemaxat
  endif
!
!  Allocate memory for matrices
!
  allocate(emat(numat,numat),stat=status)
  if (status/=0) call outofmemory('eemsplit','emat')
  allocate(ematb(neembond,neembond),stat=status)
  if (status/=0) call outofmemory('eemsplit','ematb')
!
!  Decide on total EEM fragment charge - for cluster there is no constraint so use sum of initial charges.
!  For periodic system, charge must be equal to the negative sum of the non-EEM ion charges.
!
  if (ndim.eq.0) then
    qtot = qsum
  endif
!
!  Set iteration logical if Q is multirange
!
  literate = lmultiqrange
  if (literate) then
    nitereem = nqeqitermax
  else
    nitereem = 1
  endif
!
!  Allocate local memory that depends on neembond
!
  allocate(z(max(neembond,numat)),stat=status)
  if (status/=0) call outofmemory('eemsplit','z')
  allocate(z2(numat),stat=status)
  if (status/=0) call outofmemory('eemsplit','z2')
!
!  If iterative then set up pointers to local elements
!
  if (lqiter) then
    nloc = neembond
    allocate(nlocptr(neembond),stat=status)
    if (status/=0) call outofmemory('eemsplit','nlocptr')
    do i = 1,neembond
      nlocptr(i) = i
    enddo
  endif
  if (literate) then
    allocate(oldqf(numat),stat=status)
    if (status/=0) call outofmemory('eemsplit','oldqf')
  endif
  if (lfieldcfg(ncf).or.(ntdfieldcfg(ncf).gt.0)) then
    allocate(vfield(numat),stat=status)
    if (status/=0) call outofmemory('eemsplit','vfield')
  endif
!
!  Calculate reciprocal lattice vectors
!
  if (lewald.and.ndim.gt.1) call kindex
!
!  For 1-D case set guess at charges based EEM parameters and scaled up by 1.5 to allow for increase in ionicity.
!
  if (ndim.eq.1.and.ncf.ne.ncfold) then
    ncfold = ncf
    qguesstot = 0.0_dp
    rnguess = 0.0_dp
    do i = 1,neem
      ii = neemptr(i)
      ni = nat(ii)
      if (lmultiqrange) then
        nqr = nqrnow(i)
      else
        nqr = 1
      endif
      chii = chirange(nqr,ni,neemtype)
      rmui = murange(nqr,ni,neemtype)
      q0i  = q0range(nqr,ni,neemtype)
      qf(ii) = q0i - chii/rmui
      qguesstot = qguesstot + qf(ii)*occuf(ii)
      rnguess = rnguess + occuf(ii)
    enddo
    qguesstot = (qguesstot + qtot)/rnguess
    do i = 1,neem
      ii = neemptr(i)
      qf(ii) = qf(ii) - qguesstot
      if (abs(qtot).lt.1.0d-12) qf(ii) = 1.5_dp*qf(ii)
    enddo
    do i = 1,nasym
      qa(i) = qf(nrela2f(i))
    enddo
    if (index(keyword,'debu').ne.0.and.ioproc) then
      write(ioout,'(/,''  Initial guess for 1-D variable charges :'',/)')
      write(ioout,'('' Atom        Q'')')
      do i = 1,neem
        ii = neemptr(i)
        write(ioout,'(i5,1x,f12.6)') ii,qf(ii)
      enddo
      write(ioout,'(/)')
    endif
  endif
!
!  Setup coordinates
!
  if (lsymopt) then
    do i = 1,nasym
      nr = nrela2f(i)
      xalat(i) = xclat(nr)
      yalat(i) = yclat(nr)
      zalat(i) = zclat(nr)
    enddo
  else
!
!  Avoid overwriting absolute coordinates for MD
!
    if (lmd.and.labscoany) then
      do i = 1,numat
        if (.not.labsco(i)) then
          xalat(i) = xclat(i)
          yalat(i) = yclat(i)
          zalat(i) = zclat(i)
        endif
      enddo
    else
      do i = 1,numat
        xalat(i) = xclat(i)
        yalat(i) = yclat(i)
        zalat(i) = zclat(i)
      enddo
    endif
  endif
!
!  Store charges for convergence check
!
  if (literate) then
    do i = 1,numat
      oldqf(i) = qf(i)
    enddo
  endif
!
!  Generate electric field potential
!
  if (lfieldcfg(ncf).or.(ntdfieldcfg(ncf).gt.0)) then
    vfield(1:numat) = 0.0_dp
    call electricfieldpotl(vfield)
  endif
!****************************
!  Start of iterative loop  *
!****************************
  lconverged = .false.
  niter = 0
  if (literate.and.lmain.and.ioproc) then
    write(ioout,'(''  Iterative solution of charge equilibration :'',/)')
  endif
  do while (niter.lt.nitereem.and..not.lconverged)
    niter = niter + 1
!
!  Zero right hand vector
!
    z = 0.0_dp
!*******************************************
!  Generate electrostatic potential terms  *
!*******************************************
    if (lnoqeem) then
      emat(1:numat,1:numat) = 0.0_dp
    else
      call genpot(emat,numat,z,1_i4)
    endif
!
!  Reduce emat to neem x neem form
!
    do i = 1,neem
!
!  Zero storage vector for emat array
!
      z2(1:numat) = 0.0_dp
!
!  Place i-j potential terms into emat
!
      do j = 1,numat
        jj = neemrptr(j)
!
!  Variable j charge case
!
        if (jj.eq.j) then
          z2(jj) = z2(jj) + emat(j,neemptr(i))*occuf(j)
        else
          z(i) = z(i) - qf(j)*emat(j,neemptr(i))*occuf(j)
        endif
      enddo
!
!  Copy temporary storage vector back into emat array
!
      do j = 1,neem
        emat(j,i) = z2(j)
      enddo
    enddo
!
!  Now create matrix for split bonds
!
    ematb(1:neembond,1:neembond) = 0.0_dp
    do i = 1,neembond
      ib1 = neembonded(1,i)
      ib2 = neembonded(2,i)
      do j = 1,neembond
        jb1 = neembonded(1,j)
        jb2 = neembonded(2,j)
        ematb(j,i) = ematb(j,i) + ((emat(jb1,ib1) - emat(jb2,ib1)) - (emat(jb1,ib2) - emat(jb2,ib2)))
      enddo
    enddo
    if (lSandM) then
!
!  From S & M, where z has been set without reference to neem, reduce elements to those that are needed
!
      do i = 1,neem
        z2(i) = z(neemptr(i))
      enddo
    else
!
!  For methods other than S & M, copy into z2 for z reduction
!
      do i = 1,neem
        z2(i) = z(i)
      enddo
    endif
!
!  Reduce z to split bond charge form
!
    z(1:neembond) = 0.0_dp
    do i = 1,neembond
      ib1 = neembonded(1,i)
      ib2 = neembonded(2,i)
      z(i) = z2(ib1) - z2(ib2)
    enddo
!********************************
!  Form matrix of coefficients  *
!********************************
    do i = 1,neembond
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
          ematb(j,i) = ematb(j,i) + 2.0_dp*murange(nqrib1,nat(ib1),neemtype)
        endif
        if (ib1.eq.jb2) then
          ematb(j,i) = ematb(j,i) - 2.0_dp*murange(nqrib1,nat(ib1),neemtype)
        endif
        if (ib2.eq.jb1) then
          ematb(j,i) = ematb(j,i) - 2.0_dp*murange(nqrib2,nat(ib2),neemtype)
        endif
        if (ib2.eq.jb2) then
          ematb(j,i) = ematb(j,i) + 2.0_dp*murange(nqrib2,nat(ib2),neemtype)
        endif
      enddo
    enddo
!
!  Add external potential
!
    do i = 1,neembond
      ib1 = neembonded(1,i)
      ib2 = neembonded(2,i)
      z(i) = z(i) - extpotcfg(nsft+nrelf2a(ib1))
      z(i) = z(i) + extpotcfg(nsft+nrelf2a(ib2))
    enddo
    if (lmultiqrange) then
      do i = 1,neembond
        ib1 = neembonded(1,i)
        ib2 = neembonded(2,i)
        ni = nat(ib1)
        nj = nat(ib2)
        nqrib1 = nqrnow(neemrptr(ib1))
        nqrib2 = nqrnow(neemrptr(ib2))
        z(i) = z(i) + (chirange(nqrib2,nj,neemtype) - chirange(nqrib1,ni,neemtype))
      enddo
    else
      do i = 1,neembond
        ib1 = neembonded(1,i)
        ib2 = neembonded(2,i)
        ni = nat(ib1)
        nj = nat(ib2)
        z(i) = z(i) + (chirange(1,nj,neemtype) - chirange(1,ni,neemtype))
      enddo
    endif
!
    if (lfieldcfg(ncf).or.(ntdfieldcfg(ncf).gt.0)) then
      do i = 1,neembond
        ib1 = neembonded(1,i)
        ib2 = neembonded(2,i)
        z(i) = z(i) - vfield(ib1) + vfield(ib2)
      enddo
    endif
!
!  Debugging output
!
    if (index(keyword,'debu').ne.0.and.ioproc) then
      write(ioout,'(/,''  EEM/QEq Matrix :'',/)')
      do i = 1,neembond
        write(ioout,'(10(1x,f9.5))')(ematb(j,i),j=1,neembond),z(i)
      enddo
    endif
    if (lqiter) then
!***************************
!  Iterative charge solve  *
!***************************
!
!  Initialise bond charge array
!
      qbond(1:neembond) = 0.0_dp
!
!  Solve using iterative route
!
      call dbcgsolve(0_i4,neembond,neembond,nlocptr,ematb,neembond,z,qbond,qitertol,nqitermax,iter,err)
!
      if (ioproc) then
        if (index(keyword,'verb').ne.0) then
          write(ioout,'('' Number of iterations / error in dbcgsolve = '',i4,1x,f16.14)') iter,err
        endif
      endif
!
!  Was iterative solution successful?
!
      if (iter.ge.nqitermax) then
        call outerror('iterative charge solution failed in eemsplit',0_i4)
        call stopnow('eemsplit')
      endif
    else
!******************
!  Invert matrix  *
!******************
      ifail = 0
      n = neembond
!************************
!  Symmetric inversion  *
!************************
      call matrix_inversion_library(n,1_i4,neembond,nblocksize,ematb,0_i4,ifail)
!
!  Was inversion successful?
!
      if (ifail.ne.0) then
        call outerror('matrix inversion failed in EEM/QEq',0_i4)
        call stopnow('eemsplit')
      endif
!  
!  Multiply inverse matrix and chi matrix to get bond charges
!
      do i = 1,neembond
        qbond(i) = 0.0_dp
        do j = 1,neembond
          qbond(i) = qbond(i) + z(j)*ematb(j,i)
        enddo
      enddo
    endif
!
!  Set atom charges based on split bond charges
!
    qf(1:numat) = 0.0_dp
    do i = 1,neembond
      ib1 = neembonded(1,i)
      ib2 = neembonded(2,i)
      qf(ib1) = qf(ib1) + qbond(i)
      qf(ib2) = qf(ib2) - qbond(i)
    enddo
!
!  Transfer charges to qa
!
    do i = 1,nasym
      nr = nrela2f(i)
      qa(i) = qf(nr)
    enddo
!
    if (literate) then
!
!  For multiple q ranges check whether range has changed
!
      if (lmultiqrange) then
        nqrlast(1:neem) = nqrnow(1:neem)
        lqchange = .false.
        do i = 1,neem
          ii = neemptr(i)
          ni = nat(ii)
          qi = qf(ii)
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
            nqrnow(i) = nqr
          else
            nqrnow(i) = 1
          endif
          lqchange = (nqrnow(i).ne.nqrlast(i))
        enddo
      else
        lqchange = .false.
      endif
!
!  Check for convergence based on charge differences
!
      qdiff = 0.0_dp
      do i = 1,numat
        qd = qf(i) - oldqf(i)
        qdiff = qdiff + abs(qd)
      enddo
      qdiff = qdiff/dble(numat)
      lconverged = (qdiff.lt.qeqscfcrit.and..not.lqchange)
      if (lmain.and.ioproc) then
        write(ioout,'(''  ** Cycle : '',i4,'' Qdiff : '',f10.8)') niter,qdiff
      endif
      if (.not.lconverged) then
        do i = 1,neem
          ii = neemptr(i)
          oldqf(ii) = qf(ii)
        enddo
      endif
    endif
!
!  Transfer charges to qa
!
    do i = 1,nasym
      nr = nrela2f(i)
      qa(i) = qf(nr)
    enddo
!*****************************
!  End loop over iterations  *
!*****************************
  enddo
!
!  Store charges in configurational array
!
  do i = 1,nasym
    qlcfg(nsft+i) = qa(i)
  enddo
!**************************
!  Calculate self energy  *
!**************************
  eself = 0.0_dp
  do i = 1,neembond
    ib1 = neembonded(1,i)
    ib2 = neembonded(2,i)
    ni = nat(ib1)
    nj = nat(ib2)
!
    if (lmultiqrange) then
      nqrib1 = nqrnow(neemrptr(ib1))
      nqrib2 = nqrnow(neemrptr(ib2))
    else
      nqrib1 = 1
      nqrib2 = 1
    endif
!
    eself_before = eself
    qi = qbond(i)
    eself = eself + qi*(chirange(nqrib1,ni,neemtype)-chirange(nqrib2,nj,neemtype))
!
!  Add external potential for site
!
    eself = eself + qi*extpotcfg(nsft+nrelf2a(ib1))
    eself = eself - qi*extpotcfg(nsft+nrelf2a(ib2))
!
    nregioni = nregionno(nsft+nrelf2a(ib1))
    eregion2region(nregioni,nregioni) = eregion2region(nregioni,nregioni) + eself - eself_before
    nregioni = nregionno(nsft+nrelf2a(ib2))
    eregion2region(nregioni,nregioni) = eregion2region(nregioni,nregioni) + eself - eself_before
!
    siteenergy(ib1) = siteenergy(ib1) + eself - eself_before 
    siteenergy(ib2) = siteenergy(ib2) + eself - eself_before 
  enddo
!
!  Use loop over atoms for hardness as it depends on the square of the total charge
!
  do i = 1,neem
    ii = neemptr(i)
    ni = nat(ii)
    eself_before = eself
    qi = qf(ii)**2
    if (lmultiqrange) then
      nqrib1 = nqrnow(i)
    else
      nqrib1 = 1
    endif
    eself = eself + qi*murange(nqrib1,ni,neemtype)
    eself = eself + e0range(nqrib1,ni,neemtype)
!
    nregioni = nregionno(nsft+nrelf2a(ii))
    eregion2region(nregioni,nregioni) = eregion2region(nregioni,nregioni) + eself - eself_before
!
    siteenergy(ii) = siteenergy(ii) + eself - eself_before
  enddo
!*********************************
!  Calculate charge derivatives  *
!*********************************
  if (lgrad2.or.(lgrad1.and.(lfieldcfg(ncf).or.(ntdfieldcfg(ncf).gt.0))).or.ldcharge) then
    call dchargesplit(lmain,emat,numat,ematb,neembond,.false.,.true.)
  endif
!
!  For Pacha we can also compute approximate NMR shifts for relevant nuclei
!
  if (lpacha) then
    allocate(qnmr(numat),stat=status)
    if (status/=0) call outofmemory('eemsplit','qnmr')
    call getnmr(nasym,iatn,qa,qnmr)
  endif
!*******************
!  Output results  *
!*******************
  if ((lmain.or.index(keyword,'debu').ne.0).and.ioproc) then
    if (lqeq) then
      write(ioout,'(//,''  Final charges from split-bond QEq :'',/)')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''    Atom no.            Atomic No.             Charge'')')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      do i = 1,nasym
        write(ioout,'(4x,i6,17x,i3,16x,f10.7)') i,iatn(i),qa(i)
      enddo
    elseif (lpacha) then
      write(ioout,'(//,''  Final charges from split-bond PACHA-EEM :'',/)')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''    Atom no.            Atomic No.             Charge       Chemical shift'')')
      write(ioout,'(''                                                 (e)            (ppm)     '')')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      do i = 1,nasym
        write(ioout,'(4x,i6,17x,i3,16x,f10.7,2x,f11.4)') i,iatn(i),qa(i),qnmr(i)
      enddo
    else
      write(ioout,'(//,''  Final charges from split-bond EEM :'',/)')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''    Atom no.            Atomic No.             Charge'')')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      do i = 1,nasym
        write(ioout,'(4x,i6,17x,i3,16x,f10.7)') i,iatn(i),qa(i)
      enddo
    endif
!
!  Charges for bonds
!
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    write(ioout,'(''    Bond no.            Atom Nos.              Charge'')')
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    do i = 1,neembond
      write(ioout,'(4x,i6,12x,2i6,12x,f10.7)') i,neembonded(1,i),neembonded(2,i),qbond(i)
    enddo
!
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    write(ioout,'(''  Self energy       = '',f16.6,'' eV'')') eself
    write(ioout,'(''--------------------------------------------------------------------------------'')')
  endif
!
!  Free local memory 
!
  if (lpacha) then
    deallocate(qnmr,stat=status)
    if (status/=0) call deallocate_error('eemsplit','qnmr')
  endif
  if (lfieldcfg(ncf).or.(ntdfieldcfg(ncf).gt.0)) then
    deallocate(vfield,stat=status)
    if (status/=0) call deallocate_error('eemsplit','vfield')
  endif
  if (literate) then
    deallocate(oldqf,stat=status)
    if (status/=0) call deallocate_error('eemsplit','oldqf')
  endif
  if (lqiter) then
    deallocate(nlocptr,stat=status)
    if (status/=0) call deallocate_error('eemsplit','nlocptr')
  endif
  deallocate(z2,stat=status)
  if (status/=0) call deallocate_error('eemsplit','z2')
  deallocate(z,stat=status)
  if (status/=0) call deallocate_error('eemsplit','z')
  deallocate(ematb,stat=status)
  if (status/=0) call deallocate_error('eemsplit','ematb')
  deallocate(emat,stat=status)
  if (status/=0) call deallocate_error('eemsplit','emat')
!
!  Exit point if there are no molecules
!
100 continue
  if (lmultiqrange) then
    deallocate(nqrlast,stat=status)
    if (status/=0) call deallocate_error('eemsplit','nqrlast')
  endif
  deallocate(leemfoc,stat=status)
  if (status/=0) call deallocate_error('eemsplit','leemfoc')
!
!  Timing
!
  time2 = g_cpu_time()
  teem = teem + time2 - time1
#ifdef TRACE
  call trace_out('eemsplit')
#endif
!
  return
  end
