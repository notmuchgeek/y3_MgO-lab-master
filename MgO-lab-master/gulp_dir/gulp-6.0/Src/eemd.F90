  subroutine eemd(lmain)
!
!  Subroutine for performing electronegativity equilisation calcns
!  according to work of Mortier.
!
!  Distributed memory version for parallel case - uses iterative algorithm.
!
!  lmain = .true.  => call from gulp, one off calc with output
!        = .false. => call from energy, needed for charges
!
!   2/14 Created from eem
!   3/14 Unused lgrad1 argument removed
!   4/14 Fixed parameters for iterative charge solve replaced by input values
!   3/15 lnoqeem added
!   7/17 External potential added
!   7/17 Trap added for nloc = 0
!  10/17 Modified so that absolute coordinates are not overwritten for MD
!   1/18 Trace added
!   4/18 Parallel I/O corrected
!   5/18 Correction for S&M case in globalisation of z
!   5/18 oldeem and lelementOK handling moved to seteem
!   5/18 Multiple qranges added
!   6/18 Extra flag passed to dbcgsolve
!   6/18 e0range added
!   6/18 neemrptr set
!   6/18 potl renamed to emat for consistency across routines
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
  use eemdata
  use element
  use energies
  use field,           only : lfieldcfg, ntdfieldcfg
  use iochannels
  use mdlogic,         only : lmd
  use moldyn,          only : labscoany, labsco
  use parallel
  use partial
  use symmetry
  use times
#ifdef TRACE
  use trace,           only : trace_in, trace_out
#endif
  implicit none
#ifdef MPI
  include 'mpif.h'
#endif
!
!  Passed variables
!
  logical, intent(in)                          :: lmain
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ii
  integer(i4)                                  :: il
  integer(i4)                                  :: iopr
  integer(i4)                                  :: iter
  integer(i4)                                  :: j
  integer(i4)                                  :: jj
  integer(i4)                                  :: k
  integer(i4)                                  :: kk
  integer(i4)                                  :: maxp
  integer(i4),                            save :: ncfold = 0
  integer(i4)                                  :: neemfoc
  integer(i4)                                  :: ni
  integer(i4)                                  :: niter
  integer(i4)                                  :: nitereem
  integer(i4)                                  :: nloc
  integer(i4), dimension(:), allocatable       :: nlocptr
  integer(i4)                                  :: nqr
  integer(i4), dimension(:), allocatable       :: nqrlast
  integer(i4)                                  :: nr
  integer(i4)                                  :: nregioni
  integer(i4)                                  :: status
#ifdef MPI
  integer                                      :: MPIerror
  integer                                      :: ntag
  integer                                      :: nnode
  integer                                      :: ntmp
  integer                                      :: Request
  integer(i4), dimension(:), allocatable       :: StatMPI       ! Array for status from MPI
#endif
  logical,     dimension(:), allocatable       :: leemfoc
  logical                                      :: ldamp
  logical                                      :: lconverged
  logical                                      :: lfound
  logical                                      :: literate
  logical                                      :: lqchange
  real(dp)                                     :: chii
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: enega
  real(dp)                                     :: err
  real(dp)                                     :: eself_before
  real(dp),    dimension(:,:), allocatable     :: emat
  real(dp)                                     :: q0i
  real(dp)                                     :: qdiff
  real(dp)                                     :: qd
  real(dp)                                     :: qguesstot
  real(dp)                                     :: qi
  real(dp),    dimension(:), allocatable       :: qnmr
  real(dp)                                     :: qsum
  real(dp)                                     :: qtot
  real(dp)                                     :: reqv
  real(dp)                                     :: rjfac
  real(dp)                                     :: rmui
  real(dp)                                     :: rnguess
  real(dp)                                     :: time1
  real(dp)                                     :: time2
#ifdef MPI
  real(dp),    dimension(:), allocatable       :: tmp
  real(dp),    dimension(:), allocatable       :: tmp2
#endif
  real(dp)                                     :: zetah0
  real(dp),    dimension(:), allocatable       :: oldqa
  real(dp),    dimension(:), allocatable       :: vfield
  real(dp),    dimension(:), allocatable       :: z
  real(dp),    dimension(:), allocatable       :: z2
#ifdef TRACE
  call trace_in('eemd')
#endif
!
  time1 = g_cpu_time()
!
!  Allocate local memory
!
  allocate(leemfoc(numat),stat=status)
  if (status/=0) call outofmemory('eemd','leemfoc')
  if (lmultiqrange) then
    allocate(nqrlast(numat),stat=status)
    if (status/=0) call outofmemory('eemd','nqrlast')
  endif
!
  qsum = 0.0_dp
  qtot = 0.0_dp
  neem = 0
  neemrptr(1:numat) = 0
!
!  Check elements
!
  if (lsymopt) then
    do i = 1,nasym
      ni = iatn(i)
      qi = qa(i)
      if (lelementOK(ni).and.nregionno(nsft+i).eq.1) then
        neem = neem + 1
        neemptr(neem) = i
        neemrptr(i) = neem
        qsum = qsum - neqv(i)*qi*occua(i)
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
        call stopnow('eemd')
      else
        qtot = qtot + neqv(i)*qi*occua(i)
      endif
    enddo
  else
    do i = 1,numat
      ni = nat(i)
      qi = qf(i)
      if (lelementOK(ni).and.nregionno(nsft+nrelf2a(i)).eq.1) then
        neem = neem + 1
        neemptr(neem) = i
        neemrptr(i) = neem
        qsum = qsum - qi*occuf(i)
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
        call stopnow('eemd')
      else
        qtot = qtot + qi*occuf(i)
      endif
    enddo
  endif
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
!  Check the memory for the linear arrays
!
  if (numat+1.gt.maxat) then
    maxat = numat + 1
    call changemaxat
  endif
!
!  Set the pointer to where the electronegativity should be as well
!
  neemptr(neemfoc+1) = nasym + 1
!
!  Decide on total EEM fragment charge - for cluster
!  there is no constraint so use sum of initial charges.
!  For periodic system, charge must be equal to the
!  negative sum of the non-EEM ion charges.
!
  if (ndim.eq.0) then
    qtot = qsum
  endif
!*****************************************************************
!  Is hydrogen present in QEq? If so then solution is iterative  *
!*****************************************************************
  literate = lmultiqrange
  ldamp = .false.
  if (lqeq.and..not.literate) then
    i = 0
    do while (i.lt.nasym.and..not.literate)
      i = i + 1
      literate = (iatn(i).eq.1)
      if (literate) ldamp = .true.
    enddo
  endif
  if (literate) then
    nitereem = nqeqitermax
  else
    nitereem = 1
  endif
!
!  Allocate local memory that depends on neem
!
  allocate(z(max(neem+1,numat)),stat=status)
  if (status/=0) call outofmemory('eemd','z')
  allocate(z2(numat),stat=status)
  if (status/=0) call outofmemory('eemd','z2')
  allocate(nlocptr(neem+1),stat=status)
  if (status/=0) call outofmemory('eemd','nlocptr')
!
!  Set up pointers to local elements for parallelisation
!
  if (nprocs.eq.1) then
    nloc = neem
    do i = 1,neem
      nlocptr(i) = i
    enddo
  else
    nloc = 0
    do i = procid+1,neem,nprocs
      nloc = nloc + 1
      nlocptr(nloc) = i
    enddo
  endif
!
!  Allocate the memory for the potential array
!
!  procid = 0 carries the extra row for the constraints
!
  maxp = numat + 1
  if (procid.eq.0) then
    allocate(emat(maxp,nloc+1),stat=status)
    if (status/=0) call outofmemory('eemd','emat')
  else
    allocate(emat(maxp,nloc),stat=status)
    if (status/=0) call outofmemory('eemd','emat')
  endif
!
  if (literate) then
    allocate(oldqa(nasym),stat=status)
    if (status/=0) call outofmemory('eemd','oldqa')
  endif
  if (lfieldcfg(ncf).or.(ntdfieldcfg(ncf).gt.0)) then
    allocate(vfield(numat),stat=status)
    if (status/=0) call outofmemory('eemd','vfield')
  endif
!
!  Calculate reciprocal lattice vectors
!
  if (lewald.and.ndim.gt.1) call kindex
!
!  For 1-D case set guess at charges based EEM parameters
!  and scaled up by 1.5 to allow for increase in ionicity.
!
  if (ndim.eq.1.and.ncf.ne.ncfold) then
    ncfold = ncf
    qguesstot = 0.0_dp
    rnguess = 0.0_dp
    do i = 1,neem
      ii = neemptr(i)
      ni = iatn(ii)
      if (lmultiqrange) then
        nqr = nqrnow(i)
      else
        nqr = 1
      endif
      chii = chirange(nqr,ni,neemtype)
      rmui = murange(nqr,ni,neemtype)
      q0i  = q0range(nqr,ni,neemtype)
      qa(ii) = q0i - chii/rmui
      qguesstot = qguesstot + qa(ii)*dble(neqv(ii))*occua(ii)
      rnguess = rnguess + dble(neqv(ii))*occua(ii)
    enddo
    qguesstot = (qguesstot + qtot)/rnguess
    do i = 1,neem
      ii = neemptr(i)
      qa(ii) = qa(ii) - qguesstot
      if (abs(qtot).lt.1.0d-12) qa(ii) = 1.5_dp*qa(ii)
      do j = 1,numat
        if (nrelf2a(j).eq.ii) then
          qf(j) = qa(ii)
        endif
      enddo
    enddo
    if (index(keyword,'debu').ne.0) then
      if (ioproc) then
        write(ioout,'(/,''  Initial guess for 1-D variable charges :'',/)')
        write(ioout,'('' Atom        Q'')')
      endif
#ifdef MPI
      if (lioproconly) then
        allocate(tmp(numat),stat=status)
        if (status/=0) call outofmemory('eemd','tmp')
        allocate(tmp2(numat),stat=status)
        if (status/=0) call outofmemory('eemd','tmp2')
!
        tmp(1:numat) = 0.0_dp
!
        iopr = 0
        do i = 1,neem
          if (iopr.eq.procid) then
            ii = neemptr(i)
            tmp(ii) = qa(ii)
          endif
          iopr = iopr + 1
          iopr = mod(iopr,nprocs)
        enddo
!
        call sumall(tmp,tmp2,numat,"eemd","qa")
!
        if (ioproc) then
          do i = 1,neem
            ii = neemptr(i)
            write(ioout,'(i5,1x,f12.6)') ii,tmp2(ii)
          enddo
          write(ioout,'(/)')
        endif
!
        deallocate(tmp2,stat=status)
        if (status/=0) call deallocate_error('eemd','tmp2')
        deallocate(tmp,stat=status)
        if (status/=0) call deallocate_error('eemd','tmp')
      else
#endif
        iopr = 0
        do i = 1,neem
          if (iopr.eq.procid) then
            ii = neemptr(i)
            write(ioout,'(i5,1x,f12.6)') ii,qa(ii)
          endif
          iopr = iopr + 1
          iopr = mod(iopr,nprocs)
        enddo
        if (ioproc) then
          write(ioout,'(/)')
        endif
#ifdef MPI
      endif
#endif
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
    do i = 1,nasym
      oldqa(i) = qa(i)
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
    write(ioout,'(''  Iterative solution of QEq :'',/)')
  endif
  do while (niter.lt.nitereem.and..not.lconverged)
    niter = niter + 1
!
!  Zero right hand vector
!
    z(1:neem) = 0.0_dp
!*******************************************
!  Generate electrostatic potential terms  *
!*******************************************
    if (lnoqeem) then
      emat(1:numat,1:nloc) = 0.0_dp
    else
      if (lsymopt) then
        call sympotd(nloc,nlocptr,emat,maxp,z)
      else
        call genpotd(nloc,nlocptr,emat,maxp,z)
      endif
    endif
!
!  From S & M, where z has been set without reference to neem, reduce elements to those that are needed
!
    if (lSandM) then
!
!  Globalise z
!
      call sumall(z,z2,numat,"eemd","z")
      do i = 1,neem
        z(i) = z2(neemptr(i))
      enddo
    else
      call sumall(z,z2,neem,"eemd","z")
      z(1:neem) = z2(1:neem)
    endif
!
!  Reduce to nasym x nasym form
!
    do il = 1,nloc
      i = nlocptr(il)
!
!  Zero storage vector for emat array
!
      do j = 1,numat
        z2(j) = 0.0_dp
      enddo
!
!  Place i-j potential terms into emat
!
      do j = 1,numat
        k = nrelf2a(j)
        jj = 1
        kk = neemptr(jj)
        do while (jj.lt.neem.and.kk.ne.k)
          jj = jj + 1
          kk = neemptr(jj)
        enddo
!
!  Variable j charge case
!
        if (kk.eq.k) then
          z2(jj) = z2(jj) + emat(j,il)*occuf(j)
        else
          z(i) = z(i) - qf(j)*emat(j,il)*occuf(j)
        endif
      enddo
!
!  Copy temporary storage vector back into emat array
!
      do j = 1,neem
        emat(j,il) = z2(j)
      enddo
    enddo
!********************************
!  Form matrix of coefficients  *
!********************************
    if (nloc.gt.0) then
      if (lqeq) then
        do il = 1,nloc
          i = nlocptr(il)
          ii = neemptr(i)
          ni = iatn(ii)
          if (lmultiqrange) then
            nqr = nqrnow(i)
          else
            nqr = 1
          endif
          if (ni.ne.1) then
            emat(i,il) = emat(i,il) + 2.0_dp*murange(nqr,ni,neemtype)*occua(ii)
          else
!
!  For hydrogen charge dependant factor must be introduced
!
            zetah0 = 0.529177_dp*0.75_dp/radrange(nqr,1,neemtype)
            rjfac = 1.0_dp + (qa(ii)/zetah0)
            emat(i,il) = emat(i,il) + 2.0_dp*murange(nqr,1,neemtype)*occua(ii)*rjfac
          endif
        enddo
      else
        do il = 1,nloc
          i = nlocptr(il)
          ii = neemptr(i)
          ni = iatn(ii)
          if (lmultiqrange) then
            nqr = nqrnow(i)
          else
            nqr = 1
          endif
          emat(i,il) = emat(i,il) + 2.0_dp*murange(nqr,ni,neemtype)*occua(ii)
        enddo
      endif
      emat(nasym+1,nloc) = 0.0_dp
      do il = 1,nloc
        i = nlocptr(il)
        emat(neem+1,il) = 1.0_dp
      enddo
      if (procid.eq.0) then
        do i = 1,neem
          ii = neemptr(i)
          emat(i,nloc+1) = dble(neqv(ii))*occua(ii)
        enddo
        emat(neem+1,nloc+1) = 0.0_dp
      endif
    endif
!
!  Add external potential
!
    do i = 1,neem
      ii = neemptr(i)
      z(i) = z(i) - extpotcfg(nsft+ii)
    enddo
    if (lmultiqrange) then
      do i = 1,neem
        ii = neemptr(i)
        ni = iatn(ii)
        nqr = nqrnow(i)
        z(i) = z(i) - chirange(nqr,ni,neemtype) + 2.0_dp*murange(nqr,ni,neemtype)*q0range(nqr,ni,neemtype)
      enddo
    else
      do i = 1,neem
        ii = neemptr(i)
        ni = iatn(ii)
        z(i) = z(i) - chirange(1,ni,neemtype) + 2.0_dp*murange(1,ni,neemtype)*q0range(1,ni,neemtype)
      enddo
    endif
    if (lfieldcfg(ncf).or.(ntdfieldcfg(ncf).gt.0)) then
      do i = 1,neem
        ii = neemptr(i)
        z(i) = z(i) - vfield(ii)
      enddo
    endif
    z(neem+1) = - qtot
!
!  Debugging output
!
    if (index(keyword,'debu').ne.0) then
      if (ioproc) then
        write(ioout,'(/,''  EEM/QEq Matrix :'',/)')
      endif
#ifdef MPI
      if (lioproconly) then
!
!  Allocate temporary workspace for communication
!
        ntmp = neem + 1
        ntag = 1
        allocate(tmp(ntmp),stat=status)
        if (status/=0) call outofmemory('eemd','tmp')
        allocate(StatMPI(MPI_Status_Size),stat=status)
        if (status/=0) call outofmemory('eemd','StatMPI')
!
        iopr = 0
        il = 0
        do i = 1,neem
          if (iopr.ne.0_i4) then
!
!  Post receive
!
            if (ioproc) then
              nnode = iopr
              call MPI_IRecv(tmp,ntmp,MPI_double_precision,nnode, &
                             ntag,MPI_Comm_World,Request,MPIerror)
            endif
!
!  Pass data to ioproc for writing
!
            if (procid.eq.iopr) then
              il = il + 1
              tmp(1:neem+1) = emat(1:neem+1,il)
!
!  Post send
!
              call MPI_ISend(tmp,ntmp,MPI_double_precision,0, &
                             ntag,MPI_Comm_World,Request,MPIerror)
            endif
            if (ioproc.or.procid.eq.iopr) then
              call MPI_WaitAll(1,Request,StatMPI,MPIerror)
            endif
            if (ioproc) then
!
!  Write on I/O node
!
              write(ioout,'(i5,10(1x,f9.5))') i,(tmp(j),j=1,neem+1),z(i)
            endif
          elseif (ioproc) then
            il = il + 1
            write(ioout,'(i5,10(1x,f9.5))') i,(emat(j,il),j=1,neem+1),z(i)
          endif
          iopr = iopr + 1
          iopr = mod(iopr,nprocs)
        enddo
        if (procid.eq.0) then
          write(ioout,'(i5,10(1x,f9.5))') i,(emat(j,nloc+1),j=1,neem+1),z(neem+1)
        endif
!
        deallocate(StatMPI,stat=status)
        if (status/=0) call deallocate_error('eemd','StatMPI')
        deallocate(tmp,stat=status)
        if (status/=0) call deallocate_error('eemd','tmp')
      else
#endif
        iopr = 0
        il = 0
        do i = 1,neem
          if (procid.eq.iopr) then
            il = il + 1
            write(ioout,'(i5,10(1x,f9.5))') i,(emat(j,il),j=1,neem+1),z(i)
          endif
          iopr = iopr + 1
          iopr = mod(iopr,nprocs)
        enddo
        if (procid.eq.0) then
          write(ioout,'(i5,10(1x,f9.5))') i,(emat(j,nloc+1),j=1,neem+1),z(neem+1)
        endif
#ifdef MPI
      endif
#endif
    endif
!***************************
!  Iterative charge solve  *
!***************************
    do i = 1,neem
      ii = neemptr(i)
      qf(i) = qa(ii)
    enddo
    qf(neem+1) = 0.0_dp
!
!  Solve using iterative route
!
    call dbcgsolve(1_i4,neem+1_i4,nloc,nlocptr,emat,maxp,z,qf,qitertol,nqitermax,iter,err)
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
      call outerror('iterative charge solution failed in eemd',0_i4)
      call stopnow('eemd')
    endif
    do i = 1,neem
      ii = neemptr(i)
      qa(ii) = qf(i)
    enddo
    enega = - qf(neem+1)
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
          ni = iatn(ii)
          if (lsymopt) then
            qi = qa(ii)
          else
            qi = qf(ii)
          endif
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
      do i = 1,nasym
        qd = qa(i) - oldqa(i)
        qdiff = qdiff + abs(qd)
      enddo
      qdiff = qdiff/dble(nasym)
      lconverged = (qdiff.lt.qeqscfcrit.and..not.lqchange)
      if (lmain.and.ioproc) then
        write(ioout,'(''  ** Cycle : '',i4,'' Qdiff : '',f10.8)') niter,qdiff
      endif
      if (.not.lconverged) then
        if (ldamp) then
!
!  Damp change to improve convergence if QEq has triggered iteration
!
          do i = 1,neem
            ii = neemptr(i)
            qd = qa(ii) - oldqa(ii)
            qa(ii) = qa(ii) - 0.25_dp*qd
            oldqa(ii) = qa(ii)
          enddo
        else
          do i = 1,neem
            ii = neemptr(i)
            oldqa(ii) = qa(ii)
          enddo
        endif
      endif
    endif
!
!  Transfer charges to qf
!
    do i = 1,numat
      nr = nrelf2a(i)
      qf(i) = qa(nr)
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
  do i = 1,neem
    ii = neemptr(i)
    qi = qa(ii)
    ni = iatn(ii)
    reqv = dble(neqv(ii))*occua(ii)
    eself_before = eself
!
    if (lmultiqrange) then
      nqr = nqrnow(i)
    else
      nqr = 1
    endif
!
    eself = eself + reqv*e0range(nqr,ni,neemtype)
    if (lqeq) then
      q0i = q0range(nqr,ni,neemtype)
      if (ni.ne.1) then
        eself = eself + (qi-q0i)*reqv*(chirange(nqr,ni,neemtype)+(qi-q0i)*murange(nqr,ni,neemtype))
      else
        zetah0 = 0.529177_dp*0.75_dp/radrange(nqr,1,neemtype)
        eself = eself + (qi-q0i)*reqv*(chirange(nqr,1,neemtype)+(qi-q0i)*murange(nqr,1,neemtype)* &
                        (1.0_dp+(2.0_dp*(qi-q0i)/(3.0_dp*zetah0))))
      endif
    else
      q0i = q0range(nqr,ni,neemtype)
      eself = eself + (qi-q0i)*reqv*(chirange(nqr,ni,neemtype)+(qi-q0i)*murange(nqr,ni,neemtype))
    endif
!
!  Add external potential for site
!
    eself = eself + qi*reqv*extpotcfg(nsft+ii)
!
    nregioni = nregionno(nsft+ii)
    eregion2region(nregioni,nregioni) = eregion2region(nregioni,nregioni) + eself - eself_before
!
    siteenergy(ii) = siteenergy(ii) + eself - eself_before 
  enddo
!
!  For Pacha we can also compute approximate NMR shifts for relevant nuclei
!
  if (lpacha) then
    allocate(qnmr(numat),stat=status)
    if (status/=0) call outofmemory('eemd','qnmr')
    call getnmr(nasym,iatn,qa,qnmr)
  endif
!*******************
!  Output results  *
!*******************
  if ((lmain.or.index(keyword,'debu').ne.0).and.ioproc) then
    if (lqeq) then
      write(ioout,'(//,''  Final charges from QEq :'',/)')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''    Atom no.            Atomic No.             Charge'')')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      do i = 1,nasym
        write(ioout,'(4x,i6,17x,i3,16x,f10.7)') i,iatn(i),qa(i)
      enddo
    elseif (lpacha) then
      write(ioout,'(//,''  Final charges from PACHA-EEM :'',/)')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''    Atom no.            Atomic No.             Charge       Chemical shift'')')
      write(ioout,'(''                                                 (e)            (ppm)     '')')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      do i = 1,nasym
        write(ioout,'(4x,i6,17x,i3,16x,f10.7,2x,f11.4)') i,iatn(i),qa(i),qnmr(i)
      enddo
    else
      write(ioout,'(//,''  Final charges from EEM :'',/)')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''    Atom no.            Atomic No.             Charge'')')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      do i = 1,nasym
        write(ioout,'(4x,i6,17x,i3,16x,f10.7)') i,iatn(i),qa(i)
      enddo
    endif
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    write(ioout,'(''  Electronegativity = '',f16.6,'' eV'')') enega
    write(ioout,'(''  Self energy       = '',f16.6,'' eV'')') eself
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    if (lqeq) then
      if (literate) then
        if (lconverged) then
          write(ioout,'(/,''  Charges converged in '',i3,'' iterations'',/)') niter
        else
          write(ioout,'(/,''  Failed to converged after '',i3,'' iterations'',/)') nitereem
        endif
      else
        write(ioout,'(/,''  No hydrogens present - no iteration needed'',/)')
      endif
    endif
  endif
!
!  Free local memory 
!
  if (lpacha) then
    deallocate(qnmr,stat=status)
    if (status/=0) call deallocate_error('eemd','qnmr')
  endif
  if (lfieldcfg(ncf).or.(ntdfieldcfg(ncf).gt.0)) then
    deallocate(vfield,stat=status)
    if (status/=0) call deallocate_error('eemd','vfield')
  endif
  if (literate) then
    deallocate(oldqa,stat=status)
    if (status/=0) call deallocate_error('eemd','oldqa')
  endif
  deallocate(emat,stat=status)
  if (status/=0) call deallocate_error('eemd','emat')
  deallocate(nlocptr,stat=status)
  if (status/=0) call deallocate_error('eemd','nlocptr')
  deallocate(z2,stat=status)
  if (status/=0) call deallocate_error('eemd','z2')
  deallocate(z,stat=status)
  if (status/=0) call deallocate_error('eemd','z')
  if (lmultiqrange) then
    deallocate(nqrlast,stat=status)
    if (status/=0) call deallocate_error('eemd','nqrlast')
  endif
  deallocate(leemfoc,stat=status)
  if (status/=0) call deallocate_error('eemd','leemfoc')
!
!  Timing
!
  time2 = g_cpu_time()
  teem = teem + time2 - time1
#ifdef TRACE
  call trace_out('eemd')
#endif
!
  return
  end
