  subroutine eem(lmain,lgrad1,lgrad2)
!
!  Subroutine for performing electronegativity equilisation calcns
!  according to work of Mortier
!
!  Periodic boundary condition version
!  Now uses only the asymmetric unit
!
!  lmain = .true.  => call from gulp, one off calc with output
!        = .false. => call from energy, needed for charges
!
!  If second derivatives are being used then tmat must be stored
!  on channel 54 while charges are calculated.
!
!  12/97 Modified to included QEq scheme. Note that when hydrogen
!        is present solution is iterative as interaction depends
!        on the charge in a non-linear fashion.
!  12/97 Self energy calculated and stored for inclusion in total
!        energy.
!  12/97 Modified so that only one of the second derivative arrays
!        has to be passed to make it easier to incorporate into
!        energy call.
!  12/97 Gradient logicals added to call as there are now EEM/QEq
!        contributions to the derivatives
!   1/98 Correction of derivatives removed as this is now handled
!        during main derivative calculation
!   2/98 If lgrad2 and .not.lmain then force dcharge to be used
!        to generate dq/dalpha for the full set of atoms as this 
!        is needed in the second derivative calculation.
!   6/98 Matrix inversion failure trapped
!  12/99 Modified so that unparameterised elements are treated as
!        fixed point charges.
!   2/01 Modified to handle mean field models by reducing sites
!        down to the fully occupied set.
!   5/02 Charges fixed for region 2 in surface calculation
!  10/02 ReaxFF modifications added
!  11/02 Calculation of strain derivatives of charge turned on
!  11/02 Initial charges now set for the benefit of the 1-D case
!        where the convergence of the energy is tested in order
!        to find the range of the Coulomb interaction.
!   2/03 Use of matinv removed to increase speed
!   3/03 dgetrf/i used for lsymopt case since matrix is not symmetric
!   7/05 Streitz and Mintmire modifications added
!   7/05 Array that holds the negative electronegativity now passed
!        to genpot/sympot for correction in S and M scheme
!   5/07 Partial occupancy data moved to module
!   6/07 Size of z array corrected due to requirements of genpot
!   6/07 z array reduced after genpot call for S and M method to required
!        subset of elements
!   3/08 ReaxFF charge scheme removed as this is handled elsewhere
!   6/09 Bug in dimensioning of linear arrays fixed
!   6/09 Site energy added
!   8/11 Electric field added to charge calculation
!   9/11 lgrad1 added to incoming argument list
!  11/11 Region-region energy contributions stored
!   9/12 Pacha added
!  12/12 Option for q0 added
!  12/12 Time-dependent field added
!  12/12 Modified to allow for multiple time-dependent fields
!  11/13 Output format modified to allow for large numbers
!   2/14 Iterative calculation of charges added
!   2/14 Bug in handling of case where charges are fixed corrected
!   4/14 Fixed parameters for iterative charge solve replaced by input values
!   3/15 Option to exclude the Coulomb terms added
!   7/15 External potential added
!   9/16 Matrix inversion for symmetric case moved to subroutine
!   2/17 Blocksize added to call to matrix_inversion_library
!  10/17 Modified so that absolute coordinates are not overwritten for MD
!   1/18 Trace added
!   5/18 oldeem and lelementOK handling moved to seteem
!   5/18 Multiple qranges added
!   6/18 Extra flag passed to dbcgsolve
!   6/18 e0range added
!   6/18 neemrptr set
!   6/18 EEM matrix no longer overwrites derv2/dervi to avoid issues
!        when using numerical derivatives
!   6/18 Correction to setting of emat dimensions
!   7/18 lDoChargeDerv now replaced by ldcharge
!   7/18 emat passed to dcharge routines
!   7/18 Right-hand dimension of emat corrected for case where dcharge is called
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
  integer(i4)                                  :: ifail
  integer(i4)                                  :: ii
  integer(i4)                                  :: ilaenv
  integer(i4), dimension(:), allocatable       :: ipivot
  integer(i4)                                  :: iter
  integer(i4)                                  :: j
  integer(i4)                                  :: jj
  integer(i4)                                  :: k
  integer(i4)                                  :: kk
  integer(i4)                                  :: lwrk
  integer(i4)                                  :: n
  integer(i4),                            save :: ncfold = 0
  integer(i4)                                  :: neemfoc
  integer(i4)                                  :: ni
  integer(i4)                                  :: niter
  integer(i4)                                  :: nitereem
  integer(i4)                                  :: nloc
  integer(i4), dimension(:), allocatable       :: nlocptr
  integer(i4)                                  :: nmax
  integer(i4)                                  :: nmaxu
  integer(i4)                                  :: nqr
  integer(i4), dimension(:), allocatable       :: nqrlast
  integer(i4)                                  :: nr
  integer(i4)                                  :: nregioni
  integer(i4)                                  :: status
  logical,     dimension(:), allocatable       :: leemfoc
  logical                                      :: lconverged
  logical                                      :: ldamp
  logical                                      :: lfound
  logical                                      :: literate
  logical                                      :: lqchange
  logical                                      :: lqiter
  real(dp)                                     :: chii
  real(dp)                                     :: g_cpu_time
  real(dp),    dimension(:,:), allocatable     :: emat
  real(dp)                                     :: enega
  real(dp)                                     :: err
  real(dp)                                     :: eself_before
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
  real(dp),    dimension(:), allocatable       :: wrk
  real(dp)                                     :: zetah0
  real(dp),    dimension(:), allocatable       :: oldqa
  real(dp),    dimension(:), allocatable       :: vfield
  real(dp),    dimension(:), allocatable       :: z
  real(dp),    dimension(:), allocatable       :: z2
#ifdef TRACE
  call trace_in('eem')
#endif
!
!  If this is a parallel run then call distributed memory version
!
  if (nprocs.gt.1) then
    if (lgrad2) then
      call outerror('cannot use EEM for second derivatives in parallel',0_i4)
      call stopnow('eem')
    endif
    call eemd(lmain)
#ifdef TRACE
  call trace_out('eem')
#endif
    return
  endif
!
  time1 = g_cpu_time()
!
!  Allocate local memory
!
  allocate(leemfoc(numat),stat=status)
  if (status/=0) call outofmemory('eem','leemfoc')
  if (lmultiqrange) then
    allocate(nqrlast(numat),stat=status)
    if (status/=0) call outofmemory('eem','nqrlast')
  endif
!
!  Set flag for iterative charges - can't be used for all algorithms
!
  lqiter = literativeQ
  if (lgrad2.or.(lgrad1.and.(lfieldcfg(ncf).or.(ntdfieldcfg(ncf).gt.0))).or.ldcharge) lqiter = .false.
!
  qsum = 0.0_dp
  qtot = 0.0_dp
  neem = 0
  neemrptr(1:numat) = 0
!
!  Check elements and set qrange based on initial charge
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
        call stopnow('eem')
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
        call stopnow('eem')
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
!  Allocate the memory for the main matrix
!
  if (lsymopt.and.lmain) then
    nmaxu = nasym + 1
    nmax = numat + 1
  else
    nmaxu = numat + 1
    nmax = numat + 1
  endif
  allocate(emat(nmax,nmaxu),stat=status)
  if (status/=0) call outofmemory('eem','emat')
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
  if (status/=0) call outofmemory('eem','z')
  allocate(z2(numat),stat=status)
  if (status/=0) call outofmemory('eem','z2')
!
!  If iterative then set up pointers to local elements
!
  if (lqiter) then
    nloc = neem
    allocate(nlocptr(neem+1),stat=status)
    if (status/=0) call outofmemory('eem','nlocptr')
    do i = 1,neem
      nlocptr(i) = i
    enddo
  endif
  if (literate) then
    allocate(oldqa(nasym),stat=status)
    if (status/=0) call outofmemory('eem','oldqa')
  endif
  if (lfieldcfg(ncf).or.(ntdfieldcfg(ncf).gt.0)) then
    allocate(vfield(numat),stat=status)
    if (status/=0) call outofmemory('eem','vfield')
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
    if (index(keyword,'debu').ne.0.and.ioproc) then
      write(ioout,'(/,''  Initial guess for 1-D variable charges :'',/)')
      write(ioout,'('' Atom        Q'')')
      do i = 1,neem
        ii = neemptr(i)
        write(ioout,'(i5,1x,f12.6)') ii,qa(ii)
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
    write(ioout,'(''  Iterative solution of charge equilibration :'',/)')
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
      if (lsymopt) then
        emat(1:numat,1:nasym) = 0.0_dp
      else
        emat(1:numat,1:numat) = 0.0_dp
      endif
    else
      if (lsymopt) then
        call sympot(emat,nmax,z,1_i4)
      else
        call genpot(emat,nmax,z,1_i4)
      endif
    endif
!
!  From S & M, where z has been set without reference to neem, reduce elements to those that are needed
!
    if (lSandM) then
      do i = 1,neem
        z2(i) = z(neemptr(i))
      enddo
      z(1:neem) = z2(1:neem)
    endif
!
!  Reduce to nasym x nasym form
!
    do i = 1,neem
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
!********************************
!  Form matrix of coefficients  *
!********************************
    if (lqeq) then
      do i = 1,neem
        ii = neemptr(i)
        ni = iatn(ii)
        if (lmultiqrange) then
          nqr = nqrnow(i)
        else
          nqr = 1
        endif
        if (ni.ne.1) then
          emat(i,i) = emat(i,i) + 2.0_dp*murange(nqr,ni,neemtype)*occua(ii)
        else
!
!  For hydrogen charge dependant factor must be introduced
!
          zetah0 = 0.529177_dp*0.75_dp/radrange(nqr,1,neemtype)
          rjfac = 1.0_dp + (qa(ii)/zetah0)
          emat(i,i) = emat(i,i) + 2.0_dp*murange(nqr,1,neemtype)*occua(ii)*rjfac
        endif
      enddo
    else
      do i = 1,neem
        ii = neemptr(i)
        ni = iatn(ii)
        if (lmultiqrange) then
          nqr = nqrnow(i)
        else
          nqr = 1
        endif
        emat(i,i) = emat(i,i) + 2.0_dp*murange(nqr,ni,neemtype)*occua(ii)
      enddo
    endif
    emat(nasym+1,nasym+1) = 0.0_dp
    do i = 1,neem
      ii = neemptr(i)
      emat(i,neem+1) = dble(neqv(ii))*occua(ii)
      emat(neem+1,i) = 1.0_dp
    enddo
    emat(neem+1,neem+1) = 0.0_dp
!
!  Add external potential
!
    do i = 1,neem
      ii = neemptr(i)
      z(i) = z(i) - extpotcfg(nsft+ii)
    enddo
!
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
!
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
    if (index(keyword,'debu').ne.0.and.ioproc) then
      write(ioout,'(/,''  EEM/QEq Matrix :'',/)')
      do i = 1,neem + 1
        write(ioout,'(10(1x,f9.5))')(emat(j,i),j=1,neem+1),z(i)
      enddo
    endif
    if (lqiter) then
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
      call dbcgsolve(1_i4,neem+1_i4,neem,nlocptr,emat,nmax,z,qf,qitertol,nqitermax,iter,err)
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
        call outerror('iterative charge solution failed in eem',0_i4)
        call stopnow('eem')
      endif
      do i = 1,neem
        ii = neemptr(i)
        qa(ii) = qf(i)
      enddo
      enega = - qf(neem+1)
    else
!******************
!  Invert matrix  *
!******************
      ifail = 0
      n = neem + 1
      if (lsymopt) then
!*************************
!  Asymmetric inversion  *
!*************************
!     
!  Allocate workspace for inversion
!     
        lwrk = n*ilaenv(1_i4,'DGETRI',' ',n,-1_i4,-1_i4,-1_i4)
        allocate(ipivot(n),stat=status)
        if (status/=0) call outofmemory('eem','ipivot')
        allocate(wrk(lwrk),stat=status)
        if (status/=0) call outofmemory('eem','wrk')
!     
!  Factorise matrix
!     
        call dgetrf(n,n,emat,nmax,ipivot,ifail)
        if (ifail.eq.0) then
!     
!  Form inverse
!     
          call dgetri(n,emat,nmax,ipivot,wrk,lwrk,ifail)
        endif
!     
!  Free workspace  
!     
        deallocate(wrk,stat=status)
        if (status/=0) call deallocate_error('eem','wrk')
        deallocate(ipivot,stat=status)  
        if (status/=0) call deallocate_error('eem','ipivot')
      else
!************************
!  Symmetric inversion  *
!************************
        call matrix_inversion_library(n,1_i4,nmax,nblocksize,emat,0_i4,ifail)
      endif
!
!  Was inversion successful?
!
      if (ifail.ne.0) then
        call outerror('matrix inversion failed in EEM/QEq',0_i4)
        call stopnow('eem')
      endif
!  
!  Multiply inverse matrix and chi matrix to get charges
!
      do i = 1,neem + 1
        ii = neemptr(i)
        qf(ii) = 0.0_dp
        do j = 1,neem + 1
          qf(ii) = qf(ii) + z(j)*emat(j,i)
        enddo
      enddo
      do i = 1,neem
        ii = neemptr(i)
        qa(ii) = qf(ii)
      enddo
      enega = - qf(nasym+1)
    endif
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
!*********************************
!  Calculate charge derivatives  *
!*********************************
  if (lgrad2.or.(lgrad1.and.(lfieldcfg(ncf).or.(ntdfieldcfg(ncf).gt.0))).or.ldcharge) then
    if (lsymopt.and.lmain) then
      call dcharges(lmain,emat,nmax)
    else
      call dcharge(lmain,emat,nmax,lsymopt,.true.)
    endif
  endif
!
!  For Pacha we can also compute approximate NMR shifts for relevant nuclei
!
  if (lpacha) then
    allocate(qnmr(numat),stat=status)
    if (status/=0) call outofmemory('eem','qnmr')
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
          if (niter.eq.1) then
            write(ioout,'(/,''  Charges converged in '',i3,'' iteration'',/)') niter
          else
            write(ioout,'(/,''  Charges converged in '',i3,'' iterations'',/)') niter
          endif
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
    if (status/=0) call deallocate_error('eem','qnmr')
  endif
  if (lfieldcfg(ncf).or.(ntdfieldcfg(ncf).gt.0)) then
    deallocate(vfield,stat=status)
    if (status/=0) call deallocate_error('eem','vfield')
  endif
  if (literate) then
    deallocate(oldqa,stat=status)
    if (status/=0) call deallocate_error('eem','oldqa')
  endif
  if (lqiter) then
    deallocate(nlocptr,stat=status)
    if (status/=0) call deallocate_error('eem','nlocptr')
  endif
  deallocate(z2,stat=status)
  if (status/=0) call deallocate_error('eem','z2')
  deallocate(z,stat=status)
  if (status/=0) call deallocate_error('eem','z')
  deallocate(emat,stat=status)
  if (status/=0) call deallocate_error('eem','emat')
  if (lmultiqrange) then
    deallocate(nqrlast,stat=status)
    if (status/=0) call deallocate_error('eem','nqrlast')
  endif
  deallocate(leemfoc,stat=status)
  if (status/=0) call deallocate_error('eem','leemfoc')
!
!  Timing
!
  time2 = g_cpu_time()
  teem = teem + time2 - time1
#ifdef TRACE
  call trace_out('eem')
#endif
!
  return
  end
