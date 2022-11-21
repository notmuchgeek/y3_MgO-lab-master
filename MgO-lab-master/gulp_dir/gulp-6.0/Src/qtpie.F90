  subroutine qtpie(lmain)
!
!  Subroutine for computing charges according to the QTPIE algorithm.
!
!  lmain = .true.  => call from gulp, one off calc with output
!        = .false. => call from energy, needed for charges
!
!   4/10 Created based on eem
!   7/11 dqtot variable replaced by totalcharge as dqtot was unset
!   9/11 Electric field added to charge calculation
!   9/12 Pacha added
!  12/12 Pacha missing term added
!  12/12 q0 terms added
!  12/12 Time-dependent field added
!  12/12 Modified to allow for multiple time-dependent fields
!   2/14 Same change as per eem made to correct for fixed charge case
!   3/15 lnoqeem added
!   7/15 External potential added
!  10/17 Modified so that absolute coordinates are not overwritten for MD
!   2/18 Trace added
!   5/18 oldeem and lelementOK handling moved to seteem
!   5/18 Multiple qranges added
!   6/18 e0range added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!  12/20 Atom number format increased to allow for larger values
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
!  Julian Gale, CIC, Curtin University, December 2020
!
  use control
  use configurations
  use current
  use derivatives,     only : maxd2, maxd2u, derv2, dervi
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
  logical, intent(in)                          :: lmain
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ifail
  integer(i4)                                  :: ii
  integer(i4)                                  :: j
  integer(i4)                                  :: jj
  integer(i4)                                  :: k
  integer(i4)                                  :: kk
  integer(i4),                            save :: ncfold = 0
  integer(i4)                                  :: neemfoc
  integer(i4)                                  :: ni
  integer(i4)                                  :: niter
  integer(i4)                                  :: nitereem
  integer(i4)                                  :: nqr
  integer(i4), dimension(:), allocatable       :: nqrlast
  integer(i4)                                  :: nr
  integer(i4)                                  :: status
  logical,     dimension(:), allocatable       :: leemfoc
  logical                                      :: lconverged
  logical                                      :: ldamp
  logical                                      :: lfound
  logical                                      :: literate
  logical                                      :: lqchange
  real(dp)                                     :: chii
  real(dp)                                     :: g_cpu_time
  real(dp),    dimension(:), allocatable       :: dEdq
  real(dp),    dimension(:), allocatable       :: qvar
  real(dp)                                     :: ect
  real(dp)                                     :: eself_before
  real(dp)                                     :: q0i
  real(dp)                                     :: qdiff
  real(dp)                                     :: qd
  real(dp)                                     :: qguesstot
  real(dp)                                     :: qi
  real(dp)                                     :: qsum
  real(dp)                                     :: qtot
  real(dp)                                     :: reqv
  real(dp)                                     :: rjfac
  real(dp)                                     :: rmui
  real(dp)                                     :: rnguess
  real(dp)                                     :: time1
  real(dp)                                     :: time2
  real(dp)                                     :: zetah0
  real(dp),    dimension(:), allocatable       :: oldqa
  real(dp),    dimension(:), allocatable       :: vfield
  real(dp),    dimension(:), allocatable       :: z
  real(dp),    dimension(:), allocatable       :: z2
#ifdef TRACE
  call trace_in('qtpie')
#endif
!
  time1 = g_cpu_time()
!
!  Allocate local memory
!
  allocate(leemfoc(numat),stat=status)
  if (status/=0) call outofmemory('qtpie','leemfoc')
  if (lmultiqrange) then
    allocate(nqrlast(numat),stat=status)
    if (status/=0) call outofmemory('qtpie','nqrlast')
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
        call outerror('cannot use QTPIE EEM with shells present',0_i4)
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
        call outerror('cannot use QTPIE EEM with shells present',0_i4)
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
!  Check the memory for the square arrays
!
  if (lsymopt) then
    if (nasym+1.gt.maxd2u) then
      maxd2u = nasym + 1
      call changemaxd2
    endif
    if (numat+1.gt.maxd2) then
      maxd2 = numat + 1
      call changemaxd2
    endif
  else
    if (numat+1.gt.maxd2u) then
      maxd2u = numat + 1
      call changemaxd2
    endif
    if (numat+1.gt.maxd2) then
      maxd2 = numat + 1
      call changemaxd2
    endif
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
!  Allocate memory for the charges and their derivatives
!
  allocate(qvar(neem),stat=status)
  if (status/=0) call outofmemory('qtpie','qvar')
  allocate(dEdq(neem),stat=status)
  if (status/=0) call outofmemory('qtpie','dEdq')
!
!  Allocate local memory that depends on neem
!
  allocate(z(max(neem+1,numat)),stat=status)
  if (status/=0) call outofmemory('qtpie','z')
  allocate(z2(numat),stat=status)
  if (status/=0) call outofmemory('qtpie','z2')
  if (literate) then
    allocate(oldqa(nasym),stat=status)
    if (status/=0) call outofmemory('qtpie','oldqa')
  endif
  if (lfieldcfg(ncf).or.(ntdfieldcfg(ncf).gt.0)) then
    allocate(vfield(numat),stat=status)
    if (status/=0) call outofmemory('qtpie','vfield')
  endif
!
!  Set initial charge state
!
  do i = 1,neem
    qvar(i) = qf(neemptr(i))
  enddo
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
      derv2(1:numat,1:numat) = 0.0_dp
    else
      call genpot(derv2,maxd2,z,1_i4)
    endif
    call screenct(dervi,maxd2)
!
!  From S & M, where z has been set without reference to neem, reduce
!  elements to those that are needed
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
!  Zero storage vector for derv2 array
!
      do j = 1,numat
        z2(j) = 0.0_dp
      enddo
!
!  Place i-j potential terms into derv2
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
          z2(jj) = z2(jj) + derv2(j,neemptr(i))*occuf(j)
        else
          z(i) = z(i) - qf(j)*derv2(j,neemptr(i))*occuf(j)
        endif
      enddo
!
!  Copy temporary storage vector back into derv2 array
!
      do j = 1,neem
        derv2(j,i) = z2(j)
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
          derv2(i,i) = derv2(i,i) + 2.0_dp*murange(nqr,ni,neemtype)*occua(ii)
        else
!
!  For hydrogen charge dependant factor must be introduced
!
          zetah0 = 0.529177_dp*0.75_dp/radrange(nqr,1,neemtype)
          rjfac = 1.0_dp+(qa(ii)/zetah0)
          derv2(i,i) = derv2(i,i) + 2.0_dp*murange(nqr,1,neemtype)*occua(ii)*rjfac
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
        derv2(i,i) = derv2(i,i) + 2.0_dp*murange(nqr,ni,neemtype)*occua(ii)
      enddo
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
        z(i) = z(i) + chirange(nqr,ni,neemtype) - 2.0_dp*murange(nqr,ni,neemtype)*q0range(nqr,ni,neemtype)
      enddo
    else
      do i = 1,neem
        ii = neemptr(i)
        ni = iatn(ii)
        z(i) = z(i) + chirange(1,ni,neemtype) - 2.0_dp*murange(1,ni,neemtype)*q0range(1,ni,neemtype)
      enddo
    endif
    if (lfieldcfg(ncf).or.(ntdfieldcfg(ncf).gt.0)) then
      do i = 1,neem
        ii = neemptr(i)
        z(i) = z(i) - vfield(ii)
      enddo
    endif
    z(neem+1) = totalcharge
!
!  Debugging output
!
    if (index(keyword,'debu').ne.0.and.ioproc) then
      write(ioout,'(/,''  QTPIE: EEM/QEq Matrix :'',/)')
      do i = 1,neem
        write(ioout,'(10(1x,f9.5))')(derv2(j,i),j=1,neem),z(i)
      enddo
    endif
!*********************************************
!  Solve for charge transfer matrix elements *
!*********************************************
    ifail = 0
    call qtmin(neem,qvar,ect,dEdq,maxd2,derv2,dervi,z,ifail)
!  
!  Compute the charges from the charge transfer matrix elements
!
    do i = 1,neem
      ii = neemptr(i)
      qf(ii) = qvar(i)
    enddo
    do i = 1,neem
      ii = neemptr(i)
      qa(ii) = qf(ii)
    enddo
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
    eself = eself + e0range(nqr,ni,neemtype)
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
    siteenergy(i) = siteenergy(i) + eself - eself_before 
  enddo
!*******************
!  Output results  *
!*******************
  if ((lmain.or.index(keyword,'debu').ne.0).and.ioproc) then
    if (lqeq) then
      write(ioout,'(//,''  Final charges from QTPIE QEq :'',/)')
      if (literate) then
        if (lconverged) then
          write(ioout,'(''  Charge transfer converged in '',i3,'' iterations'',/)') niter
        else
          write(ioout,'(''  Failed to converged after '',i3,'' iterations'',/)') nitereem
        endif
      else
        write(ioout,'(''  No hydrogens present - no iteration needed'',/)')
      endif
    else
      write(ioout,'(//,''  Final charges from QTPIE EEM :'',/)')
    endif
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    write(ioout,'(''    Atom no.            Atomic No.             Charge'')')
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    do i = 1,nasym
      write(ioout,'(4x,i6,18x,i2,16x,f10.7)') i,iatn(i),qa(i)
    enddo
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    write(ioout,'(''  Self energy       = '',f16.6,'' eV'')') eself
    write(ioout,'(''--------------------------------------------------------------------------------'')')
  endif
!
!  Free local memory
!
  if (lfieldcfg(ncf).or.(ntdfieldcfg(ncf).gt.0)) then
    deallocate(vfield,stat=status)
    if (status/=0) call deallocate_error('qtpie','vfield')
  endif
  if (literate) then
    deallocate(oldqa,stat=status)
    if (status/=0) call deallocate_error('qtpie','oldqa')
  endif
  deallocate(z2,stat=status)
  if (status/=0) call deallocate_error('qtpie','z2')
  deallocate(z,stat=status)
  if (status/=0) call deallocate_error('qtpie','z')
  deallocate(dEdq,stat=status)
  if (status/=0) call deallocate_error('qtpie','dEdq')
  deallocate(qvar,stat=status)
  if (status/=0) call deallocate_error('qtpie','qvar')
  if (lmultiqrange) then
    deallocate(nqrlast,stat=status)
    if (status/=0) call deallocate_error('qtpie','nqrlast')
  endif
  deallocate(leemfoc,stat=status)
  if (status/=0) call deallocate_error('qtpie','leemfoc')
!
!  Timing
!
  time2 = g_cpu_time()
  teem = teem + time2 - time1
#ifdef TRACE
  call trace_out('qtpie')
#endif
!
  return
  end
