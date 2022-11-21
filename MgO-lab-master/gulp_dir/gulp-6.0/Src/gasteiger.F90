  subroutine gasteiger(lmain)
!
!  Subroutine for calculating charges according to the scheme
!  of Gasteiger and Marsili.
!
!  lmain = .true.  => call from gulp, one off calc with output
!        = .false. => call from energy, needed for charges
!
!   1/07 Created based on eem.f
!   4/07 Searching for bonded atoms now uses nbonds number rather
!        than checking for a non-zero type
!   5/07 Partial occupancy data moved to module
!  12/07 References to nitergast replaced with niter
!  12/07 Unused variables removed
!   3/15 Modified to allow Gasteiger charges for non-neutral
!        systems
!   3/15 Search for species specific parameters added
!   3/15 Modified to allow for input of damping options
!   3/15 Modified so that charges are always initialised from the
!        same initial value to ensure consistency.
!   3/15 nbonds used directly in hybridisation check
!   1/18 Trace added
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
  use element
  use energies
  use iochannels
  use parallel
  use partial
  use species,          only : lgastinspec, gastspec
  use symmetry
  use times
#ifdef TRACE
  use trace,            only : trace_in, trace_out
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
  integer(i4)                                  :: ia
  integer(i4)                                  :: imm
  integer(i4)                                  :: ii
  integer(i4)                                  :: j
  integer(i4)                                  :: nbond
  integer(i4)                                  :: ngast
  integer(i4)                                  :: ngastfoc
  integer(i4)                                  :: ni
  integer(i4)                                  :: niter
  integer(i4)                                  :: nri
  integer(i4), dimension(:), allocatable       :: ngastptr
  integer(i4), dimension(:), allocatable       :: ngastrptr
  integer(i4)                                  :: status
  logical,     dimension(:), allocatable       :: lgastfoc
  logical                                      :: lconverged
  real(dp)                                     :: damping
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: qdiff
  real(dp)                                     :: qd
  real(dp)                                     :: qsum
  real(dp)                                     :: time1
  real(dp)                                     :: time2
  real(dp),    dimension(:), allocatable       :: chiG
  real(dp),    dimension(:), allocatable       :: chiGplus
  real(dp),    dimension(:), allocatable       :: gastA
  real(dp),    dimension(:), allocatable       :: gastB
  real(dp),    dimension(:), allocatable       :: gastC
  real(dp),    dimension(:), allocatable       :: oldqa
#ifdef TRACE
  call trace_in('gasteiger')
#endif
!
  time1 = g_cpu_time()
!
!  Allocate local memory
!
  allocate(ngastptr(numat),stat=status)
  if (status/=0) call outofmemory('gasteiger','ngastptr')
  allocate(ngastrptr(numat),stat=status)
  if (status/=0) call outofmemory('gasteiger','ngastrptr')
  allocate(lgastfoc(numat),stat=status)
  if (status/=0) call outofmemory('gasteiger','lgastfoc')
  allocate(chiG(numat),stat=status)
  if (status/=0) call outofmemory('gasteiger','chiG')
  allocate(chiGplus(numat),stat=status)
  if (status/=0) call outofmemory('gasteiger','chiGplus')
  allocate(gastA(numat),stat=status)
  if (status/=0) call outofmemory('gasteiger','gastA')
  allocate(gastB(numat),stat=status)
  if (status/=0) call outofmemory('gasteiger','gastB')
  allocate(gastC(numat),stat=status)
  if (status/=0) call outofmemory('gasteiger','gastC')
  allocate(oldqa(numat),stat=status)
  if (status/=0) call outofmemory('gasteiger','oldqa')
!
!  Set up list of allowed elements
!
  do i = 1,maxele
    if (abs(gasteigerA(i)).gt.1.0d-6.or.abs(gasteigerB(i)).gt.1.0d-6) then
      lelementOK(i) = .true.
    else
      lelementOK(i) = .false.
    endif
  enddo
!
  ngast = 0
  ngastrptr(1:nasym) = 0
!
!  Check elements
!
  if (lsymopt) then
    do i = 1,nasym
      ia = iatn(i)
      if (lelementOK(ia).and.nregionno(nsft+i).eq.1) then
        ngast = ngast + 1
        ngastptr(ngast) = i
        ngastrptr(i) = ngast
      elseif (ia.gt.maxele) then
        call outerror('cannot use Gasteiger with shells present',0_i4)
        call stopnow('gasteiger')
      endif
    enddo
  else
    do i = 1,numat
      ia = nat(i)
      if (lelementOK(ia).and.nregionno(nsft+nrelf2a(i)).eq.1) then
        ngast = ngast + 1
        ngastptr(ngast) = i
        ngastrptr(i) = ngast
      elseif (ia.gt.maxele) then
        call outerror('cannot use Gasteiger with shells present',0_i4)
        call stopnow('gasteiger')
      endif
    enddo
  endif
!
!  Now find the number of fully occupied sites for Gasteiger
!
  lgastfoc(1:numat) = .false.
  do i = 1,ngast
    ii = iocptr(ngastptr(i))
    lgastfoc(ii) = .true.
  enddo
  ngastfoc = 0
  do i = 1,ncfoc
    if (lgastfoc(i)) ngastfoc = ngastfoc + 1
  enddo
!
!  Assign parameters - take defaults from main arrays and then check hybridisation of C, N & O
!
  do i = 1,ngast
    ii = ngastptr(i)
    ni = iatn(ii)
    nri = nrela2f(ii)
    gastA(i) = gasteigerA(ni)
    gastB(i) = gasteigerB(ni)
    gastC(i) = gasteigerC(ni)
    qa(i)    = 0.0_dp
!
!  Get number of bonds for C, N or O
!
    nbond = nbonds(nri)
!
!  Modify parameters for C, N or O if not sp3 hybridised
!
    if (ni.eq.6.and.nbond.eq.3) then
      gastA(i) =  8.79_dp
      gastB(i) =  9.32_dp
      gastC(i) =  1.51_dp
    elseif (ni.eq.6.and.nbond.le.2) then
      gastA(i) = 10.39_dp
      gastB(i) =  9.45_dp
      gastC(i) =  0.73_dp
    elseif (ni.eq.7.and.nbond.eq.2) then
      gastA(i) = 12.87_dp
      gastB(i) = 11.15_dp
      gastC(i) =  0.85_dp
    elseif (ni.eq.7.and.nbond.eq.1) then
      gastA(i) = 15.68_dp
      gastB(i) = 11.70_dp
      gastC(i) = -0.27_dp
    elseif (ni.eq.8.and.nbond.eq.1) then
      gastA(i) = 17.07_dp
      gastB(i) = 13.79_dp
      gastC(i) =  0.47_dp
    endif
!
!  Check for species specific parameters
!
    if (lgastinspec(nspecptr(ii))) then
      gastA(i) = gastspec(1,nspecptr(ii))
      gastB(i) = gastspec(2,nspecptr(ii))
      gastC(i) = gastspec(3,nspecptr(ii))
      qa(i)    = gastspec(4,nspecptr(ii))
    endif
!
!  Set electronegativity of plus one ion - doesn't depend on charge!
!  Exception is hydrogen where a special value is used.
!
    if (ni.eq.1) then
      chiGplus(i) = 20.02_dp
    else
      chiGplus(i) = gastA(i) + gastB(i) + gastC(i)
    endif
  enddo
!
!  Save old charges and store for convergence check
!
  do i = 1,nasym
    oldqa(i) = qa(i)
  enddo
!****************************
!  Start of iterative loop  *
!****************************
  lconverged = .false.
  niter = 0
  if (lmain.and.ioproc) then
    write(ioout,'(''  Iterative solution of Gasteiger charges :'',/)')
  endif
  damping = gastdamp
  do while (niter.lt.ngastitermax.and..not.lconverged)
    niter = niter + 1
    if (ngastdamptype.eq.3) then
!
!  Accelerate damping
!
      damping = 1.0_dp - gastdamp**niter
    elseif (ngastdamptype.eq.1) then
!
!  Original damping
!
      damping = damping*gastdamp
    endif
!***********************************************
!  Find electronegativities at this iteration  *
!***********************************************
    do i = 1,ngast
      ii = ngastptr(i)
      chiG(i) = gastA(i) + qa(ii)*(gastB(i) + qa(ii)*gastC(i))
    enddo
!****************************************
!  Compute charges at latest iteration  *
!****************************************
    do i = 1,ngast
      ii = ngastptr(i)
      nri = nrela2f(ii)
!
!  Loop over bonds and compute contributions to charge
!
      nbond = 0
      imm = 1
      do while (imm.gt.0.and.nbond.lt.nbonds(nri))
        nbond = nbond + 1
        imm = nbonded(nbond,nri)
        if (imm.gt.0) then
          j = ngastrptr(nrelf2a(imm))
          if (j.gt.0) then
            if (chiG(j).gt.chiG(i)) then
              qa(ii) = qa(ii) + damping*(chiG(j)-chiG(i))/chiGplus(i)
            else
              qa(ii) = qa(ii) + damping*(chiG(j)-chiG(i))/chiGplus(j)
            endif
          endif
        endif
      enddo
    enddo
!**************************
!  Check for convergence  *
!**************************
    qdiff = 0.0_dp
    do i = 1,nasym
      qd = qa(i) - oldqa(i)
      qdiff = qdiff + abs(qd)*neqv(i)
      oldqa(i) = qa(i)
    enddo
    qdiff = qdiff/dble(numat)
    lconverged = (qdiff.lt.gasttol)
    if (lmain.and.ioproc) then
      write(ioout,'(''  ** Cycle : '',i4,'' Qdiff : '',f10.8)') niter,qdiff
    endif
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
!*******************
!  Output results  *
!*******************
  if (lmain.and.ioproc) then
    write(ioout,'(//,''  Final Gasteiger charges:'',/)')
    if (lconverged) then
      write(ioout,'(''  Charges converged in '',i3,'' iterations'',/)') niter
    else
      write(ioout,'(''  Failed to converged after '',i3,'' iterations'',/)') niter
    endif
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    write(ioout,'(''    Atom no.            Atomic No.             Charge'')')
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    qsum = 0.0_dp
    do i = 1,nasym
      write(ioout,'(6x,i4,18x,i2,16x,f10.7)') i,iatn(i),qa(i)
      qsum = qsum + qa(i)*dble(neqv(i))
    enddo
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    write(ioout,'(''    Total '',36x,f10.7)') qsum
    write(ioout,'(''--------------------------------------------------------------------------------'')')
  endif
!
!  Free local memory 
!
  deallocate(oldqa,stat=status)
  if (status/=0) call deallocate_error('gasteiger','oldqa')
  deallocate(gastC,stat=status)
  if (status/=0) call deallocate_error('gasteiger','gastC')
  deallocate(gastB,stat=status)
  if (status/=0) call deallocate_error('gasteiger','gastB')
  deallocate(gastA,stat=status)
  if (status/=0) call deallocate_error('gasteiger','gastA')
  deallocate(chiGplus,stat=status)
  if (status/=0) call deallocate_error('gasteiger','chiGplus')
  deallocate(chiG,stat=status)
  if (status/=0) call deallocate_error('gasteiger','chiG')
  deallocate(lgastfoc,stat=status)
  if (status/=0) call deallocate_error('gasteiger','lgastfoc')
  deallocate(ngastrptr,stat=status)
  if (status/=0) call deallocate_error('gasteiger','ngastrptr')
  deallocate(ngastptr,stat=status)
  if (status/=0) call deallocate_error('gasteiger','ngastptr')
!
!  Timing
!
  time2 = g_cpu_time()
  teem = teem + time2 - time1
#ifdef TRACE
  call trace_out('gasteiger')
#endif
!
  return
  end
