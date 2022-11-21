  subroutine shellopt
!
!  Subroutine for optimisation of shell positions only.
!
!   5/14 Created from optim
!   6/14 lfreezeok set to false
!   1/17 Call to transmat modified for distributed memory
!   2/17 nmin removed from arguments to minimise
!   3/17 hessian matrix changed to 2-D matrix to accommodate parallel case
!   3/17 Modifications made to allow for new variable order in iopt
!   7/17 Call to setoptptr added
!   2/18 Trace added
!   4/18 Allocation of hessian corrected to use nvaronnode for right-hand side
!   6/18 Parallel handling of nvar corrected
!   8/18 Modified due to changes in lstraincell algorithm
!   8/18 Adding 1 to strains 1-3 removed
!   3/19 x0 removed
!   3/19 iopt replaced by ioptindex and iopttype
!   3/19 Constraint handling now moved to subroutine
!   8/19 Arguments passed to harmonicrelax corrected
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
!  Julian Gale, CIC, Curtin University, August 2019
!
  use configurations
  use control
  use current
  use element,       only : maxele
  use energies,      only : fcsave
  use gulp_files
  use fitting
  use general
  use iochannels
  use optimisation
  use parallel
  use reallocate
  use symmetry
  use times
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  use xcgc

  implicit none
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: icount
  integer(i4)                                  :: icx
  integer(i4)                                  :: icy
  integer(i4)                                  :: icz
  integer(i4)                                  :: ierror
  integer(i4)                                  :: ifail
  integer(i4)                                  :: ind
  integer(i4)                                  :: indc
  integer(i4), dimension(:), allocatable       :: ioptindexsave
  integer(i4), dimension(:), allocatable       :: iopttypesave
  integer(i4)                                  :: j
  integer(i4)                                  :: ncellsave
  integer(i4)                                  :: ncellmaxsave
  integer(i4)                                  :: ncellminsave
  integer(i4)                                  :: node
  integer(i4), dimension(:), allocatable       :: node2varsave
  integer(i4)                                  :: nvarsave
  integer(i4)                                  :: nvaronnodesave
  integer(i4), dimension(:), allocatable       :: nvar2localsave
  integer(i4), dimension(:), allocatable       :: nvar2nodesave
  integer(i4),                            save :: ncflast = 0
  integer(i4),                            save :: nhwords = 1
  integer(i4),                            save :: nhuwords = 1
  integer(i4)                                  :: status
  logical,                                save :: lfirsttime = .true.
  logical,                                save :: lhess2D = .false.
  logical                                      :: lfreezeok
  logical                                      :: lgradloc
  logical                                      :: lgrad2loc
  logical                                      :: lharmloc
  logical                                      :: loptiloc
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: fc
  real(dp),    dimension(:,:), pointer,   save :: hess => null()
  real(dp)                                     :: pv
  real(dp)                                     :: time1
  real(dp)                                     :: xf
  real(dp)                                     :: yf
  real(dp)                                     :: zf
#ifdef TRACE
  call trace_in('shellopt')
#endif
!
!  Turn off printing during optimisation
!
  lopprt = .false. 
!
!  Save ncell since we don't want to optimise the cell here
!
  ncellsave = ncell
  ncellmaxsave = ncellmax
  ncellminsave = ncellmin
  ncell = 0
!
!  Save nvar and ioptindex/iopttype in case more than shell optimisation was specified in the input
!
  nvarsave = nvar
  nvaronnodesave = nvaronnode
  allocate(ioptindexsave(nvar),stat=status)
  if (status/=0) call outofmemory('shellopt','ioptindexsave')
  allocate(iopttypesave(nvar),stat=status)
  if (status/=0) call outofmemory('shellopt','iopttypesave')
  allocate(node2varsave(nvar),stat=status)
  if (status/=0) call outofmemory('shellopt','node2varsave')
  allocate(nvar2localsave(nvar),stat=status)
  if (status/=0) call outofmemory('shellopt','nvar2localsave')
  allocate(nvar2nodesave(nvar),stat=status)
  if (status/=0) call outofmemory('shellopt','nvar2nodesave')
!
  ioptindexsave(1:nvar) = ioptindex(1:nvar)
  iopttypesave(1:nvar) = iopttype(1:nvar)
  node2varsave(1:nvaronnode) = node2var(1:nvaronnode)
  nvar2localsave(1:nvar) = nvar2local(1:nvar)
  nvar2nodesave(1:nvar) = nvar2node(1:nvar)
!
!  Reset nvar and iopt based on shells only
!
  nvar = 0
  do i = 1,nvarsave
    ind = ioptindexsave(i)
!
!  Exclude strains
!
    if (iopttypesave(i).eq.iopt_xf.or.iopttypesave(i).eq.iopt_yf.or.iopttypesave(i).eq.iopt_zf) then
      if (iatn(ind).gt.maxele) then
        nvar = nvar + 1
        ioptindex(nvar) = ioptindexsave(i)
        iopttype(nvar) = iopttypesave(i)
      endif
    endif
  enddo
!
  if (nprocs.gt.1) then
    nvaronnode = 0
    node = 0
    icount = 0
!
!  If block size hasn't been input then choose a value based on the number of atoms versus processors
!
    if (nblocksizevar.eq.0) then
      nblocksizevar = 3
    endif
!
    do i = 1,nvar
      icount = icount + 1
      nvar2node(i) = node
      if (node.eq.procid) then
        nvaronnode = nvaronnode + 1
        node2var(nvaronnode) = i
        nvar2local(i) = nvaronnode
      else
        nvar2local(i) = 0
      endif
      if (icount.eq.nblocksizevar) then
        icount = 0
        node = node + 1
        if (node.eq.nprocs) node = 0
      endif
    enddo
  else
    nvaronnode = nvar
    do i = 1,nvar
      nvar2node(i) = 0
      node2var(i) = i
      nvar2local(i) = i
    enddo
  endif
!
!  Nullify Hessian pointer and initialise with basic size
!
  if (lfirsttime) then
    lfirsttime = .false.
    if (nprocs.gt.1) lhess2D = .true.
    if (nprocs.gt.1) then
      call realloc(hess,nhwords,nhuwords,ierror)
      if (ierror.ne.0) call outofmemory('shellopt','hess')
    else
      call realloc(hess,nhwords,nhuwords,ierror)
      if (ierror.ne.0) call outofmemory('shellopt','hess')
    endif
  endif
!
  lfreezeok = .false.
  loptiloc = lopt
  lgradloc = lgrad
  lharmloc = lharmrelax
  loptsuccess = .false.
!
!  Set flag to indicate whether second derivatives and therefore tmat will ever be needed
!
  lgrad2loc = (.not.lconj.and..not.llbfgs.and..not.lunit)
  if (lminch.and.mintype.le.2) lgrad2loc = .true.
!
  fc = 0.0_dp
  pv = 0.0_dp
!
!  Transfer all coordinates to configuration, including cores
!
  do i = 1,nasym
    call cart2frac(ndimen(ncf),xclat(i),yclat(i),zclat(i),rv,xf,yf,zf,icx,icy,icz)
    xafrac(i) = xf - icx
    yafrac(i) = yf - icy
    zafrac(i) = zf - icz
  enddo
!
!  Transfer breathing shell radii
!
  if (nbsmat.gt.0) then
    do i = 1,nasym
      rada(i) = radcfg(nsft+i)
    enddo
  endif
!***********************
!  Set freezing flags  *
!***********************
!
!  Initialise to fixed
!
  do i = 1,nasym
    lopf(i) = .false.
  enddo
  if (lbulknoopt) then
    loptiloc = .false.
  elseif (nvar.gt.0) then
!
!  Transfer configuration to variables
!
    call cfgtovar(nvar,xc)
!
    do i = 1,nvar
      ind = ioptindex(i)
      lopf(ind) = .true.
!
!  Check for constrained atoms
!
      if (ncon.gt.0) then
        do j = 1,ncon
          if (ioptindex(i).eq.ncvarind(j).and.iopttype(i).eq.ncvartyp(j)) then
            indc = ncfixind(j)
            if (ncfixtyp(j).eq.iopt_radius) then
              lopf(indc) = .true.
            elseif (ncfixtyp(j).eq.iopt_xf) then
              lopf(indc) = .true.
            elseif (ncfixtyp(j).eq.iopt_yf) then
              lopf(indc) = .true.
            elseif (ncfixtyp(j).eq.iopt_zf) then
              lopf(indc) = .true.
            endif
          endif
        enddo
      endif
    enddo
  elseif (nvar.eq.0.and.(lopt.or.lgrad.or.lharmrelax)) then
    nwarn = nwarn + 1
    loptiloc = .false.
    lgradloc = .false.
    lharmloc = .false.
  endif
  if ((loptiloc.or.lharmloc).and.(.not.lconj.or.lminch)) then
!*********************
!  Allocate hessian  *
!*********************
    if (lhess2D) then
!
!  2-D hessian
!
      if (llbfgs) then
        nhwords = nvar*(2*lmbfgsorder + 1) + 2*lmbfgsorder
!
!  Check if minimiser might change and need more memory
!
        if (lminch.and.mintype.le.4) then
          nhuwords = nvaronnode
        else
          nhuwords = 1_i4
        endif
      else
        nhwords = nvar
        nhuwords = nvaronnode
      endif
    else
!
!  1-D triangular half hessian
!
      if (llbfgs) then
        nhwords = nvar*(2*lmbfgsorder + 1) + 2*lmbfgsorder
!
!  Check if minimiser might change and need more memory
!
        if (lminch.and.mintype.le.4) then
          nhwords = max(nhwords,nvar*(nvar+1)/2)
        endif
      else
        nhwords = nvar*(nvar+1)/2
      endif
    endif
    call realloc(hess,nhwords,nhuwords,ierror)
    if (ierror.ne.0) call outofmemory('shellopt','hess')
  endif
  fcsave = fc
!********************************
!       Optimisation            *
!********************************
  if (loptiloc.or.lharmloc) then
    ifail = 1
    lfreeze = lfreezeok
!
!  Set pointers to atoms for optimisation after setting lfreeze
!
    call setoptptr(lfreeze)
!
    if (ioproc) call gflush(ioout)
!
!  Setup transformation matrix if needed
!
    if (ncf.ne.ncflast.and.lgrad2loc) then
      if (nprocs.gt.1) then
        call transmatd
      else
        call transmat
      endif
    endif
!
!  Start minimisation
!
    if (lharmloc) then
!
!  Implicit harmonic relaxation
!
      call harmonicrelax(xc,fc,gc,hess,nhwords,lhess2D,1_i4)
      ifail = 4
      loptsuccess = .true.
      if (ioproc) call outener
    else
!
!  Minimise static/free energy
!
      call minimise(xc,fc,gc,hess,nhwords,lhess2D,ifail,1_i4,.false.)
      loptsuccess = (ifail.eq.0)
    endif
!
!  Substitute parameters into place
!
    if (.not.lstraincell) then
      strain(1:nstrains) = 0.0_dp
    endif
!****************************************
!  Transfer variables to configuration  *
!****************************************
    call vartocfg(nvar,xc)
!**********************
!  Apply constraints  *
!**********************
    if (ncon.gt.0) then
      call applyconstraints
    endif
!****************************************
!  Return data to configuration arrays  *
!****************************************
!
!  Atomic positions
!
    if (ndim.ge.1.and.lmodco) then
      do i = 1,nasym
        xcfg(i+nsft) = mod(xafrac(i)+1.0_dp,1.0_dp)
      enddo
    else
      do i = 1,nasym
        xcfg(i+nsft) = xafrac(i)
      enddo
    endif
    if (ndim.ge.2.and.lmodco) then
      do i = 1,nasym
        ycfg(i+nsft) = mod(yafrac(i)+1.0_dp,1.0_dp)
      enddo
    else
      do i = 1,nasym
        ycfg(i+nsft) = yafrac(i)
      enddo
    endif
    if (ndim.eq.3.and.lmodco) then
      do i = 1,nasym
        zcfg(i+nsft) = mod(zafrac(i)+1.0_dp,1.0_dp)
      enddo
    else
      do i = 1,nasym
        zcfg(i+nsft) = zafrac(i)
      enddo
    endif
!
!  Radii
!
    if (nbsmat.gt.0) then
      do i = 1,nasym
        radcfg(i+nsft) = rada(i)
      enddo
    endif
!
!  Copy configuration coordinates back to current arrays
!
    do i = 1,nasym
      xafrac(i) = xcfg(nsft+i)
      yafrac(i) = ycfg(nsft+i)
      zafrac(i) = zcfg(nsft+i)
    enddo
    if (lsymopt) then
      call equpos(.true.,.false.)
    else
      do i = 1,numat
        xfrac(i) = xafrac(i)
        yfrac(i) = yafrac(i)
        zfrac(i) = zafrac(i)
      enddo
    endif
  endif
!************************
!  End of optimisation  *
!************************
  ncflast = ncf
!
!  Timing
!
  time1 = g_cpu_time()
  time1 = time1 - time0
!
!  Restore nvar and ioptindex/iopttype in case more than shell optimisation was specified in the input
!
  nvar = nvarsave
  nvaronnode = nvaronnodesave
!
  ioptindex(1:nvar) = ioptindexsave(1:nvar)
  iopttype(1:nvar) = iopttypesave(1:nvar)
  node2var(1:nvaronnode) = node2varsave(1:nvaronnode)
  nvar2local(1:nvar) = nvar2localsave(1:nvar)
  nvar2node(1:nvar) = nvar2nodesave(1:nvar)
!
  deallocate(nvar2nodesave,stat=status)
  if (status/=0) call deallocate_error('shellopt','nvar2nodesave')
  deallocate(nvar2localsave,stat=status)
  if (status/=0) call deallocate_error('shellopt','nvar2localsave')
  deallocate(node2varsave,stat=status)
  if (status/=0) call deallocate_error('shellopt','node2varsave')
  deallocate(iopttypesave,stat=status)
  if (status/=0) call deallocate_error('shellopt','iopttypesave')
  deallocate(ioptindexsave,stat=status)
  if (status/=0) call deallocate_error('shellopt','ioptindexsave')
!
!  Reset ncell 
!
  ncell = ncellsave
  ncellmax = ncellmaxsave
  ncellmin = ncellminsave
#ifdef TRACE
  call trace_out('shellopt')
#endif
!
  return
  end
