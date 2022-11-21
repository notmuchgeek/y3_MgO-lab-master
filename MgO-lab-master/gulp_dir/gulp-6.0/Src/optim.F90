  subroutine optim(lmain,lscan,nscanpoint)
!
!  Subroutine for optimisation of crystal structure
!  also functions in single point calculation mode
!  Second derivative option now included
!
!  lmain = is this a main call of optim?
!  lscan = is this one point of a scan?
!  nscanpoint = no. of scanpoint
!
!   6/95 Use of symmetrised second derivatives in property calc added
!   6/95 Modified for use of additive constraints
!   7/96 Correction added for effect of constraints on lopf
!  11/96 Change added to allow transmat to be called when lfreezeok
!        and number of variables is =< maxd2
!  12/96 Dynamical hessian allocation added
!   1/97 Dynamical hessian generalised for non-SG cases
!   2/98 Freezing turned off if symmetry is being used with EEM as
!        this algorithm combination doesn't work.
!   4/98 Freezing turned off altogether for EEM otherwise changes in
!        self energy get missed
!   6/98 Remove call to gaopt (SMW)
!   8/98 Analytical free energy derivative modifications added
!   2/99 Numerical gradients including hessian (SMW)
!   3/99 Option to minimise static energy before free energy added
!   5/99 Check on lfcphon to decide whether transmat should be called
!        added as fitting with phonons can overwrite tmat
!   8/99 Static call to funct added so that properties are correct
!        after a free energy minimisation
!  12/00 Modified to handle 1-D and 2-D cases
!   3/01 Freezing turned off for symmetry adapted optimisation since
!        it is hard to get the scaling of interactions correct for 
!        all cases.
!   8/01 lopf is now set even for non-freezing case since this can
!        still be used to save work
!  10/01 lfreeloc removed since free energy minimisation can be required
!        even when T = 0 due to ZPE
!  10/01 energycfg now scaled for full cell
!   2/02 Extra call to funct added to recalculate the second derivatives
!        when phonon/property is going to be called after an optimisation
!        otherwise acoustic branches can be wrong and xc contents updated
!   5/02 Array hess set to ensure that pointer is valid on calling functn
!   5/02 Order of calls to borncharge and property swapped to preserve
!        contents of dervi for defect calculations
!   6/02 Analytical calculation of second derivatives removed as precursor
!        to numerical first derivatives
!   9/02 Setting of lsymderv2 altered so that a single point calc can use
!        symmetry when no second derivatives are required
!  10/02 lopf initialised for all cases
!  11/02 scan input variables added to control outarc calls
!   6/03 XML modifications added
!   6/03 LM-BFGS modifications added
!   5/04 lmodco option introduced
!   9/04 Freezing turned off if variable charges are being used
!  10/04 Handling of finite difference cleaned up
!  11/04 Setting of inverse cell parameters added
!  12/04 Freezing turned off when bond order models are used 
!   3/06 Numerical free energy derivatives more widely enabled
!   5/06 Calls to outarc and outxyz suppressed during relaxed
!        fitting
!   5/06 Property output for clusters added
!  11/06 lfirst argument added to equpos call
!   3/07 Chemshell changes added
!   3/07 Flag for funct now passed as variable to avoid error
!  11/07 Modifications for reaxFF Coulomb term added
!  12/07 Unused variables removed
!   1/08 lreaxFFqreal removed
!   4/08 Misleading comment corrected about writing only shells
!   4/08 Chemshell output now requires that lgradloc be true
!   4/08 lphonloc now not disabled for finite difference case
!   4/08 ChemShell file handling for the lgradloc conditional
!   6/08 Initialisation of numnonsh corrected.
!   8/08 Format of chemshell write adjusted
!   9/08 Logic of output modified for ChemShell to handle loptiloc case
!   2/09 Hessian dimension passed to minimize for checking
!   2/09 Modified to accommodate new version of FoX and gulp_cml
!   3/09 Use of lfinitediff replaces testing of finite difference value
!   5/09 Option to generate second derivatives for properties by finite
!        differences added 
!   6/09 Error format updated to standard form
!   6/10 loptsuccess flag added to indicate a successful optimisation
!   7/11 Call to transmat suppressed if not needed as this can cause
!        a large memory overhead.
!  10/11 Correction to allocated memory for case where switch occurs
!        from lmbfgs to a Hessian based method. 
!   4/12 Comments relating to calling cml directly removed
!   5/12 Stress output moved here from property.
!   6/12 Dummy variables added to phonon and deffreq call
!   5/13 Freezing turned off for parallel second derivatives
!  10/13 Call to evaluate properties enabled for all dimensionalities
!        since it is now needed to compute Raman susceptibilities
!   3/14 Harmonic relaxation energy added
!   3/14 Calls to matinv renamed to matrix_inversion for benefit of ChemShell
!   3/14 Second derivatives added for shengBTE
!   3/14 Output channel for shengBTE set from module
!   6/14 Use of lnoexclude flag added
!  10/14 fc_supercell algorithm added
!   1/15 Call to outinertia changed to append frames
!   2/15 Modified to allow multiple frames to be written from
!        a translate scan.
!   2/15 Calls to outshengBTE moved into phonon routines
!   9/16 Flag now set during optimisation so the other routines
!        can adjust settings accordingly
!  11/16 Call to phonond added for the parallel case
!  12/16 Call added to phonond for parallel cluster frequencies 
!  12/16 Call changed to phonon for serial cluster frequencies 
!   1/17 Call to setoptptr added
!   1/17 Transmat call modified for distributed memory
!   2/17 nmin removed from arguments to minimise
!   2/17 hessian matrix changed to 2-D matrix to accommodate parallel case
!   3/17 Call to setoptptr moved to setup since it can be needed outside
!        of optimisation
!   4/17 Call to parallel phonon_fcd added
!   4/17 ChemShell interaction modified
!   4/17 lhess2D now access from module via an option
!   6/17 Old ChemShell restored as a compile option
!   7/17 Call to setoptptr moved here from setup since it needs lopf flags
!        to be set
!   8/17 Use of numerical Hessian for optimisation added
!  11/17 Reallocation of Hessian included for finite differences ahead of
!        call to functn
!  12/17 Calls to functnff added
!   1/18 Ghostcell algorithm added for phonons
!   1/18 Trace added
!   4/18 Allocation of hessian corrected to use nvaronnode for right-hand side
!   6/18 Strain cell option added
!   6/18 Handling of xyz files corrected for alterative scan mode
!   8/18 Changes made for lstraincell algorithm
!   8/18 Adding 1 to strains 1-3 removed
!  11/18 Turning off of lfinitestrain added when properties are being computed
!   1/19 Close for xyz changed
!   3/19 x0 removed
!   3/19 iopt replaced by ioptindex and iopttype
!   3/19 Constraint handling now moved to subroutine
!   4/19 Reseting of strains added for property calc after optimisation
!   4/19 Handling of saving of lfinitestrain for property calculation 
!        corrected to avoid potential access to uninitialised variable
!   5/19 Finite difference flag split for first and second derivatives
!   8/19 Arguments passed to harmonicrelax corrected
!  10/19 Freezing turned off for rigid molecule case
!  12/19 rvcfg cell is now explicitly saved during relaxed fitting and reset
!   2/20 Call to transmatmol added
!   3/20 Keyword added to output list of variables
!   4/20 Restart for rigid molecules added
!   7/20 Site energy call separated from optout
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
!  Scott Woodley, R.I.G.B., June 1997
!
  use bondorderdata, only : nbopot, nboQ
  use configurations
  use control
  use current
  use derivatives,   only : xdrv, ydrv, zdrv
  use derivatives,   only : lfcsupercell, lfinitestrain
  use element,       only : maxele
  use energies,      only : fcsave, fcstore
  use gulp_files
  use fitting
  use general
  use genetic,       only : lgacost
  use gulpchemsh
  use gulp_cml,      only : lcml, gulp_cml_StartModule, gulp_cml_EndModule, gulp_cml_structure
  use iochannels
  use maths,         only : lhess2D
  use molecule,      only : nmolptr, nmolatom, nmollist, molcom, molcomcfg, nmol, molQ, molQcfg
  use optimisation
  use parallel
  use progress,      only : lduring_opt
  use reallocate
  use scan,          only : ntran, ncscan
  use sutton,        only : lsuttonc
  use symmetry
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  use times
  use xcgc

  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                      :: nscanpoint
  logical,     intent(in)                      :: lmain
  logical,     intent(in)                      :: lscan
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ierror
  integer(i4)                                  :: ifail
  integer(i4)                                  :: iflag
  integer(i4)                                  :: ii
  integer(i4)                                  :: ind
  integer(i4)                                  :: indc
  integer(i4)                                  :: j
  integer(i4),                            save :: ncflast = 0
  integer(i4)                                  :: nfwords
  integer(i4)                                  :: nfuwords
  integer(i4),                            save :: nhwords = 1
  integer(i4),                            save :: nhuwords = 1
  integer(i4)                                  :: nm
  integer(i4)                                  :: ipunch      !ChemShell
  integer(i4)                                  :: numnonsh    !ChemShell
  integer(i4)                                  :: status
  logical                                      :: lappend
  logical,                                save :: lfirsttime = .true.
  logical                                      :: lfreezeok
  logical                                      :: lgradloc
  logical                                      :: lgrad2loc
  logical                                      :: lharmloc
  logical                                      :: loptiloc
  logical                                      :: lphonloc
  logical                                      :: lprintsave
  logical                                      :: lproploc
  logical                                      :: lsavefinitestrain
  logical                                      :: lstrold
  logical                                      :: lusenumericd2
  real(dp)                                     :: cost
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: fc
  real(dp),    dimension(:,:), pointer,   save :: hess => null()
  real(dp)                                     :: oldcell(6)
  real(dp)                                     :: pv
  real(dp)                                     :: rvcfgsave(3,3)
  real(dp)                                     :: strloc(6)
  real(dp)                                     :: time1
  real(dp)                                     :: face        !ChemShell
  real(dp)                                     :: facg        !ChemShell
  real(dp)                                     :: xtt         !ChemShell
  real(dp)                                     :: ytt         !ChemShell
  real(dp)                                     :: ztt         !ChemShell
  real(dp)                                     :: rvi(3,3)    !ChemShell
  real(dp)                                     :: rr(6)       !ChemShell
#ifdef TRACE
  call trace_in('optim')
#endif
!
!  Nullify Hessian pointer and initialise with basic size
!
  if (lfirsttime) then
    lfirsttime = .false.
    call realloc(hess,nhwords,nhuwords,ierror)
    if (ierror.ne.0) call outofmemory('optim','hess')
  endif
!
  lopprt = lmain
  lfreezeok = (.not.lnoexclude.and.ncell.eq.0.and..not.lfree.and. &
    .not.lsymopt.and..not.lsuttonc.and..not.lDoQDeriv1.and..not.lDoQDeriv2 &
    .and..not.lbrenner.and.(nbopot+nboQ.eq.0).and..not.lcosmo.and..not.lreaxFF &
    .and..not.lrigid)
  lgacost = (index(keyword,'cost').ne.0.and.lopt)
  loptiloc = lopt
  lgradloc = lgrad
  lharmloc = lharmrelax
  lphonloc = lphon
  lproploc = (lprop.or.lphon.or.lraman)
  loptsuccess = .false.
!
!  Set flag to indicate whether numerical second derivatives should be used for properties
!
  lusenumericd2 = (lproploc.and.lnoanald2.or.lnumerical)
!
!  Set flag to indicate whether second derivatives and therefore tmat will ever be needed
!
  lgrad2loc = (.not.lconj.and..not.llbfgs.and..not.lunit)
  if (lminch.and.mintype.le.2) lgrad2loc = .true.
!
!  If lgrad2loc and parallel then freezing can't be used
!
  if (lgrad2loc.and.nprocs.gt.1) lfreezeok = .false.
!
  fc = 0.0_dp
  pv = 0.0_dp
!
!  If relax mode then restore original cell coordinates
!
  if (lrelax) then
    do i = 1,nasym
      xstore(i) = xcfg(nsft+i)
      ystore(i) = ycfg(nsft+i)
      zstore(i) = zcfg(nsft+i)
    enddo
    if (nbsmat.gt.0) then
      do i = 1,nasym
        rstore(i) = radcfg(nsft+i)
      enddo
    endif
  elseif (lcomp) then
!
!  Store structure for comparison afterwards
!
    if (ndim.eq.3) then
      oldcell(1) = a
      oldcell(2) = b
      oldcell(3) = c
      oldcell(4) = alpha
      oldcell(5) = beta
      oldcell(6) = gamma
    elseif (ndim.eq.2) then
      oldcell(1) = a
      oldcell(2) = b
      oldcell(3) = alpha
    elseif (ndim.eq.1) then
      oldcell(1) = a
    endif
    do i = 1,nasym
      xstore(i) = xcfg(nsft+i)
      ystore(i) = ycfg(nsft+i)
      zstore(i) = zcfg(nsft+i)
    enddo
    if (nbsmat.gt.0) then
      do i = 1,nasym
        rstore(i) = radcfg(nsft+i)
      enddo
    endif
  endif
!
!  Copy configuration cell to saved copy
!
  if (ndim.gt.0) then
    rvcfgsave(1:3,1:ndim) = rvcfg(1:3,1:ndim,ncf)
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
    call cfgtovar(nvar,xc)
!
!  Set flag to indicate that atom has a variable
!
    do i = 1,nvar
      ind = ioptindex(i)
      if (iopttype(i).eq.iopt_xf) then
        lopf(ind) = .true.
      elseif (iopttype(i).eq.iopt_yf) then
        lopf(ind) = .true.
      elseif (iopttype(i).eq.iopt_zf) then
        lopf(ind) = .true.
      elseif (iopttype(i).eq.iopt_radius) then
        lopf(ind) = .true.
      elseif (iopttype(i).eq.iopt_xcom.or.iopttype(i).eq.iopt_ycom.or.iopttype(i).eq.iopt_zcom) then
        do nm = 1,nmolatom(ind)
          ii = nmollist(nmolptr(ind)+nm)
          lopf(ii) = .true.
        enddo
      elseif (iopttype(i).eq.iopt_xqtn.or.iopttype(i).eq.iopt_yqtn.or.iopttype(i).eq.iopt_zqtn) then
        do nm = 1,nmolatom(ind)
          ii = nmollist(nmolptr(ind)+nm)
          lopf(ii) = .true.
        enddo
      endif
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
    if (lopt) then
      call outwarning('No variables to optimise - single point performed',0_i4)
    else
      call outwarning('No variables for gradients - only energy calculated',0_i4)
    endif
    loptiloc = .false.
    lgradloc = .false.
    lharmloc = .false.
  endif
!
  if ((loptiloc.or.lharmloc.or.lrelax).and.(.not.lconj.or.lminch)) then
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
!
!  For rfo need storage space to avoid using disk
!
        if (lfree) then
          nfwords = 3*numat
          nhwords = max(nhwords,nfwords)
          nfuwords = 3*natomsonnode
          nhuwords = max(nhuwords,nfuwords)
        endif
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
!
!  For rfo need storage space to avoid using disk
!
        if (lfree) then
          nfwords = 3*numat*(3*numat+1)/2
          nhwords = max(nhwords,nfwords)
        endif
      endif
    endif
    call realloc(hess,nhwords,nhuwords,ierror)
    if (ierror.ne.0) call outofmemory('optim','hess')
  endif
!
!  Set flag to indicate that we are in optimisation mode
!
  lduring_opt = loptiloc
!
  if (lmain) then
    lfirst = .true.
    lfreeze = .false.
!
!  Set pointers to atoms for optimisation after setting lfreeze
!
    call setoptptr(lfreeze)
!
    iflag = 0
    if (.not.loptiloc) then
      if (lgradloc.or.lfit) then
        iflag = 1
      endif
      if (lborn.or.lproploc.or.lphon.or.lharmrelax) iflag = 2
      lstr = (ndim.gt.0)
      ifail = 4
    endif
    if (index(keyword,'sing').ne.0.or.nvar.eq.0) then
      ifail = 5
    endif
!
!  Set lfinitestrain before property calculation
!
    if (.not.lstraincellprop.and.iflag.ge.2) then
      lsavefinitestrain = lfinitestrain
      lfinitestrain = .false.
    endif
!********************************
!   Single point calculation    *
!********************************
    if (lfree) then
      lsymderv2 = .false.
      if (iflag.ge.2) then
        if (lhess2D) then
          nhwords = nvar
          nhuwords = nvaronnode
        else
          nhwords = nvar*(nvar+1)/2
        endif
        call realloc(hess,nhwords,nhuwords,ierror)
        if (ierror.ne.0) call outofmemory('optim','hess')
      endif
      if ((lfinitediff2.and.iflag.eq.2).and..not.ldefect) then
        call fefunctn(iflag,nvar,xc,fc,gc,hess,nhwords,lhess2D)
      elseif ((lfinitediff1.and.iflag.eq.1).and..not.ldefect) then
        call fefunctn(iflag,nvar,xc,fc,gc,hess,nhwords,lhess2D)
      else
        call fefunct(iflag,nvar,xc,fc,gc,hess,nhwords,lhess2D)
      endif
    else
!
!  If not doing optimisation turn off symmetrisation of second derivatives for property/phonon evaluation
!
      if (.not.loptiloc.and..not.lharmloc.and.iflag.ge.2) lsymderv2 = .false.
      if (lfinitediff2.and..not.ldefect) then
        if (iflag.ge.2) then
          if (lhess2D) then
            nhwords = nvar
            nhuwords = nvaronnode
          else
            nhwords = nvar*(nvar+1)/2
          endif
          call realloc(hess,nhwords,nhuwords,ierror)
          if (ierror.ne.0) call outofmemory('optim','hess')
        endif
        if (lnumerical.and.iflag.eq.2) then
          call functnf(iflag,nvar,xc,fc,gc,.true.)
        endif
        call functn(iflag,nvar,xc,fc,gc,hess,nhwords,lhess2D,1_i4)
      else
        if (lusenumericd2.and.iflag.eq.2) then
          if (lfastfd) then
            call functnff(iflag,nvar,xc,fc,gc,.true.)
          else
            call functnf(iflag,nvar,xc,fc,gc,.true.)
          endif
        else
          call funct(iflag,nvar,xc,fc,gc)
        endif
      endif
    endif
!
!  Reset lfinitestrain after property calculation
!
    if (.not.lstraincellprop.and.iflag.ge.2) then
      lfinitestrain = lsavefinitestrain
!
!  If second derivatives were calculated without finite strain due to lstraincellprop
!  but finite strain derivatives are needed, then re-compute first derivatives
!
      if (lfinitestrain.and.iflag.ge.1) then
        if (lfree) then
          if (lfinitediff1.and..not.ldefect) then
            call fefunctn(1_i4,nvar,xc,fc,gc,hess,nhwords,lhess2D)
          else
            call fefunct(1_i4,nvar,xc,fc,gc,hess,nhwords,lhess2D)
          endif
        else
          if (lfinitediff1.and..not.ldefect) then
            call functn(1_i4,nvar,xc,fc,gc,hess,nhwords,lhess2D,1_i4)
          else
            call funct(1_i4,nvar,xc,fc,gc)
          endif
        endif
      endif
    endif
!
    if (ioproc.and..not.lharmloc) call outener
    if (.not.lfree.and..not.loptiloc.and..not.lharmloc) then
      if ((lprop.or.lborn.or.lphon).and.(lewald.or.lwolf)) then
        call borncharge(.true.)
      endif
      if (lproploc) then
        if (ndim.eq.3) then
          call property3(.true.)
        elseif (ndim.eq.1.or.ndim.eq.2) then
          call property12(.true.)
        elseif (ndim.eq.0) then
          call property0(.true.)
        endif
      endif
      if (lphonloc) then
        if (ndim.gt.0) then
          if (lfcsupercell) then
            if (nprocs.gt.1) then
              call phonon_fcd(.true.,fc,0_i4,0_i4)
            else
              if (lghost) then
                call phonon_fcg(.true.,fc,0_i4,0_i4)
              else
                call phonon_fc(.true.,fc,0_i4,0_i4)
              endif
            endif
          else
            if (nprocs.gt.1) then
              call phonond(.true.,fc,0_i4,0_i4)
            else
              call phonon(.true.,fc,0_i4,0_i4)
            endif
          endif
        else
          if (nprocs.gt.1) then
            call phonond(.true.,fc,0_i4,0_i4)
          else
            call phonon(.true.,fc,0_i4,0_i4)
          endif
        endif
        if (llower) call setup(.true.)
      endif
      if (lproploc.and.ldefect.and.(lphonloc)) then
!
!  If a defect calculation is to be performed do silent property calculation
!  to restore matrices after phonon calculation
!
        lsymderv2 = .false.
        iflag = 2
        if (lusenumericd2) then
          if (lfastfd) then
            call functnff(iflag,nvar,xc,fc,gc,.true.)
          else
            call functnf(iflag,nvar,xc,fc,gc,.true.)
          endif
        else
          call funct(iflag,nvar,xc,fc,gc)
        endif
        call property3(.false.)
      endif
    endif
  endif
!**********************************************
!  Initialise arcfile for movie if requested  *
!**********************************************
  fcsave = fc
  if (ioproc.and..not.lrelax) then
    lappend = (lscan.and.nscanpoint.gt.0)
    if (larc.and.lmovie) then
      call outarc(16_i4,lappend,.false.)
    endif
    if (lxyz.and.lxyzmovie) then
      call outxyz(18_i4,lappend,.false.)
    endif
  endif
!********************************
!       Optimisation            *
!********************************
  if (loptiloc.or.lharmloc.or.lrelax) then
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
      if (lrigid) then
        if (ninternalmolT.gt.0) call transmatmolT
        if (ninternalmolQ.gt.0) call transmatmolQ
      endif
    endif
!
!  Output optimisation parameters, if main call
!
    if (lmain.and..not.lharmloc) then
      if (lPrintVar) call outoptvar(xc)
      call optin(lhess2D)
    endif
!
!  Start minimisation
!
    if (lstaticfirst.and.lfree) then
!
!  Minimise static energy prior to free energy minimisation
!
      lprintsave = lopprt
      lopprt = .false.
      call minimise(xc,fc,gc,hess,nhwords,lhess2D,ifail,1_i4,.false.)
      lopprt = lprintsave
      loptsuccess = (ifail.eq.0)
!
!  Reset number of frequencies to avoid false warning messages
!  at start of free energy minimisation
!
      nummode = 0
    endif
    if (lharmloc) then
!
!  Implicit harmonic relaxation
!
      call harmonicrelax(xc,fc,gc,hess(1,1),nhwords,lhess2D,1_i4)
      ifail = 4
      loptsuccess = .true.
      if (ioproc) call outener
    else
!
!  Minimise static/free energy
!
      call minimise(xc,fc,gc,hess,nhwords,lhess2D,ifail,1_i4,lfree)
      loptsuccess = (ifail.eq.0)
    endif
!
!  Set flag to indicate that we have finished with optimisation mode
!
    lduring_opt = .false.
!
!  Calculate the elastic constants at end of optimisation
!  Need to switch off lopt to enable second deriv calcn.
!
    lfreeze = .false.
!
!  Set pointers to atoms for optimisation after setting lfreeze
!
    call setoptptr(lfreeze)
!
    if ((lborn.or.lproploc.or.lphonloc).and.(loptiloc.or.lharmloc)) then
      lopt = .false.
      lstrold = lstr
      lstr = (ndim.eq.3)
      iflag = 2
      lsymderv2 = .false.
      if (lfree) then
        if (lfinitediff2.and..not.ldefect) then
          call fefunctn(iflag,nvar,xc,fc,gc,hess,nhwords,lhess2D)
        else
          call fefunct(iflag,nvar,xc,fc,gc,hess,nhwords,lhess2D)
        endif
      else
        if (lusenumericd2) then
          if (lfastfd) then
            call functnff(iflag,nvar,xc,fc,gc,.true.)
          else
            call functnf(iflag,nvar,xc,fc,gc,.true.)
          endif
        else
          call funct(iflag,nvar,xc,fc,gc)
        endif
      endif
      lopt = .true.
      lstr = lstrold
    elseif (lfreezeok) then
      iflag = 1
      if (lfree) then
        if (lfinitediff1.and..not.ldefect) then
          call fefunctn(iflag,nvar,xc,fc,gc,hess,nhwords,lhess2D)
        else
          call fefunct(iflag,nvar,xc,fc,gc,hess,nhwords,lhess2D)
        endif
      else
        call funct(iflag,nvar,xc,fc,gc)
      endif
    endif
!
!  Substitute parameters into place
!
    if (.not.loptcellpar.and..not.lstraincell) then
      do i = 1,nstrains
        strain(i) = 0.0_dp
      enddo
    endif
!**************************************
!  Transfer xc back to configuration  *
!**************************************
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
    if (ncell.gt.0) then
      if (.not.loptcellpar) then
!
!  Strains
!
        do i = 1,3
          rv(1,i) = rvcfg(1,i,ncf)
          rv(2,i) = rvcfg(2,i,ncf)
          rv(3,i) = rvcfg(3,i,ncf)
        enddo
!
!  Apply strain due to optimisation
!
        if (ndim.eq.3) then
          call strain3D(strain,rv)
          call uncell3D(rv,a,b,c,alpha,beta,gamma)
          if (a.gt.1.0d-12) then
            recipa = 1.0_dp/a
          else
            recipa = 0.0_dp
          endif
          if (b.gt.1.0d-12) then
            recipb = 1.0_dp/b
          else
            recipb = 0.0_dp
          endif
          if (c.gt.1.0d-12) then
            recipc = 1.0_dp/c
          else
            recipc = 0.0_dp
          endif
          if (.not.lrelax) then
            if (lstraincell) then
              call getstrain3D(rv,strloc)
              straincfg(1:6,ncf) = strloc(1:6)
            else
              do i = 1,3
                rvcfg(1,i,ncf) = rv(1,i)
                rvcfg(2,i,ncf) = rv(2,i)
                rvcfg(3,i,ncf) = rv(3,i)
              enddo
            endif
          endif
        elseif (ndim.eq.2) then
          call strain2D(strain,rv)
          call uncell2D(rv,a,b,alpha)
          if (a.gt.1.0d-12) then
            recipa = 1.0_dp/a
          else
            recipa = 0.0_dp
          endif
          if (b.gt.1.0d-12) then
            recipb = 1.0_dp/b
          else
            recipb = 0.0_dp
          endif
          if (.not.lrelax) then
            if (lstraincell) then
              call getstrain2D(rv,strloc)
              straincfg(1:3,ncf) = strloc(1:3)
            else
              do i = 1,2
                rvcfg(1,i,ncf) = rv(1,i)
                rvcfg(2,i,ncf) = rv(2,i)
                rvcfg(3,i,ncf) = rv(3,i)
              enddo
            endif
          endif
        elseif (ndim.eq.1) then
          call strain1D(strain,rv)
          call uncell1D(rv,a)
          if (a.gt.1.0d-12) then
            recipa = 1.0_dp/a
          else
            recipa = 0.0_dp
          endif
          if (.not.lrelax) then
            if (lstraincell) then
              call getstrain1D(rv,strloc)
              straincfg(1,ncf) = strloc(1)
            else
              rvcfg(1,1,ncf) = rv(1,1)
              rvcfg(2,1,ncf) = rv(2,1)
              rvcfg(3,1,ncf) = rv(3,1)
            endif
          endif
        endif
      endif
    endif
    if (lrelax) then
!
!  Reset configuration cell for safety during relax fitting
!
      if (ndim.gt.0) then
        rvcfg(1:3,1:ndim,ncf) = rvcfgsave(1:3,1:ndim)
      endif
    endif
!
!  Atomic positions
!
    if (.not.lrelax) then
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
!  Rigid molecule centre of mass
!
      if (lrigid) then
!
!  Convert COM to Cartesian for configuration store and copy quaternions back to configuration array
!
        do i = 1,nmol
          molQcfg(1:3,i,ncf) = molQ(1:3,i)
          if (ndim.eq.3) then
            molcomcfg(1,i,ncf) = molcom(1,i)*rvcfg(1,1,ncf) + molcom(2,i)*rvcfg(1,2,ncf) + molcom(3,i)*rvcfg(1,3,ncf)
            molcomcfg(2,i,ncf) = molcom(1,i)*rvcfg(2,1,ncf) + molcom(2,i)*rvcfg(2,2,ncf) + molcom(3,i)*rvcfg(2,3,ncf)
            molcomcfg(3,i,ncf) = molcom(1,i)*rvcfg(3,1,ncf) + molcom(2,i)*rvcfg(3,2,ncf) + molcom(3,i)*rvcfg(3,3,ncf)
          elseif (ndim.eq.2) then
            molcomcfg(1,i,ncf) = molcom(1,i)*rvcfg(1,1,ncf) + molcom(2,i)*rvcfg(1,2,ncf)
            molcomcfg(2,i,ncf) = molcom(1,i)*rvcfg(2,1,ncf) + molcom(2,i)*rvcfg(2,2,ncf)
            molcomcfg(3,i,ncf) = molcom(3,i)
          elseif (ndim.eq.1) then
            molcomcfg(1,i,ncf) = molcom(1,i)*rvcfg(1,1,ncf)
            molcomcfg(2,i,ncf) = molcom(2,i)
            molcomcfg(3,i,ncf) = molcom(3,i)
          else
            molcomcfg(1,i,ncf) = molcom(1,i)
            molcomcfg(2,i,ncf) = molcom(2,i)
            molcomcfg(3,i,ncf) = molcom(3,i)
          endif
        enddo
      endif
    endif
  endif
!**********************************
!  Output arc file or stop movie  *
!**********************************
  if (ioproc.and..not.lrelax) then
    if (larc) then
      if (lmovie) then
        if (ntran(ncf).gt.0) then
          if (nscanpoint.eq.ntran(ncf)) close(16)
        elseif (ncscan(ncf).gt.0) then
          if (nscanpoint.eq.ncscan(ncf)) close(16)
        endif
      else
        call outarc(16_i4,.false.,.false.)
        close(16)
      endif
    endif
    if (lxyz) then
      if (lxyzmovie) then
        if (ntran(ncf).gt.0) then
          if (nscanpoint.eq.ntran(ncf)) close(18)
        elseif (ncscan(ncf).gt.0) then
          if (nscanpoint.eq.ncscan(ncf)) close(18)
        else
          close(18)
        endif
      else
        lappend = (lscan.and.nscanpoint.gt.0)
        call outxyz(18_i4,.false.,.false.)
        close(18)
      endif
    endif
  endif
!*******************************************
!  Output moment of inertia for molecules  *
!*******************************************
  if (ioproc.and.linertia) then
    lappend = (lscan.and.nscanpoint.gt.0)
    call outinertia(19_i4,lappend)
  endif
!************************
!  End of optimisation  *
!************************
  ncflast = ncf
  if (lmain) then
! 
! New CML - open finalization module and dump final structure
!
    if (lcml) then
      call gulp_cml_StartModule(title='Finalization') 
      call gulp_cml_structure(ncf)
    endif 
!
    energycfg(ncf) = fc*dble(icentfct(ncbl))
    if (loptiloc.or.lharmloc.or.lgradloc) then
      call optout(ifail,ncf,fc,gc,oldcell)
      if (lstressout) call outstresses(.true.)
      if (latomicstress) call outatomicstress
    endif
    if (lsiteenergy) call outsiteenergy(fc)
!
!  Set lfinitestrain before property calculation
!
    if (.not.lstraincellprop) then
      lsavefinitestrain = lfinitestrain
      lfinitestrain = .false.
    endif
!
    if ((loptiloc.or.lharmloc.or.lfree).and.(lborn.or.lproploc.or.lphonloc)) then
!
!  If this was a free energy minimisation then perform a static calculation in order 
!  to get the correct properties at the end. It is also important to reset the values 
!  of xc so that the correct structure is used for the properties.
!
      call cfgtovar(nvar,xc)
!
!  Reset strains to zero for property calculation
!
      do i = 1,nvar
        if (iopttype(i).eq.iopt_strain.and..not.loptcellpar.and..not.lstraincell) then
          xc(i) = 0.0_dp
        endif
      enddo
!
      iflag = 2
      if (lusenumericd2) then
        if (lfastfd) then
          call functnff(iflag,nvar,xc,fc,gc,.true.)
        else
          call functnf(iflag,nvar,xc,fc,gc,.true.)
        endif
      else
        call funct(iflag,nvar,xc,fc,gc)
      endif
    endif
    if (loptiloc.or.lharmloc.or.lfree) then
      if ((lprop.or.lborn.or.lphon).and.(lewald.or.lwolf)) then
        call borncharge(.true.)
      endif
      if (lproploc) then
        if (ndim.eq.3) then
          call property3(.true.)
        elseif (ndim.eq.1.or.ndim.eq.2) then
          call property12(.true.)
        elseif (ndim.eq.0) then
          call property0(.true.)
        endif
      endif
      if (lphonloc) then
        if (ndim.gt.0) then
          if (lfcsupercell) then
            if (nprocs.gt.1) then
              call phonon_fcd(.true.,fc,0_i4,0_i4)
            else
              if (lghost) then
                call phonon_fcg(.true.,fc,0_i4,0_i4)
              else
                call phonon_fc(.true.,fc,0_i4,0_i4)
              endif
            endif
          else
            if (nprocs.gt.1) then
              call phonond(.true.,fc,0_i4,0_i4)
            else
              call phonon(.true.,fc,0_i4,0_i4)
            endif
          endif
        else
          if (nprocs.gt.1) then
            call phonond(.true.,fc,0_i4,0_i4)
          else
            call phonon(.true.,fc,0_i4,0_i4)
          endif
        endif
        if (llower) call setup(.true.)
      endif
      if (lproploc.and.ldefect.and.lphonloc) then
!
!  If a defect calculation is to be performed do silent property calculation
!  to restore matrices after phonon calculation
!
        iflag = 2
        if (lusenumericd2) then
          if (lfastfd) then
            call functnff(iflag,nvar,xc,fc,gc,.true.)
          else
            call functnf(iflag,nvar,xc,fc,gc,.true.)
          endif
        else
          call funct(iflag,nvar,xc,fc,gc)
        endif
        call property3(.false.)
      endif
    endif
    if (lcml) call gulp_cml_EndModule
!
!  Reset lfinitestrain after property calculation
!
    if (.not.lstraincellprop) then
      lfinitestrain = lsavefinitestrain
    endif
  endif
!
!  Print out value of cost function if requested keyword 'cost'
!
  if (ioproc.and.index(keyword,'cost').ne.0.and..not.lopt) then
    if (.not.lga.and..not.lpredict) then
      lgacost = .true.
      iflag = 0
      call funct(iflag,nvar,xc,cost,gc)
      write(ioout,'(/,''  The value of the cost function for this initial structure is '',f14.7,/)')cost
      write(ioout,'(/,''--------------------------------------------------------------------------------''/)')
      lgacost = .false.
    endif
  endif
!
!  Timing
!
  time1 = g_cpu_time()
  time1 = time1 - time0
  if (ioproc) then
    if (lmain.and.loptiloc) then
      write(ioout,'(/,''  Time to end of optimisation = '',f12.4,'' seconds'',/)')time1
    elseif (lmain.and.lharmloc) then
      write(ioout,'(/,''  Time to end of harmonic relaxtion = '',f12.4,'' seconds'',/)')time1
    elseif (lmain.and.(lproploc.or.lphonloc)) then
      write(ioout,'(/,''  Time to end of properties = '',f12.4,'' seconds'',/)')time1
    elseif (lmain.and.lgradloc) then
      write(ioout,'(/,''  Time to end of gradients = '',f12.4,'' seconds'',/)')time1
    endif
  endif
!********************
!  CHEMSHELL_START  *
!********************
#ifdef OLDCS
  if (ioproc .and. ichemsh_qm.ge.0) then
#else
  if (ioproc .and. ichemsh_output.eq.2) then
#endif

    ipunch = 7

    open(ipunch,file='gulp.energy',form='formatted')

    face = 0.03674896_dp
    facg = 0.0194467_dp

    write(ipunch,*) "block=matrix records=0"
    write(ipunch,*) "block=matrix_title records=1"
    write(ipunch,*) "Gulp energy"
    write(ipunch,*) "block=dense_real_matrix records=1 dimensions=1 1"

    write(ipunch,100) fcstore*face

    close(unit=ipunch)

    numnonsh = 0
    do i = 1,numat
      if (iatn(i).le.maxele) then
        numnonsh = numnonsh + 1
      endif
    enddo

    if (lgradloc.or.loptiloc.or.lharmloc) then
      open(ipunch,file='gulp.gradient',form='formatted')

      write(ipunch,*) "block=matrix records=0"
      write(ipunch,*) "block=matrix_title records=1"
      write(ipunch,*) "Gulp gradient"

      write(ipunch,101) 3*numnonsh, 3, numnonsh
!
!  Check that derivatives have been generated; if not, set them to
!  zeroes as the old code would have (GULP 1.3)
!
      if (.not. (associated(xdrv) .and. associated(ydrv) .and. associated(zdrv))) then
        call realloc(xdrv,numat,status)
        if (status.ne.0) call outofmemory('optim','xdrv')
        call realloc(ydrv,numat,status)
        if (status.ne.0) call outofmemory('optim','ydrv')
        call realloc(zdrv,numat,status)
        if (status.ne.0) call outofmemory('optim','zdrv')
        xdrv = 0.0_dp
        ydrv = 0.0_dp
        zdrv = 0.0_dp
      endif
!
!  Write only atoms, not shells
!
      if (ndimen(ncf).eq.3) then
        if (ncbl.gt.1) then
          call outerror('No centred cells',0_i4)
          call stopnow('optim')
        endif
      
        do i = 1,3
          rvi(1,i) = rv(1,i)
          rvi(2,i) = rv(2,i)
          rvi(3,i) = rv(3,i)
        enddo
        ifail = 0
        call matrix_inversion(rvi,3_i4,3_i4,rr,ifail)
      
        do i = 1,numat
          if (nat(i).le.maxele) then
            xtt = xdrv(i)*rvi(1,1) + ydrv(i)*rvi(2,1) + zdrv(i)*rvi(3,1)
            ytt = xdrv(i)*rvi(1,2) + ydrv(i)*rvi(2,2) + zdrv(i)*rvi(3,2)
            ztt = xdrv(i)*rvi(1,3) + ydrv(i)*rvi(2,3) + zdrv(i)*rvi(3,3)
            write(ipunch,100) xtt*facg
            write(ipunch,100) ytt*facg
            write(ipunch,100) ztt*facg
          endif
        enddo
      else
        do i = 1,numat
          if (nat(i).le.maxele) then
            write(ipunch,100) xdrv(i)*facg
            write(ipunch,100) ydrv(i)*facg
            write(ipunch,100) zdrv(i)*facg
          endif
        enddo
      endif

      close(unit=ipunch)

    endif

    open(ipunch,file='gulp.relshel',form='formatted')

    write(ipunch,*) "block=matrix records=0"
    write(ipunch,*) "block=matrix_title records=1"
    write(ipunch,*) "Gulp relaxed shell coordinates"

    write(ipunch,101) 3*(nasym-numnonsh), 3, nasym-numnonsh

    do i = 1,numat
      if (iatn(i).gt.maxele) then
        write(ipunch,100) xclat(i)/0.52917706_dp
        write(ipunch,100) yclat(i)/0.52917706_dp
        write(ipunch,100) zclat(i)/0.52917706_dp
      endif
    enddo

    close(unit=ipunch)

  endif ! ChemShell output

100  format(f28.14)
101  format('block=dense_real_matrix records=',i6,' dimensions=',2i6)
!******************
!  CHEMSHELL_END  *
!******************
!
#ifdef TRACE
  call trace_out('optim')
#endif
  return
  end
