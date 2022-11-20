  subroutine defopt(lmain)
!
!  Main subroutine that controls the optimisation of defects
!
!   6/95 Modified to allow for additive defect constraints
!  11/95 Defect flags modified for correct restart
!   3/97 Kindex now always called for safety
!   7/97 Initialisation of region2 densities for EAM added
!   8/02 Brenner potential added
!  12/03 Finite difference option added for defect derivatives
!  10/06 Region 1 - 2 handling modified to avoid subtraction
!  11/07 Unused variables removed
!  12/08 scrho, and defect versions, converted to 2-D array for MEAM
!   2/09 Hessian dimension passed to minimize for checking
!   3/09 Use of lfinitediff replaces testing of finite difference value
!   6/09 Module name changed from three to m_three
!  11/09 Hessian allocation for lmbfgs case corrected
!   6/12 Dummy variables added to deffreq call.
!   3/13 Use of scratch files for defect information removed
!   3/13 xdis, ydis, zdis no longer written to channel 48 since array
!        shouldn't be modified prior to use in move2a1
!   3/14 Pointer nullified on declaration line
!   8/14 eatom added to call to density routines
!  12/16 imode argument removed from deffreq calls
!   2/17 nmin removed from arguments to minimise
!   3/17 2-D hessian option added for compatibility with bulk calculation
!        and parallelisation
!   5/17 lhess2D now access from module via an option
!   5/17 Argument added to setvarnoded2 call
!   5/17 lnewdefalg now accessed from control module
!   5/17 lnewdefalg parallel calls added
!   5/17 Handling of lminch with 2-D hessian added
!   5/17 Calls to deffreqd added
!   6/17 four12d call added
!   1/18 Trace added
!   4/18 Allocation of hessian corrected to use nvaronnode for right-hand side
!   3/19 x0 removed
!   3/19 Change of idopt to idoptindex and idopttype
!   3/19 Change of defect constraints to have index and type
!   5/19 Finite difference flag split for first and second derivatives
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
  use bondorderdata, only : nbopot
  use control
  use current
  use defects
  use derivatives
  use eam,           only : lMEAMden, maxmeamcomponent
  use energies
  use gulp_files
  use four
  use general,       only : time0, nwarn, lfinitediff2
  use iochannels
  use maths,         only : lhess2D
  use optimisation
  use parallel
  use reallocate
  use region2a
  use sutton
  use symmetry
  use m_three
  use times
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  logical                                      :: lmain
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ierror
  integer(i4)                                  :: ifail
  integer(i4)                                  :: iflag
  integer(i4)                                  :: ii
  integer(i4)                                  :: ip
  integer(i4), dimension(:), allocatable       :: itmp
  integer(i4)                                  :: ix
  integer(i4)                                  :: iy
  integer(i4)                                  :: iz
  integer(i4)                                  :: ndind
  integer(i4),                            save :: ndwords
  integer(i4),                            save :: nduwords
  integer(i4)                                  :: nf
  integer(i4)                                  :: nlreg2
  integer(i4)                                  :: npsi
  integer(i4)                                  :: nv
  integer(i4)                                  :: status
  logical,                                save :: lfirsttime = .true.
  logical                                      :: lgradloc
  logical                                      :: loptiloc
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: csft(3)
  real(dp)                                     :: fc
  real(dp),    dimension(:), allocatable       :: gc
  real(dp),    dimension(:,:), pointer,   save :: hess => null()
  real(dp)                                     :: time1
  real(dp)                                     :: x0nf
  real(dp),    dimension(:), allocatable       :: xc
  real(dp),    dimension(:), allocatable       :: xdstore
  real(dp),    dimension(:), allocatable       :: ydstore
  real(dp),    dimension(:), allocatable       :: zdstore
#ifdef TRACE
  call trace_in('defopt')
#endif
!
!  Nullify Hessian pointer and initialise with basic size 
!
  if (lfirsttime) then
    lfirsttime = .false.
    ndwords = 1
    nduwords = 1
    call realloc(hess,ndwords,nduwords,ierror)
    if (ierror.ne.0) call outofmemory('defopt','hess')
  endif
!
  lopprt = lmain
  loptiloc = lopt
  lgradloc = lgrad
!
!  If ndasym  =  nreg1 turn off symmetry in calculation method
!
  ld1sym = (ndasym.le.nreg1.and.ldsym)
!
!  Set second derivative symmetry flag - only use symmetry
!  adapted algorithm if ndasym < 1/2 nreg1, otherwise the
!  memory demands are larger.
!
!  Algorithm based on symmetry is actually faster because
!  of nature of full multiply
!      ld2sym = (index(keyword,'nod2').eq.0.and.ld1sym.and.ndasym.le.(nreg1/2))
!
  ld2sym = (index(keyword,'nod2').eq.0.and.ld1sym)
!
!  Allocate local memory
!
  allocate(gc(nvar),stat=status)
  if (status/=0) call outofmemory('defopt','gc')
  allocate(xc(nvar),stat=status)
  if (status/=0) call outofmemory('defopt','xc')
!
!  Allocate module memory
!
  if (ldsym) then
    allocate(xdstore(ndasym),stat=status)
    if (status/=0) call outofmemory('defopt','xdstore')
    allocate(ydstore(ndasym),stat=status)
    if (status/=0) call outofmemory('defopt','ydstore')
    allocate(zdstore(ndasym),stat=status)
    if (status/=0) call outofmemory('defopt','zdstore')
  else
    allocate(xdstore(nreg1),stat=status)
    if (status/=0) call outofmemory('defopt','xdstore')
    allocate(ydstore(nreg1),stat=status)
    if (status/=0) call outofmemory('defopt','ydstore')
    allocate(zdstore(nreg1),stat=status)
    if (status/=0) call outofmemory('defopt','zdstore')
  endif
!
!  Set parallel distribution of variables
!
  call setvarnoded2(.true.)
!
  if (lcomp) then
    if (ldsym) then
      do i = 1,ndasym
        ii = ndsptr(i)
        xdstore(i) = xdefe(ii)
        ydstore(i) = ydefe(ii)
        zdstore(i) = zdefe(ii)
      enddo
    else
      do i = 1,nreg1
        xdstore(i) = xdefe(i)
        ydstore(i) = ydefe(i)
        zdstore(i) = zdefe(i)
      enddo
    endif
  endif
!
!  Initialise region 2 displacements
!
  nlreg2 = max(ntreg2,npreg2)
  do i = 1,nlreg2
    xdis(i) = 0.0_dp
    ydis(i) = 0.0_dp
    zdis(i) = 0.0_dp
  enddo
  if (nvar.gt.0) then
    call defcfgtovar(nvar,xc)
  elseif (nvar.eq.0.and.(lopt.or.lgrad)) then
    nwarn = nwarn + 1
    call outwarning('No variables specified for optimisation - single point performed',0_i4)
    loptiloc = .false.
    lgradloc = .false.
  endif
  if (lmain) then
    leprint = .true.
    lfirst = .true.
    iflag = 0
    if (.not.loptiloc) then
      if (lgradloc.or.lfit) then
        iflag = 1
      endif
      if (lfreq) iflag = 2
      ifail = 4
    endif
    if (index(keyword,'sing').ne.0.or.nvar.eq.0) then
      ifail = 5
    endif
  endif
!**********************************************************
!  Subtract perfect region 1 contribution - one off task  *
!**********************************************************
!
!  If a bulk phonon calculation has been performed then
!  we need to reset the K vectors before starting on
!  defect calculation
!
!  It turns out that for numerical stability it is necessary
!  to call kindex at the start of all types of defect runs
!  as just slight differences between the K vectors at the
!  end of an optimisation and here can cause problems.
!
  if (lewald) call kindex
  if (lsuttonc) then
!
!  Initialise region 2 densities to bulk values so we can
!  subtract off density from perfect region 1 here.
!
    if (lMEAMden) then
      do i = 1,npreg2
        npsi = nps(i)
        dscrhor2p(1:maxmeamcomponent,i) = scrho(1:maxmeamcomponent,nrelf2a(npsi))
      enddo
    else
      do i = 1,npreg2
        npsi = nps(i)
        dscrhor2p(1,i) = scrho(1,nrelf2a(npsi))
      enddo
    endif
  endif
!
!  Twobody potentials
!
  eatom = 0.0_dp
  ereal = 0.0_dp
  ec6 = 0.0_dp
  erecip = 0.0_dp
  if (lewald) then
    call recip12a(erecip,ec6,.false.,.false.,2_i4)
  endif
  if (lnewdefalg) then
    if (nprocs.gt.1) then
      call real12a2d(eatom,ereal,erecip,ec6,.false.,.false.,2_i4)
    else
      call real12a2(eatom,ereal,erecip,ec6,.false.,.false.,2_i4)
    endif
    if (lsuttonc.and.lMEAMden) call density12a2(2_i4,eatom)
  else
    call real12a(eatom,ereal,erecip,ec6,.false.,.false.,2_i4)
  endif
  e12aold = erecip + eatom + ereal + ec6
  ereal = 0.0_dp
  eatom = 0.0_dp
  if (lnewdefalg) then
    if (nprocs.gt.1) then
      call real112d(eatom,ereal,.false.,.false.,4_i4)
    else
      call real112(eatom,ereal,.false.,.false.,4_i4)
    endif
  endif
  e12aold = e12aold - eatom - ereal
  ereal = 0.0_dp
  eatom = 0.0_dp
  if (lnewdefalg) then
    if (nprocs.gt.1) then
      call real112d(eatom,ereal,.false.,.false.,2_i4)
    else
      call real112(eatom,ereal,.false.,.false.,2_i4)
    endif
  else
    call real11(eatom,ereal,.false.,.false.,2_i4)
  endif
  e11old = eatom + ereal
!
!  Many body potentials
!
  e12told = 0.0_dp
  e12fold = 0.0_dp
  e12mold = 0.0_dp
  e12bold = 0.0_dp
  e12boold = 0.0_dp
  if (nthb.gt.0) then
    if (nprocs.gt.1) then
      call three12d(e12told,.false.,.false.,2_i4)
    else
      call three12(e12told,.false.,.false.,2_i4)
    endif
  endif
  if (nfor.gt.0) then
    if (nprocs.gt.1) then
      call four12d(e12fold,.false.,.false.,2_i4)
    else
      call four12(e12fold,.false.,.false.,2_i4)
    endif
  endif
  if (lsuttonc) then
    if (nprocs.gt.1) then
      call many12d(e12mold,.false.,.false.,2_i4,.false.)
    else
      call many12(e12mold,.false.,.false.,2_i4,.false.)
    endif
  endif
  if (lbrenner) then
    if (nprocs.gt.1) then
      call brenner12d(e12bold,2_i4,.false.,.false.)
    else
      call brenner12(e12bold,2_i4,.false.,.false.)
    endif
  endif
  if (nbopot.gt.0) then
    if (nprocs.gt.1) then
      call bondorder12d(e12boold,2_i4,.false.,.false.)
    else
      call bondorder12(e12boold,2_i4,.false.,.false.)
    endif
  endif
!********************************
!   Single point calculation    *
!********************************
  if (lmain) then
    if (.not.loptiloc) ld2sym = .false.
    if (lfinitediff2) then
!
!  Allocate hessian
!
      if (lhess2D) then
        if (iflag.ge.2) then
          ndwords = nvar
          nduwords = nvaronnode
          call realloc(hess,ndwords,nduwords,ierror)
          if (ierror.ne.0) call outofmemory('defopt','hess')
        endif
      else
        if (iflag.ge.2) then
          ndwords = nvar*(nvar+1)/2
          call realloc(hess,ndwords,nduwords,ierror)
          if (ierror.ne.0) call outofmemory('defopt','hess')
        endif
      endif
      call functn(iflag,nvar,xc,fc,gc,hess,ndwords,lhess2D,2_i4)
    else
      call deffun(iflag,nvar,xc,fc,gc)
    endif
    leprint = .false.
    if (lfreq.and..not.loptiloc) then
      if (nprocs.gt.1) then
        call deffreqd(.true.,fc)
      else
        call deffreq(.true.,fc)
      endif
    endif
  endif
!*******************************************
!  Initialise arc file if movie requested  *
!*******************************************
  if (ioproc) then
    if (larc.and.lmovie) then
      call outarc(16_i4,.false.,.true.)
    endif
    if (lxyz.and.lxyzmovie) then
      call outxyz(18_i4,.false.,.true.)
    endif
  endif
!********************************
!       Optimisation            *
!********************************
  if (loptiloc) then
    ifail = 1
    if (ioproc) call gflush(ioout)
    if (ld2sym) then
      if (ndasym.eq.nreg1) ld2sym = .false.
      if (ldbsm) ld2sym = .false.
    endif
!
!  Set up second derivative transformation matrix
!
    if (ldsym) call deftmat
!
!  Output optimisation parameters, if main call
!
    if (lmain) call defoin
!
!  Allocate hessian
!
    if (.not.lconj.or.lminch) then
      if (lhess2D) then
        if (llbfgs) then
          ndwords = nvar*(2*lmbfgsorder + 1) + 2*lmbfgsorder
!
!  Check if minimiser might change and need more memory
!
          if (lminch.and.mintype.le.4) then
            nduwords = nvaronnode
          else
            nduwords = 1_i4
          endif
        else
          ndwords = nvar
          nduwords = nvaronnode
        endif
      else
        if (llbfgs) then
          ndwords = nvar*(2*lmbfgsorder + 1) + 2*lmbfgsorder
        else
          ndwords = 1 + nvar*(nvar+1)/2
        endif
        if ((lrfo.and.nupdate.ne.1).or.(lminch.and.mintype.eq.2)) ndwords = 2*ndwords
      endif
      allocate(hess(ndwords,nduwords),stat=status)
      if (status/=0) call outofmemory('defopt','hess')
    endif
!
!  Start minimisation
!
    call minimise(xc,fc,gc,hess,ndwords,lhess2D,ifail,2_i4,.false.)
!
!  Dellocate hessian
!
    if (.not.lconj.or..not.lminch) then
      deallocate(hess,stat=status)
      if (status/=0) call deallocate_error('defopt','hess')
    endif
!
    if (lfreq.and.loptiloc) then
      lopt = .false.
      ld2sym = .false.
      iflag = 2
      call deffun(iflag,nvar,xc,fc,gc)
      lopt = .true.
    endif
!
!  Substitute parameters into place
!
    call defvartocfg(nvar,xc)
!
!  Now apply constraints
!
    csft(1) = xdc
    csft(2) = ydc
    csft(3) = zdc
    if (ndcon.gt.0) then
      do i = 1,ndcon
        nf = ncdfixind(i)
        if (ncdfixtyp(i).eq.idopt_dx) then
          xdefe(nf) = csft(1)
        elseif (ncdfixtyp(i).eq.idopt_dy) then
          ydefe(nf) = csft(2)
        elseif (ncdfixtyp(i).eq.idopt_dz) then
          zdefe(nf) = csft(3)
        endif
      enddo
      do i = 1,ndcon
        nf = ncdfixind(i)
        if (ncdfixtyp(i).eq.idopt_dx) then
          x0nf = xdefe(nf)
        elseif (ncdfixtyp(i).eq.idopt_dy) then
          x0nf = ydefe(nf)
        elseif (ncdfixtyp(i).eq.idopt_dz) then
          x0nf = zdefe(nf)
        endif
!
        nv = ncdvarind(i)
        if (ncdvartyp(i).eq.idopt_dx) then
          x0nf = (xdefe(nv)-csft(1))*dconco(i) + x0nf
        elseif (ncdvartyp(i).eq.idopt_dy) then
          x0nf = (ydefe(nv)-csft(2))*dconco(i) + x0nf
        elseif (ncdvartyp(i).eq.idopt_dz) then
          x0nf = (zdefe(nv)-csft(3))*dconco(i) + x0nf
        endif
!
        if (ncdfixtyp(i).eq.idopt_dx) then
          xdefe(nf) = x0nf
        elseif (ncdfixtyp(i).eq.idopt_dy) then
          ydefe(nf) = x0nf
        elseif (ncdfixtyp(i).eq.idopt_dz) then
          zdefe(nf) = x0nf
        endif
      enddo
    endif
!
!  If main optimisation then save region 1 to disk for restarts
!
    if (lmain) then
      allocate(itmp(3*ndasym),stat=status)
      if (status/=0) call outofmemory('defopt','itmp')
      if (ldsym) then
        do i = 1,3*ndasym
          itmp(i) = 0
        enddo
        do i = 1,nvar
          ii = idoptindex(i)
          if (idopttype(i).eq.idopt_dx) then
            itmp(3*ii-2) = 1
          elseif (idopttype(i).eq.idopt_dy) then
            itmp(3*ii-1) = 1
          elseif (idopttype(i).eq.idopt_dz) then
            itmp(3*ii) = 1
          endif
        enddo
!
!  As defect constraints are not output we need to set flags for 
!  symmetry constrained coordinates to 1 as well.
!
        do i = 1,ndcon
          nf = ncdfixind(i)
          if (ncdfixtyp(i).eq.idopt_dx) then
            ii = itmp(3*nf-2)
          elseif (ncdfixtyp(i).eq.idopt_dy) then
            ii = itmp(3*nf-1)
          elseif (ncdfixtyp(i).eq.idopt_dz) then
            ii = itmp(3*nf)
          endif
          nv = ncdvarind(i)
          if (ncdvartyp(i).eq.idopt_dx) then
            itmp(3*nv-2) = ii
          elseif (ncdvartyp(i).eq.idopt_dy) then
            itmp(3*nv-1) = ii
          elseif (ncdvartyp(i).eq.idopt_dz) then
            itmp(3*nv) = ii
          endif
        enddo
      else
        do i = 1,3*nreg1
          itmp(i) = 0
        enddo
        do i = 1,nvar
          ii = idoptindex(i)
          if (idopttype(i).eq.idopt_dx) then
            itmp(3*ii-2) = 1
          elseif (idopttype(i).eq.idopt_dy) then
            itmp(3*ii-1) = 1
          elseif (idopttype(i).eq.idopt_dz) then
            itmp(3*ii) = 1
          endif
        enddo
      endif
!
!  Find start of atoms in defect configuration arrays
!
      ndind = 0
      do i = 1,ncf-1
        ndind = ndind + nreg1cfg(i)
      enddo
!
!  Copy information back to configuration arrays
!
      do i = 1,nreg1
        natdefecfg(ndind+i) = natdefe(i)
        ntypdefecfg(ndind+i) = ntypdefe(i)
        xdefecfg(ndind+i) = xdefe(i)
        ydefecfg(ndind+i) = ydefe(i)
        zdefecfg(ndind+i) = zdefe(i)
        qdefecfg(ndind+i) = qdefe(i)
        occdefecfg(ndind+i) = occdefe(i)
        radefecfg(ndind+i) = radefe(i)
        ndefmolcfg(ndind+i) = ndefmol(i)
        ndefindcfg(ndind+i) = ndefind(i)
        ldefbsmatcfg(ndind+i) = ldefbsmat(i)
        ldqmatomcfg(ndind+i) = ldqmatom(i)
        if (ldsym) then
!
!  For symmetry related species need to transform flags
!
          ii = ndrel(i)
          ii = 3*(ii-1)
          ix = itmp(ii+1)
          iy = itmp(ii+2)
          iz = itmp(ii+3)
          ip = ndrelop(i)
          idoptcfg(1,ndind+i) = abs(nint(ix*dsymop(1,1,ip) + iy*dsymop(1,2,ip) + iz*dsymop(1,3,ip)))
          idoptcfg(2,ndind+i) = abs(nint(ix*dsymop(2,1,ip) + iy*dsymop(2,2,ip) + iz*dsymop(2,3,ip)))
          idoptcfg(3,ndind+i) = abs(nint(ix*dsymop(3,1,ip) + iy*dsymop(3,2,ip) + iz*dsymop(3,3,ip)))
        else
          idoptcfg(1,ndind+i) = itmp(3*(i-1)+1)
          idoptcfg(2,ndind+i) = itmp(3*(i-1)+2)
          idoptcfg(3,ndind+i) = itmp(3*(i-1)+3)
        endif
      enddo
      deallocate(itmp,stat=status)
      if (status/=0) call deallocate_error('defopt','itmp')
    endif
  endif
!**********************************
!  Output arc file or stop movie  *
!**********************************
  if (ioproc) then
    if (larc) then
      if (lmovie) then
        close(16)
      else
        call outarc(16_i4,.false.,.true.)
        close(16)
      endif
    endif
    if (lxyz) then
      if (lxyzmovie) then
        close(18)
      else
        call outxyz(18_i4,.false.,.true.)
        close(18)
      endif
    endif
  endif
!************************
!  End of optimisation  *
!************************
  if (lmain) then
    call defout(ifail,fc,gc,xdstore,ydstore,zdstore)
    if (lfreq.and.loptiloc) then
      if (nprocs.gt.1) then
        call deffreqd(.true.,fc)
      else
        call deffreq(.true.,fc)
      endif
    endif
  endif
!
!  Free module memory
!
  deallocate(zdstore,stat=status)
  if (status/=0) call deallocate_error('defopt','zdstore')
  deallocate(ydstore,stat=status)
  if (status/=0) call deallocate_error('defopt','ydstore')
  deallocate(xdstore,stat=status)
  if (status/=0) call deallocate_error('defopt','xdstore')
  deallocate(xc,stat=status)
  if (status/=0) call deallocate_error('defopt','xc')
  deallocate(gc,stat=status)
  if (status/=0) call deallocate_error('defopt','gc')
  call realloc(hess,1_i4,1_i4,ierror)
  if (ierror.ne.0) call outofmemory('defopt','hess')
!
!  Timing
!
  time1 = g_cpu_time()
  time1 = time1 - time0
  if (ioproc) then
    if (lmain.and.loptiloc) then
      write(ioout,'(/,''  Time to end of optimisation  =  '',f12.4,'' seconds'',/)') time1
    elseif (lmain.and.lfreq) then
      write(ioout,'(/,''  Time to end of properties  =  '',f12.4,'' seconds'',/)') time1
    elseif (lmain.and.lgradloc) then
      write(ioout,'(/,''  Time to end of gradients  =  '',f12.4,'' seconds'',/)') time1
    endif
  endif
#ifdef TRACE
  call trace_out('defopt')
#endif
!
  return
  end
