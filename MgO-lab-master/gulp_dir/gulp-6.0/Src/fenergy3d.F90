  subroutine fenergy3d(fc,lgrad1,hessian,nhwords,lkeepd2)
!
!  Supplies the function and first derivatives of the free energy for a solid.
!  Distributed memory parallel version.
!
!  NOTE : BSM or partial occupancy not allowed for!
!
!   4/17 Created from fenergy3
!   5/17 Trap for case where no modes are local to the current node added
!   5/17 Call to transmatd added for parallel case
!   7/17 nshell and nshellonnode added to complex subroutine call
!   7/17 Copy back from eigc to derv2/dervi extended to shells to fix bug
!   7/17 Generation of Pns modified to remove bug
!   8/17 Routine wrapped in ifdef since it won't work or be called without MPI
!   1/18 Trace added
!   3/18 dynam shortened to dyna in keyword check
!   3/18 Parallel I/O corrected
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!   4/20 maxfqat changed to maxfreq since it is no longer just atom based
!   4/20 Mass arrays now use fmass and rfmass
!   7/20 Separate routine for sumall with 1 argument added
!
!  nkpt    = total number of k points across all structures
!  nlkpt   = pointer to lowest k point
!  nukpt   = pointer to upper k point
!  xkpt    = fractional x component of k point
!  ykpt    = fractional y component of k point
!  zkpt    = fractional z component of k point
!  wkpt    = weight of each k point
!  nkptcfg = configuration pointer for each k point
!  sumwkpt = sum over weights of k points
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
!
  use g_constants
  use control
  use current
  use derivatives
  use element
  use energies
  use feworkspace
  use four
  use frequencies
#ifdef MPI
  use general,     only : nwarn
#endif
  use iochannels
  use ksample
  use m_three
  use parallel
  use partial
  use shells
#ifdef MPI
  use species,     only : massspec, natspec, ntypspec
#endif
  use sutton
  use times
#ifdef TRACE
  use trace,       only : trace_in, trace_out
#endif
  implicit none
#ifdef MPI
  include 'mpif.h'
#endif
!
!  Passed variables
!
  integer(i4),  intent(in)                      :: nhwords
  logical,      intent(in)                      :: lgrad1
  logical,      intent(in)                      :: lkeepd2
  real(dp),     intent(inout)                   :: fc
  real(dp)                                      :: hessian(nhwords,*)
#ifdef MPI
!
!  Local variables
!
  character(len=5)                              :: lab1
  complex(dpc), dimension(:,:), allocatable     :: eigc
  integer(i4)                                   :: i
  integer(i4)                                   :: ifail
  integer(i4)                                   :: ii
  integer(i4)                                   :: iloc
  integer(i4)                                   :: ilower
  integer(i4)                                   :: ind
  integer(i4)                                   :: indi
  integer(i4)                                   :: indiloc
  integer(i4)                                   :: irem
  integer(i4)                                   :: ix
  integer(i4)                                   :: j
  integer(i4)                                   :: k
  integer(i4)                                   :: maxeigc
  integer(i4)                                   :: mcv
  integer(i4)                                   :: mcv4
  integer(i4)                                   :: mcvloc
  integer(i4), dimension(:), allocatable        :: mptr
  integer(i4), dimension(:), allocatable        :: mnodeptr
  integer(i4), dimension(:), allocatable        :: mlocptr
  integer(i4)                                   :: mcvmax
  integer(i4)                                   :: mcvmaxloc
  integer(i4)                                   :: mcvmin
  integer(i4)                                   :: mcvminloc
  integer(i4)                                   :: mcvmintot
  integer(i4)                                   :: mint
  integer(i4)                                   :: mintloc
  integer(i4)                                   :: msv
  integer(i4)                                   :: msvloc
  integer(i4)                                   :: nk
  integer(i4)                                   :: nlkpt
  integer(i4)                                   :: nsi
  integer(i4)                                   :: nukpt
  integer(i4)                                   :: status
!
  integer                                       :: idesc(9)
  integer                                       :: idesd(9)
  integer                                       :: ifails
  integer                                       :: ld
  integer                                       :: nb
  integer                                       :: ncs
  integer                                       :: MPIerror
  integer                                       :: ntag
  integer                                       :: nnode
  integer                                       :: ntmp
  integer                                       :: Request
  integer,      dimension(:), allocatable       :: StatMPI       ! Array for status from MPI
!
  logical                                       :: lfirstkpt
  logical                                       :: lmanybodyderv
  logical                                       :: lnozeropt
  real(dp)                                      :: cmfact
  real(dp)                                      :: g_cpu_time
  real(dp),     dimension(:),   allocatable     :: dsqh
  real(dp),     dimension(:),   allocatable     :: dfdw2
  real(dp),     dimension(:),   allocatable     :: dtmp
  real(dp),     dimension(:,:), allocatable     :: eigi
  real(dp),     dimension(:,:), allocatable     :: eigr
  real(dp)                                      :: esum
  real(dp)                                      :: exptrm
  real(dp)                                      :: fetrm
  real(dp)                                      :: fscale
  real(dp)                                      :: rkptcmfct
  real(dp)                                      :: rkptfct
  real(dp)                                      :: rmassi
  real(dp)                                      :: rnokpt
  real(dp)                                      :: sumwkpt
  real(dp)                                      :: t1i
  real(dp)                                      :: t2i
  real(dp)                                      :: wk
  real(dp),     dimension(:,:), allocatable     :: wrk2
  real(dp),     dimension(:,:), allocatable     :: wrk3
  real(dp)                                      :: zpe
#ifdef TRACE
  call trace_in('fenergy3d')
#endif
!
!  Check that there are no breathing shells
!
  if (nbsmat.gt.0) then
    call outerror('breathing shells not allowed for FEM',0_i4)
    call stopnow('fenergy3d')
  endif
!**********************************
!  Set local variables and flags  *
!**********************************
  mint = 3*numat
  mintloc = 3*natomsonnode
  msv = 3*nshell
  msvloc = 3*nshellonnode
  mcv = 3*ncore
  mcvloc = 3*ncoreonnode
!
!  Allocate arrays that depends on mcv
!
  allocate(mptr(mint),stat=status)
  if (status/=0) call outofmemory('fenergy3d','mptr')
  allocate(mnodeptr(mint),stat=status)
  if (status/=0) call outofmemory('fenergy3d','mnodeptr')
  allocate(mlocptr(mintloc),stat=status)
  if (status/=0) call outofmemory('fenergy3d','mlocptr')
  allocate(dfdw2(mcvloc),stat=status)
  if (status/=0) call outofmemory('fenergy3d','dfdw2')
!
!  Set up pointer from local modes to global ones
!
  mptr(1:mint) = 0
  do i = 1,mintloc
    ii = (i-1)/3
    irem = i - 3*ii
    ii = ii + 1
    mlocptr(i) = 3*(node2atom(ii)-1) + irem
    mptr(mlocptr(i)) = i
  enddo
!
!  Set pointer to nodes
!
  do i = 1,mint
    ii = (i-1)/3
    irem = i - 3*ii
    ii = ii + 1
    mnodeptr(i) = atom2node(ii)
  enddo
!
!  For ZSISA find the local value that is greater than or equal to 4
!
  if (nprocs.gt.1) then
    mcv4 = 0
    i = 0
    do while (i.lt.mintloc.and.mcv4.eq.0)
      i = i + 1
      if (mlocptr(i).ge.4) mcv4 = i
    enddo
  else
    mcv4 = 4
  endif
!
!  Ensure that frequency array is large enough
!
  call changemaxfreq(3_i4*ncore)
!
  fscale = sqrt(1.0d23*evtoj*avogadro)
  fscale = fscale/(2.0_dp*pi*speedl)
  if (temperature.gt.1.0d-6) then
    cmfact = planck*speedl/(boltz*temperature)
  endif
  evib = 0.0_dp
  lnozeropt = (index(keyword,'noze').ne.0)
  lmanybodyderv = ((nthb+nfor).gt.0.or.lsuttonc)
!*******************************************************
!  Store contents of second derivative arrays on disk  *
!*******************************************************
!
!  Second derivatives of U needed to approximate hessian
!
  if (lkeepd2) then
    do i = 1,mintloc
      do j = 1,mint
        hessian(j,i) = derv2(j,i)
      enddo
    enddo
  endif
!******************************************
!  Allocate memory for 3-body 3rd derivs  *
!******************************************
  if (lmanybodyderv) then
    maxmany = max(2*numat,54)
    if (lsuttonc) then
      maxmany2 = numat*(numat + 1)
    else
      maxmany2 = maxmany
    endif
    call changemaxmany
  endif
!************************************
!  Find K points for configuration  *
!************************************
  nlkpt = 0
  do i = 1,nkpt
    nk = nkptcfg(i)
    if (nlkpt.eq.0.and.nk.eq.ncf) nlkpt = i
    if (nk.eq.ncf) nukpt = i
  enddo
!
!  No k points found for current structure, so return
!
  if (nlkpt.eq.0) goto 999
!********************************************
!  Calculate inverse square root of masses  *
!********************************************
  do i = 1,ncore
    nsi = nspecptr(nrelf2a(i))
    fmass(3*(i-1)+1) = massspec(nsi)*occuf(i)
    fmass(3*(i-1)+2) = massspec(nsi)*occuf(i)
    fmass(3*(i-1)+3) = massspec(nsi)*occuf(i)
    if (abs(fmass(3*(i-1)+1)).lt.1.0d-12) then
      call label(natspec(nsi),ntypspec(nsi),lab1)
      call outerror('mass of species '//lab1//' is zero',0_i4)
      goto 999
    endif
  enddo
  do i = 1,3*ncore
    if (fmass(i).eq.0.0_dp) then
      call outerror('site has total mass of zero in phonon',0_i4)
      call stopnow('fenergy3d')
    endif
    rfmass(i) = 1.0_dp/sqrt(fmass(i))
  enddo
!
!  Check normalisation of weighting factors 
!
  sumwkpt = 0.0_dp
  do i = nlkpt,nukpt
    sumwkpt = sumwkpt + wkpt(i)
  enddo
  rnokpt = 1.0_dp/sumwkpt
  rkptfct = rnokpt*boltz*temperature/evtoj
  rkptcmfct = rnokpt*planck*speedl/evtoj
!***************************************
!  Allocate memory for complex matrix  *
!***************************************
  maxeigc = maxd2
  allocate(eigc(maxeigc,maxd2u),stat=status)
  if (status/=0) call outofmemory('fenergy3d','eigc')
!
!  Create space for real and imaginary separate copies of eigenvectors to suit old form of code
!
  allocate(eigr(maxeigc,mintloc),stat=status)
  if (status/=0) call outofmemory('fenergy3d','eigr')
  allocate(eigi(maxeigc,mintloc),stat=status)
  if (status/=0) call outofmemory('fenergy3d','eigi')
!
  if (lzsisa) then
    allocate(dsqh(nstrains*maxd2),stat=status)
    if (status/=0) call outofmemory('fenergy3d','dsqh')
  endif
  mcvmintot = 0
!*************************************************
!  Store diagonal blocks to avoid recalculation  *
!*************************************************
  do i = 1,numat
    iloc = atom2local(i)
    if (iloc.gt.0) then
!
!  Atom is local to this node
!
      indi = 3*(i-1)
      indiloc = 3*(iloc-1)
      diagblock(indiloc+1,1) = derv2(indi+1,indiloc+1)
      diagblock(indiloc+2,1) = derv2(indi+2,indiloc+1)
      diagblock(indiloc+3,1) = derv2(indi+3,indiloc+1)
      diagblock(indiloc+1,2) = derv2(indi+1,indiloc+2)
      diagblock(indiloc+2,2) = derv2(indi+2,indiloc+2)
      diagblock(indiloc+3,2) = derv2(indi+3,indiloc+2)
      diagblock(indiloc+1,3) = derv2(indi+1,indiloc+3)
      diagblock(indiloc+2,3) = derv2(indi+2,indiloc+3)
      diagblock(indiloc+3,3) = derv2(indi+3,indiloc+3)
    endif
  enddo
!****************************************************************
!  For ZSISA approximation mode build strain correction matrix  *
!****************************************************************
  if (lzsisa) then
    call setquasiharm(dsqh,mint,mintloc,mcv4)
  endif
!***********************
!  Loop over k points  *
!***********************
  do k = nlkpt,nukpt
!***************************************
!  Generate phased second derivatives  *
!***************************************
    lfirstkpt = (k.eq.nlkpt)
    call dynamic(xkpt(k),ykpt(k),zkpt(k))
!
!  Change sign of imaginary part of dynamical matrix
!
    do i = 1,mintloc
      do j = 1,mint
        dervi(j,i) = - dervi(j,i)
      enddo
    enddo
    wk = wkpt(k)
!
!  Include diagonal blocks, stored in diagblock
!
    do i = 1,numat
      indi = 3*(i-1)
      iloc = atom2local(i)
      if (iloc.gt.0) then
!
!  Atom is local to this node
!
        indi = 3*(i-1)
        indiloc = 3*(iloc-1)
        derv2(indi+1,indiloc+1) = derv2(indi+1,indiloc+1) + diagblock(indiloc+1,1)
        derv2(indi+2,indiloc+1) = derv2(indi+2,indiloc+1) + diagblock(indiloc+2,1)
        derv2(indi+3,indiloc+1) = derv2(indi+3,indiloc+1) + diagblock(indiloc+3,1)
        derv2(indi+1,indiloc+2) = derv2(indi+1,indiloc+2) + diagblock(indiloc+1,2)
        derv2(indi+2,indiloc+2) = derv2(indi+2,indiloc+2) + diagblock(indiloc+2,2)
        derv2(indi+3,indiloc+2) = derv2(indi+3,indiloc+2) + diagblock(indiloc+3,2)
        derv2(indi+1,indiloc+3) = derv2(indi+1,indiloc+3) + diagblock(indiloc+1,3)
        derv2(indi+2,indiloc+3) = derv2(indi+2,indiloc+3) + diagblock(indiloc+2,3)
        derv2(indi+3,indiloc+3) = derv2(indi+3,indiloc+3) + diagblock(indiloc+3,3)
      endif
    enddo
!**********************************
!  Eliminate shell contributions  *
!**********************************
    if (msv.gt.0) then
!*****************************
!  Complex Matrix Inversion  *
!*****************************
      t1i = g_cpu_time()
!
!  Transfer real and imaginary components to complex matrix 
!  NB: must be in the same position due to data decomposition being
!      position dependent
!
      do i = 1,msvloc
        do j = 1,msv
          eigc(mcv+j,mcvloc+i) = dcmplx(derv2(mcv+j,mcvloc+i),dervi(mcv+j,mcvloc+i))
        enddo
      enddo
!
!  Call library to invert matrix stored in eigc
!
      call cmatrix_inversion_shells(msv,mcv+1_i4,maxeigc,eigc,nshell,nshellonnode,ifail)
!
!  Check return flag
!
      if (ifail.ne.0) then
        call outerror('inversion of shell 2nd derivatives failed',0_i4)
        call stopnow('fenergy3d')
      endif
!
!  Transfer data back
!
      do i = 1,msvloc
        do j = 1,msv
          derv2(mcv+j,mcvloc+i) = dble(eigc(mcv+j,mcvloc+i))
          dervi(mcv+j,mcvloc+i) = dimag(eigc(mcv+j,mcvloc+i))
        enddo
      enddo
      t2i = g_cpu_time()
      tmati = tmati + t2i - t1i
!***********************************************
!  Corrected second derivatives = R - T*S-1*T  *
!***********************************************
!
!  Set up Blacs descriptors for matrices
!
      nb = nblocksize
      ifails = 0
      ncs = mcv + msv
      ld = maxeigc
      call descinit( idesc, ncs, ncs, 3*nb, 3*nb, 0, 0, iBlacsContext, ld, ifails )
!
      if (ifails.ne.0) then
        call outerror('initialisation in descinit failed',0_i4)
        call stopnow('fenergy3d')
      endif
!
      ld = maxd2
      call descinit( idesd, ncs, ncs, 3*nb, 3*nb, 0, 0, iBlacsContext, ld, ifails )
!
      if (ifails.ne.0) then
        call outerror('initialisation in descinit failed',0_i4)
        call stopnow('fenergy3d')
      endif
!
!  Transfer data to eigc as complex version of derv2/dervi
!
      do i = 1,mcvloc+msvloc
        do j = 1,mcv+msv
          eigc(j,i) = dcmplx(derv2(j,i),dervi(j,i))
        enddo
      enddo
!
!  First pass : S-1*T stored in diagonally opposite copy of T
!
      call pzgemm('N','C',msv,mcv,msv,(1.0d0,0.0d0),eigc,mcv+1,mcv+1,idesc,eigc,1,mcv+1, &
                  idesc,(0.0d0,0.0d0),eigc,mcv+1,1,idesc)
!
!  Second pass : T*(S-1*T)
!
      call pzgemm('N','N',mcv,mcv,msv,(-1.0d0,0.0d0),eigc,1,mcv+1,idesc,eigc,mcv+1,1, &
                  idesc,(1.0d0,0.0d0),eigc,1,1,idesc)
!
!  Transfer back from eigc to derv2/dervi - NB important to copy back shell parts too since they are used later
!
      do i = 1,mcvloc+msvloc
        do j = 1,mcv+msv
          derv2(j,i) = dble(eigc(j,i))
          dervi(j,i) = dimag(eigc(j,i))
        enddo
      enddo
    endif
!****************************
!  End of shell correction  *
!****************************
!
!  Multiply by mass-factors
!
    ind = 0
    do i = 1,ncoreonnode
      indi = 3*(node2atom(i)-1)
      do ix = 1,3
        ind  = ind + 1
        indi = indi + 1
        rmassi = rfmass(indi)
        do j = 1,3*ncore
          derv2(j,ind) = rmassi*rfmass(j)*derv2(j,ind)
          dervi(j,ind) = rmassi*rfmass(j)*dervi(j,ind)
        enddo
      enddo
    enddo
!*********************
!  Debugging output  *
!*********************
    if (index(keyword,'dyna').ne.0) then
      if (ioproc) then
        write(ioout,'(/,''  Real Dynamical matrix :'',/)')
      endif
      if (lioproconly) then
!
!  Allocate temporary workspace for communication
!
        ntmp = mcv
        ntag = 1
        allocate(dtmp(mcv),stat=status)
        if (status/=0) call outofmemory('fenergy3d','dtmp')
        allocate(StatMPI(MPI_Status_Size),stat=status)
        if (status/=0) call outofmemory('fenergy3d','StatMPI')
!
        do i = 1,mcv
          iloc = mptr(i)
          if (iloc.gt.0) then
            if (mnodeptr(i).ne.0_i4) then
!
!  Post receive
!
              if (ioproc) then
                nnode = mnodeptr(i)
                call MPI_IRecv(dtmp,ntmp,MPI_double_precision,nnode, &
                               ntag,MPI_Comm_World,Request,MPIerror)
              endif
!
!  Pass data to ioproc for writing
!
              if (iloc.gt.0) then
                dtmp(1:ntmp) = derv2(1:ntmp,iloc)
!
!  Post send
!
                call MPI_ISend(dtmp,ntmp,MPI_double_precision,0, &
                               ntag,MPI_Comm_World,Request,MPIerror)
              endif
              if (ioproc.or.iloc.gt.0) then
                call MPI_WaitAll(1,Request,StatMPI,MPIerror)
              endif
              if (ioproc) then
!
!  Write on I/O node
!
                write(ioout,'(9f12.8)')(dtmp(j),j=1,ntmp)
              endif
            else
              if (iloc.gt.0) then
                write(ioout,'(9f12.8)')(derv2(j,iloc),j=1,mcv)
              endif
            endif
          endif
        enddo
      else
        call mpbarrier
        do i = 1,mcv
          iloc = mptr(i)
          if (iloc.gt.0) then
            write(ioout,'(9f12.8)')(derv2(j,iloc),j=1,mcv)
          endif
          call mpbarrier
        enddo
      endif
      if (ioproc) then
        write(ioout,'(/,''  Imaginary Dynamical matrix :'',/)')
      endif
      if (lioproconly) then
        do i = 1,mcv
          iloc = mptr(i)
          if (iloc.gt.0) then
            if (mnodeptr(i).ne.0_i4) then
!
!  Post receive
!
              if (ioproc) then
                nnode = mnodeptr(i)
                call MPI_IRecv(dtmp,ntmp,MPI_double_precision,nnode, &
                               ntag,MPI_Comm_World,Request,MPIerror)
              endif
!
!  Pass data to ioproc for writing
!
              if (iloc.gt.0) then
                dtmp(1:ntmp) = dervi(1:ntmp,iloc)
!
!  Post send
!
                call MPI_ISend(dtmp,ntmp,MPI_double_precision,0, &
                               ntag,MPI_Comm_World,Request,MPIerror)
              endif
              if (ioproc.or.iloc.gt.0) then
                call MPI_WaitAll(1,Request,StatMPI,MPIerror)
              endif
              if (ioproc) then
!
!  Write on I/O node
!
                write(ioout,'(9f12.8)')(dtmp(j),j=1,ntmp)
              endif
            else
              if (iloc.gt.0) then
                write(ioout,'(9f12.8)')(dervi(j,iloc),j=1,mcv)
              endif
            endif
          endif
        enddo
!
        deallocate(StatMPI,stat=status)
        if (status/=0) call deallocate_error('fenergy3d','StatMPI')
        deallocate(dtmp,stat=status)
        if (status/=0) call deallocate_error('fenergy3d','dtmp')
      else
        call mpbarrier
        do i = 1,mcv
          iloc = mptr(i)
          if (iloc.gt.0) then
            write(ioout,'(9f12.8)')(dervi(j,iloc),j=1,mcv)
          endif
          call mpbarrier
        enddo
      endif
    endif
!*****************************************************************
!  Diagonalise dynamical matrix => Eigenvectors and eigenvalues  *
!*****************************************************************
    ifail = 0
    call pdiagd(mcv,mcvloc,maxd2,derv2,dervi,maxeigc,eigc,freq(1,k-nlkpt+1),fscale,.true.,.false.,ifail)
!
!  Copy eigenvectors from complex form to real and imaginary components
!
    do i = 1,mcvloc
      do j = 1,mcv
        eigr(j,i) = dble(eigc(j,i))
        eigi(j,i) = aimag(eigc(j,i))
      enddo
    enddo
!
!  Set minimum and maximum mode numbers
!
    mcvmin = minmode
    ilower = mcvmin
    do i = ilower,mcv
      if (freq(i,k-nlkpt+1).lt.1.0_dp) mcvmin = mcvmin + 1
    enddo
!
    mcvmintot = mcvmintot + mcvmin
    mcvmax = mcv
    if (minmode.ne.1) mcvmin = minmode
    if (maxmode.ne.0) mcvmax = maxmode
    if (index(keyword,'verb').ne.0.and.ioproc) then
      write(ioout,'(/,''  Minimum mode number for K point '',i5,'' = '',i5)') k-nlkpt+1,mcvmin
      write(ioout,'(''  Maximum mode number for K point '',i5,'' = '',i5,/)') k-nlkpt+1,mcvmax
    endif
!
!  Find local versions of mcvmin and mcvmax
!
    mcvmaxloc = 0
    mcvminloc = 0
    i = 0
    do while (i.lt.mcvloc.and.mcvminloc.eq.0)
      i = i + 1
      if (mlocptr(i).ge.mcvmin) then
        mcvminloc = i
      endif
    enddo
    i = mcvloc
    do while (i.gt.0.and.mcvmaxloc.eq.0)
      if (mlocptr(i).le.mcvmax) then
        if (mlocptr(i).ge.mcvmin) then
          mcvmaxloc = i
        endif
      endif
      i = i - 1
    enddo
!
!  Trap case where no modes are on the local mode
!
    if (mcvminloc.eq.0.and.mcvmaxloc.eq.0) mcvmaxloc = -1
!**************************************
!  Calculate transformation matrices  *
!  and contributions to free energy   *
!**************************************
    fetrm = 0.0_dp
    zpe = 0.0_dp
!
!  Calculate (1/2w).(dG/dw) => dfdw2
!
    if (lnozeropt) then
      do iloc = mcvminloc,mcvmaxloc
        i = mlocptr(iloc)
        if (temperature.gt.1.0d-6) then
          exptrm = exp(-freq(i,k-nlkpt+1)*cmfact)
        else
          exptrm = 0.0_dp
        endif
        if (exptrm.lt.1.0_dp) then
          fetrm = fetrm + log(1.0_dp-exptrm)
          dfdw2(iloc) = wk*rkptcmfct*(exptrm/(1.0_dp-exptrm))/(2.0_dp*freq(i,k-nlkpt+1))
!
!  Sqrt dfdw2 so that it can be convoluted into the eigenvectors 
!  to save time on the projection of the third derivatives
!
          dfdw2(iloc) = sqrt(dfdw2(iloc))*fscale
        else
          dfdw2(iloc) = 0.0_dp
        endif
      enddo
    else
      do iloc = mcvminloc,mcvmaxloc
        i = mlocptr(iloc)
        if (temperature.gt.1.0d-6) then
          exptrm = exp(-freq(i,k-nlkpt+1)*cmfact)
        else
          exptrm = 0.0_dp
        endif
        zpe = zpe + wk*freq(i,k-nlkpt+1)*rkptcmfct
        if (exptrm.lt.1.0_dp) then
          fetrm = fetrm + log(1.0_dp-exptrm)
          dfdw2(iloc) = wk*rkptcmfct*(0.5_dp+exptrm/(1.0_dp-exptrm))/(2.0_dp*freq(i,k-nlkpt+1))
!
!  Sqrt dfdw2 so that it can be convoluted into the eigenvectors 
!  to save time on the projection of the third derivatives
!
          dfdw2(iloc) = sqrt(dfdw2(iloc))*fscale
        else
          dfdw2(iloc) = 0.0_dp
        endif
      enddo
    endif
    zpe = 0.5_dp*zpe
    evib = evib + zpe + fetrm*wk*rkptfct
!
    if (lgrad1) then
!*********************
!  Debugging output  *
!*********************
      if (index(keyword,'dyna').ne.0) then
        if (ioproc) then
          write(ioout,'(/,''  Real Eigenvectors :'',/)')
        endif
        if (lioproconly) then
!
!  Allocate temporary workspace for communication
!
          ntmp = mcvmax
          ntag = 1
          allocate(dtmp(mcvmax),stat=status)
          if (status/=0) call outofmemory('fenergy3d','dtmp')
          allocate(StatMPI(MPI_Status_Size),stat=status)
          if (status/=0) call outofmemory('fenergy3d','StatMPI')
!
          do i = 1,mcv
            iloc = mptr(i)
            if (iloc.gt.0) then
              if (mnodeptr(i).ne.0_i4) then
!
!  Post receive
!
                if (ioproc) then
                  nnode = mnodeptr(i)
                  call MPI_IRecv(dtmp,ntmp,MPI_double_precision,nnode, &
                                 ntag,MPI_Comm_World,Request,MPIerror)
                endif
!
!  Pass data to ioproc for writing
!
                if (iloc.gt.0) then
                  dtmp(1:ntmp) = eigr(1:ntmp,iloc)
!
!  Post send
!
                  call MPI_ISend(dtmp,ntmp,MPI_double_precision,0, &
                                 ntag,MPI_Comm_World,Request,MPIerror)
                endif
                if (ioproc.or.iloc.gt.0) then
                  call MPI_WaitAll(1,Request,StatMPI,MPIerror)
                endif
                if (ioproc) then
!
!  Write on I/O node
!
                  write(ioout,'(12f11.6)') freq(i,k-nlkpt+1),(dtmp(j),j=1,mcvmax)
                endif
              else
                if (iloc.gt.0) then
                  write(ioout,'(12f11.6)') freq(i,k-nlkpt+1),(eigr(j,iloc),j=1,mcvmax)
                endif
              endif
            endif
          enddo
        else
          call mpbarrier
          do i = 1,mcvmax
            iloc = mptr(i)
            if (iloc.gt.0) then
              write(ioout,'(12f11.6)') freq(i,k-nlkpt+1),(eigr(j,iloc),j=1,mcvmax)
            endif
            call mpbarrier
          enddo
        endif
        if (ioproc) then
          write(ioout,'(/,''  Imag Eigenvectors :'',/)')
        endif
        if (lioproconly) then
          do i = 1,mcv
            iloc = mptr(i)
            if (iloc.gt.0) then
              if (mnodeptr(i).ne.0_i4) then
!
!  Post receive
!
                if (ioproc) then
                  nnode = mnodeptr(i)
                  call MPI_IRecv(dtmp,ntmp,MPI_double_precision,nnode, &
                                 ntag,MPI_Comm_World,Request,MPIerror)
                endif
!
!  Pass data to ioproc for writing
!
                if (iloc.gt.0) then
                  dtmp(1:ntmp) = eigi(1:ntmp,iloc)
!
!  Post send
!
                  call MPI_ISend(dtmp,ntmp,MPI_double_precision,0, &
                                 ntag,MPI_Comm_World,Request,MPIerror)
                endif
                if (ioproc.or.iloc.gt.0) then
                  call MPI_WaitAll(1,Request,StatMPI,MPIerror)
                endif
                if (ioproc) then
!
!  Write on I/O node
!
                  write(ioout,'(12f11.6)') freq(i,k-nlkpt+1),(dtmp(j),j=1,mcvmax)
                endif
              else
                if (iloc.gt.0) then
                  write(ioout,'(12f11.6)') freq(i,k-nlkpt+1),(eigi(j,iloc),j=1,mcvmax)
                endif
              endif
            endif
          enddo
!
          deallocate(StatMPI,stat=status)
          if (status/=0) call deallocate_error('fenergy3d','StatMPI')
          deallocate(dtmp,stat=status)
          if (status/=0) call deallocate_error('fenergy3d','dtmp')
        else
          call mpbarrier
          do i = 1,mcvmax
            iloc = mptr(i)
            if (iloc.gt.0) then
              write(ioout,'(12f11.6)')freq(i,k-nlkpt+1),(eigi(j,iloc),j=1,mcvmax)
            endif
            call mpbarrier
          enddo
        endif
      endif
!
!  Scale eigenvectors by sqrt(dfdw2)
!
      call scaleevec(mcv,mcvloc,eigr,maxeigc,dfdw2)
      call scaleevec(mcv,mcvloc,eigi,maxeigc,dfdw2)
!
!  Calculate shell projection matrices
!
!  Two components to store real/imaginary
!  
!  Msc => already stored in derv2(mcv+j,i)/dervi(mcv+j,i)
!  Pns => store in derv2(i,mcv+j)/dervi(i,mcv+j)
!      =  Enc.Mcs, where Enc = eigenvector for mode n
!
      if (msv.gt.0) then
!
!  Scale Msc by mass factor for use in derivatives
!
        ind = 0
        do i = 1,ncoreonnode
          indi = 3*(node2atom(i)-1)
          do ix = 1,3
            ind  = ind + 1
            indi = indi + 1
            rmassi = rfmass(indi)
            do j = 1,msv
              derv2(mcv+j,ind) = rmassi*derv2(mcv+j,ind)
              dervi(mcv+j,ind) = rmassi*dervi(mcv+j,ind)
            enddo
          enddo
        enddo
!*********************
!  Debugging output  *
!*********************
        if (index(keyword,'dyna').ne.0) then
          if (ioproc) then
            write(ioout,'(/,''  Msc (real) :'',/)')
          endif
          if (lioproconly) then
!
!  Allocate temporary workspace for communication
!
            ntmp = msv
            ntag = 1
            allocate(dtmp(msv),stat=status)
            if (status/=0) call outofmemory('fenergy3d','dtmp')
            allocate(StatMPI(MPI_Status_Size),stat=status)
            if (status/=0) call outofmemory('fenergy3d','StatMPI')
!
            do i = 1,mcv
              iloc = mptr(i)
              if (iloc.gt.0) then
                if (mnodeptr(i).ne.0_i4) then
!
!  Post receive
!
                  if (ioproc) then
                    nnode = mnodeptr(i)
                    call MPI_IRecv(dtmp,ntmp,MPI_double_precision,nnode, &
                                   ntag,MPI_Comm_World,Request,MPIerror)
                  endif
!
!  Pass data to ioproc for writing
!
                  if (iloc.gt.0) then
                    dtmp(1:ntmp) = derv2(mcv+1:mcv+ntmp,iloc)
!
!  Post send
!
                    call MPI_ISend(dtmp,ntmp,MPI_double_precision,0, &
                                   ntag,MPI_Comm_World,Request,MPIerror)
                  endif
                  if (ioproc.or.iloc.gt.0) then
                    call MPI_WaitAll(1,Request,StatMPI,MPIerror)
                  endif
                  if (ioproc) then
!
!  Write on I/O node
!
                    write(ioout,'(12f11.6)')(dtmp(j),j=1,ntmp)
                  endif
                else
                  if (iloc.gt.0) then
                    write(ioout,'(12f11.6)')(derv2(mcv+j,iloc),j=1,msv)
                  endif
                endif
              endif
            enddo
          else
            call mpbarrier
            do i = 1,mcv
              iloc = mptr(i)
              if (iloc.gt.0) then
                write(ioout,'(12f11.6)')(derv2(mcv+j,iloc),j=1,msv)
              endif
              call mpbarrier
            enddo
          endif
          call gflush(ioout)
          if (ioproc) then
            write(ioout,'(/,''  Msc (imag) :'',/)')
          endif
          if (lioproconly) then
            do i = 1,mcv
              iloc = mptr(i)
              if (iloc.gt.0) then
                if (mnodeptr(i).ne.0_i4) then
!
!  Post receive
!
                  if (ioproc) then
                    nnode = mnodeptr(i)
                    call MPI_IRecv(dtmp,ntmp,MPI_double_precision,nnode, &
                                   ntag,MPI_Comm_World,Request,MPIerror)
                  endif
!
!  Pass data to ioproc for writing
!
                  if (iloc.gt.0) then
                    dtmp(1:ntmp) = dervi(mcv+1:mcv+ntmp,iloc)
!
!  Post send
!
                    call MPI_ISend(dtmp,ntmp,MPI_double_precision,0, &
                                   ntag,MPI_Comm_World,Request,MPIerror)
                  endif
                  if (ioproc.or.iloc.gt.0) then
                    call MPI_WaitAll(1,Request,StatMPI,MPIerror)
                  endif
                  if (ioproc) then
!
!  Write on I/O node
!
                    write(ioout,'(12f11.6)')(dtmp(j),j=1,ntmp)
                  endif
                else
                  if (iloc.gt.0) then
                    write(ioout,'(12f11.6)')(dervi(mcv+j,iloc),j=1,msv)
                  endif
                endif
              endif
            enddo
!
            deallocate(StatMPI,stat=status)
            if (status/=0) call deallocate_error('fenergy3d','StatMPI')
            deallocate(dtmp,stat=status)
            if (status/=0) call deallocate_error('fenergy3d','dtmp')
          else
            call mpbarrier
            do i = 1,mcv
              iloc = mptr(i)
              if (iloc.gt.0) then
                write(ioout,'(12f11.6)')(dervi(mcv+j,iloc),j=1,msv)
              endif
              call mpbarrier
            enddo
          endif
          call gflush(ioout)
        endif
!
!  Generate Pns from Msc and eigenvectors
!  NB position in derv2 of final matrix is flipped relative to old position
!
        allocate(wrk2(maxd2,mcvloc),stat=status)
        if (status/=0) call outofmemory('fenergy3d','wrk2')
        allocate(wrk3(maxd2,mcvloc),stat=status)
        if (status/=0) call outofmemory('fenergy3d','wrk3')
!
        call pdgemm('N','N',msv,mcv,mcv,1.0d0,derv2,mcv+1,1,idesd,eigr,1,1,idesc,0.0d0,wrk2,mcv+1,1,idesd)
        call pdgemm('N','N',msv,mcv,mcv,-1.0d0,dervi,mcv+1,1,idesd,eigi,1,1,idesc,1.0d0,wrk2,mcv+1,1,idesd)
!
        call pdgemm('N','N',msv,mcv,mcv,1.0d0,derv2,mcv+1,1,idesd,eigi,1,1,idesc,0.0d0,wrk3,mcv+1,1,idesd)
        call pdgemm('N','N',msv,mcv,mcv,1.0d0,dervi,mcv+1,1,idesd,eigr,1,1,idesc,1.0d0,wrk3,mcv+1,1,idesd)
!
        do i = 1,mcvloc
          do j = 1,msv
            derv2(mcv+j,i) = wrk2(mcv+j,i)
            dervi(mcv+j,i) = wrk3(mcv+j,i)
          enddo
        enddo
!
        deallocate(wrk3,stat=status)
        if (status/=0) call deallocate_error('fenergy3d','wrk3')
        deallocate(wrk2,stat=status)
        if (status/=0) call deallocate_error('fenergy3d','wrk2')
!*********************
!  Debugging output  *
!*********************
        if (index(keyword,'dyna').ne.0) then
          if (ioproc) then
            write(ioout,'(/,''  Pns (real) :'',/)')
          endif
          if (lioproconly) then
!
!  Allocate temporary workspace for communication
!
            ntmp = msv
            ntag = 1
            allocate(dtmp(msv),stat=status)
            if (status/=0) call outofmemory('fenergy3d','dtmp')
            allocate(StatMPI(MPI_Status_Size),stat=status)
            if (status/=0) call outofmemory('fenergy3d','StatMPI')
!
            do i = 1,mcv
              iloc = mptr(i)
              if (iloc.gt.0) then
                if (mnodeptr(i).ne.0_i4) then
!
!  Post receive
!
                  if (ioproc) then
                    nnode = mnodeptr(i)
                    call MPI_IRecv(dtmp,ntmp,MPI_double_precision,nnode, &
                                   ntag,MPI_Comm_World,Request,MPIerror)
                  endif
!
!  Pass data to ioproc for writing
!
                  if (iloc.gt.0) then
                    dtmp(1:ntmp) = derv2(mcv+1:mcv+ntmp,iloc)
!
!  Post send
!
                    call MPI_ISend(dtmp,ntmp,MPI_double_precision,0, &
                                   ntag,MPI_Comm_World,Request,MPIerror)
                  endif
                  if (ioproc.or.iloc.gt.0) then
                    call MPI_WaitAll(1,Request,StatMPI,MPIerror)
                  endif
                  if (ioproc) then
!
!  Write on I/O node
!
                    write(ioout,'(12f11.6)')(dtmp(j),j=1,ntmp)
                  endif
                else
                  if (iloc.gt.0) then
                    write(ioout,'(12f11.6)')(derv2(mcv+j,iloc),j=1,msv)
                  endif
                endif
              endif
            enddo
          else
            call mpbarrier
            do i = 1,mcv
              iloc = mptr(i)
              if (iloc.gt.0) then
                write(ioout,'(12f11.6)')(derv2(mcv+j,iloc),j=1,msv)
              endif
              call mpbarrier
            enddo
          endif
          call gflush(ioout)
          if (ioproc) then
            write(ioout,'(/,''  Pns (imag) :'',/)')
          endif
          if (lioproconly) then
            do i = 1,mcv
              iloc = mptr(i)
              if (iloc.gt.0) then
                if (mnodeptr(i).ne.0_i4) then
!
!  Post receive
!
                  if (ioproc) then
                    nnode = mnodeptr(i)
                    call MPI_IRecv(dtmp,ntmp,MPI_double_precision,nnode, &
                                   ntag,MPI_Comm_World,Request,MPIerror)
                  endif
!
!  Pass data to ioproc for writing
!
                  if (iloc.gt.0) then
                    dtmp(1:ntmp) = dervi(mcv+1:mcv+ntmp,iloc)
!
!  Post send
!
                    call MPI_ISend(dtmp,ntmp,MPI_double_precision,0, &
                                   ntag,MPI_Comm_World,Request,MPIerror)
                  endif
                  if (ioproc.or.iloc.gt.0) then
                    call MPI_WaitAll(1,Request,StatMPI,MPIerror)
                  endif
                  if (ioproc) then
!
!  Write on I/O node
!
                    write(ioout,'(12f11.6)')(dtmp(j),j=1,ntmp)
                  endif
                else
                  if (iloc.gt.0) then
                    write(ioout,'(12f11.6)')(dervi(mcv+j,iloc),j=1,msv)
                  endif
                endif
              endif
            enddo
!
            deallocate(StatMPI,stat=status)
            if (status/=0) call deallocate_error('fenergy0d','StatMPI')
            deallocate(dtmp,stat=status)
            if (status/=0) call deallocate_error('fenergy0d','dtmp')
          else
            call mpbarrier
            do i = 1,mcv
              iloc = mptr(i)
              if (iloc.gt.0) then
                write(ioout,'(12f11.6)')(dervi(mcv+j,iloc),j=1,msv)
              endif
              call mpbarrier
            enddo
          endif
          call gflush(ioout)
        endif
!
!  End of shell condition
!
      endif
!********************************
!  Evaluate phonon derivatives  *
!********************************
      call realrecip3d3d(k,mcvminloc,mcvmaxloc,eigr,eigi,dsqh,maxd2,lfirstkpt)
    endif
!***************************
!  End loop over K points  *
!***************************
  enddo
  if (nummode.eq.0) then
    nummode = mcvmintot
  elseif (nummode.ne.mcvmintot) then
    if (ioproc) then
      write(ioout,'(''**** Warning - number of modes has changed to '',i8,'' ****'')')(mcv+1)*nkpt-mcvmintot
    endif
    nwarn = nwarn + 1
    nummode = mcvmintot
  endif
!
!  Global sum of vibrational energy
!
  call sumone(evib,esum,"fenergy3d","evib")
  evib = esum
!
!  Sum vibration contribution and internal energy
!
  fc = fc + evib
!***************
!  Exit point  *
!***************
999 continue
!****************
!  Free memory  *
!****************
  if (lzsisa.and.allocated(dsqh)) then
    deallocate(dsqh,stat=status)
    if (status/=0) call deallocate_error('fenergy3d','dsqh')
  endif
  deallocate(eigi,stat=status)
  if (status/=0) call deallocate_error('fenergy3d','eigi')
  deallocate(eigr,stat=status)
  if (status/=0) call deallocate_error('fenergy3d','eigr')
  deallocate(eigc,stat=status)
  if (status/=0) call deallocate_error('fenergy3d','eigc')
!
  deallocate(dfdw2,stat=status)
  if (status/=0) call deallocate_error('fenergy3d','dfdw2')
  deallocate(mlocptr,stat=status)
  if (status/=0) call deallocate_error('fenergy3d','mlocptr')
  deallocate(mnodeptr,stat=status)
  if (status/=0) call deallocate_error('fenergy3d','mnodeptr')
  deallocate(mptr,stat=status)
  if (status/=0) call deallocate_error('fenergy3d','mptr')
!*************************************************
!  Recover contents of second derivative arrays  *
!*************************************************
!
!  Re-generate transformation matrix if needed for optimisation
!  This is generally faster than saving it to disk
!
  if ((lopt.or.lrelax).and.lkeepd2) then
    call transmatd
  endif
!
!  Second derivatives of U needed to approximate hessian
!
  if (lkeepd2) then
    do i = 1,mintloc
      do j = 1,mint
        derv2(j,i) = hessian(j,i)
      enddo
    enddo
  endif
#ifdef TRACE
  call trace_out('fenergy3d')
#endif
#else
  call outerror('fenergy3d called when not compiled with MPI',0_i4)
  call stopnow('fenergy3d')
#endif
!
  return
  end
