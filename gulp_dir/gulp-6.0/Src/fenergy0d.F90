  subroutine fenergy0d(fc,lgrad1,hessian,nhwords,lkeepd2)
!
!  Supplies the function and first derivatives of the free energy
!  for a molecule.
!  Distributed memory parallel version.
!  NB: BSM or partial occupancy not currently allowed
!
!   4/17 Created from fenergy0
!   5/17 Trap for case where no modes are local to the current node added
!   5/17 Call to transmatd added for parallel case
!   5/17 nshell and nshellonnode added to subroutine call
!   8/17 Routine wrapped in ifdef since it won't work or be called without MPI
!   1/18 Trace added
!   3/18 dynam shortened to dyna in keyword check
!   3/18 Parallel I/O corrected
!   8/19 fc changed from intent in to inout
!   4/20 maxfqat changed to maxfreq since it is no longer just atom based
!   4/20 Mass arrays now use fmass and rfmass
!   7/20 Separate routine for sumall with 1 argument added
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
  use general,      only : nwarn
#endif
  use iochannels
  use m_three
  use parallel
  use partial
  use shells
#ifdef MPI
  use species,      only : massspec, natspec, ntypspec
#endif
  use sutton
  use times
#ifdef TRACE
  use trace,        only : trace_in, trace_out
#endif
  implicit none
#ifdef MPI
  include 'mpif.h'
#endif
!
!  Passed variables
!
  integer(i4),                   intent(in)    :: nhwords
  logical,                       intent(in)    :: lgrad1
  logical,                       intent(in)    :: lkeepd2
  real(dp),                      intent(inout) :: fc
  real(dp)                                     :: hessian(nhwords,*)
#ifdef MPI
!
!  Local variables
!
  character(len=5)                             :: lab1
  integer(i4)                                  :: i
  integer(i4)                                  :: ifail
  integer(i4)                                  :: ii
  integer(i4)                                  :: iloc
  integer(i4)                                  :: ind
  integer(i4)                                  :: indi
  integer(i4)                                  :: irem
  integer(i4)                                  :: ix
  integer(i4)                                  :: j
  integer(i4)                                  :: lwrk
  integer(i4)                                  :: mcv
  integer(i4)                                  :: mcvloc
  integer(i4), dimension(:), allocatable       :: mptr
  integer(i4), dimension(:), allocatable       :: mnodeptr
  integer(i4), dimension(:), allocatable       :: mlocptr
  integer(i4)                                  :: mcvmax
  integer(i4)                                  :: mcvmaxloc
  integer(i4)                                  :: mcvmin
  integer(i4)                                  :: mcvminloc
  integer(i4)                                  :: mint
  integer(i4)                                  :: mintloc
  integer(i4)                                  :: msv
  integer(i4)                                  :: msvloc
  integer(i4)                                  :: nimag
  integer(i4)                                  :: nsi
  integer(i4)                                  :: status
!
  integer                                      :: idesc(9)
  integer                                      :: idescs(9)
  integer                                      :: ifails
  integer                                      :: ld
  integer                                      :: n
  integer                                      :: nb
  integer                                      :: ncs
  integer                                      :: nloc
  integer                                      :: nsize
  integer                                      :: MPIerror
  integer                                      :: ntag
  integer                                      :: nnode
  integer                                      :: ntmp
  integer                                      :: Request
  integer,     dimension(:), allocatable       :: StatMPI       ! Array for status from MPI
!
  logical                                      :: linear
  logical                                      :: lmanybodyderv
  logical                                      :: lnozeropt
  real(dp)                                     :: cmfact
  real(dp)                                     :: g_cpu_time
  real(dp),    dimension(:), allocatable       :: dfdw2
  real(dp),    dimension(:), allocatable       :: dtmp
  real(dp)                                     :: esum
  real(dp)                                     :: exptrm
  real(dp)                                     :: fetrm
  real(dp)                                     :: fscale
  real(dp)                                     :: rkptcmfct
  real(dp)                                     :: rkptfct
  real(dp)                                     :: rmassi
  real(dp)                                     :: rt
  real(dp)                                     :: t1d
  real(dp)                                     :: t1i
  real(dp)                                     :: t2d
  real(dp)                                     :: t2i
  real(dp),    dimension(:),   allocatable     :: wrk
  real(dp),    dimension(:,:), allocatable     :: wrk2
  real(dp)                                     :: zpe
#ifdef TRACE
  call trace_in('fenergy0d')
#endif
!
!  Check that there are no breathing shells
!
  if (nbsmat.gt.0) then
    call outerror('breathing shells not allowed for FEM',0_i4)
    call stopnow('fenergy0d')
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
  if (status/=0) call outofmemory('fenergy0d','mptr')
  allocate(mnodeptr(mint),stat=status)
  if (status/=0) call outofmemory('fenergy0d','mnodeptr')
  allocate(mlocptr(mintloc),stat=status)
  if (status/=0) call outofmemory('fenergy0d','mlocptr')
  allocate(dfdw2(mcvloc),stat=status)
  if (status/=0) call outofmemory('fenergy0d','dfdw2')
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
!  Ensure that frequency array is large enough
!
  call changemaxfreq(3_i4*ncore)
!
  fscale = sqrt(1.0d23*evtoj*avogadro)
  fscale = fscale/(2.0_dp*pi*speedl)
  if (temperature.gt.1.0d-6) then
    cmfact = planck*speedl/(boltz*temperature)
  endif
  rkptfct = boltz*temperature/evtoj
  rkptcmfct = planck*speedl/evtoj
  lmanybodyderv = ((nthb+nfor).gt.0.or.lsuttonc)
  lnozeropt = (index(keyword,'noze').ne.0)
!***********************************************
!  Store contents of second derivative arrays  *
!***********************************************
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
!*********************************************
!  Allocate memory for many-body 3rd derivs  *
!*********************************************
  if (lmanybodyderv) then
    maxmany = numat
    if (lsuttonc) then
      maxmany2 = numat*(numat - 1)/2
    else
      maxmany2 = numat
    endif
    call changemaxmany
  endif
!*********************
!  Debugging output  *
!*********************
  if (index(keyword,'dyna').ne.0) then
    if (ioproc) then
      write(ioout,'(/,''  Second derivative matrix:'',/)')
    endif
    if (lioproconly) then
!
!  Allocate temporary workspace for communication
!
      ntmp = mcv + msv
      ntag = 1
      allocate(dtmp(mcv+msv),stat=status)
      if (status/=0) call outofmemory('fenergy0d','dtmp')
      allocate(StatMPI(MPI_Status_Size),stat=status)
      if (status/=0) call outofmemory('fenergy0d','StatMPI')
!
      do i = 1,mcv+msv
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
              write(ioout,'(12f11.6)')(dtmp(j),j=1,ntmp)
            endif
          else
            if (iloc.gt.0) then
              write(ioout,'(12f11.6)')(derv2(j,iloc),j=1,mcv+msv)
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
      do i = 1,mcv+msv
        iloc = mptr(i)
        if (iloc.gt.0) then
          write(ioout,'(12f11.6)')(derv2(j,iloc),j=1,mcv+msv)
        endif
        call mpbarrier
      enddo
    endif
  endif
!********************************************
!  Calculate inverse square root of masses  *
!********************************************
  do i = 1,ncore
    nsi = nspecptr(i)
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
      call stopnow('fenergy0d')
    endif
    rfmass(i) = 1.0_dp/sqrt(fmass(i))
  enddo
!**********************************
!  Eliminate shell contributions  *
!**********************************
  if (msv.gt.0) then
!*********************
!  Matrix Inversion  *
!*********************
    t1i = g_cpu_time()
!
!  Call library to invert matrix stored in eigr
!
    call matrix_inversion_shells(msv,mcv+1_i4,maxd2,derv2,nshell,nshellonnode,ifail)
!
!  Check return flag
!
    if (ifail.ne.0) then
      call outerror('inversion of shell 2nd derivatives failed',0_i4)
      call stopnow('fenergy0d')
    endif
!
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
    ld = maxd2
    call descinit( idescs, ncs, ncs, 3*nb, 3*nb, 0, 0, iBlacsContext, ld, ifails )
!
    if (ifails.ne.0) then
      call outerror('initialisation in descinit failed',0_i4)
      call stopnow('fenergy0d')
    endif
!
!  First pass : S-1*T stored in diagonally opposite copy of T
!
    call pdgemm('N','T',msv,mcv,msv,1.0d0,derv2,mcv+1,mcv+1,idescs,derv2,1,mcv+1,idescs,0.0d0,derv2,mcv+1,1,idescs)
!
!  Second pass : T*(S-1*T)
!
    call pdgemm('N','N',mcv,mcv,msv,-1.0d0,derv2,1,mcv+1,idescs,derv2,mcv+1,1,idescs,1.0d0,derv2,1,1,idescs)
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
      if (status/=0) call outofmemory('fenergy0d','dtmp')
      allocate(StatMPI(MPI_Status_Size),stat=status)
      if (status/=0) call outofmemory('fenergy0d','StatMPI')
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
              write(ioout,'(12f11.6)')(dtmp(j),j=1,ntmp)
            endif
          else
            if (iloc.gt.0) then
              write(ioout,'(12f11.6)')(derv2(j,iloc),j=1,mcv)
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
          write(ioout,'(12f11.6)')(derv2(j,iloc),j=1,mcv)
        endif
        call mpbarrier
      enddo
    endif
  endif
!*****************************************************************
!  Diagonalise dynamical matrix => eigenvectors and eigenvalues  *
!*****************************************************************
  t1d = g_cpu_time()
!
!  Set up Blacs descriptors for arrays
!
  n = mcv
  nloc = mcvloc
  nb = 3*nblocksize
  ld = maxd2
  call descinit(idesc, n, n, nb, nb, 0, 0, iBlacsContext, ld, ifails)
!
  if (ifails.ne.0) then
    call outerror('initialisation in descinit failed - idesc ',0_i4)
    call stopnow('fenergy0d')
  endif
!
!  Find size of work space array
!
  nsize = max((nb*(nb-1))/2,(n+nloc)*nb) + nb*nb
  lwrk = 5*n + n*max(1,nloc) + max(2*n,nsize) + 1
!
!  Allocate work space
!
  allocate(wrk(lwrk),stat=status)
  if (status/=0) call outofmemory('fenergy0d','wrk')
  allocate(wrk2(maxd2,mcvloc),stat=status)
  if (status/=0) call outofmemory('fenergy0d','wrk2')
!
!  Make copy of derv2 to avoid overwriting
!
  do i = 1,mcvloc
    do j = 1,mcv
      wrk2(j,i) = derv2(j,i)
    enddo
  enddo
!
!  Call Scalapack eigensolver
!
  call pdsyev( 'V', 'U', n, wrk2, 1, 1, idesc, freq, dervi, 1, 1, idesc, wrk, lwrk, ifails )
!
  ifail = ifails
!
!  Free workspace
!
  deallocate(wrk2,stat=status)
  if (status/=0) call deallocate_error('fenergy0d','wrk2')
  deallocate(wrk,stat=status)
  if (status/=0) call deallocate_error('fenergy0d','wrk')
!
  t2d = g_cpu_time()
  tdiag = tdiag + t2d - t1d
!***************************************
!  Convert eigenvalues to frequencies  *
!***************************************
!
!  For molecule exclude rotations
!
  call moltype(linear)
  if (linear) then
    mcvmin = 6
  else
    mcvmin = 7
  endif
!
!  Check for imaginary modes and exclude them
!
!  Note that rotations tend to be imaginary modes 
!
  nimag = 0
  do i = 1,mcv
    if (freq(i,1).lt.-1.0d-3) nimag = nimag + 1
  enddo
  if (nimag.gt.3) mcvmin = mcvmin + (nimag - 3)
  mcvmax = mcv
  if (minmode.ne.1) mcvmin = minmode
  if (maxmode.ne.0) mcvmax = maxmode
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
!
  do i = mcvmin,mcvmax
    rt = freq(i,1)
    if (rt.ge.0.0_dp) then
      freq(i,1) = sqrt(rt)*fscale
    else
      rt = abs(rt)
      freq(i,1) = - sqrt(rt)*fscale
    endif
  enddo
  if (nummode.eq.0) then
    nummode = mcvmin
  elseif (nummode.ne.mcvmin) then
    nwarn = nwarn + 1
    if (ioproc) then
      write(ioout,'(''**** Warning - number of modes has changed during free energy minimisation ****'')')
    endif
  endif
  if (index(keyword,'dyna').ne.0.and.ioproc) then
    write(ioout,'(/,''  Frequencies for modes included in free energy (cm^-1) :'',/)')
    do i = mcvmin,mcvmax
      write(ioout,'(i8,2x,f15.6)') i,freq(i,1)
    enddo
  endif
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
        exptrm = exp(-freq(i,1)*cmfact)   
      else
        exptrm = 0.0_dp
      endif
      if (exptrm.lt.1.0_dp) then
        fetrm = fetrm + log(1.0_dp - exptrm)
        dfdw2(iloc) = rkptcmfct*(exptrm/(1.0_dp - exptrm))/(2.0_dp*freq(i,1))
        dfdw2(iloc) = sqrt(dfdw2(iloc))*fscale
      else
        dfdw2(iloc) = 0.0_dp
      endif
    enddo
  else 
    do iloc = mcvminloc,mcvmaxloc
      i = mlocptr(iloc)
      if (temperature.gt.1.0d-6) then
        exptrm = exp(-freq(i,1)*cmfact)   
      else
        exptrm = 0.0_dp
      endif
      zpe = zpe + freq(i,1)*rkptcmfct
      if (exptrm.lt.1.0_dp) then
        fetrm = fetrm + log(1.0_dp-exptrm)
        dfdw2(iloc) = rkptcmfct*(0.5_dp + exptrm/(1.0_dp - exptrm))/(2.0_dp*freq(i,1))
        dfdw2(iloc) = sqrt(dfdw2(iloc))*fscale
      else
        dfdw2(iloc) = 0.0_dp
      endif
    enddo
  endif
  zpe = 0.5_dp*zpe
  evib = zpe + fetrm*rkptfct
!
!  Global sum of vibrational energy
!
  call sumone(evib,esum,"fenergy0d","evib")
  evib = esum
!
  fc = fc + evib
  if (lgrad1) then
!*********************
!  Debugging output  *
!*********************
    if (index(keyword,'dyna').ne.0) then
      if (ioproc) then
        write(ioout,'(/,''  Eigenvectors :'',/)')
      endif
      if (lioproconly) then
!
!  Allocate temporary workspace for communication
!
        ntmp = mcv
        ntag = 1
        allocate(dtmp(mcv),stat=status)
        if (status/=0) call outofmemory('fenergy0d','dtmp')
        allocate(StatMPI(MPI_Status_Size),stat=status)
        if (status/=0) call outofmemory('fenergy0d','StatMPI')
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
                write(ioout,'(12f11.6)')(dtmp(j),j=1,ntmp)
              endif
            else
              if (iloc.gt.0) then
                write(ioout,'(12f11.6)')(derv2(j,iloc),j=1,mcv)
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
        do i = mcvmin,mcvmax
          iloc = mptr(i)
          if (iloc.gt.0) then
            write(ioout,'(12f11.6)')(dervi(j,iloc),j=1,mcv)
          endif
          call mpbarrier
        enddo
      endif
    endif
!********************************
!  Evaluate phonon derivatives  *
!********************************
!
!  Scale eigenvectors by sqrt(dfdw2)
!
    call scaleevec(mcv,mcvloc,dervi,maxd2,dfdw2)
!
!  Calculate shell projection matrices 
!
!  Msc => already stored in derv2(mcv+j,i)
!  Pns => store in derv2(i,mcv+j)
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
          enddo
        enddo
      enddo
!*********************
!  Debugging output  *
!*********************
      if (index(keyword,'dyna').ne.0) then
        if (ioproc) then
          write(ioout,'(/,''  Msc :'',/)')
        endif
        if (lioproconly) then
!
!  Allocate temporary workspace for communication
!
          ntmp = msv
          ntag = 1
          allocate(dtmp(msv),stat=status)
          if (status/=0) call outofmemory('fenergy0d','dtmp')
          allocate(StatMPI(MPI_Status_Size),stat=status)
          if (status/=0) call outofmemory('fenergy0d','StatMPI')
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
              write(ioout,'(12f11.6)')(derv2(mcv+j,iloc),j=1,msv)
            endif
            call mpbarrier
          enddo
        endif
      endif
!
!  Generate Pns from Msc and eigenvectors
!  NB position in derv2 of final matrix is flipped relative to old position
!
      allocate(wrk2(maxd2,mcvloc),stat=status)
      if (status/=0) call outofmemory('fenergy0d','wrk2')
!
      call pdgemm('N','N',msv,mcv,mcv,1.0d0,derv2,mcv+1,1,idescs,dervi,1,1,idescs,0.0d0,wrk2,1,1,idesc)
!
!  Copy wrk2 back to derv2
!
      do i = 1,mcvloc
        do j = 1,msv
          derv2(mcv+j,i) = wrk2(j,i)
        enddo
      enddo
!
      deallocate(wrk2,stat=status)
      if (status/=0) call deallocate_error('fenergy0d','wrk2')
!*********************
!  Debugging output  *
!*********************
      if (index(keyword,'dyna').ne.0) then
        if (ioproc) then
          write(ioout,'(/,''  Pns :'',/)')
        endif
        if (lioproconly) then
!
!  Allocate temporary workspace for communication
!
          ntmp = msv
          ntag = 1
          allocate(dtmp(msv),stat=status)
          if (status/=0) call outofmemory('fenergy0d','dtmp')
          allocate(StatMPI(MPI_Status_Size),stat=status)
          if (status/=0) call outofmemory('fenergy0d','StatMPI')
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
              write(ioout,'(12f11.6)')(derv2(mcv+j,iloc),j=1,msv)
            endif
            call mpbarrier
          enddo
        endif
      endif
!
!  End of shell condition
!
    endif
!
!  Real space component
!
    call real0d3d(mcvminloc,mcvmaxloc)
  endif
!***************
!  Exit point  *
!***************
999 continue
!****************
!  Free memory  *
!****************
  deallocate(dfdw2,stat=status)
  if (status/=0) call deallocate_error('fenergy0d','dfdw2')
  deallocate(mlocptr,stat=status)
  if (status/=0) call deallocate_error('fenergy0d','mlocptr')
  deallocate(mnodeptr,stat=status)
  if (status/=0) call deallocate_error('fenergy0d','mnodeptr')
  deallocate(mptr,stat=status)
  if (status/=0) call deallocate_error('fenergy0d','mptr')
!*************************************************
!  Recover contents of second derivative arrays  *
!*************************************************
!
!  Re-generate transformation matrix if needed for optimisation
!  This is in general faster than saving to disk
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
  call trace_out('fenergy0d')
#endif
#else
  call outerror('fenergy0d called when not compiled with MPI',0_i4)
  call stopnow('fenergy0d')
#endif
!
  return
  end
