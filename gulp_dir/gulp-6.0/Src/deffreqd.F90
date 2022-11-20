  subroutine deffreqd(lprint,fc)
!
!  Calculates defect frequencies.
!  Distributed memory parallel version.
!
!  NB: BSM and partial occupancy not currently supported
!
!  The block diagonal elements of derv2 must be stored, as these are
!  not recalculated in dynamic as the summing of off diagonals is no
!  longer applicable.
!
!  leigloc = local flag to indicate whether eigenvectors are 
!            to be generated for this configuration
!  fhenergy= Helmholtz free-energy
!  rtlnz   = R*T*ln(z), where z=partition function
!
!   5/17 Created from deffreq
!   5/17 nobsmode removed since this no longer used since clusters stopped using deffreq
!   8/17 Routine wrapped in ifdef since it won't work or be called without MPI
!  10/17 fhenergy moved to energies module
!  11/17 Output format tweaked
!   2/18 Trace added
!   3/18 dynam shortened to dyna in keyword check
!   3/18 Parallel I/O corrected
!   4/18 leigloc now set here for llower and linten
!   3/19 Multiple temperature ramps added
!   5/19 Finite difference flag split for first and second derivatives
!   4/20 maxfqat changed to maxfreq since it is no longer just atom based
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
!  Julian Gale, CIC, Curtin University, April 2020
!
  use g_constants
  use control
  use current
  use defects
  use derivatives, vectors => dervi
  use element
#ifdef MPI
  use energies,    only : fhenergy
#endif
  use gulp_files
  use frequencies
#ifdef MPI
  use general,     only : lfinitediff2
#endif
  use genetic
  use iochannels
  use parallel
  use projectdos
  use shells
  use species
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
  logical,      intent(in)                     :: lprint        ! If true then output results
  real(dp),     intent(in)                     :: fc            ! Internal energy
#ifdef MPI
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4), dimension(:),     allocatable   :: ibocptr
  integer(i4)                                  :: ifail
  integer(i4)                                  :: igroup
  integer(i4)                                  :: ii
  integer(i4)                                  :: iloc
  integer(i4)                                  :: imoddb
  integer(i4)                                  :: inat
  integer(i4)                                  :: ind
  integer(i4)                                  :: indi
  integer(i4)                                  :: indj
  integer(i4), dimension(:),     allocatable   :: iocptr
  integer(i4)                                  :: iresid
  integer(i4), dimension(:),     allocatable   :: itmp
  integer(i4)                                  :: itype
  integer(i4)                                  :: j
  integer(i4)                                  :: jj
  integer(i4)                                  :: k
  integer(i4)                                  :: l
  integer(i4)                                  :: m
  integer(i4)                                  :: maxlim
  integer(i4)                                  :: mcv
  integer(i4)                                  :: mcvloc
  integer(i4),  dimension(:),     allocatable  :: mcvptr
  integer(i4),  dimension(:),     allocatable  :: mcvnptr
  integer(i4),  dimension(:),     allocatable  :: mcvrptr
  integer(i4)                                  :: mcvmax
  integer(i4)                                  :: mcvmin
  integer(i4)                                  :: mint
  integer(i4)                                  :: mintloc
  integer(i4)                                  :: mm
  integer(i4)                                  :: msv
  integer(i4)                                  :: msvloc
  integer(i4)                                  :: nati
  integer(i4)                                  :: natj
  integer(i4), dimension(:),     allocatable   :: ndbsptr
  integer(i4)                                  :: nbfoc
  integer(i4)                                  :: ndbs
  integer(i4)                                  :: ncfoc
  integer(i4)                                  :: ni
  integer(i4)                                  :: nimag
  integer(i4)                                  :: np
  integer(i4)                                  :: npc
  integer(i4)                                  :: npfirst
  integer(i4)                                  :: npi
  integer(i4)                                  :: npifirst
  integer(i4)                                  :: npilast
  integer(i4)                                  :: nplast
  integer(i4)                                  :: nproj
  integer(i4)                                  :: nr1
  integer(i4)                                  :: nsfoc
  integer(i4)                                  :: nsi
  integer(i4)                                  :: nt
  integer(i4)                                  :: nt0
  integer(i4)                                  :: ntj
  integer(i4)                                  :: ntr
  integer(i4)                                  :: status
!
  integer                                      :: MPIerror
  integer                                      :: idesc(9)
  integer                                      :: ifails
  integer                                      :: ld
  integer                                      :: nb
  integer                                      :: ncs
  integer                                      :: nnode
  integer                                      :: ntag
  integer                                      :: ntmp
  integer                                      :: Request
  integer,     dimension(:),   allocatable     :: StatMPI       ! Array for status from MPI
!
  logical                                      :: lbsmi
  logical                                      :: lbsmj
  logical                                      :: lcorei
  logical                                      :: lcorej
  logical                                      :: leigloc
  logical                                      :: lfound
  logical                                      :: lnozero
  logical                                      :: lpocc
  logical                                      :: lproj
  logical                                      :: lprinloc
  real(dp)                                     :: cmfact
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: cv
  real(dp),    dimension(:,:),   allocatable   :: dtmp
  real(dp)                                     :: cv2
  real(dp)                                     :: ent2
  real(dp)                                     :: entropy
  real(dp)                                     :: factor
  real(dp)                                     :: fe_equipartition
  real(dp)                                     :: freqmin
  real(dp)                                     :: fscale
  real(dp),    dimension(:),     allocatable   :: irx
  real(dp),    dimension(:),     allocatable   :: iry
  real(dp),    dimension(:),     allocatable   :: irz
  real(dp),    dimension(:),     allocatable   :: massdef
  real(dp)                                     :: oci
  real(dp)                                     :: qj
  real(dp)                                     :: rkt
  real(dp),    dimension(:),     allocatable   :: rmassdef
  real(dp)                                     :: rmassi
  real(dp)                                     :: rmode
  real(dp)                                     :: rsum
  real(dp),    dimension(:),     allocatable   :: rtmp2
  real(dp)                                     :: rtlnz
  real(dp)                                     :: s_equipartition
  real(dp)                                     :: t1d
  real(dp)                                     :: t1i
  real(dp)                                     :: t1t
  real(dp)                                     :: t2d
  real(dp)                                     :: t2i
  real(dp)                                     :: t2t
  real(dp)                                     :: tem
  real(dp)                                     :: trm
  real(dp)                                     :: trm1
  real(dp)                                     :: trmcv
  real(dp)                                     :: trmen
  real(dp)                                     :: trmfe
  real(dp)                                     :: trmfe_eq
  real(dp)                                     :: trmj
  real(dp)                                     :: trms_eq
  real(dp)                                     :: trmzp
  real(dp),    dimension(:),     allocatable   :: w1
  real(dp),    dimension(:),     allocatable   :: w2
  real(dp),    dimension(:),     allocatable   :: w3
  real(dp)                                     :: xd
  real(dp)                                     :: yd
  real(dp)                                     :: zd
  real(dp)                                     :: xi
  real(dp)                                     :: yi
  real(dp)                                     :: zi
  real(dp)                                     :: xir
  real(dp)                                     :: yir
  real(dp)                                     :: zir
  real(dp)                                     :: zpe
#ifdef TRACE
  call trace_in('deffreqd')
#endif
!
  t1t = g_cpu_time()
  leigloc = leigen
  if (linten) leigloc = .true.
  if (llower) leigloc = .true.
  lnozero  = (index(keyword,'noze').ne.0)
  lproj = ((nprojcfg(ncf)-nprojdef(ncf)).gt.0) 
  if (lproj) leigloc = .true.
  if (.not.lprint) leigloc = .false.
  lprinloc = (lprint.and.ioproc)
  fscale = sqrt(1.0d23*evtoj*avogadro)
  fscale = fscale/(2.0_dp*pi*speedl)
!
!  If analytic second derivatives are not available (or finite differences are requested) then try finite differences
!
  if (lnoanald2.or.lfinitediff2) then
    call dynamicn_defect
  endif
!
!  Set mode for projection of DOS
!
  imoddb = 2
!***********
!  Defect  *
!***********
!
!  Allocate pointer arrays
!
  allocate(ndbsptr(nreg1),stat=status)
  if (status/=0) call outofmemory('deffreqd','ndbsptr')
  allocate(ibocptr(nreg1),stat=status)
  if (status/=0) call outofmemory('deffreqd','ibocptr')
  allocate(iocptr(nreg1),stat=status)
  if (status/=0) call outofmemory('deffreqd','iocptr')
!
  ndbs = 0
  do i = 1,nreg1
    if (ldefbsmat(i)) then
      ndbs = ndbs + 1
      ndbsptr(ndbs) = i
    endif
  enddo
!
!  Setup partial occupancy pointer
!
  ncfoc = 0
  nsfoc = 0
  nbfoc = 0
  do i = 1,nreg1
    ibocptr(i) = 0
    if (occdefe(i).eq.1.0_dp.and..not.lallowgt1) then
!
!  Fully occupied site
!
      if (natdefe(i).gt.maxele) then
        nsfoc = nsfoc + 1
      else
        ncfoc = ncfoc + 1
      endif
      iocptr(i) = ncfoc + nsfoc
      if (ldefbsmat(i)) then
        nbfoc = nbfoc + 1
        ibocptr(i) = nbfoc
      endif
    else
!
!  Partially occupied site
!  Check to see if there is a previous atom on this site
!
      xi = xdefe(i)
      yi = ydefe(i)
      zi = zdefe(i)
      nati = natdefe(i)
      lcorei = (nati.le.maxele)
      lbsmi = (ldefbsmat(i))
      j = 1
      lfound = .false.
      do while (j.lt.i.and..not.lfound)
        lcorej = (natdefe(j).le.maxele)
        if ((lcorei.and.lcorej).or.(.not.lcorei.and..not.lcorej)) then
          xd = xdefe(j) - xi
          yd = ydefe(j) - yi
          zd = zdefe(j) - zi
          rsum = abs(xd) + abs(yd) + abs(zd)
          if (rsum.lt.1.0d-4) then
            lfound = .true.
            iocptr(i) = iocptr(j)
            lbsmj = (ldefbsmat(j))
            if (lbsmi) then
              if (lbsmj) then
                ibocptr(i) = ibocptr(j)
              else
                nbfoc = nbfoc + 1
                ibocptr(i) = nbfoc
              endif
            endif
          endif
        endif
        j = j + 1
      enddo
      if (.not.lfound) then
!
!  Must be new site
!
        if (lcorei) then
          ncfoc = ncfoc + 1
        else
          nsfoc = nsfoc + 1
        endif
        iocptr(i) = ncfoc + nsfoc
        if (ldefbsmat(i)) then
          nbfoc = nbfoc + 1
          ibocptr(i) = nbfoc
        endif
      endif
    endif
  enddo
!
!  Trap cases not currently implemented with distributed memory
!
  lpocc = (nsfoc+ncfoc.ne.nreg1)
  if (lpocc) then
    call outerror('partial occupancy not compatible with distributed 2nd derivs',0_i4)
    call stopnow('deffreqd')
  endif
  if (ndbs.gt.0) then
    call outerror('breathing shells not compatible with distributed 2nd derivs',0_i4)
    call stopnow('deffreqd')
  endif
!
!  End of partial occupancy pointer generation
!
  mint = 3*nreg1
  mintloc = 3*nreg1onnode
  maxlim = mint
  if (ndbs.gt.0) maxlim = maxlim + nreg1
  nr1 = nreg1
  msv = 3*nshreg1
  mcv = 3*ncoreg1
  msvloc = 3*nreg1onnodes
  mcvloc = 3*nreg1onnodec
!
!  Ensure that frequency array is large enough
!
  call changemaxfreq(3_i4*ncoreg1)
!
!  Allocate local memory
!
  allocate(mcvptr(mcv),stat=status)
  if (status/=0) call outofmemory('deffreqd','mcvptr')
  allocate(mcvnptr(mcv),stat=status)
  if (status/=0) call outofmemory('deffreqd','mcvnptr')
  allocate(mcvrptr(mcv),stat=status)
  if (status/=0) call outofmemory('deffreqd','mcvrptr')
  allocate(irx(mcv),stat=status)
  if (status/=0) call outofmemory('deffreqd','irx')
  allocate(iry(mcv),stat=status)
  if (status/=0) call outofmemory('deffreqd','iry')
  allocate(irz(mcv),stat=status)
  if (status/=0) call outofmemory('deffreqd','irz')
  allocate(rtmp2(mcv),stat=status)
  if (status/=0) call outofmemory('deffreqd','rtmp2')
  allocate(w1(3*mcv),stat=status)
  if (status/=0) call outofmemory('deffreqd','w1')
  allocate(w2(3*mcv),stat=status)
  if (status/=0) call outofmemory('deffreqd','w2')
  allocate(w3(mcv),stat=status)
  if (status/=0) call outofmemory('deffreqd','w3')
  allocate(massdef(ncfoc),stat=status)
  if (status/=0) call outofmemory('deffreqd','massdef')
  allocate(rmassdef(ncfoc),stat=status)
  if (status/=0) call outofmemory('deffreqd','rmassdef')
!
!
  if (nprocs.gt.1) then
    mcvrptr(1:mcv) = 0
    do iloc = 1,nreg1onnodec
      i = node2reg1(iloc)
      m = 3*(iloc-1)
      mm = 3*(i-1)
      mcvptr(m+1) = mm + 1
      mcvptr(m+2) = mm + 2
      mcvptr(m+3) = mm + 3
      mcvrptr(mm+1) = m + 1
      mcvrptr(mm+2) = m + 2
      mcvrptr(mm+3) = m + 3
    enddo
!
!  Set up pointer to node for each value of mcvrptr
!
    do i = 1,nreg1
      if (natdefe(i).le.maxele) then
        mm = 3*(i-1)
        mcvnptr(mm+1) = reg12node(i)
        mcvnptr(mm+2) = reg12node(i)
        mcvnptr(mm+3) = reg12node(i)
      endif
    enddo
  else
    do m = 1,mcvloc
      mcvptr(m) = m
      mcvnptr(m) = 0_i4
      mcvrptr(m) = m
    enddo
  endif
!
!  Calculate inversion square root of masses
!
!  Now modified to handle partial occupancies
!
  do i = 1,ncoreg1
    massdef(i) = 0.0_dp
  enddo
  do i = 1,ncoreg1
    ni = natdefe(i)
    nt = ntypdefe(i)
    oci = occdefe(i)
    lfound = .false.
    nsi = 0
    do while (nsi.le.nspec.and..not.lfound)
      nsi = nsi + 1
      lfound = (ni.eq.natspec(nsi).and.nt.eq.ntypspec(nsi))
    enddo
    if (lfound) then
      rmassi = massspec(nsi)*oci
    else
      rmassi = atmass(ni)*oci
    endif
    if (rmassi.eq.0.0) then
      call outerror('mass of element '//atsym(ni)//' is zero',0_i4)
      call stopnow('deffreqd')
    endif
    massdef(i) = massdef(i) + rmassi
  enddo
  do i = 1,ncoreg1
    if (massdef(i).eq.0.0_dp) then
      call outerror('site has total mass of zero in phonon',0_i4)
      call stopnow('deffreqd')
    endif
    rmassdef(i) = 1.0_dp/sqrt(massdef(i))
  enddo
!
!  Output frequency header
!
  if (lprinloc) then
    write(ioout,'(/,''  Vibrational Frequency Calculation : '',/)')
    write(ioout,'(''--------------------------------------------------------------------------------'')')
  endif
!**********************************
!  Eliminate shell contributions  *
!**********************************
  if (msv.gt.0) then
!
!  Invert dynamical shell-shell matrix
!
    t1i = g_cpu_time()
    ifail = 0
!
!  Call library to invert matrix stored in eigr
!
    call matrix_inversion_shells(msv,mcv+1_i4,maxd2,derv2,nshreg1,nreg1onnodes,ifail)
!
!  Check return flag
!
    if (ifail.ne.0) then
      call outerror('inversion of shell 2nd derivatives failed',0_i4)
      call stopnow('deffreqd')
    endif
!
    t2i = g_cpu_time()
    tmati = tmati + t2i - t1i
!*************************************************
!  Corrected second derivatives  =  R - T*S-1*T  *
!*************************************************
!
!  Set up Blacs descriptors for matrices
!
    nb = nblocksize
    ifails = 0
    ncs = mcv + msv
    ld = maxd2
    call descinit( idesc, ncs, ncs, 3*nb, 3*nb, 0, 0, iBlacsContext, ld, ifails )
!
    if (ifails.ne.0) then
      call outerror('initialisation in descinit failed',0_i4)
      call stopnow('deffreqd')
    endif
!
!  First pass : S-1*T stored in diagonally opposite copy of T
!
    call pdgemm('N','T',msv,mcv,msv,1.0d0,derv2,mcv+1,mcv+1,idesc,derv2,1,mcv+1,idesc,0.0d0,derv2,mcv+1,1,idesc)
!
!  Second pass : T*(S-1*T)
!
    call pdgemm('N','N',mcv,mcv,msv,-1.0d0,derv2,1,mcv+1,idesc,derv2,mcv+1,1,idesc,1.0d0,derv2,1,1,idesc)
  endif
!****************************
!  End of shell correction  *
!****************************
!
!  Multiply by mass-factors
!
  do ii = 1,3
    indi = ii - 3
    do jj = 1,3
      indj = jj - 3
      do i = 1,nreg1onnodec
        rmassi = rmassdef(node2reg1(i))
        do j = 1,ncoreg1
          derv2(indj+3*j,indi+3*i) = rmassi*rmassdef(j)*derv2(indj+3*j,indi+3*i)
        enddo
      enddo
    enddo
  enddo
!
!  If debugging print out dynamical matrix
!
  if (index(keyword,'dyna').ne.0) then
    if (ioproc) then
      write(ioout,'(/,''  Real Dynamical matrix :'',/)')
    endif
#ifdef MPI
    if (lioproconly) then
!
!  Allocate temporary workspace for communication
!
      ntmp = 3*mcv
      ntag = 1
      allocate(dtmp(mcv,3_i4),stat=status)
      if (status/=0) call outofmemory('deffreqd','dtmp')
      allocate(StatMPI(MPI_Status_Size),stat=status)
      if (status/=0) call outofmemory('deffreqd','StatMPI')
!
      do i = 1,ncoreg1
        iloc = reg12local(i)
        if (reg12node(i).gt.0) then
!
!  Post receive
!
          if (ioproc) then
            nnode = reg12node(i)
            call MPI_IRecv(dtmp,ntmp,MPI_double_precision,nnode, &
                           ntag,MPI_Comm_World,Request,MPIerror)
          endif
!
!  Pass data to ioproc for writing
!
          if (iloc.gt.0) then
            ind = 3*(iloc-1)
            do ii = 1,3
              do j = 1,mcv
                dtmp(j,ii) = derv2(j,ind+ii)
              enddo
            enddo
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
            do ii = 1,3
              write(ioout,'(12f10.5)')(dtmp(j,ii),j=1,mcv)
            enddo
          endif
        else
          if (iloc.gt.0) then
            ind = 3*(iloc-1)
            do ii = 1,3
              write(ioout,'(12f10.5)')(derv2(j,ind+ii),j = 1,mcv)
            enddo
          endif
        endif
      enddo
!
      deallocate(StatMPI,stat=status)
      if (status/=0) call deallocate_error('deffreqd','StatMPI')
      deallocate(dtmp,stat=status)
      if (status/=0) call deallocate_error('deffreqd','dtmp')
    else
#endif
      call mpbarrier
      do i = 1,ncoreg1
        iloc = reg12local(i)
        if (iloc.gt.0) then
          ind = 3*(iloc-1)
          do ii = 1,3
            write(ioout,'(12f10.5)')(derv2(j,ind+ii),j = 1,mcv)
          enddo
        endif
        call mpbarrier
      enddo
#ifdef MPI
    endif
#endif
  endif
!*********************************
!  Diagonalise dynamical matrix  *
!*********************************
  ifail = 0
  t1d = g_cpu_time()
!
  call pdiaggd(mcv,mcvloc,maxd2,derv2,maxd2,vectors,freq,fscale,leigloc,.true.,ifail)
!
  t2d = g_cpu_time()
  tdiag = tdiag + t2d - t1d
!***********************
!  Output frequencies  *
!***********************
  if (leigloc) then
!**********************************
!  Projected densities of states  *
!**********************************
    if (lproj) then
      allocate(itmp(nr1),stat=status)
      if (status/=0) call outofmemory('deffreqd','itmp')
!
!  Loop over projections
!
      nproj = nprojcfg(ncf)
      npfirst = 1
      npifirst = 1
      ii = 0
      do i = 1,ncf-1
        npc = nprojcfg(i)
        npfirst = npfirst + npc
        do j = 1,npc
          npifirst = npifirst + nprojit(ii+j)
        enddo
        ii = ii + npc
      enddo
      nplast = npfirst + nproj - 1
      npilast = npifirst
      do i = 1,nproj
        npilast = npilast + nprojit(ii+i)
      enddo
      npilast = npilast - 1
      do np = npfirst,nplast
        if (nprojdb(np).eq.imoddb) then
          do i = 1,nr1
            itmp(i) = 0
          enddo
!
!  Find atoms of projection
!
          do npi = npifirst,npilast
            if (nprojptr(npi).eq.np) then
              if (nprojtyp(npi).gt.99) then
                itmp(nprojnat(npi)) = 1
              else
                inat = nprojnat(npi)
                itype = nprojtyp(npi)
                do i = 1,nr1
                  if (inat.eq.natdefe(i).and.(itype.eq.ntypdefe(i).or.itype.eq.0)) itmp(i) = 1
                enddo
              endif
            endif
          enddo
!
!  Loop over frequencies
!
          w1(1:mcv) = 0.0_dp
          do j = 1,mcvloc
            jj = mcvptr(j)
            ind = 0
!
!  Loop over atoms
!
            do k = 1,ncfoc
              lfound = .false.
              l = 1
              do while (l.le.ncoreg1.and..not.lfound)
                if (iocptr(l).eq.k) then
                  if (itmp(l).eq.1) then
                    lfound = .true.
                    w1(jj) = w1(jj) + vectors(ind+1,j)**2 + &
                                      vectors(ind+2,j)**2 + &
                                      vectors(ind+3,j)**2
                  endif
                endif
                l = l + 1
              enddo
              ind = ind + 3
            enddo
          enddo
          call sumall(w1,w2,mcv,"deffreqd","w")
          if (ioproc) then
            do j = 1,mcv
              write(59) w2(j)
            enddo
          endif
        endif
      enddo
      deallocate(itmp,stat=status)
      if (status/=0) call deallocate_error('deffreqd','itmp')
    endif
!*********************************************
!  Evaluate infra-red and Raman intensities  *
!*********************************************
    if (linten) then
!
!  Old algorithm - no use of oscillator strengths
!
      do j = 1,ncoreg1
        w3(j) = 0.0_dp
        qj   = qa(j)
        natj = natdefe(j)
        ntj  = ntypdefe(j)
        trmj = occua(j)/sqrt(atmass(natj))
!
!  Add on any shell charge
!
        do k = 1,nspec
          if ((natspec(k)-maxele).eq.natj) then
            if (ntj.eq.ntypspec(k).or.ntypspec(k).eq.0) then
              qj = qj + qlspec(k)
            endif
          endif
        enddo
!
!  Term for IR
!
        w3(j) = qj*trmj
      enddo
!
!  Sum eigenvector components multiplied by charge
!
      rkt = boltz*temperature
      if (abs(rkt).gt.1.0d-12) then
        cmfact = planck*speedl/rkt
      endif
      do i = 1,mcvloc
        ii = mcvptr(i)
        xir = 0.0_dp
        yir = 0.0_dp
        zir = 0.0_dp
        ind = 0
        do j = 1,ncfoc
          xir = xir + w3(j)*vectors(ind+1,i)
          yir = yir + w3(j)*vectors(ind+2,i)
          zir = zir + w3(j)*vectors(ind+3,i)
          ind = ind + 3
        enddo
        irx(ii) = xir*xir
        iry(ii) = yir*yir
        irz(ii) = zir*zir
      enddo
      call sumall(irx,w1,mcv,"deffreqd","irx")
      irx(1:mcv) = w1(1:mcv)
      call sumall(iry,w1,mcv,"deffreqd","iry")
      iry(1:mcv) = w1(1:mcv)
      call sumall(irz,w1,mcv,"deffreqd","irz")
      irz(1:mcv) = w1(1:mcv)
      do i = 1,mcv
        IRintensity(i,1) = irx(i) + iry(i) + irz(i)
      enddo
    endif
!************************
!  Output eigenvectors  *
!************************
    if (index(keyword,'eige').ne.0) then
      if (ioproc) then
        write(ioout,'(/,''  Frequencies (cm-1) and Eigenvectors : '',/)')
      endif
      call mpbarrier
#ifdef MPI
      if (lioproconly) then
!
!  Allocate temporary workspace for communication
!
        ntag = 1
        allocate(dtmp(3_i4*ncfoc,3_i4),stat=status)
        if (status/=0) call outofmemory('deffreqd','dtmp')
        allocate(StatMPI(MPI_Status_Size),stat=status)
        if (status/=0) call outofmemory('deffreqd','StatMPI')
      endif
#endif
      igroup = mcv/3
      iresid = mcv - igroup*3
      indi = 0
      if (linten) then
        if (igroup.gt.0) then
          do i = 1,igroup
            if (ioproc) then
              write(ioout,'(''  Frequency   '',6f10.4)') (freq(indi+j,1),j=1,3)
              write(ioout,'(''  IR Intensity'',6f10.4)') (IRintensity(indi+j,1),j=1,3)
              write(ioout,'(''     in X     '',6f10.4)') (irx(indi+j),j=1,3)
              write(ioout,'(''     in Y     '',6f10.4)') (iry(indi+j),j=1,3)
              write(ioout,'(''     in Z     '',6f10.4,/)') (irz(indi+j),j=1,3)
            endif
            if (nprocs.gt.1) call mpbarrier
#ifdef MPI
            if (lioproconly) then
              iloc  = mcvrptr(indi+1)
              if (mcvnptr(i).ne.0_i4) then
!
!  Post receive
!
                if (ioproc) then
                  nnode = mcvnptr(i)
                  ntmp = 9*ncfoc
                  call MPI_IRecv(dtmp,ntmp,MPI_double_precision,nnode, &
                                 ntag,MPI_Comm_World,Request,MPIerror)
                endif
!
!  Pass data to ioproc for writing
!
                if (iloc.gt.0) then
                  do k = 1,3
                    dtmp(1:3*ncfoc,k) = vectors(1:3*ncfoc,mcvrptr(indi+k))
                  enddo
!
!  Post send
!
                  ntmp = 9*ncfoc
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
                  indj = 0
                  do j = 1,ncfoc
                    write(ioout,'(i6,'' x '',4x,6f10.6)') j,(dtmp(indj+1,k),k=1,3)
                    write(ioout,'(i6,'' y '',4x,6f10.6)') j,(dtmp(indj+2,k),k=1,3)
                    write(ioout,'(i6,'' z '',4x,6f10.6)') j,(dtmp(indj+3,k),k=1,3)
                    indj = indj + 3
                  enddo
                endif
              else
                if (iloc.gt.0) then
                  indj = 0
                  do j = 1,ncfoc
                    write(ioout,'(i6,'' x '',4x,6f10.6)') j,(vectors(indj+1,mcvrptr(indi+k)),k=1,3)
                    write(ioout,'(i6,'' y '',4x,6f10.6)') j,(vectors(indj+2,mcvrptr(indi+k)),k=1,3)
                    write(ioout,'(i6,'' z '',4x,6f10.6)') j,(vectors(indj+3,mcvrptr(indi+k)),k=1,3)
                    indj = indj + 3
                  enddo
                endif
              endif
            else
#endif
              if (mcvrptr(indi+1).gt.0) then
                indj = 0
                do j = 1,ncfoc
                  write(ioout,'(i6,'' x '',4x,6f10.6)') j,(vectors(indj+1,mcvrptr(indi+k)),k=1,3)
                  write(ioout,'(i6,'' y '',4x,6f10.6)') j,(vectors(indj+2,mcvrptr(indi+k)),k=1,3)
                  write(ioout,'(i6,'' z '',4x,6f10.6)') j,(vectors(indj+3,mcvrptr(indi+k)),k=1,3)
                  indj = indj + 3
                enddo
              endif
#ifdef MPI
            endif
#endif
            indi = indi + 3
            if (nprocs.gt.1) call mpbarrier
            if (ioproc) then
              write(ioout,'(/)')
            endif
          enddo
        endif
        if (iresid.gt.0) then
          if (ioproc) then
            write(ioout,'(''  Frequency   '',6f10.4)') (freq(indi+j,1),j=1,iresid)
            write(ioout,'(''  IR Intensity'',6f10.4)') (IRintensity(indi+j,1),j=1,iresid)
            write(ioout,'(''     in X     '',6f10.4)') (irx(indi+j),j=1,iresid)
            write(ioout,'(''     in Y     '',6f10.4)') (iry(indi+j),j=1,iresid)
            write(ioout,'(''     in Z     '',6f10.4)') (irz(indi+j),j=1,iresid)
            write(ioout,'(/)',advance='no')
          endif
          if (nprocs.gt.1) call mpbarrier
#ifdef MPI
          if (lioproconly) then
            iloc  = mcvrptr(indi+1)
            if (mcvnptr(i).ne.0_i4) then
!
!  Post receive
!
              if (ioproc) then
                nnode = mcvnptr(i)
                ntmp = 3*ncfoc*iresid
                call MPI_IRecv(dtmp,ntmp,MPI_double_precision,nnode, &
                               ntag,MPI_Comm_World,Request,MPIerror)
              endif
!
!  Pass data to ioproc for writing
!
              if (iloc.gt.0) then
                do k = 1,iresid
                  dtmp(1:3*ncfoc,k) = vectors(1:3*ncfoc,mcvrptr(indi+k))
                enddo
!
!  Post send
!
                ntmp = 3*ncfoc*iresid
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
                indj = 0
                do j = 1,ncfoc
                  write(ioout,'(i6,'' x '',4x,6f10.6)') j,(dtmp(indj+1,k),k=1,iresid)
                  write(ioout,'(i6,'' y '',4x,6f10.6)') j,(dtmp(indj+2,k),k=1,iresid)
                  write(ioout,'(i6,'' z '',4x,6f10.6)') j,(dtmp(indj+3,k),k=1,iresid)
                  indj = indj + 3
                enddo
              endif
            else
              if (iloc.gt.0) then
                indj = 0
                do j = 1,ncfoc
                  write(ioout,'(i6,'' x '',4x,6f10.6)') j,(vectors(indj+1,mcvrptr(indi+k)),k=1,iresid)
                  write(ioout,'(i6,'' y '',4x,6f10.6)') j,(vectors(indj+2,mcvrptr(indi+k)),k=1,iresid)
                  write(ioout,'(i6,'' z '',4x,6f10.6)') j,(vectors(indj+3,mcvrptr(indi+k)),k=1,iresid)
                  indj = indj + 3
                enddo
              endif
            endif
          else
#endif
            if (mcvrptr(indi+1).gt.0) then
              indj = 0
              do j = 1,ncfoc
                write(ioout,'(i6,'' x '',4x,6f10.6)') j,(vectors(indj+1,mcvrptr(indi+k)),k=1,iresid)
                write(ioout,'(i6,'' y '',4x,6f10.6)') j,(vectors(indj+2,mcvrptr(indi+k)),k=1,iresid)
                write(ioout,'(i6,'' z '',4x,6f10.6)') j,(vectors(indj+3,mcvrptr(indi+k)),k=1,iresid)
                indj = indj + 3
              enddo
            endif
#ifdef MPI
          endif
#endif
          if (nprocs.gt.1) call mpbarrier
        endif
        if (ioproc) then
          write(ioout,'(/)')
        endif
      else
        if (igroup.gt.0) then
          do i = 1,igroup
            if (ioproc) then
              write(ioout,'(''  Frequency   '',6f10.4,/)') (freq(indi+j,1),j = 1,3)
            endif
            if (nprocs.gt.1) call mpbarrier
#ifdef MPI
            if (lioproconly) then
              iloc  = mcvrptr(indi+1)
              if (mcvnptr(i).ne.0_i4) then
!
!  Post receive
!
                if (ioproc) then
                  nnode = mcvnptr(i)
                  ntmp = 9*ncfoc
                  call MPI_IRecv(dtmp,ntmp,MPI_double_precision,nnode, &
                                 ntag,MPI_Comm_World,Request,MPIerror)
                endif
!
!  Pass data to ioproc for writing
!
                if (iloc.gt.0) then
                  do k = 1,3
                    dtmp(1:3*ncfoc,k) = vectors(1:3*ncfoc,mcvrptr(indi+k))
                  enddo
!
!  Post send
!
                  ntmp = 9*ncfoc
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
                  indj = 0
                  do j = 1,ncfoc
                    write(ioout,'(i6,'' x '',4x,6f10.6)') j,(dtmp(indj+1,k),k=1,3)
                    write(ioout,'(i6,'' y '',4x,6f10.6)') j,(dtmp(indj+2,k),k=1,3)
                    write(ioout,'(i6,'' z '',4x,6f10.6)') j,(dtmp(indj+3,k),k=1,3)
                    indj = indj + 3
                  enddo
                endif
              else
                if (iloc.gt.0) then
                  indj = 0
                  do j = 1,ncfoc
                    write(ioout,'(i6,'' x '',4x,6f10.6)') j,(vectors(indj+1,mcvrptr(indi+k)),k=1,3)
                    write(ioout,'(i6,'' y '',4x,6f10.6)') j,(vectors(indj+2,mcvrptr(indi+k)),k=1,3)
                    write(ioout,'(i6,'' z '',4x,6f10.6)') j,(vectors(indj+3,mcvrptr(indi+k)),k=1,3)
                    indj = indj + 3
                  enddo
                endif
              endif
            else
#endif
              if (mcvrptr(indi+1).gt.0) then
                indj = 0
                do j = 1,ncfoc
                  write(ioout,'(i6,'' x '',4x,6f10.6)') j,(vectors(indj+1,mcvrptr(indi+k)),k=1,3)
                  write(ioout,'(i6,'' y '',4x,6f10.6)') j,(vectors(indj+2,mcvrptr(indi+k)),k=1,3)
                  write(ioout,'(i6,'' z '',4x,6f10.6)') j,(vectors(indj+3,mcvrptr(indi+k)),k=1,3)
                  indj = indj + 3
                enddo
              endif
#ifdef MPI
            endif
#endif
            indi = indi + 3
            if (nprocs.gt.1) call mpbarrier
            if (ioproc) then
              write(ioout,'(/)')
            endif
          enddo
        endif
        if (iresid.gt.0) then
          if (ioproc) then
            write(ioout,'(''  Frequency   '',6f10.4)') (freq(indi+j,1),j = 1,iresid)
            write(ioout,'(/)',advance='no')
          endif
          if (nprocs.gt.1) call mpbarrier
#ifdef MPI
          if (lioproconly) then
            iloc  = mcvrptr(indi+1)
            if (mcvnptr(i).ne.0_i4) then
!
!  Post receive
!
              if (ioproc) then
                nnode = mcvnptr(i)
                ntmp = 3*ncfoc*iresid
                call MPI_IRecv(dtmp,ntmp,MPI_double_precision,nnode, &
                               ntag,MPI_Comm_World,Request,MPIerror)
              endif
!
!  Pass data to ioproc for writing
!
              if (iloc.gt.0) then
                do k = 1,iresid
                  dtmp(1:3*ncfoc,k) = vectors(1:3*ncfoc,mcvrptr(indi+k))
                enddo
!
!  Post send
!
                ntmp = 3*ncfoc*iresid
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
                indj = 0
                do j = 1,ncfoc
                  write(ioout,'(i6,'' x '',4x,6f10.6)') j,(dtmp(indj+1,k),k=1,iresid)
                  write(ioout,'(i6,'' y '',4x,6f10.6)') j,(dtmp(indj+2,k),k=1,iresid)
                  write(ioout,'(i6,'' z '',4x,6f10.6)') j,(dtmp(indj+3,k),k=1,iresid)
                  indj = indj + 3
                enddo
              endif
            else
              if (iloc.gt.0) then
                indj = 0
                do j = 1,ncfoc
                  write(ioout,'(i6,'' x '',4x,6f10.6)') j,(vectors(indj+1,mcvrptr(indi+k)),k=1,iresid)
                  write(ioout,'(i6,'' y '',4x,6f10.6)') j,(vectors(indj+2,mcvrptr(indi+k)),k=1,iresid)
                  write(ioout,'(i6,'' z '',4x,6f10.6)') j,(vectors(indj+3,mcvrptr(indi+k)),k=1,iresid)
                  indj = indj + 3
                enddo
              endif
            endif
          else
#endif
            if (mcvrptr(indi+1).gt.0) then
              indj = 0
              do j = 1,ncfoc
                write(ioout,'(i6,'' x '',4x,6f10.6)') j,(vectors(indj+1,mcvrptr(indi+k)),k=1,iresid)
                write(ioout,'(i6,'' y '',4x,6f10.6)') j,(vectors(indj+2,mcvrptr(indi+k)),k=1,iresid)
                write(ioout,'(i6,'' z '',4x,6f10.6)') j,(vectors(indj+3,mcvrptr(indi+k)),k=1,iresid)
                indj = indj + 3
              enddo
            endif
#ifdef MPI
          endif
#endif
          if (nprocs.gt.1) call mpbarrier
        endif
        if (ioproc) then
          write(ioout,'(/)')
        endif
      endif
#ifdef MPI
      if (lioproconly) then
        deallocate(StatMPI,stat=status)
        if (status/=0) call deallocate_error('deffreqd','StatMPI')
        deallocate(dtmp,stat=status)
        if (status/=0) call deallocate_error('deffreqd','dtmp')
      endif
#endif
      write(ioout,'(/)')
    elseif (linten.and.ioproc) then
      write(ioout,'(/,''  Frequencies (cm-1) and IR Intensities : '',/)')
      igroup = mcv/6
      iresid = mcv - igroup*6
      indi = 0
      if (igroup.gt.0) then
        do i = 1,igroup
          write(ioout,'(''  Frequency   '',6f10.4)') (freq(indi+j,1),j = 1,6)
          write(ioout,'(''  IR Intensity'',6f10.4)') (IRintensity(indi+j,1),j = 1,6)
          write(ioout,'(''     in X     '',6f10.4)') (irx(indi+j),j = 1,6)
          write(ioout,'(''     in Y     '',6f10.4)') (iry(indi+j),j = 1,6)
          write(ioout,'(''     in Z     '',6f10.4,/)') (irz(indi+j),j = 1,6)
          indi = indi + 6
        enddo
      endif
      if (iresid.gt.0) then
        write(ioout,'(''  Frequency   '',6f10.4)') (freq(indi+j,1),j = 1,iresid)
        write(ioout,'(''  IR Intensity'',6f10.4)') (IRintensity(indi+j,1),j = 1,iresid)
        write(ioout,'(''     in X     '',6f10.4)') (irx(indi+j),j = 1,iresid)
        write(ioout,'(''     in Y     '',6f10.4)') (iry(indi+j),j = 1,iresid)
        write(ioout,'(''     in Z     '',6f10.4)') (irz(indi+j),j = 1,iresid)
        write(ioout,'(/)',advance='no')
      endif
      write(ioout,'(/)')
    elseif (lfreqout.and.lprinloc) then
      write(ioout,'(/,''  Frequencies (cm-1) : '',/)')
      write(ioout,'(9(f8.2))') (freq(j,1),j = 1,mcv)
      write(ioout,'(/)')
    endif
  else
    if (lfreqout.and.lprinloc) then
      write(ioout,'(/,''  Frequencies (cm-1) :'',/)')
      write(ioout,'(9f8.2)')(freq(i,1),i = 1,mcv)
      write(ioout,'(/)')
    endif
  endif
!
!  Output eigenvectors if requested
!
  if (lfreqout.and.lprinloc) then
    write(ioout,'(''--------------------------------------------------------------------------------'',/)')
  endif
!
!  Lower symmetry to remove imaginary modes if selected
!
  call lower(mcv,mcvrptr,freq(1,1),ncoreg1,iocptr,ncoreg1,iocptr,maxd2,vectors)
!
  mcvmin = 1
!
!  Check for imaginary modes and exclude them
!
  nimag = 0
  do i = 1,mcv
    if (freq(i,1).lt.-0.5_dp) nimag = nimag + 1
  enddo
  mcvmin = mcvmin + max(0,nimag-3)
  mcvmax = mcv
  if (minmode.ne.1) mcvmin = minmode
  if (maxmode.ne.0) mcvmax = maxmode
!*************************************
!  Output phonon related properties  *
!*************************************
  if (lprinloc) then
    if (ntemperatureramp.gt.0) then
      do ntr = 1,ntemperatureramp
!
!  For first ramp need to start from initial temperature
!  For subsequent ramps this would be the same as the end of the previous ramp
!
        if (ntr.eq.1) then
          nt0 = 0
        else
          nt0 = 1
        endif
        do nt = nt0,ntemperaturestep(ntr)
          tem = temperaturestart(ntr) + dble(nt)*temperaturestep(ntr)
!
!  Zero thermodynamic properties
!
          zpe = 0.0_dp
          entropy = 0.0_dp
          fhenergy = 0.0_dp
          rtlnz = 0.0_dp
          cv = 0.0_dp
          fe_equipartition = 0.0_dp
          s_equipartition = 0.0_dp
!***************************************
!  Evaluate phonon related properties  *
!***************************************
          if (tem.gt.1.0d-3) then
!
!  Scale frequencies to hw/kT
!
            rkt = boltz*tem
            cmfact = planck*speedl/rkt
            do i = 1,mcv
              rtmp2(i) = cmfact*freq(i,1)
            enddo
!
!  Store exp(x) in w1, exp(x)-1 in w2 and 1/(exp(x)-1) in w3
!
            do i = 1,mcv
              if (rtmp2(i).lt.12.0_dp) then
                w1(i) = exp(rtmp2(i))
                w2(i) = w1(i) - 1.0_dp
                if (abs(w2(i)).gt.0.0_dp) w3(i) = 1.0_dp/w2(i)
              else
                w3(i) = exp(-rtmp2(i))
              endif
            enddo
!
!  Zero point energy
!
            factor = 0.5_dp*rkt/evtoj
            trmzp = 0.0_dp
            if (.not.lnozero) then
              do i = mcvmin,mcvmax
                if (rtmp2(i).gt.cmfact) trmzp = trmzp + rtmp2(i)
              enddo
              trmzp = factor*trmzp
              zpe = zpe + trmzp
            endif
!
!  Entropy and free energy
!
            trmfe = 0.0_dp
            trmen = 0.0_dp
            freqmin = cmfact
            do i = mcvmin,mcvmax
              trm1 = rtmp2(i)
              if (trm1.gt.freqmin) then
                trm1 = 1.0_dp - exp(-trm1)
                trmfe = trmfe + log(trm1)
                trmen = trmen + rtmp2(i)*w3(i)
              endif
            enddo
!
!  Equipartition free energy and entropy
!
            trmfe_eq = 0.0_dp
            trms_eq = 0.0_dp
            rmode = 0.0_dp
            do i = 1,mcv
              trm1 = rtmp2(i)
              if (trm1.gt.freqmin) then
                trmfe_eq = trmfe_eq + log(trm1)
                trms_eq  = trms_eq + log(trm1) - 1.0_dp
                rmode = rmode + 1.0_dp
              endif
            enddo
!
            factor = 2.0_dp*factor
            trm = factor*trmfe
            fhenergy = fhenergy + trm + trmzp
            rtlnz = rtlnz - trm
            fe_equipartition = fe_equipartition + factor*trmfe_eq
            s_equipartition = s_equipartition + factor*trms_eq
            factor = factor/tem
            entropy = entropy + factor*trmen
!
!  Heat capacity - constant volume
!
            trmcv = 0.0_dp
            do i = mcvmin,mcvmax
              if (rtmp2(i).gt.freqmin) then
                if (rtmp2(i).lt.12.0_dp) then
                  trmcv = trmcv + rtmp2(i)*rtmp2(i)*w1(i)*w3(i)*w3(i)
                else
                  trmcv = trmcv + rtmp2(i)*rtmp2(i)*w3(i)
                endif
              endif
            enddo
            cv = cv + factor*trmcv
          else
!
!  Zero point energy
!
            factor = 0.5_dp*planck*speedl/evtoj
            if (.not.lnozero) then
              trmzp = 0.0_dp
              do i = mcvmin,mcvmax
                if (freq(i,1).gt.1.0_dp) trmzp = trmzp + freq(i,1)
              enddo
              trmzp = factor*trmzp
              zpe = zpe + trmzp
            endif
          endif
          write(ioout,'(''  Vibrational properties (for region 1): '',''Temperature  =  '',f10.3,'' K'')')tem
          write(ioout,'(''--------------------------------------------------------------------------------'')')
          write(ioout,'(''  Zero point energy             =  '',f15.6,'' eV'')')zpe
          if (tem.gt.1.0d-03) then
            trmen = fhenergy - zpe
            entropy = entropy - trmen/tem
            ent2 = entropy*evtoj*avogadro
            cv2 = cv*evtoj*avogadro
            write(ioout,'(''  Entropy                       =  '',f15.6,'' eV/K'')') entropy
            write(ioout,'(''                                =  '',f15.6,'' J/(mol.K)'')') ent2
            write(ioout,'(''  Helmholtz free-energy         =  '',f15.6,'' eV'')') fhenergy + fc
            write(ioout,'(''                                =  '',f15.6,'' kJmol-1'')') (fhenergy+fc)*evtoj*avogadro*0.001_dp
            write(ioout,'(''  Heat capacity - const volume  =  '',f15.6,'' eV/K'')') cv
            write(ioout,'(''                                =  '',f15.6,'' J/(mol.K)'')') cv2
          endif
          write(ioout,'(''--------------------------------------------------------------------------------'')')
        enddo
      enddo
    else
!
!  Single temperature
!
      tem = temperature
!
!  Zero thermodynamic properties
!
      zpe = 0.0_dp
      entropy = 0.0_dp
      fhenergy = 0.0_dp
      rtlnz = 0.0_dp
      cv = 0.0_dp
      fe_equipartition = 0.0_dp
      s_equipartition = 0.0_dp
!***************************************
!  Evaluate phonon related properties  *
!***************************************
      if (tem.gt.1.0d-3) then
!
!  Scale frequencies to hw/kT
!
        rkt = boltz*tem
        cmfact = planck*speedl/rkt
        do i = 1,mcv
          rtmp2(i) = cmfact*freq(i,1)
        enddo
!
!  Store exp(x) in w1, exp(x)-1 in w2 and 1/(exp(x)-1) in w3
!
        do i = 1,mcv
          if (rtmp2(i).lt.12.0_dp) then
            w1(i) = exp(rtmp2(i))
            w2(i) = w1(i) - 1.0_dp
            if (abs(w2(i)).gt.0.0_dp) w3(i) = 1.0_dp/w2(i)
          else
            w3(i) = exp(-rtmp2(i))
          endif
        enddo
!
!  Zero point energy
!
        factor = 0.5_dp*rkt/evtoj
        trmzp = 0.0_dp
        if (.not.lnozero) then
          do i = mcvmin,mcvmax
            if (rtmp2(i).gt.cmfact) trmzp = trmzp + rtmp2(i)
          enddo
          trmzp = factor*trmzp
          zpe = zpe + trmzp
        endif
!
!  Entropy and free energy
!
        trmfe = 0.0_dp
        trmen = 0.0_dp
        freqmin = cmfact
        do i = mcvmin,mcvmax
          trm1 = rtmp2(i)
          if (trm1.gt.freqmin) then
            trm1 = 1.0_dp - exp(-trm1)
            trmfe = trmfe + log(trm1)
            trmen = trmen + rtmp2(i)*w3(i)
          endif
        enddo
!
!  Equipartition free energy and entropy
!
        trmfe_eq = 0.0_dp
        trms_eq = 0.0_dp
        rmode = 0.0_dp
        do i = 1,mcv
          trm1 = rtmp2(i)
          if (trm1.gt.freqmin) then
            trmfe_eq = trmfe_eq + log(trm1)
            trms_eq  = trms_eq + log(trm1) - 1.0_dp
            rmode = rmode + 1.0_dp
          endif
        enddo
!
        factor = 2.0_dp*factor
        trm = factor*trmfe
        fhenergy = fhenergy + trm + trmzp
        rtlnz = rtlnz - trm
        fe_equipartition = fe_equipartition + factor*trmfe_eq
        s_equipartition = s_equipartition + factor*trms_eq
        factor = factor/tem
        entropy = entropy + factor*trmen
!
!  Heat capacity - constant volume
!
        trmcv = 0.0_dp
        do i = mcvmin,mcvmax
          if (rtmp2(i).gt.freqmin) then
            if (rtmp2(i).lt.12.0_dp) then
              trmcv = trmcv + rtmp2(i)*rtmp2(i)*w1(i)*w3(i)*w3(i)
            else
              trmcv = trmcv + rtmp2(i)*rtmp2(i)*w3(i)
            endif
          endif
        enddo
        cv = cv + factor*trmcv
      else
!
!  Zero point energy
!
        factor = 0.5_dp*planck*speedl/evtoj
        if (.not.lnozero) then
          trmzp = 0.0_dp
          do i = mcvmin,mcvmax
            if (freq(i,1).gt.1.0_dp) trmzp = trmzp + freq(i,1)
          enddo
          trmzp = factor*trmzp
          zpe = zpe + trmzp
        endif
      endif
      write(ioout,'(''  Vibrational properties (for region 1): '',''Temperature  =  '',f10.3,'' K'')')tem
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''  Zero point energy             =  '',f15.6,'' eV'')')zpe
      if (tem.gt.1.0d-03) then
        trmen = fhenergy - zpe
        entropy = entropy - trmen/tem
        ent2 = entropy*evtoj*avogadro
        cv2 = cv*evtoj*avogadro
        write(ioout,'(''  Entropy                       =  '',f15.6,'' eV/K'')') entropy
        write(ioout,'(''                                =  '',f15.6,'' J/(mol.K)'')') ent2
        write(ioout,'(''  Helmholtz free-energy         =  '',f15.6,'' eV'')') fhenergy + fc
        write(ioout,'(''                                =  '',f15.6,'' kJmol-1'')') (fhenergy+fc)*evtoj*avogadro*0.001_dp
        write(ioout,'(''  Heat capacity - const volume  =  '',f15.6,'' eV/K'')') cv
        write(ioout,'(''                                =  '',f15.6,'' J/(mol.K)'')') cv2
      endif
      write(ioout,'(''--------------------------------------------------------------------------------'')')
    endif
  endif
!
!  If output if required then call outfreq
!
  if (lprinloc) call outfreq(leigloc,1_i4,freq,IRintensity,ncoreg1)
!
!  Free local memory
!
  deallocate(rmassdef,stat=status)
  if (status/=0) call deallocate_error('deffreqd','rmassdef')
  deallocate(massdef,stat=status)
  if (status/=0) call deallocate_error('deffreqd','massdef')
  deallocate(w3,stat=status)
  if (status/=0) call deallocate_error('deffreqd','w3')
  deallocate(w2,stat=status)
  if (status/=0) call deallocate_error('deffreqd','w2')
  deallocate(w1,stat=status)
  if (status/=0) call deallocate_error('deffreqd','w1')
  deallocate(rtmp2,stat=status)
  if (status/=0) call deallocate_error('deffreqd','rtmp2')
  deallocate(irz,stat=status)
  if (status/=0) call deallocate_error('deffreqd','irz')
  deallocate(iry,stat=status)
  if (status/=0) call deallocate_error('deffreqd','iry')
  deallocate(irx,stat=status)
  if (status/=0) call deallocate_error('deffreqd','irx')
  deallocate(mcvrptr,stat=status)
  if (status/=0) call deallocate_error('deffreqd','mcvrptr')
  deallocate(mcvnptr,stat=status)
  if (status/=0) call deallocate_error('deffreqd','mcvnptr')
  deallocate(mcvptr,stat=status)
  if (status/=0) call deallocate_error('deffreqd','mcvptr')
  deallocate(iocptr,stat=status)
  if (status/=0) call deallocate_error('deffreqd','iocptr')
  deallocate(ibocptr,stat=status)
  if (status/=0) call deallocate_error('deffreqd','ibocptr')
  deallocate(ndbsptr,stat=status)
  if (status/=0) call deallocate_error('deffreqd','ndbsptr')
!
!  Timing
!
  t2t = g_cpu_time()
  tphon = t2t - t1t + tphon
#ifdef TRACE
  call trace_out('deffreqd')
#endif
#else
  call outerror('deffreqd called when not compiled with MPI',0_i4)
  call stopnow('deffreqd')
#endif
!
  return
  end
