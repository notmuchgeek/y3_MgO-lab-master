  subroutine phonond(lprint,fc,nobsmodeptr0,nobsmode)
!
!  Calculates the phonons at a given set of k points.
!  Distributed second derivative version. 
!
!  NB: Currently doesn't support partial occupancy!
!
!  The block diagonal elements of derv2 must be stored, as these are
!  not recalculated in dynamic as the summing of off diagonals is no
!  longer applicable.
!
!  nkpt    = total number of k points across all structures
!  nllkpt  = number of k points for this structure
!  nlkpt   = pointer to lowest k point
!  nukpt   = pointer to upper k point
!  xkpt    = fractional x component of k point
!  ykpt    = fractional y component of k point
!  zkpt    = fractional z component of k point
!  wkpt    = weight of each k point
!  nkptcfg = configuration pointer for each k point
!  sumwkpt = sum over weights of k points
!  leigloc = local flag to indicate whether eigenvectors are 
!            to be generated for this configuration
!
!  11/16 Created from phonon
!  12/16 Merge cluster frequency functionality into this routine
!   5/17 IO handling in parallel added for dynamical matrix
!   5/17 Calls to matrix inversion modified to use nmin.ne.1
!   5/17 Parallel thermal conductivity call added
!   5/17 Handling of output of dynamical matrix added
!   5/17 nshell and nshellonnode added to subroutine call
!   6/17 Third derivatives output for shengBTE added back
!   7/17 nshell and nshellonnode added to complex subroutine call
!   7/17 outfrc now called in parallel
!   7/17 Call to getvibmode replaced by getvibmoded
!   8/17 Routine wrapped in ifdef since it won't work or be called without MPI
!   8/17 Third order derivatives removed for ShengBTE since this
!        is not yet supported in parallel
!   8/17 Correction for region 2 phonon case
!   8/17 Modifications to accommodate Intel MPI and scalapack
!  11/17 leigloc set to true for MSD calc
!  11/17 wk now passed to peigend
!   1/18 Trace added
!   2/18 Grueneisen parameters added
!   2/18 Set up of fourbody list added for Grueneisen calculation
!   3/18 dynam shortened to dyna in keyword check
!   3/18 Parallel I/O corrected
!   4/18 leigloc now set here for llower and linten
!   3/19 Multiple temperature ramps added
!   5/19 Finite difference flag split for first and second derivatives
!   8/19 Trap added for parallel error when there is no work on a node
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!   4/20 maxfqat changed to maxfreq since it is no longer just atom based
!   4/20 Mass arrays now use fmass and rfmass
!   5/20 ncfoc removed from arguments to groupvelocities
!   6/20 Rigid molecule modifications added
!   6/20 Correction to integerKpoint logic
!   6/20 Calls to lower removed except for gamma point
!   6/20 Last argument to outphon changed
!   6/20 mtvptr passed to groupvelocitiesd
!   6/20 Criteria for gamma point made consistent
!   6/20 nmolcore changes added
!   6/20 Order of operations changed - shells now handled before
!        generating molecular terms so that eckart works with shells
!   6/20 Correction to shell contribution to group velocities
!   6/20 Wrapping of communication of group velocities with not lgamma added
!   7/20 Correction to maxeigc for case with shells
!   7/20 Correction to second dimension for transfer back to derv2/dervi when
!        removing shell contribution
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
  use configurations
  use g_constants
  use control
  use current
  use derivatives
  use dispersion
  use element
  use gulp_files
#ifdef MPI
  use four,            only : nfor
#endif
  use frequencies
  use general
#ifdef MPI
  use gulp_cml,        only : lcml
  use gulp_cml_phonon, only : gulp_cml_startPhonons, gulp_cml_outPDF
  use gulp_cml_phonon, only : gulp_cml_addKPoint, gulp_cml_endPhonons, gulp_cml_PDFstats
#endif
  use iochannels
  use ksample
  use ksample_scatter
#ifdef MPI
  use m_pdf,           only : closepdfphonon, pdfsetup
  use m_pdfneutron
  use molecule
  use observables,     only : fobsmodefreq, fobsmodeover
#endif
  use parallel
  use partial
  use phononatoms
  use projectdos
  use properties
#ifdef MPI
  use scatterdata,     only : lscattercall
  use species,         only : massspec, natspec, ntypspec
#endif
  use shells
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
  integer(i4),  intent(in)                          :: nobsmodeptr0  ! Pointer to first observable mode - 1
  integer(i4),  intent(in)                          :: nobsmode      ! Number of observable modes if fitting run
  logical,      intent(in)                          :: lprint        ! If true then output results
  real(dp),     intent(in)                          :: fc            ! Internal energy
#ifdef MPI
!
!  Local variables
!
  character(len=5)                                  :: lab1
  character(len=2)                                  :: fstring1
  character(len=2)                                  :: fstring2
  character(len=12)                                 :: fstring
  integer(i4)                                       :: i
  integer(i4)                                       :: iatm
  integer(i4)                                       :: ifail
  integer(i4)                                       :: ii
  integer(i4)                                       :: icount
  integer(i4)                                       :: iloc
  integer(i4)                                       :: indi
  integer(i4)                                       :: indif
  integer(i4)                                       :: indii
  integer(i4)                                       :: indil
  integer(i4)                                       :: indiloc
  integer(i4)                                       :: indj
  integer(i4)                                       :: indjf
  integer(i4)                                       :: indjj
  integer(i4)                                       :: indjl
  integer(i4)                                       :: ix
  integer(i4)                                       :: iy
  integer(i4)                                       :: iz
  integer(i4)                                       :: j
  integer(i4)                                       :: jj
  integer(i4)                                       :: jx
  integer(i4)                                       :: k
  integer(i4)                                       :: m
  integer(i4)                                       :: maxeigc
  integer(i4)                                       :: maxeigr
  integer(i4)                                       :: maxlim
  integer(i4)                                       :: mcnmv
  integer(i4)                                       :: mcv
  integer(i4)                                       :: mcvloc
  integer(i4),  dimension(:),     allocatable       :: mcvptr
  integer(i4)                                       :: mind
  integer(i4)                                       :: mindi
  integer(i4)                                       :: mindj
  integer(i4),  dimension(:),     allocatable       :: mtvptr
  integer(i4),  dimension(:),     allocatable       :: mtvnptr
  integer(i4),  dimension(:),     allocatable       :: mtvrptr
  integer(i4)                                       :: mint
  integer(i4)                                       :: mintloc
  integer(i4)                                       :: mis
  integer(i4)                                       :: mjs
  integer(i4)                                       :: mm
  integer(i4)                                       :: mmv
  integer(i4)                                       :: mmvloc
  integer(i4)                                       :: msv
  integer(i4)                                       :: msvloc
  integer(i4)                                       :: mtsv        ! Number of vibrations plus shell coordinates = mtv + msv
  integer(i4)                                       :: mtsvloc     ! Number of vibrations plus shell coordinates on local node
  integer(i4)                                       :: mtv         ! Number of vibrations = mcv + mvv
  integer(i4)                                       :: mtvloc      ! Number of vibrations on local node
  integer(i4)                                       :: nd
  integer(i4)                                       :: nf
  integer(i4)                                       :: nfitmode
  integer(i4)                                       :: nk
  integer(i4)                                       :: nldpt
  integer(i4)                                       :: nlkpt
  integer(i4)                                       :: nllkpt
  integer(i4)                                       :: nobm
  integer(i4)                                       :: node
  integer(i4)                                       :: np
  integer(i4)                                       :: np_nlkpt
  integer(i4)                                       :: np_procs
  integer(i4)                                       :: nri
  integer(i4)                                       :: nrj
  integer(i4)                                       :: nsi
  integer(i4)                                       :: nt
  integer(i4)                                       :: ntmax
  integer(i4)                                       :: nudpt
  integer(i4)                                       :: nukpt
  integer(i4)                                       :: status
!
  integer                                           :: MPIerror
  integer                                           :: idesc(9)
  integer                                           :: ifails
  integer                                           :: ld
  integer                                           :: nb
  integer                                           :: ncs
  integer                                           :: nsize
  integer                                           :: nsizec
  integer                                           :: ntag
  integer                                           :: nnode
  integer                                           :: ntmp
  integer                                           :: Request
  integer                                           :: Requestc
  integer                                           :: rnode
  integer                                           :: snode
  integer(i4),      dimension(:),   allocatable     :: StatMPI       ! Array for status from MPI
!
#ifdef INTEL
  logical                                           :: leigcOK
#endif
  logical                                           :: lcmlloc
  logical                                           :: lcluster
  logical                                           :: ld2loc
  logical                                           :: ldiloc
  logical                                           :: ldendisp
  logical                                           :: leigloc
  logical                                           :: lfound
  logical                                           :: lgamma
  logical                                           :: lgammaonly
  logical                                           :: lintegerKpoint
  logical                                           :: lintegerKpoint_present
  logical                                           :: lnonanal
  logical                                           :: lnoanald2loc
  logical                                           :: lnoanald3loc
  logical                                           :: lpartofdisp
  logical                                           :: lpdfloc
  logical                                           :: lprinloc
  logical                                           :: ltemprop
  complex(dpc), dimension(:,:,:), allocatable       :: ctmp
  complex(dpc), dimension(:,:),   allocatable       :: eigc
  real(dp),     dimension(:),     allocatable       :: averfreq
  real(dp)                                          :: bornkloc(3)
  real(dp)                                          :: cmfact
  real(dp),     dimension(:,:),   allocatable       :: dtmp
  real(dp)                                          :: g_cpu_time
  real(dp)                                          :: dx
  real(dp)                                          :: dy
  real(dp)                                          :: dz
  real(dp),     dimension(:,:),   allocatable       :: eigr
  real(dp)                                          :: factor
  real(dp)                                          :: frq
  real(dp)                                          :: fscale
  real(dp)                                          :: mI(3)
  real(dp),     dimension(:),     allocatable       :: meanKEperatom
  real(dp),     dimension(:,:,:), allocatable       :: oscstrength
  real(dp)                                          :: overlap
  real(dp)                                          :: phi
  real(dp)                                          :: phistep
  real(dp)                                          :: projectionfactor
  real(dp)                                          :: qfrac(3)
  real(dp)                                          :: rkt
  real(dp),     dimension(:,:,:), allocatable       :: ramstrength
  real(dp)                                          :: rmassi
  real(dp)                                          :: rnokpt
  real(dp),     dimension(:),     allocatable       :: savefreq
  real(dp),     dimension(:,:),   allocatable       :: savederv2
  real(dp),     dimension(:,:),   allocatable       :: savedervi
  real(dp),     dimension(:,:),   allocatable       :: Sij
  real(dp),     dimension(:),     allocatable       :: sum
  real(dp),     dimension(:),     allocatable       :: sum2
  real(dp)                                          :: sumwkpt
  real(dp)                                          :: t1
  real(dp)                                          :: t1i
  real(dp)                                          :: t1t
  real(dp)                                          :: t2
  real(dp)                                          :: t2i
  real(dp)                                          :: t2t
  real(dp)                                          :: theta
  real(dp)                                          :: thetastep
  real(dp)                                          :: trmke
  real(dp)                                          :: w1
  real(dp)                                          :: w2
  real(dp)                                          :: w3
  real(dp)                                          :: weightpt
  real(dp)                                          :: wk
  real(dp)                                          :: xkt
  real(dp)                                          :: ykt
  real(dp)                                          :: zkt
  real(dp)                                          :: xmod
  real(dp)                                          :: ymod
  real(dp)                                          :: zmod
  real(dp)                                          :: wrk(9)
#ifdef TRACE
  call trace_in('phonond')
#endif
!
  t1t = g_cpu_time()
!
!  Check for things that aren't implemented yet
!
  if (nbsmat.gt.0) then
    call outerror('distributed breathing shell 2nd derivs not done',0_i4)
    call stopnow('phonond')
  endif
!
!  Set logicals
!
  leigloc = (leigen.or.leig.or.lcas)
  if ((nprojcfg(ncf)-nprojdef(ncf)).gt.0) leigloc = .true.
  if (.not.lprint) leigloc = .false.
  if (lomega(ncf)) leigloc = .true.
  if (lraman) leigloc = .true.
  if (nbornstep(ncf).gt.0) leigloc = .true.
  if (lmeanke) leigloc = .true.
  if (nobsmode.gt.0) leigloc = .true.
  if (lmsd) leigloc = .true.
  if (lgrueneisen) leigloc = .true.
  if (lthermal) leigloc = .true.
  if (lgroupvelocity) leigloc = .true.
  if (linten) leigloc = .true.
  if (llower) leigloc = .true.
  lcluster = (ndim.eq.0)
  lpdfloc = (lpdf.and..not.lcluster)
  lprinloc = (lprint.and.ioproc)
  ltemprop = (temperature.gt.1.0d-6.or.ntemperatureramp.gt.0)
  lnonanal = (index(keyword,'nono').eq.0.and.ndim.eq.3.and..not.leem.and..not.lnoanald2)
  lnoanald2loc = lnoanald2
  lnoanald3loc = lnoanald3
!
!  Is keyword present to force numerical third order force constants?
!
  if (index(keyword,'num3').eq.1.or.index(keyword,' num3').ne.0) lnoanald3loc = .true.
!
!  For a thermal conductivity calculation we need to turn off the non-analytic correction
!  to avoid the force constant matrix being overwritten.
!
  if (lthermal) lnonanal = .false.
!
!  CML phonon output
!
  lcmlloc = (lcml.and.lfreqout)
  if (lcmlloc) call gulp_cml_startPhonons
!
!  Set constants
!
!  fscale is the conversion factor from (eV/Ang^2)(mol/g) in the dynamical matrix to
!  wavenumbers for the frequencies. The factor of 10^23 comes from Ang^2 -> m^2 and
!  g -> kg.
!
  fscale = sqrt(1.0d23*evtoj*avogadro)
  fscale = fscale/(2.0_dp*pi*speedl)
!
!  Set flag as to whether PDOS or dispersion curves
!  are to be produced
!
  if ((lphono.or.ndline.gt.0.or.index(keyword,'node').eq.0)) then
    ldendisp = lprint
  else
    ldendisp = .false.
  endif
!
!  Work out number of k points for this structure.
!
!  Requires k points to be sorted such that all points
!  relating a given structure are consecutive.
!  
  nlkpt = 0
  do i = 1,nkpt
    nk = nkptcfg(i)
    if (nlkpt.eq.0.and.nk.eq.ncf) nlkpt = i
    if (nk.eq.ncf) nukpt = i
  enddo
!
!  No k points found for current structure, so return
!
  if (nlkpt.eq.0) then
#ifdef TRACE
    call trace_out('phonond')
#endif
    return
  endif
  nllkpt = nukpt - nlkpt + 1
!     
!  Work out lower and upper dispersion lines for this configuration
!     
  nldpt = 0
  if (ndline.gt.0) then
    do i = 1,ndline 
      nd = ndispcfg(i)
      if (nldpt.eq.0.and.nd.eq.ncf) nldpt = i
      if (nd.eq.ncf) nudpt = i
    enddo
  else
    nudpt = -1
  endif
!
!  If gamma approach direction hasn't been input and there is no dispersion then turn off non-analytic correction
!
  if ((nldpt-nudpt).eq.1.and..not.lbornkin(ncf)) lnonanal = .false.
!************************************
!  Build potential lists if needed  *
!************************************
  if (lgrueneisen) then
    if (nfor.gt.0) call setlist4
  endif
!************************************************
!  Find number of atoms for phonon calculation  *
!************************************************
!
!  Set phonon pointers 
!
  call setphonptr
!
!  Check whether partial occupancy is present that will fail with this routine
!
  if (lrigid) then
    lpocc = (ncore+nshell.ne.numat)
  else
    lpocc = (nsfoc+ncfoc.ne.nphonat)
  endif
  if (lpocc) then
    call outerror('partial occupancy not compatible with distributed 2nd derivs',0_i4)
    call stopnow('phonond')
  endif
!
!  Allocate local pointer arrays
!
  allocate(meanKEperatom(numat),stat=status)
  if (status/=0) call outofmemory('phonond','meanKEperatom')
!
!  Calculate a few constants to do with the size of the problem
!
  mint = 3*nphonat
  mintloc = 3*nphonatonnode
!
  maxlim = mint
  if (nbsmat.gt.0) maxlim = maxlim + nphonat
  if (lrigid) then
!
!  Rigid molecules - no partial occupancy
!
    msv = 3*nshell + nbsmatnomol
    mcv = 3*ncore
    mcvloc = 3*ncoreonnode
    msvloc = 3*nshellonnode
    mcnmv = 3*ncorenomol
    mmv = 6*nphonatm
    mmvloc = 6*nphonatonnodem
  else
!
!  General case
!
    msv = 3*nphonats + nphonatb
    mcv = 3*nphonatc
    mcnmv = 3*nphonatc
    msvloc = 3*nphonatonnodes + nphonatonnodeb
    mcvloc = 3*nphonatonnodec
    mmv = 0
    mmvloc = 0
    mtvloc = mcvloc
  endif
!
  mtv = mcnmv + mmv
  mtsv = mtv + msv
  mtsvloc = mcvloc + msvloc
!
  allocate(mtvptr(mtsv),stat=status)
  if (status/=0) call outofmemory('phonond','mtvptr')
  allocate(mtvnptr(mtsv),stat=status)
  if (status/=0) call outofmemory('phonond','mtvnptr')
  allocate(mtvrptr(mtsv),stat=status)
  if (status/=0) call outofmemory('phonond','mtvrptr')
!
  if (lrigid) then
!
!  Compute distribution of atoms and molecules for rigid molecule case
!
!  => Block cyclic distribution with numatnomol followed by nmol in block multiples of 3
!
    node = 0
    mtvloc = 0
    mtsvloc = 0
    icount = 0
    do i = 1,mtsv
      icount = icount + 1
      mtvnptr(i) = node
      if (node.eq.procid) then
        mtsvloc = mtsvloc + 1
        if (i.le.mtv) then
!
!  Only increment mtvloc until counter for mtsv reaches mtv
!
          mtvloc = mtvloc + 1
        endif
        mtvptr(mtsvloc) = i
        mtvrptr(i) = mtsvloc
      else
        mtvrptr(i) = 0
      endif
      if (icount.eq.3_i4*nblocksize) then
        icount = 0
        node = node + 1
        if (node.eq.nprocs) node = 0
      endif
    enddo
  else
    mtvrptr(1:mtsv) = 0
    do iloc = 1,nphonatonnodec
      i = nphonatonnodecptr(iloc)
      m = 3*(iloc-1)
      mm = 3*(i-1)
      mtvptr(m+1) = mm + 1
      mtvptr(m+2) = mm + 2
      mtvptr(m+3) = mm + 3
      mtvrptr(mm+1) = m + 1
      mtvrptr(mm+2) = m + 2
      mtvrptr(mm+3) = m + 3
    enddo
!
!  Set up pointer to node for each value of mcvrptr
!
    do i = 1,nphonatc
      mm = 3*(nphonatcptr(i)-1)
      mtvnptr(mm+1) = atom2node(nphonatcptr(i))
      mtvnptr(mm+2) = atom2node(nphonatcptr(i))
      mtvnptr(mm+3) = atom2node(nphonatcptr(i))
    enddo
  endif
!
!  Check that maxd2 is greater than or equal to mcv
!
  if (maxd2.lt.mtv) then
    maxd2 = mtv
    call changemaxd2
  endif
  if (maxd2u.lt.mtvloc) then
    maxd2u = mtvloc
    call changemaxd2
  endif
!
  if (leckart) then
!
!  Set up pointer to pass to eckart
!
    allocate(mcvptr(mcvloc),stat=status)
    if (status/=0) call outofmemory('phonon','mcvptr')
!
    do iloc = 1,nphonatonnodec
      i = nphonatonnodecptr(iloc)
      m = 3*(iloc-1)
      mm = 3*(i-1)
      mcvptr(m+1) = mm + 1
      mcvptr(m+2) = mm + 2
      mcvptr(m+3) = mm + 3
    enddo
  endif
!
!  Ensure that frequency array is large enough
!
  call changemaxfreq(mtv)
!
!  Allocate array to hold oscillator strengths if needed
!
  if (leigloc) then
    allocate(savefreq(mtv),stat=status)
    if (status/=0) call outofmemory('phonond','savefreq')
    allocate(oscstrength(3,3,mtv),stat=status)
    if (status/=0) call outofmemory('phonond','oscstrength')
    allocate(ramstrength(3,3,mtv),stat=status)
    if (status/=0) call outofmemory('phonond','ramstrength')
  endif
!
!  Allocate array for thermal conductivity if needed
!
  if (lthermal) then
    allocate(Sij(mtv,mtvloc),stat=status)
    if (status/=0) call outofmemory('phonond','Sij')
  endif
!
  rkt = boltz*temperature
!
!  If mean KE per atom is requested then initialise array to zero
!
  if (lmeanke) then
    meanKEperatom(1:nphonatc) = 0.0_dp
  endif
!
!  Calculate inversion square root of masses for atoms and momemt of inertia for quaternions
!
!  Now modified to handle partial occupancies
!
  do i = 1,mtv
    fmass(i) = 0.0_dp
  enddo
  do i = 1,nphonatc
    ii = 3*(i - 1)
    nsi = nspecptr(nrelf2a(nphonatcptr(i)))
    rmassi = massspec(nsi)*occuf(nphonatcptr(i))
    if (natspec(nsi).le.maxele.and.abs(rmassi).lt.1.0d-12) then
      call label(natspec(nsi),ntypspec(nsi),lab1)
      call outerror('mass of species '//lab1//' is zero',0_i4)
      goto 999
    endif
    fmass(ii+1) = fmass(ii+1) + rmassi
    fmass(ii+2) = fmass(ii+2) + rmassi
    fmass(ii+3) = fmass(ii+3) + rmassi
  enddo
  if (lrigid) then
    do i = 1,nmol
      ii = mcnmv + 6*(i-1)
!
      molaxes(1:3,1:3,i) = 0.0_dp
!
!  Sum atom masses for rigid molecule translations and recompute moment of inertia orientation
!
      do j = 1,nmolcore(i)
        k = nmollist(nmolptr(i)+j)
        nsi = nspecptr(nrelf2a(k))
        rmassi = massspec(nsi)*occuf(k)
        fmass(ii+1) = fmass(ii+1) + rmassi
        fmass(ii+2) = fmass(ii+2) + rmassi
        fmass(ii+3) = fmass(ii+3) + rmassi
!
        dx = molxyz(1,j,i)
        dy = molxyz(2,j,i)
        dz = molxyz(3,j,i)
!
        molaxes(1,1,i) = molaxes(1,1,i) + rmassi*dy*dy + rmassi*dz*dz
        molaxes(2,1,i) = molaxes(2,1,i) - rmassi*dx*dy
        molaxes(3,1,i) = molaxes(3,1,i) - rmassi*dx*dz
        molaxes(1,2,i) = molaxes(1,2,i) - rmassi*dy*dx
        molaxes(2,2,i) = molaxes(2,2,i) + rmassi*dx*dx + rmassi*dz*dz
        molaxes(3,2,i) = molaxes(3,2,i) - rmassi*dy*dz
        molaxes(1,3,i) = molaxes(1,3,i) - rmassi*dz*dx
        molaxes(2,3,i) = molaxes(2,3,i) - rmassi*dz*dy
        molaxes(3,3,i) = molaxes(3,3,i) + rmassi*dx*dx + rmassi*dy*dy
      enddo
!
!  Diagonalise tensor to get axes for moment of inertia tensor
!
      call dsyev('V','U',3_i4,molaxes(1,1,i),3_i4,mI,wrk,9_i4,ifail)
!
!  Moment of inertia for rotations
!
      fmass(ii+4) = mI(1)
      fmass(ii+5) = mI(2)
      fmass(ii+6) = mI(3)
    enddo
  endif
  do i = 1,mtv
    if (fmass(i).lt.too_small) then
      if (i.le.mcnmv) then
        call outerror('site has total mass close to zero in phonon',0_i4)
        call stopnow('phonond')
      elseif (mod(i-mcnmv-1,6_i4).lt.3) then
        call outerror('molecule has total mass close to zero in phonon',0_i4)
        call stopnow('phonond')
      else
!
!  Linear molecule - set mass to be large so that frequency is close to zero
!
        fmass(i) = 10000000.0_dp
      endif
    endif
    rfmass(i) = 1.0_dp/sqrt(fmass(i))
  enddo
  if (.not.lnoanald2loc.and..not.lfinitediff2.and..not.lcluster) then
!
!  Store diagonal blocks in derv3 to avoid recalculation
!
    do i = 1,nphonat
      iloc = atom2local(nphonatptr(i))
      if (iloc.gt.0) then
!
!  Atom is local to this node
!
        indi = 3*(nphonatptr(i)-1)
        indiloc = 3*(iloc-1)
        derv3(indiloc+1,1) = derv2(indi+1,indiloc+1)
        derv3(indiloc+2,1) = derv2(indi+2,indiloc+1)
        derv3(indiloc+3,1) = derv2(indi+3,indiloc+1)
        derv3(indiloc+1,2) = derv2(indi+1,indiloc+2)
        derv3(indiloc+2,2) = derv2(indi+2,indiloc+2)
        derv3(indiloc+3,2) = derv2(indi+3,indiloc+2)
        derv3(indiloc+1,3) = derv2(indi+1,indiloc+3)
        derv3(indiloc+2,3) = derv2(indi+2,indiloc+3)
        derv3(indiloc+3,3) = derv2(indi+3,indiloc+3)
        if (nbsmat.gt.0) then
          indi = mint + i
          indiloc = mintloc + iloc
          derv3(indiloc,1) = derv2(indi,indiloc)
        endif
      endif
    enddo
  endif
!
!  Sum up k point weights and check for gamma point
!
  lgamma = .false.
  lintegerKpoint_present = .false.
  if (lcluster) then
    lgamma = .true.
    lgammaonly = .true.
    sumwkpt = 1.0_dp
  else
    sumwkpt = 0.0_dp
    do i = nlkpt,nukpt
      sumwkpt = sumwkpt + wkpt(i)
      if (abs(xkpt(i)).lt.too_small.and.abs(ykpt(i)).lt.too_small.and.abs(zkpt(i)).lt.too_small) then
        lgamma = .true.
      else
        xmod = mod(xkpt(i)+100.0_dp,1.0_dp)
        ymod = mod(ykpt(i)+100.0_dp,1.0_dp)
        zmod = mod(zkpt(i)+100.0_dp,1.0_dp)
        if (xmod+ymod+zmod.lt.1.0d-10) then
          lintegerKpoint_present = .true.
        endif
      endif
    enddo
!
!  Is the gamma point the only K point?
!
    lgammaonly = (lgamma.and.nllkpt.eq.1)
  endif
!
!  If this is a finite difference second derivative run and there is a non-gamma k-point then stop.
!
  if ((lfinitediff2.or.lnoanald2loc).and..not.lgammaonly) then
    call outerror('Finite difference phonons only allowed for gamma point',0_i4)
    call stopnow('phonond')
  endif
!
  rnokpt = 1.0_dp/sumwkpt
!
!  Ensure we have k vectors
!
  if (ndim.eq.3) then
    call kvector3D
  elseif (ndim.eq.2) then
    call kvector2D
  elseif (ndim.eq.1) then
    call kvector1D
  endif
!
!  Output phonon header
!
  if (lprinloc) then
    if (lcluster) then
      write(ioout,'(/,''  Vibrational Frequency Calculation : '',/)')
    else
      write(ioout,'(/,''  Phonon Calculation : '',/)')
      if ((lgamma.or.lintegerKpoint_present).and.lnonanal) then
        if (nbornstep(ncf).gt.0) then
          write(ioout,'(''  Number of angular steps for n-a correction at gamma = '',i6,/)') nbornstep(ncf)
        else
          lpartofdisp = (nudpt.ge.nldpt)
          if (lpartofdisp) then
            write(ioout,'(''  K direction for n-a correction at gamma ='',3(1x,f8.5))') bornk(1,ncf),bornk(2,ncf),bornk(3,ncf)
            write(ioout,'(''              or the direction of phonon dispersion'',/)')
          else
            write(ioout,'(''  K direction for n-a correction at gamma ='',3(1x,f8.5),/)') bornk(1,ncf),bornk(2,ncf),bornk(3,ncf)
          endif
        endif
      endif
    endif
    if ((lgamma.or.lintegerKpoint_present).and.lraman) then
      if (ramandirtype(ncf).eq.1) then
        write(ioout,'(''  Directions for Raman intensities: Type  =  Cartesian'')')
      else
        write(ioout,'(''  Directions for Raman intensities: Type  =  Fractional '')')
      endif
      write(ioout,'(''                                    In    ='',3(1x,f8.5))') ramandir(1,ncf),ramandir(2,ncf),ramandir(3,ncf)
      write(ioout,'(''                                    Out   ='',3(1x,f8.5),/)') ramandir(4,ncf),ramandir(5,ncf),ramandir(6,ncf)
    endif
    if (lcluster) then
      write(ioout,'(''--------------------------------------------------------------------------------'')')
    else
      write(ioout,'(''  Number of k points for this configuration = '',i8,/)')(nukpt-nlkpt+1)
      write(ioout,'(''--------------------------------------------------------------------------------'')')
    endif
  endif
!**********************************************************************************************
!  If frequencies are to be written to a permanent file then open with the correct type here  *
!**********************************************************************************************
  if (lfrq.and.ioproc) then
    if (index(freqfile,' ').ne.1) then
      if (lfrqbin) then
        open(51,file=freqfile,form='unformatted',status='unknown')
      else
        open(52,file=freqfile,status='unknown')
!
!  Set format statement
!
        fstring = ' '
        call itow(fstring1,nfreqdecimals+5_i4,2_i4)
        call itow(fstring2,nfreqdecimals,2_i4)
        fstring = adjustl(fstring1)//"."//adjustl(fstring2)
      endif
    else
      if (lfrqbin) then
        open(51,form='unformatted',status='unknown')
      else
        open(52,status='unknown')
!
!  Set format statement
!
        fstring = ' '
        call itow(fstring1,nfreqdecimals+5_i4,2_i4)
        call itow(fstring2,nfreqdecimals,2_i4)
        fstring = adjustl(fstring1)//"."//adjustl(fstring2)
      endif
    endif
  endif
!*******************************************************
!  Open eigenvector file if required and write header  *
!*******************************************************
  if (leig.and.ioproc) then
    open(53,file=eigfile,status='unknown',form='formatted')
    write(53,'(i6)') nphonatc
    do i = 1,nphonatc
      ii = nphonatptr(i)
      write(53,'(i3,1x,3(f15.6,1x))') nat(ii),xclat(ii),yclat(ii),zclat(ii)
    enddo
    write(53,'(i6)') nllkpt
    write(53,'(i6)') 3*nphonatc
  endif
!*********************************************************
!  Open CASTEP phonon file if required and write header  *
!*********************************************************
  if (lcas.and.ioproc) then
!
!  Open file
!
    open(54,file=casfile,status='unknown',form='formatted')
!
!  Open file header
!
    write(54,'(a13)') " BEGIN header"
    write(54,'(a22,i6)') " Number of ions       ", nphonatc
    write(54,'(a22,i6)') " Number of branches   ", 3*nphonatc
    write(54,'(a22,i6)') " Number of wavevectors", nllkpt
    write(54,'(a28)') " Frequencies in         cm-1"
!
!  Are intensity units correct?
!
    write(54,'(a36)') " IR intensities in      (D/A)**2/amu" ! ??
    write(54,'(a28)') " Raman intensities in   A**4"         ! ??
!
!  Write out cell vectors
!
    write(54,'(a)')  " Unit cell vectors (A)"
    do i = 1,3
      write(54,'(3f12.6)') (rv(j,i),j=1,3)
    enddo
!
!  Write out a list of fully occupied atoms in the FULL unit cell
!
    write(54,'(a)') " Fractional Co-ordinates"
    do i = 1,nphonatc
      ii = nphonatptr(i)
      call label(nat(ii),nftype(ii),lab1)
      write(54,'(i6,1x,3f12.6,3x,a5,1x,f12.6)') i,xfrac(ii),yfrac(ii),zfrac(ii),lab1,fmass(3*(ii-1)+1)
    enddo
!
!  Close file header
!
    write(54,'(a)') " END header"
  endif
!************************************************
!  Allocate second derivative workspace memory  *
!************************************************
  if (leigloc) then
    maxeigc = maxd2
    maxeigr = maxd2
    if (lgammaonly) then
      allocate(eigr(maxeigr,mtsvloc),stat=status)
      if (status/=0) call outofmemory('phonond','eigr')
      if (nbornstep(ncf).gt.0.or..not.lnonanal) then
        allocate(savederv2(maxeigr,mtvloc),stat=status)
        if (status/=0) call outofmemory('phonond','savederv2')
      endif
    else
      if (lgamma) then
        allocate(eigr(maxeigr,mtsvloc),stat=status)
        if (status/=0) call outofmemory('phonond','eigr')
      endif
      allocate(eigc(maxeigc,mtsvloc),stat=status)
      if (status/=0) call outofmemory('phonond','eigc')
      if (nbornstep(ncf).gt.0.or.lnonanal) then
        allocate(savederv2(maxeigr,mtvloc),stat=status)
        if (status/=0) call outofmemory('phonond','savederv2')
        allocate(savedervi(maxeigr,mtvloc),stat=status)
        if (status/=0) call outofmemory('phonond','savedervi')
      endif
    endif
  else
    if (msv.gt.0) then
      maxeigc = mtsv
      maxeigr = msv 
      if (lgammaonly) then
        allocate(eigr(maxeigr,msv),stat=status)
        if (status/=0) call outofmemory('phonond','eigr')
      else
        if (lgamma) then
          allocate(eigr(maxeigr,msv),stat=status)
          if (status/=0) call outofmemory('phonond','eigr')
        endif
        allocate(eigc(maxeigc,mtsvloc),stat=status)
        if (status/=0) call outofmemory('phonond','eigc')
      endif
    else
      maxeigc = 1
      maxeigr = 1
      if (lgammaonly) then
        maxeigr = mtv
        allocate(eigr(mtv,mtvloc),stat=status)
        if (status/=0) call outofmemory('phonond','eigr')
      else
        maxeigc = mtv
        if (lgamma) then
          maxeigr = mtv
          allocate(eigr(mtv,mtvloc),stat=status)
          if (status/=0) call outofmemory('phonond','eigr')
        endif
        allocate(eigc(mtv,mtvloc),stat=status)
        if (status/=0) call outofmemory('phonond','eigc')
      endif
    endif
  endif
!***********************
!  Setup PDF Variables *
!***********************
  if (lpdfloc) then
    call setuppdfphonon(nllkpt,mtv,nphonatc,nphonatptr,ncf)
  endif
!******************************
!  Initialise MSDs if needed  *
!******************************
  if (lmsd) then
    if (lrigid) then
      msdx(1:numatnomol+nmol) = 0.0_dp
      msdy(1:numatnomol+nmol) = 0.0_dp
      msdz(1:numatnomol+nmol) = 0.0_dp
    else
      msdx(1:ncfoc) = 0.0_dp
      msdy(1:ncfoc) = 0.0_dp
      msdz(1:ncfoc) = 0.0_dp
    endif
  endif
!***********************
!  Loop over k points  *
!***********************
  np_nlkpt = nlkpt
  np_procs = 1_i4
  do k = np_nlkpt,nukpt,np_procs
    if (lpdfloc) then
!
!  Store k point information in PDF module
!
      call setcurk(k,lallk=.true.)
    endif
    if (lscattercall) then
      xkt = xskpt(k)
      ykt = yskpt(k)
      zkt = zskpt(k)
      wk  = wskpt(k)*rnokpt
    else
      xkt = xkpt(k)
      ykt = ykpt(k)
      zkt = zkpt(k)
      wk  = wkpt(k)*rnokpt
    endif
    lgamma = (abs(xkt).lt.too_small.and.abs(ykt).lt.too_small.and.abs(zkt).lt.too_small)
    if (.not.lgamma) then
      xmod = mod(xkt+100.0_dp,1.0_dp)
      ymod = mod(ykt+100.0_dp,1.0_dp)
      zmod = mod(zkt+100.0_dp,1.0_dp)
      lintegerKpoint = (xmod+ymod+zmod.lt.1.0d-10)
    else
      lintegerKpoint = .false.
    endif
!
    if (lfreqout.and.lprinloc.and..not.lcluster) then
      write(ioout,'(''  K point '',i6,'' = '',3f10.6,''  Weight = '',f8.3)') k,xkt,ykt,zkt,wk
      write(ioout,'(''--------------------------------------------------------------------------------'')')
    endif
!
!  Write K point headers to files if required
!
    if (leig.and.ioproc.and..not.lcluster) then
      write(53,'(''K point at '',3f10.6,'' in BZ '')') xkt,ykt,zkt
    endif
    if (lcas.and.ioproc.and..not.lcluster) then
      write(54,'(a10,i5,3f12.6,f14.6)') "     q-pt=",k,xkt,ykt,zkt,wk
    endif
    if (lcluster) then
!
!  Generate derivatives for cluster case by finite differences if not available analytically
!
      if (lnoanald2loc.or.lfinitediff2) then
        call dynamicn
      endif
    else
!
!  Generate phased second derivatives for periodic case
!
      if (lnoanald2loc.or.lfinitediff2) then
        call dynamicn
      else
        call dynamic(xkt,ykt,zkt)
      endif
    endif
!
    if (.not.lnoanald2loc.and..not.lfinitediff2.and..not.lcluster) then
!
!  Include diagonal blocks, stored in derv3
!
      do i = 1,nphonat
        iloc = atom2local(nphonatptr(i))
        if (iloc.gt.0) then
!
!  Atom is local to this node
!
          indi = 3*(nphonatptr(i)-1)
          indiloc = 3*(iloc-1)
          derv2(indi+1,indiloc+1) = derv2(indi+1,indiloc+1) + derv3(indiloc+1,1)
          derv2(indi+2,indiloc+1) = derv2(indi+2,indiloc+1) + derv3(indiloc+2,1)
          derv2(indi+3,indiloc+1) = derv2(indi+3,indiloc+1) + derv3(indiloc+3,1)
          derv2(indi+1,indiloc+2) = derv2(indi+1,indiloc+2) + derv3(indiloc+1,2)
          derv2(indi+2,indiloc+2) = derv2(indi+2,indiloc+2) + derv3(indiloc+2,2)
          derv2(indi+3,indiloc+2) = derv2(indi+3,indiloc+2) + derv3(indiloc+3,2)
          derv2(indi+1,indiloc+3) = derv2(indi+1,indiloc+3) + derv3(indiloc+1,3)
          derv2(indi+2,indiloc+3) = derv2(indi+2,indiloc+3) + derv3(indiloc+2,3)
          derv2(indi+3,indiloc+3) = derv2(indi+3,indiloc+3) + derv3(indiloc+3,3)
          if (nbsmat.gt.0) then
            indi = mint + i
            indiloc = mintloc + iloc
            derv2(indi,indiloc) = derv2(indi,indiloc) + derv3(indiloc,1)
          endif
        endif
      enddo
    endif
!----------------------------------------------------------------------------------------------------
!  Construct second derivative matrix : Atom only case                                              |
!----------------------------------------------------------------------------------------------------
    if (nphonatonnode.gt.0) then
      if (numat.ne.nphonat.and.atom2local(nphonatonnodeptr(nphonatonnode)).ne.nphonatonnode) then
!*********************************************************************
!  Compress full second derivatives down to region 1 only if needed  *
!*********************************************************************
        do i = 1,nphonatonnode
          if (atom2local(nphonatonnodeptr(i)).ne.i) then
            indi  = 3*(i-1)
            indii = 3*(atom2local(nphonatonnodeptr(i))-1)
            do j = 1,nphonat
              indj  = 3*(j-1)
              indjj = 3*(nphonatptr(j)-1)
              derv2(indj+1,indi+1) = derv2(indjj+1,indii+1)
              derv2(indj+2,indi+1) = derv2(indjj+2,indii+1)
              derv2(indj+3,indi+1) = derv2(indjj+3,indii+1)
              derv2(indj+1,indi+2) = derv2(indjj+1,indii+2)
              derv2(indj+2,indi+2) = derv2(indjj+2,indii+2)
              derv2(indj+3,indi+2) = derv2(indjj+3,indii+2)
              derv2(indj+1,indi+3) = derv2(indjj+1,indii+3)
              derv2(indj+2,indi+3) = derv2(indjj+2,indii+3)
              derv2(indj+3,indi+3) = derv2(indjj+3,indii+3)
              if (.not.lgamma) then
                dervi(indj+1,indi+1) = dervi(indjj+1,indii+1)
                dervi(indj+2,indi+1) = dervi(indjj+2,indii+1)
                dervi(indj+3,indi+1) = dervi(indjj+3,indii+1)
                dervi(indj+1,indi+2) = dervi(indjj+1,indii+2)
                dervi(indj+2,indi+2) = dervi(indjj+2,indii+2)
                dervi(indj+3,indi+2) = dervi(indjj+3,indii+2)
                dervi(indj+1,indi+3) = dervi(indjj+1,indii+3)
                dervi(indj+2,indi+3) = dervi(indjj+2,indii+3)
                dervi(indj+3,indi+3) = dervi(indjj+3,indii+3)
              endif
            enddo
          endif
        enddo
      endif
    endif
!**************************************************************
!  Compress second derivative matrix w.r.t. breathing shells  *
!**************************************************************
    if (nbsmat.gt.0) then
!
!  Full occupancy
!
      mis = mint
      mjs = mint
      do i = 1,nbs
        nri = nbsptr(i)
        do j = 1,nphonat
          indj = 3*(j-1)
          indjj = 3*(nphonatptr(j) - 1)
          derv2(indj+1,mis+i) = derv2(indjj+1,mis+nri)
          derv2(indj+2,mis+i) = derv2(indjj+2,mis+nri)
          derv2(indj+3,mis+i) = derv2(indjj+3,mis+nri)
          derv2(mis+i,indj+1) = derv2(mis+nri,indjj+1)
          derv2(mis+i,indj+2) = derv2(mis+nri,indjj+2)
          derv2(mis+i,indj+3) = derv2(mis+nri,indjj+3)
        enddo
        do j = 1,nbs
          nrj = nbsptr(j)
          derv2(mjs+j,mis+i) = derv2(mjs+nrj,mis+nri)
        enddo
        if (.not.lgamma) then
          do j = 1,nphonat
            indj = 3*(j-1)
            indjj = 3*(nphonatptr(j) - 1)
            dervi(indj+1,mis+i) = dervi(indjj+1,mis+nri)
            dervi(indj+2,mis+i) = dervi(indjj+2,mis+nri)
            dervi(indj+3,mis+i) = dervi(indjj+3,mis+nri)
            dervi(mis+i,indj+1) = dervi(mis+nri,indjj+1)
            dervi(mis+i,indj+2) = dervi(mis+nri,indjj+2)
            dervi(mis+i,indj+3) = dervi(mis+nri,indjj+3)
          enddo
          do j = 1,nbs
            nrj = nbsptr(j)
            dervi(mjs+j,mis+i) = dervi(mjs+nrj,mis+nri)
          enddo
        endif
      enddo
    endif
!
!  Output the uncompressed second derivatives for debugging
!
    if (index(keyword,'dyna').ne.0.and.index(keyword,'debu').ne.0) then
      if (ioproc) then
        write(ioout,'(/,''  Uncompressed Real Dynamical matrix :'',/)')
      endif
      if (lioproconly) then
!
!  Allocate temporary workspace for communication
!
        ntmp = mtsv
        ntag = 1
        allocate(dtmp(mtsv,3_i4),stat=status)
        if (status/=0) call outofmemory('phonond','dtmp')
        allocate(StatMPI(MPI_Status_Size),stat=status)
        if (status/=0) call outofmemory('phonond','StatMPI')
      endif
      call mpbarrier
      do i = 1,mtsv,3
        iloc = mtvrptr(i)
        if (lioproconly.and.mtvnptr(i).ne.0_i4) then
!
!  Post receive
!
          if (ioproc) then
            nnode = mtvnptr(i)
            call MPI_IRecv(dtmp,ntmp,MPI_double_precision,nnode, &
                           ntag,MPI_Comm_World,Request,MPIerror)
          endif
!
!  Pass data to ioproc for writing
!
          if (iloc.gt.0) then
            ix = iloc
            iy = ix + 1
            iz = ix + 2
            dtmp(1:mtsv,1) = derv2(1:mtsv,ix)
            dtmp(1:mtsv,2) = derv2(1:mtsv,iy)
            dtmp(1:mtsv,3) = derv2(1:mtsv,iz)
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
            write(ioout,'(12f11.6)')(dtmp(j,1),j=1,mtsv)
            write(ioout,'(12f11.6)')(dtmp(j,2),j=1,mtsv)
            write(ioout,'(12f11.6)')(dtmp(j,3),j=1,mtsv)
          endif
        else
          if (iloc.gt.0) then
            ix = iloc
            iy = ix + 1
            iz = ix + 2
            write(ioout,'(12f11.6)')(derv2(j,ix),j=1,mtsv)
            write(ioout,'(12f11.6)')(derv2(j,iy),j=1,mtsv)
            write(ioout,'(12f11.6)')(derv2(j,iz),j=1,mtsv)
          endif
        endif
        call mpbarrier
      enddo
      if (.not.lgamma) then
        if (ioproc) then
          write(ioout,'(/,''  Uncompressed Imag Dynamical matrix :'',/)')
        endif
        call mpbarrier
        do i = 1,mtsv,3
          iloc = mtvrptr(i)
          if (lioproconly.and.mtvnptr(i).ne.0_i4) then
!
!  Post receive
!
            if (ioproc) then
              nnode = mtvnptr(i)
              call MPI_IRecv(dtmp,ntmp,MPI_double_precision,nnode, &
                             ntag,MPI_Comm_World,Request,MPIerror)
            endif
!
!  Pass data to ioproc for writing
!
            if (iloc.gt.0) then
              ix = iloc
              iy = ix + 1
              iz = ix + 2
              dtmp(1:mtsv,1) = dervi(1:mtsv,ix)
              dtmp(1:mtsv,2) = dervi(1:mtsv,iy)
              dtmp(1:mtsv,3) = dervi(1:mtsv,iz)
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
              write(ioout,'(12f11.6)')(dtmp(j,1),j=1,mtsv)
              write(ioout,'(12f11.6)')(dtmp(j,2),j=1,mtsv)
              write(ioout,'(12f11.6)')(dtmp(j,3),j=1,mtsv)
            endif
          else
            if (iloc.gt.0) then
              ix = iloc
              iy = ix + 1
              iz = ix + 2
              write(ioout,'(12f11.6)')(dervi(j,ix),j=1,mtsv)
              write(ioout,'(12f11.6)')(dervi(j,iy),j=1,mtsv)
              write(ioout,'(12f11.6)')(dervi(j,iz),j=1,mtsv)
            endif
          endif
          call mpbarrier
        enddo
      endif
      if (lioproconly) then
        deallocate(StatMPI,stat=status)
        if (status/=0) call deallocate_error('phonond','StatMPI')
        deallocate(dtmp,stat=status)
        if (status/=0) call deallocate_error('phonond','dtmp')
      endif
    endif
!**********************************
!  Eliminate shell contributions  *
!**********************************
    if (msv.gt.0) then
      if (lgamma) then
!**************************
!  Real Matrix Inversion  *
!**************************
        t1i = g_cpu_time()
        ifail = 0                   
#ifdef INTEL
        leigcOK = (allocated(eigc))
        if (.not.leigcOK) then
          allocate(eigc(maxeigc,mtsvloc),stat=status)
          if (status/=0) call outofmemory('phonond','eigc')
        endif
!
!  Copy to eigc as workspace
!
        do i = 1,msvloc
          do j = 1,msv
            eigc(mcv+j,mcvloc+i) = dcmplx(derv2(mcv+j,mcvloc+i),0.0_dp)
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
          call stopnow('phonond')
        endif
!
!  Transfer data back
!
        do i = 1,msvloc
          do j = 1,msv
            derv2(mcv+j,mcvloc+i) = dble(eigc(mcv+j,mcvloc+i))
          enddo
        enddo
!
        if (.not.leigcOK) then
          deallocate(eigc,stat=status)
          if (status/=0) call deallocate_error('phonond','eigc')
        endif
#else
!
!  Call library to invert shell-shell matrix 
!
        call matrix_inversion_shells(msv,mcv+1_i4,maxd2,derv2,nshell,nshellonnode,ifail)
!
!  Check return flag
!
        if (ifail.ne.0) then
          call outerror('inversion of shell 2nd derivatives failed',0_i4)
          call stopnow('phonond')
        endif
#endif
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
        call descinit( idesc, ncs, ncs, 3*nb, 3*nb, 0, 0, iBlacsContext, ld, ifails )
!
        if (ifails.ne.0) then
          call outerror('initialisation in descinit failed',0_i4)
          call stopnow('phonond')
        endif
!
!  First pass : S-1*T stored in diagonally opposite copy of T
!
        call pdgemm('N','T',msv,mcv,msv,1.0d0,derv2,mcv+1,mcv+1,idesc,derv2,1,mcv+1,idesc,0.0d0,derv2,mcv+1,1,idesc)
!
!  Second pass : T*(S-1*T)
!
        call pdgemm('N','N',mcv,mcv,msv,-1.0d0,derv2,1,mcv+1,idesc,derv2,mcv+1,1,idesc,1.0d0,derv2,1,1,idesc)
      else
!*****************************
!  Complex Matrix Inversion  *
!*****************************
        t1i = g_cpu_time()
        ifail = 0                   
!
!  Copy to eigc as workspace
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
          call stopnow('phonond')
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
          call stopnow('phonond')
        endif
!
!  Transfer data to eigc as complex version of derv2/dervi
!
        do i = 1,mcvloc+msvloc
          do j = 1,mcv
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
        if (lgroupvelocity) then
!
!  Correct group velocities for shells
!
          do ii = 1,3
!
!  Copy core-core component of derv2dk to eigc
!
            do i = 1,mcvloc
              do j = 1,mcv
                eigc(j,i) = derv2dk(ii,j,i)
              enddo
            enddo
!
!  Copy shell-shell component of derv2dk to eigc
!
            do i = 1,msvloc
              do j = 1,msv
                eigc(mcv+j,mcvloc+i) = derv2dk(ii,mcv+j,mcvloc+i)
              enddo
            enddo
!
!  First pass 
!
            call pzgemm('C','N',mcv,msv,msv,(1.0d0,0.0d0),eigc,mcv+1,1,idesc,eigc,mcv+1,mcv+1, &
                        idesc,(0.0d0,0.0d0),eigc,1,mcv+1,idesc)
!
!  Second pass 
!
            call pzgemm('N','N',mcv,mcv,msv,(1.0d0,0.0d0),eigc,1,mcv+1,idesc,eigc,mcv+1,1, &
                        idesc,(1.0d0,0.0d0),eigc,1,1,idesc)
!
!  Copy back core-core component of derv2dk from eigc
!
            do i = 1,mcvloc
              do j = 1,mcv
                derv2dk(ii,j,i) = eigc(j,i)
              enddo
            enddo
          enddo
!
!  Second term : core-shell x shell-core
!
!  Correct group velocities for shells
!
          do ii = 1,3
!
!  Copy shell-core component of derv2dk to eigc
!
            do i = 1,msvloc
              do j = 1,mcv
                eigc(j,mcvloc+i) = derv2dk(ii,j,mcvloc+i)
              enddo
            enddo
!
!  Multiply core-shell of eigc with shell-core of derv2dk
!
            call pzgemm('N','N',mcv,mcv,msv,(-1.0d0,0.0d0),eigc,1,mcv+1,idesc,eigc,mcv+1,1, &
                        idesc,(0.0d0,0.0d0),eigc,1,1,idesc)
            call pzgemm('C','C',mcv,mcv,msv,(-1.0d0,0.0d0),eigc,mcv+1,1,idesc,eigc,1,mcv+1, &
                        idesc,(1.0d0,0.0d0),eigc,1,1,idesc)
!
!  Copy back core-core component of derv2dk from eigc
!
            do i = 1,mcvloc
              do j = 1,mcv
                derv2dk(ii,j,i) = derv2dk(ii,j,i) + eigc(j,i)
              enddo
            enddo
          enddo
        endif
      endif
    endif
!****************************
!  End of shell correction  *
!****************************
!
!  Output the second derivatives after shell removal for debugging
!
    if (lrigid.and.index(keyword,'dyna').ne.0.and.index(keyword,'debu').ne.0) then
      if (ioproc) then
        write(ioout,'(/,''  Real Dynamical matrix for the cores :'',/)')
      endif
      if (lioproconly) then
!
!  Allocate temporary workspace for communication
!
        ntmp = mcv
        ntag = 1
        allocate(dtmp(mcv,3_i4),stat=status)
        if (status/=0) call outofmemory('phonond','dtmp')
        allocate(StatMPI(MPI_Status_Size),stat=status)
        if (status/=0) call outofmemory('phonond','StatMPI')
      endif
      call mpbarrier
      do i = 1,mcv,3
        iloc = mtvrptr(i)
        if (lioproconly.and.mtvnptr(i).ne.0_i4) then
!
!  Post receive
!
          if (ioproc) then
            nnode = mtvnptr(i)
            call MPI_IRecv(dtmp,ntmp,MPI_double_precision,nnode, &
                           ntag,MPI_Comm_World,Request,MPIerror)
          endif
!
!  Pass data to ioproc for writing
!
          if (iloc.gt.0) then
            ix = iloc
            iy = ix + 1
            iz = ix + 2
            dtmp(1:mcv,1) = derv2(1:mcv,ix)
            dtmp(1:mcv,2) = derv2(1:mcv,iy)
            dtmp(1:mcv,3) = derv2(1:mcv,iz)
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
            write(ioout,'(12f11.6)')(dtmp(j,1),j=1,mcv)
            write(ioout,'(12f11.6)')(dtmp(j,2),j=1,mcv)
            write(ioout,'(12f11.6)')(dtmp(j,3),j=1,mcv)
          endif
        else
          if (iloc.gt.0) then
            ix = iloc
            iy = ix + 1
            iz = ix + 2
            write(ioout,'(12f11.6)')(derv2(j,ix),j=1,mcv)
            write(ioout,'(12f11.6)')(derv2(j,iy),j=1,mcv)
            write(ioout,'(12f11.6)')(derv2(j,iz),j=1,mcv)
          endif
        endif
        call mpbarrier
      enddo
      if (.not.lgamma) then
        if (ioproc) then
          write(ioout,'(/,''  Uncompressed Imag Dynamical matrix for the cores :'',/)')
        endif
        call mpbarrier
        do i = 1,mcv,3
          iloc = mtvrptr(i)
          if (lioproconly.and.mtvnptr(i).ne.0_i4) then
!
!  Post receive
!
            if (ioproc) then
              nnode = mtvnptr(i)
              call MPI_IRecv(dtmp,ntmp,MPI_double_precision,nnode, &
                             ntag,MPI_Comm_World,Request,MPIerror)
            endif
!
!  Pass data to ioproc for writing
!
            if (iloc.gt.0) then
              ix = iloc
              iy = ix + 1
              iz = ix + 2
              dtmp(1:mcv,1) = dervi(1:mcv,ix)
              dtmp(1:mcv,2) = dervi(1:mcv,iy)
              dtmp(1:mcv,3) = dervi(1:mcv,iz)
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
              write(ioout,'(12f11.6)')(dtmp(j,1),j=1,mcv)
              write(ioout,'(12f11.6)')(dtmp(j,2),j=1,mcv)
              write(ioout,'(12f11.6)')(dtmp(j,3),j=1,mcv)
            endif
          else
            if (iloc.gt.0) then
              ix = iloc
              iy = ix + 1
              iz = ix + 2
              write(ioout,'(12f11.6)')(dervi(j,ix),j=1,mcv)
              write(ioout,'(12f11.6)')(dervi(j,iy),j=1,mcv)
              write(ioout,'(12f11.6)')(dervi(j,iz),j=1,mcv)
            endif
          endif
          call mpbarrier
        enddo
      endif
      if (lioproconly) then
        deallocate(StatMPI,stat=status)
        if (status/=0) call deallocate_error('phonond','StatMPI')
        deallocate(dtmp,stat=status)
        if (status/=0) call deallocate_error('phonond','dtmp')
      endif
    endif
!**********************************************
!  Option to write out a .frc file for QMPOT  *
!**********************************************
    if (lfrc.and.lgamma) call outfrc(fc,.true.,.true.)
!************************
!  Eckart purification  *
!************************
    if (leckart.and.lcluster) then
      call eckart(mcv,mcvloc,mcvptr,maxd2,derv2,maxd2,dervi)
    endif
!----------------------------------------------------------------------------------------------------
!  Construct second derivative matrix : Rigid molecule case                                         |
!----------------------------------------------------------------------------------------------------
    if (lrigid) then
!****************************************
!  Build rotational second derivatives  *
!****************************************
      call rigidmoleculephon(.not.lgamma)
      if (lgroupvelocity) call rigidmoleculegv
!***************************************************
!  Compress second derivative matrix w.r.t. cores  *
!***************************************************
      if (.not.lgamma) then
        allocate(dtmp(3*numatnomol,6),stat=status)
        if (status/=0) call outofmemory('phonond','dtmp')
        if (lgroupvelocity) then
          allocate(ctmp(3,3*numatnomol,3),stat=status)
          if (status/=0) call outofmemory('phonond','ctmp')
        endif
      else
        allocate(dtmp(3*numatnomol,3),stat=status)
        if (status/=0) call outofmemory('phonond','dtmp')
      endif
      allocate(StatMPI(MPI_Status_Size),stat=status)
      if (status/=0) call outofmemory('property3','StatMPI')
!
      do i = 1,ncorenomol
!
!  Is derv2 local to this node?
!
        snode = (atom2node(ncorenomolptr(i)))
        if (snode.eq.procid) then
          ld2loc = .true.
          indif = 3*(atom2local(ncorenomolptr(i)) - 1)
        else
          ld2loc = .false.
        endif
!
!  Is dervi local to this node?
!
        rnode = mtvnptr(3*(i-1)+1)
        if (rnode.eq.procid) then
          ldiloc = .true.
          iloc = mtvrptr(3*(i-1)+1)
          indil = iloc -1
        else
          ldiloc = .false.
        endif
!
!  If derv2 and dervi are both local then copy
!
        if (ld2loc.and.ldiloc) then
          do j = 1,ncorenomol
            indjl = 3*(j-1)
            indjf = 3*(ncorenomolptr(j)-1)
            do ix = 1,3
              do jx = 1,3
                derv2(indjl+jx,indil+ix) = derv2(indjf+jx,indif+ix)
              enddo
            enddo
          enddo
          if (.not.lgamma) then
            do j = 1,ncorenomol
              indjl = 3*(j-1)
              indjf = 3*(ncorenomolptr(j)-1)
              do ix = 1,3
                do jx = 1,3
                  dervi(indjl+jx,indil+ix) = dervi(indjf+jx,indif+ix)
                enddo
              enddo
            enddo
            if (lgroupvelocity) then
              do j = 1,ncorenomol
                indjl = 3*(j-1)
                indjf = 3*(ncorenomolptr(j)-1)
                do ix = 1,3
                  do jx = 1,3
                    derv2dk(1:3,indjl+jx,indil+ix) = derv2dk(1:3,indjf+jx,indif+ix)
                  enddo
                enddo
              enddo
            endif
          endif
        elseif (ldiloc) then
!
!  If dervi is local then receive
!
          if (.not.lgamma) then
            nsize = 18*ncorenomol
          else
            nsize = 9*ncorenomol
          endif
          ntag = i
          call MPI_IRecv(dtmp,nsize,MPI_double_precision,snode,ntag,MPI_Comm_World,Request,MPIerror)
!
          if (.not.lgamma) then
            if (lgroupvelocity) then
              ntag = ncorenomol + i
              nsizec = 27*ncorenomol
              call MPI_IRecv(ctmp,nsizec,MPI_double_complex,snode,ntag,MPI_Comm_World,Requestc,MPIerror)
            endif
          endif
        elseif (ld2loc) then
!
!  If derv2 is local then send
!
          do j = 1,ncorenomol
            indjl = 3*(j-1)
            indjf = 3*(ncorenomolptr(j)-1)
            do ix = 1,3
              do jx = 1,3
                dtmp(indjl+jx,ix) = derv2(indjf+jx,indif+ix)
              enddo
            enddo
          enddo
          if (.not.lgamma) then
            do j = 1,ncorenomol
              indjl = 3*(j-1)
              indjf = 3*(ncorenomolptr(j)-1)
              do ix = 1,3
                do jx = 1,3
                  dtmp(indjl+jx,3+ix) = dervi(indjf+jx,indif+ix)
                enddo
              enddo
            enddo
          endif
          if (.not.lgamma) then
            nsize = 18*ncorenomol
          else
            nsize = 9*ncorenomol
          endif
          ntag = i
          call MPI_ISend(dtmp,nsize,MPI_double_precision,snode,ntag,MPI_Comm_World,Request,MPIerror)
!
          if (.not.lgamma) then
            if (lgroupvelocity) then
              do j = 1,ncorenomol
                indjl = 3*(j-1)
                indjf = 3*(ncorenomolptr(j)-1)
                do ix = 1,3
                  do jx = 1,3
                    ctmp(1:3,indjl+jx,3+ix) = derv2dk(1:3,indjf+jx,indif+ix)
                  enddo
                enddo
              enddo
              ntag = ncorenomol + i
              nsizec = 27*ncorenomol
              call MPI_ISend(ctmp,nsizec,MPI_double_complex,snode,ntag,MPI_Comm_World,Requestc,MPIerror)
            endif
          endif
        endif
!
!  Wait for communication to finish
!
        if ((ld2loc.and..not.ldiloc).or.(ldiloc.and..not.ld2loc)) then
          call MPI_WaitAll(1,Request,StatMPI,MPIerror)
          if (ldiloc) then
!
!  After communication place data in dervi
!
            do j = 1,ncorenomol
              indjl = 3*(j-1)
              do ix = 1,3
                do jx = 1,3
                  derv2(indjl+jx,indil+ix) = dtmp(indjl+jx,ix)
                enddo
              enddo
            enddo
            if (.not.lgamma) then
              do j = 1,ncorenomol
                indjl = 3*(j-1)
                do ix = 1,3
                  do jx = 1,3
                    dervi(indjl+jx,indil+ix) = dtmp(indjl+jx,3+ix)
                  enddo
                enddo
              enddo
            endif
          endif
          if (.not.lgamma) then
            if (lgroupvelocity) then
              call MPI_WaitAll(1,Requestc,StatMPI,MPIerror)
              if (ldiloc) then
!
!  After communication place data in derv2dk
!
                do j = 1,ncorenomol
                  indjl = 3*(j-1)
                  do ix = 1,3
                    do jx = 1,3
                      derv2dk(1:3,indjl+jx,indil+ix) = ctmp(1:3,indjl+jx,ix)
                    enddo
                  enddo
                enddo
              endif
            endif
          endif
        endif
      enddo
!
      deallocate(StatMPI,stat=status)
      if (status/=0) call deallocate_error('phonond','StatMPI')
      if (.not.lgamma.and.lgroupvelocity) then
        deallocate(ctmp,stat=status)
        if (status/=0) call deallocate_error('phonond','ctmp')
      endif
      deallocate(dtmp,stat=status)
      if (status/=0) call deallocate_error('phonond','dtmp')
!****************************************
!  Add in rigid molecule contributions  *
!****************************************
!
!  Copy molecule terms
!
      mind = 3*ncorenomol
      do i = 1,ncorenomol
        if (mtvnptr(3*(i-1)+1).eq.procid) then
          indil = mtvrptr(3*(i-1)+1) - 1
          indif = 3*(ncorenomolptr(i)-1)
          do j = 1,nmol
            mindj = mind + 6*(j-1)
!
!  Copy rigid molecule - atom terms
!
            do ix = 1,3
              do jx = 1,3
                derv2(mindj+jx,indil+ix) = molTCdrv(indif+ix,jx,j)
              enddo
            enddo
            do ix = 1,3
              do jx = 1,3
                derv2(mindj+3+jx,indil+ix) = molQCdrv(indif+ix,jx,j)
              enddo
            enddo
            if (.not.lgamma) then
              do ix = 1,3
                do jx = 1,3
                  dervi(mindj+jx,indil+ix) = - molTCdri(indif+ix,jx,j)
                enddo
              enddo
              do ix = 1,3
                do jx = 1,3
                  dervi(mindj+3+jx,indil+ix) = - molQCdri(indif+ix,jx,j)
                enddo
              enddo
              if (lgroupvelocity) then
                do ix = 1,3
                  do jx = 1,3
                    derv2dk(1:3,mindj+jx,indil+ix) = conjg(molTCdk(1:3,indif+ix,jx,j))
                  enddo
                enddo
                do ix = 1,3
                  do jx = 1,3
                    derv2dk(1:3,mindj+3+jx,indil+ix) = conjg(molQCdk(1:3,indif+ix,jx,j))
                  enddo
                enddo
              endif
            endif
          enddo
        endif
      enddo
      do i = 1,nmol
        iloc = ncorenomol + 2*(i-1) + 1
        if (mtvnptr(3*(iloc-1)+1).eq.procid) then
          mindi = mtvrptr(3*(iloc-1)+1) - 1
!
!  Copy rigid molecule - atom terms for translation
!
          do j = 1,ncorenomol
            indjl = 3*(j-1)
            indjf = 3*(ncorenomolptr(j)-1)
            do ix = 1,3
              do jx = 1,3
                derv2(indjl+jx,mindi+ix) = molTCdrv(indjf+jx,ix,i)
              enddo
            enddo
            if (.not.lgamma) then
              do ix = 1,3
                do jx = 1,3
                  dervi(indjl+jx,mindi+ix) = molTCdri(indjf+jx,ix,i)
                enddo
              enddo
              if (lgroupvelocity) then
                do ix = 1,3
                  do jx = 1,3
                    derv2dk(1:3,indjl+jx,mindi+ix) = molTCdk(1:3,indjf+jx,ix,i)
                  enddo
                enddo
              endif
            endif
          enddo
!
!  Copy rigid molecule - rigid molecule terms
!
          do j = 1,nmol
            mindj = mind + 6*(j-1)
            do ix = 1,3
              do jx = 1,3
                derv2(mindj+jx,mindi+ix) = molTTdrv(jx,ix,j,i)
              enddo
            enddo
            do ix = 1,3
              do jx = 1,3
                derv2(mindj+3+jx,mindi+ix) = molQTdrv(jx,ix,j,i)
              enddo
            enddo
            if (.not.lgamma) then
              do ix = 1,3
                do jx = 1,3
                  dervi(mindj+jx,mindi+ix) = molTTdri(jx,ix,j,i)
                enddo
              enddo
              do ix = 1,3
                do jx = 1,3
                  dervi(mindj+3+jx,mindi+ix) = - molQTdri(jx,ix,j,i)
                enddo
              enddo
              if (lgroupvelocity) then
                do ix = 1,3
                  do jx = 1,3
                    derv2dk(1:3,mindj+jx,mindi+ix) = molTTdk(1:3,jx,ix,j,i)
                  enddo
                enddo
                do ix = 1,3
                  do jx = 1,3
                    derv2dk(1:3,mindj+3+jx,mindi+ix) = conjg(molQTdk(1:3,jx,ix,j,i))
                  enddo
                enddo
              endif
            endif
          enddo
        endif
        if (mtvnptr(3*iloc+1).eq.procid) then
          mindi = mtvrptr(3*iloc+1) - 1
!
!  Copy rigid molecule - atom terms for rotation
!
          do j = 1,ncorenomol
            indjl = 3*(j-1)
            indjf = 3*(ncorenomolptr(j)-1)
            do ix = 1,3
              do jx = 1,3
                derv2(indjl+jx,mindi+ix) = molQCdrv(indjf+jx,ix,i)
              enddo
            enddo
            if (.not.lgamma) then
              do ix = 1,3
                do jx = 1,3
                  dervi(indjl+jx,mindi+ix) = molQCdri(indjf+jx,ix,i)
                enddo
              enddo
              if (lgroupvelocity) then
                do ix = 1,3
                  do jx = 1,3
                    derv2dk(1:3,indjl+jx,mindi+ix) = molQCdk(1:3,indjf+jx,ix,i)
                  enddo
                enddo
              endif
            endif
          enddo
!
!  Copy rigid molecule - rigid molecule terms
!
          do j = 1,nmol
            mindj = mind + 6*(j-1)
            do ix = 1,3
              do jx = 1,3
                derv2(mindj+jx,mindi+ix) = molQTdrv(ix,jx,i,j)
                derv2(mindj+3+jx,mindi+ix) = molQQdrv(jx,ix,j,i)
              enddo
            enddo
            if (.not.lgamma) then
              do ix = 1,3
                do jx = 1,3
                  dervi(mindj+jx,mindi+ix) = molQTdri(ix,jx,i,j)
                  dervi(mindj+3+jx,mindi+ix) = molQQdri(jx,ix,j,i)
                enddo
              enddo
              if (lgroupvelocity) then
                do ix = 1,3
                  do jx = 1,3
                    derv2dk(1:3,mindj+jx,mindi+ix) = molQTdk(1:3,ix,jx,i,j)
                    derv2dk(1:3,mindj+3+jx,mindi+ix) = molQQdk(1:3,jx,ix,j,i)
                  enddo
                enddo
              endif
            endif
          enddo
        endif
      enddo
    endif
!*****************************
!  Multiply by mass-factors  *
!*****************************
    if (lgamma) then
      do i = 1,mtvloc
        rmassi = rfmass(mtvptr(i))
        do j = 1,mtv
          derv2(j,i) = rmassi*rfmass(j)*derv2(j,i)
        enddo
      enddo
    else
      do i = 1,mtvloc
        rmassi = rfmass(mtvptr(i))
        do j = 1,mtv
          derv2(j,i) = rmassi*rfmass(j)*derv2(j,i)
          dervi(j,i) = rmassi*rfmass(j)*dervi(j,i)
        enddo
      enddo
    endif
!
!  If debugging print out dynamical matrix
!
    if (index(keyword,'dyna').ne.0) then
      if (ioproc) then
        write(ioout,'(/,''  Real Dynamical matrix :'',/)')
      endif
      if (lioproconly) then
!
!  Allocate temporary workspace for communication
!
        ntmp = mtv
        ntag = 1
        allocate(dtmp(mtv,1_i4),stat=status)
        if (status/=0) call outofmemory('phonond','dtmp')
        allocate(StatMPI(MPI_Status_Size),stat=status)
        if (status/=0) call outofmemory('phonond','StatMPI')
      endif
      call mpbarrier
      do i = 1,mtv
        iloc  = mtvrptr(i)
        if (lioproconly.and.mtvnptr(i).ne.0_i4) then
!
!  Post receive
!
          if (ioproc) then
            nnode = mtvnptr(i)
            call MPI_IRecv(dtmp,ntmp,MPI_double_precision,nnode, &
                           ntag,MPI_Comm_World,Request,MPIerror)
          endif
!
!  Pass data to ioproc for writing
!
          if (iloc.gt.0) then
            dtmp(1:mtv,1) = derv2(1:mtv,iloc)
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
            write(ioout,'(12f11.6)')(dtmp(j,1),j=1,mtv)
          endif
        else
          if (iloc.gt.0) then
            write(ioout,'(12f11.6)')(derv2(j,iloc),j=1,mtv)
          endif
        endif
        call mpbarrier
      enddo
      if (.not.lgamma) then
        if (ioproc) then
          write(ioout,'(/,''  Imaginary Dynamical matrix :'',/)')
        endif
        call mpbarrier
        do i = 1,mtv
          iloc = mtvrptr(i)
          if (lioproconly.and.mtvnptr(i).ne.0_i4) then
!
!  Post receive
!
            if (ioproc) then
              nnode = mtvnptr(i)
              call MPI_IRecv(dtmp,ntmp,MPI_double_precision,nnode, &
                             ntag,MPI_Comm_World,Request,MPIerror)
            endif
!
!  Pass data to ioproc for writing
!
            if (iloc.gt.0) then
              dtmp(1:mtv,1) = dervi(1:mtv,iloc)
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
              write(ioout,'(12f11.6)')(dtmp(j,1),j=1,mtv)
            endif
          else
            if (iloc.gt.0) then
              write(ioout,'(12f11.6)')(dervi(j,iloc),j=1,mtv)
            endif
          endif
          call mpbarrier
        enddo
      endif
      if (lioproconly) then
        deallocate(StatMPI,stat=status)
        if (status/=0) call deallocate_error('phonond','StatMPI')
        deallocate(dtmp,stat=status)
        if (status/=0) call deallocate_error('phonond','dtmp')
      endif
    endif
!*********************************
!  Diagonalise dynamical matrix  *
!*********************************
    ifail = 0
    if (leigloc) then
!*********************************
!  Eigenvectors and eigenvalues  *
!*********************************
      if (lgamma) then
!
!  Calculate eigenvalues and eigenvectors of uncorrected dynamical matrix
!
        call pdiaggd(mtv,mtvloc,maxd2,derv2,maxeigr,eigr,freq(1,k-nlkpt+1),fscale,.true.,lprint,ifail)
!
!  Calculate the oscillator strengths
!
        call oscillatorstrengthg(mtv,mtvloc,mtvptr,nphonatc,nphonatptr,nphonatc,nphonatptr, &
                                 eigr,maxeigr,oscstrength)
!
!  Optionally call Raman strengths
!
        if (lraman) then
          call ramanstrength(mtv,mtvloc,mtvptr,nphonatc,nphonatptr,nphonatc,nphonatptr, &
                             eigr,maxeigr,ramstrength)
        endif
!
!  Compute Grueneisen parameters
!
        if (lgrueneisen) then
          call gruneisengd(k-nlkpt+1_i4,mtv,mtvloc,mtvptr,msv,nphonatc,eigr,maxeigr,derv2,maxd2,fscale,lprint)
        endif
!
!  Store uncorrected frequencies
!
        savefreq(1:mtv) = freq(1:mtv,k-nlkpt+1)
        if (nbornstep(ncf).eq.0.and.lnonanal) then
!
!  Non-analytic correction for gamma point
!
          if (lnonanal) then
!
!  Find if this point is part of a dispersion curve and if it is substitute direction of approach
!
            lpartofdisp = (nudpt.ge.nldpt)
            if (lpartofdisp) then
              if (k.ge.ndds(nldpt).and.k.le.ndde(nudpt)) then
                lfound = .false.
                nd = nldpt - 1
                do while (.not.lfound.and.nd.lt.nudpt)
                  nd = nd + 1
                  lfound = (k.ge.ndds(nd).and.k.le.ndde(nd))
                enddo
                bornkloc(1) = xkpt(ndde(nd)) - xkpt(ndds(nd))
                bornkloc(2) = ykpt(ndde(nd)) - ykpt(ndds(nd))
                bornkloc(3) = zkpt(ndde(nd)) - zkpt(ndds(nd))
              else
                lpartofdisp = .false.
              endif
            endif
            if (lpartofdisp) then
              if (lrigid) then
                call nagammarigid(mtvrptr,rfmass,bornkloc,maxd2,derv2)
              else
                call nagamma(nphonatc,nphonatc,nphonatptr,rfmass,bornkloc,maxd2,derv2)
              endif
            else
              if (lrigid) then
                call nagammarigid(mtvrptr,rfmass,bornk(1,ncf),maxd2,derv2)
              else
                call nagamma(nphonatc,nphonatc,nphonatptr,rfmass,bornk(1,ncf),maxd2,derv2)
              endif
            endif
          endif
!
!  Calculate eigenvalues and eigenvectors of corrected dynamical matrix
!
          call pdiaggd(mtv,mtvloc,maxd2,derv2,maxeigr,eigr,freq(1,k-nlkpt+1),fscale,.false.,lprint,ifail)
        elseif (nbornstep(ncf).ne.0) then
!******************************************************
!  Polycrystalline average of gamma point correction  *
!******************************************************
!
!  Allocate extra storage
!
          allocate(averfreq(mtv),stat=status)
          if (status/=0) call outofmemory('phonond','averfreq')
!
!  Store derv2
!
          savederv2(1:mtv,1:mtvloc) = derv2(1:mtv,1:mtvloc)
          averfreq(1:mtv) = 0.0_dp
!
!  Set theta and phi step sizes
!
          weightpt = 1.0_dp/dble(nbornstep(ncf))
          phistep = 0.5_dp*pi*weightpt
          thetastep = phistep
!
!  Loop over theta and phi in nbornstep divisions per angle
!
          weightpt = 0.0_dp
          phi = - phistep
          do np = 0,nbornstep(ncf)
            phi = phi + phistep
            theta = - thetastep
            if (np.eq.0) then
              ntmax = 1
            elseif (np.eq.nbornstep(ncf)) then
              ntmax = 2*nbornstep(ncf)
            else
              ntmax = 4*nbornstep(ncf)
            endif
            do nt = 1,ntmax
              theta = theta + thetastep
              weightpt = weightpt + 1.0_dp
!
!  Find current Q direction of approach
!
              qfrac(1) = sin(phi)*cos(theta)
              qfrac(2) = sin(phi)*sin(theta)
              qfrac(3) = cos(phi)
!
!  Restore derv2
!
              derv2(1:mtv,1:mtvloc) = savederv2(1:mtv,1:mtvloc)
!
!  Non-analytic correction for gamma point
!
              if (lrigid) then
                call nagammarigid(mtvrptr,rfmass,qfrac,maxd2,derv2)
              else
                call nagamma(nphonatc,nphonatc,nphonatptr,rfmass,qfrac,maxd2,derv2)
              endif
!
!  Calculate eigenvalues and eigenvectors of corrected dynamical matrix
!
              call pdiaggd(mtv,mtvloc,maxd2,derv2,maxeigr,eigr,freq(1,k-nlkpt+1),fscale,.false.,lprint,ifail)
!
!  Find correspondence of modes by projection of eigenvectors
!

!
!  Add frequencies to average for appropriate mode
!
              do nf = 1,mtv
                averfreq(nf) = averfreq(nf) + freq(nf,k-nlkpt+1)
              enddo
!
!  End loops over angles
!
            enddo
          enddo
!
!  Move average frequencies back to main frequency array
!
          freq(1:mtv,k-nlkpt+1) = averfreq(1:mtv)/weightpt
!
          deallocate(averfreq,stat=status)
          if (status/=0) call deallocate_error('phonond','averfreq')
        endif
!
!  Output frequencies / DOS / intensities
!
        call peigengd(mtv,mtvloc,mtvptr,mtvnptr,mtvrptr,freq(1,k-nlkpt+1),ncore,nphonatc,nphonats, &
                      nphonatptr,maxeigr,eigr,oscstrength,ramstrength,IRintensity(1,k-nlkpt+1)) 
        if (nobsmode.gt.0) then
!
!  Look for vibrational mode frequencies projected on to eigenvectors
!
          do nobm = 1,nobsmode
            call getvibmoded(nobsmodeptr0+nobm,maxeigr,eigr,mtvloc,mtvptr,nfitmode,overlap)
            fobsmodefreq(nobsmodeptr0+nobm) = freq(nfitmode,k-nlkpt+1)
            fobsmodeover(nobsmodeptr0+nobm) = overlap
          enddo
        endif
!
!  Thermal conductivity calculation
!
        if (lthermal) then
          call thermalconductivity_af_d(mtv,mtvloc,mtvptr,mtvrptr,derv2,eigr,Sij,freq(1,k-nlkpt+1), &
                                        nphonatc,nphonatptr,maxd2,fscale)
        endif
!
!  Lower symmetry to remove imaginary modes if selected
!
        call lower(mtv,mtvrptr,freq(1,k-nlkpt+1),nphonatc,nphonatptr,nphonatc,nphonatptr,maxeigr,eigr)
      else
!*****************************************
!  Complex eigenvalues and eigenvectors  *
!*****************************************
        if (lintegerKpoint) then
!
!  Store derv2
!
          if (nbornstep(ncf).gt.0.or.lnonanal) then
            savederv2(1:mtv,1:mtvloc) = derv2(1:mtv,1:mtvloc)
            savedervi(1:mtv,1:mtvloc) = dervi(1:mtv,1:mtvloc)
          endif
!
!  Calculate eigenvalues and eigenvectors of dynamical matrix
!
          call pdiagd(mtv,mtvloc,maxd2,derv2,dervi,maxeigc,eigc,freq(1,k-nlkpt+1),fscale,.true.,lprint,ifail)
!
!  Store uncorrected frequencies
!
          savefreq(1:mtv) = freq(1:mtv,k-nlkpt+1)
!
!  Restore derv2
!
          if (nbornstep(ncf).gt.0.or.lnonanal) then
            derv2(1:mtv,1:mtvloc) = savederv2(1:mtv,1:mtvloc)
            dervi(1:mtv,1:mtvloc) = savedervi(1:mtv,1:mtvloc)
          endif
!
          if (nbornstep(ncf).eq.0.and.lnonanal) then
!
!  Non-analytic correction for gamma point
!
            if (lnonanal) then
!
!  Find if this point is part of a dispersion curve and if it is substitute direction of approach
!
              lpartofdisp = (nudpt.ge.nldpt)
              if (lpartofdisp) then
                if (k.ge.ndds(nldpt).and.k.le.ndde(nudpt)) then
                  lfound = .false.
                  nd = nldpt - 1
                  do while (.not.lfound.and.nd.lt.nudpt)
                    nd = nd + 1
                    lfound = (k.ge.ndds(nd).and.k.le.ndde(nd))
                  enddo
                  bornkloc(1) = xkpt(ndde(nd)) - xkpt(ndds(nd))
                  bornkloc(2) = ykpt(ndde(nd)) - ykpt(ndds(nd))
                  bornkloc(3) = zkpt(ndde(nd)) - zkpt(ndds(nd))
                else
                  lpartofdisp = .false.
                endif
              endif
              if (lpartofdisp) then
                if (lrigid) then
                  call nagammacrigid(mtvrptr,rfmass,bornkloc,maxd2,derv2,dervi,xkt,ykt,zkt)
                else
                  call nagammac(nphonatc,nphonatc,nphonatptr,nphonatptr,rfmass,bornkloc,maxd2,derv2,dervi,xkt,ykt,zkt)
                endif
              else
                if (lrigid) then
                  call nagammacrigid(mtvrptr,rfmass,bornk(1,ncf),maxd2,derv2,dervi,xkt,ykt,zkt)
                else
                  call nagammac(nphonatc,nphonatc,nphonatptr,nphonatptr,rfmass,bornk(1,ncf),maxd2,derv2,dervi,xkt,ykt,zkt)
                endif
              endif
            endif
!
!  Calculate eigenvalues and eigenvectors of corrected dynamical matrix
!
            call pdiagd(mtv,mtvloc,maxd2,derv2,dervi,maxeigc,eigc,freq(1,k-nlkpt+1),fscale,.true.,lprint,ifail)
          elseif (nbornstep(ncf).ne.0) then
!******************************************************
!  Polycrystalline average of gamma point correction  *
!******************************************************
!
!  Allocate extra storage
!
            allocate(averfreq(mtv),stat=status)
            if (status/=0) call outofmemory('phonond','averfreq')
!
!  Store derv2
!
            savederv2(1:mtv,1:mtvloc) = derv2(1:mtv,1:mtvloc)
            savedervi(1:mtv,1:mtvloc) = dervi(1:mtv,1:mtvloc)
            averfreq(1:mtv) = 0.0_dp
!
!  Set theta and phi step sizes
!
            weightpt = 1.0_dp/dble(nbornstep(ncf))
            phistep = 0.5_dp*pi*weightpt
            thetastep = phistep
!
!  Loop over theta and phi in nbornstep divisions per angle
!
            weightpt = 0.0_dp
            phi = - phistep
            do np = 0,nbornstep(ncf)
              phi = phi + phistep
              theta = - thetastep
              if (np.eq.0) then
                ntmax = 1
              elseif (np.eq.nbornstep(ncf)) then
                ntmax = 2*nbornstep(ncf)
              else
                ntmax = 4*nbornstep(ncf)
              endif
              do nt = 1,ntmax
                theta = theta + thetastep
                weightpt = weightpt + 1.0_dp
!
!  Find current Q direction of approach
!
                qfrac(1) = sin(phi)*cos(theta)
                qfrac(2) = sin(phi)*sin(theta)
                qfrac(3) = cos(phi)
!
!  Restore derv2
!
                derv2(1:mtv,1:mtvloc) = savederv2(1:mtv,1:mtvloc)
                dervi(1:mtv,1:mtvloc) = savedervi(1:mtv,1:mtvloc)
!
!  Non-analytic correction for gamma point
!
                if (lrigid) then
                  call nagammacrigid(mtvrptr,rfmass,qfrac,maxd2,derv2,dervi,xkt,ykt,zkt)
                else
                  call nagammac(nphonatc,nphonatc,nphonatptr,nphonatptr,rfmass,qfrac,maxd2,derv2,dervi,xkt,ykt,zkt)
                endif
!
!  Calculate eigenvalues and eigenvectors of corrected dynamical matrix
!
                call pdiagd(mtv,mtvloc,maxd2,derv2,dervi,maxeigc,eigc,freq(1,k-nlkpt+1),fscale,.true.,lprint,ifail)
!
!  Add frequencies to average for appropriate mode
!
                do nf = 1,mtv
                  averfreq(nf) = averfreq(nf) + freq(nf,k-nlkpt+1)
                enddo
!
!  End loops over angles
!
              enddo
            enddo
!
!  Move average frequencies back to main frequency array
!
            freq(1:mtv,k-nlkpt+1) = averfreq(1:mtv)/weightpt
!
            deallocate(averfreq,stat=status)
            if (status/=0) call deallocate_error('phonond','averfreq')
          endif
!
!  Calculate the oscillator strengths
!
          call oscillatorstrengthc(mtv,mtvloc,mtvptr,nphonatc,nphonatptr,nphonatc,nphonatptr, &
                                   eigc,maxeigc,oscstrength)
!
!  Output frequencies / DOS / intensities
!
          call peigend(mtv,mtvloc,mtvptr,mtvnptr,mtvrptr,freq(1,k-nlkpt+1),ncore,nphonatc,nphonats, &
                       nphonatptr,maxeigc,eigc,oscstrength,IRintensity(1,k-nlkpt+1),wk)
!
!  Compute group velocitiies
!
          if (lgroupvelocity) then
            call groupvelocitiesd(k-nlkpt+1,mtv,mtvloc,mtvptr,eigc,maxeigc,fscale,lprint)
          endif
!
!  Compute Grueneisen parameters
!
          if (lgrueneisen) then
            call gruneisend(k-nlkpt+1_i4,mtv,mtvloc,mtvptr,msv,nphonatc,eigc,maxeigc,derv2,dervi,maxd2,fscale,lprint)
          endif
        else
!
!  Calculate eigenvalues and eigenvectors of dynamical matrix
!
          call pdiagd(mtv,mtvloc,maxd2,derv2,dervi,maxeigc,eigc,freq(1,k-nlkpt+1),fscale,.true.,lprint,ifail)
!
!  Calculate the oscillator strengths
!
          call oscillatorstrengthc(mtv,mtvloc,mtvptr,nphonatc,nphonatptr,nphonatc,nphonatptr, &
                                   eigc,maxeigc,oscstrength)
!
!  Output frequencies / DOS / intensities
!
          call peigend(mtv,mtvloc,mtvptr,mtvnptr,mtvrptr,freq(1,k-nlkpt+1),ncore,nphonatc,nphonats, &
                       nphonatptr,maxeigc,eigc,oscstrength,IRintensity(1,k-nlkpt+1),wk)
!
!  Compute group velocitiies
!
          if (lgroupvelocity) then
            call groupvelocitiesd(k-nlkpt+1,mtv,mtvloc,mtvptr,eigc,maxeigc,fscale,lprint)
          endif
!
!  Compute Grueneisen parameters
!
          if (lgrueneisen) then
            call gruneisend(k-nlkpt+1_i4,mtv,mtvloc,mtvptr,msv,nphonatc,eigc,maxeigc,derv2,dervi,maxd2,fscale,lprint)
          endif
        endif
      endif
    else
!*********************
!  Eigenvalues only  *
!*********************
      if (lgamma.or.lintegerKpoint) then
!
!  Non-analytic correction for gamma point
!
        if (lnonanal) then
!
!  Find if this point is part of a dispersion curve and if it is substitute direction of approach
!
          lpartofdisp = (nudpt.ge.nldpt)
          if (lpartofdisp) then
            if (k.ge.ndds(nldpt).and.k.le.ndde(nudpt)) then
              lfound = .false.
              nd = nldpt - 1
              do while (.not.lfound.and.nd.lt.nudpt)
                nd = nd + 1
                lfound = (k.ge.ndds(nd).and.k.le.ndde(nd))
              enddo
              bornkloc(1) = xkpt(ndde(nd)) - xkpt(ndds(nd))
              bornkloc(2) = ykpt(ndde(nd)) - ykpt(ndds(nd))
              bornkloc(3) = zkpt(ndde(nd)) - zkpt(ndds(nd))
            else
              lpartofdisp = .false. 
            endif
          endif
          if (lintegerKpoint) then
            if (lpartofdisp) then
              if (lrigid) then
                call nagammacrigid(mtvrptr,rfmass,bornkloc,maxd2,derv2,dervi,xkt,ykt,zkt)
              else
                call nagammac(nphonatc,nphonatc,nphonatptr,nphonatptr,rfmass,bornkloc,maxd2,derv2,dervi,xkt,ykt,zkt)
              endif
            else
              if (lrigid) then
                call nagammacrigid(mtvrptr,rfmass,bornk(1,ncf),maxd2,derv2,dervi,xkt,ykt,zkt)
              else
                call nagammac(nphonatc,nphonatc,nphonatptr,nphonatptr,rfmass,bornk(1,ncf),maxd2,derv2,dervi,xkt,ykt,zkt)
              endif
            endif
          else
            if (lpartofdisp) then
              if (lrigid) then
                call nagammarigid(mtvrptr,rfmass,bornkloc,maxd2,derv2)
              else
                call nagamma(nphonatc,nphonatc,nphonatptr,rfmass,bornkloc,maxd2,derv2)
              endif
            else
              if (lrigid) then
                call nagammarigid(mtvrptr,rfmass,bornk(1,ncf),maxd2,derv2)
              else
                call nagamma(nphonatc,nphonatc,nphonatptr,rfmass,bornk(1,ncf),maxd2,derv2)
              endif
            endif
          endif
        endif
!
!  Calculate eigenvalues of uncorrected dynamical matrix
!
        if (lintegerKpoint) then
          call pdiagd(mtv,mtvloc,maxd2,derv2,dervi,maxeigc,eigc,freq(1,k-nlkpt+1),fscale,.false.,lprint,ifail)
        else
          call pdiaggd(mtv,mtvloc,maxd2,derv2,maxeigr,eigr,freq(1,k-nlkpt+1),fscale,.false.,lprint,ifail)
        endif
      else
!
!  Calculate eigenvalues of dynamical matrix
!
        call pdiagd(mtv,mtvloc,maxd2,derv2,dervi,maxeigc,eigc,freq(1,k-nlkpt+1),fscale,.false.,lprint,ifail)
      endif
!
!  Output frequencies
!
      if (lfreqout.and.lprinloc) then
        write(ioout,'(/,''  Frequencies (cm-1) [NB: Negative implies an imaginary mode]:'',/)')
        write(ioout,'(9f8.2)')(freq(i,k-nlkpt+1),i=1,mtv)
        write(ioout,'(/)')
      endif
    endif
!
!  Store frequencies if needed for output options
!
    if (ldendisp.or.(lfrq.and.lfrqbin).or.ltemprop.or.(lprinloc.and.(ntemperatureramp.gt.0))) then
      t1 = g_cpu_time()
      if (lfrq.and.ioproc) then
        if (lfrqbin) then
          do i = 1,mtv
            write(51) freq(i,k-nlkpt+1)
          enddo
        else
          do i = 1,mtv
            write(52,"(f" // adjustl(fstring) // ")") freq(i,k-nlkpt+1)
          enddo
        endif
      endif
      t2 = g_cpu_time()
      tdisk = tdisk + t2 - t1
    elseif (lfrq.and.ioproc) then
      t1 = g_cpu_time()
      do i = 1,mtv
        write(52,'(f12.6)') freq(i,k-nlkpt+1)
      enddo
      t2 = g_cpu_time()
      tdisk = tdisk + t2 - t1
    endif
!***************************************
!  Evaluate phonon related properties  *
!***************************************
    if (lmeanke) then
      if (temperature.gt.1.0d-6) then
!
!  Scale frequencies to hw/kT
!
        cmfact = planck*speedl/rkt
        factor = 0.5_dp*wk*rkt/evtoj
!
!  The following is equivalent to looping over mtvloc
!
        do iloc = 1,nphonatonnodec
          iatm = nphonatonnodecptr(iloc)
          i = 3*(iatm-1)
          do ii = 1,3
            i = i + 1
            frq = cmfact*freq(i,k-nlkpt+1)
!
!  Store exp(x) in w1, exp(x)-1 in w2 and 1/(exp(x)-1) in w3
!
            if (frq.lt.12.0_dp) then
              w1 = exp(frq)
              w2 = w1 - 1.0_dp
              if (abs(w2).gt.0.0_dp) w3 = 1.0_dp/w2
            else
              w3 = exp(-frq)
            endif
!
!  Aportion kinetic energy to atoms based on eigenvector contribution
!
            trmke = factor*frq*(0.5_dp + w3)
            if (lgamma) then
              do j = 1,mtv
                jj = (j-1)/3_i4 + 1
                projectionfactor = eigr(j,iloc)**2
                meanKEperatom(jj) = meanKEperatom(jj) + trmke*projectionfactor
              enddo
            else
              do j = 1,mtv
                jj = (j-1)/3_i4 + 1
                projectionfactor = dble(conjg(eigc(j,iloc))*eigc(j,iloc))
                meanKEperatom(jj) = meanKEperatom(jj) + trmke*projectionfactor
              enddo
            endif
          enddo
        enddo
!
!  Global reduction of array
!
        allocate(sum(mtv),stat=status)
        if (status/=0) call outofmemory('phonond','sum')
!
        call sumall(meanKEperatom,sum,mtv,"phonond","meanKEperatom")
        do i = 1,mtv
          meanKEperatom(i) = sum(i)
        enddo
!
        deallocate(sum,stat=status)
        if (status/=0) call deallocate_error('phonond','sum')
      endif
    endif
!
!  Output eigenvectors if requested
!
    if (lfreqout.and.lprinloc) then
      write(ioout,'(''--------------------------------------------------------------------------------'')')
    endif
!
!  Frequency-dependent properties
!
    if (lomega(ncf).and.lgamma) then
      call omegaproperty(mtv,savefreq,fscale,oscstrength)
    endif
!***************************
!  End loop over K points  *
!***************************
  enddo
  if (lpdfloc) then
!
!  Set up variables for gulp_cml_PDFstats
!
    call pdfsetup(lcml) 
  endif
!**************************************************
!  Deallocate second derivative workspace memory  *
!**************************************************
  if (lgammaonly) then
    if (allocated(savederv2)) then
      deallocate(savederv2,stat=status)
      if (status/=0) call deallocate_error('phonond','savederv2')
    endif
    deallocate(eigr,stat=status)
    if (status/=0) call deallocate_error('phonond','eigr')
  else
    if (allocated(savederv2)) then
      deallocate(savedervi,stat=status)
      if (status/=0) call deallocate_error('phonond','savedervi')
      deallocate(savederv2,stat=status)
      if (status/=0) call deallocate_error('phonond','savederv2')
    endif
    deallocate(eigc,stat=status)
    if (status/=0) call deallocate_error('phonond','eigc')
    if (allocated(eigr)) then
      deallocate(eigr,stat=status)
      if (status/=0) call deallocate_error('phonond','eigr')
    endif
  endif
!*************************************
!  Output phonon related properties  *
!*************************************
!
!  Compute vibrational energies and output
!
  call vibenergy(mtv,nllkpt,wkpt(nlkpt),freq(1,1),maxfreq,rnokpt,lprint,fc)
!
!  Mean squared displacements
!
  if (lmsd) then
!
!  Global reduction of array
!
    allocate(sum(3*ncfoc),stat=status)
    if (status/=0) call outofmemory('phonond','sum')
    allocate(sum2(3*ncfoc),stat=status)
    if (status/=0) call outofmemory('phonond','sum2')
!
    sum(1:ncfoc) = msdx(1:ncfoc)
    sum(ncfoc+1:2*ncfoc) = msdy(1:ncfoc)
    sum(2*ncfoc+1:3*ncfoc) = msdz(1:ncfoc)
!
    call sumall(sum,sum2,3_i4*ncfoc,"phonond","msd")
!
    msdx(1:ncfoc) = sum2(1:ncfoc)
    msdy(1:ncfoc) = sum2(ncfoc+1:2*ncfoc)
    msdz(1:ncfoc) = sum2(2*ncfoc+1:3*ncfoc)
!
    deallocate(sum2,stat=status)
    if (status/=0) call deallocate_error('phonond','sum2')
    deallocate(sum,stat=status)
    if (status/=0) call deallocate_error('phonond','sum')
!
    if (lprinloc) then
      write(ioout,'(/,''  Mean-squared displacements for phonons (Ang**2) : '',/)')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''   Atom '',8x,''x'',12x,''y'',12x,''z'')')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      if (lrigid) then
        do i = 1,numatnomol
          write(ioout,'(i8,5x,3(f12.6,1x))') i,msdx(i),msdy(i),msdz(i)
        enddo
        write(ioout,'(''--------------------------------------------------------------------------------'')')
        write(ioout,'(''   Molecule '',8x,''x'',12x,''y'',12x,''z'')')
        write(ioout,'(''--------------------------------------------------------------------------------'')')
        do i = 1,nmol
          j = numatnomol + i
          write(ioout,'(i8,5x,3(f12.6,1x))') i,msdx(j),msdy(j),msdz(j)
        enddo
      else
        do i = 1,ncfoc
          write(ioout,'(i8,1x,3(f12.6,1x))') i,msdx(i),msdy(i),msdz(i)
        enddo
      endif
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(/)')
    endif
  endif
!
!  Mean kinetic energies
!
  if (lprinloc.and.lmeanke) then
    write(ioout,'(''  Mean kinetic energy per site: (eV)'')')
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    do i = 1,nphonatc
      write(ioout,'(i6,2x,f15.6)') i,meanKEperatom(i)
    enddo
    write(ioout,'(''--------------------------------------------------------------------------------'')')
  endif
!
!  If output if required then call outphon
!
  if (ldendisp.and.ioproc) then
    call outphon(nlkpt,nukpt,leigloc,mtv)
    if (lfrq.and.lfrqbin) close(51)
  elseif (lfrq.and.lfrqbin.and.ioproc) then
    close(51)
  endif
!***************
!  Exit point  *
!***************
999 continue
!
!  Close I/O channels
!
  if (leig.and.ioproc) then
    close(53)
  endif
  if (lcas.and.ioproc) then
    close(54)
  endif
!
!  CML output
!
  if (lcmlloc) call gulp_cml_endPhonons
  if (lcml.and.lpdfloc) then
    call gulp_cml_PDFstats
    call gulp_cml_outPDF
  endif
!
!  Free local memory
!
  if (lthermal) then
    deallocate(Sij,stat=status)
    if (status/=0) call deallocate_error('phonond','Sij')
  endif
  if (leigloc) then
    deallocate(ramstrength,stat=status)
    if (status/=0) call deallocate_error('phonond','ramstrength')
    deallocate(oscstrength,stat=status)
    if (status/=0) call deallocate_error('phonond','oscstrength')
    deallocate(savefreq,stat=status)
    if (status/=0) call deallocate_error('phonond','savefreq')
  endif
  deallocate(mtvrptr,stat=status)
  if (status/=0) call deallocate_error('phonond','mtvrptr')
  deallocate(mtvnptr,stat=status)
  if (status/=0) call deallocate_error('phonond','mtvnptr')
  deallocate(mtvptr,stat=status)
  if (status/=0) call deallocate_error('phonond','mtvptr')
  if (lrigid.and.leckart) then
    deallocate(mcvptr,stat=status)
    if (status/=0) call deallocate_error('phonond','mcvptr')
  endif
  deallocate(meanKEperatom,stat=status)
  if (status/=0) call deallocate_error('phonond','meanKEperatom')
!
!  Close PDF related parts of phonon 
!
  call closepdfphonon
  if (lshengBTE) then
!
!  Output shengBTE control file
!
    call outshengBTEcontrol(iout_shengBTE,.false.)
!
!  Output shengBTE second derivative file
!
    call outshengBTEfc2(iout_shengBTE)
  endif
!
!  Timing
!
  t2t = g_cpu_time()
  tphon = t2t - t1t + tphon
#ifdef TRACE
  call trace_out('phonond')
#endif
#else
  call outerror('phonond called when not compiled with MPI',0_i4)
  call stopnow('phonond')
#endif
!
  return
  end
