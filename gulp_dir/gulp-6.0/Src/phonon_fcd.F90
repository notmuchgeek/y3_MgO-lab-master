  subroutine phonon_fcd(lprint,fc,nobsmodeptr0,nobsmode)
!
!  Calculates the phonons at a given set of k points.
!  This version uses a different algorithm in which the second
!  derivatives are stored by cell image and only calculated once.
!  The phasing for k points is then applied later.
!  Distributed second derivative version.
!
!  The block diagonal elements of derv2 must be stored, as these are
!  not recalculated in dynamic as the summing of off diagonals is no
!  longer applicable.
!
!  NB: Currently doesn't support partial occupancy!
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
!   4/17 Created from phonon_fc
!   5/17 Calls to matrix inversion modified to use nmin.ne.1
!   5/17 Parallel thermal conductivity call added
!   5/17 nshell and nshellonnode added to subroutine call
!   7/17 nshell and nshellonnode added to complex subroutine call
!   7/17 outfrc now called in parallel
!   7/17 Call to getvibmode replaced by getvibmoded
!   8/17 Routine wrapped in ifdef since it won't work or be called without MPI
!   8/17 Correction for region 2 phonon case
!   8/17 Modifications to accommodate Intel MPI and scalapack
!  11/17 leigloc set to true for MSD calc
!  11/17 wk now passed to peigend
!   1/18 Trace added
!   3/18 dynam shortened to dyna in keyword check
!   3/18 Parallel I/O corrected
!   4/18 leigloc now set here for llower and linten
!   3/19 Multiple temperature ramps added
!   5/19 Finite difference flag split for first and second derivatives
!   8/19 Call to peigend corrected
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!   4/20 maxfqat changed to maxfreq since it is no longer just atom based
!   4/20 Mass arrays now use fmass and rfmass
!   5/20 ncfoc removed from arguments to groupvelocities
!   6/20 Argument list to peigen routines changed
!   6/20 Last argument to outphon changed
!   6/20 mcvptr passed to groupvelocitiesd
!  12/20 Check added for gamma point only for lused1 with EAM/MEAM
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
  use configurations
  use g_constants
  use control
  use current
  use derivatives
  use dispersion
  use element
  use gulp_files
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
#ifdef MPI
  use sutton,          only : lsuttonc
#endif
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
  integer(i4)                                       :: iloc
  integer(i4)                                       :: ind
  integer(i4)                                       :: indi
  integer(i4)                                       :: indiloc
  integer(i4)                                       :: indii
  integer(i4)                                       :: indj
  integer(i4)                                       :: indjj
  integer(i4)                                       :: j
  integer(i4)                                       :: jj
  integer(i4)                                       :: k
  integer(i4)                                       :: kk
  integer(i4),  dimension(:),     allocatable       :: kpvt
  integer(i4)                                       :: m
  integer(i4)                                       :: mm
  integer(i4)                                       :: maxeigc
  integer(i4)                                       :: maxeigr
  integer(i4)                                       :: maxlim
  integer(i4)                                       :: mcv
  integer(i4)                                       :: mcvloc
  integer(i4),  dimension(:),     allocatable       :: mcvptr
  integer(i4),  dimension(:),     allocatable       :: mcvnptr
  integer(i4),  dimension(:),     allocatable       :: mcvrptr
  integer(i4)                                       :: mint
  integer(i4)                                       :: mintloc
  integer(i4)                                       :: msv
  integer(i4)                                       :: msvloc
  integer(i4)                                       :: ncind
  integer(i4)                                       :: nd
  integer(i4)                                       :: nf
  integer(i4)                                       :: nfitmode
  integer(i4)                                       :: nk
  integer(i4)                                       :: nldpt
  integer(i4)                                       :: nlkpt
  integer(i4)                                       :: nllkpt
  integer(i4)                                       :: nobm
  integer(i4)                                       :: np
  integer(i4)                                       :: np_nlkpt
  integer(i4)                                       :: np_procs
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
  integer                                           :: ntag
  integer                                           :: nnode
  integer                                           :: ntmp
  integer                                           :: Request
  integer(i4),      dimension(:),   allocatable     :: StatMPI       ! Array for status from MPI
!
#ifdef INTEL
  logical                                           :: leigcOK
#endif
  logical                                           :: lcmlloc
  logical                                           :: ldendisp
  logical                                           :: leigloc
  logical                                           :: lfound
  logical                                           :: lgamma
  logical                                           :: lgammaonly
  logical                                           :: lintegerKpoint
  logical                                           :: lintegerKpoint_present
  logical                                           :: lnonanal
  logical                                           :: lnoanald2loc
  logical                                           :: lpartofdisp
  logical                                           :: lprinloc
  logical                                           :: ltemprop
  logical                                           :: lused1
  complex(dpc), dimension(:,:),   allocatable       :: eigc
  real(dp),     dimension(:),     allocatable       :: averfreq
  real(dp)                                          :: bornkloc(3)
  real(dp)                                          :: cmfact
  real(dp),     dimension(:,:),   allocatable       :: dtmp
  real(dp),     dimension(:,:),   allocatable       :: eigr
  real(dp)                                          :: g_cpu_time
  real(dp)                                          :: factor
  real(dp)                                          :: frq
  real(dp)                                          :: fscale
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
#ifdef TRACE
  call trace_in('phonon_fcd')
#endif
!
  t1t = g_cpu_time()
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
  if (lthermal) leigloc = .true.
  if (lgroupvelocity) leigloc = .true.
  if (linten) leigloc = .true.
  if (llower) leigloc = .true.
  lprinloc = (lprint.and.ioproc)
  ltemprop = (temperature.gt.1.0d-6.or.ntemperatureramp.gt.0)
  lnonanal = (index(keyword,'nono').eq.0.and.ndim.eq.3.and..not.leem.and..not.lnoanald2)
  lnoanald2loc = lnoanald2
!
!  Set flag as to whether to use first derivative algorithm or second
!
  lused1 = (lnoanald2loc.or.lfinitediff2)
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
  if (lshengBTE) then
!
!  Output shengBTE control file
!
    call outshengBTEcontrol(iout_shengBTE,.true.)
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
    call trace_out('phonon_fcd')
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
  lpocc = (nsfoc+ncfoc.ne.nphonat)
  if (lpocc) then
    call outerror('partial occupancy not compatible with distributed 2nd derivs',0_i4)
    call stopnow('phonon_fcd')
  endif
!
!  Allocate local pointer arrays
!
  allocate(meanKEperatom(numat),stat=status)
  if (status/=0) call outofmemory('phonon_fcd','meanKEperatom')
!
!  Calculate a few constants to do with the size of the problem
!
  mint = 3*nphonat
  mintloc = 3*nphonatonnode
!
  maxlim = mint
  if (nbsmat.gt.0) maxlim = maxlim + nphonat
!
  msv = 3*nphonats + nphonatb
  mcv = 3*nphonatc
  msvloc = 3*nphonatonnodes + nphonatonnodeb
  mcvloc = 3*nphonatonnodec
!
  allocate(mcvptr(mcvloc),stat=status)
  if (status/=0) call outofmemory('phonon_fcd','mcvptr')
  allocate(mcvnptr(mcv),stat=status)
  if (status/=0) call outofmemory('phonon_fcd','mcvnptr')
  allocate(mcvrptr(mcv),stat=status)
  if (status/=0) call outofmemory('phonon_fcd','mcvrptr')
!
  mcvrptr(1:mcv) = 0
  do iloc = 1,nphonatonnodec
    i = nphonatonnodecptr(iloc)
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
  do i = 1,nphonat
    mm = 3*(nphonatptr(i)-1)
    mcvnptr(mm+1) = atom2node(nphonatptr(i))
    mcvnptr(mm+2) = atom2node(nphonatptr(i))
    mcvnptr(mm+3) = atom2node(nphonatptr(i))
  enddo
!
!  Check that maxd2 is greater than or equal to mcv
!
  if (maxd2.lt.mcv+msv) then
    maxd2 = mcv + msv
    call changemaxd2
  endif
  if (maxd2u.lt.mcvloc+msvloc) then
    maxd2u = mcvloc + msvloc
    call changemaxd2
  endif
!
!  Check number of supercells
!
  if (maxd2cells.lt.nd2cells) then
    maxd2cells = nd2cells
    call changemaxd2cells
  endif
  if (maxd2cell.lt.nd2cell(1)) then
    maxd2cell = nd2cell(1)
    call changemaxd2cell
  endif
  if (maxd2cell.lt.nd2cell(2)) then
    maxd2cell = nd2cell(2)
    call changemaxd2cell
  endif
  if (maxd2cell.lt.nd2cell(3)) then
    maxd2cell = nd2cell(3)
    call changemaxd2cell
  endif
!
!  Build index pointer
!
  ncind = 0
  do ii = -nd2cell(1),nd2cell(1)
    do jj = -nd2cell(2),nd2cell(2)
      do kk = -nd2cell(3),nd2cell(3)
        ncind = ncind + 1
        nd2cellptr(nd2cell(1)+1+ii,nd2cell(2)+1+jj,nd2cell(3)+1+kk) = ncind
      enddo
    enddo
  enddo
!
!  Set central cell index
!
  nd2central = (nd2cells + 1)/2
!
!  Ensure that frequency array is large enough
!
  call changemaxfreq(3_i4*nphonatc)
!
!  Allocate array to hold oscillator strengths if needed
!
  if (leigloc) then
    allocate(savefreq(mcv),stat=status)
    if (status/=0) call outofmemory('phonon_fcd','savefreq')
    allocate(oscstrength(3,3,mcv),stat=status)
    if (status/=0) call outofmemory('phonon_fcd','oscstrength')
    allocate(ramstrength(3,3,mcv),stat=status)
    if (status/=0) call outofmemory('phonon_fcd','ramstrength')
  endif
!
!  Allocate array for thermal conductivity if needed
!
  if (lthermal) then
    allocate(Sij(mcv,mcvloc),stat=status)
    if (status/=0) call outofmemory('phonon_fcd','Sij')
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
  do i = 1,3*nphonatc
    fmass(i) = 0.0_dp
  enddo
  do i = 1,nphonatc
    nsi = nspecptr(nrelf2a(nphonatptr(i)))
    rmassi = massspec(nsi)*occuf(nphonatptr(i))
    if (natspec(nsi).le.maxele.and.abs(rmassi).lt.1.0d-12) then
      call label(natspec(nsi),ntypspec(nsi),lab1)
      call outerror('mass of species '//lab1//' is zero',0_i4)
      goto 999
    endif
    fmass(3*(i-1)+1) = fmass(3*(i-1)+1) + rmassi
    fmass(3*(i-1)+2) = fmass(3*(i-1)+2) + rmassi
    fmass(3*(i-1)+3) = fmass(3*(i-1)+3) + rmassi
  enddo
  do i = 1,3*nphonatc
    if (fmass(i).eq.0.0_dp) then
      call outerror('site has total mass of zero in phonon_fcd',0_i4)
      call stopnow('phonon_fcd')
    endif
    rfmass(i) = 1.0_dp/sqrt(fmass(i))
  enddo
!
!  Sum up k point weights and check for gamma point
!
  lgamma = .false.
  lintegerKpoint = .false.
  sumwkpt = 0.0_dp
  do i = nlkpt,nukpt
    sumwkpt = sumwkpt + wkpt(i)
    if (xkpt(i).eq.0.0_dp.and.ykpt(i).eq.0.0_dp.and.zkpt(i).eq.0.0_dp) then
      lgamma = .true.
    else
      xmod = mod(xkpt(i)+100.0_dp,2.0_dp)
      ymod = mod(ykpt(i)+100.0_dp,2.0_dp)
      zmod = mod(zkpt(i)+100.0_dp,2.0_dp)
      if (xmod+ymod+zmod.lt.1.0d-10) then
        lintegerKpoint_present = .true.
      endif
    endif
  enddo
!
!  Is the gamma point the only K point?
!
  lgammaonly = (lgamma.and.nllkpt.eq.1)
!
!  Check algorithm works for k points
!
  if (.not.lgammaonly.and.lused1.and.lsuttonc) then
    call outerror('phonon algorithm only available for gamma point when using EAM/MEAM',0_i4)
    call stopnow('phonon_fcd')
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
    write(ioout,'(/,''  Phonon Calculation : '',/)')
    if ((lgamma.or.lintegerKpoint).and.lnonanal) then
      if (nbornstep(ncf).gt.0) then
        write(ioout,'(''  Number of angular steps for n-a correction at gamma = '',i6,/)') nbornstep(ncf)
      else
        write(ioout,'(''  K direction for n-a correction at gamma ='',3(1x,f8.5),/)') bornk(1,ncf),bornk(2,ncf),bornk(3,ncf)
      endif
    endif
    if ((lgamma.or.lintegerKpoint).and.lraman) then
      if (ramandirtype(ncf).eq.1) then
        write(ioout,'(''  Directions for Raman intensities: Type  =  Cartesian'')')
      else
        write(ioout,'(''  Directions for Raman intensities: Type  =  Fractional '')')
      endif
      write(ioout,'(''                                    In    ='',3(1x,f8.5))') ramandir(1,ncf),ramandir(2,ncf),ramandir(3,ncf)
      write(ioout,'(''                                    Out   ='',3(1x,f8.5),/)') ramandir(4,ncf),ramandir(5,ncf),ramandir(6,ncf)
    endif
    write(ioout,'(''  Number of k points for this configuration = '',i8,/)')(nukpt-nlkpt+1)
    write(ioout,'(''--------------------------------------------------------------------------------'')')
  endif
!****************************************************
!  If frequencies are to be written to a permanent  *
!  file then open with the correct type here        *
!****************************************************
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
      allocate(eigr(maxeigr,mcvloc),stat=status)
      if (status/=0) call outofmemory('phonon_fcd','eigr')
    else
      if (lgamma.or.lintegerKpoint) then
        allocate(eigr(maxeigr,mcvloc),stat=status)
        if (status/=0) call outofmemory('phonon_fcd','eigr')
      endif
      allocate(eigc(maxeigc,mcvloc+msvloc),stat=status)
      if (status/=0) call outofmemory('phonon_fcd','eigc')
    endif
  else
    if (msv.gt.0) then
      maxeigc = mcv + msv
      maxeigr = msv
      if (lgammaonly) then
        allocate(eigr(maxeigr,msv),stat=status)
        if (status/=0) call outofmemory('phonon_fcd','eigr')
      else
        if (lgamma.or.lintegerKpoint) then
          allocate(eigr(maxeigr,msv),stat=status)
          if (status/=0) call outofmemory('phonon_fcd','eigr')
        endif
        allocate(eigc(maxeigc,mcvloc+msvloc),stat=status)
        if (status/=0) call outofmemory('phonon_fcd','eigc')
      endif
    else
      maxeigc = 1
      maxeigr = 1
      if (lgammaonly) then
        maxeigr = mcv
        allocate(eigr(mcv,mcvloc),stat=status)
        if (status/=0) call outofmemory('phonon_fcd','eigr')
      else
        maxeigc = mcv
        if (lgamma.or.lintegerKpoint) then
          maxeigr = mcv
          allocate(eigr(mcv,mcvloc),stat=status)
          if (status/=0) call outofmemory('phonon_fcd','eigr')
        endif
        allocate(eigc(mcv,mcvloc),stat=status)
        if (status/=0) call outofmemory('phonon_fcd','eigc')
      endif
    endif
  endif
!***********************
!  Setup PDF Variables *
!***********************
  if (lpdf) then
    call setuppdfphonon(nllkpt,mcv,nphonatc,nphonatptr,ncf)
  endif
!*************************************************************
!  Compute second derivatives once and store per cell image  *
!*************************************************************
  if (lused1) then
    call dynamicn_fcd
  else
    call dynamic_fc
  endif
!
!  Set diagonal elements
!
  call diagonal_fcd(lused1)
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
  if (lshengBTE) then
!
!  Output shengBTE second derivative file using supercell format
!
    call outshengBTEfc2s(iout_shengBTE)
  endif
!******************************
!  Initialise MSDs if needed  *
!******************************
  if (lmsd) then
    msdx(1:ncfoc) = 0.0_dp
    msdy(1:ncfoc) = 0.0_dp
    msdz(1:ncfoc) = 0.0_dp
  endif
!***********************
!  Loop over k points  *
!***********************
  np_nlkpt = nlkpt
  np_procs = 1_i4
  do k = np_nlkpt,nukpt,np_procs
    if (lpdf) then
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
    lgamma = (abs(xkt)+abs(ykt)+abs(zkt).lt.1.0d-10)
    if (.not.lgamma) then
      xmod = mod(xkt+100.0_dp,2.0_dp)
      ymod = mod(ykt+100.0_dp,2.0_dp)
      zmod = mod(zkt+100.0_dp,2.0_dp)
      lintegerKpoint = (xmod+ymod+zmod.lt.1.0d-10)
    else
      lintegerKpoint = .false.
    endif
!
    if (lfreqout.and.lprinloc) then
      write(ioout,'(''  K point '',i6,'' = '',3f10.6,''  Weight = '',f8.3)') k,xkt,ykt,zkt,wk
      write(ioout,'(''--------------------------------------------------------------------------------'')')
    endif
!
!  Write K point headers to files if required
!
    if (leig.and.ioproc) then
      write(53,'(''K point at '',3f10.6,'' in BZ '')') xkt,ykt,zkt
    endif
    if (lcas.and.ioproc) then
      write(54,'(a10,i5,3f12.6,f14.6)') "     q-pt=",k,xkt,ykt,zkt,wk
    endif
!
!  Generate phased second derivatives
!
    call phase_fcd(xkt,ykt,zkt)
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
      endif
    enddo
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
        ntmp = 3*(mcv + msv)
        ntag = 1
        allocate(dtmp(mcv+msv,3_i4),stat=status)
        if (status/=0) call outofmemory('phonon_fcd','dtmp')
        allocate(StatMPI(MPI_Status_Size),stat=status)
        if (status/=0) call outofmemory('phonon_fcd','StatMPI')
      endif
      call mpbarrier
      do iloc = 1,nphonat
        i = nphonatonnoderptr(iloc)
        if (lioproconly.and.mcvnptr(i).ne.0_i4) then
!
!  Post receive
!
          if (ioproc) then
            nnode = mcvnptr(i)
            call MPI_IRecv(dtmp,ntmp,MPI_double_precision,nnode, &
                           ntag,MPI_Comm_World,Request,MPIerror)
          endif
!
!  Pass data to ioproc for writing
!
          if (i.gt.0) then
            indi = 3*(i-1)
            dtmp(1:mcv+msv,1) = derv2(1:mcv+msv,indi+1)
            dtmp(1:mcv+msv,2) = derv2(1:mcv+msv,indi+2)
            dtmp(1:mcv+msv,3) = derv2(1:mcv+msv,indi+3)
!
!  Post send
!
            call MPI_ISend(dtmp,ntmp,MPI_double_precision,0, &
                           ntag,MPI_Comm_World,Request,MPIerror)
          endif
          if (ioproc.or.i.gt.0) then
            call MPI_WaitAll(1,Request,StatMPI,MPIerror)
          endif
          if (ioproc) then
!
!  Write on I/O node
!
            write(ioout,'(12f11.6)')(dtmp(j,1),j=1,mcv+msv)
            write(ioout,'(12f11.6)')(dtmp(j,2),j=1,mcv+msv)
            write(ioout,'(12f11.6)')(dtmp(j,3),j=1,mcv+msv)
          endif
        else
          if (i.gt.0) then
            indi = 3*(i-1)
            write(ioout,'(12f11.6)')(derv2(j,indi+1),j=1,mcv+msv)
            write(ioout,'(12f11.6)')(derv2(j,indi+2),j=1,mcv+msv)
            write(ioout,'(12f11.6)')(derv2(j,indi+3),j=1,mcv+msv)
          endif
        endif
        call mpbarrier
      enddo
      if (.not.lgamma) then
        if (ioproc) then
          write(ioout,'(/,''  Uncompressed Imag Dynamical matrix :'',/)')
        endif
        call mpbarrier
        do iloc = 1,nphonat
          i = nphonatonnoderptr(iloc)
          if (lioproconly.and.mcvnptr(i).ne.0_i4) then
!
!  Post receive
!
            if (ioproc) then
              nnode = mcvnptr(i)
              call MPI_IRecv(dtmp,ntmp,MPI_double_precision,nnode, &
                             ntag,MPI_Comm_World,Request,MPIerror)
            endif
!
!  Pass data to ioproc for writing
!
            if (i.gt.0) then
              indi = 3*(i-1)
              dtmp(1:mcv+msv,1) = dervi(1:mcv+msv,indi+1)
              dtmp(1:mcv+msv,2) = dervi(1:mcv+msv,indi+2)
              dtmp(1:mcv+msv,3) = dervi(1:mcv+msv,indi+3)
!
!  Post send
!
              call MPI_ISend(dtmp,ntmp,MPI_double_precision,0, &
                             ntag,MPI_Comm_World,Request,MPIerror)
            endif
            if (ioproc.or.i.gt.0) then
              call MPI_WaitAll(1,Request,StatMPI,MPIerror)
            endif
          else
            if (i.gt.0) then
              indi = 3*(i-1)
              write(ioout,'(12f11.6)')(dervi(j,indi+1),j=1,mcv+msv)
              write(ioout,'(12f11.6)')(dervi(j,indi+2),j=1,mcv+msv)
              write(ioout,'(12f11.6)')(dervi(j,indi+3),j=1,mcv+msv)
            endif
          endif
          call mpbarrier
        enddo
      endif
      if (lioproconly) then
        deallocate(StatMPI,stat=status)
        if (status/=0) call deallocate_error('phonon_fcd','StatMPI')
        deallocate(dtmp,stat=status)
        if (status/=0) call deallocate_error('phonon_fcd','dtmp')
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
!
!  Call library to invert matrix stored in eigr
!
#ifdef INTEL
        leigcOK = (allocated(eigc))
        if (.not.leigcOK) then
          allocate(eigc(maxeigc,mcvloc+msvloc),stat=status)
          if (status/=0) call outofmemory('phonon_fcd','eigc')
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
          call stopnow('phonon_fcd')
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
          if (status/=0) call deallocate_error('phonon_fcd','eigc')
        endif
#else
        call matrix_inversion_shells(msv,mcv+1_i4,maxd2,derv2,nshell,nshellonnode,ifail)
!
!  Check return flag
!
        if (ifail.ne.0) then
          call outerror('inversion of shell 2nd derivatives failed',0_i4)
          call stopnow('phonon_fcd')
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
          call stopnow('phonon_fcd')
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
          call stopnow('phonon_fcd')
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
          call stopnow('phonon_fcd')
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
!  Transfer back from eigc to derv2/dervi
!
        do i = 1,mcvloc
          do j = 1,mcv
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
!  Copy core-core component of derv2dk to eigc
!
            do i = 1,mcvloc
              do j = 1,mcv
                eigc(j,i) = 0.0_dpc
              enddo
            enddo
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
                        idesc,(-1.0d0,0.0d0),eigc,1,1,idesc)
!
!  Copy back core-core component of derv2dk from eigc
!
            do i = 1,mcvloc
              do j = 1,mcv
                derv2dk(ii,j,i) = derv2dk(ii,j,i) + eigc(j,i) + conjg(eigc(j,i))
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
!  Option to write out a .frc file for QMPOT
!
    if (lfrc.and.lgamma) call outfrc(fc,.true.,.true.)
!*****************************
!  Multiply by mass-factors  *
!*****************************
    if (lgamma) then
      ind = 0
      do i = 1,nphonatonnodec
        indi = 3*(nphonatonnodecptr(i)-1)
        do ii = 1,3
          ind  = ind + 1
          indi = indi + 1
          rmassi = rfmass(indi)
          do j = 1,3*nphonatc
            derv2(j,ind) = rmassi*rfmass(j)*derv2(j,ind)
          enddo
        enddo
      enddo
    else
      ind = 0
      do i = 1,nphonatonnodec
        indi = 3*(nphonatonnodecptr(i)-1)
        do ii = 1,3
          ind  = ind + 1
          indi = indi + 1
          rmassi = rfmass(indi)
          do j = 1,3*nphonatc
            derv2(j,ind) = rmassi*rfmass(j)*derv2(j,ind)
            dervi(j,ind) = rmassi*rfmass(j)*dervi(j,ind)
          enddo
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
        ntmp = 3*mcv
        ntag = 1
        allocate(dtmp(mcv,3_i4),stat=status)
        if (status/=0) call outofmemory('phonon_fcd','dtmp')
        allocate(StatMPI(MPI_Status_Size),stat=status)
        if (status/=0) call outofmemory('phonon_fcd','StatMPI')
      endif
      call mpbarrier
      do iloc = 1,nphonatc
        i = nphonatonnodercptr(iloc)
        if (lioproconly.and.mcvnptr(i).ne.0_i4) then
!
!  Post receive
!
          if (ioproc) then
            nnode = mcvnptr(i)
            call MPI_IRecv(dtmp,ntmp,MPI_double_precision,nnode, &
                           ntag,MPI_Comm_World,Request,MPIerror)
          endif
!
!  Pass data to ioproc for writing
!
          if (i.gt.0) then
            indi = 3*(i-1)
            dtmp(1:mcv,1) = derv2(1:mcv,indi+1)
            dtmp(1:mcv,2) = derv2(1:mcv,indi+2)
            dtmp(1:mcv,3) = derv2(1:mcv,indi+3)
!
!  Post send
!
            call MPI_ISend(dtmp,ntmp,MPI_double_precision,0, &
                           ntag,MPI_Comm_World,Request,MPIerror)
          endif
          if (ioproc.or.i.gt.0) then
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
          if (i.gt.0) then
            indi = 3*(i-1)
            write(ioout,'(12f11.6)')(derv2(j,indi+1),j=1,mcv)
            write(ioout,'(12f11.6)')(derv2(j,indi+2),j=1,mcv)
            write(ioout,'(12f11.6)')(derv2(j,indi+3),j=1,mcv)
          endif
        endif
        call mpbarrier
      enddo
      if (.not.lgamma) then
        if (ioproc) then
          write(ioout,'(/,''  Imaginary Dynamical matrix :'',/)')
        endif
        call mpbarrier
        do iloc = 1,nphonatc
          i = nphonatonnodercptr(iloc)
          if (lioproconly.and.mcvnptr(i).ne.0_i4) then
!
!  Post receive
!
            if (ioproc) then
              nnode = mcvnptr(i)
              call MPI_IRecv(dtmp,ntmp,MPI_double_precision,nnode, &
                             ntag,MPI_Comm_World,Request,MPIerror)
            endif
!
!  Pass data to ioproc for writing
!
            if (i.gt.0) then
              indi = 3*(i-1)
              dtmp(1:mcv,1) = dervi(1:mcv,indi+1)
              dtmp(1:mcv,2) = dervi(1:mcv,indi+2)
              dtmp(1:mcv,3) = dervi(1:mcv,indi+3)
!
!  Post send
!
              call MPI_ISend(dtmp,ntmp,MPI_double_precision,0, &
                             ntag,MPI_Comm_World,Request,MPIerror)
            endif
            if (ioproc.or.i.gt.0) then
              call MPI_WaitAll(1,Request,StatMPI,MPIerror)
            endif
          else
            if (i.gt.0) then
              indi = 3*(i-1)
              write(ioout,'(12f11.6)')(dervi(j,indi+1),j=1,mcv)
              write(ioout,'(12f11.6)')(dervi(j,indi+2),j=1,mcv)
              write(ioout,'(12f11.6)')(dervi(j,indi+3),j=1,mcv)
            endif
          endif
          call mpbarrier
        enddo
      endif
      if (lioproconly) then
        deallocate(StatMPI,stat=status)
        if (status/=0) call deallocate_error('phonon_fcd','StatMPI')
        deallocate(dtmp,stat=status)
        if (status/=0) call deallocate_error('phonon_fcd','dtmp')
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
        call pdiaggd(mcv,mcvloc,maxd2,derv2,maxeigr,eigr,freq(1,k-nlkpt+1),fscale,.true.,lprint,ifail)
!
!  Calculate the oscillator strengths
!
        call oscillatorstrengthg(mcv,mcvloc,mcvptr,nphonatc,nphonatptr,nphonatc,iocptr, &
                                 eigr,maxeigr,oscstrength)
!
!  Optionally call Raman strengths
!
        if (lraman) then
          call ramanstrength(mcv,mcvloc,mcvptr,nphonatc,nphonatptr,nphonatc,iocptr, &
                             eigr,maxeigr,ramstrength)
        endif
!
!  Store uncorrected frequencies
!
        savefreq(1:mcv) = freq(1:mcv,k-nlkpt+1)
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
              call nagamma(nphonatc,nphonatc,iocptr,rfmass,bornkloc,maxd2,derv2)
            else
              call nagamma(nphonatc,nphonatc,iocptr,rfmass,bornk(1,ncf),maxd2,derv2)
            endif
          endif
!
!  Calculate eigenvalues and eigenvectors of corrected dynamical matrix
!
          call pdiaggd(mcv,mcvloc,maxd2,derv2,maxeigr,eigr,freq(1,k-nlkpt+1),fscale,.false.,lprint,ifail)
        elseif (nbornstep(ncf).ne.0) then
!******************************************************
!  Polycrystalline average of gamma point correction  *
!******************************************************
!
!  Allocate extra storage
!
          allocate(averfreq(mcv),stat=status)
          if (status/=0) call outofmemory('phonon_fcd','averfreq')
!
!  Store derv2
!
          savederv2(1:mcv,1:mcvloc) = derv2(1:mcv,1:mcvloc)
          averfreq(1:mcv) = 0.0_dp
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
              derv2(1:mcv,1:mcvloc) = savederv2(1:mcv,1:mcvloc)
!
!  Non-analytic correction for gamma point
!
              call nagamma(nphonatc,nphonatc,iocptr,rfmass,qfrac,maxd2,derv2)
!
!  Calculate eigenvalues and eigenvectors of corrected dynamical matrix
!
              call pdiaggd(mcv,mcvloc,maxd2,derv2,maxeigr,eigr,freq(1,k-nlkpt+1),fscale,.true.,lprint,ifail)
!
!  Find correspondence of modes by projection of eigenvectors
!

!
!  Add frequencies to average for appropriate mode
!
              do nf = 1,mcv
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
          freq(1:mcv,k-nlkpt+1) = averfreq(1:mcv)/weightpt
!
          deallocate(averfreq,stat=status)
          if (status/=0) call deallocate_error('phonon_fcd','averfreq')
        endif
!
!  Output frequencies / DOS / intensities
!
        call peigengd(mcv,mcvloc,mcvptr,mcvnptr,mcvrptr,freq(1,k-nlkpt+1),ncore,nphonatc,nphonats, &
                      nphonatptr,maxeigr,eigr,oscstrength,ramstrength,IRintensity(1,k-nlkpt+1))
        if (nobsmode.gt.0) then
!
!  Look for vibrational mode frequencies projected on to eigenvectors
!
          do nobm = 1,nobsmode
            call getvibmoded(nobsmodeptr0+nobm,maxd2,eigr,mcvloc,mcvptr,nfitmode,overlap)
            fobsmodefreq(nobsmodeptr0+nobm) = freq(nfitmode,k-nlkpt+1)
            fobsmodeover(nobsmodeptr0+nobm) = overlap
          enddo
        endif
!
!  Thermal conductivity calculation
!
        if (lthermal) then
          call thermalconductivity_af_d(mcv,mcvloc,mcvptr,mcvrptr,derv2,eigr,Sij,freq(1,k-nlkpt+1), &
                                        nphonatc,nphonatptr,maxd2,fscale)
        endif
!
!  Lower symmetry to remove imaginary modes if selected
!
        call lower(mcv,mcvrptr,freq(1,k-nlkpt+1),nphonatc,nphonatptr,nphonatc,iocptr,maxeigr,eigr)
      else
!*****************************************
!  Complex eigenvalues and eigenvectors  *
!*****************************************
        if (lintegerKpoint) then
!
!  Calculate eigenvalues and eigenvectors of dynamical matrix
!
          call pdiagd(mcv,mcvloc,maxd2,derv2,dervi,maxeigc,eigc,freq(1,k-nlkpt+1),fscale,.true.,lprint,ifail)
!
!  Store uncorrected frequencies
!
          savefreq(1:mcv) = freq(1:mcv,k-nlkpt+1)
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
                call nagammac(nphonatc,nphonatc,nphonatptr,iocptr,rfmass,bornkloc,maxd2,derv2,dervi,xkt,ykt,zkt)
              else
                call nagammac(nphonatc,nphonatc,nphonatptr,iocptr,rfmass,bornk(1,ncf),maxd2,derv2,dervi,xkt,ykt,zkt)
              endif
            endif
!
!  Calculate eigenvalues and eigenvectors of corrected dynamical matrix
!
            call pdiagd(mcv,mcvloc,maxd2,derv2,dervi,maxeigc,eigc,freq(1,k-nlkpt+1),fscale,.true.,lprint,ifail)
          elseif (nbornstep(ncf).ne.0) then
!******************************************************
!  Polycrystalline average of gamma point correction  *
!******************************************************
!
!  Allocate extra storage
!
            allocate(averfreq(mcv),stat=status)
            if (status/=0) call outofmemory('phonon_fcd','averfreq')
!
!  Store derv2
!
            savederv2(1:mcv,1:mcvloc) = derv2(1:mcv,1:mcvloc)
            savedervi(1:mcv,1:mcvloc) = dervi(1:mcv,1:mcvloc)
            averfreq(1:mcv) = 0.0_dp
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
                derv2(1:mcv,1:mcvloc) = savederv2(1:mcv,1:mcvloc)
                dervi(1:mcv,1:mcvloc) = savedervi(1:mcv,1:mcvloc)
!
!  Non-analytic correction for gamma point
!
                call nagammac(nphonatc,nphonatc,nphonatptr,iocptr,rfmass,qfrac,maxd2,derv2,dervi,xkt,ykt,zkt)
!
!  Calculate eigenvalues and eigenvectors of corrected dynamical matrix
!
                call pdiagd(mcv,mcvloc,maxd2,derv2,dervi,maxeigc,eigc,freq(1,k-nlkpt+1),fscale,.true.,lprint,ifail)
!
!  Add frequencies to average for appropriate mode
!
                do nf = 1,mcv
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
            freq(1:mcv,k-nlkpt+1) = averfreq(1:mcv)/weightpt
!
            deallocate(averfreq,stat=status)
            if (status/=0) call deallocate_error('phonon_fcd','averfreq')
          endif
!
!  Calculate the oscillator strengths
!
          call oscillatorstrengthc(mcv,mcvloc,mcvptr,nphonatc,nphonatptr,nphonatc,nphonatptr, &
                                   eigc,maxeigc,oscstrength)
!
!  Output frequencies / DOS / intensities
!
          call peigend(mcv,mcvloc,mcvptr,mcvnptr,mcvrptr,freq(1,k-nlkpt+1),ncore,nphonatc,nphonats, &
                       nphonatptr,maxeigc,eigc,oscstrength,IRintensity(1,k-nlkpt+1),wk)
!
!  Compute group velocitiies
!
          if (lgroupvelocity) then
            call groupvelocitiesd(k-nlkpt+1,mcv,mcvloc,mcvptr,eigc,maxeigc,fscale,lprint)
          endif
!
!  Lower symmetry to remove imaginary modes if selected
!
          call lower(mcv,mcvrptr,freq(1,k-nlkpt+1),nphonatc,nphonatptr,nphonatc,iocptr,maxeigr,eigr)
        else
!
!  Calculate eigenvalues and eigenvectors of dynamical matrix
!
          call pdiagd(mcv,mcvloc,maxd2,derv2,dervi,maxeigc,eigc,freq(1,k-nlkpt+1),fscale,.true.,lprint,ifail)
!
!  Calculate the oscillator strengths
!
          call oscillatorstrengthc(mcv,mcvloc,mcvptr,nphonatc,nphonatptr,nphonatc,nphonatptr, &
                                   eigc,maxeigc,oscstrength)
!
!  Output frequencies / DOS / intensities
!
          call peigend(mcv,mcvloc,mcvptr,mcvnptr,mcvrptr,freq(1,k-nlkpt+1),ncore,nphonatc,nphonats, &
                       nphonatptr,maxeigc,eigc,oscstrength,IRintensity(1,k-nlkpt+1),wk)
!
!  Compute group velocitiies
!
          if (lgroupvelocity) then
            call groupvelocitiesd(k-nlkpt+1,mcv,mcvloc,mcvptr,eigc,maxeigc,fscale,lprint)
          endif
!
!  Lower symmetry to remove imaginary modes if selected
!
          call lower(mcv,mcvrptr,freq(1,k-nlkpt+1),nphonatc,nphonatptr,nphonatc,iocptr,maxeigr,eigr)
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
              call nagammac(nphonatc,nphonatc,nphonatptr,iocptr,rfmass,bornkloc,maxd2,derv2,dervi,xkt,ykt,zkt)
            else
              call nagammac(nphonatc,nphonatc,nphonatptr,iocptr,rfmass,bornk(1,ncf),maxd2,derv2,dervi,xkt,ykt,zkt)
            endif
          else
            if (lpartofdisp) then
              call nagamma(nphonatc,nphonatc,iocptr,rfmass,bornkloc,maxd2,derv2)
            else
              call nagamma(nphonatc,nphonatc,iocptr,rfmass,bornk(1,ncf),maxd2,derv2)
            endif
          endif
        endif
!
!  Calculate eigenvalues of uncorrected dynamical matrix
!
        if (lintegerKpoint) then
          call pdiagd(mcv,mcvloc,maxd2,derv2,dervi,maxeigc,eigc,freq(1,k-nlkpt+1),fscale,.false.,lprint,ifail)
        else
          call pdiaggd(mcv,mcvloc,maxd2,derv2,maxeigr,eigr,freq(1,k-nlkpt+1),fscale,.false.,lprint,ifail)
        endif
      else
!
!  Calculate eigenvalues of dynamical matrix
!
        call pdiagd(mcv,mcvloc,maxd2,derv2,dervi,maxeigc,eigc,freq(1,k-nlkpt+1),fscale,.false.,lprint,ifail)
      endif
!
!  Output frequencies
!
      if (lfreqout.and.lprinloc) then
        write(ioout,'(/,''  Frequencies (cm-1) [NB: Negative implies an imaginary mode]:'',/)')
        write(ioout,'(9f8.2)')(freq(i,k-nlkpt+1),i=1,mcv)
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
          do i = 1,mcv
            write(51) freq(i,k-nlkpt+1)
          enddo
        else
          do i = 1,mcv
            write(52,"(f" // adjustl(fstring) // ")") freq(i,k-nllkpt+1)
          enddo
        endif
      endif
      t2 = g_cpu_time()
      tdisk = tdisk + t2 - t1
    elseif (lfrq.and.ioproc) then
      t1 = g_cpu_time()
      do i = 1,mcv
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
!  The following is equivalent to looping over mcvloc
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
              do j = 1,mcv
                jj = (j-1)/3_i4 + 1
                projectionfactor = eigr(j,iloc)**2
                meanKEperatom(jj) = meanKEperatom(jj) + trmke*projectionfactor
              enddo
            else
              do j = 1,mcv
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
        allocate(sum(mcv),stat=status)
        if (status/=0) call outofmemory('phonon_fcd','sum')
!
        call sumall(meanKEperatom,sum,mcv,"phonon_fcd","meanKEperatom")
        do i = 1,mcv
          meanKEperatom(i) = sum(i)
        enddo
!
        deallocate(sum,stat=status)
        if (status/=0) call deallocate_error('phonon_fcd','sum')
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
      call omegaproperty(mcv,savefreq,fscale,oscstrength)
    endif
!***************************
!  End loop over K points  *
!***************************
  enddo
  if (lpdf) then
!
!  Set up variables for gulp_cml_PDFstats
!
    call pdfsetup(lcml) 
  endif
!**************************************************
!  Deallocate second derivative workspace memory  *
!**************************************************
  if (lgammaonly) then
    deallocate(eigr,stat=status)
    if (status/=0) call deallocate_error('phonon_fcd','eigr')
  else
    deallocate(eigc,stat=status)
    if (status/=0) call deallocate_error('phonon_fcd','eigc')
    if (allocated(eigr)) then
      deallocate(eigr,stat=status)
      if (status/=0) call deallocate_error('phonon_fcd','eigr')
    endif
  endif
  if (allocated(savedervi)) then
    deallocate(savedervi,stat=status)
    if (status/=0) call deallocate_error('phonon_fcd','savedervi')
  endif
  if (allocated(savederv2)) then
    deallocate(savederv2,stat=status)
    if (status/=0) call deallocate_error('phonon_fcd','savederv2')
  endif
!*************************************
!  Output phonon related properties  *
!*************************************
!
!  Compute vibrational energies and output
!
  call vibenergy(mcv,nllkpt,wkpt(nlkpt),freq(1,1),maxfreq,rnokpt,lprint,fc)
!
!  Mean squared displacements
!
  if (lmsd) then
!
!  Global reduction of array
!
    allocate(sum(3*ncfoc),stat=status)
    if (status/=0) call outofmemory('phonon_fcd','sum')
    allocate(sum2(3*ncfoc),stat=status)
    if (status/=0) call outofmemory('phonon_fcd','sum2')
!
    sum(1:ncfoc) = msdx(1:ncfoc)
    sum(ncfoc+1:2*ncfoc) = msdy(1:ncfoc)
    sum(2*ncfoc+1:3*ncfoc) = msdz(1:ncfoc)
!
    call sumall(sum,sum2,3_i4*ncfoc,"phonon_fcd","msd")
!
    msdx(1:ncfoc) = sum2(1:ncfoc)
    msdy(1:ncfoc) = sum2(ncfoc+1:2*ncfoc)
    msdz(1:ncfoc) = sum2(2*ncfoc+1:3*ncfoc)
!
    deallocate(sum2,stat=status)
    if (status/=0) call deallocate_error('phonon_fcd','sum2')
    deallocate(sum,stat=status)
    if (status/=0) call deallocate_error('phonon_fcd','sum')
!
    if (lprinloc) then
      write(ioout,'(/,''  Mean-squared displacements for phonons (Ang**2) : '',/)')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''   Atom '',8x,''x'',12x,''y'',12x,''z'')')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      do i = 1,ncfoc
        write(ioout,'(i8,1x,3(f12.6,1x))') i,msdx(i),msdy(i),msdz(i)
      enddo
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
    call outphon(nlkpt,nukpt,leigloc,3_i4*nphonatc)
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
  if (lcml.and.lpdf) then
    call gulp_cml_PDFstats
    call gulp_cml_outPDF
  endif
!
!  Free local memory
!
  if (lthermal) then
    deallocate(Sij,stat=status)
    if (status/=0) call deallocate_error('phonon_fcd','Sij')
  endif
  if (leigloc) then
    deallocate(ramstrength,stat=status)
    if (status/=0) call deallocate_error('phonon_fcd','ramstrength')
    deallocate(oscstrength,stat=status)
    if (status/=0) call deallocate_error('phonon_fcd','oscstrength')
    deallocate(savefreq,stat=status)
    if (status/=0) call deallocate_error('phonon_fcd','savefreq')
  endif
  if (allocated(kpvt)) then
    deallocate(kpvt,stat=status)
    if (status/=0) call deallocate_error('phonon_fcd','kpvt')
  endif
  deallocate(mcvrptr,stat=status)
  if (status/=0) call deallocate_error('phonon_fcd','mcvrptr')
  deallocate(mcvnptr,stat=status)
  if (status/=0) call deallocate_error('phonon_fcd','mcvnptr')
  deallocate(mcvptr,stat=status)
  if (status/=0) call deallocate_error('phonon_fcd','mcvptr')
  deallocate(meanKEperatom,stat=status)
  if (status/=0) call deallocate_error('phonon_fcd','meanKEperatom')
!
!  Close PDF related parts of phonon 
!
  call closepdfphonon
!
!  Timing
!
  t2t = g_cpu_time()
  tphon = t2t - t1t + tphon
#ifdef TRACE
  call trace_out('phonon_fcd')
#endif
#else
  call outerror('phonon_fcd called when not compiled with MPI',0_i4)
  call stopnow('phonon_fcd')
#endif
!
  return
  end
