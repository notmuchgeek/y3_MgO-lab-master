  subroutine phonon_fcg(lprint,fc,nobsmodeptr0,nobsmode)
!
!  Calculates the phonons at a given set of k points.
!  This version uses a different algorithm in which the second
!  derivatives are stored by cell image and only calculated once.
!  The phasing for k points is then applied later.
!  NB: In this version the force constants are computed for a
!  supercell and then contracted to the original ghostcell
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
!   1/18 Created from phonon_fc
!   1/18 Trace added
!   3/18 dynam shortened to dyna in keyword check
!   4/18 leigloc now set here for llower and linten
!   3/19 Multiple temperature ramps added
!   5/19 Finite difference flag split for first and second derivatives
!   8/19 Calls to peigen corrected for missing argument
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!   4/20 maxfqat changed to maxfreq since it is no longer just atom based
!   4/20 Mass arrays now use fmass and rfmass
!   5/20 ncfoc removed from arguments to groupvelocities
!   6/20 Last argument to outphon changed
!   7/20 bornk passed as local array
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
  use gulp_cml,        only : lcml
  use gulp_cml_phonon, only : gulp_cml_startPhonons
  use gulp_cml_phonon, only : gulp_cml_addKPoint, gulp_cml_endPhonons
  use iochannels
  use ksample
  use ksample_scatter
  use maths,           only : leispack_eigensolve
  use m_pdfneutron,    only : lpdf
  use observables,     only : fobsmodefreq, fobsmodeover
  use parallel
  use partial
  use phononatoms
  use projectdos
  use properties
  use scatterdata,     only : lscattercall
  use species,         only : massspec, natspec, ntypspec
  use shells
  use sutton,          only : lsuttonc
  use times
#ifdef TRACE
  use trace,           only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4),  intent(in)                          :: nobsmodeptr0  ! Pointer to first observable mode - 1
  integer(i4),  intent(in)                          :: nobsmode      ! Number of observable modes if fitting run
  logical,      intent(in)                          :: lprint        ! If true then output results
  real(dp),     intent(in)                          :: fc            ! Internal energy
!
!  Local variables
!
  character(len=5)                                  :: lab1
  character(len=2)                                  :: fstring1
  character(len=2)                                  :: fstring2
  character(len=12)                                 :: fstring
  complex(dpc), dimension(:),     allocatable       :: ctmp
  integer(i4)                                       :: i
  integer(i4)                                       :: ifail
  integer(i4)                                       :: ig
  integer(i4)                                       :: ii
  integer(i4)                                       :: indgi
  integer(i4)                                       :: indi
  integer(i4)                                       :: indii
  integer(i4)                                       :: indj
  integer(i4)                                       :: indjj
  integer(i4)                                       :: inert(3)
  integer(i4),  dimension(:),     allocatable       :: ipivot
  integer(i4)                                       :: j
  integer(i4)                                       :: jg
  integer(i4)                                       :: jj
  integer(i4)                                       :: job
  integer(i4)                                       :: k
  integer(i4)                                       :: kk
  integer(i4),  dimension(:),     allocatable       :: kpvt
  integer(i4)                                       :: l
  integer(i4)                                       :: m
  integer(i4)                                       :: maxlim
  integer(i4)                                       :: mcv
  integer(i4),  dimension(:),     allocatable       :: mcvptr
  integer(i4),  dimension(:),     allocatable       :: mcvrptr
  integer(i4)                                       :: mint
  integer(i4)                                       :: msv
  integer(i4)                                       :: nbfocg
  integer(i4)                                       :: ncfocg
  integer(i4)                                       :: nsfocg
  integer(i4)                                       :: ncind
  integer(i4)                                       :: nd
  integer(i4)                                       :: nf
  integer(i4)                                       :: nfitmode
  integer(i4)                                       :: nghostatom
  integer(i4)                                       :: nghostcell
  integer(i4)                                       :: nghostcore
  integer(i4)                                       :: nk
  integer(i4)                                       :: nldpt
  integer(i4)                                       :: nlkpt
  integer(i4)                                       :: nllkpt
  integer(i4)                                       :: nobm
  integer(i4)                                       :: np
  integer(i4)                                       :: np_nlkpt
  integer(i4)                                       :: np_procs
  integer(i4)                                       :: nphonatg
  integer(i4)                                       :: nsi
  integer(i4)                                       :: nt
  integer(i4)                                       :: ntmax
  integer(i4)                                       :: nudpt
  integer(i4)                                       :: nukpt
  integer(i4)                                       :: status
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
  real(dp),     dimension(:),     allocatable       :: averfreq
  real(dp)                                          :: bornkloc(3)
  real(dp)                                          :: bornkpt(3)
  real(dp)                                          :: cmfact
  real(dp)                                          :: g_cpu_time
  real(dp)                                          :: det(2)
  real(dp)                                          :: eigcomp
  real(dp),     dimension(:),     allocatable       :: eigr
  real(dp)                                          :: eigreal
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
  real(dp),     dimension(:),     allocatable       :: wtmp
  real(dp)                                          :: weightpt
  real(dp)                                          :: wi
  real(dp)                                          :: wk
  real(dp)                                          :: wr
  real(dp)                                          :: xkt
  real(dp)                                          :: ykt
  real(dp)                                          :: zkt
  real(dp)                                          :: xmod
  real(dp)                                          :: ymod
  real(dp)                                          :: zmod
#ifdef TRACE
  call trace_in('phonon_fcg')
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
  if (lthermal) leigloc = .true.
  if (lgroupvelocity) leigloc = .true.
  if (linten) leigloc = .true.
  if (llower) leigloc = .true.
  lprinloc = (lprint.and.ioproc)
  ltemprop = (temperature.gt.1.0d-6.or.ntemperatureramp.gt.0)
  lnonanal = (index(keyword,'nono').eq.0.and.ndim.eq.3.and..not.leem.and..not.lnoanald2)
  lnoanald2loc = (lnoanald2)
!
!  Trap features that are not compatible with this algorithm
!
  if (lpdf) then
    call outerror('PDF calculation not yet implemented for ghost cell phonons',0_i4)
    call stopnow('phonon_fcg')
  endif
  if (nbsmat.gt.0) then
    call outerror('Breathing shells not yet implemented for ghost cell phonons',0_i4)
    call stopnow('phonon_fcg')
  endif
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
    call trace_out('phonon_fcg')
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
!  Find number of ghost cells
!
  nghostcell = nsuperghost(1,ncf)*nsuperghost(2,ncf)*nsuperghost(3,ncf)
!
!  Check the number of atoms is divisable by the number of ghost cells
!
  if (mod(numat,nghostcell).ne.0) then
    call outerror('number of atoms is not compatible with ghost_supercell',0_i4)
    call stopnow('phonon_fcg')
  endif
!
  nghostatom = numat/nghostcell
  nghostcore = ncore/nghostcell
!
  nphonatg = nphonat/nghostcell
  nbfocg = nbfoc/nghostcell
  ncfocg = ncfoc/nghostcell
  nsfocg = nsfoc/nghostcell
!
!  Allocate local pointer arrays
!
  allocate(meanKEperatom(nghostcore),stat=status)
  if (status/=0) call outofmemory('phonon_fcg','meanKEperatom')
!
  lpocc = (nsfoc+ncfoc.ne.nphonat)
!
!  At present partial occupancy is not allowed with this algorithm
!
  if (lpocc) then
    call outerror('partial occupancy not yet implemented for ghost cell phonons',0_i4)
    call stopnow('phonon_fcg')
  endif
!
!  Calculate a few constants to do with the size of the problem
!
  mint = 3*nphonatg
  maxlim = mint
  if (nbsmat.gt.0) maxlim = maxlim + nphonatg
  msv = 3*nsfocg + nbfocg
  mcv = 3*ncfocg
!
  allocate(mcvptr(mcv),stat=status)
  if (status/=0) call outofmemory('phonon_fcg','mcvptr')
  allocate(mcvrptr(mcv),stat=status)
  if (status/=0) call outofmemory('phonon_fcg','mcvrptr')
!
  do m = 1,mcv
    mcvptr(m) = m 
    mcvrptr(m) = m 
  enddo
!
!  Check that maxd2 is greater than or equal to mcv
!
  if (maxd2.lt.mcv) then
    maxd2 = mcv
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
  call changemaxfreq(3_i4*ncfocg)
!
!  Allocate array to hold oscillator strengths if needed
!
  if (leigloc) then
    allocate(savefreq(mcv),stat=status)
    if (status/=0) call outofmemory('phonon_fcg','savefreq')
    allocate(oscstrength(3,3,mcv),stat=status)
    if (status/=0) call outofmemory('phonon_fcg','oscstrength')
    allocate(ramstrength(3,3,mcv),stat=status)
    if (status/=0) call outofmemory('phonon_fcg','ramstrength')
  endif
!
!  Allocate array for thermal conductivity if needed
!
  if (lthermal) then
    allocate(Sij(mcv,mcv),stat=status)
    if (status/=0) call outofmemory('phonon_fcg','Sij')
  endif
!
  rkt = boltz*temperature
!
!  If mean KE per atom is requested then initialise array to zero
!
  if (lmeanke) then
    meanKEperatom(1:ncfocg) = 0.0_dp
  endif
!
!  Calculate inversion square root of masses for atoms and momemt of inertia for quaternions
!
!  Now modified to handle partial occupancies
!
  do i = 1,3*ncfocg
    fmass(i) = 0.0_dp
  enddo
  do i = 1,nphonatc,nghostcell
    ii = iocptr(i)
    ig = (ii - 1)/nghostcell + 1
    nsi = nspecptr(nrelf2a(nphonatptr(i)))
    rmassi = massspec(nsi)*occuf(nphonatptr(i))
    if (natspec(nsi).le.maxele.and.abs(rmassi).lt.1.0d-12) then
      call label(natspec(nsi),ntypspec(nsi),lab1)
      call outerror('mass of species '//lab1//' is zero',0_i4)
      goto 999
    endif
    fmass(3*(ig-1)+1) = fmass(3*(ig-1)+1) + rmassi
    fmass(3*(ig-1)+2) = fmass(3*(ig-1)+2) + rmassi
    fmass(3*(ig-1)+3) = fmass(3*(ig-1)+3) + rmassi
  enddo
  do i = 1,3*ncfocg
    if (fmass(i).eq.0.0_dp) then
      call outerror('site has total mass of zero in phonon_fcg',0_i4)
      call stopnow('phonon_fcg')
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
    call stopnow('phonon_fcg')
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
    write(ioout,'(/,''  Phonon Calculation using ghost cell: '',/)')
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
    write(53,'(i6)') nphonatc/nghostcell
    do i = 1,nphonatc,nghostcell
      ii = nphonatptr(i)
      write(53,'(i3,1x,3(f15.6,1x))') nat(ii),xclat(ii),yclat(ii),zclat(ii)
    enddo
    write(53,'(i6)') nllkpt
    write(53,'(i6)') 3*nphonatc/nghostcell
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
    write(54,'(a22,i6)') " Number of ions       ", nphonatc/nghostcell
    write(54,'(a22,i6)') " Number of branches   ", 3*ncfocg
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
      write(54,'(3f12.6)') (rv(j,i)/dble(nsuperghost(i,ncf)),j=1,3)
    enddo
!
!  Write out a list of fully occupied atoms in the FULL unit cell
!
    write(54,'(a)') " Fractional Co-ordinates"
    ig = 0
    do i = 1,nphonatc,nghostcell
      ig = ig + 1
      ii = nphonatptr(i)
      call label(nat(ii),nftype(ii),lab1)
      write(54,'(i6,1x,3f12.6,3x,a5,1x,f12.6)') ig,xfrac(ii)*dble(nsuperghost(1,ncf)), &
                                                   yfrac(ii)*dble(nsuperghost(2,ncf)), &
                                                   zfrac(ii)*dble(nsuperghost(3,ncf)), &
                                                   lab1,fmass(3*(ii-1)+1)
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
    if (lgammaonly.and.nbornstep(ncf).eq.0) then
      allocate(eigr(maxd2*maxd2),stat=status)
    else
      allocate(eigr(2*maxd2*maxd2),stat=status)
    endif
    if (status/=0) call outofmemory('phonon_fcg','eigr')
    if (nbornstep(ncf).gt.0.and..not.lnonanal) then
      allocate(savederv2(mcv,mcv),stat=status)
      if (status/=0) call outofmemory('phonon_fcg','savederv2')
      if (lintegerKpoint_present) then
        allocate(savedervi(mcv,mcv),stat=status)
        if (status/=0) call outofmemory('phonon_fcg','savedervi')
      endif
    endif
  else
    if (leispack_eigensolve) then
      if (msv.gt.0) then
        if (lgammaonly) then
          allocate(eigr(msv*msv),stat=status)
        else
          allocate(eigr(2*msv*msv),stat=status)
        endif
      else
        allocate(eigr(1),stat=status)
      endif
    elseif (.not.lgammaonly) then
      allocate(eigr(2*maxd2*mcv),stat=status)
    elseif (msv.gt.0) then
      allocate(eigr(msv*(msv+1)/2),stat=status)
    else
      allocate(eigr(1),stat=status)
    endif
    if (status/=0) call outofmemory('phonon_fcg','eigr')
  endif
!*************************************************************
!  Compute second derivatives once and store per cell image  *
!*************************************************************
  if (lused1) then
    call dynamicn_fc
  else
    call dynamic_fc
  endif
!
!  Set diagonal elements
!
  call diagonal_fc(lused1)
!
!  Store diagonal blocks in derv3 to avoid recalculation
!
  indgi = 0
  do i = 1,nphonat,nghostcell
    indi = 3*(nphonatptr(i)-1)
    derv3(indgi+1,1) = derv2(indi+1,indi+1)
    derv3(indgi+2,1) = derv2(indi+2,indi+1)
    derv3(indgi+3,1) = derv2(indi+3,indi+1)
    derv3(indgi+1,2) = derv2(indi+1,indi+2)
    derv3(indgi+2,2) = derv2(indi+2,indi+2)
    derv3(indgi+3,2) = derv2(indi+3,indi+2)
    derv3(indgi+1,3) = derv2(indi+1,indi+3)
    derv3(indgi+2,3) = derv2(indi+2,indi+3)
    derv3(indgi+3,3) = derv2(indi+3,indi+3)
    indgi = indgi + 3
  enddo
  if (lshengBTE) then
!
!  Output shengBTE second derivative file using supercell format
!
    call outshengBTEfc2s(iout_shengBTE)
  endif
!***********************
!  Loop over k points  *
!***********************
  np_nlkpt = nlkpt
  np_procs = 1_i4
  do k = np_nlkpt,nukpt,np_procs
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
    call phase_fcg(xkt,ykt,zkt)
!
!  Include diagonal blocks, stored in derv3
!
    indi = 0
    do i = 1,nphonatg
      derv2(indi+1,indi+1) = derv2(indi+1,indi+1) + derv3(indi+1,1)
      derv2(indi+2,indi+1) = derv2(indi+2,indi+1) + derv3(indi+2,1)
      derv2(indi+3,indi+1) = derv2(indi+3,indi+1) + derv3(indi+3,1)
      derv2(indi+1,indi+2) = derv2(indi+1,indi+2) + derv3(indi+1,2)
      derv2(indi+2,indi+2) = derv2(indi+2,indi+2) + derv3(indi+2,2)
      derv2(indi+3,indi+2) = derv2(indi+3,indi+2) + derv3(indi+3,2)
      derv2(indi+1,indi+3) = derv2(indi+1,indi+3) + derv3(indi+1,3)
      derv2(indi+2,indi+3) = derv2(indi+2,indi+3) + derv3(indi+2,3)
      derv2(indi+3,indi+3) = derv2(indi+3,indi+3) + derv3(indi+3,3)
      indi = indi + 3
    enddo
!*****************************************************************
!  Compress second derivatives according to partial occupancies  *
!*****************************************************************
! DEBUG - needs handling
    if (lpocc) then
      ncsfoc = ncfoc + nsfoc
      call compressd2(derv2,maxd2,ncfoc,nsfoc,nbfoc,numat,iocptr,ibocptr)
      if (.not.lgamma) then
        call compressd2(dervi,maxd2,ncfoc,nsfoc,nbfoc,numat,iocptr,ibocptr)
      endif
    elseif (numat.ne.nphonat) then
!*********************************************************************
!  Compress full second derivatives down to region 1 only if needed  *
!*********************************************************************
      ig = 0
      do i = 1,nphonat,nghostcell
        ig = ig + 1
        indi  = 3*(ig-1)
        indii = 3*(nphonatptr(i)-1)/nghostcell
        jg = 0
        do j = 1,nphonat,nghostcell
          jg = jg + 1
          indj  = 3*(jg-1)
          indjj = 3*(nphonatptr(j)-1)/nghostcell
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
      enddo
    endif
!
!  Output the uncompressed second derivatives for debugging
!
    if (index(keyword,'dyna').ne.0.and.index(keyword,'debu').ne.0.and.ioproc) then
      write(ioout,'(/,''  Uncompressed Real Dynamical matrix :'',/)')
      do i = 1,mcv+msv
        write(ioout,'(12f11.6)')(derv2(j,i),j=1,mcv+msv)
      enddo
      write(ioout,'(/,''  Uncompressed Imag Dynamical matrix :'',/)')
      do i = 1,mcv+msv
        write(ioout,'(12f11.6)')(dervi(j,i),j=1,mcv+msv)
      enddo
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
!  Allocate workspace for inversion
!     
        allocate(ipivot(msv),stat=status)
        if (status/=0) call outofmemory('phonon_fcg','ipivot')                  
        allocate(wtmp(3*msv),stat=status)
        if (status/=0) call outofmemory('phonon_fcg','wtmp')
!
!  Transfer data to packed storage
!    
        kk = 0
        do i = 1,msv
          do j = 1,i
            kk = kk + 1
            eigr(kk) = derv2(mcv+j,mcv+i)
          enddo
        enddo                                    
!         
!  Factorise matrix
!  
        call dsptrf('U',msv,eigr,ipivot,ifail)
        if (ifail.eq.0) then
!
!  Form inverse
!
          call dsptri('U',msv,eigr,ipivot,wtmp,ifail)
!
!  Transfer data back
!
          kk = 0
          do i = 1,msv
            do j = 1,i
              kk = kk + 1
              derv2(mcv+j,mcv+i) = eigr(kk)
              derv2(mcv+i,mcv+j) = eigr(kk)
            enddo
          enddo
        endif
!
!  Free workspace  
!
        deallocate(wtmp,stat=status)
        if (status/=0) call deallocate_error('phonon_fcg','wtmp')
        deallocate(ipivot,stat=status)
        if (status/=0) call deallocate_error('phonon_fcg','ipivot')  
!
        t2i = g_cpu_time()
        tmati = tmati + t2i - t1i
        if (ifail.ne.0) then
          call outerror('inversion of shell 2nd derivatives failed',0_i4)
          call stopnow('phonon_fcg')
        endif
!***********************************************
!  Corrected second derivatives = R - T*S-1*T  *
!***********************************************
!
!  First pass : S-1*T stored in diagonally opposite copy of T
!
        do i = 1,mcv
          do j = 1,msv
            wr = 0.0_dp
            do l = 1,msv
              wr = wr + derv2(j+mcv,l+mcv)*derv2(i,l+mcv)
            enddo
            derv2(mcv+j,i) = wr
          enddo
        enddo
!
!  Second pass : T*(S-1*T)
!
        do i = 1,mcv
          do j = 1,mcv
            wr = 0.0_dp
            do l = 1,msv
              wr = wr - derv2(j,l+mcv)*derv2(mcv+l,i)
            enddo
            derv2(j,i) = derv2(j,i) + wr
          enddo
        enddo
      else
!*****************************
!  Complex Matrix Inversion  *
!*****************************
!
!  Transfer real and imaginary components to complex matrix 
!
        call phoncopy1(derv2,dervi,eigr,maxd2,mcv,mcv,msv,msv)
!
!  Invert complex matrix
!
        job = 1
        t1i = g_cpu_time()
        allocate(kpvt(msv),stat=status)
        if (status/=0) call outofmemory('phonon_fcg','kpvt')
        call zhifa(eigr,msv,msv,kpvt,ifail)
        if (ifail.ne.0) then
          call outerror('inversion of shell 2nd derivatives failed',0_i4)
          goto 999
        endif
        allocate(ctmp(msv),stat=status)
        if (status/=0) call outofmemory('phonon_fcg','ctmp')
        call zhidi(eigr,msv,msv,kpvt,det,inert,ctmp,job)
        deallocate(ctmp,stat=status)
        if (status/=0) call deallocate_error('phonon_fcg','ctmp')
        deallocate(kpvt,stat=status)
        if (status/=0) call deallocate_error('phonon_fcg','kpvt')
        t2i = g_cpu_time()
        tmati = tmati + t2i - t1i
!
!  Return inverse complex matrix to separate real and imaginary matrices and resymmetrise
!
        call phoncopy2e(0_i4,derv2,dervi,eigr,maxd2,mcv,mcv,msv,msv)
!***********************************************
!  Corrected second derivatives = R - T*S-1*T  *
!***********************************************
!
!  First pass : S-1*T stored in diagonally opposite copy of T
!
        do i = 1,mcv
          do j = 1,msv
            wr = 0.0_dp
            wi = 0.0_dp
            do l = 1,msv
              wr = wr + derv2(j+mcv,l+mcv)*derv2(i,l+mcv)
              wr = wr + dervi(j+mcv,l+mcv)*dervi(i,l+mcv)
              wi = wi + dervi(j+mcv,l+mcv)*derv2(i,l+mcv)
              wi = wi - derv2(j+mcv,l+mcv)*dervi(i,l+mcv)
            enddo
            derv2(mcv+j,i) = wr
            dervi(mcv+j,i) = wi
          enddo
        enddo
!
!  Second pass : T*(S-1*T) - for imaginary case sign has to be adjusted
!  to allow for the fact that T(conj) and T are of opposite signs
!
        do i = 1,mcv
          do j = 1,mcv
            wr = 0.0_dp
            wi = 0.0_dp
            do l = 1,msv
              wr = wr - derv2(j,l+mcv)*derv2(mcv+l,i)
              wr = wr + dervi(j,l+mcv)*dervi(mcv+l,i)
              wi = wi - derv2(j,l+mcv)*dervi(mcv+l,i)
              wi = wi - dervi(j,l+mcv)*derv2(mcv+l,i)
            enddo
            derv2(j,i) = derv2(j,i) + wr
            dervi(j,i) = dervi(j,i) + wi
          enddo
        enddo
      endif
    endif
!****************************
!  End of shell correction  *
!****************************
!
!  Option to write out a .frc file for QMPOT
!
    if (lfrc.and.ioproc.and.lgamma) call outfrc(fc,.true.,.true.)
!*****************************
!  Multiply by mass-factors  *
!*****************************
    if (lgamma) then
      do i = 1,3*ncfocg
        rmassi = rfmass(i)
        do j = 1,3*ncfocg
          derv2(j,i) = rmassi*rfmass(j)*derv2(j,i)
        enddo
      enddo
    else
      do i = 1,3*ncfocg
        rmassi = rfmass(i)
        do j = 1,3*ncfocg
          derv2(j,i) = rmassi*rfmass(j)*derv2(j,i)
          dervi(j,i) = rmassi*rfmass(j)*dervi(j,i)
        enddo
      enddo
    endif
!
!  If debugging print out dynamical matrix
!
    if (index(keyword,'dyna').ne.0.and.ioproc) then
      write(ioout,'(/,''  Real Dynamical matrix :'',/)')
      do i = 1,mcv
        write(ioout,'(12f11.6)')(derv2(j,i),j=1,mcv)
      enddo
      if (.not.lgamma) then
        write(ioout,'(/,''  Imaginary Dynamical matrix :'',/)')
        do i = 1,mcv
          write(ioout,'(12f11.6)')(dervi(j,i),j=1,mcv)
        enddo
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
        call pdiagg(mcv,maxd2,derv2,eigr,freq(1,k-nlkpt+1),fscale,.true.,lprint,ifail)
!
!  Calculate the oscillator strengths
!
        call oscillatorstrengthg(mcv,mcv,mcvptr,nphonatc,nphonatptr,ncfocg,iocptr, &
                                 eigr,maxd2,oscstrength)
!
!  Optionally call Raman strengths
!
        if (lraman) then
          call ramanstrength(mcv,mcv,mcvptr,nphonatc,nphonatptr,ncfocg,iocptr, &
                             eigr,maxd2,ramstrength)
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
              call nagamma(ncfocg,nphonatc,iocptr,rfmass,bornkloc,maxd2,derv2)
            else
              bornkpt(1:3) = bornk(1:3,ncf)
              call nagamma(ncfocg,nphonatc,iocptr,rfmass,bornkpt,maxd2,derv2)
            endif
          endif
!
!  Calculate eigenvalues and eigenvectors of corrected dynamical matrix
!
          call pdiagg(mcv,maxd2,derv2,eigr,freq(1,k-nlkpt+1),fscale,.false.,lprint,ifail)
        elseif (nbornstep(ncf).ne.0) then
!******************************************************
!  Polycrystalline average of gamma point correction  *
!******************************************************
!
!  Allocate extra storage
!
          allocate(averfreq(mcv),stat=status)
          if (status/=0) call outofmemory('phonon_fcg','averfreq')
!
!  Store derv2
!
          savederv2(1:mcv,1:mcv) = derv2(1:mcv,1:mcv)
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
              derv2(1:mcv,1:mcv) = savederv2(1:mcv,1:mcv)
!
!  Non-analytic correction for gamma point
!
              call nagamma(ncfocg,nphonatc,iocptr,rfmass,qfrac,maxd2,derv2)
!
!  Calculate eigenvalues and eigenvectors of corrected dynamical matrix
!
              call pdiagg(mcv,maxd2,derv2,eigr(maxd2*maxd2+1),freq(1,k-nlkpt+1),fscale,.true.,lprint,ifail)
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
          if (status/=0) call deallocate_error('phonon_fcg','averfreq')
        endif
!
!  Output frequencies / DOS / intensities
!
        call peigeng(mcv,freq(1,k-nlkpt+1),ncore,ncfocg,iocptr,maxd2,eigr, &
                     oscstrength,ramstrength,IRintensity(1,k-nlkpt+1))
        if (nobsmode.gt.0) then
!
!  Look for vibrational mode frequencies projected on to eigenvectors
!
          do nobm = 1,nobsmode
            call getvibmode(nobsmodeptr0+nobm,maxd2,eigr,nfitmode,overlap)
            fobsmodefreq(nobsmodeptr0+nobm) = freq(nfitmode,k-nlkpt+1)
            fobsmodeover(nobsmodeptr0+nobm) = overlap
          enddo
        endif
!
!  Thermal conductivity calculation
!
        if (lthermal) then
          call thermalconductivity_af(mcv,derv2,eigr,Sij,freq(1,k-nlkpt+1),nphonatc,ncfocg,nphonatptr,maxd2,fscale)
        endif
!
!  Lower symmetry to remove imaginary modes if selected
!
        call lower(mcv,mcvrptr,freq(1,k-nlkpt+1),nphonatc,nphonatptr,ncfocg,iocptr,maxd2,eigr)
      else
!*****************************************
!  Complex eigenvalues and eigenvectors  *
!*****************************************
        if (lintegerKpoint) then
!
!  Calculate eigenvalues and eigenvectors of dynamical matrix
!
          if (leispack_eigensolve) then
            call pdiage(mcv,maxd2,derv2,dervi,eigr(1),eigr(maxd2*maxd2+1),freq(1,k-nlkpt+1),fscale,.true.,lprint,ifail)
          else
            call pdiagl(k-nlkpt+1,mcv,maxd2,derv2,dervi,eigr(1),eigr(maxd2*maxd2+1),freq(1,k-nlkpt+1),fscale,.true.,lprint,ifail)
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
                call nagammac(ncfocg,nphonatc,nphonatptr,iocptr,rfmass,bornkloc,maxd2,derv2,dervi,xkt,ykt,zkt)
              else
                bornkpt(1:3) = bornk(1:3,ncf)
                call nagammac(ncfocg,nphonatc,nphonatptr,iocptr,rfmass,bornkpt,maxd2,derv2,dervi,xkt,ykt,zkt)
              endif
            endif
!
!  Calculate eigenvalues and eigenvectors of corrected dynamical matrix
!
            if (leispack_eigensolve) then
              call pdiage(mcv,maxd2,derv2,dervi,eigr(1),eigr(maxd2*maxd2+1),freq(1,k-nlkpt+1),fscale,.true.,lprint,ifail)
            else
              call pdiagl(k-nlkpt+1,mcv,maxd2,derv2,dervi,eigr(1),eigr(maxd2*maxd2+1),freq(1,k-nlkpt+1),fscale,.true.,lprint,ifail)
            endif
          elseif (nbornstep(ncf).ne.0) then
!******************************************************
!  Polycrystalline average of gamma point correction  *
!******************************************************
!
!  Allocate extra storage
!
            allocate(averfreq(mcv),stat=status)
            if (status/=0) call outofmemory('phonon_fcg','averfreq')
!
!  Store derv2
!
            savederv2(1:mcv,1:mcv) = derv2(1:mcv,1:mcv)
            savedervi(1:mcv,1:mcv) = dervi(1:mcv,1:mcv)
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
                derv2(1:mcv,1:mcv) = savederv2(1:mcv,1:mcv)
                dervi(1:mcv,1:mcv) = savedervi(1:mcv,1:mcv)
!
!  Non-analytic correction for gamma point
!
                call nagammac(ncfocg,nphonatc,nphonatptr,iocptr,rfmass,qfrac,maxd2,derv2,dervi,xkt,ykt,zkt)
!
!  Calculate eigenvalues and eigenvectors of corrected dynamical matrix
!
                if (leispack_eigensolve) then
                  call pdiage(mcv,maxd2,derv2,dervi,eigr(1),eigr(maxd2*maxd2+1),freq(1,k-nlkpt+1),fscale,.true.,lprint,ifail)
                else
                  call pdiagl(k-nlkpt+1,mcv,maxd2,derv2,dervi,eigr(1),eigr(maxd2*maxd2+1),freq(1,k-nlkpt+1),fscale,.true., &
                              lprint,ifail)
                endif
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
            if (status/=0) call deallocate_error('phonon_fcg','averfreq')
          endif
!
!  Calculate the oscillator strengths
!
          call oscillatorstrength(mcv,nphonatc,nphonatptr,ncfocg,iocptr,eigr(1),eigr(maxd2*maxd2+1), &
                                  maxd2,oscstrength)
!
!  Output frequencies / DOS / intensities
!
          call peigen(mcv,freq(1,k-nlkpt+1),ncore,nphonatc,ncfocg,nsfocg,iocptr,maxd2, &
                      eigr(1),eigr(maxd2*maxd2+1),oscstrength,IRintensity(1,k-nlkpt+1),wk) 
!
!  Compute group velocitiies
!
          if (lgroupvelocity) then
            call groupvelocities(k-nlkpt+1,mcv,eigr(1),eigr(maxd2*maxd2+1),maxd2,fscale,lprint)
          endif
!
!  Lower symmetry to remove imaginary modes if selected
!
          call lower(mcv,mcvrptr,freq(1,k-nlkpt+1),nphonatc,nphonatptr,ncfocg,iocptr,maxd2,eigr)
        else
!
!  Calculate eigenvalues and eigenvectors of dynamical matrix
!
          if (leispack_eigensolve) then
            call pdiage(mcv,maxd2,derv2,dervi,eigr(1),eigr(maxd2*maxd2+1),freq(1,k-nlkpt+1),fscale,.true.,lprint,ifail)
          else
            call pdiagl(k-nlkpt+1,mcv,maxd2,derv2,dervi,eigr(1),eigr(maxd2*maxd2+1),freq(1,k-nlkpt+1),fscale,.true.,lprint,ifail)
          endif
!
!  Calculate the oscillator strengths
!
          call oscillatorstrength(mcv,nphonatc,nphonatptr,ncfocg,iocptr,eigr(1),eigr(maxd2*maxd2+1), &
                                  maxd2,oscstrength)
!
!  Output frequencies / DOS / intensities
!
          call peigen(mcv,freq(1,k-nlkpt+1),ncore,nphonatc,ncfocg,nsfocg,iocptr,maxd2, &
                      eigr(1),eigr(maxd2*maxd2+1),oscstrength,IRintensity(1,k-nlkpt+1),wk) 
!
!  Compute group velocitiies
!
          if (lgroupvelocity) then
            call groupvelocities(k-nlkpt+1,mcv,eigr(1),eigr(maxd2*maxd2+1),maxd2,fscale,lprint)
          endif
!
!  Lower symmetry to remove imaginary modes if selected
!
          call lower(mcv,mcvrptr,freq(1,k-nlkpt+1),nphonatc,nphonatptr,ncfocg,iocptr,maxd2,eigr)
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
              call nagammac(ncfocg,nphonatc,nphonatptr,iocptr,rfmass,bornkloc,maxd2,derv2,dervi,xkt,ykt,zkt)
            else
              bornkpt(1:3) = bornk(1:3,ncf)
              call nagammac(ncfocg,nphonatc,nphonatptr,iocptr,rfmass,bornkpt,maxd2,derv2,dervi,xkt,ykt,zkt)
            endif
          else
            if (lpartofdisp) then
              call nagamma(ncfocg,nphonatc,iocptr,rfmass,bornkloc,maxd2,derv2)
            else
              bornkpt(1:3) = bornk(1:3,ncf)
              call nagamma(ncfocg,nphonatc,iocptr,rfmass,bornkpt,maxd2,derv2)
            endif
          endif
        endif
!
!  Calculate eigenvalues of uncorrected dynamical matrix
!
        if (lintegerKpoint) then
          if (leispack_eigensolve) then
            call pdiage(mcv,maxd2,derv2,dervi,eigr,eigr,freq(1,k-nlkpt+1),fscale,.false.,lprint,ifail)
          else
            call pdiagl(k-nlkpt+1,mcv,maxd2,derv2,dervi,eigr,eigr,freq(1,k-nlkpt+1),fscale,.false.,lprint,ifail)
          endif
        else
          call pdiagg(mcv,maxd2,derv2,eigr,freq(1,k-nlkpt+1),fscale,.false.,lprint,ifail)
        endif
      else
!
!  Calculate eigenvalues of dynamical matrix
!
        if (leispack_eigensolve) then
          call pdiage(mcv,maxd2,derv2,dervi,eigr,eigr,freq(1,k-nlkpt+1),fscale,.false.,lprint,ifail)
        else
          call pdiagl(k-nlkpt+1,mcv,maxd2,derv2,dervi,eigr,eigr,freq(1,k-nlkpt+1),fscale,.false.,lprint,ifail)
        endif
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
        do i = 1,mcv
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
!
!  Aportion kinetic energy to atoms based on eigenvector contribution
!
          trmke = factor*frq*(0.5_dp + w3)
          if (lgamma) then
            do j = 1,mcv
              jj = (j-1)/3_i4 + 1
              projectionfactor = eigr(maxd2*(i-1)+j)**2
              meanKEperatom(jj) = meanKEperatom(jj) + trmke*projectionfactor
            enddo
          else
            do j = 1,mcv
              jj = (j-1)/3_i4 + 1
              eigreal = eigr(maxd2*(i-1)+j)
              eigcomp = eigr(maxd2*maxd2+maxd2*(i-1)+j)
              projectionfactor = eigreal**2 + eigcomp**2
              meanKEperatom(jj) = meanKEperatom(jj) + trmke*projectionfactor
            enddo
          endif
        enddo
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
!**************************************************
!  Deallocate second derivative workspace memory  *
!**************************************************
  deallocate(eigr,stat=status)
  if (status/=0) call deallocate_error('phonon_fcg','eigr')
  if (allocated(savederv2)) then
    if (lintegerKpoint_present) then
      deallocate(savedervi,stat=status)
      if (status/=0) call deallocate_error('phonon_fcg','savedervi')
    endif
    deallocate(savederv2,stat=status)
    if (status/=0) call deallocate_error('phonon_fcg','savederv2')
  endif
!*************************************
!  Output phonon related properties  *
!*************************************
!
!  Compute vibrational energies and output
!
  call vibenergy(mcv,nllkpt,wkpt(nlkpt),freq(1,1),maxfreq,rnokpt,lprint,fc)
!
!  Mean kinetic energies
!
  if (lprinloc.and.lmeanke) then
    write(ioout,'(''  Mean kinetic energy per site: (eV)'')')
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    do i = 1,ncfoc
      write(ioout,'(i6,2x,f15.6)') i,meanKEperatom(i)
    enddo
    write(ioout,'(''--------------------------------------------------------------------------------'')')
  endif
!
!  If output if required then call outphon
!
  if (ldendisp.and.ioproc) then
    call outphon(nlkpt,nukpt,leigloc,3_i4*ncfoc)
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
!
!  Free local memory
!
  if (lthermal) then
    deallocate(Sij,stat=status)
    if (status/=0) call deallocate_error('phonon_fcg','Sij')
  endif
  if (leigloc) then
    deallocate(ramstrength,stat=status)
    if (status/=0) call deallocate_error('phonon_fcg','ramstrength')
    deallocate(oscstrength,stat=status)
    if (status/=0) call deallocate_error('phonon_fcg','oscstrength')
    deallocate(savefreq,stat=status)
    if (status/=0) call deallocate_error('phonon_fcg','savefreq')
  endif
  if (allocated(kpvt)) then
    deallocate(kpvt,stat=status)
    if (status/=0) call deallocate_error('phonon_fcg','kpvt')
  endif
  deallocate(mcvrptr,stat=status)
  if (status/=0) call deallocate_error('phonon_fcg','mcvrptr')
  deallocate(mcvptr,stat=status)
  if (status/=0) call deallocate_error('phonon_fcg','mcvptr')
  deallocate(meanKEperatom,stat=status)
  if (status/=0) call deallocate_error('phonon_fcg','meanKEperatom')
!
!  Timing
!
  t2t = g_cpu_time()
  tphon = t2t - t1t + tphon
#ifdef TRACE
  call trace_out('phonon_fcg')
#endif
!
  return
  end
