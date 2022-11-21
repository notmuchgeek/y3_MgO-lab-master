  subroutine genword(nru,word,lwordok,iline,line,l55,l1000,llibrary,lfflags,ncurr)
!
!  Processes input after keyword line looking for general input
!
!  nru is fortran read channel
!
!   5/00 Created from inword
!   6/01 Initialisation of line added for benefit of some compilers
!   6/01 Modification of VDW radii added
!  10/01 1-D accuracy parameters added
!   5/02 Minimum cell parameter added
!   5/02 Scaling of lower shift rename to avoid conflict
!   8/02 Files added for DLV
!   8/02 Flag for ASCII trajectory file added
!  10/02 Separate option for iseed added
!   1/03 Wolf sum parameters added
!   5/03 XML option added
!   5/03 Potential sites specification added
!   6/03 lmbfgs_order added
!   4/04 Output of pressure file added
!   4/04 Terse option added
!   9/04 Seed for random number generator forced to be negative
!  11/04 Sqrt pi accessed from constants
!   1/05 output sas option added
!   2/05 Potential interpolation option added 
!   4/05 cwolf option added
!   6/05 Fitting flags added for electronegativity input
!   6/05 llibrary flag added as input to genword
!   7/05 Input for Streitz and Mintmire parameters added
!   7/05 Input of number of table points added
!  11/05 Handling of "bin" sub-option in output option modified so that
!        it is only used as such if this is for frequency output
!  11/05 Voter taper option added to cutp
!   3/06 Output file added for oscillator strengths
!   3/06 Bug in adding of .trg to trjfile fixed
!   5/06 Reading of species masses accommodated as well as updating of
!        any previously set species masses
!   7/06 Checking added so that only masses for library species that
!        are relevant are added to list
!   8/06 Radius added to QEq parameter input - note this has been done
!        to try to maintain backwards compatibility
!   8/06 nru passed to linepro
!   8/06 Separate flag added for fitting flags
!   9/06 Literal symbol added to arguments for getpotsymbol1
!  10/06 CML mods added
!   1/07 Gasteiger options added
!   3/07 Radial force added
!   3/07 Gauss renamed to GULP_gauss
!   5/07 Exponential taper added
!   6/07 Option to control some units from the input file added
!  11/07 MDF taper added
!  11/07 Option to output geo file for ReaxFF added
!  12/07 Option to control the frequency of archive frame writes added
!  12/07 730 -> option to set second derivative size remove as it is redundant
!  12/07 Typo reference to qeqtol replaced by qeqscfcrit
!   4/08 Option to input spatial decomposition region size added
!   4/08 Option to input phondiff for phonon finite differences added
!  10/08 Option to prevent overwriting of dumpfiles added
!  10/08 Wolf/COSMIC options merged in
!  12/08 Module input renamed to gulpinput
!   2/09 Modified to accommodate new version of FoX and gulp_cml
!   2/09 Old xml calls removed
!   3/09 Anisotropic rcspatial suboption added
!   3/09 lfinitediff now set when finite option is present
!   5/09 finite difference intervals for property evaluation added
!  11/09 Bug in setting of mass for general element fixed
!  12/09 Handling of units modified for potsites input
!   3/10 Line converted to lower case before checking for end in title option
!   4/10 COSMO file type added
!   8/10 Interpretation of output option modified to allowed for numeric filenames
!  10/10 qbo file output added
!  11/10 Anisotropic pressure added
!   1/11 Force minimisation forced for anisotropic pressure case
!   3/11 Lammps potential file added
!   8/11 Plumed option added
!   3/12 Output of DCD files added
!   5/12 Logic of dump option corrected
!   6/12 Option to select the maths library added
!   6/12 File handling changed to allow for spaces in names
!   7/12 Blocksize added
!   9/12 NMR parameters added for species
!  10/12 Modifications added for lammps file output
!  11/12 Modified Wolf sum of Fennell and Gezelter added
!  11/12 Copy of floatwords to word now restricted to size of word. 
!  12/12 q0 parameters added to electronegativity methods
!   5/13 Option to set the bond order for UFF added
!   6/13 Units of bar added to pressure
!   9/13 Parallel option added
!  10/13 End of line traps added
!  10/13 index_k option added
!  11/13 NMR check on number of floating point numbers corrected
!  12/13 CASTEP phonon output added
!   3/14 Selection of iterative charge solver added
!   3/14 Output of shengBTE files added
!   3/14 derfc changed to g_derfc for benefit of ChemShell
!   3/14 Fractional bond orders added
!   4/14 Iterative charge parameters added
!   1/15 Output of moment of inertia added
!   3/15 cif_dummylattice added
!   3/15 Option added to specify default bond types for species pairs
!   3/15 Input of species specific Gasteiger parameters added
!   3/15 Input of damping parameters for Gasteiger added
!   5/15 kpt file added
!   8/15 Thresholds for third order force constants added
!   8/15 Option to switch minimiser after lower has been used added
!   9/15 Flag for input of stepmax added
!   9/15 Separate stepmax added for rfo
!   1/16 Species specific VDW radius added for cosmo
!   3/16 Plumed options changed
!   5/16 nqkgrid added
!   6/16 nBsplineorder added
!   8/16 Sub-option added to output phonon DOS with IR intensities
!  12/16 Trap for anisotropic pressure with an isotropic cell added
!   2/17 llammpsvar added
!   4/17 Reading of ChemShell option added
!   4/17 Masses and VDW radii now set to element default when adding a new species
!   4/17 Sub-options added for blocksize
!   7/17 accf1D added
!   7/17 nblocksizevar forced to be 3 x nblocksize
!   7/17 Option added to dump to force generated potentials to be output to restart file
!  10/17 Dipole output file added
!   3/18 frqtol added
!   4/18 XSF file added
!   5/18 qrange modifications made
!   6/18 e0range added
!   6/18 switch_stepmx added
!   8/18 7th order polynominal taper added
!  10/18 RFO control variables added and stepmaxrfo renamed 
!  11/18 lowersign added
!   1/19 maxlinelength change added
!   1/19 Worsy changed to use word length parameter
!   2/19 ldump80 added
!   5/19 Finite difference flag split for first and second derivatives
!   6/19 Spin added
!   8/19 Short range damping of polarisation added
!  10/19 Correction to reading of qpolspec
!  10/19 Langevin damping of dipoles added
!  10/19 Error in reading of potsites when at the end of the input fixed
!  12/19 Trap option added
!   2/20 npotsitescfg added
!   3/20 Dielectric constant added
!   3/20 Units option sets inverse_angstroms_to_ev
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
!  Julian Gale, CIC, Curtin University, March 2020
!
  use cellmultipole
  use configurations
  use g_constants
  use control
  use current
  use dump
  use element
  use eemdata
  use gulp_files
  use fitting
  use general
  use genetic,                only : iseed
  use gulpinput
  use gulp_cml,               only : lcml, lvcml, cmlfilename
  use gulpchemsh,             only : ichemsh_output
  use gulp_lengths
  use iochannels
  use kspace,                 only : tweatpi, maxindk, maxindk2, maxindk3, accf1D
  use m_pdfneutron,           only : pdffiles
  use maths,                  only : leispack_eigensolve, ldivide_and_conquer, qsolver
  use maths,                  only : lhess2D, lcosmo2D
  use molecule
  use m_plumed,               only : lplumed, lplumed_available, plumedinput, plumedlog
  use optimisation
  use parallel
  use polarise
  use potentialgrid
  use potentialinterpolation, only : nptsinterpolate, lpotlinterpolate
  use potentialsites
  use radial
  use shells
  use spatial,                only : rcspatial, lrcspatial_anisotropic
  use spatial,                only : rcspatialx, rcspatialy, rcspatialz
  use spatialbo,              only : rcspatialbo, lrcspatialBO_anisotropic
  use spatialbo,              only : rcspatialbox, rcspatialboy, rcspatialboz
  use species
  use spme,                   only : nqkgrid, nBsplineorder
  use terse
  use thresholds
  use trap
  use two
  use uffdata,                only : UFFbondorder
  use wolfcosmo
  implicit none
!
!  Passed variables
!
  character(len=maxwordlength) :: word
  character(len=maxlinelength) :: line
  character(len=maxlinelength) :: linelc
  integer(i4)                  :: iline
  integer(i4)                  :: ncurr
  integer(i4)                  :: nru
  logical                      :: l55
  logical                      :: l1000
  logical                      :: lfflags
  logical                      :: llibrary
  logical                      :: lwordok
!
!  Local variables
!
  character(len=2)             :: cha
  character(len=2)             :: s1
  character(len=5)             :: sym1
  integer(i4)                  :: i
  integer(i4)                  :: iend
  integer(i4)                  :: ii
  integer(i4)                  :: ilp
  integer(i4)                  :: inat
  integer(i4)                  :: ind
  integer(i4)                  :: ind2
  integer(i4)                  :: istart
  integer(i4)                  :: istoan
  integer(i4)                  :: itype
  integer(i4)                  :: itype1
  integer(i4)                  :: itype2
  integer(i4)                  :: iwc
  integer(i4)                  :: j
  integer(i4)                  :: mcal
  integer(i4)                  :: n
  integer(i4)                  :: n1
  integer(i4)                  :: n2
  integer(i4)                  :: n3
  integer(i4)                  :: n4
  integer(i4)                  :: na
  integer(i4)                  :: nbeg
  integer(i4)                  :: ncsword
  integer(i4)                  :: nextraword
  integer(i4)                  :: never
  integer(i4)                  :: nflo
  integer(i4)                  :: nfloatend
  integer(i4)                  :: nfloatused
  integer(i4)                  :: nonzero
  integer(i4)                  :: nqrtyp
  integer(i4)                  :: ns
  integer(i4)                  :: ntot
  integer(i4)                  :: nvar1
  integer(i4)                  :: nvar2
  integer(i4)                  :: nvarshift1
  integer(i4)                  :: nvarshift2
  integer(i4)                  :: nw
  integer(i4)                  :: nw2
  integer(i4)                  :: nwo
  logical                      :: latmassin
  logical                      :: lcart
  logical                      :: lever
  logical                      :: lfound
  logical                      :: lin
  logical                      :: llocalseq
  logical                      :: llocalsym
  logical                      :: llocalmovie
  logical                      :: lnodot
  logical                      :: lok1
  logical                      :: lout
  logical                      :: lsymbol
  logical                      :: lset
  logical                      :: ltype01
  logical                      :: lvalid
  logical                      :: lvdwin
  real(dp)                     :: g_derfc
  real(dp)                     :: rd5
  real(dp)                     :: rd7
  real(dp)                     :: ta
  real(dp)                     :: ta2
  real(dp)                     :: taperrange
  real(dp)                     :: tb
  real(dp)                     :: tb2
  real(dp)                     :: units
  real(dp)                     :: val
!
!  Options
!
  if (index(word,'accu').eq.1) goto 100
  if (index(word,'qele').eq.1) goto 110
  if (index(word,'cmm').eq.1)  goto 120
  if (index(word,'elec').eq.1) goto 130
  if (index(word,'name').eq.1) goto 140
  if (index(word,'cutb').eq.1) goto 150
  if (index(word,'cutd').eq.1) goto 150
  if (index(word,'cutp').eq.1) goto 160
  if (index(word,'site').eq.1) goto 170
  if (index(word,'seed').eq.1) goto 180
  if (index(word,'ewal').eq.1) goto 190
  if (index(word,'slow').eq.1) goto 200
  if (index(word,'pola').eq.1) goto 210
  if (index(word,'delt').eq.1) goto 220
  if (index(word,'maxc').eq.1) goto 230
  if (index(word,'step').eq.1) goto 240
  if (index(word,'xtol').eq.1) goto 250
  if (index(word,'switch_s').eq.1) goto 255
  if (index(word,'swit').eq.1) goto 260
  if (index(word,'potg').eq.1) goto 270
  if (index(word,'time ').eq.1) goto 280
  if (index(word,'qeqt').eq.1) goto 290
  if (index(word,'qeqi').eq.1) goto 300
  if (index(word,'qeqr').eq.1) goto 310
  if (index(word,'titl').eq.1) goto 320
  if (index(word,'rtol').eq.1) goto 330
  if (index(word,'cdum').eq.1) goto 340
  if (index(word,'dump').eq.1) goto 350
  if (index(word,'fini').eq.1) goto 360
  if (index(word,'cuts').eq.1) goto 370
  if (index(word,'minc').eq.1) goto 380
  if (index(word,'gtol').eq.1) goto 390
  if (index(word,'gmax').eq.1) goto 400
  if (index(word,'ters').eq.1) goto 410
  if (index(word,'ftol').eq.1) goto 420
  if (index(word,'prin').eq.1) goto 430
  if (index(word,'elem').eq.1) goto 440
  if (index(word,'mass').eq.1) goto 445
  if (index(word,'ioni').eq.1) goto 445
  if (index(word,'cova').eq.1) goto 445
  if (index(word,'vdw').eq.1)  goto 445
  if (index(word,'nmr').eq.1)  goto 445
  if (index(word,'symb').eq.1) goto 445
  if (index(word,'bbar').eq.1) goto 445
  if (index(word,'sigi').eq.1) goto 445
  if (index(word,'spin').eq.1) goto 445
  if (index(word,'smel').eq.1) goto 490
  if (index(word,'maxl').eq.1) goto 500
  if (index(word,'outp').eq.1) goto 510
  if (index(word,'rspe').eq.1) goto 520
  if (index(word,'nadd').eq.1) goto 530
  if (index(word,'gdcr').eq.1) goto 540
  if (index(word,'nobo').eq.1) goto 550
  if (index(word,'upda').eq.1) goto 560
  if (index(word,'lbfg').eq.1) goto 570
  if (index(word,'line').eq.1) goto 580
  if (index(word,'pots').eq.1) goto 590
  if (index(word,'pote').eq.1) goto 600
  if (index(word,'delf').eq.1) goto 610
  if (index(word,'marv').eq.1) goto 620
  if (index(word,'qmmm').eq.1) goto 630
  if (index(word,'unit').eq.1) goto 640
  if (index(word,'math').eq.1) goto 650
  if (index(word,'para').eq.1) goto 660
  if (index(word,'qsol').eq.1) goto 670
  if (index(word,'stre').eq.1) goto 680
  if (index(word,'maxi').eq.1) goto 690
  if (index(word,'pres').eq.1) goto 700
  if (index(word,'qwol').eq.1) goto 710
  if (index(word,'radi').eq.1) goto 720
  if (index(word,'rcsp').eq.1) goto 730
  if (index(word,'gastt').eq.1) goto 740
  if (index(word,'gasti').eq.1) goto 750
  if (index(word,'pfin').eq.1) goto 760
  if (index(word,'qwol').eq.1) goto 770
  if (index(word,'cwol').eq.1) goto 780
  if (index(word,'sfin').eq.1) goto 790
  if (index(word,'anis').eq.1) goto 800
  if (index(word,'plumed_l').eq.1) goto 805
  if (index(word,'plum').eq.1) goto 810
  if (index(word,'bloc').eq.1) goto 820
  if (index(word,'uff_b').eq.1) goto 825
  if (index(word,'inde').eq.1) goto 830
  if (index(word,'qite').eq.1) goto 840
  if (index(word,'bond').eq.1) goto 850
  if (index(word,'gastp').eq.1) goto 860
  if (index(word,'gastd').eq.1) goto 870
  if (index(word,'thres').eq.1) goto 880
  if (index(word,'qgri').eq.1) goto 890
  if (index(word,'bspl').eq.1) goto 900
  if (index(word,'chem').eq.1) goto 910
  if (index(word,'matr').eq.1) goto 920
  if (index(word,'frqt').eq.1) goto 930
  if (index(word,'rfo_e').eq.1) goto 940
  if (index(word,'rfo_g').eq.1) goto 950
  if (index(word,'diel').eq.1) goto 960
  if (index(word,'trap').eq.1) goto 970
  return
!*************************************
!  Accuracy factor for lattice sums  *
!*************************************
100 if (nfloat.gt.0) then
    accuracy = floats(1)
    if (nfloat.ge.2) then
      nemorder = nint(floats(2))
    endif
    if (nfloat.ge.3) then
      nmaxcells = nint(floats(3))
    endif
    if (nfloat.ge.4) then
      accf1D = abs(floats(4))
    endif
  else
    read(nru,*,err=99,end=99) accuracy
    iline = iline + 1
  endif
  lwordok = .true.
  return
!**************************
!  Cell multipole method  *
!**************************
120 icmm = 3
  if (nfloat.gt.0) then
    rbox = floats(1)
  endif
  if (nword.gt.1) then
    do i = 2,nword
      if (index(words(i),'mo').eq.1) icmm = 1
      if (index(words(i),'di').eq.1) icmm = 2
      if (index(words(i),'qu').eq.1) icmm = 3
      if (index(words(i),'oc').eq.1) icmm = 4
    enddo
  endif
  lwordok = .true.
  return
!*********************************
!  Electronegativity parameters  *
!*********************************
110 units = 1.0_dp
  nqrtyp = 0
  if (nword.gt.1) then
    do i = 2,nword
      call stolc(words(i),maxword)
      if (index(words(i),'au').eq.1) units = autoangs
      if (index(words(i),'qmi').eq.1) nqrtyp = 1
      if (index(words(i),'qma').eq.1) nqrtyp = 2
      if (index(words(i),'qra').eq.1) nqrtyp = 3
    enddo
  endif
!
!  Start of input loop
!
115 line = '  '
  read(nru,'(a)',end=1000) line
  iline = iline + 1
  call linepro(nru,line,iline)
!
!  Check whether input is symbol or option
!
  if ((nword+nfloat).eq.0) goto 115
!
!  Check for old fashion specification of number of atoms
!
  if (nword.gt.0) then
    word = words(1)
    call worsy(word,lsymbol,.false.)
    if (.not.lsymbol) then
      l55 = .true.
      return
    endif
  endif
!
!  Find atom number for element
!
!  Symbols used in input
!
  if (nword.gt.0) then
    call ltont(words(1),inat,itype)
    leemparaltered(inat) = .true.
    nfloatused = 0
  else
    inat = int(floats(1))
    leemparaltered(inat) = .true.
    nfloatused = 1
  endif
!
!  Set end of float range for given qrange type
!
  nfloatend = nfloat
  if (nqrtyp.eq.3) then
    if (lfit.and.lfflags) then
      nfloatend = min(nfloatend,nfloatused+10_i4)
      if (nfloatend.lt.6_i4) then
        call outerror('insufficient input parameters for electronegativity',iline)
        call stopnow('genword')
      endif
    else
      nfloatend = min(nfloatend,nfloatused+7_i4)
      if (nfloatend.lt.3_i4) then
        call outerror('insufficient input parameters for electronegativity',iline)
        call stopnow('genword')
      endif
    endif
  elseif (nqrtyp.eq.2.or.nqrtyp.eq.1) then
    if (lfit.and.lfflags) then
      nfloatend = min(nfloatend,nfloatused+9_i4)
      if (nfloatend.lt.5_i4) then
        call outerror('insufficient input parameters for electronegativity',iline)
        call stopnow('genword')
      endif
    else
      nfloatend = min(nfloatend,nfloatused+6_i4)
      if (nfloatend.lt.2_i4) then
        call outerror('insufficient input parameters for electronegativity',iline)
        call stopnow('genword')
      endif
    endif
  else
    if (lfit.and.lfflags) then
      nfloatend = min(nfloatend,nfloatused+8_i4)
      if (nfloatend.lt.4_i4) then
        call outerror('insufficient input parameters for electronegativity',iline)
        call stopnow('genword')
      endif
    else
      nfloatend = min(nfloatend,nfloatused+5_i4)
      if (nfloatend.lt.1_i4) then
        call outerror('insufficient input parameters for electronegativity',iline)
        call stopnow('genword')
      endif
    endif
  endif
!
!  Increment number of ranges
!
  nqrange(inat,2) = nqrange(inat,2) + 1
!
!  Check dimensions
!
  if (nqrange(inat,2).gt.maxqrange) then
    maxqrange = nqrange(inat,2)
    call changemaxqrange
  endif
!
!  Assign range type 
!
  nqrangetype(nqrange(inat,2),inat,2) = nqrtyp
!    
!  Fitting flags
!
  if (lfit.and.lfflags) then
    n1 = int(floats(nfloatend-2))
    n2 = int(floats(nfloatend-1))
    n3 = int(floats(nfloatend))
    if (n1.eq.1) then
      nfit = nfit + 1
      if (nfit.ge.maxfit) then
        maxfit = nfit + 10
        call changemaxfit
      endif
      nftyp(nfit) = 1
      nfpot(nfit) = 11
      nfvar(nfit) = inat
      nfvar2(nfit) = nqrange(inat,2)
      nfvar3(nfit) = 2
    endif
    if (n2.eq.1) then
      nfit = nfit + 1
      if (nfit.ge.maxfit) then
        maxfit = nfit + 10
        call changemaxfit
      endif
      nftyp(nfit) = 1
      nfpot(nfit) = 12
      nfvar(nfit) = inat
      nfvar2(nfit) = nqrange(inat,2)
      nfvar3(nfit) = 2
    endif
    if (n3.eq.1) then
      nfit = nfit + 1
      if (nfit.ge.maxfit) then
        maxfit = nfit + 10
        call changemaxfit
      endif
      nftyp(nfit) = 1
      nfpot(nfit) = 13
      nfvar(nfit) = inat
      nfvar2(nfit) = nqrange(inat,2)
      nfvar3(nfit) = 2
    endif
    nfloatend = nfloatend - 3
  endif
!
!  Assign range parameters
!
  if (nqrtyp.eq.3) then
    qrangemax(nqrange(inat,2),inat,2) = floats(nfloatend)
    qrangemin(nqrange(inat,2),inat,2) = floats(nfloatend-1)
    nfloatend = nfloatend - 2
  elseif (nqrtyp.eq.2) then
    qrangemax(nqrange(inat,2),inat,2) = floats(nfloatend)
    nfloatend = nfloatend - 1
  elseif (nqrtyp.eq.1) then
    qrangemin(nqrange(inat,2),inat,2) = floats(nfloatend)
    nfloatend = nfloatend - 1
  endif
!
!  Assign electronegativity parameters
!
  nfloatused = nfloatused + 1
  chirange(nqrange(inat,2),inat,2) = floats(nfloatused)
  if (nfloatused.lt.nfloatend) then
    nfloatused = nfloatused + 1
    murange(nqrange(inat,2),inat,2) = floats(nfloatused)
    if (nfloatused.lt.nfloatend) then
      nfloatused = nfloatused + 1
      radrange(nqrange(inat,2),inat,2) = floats(nfloatused)
      if (nfloatused.lt.nfloatend) then
        nfloatused = nfloatused + 1
        q0range(nqrange(inat,2),inat,2) = floats(nfloatused)
        if (nfloatused.lt.nfloatend) then
          nfloatused = nfloatused + 1
          e0range(nqrange(inat,2),inat,2) = floats(nfloatused)
        endif
      endif
    endif
  endif
  goto 115
!*********************************
!  Electronegativity parameters  *
!*********************************
130 units = 1.0_dp
  nqrtyp = 0
  if (nword.gt.1) then
    do i = 2,nword
      call stolc(words(i),maxword)
      if (index(words(i),'au').eq.1) units = autoangs
      if (index(words(i),'qmi').eq.1) nqrtyp = 1
      if (index(words(i),'qma').eq.1) nqrtyp = 2
      if (index(words(i),'qra').eq.1) nqrtyp = 3
    enddo
  endif
!
!  Start of input loop
!
135 line = '  '
  read(nru,'(a)',end=1000) line
  iline = iline + 1
  call linepro(nru,line,iline)
!
!  Check whether input is symbol or option
!
  if ((nword+nfloat).eq.0) goto 135
!
!  Check for old fashion specification of number of atoms
!
  if (nword.gt.0) then
    word = words(1)
    call worsy(word,lsymbol,.false.)
    if (.not.lsymbol) then
      l55 = .true.
      return
    endif
  endif
!
!  Find atom number for element
!
!  Symbols used in input
!
  if (nword.gt.0) then
    call ltont(words(1),inat,itype)
    leemparaltered(inat) = .true.
    nfloatused = 0
  else
    inat = int(floats(1))
    leemparaltered(inat) = .true.
    nfloatused = 1
  endif
!
!  Set end of float range for given qrange type
!
  nfloatend = nfloat
  if (nqrtyp.eq.3) then
    if (lfit.and.lfflags) then
      nfloatend = min(nfloatend,nfloatused+8_i4)
      if (nfloatend.lt.5_i4) then
        call outerror('insufficient input parameters for electronegativity',iline)
        call stopnow('genword')
      endif
    else
      nfloatend = min(nfloatend,nfloatused+6_i4)
      if (nfloatend.lt.3_i4) then
        call outerror('insufficient input parameters for electronegativity',iline)
        call stopnow('genword')
      endif
    endif
  elseif (nqrtyp.eq.2.or.nqrtyp.eq.1) then
    if (lfit.and.lfflags) then
      nfloatend = min(nfloatend,nfloatused+7_i4)
      if (nfloatend.lt.4_i4) then
        call outerror('insufficient input parameters for electronegativity',iline)
        call stopnow('genword')
      endif
    else
      nfloatend = min(nfloatend,nfloatused+5_i4)
      if (nfloatend.lt.2_i4) then
        call outerror('insufficient input parameters for electronegativity',iline)
        call stopnow('genword')
      endif
    endif
  else
    if (lfit.and.lfflags) then
      nfloatend = min(nfloatend,nfloatused+6_i4)
      if (nfloatend.lt.3_i4) then
        call outerror('insufficient input parameters for electronegativity',iline)
        call stopnow('genword')
      endif
    else
      nfloatend = min(nfloatend,nfloatused+4_i4)
      if (nfloatend.lt.1_i4) then
        call outerror('insufficient input parameters for electronegativity',iline)
        call stopnow('genword')
      endif
    endif
  endif
!
!  Increment number of ranges
!
  nqrange(inat,1) = nqrange(inat,1) + 1
!
!  Check dimensions
!
  if (nqrange(inat,1).gt.maxqrange) then
    maxqrange = nqrange(inat,1)
    call changemaxqrange
  endif
!
!  Assign range type 
!
  nqrangetype(nqrange(inat,1),inat,1) = nqrtyp
!    
!  Fitting flags
!
  if (lfit.and.lfflags) then
    n1 = int(floats(nfloatend-1))
    n2 = int(floats(nfloatend))
    if (n1.eq.1) then
      nfit = nfit + 1
      if (nfit.ge.maxfit) then
        maxfit = nfit + 10
        call changemaxfit
      endif
      nftyp(nfit) = 1
      nfpot(nfit) = 11
      nfvar(nfit) = inat
      nfvar2(nfit) = nqrange(inat,1)
      nfvar3(nfit) = 1
    endif
    if (n2.eq.1) then
      nfit = nfit + 1
      if (nfit.ge.maxfit) then
        maxfit = nfit + 10
        call changemaxfit
      endif
      nftyp(nfit) = 1
      nfpot(nfit) = 12
      nfvar(nfit) = inat
      nfvar2(nfit) = nqrange(inat,1)
      nfvar3(nfit) = 1
    endif
    nfloatend = nfloatend - 2
  endif
!
!  Assign range parameters
!
  if (nqrtyp.eq.3) then
    qrangemax(nqrange(inat,1),inat,1) = floats(nfloatend)
    qrangemin(nqrange(inat,1),inat,1) = floats(nfloatend-1)
    nfloatend = nfloatend - 2
  elseif (nqrtyp.eq.2) then
    qrangemax(nqrange(inat,1),inat,1) = floats(nfloatend)
    nfloatend = nfloatend - 1
  elseif (nqrtyp.eq.1) then
    qrangemin(nqrange(inat,1),inat,1) = floats(nfloatend)
    nfloatend = nfloatend - 1
  endif
!
!  Assign electronegativity parameters
!
  nfloatused = nfloatused + 1
  chirange(nqrange(inat,1),inat,1) = floats(nfloatused)
  if (nfloatused.lt.nfloatend) then
    nfloatused = nfloatused + 1
    murange(nqrange(inat,1),inat,1) = floats(nfloatused)
    if (nfloatused.lt.nfloatend) then
      nfloatused = nfloatused + 1
      q0range(nqrange(inat,1),inat,1) = floats(nfloatused)
      if (nfloatused.lt.nfloatend) then
        nfloatused = nfloatused + 1
        e0range(nqrange(inat,1),inat,1) = floats(nfloatused)
      endif
    endif
  endif
  goto 135
!*******************
!  Structure name  *
!*******************
140 if (ncfg+1.gt.maxcfg) then
    maxcfg = ncfg + 1
    call changemaxcfg
  endif
  if (nword.eq.1) then
    line = '  '
    read(nru,'(a)') line
    iline = iline + 1
    call linepro(nru,line,iline)
    names(ncfg+1) = words(1)
  else
    names(ncfg+1) = words(2)
  endif
  lwordok = .true.
  return
!*************************************
!  Cutoff for distance calculations  *
!*************************************
150 units = 1.0_dp
  if (nword.gt.1) then
    call stolc(words(2),maxword)
    if (index(words(2),'au').eq.1) units = autoangs
  endif
  if (nfloat.gt.0) then
    cutb = floats(1)*units
  else
    line = '  '
    read(nru,'(a)') line
    iline = iline + 1
    call linepro(nru,line,iline)
    cutb = floats(1)
    if (nword.ge.1) then
      call stolc(words(1),maxword)
      if (index(words(1),'au').eq.1) units = autoangs
    endif
    cutb = cutb*units
  endif
  lwordok = .true.
  return
!**************************************************
!  Cutoff for interatomic potential calculations  *
!**************************************************
160 units = 1.0_dp
  if (nword.gt.1) then
    do ii = 2,nword
      call stolc(words(ii),maxword)
      if (index(words(ii),'au').eq.1) then
        units = autoangs
      elseif (index(words(ii),'po').eq.1) then
        tapertype = 1
      elseif (index(words(ii),'co').eq.1) then
        tapertype = 2
      elseif (index(words(ii),'vo').eq.1) then
        tapertype = 3
      elseif (index(words(ii),'ex').eq.1) then
        tapertype = 4
      elseif (index(words(ii),'md').eq.1) then
        tapertype = 5
      elseif (index(words(ii),'p7').eq.1) then
        tapertype = 6
      endif
    enddo
  endif
  if (tapertype.eq.3) then
!
!  Voter style taper
!
    if (nfloat.gt.0) then
      cutp = floats(1)*units
      if (nfloat.gt.1) then
        taperm = abs(floats(2))
      endif
    else
      line = '  '
      read(nru,'(a)') line
      iline = iline + 1
      call linepro(nru,line,iline)
      cutp = floats(1)
      if (nfloat.gt.1) then
        taperm = abs(floats(2))
      endif
      if (nword.ge.1) then
        call stolc(words(1),maxword)
        if (index(words(1),'au').eq.1) units = autoangs
      endif
      cutp = cutp*units
    endif
    if (taperm.lt.1.0d-12) then
      call outerror('Voter taper m cannot be zero',iline)
      call stopnow('genword')
    endif
    taperrange = cutp
  elseif (tapertype.eq.4.or.tapertype.eq.6) then
!
!  Exponential or 7th order polynominal style taper
!
    if (nfloat.gt.0) then
      cutp = floats(1)*units
    else
      line = '  '
      read(nru,'(a)') line
      iline = iline + 1
      call linepro(nru,line,iline)
      cutp = floats(1)
      if (nword.ge.1) then
        call stolc(words(1),maxword)
        if (index(words(1),'au').eq.1) units = autoangs
      endif
      cutp = cutp*units
    endif
    taperrange = cutp
  else
    if (nfloat.gt.0) then
      cutp = floats(1)*units
      if (nfloat.gt.1) then
        taperrange = abs(floats(2))*units
      endif
    else
      line = '  '
      read(nru,'(a)') line
      iline = iline + 1
      call linepro(nru,line,iline)
      cutp = floats(1)
      if (nfloat.gt.1) then
        taperrange = abs(floats(2))
      endif
      if (nword.ge.1) then
        call stolc(words(1),maxword)
        if (index(words(1),'au').eq.1) units = autoangs
      endif
      cutp = cutp*units
      taperrange = taperrange*units
    endif
  endif
!
!  Set up taper constants
!
  if (taperrange.gt.1.0d-2) then
    tapermax = cutp
    tapermin = cutp - taperrange
    tb = cutp
    tb2 = tb*tb
    ta = tb - taperrange
    ta2 = ta*ta
    rd5 = 1.0_dp/(taperrange**5.0_dp)
    pts0 = (10.0_dp*ta2*tb+(tb-5.0_dp*ta)*tb2)*tb2*rd5
    pts1 = - 30.0_dp*ta2*tb2*rd5
    pts2 = 30.0_dp*(ta2*tb+ta*tb2)*rd5
    pts3 = - 10.0_dp*(ta2+4.0_dp*ta*tb+tb2)*rd5
    pts4 = 15.0_dp*(ta+tb)*rd5
    pts5 = - 6.0_dp*rd5
    if (tapertype.eq.6) then
      rd7 = 1.0_dp/(tb**7.0_dp)
      pts7_6 =  20.0_dp*rd7
      pts6_6 = -70.0_dp*tb*rd7
      pts5_6 =  84.0_dp*tb2*rd7
      pts4_6 = -35.0_dp*tb*tb2*rd7
      pts3_6 =   0.0_dp
      pts2_6 =   0.0_dp
      pts1_6 =   0.0_dp
      pts0_6 =   1.0_dp
    endif
  endif
  lwordok = .true.
  return
!*********************
!  Site list option  *
!*********************
170 if (ioproc) then
    write(ioout,'(''  **** Site list option has been withdrawn ****'')')
  endif
  lwordok = .true.
  return
!*************************************
!  Seed for random number generator  *
!*************************************
180 if (nfloat.ge.1) then  
    iseed = nint(floats(1))
  else  
    read(nru,*,err=99,end=99) iseed
    iline = iline + 1
  endif
  iseed = - abs(iseed)
  lwordok = .true.
  return
!**************************************************
!  Target real space radius for Ewald/Parry sums  *
!**************************************************
190 if (nfloat.ge.1) then  
    targetrradmax = abs(floats(1))
  else  
    read(nru,*,err=99,end=99) targetrradmax
    iline = iline + 1
  endif
  lwordok = .true.
  return
!*****************************************************
!  Change default value of scaling for lower option  *
!*****************************************************
200 continue
  if (nfloat.gt.0) then
    lowerscale = floats(1)
    if (nfloat.gt.1) then
      if (floats(2).lt.0.0_dp) then
        lowersign = - 1.0_dp
      endif
    endif
  else
    read(nru,*,err=99,end=99) lowerscale
    iline = iline + 1
  endif
  lwordok = .true.
  return
!************************
!  Polarisability data  *
!************************
210 lpolar = .true.
  if (nword.gt.1) then
    do ii = 2,nword
      call stolc(words(ii),maxword)
      if (index(words(ii),'da').eq.1) then
        lpoldamp = .true.
      elseif (index(words(ii),'la').eq.1) then
        lpollangevin = .true.
      endif
    enddo
  endif
  if (lpoldamp.and.nfloat.gt.0) then
    bpdamp = abs(floats(1))
  endif
215 line = '  '
  read(nru,'(a)',end=218) line
  iline = iline + 1
  call linepro(nru,line,iline)
!
!  Check whether input is symbol or option
!
  if ((nword+nfloat).eq.0) goto 215
!
!  Check for old fashion specification of number of atoms
!
  if (nword.eq.0.and.nfloat.eq.1) goto 215
  if (nword.gt.0) then
    word = words(1)
    call worsy(word,lsymbol,.false.)
    if (.not.lsymbol) then
      call stolc(word,maxwordlength)
      if (index(word,'end').eq.1) then
        lwordok = .true.
        return
      else
        l55 = .true.
        goto 218
      endif
    endif
  endif
  npolspec = npolspec + 1
  if (npolspec.gt.maxspec) then
    maxspec = npolspec + 20
    call changemaxspec
  endif
  if (lpollangevin) then
    if (nword.eq.0) then
      inat = int(floats(1))
      if (inat.gt.100) inat = inat - 100
      natpolspec(npolspec) = inat
      ntyppolspec(npolspec) = 0
      dpolspec(npolspec) = floats(2)
      if (nfloat.gt.2) then
        dpolmaxspec(npolspec) = abs(floats(3))
        if (dpolmaxspec(npolspec).lt.1.0d-12) then
          call outerror('Saturation dipole value is too small',iline)
          call stopnow('genword')
        endif
      else
        call outerror('Saturation dipole value missing',iline)
        call stopnow('genword')
      endif
    elseif (nword.ge.1) then
      call ltont(words(1),inat,itype)
      natpolspec(npolspec) = inat
      ntyppolspec(npolspec) = itype
      dpolspec(npolspec) = floats(1)
      if (nfloat.gt.1) then
        dpolmaxspec(npolspec) = abs(floats(2))
        if (dpolmaxspec(npolspec).lt.1.0d-12) then
          call outerror('Saturation dipole value is too small',iline)
          call stopnow('genword')
        endif
      else
        call outerror('Saturation dipole value missing',iline)
        call stopnow('genword')
      endif
    endif
  else
    if (nword.eq.0) then
      inat = int(floats(1))
      if (inat.gt.100) inat = inat - 100
      natpolspec(npolspec) = inat
      ntyppolspec(npolspec) = 0
      dpolspec(npolspec) = floats(2)
      if (nfloat.gt.2) then
        qpolspec(npolspec) = floats(3)
      endif
    elseif (nword.ge.1) then
      call ltont(words(1),inat,itype)
      natpolspec(npolspec) = inat
      ntyppolspec(npolspec) = itype
      dpolspec(npolspec) = floats(1)
      if (nfloat.gt.1) then
        qpolspec(npolspec) = floats(2)
      endif
    endif
  endif
  goto 215
218 if (.not.l55) l1000 = .true.
  lwordok = .true.
  return
!**********************************
!  Change default value of delta  *
!**********************************
220 if (nfloat.gt.0) then
    delta = floats(1)
  else
    read(nru,*,err=99,end=99) delta
    iline = iline + 1
  endif
!
  delta = abs(delta)
  if (delta.ge.1.0_dp) delta = 10.0**(-delta)
  lwordok = .true.
  return
!***********************************
!  Change default value of maxcal  *
!***********************************
230 if (nword.gt.1) then
    if (index(words(2),'opt').eq.1) then
!
!  Set only optimisation maxcal
!
      if (nfloat.gt.0) then
        maxcal = int(floats(1))
      else
        read(nru,*,err=99,end=99) maxcal
        iline = iline + 1
      endif
    elseif (index(words(2),'fit').eq.1) then
!
!  Set only fit maxcal
!
      if (nfloat.gt.0) then
        maxfcal = int(floats(1))
      else
        read(nru,*,err=99,end=99) maxfcal
        iline = iline + 1
      endif
    else
!
!  Set both
!
      if (nfloat.gt.0) then
        mcal = int(floats(1))
      else
        read(nru,*,err=99,end=99) mcal
        iline = iline + 1
      endif
      if (lfit) maxfcal = mcal
      if (lopt) maxcal = mcal
    endif
  else
!
!  Set both
!
    if (nfloat.gt.0) then
      maxcal = int(floats(1))
    else
      read(nru,*,err=99,end=99) maxcal
      iline = iline + 1
    endif
    maxfcal = maxcal
  endif
  lwordok = .true.
  return
!************************************
!  Change default value of stepmax  *
!************************************
240 if (nword.gt.1) then
    if (index(words(2),'opt').eq.1) then
!
!  Set only optimisation stepmax
!
      lstepmaxin = .true.
      if (nfloat.gt.0) then
        stepmax = floats(1)
      else
        read(nru,*,err=99,end=99) stepmax
        iline = iline + 1
      endif
    elseif (index(words(2),'rfo').eq.1) then
!
!  Set only optimisation stepmax for RFO 
!
      if (nfloat.gt.0) then
        rfostepmax = floats(1)
      else
        read(nru,*,err=99,end=99) rfostepmax
        iline = iline + 1
      endif
    elseif (index(words(2),'fit').eq.1) then
!
!  Set only fit stepmax
!
      if (nfloat.gt.0) then
        fstepmx = floats(1)
      else
        read(nru,*,err=99,end=99) fstepmx
        iline = iline + 1
      endif
    else
!
!  Set according to job type
!
      if (lopt) then
        lstepmaxin = .true.
        if (nfloat.gt.0) then
          stepmax = floats(1)
        else
          read(nru,*,err=99,end=99) stepmax
          iline = iline + 1
        endif
        if (lfit) fstepmx = stepmax
      elseif (lfit) then
        if (nfloat.gt.0) then
          fstepmx = floats(1)
        else
          read(nru,*,err=99,end=99) fstepmx
          iline = iline + 1
        endif
      endif
    endif
  else
!
!  Set according to job type
!
    if (lopt) then
      lstepmaxin = .true.
      if (nfloat.gt.0) then
        stepmax = floats(1)
      else
        read(nru,*,err=99,end=99) stepmax
        iline = iline + 1
      endif
      if (lfit) fstepmx = stepmax
    elseif (lfit) then
      if (nfloat.gt.0) then
        fstepmx = floats(1)
      else
        read(nru,*,err=99,end=99) fstepmx
        iline = iline + 1
      endif
    endif
  endif
  lwordok = .true.
  return
!*********************************
!  Change default value of xtol  *
!*********************************
250 if (nword.gt.1) then
    if (index(words(2),'opt').eq.1) then
!
!  Set only optimisation xtol
!
      if (nfloat.gt.0) then
        xtol = floats(1)
      else
        read(nru,*,err=99,end=99) xtol
        iline = iline + 1
      endif
    elseif (index(words(2),'fit').eq.1) then
!
!  Set only fit xtol
!
      if (nfloat.gt.0) then
        fxtol = floats(1)
      else
        read(nru,*,err=99,end=99) fxtol
        iline = iline + 1
      endif
    else
!
!  Set according to job type
!
      if (lopt) then
        if (nfloat.gt.0) then
          xtol = floats(1)
        else
          read(nru,*,err=99,end=99) xtol
          iline = iline + 1
        endif
        if (lfit) fxtol = xtol
      elseif (lfit) then
        if (nfloat.gt.0) then
          fxtol = floats(1)
        else
          read(nru,*,err=99,end=99) fxtol
          iline = iline + 1
        endif
      endif
    endif
  else
!
!  Set according to job type
!
    if (lopt) then
      if (nfloat.gt.0) then
        xtol = floats(1)
      else
        read(nru,*,err=99,end=99) xtol
        iline = iline + 1
      endif
      if (lfit) fxtol = xtol
    elseif (lfit) then
      if (nfloat.gt.0) then
        fxtol = floats(1)
      else
        read(nru,*,err=99,end=99) fxtol
        iline = iline + 1
      endif
    endif
  endif
  xtol = abs(xtol)
  if (xtol.ge.1) xtol = 10.0**(-xtol)
  fxtol = abs(fxtol)
  if (fxtol.ge.1) fxtol = 10.0**(-fxtol)
  lwordok = .true.
  return
!*******************
!  Switch stepmax  *
!*******************
!
!  lmstpch    = logical controlling whether change is allowed
!  mstpchcrit = indicator to change criteria
!              1 => no. of cycles
!              2 => gradient norm
!  stpchcrit  = change criteria (no. of cycles / gnorm)
!  stepmxch   = new step max size
!
255 if (nfloat+nword.eq.1) then
    line = '  '
    read(nru,'(a)',end=1000) line
    iline = iline + 1
    call linepro(nru,line,iline)
    istart = 1
  else
    istart = 2
  endif
  mstpchcrit = 2
  lmstpch = .true.
!
!  Find option words
!
  do i = istart,nword
    if (index(words(i),'cyc').eq.1) then
      minchcrit = 1
    elseif (index(words(i),'gno').eq.1) then
      minchcrit = 2
    endif
  enddo
  if (nfloat.lt.2) then
    call outerror('insufficient values given for switch_stepmx option',iline)
    call stopnow('genword')
  endif
  stepmxch = abs(floats(1))
  stpchcrit = abs(floats(2))
  lwordok = .true.
  return
!*********************
!  Switch minimiser  *
!*********************
!
!  lminch    = logical controlling whether change is allowed
!  mintype   = index for new minimiser
!              1 => exact hessian + BFGS
!              2 => exact hessian + RFO
!              3 => unit hessian + BFGS
!              4 => numerical (numerical diagonal)
!              5 => conjugate gradients
!              6 => limited memory BFGS
!  minchcrit = indicator to change criteria
!              1 => no. of cycles
!              2 => gradient norm
!              3 => if lower has been used
!  chcrit    = change criteria (no. of cycles / gnorm)
!
260 if (nfloat+nword.eq.1) then
    line = '  '
    read(nru,'(a)',end=1000) line
    iline = iline + 1
    call linepro(nru,line,iline)
    istart = 1
  else
    istart = 2
  endif
  minchcrit = 2
  mintype = 0
  lminch = .true.
!
!  Find option words
!
  do i = istart,nword
    if (index(words(i),'bfgs').eq.1) then
      mintype = 1
    elseif (index(words(i),'rfo').eq.1) then
      mintype = 2
    elseif (index(words(i),'unit').eq.1) then
      mintype = 3
    elseif (index(words(i),'nume').eq.1) then
      mintype = 4
    elseif (index(words(i),'conj').eq.1) then
      mintype = 5
    elseif (index(words(i),'lbfg').eq.1) then
      mintype = 6
    elseif (index(words(i),'cyc').eq.1) then
      minchcrit = 1
    elseif (index(words(i),'low').eq.1) then
      minchcrit = 3
    endif
  enddo
  if (mintype.eq.0) then
    call outerror('no new minimiser type given in switch option',iline)
    call stopnow('genword')
  endif
  if (nfloat.eq.0.and.minchcrit.ne.3) then
    call outerror('no change criterion given in switch option',iline)
    call stopnow('genword')
  endif
  chcrit = floats(1)
  lwordok = .true.
  return
!************************************
!  Potential on a grid calculation  *
!************************************
270 if (nfloat.eq.0) then
    line = '  '
    read(nru,'(a)') line
    iline = iline + 1
    call linepro(nru,line,iline)
  endif
  if (nfloat.ge.9) then
    xminpg(ncurr) = floats(1)
    xmaxpg(ncurr) = floats(2)
    yminpg(ncurr) = floats(3)
    ymaxpg(ncurr) = floats(4)
    zminpg(ncurr) = floats(5)
    zmaxpg(ncurr) = floats(6)
    nxpg(ncurr) = nint(floats(7))
    nypg(ncurr) = nint(floats(8))
    nzpg(ncurr) = nint(floats(9))
  elseif (nfloat.ge.3) then
    nxpg(ncurr) = nint(floats(1))
    nypg(ncurr) = nint(floats(2))
    nzpg(ncurr) = nint(floats(3))
  else
    call outerror('insufficient input for potgrid option',iline)
    call stopnow('genword')
  endif
  if (nxpg(ncurr).eq.0.or.nypg(ncurr).eq.0.or.nzpg(ncurr).eq.0) then
    call outerror('no. of grid points in zero for potgrid',iline)
    call stopnow('genword')
  endif
  lwordok = .true.
  return
!*******************
!  Set time limit  *
!*******************
280 units = 1.0_dp
  if (nword.gt.1) then
    call stolc(words(2),maxword)
    if (index(words(2),'min').eq.1) then
      units = 60.0_dp
    elseif (index(words(2),'hou').eq.1) then
      units = 3600.0_dp
    endif
  endif
  if (nfloat.gt.0) then
    timmax = floats(1)*units
  else
    line = '  '
    read(nru,'(a)') line
    iline = iline + 1
    call linepro(nru,line,iline)
    timmax = floats(1)
    if (nword.ge.1) then
      call stolc(words(1),maxword)
      if (index(words(1),'min').eq.1) then
        units = 60.0_dp
      elseif (index(words(1),'hou').eq.1) then
        units = 3600.0_dp
      endif
    endif
    timmax = timmax*units
  endif
  lwordok = .true.
  return
!******************************
!  QEq convergance tolerance  *
!******************************
290 if (nfloat.gt.0) then
    qeqscfcrit = floats(1)
  else
    read(nru,*,err=99,end=99) qeqscfcrit
    iline = iline + 1
  endif
  qeqscfcrit = max(abs(qeqscfcrit),1.0d-8)
  lwordok = .true.
  return
!*************************************
!  QEq maximum number of iterations  *
!*************************************
300 if (nfloat.gt.0) then
    nqeqitermax = nint(floats(1))
  else
    read(nru,*,err=99,end=99) nqeqitermax
    iline = iline + 1
  endif
  nqeqitermax = max(1,abs(nqeqitermax))
  lwordok = .true.
  return
!************************
!  QEq integral radius  *
!************************
310 if (nfloat.gt.0) then
    rqeq = floats(1)
  else
    read(nru,*,err=99,end=99) rqeq
    iline = iline + 1
  endif
  rqeq = abs(rqeq)
  lwordok = .true.
  return
!******************
!  Read in title  *
!******************
320 if (nfloat.gt.0) then
    ntitle = nint(floats(1))
    if (ntitle.gt.maxtitle) then
      maxtitle = ntitle + 5
      call changemaxtitle
    endif
    do i = 1,ntitle
      read(nru,'(a)') titleword(i)
      if (ioproc) then
        write(ioout,'(''* '',a76,'' *'')') titleword(i)(1:76)
      endif
      iline = iline + 1
    enddo
  else
325 read(nru,'(a)',end=1000) line
    iline = iline + 1
!
!  Convert line to lower case before checking for "end"
!
    linelc = line
    call stolc(linelc,maxlinelength)
    if (index(linelc,'end').eq.1) goto 328
!
    ntitle = ntitle + 1
    if (ntitle.gt.maxtitle) then
      maxtitle = ntitle + 5
      call changemaxtitle
    endif
    titleword(ntitle) = line
    goto 325
  endif
328 continue
  lwordok = .true.
  return
!**************************
!  Bond length tolerance  *
!  Default = 1.2          *
!**************************
330 if (nfloat.gt.0) then
    rtol = floats(1)
  else
    read(nru,*,err=99,end=99) rtol
    iline = iline + 1
  endif
  lwordok = .true.
  return
!**********************************
!  Fortran channel for dump file  *
!  Default = 12                   *
!**********************************
340 if (nfloat.gt.0) then
    idump = int(floats(1))
  else
    read(nru,*,err=99,end=99) idump
    iline = iline + 1
  endif
  lwordok = .true.
  return
!******************************
!  Write out dumpfile at end  *
!******************************
350 idump = 12
  ntot = nfloat + nword
  if (ntot.gt.0) then
    nwo = 1
    nflo = 0
    never = -1
    j = 2
    lever = .false.
    do while (j.le.ntot)
      if (nlorder(j).eq.1) then
        nwo = nwo + 1
        if (index(words(nwo),'ever').eq.1) then
          lever = .true.
          if ((nfloat-nflo).gt.0) then
            nflo = nflo + 1
            ncycd = nint(floats(nflo))
            never = nflo
          else
            ncycd = 1
          endif
        elseif (index(words(nwo),'duri').eq.1) then
          lever = .true.
          ncycd = 1
        elseif (index(words(nwo),'cart').eq.1) then
          ldumpcart = .true.
        elseif (index(words(nwo),'noov').eq.1) then
          ldumpnooverwrite = .true.
        elseif (index(words(nwo),'conn').eq.1) then
          ldumpconnectivity = .true.
        elseif (index(words(nwo),'gene').eq.1) then
          ldumpgenerated = .true.
        elseif (index(words(nwo),'old ').eq.1) then
          ldump80 = .true.
        else
          dfile = words(nwo)
        endif
      else
        if (nflo.eq.never) then
          never = -1
        else
          nflo = nflo + 1
          idump = nint(floats(nflo))
        endif
      endif
      j = j + 1
    enddo
  endif
!
!  If not 'dump every' turn off ldumpnooverwrite option
!
  if (.not.lever.and.ncycd.ne.1) ldumpnooverwrite = .false.
!
  lwordok = .true.
  return
!*******************************************************************
!  Finite difference interval for numerical free energy gradients  *
!*******************************************************************
360 continue
  if (nword.gt.1) then
    do n = 2,nword
      if (index(words(n),'fir').eq.1) then
        lfinitediff1 = .true.
      elseif (index(words(n),'sec').eq.1) then
        lfinitediff2 = .true.
      elseif (index(words(n),'all').eq.1) then
        lfinitediff1 = .true.
        lfinitediff2 = .true.
      endif
    enddo
  else
    lfinitediff1 = .true.
    lfinitediff2 = .true.
  endif
  if (nfloat.gt.0) then
    findiff = abs(floats(1))
  else
    findiff = 1.0d-5
  endif
  lwordok = .true.
  return
!**********************
!  Core-shell cutoff  *
!**********************
370 units = 1.0_dp
  if (nword.gt.1) then
    if (index(words(2),'au').eq.1) then
      units = autoangs
    endif
  endif
  if (nfloat.gt.0) then
    cuts = abs(floats(1))*units
  else
    line = '  '
    read(nru,'(a)') line
    iline = iline + 1
    call linepro(nru,line,iline)
    cuts = abs(floats(1))
    if (nword.ge.1) then
      if (index(words(1),'au').eq.1) then
        units = autoangs
      endif
    endif
    cuts = cuts*units
  endif
  cuts = max(cuts,1.0d-8)
  lwordok = .true.
  return
!***************************
!  Minimum cell parameter  *
!***************************
380 units = 1.0_dp
  if (nword.gt.1) then
    if (index(words(2),'au').eq.1) then
      units = autoangs
    endif
  endif
  if (nfloat.gt.0) then
    cellmin = abs(floats(1))*units
  else
    line = '  '
    read(nru,'(a)') line
    iline = iline + 1
    call linepro(nru,line,iline)
    cellmin = abs(floats(1))
    if (nword.ge.1) then
      if (index(words(1),'au').eq.1) then
        units = autoangs
      endif
    endif
    cellmin = cellmin*units
  endif
  cellmin = max(cellmin,1.0d-8)
  lwordok = .true.
  return
!*********************************
!  Change default value of gtol  *
!*********************************
390 if (nword.gt.1) then
    if (index(words(2),'opt').eq.1) then
!
!  Set only optimisation gtol
!
      if (nfloat.gt.0) then
        gtol = floats(1)
      else
        read(nru,*,err=99,end=99) gtol
        iline = iline + 1
      endif
    elseif (index(words(2),'fit').eq.1) then
!
!  Set only fit gtol
!
      if (nfloat.gt.0) then
        fgtol = floats(1)
      else
        read(nru,*,err=99,end=99) fgtol
        iline = iline + 1
      endif
    else
!
!  Set according to job type
!
      if (lopt) then
        if (nfloat.gt.0) then
          gtol = floats(1)
        else
          read(nru,*,err=99,end=99) gtol
          iline = iline + 1
        endif
        if (lfit) fgtol = gtol
      elseif (lfit) then
        if (nfloat.gt.0) then
          fgtol = floats(1)
        else
          read(nru,*,err=99,end=99) fgtol
          iline = iline + 1
        endif
      endif
    endif
  else
!
!  Set according to job type
!
    if (lopt) then
      if (nfloat.gt.0) then
        gtol = floats(1)
      else
        read(nru,*,err=99,end=99) gtol
        iline = iline + 1
      endif
      if (lfit) fgtol = gtol
    elseif (lfit) then
      if (nfloat.gt.0) then
        fgtol = floats(1)
      else
        read(nru,*,err=99,end=99) fgtol
        iline = iline + 1
      endif
    endif
  endif
  gtol = abs(gtol)
  if (gtol.ge.1.0) gtol = 10.0**(-gtol)
  fgtol = abs(fgtol)
  if (fgtol.ge.1.0) fgtol = 10.0**(-fgtol)
  lwordok = .true.
  return
!*********************************
!  Change default value of gmax  *
!*********************************
400 if (nword.gt.1) then
    if (index(words(2),'opt').eq.1) then
!
!  Set only optimisation grmax
!
      if (nfloat.gt.0) then
        grmax = floats(1)
      else
        read(nru,*,err=99,end=99) grmax
        iline = iline + 1
      endif
    elseif (index(words(2),'fit').eq.1) then
!
!  Set only fit gmax
!
      if (nfloat.gt.0) then
        fgmax = floats(1)
      else
        read(nru,*,err=99,end=99) fgmax
        iline = iline + 1
      endif
    else
!
!  Set according to job type
!
      if (lopt) then
        if (nfloat.gt.0) then
          grmax = floats(1)
        else
          read(nru,*,err=99,end=99) grmax
          iline = iline + 1
        endif
        if (lfit) fgmax = grmax
      elseif (lfit) then
        if (nfloat.gt.0) then
          fgmax = floats(1)
        else
          read(nru,*,err=99,end=99) fgmax
          iline = iline + 1
        endif
      endif
    endif
  else
!
!  Set according to job type
!
    if (lopt) then
      if (nfloat.gt.0) then
        grmax = floats(1)
      else
        read(nru,*,err=99,end=99) grmax
        iline = iline + 1
      endif
      if (lfit) fgmax = grmax
    elseif (lfit) then
      if (nfloat.gt.0) then
        fgmax = floats(1)
      else
        read(nru,*,err=99,end=99) fgmax
        iline = iline + 1
      endif
    endif
  endif
  grmax = abs(grmax)
  if (grmax.ge.1.0) grmax = 10.0**(-grmax)
  fgmax = abs(fgmax)
  if (fgmax.ge.1.0) fgmax = 10.0**(-fgmax)
  lwordok = .true.
  return
!******************************************
!  Set options to make output more terse  *
!******************************************
410 if (nword.gt.1) then
    lin = .true.
    lout = .true.
    do n = 2,nword
      if (index(words(n),'in ').eq.1) then
        lin = .true.
        lout = .false.
      elseif (index(words(n),'out').eq.1) then
        lin = .false.
        lout = .true.
      elseif (index(words(n),'ino').eq.1) then
        lin = .true.
        lout = .true.
      elseif (index(words(n),'cell').eq.1) then
        if (lin.and.lout) then
          lterseincell = .true.
          lterseoutcell = .true.
        elseif (lin) then
          lterseincell = .true.
        elseif (lout) then
          lterseoutcell = .true.
        endif
      elseif (index(words(n),'coor').eq.1) then
        if (lin.and.lout) then
          lterseincoords = .true.
          lterseoutcoords = .true.
        elseif (lin) then
          lterseincoords = .true.
        elseif (lout) then
          lterseoutcoords = .true.
        endif
      elseif (index(words(n),'stru').eq.1) then
        if (lin.and.lout) then
          lterseincell = .true.
          lterseincoords = .true.
          lterseinmol = .true.
          lterseoutcell = .true.
          lterseoutcoords = .true.
        elseif (lin) then
          lterseincell = .true.
          lterseincoords = .true.
        elseif (lout) then
          lterseoutcell = .true.
          lterseoutcoords = .true.
        endif
      elseif (index(words(n),'pot').eq.1) then
        ltersepotentials = .true.
      elseif (index(words(n),'der').eq.1) then
        ltersederivs = .true.
      endif
    enddo
  endif
  lwordok = .true.
  return
!**************************************************************
!  Change default value of ftol                               *
!  Controls accuracy in energy during optimisation with BFGS  *
!**************************************************************
420 if (nword.gt.1) then
    if (index(words(2),'opt').eq.1) then
!
!  Set only optimisation ftol
!
      if (nfloat.gt.0) then
        ftol = floats(1)
      else
        read(nru,*,err=99,end=99) ftol
        iline = iline + 1
      endif
    elseif (index(words(2),'fit').eq.1) then
!
!  Set only fit ftol
!
      if (nfloat.gt.0) then
        fftol = floats(1)
      else
        read(nru,*,err=99,end=99) fftol
        iline = iline + 1
      endif
    else
!
!  Set according to job type
!
      if (lopt) then
        if (nfloat.gt.0) then
          ftol = floats(1)
        else
          read(nru,*,err=99,end=99) ftol
          iline = iline + 1
        endif
        if (lfit) fftol = ftol
      elseif (lfit) then
        if (nfloat.gt.0) then
          fftol = floats(1)
        else
          read(nru,*,err=99,end=99) fftol
          iline = iline + 1
        endif
      endif
    endif
  else
!
!  Set according to job type
!
    if (lopt) then
      if (nfloat.gt.0) then
        ftol = floats(1)
      else
        read(nru,*,err=99,end=99) ftol
        iline = iline + 1
      endif
      if (lfit) fftol = ftol
    elseif (lfit) then
      if (nfloat.gt.0) then
        fftol = floats(1)
      else
        read(nru,*,err=99,end=99) fftol
        iline = iline + 1
      endif
    endif
  endif
  ftol = abs(ftol)
  if (ftol.ge.1.0) ftol = 10.0**(-ftol)
  fftol = abs(fftol)
  if (fftol.ge.1.0) fftol = 10.0**(-fftol)
  lwordok = .true.
  return
!************************************************************
!  Output current parameters every ncycp cycles of fitting  *
!************************************************************
430 if (nfloat.gt.0) then
    ncycp = int(floats(1))
  else
    read(nru,*,err=99,end=99) ncycp
    iline = iline + 1
  endif
  lwordok = .true.
  return
!************************
!  Change element data  *
!************************
440 line = '  '
  read(nru,'(a)',err=99,end=99) line
  iline = iline + 1
  if (index(line,'#').eq.1) goto 440
  call linepro(nru,line,iline)
445 word = words(1)
  if (index(word,'mass').eq.1) then
    latmassin = .true.
    call getpotsymbol1(iline,llibrary,inat,itype,sym1,1_i4,nbeg,lvalid)
    if (lvalid) then
      if (itype.eq.0) then
!
!  Standard atomic mass
!
        if (nfloat.ge.nbeg+1) then
          val = floats(nbeg+1)
        else
          read(nru,*,err=99,end=99) val
          iline = iline + 1
        endif
      else
!
!  Species specific mass
!
        latmassin = .false.
!
!  Find species
!
        lfound = .false.
        ns = 0
        do while (ns.lt.nspec.and..not.lfound)
          ns = ns + 1
          if (inat.eq.natspec(ns).and.itype.eq.ntypspec(ns)) then
            lfound = .true.
          endif
        enddo
        if (lfound) then
          massspec(ns) = floats(1)
          lmassinspec(ns) = .true.
        else
          nspec = nspec + 1
          if (nspec.gt.maxspec) then
            maxspec = nspec + 20 
            call changemaxspec
          endif
          linspec(nspec) = .true.
          natspec(nspec) = inat
          ntypspec(nspec) = itype
          massspec(nspec) = floats(1)
          lmassinspec(nspec) = .true.
          vdwspec(nspec) = rvdw(inat)
          lvdwinspec(nspec) = .false.
          qlspec(nspec) = 0.0_dp
          radspec(nspec) = 0.0_dp
        endif
      endif
    endif
    if (latmassin.and.lvalid) then
      atmass(inat) = val
!
!  Update species masses where not set already based on specific species
!
      do i = 1,nspec
        if (.not.lmassinspec(i)) then
          if (natspec(i).eq.inat) then
            if (natspec(i).le.maxele) then
              massspec(i) = atmass(natspec(i))
            else
              massspec(i) = 0.0_dp
            endif
          endif
        endif
      enddo
    endif
    goto 440
  elseif (index(word,'nmr').eq.1) then
    call getpotsymbol1(iline,llibrary,inat,itype,sym1,1_i4,nbeg,lvalid)
    if (lvalid) then
!
!  Species specific NMR parameters
!
!  Find species
!
      lfound = .false.
      ns = 0
      do while (ns.lt.nspec.and..not.lfound)
        ns = ns + 1
        if (inat.eq.natspec(ns).and.itype.eq.ntypspec(ns)) then
          lfound = .true.
        endif
      enddo
      if (nfloat.lt.3) then
        call outerror('insufficient input reading nmr',iline)
        call stopnow('genword')
      endif
      if (lfound) then
        nmrspec(1,ns) = floats(1)
        nmrspec(2,ns) = floats(2)
        nmrspec(3,ns) = floats(3)
        lnmrinspec(ns) = .true.
      else
        nspec = nspec + 1
        if (nspec.gt.maxspec) then
          maxspec = nspec + 20
          call changemaxspec
        endif
        linspec(nspec) = .true.
        natspec(nspec) = inat
        ntypspec(nspec) = itype
        massspec(nspec) = atmass(inat)
        lmassinspec(nspec) = .false.
        vdwspec(nspec) = rvdw(inat)
        lvdwinspec(nspec) = .false.
        nmrspec(1,nspec) = floats(1)
        nmrspec(2,nspec) = floats(2)
        nmrspec(3,nspec) = floats(3)
        lnmrinspec(nspec) = .true.
        ns = nspec
      endif
    endif
    goto 440
  elseif (index(word,'ioni').eq.1) then
    if (nword.ge.2) then
      s1 = words(2)(1:2)
      na = istoan(s1)
      if (nfloat.ge.1) then
        val = floats(1)
      else
        read(nru,*,err=99,end=99) val
        iline = iline + 1
      endif
    elseif (nfloat.ge.1) then
      na = nint(floats(1))
      if (nfloat.ge.2) then
        val = floats(2)
      else
        read(nru,*,err=99,end=99) val
        iline = iline + 1
      endif
    else
      line = '  '
      read(nru,'(a)',err=99,end=99) line
      iline = iline + 1
      call linepro(nru,line,iline)
      if (nword.gt.0.and.nfloat.gt.0) then
        s1 = words(1)(1:2)
        na = istoan(s1)
        val = floats(1)
      elseif (nfloat.ge.2) then
        na = nint(floats(1))
        val = floats(2)
      else
        call outerror('insufficient input reading mass',iline)
        call stopnow('genword')
      endif
    endif
    rion(na) = val
    goto 440
  elseif (index(word,'cova').eq.1) then
    if (nword.ge.2) then
      s1 = words(2)(1:2)
      na = istoan(s1)
      if (nfloat.ge.1) then
        val = floats(1)
      else
        read(nru,*,err=99,end=99) val
        iline = iline + 1
      endif
    elseif (nfloat.ge.1) then
      na = nint(floats(1))
      if (nfloat.ge.2) then
        val = floats(2)
      else
        read(nru,*,err=99,end=99) val
        iline = iline + 1
      endif
    else
      line = '  '
      read(nru,'(a)',err=99,end=99) line
      iline = iline + 1
      call linepro(nru,line,iline)
      if (nword.gt.0.and.nfloat.gt.0) then
        s1 = words(1)(1:2)
        na = istoan(s1)
        val = floats(1)
      elseif (nfloat.ge.2) then
        na = nint(floats(1))
        val = floats(2)
      else
        call outerror('insufficient input reading covalent radius',iline)
        call stopnow('genword')
      endif
    endif
    rcov(na) = val
    goto 440
  elseif (index(word,'spin').eq.1) then
    call getpotsymbol1(iline,llibrary,inat,itype,sym1,1_i4,nbeg,lvalid)
    if (lvalid) then
!
!  Species specific spin parameters
!
!  Find species
!
      lfound = .false.
      ns = 0
      do while (ns.lt.nspec.and..not.lfound)
        ns = ns + 1
        if (inat.eq.natspec(ns).and.itype.eq.ntypspec(ns)) then
          lfound = .true.
        endif
      enddo
      if (nfloat.lt.1) then
        call outerror('insufficient input reading spin',iline)
        call stopnow('genword')
      endif
      if (lfound) then
        spinspec(ns) = floats(1)
        lspininspec(ns) = .true.
      else
        nspec = nspec + 1
        if (nspec.gt.maxspec) then
          maxspec = nspec + 20
          call changemaxspec
        endif
        linspec(nspec) = .true.
        natspec(nspec) = inat
        ntypspec(nspec) = itype
        massspec(nspec) = atmass(inat)
        lmassinspec(nspec) = .false.
        vdwspec(nspec) = rvdw(inat)
        lvdwinspec(nspec) = .false.
        spinspec(nspec) = floats(1)
        lspininspec(nspec) = .true.
        ns = nspec
      endif
    endif
    goto 440
  elseif (index(word,'vdw').eq.1) then
    lvdwin = .true.
    call getpotsymbol1(iline,llibrary,inat,itype,sym1,1_i4,nbeg,lvalid)
    if (lvalid) then
      if (itype.eq.0) then
!
!  Standard atomic VDW radius
!
        if (nfloat.ge.nbeg+1) then
          val = floats(nbeg+1)
        else
          read(nru,*,err=99,end=99) val
          iline = iline + 1
        endif
      else
!
!  Species specific VDW radius
!
        lvdwin = .false.
!
!  Find species
!
        lfound = .false.
        ns = 0
        do while (ns.lt.nspec.and..not.lfound)
          ns = ns + 1
          if (inat.eq.natspec(ns).and.itype.eq.ntypspec(ns)) then
            lfound = .true.
          endif
        enddo
        if (lfound) then
          vdwspec(ns) = floats(1)
          lvdwinspec(ns) = .true.
        else
          nspec = nspec + 1
          if (nspec.gt.maxspec) then
            maxspec = nspec + 20 
            call changemaxspec
          endif
          linspec(nspec) = .true.
          natspec(nspec) = inat
          ntypspec(nspec) = itype
          massspec(nspec) = atmass(inat)
          lmassinspec(nspec) = .false.
          vdwspec(nspec) = floats(1)
          lvdwinspec(nspec) = .true.
          qlspec(nspec) = 0.0_dp
          radspec(nspec) = 0.0_dp
        endif
      endif
    endif
    if (lvdwin.and.lvalid) then
      rvdw(inat) = val
!
!  Update species VDW radii where not set already based on specific species
!
      do i = 1,nspec
        if (.not.lvdwinspec(i)) then
          if (natspec(i).eq.inat) then
            if (natspec(i).le.maxele) then
              vdwspec(i) = rvdw(natspec(i))
            else
              vdwspec(i) = 0.0_dp
            endif
          endif
        endif
      enddo
    endif
    goto 440
  elseif (index(word,'symb').eq.1) then
    if (nfloat.ge.1) then
      na = nint(floats(1))
      if (nword.ge.2) then
        cha = words(2)(1:2)
      else
        line = '  '
        read(nru,'(a)',err=99,end=99) line
        iline = iline + 1
        call linepro(nru,line,iline)
        cha = words(1)(1:2)
      endif
    else
      line = '  '
      read(nru,'(a)',err=99,end=99) line
      iline = iline + 1
      call linepro(nru,line,iline)
      if (nfloat.gt.0.and.nword.gt.0) then
        na = nint(floats(1))
        cha = words(1)(1:2)
      else
        call outerror('insufficient input reading element symbol',iline)
        call stopnow('genword')
      endif
    endif
    atsym(na) = cha
    goto 440
  elseif (index(word,'bbar').eq.1) then
    if (nword.ge.2) then
      s1 = words(2)(1:2)
      na = istoan(s1)
      if (nfloat.ge.1) then
        val = floats(1)
      else
        read(nru,*,err=99,end=99) val
        iline = iline + 1
      endif
    elseif (nfloat.ge.1) then
      na = nint(floats(1))
      if (nfloat.ge.2) then
        val = floats(2)
      else
        read(nru,*,err=99,end=99) val
        iline = iline + 1
      endif
    else
      line = '  '
      read(nru,'(a)',err=99,end=99) line
      iline = iline + 1
      call linepro(nru,line,iline)
      if (nword.gt.0.and.nfloat.gt.0) then
        s1 = words(1)(1:2)
        na = istoan(s1)
        val = floats(1)
      elseif (nfloat.ge.2) then
        na = nint(floats(1))
        val = floats(2)
      else
        call outerror('insufficient input reading bbar',iline)
        call stopnow('genword')
      endif
    endif
    bbar(na) = val !in Angstrom
    goto 440
  elseif (index(word,'sigi').eq.1) then
    if (nword.ge.2) then
      s1 = words(2)(1:2)
      na = istoan(s1)
      if (nfloat.ge.1) then
        val = floats(1)
      else
        read(nru,*,err=99,end=99) val
        iline = iline + 1
      endif
    elseif (nfloat.ge.1) then
      na = nint(floats(1))
      if (nfloat.ge.2) then
        val = floats(2)
      else
        read(nru,*,err=99,end=99) val
        iline = iline + 1
      endif
    else
      line = '  '
      read(nru,'(a)',err=99,end=99) line
      iline = iline + 1
      call linepro(nru,line,iline)
      if (nword.gt.0.and.nfloat.gt.0) then
        s1 = words(1)(1:2)
        na = istoan(s1)
        val = floats(1)
      elseif (nfloat.ge.2) then
        na = nint(floats(1))
        val = floats(2)
      else
        call outerror('insufficient input reading sigi',iline)
        call stopnow('genword')
      endif
    endif
    siginc(na) = val*1.0d8 !input units Angstrom squared, stored in barns
    goto 440
  elseif (index(word,'end').eq.1) then
    lwordok = .true.
    return
  else
    l55 = .true.
    return
  endif
!*****************************************
!  S and M electronegativity parameters  *
!*****************************************
490 units = 1.0_dp
  nqrtyp = 0
  if (nword.gt.1) then
    do i = 2,nword
      call stolc(words(i),maxword)
      if (index(words(i),'au').eq.1) units = autoangs
      if (index(words(i),'qmi').eq.1) nqrtyp = 1
      if (index(words(i),'qma').eq.1) nqrtyp = 2
      if (index(words(i),'qra').eq.1) nqrtyp = 3
    enddo
  endif
!
!  Start of input loop
!
495 line = '  '
  read(nru,'(a)',end=1000) line
  iline = iline + 1
  call linepro(nru,line,iline)
!
!  Check whether input is symbol or option
!
  if ((nword+nfloat).eq.0) goto 495
!
!  Check for old fashion specification of number of atoms
!
  if (nword.gt.0) then
    word = words(1)
    call worsy(word,lsymbol,.false.)
    if (.not.lsymbol) then
      l55 = .true.
      return
    endif
  endif
!
!  Find atom number for element
!
!  Symbols used in input
!
  if (nword.gt.0) then
    call ltont(words(1),inat,itype)
    leemparaltered(inat) = .true.
    nfloatused = 0
  else
    inat = int(floats(1))
    leemparaltered(inat) = .true.
    nfloatused = 1
  endif
!
!  Set end of float range for given qrange type
!
  nfloatend = nfloat
  if (nqrtyp.eq.3) then
    if (lfit.and.lfflags) then
      nfloatend = min(nfloatend,nfloatused+12_i4)
      if (nfloatend.lt.10_i4) then
        call outerror('insufficient input parameters for electronegativity',iline)
        call stopnow('genword')
      endif
    else
      nfloatend = min(nfloatend,nfloatused+8_i4)
      if (nfloatend.lt.6_i4) then
        call outerror('insufficient input parameters for electronegativity',iline)
        call stopnow('genword')
      endif
    endif
  elseif (nqrtyp.eq.2.or.nqrtyp.eq.1) then
    if (lfit.and.lfflags) then
      nfloatend = min(nfloatend,nfloatused+11_i4)
      if (nfloatend.lt.9_i4) then
        call outerror('insufficient input parameters for electronegativity',iline)
        call stopnow('genword')
      endif
    else
      nfloatend = min(nfloatend,nfloatused+7_i4)
      if (nfloatend.lt.5_i4) then
        call outerror('insufficient input parameters for electronegativity',iline)
        call stopnow('genword')
      endif
    endif
  else
    if (lfit.and.lfflags) then
      nfloatend = min(nfloatend,nfloatused+10_i4)
      if (nfloatend.lt.8_i4) then
        call outerror('insufficient input parameters for electronegativity',iline)
        call stopnow('genword')
      endif
    else
      nfloatend = min(nfloatend,nfloatused+6_i4)
      if (nfloatend.lt.4_i4) then
        call outerror('insufficient input parameters for electronegativity',iline)
        call stopnow('genword')
      endif
    endif
  endif
!
!  Increment number of ranges
!
  nqrange(inat,3) = nqrange(inat,3) + 1
!
!  Check dimensions
!
  if (nqrange(inat,3).gt.maxqrange) then
    maxqrange = nqrange(inat,3)
    call changemaxqrange
  endif
!
!  Assign range type 
!
  nqrangetype(nqrange(inat,3),inat,3) = nqrtyp
!    
!  Fitting flags
!
  if (lfit.and.lfflags) then
    n1 = int(floats(nfloatend-3))
    n2 = int(floats(nfloatend-2))
    n3 = int(floats(nfloatend-1))
    n4 = int(floats(nfloatend))
    if (n1.eq.1) then
      nfit = nfit + 1
      if (nfit.ge.maxfit) then
        maxfit = nfit + 10
        call changemaxfit
      endif
      nftyp(nfit) = 1
      nfpot(nfit) = 11
      nfvar(nfit) = inat
      nfvar2(nfit) = nqrange(inat,3)
      nfvar3(nfit) = 3
    endif
    if (n2.eq.1) then
      nfit = nfit + 1
      if (nfit.ge.maxfit) then
        maxfit = nfit + 10
        call changemaxfit
      endif
      nftyp(nfit) = 1
      nfpot(nfit) = 12
      nfvar(nfit) = inat
      nfvar2(nfit) = nqrange(inat,3)
      nfvar3(nfit) = 3
    endif
    if (n3.eq.1) then
      nfit = nfit + 1
      if (nfit.ge.maxfit) then
        maxfit = nfit + 10
        call changemaxfit
      endif
      nftyp(nfit) = 1
      nfpot(nfit) = 14
      nfvar(nfit) = inat
      nfvar2(nfit) = nqrange(inat,3)
      nfvar3(nfit) = 3
    endif
    if (n4.eq.1) then
      nfit = nfit + 1
      if (nfit.ge.maxfit) then
        maxfit = nfit + 10
        call changemaxfit
      endif
      nftyp(nfit) = 1
      nfpot(nfit) = 15
      nfvar(nfit) = inat
      nfvar2(nfit) = nqrange(inat,3)
      nfvar3(nfit) = 3
    endif
    nfloatend = nfloatend - 4
  endif
!
!  Assign range parameters
!
  if (nqrtyp.eq.3) then
    qrangemax(nqrange(inat,3),inat,3) = floats(nfloatend)
    qrangemin(nqrange(inat,3),inat,3) = floats(nfloatend-1)
    nfloatend = nfloatend - 2
  elseif (nqrtyp.eq.2) then
    qrangemax(nqrange(inat,3),inat,3) = floats(nfloatend)
    nfloatend = nfloatend - 1
  elseif (nqrtyp.eq.1) then
    qrangemin(nqrange(inat,3),inat,3) = floats(nfloatend)
    nfloatend = nfloatend - 1
  endif
!
!  Assign electronegativity parameters
!
  nfloatused = nfloatused + 1
  chirange(nqrange(inat,3),inat,3) = floats(nfloatused)
  if (nfloatused.lt.nfloatend) then
    nfloatused = nfloatused + 1
    murange(nqrange(inat,3),inat,3) = floats(nfloatused)
    if (nfloatused.lt.nfloatend) then
      nfloatused = nfloatused + 1
      zetarange(nqrange(inat,3),inat,3) = floats(nfloatused)
      if (nfloatused.lt.nfloatend) then
        nfloatused = nfloatused + 1
        znucrange(nqrange(inat,3),inat,3) = floats(nfloatused)
        if (nfloatused.lt.nfloatend) then
          nfloatused = nfloatused + 1
          q0range(nqrange(inat,3),inat,3) = floats(nfloatused)
          if (nfloatused.lt.nfloatend) then
            nfloatused = nfloatused + 1
            e0range(nqrange(inat,3),inat,3) = floats(nfloatused)
          endif
        endif
      endif
    endif
  endif
  goto 495
!*****************************************
!  Maximum number of bfgs line searches  *
!*****************************************
500 if (nfloat.gt.0) then
    maxline = int(floats(1))
  else
    read(nru,*,err=99,end=99) maxline
    iline = iline + 1
  endif
  lwordok = .true.
  return
!************************************************************************
!  Select additional output file formats - option to supply file names  *
!************************************************************************
510 if (nword.lt.2) then
    if (ioproc) then
      write(ioout,'(''  **** No output file type specified ****'')')
      write(ioout,'(''  **** Directive ignored             ****'')')
    endif
    lwordok = .true.
    return
  endif
  iwc = 0
  llocalseq = .false.
  llocalsym = .false.
  llocalmovie = .false.
  nextraword = 0
  do i = 2,nword
    word = words(i)
    lset = .false.
    lnodot = (index(word,'.').eq.0)
!
!  Check that there is no dot in the word, 
!  otherwise it's a file name not an option
!
    if (lnodot) then
      if (.not.lmarv) then
        if (index(word,'marv').eq.1) then
          lmarv = .true.
          lmarv2 = (index(word,'2').ne.0)
          iwc = 1
          lset = .true.
        endif
      endif
      if (.not.lthb) then
        if (index(word,'thb ').eq.1) then
          lthb = .true.
          iwc = 2
          lset = .true.
        endif
      endif
      if (.not.lxtl) then
        if (index(word,'xtl ').eq.1) then
          lxtl = .true.
          iwc = 3
          lset = .true.
        endif
      endif
      if (.not.lphono) then
        if (index(word,'phon').eq.1) then
          lphono = .true.
          iwc = 4
          lset = .true.
        endif
      endif
      if (.not.larc) then
        if (index(word,'arc ').eq.1) then
          larc = .true.
          iwc = 5
          lset = .true.
        endif
      endif
      if (.not.lxr) then
        if (index(word,'xr ').eq.1) then
          lxr = .true.
          iwc = 6
          lset = .true.
        endif
      endif
      if (.not.lcssr) then
        if (index(word,'cssr').eq.1) then
          lcssr = .true.
          iwc = 7
          lset = .true.
        endif
      endif
      if (.not.ltrj) then
        if (index(word,'traj').eq.1) then
          ltrj = .true.
          iwc = 8
          lset = .true.
        endif
      else
        if (index(word,'asci').eq.1) then
          ltrjascii  =  .true.
        endif
        if (index(word,'equi').eq.1) then
          ltrjequil  =  .true.
        endif
      endif
      if (.not.lfrq) then
        if (index(word,'freq').eq.1) then
          lfrq = .true.
          iwc = 9
          lset = .true.
        endif
      endif
      if (.not.lxyz) then
        if (index(word,'xyz ').eq.1) then
          lxyz = .true.
          iwc = 10
          lset = .true.
        endif
      endif
      if (.not.lhis) then
        if (index(word,'his').eq.1) then
          lhis = .true.
          iwc = 11
          lset = .true.
        endif
      endif
      if (.not.lfdf) then
        if (index(word,'fdf ').eq.1) then
          lfdf = .true.
          iwc = 12
          lset = .true.
        endif
      endif
      if (.not.ldrv) then
        if (index(word,'der').eq.1.or.index(word,'drv ').eq.1) then
          ldrv = .true.
          iwc = 13
          lset = .true.
        endif
      endif
      if (.not.lfrc) then
        if (index(word,'frc ').eq.1) then
          lfrc = .true.
          iwc = 14
          lset = .true.
        endif
      endif
      if (.not.lcif) then
        if (index(word,'cif ').eq.1) then
          lcif = .true.
          iwc = 15
          lset = .true.
        endif
      endif
      if (.not.ldlv) then
        if (index(word,'str ').eq.1) then
          ldlv = .true.
          iwc = 16
          lset = .true.
        endif
      endif
      if (.not.leig) then
        if (index(word,'eig ').eq.1) then
          leig = .true.
          iwc = 17
          lset = .true.
        endif
      endif
      if (.not.lpre) then
        if (index(word,'pres').eq.1) then
          if (ioproc) lpre = .true.
          iwc = 19
          lset = .true.
        endif
      endif
      if (.not.lsas) then
        if (index(word,'sas ').eq.1) then
          if (ioproc) lsas = .true.
          iwc = 20
          lset = .true.
        endif
      endif
      if (.not.losc) then
        if (index(word,'osc ').eq.1) then
          if (ioproc) losc = .true.
          iwc = 21
          lset = .true.
        endif
      endif
      if (.not.lbio) then
        if (index(word,'bio ').eq.1) then
          if (ioproc) lbio = .true.
          iwc = 22
          lset = .true.
        endif
      endif
      if (index(word,'pdf').eq.1) then
        iwc = 23
        lset = .true.
      endif
      if (.not.lcosmofile) then
        if (index(word,'cos').eq.1) then
          lcosmofile = .true.
          iwc = 24
          lset = .true.
        endif
      endif
      if (.not.lqbo) then
        if (index(word,'qbo ').eq.1) then
          if (ioproc) lqbo = .true.
          iwc = 25
          lset = .true.
        endif
      endif
      if (.not.llammpspots) then
        if (index(word,'lammps_p').eq.1) then
          llammpspots = .true.
          iwc = 26
          lset = .true.
        endif
      endif
      if (.not.ldcd) then
        if (index(word,'dcd').eq.1) then
          ldcd = .true.
          iwc = 27
          lset = .true.
        endif
      endif
      if (.not.llammps) then
        if (index(word,'lam').eq.1) then
          llammps = .true.
          iwc = 28
          lset = .true.
        endif
      endif
      if (.not.lcas) then
        if (index(word,'cas').eq.1) then
          lcas = .true.
          iwc = 29
          lset = .true.
        endif
      endif
      if (.not.lshengBTE) then
        if (index(word,'she').eq.1) then
          lshengBTE = .true.
          iwc = 30
          lset = .true.
        endif
      endif
      if (.not.linertia) then
        if (index(word,'ine').eq.1) then
          linertia = .true.
          iwc = 31
          lset = .true.
        endif
      endif
      if (.not.lkpt) then
        if (index(word,'kpt').eq.1) then
          lkpt = .true.
          iwc = 32
          lset = .true.
        endif
      endif
      if (.not.ldip) then
        if (index(word,'dip').eq.1) then
          ldip = .true.
          iwc = 33
          lset = .true.
        endif
      endif
      if (.not.lxsf) then
        if (index(word,'xsf ').eq.1) then
          lxsf = .true.
          iwc = 34
          lset = .true.
        endif
      endif
!
!  CML output
!
      if (.not.lcml) then
        if (index(word,'cml ').eq.1) then
          if (ioproc) lcml = .true.
          iwc = 1000 
          lset = .true.
        endif 
      endif
      if (.not.lvcml) then 
        if (index(word,'vcml ').eq.1) then
          if (ioproc) lcml = .true.
          if (ioproc) lvcml = .true.
          iwc = 1001
          lset = .true.
        endif
      endif
    endif
    if (index(word,'movi').eq.1.and.lnodot) then
      llocalmovie = .true.
      nextraword = nextraword + 1
    elseif (index(word,'sym').eq.1.and.lnodot) then
      llocalsym = .true.
      nextraword = nextraword + 1
    elseif (index(word,'seq').eq.1.and.lnodot) then
      llocalseq = .true.
      nextraword = nextraword + 1
    elseif (index(word,'shel').eq.1.and.lnodot) then
      loutshell = .true.
      nextraword = nextraword + 1
    elseif (lfrq.and.index(word,'bin').eq.1.and.lnodot) then
      lfrqbin = .true.
      nextraword = nextraword + 1
    elseif (index(word,'tex').eq.1.and.lnodot) then
      lfrqbin = .false.
      nextraword = nextraword + 1
    elseif (llammps.and.index(word,'var').eq.1.and.lnodot) then
      llammpsvar = .true.
      nextraword = nextraword + 1
    elseif (lqbo.and.index(word,'app').eq.1.and.lnodot) then
      lqboappend = .true.
      nextraword = nextraword + 1
    elseif (lphono.and.index(word,'ir').eq.1) then
      lirout = .true.
      linten = .true.
      nextraword = nextraword + 1
    elseif (.not.lset) then
!
!  Assign any output files that may have been given
!
      if (iwc.eq.1) then
        marvfile = words(i)
        if (index(marvfile,'.mvn').eq.0) then
          call endstring(marvfile,len(marvfile),iend)
          if (iend.eq.0) iend = len(marvfile) - 3
          marvfile(iend:iend+3) = '.mvn'
        endif
      elseif (iwc.eq.2) then
        thbfile = words(i)
      elseif (iwc.eq.3) then
        xtlfile = words(i)
        if (index(xtlfile,'.xtl').eq.0) then
          call endstring(xtlfile,len(xtlfile),iend)
          if (iend.eq.0) iend = len(xtlfile) - 3
          xtlfile(iend:iend+3) = '.xtl'
        endif
      elseif (iwc.eq.4) then
        phonfile = words(i)
      elseif (iwc.eq.5) then
        arcfile = words(i)
        if (index(arcfile,'.arc').eq.0.and.index(arcfile,'.car').eq.0) then
          call endstring(arcfile,len(arcfile),iend)
          if (iend.eq.0) iend = len(arcfile) - 3
          arcfile(iend:iend+3) = '.arc'
        endif
        lmovie = llocalmovie
        if (nfloat.gt.0.and.lmovie) then
          narcwrite = nint(abs(floats(1)))
          narcwrite = max(narcwrite,1_i4)
        endif
      elseif (iwc.eq.6) then
        xrfile = words(i)
        if (index(xrfile,'.xr').eq.0) then
          call endstring(xrfile,len(xrfile),iend)
          if (iend.eq.0) iend = len(xrfile) - 2
          xrfile(iend:iend+2) = '.xr'
        endif
      elseif (iwc.eq.7) then
        cssrfile = words(i)
        if (index(cssrfile,'.cssr').eq.0) then
          call endstring(cssrfile,len(cssrfile),iend)
          if (iend.eq.0) iend = len(cssrfile) - 4
          cssrfile(iend:iend+4) = '.cssr'
        endif
      elseif (iwc.eq.8) then
        trjfile = words(i)
        if (index(trjfile,'.trg').eq.0) then
          call endstring(trjfile,len(trjfile),iend)
          if (iend.eq.0) iend = len(trjfile) - 3
          trjfile(iend:iend+3) = '.trg'
        endif
      elseif (iwc.eq.9) then
        freqfile = words(i)
        if (nfloat.gt.0) then
          nfreqdecimals = abs(nint(floats(1)))
          nfreqdecimals = min(nfreqdecimals,12_i4)
          nfreqdecimals = max(nfreqdecimals,1_i4)
        endif
      elseif (iwc.eq.10) then
        xyzfile = words(i)
        if (index(xyzfile,'.xyz').eq.0) then
          call endstring(xyzfile,len(xyzfile),iend)
          if (iend.eq.0) iend = len(xyzfile) - 3
          xyzfile(iend:iend+3) = '.xyz'
        endif
        lxyzmovie = llocalmovie
      elseif (iwc.eq.11) then
        hisfile = words(i)
        if (index(hisfile,'.his').eq.0) then
          call endstring(hisfile,len(hisfile),iend)
          if (iend.eq.0) iend = len(hisfile) - 3
          hisfile(iend:iend+3) = '.his'
        endif
      elseif (iwc.eq.12) then
        fdffile = words(i)
        if (index(fdffile,'.fdf').eq.0) then
          call endstring(fdffile,len(fdffile),iend)
          if (iend.eq.0) iend = len(fdffile) - 3
          fdffile(iend:iend+3) = '.fdf'
        endif
      elseif (iwc.eq.13) then
        drvfile = words(i)
        if (index(drvfile,'.drv').eq.0) then
          call endstring(drvfile,len(drvfile),iend)
          if (iend.eq.0) iend = len(drvfile) - 3
          drvfile(iend:iend+3) = '.drv'
        endif
      elseif (iwc.eq.14) then
        frcfile = words(i)
        if (index(frcfile,'.frc').eq.0) then
          call endstring(frcfile,len(frcfile),iend)
          if (iend.eq.0) iend = len(frcfile) - 3
          frcfile(iend:iend+3) = '.frc'
        endif
      elseif (iwc.eq.15) then
        ciffile = words(i)
        if (index(ciffile,'.cif').eq.0) then
          call endstring(ciffile,len(ciffile),iend)
          if (iend.eq.0) iend = len(ciffile) - 3
          ciffile(iend:iend+3) = '.cif'
        endif
        if (nfloat.gt.0) then
          cif_dummylattice = abs(floats(1))
        endif
      elseif (iwc.eq.16) then
        dlvfile = words(i)
        if (index(dlvfile,'.str').eq.0) then
          call endstring(dlvfile,len(dlvfile),iend)
          if (iend.eq.0) iend = len(dlvfile) - 3
          dlvfile(iend:iend+3) = '.str'
        endif
      elseif (iwc.eq.17) then
        eigfile = words(i)
        if (index(eigfile,'.eig').eq.0) then
          call endstring(eigfile,len(eigfile),iend)
          if (iend.eq.0) iend = len(eigfile) - 3
          eigfile(iend:iend+3) = '.eig'
        endif
      elseif (iwc.eq.19) then
        prefile = words(i)
        if (index(prefile,'.pre').eq.0) then
          call endstring(prefile,len(prefile),iend)
          if (iend.eq.0) iend = len(prefile) - 3
          prefile(iend:iend+3) = '.pre'
        endif
      elseif (iwc.eq.20) then
        sasfile = words(i)
        if (index(sasfile,'.sas').eq.0) then
          call endstring(sasfile,len(sasfile),iend)
          if (iend.eq.0) iend = len(sasfile) - 3
          sasfile(iend:iend+3) = '.sas'
        endif
      elseif (iwc.eq.21) then
        oscfile = words(i)
        if (index(oscfile,'.osc').eq.0) then
          call endstring(oscfile,len(oscfile),iend)
          if (iend.eq.0) iend = len(oscfile) - 3
          oscfile(iend:iend+3) = '.osc'
        endif
      elseif (iwc.eq.22) then
        biofile = words(i)
        if (index(biofile,'.bio').eq.0) then
          call endstring(biofile,len(biofile),iend)
          if (iend.eq.0) iend = len(biofile) - 3
          biofile(iend:iend+3) = '.bio'
        endif
      elseif (iwc.eq.23) then
        pdffiles(ncurr) = words(i)
        if (index(pdffiles(ncurr),'.wid').eq.0) then
          iend = index(pdffiles(ncurr),' ')
          call endstring(pdffiles(ncurr),len(pdffiles(ncurr)),iend)
          if (iend.eq.0) iend = len(pdffiles(ncurr)) - 3
          pdffiles(ncurr)(iend:iend+3) = '.wid'
        endif
      elseif (iwc.eq.24) then
        cosmofile = words(i)
        if (index(cosmofile,'.cosmo').eq.0) then
          call endstring(cosmofile,len(cosmofile),iend)
          if (iend.eq.0) iend = len(cosmofile) - 5
          cosmofile(iend:iend+5) = '.cosmo'
        endif
      elseif (iwc.eq.25) then
        qbofile = words(i)
        if (index(qbofile,'.qbo').eq.0) then
          call endstring(qbofile,len(qbofile),iend)
          if (iend.eq.0) iend = len(qbofile) - 3
          qbofile(iend:iend+3) = '.qbo'
        endif
      elseif (iwc.eq.26) then
        lammpspotsfile = words(i)
        if (index(lammpspotsfile,'.tab').eq.0) then
          call endstring(lammpspotsfile,len(lammpspotsfile),iend)
          if (iend.eq.0) iend = len(lammpspotsfile) - 3
          lammpspotsfile(iend:iend+3) = '.tab'
        endif
        if (nfloat.ge.3) then
          lammps_r0   = abs(floats(1))
          lammps_rend = abs(floats(2))
          if (lammps_rend.lt.lammps_r0) then
            call outerror('lammps potentials - end point before start',iline)
            call stopnow('genword')
          endif
          nlammpspoints   = nint(abs(floats(3)))
          if (nlammpspoints.le..1) then
            call outerror('lammps potentials - number of points too small',iline)
            call stopnow('genword')
          endif
        else
          call outerror('value(s) missing from input for lammps potentials',iline)
          call stopnow('genword')
        endif
      elseif (iwc.eq.27) then
        dcdfile = words(i)
        if (index(dcdfile,'.dcd').eq.0) then
          call endstring(dcdfile,len(dcdfile),iend)
          if (iend.eq.0) iend = len(dcdfile) - 3
          dcdfile(iend:iend+3) = '.dcd'
        endif
      elseif (iwc.eq.28) then
        lammpsfile = words(i)
        if (index(lammpsfile,'.lmp').eq.0) then
          call endstring(lammpsfile,len(lammpsfile),iend)
          if (iend.eq.0) iend = len(lammpsfile) - 3
          lammpsfile(iend:iend+3) = '.lmp'
        endif
      elseif (iwc.eq.29) then
        casfile = words(i)
        if (index(casfile,'.csp').eq.0) then
          call endstring(casfile,len(casfile),iend)
          if (iend.eq.0) iend = len(casfile) - 3
          casfile(iend:iend+3) = '.csp'
        endif
      elseif (iwc.eq.30) then
        shengfile = words(i)
      elseif (iwc.eq.31) then
        inertfile = words(i)
        if (index(inertfile,'.ine').eq.0) then
          call endstring(inertfile,len(inertfile),iend)
          if (iend.eq.0) iend = len(inertfile) - 3
          inertfile(iend:iend+3) = '.ine'
        endif
      elseif (iwc.eq.32) then
        kptfile = words(i)
        if (index(kptfile,'.kpt').eq.0) then
          call endstring(kptfile,len(kptfile),iend)
          if (iend.eq.0) iend = len(kptfile) - 3
          kptfile(iend:iend+3) = '.kpt'
        endif
      elseif (iwc.eq.33) then
        dipfile = words(i)
        if (index(dipfile,'.dip').eq.0) then
          call endstring(dipfile,len(dipfile),iend)
          if (iend.eq.0) iend = len(dipfile) - 3
          dipfile(iend:iend+3) = '.dip'
        endif
      elseif (iwc.eq.34) then
        xsffile = words(i)
        if (index(xsffile,'.xsf').eq.0) then
          call endstring(xsffile,len(xsffile),iend)
          if (iend.eq.0) iend = len(xsffile) - 3
          xsffile(iend:iend+3) = '.xsf'
        endif
        lxsfsym = llocalsym
        lxsfseq = llocalseq
!
! CML output
!
      elseif (iwc.eq.1000) then
        cmlfilename = words(i)
        if (index(cmlfilename,'.xml').eq.0) then
          call endstring(cmlfilename,len(cmlfilename),iend)
          if (iend.eq.0) iend = len(cmlfilename) - 3
          cmlfilename(iend:iend+3) = '.xml'
        endif
      elseif (iwc.eq.1001) then
        cmlfilename = words(i)
        if (index(cmlfilename,'.xml').eq.0) then
          call endstring(cmlfilename,len(cmlfilename),iend)
          if (iend.eq.0) iend = len(cmlfilename) - 3
          cmlfilename(iend:iend+3) = '.xml'
        endif
      else
        nwarn = nwarn + 1
        call outwarning('invalid filetype supplied in output option',0_i4)
      endif
    endif
  enddo
  if (((nword-nextraword).eq.2).and.nfloat.ge.1) then
!
!  File name was not given but there is a numeric file name
!
    if (nfloat.ge.2) then
      word = floatwords(2)
    else
      word = floatwords(1)
    endif
    lset = .false.
    lnodot = (index(word,'.').eq.0)
!
!  Check that there is no dot in the word, 
!  otherwise it's a file name not an option
!
    if (.not.lset) then
!
!  Assign any output files that may have been given
!
      if (iwc.eq.1) then
        marvfile = word
        if (index(marvfile,'.mvn').eq.0) then
          call endstring(marvfile,len(marvfile),iend)
          if (iend.eq.0) iend = len(marvfile) - 3
          marvfile(iend:iend+3) = '.mvn'
        endif
      elseif (iwc.eq.2) then
        thbfile = word
      elseif (iwc.eq.3) then
        xtlfile = word
        if (index(xtlfile,'.xtl').eq.0) then
          call endstring(xtlfile,len(xtlfile),iend)
          if (iend.eq.0) iend = len(xtlfile) - 3
          xtlfile(iend:iend+3) = '.xtl'
        endif
      elseif (iwc.eq.4) then
        phonfile = word
      elseif (iwc.eq.5) then
        arcfile = word
        if (index(arcfile,'.arc').eq.0.and.index(arcfile,'.car').eq.0) then
          call endstring(arcfile,len(arcfile),iend)
          if (iend.eq.0) iend = len(arcfile) - 3
          arcfile(iend:iend+3) = '.arc'
        endif
        lmovie = llocalmovie
        if (nfloat.gt.1.and.lmovie) then
          narcwrite = nint(abs(floats(1)))
          narcwrite = max(narcwrite,1_i4)
        endif
      elseif (iwc.eq.6) then
        xrfile = word
        if (index(xrfile,'.xr').eq.0) then
          iend = index(xrfile,' ')
          if (iend.eq.0) iend = len(xrfile) - 2
          xrfile(iend:iend+2) = '.xr'
        endif
      elseif (iwc.eq.7) then
        cssrfile = word
        if (index(cssrfile,'.cssr').eq.0) then
          call endstring(cssrfile,len(cssrfile),iend)
          if (iend.eq.0) iend = len(cssrfile) - 4
          cssrfile(iend:iend+4) = '.cssr'
        endif
      elseif (iwc.eq.8) then
        trjfile = word
        if (index(trjfile,'.trg').eq.0) then
          call endstring(trjfile,len(trjfile),iend)
          if (iend.eq.0) iend = len(trjfile) - 3
          trjfile(iend:iend+3) = '.trg'
        endif
      elseif (iwc.eq.9) then
        freqfile = word
      elseif (iwc.eq.10) then
        xyzfile = word
        if (index(xyzfile,'.xyz').eq.0) then
          call endstring(xyzfile,len(xyzfile),iend)
          if (iend.eq.0) iend = len(xyzfile) - 3
          xyzfile(iend:iend+3) = '.xyz'
        endif
        lxyzmovie = llocalmovie
      elseif (iwc.eq.11) then
        hisfile = word
        if (index(hisfile,'.his').eq.0) then
          call endstring(hisfile,len(hisfile),iend)
          if (iend.eq.0) iend = len(hisfile) - 3
          hisfile(iend:iend+3) = '.his'
        endif
      elseif (iwc.eq.12) then
        fdffile = word
        if (index(fdffile,'.fdf').eq.0) then
          call endstring(fdffile,len(fdffile),iend)
          if (iend.eq.0) iend = len(fdffile) - 3
          fdffile(iend:iend+3) = '.fdf'
        endif
      elseif (iwc.eq.13) then
        drvfile = word
        if (index(drvfile,'.drv').eq.0) then
          call endstring(drvfile,len(drvfile),iend)
          if (iend.eq.0) iend = len(drvfile) - 3
          drvfile(iend:iend+3) = '.drv'
        endif
      elseif (iwc.eq.14) then
        frcfile = word
        if (index(frcfile,'.frc').eq.0) then
          call endstring(frcfile,len(frcfile),iend)
          if (iend.eq.0) iend = len(frcfile) - 3
          frcfile(iend:iend+3) = '.frc'
        endif
      elseif (iwc.eq.15) then
        ciffile = word
        if (index(ciffile,'.cif').eq.0) then
          call endstring(ciffile,len(ciffile),iend)
          if (iend.eq.0) iend = len(ciffile) - 3
          ciffile(iend:iend+3) = '.cif'
        endif
        if (nfloat.gt.0) then
          cif_dummylattice = abs(floats(1))
        endif
      elseif (iwc.eq.16) then
        dlvfile = word
        if (index(dlvfile,'.str').eq.0) then
          call endstring(dlvfile,len(dlvfile),iend)
          if (iend.eq.0) iend = len(dlvfile) - 3
          dlvfile(iend:iend+3) = '.str'
        endif
      elseif (iwc.eq.17) then
        eigfile = word
        if (index(eigfile,'.eig').eq.0) then
          call endstring(eigfile,len(eigfile),iend)
          if (iend.eq.0) iend = len(eigfile) - 3
          eigfile(iend:iend+3) = '.eig'
        endif
      elseif (iwc.eq.19) then
        prefile = word
        if (index(prefile,'.pre').eq.0) then
          call endstring(prefile,len(prefile),iend)
          if (iend.eq.0) iend = len(prefile) - 3
          prefile(iend:iend+3) = '.pre'
        endif
      elseif (iwc.eq.20) then
        sasfile = word
        if (index(sasfile,'.sas').eq.0) then
          call endstring(sasfile,len(sasfile),iend)
          if (iend.eq.0) iend = len(sasfile) - 3
          sasfile(iend:iend+3) = '.sas'
        endif
      elseif (iwc.eq.21) then
        oscfile = word
        if (index(oscfile,'.osc').eq.0) then
          call endstring(oscfile,len(oscfile),iend)
          if (iend.eq.0) iend = len(oscfile) - 3
          oscfile(iend:iend+3) = '.osc'
        endif
      elseif (iwc.eq.22) then
        biofile = word
        if (index(biofile,'.bio').eq.0) then
          call endstring(biofile,len(biofile),iend)
          if (iend.eq.0) iend = len(biofile) - 3
          biofile(iend:iend+3) = '.bio'
        endif
      elseif (iwc.eq.23) then
        pdffiles(ncurr) = word
        if (index(pdffiles(ncurr),'.wid').eq.0) then
          iend = index(pdffiles(ncurr),' ')
          call endstring(pdffiles(ncurr),len(pdffiles(ncurr)),iend)
          if (iend.eq.0) iend = len(pdffiles(ncurr)) - 3
          pdffiles(ncurr)(iend:iend+3) = '.wid'
        endif
      elseif (iwc.eq.24) then
        cosmofile = word
        if (index(cosmofile,'.cosmo').eq.0) then
          call endstring(cosmofile,len(cosmofile),iend)
          if (iend.eq.0) iend = len(cosmofile) - 5
          cosmofile(iend:iend+5) = '.cosmo'
        endif
      elseif (iwc.eq.25) then
        qbofile = word
        if (index(qbofile,'.qbo').eq.0) then
          call endstring(qbofile,len(qbofile),iend)
          if (iend.eq.0) iend = len(qbofile) - 3
          qbofile(iend:iend+5) = '.qbo'
        endif
      elseif (iwc.eq.27) then
        dcdfile = word
        if (index(dcdfile,'.dcd').eq.0) then
          call endstring(dcdfile,len(dcdfile),iend)
          if (iend.eq.0) iend = len(dcdfile) - 3
          dcdfile(iend:iend+5) = '.dcd'
        endif
      elseif (iwc.eq.28) then
        lammpsfile = words(i)
        if (index(lammpsfile,'.lmp').eq.0) then
          call endstring(lammpsfile,len(lammpsfile),iend)
          if (iend.eq.0) iend = len(lammpsfile) - 3
          lammpsfile(iend:iend+3) = '.lmp'
        endif
      elseif (iwc.eq.29) then
        casfile = word
        if (index(casfile,'.csp').eq.0) then
          call endstring(casfile,len(casfile),iend)
          if (iend.eq.0) iend = len(casfile) - 3
          casfile(iend:iend+3) = '.csp'
        endif
      elseif (iwc.eq.30) then
        shengfile = word
      elseif (iwc.eq.31) then
        inertfile = word
        if (index(inertfile,'.ine').eq.0) then
          call endstring(inertfile,len(inertfile),iend)
          if (iend.eq.0) iend = len(inertfile) - 3
          inertfile(iend:iend+3) = '.ine'
        endif
      elseif (iwc.eq.32) then
        kptfile = word
        if (index(kptfile,'.kpt').eq.0) then
          call endstring(kptfile,len(kptfile),iend)
          if (iend.eq.0) iend = len(kptfile) - 3
          kptfile(iend:iend+3) = '.kpt'
        endif
      elseif (iwc.eq.33) then
        dipfile = word
        if (index(dipfile,'.dip').eq.0) then
          call endstring(dipfile,len(dipfile),iend)
          if (iend.eq.0) iend = len(dipfile) - 3
          dipfile(iend:iend+3) = '.dip'
        endif
!
! CML output
!
      elseif (iwc.eq.1000) then
        cmlfilename = word
        if (index(cmlfilename,'.xml').eq.0) then
          call endstring(cmlfilename,len(cmlfilename),iend)
          if (iend.eq.0) iend = len(cmlfilename) - 3
          cmlfilename(iend:iend+3) = '.xml'
        endif
      elseif (iwc.eq.1001) then
        cmlfilename = word
        if (index(cmlfilename,'.xml').eq.0) then
          call endstring(cmlfilename,len(cmlfilename),iend)
          if (iend.eq.0) iend = len(cmlfilename) - 3
          cmlfilename(iend:iend+3) = '.xml'
        endif
      else
        nwarn = nwarn + 1
        call outwarning('invalid filetype supplied in output option',0_i4)
      endif
    endif
  endif
  lwordok = .true.
  return
!********************************************************************
!  Relative speed of reciprocal and real space electrostatic terms  *
!********************************************************************
520 if (nfloat.ge.1) then
    rspeed0 = abs(floats(1))
  else
    call outerror('rspeed value(s) missing from input',iline)
    call stopnow('genword')
  endif
  lwordok = .true.
  return
!*********
!  Nadd  *
!*********
530 if (nfloat.gt.0) then
    nadd = int(floats(1))
  else
    read(nru,*,err=99,end=99) nadd
    iline = iline + 1
  endif
  lwordok = .true.
  return
!*********************************************************
!  gdcrit - for switch between energy and force balance  *
!*********************************************************
540 if (nfloat.gt.0) then
    gdcrit = floats(1)
  else
    read(nru,*,err=99,end=99) gdcrit
    iline = iline + 1
  endif
  gdcrit = abs(gdcrit)
  lwordok = .true.
  return
!******************************************************************
!  No bond - prevent molecule option from locating bonds between  *
!  specified atom types                                           *
!  Save index as six figure number - each three digits implies    *
!  atom type                                                      *
!******************************************************************
550 nnobo = nnobo + 1
  if (nnobo.gt.maxnobo) then
    maxnobo = nnobo + 10
    call changemaxnobo
  endif
  if (nfloat.ge.2) then
    nvar1 = nint(floats(1))
    if (nvar1.gt.100) nvar1 = nvar1 - 100 + maxele
    nvar2 = nint(floats(2))
    if (nvar2.gt.100) nvar2 = nvar2 - 100 + maxele
    itype1 = 0
    itype2 = 0
  elseif (nword.eq.3) then
    call ltont(words(2),nvar1,itype1)
    call ltont(words(3),nvar2,itype2)
  elseif (nword.eq.4) then
    call ltont(words(2),nvar1,itype1)
    word = words(3)
    call stolc(word,maxwordlength)
    if (index(word,'cor').eq.1.or.index(word,'bco').eq.1) then
      call ltont(words(4),nvar2,itype2)
    elseif (index(word,'she').eq.1.or.index(word,'bsh').eq.1) then
      nvar1 = nvar1 + maxele
      call ltont(words(4),nvar2,itype2)
    else
      call ltont(words(3),nvar2,itype2)
      word = words(4)
      call stolc(word,maxwordlength)
      if (index(word,'she').eq.1.or.index(word,'bsh').eq.1) then
        nvar2 = nvar2 + maxele
      endif
    endif
  elseif (nword.ge.5) then
    call ltont(words(2),nvar1,itype1)
    call ltont(words(4),nvar2,itype2)
    if ((index(words(3),'s').eq.1).or.(index(words(3),'S').eq.1)) nvar1 = nvar1 + maxele
    if ((index(words(5),'s').eq.1).or.(index(words(5),'S').eq.1)) nvar2 = nvar2 + maxele
    if ((index(words(3),'bs').eq.1).or.(index(words(3),'BS').eq.1)) nvar1 = nvar1 + maxele
    if ((index(words(5),'bs').eq.1).or.(index(words(5),'BS').eq.1)) nvar2 = nvar2 + maxele
  else
    line = '  '
    read(nru,'(a)',err=99,end=99) line
    iline = iline + 1
    call linepro(nru,line,iline)
    if (nfloat.ge.2) then
      nvar1 = nint(floats(1))
      if (nvar1.gt.100) nvar1 = nvar1 - 100 + maxele
      nvar2 = nint(floats(2))
      if (nvar2.gt.100) nvar2 = nvar2 - 100 + maxele
      itype1 = 0
      itype2 = 0
    elseif (nword.eq.2) then
      call ltont(words(1),nvar1,itype1)
      call ltont(words(2),nvar2,itype2)
    elseif (nword.eq.3) then
      call ltont(words(1),nvar1,itype1)
      word = words(2)
      call stolc(word,maxwordlength)
      if (index(word,'cor').eq.1.or.index(word,'bco').eq.1) then
        call ltont(words(3),nvar2,itype2)
      elseif (index(word,'she').eq.1.or.index(word,'bsh').eq.1) then
        nvar1 = nvar1 + maxele
        call ltont(words(3),nvar2,itype2)
      else
        call ltont(words(2),nvar2,itype2)
        word = words(3)
        call stolc(word,maxwordlength)
        if (index(word,'she').eq.1.or.index(word,'bsh').eq.1) then
          nvar2 = nvar2 + maxele
        endif
      endif
    elseif (nword.ge.4) then
      call ltont(words(1),nvar1,itype1)
      call ltont(words(3),nvar2,itype2)
      if ((index(words(2),'s').eq.1).or.(index(words(2),'S').eq.1)) nvar1 = nvar1 + maxele
      if ((index(words(4),'s').eq.1).or.(index(words(4),'S').eq.1)) nvar2 = nvar2 + maxele
      if ((index(words(2),'bs').eq.1).or.(index(words(2),'BS').eq.1)) nvar1 = nvar1 + maxele
      if ((index(words(4),'bs').eq.1).or.(index(words(4),'BS').eq.1)) nvar2 = nvar2 + maxele
    else
      call outerror('Error in species input for nobond option',iline)
      call stopnow('genword')
    endif
  endif
  if (nvar1.eq.nvar2) then
    nvar2 = nvar2 + 1000*nvar1
    nobond(nnobo) = nvar2
    if (itype1.lt.itype2) then
      itype2 = itype2 + 1000*itype1
      nobotyp(nnobo) = itype2
    else
      itype1 = itype1 + 1000*itype2
      nobotyp(nnobo) = itype1
    endif
  elseif (nvar1.lt.nvar2) then
    nvar2 = nvar2 + 1000*nvar1
    nobond(nnobo) = nvar2
    itype2 = itype2 + 1000*itype1
    nobotyp(nnobo) = itype2
  else
    nvar1 = nvar1 + 1000*nvar2
    nobond(nnobo) = nvar1
    itype1 = itype1 + 1000*itype2
    nobotyp(nnobo) = itype1
  endif
  lwordok = .true.
  return
!****************************************************
!  Number of cycles before enforced Hessian update  *
!****************************************************
560 if (nword.gt.1) then
    if (index(words(2),'opt').eq.1) then
!
!  Set only optimisation update
!
      if (nfloat.gt.0) then
        nupdate = nint(floats(1))
      else
        read(nru,*,err=99,end=99) nupdate
        iline = iline + 1
      endif
    elseif (index(words(2),'fit').eq.1) then
!
!  Set only fit update
!
      if (nfloat.gt.0) then
        nfupdate = nint(floats(1))
      else
        read(nru,*,err=99,end=99) nfupdate
        iline = iline + 1
      endif
    else
!
!  Set according to job type
!
      if (lopt) then
        if (nfloat.gt.0) then
          nupdate = nint(floats(1))
        else
          read(nru,*,err=99,end=99) nupdate
          iline = iline + 1
        endif
        if (lfit) nfupdate = nupdate
      elseif (lfit) then
        if (nfloat.gt.0) then
          nfupdate = nint(floats(1))
        else
          read(nru,*,err=99,end=99) nfupdate
          iline = iline + 1
        endif
      endif
    endif
  else
!
!  Set according to job type
!
    if (lopt) then
      if (nfloat.gt.0) then
        nupdate = nint(floats(1))
      else
        read(nru,*,err=99,end=99) nupdate
        iline = iline + 1
      endif
      if (lfit) nfupdate = nupdate
    elseif (lfit) then
      if (nfloat.gt.0) then
        nfupdate = nint(floats(1))
      else
        read(nru,*,err=99,end=99) nfupdate
        iline = iline + 1
      endif
    endif
  endif
  if (nupdate.le.0) then
    nwarn = nwarn + 1
    if (ioproc) then
      write(ioout,'(/,''  **** Warning - updating of Hessian every zero cycles is meaningless ****'')')
      write(ioout,'(''  **** The update parameter will be reset to the default value        ****'')')
    endif
    nupdate = 10
  endif
  if (nfupdate.le.0) then
    nwarn = nwarn + 1
    if (ioproc) then
      write(ioout,'(/,''  **** Warning - updating of Hessian every zero cycles is meaningless ****'')')
      write(ioout,'(''  **** The update parameter will be reset to the default value        ****'')')
    endif
    nfupdate = 20
  endif
  lwordok = .true.
  return
!***********************************************
!  Order of the limited memory BFGS algorithm  *
!***********************************************
570 if (nfloat.ge.1) then
    lmbfgsorder = nint(floats(1))
  else
    read(nru,*,err=99,end=99) lmbfgsorder
    iline=iline+1
  endif
  lwordok = .true.    
  return
!**************************************************
!  Maximum number of points in line minimisation  *
!**************************************************
580 if (nfloat.ge.1) then
    nlinmin = nint(floats(1))
  else
    read(nru,*,err=99,end=99) nlinmin
    iline=iline+1
  endif
  lwordok = .true.
  return
!********************
!  Potential sites  *
!********************
590 units = 1.0_dp
  lcart = (ndimen(ncurr).eq.0)
  if (nword.gt.1) then
    do i = 2,nword
      if (index(words(i),'cart').eq.1) lcart = .true.
      if (index(words(i),'au').eq.1) units = autoangs
    enddo
  endif
595 line = '  '
  read(nru,'(a)',err=1000,end=1000) line
  iline = iline + 1
  call linepro(nru,line,iline)
  if (nword.gt.0) then
    l55 = .true.
    return
  elseif ((nword+nfloat).eq.0) then
    goto 595
  endif
  npotsites = npotsites + 1
  npotsitescfg(ncurr) = npotsitescfg(ncurr) + 1
  lpotsitecartcfg(ncurr) = lcart
  if (npotsites.gt.maxpotsites) then
    maxpotsites = npotsites + 20
    call changemaxpotsites
  endif
  if (nfloat.lt.3) then
    call outerror('Insufficient coordinates specified',iline)
    call stopnow('genword')
  endif
  npotsitecfg(npotsites) = ncurr
  if (lcart) then
    xpotsite(npotsites) = floats(1)*units
    ypotsite(npotsites) = floats(2)*units
    zpotsite(npotsites) = floats(3)*units
  else
    if (ndimen(ncurr).eq.1) then
      xpotsite(npotsites) = floats(1)
      ypotsite(npotsites) = floats(2)*units
      zpotsite(npotsites) = floats(3)*units
    elseif (ndimen(ncurr).eq.2) then
      xpotsite(npotsites) = floats(1)
      ypotsite(npotsites) = floats(2)
      zpotsite(npotsites) = floats(3)*units
    elseif (ndimen(ncurr).eq.3) then
      xpotsite(npotsites) = floats(1)
      ypotsite(npotsites) = floats(2)
      zpotsite(npotsites) = floats(3)
    endif
  endif
  goto 595
!****************************
!  Potential interpolation  *
!****************************
600 if (nfloat.gt.0) then
    lpotlinterpolate = .true.
    nptsinterpolate = nint(abs(floats(1)))
  endif
  lwordok = .true.
  return
!************************************************************
!  Maximum change in function before recalculating Hessian  *
!************************************************************
610 if (nfloat.ge.1) then
    delfc = floats(1)
  else
    read(nru,*,err=99,end=99) delfc
    iline = iline + 1
  endif
  lwordok = .true.
  return
!****************************************
!  Marvin insert for Marvin input file  *
!****************************************
620 if (lmarv) then
    ind = index(marvfile,'.mvn') + 1
    if (ind.eq.0) then
      marvtemp = 'marvin.gmt'
    else
      marvtemp = marvfile
      marvtemp(ind:ind+3) = 'gmt'
    endif
  else
    marvtemp = 'marvin.gmt'
  endif
  open(7,file=marvtemp,status='unknown',err=625)
623 read(nru,'(a)',end=625,err=625) line
  if (index(line,'end').eq.1) goto 625
  if (ioproc) then
    write(7,'(a)') trim(line)
  endif
  goto 623
625 close(7)
  lwordok = .true.
  return
!***********************
!  QM/MM mode control  *
!***********************
630 if (nword.ge.2) then
    call stolc(words(2),maxword)
    if (index(words(2),'me').eq.1) then
      QMMMmode(ncurr) = 1
    elseif (index(words(2),'el').eq.1) then
      QMMMmode(ncurr) = 2
    endif
  endif
  lwordok = .true.
  return
!***********************
!  Units modification  *
!***********************
640 if (nword.ge.2.and.nfloat.ge.1) then
    call stolc(words(2),maxword)
    if (index(words(2),'kcaltoev').eq.1) then
      kcaltoev = abs(floats(1))
    elseif (index(words(2),'angstoev').eq.1) then
      inverse_angstroms_to_ev = abs(floats(1))
    endif
  endif
  lwordok = .true.
  return
!*******************************
!  Selection of maths library  *
!*******************************
650 if (nword.ge.2) then
    call stolc(words(2),maxword)
    if (index(words(2),'lap').eq.1) then
      leispack_eigensolve = .false.
      if (nword.ge.3) then
        call stolc(words(3),maxword)
        if (index(words(3),'nod').eq.1) ldivide_and_conquer = .false.
      endif
    elseif (index(words(2),'eis').eq.1) then
      leispack_eigensolve = .true.
    endif
  endif
  lwordok = .true.
  return
!************************
!  Parallel algorithms  *
!************************
660 if (nword.ge.2) then
    call stolc(words(2),maxword)
    if (index(words(2),'avoi').eq.1) then
      lfastMPI = .false.
    endif
  endif
  lwordok = .true.
  return
!*****************************************
!  Selection of iterative charge solver  *
!*****************************************
!
!  qsolver = 1 => itpack
!  qsolver = 2 => lapack
!  qsolver = 3 => (d)bcgsolve
!
670 if (nword.ge.2) then
    call stolc(words(2),maxword)
    if (index(words(2),'lap').eq.1) then
      qsolver = 2
    elseif (index(words(2),'bcg').eq.1) then
      qsolver = 3
    elseif (index(words(2),'itp').eq.1) then
      qsolver = 1
    endif
  endif
  lwordok = .true.
  return
!******************
!  Stress tensor  *
!******************
680 units = 1.0_dp
  if (nword.ge.2) then
    call stolc(words(2),maxword)
    if (index(words(2),'pa ').eq.1) then
      units = 1.0d-9
    elseif (index(words(2),'kpa ').eq.1) then
      units = 1.0d-6
    elseif (index(words(2),'atm').eq.1) then
      units = 101.325d-6
    elseif (index(words(2),'nm-2').eq.1) then
      units = 1.0d-9
    elseif (index(words(2),'kbar').eq.1) then
      units = 0.1_dp
    endif
  endif
  if (nfloat.lt.6) then
    line = '  '
    read(nru,'(a)')line
    iline = iline + 1
    call linepro(nru,line,iline)
    if (nfloat.ge.6) then
      if (nword.ge.1) then
        call stolc(words(1),maxword)
        if (index(words(1),'pa ').eq.1) then
          units = 1.0d-9
        elseif (index(words(1),'kpa ').eq.1) then
          units = 1.0d-6
        elseif (index(words(1),'atm').eq.1) then
          units = 101.325d-6
        elseif (index(words(1),'nm-2').eq.1) then
          units = 1.0d-9
        elseif (index(words(1),'kbar').eq.1) then
          units = 0.1_dp
        endif
      endif
    else
      call outerror('Stress tensor is missing or incomplete',iline)
      call stopnow('genword')
    endif
  endif
  stresscfg(1,ncurr) = floats(1)*units
  stresscfg(2,ncurr) = floats(2)*units
  stresscfg(3,ncurr) = floats(3)*units
  stresscfg(4,ncurr) = floats(4)*units
  stresscfg(5,ncurr) = floats(5)*units
  stresscfg(6,ncurr) = floats(6)*units
!
!  Set pressure according to average of on diagonal elements of stress tensor
!
  presscfg(ncurr) = (stresscfg(1,ncurr)+stresscfg(2,ncurr)+stresscfg(3,ncurr))
  nonzero = 0
  if (abs(stresscfg(1,ncurr)).gt.1.0d-12) nonzero = nonzero + 1
  if (abs(stresscfg(2,ncurr)).gt.1.0d-12) nonzero = nonzero + 1
  if (abs(stresscfg(3,ncurr)).gt.1.0d-12) nonzero = nonzero + 1
  if (nonzero.gt.0) presscfg(ncurr) = presscfg(ncurr)/dble(nonzero)
  lwordok = .true.
  return
!******************************
!  Maximise using RFO method  *
!******************************
690 if (nword.gt.1) then
    if (index(words(2),'mode').eq.1) then
      if (morder.gt.1) then
        call outerror('Mode following only allowed for 1st order TS',iline)
        call stopnow('genword')
      endif
      if (nfloat.ge.1) then
        mode = nint(floats(1))
        morder = 1
      else
        nwarn = nwarn + 1
        call outwarning('no mode given for maximisation',iline)
      endif
    elseif (index(words(2),'orde').eq.1) then
      if (nfloat.ge.1) then
        morder = nint(floats(1))
      else
        nwarn = nwarn + 1
        call outwarning('no mode given for maximisation',iline)
      endif
      if (mode.gt.0.and.morder.gt.1) then
        call outerror('Mode following only allowed for 1st order TS',iline)
        call stopnow('genword')
      endif
    endif
  endif
  lwordok = .true.
  return
!******************************
!  Pressure of configuration  *
!******************************
700 if (nfloat.ge.1) then
    press = floats(1)
    if (nword.ge.2) then
      call stolc(words(2),maxword)
      if (index(words(2),'pa ').eq.1) then
        press = press*1.0d-9
      elseif (index(words(2),'kpa ').eq.1) then
        press = press*1.0d-6
      elseif (index(words(2),'mpa ').eq.1) then
        press = press*1.0d-3
      elseif (index(words(2),'atm').eq.1) then
        press = press*101.325d-6
      elseif (index(words(2),'nm-2').eq.1) then
        press = press*1.0d-9
      elseif (index(words(2),'kbar').eq.1) then
        press = press*0.1_dp
      elseif (index(words(2),'bar').eq.1) then
        press = press*0.0001_dp
      endif
    endif
  else
    line = '  '
    read(nru,'(a)') line
    iline = iline + 1
    call linepro(nru,line,iline)
    if (nfloat.ge.1) then
      press = floats(1)
      if (nword.ge.1) then
        call stolc(words(1),maxword)
        if (index(words(1),'pa ').eq.1) then
          press = press*1.0d-9
        elseif (index(words(1),'kpa ').eq.1) then
          press = press*1.0d-6
        elseif (index(words(1),'mpa ').eq.1) then
          press = press*1.0d-3
        elseif (index(words(1),'atm').eq.1) then
          press = press*101.325d-6
        elseif (index(words(1),'nm-2').eq.1) then
          press = press*1.0d-9
        elseif (index(words(1),'kbar').eq.1) then
          press = press*0.1_dp
        endif
      endif
    else
      call outerror('Pressure value is missing from input',iline)
      call stopnow('genword')
    endif
  endif
  presscfg(ncurr) = press
!
!  Set pressure in stress tensor
!
  stresscfg(1,ncurr) = press
  stresscfg(2,ncurr) = press
  stresscfg(3,ncurr) = press
  stresscfg(4,ncurr) = 0.0_dp
  stresscfg(5,ncurr) = 0.0_dp
  stresscfg(6,ncurr) = 0.0_dp
  lwordok = .true.
  return
!****************************
!  Parameters for Wolf sum  *
!****************************
710 if (nfloat.ge.2) then
    etaw = abs(floats(1))
    cutw = abs(floats(2))
    if (nword.gt.1) then
      lwolforiginal = (index(words(2),'ori').eq.1)
      lwolffennell  = (index(words(2),'fen').eq.1)
    endif
  else
    call outerror('insufficient parameters in input for qwolf',iline)
    call stopnow('genword')
  endif
!
!  Set Wolf flag
!
  lwolf = .true.
!
!  Check that this is not a defect calculation
!
  if (ldefect) then
    call outerror('Wolf sum cannot be used in a defect calculation yet',iline)
    call stopnow('genword')
  endif
!
!  Set tweatpi
!
  tweatpi = 2.0_dp*etaw/sqrtpi
  if (abs(cutw).gt.1.0d-12) then
    rkw = 1.0_dp/cutw
    selfwolf = g_derfc(etaw*cutw)*rkw 
    if (lwolffennell) then
      selfwolfforce  = selfwolf*rkw + tweatpi*exp(-(etaw*cutw)**2)*rkw
      selfwolfenergy = selfwolfforce
    elseif (lwolforiginal) then
      selfwolfforce  = selfwolf*rkw + tweatpi*exp(-(etaw*cutw)**2)*rkw
      selfwolfenergy = 0.0_dp
    else
      selfwolfforce  = 0.0_dp
      selfwolfenergy = 0.0_dp
    endif
  endif
!
  lwordok = .true.
  return
!*****************
!  Radial force  *
!*****************
720 if (nfloat.ge.4) then
    radialKcfg(ncurr) = floats(1)
    radialXYZcfg(1,ncurr) = floats(2)
    radialXYZcfg(2,ncurr) = floats(3)
    radialXYZcfg(3,ncurr) = floats(4)
  else
    call outerror('insufficient parameters in input for radial_force',iline)
    call stopnow('genword')
  endif
!
!  Check that this configuration is non-periodic
!
  if (ndimen(ncurr).gt.0) then
    call outerror('Radial force cannot be used for periodic systems',iline)
    call stopnow('genword')
  endif
!
!  Set radial force flag
!
  lradialcfg(ncurr) = .true.
!
  lwordok = .true.
  return
!************************************
!  Spatial decomposition cutoff(s)  *
!************************************
730 continue
  if (nword.gt.1) then
    lrcspatial_anisotropic = (index(words(2),'ani').eq.1) 
    lrcspatialBO_anisotropic = lrcspatial_anisotropic
  endif
  if (lrcspatial_anisotropic) then
    if (nfloat.ge.6) then  
      rcspatialx   = abs(floats(1))
      rcspatialy   = abs(floats(2))
      rcspatialz   = abs(floats(3))
      rcspatialbox = abs(floats(4))
      rcspatialboy = abs(floats(5))
      rcspatialboz = abs(floats(6))
    elseif (nfloat.ge.3) then
      rcspatialx   = abs(floats(1))
      rcspatialy   = abs(floats(2))
      rcspatialz   = abs(floats(3))
      rcspatialbox = 0.0_dp
      rcspatialboy = 0.0_dp
      rcspatialboz = 0.0_dp
    endif
  else
    if (nfloat.ge.2) then  
      rcspatial   = abs(floats(1))
      rcspatialbo = abs(floats(2))
    elseif (nfloat.eq.1) then
      rcspatial   = abs(floats(1))
      rcspatialbo = rcspatial
    endif
  endif
  lwordok = .true.
  return
!************************************
!  Gasteiger convergance tolerance  *
!************************************
740 if (nfloat.gt.0) then
    gasttol = floats(1)
  else
    read(nru,*,err=99,end=99) gasttol
    iline = iline + 1
  endif
  gasttol = max(abs(gasttol),1.0d-10)
  lwordok = .true.
  return
!*******************************************
!  Gasteiger maximum number of iterations  *
!*******************************************
750 if (nfloat.gt.0) then
    ngastitermax = nint(floats(1))
  else
    read(nru,*,err=99,end=99) ngastitermax
    iline = iline + 1
  endif
  ngastitermax = max(1,abs(ngastitermax))
  lwordok = .true.
  return
!*****************************************************
!  Finite difference interval for numerical phonons  *
!*****************************************************
760 if (nfloat.gt.0) then
    phondiff = abs(floats(1))
  else
    phondiff = 1.0d-5
  endif
  lwordok = .true.
  return
!****************************
!  Parameters for Wolf sum  *
!****************************
770 if (nfloat.ge.2) then
    etaw = abs(floats(1))
    cutw = abs(floats(2))
    if (nword.gt.1) then
      lwolforiginal = (index(words(1),'ori').eq.1)
    endif
  else
    call outerror('insufficient parameters in input for qwolf',iline)
    call stopnow('genword')
  endif
!
!  Set Wolf flag
!
  lwolf = .true.
!
!  Check that this is not a defect calculation
!
  if (ldefect) then
    call outerror('Wolf sum cannot be used in a defect calculation yet',iline)
    call stopnow('genword')
  endif
!  
!  Set tweatpi
!     
  tweatpi = 2.0_dp*etaw/sqrtpi
  if (abs(cutw).gt.1.0d-10) then
    rkw = 1.0_dp/cutw
    selfwolf = g_derfc(etaw*cutw)*rkw
  endif
!     
  lwordok = .true.
  return
!*****************************************
!  Parameters for COSMO/COSMIC Wolf sum  *
!*****************************************
780 if (nfloat.ge.2) then
    etawc = abs(floats(1))
    cutwc = abs(floats(2))
  else
    call outerror('insufficient parameters in input for cwolf',iline)
    call stopnow('genword')
  endif
!
!  Set tweatpi
!
  tweatpic = 2.0_dp*etawc/sqrtpi
  if (abs(cutwc).gt.1.0d-10) then
    selfwolfc = g_derfc(etawc*cutwc)/cutwc
  endif
!
  lwordok = .true.
  return
!*****************************************************************
!  Finite difference interval for numerical property evaluation  *
!*****************************************************************
790 if (nfloat.ge.2) then
    findiffc = abs(floats(1))
    findiffs = abs(floats(2))
  elseif (nfloat.eq.1) then
    findiffc = abs(floats(1))
  endif
  lwordok = .true.
  return
!*************************
!  Anisotropic pressure  *
!*************************
800 units = 1.0_dp
  if (nword.ge.2) then
    call stolc(words(2),maxword)
    if (index(words(2),'pa ').eq.1) then
      units = 1.0d-9
    elseif (index(words(2),'kpa ').eq.1) then
      units = 1.0d-6
    elseif (index(words(2),'atm').eq.1) then
      units = 101.325d-6
    elseif (index(words(2),'nm-2').eq.1) then
      units = 1.0d-9
    elseif (index(words(2),'kbar').eq.1) then
      units = 0.1_dp
    endif
  endif
!
!  Check that cell is not isotropic since this is incompatible
!
  if (index(keyword,'iso').ne.0) then
    call outerror('anisotropic pressure cannot be used with an isotropic cell',iline)
    call stopnow('genword')
  endif
  if (ndimen(ncurr).eq.3) then
    if (nfloat.lt.6) then
      call outerror('anisotropic pressure input is incomplete',iline)
      call stopnow('genword')
    endif
    anisotropicpresscfg(1,ncurr) = floats(1)*units
    anisotropicpresscfg(2,ncurr) = floats(2)*units
    anisotropicpresscfg(3,ncurr) = floats(3)*units
    anisotropicpresscfg(4,ncurr) = floats(4)*units
    anisotropicpresscfg(5,ncurr) = floats(5)*units
    anisotropicpresscfg(6,ncurr) = floats(6)*units
  elseif (ndimen(ncurr).eq.2) then
    if (nfloat.lt.3) then
      call outerror('anisotropic pressure input is incomplete',iline)
      call stopnow('genword')
    endif
    anisotropicpresscfg(1,ncurr) = floats(1)*units
    anisotropicpresscfg(2,ncurr) = floats(2)*units
    anisotropicpresscfg(3,ncurr) = floats(3)*units
  elseif (ndimen(ncurr).eq.1) then
    if (nfloat.lt.1) then
      call outerror('anisotropic pressure input is incomplete',iline)
      call stopnow('genword')
    endif
    anisotropicpresscfg(1,ncurr) = floats(1)*units
  endif
  lanisotropicpresscfg(ncurr) = .true.
!
!  Since an anisotropic pressure has been set then we forced to use force minimisation
!
  lforcemin = .true.
!
  lwordok = .true.
  return
!**********************
!  Plumed log option  *
!**********************
805 continue
  if (.not.lplumed_available) then
    call outerror('PLUMED run but GULP has not been built with PLUMED'//word,iline)
    call stopnow('genword')
  endif
  if (nword.gt.1) then
    plumedlog = words(2)
  else
    plumedlog = 'plumed.log'
  endif
  lwordok = .true.
  return
!******************
!  Plumed option  *
!******************
810 continue
  if (.not.lplumed_available) then
    call outerror('PLUMED run but GULP has not been built with PLUMED'//word,iline)
    call stopnow('genword')
  endif
  lplumed = .true.
  if (nword.gt.1) then
    plumedinput = words(2)
  else
    plumedinput = 'plumed.dat'
  endif
  lwordok = .true.
  return
!***************************************************************
!  Blocksize for parallel second derivative data distribution  *
!***************************************************************
820 continue
  if (nword.gt.1.and.nfloat.ge.1) then
    call stolc(words(2),maxword)
    if (index(words(2),'at').eq.1) then
      nblocksize = abs(nint(floats(1)))
      nblocksizevar = 3*nblocksize
    elseif (index(words(2),'sa').eq.1) then
      nblocksizesas = abs(nint(floats(1)))
    endif
  else
    call outerror('Blocksize input is incomplete',iline)
    call stopnow('genword')
  endif
  lwordok = .true.
  return
!*******************
!  UFF bond order  *
!*******************
825 continue
  if (nword.ge.2) then
    if (nfloat.lt.1) then
      call outerror('bond order is missing for uff_bondorder option',iline)
      call stopnow('genword')
    endif
    call stolc(words(2),maxword)
    if (index(words(2),'sin').eq.1) then
      call outerror('single bond order cannot be changed',iline)
      call stopnow('genword')
    elseif (index(words(2),'dou').eq.1) then
      UFFbondorder(2) = floats(1)
    elseif (index(words(2),'tri').eq.1) then
      UFFbondorder(3) = floats(1)
    elseif (index(words(2),'quad').eq.1) then
      UFFbondorder(4) = floats(1)
    elseif (index(words(2),'res').eq.1) then
      UFFbondorder(5) = floats(1)
    elseif (index(words(2),'ami').eq.1) then
      UFFbondorder(6) = floats(1)
    elseif (index(words(2),'cus').eq.1) then
      UFFbondorder(7) = floats(1)
    elseif (index(words(2),'hal').eq.1) then
      UFFbondorder(8) = floats(1)
    elseif (index(words(2),'quar').eq.1) then
      UFFbondorder(9) = floats(1)
    elseif (index(words(2),'thi').eq.1) then
      UFFbondorder(10) = floats(1)
    endif
  else
    call outerror('bond type missing for uff_bondorder option',iline)
    call stopnow('genword')
  endif
  lwordok = .true.
  return
!************
!  index_k  *
!************
830 continue
  if (nfloat.ge.1) then
    maxindk = nint(abs(floats(1)))
    if (maxindk.eq.0) then
      call outerror('index_k value cannot be specified to be zero',iline)
      call stopnow('genword')
    endif
!
!  Set dependent variables
!
    maxindk2 = 2*maxindk
    maxindk3 = maxindk2**2
  else
    call outerror('value is missing for index_k option',iline)
    call stopnow('genword')
  endif
  lwordok = .true.
  return
!********************************
!  Iterative charge parameters  *
!********************************
840 if (nfloat.gt.0) then
    nqitermax = abs(nint(floats(1)))
    if (nfloat.gt.1) then
      qitertol = abs((floats(2)))
    endif
  else
    read(nru,*,err=99,end=99) nqitermax
    iline = iline + 1
  endif
  lwordok = .true.
  return
!***************************************
!  Bond type default for species pair  *
!***************************************
850 nbondtype = nbondtype + 1
  if (nbondtype.gt.maxbondtype) then
    maxbondtype = nbondtype + 10
    call changemaxbondtype
  endif
!
!  Find number of core/shell words in list
!
  ncsword = 0
  nvarshift1 = 0
  nvarshift2 = 0
  nw = 0
  do while (nw.lt.nword) 
    nw = nw + 1
    word = words(nw)
    call stolc(word,maxwordlength)
    if (ncsword.eq.0) then
      if (index(word,'sh').eq.1) then
        nvarshift1 = nvarshift1 + maxele
        ncsword = ncsword + 1
      elseif (index(word,'bs').eq.1) then
        nvarshift1 = nvarshift1 + maxele
        ncsword = ncsword + 1
      elseif (index(word,'cor').eq.1) then
        ncsword = ncsword + 1
      elseif (index(word,'bco').eq.1) then
        ncsword = ncsword + 1
      endif
      if (ncsword.eq.1) then
        do nw2 = nw+1,nword
          words(nw2-1) = words(nw2)
        enddo
        nword = nword - 1
      endif
    elseif (ncsword.eq.1) then
      if (index(word,'sh').eq.1) then
        nvarshift2 = nvarshift2 + maxele
        ncsword = ncsword + 1
      elseif (index(word,'bs').eq.1) then
        nvarshift2 = nvarshift2 + maxele
        ncsword = ncsword + 1
      elseif (index(word,'cor').eq.1) then
        ncsword = ncsword + 1
      elseif (index(word,'bco').eq.1) then
        ncsword = ncsword + 1
      endif
      if (ncsword.eq.2) then
        do nw2 = nw+1,nword
          words(nw2-1) = words(nw2)
        enddo
        nword = nword - 1
      endif
    endif
  enddo
  ind2 = 0
  if (nword.ge.4) then
    call ltont(words(2),nvar1,itype1)
    nvar1 = nvar1 + nvarshift1
    call ltont(words(3),nvar2,itype2)
    nvar2 = nvar2 + nvarshift2
    ind = 4
    if (nword.ge.5) ind2 = 5
  else
    call outerror('Error in input for bondtype option',iline)
    call stopnow('genword')
  endif
  if (nvar1.eq.nvar2) then
    natbondtype(1,nbondtype) = nvar1
    natbondtype(2,nbondtype) = nvar2
    if (itype1.lt.itype2) then
      ntypbondtype(1,nbondtype) = itype2
      ntypbondtype(2,nbondtype) = itype1
    else
      ntypbondtype(1,nbondtype) = itype1
      ntypbondtype(2,nbondtype) = itype2
    endif
  elseif (nvar1.lt.nvar2) then
    natbondtype(1,nbondtype) = nvar2
    natbondtype(2,nbondtype) = nvar1
    ntypbondtype(1,nbondtype) = itype2
    ntypbondtype(2,nbondtype) = itype1
  else
    natbondtype(1,nbondtype) = nvar1
    natbondtype(2,nbondtype) = nvar2
    ntypbondtype(1,nbondtype) = itype1
    ntypbondtype(2,nbondtype) = itype2
  endif
  if (index(words(ind),'sin').eq.1) then
    nbondtypeptr(1,nbondtype) = 1
  elseif (index(words(ind),'dou').eq.1) then
    nbondtypeptr(1,nbondtype) = 2
  elseif (index(words(ind),'tri').eq.1) then
    nbondtypeptr(1,nbondtype) = 3
  elseif (index(words(ind),'quad').eq.1) then
    nbondtypeptr(1,nbondtype) = 4
  elseif (index(words(ind),'res').eq.1) then
    nbondtypeptr(1,nbondtype) = 5
  elseif (index(words(ind),'ami').eq.1) then
    nbondtypeptr(1,nbondtype) = 6
  elseif (index(words(ind),'cus').eq.1) then
    nbondtypeptr(1,nbondtype) = 7
  elseif (index(words(ind),'hal').eq.1) then
    nbondtypeptr(1,nbondtype) = 8
  elseif (index(words(ind),'quar').eq.1) then
    nbondtypeptr(1,nbondtype) = 9
  elseif (index(words(ind),'thi').eq.1) then
    nbondtypeptr(1,nbondtype) = 10
  else
    nbondtypeptr(1,nbondtype) = 0
  endif
  if (ind2.gt.0) then
    if (index(words(ind2),'cyc').eq.1) then
      nbondtypeptr(2,nbondtype) = 2
    elseif (index(words(ind2),'exo').eq.1) then
      nbondtypeptr(2,nbondtype) = 3
    else
      nbondtypeptr(2,nbondtype) = 1
    endif
  else
    nbondtypeptr(2,nbondtype) = 1
  endif
  lwordok = .true.
  return
!*************************************
!  Gasteiger parameters for species  *
!*************************************
860 continue
  if (nword.ge.2) then
    word = words(2)
  else
    call outerror('species missing for Gasteiger parameters',iline)
    call stopnow('genword')
  endif
!
!  Check whether input is symbol or option
!
  call worsy(word,lsymbol,.true.)
  if (.not.lsymbol) then
    l55 = .true.
    goto 868
  endif
!
  if (llibrary) then
    call okspec(lok1,word,ilp,.true.,ltype01)
    if (.not.lok1) then
      l55 = .true.
      goto 868
    endif
    if (ilp.gt.0.and.ilp.le.nspec) then
      inat = natspec(ilp)
      if (ltype01) then
        itype = 0
      else
        itype = ntypspec(ilp)
      endif
    elseif (ilp.eq.-1) then
      inat = maxele
      itype = 0
    endif
  else
    call ltont(word,inat,itype)
  endif
!
!  Find species
!
  lfound = .false.
  ns = 0
  do while (ns.lt.nspec.and..not.lfound)
    ns = ns + 1
    if (inat.eq.natspec(ns).and.itype.eq.ntypspec(ns)) then
      lfound = .true.
    endif
  enddo
  if (nfloat.lt.3) then
    call outerror('insufficient input reading Gasteiger parameters',iline)
    call stopnow('genword')
  endif
  if (lfound) then
    gastspec(1,ns) = floats(1)
    gastspec(2,ns) = floats(2)
    gastspec(3,ns) = floats(3)
    if (nfloat.ge.4) then
      gastspec(4,ns) = floats(4)
    else
      gastspec(4,ns) = 0.0_dp
    endif
    lgastinspec(ns) = .true.
    lgastinlibspec(ns) = llibrary
  else
    nspec = nspec + 1
    if (nspec.gt.maxspec) then
      maxspec = nspec + 20
      call changemaxspec
    endif
    linspec(nspec) = .true.
    natspec(nspec) = inat
    ntypspec(nspec) = itype
    massspec(nspec) = atmass(inat)
    lmassinspec(nspec) = .false.
    vdwspec(nspec) = rvdw(inat)
    lvdwinspec(nspec) = .false.
    gastspec(1,nspec) = floats(1)
    gastspec(2,nspec) = floats(2)
    gastspec(3,nspec) = floats(3)
    if (nfloat.ge.4) then
      gastspec(4,nspec) = floats(4)
    else
      gastspec(4,nspec) = 0.0_dp
    endif
    lgastinspec(nspec) = .true.
    lgastinlibspec(nspec) = llibrary
    ns = nspec
  endif
!     
!  Fitting flags
!       
  nfloatused = 4
  if (lfit.and.lfflags) then
    if (nfloat.ge.nfloatused+4) then
      n1 = int(floats(nfloatused+1))
      n2 = int(floats(nfloatused+2))
      n3 = int(floats(nfloatused+3))
      n4 = int(floats(nfloatused+4))
    else
      call outerror('insufficient fitting flags for Gasteiger parameters',iline)
      call stopnow('genword')
    endif
    if (n1.eq.1) then
      nfit = nfit + 1
      if (nfit.ge.maxfit) then
        maxfit = nfit + 10
        call changemaxfit
      endif
      nftyp(nfit) = 1
      nfpot(nfit) = 39
      nfvar(nfit) = ns
    endif
    if (n2.eq.1) then
      nfit = nfit + 1
      if (nfit.ge.maxfit) then
        maxfit = nfit + 10
        call changemaxfit
      endif
      nftyp(nfit) = 1
      nfpot(nfit) = 40
      nfvar(nfit) = ns
    endif
    if (n3.eq.1) then
      nfit = nfit + 1
      if (nfit.ge.maxfit) then
        maxfit = nfit + 10
        call changemaxfit
      endif
      nftyp(nfit) = 1
      nfpot(nfit) = 41
      nfvar(nfit) = ns
    endif
    if (n4.eq.1) then
      nfit = nfit + 1
      if (nfit.ge.maxfit) then
        maxfit = nfit + 10
        call changemaxfit
      endif
      nftyp(nfit) = 1
      nfpot(nfit) = 42
      nfvar(nfit) = ns
    endif
  endif
868 continue
  lwordok = .true.
  return
!*********************************
!  Gasteiger damping parameters  *
!*********************************
870 if (nfloat.gt.0) then
    gastdamp = floats(1)
  else
    read(nru,*,err=99,end=99) gastdamp
    iline = iline + 1
  endif
  if (nword.ge.2) then
    call stolc(words(2),maxword)
    if (index(words(2),'fix').eq.1) then
      ngastdamptype = 2
    elseif (index(words(2),'acc').eq.1) then
      ngastdamptype = 3
    endif
  endif
  gastdamp = min(abs(gastdamp),1.0_dp)
  lwordok = .true.
  return
!***********************************************
!  Thresholds for third order force constants  *
!***********************************************
880 if (nfloat.gt.0) then
    thresh_fc3_ind = abs(floats(1))
    if (nfloat.ge.2) then
      thresh_fc3_tot = abs(floats(2))
    endif
  endif
  lwordok = .true.
  return
!*****************************
!  Grid dimensions for SPME  *
!*****************************
890 continue
  if (nfloat.ge.3) then
    nqkgrid(1,ncurr) = abs(nint(floats(1)))
    nqkgrid(2,ncurr) = abs(nint(floats(2)))
    nqkgrid(3,ncurr) = abs(nint(floats(3)))
  else
    call outerror('insufficient grid dimensions for nqkgrid ',iline)
    call stopnow('genword')
  endif
  lwordok = .true.
  return
!****************************
!  B-spline order for SPME  *
!****************************
900 continue
  if (nfloat.ge.1) then
    nBsplineorder = abs(nint(floats(1)))
  else
    call outerror('B-spline order missing',iline)
    call stopnow('genword')
  endif
  if (nBsplineorder.lt.4) then
    call outerror('B-spline order less than 4 specified',iline)
    call stopnow('genword')
  endif
  lwordok = .true.
  return
!*****************************
!  ChemShell option setting  *
!*****************************
910 continue
  if (nword.gt.1) then
    call stolc(words(2),maxword)
    if (index(words(2),'init').eq.1) then
      ichemsh_output = 1
    elseif (index(words(2),'calc').eq.1) then
      ichemsh_output = 2
    else
      ichemsh_output = 0
    endif
  endif
  lwordok = .true.
  return
!*******************************
!  Selection of matrix format  *
!*******************************
920 if (nword.ge.3) then
    call stolc(words(2),maxword)
    call stolc(words(3),maxword)
    if (index(words(2),'hes').eq.1) then
      if (index(words(3),'tr').eq.1) then
        lhess2D = .false.
      elseif (index(words(3),'tw').eq.1) then
        lhess2D = .true.
      endif
!
!  Check that parallel choice makes sense
!
      if (nprocs.gt.1.and..not.lhess2D) then
        call outerror('triangular hessian cannot be used in parallel',iline)
        call stopnow('genword')
      endif
    elseif (index(words(2),'cos').eq.1) then
      if (index(words(3),'tr').eq.1) then
        lcosmo2D = .false.
      elseif (index(words(3),'tw').eq.1) then
        lcosmo2D = .true.
      endif
!
!  Check that parallel choice makes sense
!
      if (nprocs.gt.1.and..not.lcosmo2D) then
        call outerror('triangular cosmoA cannot be used in parallel',iline)
        call stopnow('genword')
      endif
    endif
  endif
  lwordok = .true.
  return
!**********************************************************
!  Change default value of imaginary frequency tolerance  *
!**********************************************************
930 if (nfloat.gt.0) then
    frqtol = abs(floats(1))
  else
    read(nru,*,err=99,end=99) frqtol
    frqtol = abs(frqtol)
    iline = iline + 1
  endif
  lwordok = .true.
  return
!*****************************************************************
!  Change default value of negative eigenvalue tolerance for RFO *
!*****************************************************************
940 if (nfloat.gt.0) then
    rfotoleig = abs(floats(1))
  else
    read(nru,*,err=99,end=99) rfotoleig
    rfotoleig = abs(rfotoleig)
    iline = iline + 1
  endif
  lwordok = .true.
  return
!****************************************************************
!  Change default value of projected gradient tolerance for RFO *
!****************************************************************
950 continue
  if (nfloat.gt.0) then
    rfotolgrad = abs(floats(1))
  else
    read(nru,*,err=99,end=99) rfotolgrad
    rfotolgrad = abs(rfotolgrad)
    iline = iline + 1
  endif
  if (rfotolgrad.gt.1.0_dp) then
    rfotolgrad = 10.0**(-rfotolgrad)
  endif
  lwordok = .true.
  return
!************************
!  Dielectric constant  *
!************************
960 continue
  if (nfloat.gt.0) then
    dielectriccfg(ncurr) = abs(floats(1))
  else
    read(nru,*,err=99,end=99) dielectriccfg(ncurr)
    dielectriccfg(ncurr) = abs(dielectriccfg(ncurr))
    iline = iline + 1
    if (dielectriccfg(ncurr).lt.1.0_dp) then
      call outerror('dielectric constant cannot be less than 1',iline)
      call stopnow('genword')
    endif
  endif
  lwordok = .true.
  return
!****************************************
!  Set options for trapping exceptions  *
!****************************************
970 if (nword.gt.1) then
    if (index(words(2),'fc ').eq.1) then
      ltrap_fc = .true.
      if (nfloat.gt.0) then
        trap_fc = floats(1)
      endif
    endif
  endif
  lwordok = .true.
  return
!*****************************
!  End of input for options  *
!*****************************
!
!  Error handling
!
99 call outerror('Error reading input data - check format near '//word,iline)
  call stopnow('genword')
1000 l1000 = .true.
  return
  end
