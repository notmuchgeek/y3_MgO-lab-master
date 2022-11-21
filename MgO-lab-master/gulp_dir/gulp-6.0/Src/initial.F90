  subroutine initial
!
!  Initialise arrays and variables at start of run
!
!  10/98 Fitting variable names are now initialise here
!   6/01 Fit labels for B6/B8 swapped to correct order
!  10/02 Torharm potential added
!  10/02 Lambda sums initialised
!   3/03 Fitlabel initialisation dimensions corrected
!  11/04 Initialisation of pi now done here
!  11/04 linspec removed since it is now a pointer
!   4/05 Setting up of default values for tweatpic & selfwolfc added
!   6/05 Modified to allow fitting of electronegativity parameters
!   7/05 ladiabatic added
!   7/05 Initialisation of Streitz and Mintmire parameters added
!  11/05 EAM potential shift added
!   3/06 Poly harmonic potential added
!   3/06 Fitlabel corrected for polynomial potentials
!   5/06 maxfnspec added
!   7/06 Sixbody potentials added
!   8/06 QEq radius added to fitlabel
!   8/06 Optimisation change criteria now initialised
!  10/06 Initialisation for neutron variables added ers29
!  11/06 NEB modifications added
!  12/06 Dummy line added to mark place of qoverr2 in fitlabel added
!   1/07 Gasteiger data added
!   1/07 lautobond option initialised
!   1/07 Force constant potential added
!   2/07 p4 and p5 added to EAM functional in fitlabel
!   3/07 sr_glue potential added
!   3/07 lpreserveQ added
!   4/07 UFF data added
!   5/07 More initialisations added
!   5/07 Probability of straining cell added
!   5/07 lx0centroid added
!   5/07 Strain target data added
!   5/07 UFFtor added
!   5/07 Number of bond charge increments (nbondQ) added
!   5/07 Bond charge increments added to fitlabel
!   6/07 nrotationtype added
!   7/07 reaxFFtol added
!   7/07 Plane potential added
!   8/07 Default for cutp increased to 100000.0 Ang
!  10/07 Angle-angle cross potential added
!  10/07 QM/MM initialisation call to resetcfgdefaults added
!  11/07 ReaxFFtol changed to match value in reaxff code
!  11/07 fitlabel array modified to reflect new EAM parameters
!  12/07 mdmaxtemp option added
!   1/08 lreaxFFqreal removed
!   3/07 mdmaxvolume option added
!   3/08 erferfc, erf and reperfc potentials added
!   3/08 fitlabel expanded for reaxFF parameters
!   4/08 Option to input spatial decomposition size added
!   4/08 reaxFFatol added
!   5/08 UFF oop potential added
!   5/08 New UFF generators for OOP added
!   5/08 Default number of fitting cycles set to 5000
!   6/08 ReaxFF Q shell variables added
!   7/08 pval6 added to reaxFFval3 array
!  10/08 ldumpnooverwrite and ndumpstep added
!  10/08 COSMO variables added
!  11/08 bacoscross form added
!  11/08 xcosangleangle potential added
!  11/08 torcosangle potential added
!  12/08 defaults for parameters in synchro module added
!  12/08 reaxFFgamma initialised to zero so that input values can be found
!   1/09 swap move added to Monte Carlo
!   1/09 Use of nfvar for two-body potentials modified
!   1/09 VBO_twobody potential added
!   2/09 Default value of accuracy changed to 12
!   3/08 3coulomb potential added
!   3/08 smallself added
!   3/08 lfinitediff added
!   5/09 findiffc and findiffs defaults added
!   6/09 Charge as a coordinate option added
!   6/09 Module name changed from three to m_three
!   7/09 exp2 potential added
!   7/09 exppowers potential added
!   8/09 Grimme_C6 potential added
!   9/09 Maximum order of polynomial potential increased to 8
!  10/09 ReaxFF qr12 term added
!   1/10 One-body potentials added
!   2/10 Central force model potentials added
!   3/10 Separate default weights for cell, coordinates and properties added
!   4/10 COSMO file added
!   5/10 gcoulomb potential added
!   8/10 Initialisation of p_iso and p_flx input flags added
!   8/10 lpisoinput/lpflxinput added
!   8/10 lfix1atom added
!   9/10 Neutron scattering modifications added
!   9/10 All scalar variables from modules now initialised here
!   9/10 Default MD integrator changed to stochastic
!  10/10 EDIP linear threebody modifications added
!  10/10 EDIP accuracy values set
!  10/10 qbo file output added
!  12/10 Hiding of shells added
!   1/11 Setting of alternate space group symbols moved here from settings
!   1/11 lforcemin flag added
!   2/11 lreaxFFQ added
!   3/11 lstressout added
!   8/11 reaxFFchi and reaxFFmu initialised to zero
!   8/11 Initialise of plumed variables added
!   8/11 Routine renamed to .F90 so that preprocessor options can be used
!   9/11 g3coulomb potential added to fitlabel
!   9/11 Madelung correction added
!   3/12 Output of DCD files added
!   4/12 lreaxFFunder added
!   5/12 sixth added 
!   6/12 nobsmode added
!   9/12 lpacha added
!  10/12 Support for OpenKIM added
!  10/12 Modifications for output of lammps files added
!  11/12 Modified Wolf sum of Fennell and Gezelter added
!  12/12 Initialisation of Q0 values added
!   2/13 BAGcross potential added
!   3/13 Initialisation of nreg1tot added
!   3/13 Error in labelling of terms in fitlabel for BAcross potential
!        corrected.
!   4/13 Default cuts changed back to old value of 0.8
!   5/13 Custom UFF bond order added
!   6/13 BALcross potential added
!   7/13 Extra initialisations of variables added to handle change in
!        setkeyword
!   7/13 lbroaden_scale initialised
!   8/13 Keyword to include imaginary phonons in output added
!   8/13 lfindsym added
!   9/13 Extra timing info added
!  10/13 MD iteration defaults changed
!  10/13 lsopt added
!  10/13 maxindk2 & maxindk3 initialised based on maxindk
!  10/13 lbroaden_gauss initialised
!  12/13 CASTEP phonon output added
!   1/14 Fitting label added for EDIP2p
!   2/14 Initialisation of KIM test descriptor wrapped
!   3/14 Output of ShengBTE files added
!   3/14 derfc changed to g_derfc for benefit of ChemShell
!   3/14 Fractional bond orders added
!   3/14 literativeQ logical added
!   4/14 Iterative charge parameters added
!   6/14 lnoexclude added
!   8/14 lgroupvelocity calculation added
!   8/14 fitlabel modified for Baskes a4 potential
!   8/14 lMEAMscreen removed as this is now a species specific array
!   8/14 nmeamrhotype added
!   8/14 nmeamcombotype added
!   9/14 nlibeam* added
!  10/14 fc_supercell variables added
!  12/14 Output of moment of inertia added
!   2/15 nlibnobo added
!   2/15 MM3 potentials added
!   3/15 MM3angle potential added
!   3/15 New fitlabel related arrays added
!   3/15 cif_dummylattice added
!   3/15 bond type variables added
!   3/15 Gasteiger species fitting labels added
!   3/15 gasttol changed from 1.0d-3 to 1.0d-6
!   3/15 Gasteiger damping options added
!   3/15 lnoqeem added
!   5/15 kpt file added
!   8/15 Thresholds added for third order force constants
!   8/15 Garofalini form of sw3 added
!   8/15 New lower options added
!   9/15 Fitting numbers changed for bond order parameters
!   9/15 Flag for input of stepmax added
!   9/15 Separate stepmax added for rfo
!   3/16 Default for lfastMPI changed to true
!   3/16 BO coordination parameters added to fitlabel
!   3/16 Changes to plumed file names
!   4/16 nmcswappair added
!   5/16 Multiple mcswaps added
!   6/16 SPME added
!   7/16 lshengBTE3rd added
!   8/16 IR intensity output added
!   9/16 eamXcutfactor added
!  10/16 Lennard-Jones with fixed coefficients added
!   1/17 lPrintSix added
!   2/17 llammpsvar added
!   3/17 loldvarorder added
!   4/17 Chemshell flag initialised
!   4/17 Initialisation of matrix format flags added
!   4/17 Blocksizes explicitly initialised
!   5/17 lnewdefalg initialised
!   5/17 Defect symmetry turned off for parallel case by default
!   7/17 nmaxcells increased
!   7/17 accf1D added
!   7/17 nblocksizevar forced to be 3 x nblocksize
!   7/17 ldumpgenerated added
!   8/17 lnumerical added
!   8/17 accf1D changed
!  10/17 Delta_dipole added
!  10/17 Absolute coordinate flag added
!  10/17 Dipole output file added
!  11/17 Buffered_LJ added
!  11/17 lmsd added
!  12/17 lfastfd added
!   1/18 lghost added
!   1/18 Default for stepmax increased
!   1/18 Trace added
!   2/18 lalamode added
!   2/18 Slater potential added
!   3/18 frqtol added
!   4/18 leembond added
!   4/18 XSF file added
!   5/18 nqrange added
!   5/18 Electronegativity fitting changed
!   6/18 lallbonds added
!   6/18 lstraincell added
!   6/18 switch_stepmx added
!   6/18 Coulomb subtract flag added to fit labels
!   7/18 Default for stepmax decreased back to 1
!   7/18 stpchcrit initialised
!   8/18 Thresholds for force constants decreased
!   8/18 Modified for version 2 of OpenKIM
!  10/18 RFO control variables added and stepmaxrfo renamed 
!  11/18 lowersign added
!  12/18 lusevectors added
!   2/19 ldump80 added
!   2/19 Rigid molecules added
!   3/19 Lone pair smoothing added to ReaxFF
!   3/19 Multiple temperature ramps added
!   4/19 lfinitestrain initialised
!   5/19 Finite difference flag split for first and second derivatives
!   6/19 j2 and j3 potentials added
!   7/19 lceigen added
!   8/19 Short range damping of polarisation added
!   8/19 lmolatom added
!   9/19 ala_mode added
!   8/19 Coulomb-related variables initialised to zero in case they are accessed even when not needed
!   9/19 lnotorXduplicates flag added
!   9/19 lnoshellzero added
!   9/19 Initialisation of reaxFF charge arrays moved to initmaxreaxffspecdefaults
!  10/19 Ala_mode added
!  10/19 Langevin damping of dipoles added
!  11/19 ppp3body added
!  12/19 Trap module added
!   3/20 SW2 potential generalised to include p, q and C
!   4/20 lmolcom_mass added
!   4/20 lnorotate added
!   5/20 lphasecom added
!   6/20 nmaxcells increased from 50 to 200
!   7/20 Default weight for quaternions added
!   7/20 lkcal and lkjmol added
!   8/20 sumsfac, sfac and sfac0 now initialised here
!   8/20 Correction to BO h fitting label
!   9/20 lpolarisation added
!  10/20 lterseinmol added
!  11/20 Tersoff reorganised
!   1/21 BO lambda and omega fitting variables swapped
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
!  Copyright Curtin University 2021
!
!  Julian Gale, CIC, Curtin University, January 2021
!
  use bondcharge,    only : nbondQ
  use bondorderdata
  use brennerdata,   only : lbrennersplinef, lbrennersplineh, nat2REBOspecies, nbrennertype
  use cellinputflag
  use cellmultipole
  use configurations
  use g_constants,   only : pi, sqrtpi, hbar_A2Kgpers
  use g_constants,   only : hbar_js, planck, cmtorads, thztorad, mevtorad, speedl, mevtot
  use control
  use cosmic,        only : lcosmic, isasatomoption, cosmorange, cosmormax, cosmormaxs
  use cosmic,        only : nspa, nspah, nppa, npts, nptsh, nallnearseg, ldodeca, lsegsmooth
  use costfunction
  use current
  use defects
  use derivatives,   only : nd2cells, nd2cell, lfcsupercell, lfinitestrain
  use dispersion
  use distances,     only : extracutoff, ndistancetotal, lStoreVectors
  use dump
  use eam
  use EDIPdata,      only : EDIPaccuracy1, EDIPaccuracy2
  use EDIPdata,      only : EDIPaccuracy1drmax, EDIPaccuracy2drmax
  use element
  use eemdata,       only : nqrange, maxeemtype, lmultiqrange
  use gulp_files
  use fitting
  use four
  use freeze
  use general
  use genetic
  use gulpchemsh,    only : ichemsh_output
  use kim_models,    only : lkim_model, nkimmodel
  use ksample
  use kspace,        only : maxindk, maxindk2, maxindk3, accf1D, eta, tweatpi, seta
  use library
  use m_alamode,     only : ala_mode
  use m_pdfneutron,  only : init_pdf
  use m_plumed
  use m_pr,          only : lpisoinput, lpflxinput
  use m_three
  use maths,         only : lhess2D, lcosmo2D
  use mdlogic,       only : ladiabatic, lmd
  use moldyn
  use molecule
  use montecarlo
  use g_neb
  use numbers
  use observables
  use one
  use optimisation
  use parallel,      only : lfastMPI, nprocs
  use parallel,      only : nblocksize, nblocksizesas, nblocksizevar
  use phonout
  use plane
  use polarise
  use potentialgrid
  use potentialpoints
  use potentialsites
  use randomnumbers
  use reaxFFdata
  use scan
  use scatterdata
  use shells
  use shellextrapolation
  use shifts
  use six,           only : lPrintSix
  use spatial
  use spatialbo
  use species
  use spme,          only : lspme, maxqkgrid, nBsplineorder
  use sutton
  use symmetry
  use synchro
  use terse
  use thresholds
  use times
  use trap
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  use two
  use uffdata
  use velocities
  use wolfcosmo
  use xcgc
  implicit none
!
!  Local variables
!
  integer(i4)           :: i
  integer(i4)           :: j
  integer(i4)           :: k
  real(dp)              :: g_cpu_time
  real(dp)              :: g_derfc
#ifdef TRACE
  call trace_in('initial')
#endif
!
!  QM/MM requirement 
!
  call resetcfgdefaults(1_i4)
!*************************
!  Initialise variables  *
!*************************
!
!  Logicals
!
  labscoany = .false.
  ladiabatic = .true.
  lalamode = .false.
  lallbonds = .false.
  lallowgt1 = .false.
  langle = .false.
  lanneal = .false.
  lanyMEAMscreen = .false.
  larc = .false.
  latomicstress = .false.
  laver = .false.
  lbio = .false.
  lbond = .false.
  lbrenner = .false.
  lbrennersplinef = .true.
  lbrennersplineh = .true.
  lbroad = .false.
  lbroaden_gauss = .false.
  lbroaden_scale = .false.
  lbulknoopt = .false.
  lc6 = .false.
  lc6one = .false.
  lcas = .false.
  lcelllasttime = .false.
  lcello = .false.
  lceigen = .false.
  lcif = .false.
  lcomp = .false.
  lconj = .false.
  lconp = .false.
  lconv = .false.
  lcosmo = .false.
  lcosmo2D = (nprocs.gt.1)
  lcosmofile = .false.
  lcosmic = .false.
  lcssr = .false.
  ldcd = .false.
  ldcellr = .false.
  ldcharge = .false.
  ldebug = .false.
  ldefect = .false.
  ldfp = .false.
  lddipole = .false.
  ldip = .false.
  ldipole = .false.
  ldist = .false.
  ldlv = .false.
  lDoElectrostatics = .true.
  lDoQDeriv1 = .false.
  lDoQDeriv2 = .false.
  ldump80 = .false.
  lnebdoublynudged = .true.
  ldrv = .false.
  ldsym = (nprocs.eq.1)
  ldumpcart = .false.
  ldumpconnectivity = .false.
  ldumpgenerated = .false.
  ldumpnooverwrite = .false.
  leckart = .false.
  lEDIP = .false.
  leem = .false.
  leembond = .false.
  lefg = .false.
  leig = .false.
  leigen = .false.
  leprint = .false.
  leregion = .false.
  lewald = .false.
  lextrapolateshells = .true.
  lfastfd = .false.
  lfastMPI = .true.
  lfbfgs = .false.
  lfdf = .false.
  lfindsym = .false.
  lfinitediff1 = .false.
  lfinitediff2 = .false.
  lfinitestrain = .false.
  lfirst = .true.
  lfit = .false.
  lfix1atom = .true.
  lfixtangent = .false.
  lflags = .false.
  lforcemin = .false.
  lfrc = .false.
  lfree = .false.
  lfreq = .false.
  lfreqout = .false.
  lfrq = .false.
  lfrqbin = .true.
  lga = .false.
  lgasteiger = .false.
  lghost = .false.
  lgrad = .false.
  lgroupvelocity = .false.
  lharmrelax = .false.
  lhess2D = (nprocs.gt.1)
  lhex = .false.
  lhideshells = .false.
  lhis = .false.
  lkpt = .false.
  linclude_imaginary = .false.
  linertia = .false.
  linten = .false.
  linputvelc = .false.
  linputvelxyz = .false.
  lirout = .false.
  literativeQ = .false.
  lkcal = .false.
  lkjmol = .false.
  lkfull = .false.
  llammps = .false.
  llammpspots = .false.
  llammpsvar = .false.
  linputvelc = .false.
  linputvelxyz = .false.
  lkim_model = .false.
  llbfgs = .false.
  llibdump = .true.
  llibsymdump = .false.
  llower = .false.
  lmadelung = .false.
  lmarv = .false.
  lmarv2 = .false.
  lmarvreg2 = .false.
  lmc = .false.
  lmclowestwrite = .false.
  lmcout = .false.
  lmd = .false.
  lmeanke = .false.
  lminch = .false.
  lminimage = .false.
  lmodco = .true.
  lmodeset = .false.
  lmol = .false.
  lmolatom = .false.
  lmolcom_mass = .true.
  lmolq = .false.
  lmolmec = .false.
  lmolfix = .false.
  lmolrestart = .false.
  lmolrigid = .false.
  lmovie = .false.
  lmdconstrain = .false.
  lmsd = .false.
  lmstpch = .false.
  lmultiqrange = .false.
  lneb = .false.
  lnebclimbingimage = .false.
  lnebdoublynudged = .true.
  lnebvaryspring = .false.
  lnewdefalg = .false.
  lnoautobond = .false.
  lnoanald1 = .false.
  lnoanald2 = .false.
  lnoanald3 = .false.
  lnoenergy = .false.
  lnoexclude = .false.
  lnoflags = .false.
  lnoqeem = .false.
  lnoreal = .false.
  lnorecip = .false.
  lnorotate = .false.
  lnoshellzero = .false.
  lnotorXduplicates = .false.
  lnomottlittleton = .false.
  lnudgegc = .false.
  lnumdiag = .false.
  lnumerical = .false.
  loldinten = .false.
  loldvarorder = .false.
  lopt = .false.
  loptcellpar = .false.
  loptlower = .false.
  losc = .false.
  loutshell = .false.
  lpacha = .false.
  lphasecom = .true.
  lphon = .false.
  lphono = .false.
  lplumed = .false.
#ifdef PLUMED
  lplumed_available = .true.
#else
  lplumed_available = .false.
#endif
  lpolar = .false.
  lpolarisation = .false.
  lpoldamp = .false.
  lpollangevin = .false.
  lposidef = .false.
  lpot = .false.
  lpredict = .false.
  lqbo = .false.
  lqbond = .false.
  lqboappend = .false.
  lqpolar = .false.
  lqtpie = .false.
  lpacha = .false.
  lpflxinput = .false.
  lpisoinput = .false.
  lpre = .false.
  lpreserveQ = .false.
  lprop = .false.
  lqeq = .false.
  lquicksearch = .true.
  lraman = .false.
  lrcspatial_anisotropic = .false.
  lrcspatialBO_anisotropic = .false.
  lreaxFF = .false.
  lreaxFFQ = .true.
  lreaxFFlpsmooth = .false.
  lreaxFFovsmooth = .false.
  lrelax = .false.
  lregionforce = .false.
  lrest = .false.
  lrigid = .false.
  lrfo = .false.
  lSandM = .false.
  lSandMnoZZ = .true.
  lsas = .false.
  lsave = .false.
  lshello = .false.
  lsiteenergy = .false.
  lsopt = .false.
  lspme = .false.
  lstraincell = .false.
  lstraincellprop = .false.
  lstepmaxin = .false.
  lstressout = .false.
  lthermal = .false.
!
!  Scatter
!
  lscatter = .false.
  lscattercall = .false.
  do i = 1,9
     lscatopt(i) = .false.
     lscatanopt(i) = .false.
  enddo
  lscatoptflag = .false.
  lscatanoptflag = .false.
  lscatsurfaceweight = .false.
!
!  End scatter
!
  lshengBTE3rd = .true.
  lshengBTE = .false.
  lspatial = .false.
  lspatialBOok = .false.
  lspatialok = .false.
  lstaticfirst = .false.
  lsuttonc = .false.
  lstr = .false.
  lsym = .true.
  lsymdok = .true.
  lsymoff = .false.
  lsymopt = .false.
  ltersederivs = .false.
  lterseincell = .false.
  lterseincoords = .false.
  lterseinmol = .false.
  lterseoutcell = .false.
  lterseoutcoords = .false.
  ltersepotentials = .false.
  ltethered = .false.
  ltors = .false.
  ltran = .false.
  ltrap_fc = .false.
  ltrjascii = .false.
  ltrjequil = .false.
  lthb = .false.
  ltrj = .false.
  lunit = .false.
  lusevectors = .false.
  lx0centroid = .false.
  lxr = .false.
  lxray = .false.
  lxsf = .false.
  lxsfseq = .false.
  lxsfsym = .false.
  lxtl = .false.
  lxyz = .false.
  lxyzmovie = .false.
  lwolf = .false.
  lwolffennell = .false.
  lwolforiginal = .false.
  lzsisa = .false.
!
!  Reals and integers
!
  pi = 4.0_dp*atan(1.0_dp)
  sqrtpi = sqrt(pi)
!
!  Derived constants for neutrons : ers29
!
  hbar_js = planck/(2.0_dp*pi) 
  hbar_A2Kgpers = (planck*1d20)/(2.0_dp*pi)
  cmtorads = speedl*2.0_dp*pi 
  thztorad = (2.0_dp*pi)*1d12
  mevtot = 0.241797_dp
  mevtorad = mevtot*thztorad 
!
!  Gasteiger data
!
  gasteigerA(1:maxele) = 0.0_dp
  gasteigerB(1:maxele) = 0.0_dp
  gasteigerC(1:maxele) = 0.0_dp
!  Hydrogen
  gasteigerA(1) =  7.17_dp
  gasteigerB(1) =  6.24_dp
  gasteigerC(1) = -0.56_dp
!  Carbon - sp3
  gasteigerA(6) =  7.98_dp
  gasteigerB(6) =  9.18_dp
  gasteigerC(6) =  1.88_dp
!  Nitrogen - sp3
  gasteigerA(7) = 11.54_dp
  gasteigerB(7) = 10.82_dp
  gasteigerC(7) =  1.36_dp
!  Oxygen - sp3
  gasteigerA(8) = 14.18_dp
  gasteigerB(8) = 12.92_dp
  gasteigerC(8) =  1.39_dp
!  Fluorine
  gasteigerA(9) = 14.66_dp
  gasteigerB(9) = 13.85_dp
  gasteigerC(9) =  2.31_dp
!  Phosphorous
  gasteigerA(15) =  7.40_dp
  gasteigerB(15) =  3.00_dp
  gasteigerC(15) = -1.00_dp
!  Sulphur
  gasteigerA(16) = 10.14_dp
  gasteigerB(16) =  9.13_dp
  gasteigerC(16) =  1.38_dp
!  Chlorine
  gasteigerA(17) = 11.00_dp
  gasteigerB(17) =  9.69_dp
  gasteigerC(17) =  1.35_dp
!  Bromine
  gasteigerA(35) = 10.08_dp
  gasteigerB(35) =  8.47_dp
  gasteigerC(35) =  1.16_dp
!  Iodine
  gasteigerA(53) =  9.90_dp
  gasteigerB(53) =  7.96_dp
  gasteigerC(53) =  0.96_dp
!
  do i = 1,18
    chi(i) = chiold(i)
    rmu(i) = rmuold(i)
  enddo
  do i = 19,maxele
    chi(i) = 0.0_dp
    rmu(i) = 0.0_dp
  enddo
  do i = 1,maxele
    atsymin(i) = atsym(i)
    atmassin(i) = atmass(i)
    bbarin(i) = bbar(i)
    leemparaltered(i) = .false.
    qeqq0(i) = 0.0_dp
    rcovin(i) = rcov(i)
    rionin(i) = rion(i)
    rvdwin(i) = rvdw(i)
    sigincin(i) = siginc(i)
    smchi(i) = 0.0_dp
    smmu(i) = 0.0_dp
    smq0(i) = 0.0_dp
    smzeta(i) = 0.0_dp
    smZnuc(i) = 0.0_dp
    q0(i) = 0.0_dp
    q0_pacha(i) = 0.0_dp
  enddo
  accf1D = 8.0_dp
  accuracy = 12.0_dp
  ala_mode = 'all'
  bfactor = 0.2_dp
  bpdamp = 0.0_dp
  rpdamp = 0.0_dp
  Lor_tol = 1.0d-2
!
  cellmin = 0.5_dp
  cif_dummylattice = 0.0_dp
  chcrit = 0.0_dp
  chempot = 0.0_dp
  cutb = 2.0_dp
  cutp = 100000.0_dp
  cuts = 0.8_dp
  cutw = 0.0_dp
  delfc = 10.0_dp
  delta = 1.0d-4
  delwht = 1.0_dp
  delwht_angle = 1.0_dp
  delwht_bond = 1.0_dp
  delwht_cell_length = 1000.0_dp
  delwht_cell_angle = 1000.0_dp
  delwht_coord = 1000.0_dp
  delwht_dielectric = 1.0_dp
  delwht_dipole = 10.0_dp
  delwht_elastic = 0.01_dp
  delwht_energy = 1.0_dp
  delwht_frac = 10000.0_dp
  delwht_freq = 1.0_dp
  delwht_grad = 1.0_dp
  delwht_modulus = 1.0_dp
  delwht_qtn = 1.0_dp
  delwht_radius = 1000.0_dp
  delwht_strain = 1.0_dp
  delwht_stress = 1.0_dp
  dmaxmc = 0.05_dp
  eamXcutfactor = 0.5_dp
  eta = 0.0_dp
  etaw = 0.0_dp
  extracutoff = 0.0_dp
  fbox = 0.0_dp
  flbox = 0.0_dp
  ftol = 0.00001_dp
  fftol = 0.00001_dp
  fgmax = 0.001_dp
  fgtol = 0.0001_dp
  findiff = 0.0_dp
  findiffc = 0.00001_dp
  findiffs = 0.00001_dp
  frqtol = 0.5_dp
  fxtol = 0.00001_dp
  fstepmx = 1000.0_dp
  ngastdamptype = 1_i4
  gastdamp = 0.5_dp
  gasttol = 0.000001_dp
  gdcrit = 4.0_dp
  grmax = 0.01_dp
  gtol = 0.001_dp
  lmgtol = 0.9_dp
  lowerscale = 0.001_dp
  lowersign = 1.0_dp
  mcelowest = 0.0_dp
  mcemean = 0.0_dp
  mcnmean = 0.0_dp
  mcvolume = 0.0_dp
  pcreate = 0.0_dp
  pdestroy = 0.0_dp
  phondiff = 0.00001_dp
  pmove = 1.0_dp
  protate = 0.0_dp
  pstrain = 0.0_dp
  prob(1) = 0.8_dp
  prob(3) = 2.0_dp
  prob(4) = 0.4_dp
  prob(6) = 2.0_dp
  prob(7) = -1.0_dp
  prob(9) = 2.0_dp
  pswapsum = 0.0_dp
  qeqlambda = 0.5_dp
  qeqscfcrit = 1.0d-6
  qitertol = 1.0d-8
  rbox = 6.0_dp
  rcspatial = 0.0_dp
  rcspatialx = 0.0_dp
  rcspatialy = 0.0_dp
  rcspatialz = 0.0_dp
  rcspatialbo = 0.0_dp
  rcspatialbox = 0.0_dp
  rcspatialboy = 0.0_dp
  rcspatialboz = 0.0_dp
  reaxFFcutoff = 10.0_dp
  reaxFFcutoffQ   = 10.0_dp
  reaxFFcutoffVDW = 10.0_dp
  reaxFFtol = 1.0d-3
  reaxFFatol = 1.0d-3
  reaxFFatol2 = reaxFFatol**2
  reaxFFatol3 = reaxFFatol**3
  reaxFFhtol = 0.01_dp
  reaxFFlpsmooth = 80.0_dp
  reaxFFrhtol = 7.5_dp
  reaxFFtaperscale = 2.0_dp
  reaxFFqconverged = 1.0d-6
  reaxFFqconverged1 = 1.0d-5
  reaxFFqdamp = 0.15_dp
  nreaxFFqiter = 100
!
!  Scatter
!
  phi_initial = - 0.5_dp*pi
  theta_initial = - 0.5_dp*pi
!
  rfostepmax = 100.0_dp
  rfotoleig = 0.001_dp
  rfotolgrad = 1.0d-8
  rkw = 0.0_dp
  rmaxmc = 1.0_dp
  rmdmaxtemp = 100.0_dp
  rmdmaxvol = 100.0_dp
  rpmax = 0.0_dp
  rqeq = 15.0_dp
  rspeed0 = 1.0_dp
  rtol = 1.2_dp
  scalefactor = 1.0_dp
  scmaxsearch = 2.0_dp
  selfwolf = 0.0_dp
  selfwolfenergy = 0.0_dp
  selfwolfforce = 0.0_dp
  seta = 0.0_dp
  sfac = 0.0_dp
  sfac0 = 0.0_dp
  sgtol = 1.0d-8
  sixth = 1.0_dp/6.0_dp
  smallself = 1.0d-12
  smaxmc = 0.1_dp
  stepmax = 1.0_dp
  stepmin = 1.0d-12
  stpchcrit = 0.0_dp
  sumlambdaR = 0.0_dp
  sumlambdaV = 0.0_dp
  sumsfac = 0.0_dp
  tapertype = 1
  tapermax = 0.0_dp
  tapermin = 0.0_dp
  taperm = 20.0_dp
  targetmove = 0.0_dp
  targetrota = 0.0_dp
  targetstrain = 0.0_dp
  targetrradmax = 0.0_dp
  third = 1.0_dp/3.0_dp
  thresh_fc3_ind = 0.0000001_dp
  thresh_fc3_tot = 0.00001_dp
  timesofar = 0.0_dp
  trap_fc = 1000.0_dp
  tweatpi = 0.0_dp
  xtol = 0.00001_dp
  velmax = 100.0_dp
  ichemsh_output = 0
  icmm = 0
  idump = 0
  iseed = -1
  lmbfgsorder = 10
  maxgacyc = 300
  maxindk2 = 2*maxindk
  maxindk3 = maxindk2**2
  maxsynciter = 1000
  maxsyncstep = 1000
  mgacfg = 10
  maxcal = 1000
  maxextrapol = 8
  maxfcal = 5000
  maxqkgrid(1:3) = 2
  mode = 0
  mode2a = 4
  moptit = 99
  morder = 0
  nadd = 0
  nasum = 0
  natab = 0
  nblocksize = 1
  nblocksizesas = 12
  nblocksizevar = 3*nblocksize
  nBsplineorder = 8
  nbufferx = 1
  nbuffery = 1
  nbufferz = 1
  nbufferxbo = 1
  nbufferybo = 1
  nbufferzbo = 1
  nbrennertype = 3
  nboA = 0
  nboO = 0
  nboR = 0
  nboQ = 0
  nboQ0 = 0
  nboZ = 0
  nbopot = 0
  nlibeamfnspec = 0
  nlibeamspec = 0
  nlibnboA = 0
  nlibnboO = 0
  nlibnboR = 0
  nlibnboQ = 0
  nlibnboQ0 = 0
  nlibnbopot = 0
  nlibnobo = 0
  nlibbondtype = 0
  nbondQ = 0
  ncell = 0
  ncellsearch = 1
  ncellsearchbo = 1
  ncfg = 0
  ncfmdrun = 0
  nconnect = 0
  ncontot = 0
  ncycd = 1000
  ncycp = 1000
  nd = 6
  ndi = -1
  ndmax = 6
  ndcon = 0
  ndef = 0
  ndispres = 30
  ndline = 0
  ndpoint = 0
  ndumpstep = 0
  neamfn = 1
  neampower = 2
  neamfnspec = 0
  neamspec = 0
  nebmaxdisp = 0.1_dp
  nebrandom = 0.0_dp
  nebtangent = 3
  nebtol = 1.0d-4
  nnebreplicatot = 0
  nfcon = 0
  nfgrad = 0
  nfstrain = 0
  nfit = 0
  ngabest = 2
  ngabset = 0
  ngacfg = 10
  ngacjg = 0
  ngcmcmol = 0
  ngcmcspec = 0
  nkimmodel = 0
  nmcaccepted = 0
  nmcoutfreq = 100
  nmcsample = 10
  nmcswaps = 0
  nmcstep = 1
  nmctrial = 0
  nmeamcombotype = 2
  nmeamrhotype = 2
  nqrange(1:maxele,1:maxeemtype) = 0
  ntargetfreq = 0
  ntargetfreqr = 0
  ntargetfreqs = 0
  ntrialatom = 0
  nfor = 0
  nkpt = 0
  nbox = 0
  nlbox = 0
  nlib = 0
  nlinmin = 8
  ninte = 0
  nvaca = 0
  ntotinte = 0
  ntotvaca = 0
  ndistancetotal = 1
  nmaxcells = 200
  nemorder = 4
  nmdintegrator = 4
  nmditer = 3
  nnebiter = 1000
  nbondtype = 0
  nnobo = 0
  nobs = 0
  nobsmode = 0
  none = 0
  nplanepot = 0
  npolspec = 0
  npote = 0
  npotpt = 0
  npotsites = 0
  ngastitermax = 50
  nqeqitermax = 20
  nqitermax = 200
  nratiomspec = 0
  nrotationtype = 0
  nreg1tot = 0
  nseps = 0
  nshift = 0
  nspar = 0
  nspec = 0
  ntitle = 0
  nthb = 0
  nupdate = 10
  nfupdate = 20
  nUFFspec = 0
  nwarn = 0
!
!  fc_supercell variables
!
  lfcsupercell = .false.
  nd2cells = 0
  nd2cell(1:3) = 0
!
!  Scatter
!
  nq_step_fit = 0
  nw_step_fit = 0
!
  synctol = 1.0d-4
  time0 = g_cpu_time()
  timmax = -1.0_dp
!
!  COSMO variables
!
  isasatomoption = 1
  nallnearseg = 0
  nppa = 974
  npts = 0
  nptsh = 0
  nspa = 110
  nspah = 110
  cosmorange = 0.0_dp
  cosmormax = 10.0_dp
  cosmormaxs = 1.0_dp
  ldodeca = .false.
  lsegsmooth = .true.
!
!  Cost function variables
!
  kacf = 1.0_dp
  kbcf = 1.0_dp
  kcccf = 0.1_dp
  kcacf = 0.2_dp
  kqccf = 1.0_dp
  kqacf = 1.0_dp
  kscf = 0.0_dp
!
!  EAM variables
!
  lPrintEAM  = .false.
  lMEAM    = .false.
  lMEAMfn  = .false.
  lMEAMden = .false.
!
!  Genetic variables
!
  ngabset = 0
  l2pxo = .false.
  lgabest = .false.
  lgadef = .false.
  lgaexpw = .false.
  loxygen = .false.
  lstruc = .false.
  udif = 0.0_dp
!
!  Random numbers
!
  nrandomcalls = 0
  npr_randomcalls = 0
  npr_grandomcalls = 0
  npr_randomcalls_adv = 0
  npr_grandomcalls_adv = 0
  lGaussianLast = .false.
!
!  ReaxFF
!
  nreaxFFspec = 0
  nlibreaxFFspec = 0
  nreaxFFqiter = 100
!
!  Wolf cosmo
!
  lPureCoulomb0D = .false.
  cutwc = 20.0_dp
  etawc = 0.05_dp
!
!  Optimisation flags
!
  lPrintVar = .false.
!
!  EDIP parameters
!
  EDIPaccuracy1 = 0.000001_dp
  EDIPaccuracy2 = 0.0000000001_dp
  EDIPaccuracy1drmax = 1.0_dp/log(EDIPaccuracy1)
  EDIPaccuracy2drmax = 1.0_dp/log(EDIPaccuracy2)
!
!  Odd variables
!
  narcwrite = 1
  nfreqdecimals = 6
  lmdarcunique = .false.
  llibsymdump = .false.
  lPrintFour = .false.
  lPrintSix = .false.
  lPrintThree = .false.
  lPrintTwo = .false.
  lStoreVectors = .false.
!
!  Times
!
  tatom = 0.0_dp
  tbondorder = 0.0_dp
  tbrenner = 0.0_dp
  tcorrect = 0.0_dp
  tcosmo = 0.0_dp
  tcosmoderv = 0.0_dp
  tderv3 = 0.0_dp
  tdiag = 0.0_dp
  tdisk = 0.0_dp
  tedip = 0.0_dp
  teem = 0.0_dp
  tfederiv = 0.0_dp
  tfitf = 0.0_dp
  tfun = 0.0_dp
  tfour = 0.0_dp
  thes = 0.0_dp
  tion = 0.0_dp
  tkim = 0.0_dp
  tmany = 0.0_dp
  tmati = 0.0_dp
  tmc   = 0.0_dp
  tmdinit = 0.0_dp
  tphon = 0.0_dp
  tpolar = 0.0_dp
  tpredict = 0.0_dp
  tproj = 0.0_dp
  tprop = 0.0_dp
  treaxFF = 0.0_dp
  treg1 = 0.0_dp
  treg2a = 0.0_dp
  treg2b = 0.0_dp
  treg3 = 0.0_dp
  treg4 = 0.0_dp
  tregm = 0.0_dp
  tres = 0.0_dp
  trls = 0.0_dp
  tsearch = 0.0_dp
  tsix = 0.0_dp
  tspline = 0.0_dp
  tsum = 0.0_dp
  tsym = 0.0_dp
  tthree = 0.0_dp
  ttmat = 0.0_dp
  tvelcor = 0.0_dp
!
!  Initialise neutron variables
!
  call init_pdf
  call init_ins
!
!  Set some dependent values
!
  tweatpic = 2.0_dp*etawc/sqrtpi
  if (abs(cutwc).gt.1.0d-10) then
    selfwolfc = g_derfc(etawc*cutwc)/cutwc
  endif
!  
!  Set array of pointers from atomic number to REBO species
!  
!  Default = 0 => not parameterised yet
!  
  nat2REBOspecies(1:maxele) = 0
  nat2REBOspecies(6) = 1
  nat2REBOspecies(1) = 2
  nat2REBOspecies(8) = 3
  nat2REBOspecies(9) = 4
!***********************
!  Initialise strings  *
!***********************
  arcfile = ' '
  biofile = ' '
  casfile = ' '
  ciffile = ' '
  cosmofile = ' '
  cssrfile = ' '
  dcdfile = ' '
  dfile = ' '
  dipfile = ' '
  dlvfile = ' '
  drvfile = ' '
  eigfile = ' '
  fdffile = ' '
  frcfile = ' '
  freqfile = ' '
  hisfile = ' '
  inertfile = ' '
  kptfile = ' '
  lammpsfile = ' '
  lammpspotsfile = ' '
  marvfile = ' '
  marvtemp = ' '
  mcfile = ' '
  mdafil = ' '
  oscfile = ' '
  phonfile = ' '
  plumedinput = 'plumed.dat'
  plumedlog = 'plumed.log'
  prefile = ' '
  qbofile = ' '
  sasfile = ' '
  shengfile = ' '
  thbfile = ' '
  trjfile = ' '
  xrfile = ' '
  xsffile = ' '
  xtlfile = ' '
  xyzfile = ' '
  elefile = ' '
  helpfile = ' '
  keyword = ' '
  site = ' '
!
!  Neutron related strings : ers29
!
  pdffile = ' '
!*****************************
!  Alternative space groups  *
!*****************************
  do i = 1,74
    naltgnam(i) = 1
  enddo
  do i = 1,74
    do j = 1,18
      altgnam(j,i) = ' '
    enddo
  enddo
!
!  Numbers of alterative group names
!
  naltgnam(3:4) = 3
  naltgnam(5) = 9
  naltgnam(6) = 3
  naltgnam(7:8) = 9
  naltgnam(9) = 18
  naltgnam(10:11) = 3
  naltgnam(12:14) = 9
  naltgnam(15) = 18
  naltgnam(17:18) = 3
  naltgnam(20:21) = 3
  naltgnam(25) = 3
  naltgnam(26) = 6
  naltgnam(27) = 3
  naltgnam(28:31) = 6
  naltgnam(32) = 3
  naltgnam(33) = 6
  naltgnam(34:35) = 3
  naltgnam(36) = 6
  naltgnam(37) = 3
  naltgnam(38:41) = 6
  naltgnam(42:45) = 3
  naltgnam(46) = 6
  naltgnam(49:50) = 3
  naltgnam(51:54) = 6
  naltgnam(55:56) = 3
  naltgnam(57) = 6
  naltgnam(58:59) = 3
  naltgnam(60) = 6
  naltgnam(61) = 2
  naltgnam(62:64) = 6
  naltgnam(65:66) = 3
  naltgnam(67:68) = 6
  naltgnam(72) = 3
  naltgnam(73) = 2
  naltgnam(74) = 6
!
!  Alternative group names
!
  altgnam(1,1)   = 'P 1'
  altgnam(1,2)   = 'P -1'
  altgnam(1,3)   = 'P 1 2 1'
  altgnam(2,3)   = 'P 1 1 2'
  altgnam(3,3)   = 'P 2 1 1'
  altgnam(1,4)   = 'P 1 21 1'
  altgnam(2,4)   = 'P 1 1 21'
  altgnam(3,4)   = 'P 21 1 1'
  altgnam(1,5)   = 'C 1 2 1'
  altgnam(2,5)   = 'A 1 2 1'
  altgnam(3,5)   = 'I 1 2 1'
  altgnam(4,5)   = 'A 1 1 2'
  altgnam(5,5)   = 'B 1 1 2'
  altgnam(6,5)   = 'I 1 1 2'
  altgnam(7,5)   = 'B 2 1 1'
  altgnam(8,5)   = 'C 2 1 1'
  altgnam(9,5)   = 'I 2 1 1'
  altgnam(1,6)   = 'P 1 M 1'
  altgnam(2,6)   = 'P 1 1 M'
  altgnam(3,6)   = 'P M 1 1'
  altgnam(1,7)   = 'P 1 C 1'
  altgnam(2,7)   = 'P 1 N 1'
  altgnam(3,7)   = 'P 1 A 1'
  altgnam(4,7)   = 'P 1 1 A'
  altgnam(5,7)   = 'P 1 1 N'
  altgnam(6,7)   = 'P 1 1 B'
  altgnam(7,7)   = 'P B 1 1'
  altgnam(8,7)   = 'P N 1 1'
  altgnam(9,7)   = 'P C 1 1'
  altgnam(1,8)   = 'C 1 M 1'
  altgnam(2,8)   = 'A 1 M 1'
  altgnam(3,8)   = 'I 1 M 1'
  altgnam(4,8)   = 'A 1 1 M'
  altgnam(5,8)   = 'B 1 1 M'
  altgnam(6,8)   = 'I 1 1 M'
  altgnam(7,8)   = 'B M 1 1'
  altgnam(8,8)   = 'C M 1 1'
  altgnam(9,8)   = 'I M 1 1'
  altgnam(1,9)   = 'C 1 C 1'
  altgnam(2,9)   = 'A 1 N 1'
  altgnam(3,9)   = 'I 1 A 1'
  altgnam(4,9)   = 'A 1 A 1'
  altgnam(5,9)   = 'C 1 N 1'
  altgnam(6,9)   = 'I 1 C 1'
  altgnam(7,9)   = 'A 1 1 A'
  altgnam(8,9)   = 'B 1 1 N'
  altgnam(9,9)   = 'I 1 1 B'
  altgnam(10,9)  = 'B 1 1 B'
  altgnam(11,9)  = 'A 1 1 N'
  altgnam(12,9)  = 'I 1 1 A'
  altgnam(13,9)  = 'B B 1 1'
  altgnam(14,9)  = 'C N 1 1'
  altgnam(15,9)  = 'I C 1 1'
  altgnam(16,9)  = 'C C 1 1'
  altgnam(17,9)  = 'B N 1 1'
  altgnam(18,9)  = 'I B 1 1'
  altgnam(1,10)  = 'P 1 2/M 1'
  altgnam(2,10)  = 'P 1 1 2/M'
  altgnam(3,10)  = 'P 2/M 1 1'
  altgnam(1,11)  = 'P 1 21/M 1'
  altgnam(2,11)  = 'P 1 1 21/M'
  altgnam(3,11)  = 'P 21/M 1 1'
  altgnam(1,12)  = 'C 1 2/M 1'
  altgnam(2,12)  = 'A 1 2/M 1'
  altgnam(3,12)  = 'I 1 2/M 1'
  altgnam(4,12)  = 'A 1 1 2/M'
  altgnam(5,12)  = 'B 1 1 2/M'
  altgnam(6,12)  = 'I 1 1 2/M'
  altgnam(7,12)  = 'B 2/M 1 1'
  altgnam(8,12)  = 'C 2/M 1 1'
  altgnam(9,12)  = 'I 2/M 1 1'
  altgnam(1,13)  = 'P 1 2/C 1'
  altgnam(2,13)  = 'P 1 2/N 1'
  altgnam(3,13)  = 'P 1 2/A 1'
  altgnam(4,13)  = 'P 1 1 2/A'
  altgnam(5,13)  = 'P 1 1 2/N'
  altgnam(6,13)  = 'P 1 1 2/B'
  altgnam(7,13)  = 'P 2/B 1 1'
  altgnam(8,13)  = 'P 2/N 1 1'
  altgnam(9,13)  = 'P 2/C 1 1'
  altgnam(1,14)  = 'P 1 21/C 1'
  altgnam(2,14)  = 'P 1 21/N 1'
  altgnam(3,14)  = 'P 1 21/A 1'
  altgnam(4,14)  = 'P 1 1 21/A'
  altgnam(5,14)  = 'P 1 1 21/N'
  altgnam(6,14)  = 'P 1 1 21/B'
  altgnam(7,14)  = 'P 21/B 1 1'
  altgnam(8,14)  = 'P 21/N 1 1'
  altgnam(9,14)  = 'P 21/C 1 1'
  altgnam(1,15)  = 'C 1 2/C 1'
  altgnam(2,15)  = 'A 1 2/N 1'
  altgnam(3,15)  = 'I 1 2/A 1'
  altgnam(4,15)  = 'A 1 2/A 1'
  altgnam(5,15)  = 'C 1 2/N 1'
  altgnam(6,15)  = 'I 1 2/C 1'
  altgnam(7,15)  = 'A 1 1 2/A'
  altgnam(8,15)  = 'B 1 1 2/N'
  altgnam(9,15)  = 'I 1 1 2/B'
  altgnam(10,15) = 'B 1 1 2/B'
  altgnam(11,15) = 'A 1 1 2/N'
  altgnam(12,15) = 'I 1 1 2/A'
  altgnam(13,15) = 'B 2/B 1 1'
  altgnam(14,15) = 'C 2/N 1 1'
  altgnam(15,15) = 'I 2/C 1 1'
  altgnam(16,15) = 'C 2/C 1 1'
  altgnam(17,15) = 'B 2/N 1 1'
  altgnam(18,15) = 'I 2/B 1 1'
  altgnam(1,16)  = 'P 2 2 2'
  altgnam(1,17)  = 'P 2 2 21'
  altgnam(2,17)  = 'P 21 2 2'
  altgnam(3,17)  = 'P 2 21 2'
  altgnam(1,18)  = 'P 21 21 2'
  altgnam(2,18)  = 'P 2 21 21'
  altgnam(3,18)  = 'P 21 2 21'
  altgnam(1,19)  = 'P 21 21 21'
  altgnam(1,20)  = 'C 2 2 21'
  altgnam(2,20)  = 'A 21 2 2'
  altgnam(3,20)  = 'B 2 21 2'
  altgnam(1,21)  = 'C 2 2 2'
  altgnam(2,21)  = 'A 2 2 2'
  altgnam(3,21)  = 'B 2 2 2'
  altgnam(1,22)  = 'F 2 2 2'
  altgnam(1,23)  = 'I 2 2 2'
  altgnam(1,24)  = 'I 21 21 21'
  altgnam(1,25)  = 'P M M 2'
  altgnam(2,25)  = 'P 2 M M'
  altgnam(3,25)  = 'P M 2 M'
  altgnam(1,26)  = 'P M C 21'
  altgnam(2,26)  = 'P C M 21'
  altgnam(3,26)  = 'P 21 M A'
  altgnam(4,26)  = 'P 21 A M'
  altgnam(5,26)  = 'P B 21 M'
  altgnam(6,26)  = 'P M 21 B'
  altgnam(1,27)  = 'P C C 2'
  altgnam(2,27)  = 'P 2 C C'
  altgnam(3,27)  = 'P C 2 C'
  altgnam(1,28)  = 'P M A 2'
  altgnam(2,28)  = 'P B M 2'
  altgnam(3,28)  = 'P 2 M B'
  altgnam(4,28)  = 'P 2 C M'
  altgnam(5,28)  = 'P C 2 M'
  altgnam(6,28)  = 'P M 2 A'
  altgnam(1,29)  = 'P C A 21'
  altgnam(2,29)  = 'P B C 21'
  altgnam(3,29)  = 'P 21 A B'
  altgnam(4,29)  = 'P 21 C A'
  altgnam(5,29)  = 'P C 21 B'
  altgnam(6,29)  = 'P B 21 A'
  altgnam(1,30)  = 'P N C 2'
  altgnam(2,30)  = 'P C N 2'
  altgnam(3,30)  = 'P 2 A N'
  altgnam(4,30)  = 'P 2 N A'
  altgnam(5,30)  = 'P B 2 N'
  altgnam(6,30)  = 'P N 2 B'
  altgnam(1,31)  = 'P M N 21'
  altgnam(2,31)  = 'P N M 21'
  altgnam(3,31)  = 'P 21 M N'
  altgnam(4,31)  = 'P 21 N M'
  altgnam(5,31)  = 'P N 21 M'
  altgnam(6,31)  = 'P M 21 N'
  altgnam(1,32)  = 'P B A 2'
  altgnam(2,32)  = 'P 2 C B'
  altgnam(3,32)  = 'P C 2 A'
  altgnam(1,33)  = 'P N A 21'
  altgnam(2,33)  = 'P B N 21'
  altgnam(3,33)  = 'P 21 N B'
  altgnam(4,33)  = 'P 21 C N'
  altgnam(5,33)  = 'P C 21 N'
  altgnam(6,33)  = 'P N 21 A'
  altgnam(1,34)  = 'P N N 2'
  altgnam(2,34)  = 'P 2 N N'
  altgnam(3,34)  = 'P N 2 N'
  altgnam(1,35)  = 'C M M 2'
  altgnam(2,35)  = 'C 2 M M'
  altgnam(3,35)  = 'C M 2 M'
  altgnam(1,36)  = 'C M C 21'
  altgnam(2,36)  = 'C C M 21'
  altgnam(3,36)  = 'A 21 M A'
  altgnam(4,36)  = 'A 21 A M'
  altgnam(5,36)  = 'B B 21 M'
  altgnam(6,36)  = 'B M 21 B'
  altgnam(1,37)  = 'C C C 2'
  altgnam(2,37)  = 'A 2 A A'
  altgnam(3,37)  = 'B B 2 B'
  altgnam(1,38)  = 'A M M 2'
  altgnam(2,38)  = 'B M M 2'
  altgnam(3,38)  = 'B 2 M M'
  altgnam(4,38)  = 'C 2 M M'
  altgnam(5,38)  = 'C M 2 M'
  altgnam(6,38)  = 'A M 2 M'
  altgnam(1,39)  = 'A B M 2'
  altgnam(2,39)  = 'B M A 2'
  altgnam(3,39)  = 'B 2 C M'
  altgnam(4,39)  = 'C 2 M B'
  altgnam(5,39)  = 'C M 2 A'
  altgnam(6,39)  = 'A C 2 M'
  altgnam(1,40)  = 'A M A 2'
  altgnam(2,40)  = 'B B M 2'
  altgnam(3,40)  = 'B 2 M B'
  altgnam(4,40)  = 'C 2 C M'
  altgnam(5,40)  = 'C C 2 M'
  altgnam(6,40)  = 'A M 2 A'
  altgnam(1,41)  = 'A B A 2'
  altgnam(2,41)  = 'B B A 2'
  altgnam(3,41)  = 'B 2 C B'
  altgnam(4,41)  = 'C 2 C B'
  altgnam(5,41)  = 'C C 2 A'
  altgnam(6,41)  = 'A C 2 A'
  altgnam(1,42)  = 'F M M 2'
  altgnam(2,42)  = 'F 2 M M'
  altgnam(3,42)  = 'F M 2 M'
  altgnam(1,43)  = 'F D D 2'
  altgnam(2,43)  = 'F 2 D D'
  altgnam(3,43)  = 'F D 2 D'
  altgnam(1,44)  = 'I M M 2'
  altgnam(2,44)  = 'I 2 M M'
  altgnam(3,44)  = 'I M 2 M'
  altgnam(1,45)  = 'I B A 2'
  altgnam(2,45)  = 'I 2 C B'
  altgnam(3,45)  = 'I C 2 A'
  altgnam(1,46)  = 'I M A 2'
  altgnam(2,46)  = 'I B M 2'
  altgnam(3,46)  = 'I 2 M B'
  altgnam(4,46)  = 'I 2 C M'
  altgnam(5,46)  = 'I C 2 M'
  altgnam(6,46)  = 'I M 2 A'
  altgnam(1,47)  = 'P M M M'
  altgnam(1,48)  = 'P N N N'
  altgnam(1,49)  = 'P C C M'
  altgnam(2,49)  = 'P M A A'
  altgnam(3,49)  = 'P B M B'
  altgnam(1,50)  = 'P B A N'
  altgnam(2,50)  = 'P N C B'
  altgnam(3,50)  = 'P C N A'
  altgnam(1,51)  = 'P M M A'
  altgnam(2,51)  = 'P M M B'
  altgnam(3,51)  = 'P B M M'
  altgnam(4,51)  = 'P C M M'
  altgnam(5,51)  = 'P M C M'
  altgnam(6,51)  = 'P M A M'
  altgnam(1,52)  = 'P N N A'
  altgnam(2,52)  = 'P N N B'
  altgnam(3,52)  = 'P B N N'
  altgnam(4,52)  = 'P C N N'
  altgnam(5,52)  = 'P N C N'
  altgnam(6,52)  = 'P N A N'
  altgnam(1,53)  = 'P M N A'
  altgnam(2,53)  = 'P N M B'
  altgnam(3,53)  = 'P B M N'
  altgnam(4,53)  = 'P C N M'
  altgnam(5,53)  = 'P N C M'
  altgnam(6,53)  = 'P M A N'
  altgnam(1,54)  = 'P C C A'
  altgnam(2,54)  = 'P C C B'
  altgnam(3,54)  = 'P B A A'
  altgnam(4,54)  = 'P C A A'
  altgnam(5,54)  = 'P B C B'
  altgnam(6,54)  = 'P B A B'
  altgnam(1,55)  = 'P B A M'
  altgnam(2,55)  = 'P M C B'
  altgnam(3,55)  = 'P C M A'
  altgnam(1,56)  = 'P C C N'
  altgnam(2,56)  = 'P N A A'
  altgnam(3,56)  = 'P B N B'
  altgnam(1,57)  = 'P B C M'
  altgnam(2,57)  = 'P C A M'
  altgnam(3,57)  = 'P M C A'
  altgnam(4,57)  = 'P M A B'
  altgnam(5,57)  = 'P B M A'
  altgnam(6,57)  = 'P C M B'
  altgnam(1,58)  = 'P N N M'
  altgnam(2,58)  = 'P M N N'
  altgnam(3,58)  = 'P N M N'
  altgnam(1,59)  = 'P M M N'
  altgnam(2,59)  = 'P N M M'
  altgnam(3,59)  = 'P M N M'
  altgnam(1,60)  = 'P B C N'
  altgnam(2,60)  = 'P C A N'
  altgnam(3,60)  = 'P N C A'
  altgnam(4,60)  = 'P N A B'
  altgnam(5,60)  = 'P B N A'
  altgnam(6,60)  = 'P C N B'
  altgnam(1,61)  = 'P B C A'
  altgnam(2,61)  = 'P C A B'
  altgnam(1,62)  = 'P N M A'
  altgnam(2,62)  = 'P M N B'
  altgnam(3,62)  = 'P B N M'
  altgnam(4,62)  = 'P C M N'
  altgnam(5,62)  = 'P M C N'
  altgnam(6,62)  = 'P N A M'
  altgnam(1,63)  = 'C M C M'
  altgnam(2,63)  = 'C C M M'
  altgnam(3,63)  = 'A M M A'
  altgnam(4,63)  = 'A M A M'
  altgnam(5,63)  = 'B B M M'
  altgnam(6,63)  = 'B M M B'
  altgnam(1,64)  = 'C M C A'
  altgnam(2,64)  = 'C C M B'
  altgnam(3,64)  = 'A B M A'
  altgnam(4,64)  = 'A C A M'
  altgnam(5,64)  = 'B B C M'
  altgnam(6,64)  = 'B M A B'
  altgnam(1,65)  = 'C M M M'
  altgnam(2,65)  = 'A M M M'
  altgnam(3,65)  = 'B M M M'
  altgnam(1,66)  = 'C C C M'
  altgnam(2,66)  = 'A M A A'
  altgnam(3,66)  = 'B B M B'
  altgnam(1,67)  = 'C M M A'
  altgnam(2,67)  = 'C M M B'
  altgnam(3,67)  = 'A B M M'
  altgnam(4,67)  = 'A C M M'
  altgnam(5,67)  = 'B M C M'
  altgnam(6,67)  = 'B M A M'
  altgnam(1,68)  = 'C C C A'
  altgnam(2,68)  = 'C C C B'
  altgnam(3,68)  = 'A B A A'
  altgnam(4,68)  = 'A C A A'
  altgnam(5,68)  = 'B B C B'
  altgnam(6,68)  = 'B B A B'
  altgnam(1,69)  = 'F M M M'
  altgnam(1,70)  = 'F D D D'
  altgnam(1,71)  = 'I M M M'
  altgnam(1,72)  = 'I B A M'
  altgnam(2,72)  = 'I M C B'
  altgnam(3,72)  = 'I C M A'
  altgnam(1,73)  = 'I B C A'
  altgnam(2,73)  = 'I C A B'
  altgnam(1,74)  = 'I M M A'
  altgnam(2,74)  = 'I M M B'
  altgnam(3,74)  = 'I B M M'
  altgnam(4,74)  = 'I C M M'
  altgnam(5,74)  = 'I M C M'
  altgnam(6,74)  = 'I M A M'
!****************************
!  Fitting variable labels  *
!****************************
  do i = 1,9
    do j = 1,63
      nfitlabel(j,i) = 0_i4
      do k = 1,30
        fitlabel(k,j,i) = 'Unknown type   '
        fitlabelunits(k,j,i) = 'Unknown   '
      enddo
    enddo
  enddo
!----------------------------
!  UFF bond order defaults  !
!----------------------------
  UFFbondorder(1) = 1.0_dp       ! Single bond
  UFFbondorder(2) = 2.0_dp       ! Double bond
  UFFbondorder(3) = 3.0_dp       ! Treble bond
  UFFbondorder(4) = 4.0_dp       ! Quadruple bond
  UFFbondorder(5) = 1.5_dp       ! Resonant bond
  UFFbondorder(6) = 1.41_dp      ! Amide bond
  UFFbondorder(7) = 1.0_dp       ! Custom bond -> default as single
  UFFbondorder(8) = 0.5_dp       ! Half bond
  UFFbondorder(9) = 0.25_dp      ! Quarter bond
  UFFbondorder(10) = 0.333333_dp ! Third bond
!----------------------
!  General parameters -
!----------------------
  fitlabel(1,1,1)='Shell position '
  fitlabel(1,2,1)='BSM radius     '
  fitlabel(1,3,1)='Energy shift   '
  fitlabel(1,4,1)='Charge         '
  fitlabel(1,5,1)='Polarisability '
  fitlabel(1,6,1)='Core charge    '
  fitlabel(1,7,1)='Species epsilon'
  fitlabel(1,8,1)='Species sigma  '
  fitlabel(1,9,1)='Atom A (lenn)  '
  fitlabel(1,10,1)='Atom B (lenn)  '
  fitlabel(1,11,1)='EEM chi        '
  fitlabel(1,12,1)='EEM mu         '
  fitlabel(1,13,1)='QEq radius     '
  fitlabel(1,14,1)='S and M zeta   '
  fitlabel(1,15,1)='S and M Znuc   '
  fitlabel(1,16,1)='Not in use     '
  fitlabel(1,17,1)='Not in use     '
  fitlabel(1,18,1)='Not in use     '
  fitlabel(1,19,1)='Not in use     '
  fitlabel(1,20,1)='UFF radius     '
  fitlabel(1,21,1)='UFF theta      '
  fitlabel(1,22,1)='UFF distance   '
  fitlabel(1,23,1)='UFF energy     '
  fitlabel(1,24,1)='UFF zeta       '
  fitlabel(1,25,1)='UFF Z effective'
  fitlabel(1,26,1)='UFF torsion    '
  fitlabel(1,27,1)='Bond charge    '
  fitlabel(1,28,1)='Plane pot A    '
  fitlabel(1,29,1)='Plane pot B    '
  fitlabel(1,30,1)='UFF Koop       '
  fitlabel(1,31,1)='UFF theta oop  '
  fitlabel(1,32,1)='UFF chi        '
  fitlabel(1,33,1)='EVB Q1         '
  fitlabel(1,34,1)='EVB Q2         '
  fitlabel(1,35,1)='EVB Q3         '
  fitlabel(1,36,1)='EVB Q4         '
  fitlabel(1,37,1)='EVB Q5         '
  fitlabel(1,38,1)='Ashift         '
  fitlabel(1,39,1)='Gasteiger A    '
  fitlabel(1,40,1)='Gasteiger B    '
  fitlabel(1,41,1)='Gasteiger C    '
  fitlabel(1,42,1)='Gasteiger Q0   '
!-----------------------
!  Two-body parameters -
!-----------------------
  nfitlabel(1,2) = 3
  fitlabel(1,1,2)='Buckingham A   '
  fitlabel(2,1,2)='Buckingham rho '
  fitlabel(3,1,2)='Buckingham C   '
  fitlabelunits(1,1,2)='eV        '
  fitlabelunits(2,1,2)='Ang       '
  fitlabelunits(3,1,2)='eV*Ang^6  '
!
  nfitlabel(2,2) = 2
  fitlabel(1,2,2)='Lennard-Jones A'
  fitlabel(2,2,2)='Lennard-Jones B'
  fitlabelunits(1,2,2)='eV*Ang^m  '
  fitlabelunits(2,2,2)='eV*Ang^n  '
!
  nfitlabel(3,2) = 3
  fitlabel(1,3,2)='Morse  De      '
  fitlabel(2,3,2)='Morse  a0      '
  fitlabel(3,3,2)='Morse  r0      '
  fitlabelunits(1,3,2)='eV        '
  fitlabelunits(2,3,2)='Ang^-1    '
  fitlabelunits(3,3,2)='Ang       '
!
  nfitlabel(4,2) = 3
  fitlabel(1,4,2)='Morse  De      '
  fitlabel(2,4,2)='Morse  a0      '
  fitlabel(3,4,2)='Morse  r0      '
  fitlabelunits(1,4,2)='eV        '
  fitlabelunits(2,4,2)='Ang^-1    '
  fitlabelunits(3,4,2)='Ang       '
!
  nfitlabel(5,2) = 4
  fitlabel(1,5,2)='Harmonic k2    '
  fitlabel(2,5,2)='Harmonic r0    '
  fitlabel(3,5,2)='Harmonic k3    '
  fitlabel(4,5,2)='Harmonic k4    '
  fitlabelunits(1,5,2)='eV*Ang^-2 '
  fitlabelunits(2,5,2)='Ang       '
  fitlabelunits(3,5,2)='eV*Ang^-3 '
  fitlabelunits(4,5,2)='eV*Ang^-4 '
!
  nfitlabel(6,2) = 4
  fitlabel(1,6,2)='Harmonic k2    '
  fitlabel(2,6,2)='Harmonic r0    '
  fitlabel(3,6,2)='Harmonic k3    '
  fitlabel(4,6,2)='Harmonic k4    '
  fitlabelunits(1,6,2)='eV*Ang^-2 '
  fitlabelunits(2,6,2)='Ang       '
  fitlabelunits(3,6,2)='eV*Ang^-3 '
  fitlabelunits(4,6,2)='eV*Ang^-4 '
!
  nfitlabel(7,2) = 3
  fitlabel(1,7,2)='General A      '
  fitlabel(2,7,2)='General rho    '
  fitlabel(3,7,2)='General C      '
  fitlabelunits(1,7,2)='eV*Ang^m  '
  fitlabelunits(2,7,2)='Ang       '
  fitlabelunits(3,7,2)='eV*Ang^n  '
!
  nfitlabel(8,2) = 2
  fitlabel(1,8,2)='Spring k 2     '
  fitlabel(2,8,2)='Spring k 4     '
  fitlabelunits(1,8,2)='eV*Ang^-2 '
  fitlabelunits(2,8,2)='eV*Ang^-4 '
!
  nfitlabel(9,2) = 1
  fitlabel(1,9,2)='Coulomb sub    '
  fitlabelunits(1,9,2)='None      '
!
  nfitlabel(10,2) = 4
  fitlabel(1,10,2)='Buck4 A        '
  fitlabel(2,10,2)='Buck4 rho      '
  fitlabel(3,10,2)='Buck4 C        '
  fitlabel(4,10,2)='Buck4 minimum  '
  fitlabelunits(1,10,2)='eV        '
  fitlabelunits(2,10,2)='Ang       '
  fitlabelunits(3,10,2)='eV*Ang^6  '
  fitlabelunits(4,10,2)='Ang       '
!
  nfitlabel(11,2) = 1
  fitlabel(1,11,2)='Spline shift   '
  fitlabelunits(1,11,2)='eV        '
!
  nfitlabel(12,2) = 2
  fitlabel(1,12,2)='Lennard epsilon'
  fitlabel(2,12,2)='Lennard sigma/z'
  fitlabelunits(1,12,2)='eV        '
  fitlabelunits(2,12,2)='Ang       '
!
  nfitlabel(13,2) = 2
  fitlabel(1,13,2)='Lennard epsilon'
  fitlabel(2,13,2)='Lennard sigma/m'
  fitlabelunits(1,13,2)='eV        '
  fitlabelunits(2,13,2)='Ang       '
!
  nfitlabel(14,2) = 2
  fitlabel(1,14,2)='BSM force const'
  fitlabel(2,14,2)='BSM radius     '
  fitlabelunits(1,14,2)='eV*Ang^-2 '
  fitlabelunits(2,14,2)='Ang       '
!
  nfitlabel(15,2) = 3
  fitlabel(1,15,2)='SW2 A          '
  fitlabel(2,15,2)='SW2 rho        '
  fitlabel(3,15,2)='SW2 B          '
  fitlabelunits(1,15,2)='eV        '
  fitlabelunits(2,15,2)='Ang       '
  fitlabelunits(3,15,2)='Ang^4     '
!
  nfitlabel(16,2) = 3
  fitlabel(1,16,2)='Inv Gauss K    '
  fitlabel(2,16,2)='Inv Gauss a    '
  fitlabel(3,16,2)='Inv Gauss r0   '
  fitlabelunits(1,16,2)='eV        '
  fitlabelunits(2,16,2)='Ang^-2    '
  fitlabelunits(3,16,2)='Ang       '
!
  nfitlabel(17,2) = 3
  fitlabel(1,17,2)='BSM force const'
  fitlabel(2,17,2)='BSM rho        '
  fitlabel(3,17,2)='BSM radius     '
  fitlabelunits(1,17,2)='eV        '
  fitlabelunits(2,17,2)='Ang^-1    '
  fitlabelunits(3,17,2)='Ang       '
!
  nfitlabel(18,2) = 11
  fitlabel(1,18,2)='Damped Disp C6 '
  fitlabel(2,18,2)='Damped Disp C8 '
  fitlabel(3,18,2)='Damped Disp B6 '
  fitlabel(4,18,2)='Damped Disp B8 '
  fitlabel(10,18,2)='Damped Disp C10'
  fitlabel(11,18,2)='Damped Disp B10'
  fitlabelunits(1,18,2)='eV*Ang^6  '
  fitlabelunits(2,18,2)='eV*Ang^8  '
  fitlabelunits(3,18,2)='Ang^-1    '
  fitlabelunits(4,18,2)='Ang^-1    '
  fitlabelunits(10,18,2)='eV*Ang^10 '
  fitlabelunits(11,18,2)='Ang^-1    '
!
  nfitlabel(19,2) = 0
!
  nfitlabel(20,2) = 3
  fitlabel(1,20,2)='Rydberg A      '
  fitlabel(2,20,2)='Rydberg B      '
  fitlabel(3,20,2)='Rydberg r0     '
  fitlabelunits(1,20,2)='eV        '
  fitlabelunits(2,20,2)='None      '
  fitlabelunits(3,20,2)='Ang       '
!
  nfitlabel(21,2) = 0
!
  nfitlabel(22,2) = 1
  fitlabel(1,22,2)='Coulomb taper C'
  fitlabelunits(1,22,2)='eV        '
!
  nfitlabel(23,2) = 18
  fitlabel(10,23,2)='Polynomial c0  '
  fitlabel(11,23,2)='Polynomial c1  '
  fitlabel(12,23,2)='Polynomial c2  '
  fitlabel(13,23,2)='Polynomial c3  '
  fitlabel(14,23,2)='Polynomial c4  '
  fitlabel(15,23,2)='Polynomial c5  '
  fitlabel(16,23,2)='Polynomial c6  '
  fitlabel(17,23,2)='Polynomial c7  '
  fitlabel(18,23,2)='Polynomial c8  '
  fitlabelunits(10,23,2)='eV        '
  fitlabelunits(11,23,2)='eV*Ang^-1 '
  fitlabelunits(12,23,2)='eV*Ang^-2 '
  fitlabelunits(13,23,2)='eV*Ang^-3 '
  fitlabelunits(14,23,2)='eV*Ang^-4 '
  fitlabelunits(15,23,2)='eV*Ang^-5 '
  fitlabelunits(16,23,2)='eV*Ang^-6 '
  fitlabelunits(17,23,2)='eV*Ang^-7 '
  fitlabelunits(18,23,2)='eV*Ang^-8 '
!
  nfitlabel(24,2) = 1
  fitlabel(1,24,2)='Coulomb erfc a '
  fitlabelunits(1,24,2)='Ang       '
!
  nfitlabel(25,2) = 3
  fitlabel(1,25,2)='CovExp De      '
  fitlabel(2,25,2)='CovExp a0      '
  fitlabel(3,25,2)='CovExp r0      '
  fitlabelunits(1,25,2)='eV        '
  fitlabelunits(2,25,2)='Ang^-1    '
  fitlabelunits(3,25,2)='Ang       '
!
  nfitlabel(26,2) = 3
  fitlabel(1,26,2)='Fermi-Dirac a  '
  fitlabel(2,26,2)='Fermi-Dirac b  '
  fitlabel(3,26,2)='Fermi-Dirac c  '
  fitlabelunits(1,26,2)='eV        '
  fitlabelunits(2,26,2)='Ang^-1    '
  fitlabelunits(3,26,2)='Ang       '
!
  nfitlabel(27,2) = 3
  fitlabel(1,27,2)='L-J buffered A '
  fitlabel(2,27,2)='L-J buffered B '
  fitlabel(3,27,2)='L-J buffered r0'
  fitlabelunits(1,27,2)='eV*Ang^m  '
  fitlabelunits(2,27,2)='eV*Ang^n  '
  fitlabelunits(3,27,2)='Ang       '
!
  nfitlabel(28,2) = 2
  fitlabel(1,28,2)='SqrHarmonic k2 '
  fitlabel(2,28,2)='SqrHarmonic r0 '
  fitlabelunits(1,28,2)='eV*Ang^-4 '
  fitlabelunits(2,28,2)='Ang       '
!
  nfitlabel(29,2) = 2
  fitlabel(1,29,2)='SqrHarmonic k2 '
  fitlabel(2,29,2)='SqrHarmonic r0 '
  fitlabelunits(1,29,2)='eV*Ang^-4 '
  fitlabelunits(2,29,2)='Ang       '
!
  nfitlabel(30,2) = 3
  fitlabel(1,30,2)='Tsuneyuki Q1   '
  fitlabel(2,30,2)='Tsuneyuki Q2   '
  fitlabel(3,30,2)='Tsuneyuki zeta '
  fitlabelunits(1,30,2)='a.u.      '
  fitlabelunits(2,30,2)='a.u.      '
  fitlabelunits(3,30,2)='Ang^-1    '
!
  nfitlabel(31,2) = 3
  fitlabel(1,31,2)='BSM force const'
  fitlabel(2,31,2)='BSM rho        '
  fitlabel(3,31,2)='BSM radius     '
  fitlabelunits(1,31,2)='eV        '
  fitlabelunits(2,31,2)='Ang^-1    '
  fitlabelunits(3,31,2)='Ang       '
!
  nfitlabel(32,2) = 4
  fitlabel(1,32,2)='SW2jb A        '
  fitlabel(2,32,2)='SW2jb rho      '
  fitlabel(3,32,2)='SW2jb B        '
  fitlabel(4,32,2)='SW2jb Qs       '
  fitlabelunits(1,32,2)='eV        '
  fitlabelunits(2,32,2)='Ang       '
  fitlabelunits(3,32,2)='Ang^4     '
  fitlabelunits(4,32,2)='a.u.      '
!
  nfitlabel(33,2) = 2
  fitlabel(1,33,2)='Cosh spring k  '
  fitlabel(2,33,2)='Cosh spring d  '
  fitlabelunits(1,33,2)='eV*Ang^-2 '
  fitlabelunits(2,33,2)='Ang       '
!
  nfitlabel(34,2) = 2
  fitlabel(1,34,2)='EAM pot g      '
  fitlabel(2,34,2)='EAM pot beta   '
  fitlabelunits(1,34,2)='None      '
  fitlabelunits(2,34,2)='Ang^-1    '
!
  nfitlabel(35,2) = 15
  fitlabel(1,35,2) ='Polynomial r0  '
  fitlabel(10,35,2)='Polynomial c0  '
  fitlabel(11,35,2)='Polynomial c1  '
  fitlabel(12,35,2)='Polynomial c2  '
  fitlabel(13,35,2)='Polynomial c3  '
  fitlabel(14,35,2)='Polynomial c4  '
  fitlabel(15,35,2)='Polynomial c5  '
  fitlabelunits(1,35,2) ='Ang       '
  fitlabelunits(10,35,2)='eV*Ang^-2 '
  fitlabelunits(11,35,2)='eV*Ang^-3 '
  fitlabelunits(12,35,2)='eV*Ang^-4 '
  fitlabelunits(13,35,2)='eV*Ang^-5 '
  fitlabelunits(14,35,2)='eV*Ang^-6 '
  fitlabelunits(15,35,2)='eV*Ang^-7 '
!
  nfitlabel(36,2) = 0
  fitlabel(1,36,2)='Q over R2      '
!
  nfitlabel(37,2) = 2
  fitlabel(1,37,2)='Longitudinal K '
  fitlabel(2,37,2)='Transverse   K '
  fitlabelunits(1,37,2)='eV*Ang^-2 '
  fitlabelunits(2,37,2)='eV*Ang^-2 '
!
  nfitlabel(38,2) = 21
  fitlabel(1,38,2) ='SRglue d      '
  fitlabel(10,38,2)='SRglue a14    '
  fitlabel(11,38,2)='SRglue a13    '
  fitlabel(12,38,2)='SRglue a12    '
  fitlabel(13,38,2)='SRglue a11    '
  fitlabel(14,38,2)='SRglue a10    '
  fitlabel(15,38,2)='SRglue a26    '
  fitlabel(16,38,2)='SRglue a25    '
  fitlabel(17,38,2)='SRglue a24    '
  fitlabel(18,38,2)='SRglue a23    '
  fitlabel(19,38,2)='SRglue a22    '
  fitlabel(20,38,2)='SRglue a21    '
  fitlabel(21,38,2)='SRglue a20    '
  fitlabelunits(1,38,2)='Ang       '
  fitlabelunits(10,38,2)='eV*Ang^-4 '
  fitlabelunits(11,38,2)='eV*Ang^-3 '
  fitlabelunits(12,38,2)='eV*Ang^-2 '
  fitlabelunits(13,38,2)='eV*Ang^-1 '
  fitlabelunits(14,38,2)='eV        '
  fitlabelunits(15,38,2)='eV*Ang^-6 '
  fitlabelunits(16,38,2)='eV*Ang^-5 '
  fitlabelunits(17,38,2)='eV*Ang^-4 '
  fitlabelunits(18,38,2)='eV*Ang^-3 '
  fitlabelunits(19,38,2)='eV*Ang^-2 '
  fitlabelunits(20,38,2)='eV*Ang^-1 '
  fitlabelunits(21,38,2)='eV        '
!
  nfitlabel(39,2) = 3
  fitlabel(1,39,2)='Morse-etaper De'
  fitlabel(2,39,2)='Morse-etaper a0'
  fitlabel(3,39,2)='Morse-etaper r0'
  fitlabelunits(1,39,2)='eV        '
  fitlabelunits(2,39,2)='Ang^-1    '
  fitlabelunits(3,39,2)='Ang       '
!
  nfitlabel(40,2) = 3
  fitlabel(1,40,2)='Morse-etaper De'
  fitlabel(2,40,2)='Morse-etaper a0'
  fitlabel(3,40,2)='Morse-etaper r0'
  fitlabelunits(1,40,2)='eV        '
  fitlabelunits(2,40,2)='Ang^-1    '
  fitlabelunits(3,40,2)='Ang       '
!
  nfitlabel(41,2) = 4
  fitlabel(1,41,2)='Mei-Dprt phi0  '
  fitlabel(2,41,2)='Mei-Dprt delta '
  fitlabel(3,41,2)='Mei-Dprt gamma '
  fitlabel(4,41,2)='Mei-Dprt r0    '
  fitlabelunits(1,41,2)='eV        '
  fitlabelunits(2,41,2)='None      '
  fitlabelunits(3,41,2)='None      '
  fitlabelunits(4,41,2)='Ang       '
!
  nfitlabel(42,2) = 3
  fitlabel(1,42,2)='erferfc A      '
  fitlabel(2,42,2)='erferfc alpha  '
  fitlabel(3,42,2)='erferfc beta   '
  fitlabelunits(1,42,2)='eV        '
  fitlabelunits(2,42,2)='Ang       '
  fitlabelunits(3,42,2)='Ang       '
!
  nfitlabel(43,2) = 2
  fitlabel(1,43,2)='reperfc A      '
  fitlabel(2,43,2)='reperfc beta   '
  fitlabelunits(1,43,2)='eV        '
  fitlabelunits(2,43,2)='Ang       '
!
  nfitlabel(44,2) = 2
  fitlabel(1,44,2)='erfpot A       '
  fitlabel(2,44,2)='erfpot alpha   '
  fitlabelunits(1,44,2)='eV        '
  fitlabelunits(2,44,2)='Ang       '
!
  nfitlabel(45,2) = 7
  fitlabel(1,45,2)='Baskes Ec      '
  fitlabel(2,45,2)='Baskes A       '
  fitlabel(3,45,2)='Baskes alpha   '
  fitlabel(4,45,2)='Baskes r0      '
  fitlabel(5,45,2)='Baskes Z       '
  fitlabel(6,45,2)='Baskes rho0    '
  fitlabel(7,45,2)='Baskes d       '
  fitlabelunits(1,45,2)='eV        '
  fitlabelunits(2,45,2)='None      '
  fitlabelunits(3,45,2)='None      '
  fitlabelunits(4,45,2)='Ang       '
  fitlabelunits(5,45,2)='None      '
  fitlabelunits(6,45,2)='None      '
  fitlabelunits(7,45,2)='None      '
!
  nfitlabel(46,2) = 4
  fitlabel(1,46,2)='VBO 2body c    '
  fitlabel(2,46,2)='VBO 2body gamma'
  fitlabel(3,46,2)='VBO 2body R0   '
  fitlabel(4,46,2)='VBO 2body delta'
  fitlabelunits(1,46,2)='eV        '
  fitlabelunits(2,46,2)='None      '
  fitlabelunits(3,46,2)='Ang       '
  fitlabelunits(4,46,2)='Ang       '
!
  nfitlabel(47,2) = 5
  fitlabel(1,47,2)='Exp-powers A   '
  fitlabel(2,47,2)='Exp-powers B0  '
  fitlabel(3,47,2)='Exp-powers B1  '
  fitlabel(4,47,2)='Exp-powers B2  '
  fitlabel(5,47,2)='Exp-powers B3  '
  fitlabelunits(1,47,2)='eV        '
  fitlabelunits(2,47,2)='None      '
  fitlabelunits(3,47,2)='Ang^-1    '
  fitlabelunits(4,47,2)='Ang^-2    '
  fitlabelunits(5,47,2)='Ang^-3    '
!
  nfitlabel(48,2) = 3
  fitlabel(1,48,2)='Grimme_C6  C6  '
  fitlabel(2,48,2)='Grimme_C6  d   '
  fitlabel(3,48,2)='Grimme_C6  r0  '
  fitlabelunits(1,48,2)='eV*Ang^6  '
  fitlabelunits(2,48,2)='None      '
  fitlabelunits(3,48,2)='Ang       '
!
  nfitlabel(49,2) = 4
  fitlabel(1,49,2)='CFM harm   k   '
  fitlabel(2,49,2)='CFM harm   r0  '
  fitlabel(3,49,2)='CFM harm   R   '
  fitlabel(4,49,2)='CFM harm   w   '
  fitlabelunits(1,49,2)='eV*Ang^-2 '
  fitlabelunits(2,49,2)='Ang       '
  fitlabelunits(3,49,2)='Ang       '
  fitlabelunits(4,49,2)='Ang       '
!
  nfitlabel(50,2) = 5
  fitlabel(1,50,2)='CFM gauss  k   '
  fitlabel(2,50,2)='CFM gauss  zeta'
  fitlabel(3,50,2)='CFM gauss  r0  '
  fitlabel(4,50,2)='CFM gauss  R   '
  fitlabel(5,50,2)='CFM gauss  w   '
  fitlabelunits(1,50,2)='eV        '
  fitlabelunits(2,50,2)='Ang^-2    '
  fitlabelunits(3,50,2)='Ang       '
  fitlabelunits(4,50,2)='Ang       '
  fitlabelunits(5,50,2)='Ang       '
!
  nfitlabel(51,2) = 4
  fitlabel(1,51,2)='CFM power  A   '
  fitlabel(2,51,2)='CFM power  pwr '
  fitlabel(3,51,2)='CFM power  R   '
  fitlabel(4,51,2)='CFM power  w   '
  fitlabelunits(1,51,2)='eV*Ang^pwr'
  fitlabelunits(2,51,2)='None      '
  fitlabelunits(3,51,2)='Ang       '
  fitlabelunits(4,51,2)='Ang       '
!
  nfitlabel(52,2) = 5
  fitlabel(1,52,2)='CFM fermi  k   '
  fitlabel(2,52,2)='CFM fermi  zeta'
  fitlabel(3,52,2)='CFM fermi  r0  '
  fitlabel(4,52,2)='CFM fermi  R   '
  fitlabel(5,52,2)='CFM fermi  w   '
  fitlabelunits(1,52,2)='eV        '
  fitlabelunits(2,52,2)='Ang^-1    '
  fitlabelunits(3,52,2)='Ang       '
  fitlabelunits(4,52,2)='Ang       '
  fitlabelunits(5,52,2)='Ang       '
!
  nfitlabel(53,2) = 2
  fitlabel(1,53,2)='gCoulomb A     '
  fitlabel(2,53,2)='gCoulomb gamma '
  fitlabelunits(1,53,2)='a.u.      '
  fitlabelunits(2,53,2)='Ang^3     '
!
  nfitlabel(54,2) = 2
  fitlabel(1,54,2)='BeckeJ_C6  C6  '
  fitlabel(2,54,2)='BeckeJ_C6  r0  '
  fitlabelunits(1,54,2)='eV*Ang^6  '
  fitlabelunits(2,54,2)='Ang       '
!
  nfitlabel(55,2) = 9
  fitlabel(1,55,2)='Baskes Ec      '
  fitlabel(2,55,2)='Baskes A       '
  fitlabel(3,55,2)='Baskes alpha   '
  fitlabel(4,55,2)='Baskes r0      '
  fitlabel(5,55,2)='Baskes Z       '
  fitlabel(6,55,2)='Baskes rho0    '
  fitlabel(7,55,2)='Baskes d       '
  fitlabel(8,55,2)='Baskes gamma   '
  fitlabel(9,55,2)='Baskes lambda  '
  fitlabelunits(1,55,2)='eV        '
  fitlabelunits(2,55,2)='None      '
  fitlabelunits(3,55,2)='None      '
  fitlabelunits(4,55,2)='Ang       '
  fitlabelunits(5,55,2)='None      '
  fitlabelunits(6,55,2)='None      '
  fitlabelunits(7,55,2)='None      '
  fitlabelunits(8,55,2)='Ang       '
  fitlabelunits(9,55,2)='None      '
!
  nfitlabel(56,2) = 9
  fitlabel(1,56,2)='ZBL b1         '
  fitlabel(2,56,2)='ZBL c1         '
  fitlabel(3,56,2)='ZBL b2         '
  fitlabel(4,56,2)='ZBL c2         '
  fitlabel(5,56,2)='ZBL b3         '
  fitlabel(6,56,2)='ZBL c3         '
  fitlabel(7,56,2)='ZBL b4         '
  fitlabel(8,56,2)='ZBL c4         '
  fitlabel(9,56,2)='ZBL d          '
  fitlabelunits(1,56,2)='None      '
  fitlabelunits(2,56,2)='None      '
  fitlabelunits(3,56,2)='None      '
  fitlabelunits(4,56,2)='None      '
  fitlabelunits(5,56,2)='None      '
  fitlabelunits(6,56,2)='None      '
  fitlabelunits(7,56,2)='None      '
  fitlabelunits(8,56,2)='None      '
  fitlabelunits(9,56,2)='None      '
!
  nfitlabel(57,2) = 5
  fitlabel(1,57,2)='MM3buck A      '
  fitlabel(2,57,2)='MM3buck B      '
  fitlabel(3,57,2)='MM3buck C      '
  fitlabel(4,57,2)='MM3buck epsilon'
  fitlabel(5,57,2)='MM3buck rv     '
  fitlabelunits(1,57,2)='None      '
  fitlabelunits(2,57,2)='None      '
  fitlabelunits(3,57,2)='None      '
  fitlabelunits(4,57,2)='eV        '
  fitlabelunits(5,57,2)='Ang       '
!
  nfitlabel(58,2) = 5
  fitlabel(1,58,2)='MM3stretch k2  '
  fitlabel(2,58,2)='MM3stretch r0  '
  fitlabel(3,58,2)='MM3stretch A   '
  fitlabel(4,58,2)='MM3stretch B   '
  fitlabel(5,58,2)='MM3stretch C   '
  fitlabelunits(1,58,2)='eV*Ang^-2 '
  fitlabelunits(2,58,2)='Ang       '
  fitlabelunits(3,58,2)='None      '
  fitlabelunits(4,58,2)='Ang^-1    '
  fitlabelunits(5,58,2)='Ang^-1    '
!
  nfitlabel(59,2) = 2
  fitlabel(1,59,2)='Lennard epsilon'
  fitlabel(2,59,2)='Lennard sigma/m'
  fitlabelunits(1,59,2)='eV        '
  fitlabelunits(2,59,2)='Ang       '
!
  nfitlabel(60,2) = 4
  fitlabel(1,60,2)='Buff_LJ epsilon'
  fitlabel(2,60,2)='Buff_LJ r0     '
  fitlabel(3,60,2)='Buff_LJ delta  '
  fitlabel(4,60,2)='Buff_LJ gamma  '
  fitlabelunits(1,60,2)='eV        '
  fitlabelunits(2,60,2)='Ang       '
  fitlabelunits(3,60,2)='None      '
  fitlabelunits(4,60,2)='None      '
!
  nfitlabel(61,2) = 2
  fitlabel(1,61,2)='Slater A   '
  fitlabel(2,61,2)='Slater B   '
  fitlabelunits(1,61,2)='eV        '
  fitlabelunits(2,61,2)='Ang^-1    '
!
  nfitlabel(62,2) = 1
  fitlabel(1,62,2)='j2 J2          '
  fitlabelunits(1,62,2)='eV        '
!
  nfitlabel(63,2) = 5
  fitlabel(1,63,2)='SW2 Gen A      '
  fitlabel(2,63,2)='SW2 Gen sigma  '
  fitlabel(3,63,2)='SW2 Gen B      '
  fitlabel(4,63,2)='SW2 Gen p      '
  fitlabel(5,63,2)='SW2 Gen q      '
  fitlabelunits(1,63,2)='eV        '
  fitlabelunits(2,63,2)='Ang       '
  fitlabelunits(3,63,2)='None      '
  fitlabelunits(4,63,2)='None      '
  fitlabelunits(5,63,2)='None      '
!
  nfitlabel(64,2) = 3
  fitlabel(1,64,2)='Screened Esat  '
  fitlabel(2,64,2)='Screened r_me  '
  fitlabel(3,64,2)='Screened sig_e '
  fitlabelunits(1,64,2)='None      '
  fitlabelunits(2,64,2)='Ang       '
  fitlabelunits(3,64,2)='Ang       '
!-------------------------
!  Three-body parameters -
!-------------------------
  nfitlabel(1,3) = 4
  fitlabel(1,1,3)='Three-body cnst'
  fitlabel(2,1,3)='Three-body angl'
  fitlabel(3,1,3)='Three-body k4  '
  fitlabel(4,1,3)='Three-body k3  '
  fitlabelunits(1,1,3)='eV*rad^-2 '
  fitlabelunits(2,1,3)='degrees   '
  fitlabelunits(3,1,3)='eV*rad^-4 '
  fitlabelunits(4,1,3)='eV*rad^-3 '
!
  nfitlabel(2,3) = 4
  fitlabel(1,2,3)='Three-body cnst'
  fitlabel(2,2,3)='Three-body angl'
  fitlabel(3,2,3)='Three-body rho1'
  fitlabel(4,2,3)='Three-body rho2'
  fitlabelunits(1,2,3)='eV*rad^-2 '
  fitlabelunits(2,2,3)='degrees   '
  fitlabelunits(3,2,3)='Ang       '
  fitlabelunits(4,2,3)='Ang       '
!
  nfitlabel(3,3) = 1
  fitlabel(1,3,3)='Three-body cnst'
  fitlabelunits(1,3,3)='eV*Ang^9  '
!
  nfitlabel(4,3) = 4
  fitlabel(1,4,3)='Three-body cnst'
  fitlabel(2,4,3)='Three-body rho1'
  fitlabel(3,4,3)='Three-body rho2'
  fitlabel(4,4,3)='Three-body rho3'
  fitlabelunits(1,4,3)='eV        '
  fitlabelunits(2,4,3)='Ang       '
  fitlabelunits(3,4,3)='Ang       '
  fitlabelunits(4,4,3)='Ang       '
!
  nfitlabel(5,3) = 4
  fitlabel(1,5,3)='Three-body cnst'
  fitlabel(2,5,3)='Three-body angl'
  fitlabel(3,5,3)='Three-body rho1'
  fitlabel(4,5,3)='Three-body rho2'
  fitlabelunits(1,5,3)='eV        '
  fitlabelunits(2,5,3)='degrees   '
  fitlabelunits(3,5,3)='Ang       '
  fitlabelunits(4,5,3)='Ang       '
!
  nfitlabel(6,3) = 4
  fitlabel(1,6,3)='Three-body cnst'
  fitlabel(3,6,3)='Three-body r0/1'
  fitlabel(4,6,3)='Three-body r0/2'
  fitlabelunits(1,6,3)='eV*Ang^-2 '
  fitlabelunits(2,6,3)='Ang       '
  fitlabelunits(3,6,3)='Ang       '
!
  nfitlabel(7,3) = 2
  fitlabel(1,7,3)='Three-body cnst'
  fitlabel(2,7,3)='Three-body r0  '
  fitlabelunits(1,7,3)='eV*Ang^-2 '
  fitlabelunits(2,7,3)='Ang       '
!
  nfitlabel(8,3) = 4
  fitlabel(1,8,3)='Three-body cnst'
  fitlabel(2,8,3)='Three-body angl'
  fitlabel(3,8,3)='Three-body rho1'
  fitlabel(4,8,3)='Three-body rho2'
  fitlabelunits(1,8,3)='eV*rad^-2 '
  fitlabelunits(2,8,3)='degrees   '
  fitlabelunits(3,8,3)='Ang       '
  fitlabelunits(4,8,3)='Ang       '
!
  nfitlabel(9,3) = 4
  fitlabel(1,9,3)='Three-body cnst'
  fitlabel(2,9,3)='Three-body angl'
  fitlabel(3,9,3)='Three-body k3  '
  fitlabel(4,9,3)='Three-body k4  '
  fitlabelunits(1,9,3)='eV        '
  fitlabelunits(2,9,3)='degrees   '
  fitlabelunits(3,9,3)='eV        '
  fitlabelunits(4,9,3)='eV        '
!
  nfitlabel(10,3) = 16
  fitlabel(1,10,3)='Three-body cnst'
  fitlabel(2,10,3)='Three-body rho '
  fitlabel(3,10,3)='Three-body r0/1'
  fitlabel(4,10,3)='Three-body r0/2'
  fitlabel(5,10,3)='Three-body r0/3'
  fitlabel(6,10,3)='Three-body c0  '
  fitlabel(7,10,3)='Three-body c1  '
  fitlabel(8,10,3)='Three-body c2  '
  fitlabel(9,10,3)='Three-body c3  '
  fitlabel(10,10,3)='Three-body c4  '
  fitlabel(11,10,3)='Three-body c5  '
  fitlabel(12,10,3)='Three-body c6  '
  fitlabel(13,10,3)='Three-body c7  '
  fitlabel(14,10,3)='Three-body c8  '
  fitlabel(15,10,3)='Three-body c9  '
  fitlabel(16,10,3)='Three-body c10 '
  fitlabelunits(1,10,3)='eV        '
  fitlabelunits(2,10,3)='Ang^-1    '
  fitlabelunits(3,10,3)='Ang       '
  fitlabelunits(4,10,3)='Ang       '
  fitlabelunits(5,10,3)='Ang       '
  fitlabelunits(6,10,3)='None      '
  fitlabelunits(7,10,3)='Ang^-1    '
  fitlabelunits(8,10,3)='Ang^-2    '
  fitlabelunits(9,10,3)='Ang^-2    '
  fitlabelunits(10,10,3)='Ang^-3    '
  fitlabelunits(11,10,3)='Ang^-3    '
  fitlabelunits(12,10,3)='Ang^-3    '
  fitlabelunits(13,10,3)='Ang^-4    '
  fitlabelunits(14,10,3)='Ang^-4    '
  fitlabelunits(15,10,3)='Ang^-4    '
  fitlabelunits(16,10,3)='Ang^-4    '
!
  nfitlabel(11,3) = 5
  fitlabel(1,11,3)='Three-body K1  '
  fitlabel(2,11,3)='Three-body angl'
  fitlabel(3,11,3)='Three-body r0/1'
  fitlabel(4,11,3)='Three-body r0/2'
  fitlabel(5,11,3)='Three-body K2  '
  fitlabelunits(1,11,3)='eV/Ang/rad'
  fitlabelunits(2,11,3)='degrees   '
  fitlabelunits(3,11,3)='Ang       '
  fitlabelunits(4,11,3)='Ang       '
  fitlabelunits(5,11,3)='eV/Ang/rad'
!
  nfitlabel(12,3) = 2
  fitlabel(1,12,3)='Three-body cnst'
  fitlabel(2,12,3)='Three-body sign'
  fitlabelunits(1,12,3)='eV        '
  fitlabelunits(2,12,3)='None      '
!
  nfitlabel(13,3) = 4
  fitlabel(1,13,3)='Three-body cnst'
  fitlabel(2,13,3)='Three-body b   '
  fitlabel(3,13,3)='Three-body r0/1'
  fitlabel(4,13,3)='Three-body r0/2'
  fitlabelunits(1,13,3)='eV*Ang^-2 '
  fitlabelunits(2,13,3)='None      '
  fitlabelunits(3,13,3)='Ang       '
  fitlabelunits(4,13,3)='Ang       '
!
  nfitlabel(14,3) = 5
  fitlabel(1,14,3)='Three-body cnst'
  fitlabel(2,14,3)='Three-body angl'
  fitlabel(3,14,3)='Three-body rho1'
  fitlabel(4,14,3)='Three-body rho2'
  fitlabel(5,14,3)='Three-body Q   '
  fitlabelunits(1,14,3)='eV        '
  fitlabelunits(2,14,3)='degrees   '
  fitlabelunits(3,14,3)='Ang       '
  fitlabelunits(4,14,3)='Ang       '
  fitlabelunits(5,14,3)='a.u.      '
!
  nfitlabel(15,3) = 2
  fitlabel(1,15,3)='Three-body A   '
  fitlabel(2,15,3)='Three-body B   '
  fitlabelunits(1,15,3)='eV*Ang^m  '
  fitlabelunits(2,15,3)='eV*Ang^n  '
!
  nfitlabel(16,3) = 4
  fitlabel(1,16,3)='Three-body K   '
  fitlabel(2,16,3)='Three-body n   '
  fitlabel(3,16,3)='Three-body beta'
  fitlabel(4,16,3)='Three-body r0  '
  fitlabelunits(1,16,3)='eV        '
  fitlabelunits(2,16,3)='None      '
  fitlabelunits(3,16,3)='Ang^-1    '
  fitlabelunits(4,16,3)='Ang       '
!
  nfitlabel(17,3) = 2
  fitlabel(1,17,3)='Three-body cnst'
  fitlabel(2,17,3)='Three-body angl'
  fitlabelunits(1,17,3)='eV        '
  fitlabelunits(2,17,3)='degrees   '
!
  nfitlabel(18,3) = 5
  fitlabel(1,18,3)='Three-body K1  '
  fitlabel(2,18,3)='Three-body K2  '
  fitlabel(3,18,3)='Three-body r0/1'
  fitlabel(4,18,3)='Three-body r0/2'
  fitlabel(5,18,3)='Three-body angl'
  fitlabelunits(1,18,3)='eV*Ang^-1 '
  fitlabelunits(2,18,3)='eV*Ang^-1 '
  fitlabelunits(3,18,3)='Ang       '
  fitlabelunits(4,18,3)='Ang       '
  fitlabelunits(5,18,3)='degrees   '
!
  nfitlabel(19,3) = 1
  fitlabel(1,19,3)='Three-body Qsub'
  fitlabelunits(1,19,3)='None      '
!
  nfitlabel(20,3) = 5
  fitlabel(1,20,3)='Three-body K   '
  fitlabel(2,20,3)='Three-body rho2'
  fitlabel(3,20,3)='Three-body r120'
  fitlabel(4,20,3)='Three-body rho3'
  fitlabel(5,20,3)='Three-body r130'
  fitlabelunits(1,20,3)='eV        '
  fitlabelunits(2,20,3)='Ang^-1    '
  fitlabelunits(3,20,3)='Ang       '
  fitlabelunits(4,20,3)='Ang^-1    '
  fitlabelunits(5,20,3)='Ang       '
!
  nfitlabel(21,3) = 2
  fitlabel(1,21,3)='Three-body A   '
  fitlabel(2,21,3)='Three-body gamm'
  fitlabelunits(1,21,3)='a.u.      '
  fitlabelunits(2,21,3)='Ang^3     '
!
  nfitlabel(22,3) = 4
  fitlabel(1,22,3)='Three-body K   '
  fitlabel(2,22,3)='Three-body r0/1'
  fitlabel(3,22,3)='Three-body r0/2'
  fitlabel(4,22,3)='Three-body angl'
  fitlabelunits(1,22,3)='eV/A^mnr^p'
  fitlabelunits(2,22,3)='Ang       '
  fitlabelunits(3,22,3)='Ang       '
  fitlabelunits(4,22,3)='degrees   '
!
  nfitlabel(23,3) = 3
  fitlabel(1,23,3)='Three-body K   '
  fitlabel(2,23,3)='Three-body r0/1'
  fitlabel(3,23,3)='Three-body r0/2'
  fitlabelunits(1,23,3)='eV/Ang^m+n'
  fitlabelunits(2,23,3)='Ang       '
  fitlabelunits(3,23,3)='Ang       '
!
  nfitlabel(24,3) = 7
  fitlabel(1,24,3)='Three-body cnst'
  fitlabel(2,24,3)='Three-body angl'
  fitlabel(3,24,3)='Three-body A   '
  fitlabel(4,24,3)='Three-body B   '
  fitlabel(5,24,3)='Three-body C   '
  fitlabel(6,24,3)='Three-body D   '
  fitlabel(7,24,3)='Three-body E   '
  fitlabelunits(1,24,3)='eV*rad^-2 '
  fitlabelunits(2,24,3)='degrees   '
  fitlabelunits(3,24,3)='None      '
  fitlabelunits(4,24,3)='rad^-1    '
  fitlabelunits(5,24,3)='rad^-2    '
  fitlabelunits(6,24,3)='rad^-3    '
  fitlabelunits(7,24,3)='rad^-4    '
!
  nfitlabel(25,3) = 4
  fitlabel(1,25,3)='Three-body cnst'
  fitlabel(2,25,3)='Three-body angl'
  fitlabel(3,25,3)='Three-body rho1'
  fitlabel(4,25,3)='Three-body rho2'
  fitlabelunits(1,25,3)='eV        '
  fitlabelunits(2,25,3)='degrees   '
  fitlabelunits(3,25,3)='Ang       '
  fitlabelunits(4,25,3)='Ang       '
!
  nfitlabel(26,3) = 1
  fitlabel(1,26,3)='j3 J3          '
  fitlabelunits(1,26,3)='eV        '
!
  nfitlabel(27,3) = 5
  fitlabel(1,27,3)='Three-body cnst'
  fitlabel(3,27,3)='Three-bdy rmin1'
  fitlabel(4,27,3)='Three-bdy rmin2'
  fitlabel(5,27,3)='Three-body rho '
  fitlabelunits(1,27,3)='eV        '
  fitlabelunits(3,27,3)='Ang       '
  fitlabelunits(4,27,3)='Ang       '
  fitlabelunits(5,27,3)='Ang^-1    '
!------------------------
!  Four-body parameters -
!------------------------
  nfitlabel(1,4) = 1
  fitlabel(1,1,4)='Torsional const'
  fitlabelunits(1,1,4)='eV        '
!
  nfitlabel(2,4) = 6
  fitlabel(1,2,4)='Torsional const'
  fitlabel(2,2,4)='Ryckaert coeff1'
  fitlabel(3,2,4)='Ryckaert coeff2'
  fitlabel(4,2,4)='Ryckaert coeff3'
  fitlabel(5,2,4)='Ryckaert coeff4'
  fitlabel(6,2,4)='Ryckaert coeff5'
  fitlabelunits(1,2,4)='eV        '
  fitlabelunits(2,2,4)='eV        '
  fitlabelunits(3,2,4)='eV        '
  fitlabelunits(4,2,4)='eV        '
  fitlabelunits(5,2,4)='eV        '
  fitlabelunits(6,2,4)='eV        '
!
  nfitlabel(3,4) = 2
  fitlabel(1,3,4)='Out of plane K2'
  fitlabel(2,3,4)='Out of plane K4'
  fitlabelunits(1,3,4)='eV*Ang^-2 '
  fitlabelunits(2,3,4)='eV*Ang^-4 '
!
  nfitlabel(4,4) = 2
  fitlabel(1,4,4)='Torsional K1   '
  fitlabel(2,4,4)='Torsional K2   '
  fitlabelunits(1,4,4)='eV        '
  fitlabelunits(2,4,4)='eV        '
!
  nfitlabel(5,4) = 2
  fitlabel(1,5,4)='Torsional const'
  fitlabel(2,5,4)='Torsional phi0 '
  fitlabelunits(1,5,4)='eV*rad^-2 '
  fitlabelunits(2,5,4)='degrees   '
!
  nfitlabel(6,4) = 5
  fitlabel(1,6,4)='Torsional const'
  fitlabel(2,6,4)='Torsional phi0 '
  fitlabel(3,6,4)='Torsional rho12'
  fitlabel(4,6,4)='Torsional rho23'
  fitlabel(5,6,4)='Torsional rho34'
  fitlabelunits(1,6,4)='eV        '
  fitlabelunits(2,6,4)='degrees   '
  fitlabelunits(3,6,4)='Ang       '
  fitlabelunits(4,6,4)='Ang       '
  fitlabelunits(5,6,4)='Ang       '
!
  nfitlabel(7,4) = 5
  fitlabel(1,7,4)='Torsional K1   '
  fitlabel(2,7,4)='Torsional K2   '
  fitlabel(3,7,4)='Torsional rho12'
  fitlabel(4,7,4)='Torsional rho23'
  fitlabel(5,7,4)='Torsional rho34'
  fitlabelunits(1,7,4)='eV        '
  fitlabelunits(2,7,4)='eV        '
  fitlabelunits(3,7,4)='Ang       '
  fitlabelunits(4,7,4)='Ang       '
  fitlabelunits(5,7,4)='Ang       '
!
  nfitlabel(8,4) = 1
  fitlabel(1,8,4)='Torsional const'
  fitlabelunits(1,8,4)='eV        '
!
  nfitlabel(9,4) = 2
  fitlabel(1,9,4)='Torsional K1   '
  fitlabel(2,9,4)='Torsional K2   '
  fitlabelunits(1,9,4)='eV        '
  fitlabelunits(2,9,4)='eV        '
!
  nfitlabel(10,4) = 3
  fitlabel(1,10,4)='Torsional K    '
  fitlabel(2,10,4)='Torsion theta 1'
  fitlabel(3,10,4)='Torsion theta 2'
  fitlabelunits(1,10,4)='eV*rad^-2 '
  fitlabelunits(2,10,4)='degrees   '
  fitlabelunits(3,10,4)='degrees   '
!
  nfitlabel(11,4) = 1
  fitlabel(1,11,4)='Inversion K    '
  fitlabelunits(1,11,4)='eV        '
!
  nfitlabel(12,4) = 2
  fitlabel(1,12,4)='Inversion Sq K '
  fitlabel(2,12,4)='Inversion Sq K0'
  fitlabelunits(1,12,4)='eV        '
  fitlabelunits(2,12,4)='degrees   '
!
  nfitlabel(13,4) = 1
  fitlabel(1,13,4)='Torsional const'
  fitlabelunits(1,13,4)='eV        '
!
  nfitlabel(14,4) = 6
  fitlabel(1,14,4)='Xa-a K(213/4)  '
  fitlabel(2,14,4)='Xa-a K(312/4)  '
  fitlabel(3,14,4)='Xa-a K(412/3)  '
  fitlabel(4,14,4)='Xa-a Theta0_213'
  fitlabel(5,14,4)='Xa-a Theta0_214'
  fitlabel(6,14,4)='Xa-a Theta0_314'
  fitlabelunits(1,14,4)='eV        '
  fitlabelunits(2,14,4)='eV        '
  fitlabelunits(3,14,4)='eV        '
  fitlabelunits(4,14,4)='degrees   '
  fitlabelunits(5,14,4)='degrees   '
  fitlabelunits(6,14,4)='degrees   '
!
  nfitlabel(15,4) = 4
  fitlabel(1,15,4)='UFF oop K      '
  fitlabel(2,15,4)='UFF oop C0     '
  fitlabel(3,15,4)='UFF oop C1     '
  fitlabel(4,15,4)='UFF oop C2     '
  fitlabelunits(1,15,4)='eV        '
  fitlabelunits(2,15,4)='None      '
  fitlabelunits(3,15,4)='None      '
  fitlabelunits(4,15,4)='None      '
!
  nfitlabel(16,4) = 6
  fitlabel(1,16,4)='Xca-a K(213/4) '
  fitlabel(2,16,4)='Xca-a K(312/4) '
  fitlabel(3,16,4)='Xca-a K(412/3) '
  fitlabel(4,16,4)='Xca-a Theta_213'
  fitlabel(5,16,4)='Xca-a Theta_214'
  fitlabel(6,16,4)='Xca-a Theta_314'
  fitlabelunits(1,16,4)='eV        '
  fitlabelunits(2,16,4)='eV        '
  fitlabelunits(3,16,4)='eV        '
  fitlabelunits(4,16,4)='degrees   '
  fitlabelunits(5,16,4)='degrees   '
  fitlabelunits(6,16,4)='degrees   '
!
  nfitlabel(17,4) = 3
  fitlabel(1,17,4)='Torsional K    '
  fitlabel(2,17,4)='Torsion theta 1'
  fitlabel(3,17,4)='Torsion theta 2'
  fitlabelunits(1,17,4)='eV*rad^-2 '
  fitlabelunits(2,17,4)='degrees   '
  fitlabelunits(3,17,4)='degrees   '
!------------------------
!  Many-body parameters -
!------------------------
  fitlabel(1,1,5) ='EAM density p 1'
  fitlabel(2,1,5) ='EAM density p 2'
  fitlabel(3,1,5) ='EAM density p 3'
  fitlabel(4,1,5) ='EAM function p1'
  fitlabel(5,1,5) ='EAM function p2'
  fitlabel(6,1,5) ='EAM function p3'
  fitlabel(7,1,5) ='EAM alloy scale'
  fitlabel(8,1,5) ='EAM alloy add  '
  fitlabel(9,1,5) ='EAM function p4'
  fitlabel(10,1,5)='EAM function p5'
  fitlabel(11,1,5)='EAM function p6'
  fitlabel(12,1,5)='EAM function p7'
  fitlabel(13,1,5)='EAM function p8'
  fitlabel(14,1,5)='EAM function p9'
  fitlabel(15,1,5)='EAM density p 4'
  fitlabel(16,1,5)='EAM density p 5'
  fitlabel(17,1,5)='EAM density p 6'
  fitlabel(18,1,5)='EAM density p 7'
  fitlabel(19,1,5)='MEAM function t'
!-------------------------
!  Bond-order parameters -
!-------------------------
  fitlabel(1,1,6) ='B-O twobody A '
  fitlabel(2,1,6) ='B-O twobody B '
  fitlabel(3,1,6) ='B-O twobody Za'
  fitlabel(4,1,6) ='B-O twobody Zb'
  fitlabel(5,1,6) ='B-O rep Alpha '
  fitlabel(6,1,6) ='B-O att Alpha '
  fitlabel(7,1,6) ='B-O rep N     '
  fitlabel(8,1,6) ='B-O att N     '
  fitlabel(9,1,6) ='B-O rep Omega '
  fitlabel(10,1,6)='B-O rep Lambda'
  fitlabel(11,1,6)='B-O rep C / C1'
  fitlabel(12,1,6)='B-O rep D / C2'
  fitlabel(13,1,6)='B-O rep C3    '
  fitlabel(14,1,6)='B-O rep C4    '
  fitlabel(15,1,6)='B-O rep C5    '
  fitlabel(16,1,6)='B-O rep H     '
  fitlabel(17,1,6)='B-O att Omega '
  fitlabel(18,1,6)='B-O att Lambda'
  fitlabel(19,1,6)='B-O att C / C1'
  fitlabel(20,1,6)='B-O att D / C2'
  fitlabel(21,1,6)='B-O att C3    '
  fitlabel(22,1,6)='B-O att C4    '
  fitlabel(23,1,6)='B-O att C5    '
  fitlabel(24,1,6)='B-O att H     '
  fitlabel(25,1,6)='B-O rep chi   '
  fitlabel(26,1,6)='B-O att chi   '
  fitlabel(27,1,6)='B-O cn C1     '
  fitlabel(28,1,6)='B-O cn C2     '
  fitlabel(29,1,6)='B-O cn Z0     '
  fitlabel(30,1,6)='B-O cn E0     '
  fitlabelunits(1,1,6)='eV        '
  fitlabelunits(2,1,6)='eV        '
  fitlabelunits(3,1,6)='1/Ang     '
  fitlabelunits(4,1,6)='1/Ang     '
  fitlabelunits(7,1,6)='1/Ang     '
  fitlabelunits(16,1,6)='1/Ang     '
!-----------------------
!  Six-body parameters -
!-----------------------
  fitlabel(1,1,7)='XOut of plane K'
  fitlabelunits(1,1,7)='eV*Ang^-2 '
!----------------------
!  ReaxFF  parameters -
!----------------------
  fitlabel(1,1,8)='RFF1 rad sigma '
  fitlabel(2,1,8)='RFF1 rad pi    '
  fitlabel(3,1,8)='RFF1 rad pi-pi '
  fitlabel(1,2,8)='RFF1 val norm  '
  fitlabel(2,2,8)='RFF1 val boc   '
  fitlabel(3,2,8)='RFF1 val lp    '
  fitlabel(4,2,8)='RFF1 val angle '
  fitlabel(1,3,8)='RFF1 over boc3 '
  fitlabel(2,3,8)='RFF1 over boc4 '
  fitlabel(3,3,8)='RFF1 over boc5 '
  fitlabel(4,3,8)='RFF1 over ovun2'
  fitlabel(1,4,8)='RFF1 undr ovun5'
  fitlabel(1,5,8)='RFF1 lp nlp_opt'
  fitlabel(2,5,8)='RFF1 lp p_lp2  '
  fitlabel(1,6,8)='RFF1 angle val3'
  fitlabel(2,6,8)='RFF1 angle val5'
  fitlabel(1,7,8)='RFF1 morse alph'
  fitlabel(2,7,8)='RFF1 morse Dij '
  fitlabel(3,7,8)='RFF1 morse rvdw'
  fitlabel(4,7,8)='RFF1 morse gamm'
  fitlabel(1,8,8)='RFFq chi       '
  fitlabel(2,8,8)='RFFq mu        '
  fitlabel(3,8,8)='RFFq gamma     '
  fitlabel(4,8,8)='RFFq shell dMu '
  fitlabel(5,8,8)='RFFq shell Qs  '
  fitlabel(6,8,8)='RFFq shell beta'
  fitlabel(7,8,8)='RFFq fixed Q   '
  fitlabel(1,9,8) ='RFF0 bond boc1 '
  fitlabel(2,9,8) ='RFF0 bond boc2 '
  fitlabel(1,10,8)='RFF0 over ovun3'
  fitlabel(2,10,8)='RFF0 over ovun4'
  fitlabel(3,10,8)='RFF0 over ovun6'
  fitlabel(4,10,8)='RFF0 over ovun7'
  fitlabel(5,10,8)='RFF0 over ovun8'
  fitlabel(1,11,8)='RFF0 val  val6 '
  fitlabel(2,11,8)='RFF0 val  val8 '
  fitlabel(3,11,8)='RFF0 val  val9 '
  fitlabel(4,11,8)='RFF0 val  val10'
  fitlabel(1,12,8)='RFF0 pen  pen2 '
  fitlabel(2,12,8)='RFF0 pen  pen3 '
  fitlabel(3,12,8)='RFF0 pen  pen4 '
  fitlabel(1,13,8)='RFF0 tors tor2 '
  fitlabel(2,13,8)='RFF0 tors tor3 '
  fitlabel(3,13,8)='RFF0 tors tor4 '
  fitlabel(4,13,8)='RFF0 tors cot2 '
  fitlabel(1,14,8)='RFF0 vdw  vdw1 '
  fitlabel(1,15,8)='RFF0 lp  p_lp1 '
  fitlabel(1,16,8)='RFF2 bond De_s '
  fitlabel(2,16,8)='RFF2 bond De_p '
  fitlabel(3,16,8)='RFF2 bond De_pp'
  fitlabel(4,16,8)='RFF2 bond p_be1'
  fitlabel(5,16,8)='RFF2 bond p_be2'
  fitlabel(1,17,8)='RFF2 over ovun1'
  fitlabel(1,18,8)='RFF2 bo p_bo1  '
  fitlabel(2,18,8)='RFF2 bo p_bo2  '
  fitlabel(3,18,8)='RFF2 bo p_bo3  '
  fitlabel(4,18,8)='RFF2 bo p_bo4  '
  fitlabel(5,18,8)='RFF2 bo p_bo5  '
  fitlabel(6,18,8)='RFF2 bo p_bo6  '
  fitlabel(1,19,8)='RFF2 morse De  '
  fitlabel(2,19,8)='RFF2 morse alph'
  fitlabel(3,19,8)='RFF2 morse r0  '
  fitlabel(4,19,8)='RFF2 morse r_s '
  fitlabel(5,19,8)='RFF2 morse r_p '
  fitlabel(6,19,8)='RFF2 morse r_pp'
  fitlabel(1,20,8)='RFF2 pen p_pen1'
  fitlabel(2,20,8)='RFF2 pen p_pen2'
  fitlabel(3,20,8)='RFF2 pen p_pen3'
  fitlabel(1,21,8)='RFF3 ang theta0'
  fitlabel(2,21,8)='RFF3 ang p_val1'
  fitlabel(3,21,8)='RFF3 ang p_val2'
  fitlabel(4,21,8)='RFF3 ang p_val4'
  fitlabel(5,21,8)='RFF3 ang p_val7'
  fitlabel(6,21,8)='RFF3 ang p_val6'
  fitlabel(1,22,8)='RFF3 pen p_pen1'
  fitlabel(1,23,8)='RFF3 conj coa1 '
  fitlabel(2,23,8)='RFF3 conj coa2 '
  fitlabel(3,23,8)='RFF3 conj coa3 '
  fitlabel(4,23,8)='RFF3 conj coa4 '
  fitlabel(1,24,8)='RFF3 hb r0_hb  '
  fitlabel(2,24,8)='RFF3 hb p_hb1  '
  fitlabel(3,24,8)='RFF3 hb p_hb2  '
  fitlabel(4,24,8)='RFF3 hb p_hb3  '
  fitlabel(1,25,8)='RFF4 tor V1    '
  fitlabel(2,25,8)='RFF4 tor V2    '
  fitlabel(3,25,8)='RFF4 tor V3    '
  fitlabel(4,25,8)='RFF4 tor p_tor1'
  fitlabel(5,25,8)='RFF4 tor p_cot1'
!-------------------
!  EDIP parameters -
!-------------------
  fitlabel(1,1,9)='EDIP cn alpha  '
  fitlabel(2,1,9)='EDIP cn Zdih   '
  fitlabel(3,1,9)='EDIP cn Zrep   '
  fitlabel(4,1,9)='EDIP cn c0     '
  fitlabel(1,2,9)='EDIP 2b epsilon'
  fitlabel(2,2,9)='EDIP 2b B      '
  fitlabel(3,2,9)='EDIP 2b p      '
  fitlabel(4,2,9)='EDIP 2b beta   '
  fitlabel(5,2,9)='EDIP 2b sigma  '
  fitlabel(6,2,9)='EDIP 2b a      '
  fitlabel(7,2,9)='EDIP 2b a prime'
  fitlabel(1,3,9)='EDIP 3b lambda0'
  fitlabel(2,3,9)='EDIP 3b lambdap'
  fitlabel(3,3,9)='EDIP 3b Z0     '
  fitlabel(4,3,9)='EDIP 3b gamma0 '
  fitlabel(5,3,9)='EDIP 3b gammap '
  fitlabel(6,3,9)='EDIP 3b q      '
  fitlabel(7,3,9)='EDIP 3b kq2    '
#ifdef TRACE
  call trace_out('initial')
#endif
!
  return
  end
