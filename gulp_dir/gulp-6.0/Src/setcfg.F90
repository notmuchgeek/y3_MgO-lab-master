  subroutine setcfg
!
!  One off configuration set up
!
!  11/96 Make sure that maxd2 is no larger than 3*maxat+6 to
!        avoid wasting memory.
!   6/97 Define cutoff (costcut) for Pannetier Cost Function (SMW)
!   7/97 Don't print co-ords if called from predict{maxd2<>maxd2in}(SMW)
!   5/98 Comment out cutoff defn for costcut - common genconst (SMW)
!   5/98 Output of stress tensor added
!   8/99 Testing of angles for monoclinic case corrected -
!        only affects version 1.3
!  10/99 Bug in setting k points for multiple configurations fixed
!   7/00 lflags now in control module
!  11/00 maxd2 removed as handled elsewhere now
!  12/00 2-D modifications added
!  12/00 call to rlist added before geometry measurements for safety
!   6/01 Solvation model details now output
!   8/01 More solvation model details added to output
!   9/01 Default K point for FEM changed to a symmetric point.
!   9/01 Temperature converted to g format
!  11/01 Set up of growth slice using dhkl added and nzmol
!   8/02 Output of external forces added
!  11/02 Einstein model data output added
!  11/02 If Einstein model, then no atom is forced to be fixed.
!   5/03 Checking of memory for maxvar altered to minimise number of
!        allocations and deallocations
!   6/03 XML modifications added
!   9/03 Rigid region frozen direction flags now set
!  10/03 Rhombohedral coordinates for hexagonal system now output
!   2/04 Time-dependent force added
!   3/04 Spatial decomposition version of poccon added
!   4/04 Terse options added
!  10/04 Modifications for non-standard (>230) space groups made
!  12/04 Fixed atom choice for clusters and polymers improved
!   4/05 Shells excluded from fix atom search
!   8/05 Setting of fixed directions changed for nspg = 1, nccs > 0
!   2/06 Modified so that setting of a fixed direction discounts constrained atoms
!   5/06 Fitting to stresses added
!   9/06 Hiccup with fitting flags fixed in MC where no atoms are yet present
!   9/06 Call to formula changed to ioout
!  11/06 NEB modifications added
!  11/06 Checking of user input flags for symmetry correctness
!  11/06 Approach to flags modified. Default values of ltmp set first based 
!        either on user input or GULP defaults. Symmetry correctness then 
!        enforced. 
!  11/06 Total occupancy is zero warning added
!   2/07 Electric field output added
!   3/07 Radial force added
!   3/07 Calls to mxmb renamed to GULP_mxmb
!   3/07 Call to bond changed to GULP_bond
!   5/07 Output of QM or MM region type added
!   5/07 Output of QM/MM mode added
!   6/07 Forcing of first atom to be fixed turned off for MC if lflags
!   7/07 Flag setting modified so that if a plane potential is present then 
!        there is no default fixing in the z direction. 
!  12/07 Unused variables removed
!   1/08 Code that frees up first atom for MD excluded from acting in MC case
!   4/08 Logic for setting of lphonfit modified to avoid referencing ndimen(0)
!        due to reaction energies
!   5/08 Orthorhombic cell relaxation option added
!   6/08 Fixing of a single atom turned off for NEB/Sync and unfix option modified
!   8/08 Bug in unfixing atoms for MD on clusters corrected
!  10/08 COSMO/COSMIC changes merged in
!   1/09 Arrays apot, bpot, cpot & dpot replaced by new 2-D array twopot
!   1/09 Integer datatypes all explicitly declared
!   2/09 XML calls removed
!   3/09 lkptdispersion added
!   6/09 Module name changed from three to m_three
!   3/10 Modification of default weights added 
!   6/10 Number of decimal places for charges increased to 5.
!   8/10 lfix1atom introduced to indicate whether unfix keyword is specified or not
!  12/10 Check for missing shells + automatic addition added
!  12/10 Output of anisotropic pressure tensor added
!  12/10 Hide shell option added
!   6/11 Electric field enabled for periodic systems
!   7/11 Keyword added to prevent automatic adding of shells
!   7/11 Check on cell / space group consistency added
!   7/11 Output of constraints modified to allow for larger systems
!   8/11 Call to poccons now controlled by lspatialok and not lspatial
!  10/11 Output of electric field modified
!   7/12 C 1 & C -1 changed to triclinic from monoclinic
!  10/12 Call to setup KIM potentials added - needs configuration information
!  10/12 Call to angle and torsion modified to allow for new arguments being passed
!  11/12 Number of decimal places for charge output increased to 6
!  12/12 Time-dependent field added
!  12/12 Modified to allow for multiple time-dependent fields
!   7/13 Symmetry number added
!   8/13 Option to find space group added
!   9/13 nsuper changed to be a 2-D array
!  10/13 Referencing of lbsmat corrected for symmetry adapted case
!  12/13 Strains and stresses separated for fitting
!   3/14 Modified for harmonic relaxation
!   3/14 Call to angle changed to getangles for benefit of ChemShell
!   3/14 Call to torsion changed to gettorsions for benefit of ChemShell
!   3/14 Call to super changed to build_supercell for benefit of ChemShell
!  10/14 Output of force constant supercell added
!   4/15 Ghost supercell array added
!   7/15 Output of external potential added
!   3/17 fix_atom option added
!   3/17 Order of variables in iopt changed so that internals are first
!        and cell last for benefit of parallel code
!   3/17 lphonfit now set to be true if observables are present regardless
!        of dimensionality since phonon is being used in all cases
!   3/17 nobtyp = 29 added to the triggers for lphonfit
!   8/17 lodd freezing of atoms only applied if lsymopt is true since 
!        otherwise fixing the first atom is incompatible with nfixatom
!   1/18 Position of default k point setting moved so that k points are
!        output correctly.
!   1/18 Default free energy k point changed to 1/4,1/4,1/4
!   1/18 Trace added
!   2/18 nxks, nyks, nzks converted to a single array
!   2/18 Saving of primitive cell added
!   4/18 Twist option added
!   4/18 Format for number of configurations changed to allow for larger numbers
!   6/18 Strain cell option added
!   9/18 Modified so that both reference and actual cells are output for 
!        lstraincell algorithm
!  11/18 Output corrected for reference cells
!  12/18 Shear force added
!   2/19 x0 removed
!   2/19 Format for temperature changed to satisfy Intel compiler
!   3/19 iopt replaced by ioptindex and iopttype
!   3/19 ltmp changed for lopt for each type of variable
!   3/19 Constraint arrays changed to have index and type
!   8/19 Rigid molecules added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!   9/19 rvpcfg only set if not read in
!  12/19 Corrections to flag handling for rigid molecules
!   2/20 Call to specialmol added for rigid molecule case
!   2/20 Breathing shell radii now frozen in rigid regions
!   2/20 lksorigin added
!   3/20 Configuration specific flag added to indicate whether shell gradients
!        should be added to those from the core
!   3/20 Dielectric constant added
!   3/20 Rotation of rigid molecules excluded in rigid regions
!   4/20 Rigid molecule modifications added
!   5/20 Output of rigid molecule constraints added
!   5/20 Variable order changed for rigid molecules to be all 
!        translations before all quaternions to help tmat algorithm
!   6/20 Handling of shells in setting nfixatom corrected
!   6/20 Corrections to use of nasymnomol with nasymnomolptr
!   7/20 Handling of constraints for cell corrected for ocell option
!   7/20 Algorithm for fixing molecules in rigid regions changed
!   7/20 Check on calling setregiontrans corrected
!   7/20 Missing breathing radius variables set
!   7/20 Call to symupdatemol added prior to specialcom to setup molQsym
!   7/20 Rigid molecule fitting added
!  10/20 lterseinmol added
!  12/20 Initialisation of optimisation logical arrays added to ensure
!        consistent behaviour for cases where optimisation is not invoked
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
  use control
  use cosmic
  use current
  use derivatives,   only : lfcsupercell
  use dispersion
  use element
  use field
  use fitting,       only : lsumcoreshell
  use four
  use freeze
  use general
  use g_constants,   only : pi
  use iochannels
  use kim_models,    only : lkim_model
  use ksample
  use m_three
  use mdlogic
  use moldyn,        only : lfix
  use molecule
  use parallel
  use plane,         only : nplanepot
  use observables
  use optimisation
  use radial
  use shells
  use spatial,       only : lspatialok
  use symmetry
  use terse
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  use two
  use wolfcosmo,     only : etawc, cutwc
  implicit none
!
!  Local variables
!
  character(len=4)                         :: crd4(9)
  character(len=1)                         :: crd1(3)
  character(len=1)                         :: ocha(3)
  character(len=2)                         :: cstype
  character(len=3)                         :: fixed
  character(len=3)                         :: fixstring
  character(len=5)                         :: lab
  character(len=7)                         :: systype(4)
  integer(i4)                              :: i
  integer(i4)                              :: ic
  integer(i4)                              :: icfg
  integer(i4)                              :: ii
  integer(i4)                              :: inat
  integer(i4)                              :: ind
  integer(i4)                              :: itp
  integer(i4)                              :: itype
  integer(i4)                              :: ix
  integer(i4)                              :: iy
  integer(i4)                              :: iz
  integer(i4)                              :: j
  integer(i4)                              :: jj
  integer(i4)                              :: k
  integer(i4)                              :: na
  integer(i4)                              :: nangtot
  integer(i4)                              :: nati
  integer(i4)                              :: ncm
  integer(i4)                              :: ncrf
  integer(i4)                              :: ncrv
  integer(i4)                              :: ncvi
  integer(i4)                              :: ndir
  integer(i4)                              :: ndirp(3)
  integer(i4)                              :: nfv
  integer(i4)                              :: ngroup
  integer(i4)                              :: nin
  integer(i4)                              :: nk
  integer(i4)                              :: nkp
  integer(i4)                              :: nlfgra
  integer(i4)                              :: nlfstr
  integer(i4)                              :: nlfgrad
  integer(i4)                              :: nlfstrain
  integer(i4)                              :: nlkpt
  integer(i4)                              :: nobsold
  integer(i4)                              :: np
  integer(i4)                              :: npt
  integer(i4)                              :: nr
  integer(i4)                              :: nri
  integer(i4)                              :: nrj
  integer(i4)                              :: nspg
  integer(i4)                              :: nphitot
  integer(i4)                              :: nufgra
  integer(i4)                              :: nufstr
  integer(i4)                              :: nvarl
  integer(i4)                              :: nvarlold
  integer(i4)                              :: status
  logical                                  :: loptc(6)    ! Optimisation flag - cell
  logical, dimension(:,:), allocatable     :: loptm       ! Optimisation flag - molecules
  logical, dimension(:),   allocatable     :: loptr       ! Optimisation flag - radius
  logical, dimension(:),   allocatable     :: loptx       ! Optimisation flag - x
  logical, dimension(:),   allocatable     :: lopty       ! Optimisation flag - y
  logical, dimension(:),   allocatable     :: loptz       ! Optimisation flag - z
  logical, dimension(:),   allocatable     :: loptxloc
  logical, dimension(:),   allocatable     :: loptyloc
  logical, dimension(:),   allocatable     :: loptzloc
  logical                                  :: lalleinstein
  logical                                  :: lbreathe
  logical                                  :: lcore
  logical                                  :: lfirstout
  logical                                  :: lfixeddirection
  logical                                  :: lfound
  logical                                  :: lfound1
  logical                                  :: lfound2
  logical                                  :: liso
  logical                                  :: lnoflagsloc
  logical                                  :: lnxi
  logical                                  :: lnyi
  logical                                  :: lnzi
  logical                                  :: lodd
  logical                                  :: lortho
  logical                                  :: lphonfit
  logical                                  :: lrhombo
  logical                                  :: lshrink
  logical                                  :: lslice
  logical                                  :: lsuper
  logical                                  :: lsuperghost
  real(dp)                                 :: afull
  real(dp)                                 :: alphafull
  real(dp)                                 :: alpprim
  real(dp)                                 :: aprim
  real(dp)                                 :: ara
  real(dp)                                 :: area
  real(dp)                                 :: betafull
  real(dp)                                 :: betprim
  real(dp)                                 :: bfull
  real(dp)                                 :: bprim
  real(dp)                                 :: cfull
  real(dp)                                 :: cprim
  real(dp)                                 :: dipolex
  real(dp)                                 :: dipoley
  real(dp)                                 :: dipolez
  real(dp)                                 :: fieldnorm
  real(dp)                                 :: forcenorm
  real(dp)                                 :: gamprim
  real(dp)                                 :: gammafull
  real(dp)                                 :: occtot
  real(dp)                                 :: r2
  real(dp)                                 :: r2best
  real(dp)                                 :: ra
  real(dp)                                 :: rvt(3,3)
  real(dp)                                 :: sum
  real(dp)                                 :: sumx
  real(dp)                                 :: sumy
  real(dp)                                 :: sumz
  real(dp)                                 :: vol
  real(dp)                                 :: volp
  real(dp)                                 :: volume
  real(dp)                                 :: x(3)
  real(dp)                                 :: xx(3)
  real(dp)                                 :: zbest
  real(dp)                                 :: xc
  real(dp)                                 :: yc
  real(dp)                                 :: zc
  real(dp)                                 :: xmid
  real(dp)                                 :: ymid
  real(dp)                                 :: zmid
!
  data crd1/'x','y','z'/
  data crd4/'x   ','y   ','z   ','xcom','ycom','zcom','xqtn','yqtn','zqtn'/
  data systype/'Cluster','Polymer','Surface','Bulk   '/
#ifdef TRACE
  call trace_in('setcfg')
#endif
!
  lnoflagsloc = lnoflags
  if (.not.lnoflagsloc) then
    lnoflagsloc = (.not.lconp.and..not.lconv.and..not.lcello.and..not.lshello)
  endif
  liso = (index(keyword,'iso').eq.1.or.index(keyword,' iso').ne.0) 
  lortho = (index(keyword,'ort').eq.1.or.index(keyword,' ort').ne.0) 
  lbreathe = (index(keyword,' brea').ne.0.or.index(keyword,'brea').eq.1)
!
!  Store number of constraints read in - if none have been read in then there is no need 
!  to dump constraints as they were automatically generated by GULP
!
  nconin = ncontot
!
!  Check if phonon fitting is required
!
  lphonfit = .false.
  if (lfit) then
    do i = 1,nobs
      if (nobtyp(i).eq.9.or.nobtyp(i).eq.13.or.nobtyp(i).eq.14.or.nobtyp(i).eq.29) then
        lphonfit = .true.
      endif
    enddo
  endif
!******************************************
!  Output total number of configurations  *
!******************************************
  if (ioproc) then
    write(ioout,'(/,''  Total number of configurations input = '',i6)') ncfg
  endif
!***************************************************
!  Check for missing shells that need to be added  *
!***************************************************
  if (index(keyword,' noad').eq.0.and.index(keyword,'noad').ne.1) then
    do i = 1,ncfg
      call setup(.false.)
      call addshell
    enddo
  endif
!*********************************************************
!  Transform configurations input in rhombohedral form   *
!  hexagonal setting (ifhr=1)                            *
!*********************************************************
  do i = 1,ncfg
    if (ndimen(i).eq.3) then
      ncf = i
      call setup(.false.)
      nspg = nspcg(i)
!
!  Option to find symmetry
!
!      if (lfindsym) then
!        call findspacegroup(ngroup)
!      endif
!
      if (nspg.le.2) then
        ictype = 1
      elseif (nspg.ge.3.and.nspg.le.15) then
        ictype = 2
      elseif (nspg.ge.16.and.nspg.le.74) then
        ictype = 3
      elseif (nspg.ge.75.and.nspg.le.142) then
        ictype = 4
      elseif (nspg.ge.143.and.nspg.le.194) then
        ictype = 5
        call cellfhr(icfhr,rv)
        ifhr(i) = icfhr
      elseif (nspg.ge.195.and.nspg.le.230) then
        ictype = 6
      elseif (nspg.ge.231.and.nspg.le.232) then
        ictype = 1
      endif
      call setup(.false.)
      if (nccs.ne.5) ifhr(i) = 0
      if (ifhr(i).eq.1) then
!
!  Change cell to hexgonal form
!
        call rhtohex
!
!  Change fractional coordinates
!
        do na = 1,nasym
          x(1) = xcfg(nsft+na)
          x(2) = ycfg(nsft+na)
          x(3) = zcfg(nsft+na)
          xx(1) = 0.0_dp
          xx(2) = 0.0_dp
          xx(3) = 0.0_dp
          call GULP_mxmb(w(ncbl,1,1),7_i4,21_i4,x,1_i4,1_i4,xx,1_i4,1_i4,3_i4,3_i4,1_i4)
!
!  Place fractional coordinates in range 0-1
!
          xx(1) = xx(1) + 3.0_dp
          nin = xx(1)
          xx(1) = xx(1) - nin
          xx(2) = xx(2) + 3.0_dp
          nin = xx(2)
          xx(2) = xx(2) - nin
          xx(3) = xx(3) + 3.0_dp
          nin = xx(3)
          xx(3) = xx(3) - nin
          xcfg(nsft+na) = xx(1)
          ycfg(nsft+na) = xx(2)
          zcfg(nsft+na) = xx(3)
        enddo
      endif
    endif
  enddo
!******************************
!  Twist option for polymers  *
!******************************
  do i = 1,ncfg
    if (ndimen(i).eq.1.and.ntwistcfg(i).ne.0_i4) then
      ncf = i
      call setup(.false.)
!
!  Loop over atoms to find central axis
!
      ymid = 0.0_dp
      zmid = 0.0_dp
      do na = 1,nasym
        ymid = ymid + yalat(na)
        zmid = zmid + zalat(na)
      enddo
      ymid = ymid/dble(nasym)
      zmid = zmid/dble(nasym)
!
!  Loop over atoms rotating them around the periodic x axis
!
      ra = 2.0_dp*dble(ntwistcfg(i))*pi/a
      do na = 1,nasym
        xc = xalat(na)
        yc = yalat(na) - ymid
        zc = zalat(na) - zmid
        r2 = yc*yc + zc*zc
        r2 = sqrt(r2)
        ycfg(nsft+na) = ymid + yc*cos(ra*xc) - zc*sin(ra*xc)
        zcfg(nsft+na) = zmid + yc*sin(ra*xc) + zc*cos(ra*xc)
      enddo
!
!  Set twist to zero to avoid being repeated
!
      ntwistcfg(i) = 0_i4
    endif
  enddo
!**************************************************
!  One off initial set up for each configuration  *
!**************************************************
!
!  (1) Call setup to do symmetry
!  (2) Centre cell if necessary
!  (3) Create optimisation pointer arrays ioptindexcfg and iopttypecfg
!
  nvarl = 0
  do icfg = 1,ncfg
    ncf = icfg
    nspg = max(nspcg(icfg),1)
    lsuper = (nsuper(1,icfg).gt.1.or.nsuper(2,icfg).gt.1.or.nsuper(3,icfg).gt.1)
    lsuperghost = (nsuperghost(1,icfg).gt.1.or.nsuperghost(2,icfg).gt.1.or.nsuperghost(3,icfg).gt.1)
    lsymopt = lsymset(icfg)
    call setup(.false.)
    if (ndimen(icfg).eq.3) then
      if (lsymopt) call centre
      do i = 1,3
        rv(1,i) = rvcfg(1,i,icfg)
        rv(2,i) = rvcfg(2,i,icfg)
        rv(3,i) = rvcfg(3,i,icfg)
      enddo
!
!  Save primitive cell in case needed later
!
      if (.not.lvecpin(icfg)) then
        ncorepcfg(icfg) = ncore
        do i = 1,3
          rvpcfg(1:3,i,icfg) = rvcfg(1:3,i,icfg)
        enddo
      endif
!
      do i = 1,3
        rvt(1,i) = rv(1,i)
        rvt(2,i) = rv(2,i)
        rvt(3,i) = rv(3,i)
      enddo
      call uncentre(rvt)
      call uncell3D(rvt,afull,bfull,cfull,alphafull,betafull,gammafull)
    endif
!
!  If symmetry is to be switched off then expand structure due to nosym option or supercell calcn
!
    if ((lsymopt.and..not.lsym).or.lsuper.or.lmc.or.lmd) then
!
!  Call setup to reinitialise qf()
!
      call symoff
      call setup(.false.)
      lsym = .false.
      lsymopt = .false.
      lsymset(icfg) = .false.
      ictype = 1
      if (lsuper) call build_supercell
    endif
!
!  Generate coordinates here in case molecule calculation is to be performed
!
    if (ndimen(icfg).eq.3) then
!
!  3-D
!
      do i = 1,numat
        xclat(i) = xfrac(i)*rv(1,1) + yfrac(i)*rv(1,2) + zfrac(i)*rv(1,3)
        yclat(i) = xfrac(i)*rv(2,1) + yfrac(i)*rv(2,2) + zfrac(i)*rv(2,3)
        zclat(i) = xfrac(i)*rv(3,1) + yfrac(i)*rv(3,2) + zfrac(i)*rv(3,3)
      enddo
    elseif (ndimen(icfg).eq.2) then
!
!  2-D
!
      do i = 1,numat
        xclat(i) = xfrac(i)*rv(1,1) + yfrac(i)*rv(1,2)
        yclat(i) = xfrac(i)*rv(2,1) + yfrac(i)*rv(2,2)
        zclat(i) = zfrac(i)
      enddo
    endif
    call setup(.true.)
!
!  Count number of breathing shells and set default radii
!
    nbsm = 0
    do i = 1,nasym
      if (lbsmat(nsft+i)) then
        nbsm = nbsm + 1
        if (radcfg(nsft+i).eq.0.0_dp) then
!
!  Find potential and set radius equal to equilibrium value
!
          lfound = .false.
          np = 0
          do while (np.lt.npote.and..not.lfound)
            np = np + 1
            npt = nptype(np)
            lfound = ((npt.eq.14.or.npt.eq.17.or.npt.eq.31).and.iatn(i).eq.nspec1(np).and. &
              (natype(i).eq.nptyp1(np).or.nptyp1(np).eq.0))
          enddo
          if (lfound) then
            radcfg(nsft+i) = twopot(2,np)
          else
            nati = iatn(i)
            if (nati.gt.maxele) nati = nati - maxele
            radcfg(nsft+i) = rion(nati)
          endif
        endif
      endif
    enddo
!
!  Set flag according to whether all atoms are fixed to sites by Einstein model
!
    lalleinstein = .true.
    do i = 1,nasym
      if (.not.leinsteinat(nsft + i)) lalleinstein = .false.
    enddo
!
    if (nspg.le.2) then
      ictype = 1
    elseif (nspg.ge.3.and.nspg.le.15) then
      ictype = 2
    elseif (nspg.ge.16.and.nspg.le.74) then
      ictype = 3
    elseif (nspg.ge.75.and.nspg.le.142) then
      ictype = 4
    elseif (nspg.ge.143.and.nspg.le.194) then
      ictype = 5
    elseif (nspg.ge.195.and.nspg.le.230) then
      ictype = 6
    elseif (nspg.ge.231.and.nspg.le.232) then
      ictype = 1
    endif
    lrhombo = (ifhr(ncf).eq.1.and.(.not.lhex))
!
!  Set up growth slice if necessary
!
    if (ndimen(icfg).eq.2) then
      if (dhklcfg(icfg).gt.0.0_dp) then
        call setslice(icfg)
      endif
      call setzmolslice
    endif
!
!  Allocate array for optimisation logicals
!
    allocate(loptm(6_i4,nmolasym),stat=status)
    if (status/=0) call outofmemory('setcfg','loptm')
    allocate(loptr(nasym),stat=status)
    if (status/=0) call outofmemory('setcfg','loptr')
    allocate(loptx(nasym),stat=status)
    if (status/=0) call outofmemory('setcfg','loptx')
    allocate(lopty(nasym),stat=status)
    if (status/=0) call outofmemory('setcfg','lopty')
    allocate(loptz(nasym),stat=status)
    if (status/=0) call outofmemory('setcfg','loptz')
!
!  Initial arrays 
!
    loptc(1:6) = .false.
    loptm(1:6,1:nmolasym) = .false.
    loptr(1:nasym) = .false.
    loptx(1:nasym) = .false.
    lopty(1:nasym) = .false.
    loptz(1:nasym) = .false.
!
    if ((lopt.or.lharmrelax.or.lgrad.or.lfit.or.lmc.or.lmd.or.lneb).and..not.lbulknoopt) then
!*******************
!  Unit cell flags *
!*******************
      if (ndimen(icfg).gt.0) then
!
!  Set initial flags prior to symmetry checking
!
        if (lflags) then
          ind = 6*(icfg-1)
          do j = 1,nstrains
            loptc(j) = lopfc(ind+j)
          enddo
        else
          if (lconp.or.lcello) then
            if (liso) then
!
!  Isotropic cell expansion only
!
              loptc(1) = .true.
              do j = 2,nstrains
                loptc(j) = .false.
              enddo
            elseif (lortho) then
!
!  Orthorhombic cell expansion only
!
              if (ndimen(icfg).eq.3) then
                loptc(1) = .true.
                loptc(2) = .true.
                loptc(3) = .true.
                do j = 4,nstrains
                  loptc(j) = .false.
                enddo
              elseif (ndimen(icfg).eq.2) then
                loptc(1) = .true.
                loptc(2) = .true.
                loptc(3) = .false.
              elseif (ndimen(icfg).eq.1) then
                loptc(1) = .true.
              endif
            else
!
!  Anisotropic cell expansion
!
              do j = 1,nstrains
                loptc(j) = .true.
              enddo
            endif
          else
            do j = 1,nstrains
              loptc(j) = .false.
            enddo
          endif
        endif
        if (liso) then
!
!  Isotropic cell expansion only
!
          if (ndim.eq.3) then
            if (ncontot+2.ge.maxcontot) then
              maxcontot = ncontot + 50
              call changemaxcontot
            endif 
            if (icfg.lt.ncfg) then
              do k = ncontot,n1con(icfg+1),-1
                ncvarindcfg(k+2) = ncvarindcfg(k)
                ncvartypcfg(k+2) = ncvartypcfg(k)
                ncfixindcfg(k+2) = ncfixindcfg(k)
                ncfixtypcfg(k+2) = ncfixtypcfg(k)
                concocfg(k+2) = concocfg(k)
                nconcfg(k+2) = nconcfg(k)
                conaddcfg(k+2) = conaddcfg(k)
              enddo
              do k = icfg+1,ncfg
                n1con(k) = n1con(k) + 2
              enddo
            endif
            ncontot = ncontot + 2
            ncon = ncon + 1
            ncvarindcfg(ncfst+ncon) = 1
            ncfixindcfg(ncfst+ncon) = 2
            if (loptcellpar) then
              ncvartypcfg(ncfst+ncon) = iopt_cell
              ncfixtypcfg(ncfst+ncon) = iopt_cell
            else
              ncvartypcfg(ncfst+ncon) = iopt_strain
              ncfixtypcfg(ncfst+ncon) = iopt_strain
            endif
            concocfg(ncfst+ncon) = 1.0_dp
            conaddcfg(ncfst+ncon) = 0.0_dp
            nconcfg(ncfst+ncon) = icfg
            ncon = ncon + 1
            ncvarindcfg(ncfst+ncon) = 1
            ncfixindcfg(ncfst+ncon) = 3
            if (loptcellpar) then
              ncvartypcfg(ncfst+ncon) = iopt_cell
              ncfixtypcfg(ncfst+ncon) = iopt_cell
            else
              ncvartypcfg(ncfst+ncon) = iopt_strain
              ncfixtypcfg(ncfst+ncon) = iopt_strain
            endif
            concocfg(ncfst+ncon) = 1.0_dp
            conaddcfg(ncfst+ncon) = 0.0_dp
            nconcfg(ncfst+ncon) = icfg
          elseif (ndim.eq.2) then
            if (ncontot+1.ge.maxcontot) then
              maxcontot = ncontot + 50
              call changemaxcontot
            endif 
            if (icfg.lt.ncfg) then
              do k = ncontot,n1con(icfg+1),-1
                ncvarindcfg(k+1) = ncvarindcfg(k)
                ncvartypcfg(k+1) = ncvartypcfg(k)
                ncfixindcfg(k+1) = ncfixindcfg(k)
                ncfixtypcfg(k+1) = ncfixtypcfg(k)
                concocfg(k+1) = concocfg(k)
                nconcfg(k+1) = nconcfg(k)
                conaddcfg(k+1) = conaddcfg(k)
              enddo
              do k = icfg+1,ncfg
                n1con(k) = n1con(k) + 1
              enddo
            endif
            ncontot = ncontot + 1
            ncon = ncon + 1
            ncvarindcfg(ncfst+ncon) = 1
            ncfixindcfg(ncfst+ncon) = 2
            if (loptcellpar) then
              ncvartypcfg(ncfst+ncon) = iopt_cell
              ncfixtypcfg(ncfst+ncon) = iopt_cell
            else
              ncvartypcfg(ncfst+ncon) = iopt_strain
              ncfixtypcfg(ncfst+ncon) = iopt_strain
            endif
            concocfg(ncfst+ncon) = 1.0_dp
            conaddcfg(ncfst+ncon) = 0.0_dp
            nconcfg(ncfst+ncon) = icfg
          endif
        else
!
!  Anisotropic cell expansion
!
          if (lsymopt.and.nspg.gt.1) then
!
!  For certain cell types remove redundant cell strains
!
            if (ictype.gt.2) then
!
!  Orthorhombic, tetragonal, hexagonal, trigonal and cubic
!
              loptc(4) = .false.
              loptc(5) = .false.
              loptc(6) = .false.
              if (.not.lvecin(ncf)) then
                if (ictype.eq.6) then
                  loptc(2) = .false.
                  loptc(3) = .false.
!
!  Check to see if constraint already exists
!
                  lfound1 = .false.
                  lfound2 = .false.
                  do k = 1,ncon
                    if (ncvartypcfg(k+ncfst).eq.iopt_strain.or.ncvartypcfg(k+ncfst).eq.iopt_cell) then
                      if (ncfixtypcfg(k+ncfst).eq.iopt_strain.or.ncfixtypcfg(k+ncfst).eq.iopt_cell) then
                        if (ncvarindcfg(k+ncfst).eq.1.and.ncfixindcfg(k+ncfst).eq.2) then 
                          lfound1 = .true.
                        elseif (ncvarindcfg(k+ncfst).eq.1.and.ncfixindcfg(k+ncfst).eq.3) then
                          lfound2 = .true.
                        endif
                      endif
                    endif
                  enddo
                  if (.not.lfound1.and..not.lfound2) then
                    if (ncontot+2.ge.maxcontot) then
                      maxcontot = ncontot + 50
                      call changemaxcontot
                    endif 
                    if (icfg.lt.ncfg) then
                      do k = ncontot,n1con(icfg+1),-1
                        ncvarindcfg(k+2) = ncvarindcfg(k)
                        ncvartypcfg(k+2) = ncvartypcfg(k)
                        ncfixindcfg(k+2) = ncfixindcfg(k)
                        ncfixtypcfg(k+2) = ncfixtypcfg(k)
                        concocfg(k+2) = concocfg(k)
                        conaddcfg(k+2) = conaddcfg(k)
                        nconcfg(k+2) = nconcfg(k)
                      enddo
                      do k = icfg+1,ncfg
                        n1con(k) = n1con(k) + 2
                      enddo
                    endif
                    ncontot = ncontot + 2
                    ncon = ncon + 1
                    ncvarindcfg(ncfst+ncon) = 1
                    ncfixindcfg(ncfst+ncon) = 2
                    if (loptcellpar) then
                      ncvartypcfg(ncfst+ncon) = iopt_cell
                      ncfixtypcfg(ncfst+ncon) = iopt_cell
                    else
                      ncvartypcfg(ncfst+ncon) = iopt_strain
                      ncfixtypcfg(ncfst+ncon) = iopt_strain
                    endif
                    concocfg(ncfst+ncon) = 1.0_dp
                    conaddcfg(ncfst+ncon) = 0.0_dp
                    nconcfg(ncfst+ncon) = icfg
                    ncon = ncon + 1
                    ncvarindcfg(ncfst+ncon) = 1
                    ncfixindcfg(ncfst+ncon) = 3
                    if (loptcellpar) then
                      ncvartypcfg(ncfst+ncon) = iopt_cell
                      ncfixtypcfg(ncfst+ncon) = iopt_cell
                    else
                      ncvartypcfg(ncfst+ncon) = iopt_strain
                      ncfixtypcfg(ncfst+ncon) = iopt_strain
                    endif
                    concocfg(ncfst+ncon) = 1.0_dp
                    conaddcfg(ncfst+ncon) = 0.0_dp
                    nconcfg(ncfst+ncon) = icfg
                  elseif (lfound1.and..not.lfound2) then
                    if (ncontot.ge.maxcontot) then
                      maxcontot = ncontot + 50
                      call changemaxcontot
                    endif 
                    if (icfg.lt.ncfg) then
                      do k = ncontot,n1con(icfg+1),-1
                        ncvarindcfg(k+1) = ncvarindcfg(k)
                        ncvartypcfg(k+1) = ncvartypcfg(k)
                        ncfixindcfg(k+1) = ncfixindcfg(k)
                        ncfixtypcfg(k+1) = ncfixtypcfg(k)
                        concocfg(k+1) = concocfg(k)
                        conaddcfg(k+1) = conaddcfg(k)
                        nconcfg(k+1) = nconcfg(k)
                      enddo
                      do k = icfg+1,ncfg
                        n1con(k) = n1con(k) + 1
                      enddo
                    endif
                    ncontot = ncontot + 1
                    ncon = ncon + 1
                    ncvarindcfg(ncfst+ncon) = 1
                    ncfixindcfg(ncfst+ncon) = 3
                    if (loptcellpar) then
                      ncvartypcfg(ncfst+ncon) = iopt_cell
                      ncfixtypcfg(ncfst+ncon) = iopt_cell
                    else
                      ncvartypcfg(ncfst+ncon) = iopt_strain
                      ncfixtypcfg(ncfst+ncon) = iopt_strain
                    endif
                    concocfg(ncfst+ncon) = 1.0_dp
                    conaddcfg(ncfst+ncon) = 0.0_dp
                    nconcfg(ncfst+ncon) = icfg
                  elseif (.not.lfound1.and.lfound2) then
                    if (ncontot.ge.maxcontot) then
                      maxcontot = ncontot + 50
                      call changemaxcontot
                    endif 
                    if (icfg.lt.ncfg) then
                      do k = ncontot,n1con(icfg+1),-1
                        ncvarindcfg(k+1) = ncvarindcfg(k)
                        ncvartypcfg(k+1) = ncvartypcfg(k)
                        ncfixindcfg(k+1) = ncfixindcfg(k)
                        ncfixtypcfg(k+1) = ncfixtypcfg(k)
                        concocfg(k+1) = concocfg(k)
                        nconcfg(k+1) = nconcfg(k)
                        conaddcfg(k+1) = conaddcfg(k)
                      enddo
                      do k = icfg+1,ncfg
                        n1con(k) = n1con(k) + 1
                      enddo
                    endif
                    ncontot = ncontot + 1
                    ncon = ncon + 1
                    ncvarindcfg(ncfst+ncon) = 1
                    ncfixindcfg(ncfst+ncon) = 2
                    if (loptcellpar) then
                      ncvartypcfg(ncfst+ncon) = iopt_cell
                      ncfixtypcfg(ncfst+ncon) = iopt_cell
                    else
                      ncvartypcfg(ncfst+ncon) = iopt_strain
                      ncfixtypcfg(ncfst+ncon) = iopt_strain
                    endif
                    concocfg(ncfst+ncon) = 1.0_dp
                    conaddcfg(ncfst+ncon) = 0.0_dp
                    nconcfg(ncfst+ncon) = icfg
                  endif
                elseif (ictype.eq.4.or.ictype.eq.5) then
                  loptc(2) = .false.
!
!  Check to see if constraint already exists
!
                  lfound1 = .false.
                  do k = 1,ncon
                    if (ncvartypcfg(k+ncfst).eq.iopt_strain.or.ncvartypcfg(k+ncfst).eq.iopt_cell) then
                      if (ncfixtypcfg(k+ncfst).eq.iopt_strain.or.ncfixtypcfg(k+ncfst).eq.iopt_cell) then
                        if (ncvarindcfg(k+ncfst).eq.1.and.ncfixindcfg(k+ncfst).eq.2) then 
                          lfound1 = .true.
                        endif
                      endif
                    endif
                  enddo
                  if (.not.lfound1) then
                    if (ncontot.ge.maxcontot) then
                      maxcontot = ncontot + 50
                      call changemaxcontot
                    endif 
                    if (icfg.lt.ncfg) then
                      do k = ncontot,n1con(icfg+1),-1
                        ncvarindcfg(k+1) = ncvarindcfg(k)
                        ncvartypcfg(k+1) = ncvartypcfg(k)
                        ncfixindcfg(k+1) = ncfixindcfg(k)
                        ncfixtypcfg(k+1) = ncfixtypcfg(k)
                        concocfg(k+1) = concocfg(k)
                        conaddcfg(k+1) = conaddcfg(k)
                        nconcfg(k+1) = nconcfg(k)
                      enddo
                      do k = icfg+1,ncfg
                        n1con(k) = n1con(k) + 1
                      enddo
                    endif
                    ncontot = ncontot + 1
                    ncon = ncon + 1
                    ncvarindcfg(ncfst+ncon) = 1
                    ncfixindcfg(ncfst+ncon) = 2
                    if (loptcellpar) then
                      ncvartypcfg(ncfst+ncon) = iopt_cell
                      ncfixtypcfg(ncfst+ncon) = iopt_cell
                    else
                      ncvartypcfg(ncfst+ncon) = iopt_strain
                      ncfixtypcfg(ncfst+ncon) = iopt_strain
                    endif
                    concocfg(ncfst+ncon) = 1.0_dp
                    conaddcfg(ncfst+ncon) = 0.0_dp
                    nconcfg(ncfst+ncon) = icfg
                  endif
                endif
              endif
            elseif (ictype.eq.2) then
!
!  Monoclinic
!
              if (abs(alphafull-90.0_dp).gt.1.0d-4) then
                loptc(5) = .false.
                loptc(6) = .false.
              elseif (abs(gammafull-90.0_dp).gt.1.0d-4) then
                loptc(4) = .false.
                loptc(5) = .false.
              else
                loptc(4) = .false.
                loptc(6) = .false.
              endif
            endif
          endif
        endif
      endif
!******************
!  Internal flags *
!******************
!
!  Set initial flags based on keywords or lack of them
!
      if (lflags) then
!
!  Atoms
!
        do j = 1,nasymnomol
          jj = nasymnomolptr(j)
          loptx(j) = lopfi(3*(jj+nsft-1)+1)
          lopty(j) = lopfi(3*(jj+nsft-1)+2)
          loptz(j) = lopfi(3*(jj+nsft-1)+3)
        enddo
        if (lrigid) then
!
!  Rigid molecules
!
          do j = 1,nmolasym
            jj = nmolasymptr(1,j)
            loptm(1,j) = lopfi(3*(jj+nsft-1)+1)
            loptm(2,j) = lopfi(3*(jj+nsft-1)+2)
            loptm(3,j) = lopfi(3*(jj+nsft-1)+3)
            if (nmolasymno(j).gt.1) then
              jj = nmolasymptr(2,j)
              loptm(4,j) = lopfi(3*(jj+nsft-1)+1)
              loptm(5,j) = lopfi(3*(jj+nsft-1)+2)
              loptm(6,j) = lopfi(3*(jj+nsft-1)+3)
            else
              loptm(4:6,j) = .false.
            endif
          enddo
        endif
      elseif (lcello.or.lnoflags.or.lbreathe) then
!
!  Atoms
!
        do j = 1,nasymnomol
          loptx(j) = .false.
          lopty(j) = .false.
          loptz(j) = .false.
        enddo
!
!  Molecules
!
        do j = 1,nmolasym
          loptm(1:6,j) = .false.
        enddo
      elseif (lshello) then
!
!  Atoms
!
        do j = 1,nasymnomol
          jj = nasymnomolptr(j)
          if (iatn(jj).lt.maxele) then
            loptx(j) = .false.
            lopty(j) = .false.
            loptz(j) = .false.
          else
            loptx(j) = .true.
            lopty(j) = .true.
            loptz(j) = .true.
          endif
        enddo
!
!  Molecules
!
        do j = 1,nmolasym
          loptm(1:6,j) = .false.
        enddo
      else
!
!  Atoms
!
        do j = 1,nasymnomol
          loptx(j) = .true.
          lopty(j) = .true.
          loptz(j) = .true.
        enddo
!
!  Molecules
!
        if (lnorotate) then
          do j = 1,nmolasym
            loptm(1:3,j) = .true.
          enddo
          do j = 1,nmolasym
            loptm(4:6,j) = .false.
          enddo
        else
          do j = 1,nmolasym
            loptm(1:6,j) = .true.
          enddo
        endif
      endif
!
!  Unfreezing option
!
      if (lufree(ncf)) call unfreeze(loptx,lopty,loptz)
!
!  Fix rigid region atom directions
!
      if (nregions(icfg).gt.1) then
        do j = 1,nasymnomol
          jj = nasymnomolptr(j)
          nrj = nregionno(nsft+jj)
          if (lregionrigid(nrj,ncf)) then
            if (.not.lopfreg(3*(nrj-1)+1,ncf)) then
              loptx(j) = .false.
            endif
            if (.not.lopfreg(3*(nrj-1)+2,ncf)) then
              lopty(j) = .false.
            endif
            if (.not.lopfreg(3*(nrj-1)+3,ncf)) then
              loptz(j) = .false.
            endif
          endif
        enddo
!
!  Fix rigid molecules in rigid regions
!
        if (lrigid) then
          do j = 1,nmolasym
            jj = nmolasymptr(1,j)
            nrj = nregionno(nsft+jj)
            if (lregionrigid(nrj,ncf)) then
              if (.not.lopfreg(3*(nrj-1)+1,ncf)) then
                loptm(1,j) = .false.
              endif
              if (.not.lopfreg(3*(nrj-1)+2,ncf)) then
                loptm(2,j) = .false.
              endif
              if (.not.lopfreg(3*(nrj-1)+3,ncf)) then
                loptm(3,j) = .false.
              endif
!
!  For rigid regions fix orientation of molecules
!
              loptm(4:6,j) = .false.
            endif
          enddo
        endif
      endif
!
      if (lsymopt) then
!
!  Special positions
!
        call special(loptx,lopty,loptz,.false.)
        if (lrigid) then
!
!  Symmetry generation of centres of mass and quaternions
!
          call symupdatemol
!
!  Rigid molecule centre of mass special positions
!
          call specialcom(loptm,.false.)
!
!  Rigid molecule quaternion special positions
!
          call specialqtn(loptm,.false.)
        endif
      endif
!
!  Partial occupancy constraints
!
      if (lspatialok) then
        call poccons(loptx,lopty,loptz,loptr)
      else
        call poccon(loptx,lopty,loptz,loptr)
      endif
      if (lrigid) then
!
!  Rigid molecule partial occupancy
!
! DEBUG - needs adding
      endif
!
!  Rigid region translation
!
      if (nregions(icfg).gt.1) then
        call setregiontrans(loptx,lopty,loptz,loptm)
      endif
      if (ndimen(icfg).eq.3.and.(.not.lflags.or..not.lmc).and.(lfix1atom.and..not.nfixatomtype.eq.5)) then
!
!  For certain space groups there is one arbitary coordinate 
!  in which case one variable in this direction must be removed.
!  This is true if there is no inversion symmetry in this direction.
!
        lodd = .false.
        if (nspg.eq.1.and.ngocfg(icfg).gt.1) then
!
!  General - operators used instead of space group
!       
          sumx = 0.0_dp
          sumy = 0.0_dp
          sumz = 0.0_dp
          do ngo = 1,ngocfg(icfg)
            sumx = sumx + ropcfg(1,1,ngo,icfg)
            sumy = sumy + ropcfg(2,2,ngo,icfg)
            sumz = sumz + ropcfg(3,3,ngo,icfg)       
          enddo
          ndir = 0
          if (abs(sumx).gt.0.0_dp) then
            ndir = ndir + 1
            ndirp(ndir) = 1
          endif
          if (abs(sumy).gt.0.0_dp) then
            ndir = ndir + 1
            ndirp(ndir) = 2
          endif
          if (abs(sumz).gt.0.0_dp) then
            ndir = ndir + 1
            ndirp(ndir) = 3
          endif
          lodd = (ndir.gt.0)
        elseif (nspg.eq.1.and.nccs.eq.1) then
          if (.not.lalleinstein) then
            lodd = .true.
            ndir = 3
            ndirp(1) = 1
            ndirp(2) = 2
            ndirp(3) = 3
          endif
        elseif (nspg.ge.3.and.nspg.le.5) then
!
!  Monoclinic
!
          lodd = .true.
          ndir = 1
          if (abs(alphafull-90.0_dp).gt.1.0d-4) then
            ndirp(1) = 1
          elseif (abs(gammafull-90.0_dp).gt.1.0d-4) then
            ndirp(1) = 3
          else
            ndirp(1) = 2
          endif
        elseif (nspg.ge.6.and.nspg.le.9) then
!
!  Monoclinic
!
          lodd = .true.
          ndir = 2
          if (abs(alphafull-90.0_dp).gt.1.0d-4) then
            ndirp(1) = 2
            ndirp(2) = 3
          elseif (abs(gammafull-90.0_dp).gt.1.0d-4) then
            ndirp(1) = 1
            ndirp(2) = 2
          else
            ndirp(1) = 1
            ndirp(2) = 3
          endif
        elseif (nspg.eq.231.or.nspg.eq.232) then
!
!  Triclinic - C 1 
!
          lodd = .true.
          ndir = 3
          ndirp(1) = 1
          ndirp(2) = 2
          ndirp(3) = 3
        elseif (nspg.ge.25.and.nspg.le.46) then
!
!  Orthorhombic
!
          lodd = .true.
          ndir = 1
          ndirp(1) = 3
        elseif ((nspg.ge.75.and.nspg.le.80).or.(nspg.ge.99.and.nspg.le.110)) then
!
!  Tetragonal
!
          lodd = .true.
          ndir = 1
          ndirp(1) = 3
        elseif (nspg.ge.143.and.nspg.le.146) then
!
!  Trigonal
!
          lodd = .true.
          ndir = 1
          ndirp(1) = 3
        elseif ((nspg.ge.156.and.nspg.le.159).or.(nspg.ge.168.and.nspg.le.173)) then
!
!  Trigonal/Hexagonal
!
          lodd = .true.
          ndir = 1
          ndirp(1) = 3
        elseif (nspg.ge.183.and.nspg.le.186) then
!
!  Hexagonal
!
          lodd = .true.
          ndir = 1
          ndirp(1) = 3
        endif
!
!  Correct flags for lack of inversion symmetry
!
        if (lodd.and.lsymopt) then
          if (nasymnomol.gt.0) then
            allocate(loptxloc(nasymnomol),stat=status)
            if (status/=0) call outofmemory('setcfg','loptxloc')
            allocate(loptyloc(nasymnomol),stat=status)
            if (status/=0) call outofmemory('setcfg','loptyloc')
            allocate(loptzloc(nasymnomol),stat=status)
            if (status/=0) call outofmemory('setcfg','loptzloc')
!
            loptxloc(1:nasymnomol) = loptx(1:nasymnomol)
            loptyloc(1:nasymnomol) = lopty(1:nasymnomol)
            loptzloc(1:nasymnomol) = loptz(1:nasymnomol)
!
            do j = 1,ncon
              if (ncfixtypcfg(ncfst+j).eq.iopt_xf) then
                loptxloc(ncfixindcfg(ncfst+j)) = .true.
              elseif (ncfixtypcfg(ncfst+j).eq.iopt_yf) then
                loptyloc(ncfixindcfg(ncfst+j)) = .true.
              elseif (ncfixtypcfg(ncfst+j).eq.iopt_zf) then
                loptzloc(ncfixindcfg(ncfst+j)) = .true.
              endif
            enddo
!
            lnxi = .true.
            lnyi = .true.
            lnzi = .true.
            do j = 1,nasymnomol
              if (.not.loptxloc(j)) lnxi = .false.
              if (.not.loptyloc(j)) lnyi = .false.
              if (.not.loptzloc(j)) lnzi = .false.
            enddo
            do j = 1,ndir
              if (ndirp(j).eq.1) then
                k = 0
                do while (lnxi.and.k.lt.nasymnomol)
                  k = k + 1
                  if (loptx(k)) then
                    loptx(k) = .false.
                    lnxi = .false.
                  endif
                enddo
              elseif (ndirp(j).eq.2) then
                k = 0
                do while (lnyi.and.k.lt.nasymnomol)
                  k = k + 1
                  if (lopty(k)) then
                    lopty(k) = .false.
                    lnyi = .false.
                  endif
                enddo
              else
                k = 0
                do while (lnzi.and.k.lt.nasymnomol)
                  k = k + 1
                  if (loptz(k)) then
                    loptz(k) = .false.
                    lnzi = .false.
                  endif
                enddo
              endif
            enddo
            deallocate(loptzloc,stat=status)
            if (status/=0) call deallocate_error('setcfg','loptzloc')
            deallocate(loptyloc,stat=status)
            if (status/=0) call deallocate_error('setcfg','loptyloc')
            deallocate(loptxloc,stat=status)
            if (status/=0) call deallocate_error('setcfg','loptxloc')
          endif
          if (nmolasym.gt.0) then
            allocate(loptxloc(nmolasym),stat=status)
            if (status/=0) call outofmemory('setcfg','loptxloc')
            allocate(loptyloc(nmolasym),stat=status)
            if (status/=0) call outofmemory('setcfg','loptyloc')
            allocate(loptzloc(nmolasym),stat=status)
            if (status/=0) call outofmemory('setcfg','loptzloc')
!
            loptxloc(1:nmolasym) = loptm(1,1:nmolasym)
            loptyloc(1:nmolasym) = loptm(2,1:nmolasym)
            loptzloc(1:nmolasym) = loptm(3,1:nmolasym)
!
            do j = 1,ncon
              if (ncfixtypcfg(ncfst+j).eq.iopt_xcom) then
                loptxloc(ncfixindcfg(ncfst+j)) = .true.
              elseif (ncfixtypcfg(ncfst+j).eq.iopt_ycom) then
                loptyloc(ncfixindcfg(ncfst+j)) = .true.
              elseif (ncfixtypcfg(ncfst+j).eq.iopt_zcom) then
                loptzloc(ncfixindcfg(ncfst+j)) = .true.
              endif
            enddo
!
            lnxi = .true.
            lnyi = .true.
            lnzi = .true.
            do j = 1,nmolasym
              if (.not.loptxloc(j)) lnxi = .false.
              if (.not.loptyloc(j)) lnyi = .false.
              if (.not.loptzloc(j)) lnzi = .false.
            enddo
            do j = 1,ndir
              if (ndirp(j).eq.1) then
                k = 0
                do while (lnxi.and.k.lt.nmolasym)
                  k = k + 1
                  if (loptm(1,k)) then
                    loptm(1,k) = .false.
                    lnxi = .false.
                  endif
                enddo
              elseif (ndirp(j).eq.2) then
                k = 0
                do while (lnyi.and.k.lt.nmolasym)
                  k = k + 1
                  if (loptm(2,k)) then
                    loptm(2,k) = .false.
                    lnyi = .false.
                  endif
                enddo
              else
                k = 0
                do while (lnzi.and.k.lt.nmolasym)
                  k = k + 1
                  if (loptm(3,k)) then
                    loptm(3,k) = .false.
                    lnzi = .false.
                  endif
                enddo
              endif
            enddo
            deallocate(loptzloc,stat=status)
            if (status/=0) call deallocate_error('setcfg','loptzloc')
            deallocate(loptyloc,stat=status)
            if (status/=0) call deallocate_error('setcfg','loptyloc')
            deallocate(loptxloc,stat=status)
            if (status/=0) call deallocate_error('setcfg','loptxloc')
          endif
        endif
      endif
!
!  If not symmetry optimisation or all Einstein model set derivatives of one atom to zero
!
!  For a slab, find atom in the middle to fix
!
      if (.not.lsymopt.and..not.lshello.and..not.lflags.and..not.lalleinstein.and..not.lneb &
          .and.(lfix1atom.and..not.nfixatomtype.eq.5)) then
        if (nregions(icfg).eq.1) then
          if (nfixatomtype.gt.0) then
!##########################################
!  Set fixed atom based on input options  #
!##########################################
            if (nfixatomtype.eq.1) then
!
!  First atom
!
              nfixatom = 1
            elseif (nfixatomtype.eq.2) then
!
!  Last atom - excluding shells
!
              if (nasymnomol.gt.0) then
                nfixatom = 1
                do i = 1,nasymnomol
                  ii = nasymnomolptr(i)
                  if (iatn(ii).le.maxele) then
                    nfixatom = i
                  endif
                enddo
              else
!
!  Last molecule
!
                nfixatom = nmolasym
              endif
            elseif (nfixatomtype.eq.3) then
!
!  Centre atom
!
              if (ndim.eq.2) then
!
!  Find atom close to middle of slab to fix
!
                zmid = 0.0_dp
                do i = 1,ncore
                  ic = ncoptr(i)
                  zmid = zmid + zclat(ic)
                enddo
                zmid = zmid/dble(ncore)
                zbest = 1.0d10
                do i = 1,nasymnomol
                  ii = nasymnomolptr(i)
                  if (iatn(ii).le.maxele) then
                    if (abs(zclat(ii)-zmid).lt.zbest) then
                      nfixatom = i
                      zbest = abs(zclat(ii)-zmid)
                    endif
                  endif
                enddo
              elseif (ndim.eq.1) then
!
!  Find atom close to the middle of polymer to fix
!
                ymid = 0.0_dp
                zmid = 0.0_dp
                do i = 1,numat
                  ic = ncoptr(i)
                  ymid = ymid + yclat(ic)
                  zmid = zmid + zclat(ic)
                enddo
                ymid = ymid/dble(ncore)
                zmid = zmid/dble(ncore)
                r2best = 1.0d20
                do i = 1,nasymnomol
                  ii = nasymnomolptr(i)
                  if (iatn(ii).le.maxele) then
                    r2 = (yclat(ii)-ymid)**2 + (zclat(ii)-zmid)**2
                    if (r2.lt.r2best) then
                      nfixatom = i
                      r2best = r2
                    endif
                  endif
                enddo
              else
!
!  Find atom close to the middle of cluster or solid to fix
!
                xmid = 0.0_dp
                ymid = 0.0_dp
                zmid = 0.0_dp
                do i = 1,ncore
                  ic = ncoptr(i)
                  xmid = xmid + xclat(ic)
                  ymid = ymid + yclat(ic)
                  zmid = zmid + zclat(ic)
                enddo
                xmid = xmid/dble(ncore)
                ymid = ymid/dble(ncore)
                zmid = zmid/dble(ncore)
                r2best = 1.0d20
                do i = 1,nasymnomol
                  ii = nasymnomolptr(i)
                  if (iatn(ii).le.maxele) then
                    r2 = (xclat(ii)-xmid)**2 + (yclat(ii)-ymid)**2 + (zclat(ii)-zmid)**2
                    if (r2.lt.r2best) then
                      nfixatom = i
                      r2best = r2
                    endif
                  endif
                enddo
              endif
            endif
            nfixatomcfg(ncf) = nfixatom
          else
!#####################################
!  Set fixed atom based on defaults  #
!#####################################
            if (nprocs.gt.1) then
!
!  Last atom - excluding shells
!
              if (nasymnomol.gt.0) then
                nfixatom = 1
                do i = 1,nasymnomol
                  ii = nasymnomolptr(i)
                  if (iatn(ii).le.maxele) then
                    nfixatom = i
                  endif
                enddo
              else
!
!  Last molecule
!
                nfixatom = nmolasym
              endif
            else
              nfixatom = 1
            endif
            if (ndim.eq.2) then
!
!  Find atom close to middle of slab to fix
!
              zmid = 0.0_dp
              do i = 1,ncore
                ic = ncoptr(i)
                zmid = zmid + zclat(ic)
              enddo
              zmid = zmid/dble(ncore)
              zbest = 1.0d10
              do i = 1,nasymnomol
                ii = nasymnomolptr(i)
                if (iatn(ii).le.maxele) then
                  if (abs(zclat(ii)-zmid).lt.zbest) then
                    nfixatom = i
                    zbest = abs(zclat(ii)-zmid)
                  endif
                endif
              enddo
            elseif (ndim.eq.1) then
!
!  Find atom close to the middle of polymer to fix
!
              ymid = 0.0_dp
              zmid = 0.0_dp
              do i = 1,ncore
                ic = ncoptr(i)
                ymid = ymid + yclat(ic)
                zmid = zmid + zclat(ic)
              enddo
              ymid = ymid/dble(ncore)
              zmid = zmid/dble(ncore)
              r2best = 1.0d20
              do i = 1,nasymnomol
                ii = nasymnomolptr(i)
                if (iatn(ii).le.maxele) then
                  r2 = (yclat(ii)-ymid)**2 + (zclat(ii)-zmid)**2
                  if (r2.lt.r2best) then
                    nfixatom = i
                    r2best = r2
                  endif
                endif
              enddo
            elseif (ndim.eq.0) then
!
!  Find atom close to the middle of cluster to fix
!
              xmid = 0.0_dp
              ymid = 0.0_dp
              zmid = 0.0_dp
              do i = 1,ncore
                ic = ncoptr(i)
                xmid = xmid + xclat(ic)
                ymid = ymid + yclat(ic)
                zmid = zmid + zclat(ic)
              enddo
              xmid = xmid/dble(ncore)
              ymid = ymid/dble(ncore)
              zmid = zmid/dble(ncore)
              r2best = 1.0d20
              do i = 1,nasymnomol
                ii = nasymnomolptr(i)
                if (iatn(ii).le.maxele) then
                  r2 = (xclat(ii)-xmid)**2 + (yclat(ii)-ymid)**2 + (zclat(ii)-zmid)**2
                  if (r2.lt.r2best) then
                    nfixatom = i
                    r2best = r2
                  endif
                endif
              enddo
            endif
            nfixatomcfg(ncf) = nfixatom
          endif
        endif
        if (nasymnomol.gt.0) then
          lfound = .false.
          do i = 1,nasymnomol
            if (.not.loptx(i)) lfound = .true.
          enddo
          if (.not.lfound) loptx(nfixatom) = .false.
          lfound = .false.
          do i = 1,nasymnomol
            if (.not.lopty(i)) lfound = .true.
          enddo
          if (.not.lfound) lopty(nfixatom) = .false.
          lfound = .false.
          do i = 1,nasymnomol
            if (.not.loptz(i)) lfound = .true.
          enddo
!
!  If there is a plane potential and this is 2-D then there is no need to fix
!  a coordinate in the z-direction since the plane prevents translation
!
          if ((ndim.ne.0.and.ndim.ne.2).or.nplanepot.eq.0) then
            if (.not.lfound) loptz(nfixatom) = .false.
          endif
        else
          lfound = .false.
          do i = 1,nmolasym
            if (.not.loptm(1,i)) lfound = .true.
          enddo
          if (.not.lfound) loptm(1,nfixatom) = .false.
          lfound = .false.
          do i = 1,nmolasym
            if (.not.loptm(2,i)) lfound = .true.
          enddo
          if (.not.lfound) loptm(2,nfixatom) = .false.
          lfound = .false.
          do i = 1,nmolasym
            if (.not.loptm(3,i)) lfound = .true.
          enddo
!
!  If there is a plane potential and this is 2-D then there is no need to fix
!  a coordinate in the z-direction since the plane prevents translation
!
          if ((ndim.ne.0.and.ndim.ne.2).or.nplanepot.eq.0) then
            if (.not.lfound) loptm(3,nfixatom) = .false.
          endif
        endif
      endif
!
!  For MD no need to fix any atoms
!
      if (lmd.and..not.lflags.and.numat.gt.0) then
        if (nasymnomol.gt.0) then
          loptx(nfixatom) = .true.
          lopty(nfixatom) = .true.
          loptz(nfixatom) = .true.
        else
          loptm(1:3,nfixatom) = .true.
        endif
      endif
!
!  Check constraints to remove non-allowed variables
!
      if (ncon.gt.0) then
        do i = 1,ncon
          if (ncfixtypcfg(ncfst+i).eq.iopt_xf) then
            loptx(ncfixindcfg(ncfst+i)) = .false.
          elseif (ncfixtypcfg(ncfst+i).eq.iopt_yf) then
            lopty(ncfixindcfg(ncfst+i)) = .false.
          elseif (ncfixtypcfg(ncfst+i).eq.iopt_zf) then
            loptz(ncfixindcfg(ncfst+i)) = .false.
          elseif (ncfixtypcfg(ncfst+i).eq.iopt_xcom) then
            loptm(1,ncfixindcfg(ncfst+i)) = .false.
          elseif (ncfixtypcfg(ncfst+i).eq.iopt_ycom) then
            loptm(2,ncfixindcfg(ncfst+i)) = .false.
          elseif (ncfixtypcfg(ncfst+i).eq.iopt_zcom) then
            loptm(3,ncfixindcfg(ncfst+i)) = .false.
          elseif (ncfixtypcfg(ncfst+i).eq.iopt_xqtn) then
            loptm(4,ncfixindcfg(ncfst+i)) = .false.
          elseif (ncfixtypcfg(ncfst+i).eq.iopt_yqtn) then
            loptm(5,ncfixindcfg(ncfst+i)) = .false.
          elseif (ncfixtypcfg(ncfst+i).eq.iopt_zqtn) then
            loptm(6,ncfixindcfg(ncfst+i)) = .false.
          endif
        enddo
      endif
    else
      loptm(1:6,1:nmolasym) = .false.
      loptx(1:nasymnomol) = .false.
      lopty(1:nasymnomol) = .false.
      loptz(1:nasymnomol) = .false.
    endif
!**********
!  Radii  *
!**********
    if (.not.lbulknoopt.and.index(keyword,'nobr').eq.0) then
      do j = 1,nasymnomol
        jj = nasymnomolptr(j)
        if (lbsmat(jj+nsft)) then
          nrj = nregionno(nsft+jj)
          if (lregionrigid(nrj,ncf)) then
            loptr(j) = .false.
          else
            loptr(j) = .true.
          endif
        else
          loptr(j) = .false.
        endif
      enddo
      if (ncon.gt.0) then
        do i = 1,ncon
          if (ncfixtypcfg(ncfst+i).eq.iopt_radius) then
            loptr(ncfixindcfg(ncfst+i)) = .false.
          endif
        enddo
      endif
    else
      loptr(1:nasymnomol) = .false.
    endif
!**************************************
!  Initialise variable pointer array  *
!**************************************
    n1var(icfg) = nvarl + 1
!
!  Count number of variables and check memory
!
    nvarlold = nvarl
    do j = 1,nstrains
      if (loptc(j)) then
        nvarl = nvarl + 1
      endif
    enddo
    do j = 1,nasymnomol
      if (loptx(j)) then
        nvarl = nvarl + 1
      endif
      if (lopty(j)) then
        nvarl = nvarl + 1
      endif
      if (loptz(j)) then
        nvarl = nvarl + 1
      endif
      if (loptr(j)) then
        nvarl = nvarl + 1
      endif
    enddo
    do j = 1,nmolasym
      do k = 1,6
        if (loptm(k,j)) then
          nvarl = nvarl + 1
        endif
      enddo
    enddo
    if (nvarl.ge.maxvar) then
      maxvar = nvarl + 10
      call changemaxvar
    endif
    nvarl = nvarlold
    if (loldvarorder) then
!########################################
!  Build list of variables - old order  #
!########################################
!
!  Strain or cell
!
      ncellmincfg(icfg) = 1
      do j = 1,nstrains
        if (loptc(j)) then
          nvarl = nvarl + 1
          ioptindexcfg(nvarl) = j
          if (loptcellpar) then
            iopttypecfg(nvarl) = iopt_cell
          else
            iopttypecfg(nvarl) = iopt_strain
          endif
        endif
      enddo
      ncellmaxcfg(icfg) = nvarl - nvarlold
!
!  Internals
!
      ninternalmincfg(icfg) = ncellmaxcfg(icfg) + 1
      if (lrigid) then
        do j = 1,nasymnomol
          if (loptx(j)) then
            nvarl = nvarl + 1
            ioptindexcfg(nvarl) = j
            iopttypecfg(nvarl) = iopt_xf
          endif
          if (lopty(j)) then
            nvarl = nvarl + 1
            ioptindexcfg(nvarl) = j
            iopttypecfg(nvarl) = iopt_yf
          endif
          if (loptz(j)) then
            nvarl = nvarl + 1
            ioptindexcfg(nvarl) = j
            iopttypecfg(nvarl) = iopt_zf
          endif
        enddo
!
!  Radius
!
        do j = 1,nasymnomol
          if (loptr(j)) then
            nvarl = nvarl + 1
            ioptindexcfg(nvarl) = j
            iopttypecfg(nvarl) = iopt_radius
          endif
        enddo
!
!  Rigid molecules
!
        do j = 1,nmolasym
          if (loptm(1,j)) then
            nvarl = nvarl + 1
            ioptindexcfg(nvarl) = j
            iopttypecfg(nvarl) = iopt_xcom
          endif
          if (loptm(2,j)) then
            nvarl = nvarl + 1
            ioptindexcfg(nvarl) = j
            iopttypecfg(nvarl) = iopt_ycom
          endif
          if (loptm(3,j)) then
            nvarl = nvarl + 1
            ioptindexcfg(nvarl) = j
            iopttypecfg(nvarl) = iopt_zcom
          endif
          if (loptm(4,j)) then
            nvarl = nvarl + 1
            ioptindexcfg(nvarl) = j
            iopttypecfg(nvarl) = iopt_xqtn
          endif
          if (loptm(5,j)) then
            nvarl = nvarl + 1
            ioptindexcfg(nvarl) = j
            iopttypecfg(nvarl) = iopt_yqtn
          endif
          if (loptm(6,j)) then
            nvarl = nvarl + 1
            ioptindexcfg(nvarl) = j
            iopttypecfg(nvarl) = iopt_zqtn
          endif
        enddo
      else
        do j = 1,nasym
          if (loptx(j)) then
            nvarl = nvarl + 1
            ioptindexcfg(nvarl) = j
            iopttypecfg(nvarl) = iopt_xf
          endif
          if (lopty(j)) then
            nvarl = nvarl + 1
            ioptindexcfg(nvarl) = j
            iopttypecfg(nvarl) = iopt_yf
          endif
          if (loptz(j)) then
            nvarl = nvarl + 1
            ioptindexcfg(nvarl) = j
            iopttypecfg(nvarl) = iopt_zf
          endif
        enddo
!
!  Radius
!
        do j = 1,nasym
          if (loptr(j)) then
            nvarl = nvarl + 1
            ioptindexcfg(nvarl) = j
            iopttypecfg(nvarl) = iopt_radius
          endif
        enddo
      endif
      ninternalmaxcfg(icfg) = nvarl - nvarlold
    else
!########################################
!  Build list of variables - new order  #
!########################################
!
!  Internals
!
      ninternalmincfg(icfg) = + 1
      if (lrigid) then
        do j = 1,nasymnomol
          if (loptx(j)) then
            nvarl = nvarl + 1
            ioptindexcfg(nvarl) = j
            iopttypecfg(nvarl) = iopt_xf
          endif
          if (lopty(j)) then
            nvarl = nvarl + 1
            ioptindexcfg(nvarl) = j
            iopttypecfg(nvarl) = iopt_yf
          endif
          if (loptz(j)) then
            nvarl = nvarl + 1
            ioptindexcfg(nvarl) = j
            iopttypecfg(nvarl) = iopt_zf
          endif
        enddo
!
!  Radius
!
        do j = 1,nasymnomol
          if (loptr(j)) then
            nvarl = nvarl + 1
            ioptindexcfg(nvarl) = j
            iopttypecfg(nvarl) = iopt_radius
          endif
        enddo
!
!  Rigid molecules - translations
!
        do j = 1,nmolasym
          if (loptm(1,j)) then
            nvarl = nvarl + 1
            ioptindexcfg(nvarl) = j
            iopttypecfg(nvarl) = iopt_xcom
          endif
          if (loptm(2,j)) then
            nvarl = nvarl + 1
            ioptindexcfg(nvarl) = j
            iopttypecfg(nvarl) = iopt_ycom
          endif
          if (loptm(3,j)) then
            nvarl = nvarl + 1
            ioptindexcfg(nvarl) = j
            iopttypecfg(nvarl) = iopt_zcom
          endif
        enddo
!
!  Rigid molecules - quaternions
!
        do j = 1,nmolasym
          if (loptm(4,j)) then
            nvarl = nvarl + 1
            ioptindexcfg(nvarl) = j
            iopttypecfg(nvarl) = iopt_xqtn
          endif
          if (loptm(5,j)) then
            nvarl = nvarl + 1
            ioptindexcfg(nvarl) = j
            iopttypecfg(nvarl) = iopt_yqtn
          endif
          if (loptm(6,j)) then
            nvarl = nvarl + 1
            ioptindexcfg(nvarl) = j
            iopttypecfg(nvarl) = iopt_zqtn
          endif
        enddo
      else
        do j = 1,nasym
          if (loptx(j)) then
            nvarl = nvarl + 1
            ioptindexcfg(nvarl) = j
            iopttypecfg(nvarl) = iopt_xf
          endif
          if (lopty(j)) then
            nvarl = nvarl + 1
            ioptindexcfg(nvarl) = j
            iopttypecfg(nvarl) = iopt_yf
          endif
          if (loptz(j)) then
            nvarl = nvarl + 1
            ioptindexcfg(nvarl) = j
            iopttypecfg(nvarl) = iopt_zf
          endif
        enddo
!
!  Radius
!
        do j = 1,nasym
          if (loptr(j)) then
            nvarl = nvarl + 1
            ioptindexcfg(nvarl) = j
            iopttypecfg(nvarl) = iopt_radius
          endif
        enddo
      endif
      ninternalmaxcfg(icfg) = nvarl - nvarlold
!
!  Strain
!
      ncellmincfg(icfg) = ninternalmaxcfg(icfg) + 1
      do j = 1,nstrains
        if (loptc(j)) then
          nvarl = nvarl + 1
          ioptindexcfg(nvarl) = j
          if (loptcellpar) then
            iopttypecfg(nvarl) = iopt_cell
          else
            iopttypecfg(nvarl) = iopt_strain
          endif
        endif
      enddo
      ncellmaxcfg(icfg) = nvarl - nvarlold
    endif
    nvarcfg(icfg) = nvarl - n1var(icfg) + 1
!
!  If shrinking factors have been set then set flag
!
    lshrink = ((nks(1,icfg)*nks(2,icfg)*nks(3,icfg)).gt.0)
    if (ioproc) then
!*******************************
!  Configuration based output  *
!*******************************
      write(ioout,'(/,''********************************************************************************'')')
      if (names(icfg)(1:1).eq.' ') then
        write(ioout,'(''*  Input for Configuration = '',i3,47x,''*'')') icfg
      else
        write(ioout,'(''*  Input for Configuration = '',i3,'' : '',a44,''*'')') icfg,names(icfg)(1:44)
      endif
      write(ioout,'(''********************************************************************************'')')
      call formula(ioout)
!
!  Check that total occupancy isn't zero and output warning
!
      occtot = 0.0_dp
      do i = 1,nasym
        occtot = occtot + occua(i)
      enddo
      if (occtot.lt.1.0d-12.and.nasym.gt.0) then
        nwarn = nwarn + 1
        call outwarning('Total occupancy of all atoms is zero and so energy will be zero',0_i4)
      endif
!
      write(ioout,'(/,''  Number of irreducible atoms/shells = '',i7,/)') nasym
      write(ioout,'(/,''  Total number atoms/shells = '',i7,/)') numat
      write(ioout,'(''  Dimensionality = '',i1,15x,'':'',2x,7a)') ndimen(icfg),systype(ndimen(icfg)+1)
      write(ioout,'(/)')
    endif
    sum = 0.0_dp
    do j = 1,numat
      sum = sum + qf(j)*occuf(j)
    enddo
    totalchargecfg(icfg) = sum
    if (ioproc) then
      if (ndimen(icfg).eq.0) then
        write(ioout,'(''  Charge on cluster = '',f10.6,/)') sum
      elseif (abs(sum).ge.1.0d-4) then
        if (ndimen(icfg).eq.3) then
          write(ioout,'(''  Charge on solid   = '',f10.6,'' =>neutralising background added'',/)') sum
        elseif (ndimen(icfg).eq.2) then
          write(ioout,'(''  Charge on slab    = '',f10.6,'' =>energy corrected for self-interaction'',/)') sum
        else
          write(ioout,'(''  Charge on polymer = '',f10.6,'' =>energy corrected for self-interaction'',/)') sum
          nwarn = nwarn + 1
          call outwarning('Charged polymer calculations are not reliable',0_i4)
        endif
      endif
      if (lcosmo) then
        if (lcosmic) then
          write(ioout,'(''  Solvated with COSMIC model : '',/)')
        else
          write(ioout,'(''  Solvated with COSMO model : '',/)')
        endif
        write(ioout,'(''  Dielectric constant        = '',f12.6)') cosmoepsilon(ncf)
        write(ioout,'(''  Solvent radius             = '',f12.6,'' Angstroms'')') cosmorsolv(ncf)
        write(ioout,'(''  Delta solvent radius       = '',f12.6,'' Angstroms'')') cosmodrsolv(ncf)
        write(ioout,'(''  Cutoff : point to segment  = '',f12.6,'' Angstroms'')') cosmormax
        write(ioout,'(''  Smooth : point to segment  = '',f12.6,'' Angstroms'')') cosmormaxs
        write(ioout,'(''  Smooth : point inclusion   = '',f12.6,'' Angstroms'')') cosmorange
        write(ioout,'(''  Coulomb: Wolf sum eta      = '',f12.6,'' Angstroms-1'')') etawc
        write(ioout,'(''  Coulomb: Wolf sum cutoff   = '',f12.6,'' Angstroms'')') cutwc
        write(ioout,'(''  No. of points per sphere   = '',i12)') nppa
        write(ioout,'(''  No. of segments per sphere = '',i12)') nspa
        if (isasatomoption.eq.1) then
          write(ioout,'(''  For generation of SAS use  = cores only'')')
        elseif (isasatomoption.eq.2) then
          write(ioout,'(''  For generation of SAS use  = both cores and shells'')')
        endif
        if (ldodeca) then
          write(ioout,'(''  Shape of atomic mesh       = '',''Dodecahedron'',/)')
        else
          write(ioout,'(''  Shape of atomic mesh       = '',''Octahedron'',/)')
        endif
      endif
      if (QMMMmode(icfg).eq.1) then
        write(ioout,'(''  QM/MM rules for mechanical embedding to be applied to this configuration '')')
      elseif (QMMMmode(icfg).eq.2) then
        write(ioout,'(''  QM/MM rules for electrical embedding to be applied to this configuration '')')
      endif
      if (lsymopt) call symout
!
!  If symmetry has been found then output information for potential space group
!
      if (lfindsym) then
        call findspacegroup(ngroup)
        write(ioout,'(''  Space group found : Number = '',i3,'' : Symbol = '',a16,/)') ngroup,gronam(ngroup)
      endif
!
      if (lsuper) then
        ix = nsuper(1,icfg)
        iy = nsuper(2,icfg)
        iz = nsuper(3,icfg)
        if (ndimen(icfg).eq.3) then
          write(ioout,'(/,''  Supercell dimensions :  x = '',i3,''  y = '',i3,''  z = '',i3)') ix,iy,iz
        elseif (ndimen(icfg).eq.2) then
          write(ioout,'(/,''  Supercell dimensions :  x = '',i3,''  y = '',i3)') ix,iy
        elseif (ndimen(icfg).eq.1) then
          write(ioout,'(/,''  Supercell dimensions :  x = '',i3)') ix
        endif
      elseif (lsuperghost) then
        ix = nsuperghost(1,icfg)
        iy = nsuperghost(2,icfg)
        iz = nsuperghost(3,icfg)
        if (ndimen(icfg).eq.3) then
          write(ioout,'(/,''  Ghost supercell dimensions :  x = '',i3,''  y = '',i3,''  z = '',i3)') ix,iy,iz
        elseif (ndimen(icfg).eq.2) then
          write(ioout,'(/,''  Ghost supercell dimensions :  x = '',i3,''  y = '',i3)') ix,iy
        elseif (ndimen(icfg).eq.1) then
          write(ioout,'(/,''  Ghost supercell dimensions :  x = '',i3)') ix
        endif
      endif
!
!  Force constant supercell
!
      if (lfcsupercell) then
        ix = nd2cellcfg(1,icfg)
        iy = nd2cellcfg(2,icfg)
        iz = nd2cellcfg(3,icfg)
        if (ndimen(icfg).eq.3) then
          write(ioout,'(/,''  Force constant supercell dimensions :  x = '',i3,''  y = '',i3,''  z = '',i3)') ix,iy,iz
        elseif (ndimen(icfg).eq.2) then
          write(ioout,'(/,''  Force constant supercell dimensions :  x = '',i3,''  y = '',i3)') ix,iy
        elseif (ndimen(icfg).eq.1) then
          write(ioout,'(/,''  Force constant supercell dimensions :  x = '',i3)') ix
        endif
      endif
!
      if (ndimen(icfg).eq.3) then
        if (.not.lterseincell) then
          if (lstraincell) then
            write(ioout,'(/,''  Reference Cartesian lattice vectors (Angstroms) :'',/)')
            do i = 1,3
              write(ioout,'(4x,3f12.6)')(rvcfg(j,i,icfg),j=1,3)
            enddo
            write(ioout,'(/,''  Initial strains :'',/)')
            write(ioout,'(4x,6f12.6)')(straincfg(j,icfg),j=1,6)
          endif
          write(ioout,'(/,''  Cartesian lattice vectors (Angstroms) :'',/)')
          do i = 1,3
            write(ioout,'(4x,3f12.6)')(rv(j,i),j=1,3)
          enddo
        endif
        if (ncbl.gt.1.and.ifhr(icfg).eq.0) then
          aprim = a
          bprim = b
          cprim = c
          alpprim = alpha
          betprim = beta
          gamprim = gamma
          do i = 1,3
            rvt(1,i) = rv(1,i)
            rvt(2,i) = rv(2,i)
            rvt(3,i) = rv(3,i)
          enddo
          call uncentre(rvt)
          call uncell3D(rvt,a,b,c,alpha,beta,gamma)
          if (.not.lterseincell) then
            write(ioout,'(/,''  Primitive cell parameters :'',10x,''  Full cell parameters :'',/)')
            write(ioout,'(''  a = '',f8.4,''    alpha = '',f8.4,4x,''   a = '',f8.4,''    alpha = '',f8.4)') &
              aprim,alpprim,a,alpha
            write(ioout,'(''  b = '',f8.4,''    beta  = '',f8.4,4x,''   b = '',f8.4,''    beta  = '',f8.4)') &
              bprim,betprim,b,beta
            write(ioout,'(''  c = '',f8.4,''    gamma = '',f8.4,4x,''   c = '',f8.4,''    gamma = '',f8.4)') &
              cprim,gamprim,c,gamma
          endif
          volp = volume(rvt)
          vol = volume(rv)
          if (.not.lterseincell) then
            write(ioout,'(/,''  Initial volumes (Angstroms**3):'')')
            write(ioout,'(/,''  Primitive cell = '',f18.6,2x,''  Full cell = '',f18.6)')vol,volp
          endif
        else
          vol = volume(rv)
          if (.not.lterseincell) then
            write(ioout,'(/,''  Cell parameters (Angstroms/Degrees):'',/)')
            write(ioout,'(''  a = '',f12.4,''    alpha = '',f8.4)') a,alpha
            write(ioout,'(''  b = '',f12.4,''    beta  = '',f8.4)') b,beta
            write(ioout,'(''  c = '',f12.4,''    gamma = '',f8.4)') c,gamma
            write(ioout,'(/,''  Initial cell volume = '',f18.6,'' Angs**3'')') vol
          endif
        endif
!
!  Check that cell parameters are consistent with space group
!
        if (lsym) then
          call cellcheck(nspg,a,b,c,alpha,beta,gamma)
        endif
        if (lfieldcfg(icfg)) then
          write(ioout,'(/,''  Electric field applied    = '',f13.6,'' eV/Ang.e '')') fieldcfg(icfg)
          write(ioout,'(''  Electric field direction  = '',f13.6,'' a '')') fielddirectioncfg(1,icfg)
          write(ioout,'(''                              '',f13.6,'' b '')') fielddirectioncfg(2,icfg)
          write(ioout,'(''                              '',f13.6,'' c '')') fielddirectioncfg(3,icfg)
        endif
        if (ntdfieldcfg(icfg).gt.0) then
          write(ioout,'(/,''  Number of TD-Electric fields = '',i13)') ntdfieldcfg(icfg)
          do j = 1,ntdfieldcfg(icfg)
            write(ioout,'(/,''  TD-Electric field applied    = '',f13.6,'' eV/Ang.e '')') td_fieldcfg(1,j,icfg)
            write(ioout,'(''  Time constant (rate)  B      = '',f13.6,'' ps^-1 '')') td_fieldcfg(2,j,icfg)
            write(ioout,'(''  Time constant (phase) C      = '',f13.6)') td_fieldcfg(3,j,icfg)
            write(ioout,'(/,''  TD-Electric field direction  = '',f13.6,'' a '')') td_fielddirectioncfg(1,j,icfg)
            write(ioout,'(''                                 '',f13.6,'' b '')') td_fielddirectioncfg(2,j,icfg)
            write(ioout,'(''                                 '',f13.6,'' c '')') td_fielddirectioncfg(3,j,icfg)
          enddo
        endif
      elseif (ndimen(icfg).eq.2) then
        if (.not.lterseincell) then
          if (lstraincell) then
            write(ioout,'(/,''  Reference surface Cartesian vectors (Angstroms) :'',/)')
            do i = 1,2
              write(ioout,'(4x,3f12.6)')(rvcfg(j,i,icfg),j=1,3)
            enddo
            write(ioout,'(/,''  Initial strains :'',/)')
            write(ioout,'(4x,6f12.6)')(straincfg(j,icfg),j=1,3)
          endif
          write(ioout,'(/,''  Surface Cartesian vectors (Angstroms) :'',/)')
          do i = 1,2
            write(ioout,'(4x,3f12.6)')(rv(j,i),j=1,3)
          enddo
        endif
        ara = area(rv)
        if (.not.lterseincell) then
          write(ioout,'(/,''  Surface cell parameters (Angstroms/Degrees):'',/)')
          write(ioout,'(''  a = '',f12.4,''    alpha = '',f8.4)') a,alpha
          write(ioout,'(''  b = '',f12.4)') b
          write(ioout,'(/,''  Initial surface area   = '',f13.6,'' Angs**2'')') ara
        endif
        call getdipole2D(dipolez)
        write(ioout,'(/,''  Initial surface dipole = '',f13.6,'' e.Angs'')') dipolez
        if (lfieldcfg(icfg)) then
          write(ioout,'(/,''  Electric field applied   = '',f13.6,'' eV/Ang.e'')') fieldcfg(icfg)
          write(ioout,'(''  Electric field direction = '',f13.6,'' a '')') fielddirectioncfg(1,icfg)
          write(ioout,'(''                             '',f13.6,'' b '')') fielddirectioncfg(2,icfg)
          write(ioout,'(''                             '',f13.6,'' z '')') fielddirectioncfg(3,icfg)
        endif
      elseif (ndimen(icfg).eq.1) then
        if (.not.lterseincell) then
          if (lstraincell) then
            write(ioout,'(/,''  Reference polymer Cartesian vector (Angstroms) :'',/)')
            write(ioout,'(4x,3f12.6)') (rvcfg(j,1,icfg),j=1,3)
            write(ioout,'(/,''  Initial strain :'',/)')
            write(ioout,'(4x,f12.6)') straincfg(1,icfg)
          endif
          write(ioout,'(/,''  Polymer Cartesian vector (Angstroms) :'',/)')
          write(ioout,'(4x,3f12.6)') (rv(j,1),j=1,3)
        endif
        if (.not.lterseincell) then
          write(ioout,'(/,''  Polymer cell parameter (Angstrom):'',/)')
          write(ioout,'(''  a = '',f12.4)') a
        endif
        call getdipole1D(dipoley,dipolez)
        write(ioout,'(/,''  Initial polymer dipoles : y = '',f13.6,'' e.Angs'')') dipoley
        write(ioout,'(''                            z = '',f13.6,'' e.Angs'')') dipolez
        if (lfieldcfg(icfg)) then
          write(ioout,'(/,''  Electric field applied   = '',f13.6,'' eV/Ang.e'')') fieldcfg(icfg)
          write(ioout,'(''  Electric field direction = '',f13.6,'' a '')') fielddirectioncfg(1,icfg)
          write(ioout,'(''                             '',f13.6,'' y '')') fielddirectioncfg(2,icfg)
          write(ioout,'(''                             '',f13.6,'' z '')') fielddirectioncfg(3,icfg)
        endif
      elseif (ndimen(icfg).eq.0) then
        call getdipole0D(dipolex,dipoley,dipolez)
        write(ioout,'(/,''  Initial cluster dipoles : x = '',f13.6,'' e.Angs'')') dipolex
        write(ioout,'(''                            y = '',f13.6,'' e.Angs'')') dipoley
        write(ioout,'(''                            z = '',f13.6,'' e.Angs'')') dipolez
        if (lfieldcfg(icfg)) then
          fieldnorm = fielddirectioncfg(1,icfg)**2 + fielddirectioncfg(2,icfg)**2 + fielddirectioncfg(3,icfg)**2
          fieldnorm = fieldcfg(icfg)/sqrt(fieldnorm)
          write(ioout,'(/,''  Electric field applied  : x = '',f13.6,'' eV/Ang.e'')') fielddirectioncfg(1,icfg)*fieldnorm
          write(ioout,'(  ''                          : y = '',f13.6,'' eV/Ang.e'')') fielddirectioncfg(2,icfg)*fieldnorm
          write(ioout,'(  ''                          : z = '',f13.6,'' eV/Ang.e'')') fielddirectioncfg(3,icfg)*fieldnorm
        endif
        if (lradialcfg(icfg)) then
          write(ioout,'(/,''  Radial force applied    : K = '',f13.6,'' eV/Ang**2'')') radialKcfg(icfg)
          write(ioout,'(  ''                          : x = '',f13.6,'' Ang'')') radialXYZcfg(1,icfg)
          write(ioout,'(  ''                          : y = '',f13.6,'' Ang'')') radialXYZcfg(2,icfg)
          write(ioout,'(  ''                          : z = '',f13.6,'' Ang'')') radialXYZcfg(3,icfg)
        endif
      endif
      if (lshrink) then
        if (ndimen(icfg).eq.3) then
          write(ioout,'(/,''  Shrinking factors = '',3(4x,i2))') nks(1,icfg),nks(2,icfg),nks(3,icfg)
        elseif (ndimen(icfg).eq.2) then
          write(ioout,'(/,''  Shrinking factors = '',2(4x,i2))') nks(1,icfg),nks(2,icfg)
        elseif (ndimen(icfg).eq.1) then
          write(ioout,'(/,''  Shrinking factor = '',4x,i2)') nks(1,icfg) 
        endif
        if (lksorigin(icfg)) then
          write(ioout,'(/,''  Shrinking factor grid to be centred at the origin '')')
        endif
        if (lksorigin(icfg)) then
          write(ioout,'(/,''  Shrinking factor grid to be centred at the origin '')')
        endif
      endif
      write(ioout,'(/,''  Temperature of configuration = '',g10.3,'' K '')') tempcfg(icfg)
      if (ndimen(icfg).eq.3) then
        write(ioout,'(/,''  Pressure of configuration = '',f13.3,'' GPa '')') presscfg(icfg)
        if (lanisotropicpresscfg(icfg)) then
          if (ndimen(icfg).eq.3) then
            write(ioout,'(/,''  Anisotropic Pressure Tensor (GPa) : xx = '',f13.3)') anisotropicpresscfg(1,icfg)
            write(ioout,'(''                                    : yy = '',f13.3)') anisotropicpresscfg(2,icfg)
            write(ioout,'(''                                    : zz = '',f13.3)') anisotropicpresscfg(3,icfg)
            write(ioout,'(''                                    : yz = '',f13.3)') anisotropicpresscfg(4,icfg)
            write(ioout,'(''                                    : xz = '',f13.3)') anisotropicpresscfg(5,icfg)
            write(ioout,'(''                                    : xy = '',f13.3)') anisotropicpresscfg(6,icfg)
          elseif (ndimen(icfg).eq.2) then
            write(ioout,'(/,''  Anisotropic Pressure Tensor (GPa) : xx = '',f13.3)') anisotropicpresscfg(1,icfg)
            write(ioout,'(''                                    : yy = '',f13.3)') anisotropicpresscfg(2,icfg)
            write(ioout,'(''                                    : xy = '',f13.3)') anisotropicpresscfg(3,icfg)
          elseif (ndimen(icfg).eq.1) then
            write(ioout,'(/,''  Anisotropic Pressure Tensor (GPa) : xx = '',f13.3)') anisotropicpresscfg(1,icfg)
          endif
        endif
      elseif (ndimen(icfg).eq.0) then
        write(ioout,'(/,''  Symmetry number = '',i6)') symnocfg(icfg)
      endif
      if (dielectriccfg(icfg).gt.1.0_dp) then
        write(ioout,'(/,''  Dielectric constant for configuration = '',g10.5)') dielectriccfg(icfg)
      endif
      if (.not.lterseincoords) then
        if (ndimen(icfg).eq.3) then
          if (.not.lpredict) write(ioout,'(/,''  Fractional coordinates of asymmetric unit :'',/)')
        elseif (ndimen(icfg).eq.2) then
          write(ioout,'(/,''  Mixed fractional/Cartesian coordinates of surface :'',/)')
        elseif (ndimen(icfg).eq.1) then
          write(ioout,'(/,''  Mixed fractional/Cartesian coordinates of polymer :'',/)')
        else
          write(ioout,'(/,''  Cartesian coordinates of cluster :'',/)')
        endif
      endif
      if (.not.lpredict) then
        if (.not.lterseincoords) then
          write(ioout,'(''--------------------------------------------------------------------------------'')')
          write(ioout,'(''   No.  Atomic       x           y          z         Charge      Occupancy'')')
          if (ndimen(icfg).eq.3) then
            write(ioout,'(''        Label      (Frac)      (Frac)     (Frac)        (e)         (Frac)  '')')
          elseif (ndimen(icfg).eq.2) then
            write(ioout,'(''        Label      (Frac)      (Frac)     (Angs)        (e)         (Frac)  '')')
          elseif (ndimen(icfg).eq.1) then
            write(ioout,'(''        Label      (Frac)      (Angs)     (Angs)        (e)         (Frac)  '')')
          else
            write(ioout,'(''        Label      (Angs)      (Angs)     (Angs)        (e)         (Frac)  '')')
          endif
          write(ioout,'(''--------------------------------------------------------------------------------'')')
        endif
        do nr = 1,nregions(icfg)
          if (nregions(icfg).gt.1) then
            if (lregionrigid(nr,ncf)) then
              fixstring = ' '
              lfixeddirection = .false.
              if (.not.lopfreg(3*(nr-1)+1,icfg)) then
                fixstring(1:1) = 'x'
                lfixeddirection = .true.
              endif
              if (.not.lopfreg(3*(nr-1)+2,icfg)) then
                fixstring(2:2) = 'y'
                lfixeddirection = .true.
              endif
              if (.not.lopfreg(3*(nr-1)+3,icfg)) then
                fixstring(3:3) = 'z'
                lfixeddirection = .true.
              endif
              if (.not.lterseincoords) then
                if (lfixeddirection) then
                  if (nregiontype(nr,icfg).eq.1) then
                    write(ioout,'(''  Region '',i1,'' :  QM : Rigid translation fixed in '',a3)') nr,fixstring
                  elseif (nregiontype(nr,icfg).eq.2) then
                    write(ioout,'(''  Region '',i1,'' :  MM : Rigid translation fixed in '',a3)') nr,fixstring
                  else
                    write(ioout,'(''  Region '',i1,'' :  Rigid translation fixed in '',a3)') nr,fixstring
                  endif
                else
                  if (nregiontype(nr,icfg).eq.1) then
                    write(ioout,'(''  Region '',i1,'' : QM : Rigid translation '')') nr
                  elseif (nregiontype(nr,icfg).eq.2) then
                    write(ioout,'(''  Region '',i1,'' : MM : Rigid translation '')') nr
                  else
                    write(ioout,'(''  Region '',i1,'' : Rigid translation '')') nr
                  endif
                endif
              endif
            else
              if (.not.lterseincoords) then
                if (nregiontype(nr,icfg).eq.1) then
                  write(ioout,'(''  Region '',i1,'' : QM '')') nr
                elseif (nregiontype(nr,icfg).eq.2) then
                  write(ioout,'(''  Region '',i1,'' : MM '')') nr
                else
                  write(ioout,'(''  Region '',i1,'' : '')') nr
                endif
              endif
            endif
            write(ioout,'(''--------------------------------------------------------------------------------'')')
          endif
          do i = 1,nasym
            if (nregionno(nsft+i).eq.nr) then
              inat = iatn(i)
              itype = natype(i)
!
!  Hide shells?
!
              if (inat.le.maxele.or..not.lhideshells) then
                call label(inat,itype,lab)
                if (lbsmat(i+nsft)) then
                  cstype = 'bc'
                  if (inat.gt.maxele) cstype = 'bs'
                else
                  cstype = 'c '
                  if (inat.gt.maxele) cstype = 's '
                endif
                if (nasymnomolrptr(i).ne.0) then
                  if (loptx(nasymnomolrptr(i))) then
                    ocha(1) = '*'
                  else
                    ocha(1) = ' '
                  endif
                  if (lopty(nasymnomolrptr(i))) then
                    ocha(2) = '*'
                  else
                    ocha(2) = ' '
                  endif
                  if (loptz(nasymnomolrptr(i))) then
                    ocha(3) = '*'
                  else
                    ocha(3) = ' '
                  endif
                  if (lfix(nasymnomolrptr(i))) then
                    fixed = 'fix'
                  else
                    fixed = '   '
                  endif
                else
                  ocha(1:3) = ' '
                  fixed = '   '
                endif
                if (.not.lterseincoords) then
                  if (ndimen(icfg).eq.3) then
                    if (lrhombo) then
                      nri = nrela2f(i)
                      write(ioout,'(i7,1x,a5,1x,a2,2x,f9.6,2(1x,a1,1x,f9.6),1x,a1,1x,f9.5,f12.6,1x,a3)') &
                        i,lab,cstype,xfrac(nri),ocha(1),yfrac(nri),ocha(2),zfrac(nri),ocha(3),qa(i),occua(i),fixed
                    else
                      write(ioout,'(i7,1x,a5,1x,a2,2x,f9.6,2(1x,a1,1x,f9.6),1x,a1,1x,f9.5,f12.6,1x,a3)') &
                        i,lab,cstype,xafrac(i),ocha(1),yafrac(i),ocha(2),zafrac(i),ocha(3),qa(i),occua(i),fixed
                    endif
                  elseif (ndimen(icfg).eq.2) then
                    write(ioout,'(i7,1x,a5,1x,a2,2x,2(f9.6,1x,a1,1x),f9.4,1x,a1,1x,f9.5,f12.6,1x,a3)') &
                      i,lab,cstype,xafrac(i),ocha(1),yafrac(i),ocha(2),zafrac(i),ocha(3),qa(i),occua(i),fixed
                  elseif (ndimen(icfg).eq.1) then
                    write(ioout,'(i7,1x,a5,1x,a2,2x,f9.6,2(1x,a1,1x,f9.4),1x,a1,1x,f9.5,f12.6,1x,a3)') &
                      i,lab,cstype,xafrac(i),ocha(1),yafrac(i),ocha(2),zafrac(i),ocha(3),qa(i),occua(i),fixed
                  else
                    write(ioout,'(i7,1x,a5,1x,a2,2x,f9.4,2(1x,a1,1x,f9.4),1x,a1,1x,f9.5,f12.6,1x,a3)') &
                      i,lab,cstype,xafrac(i),ocha(1),yafrac(i),ocha(2),zafrac(i),ocha(3),qa(i),occua(i),fixed
                  endif
                endif
              endif
            endif
          enddo
          write(ioout,'(''--------------------------------------------------------------------------------'')')
        enddo
      endif
      write(ioout,'(/)')
      if (ndimen(icfg).eq.2) then
!
!  Check for growth slice output
!
        lslice = .false.
        i = 0
        do while (i.lt.nasym.and..not.lslice)
          i = i + 1
          lslice = lsliceatom(nsft + i)
        enddo
        if (lslice) then
          if (.not.lterseincoords) then
            write(ioout,'(/,''  Growth Slice : '',/)')
            write(ioout,'(''  Number of formula units in slice = '',i4,/)') nzmol
            write(ioout,'(''  Mixed fractional/Cartesian coordinates'',/)')
            write(ioout,'(''--------------------------------------------------------------------------------'')')
            write(ioout,'(''   No.  Atomic       x           y          z         Charge      Occupancy'')')
            write(ioout,'(''        Label      (Frac)      (Frac)     (Angs)        (e)         (Frac)  '')')
            write(ioout,'(''--------------------------------------------------------------------------------'')')
          endif
          do i = 1,nasym
            if (lsliceatom(nsft+i)) then
              inat = iatn(i)
              itype = natype(i)
!
!  Hide shells?
!
              if (inat.le.maxele.or..not.lhideshells) then
                call label(inat,itype,lab)
                if (lbsmat(i+nsft)) then
                  cstype = 'bc'
                  if (inat.gt.maxele) cstype = 'bs'
                else
                  cstype = 'c '
                  if (inat.gt.maxele) cstype = 's '
                endif
                if (.not.lterseincoords) then
                  write(ioout,'(i7,1x,a5,1x,a2,2x,2(f9.6,3x),f9.4,3x,f9.4,f12.6,1x)') &
                    i,lab,cstype,xafrac(i),yafrac(i),zafrac(i),qa(i),occua(i)
                endif
              endif
            endif
          enddo
          write(ioout,'(''--------------------------------------------------------------------------------'')')
          write(ioout,'(/)')
        endif
      endif
      if (ndimen(icfg).gt.0.and.index(keyword,'cart').ne.0) then
        write(ioout,'(''  Initial Cartesian coordinates :'',/)')
        write(ioout,'(''--------------------------------------------------------------------------------'')')
        write(ioout,'(''   No.  Atomic        x           y           z          Charge   Occupancy'')')
        write(ioout,'(''        Label       (Angs)      (Angs)      (Angs)        (e)       (Frac)  '')')
        write(ioout,'(''--------------------------------------------------------------------------------'')')
        do i = 1,numat
          inat = nat(i)
          itype = nftype(i)
!
!  Hide shells?
!
          if (inat.le.maxele.or..not.lhideshells) then
            call label(inat,itype,lab)
            if (lbsmat(nrelf2a(i)+nsft)) then
              cstype = 'bc'
              if (inat.gt.maxele) cstype = 'bs'
            else
              cstype = 'c '
              if (inat.gt.maxele) cstype = 's '
            endif
            write(ioout,'(i7,1x,a5,1x,a2,5f12.6)')i,lab,cstype,xclat(i),yclat(i),zclat(i),qf(i),occuf(i)
          endif
        enddo
        write(ioout,'(''--------------------------------------------------------------------------------'')')
        write(ioout,'(/)')
      endif
!
!  External potential output
!
      lfirstout = .true.
      do i = 1,nasym
        if (abs(extpotcfg(nsft+i)).gt.1.0d-6) then
          if (lfirstout) then
          lfirstout = .false.
          write(ioout,'(''  External electrostatic potential on asymmetric unit:'',/)')
          write(ioout,'(''--------------------------------------------------------------------------------'')')
          write(ioout,'(''   No.  Atomic         V        '')')
          write(ioout,'(''        Label        (eV/q)     '')')
          write(ioout,'(''--------------------------------------------------------------------------------'')')
          endif
          inat = iatn(i)
          itype = natype(i)
          call label(inat,itype,lab)
          if (lbsmat(i+nsft)) then
            cstype = 'bc'
            if (inat.gt.maxele) cstype = 'bs'
          else
            cstype = 'c '
            if (inat.gt.maxele) cstype = 's '
          endif
          write(ioout,'(i7,1x,a5,1x,a2,f14.6)') i,lab,cstype,extpotcfg(nsft+i)
        endif
      enddo
      if (.not.lfirstout) then
        write(ioout,'(''--------------------------------------------------------------------------------'',/)')
      endif
!
!  External force output
!
      lfirstout = .true.
      do i = 1,nasym
        forcenorm = abs(forcecfg(1,nsft+i)) + abs(forcecfg(2,nsft+i)) + abs(forcecfg(3,nsft+i))
        if (forcenorm.gt.1.0d-6) then
          if (lfirstout) then
          lfirstout = .false.
          write(ioout,'(''  External Cartesian forces on asymmetric unit:'',/)')
          write(ioout,'(''--------------------------------------------------------------------------------'')')
          write(ioout,'(''   No.  Atomic         Fx            Fy            Fz         '')')
          write(ioout,'(''        Label       (eV/Angs)     (eV/Angs)    (eV/Angs)      '')')
          write(ioout,'(''--------------------------------------------------------------------------------'')')
          endif
          inat = iatn(i)
          itype = natype(i)
          call label(inat,itype,lab)
          if (lbsmat(i+nsft)) then
            cstype = 'bc'
            if (inat.gt.maxele) cstype = 'bs'
          else
            cstype = 'c '
            if (inat.gt.maxele) cstype = 's '
          endif
          write(ioout,'(i7,1x,a5,1x,a2,3f14.6)')  &
            i,lab,cstype,forcecfg(1,nsft+i),forcecfg(2,nsft+i),forcecfg(3,nsft+i)
        endif
      enddo
      if (.not.lfirstout) then
        write(ioout,'(''--------------------------------------------------------------------------------'',/)')
      endif
!
!  Time-dependent external force output
!
      lfirstout = .true.
      do i = 1,nasym
        do j = 1,3
          if (ltdforcecfg(j,nsft+i)) then
            if (lfirstout) then
              lfirstout = .false.
              write(ioout,'(''  Time-dependent external Cartesian forces on asymmetric unit:'',/)')
              write(ioout,'(''--------------------------------------------------------------------------------'')')
              write(ioout,'(''   No.  Atomic   Direction       FA            FB             FC         '')')
              write(ioout,'(''        Label                 (eV/Angs)      (1/ps)      (Fraction 2xpi)      '')')
              write(ioout,'(''--------------------------------------------------------------------------------'')')
            endif
            inat = iatn(i)
            itype = natype(i)
!
!  Hide shells?
!
            if (inat.le.maxele.or..not.lhideshells) then
              call label(inat,itype,lab)
              if (lbsmat(i+nsft)) then
                cstype = 'bc'
                if (inat.gt.maxele) cstype = 'bs'
              else
                cstype = 'c '
                if (inat.gt.maxele) cstype = 's '
              endif
              write(ioout,'(i7,1x,a5,1x,a2,6x,a1,3x,3f14.6)')  &
                i,lab,cstype,crd1(j),tdforcecfg(1,j,nsft+i),tdforcecfg(2,j,nsft+i),tdforcecfg(3,j,nsft+i)
            endif
          endif
        enddo
      enddo
      if (.not.lfirstout) then
        write(ioout,'(''--------------------------------------------------------------------------------'',/)')
      endif
!
!  Shear force output
!
      if (lshearforcecfg(icfg)) then
        write(ioout,'(''  Shear force :'',/)')
        write(ioout,'(''--------------------------------------------------------------------------------'')')
        write(ioout,'(''  Force constant            Cartesian force direction / '')')
        write(ioout,'(''    (eV/Ang**2)             Cartesian force normal '')')
        write(ioout,'(''--------------------------------------------------------------------------------'')')
        write(ioout,'(2x,f12.6,10x,2(f10.7,1x),f10.7)') shearforcecfg(icfg),(shearforcedircfg(i,icfg),i=1,3)
        write(ioout,'(24x,2(f10.7,1x),f10.7)') (shearforcenormcfg(i,icfg),i=1,3)
        write(ioout,'(''--------------------------------------------------------------------------------'',/)')
      endif
!
!  Einstein model output
!
      lfirstout = .true.
      do i = 1,nasym
        if (leinsteinat(nsft+i)) then
          if (lfirstout) then
            lfirstout = .false.
            write(ioout,'(''  Einstein Model sites/force constants for asymmetric unit :'',/)')
            write(ioout,'(''--------------------------------------------------------------------------------'')')
            write(ioout,'(''   No.  Atomic         x             y            z         Force constant'')')
            write(ioout,'(''        Label        (frac)        (frac)       (frac)       (eV/Angs**2)'')')
            write(ioout,'(''--------------------------------------------------------------------------------'')')
          endif
          inat = iatn(i)
          itype = natype(i)
!
!  Hide shells?
!
          if (inat.le.maxele.or..not.lhideshells) then
            call label(inat,itype,lab)
            if (lbsmat(i+nsft)) then
              cstype = 'bc'
              if (inat.gt.maxele) cstype = 'bs'
            else
              cstype = 'c '
              if (inat.gt.maxele) cstype = 's '
            endif
            write(ioout,'(i7,1x,a5,1x,a2,4f14.6)') &
              i,lab,cstype,xeinsteinat(nsft+i),yeinsteinat(nsft+i),zeinsteinat(nsft+i),keinsteinat(nsft+i)
          endif
        endif
      enddo
      if (.not.lfirstout) then
        write(ioout,'(''--------------------------------------------------------------------------------'',/)')
      endif
    endif
!********************
!  Molecule output  *
!********************
    if (lmol.and.ioproc.and..not.lterseinmol) call outmol
!**********************
!  Constraint output  *
!**********************
    if (ncon.gt.0.and.ioproc) then
      write(ioout,'(''  Constraints : '',/)')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''  Constraint no.      Unconstrained     Constrained    Coefficient    Offset'')')
      write(ioout,'(''                         Variable         Variable'')')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      ncm = 3*nasym + nstrains
      do i = 1,ncon
        ncvi = ncvarindcfg(i+ncfst)
        if (ncvartypcfg(i+ncfst).eq.iopt_strain) then
          nfv = ncfixindcfg(i+ncfst)
          write(ioout,'(8x,i4,14x,''Strain '',i1,8x,''Strain '',i1,4x,f10.5,5x,f7.4)') &
            i,ncvi,nfv,concocfg(i+ncfst),conaddcfg(i+ncfst)
        elseif (ncvartypcfg(i+ncfst).eq.iopt_cell) then
          nfv = ncfixindcfg(i+ncfst)
          write(ioout,'(8x,i4,14x,''Cell   '',i1,8x,''Cell   '',i1,4x,f10.5,5x,f7.4)') &
            i,ncvi,nfv,concocfg(i+ncfst),conaddcfg(i+ncfst)
        elseif (ncvartypcfg(i+ncfst).eq.iopt_radius) then
          nfv = ncfixindcfg(i+ncfst)
          write(ioout,'(6x,i6,14x,''Radius '',i4,5x,''Radius '',i4,1x,f10.5,5x,f7.4)') &
            i,ncvi,nfv,concocfg(i+ncfst),conaddcfg(i+ncfst)
        else
          if (ncvartypcfg(i+ncfst).eq.iopt_xf) then
            ncrv = 1
          elseif (ncvartypcfg(i+ncfst).eq.iopt_yf) then
            ncrv = 2
          elseif (ncvartypcfg(i+ncfst).eq.iopt_zf) then
            ncrv = 3
          elseif (ncvartypcfg(i+ncfst).eq.iopt_xcom) then
            ncrv = 4
          elseif (ncvartypcfg(i+ncfst).eq.iopt_ycom) then
            ncrv = 5
          elseif (ncvartypcfg(i+ncfst).eq.iopt_zcom) then
            ncrv = 6
          elseif (ncvartypcfg(i+ncfst).eq.iopt_xqtn) then
            ncrv = 7
          elseif (ncvartypcfg(i+ncfst).eq.iopt_yqtn) then
            ncrv = 8
          elseif (ncvartypcfg(i+ncfst).eq.iopt_zqtn) then
            ncrv = 9
          else
            call outerror('undefined optimisation variable type in constraint',0_i4)
            call stopnow('setcfg')
          endif
          nfv = ncfixindcfg(i+ncfst)
          if (ncfixtypcfg(i+ncfst).eq.iopt_xf) then
            ncrf = 1
          elseif (ncfixtypcfg(i+ncfst).eq.iopt_yf) then
            ncrf = 2
          elseif (ncfixtypcfg(i+ncfst).eq.iopt_zf) then
            ncrf = 3
          elseif (ncfixtypcfg(i+ncfst).eq.iopt_xcom) then
            ncrf = 4
          elseif (ncfixtypcfg(i+ncfst).eq.iopt_ycom) then
            ncrf = 5
          elseif (ncfixtypcfg(i+ncfst).eq.iopt_zcom) then
            ncrf = 6
          elseif (ncfixtypcfg(i+ncfst).eq.iopt_xqtn) then
            ncrf = 7
          elseif (ncfixtypcfg(i+ncfst).eq.iopt_yqtn) then
            ncrf = 8
          elseif (ncfixtypcfg(i+ncfst).eq.iopt_zqtn) then
            ncrf = 9
          else
            call outerror('undefined optimisation variable type in constraint',0_i4)
            call stopnow('setcfg')
          endif
          write(ioout,'(6x,i6,12x,i6,1x,a4,5x,i6,1x,a4,3x,f10.5,5x,f7.4)') &
            i,ncvi,crd4(ncrv),nfv,crd4(ncrf),concocfg(i+ncfst),conaddcfg(i+ncfst)
        endif
      enddo
      write(ioout,'(''--------------------------------------------------------------------------------'',/)')
    endif
    if (ldist.or.lbond) then
      if (ndimen(icfg).eq.3) then
        do i = 1,numat
          xclat(i) = xfrac(i)*rv(1,1) + yfrac(i)*rv(1,2) + zfrac(i)*rv(1,3)
          yclat(i) = xfrac(i)*rv(2,1) + yfrac(i)*rv(2,2) + zfrac(i)*rv(2,3)
          zclat(i) = xfrac(i)*rv(3,1) + yfrac(i)*rv(3,2) + zfrac(i)*rv(3,3)
        enddo
      elseif (ndimen(icfg).eq.2) then
        do i = 1,numat
          xclat(i) = xfrac(i)*rv(1,1) + yfrac(i)*rv(1,2)
          yclat(i) = xfrac(i)*rv(2,1) + yfrac(i)*rv(2,2)
          zclat(i) = zfrac(i)
        enddo
      elseif (ndimen(icfg).eq.1) then
        do i = 1,numat
          xclat(i) = xfrac(i)*rv(1,1)
          yclat(i) = yfrac(i)
          zclat(i) = zfrac(i)
        enddo
      endif
      do i = 1,nasym
        nr = nrela2f(i)
        xalat(i) = xclat(nr)
        yalat(i) = yclat(nr)
        zalat(i) = zclat(nr)
      enddo
      if (ioproc) then
        if (lbond) call GULP_bond
        if (ldist) call distance
      endif
    endif
!
!  Option to print out valid angles
!
    if (ioproc) then
      call rlist
      if (nthb.gt.0.and.langle) call getangles(ioout,0_i4,nangtot)
      if (nfor.gt.0.and.ltors) call gettorsions(ioout,0_i4,nphitot)
    endif
!*********************************
!  Configuration k point setups  *
!*********************************
    if (ndline.gt.0) then
      call setdisp
    endif
    if (lshrink) then
      if (ndim.eq.3) then
        call setkpt3D
      elseif (ndim.eq.2) then
        call setkpt2D
      elseif (ndim.eq.1) then
        call setkpt1D
      endif
    endif
!
!  If there is more than one k point per configuration or there
!  are too many atoms then warn that no eigenvectors will be produced.
!
    if (leigen) then
      nlkpt = 0
      do i = 1,nkpt
        nk = nkptcfg(i)
        if (nlkpt.eq.0.and.nk.eq.icfg) nlkpt = i
      enddo
    endif
!
!  Set default k points if phonon specified but no k points given
!  Modified so that each configuration is checked individually.
!
    nlkpt = 0
    do j = 1,nkpt
      if (nkptcfg(j).eq.icfg) nlkpt = nlkpt + 1
    enddo
!
!  Create space in k point arrays for any new points
!
    if (nkpt+1.gt.maxkpt) then
      maxkpt = nkpt + 1
      call changemaxkpt
    endif
    if (nlkpt.eq.0.and.(lphon.or.lphonfit).and.ndline.eq.0) then
      if (nkpt.gt.0) then
        k = 0
        nkp = 0
        do while (k.lt.nkpt.and.nkp.lt.icfg)
          k = k + 1
          nkp = nkptcfg(k)
        enddo
        if (nkp.lt.icfg) then
          k = nkpt + 1
        else
          do j = nkpt,k,-1
            xkpt(j+1) = xkpt(j)
            ykpt(j+1) = ykpt(j)
            zkpt(j+1) = zkpt(j)
            wkpt(j+1) = wkpt(j)
            nkptcfg(j+1) = nkptcfg(j)
            lkptdispersion(j+1) = lkptdispersion(j)
          enddo
        endif
      else
        k = 1
      endif
      xkpt(k) = 0.0_dp
      ykpt(k) = 0.0_dp
      zkpt(k) = 0.0_dp
      wkpt(k) = 1.0_dp
      nkptcfg(k) = icfg
      lkptdispersion(k) = .false.
      nkpt = nkpt + 1
    elseif (nlkpt.eq.0.and.lfree) then
      if (nkpt.gt.0) then
        k = 0
        nkp = 0
        do while (k.lt.nkpt.and.nkp.lt.icfg)
          k = k + 1
          nkp = nkptcfg(k)
        enddo
        if (nkp.lt.icfg) then
          k = nkpt + 1
        else
          do j = nkpt,k,-1
            xkpt(j+1) = xkpt(j)
            ykpt(j+1) = ykpt(j)
            zkpt(j+1) = zkpt(j)
            wkpt(j+1) = wkpt(j)
            nkptcfg(j+1) = nkptcfg(j)
            lkptdispersion(j+1) = lkptdispersion(j)
          enddo
        endif
      else
        k = 1
      endif
      xkpt(k) = 0.25_dp
      if (ndimen(icfg).ge.2) then
        ykpt(k) = 0.25_dp
      else
        ykpt(k) = 0.0_dp
      endif
      if (ndimen(icfg).eq.3) then
        zkpt(k) = 0.25_dp
      else
        zkpt(k) = 0.0_dp
      endif
      wkpt(k) = 1.0_dp
      nkptcfg(k) = icfg
      lkptdispersion(k) = .false.
      nkpt = nkpt + 1
    endif
    if (nkpt.gt.0.and.index(keyword,'nokp').eq.0.and.ioproc) call outkpt
!**********************************************************
!  Include optimisation variables in fitting observables  *
!**********************************************************
    if (lfit) then
      if (lrelax) then
        nfst = n1var(icfg) - 1
        do i = 1,nvarcfg(icfg)
          ind = ioptindexcfg(i+nfst)
          itp = iopttypecfg(i+nfst)
!
!  For relax fitting shell positions must not be included in observables list
!
          lcore = .true.
          if (itp.eq.iopt_xf.or.itp.eq.iopt_yf.or.itp.eq.iopt_zf) then
!
!  Coordinate variable - exclude shells
!
            if (iatn(ind).gt.maxele) lcore = .false.
          elseif (itp.eq.iopt_radius) then
!
!  Variable is a radius or charge => exclude
!
            lcore = .false.
          endif
          if (lcore) then
            nobs = nobs + 1
            if (nobs.gt.maxobs) then
              maxobs = nobs + 20
              call changemaxobs
            endif
            nobtyp(nobs) = 6
            nobcfg(nobs) = icfg
            nobptr(nobs) = i
            if (itp.eq.iopt_cell.or.itp.eq.iopt_strain) then
              if (ind.eq.1) then
                fobs(nobs) = a
              elseif (ind.eq.2) then
                fobs(nobs) = b
              elseif (ind.eq.3) then
!
!  For rhombohedral space groups, input in the rhombohedral setting the
!  second cell observable must be alpha instead of c for relax fitting.
!  Also for 2-D systems then quantity is alpha as well.
!
                if (ifhr(icfg).eq.1.or.ndimen(icfg).eq.2) then
                  fobs(nobs) = alpha
                else
                  fobs(nobs) = c
                endif
              elseif (ind.eq.4) then
                fobs(nobs) = alpha
              elseif (ind.eq.5) then
                fobs(nobs) = beta
              else
                fobs(nobs) = gamma
              endif
              if (ind.le.3) then
                if (weight(nobs).eq.0.0) weight(nobs) = delwht_cell_length
              else
                if (weight(nobs).eq.0.0) weight(nobs) = delwht_cell_angle
              endif
            elseif (itp.eq.iopt_xf) then
              fobs(nobs) = xafrac(ind)
              if (ndimen(ncf).eq.3) then
                if (weight(nobs).eq.0.0) weight(nobs) = delwht_frac
              else
                if (weight(nobs).eq.0.0) weight(nobs) = delwht_coord
              endif
            elseif (itp.eq.iopt_yf) then
              fobs(nobs) = yafrac(ind)
              if (ndimen(ncf).eq.3) then
                if (weight(nobs).eq.0.0) weight(nobs) = delwht_frac
              else
                if (weight(nobs).eq.0.0) weight(nobs) = delwht_coord
              endif
            elseif (itp.eq.iopt_zf) then
              fobs(nobs) = zafrac(ind)
              if (ndimen(ncf).eq.3) then
                if (weight(nobs).eq.0.0) weight(nobs) = delwht_frac
              else
                if (weight(nobs).eq.0.0) weight(nobs) = delwht_coord
              endif
            elseif (itp.eq.iopt_radius) then
              fobs(nobs) = rada(ind)
              if (weight(nobs).eq.0.0) weight(nobs) = delwht_radius
            elseif (itp.eq.iopt_xcom) then
              fobs(nobs) = molcoma(1,ind)
              if (ndimen(ncf).eq.3) then
                if (weight(nobs).eq.0.0) weight(nobs) = delwht_frac
              else
                if (weight(nobs).eq.0.0) weight(nobs) = delwht_coord
              endif
            elseif (itp.eq.iopt_ycom) then
              fobs(nobs) = molcoma(2,ind)
              if (ndimen(ncf).eq.3) then
                if (weight(nobs).eq.0.0) weight(nobs) = delwht_frac
              else
                if (weight(nobs).eq.0.0) weight(nobs) = delwht_coord
              endif
            elseif (itp.eq.iopt_zcom) then
              fobs(nobs) = molcoma(3,ind)
              if (ndimen(ncf).eq.3) then
                if (weight(nobs).eq.0.0) weight(nobs) = delwht_frac
              else
                if (weight(nobs).eq.0.0) weight(nobs) = delwht_coord
              endif
            elseif (itp.eq.iopt_xqtn) then
              fobs(nobs) = molQa(1,ind)
              if (weight(nobs).eq.0.0) weight(nobs) = delwht_qtn
            elseif (itp.eq.iopt_yqtn) then
              fobs(nobs) = molQa(2,ind)
              if (weight(nobs).eq.0.0) weight(nobs) = delwht_qtn
            elseif (itp.eq.iopt_zqtn) then
              fobs(nobs) = molQa(3,ind)
              if (weight(nobs).eq.0.0) weight(nobs) = delwht_qtn
            endif
          endif
        enddo
      else
        nfst = n1var(icfg) - 1
        nobsold = nobs
        do i = 1,nvarcfg(icfg)
          ind = ioptindexcfg(i+nfst)
          itp = iopttypecfg(i+nfst)
!
!  Shell positions must not be included in observables list if the summing of core and shells is specified
!
          lcore = .true.
          if ((itp.eq.iopt_xf.or.itp.eq.iopt_yf.or.itp.eq.iopt_zf).and.lsumcoreshell(icfg)) then
!
!  Coordinate variable - exclude shells
!
            if (iatn(ind).gt.maxele) lcore = .false.
          elseif (itp.eq.iopt_radius.and.lsumcoreshell(icfg)) then
!
!  Variable is a radius or charge => exclude
!
            if (iatn(ind).gt.maxele) lcore = .false.
          endif
!
          if (lcore) then
            nobs = nobs + 1
            if (nobs.gt.maxobs) then
              maxobs = nobs + 20
              call changemaxobs
            endif
            fobs(nobs) = 0.0_dp
            nobtyp(nobs) = 2
            nobcfg(nobs) = icfg
            nobptr(nobs) = i
            if (weight(nobs).eq.0.0) weight(nobs) = delwht_grad
          endif
        enddo
!
!  Look for user specified gradients to fit
!
        nlfgrad = 0
        if (nfgrad.gt.0) then
          nlfgra = 0
          nufgra = 0
          do k = 1,nfgrad
            if (nfgracfg(k).eq.icfg) then
              nufgra = k
              if (nlfgra.eq.0) nlfgra = k
            endif
          enddo
          if (nlfgra.gt.0) nlfgrad = nufgra - nlfgra + 1
        endif
        if (nlfgrad.gt.0) then
          if (nobsold+nvarcfg(icfg).gt.maxobs) then
            maxobs = nobsold + nvarcfg(icfg) + 20
            call changemaxobs
          endif
          do i = 1,nvarcfg(icfg)
            ind = ioptindexcfg(i+nfst)
            if (iopttypecfg(i+nfst).eq.iopt_xf) then
              do k = nlfgra,nufgra
                if (nfgrat(k).eq.ind) then
                  fobs(nobsold+i) = fgrad(3*(k-1)+1)
                  weight(nobsold+i) = fgradweight(k)
                endif
              enddo
            elseif (iopttypecfg(i+nfst).eq.iopt_yf) then
              do k = nlfgra,nufgra
                if (nfgrat(k).eq.ind) then
                  fobs(nobsold+i) = fgrad(3*(k-1)+2)
                  weight(nobsold+i) = fgradweight(k)
                endif
              enddo
            elseif (iopttypecfg(i+nfst).eq.iopt_zf) then
              do k = nlfgra,nufgra
                if (nfgrat(k).eq.ind) then
                  fobs(nobsold+i) = fgrad(3*(k-1)+3)
                  weight(nobsold+i) = fgradweight(k)
                endif
              enddo
            endif
          enddo
        endif
!   
!  Look for user specified strain derivatives to fit
!   
        nlfstrain = 0
        if (nfstrain.gt.0) then
          nlfstr = 0
          nufstr = 0
          do k = 1,nfstrain 
            if (nfstraincfg(k).eq.icfg) then
              nufstr = k
              if (nlfstr.eq.0) nlfstr = k
            endif
          enddo 
          if (nlfstr.gt.0) nlfstrain = nufstr - nlfstr + 1
        endif 
        if (nlfstrain.gt.0) then
          if (nobsold+nvarcfg(icfg).gt.maxobs) then 
            maxobs = nobsold + nvarcfg(icfg) + 20
            call changemaxobs 
          endif
          do i = 1,nvarcfg(icfg)
            ind = ioptindexcfg(i+nfst)
            if (iopttypecfg(i+nfst).eq.iopt_strain.or.iopttypecfg(i+nfst).eq.iopt_cell) then
              do k = nlfstr,nufstr
                if (nfstraint(k).eq.ind) then
                  fobs(nobsold+i) = fstrain(k)
                  weight(nobsold+i) = fstrainweight(k)
                endif
              enddo 
            endif
          enddo 
        endif
      endif
    endif
!
!  Free array for optimisation logicals
!
    deallocate(loptz,stat=status)
    if (status/=0) call deallocate_error('setcfg','loptz')
    deallocate(lopty,stat=status)
    if (status/=0) call deallocate_error('setcfg','lopty')
    deallocate(loptx,stat=status)
    if (status/=0) call deallocate_error('setcfg','loptx')
    deallocate(loptr,stat=status)
    if (status/=0) call deallocate_error('setcfg','loptr')
    deallocate(loptm,stat=status)
    if (status/=0) call deallocate_error('setcfg','loptm')
!*****************************
!  End of configuration loop *
!*****************************
  enddo
!
!  OpenKIM potential setup - sets configuration information relating to KIM
!
  if (lkim_model) then
    call setkim
  endif
#ifdef TRACE
  call trace_out('setcfg')
#endif
!
  return
  end
