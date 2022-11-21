  subroutine setup(lall)
!
!  Set up each configuration
!
!   4/98 lra forced to be consistent with the space group
!   5/00 call to setup polarisability data added
!  11/00 arrays for structure prediction added (cn*/ox*)
!   4/02 Pointer to cores added
!  10/02 Storing of initial coordinates for external forces moved here
!  11/02 lewald flag turned on if EEM type method is being used
!  11/02 leinstein flag set to indicate presence of Einstein model
!   1/03 definition of lsymopt know includes space group
!   1/03 lewald flag turned off if lwolf is true
!   4/03 Check on sum of charges for individual regions
!   5/03 lspatialok initialised to .false. as default for config
!   6/03 Constraint handling corrected
!   4/04 lstr set to true if pressure file is to be written during MD
!   9/04 lewald set to true if there are any bond order charge potentials
!   9/04 Requirement for ndim > 0 for lewald to be set true removed
!   7/05 Streitz and Mintmire modifications added
!   7/05 Initialisation of EAM pointer call added
!   7/05 Copy atom to species number pointer to local array for configuration
!  12/05 Symmetry adapted derivative algorithms turned off for operator input
!  11/06 lfirst argument added to equpos call
!   5/07 Partial occupancy array initialisation added
!   5/07 Call to setspatial(bo) modified
!   6/07 Computation of region 1 species counters added
!  11/07 Modifications for reaxFF Coulomb term added
!   1/08 Forcing of lewald = .true. when lreaxFF is true removed
!   1/08 lreaxFFqreal removed
!   3/08 Array containing number of atoms of each species added
!   4/08 Turning off of lewald for lwolf = .true. case added back
!   4/08 Call to spatial decomposition version of setmol added
!  10/08 COSMIC setup added
!   7/09 cutoffmax(bo) removed as this is now passed via general module
!   5/10 Spatial decomposition turned off for Monte Carlo calculations
!   8/10 Spatial decomposition enabled for ReaxFF
!  10/10 Spatial decomposition enabled for EDIP model
!  11/10 Anisotropic pressure added
!   3/11 lstr is now true if lstressout is true
!   5/11 lstr now set to be true if this MD and the cell is 3-D so that pressure is computed correctly.
!   7/11 Spatial option added for partial occupancy pointer routines
!   8/11 Check on compatibility of cell optimisation and electric field added
!   5/12 Atomic stresses causes symmetry to be turned off for derivatives
!   9/12 Pacha added
!  12/12 Time-dependent field added
!  12/12 Modified to allow for multiple time-dependent fields
!   5/13 Set logical for control C during optimisation to be reinitialised
!   9/13 Call to setatomnoded2 added
!  10/14 nd2cell added
!  11/14 Call to setup reciprocal space group operators added
!   3/15 Gasteiger charges computed as part of set up since they are geometry independent
!   3/15 lnoqeem added
!   3/15 Gasteiger only called if lall is true since the bonding is required to get a
!        correct result and setmol is only called if ladd is true
!   9/16 ncoshptr added
!   1/17 Set parallel handling of variables 
!   3/17 fix_atom option added
!   3/17 Optimisation variable order modifications added
!   3/17 Call to setoptptr added
!   5/17 lsymderv2 now set to false for parallel runs
!   5/17 Argument added to setvarnoded2 call
!   6/17 Module files renamed to gulp_files
!   7/17 Call to setoptptr moved to optim since it needs lopf to be set first
!   8/17 Parallel decomposition now only computed for call where lall is true otherwise
!        blocksize trap can be applied to cell read in instead of any supercell
!   8/17 Order of calling setmol(s) and parallel initialisation switched to ensure correct operation.
!  10/17 Initial coordinates added to restart info
!  10/17 Absolute coordinates for MD now set here if required
!   1/18 Check on dimensionality for Grueneisen parameters added
!   1/18 Trace added
!   3/18 Option to call setframe added
!   5/18 Initialisation of nqrnow added
!   8/18 Strain for configuration set 
!   9/18 Strain inverse set up for lstraincell
!   9/18 Strain module added
!   9/18 Strain applied to cell for lstraincell case
!   9/18 nstrains2 added
!  12/18 Shear force added
!   2/19 Setting of shear force origin moved to force
!   3/19 iopt replaced by ioptindex and iopttype
!   3/19 Constraint handling now moved to subroutine
!   3/19 Multiple temperature ramps added
!   7/19 Symmetry adapted second derivatives turned off for polarisation
!   8/19 molatom option added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!  10/19 Setting of nasymnomol pointers added
!   2/20 lsymderv2 set to be false for rigid molecules
!   3/20 Setup of charge pointer added
!   3/20 Dielectric constant added
!   3/20 Molecule internal counters added
!   3/20 Fractional coordinates for dump added
!   4/20 Restarting for rigid molecules added
!   4/20 Trap for partial occupancies with rigid molecules added
!   6/20 Setting of nshellr1 modified for rigid molecules to only
!        include free shells
!   7/20 lpocc now set by calls to setoccptr routines
!  10/20 Correction to setting of ncharge for variable charge case
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
!  Julian Gale, CIC, Curtin University, October 2020
!
  use bondorderdata, only : nboQ, nbopot
  use control
  use configurations
  use current
  use datatypes
  use derivatives,   only : nd2cells, nd2cell, lfcsupercell
  use element,       only : lqeq, maxele, lSandM, lpacha, lgasteiger
  use eemdata,       only : nqrnow
  use field,         only : lfieldcfg, ntdfieldcfg
  use g_constants,   only : inverse_angstroms_to_ev
  use gulp_files,    only : lpre
  use interupt,      only : controlC_opt
  use iochannels
  use m_strain,      only : getinversestrain
  use mdlogic,       only : lmd
  use moldyn,        only : xabsco, yabsco, zabsco, labsco, labscoany
  use molecule
  use optimisation
  use parallel
  use partial
  use polarise
  use reallocate
  use reaxFFdata,    only : nreaxFFspec
  use shells
  use spatial,       only : lspatialok
  use spatialbo,     only : lspatialBOok
  use species
  use sutton,        only : lsuttonc
  use symmetry
  use thresholds,    only : thresh_q, thresh_c6
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  logical, intent(in)                          :: lall
!
!  Local variables
!
  character(len=1)                             :: cs
  character(len=5)                             :: lab
  integer(i4)                                  :: i
  integer(i4)                                  :: ii
  integer(i4)                                  :: inat
  integer(i4)                                  :: itype
  integer(i4)                                  :: j
  integer(i4)                                  :: nr
  integer(i4)                                  :: nspg
  logical                                      :: lp1
  logical                                      :: lqok
  logical                                      :: lvariablecharge
  real(dp)                                     :: asum
  real(dp)                                     :: c6trm
  real(dp)                                     :: qregion(2)
  real(dp)                                     :: qtrm
  real(dp)                                     :: sum
#ifdef TRACE
  call trace_in('setup')
#endif
!
!  Initialise flags
!
  lspatialok = .false.
  lspatialBOok = .false.
  lstr = .false.
  lp1 = (hmssg(1,ncf).eq.'P'.and.hmssg(3,ncf).eq.'1'.and.hmssg(4,ncf).eq.' '.and.hmssg(5,ncf).eq.' ')
  lsymopt = (lsymset(ncf).and.(nspcg(ncf).gt.1.or..not.lp1.or.ngocfg(ncf).gt.1)) 
  lqok = (index(keyword,'qok').ne.0)
  lvariablecharge = (((leem.or.lqeq.or.lSandM.or.lpacha).and.(.not.lnoqeem)).or.(nboQ.gt.0)) 
!
!  We can only use symmetry based derivative algorithms if the space group was input so that the
!  orientation of the cell is correctly handled.
!
  if (lsymopt.and.(lsymdok.and.nspcg(ncf).gt.1)) then
    lsymderv = .true.
    lsymderv2 = (nprocs.eq.1)
  else
    lsymderv = .false.
    lsymderv2 = .false.
  endif
!
!  Symmetry adapted second derivative algorithm cannot be used with variable charges or polarisation or rigid molecules
!
  if (index(keyword,'nod2').ne.0.or.leem.or.lpolar.or.lrigid) lsymderv2 = .false.
!
!  Free energy cannot use symmetry adapted derivatives yet
!
  if (lfree) then
    lsymderv = .false.
    lsymderv2 = .false.
  endif
!
!  Atomic stresses cannot use symmetry either
!
  if (latomicstress) then
    lsymderv = .false.
    lsymderv2 = .false.
  endif
!
!  Reinitialise the control C trap flag
!
  controlC_opt = .false.
!
!  Find first atom shift
!
  nsft = 0
  if (ncf.gt.1) then
    do i = 1,ncf-1
      nsft = nsft + nascfg(i)
    enddo
  endif
!
!  Set dimensionality
!
  ndim = ndimen(ncf)
!
!  Set number of strains according to the dimensionality and pointer
!
  if (ndim.eq.3) then
    nstrains = 6
    nstrptr(1) = 1
    nstrptr(2) = 2
    nstrptr(3) = 3
    nstrptr(4) = 4
    nstrptr(5) = 5
    nstrptr(6) = 6
  elseif (ndim.eq.2) then
    nstrains = 3
    nstrptr(1) = 1
    nstrptr(2) = 2
    nstrptr(3) = 6
  elseif (ndim.eq.1) then
    nstrains = 1
    nstrptr(1) = 1
  else
    nstrains = 0
  endif
  nstrains2 = nstrains*(nstrains + 1)/2
!
!  Set total charge
!
  totalcharge = totalchargecfg(ncf)
!
!  Set pressure and temperature
!
  press = presscfg(ncf)
  lanisotropicpress = lanisotropicpresscfg(ncf)
  if (lanisotropicpresscfg(ncf)) then
    anisotropicpress(1:6) = anisotropicpresscfg(1:6,ncf)
  else
    anisotropicpress(1:6) = 0.0_dp
  endif
  temperature = tempcfg(ncf)
  ntemperatureramp = ntempramp(ncf)
  temperaturestep(1:ntemperatureramp) = tempstp(1:ntemperatureramp,ncf)
  ntemperaturestep(1:ntemperatureramp) = ntempstp(1:ntemperatureramp,ncf)
!
!  Set start and stop steps for temperature ramps based on cumulative step count
!
  temperaturestart(1) = temperature
  ntemperaturestepstart(1) = ntempstpstart(1,ncf)
  ntemperaturestepstop(1)  = ntempstpstart(1,ncf) + ntemperaturestep(1)
  do i = 2,ntemperatureramp
    temperaturestart(i) = temperaturestart(i-1) + temperaturestep(i-1)*dble(ntemperaturestep(i-1))
    ntemperaturestepstart(i) = ntemperaturestepstop(i-1) + ntempstpstart(i,ncf)
    ntemperaturestepstop(i)  = ntemperaturestepstart(i) + ntemperaturestep(i)
  enddo
!
!  Set dielectric constant
!
  dielectric = dielectriccfg(ncf)
!
!  Check that this is not a defect run if dielectric constant is not 1
!
  if (ldefect.and.dielectric.gt.1.0_dp) then
    call outerror('dielectric constant must be one for defect calculations with Mott-Littleton',0_i4)
    call stopnow('setup')
  endif
!
!  Set Coulomb conversion factor for this configuration including the dielectric constant (usually just 1)
!
  angstoev = inverse_angstroms_to_ev/dielectric
!
!  Set stress and strain
!
  stress(1:6) = stresscfg(1:6,ncf)
  strain(1:6) = straincfg(1:6,ncf)
!
!  For lstraincell algorithm set up inverse of strain matrix
!
  if (lstraincell) then
    call getinversestrain(strain)
  endif
!
!  Set mode restrictions for free energy
!
  maxmode = maxmodecfg(ncf)
  minmode = minmodecfg(ncf)
  nummode = nummodecfg(ncf)
!
!  Set force constant supercell values
!
  nd2cell(1:3) = nd2cellcfg(1:3,ncf)
  nd2cells = (2*nd2cell(1) + 1)*(2*nd2cell(2) + 1)*(2*nd2cell(3) + 1)
  lfcsupercell = (nd2cells.gt.1) 
!
!  Transfer stored configuration data into working arrays
!
  nasym = nascfg(ncf)
  nbsmat = 0
  do i = 1,nasym
    iatn(i) = natcfg(nsft+i)
    natype(i) = ntypcfg(nsft+i)
    nspecptr(i) = nspecptrcfg(nsft+i)
    xafrac(i) = xcfg(nsft+i)
    yafrac(i) = ycfg(nsft+i)
    zafrac(i) = zcfg(nsft+i)
    qa(i) = qlcfg(nsft+i)
    occua(i) = occucfg(nsft+i)
    rada(i) = radcfg(nsft+i)
    if (lbsmat(nsft+i)) nbsmat = nbsmat + 1
    oxa(i) = oxcfg(nsft+i)
    cna(i) = cncfg(nsft+i)
  enddo
!
!  Set up constraint pointers
!
  if (ncf.eq.ncfg) then
    if (ncontot.gt.0) ncon = ncontot + 1 - n1con(ncf)
  else
    ncon = n1con(ncf+1) - n1con(ncf)
  endif
  if (ncon.gt.maxcon) then
    maxcon = ncon
    call changemaxcon
  endif
  ncfst = n1con(ncf) - 1
!**********************
!  Apply constraints  *
!**********************
  if (ncon.gt.0) then
    do j = 1,ncon
      ncfixind(j) = ncfixindcfg(ncfst+j)
      ncfixtyp(j) = ncfixtypcfg(ncfst+j)
      ncvarind(j) = ncvarindcfg(ncfst+j)
      ncvartyp(j) = ncvartypcfg(ncfst+j)
      conco(j) = concocfg(ncfst+j)
      conadd(j) = conaddcfg(ncfst+j)
    enddo
    call applyconstraints
  endif
!********************
!  Symmetry set up  *
!********************
  if (lsymopt) then
    call symmet
    call equpos(lall,.true.)
    call symmetp
  else
!***********************
!  No symmetry set up  *
!***********************
    ncbl = 1
    nccs = 1
    numat = nasym
    do i = 1,nasym
      qf(i) = qa(i)
      occuf(i) = occua(i)
      radf(i) = rada(i)
      nat(i) = iatn(i)
      nftype(i) = natype(i)
      neqv(i) = 1
      nrelf2a(i) = i
      nrela2f(i) = i
      xfrac(i) = xafrac(i)
      yfrac(i) = yafrac(i)
      zfrac(i) = zafrac(i)
      cnf(i) = cna(i)
      oxf(i) = oxa(i)
    enddo
  endif
!
!  Initialise fractional coordinates for dump to handle rigid molecules with straincell
!
  if (lrigid) then
    do i = 1,nasym
      xfdmp(i) = xfrac(i)
      yfdmp(i) = yfrac(i)
      zfdmp(i) = zfrac(i)
    enddo
  endif
!***************************
!  Multiple charge ranges  *
!***************************
  nqrnow(1:numat) = 1
!********************************
!  Count atoms of each species  *
!********************************
  numofspec(1:nspec) = 0
  do i = 1,nasym
    ii = nspecptr(i)
    numofspec(ii) = numofspec(ii) + neqv(i)
  enddo
!*********************
!  Set up unit cell  *
!*********************
  if (ndim.eq.3) then
    do i = 1,3
      rv(1,i) = rvcfg(1,i,ncf)
      rv(2,i) = rvcfg(2,i,ncf)
      rv(3,i) = rvcfg(3,i,ncf)
    enddo
    if (lstraincell) then
      call strain3D(strain,rv)
    endif
    call uncell3D(rv,a,b,c,alpha,beta,gamma)
    if (a.gt.1.0d-12) then
      recipa = 1.0_dp/a
    else
      recipa = 0.0_dp
    endif
    if (b.gt.1.0d-12) then
      recipb = 1.0_dp/b
    else
      recipb = 0.0_dp
    endif
    if (c.gt.1.0d-12) then
      recipc = 1.0_dp/c
    else
      recipc = 0.0_dp
    endif
    r1x = rv(1,1)
    r1y = rv(2,1)
    r1z = rv(3,1)
    r2x = rv(1,2)
    r2y = rv(2,2)
    r2z = rv(3,2)
    r3x = rv(1,3)
    r3y = rv(2,3)
    r3z = rv(3,3)
    sum = abs(r2x) + abs(r1y) + abs(r3x) + abs(r1z) + abs(r3y) + abs(r2z)
    lra = (sum.lt.1.0d-6)
!
!  Make sure lra is consistent with the space group
!
    nspg = nspcg(ncf)
    if (nspg.le.15.or.(nspg.ge.143.and.nspg.le.194)) then
      lra = .false.
    endif
  elseif (ndim.eq.2) then
    do i = 1,2
      rv(1,i) = rvcfg(1,i,ncf)
      rv(2,i) = rvcfg(2,i,ncf)
    enddo
    if (lstraincell) then
      call strain2D(strain,rv)
    endif
    call uncell2D(rv,a,b,alpha)
    if (a.gt.1.0d-12) then
      recipa = 1.0_dp/a
    else
      recipa = 0.0_dp
    endif
    if (b.gt.1.0d-12) then
      recipb = 1.0_dp/b
    else
      recipb = 0.0_dp
    endif
    r1x = rv(1,1)
    r1y = rv(2,1)
    r1z = 0.0_dp
    r2x = rv(1,2)
    r2y = rv(2,2)
    r2z = 0.0_dp
    r3x = 0.0_dp
    r3y = 0.0_dp
    r3z = 0.0_dp
    lra = (abs(alpha-90.0_dp).lt.1.0d-6)
  elseif (ndim.eq.1) then
    rv(1,1) = rvcfg(1,1,ncf)
    if (lstraincell) then
      call strain1D(strain,rv)
    endif
    call uncell1D(rv,a)
    if (a.gt.1.0d-12) then
      recipa = 1.0_dp/a
    else
      recipa = 0.0_dp
    endif
    r1x = rv(1,1)
    r1y = 0.0_dp
    r1z = 0.0_dp
    r2x = 0.0_dp
    r2y = 0.0_dp
    r2z = 0.0_dp
    r3x = 0.0_dp
    r3y = 0.0_dp
    r3z = 0.0_dp
    lra = .false.
  elseif (ndim.eq.0) then
    r1x = 0.0_dp
    r1y = 0.0_dp
    r1z = 0.0_dp
    r2x = 0.0_dp
    r2y = 0.0_dp
    r2z = 0.0_dp
    r3x = 0.0_dp
    r3y = 0.0_dp
    r3z = 0.0_dp
  endif
!
!  Setup cell vectors for neighbouring cells
!
  call rlist
!***********************************
!  Generate cartesian coordinates  *
!***********************************
  if (ndim.eq.3) then
    do i = 1,numat
      xclat(i) = xfrac(i)*r1x + yfrac(i)*r2x + zfrac(i)*r3x
      yclat(i) = xfrac(i)*r1y + yfrac(i)*r2y + zfrac(i)*r3y
      zclat(i) = xfrac(i)*r1z + yfrac(i)*r2z + zfrac(i)*r3z
    enddo
  elseif (ndim.eq.2) then
    do i = 1,numat
      xclat(i) = xfrac(i)*r1x + yfrac(i)*r2x
      yclat(i) = xfrac(i)*r1y + yfrac(i)*r2y
      zclat(i) = zfrac(i)
    enddo
  elseif (ndim.eq.1) then
    do i = 1,numat
      xclat(i) = xfrac(i)*r1x
      yclat(i) = yfrac(i)
      zclat(i) = zfrac(i)
    enddo
  else
    do i = 1,numat
      xclat(i) = xfrac(i)
      yclat(i) = yfrac(i)
      zclat(i) = zfrac(i)
    enddo
    if (lframe.and.lall) call setframe
  endif
  if (lsymopt) then
    do i = 1,nasym
      nr = nrela2f(i)
      xalat(i) = xclat(nr)
      yalat(i) = yclat(nr)
      zalat(i) = zclat(nr)
    enddo
  else
    do i = 1,nasym
      xalat(i) = xclat(i)
      yalat(i) = yclat(i)
      zalat(i) = zclat(i)
    enddo
  endif
!
!  For MD runs check whether absolute coordinates need to be substituted in
!
  if (lmd.and.labscoany) then
    do i = 1,nasym
      if (labsco(i)) then
        xalat(i) = xabsco(i)
        yalat(i) = yabsco(i)
        zalat(i) = zabsco(i)
      endif
    enddo
  endif
  if (linitcfg(ncf)) then
!
!  If initial coordinates were read in then transfer to current arrays
!
    do i = 1,nasym
      xinitial(i) = xinitcfg(nsft+i)
      yinitial(i) = yinitcfg(nsft+i)
      zinitial(i) = zinitcfg(nsft+i)
    enddo
  else
!
!  Save initial Cartesian coordinates in case external forces are applied
!
    do i = 1,nasym
      xinitial(i) = xalat(i)
      yinitial(i) = yalat(i)
      zinitial(i) = zalat(i)
    enddo
!
!  If field or delta_dipole is specified then output initial coordinates
!
    if (lfieldcfg(ncf).or.(ntdfieldcfg(ncf).gt.0).or.lddipole.or.lshearforcecfg(ncf)) then
      linitcfg(ncf) = .true.
      do i = 1,nasym
        xinitcfg(nsft+i) = xinitial(i)
        yinitcfg(nsft+i) = yinitial(i)
        zinitcfg(nsft+i) = zinitial(i)
      enddo
    endif
  endif
!
!  COSMO parameters
!
  if (lcosmo) call setcosmo
!******************************************************************
!  Check charge and Ewald summation flag : create charge pointer  *
!******************************************************************
  ncharge = 0
  sum = 0.0_dp
  asum = 0.0_dp
  qregion(1:2) = 0.0_dp
  do i = 1,numat
    qtrm = qf(i)*occuf(i)
    sum = sum + qtrm
    asum = asum + abs(qtrm)
    if (abs(qtrm).gt.thresh_q.or.lvariablecharge) then
      ncharge = ncharge + 1
      nchargeptr(ncharge) = i
    endif
    if (nregionno(nsft+nrelf2a(i)).eq.1) then
      qregion(1) = qregion(1) + qtrm
    else
      qregion(2) = qregion(2) + qtrm
    endif
  enddo
  if (.not.lqok.and.abs(sum).gt.1d-4.and.lall.and.ndim.ne.0) then
    if (ioproc) then
      write(ioout,'(/,''  Configuration number = '',i3,/)')ncf
      write(ioout,'(''  **** Unit cell is not charge neutral    ****'')')
      write(ioout,'(''  **** Sum of charges = '',f15.10,''   ****'')')sum
      write(ioout,'(''  **** Check that a special position atom ****'')')
      write(ioout,'(''  **** coordinate has not been varied     ****'')')
      write(ioout,'(/)')
      write(ioout,'(''  Current coordinates : '',/)')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''     No.  Atomic         x            y            z          Charge  Occupancy'')')
      if (ndim.eq.3) then
        write(ioout,'(''          Number       (frac)       (frac)       (frac)         (e)  '')')
      elseif (ndim.eq.2) then
        write(ioout,'(''          Number       (frac)       (frac)       (Angs)         (e)  '')')
      elseif (ndim.eq.1) then
        write(ioout,'(''          Number       (frac)       (Angs)       (Angs)         (e)  '')')
      endif
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      do i = 1,numat
        inat = nat(i)
        itype = nftype(i)
        cs = 'c'
        if (inat.gt.maxele) cs = 's'
        call label(inat,itype,lab)
        if (ndim.eq.3) then
          write(ioout,'(1x,i6,2x,a5,1x,a1,1x,3(4x,f9.6),2x,f12.6,2x,f6.4)') &
            i,lab,cs,xfrac(i),yfrac(i),zfrac(i),qf(i),occuf(i)
        elseif (ndim.eq.2) then
          write(ioout,'(1x,i6,2x,a5,1x,a1,1x,2(4x,f9.6),4x,f9.4,2x,f12.6,2x,f6.4)') &
            i,lab,cs,xfrac(i),yfrac(i),zclat(i),qf(i),occuf(i)
        elseif (ndim.eq.1) then
          write(ioout,'(1x,i6,2x,a5,1x,a1,1x,4x,f9.6,2(4x,f9.4),2x,f12.6,2x,f6.4)') &
            i,lab,cs,xfrac(i),yfrac(i),zclat(i),qf(i),occuf(i)
        endif
      enddo
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(/)')
    endif
    call stopnow('setup')
  endif
!
!  Set flag as to whether Ewald sum is needed
!
  lewald = ((asum.ne.0.0_dp.or.(lc6.and.ndim.eq.3)).and.index(keyword,'noel').eq.0)
  if (lvariablecharge) lewald = .true.
  if (lwolf) lewald = .false.
!
!  Set Ewald parameters for system if required
!
  if (lewald.and.ndim.gt.1) call setewald
!
!  Check on region charge sum
!
  if (.not.lqok.and.abs(qregion(1)).gt.1d-4.and.lall.and.ndim.ne.0) then
    call outerror('region 1 is not charge neutral',0_i4)
    call stopnow('setup')
  elseif (.not.lqok.and.abs(qregion(2)).gt.1d-4.and.lall.and.ndim.ne.0) then
    call outerror('region 2 is not charge neutral',0_i4)
    call stopnow('setup')
  endif
!
!  Set up spatial decomposition if needed
!
  lspatialok = (lspatial.and..not.lmc)
  lspatialBOok = (lspatial.and.(lbrenner.or.lEDIP.or.(nboQ+nbopot).gt.0.or.nreaxFFspec.gt.0).and..not.lmc)
  if (lspatialok) then
    call setcutoffmax
    call setspatial(.true.)
  endif
  if (lspatialBOok) then
    call setcutoffmaxbo
    call setspatialbo(.true.)
  endif
!
!  Set flag for Einstein model
!
  leinstein = .false.
  do i = 1,nasym
    if (leinsteinat(nsft+i)) leinstein = .true.
  enddo
!
!  Set up one centre C terms
!
  if (lall.and.lc6one) then
    do i = 1,nasym
      do j = 1,nspec
        if (iatn(i).eq.natspec(j).and.(natype(i).eq.ntypspec(j).or.ntypspec(j).eq.0)) then
          c6a(i) = c6spec(j)
        endif
      enddo
    enddo
    nchargec6 = 0
    do i = 1,numat
      qtrm = abs(qf(i)*occuf(i))
      do j = 1,nspec
        if (nat(i).eq.natspec(j).and.(nftype(i).eq.ntypspec(j).or.ntypspec(j).eq.0)) then
          c6f(i) = c6spec(j)
        endif
      enddo
      c6trm = abs(c6f(i)*occuf(i))
      if (qtrm.gt.thresh_q.or.c6trm.gt.thresh_c6) then
        nchargec6 = nchargec6 + 1
        nchargec6ptr(nchargec6) = i
      endif
    enddo
  else
    nchargec6 = ncharge
    nchargec6ptr(1:ncharge) = nchargeptr(1:ncharge)
  endif
!
!  Set up shell pointer array 
!
  ncore  = 0
  nshell = 0
  do i = 1,numat
    if (nat(i).gt.maxele) then
      nshell = nshell + 1
      nshptr(nshell) = i
      ncoshptr(i) = nshell
    else
      ncore = ncore + 1
      ncoptr(ncore) = i
      ncoshptr(i) = ncore
    endif
  enddo
!
!  EAM set up if needed
!
  if (lsuttonc) call seteam
!
!  Polarisability set up if needed
!
  if (lpolar) call setpolar
!
!  Core-shell pair check
!
  if (lall) call cscheck
! 
!  Set breathing shell pointer
! 
  call setbsmptr(nbs,nbss,nbsptr,(ndim.eq.2))
! 
!  Set partial occupancy pointer
!
  if (lspatialok) then
    call setoccshptrs(nsfoc,nbsfoc,iocshptr,ibocshptr,(ndim.eq.2))
  else
    call setoccshptr(nsfoc,nbsfoc,iocshptr,ibocshptr,(ndim.eq.2))
  endif
  if (lspatialok) then
    call setoccptrs(ncfoc,nsfoc,nbfoc,iocptr,ibocptr,(ndim.eq.2),lpocc)
  else
    call setoccptr(ncfoc,nsfoc,nbfoc,iocptr,ibocptr,(ndim.eq.2),lpocc)
  endif
  ncsfoc = ncfoc + nsfoc
!
!  Rigid molecules currently don't work with partial occupancy and so check for this
!
  if (lpocc.and.lrigid) then
    call outerror('partial occupancies not current supported with rigid molecules',0_i4)
    call stopnow('setup')
  endif
!
  if (lall) then
!
!  Calculate parallel division of work :
!   
!  Spatial - divide cells over processors
!  Non-spatial - divide atoms over processors
!
    call setatomnoded2
    call setatomdistribution('a')
    call setatomnodes(numat,nprocs,procid,lspatialok)
    call setatomnodesbo(numat,nprocs,procid,lspatialBOok)
  endif
!
!  If molecule option is selected call molecule setup routine - must come after parallel setup
!
  if (lall.and.index(keyword,'mol').ne.0) then
    if (lmolatom.or.lmolrigid) then
      call setmola
    elseif (lspatialok) then
      call setmols
    else
      call setmol
    endif
  else
!
!  If molecule is not being called then set nasymnomol/numatnomol/ncorenomol
!
    nasymnomol = nasym
    numatnomol = numat
    ncorenomol = ncore
    nbsmatnomol = nbsmat
    do i = 1,nasymnomol
      nasymnomolptr(i) = i
      nasymnomolrptr(i) = i
    enddo
    do i = 1,numatnomol
      numatnomolptr(i) = i
      numatnomolrptr(i) = i
    enddo
    do i = 1,ncorenomol
      ncorenomolptr(i) = i
      ncorenomolrptr(i) = i
    enddo
  endif
!
!  Set up shell pointer array for region 1 for properties section
!
  ncorer1  = 0
  nshellr1 = 0
  if (lrigid) then
    do i = 1,numatnomol
      ii = numatnomolptr(i)
      if (nregionno(nsft+nrelf2a(ii)).eq.1) then
        if (nat(ii).gt.maxele) then
          nshellr1 = nshellr1 + 1
        else
          ncorer1 = ncorer1 + 1
        endif
      endif
    enddo
  else
    do i = 1,numat
      if (nregionno(nsft+nrelf2a(i)).eq.1) then
        if (nat(i).gt.maxele) then
          nshellr1 = nshellr1 + 1
        else
          ncorer1 = ncorer1 + 1
        endif
      endif
    enddo
  endif
!
!  The values of ncorer1 and nshellr1 need to contain the number of species in region 1 only
!  for the 2-D case. For other dimensionalities, set equal to the total number of species
!
  if (ndim.ne.2) then
    if (lrigid) then
      ncorer1 = ncorenomol
    else
      ncorer1 = ncore
    endif
    nshellr1 = nshell
  endif
!
!  Set up variables for optimisation and fitting
!
  nvar = nvarcfg(ncf)
  nfst = n1var(ncf) - 1
  nbsm = 0
  nfixatom = nfixatomcfg(ncf)
  nfixatomtype = nfixatomtypecfg(ncf)
!
  ncellmax = ncellmaxcfg(ncf)
  ncellmin = ncellmincfg(ncf)
  if (ncellmax.gt.0) then
    ncell = ncellmax - ncellmin + 1
  else
    ncell = 0
  endif
!
  ninternalmax = ninternalmaxcfg(ncf)
  ninternalmin = ninternalmincfg(ncf)
  if (ninternalmax.gt.0) then
    ninternal = ninternalmax - ninternalmin + 1
  else
    ninternal = 0
  endif
!
  do i = 1,nvar
    ioptindex(i) = ioptindexcfg(i+nfst)
    iopttype(i) = iopttypecfg(i+nfst)
!
    if (iopttype(i).eq.iopt_radius) then
      nbsm = nbsm + 1
    endif
  enddo
!
!  Set number of variables of each type if rigid molecules are present
!
  if (lrigid) then
    ninternalatm = 0
    ninternalmolT = 0
    ninternalmolQ = 0
    do i = 1,nvar
      if (iopttype(i).eq.iopt_xcom) then
        ninternalmolT = ninternalmolT + 1
      elseif (iopttype(i).eq.iopt_ycom) then
        ninternalmolT = ninternalmolT + 1
      elseif (iopttype(i).eq.iopt_zcom) then
        ninternalmolT = ninternalmolT + 1
      elseif (iopttype(i).eq.iopt_xqtn) then
        ninternalmolQ = ninternalmolQ + 1
      elseif (iopttype(i).eq.iopt_yqtn) then
        ninternalmolQ = ninternalmolQ + 1
      elseif (iopttype(i).eq.iopt_zqtn) then
        ninternalmolQ = ninternalmolQ + 1
      endif
    enddo
    ninternalmol = ninternalmolT + ninternalmolQ
    ninternalatm = ninternal - ninternalmol
  else
    ninternalatm = ninternal
    ninternalmol = 0
    ninternalmolT = 0
    ninternalmolQ = 0
  endif
!
  if (lall) then
!
!  Set up parallel distribution of variables
!
    call setvarnoded2(.false.)
  endif
!
!  Check that no cell variables are specified with Einstein model
!
  if (leinstein.and.ncell.gt.0) then
    call outerror('cell must be fixed with Einstein model',0_i4)
    call stopnow('setup')
  endif
!
!  Check that no cell variables are specified with electric field
!
  if ((lfieldcfg(ncf).or.(ntdfieldcfg(ncf)).gt.0).and.ncell.gt.0) then
    call outerror('cell must be fixed with an electric field',0_i4)
    call stopnow('setup')
  endif
!
!  Check that Grueneisen is compatible with cell
!
  if (lgrueneisen.and.ndim.ne.3) then
    call outerror('Grueneisen parameters can only be calculated for 3-D systems',0_i4)
    call stopnow('setup')
  endif
!
!  Setup Gasteiger charges if needed since they are not geometry dependent
!
  if (lall.and.lgasteiger) then
    call gasteiger(.false.)
  endif
!
!  If any cell variables are to be optimised or a property calculation performed
!  then lstr must be true. Alternatively, if the pressure file is to be generated
!  in an MD run. Also, if the stress tensor is to be output then we must do the
!  strain derivatives. Furthermore, if this is a 3-D system and MD then set lstr
!  to be true so that the pressure is computed correctly.
!
  lstr = (ncell.gt.0.or.lprop.or.(lmd.and.(lpre.or.ndim.eq.3)).or.lstressout)
!
!  Surface energy flag - if cell is being optimised, don't do surface energy
!
  lseok = ((ndim.eq.1.or.ndim.eq.2).and.sbulkecfg(ncf).ne.0.0_dp.and.ncell.eq.0)
  nzmol = nzmolcfg(ncf)
#ifdef TRACE
  call trace_out('setup')
#endif
!
  return
  end
