  subroutine energy(etot,lgrad1,lgrad2)
!
!  Subroutine for calculating lattice energy
!
!   5/95 Modifications added for symmetrisation of second derivatives
!   6/95 Calls of real space routines stopped if there are no charges
!        or potentials
!   4/97 Sutton-Chen many-body potentials added
!   7/97 Neutralising background added
!  10/97 Modification added to allow for frozen atoms during zeroing 
!        of derv2
!  11/97 Bug in zeroing of derv2 for frozen atoms fixed
!  12/97 Storage of energy components added
!  12/97 Calls to kindex modified for variable rspeed
!  12/97 Self energy term added to total for EEM/QEq
!  12/97 Zeroing of first derivatives moved to before call to eem
!        so that contributions can be set in eem if needed
!   8/98 Parts made into subroutines
!   8/98 Strfin call moved outside energy for benefit of FEM
!   3/99 Parallel modifications added
!   7/99 Option to use minimum image for large calculations added
!   5/00 Dipolar polarisation energy from point ion added
!   6/00 recipsd removed to reduce size of code
!   3/01 Calculation of surface energy added
!   9/01 Modified for 1-D systems
!  11/01 Attachment energy added
!   5/02 Scaling of shift added
!   5/02 Brenner potential added
!   8/02 Surface energy calculation removed since this is now done
!        in a separate subroutine
!   8/02 External force added
!  10/02 Interaction energy between region 1/2 corrected
!  11/02 Einstein energy added
!  11/02 Call to psumall added
!   1/03 Wolf sum modifications added
!   5/03 Spatial decomposition option introduced
!   6/03 XML modifications added
!  11/03 Bond order potentials added
!   9/04 Dual spatial decomposition introduced
!  11/04 Six-body potentials added
!   5/06 Call to real1Dmd added
!   7/06 Six-body potentials added
!   8/06 Vibrational energy initialised though not used
!  11/06 Celltype called to ensure that ictype is correct
!   3/07 Electric field option added
!   3/07 Radial force added
!   3/07 Chemshell changes added
!   5/07 MC call option added for x0tostr
!   5/07 Call to setspatial(bo) modified
!   7/07 emeta added as dummy
!   7/07 ReaxFF calls added
!   7/07 Plane potential added for 2-D case
!   8/07 Plane potential enabled for O-D case
!  11/07 Unused variables removed
!  11/07 Modifications for reaxFF Coulomb term added
!   1/08 lreaxFFqreal removed
!  10/08 COSMO/COSMIC changes merged
!  12/08 scrho, and defect versions, converted to 2-D array for MEAM
!   2/09 Old xml calls removed
!   3/09 Call to sumderv1 added to handle non-radial derivative arrays
!   4/09 Minimum image option added for many in non-symmetry case
!   4/09 Separate call to generate MEAM densities added
!   6/09 Initialisation of site energies added
!   6/09 Module name changed from three to m_three
!   7/09 Experimental call to EVB added
!   7/09 cutoffmax/cutoffmaxbo moved to general module
!   7/09 Call to setcutoffmax(bo) removed
!   1/10 One-body potentials added
!   6/10 Symmetrised calculation of EAM/MEAM density turned off when lsymderv is false
!   9/10 EDIP energy added
!  10/10 EDIP enabled for spatial decomposition
!  11/10 Anisotropic pressure added
!   7/11 esregion2 added to call for wolfself
!   9/11 lgrad1 argument added to eem call
!   9/11 Metadynamics internal code replaced with Plumed
!  11/11 Initialisation of eregion2region added
!  11/11 Summing of off-diagonal region-region energies added
!  10/12 Support for OpenKIM models added
!   7/13 Improper torsion energy added
!   9/13 Calls to parallel second derivative subroutines added
!   9/13 Use of recip3Dmd introduced
!   8/14 eatom now passed to density routines
!   2/15 Criterion for calling getBOcharge changed by adding nboQ0
!   9/15 Spatial option added for 0-D case
!  12/15 Argument to setspatial subroutines changed to false otherwise finite difference calls fail.
!   1/16 Handling of symmetry logic in calling bondorder / brenner routines
!        corrected as numerical second derivatives were wrong with symmetry
!   3/16 eplumed added
!   6/16 SPME added
!   9/16 lgrad2 now passed to psumall
!   1/17 Call to fournosd added
!   1/17 Call to sixnosd added
!   1/17 Call to many0dd added
!   2/17 Call to getBOcharged added
!   4/17 ChemShell interaction modified
!   5/17 Modified to allow use of symmetry when lsymderv2 is false
!   6/17 Old ChemShell restored as a compile option
!   6/17 Call to psumall all modified so that lsymderv is only passed as true
!        if second derivatives were not calculated
!   7/17 Calling sequence for real/recip routines changed so that distributed
!        memory parallel versions are called for first derivatives too if
!        lDoQDeriv1 is true.
!   7/17 Parallel distributed memory calls added for minimium image
!   7/17 Use of distributed memory charge derivatives disabled for spatial
!        and cell multipole moment algorithms
!   7/17 threenosd now called if charge derivatives are needed with distributed
!        memory parallel algorithm
!   8/17 echemsh explicitly initialised in all runs
!  10/17 fhenergy initialised
!  10/17 Delta_dipole added
!   1/18 Trace added
!   5/18 Split bond EEM added
!   7/18 Call to sumderv1 corrected for lsymderv 
!   8/18 Modified for version 2 of OpenKIM
!   9/18 Change for multiple models in OpenKIM
!   9/18 Region and attachment energies added to kimmd arguments
!  11/18 Call to sumderv1 removed as no longer need
!  12/18 lgrad2 added as an argument to force
!   2/19 Rigid molecules added
!   4/19 Second derivatives of polarisation added
!   5/19 Calls to density3 routines modified
!   7/19 Symmetry handling changed for polarisation
!   7/19 lspatialBOok set to false for second derivatives in parallel
!        as algorithms are currently incompatible
!  11/19 Rigid molecule derivatives only evaluated if lrigid is true
!   2/20 SPME turned off for rigid molecules as this is not implemented yet
!   2/20 Spatial algorithms now turned off if symmetry is being used
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
!  Julian Gale, CIC, Curtin University, February 2020
!
  use bondorderdata, only : nbopot, nboQ, nboQ0
  use cellmultipole
  use configurations
  use g_constants
  use control
  use current
  use eam,           only : lMEAM, lMEAMden, maxmeamcomponent
  use energies
  use field,         only : lfieldcfg
  use four
#ifdef OLDCS
  use gulpchemsh,    only : ichemsh_qm
#endif
  use iochannels
  use kim_models,    only : lkim_model, nkimmodel
  use kspace
  use m_three
  use molecule
  use one,           only : none
  use optimisation
  use plane,         only : nplanepot
  use parallel
  use polarise
  use potentialxyz
  use reaxFFdata,    only : nreaxFFspec
  use shifts
  use six
  use spatial,       only : lspatialok
  use spatialbo,     only : lspatialBOok
  use spme,          only : lspme
  use sutton
  use symmetry
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  use two
  implicit none
!
!  Passed variables
!
  real(dp), intent(out) :: etot
  logical,  intent(in)  :: lgrad1
  logical,  intent(in)  :: lgrad2
!
!  Local variables
!
  integer(i4)           :: i
  integer(i4)           :: j
  integer(i4),     save :: ncflast = 0
  logical               :: lcmm
  logical               :: lpress
  logical               :: lspmeloc            ! If true, use of SPME is OK
  logical               :: lsymoffloc
  logical               :: lsymnotused
  logical               :: luseDMalg           ! If true, use distributed memory algorithm
  real(dp)              :: esum
  real(dp)              :: sft
#ifdef TRACE
  call trace_in('energy')
#endif
!
  lcmm = (icmm.gt.0.and.ndim.eq.0)
  lpress = (abs(press).gt.0.0_dp.and.ndim.gt.0)
  lsymoffloc = lsymoff
  lspmeloc = lspme
!
!  Initialise energy components
!
  etot = 0.0_dp
  evib = 0.0_dp
  ereal = 0.0_dp
  erecip = 0.0_dp
  erelax = 0.0_dp
  eatom = 0.0_dp
  ethb = 0.0_dp
  eedip = 0.0_dp
  eeinstein = 0.0_dp
  efield = 0.0_dp
  efor = 0.0_dp
  eforce = 0.0_dp
  ekim = 0.0_dp
  emany = 0.0_dp
  eimp = 0.0_dp
  eoop = 0.0_dp
  epolar = 0.0_dp
  epv = 0.0_dp
  ebrenner = 0.0_dp
  ecmm = 0.0_dp
  ec6 = 0.0_dp
  edipole = 0.0_dp
  ebgd = 0.0_dp
  emad = 0.0_dp
  eself = 0.0_dp
  esix = 0.0_dp
  eone = 0.0_dp
  eqeq = 0.0_dp
  eattach = 0.0_dp
  esurface = 0.0_dp
  esregion12 = 0.0_dp
  esregion2 = 0.0_dp
  ewolfself = 0.0_dp
  ebondorder = 0.0_dp
  eboQself = 0.0_dp
  echargecoupled = 0.0_dp
  eradial = 0.0_dp
  ereaxFF = 0.0_dp
  eplane  = 0.0_dp
  ecosmo  = 0.0_dp
  echemsh = 0.0_dp
  eplumed = 0.0_dp
!
  fhenergy = 0.0_dp
!
!  Initialise region-region interaction energy components and site energies
!
  eregion2region(1:nregions(ncf),1:nregions(ncf)) = 0.0_dp
  siteenergy(1:numat) = 0.0_dp
!************************************
!  Initialise density for MEAM/EAM  *
!************************************
  if (lsuttonc) then
    if (lMEAM) then
      do i = 1,numat
        scrho(1:maxmeamcomponent,i) = 0.0_dp
        scrho12(1:maxmeamcomponent,i) = 0.0_dp
      enddo
    else
      do i = 1,numat
        scrho(1,i) = 0.0_dp
        scrho12(1,i) = 0.0_dp
      enddo
    endif
  endif
!*******************************************
!  Set up spatial decomposition if needed  *
!*******************************************
  lspatialok = lspatial
  lspatialBOok = (lspatial.and.(lbrenner.or.lEDIP.or.(nbopot+nboQ).gt.0.or.nreaxFFspec.gt.0))
!
!  Symmetry algorithms don't use spatial and so turn off
!
  if (lsymopt) then
    lspatialok = .false.
    lspatialBOok = .false.
  endif
!
!  Second derivative algorithms are currently incompatible with spatial in parallel for bond order potentials
!
  if (lgrad2.and.nprocs.gt.1) lspatialBOok = .false.
!
  if (lspatialok) then
    call setspatial(.false.)    ! Argument has to be false otherwise fails for phonon finite difference
  endif
  if (lspatialBOok) then
    call setspatialbo(.false.)  ! Argument has to be false otherwise fails for phonon finite difference
  endif
!
!  Calculate parallel division of work :
!  
!  Spatial - divide cells over processors     
!  Non-spatial - divide atoms over processors
!   
  call setatomnodes(numat,nprocs,procid,lspatialok)
  call setatomnodesbo(numat,nprocs,procid,lspatialBOok)
!
!  Set up local variables
!
  if (ndim.eq.3) then
!
!  Symmetry can only be used if the cell parameters conform
!  to the correct space group - perform continuous check here
!
    if (lsymopt.and..not.lsymoff) then
      call celltype(ictype,icfhr)
      lsymoffloc = (nccs.gt.ictype)
      if (index(keyword,'verb').ne.0.and.lsymoffloc) then
        if (ioproc) then
          write(ioout,'(''  ** Symmetry turned off due to cell **'')')
        endif
      endif
    endif
  endif
!
!  Set flag as to whether symmetry should not be used even if present
!
  lsymnotused = ((lgrad2.and..not.lsymderv2).or.lsymoffloc) 
!***********************************************
!  Select choice of distributed memory or not  *
!***********************************************
!
!  Currently the following algorithms are not supported with distributed memory
!
!  - spatial decomposition
!  - cell multipole moments
!  - variable charges with second derivatives
!
  luseDMalg = (nprocs.gt.1.and.(lgrad2.or.(lDoQDeriv1.and..not.lspatialok.and..not.lcmm)))
!*************************
!  SPME algorithm check  *
!*************************
  if (lDoQDeriv1.or.latomicstress.or.lsiteenergy.or.lrigid) lspmeloc = .false.
!**************************************
!  Call EEM/QEq to calculate charges  *
!**************************************
  if (leem) then
    if (luseDMalg) then
      if (leembond) then
        call eemsplitd(.false.,lgrad1,lgrad2)
      else
        call eemdm(.false.,lgrad1,lgrad2)
      endif
    else
      if (leembond) then
        call eemsplit(.false.,lgrad1,lgrad2)
      else
        call eem(.false.,lgrad1,lgrad2)
      endif
    endif
  endif
!****************************************************
!  Calculate charges according to bond order model  *
!****************************************************
  if ((nboQ+nboQ0).gt.0) then
    if (luseDMalg) then
      call getBOcharged(lgrad1,lgrad2)
    else
      call getBOcharge(lgrad1,lgrad2)
    endif
  endif
!*********************
!  Solvation set up  *
!*********************
  if (lcosmo) then
    if (ncf.ne.ncflast.or.index(keyword,'nosa').eq.0) then
      call setsas
    else
      call updatesas
    endif
  endif
!*************************************************
!  Zero derivatives                              *
!  In ChemShell case the QM force is added here  *
!*************************************************
  call initdervs(lgrad1,lgrad2)
!**************************************
!  Initialise electric field to zero  *
!**************************************
  if (lpolar) then
!
!  Initialise field areas based on full cell since use of symmetry can vary with algorithm
!
    do i = 1,numat
      vx(i) = 0.0_dp
      vy(i) = 0.0_dp
      vz(i) = 0.0_dp
      vx12(i) = 0.0_dp
      vy12(i) = 0.0_dp
      vz12(i) = 0.0_dp
    enddo
    if (lqpolar) then
      do i = 1,numat
        v2xyz(1:6,i) = 0.0_dp
        v2xyz12(1:6,i) = 0.0_dp
      enddo
    endif
  endif
!***************************************
!  Cell Multipole Method for clusters  *
!***************************************
  if (lcmm) call setcmm
!
!  Redetermine cell indices for molecule atoms in case one has moved across cell boundary.
!
  if (nmol.gt.0.and.ndim.gt.0) then
    if (lmolfix) then
      call molindfix
    else
      call molind
!
!  Recalculate bond increment charges since they depend on the connectivity
!
      if (lqbond) call bondq(.false.)
    endif
  endif
!*******************************
!  Electrostatic contribution  *
!*******************************
  if (lewald.and.ndim.gt.1) then
!
!  Reciprocal space component for 2-D and 3-D cases  
!
    call kindex
    if (lsymopt.and.lsymderv) then
      if (lsymnotused) then
        if (nprocs.gt.1) then
          call recip3Dd(erecip,ec6,lgrad1,lgrad2)
        else
          call recip3D(erecip,ec6,lgrad1,lgrad2)
        endif
      elseif (lgrad2) then
        call recipsd2(erecip,ec6,lgrad1,lgrad2)
      else
        call recipsd(erecip,ec6,lgrad1)
      endif
    else
      if (ndim.eq.3) then
        if (lgrad2) then
          if (nprocs.gt.1) then
            call recip3Dd(erecip,ec6,lgrad1,lgrad2)
          else
            call recip3D(erecip,ec6,lgrad1,lgrad2)
          endif
        else
          if (lspmeloc) then
            call recip3Dspme(erecip,lgrad1)
          else
            if (luseDMalg) then
              call recip3Dd(erecip,ec6,lgrad1,lgrad2)
            else
              call recip3Dmd(erecip,ec6,lgrad1)
            endif
          endif
        endif
      elseif (ndim.eq.2) then
        if (luseDMalg) then
          call recip2Dd(erecip,esregion12,esregion2,eattach,lgrad1,lgrad2)
        else
          call recip2D(erecip,esregion12,esregion2,eattach,lgrad1,lgrad2)
        endif
      endif
    endif
  elseif (lewald.and.ndim.eq.0) then
    rmx2 = 1.0d10
  else
    rmx2 = 0.0_dp
  endif
!**********************
!  Real space energy  *
!**********************
  if (lewald.or.lwolf.or.npote.ne.0) then
    if (lsymopt) then
      if (lsymderv.and.(lgrad1.or.lgrad2)) then
        if (lsymnotused) then
          if (.not.lgrad2) then
            if (lminimage) then
              if (luseDMalg) then
                call realmi3d(eatom,ereal,erecip,ec6,eqeq,esregion12,esregion2,eattach,lgrad1)
              else
                call realmi3(eatom,ereal,erecip,ec6,eqeq,esregion12,esregion2,eattach,lgrad1)
              endif
            else
              if (luseDMalg) then
                call realed(eatom,ereal,erecip,ec6,eqeq,esregion12,esregion2,eattach,lgrad1,lgrad2)
              else
                call realmd3(eatom,ereal,erecip,ec6,eqeq,esregion12,esregion2,eattach,lgrad1)
              endif
            endif
          else
            if (nprocs.gt.1) then
              call realed(eatom,ereal,erecip,ec6,eqeq,esregion12,esregion2,eattach,lgrad1,lgrad2)
            else
              call reale(eatom,ereal,erecip,ec6,eqeq,esregion12,esregion2,eattach,lgrad1,lgrad2)
            endif
          endif
        elseif (lgrad2) then
          call realsd2(eatom,ereal,erecip,ec6,eqeq,lgrad1,lgrad2)
        else
          call realsd(eatom,ereal,erecip,ec6,eqeq,lgrad1)
        endif
      else
        if (lgrad1.or.lgrad2.or.lsymoffloc) then
          if (.not.lgrad2) then
            if (lminimage) then
              if (luseDMalg) then
                call realmi3d(eatom,ereal,erecip,ec6,eqeq,esregion12,esregion2,eattach,lgrad1)
              else
                call realmi3(eatom,ereal,erecip,ec6,eqeq,esregion12,esregion2,eattach,lgrad1)
              endif
            else
              if (luseDMalg) then
                call realed(eatom,ereal,erecip,ec6,eqeq,esregion12,esregion2,eattach,lgrad1,lgrad2)
              else
                call realmd3(eatom,ereal,erecip,ec6,eqeq,esregion12,esregion2,eattach,lgrad1)
              endif
            endif
          else
            if (luseDMalg) then
              call realed(eatom,ereal,erecip,ec6,eqeq,esregion12,esregion2,eattach,lgrad1,lgrad2)
            else
              call reale(eatom,ereal,erecip,ec6,eqeq,esregion12,esregion2,eattach,lgrad1,lgrad2)
            endif
          endif
        else
          call realsd(eatom,ereal,erecip,ec6,eqeq,lgrad1)
        endif
      endif
    else
      if (ndim.gt.0) then
        if (.not.lgrad2) then
          if (lminimage) then
            if (luseDMalg) then
              call realmi3d(eatom,ereal,erecip,ec6,eqeq,esregion12,esregion2,eattach,lgrad1)
            else
              call realmi3(eatom,ereal,erecip,ec6,eqeq,esregion12,esregion2,eattach,lgrad1)
            endif
          else
            if (lspatialok) then
              call realmd3s(eatom,ereal,erecip,ec6,eqeq,esregion12,esregion2,eattach,lgrad1)
            else
              if (luseDMalg) then
                call realed(eatom,ereal,erecip,ec6,eqeq,esregion12,esregion2,eattach,lgrad1,lgrad2)
              else
                call realmd3(eatom,ereal,erecip,ec6,eqeq,esregion12,esregion2,eattach,lgrad1)
              endif
            endif
          endif
        else
          if (nprocs.gt.1) then
            call realed(eatom,ereal,erecip,ec6,eqeq,esregion12,esregion2,eattach,lgrad1,lgrad2)
          else
            call reale(eatom,ereal,erecip,ec6,eqeq,esregion12,esregion2,eattach,lgrad1,lgrad2)
          endif
        endif
      else
        if (lcmm) then
          call realcmm(eatom,ereal,ecmm,eqeq,lgrad1)
        else
          if (.not.lgrad2) then
            if (lspatialok) then
              call realmd0s(eatom,ereal,eqeq,lgrad1)
            else
              if (luseDMalg) then
                call real0dd(eatom,ereal,eqeq,lgrad1,lgrad2)
              else
                call realmd0(eatom,ereal,eqeq,lgrad1)
              endif
            endif
          else
            if (nprocs.gt.1) then
              call real0dd(eatom,ereal,eqeq,lgrad1,lgrad2)
            else
              call real0d(eatom,ereal,eqeq,lgrad1,lgrad2)
            endif
          endif
        endif
      endif
    endif
    if (ndim.eq.1.and..not.lwolf) then
!
!  Real space Coulomb contribution from beyond potential cut-off in 1-D
!
      if (.not.lgrad2) then
        if (luseDMalg) then
          call real1Dd(ereal,esregion12,esregion2,lgrad1,lgrad2)
        else
          call real1Dmd(ereal,esregion12,esregion2,lgrad1)
        endif
      else
        if (nprocs.gt.1) then
          call real1Dd(ereal,esregion12,esregion2,lgrad1,lgrad2)
        else
          call real1D(ereal,esregion12,esregion2,lgrad1,lgrad2)
        endif
      endif
    endif
  endif
!*****************************
!  Point-ion polarisability  *
!*****************************
  if (lpolar) then
    call polarisation(epolar,esregion12,esregion2,eattach,lsymnotused,lgrad1,lgrad2)
  endif
!**********************************
!  Bond order charge self-energy  *
!**********************************
  if (nboQ0.gt.0) then
    if (luseDMalg) then
      call BOselfd(eboQself,lgrad1,lgrad2,.false.)
    else
      call BOself(eboQself,lgrad1,lgrad2,.false.)
    endif
  endif
!*********************
!  Solvation energy  *
!*********************
  if (lcosmo) then
    call solvation(ecosmo,lgrad1,lgrad2)
  endif
!**********************
!  Three-body energy  *
!**********************
  if (nthb.gt.0) then
    if (lsymderv2) then
      if (lgrad2) then
        call threesd2(ethb,lgrad1,lgrad2)
      else
        call threesd(ethb,lgrad1)
      endif
    else
      if (lgrad2) then
        if (nprocs.gt.1) then
          call threenosd(ethb,esregion12,esregion2,eattach,lgrad1,lgrad2)
        else
          call threenos(ethb,esregion12,esregion2,eattach,lgrad1,lgrad2)
        endif
      else
        if (lspatialok) then
          call threemds(ethb,esregion12,esregion2,eattach,lgrad1)
        else
          if (luseDMalg) then
            call threenosd(ethb,esregion12,esregion2,eattach,lgrad1,lgrad2)
          else
            call threemd(ethb,esregion12,esregion2,eattach,lgrad1)
          endif
        endif
      endif
    endif
  endif
!*********************
!  Four-body energy  *
!*********************
  if (nfor.gt.0) then
    if (lsymderv2) then
      if (lgrad2) then
        call foursd2(efor,eoop,eimp,lgrad1,lgrad2)
      else
        call foursd(efor,eoop,eimp,lgrad1)
      endif
    else
      if (lgrad2) then
        if (nprocs.gt.1) then
          call fournosd(efor,eoop,eimp,esregion12,esregion2,eattach,lgrad1,lgrad2)
        else
          call fournos(efor,eoop,eimp,esregion12,esregion2,eattach,lgrad1,lgrad2)
        endif
      else
        if (lspatialok) then
          call fourmds(efor,eoop,eimp,esregion12,esregion2,eattach,lgrad1)
        else
          call fourmd(efor,eoop,eimp,esregion12,esregion2,eattach,lgrad1)
        endif
      endif
    endif
  endif
!*********************
!  Six-body energy  *
!*********************
  if (nsix.gt.0) then
    if (lsymderv2) then
      if (lgrad2) then
        call sixsd2(esix,lgrad1,lgrad2)
      else
        call sixsd(esix,lgrad1)
      endif
    else
      if (lgrad2) then
        if (nprocs.gt.1) then
          call sixnosd(esix,esregion12,esregion2,eattach,lgrad1,lgrad2)
        else
          call sixnos(esix,esregion12,esregion2,eattach,lgrad1,lgrad2)
        endif
      else
        if (lspatialok) then
          call sixmds(esix,esregion12,esregion2,eattach,lgrad1)
        else
          call sixmd(esix,esregion12,esregion2,eattach,lgrad1)
        endif
      endif
    endif
  endif
!*********************
!  Many-body energy  *
!*********************
  if (lsuttonc) then
!---------------------------
!  Compute density : MEAM  !
!---------------------------
    if (lMEAMden) then
      if (lsymopt.and.(lsymderv.or..not.lgrad1)) then
        call densitysd(eatom)
      elseif (ndim.eq.0) then
        call density0(eatom)
      else
        if (lspatialok) then
          call density3s(eatom,esregion12,esregion2,eattach)
        else
          call density3(eatom,esregion12,esregion2,eattach)
        endif
      endif
    endif
!--------------------------------
!  Compute energy : EAM & MEAM  !
!--------------------------------
    if (lsymopt) then
      if (lsymderv.and.(lgrad1.or.lgrad2)) then
        if ((lgrad2.and..not.lsymderv2).or.lsymoffloc) then
          if (lgrad2) then
            if (nprocs.gt.1) then
              call manyd(emany,esregion12,esregion2,eattach,lgrad1,lgrad2)
            else
              call many(emany,esregion12,esregion2,eattach,lgrad1,lgrad2)
            endif
          else
            if (lminimage) then
              call manymi3(emany,esregion12,esregion2,eattach,lgrad1)
            else
              call manymd3(emany,esregion12,esregion2,eattach,lgrad1)
            endif
          endif
        else
          if (lgrad2) then
            call manysd2(emany,lgrad1,lgrad2)
          else
            call manysd(emany,lgrad1)
          endif
        endif
      else
        if (lgrad1.or.lgrad2.or.lsymoffloc) then
          if (lgrad2) then
            call many(emany,esregion12,esregion2,eattach,lgrad1,lgrad2)
          else
            if (lminimage) then
              call manymi3(emany,esregion12,esregion2,eattach,lgrad1)
            else
              if (lspatialok) then
                call manymd3s(emany,esregion12,esregion2,eattach,lgrad1)
              else
                call manymd3(emany,esregion12,esregion2,eattach,lgrad1)
              endif
            endif
          endif
        else
          if (lgrad2) then
            call manysd2(emany,lgrad1,lgrad2)
          else
            call manysd(emany,lgrad1)
          endif
        endif
      endif
    else
      if (ndim.gt.0) then
        if (lgrad2) then
          if (nprocs.gt.1) then
            call manyd(emany,esregion12,esregion2,eattach,lgrad1,lgrad2)
          else
            call many(emany,esregion12,esregion2,eattach,lgrad1,lgrad2)
          endif
        else
          if (lspatialok) then
            call manymd3s(emany,esregion12,esregion2,eattach,lgrad1)
          elseif (lminimage) then
            call manymi3(emany,esregion12,esregion2,eattach,lgrad1)
          else
            call manymd3(emany,esregion12,esregion2,eattach,lgrad1)
          endif
        endif
      else
        if (lgrad2) then
          if (nprocs.gt.1) then
            call many0dd(emany,lgrad1,lgrad2)
          else
            call many0d(emany,lgrad1,lgrad2)
          endif
        else
          call manymd0(emany,lgrad1)
        endif
      endif
    endif
  endif
!**********************
!  Brenner potential  *
!**********************
  if (lbrenner) then
    if (lsymopt.and.lsymderv) then
      if ((lgrad2.and..not.lsymderv2).or.lsymoffloc) then
        if (nprocs.gt.1) then
          call brennerd(ebrenner,0.0_dp,0.0_dp,0.0_dp,lgrad1,lgrad2,.false.)
        else
          call brenner(ebrenner,0.0_dp,0.0_dp,0.0_dp,lgrad1,lgrad2,.false.)
        endif
      else
        call brennersd2(ebrenner,lgrad1,lgrad2)
      endif
    else
      if (lgrad2) then
        if (nprocs.gt.1) then
          call brennerd(ebrenner,0.0_dp,0.0_dp,0.0_dp,lgrad1,lgrad2,.false.)
        else
          call brenner(ebrenner,0.0_dp,0.0_dp,0.0_dp,lgrad1,lgrad2,.false.)
        endif
      else
        call brennermd(ebrenner,lgrad1)
      endif
    endif
  endif
!**************************
!  Bond order potentials  *
!**************************
  if (nbopot.gt.0) then
    if (lsymopt.and.lsymderv) then
      if ((lgrad2.and..not.lsymderv2).or.lsymoffloc) then
        if (nprocs.gt.1) then
          call bondorderd(ebondorder,0.0_dp,0.0_dp,0.0_dp,lgrad1,lgrad2,.false.)
        else
          call bondorder(ebondorder,0.0_dp,0.0_dp,0.0_dp,lgrad1,lgrad2,.false.)
        endif
      else
        call bondordersd2(ebondorder,lgrad1,lgrad2)
      endif
    else
      if (lgrad2) then
        if (nprocs.gt.1) then
          call bondorderd(ebondorder,0.0_dp,0.0_dp,0.0_dp,lgrad1,lgrad2,.false.)
        else
          call bondorder(ebondorder,0.0_dp,0.0_dp,0.0_dp,lgrad1,lgrad2,.false.)
        endif
      else
        call bondordermd(ebondorder,lgrad1)
      endif
    endif
  endif
!*******************
!  OpenKIM models  *
!*******************
  if (lkim_model) then
    if (nkimmodel.gt.0) then
      call kimmd(ekim,esregion12,esregion2,eattach,lgrad1)
    endif
  endif
!**********************
!  ReaxFF forcefield  *
!**********************
  if (lreaxFF) then
!    if (lsymopt.and.lsymderv) then
!      if ((lgrad2.and..not.lsymderv2).or.lsymoffloc) then
!        call reaxFF(ereaxFF,0.0_dp,0.0_dp,0.0_dp,lgrad1,lgrad2,.false.)
!      else
!        call reaxFFsd2(ereaxFF,lgrad1,lgrad2)
!      endif
!    else
!      if (lgrad2) then
!        call reaxFF(ereaxFF,0.0_dp,0.0_dp,0.0_dp,lgrad1,lgrad2,.false.)
!      else
        call reaxFFmd(ereaxFF,lgrad1)
!      endif
!    endif
  endif
!********************
!  EDIP forcefield  *
!********************
  if (lEDIP) then
    call EDIPmd(eEDIP,lgrad1)
  endif
!*************************
!  Pressure-volume term  *
!*************************
  if (lpress.or.lanisotropicpress) call pressure(epv,lgrad1,lgrad2)
!*******************
!  External force  *
!*******************
  call force(eforce,lgrad1,lgrad2)
!*****************
!  Radial force  *
!*****************
  call radialforce(eradial,lgrad1,lgrad2)
!*******************
!  Electric field  *
!*******************
  if (lfieldcfg(ncf)) then
    call electricfield(efield,lgrad1)
  endif
!*******************
!  Einstein model  *
!*******************
  if (leinstein) then
    if (lsymopt.and.lsymderv) then
      if ((lgrad2.and..not.lsymderv2).or.lsymoffloc) then
        call einstein(eeinstein,lgrad1,lgrad2)
      else
        call einsteinsd2(eeinstein,lgrad1,lgrad2)
      endif
    else
      call einstein(eeinstein,lgrad1,lgrad2)
    endif
  endif
!************************
!  One-body potentials  *
!************************
  if (none.gt.0) then
    call onebody(eone,esregion2,eattach)
  endif
!*********************
!  Plane potentials  *
!*********************
  if ((ndim.eq.0.or.ndim.eq.2).and.nplanepot.gt.0) then
    if (lgrad2) then
      call planepot2D(eplane,esregion2,eattach,lgrad1,lgrad2)
    else
      call planepotmd(eplane,esregion2,eattach,lgrad1)
    endif
  endif
!***********************
!  Wolf sum self term  *
!***********************
  if (lwolf) then
    call wolfself(ewolfself,esregion2)
  endif
!****************************
!  Neutralising background  *
!****************************
  if (ndim.gt.0.and.abs(totalcharge).gt.1.0d-4) call background(ebgd,emad,lgrad1,lgrad2)
!*****************************
!  Dipole correction energy  *
!*****************************
  if ((ldipole.or.lddipole).and.(lewald.or.lwolf).and.ndim.gt.0) then
    if (lsymopt.and.lsymderv) then
      if ((lgrad2.and..not.lsymderv2).or.lsymoffloc) then
        if (nprocs.gt.1) then
          call dipole3Dd(edipole,lgrad1,lgrad2)
        else
          call dipole3D(edipole,lgrad1,lgrad2)
        endif
      else
        call dipolesd(edipole,lgrad1,lgrad2)
      endif
    else
      if (lgrad2) then
        if (nprocs.gt.1) then
          call dipole3Dd(edipole,lgrad1,lgrad2)
        else
          call dipole3D(edipole,lgrad1,lgrad2)
        endif
      else
        call dipole3D(edipole,lgrad1,lgrad2)
      endif
    endif
  endif
#ifdef OLDCS
!**********************************
!  Chemshell energy contribution  *
!**********************************
  if (ichemsh_qm.eq.1) then
    call chemshellenergy(echemsh)
  endif
#endif
!******************************************
!  Symmetrise second derivative matrices  *
!******************************************
  if (lgrad2.and..not.lsymderv2.and.nprocs.eq.1) call symderv2(lgrad2)
!********************************
!  Parallel sum of derivatives  *
!********************************
  if (nprocs.gt.1) then
    call psumall(eatom,ereal,erecip,ec6,eqeq,eattach,esregion12,esregion2,ethb,efor, &
      eoop,emany,ecmm,ebrenner,epolar,eeinstein,ewolfself,ebondorder,eforce,esix, &
      efield,eradial,ereaxFF,eplane,ecosmo,eone,eedip,eimp,ekim,eplumed,eboQself, &
      lgrad1,lgrad2,(lsymderv.and..not.lgrad2))
  endif
!********************************
!  Complete second derivatives  *
!********************************
  if (lgrad2) then
    if (lfreeze) then
      if (lsymderv2) then
        call sumderv2f(numat,nasym,lsymderv2)
      else
        call sumderv2f(numat,numat,lsymderv2)
      endif
    else
      if (lsymderv2) then
        call sumderv2s(numat,nasym,.false.,.true.)
      else
        if (nprocs.gt.1) then
          call sumderv2p(natomsonnode,node2atom,numat,.false.)
        else
          call sumderv2(natomsonnode,.false.)
        endif
      endif
    endif
  endif
!*****************************************
!  Construct rigid molecule derivatives  *
!*****************************************
  if (lrigid.and.lgrad1) then
    call rigidmoleculedrv(lgrad2)
  endif
!
!  Sum components of total energy
!
  sft = shift(nshcfg(ncf))*shscalecfg(ncf)
  etot = erecip + eatom + sft + ethb + ereal + efor + epv + ecmm + ec6 + eoop + emany + &
         edipole + ebgd + eself + eqeq + epolar + ebrenner + eforce + eimp + &
         esregion12 + esregion2 + eeinstein + ewolfself + ebondorder + eboQself + esix + &
         echargecoupled + efield + eradial + ereaxFF + eplane + ecosmo + eone + eedip + emad + &
         ekim + eplumed
!
!  Sum off diagonal region-region energies
!
  do i = 2,nregions(ncf)
    do j = 1,i-1
      esum = eregion2region(j,i) + eregion2region(i,j)
      eregion2region(j,i) = esum
      eregion2region(i,j) = esum
    enddo
  enddo
#ifdef OLDCS
!
!  Add Chemshell energy
!
  etot = etot + echemsh
#endif
!
  fcstore = etot
  ncflast = ncf
#ifdef TRACE
  call trace_out('energy')
#endif
!
  return
  end
