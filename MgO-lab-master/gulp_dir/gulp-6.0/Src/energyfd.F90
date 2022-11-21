  subroutine energyfd(matom,mcrd,step,etot,lgrad1)
!
!  Subroutine for calculating lattice energy - computes subset for matom
!
!  12/17 Created from energy
!   1/18 Trace added
!   5/18 Split bond EEM added
!   8/18 Modified for version 2 of OpenKIM
!   9/18 Change for multiple models in OpenKIM
!   9/18 Region and attachment energies added to kimmd arguments
!  11/18 Call to sumderv1 removed as no longer needed
!  12/18 lgrad2 added as an argument to force
!   8/19 Corrections to call of density3 to include all arguments
!   2/20 SPME turned off for rigid molecules as this is not implemented yet
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
  use m_bondorder,   only : bondorderfd
  use m_brenner,     only : brennerfd
  use m_edip,        only : edipfd
  use m_reaxff,      only : reaxFFfd
  use m_three
  use molecule
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
  integer(i4), intent(in)  :: matom              ! Atom whose coordinates are being changed
  integer(i4), intent(in)  :: mcrd               ! Coordinate number to change (1=x,2=y,3=z)
  real(dp),    intent(in)  :: step               ! Step size for coordinate
  real(dp),    intent(out) :: etot
  logical,     intent(in)  :: lgrad1
!
!  Local variables
!
  integer(i4)              :: i
  integer(i4),        save :: ncflast = 0
  logical                  :: lpress
  logical                  :: lspmeloc            ! If true, use of SPME is OK
  real(dp)                 :: sft
#ifdef TRACE
  call trace_in('energyfd')
#endif
!
  lpress = (abs(press).gt.0.0_dp.and.ndim.gt.0)
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
  if (lspatialok) then
    call setspatial(.false.)    ! Argument has to be false otherwise fails for phonon finite difference
  endif
  if (lspatialBOok) then
    call setspatialbo(.false.)  ! Argument has to be false otherwise fails for phonon finite difference
  endif
!
  call setatomnodes(numat,1_i4,0_i4,lspatialok)
  call setatomnodesbo(numat,1_i4,0_i4,lspatialBOok)
!*************************
!  SPME algorithm check  *
!*************************
  if (lDoQDeriv1.or.latomicstress.or.lsiteenergy.or.lrigid) lspmeloc = .false.
!**************************************
!  Call EEM/QEq to calculate charges  *
!**************************************
  if (leem) then
    if (leembond) then
      call eemsplit(.false.,lgrad1,.false.)
    else
      call eem(.false.,lgrad1,.false.)
    endif
  endif
!****************************************************
!  Calculate charges according to bond order model  *
!****************************************************
  if ((nboQ+nboQ0).gt.0) then
    call getBOcharge(lgrad1,.false.)
  endif
!*************************************************
!  Zero derivatives                              *
!  In ChemShell case the QM force is added here  *
!*************************************************
  call initdervs(lgrad1,.false.)
!
!  Redetermine cell indices for molecule atoms
!  incase one has moved across cell boundary.
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
    if (ndim.eq.3) then
      call recip3Dfd(matom,erecip,ec6,lgrad1)
    elseif (ndim.eq.2) then
      call recip2Dfd(matom,erecip,lgrad1)
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
    if (ndim.gt.0) then
      if (lminimage) then
        call realmi3fd(matom,eatom,ereal,erecip,ec6,eqeq,lgrad1)
      else
        call realfd3(matom,eatom,ereal,erecip,ec6,eqeq,lgrad1)
      endif
    else
      call realfd0(matom,eatom,ereal,eqeq,lgrad1)
    endif
    if (ndim.eq.1.and..not.lwolf) then
!
!  Real space Coulomb contribution from beyond potential cut-off in 1-D
!
      call real1Dfd(matom,ereal,lgrad1)
    endif
  endif
!**********************
!  Three-body energy  *
!**********************
  if (nthb.gt.0) then
    call threefd(matom,ethb,lgrad1)
  endif
!*********************
!  Four-body energy  *
!*********************
  if (nfor.gt.0) then
    call fourfd(matom,efor,eoop,eimp,lgrad1)
  endif
!*********************
!  Six-body energy  *
!*********************
  if (nsix.gt.0) then
    call sixfd(matom,esix,lgrad1)
  endif
!*********************
!  Many-body energy  *
!*********************
  if (lsuttonc) then
!---------------------------
!  Compute density : MEAM  !
!---------------------------
    if (lMEAMden) then
      if (ndim.eq.0) then
        call density0(eatom)
      else
        call density3(eatom,esregion12,esregion2,eattach)
      endif
    endif
!--------------------------------
!  Compute energy : EAM & MEAM  !
!--------------------------------
!    if (ndim.gt.0) then
!      if (lminimage) then
!        call manymi3fd(matom,emany,lgrad1)
!      else
!        call manyfd3(matom,emany,lgrad1)
!      endif
!    else
!      call manyfd0(matom,emany,lgrad1)
!    endif
  endif
!**********************
!  Brenner potential  *
!**********************
  if (lbrenner) then
    call brennerfd(matom,mcrd,step,ebrenner,lgrad1)
  endif
!**************************
!  Bond order potentials  *
!**************************
  if (nbopot.gt.0) then
    call bondorderfd(matom,mcrd,step,ebondorder,lgrad1)
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
    call reaxFFfd(matom,mcrd,step,ereaxFF,lgrad1)
  endif
!********************
!  EDIP forcefield  *
!********************
  if (lEDIP) then
    call EDIPfd(matom,mcrd,step,eEDIP,lgrad1)
  endif
!*************************
!  Pressure-volume term  *
!*************************
  if (lpress.or.lanisotropicpress) call pressure(epv,lgrad1,.false.)
!*******************
!  External force  *
!*******************
  call force(eforce,lgrad1,.false.)
!*****************
!  Radial force  *
!*****************
  call radialforce(eradial,lgrad1,.false.)
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
    call einstein(eeinstein,lgrad1,.false.)
  endif
!*********************
!  Plane potentials  *
!*********************
  if ((ndim.eq.0.or.ndim.eq.2).and.nplanepot.gt.0) then
    call planepotfd(matom,eplane,lgrad1)
  endif
!****************************
!  Neutralising background  *
!****************************
  if (ndim.gt.0.and.abs(totalcharge).gt.1.0d-4) call background(ebgd,emad,lgrad1,.false.)
!*****************************
!  Dipole correction energy  *
!*****************************
  if (ldipole.and.(lewald.or.lwolf).and.ndim.gt.0) then
    call dipole3D(edipole,lgrad1,.false.)
  endif
#ifdef OLDCS
!**********************************
!  Chemshell energy contribution  *
!**********************************
  if (ichemsh_qm.eq.1) then
    call chemshellenergy(echemsh)
  endif
#endif
!
!  Sum components of total energy
!
  sft = shift(nshcfg(ncf))*shscalecfg(ncf)
  etot = erecip + eatom + sft + ethb + ereal + efor + epv + ecmm + ec6 + eoop + emany + &
         edipole + ebgd + eself + eqeq + epolar + ebrenner + eforce + eimp + &
         eeinstein + ewolfself + ebondorder + eboQself + esix + &
         echargecoupled + efield + eradial + ereaxFF + eplane + ecosmo + eone + eedip + emad + &
         ekim + eplumed
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
  call trace_out('energyfd')
#endif
!
  return
  end
