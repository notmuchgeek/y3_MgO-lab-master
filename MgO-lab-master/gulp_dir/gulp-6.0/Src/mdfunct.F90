  subroutine mdfunct(istep,etot,lsetmol,lreset,lgrad1)
!
!  Subroutine for calculating energy and forces for MD.
!
!  Unlike general routines assumes no symmetry and works in
!  cartesian space for molecular dynamics.
!
!   6/95 Constraints added 
!   2/97 Dummy argument passed to energy routines in place
!        of derv2 as second derivatives are not needed
!   3/97 Modified for breathing shell MD
!   4/97 Sutton-Chen potential added
!   7/97 Neutralising background added
!  12/97 Possible second derivative matrix now passed as an
!        argument for use in variable charge calculation.
!   7/99 Minimum image option added for speed for big systems
!  12/00 Modifications for 1-D/2-D systems made
!   4/01 Point ion polarisability added
!   9/01 Modifications for 1-D made
!   5/02 Scaling of shift added
!   5/02 Brenner potential added
!  11/02 External forces and Einstein model added
!  11/02 Parallel sum of derivatives call added
!  12/02 Conversion of Cartesian to fraction in 2-D corrected
!   1/03 Wolf sum modifications made
!   2/03 eattach now passed correctly to subroutines
!   5/03 Spatial decomposition option added
!   6/03 Error in argument list to setatomnodes fixed
!   9/03 Spatial algorithms added for real space and many-body case
!  11/03 Bond order potentials added
!   4/04 Logical controlling initialisation of strain derivatives 
!        changed to lstr from lconp
!   9/04 eforce passed to psumall
!   1/05 Option to use tabulate distance vectors added
!   1/05 Logical as to whether interatomic vector table is to be reset
!        added as an input argument
!   5/06 Specific real1Dmd routine called
!   7/06 Sixbody potentials added
!   8/06 Vibrational energy initialised though not used
!   3/07 Electric field option added
!   3/07 Radial force added
!   3/07 Chemshell energy added
!   3/07 Gauss renamed to GULP_gauss
!   5/07 Call to setspatial(bo) modified
!   7/07 Metadynamics modifications added
!   7/07 ereaxFF added
!   7/07 eplane added
!  11/07 Modifications for reaxFF Coulomb term added
!  12/07 Unused variables removed
!  12/07 Call to reaxffmd added
!   1/08 lreaxFFqreal removed
!   3/08 Modified for new MD integrator
!   4/08 Call to spatial decomposition version of setmol added
!  10/08 Error in conversion of 1-D coordinates to fractional corrected
!  12/08 scrho, and defect versions, converted to 2-D array for MEAM
!   4/09 Separate call to generate MEAM densities added
!   6/09 Site energy initialisation added
!   6/09 Module name changed from three to m_three
!   7/09 cutoffmax(bo) now passed via general module
!   7/09 Call to setcutoffmax(bo) removed
!  11/09 Region forces added
!   1/10 One-body potentials added
!   8/10 Spatial decomposition now allowed for ReaxFF
!   9/10 EDIP energy added
!  10/10 EDIP enabled for spatial decomposition
!  11/10 Anisotropic pressure added
!   8/11 Step number passed in for use with plumed
!   8/11 Subroutine name changed from .f90 to .F90 to allow for plumed
!   9/11 lgrad1 argument added to eem call
!   9/11 Metadynamics internal code replaced with Plumed
!   9/11 Madelung correction added
!  11/11 Initialisation of eregion2region added
!  11/11 Summing of off-diagonal region-region energies added
!   4/12 Virial now set using diagonal elements of strain derivative tensor
!   4/12 xvir, yvir and zvir removed
!   5/12 Atomic stress initialisation added
!  10/12 Support for OpenKIM models added
!  12/12 Time-dependent field added
!  12/12 Modified to allow for multiple time-dependent fields
!   7/13 Improper energy added
!   9/13 Call to recip3Dmd added back with new version of code
!   9/13 Extra timing added
!   9/13 Use of GULP_gauss replaced by matinv/blas
!   1/14 Initialisation of rstrdnr added
!   1/14 Non-radial strain derivatives added 
!   3/14 Calls to matinv renamed to matrix_inversion for benefit of ChemShell
!   4/14 Call to wolfself corrected by adding esregion2 to argument list
!   8/14 eatom added to call to density routines
!   2/15 Criterion for calling getBOcharge changed by adding nboQ0
!   3/15 Initialisation of xdrvnr/ydrvnr/zdrvnr added
!   3/15 Adding in of non-radial forces included
!   9/15 Spatial option added for 0-D case
!   3/16 Condition added for calling real space routines
!   3/16 Plumed modified for version 2
!   3/16 eplumed added
!   8/16 SPME added 
!   8/16 Use of lconp replaced with ncell
!   9/16 lgrad2 now passed to psumall
!   4/17 ChemShell interaction modified
!   6/17 Old ChemShell restored as a compile option
!   7/17 eboQself added to call to psumall
!   7/17 Use of distributed memory charge derivatives disabled for spatial
!        and cell multipole moment algorithms
!  10/17 Delta_dipole added
!   1/18 Trace added
!   5/18 Split bond EEM added
!   8/18 Modified for version 2 of OpenKIM
!   9/18 Change for multiple models in OpenKIM
!   9/18 Region and attachment energies added to kimmd arguments
!  11/18 Non-radial arrays removed since they are no longer needed
!  12/18 lgrad2 added as an argument to force
!   3/19 Constraint arrays changed to have index and type
!   4/19 Second derivative argument added to polarisation call
!   7/19 Extra argument added to polarisation call for symmetry
!   8/19 Argument list to density3 corrected
!   8/19 molatom option added
!   4/20 Handling of lmolrigid added
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
  use bondorderdata,  only : nbopot, nboQ, nboQ0
  use cellmultipole
  use configurations, only : maxregion, nregions
  use control
  use current
  use derivatives
  use distances,      only : lStoreVectors
  use eam,            only : lMEAM, lMEAMden, maxmeamcomponent
  use energies
  use field,          only : lfieldcfg, ntdfieldcfg
  use four
  use kim_models,     only : lkim_model, nkimmodel
  use kspace
  use mdlogic
  use moldyn
  use molecule
#ifdef PLUMED
  use m_plumed,       only : lplumed, plumed_energy_force
#endif
  use m_pr,           only : virial_m
  use m_three
  use one,            only : none
  use optimisation
  use parallel,       only : nprocs, procid
  use plane,          only : nplanepot
  use polarise
  use potentialxyz
  use shifts
  use six
  use spatial,        only : lspatialok
  use spatialbo,      only : lspatialBOok
  use spme,           only : lspme
  use sutton
  use symmetry
  use times,          only : tmdinit
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  use two
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)  :: istep
  logical,     intent(in)  :: lgrad1
  logical,     intent(in)  :: lreset
  logical,     intent(in)  :: lsetmol
  real(dp),    intent(out) :: etot
!
!  Local variables
!
  integer(i4)              :: i
  integer(i4)              :: icos1
  integer(i4)              :: icos2
  integer(i4)              :: icos3
  integer(i4)              :: ifail
  integer(i4)              :: j
  integer(i4)              :: nf
  integer(i4)              :: nv
  integer(i4),        save :: ncflast = 0
  logical                  :: lcmm
  logical                  :: lreindex
  logical                  :: luseDMalg           ! If true, use distributed memory algorithm
  real(dp)                 :: cos1
  real(dp)                 :: g_cpu_time
  real(dp)                 :: esum
  real(dp)                 :: ri
  real(dp)                 :: rmat(3,3)
  real(dp)                 :: sft
  real(dp)                 :: sum
  real(dp)                 :: t1
  real(dp)                 :: t2
  real(dp)                 :: vcrd
  real(dp)                 :: v1(3)
  real(dp)                 :: v2(3)
  real(dp)                 :: wrk(6)
  real(dp)                 :: xal
  real(dp)                 :: yal
  real(dp)                 :: zal
#ifdef TRACE
  call trace_in('mdfunct')
#endif
!
!  Initialise logicals
!
  lreindex = ((nmol.gt.0.or.(llist3.or.llist4.or.llist6)).and.ndim.gt.0)
  if ((llist3.and.llist4.or.llist6).and.(nthb.gt.0.or.nfor.gt.0).and.ndim.gt.0) lreindex = .true.
  lcmm = (icmm.gt.0.and.ndim.eq.0)
!
!  Zero energy contributions
!
  etot = 0.0_dp
  evib = 0.0_dp
  ereal = 0.0_dp
  erecip = 0.0_dp
  eatom = 0.0_dp
  eedip = 0.0_dp
  ethb = 0.0_dp
  eeinstein = 0.0_dp
  efield = 0.0_dp
  efor = 0.0_dp
  eforce = 0.0_dp
  eimp = 0.0_dp
  ekim = 0.0_dp
  emany = 0.0_dp
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
  ewolfself = 0.0_dp
  esregion12 = 0.0_dp
  esregion2 = 0.0_dp
  eattach = 0.0_dp
  ebondorder = 0.0_dp
  eboQself = 0.0_dp
  eradial = 0.0_dp
  ereaxFF = 0.0_dp
  eplane = 0.0_dp
  ecosmo = 0.0_dp
  eplumed = 0.0_dp
  echemsh = 0.0_dp
  virial = 0.0_dp
!
  t1 = g_cpu_time()
!
!  Initialise region-region interaction energies
!
  eregion2region(1:nregions(ncf),1:nregions(ncf)) = 0.0_dp
!
!  Site energy initialisation
!
  siteenergy(1:numat) = 0.0_dp
!
!  Zero many-body terms
!
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
!**********************
!  Apply constraints  *
!**********************
  if (ncon.gt.0) then
    do j = 1,ncon
      nf = ncfixind(j)
      if (ncfixtyp(j).eq.iopt_xf) then
        xalat(nf) = 0.0_dp
      elseif (ncfixtyp(j).eq.iopt_yf) then
        yalat(nf) = 0.0_dp
      elseif (ncfixtyp(j).eq.iopt_zf) then
        zalat(nf) = 0.0_dp
      endif
    enddo
    do j = 1,ncon
      nv = ncvarind(j)
      if (ncvartyp(j).eq.iopt_xf) then
        vcrd = xalat(nv)
      elseif (ncvartyp(j).eq.iopt_yf) then
        vcrd = yalat(nv)
      elseif (ncvartyp(j).eq.iopt_zf) then
        vcrd = zalat(nv)
      endif
      nf = ncfixind(j)
      if (ncfixtyp(j).eq.iopt_xf) then
        xalat(nf) = vcrd*conco(j) + conadd(j) + xalat(nf)
      elseif (ncfixtyp(j).eq.iopt_yf) then
        yalat(nf) = vcrd*conco(j) + conadd(j) + yalat(nf)
      elseif (ncfixtyp(j).eq.iopt_zf) then
        zalat(nf) = vcrd*conco(j) + conadd(j) + zalat(nf)
      endif
    enddo
  endif
  if (ndim.eq.3) then
!
!  If conp change cell parameters
!
    if (ncell.gt.0) then
      rv(1,1) = xcell(1)
      rv(2,1) = xcell(2)
      rv(3,1) = xcell(3)
      rv(1,2) = xcell(4)
      rv(2,2) = xcell(5)
      rv(3,2) = xcell(6)
      rv(1,3) = xcell(7)
      rv(2,3) = xcell(8)
      rv(3,3) = xcell(9)
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
      call rlist
    endif
!
!  Invert cell matrix
!
    rmat(1,1) = r1x
    rmat(2,1) = r1y
    rmat(3,1) = r1z
    rmat(1,2) = r2x
    rmat(2,2) = r2y
    rmat(3,2) = r2z
    rmat(1,3) = r3x
    rmat(2,3) = r3y
    rmat(3,3) = r3z
    call matrix_inversion(rmat,3_i4,3_i4,wrk,ifail)
!
!  Place coordinates back in main unit cell
!
    do i = 1,numat
      if (lopf(i)) then
        xal = xalat(i)
        yal = yalat(i)
        zal = zalat(i)
        ri = xal*xal + yal*yal + zal*zal
        if (abs(ri).lt.1.0d-12) then
          icos1 = 0
          icos2 = 0
          icos3 = 0
          xclat(i) = xal
          yclat(i) = yal
          zclat(i) = zal
        else
          v1(1) = xal
          v1(2) = yal
          v1(3) = zal
          call dgemv('N',3_i4,3_i4,1.0_dp,rmat,3_i4,v1,1_i4,0.0_dp,v2,1_i4)
!
!  If conp run use fractional coordinates to correct for change in cell parameters
!
          icos1 = int(v2(1)+100.0_dp) - 100
          icos2 = int(v2(2)+100.0_dp) - 100
          icos3 = int(v2(3)+100.0_dp) - 100
          xclat(i) = xal - icos1*r1x - icos2*r2x - icos3*r3x
          yclat(i) = yal - icos1*r1y - icos2*r2y - icos3*r3y
          zclat(i) = zal - icos1*r1z - icos2*r2z - icos3*r3z
        endif
        if (lreindex) call mdcellind(i,icos1,icos2,icos3)
      endif
    enddo
  elseif (ndim.eq.2) then
!
!  If conp change cell parameters
!
    if (ncell.gt.0) then
      rv(1,1) = xcell(1)
      rv(2,1) = xcell(2)
      rv(1,2) = xcell(3)
      rv(2,2) = xcell(4)
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
      r2x = rv(1,2)
      r2y = rv(2,2)
      sum = abs(r2x) + abs(r1y)
      lra = (sum.lt.1.0d-6)
      call rlist
    endif
!
!  Invert cell matrix
!
    rmat(1,1) = r1x
    rmat(2,1) = r1y
    rmat(1,2) = r2x
    rmat(2,2) = r2y
    call matrix_inversion(rmat,3_i4,2_i4,wrk,ifail)
!
!  Place coordinates back in main unit cell
!
    do i = 1,numat
      if (lopf(i)) then
        xal = xalat(i)
        yal = yalat(i)
        zal = zalat(i)
        ri = xal*xal + yal*yal + zal*zal
        if (abs(ri).lt.1.0d-12) then
          icos1 = 0
          icos2 = 0
          icos3 = 0
          xclat(i) = xal
          yclat(i) = yal
          zclat(i) = zal
        else
          v1(1) = xal
          v1(2) = yal
          call dgemv('N',2_i4,2_i4,1.0_dp,rmat,3_i4,v1,1_i4,0.0_dp,v2,1_i4)
!
!  If conp run use fractional coordinates to correct for change in cell parameters
!
          icos1 = int(v2(1)+100.0_dp) - 100
          icos2 = int(v2(2)+100.0_dp) - 100
          icos3 = 0
          xclat(i) = xal - icos1*r1x - icos2*r2x
          yclat(i) = yal - icos1*r1y - icos2*r2y
          zclat(i) = zal
        endif
        if (lreindex) call mdcellind(i,icos1,icos2,icos3)
      endif
    enddo
  elseif (ndim.eq.1) then
!
!  If conp change cell parameters
!
    if (ncell.gt.0) then
      rv(1,1) = xcell(1)
      call uncell1D(rv,a)
      if (a.gt.1.0d-12) then
        recipa = 1.0_dp/a
      else
        recipa = 0.0_dp
      endif
      r1x = rv(1,1)
      lra = .true.
      call rlist
    endif
!
!  Place coordinates back in main unit cell
!
    do i = 1,numat
      if (lopf(i)) then
        xal = xalat(i)
        yal = yalat(i)
        zal = zalat(i)
        ri = xal*xal + yal*yal + zal*zal
        if (abs(ri).lt.1.0d-12) then
          icos1 = 0
          icos2 = 0
          icos3 = 0
          xclat(i) = xal
          yclat(i) = yal
          zclat(i) = zal
        else
          cos1 = xal/r1x
!
!  If conp run use fractional coordinates to correct for change
!  in cell parameters
!
          icos1 = int(cos1+100.0_dp) - 100
          icos2 = 0
          icos3 = 0
          xclat(i) = xal - icos1*r1x
          yclat(i) = yal
          zclat(i) = zal
        endif
        if (lreindex) call mdcellind(i,icos1,icos2,icos3)
      endif
    enddo
  else
    do i = 1,numat
      xclat(i) = xalat(i)
      yclat(i) = yalat(i)
      zclat(i) = zalat(i)
    enddo
  endif
!*******************************************
!  Spatial decomposition set up if needed  *
!*******************************************
  lspatialok = lspatial
  lspatialBOok = (lspatial.and.(lbrenner.or.lreaxFF.or.lEDIP.or.(nbopot+nboQ).gt.0))
  if (lspatialok) then
    call setspatial(.false.)
  endif
  if (lspatialBOok) then
    call setspatialbo(.false.)
  endif
!  
!  Calculate parallel division of work :
!        
!  Spatial - divide cells over processors
!  Non-spatial - divide atoms over processors
!   
  call setatomnodes(numat,nprocs,procid,lspatialok)
  call setatomnodesbo(numat,nprocs,procid,lspatialBOok)
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
  luseDMalg = (nprocs.gt.1.and.(lDoQDeriv1.and..not.lspatialok.and..not.lcmm))
!**********************************
!  Set up of interatomic vectors  *
!**********************************
  if (lStoreVectors) then
    call setdistance(lreset)
  endif
!*******************************************
!  Set up molecules if needed (GCMD only)  *
!*******************************************
  if (lsetmol) then
    if (lmolatom.or.lmolrigid) then
      call setmola
    elseif (lspatialOK) then
      call setmols
    else
      call setmol
    endif
    if (nthb.gt.0.and.llist3) call setlist3
    if (nfor.gt.0.and.llist4) call setlist4
    if (nsix.gt.0.and.llist6) call setlist6
  endif
!***************************
!  Zero first derivatives  *
!***************************
  if (lgrad1) then
    if (nbsmat.gt.0) then
      do i = 1,numat
        raderv(i) = 0.0_dp
      enddo
    endif
    do i = 1,numat
      xdrv(i) = 0.0_dp
      ydrv(i) = 0.0_dp
      zdrv(i) = 0.0_dp
    enddo
    xregdrv(1:maxregion) = 0.0_dp
    yregdrv(1:maxregion) = 0.0_dp
    zregdrv(1:maxregion) = 0.0_dp
    if (lstr) then
      do i = 1,nstrains
        strderv(i) = 0.0_dp
        rstrd(i) = 0.0_dp
      enddo
      if (latomicstress) then
        do i = 1,numat
          atomicstress(1:nstrains,i) = 0.0_dp
        enddo
      endif
    endif
  endif
!**************************************
!  Initialise electric field to zero  *
!**************************************
  if (lpolar) then
    do i = 1,numat
      vx(i) = 0.0_dp
      vy(i) = 0.0_dp
      vz(i) = 0.0_dp
      vx12(i) = 0.0_dp
      vy12(i) = 0.0_dp
      vz12(i) = 0.0_dp
    enddo
    if (lqpolar) then
      do i = 1,nasym
        v2xyz(1:6,i) = 0.0_dp
        v2xyz12(1:6,i) = 0.0_dp
      enddo
    endif
  endif
!
  t2 = g_cpu_time()
  tmdinit = tmdinit + t2 - t1
!**************************************
!  Call EEM/QEq to calculate charges  *
!**************************************
  if (leem) then
    if (luseDMalg) then
      if (leembond) then
        call eemsplitd(.false.,lgrad1,.false.)
      else
        call eemdm(.false.,lgrad1,.false.)
      endif
    else
      if (leembond) then
        call eemsplit(.false.,lgrad1,.false.)
      else
        call eem(.false.,lgrad1,.false.)
      endif
    endif
  endif
!****************************************************
!  Calculate charges according to bond order model  *
!****************************************************
  if ((nboQ+nboQ0).gt.0) then
    if (luseDMalg) then
      call getBOcharged(lgrad1,.false.)
    else
      call getBOcharge(lgrad1,.false.)
    endif
  endif
!***************
!  Set up CMM  *
!***************
  if (lcmm) call setcmm
!*******************************
!  Electrostatic contribution  *
!*******************************
  if (lewald.and.ndim.gt.1) then
    if (ncf.ne.ncflast.or.ncell.gt.0) call kindex
    if (ndim.eq.3) then
      if (lspme) then
        call recip3Dspme(erecip,lgrad1)
      else
        if (luseDMalg) then
          call recip3Dd(erecip,ec6,lgrad1,.false.)
        else
          call recip3Dmd(erecip,ec6,lgrad1)
        endif
      endif
    elseif (ndim.eq.2) then
      if (luseDMalg) then
        call recip2Dd(erecip,esregion12,esregion2,eattach,lgrad1,.false.)
      else
        call recip2D(erecip,esregion12,esregion2,eattach,lgrad1,.false.)
      endif
    endif
  elseif (lewald.and.ndim.eq.0) then
    rmx2 = 1.0d10
  else
    rmx2 = 0.0_dp
  endif
!*************************
!  Real space component  *
!*************************
  if (lewald.or.lwolf.or.npote.ne.0) then
    if (ndim.gt.0) then
      if (lspatialok) then
        call realmd3s(eatom,ereal,erecip,ec6,eqeq,esregion12,esregion2,eattach,lgrad1)
      elseif (lminimage) then
        if (luseDMalg) then
          call realmi3d(eatom,ereal,erecip,ec6,eqeq,esregion12,esregion2,eattach,lgrad1)
        else
          call realmi3(eatom,ereal,erecip,ec6,eqeq,esregion12,esregion2,eattach,lgrad1)
        endif
      elseif (lStoreVectors) then
        call realmd3t(eatom,ereal,erecip,ec6,eqeq,esregion12,esregion2,eattach,lgrad1)
      else
        if (luseDMalg) then
          call realed(eatom,ereal,erecip,ec6,eqeq,esregion12,esregion2,eattach,lgrad1,.false.)
        else
          call realmd3(eatom,ereal,erecip,ec6,eqeq,esregion12,esregion2,eattach,lgrad1)
        endif
      endif
      if (ndim.eq.1) then
!  
!  Real space Coulomb contribution from beyond potential cut-off in 1-D
!  
        if (luseDMalg) then
          call real1Dd(ereal,esregion12,esregion2,lgrad1,.false.)
        else
          call real1Dmd(ereal,esregion12,esregion2,lgrad1)
        endif
      endif
    else
      if (lspatialok) then
        call realmd0s(eatom,ereal,eqeq,lgrad1)
      elseif (lcmm) then
        call realcmm(eatom,ereal,ecmm,eqeq,lgrad1)
      else
        if (luseDMalg) then
          call real0dd(eatom,ereal,eqeq,lgrad1,.false.)
        else
          call realmd0(eatom,ereal,eqeq,lgrad1)
        endif
      endif
    endif
  endif
!*****************************
!  Point ion polarisability  *
!*****************************
  if (lpolar) then
    call polarisation(epolar,esregion12,esregion2,eattach,.true.,lgrad1,.false.)
  endif
!**********************************
!  Bond order charge self-energy  *
!**********************************
  if (nboQ0.gt.0) then
    if (luseDMalg) then
      call BOselfd(eboQself,lgrad1,.false.,.false.)
    else
      call BOself(eboQself,lgrad1,.false.,.false.)
    endif
  endif
!*************************
!  Three-body component  *
!*************************
  if (nthb.gt.0) then
    if (llist3) then
      call threelist(ethb,lgrad1)
    else
      if (lspatialok) then
        call threemds(ethb,esregion12,esregion2,eattach,lgrad1)
      else
        call threemd(ethb,esregion12,esregion2,eattach,lgrad1)
      endif
    endif
  endif
!************************
!  Four-body component  *
!************************
  if (nfor.gt.0) then
    if (llist4) then
      call fourlist(efor,eoop,eimp,lgrad1)
    else
      if (lspatialok) then
        call fourmds(efor,eoop,eimp,esregion12,esregion2,eattach,lgrad1)
      else
        call fourmd(efor,eoop,eimp,esregion12,esregion2,eattach,lgrad1)
      endif
    endif
  endif
!***********************
!  Six-body component  *
!***********************
  if (nsix.gt.0) then
    if (llist6) then
      call sixlist(esix,lgrad1)
    else
      if (lspatialok) then
        call sixmds(esix,esregion12,esregion2,eattach,lgrad1)
      else
        call sixmd(esix,esregion12,esregion2,eattach,lgrad1)
      endif
    endif
  endif
!************************
!  Many-body component  *
!************************
  if (lsuttonc) then
!---------------------------
!  Compute density : MEAM  !
!---------------------------
    if (lMEAMden) then
      if (ndim.gt.0) then
        if (lspatialok) then
          call density3s(eatom,esregion12,esregion2,eattach)
        else
          call density3(eatom,esregion12,esregion2,eattach)
        endif
      else
        call density0(eatom)
      endif
    endif
!--------------------------------
!  Compute energy : EAM & MEAM  !
!--------------------------------
    if (ndim.gt.0) then
      if (lspatialok) then
        call manymd3s(emany,esregion12,esregion2,eattach,lgrad1)
      else
        call manymd3(emany,esregion12,esregion2,eattach,lgrad1)
      endif
    else
      call manymd0(emany,lgrad1)
    endif
  endif
!***********************
!  Brenner potentials  *
!***********************
  if (lbrenner) then
    call brennermd(ebrenner,lgrad1)
  endif
!**************************
!  Bond order potentials  *
!**************************
  if (nbopot.gt.0) then
    call bondordermd(ebondorder,lgrad1)
  endif
!*******************
!  OpenKIM models  *
!*******************
  if (lkim_model) then
    if (nkimmodel.gt.0) then
      call kimmd(ekim,esregion12,esregion2,eattach,lgrad1)
    endif
  endif
!***********************
!  ReaxFF force field  *
!***********************
  if (lreaxFF) then
    call reaxFFmd(ereaxFF,lgrad1)
  endif
!*********************
!  EDIP force field  *
!*********************
  if (lEDIP) then
    call EDIPmd(eEDIP,lgrad1)
  endif
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
  if (lfieldcfg(ncf).or.(ntdfieldcfg(ncf).gt.0)) then
    call electricfield(efield,lgrad1)
  endif
!*******************
!  Einstein model  *
!*******************
  if (leinstein) call einstein(eeinstein,lgrad1,.false.)
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
    call planepotmd(eplane,esregion2,eattach,.true.)
  endif
!***********************
!  Wolf sum self term  *
!***********************
  if (lwolf) call wolfself(ewolfself,esregion2)
!****************************
!  Neutralising background  *
!****************************
  if (ndim.gt.0.and.abs(totalcharge).gt.1.0d-6) call background(ebgd,emad,lgrad1,.false.)
!*****************************
!  Dipole correction energy  *
!*****************************
  if ((ldipole.or.lddipole).and.(lewald.or.lwolf).and.ndim.gt.0) then
    call dipole3D(edipole,lgrad1,.false.)
  endif
!*****************************
!  Plumed correction energy  *
!*****************************
#ifdef PLUMED
  if (lplumed.and.istep.ne.0) then
!
! PluMeD modifications
!
    call plumed_energy_force(istep,eplumed,xdrv,ydrv,zdrv)
  endif
#endif
!********************************
!  Parallel sum of derivatives  *
!********************************
  if (nprocs.gt.1) then
    call psumall(eatom,ereal,erecip,ec6,eqeq,eattach,esregion12,esregion2,ethb,efor, &
                 eoop,emany,ecmm,ebrenner,epolar,eeinstein,ewolfself,ebondorder, &
                 eforce,esix,efield,eradial,ereaxFF,eplane,ecosmo,eone,eedip,eimp,ekim, &
                 eplumed,eboQself,lgrad1,.false.,.false.)
  endif
!
!  Complete derivatives by adding non-radial terms 
!
  if (lgrad1) then
!
!  Sum components of strain if required
!
    if (lstr) then
      do i = 1,nstrains
        strderv(i) = strderv(i) + rstrd(i)
      enddo
      select case(ndim)
        case(1)
          virial = strderv(1)
        case(2)
          virial = strderv(1) + strderv(2)
        case(3)
          virial = strderv(1) + strderv(2) + strderv(3)
      end select
      if (nmdintegrator.eq.4) then
        virial_m(1,1) = strderv(1)
        virial_m(2,1) = strderv(6)
        virial_m(3,1) = strderv(5)
        virial_m(1,2) = strderv(6)
        virial_m(2,2) = strderv(2)
        virial_m(3,2) = strderv(4)
        virial_m(1,3) = strderv(5)
        virial_m(2,3) = strderv(4)
        virial_m(3,3) = strderv(3)
      endif
    endif
  endif
!
!  Sum components of total energy
!
  sft = shift(nshcfg(ncf))*shscalecfg(ncf)
  epv = 0.0_dp
  etot = erecip + eatom + sft + ethb + ereal + efor + epv + ecmm + ec6 + eoop + emany + &
         edipole + ebgd + eself + eqeq + epolar + ebrenner + eforce + eeinstein + ewolfself + &
         ebondorder + eboQself + esix + efield + eradial + ereaxFF + eplane + ecosmo + &
         eone + eedip + emad + eimp + ekim + eplumed
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
!
  if (leprint) call outener
  ncflast = ncf
#ifdef TRACE
  call trace_out('mdfunct')
#endif
!
  return
  end
