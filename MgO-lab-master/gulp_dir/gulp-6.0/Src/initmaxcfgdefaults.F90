  subroutine initmaxcfgdefaults(icfg)
!
!  Reinitialises defaults associated with configuration icfg
!
!   9/10 Created from changemaxcfg
!  11/10 Anisotropic pressure added
!   9/11 Metadynamics internal code replaced with Plumed
!   6/12 nobsmodecfg added
!  10/12 OpenKIM arrays added
!  12/12 Time-dependent field added
!  12/12 Modified to allow for multiple time-dependent fields
!   3/13 nreg1cfg, nvacacfg and nintecfg added
!   3/13 lindvacptr added
!   7/13 Symmetry number added
!   8/13 Thermal conductivity arrays added
!   9/13 nsuper changed to be a 2-D array
!   9/13 Raman directions added
!  10/14 nd2cellcfg added
!   4/15 Ghost supercell array added
!   8/15 New lower options added
!   9/15 translate noise added
!   5/16 SPME added
!   8/16 OpenKIM changes merged in
!   8/16 Default for nqkgrid increased
!   3/17 fix_atom option added
!   4/17 qonsas option added
!  10/17 Initial coordinates added to restart info
!   1/18 Initialisation of pkim_model removed for F03 version
!   2/18 ltrantherm added
!   2/18 nxks, nyks, nzks converted to a single array
!   2/18 Arrays for handling primitive cell added
!   3/18 Sign option added to translate
!   4/18 Arrays for scan_cell added
!   4/18 Twist option added
!   5/18 lcscanstrain added
!   6/18 straincfg added
!   8/18 Modified for version 2 of OpenKIM
!   9/18 Change for multiple models in OpenKIM
!   9/18 Correction to initialisation of lkim_model_cfg_OK
!  12/18 Shear force added
!   3/19 Multiple temperature ramps added
!   8/19 molatom option added
!   9/19 lvecpin added
!   2/20 npotsitescfg added
!   2/20 lksorigin added
!   3/20 lsumcoreshell added
!   3/20 Dielectric constant added
!   4/20 Restart arrays added for rigid molecules
!   6/20 nmolcore arrays added
!   8/20 ordersuper added
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
!  Julian Gale, CIC, Curtin University, August 2020
!
  use configurations
  use cosmic
  use current,        only : maxat
  use defects
  use distances,      only : ndistancereset
  use field
  use fitting,        only : lsumcoreshell
  use freeze
  use genetic,        only : xmaxcfg, xmincfg
  use kim_models,     only : lkim_model_cfg_OK, maxkimmodel
  use ksample
  use moldyn
  use molecule,       only : maxmol, nmolcfg, nmolatomcfg, nmolatomtotcfg
  use molecule,       only : nmolcorecfg
  use molecule,       only : molQcfg, molQxyzcfg, molcomcfg
  use m_pr
  use g_neb
  use observables,    only : freaction, maxobs, nobsmodecfg
  use potentialgrid
  use potentialpoints
  use potentialsites
  use projectdos
  use radial
  use reallocate
  use scan
  use shifts
  use spme,           only : nqkgrid
  use symmetry
  use thermalcond
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in) :: icfg
!
!  Local variables
!
  integer(i4)             :: j
!
!  Initialise defaults for new part of array
!
  if (icfg.le.maxcfg.and.icfg.ge.1) then
    anisotropicpresscfg(1:6,icfg) = 0.0_dp
    bornk(1,icfg) = 1.0_dp
    bornk(2,icfg) = 1.0_dp
    bornk(3,icfg) = 1.0_dp
    cosmoepsilon(icfg) = 1.0_dp
    cosmodrsolv(icfg) = 0.1_dp
    cosmorsolv(icfg) = 1.0_dp
    dhklcfg(icfg) = 0.0_dp
    dielectriccfg(icfg) = 1.0_dp
    energycfg(icfg) = 0.0_dp
    freaction(icfg,1:maxobs) = 0.0_dp
    lanisotropicpresscfg(icfg) = .false.
    lcosmoeigin(icfg) = .false.
    lcscanstrain(icfg) = .false.
    ldeflin(icfg) = .false.
    lfcborn(icfg) = .false.
    lfcphon(icfg) = .false.
    lfcprop(icfg) = .false.
    lfcscatter(icfg) = .false.
    lindintptr(icfg) = .false.
    lindvacptr(icfg) = .false.
    linitcfg(icfg) = .false.
    lksorigin(icfg) = .false.
    llowered(icfg) = .false.
    lmdconstrain(icfg) = .false.
    lnebvaryspring(icfg) = .false.
    lomega(icfg) = .false.
    lopfc(icfg) = .false.
    lopfreg(1:3*maxregion,icfg) = .false.
    lpotsitecartcfg(icfg) = .false.
    lreldin(icfg) = .false.
    lsumcoreshell(icfg) = .false.
    lsymset(icfg) = .false.
    ltrannoise(icfg) = .false.
    ltrantherm(icfg) = .false.
    lufree(icfg) = .false.
    lvecin(icfg) = .false.
    lvecpin(icfg) = .false.
    ifso(icfg) = 0
    ifhr(icfg) = 0
    iflags(icfg) = 0
    iperm(icfg) = 1
    iufree(icfg) = 0
    ivso(1:3,icfg) = 0
    maxmodecfg(icfg) = 0
    minmodecfg(icfg) = 1
    n1con(icfg) = 1
    n1var(icfg) = 0
    nascfg(icfg) = 0
    ncellmaxcfg(icfg) = 0
    ncellmincfg(icfg) = 0
    ncorepcfg(icfg) = 0
    nbornstep(icfg) = 0
    ndcentyp(icfg) = 0
    nd2cellcfg(1:3,icfg) = 0
    ndimen(icfg) = 0
    ndistancereset(icfg) = 1
    nebspring(icfg) = 0.00005_dp
    nebspringmin(icfg) = 0.00005_dp
    nebfinalcell(1:6,icfg) = 0.0_dp
    nebfinalradius(1:maxat,icfg) = 0.0_dp
    nebfinalxyz(1:3,1:maxat,icfg) = 0.0_dp
    neiglow(icfg) = 0
    neighigh(icfg) = 0
    nensemble(icfg) = 1
    nfixatomtypecfg(icfg) = 0
    nfixatomcfg(icfg) = 0
    ngocfg(icfg) = 1
    ninternalmaxcfg(icfg) = 0
    ninternalmincfg(icfg) = 0
    nmdconstrainatom(1:2,icfg) = 0
    nmdeq(icfg) = 0
    nmdprod(icfg) = 0
    nmdsamp(icfg) = 0
    nmdvelmode(icfg) = - 1
    nmdvelmodp(icfg) = - 1
    nmdwrite(icfg) = 0
    nmolcfg(icfg) = 0
    nmolatomcfg(1:maxmol,icfg) = 0
    nmolatomtotcfg(icfg) = 0
    nmolcorecfg(1:maxmol,icfg) = 0
    nnebreplica(icfg) = 0
    nomegastep(icfg) = 0
    nobsmodecfg(icfg) = 0
    norigkpt(icfg) = 0
    npotptcfg(icfg) = 0
    npotsitescfg(icfg) = 0
    nprojcfg(icfg) = 0
    nprojdef(icfg) = 0
    nqkgrid(1:3,icfg) = 24
    nreg1cfg(icfg) = 0
    nintecfg(icfg) = 0
    nvacacfg(icfg) = 0
    nregions(icfg) = 1
    nsasexcludemax(icfg) = - 1
    nsasexcludemin(icfg) = - 1
    nccscfg(icfg) = 1
    nshcfg(icfg) = 1
    symnocfg(icfg) = 1
    nspcg(icfg) = 1
    nsregion2(icfg) = 0
    nsuper(1,icfg) = 1
    nsuper(2,icfg) = 1
    nsuper(3,icfg) = 1
    nsuperghost(1,icfg) = 1
    nsuperghost(2,icfg) = 1
    nsuperghost(3,icfg) = 1
    ntempramp(icfg) = 0
    ntempstp(1:maxtempramp,icfg) = 0
    ntempstpstart(1:maxtempramp,icfg) = 0
    ncscan(icfg) = 0
    ntran(icfg) = 0
    ntwistcfg(icfg) = 0
    nummodecfg(icfg) = 0
    nzmolcfg(icfg) = 1
    nks(1:3,icfg) = 0
    nksala(1:3,icfg) = 0
    nvarcfg(icfg) = 0
    ordersuper(icfg) = 1
    QMMMmode(icfg) = 0
    nmdconstraindist(icfg) = 0.0_dp
    presscfg(icfg) = 0.0_dp
    pr_conscfg(icfg) = 0.0_dp
    qonsastarget(icfg) = 0.0_dp
    omega(icfg) = 0.0_dp
    omegadamping(icfg) = 5.0_dp
    omegadir(1:6,icfg) = 0.0_dp
    omegadirtype(icfg) = 1
    omegastep(icfg) = 0.0_dp
    qpres(icfg) = 0.1_dp
    qtemp(icfg) = 0.1_dp
    ramandir(1:6,icfg) = 1.0_dp
    ramandirtype(icfg) = 1
    reg1(icfg) = 0.0_dp
    reg2(icfg) = 0.0_dp
    reg1last(icfg) = 0.0_dp
    reg2a1(icfg) = 0.0_dp
    rufree(icfg) = 0.0_dp
    rvcfg(1:3,1:3,icfg) = 0.0_dp
    rvpcfg(1:3,1:3,icfg) = 0.0_dp
    sbulkecfg(icfg) = 0.0_dp
    shift(icfg) = 0.0_dp
    shscalecfg(icfg) = 1.0_dp
    straincfg(1:6,icfg) = 0.0_dp
    stresscfg(1:6,icfg) = 0.0_dp
    taubcfg(icfg) = 1.0_dp
    tautcfg(icfg) = 1.0_dp
    tempcfg(icfg) = 0.0_dp
    tempstp(1:maxtempramp,icfg) = 0.0_dp
    tmdeq(icfg) = 0.0_dp
    tmdforcestart(icfg) = 0.0_dp
    tmdforcestop(icfg) = 0.0_dp
    tmdprod(icfg) = 0.0_dp
    tmdsamp(icfg) = 0.0_dp
    tmdscale(icfg) = 0.0_dp
    tmdscint(icfg) = 0.0_dp
    tmdwrite(icfg) = 0.0_dp
    totalchargecfg(icfg) = 0.0_dp
    trannoise(icfg) = 0.0_dp
    trantherm(icfg) = 0.0_dp
    tstep(icfg) = 0.0_dp
    xdcent(icfg) = 0.0_dp
    ydcent(icfg) = 0.0_dp
    zdcent(icfg) = 0.0_dp
    cscan(1:6,icfg) = 0.0_dp
    xtran(icfg) = 0.0_dp
    ytran(icfg) = 0.0_dp
    ztran(icfg) = 0.0_dp
    xufree(1:3,icfg) = 0.0_dp
    nxpg(icfg) = 0
    nypg(icfg) = 0
    nzpg(icfg) = 0
    xminpg(icfg) = 0.0_dp
    yminpg(icfg) = 0.0_dp
    zminpg(icfg) = 0.0_dp
    xmaxpg(icfg) = 1.0_dp
    ymaxpg(icfg) = 1.0_dp
    zmaxpg(icfg) = 1.0_dp
    xmaxcfg(1:3,icfg) = 1.0_dp
    xmincfg(1:3,icfg) = 0.0_dp
!
    molcomcfg(1:3,1:maxmol,icfg) = 0.0_dp
    molQcfg(1:3,1:maxmol,icfg) = 0.0_dp
    molQxyzcfg(1:3,1:maxat,icfg) = 0.0_dp
!
    hmssg(1,icfg) = '('
    hmssg(2,icfg) = 'u'
    hmssg(3,icfg) = 'n'
    hmssg(4,icfg) = 'k'
    hmssg(5,icfg) = 'n'
    hmssg(6,icfg) = 'o'
    hmssg(7,icfg) = 'w'
    hmssg(8,icfg) = 'n'
    hmssg(9,icfg) = ')'
    do j = 10,16
      hmssg(j,icfg) = ' '
    enddo
    names(icfg) = ' '
    do j = 1,maxregion
      nregiontype(j,icfg) = 0
      if (j.eq.2) then
        lregionrigid(j,icfg) = .true. 
      else
        lregionrigid(j,icfg) = .false.
      endif
    enddo
!
    ropcfg(1:3,1:3,1,icfg) = 0.0_dp
    ropcfg(1,1,1,icfg) = 1.0_dp
    ropcfg(2,2,1,icfg) = 1.0_dp
    ropcfg(3,3,1,icfg) = 1.0_dp
    vitcfg(1:3,1,icfg) = 0.0_dp
!
    lfieldcfg(icfg) = .false.
    ntdfieldcfg(icfg) = 0
    fieldcfg(icfg) = 0.0_dp
    fielddirectioncfg(1,icfg) = 0.0_dp
    fielddirectioncfg(2,icfg) = 0.0_dp
    fielddirectioncfg(3,icfg) = 1.0_dp
    td_fieldcfg(1,1:maxtdfield,icfg) = 0.0_dp
    td_fieldcfg(2,1:maxtdfield,icfg) = 0.0_dp
    td_fieldcfg(3,1:maxtdfield,icfg) = 0.0_dp
    td_fielddirectioncfg(1,1:maxtdfield,icfg) = 0.0_dp
    td_fielddirectioncfg(2,1:maxtdfield,icfg) = 0.0_dp
    td_fielddirectioncfg(3,1:maxtdfield,icfg) = 1.0_dp
    tmdfieldstart(icfg) = 0.0_dp
    tmdfieldstop(icfg) = 0.0_dp
!
    lshearforcecfg(icfg) = .false.
!
    lradialcfg(icfg) = .false.
    radialKcfg(icfg) = 0.0_dp
    radialXYZcfg(1:3,icfg) = 0.0_dp
!
!  KIM arrays
!
    lkim_model_cfg_OK(icfg,1:maxkimmodel) = .true.
!
    lomega_af_in(icfg) = .false.
    omega_af(icfg) = 0.0_dp
    B_pr_cfg(icfg) = 0.0_dp
    v_s_cfg(icfg) = 0.0_dp
    v_p_cfg(icfg) = 0.0_dp
  endif
!
  return
  end
