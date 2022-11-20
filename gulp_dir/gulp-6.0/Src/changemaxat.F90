  subroutine changemaxat
!
!  Alters the size of the arrays associated with maxat
!
!  11/06 NEB modifications added
!  11/06 x/y/zfracimage arrays added
!   2/07 nbondedtype added
!   3/07 Chemshell changes added
!   5/07 Arrays for saving fractional coordinates added
!   5/07 Partial occupancy arrays added
!   6/07 Arrays for spatial decomposition added
!   7/07 Arrays for spatial decomposition checked against nspcellattot to
!        avoid incorrect lowering of size
!   1/08 ltrialatom added
!   3/08 qreaxFF array added for reaxFF charges
!  10/08 COSMO modifications added
!  11/08 scrho changed to be a 2-D array for benefit of MEAM
!   1/09 Integer datatypes all explicitly declared
!   1/09 Core-shell vector array added
!   3/09 neamfnspecptr added
!   6/09 site energy array added
!   6/09 New molecule indexing arrays added
!   9/10 Neutron scattering modifications added
!   9/10 Initialisations now performed in a subroutine
!   4/11 Frequency array now set in separate subroutine
!   5/12 Atomic stress array added
!   6/12 fobsmode added
!   7/11 nspcell2atptr added
!   7/12 nphonatptr and nphonatrptr added
!  10/12 Modified for OpenKIM arrays
!   9/13 Parallel arrays for cores and shells on node added
!   9/13 Raman susceptibilities added
!   8/14 Call to changemaxfreqat removed
!  12/14 nbs, nbss and nbsptr moved to shells module
!   2/15 nqatomcell and qatomxyz added
!   9/16 ncoshptr added
!  11/16 phonon pointers added for local node
!   1/17 optimisation atom pointer arrays added
!   1/17 neemlocptr added
!   2/17 Some charge arrays moved to changemaxatloc
!   2/17 Variable parallel arrays added
!   4/17 ChemShell interaction modified
!   6/17 Old ChemShell restored as a compile option
!   8/17 Modified so that changemaxr1at and changemaxat don't interfere
!        with the size of common arrays
!  10/17 Absolute coordinate flags added
!  11/17 msd arrays added
!   1/18 icosxsp, icosysp, icoszsp added
!   5/18 nbondqb added
!   5/18 nqrnow added
!   8/18 Modified for version 2 of OpenKIM
!   2/19 x0 removed
!   3/19 iopt changed to ioptindex and iopttype
!   7/19 Fractional in box arrays added
!   8/19 molatom option added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!  10/19 nasymnomolptr added
!  10/19 Langevin damping of dipoles added
!  11/19 natinmol added
!   3/20 nchargeptr added
!   3/20 Fractional coordinates for dump added
!   4/20 Pointer arrays added for numatnomol
!   4/20 Restart arrays added for rigid molecules
!   4/20 maxfqat changed to maxfreq since it is no longer just atom based
!   6/20 numatnomolonnodeptr added
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
!  Julian Gale, CIC, Curtin University, June 2020
!
  use cellmultipole,  only : nboxat
  use configurations, only : maxcfg
  use control,        only : lcosmo, lscatter
  use cosmic
  use current
  use derivatives,    only : atomicstress, sumatomicstress
  use distances,      only : lStoreVectors, ndistance, ndistanceind, ndistancemolonly, distlself, ndistanceij
  use distances,      only : icosxs, icosys, icoszs
  use eam,            only : maxmeamcomponent
  use eemdata,        only : nqrnow
  use energies,       only : siteenergy
#ifdef OLDCS
  use gulpchemsh,     only : shell_force,ichemsh_qm
#endif
  use ksample,        only : maxkpt
  use moldyn,         only : xabsco, yabsco, zabsco, labsco
  use molecule,       only : natmol, nmolind, xfsave, yfsave, zfsave, ixshift, iyshift, izshift
  use molecule,       only : nmollist, nmollistcfg, natinmol, molQxyzcfg
  use montecarlo,     only : ltrialatom
  use g_neb
  use observables,    only : fobsmode, maxobsmode
  use optimisation,   only : lopf, noptatptr, noptatrptr, noptatlocptr, noptatlocrptr
  use parallel,       only : atom2local, atom2node, node2atom, ncoonnodeptr, nshonnodeptr
  use parallel,       only : atom2locala, atom2nodea, node2atoma
  use parallel,       only : atom2localv, atom2nodev, node2atomv
  use parallel,       only : node2var, nvar2node, node2atomv, nvar2local
  use partial
  use phononatoms
  use polarise,       only : dpolar, qpolar, dpolarmax
  use potentialxyz
  use properties,     only : ramanasus
  use reallocate
  use reaxFFdata,     only : qreaxFF
  use scatterdata,    only : scatlencoh, scatleninc, sofomega
  use shells,         only : ncsptr, ncoptr, ncoshptr, nshptr, ratiom, csvector, nbsptr
  use shellextrapolation
  use spatial,        only : natomcell, xinbox, yinbox, zinbox, nspcellatptr, nspcellatptrcell, nspcell2atptr
  use spatial,        only : maxnspcellattot, xfinbox, yfinbox, zfinbox
  use spatialbo,      only : natomcellbo, xinboxbo, yinboxbo, zinboxbo, nspcellatptrbo, nspcellatptrcellbo
  use spatialbo,      only : maxnspcellattotbo, xfinboxbo, yfinboxbo, zfinboxbo
  use sutton,         only : scrho, scrho12
  use velocities 
  implicit none
!
!  Local variables
!
  integer(i4)       :: ierror, i
  integer(i4), save :: oldmaxat = 0
!
  call realloc(atom2local,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','atom2local')
  call realloc(atom2locala,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','atom2locala')
  call realloc(atom2localv,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','atom2localv')
  call realloc(atom2node,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','atom2node')
  call realloc(atom2nodea,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','atom2nodea')
  call realloc(atom2nodev,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','atom2nodev')
  call realloc(atomicstress,6_i4,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','atomicstress')
  call realloc(c6a,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','c6a')
  call realloc(c6f,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','c6f')
  call realloc(cna,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','cna')
  call realloc(cnf,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','cnf')
  call realloc(iatn,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','iatn')
  call realloc(icosx,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','icosx')
  call realloc(icosy,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','icosy')
  call realloc(icosz,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','icosz')
  call realloc(icosxsp,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','icosxsp')
  call realloc(icosysp,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','icosysp')
  call realloc(icoszsp,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','icoszsp')
  call realloc(ixshift,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','ixshift')
  call realloc(iyshift,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','iyshift')
  call realloc(izshift,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','izshift')
  call realloc(ioptindex,4_i4*maxat+6_i4,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','ioptindex')
  call realloc(iopttype,4_i4*maxat+6_i4,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','iopttype')
!
  if (4_i4*maxat+6_i4.gt.size(node2var)) then
    call realloc(node2var,4_i4*maxat+6_i4,ierror)
    if (ierror.ne.0) call outofmemory('changemaxat','node2var')
    call realloc(nvar2local,4_i4*maxat+6_i4,ierror)
    if (ierror.ne.0) call outofmemory('changemaxat','nvar2local')
    call realloc(nvar2node,4_i4*maxat+6_i4,ierror)
    if (ierror.ne.0) call outofmemory('changemaxat','nvar2node')
  endif
!
  call realloc(labsco,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','labsco')
  call realloc(lopf,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','lopf')
  call realloc(ltrialatom,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','ltrialatom')
  call realloc(mass,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','mass')
  call realloc(molQxyzcfg,3_i4,maxat,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','molQxyzcfg')
  call realloc(msdx,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','msdx')
  call realloc(msdy,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','msdy')
  call realloc(msdz,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','msdz')
  call realloc(nasymnomolptr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','nasymnomolptr')
  call realloc(nasymnomolrptr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','nasymnomolrptr')
  call realloc(ncorenomolptr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','ncorenomolptr')
  call realloc(ncorenomolrptr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','ncorenomolrptr')
  call realloc(numatnomolonnodeptr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','numatnomolonnodeptr')
  call realloc(numatnomolptr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','numatnomolptr')
  call realloc(numatnomolrptr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','numatnomolrptr')
  call realloc(nat,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','nat')
  call realloc(natmol,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','natmol')
  call realloc(natinmol,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','natinmol')
  call realloc(natomcell,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','natomcell')
  call realloc(natomcellbo,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','natomcellbo')
  call realloc(natype,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','natype')
  call realloc(nbonds,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','nbonds')
  call realloc(nbonded,maxbond,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','nbonded')
  call realloc(nbondedtype,2_i4,maxbond,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','nbondedtype')
  call realloc(nbondind,maxbond,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','nbondind')
  call realloc(nbondqb,maxbond,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','nbondqb')
  call realloc(nboxat,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','nboxat')
  call realloc(nchargeptr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','nchargeptr')
  call realloc(nchargec6ptr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','nchargec6ptr')
  call realloc(ncsptr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','ncsptr')
  call realloc(nebfinalradius,maxat,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','nebfinalradius')
  call realloc(nebfinalxyz,3_i4,maxat,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','nebfinalxyz')
  call realloc(nebreplicaradius,maxat,maxnebreplicatot,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','nebreplicaradius')
  call realloc(nebreplicaxyz,3_i4,maxat,maxnebreplicatot,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','nebreplicaxyz')
  call realloc(neamspecptr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','neamspecptr')
  call realloc(neamfnspecptr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','neamfnspecptr')
  call realloc(neemptr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','neemptr')
  call realloc(neemrptr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','neemrptr')
  call realloc(neemlocptr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','neemlocptr')
  call realloc(neemlocrptr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','neemlocrptr')
  call realloc(neqv,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','neqv')
  call realloc(nftype,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','nftype')
  call realloc(nmolind,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','nmolind')
  call realloc(nmollist,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','nmollist')
  call realloc(nmollistcfg,maxat,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','nmollistcfg')
  call realloc(nqrnow,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','nqrnow')
  call realloc(node2atom,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','node2atom')
  call realloc(node2atoma,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','node2atoma')
  call realloc(node2atomv,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','node2atomv')
  call realloc(ncoonnodeptr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','ncoonnodeptr')
  call realloc(nshonnodeptr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','nshonnodeptr')
  call realloc(noptatptr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','noptatptr')
  call realloc(noptatrptr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','noptatrptr')
  call realloc(noptatlocptr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','noptatlocptr')
  call realloc(noptatlocrptr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','noptatlocrptr')
  call realloc(nrela2f,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','nrela2f')
  call realloc(nrelf2a,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','nrelf2a')
  call realloc(nmolrela2f,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','nmolrela2f')
  call realloc(nmolrelf2a,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','nmolrelf2a')
  call realloc(nrotop,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','nrotop')
  call realloc(ncoptr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','ncoptr')
  call realloc(ncoshptr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','ncoshptr')
  call realloc(nshptr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','nshptr')
  call realloc(nspecptr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','nspecptr')
  call realloc(nspcell2atptr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','nspcell2atptr')
  call realloc(nspcellatptr,max(maxat,maxnspcellattot),ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','nspcellatptr')
  call realloc(nspcellatptrcell,max(maxat,maxnspcellattot),ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','nspcellatptrcell')
  call realloc(nspcellatptrbo,max(maxat,maxnspcellattotbo),ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','nspcellatptrbo')
  call realloc(nspcellatptrcellbo,max(maxat,maxnspcellattotbo),ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','nspcellatptrcellbo')
  call realloc(occua,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','occua')
  call realloc(occuf,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','occuf')
  call realloc(oxa,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','oxa')
  call realloc(oxf,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','oxf')
  call realloc(dpolar,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','dpolar')
  call realloc(dpolarmax,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','dpolarmax')
  call realloc(qpolar,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','qpolar')
  call realloc(qa,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','qa')
  call realloc(qf,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','qf')
  call realloc(qreaxFF,maxat+1_i4,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','qreaxFF')
  call realloc(rada,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','rada')
  call realloc(radf,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','radf')
  call realloc(ratiom,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','ratiom')
  call realloc(rmass,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','rmass')
  call realloc(scrho,maxmeamcomponent,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','scrho')
  call realloc(scrho12,maxmeamcomponent,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','scrho12')
  call realloc(siteenergy,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','siteenergy')
  call realloc(sumatomicstress,6_i4,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','sumatomicstress')
  call realloc(xabsco,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','xabsco')
  call realloc(yabsco,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','yabsco')
  call realloc(zabsco,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','zabsco')
  call realloc(xalat,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','xalat')
  call realloc(yalat,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','yalat')
  call realloc(zalat,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','zalat')
  call realloc(xclat,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','xclat')
  call realloc(yclat,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','yclat')
  call realloc(zclat,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','zclat')
  call realloc(xafrac,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','xafrac')
  call realloc(yafrac,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','yafrac')
  call realloc(zafrac,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','zafrac')
  call realloc(xfdmp,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','xfdmp')
  call realloc(yfdmp,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','yfdmp')
  call realloc(zfdmp,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','zfdmp')
  call realloc(xfrac,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','xfrac')
  call realloc(yfrac,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','yfrac')
  call realloc(zfrac,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','zfrac')
  call realloc(xfracimage,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','xfracimage')
  call realloc(yfracimage,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','yfracimage')
  call realloc(zfracimage,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','zfracimage')
  call realloc(xfsave,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','xfsave')
  call realloc(yfsave,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','yfsave')
  call realloc(zfsave,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','zfsave')
  call realloc(xinitial,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','xinitial')
  call realloc(yinitial,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','yinitial')
  call realloc(zinitial,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','zinitial')
  call realloc(xstore,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','xstore')
  call realloc(ystore,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','ystore')
  call realloc(zstore,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','zstore')
  call realloc(rstore,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','rstore')
  call realloc(velx,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','velx')
  call realloc(vely,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','vely')
  call realloc(velz,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','velz')
  call realloc(v2xyz,6_i4,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','v2xyz')
  call realloc(v2xyz12,6_i4,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','v2xyz12')
  call realloc(vx,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','vx')
  call realloc(vy,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','vy')
  call realloc(vz,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','vz')
  call realloc(vx12,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','vx12')
  call realloc(vy12,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','vy12')
  call realloc(vz12,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','vz12')
  call realloc(x2,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','x2')
  call realloc(y2,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','y2')
  call realloc(z2,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','z2')
  call realloc(x3,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','x3')
  call realloc(y3,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','y3')
  call realloc(z3,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','z3')
  call realloc(x4,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','x4')
  call realloc(y4,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','y4')
  call realloc(z4,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','z4')
  call realloc(x5,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','x5')
  call realloc(y5,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','y5')
  call realloc(z5,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','z5')
  call realloc(bornq,3_i4,3_i4,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','bornq')
  call realloc(xinbox,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','xinbox')
  call realloc(yinbox,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','yinbox')
  call realloc(zinbox,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','zinbox')
  call realloc(xfinbox,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','xfinbox')
  call realloc(yfinbox,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','yfinbox')
  call realloc(zfinbox,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','zfinbox')
  call realloc(xinboxbo,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','xinboxbo')
  call realloc(yinboxbo,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','yinboxbo')
  call realloc(zinboxbo,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','zinboxbo')
  call realloc(xfinboxbo,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','xfinboxbo')
  call realloc(yfinboxbo,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','yfinboxbo')
  call realloc(zfinboxbo,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','zfinboxbo')
  call realloc(xshellsave,maxextrapol,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','xshellsave')
  call realloc(yshellsave,maxextrapol,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','yshellsave')
  call realloc(zshellsave,maxextrapol,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','zshellsave')
  call realloc(csvector,3_i4,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','csvector')
!
  call realloc(fobsmode,3_i4,maxat,maxobsmode,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','fobsmode')
!
  call realloc(nbsptr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','nbsptr')
  call realloc(ibocptr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','ibocptr')
  call realloc(ibocshptr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','ibocshptr')
  call realloc(iocptr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','iocptr')
  call realloc(iocshptr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','iocshptr')
!
!  Phonon atom pointers
!
  call realloc(nphonatptr,maxat,ierror)
  if (ierror.gt.0) call outofmemory('changemaxat','nphonatptr')
  call realloc(nphonatcptr,maxat,ierror)
  if (ierror.gt.0) call outofmemory('changemaxat','nphonatcptr')
  call realloc(nphonatsptr,maxat,ierror)
  if (ierror.gt.0) call outofmemory('changemaxat','nphonatsptr')
  call realloc(nphonatmptr,maxat,ierror)
  if (ierror.gt.0) call outofmemory('changemaxat','nphonatmptr')
  call realloc(nphonatrptr,maxat,ierror)
  if (ierror.gt.0) call outofmemory('changemaxat','nphonatrptr')
  call realloc(nphonatrcptr,maxat,ierror)
  if (ierror.gt.0) call outofmemory('changemaxat','nphonatrcptr')
  call realloc(nphonatrsptr,maxat,ierror)
  if (ierror.gt.0) call outofmemory('changemaxat','nphonatrsptr')
  call realloc(nphonatrmptr,maxat,ierror)
  if (ierror.gt.0) call outofmemory('changemaxat','nphonatrmptr')
!
!  Phonon atom pointers on local node
!
  call realloc(nphonatonnodeptr,maxat,ierror)
  if (ierror.gt.0) call outofmemory('changemaxat','nphonatonnodeptr')
  call realloc(nphonatonnodecptr,maxat,ierror)
  if (ierror.gt.0) call outofmemory('changemaxat','nphonatonnodecptr')
  call realloc(nphonatonnodesptr,maxat,ierror)
  if (ierror.gt.0) call outofmemory('changemaxat','nphonatonnodesptr')
  call realloc(nphonatonnodemptr,maxat,ierror)
  if (ierror.gt.0) call outofmemory('changemaxat','nphonatonnodemptr')
  call realloc(nphonatonnoderptr,maxat,ierror)
  if (ierror.gt.0) call outofmemory('changemaxat','nphonatonnoderptr')
  call realloc(nphonatonnodercptr,maxat,ierror)
  if (ierror.gt.0) call outofmemory('changemaxat','nphonatonnodercptr')
  call realloc(nphonatonnodersptr,maxat,ierror)
  if (ierror.gt.0) call outofmemory('changemaxat','nphonatonnodersptr')
  call realloc(nphonatonnodermptr,maxat,ierror)
  if (ierror.gt.0) call outofmemory('changemaxat','nphonatonnodermptr')
!
!  Frequencies
!
  call changemaxfreq(3_i4*maxat)
!
  call realloc(ramanasus,3_i4,3_i4,3_i4,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','ramanasus')
!
!  Neutron scattering
!
  call realloc(scatlencoh,maxat,ierror)
  if (ierror.gt.0) call outofmemory('changemaxat','scatlencoh')
  call realloc(scatleninc,maxat,ierror)
  if (ierror.gt.0) call outofmemory('changemaxat','scatleninc')
  if (lscatter) then
    call realloc(sofomega,3_i4*maxat,maxkpt,ierror)
    if (ierror.gt.0) call outofmemory('changemaxat','sofomega')
  endif
#ifdef OLDCS
  if (ichemsh_qm==1) then
    call realloc(shell_force,3_i4,maxat,ierror)
    if (ierror.ne.0) call outofmemory('changemaxat','shell_force')
  endif
#endif
!
!  Change dimensions of interatomic distance storage vectors if not direct algorithm
!
  if (lStoreVectors) then
    call realloc(icosxs,maxat,ierror)
    if (ierror.ne.0) call outofmemory('changemaxat','icosxs')
    call realloc(icosys,maxat,ierror)
    if (ierror.ne.0) call outofmemory('changemaxat','icosys')
    call realloc(icoszs,maxat,ierror)
    if (ierror.ne.0) call outofmemory('changemaxat','icoszs')
    call realloc(ndistance,maxat,ierror)
    if (ierror.ne.0) call outofmemory('changemaxat','ndistance')
    call realloc(ndistanceij,maxat,maxat,ierror)
    if (ierror.ne.0) call outofmemory('changemaxat','ndistanceij')
    call realloc(ndistanceind,maxat,ierror)
    if (ierror.ne.0) call outofmemory('changemaxat','ndistanceind')
    call realloc(ndistancemolonly,maxat,maxat,ierror)
    if (ierror.ne.0) call outofmemory('changemaxat','ndistancemolonly')
    call realloc(distlself,maxat,maxat,ierror)
    if (ierror.ne.0) call outofmemory('changemaxat','distlself')
  endif
!
!  Change maxat dimensions in COSMO
!
  if (lcosmo) call changemaxatcosmo
!
!  Initialise new parts of data arrays
!
  if (maxat.gt.oldmaxat) then
    do i = oldmaxat+1,maxat
      call initmaxatdefaults(i)
    enddo
  endif
!
!  Save current value of maxat for next call
!
  oldmaxat = maxat
!
  return
  end
