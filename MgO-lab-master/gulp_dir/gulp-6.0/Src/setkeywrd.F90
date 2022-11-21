  subroutine setkeyword
!
!  Sets general flags based on keywords
!
!   2/04 Created from getkeyword
!   7/05 lSandM added
!   7/05 Setting of linten flag improved to remove ambiguity
!   9/06 libdump keyword added
!  11/06 NEB modifications added
!   1/07 Gasteiger charges added
!   3/07 lPrintEAM keyword added
!   3/07 lPrintTwo keyword added
!   3/07 More robust checking of keywords added by ensuring
!        string comes at the start of the word
!   3/07 lpreserveQ added
!   4/07 Conjugate gradients set if running in parallel
!   5/07 Mean KE option added
!   5/07 qbond keyword added
!   7/07 lmeta added
!  10/08 COSMO/COSMIC keywords merged in 
!  10/08 Error in logic for lmodco setting corrected
!  10/08 Error in logic for other keywords also corrected
!  11/08 lPrintFour keyword added
!  11/08 lPrintThree keyword added
!   2/09 lconj only set to be true for > one processor if not LM-BFGS
!   6/09 Module name changed from three to m_three
!   6/09 PDF keywords added
!   7/09 Symmetry turned off for derivatives in reaxFF case
!  12/09 pregionforce keyword added
!   4/10 qtpie keyword added
!   6/10 Hopping keyword added
!   8/10 lconvert, lphase, lcutbelow keywords removed
!   8/10 Keyword settings adjusted for PDF case
!   8/10 lfix1atom added
!  10/10 Symmetry turned off for EDIP potentials
!  12/10 Hiding of shells added
!   1/11 Force minimisation option added
!   2/11 Keyword added to turn off ReaxFF charges
!   3/11 lstressout added
!   9/11 Metadynamics internal code replaced with Plumed
!   9/11 Madelung correction added
!  11/11 eregion keyword added
!   5/12 Atomic stress keyword added
!   6/12 Thermal conductivity added
!   7/12 Oldinten keyword added
!   9/12 Pacha added
!   9/12 Eckart transformation added 
!  10/12 Option added to allow two atoms co-exist on a site with an
!        occupancy greater than one.
!  12/12 site_energy keyword added
!   7/13 Modified so that multiple calls do not overwrite flags set by options
!   8/13 linclude_imaginary added
!   8/13 lfindsym added
!   8/13 Raman keyword added
!  10/13 lsopt keyword added
!   3/14 Harmonic relaxation keyword added
!   3/14 qiter added
!   6/14 lnoexclude flag set
!   8/14 Calculation of group velocities added
!   9/14 Typo in setting of flags for prt_ keywords fixed
!   1/15 nowrap added as a pseudonym for nomod
!   3/15 lnoqeem added
!   8/15 New lower options added
!   5/16 lspme added
!   7/16 lShengBTE3rd added 
!   8/16 Flags set for SPME to indicate that second and third
!        derivatives are not available yet
!   1/17 lPrintSix added
!   1/17 Condition on using conjugate gradients in parallel removed
!   3/17 loldvarorder added
!   5/17 lnewdefalg added
!   6/17 Module files renamed to gulp_files
!   8/17 lnumerical added
!  10/17 Delta_dipole added
!  11/17 lmsd added
!  12/17 Fastfd keyword added
!   1/18 ghostcell keyword added
!   1/18 Grueneisen parameters added
!   1/18 Trace added
!   2/18 lalamode added
!   3/18 Frame keyword added
!   4/18 Split bond EEM added
!   6/18 lallbonds added
!   6/18 Strain cell added
!   7/18 Symmetry turned off derivatives when using split bond charges
!  11/18 Finite strain derivative flag added
!  12/18 lusevectors added
!   2/19 Rigid molecules added
!   7/19 lceigen added
!   9/19 lnotorXduplicates flag added
!   9/19 lnoshellzero added
!  11/19 lxray added
!   3/20 Modifications for rigid molecules added
!   3/20 lPrintVar added
!   4/20 lmolcom_mass added
!   4/20 lnorotate added
!   4/20 For rigid molecules set lusevectors to be true for restarting
!   5/20 lphasecom added
!   5/20 lfix forced to be true for rigid molecules
!   7/20 lkcal and lkjmol added
!   9/20 lpolarisation added
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
!  Julian Gale, CIC, Curtin University, September 2020
!
  use bondorderdata,  only : nboQ
  use control
  use cosmic,         only : lcosmic, lsegsmooth
  use derivatives,    only : lfinitestrain
  use distances,      only : lStoreVectors
  use eam,            only : lPrintEAM
  use element
  use gulp_files,     only : lShengBTE3rd
  use fitting
  use four,           only : lPrintFour
  use library,        only : llibsymdump
  use m_pdfneutron,   only : lmakeeigarray, lcoreinfo, lpdf
  use m_pdfneutron,   only : lnowidth, lpartial, lfreqcut, lkeepcut, lnoksym
  use m_pdfneutron,   only : lpdfout
  use m_three,        only : lPrintThree
  use mdlogic
  use molecule
  use g_neb,          only : lnebclimbingimage, lnebdoublynudged
  use optimisation
  use phonout,        only : linclude_imaginary
  use six,            only : lPrintSix
  use spme,           only : lspme
  use symmetry
  use synchro,        only : lfixtangent
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  use two,            only : lPrintTwo
  use wolfcosmo,      only : lPureCoulomb0D
  implicit none
#ifdef TRACE
  call trace_in('setkeywrd')
#endif
!****************************
!  Set flags from keywords  *
!****************************
!
!  Keywords with default = .false.
!
  if (.not.lalamode) lalamode = (index(keyword,' ala').ne.0.or.index(keyword,'ala').eq.1)
  if (.not.lallbonds) lallbonds = (index(keyword,' allb').ne.0.or.index(keyword,'allb').eq.1)
  if (.not.lallowgt1) lallowgt1 = (index(keyword,' allo').ne.0.or.index(keyword,'allo').eq.1)
  if (.not.langle) langle = (index(keyword,' angl').ne.0.or.index(keyword,'angl').eq.1)
  if (.not.lanneal) lanneal = (index(keyword,' anne').ne.0.or.index(keyword,'anne').eq.1)
  if (.not.latomicstress) latomicstress = (index(keyword,' atom').ne.0.or.index(keyword,'atom').eq.1)
  if (.not.laver) laver = (index(keyword,' aver').ne.0.or.index(keyword,'aver').eq.1)
  if (.not.lbond) lbond = (index(keyword,' bond').ne.0.or.index(keyword,'bond').eq.1)
  if (.not.lbroad) lbroad = (index(keyword,' broa').ne.0.or.index(keyword,'broa').eq.1)
  if (.not.lbulknoopt) lbulknoopt = (index(keyword,' bulk').ne.0.or.index(keyword,'bulk').eq.1)
  if (.not.lc6) lc6 = ((index(keyword,' c6').ne.0.or.index(keyword,'c6').eq.1).and.index(keyword,'noel').eq.0)
  if (.not.lceigen) lceigen = (index(keyword,' ceig').ne.0.or.index(keyword,'ceig').eq.1)
  if (.not.lcello) lcello = (index(keyword,' cell').ne.0.or.index(keyword,'cell').eq.1)
  if (.not.lcomp) lcomp = (index(keyword,' comp').ne.0.or.index(keyword,'comp').eq.1)
  if (.not.lconj) lconj = (index(keyword,' conj').ne.0.or.index(keyword,'conj').eq.1)
  if (.not.lconp) lconp = (index(keyword,' conp').ne.0.or.index(keyword,'conp').eq.1)
  if (.not.lconv) lconv = (index(keyword,' conv').ne.0.or.index(keyword,'conv').eq.1)
  if (.not.lcosmic) lcosmic = (index(keyword,'cosmi').ne.0)
  if (.not.lcosmo) lcosmo = (index(keyword,'cosmo').ne.0)
  if (.not.ldcharge) ldcharge = (index(keyword,' dcha').ne.0.or.index(keyword,'dcha').eq.1)
  if (.not.lddipole) lddipole = (index(keyword,' delt').ne.0.or.index(keyword,'delt').eq.1)
  if (.not.ldebug) ldebug = (index(keyword,' debu').ne.0.or.index(keyword,'debu').eq.1)
  if (.not.ldefect) ldefect = (index(keyword,' defe').ne.0.or.index(keyword,'defe').eq.1)
  if (.not.ldfp) ldfp = ((index(keyword,' dfp').ne.0.or.index(keyword,'dfp').eq.1).and.index(keyword,'bfgs').eq.0)
  if (.not.ldipole) ldipole = (index(keyword,' dipo').ne.0.or.index(keyword,'dipo').eq.1)
  if (.not.lStoreVectors) lStoreVectors = (index(keyword,' stor').ne.0.or.index(keyword,'stor').eq.1)
  if (.not.ldist) ldist = (index(keyword,' dist').ne.0.or.index(keyword,'dist').eq.1)
  if (.not.leckart) leckart = (index(keyword,' eck').ne.0.or.index(keyword,'eck').eq.1)
  if (.not.leem) leem = (index(keyword,' eem').ne.0.or.index(keyword,'eem').eq.1)
  if (.not.leembond) leembond = (index(keyword,' eemb').ne.0.or.index(keyword,'eemb').eq.1)
  if (.not.lefg) lefg = (index(keyword,' efg').ne.0.or.index(keyword,'efg').eq.1)
  if (.not.leigen) leigen = (index(keyword,' eige').ne.0.or.index(keyword,'eige').eq.1)
  if (.not.leregion) leregion = (index(keyword,' ereg').ne.0.or.index(keyword,'ereg').eq.1)
  if (.not.lfastfd) lfastfd = (index(keyword,' fast').ne.0.or.index(keyword,'fast').eq.1)
  if (.not.lfbfgs) lfbfgs = (index(keyword,' fbfg').ne.0.or.index(keyword,'fbfg').eq.1)
  if (.not.lfindsym) lfindsym = (index(keyword,' find').ne.0.or.index(keyword,'find').eq.1)
  if (.not.lfit) lfit = (index(keyword,' fit').ne.0.or.index(keyword,'fit').eq.1)
  if (.not.lforcemin) lforcemin = (index(keyword,' forc').ne.0.or.index(keyword,'forc').eq.1)
  if (.not.lframe) lframe = (index(keyword,' fram').ne.0.or.index(keyword,'fram').eq.1)
  if (.not.lfree) lfree = (index(keyword,' free').ne.0.or.index(keyword,'free').eq.1)
  if (.not.lfreq) lfreq = (index(keyword,' freq').ne.0.or.index(keyword,'freq').eq.1)
  if (.not.lfreqout) lfreqout = (index(keyword,' nofr').eq.0.and.index(keyword,'nofr').ne.1)
  if (.not.lfixtangent) lfixtangent = (index(keyword,' tfix').ne.0.and.index(keyword,'tfix').eq.1)
  if (.not.lga) lga = (index(keyword,' gene').ne.0.or.index(keyword,'gene').eq.1)
  if (.not.lgasteiger) lgasteiger = (index(keyword,' gast').ne.0.or.index(keyword,'gast').eq.1)
  if (.not.lghost) lghost = (index(keyword,' ghos').ne.0.or.index(keyword,'ghos').eq.1)
  if (.not.lgrad) lgrad = (index(keyword,' grad').ne.0.or.index(keyword,'grad').eq.1)
  if (.not.lgroupvelocity) lgroupvelocity = (index(keyword,' gro').ne.0.or.index(keyword,'gro').eq.1)
  if (.not.lgrueneisen) lgrueneisen = (index(keyword,' gru').ne.0.or.index(keyword,'gru').eq.1)
  if (.not.lharmrelax) lharmrelax = (index(keyword,' ehar').ne.0.or.index(keyword,'ehar').eq.1)
  if (.not.lhex) lhex = (index(keyword,' hex').ne.0.or.index(keyword,'hex').eq.1)
  if (.not.lhideshells) lhideshells = (index(keyword,' hide').ne.0.or.index(keyword,'hide').eq.1)
  if (.not.linclude_imaginary) linclude_imaginary = (index(keyword,' incl').ne.0.or.index(keyword,'incl').eq.1)
  if (.not.linten) linten = (index(keyword,' inte').ne.0.or.index(keyword,'inte').eq.1)
  if (.not.lmsd) lmsd = (index(keyword,' msd').ne.0.or.index(keyword,'msd').eq.1)
  if (.not.literativeQ) literativeQ = (index(keyword,' qite').ne.0.or.index(keyword,'qite').eq.1)
  if (.not.lkcal) lkcal = (index(keyword,' kcal').ne.0.or.index(keyword,'kcal').eq.1)
  if (.not.lkjmol) lkjmol = (index(keyword,' kjm').ne.0.or.index(keyword,'kjm').eq.1)
  if (.not.lkfull) lkfull = (index(keyword,' kful').ne.0.or.index(keyword,'kful').eq.1)
  if (.not.llbfgs) llbfgs = (index(keyword,' lbfg').ne.0.or.index(keyword,'lbfg').eq.1)
  if (.not.llibsymdump) llibsymdump = (index(keyword,' libd').ne.0.or.index(keyword,'libd').eq.1)
  if (.not.llower) llower = (index(keyword,' lowe').ne.0.or.index(keyword,'lowe').eq.1)
  if (.not.lmadelung) lmadelung = (index(keyword,' made').ne.0.or.index(keyword,'made').eq.1)
  if (.not.lmarvreg2) lmarvreg2 = (index(keyword,' marv').ne.0.or.index(keyword,'marv').eq.1)
  if (.not.lmc) lmc = (index(keyword,' mont').ne.0.or.index(keyword,'mont').eq.1)
  if (.not.lmd) lmd = (index(keyword,' md').ne.0.or.index(keyword,'md').eq.1)
  if (.not.lmeanke) lmeanke = (index(keyword,' mean').ne.0.or.index(keyword,'mean').eq.1)
  if (.not.lminimage) lminimage = (index(keyword,' mini').ne.0.or.index(keyword,'mini').eq.1)
  if (.not.lmol) lmol = (index(keyword,' mol').ne.0.or.index(keyword,'mol').eq.1)
  if (.not.lmolq) lmolq = (index(keyword,' molq').ne.0.or.index(keyword,'molq').eq.1)
  if (.not.lmolmec) lmolmec = (index(keyword,' molm').ne.0.or.index(keyword,'molm').eq.1)
  if (.not.lmolfix) lmolfix = (index(keyword,' fix').ne.0.or.index(keyword,'fix').eq.1)
  if (.not.lneb) lneb = (index(keyword,' neb').ne.0.or.index(keyword,'neb').eq.1)
  if (.not.lnebclimbingimage) lnebclimbingimage = (index(keyword,' cineb').ne.0.or.index(keyword,'cineb').eq.1)
  if (.not.lnewdefalg) lnewdefalg = (index(keyword,' newda').ne.0.or.index(keyword,'newda').eq.1)
  if (.not.lnoautobond) lnoautobond = (index(keyword,' noau').ne.0.or.index(keyword,'noau').eq.1)
  if (.not.lnoenergy) lnoenergy = (index(keyword,' noen').ne.0.or.index(keyword,'noen').eq.1)
  if (.not.lnoexclude) lnoexclude = (index(keyword,' noex').ne.0.or.index(keyword,'noex').eq.1)
  if (.not.lnoflags) lnoflags = (index(keyword,' nofl').ne.0.or.index(keyword,'nofl').eq.1)
  if (.not.lnoqeem) lnoqeem = (index(keyword,' noqe').ne.0.or.index(keyword,'noqe').eq.1)
  if (.not.lnoreal) lnoreal = (index(keyword,' noreal').ne.0.or.index(keyword,'noreal').eq.1)
  if (.not.lnorecip) lnorecip = (index(keyword,' noreci').ne.0.or.index(keyword,'noreci').eq.1)
  if (.not.lnorotate) lnorotate = (index(keyword,' noro').ne.0.or.index(keyword,'noro').eq.1)
  if (.not.lnoshellzero) lnoshellzero = (index(keyword,' nosh').ne.0.or.index(keyword,'nosh').eq.1)
  if (.not.lnotorXduplicates) lnotorXduplicates = (index(keyword,' no4d').ne.0.or.index(keyword,'no4d').eq.1)
  if (.not.lnumdiag) lnumdiag = (index(keyword,' numd').ne.0.or.index(keyword,'numd').eq.1)
  if (.not.lnumerical) lnumerical = (index(keyword,' nume').ne.0.or.index(keyword,'nume').eq.1)
  if (.not.loldinten) loldinten = (index(keyword,' oldi').ne.0.or.index(keyword,'oldi').eq.1)
  if (.not.loldvarorder) loldvarorder = (index(keyword,' oldv').ne.0.or.index(keyword,'oldv').eq.1)
  if (.not.loptlower) loptlower = (index(keyword,' optl').ne.0.or.index(keyword,'optl').eq.1)
  if (.not.lopt) lopt = (index(keyword,' opti').ne.0.or.index(keyword,'opti').eq.1)
  if (.not.loptcellpar) loptcellpar = (index(keyword,' ocel').ne.0.or.index(keyword,'ocel').eq.1)
  if (.not.lpacha) lpacha = (index(keyword,' pac').ne.0.or.index(keyword,'pac').eq.1)
  if (.not.lphon) lphon = (index(keyword,' phon').ne.0.or.index(keyword,'phon').eq.1)
  if (.not.lpolarisation) lpolarisation = (index(keyword,' pol').ne.0.or.index(keyword,'pol').eq.1)
  if (.not.lposidef) lposidef = (index(keyword,' posi').ne.0.or.index(keyword,'posi').eq.1)
  if (.not.lpot) lpot = (index(keyword,' pot').ne.0.or.index(keyword,'pot').eq.1)
  if (.not.lpredict) lpredict = (index(keyword,' pred').ne.0.or.index(keyword,'pred').eq.1)
  if (.not.lpreserveQ) lpreserveQ = (index(keyword,' pres').ne.0.or.index(keyword,'pres').eq.1)
  if (.not.lPrintEAM) lPrintEAM = (index(keyword,' prt_eam').ne.0.or.index(keyword,'prt_eam').eq.1)
  if (.not.lPrintFour) lPrintFour = (index(keyword,' prt_fo').ne.0.or.index(keyword,'prt_fo').eq.1)
  if (.not.lPrintSix) lPrintSix = (index(keyword,' prt_s').ne.0.or.index(keyword,'prt_s').eq.1)
  if (.not.lPrintThree) lPrintThree = (index(keyword,' prt_th').ne.0.or.index(keyword,'prt_th').eq.1)
  if (.not.lPrintTwo) lPrintTwo = (index(keyword,' prt_two').ne.0.or.index(keyword,'prt_two').eq.1)
  if (.not.lPrintVar) lPrintVar = (index(keyword,' prt_var').ne.0.or.index(keyword,'prt_var').eq.1)
  if (.not.lprop) lprop = (index(keyword,' prop').ne.0.or.index(keyword,'prop').eq.1)
  if (.not.lPureCoulomb0D) lPureCoulomb0D = (index(keyword,'pure').ne.0)
  if (.not.lqbond) lqbond = (index(keyword,' qbon').ne.0.or.index(keyword,'qbon').eq.1)
  if (.not.lqeq) lqeq = (index(keyword,' qeq').ne.0.or.index(keyword,'qeq').eq.1)
  if (.not.lqtpie) lqtpie = (index(keyword,' qtp').ne.0.or.index(keyword,'qtp').eq.1)
  if (.not.lquicksearch) lquicksearch = (index(keyword,' quic').ne.0.or.index(keyword,'quic').eq.1)
  if (.not.lraman) lraman = (index(keyword,' rama').ne.0.or.index(keyword,'rama').eq.1)
  if (.not.lrelax) lrelax = ((index(keyword,' rela').ne.0.or.index(keyword,'rela').eq.1).and.lfit)
  if (.not.lregionforce) lregionforce = (index(keyword,' preg').ne.0.or.index(keyword,'preg').eq.1)
  if (.not.lrest) lrest = (index(keyword,' rest').ne.0.or.index(keyword,'rest').eq.1)
  if (.not.lrigid) lrigid = (index(keyword,' rigi').ne.0.or.index(keyword,'rigi').eq.1)
  if (.not.lrfo) lrfo = (index(keyword,' rfo').ne.0.or.index(keyword,'rfo').eq.1.or.index(keyword,'tran').ne.0)
  if (.not.lSandM) lSandM = (index(keyword,' sm').ne.0.or.index(keyword,'sm ').eq.1)
  if (.not.lsave) lsave = (index(keyword,' save').ne.0.or.index(keyword,'save').eq.1)
  if (.not.lshello) lshello = (index(keyword,' shel').ne.0.or.index(keyword,'shel').eq.1)
  if (.not.lsiteenergy) lsiteenergy = (index(keyword,' site').ne.0.or.index(keyword,'site').eq.1)
  if (.not.lsopt) lsopt = (index(keyword,' sopt').ne.0.or.index(keyword,'sopt').eq.1)
  if (.not.lspatial) lspatial = (index(keyword,' spat').ne.0.or.index(keyword,'spat').eq.1)
  if (.not.lspme) lspme = (index(keyword,' spme').ne.0.or.index(keyword,'spme').eq.1)
  if (.not.lstaticfirst) lstaticfirst = (index(keyword,' stat').ne.0.or.index(keyword,'stat').eq.1)
  if (.not.lstraincell) lstraincell = (index(keyword,' stra').ne.0.or.index(keyword,'stra').eq.1)
  if (.not.lstraincellprop) lstraincellprop = (index(keyword,' fspr').ne.0.or.index(keyword,'fspr').eq.1)
  if (.not.lstressout) lstressout = (index(keyword,' stre').ne.0.or.index(keyword,'stre').eq.1)
  if (.not.lsymoff) lsymoff = (index(keyword,' symoff').ne.0.or.index(keyword,'symoff').eq.1)
  if (.not.lthermal) lthermal = (index(keyword,' ther').ne.0.or.index(keyword,'ther').eq.1)
  if (.not.ltors) ltors = (index(keyword,' tors').ne.0.or.index(keyword,'tors').eq.1)
  if (.not.ltran) ltran = (index(keyword,' tran').ne.0.or.index(keyword,'tran').eq.1)
  if (.not.lunit) lunit = (index(keyword,' unit').ne.0.or.index(keyword,'unit').eq.1)
  if (.not.lxray) lxray = (index(keyword,' xray').ne.0.or.index(keyword,'xray').eq.1)
  if (.not.lzsisa) lzsisa = (index(keyword,' zsis').ne.0.or.index(keyword,'zsis').eq.1)
!
!  Keywords with default = .true.
!
  if (ldsym) ldsym = (index(keyword,' nodsym').eq.0.and.index(keyword,'nodsym').ne.1)
  if (lfix1atom) lfix1atom = (index(keyword,' unfi').eq.0.and.index(keyword,'unfi').ne.1)
  if (lDoElectrostatics) lDoElectrostatics = (index(keyword,' noel').eq.0.and.index(keyword,'noel').ne.1)
  if (lmodco) lmodco = (index(keyword,' nomod').eq.0.and.index(keyword,'nomod').ne.1)
  if (lmodco) lmodco = (index(keyword,' nowr').eq.0.and.index(keyword,'nowr').ne.1)
  if (lmolcom_mass) lmolcom_mass = (index(keyword,' com_n').eq.0.and.index(keyword,'com_n').ne.1)
  if (lnebdoublynudged) lnebdoublynudged = (index(keyword,' nodn').eq.0.and.index(keyword,'nodn').ne.1)
  if (lphasecom) lphasecom = (index(keyword,' noco').eq.0.and.index(keyword,'noco').ne.1)
  if (lreaxFFQ) lreaxFFQ = (index(keyword,' norx').eq.0.and.index(keyword,'norx').ne.1)
  if (lSandMnozz) lSandMnozz = (index(keyword,' smzz').eq.0.and.index(keyword,'smzz').ne.1)
  if (lsegsmooth) lsegsmooth = (index(keyword,' nosmo').eq.0.and.index(keyword,'nosmo').ne.1)
  if (lShengBTE3rd) lShengBTE3rd = (index(keyword,' nod3').eq.0.and.index(keyword,'nod3').ne.1)
  if (lsym) lsym = (index(keyword,' nosy').eq.0.and.index(keyword,'nosy').ne.1)
  if (lsymdok) lsymdok = (index(keyword,' nosd').eq.0.and.index(keyword,'nosd').ne.1)
!
!  Neutron Keywords (ers29)
!
  if (.not.lmakeeigarray) lmakeeigarray = (index(keyword,' make').ne.0.or.index(keyword,'make').eq.1)
  if (.not.lcoreinfo) lcoreinfo = (index(keyword,' core').ne.0.or.index(keyword,'core').eq.1)
  if (.not.lnoksym) lnoksym = (index(keyword,' noks').ne.0.or.index(keyword,'noks').eq.1)
  if (.not.lpdf) lpdf = (index(keyword,' pdf').ne.0.or.index(keyword,'pdf').eq.1)
  if (.not.lfreqcut) lfreqcut = (index(keyword,' pdfc').ne.0.or.index(keyword,'pdfc').eq.1)
  if (.not.lkeepcut) lkeepcut = (index(keyword,' pdfk').ne.0.or.index(keyword,'pdfk').eq.1)
  if (.not.lnowidth) lnowidth = (index(keyword,' nowi').ne.0.or.index(keyword,'nowi').eq.1)
  if (.not.lpartial) lpartial = (index(keyword,' nopa').eq.0.or.index(keyword,'nopa').eq.1)
!*******************************
!  Handle keyword dependances  *
!*******************************
!
!  lflags depends on other keywords
!
  if (lnoflags) then
    lflags = .false.
  else
    if (lopt.or.lgrad.or.lharmrelax.or.lfit.or.lrfo.or.lmc.or.lmd.or.lneb) then
      if ((.not.lconp).and.(.not.lconv).and.(.not.lcello)) lflags = .true.
    endif
  endif
!
!  PDF dependances: (ers29)
!
  if (lkeepcut) lfreqcut = .true.
  if (lfreqcut.or.lkeepcut) lpdf = .true.
  if (lpdf) then
    lmakeeigarray = .true.
    lsym = .false.
    if (.not.(index(keyword,' full').ne.0.or.index(keyword,'full').eq.1)) then
      write(keyword,*) trim(adjustl(keyword)), " full"
    endif
  endif
  if (lcoreinfo) lphon = .true.
  if (lmakeeigarray) then
    lphon = .true.
    leigen = .true.
    lnoksym = .true.
    lpdfout = .true.
  endif
!
!  If optlower specified then lower is needed
!
  if (loptlower) llower = .true.
!
!  If CI-NEB is requested then turn on NEB too
!
  if (lnebclimbingimage) lneb = .true.
!
!  If defect calc, then we need properties
!
  if (ldefect.and..not.lrest) lprop = .true.
!
!  If prop calc, then strain derivatives are needed
!
  if (lprop) then
    lstr = .true.
  endif
!
!  If prop calc, set flag for Born effective charges to be true
!
  if (lprop) then
    lborn = .true.
  endif
!
!  If atomic stresses are requested then strain derivatives are needed
!
  if (latomicstress) then
    lstr = .true.
  endif
!
!  If symmetry is to be lowered from imaginary modes then we need nosym and phonon calc
!
  if (llower) then
    lsym = .false.
    lphon = .true.
  endif
!
!  If IR intensities or eigenvectors or MSD or Grueneisen parameters are requested, then must be a phonon calc
!
  if (linten.or.leigen.or.lmsd.or.lgrueneisen) lphon = .true.
!
!  If Raman susceptibility requested, then properties must be computed
!
  if (lraman) lprop = .true.
!
!  If this is a thermal conductivity calculation then it also must be a phonon calc
!
  if (lthermal) lphon = .true.
!
!  If transition state calc, then lopt and lrfo must be true
!
  if (ltran) then
    morder = 1
    lopt = .true.
    lrfo = .true.
  endif
!
!  Change default update for unit hessian
!
  if (lunit.and.lfit) then
    nfupdate = 100
  endif
!
!  If QEq or SM or Pacha, then leem must be true
!
  if (lqeq.or.lSandM.or.lqtpie.or.lpacha) leem = .true.
!
!  If variable charge then there are no third derivatives
!  but we do need charge second derivatives
!
  if (leem) then
    lnoanald3 = .true.
    lDoQDeriv1 = .false.
    lDoQDeriv2 = .true.
  endif
!
!  If split bond variable charges are being used with all bonds then second derivatives can't be used
!
  if (leembond.and.(lallbonds.or.literativeQ)) then
    lnoanald2 = .true.
    lDoQDeriv1 = .false.
    lDoQDeriv2 = .false.
  endif
!
!  If split bond variable charges are being used then turn off symmetry for derivatives
!
  if (leembond) then
    lsymdok = .false.
  endif
!
!  If bond order charge potentials or ReaxFF are present turn off symmetry for derivatives
!
  if (nboQ.gt.0.or.lreaxFF.or.lEDIP) then
    lsymdok = .false.
  endif
!
!  If COSMIC then set cosmo to be true
!
  if (lcosmic) lcosmo = .true.
!
!  If COSMO then there are no third derivatives and symmetry
!  must not be enabled for derivatives.
!
  if (lcosmo) then
    lnoanald3 = .true.
    lsymdok = .false.
  endif
!
!  If SPME is specified then only first derivatives are currently available
!
  if (lspme) then
    lnoanald2 = .true.
    lnoanald3 = .true.
    lsymdok = .false.
  endif
!
!  If rigid molecules are being used then compute derivatives without symmetry and trap third derivatives
!
  if (lrigid) then
    lnoanald3 = .true.
    lsymdok = .false.
!
!  Force restart files to use cell vectors to preserve Cartesian orientation for restarting
!
    lusevectors  = .true.
!
!  Use fixed molecule connectivity since molecules are rigid
!
    lmolfix = .true.
  endif
!
!  If lstraincell is specified then select finite strain derivatives by default and turn off symmetry
!
  if (lstraincell) then
    lfinitestrain = .true.
    lusevectors = .true.
    lsym = .false.
  endif
!
!  If lstraincellprop is specified then turn on property calculation 
!
  if (lstraincellprop) lprop = .true.
!
!  If option to store vectors has been set then reset maxat arrays
!
  if (lStoreVectors) then
    call changemaxat
  endif
#ifdef TRACE
  call trace_out('setkeywrd')
#endif
!
  return
  end
