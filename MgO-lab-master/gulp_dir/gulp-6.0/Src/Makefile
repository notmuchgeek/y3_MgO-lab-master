################################################################################
#
#  GULP6.0
#
# Makefile for GULP6.0
# Use ./mkgulp utility to compile : See ./mkgulp -help
#
#  Original makefile by J-R. Hill, adapted for later versions of GULP by JDG
#
#  MPI additions by Ivan Walton, HPCI Southampton
#
################################################################################

default: gulp

#
#  Combine definitions
#
DEFS   = $(UDEFS) $(MDEFS) $(PDEFS) $(KDEFS) $(FDEFS)
LIBS   = $(SLIBS) $(FLIBS) $(KLIBS) $(PLIBS) $(MLIBS)
CFLAGS = $(INC) $(COPT)
FFLAGS = $(INC) $(FOPT)

F90    = $(RUNF90) $(OPT) $(FFLAGS) $(KFLAGS) $(BAGGER) $(DEFS) -c
F90FOX = $(RUNF90) $(OPT) $(FFLAGS) $(BAGGER) $(DEFS) -DFC_HAVE_ABORT -DFC_HAVE_FLUSH -c
F90O1  = $(RUNF90) $(OPT1) $(FFLAGS) $(KFLAGS) $(BAGGER) $(DEFS) -c
F90O2  = $(RUNF90) $(OPT2) $(FFLAGS) $(KFLAGS) $(BAGGER) $(DEFS) -c
F90O0  = $(RUNF90) $(FFLAGS) $(KFLAGS) $(BAGGER) $(DEFS) -c
CC     = $(RUNCC)  $(CFLAGS) $(KFLAGS) $(BAGGER) $(DEFS) -c

COBJECTS= csignal.o
LOBJECTS= lapack.o blas.o
FOBJECTS= datatypes.o iomod.o modules.o modules_kim.o itcom.o reallocate.o pr_random.o \
	  m_pr.o m_simplex.o m_conjgr.o m_plumed.o m_pdfvariables.o m_pdfneutron.o m_pdf.o m_fft3D.o \
          m_bondorder.o m_edip.o m_reaxff.o m_brenner.o m_trace.o m_alamode.o m_strain.o \
	  gulpcml.o gulpcml_md.o gulpcml_phonon.o gulpcml_props.o addframe2arc.o addshell.o getangles.o \
          applyconstraints.o \
	  aocn2.o applyconcfg.o banner.o bcgsolve.o binary.o bond.o bonded.o baskes.o bonded2donor.o bonded3.o bondpot.o \
	  boword.o bondorder.o bondorderfc.o cellcheck.o bondordermd.o bondordersd2.o bondorder12.o bondorder12d.o \
	  brennerd.o bondorder1fc.o bonded3c.o boself.o boselfd.o boselfp.o boselffc.o \
	  boselftrial.o bondq.o bondqtrial.o bspline.o bondorderd.o brennerfcd.o bondorderfcd.o brenner12d.o \
          bondint.o brenner.o brennerfc.o brenner1fc.o brennermd.o brennersd2.o brenner12.o buck4.o cell.o cellfhr.o celltype.o centre.o \
          background.o borncharge.o calcpht.o celldrv.o channels.o changemaxat.o changemaxatot.o changemaxtotr1at.o \
          changemaxatompernode.o changemaxcellpernode.o changemaxregion.o changemaxbondvec.o changemaxmcswapspec.o changemaxmcswaps.o \
          changemaxnboa.o changemaxnbor.o changemaxbond.o changemaxcfg.o changemaxcon.o changemaxcontot.o changemaxpotsites.o changemaxconnectpermol.o \
          changemaxnbopot.o changemaxnboo.o changemaxnboq.o changemaxconnect.o changemaxneamfnnumeric.o changemaxnebreplicatot.o changemaxnboz.o \
          changemaxdcon.o changemaxdef.o changemaxdis.o changemaxd2cells.o changemaxfcon.o changemaxsymop.o changemaxccspec.o \
          changemaxatloc.o changemaxfgrad.o changemaxfit.o changemaxfor.o changemaxkpt.o changemaxtrialatom.o changemaxbondq.o \
          changemaxeembond.o changemaxfstrain.o changemaxkvec.o changemaxlib.o changemaxlist3.o changemaxlist4.o changemaxlist6.o \
          changemaxmany.o changemaxmol.o changemaxmoleqv.o changemaxnobo.o changemaxobs.o changemaxiltor.o changemaxfkpt.o changemaxskpt.o \
          changemaxone.o changemaxfreq.o changemaxobsmode.o changemaxtdfield.o changemaxd2cell.o changemaxbondtype.o \
          changemaxpot.o changemaxppt.o changemaxproj.o changemaxproji.o changemaxqatoms.o changemaxqatoms2.o changemaxqrange.o \
          changemaxpts.o changemaxr1at.o changemaxr2at.o changemaxspcell.o changemaxspec.o changemaxuffspec.o changemaxmolat.o \
          changemaxeamspec.o changemaxeamfnspec.o changemaxthb.o changemaxtitle.o changemaxvac.o changemaxint.o changemaxvar.o \
          changemaxword.o changemaxd1.o changemaxd2.o changemaxnboq0.o changemaxsix.o changemaxtotint.o changemaxtotvac.o \
          changemaxn3a.o changemaxd2q.o changemaxn3f.o changemaxeamden.o changemaxnspcellattot.o changemaxnspcellattotbo.o \
          changemaxgbcfg.o changemaxgcmcmol.o changemaxreaxffspec.o changemaxreaxffval3.o cart2frac.o closestimage.o changemaxreaxfffixqspec.o \
          changemaxatompernodebo.o changemaxcellpernodebo.o changemaxspcellbo.o changemaxplanepot.o changemaxvectorpair.o \
          changemaxndistancetotal.o changemaxextrapol.o checkdefshspec.o changemaxn3bondnono.o changemaxsa.o \
          changemaxhold.o changemaxqvector.o changemaxnqstepfit.o changemaxnwstepfit.o changemaxedipspec.o \
          changemaxtempramp.o changemaxn3ma.o changemaxn3mf.o calcf2D.o \
          copycfg.o chvo.o collect.o costfn.o cputime.o create.o cscheck.o defect.o defener.o \
          checkpot.o chemshellenergy.o connectwraptrial.o cutscheck.o checkquaternions.o \
          compressd1.o compressd2.o compressd3.o connectwrap.o connectunwrap.o datetime.o dchargesplit.o dchargesplitd.o \
          defequ.o defequ2a.o deffreq.o deffreqd.o deffun.o defoin.o defopt.o d2charge3p.o d2charge3fc.o density1fc.o \
          defvartocfg.o \
          defout.o defpot.o defsec.o defsym.o deftmat.o defword.o derfc.o destroy.o d1charge3.o d2charge3.o dcharged.o \
          density0.o density3.o density3s.o densitysd.o density2a.o density11s2.o density112.o diagonal_fc.o diagonal_fcd.o \
          density12a2.o density12as2.o dynamic_fc.o dynamicn_fc.o dynamicn_fcd.o dynamicn_defect.o dbcgsolve.o endstring.o \
          dcharge.o dcharges.o d1charge.o d1charged.o d1charges.o d2charge.o d2chargefc.o d2chargep.o d2charge2D.o d2chargep2D.o \
          diagdyn.o dipole3D.o dipole3Dd.o dipole3Dtrial.o dipolesd.o distance.o dpoccon.o dspecial.o dump34.o dlamch.o dumpbo.o \
          dumpdur.o dumpgen.o dumppot.o dutils.o dynamic.o dynamicn.o eamfnderv.o eamscalescrho.o electricfieldtrial.o eckart.o energyfd.o \
          eem.o eemd.o eemdm.o eemsplit.o eemsplitd.o electricfield.o einstein.o einsteintrial.o einsteinsd2.o eispack.o emfunc.o energy.o epot.o epsilon.o equpos.o exchng.o \
          eamrho.o edipmd.o edip1fc.o edipcn.o edip2body.o edip3body.o electricfieldpotl.o electricfieldparts.o einsteinpd.o \
          fefunct.o fefunctn.o fitfunct.o fenergy0.o fenergy0d.o fenergy3.o fenergy3d.o force.o forcetrial.o four0d3.o \
          fourfc3.o four3d3.o four3d3dV.o four3d3g.o four3d3dVg.o fourbody.o fourmc.o fourmd.o fourmds.o finddef.o finish.o fit.o \
          fitbfgs.o fitcon.o firstpass.o fitfun.o fitgrad.o fitmin.o fitword.o findfractionword.o formula.o fourooppd.o fournosd.o \
          fournos.o four12.o four12d.o fourlist.o fourfc.o fourfcd.o four1fc.o fourp.o fourpd.o foursd.o forcebycell.o fouroopfc3.o fourimpfc3.o \
          fouroop.o fouroopd.o fouroopfc.o fouroopfcd.o fouroopp.o fouroopmc.o fouroopmd.o fouroopmds.o fouroopsd.o fouroopsd2.o fouroop12.o fouroop12d.o \
          fourimp.o fourimpd.o fourimpfc.o fourimpfcd.o fourimpp.o fourimpmc.o fourimpmd.o fourimpmds.o fourimpsd.o fourimpsd2.o fourimp12.o fourimp12d.o \
          fouroop1fc.o fourimp1fc.o fourimppd.o foursd2.o ftow.o functn.o functneb.o functsync.o functnf.o functnff.o \
          funct.o gabcfg.o gacreate.o gacross.o gaexpd.o gafit.o gamaxmin.o gamuta.o gamma.o gammasm.o getvibmode.o getvibmoded.o \
          gaopt.o garepro.o gasort.o gasteiger.o gasubs.o gauss.o gaword.o genpot.o genpotd.o genpotsplitd.o genpotdm.o genword.o getkeywrd.o getcoordno.o \
          getbocharge.o getbocharged.o getbochargetrial.o getbonds.o getderv1.o getdipole.o gtheta.o gthetabo.o gettraceinertia.o \
          getnmr.o getphi.o getneighbour.o getpotsymbol.o groupvelocities.o groupvelocitiesd.o getshellmatrices.o \
          getgcindex.o getstrain.o grueneisen.o grueneiseng.o grueneisend.o grueneisengd.o \
          harmonicrelax.o help.o hfunc.o lmatch.o lmatchany.o lmatchpair.o ind1toijk.o langevinpol.o \
          indsft.o indtoijk.o initial.o init_ins.o inword.o initarc.o initdervs.o initktrm.o init2bodydefaults.o \
          initqmatrix.o istoan.o itow.o kindex.o kvector.o kvectorf.o label.o labellc.o lambda.o init1bodydefaults.o \
          initmaxcfgdefaults.o  initmaxatdefaults.o initmaxatotdefaults.o initmaxbonddefaults.o improper.o \
          initmaxbondqdefaults.o initmaxccspecdefaults.o initmaxconnectdefaults.o initmaxdefdefaults.o \
          initmaxeamdendefaults.o initmaxeamfnspecdefaults.o initmaxeamspecdefaults.o initmaxfgraddefaults.o initmaxnboodefaults.o \
          initmaxfitdefaults.o initmaxfordefaults.o initmaxgcmcmoldefaults.o initmaxnbozdefaults.o initmaxmcswaps.o \
          initmaxlibdefaults.o initmaxmcswapspec.o initmaxmoldefaults.o initmaxnboadefaults.o initmaxnbopotdefaults.o \
          initmaxnboq0defaults.o initmaxnboqdefaults.o initmaxnbordefaults.o initmaxnebreplicatotdefaults.o \
          initmaxnppadefaults.o initmaxnptsdefaults.o initmaxobsdefaults.o initmaxplanepotdefaults.o \
          initmaxr1atdefaults.o initmaxreaxffspecdefaults.o initmaxreaxffval3defaults.o initmaxregiondefaults.o \
          initmaxsixdefaults.o initmaxspcellbodefaults.o initmaxspcelldefaults.o initmaxspecdefaults.o \
          initmaxqrangedefaults.o initmaxthbdefaults.o initmaxtitledefaults.o initmaxedipspecdefaults.o insertatom.o \
          initmaxfstraindefaults.o latpot.o lmbfgs.o linepro.o linpack.o lintoijk.o local.o lower.o \
          ltont.o cmatinv.o matinv.o imaxall.o maxall.o mdcellind.o mdcorrect.o mdfunct.o mdpredict.o kinthewedge.o \
          many.o manysd.o many0d.o many0dd.o many12.o many12d.o manyfc.o manyfcd.o manyp.o manypd.o mdke.o mdkeatomicstress.o mdshscale.o \
          manyd.o manysd2.o manymd0.o manymd3.o manymd3s.o manymi3.o many1fc.o mcstorecfg.o mcenergy.o \
          mccreate.o mctrialfunction.o mcdestroy.o mcmolinsert.o mcmove.o mcrotate.o mcstrain.o mcswap.o mddistconpost.o \
          mdgcorrect.o mdvvcorrect.o mdlfvcorrect.o mdkestrain.o mdtemperature.o \
          mdprop.o mdreset.o mdrestart.o mdscale.o mdstop.o mdvelcor.o \
          mcword.o mdword.o mdwrite.o methodok.o mindtoijk.o minimise.o molind.o molindfix.o meamsymdensity.o meamscreen.o \
          meamfnderv.o meamrho.o meamtotalrho.o meamtotalrhoderv.o meamtotalrhocrossderv.o molindtrial.o molindfixtrial.o \
          meamtotalrhoscreenderv.o moltype.o move2a1.o mxmb.o nagamma.o nearest.o nebforce.o nebword.o nrhess.o nrhessn.o nrstep.o numhess.o \
          nagammac.o nagammarigid.o nagammacrigid.o nudge.o okspec.o optim.o optin.o options.o optout.o optshe.o outerror.o outadvice.o \
          omegaproperty.o oscillatorstrength.o oscillatorstrengthg.o oscillatorstrengthc.o olinmin.o outcell.o outdcdheader.o \
          outarc.o outdefin.o outdener.o outderv.o outdlv.o outdrv.o outed.o outener.o outfdf.o outfile.o \
          outfrc.o outfit.o outfitcon.o outfreq.o outgafit.o outgaopt.o outgapar.o outneb.o outnebin.o outeem.o outoptvar.o \
          outsync.o outsyncin.o outspec_pdf.o onebody.o outlammps.o outlammpspots.o outatomicstress.o outinertia.o outgen.o outkey.o \
          outkpt.o outkptfile.o outmarv.o outmarv2.o outmc.o outmd.o outmdin.o outmol.o outshengBTE.o outfour.o outobs.o outphon.o \
          outpolar.o outpot.o outsiteenergy.o outspec.o outthree.o outtwo.o outthb.o outtime.o outbiograf.o outqbo.o phonon_fc.o \
          outvar.o outvarg.o outwarning.o outxr.o outxsf.o outxtl.o outxyz.o peigen.o peigend.o phonon.o phonond.o phonword.o poccon.o poccons.o \
          outcif.o outcssr.o outofp.o pireal0d.o pirealrecip.o pirealrecips.o outstructure.o outstresses.o phase_fc.o phase_fcd.o phonon_fcd.o \
          pdiagg.o pdiaggd.o pdiagd.o pdiage.o pdiagl.o psumall.o polarisation.o potgrid.o phoncopy1.o phoncopy2.o potwordm.o phase_fcg.o \
          potcut.o potential.o potword21.o potword3.o potword4.o potword6.o pressure.o planepot2D.o planepotfd.o potwordk.o phonon_fcg.o projd3dV.o projd3dVd.o projd3dVg.o projd3dVgd.o \
          potword22.o peigeng.o peigengd.o projd0.o projd0d.o projd0r.o projd3a.o projd3ad.o projthbk3ds.o projthbk3dV.o \
          projthbk3dVg.o projthbk3ddV.o projfork0.o projfork0r.o planepotmc.o planepotmd.o \
	  pirealrecipd.o pireal0dd.o pirealrecipp.o pirealrecippd.o qmatrixelemente.o \
          predict.o projthbk0.o projthbk0r.o projthbk3.o projthbk3dm.o projthbk3s.o projthbk3d.o pr_integrator.o psiscreenderv.o psiscreendervs.o \
          projfork3.o projfork3dm.o projfork3s.o projfork3d.o projfork3dd.o projfork3ds.o projfork3dV.o projfork3dVg.o \
          projfork3ddV.o property3.o property12.o property0.o psibaskes.o putpar.o pythag.o \
          psirho.o qeqbody.o qmatrix1D.o qmatrixelement.o qupdate.o radialforce.o radialforcetrial.o real1Dmc.o real1Dd.o real12ad.o \
          random.o ratfn.o real0d.o real0dd.o real1D.o real1Dmd.o real1Dp.o real1Dpd.o real1Dq.o real11.o real11d.o real11pot.o real11s.o real12a.o \
          real12apot.o real12as.o real2a.o real2ad.o real2a1.o real2a1d.o real2a1f.o real2a1fd.o real2adef.o realfc.o realfcd.o real1fc.o realselffc.o realselffcd.o \
          real2amany.o real12a2.o real12a2d.o real12as2.o real112.o real112d.o real11s2.o realmc0.o ramanstrength.o realmd0s.o realrecipfc3.o \
          real2aef.o real2b.o real4.o realcmm.o realmc3.o realmd3.o realmd3s.o realmd3t.o realmd0.o reale.o realp.o realpd.o realrecip3d3d.o \
          realrecip3d3dV.o realrecip3d3dVd.o realrecip3d3dVg.o realrecip3d3dVgd.o realfd3.o realfd0.o realmi3fd.o real1Dfd.o recip2Dfd.o recip3Dfd.o \
          fourfd.o fouroopfd.o fourimpfd.o sixfd.o threefd.o resetrigidmol.o \
          realmi3.o realmi3d.o realmimc3.o realrecip3d3.o real0d3.o real0d3d.o real1fc3.o real1d3.o real1d3dVg.o real1d3g.o raman0.o raman3.o ramantrm.o \
          real1d3dV.o reciptrmd3.o reciptrmd3g.o reciptrmdp.o reciptrmd3dV.o reciptrmd3dVg.o \
          realsd.o realsd2.o recip2D.o recip2Dtrial.o recip3D.o recip3Dmd.o recip3Dtrial.o recip3Dspme.o \
          realed.o recip3Dd.o recip3Dpd.o recip12a.o recip12ad.o recip12apot.o recip12as.o recip2Dd.o recip2Dpd.o \
          real2Dq.o recip2Dq.o reaxff1fc.o reaxffmd.o reaxffcharges.o reaxff_bo.o reaxff_funcs.o rebospline.o rfind.o rfindboth.o \
          recip4.o recip2Dp.o recip3Dp.o recipsd.o recipsd2.o recipfc3.o remove.o reorder.o resetcfgdefaults.o \
          rfohess.o rfostep.o rhoinv.o rhtohex.o rlist.o rpair.o rpair2.o runsync.o reinitialise.o resetx0rv.o \
          rigidmoleculedrv.o rigidmoleculephon.o rigidmoleculeprop.o rigidmoleculegv.o \
          reaxffpot.o rsearch1D.o rsearch2D.o rsearch3D.o rtlist.o rtclist.o runmc.o runmd.o runneb.o samemol.o \
          setcutoffmax.o shmin.o selftermfc.o selftermfcd.o shellgroupvelocities.o shellopt.o setphonptr.o setatomnoded2.o \
          setbelashchenko.o setcutoffmaxbo.o setspatialbo.o setatomnodesbo.o setspatialcells.o setedip.o seteem.o \
          sec0.o sec0d.o sec0f.o sec3.o sec3d.o sec3f.o selfterm.o set2a.o setbrenner1.o setbrenner3.o setc6.o setcfg.o setcmm.o \
          setcellderv.o seteam.o setratiom.o setreaxffQ.o setreaxffQ1fc.o setreaxffQiter.o setreaxffQiter1fc.o setregiontrans.o \
          seteembond.o setewald.o setdef.o setdisp.o setframe.o selftermp.o \
          setele.o setfour.o setkpt.o setsix.o setuff.o setlistimp.o setbondtype.o secondorderfc.o secondorderfcd.o \
          setatomnodes.o setbondorder.o setktrm4.o setktrmdp.o setktrmfc3.o setcombine.o setlib.o setlist3.o setlist4.o setlist6.o \
          setmaxcell1D.o setmd.o setmol.o setmols.o setmoldef.o setmoldef12a.o settorcubic.o settaper.o secondpass.o \
          setpote.o setreaxFF.o setslice.o setspec.o setspring.o settricubic.o settricubic2.o sumderv2f.o sumderv2p.o \
          setbsmptr.o setdr2mat.o setktrmd3.o setlistoop.o setoccptr.o setoccptrs.o setoccshptr.o setpolar.o setspatial.o \
          getdistances.o setdistance.o smbody.o setvbon.o symdervrot.o setoccshptrs.o symgenkpt.o setoptptr.o \
          setatomdistribution.o setvarnoded2.o setptsnoded2.o setrigidmol.o setmola.o setquaternion.o \
          setthree.o settings.o setup.o setzmolslice.o sort.o sort1.o sort1sym.o sort2a.o setstrainptr.o \
          setfeshtmat.o setquasiharm.o scaleevec.o setdestroyptr.o setmoveptr.o setrotateptr.o setkeywrd.o \
          setswapptr.o sitepot.o special.o specialcom.o specialqtn.o speclibw.o specword.o spline.o stolc.o stouc.o surfword.o sixfc3.o \
          sixbody.o sixfc.o sixfcd.o six1fc.o sixlist.o sixmc.o sixmd.o sixmds.o sixnos.o sixnosd.o sixp.o sixpd.o sixsd.o sixsd2.o \
          strain.o strfin.o strsym.o strword.o sumderv2.o sumderv2s.o surfaceenergy.o syncforce.o threefc.o threefcd.o three1fc.o \
          super.o symbind.o symderv2.o symmet.o symmetp.o symoff.o symout.o sympot.o sympotd.o symprod.o tfunc.o thermalconductivity.o thermalconductivityd.o \
          taper.o threenos.o threenosd.o three12.o three12d.o threelist.o threep.o threepd.o threemc.o threemd.o threemds.o threesd.o threefc3.o \
          threesd2.o threebody.o three0d3.o three3d3.o three3d3dV.o three3d3g.o three3d3dVg.o gettorsions.o transmat.o transmatd.o \
          twobody.o twobodymd.o twoden.o vartocfg.o tangtoennies.o transmatmol.o \
          uncell.o uncentre.o unfreeze.o validtriad.o volume.o wcoord.o wolfself.o wolfselftrial.o kweight.o worsy.o wtof.o \
          thirdordern.o thirdorderfc3.o vibenergy.o x0tostrtrial.o x0strain.o x0tostr.o x0tostrcentroid.o xctox0.o $(CDABS) \
	  getline.o initcomms.o sumall.o sendall.o stopnow.o mpfinish.o mpbarrier.o \
          initmemory.o nullpointer.o outofmemory.o setrotation.o \
          iterinv.o itpack.o findcelltype.o findspacegroup.o testspacegroup.o \
          gulp.o gulpmain.o gulpsetup.o gulpfinish.o \
          changemaxatcosmo.o cosmoaiiderv.o cosmoaijderv.o cosmoaijzero.o cosmoamatdadd.o cosmoderv.o cosmoderv2D.o \
          dcosmowt.o dtotcosmowt.o setcosmo.o setcosmoainv.o setcosmoainv2D.o setcosmoamat.o setcosmoamat2D.o setcosmobmat.o \
          setdcosmoamat.o setdcosmobmat.o changemaxnearseg.o changemaxnppa.o changemaxnpts.o \
          changemaxnptsh.o changemaxnptstot.o changemaxnpwt.o changemaxsasparticles.o cosmoamatdaddd.o \
          changemaxsasparticlespart.o setsegweight.o setdsegweight.o setqonsas.o setqonsas2D.o setsas.o \
          updatesas.o solvation.o switch.o cosmoaiizero.o changemaxnset.o createsphere.o cosmoamatdadd1p.o \
          initqmatrixc.o qmatrixelementc.o changemaxallnearseg.o setdcosmoamat1p.o setdcosmobmat1p.o \
          qtpie.o ctfunct.o qtmin.o screenct.o reaxffctfunct.o \
          kim_functions.o setkim.o kimmd.o \
	  trap_sigint.o 

CSDUM = CS_dummy.o

FOBJ = $(FOBJECTS) 
OBJ  = $(FOBJ) $(COBJECTS) 

gulp: $(OBJ) $(CSDUM) $(LOBJECTS) ; $(RUNF90) $(LDFLAGS) $(BAGGER) $(OBJ) $(CSDUM) $(LIBS) $(FoX_libs) -o gulp

lib:    $(OBJ)
	$(ARCHIVE) libgulp.a $(OBJ)
	$(RANLIB) libgulp.a

clean:           ; rm -f $(OBJ) $(CSDUM) *.mod ; rm -f ../gulp
libclean:        ; rm -f $(OBJ) *.mod ; rm -f ../../libgulp.a
html:            ; ruby ../Docs/gulp2html.rb ../Docs/help.txt > ../Docs/help.html
tar  :           ; cd .. ; tar -cvf gulp-6.0.tar README Src/*.F90 Src/*.c \
                   Src/Makefile Examples/README Examples/example*.gin Examples/*.dbout Libraries/eledata Src/mkgulp \
                   Libraries/*.lib \
                   Utils/MC2arc/*.f Utils/MC2arc/*.F Utils/MC2arc/README \
                   Utils/MC2arc/Makefile Utils/MC2arc/getmachine Utils/ReaxFF/reaxff2lib.f90 \
                   Utils/MDarc2his/arc2his.f90 Utils/MDarc2his/cell.f90 Utils/MDarc2his/README \
                   Utils/MDbin2ascii/README Utils/MDbin2ascii/mdbin2ascii.f90 \
                   Utils/Alamode/README Utils/Alamode/example67.gin Utils/Alamode/runalamode \
                   Utils/Morphology/README Utils/Morphology/emorph.f90 Utils/Morphology/matinv.f90 Utils/Morphology/nacl.emf \
                   Docs/gulp.intro Docs/gulp2html.rb Docs/help.txt Docs/gulp6.0_manual.pdf \
                   Src/FoX-4.1.2/LICENSE Src/FoX-4.1.2/configure Src/FoX-4.1.2/arch.make.in Src/FoX-4.1.2/Makefile Src/FoX-4.1.2/README.FoX.txt \
                   Src/FoX-4.1.2/FoX-config.in Src/FoX-4.1.2/config/config.* Src/FoX-4.1.2/config/install-sh Src/FoX-4.1.2/utils/*.f90 \
                   Src/FoX-4.1.2/utils/*.F90 Src/FoX-4.1.2/utils/makefile \
                   Src/FoX-4.1.2/fsys/*.f90 Src/FoX-4.1.2/fsys/*.F90 Src/FoX-4.1.2/fsys/makefile Src/FoX-4.1.2/common/*.F90 Src/FoX-4.1.2/common/makefile \
                   Src/FoX-4.1.2/wxml/*.f90 Src/FoX-4.1.2/wxml/*.F90 Src/FoX-4.1.2/wxml/makefile Src/FoX-4.1.2/wcml/*.f90 Src/FoX-4.1.2/wcml/*.F90 \
                   Src/FoX-4.1.2/wcml/makefile; \
                   gzip -9 gulp-6.0.tar

## FoX config and making
FOX_PATH = ../FoX-4.1.2
FOX_CONFIG = $(FOX_PATH)/FoX-config
FOX_OBJ = $(FOX_PATH)/wcml/FoX_wcml.o

$(FOX_CONFIG):
ifeq (,$(findstring DNOFOX,$(DEFS)))
	(cd $(FOX_PATH); ./configure --enable-wcml --enable-wxml FC='$(RUNF90)' FCFLAGS='$(OPT) $(FFLAGS) $(BAGGER)' LDFLAGS='$(LDFLAGS) $(BAGGER) $(LIBS)')
else
	@echo 'Building with FoX disabled'
endif

$(FOX_OBJ): $(FOX_CONFIG)
ifeq (,$(findstring DNOFOX,$(DEFS)))
	(cd $(FOX_PATH);make wxml_lib wcml_lib)
FoX_fcflags=`$(FOX_PATH)/FoX-config --fcflags --wcml --wxml`
FoX_libs=`$(FOX_PATH)/FoX-config --libs --wcml --wxml`
else
FoX_fcflags=
FoX_libs=
endif

## End of FoX stuff

# .F90 files
banner.o                  : ../banner.F90                  ; $(F90) ../banner.F90
channels.o                : ../channels.F90                ; $(F90) ../channels.F90
changemaxat.o             : ../changemaxat.F90             ; $(F90) ../changemaxat.F90
changemaxcfg.o            : ../changemaxcfg.F90            ; $(F90) ../changemaxcfg.F90
chemshellenergy.o         : ../chemshellenergy.F90         ; $(F90) ../chemshellenergy.F90
cputime.o                 : ../cputime.F90                 ; $(F90) ../cputime.F90
energy.o                  : ../energy.F90                  ; $(F90) ../energy.F90
energyfd.o                : ../energyfd.F90                ; $(F90) ../energyfd.F90
getline.o                 : ../getline.F90                 ; $(F90) ../getline.F90
datatypes.o               : ../datatypes.F90               ; $(F90) ../datatypes.F90
dcharged.o                : ../dcharged.F90                ; $(F90) ../dcharged.F90
dchargesplit.o            : ../dchargesplit.F90            ; $(F90) ../dchargesplit.F90
dchargesplitd.o           : ../dchargesplitd.F90           ; $(F90) ../dchargesplitd.F90
deffreqd.o                : ../deffreqd.F90                ; $(F90) ../deffreqd.F90
defsec.o                  : ../defsec.F90                  ; $(F90) ../defsec.F90
defvartocfg.o             : ../defvartocfg.F90             ; $(F90) ../defvartocfg.F90
eckart.o                  : ../eckart.F90                  ; $(F90) ../eckart.F90
fenergy0d.o               : ../fenergy0d.F90               ; $(F90) ../fenergy0d.F90
fenergy3d.o               : ../fenergy3d.F90               ; $(F90) ../fenergy3d.F90
finish.o                  : ../finish.F90                  ; $(F90) ../finish.F90
fitbfgs.o                 : ../fitbfgs.F90                 ; $(F90) ../fitbfgs.F90 
groupvelocitiesd.o        : ../groupvelocitiesd.F90        ; $(F90) ../groupvelocitiesd.F90
gulp.o                    : ../gulp.F90                    ; $(F90) ../gulp.F90
gulpmain.o                : ../gulpmain.F90                ; $(F90) $(FoX_fcflags) ../gulpmain.F90
gulpsetup.o               : ../gulpsetup.F90               ; $(F90) $(FoX_fcflags) ../gulpsetup.F90
harmonicrelax.o           : ../harmonicrelax.F90           ; $(F90) ../harmonicrelax.F90
initcomms.o               : ../initcomms.F90               ; $(F90) ../initcomms.F90
initdervs.o               : ../initdervs.F90               ; $(F90) ../initdervs.F90
initial.o                 : ../initial.F90                 ; $(F90) ../initial.F90
initmaxcfgdefaults.o      : ../initmaxcfgdefaults.F90      ; $(F90) ../initmaxcfgdefaults.F90
kim_functions.o           : ../kim_functions.F90           ; $(F90) ../kim_functions.F90
kimmd.o                   : ../kimmd.F90                   ; $(F90) ../kimmd.F90
local.o                   : ../local.F90                   ; $(F90) ../local.F90
m_fft3D.o                 : ../m_fft3D.F90                 ; $(F90) ../m_fft3D.F90
m_plumed.o                : ../m_plumed.F90                ; $(F90) ../m_plumed.F90
m_simplex.o               : ../m_simplex.F90               ; $(F90) ../m_simplex.F90 
m_strain.o                : ../m_strain.F90                ; $(F90) ../m_strain.F90 
cmatinv.o                 : ../cmatinv.F90                 ; $(F90) ../cmatinv.F90
matinv.o                  : ../matinv.F90                  ; $(F90) ../matinv.F90
imaxall.o                 : ../imaxall.F90                 ; $(F90) ../imaxall.F90 
maxall.o                  : ../maxall.F90                  ; $(F90) ../maxall.F90 
mcenergy.o                : ../mcenergy.F90                ; $(F90) ../mcenergy.F90
mdfunct.o                 : ../mdfunct.F90                 ; $(F90) ../mdfunct.F90
minimise.o                : ../minimise.F90                ; $(F90) $(FoX_fcflags) ../minimise.F90
modules_kim.o             : ../modules_kim.F90             ; $(F90) ../modules_kim.F90
mpbarrier.o               : ../mpbarrier.F90               ; $(F90) ../mpbarrier.F90
mpfinish.o                : ../mpfinish.F90                ; $(F90) ../mpfinish.F90
nrhess.o                  : ../nrhess.F90                  ; $(F90) ../nrhess.F90
nrhessn.o                 : ../nrhessn.F90                 ; $(F90) ../nrhessn.F90
nrstep.o                  : ../nrstep.F90                  ; $(F90) ../nrstep.F90
nullpointer.o             : ../nullpointer.F90             ; $(F90) ../nullpointer.F90
optim.o                   : ../optim.F90                   ; $(F90) ../optim.F90
pdiagd.o                  : ../pdiagd.F90                  ; $(F90) ../pdiagd.F90
pdiaggd.o                 : ../pdiaggd.F90                 ; $(F90) ../pdiaggd.F90
phonond.o                 : ../phonond.F90                 ; $(F90) ../phonond.F90
phonon_fcd.o              : ../phonon_fcd.F90              ; $(F90) ../phonon_fcd.F90
rfostep.o                 : ../rfostep.F90                 ; $(F90) ../rfostep.F90
runmd.o                   : ../runmd.F90                   ; $(F90) $(FoX_fcflags) ../runmd.F90
sec0d.o                   : ../sec0d.F90                   ; $(F90) ../sec0d.F90
sec3d.o                   : ../sec3d.F90                   ; $(F90) ../sec3d.F90
secondorderfcd.o          : ../secondorderfcd.F90          ; $(F90) ../secondorderfcd.F90
sendall.o                 : ../sendall.F90                 ; $(F90) ../sendall.F90
setcosmoainv2D.o          : ../setcosmoainv2D.F90          ; $(F90) ../setcosmoainv2D.F90
setcutoffmax.o            : ../setcutoffmax.F90            ; $(F90) ../setcutoffmax.F90
setkim.o                  : ../setkim.F90                  ; $(F90) ../setkim.F90
setqonsas2D.o             : ../setqonsas2D.F90             ; $(F90) ../setqonsas2D.F90
setquaternion.o           : ../setquaternion.F90           ; $(F90) ../setquaternion.F90
setquasiharm.o            : ../setquasiharm.F90            ; $(F90) ../setquasiharm.F90
setrotation.o             : ../setrotation.F90             ; $(F90) ../setrotation.F90
stopnow.o                 : ../stopnow.F90                 ; $(F90) ../stopnow.F90
sumall.o                  : ../sumall.F90                  ; $(F90) ../sumall.F90
thermalconductivityd.o    : ../thermalconductivityd.F90    ; $(F90) ../thermalconductivityd.F90
CS_dummy.o                : ../CS_dummy.F90                ; $(F90) ../CS_dummy.F90
# Former .f90 files
modules.o                 : ../modules.F90                 ; $(F90) ../modules.F90
iomod.o                   : ../iomod.F90                   ; $(F90) ../iomod.F90
pr_random.o               : ../pr_random.F90               ; $(F90) ../pr_random.F90
m_alamode.o               : ../m_alamode.F90               ; $(F90) ../m_alamode.F90
m_bondorder.o             : ../m_bondorder.F90             ; $(F90) ../m_bondorder.F90
m_brenner.o               : ../m_brenner.F90               ; $(F90) ../m_brenner.F90
m_edip.o                  : ../m_edip.F90                  ; $(F90) ../m_edip.F90
m_reaxff.o                : ../m_reaxff.F90                ; $(F90) ../m_reaxff.F90
m_pr.o                    : ../m_pr.F90                    ; $(F90) ../m_pr.F90
m_trace.o                 : ../m_trace.F90                 ; $(F90) ../m_trace.F90
addframe2arc.o            : ../addframe2arc.F90            ; $(F90) ../addframe2arc.F90
addshell.o                : ../addshell.F90                ; $(F90) ../addshell.F90
aocn2.o                   : ../aocn2.F90                   ; $(F90) ../aocn2.F90
applyconcfg.o             : ../applyconcfg.F90             ; $(F90) ../applyconcfg.F90
applyconstraints.o        : ../applyconstraints.F90        ; $(F90) ../applyconstraints.F90
background.o              : ../background.F90              ; $(F90) ../background.F90
baskes.o                  : ../baskes.F90                  ; $(F90) ../baskes.F90
bcgsolve.o                : ../bcgsolve.F90                ; $(F90) ../bcgsolve.F90
binary.o                  : ../binary.F90                  ; $(F90) ../binary.F90
blas.o                    : ../blas.F90                    ; $(F90) ../blas.F90
bond.o                    : ../bond.F90                    ; $(F90) ../bond.F90
bonded.o                  : ../bonded.F90                  ; $(F90) ../bonded.F90
bonded2donor.o            : ../bonded2donor.F90            ; $(F90) ../bonded2donor.F90
bonded3.o                 : ../bonded3.F90                 ; $(F90) ../bonded3.F90
bonded3c.o                : ../bonded3c.F90                ; $(F90) ../bonded3c.F90
bondint.o                 : ../bondint.F90                 ; $(F90) ../bondint.F90
bondorder.o               : ../bondorder.F90               ; $(F90) ../bondorder.F90
bondorder12.o             : ../bondorder12.F90             ; $(F90) ../bondorder12.F90
bondorder12d.o            : ../bondorder12d.F90            ; $(F90) ../bondorder12d.F90
bondorder1fc.o            : ../bondorder1fc.F90            ; $(F90) ../bondorder1fc.F90
bondorderd.o              : ../bondorderd.F90              ; $(F90) ../bondorderd.F90
bondorderfc.o             : ../bondorderfc.F90             ; $(F90) ../bondorderfc.F90
bondorderfcd.o            : ../bondorderfcd.F90            ; $(F90) ../bondorderfcd.F90
bondordermd.o             : ../bondordermd.F90             ; $(F90) ../bondordermd.F90
bondordersd2.o            : ../bondordersd2.F90            ; $(F90) ../bondordersd2.F90
bondpot.o                 : ../bondpot.F90                 ; $(F90) ../bondpot.F90
bondq.o                   : ../bondq.F90                   ; $(F90) ../bondq.F90
bondqtrial.o              : ../bondqtrial.F90              ; $(F90) ../bondqtrial.F90
borncharge.o              : ../borncharge.F90              ; $(F90) ../borncharge.F90
boself.o                  : ../boself.F90                  ; $(F90) ../boself.F90
boselfd.o                 : ../boselfd.F90                 ; $(F90) ../boselfd.F90
boselfp.o                 : ../boselfp.F90                 ; $(F90) ../boselfp.F90
boselfpd.o                : ../boselfpd.F90                ; $(F90) ../boselfpd.F90
boselffc.o                : ../boselffc.F90                ; $(F90) ../boselffc.F90
boselftrial.o             : ../boselftrial.F90             ; $(F90) ../boselftrial.F90
boword.o                  : ../boword.F90                  ; $(F90) ../boword.F90
brenner.o                 : ../brenner.F90                 ; $(F90) ../brenner.F90
brennerd.o                : ../brennerd.F90                ; $(F90) ../brennerd.F90
brenner12.o               : ../brenner12.F90               ; $(F90) ../brenner12.F90
brenner12d.o              : ../brenner12d.F90              ; $(F90) ../brenner12d.F90
brenner1fc.o              : ../brenner1fc.F90              ; $(F90) ../brenner1fc.F90
brennerfc.o               : ../brennerfc.F90               ; $(F90) ../brennerfc.F90
brennerfcd.o              : ../brennerfcd.F90              ; $(F90) ../brennerfcd.F90
brennermd.o               : ../brennermd.F90               ; $(F90) ../brennermd.F90
brennersd2.o              : ../brennersd2.F90              ; $(F90) ../brennersd2.F90
bspline.o                 : ../bspline.F90                 ; $(F90) ../bspline.F90
buck4.o                   : ../buck4.F90                   ; $(F90) ../buck4.F90
calcf2D.o                 : ../calcf2D.F90                 ; $(F90) ../calcf2D.F90
calcpht.o                 : ../calcpht.F90                 ; $(F90) ../calcpht.F90
cdabs.o                   :                                ; $(F90) ../cdabs.F90
cell.o                    : ../cell.F90                    ; $(F90) ../cell.F90
cellcheck.o               : ../cellcheck.F90               ; $(F90) ../cellcheck.F90
celldrv.o                 : ../celldrv.F90                 ; $(F90) ../celldrv.F90
cellfhr.o                 : ../cellfhr.F90                 ; $(F90) ../cellfhr.F90
celltype.o                : ../celltype.F90                ; $(F90) ../celltype.F90
centre.o                  : ../centre.F90                  ; $(F90) ../centre.F90
cart2frac.o               : ../cart2frac.F90               ; $(F90) ../cart2frac.F90
changemaxatcosmo.o        : ../changemaxatcosmo.F90        ; $(F90) ../changemaxatcosmo.F90
changemaxatloc.o          : ../changemaxatloc.F90          ; $(F90) ../changemaxatloc.F90
changemaxatot.o           : ../changemaxatot.F90           ; $(F90) ../changemaxatot.F90
changemaxatompernode.o    : ../changemaxatompernode.F90    ; $(F90) ../changemaxatompernode.F90
changemaxatompernodebo.o  : ../changemaxatompernodebo.F90  ; $(F90) ../changemaxatompernodebo.F90
changemaxbond.o           : ../changemaxbond.F90           ; $(F90) ../changemaxbond.F90
changemaxbondq.o          : ../changemaxbondq.F90          ; $(F90) ../changemaxbondq.F90
changemaxbondtype.o       : ../changemaxbondtype.F90       ; $(F90) ../changemaxbondtype.F90
changemaxbondvec.o        : ../changemaxbondvec.F90        ; $(F90) ../changemaxbondvec.F90
changemaxccspec.o         : ../changemaxccspec.F90         ; $(F90) ../changemaxccspec.F90
changemaxcellpernode.o    : ../changemaxcellpernode.F90    ; $(F90) ../changemaxcellpernode.F90
changemaxcellpernodebo.o  : ../changemaxcellpernodebo.F90  ; $(F90) ../changemaxcellpernodebo.F90
changemaxcon.o            : ../changemaxcon.F90            ; $(F90) ../changemaxcon.F90
changemaxconnect.o        : ../changemaxconnect.F90        ; $(F90) ../changemaxconnect.F90
changemaxconnectpermol.o  : ../changemaxconnectpermol.F90  ; $(F90) ../changemaxconnectpermol.F90
changemaxcontot.o         : ../changemaxcontot.F90         ; $(F90) ../changemaxcontot.F90
changemaxd1.o             : ../changemaxd1.F90             ; $(F90) ../changemaxd1.F90
changemaxd2.o             : ../changemaxd2.F90             ; $(F90) ../changemaxd2.F90
changemaxd2cell.o         : ../changemaxd2cell.F90         ; $(F90) ../changemaxd2cell.F90
changemaxd2cells.o        : ../changemaxd2cells.F90        ; $(F90) ../changemaxd2cells.F90
changemaxd2q.o            : ../changemaxd2q.F90            ; $(F90) ../changemaxd2q.F90
changemaxdcon.o           : ../changemaxdcon.F90           ; $(F90) ../changemaxdcon.F90
changemaxdef.o            : ../changemaxdef.F90            ; $(F90) ../changemaxdef.F90
changemaxdis.o            : ../changemaxdis.F90            ; $(F90) ../changemaxdis.F90
changemaxeamden.o         : ../changemaxeamden.F90         ; $(F90) ../changemaxeamden.F90
changemaxeamfnspec.o      : ../changemaxeamfnspec.F90      ; $(F90O1) ../changemaxeamfnspec.F90
changemaxeamspec.o        : ../changemaxeamspec.F90        ; $(F90) ../changemaxeamspec.F90
changemaxedipspec.o       : ../changemaxedipspec.F90       ; $(F90) ../changemaxedipspec.F90
changemaxeembond.o        : ../changemaxeembond.F90        ; $(F90) ../changemaxeembond.F90
changemaxextrapol.o       : ../changemaxextrapol.F90       ; $(F90) ../changemaxextrapol.F90
changemaxfcon.o           : ../changemaxfcon.F90           ; $(F90) ../changemaxfcon.F90
changemaxfgrad.o          : ../changemaxfgrad.F90          ; $(F90) ../changemaxfgrad.F90
changemaxfstrain.o        : ../changemaxfstrain.F90        ; $(F90) ../changemaxfstrain.F90
changemaxfit.o            : ../changemaxfit.F90            ; $(F90) ../changemaxfit.F90
changemaxfor.o            : ../changemaxfor.F90            ; $(F90) ../changemaxfor.F90
changemaxfreq.o           : ../changemaxfreq.F90           ; $(F90) ../changemaxfreq.F90
changemaxgbcfg.o          : ../changemaxgbcfg.F90          ; $(F90) ../changemaxgbcfg.F90
changemaxgcmcmol.o        : ../changemaxgcmcmol.F90        ; $(F90) ../changemaxgcmcmol.F90
changemaxhold.o           : ../changemaxhold.F90           ; $(F90) ../changemaxhold.F90
changemaxiltor.o          : ../changemaxiltor.F90          ; $(F90) ../changemaxiltor.F90
changemaxkpt.o            : ../changemaxkpt.F90            ; $(F90) ../changemaxkpt.F90
changemaxfkpt.o           : ../changemaxfkpt.F90           ; $(F90) ../changemaxfkpt.F90
changemaxskpt.o           : ../changemaxskpt.F90           ; $(F90) ../changemaxskpt.F90
changemaxkvec.o           : ../changemaxkvec.F90           ; $(F90) ../changemaxkvec.F90
changemaxlib.o            : ../changemaxlib.F90            ; $(F90) ../changemaxlib.F90
changemaxlist3.o          : ../changemaxlist3.F90          ; $(F90) ../changemaxlist3.F90
changemaxlist4.o          : ../changemaxlist4.F90          ; $(F90) ../changemaxlist4.F90
changemaxlist6.o          : ../changemaxlist6.F90          ; $(F90) ../changemaxlist6.F90
changemaxmany.o           : ../changemaxmany.F90           ; $(F90) ../changemaxmany.F90
changemaxmcswaps.o        : ../changemaxmcswaps.F90        ; $(F90) ../changemaxmcswaps.F90
changemaxmcswapspec.o     : ../changemaxmcswapspec.F90     ; $(F90) ../changemaxmcswapspec.F90
changemaxmol.o            : ../changemaxmol.F90            ; $(F90) ../changemaxmol.F90
changemaxmolat.o          : ../changemaxmolat.F90          ; $(F90) ../changemaxmolat.F90
changemaxmoleqv.o         : ../changemaxmoleqv.F90         ; $(F90) ../changemaxmoleqv.F90
changemaxn3a.o            : ../changemaxn3a.F90            ; $(F90) ../changemaxn3a.F90
changemaxn3f.o            : ../changemaxn3f.F90            ; $(F90) ../changemaxn3f.F90
changemaxn3ma.o           : ../changemaxn3ma.F90           ; $(F90) ../changemaxn3ma.F90
changemaxn3mf.o           : ../changemaxn3mf.F90           ; $(F90) ../changemaxn3mf.F90
changemaxn3bondnono.o     : ../changemaxn3bondnono.F90     ; $(F90) ../changemaxn3bondnono.F90
changemaxallnearseg.o     : ../changemaxallnearseg.F90     ; $(F90) ../changemaxallnearseg.F90
changemaxnbopot.o         : ../changemaxnbopot.F90         ; $(F90) ../changemaxnbopot.F90
changemaxnboa.o           : ../changemaxnboa.F90           ; $(F90) ../changemaxnboa.F90
changemaxnboo.o           : ../changemaxnboo.F90           ; $(F90) ../changemaxnboo.F90
changemaxnboq.o           : ../changemaxnboq.F90           ; $(F90) ../changemaxnboq.F90
changemaxnboq0.o          : ../changemaxnboq0.F90          ; $(F90) ../changemaxnboq0.F90
changemaxnbor.o           : ../changemaxnbor.F90           ; $(F90) ../changemaxnbor.F90
changemaxnboz.o           : ../changemaxnboz.F90           ; $(F90) ../changemaxnboz.F90
changemaxndistancetotal.o : ../changemaxndistancetotal.F90 ; $(F90) ../changemaxndistancetotal.F90
changemaxneamfnnumeric.o  : ../changemaxneamfnnumeric.F90  ; $(F90) ../changemaxneamfnnumeric.F90
changemaxnearseg.o        : ../changemaxnearseg.F90        ; $(F90) ../changemaxnearseg.F90
changemaxnebreplicatot.o  : ../changemaxnebreplicatot.F90  ; $(F90) ../changemaxnebreplicatot.F90
changemaxnset.o           : ../changemaxnset.F90           ; $(F90) ../changemaxnset.F90
changemaxnspcellattot.o   : ../changemaxnspcellattot.F90   ; $(F90) ../changemaxnspcellattot.F90
changemaxnspcellattotbo.o : ../changemaxnspcellattotbo.F90 ; $(F90) ../changemaxnspcellattotbo.F90
changemaxnobo.o           : ../changemaxnobo.F90           ; $(F90) ../changemaxnobo.F90
changemaxnppa.o           : ../changemaxnppa.F90           ; $(F90) ../changemaxnppa.F90
changemaxnpwt.o           : ../changemaxnpwt.F90           ; $(F90) ../changemaxnpwt.F90
changemaxnpts.o           : ../changemaxnpts.F90           ; $(F90) ../changemaxnpts.F90
changemaxnptsh.o          : ../changemaxnptsh.F90          ; $(F90) ../changemaxnptsh.F90
changemaxnptstot.o        : ../changemaxnptstot.F90        ; $(F90) ../changemaxnptstot.F90
changemaxnqstepfit.o      : ../changemaxnqstepfit.F90      ; $(F90) ../changemaxnqstepfit.F90
changemaxnwstepfit.o      : ../changemaxnwstepfit.F90      ; $(F90) ../changemaxnwstepfit.F90
changemaxobs.o            : ../changemaxobs.F90            ; $(F90) ../changemaxobs.F90
changemaxobsmode.o        : ../changemaxobsmode.F90        ; $(F90) ../changemaxobsmode.F90
changemaxone.o            : ../changemaxone.F90            ; $(F90) ../changemaxone.F90
changemaxplanepot.o       : ../changemaxplanepot.F90       ; $(F90) ../changemaxplanepot.F90
changemaxpot.o            : ../changemaxpot.F90            ; $(F90) ../changemaxpot.F90
changemaxpotsites.o       : ../changemaxpotsites.F90       ; $(F90) ../changemaxpotsites.F90
changemaxppt.o            : ../changemaxppt.F90            ; $(F90) ../changemaxppt.F90
changemaxproj.o           : ../changemaxproj.F90           ; $(F90) ../changemaxproj.F90
changemaxproji.o          : ../changemaxproji.F90          ; $(F90) ../changemaxproji.F90
changemaxpts.o            : ../changemaxpts.F90            ; $(F90) ../changemaxpts.F90
changemaxqrange.o         : ../changemaxqrange.F90         ; $(F90) ../changemaxqrange.F90
changemaxqatoms.o         : ../changemaxqatoms.F90         ; $(F90) ../changemaxqatoms.F90
changemaxqatoms2.o        : ../changemaxqatoms2.F90        ; $(F90) ../changemaxqatoms2.F90
changemaxqvector.o        : ../changemaxqvector.F90        ; $(F90) ../changemaxqvector.F90
changemaxr1at.o           : ../changemaxr1at.F90           ; $(F90) ../changemaxr1at.F90
changemaxr2at.o           : ../changemaxr2at.F90           ; $(F90) ../changemaxr2at.F90
changemaxreaxfffixqspec.o : ../changemaxreaxfffixqspec.F90 ; $(F90) ../changemaxreaxfffixqspec.F90
changemaxreaxffspec.o     : ../changemaxreaxffspec.F90     ; $(F90) ../changemaxreaxffspec.F90
changemaxreaxffval3.o     : ../changemaxreaxffval3.F90     ; $(F90) ../changemaxreaxffval3.F90
changemaxregion.o         : ../changemaxregion.F90         ; $(F90) ../changemaxregion.F90
changemaxsa.o             : ../changemaxsa.F90             ; $(F90) ../changemaxsa.F90
changemaxsasparticles.o   : ../changemaxsasparticles.F90   ; $(F90) ../changemaxsasparticles.F90
changemaxsasparticlespart.o : ../changemaxsasparticlespart.F90 ; $(F90) ../changemaxsasparticlespart.F90
changemaxsix.o            : ../changemaxsix.F90            ; $(F90) ../changemaxsix.F90
changemaxspcell.o         : ../changemaxspcell.F90         ; $(F90) ../changemaxspcell.F90
changemaxspcellbo.o       : ../changemaxspcellbo.F90       ; $(F90) ../changemaxspcellbo.F90
changemaxspec.o           : ../changemaxspec.F90           ; $(F90) ../changemaxspec.F90
changemaxsymop.o          : ../changemaxsymop.F90          ; $(F90) ../changemaxsymop.F90
changemaxtdfield.o        : ../changemaxtdfield.F90        ; $(F90) ../changemaxtdfield.F90
changemaxtempramp.o       : ../changemaxtempramp.F90       ; $(F90) ../changemaxtempramp.F90
changemaxthb.o            : ../changemaxthb.F90            ; $(F90) ../changemaxthb.F90
changemaxtitle.o          : ../changemaxtitle.F90          ; $(F90) ../changemaxtitle.F90
changemaxtotr1at.o        : ../changemaxtotr1at.F90        ; $(F90) ../changemaxtotr1at.F90
changemaxtotint.o         : ../changemaxtotint.F90         ; $(F90) ../changemaxtotint.F90
changemaxtotvac.o         : ../changemaxtotvac.F90         ; $(F90) ../changemaxtotvac.F90
changemaxtrialatom.o      : ../changemaxtrialatom.F90      ; $(F90) ../changemaxtrialatom.F90
changemaxuffspec.o        : ../changemaxuffspec.F90        ; $(F90) ../changemaxuffspec.F90
changemaxint.o            : ../changemaxint.F90            ; $(F90) ../changemaxint.F90
changemaxvac.o            : ../changemaxvac.F90            ; $(F90) ../changemaxvac.F90
changemaxvar.o            : ../changemaxvar.F90            ; $(F90) ../changemaxvar.F90
changemaxvectorpair.o     : ../changemaxvectorpair.F90     ; $(F90) ../changemaxvectorpair.F90
changemaxword.o           : ../changemaxword.F90           ; $(F90) ../changemaxword.F90
checkdefshspec.o          : ../checkdefshspec.F90          ; $(F90) ../checkdefshspec.F90
checkpot.o                : ../checkpot.F90                ; $(F90) ../checkpot.F90
checkquaternions.o        : ../checkquaternions.F90        ; $(F90) ../checkquaternions.F90
chvo.o                    : ../chvo.F90                    ; $(F90) ../chvo.F90
closestimage.o            : ../closestimage.F90            ; $(F90) ../closestimage.F90
collect.o                 : ../collect.F90                 ; $(F90) ../collect.F90
compressd1.o              : ../compressd1.F90              ; $(F90) ../compressd1.F90
compressd2.o              : ../compressd2.F90              ; $(F90) ../compressd2.F90
compressd3.o              : ../compressd3.F90              ; $(F90) ../compressd3.F90
connectwrap.o             : ../connectwrap.F90             ; $(F90) ../connectwrap.F90
connectwraptrial.o        : ../connectwraptrial.F90        ; $(F90) ../connectwraptrial.F90
connectunwrap.o           : ../connectunwrap.F90           ; $(F90) ../connectunwrap.F90
copycfg.o                 : ../copycfg.F90                 ; $(F90) ../copycfg.F90
cosmoaiiderv.o            : ../cosmoaiiderv.F90            ; $(F90) ../cosmoaiiderv.F90
cosmoaijderv.o            : ../cosmoaijderv.F90            ; $(F90) ../cosmoaijderv.F90
cosmoaiizero.o            : ../cosmoaiizero.F90            ; $(F90) ../cosmoaiizero.F90
cosmoaijzero.o            : ../cosmoaijzero.F90            ; $(F90) ../cosmoaijzero.F90
cosmoamatdadd.o           : ../cosmoamatdadd.F90           ; $(F90) ../cosmoamatdadd.F90
cosmoamatdadd1p.o         : ../cosmoamatdadd1p.F90         ; $(F90) ../cosmoamatdadd1p.F90
cosmoamatdaddd.o          : ../cosmoamatdaddd.F90          ; $(F90) ../cosmoamatdaddd.F90
cosmoderv.o               : ../cosmoderv.F90               ; $(F90) ../cosmoderv.F90
cosmoderv2D.o             : ../cosmoderv2D.F90             ; $(F90) ../cosmoderv2D.F90
costfn.o                  : ../costfn.F90                  ; $(F90) ../costfn.F90
create.o                  : ../create.F90                  ; $(F90) ../create.F90
createsphere.o            : ../createsphere.F90            ; $(F90) ../createsphere.F90
cscheck.o                 : ../cscheck.F90                 ; $(F90) ../cscheck.F90
cutscheck.o               : ../cutscheck.F90               ; $(F90) ../cutscheck.F90
datetime.o                : ../datetime.F90                ; $(F90) ../datetime.F90
dbcgsolve.o               : ../dbcgsolve.F90               ; $(F90) ../dbcgsolve.F90
dcharge.o                 : ../dcharge.F90                 ; $(F90) ../dcharge.F90
dcharges.o                : ../dcharges.F90                ; $(F90) ../dcharges.F90
d1charge.o                : ../d1charge.F90                ; $(F90) ../d1charge.F90
d1charge3.o               : ../d1charge3.F90               ; $(F90) ../d1charge3.F90
d1charged.o               : ../d1charged.F90               ; $(F90) ../d1charged.F90
d1charges.o               : ../d1charges.F90               ; $(F90) ../d1charges.F90
d2charge.o                : ../d2charge.F90                ; $(F90) ../d2charge.F90
d2charge3.o               : ../d2charge3.F90               ; $(F90) ../d2charge3.F90
d2charge3fc.o             : ../d2charge3fc.F90             ; $(F90) ../d2charge3fc.F90
d2charge3p.o              : ../d2charge3p.F90              ; $(F90) ../d2charge3p.F90
d2chargefc.o              : ../d2chargefc.F90              ; $(F90) ../d2chargefc.F90
d2chargep.o               : ../d2chargep.F90               ; $(F90) ../d2chargep.F90
d2charge2D.o              : ../d2charge2D.F90              ; $(F90) ../d2charge2D.F90
d2chargep2D.o             : ../d2chargep2D.F90             ; $(F90) ../d2chargep2D.F90
dcosmowt.o                : ../dcosmowt.F90                ; $(F90) ../dcosmowt.F90
dtotcosmowt.o             : ../dtotcosmowt.F90             ; $(F90) ../dtotcosmowt.F90
defect.o                  : ../defect.F90                  ; $(F90) ../defect.F90
defener.o                 : ../defener.F90                 ; $(F90) ../defener.F90
defequ.o                  : ../defequ.F90                  ; $(F90) ../defequ.F90
defequ2a.o                : ../defequ2a.F90                ; $(F90) ../defequ2a.F90
deffreq.o                 : ../deffreq.F90                 ; $(F90) ../deffreq.F90
deffun.o                  : ../deffun.F90                  ; $(F90) ../deffun.F90
defoin.o                  : ../defoin.F90                  ; $(F90) ../defoin.F90
defopt.o                  : ../defopt.F90                  ; $(F90) ../defopt.F90
defout.o                  : ../defout.F90                  ; $(F90) ../defout.F90
defpot.o                  : ../defpot.F90                  ; $(F90) ../defpot.F90
defsym.o                  : ../defsym.F90                  ; $(F90) ../defsym.F90
deftmat.o                 : ../deftmat.F90                 ; $(F90) ../deftmat.F90
defword.o                 : ../defword.F90                 ; $(F90) ../defword.F90
density0.o                : ../density0.F90                ; $(F90) ../density0.F90
density112.o              : ../density112.F90              ; $(F90) ../density112.F90
density11s2.o             : ../density11s2.F90             ; $(F90) ../density11s2.F90
density1fc.o              : ../density1fc.F90              ; $(F90) ../density1fc.F90
density12a2.o             : ../density12a2.F90             ; $(F90) ../density12a2.F90
density12as2.o            : ../density12as2.F90            ; $(F90) ../density12as2.F90
density2a.o               : ../density2a.F90               ; $(F90) ../density2a.F90
density3.o                : ../density3.F90                ; $(F90) ../density3.F90
density3s.o               : ../density3s.F90               ; $(F90) ../density3s.F90
densitysd.o               : ../densitysd.F90               ; $(F90) ../densitysd.F90
derfc.o                   : ../derfc.F90                   ; $(F90) ../derfc.F90
destroy.o                 : ../destroy.F90                 ; $(F90) ../destroy.F90
diagdyn.o                 : ../diagdyn.F90                 ; $(F90) ../diagdyn.F90
diagonal_fc.o             : ../diagonal_fc.F90             ; $(F90) ../diagonal_fc.F90
diagonal_fcd.o            : ../diagonal_fcd.F90            ; $(F90) ../diagonal_fcd.F90
dipole3D.o                : ../dipole3D.F90                ; $(F90) ../dipole3D.F90
dipole3Dd.o               : ../dipole3Dd.F90               ; $(F90) ../dipole3Dd.F90
dipole3Dtrial.o           : ../dipole3Dtrial.F90           ; $(F90) ../dipole3Dtrial.F90
dipolesd.o                : ../dipolesd.F90                ; $(F90) ../dipolesd.F90
distance.o                : ../distance.F90                ; $(F90) ../distance.F90
dlamch.o                  : ../dlamch.F90                  ; $(F90O0) ../dlamch.F90
dpoccon.o                 : ../dpoccon.F90                 ; $(F90) ../dpoccon.F90
dspecial.o                : ../dspecial.F90                ; $(F90) ../dspecial.F90
dump34.o                  : ../dump34.F90                  ; $(F90O1) ../dump34.F90
dumpbo.o                  : ../dumpbo.F90                  ; $(F90O1) ../dumpbo.F90
dumpdur.o                 : ../dumpdur.F90                 ; $(F90O1) ../dumpdur.F90
dumpgen.o                 : ../dumpgen.F90                 ; $(F90O1) ../dumpgen.F90
dumppot.o                 : ../dumppot.F90                 ; $(F90O1) ../dumppot.F90
dutils.o                  : ../dutils.F90                  ; $(F90) ../dutils.F90
dynamic.o                 : ../dynamic.F90                 ; $(F90) ../dynamic.F90
dynamic_fc.o              : ../dynamic_fc.F90              ; $(F90) ../dynamic_fc.F90
dynamicn.o                : ../dynamicn.F90                ; $(F90) ../dynamicn.F90
dynamicn_fc.o             : ../dynamicn_fc.F90             ; $(F90) ../dynamicn_fc.F90
dynamicn_fcd.o            : ../dynamicn_fcd.F90            ; $(F90) ../dynamicn_fcd.F90
dynamicn_defect.o         : ../dynamicn_defect.F90         ; $(F90) ../dynamicn_defect.F90
eamfnderv.o               : ../eamfnderv.F90               ; $(F90O1) ../eamfnderv.F90
eamrho.o                  : ../eamrho.F90                  ; $(F90) ../eamrho.F90
eamscalescrho.o           : ../eamscalescrho.F90           ; $(F90) ../eamscalescrho.F90
edip1fc.o                 : ../edip1fc.F90                 ; $(F90) ../edip1fc.F90
edip2body.o               : ../edip2body.F90               ; $(F90) ../edip2body.F90
edip3body.o               : ../edip3body.F90               ; $(F90) ../edip3body.F90
edipcn.o                  : ../edipcn.F90                  ; $(F90) ../edipcn.F90
edipmd.o                  : ../edipmd.F90                  ; $(F90) ../edipmd.F90
eem.o                     : ../eem.F90                     ; $(F90) ../eem.F90
eemd.o                    : ../eemd.F90                    ; $(F90) ../eemd.F90
eemdm.o                   : ../eemdm.F90                   ; $(F90) ../eemdm.F90
eemsplit.o                : ../eemsplit.F90                ; $(F90) ../eemsplit.F90
eemsplitd.o               : ../eemsplitd.F90               ; $(F90) ../eemsplitd.F90
electricfield.o           : ../electricfield.F90           ; $(F90) ../electricfield.F90
electricfieldparts.o      : ../electricfieldparts.F90      ; $(F90) ../electricfieldparts.F90
electricfieldpotl.o       : ../electricfieldpotl.F90       ; $(F90) ../electricfieldpotl.F90
electricfieldtrial.o      : ../electricfieldtrial.F90      ; $(F90) ../electricfieldtrial.F90
einstein.o                : ../einstein.F90                ; $(F90) ../einstein.F90
einsteinpd.o              : ../einsteinpd.F90              ; $(F90) ../einsteinpd.F90
einsteinsd2.o             : ../einsteinsd2.F90             ; $(F90) ../einsteinsd2.F90
einsteintrial.o           : ../einsteintrial.F90           ; $(F90) ../einsteintrial.F90
eispack.o                 : ../eispack.F90                 ; $(F90) ../eispack.F90
emfunc.o                  : ../emfunc.F90                  ; $(F90) ../emfunc.F90
endstring.o               : ../endstring.F90               ; $(F90) ../endstring.F90
epot.o                    : ../epot.F90                    ; $(F90) ../epot.F90
epsilon.o                 : ../epsilon.F90                 ; $(F90) ../epsilon.F90
equpos.o                  : ../equpos.F90                  ; $(F90) ../equpos.F90
exchng.o                  : ../exchng.F90                  ; $(F90) ../exchng.F90
fefunct.o                 : ../fefunct.F90                 ; $(F90) ../fefunct.F90
fefunctn.o                : ../fefunctn.F90                ; $(F90) ../fefunctn.F90
fenergy0.o                : ../fenergy0.F90                ; $(F90) ../fenergy0.F90
fenergy3.o                : ../fenergy3.F90                ; $(F90) ../fenergy3.F90
findcelltype.o            : ../findcelltype.F90            ; $(F90) ../findcelltype.F90
finddef.o                 : ../finddef.F90                 ; $(F90) ../finddef.F90
findfractionword.o        : ../findfractionword.F90        ; $(F90) ../findfractionword.F90
findspacegroup.o          : ../findspacegroup.F90          ; $(F90) ../findspacegroup.F90
firstpass.o               : ../firstpass.F90               ; $(F90) ../firstpass.F90
fit.o                     : ../fit.F90                     ; $(F90) ../fit.F90
fitcon.o                  : ../fitcon.F90                  ; $(F90) ../fitcon.F90
fitfun.o                  : ../fitfun.F90                  ; $(F90) ../fitfun.F90
fitfunct.o                : ../fitfunct.F90                ; $(F90) ../fitfunct.F90
fitgrad.o                 : ../fitgrad.F90                 ; $(F90) ../fitgrad.F90
fitmin.o                  : ../fitmin.F90                  ; $(F90) ../fitmin.F90
fitword.o                 : ../fitword.F90                 ; $(F90) ../fitword.F90
force.o                   : ../force.F90                   ; $(F90) ../force.F90
forcebycell.o             : ../forcebycell.F90             ; $(F90) ../forcebycell.F90
forcetrial.o              : ../forcetrial.F90              ; $(F90) ../forcetrial.F90
formula.o                 : ../formula.F90                 ; $(F90) ../formula.F90
fournos.o                 : ../fournos.F90                 ; $(F90) ../fournos.F90
fournosd.o                : ../fournosd.F90                ; $(F90) ../fournosd.F90
four12.o                  : ../four12.F90                  ; $(F90) ../four12.F90
four12d.o                 : ../four12d.F90                 ; $(F90) ../four12d.F90
fourbody.o                : ../fourbody.F90                ; $(F90) ../fourbody.F90
fourlist.o                : ../fourlist.F90                ; $(F90) ../fourlist.F90
fourimp.o                 : ../fourimp.F90                 ; $(F90) ../fourimp.F90
fourimpd.o                : ../fourimpd.F90                ; $(F90) ../fourimpd.F90
fourimp12.o               : ../fourimp12.F90               ; $(F90) ../fourimp12.F90
fourimp12d.o              : ../fourimp12d.F90              ; $(F90) ../fourimp12d.F90
fourimp1fc.o              : ../fourimp1fc.F90              ; $(F90) ../fourimp1fc.F90
fourimpfc.o               : ../fourimpfc.F90               ; $(F90) ../fourimpfc.F90
fourimpfc3.o              : ../fourimpfc3.F90              ; $(F90) ../fourimpfc3.F90
fourimpfcd.o              : ../fourimpfcd.F90              ; $(F90) ../fourimpfcd.F90
fourimpfd.o               : ../fourimpfd.F90               ; $(F90) ../fourimpfd.F90
fourimpmc.o               : ../fourimpmc.F90               ; $(F90) ../fourimpmc.F90
fourimpmd.o               : ../fourimpmd.F90               ; $(F90) ../fourimpmd.F90
fourimpmds.o              : ../fourimpmds.F90              ; $(F90) ../fourimpmds.F90
fourimpp.o                : ../fourimpp.F90                ; $(F90) ../fourimpp.F90
fourimppd.o               : ../fourimppd.F90               ; $(F90) ../fourimppd.F90
fourimpsd.o               : ../fourimpsd.F90               ; $(F90) ../fourimpsd.F90
fourimpsd2.o              : ../fourimpsd2.F90              ; $(F90) ../fourimpsd2.F90
fouroop.o                 : ../fouroop.F90                 ; $(F90) ../fouroop.F90
fouroopd.o                : ../fouroopd.F90                ; $(F90) ../fouroopd.F90
fouroop1fc.o              : ../fouroop1fc.F90              ; $(F90) ../fouroop1fc.F90
fouroopfc.o               : ../fouroopfc.F90               ; $(F90) ../fouroopfc.F90
fouroopfc3.o              : ../fouroopfc3.F90              ; $(F90) ../fouroopfc3.F90
fouroopfcd.o              : ../fouroopfcd.F90              ; $(F90) ../fouroopfcd.F90
fouroopfd.o               : ../fouroopfd.F90               ; $(F90) ../fouroopfd.F90
fouroopp.o                : ../fouroopp.F90                ; $(F90) ../fouroopp.F90
fourooppd.o               : ../fourooppd.F90               ; $(F90) ../fourooppd.F90
fouroopmc.o               : ../fouroopmc.F90               ; $(F90) ../fouroopmc.F90
fouroopmd.o               : ../fouroopmd.F90               ; $(F90) ../fouroopmd.F90
fouroopmds.o              : ../fouroopmds.F90              ; $(F90) ../fouroopmds.F90
fouroopsd.o               : ../fouroopsd.F90               ; $(F90) ../fouroopsd.F90
fouroop12.o               : ../fouroop12.F90               ; $(F90) ../fouroop12.F90
fouroop12d.o              : ../fouroop12d.F90              ; $(F90) ../fouroop12d.F90
fouroopsd2.o              : ../fouroopsd2.F90              ; $(F90) ../fouroopsd2.F90
four1fc.o                 : ../four1fc.F90                 ; $(F90) ../four1fc.F90
fourfc.o                  : ../fourfc.F90                  ; $(F90) ../fourfc.F90
fourfd.o                  : ../fourfd.F90                  ; $(F90) ../fourfd.F90
fourfcd.o                 : ../fourfcd.F90                 ; $(F90) ../fourfcd.F90
fourp.o                   : ../fourp.F90                   ; $(F90O1) ../fourp.F90
fourpd.o                  : ../fourpd.F90                  ; $(F90O1) ../fourpd.F90
fourmc.o                  : ../fourmc.F90                  ; $(F90) ../fourmc.F90
fourmd.o                  : ../fourmd.F90                  ; $(F90) ../fourmd.F90
fourmds.o                 : ../fourmds.F90                 ; $(F90) ../fourmds.F90
foursd.o                  : ../foursd.F90                  ; $(F90) ../foursd.F90
foursd2.o                 : ../foursd2.F90                 ; $(F90) ../foursd2.F90
four0d3.o                 : ../four0d3.F90                 ; $(F90) ../four0d3.F90
four3d3.o                 : ../four3d3.F90                 ; $(F90) ../four3d3.F90
four3d3dV.o               : ../four3d3dV.F90               ; $(F90) ../four3d3dV.F90
four3d3dVg.o              : ../four3d3dVg.F90              ; $(F90) ../four3d3dVg.F90
four3d3g.o                : ../four3d3g.F90                ; $(F90) ../four3d3g.F90
fourfc3.o                 : ../fourfc3.F90                 ; $(F90) ../fourfc3.F90
ftow.o                    : ../ftow.F90                    ; $(F90) ../ftow.F90
funct.o                   : ../funct.F90                   ; $(F90) ../funct.F90
functn.o                  : ../functn.F90                  ; $(F90) ../functn.F90
functnf.o                 : ../functnf.F90                 ; $(F90) ../functnf.F90
functnff.o                : ../functnff.F90                ; $(F90) ../functnff.F90
functneb.o                : ../functneb.F90                ; $(F90) ../functneb.F90
functsync.o               : ../functsync.F90               ; $(F90) ../functsync.F90
gabcfg.o                  : ../gabcfg.F90                  ; $(F90) ../gabcfg.F90
gacreate.o                : ../gacreate.F90                ; $(F90) ../gacreate.F90
gacross.o                 : ../gacross.F90                 ; $(F90) ../gacross.F90
gaexpd.o                  : ../gaexpd.F90                  ; $(F90) ../gaexpd.F90
gafit.o                   : ../gafit.F90                   ; $(F90) ../gafit.F90
gamaxmin.o                : ../gamaxmin.F90                ; $(F90) ../gamaxmin.F90
gamma.o                   : ../gamma.F90                   ; $(F90) ../gamma.F90
gammasm.o                 : ../gammasm.F90                 ; $(F90) ../gammasm.F90
gamuta.o                  : ../gamuta.F90                  ; $(F90) ../gamuta.F90
gaopt.o                   : ../gaopt.F90                   ; $(F90) ../gaopt.F90
garepro.o                 : ../garepro.F90                 ; $(F90) ../garepro.F90
gasort.o                  : ../gasort.F90                  ; $(F90) ../gasort.F90
gasteiger.o               : ../gasteiger.F90               ; $(F90) ../gasteiger.F90
gasubs.o                  : ../gasubs.F90                  ; $(F90) ../gasubs.F90
gauss.o                   : ../gauss.F90                   ; $(F90) ../gauss.F90
gaword.o                  : ../gaword.F90                  ; $(F90) ../gaword.F90
genpot.o                  : ../genpot.F90                  ; $(F90) ../genpot.F90
genpotd.o                 : ../genpotd.F90                 ; $(F90) ../genpotd.F90
genpotdm.o                : ../genpotdm.F90                ; $(F90) ../genpotdm.F90
genpotsplitd.o            : ../genpotsplitd.F90            ; $(F90) ../genpotsplitd.F90
genword.o                 : ../genword.F90                 ; $(F90) ../genword.F90
getangles.o               : ../getangles.F90               ; $(F90) ../getangles.F90
getbocharge.o             : ../getbocharge.F90             ; $(F90) ../getbocharge.F90
getbocharged.o            : ../getbocharged.F90            ; $(F90) ../getbocharged.F90
getbochargetrial.o        : ../getbochargetrial.F90        ; $(F90) ../getbochargetrial.F90
getbonds.o                : ../getbonds.F90                ; $(F90) ../getbonds.F90
getcoordno.o              : ../getcoordno.F90              ; $(F90) ../getcoordno.F90
getderv1.o                : ../getderv1.F90                ; $(F90) ../getderv1.F90
getdipole.o               : ../getdipole.F90               ; $(F90) ../getdipole.F90
getdistances.o            : ../getdistances.F90            ; $(F90) ../getdistances.F90
getgcindex.o              : ../getgcindex.F90              ; $(F90) ../getgcindex.F90
getkeywrd.o               : ../getkeywrd.F90               ; $(F90) ../getkeywrd.F90
getneighbour.o            : ../getneighbour.F90            ; $(F90) ../getneighbour.F90
getnmr.o                  : ../getnmr.F90                  ; $(F90) ../getnmr.F90
getphi.o                  : ../getphi.F90                  ; $(F90) ../getphi.F90
getpotsymbol.o            : ../getpotsymbol.F90            ; $(F90) ../getpotsymbol.F90
getshellmatrices.o        : ../getshellmatrices.F90        ; $(F90) ../getshellmatrices.F90
getstrain.o               : ../getstrain.F90               ; $(F90) ../getstrain.F90
gettorsions.o             : ../gettorsions.F90             ; $(F90) ../gettorsions.F90
gettraceinertia.o         : ../gettraceinertia.F90         ; $(F90) ../gettraceinertia.F90
getvibmode.o              : ../getvibmode.F90              ; $(F90) ../getvibmode.F90
getvibmoded.o             : ../getvibmoded.F90             ; $(F90) ../getvibmoded.F90
groupvelocities.o         : ../groupvelocities.F90         ; $(F90) ../groupvelocities.F90
grueneisen.o              : ../grueneisen.F90              ; $(F90) ../grueneisen.F90
grueneisend.o             : ../grueneisend.F90             ; $(F90) ../grueneisend.F90
grueneiseng.o             : ../grueneiseng.F90             ; $(F90) ../grueneiseng.F90
grueneisengd.o            : ../grueneisengd.F90            ; $(F90) ../grueneisengd.F90
gtheta.o                  : ../gtheta.F90                  ; $(F90) ../gtheta.F90
gthetabo.o                : ../gthetabo.F90                ; $(F90) ../gthetabo.F90
gulpfinish.o              : ../gulpfinish.F90              ; $(F90) ../gulpfinish.F90
help.o                    : ../help.F90                    ; $(F90) ../help.F90
hfunc.o                   : ../hfunc.F90                   ; $(F90) ../hfunc.F90
improper.o                : ../improper.F90                ; $(F90) ../improper.F90
indsft.o                  : ../indsft.F90                  ; $(F90) ../indsft.F90
ind1toijk.o               : ../ind1toijk.F90               ; $(F90) ../ind1toijk.F90
indtoijk.o                : ../indtoijk.F90                ; $(F90) ../indtoijk.F90
initarc.o                 : ../initarc.F90                 ; $(F90) ../initarc.F90
init1bodydefaults.o       : ../init1bodydefaults.F90       ; $(F90) ../init1bodydefaults.F90
init2bodydefaults.o       : ../init2bodydefaults.F90       ; $(F90) ../init2bodydefaults.F90
initmaxatdefaults.o       : ../initmaxatdefaults.F90       ; $(F90) ../initmaxatdefaults.F90
initmaxatotdefaults.o     : ../initmaxatotdefaults.F90     ; $(F90) ../initmaxatotdefaults.F90
initmaxbonddefaults.o     : ../initmaxbonddefaults.F90     ; $(F90) ../initmaxbonddefaults.F90
initmaxbondqdefaults.o    : ../initmaxbondqdefaults.F90    ; $(F90) ../initmaxbondqdefaults.F90
initmaxccspecdefaults.o   : ../initmaxccspecdefaults.F90   ; $(F90) ../initmaxccspecdefaults.F90
initmaxconnectdefaults.o  : ../initmaxconnectdefaults.F90  ; $(F90) ../initmaxconnectdefaults.F90
initmaxdefdefaults.o      : ../initmaxdefdefaults.F90      ; $(F90) ../initmaxdefdefaults.F90
initmaxeamdendefaults.o   : ../initmaxeamdendefaults.F90   ; $(F90) ../initmaxeamdendefaults.F90
initmaxeamfnspecdefaults.o: ../initmaxeamfnspecdefaults.F90; $(F90) ../initmaxeamfnspecdefaults.F90
initmaxeamspecdefaults.o  : ../initmaxeamspecdefaults.F90  ; $(F90) ../initmaxeamspecdefaults.F90
initmaxedipspecdefaults.o : ../initmaxedipspecdefaults.F90 ; $(F90) ../initmaxedipspecdefaults.F90
initmaxfgraddefaults.o    : ../initmaxfgraddefaults.F90    ; $(F90) ../initmaxfgraddefaults.F90
initmaxfitdefaults.o      : ../initmaxfitdefaults.F90      ; $(F90) ../initmaxfitdefaults.F90
initmaxfordefaults.o      : ../initmaxfordefaults.F90      ; $(F90) ../initmaxfordefaults.F90
initmaxfstraindefaults.o  : ../initmaxfstraindefaults.F90  ; $(F90) ../initmaxfstraindefaults.F90
initmaxgcmcmoldefaults.o  : ../initmaxgcmcmoldefaults.F90  ; $(F90) ../initmaxgcmcmoldefaults.F90
initmaxlibdefaults.o      : ../initmaxlibdefaults.F90      ; $(F90) ../initmaxlibdefaults.F90
initmaxqrangedefaults.o   : ../initmaxqrangedefaults.F90   ; $(F90) ../initmaxqrangedefaults.F90
initmaxmcswaps.o          : ../initmaxmcswaps.F90          ; $(F90) ../initmaxmcswaps.F90
initmaxmcswapspec.o       : ../initmaxmcswapspec.F90       ; $(F90) ../initmaxmcswapspec.F90
initmaxmoldefaults.o      : ../initmaxmoldefaults.F90      ; $(F90) ../initmaxmoldefaults.F90
initmaxnboadefaults.o     : ../initmaxnboadefaults.F90     ; $(F90) ../initmaxnboadefaults.F90
initmaxnboodefaults.o     : ../initmaxnboodefaults.F90     ; $(F90) ../initmaxnboodefaults.F90
initmaxnbopotdefaults.o   : ../initmaxnbopotdefaults.F90   ; $(F90) ../initmaxnbopotdefaults.F90
initmaxnboq0defaults.o    : ../initmaxnboq0defaults.F90    ; $(F90) ../initmaxnboq0defaults.F90
initmaxnboqdefaults.o     : ../initmaxnboqdefaults.F90     ; $(F90) ../initmaxnboqdefaults.F90
initmaxnbordefaults.o     : ../initmaxnbordefaults.F90     ; $(F90) ../initmaxnbordefaults.F90
initmaxnbozdefaults.o     : ../initmaxnbozdefaults.F90     ; $(F90) ../initmaxnbozdefaults.F90
initmaxnebreplicatotdefaults.o : ../initmaxnebreplicatotdefaults.F90 ; $(F90) ../initmaxnebreplicatotdefaults.F90
initmaxnppadefaults.o     : ../initmaxnppadefaults.F90     ; $(F90) ../initmaxnppadefaults.F90
initmaxnptsdefaults.o     : ../initmaxnptsdefaults.F90     ; $(F90) ../initmaxnptsdefaults.F90
initmaxobsdefaults.o      : ../initmaxobsdefaults.F90      ; $(F90) ../initmaxobsdefaults.F90
initmaxplanepotdefaults.o : ../initmaxplanepotdefaults.F90 ; $(F90) ../initmaxplanepotdefaults.F90
initmaxr1atdefaults.o     : ../initmaxr1atdefaults.F90     ; $(F90) ../initmaxr1atdefaults.F90
initmaxreaxffspecdefaults.o : ../initmaxreaxffspecdefaults.F90 ; $(F90) ../initmaxreaxffspecdefaults.F90
initmaxreaxffval3defaults.o : ../initmaxreaxffval3defaults.F90 ; $(F90) ../initmaxreaxffval3defaults.F90
initmaxregiondefaults.o   : ../initmaxregiondefaults.F90   ; $(F90) ../initmaxregiondefaults.F90
initmaxsixdefaults.o      : ../initmaxsixdefaults.F90      ; $(F90) ../initmaxsixdefaults.F90
initmaxspcellbodefaults.o : ../initmaxspcellbodefaults.F90 ; $(F90) ../initmaxspcellbodefaults.F90
initmaxspcelldefaults.o   : ../initmaxspcelldefaults.F90   ; $(F90) ../initmaxspcelldefaults.F90
initmaxspecdefaults.o     : ../initmaxspecdefaults.F90     ; $(F90) ../initmaxspecdefaults.F90
initmaxthbdefaults.o      : ../initmaxthbdefaults.F90      ; $(F90) ../initmaxthbdefaults.F90
initmaxtitledefaults.o    : ../initmaxtitledefaults.F90    ; $(F90) ../initmaxtitledefaults.F90
init_ins.o                : ../init_ins.F90                ; $(F90) ../init_ins.F90
initktrm.o                : ../initktrm.F90                ; $(F90) ../initktrm.F90
initmemory.o              : ../initmemory.F90              ; $(F90) ../initmemory.F90
initqmatrix.o             : ../initqmatrix.F90             ; $(F90) ../initqmatrix.F90
initqmatrixc.o            : ../initqmatrixc.F90            ; $(F90) ../initqmatrixc.F90
insertatom.o              : ../insertatom.F90              ; $(F90) ../insertatom.F90
inword.o                  : ../inword.F90                  ; $(F90O1) ../inword.F90
istoan.o                  : ../istoan.F90                  ; $(F90) ../istoan.F90
itcom.o                   : ../itcom.F90                   ; $(F90) ../itcom.F90
iterinv.o                 : ../iterinv.F90                 ; $(F90) ../iterinv.F90
itow.o                    : ../itow.F90                    ; $(F90) ../itow.F90
itpack.o                  : ../itpack.F90                  ; $(F90) ../itpack.F90
kindex.o                  : ../kindex.F90                  ; $(F90) ../kindex.F90
kinthewedge.o             : ../kinthewedge.F90             ; $(F90) ../kinthewedge.F90
kvector.o                 : ../kvector.F90                 ; $(F90) ../kvector.F90
kvectorf.o                : ../kvectorf.F90                ; $(F90) ../kvectorf.F90
kweight.o                 : ../kweight.F90                 ; $(F90) ../kweight.F90
label.o                   : ../label.F90                   ; $(F90) ../label.F90
labellc.o                 : ../labellc.F90                 ; $(F90) ../labellc.F90
lambda.o                  : ../lambda.F90                  ; $(F90O1) ../lambda.F90
langevinpol.o             : ../langevinpol.F90             ; $(F90) ../langevinpol.F90
lapack.o                  : ../lapack.F90                  ; $(F90O2) ../lapack.F90
latpot.o                  : ../latpot.F90                  ; $(F90) ../latpot.F90
lmbfgs.o                  : ../lmbfgs.F90                  ; $(F90) ../lmbfgs.F90
linepro.o                 : ../linepro.F90                 ; $(F90) ../linepro.F90
linpack.o                 : ../linpack.F90                 ; $(F90) ../linpack.F90
lintoijk.o                : ../lintoijk.F90                ; $(F90) ../lintoijk.F90
lmatch.o                  : ../lmatch.F90                  ; $(F90) ../lmatch.F90
lmatchany.o               : ../lmatchany.F90               ; $(F90) ../lmatchany.F90
lmatchpair.o              : ../lmatchpair.F90              ; $(F90) ../lmatchpair.F90
lower.o                   : ../lower.F90                   ; $(F90) ../lower.F90
ltont.o                   : ../ltont.F90                   ; $(F90) ../ltont.F90
m_pdf.o                   : ../m_pdf.F90                   ; $(F90) ../m_pdf.F90
m_pdfneutron.o            : ../m_pdfneutron.F90            ; $(F90) ../m_pdfneutron.F90
m_pdfvariables.o          : ../m_pdfvariables.F90          ; $(F90) ../m_pdfvariables.F90
many.o                    : ../many.F90                    ; $(F90) ../many.F90
manyd.o                   : ../manyd.F90                   ; $(F90) ../manyd.F90
many1fc.o                 : ../many1fc.F90                 ; $(F90) ../many1fc.F90
manyfc.o                  : ../manyfc.F90                  ; $(F90) ../manyfc.F90
manyfcd.o                 : ../manyfcd.F90                 ; $(F90) ../manyfcd.F90
manyp.o                   : ../manyp.F90                   ; $(F90) ../manyp.F90
manypd.o                  : ../manypd.F90                  ; $(F90) ../manypd.F90
manysd.o                  : ../manysd.F90                  ; $(F90) ../manysd.F90
many0d.o                  : ../many0d.F90                  ; $(F90) ../many0d.F90
many0dd.o                 : ../many0dd.F90                 ; $(F90) ../many0dd.F90
many12.o                  : ../many12.F90                  ; $(F90) ../many12.F90
many12d.o                 : ../many12d.F90                 ; $(F90) ../many12d.F90
manymd0.o                 : ../manymd0.F90                 ; $(F90) ../manymd0.F90
manymd3.o                 : ../manymd3.F90                 ; $(F90) ../manymd3.F90
manymd3s.o                : ../manymd3s.F90                ; $(F90) ../manymd3s.F90
manymi3.o                 : ../manymi3.F90                 ; $(F90) ../manymi3.F90
manysd2.o                 : ../manysd2.F90                 ; $(F90) ../manysd2.F90
mccreate.o                : ../mccreate.F90                ; $(F90) ../mccreate.F90
mctrialfunction.o         : ../mctrialfunction.F90         ; $(F90) ../mctrialfunction.F90
mcdestroy.o               : ../mcdestroy.F90               ; $(F90) ../mcdestroy.F90
mcmolinsert.o             : ../mcmolinsert.F90             ; $(F90) ../mcmolinsert.F90
mcmove.o                  : ../mcmove.F90                  ; $(F90) ../mcmove.F90
mcrotate.o                : ../mcrotate.F90                ; $(F90) ../mcrotate.F90
mcstrain.o                : ../mcstrain.F90                ; $(F90) ../mcstrain.F90
mcstorecfg.o              : ../mcstorecfg.F90              ; $(F90) ../mcstorecfg.F90
mcswap.o                  : ../mcswap.F90                  ; $(F90) ../mcswap.F90
mcword.o                  : ../mcword.F90                  ; $(F90) ../mcword.F90
mdcellind.o               : ../mdcellind.F90               ; $(F90) ../mdcellind.F90
mdcorrect.o               : ../mdcorrect.F90               ; $(F90) ../mdcorrect.F90
mddistconpost.o           : ../mddistconpost.F90           ; $(F90) ../mddistconpost.F90
mdgcorrect.o              : ../mdgcorrect.F90              ; $(F90) ../mdgcorrect.F90
mdvvcorrect.o             : ../mdvvcorrect.F90             ; $(F90) ../mdvvcorrect.F90
mdlfvcorrect.o            : ../mdlfvcorrect.F90            ; $(F90) ../mdlfvcorrect.F90
mdke.o                    : ../mdke.F90                    ; $(F90) ../mdke.F90
mdkeatomicstress.o        : ../mdkeatomicstress.F90        ; $(F90) ../mdkeatomicstress.F90
mdkestrain.o              : ../mdkestrain.F90              ; $(F90) ../mdkestrain.F90
mdpredict.o               : ../mdpredict.F90               ; $(F90) ../mdpredict.F90
mdprop.o                  : ../mdprop.F90                  ; $(F90) ../mdprop.F90
mdreset.o                 : ../mdreset.F90                 ; $(F90) ../mdreset.F90
mdrestart.o               : ../mdrestart.F90               ; $(F90) ../mdrestart.F90
mdscale.o                 : ../mdscale.F90                 ; $(F90) ../mdscale.F90
mdshscale.o               : ../mdshscale.F90               ; $(F90) ../mdshscale.F90
mdstop.o                  : ../mdstop.F90                  ; $(F90) ../mdstop.F90
mdtemperature.o           : ../mdtemperature.F90           ; $(F90) ../mdtemperature.F90
mdvelcor.o                : ../mdvelcor.F90                ; $(F90) ../mdvelcor.F90
mdword.o                  : ../mdword.F90                  ; $(F90) ../mdword.F90
mdwrite.o                 : ../mdwrite.F90                 ; $(F90) ../mdwrite.F90
meamfnderv.o              : ../meamfnderv.F90              ; $(F90) ../meamfnderv.F90
meamrho.o                 : ../meamrho.F90                 ; $(F90) ../meamrho.F90
meamscreen.o              : ../meamscreen.F90              ; $(F90) ../meamscreen.F90
meamsymdensity.o          : ../meamsymdensity.F90          ; $(F90) ../meamsymdensity.F90
meamtotalrho.o            : ../meamtotalrho.F90            ; $(F90) ../meamtotalrho.F90
meamtotalrhoderv.o        : ../meamtotalrhoderv.F90        ; $(F90) ../meamtotalrhoderv.F90
meamtotalrhocrossderv.o   : ../meamtotalrhocrossderv.F90   ; $(F90) ../meamtotalrhocrossderv.F90
meamtotalrhoscreenderv.o  : ../meamtotalrhoscreenderv.F90  ; $(F90) ../meamtotalrhoscreenderv.F90
methodok.o                : ../methodok.F90                ; $(F90) ../methodok.F90
mindtoijk.o               : ../mindtoijk.F90               ; $(F90) ../mindtoijk.F90
molind.o                  : ../molind.F90                  ; $(F90) ../molind.F90
molindfix.o               : ../molindfix.F90               ; $(F90) ../molindfix.F90
molindtrial.o             : ../molindtrial.F90             ; $(F90) ../molindtrial.F90
molindfixtrial.o          : ../molindfixtrial.F90          ; $(F90) ../molindfixtrial.F90
moltype.o                 : ../moltype.F90                 ; $(F90) ../moltype.F90
move2a1.o                 : ../move2a1.F90                 ; $(F90) ../move2a1.F90
mxmb.o                    : ../mxmb.F90                    ; $(F90) ../mxmb.F90
nagamma.o                 : ../nagamma.F90                 ; $(F90) ../nagamma.F90
nagammac.o                : ../nagammac.F90                ; $(F90) ../nagammac.F90
nagammacrigid.o           : ../nagammacrigid.F90           ; $(F90) ../nagammacrigid.F90
nagammarigid.o            : ../nagammarigid.F90            ; $(F90) ../nagammarigid.F90
nearest.o                 : ../nearest.F90                 ; $(F90) ../nearest.F90
nebforce.o                : ../nebforce.F90                ; $(F90) ../nebforce.F90
nebword.o                 : ../nebword.F90                 ; $(F90) ../nebword.F90
nudge.o                   : ../nudge.F90                   ; $(F90) ../nudge.F90
numhess.o                 : ../numhess.F90                 ; $(F90) ../numhess.F90
okspec.o                  : ../okspec.F90                  ; $(F90) ../okspec.F90
olinmin.o                 : ../olinmin.F90                 ; $(F90) ../olinmin.F90
omegaproperty.o           : ../omegaproperty.F90           ; $(F90) ../omegaproperty.F90
onebody.o                 : ../onebody.F90                 ; $(F90) ../onebody.F90
optin.o                   : ../optin.F90                   ; $(F90) ../optin.F90
options.o                 : ../options.F90                 ; $(F90) ../options.F90
optout.o                  : ../optout.F90                  ; $(F90) ../optout.F90
optshe.o                  : ../optshe.F90                  ; $(F90) ../optshe.F90
oscillatorstrength.o      : ../oscillatorstrength.F90      ; $(F90) ../oscillatorstrength.F90
oscillatorstrengthc.o     : ../oscillatorstrengthc.F90     ; $(F90) ../oscillatorstrengthc.F90
oscillatorstrengthg.o     : ../oscillatorstrengthg.F90     ; $(F90) ../oscillatorstrengthg.F90
outadvice.o               : ../outadvice.F90               ; $(F90) ../outadvice.F90
outarc.o                  : ../outarc.F90                  ; $(F90) ../outarc.F90
outatomicstress.o         : ../outatomicstress.F90         ; $(F90) ../outatomicstress.F90
outbiograf.o              : ../outbiograf.F90              ; $(F90) ../outbiograf.F90
outcell.o                 : ../outcell.F90                 ; $(F90) ../outcell.F90
outcif.o                  : ../outcif.F90                  ; $(F90) ../outcif.F90
outcssr.o                 : ../outcssr.F90                 ; $(F90) ../outcssr.F90
outdcdheader.o            : ../outdcdheader.F90            ; $(F90) ../outdcdheader.F90
outdefin.o                : ../outdefin.F90                ; $(F90) ../outdefin.F90
outdener.o                : ../outdener.F90                ; $(F90) ../outdener.F90
outderv.o                 : ../outderv.F90                 ; $(F90) ../outderv.F90
outdlv.o                  : ../outdlv.F90                  ; $(F90) ../outdlv.F90
outdrv.o                  : ../outdrv.F90                  ; $(F90) ../outdrv.F90
outed.o                   : ../outed.F90                   ; $(F90) ../outed.F90
outeem.o                  : ../outeem.F90                  ; $(F90) ../outeem.F90
outener.o                 : ../outener.F90                 ; $(F90) ../outener.F90
outerror.o                : ../outerror.F90                ; $(F90) ../outerror.F90
outfdf.o                  : ../outfdf.F90                  ; $(F90) ../outfdf.F90
outfile.o                 : ../outfile.F90                 ; $(F90) ../outfile.F90
outfit.o                  : ../outfit.F90                  ; $(F90) ../outfit.F90
outfitcon.o               : ../outfitcon.F90               ; $(F90) ../outfitcon.F90
outfour.o                 : ../outfour.F90                 ; $(F90) ../outfour.F90
outfrc.o                  : ../outfrc.F90                  ; $(F90) ../outfrc.F90
outfreq.o                 : ../outfreq.F90                 ; $(F90) ../outfreq.F90
outgafit.o                : ../outgafit.F90                ; $(F90) ../outgafit.F90
outgaopt.o                : ../outgaopt.F90                ; $(F90) ../outgaopt.F90
outgapar.o                : ../outgapar.F90                ; $(F90) ../outgapar.F90
outgen.o                  : ../outgen.F90                  ; $(F90) ../outgen.F90
outinertia.o              : ../outinertia.F90              ; $(F90) ../outinertia.F90
outkey.o                  : ../outkey.F90                  ; $(F90) ../outkey.F90
outkpt.o                  : ../outkpt.F90                  ; $(F90) ../outkpt.F90
outkptfile.o              : ../outkptfile.F90              ; $(F90) ../outkptfile.F90
outlammps.o               : ../outlammps.F90               ; $(F90) ../outlammps.F90
outlammpspots.o           : ../outlammpspots.F90           ; $(F90) ../outlammpspots.F90
outmarv.o                 : ../outmarv.F90                 ; $(F90) ../outmarv.F90
outmarv2.o                : ../outmarv2.F90                ; $(F90) ../outmarv2.F90
outmc.o                   : ../outmc.F90                   ; $(F90) ../outmc.F90
outmd.o                   : ../outmd.F90                   ; $(F90) ../outmd.F90
outmdin.o                 : ../outmdin.F90                 ; $(F90) ../outmdin.F90
outmol.o                  : ../outmol.F90                  ; $(F90) ../outmol.F90
outneb.o                  : ../outneb.F90                  ; $(F90) ../outneb.F90
outnebin.o                : ../outnebin.F90                ; $(F90) ../outnebin.F90
outobs.o                  : ../outobs.F90                  ; $(F90) ../outobs.F90
outofmemory.o             : ../outofmemory.F90             ; $(F90) ../outofmemory.F90
outofp.o                  : ../outofp.F90                  ; $(F90) ../outofp.F90
outoptvar.o               : ../outoptvar.F90               ; $(F90) ../outoptvar.F90
outphon.o                 : ../outphon.F90                 ; $(F90) ../outphon.F90
outpolar.o                : ../outpolar.F90                ; $(F90) ../outpolar.F90
outpot.o                  : ../outpot.F90                  ; $(F90) ../outpot.F90
outqbo.o                  : ../outqbo.F90                  ; $(F90) ../outqbo.F90
outshengBTE.o             : ../outshengBTE.F90             ; $(F90) ../outshengBTE.F90
outsiteenergy.o           : ../outsiteenergy.F90           ; $(F90) ../outsiteenergy.F90
outspec.o                 : ../outspec.F90                 ; $(F90) ../outspec.F90
outspec_pdf.o             : ../outspec_pdf.F90             ; $(F90) ../outspec_pdf.F90
outstresses.o             : ../outstresses.F90             ; $(F90) ../outstresses.F90
outstructure.o            : ../outstructure.F90            ; $(F90) ../outstructure.F90
outsync.o                 : ../outsync.F90                 ; $(F90) ../outsync.F90
outsyncin.o               : ../outsyncin.F90               ; $(F90) ../outsyncin.F90
outthb.o                  : ../outthb.F90                  ; $(F90) ../outthb.F90
outthree.o                : ../outthree.F90                ; $(F90) ../outthree.F90
outtime.o                 : ../outtime.F90                 ; $(F90) ../outtime.F90
outtwo.o                  : ../outtwo.F90                  ; $(F90) ../outtwo.F90
outvar.o                  : ../outvar.F90                  ; $(F90) ../outvar.F90
outvarg.o                 : ../outvarg.F90                 ; $(F90) ../outvarg.F90
outwarning.o              : ../outwarning.F90              ; $(F90) ../outwarning.F90
outxr.o                   : ../outxr.F90                   ; $(F90) ../outxr.F90
outxsf.o                  : ../outxsf.F90                  ; $(F90) ../outxsf.F90
outxtl.o                  : ../outxtl.F90                  ; $(F90) ../outxtl.F90
outxyz.o                  : ../outxyz.F90                  ; $(F90) ../outxyz.F90
pdiage.o                  : ../pdiage.F90                  ; $(F90) ../pdiage.F90
pdiagg.o                  : ../pdiagg.F90                  ; $(F90) ../pdiagg.F90
pdiagl.o                  : ../pdiagl.F90                  ; $(F90) ../pdiagl.F90
peigen.o                  : ../peigen.F90                  ; $(F90) ../peigen.F90
peigend.o                 : ../peigend.F90                 ; $(F90) ../peigend.F90
peigeng.o                 : ../peigeng.F90                 ; $(F90) ../peigeng.F90
peigengd.o                : ../peigengd.F90                ; $(F90) ../peigengd.F90
phase_fc.o                : ../phase_fc.F90                ; $(F90) ../phase_fc.F90
phase_fcd.o               : ../phase_fcd.F90               ; $(F90) ../phase_fcd.F90
phase_fcg.o               : ../phase_fcg.F90               ; $(F90) ../phase_fcg.F90
phoncopy1.o               : ../phoncopy1.F90               ; $(F90) ../phoncopy1.F90
phoncopy2.o               : ../phoncopy2.F90               ; $(F90) ../phoncopy2.F90
phonon.o                  : ../phonon.F90                  ; $(F90) ../phonon.F90
phonon_fc.o               : ../phonon_fc.F90               ; $(F90) ../phonon_fc.F90
phonon_fcg.o              : ../phonon_fcg.F90              ; $(F90) ../phonon_fcg.F90
phonword.o                : ../phonword.F90                ; $(F90) ../phonword.F90
pireal0d.o                : ../pireal0d.F90                ; $(F90) ../pireal0d.F90
pireal0dd.o               : ../pireal0dd.F90               ; $(F90) ../pireal0dd.F90
pirealrecip.o             : ../pirealrecip.F90             ; $(F90) ../pirealrecip.F90
pirealrecipd.o            : ../pirealrecipd.F90            ; $(F90) ../pirealrecipd.F90
pirealrecipp.o            : ../pirealrecipp.F90            ; $(F90) ../pirealrecipp.F90
pirealrecippd.o           : ../pirealrecippd.F90           ; $(F90) ../pirealrecippd.F90
pirealrecips.o            : ../pirealrecips.F90            ; $(F90) ../pirealrecips.F90
planepot2D.o              : ../planepot2D.F90              ; $(F90) ../planepot2D.F90
planepotfd.o              : ../planepotfd.F90              ; $(F90) ../planepotfd.F90
planepotmc.o              : ../planepotmc.F90              ; $(F90) ../planepotmc.F90
planepotmd.o              : ../planepotmd.F90              ; $(F90) ../planepotmd.F90
poccon.o                  : ../poccon.F90                  ; $(F90) ../poccon.F90
poccons.o                 : ../poccons.F90                 ; $(F90) ../poccons.F90
polarisation.o            : ../polarisation.F90            ; $(F90) ../polarisation.F90
potcut.o                  : ../potcut.F90                  ; $(F90O1) ../potcut.F90
potential.o               : ../potential.F90               ; $(F90) ../potential.F90
potgrid.o                 : ../potgrid.F90                 ; $(F90) ../potgrid.F90
potword21.o               : ../potword21.F90               ; $(F90O1) ../potword21.F90
potword22.o               : ../potword22.F90               ; $(F90O1) ../potword22.F90
potword3.o                : ../potword3.F90                ; $(F90O1) ../potword3.F90
potword4.o                : ../potword4.F90                ; $(F90O1) ../potword4.F90
potword6.o                : ../potword6.F90                ; $(F90O1) ../potword6.F90
potwordk.o                : ../potwordk.F90                ; $(F90O1) ../potwordk.F90
potwordm.o                : ../potwordm.F90                ; $(F90O1) ../potwordm.F90
pr_integrator.o           : ../pr_integrator.F90           ; $(F90) ../pr_integrator.F90
predict.o                 : ../predict.F90                 ; $(F90) ../predict.F90
pressure.o                : ../pressure.F90                ; $(F90) ../pressure.F90
projd0.o                  : ../projd0.F90                  ; $(F90) ../projd0.F90
projd0d.o                 : ../projd0d.F90                 ; $(F90) ../projd0d.F90
projd0r.o                 : ../projd0r.F90                 ; $(F90) ../projd0r.F90
projd3a.o                 : ../projd3a.F90                 ; $(F90) ../projd3a.F90
projd3ad.o                : ../projd3ad.F90                ; $(F90) ../projd3ad.F90
projd3dV.o                : ../projd3dV.F90                ; $(F90) ../projd3dV.F90
projd3dVd.o               : ../projd3dVd.F90               ; $(F90) ../projd3dVd.F90
projd3dVg.o               : ../projd3dVg.F90               ; $(F90) ../projd3dVg.F90
projd3dVgd.o              : ../projd3dVgd.F90              ; $(F90) ../projd3dVgd.F90
projfork0.o               : ../projfork0.F90               ; $(F90) ../projfork0.F90
projfork0r.o              : ../projfork0r.F90              ; $(F90) ../projfork0r.F90
projfork3.o               : ../projfork3.F90               ; $(F90) ../projfork3.F90
projfork3d.o              : ../projfork3d.F90              ; $(F90) ../projfork3d.F90
projfork3ddV.o            : ../projfork3ddV.F90            ; $(F90) ../projfork3ddV.F90
projfork3dm.o             : ../projfork3dm.F90             ; $(F90) ../projfork3dm.F90
projfork3dd.o             : ../projfork3dd.F90             ; $(F90) ../projfork3dd.F90
projfork3ds.o             : ../projfork3ds.F90             ; $(F90) ../projfork3ds.F90
projfork3dV.o             : ../projfork3dV.F90             ; $(F90) ../projfork3dV.F90
projfork3dVg.o            : ../projfork3dVg.F90            ; $(F90) ../projfork3dVg.F90
projfork3s.o              : ../projfork3s.F90              ; $(F90) ../projfork3s.F90
projthbk0.o               : ../projthbk0.F90               ; $(F90) ../projthbk0.F90
projthbk0r.o              : ../projthbk0r.F90              ; $(F90) ../projthbk0r.F90
projthbk3.o               : ../projthbk3.F90               ; $(F90) ../projthbk3.F90
projthbk3d.o              : ../projthbk3d.F90              ; $(F90) ../projthbk3d.F90
projthbk3ddV.o            : ../projthbk3ddV.F90            ; $(F90) ../projthbk3ddV.F90
projthbk3dm.o             : ../projthbk3dm.F90             ; $(F90) ../projthbk3dm.F90
projthbk3ds.o             : ../projthbk3ds.F90             ; $(F90) ../projthbk3ds.F90
projthbk3dV.o             : ../projthbk3dV.F90             ; $(F90) ../projthbk3dV.F90
projthbk3dVg.o            : ../projthbk3dVg.F90            ; $(F90) ../projthbk3dVg.F90
projthbk3s.o              : ../projthbk3s.F90              ; $(F90) ../projthbk3s.F90
property0.o               : ../property0.F90               ; $(F90) ../property0.F90
property12.o              : ../property12.F90              ; $(F90) ../property12.F90
property3.o               : ../property3.F90               ; $(F90) ../property3.F90
psibaskes.o               : ../psibaskes.F90               ; $(F90) ../psibaskes.F90
psirho.o                  : ../psirho.F90                  ; $(F90) ../psirho.F90
psiscreenderv.o           : ../psiscreenderv.F90           ; $(F90) ../psiscreenderv.F90
psiscreendervs.o          : ../psiscreendervs.F90          ; $(F90) ../psiscreendervs.F90
psumall.o                 : ../psumall.F90                 ; $(F90) ../psumall.F90
putpar.o                  : ../putpar.F90                  ; $(F90) ../putpar.F90
pythag.o                  : ../pythag.F90                  ; $(F90) ../pythag.F90
qeqbody.o                 : ../qeqbody.F90                 ; $(F90) ../qeqbody.F90
qmatrix1D.o               : ../qmatrix1D.F90               ; $(F90) ../qmatrix1D.F90
qmatrixelement.o          : ../qmatrixelement.F90          ; $(F90) ../qmatrixelement.F90
qmatrixelementc.o         : ../qmatrixelementc.F90         ; $(F90) ../qmatrixelementc.F90
qmatrixelemente.o         : ../qmatrixelemente.F90         ; $(F90) ../qmatrixelemente.F90
qupdate.o                 : ../qupdate.F90                 ; $(F90) ../qupdate.F90
radialforce.o             : ../radialforce.F90             ; $(F90) ../radialforce.F90
radialforcetrial.o        : ../radialforcetrial.F90        ; $(F90) ../radialforcetrial.F90
raman0.o                  : ../raman0.F90                  ; $(F90) ../raman0.F90
raman3.o                  : ../raman3.F90                  ; $(F90) ../raman3.F90
ramanstrength.o           : ../ramanstrength.F90           ; $(F90) ../ramanstrength.F90
ramantrm.o                : ../ramantrm.F90                ; $(F90) ../ramantrm.F90
random.o                  : ../random.F90                  ; $(F90) ../random.F90
ratfn.o                   : ../ratfn.F90                   ; $(F90) ../ratfn.F90
real0d.o                  : ../real0d.F90                  ; $(F90) ../real0d.F90
real0dd.o                 : ../real0dd.F90                 ; $(F90) ../real0dd.F90
real0d3.o                 : ../real0d3.F90                 ; $(F90) ../real0d3.F90
real0d3d.o                : ../real0d3d.F90                ; $(F90) ../real0d3d.F90
real1d3.o                 : ../real1d3.F90                 ; $(F90) ../real1d3.F90
real1d3dV.o               : ../real1d3dV.F90               ; $(F90) ../real1d3dV.F90
real1d3dVg.o              : ../real1d3dVg.F90              ; $(F90) ../real1d3dVg.F90
real1d3g.o                : ../real1d3g.F90                ; $(F90) ../real1d3g.F90
real1D.o                  : ../real1D.F90                  ; $(F90) ../real1D.F90
real1Dd.o                 : ../real1Dd.F90                 ; $(F90) ../real1Dd.F90
real1Dfd.o                : ../real1Dfd.F90                ; $(F90) ../real1Dfd.F90
real1Dmc.o                : ../real1Dmc.F90                ; $(F90) ../real1Dmc.F90
real1Dmd.o                : ../real1Dmd.F90                ; $(F90) ../real1Dmd.F90
real1Dp.o                 : ../real1Dp.F90                 ; $(F90) ../real1Dp.F90
real1Dpd.o                : ../real1Dpd.F90                ; $(F90) ../real1Dpd.F90
real1Dq.o                 : ../real1Dq.F90                 ; $(F90) ../real1Dq.F90
real1fc3.o                : ../real1fc3.F90                ; $(F90) ../real1fc3.F90
real11.o                  : ../real11.F90                  ; $(F90) ../real11.F90
real112.o                 : ../real112.F90                 ; $(F90) ../real112.F90
real112d.o                : ../real112d.F90                ; $(F90) ../real112d.F90
real11d.o                 : ../real11d.F90                 ; $(F90) ../real11d.F90
real11pot.o               : ../real11pot.F90               ; $(F90) ../real11pot.F90
real11s.o                 : ../real11s.F90                 ; $(F90) ../real11s.F90
real11s2.o                : ../real11s2.F90                ; $(F90) ../real11s2.F90
real12a.o                 : ../real12a.F90                 ; $(F90) ../real12a.F90
real12a2.o                : ../real12a2.F90                ; $(F90) ../real12a2.F90
real12a2d.o               : ../real12a2d.F90               ; $(F90) ../real12a2d.F90
real12ad.o                : ../real12ad.F90                ; $(F90) ../real12ad.F90
real12apot.o              : ../real12apot.F90              ; $(F90) ../real12apot.F90
real12as.o                : ../real12as.F90                ; $(F90) ../real12as.F90
real12as2.o               : ../real12as2.F90               ; $(F90) ../real12as2.F90
real2a.o                  : ../real2a.F90                  ; $(F90O1) ../real2a.F90
real2a1.o                 : ../real2a1.F90                 ; $(F90) ../real2a1.F90
real2a1d.o                : ../real2a1d.F90                ; $(F90) ../real2a1d.F90
real2a1f.o                : ../real2a1f.F90                ; $(F90) ../real2a1f.F90
real2a1fd.o               : ../real2a1fd.F90               ; $(F90) ../real2a1fd.F90
real2ad.o                 : ../real2ad.F90                 ; $(F90O1) ../real2ad.F90
real2adef.o               : ../real2adef.F90               ; $(F90) ../real2adef.F90
real2aef.o                : ../real2aef.F90                ; $(F90) ../real2aef.F90
real2amany.o              : ../real2amany.F90              ; $(F90) ../real2amany.F90
real2b.o                  : ../real2b.F90                  ; $(F90) ../real2b.F90
real2Dq.o                 : ../real2Dq.F90                 ; $(F90) ../real2Dq.F90
real4.o                   : ../real4.F90                   ; $(F90) ../real4.F90
realcmm.o                 : ../realcmm.F90                 ; $(F90) ../realcmm.F90
reale.o                   : ../reale.F90                   ; $(F90) ../reale.F90
realed.o                  : ../realed.F90                  ; $(F90) ../realed.F90
real1fc.o                 : ../real1fc.F90                 ; $(F90) ../real1fc.F90
realfc.o                  : ../realfc.F90                  ; $(F90) ../realfc.F90
realfcd.o                 : ../realfcd.F90                 ; $(F90) ../realfcd.F90
realfd0.o                 : ../realfd0.F90                 ; $(F90) ../realfd0.F90
realfd3.o                 : ../realfd3.F90                 ; $(F90) ../realfd3.F90
reallocate.o              : ../reallocate.F90              ; $(F90) ../reallocate.F90
realmc0.o                 : ../realmc0.F90                 ; $(F90) ../realmc0.F90
realmc3.o                 : ../realmc3.F90                 ; $(F90) ../realmc3.F90
realmd0.o                 : ../realmd0.F90                 ; $(F90) ../realmd0.F90
realmd0s.o                : ../realmd0s.F90                ; $(F90) ../realmd0s.F90
realmd3.o                 : ../realmd3.F90                 ; $(F90) ../realmd3.F90
realmd3s.o                : ../realmd3s.F90                ; $(F90) ../realmd3s.F90
realmd3t.o                : ../realmd3t.F90                ; $(F90) ../realmd3t.F90
realmi3.o                 : ../realmi3.F90                 ; $(F90) ../realmi3.F90
realmi3d.o                : ../realmi3d.F90                ; $(F90) ../realmi3d.F90
realmi3fd.o               : ../realmi3fd.F90               ; $(F90) ../realmi3fd.F90
realmimc3.o               : ../realmimc3.F90               ; $(F90) ../realmimc3.F90
realp.o                   : ../realp.F90                   ; $(F90) ../realp.F90
realpd.o                  : ../realpd.F90                  ; $(F90) ../realpd.F90
realrecip3d3.o            : ../realrecip3d3.F90            ; $(F90) ../realrecip3d3.F90
realrecip3d3d.o           : ../realrecip3d3d.F90           ; $(F90) ../realrecip3d3d.F90
realrecip3d3dV.o          : ../realrecip3d3dV.F90          ; $(F90) ../realrecip3d3dV.F90
realrecip3d3dVd.o         : ../realrecip3d3dVd.F90         ; $(F90) ../realrecip3d3dVd.F90
realrecip3d3dVg.o         : ../realrecip3d3dVg.F90         ; $(F90) ../realrecip3d3dVg.F90
realrecip3d3dVgd.o        : ../realrecip3d3dVgd.F90        ; $(F90) ../realrecip3d3dVgd.F90
realrecipfc3.o            : ../realrecipfc3.F90            ; $(F90) ../realrecipfc3.F90
realsd.o                  : ../realsd.F90                  ; $(F90) ../realsd.F90
realsd2.o                 : ../realsd2.F90                 ; $(F90) ../realsd2.F90
realselffc.o              : ../realselffc.F90              ; $(F90) ../realselffc.F90
realselffcd.o             : ../realselffcd.F90             ; $(F90) ../realselffcd.F90
reaxff1fc.o               : ../reaxff1fc.F90               ; $(F90) ../reaxff1fc.F90
reaxffcharges.o           : ../reaxffcharges.F90           ; $(F90) ../reaxffcharges.F90
reaxffmd.o                : ../reaxffmd.F90                ; $(F90) ../reaxffmd.F90
reaxffpot.o               : ../reaxffpot.F90               ; $(F90) ../reaxffpot.F90
reaxff_bo.o               : ../reaxff_bo.F90               ; $(F90) ../reaxff_bo.F90
reaxff_funcs.o            : ../reaxff_funcs.F90            ; $(F90) ../reaxff_funcs.F90
rebospline.o              : ../rebospline.F90              ; $(F90) ../rebospline.F90
recip2D.o                 : ../recip2D.F90                 ; $(F90) ../recip2D.F90
recip2Dd.o                : ../recip2Dd.F90                ; $(F90) ../recip2Dd.F90
recip2Dfd.o               : ../recip2Dfd.F90               ; $(F90) ../recip2Dfd.F90
recip2Dtrial.o            : ../recip2Dtrial.F90            ; $(F90) ../recip2Dtrial.F90
recip3D.o                 : ../recip3D.F90                 ; $(F90) ../recip3D.F90
recip3Dd.o                : ../recip3Dd.F90                ; $(F90) ../recip3Dd.F90
recip3Dfd.o               : ../recip3Dfd.F90               ; $(F90) ../recip3Dfd.F90
recip3Dmd.o               : ../recip3Dmd.F90               ; $(F90) ../recip3Dmd.F90
recip3Dp.o                : ../recip3Dp.F90                ; $(F90) ../recip3Dp.F90
recip3Dpd.o               : ../recip3Dpd.F90               ; $(F90) ../recip3Dpd.F90
recip3Dspme.o             : ../recip3Dspme.F90             ; $(F90) ../recip3Dspme.F90
recip3Dtrial.o            : ../recip3Dtrial.F90            ; $(F90) ../recip3Dtrial.F90
recip12a.o                : ../recip12a.F90                ; $(F90) ../recip12a.F90
recip12ad.o               : ../recip12ad.F90               ; $(F90) ../recip12ad.F90
recip12apot.o             : ../recip12apot.F90             ; $(F90) ../recip12apot.F90
recip12as.o               : ../recip12as.F90               ; $(F90) ../recip12as.F90
recip4.o                  : ../recip4.F90                  ; $(F90) ../recip4.F90
recip2Dp.o                : ../recip2Dp.F90                ; $(F90) ../recip2Dp.F90
recip2Dpd.o               : ../recip2Dpd.F90               ; $(F90) ../recip2Dpd.F90
recip2Dq.o                : ../recip2Dq.F90                ; $(F90) ../recip2Dq.F90
recipfc3.o                : ../recipfc3.F90                ; $(F90) ../recipfc3.F90
recipsd.o                 : ../recipsd.F90                 ; $(F90) ../recipsd.F90
recipsd2.o                : ../recipsd2.F90                ; $(F90) ../recipsd2.F90
reciptrmd3.o              : ../reciptrmd3.F90              ; $(F90) ../reciptrmd3.F90
reciptrmd3dV.o            : ../reciptrmd3dV.F90            ; $(F90) ../reciptrmd3dV.F90
reciptrmd3dVg.o           : ../reciptrmd3dVg.F90           ; $(F90) ../reciptrmd3dVg.F90
reciptrmd3g.o             : ../reciptrmd3g.F90             ; $(F90) ../reciptrmd3g.F90
reciptrmdp.o              : ../reciptrmdp.F90              ; $(F90) ../reciptrmdp.F90
reinitialise.o            : ../reinitialise.F90            ; $(F90) ../reinitialise.F90
remove.o                  : ../remove.F90                  ; $(F90) ../remove.F90
reorder.o                 : ../reorder.F90                 ; $(F90) ../reorder.F90
resetcfgdefaults.o        : ../resetcfgdefaults.F90        ; $(F90) ../resetcfgdefaults.F90
resetrigidmol.o           : ../resetrigidmol.F90           ; $(F90) ../resetrigidmol.F90
resetx0rv.o               : ../resetx0rv.F90               ; $(F90) ../resetx0rv.F90
rfind.o                   : ../rfind.F90                   ; $(F90) ../rfind.F90
rfindboth.o               : ../rfindboth.F90               ; $(F90) ../rfindboth.F90
rfohess.o                 : ../rfohess.F90                 ; $(F90) ../rfohess.F90
rhoinv.o                  : ../rhoinv.F90                  ; $(F90) ../rhoinv.F90
rhtohex.o                 : ../rhtohex.F90                 ; $(F90) ../rhtohex.F90
rigidmoleculedrv.o        : ../rigidmoleculedrv.F90        ; $(F90) ../rigidmoleculedrv.F90
rigidmoleculegv.o         : ../rigidmoleculegv.F90         ; $(F90) ../rigidmoleculegv.F90
rigidmoleculephon.o       : ../rigidmoleculephon.F90       ; $(F90) ../rigidmoleculephon.F90
rigidmoleculeprop.o       : ../rigidmoleculeprop.F90       ; $(F90) ../rigidmoleculeprop.F90
rlist.o                   : ../rlist.F90                   ; $(F90) ../rlist.F90
rpair.o                   : ../rpair.F90                   ; $(F90) ../rpair.F90
rpair2.o                  : ../rpair2.F90                  ; $(F90) ../rpair2.F90
rsearch1D.o               : ../rsearch1D.F90               ; $(F90) ../rsearch1D.F90
rsearch2D.o               : ../rsearch2D.F90               ; $(F90) ../rsearch2D.F90
rsearch3D.o               : ../rsearch3D.F90               ; $(F90) ../rsearch3D.F90
rtclist.o                 : ../rtclist.F90                 ; $(F90) ../rtclist.F90
rtlist.o                  : ../rtlist.F90                  ; $(F90) ../rtlist.F90
runmc.o                   : ../runmc.F90                   ; $(F90) ../runmc.F90
runneb.o                  : ../runneb.F90                  ; $(F90) ../runneb.F90
runsync.o                 : ../runsync.F90                 ; $(F90) ../runsync.F90
samemol.o                 : ../samemol.F90                 ; $(F90) ../samemol.F90
scaleevec.o               : ../scaleevec.F90               ; $(F90) ../scaleevec.F90
sec0.o                    : ../sec0.F90                    ; $(F90) ../sec0.F90
sec0f.o                   : ../sec0f.F90                   ; $(F90) ../sec0f.F90
sec3.o                    : ../sec3.F90                    ; $(F90) ../sec3.F90
sec3f.o                   : ../sec3f.F90                   ; $(F90) ../sec3f.F90
secondorderfc.o           : ../secondorderfc.F90           ; $(F90) ../secondorderfc.F90
secondpass.o              : ../secondpass.F90              ; $(F90) ../secondpass.F90
selfterm.o                : ../selfterm.F90                ; $(F90) ../selfterm.F90
selftermfc.o              : ../selftermfc.F90              ; $(F90) ../selftermfc.F90
selftermp.o               : ../selftermp.F90               ; $(F90) ../selftermp.F90
selftermfcd.o             : ../selftermfcd.F90             ; $(F90) ../selftermfcd.F90
set2a.o                   : ../set2a.F90                   ; $(F90) ../set2a.F90
setatomdistribution.o     : ../setatomdistribution.F90     ; $(F90) ../setatomdistribution.F90
setatomnoded2.o           : ../setatomnoded2.F90           ; $(F90) ../setatomnoded2.F90
setatomnodes.o            : ../setatomnodes.F90            ; $(F90) ../setatomnodes.F90
setatomnodesbo.o          : ../setatomnodesbo.F90          ; $(F90) ../setatomnodesbo.F90
setbelashchenko.o         : ../setbelashchenko.F90         ; $(F90) ../setbelashchenko.F90
setbondorder.o            : ../setbondorder.F90            ; $(F90) ../setbondorder.F90
setbondtype.o             : ../setbondtype.F90             ; $(F90) ../setbondtype.F90
setbrenner1.o             : ../setbrenner1.F90             ; $(F90) ../setbrenner1.F90
setbrenner3.o             : ../setbrenner3.F90             ; $(F90) ../setbrenner3.F90
setbsmptr.o               : ../setbsmptr.F90               ; $(F90) ../setbsmptr.F90
setc6.o                   : ../setc6.F90                   ; $(F90) ../setc6.F90
setcellderv.o             : ../setcellderv.F90             ; $(F90) ../setcellderv.F90
setcfg.o                  : ../setcfg.F90                  ; $(F90) ../setcfg.F90
setcmm.o                  : ../setcmm.F90                  ; $(F90) ../setcmm.F90
setcosmo.o                : ../setcosmo.F90                ; $(F90) ../setcosmo.F90
setcosmoainv.o            : ../setcosmoainv.F90            ; $(F90) ../setcosmoainv.F90
setcosmoamat.o            : ../setcosmoamat.F90            ; $(F90) ../setcosmoamat.F90
setcosmoamat2D.o          : ../setcosmoamat2D.F90          ; $(F90) ../setcosmoamat2D.F90
setcosmobmat.o            : ../setcosmobmat.F90            ; $(F90) ../setcosmobmat.F90
setdcosmoamat.o           : ../setdcosmoamat.F90           ; $(F90) ../setdcosmoamat.F90
setdcosmoamat1p.o         : ../setdcosmoamat1p.F90         ; $(F90) ../setdcosmoamat1p.F90
setdcosmobmat.o           : ../setdcosmobmat.F90           ; $(F90) ../setdcosmobmat.F90
setdcosmobmat1p.o         : ../setdcosmobmat1p.F90         ; $(F90) ../setdcosmobmat1p.F90
setcutoffmaxbo.o          : ../setcutoffmaxbo.F90          ; $(F90) ../setcutoffmaxbo.F90
setdef.o                  : ../setdef.F90                  ; $(F90) ../setdef.F90
setdestroyptr.o           : ../setdestroyptr.F90           ; $(F90) ../setdestroyptr.F90
setdisp.o                 : ../setdisp.F90                 ; $(F90) ../setdisp.F90
setdistance.o             : ../setdistance.F90             ; $(F90) ../setdistance.F90
setdr2mat.o               : ../setdr2mat.F90               ; $(F90O1) ../setdr2mat.F90
seteam.o                  : ../seteam.F90                  ; $(F90) ../seteam.F90
setedip.o                 : ../setedip.F90                 ; $(F90) ../setedip.F90
seteem.o                  : ../seteem.F90                  ; $(F90) ../seteem.F90
seteembond.o              : ../seteembond.F90              ; $(F90) ../seteembond.F90
setele.o                  : ../setele.F90                  ; $(F90) ../setele.F90
setewald.o                : ../setewald.F90                ; $(F90) ../setewald.F90
setfeshtmat.o             : ../setfeshtmat.F90             ; $(F90) ../setfeshtmat.F90
setfour.o                 : ../setfour.F90                 ; $(F90) ../setfour.F90
setframe.o                : ../setframe.F90                ; $(F90) ../setframe.F90
setkeywrd.o               : ../setkeywrd.F90               ; $(F90) ../setkeywrd.F90
setkpt.o                  : ../setkpt.F90                  ; $(F90) ../setkpt.F90
setktrmd3.o               : ../setktrmd3.F90               ; $(F90) ../setktrmd3.F90
setktrmdp.o               : ../setktrmdp.F90               ; $(F90) ../setktrmdp.F90
setktrm4.o                : ../setktrm4.F90                ; $(F90) ../setktrm4.F90
setktrmfc3.o              : ../setktrmfc3.F90              ; $(F90) ../setktrmfc3.F90
setcombine.o              : ../setcombine.F90              ; $(F90) ../setcombine.F90
setlib.o                  : ../setlib.F90                  ; $(F90) ../setlib.F90
setlist3.o                : ../setlist3.F90                ; $(F90) ../setlist3.F90
setlist4.o                : ../setlist4.F90                ; $(F90) ../setlist4.F90
setlist6.o                : ../setlist6.F90                ; $(F90) ../setlist6.F90
setlistimp.o              : ../setlistimp.F90              ; $(F90) ../setlistimp.F90
setlistoop.o              : ../setlistoop.F90              ; $(F90) ../setlistoop.F90
setmaxcell1D.o            : ../setmaxcell1D.F90            ; $(F90) ../setmaxcell1D.F90
setmd.o                   : ../setmd.F90                   ; $(F90) ../setmd.F90
setmol.o                  : ../setmol.F90                  ; $(F90) ../setmol.F90
setmola.o                 : ../setmola.F90                 ; $(F90) ../setmola.F90
setmols.o                 : ../setmols.F90                 ; $(F90) ../setmols.F90
setmoldef.o               : ../setmoldef.F90               ; $(F90) ../setmoldef.F90
setmoldef12a.o            : ../setmoldef12a.F90            ; $(F90) ../setmoldef12a.F90
setmoveptr.o              : ../setmoveptr.F90              ; $(F90) ../setmoveptr.F90
setoccptr.o               : ../setoccptr.F90               ; $(F90) ../setoccptr.F90
setoccptrs.o              : ../setoccptrs.F90              ; $(F90) ../setoccptrs.F90
setoccshptr.o             : ../setoccshptr.F90             ; $(F90) ../setoccshptr.F90
setoccshptrs.o            : ../setoccshptrs.F90            ; $(F90) ../setoccshptrs.F90
setoptptr.o               : ../setoptptr.F90               ; $(F90) ../setoptptr.F90
setphonptr.o              : ../setphonptr.F90              ; $(F90) ../setphonptr.F90
setpolar.o                : ../setpolar.F90                ; $(F90) ../setpolar.F90
setpote.o                 : ../setpote.F90                 ; $(F90) ../setpote.F90
setratiom.o               : ../setratiom.F90               ; $(F90) ../setratiom.F90
setreaxFF.o               : ../setreaxFF.F90               ; $(F90) ../setreaxFF.F90
setreaxffQ.o              : ../setreaxffQ.F90              ; $(F90) ../setreaxffQ.F90
setreaxffQ1fc.o           : ../setreaxffQ1fc.F90           ; $(F90) ../setreaxffQ1fc.F90
setreaxffQiter.o          : ../setreaxffQiter.F90          ; $(F90) ../setreaxffQiter.F90
setreaxffQiter1fc.o       : ../setreaxffQiter1fc.F90       ; $(F90) ../setreaxffQiter1fc.F90
setregiontrans.o          : ../setregiontrans.F90          ; $(F90) ../setregiontrans.F90
setrigidmol.o             : ../setrigidmol.F90             ; $(F90) ../setrigidmol.F90
setrotateptr.o            : ../setrotateptr.F90            ; $(F90) ../setrotateptr.F90
setsas.o                  : ../setsas.F90                  ; $(F90) ../setsas.F90
setsegweight.o            : ../setsegweight.F90            ; $(F90) ../setsegweight.F90
setdsegweight.o           : ../setdsegweight.F90           ; $(F90) ../setdsegweight.F90
setptsnoded2.o            : ../setptsnoded2.F90            ; $(F90) ../setptsnoded2.F90
setqonsas.o               : ../setqonsas.F90               ; $(F90) ../setqonsas.F90
setstrainptr.o            : ../setstrainptr.F90            ; $(F90) ../setstrainptr.F90
setsix.o                  : ../setsix.F90                  ; $(F90) ../setsix.F90
setslice.o                : ../setslice.F90                ; $(F90) ../setslice.F90
setspatial.o              : ../setspatial.F90              ; $(F90) ../setspatial.F90
setspatialbo.o            : ../setspatialbo.F90            ; $(F90) ../setspatialbo.F90
setspatialcells.o         : ../setspatialcells.F90         ; $(F90) ../setspatialcells.F90
setspec.o                 : ../setspec.F90                 ; $(F90) ../setspec.F90
setspring.o               : ../setspring.F90               ; $(F90) ../setspring.F90
setswapptr.o              : ../setswapptr.F90              ; $(F90) ../setswapptr.F90
settaper.o                : ../settaper.F90                ; $(F90) ../settaper.F90
setthree.o                : ../setthree.F90                ; $(F90) ../setthree.F90
settings.o                : ../settings.F90                ; $(F90) ../settings.F90
settorcubic.o             : ../settorcubic.F90             ; $(F90O1) ../settorcubic.F90
settricubic.o             : ../settricubic.F90             ; $(F90O1) ../settricubic.F90
settricubic2.o            : ../settricubic2.F90            ; $(F90O1) ../settricubic2.F90
setuff.o                  : ../setuff.F90                  ; $(F90) ../setuff.F90
setup.o                   : ../setup.F90                   ; $(F90) ../setup.F90
setvarnoded2.o            : ../setvarnoded2.F90            ; $(F90) ../setvarnoded2.F90
setvbon.o                 : ../setvbon.F90                 ; $(F90) ../setvbon.F90
setzmolslice.o            : ../setzmolslice.F90            ; $(F90) ../setzmolslice.F90
shellgroupvelocities.o    : ../shellgroupvelocities.F90    ; $(F90) ../shellgroupvelocities.F90
shellopt.o                : ../shellopt.F90                ; $(F90) ../shellopt.F90
shmin.o                   : ../shmin.F90                   ; $(F90) ../shmin.F90
sitepot.o                 : ../sitepot.F90                 ; $(F90) ../sitepot.F90
sixbody.o                 : ../sixbody.F90                 ; $(F90) ../sixbody.F90
six1fc.o                  : ../six1fc.F90                  ; $(F90) ../six1fc.F90
sixfc.o                   : ../sixfc.F90                   ; $(F90) ../sixfc.F90
sixfc3.o                  : ../sixfc3.F90                  ; $(F90) ../sixfc3.F90
sixfcd.o                  : ../sixfcd.F90                  ; $(F90) ../sixfcd.F90
sixfd.o                   : ../sixfd.F90                   ; $(F90) ../sixfd.F90
sixlist.o                 : ../sixlist.F90                 ; $(F90) ../sixlist.F90
sixmc.o                   : ../sixmc.F90                   ; $(F90) ../sixmc.F90
sixmd.o                   : ../sixmd.F90                   ; $(F90) ../sixmd.F90
sixmds.o                  : ../sixmds.F90                  ; $(F90) ../sixmds.F90
sixnos.o                  : ../sixnos.F90                  ; $(F90) ../sixnos.F90
sixnosd.o                 : ../sixnosd.F90                 ; $(F90) ../sixnosd.F90
sixp.o                    : ../sixp.F90                    ; $(F90) ../sixp.F90
sixpd.o                   : ../sixpd.F90                   ; $(F90) ../sixpd.F90
sixsd.o                   : ../sixsd.F90                   ; $(F90) ../sixsd.F90
sixsd2.o                  : ../sixsd2.F90                  ; $(F90) ../sixsd2.F90
smbody.o                  : ../smbody.F90                  ; $(F90) ../smbody.F90
solvation.o               : ../solvation.F90               ; $(F90) ../solvation.F90
sort.o                    : ../sort.F90                    ; $(F90) ../sort.F90
sort1.o                   : ../sort1.F90                   ; $(F90) ../sort1.F90
sort1sym.o                : ../sort1sym.F90                ; $(F90) ../sort1sym.F90
sort2a.o                  : ../sort2a.F90                  ; $(F90) ../sort2a.F90
special.o                 : ../special.F90                 ; $(F90) ../special.F90
specialcom.o              : ../specialcom.F90              ; $(F90) ../specialcom.F90
specialqtn.o              : ../specialqtn.F90              ; $(F90) ../specialqtn.F90
speclibw.o                : ../speclibw.F90                ; $(F90) ../speclibw.F90
specword.o                : ../specword.F90                ; $(F90) ../specword.F90
spline.o                  : ../spline.F90                  ; $(F90) ../spline.F90
stolc.o                   : ../stolc.F90                   ; $(F90) ../stolc.F90
stouc.o                   : ../stouc.F90                   ; $(F90) ../stouc.F90
strain.o                  : ../strain.F90                  ; $(F90) ../strain.F90
strfin.o                  : ../strfin.F90                  ; $(F90) ../strfin.F90
strsym.o                  : ../strsym.F90                  ; $(F90) ../strsym.F90
strword.o                 : ../strword.F90                 ; $(F90) ../strword.F90
sumderv2.o                : ../sumderv2.F90                ; $(F90) ../sumderv2.F90
sumderv2f.o               : ../sumderv2f.F90               ; $(F90) ../sumderv2f.F90
sumderv2p.o               : ../sumderv2p.F90               ; $(F90) ../sumderv2p.F90
sumderv2s.o               : ../sumderv2s.F90               ; $(F90) ../sumderv2s.F90
super.o                   : ../super.F90                   ; $(F90) ../super.F90
surfaceenergy.o           : ../surfaceenergy.F90           ; $(F90) ../surfaceenergy.F90
surfword.o                : ../surfword.F90                ; $(F90) ../surfword.F90
symbind.o                 : ../symbind.F90                 ; $(F90) ../symbind.F90
symderv2.o                : ../symderv2.F90                ; $(F90) ../symderv2.F90
symdervrot.o              : ../symdervrot.F90              ; $(F90) ../symdervrot.F90
symgenkpt.o               : ../symgenkpt.F90               ; $(F90) ../symgenkpt.F90
symmet.o                  : ../symmet.F90                  ; $(F90) ../symmet.F90
symmetp.o                 : ../symmetp.F90                 ; $(F90) ../symmetp.F90
symoff.o                  : ../symoff.F90                  ; $(F90) ../symoff.F90
symout.o                  : ../symout.F90                  ; $(F90) ../symout.F90
sympot.o                  : ../sympot.F90                  ; $(F90) ../sympot.F90
sympotd.o                 : ../sympotd.F90                 ; $(F90) ../sympotd.F90
symprod.o                 : ../symprod.F90                 ; $(F90) ../symprod.F90
syncforce.o               : ../syncforce.F90               ; $(F90) ../syncforce.F90
switch.o                  : ../switch.F90                  ; $(F90) ../switch.F90
tangtoennies.o            : ../tangtoennies.F90            ; $(F90) ../tangtoennies.F90
taper.o                   : ../taper.F90                   ; $(F90) ../taper.F90
testspacegroup.o          : ../testspacegroup.F90          ; $(F90) ../testspacegroup.F90
tfunc.o                   : ../tfunc.F90                   ; $(F90) ../tfunc.F90
thermalconductivity.o     : ../thermalconductivity.F90     ; $(F90) ../thermalconductivity.F90
thirdorderfc3.o           : ../thirdorderfc3.F90           ; $(F90) ../thirdorderfc3.F90
thirdordern.o             : ../thirdordern.F90             ; $(F90) ../thirdordern.F90
three1fc.o                : ../three1fc.F90                ; $(F90) ../three1fc.F90
threefc.o                 : ../threefc.F90                 ; $(F90) ../threefc.F90
threefd.o                 : ../threefd.F90                 ; $(F90) ../threefd.F90
threefcd.o                : ../threefcd.F90                ; $(F90) ../threefcd.F90
threenos.o                : ../threenos.F90                ; $(F90) ../threenos.F90
threenosd.o               : ../threenosd.F90               ; $(F90) ../threenosd.F90
three12.o                 : ../three12.F90                 ; $(F90) ../three12.F90
three12d.o                : ../three12d.F90                ; $(F90) ../three12d.F90
threebody.o               : ../threebody.F90               ; $(F90O1) ../threebody.F90
threelist.o               : ../threelist.F90               ; $(F90) ../threelist.F90
threep.o                  : ../threep.F90                  ; $(F90) ../threep.F90
threepd.o                 : ../threepd.F90                 ; $(F90) ../threepd.F90
threemc.o                 : ../threemc.F90                 ; $(F90) ../threemc.F90
threemd.o                 : ../threemd.F90                 ; $(F90) ../threemd.F90
threemds.o                : ../threemds.F90                ; $(F90) ../threemds.F90
threesd.o                 : ../threesd.F90                 ; $(F90) ../threesd.F90
threesd2.o                : ../threesd2.F90                ; $(F90) ../threesd2.F90
three0d3.o                : ../three0d3.F90                ; $(F90) ../three0d3.F90
three3d3.o                : ../three3d3.F90                ; $(F90) ../three3d3.F90
three3d3dV.o              : ../three3d3dV.F90              ; $(F90) ../three3d3dV.F90
three3d3dVg.o             : ../three3d3dVg.F90             ; $(F90) ../three3d3dVg.F90
three3d3g.o               : ../three3d3g.F90               ; $(F90) ../three3d3g.F90
threefc3.o                : ../threefc3.F90                ; $(F90) ../threefc3.F90
transmat.o                : ../transmat.F90                ; $(F90) ../transmat.F90
transmatd.o               : ../transmatd.F90               ; $(F90) ../transmatd.F90
transmatmol.o             : ../transmatmol.F90             ; $(F90) ../transmatmol.F90
trap_sigint.o             : ../trap_sigint.F90             ; $(F90) ../trap_sigint.F90
twobody.o                 : ../twobody.F90                 ; $(F90) ../twobody.F90
twobodymd.o               : ../twobodymd.F90               ; $(F90) ../twobodymd.F90
twoden.o                  : ../twoden.F90                  ; $(F90) ../twoden.F90
uncell.o                  : ../uncell.F90                  ; $(F90) ../uncell.F90
uncentre.o                : ../uncentre.F90                ; $(F90) ../uncentre.F90
unfreeze.o                : ../unfreeze.F90                ; $(F90) ../unfreeze.F90
updatesas.o               : ../updatesas.F90               ; $(F90) ../updatesas.F90
validtriad.o              : ../validtriad.F90              ; $(F90) ../validtriad.F90
vartocfg.o                : ../vartocfg.F90                ; $(F90) ../vartocfg.F90
vibenergy.o               : ../vibenergy.F90               ; $(F90) ../vibenergy.F90
volume.o                  : ../volume.F90                  ; $(F90) ../volume.F90
wcoord.o                  : ../wcoord.F90                  ; $(F90) ../wcoord.F90
wolfself.o                : ../wolfself.F90                ; $(F90) ../wolfself.F90
wolfselftrial.o           : ../wolfselftrial.F90           ; $(F90) ../wolfselftrial.F90
worsy.o                   : ../worsy.F90                   ; $(F90) ../worsy.F90
wtof.o                    : ../wtof.F90                    ; $(F90) ../wtof.F90
x0strain.o                : ../x0strain.F90                ; $(F90) ../x0strain.F90
x0tostr.o                 : ../x0tostr.F90                 ; $(F90) ../x0tostr.F90
x0tostrcentroid.o         : ../x0tostrcentroid.F90         ; $(F90) ../x0tostrcentroid.F90
x0tostrtrial.o            : ../x0tostrtrial.F90            ; $(F90) ../x0tostrtrial.F90
xctox0.o                  : ../xctox0.F90                  ; $(F90) ../xctox0.F90
#
qtpie.o                   : ../qtpie.F90                   ; $(F90) ../qtpie.F90
qtmin.o                   : ../qtmin.F90                   ; $(F90) ../qtmin.F90
ctfunct.o                 : ../ctfunct.F90                 ; $(F90) ../ctfunct.F90
reaxffctfunct.o           : ../reaxffctfunct.F90           ; $(F90) ../reaxffctfunct.F90
m_conjgr.o                : ../m_conjgr.F90                ; $(F90) ../m_conjgr.F90
screenct.o                : ../screenct.F90                ; $(F90) ../screenct.F90
#
gulpcml.o                 : ../gulpcml.F90 $(FOX_OBJ)                  ; $(F90) $(FoX_fcflags) ../gulpcml.F90
gulpcml_md.o              : gulpcml.o ../gulpcml_md.F90 $(FOX_OBJ)     ; $(F90) $(FoX_fcflags) ../gulpcml_md.F90
gulpcml_phonon.o          : gulpcml.o ../gulpcml_phonon.F90 $(FOX_OBJ) ; $(F90) $(FoX_fcflags) ../gulpcml_phonon.F90
gulpcml_props.o           : gulpcml.o ../gulpcml_props.F90 $(FOX_OBJ)  ; $(F90) $(FoX_fcflags) ../gulpcml_props.F90
#
#  C routines
#
csignal.o                 : ../csignal.c       	           ; $(CC) ../csignal.c
#
#  Automatically generated dependences from mkmf
#
CS_dummy.o: ../CS_dummy.F90 datatypes.o
addframe2arc.o: ../addframe2arc.F90 modules.o m_trace.o
addshell.o: ../addshell.F90 modules.o iomod.o m_trace.o
aocn2.o: ../aocn2.F90 modules.o iomod.o m_trace.o
applyconcfg.o: ../applyconcfg.F90 modules.o datatypes.o reallocate.o m_trace.o
applyconstraints.o: ../applyconstraints.F90 modules.o datatypes.o reallocate.o m_trace.o
background.o: ../background.F90 modules.o m_trace.o m_strain.o
banner.o: ../banner.F90 iomod.o modules.o m_trace.o
baskes.o: ../baskes.F90 modules.o m_trace.o
bcgsolve.o: ../bcgsolve.F90 datatypes.o iomod.o modules.o m_trace.o
binary.o: ../binary.F90 modules.o m_trace.o datatypes.o
blas.o: ../blas.F90 datatypes.o
bond.o: ../bond.F90 modules.o iomod.o m_trace.o
bonded.o: ../bonded.F90 modules.o m_trace.o
bonded2donor.o: ../bonded2donor.F90 modules.o m_trace.o
bonded3.o: ../bonded3.F90 modules.o m_trace.o
bonded3c.o: ../bonded3c.F90 modules.o m_trace.o
bondint.o: ../bondint.F90 modules.o m_trace.o
bondorder.o: ../bondorder.F90 datatypes.o modules.o iomod.o m_trace.o
bondorder12.o: ../bondorder12.F90 datatypes.o modules.o iomod.o m_trace.o
bondorder12d.o: ../bondorder12d.F90 datatypes.o modules.o iomod.o m_trace.o
bondorder1fc.o: ../bondorder1fc.F90 datatypes.o modules.o iomod.o m_trace.o
bondorderd.o: ../bondorderd.F90 datatypes.o modules.o iomod.o m_trace.o
bondorderfc.o: ../bondorderfc.F90 datatypes.o modules.o iomod.o m_trace.o
bondorderfcd.o: ../bondorderfcd.F90 datatypes.o modules.o iomod.o m_trace.o
bondordermd.o: ../bondordermd.F90 datatypes.o modules.o iomod.o m_trace.o
bondordersd2.o: ../bondordersd2.F90 datatypes.o modules.o iomod.o m_trace.o
bondpot.o: ../bondpot.F90 modules.o m_trace.o
bondq.o: ../bondq.F90 modules.o iomod.o m_trace.o
bondqtrial.o: ../bondqtrial.F90 modules.o iomod.o m_trace.o
borncharge.o: ../borncharge.F90 modules.o gulpcml.o gulpcml_props.o iomod.o m_trace.o
boself.o: ../boself.F90 datatypes.o modules.o m_trace.o
boselfd.o: ../boselfd.F90 datatypes.o modules.o m_trace.o
boselffc.o: ../boselffc.F90 datatypes.o modules.o m_trace.o
boselfp.o: ../boselfp.F90 datatypes.o modules.o m_trace.o
boselftrial.o: ../boselftrial.F90 datatypes.o modules.o m_trace.o
boword.o: ../boword.F90 modules.o iomod.o
brenner.o: ../brenner.F90 datatypes.o modules.o iomod.o m_trace.o
brenner12.o: ../brenner12.F90 datatypes.o modules.o iomod.o m_trace.o
brenner12d.o: ../brenner12d.F90 datatypes.o modules.o iomod.o m_trace.o
brenner1fc.o: ../brenner1fc.F90 datatypes.o modules.o iomod.o m_trace.o
brennerd.o: ../brennerd.F90 datatypes.o modules.o iomod.o m_trace.o
brennerfc.o: ../brennerfc.F90 datatypes.o modules.o iomod.o m_trace.o
brennerfcd.o: ../brennerfcd.F90 datatypes.o modules.o iomod.o m_trace.o
brennermd.o: ../brennermd.F90 datatypes.o modules.o iomod.o m_trace.o
brennersd2.o: ../brennersd2.F90 datatypes.o modules.o iomod.o m_trace.o
bspline.o: ../bspline.F90 datatypes.o m_trace.o
buck4.o: ../buck4.F90 modules.o iomod.o m_trace.o
calcf2D.o: ../calcf2D.F90 modules.o m_trace.o datatypes.o
calcpht.o: ../calcpht.F90 modules.o m_trace.o datatypes.o
cart2frac.o: ../cart2frac.F90 datatypes.o modules.o m_trace.o
cdabs.o: ../cdabs.F90 datatypes.o
cell.o: ../cell.F90 modules.o m_trace.o
cellcheck.o: ../cellcheck.F90 datatypes.o m_trace.o
celldrv.o: ../celldrv.F90 modules.o m_trace.o
cellfhr.o: ../cellfhr.F90 datatypes.o m_trace.o
celltype.o: ../celltype.F90 datatypes.o modules.o m_trace.o
centre.o: ../centre.F90 modules.o m_trace.o
changemaxallnearseg.o: ../changemaxallnearseg.F90 modules.o reallocate.o
changemaxat.o: ../changemaxat.F90 modules.o reallocate.o
changemaxatcosmo.o: ../changemaxatcosmo.F90 modules.o reallocate.o
changemaxatloc.o: ../changemaxatloc.F90 modules.o reallocate.o
changemaxatompernode.o: ../changemaxatompernode.F90 modules.o reallocate.o
changemaxatompernodebo.o: ../changemaxatompernodebo.F90 modules.o reallocate.o
changemaxatot.o: ../changemaxatot.F90 modules.o reallocate.o
changemaxbond.o: ../changemaxbond.F90 modules.o reallocate.o
changemaxbondq.o: ../changemaxbondq.F90 modules.o reallocate.o
changemaxbondtype.o: ../changemaxbondtype.F90 modules.o reallocate.o
changemaxbondvec.o: ../changemaxbondvec.F90 modules.o reallocate.o
changemaxccspec.o: ../changemaxccspec.F90 modules.o reallocate.o
changemaxcellpernode.o: ../changemaxcellpernode.F90 modules.o reallocate.o
changemaxcellpernodebo.o: ../changemaxcellpernodebo.F90 modules.o reallocate.o
changemaxcfg.o: ../changemaxcfg.F90 modules.o modules_kim.o m_pr.o m_pdfneutron.o reallocate.o iomod.o
changemaxcon.o: ../changemaxcon.F90 modules.o reallocate.o
changemaxconnect.o: ../changemaxconnect.F90 modules.o reallocate.o
changemaxconnectpermol.o: ../changemaxconnectpermol.F90 modules.o reallocate.o
changemaxcontot.o: ../changemaxcontot.F90 modules.o reallocate.o
changemaxd1.o: ../changemaxd1.F90 modules.o reallocate.o
changemaxd2.o: ../changemaxd2.F90 modules.o reallocate.o
changemaxd2cell.o: ../changemaxd2cell.F90 modules.o reallocate.o
changemaxd2cells.o: ../changemaxd2cells.F90 modules.o reallocate.o
changemaxd2q.o: ../changemaxd2q.F90 modules.o reallocate.o
changemaxdcon.o: ../changemaxdcon.F90 modules.o reallocate.o
changemaxdef.o: ../changemaxdef.F90 modules.o reallocate.o
changemaxdis.o: ../changemaxdis.F90 modules.o reallocate.o
changemaxeamden.o: ../changemaxeamden.F90 modules.o reallocate.o
changemaxeamfnspec.o: ../changemaxeamfnspec.F90 modules.o reallocate.o iomod.o
changemaxeamspec.o: ../changemaxeamspec.F90 modules.o reallocate.o iomod.o
changemaxedipspec.o: ../changemaxedipspec.F90 modules.o reallocate.o
changemaxeembond.o: ../changemaxeembond.F90 modules.o reallocate.o
changemaxextrapol.o: ../changemaxextrapol.F90 modules.o reallocate.o
changemaxfcon.o: ../changemaxfcon.F90 modules.o reallocate.o
changemaxfgrad.o: ../changemaxfgrad.F90 modules.o reallocate.o
changemaxfit.o: ../changemaxfit.F90 modules.o reallocate.o
changemaxfkpt.o: ../changemaxfkpt.F90 modules.o reallocate.o
changemaxfor.o: ../changemaxfor.F90 modules.o reallocate.o
changemaxfreq.o: ../changemaxfreq.F90 datatypes.o modules.o reallocate.o
changemaxfstrain.o: ../changemaxfstrain.F90 modules.o reallocate.o
changemaxgbcfg.o: ../changemaxgbcfg.F90 modules.o reallocate.o
changemaxgcmcmol.o: ../changemaxgcmcmol.F90 modules.o reallocate.o
changemaxhold.o: ../changemaxhold.F90 reallocate.o modules.o
changemaxiltor.o: ../changemaxiltor.F90 modules.o reallocate.o
changemaxint.o: ../changemaxint.F90 modules.o reallocate.o
changemaxkpt.o: ../changemaxkpt.F90 modules.o reallocate.o
changemaxkvec.o: ../changemaxkvec.F90 modules.o reallocate.o
changemaxlib.o: ../changemaxlib.F90 modules.o reallocate.o iomod.o
changemaxlist3.o: ../changemaxlist3.F90 reallocate.o modules.o
changemaxlist4.o: ../changemaxlist4.F90 modules.o reallocate.o
changemaxlist6.o: ../changemaxlist6.F90 modules.o reallocate.o
changemaxmany.o: ../changemaxmany.F90 modules.o reallocate.o
changemaxmcswaps.o: ../changemaxmcswaps.F90 modules.o reallocate.o
changemaxmcswapspec.o: ../changemaxmcswapspec.F90 modules.o reallocate.o
changemaxmol.o: ../changemaxmol.F90 modules.o reallocate.o
changemaxmolat.o: ../changemaxmolat.F90 modules.o reallocate.o
changemaxmoleqv.o: ../changemaxmoleqv.F90 modules.o reallocate.o
changemaxn3a.o: ../changemaxn3a.F90 reallocate.o modules.o
changemaxn3bondnono.o: ../changemaxn3bondnono.F90 reallocate.o modules.o
changemaxn3f.o: ../changemaxn3f.F90 reallocate.o modules.o
changemaxn3ma.o: ../changemaxn3ma.F90 reallocate.o modules.o
changemaxn3mf.o: ../changemaxn3mf.F90 reallocate.o modules.o
changemaxnboa.o: ../changemaxnboa.F90 modules.o reallocate.o
changemaxnboo.o: ../changemaxnboo.F90 modules.o reallocate.o
changemaxnbopot.o: ../changemaxnbopot.F90 modules.o reallocate.o
changemaxnboq.o: ../changemaxnboq.F90 modules.o reallocate.o
changemaxnboq0.o: ../changemaxnboq0.F90 modules.o reallocate.o
changemaxnbor.o: ../changemaxnbor.F90 modules.o reallocate.o
changemaxnboz.o: ../changemaxnboz.F90 modules.o reallocate.o
changemaxndistancetotal.o: ../changemaxndistancetotal.F90 modules.o reallocate.o
changemaxneamfnnumeric.o: ../changemaxneamfnnumeric.F90 modules.o reallocate.o
changemaxnearseg.o: ../changemaxnearseg.F90 modules.o reallocate.o
changemaxnebreplicatot.o: ../changemaxnebreplicatot.F90 modules.o reallocate.o
changemaxnobo.o: ../changemaxnobo.F90 modules.o reallocate.o
changemaxnppa.o: ../changemaxnppa.F90 modules.o reallocate.o
changemaxnpts.o: ../changemaxnpts.F90 modules.o reallocate.o
changemaxnptsh.o: ../changemaxnptsh.F90 modules.o reallocate.o
changemaxnptstot.o: ../changemaxnptstot.F90 modules.o reallocate.o
changemaxnpwt.o: ../changemaxnpwt.F90 modules.o reallocate.o
changemaxnqstepfit.o: ../changemaxnqstepfit.F90 reallocate.o modules.o
changemaxnset.o: ../changemaxnset.F90 modules.o reallocate.o
changemaxnspcellattot.o: ../changemaxnspcellattot.F90 reallocate.o modules.o
changemaxnspcellattotbo.o: ../changemaxnspcellattotbo.F90 reallocate.o modules.o
changemaxnwstepfit.o: ../changemaxnwstepfit.F90 reallocate.o modules.o
changemaxobs.o: ../changemaxobs.F90 modules.o reallocate.o
changemaxobsmode.o: ../changemaxobsmode.F90 modules.o reallocate.o
changemaxone.o: ../changemaxone.F90 reallocate.o modules.o
changemaxplanepot.o: ../changemaxplanepot.F90 modules.o reallocate.o
changemaxpot.o: ../changemaxpot.F90 modules.o reallocate.o
changemaxpotsites.o: ../changemaxpotsites.F90 modules.o reallocate.o
changemaxppt.o: ../changemaxppt.F90 modules.o reallocate.o
changemaxproj.o: ../changemaxproj.F90 modules.o reallocate.o
changemaxproji.o: ../changemaxproji.F90 modules.o reallocate.o
changemaxpts.o: ../changemaxpts.F90 reallocate.o modules.o
changemaxqrange.o: ../changemaxqrange.F90 reallocate.o modules.o
changemaxqatoms.o: ../changemaxqatoms.F90 modules.o reallocate.o
changemaxqatoms2.o: ../changemaxqatoms2.F90 modules.o reallocate.o
changemaxqvector.o: ../changemaxqvector.F90 reallocate.o modules.o
changemaxr1at.o: ../changemaxr1at.F90 modules.o reallocate.o
changemaxr2at.o: ../changemaxr2at.F90 modules.o reallocate.o
changemaxreaxfffixqspec.o: ../changemaxreaxfffixqspec.F90 modules.o reallocate.o
changemaxreaxffspec.o: ../changemaxreaxffspec.F90 modules.o reallocate.o
changemaxreaxffval3.o: ../changemaxreaxffval3.F90 modules.o reallocate.o
changemaxregion.o: ../changemaxregion.F90 modules.o reallocate.o
changemaxsa.o: ../changemaxsa.F90 modules.o reallocate.o
changemaxsasparticles.o: ../changemaxsasparticles.F90 modules.o reallocate.o
changemaxsasparticlespart.o: ../changemaxsasparticlespart.F90 modules.o reallocate.o
changemaxsix.o: ../changemaxsix.F90 reallocate.o modules.o
changemaxskpt.o: ../changemaxskpt.F90 modules.o reallocate.o
changemaxspcell.o: ../changemaxspcell.F90 modules.o reallocate.o
changemaxspcellbo.o: ../changemaxspcellbo.F90 modules.o reallocate.o
changemaxspec.o: ../changemaxspec.F90 modules.o m_pr.o reallocate.o
changemaxsymop.o: ../changemaxsymop.F90 modules.o reallocate.o
changemaxtempramp.o: ../changemaxtempramp.F90 modules.o reallocate.o
changemaxtdfield.o: ../changemaxtdfield.F90 modules.o reallocate.o
changemaxthb.o: ../changemaxthb.F90 modules.o reallocate.o
changemaxtitle.o: ../changemaxtitle.F90 modules.o reallocate.o iomod.o
changemaxtotint.o: ../changemaxtotint.F90 modules.o reallocate.o
changemaxtotr1at.o: ../changemaxtotr1at.F90 modules.o reallocate.o
changemaxtotvac.o: ../changemaxtotvac.F90 modules.o reallocate.o
changemaxtrialatom.o: ../changemaxtrialatom.F90 modules.o reallocate.o
changemaxuffspec.o: ../changemaxuffspec.F90 modules.o reallocate.o
changemaxvac.o: ../changemaxvac.F90 modules.o reallocate.o
changemaxvar.o: ../changemaxvar.F90 modules.o reallocate.o
changemaxvectorpair.o: ../changemaxvectorpair.F90 modules.o reallocate.o
changemaxword.o: ../changemaxword.F90 modules.o reallocate.o iomod.o
channels.o: ../channels.F90 iomod.o modules.o
checkdefshspec.o: ../checkdefshspec.F90 modules.o m_trace.o
checkquaternions.o: ../checkquaternions.F90 modules.o m_trace.o
checkpot.o: ../checkpot.F90 modules.o iomod.o m_trace.o
chemshellenergy.o: ../chemshellenergy.F90 modules.o iomod.o m_trace.o
chvo.o: ../chvo.F90 datatypes.o m_trace.o
closestimage.o: ../closestimage.F90 datatypes.o modules.o m_trace.o
collect.o: ../collect.F90 datatypes.o m_trace.o
compressd1.o: ../compressd1.F90 datatypes.o m_trace.o
compressd2.o: ../compressd2.F90 datatypes.o m_trace.o
compressd3.o: ../compressd3.F90 datatypes.o m_trace.o
connectunwrap.o: ../connectunwrap.F90 modules.o m_trace.o
connectwrap.o: ../connectwrap.F90 modules.o m_trace.o
connectwraptrial.o: ../connectwraptrial.F90 modules.o m_trace.o
copycfg.o: ../copycfg.F90 modules.o m_pr.o reallocate.o m_trace.o iomod.o
cosmoaiiderv.o: ../cosmoaiiderv.F90 datatypes.o modules.o m_trace.o
cosmoaiizero.o: ../cosmoaiizero.F90 datatypes.o modules.o m_trace.o
cosmoaijderv.o: ../cosmoaijderv.F90 datatypes.o modules.o m_trace.o
cosmoaijzero.o: ../cosmoaijzero.F90 datatypes.o modules.o m_trace.o
cosmoamatdadd.o: ../cosmoamatdadd.F90 modules.o m_trace.o
cosmoamatdadd1p.o: ../cosmoamatdadd1p.F90 modules.o m_trace.o
cosmoamatdaddd.o: ../cosmoamatdaddd.F90 modules.o m_trace.o
cosmoderv.o: ../cosmoderv.F90 modules.o iomod.o reallocate.o m_trace.o
cosmoderv2D.o: ../cosmoderv2D.F90 modules.o iomod.o reallocate.o m_trace.o
costfn.o: ../costfn.F90 modules.o datatypes.o iomod.o m_trace.o
cputime.o: ../cputime.F90 datatypes.o
create.o: ../create.F90 modules.o m_trace.o
createsphere.o: ../createsphere.F90 datatypes.o modules.o m_trace.o
cscheck.o: ../cscheck.F90 modules.o iomod.o m_trace.o
csignal.o: ../csignal.c
ctfunct.o: ../ctfunct.F90 iomod.o modules.o m_trace.o
cutscheck.o: ../cutscheck.F90 modules.o iomod.o m_trace.o
d1charge.o: ../d1charge.F90 modules.o m_trace.o
d1charge3.o: ../d1charge3.F90 modules.o m_trace.o
d1charged.o: ../d1charged.F90 modules.o m_trace.o
d1charges.o: ../d1charges.F90 modules.o m_trace.o
d2charge.o: ../d2charge.F90 modules.o m_trace.o m_strain.o
d2charge2D.o: ../d2charge2D.F90 modules.o m_trace.o m_strain.o
d2charge3.o: ../d2charge3.F90 modules.o m_trace.o
d2charge3fc.o: ../d2charge3fc.F90 modules.o m_trace.o
d2charge3p.o: ../d2charge3p.F90 modules.o m_trace.o
d2chargefc.o: ../d2chargefc.F90 modules.o m_trace.o
d2chargep.o: ../d2chargep.F90 modules.o m_trace.o
d2chargep2D.o: ../d2chargep2D.F90 modules.o m_trace.o
datatypes.o: ../datatypes.F90
datetime.o: ../datetime.F90 iomod.o modules.o m_trace.o
dbcgsolve.o: ../dbcgsolve.F90 datatypes.o iomod.o modules.o m_trace.o
dcharge.o: ../dcharge.F90 modules.o iomod.o m_trace.o m_strain.o
dcharged.o: ../dcharged.F90 modules.o iomod.o m_trace.o m_strain.o
dcharges.o: ../dcharges.F90 modules.o iomod.o m_trace.o m_strain.o
dchargesplit.o: ../dchargesplit.F90 modules.o iomod.o m_trace.o m_strain.o
dchargesplitd.o: ../dchargesplitd.F90 modules.o iomod.o m_trace.o m_strain.o
dcosmowt.o: ../dcosmowt.F90 modules.o m_trace.o
defect.o: ../defect.F90 modules.o iomod.o m_trace.o
defener.o: ../defener.F90 modules.o m_trace.o
defequ.o: ../defequ.F90 modules.o m_trace.o
defequ2a.o: ../defequ2a.F90 modules.o m_trace.o
deffreq.o: ../deffreq.F90 modules.o iomod.o m_trace.o
deffreqd.o: ../deffreqd.F90 modules.o iomod.o m_trace.o
deffun.o: ../deffun.F90 modules.o m_trace.o
defoin.o: ../defoin.F90 modules.o iomod.o m_trace.o
defopt.o: ../defopt.F90 modules.o iomod.o reallocate.o m_trace.o
defout.o: ../defout.F90 modules.o iomod.o m_trace.o
defpot.o: ../defpot.F90 modules.o iomod.o m_trace.o
defsec.o: ../defsec.F90 modules.o iomod.o m_trace.o
defsym.o: ../defsym.F90 modules.o iomod.o m_trace.o
deftmat.o: ../deftmat.F90 modules.o iomod.o m_trace.o
defvartocfg.o: ../defvartocfg.F90 modules.o m_trace.o
defword.o: ../defword.F90 datatypes.o modules.o iomod.o reallocate.o
density0.o: ../density0.F90 modules.o m_trace.o
density112.o: ../density112.F90 modules.o m_trace.o
density11s2.o: ../density11s2.F90 modules.o m_trace.o
density12a2.o: ../density12a2.F90 modules.o m_trace.o
density12as2.o: ../density12as2.F90 modules.o m_trace.o
density1fc.o: ../density1fc.F90 modules.o datatypes.o m_trace.o
density2a.o: ../density2a.F90 modules.o m_trace.o
density3.o: ../density3.F90 modules.o datatypes.o m_trace.o
density3s.o: ../density3s.F90 modules.o datatypes.o m_trace.o
densitysd.o: ../densitysd.F90 modules.o datatypes.o m_trace.o
derfc.o: ../derfc.F90 datatypes.o
destroy.o: ../destroy.F90 modules.o m_trace.o
diagdyn.o: ../diagdyn.F90 modules.o m_trace.o
diagonal_fc.o: ../diagonal_fc.F90 modules.o m_trace.o
diagonal_fcd.o: ../diagonal_fcd.F90 modules.o m_trace.o
dipole3D.o: ../dipole3D.F90 modules.o m_trace.o m_strain.o
dipole3Dd.o: ../dipole3Dd.F90 modules.o m_trace.o m_strain.o
dipole3Dtrial.o: ../dipole3Dtrial.F90 modules.o m_trace.o
dipolesd.o: ../dipolesd.F90 modules.o m_trace.o m_strain.o
distance.o: ../distance.F90 modules.o iomod.o reallocate.o m_trace.o
dlamch.o: ../dlamch.F90
dpoccon.o: ../dpoccon.F90 modules.o m_trace.o
dspecial.o: ../dspecial.F90 modules.o iomod.o m_trace.o
dtotcosmowt.o: ../dtotcosmowt.F90 modules.o m_trace.o
dump34.o: ../dump34.F90 modules.o m_trace.o
dumpbo.o: ../dumpbo.F90 modules.o m_trace.o
dumpdur.o: ../dumpdur.F90 modules.o m_pr.o m_trace.o iomod.o
dumpgen.o: ../dumpgen.F90 modules.o gulpcml.o m_plumed.o m_trace.o m_alamode.o iomod.o
dumppot.o: ../dumppot.F90 modules.o modules_kim.o m_trace.o 
dutils.o: ../dutils.F90 modules.o m_trace.o m_strain.o
dynamic.o: ../dynamic.F90 modules.o m_trace.o
dynamic_fc.o: ../dynamic_fc.F90 modules.o m_trace.o
dynamicn.o: ../dynamicn.F90 modules.o iomod.o m_trace.o
dynamicn_defect.o: ../dynamicn_defect.F90 modules.o m_trace.o
dynamicn_fc.o: ../dynamicn_fc.F90 modules.o iomod.o m_bondorder.o m_edip.o m_trace.o
dynamicn_fcd.o: ../dynamicn_fcd.F90 modules.o iomod.o m_bondorder.o m_edip.o m_trace.o
eamfnderv.o: ../eamfnderv.F90 datatypes.o modules.o m_trace.o
eamrho.o: ../eamrho.F90 modules.o m_trace.o m_strain.o
eamscalescrho.o: ../eamscalescrho.F90 modules.o m_trace.o
eckart.o: ../eckart.F90 modules.o iomod.o m_trace.o
edip1fc.o: ../edip1fc.F90 datatypes.o modules.o iomod.o m_trace.o
edip2body.o: ../edip2body.F90 modules.o m_trace.o
edip3body.o: ../edip3body.F90 modules.o m_trace.o datatypes.o
edipcn.o: ../edipcn.F90 modules.o m_trace.o datatypes.o
edipmd.o: ../edipmd.F90 datatypes.o modules.o iomod.o m_trace.o m_strain.o
eem.o: ../eem.F90 modules.o iomod.o m_trace.o
eemd.o: ../eemd.F90 modules.o iomod.o m_trace.o
eemdm.o: ../eemdm.F90 modules.o iomod.o m_trace.o
eemsplit.o: ../eemsplit.F90 modules.o iomod.o m_trace.o
eemsplitd.o: ../eemsplitd.F90 modules.o iomod.o m_trace.o
einstein.o: ../einstein.F90 modules.o m_trace.o
einsteinpd.o: ../einsteinpd.F90 modules.o m_trace.o
einsteinsd2.o: ../einsteinsd2.F90 modules.o m_trace.o
einsteintrial.o: ../einsteintrial.F90 modules.o m_trace.o
eispack.o: ../eispack.F90 datatypes.o
electricfield.o: ../electricfield.F90 modules.o m_trace.o
electricfieldparts.o: ../electricfieldparts.F90 modules.o m_trace.o
electricfieldpotl.o: ../electricfieldpotl.F90 modules.o m_trace.o
electricfieldtrial.o: ../electricfieldtrial.F90 modules.o m_trace.o
emfunc.o: ../emfunc.F90 datatypes.o m_trace.o
endstring.o: ../endstring.F90 datatypes.o
energy.o: ../energy.F90 modules.o iomod.o modules_kim.o m_trace.o
energyfd.o: ../energyfd.F90 modules.o iomod.o modules_kim.o m_bondorder.o m_brenner.o m_edip.o m_reaxff.o m_trace.o
epot.o: ../epot.F90 modules.o datatypes.o m_trace.o
epsilon.o: ../epsilon.F90 datatypes.o
equpos.o: ../equpos.F90 modules.o m_trace.o iomod.o
exchng.o: ../exchng.F90 datatypes.o m_trace.o
fefunct.o: ../fefunct.F90 modules.o m_trace.o
fefunctn.o: ../fefunctn.F90 modules.o m_trace.o
fenergy0.o: ../fenergy0.F90 modules.o iomod.o m_trace.o
fenergy0d.o: ../fenergy0d.F90 modules.o iomod.o m_trace.o
fenergy3.o: ../fenergy3.F90 modules.o iomod.o m_trace.o
fenergy3d.o: ../fenergy3d.F90 modules.o iomod.o m_trace.o
findcelltype.o: ../findcelltype.F90 datatypes.o m_trace.o
finddef.o: ../finddef.F90 modules.o m_trace.o
findfractionword.o: ../findfractionword.F90 datatypes.o m_trace.o
findspacegroup.o: ../findspacegroup.F90 modules.o m_trace.o
finish.o: ../finish.F90 modules.o iomod.o reallocate.o m_trace.o
firstpass.o: ../firstpass.F90 modules.o iomod.o m_trace.o
fit.o: ../fit.F90 modules.o iomod.o m_simplex.o m_trace.o
fitbfgs.o: ../fitbfgs.F90 modules.o iomod.o m_trace.o
fitcon.o: ../fitcon.F90 modules.o m_trace.o
fitfun.o: ../fitfun.F90 modules.o m_trace.o
fitfunct.o: ../fitfunct.F90 modules.o iomod.o m_trace.o
fitgrad.o: ../fitgrad.F90 modules.o m_trace.o
fitmin.o: ../fitmin.F90 modules.o iomod.o m_trace.o
fitword.o: ../fitword.F90 modules.o iomod.o
force.o: ../force.F90 modules.o m_trace.o m_strain.o
forcebycell.o: ../forcebycell.F90 modules.o m_bondorder.o m_edip.o m_trace.o
forcetrial.o: ../forcetrial.F90 modules.o m_trace.o
formula.o: ../formula.F90 datatypes.o m_trace.o modules.o
four0d3.o: ../four0d3.F90 modules.o m_trace.o
four12.o: ../four12.F90 modules.o m_trace.o
four12d.o: ../four12d.F90 modules.o m_trace.o
four1fc.o: ../four1fc.F90 modules.o m_trace.o
four3d3.o: ../four3d3.F90 modules.o m_trace.o
four3d3dV.o: ../four3d3dV.F90 modules.o m_trace.o
four3d3dVg.o: ../four3d3dVg.F90 modules.o m_trace.o
four3d3g.o: ../four3d3g.F90 modules.o m_trace.o
fourbody.o: ../fourbody.F90 modules.o iomod.o m_trace.o datatypes.o
fourfc.o: ../fourfc.F90 modules.o m_trace.o
fourfc3.o: ../fourfc3.F90 modules.o m_trace.o
fourfcd.o: ../fourfcd.F90 modules.o m_trace.o
fourfd.o: ../fourfd.F90 modules.o datatypes.o m_trace.o m_strain.o
fourimp.o: ../fourimp.F90 modules.o iomod.o m_trace.o m_strain.o
fourimp12.o: ../fourimp12.F90 modules.o m_trace.o
fourimp12d.o: ../fourimp12d.F90 modules.o m_trace.o
fourimp1fc.o: ../fourimp1fc.F90 modules.o m_trace.o
fourimpd.o: ../fourimpd.F90 modules.o iomod.o m_trace.o m_strain.o
fourimpfc.o: ../fourimpfc.F90 modules.o m_trace.o
fourimpfc3.o: ../fourimpfc3.F90 modules.o m_trace.o
fourimpfcd.o: ../fourimpfcd.F90 modules.o m_trace.o
fourimpfd.o: ../fourimpfd.F90 modules.o datatypes.o m_trace.o m_strain.o
fourimpmc.o: ../fourimpmc.F90 modules.o iomod.o m_trace.o
fourimpmd.o: ../fourimpmd.F90 modules.o datatypes.o iomod.o m_trace.o m_strain.o
fourimpmds.o: ../fourimpmds.F90 modules.o iomod.o m_trace.o m_strain.o
fourimpp.o: ../fourimpp.F90 modules.o m_trace.o
fourimppd.o: ../fourimppd.F90 modules.o m_trace.o
fourimpsd.o: ../fourimpsd.F90 modules.o datatypes.o iomod.o m_trace.o m_strain.o
fourimpsd2.o: ../fourimpsd2.F90 modules.o iomod.o m_trace.o m_strain.o
fourlist.o: ../fourlist.F90 modules.o iomod.o m_trace.o m_strain.o
fourmc.o: ../fourmc.F90 modules.o iomod.o m_trace.o
fourmd.o: ../fourmd.F90 modules.o datatypes.o iomod.o m_trace.o m_strain.o
fourmds.o: ../fourmds.F90 modules.o iomod.o m_trace.o m_strain.o
fournos.o: ../fournos.F90 modules.o iomod.o m_trace.o m_strain.o
fournosd.o: ../fournosd.F90 modules.o iomod.o m_trace.o m_strain.o
fouroop.o: ../fouroop.F90 modules.o iomod.o m_trace.o m_strain.o
fouroop12.o: ../fouroop12.F90 modules.o m_trace.o
fouroop12d.o: ../fouroop12d.F90 modules.o m_trace.o
fouroop1fc.o: ../fouroop1fc.F90 modules.o m_trace.o
fouroopd.o: ../fouroopd.F90 modules.o iomod.o m_trace.o m_strain.o
fouroopfc.o: ../fouroopfc.F90 modules.o m_trace.o
fouroopfc3.o: ../fouroopfc3.F90 modules.o m_trace.o
fouroopfcd.o: ../fouroopfcd.F90 modules.o m_trace.o
fouroopfd.o: ../fouroopfd.F90 modules.o m_trace.o m_strain.o
fouroopmc.o: ../fouroopmc.F90 modules.o iomod.o m_trace.o
fouroopmd.o: ../fouroopmd.F90 modules.o iomod.o m_trace.o m_strain.o
fouroopmds.o: ../fouroopmds.F90 modules.o iomod.o m_trace.o m_strain.o
fouroopp.o: ../fouroopp.F90 modules.o m_trace.o
fourooppd.o: ../fourooppd.F90 modules.o m_trace.o
fouroopsd.o: ../fouroopsd.F90 modules.o datatypes.o iomod.o m_trace.o m_strain.o
fouroopsd2.o: ../fouroopsd2.F90 modules.o iomod.o m_trace.o m_strain.o
fourp.o: ../fourp.F90 modules.o m_trace.o
fourpd.o: ../fourpd.F90 modules.o m_trace.o
foursd.o: ../foursd.F90 modules.o datatypes.o iomod.o m_trace.o m_strain.o
foursd2.o: ../foursd2.F90 modules.o iomod.o m_trace.o m_strain.o
ftow.o: ../ftow.F90 datatypes.o modules.o
funct.o: ../funct.F90 modules.o m_trace.o
functn.o: ../functn.F90 modules.o iomod.o m_trace.o
functneb.o: ../functneb.F90 modules.o m_trace.o
functnf.o: ../functnf.F90 modules.o m_trace.o
functnff.o: ../functnff.F90 modules.o m_brenner.o m_bondorder.o m_edip.o m_reaxff.o m_trace.o
functsync.o: ../functsync.F90 modules.o m_trace.o
gabcfg.o: ../gabcfg.F90 modules.o iomod.o m_trace.o
gacreate.o: ../gacreate.F90 modules.o iomod.o m_trace.o
gacross.o: ../gacross.F90 datatypes.o modules.o m_trace.o
gaexpd.o: ../gaexpd.F90 modules.o iomod.o m_trace.o
gafit.o: ../gafit.F90 modules.o iomod.o m_trace.o datatypes.o
gamaxmin.o: ../gamaxmin.F90 modules.o m_trace.o
gamma.o: ../gamma.F90 datatypes.o m_trace.o modules.o
gammasm.o: ../gammasm.F90 datatypes.o m_trace.o
gamuta.o: ../gamuta.F90 datatypes.o modules.o m_trace.o
gaopt.o: ../gaopt.F90 modules.o iomod.o m_trace.o
garepro.o: ../garepro.F90 modules.o iomod.o m_trace.o
gasort.o: ../gasort.F90 datatypes.o iomod.o modules.o m_trace.o
gasteiger.o: ../gasteiger.F90 modules.o iomod.o m_trace.o
gasubs.o: ../gasubs.F90 modules.o iomod.o datatypes.o
gauss.o: ../gauss.F90 datatypes.o modules.o iomod.o m_trace.o
gaword.o: ../gaword.F90 modules.o iomod.o
genpot.o: ../genpot.F90 modules.o m_trace.o
genpotd.o: ../genpotd.F90 modules.o m_trace.o
genpotdm.o: ../genpotdm.F90 modules.o m_trace.o
genpotsplitd.o: ../genpotsplitd.F90 modules.o m_trace.o
genword.o: ../genword.F90 modules.o gulpcml.o iomod.o m_pdfneutron.o m_plumed.o
getangles.o: ../getangles.F90 modules.o iomod.o m_trace.o
getbocharge.o: ../getbocharge.F90 datatypes.o modules.o iomod.o m_trace.o m_strain.o
getbocharged.o: ../getbocharged.F90 datatypes.o modules.o iomod.o m_trace.o m_strain.o
getbochargetrial.o: ../getbochargetrial.F90 datatypes.o modules.o iomod.o m_trace.o
getbonds.o: ../getbonds.F90 modules.o m_trace.o
getcoordno.o: ../getcoordno.F90 modules.o iomod.o m_trace.o
getderv1.o: ../getderv1.F90 modules.o m_trace.o
getdipole.o: ../getdipole.F90 modules.o m_trace.o
getdistances.o: ../getdistances.F90 modules.o m_trace.o
getgcindex.o: ../getgcindex.F90 modules.o m_trace.o
getkeywrd.o: ../getkeywrd.F90 modules.o iomod.o m_trace.o
getline.o: ../getline.F90 modules.o iomod.o
getneighbour.o: ../getneighbour.F90 datatypes.o modules.o m_trace.o
getnmr.o: ../getnmr.F90 datatypes.o modules.o reallocate.o m_trace.o
getphi.o: ../getphi.F90 datatypes.o m_trace.o
getpotsymbol.o: ../getpotsymbol.F90 datatypes.o modules.o m_trace.o iomod.o
getshellmatrices.o: ../getshellmatrices.F90 modules.o m_trace.o
getstrain.o: ../getstrain.F90 modules.o m_trace.o
gettorsions.o: ../gettorsions.F90 modules.o iomod.o m_trace.o
gettraceinertia.o: ../gettraceinertia.F90 modules.o iomod.o m_trace.o
getvibmode.o: ../getvibmode.F90 modules.o m_trace.o
getvibmoded.o: ../getvibmoded.F90 modules.o m_trace.o
groupvelocities.o: ../groupvelocities.F90 modules.o iomod.o m_trace.o
groupvelocitiesd.o: ../groupvelocitiesd.F90 modules.o iomod.o m_trace.o
grueneisen.o: ../grueneisen.F90 modules.o iomod.o m_trace.o
grueneisend.o: ../grueneisend.F90 modules.o iomod.o m_trace.o
grueneiseng.o: ../grueneiseng.F90 modules.o iomod.o m_trace.o
grueneisengd.o: ../grueneisengd.F90 modules.o iomod.o m_trace.o
gtheta.o: ../gtheta.F90 modules.o iomod.o m_trace.o
gthetabo.o: ../gthetabo.F90 modules.o iomod.o m_trace.o
gulp.o: ../gulp.F90 datatypes.o
gulpcml.o: ../gulpcml.F90 datatypes.o modules.o m_pdfneutron.o iomod.o m_pdfvariables.o
gulpcml_md.o: ../gulpcml_md.F90 gulpcml.o datatypes.o modules.o
gulpcml_phonon.o: ../gulpcml_phonon.F90 gulpcml.o datatypes.o modules.o iomod.o m_pdfneutron.o m_pdfvariables.o
gulpcml_props.o: ../gulpcml_props.F90 datatypes.o gulpcml.o modules.o
gulpfinish.o: ../gulpfinish.F90 modules.o gulpcml.o
gulpmain.o: ../gulpmain.F90 modules.o iomod.o license.o m_trace.o
gulpsetup.o: ../gulpsetup.F90 modules.o gulpcml.o gulpcml_phonon.o iomod.o m_pdfneutron.o m_trace.o license.o
harmonicrelax.o: ../harmonicrelax.F90 modules.o iomod.o m_trace.o
help.o: ../help.F90 modules.o iomod.o
hfunc.o: ../hfunc.F90 datatypes.o m_trace.o
improper.o: ../improper.F90 modules.o iomod.o m_trace.o
ind1toijk.o: ../ind1toijk.F90 datatypes.o m_trace.o
indsft.o: ../indsft.F90 datatypes.o m_trace.o
indtoijk.o: ../indtoijk.F90 datatypes.o m_trace.o
init1bodydefaults.o: ../init1bodydefaults.F90 modules.o m_trace.o
init2bodydefaults.o: ../init2bodydefaults.F90 m_trace.o modules.o
init_ins.o: ../init_ins.F90 modules.o m_trace.o
initarc.o: ../initarc.F90 datatypes.o modules.o m_trace.o iomod.o
initcomms.o: ../initcomms.F90 iomod.o modules.o
initdervs.o: ../initdervs.F90 modules.o iomod.o m_trace.o
initial.o: ../initial.F90 modules.o modules_kim.o m_pdfneutron.o m_plumed.o m_pr.o m_trace.o
initktrm.o: ../initktrm.F90 modules.o m_trace.o
initmaxatdefaults.o: ../initmaxatdefaults.F90 modules.o
initmaxatotdefaults.o: ../initmaxatotdefaults.F90 modules.o
initmaxbonddefaults.o: ../initmaxbonddefaults.F90 modules.o
initmaxbondqdefaults.o: ../initmaxbondqdefaults.F90 modules.o
initmaxccspecdefaults.o: ../initmaxccspecdefaults.F90 modules.o
initmaxcfgdefaults.o: ../initmaxcfgdefaults.F90 modules.o modules_kim.o m_pr.o reallocate.o
initmaxconnectdefaults.o: ../initmaxconnectdefaults.F90 modules.o
initmaxdefdefaults.o: ../initmaxdefdefaults.F90 modules.o
initmaxeamdendefaults.o: ../initmaxeamdendefaults.F90 modules.o
initmaxeamfnspecdefaults.o: ../initmaxeamfnspecdefaults.F90 modules.o
initmaxeamspecdefaults.o: ../initmaxeamspecdefaults.F90 modules.o
initmaxedipspecdefaults.o: ../initmaxedipspecdefaults.F90 modules.o
initmaxfgraddefaults.o: ../initmaxfgraddefaults.F90 modules.o
initmaxfitdefaults.o: ../initmaxfitdefaults.F90 modules.o
initmaxfordefaults.o: ../initmaxfordefaults.F90 modules.o
initmaxfstraindefaults.o: ../initmaxfstraindefaults.F90 modules.o
initmaxgcmcmoldefaults.o: ../initmaxgcmcmoldefaults.F90 modules.o
initmaxlibdefaults.o: ../initmaxlibdefaults.F90 modules.o
initmaxqrangedefaults.o: ../initmaxqrangedefaults.F90 modules.o
initmaxmcswaps.o: ../initmaxmcswaps.F90 modules.o
initmaxmcswapspec.o: ../initmaxmcswapspec.F90 modules.o
initmaxmoldefaults.o: ../initmaxmoldefaults.F90 modules.o
initmaxnboadefaults.o: ../initmaxnboadefaults.F90 modules.o
initmaxnboodefaults.o: ../initmaxnboodefaults.F90 modules.o
initmaxnbopotdefaults.o: ../initmaxnbopotdefaults.F90 modules.o
initmaxnboq0defaults.o: ../initmaxnboq0defaults.F90 modules.o
initmaxnboqdefaults.o: ../initmaxnboqdefaults.F90 modules.o
initmaxnbordefaults.o: ../initmaxnbordefaults.F90 modules.o
initmaxnbozdefaults.o: ../initmaxnbozdefaults.F90 modules.o
initmaxnebreplicatotdefaults.o: ../initmaxnebreplicatotdefaults.F90 modules.o
initmaxnppadefaults.o: ../initmaxnppadefaults.F90 datatypes.o modules.o
initmaxnptsdefaults.o: ../initmaxnptsdefaults.F90 modules.o
initmaxobsdefaults.o: ../initmaxobsdefaults.F90 modules.o
initmaxplanepotdefaults.o: ../initmaxplanepotdefaults.F90 modules.o
initmaxr1atdefaults.o: ../initmaxr1atdefaults.F90 modules.o
initmaxreaxffspecdefaults.o: ../initmaxreaxffspecdefaults.F90 modules.o
initmaxreaxffval3defaults.o: ../initmaxreaxffval3defaults.F90 modules.o
initmaxregiondefaults.o: ../initmaxregiondefaults.F90 modules.o
initmaxsixdefaults.o: ../initmaxsixdefaults.F90 modules.o
initmaxspcellbodefaults.o: ../initmaxspcellbodefaults.F90 modules.o
initmaxspcelldefaults.o: ../initmaxspcelldefaults.F90 modules.o
initmaxspecdefaults.o: ../initmaxspecdefaults.F90 modules.o
initmaxthbdefaults.o: ../initmaxthbdefaults.F90 modules.o
initmaxtitledefaults.o: ../initmaxtitledefaults.F90 modules.o
initmemory.o: ../initmemory.F90 modules.o m_trace.o
initqmatrix.o: ../initqmatrix.F90 modules.o m_trace.o
initqmatrixc.o: ../initqmatrixc.F90 datatypes.o modules.o m_trace.o
insertatom.o: ../insertatom.F90 modules.o m_trace.o
inword.o: ../inword.F90 modules.o iomod.o modules_kim.o m_pdfneutron.o m_trace.o
iomod.o: ../iomod.F90 datatypes.o
istoan.o: ../istoan.F90 modules.o iomod.o
itcom.o: ../itcom.F90 datatypes.o
iterinv.o: ../iterinv.F90 datatypes.o m_trace.o
itow.o: ../itow.F90 datatypes.o
itpack.o: ../itpack.F90 datatypes.o itcom.o m_trace.o
kim_functions.o: ../kim_functions.F90 datatypes.o m_trace.o modules.o iomod.o modules_kim.o reallocate.o
kimmd.o: ../kimmd.F90 datatypes.o modules.o iomod.o modules_kim.o kim_functions.o m_trace.o
kindex.o: ../kindex.F90 modules.o iomod.o m_trace.o
kinthewedge.o: ../kinthewedge.F90 modules.o m_trace.o
kvector.o: ../kvector.F90 modules.o iomod.o m_trace.o
kvectorf.o: ../kvectorf.F90 modules.o iomod.o m_trace.o
kweight.o: ../kweight.F90 datatypes.o m_trace.o
label.o: ../label.F90 modules.o
labellc.o: ../labellc.F90 modules.o
lambda.o: ../lambda.F90 modules.o iomod.o m_trace.o
langevinpol.o: ../langevinpol.F90 datatypes.o
lapack.o: ../lapack.F90 datatypes.o
latpot.o: ../latpot.F90 modules.o gulpcml.o iomod.o m_trace.o
license.o: ../license.F90 modules.o iomod.o
linepro.o: ../linepro.F90 modules.o iomod.o m_trace.o
linpack.o: ../linpack.F90 datatypes.o
lintoijk.o: ../lintoijk.F90 datatypes.o m_trace.o
lmatch.o: ../lmatch.F90 datatypes.o modules.o
lmatchany.o: ../lmatchany.F90 datatypes.o modules.o
lmatchpair.o: ../lmatchpair.F90 datatypes.o modules.o
lmbfgs.o: ../lmbfgs.F90 datatypes.o iomod.o m_trace.o modules.o
local.o: ../local.F90 modules.o
lower.o: ../lower.F90 modules.o iomod.o m_trace.o
ltont.o: ../ltont.F90 datatypes.o
m_alamode.o: ../m_alamode.F90 datatypes.o modules.o iomod.o
m_bondorder.o: ../m_bondorder.F90 datatypes.o modules.o iomod.o
m_brenner.o: ../m_brenner.F90 datatypes.o modules.o iomod.o
m_conjgr.o: ../m_conjgr.F90 datatypes.o
m_edip.o: ../m_edip.F90 datatypes.o modules.o iomod.o
m_fft3D.o: ../m_fft3D.F90 datatypes.o modules.o
m_pdf.o: ../m_pdf.F90 datatypes.o m_pdfneutron.o modules.o m_pdfvariables.o iomod.o
m_pdfneutron.o: ../m_pdfneutron.F90 datatypes.o iomod.o modules.o reallocate.o iomod.o
m_pdfvariables.o: ../m_pdfvariables.F90 datatypes.o
m_plumed.o: ../m_plumed.F90 datatypes.o modules.o iomod.o m_pr.o
m_pr.o: ../m_pr.F90 datatypes.o
m_reaxff.o: ../m_reaxff.F90 datatypes.o modules.o iomod.o
m_simplex.o: ../m_simplex.F90 datatypes.o modules.o iomod.o
m_strain.o: ../m_strain.F90 datatypes.o modules.o m_trace.o
m_trace.o: ../m_trace.F90 datatypes.o modules.o iomod.o
many.o: ../many.F90 modules.o iomod.o m_trace.o m_strain.o
many0d.o: ../many0d.F90 modules.o iomod.o m_trace.o m_strain.o
many0dd.o: ../many0dd.F90 modules.o iomod.o m_trace.o m_strain.o
many12.o: ../many12.F90 modules.o iomod.o m_trace.o
many12d.o: ../many12d.F90 modules.o iomod.o m_trace.o
many1fc.o: ../many1fc.F90 modules.o m_trace.o
manyd.o: ../manyd.F90 modules.o iomod.o m_trace.o m_strain.o
manyfc.o: ../manyfc.F90 modules.o m_trace.o m_strain.o
manyfcd.o: ../manyfcd.F90 modules.o m_trace.o m_strain.o
manymd0.o: ../manymd0.F90 modules.o iomod.o m_trace.o
manymd3.o: ../manymd3.F90 modules.o iomod.o m_trace.o m_strain.o
manymd3s.o: ../manymd3s.F90 modules.o iomod.o m_trace.o m_strain.o
manymi3.o: ../manymi3.F90 modules.o iomod.o m_trace.o m_strain.o
manyp.o: ../manyp.F90 modules.o m_trace.o m_strain.o
manypd.o: ../manypd.F90 modules.o m_trace.o m_strain.o
manysd.o: ../manysd.F90 modules.o datatypes.o iomod.o m_trace.o m_strain.o
manysd2.o: ../manysd2.F90 modules.o iomod.o m_trace.o m_strain.o
cmatinv.o: ../cmatinv.F90 modules.o m_trace.o datatypes.o iomod.o
matinv.o: ../matinv.F90 modules.o m_trace.o datatypes.o iomod.o
imaxall.o: ../imaxall.F90 modules.o
maxall.o: ../maxall.F90 modules.o
mccreate.o: ../mccreate.F90 modules.o m_trace.o
mcdestroy.o: ../mcdestroy.F90 modules.o m_trace.o
mcenergy.o: ../mcenergy.F90 modules.o iomod.o modules_kim.o m_trace.o
mcmolinsert.o: ../mcmolinsert.F90 modules.o m_trace.o
mcmove.o: ../mcmove.F90 modules.o reallocate.o m_trace.o
mcrotate.o: ../mcrotate.F90 modules.o reallocate.o m_trace.o
mcstorecfg.o: ../mcstorecfg.F90 modules.o m_trace.o
mcstrain.o: ../mcstrain.F90 modules.o m_trace.o
mcswap.o: ../mcswap.F90 modules.o reallocate.o m_trace.o
mctrialfunction.o: ../mctrialfunction.F90 datatypes.o modules.o m_trace.o
mcword.o: ../mcword.F90 modules.o iomod.o
mdcellind.o: ../mdcellind.F90 modules.o m_trace.o
mdcorrect.o: ../mdcorrect.F90 modules.o m_pr.o m_trace.o
mddistconpost.o: ../mddistconpost.F90 modules.o m_trace.o
mdfunct.o: ../mdfunct.F90 modules.o modules_kim.o m_plumed.o m_pr.o m_trace.o
mdgcorrect.o: ../mdgcorrect.F90 modules.o m_trace.o
mdke.o: ../mdke.F90 modules.o m_trace.o
mdkeatomicstress.o: ../mdkeatomicstress.F90 modules.o m_trace.o
mdkestrain.o: ../mdkestrain.F90 modules.o m_trace.o
mdlfvcorrect.o: ../mdlfvcorrect.F90 modules.o m_trace.o
mdpredict.o: ../mdpredict.F90 modules.o m_pr.o m_trace.o
mdprop.o: ../mdprop.F90 modules.o gulpcml.o gulpcml_md.o iomod.o m_pr.o m_trace.o
mdreset.o: ../mdreset.F90 modules.o m_trace.o
mdrestart.o: ../mdrestart.F90 modules.o iomod.o m_plumed.o m_pr.o pr_random.o m_trace.o
mdscale.o: ../mdscale.F90 modules.o iomod.o m_trace.o
mdshscale.o: ../mdshscale.F90 modules.o m_trace.o
mdstop.o: ../mdstop.F90 iomod.o modules.o
mdtemperature.o: ../mdtemperature.F90 modules.o m_pr.o pr_random.o m_trace.o
mdvelcor.o: ../mdvelcor.F90 modules.o m_trace.o
mdvvcorrect.o: ../mdvvcorrect.F90 modules.o m_trace.o
mdword.o: ../mdword.F90 modules.o iomod.o m_pr.o
mdwrite.o: ../mdwrite.F90 modules.o m_trace.o
meamfnderv.o: ../meamfnderv.F90 datatypes.o modules.o m_trace.o
meamrho.o: ../meamrho.F90 modules.o m_trace.o m_strain.o
meamscreen.o: ../meamscreen.F90 modules.o m_trace.o
meamsymdensity.o: ../meamsymdensity.F90 modules.o m_trace.o
meamtotalrho.o: ../meamtotalrho.F90 modules.o m_trace.o
meamtotalrhocrossderv.o: ../meamtotalrhocrossderv.F90 modules.o m_trace.o
meamtotalrhoderv.o: ../meamtotalrhoderv.F90 modules.o m_trace.o
meamtotalrhoscreenderv.o: ../meamtotalrhoscreenderv.F90 modules.o m_trace.o m_strain.o
methodok.o: ../methodok.F90 modules.o iomod.o m_fft3D.o m_pdfneutron.o m_trace.o
mindtoijk.o: ../mindtoijk.F90 datatypes.o m_trace.o
minimise.o: ../minimise.F90 modules.o gulpcml.o iomod.o license.o m_trace.o
modules.o: ../modules.F90 datatypes.o iomod.o
modules_kim.o: ../modules_kim.F90 datatypes.o
molind.o: ../molind.F90 modules.o iomod.o reallocate.o m_trace.o
molindfix.o: ../molindfix.F90 modules.o reallocate.o m_trace.o
molindfixtrial.o: ../molindfixtrial.F90 modules.o reallocate.o m_trace.o
molindtrial.o: ../molindtrial.F90 modules.o iomod.o reallocate.o m_trace.o
moltype.o: ../moltype.F90 modules.o m_trace.o
move2a1.o: ../move2a1.F90 modules.o m_trace.o
mpbarrier.o: ../mpbarrier.F90 modules.o
mpfinish.o: ../mpfinish.F90
mxmb.o: ../mxmb.F90 datatypes.o m_trace.o
nagamma.o: ../nagamma.F90 modules.o m_trace.o
nagammac.o: ../nagammac.F90 modules.o m_trace.o
nagammacrigid.o: ../nagammacrigid.F90 modules.o m_trace.o
nagammarigid.o: ../nagammarigid.F90 modules.o m_trace.o
nearest.o: ../nearest.F90 datatypes.o m_trace.o
nebforce.o: ../nebforce.F90 modules.o m_trace.o
nebword.o: ../nebword.F90 modules.o reallocate.o iomod.o
nrhess.o: ../nrhess.F90 modules.o iomod.o m_trace.o
nrhessn.o: ../nrhessn.F90 modules.o iomod.o m_trace.o
nrstep.o: ../nrstep.F90 datatypes.o modules.o m_trace.o
nudge.o: ../nudge.F90 datatypes.o m_trace.o modules.o
nullpointer.o: ../nullpointer.F90 modules.o modules_kim.o m_pr.o reallocate.o m_trace.o
numhess.o: ../numhess.F90 iomod.o modules.o m_trace.o
okspec.o: ../okspec.F90 modules.o iomod.o
olinmin.o: ../olinmin.F90 modules.o iomod.o m_trace.o
omegaproperty.o: ../omegaproperty.F90 modules.o iomod.o m_trace.o
onebody.o: ../onebody.F90 modules.o m_trace.o
optim.o: ../optim.F90 modules.o gulpcml.o iomod.o reallocate.o m_trace.o
optin.o: ../optin.F90 modules.o iomod.o m_trace.o
options.o: ../options.F90 modules.o gulpcml.o iomod.o m_trace.o
optout.o: ../optout.F90 modules.o iomod.o m_trace.o
optshe.o: ../optshe.F90 modules.o iomod.o m_trace.o
oscillatorstrength.o: ../oscillatorstrength.F90 modules.o iomod.o m_trace.o
oscillatorstrengthc.o: ../oscillatorstrengthc.F90 modules.o iomod.o m_trace.o
oscillatorstrengthg.o: ../oscillatorstrengthg.F90 modules.o iomod.o m_trace.o
outadvice.o: ../outadvice.F90 iomod.o modules.o
outarc.o: ../outarc.F90 modules.o iomod.o
outatomicstress.o: ../outatomicstress.F90 modules.o iomod.o
outbiograf.o: ../outbiograf.F90 modules.o iomod.o
outcell.o: ../outcell.F90 modules.o iomod.o
outcif.o: ../outcif.F90 modules.o iomod.o
outcosmo.o: ../outcosmo.F90 modules.o iomod.o
outcssr.o: ../outcssr.F90 modules.o iomod.o
outdcdheader.o: ../outdcdheader.F90 modules.o
outdefin.o: ../outdefin.F90 modules.o iomod.o
outdener.o: ../outdener.F90 modules.o iomod.o
outderv.o: ../outderv.F90 modules.o iomod.o
outdlv.o: ../outdlv.F90 modules.o iomod.o
outdrv.o: ../outdrv.F90 modules.o
outed.o: ../outed.F90 modules.o iomod.o
outeem.o: ../outeem.F90 modules.o iomod.o
outener.o: ../outener.F90 modules.o gulpcml.o iomod.o modules_kim.o m_plumed.o
outerror.o: ../outerror.F90 iomod.o modules.o
outfdf.o: ../outfdf.F90 modules.o
outfile.o: ../outfile.F90 modules.o iomod.o gulpcml.o
outfit.o: ../outfit.F90 modules.o iomod.o
outfitcon.o: ../outfitcon.F90 modules.o iomod.o
outfour.o: ../outfour.F90 modules.o iomod.o
outfrc.o: ../outfrc.F90 modules.o iomod.o
outfreq.o: ../outfreq.F90 modules.o iomod.o
outgafit.o: ../outgafit.F90 modules.o iomod.o
outgaopt.o: ../outgaopt.F90 modules.o iomod.o
outgapar.o: ../outgapar.F90 modules.o iomod.o
outgen.o: ../outgen.F90 modules.o iomod.o
outinertia.o: ../outinertia.F90 modules.o iomod.o
outkey.o: ../outkey.F90 modules.o iomod.o m_pdfneutron.o
outkpt.o: ../outkpt.F90 modules.o iomod.o
outkptfile.o: ../outkptfile.F90 modules.o
outlammps.o: ../outlammps.F90 modules.o iomod.o
outlammpspots.o: ../outlammpspots.F90 modules.o iomod.o
outmarv.o: ../outmarv.F90 modules.o iomod.o
outmarv2.o: ../outmarv2.F90 modules.o iomod.o
outmc.o: ../outmc.F90 modules.o
outmd.o: ../outmd.F90 modules.o
outmdin.o: ../outmdin.F90 modules.o iomod.o m_pr.o
outmol.o: ../outmol.F90 modules.o iomod.o
outneb.o: ../outneb.F90 modules.o iomod.o
outnebin.o: ../outnebin.F90 modules.o iomod.o
outobs.o: ../outobs.F90 modules.o iomod.o
outofmemory.o: ../outofmemory.F90 iomod.o modules.o
outofp.o: ../outofp.F90 modules.o iomod.o
outoptvar.o: ../outoptvar.F90 modules.o iomod.o
outphon.o: ../outphon.F90 modules.o iomod.o
outpolar.o: ../outpolar.F90 iomod.o modules.o
outpot.o: ../outpot.F90 modules.o iomod.o modules_kim.o
outqbo.o: ../outqbo.F90 modules.o
outshengBTE.o: ../outshengBTE.F90 modules.o m_pdfneutron.o iomod.o
outsiteenergy.o: ../outsiteenergy.F90 modules.o iomod.o
outspec.o: ../outspec.F90 modules.o iomod.o
outspec_pdf.o: ../outspec_pdf.F90 modules.o iomod.o m_pdfneutron.o
outstresses.o: ../outstresses.F90 modules.o iomod.o
outstructure.o: ../outstructure.F90 modules.o iomod.o
outsync.o: ../outsync.F90 modules.o iomod.o
outsyncin.o: ../outsyncin.F90 modules.o iomod.o
outthb.o: ../outthb.F90 modules.o
outthree.o: ../outthree.F90 modules.o iomod.o
outtime.o: ../outtime.F90 iomod.o modules.o m_plumed.o
outtwo.o: ../outtwo.F90 modules.o iomod.o
outvar.o: ../outvar.F90 modules.o iomod.o
outvarg.o: ../outvarg.F90 modules.o iomod.o
outwarning.o: ../outwarning.F90 iomod.o modules.o
outxr.o: ../outxr.F90 modules.o iomod.o
outxsf.o: ../outxsf.F90 modules.o iomod.o
outxtl.o: ../outxtl.F90 modules.o iomod.o
outxyz.o: ../outxyz.F90 modules.o iomod.o
pdiagd.o: ../pdiagd.F90 modules.o m_trace.o
pdiage.o: ../pdiage.F90 modules.o m_pdfneutron.o m_trace.o
pdiagg.o: ../pdiagg.F90 modules.o m_trace.o
pdiaggd.o: ../pdiaggd.F90 modules.o m_trace.o
pdiagl.o: ../pdiagl.F90 modules.o m_pdfneutron.o m_trace.o
peigen.o: ../peigen.F90 modules.o iomod.o m_pdfneutron.o m_trace.o
peigend.o: ../peigend.F90 modules.o iomod.o m_pdfneutron.o m_trace.o
peigeng.o: ../peigeng.F90 modules.o iomod.o m_pdfneutron.o m_trace.o
peigengd.o: ../peigengd.F90 modules.o iomod.o m_pdfneutron.o m_trace.o
phase_fc.o: ../phase_fc.F90 modules.o datatypes.o m_trace.o
phase_fcd.o: ../phase_fcd.F90 modules.o datatypes.o m_trace.o
phase_fcg.o: ../phase_fcg.F90 modules.o datatypes.o m_trace.o
phoncopy1.o: ../phoncopy1.F90 datatypes.o m_trace.o
phoncopy2.o: ../phoncopy2.F90 datatypes.o modules.o m_trace.o
phonon.o: ../phonon.F90 modules.o gulpcml.o gulpcml_phonon.o iomod.o m_pdf.o m_pdfneutron.o m_trace.o
phonon_fc.o: ../phonon_fc.F90 modules.o gulpcml.o gulpcml_phonon.o iomod.o m_pdf.o m_pdfneutron.o m_trace.o
phonon_fcd.o: ../phonon_fcd.F90 modules.o gulpcml.o gulpcml_phonon.o iomod.o m_pdf.o m_pdfneutron.o m_trace.o
phonon_fcg.o: ../phonon_fcg.F90 modules.o gulpcml.o gulpcml_phonon.o iomod.o m_pdfneutron.o m_trace.o
phonond.o: ../phonond.F90 modules.o gulpcml.o gulpcml_phonon.o iomod.o m_pdf.o m_pdfneutron.o m_trace.o
phonword.o: ../phonword.F90 modules.o iomod.o m_pdfneutron.o m_alamode.o
pireal0d.o: ../pireal0d.F90 modules.o m_trace.o
pireal0dd.o: ../pireal0dd.F90 modules.o m_trace.o
pirealrecip.o: ../pirealrecip.F90 modules.o datatypes.o m_trace.o m_strain.o
pirealrecipd.o: ../pirealrecipd.F90 modules.o datatypes.o m_trace.o m_strain.o
pirealrecipp.o: ../pirealrecipp.F90 modules.o datatypes.o m_trace.o m_strain.o
pirealrecippd.o: ../pirealrecippd.F90 modules.o datatypes.o m_trace.o m_strain.o
pirealrecips.o: ../pirealrecips.F90 modules.o datatypes.o m_trace.o m_strain.o
planepot2D.o: ../planepot2D.F90 modules.o m_trace.o
planepotfd.o: ../planepotfd.F90 modules.o m_trace.o
planepotmc.o: ../planepotmc.F90 modules.o m_trace.o
planepotmd.o: ../planepotmd.F90 modules.o m_trace.o
poccon.o: ../poccon.F90 modules.o m_trace.o
poccons.o: ../poccons.F90 modules.o m_trace.o
polarisation.o: ../polarisation.F90 modules.o m_trace.o
potcut.o: ../potcut.F90 modules.o iomod.o m_trace.o
potential.o: ../potential.F90 modules.o m_trace.o
potgrid.o: ../potgrid.F90 modules.o iomod.o m_trace.o
potword21.o: ../potword21.F90 modules.o iomod.o
potword22.o: ../potword22.F90 modules.o iomod.o
potword3.o: ../potword3.F90 modules.o iomod.o
potword4.o: ../potword4.F90 modules.o iomod.o
potword6.o: ../potword6.F90 modules.o iomod.o
potwordk.o: ../potwordk.F90 modules.o modules_kim.o iomod.o
potwordm.o: ../potwordm.F90 modules.o iomod.o
pr_integrator.o: ../pr_integrator.F90 modules.o m_pr.o m_trace.o pr_random.o datatypes.o
pr_random.o: ../pr_random.F90 datatypes.o modules.o
predict.o: ../predict.F90 modules.o iomod.o m_trace.o
pressure.o: ../pressure.F90 modules.o m_trace.o
projd0.o: ../projd0.F90 datatypes.o m_trace.o
projd0d.o: ../projd0d.F90 datatypes.o m_trace.o
projd0r.o: ../projd0r.F90 datatypes.o m_trace.o
projd3a.o: ../projd3a.F90 modules.o datatypes.o m_trace.o
projd3ad.o: ../projd3ad.F90 modules.o datatypes.o m_trace.o
projd3dV.o: ../projd3dV.F90 modules.o datatypes.o m_trace.o
projd3dVd.o: ../projd3dVd.F90 modules.o datatypes.o m_trace.o
projd3dVg.o: ../projd3dVg.F90 modules.o datatypes.o m_trace.o
projd3dVgd.o: ../projd3dVgd.F90 modules.o datatypes.o m_trace.o
projfork0.o: ../projfork0.F90 modules.o m_trace.o
projfork0r.o: ../projfork0r.F90 datatypes.o modules.o m_trace.o
projfork3.o: ../projfork3.F90 modules.o m_trace.o
projfork3d.o: ../projfork3d.F90 modules.o m_trace.o
projfork3dV.o: ../projfork3dV.F90 modules.o datatypes.o m_trace.o
projfork3dVg.o: ../projfork3dVg.F90 modules.o datatypes.o m_trace.o
projfork3dd.o: ../projfork3dd.F90 modules.o m_trace.o
projfork3ddV.o: ../projfork3ddV.F90 modules.o m_trace.o
projfork3dm.o: ../projfork3dm.F90 modules.o m_trace.o
projfork3ds.o: ../projfork3ds.F90 modules.o m_trace.o
projfork3s.o: ../projfork3s.F90 modules.o datatypes.o m_trace.o
projthbk0.o: ../projthbk0.F90 modules.o m_trace.o
projthbk0r.o: ../projthbk0r.F90 datatypes.o modules.o m_trace.o
projthbk3.o: ../projthbk3.F90 modules.o m_trace.o
projthbk3d.o: ../projthbk3d.F90 modules.o m_trace.o
projthbk3dV.o: ../projthbk3dV.F90 modules.o datatypes.o m_trace.o
projthbk3dVg.o: ../projthbk3dVg.F90 modules.o datatypes.o m_trace.o
projthbk3ddV.o: ../projthbk3ddV.F90 modules.o m_trace.o
projthbk3dm.o: ../projthbk3dm.F90 modules.o m_trace.o
projthbk3ds.o: ../projthbk3ds.F90 modules.o m_trace.o
projthbk3s.o: ../projthbk3s.F90 modules.o datatypes.o m_trace.o
property0.o: ../property0.F90 modules.o gulpcml.o gulpcml_props.o iomod.o m_trace.o
property12.o: ../property12.F90 modules.o iomod.o m_trace.o
property3.o: ../property3.F90 modules.o gulpcml.o gulpcml_props.o iomod.o m_trace.o
psibaskes.o: ../psibaskes.F90 modules.o m_trace.o
psirho.o: ../psirho.F90 modules.o m_trace.o
psiscreenderv.o: ../psiscreenderv.F90 modules.o m_trace.o m_strain.o
psiscreendervs.o: ../psiscreendervs.F90 modules.o m_trace.o m_strain.o
psumall.o: ../psumall.F90 modules.o m_trace.o
putpar.o: ../putpar.F90 modules.o m_trace.o
pythag.o: ../pythag.F90 datatypes.o
qeqbody.o: ../qeqbody.F90 modules.o m_trace.o
qmatrix1D.o: ../qmatrix1D.F90 modules.o m_trace.o
qmatrixelement.o: ../qmatrixelement.F90 modules.o m_trace.o
qmatrixelementc.o: ../qmatrixelementc.F90 datatypes.o modules.o m_trace.o
qmatrixelemente.o: ../qmatrixelemente.F90 datatypes.o modules.o m_trace.o
qtmin.o: ../qtmin.F90 modules.o iomod.o m_conjgr.o m_trace.o
qtpie.o: ../qtpie.F90 modules.o iomod.o m_trace.o
qupdate.o: ../qupdate.F90 modules.o m_trace.o
radialforce.o: ../radialforce.F90 modules.o m_trace.o
radialforcetrial.o: ../radialforcetrial.F90 modules.o m_trace.o
raman0.o: ../raman0.F90 modules.o datatypes.o iomod.o m_trace.o
raman3.o: ../raman3.F90 modules.o datatypes.o iomod.o m_trace.o
ramanstrength.o: ../ramanstrength.F90 modules.o iomod.o m_trace.o
ramantrm.o: ../ramantrm.F90 modules.o m_trace.o
random.o: ../random.F90 datatypes.o modules.o
ratfn.o: ../ratfn.F90 datatypes.o m_trace.o
real0d.o: ../real0d.F90 modules.o iomod.o m_trace.o
real0d3.o: ../real0d3.F90 modules.o m_trace.o
real0d3d.o: ../real0d3d.F90 modules.o m_trace.o
real0dd.o: ../real0dd.F90 modules.o iomod.o m_trace.o
real11.o: ../real11.F90 modules.o m_trace.o
real112.o: ../real112.F90 modules.o m_trace.o
real112d.o: ../real112d.F90 modules.o m_trace.o
real11d.o: ../real11d.F90 modules.o m_trace.o
real11pot.o: ../real11pot.F90 modules.o m_trace.o
real11s.o: ../real11s.F90 modules.o m_trace.o
real11s2.o: ../real11s2.F90 modules.o m_trace.o
real12a.o: ../real12a.F90 modules.o m_trace.o
real12a2.o: ../real12a2.F90 modules.o m_trace.o
real12a2d.o: ../real12a2d.F90 modules.o m_trace.o
real12ad.o: ../real12ad.F90 modules.o m_trace.o
real12apot.o: ../real12apot.F90 modules.o m_trace.o
real12as.o: ../real12as.F90 modules.o m_trace.o
real12as2.o: ../real12as2.F90 modules.o m_trace.o
real1D.o: ../real1D.F90 modules.o iomod.o m_trace.o m_strain.o
real1Dd.o: ../real1Dd.F90 modules.o iomod.o m_trace.o m_strain.o
real1Dfd.o: ../real1Dfd.F90 modules.o m_trace.o m_strain.o
real1Dmc.o: ../real1Dmc.F90 modules.o iomod.o m_trace.o
real1Dmd.o: ../real1Dmd.F90 modules.o iomod.o m_trace.o m_strain.o
real1Dp.o: ../real1Dp.F90 modules.o m_trace.o
real1Dpd.o: ../real1Dpd.F90 modules.o m_trace.o
real1Dq.o: ../real1Dq.F90 modules.o m_trace.o
real1d3.o: ../real1d3.F90 modules.o m_trace.o
real1d3dV.o: ../real1d3dV.F90 modules.o m_trace.o
real1d3dVg.o: ../real1d3dVg.F90 modules.o m_trace.o
real1d3g.o: ../real1d3g.F90 modules.o m_trace.o
real1fc.o: ../real1fc.F90 modules.o datatypes.o m_trace.o
real1fc3.o: ../real1fc3.F90 modules.o m_trace.o
real2Dq.o: ../real2Dq.F90 modules.o m_trace.o m_strain.o
real2a.o: ../real2a.F90 modules.o iomod.o m_trace.o
real2a1.o: ../real2a1.F90 modules.o m_trace.o
real2a1d.o: ../real2a1d.F90 modules.o m_trace.o
real2a1f.o: ../real2a1f.F90 modules.o m_trace.o
real2a1fd.o: ../real2a1fd.F90 modules.o m_trace.o
real2ad.o: ../real2ad.F90 modules.o iomod.o m_trace.o
real2adef.o: ../real2adef.F90 modules.o m_trace.o
real2aef.o: ../real2aef.F90 modules.o m_trace.o
real2amany.o: ../real2amany.F90 modules.o m_trace.o
real2b.o: ../real2b.F90 modules.o m_trace.o
real4.o: ../real4.F90 modules.o datatypes.o m_trace.o
realcmm.o: ../realcmm.F90 modules.o iomod.o m_trace.o
reale.o: ../reale.F90 modules.o datatypes.o iomod.o m_trace.o m_strain.o
realed.o: ../realed.F90 modules.o datatypes.o iomod.o m_trace.o m_strain.o
realfc.o: ../realfc.F90 modules.o datatypes.o m_trace.o
realfcd.o: ../realfcd.F90 modules.o datatypes.o m_trace.o
realfd0.o: ../realfd0.F90 modules.o m_trace.o
realfd3.o: ../realfd3.F90 modules.o datatypes.o m_trace.o m_strain.o
reallocate.o: ../reallocate.F90 datatypes.o iomod.o
realmc0.o: ../realmc0.F90 modules.o m_trace.o
realmc3.o: ../realmc3.F90 modules.o datatypes.o m_trace.o
realmd0.o: ../realmd0.F90 modules.o iomod.o m_trace.o
realmd0s.o: ../realmd0s.F90 modules.o datatypes.o iomod.o m_trace.o
realmd3.o: ../realmd3.F90 modules.o datatypes.o iomod.o m_trace.o m_strain.o
realmd3s.o: ../realmd3s.F90 modules.o datatypes.o iomod.o m_trace.o m_strain.o
realmd3t.o: ../realmd3t.F90 modules.o datatypes.o iomod.o m_trace.o m_strain.o
realmi3.o: ../realmi3.F90 modules.o iomod.o m_trace.o m_strain.o
realmi3d.o: ../realmi3d.F90 modules.o iomod.o m_trace.o m_strain.o
realmi3fd.o: ../realmi3fd.F90 modules.o m_trace.o m_strain.o
realmimc3.o: ../realmimc3.F90 modules.o m_trace.o
realp.o: ../realp.F90 modules.o datatypes.o m_trace.o
realpd.o: ../realpd.F90 modules.o datatypes.o m_trace.o
realrecip3d3.o: ../realrecip3d3.F90 modules.o datatypes.o iomod.o m_trace.o
realrecip3d3d.o: ../realrecip3d3d.F90 modules.o datatypes.o iomod.o m_trace.o
realrecip3d3dV.o: ../realrecip3d3dV.F90 modules.o datatypes.o iomod.o m_trace.o
realrecip3d3dVd.o: ../realrecip3d3dVd.F90 modules.o datatypes.o iomod.o m_trace.o
realrecip3d3dVg.o: ../realrecip3d3dVg.F90 modules.o datatypes.o iomod.o m_trace.o
realrecip3d3dVgd.o: ../realrecip3d3dVgd.F90 modules.o datatypes.o iomod.o m_trace.o
realrecipfc3.o: ../realrecipfc3.F90 modules.o datatypes.o m_trace.o
realsd.o: ../realsd.F90 modules.o datatypes.o iomod.o m_trace.o m_strain.o
realsd2.o: ../realsd2.F90 modules.o datatypes.o iomod.o m_trace.o m_strain.o
realselffc.o: ../realselffc.F90 modules.o datatypes.o m_trace.o
realselffcd.o: ../realselffcd.F90 modules.o datatypes.o m_trace.o
reaxff1fc.o: ../reaxff1fc.F90 datatypes.o modules.o iomod.o m_trace.o
reaxff_bo.o: ../reaxff_bo.F90 datatypes.o modules.o iomod.o m_trace.o
reaxff_funcs.o: ../reaxff_funcs.F90 modules.o m_trace.o datatypes.o
reaxffcharges.o: ../reaxffcharges.F90 datatypes.o modules.o m_trace.o
reaxffctfunct.o: ../reaxffctfunct.F90 iomod.o modules.o m_trace.o
reaxffmd.o: ../reaxffmd.F90 datatypes.o modules.o iomod.o m_trace.o m_strain.o
reaxffpot.o: ../reaxffpot.F90 datatypes.o modules.o m_trace.o
rebospline.o: ../rebospline.F90 datatypes.o modules.o
recip12a.o: ../recip12a.F90 modules.o m_trace.o
recip12ad.o: ../recip12ad.F90 modules.o m_trace.o
recip12apot.o: ../recip12apot.F90 modules.o m_trace.o
recip12as.o: ../recip12as.F90 modules.o m_trace.o
recip2D.o: ../recip2D.F90 modules.o m_trace.o m_strain.o
recip2Dd.o: ../recip2Dd.F90 modules.o m_trace.o m_strain.o
recip2Dfd.o: ../recip2Dfd.F90 modules.o m_trace.o m_strain.o
recip2Dp.o: ../recip2Dp.F90 modules.o m_trace.o
recip2Dpd.o: ../recip2Dpd.F90 modules.o m_trace.o
recip2Dq.o: ../recip2Dq.F90 modules.o m_trace.o m_strain.o
recip2Dtrial.o: ../recip2Dtrial.F90 modules.o m_trace.o
recip3D.o: ../recip3D.F90 modules.o m_trace.o m_strain.o
recip3Dd.o: ../recip3Dd.F90 modules.o m_trace.o m_strain.o
recip3Dfd.o: ../recip3Dfd.F90 modules.o m_trace.o m_strain.o
recip3Dmd.o: ../recip3Dmd.F90 modules.o m_trace.o m_strain.o
recip3Dp.o: ../recip3Dp.F90 modules.o m_trace.o
recip3Dpd.o: ../recip3Dpd.F90 modules.o m_trace.o
recip3Dspme.o: ../recip3Dspme.F90 modules.o m_fft3D.o m_trace.o m_strain.o
recip3Dtrial.o: ../recip3Dtrial.F90 modules.o m_trace.o
recip4.o: ../recip4.F90 modules.o m_trace.o
recipfc3.o: ../recipfc3.F90 modules.o m_trace.o
recipsd.o: ../recipsd.F90 modules.o m_trace.o m_strain.o
recipsd2.o: ../recipsd2.F90 modules.o m_trace.o m_strain.o
reciptrmd3.o: ../reciptrmd3.F90 modules.o m_trace.o
reciptrmd3dV.o: ../reciptrmd3dV.F90 modules.o m_trace.o
reciptrmd3dVg.o: ../reciptrmd3dVg.F90 modules.o m_trace.o
reciptrmd3g.o: ../reciptrmd3g.F90 modules.o m_trace.o
reciptrmdp.o: ../reciptrmdp.F90 modules.o m_trace.o m_strain.o
reinitialise.o: ../reinitialise.F90 datatypes.o modules.o m_trace.o
remove.o: ../remove.F90 modules.o m_trace.o
reorder.o: ../reorder.F90 modules.o m_trace.o
resetcfgdefaults.o: ../resetcfgdefaults.F90 modules.o reallocate.o m_trace.o
resetrigidmol.o: ../resetrigidmol.F90 modules.o 
resetx0rv.o: ../resetx0rv.F90 modules.o m_trace.o
rfind.o: ../rfind.F90 modules.o m_trace.o
rfindboth.o: ../rfindboth.F90 modules.o m_trace.o
rfohess.o: ../rfohess.F90 modules.o m_trace.o
rfostep.o: ../rfostep.F90 modules.o iomod.o
rhoinv.o: ../rhoinv.F90 m_trace.o modules.o
rhtohex.o: ../rhtohex.F90 modules.o m_trace.o
rigidmoleculedrv.o: ../rigidmoleculedrv.F90 modules.o
rigidmoleculegv.o: ../rigidmoleculegv.F90 modules.o
rigidmoleculephon.o: ../rigidmoleculephon.F90 modules.o
rigidmoleculeprop.o: ../rigidmoleculeprop.F90 modules.o
rlist.o: ../rlist.F90 modules.o m_trace.o
rpair.o: ../rpair.F90 modules.o m_trace.o
rpair2.o: ../rpair2.F90 modules.o m_trace.o
rsearch1D.o: ../rsearch1D.F90 modules.o m_trace.o
rsearch2D.o: ../rsearch2D.F90 modules.o m_trace.o
rsearch3D.o: ../rsearch3D.F90 modules.o m_trace.o
rtclist.o: ../rtclist.F90 modules.o m_trace.o
rtlist.o: ../rtlist.F90 modules.o m_trace.o
runmc.o: ../runmc.F90 modules.o iomod.o m_trace.o
runmd.o: ../runmd.F90 modules.o gulpcml.o iomod.o m_plumed.o m_pr.o license.o m_trace.o
runneb.o: ../runneb.F90 modules.o iomod.o m_trace.o
runsync.o: ../runsync.F90 modules.o iomod.o m_trace.o
samemol.o: ../samemol.F90 modules.o m_trace.o
scaleevec.o: ../scaleevec.F90 datatypes.o m_trace.o
screenct.o: ../screenct.F90 modules.o m_trace.o
sec0.o: ../sec0.F90 modules.o iomod.o m_trace.o
sec0d.o: ../sec0d.F90 modules.o iomod.o m_trace.o
sec0f.o: ../sec0f.F90 modules.o iomod.o m_trace.o
sec3.o: ../sec3.F90 modules.o iomod.o m_trace.o
sec3d.o: ../sec3d.F90 modules.o iomod.o m_trace.o
sec3f.o: ../sec3f.F90 modules.o iomod.o m_trace.o
secondorderfc.o: ../secondorderfc.F90 modules.o iomod.o m_trace.o
secondorderfcd.o: ../secondorderfcd.F90 modules.o iomod.o m_trace.o
secondpass.o: ../secondpass.F90 modules.o iomod.o m_trace.o
selfterm.o: ../selfterm.F90 modules.o m_trace.o
selftermfc.o: ../selftermfc.F90 modules.o m_trace.o
selftermfcd.o: ../selftermfcd.F90 modules.o m_trace.o
selftermp.o: ../selftermp.F90 modules.o m_trace.o
sendall.o: ../sendall.F90 modules.o
set2a.o: ../set2a.F90 modules.o datatypes.o m_trace.o
setatomdistribution.o: ../setatomdistribution.F90 modules.o
setatomnoded2.o: ../setatomnoded2.F90 modules.o
setatomnodes.o: ../setatomnodes.F90 modules.o iomod.o
setatomnodesbo.o: ../setatomnodesbo.F90 datatypes.o modules.o iomod.o
setbelashchenko.o: ../setbelashchenko.F90 modules.o
setbondorder.o: ../setbondorder.F90 modules.o
setbondtype.o: ../setbondtype.F90 modules.o
setbrenner1.o: ../setbrenner1.F90 datatypes.o modules.o
setbrenner3.o: ../setbrenner3.F90 datatypes.o modules.o
setbsmptr.o: ../setbsmptr.F90 modules.o
setc6.o: ../setc6.F90 modules.o
setcellderv.o: ../setcellderv.F90 modules.o
setcfg.o: ../setcfg.F90 modules.o iomod.o modules_kim.o m_trace.o
setcmm.o: ../setcmm.F90 modules.o iomod.o
setcombine.o: ../setcombine.F90 modules.o
setcosmo.o: ../setcosmo.F90 modules.o
setcosmoainv.o: ../setcosmoainv.F90 modules.o
setcosmoainv2D.o: ../setcosmoainv2D.F90 modules.o
setcosmoamat.o: ../setcosmoamat.F90 modules.o iomod.o reallocate.o
setcosmoamat2D.o: ../setcosmoamat2D.F90 modules.o iomod.o reallocate.o
setcosmobmat.o: ../setcosmobmat.F90 modules.o iomod.o
setcutoffmax.o: ../setcutoffmax.F90 modules.o modules_kim.o
setcutoffmaxbo.o: ../setcutoffmaxbo.F90 modules.o
setdcosmoamat.o: ../setdcosmoamat.F90 modules.o iomod.o reallocate.o
setdcosmoamat1p.o: ../setdcosmoamat1p.F90 modules.o iomod.o reallocate.o
setdcosmobmat.o: ../setdcosmobmat.F90 modules.o reallocate.o
setdcosmobmat1p.o: ../setdcosmobmat1p.F90 modules.o reallocate.o
setdef.o: ../setdef.F90 modules.o reallocate.o m_trace.o
setdestroyptr.o: ../setdestroyptr.F90 modules.o reallocate.o
setdisp.o: ../setdisp.F90 modules.o
setdistance.o: ../setdistance.F90 modules.o
setdr2mat.o: ../setdr2mat.F90 modules.o iomod.o m_trace.o
setdsegweight.o: ../setdsegweight.F90 modules.o
seteam.o: ../seteam.F90 modules.o
setedip.o: ../setedip.F90 datatypes.o modules.o iomod.o
seteem.o: ../seteem.F90 modules.o m_trace.o
seteembond.o: ../seteembond.F90 modules.o iomod.o 
setele.o: ../setele.F90 modules.o m_pdfneutron.o m_trace.o
setewald.o: ../setewald.F90 modules.o iomod.o
setfeshtmat.o: ../setfeshtmat.F90 datatypes.o m_trace.o
setfour.o: ../setfour.F90 modules.o
setframe.o: ../setframe.F90 modules.o
setkeywrd.o: ../setkeywrd.F90 modules.o m_pdfneutron.o m_trace.o
setkim.o: ../setkim.F90 modules.o iomod.o kim_functions.o reallocate.o modules_kim.o datatypes.o
setkpt.o: ../setkpt.F90 modules.o iomod.o m_pdfneutron.o
setktrm4.o: ../setktrm4.F90 modules.o
setktrmd3.o: ../setktrmd3.F90 modules.o
setktrmdp.o: ../setktrmdp.F90 modules.o
setktrmfc3.o: ../setktrmfc3.F90 modules.o
setlib.o: ../setlib.F90 modules.o m_trace.o iomod.o
setlist3.o: ../setlist3.F90 modules.o
setlist4.o: ../setlist4.F90 modules.o
setlist6.o: ../setlist6.F90 modules.o
setlistimp.o: ../setlistimp.F90 modules.o
setlistoop.o: ../setlistoop.F90 modules.o
setmaxcell1D.o: ../setmaxcell1D.F90 modules.o
setmd.o: ../setmd.F90 modules.o iomod.o m_plumed.o m_pr.o pr_random.o m_trace.o
setmol.o: ../setmol.F90 modules.o reallocate.o
setmola.o: ../setmola.F90 modules.o reallocate.o
setmoldef.o: ../setmoldef.F90 modules.o reallocate.o
setmoldef12a.o: ../setmoldef12a.F90 modules.o
setmols.o: ../setmols.F90 modules.o reallocate.o
setmoveptr.o: ../setmoveptr.F90 modules.o reallocate.o
setoccptr.o: ../setoccptr.F90 modules.o
setoccptrs.o: ../setoccptrs.F90 modules.o
setoccshptr.o: ../setoccshptr.F90 modules.o
setoccshptrs.o: ../setoccshptrs.F90 modules.o
setoptptr.o: ../setoptptr.F90 modules.o
setphonptr.o: ../setphonptr.F90 modules.o
setpointgroup.o: ../setpointgroup.F90 datatypes.o
setpolar.o: ../setpolar.F90 modules.o
setpote.o: ../setpote.F90 modules.o
setptsnoded2.o: ../setptsnoded2.F90 modules.o reallocate.o
setqonsas.o: ../setqonsas.F90 modules.o iomod.o
setqonsas2D.o: ../setqonsas2D.F90 modules.o iomod.o
setquaternion.o: ../setquaternion.F90 modules.o m_trace.o
setquasiharm.o: ../setquasiharm.F90 modules.o
setratiom.o: ../setratiom.F90 modules.o
setreaxFF.o: ../setreaxFF.F90 datatypes.o modules.o iomod.o
setreaxffQ.o: ../setreaxffQ.F90 modules.o iomod.o
setreaxffQ1fc.o: ../setreaxffQ1fc.F90 modules.o iomod.o
setreaxffQiter.o: ../setreaxffQiter.F90 modules.o iomod.o reallocate.o
setreaxffQiter1fc.o: ../setreaxffQiter1fc.F90 modules.o iomod.o reallocate.o
setregiontrans.o: ../setregiontrans.F90 modules.o
setrigidmol.o: ../setrigidmol.F90 modules.o 
setrotateptr.o: ../setrotateptr.F90 modules.o reallocate.o
setrotation.o: ../setrotation.F90 modules.o m_trace.o
setsas.o: ../setsas.F90 modules.o iomod.o m_trace.o
setsegweight.o: ../setsegweight.F90 modules.o
setsix.o: ../setsix.F90 modules.o
setslice.o: ../setslice.F90 modules.o
setspatial.o: ../setspatial.F90 modules.o iomod.o
setspatialbo.o: ../setspatialbo.F90 modules.o iomod.o
setspatialcells.o: ../setspatialcells.F90 datatypes.o
setspec.o: ../setspec.F90 modules.o
setspring.o: ../setspring.F90 modules.o
setstrainptr.o: ../setstrainptr.F90 modules.o reallocate.o
setswapptr.o: ../setswapptr.F90 modules.o reallocate.o
settaper.o: ../settaper.F90 modules.o
setthree.o: ../setthree.F90 modules.o
settings.o: ../settings.F90 datatypes.o modules.o
settorcubic.o: ../settorcubic.F90 datatypes.o modules.o
settricubic.o: ../settricubic.F90 datatypes.o modules.o
settricubic2.o: ../settricubic2.F90 datatypes.o modules.o
setuff.o: ../setuff.F90 modules.o
setup.o: ../setup.F90 modules.o datatypes.o iomod.o reallocate.o m_trace.o m_strain.o
setvarnoded2.o: ../setvarnoded2.F90 modules.o reallocate.o
setvbon.o: ../setvbon.F90 modules.o iomod.o
setzmolslice.o: ../setzmolslice.F90 modules.o
shellgroupvelocities.o: ../shellgroupvelocities.F90 modules.o iomod.o m_trace.o
shellopt.o: ../shellopt.F90 modules.o iomod.o reallocate.o m_trace.o
shmin.o: ../shmin.F90 modules.o iomod.o m_conjgr.o m_trace.o
sitepot.o: ../sitepot.F90 modules.o iomod.o m_trace.o
six1fc.o: ../six1fc.F90 modules.o m_trace.o
sixbody.o: ../sixbody.F90 modules.o iomod.o m_trace.o
sixfc.o: ../sixfc.F90 modules.o m_trace.o
sixfc3.o: ../sixfc3.F90 modules.o m_trace.o
sixfcd.o: ../sixfcd.F90 modules.o m_trace.o
sixfd.o: ../sixfd.F90 modules.o m_trace.o m_strain.o
sixlist.o: ../sixlist.F90 modules.o iomod.o m_trace.o m_strain.o
sixmc.o: ../sixmc.F90 modules.o iomod.o m_trace.o
sixmd.o: ../sixmd.F90 modules.o iomod.o m_trace.o m_strain.o
sixmds.o: ../sixmds.F90 modules.o iomod.o m_trace.o m_strain.o
sixnos.o: ../sixnos.F90 modules.o iomod.o m_trace.o m_strain.o
sixnosd.o: ../sixnosd.F90 modules.o iomod.o m_trace.o m_strain.o
sixp.o: ../sixp.F90 modules.o m_trace.o
sixpd.o: ../sixpd.F90 modules.o m_trace.o
sixsd.o: ../sixsd.F90 modules.o iomod.o m_trace.o m_strain.o
sixsd2.o: ../sixsd2.F90 modules.o iomod.o m_trace.o m_strain.o
smbody.o: ../smbody.F90 modules.o m_trace.o
solvation.o: ../solvation.F90 modules.o m_trace.o
sort.o: ../sort.F90 modules.o m_trace.o
sort1.o: ../sort1.F90 modules.o reallocate.o m_trace.o
sort1sym.o: ../sort1sym.F90 modules.o m_trace.o
sort2a.o: ../sort2a.F90 modules.o m_trace.o
special.o: ../special.F90 modules.o iomod.o m_trace.o
specialcom.o: ../specialcom.F90 modules.o iomod.o m_trace.o
specialqtn.o: ../specialqtn.F90 modules.o iomod.o m_trace.o
speclibw.o: ../speclibw.F90 modules.o m_trace.o iomod.o
specword.o: ../specword.F90 modules.o iomod.o
spline.o: ../spline.F90 modules.o datatypes.o iomod.o
stolc.o: ../stolc.F90 datatypes.o
stopnow.o: ../stopnow.F90 iomod.o modules.o
stouc.o: ../stouc.F90 datatypes.o
strain.o: ../strain.F90 datatypes.o m_trace.o
strfin.o: ../strfin.F90 modules.o iomod.o m_trace.o
strsym.o: ../strsym.F90 modules.o m_trace.o
strword.o: ../strword.F90 modules.o iomod.o reallocate.o
sumall.o: ../sumall.F90 modules.o m_trace.o
sumderv2.o: ../sumderv2.F90 modules.o m_trace.o
sumderv2f.o: ../sumderv2f.F90 modules.o m_trace.o
sumderv2p.o: ../sumderv2p.F90 modules.o m_trace.o
sumderv2s.o: ../sumderv2s.F90 modules.o m_trace.o
super.o: ../super.F90 modules.o m_trace.o
surfaceenergy.o: ../surfaceenergy.F90 modules.o m_trace.o
surfword.o: ../surfword.F90 modules.o reallocate.o iomod.o
switch.o: ../switch.F90 modules.o m_trace.o
symbind.o: ../symbind.F90 datatypes.o m_trace.o
symderv2.o: ../symderv2.F90 modules.o m_trace.o
symdervrot.o: ../symdervrot.F90 modules.o m_trace.o
symgenkpt.o: ../symgenkpt.F90 modules.o m_trace.o
symmet.o: ../symmet.F90 modules.o iomod.o m_trace.o
symmetp.o: ../symmetp.F90 modules.o iomod.o m_trace.o
symoff.o: ../symoff.F90 modules.o m_trace.o
symout.o: ../symout.F90 modules.o iomod.o m_trace.o
sympot.o: ../sympot.F90 modules.o m_trace.o
sympotd.o: ../sympotd.F90 modules.o m_trace.o
symprod.o: ../symprod.F90 iomod.o modules.o m_trace.o
syncforce.o: ../syncforce.F90 modules.o m_trace.o
tangtoennies.o: ../tangtoennies.F90 datatypes.o modules.o m_trace.o
taper.o: ../taper.F90 datatypes.o modules.o m_trace.o
testspacegroup.o: ../testspacegroup.F90 modules.o m_trace.o
tfunc.o: ../tfunc.F90 datatypes.o m_trace.o
thermalconductivity.o: ../thermalconductivity.F90 modules.o iomod.o m_trace.o
thermalconductivityd.o: ../thermalconductivityd.F90 modules.o iomod.o m_trace.o
thirdorderfc3.o: ../thirdorderfc3.F90 modules.o datatypes.o iomod.o m_trace.o
thirdordern.o: ../thirdordern.F90 modules.o iomod.o m_trace.o
three0d3.o: ../three0d3.F90 modules.o m_trace.o
three12.o: ../three12.F90 modules.o m_trace.o
three12d.o: ../three12d.F90 modules.o m_trace.o
three1fc.o: ../three1fc.F90 modules.o m_trace.o
three3d3.o: ../three3d3.F90 modules.o m_trace.o
three3d3dV.o: ../three3d3dV.F90 modules.o m_trace.o
three3d3dVg.o: ../three3d3dVg.F90 modules.o m_trace.o
three3d3g.o: ../three3d3g.F90 modules.o m_trace.o
threebody.o: ../threebody.F90 modules.o m_trace.o
threefc.o: ../threefc.F90 modules.o m_trace.o
threefc3.o: ../threefc3.F90 modules.o m_trace.o
threefcd.o: ../threefcd.F90 modules.o m_trace.o
threefd.o: ../threefd.F90 modules.o m_trace.o m_strain.o
threelist.o: ../threelist.F90 modules.o iomod.o m_trace.o m_strain.o
threemc.o: ../threemc.F90 modules.o iomod.o m_trace.o
threemd.o: ../threemd.F90 modules.o iomod.o m_trace.o m_strain.o
threemds.o: ../threemds.F90 modules.o iomod.o m_trace.o m_strain.o
threenos.o: ../threenos.F90 modules.o iomod.o m_trace.o m_strain.o
threenosd.o: ../threenosd.F90 modules.o iomod.o m_trace.o m_strain.o
threep.o: ../threep.F90 modules.o m_trace.o
threepd.o: ../threepd.F90 modules.o m_trace.o
threesd.o: ../threesd.F90 modules.o datatypes.o iomod.o m_trace.o m_strain.o
threesd2.o: ../threesd2.F90 modules.o iomod.o m_trace.o m_strain.o
transmat.o: ../transmat.F90 modules.o iomod.o m_trace.o
transmatd.o: ../transmatd.F90 modules.o iomod.o m_trace.o
transmatmol.o: ../transmatmol.F90 modules.o iomod.o m_trace.o
trap_sigint.o: ../trap_sigint.F90 datatypes.o modules.o iomod.o
twobody.o: ../twobody.F90 modules.o m_trace.o
twobodymd.o: ../twobodymd.F90 modules.o m_trace.o
twoden.o: ../twoden.F90 modules.o m_trace.o
uncell.o: ../uncell.F90 modules.o m_trace.o
uncentre.o: ../uncentre.F90 modules.o datatypes.o m_trace.o
unfreeze.o: ../unfreeze.F90 modules.o m_trace.o
updatesas.o: ../updatesas.F90 modules.o iomod.o m_trace.o
validtriad.o: ../validtriad.F90 modules.o m_trace.o
vartocfg.o: ../vartocfg.F90 modules.o m_trace.o
vibenergy.o: ../vibenergy.F90 modules.o iomod.o m_trace.o
volume.o: ../volume.F90 datatypes.o
wcoord.o: ../wcoord.F90 datatypes.o modules.o m_trace.o iomod.o
wolfself.o: ../wolfself.F90 modules.o m_trace.o
wolfselftrial.o: ../wolfselftrial.F90 modules.o m_trace.o
worsy.o: ../worsy.F90 modules.o iomod.o
wtof.o: ../wtof.F90 datatypes.o
x0strain.o: ../x0strain.F90 modules.o m_trace.o m_strain.o
x0tostr.o: ../x0tostr.F90 modules.o m_trace.o
x0tostrcentroid.o: ../x0tostrcentroid.F90 modules.o m_trace.o
x0tostrtrial.o: ../x0tostrtrial.F90 modules.o m_trace.o
xctox0.o: ../xctox0.F90 modules.o m_trace.o
