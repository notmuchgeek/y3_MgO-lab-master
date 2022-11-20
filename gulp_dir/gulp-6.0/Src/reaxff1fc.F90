  subroutine reaxFF1fc(maxlhs,d1cell,matom,vmatom)
!
!  Routine for calculating the first derivatives of reaxFF potentials
!  for use in a finite difference unphased phonon calculation.
!
!  10/15 Created from reaxffmd
!   6/17 nboatomptr passed to setreaxffQiter
!   6/17 Dummy pointer added to pass to derivative calls
!   6/17 Module files renamed to gulp_files
!  11/17 Calculation of energy terms removed as they are not needed
!  12/17 Parallelisation corrected
!  12/17 Bug due to setting of xdrv/ydrv/zdrv removed
!  12/17 Spatial decomposition removed for Hbonding/Coulomb/VDW for 
!        now until algorithm is checked.
!  12/17 Coulomb derivatives now fixed
!   2/18 Trace added
!   3/19 Lone pair smoothing added to ReaxFF
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!   9/19 ReaxFF charge arrays moved to reaxFFdata
!   3/20 Location of angstoev changed to current
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
!  Julian Gale, CIC, Curtin University, March 2020
!
  use datatypes
  use configurations, only : nregionno, nregiontype, QMMMmode
  use control,        only : keyword, lreaxFFQ, literativeQ
  use current
  use derivatives,    only : nd2central, nd2cell, nd2cellptr
  use gulp_files,     only : lqbo
  use iochannels
  use neighbours
  use numbers,        only : third
  use parallel,       only : ioproc
  use reaxFFdata
  use realvectors,    only : dist, xtmp, ytmp, ztmp, cellindex
  use spatialbo
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4),                       intent(in)    :: maxlhs
  integer(i4),                       intent(in)    :: matom
  real(dp),                          intent(inout) :: d1cell(4,maxlhs,*)
  real(dp),                          intent(in)    :: vmatom(4)
!
!  Local variables
!
  integer(i4)                                      :: i
  integer(i4)                                      :: ii
  integer(i4)                                      :: imin
  integer(i4)                                      :: ind
  integer(i4)                                      :: ind1
  integer(i4)                                      :: ind2
  integer(i4)                                      :: indil
  integer(i4)                                      :: indjk
  integer(i4)                                      :: indjl
  integer(i4)                                      :: indkl
  integer(i4), dimension(:,:,:), allocatable, save :: ineigh         ! Cell indices for vector
  integer(i4)                                      :: itmp
  integer(i4)                                      :: j
  integer(i4)                                      :: jj
  integer(i4)                                      :: jmax
  integer(i4)                                      :: k
  integer(i4)                                      :: kmax
  integer(i4)                                      :: kmin
  integer(i4)                                      :: l
  integer(i4)                                      :: m
  integer(i4)                                      :: maxneigh1
  integer(i4)                                      :: maxneigh2
  integer(i4)                                      :: maxneigh2i
  integer(i4)                                      :: mneigh
  integer(i4)                                      :: n
  integer(i4)                                      :: nati
  integer(i4)                                      :: nboij
  integer(i4)                                      :: nbojk
  integer(i4), dimension(:),     allocatable, save :: nbos
  integer(i4), dimension(:),     allocatable, save :: nbosptr
  integer(i4)                                      :: nboatom
  integer(i4), dimension(:),     allocatable, save :: nboatomptr
  integer(i4), dimension(:),     allocatable, save :: nboatomRptr
  integer(i4), dimension(:),     allocatable, save :: nbodummyptr
  integer(i4)                                      :: ncindm
  integer(i4)                                      :: ncindmij
  integer(i4)                                      :: ncindmik
  integer(i4)                                      :: ncindmjk
  integer(i4)                                      :: ncindp
  integer(i4)                                      :: ncindpij
  integer(i4)                                      :: ncindpik
  integer(i4)                                      :: ncindpjk
  integer(i4)                                      :: ni
  integer(i4)                                      :: nivalid
  integer(i4)                                      :: nj
  integer(i4)                                      :: nk
  integer(i4)                                      :: nl
  integer(i4)                                      :: nmin
  integer(i4)                                      :: nmolonly
  integer(i4)                                      :: nn
  integer(i4)                                      :: nn2
  integer(i4)                                      :: nor
  integer(i4)                                      :: nptr
  integer(i4), dimension(:,:),   allocatable, save :: neighno
  integer(i4), dimension(:,:),   allocatable, save :: neighnoRptr
  integer(i4), dimension(:),     allocatable, save :: nneigh
  integer(i4)                                      :: nneighi1
  integer(i4)                                      :: nneighj1
  integer(i4)                                      :: nneighi2
  integer(i4)                                      :: nneighj2
  integer(i4)                                      :: nneighi2i
  integer(i4)                                      :: nregioni
  integer(i4)                                      :: nregionj
  integer(i4)                                      :: nregionk
  integer(i4)                                      :: nregiontypi
  integer(i4)                                      :: nregiontypj
  integer(i4)                                      :: nspeci
  integer(i4)                                      :: nspecj
  integer(i4)                                      :: nspeck
  integer(i4)                                      :: nspecl
  integer(i4)                                      :: ntypi
  integer(i4)                                      :: numattodo
  integer(i4)                                      :: numattodo_main
  integer(i4), dimension(:),     allocatable, save :: numattodoptr
  integer(i4), dimension(:),     allocatable, save :: numattodoRptr
  integer(i4)                                      :: nval3
  integer(i4)                                      :: status
  logical                                          :: lAnyHB
  logical,     dimension(:),     allocatable, save :: latomdone
  logical                                          :: lanyBOnonzero
  logical                                          :: lBOijOK
  logical                                          :: lBOikOK
  logical                                          :: lBOjkOK
  logical                                          :: lBOnonzero
  logical                                          :: lDoEpen
  logical                                          :: lfixQi
  logical                                          :: lfixQj
  logical                                          :: lfound
  logical                                          :: lmaxneighok
  logical                                          :: lQMMMelectro
  logical                                          :: lQMMMok
  logical                                          :: lreactivei
  logical                                          :: lreactivej
  logical                                          :: lreg2one
  logical                                          :: lreg2pair
  logical                                          :: lself
  logical                                          :: lStrict
  logical                                          :: ltoriGTj
  real(dp)                                         :: bo8
  real(dp)                                         :: BOij
  real(dp)                                         :: BOij_s
  real(dp)                                         :: BOij_pi
  real(dp)                                         :: BOij_pipi
  real(dp)                                         :: BOijpbe2
  real(dp)                                         :: BOik
  real(dp)                                         :: BOjk
  real(dp)                                         :: BOjl
  real(dp)                                         :: BOpij
  real(dp)                                         :: conjtrm1
  real(dp)                                         :: conjtrm2
  real(dp)                                         :: coshalftheta
  real(dp)                                         :: cosp1d(6)
  real(dp)                                         :: cosp2d(21)
  real(dp)                                         :: cosphi
  real(dp)                                         :: cos2phi
  real(dp)                                         :: cos3phi
  real(dp)                                         :: g_cpu_time
  real(dp)                                         :: cut2
  real(dp)                                         :: cutq2
  real(dp)                                         :: cutvdw2
  real(dp)                                         :: dcos2phidcosphi
  real(dp)                                         :: dcos3phidcosphi
  real(dp)                                         :: delta_angle
  real(dp)                                         :: ddeltalpddeltae
  real(dp)                                         :: ddeltalpddeltaelp
  real(dp)                                         :: ddeltalpcorrddeltai
  real(dp)                                         :: ddeltalpcorrjddeltalpj
  real(dp)                                         :: ddeltalpcorrjddeltaj
  real(dp)                                         :: deoverddeltalpcorr_j
  real(dp)                                         :: dehbdBOij
  real(dp)                                         :: dehbdrik
  real(dp)                                         :: dehbdrjk
  real(dp)                                         :: dehbdtheta
  real(dp)                                         :: debpenddeltai
  real(dp)                                         :: deloneddeltai
  real(dp)                                         :: deloneddeltalp
  real(dp)                                         :: deoverdBOij
  real(dp)                                         :: deoverdBOpi
  real(dp)                                         :: deoverddeltai
  real(dp)                                         :: deoverddeltalpcorr
  real(dp)                                         :: deconjdcosphi
  real(dp)                                         :: deconjdf12
  real(dp)                                         :: deconjdsinkij
  real(dp)                                         :: deconjdsinijl
  real(dp)                                         :: detorsdcosphi
  real(dp)                                         :: detorsdcos1phi
  real(dp)                                         :: detorsdcos2phi
  real(dp)                                         :: detorsdcos3phi
  real(dp)                                         :: detorsddeltai
  real(dp)                                         :: detorsdexpV2
  real(dp)                                         :: detorsdf10
  real(dp)                                         :: detorsdsinkij
  real(dp)                                         :: detorsdsinijl
  real(dp)                                         :: deunderddeltalpcorr
  real(dp)                                         :: deunderddeltalpcorr_j
  real(dp)                                         :: deunderddeltai
  real(dp)                                         :: deunderutrm2deltaBO
  real(dp)                                         :: deunderutrm2ddeltaj
  real(dp)                                         :: deunderdBOpi
  real(dp)                                         :: devalddeltai
  real(dp)                                         :: devaldsbo
  real(dp)                                         :: devdwdrij
  real(dp)                                         :: diffBOdelta
  real(dp)                                         :: dtheta0dsbo
  real(dp)                                         :: dthetadrij
  real(dp)                                         :: dthetadrik
  real(dp)                                         :: dthetadrjk
  real(dp)                                         :: devaldtheta
  real(dp)                                         :: dsboproddbo
  real(dp)                                         :: dsboddeltai
  real(dp)                                         :: ebpen
  real(dp)                                         :: ebpenij
  real(dp)                                         :: ecoatrm
  real(dp)                                         :: epentrm
  real(dp)                                         :: eself
  real(dp)                                         :: etorsconjtrm
  real(dp)                                         :: etrm
  real(dp)                                         :: exptheta
  real(dp),    dimension(:,:),   allocatable, save :: BO
  real(dp),    dimension(:,:),   allocatable, save :: BO_pi
  real(dp),    dimension(:,:),   allocatable, save :: BO_pipi
  real(dp),    dimension(:,:),   allocatable, save :: BOp
  real(dp),    dimension(:,:),   allocatable, save :: BOp_pi
  real(dp),    dimension(:,:),   allocatable, save :: BOp_pipi
  real(dp),    dimension(:,:),   allocatable, save :: d1BOp
  real(dp),    dimension(:,:),   allocatable, save :: d1BOp_pi
  real(dp),    dimension(:,:),   allocatable, save :: d1BOp_pipi
  real(dp),    dimension(:),     allocatable, save :: d1i
  real(dp),    dimension(:),     allocatable, save :: d1j
  real(dp),    dimension(:),     allocatable, save :: d1BOi
  real(dp),    dimension(:),     allocatable, save :: d1BOi_s
  real(dp),    dimension(:),     allocatable, save :: d1BOi_pi
  real(dp),    dimension(:),     allocatable, save :: d1BOi_pipi
  real(dp),    dimension(:),     allocatable, save :: d1BOj
  real(dp),    dimension(:),     allocatable, save :: d1BOj_s
  real(dp),    dimension(:),     allocatable, save :: d1BOj_pi
  real(dp),    dimension(:),     allocatable, save :: d1BOj_pipi
  real(dp),    dimension(:),     allocatable, save :: delta
  real(dp),    dimension(:),     allocatable, save :: deltap
  real(dp),    dimension(:),     allocatable, save :: deltalplp ! delta lp for lone pairs
  real(dp),    dimension(:),     allocatable, save :: deltalp   ! delta lp for over coordination
  real(dp),    dimension(:),     allocatable, save :: deoverunderddeltaj_neigh
  real(dp),    dimension(:),     allocatable, save :: debpendBO_neigh
  real(dp),    dimension(:),     allocatable, save :: devaldBO_neigh
  real(dp),    dimension(:),     allocatable, save :: detorsdBO_neigh
  real(dp),    dimension(:),     allocatable, save :: detorsdBOpi_neigh
  real(dp),    dimension(:),     allocatable, save :: devaldsbo_neigh
  real(dp)                                         :: dBOpijdr
  real(dp)                                         :: deltapi
  real(dp)                                         :: deltapj
  real(dp)                                         :: delta_coa
  real(dp)                                         :: delta_e
  real(dp)                                         :: delta_e0
  real(dp)                                         :: deltalp_corr_i
  real(dp)                                         :: deltalp_corr_j
  real(dp)                                         :: deijdBO_s
  real(dp)                                         :: deijdBO_pi
  real(dp)                                         :: deijdBO_pipi
  real(dp)                                         :: dexpV2dBO_pi
  real(dp)                                         :: dexpV2df11
  real(dp)                                         :: df7dBOij
  real(dp)                                         :: df7dBOik
  real(dp)                                         :: d2f7dBOij2
  real(dp)                                         :: d2f7dBOijdBOik
  real(dp)                                         :: d2f7dBOik2
  real(dp)                                         :: df8ddeltai
  real(dp)                                         :: d2f8ddeltai2
  real(dp)                                         :: df9ddeltai
  real(dp)                                         :: d2f9ddeltai2
  real(dp)                                         :: df10dBOij
  real(dp)                                         :: df10dBOik
  real(dp)                                         :: df10dBOjl
  real(dp)                                         :: d2f10dBOik2
  real(dp)                                         :: d2f10dBOikdBOij
  real(dp)                                         :: d2f10dBOikdBOjl
  real(dp)                                         :: d2f10dBOjl2
  real(dp)                                         :: df11ddeltaij
  real(dp)                                         :: d2f11ddeltaij2
  real(dp)                                         :: df12dBOij
  real(dp)                                         :: df12dBOik
  real(dp)                                         :: df12dBOjl
  real(dp)                                         :: d2f12dBOik2
  real(dp)                                         :: d2f12dBOikdBOij
  real(dp)                                         :: d2f12dBOikdBOjl
  real(dp)                                         :: d2f12dBOij2
  real(dp)                                         :: d2f12dBOijdBOjl
  real(dp)                                         :: d2f12dBOjl2
  real(dp)                                         :: df13drij
  real(dp)                                         :: d2f13drij2
  real(dp)                                         :: fij
  real(dp)                                         :: dfijdBOij
  real(dp)                                         :: d2fijdBOij2
  real(dp)                                         :: d3fijdBOij3
  real(dp)                                         :: fik
  real(dp)                                         :: dfikdBOik
  real(dp)                                         :: d2fikdBOik2
  real(dp)                                         :: d3fikdBOik3
  real(dp)                                         :: fjk
  real(dp)                                         :: dfjkdBOjk
  real(dp)                                         :: d2fjkdBOjk2
  real(dp)                                         :: d3fjkdBOjk3
  real(dp)                                         :: fjl
  real(dp)                                         :: dfjldBOjl
  real(dp)                                         :: d2fjldBOjl2
  real(dp)                                         :: d3fjldBOjl3
  real(dp)                                         :: fHB
  real(dp)                                         :: dfHBdBO
  real(dp)                                         :: d2fHBdBO2
  real(dp)                                         :: d3fHBdBO3
  real(dp)                                         :: frHB
  real(dp)                                         :: dfrHBdr
  real(dp)                                         :: d2frHBdr2
  real(dp)                                         :: d3frHBdr3
  real(dp)                                         :: dsinkijdrij
  real(dp)                                         :: dsinkijdrik
  real(dp)                                         :: dsinkijdrjk
  real(dp)                                         :: dsinijldrij
  real(dp)                                         :: dsinijldril
  real(dp)                                         :: dsinijldrjl
  real(dp)                                         :: d2sinkijdr2
  real(dp)                                         :: d2sinijldr2
  real(dp)                                         :: d2thetadr2
  real(dp)                                         :: d2theta0dsbo2
  real(dp)                                         :: dVDWtrmdf13
  real(dp)                                         :: dVDWtrmdrij
  real(dp)                                         :: expbo8
  real(dp)                                         :: expcoa2
  real(dp)                                         :: expcoa3a
  real(dp)                                         :: expcoa3b
  real(dp)                                         :: expcoa4a
  real(dp)                                         :: expcoa4b
  real(dp)                                         :: expcoaprod
  real(dp)                                         :: expcoatrm
  real(dp)                                         :: expdi
  real(dp)                                         :: expdi_j
  real(dp)                                         :: exphb1
  real(dp)                                         :: expH
  real(dp)                                         :: exphb2
  real(dp)                                         :: expij
  real(dp)                                         :: explp
  real(dp)                                         :: explp2
  real(dp)                                         :: expoc
  real(dp)                                         :: expoc_j
  real(dp)                                         :: exppen1
  real(dp)                                         :: exppen2
  real(dp)                                         :: expuc1
  real(dp)                                         :: expuc1_j
  real(dp)                                         :: expuc2
  real(dp)                                         :: expuc2_j
  real(dp)                                         :: expuc3
  real(dp)                                         :: expuc3_j
  real(dp)                                         :: expV2
  real(dp)                                         :: expVDW1
  real(dp)                                         :: expVDW2
  real(dp)                                         :: f7
  real(dp)                                         :: f8
  real(dp)                                         :: f9
  real(dp)                                         :: f10
  real(dp)                                         :: f11
  real(dp)                                         :: f12
  real(dp)                                         :: f13
  real(dp)                                         :: decouldrij
  real(dp)                                         :: dgamdrij
  real(dp)                                         :: gam
  real(dp)                                         :: gammai
  real(dp)                                         :: gammaj
  real(dp)                                         :: gammaij
  real(dp)                                         :: otrm1
  real(dp)                                         :: otrm1_j
  real(dp)                                         :: otrm2
  real(dp)                                         :: otrm2_j
  real(dp)                                         :: otrm3
  real(dp)                                         :: otrm3_j
  real(dp)                                         :: otrm4
  real(dp)                                         :: otrm4_j
  real(dp)                                         :: ppen2
  real(dp)                                         :: pcoa1
  real(dp)                                         :: pcoa2
  real(dp)                                         :: pcoa3
  real(dp)                                         :: pcoa4
  real(dp)                                         :: pcot1
  real(dp)                                         :: phb1
  real(dp)                                         :: phb2
  real(dp)                                         :: phb3
  real(dp)                                         :: ptor1
  real(dp)                                         :: pval1
  real(dp)                                         :: pval2
  real(dp)                                         :: qij
  real(dp)                                         :: r0hb
  real(dp)                                         :: reaxFFrhtollower
  real(dp)                                         :: utrm1
  real(dp)                                         :: utrm1_j
  real(dp)                                         :: utrm2
  real(dp)                                         :: utrm2_j
  real(dp)                                         :: rintde
  real(dp)                                         :: drintde
  real(dp)                                         :: rij
  real(dp)                                         :: rij2
  real(dp)                                         :: ril
  real(dp)                                         :: ril2
  real(dp)                                         :: rkl
  real(dp)                                         :: rkl2
  real(dp)                                         :: rji
  real(dp)                                         :: rjk
  real(dp)                                         :: rjk2
  real(dp)                                         :: rki
  real(dp)                                         :: rki2
  real(dp)                                         :: rkj
  real(dp)                                         :: rkj2
  real(dp)                                         :: rrij
  real(dp)                                         :: rrik
  real(dp)                                         :: rrjk
  real(dp)                                         :: rnlp
  real(dp)                                         :: rtmp
  real(dp)                                         :: sbo
  real(dp)                                         :: sboprod
  real(dp)                                         :: scale
  real(dp)                                         :: sinkij
  real(dp)                                         :: sinijl
  real(dp)                                         :: sinhalftheta
  real(dp)                                         :: sin4halftheta
  real(dp)                                         :: smoothH
  real(dp)                                         :: sumBOcoa_ij
  real(dp)                                         :: sumBOcoa_ik
  real(dp)                                         :: sumBOcoa_jk
  real(dp)                                         :: sumdeltabopi
  real(dp)                                         :: sumdeltabopi_j
  real(dp)                                         :: sumover
  real(dp)                                         :: sumover_j
  real(dp)                                         :: sum_deoverunderddeltaj_neigh
  real(dp)                                         :: t1
  real(dp)                                         :: t2
  real(dp)                                         :: t3
  real(dp)                                         :: t4
  real(dp)                                         :: theta
  real(dp)                                         :: theta0
  real(dp)                                         :: tp
  real(dp)                                         :: dtpdr
  real(dp)                                         :: d2tpdr2
  real(dp)                                         :: tpQ
  real(dp)                                         :: dtpQdr
  real(dp)                                         :: d2tpQdr2
  real(dp)                                         :: trm1lp
  real(dp)                                         :: trm2lp
  real(dp)                                         :: V1
  real(dp)                                         :: V1trm
  real(dp)                                         :: V2
  real(dp)                                         :: V2trm
  real(dp)                                         :: V3
  real(dp)                                         :: V3trm
  real(dp)                                         :: VDWtrm
  real(dp),    dimension(:,:),   allocatable, save :: rneigh
  real(dp),    dimension(:,:),   allocatable, save :: xneigh
  real(dp),    dimension(:,:),   allocatable, save :: yneigh
  real(dp),    dimension(:,:),   allocatable, save :: zneigh
  real(dp),    dimension(:),     allocatable, save :: sum
  real(dp),    dimension(:),     allocatable, save :: sum0
  real(dp),    dimension(:),     allocatable, save :: sumsij
  real(dp)                                         :: xd
  real(dp)                                         :: yd
  real(dp)                                         :: zd
  real(dp)                                         :: xdiff
  real(dp)                                         :: ydiff
  real(dp)                                         :: zdiff
  real(dp)                                         :: xil
  real(dp)                                         :: yil
  real(dp)                                         :: zil
  real(dp)                                         :: xkl
  real(dp)                                         :: ykl
  real(dp)                                         :: zkl
  real(dp)                                         :: xji
  real(dp)                                         :: yji
  real(dp)                                         :: zji
  real(dp)                                         :: xjk
  real(dp)                                         :: yjk
  real(dp)                                         :: zjk
  real(dp)                                         :: xki
  real(dp)                                         :: yki
  real(dp)                                         :: zki
  real(dp)                                         :: xkj
  real(dp)                                         :: ykj
  real(dp)                                         :: zkj
  real(dp)                                         :: xshift
  real(dp)                                         :: yshift
  real(dp)                                         :: zshift
#ifdef TRACE
  call trace_in('reaxff1fc')
#endif
!
  t1 = g_cpu_time()
!
!  If lStrict is true then we try to set the behaviour to be the same as the
!  original ReaxFF code and turn off the modifications implemented for smoothness
!
  lStrict = (index(keyword,'stri').eq.1.or.index(keyword,' stri').ne.0)
!
!  Set lower bounds for HB taper to be 0.9 of upper bound for now
!
  reaxFFrhtollower = 0.9_dp*reaxFFrhtol
!
  allocate(nboatomptr(numat),stat=status)
  if (status/=0) call outofmemory('reaxFF1fc','nboatomptr')
  allocate(nboatomRptr(numat),stat=status)
  if (status/=0) call outofmemory('reaxFF1fc','nboatomRptr')
  allocate(nbodummyptr(numat),stat=status)
  if (status/=0) call outofmemory('reaxFF1fc','nbodummyptr')
!
!  Allocate memory that does not depend on maxneigh
!
  allocate(numattodoptr(numat),stat=status)
  if (status/=0) call outofmemory('reaxFF1fc','numattodoptr')
  allocate(numattodoRptr(numat),stat=status)
  if (status/=0) call outofmemory('reaxFF1fc','numattodoRptr')
  allocate(nbos(numat),stat=status)
  if (status/=0) call outofmemory('reaxFF1fc','nbos')
  allocate(nbosptr(numat),stat=status)
  if (status/=0) call outofmemory('reaxFF1fc','nbosptr')
  allocate(latomdone(numat),stat=status)
  if (status/=0) call outofmemory('reaxFF1fc','latomdone')
  allocate(nneigh(numat),stat=status)
  if (status/=0) call outofmemory('reaxFF1fc','nneigh')
  allocate(sumsij(numat),stat=status)
  if (status/=0) call outofmemory('reaxFF1fc','sumsij')
  allocate(sum(max(2*numat,14_i4)),stat=status)
  if (status/=0) call outofmemory('reaxFF1fc','sum')
  allocate(sum0(max(2*numat,14_i4)),stat=status)
  if (status/=0) call outofmemory('reaxFF1fc','sum0')
!
!  Reinitialisation point should maxneigh be increased
!
100 continue
  lmaxneighok = .true.
  if (allocated(rneigh)) then
    deallocate(d1BOj_pipi,stat=status)
    if (status/=0) call deallocate_error('reaxFF1fc','d1BOj_pipi')
    deallocate(d1BOj_pi,stat=status)
    if (status/=0) call deallocate_error('reaxFF1fc','d1BOj_pi')
    deallocate(d1BOj_s,stat=status)
    if (status/=0) call deallocate_error('reaxFF1fc','d1BOj_s')
    deallocate(d1BOj,stat=status)
    if (status/=0) call deallocate_error('reaxFF1fc','d1BOj')
    deallocate(d1BOi_pipi,stat=status)
    if (status/=0) call deallocate_error('reaxFF1fc','d1BOi_pipi')
    deallocate(d1BOi_pi,stat=status)
    if (status/=0) call deallocate_error('reaxFF1fc','d1BOi_pi')
    deallocate(d1BOi_s,stat=status)
    if (status/=0) call deallocate_error('reaxFF1fc','d1BOi_s')
    deallocate(d1BOi,stat=status)
    if (status/=0) call deallocate_error('reaxFF1fc','d1BOi')
    deallocate(d1j,stat=status)
    if (status/=0) call deallocate_error('reaxFF1fc','d1j')
    deallocate(d1i,stat=status)
    if (status/=0) call deallocate_error('reaxFF1fc','d1i')
    deallocate(devaldsbo_neigh,stat=status)
    if (status/=0) call deallocate_error('reaxFF1fc','devaldsbo_neigh')
    deallocate(detorsdBOpi_neigh,stat=status)
    if (status/=0) call deallocate_error('reaxFF1fc','detorsdBOpi_neigh')
    deallocate(detorsdBO_neigh,stat=status)
    if (status/=0) call deallocate_error('reaxFF1fc','detorsdBO_neigh')
    deallocate(devaldBO_neigh,stat=status)
    if (status/=0) call deallocate_error('reaxFF1fc','devaldBO_neigh')
    deallocate(debpendBO_neigh,stat=status)
    if (status/=0) call deallocate_error('reaxFF1fc','debpendBO_neigh')
    deallocate(deoverunderddeltaj_neigh,stat=status)
    if (status/=0) call deallocate_error('reaxFF1fc','deoverunderddeltaj_neigh')
    deallocate(d1BOp_pipi,stat=status)
    if (status/=0) call deallocate_error('reaxFF1fc','d1BOp_pipi')
    deallocate(d1BOp_pi,stat=status)
    if (status/=0) call deallocate_error('reaxFF1fc','d1BOp_pi')
    deallocate(d1BOp,stat=status)
    if (status/=0) call deallocate_error('reaxFF1fc','d1BOp')
    deallocate(deltalp,stat=status)
    if (status/=0) call deallocate_error('reaxFF1fc','deltalp')
    deallocate(deltap,stat=status)
    if (status/=0) call deallocate_error('reaxFF1fc','deltap')
    deallocate(delta,stat=status)
    if (status/=0) call deallocate_error('reaxFF1fc','delta')
    deallocate(BOp_pipi,stat=status)
    if (status/=0) call deallocate_error('reaxFF1fc','BOp_pipi')
    deallocate(BOp_pi,stat=status)
    if (status/=0) call deallocate_error('reaxFF1fc','BOp_pi')
    deallocate(BOp,stat=status)
    if (status/=0) call deallocate_error('reaxFF1fc','BOp')
    deallocate(BO_pipi,stat=status)
    if (status/=0) call deallocate_error('reaxFF1fc','BO_pipi')
    deallocate(BO_pi,stat=status)
    if (status/=0) call deallocate_error('reaxFF1fc','BO_pi')
    deallocate(BO,stat=status)
    if (status/=0) call deallocate_error('reaxFF1fc','BO')
    deallocate(zneigh,stat=status)
    if (status/=0) call deallocate_error('reaxFF1fc','zneigh')
    deallocate(yneigh,stat=status)
    if (status/=0) call deallocate_error('reaxFF1fc','yneigh')
    deallocate(xneigh,stat=status)
    if (status/=0) call deallocate_error('reaxFF1fc','xneigh')
    deallocate(rneigh,stat=status)
    if (status/=0) call deallocate_error('reaxFF1fc','rneigh')
    deallocate(ineigh,stat=status)
    if (status/=0) call deallocate_error('reaxFF1fc','ineigh')
    deallocate(neighnoRptr,stat=status)
    if (status/=0) call deallocate_error('reaxFF1fc','neighnoRptr')
    deallocate(neighno,stat=status)
    if (status/=0) call deallocate_error('reaxFF1fc','neighno')
  endif
!
!  Set parameter for pairwise storage memory
!
!  maxneigh2  is the size of first derivative arrays except d1i
!  maxneigh2i is the size of first derivative array d1i - this is larger to hold torsional derivatives
!
  maxneigh1   = maxneigh*(maxneigh + 1)/2
  maxneigh2   = maxneigh + maxneigh1
  maxneigh2i  = maxneigh + maxneigh1 + maxneigh*(maxneigh + 1)*maxneigh
!
!  Allocate local memory
!
  allocate(neighno(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('reaxFF1fc','neighno')
  allocate(neighnoRptr(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('reaxFF1fc','neighnoRptr')
  allocate(ineigh(3_i4,maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('reaxFF1fc','ineigh')
  allocate(rneigh(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('reaxFF1fc','rneigh')
  allocate(xneigh(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('reaxFF1fc','xneigh')
  allocate(yneigh(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('reaxFF1fc','yneigh')
  allocate(zneigh(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('reaxFF1fc','zneigh')
  allocate(BO(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('reaxFF1fc','BO')
  allocate(BO_pi(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('reaxFF1fc','BO_pi')
  allocate(BO_pipi(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('reaxFF1fc','BO_pipi')
  allocate(BOp(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('reaxFF1fc','BOp')
  allocate(BOp_pi(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('reaxFF1fc','BOp_pi')
  allocate(BOp_pipi(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('reaxFF1fc','BOp_pipi')
  allocate(delta(numat),stat=status)
  if (status/=0) call outofmemory('reaxFF1fc','delta')
  allocate(deltap(numat),stat=status)
  if (status/=0) call outofmemory('reaxFF1fc','deltap')
  allocate(deltalplp(numat),stat=status)
  if (status/=0) call outofmemory('reaxFF1fc','deltalplp')
  allocate(deltalp(numat),stat=status)
  if (status/=0) call outofmemory('reaxFF1fc','deltalp')
  allocate(d1BOp(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('reaxFF1fc','d1BOp')
  allocate(d1BOp_pi(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('reaxFF1fc','d1BOp_pi')
  allocate(d1BOp_pipi(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('reaxFF1fc','d1BOp_pipi')
  allocate(deoverunderddeltaj_neigh(maxneigh),stat=status)
  if (status/=0) call outofmemory('reaxFF1fc','deoverunderddeltaj_neigh')
  allocate(debpendBO_neigh(maxneigh),stat=status)
  if (status/=0) call outofmemory('reaxFF1fc','debpendBO_neigh')
  allocate(devaldBO_neigh(maxneigh),stat=status)
  if (status/=0) call outofmemory('reaxFF1fc','devaldBO_neigh')
  allocate(detorsdBO_neigh(maxneigh),stat=status)
  if (status/=0) call outofmemory('reaxFF1fc','detorsdBO_neigh')
  allocate(detorsdBOpi_neigh(maxneigh),stat=status)
  if (status/=0) call outofmemory('reaxFF1fc','detorsdBOpi_neigh')
  allocate(devaldsbo_neigh(maxneigh),stat=status)
  if (status/=0) call outofmemory('reaxFF1fc','devaldsbo_neigh')
  allocate(d1i(maxneigh2i),stat=status)
  if (status/=0) call outofmemory('reaxFF1fc','d1i')
  allocate(d1j(maxneigh2),stat=status)
  if (status/=0) call outofmemory('reaxFF1fc','d1j')
  allocate(d1BOi(maxneigh2),stat=status)
  if (status/=0) call outofmemory('reaxFF1fc','d1BOi')
  allocate(d1BOi_s(maxneigh2),stat=status)
  if (status/=0) call outofmemory('reaxFF1fc','d1BOi_s')
  allocate(d1BOi_pi(maxneigh2),stat=status)
  if (status/=0) call outofmemory('reaxFF1fc','d1BOi_pi')
  allocate(d1BOi_pipi(maxneigh2),stat=status)
  if (status/=0) call outofmemory('reaxFF1fc','d1BOi_pipi')
  allocate(d1BOj(maxneigh2),stat=status)
  if (status/=0) call outofmemory('reaxFF1fc','d1BOj')
  allocate(d1BOj_s(maxneigh2),stat=status)
  if (status/=0) call outofmemory('reaxFF1fc','d1BOj_s')
  allocate(d1BOj_pi(maxneigh2),stat=status)
  if (status/=0) call outofmemory('reaxFF1fc','d1BOj_pi')
  allocate(d1BOj_pipi(maxneigh2),stat=status)
  if (status/=0) call outofmemory('reaxFF1fc','d1BOj_pipi')
!*************************************
!  Find cut-off radii for all atoms  *
!*************************************
  nbosptr(1:numat) = 0
  nboatom = 0
  do i = 1,numat
    nati = nat(i)
    ntypi = nftype(i)
!
!  Set dummy pointer for derivative routines since neighbour list is not compressed yet for ReaxFF
!
    nbodummyptr(i) = i
!
!  Check repulsive terms and build atom pointers to species
!
    nbos(i) = 0
    do j = 1,nreaxFFspec
      if (nati.eq.natreaxFFspec(j).and.(ntypi.eq.ntypreaxFFspec(j).or.ntypreaxFFspec(j).eq.0)) then
        nbos(i) = nbos(i) + 1
        nbosptr(i) = j
      endif
    enddo
!
!  Check number of species for now
!
    if (nbos(i).gt.1) then
      call outerror('Multiple species per atom not yet allowed for in reaxFF',0_i4)
      call stopnow('reaxFF')
    elseif (nbos(i).eq.1) then
      nboatom = nboatom + 1
      nboatomptr(i) = nboatom
      nboatomRptr(nboatom) = i
    endif
  enddo
!****************************************
!  Calculate neighbour lists for atoms  *
!****************************************
  call getReaxFFneighbour(maxneigh,nbosptr,nneigh,neighno,rneigh, &
                          xneigh,yneigh,zneigh,ineigh,latomdone,lmaxneighok)
!**********************************************************************
!  If maxneigh has been exceeded, increase value and return to start  *
!**********************************************************************
  if (.not.lmaxneighok) then
    do i = 1,numat
      if (nneigh(i).gt.maxneigh) maxneigh = nneigh(i)
    enddo
    if (index(keyword,'verb').ne.0) then
      write(ioout,'(/,''  Increasing maxneigh to '',i6)') maxneigh
    endif
    goto 100
  endif
!
!  Shift coordinates of matom in neighbour list arrays and correct distances
!
  do i = 1,numat
    if (i.eq.matom) then
      do nn = 1,nneigh(i)
        if (neighno(nn,i).ne.matom) then
          xneigh(nn,i) = xneigh(nn,i) - (vmatom(1) - xclat(i))
          yneigh(nn,i) = yneigh(nn,i) - (vmatom(2) - yclat(i))
          zneigh(nn,i) = zneigh(nn,i) - (vmatom(3) - zclat(i))
          rtmp = xneigh(nn,i)**2 + yneigh(nn,i)**2 + zneigh(nn,i)**2
          rneigh(nn,i) = sqrt(rtmp)
        endif
      enddo
    else
      do nn = 1,nneigh(i)
        j = neighno(nn,i)
        if (j.eq.matom) then
          xneigh(nn,i) = xneigh(nn,i) + (vmatom(1) - xclat(j))
          yneigh(nn,i) = yneigh(nn,i) + (vmatom(2) - yclat(j))
          zneigh(nn,i) = zneigh(nn,i) + (vmatom(3) - zclat(j))
          rtmp = xneigh(nn,i)**2 + yneigh(nn,i)**2 + zneigh(nn,i)**2
          rneigh(nn,i) = sqrt(rtmp)
        endif
      enddo
    endif
  enddo
!
!  Set pointer to atoms that are needed 
!
!  The code below is the safe option that does all atoms
!
  numattodo = numat
  numattodo_main = numat
  do i = 1,numat
    numattodoptr(i) = i
    numattodoRptr(i) = i
  enddo
!*******************************
!  Sort neighbours into order  *
!******************************* 
  if (lspatialboOK) then
    do i = 1,numat
!               
!  Build pointer
!               
      do nn = 1,nneigh(i)
        nmin = numat + 1 
        do nn2 = nn,nneigh(i) 
          if (neighno(nn2,i).lt.nmin) then
            nmin = neighno(nn2,i)
            nptr = nn2  
          endif
        enddo       
!         
!  Sort quantities
!
        if (nptr.ne.nn) then
          itmp = neighno(nptr,i)
          neighno(nptr,i) = neighno(nn,i)
          neighno(nn,i)  = itmp
          itmp = ineigh(1,nptr,i)
          ineigh(1,nptr,i) = ineigh(1,nn,i)
          ineigh(1,nn,i)  = itmp
          itmp = ineigh(2,nptr,i)
          ineigh(2,nptr,i) = ineigh(2,nn,i)
          ineigh(2,nn,i)  = itmp
          itmp = ineigh(3,nptr,i)
          ineigh(3,nptr,i) = ineigh(3,nn,i)
          ineigh(3,nn,i)  = itmp
          rtmp = rneigh(nptr,i)
          rneigh(nptr,i) = rneigh(nn,i)
          rneigh(nn,i)  = rtmp
          rtmp = xneigh(nptr,i)
          xneigh(nptr,i) = xneigh(nn,i)
          xneigh(nn,i)  = rtmp
          rtmp = yneigh(nptr,i)
          yneigh(nptr,i) = yneigh(nn,i)
          yneigh(nn,i)  = rtmp
          rtmp = zneigh(nptr,i)
          zneigh(nptr,i) = zneigh(nn,i)
          zneigh(nn,i)  = rtmp
        endif  
      enddo         
    enddo
  endif
!*********************************************
!  Calculate numbers of neighbours for each  *
!*********************************************
  if (index(keyword,'debu').ne.0) then
    call mpbarrier
    if (ioproc) then
      write(ioout,'(/,''  Number of neighbours for atoms :'',/)')
      write(ioout,'(''  i'',2x,''Nat'',3x,''N'')')
    endif
    call mpbarrier
    do ii = 1,numattodo_main
      i = numattodoptr(ii)
      write(ioout,'(i4,1x,i3,1x,i6)') i,nat(i),nneigh(i)
    enddo
    call mpbarrier
    if (ioproc) then
      write(ioout,'(/,''  Neighbours of atoms :'',/)')
    endif
    call mpbarrier
    do ii = 1,numattodo_main
      i = numattodoptr(ii)
      write(ioout,'(i4,8(1x,i4))') i,(neighno(nn,i),nn=1,nneigh(i))
    enddo
    call mpbarrier
  endif
!*************************************************************************
!  Loop over pairs of atoms to compute bond order prime and delta prime  *
!*************************************************************************
!
!  mneigh contains the real value of maxneigh need after removing neighbours
!  whose bond order tolerance is below the allowed threshold
!
  mneigh = 0
  do ii = 1,numattodo
    i = numattodoptr(ii)
!
!  Set variables relating to i
!
    if (nbos(i).gt.0) then
      nspeci = nbosptr(i)
!
!  Initialise deltap
!
      deltap(i) = - reaxFFval(1,nspeci)
!
!  Loop over neighbours of i 
!
      nivalid = 0
      do ni = 1,nneigh(i)
        j = neighno(ni,i)
!
!  Does j have a bond order species type?
!
        if (nbos(j).gt.0) then
!
!  Set variables relating to j
!
          nspecj = nbosptr(j)
!
!  Set up i-j quantities
!
          rij = rneigh(ni,i)
!
!  Set index to parameters for pairwise interaction of species
!
          if (nspeci.ge.nspecj) then
            nboij = nspeci*(nspeci - 1)/2 + nspecj
          else
            nboij = nspecj*(nspecj - 1)/2 + nspeci
          endif
          lanyBOnonzero = .false.
          if (rij.lt.reaxFFrmaxpair(nboij)) then
!**********************************
!  Valid bond order contribution  *
!**********************************
!
!  Sigma
!
            call reaxFFbo_sigma(nboij,nspeci,nspecj,rij,BOpij,dBOpijdr,.true.,.not.lStrict,lBOnonzero)
            if (lBOnonzero) then
              lanyBOnonzero = .true.
              BOp(ni,i) = BOpij
              deltap(i) = deltap(i) + BOpij
              rrij = 1.0_dp/rij
              d1BOp(ni,i) = rrij*dBOpijdr
!
!  Pi
!
              call reaxFFbo_pi(nboij,nspeci,nspecj,rij,BOpij,dBOpijdr,.true.,.not.lStrict)
              BOp_pi(ni,i) = BOpij
              deltap(i) = deltap(i) + BOpij
              d1BOp_pi(ni,i) = rrij*dBOpijdr
!
!  Pi_pi
!
              call reaxFFbo_pipi(nboij,nspeci,nspecj,rij,BOpij,dBOpijdr,.true.,.not.lStrict)
              BOp_pipi(ni,i) = BOpij
              deltap(i) = deltap(i) + BOpij
              d1BOp_pipi(ni,i) = rrij*dBOpijdr
            else
              BOp(ni,i) = 0.0_dp
              BOp_pi(ni,i) = 0.0_dp
              BOp_pipi(ni,i) = 0.0_dp
              d1BOp(ni,i) = 0.0_dp
              d1BOp_pi(ni,i) = 0.0_dp
              d1BOp_pipi(ni,i) = 0.0_dp
            endif
          else
            BOp(ni,i) = 0.0_dp
            BOp_pi(ni,i) = 0.0_dp
            BOp_pipi(ni,i) = 0.0_dp
            d1BOp(ni,i) = 0.0_dp
            d1BOp_pi(ni,i) = 0.0_dp
            d1BOp_pipi(ni,i) = 0.0_dp
          endif
          if (lanyBOnonzero) then
!
!  Valid terms found - move terms to correct location
!
            nivalid = nivalid + 1
            if (nivalid.ne.ni) then
              neighno(nivalid,i) = neighno(ni,i)
              ineigh(1:3,nivalid,i) = ineigh(1:3,ni,i)
              rneigh(nivalid,i) = rneigh(ni,i)
              xneigh(nivalid,i) = xneigh(ni,i)
              yneigh(nivalid,i) = yneigh(ni,i)
              zneigh(nivalid,i) = zneigh(ni,i)
              BOp(nivalid,i) = BOp(ni,i)
              BOp_pi(nivalid,i) = BOp_pi(ni,i)
              BOp_pipi(nivalid,i) = BOp_pipi(ni,i)
              d1BOp(nivalid,i) = d1BOp(ni,i)
              d1BOp_pi(nivalid,i) = d1BOp_pi(ni,i)
              d1BOp_pipi(nivalid,i) = d1BOp_pipi(ni,i)
            endif
          endif
        endif
      enddo
!  
!  Now reset number of neighbours to reflect the number that have non-zero bond orders
!           
      nneigh(i) = nivalid
      mneigh = max(mneigh,nivalid)
    endif
!
!  End loop over atoms i
!
  enddo
!
!  Set neighnoRptr
!
  do ii = 1,numattodo
    i = numattodoptr(ii)
    if (nbos(i).gt.0) then
!
!  Loop over neighbours of i 
!
      nivalid = 0
      do ni = 1,nneigh(i)
        j = neighno(ni,i)
!
!  Does j have a bond order species type?
!
        if (nbos(j).gt.0) then
!
!  Set up i-j quantities
!
          rij = rneigh(ni,i)
          xji = xneigh(ni,i)
          yji = yneigh(ni,i)
          zji = zneigh(ni,i)
!
!  Find i in neighbour list for j
!
          nj = 1
          lfound = .false.
          do while (nj.le.nneigh(j).and..not.lfound)
            if (neighno(nj,j).eq.i) then
              xdiff = xneigh(nj,j) + xji
              ydiff = yneigh(nj,j) + yji
              zdiff = zneigh(nj,j) + zji
              lfound = ((abs(xdiff)+abs(ydiff)+abs(zdiff)).lt.1.0d-6)
            endif
            if (.not.lfound) nj = nj + 1
          enddo
          if (lfound) then
            neighnoRptr(ni,i) = nj
          else
            call outerror('neighbour lists are inconsistent in reaxFF',0_i4)
            call stopnow('reaxff')
          endif
        endif
      enddo
    endif
!
!  End loop over atoms i
!
  enddo
  if (index(keyword,'debu').ne.0) then
    call mpbarrier
    if (ioproc) then
      write(ioout,'(/,''  Number of neighbours for atoms after bond order screening:'',/)')
      write(ioout,'(''  i'',2x,''Nat'',3x,''N'')')
    endif
    call mpbarrier
    do ii = 1,numattodo_main
      i = numattodoptr(ii)
      write(ioout,'(i4,1x,i3,1x,i6)') i,nat(i),nneigh(i)
    enddo
    call mpbarrier
  endif
  if (index(keyword,'verb').ne.0) then
    call mpbarrier
    if (ioproc) then
      write(ioout,'(/,''  ReaxFF: Delta and Bond Order prime: '',/)')
    endif
    call mpbarrier
    do ii = 1,numattodo_main
      i = numattodoptr(ii)
      write(ioout,'(i4,f10.6,6(1x,i4,1x,f7.4))') i,deltap(i),(neighno(j,i),(BOp(j,i)+BOp_pi(j,i)+BOp_pipi(j,i)),j=1,nneigh(i))
    enddo
    call mpbarrier
    if (ioproc) then
      write(ioout,'(/,''  ReaxFF: Bond Order sigma prime: '',/)')
    endif
    call mpbarrier
    do ii = 1,numattodo_main
      i = numattodoptr(ii)
      write(ioout,'(i4,6(1x,i4,1x,f7.4))') i,(neighno(j,i),(BOp(j,i)),j=1,nneigh(i))
    enddo
    call mpbarrier
    if (ioproc) then
      write(ioout,'(/,''  ReaxFF: Bond Order pi prime: '',/)')
    endif
    call mpbarrier
    do ii = 1,numattodo_main
      i = numattodoptr(ii)
      write(ioout,'(i4,6(1x,i4,1x,f7.4))') i,(neighno(j,i),(BOp_pi(j,i)),j=1,nneigh(i))
    enddo
    call mpbarrier
    if (ioproc) then
      write(ioout,'(/,''  ReaxFF: Bond Order pi-pi prime: '',/)')
    endif
    call mpbarrier
    do ii = 1,numattodo_main
      i = numattodoptr(ii)
      write(ioout,'(i4,6(1x,i4,1x,f7.4))') i,(neighno(j,i),(BOp_pipi(j,i)),j=1,nneigh(i))
    enddo
    call mpbarrier
    if (ioproc) then
      write(ioout,'(/)')
    endif
    call mpbarrier
  endif
!*****************************
!  Loop over pairs of atoms  *
!*****************************
  iloop: do ii = 1,numattodo
    i = numattodoptr(ii)
!
!  If this atom has no ReaxFF terms then there is nothing to do...
!
    if (nbos(i).eq.0) cycle iloop
!
!  Set variables relating to i
!
    nspeci = nbosptr(i)
    nregioni = nregionno(nsft + nrelf2a(i))
    nregiontypi = nregiontype(nregioni,ncf)
!
!  Set total number of distances for neighbours of i
!
    nneighi1   = nneigh(i)*(nneigh(i) + 1)/2 
    nneighi2   = nneigh(i) + nneighi1
!
    nneighi2i  = nneigh(i) + nneighi1 + nneigh(i)*(mneigh+1)*mneigh
!
!  Initialise derivative storage for neighbours of i
!
    d1i(1:nneighi2) = 0.0_dp
!
!  Loop over neighbours of i (=> j)
!
    do ni = 1,nneigh(i)
      j = neighno(ni,i)
!
!  If j > i then exit loop
!
      if (j.gt.i) exit
!
      nregionj = nregionno(nsft + nrelf2a(j))
      nregiontypj = nregiontype(nregionj,ncf)
!
!  QM/MM handling : i & j are both QM atoms => exclude
!
      lQMMMok = .true.
      if (QMMMmode(ncf).gt.0) then
        if (nregiontypi.eq.1.and.nregiontypj.eq.1) lQMMMok = .false.
      endif
!
!  Do we need to do this pair of atoms
!
      if (nbos(j).gt.0) then
!
!  If i = j set scale to be half to correct for double counting
!
        if (i.eq.j) then
          scale = 0.5_dp
        else
          scale = 1.0_dp
        endif
!
!  Set variables relating to j
!
        nspecj = nbosptr(j)
!
!  Set up i-j quantities
!
        lreg2one  = .false.
        lreg2pair = .false.
!
!  Find i in neighbour list for j
!
        nj = neighnoRptr(ni,i)
!
!  Set total number of distances for neighbours of j
!
        nneighj1   = nneigh(j)*(nneigh(j) + 1)/2
        nneighj2   = nneigh(j) + nneighj1
!
!  Initialise derivative storage for neighbours of i
!
        d1j(1:nneighj2) = 0.0_dp
!
!  Set index to parameters for pairwise interaction of species
!
        if (nspeci.ge.nspecj) then
          nboij = nspeci*(nspeci - 1)/2 + nspecj
        else
          nboij = nspecj*(nspecj - 1)/2 + nspeci
        endif
!
!  Set delta' 
!
        deltapi = deltap(i)
        deltapj = deltap(j)
!***************************************
!  Valid pairwise reaxFF contribution  *
!***************************************
        call reaxFF_bo(i,j,nspeci,nspecj,ni,nj,deltapi,deltapj,nneigh,BOp,BOp_pi,BOp_pipi, &
                       d1BOp,d1BOp_pi,d1BOp_pipi,BOij_s,BOij_pi,BOij_pipi, &
                       d1BOi_s,d1BOi_pi,d1BOi_pipi,d1BOj_s,d1BOj_pi,d1BOj_pipi,.true.)
!
!  Save corrected bond orders
!
        BO(ni,i) = BOij_s + BOij_pi + BOij_pipi
        BO(nj,j) = BOij_s + BOij_pi + BOij_pipi
        BO_pi(ni,i) = BOij_pi
        BO_pi(nj,j) = BOij_pi
        BO_pipi(ni,i) = BOij_pipi
        BO_pipi(nj,j) = BOij_pipi
!
!  If lQMMMok then compute energy terms
!
        if (lQMMMok) then
!
!  Raise BOij to the power of Pbe,2 (in paper says Pbe,1) but I think this is a typo!
!
          BOijpbe2 = BOij_s**(reaxFFpbe(2,nboij))
!
!  Compute exponential factor
!
          expij = exp(reaxFFpbe(1,nboij)*(1.0_dp - BOijpbe2))
!
!  Derivatives of Bond Order potential energy
!
          deijdBO_s    = - scale*reaxFFDe(1,nboij)*expij*(1.0_dp - reaxFFpbe(1,nboij)*reaxFFpbe(2,nboij)*BOijpbe2)
          deijdBO_pi   = - scale*reaxFFDe(2,nboij) 
          deijdBO_pipi = - scale*reaxFFDe(3,nboij) 
          do k = 1,nneigh(i)
            d1i(k) = d1i(k) + deijdBO_s*d1BOi_s(k) + deijdBO_pi*d1BOi_pi(k) + deijdBO_pipi*d1BOi_pipi(k)
          enddo
          do k = 1,nneigh(j)
            d1j(k) = d1j(k) + deijdBO_s*d1BOj_s(k) + deijdBO_pi*d1BOj_pi(k) + deijdBO_pipi*d1BOj_pipi(k)
          enddo
!
!  Add derivatives due to neighbours of j 
!
          call d1addfc(j,maxneigh,mneigh,nneigh,neighno,xneigh,yneigh,zneigh,ineigh,maxlhs, &
                       d1cell,matom,nbodummyptr,d1j,.false.)
        endif
!
!  End condition section on i or j being associated with moving atom
!
      endif
    enddo
!
!  Add derivatives due to neighbours of i
!
    if (ii.le.numattodo_main) then
      call d1addfc(i,maxneigh,mneigh,nneigh,neighno,xneigh,yneigh,zneigh,ineigh,maxlhs, &
                   d1cell,matom,nbodummyptr,d1i,.false.)
    endif
  enddo iloop
!***************************************************
!  Loop over neighbours of atoms to compute delta  *
!***************************************************
  delta(1:numat) = 0.0_dp
  deltalp(1:numat) = 0.0_dp
  deltalplp(1:numat) = 0.0_dp
  do ii = 1,numattodo_main
    i = numattodoptr(ii)
!
!  Set variables relating to i
!
    if (nbos(i).gt.0) then
      nspeci = nbosptr(i)
!
!  Compute delta from refined bond orders
!
      delta(i) = - reaxFFval(1,nspeci)
      do ni = 1,nneigh(i)
        delta(i) = delta(i) + BO(ni,i)
      enddo
!
!  Compute 1/2 x delta_e
!
      delta_e = 0.5_dp*(delta(i) + reaxFFval(1,nspeci) - reaxFFval(3,nspeci))
!
!  Compute deltalp for lone pairs
!
      if (lreaxFFlpsmooth) then
!
!  New smooth form
!
        delta_e0 = nint(delta_e)
        smoothH = 1.0_dp/(1.0_dp + exp(-reaxFFlpsmooth*(delta_e - delta_e0)))
        rintde = delta_e0 - 1.0_dp + smoothH
      else
!
!  Old discontinuous form
!
        rintde = dble(int(delta_e))
      endif
      rnlp = - rintde + exp(-4.0_dp*reaxFFlam(29)*(1.0_dp + delta_e - rintde)**2)
      deltalplp(i) = reaxFFlp(1,nspeci) - rnlp
!
!  Compute deltalp for lone pairs for overcoordination term
!
      if (lreaxFFovsmooth) then
!
!  New smooth form
!
        delta_e0 = nint(delta_e)
        smoothH = 1.0_dp/(1.0_dp + exp(-reaxFFlpsmooth*(delta_e - delta_e0)))
        rintde = delta_e0 - 1.0_dp + smoothH
      else
!
!  Old discontinuous form
!
        rintde = dble(int(delta_e))
      endif
      rnlp = - rintde + exp(-4.0_dp*reaxFFlam(29)*(1.0_dp + delta_e - rintde)**2)
      deltalp(i) = reaxFFlp(1,nspeci) - rnlp
    endif
  enddo
  if (index(keyword,'verb').ne.0) then
    call mpbarrier
    if (ioproc) then
      write(ioout,'(/,''  ReaxFF: Delta and Bond Order: '',/)')
    endif
    call mpbarrier
    do ii = 1,numattodo_main
      i = numattodoptr(ii)
      write(ioout,'(i4,f10.6,6(1x,i4,1x,f7.4))') i,delta(i),(neighno(j,i),BO(j,i),j=1,nneigh(i))
    enddo
    call mpbarrier
    if (ioproc) then
      write(ioout,'(/,''  ReaxFF: Delta lone pair: '',/)')
    endif
    call mpbarrier
    do ii = 1,numattodo_main
      i = numattodoptr(ii)
      write(ioout,'(i4,f10.6)') i,deltalplp(i)
    enddo
    call mpbarrier
    if (ioproc) then
      write(ioout,'(/,''  ReaxFF: Delta lone pair for overcoordination: '',/)')
      do i = 1,numat
        write(ioout,'(i4,f10.6)') i,deltalp(i)
      enddo
      write(ioout,'(/)')
    endif
    call mpbarrier
    if (ioproc) then
      write(ioout,'(/)')
    endif
    call mpbarrier
  endif
!***********************************************************
!  Loop over atoms to compute over and under coordination  *
!***********************************************************
  do ii = 1,numattodo_main
    i = numattodoptr(ii)
!
!  Set variables relating to i
!
    if (nbos(i).gt.0) then
      nspeci = nbosptr(i)
      nregioni = nregionno(nsft + nrelf2a(i))
      nregiontypi = nregiontype(nregioni,ncf)
!   
!  QM/MM handling : i is a QM atom => exclude
!   
      lQMMMok = .true.                           
      if (QMMMmode(ncf).gt.0) then               
        if (nregiontypi.eq.1) lQMMMok = .false.  
      endif
!
!  Do we need to do this atom based on QMMM scheme?
!
      if (lQMMMok) then
!
!  Set total number of distances for neighbours of i
!
        nneighi1   = nneigh(i)*(nneigh(i) + 1)/2 
        nneighi2   = nneigh(i) + nneighi1
        nneighi2i  = nneigh(i) + nneighi1 + nneigh(i)*(mneigh+1)*mneigh
!
!  Initialise derivative storage for neighbours of i
!
        d1i(1:nneighi2i) = 0.0_dp
!
!  Lone pair energy
!
        if (abs(reaxFFlp(2,nspeci)).gt.1.0d-12) then
          explp = exp(-reaxFFlam(30)*deltalp(i))
          trm1lp = 1.0_dp/(1.0_dp + explp)
        endif
!
!  Bond penalty energy
!
        do ni = 1,nneigh(i)
          j = neighno(ni,i)
          if (nbos(j).gt.0) then
            nspecj = nbosptr(j)
            if (nspeci.ge.nspecj) then
              nboij = nspeci*(nspeci - 1)/2 + nspecj
            else
              nboij = nspecj*(nspecj - 1)/2 + nspeci
            endif
            diffBOdelta = BO(ni,i) - delta(i) - reaxFFpen2(2,nboij)*delta(i)**4 - reaxFFpen2(3,nboij)
            if (diffBOdelta.gt.0.0_dp) then
              ebpenij = reaxFFpen2(1,nboij)*diffBOdelta**2
              ebpen = ebpen + ebpenij
            endif
          endif
        enddo
!
!  Compute sum over neighbours for over/under-coordination energy terms
!
        sumdeltabopi = 0.0_dp
        sumover = 0.0_dp
        do ni = 1,nneigh(i)
          j = neighno(ni,i)
          if (nbos(j).gt.0) then
            nspecj = nbosptr(j)
            if (nspeci.ge.nspecj) then
              nboij = nspeci*(nspeci - 1)/2 + nspecj
            else
              nboij = nspecj*(nspecj - 1)/2 + nspeci
            endif
            if (.not.lreaxFFunder(nspeci).or..not.lreaxFFunder(nspecj)) then
              sumdeltabopi = sumdeltabopi + delta(j)*(BO_pi(ni,i) + BO_pipi(ni,i))
            else
              sumdeltabopi = sumdeltabopi + (delta(j) - deltalp(j))*(BO_pi(ni,i) + BO_pipi(ni,i))
            endif
            sumover = sumover + reaxFFoc2(nboij)*reaxFFDe(1,nboij)*BO(ni,i)
          endif
        enddo
!
!  Compute delta_lpcorr_i - expression is different for first row elements
!
        if (.not.lreaxFFunder(nspeci)) then
          expdi = 0.0_dp
          otrm1 = 0.0_dp
          deltalp_corr_i = delta(i)
        else
          expdi = exp(reaxFFlam(31)*sumdeltabopi)
          otrm1 = 1.0_dp/(1.0_dp + reaxFFlam(6)*expdi)
          deltalp_corr_i = delta(i) - deltalp(i)*otrm1
        endif
!
!  Over coordination term
!
        otrm2 = 1.0_dp/(deltalp_corr_i + reaxFFval(1,nspeci) + 1.0d-8)
        expoc = exp(reaxFFoc1(nspeci)*deltalp_corr_i)
        otrm3 = 1.0_dp/(1.0_dp + expoc)
!
!  Under coordination term
!
        expuc1 = exp(reaxFFlam(7)*deltalp_corr_i)
        expuc2 = exp(-reaxFFoc1(nspeci)*deltalp_corr_i)
        utrm1  = 1.0_dp/(1.0_dp + expuc2)
        expuc3 = exp(reaxFFlam(9)*sumdeltabopi)
        utrm2  = 1.0_dp/(1.0_dp + reaxFFlam(8)*expuc3)
!
!  Compute SBO term needed for valence angle energy
!
!  During pass over neighbours, sbo stores sum of pi bond orders & sboprod is product of exponentiated bond orders
!
        sbo = 0.0_dp
        sboprod = 1.0_dp
        do ni = 1,nneigh(i)
          sbo = sbo + BO_pi(ni,i) + BO_pipi(ni,i)
          bo8 = BO(ni,i)**8
          expbo8 = exp(-bo8)
          sboprod = sboprod*expbo8
        enddo
!
!  Combine sbo contributions
!
        delta_angle = delta(i) + reaxFFval(1,nspeci) - reaxFFval(4,nspeci)
        rnlp = reaxFFlp(1,nspeci) - deltalp(i)
        sbo = sbo - (1.0_dp - sboprod)*(delta_angle + reaxFFlam(16)*rnlp)
!
!  Valence angle and penalty term
!
        devaldsbo = 0.0_dp
        devalddeltai = 0.0_dp
        devaldBO_neigh(1:nneigh(i)) = 0.0_dp
        devaldsbo_neigh(1:nneigh(i)) = 0.0_dp
!
!  Set up derivative of delta_lp with respect to delta_e
!
        delta_e = 0.5_dp*(delta(i) + reaxFFval(1,nspeci) - reaxFFval(3,nspeci))
!
        if (lreaxFFlpsmooth) then
!
!  New smooth form
!
          delta_e0 = nint(delta_e)
          expH = exp(-reaxFFlpsmooth*(delta_e - delta_e0))
          smoothH = 1.0_dp/(1.0_dp + expH)
          rintde = delta_e0 - 1.0_dp + smoothH
          drintde = 0.5_dp*reaxFFlpsmooth*expH*smoothH**2
        else
!
!  Old discontinuous form
!
          rintde = dble(int(delta_e))
          drintde = 0.0_dp
        endif
!
        trm2lp = 2.0_dp*(1.0_dp + delta_e - rintde)
        explp2 = exp(-reaxFFlam(29)*trm2lp**2)
        ddeltalpddeltaelp = drintde + 4.0_dp*reaxFFlam(29)*trm2lp*explp2*(0.5_dp - drintde)
!
        if (lreaxFFovsmooth) then
!
!  New smooth form
!
          delta_e0 = nint(delta_e)
          expH = exp(-reaxFFlpsmooth*(delta_e - delta_e0))
          smoothH = 1.0_dp/(1.0_dp + expH)
          rintde = delta_e0 - 1.0_dp + smoothH
          drintde = 0.5_dp*reaxFFlpsmooth*expH*smoothH**2
        else
!
!  Old discontinuous form
!
          rintde = dble(int(delta_e))
          drintde = 0.0_dp
        endif
!
        trm2lp = 2.0_dp*(1.0_dp + delta_e - rintde)
        explp2 = exp(-reaxFFlam(29)*trm2lp**2)
        ddeltalpddeltae = drintde + 4.0_dp*reaxFFlam(29)*trm2lp*explp2*(0.5_dp - drintde)
!
!  Compute delta for 3-body conjugation
!
        delta_coa = delta(i) + reaxFFval(1,nspeci) - reaxFFval(2,nspeci)
!******************
!  Angular terms  *
!******************
!
!  NB: The use of cutoffs is not applied to the hydrogen bond term, as per the original program
!
        do ni = 2,nneigh(i)
          j = neighno(ni,i)
          if (nbos(j).gt.0) then
            nspecj = nbosptr(j)
            rij = rneigh(ni,i)
!
!  Bond order check for i-j
!
            BOij = BO(ni,i) - reaxFFatol
            lBOijOK = (BOij.gt.0.0_dp) 
!
!  Calculate taper function to smooth cutoff for BOij
!
            if (lStrict) then
              fij = 1.0_dp
              dfijdBOij = 0.0_dp
            else
              call ataper(.false.,BO(ni,i),reaxFFatol,reaxFFtaperscale*reaxFFatol,fij,dfijdBOij,d2fijdBOij2,d3fijdBOij3, &
                          .true.,.false.,.false.)
            endif
            do nk = 1,ni-1
              k = neighno(nk,i)
              if (nbos(k).gt.0) then
                nspeck = nbosptr(k)
!
!  Bond order check for i-k
!
                BOik = BO(nk,i) - reaxFFatol
                lBOikOK = (BOik.gt.0.0_dp)
!
!  Calculate taper function to smooth cutoff for BOij
!
                if (lStrict) then
                  fik = 1.0_dp
                  dfikdBOik = 0.0_dp
                else
                  call ataper(.false.,BO(nk,i),reaxFFatol,reaxFFtaperscale*reaxFFatol,fik,dfikdBOik,d2fikdBOik2,d3fikdBOik3, &
                              .true.,.false.,.false.)
                endif
!
!  Find triad index
!
                if (nspecj.ge.nspeck) then
                  ind = nspecj*(nspecj - 1)/2 + nspeck
                else
                  ind = nspeck*(nspeck - 1)/2 + nspecj
                endif
!
!  Compute distance squared between j - k
!
                xjk = xneigh(ni,i) - xneigh(nk,i)
                yjk = yneigh(ni,i) - yneigh(nk,i)
                zjk = zneigh(ni,i) - zneigh(nk,i)
                rjk2 = xjk**2 + yjk**2 + zjk**2
!
!  Compute theta for j - i - k
!
                rjk = sqrt(rjk2)
                call reaxff_theta(rneigh(ni,i),rneigh(nk,i),rjk,theta,dthetadrij,dthetadrik,dthetadrjk, &
                                  d2thetadr2,.true.,.false.)
                if (lBOijOK.and.lBOikOK.and.(BO(ni,i)*BO(nk,i).gt.reaxFFatol2)) then
!------------------
!  Valence energy |
!------------------
!
!  Loop over three-body valence potentials
!
                  do nval3 = 1,nreaxFFval3(ind,nspeci)
                    pval1 = reaxFFval3(2,nval3,ind,nspeci)
                    pval2 = reaxFFval3(3,nval3,ind,nspeci)
!
!  Compute f7, f8 and theta0 for this triad
!
                    call reaxFF_f7(nspeci,nspecj,nspeck,nval3,BOij,BOik,f7,df7dBOij,df7dBOik, &
                                   d2f7dBOij2,d2f7dBOijdBOik,d2f7dBOik2,.true.,.false.)
                    call reaxFF_f8(nspeci,nspecj,nspeck,nval3,delta(i),f8,df8ddeltai,d2f8ddeltai2,.true.,.false.)
                    call reaxFF_theta0(nspeci,nspecj,nspeck,nval3,sbo,theta0,dtheta0dsbo,d2theta0dsbo2,.true.,.false.)
!
!  Compute theta dependent terms
!
                    exptheta = exp(-pval2*(theta0 - theta)**2)
!
!  Derivatives of valence energy
!
!  
!  Compute derivative terms from f7 and f8
!                 
                    devalddeltai = devalddeltai + pval1*fij*fik*f7*df8ddeltai*(1.0_dp - exptheta)
                    devaldBO_neigh(ni) = devaldBO_neigh(ni) + pval1*fij*fik*f8*df7dBOij*(1.0_dp - exptheta)
                    devaldBO_neigh(ni) = devaldBO_neigh(ni) + pval1*dfijdBOij*fik*f8*f7*(1.0_dp - exptheta)
                    devaldBO_neigh(nk) = devaldBO_neigh(nk) + pval1*fij*fik*f8*df7dBOik*(1.0_dp - exptheta)
                    devaldBO_neigh(nk) = devaldBO_neigh(nk) + pval1*fij*dfikdBOik*f8*f7*(1.0_dp - exptheta)
!                 
!  Compute derivative terms from (1 - exp(-pval2*(theta0 - theta)**2) - Part 1 derivatives for theta
!                 
                    devaldtheta = pval1*fij*fik*f7*f8*(2.0_dp*pval2*exptheta*(theta0 - theta))
                    d1i(ni) = d1i(ni) - devaldtheta*dthetadrij
                    d1i(nk) = d1i(nk) - devaldtheta*dthetadrik
                    indjk = nneigh(i) + ni*(ni - 1)/2 + nk
                    d1i(indjk) = d1i(indjk) - devaldtheta*dthetadrjk
!
!  Compute derivative terms from (1 - exp(-pval2*(theta0 - theta)**2) - Part 2 derivatives for theta0
!               
                    devaldsbo = devaldsbo + devaldtheta*dtheta0dsbo
                    dsboddeltai = - (1.0_dp - sboprod)
                    devalddeltai = devalddeltai + devaldtheta*dtheta0dsbo*dsboddeltai*(1.0_dp - reaxFFlam(16)*ddeltalpddeltae)
                    do nl = 1,nneigh(i)
                      dsboproddbo = - (8.0_dp*BO(nl,i)**7)*sboprod*(delta_angle + reaxFFlam(16)*rnlp)
                      devaldsbo_neigh(nl) = devaldsbo_neigh(nl) + devaldtheta*dtheta0dsbo*dsboproddbo
                    enddo
!
!  End of loop over three-body valence potentials
!
                  enddo
!------------------
!  Penalty energy |
!------------------
                  lDoEpen = (abs(reaxFFpen3(ind,nspeci)).gt.1.0d-12)
                  if (lDoEpen) then
!
!  Compute f9 for this triad
!
                    call reaxFF_f9(delta(i),f9,df9ddeltai,d2f9ddeltai2,.true.,.false.)
!
!  Evaluate penalty energy
!
                    ppen2 = reaxFFlam(20)
                    exppen1 = exp(-ppen2*(BOij - 2.0_dp)**2)
                    exppen2 = exp(-ppen2*(BOik - 2.0_dp)**2)
                    epentrm = reaxFFpen3(ind,nspeci)*f9*exppen1*exppen2
!  
!  Derivatives of penalty energy with respect to deltai & bond orders -> add to eval terms
!
                    devalddeltai = devalddeltai + fij*fik*reaxFFpen3(ind,nspeci)*df9ddeltai*exppen1*exppen2
                    devaldBO_neigh(ni) = devaldBO_neigh(ni) - 2.0_dp*fij*fik*ppen2*epentrm*(BOij - 2.0_dp)
                    devaldBO_neigh(ni) = devaldBO_neigh(ni) + dfijdBOij*fik*epentrm
                    devaldBO_neigh(nk) = devaldBO_neigh(nk) - 2.0_dp*fij*fik*ppen2*epentrm*(BOik - 2.0_dp)
                    devaldBO_neigh(nk) = devaldBO_neigh(nk) + dfikdBOik*fij*epentrm
                  endif
!-----------------------------
!  3-body conjugation energy |
!-----------------------------
!
!  Set parameters for ecoa
!
                  pcoa1 = reaxFFconj3(1,ind,nspeci)
!
!  If pcoa1 is zero then there is no potential to do
!
                  if (abs(pcoa1).gt.1.0d-12) then
                    pcoa2 = reaxFFconj3(2,ind,nspeci)
                    pcoa3 = reaxFFconj3(3,ind,nspeci)
                    pcoa4 = reaxFFconj3(4,ind,nspeci)
!
!  Compute sum of bond orders for terminal atoms excluding i-j / i-k contribution
!  => use the fact that delta is already the sum of all bond orders
!
                    sumBOcoa_ij = delta(j) + reaxFFval(1,nspecj) - BOij
                    sumBOcoa_ik = delta(k) + reaxFFval(1,nspeck) - BOik
!
!  Evaluate 3-body conjugation energy
!
                    expcoa2 = exp(pcoa2*delta_coa)
                    expcoatrm = 1.0_dp/(1.0_dp + expcoa2)
                    expcoa3a = exp(-pcoa3*sumBOcoa_ij**2)
                    expcoa3b = exp(-pcoa3*sumBOcoa_ik**2)
                    expcoa4a = exp(-pcoa4*(BOij - 1.5_dp)**2)
                    expcoa4b = exp(-pcoa4*(BOik - 1.5_dp)**2)
                    expcoaprod = pcoa1*expcoa3a*expcoa3b*expcoa4a*expcoa4b
                    ecoatrm = expcoaprod*expcoatrm
!
!  Derivatives of 3-body conjugation energy with respect to deltai & bond orders -> add to eval terms
!
                    devalddeltai = devalddeltai - fij*fik*pcoa2*expcoatrm*ecoatrm*expcoa2
                    devaldBO_neigh(ni) = devaldBO_neigh(ni) - 2.0_dp*fij*fik*pcoa4*(BOij - 1.5_dp)*ecoatrm
                    devaldBO_neigh(nk) = devaldBO_neigh(nk) - 2.0_dp*fij*fik*pcoa4*(BOik - 1.5_dp)*ecoatrm
                    if (.not.lStrict) then
                      devaldBO_neigh(ni) = devaldBO_neigh(ni) + dfijdBOij*fik*ecoatrm
                      devaldBO_neigh(nk) = devaldBO_neigh(nk) + dfikdBOik*fij*ecoatrm
                    endif
                  endif
                endif
!
              endif
            enddo
!
          endif
        enddo
!***************************
!  Angular terms - part 2  *
!***************************
!
!  If doing gradients then we need to search for cases where i is a terminal atom
!  in order to compute the derivatives of the terminal atom bond orders for the
!  three-body conjugation energy.
!
        do ni = 1,nneigh(i)
          j = neighno(ni,i)
          if (nbos(j).gt.0) then
            nspecj = nbosptr(j)
!
!  Bond order check for i-j
!
            BOij = BO(ni,i) - reaxFFatol
            lBOijOK = (BOij.gt.0.0_dp) 
            if (lBOijOK) then
!
!  Calculate taper function to smooth cutoff for BOij
!
              if (lStrict) then
                fij = 1.0_dp
                dfijdBOij = 0.0_dp
              else
                call ataper(.false.,BO(ni,i),reaxFFatol,reaxFFtaperscale*reaxFFatol,fij,dfijdBOij,d2fijdBOij2,d3fijdBOij3, &
                            .true.,.false.,.false.)
              endif
!
!  Compute delta for j for 3-body conjugation
!
              delta_coa = delta(j) + reaxFFval(1,nspecj) - reaxFFval(2,nspecj)
!
!  Compute terms that depend on i-j only
!
              sumBOcoa_ij = delta(i) + reaxFFval(1,nspeci) - BOij
            endif
!
            do nk = 1,nneigh(j)
              k = neighno(nk,j)
              if (nbos(k).gt.0.and.k.ne.i) then
                nspeck = nbosptr(k)
!
!  Bond order check for j-k
!
                BOjk = BO(nk,j) - reaxFFatol
                lBOjkOK = (BOjk.gt.0.0_dp) 
!
!  Calculate taper function to smooth cutoff for BOij
!
                if (lStrict) then
                  fjk = 1.0_dp
                  dfjkdBOjk = 0.0_dp
                else
                  call ataper(.false.,BO(nk,j),reaxFFatol,reaxFFtaperscale*reaxFFatol,fjk,dfjkdBOjk,d2fjkdBOjk2,d3fjkdBOjk3, &
                              .true.,.false.,.false.)
                endif
!
!  Find triad indices for j being pivot atom
!
                if (nspeci.ge.nspeck) then
                  ind = nspeci*(nspeci - 1)/2 + nspeck
                else
                  ind = nspeck*(nspeck - 1)/2 + nspeci
                endif
                if (lBOijOK.and.lBOjkOK.and.(BO(ni,i)*BO(nk,j).gt.reaxFFatol2)) then
!
!  Set parameters
!
                  pcoa1 = reaxFFconj3(1,ind,nspecj)
!
!  If pcoa1 is zero then there is no potential to do
!
                  if (abs(pcoa1).gt.1.0d-12) then
                    pcoa2 = reaxFFconj3(2,ind,nspecj)
                    pcoa3 = reaxFFconj3(3,ind,nspecj)
                    pcoa4 = reaxFFconj3(4,ind,nspecj)
!
!  Compute sum of bond orders for terminal atoms excluding i-j / j-k contribution
!  => use the fact that delta is already the sum of all bond orders
!
                    sumBOcoa_jk = delta(k) + reaxFFval(1,nspeck) - BOjk
!
!  Evaluate 3-body conjugation terms
!
                    expcoa2 = exp(pcoa2*delta_coa)
                    expcoatrm = 1.0_dp/(1.0_dp + expcoa2)
                    expcoa3a = exp(-pcoa3*sumBOcoa_ij**2)
                    expcoa3b = exp(-pcoa3*sumBOcoa_jk**2)
                    expcoa4a = exp(-pcoa4*(BOij - 1.5_dp)**2)
                    expcoa4b = exp(-pcoa4*(BOjk - 1.5_dp)**2)
                    expcoaprod = pcoa1*expcoa3a*expcoa3b*expcoa4a*expcoa4b
                    ecoatrm = - 2.0_dp*pcoa3*sumBOcoa_ij*expcoaprod*expcoatrm
!
!  Compute derivatives of 3-body conjugation energy with respect to bond orders of i -> add to eval terms
!
                    do m = 1,nneigh(i)
                      if (m.ne.ni) then
                        devaldBO_neigh(m) = devaldBO_neigh(m) + fij*fjk*ecoatrm
                      endif
                    enddo
                  endif
!
!  End bond order check for j-k
!
                endif
              endif
            enddo
!
          endif
        enddo
!********************
!  Torsional terms  *
!********************
!
!  Initialise arrays for summing the derivatives of torsional energy with respect to bond orders of neighbours of i
!
        detorsdBO_neigh(1:nneigh(i)) = 0.0_dp
        detorsdBOpi_neigh(1:nneigh(i)) = 0.0_dp
        detorsddeltai = 0.0_dp
!
!  Loop over neighbours of j to find central bonds 
!
        do ni = 1,nneigh(i)
          j = neighno(ni,i)
!
!  Impose conditions that j must be a reaxFF atom 
!
!  Ideally we would only do the case where j < i. However, in order to get the bond order derivatives
!  we have to do both cases. When j < i, the energy and derivatives due to all terms bar one are 
!  computed. When j > i, the only thing computed is the derivative of the torsional energy with
!  respect to the i-k bond order.
!
          if (nbos(j).gt.0) then
!
!  Check bond order for i-j
!
            BOij = BO(ni,i) - reaxFFatol
            if (BOij.gt.0.0_dp) then
              if (lStrict) then                    
                fij = 1.0_dp                       
                dfijdBOij = 0.0_dp                 
              else                                 
                call ataper(.false.,BO(ni,i),reaxFFatol,reaxFFtaperscale*reaxFFatol,fij,dfijdBOij,d2fijdBOij2,d3fijdBOij3, &
                            .true.,.false.,.false.)
              endif
              if (i.gt.j) then
                ltoriGTj = .true.
              else
                ltoriGTj = .false.
              endif
              nspecj = nbosptr(j)
              if (nspeci.ge.nspecj) then
                ind1 = nspeci*(nspeci - 1)/2 + nspecj
              else
                ind1 = nspecj*(nspecj - 1)/2 + nspeci
              endif
!
!  Loop over neighbours of i to find terminal atom k while excluding k = j
!
              do nk = 1,nneigh(i)
                k = neighno(nk,i)
!  The condition below modified from k.ne.i to handle small cell systems
                if (nbos(k).gt.0.and.nk.ne.ni) then
!
!  Check bond order for i-k
!
                  BOik = BO(nk,i) - reaxFFatol
                  if (BOik.gt.0.0_dp) then
                    if (lStrict) then
                      fik = 1.0_dp
                      dfikdBOik = 0.0_dp
                    else
                      call ataper(.false.,BO(nk,i),reaxFFatol,reaxFFtaperscale*reaxFFatol,fik,dfikdBOik,d2fikdBOik2, &
                                  d3fikdBOik3,.true.,.false.,.false.)
                    endif
                    nspeck = nbosptr(k)
!
!  Loop over neighbours of j to find terminal atom l while excluding i = l and k = l
!
                    do nl = 1,nneigh(j)
                      l = neighno(nl,j)
!
!  Check bond order for j-l
!
                      BOjl = BO(nl,j) - reaxFFatol
                      if (BOjl.gt.0.0_dp.and.BO(ni,i)*BO(nk,i)*BO(nl,j).gt.reaxFFatol3) then
                        if (lStrict) then
                          fjl = 1.0_dp
                          dfjldBOjl = 0.0_dp
                        else
                          call ataper(.false.,BO(nl,j),reaxFFatol,reaxFFtaperscale*reaxFFatol,fjl,dfjldBOjl,d2fjldBOjl2, &
                                      d3fjldBOjl3,.true.,.false.,.false.)
                        endif
!
                        xkl = xneigh(ni,i) + xneigh(nl,j) - xneigh(nk,i)
                        ykl = yneigh(ni,i) + yneigh(nl,j) - yneigh(nk,i)
                        zkl = zneigh(ni,i) + zneigh(nl,j) - zneigh(nk,i)
                        rkl2 = xkl*xkl + ykl*ykl + zkl*zkl
!  The condition below modified from l.ne.i & l.ne.k to handle small cell systems
                        if (nbos(l).gt.0.and.nl.ne.neighnoRptr(ni,i).and.rkl2.gt.1.0d-8) then
                          nspecl = nbosptr(l)
!
!  We now have a quartet of valid reaxFF species
!
                          if (nspeck.ge.nspecl) then
                            ind2 = nspeck*(nspeck - 1)/2 + nspecl + 1
                          else
                            ind2 = nspecl*(nspecl - 1)/2 + nspeck + 1
                          endif
!
!  Assign local scalars
!
                          V1    = 0.5_dp*reaxFFtor4(1,ind2,ind1)
                          V2    = 0.5_dp*reaxFFtor4(2,ind2,ind1)
                          V3    = 0.5_dp*reaxFFtor4(3,ind2,ind1)
                          ptor1 = reaxFFtor4(4,ind2,ind1)
                          pcot1 = reaxFFtor4(5,ind2,ind1)
                          if (V1.eq.0.0_dp.and.V2.eq.0.0_dp.and.V3.eq.0.0_dp) then
!
!  No specific term so try wildcard
!
                            V1 = 0.5_dp*reaxFFtor4(1,1,ind1)
                            V2 = 0.5_dp*reaxFFtor4(2,1,ind1)
                            V3 = 0.5_dp*reaxFFtor4(3,1,ind1)
                            ptor1 = reaxFFtor4(4,1,ind1)
                            pcot1 = reaxFFtor4(5,1,ind1)
                          endif
!
!  Check whether the potentials are non-zero to avoid wasted calculation
!
                          if (abs(V1)+abs(V2)+abs(V3)+abs(pcot1).gt.1.0d-12) then
!
!  Compute unknown distances, j-k and i-l
!
                            xjk =   xneigh(ni,i) - xneigh(nk,i)
                            yjk =   yneigh(ni,i) - yneigh(nk,i)
                            zjk =   zneigh(ni,i) - zneigh(nk,i)
                            rjk2 = xjk**2 + yjk**2 + zjk**2
                            xil = - xneigh(ni,i) - xneigh(nl,j)
                            yil = - yneigh(ni,i) - yneigh(nl,j)
                            zil = - zneigh(ni,i) - zneigh(nl,j)
                            ril2 = xil**2 + yil**2 + zil**2
                            xkl = xil + xneigh(nk,i) 
                            ykl = yil + yneigh(nk,i) 
                            zkl = zil + zneigh(nk,i) 
                            rkl2 = xkl**2 + ykl**2 + zkl**2
!
!  Calculate angles k-i-j and i-j-l
!
                            rjk = sqrt(rjk2)
                            ril = sqrt(ril2)
                            rkl = sqrt(rkl2)
                            call reaxff_sintheta(rneigh(ni,i),rneigh(nk,i),rjk,sinkij,dsinkijdrij,dsinkijdrik,dsinkijdrjk, &
                                                 d2sinkijdr2,.true.,.false.)
                            call reaxff_sintheta(rneigh(ni,i),rneigh(nl,j),ril,sinijl,dsinijldrij,dsinijldrjl,dsinijldril, &
                                                 d2sinijldr2,.true.,.false.)
!
!  Compute torsion angle
!
                            call reaxff_cosphi(rneigh(nk,i),rjk,rkl,rneigh(ni,i),ril,rneigh(nl,j),cosphi,cosp1d,cosp2d, &
                                               .true.,.false.)
!
!  Compute multiple angle terms
!
                            cos2phi = 2.0_dp*cosphi*cosphi - 1.0_dp
                            cos3phi = cosphi*(4.0_dp*cosphi*cosphi - 3.0_dp)
!
!  Compute functions f10, f11 and f12
!
                            call reaxFF_f10(BOik,BOij,BOjl,f10,df10dBOik,df10dBOij,df10dBOjl, &
                                            d2f10dBOik2,d2f10dBOikdBOij,d2f10dBOikdBOjl,d2f10dBOik2,d2f10dBOikdBOjl, &
                                            d2f10dBOjl2,.true.,.false.)
                            call reaxFF_f11(nspeci,nspecj,delta(i),delta(j),f11,df11ddeltaij,d2f11ddeltaij2,.true.,.false.)
                            call reaxFF_f12(BOik,BOij,BOjl,f12,df12dBOik,df12dBOij,df12dBOjl, &
                                            d2f12dBOik2,d2f12dBOikdBOij,d2f12dBOikdBOjl,d2f12dBOij2, &
                                            d2f12dBOijdBOjl,d2f12dBOjl2,lStrict,.true.,.false.)
!
!  Compute phi related terms
!
                            expV2 = exp(ptor1*(2.0_dp - BO_pi(ni,i) - f11)**2)
                            V1trm = V1*(1.0_dp + cosphi)
                            V2trm = V2*expV2*(1.0_dp - cos2phi)
                            V3trm = V3*(1.0_dp + cos3phi)
                            conjtrm1 = (cosphi**2 - 1.0_dp)
                            conjtrm2 = (1.0_dp + conjtrm1*sinkij*sinijl)
                            etorsconjtrm = f10*sinkij*sinijl*(V1trm + V2trm + V3trm) + f12*pcot1*conjtrm2
!
!  Compute torsional and conjugation energies
!
                            if (ltoriGTj) then
!
!  Derivatives of torsional angles
!
                              dcos2phidcosphi = 4.0_dp*cosphi
                              dcos3phidcosphi = 12.0_dp*cosphi*cosphi - 3.0_dp
!
!  Derivatives of energy with respect to multiple angles
!
                              detorsdcos1phi = f10*sinkij*sinijl*V1
                              detorsdcos2phi = - f10*sinkij*sinijl*V2*expV2
                              detorsdcos3phi = f10*sinkij*sinijl*V3
!
!  Finish torsional angle derivatives + taper derivatives where needed
!
                              detorsdcosphi = (detorsdcos1phi + detorsdcos2phi*dcos2phidcosphi + detorsdcos3phi*dcos3phidcosphi)
                              deconjdcosphi = 2.0_dp*f12*pcot1*cosphi*sinkij*sinijl
                              detorsdcosphi = detorsdcosphi + deconjdcosphi
                              detorsdcosphi = fij*fik*fjl*detorsdcosphi
!  rik
                              d1i(nk) = d1i(nk) + detorsdcosphi*cosp1d(1)
!  rjk
                              if (ni.ge.nk) then
                                indjk = nneigh(i) + ni*(ni - 1)/2 + nk
                              else
                                indjk = nneigh(i) + nk*(nk - 1)/2 + ni
                              endif
                              d1i(indjk) = d1i(indjk) + detorsdcosphi*cosp1d(2)
!  rkl
                              indkl = nneigh(i) + nneigh(i)*(nneigh(i) + 1)/2 + (ni - 1)*(mneigh+1)*mneigh + nk*mneigh + nl
                              d1i(indkl) = d1i(indkl) + detorsdcosphi*cosp1d(3)
!  rij
                              d1i(ni) = d1i(ni) + detorsdcosphi*cosp1d(4)
!  ril
                              indil = nneigh(i) + nneigh(i)*(nneigh(i) + 1)/2 + (ni - 1)*(mneigh+1)*mneigh + nl
                              d1i(indil) = d1i(indil) + detorsdcosphi*cosp1d(5)
!  rjl
                              indjl = nneigh(i) + nneigh(i)*(nneigh(i) + 1)/2 + (ni - 1)*(mneigh+1)*mneigh + ni*mneigh + nl
                              d1i(indjl) = d1i(indjl) + detorsdcosphi*cosp1d(6) 
!
!  Derivatives with respect to sine of angles
!
                              detorsdsinkij = f10*sinijl*(V1trm + V2trm + V3trm)
                              detorsdsinijl = f10*sinkij*(V1trm + V2trm + V3trm)
                              deconjdsinkij = f12*pcot1*conjtrm1*sinijl
                              deconjdsinijl = f12*pcot1*conjtrm1*sinkij
                              detorsdsinkij = detorsdsinkij + deconjdsinkij
                              detorsdsinijl = detorsdsinijl + deconjdsinijl
!
                              detorsdsinkij = fij*fik*fjl*detorsdsinkij
                              detorsdsinijl = fij*fik*fjl*detorsdsinijl
!  rij
                              d1i(ni) = d1i(ni) + detorsdsinkij*dsinkijdrij + detorsdsinijl*dsinijldrij
!  rik
                              d1i(nk) = d1i(nk) + detorsdsinkij*dsinkijdrik
!  ril
                              d1i(indil) = d1i(indil) + detorsdsinijl*dsinijldril
!  rjk
                              d1i(indjk) = d1i(indjk) + detorsdsinkij*dsinkijdrjk
!  rjl
                              d1i(indjl) = d1i(indjl) + detorsdsinijl*dsinijldrjl
!
!  Add on contributions to bond order derivatives for i-j & i-k
!
                              detorsdf10 = fij*fik*fjl*sinkij*sinijl*(V1trm + V2trm + V3trm)
                              detorsdBO_neigh(ni) = detorsdBO_neigh(ni) + detorsdf10*df10dBOij
                              detorsdBO_neigh(nk) = detorsdBO_neigh(nk) + detorsdf10*df10dBOik
                              deconjdf12 = fij*fik*fjl*pcot1*conjtrm2
                              detorsdBO_neigh(ni) = detorsdBO_neigh(ni) + deconjdf12*df12dBOij
                              detorsdBO_neigh(nk) = detorsdBO_neigh(nk) + deconjdf12*df12dBOik
!
!  Derivatives of taper functions with respect to bond orders
!
                              detorsdBO_neigh(ni) = detorsdBO_neigh(ni) + etorsconjtrm*fik*fjl*dfijdBOij
                              detorsdBO_neigh(nk) = detorsdBO_neigh(nk) + etorsconjtrm*fij*fjl*dfikdBOik
!
!  Add on contributions to pi bond order derivatives for i-j from expV2
!
                              detorsdexpV2 = fij*fik*fjl*f10*sinkij*sinijl*V2*(1.0_dp - cos2phi)
                              dexpV2dBO_pi = - 2.0_dp*ptor1*expV2*(2.0_dp - BO_pi(ni,i) - f11)
                              detorsdBOpi_neigh(ni) = detorsdBOpi_neigh(ni) + detorsdexpV2*dexpV2dBO_pi
!
!  Add on contributions to derivatives with respect to delta_i
!
                              dexpV2df11 = - 2.0_dp*ptor1*expV2*(2.0_dp - BO_pi(ni,i) - f11)
                              detorsddeltai = detorsddeltai + detorsdexpV2*dexpV2df11*df11ddeltaij
                            else
!
!  Add on contributions to bond order derivatives that are effectively for j-l but here are i-k
!
                              detorsdf10 = fij*fik*fjl*sinkij*sinijl*(V1trm + V2trm + V3trm)
                              detorsdBO_neigh(nk) = detorsdBO_neigh(nk) + detorsdf10*df10dBOik
                              deconjdf12 = fij*fik*fjl*pcot1*conjtrm2
                              detorsdBO_neigh(nk) = detorsdBO_neigh(nk) + deconjdf12*df12dBOik
!
!  Add on taper function derivatives with respect to bond order
!
                              detorsdBO_neigh(nk) = detorsdBO_neigh(nk) + etorsconjtrm*fij*fjl*dfikdBOik
!
!  Add on contributions to derivatives with respect to delta_j
!
                              detorsdexpV2 = fij*fik*fjl*f10*sinkij*sinijl*V2*(1.0_dp - cos2phi)
                              dexpV2df11 = - 2.0_dp*ptor1*expV2*(2.0_dp - BO_pi(ni,i) - f11)
                              detorsddeltai = detorsddeltai + detorsdexpV2*dexpV2df11*df11ddeltaij
                            endif
                          endif
                        endif
!
!  End of check on bond order for j-l
!
                      endif
                    enddo
!
!  End of check on bond order for i-k
!
                  endif
!
                endif
              enddo
!
!  End of check on bond order for i-j
!
            endif
!
          endif
        enddo
!
!  Now for the derivatives......
!
!  Set up derivative of deltalp_corr_i
!
        ddeltalpcorrddeltai = 1.0_dp - otrm1*ddeltalpddeltae
!
!  Set up term needed later for derivative of deltalp_corr_i with respect to pi bond orders
!
        otrm4 = (otrm1**2)*expdi*reaxFFlam(6)*reaxFFlam(31)*deltalp(i)
!
!  Loop over neighbours of i to compute the delta_lpcorr values for each neighbour
!
        do ni = 1,nneigh(i)
          j = neighno(ni,i)
          if (nbos(j).gt.0) then
            nspecj = nbosptr(j)
!
!  Loop over neighbours of j to sum the products of terms required to compute contribution of j to eover
!
            sumdeltabopi_j = 0.0_dp
            sumover_j = 0.0_dp
            do nj = 1,nneigh(j)
              k = neighno(nj,j)
              if (nbos(k).gt.0) then
                nspeck = nbosptr(k)
                if (nspecj.ge.nspeck) then
                  nbojk = nspecj*(nspecj - 1)/2 + nspeck
                else
                  nbojk = nspeck*(nspeck - 1)/2 + nspecj
                endif
                if (.not.lreaxFFunder(nspecj).or..not.lreaxFFunder(nspeck)) then
                  sumdeltabopi_j = sumdeltabopi_j + delta(k)*(BO_pi(nj,j) + BO_pipi(nj,j))
                else
                  sumdeltabopi_j = sumdeltabopi_j + (delta(k) - deltalp(k))*(BO_pi(nj,j) + BO_pipi(nj,j))
                endif
                sumover_j = sumover_j + reaxFFoc2(nbojk)*reaxFFDe(1,nbojk)*BO(nj,j)
              endif
            enddo
!
!  Compute delta_lpcorr for j
!
            if (.not.lreaxFFunder(nspecj)) then
              expdi_j = 0.0_dp
              otrm1_j = 0.0_dp
              deltalp_corr_j = delta(j) 
            else
              expdi_j = exp(reaxFFlam(31)*sumdeltabopi_j)
              otrm1_j = 1.0_dp/(1.0_dp + reaxFFlam(6)*expdi_j)
              deltalp_corr_j = delta(j) - deltalp(j)*otrm1_j
            endif
!
!  Compute derivative of eover term for j with respect to deltalp_corr_j
!
            otrm2_j = 1.0_dp/(deltalp_corr_j + reaxFFval(1,nspecj) + 1.0d-8)
            expoc_j = exp(reaxFFoc1(nspecj)*deltalp_corr_j)
            otrm3_j = 1.0_dp/(1.0_dp + expoc_j)
            deoverddeltalpcorr_j = sumover_j*otrm2_j*otrm3_j*(1.0_dp - deltalp_corr_j* &
              (otrm2_j + reaxFFoc1(nspecj)*otrm3_j*expoc_j))
!
!  Compute derivative of deltalp_corr_j with respect to delta_j and deltalp_j
!
            otrm4_j = (otrm1_j**2)*expdi_j*reaxFFlam(6)*reaxFFlam(31)*deltalp(j)
            ddeltalpcorrjddeltaj = otrm4_j*(BO_pi(ni,i) + BO_pipi(ni,i))
            if (.not.lreaxFFunder(nspecj)) then
              ddeltalpcorrjddeltalpj = 0.0_dp
            else
              ddeltalpcorrjddeltalpj = - otrm4_j*(BO_pi(ni,i) + BO_pipi(ni,i))
            endif
!
!  Compute derivatives of eunder term for j with respect to deltalp_corr_j
!
            expuc1_j = exp(reaxFFlam(7)*deltalp_corr_j)
            expuc2_j = exp(-reaxFFoc1(nspecj)*deltalp_corr_j)
            utrm1_j  = 1.0_dp/(1.0_dp + expuc2_j)
            expuc3_j = exp(reaxFFlam(9)*sumdeltabopi_j)
            utrm2_j  = 1.0_dp/(1.0_dp + reaxFFlam(8)*expuc3_j)
            deunderddeltalpcorr_j = - reaxFFuc1(nspecj)*utrm1_j*utrm2_j* &
              ((1.0_dp - expuc1_j)*utrm1_j*reaxFFoc1(nspecj)*expuc2_j - reaxFFlam(7)*expuc1_j)
!
!  Compute derivatives of eunder term for j with respect to delta_j and deltalp_j
!
            deunderutrm2ddeltaj = reaxFFuc1(nspecj)*(1.0_dp - expuc1_j)*utrm1_j*(utrm2_j**2)* &
              reaxFFlam(8)*reaxFFlam(9)*expuc3_j*(BO_pi(ni,i) + BO_pipi(ni,i))
!
!  Store derivative terms for neighbours of i - eover contribution first...
!
            deoverunderddeltaj_neigh(ni) = deoverddeltalpcorr_j*(ddeltalpcorrjddeltaj + ddeltalpcorrjddeltalpj*ddeltalpddeltae)
!
!  ... now add on eunder contribution.
!
            if (.not.lreaxFFunder(nspecj)) then
              deoverunderddeltaj_neigh(ni) = deoverunderddeltaj_neigh(ni) + &
                deunderddeltalpcorr_j*(ddeltalpcorrjddeltaj + ddeltalpcorrjddeltalpj*ddeltalpddeltae) + &
                deunderutrm2ddeltaj
            else
              deoverunderddeltaj_neigh(ni) = deoverunderddeltaj_neigh(ni) + &
                deunderddeltalpcorr_j*(ddeltalpcorrjddeltaj + ddeltalpcorrjddeltalpj*ddeltalpddeltae) + &
                deunderutrm2ddeltaj*(1.0_dp - ddeltalpddeltae)
            endif
          endif
        enddo
!
!  Set up terms for derivatives of lone pair energy
!
        if (abs(reaxFFlp(2,nspeci)).gt.1.0d-12) then
          deloneddeltalp = reaxFFlp(2,nspeci)*trm1lp*(trm1lp*deltalp(i)*explp*reaxFFlam(30) + 1.0_dp)
          deloneddeltai = deloneddeltalp*ddeltalpddeltae
        else
          deloneddeltai = 0.0_dp
        endif
!
!  Derivatives of Bond penalty energy
!
        debpenddeltai = 0.0_dp
        do ni = 1,nneigh(i)
          j = neighno(ni,i)
          debpendBO_neigh(ni) = 0.0_dp
          if (nbos(j).gt.0) then
            nspecj = nbosptr(j)
            if (nspeci.ge.nspecj) then
              nboij = nspeci*(nspeci - 1)/2 + nspecj
            else
              nboij = nspecj*(nspecj - 1)/2 + nspeci
            endif
            diffBOdelta = BO(ni,i) - delta(i) - reaxFFpen2(2,nboij)*delta(i)**4 - reaxFFpen2(3,nboij)
            if (diffBOdelta.gt.0.0_dp) then
              debpenddeltai = debpenddeltai - 2.0_dp*reaxFFpen2(1,nboij)*diffBOdelta* &
                (1.0_dp + 4.0_dp*reaxFFpen2(2,nboij)*delta(i)**3)
              debpendBO_neigh(ni) = debpendBO_neigh(ni) + 2.0_dp*reaxFFpen2(1,nboij)*diffBOdelta
            endif
          endif
        enddo
!
!  Set up terms for derivatives of overcoordination energy
!
!  NB: Factor of p_ovun1.De_sigma not included in deoverdBOij since it's bond specific
!
        deoverdBOij = otrm2*deltalp_corr_i*otrm3
        deoverddeltalpcorr = sumover*otrm2*otrm3*(1.0_dp - deltalp_corr_i*(otrm2 + reaxFFoc1(nspeci)*otrm3*expoc))
        deoverddeltai = deoverddeltalpcorr*ddeltalpcorrddeltai
!
!  Set up terms for derivatives of eunder
!
        deunderddeltalpcorr = - reaxFFuc1(nspeci)*utrm1*utrm2* &
          ((1.0_dp - expuc1)*utrm1*reaxFFoc1(nspeci)*expuc2 - reaxFFlam(7)*expuc1)
        deunderddeltai = deunderddeltalpcorr*ddeltalpcorrddeltai
        deunderutrm2deltaBO = reaxFFuc1(nspeci)*(1.0_dp - expuc1)*utrm1*utrm2*utrm2*reaxFFlam(8)*reaxFFlam(9)*expuc3
!
!  Loop over neighbours of i (=> j)
!
        do ni = 1,nneigh(i)
          j = neighno(ni,i)
          nregionj = nregionno(nsft + nrelf2a(j))
          nregiontypj = nregiontype(nregionj,ncf)
!
!  Do we need to do this pair of atoms
!
          if (nbos(j).gt.0) then
!
!  Set variables relating to j
!
            nspecj = nbosptr(j)
!
!  Set up i-j quantities
!
            lreg2one  = .false.
            lreg2pair = .false.
!
!  Find i in neighbour list for j
!
            nj = neighnoRptr(ni,i)
!
!  Set total number of distances for neighbours of j
!
            nneighj1 = nneigh(j)*(nneigh(j) + 1)/2
            nneighj2 = nneigh(j) + nneighj1
!
!  Initialise derivative storage for neighbours of i
!
            d1j(1:nneighj2) = 0.0_dp
!
!  Set index to parameters for pairwise interaction of species
!
            if (nspeci.ge.nspecj) then
              nboij = nspeci*(nspeci - 1)/2 + nspecj
            else
              nboij = nspecj*(nspecj - 1)/2 + nspeci
            endif
!
!  Set delta' 
!
            deltapi = deltap(i)
            deltapj = deltap(j)
!***************************************
!  Valid pairwise reaxFF contribution  *
!***************************************
            call reaxFF_bo(i,j,nspeci,nspecj,ni,nj,deltapi,deltapj,nneigh,BOp,BOp_pi,BOp_pipi, &
                           d1BOp,d1BOp_pi,d1BOp_pipi,BOij_s,BOij_pi,BOij_pipi, &
                             d1BOi_s,d1BOi_pi,d1BOi_pipi,d1BOj_s,d1BOj_pi,d1BOj_pipi,.true.)
!
!  Calculate total bond order derivatives from sum of components
!
            do k = 1,nneighi2
              d1BOi(k) = d1BOi_s(k) + d1BOi_pi(k) + d1BOi_pipi(k)
            enddo
            do k = 1,nneighj2
              d1BOj(k) = d1BOj_s(k) + d1BOj_pi(k) + d1BOj_pipi(k)
            enddo
!
!  Derivatives of over coordination energy with respect to BOij
!
            do k = 1,nneighi2
              d1i(k) = d1i(k) + deoverdBOij*reaxFFoc2(nboij)*reaxFFDe(1,nboij)*d1BOi(k)
            enddo
            do k = 1,nneighj2
              d1j(k) = d1j(k) + deoverdBOij*reaxFFoc2(nboij)*reaxFFDe(1,nboij)*d1BOj(k)
            enddo
!
!  Derivatives of over-/under-coordination energy with respect to terms summed over neighbours (delta_j - deltalp_j)*(BO_pi + BO_pipi)
!
            sum_deoverunderddeltaj_neigh = 0.0_dp
            do l = 1,nneigh(i)
              sum_deoverunderddeltaj_neigh = sum_deoverunderddeltaj_neigh + deoverunderddeltaj_neigh(l)
            enddo
            do k = 1,nneighi2
              d1i(k) = d1i(k) + sum_deoverunderddeltaj_neigh*d1BOi(k)
            enddo
            do k = 1,nneighj2
              d1j(k) = d1j(k) + sum_deoverunderddeltaj_neigh*d1BOj(k)
            enddo
!
!  Derivatives of deltalp_corr_i with respect to pi & pi-pi bond orders
!
            if (.not.lreaxFFunder(nspeci)) then
              deoverdBOpi   = deoverddeltalpcorr*otrm4*delta(j)
            else
              deoverdBOpi   = deoverddeltalpcorr*otrm4*(delta(j) - deltalp(j))
            endif
            do k = 1,nneighi2
              d1i(k) = d1i(k) + deoverdBOpi*(d1BOi_pi(k) + d1BOi_pipi(k))
            enddo
            do k = 1,nneighj2
              d1j(k) = d1j(k) + deoverdBOpi*(d1BOj_pi(k) + d1BOj_pipi(k))
            enddo
            if (.not.lreaxFFunder(nspeci).or..not.lreaxFFunder(nspecj)) then
              deunderdBOpi   = (deunderddeltalpcorr*otrm4 + deunderutrm2deltaBO)*delta(j)
            else
              deunderdBOpi   = (deunderddeltalpcorr*otrm4 + deunderutrm2deltaBO)*(delta(j) - deltalp(j))
            endif
            do k = 1,nneighi2
              d1i(k) = d1i(k) + deunderdBOpi*(d1BOi_pi(k) + d1BOi_pipi(k))
            enddo
            do k = 1,nneighj2
              d1j(k) = d1j(k) + deunderdBOpi*(d1BOj_pi(k) + d1BOj_pipi(k))
            enddo
!
!  Derivatives of lonepair, bond penalty, over- and under-coordination energies with respect to delta(i)
!
            do k = 1,nneighi2
              d1i(k) = d1i(k) + (deloneddeltai + deoverddeltai + deunderddeltai + debpenddeltai)*d1BOi(k)
            enddo
            do k = 1,nneigh(j)
              ind = nneigh(i) + nneigh(i)*(nneigh(i) + 1)/2 + (ni - 1)*(mneigh+1)*mneigh + ni*mneigh + k
              d1i(ind) = d1i(ind) + (deloneddeltai + deoverddeltai + deunderddeltai + debpenddeltai)*d1BOj(k)
            enddo
!
!  Derivatives of ebpen with respect to BO_ij
!
            do k = 1,nneighi2
              d1i(k) = d1i(k) + debpendBO_neigh(ni)*d1BOi(k)
            enddo
            do k = 1,nneighj2
              d1j(k) = d1j(k) + debpendBO_neigh(ni)*d1BOj(k)
            enddo
!
!  Derivatives of eval with respect to delta(i)
!
            do k = 1,nneighi2
              d1i(k) = d1i(k) + devalddeltai*d1BOi(k)
            enddo
            do k = 1,nneighj2
              d1j(k) = d1j(k) + devalddeltai*d1BOj(k)
            enddo
!
!  Derivatives of eval with respect to BO_ij
!
            do k = 1,nneighi2
              d1i(k) = d1i(k) + devaldBO_neigh(ni)*d1BOi(k)
            enddo
            do k = 1,nneighj2
              d1j(k) = d1j(k) + devaldBO_neigh(ni)*d1BOj(k)
            enddo
!
!  Derivatives of eval with respect sbo multiplied by derivatives of pi bond orders 
!
            do k = 1,nneighi2
              d1i(k) = d1i(k) + devaldsbo*(d1BOi_pi(k) + d1BOi_pipi(k))
            enddo
            do k = 1,nneighj2
              d1j(k) = d1j(k) + devaldsbo*(d1BOj_pi(k) + d1BOj_pipi(k))
            enddo
!
!  Derivatives of eval with respect sbo multiplied by derivatives of bond orders 
!
            do k = 1,nneighi2
              d1i(k) = d1i(k) + devaldsbo_neigh(ni)*d1BOi(k)
            enddo
            do k = 1,nneighj2
              d1j(k) = d1j(k) + devaldsbo_neigh(ni)*d1BOj(k)
            enddo
!
!  Derivatives of etors with respect to BO_ij
!
            do k = 1,nneighi2
              d1i(k) = d1i(k) + detorsdBO_neigh(ni)*d1BOi(k)
            enddo
            do k = 1,nneighj2
              d1j(k) = d1j(k) + detorsdBO_neigh(ni)*d1BOj(k)
            enddo
!
!  Derivatives of etors with respect to BOpi_ij
!
            do k = 1,nneighi2
              d1i(k) = d1i(k) + detorsdBOpi_neigh(ni)*d1BOi_pi(k)
            enddo
            do k = 1,nneighj2
              d1j(k) = d1j(k) + detorsdBOpi_neigh(ni)*d1BOj_pi(k)
            enddo
!
!  Derivatives of etors with respect to delta(i)
!
            do k = 1,nneighi2
              d1i(k) = d1i(k) + detorsddeltai*d1BOi(k)
            enddo
            do k = 1,nneighj2
              d1j(k) = d1j(k) + detorsddeltai*d1BOj(k)
            enddo
!
!  Add derivatives due to neighbours of j 
!
            call d1addfc(j,maxneigh,mneigh,nneigh,neighno,xneigh,yneigh,zneigh,ineigh,maxlhs, &
                         d1cell,matom,nbodummyptr,d1j,.false.)
!
!  End condition section on i or j being associated with moving atom
!
          endif
!
!  End of loop over neighbours of i
!
        enddo
!
!  Add derivatives due to neighbours of i - here we include torsional derivatives
!
        call d1addfc(i,maxneigh,mneigh,nneigh,neighno,xneigh,yneigh,zneigh,ineigh,maxlhs, &
                     d1cell,matom,nbodummyptr,d1i,.true.)
!
!  End of QMMM condition
!
      endif
!
!  End of condition on i having a ReaxFF species type
!
    endif
  enddo
!**********************************
!  Hydrogen bonding contribution  *
!**********************************
  cut2 = reaxFFrhtol**2
!
!  Loop over atoms looking for hydrogen
!
  do i = 1,numat
    if (nbos(i).gt.0.and.nat(i).eq.1) then
!
!  Set variables relating to i
!
      nspeci = nbosptr(i)
      nregioni = nregionno(nsft + nrelf2a(i))
      nregiontypi = nregiontype(nregioni,ncf)
      deltapi = deltap(i)
!
!  QM/MM handling : i is a QM atom => exclude
!
      lQMMMok = .true.
      if (QMMMmode(ncf).gt.0) then
        if (nregiontypi.eq.1) lQMMMok = .false.
      endif
!
!  Do we need to do this atom based on QMMM scheme?
!
      if (lQMMMok) then
!           
!  Set total number of distances for neighbours of i
!             
        nneighi1   = nneigh(i)*(nneigh(i) + 1)/2
        nneighi2   = nneigh(i) + nneighi1 
        nneighi2i  = nneigh(i) + nneighi1 + nneigh(i)*(mneigh+1)*mneigh
!           
!  Initialise derivative storage for neighbours of i
!       
        d1i(1:nneighi2i) = 0.0_dp
!
!  Loop over atoms that are connected by bond order to i
!
        do ni = 1,nneigh(i)
          j = neighno(ni,i)
          if (nbos(j).gt.0) then
            nspecj = nbosptr(j)
!
!  Is the bond order above the threshold for a hydrogen bond?
!
            if (BO(ni,i).gt.reaxFFhtol) then
!
!  Is there any H-bonding term for this i-j pair?
!
              lAnyHB = .false.
              do nspeck = 1,nreaxFFspec
                phb1 = reaxFFhb3(2,nspeck,nspecj,nspeci)
                if (abs(phb1).gt.1.0d-12) lAnyHB = .true.
              enddo
              if (lAnyHB) then
                xji = xneigh(ni,i)
                yji = yneigh(ni,i)
                zji = zneigh(ni,i)
                rji = rneigh(ni,i)
!                 
!  Calculate taper function to smooth cutoff
!                 
                if (lStrict) then
                  fHB = 1.0_dp
                  dfHBdBO = 0.0_dp
                else 
                  call ataper(.false.,BO(ni,i),reaxFFhtol,reaxFFtaperscale*reaxFFhtol,fHB,dfHBdBO,d2fHBdBO2,d3fHBdBO3, &
                              .true.,.false.,.false.)
                endif
!
!  Set region for j
!
                nregionj = nregionno(nsft + nrelf2a(j))
!
!  Find i in neighbour list for j
!
                nj = neighnoRptr(ni,i)
!
!  Set total number of distances for neighbours of j
!
                nneighj1 = nneigh(j)*(nneigh(j) + 1)/2
                nneighj2 = nneigh(j) + nneighj1
!
!  Initialise derivative storage for neighbours of i
!
                d1j(1:nneighj2) = 0.0_dp
!
!  Set delta'
!
                deltapj = deltap(j)
!
!  Compute bond order derivatives
!
                call reaxFF_bo(i,j,nspeci,nspecj,ni,nj,deltapi,deltapj,nneigh,BOp,BOp_pi,BOp_pipi, &
                               d1BOp,d1BOp_pi,d1BOp_pipi,BOij_s,BOij_pi,BOij_pipi, &
                               d1BOi_s,d1BOi_pi,d1BOi_pipi,d1BOj_s,d1BOj_pi,d1BOj_pipi,.true.)
!
!  Calculate total bond order derivatives from sum of components
!
                do k = 1,nneighi2
                  d1BOi(k) = d1BOi_s(k) + d1BOi_pi(k) + d1BOi_pipi(k)
                enddo
                do k = 1,nneighj2
                  d1BOj(k) = d1BOj_s(k) + d1BOj_pi(k) + d1BOj_pipi(k)
                enddo
!
!  One of the atoms must be matom & so set range for k accordingly
!
                if (i.ne.matom.and.j.ne.matom) then
                  kmin = 1
                  kmax = numat
                else
                  kmin = matom
                  kmax = matom
                endif
!
!  Loop over atoms looking for an acceptor
!
                kloop: do k = kmin,kmax
!
!  Does k have a reaxFF species?
!
                  if (nbos(k).gt.0) then
                    nspeck = nbosptr(k)
!
!  Is this triad valid for a hydrogen bond?
!
                    phb1 = reaxFFhb3(2,nspeck,nspecj,nspeci)
                    if (abs(phb1).lt.1.0d-12) cycle kloop
!
!  Compute basic interatomic vector
!
                    xkj = xclat(k) - xclat(j)
                    ykj = yclat(k) - yclat(j)
                    zkj = zclat(k) - zclat(j)
!
!  Compute shift if relevant
!
                    if (j.eq.matom) then
                      xshift = vmatom(1) - xclat(j)
                      yshift = vmatom(2) - yclat(j)
                      zshift = vmatom(3) - zclat(j)
                    elseif (k.eq.matom) then
                      xshift = vmatom(1) - xclat(k)
                      yshift = vmatom(2) - yclat(k)
                      zshift = vmatom(3) - zclat(k)
                    endif
!
!  Find valid vectors
!
                    nor = 0
                    if (ndim.eq.3) then
                      call rsearch3D(xkj,ykj,zkj,.false.,.false.,j,k,0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2)
                    elseif (ndim.eq.2) then
                      call rsearch2D(xkj,ykj,zkj,.false.,.false.,j,k,0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2)
                    elseif (ndim.eq.1) then
                      call rsearch1D(xkj,ykj,zkj,.false.,.false.,j,k,0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2)
                    elseif (ndim.eq.0) then
                      rkj2 = xkj*xkj + ykj*ykj + zkj*zkj
                      if (rkj2.lt.cut2) then
                        nor = nor + 1
                      endif
                    endif
!
!  If no distances then cycle
!
                    if (nor.eq.0) cycle kloop
!
!  Correct distances for displacement
!
                    if (ndim.gt.0) then
                      if (j.eq.matom) then
                        do n = 1,nor
                          xtmp(n) = xtmp(n) - xshift
                          ytmp(n) = ytmp(n) - yshift
                          ztmp(n) = ztmp(n) - zshift
                          dist(n) = xtmp(n)*xtmp(n) + ytmp(n)*ytmp(n) + ztmp(n)*ztmp(n)
                        enddo
                      elseif (k.eq.matom) then
                        do n = 1,nor
                          xtmp(n) = xtmp(n) + xshift
                          ytmp(n) = ytmp(n) + yshift
                          ztmp(n) = ztmp(n) + zshift
                          dist(n) = xtmp(n)*xtmp(n) + ytmp(n)*ytmp(n) + ztmp(n)*ztmp(n)
                        enddo
                      endif
                    else
                      if (j.eq.matom) then
                        xkj = xkj - xshift
                        ykj = ykj - yshift
                        zkj = zkj - zshift
                      elseif (k.eq.matom) then
                        xkj = xkj + xshift
                        ykj = ykj + yshift
                        zkj = zkj + zshift
                      endif
                      rkj2 = xkj*xkj + ykj*ykj + zkj*zkj
                    endif
!
!  Set remaining parameters
!
                    r0hb = reaxFFhb3(1,nspeck,nspecj,nspeci)
                    phb2 = reaxFFhb3(3,nspeck,nspecj,nspeci)
                    phb3 = reaxFFhb3(4,nspeck,nspecj,nspeci)
!
!  Set region for k
!
                    nregionk = nregionno(nsft + nrelf2a(k))
!
!  Loop over valid distances and calculate contributions
!
                    nloop: do n = 1,nor
                      if (ndim.gt.0) then
                        rkj2 = dist(n)
                        xkj = xtmp(n)
                        ykj = ytmp(n)
                        zkj = ztmp(n)
                      endif
!
!  Find cell index for i-j
!
                      if (abs(ineigh(1,ni,i)).gt.nd2cell(1).or. &
                          abs(ineigh(2,ni,i)).gt.nd2cell(2).or. &
                          abs(ineigh(3,ni,i)).gt.nd2cell(3)) then
!
!  Find cell index - if outside user range then assign to central cell
!
                        ncindmij = nd2central
                        ncindpij = nd2central
                      else
!
!  Compute index
!
                        ncindpij = nd2cellptr(nd2cell(1)+1+ineigh(1,ni,i), &
                                              nd2cell(2)+1+ineigh(2,ni,i), &
                                              nd2cell(3)+1+ineigh(3,ni,i))
                        ncindmij = nd2cellptr(nd2cell(1)+1-ineigh(1,ni,i), &
                                              nd2cell(2)+1-ineigh(2,ni,i), &
                                              nd2cell(3)+1-ineigh(3,ni,i))
                      endif
!
!  Find cell index for j-k
!
                      if (abs(cellindex(1,n)).gt.nd2cell(1).or. &
                          abs(cellindex(2,n)).gt.nd2cell(2).or. &
                          abs(cellindex(3,n)).gt.nd2cell(3)) then
!
!  Find cell index - if outside user range then assign to central cell
!
                        ncindmjk = nd2central
                        ncindpjk = nd2central
                      else
!
!  Compute index
!
                        ncindpjk = nd2cellptr(nd2cell(1)+1+cellindex(1,n),nd2cell(2)+1+cellindex(2,n),nd2cell(3)+1+cellindex(3,n))
                        ncindmjk = nd2cellptr(nd2cell(1)+1-cellindex(1,n),nd2cell(2)+1-cellindex(2,n),nd2cell(3)+1-cellindex(3,n))
                      endif
!
!  Find cell index for i-k
!
                      if (abs(cellindex(1,n)+ineigh(1,ni,i)).gt.nd2cell(1).or. &
                          abs(cellindex(2,n)+ineigh(2,ni,i)).gt.nd2cell(2).or. &
                          abs(cellindex(3,n)+ineigh(3,ni,i)).gt.nd2cell(3)) then
!
!  Find cell index - if outside user range then assign to central cell
!
                        ncindmik = nd2central
                        ncindpik = nd2central
                      else
!
!  Compute index
!
                        ncindpik = nd2cellptr(nd2cell(1)+1+cellindex(1,n)+ineigh(1,ni,i), &
                                              nd2cell(2)+1+cellindex(2,n)+ineigh(2,ni,i), &
                                              nd2cell(3)+1+cellindex(3,n)+ineigh(3,ni,i))
                        ncindmik = nd2cellptr(nd2cell(1)+1-cellindex(1,n)-ineigh(1,ni,i), &
                                              nd2cell(2)+1-cellindex(2,n)-ineigh(2,ni,i), &
                                              nd2cell(3)+1-cellindex(3,n)-ineigh(3,ni,i))
                      endif
!
!  Calculate distances for i - k
!
                      xki = xkj + xji
                      yki = ykj + yji
                      zki = zkj + zji
                      rki2 = xki*xki + yki*yki + zki*zki
!
!  Exclude i = k case and j = k case
!
                      if (i.eq.k.and.rki2.lt.1.0d-12) cycle nloop
                      if (j.eq.k.and.rkj2.lt.1.0d-12) cycle nloop
!
                      rki = sqrt(rki2)
                      rkj = sqrt(rkj2)
!                         
!  Compute distance taper function 
!
                      if (lStrict) then
                        frHB = 1.0_dp
                        dfrHBdr = 0.0_dp
                      else
                        call ataper(.true.,rkj,reaxFFrhtollower,reaxFFrhtol,frHB,dfrHBdr,d2frHBdr2,d3frHBdr3, &
                                    .true.,.false.,.false.)
                      endif
!  
!  Compute theta for j - i - k
!                     
                      call reaxff_theta(rji,rki,rkj,theta,dthetadrij,dthetadrik,dthetadrjk, &
                                        d2thetadr2,.true.,.false.)
!
!  Compute remaining terms for hydrogen bond energy
!
                      rrik = 1.0_dp/rki
                      rrij = 1.0_dp/rji
                      sinhalftheta  = sin(0.5_dp*theta)
                      sin4halftheta = sinhalftheta**4
                      exphb1 = exp(-phb2*BO(ni,i))
                      exphb2 = exp(-phb3*((r0hb*rrik) + (rki/r0hb) - 2.0_dp))
!
!  Hydrogen bond derivatives with respect to the distance i-k
!
                      dehbdrik = - rrik*fHB*frHB*phb1*phb3*(1.0_dp - exphb1)*exphb2*sin4halftheta* &
                                   (1.0_dp/r0hb - r0hb*(rrik**2))
!
                      if (k.eq.matom) then
                        d1cell(1,i,ncindmik) = d1cell(1,i,ncindmik) - dehbdrik*xki
                        d1cell(2,i,ncindmik) = d1cell(2,i,ncindmik) - dehbdrik*yki
                        d1cell(3,i,ncindmik) = d1cell(3,i,ncindmik) - dehbdrik*zki
                      endif
                      if (i.eq.matom) then
                        d1cell(1,k,ncindpik) = d1cell(1,k,ncindpik) + dehbdrik*xki
                        d1cell(2,k,ncindpik) = d1cell(2,k,ncindpik) + dehbdrik*yki
                        d1cell(3,k,ncindpik) = d1cell(3,k,ncindpik) + dehbdrik*zki
                      endif
!
                      if (.not.lStrict) then
!                           
!  Hydrogen bond derivatives with respect to the distance j-k
!                           
                        rrjk = 1.0_dp/rkj
                        dehbdrjk = fHB*(rrjk*dfrHBdr)*phb1*(1.0_dp - exphb1)*exphb2*sin4halftheta
!
                        if (k.eq.matom) then
                          d1cell(1,j,ncindmjk) = d1cell(1,j,ncindmjk) - dehbdrjk*xkj
                          d1cell(2,j,ncindmjk) = d1cell(2,j,ncindmjk) - dehbdrjk*ykj
                          d1cell(3,j,ncindmjk) = d1cell(3,j,ncindmjk) - dehbdrjk*zkj
                        endif
                        if (j.eq.matom) then
                          d1cell(1,k,ncindpjk) = d1cell(1,k,ncindpjk) + dehbdrjk*xkj
                          d1cell(2,k,ncindpjk) = d1cell(2,k,ncindpjk) + dehbdrjk*ykj
                          d1cell(3,k,ncindpjk) = d1cell(3,k,ncindpjk) + dehbdrjk*zkj
                        endif
                      endif
!
!  Hydrogen bond derivatives with respect to the angle
!
                      coshalftheta = cos(0.5_dp*theta)
                      dehbdtheta = fHB*frHB*(phb1*(1.0_dp - exphb1)*exphb2)*(2.0_dp*coshalftheta*sinhalftheta**3)
!
                      if (j.eq.matom) then
                        d1cell(1,i,ncindmij) = d1cell(1,i,ncindmij) - dehbdtheta*dthetadrij*xji
                        d1cell(2,i,ncindmij) = d1cell(2,i,ncindmij) - dehbdtheta*dthetadrij*yji
                        d1cell(3,i,ncindmij) = d1cell(3,i,ncindmij) - dehbdtheta*dthetadrij*zji
                      endif
                      if (i.eq.matom) then
                        d1cell(1,j,ncindpij) = d1cell(1,j,ncindpij) + dehbdtheta*dthetadrij*xji
                        d1cell(2,j,ncindpij) = d1cell(2,j,ncindpij) + dehbdtheta*dthetadrij*yji
                        d1cell(3,j,ncindpij) = d1cell(3,j,ncindpij) + dehbdtheta*dthetadrij*zji
                      endif
!
                      if (k.eq.matom) then
                        d1cell(1,i,ncindmik) = d1cell(1,i,ncindmik) - dehbdtheta*dthetadrij*xki
                        d1cell(2,i,ncindmik) = d1cell(2,i,ncindmik) - dehbdtheta*dthetadrij*yki
                        d1cell(3,i,ncindmik) = d1cell(3,i,ncindmik) - dehbdtheta*dthetadrij*zki
                      endif
                      if (i.eq.matom) then
                        d1cell(1,k,ncindpik) = d1cell(1,k,ncindpik) + dehbdtheta*dthetadrij*xki
                        d1cell(2,k,ncindpik) = d1cell(2,k,ncindpik) + dehbdtheta*dthetadrij*yki
                        d1cell(3,k,ncindpik) = d1cell(3,k,ncindpik) + dehbdtheta*dthetadrij*zki
                      endif
!
                      if (k.eq.matom) then
                        d1cell(1,j,ncindmjk) = d1cell(1,j,ncindmjk) - dehbdtheta*dthetadrij*xkj
                        d1cell(2,j,ncindmjk) = d1cell(2,j,ncindmjk) - dehbdtheta*dthetadrij*ykj
                        d1cell(3,j,ncindmjk) = d1cell(3,j,ncindmjk) - dehbdtheta*dthetadrij*zkj
                      endif
                      if (j.eq.matom) then
                        d1cell(1,k,ncindpjk) = d1cell(1,k,ncindpjk) + dehbdtheta*dthetadrij*xkj
                        d1cell(2,k,ncindpjk) = d1cell(2,k,ncindpjk) + dehbdtheta*dthetadrij*ykj
                        d1cell(3,k,ncindpjk) = d1cell(3,k,ncindpjk) + dehbdtheta*dthetadrij*zkj
                      endif
!
!  Derivatives of hydrogen bond energy with respect to bond orders -> add to eval terms
!
                      dehbdBOij = phb1*frHB*(fHB*phb2*exphb1 + dfHBdBO*(1.0_dp - exphb1))*exphb2*sin4halftheta
                      do l = 1,nneighi2
                        d1i(l) = d1i(l) + dehbdBOij*d1BOi(l)
                      enddo
                      do l = 1,nneighj2
                        d1j(l) = d1j(l) + dehbdBOij*d1BOj(l)
                      enddo
!
                    enddo nloop
!
!  End check on whether k is ReaxFF species
!
                  endif
!
!  End of loop over k
!
                enddo kloop
!
!  Add derivatives due to neighbours of j 
!
                call d1addfc(j,maxneigh,mneigh,nneigh,neighno,xneigh,yneigh,zneigh,ineigh,maxlhs, &
                             d1cell,matom,nbodummyptr,d1j,.false.)
!
!  End of check as to whether there are any H-bonding terms for this i-j pair
!
              endif
!
!  End of check on bond order threshold
!
            endif
!
!  End of check on whether j is a ReaxFF atom
!
          endif
        enddo
!
!  Add derivatives due to neighbours of i 
!
        call d1addfc(i,maxneigh,mneigh,nneigh,neighno,xneigh,yneigh,zneigh,ineigh,maxlhs, &
                     d1cell,matom,nbodummyptr,d1i,.false.)
!
!  End of QMMM condition
!
      endif
    endif
  enddo
!***************************************
!  Twobody VDW & Coulomb contribution  *
!***************************************
  t3 = g_cpu_time()
  if (lreaxFFQ) then
    eself = 0.0_dp
    if (literativeQ) then
      call setreaxffQiter1fc(nboatom,nboatomptr,nboatomRptr,nbos,nbosptr,qreaxFF,eself,.false., &
                             matom,vmatom)
    else
      call setreaxffQ1fc(nboatom,nboatomRptr,nbos,nbosptr,qreaxFF,eself,.false.,matom,vmatom)
    endif
  else
    eself = 0.0_dp
    qreaxFF(1:nboatom) = 0.0_dp
  endif
  t4 = g_cpu_time()
!
!  Set cutoffs
!
  cutvdw2 = reaxFFcutoffVDW**2
  cutq2   = reaxFFcutoffQ**2
  cut2    = max(cutvdw2,cutq2)
!
  xshift = vmatom(1) - xclat(matom)
  yshift = vmatom(2) - yclat(matom)
  zshift = vmatom(3) - zclat(matom)
!
!  Set lower bound for i loop
!
  if (ndim.gt.0) then
    imin = 1
  else
    imin = 2
  endif
!
!  Loop over pairs
!
  do ii = imin,nboatom
    i = nboatomRptr(ii)
!
!  Does i have a reaxFF species or a fixed charge?
!
    nspeci = nbosptr(i)
    lfixQi = lreaxFFqfix(nspeci)
    lreactivei = (nbos(i).gt.0.and..not.lfixQi)
!
    if (lreactivei.or.lfixQi) then
      gammai = reaxFFgamma(nspeci)
!
!  Set region for i
!
      nregioni = nregionno(nsft + nrelf2a(i))
      nregiontypi = nregiontype(nregioni,ncf)
!
!  Set upper bound for j loop
!
      if (ndim.gt.0) then
        jmax = ii
      else
        jmax = ii - 1
      endif
!
!  Loop over second atom
!
      jloop: do jj = 1,jmax
        j = nboatomRptr(jj)
!
!  Does j have a reaxFF species or a fixed charge?
!
        nspecj = nbosptr(j)
        lfixQj = lreaxFFqfix(nspecj)
        lreactivej = (nbos(j).gt.0.and..not.lfixQj)
!
        if (lreactivej.or.lfixQj) then
          gammaj = reaxFFgamma(nspecj)
!
!  Compute basic interatomic vector
!
          xji = xclat(j) - xclat(i)
          yji = yclat(j) - yclat(i)
          zji = zclat(j) - zclat(i)
!
!  Find valid vectors
!
          nor = 0
          if (ndim.eq.3) then
            call rsearch3D(xji,yji,zji,.false.,.false.,i,j,0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2)
          elseif (ndim.eq.2) then
            call rsearch2D(xji,yji,zji,.false.,.false.,i,j,0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2)
          elseif (ndim.eq.1) then
            call rsearch1D(xji,yji,zji,.false.,.false.,i,j,0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2)
          elseif (ndim.eq.0) then
            rij2 = xji*xji + yji*yji + zji*zji
            if (rij2.lt.cut2) then
              nor = nor + 1
            endif
          endif
!
!  If no distances then cycle
!
          if (nor.eq.0) cycle jloop
!
!  Correct distances for displacement
!
          if (ndim.gt.0) then
            if (i.eq.matom) then
              do n = 1,nor
                xtmp(n) = xtmp(n) - xshift
                ytmp(n) = ytmp(n) - yshift
                ztmp(n) = ztmp(n) - zshift
                dist(n) = xtmp(n)*xtmp(n) + ytmp(n)*ytmp(n) + ztmp(n)*ztmp(n)
              enddo
            elseif (j.eq.matom) then
              do n = 1,nor
                xtmp(n) = xtmp(n) + xshift
                ytmp(n) = ytmp(n) + yshift
                ztmp(n) = ztmp(n) + zshift
                dist(n) = xtmp(n)*xtmp(n) + ytmp(n)*ytmp(n) + ztmp(n)*ztmp(n)
              enddo
            endif
          else
            if (i.eq.matom) then
              xji = xji - xshift
              yji = yji - yshift
              zji = zji - zshift
            elseif (j.eq.matom) then
              xji = xji + xshift
              yji = yji + yshift
              zji = zji + zshift
            endif
            rij2 = xji*xji + yji*yji + zji*zji
          endif
!
          if (i.eq.j) then
            scale = 0.5_dp
          else
            scale = 1.0_dp
          endif
!
!  Set index for pairwise parameters
!
          if (nspeci.ge.nspecj) then
            ind = nspeci*(nspeci - 1)/2 + nspecj
          else
            ind = nspecj*(nspecj - 1)/2 + nspeci
          endif
!
!  Set region for j
!
          nregionj = nregionno(nsft + nrelf2a(j))
          nregiontypj = nregiontype(nregionj,ncf)
!
!  QM/MM handling : i & j are both QM atoms => exclude
!
          lQMMMok = .true.
          if (QMMMmode(ncf).gt.0) then
            if (nregiontypi.eq.1.and.nregiontypj.eq.1) lQMMMok = .false.
          endif
          if (lQMMMok) then
!
!  QM/MM : Set electrostatic embedding flag : If either i or j are QM atoms => exclude electrostatics
!
            lQMMMelectro = (QMMMmode(ncf).eq.2.and.(nregiontypi.eq.1.or.nregiontypj.eq.1))
!
            qij = scale*qreaxFF(i)*qreaxFF(j)*angstoev
!
!  Loop over valid distances and calculate contributions
!
            do n = 1,nor
              if (ndim.gt.0) rij2 = dist(n)
              rij = sqrt(rij2)
!
!  Find cell index
!
              if (abs(cellindex(1,n)).gt.nd2cell(1).or. &
                  abs(cellindex(2,n)).gt.nd2cell(2).or. &
                  abs(cellindex(3,n)).gt.nd2cell(3)) then
!
!  Find cell index - if outside user range then assign to central cell
!
                ncindp = nd2central
                ncindm = nd2central
              else
!
!  Compute index
!
                ncindp = nd2cellptr(nd2cell(1)+1+cellindex(1,n),nd2cell(2)+1+cellindex(2,n),nd2cell(3)+1+cellindex(3,n))
                ncindm = nd2cellptr(nd2cell(1)+1-cellindex(1,n),nd2cell(2)+1-cellindex(2,n),nd2cell(3)+1-cellindex(3,n))
              endif
!
!  Compute taper function(s)
!
              call p7reaxFFvdwtaper(rij,reaxFFcutoffVDW,tp,dtpdr,d2tpdr2,.true.,.false.)
              if (reaxFFcutoffQ.ne.reaxFFcutoffVDW) then
                call p7reaxFFqtaper(rij,reaxFFcutoffQ,tpQ,dtpQdr,d2tpQdr2,.true.,.false.)
              else
                tpQ      = tp
                dtpQdr   = dtpdr
              endif
              if (lreactivei.and.lreactivej) then
!
!  Compute f13
!
                call reaxFF_f13(rij,reaxFFgammaVDW(ind),f13,df13drij,d2f13drij2,.true.,.false.)
!
!  Energy terms
!
                VDWtrm  = reaxFFalphaVDW(ind)*(1.0_dp - f13/reaxFFr0VDW(ind))
                expVDW1 = exp(0.5_dp*VDWtrm)
                expVDW2 = expVDW1*expVDW1
!
!  Add VDW energy
!
                etrm = scale*reaxFFDeVDW(ind)*(expVDW2 - 2.0_dp*expVDW1)
              endif
              if (.not.lQMMMelectro) then
!
!  Compute Coulomb shielding parameters
!
                gammaij  = sqrt(gammai*gammaj)
                gammaij  = 1.0_dp/(gammaij**3)
!                     
!  Pure real space tapered lreaxFF case
!                     
                gam = qij/(rij*rij2 + gammaij)**third
              endif
!
              if (lreactivei.and.lreactivej) then
!
!  Compute derivative of VDW energy
!
                dVDWtrmdf13 = - reaxFFalphaVDW(ind)/reaxFFr0VDW(ind)
                dVDWtrmdrij = dVDWtrmdf13*df13drij
                devdwdrij = scale*reaxFFDeVDW(ind)*(expVDW2 - expVDW1)*dVDWtrmdrij*tp + etrm*dtpdr
              else
                devdwdrij = 0.0_dp
              endif
              if (.not.lQMMMelectro) then
!
!  Compute derivative of Coulomb energy
!           
                dgamdrij = - (gam*rij)/(rij*rij2 + gammaij)
                decouldrij = gam*dtpQdr + dgamdrij*tpQ
              else
                decouldrij = 0.0_dp
              endif
!
!  Get Cartesian components of vector and derivative
!
              if (ndim.gt.0) then
                xji = xtmp(n)
                yji = ytmp(n)
                zji = ztmp(n)
              endif
              xd = xji*(devdwdrij + decouldrij)
              yd = yji*(devdwdrij + decouldrij)
              zd = zji*(devdwdrij + decouldrij)
!
!  Add derivatives to main arrays
!
              if (i.ne.matom) then
                d1cell(1,i,ncindm) = d1cell(1,i,ncindm) - xd
                d1cell(2,i,ncindm) = d1cell(2,i,ncindm) - yd
                d1cell(3,i,ncindm) = d1cell(3,i,ncindm) - zd
              endif
              if (j.ne.matom) then
                d1cell(1,j,ncindp) = d1cell(1,j,ncindp) + xd
                d1cell(2,j,ncindp) = d1cell(2,j,ncindp) + yd
                d1cell(3,j,ncindp) = d1cell(3,j,ncindp) + zd
              endif
!
!  End loop over distances
!
            enddo
          endif
        endif
!
!  End of loop over j
!
      enddo jloop
    endif
  enddo
!
!  Output of charge / bond orders
!
  if (lqbo) then
    call outqbo(1_i4,numat,qreaxFF,maxneigh,nneigh,neighno,BO)
  endif
!
!  Free local memory
!
  deallocate(d1BOj_pipi,stat=status)
  if (status/=0) call deallocate_error('reaxFF1fc','d1BOj_pipi')
  deallocate(d1BOj_pi,stat=status)
  if (status/=0) call deallocate_error('reaxFF1fc','d1BOj_pi')
  deallocate(d1BOj_s,stat=status)
  if (status/=0) call deallocate_error('reaxFF1fc','d1BOj_s')
  deallocate(d1BOj,stat=status)
  if (status/=0) call deallocate_error('reaxFF1fc','d1BOj')
  deallocate(d1BOi_pipi,stat=status)
  if (status/=0) call deallocate_error('reaxFF1fc','d1BOi_pipi')
  deallocate(d1BOi_pi,stat=status)
  if (status/=0) call deallocate_error('reaxFF1fc','d1BOi_pi')
  deallocate(d1BOi_s,stat=status)
  if (status/=0) call deallocate_error('reaxFF1fc','d1BOi_s')
  deallocate(d1BOi,stat=status)
  if (status/=0) call deallocate_error('reaxFF1fc','d1BOi')
  deallocate(d1j,stat=status)
  if (status/=0) call deallocate_error('reaxFF1fc','d1j')
  deallocate(d1i,stat=status)
  if (status/=0) call deallocate_error('reaxFF1fc','d1i')
  deallocate(devaldsbo_neigh,stat=status)
  if (status/=0) call deallocate_error('reaxFF1fc','devaldsbo_neigh')
  deallocate(detorsdBOpi_neigh,stat=status)
  if (status/=0) call deallocate_error('reaxFF1fc','detorsdBOpi_neigh')
  deallocate(detorsdBO_neigh,stat=status)
  if (status/=0) call deallocate_error('reaxFF1fc','detorsdBO_neigh')
  deallocate(devaldBO_neigh,stat=status)
  if (status/=0) call deallocate_error('reaxFF1fc','devaldBO_neigh')
  deallocate(debpendBO_neigh,stat=status)
  if (status/=0) call deallocate_error('reaxFF1fc','debpendBO_neigh')
  deallocate(deoverunderddeltaj_neigh,stat=status)
  if (status/=0) call deallocate_error('reaxFF1fc','deoverunderddeltaj_neigh')
  deallocate(d1BOp_pipi,stat=status)
  if (status/=0) call deallocate_error('reaxFF1fc','d1BOp_pipi')
  deallocate(d1BOp_pi,stat=status)
  if (status/=0) call deallocate_error('reaxFF1fc','d1BOp_pi')
  deallocate(d1BOp,stat=status)
  if (status/=0) call deallocate_error('reaxFF1fc','d1BOp')
  deallocate(deltalplp,stat=status)
  if (status/=0) call deallocate_error('reaxFF1fc','deltalplp')
  deallocate(deltalp,stat=status)
  if (status/=0) call deallocate_error('reaxFF1fc','deltalp')
  deallocate(deltap,stat=status)
  if (status/=0) call deallocate_error('reaxFF1fc','deltap')
  deallocate(delta,stat=status)
  if (status/=0) call deallocate_error('reaxFF1fc','delta')
  deallocate(BOp_pipi,stat=status)
  if (status/=0) call deallocate_error('reaxFF1fc','BOp_pipi')
  deallocate(BOp_pi,stat=status)
  if (status/=0) call deallocate_error('reaxFF1fc','BOp_pi')
  deallocate(BOp,stat=status)
  if (status/=0) call deallocate_error('reaxFF1fc','BOp')
  deallocate(BO_pipi,stat=status)
  if (status/=0) call deallocate_error('reaxFF1fc','BO_pipi')
  deallocate(BO_pi,stat=status)
  if (status/=0) call deallocate_error('reaxFF1fc','BO_pi')
  deallocate(BO,stat=status)
  if (status/=0) call deallocate_error('reaxFF1fc','BO')
  deallocate(zneigh,stat=status)
  if (status/=0) call deallocate_error('reaxFF1fc','zneigh')
  deallocate(yneigh,stat=status)
  if (status/=0) call deallocate_error('reaxFF1fc','yneigh')
  deallocate(xneigh,stat=status)
  if (status/=0) call deallocate_error('reaxFF1fc','xneigh')
  deallocate(rneigh,stat=status)
  if (status/=0) call deallocate_error('reaxFF1fc','rneigh')
  deallocate(ineigh,stat=status)
  if (status/=0) call deallocate_error('reaxFF1fc','ineigh')
  deallocate(neighnoRptr,stat=status)
  if (status/=0) call deallocate_error('reaxFF1fc','neighnoRptr')
  deallocate(neighno,stat=status)
  if (status/=0) call deallocate_error('reaxFF1fc','neighno')
  deallocate(sum0,stat=status)
  if (status/=0) call deallocate_error('reaxFF1fc','sum0')
  deallocate(sum,stat=status)
  if (status/=0) call deallocate_error('reaxFF1fc','sum')
  deallocate(sumsij,stat=status)
  if (status/=0) call deallocate_error('reaxFF1fc','sumsij')
  deallocate(nneigh,stat=status)
  if (status/=0) call deallocate_error('reaxFF1fc','nneigh')
  deallocate(latomdone,stat=status)
  if (status/=0) call deallocate_error('reaxFF1fc','latomdone')
  deallocate(nbosptr,stat=status)
  if (status/=0) call deallocate_error('reaxFF1fc','nbosptr')
  deallocate(nbos,stat=status)
  if (status/=0) call deallocate_error('reaxFF1fc','nbos')
  deallocate(numattodoRptr,stat=status)
  if (status/=0) call deallocate_error('reaxFF1fc','numattodoRptr')
  deallocate(numattodoptr,stat=status)
  if (status/=0) call deallocate_error('reaxFF1fc','numattodoptr')
  deallocate(nbodummyptr,stat=status)
  if (status/=0) call deallocate_error('reaxFF1fc','nbodummyptr')
  deallocate(nboatomRptr,stat=status)
  if (status/=0) call deallocate_error('reaxFF1fc','nboatomRptr')
  deallocate(nboatomptr,stat=status)
  if (status/=0) call deallocate_error('reaxFF1fc','nboatomptr')
!
  t2 = g_cpu_time()
  treaxFF = treaxFF + t2 - t1 - (t4 - t3)
  teem    = teem + t4 - t3
#ifdef TRACE
  call trace_out('reaxff1fc')
#endif
!
  return
  end
