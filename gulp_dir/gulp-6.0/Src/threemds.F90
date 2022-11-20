  subroutine threemds(ethb,esregion12,esregion2,eattach,lgrad1)
!
!  Subroutine for three-body energy and forces using spatial decomposition.
!  Modified for distance dependent three body terms.
!
!  Strategy - sift by potential first, then cutoffs
!
!  lsymijk = if .true. then potential is symmetric w.r.t. i,j and k
!            and is not of angle centred form
!
!   5/03 Created from threemd
!   6/03 Global sums removed
!   7/03 Rstdl removed
!  10/03 Modifications for new spatial algorithm added
!   6/04 Sign of virial corrected
!   6/04 Virial now added on to total value
!   9/04 New arguments added to threebody
!  11/04 Bug fixed - esregion2l added to total
!  10/05 Hydrogen-bond potential added
!  12/05 ESFF equatorial potential added
!   9/06 Theta tapering added
!   1/07 UFF3 potential added
!   2/07 Bonding types added
!   5/07 QMMM schemes added
!   5/07 Dreiding option added
!   6/07 Structure of arrays for storing spatial distribution changed to 1-D
!   6/07 Dreiding option bonded2donorJK check added
!   7/07 Checking of bond orders added 
!  10/07 Error in checking of exocyclic attribute for bonds corrected
!  12/07 Unused variables removed
!   4/08 Modified for variable domain size
!   4/08 xvec1cell replaced by xvec2cell etc for spatial algorithm
!   4/08 ind1toijk replaced by ind2toijk
!   5/08 Handling of asymmetric bond orders corrected
!   6/08 Checking of bond numbers added
!  11/08 BAcoscross form added
!  11/08 Option to output energy terms added
!   3/08 3coulomb potential added
!   6/09 Site energy and virials added
!   6/09 Module name changed from three to m_three
!   7/09 Modifications for exp2 potential added
!  11/09 Region derivatives added
!   5/10 g3coulomb potential added
!  11/11 Region-region energy contributions stored
!   4/12 Explicit virial calculation removed as no longer needed
!   4/12 xvir, yvir and zvir removed
!   5/12 Atomic stresses added
!   2/13 BAGcross potential added
!   6/13 BALcross potential added
!   1/15 k3 and k4 terms added for cosine form
!   1/15 Modified Stillinger-Weber 3-body added
!   2/15 Correction added to include lswapk in initial setting of k values
!   3/15 MM3angle added
!   8/15 Garofalini form of sw3 added
!   1/18 Cell index shifts handled to allow for nomod
!   2/18 Trace added
!   8/18 Bug in swapping of k values corrected
!   9/18 Strain module added
!   9/18 Call to threestrterms modified to allow for second derivatives
!   9/18 Call to threestrterms replaced with more general realstrterms
!   6/19 j3 potential added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!  11/19 ppp3body potential added
!   2/20 Correction to cell index handling for bonding
!   3/20 lsymijk for ppp3body corrected
!   4/20 Rigid molecule modifications added
!   7/20 Setting of third com vector changed
!   7/20 lneedmol set to true for rigid molecule case
!   8/20 if statement for not cell buffer moved to avoid unnecessary calculation
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
  use configurations, only : lsliceatom, nregions, nregionno, nregiontype, QMMMmode
  use g_constants
  use control,        only : lmarvreg2, lseok, lDoQDeriv1, latomicstress, lrigid
  use current
  use derivatives
  use energies,       only : siteenergy, eregion2region
  use iochannels,     only : ioout
  use m_strain,       only : realstrterms
  use m_three
  use mdlogic
  use molecule
  use numbers,        only : third
  use optimisation
  use parallel
  use species,        only : spinspec
  use spatial
  use symmetry
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  real(dp),   intent(inout)                 :: ethb
  real(dp),   intent(inout)                 :: esregion12
  real(dp),   intent(inout)                 :: esregion2 
  real(dp),   intent(inout)                 :: eattach
  logical,    intent(in)                    :: lgrad1
!
!  Local variables
!
  integer(i4)                               :: i
  integer(i4)                               :: ic
  integer(i4)                               :: icx
  integer(i4)                               :: icy
  integer(i4)                               :: icz
  integer(i4)                               :: ii
  integer(i4)                               :: ijcx
  integer(i4)                               :: ijcy
  integer(i4)                               :: ijcz
  integer(i4)                               :: ikcx
  integer(i4)                               :: ikcy
  integer(i4)                               :: ikcz
  integer(i4)                               :: imx
  integer(i4)                               :: imy
  integer(i4)                               :: imz
  integer(i4)                               :: in3
  integer(i4)                               :: ind
  integer(i4)                               :: ind2
  integer(i4)                               :: indm
  integer(i4)                               :: indmj
  integer(i4)                               :: indmk
  integer(i4)                               :: indn
  integer(i4)                               :: ix
  integer(i4)                               :: iy
  integer(i4)                               :: iz
  integer(i4)                               :: ixi
  integer(i4)                               :: iyi
  integer(i4)                               :: izi
  integer(i4)                               :: ixj
  integer(i4)                               :: iyj
  integer(i4)                               :: izj
  integer(i4)                               :: ixk
  integer(i4)                               :: iyk
  integer(i4)                               :: izk
  integer(i4)                               :: ixyz
  integer(i4)                               :: j
  integer(i4)                               :: jc
  integer(i4)                               :: jcx
  integer(i4)                               :: jcy
  integer(i4)                               :: jcz
  integer(i4)                               :: jj
  integer(i4)                               :: jmx
  integer(i4)                               :: jmy
  integer(i4)                               :: jmz
  integer(i4)                               :: jndn
  integer(i4)                               :: k
  integer(i4)                               :: kc
  integer(i4)                               :: kcx
  integer(i4)                               :: kcy
  integer(i4)                               :: kcz
  integer(i4)                               :: kk
  integer(i4)                               :: kl
  integer(i4)                               :: ks
  integer(i4)                               :: maxx
  integer(i4)                               :: maxxy
  integer(i4)                               :: n
  integer(i4)                               :: n1i
  integer(i4)                               :: n1j
  integer(i4)                               :: n1k
  integer(i4)                               :: n3ty
  integer(i4)                               :: nati
  integer(i4)                               :: natj
  integer(i4)                               :: natk
  integer(i4)                               :: nbotyp11
  integer(i4)                               :: nbotyp12
  integer(i4)                               :: nbotyp21
  integer(i4)                               :: nbotyp22
  integer(i4)                               :: nbtyp11
  integer(i4)                               :: nbtyp12
  integer(i4)                               :: nbtyp21
  integer(i4)                               :: nbtyp22
  integer(i4)                               :: nbtypeij
  integer(i4)                               :: nbtypeik
  integer(i4)                               :: nbtypejk
  integer(i4)                               :: nbtypeij2
  integer(i4)                               :: nbtypeik2
  integer(i4)                               :: nbtypejk2
  integer(i4)                               :: neq
  integer(i4)                               :: ni
  integer(i4)                               :: nj
  integer(i4)                               :: nk
  integer(i4)                               :: nmi
  integer(i4)                               :: nmj
  integer(i4)                               :: nmk
  integer(i4)                               :: nr
  integer(i4)                               :: nregioni
  integer(i4)                               :: nregionj
  integer(i4)                               :: nregionk
  integer(i4)                               :: nregiontypi
  integer(i4)                               :: nregiontypj
  integer(i4)                               :: nregiontypk
  integer(i4)                               :: nsame
  integer(i4)                               :: nspeci
  integer(i4)                               :: nspecj
  integer(i4)                               :: nspeck
  integer(i4)                               :: nsplower(3)
  integer(i4)                               :: nspupper(3)
  integer(i4)                               :: nt1
  integer(i4)                               :: nt2
  integer(i4)                               :: nt3
  integer(i4)                               :: ntmp
  integer(i4)                               :: nto
  integer(i4)                               :: ntyp1
  integer(i4)                               :: ntyp2
  integer(i4)                               :: ntyp3
  integer(i4)                               :: ntypi
  integer(i4)                               :: ntypj
  integer(i4)                               :: ntypk
  integer(i4)                               :: ntypo
  integer(i4)                               :: status
  logical                                   :: bonded2donor
  logical                                   :: bonded2donorJK
  logical                                   :: l2bonds
  logical                                   :: lattach
  logical                                   :: lbonded
  logical                                   :: lbondnoOK
  logical                                   :: lbondtypeOK
  logical                                   :: lbtyp
  logical                                   :: ldiff23typ
  logical                                   :: ldiff23cut
  logical                                   :: ldiff23bo
  logical                                   :: linter_only
  logical                                   :: lintra_only
  logical                                   :: lmatch
  logical                                   :: lmolok
  logical                                   :: lmolok2
  logical                                   :: lneedmol
  logical                                   :: lopi
  logical                                   :: lopj
  logical                                   :: lopk
  logical                                   :: lreg12
  logical                                   :: lreg2trio
  logical                                   :: lsamemol
  logical                                   :: lsg1
  logical                                   :: lslicei
  logical                                   :: lslicej
  logical                                   :: lslicek
  logical                                   :: lsymijk
  logical                                   :: lswaprho
  logical                                   :: lswapk
  real(dp)                                  :: ang
  real(dp)                                  :: g_cpu_time
  real(dp)                                  :: d0i
  real(dp)                                  :: d0j
  real(dp)                                  :: d0k
  real(dp)                                  :: d1q(3,3)
  real(dp)                                  :: d2q(6)
  real(dp)                                  :: dot
  real(dp)                                  :: e2d(6)
  real(dp)                                  :: e3d(1)
  real(dp)                                  :: eattachl
  real(dp)                                  :: ed11
  real(dp)                                  :: ed12
  real(dp)                                  :: ed13
  real(dp)                                  :: esregion12l
  real(dp)                                  :: esregion2l
  real(dp)                                  :: ethb1
  real(dp)                                  :: ofct
  real(dp)                                  :: oci
  real(dp)                                  :: ocj
  real(dp)                                  :: ock
  real(dp)                                  :: one
  real(dp)                                  :: qli
  real(dp)                                  :: qlj
  real(dp)                                  :: qlk
  real(dp)                                  :: dr2ds(6,3)
  real(dp)                                  :: d2r2dx2(3,3,3)
  real(dp)                                  :: d2r2ds2(6,6,3)
  real(dp)                                  :: d2r2dsdx(6,3,3)
  real(dp)                                  :: r12
  real(dp)                                  :: r122
  real(dp)                                  :: r13
  real(dp)                                  :: r132
  real(dp)                                  :: r23
  real(dp)                                  :: r232
  real(dp)                                  :: rho1
  real(dp)                                  :: rho2
  real(dp)                                  :: rho3
  real(dp)                                  :: rho4
  real(dp)                                  :: rho5
  real(dp)                                  :: rk32
  real(dp)                                  :: rk33
  real(dp)                                  :: rk34
  real(dp)                                  :: rkt2
  real(dp)                                  :: rkt3
  real(dp)                                  :: rkt4
  real(dp)                                  :: rkthb
  real(dp)                                  :: rkthb3
  real(dp)                                  :: rkthb4
  real(dp)                                  :: rktmp
  real(dp)                                  :: rstrdloc(6)
  real(dp)                                  :: ro1
  real(dp)                                  :: ro2
  real(dp)                                  :: ro3
  real(dp)                                  :: ro4
  real(dp)                                  :: ro5
  real(dp)                                  :: rro
  real(dp)                                  :: symfct
  real(dp)                                  :: the0
  real(dp)                                  :: time1
  real(dp)                                  :: time2
  real(dp)                                  :: tr1
  real(dp)                                  :: tr11
  real(dp)                                  :: tr11m
  real(dp)                                  :: tr1m
  real(dp)                                  :: tr2
  real(dp)                                  :: tr21
  real(dp)                                  :: tr21m
  real(dp)                                  :: tr2m
  real(dp)                                  :: tr3
  real(dp)                                  :: tr31
  real(dp)                                  :: tr31m
  real(dp)                                  :: tr3m
  real(dp)                                  :: ttmp
  real(dp)                                  :: ttr1
  real(dp)                                  :: ttr1m
  real(dp)                                  :: ttr11
  real(dp)                                  :: ttr2
  real(dp)                                  :: ttr2m
  real(dp)                                  :: ttr21
  real(dp)                                  :: x23
  real(dp)                                  :: y23
  real(dp)                                  :: z23
  real(dp)                                  :: x21
  real(dp)                                  :: y21
  real(dp)                                  :: z21
  real(dp)                                  :: x31
  real(dp)                                  :: y31
  real(dp)                                  :: z31
  real(dp)                                  :: xcom(3)
  real(dp)                                  :: ycom(3)
  real(dp)                                  :: zcom(3)
  real(dp)                                  :: xcomi
  real(dp)                                  :: ycomi
  real(dp)                                  :: zcomi
  real(dp), dimension(:), allocatable, save :: xderv
  real(dp), dimension(:), allocatable, save :: yderv
  real(dp), dimension(:), allocatable, save :: zderv
  real(dp)                                  :: xi
  real(dp)                                  :: yi
  real(dp)                                  :: zi
  real(dp)                                  :: xv3(3)
  real(dp)                                  :: yv3(3)
  real(dp)                                  :: zv3(3)
#ifdef TRACE
  call trace_in('threemds')
#endif
!
  time1 = g_cpu_time()
  lsg1 = (lstr.and.lgrad1)
!
!  Opening banner for energy decomposition
!
  if (lPrintThree) then
    call mpbarrier
    if (ioproc) then
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''  Three: Atom No. 1  Atom No. 2  Atom No. 3              Threebody energy (eV)  '')')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
    endif
  endif
!
!  Initialisation
!
  ethb = 0.0_dp
  esregion12l = 0.0_dp
  esregion2l = 0.0_dp
  eattachl = 0.0_dp
  one = 1.0_dp
!
!  Allocate local memory
!
  allocate(xderv(numat),stat=status)
  if (status/=0) call outofmemory('threemds','xderv')
  allocate(yderv(numat),stat=status)
  if (status/=0) call outofmemory('threemds','yderv')
  allocate(zderv(numat),stat=status)
  if (status/=0) call outofmemory('threemds','zderv')
!
  if (lgrad1) then
    do i = 1,numat
      xderv(i) = 0.0_dp
      yderv(i) = 0.0_dp
      zderv(i) = 0.0_dp
    enddo
  endif
!
!  Set variables for block distribution
!
  maxxy = nspcell(1)*nspcell(2)
  maxx  = nspcell(1)
!
!  Loop over potentials
!
  potentials: do n = 1,nthb
    n3ty = nthrty(n)
    nt1 = ntspec1(n)
    nt2 = ntspec2(n)
    nt3 = ntspec3(n)
    ntyp1 = ntptyp1(n)
    ntyp2 = ntptyp2(n)
    ntyp3 = ntptyp3(n)
    nbtyp11 = n3botype(1,1,n)
    nbtyp12 = n3botype(2,1,n)
    nbtyp21 = n3botype(1,2,n)
    nbtyp22 = n3botype(2,2,n)
    lbtyp = (mmtexc(n).eq.1)
    tr11m = thr1min(n)
    tr21m = thr2min(n)
    tr31m = thr3min(n)
    tr11 = thr1(n)
    tr21 = thr2(n)
    tr31 = thr3(n)
!
!  Check that atomic numbers match
!
    ldiff23typ = (ntyp2.eq.ntyp3.and.nt2.eq.nt3)
!
!  ..and the cutoffs..
!
    ldiff23cut = (tr11.eq.tr21.and.tr11m.eq.tr21m)
!
!  ..and the bond orders
!
    ldiff23bo = (nbtyp11.ne.nbtyp21.or.nbtyp12.ne.nbtyp22)
!
    tr1m = tr11m*tr11m
    tr2m = tr21m*tr21m
    tr3m = tr31m*tr31m
    tr1 = tr11*tr11
    tr2 = tr21*tr21
    tr3 = tr31*tr31
    rkt2 = thbk(n)
    rkt3 = 0.0_dp
    rkt4 = 0.0_dp
    ro1 = 0.0_dp
    ro2 = 0.0_dp
    ro3 = 0.0_dp
    ro4 = 0.0_dp
    ro5 = 0.0_dp
    if (n3ty.eq.2) then
      the0 = theta(n)*degtorad
      ro1 = thrho1(n)
      ro2 = thrho2(n)
      if (ro1.ne.0.0_dp) ro1 = 1.0_dp/ro1
      if (ro2.ne.0.0_dp) ro2 = 1.0_dp/ro2
      lsymijk = .false.
      lswaprho = .true.
      lswapk = .false.
    elseif (n3ty.eq.1) then
      the0 = theta(n)*degtorad
      rkt3 = thrho2(n)/6.0_dp
      rkt4 = thrho1(n)/24.0_dp
      lsymijk = .false.
      lswaprho = .false.
      lswapk = .false.
    elseif (n3ty.eq.3) then
      the0 = 0.0_dp
      lsymijk = .true.
      lswaprho = .false.
      lswapk = .false.
    elseif (n3ty.eq.4) then
      the0 = 0.0_dp
      ro1 = theta(n)
      ro2 = thrho1(n)
      ro3 = thrho2(n)
      lsymijk = .true.
      lswaprho = .true.
      lswapk = .false.
    elseif (n3ty.eq.5.or.n3ty.eq.25) then
      the0 = cos(theta(n)*degtorad)
      ro1 = thrho1(n)
      ro2 = thrho2(n)
      ro3 = thrho3(n)
      lsymijk = .false.
      lswaprho = .true.
      lswapk = .false.
    elseif (n3ty.eq.6) then
      the0 = 0.0_dp
      ro1 = thrho1(n)
      ro2 = thrho2(n)
      lsymijk = .false.
      lswaprho = .true.
      lswapk = .false.
    elseif (n3ty.eq.7) then
      the0 = theta(n)
      lsymijk = .false.
      lswaprho = .false.
      lswapk = .false.
    elseif (n3ty.eq.8) then
      the0 = theta(n)*degtorad
      ro1 = thrho1(n)
      ro2 = thrho2(n)
      if (ro1.ne.0.0_dp) ro1 = 1.0_dp/ro1
      if (ro2.ne.0.0_dp) ro2 = 1.0_dp/ro2
      the0 = the0 - pi
      the0 = the0*the0
      rkt2 = 0.25_dp*rkt2/the0
      lsymijk = .false.
      lswaprho = .true.
      lswapk = .false.
    elseif (n3ty.eq.9) then
      the0 = cos(theta(n)*degtorad)
      rkt3 = thrho2(n)/6.0_dp
      rkt4 = thrho1(n)/24.0_dp
      lsymijk = .false.
      lswaprho = .false.
      lswapk = .false.
    elseif (n3ty.eq.10) then
      the0 = theta(n)
      ro1 = thrho1(n)
      ro2 = thrho2(n)
      ro3 = thrho3(n)
      lsymijk = .true.
      lswaprho = .false.
      lswapk = .false.
    elseif (n3ty.eq.11) then
      rkt3 = thrho3(n)
      ro1 = thrho1(n)
      ro2 = thrho2(n)
      the0 = theta(n)*degtorad
      lsymijk = .false.
      lswaprho = .true.
      lswapk = .true.
    elseif (n3ty.eq.12) then
      the0 = theta(n)
      rkt3 = nint(thrho1(n))
      lsymijk = .false.
      lswaprho = .false.
      lswapk = .false.
    elseif (n3ty.eq.13) then
      the0 = theta(n)
      rkt3 = thrho3(n)
      ro1 = thrho1(n)
      ro2 = thrho2(n)  
      lsymijk = .false.
      lswaprho = .true.
      lswapk = .false.
    elseif (n3ty.eq.14) then
      the0 = cos(theta(n)*degtorad)
      ro1 = thrho1(n)
      ro2 = thrho2(n)
      ro3 = thrho3(n)
      lsymijk = .false.
      lswaprho = .true.
      lswapk = .false.
    elseif (n3ty.eq.15) then
      rkt3 = theta(n)
      ro1 = thrho1(n)
      ro2 = thrho2(n)
      ro3 = thrho3(n)
      lsymijk = .false.
      lswaprho = .true.
      lswapk = .false.
    elseif (n3ty.eq.16) then
      the0 = theta(n)
      ro1 = thrho1(n)
      ro2 = thrho2(n)
      lsymijk = .false.
      lswaprho = .true.
      lswapk = .false.
    elseif (n3ty.eq.17) then
      the0 = theta(n)*degtorad
      rkt4 = 1.0_dp/(2.0_dp*sin(the0))**2
      rkt3 = - 4.0_dp*rkt4*cos(the0)
      the0 = rkt4*(2.0_dp*cos(the0)**2 + 1.0_dp)
      lsymijk = .false.
      lswaprho = .false.
      lswapk = .false.
    elseif (n3ty.eq.18) then
      rkt3 = thrho3(n)
      ro1 = thrho1(n)
      ro2 = thrho2(n)
      the0 = cos(theta(n)*degtorad)
      lsymijk = .false.
      lswaprho = .true.
      lswapk = .true.
    elseif (n3ty.eq.19) then
      lsymijk = .false.
      lswaprho = .false.
      lswapk = .false.
    elseif (n3ty.eq.20) then
      the0 = 0.0_dp
      ro1 = theta(n)
      ro2 = thrho2(n)
      ro4 = thrho1(n)
      ro5 = thrho3(n)
      lsymijk = .false.
      lswaprho = .true.
      lswapk = .false.
    elseif (n3ty.eq.21) then
      the0 = theta(n)
      lsymijk = .false.
      lswaprho = .false.
    elseif (n3ty.eq.22) then
      ro1 = thrho1(n)
      ro2 = thrho2(n)
      the0 = theta(n)*degtorad
      lsymijk = .false.
      lswaprho = .true.
      lswapk = .true.
    elseif (n3ty.eq.23) then
      ro1 = thrho1(n)
      ro2 = thrho2(n)
      the0 = theta(n)
      lsymijk = .false.
      lswaprho = .true.
      lswapk = .true.
    elseif (n3ty.eq.24) then
      the0 = theta(n)*degtorad
      lsymijk = .false.
      lswaprho = .true.
      lswapk = .false.
    elseif (n3ty.eq.26) then
      the0 = 0.0_dp
      lsymijk = .true.
      lswaprho = .false.
      lswapk = .false.
    elseif (n3ty.eq.27) then
      the0 = 0.0_dp
      ro1 = thrho1(n)
      ro2 = thrho2(n)
      ro3 = thrho3(n)
      lsymijk = .false.
      lswaprho = .true.
      lswapk = .false.
    endif
    lintra_only = (ltintra(n).and..not.ltinter(n))
    linter_only = (ltinter(n).and..not.ltintra(n))
    lneedmol = (lintra_only.or.linter_only.or.lbtyp.or.lrigid)
!
!  Work out symmetry factor for symmetric potentials
!
    if (lsymijk) then
      nsame = 0
      if (nt1.eq.nt2.and.ntyp1.eq.ntyp2) nsame = nsame + 1
      if (nt1.eq.nt3.and.ntyp1.eq.ntyp3) nsame = nsame + 1
      if (nsame.eq.0) then
        symfct = 1.0_dp
      elseif (nsame.eq.1) then
        symfct = 0.5_dp
      elseif (nsame.eq.2) then
        symfct = 1.0_dp/3.0_dp
      endif
    endif
!
!  Loop over all local spatial cells except buffer regions
!
    do ixyz = 1,ncellpernode
      if (.not.lbuffercell(ixyz)) then
        ind = ncellnodeptr(ixyz)
        ind2 = ind - 1
        iz = ind2/maxxy
        ind2 = ind2 - maxxy*iz
        iy = ind2/maxx
        ix = ind2 - maxx*iy + 1
        iy = iy + 1
        iz = iz + 1
!
!  Set cell search bounds
!
        nspupper(1) = min(ix+ncellsearch(1),nspcell(1))
        nspupper(2) = min(iy+ncellsearch(2),nspcell(2))
        nspupper(3) = min(iz+ncellsearch(3),nspcell(3))
        nsplower(1) = max(ix-ncellsearch(1),1)
        nsplower(2) = max(iy-ncellsearch(2),1)
        nsplower(3) = max(iz-ncellsearch(3),1)
!
!  Get number of atoms in this cell
!  
        ni = nspcellat(ind)
        n1i = nspcellat1ptr(ind)
!
!  Outer loop over atoms within this cell
!
        iloop: do ii = 1,ni
          i = nspcellatptr(n1i+ii)
!
!  Check i is allowed for n
!
          nati = nat(i)
          ntypi = nftype(i)   
          nspeci = nspecptr(nrelf2a(i))
          if (.not.lmatch(nati,ntypi,nt1,ntyp1,.true.)) cycle iloop
!
!  Set properties of atom i
!  
          ic = nspcellatptrcell(n1i+ii)
          call ind2toijk(ic,icx,icy,icz)
          oci = occuf(i)
          qli = qf(i)
          nregioni = nregionno(nsft+nrelf2a(i))
          nregiontypi = nregiontype(nregioni,ncf)
          lopi = (.not.lfreeze.or.lopf(nrelf2a(i)))
          lslicei = lsliceatom(nsft + nrelf2a(i))
!
!  QM/MM handling
!
          if (QMMMmode(ncf).gt.0) then
            if (nregiontypi.eq.1.and.lbtyp) cycle iloop
          endif
! 
!  Dreiding option handling                 
! 
          if (ltdreiding(n)) then               
            if (.not.bonded2donor(i)) cycle iloop
          endif
!
!  Set coordinates of atom i
!
          xi = xinbox(i) + xvec2cell(ic)
          yi = yinbox(i) + yvec2cell(ic)
          zi = zinbox(i) + zvec2cell(ic)
!
!  Molecule handling
!
          if (lmol.and.lneedmol) then
            nmi = natmol(i)
            if (ndim.gt.0) then
              indm = nmolind(i)
              call mindtoijk(indm,ixi,iyi,izi)
            endif
!
!  Set COM coordinates
!
            if (lrigid.and.nmi.gt.0) then
              xcomi = molxyz(1,natinmol(i),nmi)
              ycomi = molxyz(2,natinmol(i),nmi)
              zcomi = molxyz(3,natinmol(i),nmi)
            else
              xcomi = 0.0_dp
              ycomi = 0.0_dp
              zcomi = 0.0_dp
            endif
          else
            xcomi = 0.0_dp
            ycomi = 0.0_dp
            zcomi = 0.0_dp
          endif
!     
!  Check number of bonds if necessary
!  
          if (n3bondnono(1,n).gt.0) then
            lbondnoOK = .false.
            do in3 = 1,n3bondnono(1,n)
              if (nbonds(i).eq.n3bondno(in3,1,n)) lbondnoOK = .true.
            enddo
            if (.not.lbondnoOK) cycle iloop
          endif
          if (n3bondnono(2,n).gt.0) then
            lbondnoOK = .true.
            do in3 = 1,n3bondnono(2,n)
              if (nbonds(i).eq.n3bondno(in3,2,n)) lbondnoOK = .false.
            enddo
            if (.not.lbondnoOK) cycle iloop
          endif
!
!  Loop over neighbouring cells
!  
          do imz = nsplower(3),nspupper(3)
            do imy = nsplower(2),nspupper(2)
              do imx = nsplower(1),nspupper(1)
                indn = (imz-1)*maxxy + (imy-1)*maxx + imx
!
!  Loop over atoms within neighbouring cells
!
                nj = nspcellat(indn)
                n1j = nspcellat1ptr(indn)
                jloop: do jj = 1,nj
                  j = nspcellatptr(n1j+jj)       
!
!  Set species type parameters for atom j
!
                  natj = nat(j)
                  ntypj = nftype(j)
                  nspecj = nspecptr(nrelf2a(j))
!
!  Check j is allowed for n
!
                  if (lmatch(natj,ntypj,nt2,ntyp2,.false.)) then
                    nto = nt3
                    ntypo = ntyp3
                    nbotyp11 = nbtyp11
                    nbotyp12 = nbtyp12
                    nbotyp21 = nbtyp21
                    nbotyp22 = nbtyp22
                    ttr1m = tr1m
                    ttr2m = tr2m
                    ttr1 = tr1
                    ttr2 = tr2
                    ttr11 = tr11
                    ttr21 = tr21
                    rho1 = ro1
                    rho2 = ro2
                    rho3 = ro3
                    rho4 = ro4
                    rho5 = ro5
                    rkthb  = rkt2
                    rkthb3 = rkt3
                    rkthb4 = rkt4
                  elseif (lmatch(natj,ntypj,nt3,ntyp3,.true.)) then
                    nto = nt2
                    ntypo = ntyp2
                    nbotyp11 = nbtyp21
                    nbotyp12 = nbtyp22
                    nbotyp21 = nbtyp11
                    nbotyp22 = nbtyp12
                    ttr1m = tr2m
                    ttr2m = tr1m
                    ttr1 = tr2
                    ttr2 = tr1
                    ttr11 = tr21
                    ttr21 = tr11
                    rho1 = ro2
                    rho2 = ro1
                    rho3 = ro3
                    rho4 = ro5
                    rho5 = ro4
                    if (lswapk) then
                      rkthb  = rkt3
                      rkthb3 = rkt2
                    else
                      rkthb  = rkt2
                      rkthb3 = rkt3
                    endif
                    rkthb4 = rkt4
                  else
                    cycle jloop
                  endif
                  if (lmol.and.lneedmol) then
!
!  Molecule handling
!
                    nmj = natmol(j)
                    if (ndim.gt.0) then
                      indmj = nmolind(j)
                      call mindtoijk(indmj,ixj,iyj,izj)
                      ixj = ixj - ixi
                      iyj = iyj - iyi
                      izj = izj - izi
                    endif
                    lmolok = (nmi.eq.nmj.and.nmi.gt.0)
!
!  Set COM coordinates
!
                    if (lrigid.and.nmj.gt.0) then
                      xcom(1) = molxyz(1,natinmol(j),nmj) - xcomi
                      ycom(1) = molxyz(2,natinmol(j),nmj) - ycomi
                      zcom(1) = molxyz(3,natinmol(j),nmj) - zcomi
                    else
                      xcom(1) = - xcomi
                      ycom(1) = - ycomi
                      zcom(1) = - zcomi
                    endif
                  else
                    lmolok = .false.
                    xcom(1) = - xcomi
                    ycom(1) = - ycomi
                    zcom(1) = - zcomi
                  endif
!
!  Check for intra and but not in same molecule
!
                  if (lintra_only.and..not.lmolok) cycle jloop
                  if (lbtyp.and..not.lmolok) cycle jloop
!
!  Set remaining properties for j
!
                  jc = nspcellatptrcell(n1j+jj)
                  call ind2toijk(jc,jcx,jcy,jcz)
                  nregionj = nregionno(nsft+nrelf2a(j))
                  nregiontypj = nregiontype(nregionj,ncf)
                  ocj = occuf(j)
                  qlj = qf(j)
                  lopj = (lopf(nrelf2a(j)).or..not.lfreeze)
                  lslicej = lsliceatom(nsft + nrelf2a(j))
!
!  Calculate vector from atom 1 to atom 2
!
                  x21 = xvec2cell(jc) + xinbox(j) - xi
                  y21 = yvec2cell(jc) + yinbox(j) - yi
                  z21 = zvec2cell(jc) + zinbox(j) - zi
!
!  Calculate and check r12 is OK
!
                  r122 = x21*x21 + y21*y21 + z21*z21
                  if (r122.lt.1.0d-10) cycle jloop
!
!  Molecule checking
!
                  lbonded = .false.
                  if (lmolok) then
                    if (ndim.eq.0) then
                      if (linter_only) cycle jloop
                      if (lbtyp) then
                        call bonded(lbonded,l2bonds,nbtypeij,nbtypeij2,i,j,0_i4,0_i4,0_i4)
                        if (.not.lbonded) cycle jloop
                      endif
                    else
                      ijcx = jcx - icx + icosxsp(j) - icosxsp(i)
                      ijcy = jcy - icy + icosysp(j) - icosysp(i)
                      ijcz = jcz - icz + icoszsp(j) - icoszsp(i)
!
                      if (lbtyp) then
                        call bonded(lbonded,l2bonds,nbtypeij,nbtypeij2,i,j,ijcx,ijcy,ijcz)
                        if (.not.lbonded) cycle jloop
                        lsamemol = (lbonded.or.l2bonds)
                      else
                        lsamemol = .false.
                      endif
                      if (.not.lsamemol) then
                        call samemol(lsamemol,nmi,ijcx,ijcy,ijcz,ixj,iyj,izj)
                      endif
                      if (lintra_only.and..not.lsamemol) cycle jloop
                      if (linter_only.and.lsamemol) cycle jloop
                    endif
                  endif
!
!  Distance checking
!
                  if ((r122.gt.ttr1.or.r122.lt.ttr1m).and.(.not.lbtyp.or..not.lbonded)) then
                    if (ldiff23typ.or..not.ldiff23cut) cycle jloop
                    if (r122.lt.ttr2.and.r122.gt.ttr2m) then
                      ttmp = ttr2m
                      ttr2m = ttr1m
                      ttr1m = ttmp
                      ttmp = ttr2
                      ttr2 = ttr1
                      ttr1 = ttmp
                      if (lswaprho) then
                        rro = rho1
                        rho1 = rho2
                        rho2 = rro
                        rro = rho4
                        rho4 = rho5
                        rho5 = rro
                        if (lswapk) then
                          rktmp = rkthb
                          rkthb = rkthb3
                          rkthb3 = rktmp
                        endif
                      endif
                    else
                      cycle jloop
                    endif
                  endif
                  r12 = sqrt(r122)
!
!  Loop over neighbouring cells for k
!
                  do jmz = nsplower(3),nspupper(3)
                    do jmy = nsplower(2),nspupper(2)
                      do jmx = nsplower(1),nspupper(1)
                        jndn = (jmz-1)*maxxy + (jmy-1)*maxx + jmx
!
!  Loop over atoms within neighbouring cells
!
                        nk = nspcellat(jndn)
                        n1k = nspcellat1ptr(jndn)
                        kloop: do kk = 1,nk
                          k = nspcellatptr(n1k+kk)       
!
!  Only calculate lower-half triangular interactions
!
                          if (k.gt.j) then
                            natk = nat(k)
                            ntypk = nftype(k)
                            nspeck = nspecptr(nrelf2a(k))
!
!  Check k is allowed for n, and not equivalent to j
!
                            if (.not.lmatch(natk,ntypk,nto,ntypo,.true.)) cycle kloop
!
!  If lfreeze and all atoms are fixed then skip this
!  three body term
!
                            lopk = (lopf(nrelf2a(k)).or..not.lfreeze)
                            if (.not.lopi.and..not.lopj.and..not.lopk) cycle kloop
!
!  Dreiding option handling
!
                            if (ltdreiding(n)) then
                              if (.not.bonded2donorJK(i,j,k)) cycle kloop
                            endif
!
!  Set remaining properties of atom k
!
                            kc = nspcellatptrcell(n1k+kk)
                            call ind2toijk(kc,kcx,kcy,kcz)
                            nregionk = nregionno(nsft+nrelf2a(k))
                            nregiontypk = nregiontype(nregionk,ncf)
                            ock = occuf(k)
                            qlk = qf(k)
!  
!  QM/MM handling
!           
                            if (QMMMmode(ncf).gt.0) then
                              if (nregiontypi.eq.1.and.nregiontypj.eq.1.and.nregiontypk.eq.1) cycle kloop
                            endif
!
                            if (lmol.and.lneedmol) then
!
!  Molecule handling
!
                              nmk = natmol(k)
                              if (ndim.gt.0) then
                                indmk = nmolind(k)
                                call mindtoijk(indmk,ixk,iyk,izk)
                                ixk = ixk - ixi
                                iyk = iyk - iyi
                                izk = izk - izi
                              endif
                              lmolok2 = (nmi.eq.nmk.and.nmi.gt.0)
!
!  Set COM coordinates
!
                              if (lrigid.and.nmk.gt.0) then
                                xcom(2) = molxyz(1,natinmol(k),nmk) - xcomi
                                ycom(2) = molxyz(2,natinmol(k),nmk) - ycomi
                                zcom(2) = molxyz(3,natinmol(k),nmk) - zcomi
                              else
                                xcom(2) = - xcomi
                                ycom(2) = - ycomi
                                zcom(2) = - zcomi
                              endif
                            else
                              lmolok2 = .false.
                              xcom(2) = - xcomi
                              ycom(2) = - ycomi
                              zcom(2) = - zcomi
                            endif
                            xcom(3) = xcom(2) - xcom(1)
                            ycom(3) = ycom(2) - ycom(1)
                            zcom(3) = zcom(2) - zcom(1)
!
!  Check for intra and but not in same molecule
!
                            if (lintra_only.and..not.lmolok2) cycle kloop
                            if (lbtyp.and..not.lmolok2) cycle kloop
!
!  Calculate vector from i to k
!
                            x31 = xvec2cell(kc) + xinbox(k) - xi
                            y31 = yvec2cell(kc) + yinbox(k) - yi
                            z31 = zvec2cell(kc) + zinbox(k) - zi
!
!  Check r13 is OK
!
                            r132 = x31*x31 + y31*y31 + z31*z31
                            if (r132.lt.1.0d-10) cycle kloop
!
!  Molecule checking
!
                            lbonded = .false.
                            if (lmolok2) then
                              if (ndim.eq.0) then
                                if (linter_only) cycle kloop
                                if (lbtyp) then
                                  call bonded(lbonded,l2bonds,nbtypeik,nbtypeik2,i,k,0_i4,0_i4,0_i4)
                                  if (.not.lbonded) cycle kloop
                                  if (lsymijk) then
                                    call bonded(lbonded,l2bonds,nbtypejk,nbtypejk2,j,k,0_i4,0_i4,0_i4)
                                    if (.not.lbonded) cycle kloop
                                  endif
                                endif
                              else
                                ikcx = kcx - icx + icosxsp(k) - icosxsp(i)
                                ikcy = kcy - icy + icosysp(k) - icosysp(i)
                                ikcz = kcz - icz + icoszsp(k) - icoszsp(i)
!
                                if (lbtyp) then
                                  call bonded(lbonded,l2bonds,nbtypeik,nbtypeik2,i,k,ikcx,ikcy,ikcz)
                                  if (.not.lbonded) cycle kloop
                                  if (lsymijk) then
                                    call bonded(lbonded,l2bonds,nbtypejk,nbtypejk2,j,k,ikcx-ijcx,ikcy-ijcy,ikcz-ijcz)
                                    if (.not.lbonded) cycle kloop
                                  endif
                                  lsamemol = (lbonded.or.l2bonds)
                                else
                                  lsamemol = .false.
                                endif
                                if (.not.lsamemol) then
                                  call samemol(lsamemol,nmi,ikcx,ikcy,ikcz,ixk,iyk,izk)
                                endif
                                if (lintra_only.and..not.lsamemol) cycle kloop
                                if (linter_only.and.lsamemol) cycle kloop
                              endif
                            endif
!
!  Distance checking
!
!  Modification to handle case where species for 2 and 3 are the same
!  but cutoffs are different
!
                            if ((r132.gt.ttr2.or.r132.lt.ttr2m).and.(.not.lbtyp.or..not.lbonded)) then
                              if (ldiff23typ.or..not.ldiff23cut) cycle kloop
                              if (r122.gt.ttr2.or.r132.gt.ttr1) cycle kloop
                              if (r122.lt.ttr2m.or.r132.lt.ttr1m) cycle kloop
                              if (lswaprho) then
                                rro = rho1
                                rho1 = rho2
                                rho2 = rro
                                rro = rho4
                                rho4 = rho5
                                rho5 = rro
                                if (lswapk) then
                                  rktmp = rkthb
                                  rkthb = rkthb3
                                  rkthb3 = rktmp
                                endif
                              endif
                            endif
                            if (lbtyp) then
!
!  Bond type checking
!
!  If we've made it this far then the atoms must be bonded. We just have
!  to check if the bond orders match.
!
                              lbondtypeOK = .true.
!
!  Check i-j bond for correct order
!
                              if (nbotyp11.gt.0.and.nbotyp11.ne.nbtypeij) lbondtypeOK = .false.
                              if (nbotyp12.gt.1.and.nbotyp12.ne.nbtypeij2) lbondtypeOK = .false.
!
!  Check i-k bond for correct order
!
                              if (nbotyp21.gt.0.and.nbotyp21.ne.nbtypeik) lbondtypeOK = .false.
                              if (nbotyp22.gt.1.and.nbotyp22.ne.nbtypeik2) lbondtypeOK = .false.
!
                              if (.not.lbondtypeOK) then
!
!  If bond types don't match, but atom types are symmetric then try other permutation
!
                                if (.not.ldiff23typ.or..not.ldiff23bo) cycle kloop
!
!  Check i-j bond for correct order
! 
                                if (nbotyp21.gt.0.and.nbotyp21.ne.nbtypeij) cycle kloop
                                if (nbotyp22.gt.1.and.nbotyp22.ne.nbtypeij2) cycle kloop
!
!  Check i-k bond for correct order
! 
                                if (nbotyp11.gt.0.and.nbotyp11.ne.nbtypeik) cycle kloop
                                if (nbotyp12.gt.1.and.nbotyp12.ne.nbtypeik2) cycle kloop
!
!  If we make it to here then bond orders are the wrong way wrong and i-j/i-k terms should be swapped
!
                                if (lswaprho) then
                                  ntmp = nbotyp11
                                  nbotyp11 = nbotyp21
                                  nbotyp21 = ntmp
                                  ntmp = nbotyp12
                                  nbotyp12 = nbotyp22
                                  nbotyp22 = ntmp
                                  rro = rho1
                                  rho1 = rho2
                                  rho2 = rro
                                  rro = rho4
                                  rho4 = rho5
                                  rho5 = rro
                                  if (lswapk) then
                                    rktmp = rkthb
                                    rkthb = rkthb3
                                    rkthb3 = rktmp
                                  endif
                                endif
                              endif
                            endif
!
                            r13 = sqrt(r132)
!
!  Check r23 is OK
!
                            x23 = x31 - x21
                            y23 = y31 - y21
                            z23 = z31 - z21
                            r232 = x23**2 + y23**2 + z23**2
                            if (r232.gt.tr3.and..not.lbtyp) cycle kloop
                            if (r232.lt.tr3m.or.r232.lt.1d-10) cycle kloop
!
!  Valid three-body term => calculate potential
!
                            if (n3ty.ne.3.and.n3ty.ne.4.and.n3ty.ne.6.and.n3ty.ne.7.and.n3ty.ne.19.and.n3ty.ne.26 &
                                .and.n3ty.ne.27) then
                              dot = x21*x31 + y21*y31 + z21*z31
                              dot = dot/(r12*r13)
                              if (abs(dot).gt.0.999999999999_dp) dot = sign(one,dot)
                              if (n3ty.eq.9) then
                                ang = dot
                              else
                                ang = acos(dot)
                              endif
                            else
                              dot = 0.0_dp
                              ang = 0.0_dp
                            endif
                            r23 = sqrt(r232)
                            ofct = oci*ocj*ock
                            if (lsymijk) ofct = ofct*symfct
                            if (n3ty.eq.19) then
                              rk32 = rkthb*ofct*qf(j)*qf(k)
                            else
                              rk32 = rkthb*ofct
                            endif
                            if (n3ty.eq.12.or.n3ty.eq.17) then
                              rk33 = rkthb3
                            elseif (n3ty.eq.26) then
                              rk33 = spinspec(nspeci)*spinspec(nspecj)*spinspec(nspeck)
                            else
                              rk33 = rkthb3*ofct
                            endif
                            if (n3ty.eq.17) then
                              rk34 = rkthb4
                            else
                              rk34 = rkthb4*ofct
                            endif
                            if (n3ty.eq.15) then
                              rho1 = thrho1(n)
                              rho2 = thrho2(n)
                              rho3 = thrho3(n)
                            elseif (n3ty.eq.16) then
                              rho1 = thrho1(n)
                              rho2 = thrho2(n)
                            endif
!
!  Set region 2 trio flag
!
                            lreg12    = .false.
                            lreg2trio = .false.
                            if (lseok.and.nregions(ncf).gt.1) then
                              lreg2trio = (nregioni.gt.1.and.nregionj.gt.1.and.nregionk.gt.1)
                              if (.not.lreg2trio) lreg12 = (nregioni.gt.1.or.nregionj.gt.1.or.nregionk.gt.1)
                            endif
                            lslicek = lsliceatom(nsft + nrelf2a(k))
                            lattach = .true.
                            if (lslicei.and.lslicej.and.lslicek) lattach = .false.
                            if (.not.lslicei.and..not.lslicej.and..not.lslicek) lattach = .false.
!*****************************************************
!  Calculate derivatives with respect to potentials  *
!*****************************************************
                            call threebody(1_i4,n3ty,r12,r13,r23,ed11,ed12,ed13,ethb1,e2d,e3d,ttr11,ttr21,tr31, &
                                           rho1,rho2,rho3,rho4,rho5,rk32,rk33,rk34,the0,ang,dot,lgrad1,.false., &
                                           .false.,n,qli,qlj,qlk,d0i,d0j,d0k,d1q,d2q,lthetataper(n), &
                                           thetatapermin(n),thetatapermax(n))
                            if (lreg2trio) then
                              esregion2l = esregion2l + ethb1
                            elseif (lreg12) then
                              esregion12l = esregion12l + ethb1
                            else
                              ethb = ethb + ethb1
                            endif
                            if (lattach) eattachl = eattachl + ethb1
!
                            eregion2region(nregionj,nregioni) = eregion2region(nregionj,nregioni) + third*ethb1
                            eregion2region(nregionk,nregioni) = eregion2region(nregionk,nregioni) + third*ethb1
                            eregion2region(nregionk,nregionj) = eregion2region(nregionk,nregionj) + third*ethb
!
                            siteenergy(i) = siteenergy(i) + third*ethb1
                            siteenergy(j) = siteenergy(j) + third*ethb1
                            siteenergy(k) = siteenergy(k) + third*ethb1
!
!  Output energy contribution
!
                            if (lPrintThree) then
                              write(ioout,'(4x,3i12,8x,f27.10)') i,j,k,ethb1
                            endif
!*************************
!  Start of derivatives  *
!*************************
                            if (lgrad1) then
!
!  Set up strain products
!
                              if (lsg1) then
                                xv3(1) = x21
                                xv3(2) = x31
                                xv3(3) = x23
                                yv3(1) = y21
                                yv3(2) = y31
                                yv3(3) = y23
                                zv3(1) = z21
                                zv3(2) = z31
                                zv3(3) = z23
                                call realstrterms(ndim,3_i4,3_i4,xv3,yv3,zv3,xcom,ycom,zcom,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
                              endif
!   
!  Charge derivatives
!   
                              if (lDoQDeriv1) then
                                call d1charge3(i,j,k,lopi,lopj,lopk,1_i4,d0i,d0j,d0k)
                              endif
                            endif
!***********************
!  Strain derivatives  *
!***********************
                            if (lsg1) then
!
!  First strain derivatives
!
                              rstrdloc(1:nstrains) = 0.0_dp
                              do kl = 1,nstrains
                                ks = nstrptr(kl)
                                rstrdloc(kl) = rstrdloc(kl) + ed11*dr2ds(ks,1)
                                rstrdloc(kl) = rstrdloc(kl) + ed12*dr2ds(ks,2)
                                rstrdloc(kl) = rstrdloc(kl) + ed13*dr2ds(ks,3)
                                rstrd(kl) = rstrd(kl) + rstrdloc(kl)
                              enddo
                              if (latomicstress) then
                                do kl = 1,nstrains
                                  atomicstress(kl,i) = atomicstress(kl,i) + third*rstrdloc(kl)
                                  atomicstress(kl,j) = atomicstress(kl,j) + third*rstrdloc(kl)
                                  atomicstress(kl,k) = atomicstress(kl,k) + third*rstrdloc(kl)
                                enddo
                              endif
                            endif
!*************************
!  Internal derivatives  *
!*************************
                            if (lgrad1) then
                              xderv(i) = xderv(i) - x21*ed11 - x31*ed12
                              yderv(i) = yderv(i) - y21*ed11 - y31*ed12
                              zderv(i) = zderv(i) - z21*ed11 - z31*ed12
                              xderv(j) = xderv(j) + x21*ed11 - x23*ed13
                              yderv(j) = yderv(j) + y21*ed11 - y23*ed13
                              zderv(j) = zderv(j) + z21*ed11 - z23*ed13
                              xderv(k) = xderv(k) + x31*ed12 + x23*ed13
                              yderv(k) = yderv(k) + y31*ed12 + y23*ed13
                              zderv(k) = zderv(k) + z31*ed12 + z23*ed13
                            endif
!
!  End of condition on k being greater than j
!
                          endif
!
!  End of loop over atom k
!
                        enddo kloop
!
!  End of loops over neighbouring cells for k
!
                      enddo
                    enddo
                  enddo
!
!  End of loop over atom j
!
                enddo jloop
!
!  End of loops over neighbouring cells for j
!
              enddo
            enddo
          enddo
!
!  End of loop over atom i
!
        enddo iloop
!
!  End if for non-buffer cell
!
      endif
!
!  End of loop over central cell
!
    enddo
!
!  End loop over threebody terms
!
  enddo potentials
!
!  Closing banner for energy decomposition
!
  if (lPrintThree) then
    call mpbarrier
    if (ioproc) then
      write(ioout,'(''--------------------------------------------------------------------------------'',/)')
    endif
  endif
!
!  Marvin compatibility option -> all three body terms are in region 1
!
  if (lmarvreg2) then
    ethb = ethb + esregion12l
  else
    esregion12 = esregion12 + esregion12l
    esregion2 = esregion2 + esregion2l
  endif
!
!  If symmetry adapted derivatives have been calculated elsewhere
!  then add derivatives of related atoms
!
  if (lgrad1) then
    if (lsymderv) then
      do i = 1,nasym
        nr = nrela2f(i)
        neq = neqv(i)
        xdrv(i) = xdrv(i) + neq*xderv(nr)
        ydrv(i) = ydrv(i) + neq*yderv(nr)
        zdrv(i) = zdrv(i) + neq*zderv(nr)
!
        nregioni = nregionno(nsft+i)
        xregdrv(nregioni) = xregdrv(nregioni) + neq*xderv(nr)
        yregdrv(nregioni) = yregdrv(nregioni) + neq*yderv(nr)
        zregdrv(nregioni) = zregdrv(nregioni) + neq*zderv(nr)
      enddo
    else
      do i = 1,numat
        xdrv(i) = xdrv(i) + xderv(i)
        ydrv(i) = ydrv(i) + yderv(i)
        zdrv(i) = zdrv(i) + zderv(i)
!
        nregioni = nregionno(nsft+nrelf2a(i))
        xregdrv(nregioni) = xregdrv(nregioni) + xderv(i)
        yregdrv(nregioni) = yregdrv(nregioni) + yderv(i)
        zregdrv(nregioni) = zregdrv(nregioni) + zderv(i)
      enddo
    endif
  endif
!
!  Free local memory
!
  deallocate(zderv,stat=status)
  if (status/=0) call deallocate_error('threemds','zderv')
  deallocate(yderv,stat=status)
  if (status/=0) call deallocate_error('threemds','yderv')
  deallocate(xderv,stat=status)
  if (status/=0) call deallocate_error('threemds','xderv')
!
!  Timing
!
  time2 = g_cpu_time()
  tthree = tthree + time2 - time1
#ifdef TRACE
  call trace_out('threemds')
#endif
!
  return
  end
