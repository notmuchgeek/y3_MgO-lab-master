  subroutine realmd3s(eatom,ereal,erecip,ec6,eqeq,esregion12,esregion2,eattach,lgrad1)
!
!  Subroutine for calculating real space energy and gradients using a spatial decomposition
!  algorithm
!
!  lfreeze = if .true. only do derivatives for atoms with a degree of freedom
!
!   5/03 Created from realmd3
!   5/03 Region 3 modifications added
!   9/03 Flag for periodicity / molecule handling corrected
!  10/03 Modified parallel algorithm introduced to remove poor scaling
!        of nprocs > 1 section in brennermd
!   9/04 Charge first derivatives added
!   1/05 r2 sqrt'd before passing to twobody1
!   4/05 Mods for cosh-spring added
!   7/05 Streitz and Mintmire modifications added
!   2/07 Bonding types added
!   3/07 Printing of twobody energies added as an option
!   3/07 Bonding types modified
!   5/07 QM/MM scheme added
!   5/07 Argument list for twobody call modified
!   6/07 Structure of arrays for storing spatial distribution changed to 1-D
!  11/07 Unused variables cleaned up
!  11/07 Modifications for reaxFF Coulomb term added
!   1/08 lreaxFFqreal removed
!   3/08 ReaxFF charge scheme removed as this is handled elsewhere
!   4/08 Modified for variable domain size
!   4/08 xvec1cell replaced by xvec2cell etc for spatial algorithm
!   4/08 ind1toijk replaced by ind2toijk
!  11/08 x/y/z components passed to twobody1
!  12/08 scrho, and defect versions, converted to 2-D array for MEAM
!   1/09 Arrays apot, bpot, cpot & dpot replaced by new 2-D array twopot
!   2/09 Expanding and contracting of maxdis arrays removed
!   3/09 lorder12loc added to twobody1 argument list
!   3/09 Value of distance check for lself replaced by global value smallself from general module
!   3/09 Breathing shell self term added
!   4/09 Globalisation of scrho values now only done if nprocs=1
!   4/09 MEAM density modifications removed since these are now handled in 
!        a separate routine.
!   6/09 Site energy and virials added
!  11/09 Region derivatives added
!  11/11 Region-region energy contributions stored
!   4/12 Explicit virial calculation removed as no longer needed
!   4/12 xvir, yvir and zvir removed
!   5/12 Atomic stresses added
!   7/13 Call to twobody modified to include core-shell occupancy factor
!   3/14 lorder12loc removed twobody argument list
!   2/15 MM3buck added
!   1/18 Cell index shifts handled to allow for nomod
!   2/18 Trace added
!   5/18 Modified for q ranges
!   9/18 Modified for changes due to lstraincell algorithm
!   9/18 Strain module introduced
!   5/19 Skipping of ij loop now checks for polarisation
!   8/19 Short range damping of polarisation added
!   8/19 Trap for neemrptr being zero added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!  11/19 Rigid molecule modifications added
!   2/20 Correction to cell index handling for bonding
!   2/20 Breathing shell contribution now made conditional on lopi
!   3/20 Electrostatic cutoff only included where the either charge exceeds the threshold
!   4/20 derv3c and d2r2dsdc added for benefit of rigid molecules
!   4/20 derv3c changes reversed as they are no longer required
!   8/20 if statement for not cell buffer moved to avoid unnecessary calculation
!  10/20 eqeq now included in site energy and eregion2region
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
  use configurations, only : lbsmat,lsliceatom,nregions,nregionno, nregiontype, QMMMmode
  use g_constants
  use control
  use current
  use datatypes
  use derivatives,    only : xdrv, ydrv, zdrv, rstrd, raderv
  use derivatives,    only : xregdrv, yregdrv, zregdrv, atomicstress
  use eam,            only : lMEAMden
  use eemdata
  use element
  use energies,       only : siteenergy, eregion2region
  use general,        only : cutw, smallself
  use iochannels,     only : ioout
  use kspace
  use m_strain,       only : twostrterms
  use mdlogic
  use molecule
  use optimisation
  use parallel
  use polarise
  use potentialxyz
  use realvectors,    only : dr2ds, d2r2dx2, d2r2ds2, d2r2dsdx, derivqd
  use shells
  use spatial
  use sutton
  use symmetry
  use thresholds,     only : thresh_q
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  use two
  implicit none
!
!  Passed variables
!
  logical,     intent(in)                      :: lgrad1
  real(dp),    intent(inout)                   :: eatom
  real(dp),    intent(inout)                   :: eattach
  real(dp),    intent(inout)                   :: ec6
  real(dp),    intent(inout)                   :: eqeq
  real(dp),    intent(inout)                   :: ereal
  real(dp),    intent(inout)                   :: erecip
  real(dp),    intent(inout)                   :: esregion12
  real(dp),    intent(inout)                   :: esregion2
!
!  Local variables
!
  integer(i4), dimension(:), allocatable       :: npotl
  integer(i4)                                  :: i
  integer(i4)                                  :: ic
  integer(i4)                                  :: icx
  integer(i4)                                  :: icy
  integer(i4)                                  :: icz
  integer(i4)                                  :: ii
  integer(i4)                                  :: ijcx
  integer(i4)                                  :: ijcy
  integer(i4)                                  :: ijcz
  integer(i4)                                  :: imx
  integer(i4)                                  :: imy
  integer(i4)                                  :: imz
  integer(i4)                                  :: indm
  integer(i4)                                  :: indmj
  integer(i4)                                  :: ind
  integer(i4)                                  :: ind2
  integer(i4)                                  :: indn
  integer(i4)                                  :: ix
  integer(i4)                                  :: iy
  integer(i4)                                  :: iz
  integer(i4)                                  :: ixi
  integer(i4)                                  :: iyi
  integer(i4)                                  :: izi
  integer(i4)                                  :: ixj
  integer(i4)                                  :: iyj
  integer(i4)                                  :: izj
  integer(i4)                                  :: ixyz
  integer(i4)                                  :: j
  integer(i4)                                  :: jc
  integer(i4)                                  :: jcx
  integer(i4)                                  :: jcy
  integer(i4)                                  :: jcz
  integer(i4)                                  :: jj
  integer(i4)                                  :: kl
  integer(i4)                                  :: ks
  integer(i4)                                  :: maxxy
  integer(i4)                                  :: maxx
  integer(i4)                                  :: n
  integer(i4)                                  :: n1i
  integer(i4)                                  :: n1j
  integer(i4)                                  :: nat1
  integer(i4)                                  :: nat2
  integer(i4)                                  :: nati
  integer(i4)                                  :: natj
  integer(i4)                                  :: nbtypeij
  integer(i4)                                  :: nbtypeij2
  integer(i4)                                  :: ni
  integer(i4)                                  :: nj
  integer(i4)                                  :: nmi
  integer(i4)                                  :: nmj
  integer(i4)                                  :: nmolonly
  integer(i4)                                  :: npots
  integer(i4)                                  :: nqri
  integer(i4)                                  :: nqrj
  integer(i4)                                  :: nregioni
  integer(i4)                                  :: nregionj
  integer(i4)                                  :: nregiontypi
  integer(i4)                                  :: nregiontypj
  integer(i4)                                  :: nsplower(3)
  integer(i4)                                  :: nspupper(3)
  integer(i4)                                  :: ntyp1
  integer(i4)                                  :: ntyp2
  integer(i4)                                  :: ntypi
  integer(i4)                                  :: ntypj
  integer(i4)                                  :: status
  logical                                      :: l2bonds
  logical                                      :: l3bonds
  logical                                      :: lattach
  logical                                      :: lbonded
  logical                                      :: lc6loc
  logical                                      :: lcspair
  logical                                      :: lewaldtype
  logical                                      :: lgrad1p
  logical                                      :: lmatch
  logical                                      :: lmdl
  logical                                      :: lmolok
  logical                                      :: lneedmol
  logical                                      :: lopi
  logical                                      :: lopj
  logical                                      :: lorder12loc
  logical                                      :: lptrmol
  logical                                      :: lQMMMelectro
  logical                                      :: lreg2one
  logical                                      :: lreg2pair
  logical                                      :: lsamemol
  logical                                      :: lself
  logical                                      :: lsg1  
  logical                                      :: lslicei
  logical                                      :: lslicej
  real(dp)                                     :: apt
  real(dp)                                     :: bpt
  real(dp)                                     :: c6tot
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: cut2
  real(dp)                                     :: cut2e
  real(dp)                                     :: cut2p
  real(dp)                                     :: cut2q
  real(dp)                                     :: cut2r
  real(dp)                                     :: cut2s
  real(dp)                                     :: cut2w
  real(dp)                                     :: d0i
  real(dp)                                     :: d0j
  real(dp)                                     :: d1i
  real(dp)                                     :: d1j
  real(dp)                                     :: d2i2
  real(dp)                                     :: d2ij
  real(dp)                                     :: d2j2
  real(dp)                                     :: d2self
  real(dp)                                     :: deriv
  real(dp)                                     :: deriv2
  real(dp)                                     :: deriv3
  real(dp)                                     :: derive0self
  real(dp)                                     :: derive0
  real(dp)                                     :: derive
  real(dp)                                     :: derive2
  real(dp)                                     :: derive3
  real(dp)                                     :: dist
  real(dp)                                     :: eatom2
  real(dp)                                     :: ec62 
  real(dp)                                     :: eqeq2  
  real(dp)                                     :: ereal2
  real(dp)                                     :: erecip2
  real(dp)                                     :: esum 
  real(dp)                                     :: etrm
  real(dp)                                     :: etrm1
  real(dp)                                     :: etrm2
  real(dp)                                     :: factor
  real(dp)                                     :: fct
  real(dp)                                     :: oci
  real(dp)                                     :: ocj
  real(dp)                                     :: ofct
  real(dp)                                     :: ospfct
  real(dp)                                     :: polfct
  real(dp)                                     :: qli
  real(dp)                                     :: qlj
  real(dp)                                     :: r2
  real(dp)                                     :: radi
  real(dp)                                     :: radj
  real(dp)                                     :: radsum
  real(dp)                                     :: rderiv
  real(dp)                                     :: rdiff
  real(dp)                                     :: rp
  real(dp)                                     :: rqeq2
  real(dp)                                     :: rstrdloc(6)
  real(dp)                                     :: rtrm1
  real(dp)                                     :: rtrm2
  real(dp)                                     :: rtrm3
  real(dp)                                     :: rtrm32
  real(dp)                                     :: sctrm1
  real(dp)                                     :: sctrm2
  real(dp),    dimension(:), allocatable       :: sum
  real(dp),    dimension(:), allocatable       :: sum2
  real(dp)                                     :: time1
  real(dp)                                     :: time2
  real(dp)                                     :: tsum0
  real(dp)                                     :: tsuml
  real(dp)                                     :: xcom
  real(dp)                                     :: ycom
  real(dp)                                     :: zcom
  real(dp)                                     :: xcomi
  real(dp)                                     :: ycomi
  real(dp)                                     :: zcomi
  real(dp)                                     :: xi
  real(dp)                                     :: yi
  real(dp)                                     :: zi
  real(dp)                                     :: xji
  real(dp)                                     :: yji
  real(dp)                                     :: zji
  real(dp)                                     :: x1(1,1)
  real(dp)                                     :: y1(1,1)
  real(dp)                                     :: z1(1,1)
#ifdef TRACE
  call trace_in('realmd3s')
#endif
!
  time1 = g_cpu_time()
!
!  Local variables
!
  lmdl = lmd
  if (.not.lgrad1) lmdl = .false.
  lgrad1p = (lgrad1.or.lpolar)
  lc6loc = (lc6.and.ndim.eq.3)
  rqeq2 = rqeq*rqeq
!
!  Set the Coulomb term type based on dimensionality :
!
!  1-D         => 1/r
!  2-D and 3-D => erfc(seta*r)/r
!
  lewaldtype = (ndim.ne.1.or.lwolf)
!
  tsuml = 0.0_dp
  lsg1 = (lgrad1.and.lstr)
!
!  Set up cutoffs
!
  cut2p = cutp*cutp
  cut2e = rmx2
  cut2s = cuts*cuts
  cut2w = cutw*cutw
  if (lwolf) then
    cut2q = cut2w
  else
    cut2q = cut2e
  endif
!
!  Allocate local memory
!
  allocate(npotl(npote),stat=status)
  if (status/=0) call outofmemory('realmd3s','npotl')
  allocate(sum(numat),stat=status)
  if (status/=0) call outofmemory('realmd3s','sum')
  allocate(sum2(numat),stat=status)
  if (status/=0) call outofmemory('realmd3s','sum2')
!
  if (.not.lnoreal) then
!
!  Opening banner for energy decomposition
!
    if (lPrintTwo) then
      call mpbarrier
      if (ioproc) then
        write(ioout,'(''--------------------------------------------------------------------------------'')')
        write(ioout,'(''  Two : Atom No. 1  Atom No. 2    Short-range energy (eV)   Coulomb energy (eV) '')')
        write(ioout,'(''--------------------------------------------------------------------------------'')')
      endif
    endif
!***************************************************************
!  Atomistic and real space electrostatic component of energy  *
!***************************************************************
    maxxy = nspcell(1)*nspcell(2)
    maxx  = nspcell(1)
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
        do ii = 1,ni
          i = nspcellatptr(n1i+ii)
          ic = nspcellatptrcell(n1i+ii)
          call ind2toijk(ic,icx,icy,icz)
!
!  Set coordinates of atom i
!
          xi = xinbox(i) + xvec2cell(ic)
          yi = yinbox(i) + yvec2cell(ic)
          zi = zinbox(i) + zvec2cell(ic)
!
!  Set other properties of atom i
!
          nati = nat(i)
          ntypi = nftype(i)
          qli = qf(i)
          oci = occuf(i)
          nregioni = nregionno(nsft+nrelf2a(i))
          nregiontypi = nregiontype(nregioni,ncf)
          lopi = (.not.lfreeze.or.lopf(nrelf2a(i)))
          lslicei = lsliceatom(nsft + nrelf2a(i))
          if (lbsmat(nrelf2a(i)+nsft)) then
            radi = radf(i)
          else
            radi = 0.0_dp
          endif
          if (leem.and.lmultiqrange.and.neemrptr(i).ne.0) then
            nqri = nqrnow(neemrptr(i))
          else
            nqri = 1
          endif
!
!  Molecule handling for atom i
!
          if (lmol) then
            nmi = natmol(i)
            indm = nmolind(i)
            call mindtoijk(indm,ixi,iyi,izi)
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
                  jc = nspcellatptrcell(n1j+jj)
                  call ind2toijk(jc,jcx,jcy,jcz)
!
!  Only calculate lower-half triangular interactions
!
                  if (j.le.i) then
!
!  Freezing flag
!
                    lopj = (.not.lfreeze.or.lopf(nrelf2a(j)))
                    if (.not.lopi.and..not.lopj) cycle jloop
!
!  Set coordinate differences and calculate square of distance
!  
                    xji = xvec2cell(jc) + xinbox(j) - xi
                    yji = yvec2cell(jc) + yinbox(j) - yi
                    zji = zvec2cell(jc) + zinbox(j) - zi
                    r2 = xji*xji + yji*yji + zji*zji
!
!  Set species type parameters for atom j
!
                    natj = nat(j)
                    ntypj = nftype(j)
!
!  Arrange nat1 to be lower than nat2. If they are equal
!  then arrange that ntyp1 is lower than ntyp2
!
                    if (nati.eq.natj) then
                      nat1 = nati
                      nat2 = natj
                      if (ntypi.lt.ntypj) then
                        lorder12loc = .true.
                        ntyp1 = ntypi
                        ntyp2 = ntypj
                      else
                        lorder12loc = .false.
                        ntyp1 = ntypj
                        ntyp2 = ntypi
                      endif
                    elseif (nati.lt.natj) then
                      lorder12loc = .true.
                      nat1 = nati
                      nat2 = nat(j)
                      ntyp1 = ntypi
                      ntyp2 = nftype(j)
                    else
                      lorder12loc = .false.
                      nat1 = nat(j)
                      nat2 = nati
                      ntyp1 = nftype(j)
                      ntyp2 = ntypi
                    endif
                    nregionj = nregionno(nsft+nrelf2a(j))
                    nregiontypj = nregiontype(nregionj,ncf)
!
!  Set region 2 pair flag
!
                    lreg2one  = .false.
                    lreg2pair = .false.
                    if (lseok.and.nregions(ncf).ge.2) then
                      lreg2pair = (nregioni.eq.2.and.nregionj.eq.2)
                      if (.not.lreg2pair) lreg2one = (nregioni.eq.2.or.nregionj.eq.2)
                    endif
                    lslicej = lsliceatom(nsft + nrelf2a(j))
                    lattach = (lslicei.and..not.lslicej.or.lslicej.and..not.lslicei)
!
!  QM/MM handling : i & j are both QM atoms => exclude
!
                    if (QMMMmode(ncf).gt.0) then
                      if (nregiontypi.eq.1.and.nregiontypj.eq.1) cycle jloop
                    endif
!           
!  QM/MM : Set electrostatic embedding flag : If either i or j are QM atoms => exclude electrostatics
!       
                    lQMMMelectro = (QMMMmode(ncf).eq.2.and.(nregiontypi.eq.1.or.nregiontypj.eq.1))
!
!  Set remaining properties for atom type j
!
                    qlj = qf(j)
                    ocj = occuf(j)
                    if (lbsmat(nsft+nrelf2a(j))) then
                      radj = radf(j)
                    else
                      radj = 0.0_dp
                    endif
                    radsum = radi + radj
                    if (i.eq.j) then
                      ofct = 0.5_dp*oci*ocj
                    else
                      ofct = oci*ocj
                    endif
                    fct = ofct*angstoev
                    factor = qli*qlj*fct
                    if (lpolar) then
                      polfct = abs(qli*dpolar(j)) + abs(qlj*dpolar(i))
                    else
                      polfct = 0.0_dp
                    endif
                    if (leem.and.lmultiqrange.and.neemrptr(j).ne.0) then
                      nqrj = nqrnow(neemrptr(j))
                    else
                      nqrj = 1
                    endif
!
!  Possible core-shell flag
!
                    lcspair = (abs(nat1-nat2).eq.maxele.or.(oci+ocj).lt.1.0001_dp)
                    if (lcspair) then
                      ospfct = oci
                    else
                      ospfct = ofct
                    endif
!
!  Locate potential number
!  Check whether potential requires specific types
!  Calculate sum of all dispersion terms for pair
!
                    rp = 0.0_dp
                    npots = 0
                    c6tot = 0.0_dp
                    lneedmol = (lmol.and..not.lmolq)
                    do n = 1,npote
                      if (lmatch(nat1,ntyp1,nspec1(n),nptyp1(n),.true.)) then
                        if (lmatch(nat2,ntyp2,nspec2(n),nptyp2(n),.true.)) then
                          npots = npots + 1
                          npotl(npots) = n
                          if (mmexc(n).gt.0.or.(lintra(n).and..not.linter(n).or..not.lintra(n).and.linter(n)))  &
                            lneedmol = .true.
                          if (nptype(n).eq.8.or.nptype(n).eq.33) then
                            if (cuts.gt.rp) rp = cuts
                          elseif (lc6loc) then
                            if (nptype(n).eq.1.or.nptype(n).eq.7) then
                              c6tot = c6tot + twopot(3,n)
                              if (repcut(n).gt.rp) rp = repcut(n)
                            elseif (nptype(n).eq.2.or.nptype(n).eq.21) then
                              c6tot = c6tot + twopot(2,n)
                              if (repcut(n).gt.rp) rp = repcut(n)
                            elseif (nptype(n).eq.57) then
                              c6tot = c6tot + twopot(3,n)*twopot(4,n)*twopot(5,n)**6
                              if (repcut(n).gt.rp) rp = repcut(n)
                            else
                              if (rpot(n).gt.rp) rp = rpot(n)
                            endif
                          else
                            if (rpot(n).gt.rp) rp = rpot(n)
                          endif
                        endif
                      endif
                    enddo
!
!  If no valid potentials, charge product is zero, and no polarisation then skip loop
!
                    if (npots.eq.0.and.abs(factor).lt.1.0d-8.and.polfct.lt.1.0d-8) cycle jloop
                    cut2r = rp*rp
                    if (cut2r.gt.cut2p) cut2r = cut2p
                    cut2 = cut2r
!
!  If both charges are less than threshold then exclude electrostatics from cutoff
!
                    if (abs(qli*oci)+abs(qlj*ocj).gt.thresh_q) then
                      if (cut2e.gt.cut2.and.lewald) cut2 = cut2e
                      if (cut2w.gt.cut2.and.lwolf) cut2 = cut2w
                    endif
                    if (lqeq.or.lSandM) cut2 = max(cut2,rqeq2)
!
!  Molecule handling
!
                    if (lmol) then
                      nmj = natmol(j)
                      indmj = nmolind(j)
                      call mindtoijk(indmj,ixj,iyj,izj)
!
                      ixj = ixj - ixi
                      iyj = iyj - iyi
                      izj = izj - izi
                      lmolok = (nmi.eq.nmj.and.nmi.gt.0)
!
!  Set COM coordinates
!
                      if (lrigid.and.nmj.gt.0) then
                        xcom = molxyz(1,natinmol(j),nmj) - xcomi
                        ycom = molxyz(2,natinmol(j),nmj) - ycomi
                        zcom = molxyz(3,natinmol(j),nmj) - zcomi
                      else
                        xcom = - xcomi
                        ycom = - ycomi
                        zcom = - zcomi
                      endif
                    else
                      lmolok = .false.
                      xcom = - xcomi
                      ycom = - ycomi
                      zcom = - zcomi
                    endif
!
!  If the Coulomb term and the potentials don't need molecule info for
!  for this pair of atoms, then turn off molecule checking to save time
!
                    if (.not.lneedmol) lmolok = .false.
                    nmolonly = 0
!
!  Molecule - check index
!
                    if (lmolok.and.(r2.gt.cut2s.or..not.lcspair)) then
                      ijcx = jcx - icx + icosxsp(j) - icosxsp(i)
                      ijcy = jcy - icy + icosysp(j) - icosysp(i)
                      ijcz = jcz - icz + icoszsp(j) - icoszsp(i)
!
                      call bonded3(lbonded,l2bonds,l3bonds,nbtypeij,nbtypeij2,i,j,ijcx,ijcy,ijcz)
                      lsamemol = (lbonded.or.l2bonds.or.l3bonds)
                      if (.not.lsamemol) then
                        call samemol(lsamemol,nmi,ijcx,ijcy,ijcz,ixj,iyj,izj)
                      endif
                      lptrmol = lsamemol
                      if (lsamemol) then
                        if (r2.gt.cut2e) nmolonly = 1
                      else
                        lbonded   = .false.
                        l2bonds   = .false.
                        l3bonds   = .false.
                        nbtypeij  = 0
                        nbtypeij2 = 0
                      endif
                    else
                      lptrmol   = .false.
                      lbonded   = .false.
                      l2bonds   = .false.
                      l3bonds   = .false.
                      nbtypeij  = 0
                      nbtypeij2 = 0
                    endif
!
!  Check distance against potential cutoff
!
                    if (r2.gt.cut2.and..not.lptrmol) cycle jloop
!
!  Set self term flag
!
                    lself = (r2.lt.smallself)
!
                    if (lself) then
                      d2self = 0.0_dp
                      ereal2 = 0.0_dp
                      erecip2 = 0.0_dp
                      ec62 = 0.0_dp
                      derive0self = 0.0_dp
                      call selfterm(ereal2,erecip2,ec62,derive0self,factor,fct,ofct,ospfct,1.0_dp,npotl,npots, &
                                    c6tot,d2self,lgrad1,.false.,i,j,0,0,1.0_dp,0.5_dp,lewaldtype,qli,qlj)
                      esum = ereal2 + erecip2 + ec62
                      if (lreg2one) then
                        esregion12 = esregion12 + esum
                      elseif (lreg2pair) then
                        esregion2 = esregion2 + esum
                      else
                        ereal = ereal + ereal2
                        erecip = erecip + erecip2
                        ec6 = ec6 + ec62
                      endif
                      if (lattach) eattach = eattach + esum
!
                      eregion2region(nregionj,nregioni) = eregion2region(nregionj,nregioni) + esum
!
                      siteenergy(i) = siteenergy(i) + 0.5_dp*esum
                      siteenergy(j) = siteenergy(j) + 0.5_dp*esum
                    endif
!*******************************************
!  Charge first derivatives for self-term  *
!*******************************************
                    if (lself.and.lgrad1.and.lDoQDeriv1) then
                      call d1charge(i,j,lopi,lopj,1_i4,derive0self*qlj,derive0self*qli)
                    endif
!
!  If self term, skip rest of working
!
                    if (lself) cycle jloop
!*************************************************
!  Calculate twobody contribution in real space  *
!*************************************************
                    eatom2 = 0.0_dp
                    ereal2 = 0.0_dp
                    ec62 = 0.0_dp
                    dist = sqrt(r2)
                    call twobody1(eatom2,ereal2,ec62,lgrad1p,.false.,.false.,1_i4,1_i4,dist,xji,yji,zji, &
                                  d0i,d0j,deriv,deriv2,deriv3,derive0,derive,derive2,derive3,rderiv,npots, &
                                  npotl,cut2r,cut2q,cut2s,lptrmol,nmolonly,factor,ofct,ospfct,radsum,rtrm1, &
                                  rtrm2,rtrm3,rtrm32,sctrm1,sctrm2,qli,qlj,lcspair,lewaldtype,.false., &
                                  lbonded,l2bonds,l3bonds,nbtypeij,nbtypeij2,.false.,lQMMMelectro, &
                                  d1i,d1j,d2i2,d2ij,d2j2)
                    esum = eatom2 + ereal2 + ec62
                    if (lreg2one) then
                      esregion12 = esregion12 + esum
                    elseif (lreg2pair) then
                      esregion2 = esregion2 + esum
                    else
                      eatom = eatom + eatom2
                      ereal = ereal + ereal2
                      ec6 = ec6 + ec62
                    endif
                    if (lattach) eattach = eattach + esum
!
                    eregion2region(nregionj,nregioni) = eregion2region(nregionj,nregioni) + esum
!
                    siteenergy(i) = siteenergy(i) + 0.5_dp*esum
                    siteenergy(j) = siteenergy(j) + 0.5_dp*esum
!
                    if (lPrintTwo) then
                      write(ioout,'(4x,i12,i12,f22.10,1x,f22.10)') i,j,eatom2+ec62,ereal2
                    endif
!
                    if (lDoQDeriv1.and.lgrad1) then
                      d0i = d0i + derive0*qlj
                      d0j = d0j + derive0*qli
                    endif
                    if (leem) then
                      if (lqeq.or.lSandM) then
                        eqeq2 = 0.0_dp
                        if (lqeq) then
                          call qeqbody(eqeq2,lgrad1,.false.,1_i4,1_i4,fct,qli,qlj,nati,natj,nqri,nqrj)
                        elseif (lSandM) then
                          call smbody(eqeq2,lgrad1,.false.,1_i4,1_i4,fct,qli,qlj,nati,natj,nqri,nqrj)
                        endif
                        if (lreg2one) then
                          esregion12 = esregion12 + eqeq2
                        elseif (lreg2pair) then
                          esregion2 = esregion2 + eqeq2
                        else
                          eqeq = eqeq + eqeq2
                        endif
                        if (lattach) eattach = eattach + eqeq2
!
                        eregion2region(nregionj,nregioni) = eregion2region(nregionj,nregioni) + eqeq2
!
                        siteenergy(i) = siteenergy(i) + 0.5_dp*eqeq2
                        siteenergy(j) = siteenergy(j) + 0.5_dp*eqeq2
                      endif
                    endif
                    if (lsuttonc) then
                      if (.not.lMEAMden) then
                        if (lorder12loc) then
                          scrho(1,i) = scrho(1,i) + sctrm1*ocj
                          scrho(1,j) = scrho(1,j) + sctrm2*oci
                          if (lattach) then
                            scrho12(1,i) = scrho12(1,i) + sctrm1*ocj
                            scrho12(1,j) = scrho12(1,j) + sctrm2*oci
                          endif
                        else
                          scrho(1,i) = scrho(1,i) + sctrm2*ocj
                          scrho(1,j) = scrho(1,j) + sctrm1*oci
                          if (lattach) then
                            scrho12(1,i) = scrho12(1,i) + sctrm2*ocj
                            scrho12(1,j) = scrho12(1,j) + sctrm1*oci
                          endif
                        endif
                      endif
                    endif
!
!  Generate products for derivatives
!
                    if (lmdl.or.lsg1) then
                      x1(1,1) = xji
                      y1(1,1) = yji
                      z1(1,1) = zji
                      call twostrterms(ndim,maxdis,1_i4,x1,y1,z1,xcom,ycom,zcom,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
                    endif
!*****************************
!  Charge first derivatives  *
!*****************************
                    if (lgrad1.and.lDoQDeriv1) then
                      call d1charge(i,j,lopi,lopj,1_i4,d0i,d0j)
                    endif
!*************************************
!  Electrostatic potential on-sites  *
!*************************************
                    if (lpolar) then
                      vx(i) = vx(i) - qlj*derivqd(1)*xji
                      vy(i) = vy(i) - qlj*derivqd(1)*yji
                      vz(i) = vz(i) - qlj*derivqd(1)*zji
                      vx(j) = vx(j) + qli*derivqd(1)*xji
                      vy(j) = vy(j) + qli*derivqd(1)*yji
                      vz(j) = vz(j) + qli*derivqd(1)*zji
                      if (lattach) then
                        vx12(i) = vx12(i) - qlj*derivqd(1)*xji
                        vy12(i) = vy12(i) - qlj*derivqd(1)*yji
                        vz12(i) = vz12(i) - qlj*derivqd(1)*zji
                        vx12(j) = vx12(j) + qli*derivqd(1)*xji
                        vy12(j) = vy12(j) + qli*derivqd(1)*yji
                        vz12(j) = vz12(j) + qli*derivqd(1)*zji
                      endif
                    endif
!***********************
!  Radial derivatives  *
!***********************
                    if (lgrad1) then
                      if (radi.gt.0.0_dp.and.lopi) then
                        raderv(i) = raderv(i) + rtrm1
                      endif
                      if (radj.gt.0.0_dp.and.lopj) then
                        raderv(j) = raderv(j) + rtrm1
                      endif
                    endif
!************************
!  Internal Derivatives *
!************************
!
!  First derivatives
!
                    if (lgrad1) then
                      if (lopi) then
                        xdrv(i) = xdrv(i) - deriv*xji
                        ydrv(i) = ydrv(i) - deriv*yji
                        zdrv(i) = zdrv(i) - deriv*zji
                      endif
                      if (lopj) then
                        xdrv(j) = xdrv(j) + deriv*xji
                        ydrv(j) = ydrv(j) + deriv*yji
                        zdrv(j) = zdrv(j) + deriv*zji
                      endif
                      if (nregioni.ne.nregionj) then
                        xregdrv(nregioni) = xregdrv(nregioni) - deriv*xji
                        yregdrv(nregioni) = yregdrv(nregioni) - deriv*yji
                        zregdrv(nregioni) = zregdrv(nregioni) - deriv*zji
                        xregdrv(nregionj) = xregdrv(nregionj) + deriv*xji
                        yregdrv(nregionj) = yregdrv(nregionj) + deriv*yji
                        zregdrv(nregionj) = zregdrv(nregionj) + deriv*zji
                      endif
                    endif
!***********************
!  Strain derivatives  *
!***********************
!
!  Strain only terms 
!
!  First derivatives 
!
                    if (lsg1) then
                      rstrdloc(1:nstrains) = 0.0_dp
                      do kl = 1,nstrains
                        ks = nstrptr(kl)
                        rstrdloc(kl) = rstrdloc(kl) + deriv*dr2ds(1,ks)
                        rstrd(kl) = rstrd(kl) + deriv*dr2ds(1,ks)
                      enddo
                      if (latomicstress) then
                        do kl = 1,nstrains
                          atomicstress(kl,i) = atomicstress(kl,i) + 0.5_dp*rstrdloc(kl)
                          atomicstress(kl,j) = atomicstress(kl,j) + 0.5_dp*rstrdloc(kl)
                        enddo
                      endif
                    endif
                  endif
!
!  End loop over atom j
!
                enddo jloop
!
!  End loops over neighbouring cells
!
              enddo
            enddo
          enddo
!*********************
!  Radial self term  *
!*********************
          if (lbsmat(nsft+nrelf2a(i)).and.lopi) then
!
!  Find self term
!
            eatom2 = 0.0_dp
            do n = 1,npote
              if (nptype(n).eq.14) then
                if (lmatch(nati,ntypi,nspec1(n),nptyp1(n),.true.)) then
                  rdiff = radi - twopot(2,n)
                  apt = twopot(1,n)*oci
                  eatom2 = eatom2 + 0.5_dp*apt*rdiff*rdiff
                  if (lgrad1) then
                    raderv(i) = raderv(i) + apt*rdiff
                  endif
                endif
              elseif (nptype(n).eq.17) then
                if (lmatch(nati,ntypi,nspec1(n),nptyp1(n),.true.)) then
                  rdiff = radi - twopot(3,n)
                  apt = twopot(1,n)*oci
                  bpt = twopot(2,n)
                  etrm1 = exp(bpt*rdiff)
                  etrm2 = 1.0_dp/etrm1
                  etrm = apt*(etrm1 + etrm2)
                  eatom2 = eatom2 + etrm
                  if (lgrad1) then
                    raderv(i) = raderv(i) + apt*bpt*(etrm1 - etrm2)
                  endif
                endif
              elseif (nptype(n).eq.31) then
                if (lmatch(nati,ntypi,nspec1(n),nptyp1(n),.true.)) then
                  rdiff = radi - twopot(3,n)
                  apt = twopot(1,n)*oci
                  bpt = twopot(2,n)
                  etrm1 = exp(bpt*rdiff)
                  etrm = apt*etrm1
                  eatom2 = eatom2 + etrm
                  if (lgrad1) then
                    raderv(i) = raderv(i) + apt*bpt*etrm1
                  endif
                endif
              endif
            enddo
!
            eregion2region(nregioni,nregioni) = eregion2region(nregioni,nregioni) + eatom2
!
            siteenergy(i) = siteenergy(i) + eatom2
!
!  Set region 2 flag
!
            lreg2pair = .false.
            if (lseok.and.nregions(ncf).ge.2) then
              lreg2pair = (nregionno(nsft+nrelf2a(i)).eq.2)
            endif
            if (lreg2pair) then
              esregion2 = esregion2 + eatom2
            else
              eatom = eatom + eatom2
            endif
          endif
!
!  End loop over atom i
!
        enddo
!
!  End checks on whether cell is required
!
      endif
!
!  End loop over cells on node
!
    enddo
!****************
!  Global sums  *
!****************
    tsum0 = g_cpu_time()
    if (lsuttonc.and.nprocs.gt.1) then
      if (.not.lMEAMden) then
        do i = 1,numat
          sum2(i) = scrho(1,i) 
        enddo
        call sumall(sum2,sum,numat,"realmd3","scrho")
        do i = 1,numat
          scrho(1,i) = sum(i)
        enddo
        do i = 1,numat
          sum2(i) = scrho12(1,i) 
        enddo
        call sumall(sum2,sum,numat,"realmd3","scrho12")
        do i = 1,numat
          scrho12(1,i) = sum(i)
        enddo
      endif
    endif
    tsuml = g_cpu_time() - tsum0
    tsum = tsum + tsuml
  endif
!
!  Closing banner for energy decomposition
!
  if (lPrintTwo) then
    call mpbarrier
    if (ioproc) then
      write(ioout,'(''--------------------------------------------------------------------------------'')')
    endif
  endif
!
!  Free local memory
!
  deallocate(sum2,stat=status)
  if (status/=0) call deallocate_error('realmd3s','sum2')
  deallocate(sum,stat=status)
  if (status/=0) call deallocate_error('realmd3s','sum')
  deallocate(npotl,stat=status)
  if (status/=0) call deallocate_error('realmd3s','npotl')
!
!  Timing
!
  time2 = g_cpu_time()
  tatom = tatom + time2 - time1 - tsuml
#ifdef TRACE
  call trace_out('realmd3s')
#endif
!
  return
  end
