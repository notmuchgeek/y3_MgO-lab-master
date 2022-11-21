  subroutine raman0(qD,maxqD)
!
!  Computes the atomic Raman susceptibilites using third derivatives: 0-D
!
!  NOTE : Needs finishing for BSM 
!
!   9/13 Created from raman3
!   8/17 Parallelisation added
!   2/18 Trace added
!   8/19 Corrections to call to twobody1
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!
!  On entry:
!
!  qD = matrix containing inverse shell-shell second derivative matrix multiplied by charges of shells
!
!  nmanyk = number of three-/four-body term K atoms whose derivatives
!          depend on the i-j block
!  nptrmanyk = pointer to the atoms for each nmanyk
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
!  Copyright Curtin University 2019
!
!  Julian Gale, CIC, Curtin University, September 2019
!
  use configurations, only : lbsmat
  use g_constants
  use control
  use current
  use datatypes
  use derivatives
  use eam,            only : maxmeamcomponent
  use element,        only : maxele
  use feworkspace
  use four
  use general
  use iochannels
  use ksample
  use kspace
  use m_three
  use molecule
  use parallel
  use partial,        only : iocshptr
  use properties,     only : ramanasus
  use shells
  use symmetry
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  use two
  implicit none
!
!  Passed variables
!
  integer(i4)                                  :: maxqD
  real(dp)                                     :: qD(maxqD,*)
!
!  Local variables
!
  character(len=4)                             :: cstype
  character(len=5)                             :: lab
  integer(i4)                                  :: i
  integer(i4)                                  :: inat
  integer(i4)                                  :: ind
  integer(i4)                                  :: indi
  integer(i4)                                  :: indj
  integer(i4)                                  :: indm
  integer(i4)                                  :: indmi
  integer(i4)                                  :: indmj
  integer(i4)                                  :: ioc
  integer(i4)                                  :: ish
  integer(i4)                                  :: itype
  integer(i4)                                  :: ix
  integer(i4)                                  :: iy
  integer(i4)                                  :: iz
  integer(i4)                                  :: ixx
  integer(i4)                                  :: iyy
  integer(i4)                                  :: izz
  integer(i4)                                  :: j
  integer(i4)                                  :: joc
  integer(i4)                                  :: jsh
  integer(i4)                                  :: jx
  integer(i4)                                  :: jy
  integer(i4)                                  :: jz
  integer(i4)                                  :: jxx
  integer(i4)                                  :: jyy
  integer(i4)                                  :: jzz
  integer(i4)                                  :: maxlim
  integer(i4)                                  :: mint
  integer(i4)                                  :: n
  integer(i4)                                  :: nat1
  integer(i4)                                  :: nat2
  integer(i4)                                  :: nati
  integer(i4)                                  :: natj
  integer(i4)                                  :: nbtypeij
  integer(i4)                                  :: nbtypeij2
  integer(i4)                                  :: nforkl
  integer(i4)                                  :: nmanyk
  integer(i4)                                  :: nmi
  integer(i4)                                  :: nmj
  integer(i4)                                  :: nor
  integer(i4), dimension(:), allocatable       :: npotl
  integer(i4)                                  :: npots
  integer(i4)                                  :: ntyp1
  integer(i4)                                  :: ntyp2
  integer(i4)                                  :: ntypi
  integer(i4)                                  :: ntypj
  integer(i4)                                  :: status
  logical                                      :: l2bonds
  logical                                      :: l3bonds
  logical                                      :: lbonded
  logical                                      :: lcspair
  logical                                      :: lcorei
  logical                                      :: lcorej
  logical                                      :: lfor
  logical                                      :: lthb
  logical                                      :: lmany
  logical                                      :: lmatch
  logical                                      :: lmolok
  logical                                      :: lneedmol 
  logical                                      :: lptrmol 
  logical                                      :: lorder12loc
  real(dp)                                     :: g_cpu_time
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
  real(dp)                                     :: d3(3,3,3)
  real(dp)                                     :: d2dx(3,3)
  real(dp)                                     :: d2dy(3,3)
  real(dp)                                     :: d2dz(3,3)
  real(dp)                                     :: deriv
  real(dp)                                     :: deriv2
  real(dp)                                     :: deriv3
  real(dp)                                     :: derive0
  real(dp)                                     :: derive
  real(dp)                                     :: derive2
  real(dp)                                     :: derive3
  real(dp)                                     :: eatom
  real(dp)                                     :: ec6
  real(dp)                                     :: ereal
  real(dp)                                     :: factor
  real(dp)                                     :: oci
  real(dp)                                     :: ocj
  real(dp)                                     :: ofct
  real(dp)                                     :: ospfct
  real(dp)                                     :: qli
  real(dp)                                     :: qlj
  real(dp)                                     :: r
  real(dp)                                     :: radi
  real(dp)                                     :: radj
  real(dp)                                     :: radsum
  real(dp)                                     :: rderiv
  real(dp)                                     :: rpd1
  real(dp)                                     :: rpd2
  real(dp)                                     :: rpd3
  real(dp)                                     :: rpd4
  real(dp)                                     :: rpd5
  real(dp)                                     :: rpd6
  real(dp)                                     :: rp
  real(dp)                                     :: rtrm1
  real(dp)                                     :: rtrm2
  real(dp)                                     :: rtrm3
  real(dp)                                     :: rtrm32
  real(dp)                                     :: sctrm1(maxmeamcomponent)
  real(dp)                                     :: sctrm2(maxmeamcomponent)
  real(dp)                                     :: small
  real(dp)                                     :: small2
  real(dp),    dimension(:,:,:,:), allocatable :: sum4
  real(dp)                                     :: t1
  real(dp)                                     :: t2
  real(dp)                                     :: time1
  real(dp)                                     :: time2
  real(dp)                                     :: tproj0
  real(dp)                                     :: xal
  real(dp)                                     :: yal
  real(dp)                                     :: zal
  real(dp)                                     :: xcrd
  real(dp)                                     :: ycrd
  real(dp)                                     :: zcrd
  real(dp)                                     :: xd2
  real(dp)                                     :: yd2
  real(dp)                                     :: zd2
  real(dp)                                     :: xd3
  real(dp)                                     :: yd3
  real(dp)                                     :: zd3
#ifdef TRACE
  call trace_in('raman0')
#endif
!
  time1 = g_cpu_time()
  tproj0 = tproj
!
!  Initialise tensors
!
  ramanasus(1:3,1:3,1:3,1:numat) = 0.0_dp
!
!  Local variables
!
  small = 1.0d-12
  small2 = 1.0d-2
  lthb = (nthb.gt.0)
  lfor = (nfor.gt.0)
  lmany = (lthb.or.lfor)
!
!  Build potential lists if needed
!
  if (lfor) call setlist4
!  
!  Set the Coulomb term type based on dimensionality :
!     
!  1-D         => 1/r
!  2-D and 3-D => erfc(seta*r)/r 
!
  mint = 3*numat
  maxlim = mint
  if (nbsmat.gt.0) maxlim = maxlim + numat
!
!  Zero energies although not needed to avoid overflow
!
  eatom = 0.0_dp
  ereal = 0.0_dp
  ec6 = 0.0_dp
!
!  Set up cutoffs
!
  cut2p = cutp*cutp
  cut2s = cuts*cuts
  cut2w = cutw*cutw
  if (lwolf) then
    cut2q = cut2w
  else
    cut2q = 1.0d12
  endif
!
!  Allocate local memory
!
  allocate(npotl(npote),stat=status)
  if (status/=0) call outofmemory('raman0','npotl')
!***************************************************************
!  Atomistic and real space electrostatic second derivatives   *
!***************************************************************
!
!  Outer loop over sites
!
  do ish = 1,nshellonnode
    i = node2atom(nshonnodeptr(ish))
!
!  Inner loop over second site
!
    xal = xclat(i)
    yal = yclat(i)
    zal = zclat(i)
    nati = nat(i)
    ntypi = nftype(i)
    lcorei = (nati.le.maxele)
    qli = qf(i)
    oci = occuf(i)
    if (nprocs.gt.1) then
      ioc = i - ncore
    else
      ioc = iocshptr(ish)
    endif
    if (lbsmat(nsft+nrelf2a(i))) then
      radi = radf(i)
    else
      radi = 0.0_dp
    endif
    indi = 3*(ioc-1)
    ix = indi + 1
    iy = indi + 2
    iz = indi + 3
!
!  Molecule handling
!
    if (lmol) then
      nmi = natmol(i)
      indm = nmolind(i)
    endif
!
!  Start of second atom loop
!
    do j = 1,numat
      natj = nat(j)
      ntypj = nftype(j)
      lcorej = (natj.le.maxele)
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
      xcrd = xclat(j) - xal
      ycrd = yclat(j) - yal
      zcrd = zclat(j) - zal
      qlj = qf(j)
      ocj = occuf(j)
      if (lbsmat(nsft+nrelf2a(j))) then
        radj = radf(j)
      else
        radj = 0.0_dp
      endif
      radsum = radi + radj
      ofct = oci*ocj
      if (.not.lcorej) then
        jsh = j - ncore
        joc = iocshptr(jsh)
        indj = 3*(joc - 1)
        jx = indj + 1
        jy = indj + 2
        jz = indj + 3
      endif
      factor = qli*qlj*ofct*angstoev
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
!  Molecule handling
!
      if (lmol) then
        nmj = natmol(j)
        indmj = nmolind(j)
      endif
!
!  Locate potential number
!  Check whether potential requires specific types
!
      rp = 0.0_dp
      npots = 0
      lneedmol = (lmol.and..not.lmolq)
      do n = 1,npote
        if (lmatch(nat1,ntyp1,nspec1(n),nptyp1(n),.true.)) then
          if (lmatch(nat2,ntyp2,nspec2(n),nptyp2(n),.true.)) then
            npots = npots + 1
            npotl(npots) = n
            if (mmexc(n).gt.0.or.(lintra(n).and..not.linter(n).or..not.lintra(n).and.linter(n))) lneedmol = .true.
            if (nptype(n).eq.8.or.nptype(n).eq.33) then
              if (cuts.gt.rp) rp = cuts
            else
              if (rpot(n).gt.rp) rp = rpot(n)
            endif
          endif
        endif
      enddo
      cut2r = rp*rp
      if (cut2r.gt.cut2p) cut2r = cut2p
      r = xcrd*xcrd + ycrd*ycrd + zcrd*zcrd
!
!  Molecule and bonding checks
!
      if (lmol) then
        lmolok = (nmi.eq.nmj.and.nmi.ne.0)
      else
        lmolok = .false.
      endif
!  
!  If the Coulomb term and the potentials don't need molecule info for
!  for this pair of atoms, then turn off molecule checking to save time
!
      if (.not.lneedmol) lmolok = .false.
!
      if (lmolok.and.(r.gt.cut2s.or..not.lcspair)) then
        ind = indmj - indmi
        lptrmol = (ind.eq.0)
        if (.not.lptrmol) then
          call mindtoijk(indmj,jxx,jyy,jzz)
          call mindtoijk(indmi,ixx,iyy,izz)
          jxx = jxx - ixx
          jyy = jyy - iyy
          jzz = jzz - izz
          call samemol(lptrmol,nmi,jxx,jyy,jzz,0_i4,0_i4,0_i4)
        endif
        if (lptrmol) then
          call bonded3(lbonded,l2bonds,l3bonds,nbtypeij,nbtypeij2,i,j,0_i4,0_i4,0_i4)
        else
          lbonded = .false.
          l2bonds = .false.
          l3bonds = .false.
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
!  Zero third derivative arrays
!
      d3(1:3,1:3,1:3) = 0.0_dp
!
      if (abs(r-small2).lt.1.0d-12) r = small2
      if (r.lt.small) then
!
!  Core-shell spring constant makes no contribution at small distances to third derivatives
!
        goto 1000
      else
!
!  Store vector
!
        nor = 1
        r = sqrt(r)
      endif
!*************************************************
!  Calculate twobody contribution in real space  *
!*************************************************
      call twobody1(eatom,ereal,ec6,.true.,.true.,.true.,nor,1_i4,r,xcrd,ycrd,zcrd,d0i,d0j, &
                    deriv,deriv2,deriv3,derive0,derive,derive2,derive3,rderiv,npots,npotl, &
                    cut2r,cut2q,cut2s,lptrmol,0_i4,factor,ofct,ospfct,radsum,rtrm1,rtrm2,rtrm3,rtrm32, &
                    sctrm1,sctrm2,qli,qlj,lcspair,.false.,.false.,lbonded,l2bonds,l3bonds, &
                    nbtypeij,nbtypeij2,.false.,.false.,d1i,d1j,d2i2,d2ij,d2j2)
!
!  Generate products for derivatives
!
      rpd1 = xcrd*xcrd
      rpd2 = ycrd*ycrd
      rpd3 = zcrd*zcrd
      rpd4 = ycrd*zcrd
      rpd5 = xcrd*zcrd
      rpd6 = xcrd*ycrd
      xd2 = xcrd*deriv2
      yd2 = ycrd*deriv2
      zd2 = zcrd*deriv2
      xd3 = xcrd*deriv3
      yd3 = ycrd*deriv3
      zd3 = zcrd*deriv3
!
!  Calculate third derivative matrix - first term
!
      d3(1,1,1) = d3(1,1,1) + rpd1*xd3
      d3(2,1,1) = d3(2,1,1) + rpd6*xd3
      d3(3,1,1) = d3(3,1,1) + rpd5*xd3
      d3(2,2,1) = d3(2,2,1) + rpd2*xd3
      d3(3,2,1) = d3(3,2,1) + rpd4*xd3
      d3(3,3,1) = d3(3,3,1) + rpd3*xd3
      d3(2,2,2) = d3(2,2,2) + rpd2*yd3
      d3(3,2,2) = d3(3,2,2) + rpd4*yd3
      d3(3,3,2) = d3(3,3,2) + rpd3*yd3
      d3(3,3,3) = d3(3,3,3) + rpd3*zd3
!
!  Add second term to third derivative matrix
!
      d3(1,1,1) = d3(1,1,1) + 3.0_dp*xd2
      d3(2,1,1) = d3(2,1,1) + yd2
      d3(3,1,1) = d3(3,1,1) + zd2
      d3(2,2,1) = d3(2,2,1) + xd2
      d3(3,3,1) = d3(3,3,1) + xd2
      d3(2,2,2) = d3(2,2,2) + 3.0_dp*yd2
      d3(3,2,2) = d3(3,2,2) + zd2
      d3(3,3,2) = d3(3,3,2) + yd2
      d3(3,3,3) = d3(3,3,3) + 3.0_dp*zd2
!
!  Symmetrise d3 matrices
!
      d3(1,2,1) = d3(2,1,1)
      d3(1,3,1) = d3(3,1,1)
      d3(2,3,1) = d3(3,2,1)
      d3(2,3,2) = d3(3,2,2)
      d3(1,1,2) = d3(2,1,1)
      d3(2,1,2) = d3(2,2,1)
      d3(3,1,2) = d3(3,2,1)
      d3(1,2,2) = d3(2,2,1)
      d3(1,3,2) = d3(3,2,1)
      d3(2,3,2) = d3(3,2,2)
      d3(1,1,3) = d3(3,1,1)
      d3(2,1,3) = d3(3,2,1)
      d3(3,1,3) = d3(3,3,1)
      d3(1,2,3) = d3(3,2,1)
      d3(2,2,3) = d3(3,2,2)
      d3(3,2,3) = d3(3,3,2)
      d3(1,3,3) = d3(3,3,1)
      d3(2,3,3) = d3(3,3,2)
!****************************
!  Three-body contribution  *
!****************************
      nmanyk = 0
      nforkl = 0
      if (lthb) then
        call three0d3(i,j,nati,ntypi,natj,ntypj,d3,xal,yal,zal,xcrd,ycrd,zcrd,nmanyk)
      endif
!***************************
!  Four-body contribution  *
!***************************
      if (lfor) then
        call four0d3(i,j,d3,nmanyk,nforkl)
      endif
!**************************************************
!  Project contributions on to shell coordinates  *
!**************************************************
      t1 = g_cpu_time()
      call projd0r(d3,i,ix,iy,iz,j,jx,jy,jz,lcorei,lcorej,d2dx,d2dy,d2dz,qD,maxqD, &
                   lmany,nmanyk,nptrmanyk,d33,nforkl,nptrfork,nptrforl,d34)
      t2 = g_cpu_time()
      tproj = tproj + t2 - t1
!
!  Add on derivative contributions
!
      ramanasus(1:3,1:3,1,i) = ramanasus(1:3,1:3,1,i) - d2dx(1:3,1:3)
      ramanasus(1:3,1:3,2,i) = ramanasus(1:3,1:3,2,i) - d2dy(1:3,1:3)
      ramanasus(1:3,1:3,3,i) = ramanasus(1:3,1:3,3,i) - d2dz(1:3,1:3)
!
      ramanasus(1:3,1:3,1,j) = ramanasus(1:3,1:3,1,j) + d2dx(1:3,1:3)
      ramanasus(1:3,1:3,2,j) = ramanasus(1:3,1:3,2,j) + d2dy(1:3,1:3)
      ramanasus(1:3,1:3,3,j) = ramanasus(1:3,1:3,3,j) + d2dz(1:3,1:3)
!
!  Skip to here if nothing to be done for this pair
!
1000  continue
    enddo
  enddo
!
!  Globalise ramanasus
!
  if (nprocs.gt.1) then
    allocate(sum4(3,3,3,numat),stat=status)
    if (status/=0) call outofmemory('raman0','sum4')
!
    call sumall(ramanasus,sum4,27_i4*numat,"raman0","ramanasus")
    ramanasus(1:3,1:3,1:3,1:numat) = angstoev*sum4(1:3,1:3,1:3,1:numat)
!
    deallocate(sum4,stat=status)
    if (status/=0) call deallocate_error('raman0','sum4')
  else
!
!  Multiply tensors by units conversion to make them Ang**-1
!
    ramanasus(1:3,1:3,1:3,1:numat) = angstoev*ramanasus(1:3,1:3,1:3,1:numat)
  endif
!*******************
!  Output tensors  *
!*******************
  if (ioproc) then
    write(ioout,'(/,''  Core and shell Raman susceptibility tensors (Ang**-1): '',/)')
    write(ioout,'(''-------------------------------------------------------------------------------'')')
    write(ioout,'(''  Atom                   x           y           z'')')
    write(ioout,'(''-------------------------------------------------------------------------------'')')
    do i = 1,numat
      inat = nat(i)
      itype = nftype(i)
      call label(inat,itype,lab)
      if (lbsmat(nrelf2a(i)+nsft)) then
        cstype = 'bcor'
        if (inat.gt.maxele) cstype = 'bshe'
      else
        cstype = 'core'
        if (inat.gt.maxele) cstype = 'shel'
      endif
      write(ioout,'(i4,1x,a5,1x,a4,1x,''x '',3(2x,f10.4))') i,lab,cstype,ramanasus(1,1,1,i),ramanasus(2,1,1,i),ramanasus(3,1,1,i)
      write(ioout,'(5x,''X'',10x,''y '',3(2x,f10.4))') ramanasus(1,2,1,i),ramanasus(2,2,1,i),ramanasus(3,2,1,i)
      write(ioout,'(16x,''z '',3(2x,f10.4))') ramanasus(1,3,1,i),ramanasus(2,3,1,i),ramanasus(3,3,1,i)
      write(ioout,'(16x,''x '',3(2x,f10.4))') ramanasus(1,1,2,i),ramanasus(2,1,2,i),ramanasus(3,1,2,i)
      write(ioout,'(5x,''Y'',10x,''y '',3(2x,f10.4))') ramanasus(1,2,2,i),ramanasus(2,2,2,i),ramanasus(3,2,2,i)
      write(ioout,'(16x,''z '',3(2x,f10.4))') ramanasus(1,3,2,i),ramanasus(2,3,2,i),ramanasus(3,3,2,i)
      write(ioout,'(16x,''x '',3(2x,f10.4))') ramanasus(1,1,3,i),ramanasus(2,1,3,i),ramanasus(3,1,3,i)
      write(ioout,'(5x,''Z'',10x,''y '',3(2x,f10.4))') ramanasus(1,2,3,i),ramanasus(2,2,3,i),ramanasus(3,2,3,i)
      write(ioout,'(16x,''z '',3(2x,f10.4))') ramanasus(1,3,3,i),ramanasus(2,3,3,i),ramanasus(3,3,3,i)
      write(ioout,'(''-------------------------------------------------------------------------------'')')
    enddo
    write(ioout,'(/)')
    call gflush(ioout)
  endif
!
!  Having output tensors we can now condense to cores only as this is what is needed for use in phonon calculations
!
!  Start by adding shells to their corresponding cores
!
  do ish = 1,nshell
    i = nshptr(ish)
    j = ncsptr(i)
    ramanasus(1:3,1:3,1,j) = ramanasus(1:3,1:3,1,j) + ramanasus(1:3,1:3,1,i)
    ramanasus(1:3,1:3,2,j) = ramanasus(1:3,1:3,2,j) + ramanasus(1:3,1:3,2,i)
    ramanasus(1:3,1:3,3,j) = ramanasus(1:3,1:3,3,j) + ramanasus(1:3,1:3,3,i)
  enddo
!
!  Now compact tensors into 1 -> ncore
!
  do i = 1,ncore
    j = ncoptr(i)
    ramanasus(1:3,1:3,1,i) = ramanasus(1:3,1:3,1,j)
    ramanasus(1:3,1:3,2,i) = ramanasus(1:3,1:3,2,j)
    ramanasus(1:3,1:3,3,i) = ramanasus(1:3,1:3,3,j)
  enddo
!
  if (ioproc) then
    write(ioout,'(/,''  Atomic Raman susceptibility tensors (Ang**-1): '',/)')
    write(ioout,'(''-------------------------------------------------------------------------------'')')
    write(ioout,'(''  Atom                   x           y           z'')')
    write(ioout,'(''-------------------------------------------------------------------------------'')')
    do i = 1,ncore
      j = ncoptr(i)
      inat = nat(j)
      itype = nftype(j)
      call label(inat,itype,lab)
      if (lbsmat(nrelf2a(j)+nsft)) then
        cstype = 'bcor'
        if (inat.gt.maxele) cstype = 'bshe'
      else
        cstype = 'core'
        if (inat.gt.maxele) cstype = 'shel'
      endif
      write(ioout,'(i4,1x,a5,1x,a4,1x,''x '',3(2x,f10.4))') i,lab,cstype,ramanasus(1,1,1,i),ramanasus(2,1,1,i),ramanasus(3,1,1,i)
      write(ioout,'(5x,''X'',10x,''y '',3(2x,f10.4))') ramanasus(1,2,1,i),ramanasus(2,2,1,i),ramanasus(3,2,1,i)
      write(ioout,'(16x,''z '',3(2x,f10.4))') ramanasus(1,3,1,i),ramanasus(2,3,1,i),ramanasus(3,3,1,i)
      write(ioout,'(16x,''x '',3(2x,f10.4))') ramanasus(1,1,2,i),ramanasus(2,1,2,i),ramanasus(3,1,2,i)
      write(ioout,'(5x,''Y'',10x,''y '',3(2x,f10.4))') ramanasus(1,2,2,i),ramanasus(2,2,2,i),ramanasus(3,2,2,i)
      write(ioout,'(16x,''z '',3(2x,f10.4))') ramanasus(1,3,2,i),ramanasus(2,3,2,i),ramanasus(3,3,2,i)
      write(ioout,'(16x,''x '',3(2x,f10.4))') ramanasus(1,1,3,i),ramanasus(2,1,3,i),ramanasus(3,1,3,i)
      write(ioout,'(5x,''Z'',10x,''y '',3(2x,f10.4))') ramanasus(1,2,3,i),ramanasus(2,2,3,i),ramanasus(3,2,3,i)
      write(ioout,'(16x,''z '',3(2x,f10.4))') ramanasus(1,3,3,i),ramanasus(2,3,3,i),ramanasus(3,3,3,i)
      write(ioout,'(''-------------------------------------------------------------------------------'')')
    enddo
    write(ioout,'(/)')
    call gflush(ioout)
  endif
!************************
!  Free dynamic memory  *
!************************
!
!  Free local memory
!
  deallocate(npotl,stat=status)
  if (status/=0) call deallocate_error('raman0','npotl')
!
!  Add on CPU time less time spent projecting derivatives
!
  time2 = g_cpu_time()
  tderv3 = tderv3 + time2 - time1 - tproj + tproj0
#ifdef TRACE
  call trace_out('raman0')
#endif
!
  return
  end
