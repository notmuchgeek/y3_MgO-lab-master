  subroutine realfc
!
!  Routine for calculating the second derivatives for
!  use in a phonon calculation.
!
!  10/14 Created from realp
!  12/14 rtrm1 changed from scalar to array
!   2/15 MM3buck added
!   2/18 Trace added
!   5/18 Modified for q ranges
!   8/19 Modifications for Intel compiler
!   8/19 Trap for neemrptr being zero added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!   3/20 Electrostatic cutoff only included where the either charge exceeds the threshold
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
  use configurations, only : lbsmat
  use g_constants
  use control
  use current
  use datatypes
  use derivatives
  use eemdata
  use element
  use general,        only : cutw
  use kspace
  use molecule
  use parallel
  use realvectors
  use shells
  use symmetry
  use thresholds,     only : thresh_q
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  use two
  implicit none
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: iis
  integer(i4)                                  :: indi
  integer(i4)                                  :: indj
  integer(i4)                                  :: indm
  integer(i4)                                  :: indmj
  integer(i4)                                  :: indri
  integer(i4)                                  :: indrj
  integer(i4)                                  :: ix
  integer(i4)                                  :: iy
  integer(i4)                                  :: iz
  integer(i4)                                  :: ixi
  integer(i4)                                  :: iyi
  integer(i4)                                  :: izi
  integer(i4)                                  :: ixj
  integer(i4)                                  :: iyj
  integer(i4)                                  :: izj
  integer(i4)                                  :: j
  integer(i4)                                  :: jjs
  integer(i4)                                  :: jx
  integer(i4)                                  :: jy
  integer(i4)                                  :: jz
  integer(i4)                                  :: k
  integer(i4)                                  :: m
  integer(i4)                                  :: maxlim
  integer(i4)                                  :: mint
  integer(i4)                                  :: n
  integer(i4)                                  :: nat1
  integer(i4)                                  :: nat2
  integer(i4)                                  :: nati
  integer(i4)                                  :: natj
  integer(i4)                                  :: nbsi
  integer(i4)                                  :: nbsj
  integer(i4)                                  :: ncindm
  integer(i4)                                  :: ncindp
  integer(i4)                                  :: nmi
  integer(i4)                                  :: nmj
  integer(i4)                                  :: nmolonly
  integer(i4)                                  :: nor
  integer(i4), dimension(:), allocatable       :: npotl
  integer(i4)                                  :: npots
  integer(i4)                                  :: nqri
  integer(i4)                                  :: nqrj
  integer(i4)                                  :: ntyp1
  integer(i4)                                  :: ntyp2
  integer(i4)                                  :: ntypi
  integer(i4)                                  :: ntypj
  integer(i4)                                  :: status
  logical                                      :: lmatch
  logical                                      :: lc6loc
  logical                                      :: lcspair
  logical                                      :: lewaldtype
  logical                                      :: lmolok
  logical                                      :: lneedmol
  logical                                      :: lorder12loc
  logical                                      :: lself
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
  real(dp)                                     :: d1ix
  real(dp)                                     :: d1iy
  real(dp)                                     :: d1iz
  real(dp)                                     :: d1jx
  real(dp)                                     :: d1jy
  real(dp)                                     :: d1jz
  real(dp)                                     :: d2k
  real(dp)                                     :: d2self
  real(dp)                                     :: derive0self
  real(dp)                                     :: derive0selfi
  real(dp)                                     :: derive0selfj
  real(dp)                                     :: dk
  real(dp)                                     :: eatom
  real(dp)                                     :: ec6
  real(dp)                                     :: eqeq  
  real(dp)                                     :: ereal
  real(dp)                                     :: erecip
  real(dp)                                     :: etrm
  real(dp)                                     :: etrm1
  real(dp)                                     :: etrm2
  real(dp)                                     :: factor
  real(dp)                                     :: fct
  real(dp)                                     :: hfactor
  real(dp)                                     :: oci
  real(dp)                                     :: ocj
  real(dp)                                     :: ofct
  real(dp)                                     :: ospfct
  real(dp)                                     :: qli
  real(dp)                                     :: qlj
  real(dp)                                     :: radi
  real(dp)                                     :: radj
  real(dp)                                     :: radsum
  real(dp)                                     :: rdiff
  real(dp)                                     :: rp
  real(dp)                                     :: rqeq2
  real(dp)                                     :: sctrm1
  real(dp)                                     :: sctrm2
  real(dp)                                     :: time1
  real(dp)                                     :: time2
  real(dp)                                     :: xal
  real(dp)                                     :: yal
  real(dp)                                     :: zal
  real(dp)                                     :: xcrd
  real(dp)                                     :: ycrd
  real(dp)                                     :: zcrd
#ifdef TRACE
  call trace_in('realfc')
#endif
!
  time1 = g_cpu_time()
!
!  Local variables
!
  lc6loc = (lc6.and.ndim.eq.3)
  rqeq2 = rqeq*rqeq
  mint = 3*numat
  maxlim = mint
  if (nbsmat.gt.0) maxlim = maxlim + numat
  eatom = 0.0_dp
  ereal = 0.0_dp
  erecip = 0.0_dp
  ec6 = 0.0_dp
  eqeq = 0.0_dp
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
!  Set the Coulomb term type based on dimensionality :
!
!  1-D         => 1/r
!  2-D and 3-D => erfc(seta*r)/r
!
  lewaldtype = (ndim.ne.1.or.lwolf)
!
!  Allocate local memory
!
  allocate(npotl(npote),stat=status)
  if (status/=0) call outofmemory('realfc','npotl')
!
  if (.not.lnoreal) then
!***************************************************************
!  Atomistic and real space electrostatic second derivatives   *
!***************************************************************
!
!  Outer loop over sites
!
    nbsi = 0
    do i = 1,numat
!
!  Inner loop over second site
!
      xal = xclat(i)
      yal = yclat(i)
      zal = zclat(i)
      nati = nat(i)
      ntypi = nftype(i)
      qli = qf(i)
      oci = occuf(i)
      if (lbsmat(nsft+nrelf2a(i))) then
        radi = radf(i)
        nbsi = nbsi + 1
        indri = 3*numat + nbsi
!
!  Compute breathing shell spring constant contribution
!
        do m = 1,npote
          if (nptype(m).eq.14) then
            if (lmatch(nati,ntypi,nspec1(m),nptyp1(m),.true.)) then
              rdiff = radi - twopot(2,m)
              apt = twopot(1,m)*oci
              d2cell(indri,indri,nd2central) = d2cell(indri,indri,nd2central) + apt 
            endif
          elseif (nptype(m).eq.17) then
            if (lmatch(nati,ntypi,nspec1(m),nptyp1(m),.true.)) then
              rdiff = radi - twopot(3,m)
              apt = twopot(1,m)*oci
              bpt = twopot(2,m)
              etrm1 = exp(bpt*rdiff)
              etrm2 = 1.0_dp/etrm1
              etrm = apt*(etrm1 + etrm2)
              d2cell(indri,indri,nd2central) = d2cell(indri,indri,nd2central) + bpt*bpt*etrm
            endif
          elseif (nptype(m).eq.31) then
            if (lmatch(nati,ntypi,nspec1(m),nptyp1(m),.true.)) then
              rdiff = radi - twopot(3,m)
              apt = twopot(1,m)*oci
              bpt = twopot(2,m)
              etrm1 = exp(bpt*rdiff)
              etrm = apt*etrm1
              d2cell(indri,indri,nd2central) = d2cell(indri,indri,nd2central) + bpt*bpt*etrm
            endif
          endif
        enddo
      else
        radi = 0.0_dp
      endif
!
      if (leem.and.lmultiqrange.and.neemrptr(i).ne.0) then
        nqri = nqrnow(neemrptr(i))
      else
        nqri = 1
      endif
!
      indi = 3*(i - 1)
      ix = indi + 1
      iy = indi + 2
      iz = indi + 3
!
!  Molecule handling
!
      if (lmol) then
        nmi = natmol(i)
        indm = nmolind(i)
        call mindtoijk(indm,ixi,iyi,izi)
      endif
!
!  Start of second atom loop
!
      nbsj = 0
      do j = 1,i
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
        xcrd = xclat(j) - xal
        ycrd = yclat(j) - yal
        zcrd = zclat(j) - zal
        qlj = qf(j)
        ocj = occuf(j)
        if (lbsmat(nsft+nrelf2a(j))) then
          radj = radf(j)
          nbsj = nbsj + 1
          indrj = 3*numat + nbsj
        else
          radj = 0.0_dp
        endif
        radsum = radi + radj
        ofct = oci*ocj
        indj = 3*(j - 1)
        jx = indj + 1
        jy = indj + 2
        jz = indj + 3
        fct = ofct*angstoev
        factor = qli*qlj*fct
!
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
!  Molecule handling
!
        if (lmol) then
          nmj = natmol(j)
          indmj = nmolind(j)
          call mindtoijk(indmj,ixj,iyj,izj)
          ixj = ixj - ixi
          iyj = iyj - iyi
          izj = izj - izi
          lmolok = (nmi.eq.nmj.and.nmi.gt.0)
        else
          lmolok = .false.
        endif
!
!  Locate potential number
!  Check whether potential requires specific types
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
              if (mmexc(n).gt.0.or.(lintra(n).and..not.linter(n).or..not.lintra(n).and.linter(n))) lneedmol = .true.
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
!  If no valid potentials and charge product is zero
!  then no need to search for distances
!
        if (npots.eq.0.and.abs(factor).lt.1.0d-8) goto 1120
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
!  If the Coulomb term and the potentials don't need molecule info for
!  for this pair of atoms, then turn off molecule checking to save time
!
        if (.not.lneedmol) lmolok = .false.
!
        d2self = 0.0_dp
!***********************
!  Find valid vectors  *
!***********************
        if (ndim.eq.3) then
          call rsearch3D(xcrd,ycrd,zcrd,lmolok,lcspair,i,j,nmi,ixj,iyj,izj,nor,nmolonly,lself,cut2)
        elseif (ndim.eq.2) then
          call rsearch2D(xcrd,ycrd,zcrd,lmolok,lcspair,i,j,nmi,ixj,iyj,izj,nor,nmolonly,lself,cut2)
        elseif (ndim.eq.1) then
          call rsearch1D(xcrd,ycrd,zcrd,lmolok,lcspair,i,j,nmi,ixj,iyj,izj,nor,nmolonly,lself,cut2)
        endif
!
        derive0self = 0.0_dp
        if (lself) call selftermfc(erecip,ec6,derive0self,factor,fct,ofct,ospfct,1.0_dp,npotl,npots,c6tot, &
                                   d2self,.true.,.true.,i,j,ix,jx,1.0_dp,0.5_dp,lewaldtype)
!
        if (nor.eq.0) goto 1110
!
        do k = 1,nor
          deriv2(k) = 0.0_dp
        enddo
!
!  Sqrt distances
!
        do k = 1,nor
          dist(k) = sqrt(dist(k))
        enddo
!
        call twobody(eatom,ereal,ec6,.true.,.true.,.false.,nor,1_i4,npots,npotl,cut2r,cut2q,cut2s, &
                     nmolonly,factor,ofct,ospfct,radsum,sctrm1,sctrm2,qli,qlj,lcspair,lewaldtype, &
                     .false.,.false.,.false.)
        if (leem.or.lDoQDeriv2) then
          do k = 1,nor
            d0i(k) = d0i(k) + derive0(k)*qlj
            d0j(k) = d0j(k) + derive0(k)*qli
          enddo
        endif
        if (leem) then
          if (lqeq) then
            call qeqbody(eqeq,.true.,.true.,nor,1_i4,fct,qli,qlj,nati,natj,nqri,nqrj)
          elseif (lSandM) then
            call smbody(eqeq,.true.,.true.,nor,1_i4,fct,qli,qlj,nati,natj,nqri,nqrj)
          endif
        endif
!
!  Generate products for derivatives
!
        do k = 1,nor
          rpd(k,1) = xtmp(k)*xtmp(k)
          rpd(k,2) = ytmp(k)*ytmp(k)
          rpd(k,3) = ztmp(k)*ztmp(k)
          rpd(k,4) = ytmp(k)*ztmp(k)
          rpd(k,5) = xtmp(k)*ztmp(k)
          rpd(k,6) = xtmp(k)*ytmp(k)
        enddo
!******************************************
!  Generate Cartesian second derivatives  *
!******************************************
        do k = 1,nor
!
!  Find cell index
!
          if (abs(cellindex(1,k)).gt.nd2cell(1).or. &
              abs(cellindex(2,k)).gt.nd2cell(2).or. &
              abs(cellindex(3,k)).gt.nd2cell(3)) then
!
!  Find cell index - if outside user range then assign to central cell
!
            ncindm = nd2central
            ncindp = nd2central
          else
!
!  Compute index
!
            ncindp = nd2cellptr(nd2cell(1)+1+cellindex(1,k),nd2cell(2)+1+cellindex(2,k),nd2cell(3)+1+cellindex(3,k))
            ncindm = nd2cellptr(nd2cell(1)+1-cellindex(1,k),nd2cell(2)+1-cellindex(2,k),nd2cell(3)+1-cellindex(3,k))
          endif
!
          dk = deriv(k)
          d2k = deriv2(k)
!
          d2cell(jx,ix,ncindp) = d2cell(jx,ix,ncindp) - d2k*rpd(k,1)
          d2cell(jy,ix,ncindp) = d2cell(jy,ix,ncindp) - d2k*rpd(k,6)
          d2cell(jz,ix,ncindp) = d2cell(jz,ix,ncindp) - d2k*rpd(k,5)
          d2cell(jx,iy,ncindp) = d2cell(jx,iy,ncindp) - d2k*rpd(k,6)
          d2cell(jy,iy,ncindp) = d2cell(jy,iy,ncindp) - d2k*rpd(k,2)
          d2cell(jz,iy,ncindp) = d2cell(jz,iy,ncindp) - d2k*rpd(k,4)
          d2cell(jx,iz,ncindp) = d2cell(jx,iz,ncindp) - d2k*rpd(k,5)
          d2cell(jy,iz,ncindp) = d2cell(jy,iz,ncindp) - d2k*rpd(k,4)
          d2cell(jz,iz,ncindp) = d2cell(jz,iz,ncindp) - d2k*rpd(k,3)
          d2cell(jx,ix,ncindp) = d2cell(jx,ix,ncindp) - dk
          d2cell(jy,iy,ncindp) = d2cell(jy,iy,ncindp) - dk
          d2cell(jz,iz,ncindp) = d2cell(jz,iz,ncindp) - dk
!
          if (i.ne.j) then
            d2cell(ix,jx,ncindm) = d2cell(ix,jx,ncindm) - d2k*rpd(k,1)
            d2cell(iy,jx,ncindm) = d2cell(iy,jx,ncindm) - d2k*rpd(k,6)
            d2cell(iz,jx,ncindm) = d2cell(iz,jx,ncindm) - d2k*rpd(k,5)
            d2cell(ix,jy,ncindm) = d2cell(ix,jy,ncindm) - d2k*rpd(k,6)
            d2cell(iy,jy,ncindm) = d2cell(iy,jy,ncindm) - d2k*rpd(k,2)
            d2cell(iz,jy,ncindm) = d2cell(iz,jy,ncindm) - d2k*rpd(k,4)
            d2cell(ix,jz,ncindm) = d2cell(ix,jz,ncindm) - d2k*rpd(k,5)
            d2cell(iy,jz,ncindm) = d2cell(iy,jz,ncindm) - d2k*rpd(k,4)
            d2cell(iz,jz,ncindm) = d2cell(iz,jz,ncindm) - d2k*rpd(k,3)
            d2cell(ix,jx,ncindm) = d2cell(ix,jx,ncindm) - dk
            d2cell(iy,jy,ncindm) = d2cell(iy,jy,ncindm) - dk
            d2cell(iz,jz,ncindm) = d2cell(iz,jz,ncindm) - dk
          endif
        enddo
        if (radsum.gt.0.0_dp.and.i.ne.j) then
!
!  Radial components
!
          do k = 1,nor
!
!  Find cell index
!
            if (abs(cellindex(1,k)).gt.nd2cell(1).or. &
                abs(cellindex(2,k)).gt.nd2cell(2).or. &
                abs(cellindex(3,k)).gt.nd2cell(3)) then
!
!  Find cell index - if outside user range then assign to central cell
!
              ncindm = nd2central
              ncindp = nd2central
            else
!
!  Compute index
!
              ncindp = nd2cellptr(nd2cell(1)+1+cellindex(1,k),nd2cell(2)+1+cellindex(2,k),nd2cell(3)+1+cellindex(3,k))
              ncindm = nd2cellptr(nd2cell(1)+1-cellindex(1,k),nd2cell(2)+1-cellindex(2,k),nd2cell(3)+1-cellindex(3,k))
            endif
            if (radi.gt.0.0_dp) then
              iis = indri
              d2cell(ix,iis,nd2central) = d2cell(ix,iis,nd2central) - rderiv(k)*xtmp(k)
              d2cell(iy,iis,nd2central) = d2cell(iy,iis,nd2central) - rderiv(k)*ytmp(k)
              d2cell(iz,iis,nd2central) = d2cell(iz,iis,nd2central) - rderiv(k)*ztmp(k)
!
              d2cell(iis,ix,nd2central) = d2cell(iis,ix,nd2central) - rderiv(k)*xtmp(k)
              d2cell(iis,iy,nd2central) = d2cell(iis,iy,nd2central) - rderiv(k)*ytmp(k)
              d2cell(iis,iz,nd2central) = d2cell(iis,iz,nd2central) - rderiv(k)*ztmp(k)
!
              d2cell(jx,iis,ncindp) = d2cell(jx,iis,ncindp) + rderiv(k)*xtmp(k)
              d2cell(jy,iis,ncindp) = d2cell(jy,iis,ncindp) + rderiv(k)*ytmp(k)
              d2cell(jz,iis,ncindp) = d2cell(jz,iis,ncindp) + rderiv(k)*ztmp(k)
!
              d2cell(iis,jx,ncindm) = d2cell(iis,jx,ncindm) + rderiv(k)*xtmp(k)
              d2cell(iis,jy,ncindm) = d2cell(iis,jy,ncindm) + rderiv(k)*ytmp(k)
              d2cell(iis,jz,ncindm) = d2cell(iis,jz,ncindm) + rderiv(k)*ztmp(k)
!
              d2cell(iis,iis,nd2central) = d2cell(iis,iis,nd2central) + rtrm2(k)
            endif
            if (radj.gt.0.0_dp) then
              jjs = indrj
              d2cell(jx,jjs,nd2central) = d2cell(jx,jjs,nd2central) + rderiv(k)*xtmp(k)
              d2cell(jy,jjs,nd2central) = d2cell(jy,jjs,nd2central) + rderiv(k)*ytmp(k)
              d2cell(jz,jjs,nd2central) = d2cell(jz,jjs,nd2central) + rderiv(k)*ztmp(k)
!
              d2cell(jjs,jx,nd2central) = d2cell(jjs,jx,nd2central) + rderiv(k)*xtmp(k)
              d2cell(jjs,jy,nd2central) = d2cell(jjs,jy,nd2central) + rderiv(k)*ytmp(k)
              d2cell(jjs,jz,nd2central) = d2cell(jjs,jz,nd2central) + rderiv(k)*ztmp(k)
!
              d2cell(ix,jjs,ncindm) = d2cell(ix,jjs,ncindm) - rderiv(k)*xtmp(k)
              d2cell(iy,jjs,ncindm) = d2cell(iy,jjs,ncindm) - rderiv(k)*ytmp(k)
              d2cell(iz,jjs,ncindm) = d2cell(iz,jjs,ncindm) - rderiv(k)*ztmp(k)
!
              d2cell(jjs,ix,ncindp) = d2cell(jjs,ix,ncindp) - rderiv(k)*xtmp(k)
              d2cell(jjs,iy,ncindp) = d2cell(jjs,iy,ncindp) - rderiv(k)*ytmp(k)
              d2cell(jjs,iz,ncindp) = d2cell(jjs,iz,ncindp) - rderiv(k)*ztmp(k)
!
              d2cell(jjs,jjs,nd2central) = d2cell(jjs,jjs,nd2central) + rtrm2(k)
              if (radi.gt.0.0_dp) then
                d2cell(jjs,iis,ncindp) = d2cell(jjs,iis,ncindp) + rtrm2(k)
                d2cell(iis,jjs,ncindm) = d2cell(iis,jjs,ncindm) + rtrm2(k)
              endif
            endif
          enddo
        endif
!
!  If nor = 0, then rejoin here since there can still be a contribution to 
!  the second derivatives from the variable charges and the self term.
!
1110    continue
!********************************************
!  Variable charge contribution to phonons  *
!********************************************
        if (lDoQDeriv2) then
          hfactor = 1.0_dp
          if (i.eq.j) hfactor = 0.5_dp
          d2self = hfactor*d2self
          derive0selfi = hfactor*derive0self*qlj
          derive0selfj = hfactor*derive0self*qli
          d1ix = 0.0_dp
          d1iy = 0.0_dp
          d1iz = 0.0_dp
          d1jx = 0.0_dp
          d1jy = 0.0_dp
          d1jz = 0.0_dp
          do k = 1,nor
            d0i(k) = d0i(k)*hfactor
            d0j(k) = d0j(k)*hfactor
            d1i(k) = d1i(k)*hfactor
            d1j(k) = d1j(k)*hfactor
            d2i2(k) = d2i2(k)*hfactor
            d2ij(k) = d2ij(k)*hfactor
            d2j2(k) = d2j2(k)*hfactor
            d1ix = d1ix + d1i(k)*xtmp(k)
            d1iy = d1iy + d1i(k)*ytmp(k)
            d1iz = d1iz + d1i(k)*ztmp(k)
            d1jx = d1jx + d1j(k)*xtmp(k)
            d1jy = d1jy + d1j(k)*ytmp(k)
            d1jz = d1jz + d1j(k)*ztmp(k)
          enddo
!
!  Apply variable charge correction to second derivatives
!
          call d2chargefc(i,j,nor,ix,iy,iz,jx,jy,jz,d0i,d0j,d1ix,d1iy,d1iz, &
                         d1jx,d1jy,d1jz,d2i2,d2ij,d2j2,d2self,derive0selfi,derive0selfj,.true.)
        endif
1120    continue
      enddo
    enddo
!
!  End of real space part - perform general tasks
!
  endif
!
!  Free local memory
!
  deallocate(npotl,stat=status)
  if (status/=0) call deallocate_error('realfc','npotl')
!
!  Timing
!
  time2 = g_cpu_time()
  tatom = tatom + time2 - time1
#ifdef TRACE
  call trace_out('realfc')
#endif
!
  return
  end
