  subroutine real1fc(maxlhs,d1cell,matom,vmatom)
!
!  Routine for calculating the first derivatives for
!  use in a finite difference unphased phonon calculation.
!
!  11/14 Created from realfc
!  12/14 Displaced atom coordinates/radius passed in separately
!        Strategy now finds atoms based on unperturbed coordinates
!        to ensure cell is correct. The coordinates and distances
!        are then modified afterwards. 
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
!  Passed variables
!
  integer(i4),                     intent(in)  :: maxlhs
  integer(i4),                     intent(in)  :: matom
  real(dp),                        intent(out) :: d1cell(4,maxlhs,*)
  real(dp),                        intent(in)  :: vmatom(4)
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: indm
  integer(i4)                                  :: indmj
  integer(i4)                                  :: ixi
  integer(i4)                                  :: iyi
  integer(i4)                                  :: izi
  integer(i4)                                  :: ixj
  integer(i4)                                  :: iyj
  integer(i4)                                  :: izj
  integer(i4)                                  :: j
  integer(i4)                                  :: k
  integer(i4)                                  :: m
  integer(i4)                                  :: maxlim
  integer(i4)                                  :: mint
  integer(i4)                                  :: n
  integer(i4)                                  :: nat1
  integer(i4)                                  :: nat2
  integer(i4)                                  :: nati
  integer(i4)                                  :: natj
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
  real(dp)                                     :: d2self
  real(dp)                                     :: derive0self
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
  real(dp)                                     :: xshift
  real(dp)                                     :: yshift
  real(dp)                                     :: zshift
#ifdef TRACE
  call trace_in('real1fc')
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
  if (status/=0) call outofmemory('real1fc','npotl')
!
  if (.not.lnoreal) then
!***************************************************************
!  Atomistic and real space electrostatic second derivatives   *
!***************************************************************
!
!  Only do i = matom case
!
    i = matom
!
!  Inner loop over second site
!
    xal = xclat(i)
    yal = yclat(i)
    zal = zclat(i)
    xshift = vmatom(1) - xclat(i)
    yshift = vmatom(2) - yclat(i)
    zshift = vmatom(3) - zclat(i)
    nati = nat(i)
    ntypi = nftype(i)
    qli = qf(i)
    oci = occuf(i)
    if (lbsmat(nsft+nrelf2a(i))) then
      radi = vmatom(4)
!
!  Compute breathing shell spring constant contribution
!
      do m = 1,npote
        if (nptype(m).eq.14) then
          if (lmatch(nati,ntypi,nspec1(m),nptyp1(m),.true.)) then
            rdiff = radi - twopot(2,m)
            apt = twopot(1,m)*oci
            d1cell(4,i,nd2central) = d1cell(4,i,nd2central) + apt*rdiff
          endif
        elseif (nptype(m).eq.17) then
          if (lmatch(nati,ntypi,nspec1(m),nptyp1(m),.true.)) then
            rdiff = radi - twopot(3,m)
            apt = twopot(1,m)*oci
            bpt = twopot(2,m)
            etrm1 = exp(bpt*rdiff)
            etrm2 = 1.0_dp/etrm1
            etrm = apt*(etrm1 + etrm2)
            d1cell(4,i,nd2central) = d1cell(4,i,nd2central) + apt*bpt*(etrm1 - etrm2)
          endif
        elseif (nptype(m).eq.31) then
          if (lmatch(nati,ntypi,nspec1(m),nptyp1(m),.true.)) then
            rdiff = radi - twopot(3,m)
            apt = twopot(1,m)*oci
            bpt = twopot(2,m)
            etrm1 = exp(bpt*rdiff)
            etrm = apt*etrm1
            d1cell(4,i,nd2central) = d1cell(4,i,nd2central) + apt*bpt*etrm1
          endif
        endif
      enddo
    else
      radi = 0.0_dp
    endif
    if (leem.and.lmultiqrange.and.neemrptr(i).ne.0) then
      nqri = nqrnow(neemrptr(i))
    else
      nqri = 1
    endif
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
    do j = 1,numat
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
      else
        radj = 0.0_dp
      endif
      radsum = radi + radj
      ofct = oci*ocj
      fct = ofct*angstoev
      factor = qli*qlj*fct
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
      if (npots.eq.0.and.abs(factor).lt.1.0d-8) goto 1110
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
                                 d2self,.true.,.false.,i,j,1_i4,1_i4,1.0_dp,0.5_dp,lewaldtype)
!
      if (nor.eq.0) goto 1110
!
!  Correct distances for displacement
!
      do k = 1,nor
        xtmp(k) = xtmp(k) - xshift
        ytmp(k) = ytmp(k) - yshift
        ztmp(k) = ztmp(k) - zshift
        dist(k) = xtmp(k)*xtmp(k) + ytmp(k)*ytmp(k) + ztmp(k)*ztmp(k)
      enddo
!
!  Sqrt distances
!
      do k = 1,nor
        dist(k) = sqrt(dist(k))
      enddo
!
      call twobody(eatom,ereal,ec6,.true.,.false.,.false.,nor,1_i4,npots,npotl,cut2r,cut2q,cut2s, &
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
          call qeqbody(eqeq,.true.,.false.,nor,1_i4,fct,qli,qlj,nati,natj,nqri,nqrj)
        elseif (lSandM) then
          call smbody(eqeq,.true.,.false.,nor,1_i4,fct,qli,qlj,nati,natj,nqri,nqrj)
        endif
      endif
!*****************************************
!  Generate Cartesian first derivatives  *
!*****************************************
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
!
        d1cell(1,j,ncindp) = d1cell(1,j,ncindp) + dk*xtmp(k)
        d1cell(2,j,ncindp) = d1cell(2,j,ncindp) + dk*ytmp(k)
        d1cell(3,j,ncindp) = d1cell(3,j,ncindp) + dk*ztmp(k)
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
          if (radi.gt.0.0_dp.or.radj.gt.0.0_dp) then
            d1cell(4,j,ncindp) = d1cell(4,j,ncindp) + rtrm1(k)
          endif
          if (radi.gt.0.0_dp) then
            d1cell(4,i,nd2central) = d1cell(4,i,nd2central) + rtrm1(k)
          endif
        enddo
      endif
1110  continue
    enddo
!
!  End of real space part - perform general tasks
!
  endif
!
!  Free local memory
!
  deallocate(npotl,stat=status)
  if (status/=0) call deallocate_error('real1fc','npotl')
!
!  Timing
!
  time2 = g_cpu_time()
  tatom = tatom + time2 - time1
#ifdef TRACE
  call trace_out('real1fc')
#endif
!
  return
  end
