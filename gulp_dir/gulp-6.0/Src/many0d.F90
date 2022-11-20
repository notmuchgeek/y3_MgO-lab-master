  subroutine many0d(emany,lgrad1,lgrad2)
!
!  Subroutine for calculating the many-body energy from the
!  Sutton-Chen potential(s). Requires real space routine to
!  be called first to evaluate the repulsive contribution
!  and the rho parameters.
!
!  Finite cluster only version
!
!  On entry the array scrho must contain the density at each atomic site.
!
!  lfreeze = if .true. only do derivatives for atoms with a degree of freedom
!
!  10/18 Created from many - new version that allows for screening
!  11/18 Non-radial arrays removed since they are no longer needed
!  11/18 Cutoff check added for baskes
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!  12/19 Rigid molecule modifications added
!   7/20 scrho now passed to meamfnderv as frho
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
!  Julian Gale, CIC, Curtin University, July 2020
!
  use configurations, only : nregionno, nregiontype, QMMMmode
  use control
  use current
  use derivatives
  use eam
  use energies,       only : siteenergy
  use general
  use iochannels,     only : ioout
  use m_strain,       only : real1strterm
  use mdlogic
  use optimisation
  use progress,       only : lduring_opt
  use sutton
  use symmetry
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  use two
  use vectors,        only : vector_pair
  implicit none
!
!  Passed variables
!
  real(dp),    intent(inout)                   :: emany
  logical,     intent(in)                      :: lgrad1
  logical,     intent(in)                      :: lgrad2
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ind
  integer(i4)                                  :: indij
  integer(i4)                                  :: indik
  integer(i4)                                  :: indjk
  integer(i4)                                  :: ix
  integer(i4)                                  :: iy
  integer(i4)                                  :: iz
  integer(i4)                                  :: ixc
  integer(i4)                                  :: iyc
  integer(i4)                                  :: izc
  integer(i4)                                  :: j
  integer(i4)                                  :: jx
  integer(i4)                                  :: jy
  integer(i4)                                  :: jz
  integer(i4)                                  :: jxc
  integer(i4)                                  :: jyc
  integer(i4)                                  :: jzc
  integer(i4)                                  :: k
  integer(i4)                                  :: kl
  integer(i4)                                  :: km
  integer(i4)                                  :: l
  integer(i4)                                  :: lvec
  integer(i4)                                  :: lxc
  integer(i4)                                  :: lyc
  integer(i4)                                  :: lzc
  integer(i4)                                  :: mp
  integer(i4)                                  :: n
  integer(i4)                                  :: nat1
  integer(i4)                                  :: nat2
  integer(i4)                                  :: nati
  integer(i4)                                  :: natj
  integer(i4)                                  :: natk
  integer(i4)                                  :: neamfnspec2
  integer(i4)                                  :: neamspeci
  integer(i4)                                  :: neamspecj
  integer(i4)                                  :: neamspeck
  integer(i4)                                  :: neamspecl
  integer(i4)                                  :: np1
  integer(i4)                                  :: npartial
  integer(i4)                                  :: npartialik
  integer(i4)                                  :: npartialjk
  integer(i4)                                  :: npot
  integer(i4)                                  :: npotik
  integer(i4)                                  :: npotjk
  integer(i4)                                  :: npots
  integer(i4)                                  :: nregioni
  integer(i4)                                  :: nregionj
  integer(i4)                                  :: nregionl
  integer(i4)                                  :: nregiontypi
  integer(i4)                                  :: nregiontypj
  integer(i4)                                  :: ntyp1  
  integer(i4)                                  :: ntyp2
  integer(i4)                                  :: ntypi
  integer(i4)                                  :: ntypj    
  integer(i4)                                  :: ntypk
  integer(i4), dimension(:), allocatable       :: npotikptr
  integer(i4), dimension(:), allocatable       :: npotjkptr
  integer(i4), dimension(:), allocatable       :: npotl
  integer(i4)                                  :: nvec
  integer(i4)                                  :: nvec0
  integer(i4)                                  :: status
  logical                                      :: lanyvalidik
  logical                                      :: lanyvalidjk
  logical,                                save :: lfirstcall = .true.
  logical                                      :: lnonzeroSij
  logical                                      :: lnonzeroSik
  logical                                      :: lnonzeroSjk
  logical                                      :: lopi
  logical                                      :: lopj
  logical                                      :: lopl
  logical                                      :: lpartial
  logical                                      :: lQMMMok
  logical                                      :: lvalidij
  logical                                      :: lvalidik
  logical                                      :: lvalidjk
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: cut2
  real(dp)                                     :: cut2Xscale
  real(dp)                                     :: cut2k
  real(dp)                                     :: cut2p
  real(dp)                                     :: cut2r
  real(dp)                                     :: cut2rk
  real(dp)                                     :: deriv(3)
  real(dp)                                     :: deriv2(6)
  real(dp)                                     :: drhoij(3,maxmeamcomponent)
  real(dp)                                     :: drhoijs(6,maxmeamcomponent)
  real(dp)                                     :: drhoij2(6,maxmeamcomponent)
  real(dp)                                     :: drhoij2s(21,maxmeamcomponent)
  real(dp)                                     :: drhoij2m(6,3,maxmeamcomponent)
  real(dp)                                     :: drhoik(3,maxmeamcomponent)
  real(dp)                                     :: drhoiks(6,maxmeamcomponent)
  real(dp)                                     :: drhoik2(6,maxmeamcomponent)
  real(dp)                                     :: drhoik2s(21,maxmeamcomponent)
  real(dp)                                     :: drhoik2m(6,3,maxmeamcomponent)
  real(dp)                                     :: drhoji(3,maxmeamcomponent)
  real(dp)                                     :: drhojis(6,maxmeamcomponent)
  real(dp)                                     :: drhoji2(6,maxmeamcomponent)
  real(dp)                                     :: drhoji2s(21,maxmeamcomponent)
  real(dp)                                     :: drhoji2m(6,3,maxmeamcomponent)
  real(dp)                                     :: drhojk(3,maxmeamcomponent)
  real(dp)                                     :: drhojks(6,maxmeamcomponent)
  real(dp)                                     :: drhojk2(6,maxmeamcomponent)
  real(dp)                                     :: drhojk2s(21,maxmeamcomponent)
  real(dp)                                     :: drhojk2m(6,3,maxmeamcomponent)
  real(dp)                                     :: drhoki(3,maxmeamcomponent)
  real(dp)                                     :: drhokis(6,maxmeamcomponent)
  real(dp)                                     :: drhoki2(6,maxmeamcomponent)
  real(dp)                                     :: drhoki2s(21,maxmeamcomponent)
  real(dp)                                     :: drhoki2m(6,3,maxmeamcomponent)
  real(dp)                                     :: drhokj(3,maxmeamcomponent)
  real(dp)                                     :: drhokjs(6,maxmeamcomponent)
  real(dp)                                     :: drhokj2(6,maxmeamcomponent)
  real(dp)                                     :: drhokj2s(21,maxmeamcomponent)
  real(dp)                                     :: drhokj2m(6,3,maxmeamcomponent)
  real(dp)                                     :: drhototij(3)
  real(dp)                                     :: drhototijs(6)
  real(dp)                                     :: drhototij2(6)
  real(dp)                                     :: drhototij2s(21)
  real(dp)                                     :: drhototij2m(6,3)
  real(dp)                                     :: drhototij3(10)
  real(dp)                                     :: drhototik(3)
  real(dp)                                     :: drhototiksum(3)
  real(dp)                                     :: drhototiks(6)
  real(dp)                                     :: drhototik2(6)
  real(dp)                                     :: drhototik2s(21)
  real(dp)                                     :: drhototik2m(6,3)
  real(dp)                                     :: drhototik3(10)
  real(dp)                                     :: drhototji(3)
  real(dp)                                     :: drhototjis(6)
  real(dp)                                     :: drhototji2(6)
  real(dp)                                     :: drhototji2s(21)
  real(dp)                                     :: drhototji2m(6,3)
  real(dp)                                     :: drhototji3(10)
  real(dp)                                     :: drhototjk(3)
  real(dp)                                     :: drhototjksum(3)
  real(dp)                                     :: drhototjks(6)
  real(dp)                                     :: drhototjk2(6)
  real(dp)                                     :: drhototjk2s(21)
  real(dp)                                     :: drhototjk2m(6,3)
  real(dp)                                     :: drhototjk3(10)
  real(dp)                                     :: drhototki(3)
  real(dp)                                     :: drhototkis(6)
  real(dp)                                     :: drhototki2(6)
  real(dp)                                     :: drhototki2s(21)
  real(dp)                                     :: drhototki2m(6,3)
  real(dp)                                     :: drhototki3(10)
  real(dp)                                     :: drhototkj(3)
  real(dp)                                     :: drhototkjs(6)
  real(dp)                                     :: drhototkj2(6)
  real(dp)                                     :: drhototkj2s(21)
  real(dp)                                     :: drhototkj2m(6,3)
  real(dp)                                     :: drhototkj3(10)
  real(dp)                                     :: drhototijk2(3,3)
  real(dp)                                     :: drhototijk2sum(3,3)
  real(dp)                                     :: drhototjik2(3,3)
  real(dp)                                     :: drhototjik2sum(3,3)
  real(dp)                                     :: drhototkij2(3,3)
  real(dp)                                     :: drhototijk2s(6,6)
  real(dp)                                     :: drhototjik2s(6,6)
  real(dp)                                     :: drhototkij2s(6,6)
  real(dp)                                     :: drhototijk2m(6,3)
  real(dp)                                     :: drhototjik2m(6,3)
  real(dp)                                     :: drhototkij2m(6,3)
  real(dp)                                     :: dr2ds(6)
  real(dp)                                     :: d2r2dx2(6)
  real(dp)                                     :: d2r2ds2(6,6)
  real(dp)                                     :: d2r2dsdx(6,3)
  real(dp)                                     :: dt1
  real(dp)                                     :: dt2
  real(dp)                                     :: ebas
  real(dp)                                     :: d1bas
  real(dp)                                     :: d2bas
  real(dp)                                     :: eeam
  real(dp)                                     :: emanytrm
  real(dp)                                     :: oci      
  real(dp)                                     :: ocj  
  real(dp)                                     :: ock
  real(dp)                                     :: ofct
  real(dp)                                     :: ofctijk
  real(dp)                                     :: r
  real(dp)                                     :: r2
  real(dp)                                     :: rcut2
  real(dp)                                     :: rcutfactor
  real(dp)                                     :: rhoi
  real(dp)                                     :: rhoj
  real(dp)                                     :: rhok
  real(dp)                                     :: rhoij(maxmeamcomponent)
  real(dp)                                     :: rhoji(maxmeamcomponent)
  real(dp)                                     :: rhoik(maxmeamcomponent)
  real(dp)                                     :: rhoki(maxmeamcomponent)
  real(dp)                                     :: rhojk(maxmeamcomponent)
  real(dp)                                     :: rhokj(maxmeamcomponent)
  real(dp)                                     :: rik
  real(dp)                                     :: rjk
  real(dp)                                     :: rik2
  real(dp)                                     :: rjk2
  real(dp)                                     :: rk
  real(dp)                                     :: rp
  real(dp)                                     :: rpik
  real(dp)                                     :: rpjk
  real(dp)                                     :: rpijk
  real(dp)                                     :: rscrhoi
  real(dp)                                     :: rscrhoi3
  real(dp)                                     :: rscrhoi5
  real(dp)                                     :: rscrhoj
  real(dp)                                     :: rscrhoj3
  real(dp)                                     :: rscrhoj5
  real(dp)                                     :: rscrhok
  real(dp)                                     :: rscrhok3
  real(dp)                                     :: rscrhok5
  real(dp)                                     :: scmax
  real(dp)                                     :: Sij
  real(dp)                                     :: Sik
  real(dp)                                     :: Sjk
  real(dp)                                     :: Silj
  real(dp)                                     :: Silk
  real(dp)                                     :: Sjlk
  real(dp)                                     :: dSiljdr(3)
  real(dp)                                     :: dSilkdr(3)
  real(dp)                                     :: dSjlkdr(3)
  real(dp)                                     :: time1
  real(dp)                                     :: time2  
  real(dp)                                     :: xal 
  real(dp)                                     :: yal    
  real(dp)                                     :: zal
  real(dp)                                     :: xcd 
  real(dp)                                     :: ycd    
  real(dp)                                     :: zcd
  real(dp)                                     :: xcd1
  real(dp)                                     :: ycd1   
  real(dp)                                     :: zcd1
  real(dp)                                     :: xcd2
  real(dp)                                     :: ycd2   
  real(dp)                                     :: zcd2
  real(dp)                                     :: xil0
  real(dp)                                     :: yil0
  real(dp)                                     :: zil0
  real(dp)                                     :: xjl0
  real(dp)                                     :: yjl0
  real(dp)                                     :: zjl0
  real(dp)                                     :: xkl0
  real(dp)                                     :: ykl0
  real(dp)                                     :: zkl0
  type(screening_atoms)                        :: partial
  type(screening_atoms)                        :: partialik
  type(screening_atoms)                        :: partialjk
  type(vector_pair),                      save :: vectorpair
#ifdef TRACE
  call trace_in('many0d')
#endif
!
  time1 = g_cpu_time()
!
!  If this is the first call then initialise vectorpair
!
  if (lfirstcall) then
    lfirstcall = .false.
    call changemaxvectorpair(vectorpair,0_i4)
  endif
!
!  For screened MEAM, set scale factor that determines searching cutoff based on which axis of the ellipse is largest.
!  Note: the cutoff is applied to the mid point of the i-j vector and so this is guaranteed to find all distances.
!
  rcutfactor = 1.0_dp
  if (lanyMEAMscreen) then
    neamfnspec2 = neamfnspec*(neamfnspec+1)/2
    do i = 1,neamfnspec
      do j = 1,neamfnspec2
        if (lMEAMscreen(j,i)) then
          rcutfactor = max(1.0_dp,0.25_dp*(1.0_dp + meam_Cmax(j,i)))
        endif
      enddo
    enddo
  endif
!
!  Scale density
!
  call eamscalescrho(1_i4)
  if (lPrintEAM) then
!
!  Opening banner for energy decomposition
!
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    write(ioout,'(''  EAM : Atom No.                Density                 Atom energy (eV) '')')
    write(ioout,'(''--------------------------------------------------------------------------------'')')
  endif
!
!  Energy calculation
!
  emany = 0.0_dp
  do i = 1,numat
    neamspeci = neamfnspecptr(i)
    rhoi = scrho(1,i)
    nregioni = nregionno(nsft+nrelf2a(i))
    nregiontypi = nregiontype(nregioni,ncf)
!     
!  QM/MM handling : i is a QM atom => exclude
!     
    lQMMMok = .true.
    if (QMMMmode(ncf).gt.0) then
      if (nregiontypi.eq.1) lQMMMok = .false.
    endif
    if (neamspeci.gt.0.and.rhoi.gt.1.0d-12.and.lQMMMok) then
      if (lMEAMfn) then
        call meamfnderv(neamfn,neamspeci,scrho(1,i),rhoi,eeam,rscrhoi,rscrhoi3,rscrhoi5,.false.,.false.,.false.)
      else
        call eamfnderv(neamfn,neamspeci,scrho(1,i),eeam,rscrhoi,rscrhoi3,rscrhoi5,.false.,.false.,.false.)
      endif
      emanytrm = occuf(i)*eeam
      emany = emany + emanytrm
      siteenergy(i) = siteenergy(i) + emanytrm
!
      if (lPrintEAM) then
        write(ioout,'(7x,i8,4x,f24.10,4x,f24.12)') i,rhoi,emanytrm
      endif
    endif
  enddo
  if (lPrintEAM) then
!
!  Closing banner for energy decomposition
!
    write(ioout,'(''--------------------------------------------------------------------------------'')')
  endif
!
!  If no forces are needed then we don't need to do loops over atoms, so just return
!
  if (.not.lgrad1) goto 1000
!
  if (lduring_opt) then
    cut2Xscale = eamXcutfactor**2
  else
    cut2Xscale = 1.0_dp
  endif
!
!  From here on we can assume that lgrad1 = .true.
!
!***************************
!  Set up local variables  *
!***************************
!
!  Find maximum cut-off radius
!
  scmax = 0.0_dp
  do i = 1,npote
    if (nptype(i).eq.19.or.nptype(i).eq.45.or.nptype(i).eq.55) then
      if (rpot(i).gt.scmax) scmax = rpot(i)
    endif
  enddo
!
!  Allocate local memory
!
  allocate(npotl(npote),stat=status)
  if (status/=0) call outofmemory('many','npotl')
  if (lgrad2) then
    allocate(npotikptr(npote),stat=status)
    if (status/=0) call outofmemory('many','npotikptr')
    allocate(npotjkptr(npote),stat=status)
    if (status/=0) call outofmemory('many','npotjkptr')
  endif
!
!  Set up cutoffs
!
  cut2p = cutp*cutp
  if (lnoreal) goto 999
!
!  Outer loop over sites
!
  ixc = - 2
  iyc = - 1
  izc = 0
  iloop: do i = 1,numat
!
!  Handle ixc-izc before rho check otherwise
!  values go wrong
!
    lopi = (.not.lfreeze.or.lopf(nrelf2a(i)))
    if (lopi) then
      ixc = ixc + 3
      iyc = iyc + 3
      izc = izc + 3
    endif
!
    xal = xclat(i)
    yal = yclat(i)
    zal = zclat(i)
    nati = nat(i)
    ntypi = nftype(i)
    nregioni = nregionno(nsft+nrelf2a(i))
    nregiontypi = nregiontype(nregioni,ncf)
    oci = occuf(i)
!     
!  Find EAM species for i
!  
    neamspeci = neamfnspecptr(i)
!
!  Evaluate functional derivatives
!
    if (lMEAMfn) then
      call meamfnderv(neamfn,neamspeci,scrho(1,i),rhoi,eeam,rscrhoi,rscrhoi3,rscrhoi5,.true.,lgrad2,.false.)
    else
      call eamfnderv(neamfn,neamspeci,scrho(1,i),eeam,rscrhoi,rscrhoi3,rscrhoi5,.true.,lgrad2,.false.)
      rhoi = scrho(1,i)
    endif
!
!  Start of second atom loop
!
    jxc = - 2
    jyc = - 1
    jzc = 0
    jloop: do j = 1,i
!
!  Freezing flag - need to do this before checking rho values otherwise jxc-jzc go wrong
!
      lopj = (.not.lfreeze.or.lopf(nrelf2a(j)))
      if (.not.lopi.and..not.lopj) cycle jloop
      if (lopj) then
        jxc = jxc + 3
        jyc = jyc + 3
        jzc = jzc + 3
      endif
!
      natj = nat(j)
      ntypj = nftype(j)
      nregionj = nregionno(nsft+nrelf2a(j))
      nregiontypj = nregiontype(nregionj,ncf)
!     
!  QM/MM handling : i & j are both QM atoms => no forces to compute
!     
      if (QMMMmode(ncf).gt.0) then
        if (nregiontypi.eq.1.and.nregiontypj.eq.1) cycle jloop
      endif
!     
!  Find EAM species for j
!  
      neamspecj = neamfnspecptr(j)
!
!  Find index for i-j
!
      if (neamspeci.gt.neamspecj) then
        indij = neamspeci*(neamspeci - 1)/2 + neamspecj
      else
        indij = neamspecj*(neamspecj - 1)/2 + neamspeci
      endif
!
!  Evaluate functional derivatives
!
      if (lMEAMfn) then
        call meamfnderv(neamfn,neamspecj,scrho(1,j),rhoj,eeam,rscrhoj,rscrhoj3,rscrhoj5,.true.,lgrad2,.false.)
      else
        call eamfnderv(neamfn,neamspecj,scrho(1,j),eeam,rscrhoj,rscrhoj3,rscrhoj5,.true.,lgrad2,.false.)
        rhoj = scrho(1,j)
      endif
!
!  If there is no density at either centre and this is not a second derivative run then cycle
!
      if (.not.lMEAMden.and.rhoi.eq.0.0_dp.and.rhoj.eq.0.0_dp.and..not.lgrad2) cycle jloop
!
!  Arrange nat1 to be lower than nat2. If they are equal
!  then arrange that ntyp1 is lower than ntyp2
!
      if (nati.eq.natj) then
        nat1 = nati
        nat2 = natj
        if (ntypi.lt.ntypj) then
          ntyp1 = ntypi
          ntyp2 = ntypj
        else
          ntyp1 = ntypj
          ntyp2 = ntypi
        endif
      elseif (nati.lt.natj) then
        nat1 = nati
        nat2 = nat(j)
        ntyp1 = ntypi
        ntyp2 = nftype(j)
      else
        nat1 = nat(j)
        nat2 = nati
        ntyp1 = nftype(j)
        ntyp2 = ntypi
      endif
      xcd = xclat(j) - xal
      ycd = yclat(j) - yal
      zcd = zclat(j) - zal
      ocj = occuf(j)
!
!  Set flags such that if one atom is not being optimised place 
!  3 x 3 second derivative matrix in the on-diagonal block
!
      if (lopi.and.lopj) then
        ix = ixc
        iy = iyc
        iz = izc
        jx = jxc
        jy = jyc
        jz = jzc
      elseif (lopi) then
        ix = ixc
        iy = iyc
        iz = izc
        jx = ixc
        jy = iyc
        jz = izc
      elseif (lopj) then
        ix = jxc
        iy = jyc
        iz = jzc
        jx = jxc
        jy = jyc
        jz = jzc
      endif
      ofct = oci*ocj
!
!  Locate potential number
!  Check whether potential requires specific types
!  Calculate sum of all dispersion terms for pair
!
      npots = 0
      rp = 0.0_dp
      do n = 1,npote
        if (nptype(n).eq.19.or.nptype(n).eq.45.or.nptype(n).eq.55) then
          if (nat1.eq.nspec1(n).and.nat2.eq.nspec2(n)) then
            if ((ntyp1.eq.nptyp1(n).or.nptyp1(n).eq.0).and.(ntyp2.eq.nptyp2(n).or.nptyp2(n).eq.0)) then
              npots = npots + 1
              npotl(npots) = n
              if (rpot(n).gt.rp) rp = rpot(n)
            endif
          endif
        endif
      enddo
!
!  If no valid potentials then there is no need
!  to continue with this pair, unless this is
!  a second derivative calculation, in which
!  case there may be a contribution from triangles
!  of interactions.
!
      if (npots.eq.0.and..not.lgrad2) cycle jloop
      if (lgrad2) then
!
!  Need to make cut-off equal to double the maximum
!  to ensure all triangles are included
!
        rp = 2.0_dp*scmax
        cut2r = rp*rp
        if (cut2r.gt.4.0_dp*cut2p) cut2r = cut2p
      else
        cut2r = rp*rp
        if (cut2r.gt.cut2p) cut2r = cut2p
      endif
      cut2 = cut2r
      rp = sqrt(cut2)
      r2 = xcd*xcd + ycd*ycd + zcd*zcd
      if (r2.gt.smallself.and.r2.le.cut2) then
!***************************************************
!  Calculate many-body contribution in real space  *
!***************************************************
        deriv(1:3) = 0.0_dp
        if (lgrad2) then
          deriv2(1:6) = 0.0_dp
        endif
        r = sqrt(r2)
!***************************************
!  Valid many-body potentials for i-j  *
!***************************************
        if (lMEAM) then
          rhoij(1:maxmeamcomponent) = 0.0_dp
          rhoji(1:maxmeamcomponent) = 0.0_dp
          drhoij(1:3,1:maxmeamcomponent) = 0.0_dp
          drhoji(1:3,1:maxmeamcomponent) = 0.0_dp
          if (lgrad2) then
            drhoij2(1:6,1:maxmeamcomponent) = 0.0_dp
            drhoji2(1:6,1:maxmeamcomponent) = 0.0_dp
          endif
        else
          rhoij(1) = 0.0_dp
          rhoji(1) = 0.0_dp
        endif
        drhototij(1:3) = 0.0_dp
        drhototji(1:3) = 0.0_dp
        if (lgrad2) then
          drhototij2(1:6) = 0.0_dp
          drhototji2(1:6) = 0.0_dp
        endif
        lvalidij = .false.
        if (npots.gt.0) then
          if (lMEAMden) then
            do mp = 1,npots
              npot = npotl(mp)
              if (nptype(npot).eq.19) then
                if (r.gt.rpot2(npot).and.r.le.rpot(npot).and.r.le.rpmax) then
                  lvalidij = .true.
!**********************************
!  Calculate density derivatives  *
!**********************************
                  call meamrho(nati,ntypi,natj,ntypj,r,rpot(npot),xcd,ycd,zcd,rhoij,rhoji,drhoij,drhoji, &
                               drhoijs,drhojis,drhoij2,drhoji2,drhoij2s,drhoji2s,drhoij2m,drhoji2m, &
                               1.0_dp,1.0_dp,.true.,.false.,.true.,lgrad2,twopot(1,npot))
                endif
              endif
            enddo
!*******************************
!  Compute screening function  *
!*******************************
!
!  Initialise screening function to 1
!
            lnonzeroSij = .true.
            Sij = 1.0_dp
!
            if (lanyMEAMscreen) then
!
!  Set cutoffs for possible screening atoms
!
              rcut2 = rcutfactor*r2
!
!  Loop over atoms to search for images that may contribute to the screening
!
              l = 0
              npartial = 0
              lxc = - 2
              lyc = - 1
              lzc =   0
              do while (l.lt.numat.and.lnonzeroSij)
                l = l + 1
                neamspecl = neamfnspecptr(l)
                lopl = (.not.lfreeze.or.lopf(nrelf2a(l)))
                if (lopl) then
                  lxc = lxc + 3
                  lyc = lyc + 3
                  lzc = lzc + 3
                endif
                if (lMEAMscreen(indij,neamspecl)) then
!
!  Set basic vectors between atoms
!
                  xil0 = xclat(l) - xal
                  yil0 = yclat(l) - yal
                  zil0 = zclat(l) - zal
                  xjl0 = xil0 - xcd
                  yjl0 = yil0 - ycd
                  zjl0 = zil0 - zcd
!
!  Find images within cutoffs of both atoms - excluding self images
!
                  nvec0 = 0
                  call rfindmid(xil0,yil0,zil0,xjl0,yjl0,zjl0,rcut2,.false.,nvec0,nvec,vectorpair)
!
!  Loop over results of search
!
                  lvec = 0
                  do while (lvec.lt.nvec.and.lnonzeroSij)
                    lvec = lvec + 1
!
!  Compute screening function
!
                    call meamscreen(neamspeci,neamspecj,neamspecl,r2,vectorpair%distance_pair1(lvec), &
                                    vectorpair%distance_pair2(lvec),Silj,dSiljdr,lpartial,.true.)
!
!  If screening function contribution is 0, then no need to continue for this pair
!
                    if (Silj.eq.0.0_dp) then
                      lnonzeroSij = .false.
                      Sij = 0.0_dp
                    else
!
!  Multiply total screening product
!
                      Sij = Sij*Silj
                      if (lpartial) then
!
!  If this atom has a screening factor between 0 and 1, we need to keep track of it since it will generate non-zero derivatives
!
                        npartial = npartial + 1
                        if (npartial.gt.partial%sa_maxdim) then
                          call changemaxsa(partial,npartial)
                        endif
                        partial%sa_atom(npartial) = l
                        partial%sa_kxc(npartial) = lxc
                        partial%sa_xik(npartial) = vectorpair%xvector_pair1(lvec)
                        partial%sa_yik(npartial) = vectorpair%yvector_pair1(lvec)
                        partial%sa_zik(npartial) = vectorpair%zvector_pair1(lvec)
                        partial%sa_xjk(npartial) = vectorpair%xvector_pair2(lvec)
                        partial%sa_yjk(npartial) = vectorpair%yvector_pair2(lvec)
                        partial%sa_zjk(npartial) = vectorpair%zvector_pair2(lvec)
                        partial%sa_Sikj(npartial) = Silj
                        partial%sa_dSikjdr(1:3,npartial)   = dSiljdr(1:3)
                      endif
                    endif
!
!  End loop over images of possible screening atoms
!
                  enddo
                endif
!
!  End loop over possible screening atoms
!
              enddo
!
!  End of screening function
!
            endif
!
!  Only do remainder of work if the screening factor is non-zero
!
            if (lnonzeroSij) then
              do mp = 1,npots
                npot = npotl(mp)
!
!  Pair potential contribution
!
                if (nptype(npot).eq.45.or.nptype(npot).eq.55) then
                  if (r.gt.rpot2(npot).and.r.le.rpot(npot)) then
                    lvalidij = .true.
                    rk = 1.0_dp/r
                    call baskes(nati,ntypi,neamspeci,natj,ntypj,neamspecj,npot,r,rk,1.0_dp,ebas,d1bas,d2bas,.true.,lgrad2)
!
                    call real1strterm(ndim,xcd,ycd,zcd,0.0_dp,0.0_dp,0.0_dp,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,lgrad2)
!
                    ebas  = ebas*ofct
                    d1bas = rk*d1bas*ofct
                    d2bas = d2bas*ofct
                    d2bas = rk*rk*(d2bas - d1bas)
!
                    deriv(1) = deriv(1) + d1bas*xcd*Sij
                    deriv(2) = deriv(2) + d1bas*ycd*Sij
                    deriv(3) = deriv(3) + d1bas*zcd*Sij
!
                    if (lgrad2) then
                      deriv2(1) = deriv2(1) + d2bas*Sij*d2r2dx2(1) + d1bas*Sij
                      deriv2(2) = deriv2(2) + d2bas*Sij*d2r2dx2(6)
                      deriv2(3) = deriv2(3) + d2bas*Sij*d2r2dx2(5)
                      deriv2(4) = deriv2(4) + d2bas*Sij*d2r2dx2(2) + d1bas*Sij
                      deriv2(5) = deriv2(5) + d2bas*Sij*d2r2dx2(4)
                      deriv2(6) = deriv2(6) + d2bas*Sij*d2r2dx2(3) + d1bas*Sij
                    endif
!
                    call psiscreenderv(i,j,npartial,partial,xcd,ycd,zcd,ebas,Sij,.false.)
                  endif
                endif
!
!  Density contribution
!
                if (nptype(npot).eq.19) then
                  if (lanyMEAMscreen) then
                    if (npartial.gt.0) then
!
!  Compute derivative contributions of the screening function w.r.t. total density of atom i
!
                      call meamtotalrhoscreenderv(neamspeci,npartial,partial,xcd,ycd,zcd,scrho(1,i), &
                                                  rhoi,rscrhoi,rhoij,Sij,.false.)
!
                      do np1 = 1,npartial
                        l = partial%sa_atom(np1)
                        lopl = (.not.lfreeze.or.lopf(nrelf2a(l)))
                        nregionl = nregionno(nsft+nrelf2a(l))
!
!  i-j contribution
!
                        if (lopi) then
                          xdrv(i) = xdrv(i) - partial%sa_drhototij(1,np1)*ofct
                          ydrv(i) = ydrv(i) - partial%sa_drhototij(2,np1)*ofct
                          zdrv(i) = zdrv(i) - partial%sa_drhototij(3,np1)*ofct
                        endif
                        if (lopj) then
                          xdrv(j) = xdrv(j) + partial%sa_drhototij(1,np1)*ofct
                          ydrv(j) = ydrv(j) + partial%sa_drhototij(2,np1)*ofct
                          zdrv(j) = zdrv(j) + partial%sa_drhototij(3,np1)*ofct
                        endif
                        if (nregioni.ne.nregionj) then
                          xregdrv(nregioni) = xregdrv(nregioni) - partial%sa_drhototij(1,np1)*ofct
                          yregdrv(nregioni) = yregdrv(nregioni) - partial%sa_drhototij(2,np1)*ofct
                          zregdrv(nregioni) = zregdrv(nregioni) - partial%sa_drhototij(3,np1)*ofct
                          xregdrv(nregionj) = xregdrv(nregionj) + partial%sa_drhototij(1,np1)*ofct
                          yregdrv(nregionj) = yregdrv(nregionj) + partial%sa_drhototij(2,np1)*ofct
                          zregdrv(nregionj) = zregdrv(nregionj) + partial%sa_drhototij(3,np1)*ofct
                        endif
!
!  i-l contribution
!
                        if (lopi) then
                          xdrv(i) = xdrv(i) - partial%sa_drhototik(1,np1)*ofct
                          ydrv(i) = ydrv(i) - partial%sa_drhototik(2,np1)*ofct
                          zdrv(i) = zdrv(i) - partial%sa_drhototik(3,np1)*ofct
                        endif
                        if (lopl) then
                          xdrv(l) = xdrv(l) + partial%sa_drhototik(1,np1)*ofct
                          ydrv(l) = ydrv(l) + partial%sa_drhototik(2,np1)*ofct
                          zdrv(l) = zdrv(l) + partial%sa_drhototik(3,np1)*ofct
                        endif
                        if (nregioni.ne.nregionl) then
                          xregdrv(nregioni) = xregdrv(nregioni) - partial%sa_drhototik(1,np1)*ofct
                          yregdrv(nregioni) = yregdrv(nregioni) - partial%sa_drhototik(2,np1)*ofct
                          zregdrv(nregioni) = zregdrv(nregioni) - partial%sa_drhototik(3,np1)*ofct
                          xregdrv(nregionl) = xregdrv(nregionl) + partial%sa_drhototik(1,np1)*ofct
                          yregdrv(nregionl) = yregdrv(nregionl) + partial%sa_drhototik(2,np1)*ofct
                          zregdrv(nregionl) = zregdrv(nregionl) + partial%sa_drhototik(3,np1)*ofct
                        endif
!
!  j-l contribution
!
                        if (lopj) then
                          xdrv(j) = xdrv(j) - partial%sa_drhototjk(1,np1)*ofct
                          ydrv(j) = ydrv(j) - partial%sa_drhototjk(2,np1)*ofct
                          zdrv(j) = zdrv(j) - partial%sa_drhototjk(3,np1)*ofct
                        endif
                        if (lopl) then
                          xdrv(l) = xdrv(l) + partial%sa_drhototjk(1,np1)*ofct
                          ydrv(l) = ydrv(l) + partial%sa_drhototjk(2,np1)*ofct
                          zdrv(l) = zdrv(l) + partial%sa_drhototjk(3,np1)*ofct
                        endif
                        if (nregionj.ne.nregionl) then
                          xregdrv(nregionj) = xregdrv(nregionj) - partial%sa_drhototjk(1,np1)*ofct
                          yregdrv(nregionj) = yregdrv(nregionj) - partial%sa_drhototjk(2,np1)*ofct
                          zregdrv(nregionj) = zregdrv(nregionj) - partial%sa_drhototjk(3,np1)*ofct
                          xregdrv(nregionl) = xregdrv(nregionl) + partial%sa_drhototjk(1,np1)*ofct
                          yregdrv(nregionl) = yregdrv(nregionl) + partial%sa_drhototjk(2,np1)*ofct
                          zregdrv(nregionl) = zregdrv(nregionl) + partial%sa_drhototjk(3,np1)*ofct
                        endif
                      enddo
!
!  Compute derivative contributions of the screening function w.r.t. total density of atom j
!
                      call meamtotalrhoscreenderv(neamspecj,npartial,partial,xcd,ycd,zcd,scrho(1,j), &
                                                  rhoj,rscrhoj,rhoji,Sij,.false.)
!
                      do np1 = 1,npartial
                        l = partial%sa_atom(np1)
                        lopl = (.not.lfreeze.or.lopf(nrelf2a(l)))
                        nregionl = nregionno(nsft+nrelf2a(l))
!
!  i-j contribution
!
                        if (lopi) then
                          xdrv(i) = xdrv(i) - partial%sa_drhototij(1,np1)*ofct
                          ydrv(i) = ydrv(i) - partial%sa_drhototij(2,np1)*ofct
                          zdrv(i) = zdrv(i) - partial%sa_drhototij(3,np1)*ofct
                        endif
                        if (lopj) then
                          xdrv(j) = xdrv(j) + partial%sa_drhototij(1,np1)*ofct
                          ydrv(j) = ydrv(j) + partial%sa_drhototij(2,np1)*ofct
                          zdrv(j) = zdrv(j) + partial%sa_drhototij(3,np1)*ofct
                        endif
                        if (nregioni.ne.nregionj) then
                          xregdrv(nregioni) = xregdrv(nregioni) - partial%sa_drhototij(1,np1)*ofct
                          yregdrv(nregioni) = yregdrv(nregioni) - partial%sa_drhototij(2,np1)*ofct
                          zregdrv(nregioni) = zregdrv(nregioni) - partial%sa_drhototij(3,np1)*ofct
                          xregdrv(nregionj) = xregdrv(nregionj) + partial%sa_drhototij(1,np1)*ofct
                          yregdrv(nregionj) = yregdrv(nregionj) + partial%sa_drhototij(2,np1)*ofct
                          zregdrv(nregionj) = zregdrv(nregionj) + partial%sa_drhototij(3,np1)*ofct
                        endif
!
!  i-l contribution
!
                        if (lopi) then
                          xdrv(i) = xdrv(i) - partial%sa_drhototik(1,np1)*ofct
                          ydrv(i) = ydrv(i) - partial%sa_drhototik(2,np1)*ofct
                          zdrv(i) = zdrv(i) - partial%sa_drhototik(3,np1)*ofct
                        endif
                        if (lopl) then
                          xdrv(l) = xdrv(l) + partial%sa_drhototik(1,np1)*ofct
                          ydrv(l) = ydrv(l) + partial%sa_drhototik(2,np1)*ofct
                          zdrv(l) = zdrv(l) + partial%sa_drhototik(3,np1)*ofct
                        endif
                        if (nregioni.ne.nregionl) then
                          xregdrv(nregioni) = xregdrv(nregioni) - partial%sa_drhototik(1,np1)*ofct
                          yregdrv(nregioni) = yregdrv(nregioni) - partial%sa_drhototik(2,np1)*ofct
                          zregdrv(nregioni) = zregdrv(nregioni) - partial%sa_drhototik(3,np1)*ofct
                          xregdrv(nregionl) = xregdrv(nregionl) + partial%sa_drhototik(1,np1)*ofct
                          yregdrv(nregionl) = yregdrv(nregionl) + partial%sa_drhototik(2,np1)*ofct
                          zregdrv(nregionl) = zregdrv(nregionl) + partial%sa_drhototik(3,np1)*ofct
                        endif
!
!  j-l contribution
!
                        if (lopj) then
                          xdrv(j) = xdrv(j) - partial%sa_drhototjk(1,np1)*ofct
                          ydrv(j) = ydrv(j) - partial%sa_drhototjk(2,np1)*ofct
                          zdrv(j) = zdrv(j) - partial%sa_drhototjk(3,np1)*ofct
                        endif
                        if (lopl) then
                          xdrv(l) = xdrv(l) + partial%sa_drhototjk(1,np1)*ofct
                          ydrv(l) = ydrv(l) + partial%sa_drhototjk(2,np1)*ofct
                          zdrv(l) = zdrv(l) + partial%sa_drhototjk(3,np1)*ofct
                        endif
                        if (nregionj.ne.nregionl) then
                          xregdrv(nregionj) = xregdrv(nregionj) - partial%sa_drhototjk(1,np1)*ofct
                          yregdrv(nregionj) = yregdrv(nregionj) - partial%sa_drhototjk(2,np1)*ofct
                          zregdrv(nregionj) = zregdrv(nregionj) - partial%sa_drhototjk(3,np1)*ofct
                          xregdrv(nregionl) = xregdrv(nregionl) + partial%sa_drhototjk(1,np1)*ofct
                          yregdrv(nregionl) = yregdrv(nregionl) + partial%sa_drhototjk(2,np1)*ofct
                          zregdrv(nregionl) = zregdrv(nregionl) + partial%sa_drhototjk(3,np1)*ofct
                        endif
                      enddo
                    endif
                  endif
!
!  Scale density and derivatives by screening factor
!
                  drhoij(1:3,1:maxmeamcomponent) = Sij*drhoij(1:3,1:maxmeamcomponent)
                  drhoji(1:3,1:maxmeamcomponent) = Sij*drhoji(1:3,1:maxmeamcomponent)
                  if (lgrad2) then
                    drhoij2(1:6,1:maxmeamcomponent) = Sij*drhoij2(1:6,1:maxmeamcomponent)
                    drhoji2(1:6,1:maxmeamcomponent) = Sij*drhoji2(1:6,1:maxmeamcomponent)
                  endif
!
!  Compute total derivatives of MEAM density
!
                  call meamtotalrhoderv(neamspeci,scrho(1,i),rhoi,drhoij,drhototij,drhoijs,drhototijs, &
                                        drhoij2,drhototij2,drhoij2s,drhototij2s,drhoij2m,drhototij2m, &
                                        .false.,lgrad2)
                  call meamtotalrhoderv(neamspecj,scrho(1,j),rhoj,drhoji,drhototji,drhojis,drhototjis, &
                                        drhoji2,drhototji2,drhoji2s,drhototji2s,drhoji2m,drhototji2m, &
                                        .false.,lgrad2)
                endif
              enddo
            endif
          else
            do mp = 1,npots
              npot = npotl(mp)
              if (nptype(npot).eq.19) then
                if (r.gt.rpot2(npot).and.r.le.rpot(npot)) then
                  lvalidij = .true.
                  call eamrho(nati,ntypi,natj,ntypj,r,rpot(npot),xcd,ycd,zcd,rhoij,rhoji,drhototij,drhototji, &
                              drhototijs,drhototjis,drhototij2,drhototji2,drhototij2s,drhototji2s, &
                              drhototij2m,drhototji2m,drhototij3,drhototji3,1.0_dp,1.0_dp,.false.,.true.,lgrad2,.false., &
                              twopot(1,npot))
                endif
              endif
            enddo
          endif
!
!  Combine derivative terms
!
          if (QMMMmode(ncf).gt.0) then
            if (nregiontypi.ne.1) then
              deriv(1:3) = deriv(1:3) + rscrhoi*drhototij(1:3)*ofct
            endif
            if (nregiontypj.ne.1) then
              deriv(1:3) = deriv(1:3) + rscrhoj*drhototji(1:3)*ofct
            endif
          else
            deriv(1:3) = deriv(1:3) + (rscrhoi*drhototij(1:3) + rscrhoj*drhototji(1:3))*ofct
          endif
          if (lgrad2) then
            ind = 0
            do kl = 1,3
              do km = kl,3
                ind = ind + 1
                deriv2(ind) = deriv2(ind) + (rscrhoi*drhototij2(ind) + rscrhoj*drhototji2(ind))*ofct
                deriv2(ind) = deriv2(ind) + ocj*rscrhoi3*drhototij(kl)*drhototij(km)*ofct
                deriv2(ind) = deriv2(ind) + oci*rscrhoj3*drhototji(kl)*drhototji(km)*ofct
              enddo
            enddo
          endif
        endif
        if (lvalidij) then
!******************************
!  Internal first derivatives *
!******************************
          if (lopi) then
            xdrv(i) = xdrv(i) - deriv(1)
            ydrv(i) = ydrv(i) - deriv(2)
            zdrv(i) = zdrv(i) - deriv(3)
          endif
          if (lopj) then
            xdrv(j) = xdrv(j) + deriv(1)
            ydrv(j) = ydrv(j) + deriv(2)
            zdrv(j) = zdrv(j) + deriv(3)
          endif
          if (nregioni.ne.nregionj) then
            xregdrv(nregioni) = xregdrv(nregioni) - deriv(1)
            yregdrv(nregioni) = yregdrv(nregioni) - deriv(2)
            zregdrv(nregioni) = zregdrv(nregioni) - deriv(3)
            xregdrv(nregionj) = xregdrv(nregionj) + deriv(1)
            yregdrv(nregionj) = yregdrv(nregionj) + deriv(2)
            zregdrv(nregionj) = zregdrv(nregionj) + deriv(3)
          endif
!********************************
!  Internal second derivatives  *
!********************************
          if (lgrad2) then
            derv2(jx,ix) = derv2(jx,ix) - deriv2(1)
            derv2(jy,ix) = derv2(jy,ix) - deriv2(2)
            derv2(jz,ix) = derv2(jz,ix) - deriv2(3)
            derv2(jx,iy) = derv2(jx,iy) - deriv2(2)
            derv2(jy,iy) = derv2(jy,iy) - deriv2(4)
            derv2(jz,iy) = derv2(jz,iy) - deriv2(5)
            derv2(jx,iz) = derv2(jx,iz) - deriv2(3)
            derv2(jy,iz) = derv2(jy,iz) - deriv2(5)
            derv2(jz,iz) = derv2(jz,iz) - deriv2(6)
          endif
        endif
        if (lgrad2) then
!******************************************************************
!  Start of third atom loop - only needed for second derivatives  *
!******************************************************************
          kloop: do k = 1,numat
!
!  If k  =  i or j then skip
!
            if (k.eq.i.or.k.eq.j) cycle kloop
!
            natk = nat(k)
            ntypk = nftype(k)
            xcd1 = xclat(k) - xal
            ycd1 = yclat(k) - yal
            zcd1 = zclat(k) - zal
            xcd2 = xcd1 - xcd
            ycd2 = ycd1 - ycd
            zcd2 = zcd1 - zcd
            ock = occuf(k)
            ofctijk = oci*ocj*ock
!
!  Check whether there are any potentials between i-k or j-k
!
            npotik = 0
            npotjk = 0
            rpik = 0.0_dp
            rpjk = 0.0_dp
            do n = 1,npote
              if (nptype(n).eq.19) then
                lvalidik = .false.
                if (nati.eq.nspec1(n).and.natk.eq.nspec2(n)) then
                  if (ntypi.eq.nptyp1(n).or.nptyp1(n).eq.0) then
                    if (ntypk.eq.nptyp2(n).or.nptyp2(n).eq.0) lvalidik = .true.
                  endif
                elseif (nati.eq.nspec2(n).and.natk.eq.nspec1(n)) then
                  if (ntypi.eq.nptyp2(n).or.nptyp2(n).eq.0) then
                    if (ntypk.eq.nptyp1(n).or.nptyp1(n).eq.0) lvalidik = .true.
                  endif
                endif
                if (lvalidik) then
                  npotik = npotik + 1
                  npotikptr(npotik) = n
                  if (rpot(n).gt.rpik) rpik = rpot(n)
                endif
!
                lvalidjk = .false.
                if (natj.eq.nspec1(n).and.natk.eq.nspec2(n)) then
                  if (ntypj.eq.nptyp1(n).or.nptyp1(n).eq.0) then
                    if (ntypk.eq.nptyp2(n).or.nptyp2(n).eq.0) lvalidjk = .true.
                  endif
                elseif (natj.eq.nspec2(n).and.natk.eq.nspec1(n)) then
                  if (ntypj.eq.nptyp2(n).or.nptyp2(n).eq.0) then
                    if (ntypk.eq.nptyp1(n).or.nptyp1(n).eq.0) lvalidjk = .true.
                  endif
                endif
                if (lvalidjk) then
                  npotjk = npotjk + 1
                  npotjkptr(npotjk) = n
                  if (rpot(n).gt.rpjk) rpjk = rpot(n)
                endif
              endif
            enddo
            rpijk = max(rpik,rpjk)
!
!  If no valid potentials for i-k or j-k then skip
!
            if ((npotik+npotjk).eq.0) cycle kloop
            cut2rk = rpijk*rpijk
            if (cut2rk.gt.cut2p) cut2rk = cut2p
            cut2k = cut2Xscale*cut2rk
!     
!  Find EAM species for k
!  
            neamspeck = neamfnspecptr(k)
!
!  Find index for i-k
!
            if (neamspeci.gt.neamspeck) then
              indik = neamspeci*(neamspeci - 1)/2 + neamspeck
            else
              indik = neamspeck*(neamspeck - 1)/2 + neamspeci
            endif
!
!  Find index for j-k
!
            if (neamspecj.gt.neamspeck) then
              indjk = neamspecj*(neamspecj - 1)/2 + neamspeck
            else
              indjk = neamspeck*(neamspeck - 1)/2 + neamspecj
            endif
!
!  Evaluate functional derivatives
!
            if (lMEAMfn) then
              call meamfnderv(neamfn,neamspeck,scrho(1,k),rhok,eeam,rscrhok,rscrhok3,rscrhok5,.true.,lgrad2,.false.)
            else
              call eamfnderv(neamfn,neamspeck,scrho(1,k),eeam,rscrhok,rscrhok3,rscrhok5,.true.,lgrad2,.false.)
              rhok = scrho(1,k)
            endif
!
!  If no rho then skip
!
            if (abs(rhoi+rhoj+rhok).eq.0.0_dp) cycle kloop
!
!  Set up constants for k
!
            rik2 = xcd1*xcd1 + ycd1*ycd1 + zcd1*zcd1
            rjk2 = xcd2*xcd2 + ycd2*ycd2 + zcd2*zcd2
!
            drhototiksum(1:3) = 0.0_dp
            drhototjksum(1:3) = 0.0_dp
            drhototijk2sum(1:3,1:3) = 0.0_dp
            drhototjik2sum(1:3,1:3) = 0.0_dp
!************************************************************
!  Calculate triangular contribution to second derivatives  *
!************************************************************
            lanyvalidik = .false.
            lanyvalidjk = .false.
            if (rik2.gt.smallself.and.rik2.le.cut2k.and.npotik.gt.0) then
!*********************
!  i-k contribution  *
!*********************
!
!  Zero terms
!
              if (lMEAM) then
                rhoik(1:maxmeamcomponent) = 0.0_dp
                rhoki(1:maxmeamcomponent) = 0.0_dp
                drhoik(1:3,1:maxmeamcomponent) = 0.0_dp
                drhoki(1:3,1:maxmeamcomponent) = 0.0_dp
              else
                rhoik(1) = 0.0_dp
                rhoki(1) = 0.0_dp
              endif
              drhototik(1:3) = 0.0_dp
              drhototki(1:3) = 0.0_dp
              drhototijk2(1:3,1:3) = 0.0_dp
!
              rik = sqrt(rik2)
!
!  Loop over potentials to find many-body ones
!
              do n = 1,npotik
                mp = npotikptr(n)
                lvalidik = .false.
                if (rik.gt.rpot2(mp).and.rik.le.rpot(mp)) then
                  lvalidik = .true.
                  if (lvalidik) then
!**********************************
!  Calculate density derivatives  *
!**********************************
                    if (lMEAMden) then
                      call meamrho(nati,ntypi,natk,ntypk,rik,rpot(mp),xcd1,ycd1,zcd1,rhoik,rhoki,drhoik,drhoki, &
                                   drhoiks,drhokis,drhoik2,drhoki2,drhoik2s,drhoki2s,drhoik2m,drhoki2m, &
                                   1.0_dp,1.0_dp,.true.,.false.,.true.,.false.,twopot(1,mp))
!*******************************
!  Compute screening function  *
!*******************************
!
!  Initialise screening function to 1
!
                      lnonzeroSik = .true.
                      Sik = 1.0_dp
!
                      if (lanyMEAMscreen) then
!
!  Set cutoffs for possible screening atoms
!
                        rcut2 = rcutfactor*rik2
!
!  Loop over atoms to search for images that may contribute to the screening
!
                        l = 0
                        npartialik = 0
                        do while (l.lt.numat.and.lnonzeroSik)
                          l = l + 1
                          neamspecl = neamfnspecptr(l)
                          if (lMEAMscreen(indik,neamspecl)) then
!
!  Set basic vectors between atoms
!
                            xil0 = xclat(l) - xal
                            yil0 = yclat(l) - yal
                            zil0 = zclat(l) - zal
                            xkl0 = xil0 - xcd1
                            ykl0 = yil0 - ycd1
                            zkl0 = zil0 - zcd1
!
!  Find images within cutoffs of both atoms - excluding self images
!
                            nvec0 = 0
                            call rfindmid(xil0,yil0,zil0,xkl0,ykl0,zkl0,rcut2,.false.,nvec0,nvec,vectorpair)
!
!  Loop over results of search
!
                            lvec = 0
                            do while (lvec.lt.nvec.and.lnonzeroSij)
                              lvec = lvec + 1
!
!  Compute screening function
!
                              call meamscreen(neamspeci,neamspeck,neamspecl,rik2,vectorpair%distance_pair1(lvec), &
                                              vectorpair%distance_pair2(lvec),Silk,dSilkdr,lpartial,.true.)
!
!  If screening function contribution is 0, then no need to continue for this pair
!
                              if (Silk.eq.0.0_dp) then
                                lnonzeroSik = .false.
                                Sik = 0.0_dp
                              else
!
!  Multiply total screening product
!
                                Sik = Sik*Silk
                                if (lpartial) then
!
!  If this atom has a screening factor between 0 and 1, we need to keep track of it since it will generate non-zero derivatives
!
                                  npartialik = npartialik + 1
                                  if (npartialik.gt.partialik%sa_maxdim) then
                                    call changemaxsa(partialik,npartialik)
                                  endif
                                  partialik%sa_atom(npartialik) = l
                                  partialik%sa_xik(npartialik) = vectorpair%xvector_pair1(lvec)
                                  partialik%sa_yik(npartialik) = vectorpair%yvector_pair1(lvec)
                                  partialik%sa_zik(npartialik) = vectorpair%zvector_pair1(lvec)
                                  partialik%sa_xjk(npartialik) = vectorpair%xvector_pair2(lvec)
                                  partialik%sa_yjk(npartialik) = vectorpair%yvector_pair2(lvec)
                                  partialik%sa_zjk(npartialik) = vectorpair%zvector_pair2(lvec)
                                  partialik%sa_Sikj(npartialik) = Silk
                                  partialik%sa_dSikjdr(1:3,npartialik) = dSilkdr(1:3)
                                endif
                              endif
!
!  End loop over images of possible screening atoms
!
                            enddo
                          endif
!
!  End loop over possible screening atoms
!
                        enddo
!
!  End of screening function
!
                      endif
!
!  Only do remainder of work if the screening factor is non-zero
!
                      if (lnonzeroSik) then
!
!  Scale density and derivatives by screening factor
!
                        drhoik(1:3,1:maxmeamcomponent) = Sik*drhoik(1:3,1:maxmeamcomponent)
                        drhoki(1:3,1:maxmeamcomponent) = Sik*drhoki(1:3,1:maxmeamcomponent)
!
                        call meamtotalrhoderv(neamspeci,scrho(1,i),rhoi,drhoik,drhototik,drhoiks,drhototiks, &
                                              drhoik2,drhototik2,drhoik2s,drhototik2s,drhoik2m,drhototik2m, &
                                              .false.,.false.)
                        call meamtotalrhoderv(neamspeck,scrho(1,k),rhok,drhoki,drhototki,drhokis,drhototkis, &
                                              drhoki2,drhototki2,drhoki2s,drhototki2s,drhoki2m,drhototki2m, &
                                              .false.,.false.)
                        call meamtotalrhocrossderv(neamspeci,scrho(1,i),rhoi,drhoij,drhototij,drhoik,drhototik, &
                                                   drhoijs,drhoiks,drhototijk2,drhototijk2s,drhototijk2m,.false.)
                      endif
                    else
                      call eamrho(nati,ntypi,natk,ntypk,rik,rpot(mp),xcd1,ycd1,zcd1,rhoik,rhoki,drhototik,drhototki, &
                                  drhototiks,drhototkis,drhototik2,drhototki2,drhototik2s,drhototki2s, &
                                  drhototik2m,drhototki2m,drhototik3,drhototki3,1.0_dp,1.0_dp,.false.,.true.,.false.,.false., &
                                  twopot(1,mp))
                    endif
                    lanyvalidik = .true.
                  endif
                endif
              enddo
            endif
            if (rjk2.gt.smallself.and.rjk2.le.cut2.and.npotjk.gt.0) then
!*********************
!  j-k contribution  *
!*********************
!
!  Zero terms
!
              if (lMEAM) then
                rhojk(1:maxmeamcomponent) = 0.0_dp
                rhokj(1:maxmeamcomponent) = 0.0_dp
                drhojk(1:3,1:maxmeamcomponent) = 0.0_dp
                drhokj(1:3,1:maxmeamcomponent) = 0.0_dp
              else
                rhojk(1) = 0.0_dp
                rhokj(1) = 0.0_dp
              endif
              drhototjk(1:3) = 0.0_dp
              drhototkj(1:3) = 0.0_dp
              drhototjik2(1:3,1:3) = 0.0_dp
!
              rjk = sqrt(rjk2)
!
!  Loop over potentials to find many-body ones
!
              do n = 1,npotjk
                mp = npotjkptr(n)
                lvalidjk = .false.
                if (rjk.gt.rpot2(mp).and.rjk.le.rpot(mp)) then
                  lvalidjk = .true.
                  if (lvalidjk) then
!**********************************
!  Calculate density derivatives  *
!**********************************
                    if (lMEAMden) then
                      call meamrho(natj,ntypj,natk,ntypk,rjk,rpot(mp),xcd2,ycd2,zcd2,rhojk,rhokj,drhojk,drhokj, &
                                   drhojks,drhokjs,drhojk2,drhokj2,drhojk2s,drhokj2s,drhojk2m,drhokj2m, &
                                   1.0_dp,1.0_dp,.true.,.false.,.true.,.false.,twopot(1,mp))
!*******************************
!  Compute screening function  *
!*******************************
!
!  Initialise screening function to 1
!
                      lnonzeroSjk = .true.
                      Sjk = 1.0_dp
!
                      if (lanyMEAMscreen) then
!
!  Set cutoffs for possible screening atoms
!
                        rcut2 = rcutfactor*rjk2
!
!  Loop over atoms to search for images that may contribute to the screening
!
                        l = 0
                        npartialjk = 0
                        do while (l.lt.numat.and.lnonzeroSjk)
                          l = l + 1
                          neamspecl = neamfnspecptr(l)
                          if (lMEAMscreen(indjk,neamspecl)) then
!
!  Set basic vectors between atoms
!
                            xjl0 = xclat(l) - xcd + xal
                            yjl0 = yclat(l) - ycd + yal
                            zjl0 = zclat(l) - zcd + zal
                            xkl0 = xjl0 - xcd2
                            ykl0 = yjl0 - ycd2
                            zkl0 = zjl0 - zcd2
!
!  Find images within cutoffs of both atoms - excluding self images
!
                            nvec0 = 0
                            call rfindmid(xjl0,yjl0,zjl0,xkl0,ykl0,zkl0,rcut2,.false.,nvec0,nvec,vectorpair)
!
!  Loop over results of search
!
                            lvec = 0
                            do while (lvec.lt.nvec.and.lnonzeroSjk)
                              lvec = lvec + 1
!
!  Compute screening function
!
                              call meamscreen(neamspecj,neamspeck,neamspecl,rjk2,vectorpair%distance_pair1(lvec), &
                                              vectorpair%distance_pair2(lvec),Sjlk,dSjlkdr,lpartial,.true.)
!
!  If screening function contribution is 0, then no need to continue for this pair
!
                              if (Sjlk.eq.0.0_dp) then
                                lnonzeroSjk = .false.
                                Sjk = 0.0_dp
                              else
!
!  Multiply total screening product
!
                                Sjk = Sjk*Sjlk
                                if (lpartial) then
!
!  If this atom has a screening factor between 0 and 1, we need to keep track of it since it will generate non-zero derivatives
!
                                  npartialjk = npartialjk + 1
                                  if (npartialjk.gt.partialjk%sa_maxdim) then
                                    call changemaxsa(partialjk,npartialjk)
                                  endif
                                  partialjk%sa_atom(npartialjk) = l
                                  partialjk%sa_xik(npartialjk) = vectorpair%xvector_pair1(lvec)
                                  partialjk%sa_yik(npartialjk) = vectorpair%yvector_pair1(lvec)
                                  partialjk%sa_zik(npartialjk) = vectorpair%zvector_pair1(lvec)
                                  partialjk%sa_xjk(npartialjk) = vectorpair%xvector_pair2(lvec)
                                  partialjk%sa_yjk(npartialjk) = vectorpair%yvector_pair2(lvec)
                                  partialjk%sa_zjk(npartialjk) = vectorpair%zvector_pair2(lvec)
                                  partialjk%sa_Sikj(npartialjk) = Sjlk
                                  partialjk%sa_dSikjdr(1:3,npartialjk) = dSjlkdr(1:3)
                                endif
                              endif
!
!  End loop over images of possible screening atoms
!
                            enddo
                          endif
!
!  End loop over possible screening atoms
!
                        enddo
!
!  End of screening function
!
                      endif
!
!  Only do remainder of work if the screening factor is non-zero
!
                      if (lnonzeroSjk) then
!
!  Scale density and derivatives by screening factor
!
                        drhojk(1:3,1:maxmeamcomponent) = Sjk*drhojk(1:3,1:maxmeamcomponent)
                        drhokj(1:3,1:maxmeamcomponent) = Sjk*drhokj(1:3,1:maxmeamcomponent)
                        call meamtotalrhoderv(neamspecj,scrho(1,j),rhoj,drhojk,drhototjk,drhojks,drhototjks, &
                                              drhojk2,drhototjk2,drhojk2s,drhototjk2s,drhojk2m,drhototjk2m, &
                                              .false.,.false.)
                        call meamtotalrhoderv(neamspeck,scrho(1,k),rhok,drhokj,drhototkj,drhokjs,drhototkjs, &
                                              drhokj2,drhototkj2,drhokj2s,drhototkj2s,drhokj2m,drhototkj2m, &
                                              .false.,.false.)
                        call meamtotalrhocrossderv(neamspecj,scrho(1,j),rhoj,drhoji,drhototji,drhojk,drhototjk, &
                                                   drhojis,drhojks,drhototjik2,drhototjik2s,drhototjik2m,.false.)
                      endif
                    else
                      call eamrho(natj,ntypj,natk,ntypk,rjk,rpot(mp),xcd2,ycd2,zcd2,rhojk,rhokj,drhototjk,drhototkj, &
                                  drhototjks,drhototkjs,drhototjk2,drhototkj2,drhototjk2s,drhototkj2s, &
                                  drhototjk2m,drhototkj2m,drhototjk3,drhototkj3,1.0_dp,1.0_dp,.false.,.true.,.false.,.false., &
                                  twopot(1,mp))
                    endif
                    lanyvalidjk = .true.
                  endif
                endif
              enddo
            endif
!
!  Cross term derivative for k-i / k-j
!
            if (lanyvalidik.and.lanyvalidjk.and.lMEAMden) then
              drhototkij2(1:3,1:3) = 0.0_dp
              call meamtotalrhocrossderv(neamspeck,scrho(1,k),rhok,drhoki,drhototki,drhokj,drhototkj,drhokis,drhokjs, &
                                         drhototkij2,drhototkij2s,drhototkij2m,.false.)
            endif
!**********************
!  Add terms to sums  *
!**********************
            if (lanyvalidik) then
              drhototiksum(1:3) = drhototiksum(1:3) + drhototik(1:3)
              drhototijk2sum(1:3,1:3) = drhototijk2sum(1:3,1:3) + drhototijk2(1:3,1:3)
            endif
            if (lanyvalidjk) then
              drhototjksum(1:3) = drhototjksum(1:3) + drhototjk(1:3)
              drhototjik2sum(1:3,1:3) = drhototjik2sum(1:3,1:3) + drhototjik2(1:3,1:3)
            endif
!****************************************************************************************************************
!  Calculate second derivatives for i-k/j-k terms - terms that have to be summed for each distance combination  *
!****************************************************************************************************************
!
!  i-k/j-k
!
            if (lanyvalidik.and.lanyvalidjk) then
              if (lMEAM) then
                dt1 = rscrhok3*ofctijk
                dt2 = rscrhok*ofctijk
                derv2(jx,ix) = derv2(jx,ix) + dt1*drhototki(1)*drhototkj(1) + dt2*drhototkij2(1,1)
                derv2(jy,ix) = derv2(jy,ix) + dt1*drhototki(1)*drhototkj(2) + dt2*drhototkij2(1,2)
                derv2(jz,ix) = derv2(jz,ix) + dt1*drhototki(1)*drhototkj(3) + dt2*drhototkij2(1,3)
                derv2(jx,iy) = derv2(jx,iy) + dt1*drhototki(2)*drhototkj(1) + dt2*drhototkij2(2,1)
                derv2(jy,iy) = derv2(jy,iy) + dt1*drhototki(2)*drhototkj(2) + dt2*drhototkij2(2,2)
                derv2(jz,iy) = derv2(jz,iy) + dt1*drhototki(2)*drhototkj(3) + dt2*drhototkij2(2,3)
                derv2(jx,iz) = derv2(jx,iz) + dt1*drhototki(3)*drhototkj(1) + dt2*drhototkij2(3,1)
                derv2(jy,iz) = derv2(jy,iz) + dt1*drhototki(3)*drhototkj(2) + dt2*drhototkij2(3,2)
                derv2(jz,iz) = derv2(jz,iz) + dt1*drhototki(3)*drhototkj(3) + dt2*drhototkij2(3,3)
              else
                dt1 = rscrhok3*ofctijk
                derv2(jx,ix) = derv2(jx,ix) + dt1*drhototki(1)*drhototkj(1)
                derv2(jy,ix) = derv2(jy,ix) + dt1*drhototki(1)*drhototkj(2)
                derv2(jz,ix) = derv2(jz,ix) + dt1*drhototki(1)*drhototkj(3)
                derv2(jx,iy) = derv2(jx,iy) + dt1*drhototki(2)*drhototkj(1)
                derv2(jy,iy) = derv2(jy,iy) + dt1*drhototki(2)*drhototkj(2)
                derv2(jz,iy) = derv2(jz,iy) + dt1*drhototki(2)*drhototkj(3)
                derv2(jx,iz) = derv2(jx,iz) + dt1*drhototki(3)*drhototkj(1)
                derv2(jy,iz) = derv2(jy,iz) + dt1*drhototki(3)*drhototkj(2)
                derv2(jz,iz) = derv2(jz,iz) + dt1*drhototki(3)*drhototkj(3)
              endif
            endif
!*****************************************************************************************
!  Calculate second derivatives for i-k/j-k terms that can be summed over all distances  *
!*****************************************************************************************
!
!  i-k
!
            if (lMEAM) then
              dt1 = rscrhoi3*ofctijk
              dt2 = rscrhoi*ofctijk
              derv2(jx,ix) = derv2(jx,ix) - dt1*drhototiksum(1)*drhototij(1) - dt2*drhototijk2sum(1,1)
              derv2(jy,ix) = derv2(jy,ix) - dt1*drhototiksum(1)*drhototij(2) - dt2*drhototijk2sum(2,1)
              derv2(jz,ix) = derv2(jz,ix) - dt1*drhototiksum(1)*drhototij(3) - dt2*drhototijk2sum(3,1)
              derv2(jx,iy) = derv2(jx,iy) - dt1*drhototiksum(2)*drhototij(1) - dt2*drhototijk2sum(1,2)
              derv2(jy,iy) = derv2(jy,iy) - dt1*drhototiksum(2)*drhototij(2) - dt2*drhototijk2sum(2,2)
              derv2(jz,iy) = derv2(jz,iy) - dt1*drhototiksum(2)*drhototij(3) - dt2*drhototijk2sum(3,2)
              derv2(jx,iz) = derv2(jx,iz) - dt1*drhototiksum(3)*drhototij(1) - dt2*drhototijk2sum(1,3)
              derv2(jy,iz) = derv2(jy,iz) - dt1*drhototiksum(3)*drhototij(2) - dt2*drhototijk2sum(2,3)
              derv2(jz,iz) = derv2(jz,iz) - dt1*drhototiksum(3)*drhototij(3) - dt2*drhototijk2sum(3,3)
            else
              dt1 = rscrhoi3*ofctijk
              derv2(jx,ix) = derv2(jx,ix) - dt1*drhototiksum(1)*drhototij(1)
              derv2(jy,ix) = derv2(jy,ix) - dt1*drhototiksum(1)*drhototij(2)
              derv2(jz,ix) = derv2(jz,ix) - dt1*drhototiksum(1)*drhototij(3)
              derv2(jx,iy) = derv2(jx,iy) - dt1*drhototiksum(2)*drhototij(1)
              derv2(jy,iy) = derv2(jy,iy) - dt1*drhototiksum(2)*drhototij(2)
              derv2(jz,iy) = derv2(jz,iy) - dt1*drhototiksum(2)*drhototij(3)
              derv2(jx,iz) = derv2(jx,iz) - dt1*drhototiksum(3)*drhototij(1)
              derv2(jy,iz) = derv2(jy,iz) - dt1*drhototiksum(3)*drhototij(2)
              derv2(jz,iz) = derv2(jz,iz) - dt1*drhototiksum(3)*drhototij(3)
            endif
!
!  j-k
!
            if (lMEAM) then
              dt1 = rscrhoj3*ofctijk
              dt2 = rscrhoj*ofctijk
              derv2(jx,ix) = derv2(jx,ix) + dt1*drhototjksum(1)*drhototji(1) + dt2*drhototjik2sum(1,1)
              derv2(jy,ix) = derv2(jy,ix) + dt1*drhototjksum(2)*drhototji(1) + dt2*drhototjik2sum(1,2)
              derv2(jz,ix) = derv2(jz,ix) + dt1*drhototjksum(3)*drhototji(1) + dt2*drhototjik2sum(1,3)
              derv2(jx,iy) = derv2(jx,iy) + dt1*drhototjksum(1)*drhototji(2) + dt2*drhototjik2sum(2,1)
              derv2(jy,iy) = derv2(jy,iy) + dt1*drhototjksum(2)*drhototji(2) + dt2*drhototjik2sum(2,2)
              derv2(jz,iy) = derv2(jz,iy) + dt1*drhototjksum(3)*drhototji(2) + dt2*drhototjik2sum(2,3)
              derv2(jx,iz) = derv2(jx,iz) + dt1*drhototjksum(1)*drhototji(3) + dt2*drhototjik2sum(3,1)
              derv2(jy,iz) = derv2(jy,iz) + dt1*drhototjksum(2)*drhototji(3) + dt2*drhototjik2sum(3,2)
              derv2(jz,iz) = derv2(jz,iz) + dt1*drhototjksum(3)*drhototji(3) + dt2*drhototjik2sum(3,3)
            else
              dt1 = rscrhoj3*ofctijk
              derv2(jx,ix) = derv2(jx,ix) + dt1*drhototjksum(1)*drhototji(1)
              derv2(jy,ix) = derv2(jy,ix) + dt1*drhototjksum(2)*drhototji(1)
              derv2(jz,ix) = derv2(jz,ix) + dt1*drhototjksum(3)*drhototji(1)
              derv2(jx,iy) = derv2(jx,iy) + dt1*drhototjksum(1)*drhototji(2)
              derv2(jy,iy) = derv2(jy,iy) + dt1*drhototjksum(2)*drhototji(2)
              derv2(jz,iy) = derv2(jz,iy) + dt1*drhototjksum(3)*drhototji(2)
              derv2(jx,iz) = derv2(jx,iz) + dt1*drhototjksum(1)*drhototji(3)
              derv2(jy,iz) = derv2(jy,iz) + dt1*drhototjksum(2)*drhototji(3)
              derv2(jz,iz) = derv2(jz,iz) + dt1*drhototjksum(3)*drhototji(3)
            endif
!***************************
!  End of third atom loop  *
!***************************
          enddo kloop
        endif
!**************************************
!  End of valid distance i-j section  *
!**************************************
      endif
    enddo jloop
  enddo iloop
!
!  End of real space part - perform general tasks
!
999 continue
!
!  Free local memory
!
  if (lgrad2) then
    deallocate(npotjkptr,stat=status)
    if (status/=0) call deallocate_error('many','npotjkptr')
    deallocate(npotikptr,stat=status)
    if (status/=0) call deallocate_error('many','npotikptr')
  endif
  deallocate(npotl,stat=status)
  if (status/=0) call deallocate_error('many','npotl')
!
!  Exit point
!
1000 continue
!
!  Unscale density
!
  call eamscalescrho(-1_i4)
!
!  Timing
!
  time2 = g_cpu_time()
  tmany = tmany + time2 - time1
#ifdef TRACE
  call trace_out('many0d')
#endif
!
  return
  end
