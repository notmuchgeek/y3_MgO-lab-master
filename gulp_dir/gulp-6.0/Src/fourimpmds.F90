  subroutine fourimpmds(eimp,esregion12,esregion2,eattach,lgrad1,xderv,yderv,zderv,rstrd)
!
!  Subroutine for four-body improper torsion energy and forces using spatial decomposition
!
!   7/13 Created from foundmds
!  12/14 Corrected so that ll vector is added to 42 and not 43
!   1/18 Cell index shifts handled to allow for nomod
!   2/18 Trace added
!   9/18 Strain module introduced
!   9/18 Call to fourstrterms replaced with more general realstrterms
!   4/19 Call to lmatchpair modified to include wildcard argument
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!   2/20 Correction to cell index handling for bonding
!   4/20 Rigid molecule modifications added
!   5/20 Centre of mass handling added
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
  use control,        only : lseok, latomicstress, lrigid
  use current
  use derivatives,    only : atomicstress
  use energies,       only : siteenergy, eregion2region
  use four
  use iochannels,     only : ioout
  use m_strain,       only : realstrterms
  use mdlogic
  use molecule
  use optimisation
  use parallel
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
  real(dp), intent(inout)                      :: eimp
  real(dp), intent(inout)                      :: esregion12
  real(dp), intent(inout)                      :: esregion2
  real(dp), intent(inout)                      :: eattach
  real(dp), intent(inout)                      :: rstrd(*)
  real(dp), intent(inout)                      :: xderv(*)
  real(dp), intent(inout)                      :: yderv(*)
  real(dp), intent(inout)                      :: zderv(*)
  logical,  intent(in)                         :: lgrad1
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ic
  integer(i4)                                  :: icx
  integer(i4)                                  :: icy
  integer(i4)                                  :: icz
  integer(i4)                                  :: ii
  integer(i4)                                  :: imx
  integer(i4)                                  :: imy
  integer(i4)                                  :: imz
  integer(i4)                                  :: isgn
  integer(i4)                                  :: ind
  integer(i4)                                  :: ind2
  integer(i4)                                  :: indmi
  integer(i4)                                  :: indmj
  integer(i4)                                  :: indmk
  integer(i4)                                  :: indml
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
  integer(i4)                                  :: ixk
  integer(i4)                                  :: iyk
  integer(i4)                                  :: izk
  integer(i4)                                  :: ixl
  integer(i4)                                  :: iyl
  integer(i4)                                  :: izl
  integer(i4)                                  :: ixyz
  integer(i4)                                  :: j
  integer(i4)                                  :: jc
  integer(i4)                                  :: jcx
  integer(i4)                                  :: jcy
  integer(i4)                                  :: jcz
  integer(i4)                                  :: jicx
  integer(i4)                                  :: jicy
  integer(i4)                                  :: jicz
  integer(i4)                                  :: jj
  integer(i4)                                  :: jkcx
  integer(i4)                                  :: jkcy
  integer(i4)                                  :: jkcz
  integer(i4)                                  :: jmx
  integer(i4)                                  :: jmy
  integer(i4)                                  :: jmz
  integer(i4)                                  :: jndn
  integer(i4)                                  :: k
  integer(i4)                                  :: kc
  integer(i4)                                  :: kcx
  integer(i4)                                  :: kcy
  integer(i4)                                  :: kcz
  integer(i4)                                  :: kk
  integer(i4)                                  :: kl
  integer(i4)                                  :: klcx
  integer(i4)                                  :: klcy
  integer(i4)                                  :: klcz
  integer(i4)                                  :: kmx
  integer(i4)                                  :: kmy
  integer(i4)                                  :: kmz
  integer(i4)                                  :: kndn
  integer(i4)                                  :: ks
  integer(i4)                                  :: l
  integer(i4)                                  :: lc
  integer(i4)                                  :: lcx
  integer(i4)                                  :: lcy
  integer(i4)                                  :: lcz
  integer(i4)                                  :: ll
  integer(i4)                                  :: maxx
  integer(i4)                                  :: maxxy
  integer(i4)                                  :: n
  integer(i4)                                  :: n1i
  integer(i4)                                  :: n1j
  integer(i4)                                  :: n1k
  integer(i4)                                  :: n1l
  integer(i4), dimension(:), allocatable       :: natmiddle
  integer(i4), dimension(:), allocatable       :: ntypmiddle
  integer(i4)                                  :: nbtypeji
  integer(i4)                                  :: nbtypejk
  integer(i4)                                  :: nbtypekl
  integer(i4)                                  :: nbtypeji2
  integer(i4)                                  :: nbtypejk2
  integer(i4)                                  :: nbtypekl2
  integer(i4)                                  :: nfortype
  integer(i4)                                  :: nforimp
  integer(i4)                                  :: ni
  integer(i4)                                  :: nil
  integer(i4)                                  :: niltor
  integer(i4)                                  :: nj
  integer(i4)                                  :: nk
  integer(i4)                                  :: nl
  integer(i4)                                  :: nmi 
  integer(i4)                                  :: nmiddle
  integer(i4)                                  :: nmj 
  integer(i4)                                  :: nmk 
  integer(i4)                                  :: nml 
  integer(i4)                                  :: nn
  integer(i4)                                  :: nnl
  integer(i4)                                  :: npha
  integer(i4), dimension(:), allocatable       :: nptrnforimp
  integer(i4)                                  :: nregioni
  integer(i4)                                  :: nregionj
  integer(i4)                                  :: nregionk
  integer(i4)                                  :: nregionl
  integer(i4)                                  :: nregiontypi
  integer(i4)                                  :: nregiontypj
  integer(i4)                                  :: nregiontypk
  integer(i4)                                  :: nregiontypl
  integer(i4)                                  :: nsplower(3)
  integer(i4)                                  :: nspupper(3)
  integer(i4)                                  :: nt1
  integer(i4)                                  :: nt2
  integer(i4)                                  :: nt3
  integer(i4)                                  :: nt4
  integer(i4)                                  :: ntyp1
  integer(i4)                                  :: ntyp2
  integer(i4)                                  :: ntyp3
  integer(i4)                                  :: ntyp4
  integer(i4)                                  :: ntypi
  integer(i4)                                  :: ntypj
  integer(i4)                                  :: ntypk
  integer(i4)                                  :: ntypl
  integer(i4)                                  :: status
  logical                                      :: l2bondsij
  logical                                      :: l2bondsjk
  logical                                      :: l2bondskl
  logical                                      :: lanybtyp
  logical                                      :: lanyneedmol
  logical                                      :: lattach
  logical                                      :: lbondedij
  logical                                      :: lbondedjk
  logical                                      :: lbondedkl
  logical                                      :: lbtyp
  logical                                      :: liok
  logical                                      :: limatch1
  logical                                      :: lintra_only
  logical                                      :: linter_only
  logical                                      :: ljmatch2
  logical                                      :: ljkmatch
  logical                                      :: lkmatch3
  logical                                      :: lmatch
  logical                                      :: lmatchany
  logical                                      :: lmatchpair
  logical                                      :: lmolloc
  logical                                      :: lmolok
  logical                                      :: lmolokjk
  logical                                      :: lneedmol
  logical                                      :: lopi
  logical                                      :: lopj
  logical                                      :: lopk
  logical                                      :: lopl
  logical                                      :: lreg12
  logical                                      :: lreg2qtet
  logical                                      :: lsamemolij
  logical                                      :: lsamemoljk
  logical                                      :: lsamemolkl
  logical                                      :: lsg1
  logical                                      :: lslicei
  logical                                      :: lslicej
  logical                                      :: lslicek
  logical                                      :: lslicel
  real(dp)                                     :: cut
  real(dp)                                     :: cutmax
  real(dp)                                     :: e1d(6)
  real(dp)                                     :: e2d(21)
  real(dp)                                     :: e3d(1)
  real(dp)                                     :: eterm
  real(dp)                                     :: eterm4th
  real(dp)                                     :: eterm6th
  real(dp)                                     :: fpoly(5)
  real(dp)                                     :: oci
  real(dp)                                     :: ocj
  real(dp)                                     :: ock
  real(dp)                                     :: ocl 
  real(dp)                                     :: ofct 
  real(dp)                                     :: phi0 
  real(dp)                                     :: phi0o
  real(dp)                                     :: r21
  real(dp)                                     :: r212
  real(dp)                                     :: r31
  real(dp)                                     :: r312
  real(dp)                                     :: r32
  real(dp)                                     :: r322
  real(dp)                                     :: r41
  real(dp)                                     :: r412
  real(dp)                                     :: r42
  real(dp)                                     :: r422
  real(dp)                                     :: r43
  real(dp)                                     :: r432
  real(dp)                                     :: rkforloc
  real(dp)                                     :: rko
  real(dp)                                     :: rstrdloc(6)
  real(dp)                                     :: rn
  real(dp)                                     :: dr2ds(6,6)
  real(dp)                                     :: d2r2dx2(3,3,6)
  real(dp)                                     :: d2r2ds2(6,6,6)
  real(dp)                                     :: d2r2dsdx(6,3,6)
  real(dp)                                     :: tr1
  real(dp)                                     :: tr2
  real(dp)                                     :: tr2max
  real(dp)                                     :: tr3
  real(dp)                                     :: tr4
  real(dp)                                     :: x21
  real(dp)                                     :: y21
  real(dp)                                     :: z21
  real(dp)                                     :: x31
  real(dp)                                     :: y31
  real(dp)                                     :: z31
  real(dp)                                     :: x41
  real(dp)                                     :: y41
  real(dp)                                     :: z41
  real(dp)                                     :: x32
  real(dp)                                     :: y32
  real(dp)                                     :: z32
  real(dp)                                     :: x42
  real(dp)                                     :: y42
  real(dp)                                     :: z42
  real(dp)                                     :: x43
  real(dp)                                     :: y43
  real(dp)                                     :: z43
  real(dp)                                     :: xc1
  real(dp)                                     :: yc1
  real(dp)                                     :: zc1
  real(dp)                                     :: xc2
  real(dp)                                     :: yc2
  real(dp)                                     :: zc2
  real(dp)                                     :: xc3
  real(dp)                                     :: yc3
  real(dp)                                     :: zc3
  real(dp)                                     :: xc4
  real(dp)                                     :: yc4
  real(dp)                                     :: zc4
  real(dp)                                     :: xcom(6)
  real(dp)                                     :: ycom(6)
  real(dp)                                     :: zcom(6)
  real(dp)                                     :: xcomj
  real(dp)                                     :: ycomj
  real(dp)                                     :: zcomj
  real(dp)                                     :: xv4(6)
  real(dp)                                     :: yv4(6)
  real(dp)                                     :: zv4(6)
#ifdef TRACE
  call trace_in('fourimpmds')
#endif
!
  lsg1 = (lstr.and.lgrad1)
  lmolloc = (nmol.gt.0)
!
!  Allocate local memory
!
  allocate(natmiddle(2*nfor),stat=status)
  if (status/=0) call outofmemory('fourimpmds','natmiddle')
  allocate(ntypmiddle(2*nfor),stat=status)
  if (status/=0) call outofmemory('fourimpmds','ntypmiddle')
  allocate(nptrnforimp(nfor),stat=status)
  if (status/=0) call outofmemory('fourimpmds','nptrnforimp')
!
!  Initialise centre of mass arrays to zero
!
  xcom(1:6) = 0.0_dp
  ycom(1:6) = 0.0_dp
  zcom(1:6) = 0.0_dp
!
!  Set variables for cell distribution
!
  maxxy = nspcell(1)*nspcell(2)
  maxx  = nspcell(1)
!
!  Check how many four-body potentials are of out of plane type and how many aren't
!
  nforimp = 0
  do n = 1,nfor
    if (mmfexc(n).eq.2) then
      nforimp = nforimp + 1
      nptrnforimp(nforimp) = n
    endif
  enddo
!
!  Find out if any require molecule information and whether any potential is of bonded type
!
  lanybtyp = .false.
  lanyneedmol = .false.
  do n = 1,nfor
    if (mmfexc(n).eq.2) then
      lbtyp = (mmfexc(n).ge.1)
      lintra_only = (lfintra(n).and..not.lfinter(n))
      linter_only = (lfinter(n).and..not.lfintra(n))
      lneedmol = (lintra_only.or.linter_only.or.lbtyp.or.lrigid)
      if (lneedmol) lanyneedmol = .true.
      if (lbtyp) lanybtyp = .true.
    endif
  enddo
!
!  Build a list of middle atom species types for potentials
!
  nmiddle = 0
  do n = 1,nfor
    if (mmfexc(n).eq.2) then
      if (.not.lmatchany(nfspec2(n),nfptyp2(n),nmiddle,natmiddle,ntypmiddle,.true.)) then
        nmiddle = nmiddle + 1
        natmiddle(nmiddle) = nfspec2(n)
        ntypmiddle(nmiddle) = nfptyp2(n)
      endif
      if (.not.lmatchany(nfspec3(n),nfptyp3(n),nmiddle,natmiddle,ntypmiddle,.true.)) then
        nmiddle = nmiddle + 1
        natmiddle(nmiddle) = nfspec3(n)
        ntypmiddle(nmiddle) = nfptyp3(n)
      endif
    endif
  enddo
!
!  Find maximum cutoff distance for ends and middle atoms
!
  cutmax = 0.0_dp
  tr2max = 0.0_dp
  do n = 1,nfor
    if (mmfexc(n).eq.2) then
      cut = max(for1(n)+for2(n),for2(n)+for3(n),for1(n)+for3(n),for4(n))
      cutmax = max(cut,cutmax)
      tr2 = for2(n)**2
      tr2max = max(tr2,tr2max)
    endif
  enddo
!
!  Opening banner for energy decomposition
!
  if (lPrintFour) then
    call mpbarrier
    if (ioproc) then
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''  IMP  : Atom No. 1  Atom No. 2  Atom No. 3  Atom No. 4    Torsion energy (eV)  '')')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
    endif
  endif
!************************************************************
!  Loop over all local spatial cells except buffer regions  *
!************************************************************
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
      nj = nspcellat(ind)
      n1j = nspcellat1ptr(ind)
!
!  Outer loop over atoms within this cell
!
      jloop: do jj = 1,nj
        j = nspcellatptr(n1j+jj)
!
!  Set properties of atom j
!
        nj = nat(j)
        ntypj = nftype(j)
!
!  Check whether species may be valid
!
        if (.not.lmatchany(nj,ntypj,nmiddle,natmiddle,ntypmiddle,.true.)) cycle jloop
!
        jc = nspcellatptrcell(n1j+jj)
        call ind2toijk(jc,jcx,jcy,jcz)
        nregionj = nregionno(nsft+nrelf2a(j))
        nregiontypj = nregiontype(nregionj,ncf)
        ocj = occuf(j)
        lopj = (lopf(nrelf2a(j)).or..not.lfreeze)
        lslicej = lsliceatom(nsft + nrelf2a(j))
!
!  Set coordinates of atom j
!
        xc2 = xinbox(j) + xvec2cell(jc)
        yc2 = yinbox(j) + yvec2cell(jc)
        zc2 = zinbox(j) + zvec2cell(jc)
!
!  Molecule handling
!
        if (lmolloc.and.lanyneedmol) then
          nmj = natmol(j)
          if (ndim.gt.0) then
            indmj = nmolind(j)
            call mindtoijk(indmj,ixj,iyj,izj)
          endif
!
!  Set COM coordinates
!
          if (lrigid.and.nmj.gt.0) then
            xcomj = molxyz(1,natinmol(j),nmj)
            ycomj = molxyz(2,natinmol(j),nmj)
            zcomj = molxyz(3,natinmol(j),nmj)
          else
            xcomj = 0.0_dp
            ycomj = 0.0_dp
            zcomj = 0.0_dp
          endif
        else
          xcomj = 0.0_dp
          ycomj = 0.0_dp
          zcomj = 0.0_dp
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
              nk = nspcellat(indn)
              n1k = nspcellat1ptr(indn)
              kloop: do kk = 1,nk
                k = nspcellatptr(n1k+kk)
!  
!  Prevent atoms j and k being the same atom and make sure that k < j
!               
                if (k.ge.j) cycle kloop
!  
!  Set properties for k
!
                nk = nat(k)
                ntypk = nftype(k)
!
!  Check whether species may be valid
!
                if (.not.lmatchany(nk,ntypk,nmiddle,natmiddle,ntypmiddle,.true.)) cycle kloop
                if (.not.lmatchpair(nj,ntypj,nk,ntypk,nfor,nfspec2,nfptyp2,nfspec3,nfptyp3,.true.)) cycle kloop
!
                kc = nspcellatptrcell(n1k+kk)
                call ind2toijk(kc,kcx,kcy,kcz)
!
                if (lmolloc.and.lanyneedmol) then
!
!  Molecule handling
!
                  nmk = natmol(k)
                  if (ndim.gt.0) then
                    indmk = nmolind(k)
                    call mindtoijk(indmk,ixk,iyk,izk)
!
                    ixk = ixk - ixj
                    iyk = iyk - iyj
                    izk = izk - izj
                  endif
                  lmolokjk = (nmj.eq.nmk.and.nmj.gt.0)
!
!  Set COM coordinates
!
                  if (lrigid.and.nmk.gt.0.and.lmolokjk) then
                    xcom(4) = molxyz(1,natinmol(k),nmk) - xcomj
                    ycom(4) = molxyz(2,natinmol(k),nmk) - ycomj
                    zcom(4) = molxyz(3,natinmol(k),nmk) - zcomj
                  else
                    xcom(4) = - xcomj
                    ycom(4) = - ycomj
                    zcom(4) = - zcomj
                  endif
                else
                  lmolokjk = .false.
                  xcom(4) = - xcomj
                  ycom(4) = - ycomj
                  zcom(4) = - zcomj
                endif
!
!  Calculate vector from atom 1 to atom 2
!
                xc3 = xvec2cell(kc) + xinbox(k)
                yc3 = yvec2cell(kc) + yinbox(k)
                zc3 = zvec2cell(kc) + zinbox(k)
                x32 = xc3 - xc2
                y32 = yc3 - yc2
                z32 = zc3 - zc2
!
!  Check r32 is OK
!
                r322 = x32*x32 + y32*y32 + z32*z32
                if (r322.lt.1.0d-12) cycle kloop
!
!  Molecule checking
!
                lbondedjk = .false.
                if (lmolokjk) then
                  if (ndim.eq.0) then
                    if (lanybtyp) then
                      call bonded(lbondedjk,l2bondsjk,nbtypejk,nbtypejk2,j,k,0_i4,0_i4,0_i4)
                    endif
                  else
                    jkcx = kcx - jcx + icosxsp(k) - icosxsp(j)
                    jkcy = kcy - jcy + icosysp(k) - icosysp(j)
                    jkcz = kcz - jcz + icoszsp(k) - icoszsp(j)
                    if (lanybtyp) then
                      call bonded(lbondedjk,l2bondsjk,nbtypejk,nbtypejk2,j,k,jkcx,jkcy,jkcz)
                      lsamemoljk = (lbondedjk.or.l2bondsjk)
                    else
                      lsamemoljk = .false.
                    endif
                    if (.not.lsamemoljk) then
                      call samemol(lsamemoljk,nmj,jkcx,jkcy,jkcz,ixk,iyk,izk)
                    endif
                  endif
                endif
!
!  Distance checking
!
                if (r322.gt.tr2max.and.(.not.lanybtyp.or..not.lbondedjk)) cycle kloop
!  
!  Set remaining properties for k
!               
                nregionk = nregionno(nsft+nrelf2a(k))
                nregiontypk = nregiontype(nregionk,ncf)
                ock = occuf(k)
                lopk = (lopf(nrelf2a(k)).or..not.lfreeze)
                lslicek = lsliceatom(nsft + nrelf2a(k))
!
!  Set counter for number of valid i/l end atom combinations
!
                niltor = 0
!***********************************
!  Loop over four-body potentials  *
!***********************************
                pots: do nn = 1,nforimp
                  n = nptrnforimp(nn)
                  nfortype = nforty(n)
                  nt1 = nfspec1(n)
                  nt2 = nfspec2(n)
                  nt3 = nfspec3(n)
                  nt4 = nfspec4(n)
                  ntyp1 = nfptyp1(n)
                  ntyp2 = nfptyp2(n)
                  ntyp3 = nfptyp3(n)
                  ntyp4 = nfptyp4(n)
                  tr1 = for1(n)**2
                  tr2 = for2(n)**2
                  tr3 = for3(n)**2
                  tr4 = for4(n)**2
                  lbtyp = (mmfexc(n).ge.1)
                  lintra_only = (lfintra(n).and..not.lfinter(n))
                  linter_only = (lfinter(n).and..not.lfintra(n))
                  lneedmol = (lintra_only.or.linter_only.or.lbtyp.or.lrigid)
!
!  Reset lmolok to initial state for j-k pair for each potential
!
                  lmolok = lmolokjk
!
!  QM/MM handling : j & k are both QM atoms and potential is of bonded type => exclude
!
                  if (QMMMmode(ncf).gt.0) then
                    if (nregiontypj.eq.1.and.nregiontypk.eq.1.and.lbtyp) cycle pots
                  endif
!************************************
!  Validate potential for j-k pair  *
!************************************
!
!  Check whether j and k are allowed for n
!
                  ljmatch2 = lmatch(nj,ntypj,nt2,ntyp2,.true.)
                  lkmatch3 = lmatch(nk,ntypk,nt3,ntyp3,.true.)
!
!  Check whether j-k are OK for 2-3
!
                  ljkmatch = (ljmatch2.and.lkmatch3)
!
!  If no pair of matches can be found then cycle
!
                  if (.not.ljkmatch) cycle pots
!
!  Distance checking for j-k
!
                  if (lbtyp) then
                    if (.not.lbondedjk) cycle pots
                  else
                    if (r322.gt.tr2) cycle pots
                  endif
!       
!  Check for intra and but not in same molecule
!       
                  if (lintra_only.and..not.lmolok) cycle pots
                  if (lbtyp.and..not.lmolok) cycle pots
!                 
!  Molecule checking
!
                  if (lmolok) then
                    if (ndim.eq.0) then
                      if (linter_only) cycle pots
                      if (lbtyp) then
                        if (.not.lbondedjk) cycle pots
!               
!  Check central bond type for correct order
!
                        if (n4botype(1,n).gt.0) then
                          if (n4botype(1,n).ne.nbtypejk) cycle pots
                        endif
                        if (n4botype(2,n).ne.nbtypejk2) cycle pots
                      endif
                    else
                      if (lbtyp) then
                        if (.not.lbondedjk) cycle pots
!               
!  Check central bond type for correct order
!
                        if (n4botype(1,n).gt.0) then
                          if (n4botype(1,n).ne.nbtypejk) cycle pots
                        endif
                        if (n4botype(2,n).ne.nbtypejk2) cycle pots
                      endif
                      if (lintra_only.and..not.lsamemoljk) cycle pots
                      if (linter_only.and.lsamemoljk) cycle pots
                    endif
                  endif
!  
!  Loop over neighbouring cells for i
!               
                  do jmz = nsplower(3),nspupper(3)
                    do jmy = nsplower(2),nspupper(2)
                      do jmx = nsplower(1),nspupper(1)
                        jndn = (jmz-1)*maxxy + (jmy-1)*maxx + jmx
!  
!  Loop over atoms within neighbouring cells
!                     
                        ni = nspcellat(jndn)
                        n1i = nspcellat1ptr(jndn)
                        iloop: do ii = 1,ni
                          i = nspcellatptr(n1i+ii)
!  
!  Prevent atoms i and k being the same atom
!                       
                          if (k.eq.i) cycle iloop
!  
!  Set properties of atom i
!                       
                          ic = nspcellatptrcell(n1i+ii)
                          call ind2toijk(ic,icx,icy,icz)
                          ni = nat(i)
                          ntypi = nftype(i)
!
!  Check whether i matches either of types 1 and 4
!
                          limatch1 = lmatch(ni,ntypi,nt1,ntyp1,.true.)
!
!  Is i allowed for type 1?
!
                          liok = limatch1
                          if (.not.liok) cycle iloop
!
!  Molecule handling
!
                          if (lmolloc.and.lneedmol) then
                            nmi = natmol(i)
                            if (ndim.gt.0) then
                              indmi = nmolind(i)
                              call mindtoijk(indmi,ixi,iyi,izi)
!
                              ixi = ixi - ixj
                              iyi = iyi - iyj
                              izi = izi - izj
                            endif
                            lmolok = (nmj.eq.nmi.and.nmj.gt.0)
!
!  Set COM coordinates
!
                            if (lrigid.and.nmi.gt.0.and.lmolok) then
                              xcom(1) = xcomj - molxyz(1,natinmol(i),nmi)
                              ycom(1) = ycomj - molxyz(2,natinmol(i),nmi)
                              zcom(1) = zcomj - molxyz(3,natinmol(i),nmi)
                            else
                              xcom(1) = xcomj
                              ycom(1) = ycomj
                              zcom(1) = zcomj
                            endif
                            xcom(2) = xcom(4) + xcom(1)
                            ycom(2) = ycom(4) + ycom(1)
                            zcom(2) = zcom(4) + zcom(1)
                          else
                            lmolok = .false.
                            xcom(1) = xcomj
                            ycom(1) = ycomj
                            zcom(1) = zcomj
                            xcom(2) = xcom(4) + xcom(1)
                            ycom(2) = ycom(4) + ycom(1)
                            zcom(2) = zcom(4) + zcom(1)
                          endif
!
!  Check for intra and but not in same molecule
!
                          if (lintra_only.and..not.lmolok) cycle iloop
                          if (lbtyp.and..not.lmolok) cycle iloop
!
!  Calculate vectors
!
                          xc1 = xvec2cell(ic) + xinbox(i)
                          yc1 = yvec2cell(ic) + yinbox(i)
                          zc1 = zvec2cell(ic) + zinbox(i)
!
!  Check r21 is OK
!
                          x21 = xc2 - xc1
                          y21 = yc2 - yc1
                          z21 = zc2 - zc1
                          r212 = x21*x21 + y21*y21 + z21*z21
                          if (r212.lt.1d-12) cycle iloop
!
!  Check r32 is OK
!
                          x31 = xc3 - xc1
                          y31 = yc3 - yc1
                          z31 = zc3 - zc1
                          r312 = x31*x31 + y31*y31 + z31*z31
                          if (r312.lt.1.0d-12) cycle iloop
!
!  Molecule checking
!
                          lbondedij = .false.
                          if (lmolok) then
                            if (ndim.eq.0) then
                              if (linter_only) cycle iloop
                              if (lbtyp) then
                                call bonded(lbondedij,l2bondsij,nbtypeji,nbtypeji2,j,i,0_i4,0_i4,0_i4)
                                if (.not.lbondedij) cycle iloop
                              endif
                            else
                              jicx = icx - jcx + icosxsp(i) - icosxsp(j)
                              jicy = icy - jcy + icosysp(i) - icosysp(j)
                              jicz = icz - jcz + icoszsp(i) - icoszsp(j)
                              if (lbtyp) then
                                call bonded(lbondedij,l2bondsij,nbtypeji,nbtypeji2,j,i,jicx,jicy,jicz)
                                if (.not.lbondedij) cycle iloop
                                lsamemolij = (lbondedij.or.l2bondsij)
                              else
                                lsamemolij = .false.
                              endif
                              if (.not.lsamemolij) then
                                call samemol(lsamemolij,nmj,jicx,jicy,jicz,ixi,iyi,izi)
                              endif
                              if (lintra_only.and..not.lsamemolij) cycle iloop
                              if (linter_only.and.lsamemolij) cycle iloop
                            endif
                          endif
!
!  Distance checking
!
                          if (lbtyp) then
                            if (.not.lbondedij) cycle iloop
                          else
                            if (r212.gt.tr1) cycle iloop
                          endif
! 
!  Set remaining properties of atom i
!                          
                          nregioni = nregionno(nsft+nrelf2a(i))
                          nregiontypi = nregiontype(nregioni,ncf)
                          oci = occuf(i)
                          lopi = (lopf(nrelf2a(i)).or..not.lfreeze)
                          lslicei = lsliceatom(nsft + nrelf2a(i))
!
!  Loop over neighbouring cells for l
!
                          do kmz = nsplower(3),nspupper(3)
                            do kmy = nsplower(2),nspupper(2)
                              do kmx = nsplower(1),nspupper(1)
                                kndn = (kmz-1)*maxxy + (kmy-1)*maxx + kmx
!
!  Loop over atoms within neighbouring cells
!                                
                                nnl = nspcellat(kndn)
                                n1l = nspcellat1ptr(kndn)
                                lloop: do ll = 1,nnl
                                  l = nspcellatptr(n1l+ll)
!  
!  Prevent l and j from being the same atom
!                                  
                                  if (l.eq.j) cycle lloop
!  
!  If lfreeze=.true. and no atoms have any variables
!  then skip this four body term
!                               
                                  lopl = (lopf(nrelf2a(l)).or..not.lfreeze)
                                  if (.not.lopi.and..not.lopj.and..not.lopk.and..not.lopl) cycle lloop
!
!  Check l is allowed for n
!             
                                  nl = nat(l)
                                  ntypl = nftype(l)
                                  if (.not.lmatch(nl,ntypl,nt4,ntyp4,.true.)) cycle lloop
!  
!  Set remaining properties of l
!                               
                                  nregionl = nregionno(nsft+nrelf2a(l))
                                  nregiontypl = nregiontype(nregionl,ncf)
                                  ocl = occuf(l)
                                  lc = nspcellatptrcell(n1l+ll)
                                  call ind2toijk(lc,lcx,lcy,lcz)
!
                                  if (lmolloc.and.lanyneedmol) then
!
!  Molecule handling
!
                                    nml = natmol(l)
                                    if (ndim.gt.0) then
                                      indml = nmolind(l)
                                      call mindtoijk(indml,ixl,iyl,izl)
!
                                      ixl = ixl - ixj
                                      iyl = iyl - iyj
                                      izl = izl - izj
                                    endif
                                    lmolok = (nmj.eq.nml.and.nmj.gt.0)
!
!  Set COM coordinates
!
                                    if (lrigid.and.nml.gt.0.and.lmolok) then
                                      xcom(5) = molxyz(1,natinmol(l),nml) - xcomj
                                      ycom(5) = molxyz(2,natinmol(l),nml) - ycomj
                                      zcom(5) = molxyz(3,natinmol(l),nml) - zcomj
                                    else
                                      xcom(5) = - xcomj
                                      ycom(5) = - ycomj
                                      zcom(5) = - zcomj
                                    endif
                                    xcom(3) = xcom(5) + xcom(1)
                                    ycom(3) = ycom(5) + ycom(1)
                                    zcom(3) = zcom(5) + zcom(1)
                                    xcom(6) = xcom(5) - xcom(4)
                                    ycom(6) = ycom(5) - ycom(4)
                                    zcom(6) = zcom(5) - zcom(4)
                                  else
                                    lmolok = .false.
                                    xcom(5) = - xcomj
                                    ycom(5) = - ycomj
                                    zcom(5) = - zcomj
                                    xcom(3) = xcom(5) + xcom(1)
                                    ycom(3) = ycom(5) + ycom(1)
                                    zcom(3) = zcom(5) + zcom(1)
                                    xcom(6) = xcom(5) - xcom(4)
                                    ycom(6) = ycom(5) - ycom(4)
                                    zcom(6) = zcom(5) - zcom(4)
                                  endif
!
!  Check for intra and but not in same molecule
!
                                  if (lintra_only.and..not.lmolok) cycle lloop
                                  if (lbtyp.and..not.lmolok) cycle lloop
!                                 
!  Calculate vectors            
!
                                  xc4 = xinbox(l) + xvec2cell(lc)
                                  yc4 = yinbox(l) + yvec2cell(lc)
                                  zc4 = zinbox(l) + zvec2cell(lc)
                                  x43 = xc4 - xc2
                                  y43 = yc4 - yc2
                                  z43 = zc4 - zc2
!                               
!  Check r42 is OK
!
                                  r422 = x42*x42 + y42*y42 + z42*z42
                                  if (r422.lt.1d-12) cycle lloop
!
!  Set region 2 quartet flag
!
                                  lreg12    = .false.
                                  lreg2qtet = .false.
                                  if (lseok.and.nregions(ncf).gt.1) then
                                    lreg2qtet = (nregioni.gt.1.and.nregionj.gt.1.and.nregionk.gt.1.and.nregionl.gt.1)
                                    if (.not.lreg2qtet) lreg12 = (nregioni.gt.1.or.nregionj.gt.1.or.nregionk.gt.1.or.nregionl.gt.1)
                                  endif
                                  lslicel = lsliceatom(nsft + nrelf2a(l))
                                  lattach = .true.
                                  if (lslicei.and.lslicej.and.lslicek.and.lslicel) lattach = .false.  
                                  if (.not.lslicei.and..not.lslicej.and..not.lslicek.and..not.lslicel) lattach = .false.
!
!  QM/MM handling : i, j, k & l are all QM atoms => exclude
!
                                  if (QMMMmode(ncf).gt.0) then
                                    if (nregiontypi.eq.1.and.nregiontypj.eq.1.and.nregiontypk.eq.1.and.nregiontypl.eq.1) cycle lloop
                                  endif
!
!  Molecule checking
!
                                  lbondedkl = .false.
                                  if (lmolok) then
                                    if (ndim.eq.0) then
                                      if (linter_only) cycle lloop
                                      if (lbtyp) then
                                        call bonded(lbondedkl,l2bondskl,nbtypekl,nbtypekl2,j,l,0_i4,0_i4,0_i4)
                                        if (.not.lbondedkl) cycle lloop
                                      endif
                                    else
                                      klcx = lcx - jcx + icosxsp(l) - icosxsp(j)
                                      klcy = lcy - jcy + icosysp(l) - icosysp(j)
                                      klcz = lcz - jcz + icoszsp(l) - icoszsp(j)
                                      if (lbtyp) then
                                        call bonded(lbondedkl,l2bondskl,nbtypekl,nbtypekl2,j,l,klcx,klcy,klcz)
                                        if (.not.lbondedkl) cycle lloop
                                        lsamemolkl = (lbondedkl.or.l2bondskl)
                                      else
                                        lsamemolkl = .false.
                                      endif
                                      if (.not.lsamemolkl) then
                                        call samemol(lsamemolkl,nmj,klcx,klcy,klcz,ixl,iyl,izl)
                                      endif
                                      if (lintra_only.and..not.lsamemolkl) cycle lloop
                                      if (linter_only.and.lsamemolkl) cycle lloop
                                    endif
                                  endif
!
!  Distance checking
!
                                  if (lbtyp) then
                                    if (.not.lbondedkl) cycle lloop
                                  else
                                    if (r422.gt.tr3) cycle lloop
                                  endif
!
!  Check r41 is OK
!
                                  x41 = x43 + x32 + x21
                                  y41 = y43 + y32 + y21
                                  z41 = z43 + z32 + z21
                                  r412 = x41*x41 + y41*y41 + z41*z41
                                  if (r412.lt.1.0d-12) cycle lloop
!
!  Check r43 is OK
!
                                  x43 = x42 - x32
                                  y43 = y42 - y32
                                  z43 = z42 - z32
                                  r432 = x43*x43 + y43*y43 + z43*z43
                                  if (r432.lt.1.0d-12) cycle lloop
!*********************************
!  Valid four-body term located  *
!*********************************
!
!  Finish calculating distances
!
                                  r21 = sqrt(r212)
                                  r31 = sqrt(r312)
                                  r32 = sqrt(r322)
                                  r41 = sqrt(r412)
                                  r42 = sqrt(r422)
                                  r43 = sqrt(r432)
!
!  Store information into iltor arrays
!
                                  niltor = niltor + 1
                                  if (niltor.gt.maxiltor) then
                                    maxiltor = niltor + 3
                                    call changemaxiltor
                                  endif
!
                                  nfortor(niltor) = n
!
                                  iltor(1,niltor) = i
                                  iltor(2,niltor) = l
!
                                  riltor(1,niltor) = r21
                                  riltor(2,niltor) = r31
                                  riltor(3,niltor) = r41
                                  riltor(4,niltor) = r42
                                  riltor(5,niltor) = r43
!
                                  xiltor(1,niltor) = x21
                                  yiltor(1,niltor) = y21
                                  ziltor(1,niltor) = z21
                                  xiltor(2,niltor) = x31
                                  yiltor(2,niltor) = y31
                                  ziltor(2,niltor) = z31
                                  xiltor(3,niltor) = x41
                                  yiltor(3,niltor) = y41
                                  ziltor(3,niltor) = z41
                                  xiltor(4,niltor) = x42
                                  yiltor(4,niltor) = y42
                                  ziltor(4,niltor) = z42
                                  xiltor(5,niltor) = x43
                                  yiltor(5,niltor) = y43
                                  ziltor(5,niltor) = z43
!
                                  if (lrigid) then
                                    xcomtor(1,niltor) = xcom(1)
                                    ycomtor(1,niltor) = ycom(1)
                                    zcomtor(1,niltor) = zcom(1)
                                    xcomtor(2,niltor) = xcom(2)
                                    ycomtor(2,niltor) = ycom(2)
                                    zcomtor(2,niltor) = zcom(2)
                                    xcomtor(3,niltor) = xcom(3)
                                    ycomtor(3,niltor) = ycom(3)
                                    zcomtor(3,niltor) = zcom(3)
                                    xcomtor(4,niltor) = xcom(5)
                                    ycomtor(4,niltor) = ycom(5)
                                    zcomtor(4,niltor) = zcom(5)
                                    xcomtor(5,niltor) = xcom(6)
                                    ycomtor(5,niltor) = ycom(6)
                                    zcomtor(5,niltor) = zcom(6)
                                  endif
!
                                  oiltor(niltor) = oci*ocj*ock*ocl
                                  lopiltor(1,niltor) = lopi
                                  lopiltor(2,niltor) = lopl
                                  lsurfiltor(1,niltor) = lreg12
                                  lsurfiltor(2,niltor) = lreg2qtet
                                  lsurfiltor(3,niltor) = lattach
!
!  End of loop over l             
!                                 
                                enddo lloop
!                                 
!  End of loops over cell images for l
!                                 
                              enddo 
                            enddo   
                          enddo
!                                 
!  End of loop over i             
!                                 
                        enddo iloop
!  
!  End of loops over cell images
!                   
                      enddo
                    enddo
                  enddo
!
!  End loop over potentials
!
                enddo pots
!*******************************
!  Loop over i/l combinations  *
!*******************************
                do nil = 1,niltor
!
!  Return values to local variables
!
                  n = nfortor(nil)
!
                  i = iltor(1,nil)
                  l = iltor(2,nil)
!
                  r21 = riltor(1,nil)
                  r31 = riltor(2,nil)
                  r41 = riltor(3,nil)
                  r42 = riltor(4,nil)
                  r43 = riltor(5,nil)
!
                  x21 = xiltor(1,nil)
                  y21 = yiltor(1,nil)
                  z21 = ziltor(1,nil)
                  x31 = xiltor(2,nil)
                  y31 = yiltor(2,nil)
                  z31 = ziltor(2,nil)
                  x41 = xiltor(3,nil)
                  y41 = yiltor(3,nil)
                  z41 = ziltor(3,nil)
                  x42 = xiltor(4,nil)
                  y42 = yiltor(4,nil)
                  z42 = ziltor(4,nil)
                  x43 = xiltor(5,nil)
                  y43 = yiltor(5,nil)
                  z43 = ziltor(5,nil)
!
                  if (lrigid) then
                    xcom(1) = xcomtor(1,nil)
                    ycom(1) = ycomtor(1,nil)
                    zcom(1) = zcomtor(1,nil)
                    xcom(2) = xcomtor(2,nil)
                    ycom(2) = ycomtor(2,nil)
                    zcom(2) = zcomtor(2,nil)
                    xcom(3) = xcomtor(3,nil)
                    ycom(3) = ycomtor(3,nil)
                    zcom(3) = zcomtor(3,nil)
                    xcom(5) = xcomtor(4,nil)
                    ycom(5) = ycomtor(4,nil)
                    zcom(5) = zcomtor(4,nil)
                    xcom(6) = xcomtor(5,nil)
                    ycom(6) = ycomtor(5,nil)
                    zcom(6) = zcomtor(5,nil)
                  endif
!
                  ofct = oiltor(nil)
                  lopi = lopiltor(1,nil)
                  lopl = lopiltor(2,nil)
                  lreg12 = lsurfiltor(1,nil) 
                  lreg2qtet = lsurfiltor(2,nil) 
                  lattach = lsurfiltor(3,nil) 
!
!  Set terms for potentials
!
                  rkforloc = fork(n)
!
!  If this is Dreiding mode then divide force constant by number of torsions
!
                  if (lfdreiding(n)) then
                    rkforloc = rkforloc/dble(niltor)
                  endif
                  npha = 0
                  nfortype = nforty(n)
                  if (nfortype.eq.1) then
                    npha = npfor(n)
                    if (npha.gt.0) then
                      isgn = 1
                    else
                      isgn = - 1
                    endif
                    npha = abs(npha)
                    phi0 = forpoly(1,n)*degtorad
                  elseif (nfortype.eq.4.or.nfortype.eq.6.or.nfortype.eq.7) then
                    npha = npfor(n)
                    if (npha.gt.0) then
                      isgn = 1
                    else
                      isgn = - 1
                    endif
                    npha = abs(npha)
                    if (nfortype.eq.6) then
                      phi0 = forpoly(1,n)*degtorad
                    else
                      phi0 = forpoly(1,n)
                    endif
                    if (nfortype.eq.6.or.nfortype.eq.7) then
                      fpoly(2:4) = forpoly(2:4,n)
                    endif
                  elseif (nfortype.eq.8.or.nfortype.eq.9) then
                    npha = npfor(n)
                    if (npha.gt.0) then
                      isgn = 1
                    else
                      isgn = - 1
                    endif
                    npha = abs(npha)
                    if (nfortype.eq.8) then
                      phi0 = forpoly(1,n)*degtorad
                    else
                      phi0 = forpoly(1,n)
                    endif
                    fpoly(2) = forpoly(2,n)
                    fpoly(3) = for1(n)
                    fpoly(4) = for2(n)
                    fpoly(5) = for3(n)
                  elseif (nfortype.eq.2) then
                    npha = npfor(n)
                  elseif (nfortype.eq.5) then
                    phi0 = forpoly(1,n)*degtorad
                  elseif (nfortype.eq.10.or.nfortype.eq.17) then
                    fpoly(1) = forpoly(1,n)*degtorad
                    fpoly(2) = forpoly(2,n)*degtorad
                  elseif (nfortype.eq.13) then
                    npha = abs(npfor(n))
                    phi0 = forpoly(1,n)*degtorad
                  endif
                  rn = dble(npha)
!
!  Scaling of terms
!
                  rko = rkforloc*ofct
                  phi0o = phi0
                  if (nfortype.eq.2) then
                    do kl = 1,npha
                      fpoly(kl) = forpoly(kl,n)*ofct
                    enddo
                  elseif (nfortype.eq.4.or.nfortype.eq.7.or.nfortype.eq.9) then
                    phi0o = phi0*ofct
                  endif
!
!  Call subroutine to calculate energy and derivatives
!
                  call fourbody(n,nfortype,r21,r31,r41,r32,r42,r43,eterm,e1d,e2d,e3d, &
                                rko,rn,phi0o,isgn,fpoly,lgrad1,.false.,.false.)
                  if (lreg2qtet) then
                    esregion2 = esregion2 + eterm
                  elseif (lreg12) then
                    esregion12 = esregion12 + eterm
                  else
                    eimp = eimp + eterm
                  endif
                  if (lattach) eattach = eattach + eterm
!
                  eterm6th = eterm/6.0_dp
                  eregion2region(nregionj,nregioni) = eregion2region(nregionj,nregioni) + eterm6th
                  eregion2region(nregionk,nregioni) = eregion2region(nregionk,nregioni) + eterm6th
                  eregion2region(nregionl,nregioni) = eregion2region(nregionl,nregioni) + eterm6th
                  eregion2region(nregionk,nregionj) = eregion2region(nregionk,nregionj) + eterm6th
                  eregion2region(nregionl,nregionj) = eregion2region(nregionl,nregionj) + eterm6th
                  eregion2region(nregionl,nregionk) = eregion2region(nregionl,nregionk) + eterm6th
!
                  eterm4th = 0.25_dp*eterm
                  siteenergy(i) = siteenergy(i) + eterm4th
                  siteenergy(j) = siteenergy(j) + eterm4th
                  siteenergy(k) = siteenergy(k) + eterm4th
                  siteenergy(l) = siteenergy(l) + eterm4th
!
!  Output energy contribution
!
                  if (lPrintFour) then
                    write(ioout,'(4x,4i12,1x,f22.10)') i,j,k,l,eterm
                  endif
!**************************
!  Torsional derivatives  *
!**************************
                  if (lgrad1) then
!
!  Set up strain products
!
                    if (lsg1) then
                      xv4(1) = x21
                      yv4(1) = y21
                      zv4(1) = z21
                      xv4(2) = x31
                      yv4(2) = y31
                      zv4(2) = z31
                      xv4(3) = x41
                      yv4(3) = y41
                      zv4(3) = z41
                      xv4(4) = x32
                      yv4(4) = y32
                      zv4(4) = z32
                      xv4(5) = x42
                      yv4(5) = y42
                      zv4(5) = z42
                      xv4(6) = x43
                      yv4(6) = y43
                      zv4(6) = z43
                      call realstrterms(ndim,6_i4,6_i4,xv4,yv4,zv4,xcom,ycom,zcom,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
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
                      rstrdloc(kl) = rstrdloc(kl) + e1d(1)*dr2ds(ks,1)
                      rstrdloc(kl) = rstrdloc(kl) + e1d(2)*dr2ds(ks,2)
                      rstrdloc(kl) = rstrdloc(kl) + e1d(3)*dr2ds(ks,3)
                      rstrdloc(kl) = rstrdloc(kl) + e1d(4)*dr2ds(ks,4)
                      rstrdloc(kl) = rstrdloc(kl) + e1d(5)*dr2ds(ks,5)
                      rstrdloc(kl) = rstrdloc(kl) + e1d(6)*dr2ds(ks,6)
                      rstrd(kl) = rstrd(kl) + rstrdloc(kl)
                    enddo
                    if (latomicstress) then
                      do kl = 1,nstrains
                        atomicstress(kl,i) = atomicstress(kl,i) + 0.25_dp*rstrdloc(kl)
                        atomicstress(kl,j) = atomicstress(kl,j) + 0.25_dp*rstrdloc(kl)
                        atomicstress(kl,k) = atomicstress(kl,k) + 0.25_dp*rstrdloc(kl)
                        atomicstress(kl,l) = atomicstress(kl,l) + 0.25_dp*rstrdloc(kl)
                      enddo
                    endif
                  endif
!*************************
!  Internal derivatives  *
!*************************
                  if (lgrad1) then
                    xderv(i) = xderv(i) - x21*e1d(1) - x31*e1d(2) - x41*e1d(3)
                    yderv(i) = yderv(i) - y21*e1d(1) - y31*e1d(2) - y41*e1d(3)
                    zderv(i) = zderv(i) - z21*e1d(1) - z31*e1d(2) - z41*e1d(3)
                    xderv(j) = xderv(j) - x32*e1d(4) + x21*e1d(1) - x42*e1d(5)
                    yderv(j) = yderv(j) - y32*e1d(4) + y21*e1d(1) - y42*e1d(5)
                    zderv(j) = zderv(j) - z32*e1d(4) + z21*e1d(1) - z42*e1d(5)
                    xderv(k) = xderv(k) + x32*e1d(4) - x43*e1d(6) + x31*e1d(2)
                    yderv(k) = yderv(k) + y32*e1d(4) - y43*e1d(6) + y31*e1d(2)
                    zderv(k) = zderv(k) + z32*e1d(4) - z43*e1d(6) + z31*e1d(2)
                    xderv(l) = xderv(l) + x43*e1d(6) + x42*e1d(5) + x41*e1d(3)
                    yderv(l) = yderv(l) + y43*e1d(6) + y42*e1d(5) + y41*e1d(3)
                    zderv(l) = zderv(l) + z43*e1d(6) + z42*e1d(5) + z41*e1d(3)
                  endif
                enddo
!  
!  End of loop over k
!             
              enddo kloop
!  
!  End of loops over cell images for k
!           
            enddo
          enddo
        enddo
!  
!  End of loop over j
!     
      enddo jloop
    endif
!  
!  End of ixyz loop
! 
  enddo
!
!  Closing banner for energy decomposition
!
  if (lPrintFour) then
    call mpbarrier
    if (ioproc) then
      write(ioout,'(''--------------------------------------------------------------------------------'',/)')
    endif
  endif
!
!  Free local memory
!
  deallocate(nptrnforimp,stat=status)
  if (status/=0) call deallocate_error('fourimpmds','nptrnforimp')
  deallocate(ntypmiddle,stat=status)
  if (status/=0) call deallocate_error('fourimpmds','ntypmiddle')
  deallocate(natmiddle,stat=status)
  if (status/=0) call deallocate_error('fourimpmds','natmiddle')
#ifdef TRACE
  call trace_out('fourimpmds')
#endif
!
!  All tidying up of derivatives is handled by four so we can just return here
!
  return
  end
