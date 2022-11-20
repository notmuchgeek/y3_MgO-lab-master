  subroutine sixmds(esix,esregion12,esregion2,eattach,lgrad1)
!
!  Subroutine for six-body potentials - first derivatives only 
!  using spatial decomposition
!
!   7/06 Created from sixmd.f
!   2/07 Bonding types and test added
!   5/07 QM/MM schemes added
!   6/07 Structure of arrays for storing spatial distribution changed to 1-D
!  12/07 Unused variables removed
!   4/08 Modified for variable domain size
!   4/08 xvec1cell replaced by xvec2cell etc for spatial algorithm
!   4/08 ind1toijk replaced by ind2toijk
!   6/09 Site energy and virials added
!  11/09 Region derivatives added
!  11/11 Region-region energy contributions stored
!   4/12 Explicit virial calculation removed as no longer needed
!   4/12 xvir, yvir and zvir removed
!   5/12 Atomic stresses added
!   1/17 lPrintSix added
!   1/18 Cell index shifts handled to allow for nomod
!   2/18 Trace added
!   9/18 Strain module introduced
!   9/18 Call to sixstrterms replaced with more general realstrterms
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
  use control,        only : lseok, latomicstress, lrigid
  use current
  use derivatives
  use energies,       only : siteenergy, eregion2region
  use iochannels,     only : ioout
  use m_strain,       only : realstrterms
  use mdlogic
  use molecule
  use numbers,        only : sixth, third
  use optimisation
  use parallel,       only : ioproc
  use six
  use spatial
  use symmetry
  use times,          only : tsix
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  real(dp), intent(inout)                      :: esix
  real(dp), intent(inout)                      :: esregion12
  real(dp), intent(inout)                      :: esregion2
  real(dp), intent(inout)                      :: eattach
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
  integer(i4)                                  :: ijcx
  integer(i4)                                  :: ijcy
  integer(i4)                                  :: ijcz
  integer(i4)                                  :: ikcx
  integer(i4)                                  :: ikcy
  integer(i4)                                  :: ikcz
  integer(i4)                                  :: ilcx
  integer(i4)                                  :: ilcy
  integer(i4)                                  :: ilcz
  integer(i4)                                  :: imx
  integer(i4)                                  :: imy
  integer(i4)                                  :: imz
  integer(i4)                                  :: ind
  integer(i4)                                  :: ind2
  integer(i4)                                  :: indm
  integer(i4)                                  :: indmj
  integer(i4)                                  :: indmk
  integer(i4)                                  :: indml
  integer(i4)                                  :: indmm
  integer(i4)                                  :: indmn
  integer(i4)                                  :: indn
  integer(i4)                                  :: indvec
  integer(i4)                                  :: isatom(6)
  integer(i4)                                  :: ivec
  integer(i4)                                  :: ix
  integer(i4)                                  :: ixm
  integer(i4)                                  :: ixn
  integer(i4)                                  :: iy
  integer(i4)                                  :: iym
  integer(i4)                                  :: iyn
  integer(i4)                                  :: iz
  integer(i4)                                  :: izm
  integer(i4)                                  :: izn
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
  integer(i4)                                  :: jj
  integer(i4)                                  :: jmcx
  integer(i4)                                  :: jmcy
  integer(i4)                                  :: jmcz
  integer(i4)                                  :: jncx
  integer(i4)                                  :: jncy
  integer(i4)                                  :: jncz
  integer(i4)                                  :: jmx
  integer(i4)                                  :: jmy
  integer(i4)                                  :: jmz
  integer(i4)                                  :: jndn
  integer(i4)                                  :: jvec
  integer(i4)                                  :: k
  integer(i4)                                  :: kc
  integer(i4)                                  :: kcx
  integer(i4)                                  :: kcy
  integer(i4)                                  :: kcz
  integer(i4)                                  :: kk
  integer(i4)                                  :: kl
  integer(i4)                                  :: ks
  integer(i4)                                  :: kmx
  integer(i4)                                  :: kmy
  integer(i4)                                  :: kmz
  integer(i4)                                  :: kndn
  integer(i4)                                  :: l
  integer(i4)                                  :: lc
  integer(i4)                                  :: lcx
  integer(i4)                                  :: lcy
  integer(i4)                                  :: lcz
  integer(i4)                                  :: ll
  integer(i4)                                  :: lmx
  integer(i4)                                  :: lmy
  integer(i4)                                  :: lmz
  integer(i4)                                  :: lndn
  integer(i4)                                  :: m
  integer(i4)                                  :: maxx
  integer(i4)                                  :: maxxy
  integer(i4)                                  :: mc
  integer(i4)                                  :: mcx
  integer(i4)                                  :: mcy
  integer(i4)                                  :: mcz
  integer(i4)                                  :: mm
  integer(i4)                                  :: mmx
  integer(i4)                                  :: mmy
  integer(i4)                                  :: mmz
  integer(i4)                                  :: mndn
  integer(i4)                                  :: n
  integer(i4)                                  :: n1i
  integer(i4)                                  :: n1j
  integer(i4)                                  :: n1k
  integer(i4)                                  :: n1l
  integer(i4)                                  :: n1m
  integer(i4)                                  :: n1n
  integer(i4)                                  :: nbtypeij
  integer(i4)                                  :: nbtypeik
  integer(i4)                                  :: nbtypeil
  integer(i4)                                  :: nbtypejm
  integer(i4)                                  :: nbtypejn
  integer(i4)                                  :: nbtypeij2
  integer(i4)                                  :: nbtypeik2
  integer(i4)                                  :: nbtypeil2
  integer(i4)                                  :: nbtypejm2
  integer(i4)                                  :: nbtypejn2
  integer(i4)                                  :: nc
  integer(i4)                                  :: ncx
  integer(i4)                                  :: ncy
  integer(i4)                                  :: ncz
  integer(i4)                                  :: neq
  integer(i4)                                  :: nsixtype
  integer(i4)                                  :: ni
  integer(i4)                                  :: nj
  integer(i4)                                  :: nk
  integer(i4)                                  :: nl 
  integer(i4)                                  :: nm 
  integer(i4)                                  :: nn 
  integer(i4)                                  :: nnn
  integer(i4)                                  :: nmi 
  integer(i4)                                  :: nmj 
  integer(i4)                                  :: nmk 
  integer(i4)                                  :: nml 
  integer(i4)                                  :: nmm 
  integer(i4)                                  :: nmn 
  integer(i4)                                  :: np
  integer(i4)                                  :: nr
  integer(i4)                                  :: nregioni
  integer(i4)                                  :: nregionj
  integer(i4)                                  :: nregionk
  integer(i4)                                  :: nregionl
  integer(i4)                                  :: nregionm
  integer(i4)                                  :: nregionn
  integer(i4)                                  :: nregiontypi
  integer(i4)                                  :: nregiontypj
  integer(i4)                                  :: nregiontypk
  integer(i4)                                  :: nregiontypl
  integer(i4)                                  :: nregiontypm
  integer(i4)                                  :: nregiontypn
  integer(i4)                                  :: nsplower(3)
  integer(i4)                                  :: nspupper(3)
  integer(i4)                                  :: nt2
  integer(i4)                                  :: nt3
  integer(i4)                                  :: nt4
  integer(i4)                                  :: nt5
  integer(i4)                                  :: nt6
  integer(i4)                                  :: ntyp2
  integer(i4)                                  :: ntyp3
  integer(i4)                                  :: ntyp4
  integer(i4)                                  :: ntyp5
  integer(i4)                                  :: ntyp6
  integer(i4)                                  :: ntypi
  integer(i4)                                  :: ntypj
  integer(i4)                                  :: ntypk
  integer(i4)                                  :: ntypl
  integer(i4)                                  :: ntypm
  integer(i4)                                  :: ntypn
  integer(i4)                                  :: status
  logical                                      :: l2bonds
  logical                                      :: lattach
  logical                                      :: lbonded
  logical                                      :: lbtyp
  logical                                      :: linter_only
  logical                                      :: lintra_only
  logical                                      :: ljok
  logical                                      :: llok
  logical                                      :: lnok
  logical                                      :: lmatch
  logical                                      :: lmolok
  logical                                      :: lneedmol
  logical                                      :: lopi
  logical                                      :: lopj
  logical                                      :: lopk
  logical                                      :: lopl
  logical                                      :: lopm
  logical                                      :: lopn
  logical                                      :: lreg12
  logical                                      :: lreg2stet
  logical                                      :: lsamemol
  logical                                      :: lsg1
  logical                                      :: lslicei
  logical                                      :: lslicej
  logical                                      :: lslicek
  logical                                      :: lslicel
  logical                                      :: lslicem
  logical                                      :: lslicen
  logical                                      :: ltsym12
  logical                                      :: ltsym34
  logical                                      :: ltsym56
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: e1d(15)
  real(dp)                                     :: e2d(120)
  real(dp)                                     :: e3d(1)
  real(dp)                                     :: eterm
  real(dp)                                     :: oci
  real(dp)                                     :: ocj
  real(dp)                                     :: ock
  real(dp)                                     :: ocl 
  real(dp)                                     :: ocm 
  real(dp)                                     :: ocn 
  real(dp)                                     :: ofct 
  real(dp)                                     :: r212
  real(dp)                                     :: r312
  real(dp)                                     :: r412
  real(dp)                                     :: r522
  real(dp)                                     :: r622
  real(dp)                                     :: rksix
  real(dp)                                     :: rko
  real(dp)                                     :: dr2ds(6,15)
  real(dp)                                     :: d2r2dx2(3,3,15)
  real(dp)                                     :: d2r2ds2(6,6,15)
  real(dp)                                     :: d2r2dsdx(6,3,15)
  real(dp)                                     :: rstrdloc(6)
  real(dp)                                     :: time1
  real(dp)                                     :: time2
  real(dp)                                     :: tr1
  real(dp)                                     :: tr2
  real(dp)                                     :: tr3
  real(dp)                                     :: tr4
  real(dp)                                     :: tr5
  real(dp)                                     :: ttr2
  real(dp)                                     :: ttr3
  real(dp)                                     :: ttr4
  real(dp)                                     :: ttr5
  real(dp)                                     :: x21
  real(dp)                                     :: y21
  real(dp)                                     :: z21
  real(dp)                                     :: x31
  real(dp)                                     :: y31
  real(dp)                                     :: z31
  real(dp)                                     :: x41
  real(dp)                                     :: y41
  real(dp)                                     :: z41
  real(dp)                                     :: x52
  real(dp)                                     :: y52
  real(dp)                                     :: z52
  real(dp)                                     :: x62
  real(dp)                                     :: y62
  real(dp)                                     :: z62
  real(dp)                                     :: sdist(15)
  real(dp)                                     :: svec(3,6,6)
  real(dp)                                     :: sxyz(3,6)
  real(dp)                                     :: xcom(15)
  real(dp)                                     :: ycom(15)
  real(dp)                                     :: zcom(15)
  real(dp)                                     :: xcomi
  real(dp)                                     :: ycomi
  real(dp)                                     :: zcomi
  real(dp)                                     :: xv6(15)
  real(dp)                                     :: yv6(15)
  real(dp)                                     :: zv6(15)
  real(dp), dimension(:), allocatable          :: xderv
  real(dp), dimension(:), allocatable          :: yderv
  real(dp), dimension(:), allocatable          :: zderv
#ifdef TRACE
  call trace_in('sixmds')
#endif
!
  time1 = g_cpu_time()
  lsg1 = (lstr.and.lgrad1)
!     
!  Allocate local memory                           
!
  allocate(xderv(numat),stat=status)                       
  if (status/=0) call outofmemory('sixmds','xderv')
  allocate(yderv(numat),stat=status)                       
  if (status/=0) call outofmemory('sixmds','yderv')
  allocate(zderv(numat),stat=status)
  if (status/=0) call outofmemory('sixmds','zderv')
!
!  Zero derivatives
!
  do i = 1,numat
    xderv(i) = 0.0_dp
    yderv(i) = 0.0_dp
    zderv(i) = 0.0_dp
  enddo
!
!  Initialise centre of mass arrays to zero since there should be no intermolecular terms here
!
  xcom(1:15) = 0.0_dp
  ycom(1:15) = 0.0_dp
  zcom(1:15) = 0.0_dp
!
!  Opening banner for energy decomposition
!
  if (lPrintSix) then
    call mpbarrier
    if (ioproc) then
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''  Six  :  Mid Atom 1 + End 3/4    Mid Atom 2 + End 5/6     Sixbody energy (eV)  '')')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
    endif
  endif
!
!  Set variables for cell distribution
!
  maxxy = nspcell(1)*nspcell(2)
  maxx  = nspcell(1)
!*************************
!  Loop over potentials  *
!*************************
  do np = 1,nsix
    nsixtype = nsixty(np)
    tr1 = six1(np)**2
    ttr2 = six2(np)**2
    ttr3 = six3(np)**2
    ttr4 = six4(np)**2
    ttr5 = six5(np)**2
    ltsym12 = (lmatch(nsspec1(np),nsptyp1(np),nsspec2(np),nsptyp2(np),.true.).or. &
               lmatch(nsspec2(np),nsptyp2(np),nsspec1(np),nsptyp1(np),.true.))
    lbtyp = (mmsexc(np).eq.1)
    rksix = sixk(np)
    lintra_only = (lsintra(np).and..not.lsinter(np))
    linter_only = (lsinter(np).and..not.lsintra(np))
    lneedmol = (lintra_only.or.linter_only.or.lbtyp.or.lrigid)
!********************************
!  Loop over middle site 1 / i  *
!********************************
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
!  Outer loop over atoms within this cell
!           
        n1i = nspcellat1ptr(ind)
        iloop: do ii = 1,nspcellat(ind)
          i = nspcellatptr(n1i+ii)
          ni = nat(i)
          ntypi = nftype(i)
          nregioni = nregionno(nsft+nrelf2a(i))
          nregiontypi = nregiontype(nregioni,ncf)
          oci = occuf(i)
          lopi = (lopf(nrelf2a(i)).or..not.lfreeze)
          lslicei = lsliceatom(nsft + nrelf2a(i))
!
!  Check i is allowed for n
!
          if (lmatch(ni,ntypi,nsspec1(np),nsptyp1(np),.true.)) then
            nt2 = nsspec2(np)
            ntyp2 = nsptyp2(np)
            nt3 = nsspec3(np)
            nt4 = nsspec4(np)
            nt5 = nsspec5(np)
            nt6 = nsspec6(np)
            ntyp3 = nsptyp3(np)
            ntyp4 = nsptyp4(np)
            ntyp5 = nsptyp5(np)
            ntyp6 = nsptyp6(np)
            tr2 = ttr2
            tr3 = ttr3
            tr4 = ttr4
            tr5 = ttr5
          elseif (lmatch(ni,ntypi,nsspec2(np),nsptyp2(np),.true.)) then
            nt2 = nsspec1(np)
            ntyp2 = nsptyp1(np)
            nt3 = nsspec5(np)
            nt4 = nsspec6(np)
            nt5 = nsspec3(np)
            nt6 = nsspec4(np)
            ntyp3 = nsptyp5(np)
            ntyp4 = nsptyp6(np)
            ntyp5 = nsptyp3(np)
            ntyp6 = nsptyp4(np)
            tr2 = ttr4
            tr3 = ttr5
            tr4 = ttr2
            tr5 = ttr3
          else
            cycle iloop
          endif
          ic = nspcellatptrcell(n1i+ii)
          call ind2toijk(ic,icx,icy,icz)
!
!  Set flags as to whether potentials are symmetric at the two ends
!
          ltsym34 = (lmatch(nt3,ntyp3,nt4,ntyp4,.true.).or.lmatch(nt4,ntyp4,nt3,ntyp3,.true.))
          ltsym56 = (lmatch(nt5,ntyp5,nt6,ntyp6,.true.).or.lmatch(nt6,ntyp6,nt5,ntyp5,.true.))
!
!  i has been accepted
!
          isatom(1) = i
          sxyz(1,1) = xclat(i)
          sxyz(2,1) = yclat(i)
          sxyz(3,1) = zclat(i)
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
!  Loop over neighbouring cells
!
          do imz = nsplower(3),nspupper(3)
            do imy = nsplower(2),nspupper(2)
              do imx = nsplower(1),nspupper(1)
                indn = (imz-1)*maxxy + (imy-1)*maxx + imx
!********************************
!  Loop over middle site 2 / j  *
!********************************
                n1j = nspcellat1ptr(indn)
                jloop: do jj = 1,nspcellat(indn)
                  j = nspcellatptr(n1j+jj)
                  if (ltsym12) then
                    ljok = (j.lt.i)
                  else
                    ljok = (j.ne.i)
                  endif
!
!  If j is equal to i or excluded because the potential is symmetric then skip
!
                  if (ljok) then
                    nj = nat(j)
                    ntypj = nftype(j)
!
!  Check j is allowed for n
!
                    if (.not.lmatch(nj,ntypj,nt2,ntyp2,.true.)) cycle jloop
                    if (lmol.and.lneedmol) then
!
!  Molecule handling
!
                      nmj = natmol(j)
                      if (ndim.gt.0) then
                        indmj = nmolind(j)
                        call mindtoijk(indmj,ixj,iyj,izj)
!
                        ixj = ixj - ixi
                        iyj = iyj - iyi
                        izj = izj - izi
                      endif
                      lmolok = (nmi.eq.nmj.and.nmi.gt.0)
!
!  Set COM coordinates
!
                      if (lrigid.and.nmj.gt.0.and.lmolok) then
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
!  Set properties for j
!  
                    jc = nspcellatptrcell(n1j+jj)
                    call ind2toijk(jc,jcx,jcy,jcz)
                    sxyz(1,2) = xvec2cell(jc) + xinbox(j)
                    sxyz(2,2) = yvec2cell(jc) + yinbox(j)
                    sxyz(3,2) = zvec2cell(jc) + zinbox(j)
                    x21 = sxyz(1,2) - sxyz(1,1)
                    y21 = sxyz(2,2) - sxyz(2,1)
                    z21 = sxyz(3,2) - sxyz(3,1)
!
!  Check r21 is OK
!
                    r212 = x21*x21 + y21*y21 + z21*z21
                    if (r212.lt.1d-12) cycle jloop
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
!
!  Check central bond type for correct order
!
                          if (n6botype(1,np).gt.0) then
                            if (n6botype(1,np).ne.nbtypeij) cycle jloop
                          endif
                          if (n6botype(2,np).gt.0) then
                            if (n6botype(2,np).ne.nbtypeij2) cycle jloop
                          endif
                        endif
                      else
                        ijcx = jcx - icx + icosxsp(j) - icosxsp(i)
                        ijcy = jcy - icy + icosysp(j) - icosysp(i)
                        ijcz = jcz - icz + icoszsp(j) - icoszsp(i)
!
                        if (lbtyp) then
                          call bonded(lbonded,l2bonds,nbtypeij,nbtypeij2,i,j,ijcx,ijcy,ijcz)
                          if (.not.lbonded) cycle jloop
!
!  Check central bond type for correct order
!
                          if (n6botype(1,np).gt.0) then
                            if (n6botype(1,np).ne.nbtypeij) cycle jloop
                          endif
                          if (n6botype(2,np).gt.0) then
                            if (n6botype(2,np).ne.nbtypeij2) cycle jloop
                          endif
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
                    if (r212.gt.tr1.and.(.not.lbtyp.or..not.lbonded)) cycle jloop
!
!  j has been accepted
!
                    isatom(2) = j
!
!  Set remaining properties for j
!
                    nregionj = nregionno(nsft+nrelf2a(j))
                    nregiontypj = nregiontype(nregionj,ncf)
                    ocj = occuf(j)
                    lopj = (lopf(nrelf2a(j)).or..not.lfreeze)
                    lslicej = lsliceatom(nsft + nrelf2a(j))
!
!  QM/MM handling : i & j are both QM atoms and potential is of bonded type => exclude
!
                    if (QMMMmode(ncf).gt.0) then
                      if (nregiontypi.eq.1.and.nregiontypj.eq.1.and.lbtyp) cycle jloop
                    endif
!*****************************************
!  Loop over first end site bonded to i  *
!*****************************************
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
                          n1k = nspcellat1ptr(jndn)
                          kloop: do kk = 1,nspcellat(jndn)
                            k = nspcellatptr(n1k+kk)
!                     
!  Check k is allowed for n
!                         
                            if (k.eq.j) cycle kloop
!                           
!  Prevent atoms i and k being the same atom
!                         
                            if (k.eq.i) cycle kloop
                            nk = nat(k)
                            ntypk = nftype(k)
!
!  Check k is allowed for n
!
                            if (.not.lmatch(nk,ntypk,nt3,ntyp3,.true.)) cycle kloop
                            if (lmol.and.lneedmol) then
!
!  Molecule handling
!
                              nmk = natmol(k)
                              if (ndim.gt.0) then
                                indmk = nmolind(k)
                                call mindtoijk(indmk,ixk,iyk,izk)
!
                                ixk = ixk - ixi
                                iyk = iyk - iyi
                                izk = izk - izi
                              endif
                              lmolok = (nmi.eq.nmk.and.nmi.gt.0)
!
!  Set COM coordinates
!
                              if (lrigid.and.nmk.gt.0.and.lmolok) then
                                xcom(2) = molxyz(1,natinmol(k),nmk) - xcomi
                                ycom(2) = molxyz(2,natinmol(k),nmk) - ycomi
                                zcom(2) = molxyz(3,natinmol(k),nmk) - zcomi
                              else
                                xcom(2) = - xcomi
                                ycom(2) = - ycomi
                                zcom(2) = - zcomi
                              endif
                              xcom(6) = xcom(2) - xcom(1)
                              ycom(6) = ycom(2) - ycom(1)
                              zcom(6) = zcom(2) - zcom(1)
                            else
                              lmolok = .false.
                              xcom(2) = - xcomi
                              ycom(2) = - ycomi
                              zcom(2) = - zcomi
                              xcom(6) = xcom(2) - xcom(1)
                              ycom(6) = ycom(2) - ycom(1)
                              zcom(6) = zcom(2) - zcom(1)
                            endif
!
!  Check for intra and but not in same molecule
!
                            if (lintra_only.and..not.lmolok) cycle kloop
                            if (lbtyp.and..not.lmolok) cycle kloop
!
!  Set properties of atom k
!
                            kc = nspcellatptrcell(n1k+kk)
                            call ind2toijk(kc,kcx,kcy,kcz)
                            sxyz(1,3) = xvec2cell(kc) + xinbox(k)
                            sxyz(2,3) = yvec2cell(kc) + yinbox(k)
                            sxyz(3,3) = zvec2cell(kc) + zinbox(k)
                            x31 = sxyz(1,3) - sxyz(1,1)
                            y31 = sxyz(2,3) - sxyz(2,1)
                            z31 = sxyz(3,3) - sxyz(3,1)
!
!  Check r31 is OK
!
                            r312 = x31*x31 + y31*y31 + z31*z31
                            if (r312.lt.1d-12) cycle kloop
!
!  Molecule checking
!
                            lbonded = .false.
                            if (lmolok) then
                              if (ndim.eq.0) then
                                if (linter_only) cycle kloop
                                if (lbtyp) then
                                  call bonded(lbonded,l2bonds,nbtypeik,nbtypeik2,i,k,0_i4,0_i4,0_i4)
                                  if (.not.lbonded) cycle kloop
                                endif
                              else
                                ikcx = kcx - icx + icosxsp(k) - icosxsp(i)
                                ikcy = kcy - icy + icosysp(k) - icosysp(i)
                                ikcz = kcz - icz + icoszsp(k) - icoszsp(i)
!
                                if (lbtyp) then
                                  call bonded(lbonded,l2bonds,nbtypeik,nbtypeik2,i,k,ikcx,ikcy,ikcz)
                                  if (.not.lbonded) cycle kloop
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
                            if (r312.gt.tr2.and.(.not.lbtyp.or..not.lbonded)) cycle kloop
!
!  k has been accepted
!
                            isatom(3) = k
!
!  Set remaining properties of atom k
!
                            nregionk = nregionno(nsft+nrelf2a(k))
                            nregiontypk = nregiontype(nregionk,ncf)
                            ock = occuf(k)
                            lopk = (lopf(nrelf2a(k)).or..not.lfreeze)
                            lslicek = lsliceatom(nsft + nrelf2a(k))
!************************************
!  Loop over second end site for i  *
!************************************
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
                                  n1l = nspcellat1ptr(kndn)
                                  lloop: do ll = 1,nspcellat(kndn)
                                    l = nspcellatptr(n1l+ll)
!
!  Check whether l is allowed based on symmetry of pair of atoms attached to i
!
                                    if (ltsym34) then
                                      llok = (l.lt.k)
                                    else
                                      llok = .true.
                                    endif
                                    if (llok) then
                                      nl = nat(l)
                                      ntypl = nftype(l)
!
!  Prevent atoms i and l being the same atom
!
                                      if (l.eq.i) cycle lloop
!
!  Prevent atoms j and l being the same atom
!
                                      if (l.eq.j) cycle lloop
!
!  Prevent atoms k and l being the same atom
!
                                      if (l.eq.k) cycle lloop
!
!  Check l is allowed for n
!
                                      if (.not.lmatch(nl,ntypl,nt4,ntyp4,.true.)) cycle lloop
!
!  Molecule handling
!
                                      if (lmol.and.lneedmol) then
                                        nml = natmol(l)
                                        if (ndim.gt.0) then
                                          indml = nmolind(l)
                                          call mindtoijk(indml,ixl,iyl,izl)
!
                                          ixl = ixl - ixi
                                          iyl = iyl - iyi
                                          izl = izl - izi
                                        endif
                                        lmolok = (nmi.eq.nml.and.nmi.gt.0)
!
!  Set COM coordinates
!
                                        if (lrigid.and.nml.gt.0.and.lmolok) then
                                          xcom(3) = molxyz(1,natinmol(l),nml) - xcomi
                                          ycom(3) = molxyz(2,natinmol(l),nml) - ycomi
                                          zcom(3) = molxyz(3,natinmol(l),nml) - zcomi
                                        else
                                          xcom(3) = - xcomi
                                          ycom(3) = - ycomi
                                          zcom(3) = - zcomi
                                        endif
                                        xcom(7) = xcom(3) - xcom(1)
                                        ycom(7) = ycom(3) - ycom(1)
                                        zcom(7) = zcom(3) - zcom(1)
                                        xcom(10) = xcom(3) - xcom(2)
                                        ycom(10) = ycom(3) - ycom(2)
                                        zcom(10) = zcom(3) - zcom(2)
                                      else
                                        lmolok = .false.
                                        xcom(3) = - xcomi
                                        ycom(3) = - ycomi
                                        zcom(3) = - zcomi
                                        xcom(7) = xcom(3) - xcom(1)
                                        ycom(7) = ycom(3) - ycom(1)
                                        zcom(7) = zcom(3) - zcom(1)
                                        xcom(10) = xcom(3) - xcom(2)
                                        ycom(10) = ycom(3) - ycom(2)
                                        zcom(10) = zcom(3) - zcom(2)
                                      endif
!
!  Check for intra and but not in same molecule
!
                                      if (lintra_only.and..not.lmolok) cycle lloop
                                      if (lbtyp.and..not.lmolok) cycle lloop
!
!  Set properties of l atom
!
                                      lc = nspcellatptrcell(n1l+ll)
                                      call ind2toijk(lc,lcx,lcy,lcz)
                                      sxyz(1,4) = xvec2cell(lc) + xinbox(l)
                                      sxyz(2,4) = yvec2cell(lc) + yinbox(l)
                                      sxyz(3,4) = zvec2cell(lc) + zinbox(l)
                                      x41 = sxyz(1,4) - sxyz(1,1)
                                      y41 = sxyz(2,4) - sxyz(2,1)
                                      z41 = sxyz(3,4) - sxyz(3,1)
!
!  Check r41 is OK
!
                                      r412 = x41*x41 + y41*y41 + z41*z41
                                      if (r412.lt.1d-12) cycle lloop
!
!  Molecule checking
!
                                      lbonded = .false.
                                      if (lmolok) then
                                        if (ndim.eq.0) then
                                          if (linter_only) cycle lloop
                                          if (lbtyp) then
                                            call bonded(lbonded,l2bonds,nbtypeil,nbtypeil2,i,l,0_i4,0_i4,0_i4)
                                            if (.not.lbonded) cycle lloop
                                          endif
                                        else
                                          ilcx = lcx - icx + icosxsp(l) - icosxsp(i)
                                          ilcy = lcy - icy + icosysp(l) - icosysp(i)
                                          ilcz = lcz - icz + icoszsp(l) - icoszsp(i)
!
                                          if (lbtyp) then
                                            call bonded(lbonded,l2bonds,nbtypeil,nbtypeil2,i,l,ilcx,ilcy,ilcz)
                                            if (.not.lbonded) cycle lloop
                                            lsamemol = (lbonded.or.l2bonds)
                                          else
                                            lsamemol = .false.
                                          endif
                                          if (.not.lsamemol) then
                                            call samemol(lsamemol,nmi,ilcx,ilcy,ilcz,ixl,iyl,izl)
                                          endif
                                          if (lintra_only.and..not.lsamemol) cycle lloop
                                          if (linter_only.and.lsamemol) cycle lloop
                                        endif
                                      endif
!
!  Distance checking
!
                                      if (r412.gt.tr3.and.(.not.lbtyp.or..not.lbonded)) cycle lloop
!
!  l has been accepted
!
                                      isatom(4) = l
!
!  Set remaining properties of l
!
                                      nregionl = nregionno(nsft+nrelf2a(l))
                                      nregiontypl = nregiontype(nregionl,ncf)
                                      ocl = occuf(l)
                                      lopl = (lopf(nrelf2a(l)).or..not.lfreeze)
                                      lslicel = lsliceatom(nsft + nrelf2a(l))
!*****************************************
!  Loop over first end site bonded to j  *
!*****************************************
!
!  Loop over neighbouring cells for l
!
                                      do lmz = nsplower(3),nspupper(3)
                                        do lmy = nsplower(2),nspupper(2)
                                          do lmx = nsplower(1),nspupper(1)
                                            lndn = (lmz-1)*maxxy + (lmy-1)*maxx + lmx
!
!  Loop over atoms within neighbouring cells
!
                                            n1m = nspcellat1ptr(lndn)
                                            mloop: do mm = 1,nspcellat(lndn)
                                              m = nspcellatptr(n1m+mm)
                                              nm = nat(m)
                                              ntypm = nftype(m)
!
!  Prevent atoms i and m being the same atom
!
                                              if (m.eq.i) cycle mloop
!
!  Prevent atoms j and m being the same atom
!
                                              if (m.eq.j) cycle mloop
!
!  Check m is allowed for n
!
                                              if (.not.lmatch(nm,ntypm,nt4,ntyp4,.true.)) cycle mloop
                                              if (lmol.and.lneedmol) then
!
!  Molecule handling
!
                                                nmm = natmol(m)
                                                if (ndim.gt.0) then
                                                  indmm = nmolind(m)
                                                  call mindtoijk(indmm,ixm,iym,izm)
!
                                                  ixm = ixm - ixj
                                                  iym = iym - iyj
                                                  izm = izm - izj
                                                endif
                                                lmolok = (nmi.eq.nmm.and.nmi.gt.0)
!
!  Set COM coordinates
!
                                                if (lrigid.and.nmm.gt.0.and.lmolok) then
                                                  xcom(4) = molxyz(1,natinmol(m),nmm) - xcomi
                                                  ycom(4) = molxyz(2,natinmol(m),nmm) - ycomi
                                                  zcom(4) = molxyz(3,natinmol(m),nmm) - zcomi
                                                else
                                                  xcom(4) = - xcomi
                                                  ycom(4) = - ycomi
                                                  zcom(4) = - zcomi
                                                endif
                                                xcom(8) = xcom(4) - xcom(1)
                                                ycom(8) = ycom(4) - ycom(1)
                                                zcom(8) = zcom(4) - zcom(1)
                                                xcom(11) = xcom(4) - xcom(2)
                                                ycom(11) = ycom(4) - ycom(2)
                                                zcom(11) = zcom(4) - zcom(2)
                                                xcom(13) = xcom(4) - xcom(3)
                                                ycom(13) = ycom(4) - ycom(3)
                                                zcom(13) = zcom(4) - zcom(3)
                                              else
                                                lmolok = .false.
                                                xcom(4) = - xcomi
                                                ycom(4) = - ycomi
                                                zcom(4) = - zcomi
                                                xcom(8) = xcom(4) - xcom(1)
                                                ycom(8) = ycom(4) - ycom(1)
                                                zcom(8) = zcom(4) - zcom(1)
                                                xcom(11) = xcom(4) - xcom(2)
                                                ycom(11) = ycom(4) - ycom(2)
                                                zcom(11) = zcom(4) - zcom(2)
                                                xcom(13) = xcom(4) - xcom(3)
                                                ycom(13) = ycom(4) - ycom(3)
                                                zcom(13) = zcom(4) - zcom(3)
                                              endif
!
!  Check for intra and but not in same molecule
!
                                              if (lintra_only.and..not.lmolok) cycle mloop
                                              if (lbtyp.and..not.lmolok) cycle mloop
!
!  Set properties of l atom
!
                                              mc = nspcellatptrcell(n1m+mm)
                                              call ind2toijk(mc,mcx,mcy,mcz)
                                              sxyz(1,5) = xvec2cell(mc) + xinbox(m)
                                              sxyz(2,5) = yvec2cell(mc) + yinbox(m)
                                              sxyz(3,5) = zvec2cell(mc) + zinbox(m)
                                              x52 = sxyz(1,5) - sxyz(1,2)
                                              y52 = sxyz(2,5) - sxyz(2,2)
                                              z52 = sxyz(3,5) - sxyz(3,2)
!
!  Check r52 is OK
!
                                              r522 = x52*x52 + y52*y52 + z52*z52
                                              if (r522.lt.1d-12) cycle mloop
!
!  Molecule checking
!
                                              lbonded = .false.
                                              if (lmolok) then
                                                if (ndim.eq.0) then
                                                  if (linter_only) cycle mloop
                                                  if (lbtyp) then
                                                    call bonded(lbonded,l2bonds,nbtypejm,nbtypejm2,j,m,0_i4,0_i4,0_i4)
                                                    if (.not.lbonded) cycle mloop
                                                  endif
                                                else
                                                  jmcx = mcx - jcx + icosxsp(m) - icosxsp(j)
                                                  jmcy = mcy - jcy + icosysp(m) - icosysp(j)
                                                  jmcz = mcz - jcz + icoszsp(m) - icoszsp(j)
!
                                                  if (lbtyp) then
                                                    call bonded(lbonded,l2bonds,nbtypejm,nbtypejm2,j,m,jmcx,jmcy,jmcz)
                                                    if (.not.lbonded) cycle mloop
                                                    lsamemol = (lbonded.or.l2bonds)
                                                  else
                                                    lsamemol = .false.
                                                  endif
                                                  if (.not.lsamemol) then
                                                    call samemol(lsamemol,nmi,jmcx,jmcy,jmcz,ixm,iym,izm)
                                                  endif
                                                  if (lintra_only.and..not.lsamemol) cycle mloop
                                                  if (linter_only.and.lsamemol) cycle mloop
                                                endif
                                              endif
!
!  Distance checking
!
                                              if (r522.gt.tr4.and.(.not.lbtyp.or..not.lbonded)) cycle mloop
!
!  m has been accepted
!
                                              isatom(5) = m
                                              nregionm = nregionno(nsft+nrelf2a(m))
                                              nregiontypm = nregiontype(nregionm,ncf)
                                              ocm = occuf(m)
                                              lopm = (lopf(nrelf2a(m)).or..not.lfreeze)
                                              lslicem = lsliceatom(nsft + nrelf2a(m))
!************************************
!  Loop over second end site for j  *
!************************************
!
!  Loop over neighbouring cells for l
!
                                              do mmz = nsplower(3),nspupper(3)
                                                do mmy = nsplower(2),nspupper(2)
                                                  do mmx = nsplower(1),nspupper(1)
                                                    mndn = (mmz-1)*maxxy + (mmy-1)*maxx + mmx
!
!  Loop over atoms within neighbouring cells
!
                                                    n1n = nspcellat1ptr(lndn)
                                                    nloop: do nn = 1,nspcellat(mndn)
                                                      n = nspcellatptr(n1n+nn)
                                                      if (ltsym56) then
                                                        lnok = (n.lt.m)
                                                      else
                                                        lnok = .true.
                                                      endif
                                                      if (lnok) then
                                                        nnn = nat(n)
                                                        ntypn = nftype(n)
!
!  Prevent atoms i and n being the same atom
!
                                                        if (n.eq.i) cycle nloop
!
!  Prevent atoms j and n being the same atom
!
                                                        if (n.eq.j) cycle nloop
!
!  Prevent atoms m and n being the same atom
!
                                                        if (n.eq.m) cycle nloop
!
!  If lfreeze=.true. and no atoms have any variables then skip this four body term
!
                                                        if (.not.lopi.and..not.lopj.and..not.lopk.and..not.lopl.and. &
                                                            .not.lopm.and..not.lopn) cycle nloop
!
!  QM/MM handling : i, j, k, l, m and n are all QM atoms => exclude
!
                                                        if (QMMMmode(ncf).gt.0) then
                                                          if (nregiontypi.eq.1.and.nregiontypj.eq.1.and.nregiontypk.eq.1.and. &
                                                              nregiontypl.eq.1.and.nregiontypm.eq.1.and.nregiontypn.eq.1) &
                                                            cycle nloop
                                                        endif
!
!  Check n is allowed for n
!
                                                        if (.not.lmatch(nnn,ntypn,nt6,ntyp6,.true.)) cycle nloop
!
!  Molecule handling
!
                                                        if (lmol.and.lneedmol) then
                                                          nmn = natmol(n)
                                                          if (ndim.gt.0) then
                                                            indmn = nmolind(n)
                                                            call mindtoijk(indmn,ixn,iyn,izn)
!
                                                            ixn = ixn - ixj
                                                            iyn = iyn - iyj
                                                            izn = izn - izj
                                                          endif
                                                          lmolok = (nmi.eq.nmn.and.nmi.gt.0)
!
!  Set COM coordinates
!
                                                          if (lrigid.and.nmn.gt.0.and.lmolok) then
                                                            xcom(5) = molxyz(1,natinmol(n),nmn) - xcomi
                                                            ycom(5) = molxyz(2,natinmol(n),nmn) - ycomi
                                                            zcom(5) = molxyz(3,natinmol(n),nmn) - zcomi
                                                          else
                                                            xcom(5) = - xcomi
                                                            ycom(5) = - ycomi
                                                            zcom(5) = - zcomi
                                                          endif
                                                          xcom(9) = xcom(5) - xcom(1)
                                                          ycom(9) = ycom(5) - ycom(1)
                                                          zcom(9) = zcom(5) - zcom(1)
                                                          xcom(12) = xcom(5) - xcom(2)
                                                          ycom(12) = ycom(5) - ycom(2)
                                                          zcom(12) = zcom(5) - zcom(2)
                                                          xcom(14) = xcom(5) - xcom(3)
                                                          ycom(14) = ycom(5) - ycom(3)
                                                          zcom(14) = zcom(5) - zcom(3)
                                                          xcom(15) = xcom(5) - xcom(4)
                                                          ycom(15) = ycom(5) - ycom(4)
                                                          zcom(15) = zcom(5) - zcom(4)
                                                        else
                                                          lmolok = .false.
                                                          xcom(5) = - xcomi
                                                          ycom(5) = - ycomi
                                                          zcom(5) = - zcomi
                                                          xcom(9) = xcom(5) - xcom(1)
                                                          ycom(9) = ycom(5) - ycom(1)
                                                          zcom(9) = zcom(5) - zcom(1)
                                                          xcom(12) = xcom(5) - xcom(2)
                                                          ycom(12) = ycom(5) - ycom(2)
                                                          zcom(12) = zcom(5) - zcom(2)
                                                          xcom(14) = xcom(5) - xcom(3)
                                                          ycom(14) = ycom(5) - ycom(3)
                                                          zcom(14) = zcom(5) - zcom(3)
                                                          xcom(15) = xcom(5) - xcom(4)
                                                          ycom(15) = ycom(5) - ycom(4)
                                                          zcom(15) = zcom(5) - zcom(4)
                                                        endif
!
!  Check for intra and but not in same molecule
!
                                                        if (lintra_only.and..not.lmolok) cycle nloop
                                                        if (lbtyp.and..not.lmolok) cycle nloop
                                                        nc = nspcellatptrcell(n1n+nn)
                                                        call ind2toijk(nc,ncx,ncy,ncz)
                                                        sxyz(1,6) = xvec2cell(nc) + xinbox(n)
                                                        sxyz(2,6) = yvec2cell(nc) + yinbox(n)
                                                        sxyz(3,6) = zvec2cell(nc) + zinbox(n)
                                                        x62 = sxyz(1,6) - sxyz(1,2)
                                                        y62 = sxyz(2,6) - sxyz(2,2)
                                                        z62 = sxyz(3,6) - sxyz(3,2)
!
!  Check r62 is OK
!
                                                        r622 = x62*x62 + y62*y62 + z62*z62
                                                        if (r622.lt.1d-12) cycle nloop
!
!  Molecule checking
!
                                                        lbonded = .false.
                                                        if (lmolok) then
                                                          if (ndim.eq.0) then
                                                            if (linter_only) cycle nloop
                                                            if (lbtyp) then
                                                              call bonded(lbonded,l2bonds,nbtypejn,nbtypejn2,j,n,0_i4,0_i4,0_i4)
                                                              if (.not.lbonded) cycle nloop
                                                            endif
                                                          else
                                                            jncx = ncx - jcx + icosxsp(n) - icosxsp(j)
                                                            jncy = ncy - jcy + icosysp(n) - icosysp(j)
                                                            jncz = ncz - jcz + icoszsp(n) - icoszsp(j)
!
                                                            if (lbtyp) then
                                                              call bonded(lbonded,l2bonds,nbtypejn,nbtypejn2,j,n,jncx,jncy,jncz)
                                                              if (.not.lbonded) cycle nloop
                                                              lsamemol = (lbonded.or.l2bonds)
                                                            else
                                                              lsamemol = .false.
                                                            endif
                                                            if (.not.lsamemol) then
                                                              call samemol(lsamemol,nmi,jncx,jncy,jncz,ixn,iyn,izn)
                                                            endif
                                                            if (lintra_only.and..not.lsamemol) cycle nloop
                                                            if (linter_only.and.lsamemol) cycle nloop
                                                          endif
                                                        endif
!
!  Distance checking
!
                                                        if (r622.gt.tr5.and.(.not.lbtyp.or..not.lbonded)) cycle nloop
!
!  n has been accepted
!
                                                        isatom(6) = n
                                                        nregionn = nregionno(nsft+nrelf2a(n))
                                                        nregiontypn = nregiontype(nregionn,ncf)
                                                        ocn = occuf(n)
                                                        lopn = (lopf(nrelf2a(n)).or..not.lfreeze)
                                                        lslicen = lsliceatom(nsft + nrelf2a(n))
!     
!  Set region 2 sextet flag
!     
                                                        lreg12    = .false.
                                                        lreg2stet = .false.
                                                        if (lseok.and.nregions(ncf).gt.1) then
                                                          lreg2stet = (nregioni.gt.1.and.nregionj.gt.1.and.nregionk.gt.1.and. &
                                                                       nregionl.gt.1.and.nregionm.gt.1.and.nregionn.gt.1)
                                                          if (.not.lreg2stet) lreg12 = (nregioni.gt.1.or.nregionj.gt.1.or. &
                                                                                        nregionk.gt.1.or.nregionl.gt.1.or. &
                                                                                        nregionl.gt.1.or.nregionn.gt.1)
                                                        endif
                                                        lattach = .true.
                                                        if (lslicei.and.lslicej.and.lslicek.and.lslicel.and.lslicem.and.lslicen)  &
                                                          lattach = .false.
                                                        if (.not.lslicei.and..not.lslicej.and..not.lslicek.and..not.lslicel.and. &
                                                            .not.lslicem.and..not.lslicen) lattach = .false.
!********************************
!  Valid six-body term located  *
!********************************
!
!  Calculate vectors and remaining distances
!
                                                        do ivec = 1,6
                                                          do jvec = 1,6
                                                            svec(1,jvec,ivec) = sxyz(1,ivec) - sxyz(1,jvec)
                                                            svec(2,jvec,ivec) = sxyz(2,ivec) - sxyz(2,jvec)
                                                            svec(3,jvec,ivec) = sxyz(3,ivec) - sxyz(3,jvec)
                                                          enddo
                                                        enddo
                                                        indvec = 0
                                                        do ivec = 1,5
                                                          do jvec = ivec+1,6
                                                            indvec = indvec + 1
                                                            xv6(indvec) = svec(1,jvec,ivec)
                                                            yv6(indvec) = svec(2,jvec,ivec)
                                                            zv6(indvec) = svec(3,jvec,ivec)
                                                            sdist(indvec) = svec(1,jvec,ivec)**2 +  &
                                                                            svec(2,jvec,ivec)**2 +  &
                                                                            svec(3,jvec,ivec)**2
                                                            sdist(indvec) = sqrt(sdist(indvec))
                                                          enddo
                                                        enddo
!
                                                        ofct = oci*ocj*ock*ocl*ocm*ocn
                                                        rko = rksix*ofct
                                                        call sixbody(np,nsixtype,sdist,eterm,e1d,e2d,e3d,rko,lgrad1, &
                                                                     .false.,.false.)
                                                        if (lreg2stet) then
                                                          esregion2 = esregion2 + eterm
                                                        elseif (lreg12) then
                                                          esregion12 = esregion12 + eterm
                                                        else
                                                          esix = esix + eterm
                                                        endif
                                                        if (lattach) eattach = eattach + eterm
!
!  Assign inter-region energy to i-j pair only - sixbody terms shouldn't cross regions anyway
!
                                                        eregion2region(nregionj,nregioni) = eregion2region(nregionj,nregioni) + &
                                                          eterm
!
                                                        do ivec = 1,6
                                                          siteenergy(isatom(ivec)) = siteenergy(isatom(ivec)) + 0.5_dp*third*eterm
                                                        enddo
!
!  Output energy contribution
!
                                                        if (lPrintSix) then
                                                          write(ioout,'(4x,6i8,1x,f22.10)') i,k,l,j,m,n,eterm
                                                        endif
!***********************************
!  Cross out of plane derivatives  *
!***********************************
!
!  Set up strain products
!
                                                        if (lsg1) then
                                                          call realstrterms(ndim,15_i4,15_i4,xv6,yv6,zv6,xcom,ycom,zcom, &
                                                            dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
                                                        endif
!***********************
!  Strain derivatives  *
!***********************
                                                        if (lsg1) then
!
!  First strain derivatives
!
                                                          rstrdloc(1:nstrains) = 0.0_dp
                                                          do ivec = 1,15
                                                            do kl = 1,nstrains
                                                              ks = nstrptr(kl)
                                                              rstrdloc(kl) = rstrdloc(kl) + e1d(ivec)*dr2ds(ks,ivec)
                                                            enddo
                                                          enddo
                                                          do kl = 1,nstrains
                                                            rstrd(kl) = rstrd(kl) + rstrdloc(kl)
                                                          enddo
                                                          if (latomicstress) then
                                                            do ivec = 1,6
                                                              do kl = 1,nstrains
                                                                atomicstress(kl,isatom(ivec)) = atomicstress(kl,isatom(ivec)) + &
                                                                  sixth*rstrdloc(kl)
                                                              enddo
                                                            enddo
                                                          endif
                                                        endif
!*************************
!  Internal derivatives  *
!*************************
                                                        if (lgrad1) then
                                                          indvec = 0
                                                          do ivec = 1,5
                                                            do jvec = ivec+1,6
                                                              indvec = indvec + 1
                                                              xderv(isatom(ivec)) = xderv(isatom(ivec)) +  &
                                                                svec(1,jvec,ivec)*e1d(indvec)
                                                              yderv(isatom(ivec)) = yderv(isatom(ivec)) +  &
                                                                svec(2,jvec,ivec)*e1d(indvec)
                                                              zderv(isatom(ivec)) = zderv(isatom(ivec)) +  &
                                                                svec(3,jvec,ivec)*e1d(indvec)
                                                              xderv(isatom(jvec)) = xderv(isatom(jvec)) +  &
                                                                svec(1,ivec,jvec)*e1d(indvec)
                                                              yderv(isatom(jvec)) = yderv(isatom(jvec)) +  &
                                                                svec(2,ivec,jvec)*e1d(indvec)
                                                              zderv(isatom(jvec)) = zderv(isatom(jvec)) +  &
                                                                svec(3,ivec,jvec)*e1d(indvec)
                                                            enddo
                                                          enddo
                                                        endif
!
!  Loops / ifs for n
!
                                                      endif
                                                    enddo nloop
                                                  enddo
                                                enddo
                                              enddo
!
!  Loops / ifs for m
!
                                            enddo mloop
                                          enddo
                                        enddo
                                      enddo
!
!  Loops / ifs for l
!
                                    endif
                                  enddo lloop
                                enddo
                              enddo
                            enddo
!
!  Loops / ifs for k
!
                          enddo kloop
                        enddo
                      enddo
                    enddo
!
!  Loops / ifs for j
!
                  endif
                enddo jloop
              enddo
            enddo
          enddo
        enddo iloop
      endif
!
!  End of loop over cells for i
!
    enddo
!
!  End of outer loop over potentials
!
  enddo
!
!  Closing banner for energy decomposition
!
  if (lPrintSix) then
    call mpbarrier
    if (ioproc) then
      write(ioout,'(''--------------------------------------------------------------------------------'',/)')
    endif
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
  if (status/=0) call deallocate_error('sixmds','zderv')
  deallocate(yderv,stat=status)
  if (status/=0) call deallocate_error('sixmds','yderv')
  deallocate(xderv,stat=status)
  if (status/=0) call deallocate_error('sixmds','xderv')
!
!  Timing
!
  time2 = g_cpu_time()
  tsix = tsix + time2 - time1
#ifdef TRACE
  call trace_out('sixmds')
#endif
!
  return
  end
