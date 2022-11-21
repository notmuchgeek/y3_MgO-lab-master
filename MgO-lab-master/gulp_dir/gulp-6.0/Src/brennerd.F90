  subroutine brennerd(ebrenner,xkv,ykv,zkv,lgrad1,lgrad2,ldophonon)
!
!  Calculates the energy and derivatives for the Brenner potential.
!  It is assumed that all atoms are C or H and so this should be
!  checked during the setup phase.
!
!  Distributed memory parallel version.
!
!  On entry : 
!
!  lgrad1          = if .true. calculate the first derivatives
!  lgrad2          = if .true. calculate the second derivatives
!  ldophonon       = if .true. then this is a phonon run and so
!                    second derivatives should be phased
!  xkv             = if ldophonon this is the x kvector component
!  ykv             = if ldophonon this is the y kvector component
!  zkv             = if ldophonon this is the z kvector component
!
!  On exit :
!
!  ebrenner        = the value of the energy contribution
!
!   1/17 Created from brenner
!   4/17 Debug printing now only for ioproc
!   7/17 Sparsity added to increase speed for derivatives
!   2/18 Trace added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!  10/20 Site energies added for bond order potentials
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
  use datatypes
  use brennerdata
  use configurations, only : nregionno, nregions, lsliceatom, nregiontype, QMMMmode
  use control,        only : keyword, lseok
  use current
  use energies,       only : eattach, esregion12, esregion2
  use energies,       only : eregion2region, siteenergy
  use iochannels
  use neighbours
  use optimisation,   only : lfreeze, lopf
  use parallel
  use spatialbo,      only : lspatialok => lspatialBOok
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  real(dp),    intent(out)                         :: ebrenner
  real(dp),    intent(in)                          :: xkv
  real(dp),    intent(in)                          :: ykv
  real(dp),    intent(in)                          :: zkv
  logical,     intent(in)                          :: lgrad1
  logical,     intent(in)                          :: lgrad2
  logical,     intent(in)                          :: ldophonon
!
!  Local variables
!
  integer(i4)                                      :: i
  integer(i4)                                      :: ii
  integer(i4)                                      :: ind
  integer(i4)                                      :: indns
  integer(i4), dimension(:,:,:), allocatable, save :: ineigh
  integer(i4)                                      :: itmp
  integer(i4)                                      :: j
  integer(i4)                                      :: k
  integer(i4)                                      :: l
  integer(i4)                                      :: maxneigh2
  integer(i4)                                      :: maxneigh22
  integer(i4)                                      :: m
  integer(i4)                                      :: n
  integer(i4)                                      :: nati
  integer(i4)                                      :: natj
  integer(i4)                                      :: natk
  integer(i4)                                      :: ni
  integer(i4)                                      :: nj
  integer(i4)                                      :: nik
  integer(i4)                                      :: njk
  integer(i4)                                      :: nmin
  integer(i4)                                      :: nn
  integer(i4)                                      :: nn1
  integer(i4)                                      :: nn2
  integer(i4)                                      :: nptr
  integer(i4)                                      :: nri
  integer(i4)                                      :: nrj
  integer(i4)                                      :: nREBObo
  integer(i4)                                      :: nREBObo3
  integer(i4)                                      :: nREBOsi
  integer(i4)                                      :: nREBOsj
  integer(i4)                                      :: nREBOsk
  integer(i4)                                      :: nREBOatom
  integer(i4), dimension(:),     allocatable, save :: nREBOatomptr
  integer(i4), dimension(:),     allocatable, save :: nREBOatomRptr
  integer(i4), dimension(:,:),   allocatable, save :: nREBObond
  integer(i4), dimension(:,:),   allocatable, save :: neighno
  integer(i4), dimension(:),     allocatable, save :: nfreeatom
  integer(i4), dimension(:),     allocatable, save :: nneigh
  integer(i4)                                      :: nneighi2
  integer(i4)                                      :: nneighi2not
  integer(i4)                                      :: nneighj2
  integer(i4)                                      :: nneighj2not
  integer(i4)                                      :: nneighi22
  integer(i4)                                      :: nneighi22not
  integer(i4)                                      :: nneighj22
  integer(i4)                                      :: nneighj22not
  integer(i4)                                      :: nregioni
  integer(i4)                                      :: nregionj
  integer(i4)                                      :: nregiontypi
  integer(i4)                                      :: nregiontypj
  integer(i4)                                      :: ns1
  integer(i4)                                      :: ns2
  integer(i4)                                      :: status
  logical,     dimension(:),     allocatable, save :: latomdone
  logical                                          :: lattach
  logical                                          :: lfound
  logical                                          :: lilocal
  logical                                          :: lmaxneighok
  logical                                          :: lQMMMok
  logical                                          :: lreg2one
  logical                                          :: lreg2pair
  logical                                          :: lslicei
  logical                                          :: lslicej
  logical,     dimension(:),     allocatable, save :: lopanyneigh
  real(dp)                                         :: acoeff
  real(dp)                                         :: azeta
  real(dp)                                         :: bij
  real(dp)                                         :: bji
  real(dp)                                         :: bijsum
  real(dp)                                         :: bjisum
  real(dp)                                         :: btot
  real(dp)                                         :: btrmi
  real(dp)                                         :: btrmj
  real(dp)                                         :: g_cpu_time
  real(dp),    dimension(:),     allocatable, save :: d1i
  real(dp),    dimension(:),     allocatable, save :: d1j
  real(dp),    dimension(:),     allocatable, save :: d2i
  real(dp),    dimension(:),     allocatable, save :: d2j
  real(dp),    dimension(:),     allocatable, save :: d1Btoti
  real(dp),    dimension(:),     allocatable, save :: d1Btotj
  real(dp),    dimension(:),     allocatable, save :: d2Btoti
  real(dp),    dimension(:),     allocatable, save :: d2Btotj
  real(dp)                                         :: dfdr
  real(dp)                                         :: dfikdr
  real(dp)                                         :: dfjkdr
  real(dp)                                         :: d2fdr2
  real(dp)                                         :: d2fikdr2
  real(dp)                                         :: d2fjkdr2
  real(dp)                                         :: d3fdr3
  real(dp)                                         :: d3fikdr3
  real(dp)                                         :: d3fjkdr3
  real(dp)                                         :: dFxikdr
  real(dp)                                         :: dFxjkdr
  real(dp)                                         :: d2Fxikdr2
  real(dp)                                         :: d2Fxjkdr2
  real(dp)                                         :: d3Fxikdr3
  real(dp)                                         :: d3Fxjkdr3
  real(dp)                                         :: dGijkdr(3)
  real(dp)                                         :: dGjikdr(3)
  real(dp)                                         :: d2Gijkdr2(6)
  real(dp)                                         :: d2Gjikdr2(6)
  real(dp)                                         :: d2GijkdrdNti(3)
  real(dp)                                         :: d2GjikdrdNtj(3)
  real(dp)                                         :: d3Gijkdr3(10)
  real(dp)                                         :: d3Gjikdr3(10)
  real(dp)                                         :: dGijkdNti
  real(dp)                                         :: dGjikdNtj
  real(dp)                                         :: d2GijkdNti2
  real(dp)                                         :: d2GjikdNtj2
  real(dp)                                         :: eij
  real(dp)                                         :: expa
  real(dp)                                         :: expr
  real(dp)                                         :: expijk
  real(dp)                                         :: dexpijkdrij
  real(dp)                                         :: dexpijkdrik
  real(dp)                                         :: d2expijkdrij2
  real(dp)                                         :: d2expijkdrijdrik
  real(dp)                                         :: d2expijkdrik2
  real(dp)                                         :: expjik
  real(dp)                                         :: dexpjikdrji
  real(dp)                                         :: dexpjikdrjk
  real(dp)                                         :: d2expjikdrji2
  real(dp)                                         :: d2expjikdrjidrjk
  real(dp)                                         :: d2expjikdrjk2
  real(dp)                                         :: f
  real(dp)                                         :: Fij
  real(dp)                                         :: fik
  real(dp)                                         :: fjk
  real(dp)                                         :: Fxik
  real(dp)                                         :: Fxjk
  real(dp)                                         :: Gijk
  real(dp)                                         :: Gjik
  real(dp)                                         :: Pij
  real(dp)                                         :: Pji
  real(dp)                                         :: dPijdNi(nREBOspecies)
  real(dp)                                         :: dPjidNj(nREBOspecies)
  real(dp)                                         :: d2PijdN2i(nREBOspecies*(nREBOspecies+1)/2)
  real(dp)                                         :: d2PjidN2j(nREBOspecies*(nREBOspecies+1)/2)
  real(dp)                                         :: Nci
  real(dp)                                         :: Ncj
  real(dp)                                         :: Nhi
  real(dp)                                         :: Nhj
  real(dp),    dimension(:,:),   allocatable, save :: dNdr
  real(dp),    dimension(:,:),   allocatable, save :: d2Ndr2
  real(dp)                                         :: Nti
  real(dp)                                         :: Ntj
  real(dp)                                         :: Nconj
  real(dp)                                         :: Nconji
  real(dp)                                         :: Nconjj
  real(dp),    dimension(:),     allocatable, save :: dNconjdi
  real(dp),    dimension(:),     allocatable, save :: dNconjdj
  real(dp),    dimension(:),     allocatable, save :: d2Nconjdi2
  real(dp),    dimension(:),     allocatable, save :: d2Nconjdj2
  real(dp),    dimension(:,:),   allocatable, save :: nsneigh
  real(dp)                                         :: nsneighmfi(nREBOspecies)
  real(dp)                                         :: nsneighmfj(nREBOspecies)
  real(dp)                                         :: rcoeff
  real(dp)                                         :: rij
  real(dp)                                         :: rik
  real(dp)                                         :: rjk
  real(dp)                                         :: rrij
  real(dp)                                         :: rrik
  real(dp)                                         :: rrjk
  real(dp)                                         :: rtmp
  real(dp)                                         :: rzeta
  real(dp)                                         :: scale
  real(dp)                                         :: t1
  real(dp)                                         :: t2
  real(dp)                                         :: Tij
  real(dp),    dimension(:,:),   allocatable, save :: rneigh
  real(dp),    dimension(:,:),   allocatable, save :: xneigh
  real(dp),    dimension(:,:),   allocatable, save :: yneigh
  real(dp),    dimension(:,:),   allocatable, save :: zneigh
  real(dp)                                         :: Va
  real(dp)                                         :: Vr
  real(dp)                                         :: dVadr
  real(dp)                                         :: dVrdr
  real(dp)                                         :: d2Vadr2
  real(dp)                                         :: d2Vrdr2
  real(dp)                                         :: xdiff
  real(dp)                                         :: ydiff
  real(dp)                                         :: zdiff
  real(dp)                                         :: xik
  real(dp)                                         :: xjk
  real(dp)                                         :: xji
  real(dp)                                         :: yji
  real(dp)                                         :: zji
  real(dp)                                         :: xki
  real(dp)                                         :: yki
  real(dp)                                         :: zki
  real(dp)                                         :: xkj
  real(dp)                                         :: ykj
  real(dp)                                         :: zkj
#ifdef TRACE
  call trace_in('brennerd')
#endif
!
  t1 = g_cpu_time()
!
  allocate(nREBOatomptr(numat),stat=status)
  if (status/=0) call outofmemory('brennerd','nREBOatomptr')
  allocate(nREBOatomRptr(numat),stat=status)
  if (status/=0) call outofmemory('brennerd','nREBOatomRptr')
!
!  Set up a pointer to atoms that have a REBO type
!
  nREBOatom = 0
  do i = 1,numat
    if (nat2REBOspecies(nat(i)).gt.0) then
      nREBOatom = nREBOatom + 1
      nREBOatomptr(nREBOatom) = i
      nREBOatomRptr(i) = nREBOatom
    else
      nREBOatomRptr(i) = 0
    endif
  enddo
!
!  Allocate arrays that do not depend on maxneigh
!
  allocate(latomdone(nREBOatom),stat=status)
  if (status/=0) call outofmemory('brennerd','latomdone')
  allocate(lopanyneigh(nREBOatom),stat=status)
  if (status/=0) call outofmemory('brennerd','lopanyneigh')
  allocate(nfreeatom(numat),stat=status)
  if (status/=0) call outofmemory('brennerd','nfreeatom')
  allocate(nneigh(nREBOatom),stat=status)
  if (status/=0) call outofmemory('brennerd','nneigh')
  allocate(nsneigh(nREBOspecies,nREBOatom),stat=status)
  if (status/=0) call outofmemory('brennerd','nsneigh')
!
!  Initialise Brenner energy
!
  ebrenner = 0.0_dp
!
!  Reinitialisation point should maxneigh be increased
!
100 continue
  lmaxneighok = .true.
  if (allocated(rneigh)) then
    deallocate(d2Ndr2,stat=status)
    if (status/=0) call deallocate_error('brennerd','d2Ndr2')
    deallocate(d2Nconjdj2,stat=status)
    if (status/=0) call deallocate_error('brennerd','d2Nconjdj2')
    deallocate(d2Nconjdi2,stat=status)
    if (status/=0) call deallocate_error('brennerd','d2Nconjdi2')
    deallocate(d2Btotj,stat=status)
    if (status/=0) call deallocate_error('brennerd','d2Btotj')
    deallocate(d2Btoti,stat=status)
    if (status/=0) call deallocate_error('brennerd','d2Btoti')
    deallocate(d2j,stat=status)   
    if (status/=0) call deallocate_error('brennerd','d2j')
    deallocate(d2i,stat=status)   
    if (status/=0) call deallocate_error('brennerd','d2i')
!
    deallocate(dNdr,stat=status)
    if (status/=0) call deallocate_error('brennerd','dNdr')
    deallocate(dNconjdj,stat=status)
    if (status/=0) call deallocate_error('brennerd','dNconjdj')
    deallocate(dNconjdi,stat=status)
    if (status/=0) call deallocate_error('brennerd','dNconjdi')
    deallocate(d1Btotj,stat=status)
    if (status/=0) call deallocate_error('brennerd','d1Btotj')
    deallocate(d1Btoti,stat=status)
    if (status/=0) call deallocate_error('brennerd','d1Btoti')
    deallocate(d1j,stat=status)
    if (status/=0) call deallocate_error('brennerd','d1j')
    deallocate(d1i,stat=status)
    if (status/=0) call deallocate_error('brennerd','d1i')
!
    deallocate(zneigh,stat=status)
    if (status/=0) call deallocate_error('brennerd','zneigh')
    deallocate(yneigh,stat=status)
    if (status/=0) call deallocate_error('brennerd','yneigh')
    deallocate(xneigh,stat=status)
    if (status/=0) call deallocate_error('brennerd','xneigh')
    deallocate(rneigh,stat=status)
    if (status/=0) call deallocate_error('brennerd','rneigh')
    deallocate(ineigh,stat=status)
    if (status/=0) call deallocate_error('brennerd','ineigh')
    deallocate(neighno,stat=status)
    if (status/=0) call deallocate_error('brennerd','neighno')
    deallocate(nREBObond,stat=status)
    if (status/=0) call deallocate_error('brennerd','nREBObond')
  endif
!
!  Set parameter for pairwise storage memory
!
  maxneigh2 = maxneigh + maxneigh*(maxneigh + 1)/2 + maxneigh*(maxneigh+1)*maxneigh
  maxneigh22 = maxneigh2*(maxneigh2 + 1)/2
!
!  Allocate local memory
!
  allocate(nREBObond(maxneigh,nREBOatom),stat=status)
  if (status/=0) call outofmemory('brennerd','nREBObond')
  allocate(neighno(maxneigh,nREBOatom),stat=status)
  if (status/=0) call outofmemory('brennerd','neighno')
  allocate(ineigh(3_i4,maxneigh,nREBOatom),stat=status)
  if (status/=0) call outofmemory('brennerd','ineigh')
  allocate(rneigh(maxneigh,nREBOatom),stat=status)
  if (status/=0) call outofmemory('brennerd','rneigh')
  allocate(xneigh(maxneigh,nREBOatom),stat=status)
  if (status/=0) call outofmemory('brennerd','xneigh')
  allocate(yneigh(maxneigh,nREBOatom),stat=status)
  if (status/=0) call outofmemory('brennerd','yneigh')
  allocate(zneigh(maxneigh,nREBOatom),stat=status)
  if (status/=0) call outofmemory('brennerd','zneigh')
  if (lgrad1) then
    allocate(d1i(maxneigh2),stat=status)
    if (status/=0) call outofmemory('brennerd','d1i')
    allocate(d1j(maxneigh2),stat=status)
    if (status/=0) call outofmemory('brennerd','d1j')
    allocate(d1Btoti(maxneigh2),stat=status)
    if (status/=0) call outofmemory('brennerd','d1Btoti')
    allocate(d1Btotj(maxneigh2),stat=status)
    if (status/=0) call outofmemory('brennerd','d1Btotj')
    allocate(dNconjdi(maxneigh),stat=status)
    if (status/=0) call outofmemory('brennerd','dNconjdi')
    allocate(dNconjdj(maxneigh),stat=status)
    if (status/=0) call outofmemory('brennerd','dNconjdj')
    allocate(dNdr(maxneigh,nREBOatom),stat=status)
    if (status/=0) call outofmemory('brennerd','dNdr')
  else
    allocate(d1i(1),stat=status)
    if (status/=0) call outofmemory('brennerd','d1i')
    allocate(d1j(1),stat=status)
    if (status/=0) call outofmemory('brennerd','d1j')
    allocate(d1Btoti(1),stat=status)
    if (status/=0) call outofmemory('brennerd','d1Btoti')
    allocate(d1Btotj(1),stat=status)
    if (status/=0) call outofmemory('brennerd','d1Btotj')
    allocate(dNconjdi(1),stat=status)
    if (status/=0) call outofmemory('brennerd','dNconjdi')
    allocate(dNconjdj(1),stat=status)
    if (status/=0) call outofmemory('brennerd','dNconjdj')
    allocate(dNdr(1,1),stat=status)
    if (status/=0) call outofmemory('brennerd','dNdr')
  endif
  if (lgrad2) then
    allocate(d2i(maxneigh22),stat=status)   
    if (status/=0) call outofmemory('brennerd','d2i')
    allocate(d2j(maxneigh22),stat=status)   
    if (status/=0) call outofmemory('brennerd','d2j')
    allocate(d2Btoti(maxneigh22),stat=status)
    if (status/=0) call outofmemory('brennerd','d2Btoti')
    allocate(d2Btotj(maxneigh22),stat=status)
    if (status/=0) call outofmemory('brennerd','d2Btotj')
    allocate(d2Nconjdi2(maxneigh*(maxneigh+1)/2),stat=status)
    if (status/=0) call outofmemory('brennerd','d2Nconjdi2')
    allocate(d2Nconjdj2(maxneigh*(maxneigh+1)/2),stat=status)
    if (status/=0) call outofmemory('brennerd','d2Nconjdj2')
    allocate(d2Ndr2(maxneigh,nREBOatom),stat=status)
    if (status/=0) call outofmemory('brennerd','d2Ndr2')
  else
    allocate(d2i(1),stat=status)   
    if (status/=0) call outofmemory('brennerd','d2i')
    allocate(d2j(1),stat=status)   
    if (status/=0) call outofmemory('brennerd','d2j')
    allocate(d2Btoti(1),stat=status)
    if (status/=0) call outofmemory('brennerd','d2Btoti')
    allocate(d2Btotj(1),stat=status)
    if (status/=0) call outofmemory('brennerd','d2Btotj')
    allocate(d2Nconjdi2(1),stat=status)
    if (status/=0) call outofmemory('brennerd','d2Nconjdi2')
    allocate(d2Nconjdj2(1),stat=status)
    if (status/=0) call outofmemory('brennerd','d2Nconjdj2')
    allocate(d2Ndr2(1,1),stat=status)
    if (status/=0) call outofmemory('brennerd','d2Ndr2')
  endif
!****************************
!  Find list of free atoms  *
!****************************
  if (lfreeze.and..not.ldophonon) then
    ii = 0
    do i = 1,numat
      if (lopf(nrelf2a(i))) then
        ii = ii + 1
        nfreeatom(i) = ii
      else
        nfreeatom(i) = 0
      endif
    enddo
  else
    do i = 1,numat
      nfreeatom(i) = i
    enddo
  endif
!****************************************
!  Calculate neighbour lists for atoms  *
!****************************************
  call getREBOneighbour(maxneigh,bR22,nREBObond,nneigh,neighno,rneigh, &
                        xneigh,yneigh,zneigh,ineigh,latomdone,lmaxneighok, &
                        nREBOatom,nREBOatomptr,nREBOatomRptr)
!**********************************************************************
!  If maxneigh has been exceeded, increase value and return to start  *
!**********************************************************************
  if (.not.lmaxneighok) then
    do i = 1,nREBOatom 
      if (nneigh(i).gt.maxneigh) maxneigh = nneigh(i)
    enddo
    if (index(keyword,'verb').ne.0.and.ioproc) then
      write(ioout,'(/,''  Increasing maxneigh to '',i6)') maxneigh
    endif
    goto 100
  endif
!*******************************
!  Sort neighbours into order  *
!******************************* 
  if (lspatialok) then
    do i = 1,nREBOatom
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
          itmp = nREBObond(nptr,i)
          nREBObond(nptr,i) = nREBObond(nn,i)
          nREBObond(nn,i)  = itmp
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
  if (lgrad1) then
    dNdr(1:maxneigh,1:nREBOatom) = 0.0_dp
    if (lgrad2) then
      d2Ndr2(1:maxneigh,1:nREBOatom) = 0.0_dp
    endif
  endif
  do nri = 1,nREBOatom
    i = nREBOatomptr(nri)
    nsneigh(1:nREBOspecies,nri) = 0.0_dp
!
!  Set initial value for lopanyneigh - this
!  variable indicates whether an atom has 
!  any neighbours for which derivatives are
!  required
!
    if (.not.lfreeze.or.ldophonon) then
      lopanyneigh(nri) = .true.
    else
      lopanyneigh(nri) = lopf(nrelf2a(i))
    endif
    do n = 1,nneigh(nri)
      j = neighno(n,nri)
      nREBObo = nREBObond(n,nri)
      nREBOsj = nat2REBOspecies(nat(j))
      rij = rneigh(n,nri)
!
!  Check whether atom is free to optimise
!
      if (lopf(nrelf2a(j))) then
        lopanyneigh(nri) = .true.
      endif
!
!  Calculate function
!
      call ctaper(rij,bR1(nREBObo),bR2(nREBObo),f,dfdr,d2fdr2,d3fdr3,lgrad1,lgrad2,.false.)
      nsneigh(nREBOsj,nri) = nsneigh(nREBOsj,nri) + f
      if (lgrad1) then
        rrij = 1.0_dp/rij
        dNdr(n,nri) = dNdr(n,nri) + rrij*dfdr
        if (lgrad2) then
          d2Ndr2(n,nri) = d2Ndr2(n,nri) + rrij*rrij*(d2fdr2 - rrij*dfdr)
        endif
      endif
    enddo
  enddo
  if (index(keyword,'debu').ne.0.and.ioproc) then
    write(ioout,'(/,''  Number of neighbours for atoms :'',/)')
    write(ioout,'(''    i'',4x,''Nat'',3x,''N_C'',4x,''N_H'',4x,''N_O'',4x,''N_F'')')
    do nri = 1,nREBOatom
      i = nREBOatomptr(nri)
      write(ioout,'(i8,1x,i3,4(1x,f6.4))') i,nat(i),(nsneigh(j,nri),j=1,nREBOspecies)
    enddo
    write(ioout,'(/,''  Neighbours of atoms :'',/)')
    do nri = 1,nREBOatom
      i = nREBOatomptr(nri)
      write(ioout,'(i8,8(1x,i8))') i,(neighno(nn,nri),nn=1,nneigh(nri))
    enddo
  endif
!*****************************
!  Loop over pairs of atoms  *
!*****************************
  do nri = 1,nREBOatom
    i = nREBOatomptr(nri)
!
!  Set variables relating to i
!
    nregioni = nregionno(nsft + nrelf2a(i))
    nregiontypi = nregiontype(nregioni,ncf)
    lslicei = lsliceatom(nsft + nrelf2a(i))
    nati = nat(i)
    nREBOsi = nat2REBOspecies(nati)
!
    lilocal = (atom2local(i).gt.0)
!
!  Set total number of distances for neighbours of i
!
    nneighi2 = nneigh(nri) + nneigh(nri)*(nneigh(nri) + 1)/2 + nneigh(nri)*(maxneigh+1)*maxneigh
    nneighi22 = nneighi2*(nneighi2 + 1)/2
    nneighi2not = nneigh(nri) + nneigh(nri)*(nneigh(nri) + 1)/2
    nneighi22not = nneighi2not*(nneighi2not + 1)/2
!
!  Initialise derivative storage for neighbours of i
!
    if (lgrad1) then
      d1i(1:nneighi2) = 0.0_dp
      if (lgrad2) then
        d2i(1:nneighi22) = 0.0_dp
      endif
    endif
!
!  Loop over neighbours of i (=> j)
!
    ni = 1
    do while (ni.le.nneigh(nri).and.neighno(ni,nri).le.i)
!
      j = neighno(ni,nri)
      nrj = nREBOatomRptr(j)
!
!  Do we need to do this pair of atoms
!
      if (lopanyneigh(nri).or.lopanyneigh(nrj)) then
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
        nregionj = nregionno(nsft + nrelf2a(j))
        nregiontypj = nregiontype(nregionj,ncf)
        lslicej = lsliceatom(nsft + nrelf2a(j))
        natj = nat(j)
        nREBOsj = nat2REBOspecies(natj)
!
!  QM/MM handling : i & j are both QM atoms => exclude
!
        lQMMMok = .true.
        if (QMMMmode(ncf).gt.0) then
          if (nregiontypi.eq.1.and.nregiontypj.eq.1) lQMMMok = .false.
        endif
        if (lQMMMok) then
!
!  Set up i-j quantities
!
          nREBObo = nREBObond(ni,nri)
          rij = rneigh(ni,nri)
          xji = xneigh(ni,nri)
          yji = yneigh(ni,nri)
          zji = zneigh(ni,nri)
          rrij = 1.0_dp/rij
!
          lreg2one  = .false.
          lreg2pair = .false.
          if (lseok.and.nregions(ncf).gt.1) then
            lreg2pair = (nregioni.gt.1.and.nregionj.gt.1)
            if (.not.lreg2pair) lreg2one = (nregioni.gt.1.or.nregionj.gt.1)
          endif
          lattach = (lslicei.and..not.lslicej.or.lslicej.and..not.lslicei)
!
!  Find i in neighbour list for j
!
          nj = 1
          lfound = .false.
          do while (nj.le.nneigh(nrj).and..not.lfound)
            if (neighno(nj,nrj).eq.i) then
              xdiff = xneigh(nj,nrj) + xji
              ydiff = yneigh(nj,nrj) + yji
              zdiff = zneigh(nj,nrj) + zji
              lfound = ((abs(xdiff)+abs(ydiff)+abs(zdiff)).lt.1.0d-6)
            endif
            if (.not.lfound) nj = nj + 1
          enddo
!
!  Set total number of distances for neighbours of j
!
          nneighj2 = nneigh(nrj) + nneigh(nrj)*(nneigh(nrj) + 1)/2 + nneigh(nrj)*(maxneigh+1)*maxneigh
          nneighj22 = nneighj2*(nneighj2 + 1)/2
          nneighj2not = nneigh(nrj) + nneigh(nrj)*(nneigh(nrj) + 1)/2
          nneighj22not = nneighj2not*(nneighj2not + 1)/2
!
!  Initialise derivative storage for neighbours of i
!
          if (lgrad1) then
            d1j(1:nneighj2) = 0.0_dp
            if (lgrad2) then
              d2j(1:nneighj22) = 0.0_dp
            endif
          endif
!
!  Calculate fij
!
          call ctaper(rij,bR1(nREBObo),bR2(nREBObo),f,dfdr,d2fdr2,d3fdr3,lgrad1,lgrad2,.false.)
!
!  Calculate number of neighbours excluding i/j
!
          nsneighmfi(1:nREBOspecies) = nsneigh(1:nREBOspecies,nri)
          nsneighmfj(1:nREBOspecies) = nsneigh(1:nREBOspecies,nrj)
          nsneighmfi(nREBOsj) = nsneighmfi(nREBOsj) - f
          nsneighmfj(nREBOsi) = nsneighmfj(nREBOsi) - f
          Nci = nsneigh(1,nri)
          Nhi = nsneigh(2,nri)
          Ncj = nsneigh(1,nrj)
          Nhj = nsneigh(2,nrj)
          if (nREBOsi.eq.1) then
            Ncj = Ncj - f
          elseif (nREBOsi.eq.2) then
            Nhj = Nhj - f
          endif
          if (nREBOsj.eq.1) then
            Nci = Nci - f
          elseif (nREBOsj.eq.2) then
            Nhi = Nhi - f
          endif
          Nti = Nci + Nhi
          Ntj = Ncj + Nhj
!
!  Calculate Bij and Bji - loop over all other neighbours
!
          bijsum = 1.0_dp
          bjisum = 1.0_dp
          if (lgrad1) then
            d1Btoti(1:nneighi2) = 0.0_dp
            d1Btotj(1:nneighj2) = 0.0_dp
            if (lgrad2) then
              d2Btoti(1:nneighi22) = 0.0_dp
              d2Btotj(1:nneighj22) = 0.0_dp
            endif
          endif
!
!  Loop over neighbours of i .ne. j 
!
          do k = 1,nneigh(nri)
            natk = nat(neighno(k,nri))
            nREBOsk = nat2REBOspecies(natk)
            if (k.ne.ni) then
              rik = rneigh(k,nri)
              if (lgrad1) rrik = 1.0_dp/rik
              xki = xneigh(k,nri)
              yki = yneigh(k,nri)
              zki = zneigh(k,nri)
!
!  Set triad type indicator
!
              nREBObo3 = nREBOspecies2*(nREBOsi - 1) + nREBOspecies*(nREBOsj - 1) + nREBOsk
!
!  Calculate fik
!
              call ctaper(rik,bR1(nREBObond(k,nri)),bR2(nREBObond(k,nri)),fik,dfikdr,d2fikdr2, &
                          d3fikdr3,lgrad1,lgrad2,.false.)
!
!  Calculate Gijk
!
              call Gtheta(nREBOsi,Nti,xji,yji,zji,xki,yki,zki,Gijk,dGijkdr,d2Gijkdr2,d3Gijkdr3, &
                          dGijkdNti,d2GijkdNti2,d2GijkdrdNti,lgrad1,lgrad2,.false.)
!
!  Calculate exponential factor
!
              if (balpha(nREBObo3).ne.0.0_dp) then
                expijk = bexpco(nREBObo3)*exp(balpha(nREBObo3)*(rij - rik))
                if (lgrad1) then
                  dexpijkdrij = balpha(nREBObo3)*expijk*rrij
                  dexpijkdrik = - balpha(nREBObo3)*expijk*rrik
                  if (lgrad2) then
                    d2expijkdrij2 = rrij*dexpijkdrij*(balpha(nREBObo3) - rrij)
                    d2expijkdrijdrik = rrij*dexpijkdrik*balpha(nREBObo3)
                    d2expijkdrik2 = - rrik*dexpijkdrik*(balpha(nREBObo3) + rrik)
                  endif
                endif
              else
                expijk = bexpco(nREBObo3)
                if (lgrad1) then
                  dexpijkdrij = 0.0_dp
                  dexpijkdrik = 0.0_dp
                  if (lgrad2) then
                    d2expijkdrij2 = 0.0_dp
                    d2expijkdrijdrik = 0.0_dp
                    d2expijkdrik2 = 0.0_dp
                  endif
                endif
              endif
!
!  Combine terms
!
              bijsum = bijsum + Gijk*fik*expijk
              if (lgrad1) then
!
!  Derivatives
!
!  Find index for j-k 
!
                if (ni.ge.k) then
                  njk = nneigh(nri) + ni*(ni-1)/2 + k
                else
                  njk = nneigh(nri) + k*(k-1)/2 + ni
                endif
!
                dfikdr = rrik*dfikdr
!
                d1Btoti(k)   = d1Btoti(k)   + Gijk*dfikdr*expijk
!
                d1Btoti(ni)  = d1Btoti(ni)  + dGijkdr(1)*fik*expijk
                d1Btoti(k)   = d1Btoti(k)   + dGijkdr(2)*fik*expijk
                d1Btoti(njk) = d1Btoti(njk) + dGijkdr(3)*fik*expijk
!
                do l = 1,nneigh(nri)
                  if (l.ne.ni) then
                    d1Btoti(l) = d1Btoti(l) + dGijkdNti*dNdr(l,nri)*fik*expijk
                  endif
                enddo
!
                d1Btoti(ni) = d1Btoti(ni) + Gijk*fik*dexpijkdrij
                d1Btoti(k)  = d1Btoti(k)  + Gijk*fik*dexpijkdrik
!
                if (lgrad2) then
                  d2fikdr2 = rrik*rrik*(d2fikdr2 - dfikdr)
!
                  nn = ni*(ni + 1)/2
                  d2Btoti(nn) = d2Btoti(nn) + Gijk*fik*d2expijkdrij2
                  d2Btoti(nn) = d2Btoti(nn) + d2Gijkdr2(1)*fik*expijk
                  d2Btoti(nn) = d2Btoti(nn) + 2.0_dp*dGijkdr(1)*fik*dexpijkdrij
!
                  if (ni.ge.k) then
                    nn = ni*(ni - 1)/2 + k
                  else
                    nn = k*(k - 1)/2 + ni
                  endif
                  d2Btoti(nn) = d2Btoti(nn) + Gijk*fik*d2expijkdrijdrik
                  d2Btoti(nn) = d2Btoti(nn) + Gijk*dfikdr*dexpijkdrij
                  d2Btoti(nn) = d2Btoti(nn) + d2Gijkdr2(2)*fik*expijk
                  d2Btoti(nn) = d2Btoti(nn) + dGijkdr(1)*dfikdr*expijk
                  d2Btoti(nn) = d2Btoti(nn) + dGijkdr(1)*fik*dexpijkdrik
                  d2Btoti(nn) = d2Btoti(nn) + dGijkdr(2)*fik*dexpijkdrij
!
                  if (ni.ge.njk) then
                    nn = ni*(ni - 1)/2 + njk
                  else
                    nn = njk*(njk - 1)/2 + ni
                  endif
                  d2Btoti(nn) = d2Btoti(nn) + d2Gijkdr2(3)*fik*expijk
                  d2Btoti(nn) = d2Btoti(nn) + dGijkdr(3)*fik*dexpijkdrij
!
                  nn = k*(k + 1)/2
                  d2Btoti(nn) = d2Btoti(nn) + Gijk*d2fikdr2*expijk
                  d2Btoti(nn) = d2Btoti(nn) + 2.0_dp*Gijk*dfikdr*dexpijkdrik
                  d2Btoti(nn) = d2Btoti(nn) + Gijk*fik*d2expijkdrik2
                  d2Btoti(nn) = d2Btoti(nn) + d2Gijkdr2(4)*fik*expijk
                  d2Btoti(nn) = d2Btoti(nn) + 2.0_dp*dGijkdr(2)*dfikdr*expijk
                  d2Btoti(nn) = d2Btoti(nn) + 2.0_dp*dGijkdr(2)*fik*dexpijkdrik
!
                  if (k.ge.njk) then
                    nn = k*(k - 1)/2 + njk
                  else
                    nn = njk*(njk - 1)/2 + k
                  endif
                  d2Btoti(nn) = d2Btoti(nn) + d2Gijkdr2(5)*fik*expijk
                  d2Btoti(nn) = d2Btoti(nn) + dGijkdr(3)*dfikdr*expijk
                  d2Btoti(nn) = d2Btoti(nn) + dGijkdr(3)*fik*dexpijkdrik
!
                  nn = njk*(njk + 1)/2
                  d2Btoti(nn) = d2Btoti(nn) + d2Gijkdr2(6)*fik*expijk
!
                  do l = 1,nneigh(nri)
                    if (l.ne.ni) then
                      nn = l*(l + 1)/2
                      d2Btoti(nn) = d2Btoti(nn) + dGijkdNti*d2Ndr2(l,nri)*fik*expijk
!
                      if (k.ge.l) then
                        nn = k*(k - 1)/2 + l
                      else
                        nn = l*(l - 1)/2 + k
                      endif
                      d2Btoti(nn) = d2Btoti(nn) + dGijkdNti*dNdr(l,nri)*dfikdr*expijk
                      d2Btoti(nn) = d2Btoti(nn) + dGijkdNti*dNdr(l,nri)*fik*dexpijkdrik
                      d2Btoti(nn) = d2Btoti(nn) + 2.0_dp*d2GijkdrdNti(2)*dNdr(l,nri)*fik*expijk
!
                      if (ni.ge.l) then
                        nn = ni*(ni - 1)/2 + l
                      else
                        nn = l*(l - 1)/2 + ni
                      endif
                      d2Btoti(nn) = d2Btoti(nn) + dGijkdNti*dNdr(l,nri)*fik*dexpijkdrij
                      d2Btoti(nn) = d2Btoti(nn) + 2.0_dp*d2GijkdrdNti(1)*dNdr(l,nri)*fik*expijk
!
                      if (njk.ge.l) then
                        nn = njk*(njk - 1)/2 + l
                      else
                        nn = l*(l - 1)/2 + njk
                      endif
                      d2Btoti(nn) = d2Btoti(nn) + 2.0*d2GijkdrdNti(3)*dNdr(l,nri)*fik*expijk
!
                      do m = 1,l
                        if (m.ne.ni) then
                          nn = l*(l - 1)/2 + m
                          d2Btoti(nn) = d2Btoti(nn) + dGijkdNti*dNdr(l,nri)*dNdr(m,nri)*fik*expijk
                        endif
                      enddo
!
                    endif
                  enddo
                endif
              endif
            endif
          enddo
!
!  Loop over neighbours of j .ne. i 
!
          do k = 1,nneigh(nrj)
            natk = nat(neighno(k,nrj))
            nREBOsk = nat2REBOspecies(natk)
            if (k.ne.nj) then
              rjk = rneigh(k,nrj)
              if (lgrad1) rrjk = 1.0_dp/rjk
              xkj = xneigh(k,nrj)
              ykj = yneigh(k,nrj)
              zkj = zneigh(k,nrj)
!
!  Set triad type indicator
!           
              nREBObo3 = nREBOspecies2*(nREBOsj - 1) + nREBOspecies*(nREBOsi - 1) + nREBOsk
!
!  Calculate fik
!
              call ctaper(rjk,bR1(nREBObond(k,nrj)),bR2(nREBObond(k,nrj)),fjk,dfjkdr,d2fjkdr2, &
                          d3fjkdr3,lgrad1,lgrad2,.false.)
!
!  Calculate Gijk
!
              call Gtheta(nREBOsj,Ntj,-xji,-yji,-zji,xkj,ykj,zkj,Gjik,dGjikdr,d2Gjikdr2,d3Gjikdr3, &
                          dGjikdNtj,d2GjikdNtj2,d2GjikdrdNtj,lgrad1,lgrad2,.false.)
!
!  Calculate exponential factor
!
              if (balpha(nREBObo3).ne.0.0_dp) then
                expjik = bexpco(nREBObo3)*exp(balpha(nREBObo3)*(rij - rjk))
                if (lgrad1) then
                  dexpjikdrji = balpha(nREBObo3)*expjik*rrij
                  dexpjikdrjk = - balpha(nREBObo3)*expjik*rrjk
                  if (lgrad2) then
                    d2expjikdrji2 = rrij*dexpjikdrji*(balpha(nREBObo3) - rrij)
                    d2expjikdrjidrjk = rrij*dexpjikdrjk*balpha(nREBObo3)
                    d2expjikdrjk2 = - rrjk*dexpjikdrjk*(balpha(nREBObo3) + rrjk)
                  endif
                endif
              else
                expjik = bexpco(nREBObo3)
                if (lgrad1) then
                  dexpjikdrji = 0.0_dp
                  dexpjikdrjk = 0.0_dp
                  if (lgrad2) then
                    d2expjikdrji2 = 0.0_dp
                    d2expjikdrjidrjk = 0.0_dp
                    d2expjikdrjk2 = 0.0_dp
                  endif
                endif
              endif
!
!  Combine terms
!
              bjisum = bjisum + Gjik*fjk*expjik
              if (lgrad1) then
!
!  Derivatives
!
!  Find index for i-k
!
                if (nj.ge.k) then
                  nik = nneigh(nrj) + nj*(nj-1)/2 + k
                else
                  nik = nneigh(nrj) + k*(k-1)/2 + nj
                endif
!
                dfjkdr = rrjk*dfjkdr
!
                d1Btotj(k)   = d1Btotj(k)   + Gjik*dfjkdr*expjik
!
                d1Btotj(nj)  = d1Btotj(nj)  + dGjikdr(1)*fjk*expjik
                d1Btotj(k)   = d1Btotj(k)   + dGjikdr(2)*fjk*expjik
                d1Btotj(nik) = d1Btotj(nik) + dGjikdr(3)*fjk*expjik
!
                do l = 1,nneigh(nrj)
                  if (l.ne.nj) then
                    d1Btotj(l) = d1Btotj(l) + dGjikdNtj*dNdr(l,nrj)*fjk*expjik
                  endif
                enddo
!
                d1Btotj(nj) = d1Btotj(nj) + Gjik*fjk*dexpjikdrji
                d1Btotj(k)  = d1Btotj(k)  + Gjik*fjk*dexpjikdrjk
                if (lgrad2) then
                  d2fjkdr2 = rrjk*rrjk*(d2fjkdr2 - dfjkdr)
!
                  nn = nj*(nj + 1)/2
                  d2Btotj(nn) = d2Btotj(nn) + Gjik*fjk*d2expjikdrji2
                  d2Btotj(nn) = d2Btotj(nn) + d2Gjikdr2(1)*fjk*expjik
                  d2Btotj(nn) = d2Btotj(nn) + 2.0_dp*dGjikdr(1)*fjk*dexpjikdrji
!
                  if (nj.ge.k) then
                    nn = nj*(nj - 1)/2 + k
                  else
                    nn = k*(k - 1)/2 + nj
                  endif
                  d2Btotj(nn) = d2Btotj(nn) + Gjik*fjk*d2expjikdrjidrjk
                  d2Btotj(nn) = d2Btotj(nn) + Gjik*dfjkdr*dexpjikdrji
                  d2Btotj(nn) = d2Btotj(nn) + d2Gjikdr2(2)*fjk*expjik
                  d2Btotj(nn) = d2Btotj(nn) + dGjikdr(1)*dfjkdr*expjik
                  d2Btotj(nn) = d2Btotj(nn) + dGjikdr(1)*fjk*dexpjikdrjk
                  d2Btotj(nn) = d2Btotj(nn) + dGjikdr(2)*fjk*dexpjikdrji
!
                  if (nj.ge.nik) then
                    nn = nj*(nj - 1)/2 + nik
                  else
                    nn = nik*(nik - 1)/2 + nj
                  endif
                  d2Btotj(nn) = d2Btotj(nn) + d2Gjikdr2(3)*fjk*expjik
                  d2Btotj(nn) = d2Btotj(nn) + dGjikdr(3)*fjk*dexpjikdrji
!
                  nn = k*(k + 1)/2
                  d2Btotj(nn) = d2Btotj(nn) + Gjik*d2fjkdr2*expjik
                  d2Btotj(nn) = d2Btotj(nn) + 2.0_dp*Gjik*dfjkdr*dexpjikdrjk
                  d2Btotj(nn) = d2Btotj(nn) + Gjik*fjk*d2expjikdrjk2
                  d2Btotj(nn) = d2Btotj(nn) + d2Gjikdr2(4)*fjk*expjik
                  d2Btotj(nn) = d2Btotj(nn) + 2.0_dp*dGjikdr(2)*dfjkdr*expjik
                  d2Btotj(nn) = d2Btotj(nn) + 2.0_dp*dGjikdr(2)*fjk*dexpjikdrjk
!
                  if (k.ge.nik) then
                    nn = k*(k - 1)/2 + nik
                  else
                    nn = nik*(nik - 1)/2 + k
                  endif
                  d2Btotj(nn) = d2Btotj(nn) + d2Gjikdr2(5)*fjk*expjik
                  d2Btotj(nn) = d2Btotj(nn) + dGjikdr(3)*dfjkdr*expjik
                  d2Btotj(nn) = d2Btotj(nn) + dGjikdr(3)*fjk*dexpjikdrjk
!
                  nn = nik*(nik + 1)/2
                  d2Btotj(nn) = d2Btotj(nn) + d2Gjikdr2(6)*fjk*expjik
!
                  do l = 1,nneigh(nrj)
                    if (l.ne.nj) then
                      nn = l*(l + 1)/2
                      d2Btotj(nn) = d2Btotj(nn) + dGjikdNtj*d2Ndr2(l,nrj)*fjk*expjik
!
                      if (k.ge.l) then
                        nn = k*(k - 1)/2 + l
                      else
                        nn = l*(l - 1)/2 + k
                      endif
                      d2Btotj(nn) = d2Btotj(nn) + dGjikdNtj*dNdr(l,nrj)*dfjkdr*expjik
                      d2Btotj(nn) = d2Btotj(nn) + dGjikdNtj*dNdr(l,nrj)*fjk*dexpjikdrjk
                      d2Btotj(nn) = d2Btotj(nn) + 2.0_dp*d2GjikdrdNtj(2)*dNdr(l,nrj)*fjk*expjik
!
                      if (nj.ge.l) then
                        nn = nj*(nj - 1)/2 + l
                      else
                        nn = l*(l - 1)/2 + nj
                      endif
                      d2Btotj(nn) = d2Btotj(nn) + dGjikdNtj*dNdr(l,nrj)*fjk*dexpjikdrji
                      d2Btotj(nn) = d2Btotj(nn) + 2.0_dp*d2GjikdrdNtj(1)*dNdr(l,nrj)*fjk*expjik
!
                      if (nik.ge.l) then
                        nn = nik*(nik - 1)/2 + l
                      else
                        nn = l*(l - 1)/2 + nik
                      endif
                      d2Btotj(nn) = d2Btotj(nn) + 2.0_dp*d2GjikdrdNtj(3)*dNdr(l,nrj)*fjk*expjik
!
                      do m = 1,l
                        if (m.ne.nj) then
                          nn = l*(l - 1)/2 + m
                          d2Btotj(nn) = d2Btotj(nn) + dGjikdNtj*dNdr(l,nrj)*dNdr(m,nrj)*fjk*expjik
                        endif
                      enddo
!
                    endif
                  enddo
                endif
              endif
            endif
          enddo
!
!  Calculate and add Pij
!
          call calcP(nREBOsi,nsneighmfi,nREBObo,Pij,dPijdNi,d2PijdN2i,lgrad1,lgrad2)
          bijsum = bijsum + Pij
          if (lgrad1) then
            do nn = 1,nneigh(nri)
              if (nn.ne.ni) then
                ns1 = nat2REBOspecies(nat(neighno(nn,nri)))
                d1Btoti(nn) = d1Btoti(nn) + dPijdNi(ns1)*dNdr(nn,nri)
              endif
            enddo
            if (lgrad2) then
              do nn1 = 1,nneigh(nri)
                if (nn1.ne.ni) then
                  ns1 = nat2REBOspecies(nat(neighno(nn1,nri)))
                  do nn2 = 1,nn1
                    if (nn2.ne.ni) then
                      nn = nn1*(nn1 - 1)/2 + nn2
                      ns2 = nat2REBOspecies(nat(neighno(nn2,nri)))
                      if (ns1.gt.ns2) then
                        indns = ns1*(ns1 - 1)/2 + ns2
                      else
                        indns = ns2*(ns2 - 1)/2 + ns1
                      endif
                      d2Btoti(nn) = d2Btoti(nn) + d2PijdN2i(indns)*dNdr(nn1,nri)*dNdr(nn2,nri)
                    endif
                  enddo
!
                  nn = nn1*(nn1 + 1)/2
                  d2Btoti(nn) = d2Btoti(nn) + dPijdNi(ns1)*d2Ndr2(nn1,nri)
                endif
              enddo
            endif
          endif
!
!  Calculate and add Pji
!
          call calcP(nREBOsj,nsneighmfj,nREBObo,Pji,dPjidNj,d2PjidN2j,lgrad1,lgrad2)
          bjisum = bjisum + Pji
          if (lgrad1) then
            do nn = 1,nneigh(nrj)
              if (nn.ne.nj) then
                ns1 = nat2REBOspecies(nat(neighno(nn,nrj)))
                d1Btotj(nn) = d1Btotj(nn) + dPjidNj(ns1)*dNdr(nn,nrj)
              endif
            enddo
            if (lgrad2) then
              do nn1 = 1,nneigh(nrj)
                if (nn1.ne.nj) then
                  ns1 = nat2REBOspecies(nat(neighno(nn1,nrj)))
                  do nn2 = 1,nn1
                    if (nn2.ne.nj) then
                      nn = nn1*(nn1 - 1)/2 + nn2
                      ns2 = nat2REBOspecies(nat(neighno(nn2,nrj)))
                      if (ns1.gt.ns2) then
                        indns = ns1*(ns1 - 1)/2 + ns2
                      else
                        indns = ns2*(ns2 - 1)/2 + ns1
                      endif
                      d2Btotj(nn) = d2Btotj(nn) + d2PjidN2j(indns)*dNdr(nn1,nrj)*dNdr(nn2,nrj)
                    endif
                  enddo
!
                  nn = nn1*(nn1 + 1)/2
                  d2Btotj(nn) = d2Btotj(nn) + dPjidNj(ns1)*d2Ndr2(nn1,nrj)
                endif
              enddo
            endif
          endif
!
!  Raise terms to the power of -delta
!
          bij = bijsum**(-bdelta(nREBOsi))
          bji = bjisum**(-bdelta(nREBOsj))
!
!  Scale derivatives by bijsum/bjisum factors
!
          if (lgrad1) then
            if (lgrad2) then
!
!  Second derivatives
!
              btrmi = - 0.5_dp*bdelta(nREBOsi)*bij/bijsum
              do nn = 1,nneighi22not
                d2Btoti(nn) = btrmi*d2Btoti(nn)
              enddo
              btrmi = - btrmi*(bdelta(nREBOsi) + 1.0_dp)/bijsum
              nn = 0 
              do nn1 = 1,nneighi2not
                if (d1Btoti(nn1).ne.0.0_dp) then
                  do nn2 = 1,nn1
                    nn = nn + 1
                    d2Btoti(nn) = d2Btoti(nn) + btrmi*d1Btoti(nn2)*d1Btoti(nn1)
                  enddo
                else
                  nn = nn + nn1
                endif
              enddo
              btrmj = - 0.5_dp*bdelta(nREBOsj)*bji/bjisum
              do nn = 1,nneighj22not
                d2Btotj(nn) = btrmj*d2Btotj(nn)
              enddo
              btrmj = - btrmj*(bdelta(nREBOsj) + 1.0_dp)/bjisum
              nn = 0 
              do nn1 = 1,nneighj2not
                if (d1Btotj(nn1).ne.0.0_dp) then
                  do nn2 = 1,nn1
                    nn = nn + 1
                    d2Btotj(nn) = d2Btotj(nn) + btrmj*d1Btotj(nn2)*d1Btotj(nn1)
                  enddo
                else
                  nn = nn + nn1
                endif
              enddo
            endif
!
!  First derivatives
!
            do nn = 1,nneighi2not
              d1Btoti(nn) = -0.5_dp*bdelta(nREBOsi)*d1Btoti(nn)*bij/bijsum
            enddo
            do nn = 1,nneighj2not
              d1Btotj(nn) = -0.5_dp*bdelta(nREBOsj)*d1Btotj(nn)*bji/bjisum
            enddo
          endif
!
!  Calculate Nconj / Fij / Tij - only for C-C bonds
!
!  Loop over neighbouring carbon atoms of i and j
!
          Nconji = 0.0_dp
          Nconjj = 0.0_dp
          if (lgrad1) then
            dNconjdi(1:nneigh(nri)) = 0.0_dp
            dNconjdj(1:nneigh(nrj)) = 0.0_dp
            if (lgrad2) then
              d2Nconjdi2(1:nneigh(nri)*(nneigh(nri)+1)/2) = 0.0_dp
              d2Nconjdj2(1:nneigh(nrj)*(nneigh(nrj)+1)/2) = 0.0_dp
            endif
          endif
          do k = 1,nneigh(nri)
            if (k.ne.ni.and.nat(neighno(k,nri)).eq.6) then
              call ctaper(rneigh(k,nri),bR1(nREBObond(k,nri)),bR2(nREBObond(k,nri)),fik,dfikdr,d2fikdr2, &
                          d3fikdr3,lgrad1,lgrad2,.false.)
! Hard coded for C/H
              xik = nsneigh(1,nREBOatomRptr(neighno(k,nri))) + nsneigh(2,nREBOatomRptr(neighno(k,nri))) - fik
! Hard coded for C/H
              call ctaper(xik,2.0_dp,3.0_dp,Fxik,dFxikdr,d2Fxikdr2,d3Fxikdr3,lgrad1,lgrad2,.false.)
              Nconji = Nconji + fik*Fxik
              if (lgrad1) then
                rrik = 1.0_dp/rneigh(k,nri)
                dfikdr = rrik*dfikdr
                dNconjdi(k) = dNconjdi(k) + dfikdr*Fxik
                do nn = 1,nneigh(nri)
                  if (nn.ne.k) then
                    dNconjdi(nn) = dNconjdi(nn) + fik*dFxikdr*dNdr(nn,nri)
                  endif
                enddo
                if (lgrad2) then
                  d2fikdr2 = rrik*rrik*(d2fikdr2 - dfikdr)
                  nn = k*(k + 1)/2
                  d2Nconjdi2(nn) = d2Nconjdi2(nn) + d2fikdr2*Fxik
                  do nn1 = 1,nneigh(nri)
                    if (nn1.ne.k) then
                      nn = nn1*(nn1 + 1)/2
                      d2Nconjdi2(nn) = d2Nconjdi2(nn) + 2.0_dp*dfikdr*dFxikdr*dNdr(nn1,nri)
                      d2Nconjdi2(nn) = d2Nconjdi2(nn) + fik*dFxikdr*d2Ndr2(nn1,nri)
                      do nn2 = 1,nn1
                        if (nn2.ne.k) then
                          d2Nconjdi2(nn) = d2Nconjdi2(nn) + fik*d2Fxikdr2*dNdr(nn1,nri)*dNdr(nn2,nri)
                        endif
                      enddo
                    endif
                  enddo
                endif
              endif
            endif
          enddo
          do k = 1,nneigh(nrj)
            if (k.ne.nj.and.nat(neighno(k,nrj)).eq.6) then
              call ctaper(rneigh(k,nrj),bR1(nREBObond(k,nrj)),bR2(nREBObond(k,nrj)),fjk,dfjkdr,d2fjkdr2, &
                          d3fjkdr3,lgrad1,lgrad2,.false.)
! Hard coded for C/H
              xjk = nsneigh(1,nREBOatomRptr(neighno(k,nrj))) + nsneigh(2,nREBOatomRptr(neighno(k,nrj))) - fjk
! Hard coded for C/H
              call ctaper(xjk,2.0_dp,3.0_dp,Fxjk,dFxjkdr,d2Fxjkdr2,d3Fxjkdr3,lgrad1,lgrad2,.false.)
              Nconjj = Nconjj + fjk*Fxjk
              if (lgrad1) then
                rrjk = 1.0_dp/rneigh(k,nrj)
                dfjkdr = rrjk*dfjkdr
                dNconjdj(k) = dNconjdj(k) + dfjkdr*Fxjk
                do nn = 1,nneigh(nrj)
                  if (nn.ne.k) then
                    dNconjdj(nn) = dNconjdj(nn) + fjk*dFxjkdr*dNdr(nn,nrj)
                  endif
                enddo
                if (lgrad2) then
                  d2fjkdr2 = rrjk*rrjk*(d2fjkdr2 - dfjkdr)
                  nn = k*(k + 1)/2
                  d2Nconjdj2(nn) = d2Nconjdj2(nn) + d2fjkdr2*Fxjk
                  do nn1 = 1,nneigh(nrj)
                    if (nn1.ne.k) then
                      nn = nn1*(nn1 + 1)/2
                      d2Nconjdj2(nn) = d2Nconjdj2(nn) + 2.0_dp*dfjkdr*dFxjkdr*dNdr(nn1,nrj)
                      d2Nconjdj2(nn) = d2Nconjdj2(nn) + fjk*dFxjkdr*d2Ndr2(nn1,nrj)
                      do nn2 = 1,nn1
                        if (nn2.ne.k) then
                          d2Nconjdj2(nn) = d2Nconjdj2(nn) + fjk*d2Fxjkdr2*dNdr(nn1,nrj)*dNdr(nn2,nrj)
                        endif
                      enddo
                    endif
                  enddo
                endif
              endif
            endif
          enddo
          Nconj = 1.0_dp + Nconji**2 + Nconjj**2
!
          call calcF(nREBObo,nri,nrj,ni,nj,Nci+Nhi,Ncj+Nhj,Nconj,Nconji,Nconjj,maxneigh,nneigh,d1Btoti,d1Btotj, &
                     d2Btoti,d2Btotj,dNdr,d2Ndr2,dNconjdi,dNconjdj,d2Nconjdi2,d2Nconjdj2,Fij,lgrad1,lgrad2)
!
          call calcT(nri,nrj,ni,nj,Nci+Nhi,Ncj+Nhj,Nconj,Nconji,Nconjj,rij,xji,yji,zji,maxneigh,nneigh,rneigh, &
                     xneigh,yneigh,zneigh,nREBObond,Tij,d1Btoti,d1Btotj,d2Btoti,d2Btotj,dNdr,d2Ndr2,dNconjdi, &
                     dNconjdj,d2Nconjdi2,d2Nconjdj2,lgrad1,lgrad2)
!
!  Calculate two-body component of potential
!
          if (nbrennertype.eq.3) then
            expr = exp(-brepAlpha(nREBObo)*rij)
            Vr = brepA(nREBObo)*(1.0_dp + brepQ(nREBObo)*rrij)*expr
            Va = 0.0_dp
            do m = 1,3
              Va = Va + battB(m,nREBObo)*exp(-battBeta(m,nREBObo)*rij)
            enddo
            if (lgrad1) then
              dVrdr = - brepA(nREBObo)*brepQ(nREBObo)*(rrij**2)*expr - brepAlpha(nREBObo)*Vr
              dVadr = 0.0_dp
              do m = 1,3
                dVadr = dVadr - battB(m,nREBObo)*battBeta(m,nREBObo)*exp(-battBeta(m,nREBObo)*rij)
              enddo
              dVrdr = rrij*dVrdr
              dVadr = rrij*dVadr
              if (lgrad2) then
                d2Vrdr2 = 2.0_dp*brepA(nREBObo)*brepQ(nREBObo)*(rrij**2)*expr*(rrij + brepAlpha(nREBObo)) +  &
                          (brepAlpha(nREBObo)**2)*Vr
                d2Vadr2 = 0.0_dp
                do m = 1,3
                  d2Vadr2 = d2Vadr2 + battB(m,nREBObo)*(battBeta(m,nREBObo)**2)*exp(-battBeta(m,nREBObo)*rij)
                enddo
                d2Vrdr2 = rrij*rrij*(d2Vrdr2 - dVrdr)
                d2Vadr2 = rrij*rrij*(d2Vadr2 - dVadr)
              endif
            endif
          elseif (nbrennertype.eq.1) then
            rcoeff = brepA(nREBObo)/(brepQ(nREBObo) - 1.0_dp)
            acoeff = rcoeff*brepQ(nREBObo)
            rzeta = sqrt(2.0_dp*brepQ(nREBObo))*brepAlpha(nREBObo)
            azeta = sqrt(2.0_dp/brepQ(nREBObo))*brepAlpha(nREBObo)
            expr = exp(-rzeta*(rij - battB(1,nREBObo)))
            expa = exp(-azeta*(rij - battB(1,nREBObo)))
            Vr = rcoeff*expr
            Va = acoeff*expa
            if (lgrad1) then
              dVrdr = - rzeta*Vr
              dVadr = - azeta*Va
              dVrdr = rrij*dVrdr
              dVadr = rrij*dVadr
              if (lgrad2) then
                d2Vrdr2 = rzeta*rzeta*Vr
                d2Vadr2 = azeta*azeta*Va
                d2Vrdr2 = rrij*rrij*(d2Vrdr2 - dVrdr)
                d2Vadr2 = rrij*rrij*(d2Vadr2 - dVadr)
              endif
            endif
          endif
!
!  Calculate total i-j potential
!
          Btot = 0.5_dp*(bij + bji) + Fij + Tij
          eij = scale*f*(Vr - Btot*Va)
          if (lilocal) then
!
!  Add to surface energy totals if appropriate
!
            if (lseok) then
              if (lreg2one) then
                esregion12 = esregion12 + eij
              elseif (lreg2pair) then
                esregion2 = esregion2 + eij
              else
                ebrenner = ebrenner + eij
              endif
            else
              ebrenner = ebrenner + eij
            endif
            if (lattach) eattach = eattach + eij
!
            eregion2region(nregionj,nregioni) = eregion2region(nregionj,nregioni) + eij
!
            siteenergy(i) = siteenergy(i) + 0.5_dp*eij
            siteenergy(j) = siteenergy(j) + 0.5_dp*eij
          endif
!
!  Derivatives of Brenner potential energy
!
          if (lgrad1) then
            dfdr = rrij*dfdr
            d1i(ni) = d1i(ni) + scale*dfdr*(Vr - Btot*Va)
            d1i(ni) = d1i(ni) + scale*f*(dVrdr - Btot*dVadr)
            do nn = 1,nneighi2
              d1i(nn) = d1i(nn) - scale*f*Va*d1Btoti(nn)
            enddo
            do nn = 1,nneighj2
              d1j(nn) = d1j(nn) - scale*f*Va*d1Btotj(nn)
            enddo
            if (lgrad2) then
              ind = ni*(ni + 1)/2
              d2fdr2 = rrij*rrij*(d2fdr2 - dfdr)
              d2i(ind) = d2i(ind) + scale*d2fdr2*(Vr - Btot*Va)
              d2i(ind) = d2i(ind) + scale*2.0_dp*dfdr*(dVrdr - Btot*dVadr)
              d2i(ind) = d2i(ind) + scale*f*(d2Vrdr2 - Btot*d2Vadr2)
!
!  Use sparsity of first derivatives to accelerate second derivatives:
!
              nn = 0
              do nn1 = 1,nneighi2
                if (d1Btoti(nn1).ne.0.0_dp) then
                  if (ni.ge.nn1) then
                    ind = ni*(ni - 1)/2 + nn1
                  else
                    ind = nn1*(nn1 - 1)/2 + ni
                  endif
                  d2i(ind) = d2i(ind) - (dfdr*Va + f*dVadr)*d1Btoti(nn1)
                  if (ni.eq.nn1) then
                    d2i(ind) = d2i(ind) - (dfdr*Va + f*dVadr)*d1Btoti(nn1)
                  endif
                  do nn2 = 1,nn1
                    nn = nn + 1
                    d2i(nn) = d2i(nn) - f*Va*d2Btoti(nn)
                  enddo
                else
                  nn = nn + nn1
                endif
              enddo
              nn = 0
              do nn1 = 1,nneighj2not
                if (d1Btotj(nn1).ne.0.0_dp) then
                  if (nj.ge.nn1) then
                    ind = nj*(nj - 1)/2 + nn1
                  else
                    ind = nn1*(nn1 - 1)/2 + nj
                  endif
                  d2j(ind) = d2j(ind) - (dfdr*Va + f*dVadr)*d1Btotj(nn1)
                  if (nj.eq.nn1) then
                    d2j(ind) = d2j(ind) - (dfdr*Va + f*dVadr)*d1Btotj(nn1)
                  endif
                  do nn2 = 1,nn1
                    nn = nn + 1
                    d2j(nn) = d2j(nn) - f*Va*d2Btotj(nn)
                  enddo
                else
                  nn = nn + nn1
                endif
              enddo
            endif
          endif
!
!  Add derivatives due to neighbours of j 
!
          if (lgrad1) then
            if (lilocal) then
              call d1add(j,maxneigh,maxneigh,nneigh,neighno,xneigh,yneigh,zneigh,nfreeatom,nREBOatomRptr,d1j,.true.,.true.)
            endif
            if (lgrad2) then
              if (ldophonon) then
                call d2addpdm(j,maxneigh,nneigh,neighno,xneigh,yneigh,zneigh,nREBOatomRptr,d1j,d2j,xkv,ykv,zkv,.true.)
              else
                call d2adddm(j,maxneigh,nneigh,neighno,xneigh,yneigh,zneigh,nfreeatom,nREBOatomRptr,d1j,d2j,.true.)
              endif
            endif
          endif
!
!  End condition on QM/MM
!
        endif
!
!  End condition section on i or j being associated with moving atom
!
      endif
!
!  End of loop over neighbours of i
!
      ni = ni + 1
    enddo
!
!  Add derivatives due to neighbours of i
!
    if (lgrad1) then
      if (lilocal) then
        call d1add(i,maxneigh,maxneigh,nneigh,neighno,xneigh,yneigh,zneigh,nfreeatom,nREBOatomRptr,d1i,.true.,.true.)
      endif
      if (lgrad2) then
        if (ldophonon) then
          call d2addpdm(i,maxneigh,nneigh,neighno,xneigh,yneigh,zneigh,nREBOatomRptr,d1i,d2i,xkv,ykv,zkv,.true.)
        else
          call d2adddm(i,maxneigh,nneigh,neighno,xneigh,yneigh,zneigh,nfreeatom,nREBOatomRptr,d1i,d2i,.true.)
        endif
      endif
    endif
  enddo
!
!  Free local memory
!
  deallocate(d2Ndr2,stat=status)           
  if (status/=0) call deallocate_error('brennerd','d2Ndr2')
  deallocate(d2Nconjdj2,stat=status)             
  if (status/=0) call deallocate_error('brennerd','d2Nconjdj2')
  deallocate(d2Nconjdi2,stat=status)
  if (status/=0) call deallocate_error('brennerd','d2Nconjdi2')
  deallocate(d2Btotj,stat=status)              
  if (status/=0) call deallocate_error('brennerd','d2Btotj')
  deallocate(d2Btoti,stat=status)
  if (status/=0) call deallocate_error('brennerd','d2Btoti')
  deallocate(d2j,stat=status)
  if (status/=0) call deallocate_error('brennerd','d2j')
  deallocate(d2i,stat=status)
  if (status/=0) call deallocate_error('brennerd','d2i')
  deallocate(dNdr,stat=status)
  if (status/=0) call deallocate_error('brennerd','dNdr')
  deallocate(dNconjdj,stat=status)
  if (status/=0) call deallocate_error('brennerd','dNconjdj')
  deallocate(dNconjdi,stat=status)
  if (status/=0) call deallocate_error('brennerd','dNconjdi')
  deallocate(d1Btotj,stat=status)
  if (status/=0) call deallocate_error('brennerd','d1Btotj')
  deallocate(d1Btoti,stat=status)
  if (status/=0) call deallocate_error('brennerd','d1Btoti')
  deallocate(d1j,stat=status)
  if (status/=0) call deallocate_error('brennerd','d1j')
  deallocate(d1i,stat=status)
  if (status/=0) call deallocate_error('brennerd','d1i')
  deallocate(zneigh,stat=status)
  if (status/=0) call deallocate_error('brennerd','zneigh')
  deallocate(yneigh,stat=status)
  if (status/=0) call deallocate_error('brennerd','yneigh')
  deallocate(xneigh,stat=status)
  if (status/=0) call deallocate_error('brennerd','xneigh')
  deallocate(rneigh,stat=status)
  if (status/=0) call deallocate_error('brennerd','rneigh')
  deallocate(ineigh,stat=status)
  if (status/=0) call deallocate_error('brennerd','ineigh')
  deallocate(neighno,stat=status)
  if (status/=0) call deallocate_error('brennerd','neighno')
  deallocate(nREBObond,stat=status)
  if (status/=0) call deallocate_error('brennerd','nREBObond')
  deallocate(nsneigh,stat=status)
  if (status/=0) call deallocate_error('brennerd','nsneigh')
  deallocate(nneigh,stat=status)
  if (status/=0) call deallocate_error('brennerd','nneigh')
  deallocate(nfreeatom,stat=status)
  if (status/=0) call deallocate_error('brennerd','nfreeatom')
  deallocate(lopanyneigh,stat=status)
  if (status/=0) call deallocate_error('brennerd','lopanyneigh')
  deallocate(latomdone,stat=status)
  if (status/=0) call deallocate_error('brennerd','latomdone')
  deallocate(nREBOatomRptr,stat=status)
  if (status/=0) call deallocate_error('brennerd','nREBOatomRptr')
  deallocate(nREBOatomptr,stat=status)
  if (status/=0) call deallocate_error('brennerd','nREBOatomptr')
!
  t2 = g_cpu_time()
  tbrenner = tbrenner + t2 - t1
#ifdef TRACE
  call trace_out('brennerd')
#endif
!
  return
  end
