module m_brenner
  use datatypes
  integer(i4),                                save :: maxneigh2
  integer(i4),                                save :: nREBOatom
  integer(i4), dimension(:,:,:), allocatable, save :: ineigh
  integer(i4), dimension(:),     allocatable, save :: nREBOatomptr
  integer(i4), dimension(:),     allocatable, save :: nREBOatomRptr
  integer(i4), dimension(:,:),   allocatable, save :: nREBObond
  integer(i4), dimension(:,:),   allocatable, save :: neighno
  integer(i4), dimension(:),     allocatable, save :: nfreeatom
  integer(i4), dimension(:),     allocatable, save :: nneigh
  logical,     dimension(:),     allocatable, save :: lalreadydone
  logical,     dimension(:),     allocatable, save :: latomdone
  logical,     dimension(:),     allocatable, save :: lopanyneigh
  real(dp),    dimension(:,:),   allocatable, save :: nsneigh
  real(dp),    dimension(:,:),   allocatable, save :: rneigh
  real(dp),    dimension(:,:),   allocatable, save :: xneigh
  real(dp),    dimension(:,:),   allocatable, save :: yneigh
  real(dp),    dimension(:,:),   allocatable, save :: zneigh

contains

  subroutine setbrennerneigh
!
!  Calculates the neighbour list for the Brenner potential.
!
!  12/17 Created from brennermd.f
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
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
  use datatypes
  use brennerdata
  use control,        only : keyword
  use current
  use iochannels
  use neighbours
  use optimisation,   only : lfreeze, lopf
  use parallel
  use spatialbo,      only : lspatialok => lspatialBOok
  use times
  implicit none
!
!  Local variables
!
  integer(i4)                                      :: i
  integer(i4)                                      :: ii
  integer(i4)                                      :: itmp
  integer(i4)                                      :: j
  integer(i4)                                      :: k
  integer(i4)                                      :: l
  integer(i4)                                      :: n
  integer(i4)                                      :: n1
  integer(i4)                                      :: n2
  integer(i4)                                      :: nmin
  integer(i4)                                      :: nn
  integer(i4)                                      :: nn2
  integer(i4)                                      :: nptr
  integer(i4)                                      :: nri
  integer(i4)                                      :: nrj
  integer(i4)                                      :: nrk
  integer(i4)                                      :: nREBObo
  integer(i4)                                      :: nREBOsj
  integer(i4)                                      :: status
  logical                                          :: lmaxneighok
  real(dp)                                         :: g_cpu_time
  real(dp)                                         :: dfdr
  real(dp)                                         :: d2fdr2
  real(dp)                                         :: d3fdr3
  real(dp)                                         :: f
  real(dp)                                         :: rij
  real(dp)                                         :: rtmp
  real(dp)                                         :: t1
  real(dp)                                         :: t2
!
  t1 = g_cpu_time()
! 
  allocate(nREBOatomptr(numat),stat=status)
  if (status/=0) call outofmemory('setbrennerneigh','nREBOatomptr')
  allocate(nREBOatomRptr(numat),stat=status)
  if (status/=0) call outofmemory('setbrennerneigh','nREBOatomRptr')
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
!  Allocate local memory
!
  allocate(lalreadydone(nREBOatom),stat=status)
  if (status/=0) call outofmemory('setbrennerneigh','lalreadydone')
  allocate(latomdone(nREBOatom),stat=status)
  if (status/=0) call outofmemory('setbrennerneigh','latomdone')
  allocate(lopanyneigh(nREBOatom),stat=status)
  if (status/=0) call outofmemory('setbrennerneigh','lopanyneigh')
  allocate(nfreeatom(numat),stat=status)
  if (status/=0) call outofmemory('setbrennerneigh','nfreeatom')
  allocate(nneigh(nREBOatom),stat=status)
  if (status/=0) call outofmemory('setbrennerneigh','nneigh')
  allocate(nsneigh(nREBOspecies,nREBOatom),stat=status)
  if (status/=0) call outofmemory('setbrennerneigh','nsneigh')
!
!  Reinitialisation point should maxneigh be increased
!
100 continue
  lmaxneighok = .true.
  if (allocated(rneigh)) then
    deallocate(zneigh,stat=status)
    if (status/=0) call deallocate_error('setbrennerneigh','zneigh')
    deallocate(yneigh,stat=status)
    if (status/=0) call deallocate_error('setbrennerneigh','yneigh')
    deallocate(xneigh,stat=status)
    if (status/=0) call deallocate_error('setbrennerneigh','xneigh')
    deallocate(rneigh,stat=status)
    if (status/=0) call deallocate_error('setbrennerneigh','rneigh')
    deallocate(ineigh,stat=status)
    if (status/=0) call deallocate_error('setbrennerneigh','ineigh')
    deallocate(neighno,stat=status)
    if (status/=0) call deallocate_error('setbrennerneigh','neighno')
    deallocate(nREBObond,stat=status)
    if (status/=0) call deallocate_error('setbrennerneigh','nREBObond')
  endif
!
!  Set parameter for pairwise storage memory
!
  maxneigh2 = maxneigh + maxneigh*(maxneigh + 1)/2 + maxneigh*(maxneigh+1)*maxneigh
!
  allocate(nREBObond(maxneigh,nREBOatom),stat=status)
  if (status/=0) call outofmemory('setbrennerneigh','nREBObond')
  allocate(neighno(maxneigh,nREBOatom),stat=status)
  if (status/=0) call outofmemory('setbrennerneigh','neighno')
  allocate(ineigh(3_i4,maxneigh,nREBOatom),stat=status)
  if (status/=0) call outofmemory('setbrennerneigh','ineigh')
  allocate(rneigh(maxneigh,nREBOatom),stat=status)
  if (status/=0) call outofmemory('setbrennerneigh','rneigh')
  allocate(xneigh(maxneigh,nREBOatom),stat=status)
  if (status/=0) call outofmemory('setbrennerneigh','xneigh')
  allocate(yneigh(maxneigh,nREBOatom),stat=status)
  if (status/=0) call outofmemory('setbrennerneigh','yneigh')
  allocate(zneigh(maxneigh,nREBOatom),stat=status)
  if (status/=0) call outofmemory('setbrennerneigh','zneigh')
!****************************
!  Find list of free atoms  *
!****************************
  if (lfreeze) then
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
!
!  Set up logical array of atoms done, so that only those needed are done in parallel
!
  latomdone(1:nREBOatom) = .false.
!
!  Compute neighbour lists
!
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
    if (ioproc.and.index(keyword,'verb').ne.0) then
      write(ioout,'(/,''  Increasing maxneigh to '',i6)') maxneigh
    endif
    goto 100
  endif
!*******************************
!  Sort neighbours into order  *
!*******************************
  if (lspatialok) then
    do nri = 1,nREBOatom
      i = nREBOatomptr(nri)
      if (latomdone(nri)) then
!
!  Build pointer
!
        do nn = 1,nneigh(nri)
          nmin = numat + 1
          do nn2 = nn,nneigh(nri)
            if (neighno(nn2,nri).lt.nmin) then
              nmin = neighno(nn2,nri)
              nptr = nn2
            endif
          enddo
!
!  Sort quantities
!
          if (nptr.ne.nn) then
            itmp = neighno(nptr,nri)
            neighno(nptr,nri) = neighno(nn,nri)
            neighno(nn,nri)   = itmp
            itmp = nREBObond(nptr,nri)
            nREBObond(nptr,nri) = nREBObond(nn,nri)
            nREBObond(nn,nri)   = itmp
            itmp = ineigh(1,nptr,i)
            ineigh(1,nptr,i) = ineigh(1,nn,i)
            ineigh(1,nn,i)  = itmp
            itmp = ineigh(2,nptr,i)
            ineigh(2,nptr,i) = ineigh(2,nn,i)
            ineigh(2,nn,i)  = itmp
            itmp = ineigh(3,nptr,i)
            ineigh(3,nptr,i) = ineigh(3,nn,i)
            ineigh(3,nn,i)  = itmp
            rtmp = rneigh(nptr,nri)
            rneigh(nptr,nri) = rneigh(nn,nri)
            rneigh(nn,nri)   = rtmp
            rtmp = xneigh(nptr,nri)
            xneigh(nptr,nri) = xneigh(nn,nri)
            xneigh(nn,nri)   = rtmp
            rtmp = yneigh(nptr,nri)
            yneigh(nptr,nri) = yneigh(nn,nri)
            yneigh(nn,nri)   = rtmp
            rtmp = zneigh(nptr,nri)
            zneigh(nptr,nri) = zneigh(nn,nri)
            zneigh(nn,nri)   = rtmp
          endif
        enddo
      endif
    enddo
  endif
!*********************************************
!  Calculate numbers of neighbours for each  *
!*********************************************
  do nri = 1,nREBOatom
    i = nREBOatomptr(nri)
    if (latomdone(nri)) then
      nsneigh(1:nREBOspecies,nri) = 0.0_dp
!
!  Set initial value for lopanyneigh - this variable indicates whether an atom has 
!  any neighbours for which derivatives are required
!
      if (.not.lfreeze) then
        lopanyneigh(nri) = .true.
      else
        lopanyneigh(nri) = lopf(nrelf2a(i))
      endif
      do n = 1,nneigh(nri)
        j = neighno(n,nri)
        nrj = nREBOatomRptr(j)
        nREBObo = nREBObond(n,nri)
        nREBOsj = nat2REBOspecies(nat(j))
        rij = rneigh(n,nri)
!
!  Check whether energy of atom will be affected by those atoms free to optimise
!
        if (lopf(nrelf2a(j))) then
          lopanyneigh(nri) = .true.
        else
          do n1 = 1,nneigh(nrj)
            k = neighno(n1,nrj)
            nrk = nREBOatomRptr(k)
            if (lopf(nrelf2a(k))) then
              lopanyneigh(nri) = .true.
            else
              do n2 = 1,nneigh(nrk)
                l = neighno(n2,nrk)
                if (lopf(nrelf2a(l))) then
                  lopanyneigh(nri) = .true.
                endif
              enddo
            endif
          enddo
        endif
!
!  Calculate function
!
        call ctaper(rij,bR1(nREBObo),bR2(nREBObo),f,dfdr,d2fdr2,d3fdr3,.false.,.false.,.false.)
        nsneigh(nREBOsj,nri) = nsneigh(nREBOsj,nri) + f
      enddo
    endif
  enddo
  if (ioproc.and.index(keyword,'debu').ne.0) then
    write(ioout,'(/,''  Number of neighbours for atoms :'',/)')
    write(ioout,'(''    i'',4x,''Nat'',3x,''N_C'',4x,''N_H'',4x,''N_O'',4x,''N_F'')')
    do i = 1,numat
      nri = nREBOatomRptr(i)
      if (nri.gt.0) then
        write(ioout,'(i8,1x,i3,4(1x,f6.4))') i,nat(i),(nsneigh(j,nri),j=1,nREBOspecies)
      endif
    enddo
    write(ioout,'(/,''  Neighbours of atoms :'',/)')
    do i = 1,numat
      nri = nREBOatomRptr(i)
      if (nri.gt.0) then
        write(ioout,'(i8,8(1x,i8))') i,(neighno(nn,nri),nn=1,nneigh(nri))
      endif
    enddo
  endif
!
  t2 = g_cpu_time()
  tbrenner = tbrenner + t2 - t1
!
  return
  end subroutine setbrennerneigh

  subroutine brennerfd(matom,mcrd,step,ebrenner,lgrad1)
!
!  Calculates the energy and derivatives for the Brenner potential.
!  It is assumed that all atoms are C or H and so this should be
!  checked during the setup phase. This routine only does up to first
!  derivatives.
!
!  On entry : 
!
!  lgrad1          = if .true. calculate the first derivatives
!  matom           = pointer to atom whose derivatives are being computed
!  mcrd            = pointer to coordinate that has been changed
!  step            = change in coordinate in Angstroms
!
!  On exit :
!
!  ebrenner        = the value of the energy contribution
!
!  12/17 Created from brennermd.f90
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
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
  use datatypes
  use brennerdata
  use configurations, only : nregionno, nregiontype, QMMMmode
  use current
  use neighbours
  use times
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: matom
  integer(i4), intent(in)                          :: mcrd
  real(dp),    intent(in)                          :: step
  real(dp),    intent(inout)                       :: ebrenner
  logical,     intent(in)                          :: lgrad1
!
!  Local variables
!
  integer(i4)                                      :: i
  integer(i4)                                      :: j
  integer(i4)                                      :: k
  integer(i4)                                      :: l
  integer(i4)                                      :: m
  integer(i4)                                      :: n
  integer(i4)                                      :: natk
  integer(i4)                                      :: ni
  integer(i4)                                      :: nj
  integer(i4)                                      :: nik
  integer(i4)                                      :: njk
  integer(i4)                                      :: nn
  integer(i4)                                      :: nri
  integer(i4)                                      :: nrj
  integer(i4)                                      :: nREBObo
  integer(i4)                                      :: nREBObo3
  integer(i4)                                      :: nREBOsi
  integer(i4)                                      :: nREBOsj
  integer(i4)                                      :: nREBOsk
  integer(i4)                                      :: nneighi2
  integer(i4)                                      :: nneighi2not
  integer(i4)                                      :: nneighj2
  integer(i4)                                      :: nneighj2not
  integer(i4)                                      :: nregioni
  integer(i4)                                      :: nregionj
  integer(i4)                                      :: nregiontypi
  integer(i4)                                      :: nregiontypj
  integer(i4)                                      :: ns1
  integer(i4)                                      :: status
  logical                                          :: lfound
  logical                                          :: lQMMMok
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
  real(dp)                                         :: expjik
  real(dp)                                         :: dexpjikdrji
  real(dp)                                         :: dexpjikdrjk
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
  real(dp)                                         :: d2PijdN2i(1)
  real(dp)                                         :: d2PjidN2j(1)
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
  real(dp)                                         :: Va
  real(dp)                                         :: Vr
  real(dp)                                         :: dVadr
  real(dp)                                         :: dVrdr
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
!
  t1 = g_cpu_time()
!
!  Allocate local memory
!
  allocate(d2Btoti(1),stat=status)
  if (status/=0) call outofmemory('brennerfd','d2Btoti')
  allocate(d2Btotj(1),stat=status)
  if (status/=0) call outofmemory('brennerfd','d2Btotj')
  allocate(d2Nconjdi2(1),stat=status)
  if (status/=0) call outofmemory('brennerfd','d2Nconjdi2')
  allocate(d2Nconjdj2(1),stat=status)
  if (status/=0) call outofmemory('brennerfd','d2Nconjdj2')
  allocate(d2Ndr2(1,1),stat=status)
  if (status/=0) call outofmemory('brennerfd','d2Ndr2')
!
  if (lgrad1) then
    allocate(d1i(maxneigh2),stat=status)
    if (status/=0) call outofmemory('brennerfd','d1i')
    allocate(d1j(maxneigh2),stat=status)
    if (status/=0) call outofmemory('brennerfd','d1j')
    allocate(d1Btoti(maxneigh2),stat=status)
    if (status/=0) call outofmemory('brennerfd','d1Btoti')
    allocate(d1Btotj(maxneigh2),stat=status)
    if (status/=0) call outofmemory('brennerfd','d1Btotj')
    allocate(dNconjdi(maxneigh),stat=status)
    if (status/=0) call outofmemory('brennerfd','dNconjdi')
    allocate(dNconjdj(maxneigh),stat=status)
    if (status/=0) call outofmemory('brennerfd','dNconjdj')
    allocate(dNdr(maxneigh,numat),stat=status)
    if (status/=0) call outofmemory('brennerfd','dNdr')
  else
    allocate(d1i(1),stat=status)
    if (status/=0) call outofmemory('brennerfd','d1i')
    allocate(d1j(1),stat=status)
    if (status/=0) call outofmemory('brennerfd','d1j')
    allocate(d1Btoti(1),stat=status)
    if (status/=0) call outofmemory('brennerfd','d1Btoti')
    allocate(d1Btotj(1),stat=status)
    if (status/=0) call outofmemory('brennerfd','d1Btotj')
    allocate(dNconjdi(1),stat=status)
    if (status/=0) call outofmemory('brennerfd','dNconjdi')
    allocate(dNconjdj(1),stat=status)
    if (status/=0) call outofmemory('brennerfd','dNconjdj')
    allocate(dNdr(1,1),stat=status)
    if (status/=0) call outofmemory('brennerfd','dNdr')
  endif
!*****************************************
!  Add shift to vectors involving matom  *
!*****************************************
  do i = 1,numat
    if (i.eq.matom) then
      do nn = 1,nneigh(i)
        if (neighno(nn,i).ne.matom) then
          if (mcrd.eq.1) then
            xneigh(nn,i) = xneigh(nn,i) - step
          elseif (mcrd.eq.2) then
            yneigh(nn,i) = yneigh(nn,i) - step
          elseif (mcrd.eq.3) then
            zneigh(nn,i) = zneigh(nn,i) - step
          endif
          rtmp = xneigh(nn,i)**2 + yneigh(nn,i)**2 + zneigh(nn,i)**2
          rneigh(nn,i) = sqrt(rtmp)
        endif
      enddo
    else
      do nn = 1,nneigh(i)
        j = neighno(nn,i)
        if (j.eq.matom) then
          if (mcrd.eq.1) then
            xneigh(nn,i) = xneigh(nn,i) + step
          elseif (mcrd.eq.2) then
            yneigh(nn,i) = yneigh(nn,i) + step
          elseif (mcrd.eq.3) then
            zneigh(nn,i) = zneigh(nn,i) + step
          endif
          rtmp = xneigh(nn,i)**2 + yneigh(nn,i)**2 + zneigh(nn,i)**2
          rneigh(nn,i) = sqrt(rtmp)
        endif
      enddo
    endif
  enddo
!*********************************************
!  Calculate numbers of neighbours for each  *
!*********************************************
  if (lgrad1) then
    dNdr(1:maxneigh,1:nREBOatom) = 0.0_dp
  endif
  do nri = 1,nREBOatom
    i = nREBOatomptr(nri)
    if (latomdone(nri)) then
      nsneigh(1:nREBOspecies,nri) = 0.0_dp
      do n = 1,nneigh(nri)
        j = neighno(n,nri)
        nrj = nREBOatomRptr(j)
        nREBObo = nREBObond(n,nri)
        nREBOsj = nat2REBOspecies(nat(j))
        rij = rneigh(n,nri)
!
!  Calculate function
!
        call ctaper(rij,bR1(nREBObo),bR2(nREBObo),f,dfdr,d2fdr2,d3fdr3,lgrad1,.false.,.false.)
        nsneigh(nREBOsj,nri) = nsneigh(nREBOsj,nri) + f
        if (lgrad1) then
          rrij = 1.0_dp/rij
          dNdr(n,nri) = dNdr(n,nri) + rrij*dfdr
        endif
      enddo
    endif
  enddo
!*****************************
!  Loop over pairs of atoms  *
!*****************************
  icloop: do i = 1,numat
!
!  Check whether atom is a REBO species and if not skip
!
    nri = nREBOatomRptr(i)
    if (nri.eq.0) cycle icloop
!orig do i = 1,numat
!
!  Set variables relating to i
!
    nregioni = nregionno(nsft + nrelf2a(i))
    nregiontypi = nregiontype(nregioni,ncf)
    nREBOsi = nat2REBOspecies(nat(i))
!
!  Set total number of distances for neighbours of i
!
    nneighi2 = nneigh(nri) + nneigh(nri)*(nneigh(nri) + 1)/2 + nneigh(nri)*(maxneigh+1)*maxneigh
    nneighi2not = nneigh(nri) + nneigh(nri)*(nneigh(nri) + 1)/2
!
!  Initialise derivative storage for neighbours of i
!
    if (lgrad1) then
      d1i(1:nneighi2) = 0.0_dp
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
        nREBOsj = nat2REBOspecies(nat(j))
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
          nneighj2not = nneigh(nrj) + nneigh(nrj)*(nneigh(nrj) + 1)/2
!
!  Initialise derivative storage for neighbours of i
!
          if (lgrad1) then
            d1j(1:nneighj2) = 0.0_dp
          endif
!
!  Calculate fij
!
          call ctaper(rij,bR1(nREBObo),bR2(nREBObo),f,dfdr,d2fdr2,d3fdr3,lgrad1,.false.,.false.)
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
                          d3fikdr3,lgrad1,.false.,.false.)
!
!  Calculate Gijk
!
              call Gtheta(nREBOsi,Nti,xji,yji,zji,xki,yki,zki,Gijk,dGijkdr,d2Gijkdr2,d3Gijkdr3, &
                          dGijkdNti,d2GijkdNti2,d2GijkdrdNti,lgrad1,.false.,.false.)
!
!  Calculate exponential factor
!
              if (balpha(nREBObo3).ne.0.0_dp) then
                expijk = bexpco(nREBObo3)*exp(balpha(nREBObo3)*(rij - rik))
                if (lgrad1) then
                  dexpijkdrij = balpha(nREBObo3)*expijk*rrij
                  dexpijkdrik = - balpha(nREBObo3)*expijk*rrik
                endif
              else
                expijk = bexpco(nREBObo3)
                if (lgrad1) then
                  dexpijkdrij = 0.0_dp
                  dexpijkdrik = 0.0_dp
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
                          d3fjkdr3,lgrad1,.false.,.false.)
!
!  Calculate Gijk
!
              call Gtheta(nREBOsj,Ntj,-xji,-yji,-zji,xkj,ykj,zkj,Gjik,dGjikdr,d2Gjikdr2, &
                          d3Gjikdr3,dGjikdNtj,d2GjikdNtj2,d2GjikdrdNtj,lgrad1,.false.,.false.)
!
!  Calculate exponential factor
!
              if (balpha(nREBObo3).ne.0.0_dp) then
                expjik = bexpco(nREBObo3)*exp(balpha(nREBObo3)*(rij - rjk))
                if (lgrad1) then
                  dexpjikdrji = balpha(nREBObo3)*expjik*rrij
                  dexpjikdrjk = - balpha(nREBObo3)*expjik*rrjk
                endif
              else
                expjik = bexpco(nREBObo3)
                if (lgrad1) then
                  dexpjikdrji = 0.0_dp
                  dexpjikdrjk = 0.0_dp
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
              endif
            endif
          enddo
!
!  Calculate and add Pij
!
          call calcP(nREBOsi,nsneighmfi,nREBObo,Pij,dPijdNi,d2PijdN2i,lgrad1,.false.)
          bijsum = bijsum + Pij
          if (lgrad1) then
            do nn = 1,nneigh(nri)
              if (nn.ne.ni) then
                ns1 = nat2REBOspecies(nat(neighno(nn,nri)))
                d1Btoti(nn) = d1Btoti(nn) + dPijdNi(ns1)*dNdr(nn,nri)
              endif
            enddo
          endif
!
!  Calculate and add Pji
!
          call calcP(nREBOsj,nsneighmfj,nREBObo,Pji,dPjidNj,d2PjidN2j,lgrad1,.false.)
          bjisum = bjisum + Pji
          if (lgrad1) then
            do nn = 1,nneigh(nrj)
              if (nn.ne.nj) then
                ns1 = nat2REBOspecies(nat(neighno(nn,nrj)))
                d1Btotj(nn) = d1Btotj(nn) + dPjidNj(ns1)*dNdr(nn,nrj)
              endif
            enddo
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
!
!  First derivatives
!
            btrmi = - 0.5_dp*bdelta(nREBOsi)*bij/bijsum
            do nn = 1,nneighi2not
              d1Btoti(nn) = btrmi*d1Btoti(nn)
            enddo
            btrmj = - 0.5_dp*bdelta(nREBOsj)*bji/bjisum
            do nn = 1,nneighj2not
              d1Btotj(nn) = btrmj*d1Btotj(nn)
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
          endif
          do k = 1,nneigh(nri)
            if (k.ne.ni.and.nat(neighno(k,nri)).eq.6) then
              call ctaper(rneigh(k,nri),bR1(nREBObond(k,nri)),bR2(nREBObond(k,nri)),fik,dfikdr,d2fikdr2, &
                          d3fikdr3,lgrad1,.false.,.false.)
! Hard coded for C/H
              xik = nsneigh(1,nREBOatomRptr(neighno(k,nri))) + nsneigh(2,nREBOatomRptr(neighno(k,nri))) - fik
! Hard coded for C/H
              call ctaper(xik,2.0_dp,3.0_dp,Fxik,dFxikdr,d2Fxikdr2,d3Fxikdr3,lgrad1,.false.,.false.)
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
              endif
            endif
          enddo
          do k = 1,nneigh(nrj)
            if (k.ne.nj.and.nat(neighno(k,nrj)).eq.6) then
              call ctaper(rneigh(k,nrj),bR1(nREBObond(k,nrj)),bR2(nREBObond(k,nrj)),fjk,dfjkdr,d2fjkdr2, &
                          d3fjkdr3,lgrad1,.false.,.false.)
! Hard coded for C/H
              xjk = nsneigh(1,nREBOatomRptr(neighno(k,nrj))) + nsneigh(2,nREBOatomRptr(neighno(k,nrj))) - fjk
! Hard coded for C/H
              call ctaper(xjk,2.0_dp,3.0_dp,Fxjk,dFxjkdr,d2Fxjkdr2,d3Fxjkdr3,lgrad1,.false.,.false.)
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
              endif
            endif
          enddo
          Nconj = 1.0_dp + Nconji**2 + Nconjj**2
!
          call calcF(nREBObo,nri,nrj,ni,nj,Nci+Nhi,Ncj+Nhj,Nconj,Nconji,Nconjj,maxneigh,nneigh,d1Btoti, &
                     d1Btotj,d2Btoti,d2Btotj,dNdr,d2Ndr2,dNconjdi,dNconjdj,d2Nconjdi2,d2Nconjdj2, &
                     Fij,lgrad1,.false.)
!
          call calcT(nri,nrj,ni,nj,Nci+Nhi,Ncj+Nhj,Nconj,Nconji,Nconjj,rij,xji,yji,zji,maxneigh,nneigh, &
                     rneigh,xneigh,yneigh,zneigh,nREBObond,Tij,d1Btoti,d1Btotj,d2Btoti,d2Btotj,dNdr, &
                     d2Ndr2,dNconjdi,dNconjdj,d2Nconjdi2,d2Nconjdj2,lgrad1,.false.)
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
            endif
          endif
!
!  Calculate total i-j potential
!
          Btot = 0.5_dp*(bij + bji) + Fij + Tij
          eij = scale*f*(Vr - Btot*Va)
!
          ebrenner = ebrenner + eij
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
            do nn = 1,nneighj2not
              d1j(nn) = d1j(nn) - scale*f*Va*d1Btotj(nn)
            enddo
          endif
!
!  Add derivatives due to neighbours of j 
!
          if (lgrad1) then
            call d1add(j,maxneigh,maxneigh,nneigh,neighno,xneigh,yneigh,zneigh,nfreeatom,nREBOatomRptr,d1j,.true.,.true.)
          endif
!
!  End of QM/MM scheme condition
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
      call d1add(i,maxneigh,maxneigh,nneigh,neighno,xneigh,yneigh,zneigh,nfreeatom,nREBOatomRptr,d1i,.true.,.true.)
    endif
  enddo icloop
!**********************************************
!  Subtract shift to vectors involving matom  *
!**********************************************
  do i = 1,numat
    if (i.eq.matom) then
      do nn = 1,nneigh(i)
        if (neighno(nn,i).ne.matom) then
          if (mcrd.eq.1) then
            xneigh(nn,i) = xneigh(nn,i) + step
          elseif (mcrd.eq.2) then
            yneigh(nn,i) = yneigh(nn,i) + step
          elseif (mcrd.eq.3) then
            zneigh(nn,i) = zneigh(nn,i) + step
          endif
          rtmp = xneigh(nn,i)**2 + yneigh(nn,i)**2 + zneigh(nn,i)**2
          rneigh(nn,i) = sqrt(rtmp)
        endif
      enddo
    else
      do nn = 1,nneigh(i)
        j = neighno(nn,i)
        if (j.eq.matom) then
          if (mcrd.eq.1) then
            xneigh(nn,i) = xneigh(nn,i) - step
          elseif (mcrd.eq.2) then
            yneigh(nn,i) = yneigh(nn,i) - step
          elseif (mcrd.eq.3) then
            zneigh(nn,i) = zneigh(nn,i) - step
          endif
          rtmp = xneigh(nn,i)**2 + yneigh(nn,i)**2 + zneigh(nn,i)**2
          rneigh(nn,i) = sqrt(rtmp)
        endif
      enddo
    endif
  enddo
!
!  Free local memory
!
  deallocate(dNdr,stat=status)
  if (status/=0) call deallocate_error('brennerfd','dNdr')
  deallocate(dNconjdj,stat=status)
  if (status/=0) call deallocate_error('brennerfd','dNconjdj')
  deallocate(dNconjdi,stat=status)
  if (status/=0) call deallocate_error('brennerfd','dNconjdi')
  deallocate(d1Btotj,stat=status)
  if (status/=0) call deallocate_error('brennerfd','d1Btotj')
  deallocate(d1Btoti,stat=status)
  if (status/=0) call deallocate_error('brennerfd','d1Btoti')
  deallocate(d1j,stat=status)
  if (status/=0) call deallocate_error('brennerfd','d1j')
  deallocate(d1i,stat=status)
  if (status/=0) call deallocate_error('brennerfd','d1i')
  deallocate(d2Ndr2,stat=status)           
  if (status/=0) call deallocate_error('brennerfd','d2Ndr2')
  deallocate(d2Nconjdj2,stat=status)             
  if (status/=0) call deallocate_error('brennerfd','d2Nconjdj2')
  deallocate(d2Nconjdi2,stat=status)
  if (status/=0) call deallocate_error('brennerfd','d2Nconjdi2')
  deallocate(d2Btotj,stat=status)              
  if (status/=0) call deallocate_error('brennerfd','d2Btotj')
  deallocate(d2Btoti,stat=status)
  if (status/=0) call deallocate_error('brennerfd','d2Btoti')
!
  t2 = g_cpu_time()
  tbrenner = tbrenner + t2 - t1
!
  return
  end subroutine brennerfd

  subroutine unsetbrennerneigh
!
!  Frees memory for Brenner neighbour list
!
!  12/17 Created from brennermd.f90
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
!  Copyright Curtin University 2017
!
!  Julian Gale, CIC, Curtin University, December 2017
!
  use datatypes
  use times
  implicit none
!
!  Local variables
!
  integer(i4)                                      :: status
  real(dp)                                         :: g_cpu_time
  real(dp)                                         :: t1
  real(dp)                                         :: t2
!
  t1 = g_cpu_time()
!
!  Free local memory
!
  deallocate(zneigh,stat=status)
  if (status/=0) call deallocate_error('brennermd','zneigh')
  deallocate(yneigh,stat=status)
  if (status/=0) call deallocate_error('brennermd','yneigh')
  deallocate(xneigh,stat=status)
  if (status/=0) call deallocate_error('brennermd','xneigh')
  deallocate(rneigh,stat=status)
  if (status/=0) call deallocate_error('brennermd','rneigh')
  deallocate(ineigh,stat=status)
  if (status/=0) call deallocate_error('brennermd','ineigh')
  deallocate(neighno,stat=status)
  if (status/=0) call deallocate_error('brennermd','neighno')
  deallocate(nREBObond,stat=status)
  if (status/=0) call deallocate_error('brennermd','nREBObond')
  deallocate(nsneigh,stat=status)
  if (status/=0) call deallocate_error('brennermd','nsneigh')
  deallocate(nneigh,stat=status)
  if (status/=0) call deallocate_error('brennermd','nneigh')
  deallocate(nfreeatom,stat=status)
  if (status/=0) call deallocate_error('brennermd','nfreeatom')
  deallocate(lopanyneigh,stat=status)
  if (status/=0) call deallocate_error('brennermd','lopanyneigh')
  deallocate(latomdone,stat=status)
  if (status/=0) call deallocate_error('brennermd','latomdone')
  deallocate(lalreadydone,stat=status)
  if (status/=0) call deallocate_error('brennermd','lalreadydone')
  deallocate(nREBOatomRptr,stat=status)
  if (status/=0) call deallocate_error('brennermd','nREBOatomRptr')
  deallocate(nREBOatomptr,stat=status)
  if (status/=0) call deallocate_error('brennermd','nREBOatomptr')
!
  t2 = g_cpu_time()
  tbrenner = tbrenner + t2 - t1
!
  return
  end subroutine unsetbrennerneigh

end module m_brenner
