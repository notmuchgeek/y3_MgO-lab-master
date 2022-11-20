  module m_bondorder

    use datatypes

    integer(i4), dimension(:,:,:), allocatable, save :: ineigh         ! Cell indices for vector
    integer(i4)                                      :: maxneigh2
    integer(i4)                                      :: nboatom
    integer(i4), dimension(:,:),   allocatable, save :: nbopotptr
    integer(i4), dimension(:),     allocatable, save :: nboatomRptr
    integer(i4), dimension(:,:),   allocatable, save :: neighno
    integer(i4), dimension(:),     allocatable, save :: nfreeatom
    integer(i4), dimension(:),     allocatable, save :: nneigh
    logical,     dimension(:),     allocatable, save :: latomdone
    real(dp),    dimension(:),     allocatable, save :: rBOcutmax
    real(dp),    dimension(:,:),   allocatable, save :: rneigh
    real(dp),    dimension(:,:),   allocatable, save :: xneigh
    real(dp),    dimension(:,:),   allocatable, save :: yneigh
    real(dp),    dimension(:,:),   allocatable, save :: zneigh
    real(dp),    dimension(:),     allocatable, save :: Zcn

  contains

  subroutine setbondorderneigh
!
!  Calculates the neighbour list for bond order potentials
!
!  12/17 Created from bondordermd.f
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
  use bondorderdata
  use control,        only : keyword
  use current
  use iochannels
  use neighbours
  use parallel
  use spatialbo,      only : lspatialok => lspatialBOok
  use times
  implicit none
!
!  Local variables
!
  integer(i4)                                      :: i
  integer(i4)                                      :: itmp
  integer(i4)                                      :: j
  integer(i4)                                      :: nati
  integer(i4)                                      :: nmin
  integer(i4)                                      :: nn
  integer(i4)                                      :: nn2
  integer(i4)                                      :: nptr
  integer(i4)                                      :: ntypi
  integer(i4)                                      :: status
  logical                                          :: lmaxneighok
  real(dp)                                         :: g_cpu_time
  real(dp)                                         :: rtmp
  real(dp)                                         :: t1
  real(dp)                                         :: t2
!
  t1 = g_cpu_time()
!
  allocate(nboatomRptr(numat),stat=status)
  if (status/=0) call outofmemory('setbondorderneigh','nboatomRptr')
!
!  Set up a dummy pointer for derivative array calls
!
  nboatom = 0
  do i = 1,numat
    nboatom = nboatom + 1
    nboatomRptr(i) = nboatom
  enddo
!
!  Allocate memory that does not depend on maxneigh             
!
  allocate(latomdone(numat),stat=status)
  if (status/=0) call outofmemory('setbondorderneigh','latomdone')
  allocate(nfreeatom(numat),stat=status)
  if (status/=0) call outofmemory('setbondorderneigh','nfreeatom')
  allocate(nneigh(numat),stat=status)
  if (status/=0) call outofmemory('setbondorderneigh','nneigh')
  allocate(rBOcutmax(numat),stat=status)
  if (status/=0) call outofmemory('setbondorderneigh','rBOcutmax')
  allocate(Zcn(numat),stat=status)
  if (status/=0) call outofmemory('setbondorderneigh','Zcn')
!
!  Reinitialisation point should maxneigh be increased             
!
100 continue
  lmaxneighok = .true.
  if (allocated(rneigh)) then
    deallocate(zneigh,stat=status)
    if (status/=0) call deallocate_error('setbondorderneigh','zneigh')
    deallocate(yneigh,stat=status)
    if (status/=0) call deallocate_error('setbondorderneigh','yneigh')
    deallocate(xneigh,stat=status)
    if (status/=0) call deallocate_error('setbondorderneigh','xneigh')
    deallocate(rneigh,stat=status)
    if (status/=0) call deallocate_error('setbondorderneigh','rneigh')
    deallocate(ineigh,stat=status)
    if (status/=0) call deallocate_error('setbondorderneigh','ineigh')
    deallocate(nbopotptr,stat=status)
    if (status/=0) call deallocate_error('setbondorderneigh','nbopotptr')
    deallocate(neighno,stat=status)
    if (status/=0) call deallocate_error('setbondorderneigh','neighno')
  endif
!
!  Set parameter for pairwise storage memory
!
  maxneigh2 = maxneigh + maxneigh*(maxneigh + 1)/2
!
!  Allocate local memory
!
  allocate(neighno(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('setbondorderneigh','neighno')
  allocate(nbopotptr(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('setbondorderneigh','nbopotptr')
  allocate(ineigh(3_i4,maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('setbondorderneigh','ineigh')
  allocate(rneigh(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('setbondorderneigh','rneigh')
  allocate(xneigh(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('setbondorderneigh','xneigh')
  allocate(yneigh(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('setbondorderneigh','yneigh')
  allocate(zneigh(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('setbondorderneigh','zneigh')
!****************************
!  Find list of free atoms  *
!****************************
  do i = 1,numat
    nfreeatom(i) = i
  enddo
!*************************************
!  Find cut-off radii for all atoms  *
!*************************************
  do i = 1,numat
    nati = nat(i)
    ntypi = nftype(i)
    rBOcutmax(i) = 0.0_dp
!
!  Check twobody potentials
!
    do j = 1,nbopot
      if (nati.eq.nBOspec1(j).and.(ntypi.eq.nBOtyp1(j).or.nBOtyp1(j).eq.0)) then
        rBOcutmax(i) = max(rBOcutmax(i),rBOmax(j))
      endif
      if (nati.eq.nBOspec2(j).and.(ntypi.eq.nBOtyp2(j).or.nBOtyp2(j).eq.0)) then
        rBOcutmax(i) = max(rBOcutmax(i),rBOmax(j))
      endif
    enddo
  enddo
!****************************************
!  Calculate neighbour lists for atoms  *
!****************************************
!
!  Set up logical array of atoms done, so that only those needed are done in parallel
!
  latomdone(1:numat) = .false.
!
!  Compute neighbour list
!
  call getBOneighbour(maxneigh,rBOcutmax,nBOpotptr,nneigh,neighno,rneigh, &
                      xneigh,yneigh,zneigh,ineigh,latomdone,lmaxneighok)
!**********************************************************************
!  If maxneigh has been exceeded, increase value and return to start  *
!**********************************************************************
  if (.not.lmaxneighok) then
    do i = 1,numat
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
    do i = 1,numat
      if (latomdone(i)) then
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
            itmp = nbopotptr(nptr,i)
            nbopotptr(nptr,i) = nbopotptr(nn,i)
            nbopotptr(nn,i)  = itmp
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
      endif
    enddo
  endif
!*********************************************
!  Calculate numbers of neighbours for each  *
!*********************************************
  if (ioproc.and.index(keyword,'debu').ne.0) then
    write(ioout,'(/,''  Number of neighbours for atoms :'',/)')
    write(ioout,'(''  i'',2x,''Nat'',3x,''N'')')
    do i = 1,numat
      write(ioout,'(i4,1x,i3,1x,i6)') i,nat(i),nneigh(i)
    enddo
    write(ioout,'(/,''  Neighbours of atoms :'',/)')
    do i = 1,numat
      write(ioout,'(i4,8(1x,i4))') i,(neighno(nn,i),nn=1,nneigh(i))
    enddo
  endif
!
  t2 = g_cpu_time()
  tbondorder = tbondorder + t2 - t1
!
  return
  end subroutine setbondorderneigh

  subroutine bondorderfd(matom,mcrd,step,ebondorder,lgrad1)
!
!  Calculates the energy and up to first derivatives for the Bond Order potentials.
!  Finite difference version that focuses on derivatives of matom.
!
!  On entry : 
!
!  lgrad1          = if .true. calculate the first derivatives
!
!  On exit :
!
!  ebondorder      = the value of the energy contribution
!
!  12/17 Created from bondordermd.f
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!  11/20 Tersoff reorganised
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
!  Julian Gale, CIC, Curtin University, November 2020
!
  use datatypes
  use bondorderdata
  use configurations, only : nregionno, nregions, nregiontype, QMMMmode
  use control,        only : lseok
  use current
  use energies,       only : esregion12, esregion2
  use iochannels
  use neighbours
  use parallel
  use times
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: matom
  integer(i4), intent(in)                          :: mcrd
  real(dp),    intent(in)                          :: step
  real(dp),    intent(out)                         :: ebondorder
  logical,     intent(in)                          :: lgrad1
!
!  Local variables
!
  integer(i4)                                      :: i
  integer(i4)                                      :: j
  integer(i4)                                      :: k
  integer(i4)                                      :: kk
  integer(i4)                                      :: mA
  integer(i4)                                      :: mR
  integer(i4)                                      :: nati
  integer(i4)                                      :: natj
  integer(i4)                                      :: natk
  integer(i4)                                      :: nboij
  integer(i4)                                      :: nboAij
  integer(i4)                                      :: nboAik
  integer(i4)                                      :: nboOi
  integer(i4)                                      :: nboRij
  integer(i4)                                      :: nboRik
  integer(i4)                                      :: nboZi
  integer(i4)                                      :: ni
  integer(i4)                                      :: nj
  integer(i4)                                      :: njk
  integer(i4)                                      :: nn
  integer(i4)                                      :: npki   
  integer(i4)                                      :: nneighi2
  integer(i4)                                      :: nneighj2
  integer(i4)                                      :: nregioni
  integer(i4)                                      :: nregionj
  integer(i4)                                      :: nregiontypi
  integer(i4)                                      :: nregiontypj
  integer(i4)                                      :: ntypi
  integer(i4)                                      :: ntypj
  integer(i4)                                      :: ntypk
  integer(i4)                                      :: status
  logical                                          :: lfound
  logical                                          :: lneedBOcnderiv
  logical                                          :: lQMMMok
  logical                                          :: lreg2one
  logical                                          :: lreg2pair
  real(dp)                                         :: bijA
  real(dp)                                         :: bijR
  real(dp)                                         :: bijsumA
  real(dp)                                         :: bijsumA1
  real(dp)                                         :: bijsumAn1
  real(dp)                                         :: bijsumR
  real(dp)                                         :: bijsumR1
  real(dp)                                         :: bijsumRn1
  real(dp)                                         :: BOalpAi
  real(dp)                                         :: BOalpRi
  real(dp)                                         :: BOncoAi
  real(dp)                                         :: BOncoRi
  real(dp)                                         :: btotA
  real(dp)                                         :: btotR
  real(dp)                                         :: g_cpu_time
  real(dp),    dimension(:),     allocatable, save :: d1i
  real(dp),    dimension(:),     allocatable, save :: d1BtotiA
  real(dp),    dimension(:),     allocatable, save :: d1BtotiR
  real(dp)                                         :: dedZ
  real(dp)                                         :: dexpijkdr
  real(dp)                                         :: dfdr
  real(dp)                                         :: dfikdr
  real(dp)                                         :: d2fdr2
  real(dp)                                         :: d2fikdr2
  real(dp)                                         :: d3fdr3
  real(dp)                                         :: d3fikdr3
  real(dp)                                         :: dGijkdr(3)
  real(dp)                                         :: d2Gijkdr2(6)
  real(dp)                                         :: d3Gijkdr3(10)
  real(dp)                                         :: dZ
  real(dp)                                         :: deltaZ
  real(dp)                                         :: eij
  real(dp)                                         :: expijk
  real(dp)                                         :: f
  real(dp)                                         :: fik
  real(dp)                                         :: dfzdz
  real(dp)                                         :: d2fzdz2
  real(dp)                                         :: d3fzdz3
  real(dp)                                         :: fz
  real(dp)                                         :: Gijk
  real(dp)                                         :: omegaik
  real(dp)                                         :: rbijsumA1
  real(dp)                                         :: rbijsumR1
  real(dp)                                         :: rij
  real(dp)                                         :: rik
  real(dp)                                         :: rlambda
  real(dp)                                         :: drlambdadrij
  real(dp)                                         :: drlambdadrik
  real(dp)                                         :: RmA
  real(dp)                                         :: RmR
  real(dp)                                         :: rrij
  real(dp)                                         :: rrik
  real(dp)                                         :: rtmp
  real(dp)                                         :: t1
  real(dp)                                         :: t2
  real(dp)                                         :: Va
  real(dp)                                         :: Vr
  real(dp)                                         :: dVadr
  real(dp)                                         :: dVrdr
  real(dp)                                         :: xdiff
  real(dp)                                         :: ydiff
  real(dp)                                         :: zdiff
  real(dp)                                         :: xji
  real(dp)                                         :: yji
  real(dp)                                         :: zji
  real(dp)                                         :: xki
  real(dp)                                         :: yki
  real(dp)                                         :: zki
  real(dp)                                         :: zAi3
  real(dp)                                         :: zRi3
  real(dp)                                         :: zsign
!
  t1 = g_cpu_time()
!
!  Initialise Bond Order energy
!
  ebondorder = 0.0_dp
!
!  Allocate local memory
!
  if (lgrad1) then
    allocate(d1i(maxneigh2),stat=status)
    if (status/=0) call outofmemory('bondorderfd','d1i')
    allocate(d1BtotiA(maxneigh2),stat=status)
    if (status/=0) call outofmemory('bondorderfd','d1BtotiA')
    allocate(d1BtotiR(maxneigh2),stat=status)
    if (status/=0) call outofmemory('bondorderfd','d1BtotiR')
  else
    allocate(d1i(1),stat=status)
    if (status/=0) call outofmemory('bondorderfd','d1i')
    allocate(d1BtotiA(1),stat=status)
    if (status/=0) call outofmemory('bondorderfd','d1BtotiA')
    allocate(d1BtotiR(1),stat=status)
    if (status/=0) call outofmemory('bondorderfd','d1BtotiR')
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
!
!  Initialise coordination number
!
  Zcn(1:numat) = 0.0_dp
!*****************************
!  Loop over pairs of atoms  *
!*****************************
  do i = 1,numat
!
!  Set variables relating to i
!
    nati = nat(i)
    ntypi = nftype(i)
    nregioni = nregionno(nsft + nrelf2a(i))
    nregiontypi = nregiontype(nregioni,ncf)
!
!  Set total number of distances for neighbours of i
!
    nneighi2 = nneigh(i) + nneigh(i)*(nneigh(i) + 1)/2
!
!  Initialise derivative storage for neighbours of i
!
    if (lgrad1) then
      d1i(1:nneighi2) = 0.0_dp
    endif
!
!  Check for self energy terms for atom i
!
    if (nboZ.gt.0) then
      lfound = .false.
      nboZi = 0
      do while (.not.lfound.and.nboZi.lt.nboZ)
        nboZi = nboZi + 1
        if (nBOspecZ(nboZi).eq.nati) then
          if (nBOtypZ(nboZi).eq.ntypi.or.nBOtypZ(nboZi).eq.0) then
            lfound = .true.
          endif
        endif
      enddo
      if (lfound) then
        ebondorder = ebondorder + BOecoeffZ(nboZi)
      endif
    endif
!
!  Find one-body parameters for i
!
    lfound = .false.
    nboOi = 0
    do while (.not.lfound.and.nboOi.lt.nboO)
      nboOi = nboOi + 1
      if (nBOspec0(nboOi).eq.nati) then
        if (nBOtyp0(nboOi).eq.ntypi.or.nBOtyp0(nboOi).eq.0) then
          lfound = .true.
        endif
      endif
    enddo
    if (lfound) then
      BOncoAi = BOncoeffA(nboOi)
      BOalpAi = BOecoeffA(nboOi)
      BOncoRi = BOncoeffR(nboOi)
      BOalpRi = BOecoeffR(nboOi)
    else
      BOncoAi = 1.0_dp
      BOalpAi = 0.0_dp
      BOncoRi = 1.0_dp
      BOalpRi = 0.0_dp
    endif
    zAi3 = BOalpAi**BOncoAi
    zRi3 = BOalpRi**BOncoRi
!
!  Loop over neighbours of i
!
    ni = 1
    niloop: do while (ni.le.nneigh(i))
      j = neighno(ni,i)
!
!  Set variables relating to j
!
      natj = nat(j)
      ntypj = nftype(j)
      nregionj = nregionno(nsft + nrelf2a(j))
      nregiontypj = nregiontype(nregionj,ncf)
!
!  Set up i-j quantities
!
      rij = rneigh(ni,i)
      xji = xneigh(ni,i)
      yji = yneigh(ni,i)
      zji = zneigh(ni,i)
      rrij = 1.0_dp/rij
!
      lreg2one  = .false.
      lreg2pair = .false.
      if (lseok.and.nregions(ncf).gt.1) then
        lreg2pair = (nregioni.gt.1.and.nregionj.gt.1)
        if (.not.lreg2pair) lreg2one = (nregioni.gt.1.or.nregionj.gt.1)
      endif
!
!  Find i in neighbour list for j
!
      nj = 1
      lfound = .false.
      do while (nj.lt.nneigh(j).and..not.lfound)
        if (neighno(nj,j).eq.i) then
          xdiff = xneigh(nj,j) + xji
          ydiff = yneigh(nj,j) + yji
          zdiff = zneigh(nj,j) + zji
          lfound = ((abs(xdiff)+abs(ydiff)+abs(zdiff)).lt.1.0d-6)
        endif
        if (.not.lfound) nj = nj + 1
      enddo
!
!  Set total number of distances for neighbours of j
!
      nneighj2 = nneigh(j) + nneigh(j)*(nneigh(j) + 1)/2
!
!  Find repulsive bond order potential from j to i
!
      lfound = .false.
      nboRij = 0
      do while (.not.lfound.and.nboRij.lt.nboR)
        nboRij = nboRij + 1
        if (nBOspecR1(nboRij).eq.nati.and.nBOspecR2(nboRij).eq.natj) then
          if ((nBOtypR1(nboRij).eq.ntypi.or.nBOtypR1(nboRij).eq.0).and. &
              (nBOtypR2(nboRij).eq.ntypj.or.nBOtypR2(nboRij).eq.0)) then
            lfound = .true.
          elseif (nati.eq.natj) then
            if ((nBOtypR2(nboRij).eq.ntypi.or.nBOtypR2(nboRij).eq.0).and. &
                (nBOtypR1(nboRij).eq.ntypj.or.nBOtypR1(nboRij).eq.0)) then
              lfound = .true.
            endif
          endif
        endif
      enddo
      if (.not.lfound) nboRij = 0
!
!  Find attractive bond order potential from j to i
!
      lfound = .false.
      nboAij = 0
      do while (.not.lfound.and.nboAij.lt.nboA)
        nboAij = nboAij + 1
        if (nBOspecA1(nboAij).eq.nati.and.nBOspecA2(nboAij).eq.natj) then
          if ((nBOtypA1(nboAij).eq.ntypi.or.nBOtypA1(nboAij).eq.0).and. &
              (nBOtypA2(nboAij).eq.ntypj.or.nBOtypA2(nboAij).eq.0)) then
            lfound = .true.
          elseif (nati.eq.natj) then
            if ((nBOtypA2(nboAij).eq.ntypi.or.nBOtypA2(nboAij).eq.0).and. &
                (nBOtypA1(nboAij).eq.ntypj.or.nBOtypA1(nboAij).eq.0)) then
              lfound = .true.
            endif
          endif
        endif
      enddo
      if (.not.lfound) nboAij = 0
!
!  Find two-body bond order potential between i and j
!
      lfound = .false.
      nboij = 0
      do while (.not.lfound.and.nboij.lt.nbopot) 
        nboij = nboij + 1
        if (nBOspec1(nboij).eq.nati.and.nBOspec2(nboij).eq.natj) then
          if ((nBOtyp1(nboij).eq.ntypi.or.nBOtyp1(nboij).eq.0).and. &
              (nBOtyp2(nboij).eq.ntypj.or.nBOtyp2(nboij).eq.0)) then
            lfound = .true.
          endif
        elseif (nBOspec1(nboij).eq.natj.and.nBOspec2(nboij).eq.nati) then
          if ((nBOtyp1(nboij).eq.ntypj.or.nBOtyp1(nboij).eq.0).and. &
              (nBOtyp2(nboij).eq.ntypi.or.nBOtyp2(nboij).eq.0)) then
            lfound = .true.
          endif
        endif
      enddo
      if (.not.lfound) nboij = 0
!
!  QM/MM handling : i & j are both QM atoms => exclude
!
      lQMMMok = .true.
      if (QMMMmode(ncf).gt.0) then
        if (nregiontypi.eq.1.and.nregiontypj.eq.1) lQMMMok = .false.
      endif
      if (nboij.gt.0.and.lQMMMok) then
!****************************************
!  Valid two-body bond order potential  *
!****************************************
!
!  Calculate fij
!
        if (nBOtapertype(nboij).eq.2) then
          call ctaper(rij,rBOmin(nboij),rBOmax(nboij),f,dfdr,d2fdr2,d3fdr3,lgrad1,.false.,.false.)
        elseif (nBOtapertype(nboij).eq.3) then
          call murtytaper(rij,rBOmin(nboij),rBOmax(nboij),f,dfdr,d2fdr2,d3fdr3,lgrad1,.false.,.false.)
        else
          call mdftaper(rij,rBOmin(nboij),rBOmax(nboij),f,dfdr,d2fdr2,d3fdr3,lgrad1,.false.,.false.)
        endif
!
!  Calculate Bij and Bji - loop over all other neighbours
!
        bijsumA = 0.0_dp
        bijsumR = 0.0_dp
        if (lgrad1) then
          d1BtotiA(1:nneighi2) = 0.0_dp
          d1BtotiR(1:nneighi2) = 0.0_dp
        endif
!
!  Loop over neighbours of i .ne. j 
!
        kloop: do k = 1,nneigh(i)
          npki = nbopotptr(k,i)
          if (k.ne.ni) then
            rik = rneigh(k,i)
            xki = xneigh(k,i)
            yki = yneigh(k,i)
            zki = zneigh(k,i)
!
!  Repulsive component
!
            if (rik.lt.rBOmax(npki)) then
!
!  Find nboRik
!
              lfound = .false.
              nboRik = 0
              kk = neighno(k,i)
              natk = nat(kk)
              ntypk = nftype(kk)
              do while (.not.lfound.and.nboRik.lt.nboR)
                nboRik = nboRik + 1
                if (nBOspecR1(nboRik).eq.nati.and.nBOspecR2(nboRik).eq.natk) then
                  if ((nBOtypR1(nboRik).eq.ntypi.or.nBOtypR1(nboRik).eq.0).and. &
                      (nBOtypR2(nboRik).eq.ntypk.or.nBOtypR2(nboRik).eq.0)) then
                    lfound = .true.
                  endif
                endif
              enddo
              if (nboRik.gt.0) then
!
!  Calculate fik
!
                if (nBOtapertype(npki).eq.2) then
                  call ctaper(rik,rBOmin(npki),rBOmax(npki),fik,dfikdr,d2fikdr2,d3fikdr3,lgrad1,.false.,.false.)
                elseif (nBOtapertype(npki).eq.3) then
                  call murtytaper(rik,rBOmin(npki),rBOmax(npki),fik,dfikdr,d2fikdr2,d3fikdr3,lgrad1,.false.,.false.)
                else
                  call mdftaper(rik,rBOmin(npki),rBOmax(npki),fik,dfikdr,d2fikdr2,d3fikdr3,lgrad1,.false.,.false.)
                endif
! 
!  Set omega for i-k
!
                omegaik = BOocoeffR(nboRik)
                fik = fik*omegaik
                if (lgrad1) then
                  dfikdr = dfikdr*omegaik
                endif
!
!  Calculate Gijk
!
                if (nBOtypeR(nboRik).ne.1) then
                  call GthetaBO(xji,yji,zji,xki,yki,zki,nBOtypeR(nboRik),BOccoeffR(1,nboRik), &
                    BOhcoeffR(nboRik),Gijk,dGijkdr,d2Gijkdr2,d3Gijkdr3,lgrad1,.false.,.false.)
                else
                  Gijk = 1.0_dp
                  dGijkdr = 0.0_dp
                  d2Gijkdr2 = 0.0_dp
                endif
!
!  Calculate exponential factor
!
                mR = nint(BOmcoeffR(nboRik))
                if (nboRij.gt.0) then
                  if (lBOzrlR(nboRij)) then
                    rlambda = (BOlcoeffR(nboRij)*rij - BOlcoeffR(nboRik)*rik)
                  else
                    rlambda = BOlcoeffR(nboRik)*(rij - rik)
                  endif
                else
                  rlambda = BOlcoeffR(nboRik)*(rij - rik)
                endif
                expijk = exp(rlambda**mR)
!
!  Combine terms
!
                bijsumR = bijsumR + Gijk*fik*expijk
                if (lgrad1) then
!
!  Derivatives
!
!  Find index for j-k 
!
                  if (ni.ge.k) then
                    njk = nneigh(i) + ni*(ni-1)/2 + k
                  else
                    njk = nneigh(i) + k*(k-1)/2 + ni
                  endif
!
                  rrik = 1.0_dp/rik
                  dfikdr = rrik*dfikdr
                  RmR = dble(mR)
                  if (mR.ge.1) then
                    if (nboRij.gt.0) then
                      if (lBOzrlR(nboRij)) then
                        drlambdadrij = BOlcoeffR(nboRij)
                        drlambdadrik = BOlcoeffR(nboRik)
                      else
                        drlambdadrij = BOlcoeffR(nboRik)
                        drlambdadrik = BOlcoeffR(nboRik)
                      endif
                    else
                      drlambdadrij = BOlcoeffR(nboRik)
                      drlambdadrik = BOlcoeffR(nboRik)
                    endif
                    dexpijkdr = RmR*(rlambda**(mR-1))*expijk
                  else
                    dexpijkdr = 0.0_dp
                  endif
!
                  d1BtotiR(k)   = d1BtotiR(k)   + Gijk*dfikdr*expijk
!
                  d1BtotiR(ni)  = d1BtotiR(ni)  + dGijkdr(1)*fik*expijk
                  d1BtotiR(k)   = d1BtotiR(k)   + dGijkdr(2)*fik*expijk
                  d1BtotiR(njk) = d1BtotiR(njk) + dGijkdr(3)*fik*expijk
!
                  d1BtotiR(ni) = d1BtotiR(ni) + Gijk*fik*dexpijkdr*rrij*drlambdadrij
                  d1BtotiR(k)  = d1BtotiR(k)  - Gijk*fik*dexpijkdr*rrik*drlambdadrik
                endif
              endif
            endif
!
!  Attractive component
!
            if (rik.lt.rBOmax(npki)) then
!
!  Find nboAik
!
              lfound = .false.
              nboAik = 0
              kk = neighno(k,i)
              natk = nat(kk)
              ntypk = nftype(kk)
              do while (.not.lfound.and.nboAik.lt.nboA)
                nboAik = nboAik + 1
                if (nBOspecA1(nboAik).eq.nati.and.nBOspecA2(nboAik).eq.natk) then
                  if ((nBOtypA1(nboAik).eq.ntypi.or.nBOtypA1(nboAik).eq.0).and. &
                      (nBOtypA2(nboAik).eq.ntypk.or.nBOtypA2(nboAik).eq.0)) then
                    lfound = .true.
                  endif
                endif
              enddo
              if (nboAik.gt.0) then
!
!  Calculate fik
!
                if (nBOtapertype(npki).eq.2) then
                  call ctaper(rik,rBOmin(npki),rBOmax(npki),fik,dfikdr,d2fikdr2,d3fikdr3,lgrad1,.false.,.false.)
                elseif (nBOtapertype(npki).eq.3) then
                  call murtytaper(rik,rBOmin(npki),rBOmax(npki),fik,dfikdr,d2fikdr2,d3fikdr3,lgrad1,.false.,.false.)
                else
                  call mdftaper(rik,rBOmin(npki),rBOmax(npki),fik,dfikdr,d2fikdr2,d3fikdr3,lgrad1,.false.,.false.)
                endif
! 
!  Set omega for i-k
!
                omegaik = BOocoeffA(nboAik)
                fik = fik*omegaik
                if (lgrad1) then
                  dfikdr = dfikdr*omegaik
                endif
!
!  Calculate Gijk
!
                if (nBOtypeA(nboAik).ne.1) then
                  call GthetaBO(xji,yji,zji,xki,yki,zki,nBOtypeA(nboAik),BOccoeffA(1,nboAik), &
                    BOhcoeffA(nboAik),Gijk,dGijkdr,d2Gijkdr2,d3Gijkdr3,lgrad1,.false.,.false.)
                else
                  Gijk = 1.0_dp
                  dGijkdr = 0.0_dp
                  d2Gijkdr2 = 0.0_dp
                endif
!
!  Calculate exponential factor
!
                mA = nint(BOmcoeffA(nboAik))
                if (nboAij.gt.0) then
                  if (lBOzrlA(nboAij)) then
                    rlambda = (BOlcoeffA(nboAij)*rij - BOlcoeffA(nboAik)*rik)
                  else
                    rlambda = BOlcoeffA(nboAik)*(rij - rik)
                  endif
                else
                  rlambda = BOlcoeffA(nboAik)*(rij - rik)
                endif
                expijk = exp(rlambda**mA)
!
!  Combine terms
!
                bijsumA = bijsumA + Gijk*fik*expijk
                if (lgrad1) then
!
!  Derivatives
!
!  Find index for j-k 
!
                  if (ni.ge.k) then
                    njk = nneigh(i) + ni*(ni-1)/2 + k
                  else
                    njk = nneigh(i) + k*(k-1)/2 + ni
                  endif
!
                  rrik = 1.0_dp/rik
                  dfikdr = rrik*dfikdr
                  RmA = dble(mA)
                  if (mA.ge.1) then
                    if (nboAij.gt.0) then
                      if (lBOzrlA(nboAij)) then
                        drlambdadrij = BOlcoeffA(nboAij)
                        drlambdadrik = BOlcoeffA(nboAik)
                      else
                        drlambdadrij = BOlcoeffA(nboAik)
                        drlambdadrik = BOlcoeffA(nboAik)
                      endif
                    else
                      drlambdadrij = BOlcoeffA(nboAik)
                      drlambdadrik = BOlcoeffA(nboAik)
                    endif
                    dexpijkdr = RmA*(rlambda**(mA-1))*expijk
                  else
                    dexpijkdr = 0.0_dp
                  endif
!
                  d1BtotiA(k)   = d1BtotiA(k)   + Gijk*dfikdr*expijk
!
                  d1BtotiA(ni)  = d1BtotiA(ni)  + dGijkdr(1)*fik*expijk
                  d1BtotiA(k)   = d1BtotiA(k)   + dGijkdr(2)*fik*expijk
                  d1BtotiA(njk) = d1BtotiA(njk) + dGijkdr(3)*fik*expijk
!
                  d1BtotiA(ni) = d1BtotiA(ni) + Gijk*fik*dexpijkdr*rrij*drlambdadrij
                  d1BtotiA(k)  = d1BtotiA(k)  - Gijk*fik*dexpijkdr*rrik*drlambdadrik
                endif
              endif
            endif
          endif
        enddo kloop
!
!  Raise terms to the power of n, add 1, and then raise to -2*n
!
        if (abs(bijsumA).gt.1.0d-12) then
          bijsumAn1 = bijsumA**(BOncoAi - 1.0_dp)
        else
          bijsumAn1 = 0.0_dp
        endif
        if (abs(bijsumR).gt.1.0d-12) then
          bijsumRn1 = bijsumR**(BOncoRi - 1.0_dp)
        else
          bijsumRn1 = 0.0_dp
        endif
!
        bijsumA1 = 1.0_dp + zAi3*bijsumA*bijsumAn1
        bijsumR1 = 1.0_dp + zRi3*bijsumR*bijsumRn1
!
        rbijsumA1 = 1.0_dp/bijsumA1
        rbijsumR1 = 1.0_dp/bijsumR1
!
        bijA = bijsumA1**(-0.5_dp/BOncoAi)
        bijR = bijsumR1**(-0.5_dp/BOncoRi)
!
!  Scale derivatives by bijsum factors
!
        if (lgrad1) then
!
!  First derivatives
!
          if (bijsumA.gt.0.0_dp) then
            rtmp = 0.25_dp*zAi3*bijsumAn1*bijA*rbijsumA1
            do nn = 1,nneighi2
              d1BtotiA(nn) = - rtmp*d1BtotiA(nn)
            enddo
          endif
          if (bijsumR.gt.0.0_dp) then
            rtmp = 0.25_dp*zRi3*bijsumRn1*bijR*rbijsumR1
            do nn = 1,nneighi2
              d1BtotiR(nn) = - rtmp*d1BtotiR(nn)
            enddo
          endif
        endif
!
!  Calculate two-body component of potential
!
        Vr = BOacoeff(nboij)*exp(-BOzacoeff(nboij)*rij)
        Va = BObcoeff(nboij)*exp(-BOzbcoeff(nboij)*rij)
        if (lgrad1) then
          dVrdr = - BOzacoeff(nboij)*Vr
          dVadr = - BOzbcoeff(nboij)*Va
!
          dVrdr = rrij*dVrdr
          dVadr = rrij*dVadr
        endif
!
!  Calculate total i-j potential
!
        BtotA = 0.5_dp*bijA
        BtotR = 0.5_dp*bijR
        eij = f*(BtotR*Vr - BtotA*Va)
!
!  Add to surface energy totals if appropriate
!
        if (lseok) then
          if (lreg2one) then
            esregion12 = esregion12 + eij
          elseif (lreg2pair) then
            esregion2 = esregion2 + eij
          else
            ebondorder = ebondorder + eij
          endif
        else
          ebondorder = ebondorder + eij
        endif
!
!  Add contribution to coordination number
!
        Zcn(i) = Zcn(i) + f*bijA
!
!  Derivatives of Bond Order potential energy
!
        if (lgrad1) then
          dfdr = rrij*dfdr
          d1i(ni) = d1i(ni) + dfdr*(BtotR*Vr - BtotA*Va)
          d1i(ni) = d1i(ni) + f*(BtotR*dVrdr - BtotA*dVadr)
          do nn = 1,nneighi2
            d1i(nn) = d1i(nn) + f*(Vr*d1BtotiR(nn) - Va*d1BtotiA(nn))
          enddo
        endif
      endif
!
!  End of loop over neighbours of i
!
      ni = ni + 1
    enddo niloop
!
!  Add derivatives due to neighbours of i
!
    if (lgrad1) then
      call d1add(i,maxneigh,maxneigh,nneigh,neighno,xneigh,yneigh,zneigh,nfreeatom,nboatomRptr,d1i,.false.,.true.)
    endif
  enddo
  if (nboZ.gt.0) then
!**************************************
!  Compute coordination contribution  *
!**************************************
    do i = 1,numat
!
!  Set variables relating to i
!
      nati = nat(i)
      ntypi = nftype(i)
!
!  Check for coordination terms for atom i
!
      lfound = .false.
      nboZi = 0
      do while (.not.lfound.and.nboZi.lt.nboZ)
        nboZi = nboZi + 1
        if (nBOspecZ(nboZi).eq.nati) then
          if (nBOtypZ(nboZi).eq.ntypi.or.nBOtypZ(nboZi).eq.0) then
            lfound = .true.
          endif
        endif
      enddo
      if (lfound) then
        dZ = Zcn(i) - BOzcoeffZ(nboZi)
        zsign = sign(1.0_dp,dZ)
        call botaper(dZ,fz,dfzdz,d2fzdz2,d3fzdz3,lgrad1,.false.,.false.)
        deltaZ = zsign*fz
        lneedBOcnderiv = (abs(dble(nint(deltaZ)) - deltaZ).gt.BOcntol)
        ebondorder = ebondorder + deltaZ*(BOccoeffZ(1,nboZi) + BOccoeffZ(2,nboZi)*deltaZ)
!
!  If deltaZ is not an integer then compute derivatives
!
        if (lneedBOcnderiv.and.lgrad1) then
          dedZ = (BOccoeffZ(1,nboZi) + 2.0_dp*BOccoeffZ(2,nboZi)*deltaZ)*zsign*dfzdz
!
          nregioni = nregionno(nsft + nrelf2a(i))
          nregiontypi = nregiontype(nregioni,ncf)
!
!  Set total number of distances for neighbours of i
!
          nneighi2 = nneigh(i) + nneigh(i)*(nneigh(i) + 1)/2 
!
!  Initialise derivative storage for neighbours of i
!
          d1i(1:nneighi2) = 0.0_dp
!
!  Find one-body parameters for i
!
          lfound = .false.
          nboOi = 0
          do while (.not.lfound.and.nboOi.lt.nboO)
            nboOi = nboOi + 1
            if (nBOspec0(nboOi).eq.nati) then
              if (nBOtyp0(nboOi).eq.ntypi.or.nBOtyp0(nboOi).eq.0) then
                lfound = .true.
              endif
            endif
          enddo
          if (lfound) then
            BOncoAi = BOncoeffA(nboOi)
            BOalpAi = BOecoeffA(nboOi)
          else
            BOncoAi = 1.0_dp
            BOalpAi = 0.0_dp
          endif
          zAi3 = BOalpAi**BOncoAi
!
!  Loop over neighbours of i 
!
          ni = 1
          niloopnc: do while (ni.le.nneigh(i))
            j = neighno(ni,i)
!
!  Set variables relating to j
!
            natj = nat(j)
            ntypj = nftype(j)
            nregionj = nregionno(nsft + nrelf2a(j))
            nregiontypj = nregiontype(nregionj,ncf)
!
!  Set up i-j quantities
!
            rij = rneigh(ni,i)
            xji = xneigh(ni,i)
            yji = yneigh(ni,i)
            zji = zneigh(ni,i)
            rrij = 1.0_dp/rij
!
!  Find i in neighbour list for j
!
            nj = 1
            lfound = .false.
            do while (nj.lt.nneigh(j).and..not.lfound)
              if (neighno(nj,j).eq.i) then
                xdiff = xneigh(nj,j) + xji
                ydiff = yneigh(nj,j) + yji
                zdiff = zneigh(nj,j) + zji
                lfound = ((abs(xdiff)+abs(ydiff)+abs(zdiff)).lt.1.0d-6)
              endif
              if (.not.lfound) nj = nj + 1
            enddo
!
!  Find two-body bond order potential between i and j
!
            lfound = .false.
            nboij = 0
            do while (.not.lfound.and.nboij.lt.nbopot) 
              nboij = nboij + 1
              if (nBOspec1(nboij).eq.nati.and.nBOspec2(nboij).eq.natj) then
                if ((nBOtyp1(nboij).eq.ntypi.or.nBOtyp1(nboij).eq.0).and. &
                    (nBOtyp2(nboij).eq.ntypj.or.nBOtyp2(nboij).eq.0)) then
                  lfound = .true.
                endif
              elseif (nBOspec1(nboij).eq.natj.and.nBOspec2(nboij).eq.nati) then
                if ((nBOtyp1(nboij).eq.ntypj.or.nBOtyp1(nboij).eq.0).and. &
                    (nBOtyp2(nboij).eq.ntypi.or.nBOtyp2(nboij).eq.0)) then
                  lfound = .true.
                endif
              endif
            enddo
            if (.not.lfound) nboij = 0
!
!  Find attractive bond order potential from j to i
!
            lfound = .false.
            nboAij = 0
            do while (.not.lfound.and.nboAij.lt.nboA)
              nboAij = nboAij + 1
              if (nBOspecA1(nboAij).eq.nati.and.nBOspecA2(nboAij).eq.natj) then
                if ((nBOtypA1(nboAij).eq.ntypi.or.nBOtypA1(nboAij).eq.0).and. &
                    (nBOtypA2(nboAij).eq.ntypj.or.nBOtypA2(nboAij).eq.0)) then
                  lfound = .true.
                elseif (nati.eq.natj) then
                  if ((nBOtypA2(nboAij).eq.ntypi.or.nBOtypA2(nboAij).eq.0).and. &
                      (nBOtypA1(nboAij).eq.ntypj.or.nBOtypA1(nboAij).eq.0)) then
                    lfound = .true.
                  endif
                endif
              endif
            enddo
            if (.not.lfound) nboAij = 0
!
!  QM/MM handling : i & j are both QM atoms => exclude
!
            lQMMMok = .true.
            if (QMMMmode(ncf).gt.0) then
              if (nregiontypi.eq.1.and.nregiontypj.eq.1) lQMMMok = .false.
            endif
            if (nboij.gt.0.and.lQMMMok) then
!****************************************
!  Valid two-body bond order potential  *
!****************************************
!
!  Calculate fij
!
              if (nBOtapertype(nboij).eq.2) then
                call ctaper(rij,rBOmin(nboij),rBOmax(nboij),f,dfdr,d2fdr2,d3fdr3,lgrad1,.false.,.false.)
              elseif (nBOtapertype(nboij).eq.3) then
                call murtytaper(rij,rBOmin(nboij),rBOmax(nboij),f,dfdr,d2fdr2,d3fdr3,lgrad1,.false.,.false.)
              else
                call mdftaper(rij,rBOmin(nboij),rBOmax(nboij),f,dfdr,d2fdr2,d3fdr3,lgrad1,.false.,.false.)
              endif
!
!  Calculate Bij and Bji - loop over all other neighbours
!
              bijsumA = 0.0_dp
              d1BtotiA(1:nneighi2) = 0.0_dp
!
!  Loop over neighbours of i .ne. j 
!
              do k = 1,nneigh(i)
                npki = nbopotptr(k,i)
                if (k.ne.ni) then
                  rik = rneigh(k,i)
                  xki = xneigh(k,i)
                  yki = yneigh(k,i)
                  zki = zneigh(k,i)
!
!  Attractive component
!
                    if (rik.lt.rBOmax(npki)) then
!
!  Find nboAik
!
                      lfound = .false.
                      nboAik = 0
                      kk = neighno(k,i)
                      natk = nat(kk)
                      ntypk = nftype(kk)
                      do while (.not.lfound.and.nboAik.lt.nboA)
                        nboAik = nboAik + 1
                        if (nBOspecA1(nboAik).eq.nati.and.nBOspecA2(nboAik).eq.natk) then
                          if ((nBOtypA1(nboAik).eq.ntypi.or.nBOtypA1(nboAik).eq.0).and. &
                              (nBOtypA2(nboAik).eq.ntypk.or.nBOtypA2(nboAik).eq.0)) then
                            lfound = .true.
                          endif
                        endif
                      enddo
                      if (nboAik.gt.0) then
!
!  Calculate fik
!
                      if (nBOtapertype(npki).eq.2) then
                        call ctaper(rik,rBOmin(npki),rBOmax(npki),fik,dfikdr,d2fikdr2,d3fikdr3,lgrad1,.false.,.false.)
                      elseif (nBOtapertype(npki).eq.3) then
                        call murtytaper(rik,rBOmin(npki),rBOmax(npki),fik,dfikdr,d2fikdr2,d3fikdr3,lgrad1,.false.,.false.)
                      else
                        call mdftaper(rik,rBOmin(npki),rBOmax(npki),fik,dfikdr,d2fikdr2,d3fikdr3,lgrad1,.false.,.false.)
                      endif
! 
!  Set omega for i-k
!
                      omegaik = BOocoeffA(nboAik)
                      fik = fik*omegaik
                      if (lgrad1) then
                        dfikdr = dfikdr*omegaik
                      endif
!
!  Calculate Gijk
!
                      if (nBOtypeA(nboAik).ne.1) then
                        call GthetaBO(xji,yji,zji,xki,yki,zki,nBOtypeA(nboAik),BOccoeffA(1,nboAik), &
                                      BOhcoeffA(nboAik),Gijk,dGijkdr,d2Gijkdr2,d3Gijkdr3,lgrad1,.false.,.false.)
                      else
                        Gijk = 1.0_dp
                        dGijkdr = 0.0_dp
                        d2Gijkdr2 = 0.0_dp
                      endif
!
!  Calculate exponential factor
!
                      mA = nint(BOmcoeffA(nboAik))
                      if (nboAij.gt.0) then
                        if (lBOzrlA(nboAij)) then
                          rlambda = (BOlcoeffA(nboAij)*rij - BOlcoeffA(nboAik)*rik)
                        else
                          rlambda = BOlcoeffA(nboAik)*(rij - rik)
                        endif
                      else
                        rlambda = BOlcoeffA(nboAik)*(rij - rik)
                      endif
                      expijk = exp(rlambda**mA)
!
!  Combine terms
!
                      bijsumA = bijsumA + Gijk*fik*expijk
!
!  Derivatives
!
!  Find index for j-k 
!
                      if (ni.ge.k) then
                        njk = nneigh(i) + ni*(ni-1)/2 + k
                      else
                        njk = nneigh(i) + k*(k-1)/2 + ni
                      endif
!
                      rrik = 1.0_dp/rik
                      dfikdr = rrik*dfikdr
                      RmA = dble(mA)
                      if (mA.ge.1) then
                        if (nboAij.gt.0) then
                          if (lBOzrlA(nboAij)) then
                            drlambdadrij = BOlcoeffA(nboAij)
                            drlambdadrik = BOlcoeffA(nboAik)
                          else
                            drlambdadrij = BOlcoeffA(nboAik)
                            drlambdadrik = BOlcoeffA(nboAik)
                          endif
                        else
                          drlambdadrij = BOlcoeffA(nboAik)
                          drlambdadrik = BOlcoeffA(nboAik)
                        endif
                        dexpijkdr = RmA*(rlambda**(mA-1))*expijk
                      else
                        dexpijkdr = 0.0_dp
                      endif
!
                      d1BtotiA(k)   = d1BtotiA(k)   + Gijk*dfikdr*expijk
!
                      d1BtotiA(ni)  = d1BtotiA(ni)  + dGijkdr(1)*fik*expijk
                      d1BtotiA(k)   = d1BtotiA(k)   + dGijkdr(2)*fik*expijk
                      d1BtotiA(njk) = d1BtotiA(njk) + dGijkdr(3)*fik*expijk
!
                      d1BtotiA(ni) = d1BtotiA(ni) + Gijk*fik*dexpijkdr*rrij*drlambdadrij
                      d1BtotiA(k)  = d1BtotiA(k)  - Gijk*fik*dexpijkdr*rrik*drlambdadrik
                    endif
                  endif
                endif
              enddo
!
              if (abs(bijsumA).gt.1.0d-12) then
                bijsumAn1 = bijsumA**(BOncoAi - 1.0_dp)
              else
                bijsumAn1 = 0.0_dp
              endif
!
              bijsumA1 = 1.0_dp + zAi3*bijsumA*bijsumAn1
!
              rbijsumA1 = 1.0_dp/bijsumA1
!
              bijA = bijsumA1**(-0.5_dp/BOncoAi)
!
!  Scale derivatives by bijsum factors
!
!  First derivatives
!
              if (bijsumA.gt.0.0_dp) then
                rtmp = 0.25_dp*zAi3*bijsumAn1*bijA*rbijsumA1
                do nn = 1,nneighi2
                  d1BtotiA(nn) = - rtmp*d1BtotiA(nn)
                enddo
              endif
!
!  Derivatives of Bond Order coordination energy
!
              dfdr = rrij*dfdr
              d1i(ni) = d1i(ni) + dedZ*dfdr*bijA
              do nn = 1,nneighi2
                d1i(nn) = d1i(nn) + 2.0_dp*dedZ*f*d1BtotiA(nn)
              enddo
            endif
!
!  End of loop over neighbours of i
!
            ni = ni + 1
          enddo niloopnc
!
!  Add derivatives due to neighbours of i
!
          call d1add(i,maxneigh,maxneigh,nneigh,neighno,xneigh,yneigh,zneigh,nfreeatom,nboatomRptr,d1i,.false.,.true.)
!
!  End test over whether derivatives are needed
!
        endif
!
!  End test over whether a coordination potential was found
!
      endif
    enddo
  endif
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
  deallocate(d1BtotiR,stat=status)
  if (status/=0) call deallocate_error('bondorderfd','d1BtotiR')
  deallocate(d1BtotiA,stat=status)
  if (status/=0) call deallocate_error('bondorderfd','d1BtotiA')
  deallocate(d1i,stat=status)
  if (status/=0) call deallocate_error('bondorderfd','d1i')
!
  t2 = g_cpu_time()
  tbondorder = tbondorder + t2 - t1
!
  return
  end subroutine bondorderfd

  subroutine bondorder1fc2(maxlhs,d1cell,matom,vmatom)
!
!  Routine for calculating the first derivatives of bond order potentials
!  for use in a finite difference unphased phonon calculation.
!  This version uses a pre-computed neighbour list.
!
!  NB: Assumes that setedipneigh has been previously called to
!      set up the neighbour lists and perform the screening
!      based on coordination number.
!
!  12/17 Created from bondorder1fc
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!  11/20 Tersoff reorganised
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
!  Julian Gale, CIC, Curtin University, November 2020
!
  use datatypes
  use bondorderdata
  use configurations, only : nregionno, nregiontype, QMMMmode
  use current
  use iochannels
  use neighbours
  use times
  implicit none
!
!  Passed variables
!
  integer(i4),                       intent(in)    :: maxlhs
  integer(i4),                       intent(in)    :: matom
  real(dp),                          intent(inout) :: d1cell(4,maxlhs,*)
  real(dp),                          intent(in)    :: vmatom(4)
!
!  Local variables
!
  integer(i4)                                      :: i
  integer(i4)                                      :: j
  integer(i4)                                      :: k
  integer(i4)                                      :: kk
  integer(i4)                                      :: mA
  integer(i4)                                      :: mR
  integer(i4)                                      :: nati
  integer(i4)                                      :: natj
  integer(i4)                                      :: natk
  integer(i4)                                      :: nboij
  integer(i4)                                      :: nboAij
  integer(i4)                                      :: nboAik
  integer(i4)                                      :: nboOi
  integer(i4)                                      :: nboRij
  integer(i4)                                      :: nboRik
  integer(i4)                                      :: nboZi
  integer(i4)                                      :: ni
  integer(i4)                                      :: nj
  integer(i4)                                      :: njk
  integer(i4)                                      :: nn
  integer(i4)                                      :: npki   
  integer(i4)                                      :: nneighi2
  integer(i4)                                      :: nregioni
  integer(i4)                                      :: nregionj
  integer(i4)                                      :: nregiontypi
  integer(i4)                                      :: nregiontypj
  integer(i4)                                      :: ntypi
  integer(i4)                                      :: ntypj
  integer(i4)                                      :: ntypk
  integer(i4)                                      :: status
  logical                                          :: lfound
  logical                                          :: lneedBOcnderiv
  logical                                          :: lQMMMok
  real(dp)                                         :: bijA
  real(dp)                                         :: bijR
  real(dp)                                         :: bijsumA
  real(dp)                                         :: bijsumA1
  real(dp)                                         :: bijsumAn1
  real(dp)                                         :: bijsumAn2
  real(dp)                                         :: bijsumR
  real(dp)                                         :: bijsumR1
  real(dp)                                         :: bijsumRn1
  real(dp)                                         :: bijsumRn2
  real(dp)                                         :: BOalpAi
  real(dp)                                         :: BOalpRi
  real(dp)                                         :: BOncoAi
  real(dp)                                         :: BOncoRi
  real(dp)                                         :: btotA
  real(dp)                                         :: btotR
  real(dp)                                         :: g_cpu_time
  real(dp),    dimension(:),     allocatable, save :: d1i
  real(dp),    dimension(:),     allocatable, save :: d1BtotiA
  real(dp),    dimension(:),     allocatable, save :: d1BtotiR
  real(dp)                                         :: dexpijkdr
  real(dp)                                         :: dfdr
  real(dp)                                         :: dfikdr
  real(dp)                                         :: dfzdz
  real(dp)                                         :: d2fdr2
  real(dp)                                         :: d2fikdr2
  real(dp)                                         :: d2fzdz2
  real(dp)                                         :: d3fdr3
  real(dp)                                         :: d3fikdr3
  real(dp)                                         :: d3fzdz3
  real(dp)                                         :: dGijkdr(3)
  real(dp)                                         :: d2Gijkdr2(6)
  real(dp)                                         :: d3Gijkdr3(10)
  real(dp)                                         :: dZ
  real(dp)                                         :: deltaZ
  real(dp)                                         :: dedZ
  real(dp)                                         :: expijk
  real(dp)                                         :: f
  real(dp)                                         :: fik
  real(dp)                                         :: fz
  real(dp)                                         :: Gijk
  real(dp)                                         :: omegaik
  real(dp)                                         :: RmA
  real(dp)                                         :: RmR
  real(dp)                                         :: rbijsumA1
  real(dp)                                         :: rbijsumR1
  real(dp)                                         :: rij
  real(dp)                                         :: rik
  real(dp)                                         :: rlambda
  real(dp)                                         :: drlambdadrij
  real(dp)                                         :: drlambdadrik
  real(dp)                                         :: rrij
  real(dp)                                         :: rrik
  real(dp)                                         :: rtmp
  real(dp)                                         :: t1
  real(dp)                                         :: t2
  real(dp)                                         :: Va
  real(dp)                                         :: Vr
  real(dp)                                         :: dVadr
  real(dp)                                         :: dVrdr
  real(dp)                                         :: xdiff
  real(dp)                                         :: ydiff
  real(dp)                                         :: zdiff
  real(dp)                                         :: xji
  real(dp)                                         :: yji
  real(dp)                                         :: zji
  real(dp)                                         :: xki
  real(dp)                                         :: yki
  real(dp)                                         :: zki
  real(dp)                                         :: zAi3
  real(dp)                                         :: zRi3
  real(dp)                                         :: zsign
!
  t1 = g_cpu_time()
!
!  Allocate local memory
!
  allocate(d1i(maxneigh2),stat=status)
  if (status/=0) call outofmemory('bondorder1fc2','d1i')
  allocate(d1BtotiA(maxneigh2),stat=status)
  if (status/=0) call outofmemory('bondorder1fc2','d1BtotiA')
  allocate(d1BtotiR(maxneigh2),stat=status)
  if (status/=0) call outofmemory('bondorder1fc2','d1BtotiR')
!
!  Shift coordinates of matom in neighbour list arrays and correct distances
!
  do i = 1,numat
    if (i.eq.matom) then
      do nn = 1,nneigh(i)
        if (neighno(nn,i).ne.matom) then
          xneigh(nn,i) = xneigh(nn,i) - (vmatom(1) - xclat(i))
          yneigh(nn,i) = yneigh(nn,i) - (vmatom(2) - yclat(i))
          zneigh(nn,i) = zneigh(nn,i) - (vmatom(3) - zclat(i))
          rtmp = xneigh(nn,i)**2 + yneigh(nn,i)**2 + zneigh(nn,i)**2
          rneigh(nn,i) = sqrt(rtmp)
        endif
      enddo
    else
      do nn = 1,nneigh(i)
        j = neighno(nn,i)
        if (j.eq.matom) then
          xneigh(nn,i) = xneigh(nn,i) + (vmatom(1) - xclat(j))
          yneigh(nn,i) = yneigh(nn,i) + (vmatom(2) - yclat(j))
          zneigh(nn,i) = zneigh(nn,i) + (vmatom(3) - zclat(j))
          rtmp = xneigh(nn,i)**2 + yneigh(nn,i)**2 + zneigh(nn,i)**2
          rneigh(nn,i) = sqrt(rtmp)
        endif
      enddo
    endif
  enddo
!
!  Initialise coordination number
!
  Zcn(1:numat) = 0.0_dp
!*****************************
!  Loop over pairs of atoms  *
!*****************************
  do i = 1,numat
!
!  Set variables relating to i
!
    nregioni = nregionno(nsft + nrelf2a(i))
    nregiontypi = nregiontype(nregioni,ncf)
    nati = nat(i)
    ntypi = nftype(i)
!
!  Set total number of distances for neighbours of i
!
    nneighi2 = nneigh(i) + nneigh(i)*(nneigh(i) + 1)/2 
!
!  Initialise derivative storage for neighbours of i
!
    d1i(1:nneighi2) = 0.0_dp
!
!  Find one-body parameters for i
!
    lfound = .false.
    nboOi = 0
    do while (.not.lfound.and.nboOi.lt.nboO)
      nboOi = nboOi + 1
      if (nBOspec0(nboOi).eq.nati) then
        if (nBOtyp0(nboOi).eq.ntypi.or.nBOtyp0(nboOi).eq.0) then
          lfound = .true.
        endif
      endif
    enddo
    if (lfound) then
      BOncoAi = BOncoeffA(nboOi)
      BOalpAi = BOecoeffA(nboOi)
      BOncoRi = BOncoeffR(nboOi)
      BOalpRi = BOecoeffR(nboOi)
    else
      BOncoAi = 1.0_dp
      BOalpAi = 0.0_dp
      BOncoRi = 1.0_dp
      BOalpRi = 0.0_dp
    endif
    zAi3 = BOalpAi**BOncoAi
    zRi3 = BOalpRi**BOncoRi
!
!  Loop over neighbours of i (=> j)
!
    ni = 1
    niloop: do while (ni.le.nneigh(i))
      j = neighno(ni,i)
!
!  Set variables relating to j
!
      nregionj = nregionno(nsft + nrelf2a(j))
      nregiontypj = nregiontype(nregionj,ncf)
      natj = nat(j)
      ntypj = nftype(j)
!
!  Set up i-j quantities
!
      rij = rneigh(ni,i)
      xji = xneigh(ni,i)
      yji = yneigh(ni,i)
      zji = zneigh(ni,i)
      rrij = 1.0_dp/rij
!
!  Find i in neighbour list for j
!
      nj = 1
      lfound = .false.
      do while (nj.lt.nneigh(j).and..not.lfound)
        if (neighno(nj,j).eq.i) then
          xdiff = xneigh(nj,j) + xji
          ydiff = yneigh(nj,j) + yji
          zdiff = zneigh(nj,j) + zji
          lfound = ((abs(xdiff)+abs(ydiff)+abs(zdiff)).lt.1.0d-6)
        endif
        if (.not.lfound) nj = nj + 1
      enddo
!
!  Find two-body bond order potential between i and j
!
      lfound = .false.
      nboij = 0
      do while (.not.lfound.and.nboij.lt.nbopot) 
        nboij = nboij + 1
        if (nBOspec1(nboij).eq.nati.and.nBOspec2(nboij).eq.natj) then
          if ((nBOtyp1(nboij).eq.ntypi.or.nBOtyp1(nboij).eq.0).and. &
              (nBOtyp2(nboij).eq.ntypj.or.nBOtyp2(nboij).eq.0)) then
            lfound = .true.
          endif
        elseif (nBOspec1(nboij).eq.natj.and.nBOspec2(nboij).eq.nati) then
          if ((nBOtyp1(nboij).eq.ntypj.or.nBOtyp1(nboij).eq.0).and. &
              (nBOtyp2(nboij).eq.ntypi.or.nBOtyp2(nboij).eq.0)) then
            lfound = .true.
          endif
        endif
      enddo
      if (.not.lfound) nboij = 0
!
!  Find repulsive bond order potential from j to i
!
      lfound = .false.
      nboRij = 0
      do while (.not.lfound.and.nboRij.lt.nboR) 
        nboRij = nboRij + 1
        if (nBOspecR1(nboRij).eq.nati.and.nBOspecR2(nboRij).eq.natj) then
          if ((nBOtypR1(nboRij).eq.ntypi.or.nBOtypR1(nboRij).eq.0).and. &
              (nBOtypR2(nboRij).eq.ntypj.or.nBOtypR2(nboRij).eq.0)) then
            lfound = .true.
          elseif (nati.eq.natj) then
            if ((nBOtypR2(nboRij).eq.ntypi.or.nBOtypR2(nboRij).eq.0).and. &
                (nBOtypR1(nboRij).eq.ntypj.or.nBOtypR1(nboRij).eq.0)) then
              lfound = .true.
            endif
          endif
        endif
      enddo
      if (.not.lfound) nboRij = 0
!
!  Find attractive bond order potential from j to i
!
      lfound = .false.
      nboAij = 0
      do while (.not.lfound.and.nboAij.lt.nboA)
        nboAij = nboAij + 1
        if (nBOspecA1(nboAij).eq.nati.and.nBOspecA2(nboAij).eq.natj) then
          if ((nBOtypA1(nboAij).eq.ntypi.or.nBOtypA1(nboAij).eq.0).and. &
              (nBOtypA2(nboAij).eq.ntypj.or.nBOtypA2(nboAij).eq.0)) then
            lfound = .true.
          elseif (nati.eq.natj) then
            if ((nBOtypA2(nboAij).eq.ntypi.or.nBOtypA2(nboAij).eq.0).and. &
                (nBOtypA1(nboAij).eq.ntypj.or.nBOtypA1(nboAij).eq.0)) then
              lfound = .true.
            endif
          endif
        endif
      enddo
      if (.not.lfound) nboAij = 0
!
!  QM/MM handling : i & j are both QM atoms => exclude
!
      lQMMMok = .true.
      if (QMMMmode(ncf).gt.0) then
        if (nregiontypi.eq.1.and.nregiontypj.eq.1) lQMMMok = .false.
      endif
      if (nboij.gt.0.and.lQMMMok) then
!****************************************
!  Valid two-body bond order potential  *
!****************************************
!
!  Calculate fij
!
        if (nBOtapertype(nboij).eq.2) then
          call ctaper(rij,rBOmin(nboij),rBOmax(nboij),f,dfdr,d2fdr2,d3fdr3,.true.,.false.,.false.)
        else
          call mdftaper(rij,rBOmin(nboij),rBOmax(nboij),f,dfdr,d2fdr2,d3fdr3,.true.,.false.,.false.)
        endif
!
!  Calculate Bij and Bji - loop over all other neighbours
!
        bijsumA = 0.0_dp
        bijsumR = 0.0_dp
!
        d1BtotiA(1:nneighi2) = 0.0_dp
        d1BtotiR(1:nneighi2) = 0.0_dp
!
!  Loop over neighbours of i .ne. j 
!
        do k = 1,nneigh(i)
          npki = nbopotptr(k,i)
          if (k.ne.ni) then
            rik = rneigh(k,i)
            xki = xneigh(k,i)
            yki = yneigh(k,i)
            zki = zneigh(k,i)
!
!  Repulsive component
!
            if (rik.lt.rBOmax(npki)) then
!
!  Find nboRik
!
              lfound = .false.
              nboRik = 0
              kk = neighno(k,i)
              natk = nat(kk)
              ntypk = nftype(kk)
              do while (.not.lfound.and.nboRik.lt.nboR)
                nboRik = nboRik + 1
                if (nBOspecR1(nboRik).eq.nati.and.nBOspecR2(nboRik).eq.natk) then
                  if ((nBOtypR1(nboRik).eq.ntypi.or.nBOtypR1(nboRik).eq.0).and. &
                      (nBOtypR2(nboRik).eq.ntypk.or.nBOtypR2(nboRik).eq.0)) then
                    lfound = .true.
                  endif
                endif
              enddo
              if (nboRik.gt.0) then
!
!  Calculate fik
!
                if (nBOtapertype(npki).eq.2) then
                  call ctaper(rik,rBOmin(npki),rBOmax(npki),fik,dfikdr,d2fikdr2,d3fikdr3,.true.,.false.,.false.)
                else
                  call mdftaper(rik,rBOmin(npki),rBOmax(npki),fik,dfikdr,d2fikdr2,d3fikdr3,.true.,.false.,.false.)
                endif
! 
!  Set omega for i-k
!
                omegaik = BOocoeffR(nboRik)
                fik = fik*omegaik
                dfikdr = dfikdr*omegaik
!
!  Calculate Gijk
!
                if (nBOtypeR(nboRik).ne.1) then
                  call GthetaBO(xji,yji,zji,xki,yki,zki,nBOtypeR(nboRik),BOccoeffR(1,nboRik), &
                    BOhcoeffR(nboRik),Gijk,dGijkdr,d2Gijkdr2,d3Gijkdr3,.true.,.false.,.false.)
                else
                  Gijk = 1.0_dp
                  dGijkdr = 0.0_dp
                  d2Gijkdr2 = 0.0_dp
                endif
!
!  Calculate exponential factor
!
                mR = nint(BOmcoeffR(nboRik))
                if (nboRij.gt.0) then
                  if (lBOzrlR(nboRij)) then
                    rlambda = (BOlcoeffR(nboRij)*rij - BOlcoeffR(nboRik)*rik)
                  else
                    rlambda = BOlcoeffR(nboRik)*(rij - rik)
                  endif
                else
                  rlambda = BOlcoeffR(nboRik)*(rij - rik)
                endif
                expijk = exp(rlambda**mR)
!
!  Combine terms
!
                bijsumR = bijsumR + Gijk*fik*expijk
!
!  Derivatives
!
!  Find index for j-k 
!
                if (ni.ge.k) then
                  njk = nneigh(i) + ni*(ni-1)/2 + k
                else
                  njk = nneigh(i) + k*(k-1)/2 + ni
                endif
!
                rrik = 1.0_dp/rik
                dfikdr = rrik*dfikdr
                RmR = dble(mR)
                if (mR.ge.1) then
                  if (nboRij.gt.0) then
                    if (lBOzrlR(nboRij)) then
                      drlambdadrij = BOlcoeffR(nboRij)
                      drlambdadrik = BOlcoeffR(nboRik)
                    else
                      drlambdadrij = BOlcoeffR(nboRik)
                      drlambdadrik = BOlcoeffR(nboRik)
                    endif
                  else
                    drlambdadrij = BOlcoeffR(nboRik)
                    drlambdadrik = BOlcoeffR(nboRik)
                  endif
                  dexpijkdr = RmR*(rlambda**(mR-1))*expijk
                else
                  dexpijkdr = 0.0_dp
                endif
!
                d1BtotiR(k)   = d1BtotiR(k)   + Gijk*dfikdr*expijk
!
                d1BtotiR(ni)  = d1BtotiR(ni)  + dGijkdr(1)*fik*expijk
                d1BtotiR(k)   = d1BtotiR(k)   + dGijkdr(2)*fik*expijk
                d1BtotiR(njk) = d1BtotiR(njk) + dGijkdr(3)*fik*expijk
!
                d1BtotiR(ni) = d1BtotiR(ni) + Gijk*fik*dexpijkdr*rrij*drlambdadrij
                d1BtotiR(k)  = d1BtotiR(k)  - Gijk*fik*dexpijkdr*rrik*drlambdadrik
              endif
            endif
!
!  Attractive component
!
            if (rik.lt.rBOmax(npki)) then
!
!  Find nboAik
!
              lfound = .false.
              nboAik = 0
              kk = neighno(k,i)
              natk = nat(kk)
              ntypk = nftype(kk)
              do while (.not.lfound.and.nboAik.lt.nboA)
                nboAik = nboAik + 1
                if (nBOspecA1(nboAik).eq.nati.and.nBOspecA2(nboAik).eq.natk) then
                  if ((nBOtypA1(nboAik).eq.ntypi.or.nBOtypA1(nboAik).eq.0).and. &
                      (nBOtypA2(nboAik).eq.ntypk.or.nBOtypA2(nboAik).eq.0)) then
                    lfound = .true.
                  endif
                endif
              enddo
              if (nboAik.gt.0) then
!
!  Calculate fik
!
                if (nBOtapertype(npki).eq.2) then
                  call ctaper(rik,rBOmin(npki),rBOmax(npki),fik,dfikdr,d2fikdr2,d3fikdr3,.true.,.false.,.false.)
                else
                  call mdftaper(rik,rBOmin(npki),rBOmax(npki),fik,dfikdr,d2fikdr2,d3fikdr3,.true.,.false.,.false.)
                endif
! 
!  Set omega for i-k
!
                omegaik = BOocoeffA(nboAik)
                fik = fik*omegaik
                dfikdr = dfikdr*omegaik
!
!  Calculate Gijk
!
                if (nBOtypeA(nboAik).ne.1) then
                  call GthetaBO(xji,yji,zji,xki,yki,zki,nBOtypeA(nboAik),BOccoeffA(1,nboAik), &
                    BOhcoeffA(nboAik),Gijk,dGijkdr,d2Gijkdr2,d3Gijkdr3,.true.,.false.,.false.)
                else
                  Gijk = 1.0_dp
                  dGijkdr = 0.0_dp
                  d2Gijkdr2 = 0.0_dp
                endif
!
!  Calculate exponential factor
!
                mA = nint(BOmcoeffA(nboAik))
                if (nboAij.gt.0) then
                  if (lBOzrlA(nboAij)) then
                    rlambda = (BOlcoeffA(nboAij)*rij - BOlcoeffA(nboAik)*rik)
                  else
                    rlambda = BOlcoeffA(nboAik)*(rij - rik)
                  endif
                else
                  rlambda = BOlcoeffA(nboAik)*(rij - rik)
                endif
                expijk = exp(rlambda**mA)
!
!  Combine terms
!
                bijsumA = bijsumA + Gijk*fik*expijk
!
!  Derivatives
!
!  Find index for j-k 
!
                if (ni.ge.k) then
                  njk = nneigh(i) + ni*(ni-1)/2 + k
                else
                  njk = nneigh(i) + k*(k-1)/2 + ni
                endif
!
                rrik = 1.0_dp/rik
                dfikdr = rrik*dfikdr
                RmA = dble(mA)
                if (mA.ge.1) then
                  if (nboAij.gt.0) then
                    if (lBOzrlA(nboAij)) then
                      drlambdadrij = BOlcoeffA(nboAij)
                      drlambdadrik = BOlcoeffA(nboAik)
                    else
                      drlambdadrij = BOlcoeffA(nboAik)
                      drlambdadrik = BOlcoeffA(nboAik)
                    endif
                  else
                    drlambdadrij = BOlcoeffA(nboAik)
                    drlambdadrik = BOlcoeffA(nboAik)
                  endif
                  dexpijkdr = RmA*(rlambda**(mA-1))*expijk
                else
                  dexpijkdr = 0.0_dp
                endif
!
                d1BtotiA(k)   = d1BtotiA(k)   + Gijk*dfikdr*expijk
!
                d1BtotiA(ni)  = d1BtotiA(ni)  + dGijkdr(1)*fik*expijk
                d1BtotiA(k)   = d1BtotiA(k)   + dGijkdr(2)*fik*expijk
                d1BtotiA(njk) = d1BtotiA(njk) + dGijkdr(3)*fik*expijk
!
                d1BtotiA(ni) = d1BtotiA(ni) + Gijk*fik*dexpijkdr*rrij*drlambdadrij
                d1BtotiA(k)  = d1BtotiA(k)  - Gijk*fik*dexpijkdr*rrik*drlambdadrik
              endif
            endif
          endif
        enddo
!
!  Raise terms to the power of n, add 1, and then raise to -2*n
!
        if (abs(bijsumA).gt.1.0d-12) then
          bijsumAn2 = bijsumA**(BOncoAi - 2.0_dp)
        else
          bijsumAn2 = 0.0_dp
        endif
        if (abs(bijsumR).gt.1.0d-12) then
          bijsumRn2 = bijsumR**(BOncoRi - 2.0_dp)
        else
          bijsumRn2 = 0.0_dp
        endif
!
        bijsumAn1 = bijsumAn2*bijsumA
        bijsumRn1 = bijsumRn2*bijsumR
!
        bijsumA1 = 1.0_dp + zAi3*bijsumA*bijsumAn1
        bijsumR1 = 1.0_dp + zRi3*bijsumR*bijsumRn1
!
        rbijsumA1 = 1.0_dp/bijsumA1
        rbijsumR1 = 1.0_dp/bijsumR1
!
        bijA = bijsumA1**(-0.5_dp/BOncoAi)
        bijR = bijsumR1**(-0.5_dp/BOncoRi)
!
!  Scale derivatives by bijsum/bjisum factors
!
!
!  First derivatives
!
        if (bijsumA.gt.0.0_dp) then
          rtmp = 0.25_dp*zAi3*bijsumAn1*bijA*rbijsumA1
          do nn = 1,nneighi2
            d1BtotiA(nn) = - rtmp*d1BtotiA(nn)
          enddo
        endif
        if (bijsumR.gt.0.0_dp) then
          rtmp = 0.25_dp*zRi3*bijsumRn1*bijR*rbijsumR1
          do nn = 1,nneighi2
            d1BtotiR(nn) = - rtmp*d1BtotiR(nn)
          enddo
        endif
!
!  Calculate two-body component of potential
!
        Vr = BOacoeff(nboij)*exp(-BOzacoeff(nboij)*rij)
        Va = BObcoeff(nboij)*exp(-BOzbcoeff(nboij)*rij)
        dVrdr = - BOzacoeff(nboij)*Vr
        dVadr = - BOzbcoeff(nboij)*Va
!
        dVrdr = rrij*dVrdr
        dVadr = rrij*dVadr
!
!  Calculate total i-j potential
!
        BtotA = 0.5_dp*bijA
        BtotR = 0.5_dp*bijR
!
!  Add contribution to coordination number
!
        Zcn(i) = Zcn(i) + f*bijA
!
!  Derivatives of Bond Order potential energy
!
        dfdr = rrij*dfdr
        d1i(ni) = d1i(ni) + dfdr*(BtotR*Vr - BtotA*Va)
        d1i(ni) = d1i(ni) + f*(BtotR*dVrdr - BtotA*dVadr)
        do nn = 1,nneighi2
          d1i(nn) = d1i(nn) + f*(Vr*d1BtotiR(nn) - Va*d1BtotiA(nn))
        enddo
      endif
!
!  End of loop over neighbours of i
!
      ni = ni + 1
    enddo niloop
!
!  Add derivatives due to neighbours of i
!
    call d1addfc(i,maxneigh,maxneigh,nneigh,neighno,xneigh,yneigh,zneigh,ineigh,maxlhs,d1cell,matom,nboatomRptr,d1i,.false.)
  enddo
  if (nboZ.gt.0) then
!**************************************
!  Compute coordination contribution  *
!**************************************
    do i = 1,numat
!
!  Set variables relating to i
!
      nati = nat(i)
      ntypi = nftype(i)
!
!  Check for coordination terms for atom i
!
      lfound = .false.
      nboZi = 0
      do while (.not.lfound.and.nboZi.lt.nboZ)
        nboZi = nboZi + 1
        if (nBOspecZ(nboZi).eq.nati) then
          if (nBOtypZ(nboZi).eq.ntypi.or.nBOtypZ(nboZi).eq.0) then
            lfound = .true.
          endif
        endif
      enddo
      if (lfound) then
        dZ = Zcn(i) - BOzcoeffZ(nboZi)
        zsign = sign(1.0_dp,dZ)
        call botaper(dZ,fz,dfzdz,d2fzdz2,d3fzdz3,.true.,.false.,.false.)
        deltaZ = zsign*fz
        lneedBOcnderiv = (abs(dble(nint(deltaZ)) - deltaZ).gt.BOcntol)
!
!  If deltaZ is not an integer then compute derivatives
!
        if (lneedBOcnderiv) then
          dedZ = (BOccoeffZ(1,nboZi) + 2.0_dp*BOccoeffZ(2,nboZi)*deltaZ)*zsign*dfzdz
!
          nregioni = nregionno(nsft + nrelf2a(i))
          nregiontypi = nregiontype(nregioni,ncf)
!
!  Set total number of distances for neighbours of i
!
          nneighi2 = nneigh(i) + nneigh(i)*(nneigh(i) + 1)/2 
!
!  Initialise derivative storage for neighbours of i
!
          d1i(1:nneighi2) = 0.0_dp
!
!  Find one-body parameters for i
!
          lfound = .false.
          nboOi = 0
          do while (.not.lfound.and.nboOi.lt.nboO)
            nboOi = nboOi + 1
            if (nBOspec0(nboOi).eq.nati) then
              if (nBOtyp0(nboOi).eq.ntypi.or.nBOtyp0(nboOi).eq.0) then
                lfound = .true.
              endif
            endif
          enddo
          if (lfound) then
            BOncoAi = BOncoeffA(nboOi)
            BOalpAi = BOecoeffA(nboOi)
          else
            BOncoAi = 1.0_dp
            BOalpAi = 0.0_dp
          endif
          zAi3 = BOalpAi**BOncoAi
!
!  Loop over neighbours of i 
!
          ni = 1
          niloopnc: do while (ni.le.nneigh(i))
            j = neighno(ni,i)
!
!  Set variables relating to j
!
            natj = nat(j)
            ntypj = nftype(j)
            nregionj = nregionno(nsft + nrelf2a(j))
            nregiontypj = nregiontype(nregionj,ncf)
!
!  Set up i-j quantities
!
            rij = rneigh(ni,i)
            xji = xneigh(ni,i)
            yji = yneigh(ni,i)
            zji = zneigh(ni,i)
            rrij = 1.0_dp/rij
!
!  Find i in neighbour list for j
!
            nj = 1
            lfound = .false.
            do while (nj.lt.nneigh(j).and..not.lfound)
              if (neighno(nj,j).eq.i) then
                xdiff = xneigh(nj,j) + xji
                ydiff = yneigh(nj,j) + yji
                zdiff = zneigh(nj,j) + zji
                lfound = ((abs(xdiff)+abs(ydiff)+abs(zdiff)).lt.1.0d-6)
              endif
              if (.not.lfound) nj = nj + 1
            enddo
!
!  Find two-body bond order potential between i and j
!
            lfound = .false.
            nboij = 0
            do while (.not.lfound.and.nboij.lt.nbopot) 
              nboij = nboij + 1
              if (nBOspec1(nboij).eq.nati.and.nBOspec2(nboij).eq.natj) then
                if ((nBOtyp1(nboij).eq.ntypi.or.nBOtyp1(nboij).eq.0).and. &
                    (nBOtyp2(nboij).eq.ntypj.or.nBOtyp2(nboij).eq.0)) then
                  lfound = .true.
                endif
              elseif (nBOspec1(nboij).eq.natj.and.nBOspec2(nboij).eq.nati) then
                if ((nBOtyp1(nboij).eq.ntypj.or.nBOtyp1(nboij).eq.0).and. &
                    (nBOtyp2(nboij).eq.ntypi.or.nBOtyp2(nboij).eq.0)) then
                  lfound = .true.
                endif
              endif
            enddo
            if (.not.lfound) nboij = 0
!
!  Find attractive bond order potential from j to i
!
            lfound = .false.
            nboAij = 0
            do while (.not.lfound.and.nboAij.lt.nboA)
              nboAij = nboAij + 1
              if (nBOspecA1(nboAij).eq.nati.and.nBOspecA2(nboAij).eq.natj) then
                if ((nBOtypA1(nboAij).eq.ntypi.or.nBOtypA1(nboAij).eq.0).and. &
                    (nBOtypA2(nboAij).eq.ntypj.or.nBOtypA2(nboAij).eq.0)) then
                  lfound = .true.
                elseif (nati.eq.natj) then
                  if ((nBOtypA2(nboAij).eq.ntypi.or.nBOtypA2(nboAij).eq.0).and. &
                      (nBOtypA1(nboAij).eq.ntypj.or.nBOtypA1(nboAij).eq.0)) then
                    lfound = .true.
                  endif
                endif
              endif
            enddo
            if (.not.lfound) nboAij = 0
!
!  QM/MM handling : i & j are both QM atoms => exclude
!
            lQMMMok = .true.
            if (QMMMmode(ncf).gt.0) then
              if (nregiontypi.eq.1.and.nregiontypj.eq.1) lQMMMok = .false.
            endif
            if (nboij.gt.0.and.lQMMMok) then
!****************************************
!  Valid two-body bond order potential  *
!****************************************
!
!  Calculate fij
!
              if (nBOtapertype(nboij).eq.2) then
                call ctaper(rij,rBOmin(nboij),rBOmax(nboij),f,dfdr,d2fdr2,d3fdr3,.true.,.false.,.false.)
              else
                call mdftaper(rij,rBOmin(nboij),rBOmax(nboij),f,dfdr,d2fdr2,d3fdr3,.true.,.false.,.false.)
              endif
!
!  Calculate Bij and Bji - loop over all other neighbours
!
              bijsumA = 0.0_dp
              d1BtotiA(1:nneighi2) = 0.0_dp
!
!  Loop over neighbours of i .ne. j 
!
              do k = 1,nneigh(i)
                npki = nbopotptr(k,i)
                if (k.ne.ni) then
                  rik = rneigh(k,i)
                  xki = xneigh(k,i)
                  yki = yneigh(k,i)
                  zki = zneigh(k,i)
!
!  Attractive component
!
                  if (rik.lt.rBOmax(npki)) then
!
!  Find nboAik
!
                    lfound = .false.
                    nboAik = 0
                    kk = neighno(k,i)
                    natk = nat(kk)
                    ntypk = nftype(kk)
                    do while (.not.lfound.and.nboAik.lt.nboA)
                      nboAik = nboAik + 1
                      if (nBOspecA1(nboAik).eq.nati.and.nBOspecA2(nboAik).eq.natk) then
                        if ((nBOtypA1(nboAik).eq.ntypi.or.nBOtypA1(nboAik).eq.0).and. &
                            (nBOtypA2(nboAik).eq.ntypk.or.nBOtypA2(nboAik).eq.0)) then
                          lfound = .true.
                        endif
                      endif
                    enddo
                    if (nboAik.gt.0) then
!
!  Calculate fik
!
                      if (nBOtapertype(npki).eq.2) then
                        call ctaper(rik,rBOmin(npki),rBOmax(npki),fik,dfikdr,d2fikdr2,d3fikdr3,.true.,.false.,.false.)
                      else
                        call mdftaper(rik,rBOmin(npki),rBOmax(npki),fik,dfikdr,d2fikdr2,d3fikdr3,.true.,.false.,.false.)
                      endif
! 
!  Set omega for i-k
!
                      omegaik = BOocoeffA(nboAik)
                      fik = fik*omegaik
                      dfikdr = dfikdr*omegaik
!
!  Calculate Gijk
!
                      if (nBOtypeA(nboAik).ne.1) then
                        call GthetaBO(xji,yji,zji,xki,yki,zki,nBOtypeA(nboAik),BOccoeffA(1,nboAik), &
                                      BOhcoeffA(nboAik),Gijk,dGijkdr,d2Gijkdr2,d3Gijkdr3,.true.,.false.,.false.)
                      else
                        Gijk = 1.0_dp
                        dGijkdr = 0.0_dp
                        d2Gijkdr2 = 0.0_dp
                      endif
!
!  Calculate exponential factor
!
                      mA = nint(BOmcoeffA(nboAik))
                      if (nboAij.gt.0) then
                        if (lBOzrlA(nboAij)) then
                          rlambda = (BOlcoeffA(nboAij)*rij - BOlcoeffA(nboAik)*rik)
                        else
                          rlambda = BOlcoeffA(nboAik)*(rij - rik)
                        endif
                      else
                        rlambda = BOlcoeffA(nboAik)*(rij - rik)
                      endif
                      expijk = exp(rlambda**mA)
!
!  Combine terms
!
                      bijsumA = bijsumA + Gijk*fik*expijk
!
!  Derivatives
!
!  Find index for j-k 
!
                      if (ni.ge.k) then
                        njk = nneigh(i) + ni*(ni-1)/2 + k
                      else
                        njk = nneigh(i) + k*(k-1)/2 + ni
                      endif
!
                      rrik = 1.0_dp/rik
                      dfikdr = rrik*dfikdr
                      RmA = dble(mA)
                      if (mA.ge.1) then
                        if (nboAij.gt.0) then
                          if (lBOzrlA(nboAij)) then
                            drlambdadrij = BOlcoeffA(nboAij)
                            drlambdadrik = BOlcoeffA(nboAik)
                          else
                            drlambdadrij = BOlcoeffA(nboAik)
                            drlambdadrik = BOlcoeffA(nboAik)
                          endif
                        else
                          drlambdadrij = BOlcoeffA(nboAik)
                          drlambdadrik = BOlcoeffA(nboAik)
                        endif
                        dexpijkdr = RmA*(rlambda**(mA-1))*expijk
                      else
                        dexpijkdr = 0.0_dp
                      endif
!
                      d1BtotiA(k)   = d1BtotiA(k)   + Gijk*dfikdr*expijk
!
                      d1BtotiA(ni)  = d1BtotiA(ni)  + dGijkdr(1)*fik*expijk
                      d1BtotiA(k)   = d1BtotiA(k)   + dGijkdr(2)*fik*expijk
                      d1BtotiA(njk) = d1BtotiA(njk) + dGijkdr(3)*fik*expijk
!
                      d1BtotiA(ni) = d1BtotiA(ni) + Gijk*fik*dexpijkdr*rrij*drlambdadrij
                      d1BtotiA(k)  = d1BtotiA(k)  - Gijk*fik*dexpijkdr*rrik*drlambdadrik
                    endif
                  endif
                endif
              enddo
!
              if (abs(bijsumA).gt.1.0d-12) then
                bijsumAn2 = bijsumA**(BOncoAi - 2.0_dp)
              else
                bijsumAn2 = 0.0_dp
              endif
!
              bijsumAn1 = bijsumAn2*bijsumA
!
              bijsumA1 = 1.0_dp + zAi3*bijsumA*bijsumAn1
!
              rbijsumA1 = 1.0_dp/bijsumA1
!
              bijA = bijsumA1**(-0.5_dp/BOncoAi)
!
!  Scale derivatives by bijsum factors
!
!
!  First derivatives
!
              if (bijsumA.gt.0.0_dp) then
                rtmp = 0.25_dp*zAi3*bijsumAn1*bijA*rbijsumA1
                do nn = 1,nneighi2
                  d1BtotiA(nn) = - rtmp*d1BtotiA(nn)
                enddo
              endif
!
!  Derivatives of Bond Order coordination energy
!
              dfdr = rrij*dfdr
              d1i(ni) = d1i(ni) + dedZ*dfdr*bijA
              do nn = 1,nneighi2
                d1i(nn) = d1i(nn) + 2.0_dp*dedZ*f*d1BtotiA(nn)
              enddo
            endif
!
!  End of loop over neighbours of i
!
            ni = ni + 1
          enddo niloopnc
!
!  Add derivatives due to neighbours of i
!
          call d1addfc(i,maxneigh,maxneigh,nneigh,neighno,xneigh,yneigh,zneigh,ineigh, &
                       maxlhs,d1cell,matom,nboatomRptr,d1i,.false.)
!
!  End test over whether derivatives are needed
!
        endif
!
!  End test over whether a coordination potential was found
!
      endif
    enddo
  endif
!
!  Shift back coordinates of matom in neighbour list arrays and correct distances
!
  do i = 1,numat
    if (i.eq.matom) then
      do nn = 1,nneigh(i)
        if (neighno(nn,i).ne.matom) then
          xneigh(nn,i) = xneigh(nn,i) + (vmatom(1) - xclat(i))
          yneigh(nn,i) = yneigh(nn,i) + (vmatom(2) - yclat(i))
          zneigh(nn,i) = zneigh(nn,i) + (vmatom(3) - zclat(i))
          rtmp = xneigh(nn,i)**2 + yneigh(nn,i)**2 + zneigh(nn,i)**2
          rneigh(nn,i) = sqrt(rtmp)
        endif
      enddo
    else
      do nn = 1,nneigh(i)
        j = neighno(nn,i)
        if (j.eq.matom) then
          xneigh(nn,i) = xneigh(nn,i) - (vmatom(1) - xclat(j))
          yneigh(nn,i) = yneigh(nn,i) - (vmatom(2) - yclat(j))
          zneigh(nn,i) = zneigh(nn,i) - (vmatom(3) - zclat(j))
          rtmp = xneigh(nn,i)**2 + yneigh(nn,i)**2 + zneigh(nn,i)**2
          rneigh(nn,i) = sqrt(rtmp)
        endif
      enddo
    endif
  enddo
!
!  Free local memory
!
  deallocate(d1BtotiR,stat=status)
  if (status/=0) call deallocate_error('bondorder1fc2','d1BtotiR')
  deallocate(d1BtotiA,stat=status)
  if (status/=0) call deallocate_error('bondorder1fc2','d1BtotiA')
  deallocate(d1i,stat=status)
  if (status/=0) call deallocate_error('bondorder1fc2','d1i')
!
  t2 = g_cpu_time()
  tbondorder = tbondorder + t2 - t1
!
  return
  end subroutine bondorder1fc2

  subroutine unsetbondorderneigh
!
!  Deallocate neighbour lists for bond order potentials
!
!  12/17 Created from bondordermd.f
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
  use iochannels
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
!  Free memory
!
  deallocate(zneigh,stat=status)
  if (status/=0) call deallocate_error('unsetbondorderneigh','zneigh')
  deallocate(yneigh,stat=status)
  if (status/=0) call deallocate_error('unsetbondorderneigh','yneigh')
  deallocate(xneigh,stat=status)
  if (status/=0) call deallocate_error('unsetbondorderneigh','xneigh')
  deallocate(rneigh,stat=status)
  if (status/=0) call deallocate_error('unsetbondorderneigh','rneigh')
  deallocate(ineigh,stat=status)
  if (status/=0) call deallocate_error('unsetbondorderneigh','ineigh')
  deallocate(nbopotptr,stat=status)
  if (status/=0) call deallocate_error('unsetbondorderneigh','nbopotptr')
  deallocate(neighno,stat=status)
  if (status/=0) call deallocate_error('unsetbondorderneigh','neighno')
  deallocate(Zcn,stat=status)
  if (status/=0) call deallocate_error('unsetbondorderneigh','Zcn')
  deallocate(rBOcutmax,stat=status)
  if (status/=0) call deallocate_error('unsetbondorderneigh','rBOcutmax')
  deallocate(nneigh,stat=status)
  if (status/=0) call deallocate_error('unsetbondorderneigh','nneigh')
  deallocate(nfreeatom,stat=status)
  if (status/=0) call deallocate_error('unsetbondorderneigh','nfreeatom')
  deallocate(latomdone,stat=status)
  if (status/=0) call deallocate_error('unsetbondorderneigh','latomdone')
  deallocate(nboatomRptr,stat=status)
  if (status/=0) call deallocate_error('unsetbondorderneigh','nboatomRptr')
!
  t2 = g_cpu_time()
  tbondorder = tbondorder + t2 - t1
!
  return
  end subroutine unsetbondorderneigh

  end module m_bondorder
