  subroutine bondorder1fc(maxlhs,d1cell,matom,vmatom)
!
!  Routine for calculating the first derivatives of bond order potentials
!  for use in a finite difference unphased phonon calculation.
!
!   1/15 Created from bondorderfc
!   1/15 Typo in variable name corrected
!   1/15 do while loop changed to avoid out of bounds possible issue
!   2/15 Trap for out of bounds corrected by specifying niloop to exit
!   9/15 BOdcoeff replaced by extended BOccoeff array
!  12/15 nBOtapertype added
!   4/16 Algorithm changed slightly so that instead of using lower half
!        triangular, the full neighbour list of each atom is used, but
!        only for half the work.
!   4/16 BO coordination potential added
!   2/18 Trace added
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
  use control,        only : keyword
  use current
  use iochannels
  use neighbours
  use spatialbo,      only : lspatialok => lspatialBOok
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
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
  integer(i4), dimension(:,:,:), allocatable, save :: ineigh         ! Cell indices for vector
  integer(i4)                                      :: itmp
  integer(i4)                                      :: j
  integer(i4)                                      :: k
  integer(i4)                                      :: kk
  integer(i4)                                      :: maxneigh2
  integer(i4)                                      :: maxneigh22
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
  integer(i4), dimension(:,:),   allocatable, save :: nbopotptr
  integer(i4)                                      :: nboatom
  integer(i4), dimension(:),     allocatable, save :: nboatomRptr
  integer(i4)                                      :: nboZi
  integer(i4)                                      :: ni
  integer(i4)                                      :: nj
  integer(i4)                                      :: njk
  integer(i4)                                      :: nmin
  integer(i4)                                      :: nn
  integer(i4)                                      :: nn2
  integer(i4)                                      :: npki   
  integer(i4)                                      :: nptr
  integer(i4), dimension(:,:),   allocatable, save :: neighno
  integer(i4), dimension(:),     allocatable, save :: nneigh
  integer(i4)                                      :: nneighi2
  integer(i4)                                      :: nregioni
  integer(i4)                                      :: nregionj
  integer(i4)                                      :: nregiontypi
  integer(i4)                                      :: nregiontypj
  integer(i4)                                      :: ntypi
  integer(i4)                                      :: ntypj
  integer(i4)                                      :: ntypk
  integer(i4)                                      :: status
  logical,     dimension(:),     allocatable, save :: latomdone
  logical                                          :: lfound
  logical                                          :: lmaxneighok
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
  real(dp),    dimension(:),     allocatable, save :: rBOcutmax
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
  real(dp),    dimension(:,:),   allocatable, save :: rneigh
  real(dp),    dimension(:,:),   allocatable, save :: xneigh
  real(dp),    dimension(:,:),   allocatable, save :: yneigh
  real(dp),    dimension(:,:),   allocatable, save :: zneigh
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
  real(dp),    dimension(:),     allocatable, save :: Zcn
#ifdef TRACE
  call trace_in('bondorder1fc')
#endif
!
  t1 = g_cpu_time()
  allocate(nboatomRptr(numat),stat=status)
  if (status/=0) call outofmemory('bondorder1fc','nboatomRptr')
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
  if (status/=0) call outofmemory('bondorder1fc','latomdone')
  allocate(nneigh(numat),stat=status)
  if (status/=0) call outofmemory('bondorder1fc','nneigh')
  allocate(rBOcutmax(numat),stat=status)
  if (status/=0) call outofmemory('bondorder1fc','rBOcutmax')
  allocate(Zcn(numat),stat=status)
  if (status/=0) call outofmemory('bondorder1fc','Zcn')
!
!  Reinitialisation point should maxneigh be increased
!
100 continue
  lmaxneighok = .true.
  if (allocated(rneigh)) then
    deallocate(d1BtotiR,stat=status)
    if (status/=0) call deallocate_error('bondorder1fc','d1BtotiR')
    deallocate(d1BtotiA,stat=status)
    if (status/=0) call deallocate_error('bondorder1fc','d1BtotiA')
    deallocate(d1i,stat=status)
    if (status/=0) call deallocate_error('bondorder1fc','d1i')
    deallocate(zneigh,stat=status)
    if (status/=0) call deallocate_error('bondorder1fc','zneigh')
    deallocate(yneigh,stat=status)
    if (status/=0) call deallocate_error('bondorder1fc','yneigh')
    deallocate(xneigh,stat=status)
    if (status/=0) call deallocate_error('bondorder1fc','xneigh')
    deallocate(rneigh,stat=status)
    if (status/=0) call deallocate_error('bondorder1fc','rneigh')
    deallocate(ineigh,stat=status)
    if (status/=0) call deallocate_error('bondorder1fc','ineigh')
    deallocate(nbopotptr,stat=status)
    if (status/=0) call deallocate_error('bondorder1fc','nbopotptr')
    deallocate(neighno,stat=status)
    if (status/=0) call deallocate_error('bondorder1fc','neighno')
  endif
!
!  Set parameter for pairwise storage memory
!
  maxneigh2 = maxneigh + maxneigh*(maxneigh + 1)/2
  maxneigh22 = maxneigh2*(maxneigh2 + 1)/2
!
!  Allocate local memory
!
  allocate(neighno(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('bondorder1fc','neighno')
  allocate(nbopotptr(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('bondorder1fc','nbopotptr')
  allocate(ineigh(3_i4,maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('bondorder1fc','ineigh')
  allocate(rneigh(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('bondorder1fc','rneigh')
  allocate(xneigh(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('bondorder1fc','xneigh')
  allocate(yneigh(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('bondorder1fc','yneigh')
  allocate(zneigh(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('bondorder1fc','zneigh')
  allocate(d1i(maxneigh2),stat=status)
  if (status/=0) call outofmemory('bondorder1fc','d1i')
  allocate(d1BtotiA(maxneigh2),stat=status)
  if (status/=0) call outofmemory('bondorder1fc','d1BtotiA')
  allocate(d1BtotiR(maxneigh2),stat=status)
  if (status/=0) call outofmemory('bondorder1fc','d1BtotiR')
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
  call getBOneighbour(maxneigh,rBOcutmax,nBOpotptr,nneigh,neighno,rneigh, &
                      xneigh,yneigh,zneigh,ineigh,latomdone,lmaxneighok)
!**********************************************************************
!  If maxneigh has been exceeded, increase value and return to start  *
!**********************************************************************
  if (.not.lmaxneighok) then
    do i = 1,numat
      if (nneigh(i).gt.maxneigh) maxneigh = nneigh(i)
    enddo
    if (index(keyword,'verb').ne.0) then
      write(ioout,'(/,''  Increasing maxneigh to '',i6)') maxneigh
    endif
    goto 100
  endif
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
!*******************************
!  Sort neighbours into order  *
!******************************* 
  if (lspatialok) then
    do i = 1,numat
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
    enddo
  endif
!
  if (index(keyword,'debu').ne.0) then
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
!  Free local memory
!
  deallocate(d1BtotiR,stat=status)
  if (status/=0) call deallocate_error('bondorder1fc','d1BtotiR')
  deallocate(d1BtotiA,stat=status)
  if (status/=0) call deallocate_error('bondorder1fc','d1BtotiA')
  deallocate(d1i,stat=status)
  if (status/=0) call deallocate_error('bondorder1fc','d1i')
  deallocate(zneigh,stat=status)
  if (status/=0) call deallocate_error('bondorder1fc','zneigh')
  deallocate(yneigh,stat=status)
  if (status/=0) call deallocate_error('bondorder1fc','yneigh')
  deallocate(xneigh,stat=status)
  if (status/=0) call deallocate_error('bondorder1fc','xneigh')
  deallocate(rneigh,stat=status)
  if (status/=0) call deallocate_error('bondorder1fc','rneigh')
  deallocate(ineigh,stat=status)
  if (status/=0) call deallocate_error('bondorder1fc','ineigh')
  deallocate(nbopotptr,stat=status)
  if (status/=0) call deallocate_error('bondorder1fc','nbopotptr')
  deallocate(neighno,stat=status)
  if (status/=0) call deallocate_error('bondorder1fc','neighno')
  deallocate(Zcn,stat=status)
  if (status/=0) call deallocate_error('bondorder1fc','Zcn')
  deallocate(rBOcutmax,stat=status)
  if (status/=0) call deallocate_error('bondorder1fc','rBOcutmax')
  deallocate(nneigh,stat=status)
  if (status/=0) call deallocate_error('bondorder1fc','nneigh')
  deallocate(latomdone,stat=status)
  if (status/=0) call deallocate_error('bondorder1fc','latomdone')
  deallocate(nboatomRptr,stat=status)
  if (status/=0) call deallocate_error('bondorder1fc','nboatomRptr')
!
  t2 = g_cpu_time()
  tbondorder = tbondorder + t2 - t1
#ifdef TRACE
  call trace_out('bondorder1fc')
#endif
!
  return
  end
