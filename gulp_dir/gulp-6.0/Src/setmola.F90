  subroutine setmola
!
!  Sets up molecules using user input molecule information
!
!  Assignment of nmolind is made here for clusters - for bulk
!  the values are set in molind
!  
!  moldim  = dimensionality of molecule (0,1,2,3)
!  moldimi = index to periodic directions of molecule
!
!  Guide to molecule dimensionality indices:
!
!    if moldim = 0: moldimi = 0
!    if moldim = 1: moldimi = 1 => x
!                           = 2 => y
!                           = 3 => z
!    if moldim = 2: moldimi = 1 => xy
!                           = 2 => xz
!                           = 3 => yz
!    if moldim = 3: moldimi = 0
!
!   2/17 Created from setmol
!   2/17 Parallel globalisation of bonding corrected
!   8/19 Added to main branch
!   8/19 Rigid molecules added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!  10/19 nmola2f added
!  11/19 natinmol added
!  12/19 nmolasymptr added
!   4/20 Pointer arrays added for numatnomol
!   4/20 Rigid molecule restarting added
!   4/20 ncorenomol and pointer arrays added
!   6/20 Parallel modifications added
!   6/20 nmolcore changes added
!   7/20 Call to changemaxmolat added
!  12/20 Setup of nmolconnect and pointer added
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
!  Julian Gale, CIC, Curtin University, December 2020
!
  use configurations, only : lbsmat
  use control
  use current
  use element
  use molecule
  use parallel
  use reallocate
  use shells
  use spatial,        only : lspatialok, xinbox, yinbox, zinbox
  use times,          only : tmol
  implicit none
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ibt1
  integer(i4)                                  :: ibt2
  integer(i4)                                  :: ic
  integer(i4)                                  :: ii
  integer(i4)                                  :: iii
  integer(i4)                                  :: iimin
  integer(i4)                                  :: im
  integer(i4)                                  :: im1
  integer(i4)                                  :: ind
  integer(i4)                                  :: indb
  integer(i4)                                  :: ixx
  integer(i4)                                  :: iyy
  integer(i4)                                  :: izz
  integer(i4)                                  :: j
  integer(i4)                                  :: jj
  integer(i4)                                  :: jm
  integer(i4)                                  :: kk
  integer(i4)                                  :: nb1
  integer(i4)                                  :: nb2
  integer(i4)                                  :: nbond
  integer(i4)                                  :: nbondcheck
  integer(i4), dimension(:),     allocatable   :: nexist
  integer(i4)                                  :: ni
  integer(i4)                                  :: nj
  integer(i4)                                  :: ni1
  integer(i4)                                  :: nj1
  integer(i4)                                  :: nm
  integer(i4)                                  :: nmmin
  integer(i4)                                  :: nmp
  integer(i4)                                  :: nti
  integer(i4)                                  :: ntj
  integer(i4)                                  :: nti1
  integer(i4)                                  :: ntj1
  integer(i4), dimension(:),     allocatable   :: ntmp
  integer(i4), dimension(:,:),   allocatable   :: ntmp2
  integer(i4), dimension(:,:,:), allocatable   :: ntmp3
  integer(i4)                                  :: status
  logical                                      :: lbondok
  logical                                      :: lduplicate
  logical                                      :: ljfound
  logical,     dimension(:),     allocatable   :: lmoldone
  logical,     dimension(:),     allocatable   :: lmolunique
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: r2
  real(dp)                                     :: r2min
  real(dp)                                     :: rcut
  real(dp)                                     :: ri
  real(dp)                                     :: rij
  real(dp)                                     :: rj
  real(dp)                                     :: t1
  real(dp)                                     :: t2
  real(dp)                                     :: t3
  real(dp)                                     :: t4
  real(dp)                                     :: xal
  real(dp)                                     :: yal
  real(dp)                                     :: zal
  real(dp)                                     :: xcd
  real(dp)                                     :: ycd
  real(dp)                                     :: zcd
  real(dp)                                     :: xcdi
  real(dp)                                     :: ycdi
  real(dp)                                     :: zcdi
  real(dp)                                     :: xcrd
  real(dp)                                     :: ycrd
  real(dp)                                     :: zcrd
!
  t1 = g_cpu_time()
  lbondok = .true.
!
!  Allocate local memory
!
  allocate(nexist(numat),stat=status)
  if (status/=0) call outofmemory('setmol','nexist')
!
!  Return point for re-initialisation
!
 10 continue
!
!  Initialise arrays
!
  do i = 1,numat
    natmol(i) = 0
    natinmol(i) = 0
    nmolind(i) = 0
    nbonds(i) = 0
    do j = 1,maxbond
      nbonded(j,i) = 0
    enddo
  enddo
  do i = 1,maxmol
    moldim(i) = 0
    moldimi(i) = 0
  enddo
!******************************
! (1) Set up molecule lists   *
!******************************
  nmol = nmolcfg(ncf)
!******************************
! (2) Generate bonding lists  *
!******************************
  if (nprocs.gt.1.and.lfastMPI) then
    nmmin = procid + 1
    nmp = nprocs
  else
    nmmin = 1
    nmp = 1
  endif
!
!  Loop over molecules
!
  do nm = nmmin,nmol,nmp
!
!  Find first atom for this molecule in list
!
    im1 = 0
    do im = 1,nm-1
      im1 = im1 + nmolatomcfg(im,ncf)
    enddo
!
!  Loop over atoms in molecules
!
    do im = 1,nmolatomcfg(nm,ncf)
      i = nmollistcfg(im1+im,ncf)
      if (lspatialok) then
        xal = xinbox(i)
        yal = yinbox(i)
        zal = zinbox(i)
      else
        xal = xclat(i)
        yal = yclat(i)
        zal = zclat(i)
      endif
      ni1 = nat(i)
      nti = nftype(i)
      ni = ni1
      if (ni.gt.maxele) ni = ni - maxele
      ri = rcov(ni)
!
!  Find all atoms bonded to atom i
!
      nbond = 0
      if (nconnect.gt.0) then
        do ic = 1,nconnect
          if (nconnectcfg(ic).eq.ncf) then
            if (n1connect(ic).eq.i) then
              j = n2connect(ic)
!
!  Check that j is in molecule
!
              ljfound = .false.
              jm = 0
              do while (jm.lt.nmolatomcfg(nm,ncf).and..not.ljfound)
                jm = jm + 1
                ljfound = (nmollistcfg(im1+jm,ncf).eq.j)
              enddo
!
              if (.not.ljfound) then
                call outerror('molatom input is inconsistent with connect input',0_i4)
                call stopnow('setmola')
              endif
!
              nbond = nbond + 1
              nbonds(i) = nbonds(i) + 1
              if (nbond.gt.maxbond) then
                maxbond = nbond + 2
                call changemaxbond
              endif
              nbonded(nbond,i) = j
              nbondedtype(1,nbond,i) = nconnecttype(1,ic)
              nbondedtype(2,nbond,i) = nconnecttype(2,ic)
              if (nconnectind(ic).gt.0) then
                nbondind(nbond,i) = nconnectind(ic)
              else
!
!  Find nearest image
!
                if (lspatialok) then
                  xcdi = xinbox(j) - xal
                  ycdi = yinbox(j) - yal
                  zcdi = zinbox(j) - zal
                else
                  xcdi = xclat(j) - xal
                  ycdi = yclat(j) - yal
                  zcdi = zclat(j) - zal
                endif
                r2min = 1000000.0_dp
                iimin = 0
                do iii = 1,iimax
                  xcrd = xcdi + xvec1cell(iii)
                  ycrd = ycdi + yvec1cell(iii)
                  zcrd = zcdi + zvec1cell(iii)
                  r2 = xcrd*xcrd + ycrd*ycrd + zcrd*zcrd
                  if (r2.lt.r2min) then
                    r2min = r2
                    iimin = iii
                  endif
                enddo
                call lintoijk(ii,jj,kk,iimin,imaxl,jmaxl,kmaxl)
                nbondind(nbond,i) = (ii+5)+10*(jj+5)+100*(kk+5)
              endif
!
!  If i = j, then add bond in opposite direction
!
              if (i.eq.j) then
                nbond = nbond + 1
                nbonds(i) = nbonds(i) + 1
                if (nbond.gt.maxbond) then
                  maxbond = nbond + 2
                  call changemaxbond
                endif
                nbonded(nbond,i) = j
                nbondedtype(1,nbond,i) = nconnecttype(1,ic)
                nbondedtype(2,nbond,i) = nconnecttype(2,ic)
                if (nconnectind(ic).gt.0) then
                  nbondind(nbond,i) = 1110 - nconnectind(ic)
                endif
              endif
            elseif (n2connect(ic).eq.i) then
              nbond = nbond + 1
              nbonds(i) = nbonds(i) + 1
              if (nbond.gt.maxbond) then
                maxbond = nbond + 2
                call changemaxbond
              endif
              j = n1connect(ic)
              nbonded(nbond,i) = j
              nbondedtype(1,nbond,i) = nconnecttype(1,ic)
              nbondedtype(2,nbond,i) = nconnecttype(2,ic)
              if (nconnectind(ic).gt.0) then
                nbondind(nbond,i) = 1110 - nconnectind(ic)
              else
!
!  Find nearest image
!
                if (lspatialok) then
                  xcdi = xinbox(j) - xal
                  ycdi = yinbox(j) - yal
                  zcdi = zinbox(j) - zal
                else
                  xcdi = xclat(j) - xal
                  ycdi = yclat(j) - yal
                  zcdi = zclat(j) - zal
                endif
                r2min = 1000000.0_dp
                iimin = 0
                do iii = 1,iimax
                  xcrd = xcdi + xvec1cell(iii)
                  ycrd = ycdi + yvec1cell(iii)
                  zcrd = zcdi + zvec1cell(iii)
                  r2 = xcrd*xcrd + ycrd*ycrd + zcrd*zcrd
                  if (r2.lt.r2min) then
                    r2min = r2
                    iimin = iii
                  endif
                enddo
                call lintoijk(ii,jj,kk,iimin,imaxl,jmaxl,kmaxl)
                nbondind(nbond,i) = (ii+5) + 10*(jj+5) + 100*(kk+5)
              endif
            endif
          endif
        enddo
      endif
      if (.not.lnoautobond.and.ri.ne.0.0_dp) then
        jmloop: do jm = 1,nmolatomcfg(nm,ncf)
!
!  Exclude self term
!
          if (jm.eq.im) cycle jmloop
!
          j = nmollistcfg(im1+jm,ncf)
          if (lspatialok) then
            xcd = xinbox(j)
            ycd = yinbox(j)
            zcd = zinbox(j)
          else
            xcd = xclat(j)
            ycd = yclat(j)
            zcd = zclat(j)
          endif
          nj1 = nat(j)
          ntj = nftype(j)
          nj = nj1
          if (nj.gt.maxele) nj = nj - maxele
          rj = rcov(nj)
!
!  Check whether bond type is excluded
!
          if (nnobo.gt.0) then
            if (ni1.eq.nj1) then
              indb = nj1 + 1000*ni1
              if (nti.lt.ntj) then
                nti1 = nti
                ntj1 = ntj
              else
                nti1 = ntj
                ntj1 = nti
              endif
            elseif (ni1.lt.nj1) then
              indb = nj1 + 1000*ni1
              nti1 = nti
              ntj1 = ntj
            else
              indb = ni1 + 1000*nj1
              nti1 = ntj
              ntj1 = nti
            endif
            lbondok = .true.
            ii = 1
            do while (lbondok.and.(ii.le.nnobo))
              if (indb.eq.nobond(ii)) then
                nb1 = nobotyp(ii)/1000
                nb2 = nobotyp(ii) - 1000*nb1
                if ((nb1.eq.nti1.or.nb1.eq.0).and.(nb2.eq.ntj1.or.nb2.eq.0)) lbondok = .false.
                if (ni1.eq.nj1.and.lbondok) then
                  if ((nb1.eq.ntj1.or.nb1.eq.0).and.(nb2.eq.nti1.or.nb2.eq.0)) lbondok = .false.
                endif
              endif
              ii = ii + 1
            enddo
          endif
!
!  Find default bond types
!
          call setbondtype(ni1,nti,nj1,ntj,ibt1,ibt2)
!
!  Distance check
!
          if (rj.ne.0.0_dp.and.lbondok) then
            rcut = rtol*(ri+rj)
            rcut = rcut*rcut
            if (rj.ne.0.0_dp) then
              xcdi = xcd - xal
              ycdi = ycd - yal
              zcdi = zcd - zal
              do ii = 1,iimax
                xcrd = xcdi + xvec1cell(ii)
                ycrd = ycdi + yvec1cell(ii)
                zcrd = zcdi + zvec1cell(ii)
                rij = xcrd*xcrd + ycrd*ycrd + zcrd*zcrd
                if (rij.le.rcut.and.rij.gt.0.0_dp.and.(i.ne.j.or.ii.ne.iimid)) then
!
!  Valid bond
!
                  nbond = nbond + 1
                  nbonds(i) = nbonds(i) + 1
                  if (nbond.gt.maxbond) then
                    maxbond = nbond + 2
                    call changemaxbond
                  endif
                  nbonded(nbond,i) = j
                  nbondedtype(1,nbond,i) = ibt1
                  nbondedtype(2,nbond,i) = ibt2
                  if (ii.eq.iimid) then
                    nbondind(nbond,i) = 555
                  else
                    ind = ii - 1
                    ixx = (ind/((2*jmaxl+1)*(2*kmaxl+1)))
                    ind = ind - ixx*(2*jmaxl+1)*(2*kmaxl+1)
                    iyy = (ind/(2*kmaxl+1))
                    ind = ind - iyy*(2*kmaxl+1)
                    izz = ind - kmaxl
                    iyy = iyy - jmaxl
                    ixx = ixx - imaxl
                    nbondind(nbond,i) = ixx + 5 + 10*(iyy + 5) + 100*(izz + 5)
                  endif
!
!  Check that bond doesn't duplicate one already specified
!
                  lduplicate = .false.
                  nbondcheck = 1
                  do while (.not.lduplicate.and.nbondcheck.lt.nbond)
                    lduplicate = (nbonded(nbondcheck,i).eq.nbonded(nbond,i).and. &
                                  nbondind(nbondcheck,i).eq.nbondind(nbond,i))
                    nbondcheck = nbondcheck + 1
                  enddo
                  if (lduplicate) then
                    nbond = nbond - 1
                    nbonds(i) = nbonds(i) - 1
                  endif
                endif
              enddo
            endif
          endif
        enddo jmloop
      endif
    enddo
!
!  End of loop over molecules
!
  enddo
  if (nprocs.gt.1.and.lfastMPI) then
!
!  Globalise bonding information
!
!  First find global maxbond value
!
    call imaxall(maxbond,nbond,1_i4,"setmola","maxbond")
!
!  If maxbond has been exceeded redimension and then go back to the beginning
!
    if (nbond.gt.maxbond) then
      maxbond = nbond
!
!  Now update sizes to ensure that all arrays are large enough
!
      call changemaxbond
!
      goto 10
    endif
!
!  Start communicating bonding data
!
    allocate(ntmp(numat),stat=status)
    if (status/=0) call outofmemory('setmola','ntmp')
!
    ntmp(1:numat) = 0
    do nm = nmmin,nmol,nmp
      im1 = 0
      do im = 1,nm-1
        im1 = im1 + nmolatomcfg(im,ncf)
      enddo
      do im = 1,nmolatomcfg(nm,ncf)
        i = nmollistcfg(im1+im,ncf)
        ntmp(i) = nbonds(i)
      enddo
    enddo
    call isumall(ntmp,nbonds,numat,"setmola","nbonds")
!
    deallocate(ntmp,stat=status)
    if (status/=0) call deallocate_error('setmola','ntmp')
    allocate(ntmp2(maxbond,numat),stat=status)
    if (status/=0) call outofmemory('setmola','ntmp2')
!
    ntmp2(1:maxbond,1:numat) = 0
    do nm = nmmin,nmol,nmp
      im1 = 0
      do im = 1,nm-1
        im1 = im1 + nmolatomcfg(im,ncf)
      enddo
      do im = 1,nmolatomcfg(nm,ncf)
        i = nmollistcfg(im1+im,ncf)
        ntmp2(1:nbonds(i),i) = nbonded(1:nbonds(i),i)
      enddo
    enddo
    call isumall(ntmp2,nbonded,maxbond*numat,"setmola","nbonded")
!
    ntmp2(1:maxbond,1:numat) = 0
    do nm = nmmin,nmol,nmp
      im1 = 0
      do im = 1,nm-1
        im1 = im1 + nmolatomcfg(im,ncf)
      enddo
      do im = 1,nmolatomcfg(nm,ncf)
        i = nmollistcfg(im1+im,ncf)
        ntmp2(1:nbonds(i),i) = nbondind(1:nbonds(i),i)
      enddo
    enddo
    call isumall(ntmp2,nbondind,maxbond*numat,"setmola","nbondind")
!
    deallocate(ntmp2,stat=status)
    if (status/=0) call deallocate_error('setmola','ntmp2')
    allocate(ntmp3(2,maxbond,numat),stat=status)
    if (status/=0) call outofmemory('setmola','ntmp3')
!
    ntmp3(1:2,1:maxbond,1:numat) = 0
    do nm = nmmin,nmol,nmp
      im1 = 0
      do im = 1,nm-1
        im1 = im1 + nmolatomcfg(im,ncf)
      enddo
      do im = 1,nmolatomcfg(nm,ncf)
        i = nmollistcfg(im1+im,ncf)
        ntmp3(1:2,1:nbonds(i),i) = nbondedtype(1:2,1:nbonds(i),i)
      enddo
    enddo
    call isumall(ntmp3,nbondedtype,2_i4*maxbond*numat,"setmola","nbondedtype")
!
    deallocate(ntmp3,stat=status)
    if (status/=0) call deallocate_error('setmola','ntmp3')
  endif
!***************************
!  Set up molecule values  *
!***************************
  im1 = 0
  do nm = 1,nmol
!
!  Loop over atoms in molecules
!
    nmolatom(nm) = nmolatomcfg(nm,ncf)
    nmolcore(nm) = nmolcorecfg(nm,ncf)
    if (nmolatom(nm).gt.maxmolat) then
      maxmolat = nmolatom(nm) + 10_i4
      call changemaxmolat
    endif
    nmolptr(nm) = im1
    do im = 1,nmolatomcfg(nm,ncf)
      i = nmollistcfg(im1+im,ncf)
      natmol(i) = nm
      natinmol(i) = im
      nmolind(i) = 555
      nmollist(im1+im) = nmollistcfg(im1+im,ncf)
    enddo
!
!  Rigid molecule restarting
!
    if (lmolrigid) then
      molQ(1:3,nm) = molQcfg(1:3,nm,ncf)
      molcom(1:3,nm) = molcomcfg(1:3,nm,ncf)
      do im = 1,nmolcorecfg(nm,ncf)
        molQxyz(1:3,im,nm) = molQxyzcfg(1:3,im1+im,ncf)
      enddo
    endif
    im1 = im1 + nmolatomcfg(nm,ncf)
  enddo
!***************************************************
!  Build pointers from molecules to nconnect list  *
!***************************************************
  if (nconnect.gt.0) then
!
!  Initialise number of connections per molecule
!
    nmolconnect(1:nmol) = 0
!
!  Loop over connections and assign to molecules in this configuration
!
    do ic = 1,nconnect
      if (nconnectcfg(ic).eq.ncf) then
        nm = natmol(n1connect(ic))
        nmolconnect(nm) = nmolconnect(nm) + 1
        if (nmolconnect(nm).gt.maxconnectpermol) then
          maxconnectpermol = nmolconnect(nm) + 10
          call changemaxconnectpermol
        endif
        nmolconnectptr(nmolconnect(nm),nm) = ic
      endif
    enddo
  endif
!**************************************
!  Obtain cell indices for each atom  *
!**************************************
  if (ndim.ne.0) then
    t3 = g_cpu_time()
    call molind
    t4 = g_cpu_time()
    tmol = tmol - t4 + t3
  endif
!
!  Check that coulomb subtraction is not present with periodic
!  molecules as this would cause an error.
!
  if (index(keyword,'mole').ne.0) then
    do i = 1,nmol
      if (moldim(i).gt.0) then
        call outerror('Coulomb subtraction within periodic molecules is not allowed',0_i4)
        call stopnow('setmola')
      endif
    enddo
  endif
!
!  Check that the dimensionality of the molecule doesn't
!  exceed that of the system which would clearly be an error!
!
  do i = 1,nmol
    if (moldim(i).gt.ndim) then
      call outerror('dimensionality of molecule is too high',0_i4)
      call stopnow('setmola')
    endif
  enddo
!
!  Free local memory
!
  deallocate(nexist,stat=status)
  if (status/=0) call deallocate_error('setmola','nexist')
  if (lrigid) then
!
!  Set number of atoms not in a rigid molecule in the asymmetric unit
!
    nasymnomol = 0
    do i = 1,nasym
      if (natmol(nrela2f(i)).eq.0.or.iatn(i).gt.maxele) then
        nasymnomol = nasymnomol + 1
        nasymnomolptr(nasymnomol) = i
        nasymnomolrptr(i) = nasymnomol
      else
        nasymnomolrptr(i) = 0
      endif
    enddo
!
!  Set number of atoms not in a rigid molecule in the full cell
!
    numatnomol = 0
    nbsmatnomol = 0
    do i = 1,numat
      if (natmol(i).eq.0.or.nat(i).gt.maxele) then
        numatnomol = numatnomol + 1
        numatnomolptr(numatnomol) = i
        numatnomolrptr(i) = numatnomol
        if (lbsmat(nrelf2a(i)+nsft)) nbsmatnomol = nbsmatnomol + 1
      else
        numatnomolrptr(i) = 0
      endif
    enddo
!
!  Set number of cores not in a rigid molecule in the full cell
!
    ncorenomol = 0
    do ii = 1,ncore
      i = ncoptr(ii)
      if (natmol(i).eq.0) then
        ncorenomol = ncorenomol + 1
        ncorenomolptr(ncorenomol) = i
        ncorenomolrptr(i) = ncorenomol
      else
        ncorenomolrptr(i) = 0
      endif
    enddo
!
!  Parallel counters for rigid molecules
!
    if (nprocs.gt.1) then
      numatnomolonnode = 0
      do i = 1,numatnomol
        if (atom2node(numatnomolptr(i)).eq.procid) then
          numatnomolonnode = numatnomolonnode + 1
          numatnomolonnodeptr(numatnomolonnode) = numatnomolptr(i)
        endif
      enddo
    else
      numatnomolonnode = numatnomol
      numatnomolonnodeptr(1:numatnomol) = numatnomolptr(1:numatnomol)
    endif
!
!  Check that numatnomol is sensible
!
    if (numatnomol.lt.0) then
      call outerror('numatnomol below zero',0_i4)
      call stopnow('setmola')
    endif
!
!  Find number of symmetry unique molecules in cell
!
    allocate(lmoldone(nmol),stat=status)
    if (status/=0) call outofmemory('setmol','lmoldone')
    allocate(lmolunique(nmol),stat=status)
    if (status/=0) call outofmemory('setmol','lmolunique')
!
    lmoldone(1:nmol) = .false.
    lmolunique(1:nmol) = .true.
!
    nmolasym = 0
    do i = 1,nasym
      ii = nrela2f(i)
      if (natmol(ii).ne.0) then
        if (lmolunique(natmol(ii)).and..not.lmoldone(natmol(ii))) then
          nmolasym = nmolasym + 1
          nmolasymno(nmolasym) = 1
          nmolasymptr(1,nmolasym) = i
          nmola2f(nmolasym) = natmol(ii)
          nmolf2a(natmol(ii)) = nmolasym
          lmoldone(natmol(ii)) = .true.
        elseif (lmolunique(natmol(ii))) then
          nmolasymno(nmolf2a(natmol(ii))) = nmolasymno(nmolf2a(natmol(ii))) + 1
          nmolasymptr(nmolasymno(nmolf2a(natmol(ii))),nmolf2a(natmol(ii))) = i
        endif
!
!  Loop over atoms looking for those equivalent in a different molecule
!
        do j = 1,numat
          if (nrelf2a(j).eq.i.and.ii.ne.j) then
            if (natmol(j).ne.natmol(ii)) then
              lmolunique(natmol(j)) = .false.
              nmolf2a(natmol(j)) = nmolf2a(natmol(ii))
            endif
          endif
        enddo
      endif
    enddo
!
!  Count the number of molecules left
!
    nmolasym = 0
    do i = 1,nmol
      if (lmolunique(i)) then
        nmolasym = nmolasym + 1
        nmola2f(nmolasym) = i
      endif
    enddo
!
!  Find number of symmetry equivalent molecules 
!
    nmoleqv(1:nmolasym) = 0
    do i = 1,nmol
      nmoleqv(nmolf2a(i)) = nmoleqv(nmolf2a(i)) + 1
      if (nmoleqv(nmolf2a(i)).gt.maxmoleqv) then
        maxmoleqv = nmoleqv(nmolf2a(i)) + 1
        call changemaxmoleqv
      endif
      nmolasymeqvptr(nmoleqv(nmolf2a(i)),nmolf2a(i)) = i
    enddo
!
    deallocate(lmolunique,stat=status)
    if (status/=0) call deallocate_error('setmol','lmolunique')
    deallocate(lmoldone,stat=status)
    if (status/=0) call deallocate_error('setmol','lmoldone')
!
!  For rigid molecules, perform extra set up
!
    call setrigidmol
  else
    nasymnomol = nasym
    do i = 1,nasymnomol
      nasymnomolptr(i) = i
      nasymnomolrptr(i) = i
    enddo
    nmolasym = nmol
    do i = 1,nmol
      nmola2f(i) = i
      nmolf2a(i) = i
      nmoleqv(i) = 1
      nmolasymeqvptr(1,i) = i
    enddo
!
    numatnomol = numat
    nbsmatnomol = nbsmat
    ncorenomol = ncore
    do i = 1,numat
      numatnomolptr(i) = i
      numatnomolrptr(i) = i
    enddo
    do i = 1,ncore
      ncorenomolptr(i) = i
      ncorenomolrptr(i) = i
    enddo
  endif
!
  t2 = g_cpu_time()
  tmol = tmol + t2 - t1
!
  return
  end
