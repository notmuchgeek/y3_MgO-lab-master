  subroutine seteembond(lall)
!
!  Tries to determine the optimal set of bond for EEM split bond
!  algorithm such that there are no redundancies. If input flag
!  is true, then all bonds are found regardless of some degree of
!  redundancy.
!
!  Approach is as follows:
!  1) Find unique set of bonds amoungst bonds of atoms in molecule
!  2) For periodic molecules only include the bonds in the cell
!     and neglect those across boundaries
!
!   6/18 Created
!   8/19 Correction to setting of mi/mj 
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
  use configurations, only : nregionno
  use control
  use current
  use eembonds
  use element
  use iochannels
  use molecule
  use parallel
  implicit none
!
!  Passed variables
!
  logical,                         intent(in)  :: lall
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: j
  integer(i4)                                  :: m
  integer(i4)                                  :: mc
  integer(i4)                                  :: mi
  integer(i4)                                  :: mj
  integer(i4)                                  :: nadone
  integer(i4)                                  :: nb
  integer(i4)                                  :: nb2
  integer(i4)                                  :: nbc
  integer(i4)                                  :: nbdone
  integer(i4)                                  :: nbi
  integer(i4)                                  :: nbj
  integer(i4)                                  :: nbok
  integer(i4)                                  :: nc
  integer(i4)                                  :: nconnat
  integer(i4)                                  :: nconnatnow
  integer(i4)                                  :: neb
  integer(i4)                                  :: nebmax
  integer(i4)                                  :: nebtot
  integer(i4), dimension(:),     allocatable   :: nbokptr
  integer(i4), dimension(:),     allocatable   :: nconnatptr
  integer(i4), dimension(:),     allocatable   :: nconnatrptr
  integer(i4), dimension(:,:),   allocatable   :: nebonded
  integer(i4), dimension(:),     allocatable   :: nmolrlist
  integer(i4)                                  :: status
  logical,     dimension(:),     allocatable   :: ladone
  logical,     dimension(:),     allocatable   :: lbadded
  logical,     dimension(:),     allocatable   :: lbdone
  logical,     dimension(:),     allocatable   :: lbok
  logical                                      :: lamisave
  logical                                      :: lamjsave
  logical                                      :: lconnected
  logical                                      :: lfinished
  logical                                      :: lfound
  logical                                      :: lloop
  logical                                      :: lone
  logical                                      :: ltwo
!
!  Initialise neembond
!
  neembond = 0
!
!  Allocate workspace 
!
  allocate(nmolrlist(numat),stat=status)
  if (status/=0) call outofmemory('seteembond','nmolrlist')
!
!  Set up reverse pointer to nmolrlist : atom to atom number in its molecule
!
  nmolrlist(1:numat) = 0
  do m = 1,nmol
    do mi = 1,nmolatom(m)
      i = nmollist(nmolptr(m)+mi)
      nmolrlist(i) = mi
    enddo
  enddo
!
!  Loop over molecules in the system
!
  mloop: do m = 1,nmol
!
!  Find total number of bonds for molecule
!
    nebmax = 0
    do mi = 1,nmolatom(m)
      i = nmollist(nmolptr(m)+mi)
      nebmax = nebmax + nbonds(i)
    enddo
    nebmax = nebmax/2
!
!  Allocate workspace for this molecule
!
    allocate(nebonded(4,nebmax),stat=status)
    if (status/=0) call outofmemory('seteembond','nebonded')
    allocate(nbokptr(nebmax),stat=status)
    if (status/=0) call outofmemory('seteembond','nbokptr')
    allocate(lbok(nebmax),stat=status)
    if (status/=0) call outofmemory('seteembond','lbok')
!
    nbokptr(1:nebmax) = 0
!
!  Build a list of potential bonds for the system excluding duplicates for
!  the same pair and the same bond from two directions
!
    nebtot = 0
    do mi = 1,nmolatom(m)
      i = nmollist(nmolptr(m)+mi)
      if (lelementOK(nat(i)).and.nregionno(nsft+nrelf2a(i)).eq.1) then
        do nbi = 1,nbonds(i)
          j = nbonded(nbi,i)
          mj = nmolrlist(j)
          if (j.lt.i.and.lelementOK(nat(j)).and.nregionno(nsft+nrelf2a(j)).eq.1) then
!
!  Check whether atom pair already has a bond
!
            lfound = .false.
            nbj = 0
            do while (nbj.lt.nebtot.and..not.lfound)
              nbj = nbj + 1
              lfound = ((i.eq.nebonded(1,nbj).and.j.eq.nebonded(2,nbj)).or. &
                        (j.eq.nebonded(1,nbj).and.i.eq.nebonded(2,nbj)))
            enddo
            if (.not.lfound) then
              lfound = .false.
              nbj = 0
              do while (nbj.lt.nbonds(j).and..not.lfound)
                nbj = nbj + 1
                lfound = (i.eq.nbonded(nbj,j))
              enddo
!
              nebtot = nebtot + 1
              nebonded(1,nebtot) = mi
              nebonded(2,nebtot) = mj
              nebonded(3,nebtot) = nbi
              nebonded(4,nebtot) = nbj
            endif
          endif
        enddo
      endif
    enddo
!
!  Initialise flag that indicates bonds should be accepted
!
    lbok(1:nebtot) = .true.
!
!  If number of potential bonds is greater than or equal to the number of atoms 
!  then need to decide which ones to remove
!
    if (nebtot.ge.nmolatom(m).and..not.lall) then
      allocate(ladone(nmolatom(m)),stat=status)
      if (status/=0) call outofmemory('seteembond','ladone')
      allocate(lbdone(nebtot),stat=status)
      if (status/=0) call outofmemory('seteembond','lbdone')
      allocate(lbadded(nebtot),stat=status)
      if (status/=0) call outofmemory('seteembond','lbadded')
      allocate(nconnatptr(nmolatom(m)),stat=status)
      if (status/=0) call outofmemory('seteembond','nconnatptr')
      allocate(nconnatrptr(nmolatom(m)),stat=status)
      if (status/=0) call outofmemory('seteembond','nconnatrptr')
!
      ladone(1:nmolatom(m)) = .false.
      lbdone(1:nebtot) = .false.
      lbok(1:nebtot) = .false.
      nadone = 0
      nbdone = 0
      nbok = 0
!
!  Iterative loop trying to find bonds to connect all atoms
!
      lconnected = .false.      ! Flag that indicates whether there is a path that connects all atoms
      ltwo = .true.             ! Flag that indicates whether there are bonds with 2 atoms yet to be done
      lone = .true.             ! Flag that indicates whether there are bonds with 1 atom yet to be done
      bsearch: do while (nbok.lt.(nmolatom(m)-1))
!
!  Try to find a bond that connects 2 atoms that have not been joined yet
!
        if (ltwo) then
          nb = 0
          lfound = .false.
          do while (nb.lt.nebtot.and..not.lfound)
            nb = nb + 1
            lfound = (.not.lbdone(nb).and.(.not.ladone(nebonded(1,nb)).and..not.ladone(nebonded(2,nb))))
          enddo
          if (.not.lfound) ltwo = .false.
        endif
!
!  If there are no bond that connects 2 atoms that have not been joined then look for 1
!
        if (lone.and..not.ltwo) then
          nb = 0
          lfound = .false.
          do while (nb.lt.nebtot.and..not.lfound)
            nb = nb + 1
            lfound = (.not.lbdone(nb).and.(.not.ladone(nebonded(1,nb)).or..not.ladone(nebonded(2,nb))))
          enddo
          if (.not.lfound) lone = .false.
        endif
!
!  If all atoms have been done try first bond that hasn't been done
!
        if (.not.lone.and..not.ltwo) then
          nb = 0
          lfound = .false.
          do while (nb.lt.nebtot.and..not.lfound)
            nb = nb + 1
            lfound = (.not.lbdone(nb))
          enddo
          if (.not.lfound) then
            call outerror('error in bond search algorithm for EEM',0_i4)
            call stopnow('seteembond')
          endif
        endif
!
!  Set variables for selected bond
!
        mi = nebonded(1,nb)
        mj = nebonded(2,nb)
        i = nmollist(nmolptr(m)+mi)
        j = nmollist(nmolptr(m)+mj)
!
!  Flag atoms and bonds that are being done
!
        lamisave = ladone(mi)
        lamjsave = ladone(mj)
        ladone(mi) = .true.
        ladone(mj) = .true.
        lbdone(nb) = .true.
!
!  Increment number of accepted bonds and set pointer
!
        nbok = nbok + 1
        nbokptr(nbok) = nb
        lbok(nb) = .true.
        if (nbok.gt.1) then
!
!  Find out whether including this bond would create a loop or finish the connection of all atoms
!
          lloop = .false.
          nconnat = 2
          nconnatnow = 2
          nconnatptr(1) = mi
          nconnatptr(2) = mj
          nconnatrptr(1:nmolatom(m)) = 0
          nconnatrptr(mi) = 1
          nconnatrptr(mj) = 2
          lfinished = .false.
          lbadded(1:nbok-1) = .false.
          bloop: do while (.not.lfinished)
!
!  Loop over bonds looking for connected atoms to add 
!
            do nc = 1,nconnatnow
              mc = nconnatptr(nc)
              do nb2 = 1,nbok-1
                if (.not.lbadded(nb2)) then
                  nbc = nbokptr(nb2)
                  if (nebonded(1,nbc).eq.mc) then
                    lbadded(nb2) = .true.
                    if (nebonded(2,nbc).eq.mi.or.nebonded(2,nbc).eq.mj) lloop = .true.
                    if (lloop) exit bloop
                    nconnat = nconnat + 1
                    nconnatptr(nconnat) = nebonded(2,nbc)
                    nconnatrptr(nebonded(2,nbc)) = nconnat
                  elseif (nebonded(2,nbc).eq.mc) then
                    lbadded(nb2) = .true.
                    if (nebonded(1,nbc).eq.mi.or.nebonded(1,nbc).eq.mj) lloop = .true.
                    if (lloop) exit bloop
                    nconnat = nconnat + 1
                    nconnatptr(nconnat) = nebonded(1,nbc)
                    nconnatrptr(nebonded(1,nbc)) = nconnat
                  endif
                endif
              enddo
            enddo
!
!  If number of nconnected atoms hasn't increased then this search is finished
!
            lfinished = (nconnat.eq.nconnatnow)
!
!  Update counters for next pass
!
            nconnatnow = nconnat
          enddo bloop
!
          if (lloop) then
!
!  Remove bond from the list
!
            nbok = nbok - 1
            lbok(nb) = .false.
            ladone(mi) = lamisave
            ladone(mj) = lamjsave
            cycle bsearch
          endif
!
!  Check whether all atoms are connected
!
          lconnected = (nconnat.eq.nmolatom(m))
        else
!
!  If nbok is 1 then test to finish is that the number of atoms in the molecule is 2
!
          lconnected = (nmolatom(m).eq.2)
        endif
!
!  If connectivity is complete exit the search for bonds
!
        if (lconnected) exit bsearch
      enddo bsearch
!
      deallocate(nconnatrptr,stat=status)
      if (status/=0) call deallocate_error('seteembond','nconnatrptr')
      deallocate(nconnatptr,stat=status)
      if (status/=0) call deallocate_error('seteembond','nconnatptr')
      deallocate(lbadded,stat=status)
      if (status/=0) call deallocate_error('seteembond','lbadded')
      deallocate(lbdone,stat=status)
      if (status/=0) call deallocate_error('seteembond','lbdone')
      deallocate(ladone,stat=status)
      if (status/=0) call deallocate_error('seteembond','ladone')
    endif
!
!  Add bonds to global arrays and finish set up tasks
!
    do neb = 1,nebtot
      if (lbok(neb)) then
        neembond = neembond + 1
        if (neembond.gt.maxeembond) then
          maxeembond = maxeembond + 40
          call changemaxeembond
        endif
        neembonded(1:2,neembond) = nebonded(1:2,neb)
!
!  Set pointers from bonds to neembond
!
        mi = nebonded(1,neb)
        mj = nebonded(2,neb)
        i = nmollist(nmolptr(m)+mi)
        j = nmollist(nmolptr(m)+mj)
        nbi = nebonded(3,neb)
        nbj = nebonded(4,neb)
        nbondqb(nbi,i) = neembond
        nbondqb(nbj,j) = - neembond
      endif
    enddo
!
!  Deallocate workspace for this molecule
!
    deallocate(lbok,stat=status)
    if (status/=0) call deallocate_error('seteembond','lbok')
    deallocate(nbokptr,stat=status)
    if (status/=0) call deallocate_error('seteembond','nbokptr')
    deallocate(nebonded,stat=status)
    if (status/=0) call deallocate_error('seteembond','nebonded')
!
!  End of loop over molecules in the system
!
  enddo mloop
!
!  Deallocate workspace 
!
  deallocate(nmolrlist,stat=status)
  if (status/=0) call deallocate_error('seteembond','nmolrlist')
!
!  Debugging info
!
  if (ldebug.and.ioproc) then
    write(ioout,'(/,''  Bonds for EEM split charges: '',/)')
    do m = 1,neembond
      write(ioout,'(i6,'' : '',i7,1x,i7)') m,(neembonded(j,m),j=1,2)
    enddo
  endif
!
  return
  end
