  subroutine setatomnodesbo(numat,nprocs,procid,lspatial)
!
!  Sets up the distribution of atoms over nodes for a 
!  spatial decomposition or otherwise in parallel. Version specifically
!  for bond order potentials.
!
!   9/04 Created from setatomnodes
!  11/06 Modified in handling of npgridxptr setting
!   6/07 Structure of arrays for storing distribution changed to 1-D
!  12/07 Unused variables removed
!   7/09 y/z interchange error in parallel spatial decomposition fixed
!   9/15 Trap added for case where size of processor grid matches the
!        number of cells to ensure that cells are evenly distributed
!   8/20 Output of debugging info improved so that ioproc writes everything
!   8/20 Parallel distribution of cells for spatial simplified by using
!        atom number distribution rather than trying to localise in space
!
!  On entry :
!
!  numat         = total number of atoms
!  nprocs        = total number of nodes
!  procid        = local node number
!  lspatial      = if .true. then use a spatial decomposition
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
  use datatypes
  use control,   only : keyword
  use iochannels
  use parallel,  only : ioproc
  use spatialbo, only : lbuffercell => lbuffercellbo
  use spatialbo, only : maxatompernode => maxatompernodebo
  use spatialbo, only : maxcellpernode => maxcellpernodebo
  use spatialbo, only : natomcell => natomcellbo
  use spatialbo, only : natomnodeptr => natomnodeptrbo
  use spatialbo, only : natompernode => natompernodebo
  use spatialbo, only : nbufferx => nbufferxbo
  use spatialbo, only : nbuffery => nbufferybo
  use spatialbo, only : nbufferz => nbufferzbo
  use spatialbo, only : ncellpernode => ncellpernodebo
  use spatialbo, only : ncellnodeptr => ncellnodeptrbo
  use spatialbo, only : nspcell => nspcellbo
  use spatialbo, only : nspcellat => nspcellatbo
  use spatialbo, only : nspcellatptr => nspcellatptrbo
  use spatialbo, only : nspcellat1ptr => nspcellat1ptrbo
  use spatialbo, only : nspmax => nspmaxbo
  use spatialbo, only : nspmin => nspminbo
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                      :: nprocs
  integer(i4), intent(in)                      :: numat
  integer(i4), intent(in)                      :: procid
  logical,     intent(in)                      :: lspatial
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ii
  integer(i4)                                  :: ind
  integer(i4)                                  :: iproc
  integer(i4)                                  :: iprocmin
  integer(i4)                                  :: ix
  integer(i4)                                  :: iy
  integer(i4)                                  :: iz
  integer(i4)                                  :: maxxy
  integer(i4)                                  :: maxx
  integer(i4)                                  :: n
  integer(i4)                                  :: nmax
  integer(i4)                                  :: nmin
  integer(i4), dimension(:), allocatable       :: natompernodes
  integer(i4), dimension(:), allocatable       :: ncellpernodes
  integer(i4)                                  :: nrem
  integer(i4)                                  :: nstep
  integer(i4)                                  :: status
  logical                                      :: lokx
  logical                                      :: loky
  logical                                      :: lokz
!
  if (nprocs.eq.1) then
!**********************
!  Non-parallel case  *
!**********************
    natompernode = numat
    if (natompernode.gt.maxatompernode) then
      maxatompernode = natompernode
      call changemaxatompernodebo
    endif
    do i = 1,natompernode
       natomnodeptr(i) = i
    enddo
    if (lspatial) then
!
!  For spatial algorithm set pointers to cells
!
      ncellpernode = nspcell(1)*nspcell(2)*nspcell(3)
      if (ncellpernode.gt.maxcellpernode) then
        maxcellpernode = ncellpernode
        call changemaxcellpernodebo
      endif
      nspmax(1) = nspcell(1) - nbufferx
      nspmin(1) = nbufferx
      nspmax(2) = nspcell(2) - nbuffery
      nspmin(2) = nbuffery
      nspmax(3) = nspcell(3) - nbufferz
      nspmin(3) = nbufferz
      maxxy = nspcell(1)*nspcell(2)
      maxx  = nspcell(1)
      n = 0
      do iz = 1,nspcell(3)
        lokz = (iz.gt.nbufferz.and.iz.lt.(nspcell(3)-nbufferz+1))
        do iy = 1,nspcell(2)
          loky = (iy.gt.nbuffery.and.iy.lt.(nspcell(2)-nbuffery+1))
          do ix = 1,nspcell(1)
            lokx = (ix.gt.nbufferx.and.ix.lt.(nspcell(1)-nbufferx+1))
            ind = (iz-1)*maxxy + (iy-1)*maxx + ix
            n = n + 1
            ncellnodeptr(n) = ind
            if (lokx.and.loky.and.lokz) then
              lbuffercell(n) = .false.
            else
              lbuffercell(n) = .true.
            endif
          enddo
        enddo
      enddo
    endif
  elseif (.not.lspatial) then
!******************************
!  Non-spatial parallel case  *
!******************************
    nstep = numat/nprocs
    nrem  = numat - nprocs*nstep
    nmax = (procid + 1)*nstep + min(nrem,procid+1)
    nmin = procid*nstep + 1 + min(nrem,procid)
    natompernode = nmax - nmin + 1
    if (natompernode.gt.maxatompernode) then
      maxatompernode = natompernode
      call changemaxatompernodebo
    endif
    do i = 1,natompernode
       natomnodeptr(i) = i + nmin - 1
    enddo
  else
!**************************
!  Spatial parallel case  *
!**************************
!
!  Allocate workspace memory
!
    allocate(natompernodes(nprocs),stat=status)
    if (status/=0) call outofmemory('setatomnodesbo','natompernodes')
    allocate(ncellpernodes(nprocs),stat=status)
    if (status/=0) call outofmemory('setatomnodesbo','ncellpernodes')
!
    maxxy = nspcell(1)*nspcell(2)
    maxx  = nspcell(1)
!
!  Initialise counters
!
    natompernode = 0
    ncellpernode = 0
    natompernodes(1:nprocs) = 0
    ncellpernodes(1:nprocs) = 0
!
!  Loop over cells and assign to nodes based on balancing number of atoms per node
!
    do ix = nbufferx+1,nspcell(1) - nbufferx
      do iy = nbuffery+1,nspcell(2) - nbuffery
        do iz = nbufferz+1,nspcell(3) - nbufferz
          ind = (iz-1)*maxxy + (iy-1)*maxx + ix
!
!  Only include cells that have atoms
!
          if (nspcellat(ind).gt.0) then
!
!  Find node with the least atoms
!
!  Perform minloc operation manually to avoid issues with Intel compiler
!
            iprocmin = 0
            nmin = numat + 1
            do iproc = 1,nprocs
              if (natompernodes(iproc).lt.nmin) then
                nmin = natompernodes(iproc)
                iprocmin = iproc
              endif
            enddo
!
            natompernodes(iprocmin) = natompernodes(iprocmin) + nspcellat(ind)
            ncellpernodes(iprocmin) = ncellpernodes(iprocmin) + 1
!
            if (procid.eq.iprocmin-1) then
!
!  Node specific set up - cell
!
              ncellpernode = ncellpernode + 1
!
              if (ncellpernode.gt.maxcellpernode) then
                maxcellpernode = ncellpernode + 8
                call changemaxcellpernodebo
              endif
!
              ncellnodeptr(ncellpernode) = ind
              lbuffercell(ncellpernode) = .false.
!
!  Node specific set up - atoms
!
              if (natompernode+nspcellat(ind).gt.maxatompernode) then
                maxatompernode = natompernode + nspcellat(ind)
                call changemaxatompernodebo
              endif
              do n = 1,nspcellat(ind)
                natompernode = natompernode + 1
                natomnodeptr(natompernode) = nspcellatptr(nspcellat1ptr(ind)+n)
              enddo
            endif
          endif
        enddo
      enddo
    enddo
!
    if (index(keyword,'debu').ne.0) then
      if (ioproc) then
        do ii = 1,nprocs
          write(ioout,'(''  Cells per Processor = '',2(i8,1x))') ii-1,ncellpernodes(ii)
          write(ioout,'(''  Atoms per Processor = '',2(i8,1x))') ii-1,natompernodes(ii)
        enddo
      endif
    endif
!
!  Free workspace memory
!
    deallocate(ncellpernodes,stat=status)
    if (status/=0) call deallocate_error('setatomnodesbo','ncellpernodes')
    deallocate(natompernodes,stat=status)
    if (status/=0) call deallocate_error('setatomnodesbo','natompernodes')
  endif
!
!  Set pointer from atom to cell containing image in non-buffer region
!
  if (lspatial) then
    do iz = nbufferz+1,nspcell(3) - nbufferz
      do iy = nbuffery+1,nspcell(2) - nbuffery
        do ix = nbufferx+1,nspcell(1) - nbufferx
          ind = (iz-1)*maxxy + (iy-1)*maxx + ix
          do n = 1,nspcellat(ind)
            natomcell(nspcellatptr(nspcellat1ptr(ind)+n)) = ind
          enddo
        enddo
      enddo
    enddo
  endif
!
  return
  end
