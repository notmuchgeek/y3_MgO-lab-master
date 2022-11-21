  module kim_functions
!
!  Module that contains the data structures and associated KIM functions within GULP
!  including the neighbour lists.
!
!   7/16 Modified for OpenKIM version 1.7.3
!   2/18 Trace added
!   8/18 Modified for version 2 of OpenKIM
!   8/18 Variables for neighbour list renamed to avoid conflicts
!
!  Julian Gale, Curtin University, August 2018
!
    use datatypes
#ifdef TRACE
    use trace,         only : trace_in, trace_out
#endif
    integer(i4),                            private, save :: maxneighk = 12
    integer(i4), dimension(:),     pointer, private, save :: nneighk => null()
    integer(i4), dimension(:,:),   pointer, private, save :: neighkno => null()
    integer(i4), dimension(:,:,:), pointer, private, save :: neighkcell => null()
    real(dp),    dimension(:,:),   pointer, private, save :: rneighk => null()
    real(dp),    dimension(:,:,:), pointer, private, save :: xyzneighk => null()

  contains

  subroutine set_kim_cellimages
!
!  Store linear array of lattice vectors for cell images that are 
!  needed for building the OpenKIM neighbour list
!
!   8/18 Created from rlist
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
!  Copyright Curtin University 2018
!
!  Julian Gale, CIC, Curtin University, August 2018
!
  use current
  use kim_cellimages
  use kim_models,    only : kim_influence, kim_cutoff
  use reallocate
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Local variables
!
  integer(i4)        :: ierror
  integer(i4)        :: ii
  integer(i4)        :: iim
  integer(i4)        :: jj
  integer(i4)        :: kk
  real(dp)           :: xcdi
  real(dp)           :: ycdi
  real(dp)           :: zcdi
  real(dp)           :: xcdj
  real(dp)           :: ycdj
  real(dp)           :: zcdj
  real(dp)           :: xcrd
  real(dp)           :: ycrd
  real(dp)           :: zcrd
!
#ifdef TRACE
  call trace_in('set_kim_cellimages')
#endif
!
!  Compute number of cell images required in each direction for cutoff
!
  kim_ncells(1:3) = 0
  if (ndim.eq.3) then
    call uncell3D(rv,a,b,c,alpha,beta,gamma)
    kim_ncells(1) = (kim_cutoff/a) + 1
    kim_ncells(2) = (kim_cutoff/b) + 1
    kim_ncells(3) = (kim_cutoff/c) + 1
  elseif (ndim.eq.2) then
    call uncell2D(rv,a,b,alpha)
    kim_ncells(1) = (kim_cutoff/a) + 1
    kim_ncells(2) = (kim_cutoff/b) + 1
  elseif (ndim.eq.1) then
    call uncell1D(rv,a)
    kim_ncells(1) = (kim_cutoff/a) + 1
  endif
!
!  Compute number of cell images required in each direction for sphere of influence
!
  kim_nicells(1:3) = 0
  if (ndim.eq.3) then
    kim_nicells(1) = (kim_influence/a) + 1
    kim_nicells(2) = (kim_influence/b) + 1
    kim_nicells(3) = (kim_influence/c) + 1
  elseif (ndim.eq.2) then
    kim_nicells(1) = (kim_influence/a) + 1
    kim_nicells(2) = (kim_influence/b) + 1
  elseif (ndim.eq.1) then
    kim_nicells(1) = (kim_influence/a) + 1
  endif
!
!  Set the total number of cells and check memory for cutoff cell vectors
!
  kim_nicell = (2*kim_nicells(1) + 1)*(2*kim_nicells(2) + 1)*(2*kim_nicells(3) + 1)
  kim_ncell  = (2*kim_ncells(1)  + 1)*(2*kim_ncells(2)  + 1)*(2*kim_ncells(3)  + 1)
  if (kim_ncell.gt.kim_maxncell) then
    kim_maxncell = kim_ncell
    call realloc(kim_cellvector,3_i4,kim_maxncell,ierror)
    if (ierror.ne.0) call outofmemory('set_kim_cellimages','kim_cellvector')
    call realloc(kim_icell,3_i4,kim_maxncell,ierror)
    if (ierror.ne.0) call outofmemory('set_kim_cellimages','kim_icell')
  endif
!
!  Set numbers of central cells
!
  kim_midcell  = (kim_ncell + 1)/2
  kim_midicell = (kim_nicell + 1)/2
!
!  Build cell vectors
!
  xcdi = - (kim_ncells(1)+1)*r1x
  ycdi = - (kim_ncells(1)+1)*r1y
  zcdi = - (kim_ncells(1)+1)*r1z
  iim = 0
!
!  Loop over unit cells
!
  do ii = - kim_ncells(1),kim_ncells(1)
    xcdi = xcdi + r1x
    ycdi = ycdi + r1y
    zcdi = zcdi + r1z
    xcdj = xcdi - (kim_ncells(2)+1)*r2x
    ycdj = ycdi - (kim_ncells(2)+1)*r2y
    zcdj = zcdi - (kim_ncells(2)+1)*r2z
    do jj = - kim_ncells(2),kim_ncells(2)
      xcdj = xcdj + r2x
      ycdj = ycdj + r2y
      zcdj = zcdj + r2z
      xcrd = xcdj - (kim_ncells(3)+1)*r3x
      ycrd = ycdj - (kim_ncells(3)+1)*r3y
      zcrd = zcdj - (kim_ncells(3)+1)*r3z
      do kk = - kim_ncells(3),kim_ncells(3)
        iim = iim + 1
        xcrd = xcrd + r3x
        ycrd = ycrd + r3y
        zcrd = zcrd + r3z
        kim_cellvector(1,iim) = xcrd
        kim_cellvector(2,iim) = ycrd
        kim_cellvector(3,iim) = zcrd
        kim_icell(1,iim) = ii
        kim_icell(2,iim) = jj
        kim_icell(3,iim) = kk
      enddo
    enddo
  enddo
!
#ifdef TRACE
  call trace_out('set_kim_cellimages')
#endif
  return
  end
!
  subroutine set_kim_neighbours
!
!  Sets up the neighbour list for a KIM model calculation
!
!  10/12 Created from bondordermd.f90
!   8/18 Modified for version 2 of OpenKIM
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
!  Copyright Curtin University 2018
!
!  Julian Gale, CIC, Curtin University, August 2018
!
  use datatypes
  use control,        only : keyword
  use current
  use iochannels
  use kim_cellimages
  use kim_models
  use neighbours
  use parallel
  use spatial,        only : lspatialok
  use spatial,        only : natomcell
  use spatial,        only : natomnodeptr
  use spatial,        only : natompernode
  use spatial,        only : ncellsearch
  use spatial,        only : nspcell
  use spatial,        only : nspcellat
  use spatial,        only : nspcellatptr
  use spatial,        only : nspcellat1ptr
  use spatial,        only : nspcellatptrcell
  use spatial,        only : xinbox
  use spatial,        only : yinbox
  use spatial,        only : zinbox
  use times
!
  implicit none
!
!  Local variables
!
  integer(i4)                                      :: i
  integer(i4)                                      :: ic
  integer(i4)                                      :: ii
  integer(i4)                                      :: imx
  integer(i4)                                      :: imy
  integer(i4)                                      :: imz
  integer(i4)                                      :: ind
  integer(i4)                                      :: ind2
  integer(i4)                                      :: indn
  integer(i4)                                      :: ix
  integer(i4)                                      :: iy
  integer(i4)                                      :: iz
  integer(i4)                                      :: j
  integer(i4)                                      :: jc
  integer(i4)                                      :: jj
  integer(i4)                                      :: k
  integer(i4)                                      :: kc
  integer(i4)                                      :: kmax
  integer(i4)                                      :: l
  integer(i4)                                      :: maxxy
  integer(i4)                                      :: maxx
  integer(i4)                                      :: n1j
  integer(i4)                                      :: nati
  integer(i4)                                      :: natj
  integer(i4)                                      :: ndone
  integer(i4), dimension(:),     allocatable, save :: ndoneptr
  integer(i4)                                      :: nj
  integer(i4)                                      :: nn
  integer(i4)                                      :: nnshell
  integer(i4)                                      :: nsplower(3)
  integer(i4)                                      :: nspupper(3)
  integer(i4)                                      :: ntypi
  integer(i4)                                      :: ntypj
  integer(i4)                                      :: status
  logical,     dimension(:),     allocatable, save :: lalreadydone
  logical,     dimension(:),     allocatable, save :: latomdone
  logical,     dimension(:),     allocatable, save :: latomdone2
  real(dp)                                         :: cut2
  real(dp)                                         :: rij
  real(dp)                                         :: r2
  real(dp)                                         :: xi
  real(dp)                                         :: yi     
  real(dp)                                         :: zi
  real(dp)                                         :: xji
  real(dp)                                         :: yji
  real(dp)                                         :: zji
  real(dp)                                         :: xji0
  real(dp)                                         :: yji0
  real(dp)                                         :: zji0
#ifdef TRACE
  call trace_in('set_kim_neighbours')
#endif
!
!  Allocate local memory that does not depend on maxneighk
!
  allocate(lalreadydone(numat),stat=status)
  if (status/=0) call outofmemory('kimmd','lalreadydone')
  allocate(latomdone(numat),stat=status)
  if (status/=0) call outofmemory('kimmd','latomdone')
  allocate(latomdone2(numat),stat=status)
  if (status/=0) call outofmemory('kimmd','latomdone2')
  allocate(ndoneptr(numat),stat=status)
  if (status/=0) call outofmemory('kimmd','ndoneptr')
!
!  Call memory allocation routine to ensure that neighbour list is initialise with the right size
!
  call changemaxneighk
!********************************
!  Set square of cutoff radius  *
!********************************
  cut2 = kim_cutoff**2
!****************************************
!  Calculate neighbour lists for atoms  *
!****************************************
!
!  Set up logical array of atoms done, so that only those needed are done in parallel
!
  latomdone(1:numat) = .false.
!
!  Setup cell images 
!
  call set_kim_cellimages
!
!  Compute neighbour list
!
  call get_kim_neighbours(cut2,latomdone)
  if (nprocs.gt.1) then
!******************************
!  Parallel additional setup  *
!******************************
!
!  Loop over atoms again to do those that are neighbours of the main atoms
!  since these will be needed in the energy evaluation. This process has to
!  be done once (in contrast to Brenner potential) since there are no
!  torsional terms. However, construct is left in here in case we want to 
!  add torsions later!
!
    if (lspatialok) then 
!******************** 
!  Spatial version  *
!********************
      maxxy = nspcell(1)*nspcell(2)
      maxx  = nspcell(1)
!                 
!  Initialise atom done once pointer
!                   
      ndone = 0
      do i = 1,numat
        lalreadydone(i) = .false.
      enddo
!
!  Loop over shells
!
      do nnshell = 1,1
        latomdone2(1:numat) = latomdone(1:numat)
        if (nnshell.eq.1) then
          kmax = natompernode
        else
          kmax = numat
        endif
        do kc = 1,kmax
          if (nnshell.eq.1) then
            k = natomnodeptr(kc)
          else
            k = kc
          endif
          if (latomdone2(k)) then
            do l = 1,nneighk(k)
              i = neighkno(l,k)
              if (.not.latomdone(i)) then
                nati = nat(i)
                ntypi = nftype(i)
                nneighk(i) = 0
!
!  Find cell containing central image of i 
!
                ind = natomcell(i)
                ind2 = ind - 1
                iz = ind2/maxxy
                ind2 = ind2 - maxxy*iz
                iy = ind2/maxx
                ix = ind2 - maxx*iy + 1
                iy = iy + 1 
                iz = iz + 1
!
                xi = xinbox(i)
                yi = yinbox(i)
                zi = zinbox(i)
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
!  Loop over neighbouring cells
!
                do imz = nsplower(3),nspupper(3)
                  do imy = nsplower(2),nspupper(2)
                    do imx = nsplower(1),nspupper(1)
                      indn = (imz-1)*maxxy + (imy-1)*maxx + imx
!
!  Loop over atoms within neighbouring cells
!
                      nj = nspcellat(indn)
                      n1j = nspcellat1ptr(indn)
                      do jj = 1,nj
                        j = nspcellatptr(n1j+jj)
                        jc = nspcellatptrcell(n1j+jj)
!
!  Exclude self term
!               
                        if (.not.lalreadydone(j)) then
                          if (i.ne.j.or.ind.ne.indn) then
                            if (latomdone(j)) then
!
!  Atom has already been done and therefore information can be copied
!
                              do ii = 1,nneighk(j)
                                if (neighkno(ii,j).eq.i) then
!
!  Check size of neighbour list arrays
!
                                  if (nneighk(i).ge.maxneighk) then
                                    maxneighk = nneighk(i) + 10
                                    call changemaxneighk
                                  endif
!
!  Add neighbour to list
!
                                  nneighk(i) = nneighk(i) + 1
                                  neighkno(nneighk(i),i) = j
                                  rneighk(nneighk(i),i) = rneighk(ii,j)
                                  neighkcell(1,nneighk(i),i) = ivec2cell(1,jc)
                                  neighkcell(2,nneighk(i),i) = ivec2cell(2,jc)
                                  neighkcell(3,nneighk(i),i) = ivec2cell(3,jc)
                                  xyzneighk(1,nneighk(i),i) = - xyzneighk(1,ii,j)
                                  xyzneighk(2,nneighk(i),i) = - xyzneighk(2,ii,j)
                                  xyzneighk(3,nneighk(i),i) = - xyzneighk(3,ii,j)
                                endif
                              enddo
!
!  Set pointer to avoid repetition of this atom in the copy phase
!
                              ndone = ndone + 1
                              ndoneptr(ndone) = j
                              lalreadydone(j) = .true.
                            else
                              natj = nat(j)
                              ntypj = nftype(j)
!
!  Set centre cell coordinate differences
!
                              xji = xvec2cell(jc) + xinbox(j) - xi
                              yji = yvec2cell(jc) + yinbox(j) - yi
                              zji = zvec2cell(jc) + zinbox(j) - zi
                              r2 = xji*xji + yji*yji + zji*zji
!
!  Check distance squared against cutoff
!
                              if (r2.lt.cut2) then
!
!  Check size of neighbour list arrays
!
                                if (nneighk(i).ge.maxneighk) then
                                  maxneighk = nneighk(i) + 10
                                  call changemaxneighk
                                endif
!
!  Add neighbour to list
!
                                rij = sqrt(r2)
                                nneighk(i) = nneighk(i) + 1
                                neighkno(nneighk(i),i) = j
                                rneighk(nneighk(i),i) = rij
                                neighkcell(1,nneighk(i),i) = ivec2cell(1,jc)
                                neighkcell(2,nneighk(i),i) = ivec2cell(2,jc)
                                neighkcell(3,nneighk(i),i) = ivec2cell(3,jc)
                                xyzneighk(1,nneighk(i),i) = xji
                                xyzneighk(2,nneighk(i),i) = yji
                                xyzneighk(3,nneighk(i),i) = zji
!
!  Set pointer to avoid repetition of this atom in the search phase
!
                                ndone = ndone + 1
                                ndoneptr(ndone) = j
                                lalreadydone(j) = .true.
                              endif
                            endif
                          endif
                        endif
                      enddo
                    enddo
                  enddo
                enddo
!
!  Set flag for this atom to indicate that it has been done
!
                latomdone(i) = .true.
!
!  Clear already done pointers for this atom
!
                do j = 1,ndone
                  lalreadydone(ndoneptr(j)) = .false.
                enddo
                ndone = 0
              endif
            enddo
          endif
        enddo
      enddo
    else
!************************
!  Non-spatial version  *
!************************
      do nnshell = 1,1
        latomdone2(1:numat) = latomdone(1:numat)
        if (nnshell.eq.1) then
          kmax = natompernode   
        else
          kmax = numat
        endif
        do kc = 1,kmax
          if (nnshell.eq.1) then
            k = natomnodeptr(kc)
          else
            k = kc   
          endif
          if (latomdone2(k)) then
            do l = 1,nneighk(k)
              i = neighkno(l,k)
!
              if (.not.latomdone(i)) then
                nneighk(i) = 0
                nati = nat(i)
                ntypi = nftype(i)
!
!  Loop over atoms
!  
                do j = 1,numat
                  if (latomdone(j)) then
!  
!  Atom has already been done and therefore information can be copied      
!  
                    do ii = 1,nneighk(j)
                      if (neighkno(ii,j).eq.i) then
!
!  Check size of neighbour list arrays
!
                        if (nneighk(i).ge.maxneighk) then
                          maxneighk = nneighk(i) + 10
                          call changemaxneighk
                        endif
!
!  Add neighbour to list
!
                        nneighk(i) = nneighk(i) + 1
                        neighkno(nneighk(i),i) = j
                        neighkcell(1,nneighk(i),i) = - neighkcell(1,ii,j)
                        neighkcell(2,nneighk(i),i) = - neighkcell(2,ii,j)
                        neighkcell(3,nneighk(i),i) = - neighkcell(3,ii,j)
                        rneighk(nneighk(i),i) = rneighk(ii,j)
                        xyzneighk(1,nneighk(i),i) = - xyzneighk(1,ii,j)
                        xyzneighk(2,nneighk(i),i) = - xyzneighk(2,ii,j)
                        xyzneighk(3,nneighk(i),i) = - xyzneighk(3,ii,j)
                      endif     
                    enddo       
                  else
                    natj = nat(j)  
                    ntypj = nftype(j)
!
!  Set centre cell coordinate differences
!
                    xji0 = xclat(j) - xclat(i)
                    yji0 = yclat(j) - yclat(i)
                    zji0 = zclat(j) - zclat(i)
!
!  Loop over unit cells
!
                    do ii = 1,kim_ncell
!
!  Exclude self term
!
                      if (i.ne.j.or.ii.ne.kim_midcell) then
                        xji = xji0 + kim_cellvector(1,ii)
                        yji = yji0 + kim_cellvector(2,ii)
                        zji = zji0 + kim_cellvector(3,ii)
                        r2 = xji*xji + yji*yji + zji*zji
!
!  Check distance squared against cutoff
!
                        if (r2.lt.cut2) then
!
!  Check size of neighbour list arrays
!
                          if (nneighk(i).ge.maxneighk) then
                            maxneighk = nneighk(i) + 10
                            call changemaxneighk
                          endif
!
!  Add neighbour to list
!
                          rij = sqrt(r2)
                          nneighk(i) = nneighk(i) + 1
                          neighkno(nneighk(i),i) = j
                          neighkcell(1,nneighk(i),i) = kim_icell(1,ii)
                          neighkcell(2,nneighk(i),i) = kim_icell(2,ii)
                          neighkcell(3,nneighk(i),i) = kim_icell(3,ii)
                          rneighk(nneighk(i),i) = rij
                          xyzneighk(1,nneighk(i),i) = xji
                          xyzneighk(2,nneighk(i),i) = yji
                          xyzneighk(3,nneighk(i),i) = zji
                        endif
                      endif
                    enddo
                  endif
                enddo
!
!  Set flag for this atom to indicate that it has been done
!
                latomdone(i) = .true.
              endif
            enddo
          endif
        enddo
      enddo
    endif
  endif
!
!  Set maxneighk in KIM location
!
  kim_maxneighbour = maxneighk + 1
!*********************************************
!  Print debugging information if requested  *
!*********************************************
  if (ioproc.and.index(keyword,'debu').ne.0) then
    write(ioout,'(/,''  Number of neighbours for atoms :'',/)')
    write(ioout,'(''  i'',2x,''Nat'',3x,''N'')')
    do ic = 1,natompernode
      i = natomnodeptr(ic)
      write(ioout,'(i4,1x,i3,1x,i6)') i,nat(i),nneighk(i)
    enddo
    write(ioout,'(/,''  Neighbours of atoms :'',/)')
    do ic = 1,natompernode
      i = natomnodeptr(ic)
      write(ioout,'(i4,8(1x,i4))') i,(neighkno(nn,i),nn=1,nneighk(i))
    enddo
  endif
!
!  Free local memory
!
  deallocate(ndoneptr,stat=status)
  if (status/=0) call deallocate_error('kimmd','ndoneptr')
  deallocate(latomdone2,stat=status)
  if (status/=0) call deallocate_error('kimmd','latomdone2')
  deallocate(latomdone,stat=status)
  if (status/=0) call deallocate_error('kimmd','latomdone')
  deallocate(lalreadydone,stat=status)
  if (status/=0) call deallocate_error('kimmd','lalreadydone')
#ifdef TRACE
  call trace_in('set_kim_neighbours')
#endif
!
  end subroutine set_kim_neighbours
!***********************
!  get_kim_neighbours  *
!***********************
  subroutine get_kim_neighbours(cut2,latomdone)
!
!  Finds neighbour list for a KIM potential
!
!  On entry : 
!
!  cut2          = maximum potential cutoff from KIM squared
!
!  On exit :
!
!  nneighk       = number of neighbours for each atom
!  neighkno      = pointer to atom numbers of neighbours for each atom
!  neighkcell    = pointer to cell index for neighbours of each atom
!  rneighk       = distances of neighbours for each atom
!  xyzneighk     = x/y/z components of distances of neighbours for each atom
!  latomdone     = logical indicating whether atom was done locally
!
!  10/12 Created from getBOneighbour
!   8/18 Cell search changed to use range set based on kim_influence
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
!  Copyright Curtin University 2018
!
!  Julian Gale, CIC, Curtin University, August 2018
!
  use datatypes
  use current,        only : numat
  use current,        only : ivec2cell
  use current,        only : xclat, yclat, zclat
  use current,        only : xvec2cell, yvec2cell, zvec2cell
  use kim_cellimages
  use parallel
  use spatial,        only : lbuffercell
  use spatial,        only : lspatialok
  use spatial,        only : natompernode
  use spatial,        only : natomnodeptr
  use spatial,        only : ncellsearch
  use spatial,        only : ncellpernode
  use spatial,        only : ncellnodeptr
  use spatial,        only : nspcell
  use spatial,        only : nspcellat
  use spatial,        only : nspcellatptr
  use spatial,        only : nspcellat1ptr
  use spatial,        only : nspcellatptrcell
  use spatial,        only : xinbox
  use spatial,        only : yinbox
  use spatial,        only : zinbox
  implicit none
!
!  Passed variables
!
  real(dp),    intent(in)                        :: cut2
  logical,     intent(out)                       :: latomdone(numat)
!
!  Local variables
!
  integer(i4)                                    :: i
  integer(i4)                                    :: ia
  integer(i4)                                    :: ic
  integer(i4)                                    :: ii
  integer(i4)                                    :: imx
  integer(i4)                                    :: imy
  integer(i4)                                    :: imz
  integer(i4)                                    :: ind1
  integer(i4)                                    :: ind2
  integer(i4)                                    :: ix
  integer(i4)                                    :: ixyz
  integer(i4)                                    :: iy
  integer(i4)                                    :: iz
  integer(i4)                                    :: j
  integer(i4)                                    :: jc
  integer(i4)                                    :: jj
  integer(i4)                                    :: maxxy
  integer(i4)                                    :: maxx
  integer(i4)                                    :: n1
  integer(i4)                                    :: n1j
  integer(i4)                                    :: ni
  integer(i4)                                    :: nj
  integer(i4)                                    :: nsplower(3)
  integer(i4)                                    :: nspupper(3)
  real(dp)                                       :: rij
  real(dp)                                       :: r2
  real(dp)                                       :: xi
  real(dp)                                       :: yi     
  real(dp)                                       :: zi
  real(dp)                                       :: xji
  real(dp)                                       :: yji
  real(dp)                                       :: zji
  real(dp)                                       :: xji0
  real(dp)                                       :: yji0
  real(dp)                                       :: zji0
#ifdef TRACE
  call trace_in('get_kim_neighbours')
#endif
!
!  Initialise nneighk
!
  nneighk(1:numat) = 0
!****************************************
!  Calculate neighbour lists for atoms  *
!****************************************
  if (lspatialok) then
    maxxy = nspcell(1)*nspcell(2)
    maxx  = nspcell(1)
!  
!  Loop over all spatial cells looking for non-buffer cells
!     
    do ixyz = 1,ncellpernode
      if (.not.lbuffercell(ixyz)) then
        ind1 = ncellnodeptr(ixyz)
        ind2 = ind1 - 1
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
        ni = nspcellat(ind1)
        n1 = nspcellat1ptr(ind1)
!  
!  Loop over atoms in the cell finding neighbours
!     
        do ii = 1,ni
          i = nspcellatptr(n1+ii)
          ic = nspcellatptrcell(n1+ii)
          latomdone(i) = .true.
!       
!  Set coordinates
!
          xi = xinbox(i) + xvec2cell(ic)
          yi = yinbox(i) + yvec2cell(ic)
          zi = zinbox(i) + zvec2cell(ic)
!
!  Loop over neighbouring cells
!
          do imz = nsplower(3),nspupper(3)
            do imy = nsplower(2),nspupper(2)
              do imx = nsplower(1),nspupper(1)
                ind2 = (imz-1)*maxxy + (imy-1)*maxx + imx
!                       
!  Loop over atoms within neighbouring cells  
!                         
                nj = nspcellat(ind2)
                n1j = nspcellat1ptr(ind2)
                do jj = 1,nj
                  j = nspcellatptr(n1j+jj)
!                     
!  Exclude self term    
!                         
                  if (i.ne.j.or.ind1.ne.ind2) then
                    jc = nspcellatptrcell(n1j+jj)
!                             
!  Set centre cell coordinate differences
!  
                    xji = xvec2cell(jc) + xinbox(j) - xi
                    yji = yvec2cell(jc) + yinbox(j) - yi
                    zji = zvec2cell(jc) + zinbox(j) - zi
!  
                    r2 = xji*xji + yji*yji + zji*zji
!
!  Check whether j is within two body cut-off
!
                    if (r2.lt.cut2) then
!
!  Check whether neighbour list arrays are large enough
!
                      if (nneighk(i).ge.maxneighk) then
                        maxneighk = nneighk(i) + 10
                        call changemaxneighk
                      endif
!
!  Save information for new neighbour
!
                      rij = sqrt(r2)
                      nneighk(i) = nneighk(i) + 1
                      neighkno(nneighk(i),i) = j
                      rneighk(nneighk(i),i) = rij
                      neighkcell(1,nneighk(i),i) = ivec2cell(1,jc)
                      neighkcell(2,nneighk(i),i) = ivec2cell(2,jc)
                      neighkcell(3,nneighk(i),i) = ivec2cell(3,jc)
                      xyzneighk(1,nneighk(i),i) = xji
                      xyzneighk(2,nneighk(i),i) = yji
                      xyzneighk(3,nneighk(i),i) = zji
                    endif
                  endif
!
                enddo
              enddo
            enddo
          enddo             
!
        enddo                 
!
      endif
    enddo                       
  else
!
    do ia = 1,natompernode
      i = natomnodeptr(ia)
      latomdone(i) = .true.
!
!  Loop over atoms
!
      do j = 1,numat
!
!  Set centre cell coordinate differences
!
        xji0 = xclat(j) - xclat(i)
        yji0 = yclat(j) - yclat(i)
        zji0 = zclat(j) - zclat(i)
!
!  Loop over unit cells
!
        do ii = 1,kim_ncell
!
!  Exclude self term
!
          if (i.ne.j.or.ii.ne.kim_midcell) then
            xji = xji0 + kim_cellvector(1,ii)
            yji = yji0 + kim_cellvector(2,ii)
            zji = zji0 + kim_cellvector(3,ii)
            r2 = xji*xji + yji*yji + zji*zji
!
!  Check distance squared against cutoff
!
            if (r2.lt.cut2) then
!
!  Check whether neighbour list arrays are large enough
!
              if (nneighk(i).ge.maxneighk) then
                maxneighk = nneighk(i) + 10
                call changemaxneighk
              endif
!
!  Add neighbour to lists
!
              rij = sqrt(r2)
              nneighk(i) = nneighk(i) + 1
              neighkno(nneighk(i),i) = j
              neighkcell(1,nneighk(i),i) = kim_icell(1,ii)
              neighkcell(2,nneighk(i),i) = kim_icell(2,ii)
              neighkcell(3,nneighk(i),i) = kim_icell(3,ii)
              rneighk(nneighk(i),i) = rij
              xyzneighk(1,nneighk(i),i) = xji
              xyzneighk(2,nneighk(i),i) = yji
              xyzneighk(3,nneighk(i),i) = zji
            endif
          endif
        enddo
      enddo
    enddo
  endif
#ifdef TRACE
  call trace_out('get_kim_neighbours')
#endif
!
  end subroutine get_kim_neighbours
!
!  changemaxneighk - reallocates memory for neighbour lists
!
  subroutine changemaxneighk
!
!  Alters the size of the arrays associated with maxneighk and maxat relevant to neighbour list
!
!  10/12 Created
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
!  Copyright Curtin University 2018
!
!  Julian Gale, CIC, Curtin University, August 2018
!
  use datatypes
  use current,        only : maxat
  use reallocate
  implicit none
!
!  Local variables
!
  integer(i4)       :: ierror, i
  integer(i4), save :: oldmaxat = 0
#ifdef TRACE
  call trace_in('changemaxneighk')
#endif
!
  call realloc(nneighk,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxneighk','nneighk')
  call realloc(neighkno,maxneighk,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxneighk','neighkno')
  call realloc(neighkcell,3_i4,maxneighk,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxneighk','neighkcell')
  call realloc(rneighk,maxneighk,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxneighk','rneighk')
  call realloc(xyzneighk,3_i4,maxneighk,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxneighk','xyzneighk')
!
!  Initialise new parts of data arrays
!
  if (maxat.gt.oldmaxat) then
    do i = oldmaxat+1,maxat
      nneighk(i) = 0
    enddo
  endif
!
!  Save current value of maxat for next call
!
  oldmaxat = maxat
#ifdef TRACE
  call trace_out('changemaxneighk')
#endif
!
  return
  end subroutine changemaxneighk

  subroutine set_kim_neighbour_list(nm)
!
!  Transfer data from GULP neighbour list into arrays for KIM
!
!  For now the neighbour list for each atom is just added
!  to the list of atoms in the unit cell.
!
!   8/18 Created
!   9/18 Modified for multiple models
!
!  On entry:
!
!  nm = model number for OpenKIM
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
!  Copyright Curtin University 2018
!
!  Julian Gale, CIC, Curtin University, September 2018
!
  use, intrinsic :: iso_c_binding
  use datatypes
  use current
  use kim_cellimages
  use kim_models
  use neighbours
  use parallel
  use spatial,        only : natomnodeptr
  use spatial,        only : natompernode
!
  implicit none
!
!  Passed variables
!
  integer(i4),                        intent(in)   :: nm
!
!  Local variables
!
  integer(i4)                                      :: i
  integer(i4)                                      :: ia
  integer(i4)                                      :: ic1
  integer(i4)                                      :: ic2
  integer(i4)                                      :: ic3
  integer(i4)                                      :: ii
  integer(i4)                                      :: iim
  integer(i4)                                      :: imid
  integer(i4)                                      :: ind
  integer(i4)                                      :: j
  integer(i4)                                      :: jc1
  integer(i4)                                      :: jc2
  integer(i4)                                      :: jc3
  integer(i4)                                      :: ji
  integer(i4)                                      :: jj
  integer(i4)                                      :: k
  integer(i4)                                      :: kk
  integer(i4)                                      :: kj
  integer(i4)                                      :: nc
  real(dp)                                         :: conversion
  real(dp)                                         :: cut
  real(dp)                                         :: xcdi
  real(dp)                                         :: ycdi
  real(dp)                                         :: zcdi
  real(dp)                                         :: xcdj
  real(dp)                                         :: ycdj
  real(dp)                                         :: zcdj
  real(dp)                                         :: xcrd
  real(dp)                                         :: ycrd
  real(dp)                                         :: zcrd
!
!  Loop over atoms in cell
!
  conversion = 1.0_dp/kim_length_conversion(nm)
!
!  Build coordinates and species based on cell images
!
  iim = 0
  xcdi = - (kim_nicells(1)+1)*r1x
  ycdi = - (kim_nicells(1)+1)*r1y
  zcdi = - (kim_nicells(1)+1)*r1z
  do ii = - kim_nicells(1),kim_nicells(1)
    xcdi = xcdi + r1x
    ycdi = ycdi + r1y
    zcdi = zcdi + r1z
    xcdj = xcdi - (kim_nicells(2)+1)*r2x
    ycdj = ycdi - (kim_nicells(2)+1)*r2y
    zcdj = zcdi - (kim_nicells(2)+1)*r2z
    do jj = - kim_nicells(2),kim_nicells(2)
      xcdj = xcdj + r2x
      ycdj = ycdj + r2y
      zcdj = zcdj + r2z
      xcrd = xcdj - (kim_nicells(3)+1)*r3x
      ycrd = ycdj - (kim_nicells(3)+1)*r3y
      zcrd = zcdj - (kim_nicells(3)+1)*r3z
      do kk = - kim_nicells(3),kim_nicells(3)
        xcrd = xcrd + r3x
        ycrd = ycrd + r3y
        zcrd = zcrd + r3z
        do i = 1,numat
          iim = iim + 1
          kim_nspecies(iim) = kim_nspec(nspecptr(i),nm)
          kim_coord(1,iim) = (xclat(i) + xcrd)*conversion
          kim_coord(2,iim) = (yclat(i) + ycrd)*conversion
          kim_coord(3,iim) = (zclat(i) + zcrd)*conversion
        enddo
      enddo
    enddo
  enddo
!
!  Initialise contributing flag and number of neighbours to zero as not all atoms are needed
!
  kim_ncontributing(1:kim_numat,1:kim_num_nlist(nm)) = 0
  kim_neighbourlist(1,1:kim_numat,1:kim_num_nlist(nm)) = 0
!
!  Loop over central atoms to set up neighbour lists
!
  imid = (kim_midicell - 1)*numat
!
!  If there is a single neighbour list then we can just copy information
!
  if (kim_num_nlist(nm).eq.1) then
    do ia = 1,natompernode
      i = natomnodeptr(ia)
      iim = imid + i    ! This is the central cell atom in the influence list
!
!  Set atom to be contributing
!
      kim_ncontributing(iim,1) = 1
!
!  Set number of neighbours
!
      kim_neighbourlist(1,iim,1) = nneighk(i)
!
!  Loop over neighbours
!
      do jj = 1,nneighk(i)
        j = neighkno(jj,i)
        ic1 = neighkcell(1,jj,i)
        ic2 = neighkcell(2,jj,i)
        ic3 = neighkcell(3,jj,i)
!
!  Find index for cell of j and generate atom number
!
        ind = (ic1 + kim_nicells(1))*(2*kim_nicells(2) + 1)*(2*kim_nicells(3) + 1) + &
              (ic2 + kim_nicells(2))*(2*kim_nicells(3) + 1) + &
              (ic3 + kim_nicells(3))
        ji = ind*numat + j
        if (ji.lt.1.or.ji.gt.kim_numat) then
          call outerror('atom in neighbour list is out of range of influence cells',0_i4)
          call stopnow('set_kim_neighbour_list')
        endif
        kim_neighbourlist(1+jj,iim,1) = ji
!
        if (kim_ncontributing_only(1,nm).eq.0) then
          if (ji.le.imid.or.ji.gt.(imid+numat)) then
            if (kim_neighbourlist(1,ji,1).eq.0) then
!
!  If neighbour has been handled already then add its neighbour list if required
!
              kim_ncontributing(ji,1) = 0
              kim_neighbourlist(1,ji,1) = nneighk(j)
              do kk = 1,nneighk(j)
                k = neighkno(kk,j)
                jc1 = neighkcell(1,kk,j)
                jc2 = neighkcell(2,kk,j)
                jc3 = neighkcell(3,kk,j)
                ind = (ic1 + jc1 + kim_nicells(1))*(2*kim_nicells(2) + 1)*(2*kim_nicells(3) + 1) + &
                      (ic2 + jc2 + kim_nicells(2))*(2*kim_nicells(3) + 1) + &
                      (ic3 + jc3 + kim_nicells(3))
                kj = ind*numat + k
                if (kj.lt.1.or.kj.gt.kim_numat) then
                  call outerror('atom in neighbour list is out of range of influence cells',0_i4)
                  call stopnow('set_kim_neighbour_list')
                endif
                kim_neighbourlist(1+kk,ji,1) = kj
              enddo
            endif
          endif
        else
!
!  Add connection back to central cell if atom is not in the cell already
!
          if (ji.le.imid.or.ji.gt.(imid+numat)) then
            kim_ncontributing(ji,1) = 0
            kim_neighbourlist(1,ji,1) = kim_neighbourlist(1,ji,1) + 1
            kim_neighbourlist(1+kim_neighbourlist(1,ji,1),ji,1) = iim
          endif
        endif
      enddo
    enddo
  else
!
!  In order to handle parallel case set flags from contributing separately to neighbour lists
!
    do ia = 1,natompernode
      i = natomnodeptr(ia)
      iim = imid + i    ! This is the central cell atom in the influence list
      do nc = 1,kim_num_nlist(nm)
        kim_ncontributing(iim,nc) = 1
      enddo
    enddo
!
!  Loop over number of cutoffs
!
    do nc = 1,kim_num_nlist(nm)
      cut = kim_cutoffs(nc,nm)
      do i = 1,numat
        iim = imid + i    ! This is the central cell atom in the influence list
!
!  Loop over neighbours
!
        do jj = 1,nneighk(i)
!
!  Check cutoff
!
          if (rneighk(jj,i).le.cut.or.kim_any_ncontributing_only(nm)) then
            j = neighkno(jj,i)
            ic1 = neighkcell(1,jj,i)
            ic2 = neighkcell(2,jj,i)
            ic3 = neighkcell(3,jj,i)
!
!  Find index for cell of j and generate atom number
!
            ind = (ic1 + kim_nicells(1))*(2*kim_nicells(2) + 1)*(2*kim_nicells(3) + 1) + &
                  (ic2 + kim_nicells(2))*(2*kim_nicells(3) + 1) + &
                  (ic3 + kim_nicells(3))
            ji = ind*numat + j
            if (ji.lt.1.or.ji.gt.kim_numat) then
              call outerror('atom in neighbour list is out of range of influence cells',0_i4)
              call stopnow('set_kim_neighbour_list')
            endif
            if (rneighk(jj,i).le.cut) then
!
!  If below cutoff then include in the current neighbour list
!
              kim_neighbourlist(1,iim,nc) = kim_neighbourlist(1,iim,nc) + 1
              kim_neighbourlist(1+kim_neighbourlist(1,iim,nc),iim,nc) = ji
            endif
!
            if (kim_any_ncontributing_only(nm)) then
              if (ji.le.imid.or.ji.gt.(imid+numat)) then
                if (kim_neighbourlist(1,ji,nc).eq.0) then
!
!  If neighbour has been handled already then add its neighbour list if required
!
                  kim_ncontributing(ji,nc) = 0
                  do kk = 1,nneighk(j)
!
!  Is atom within cutoff for this neighbour list?
!
                    if (rneighk(kk,j).le.cut) then
                      kim_neighbourlist(1,ji,nc) = kim_neighbourlist(1,ji,nc) + 1
!
                      k = neighkno(kk,j)
                      jc1 = neighkcell(1,kk,j)
                      jc2 = neighkcell(2,kk,j)
                      jc3 = neighkcell(3,kk,j)
!
                      ind = (ic1 + jc1 + kim_nicells(1))*(2*kim_nicells(2) + 1)*(2*kim_nicells(3) + 1) + &
                            (ic2 + jc2 + kim_nicells(2))*(2*kim_nicells(3) + 1) + &
                            (ic3 + jc3 + kim_nicells(3))
                      kj = ind*numat + k
                      if (kj.lt.1.or.kj.gt.kim_numat) then
                        call outerror('atom in neighbour list is out of range of influence cells',0_i4)
                        call stopnow('set_kim_neighbour_list')
                      endif
                      kim_neighbourlist(1+kim_neighbourlist(1,ji,nc),ji,nc) = kj
                    endif
                  enddo
                endif
              endif
            else
!
!  Add connection back to central cell if atom is not in the cell already
!
              if (ji.le.imid.or.ji.gt.(imid+numat)) then
                kim_ncontributing(ji,nc) = 0
                kim_neighbourlist(1,ji,nc) = kim_neighbourlist(1,ji,nc) + 1
                kim_neighbourlist(1+kim_neighbourlist(1,ji,nc),ji,nc) = iim
              endif
            endif
          endif
        enddo
      enddo
    enddo
  endif

  end subroutine set_kim_neighbour_list

  subroutine get_neigh(data_object, number_of_neighbor_lists, cutoffs, &
  neighbor_list_index, request, numnei, pnei1part, ierr) bind(c)

  use kim_models
  implicit none

  !-- Transferred variables
  type(c_ptr), value, intent(in)  :: data_object
  integer(i4), value, intent(in)  :: number_of_neighbor_lists
  real(dp),           intent(in)  :: cutoffs(number_of_neighbor_lists)
  integer(i4), value, intent(in)  :: neighbor_list_index
  integer(i4), value, intent(in)  :: request
  integer(i4),        intent(out) :: numnei
  type(c_ptr),        intent(out) :: pnei1part
  integer(i4),        intent(out) :: ierr

  if (cutoffs(neighbor_list_index).gt.kim_cutoff) then
    call outerror('cutoff exceeds that of neighbour list',0_i4)
    call stopnow('get_neigh')
  endif

  if ((request.gt.kim_numat) .or. (request.lt.1)) then
    call outerror('KIM atom is out of range in get_neigh',0_i4)
    call stopnow('get_neigh')
  endif

  ! set the returned number of neighbors for the returned part
  numnei = kim_neighbourlist(1,request,neighbor_list_index)

  ! set the location for the returned neighbor list
  pnei1part = c_loc(kim_neighbourlist(2,request,neighbor_list_index))

  ierr = 0
  return

  end subroutine get_neigh

  end module kim_functions
