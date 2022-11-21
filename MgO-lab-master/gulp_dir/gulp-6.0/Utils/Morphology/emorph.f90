program emorph
!
!  Program to compute the energy of a finite crystal based on surface area and volume
!
!  Conditions of use:
!
!  Nanomorph is available free of charge to academic institutions
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
!  Julian Gale, Curtin University, October 2018
!
  implicit none
!
!  The following variables define the limits of the code in terms of problem size and
!  the amount of memory that will be used.
!
  integer*4                          :: maxpav = 6       ! Maximum number of Planes At a Vertex
  integer*4                          :: maxhkl = 40      ! Maximum number of surfaces
  integer*4                          :: maxvertex = 100  ! Maximum number of vertices
!
  character(len=132)                 :: line             ! Character string used for reading lines in
  logical                            :: leof             ! Logical used to test whether end of input has been reached
  logical                            :: linputrv         ! Flag as to whether lattice vectors have been input
  logical                            :: lvertexok        ! Flag as to whether vertex is part of morphology
  integer*4                          :: i                ! Loop counter
  integer*4                          :: il               ! Local scalar
  integer*4                          :: ind              ! Index for convoluting 2 integers into a lower-half triangular form
  integer*4                          :: iv               ! Local scalar
  integer*4                          :: iv1              ! Local scalar
  integer*4                          :: iv2              ! Local scalar
  integer*4                          :: ivv1             ! Local scalar
  integer*4                          :: ivv2             ! Local scalar
  integer*4                          :: i3               ! Integer value of 3 in integer*4 for passing to matinv
  integer*4                          :: info             ! Error flag for Lapack
  integer*4                          :: itmp(3,2)        ! Temporary array for holding edge Miller indices
  integer*4                          :: j                ! Loop counter
  integer*4                          :: k                ! Loop counter
  integer*4                          :: l                ! Loop counter
  logical                            :: lfound1          ! Logical test flag for searching in a loop
  logical                            :: lfound2          ! Logical test flag for searching in a loop
  integer*4                          :: nfaces           ! Number of visible surface planes
  integer*4                          :: nhkl             ! Number of surfaces
  integer*4                          :: nvertex          ! Number of vertices
  integer*4                          :: nzero            ! Number of zero Miller indices
  integer*4,         allocatable     :: ihkl(:,:)        ! Miller indices of surfaces
  integer*4,         allocatable     :: ifacevertex(:,:) ! Point to vertices for a visible surface
  integer*4,         allocatable     :: nfaceptr(:)      ! Pointer from nfaces to nhkl index
  integer*4,         allocatable     :: nadjfaceptr(:,:) ! Pointer from edge of nfaces to nhkl index of adjacent face
  integer*4,         allocatable     :: nfacevertex(:)   ! Number of vertices for a visible surface
  integer*4,         allocatable     :: nhklvertex(:)    ! Number of vertices for a plane
  integer*4,         allocatable     :: npavertex(:)     ! Number of planes at a vertex
  integer*4,         allocatable     :: ivertex(:,:)     ! Pointer to surfaces that meet at the vertex
  integer*4,         allocatable     :: ihklvertex(:,:)  ! Pointer to vertices for a plane
  integer*4,         allocatable     :: ifvordered(:,:)  ! Ordered list of vertices for a visible surface
  logical,           allocatable     :: lkeepvertex(:)   ! Is a vertex redundant due to being equal to another one?
  logical,           allocatable     :: lfvdone(:)       ! Flag as to whether a face vertex has already been found
  real*8,            allocatable     :: aface(:)         ! Area of visible surface
  real*8,            allocatable     :: lface(:,:)       ! Length of surface edges
  real*8,            allocatable     :: dotf(:)          ! Dot products of vectors from mid point to vertices
  real*8,            allocatable     :: ehkl(:)          ! Energy of surface
  real*8,            allocatable     :: lhkl(:)          ! Energy for edge between two surfaces
  real*8,            allocatable     :: shkl(:,:)        ! Surface normal vector from origin
  real*8,            allocatable     :: snhkl(:,:)       ! Surface normal vector from origin (normalised)
  real*8,            allocatable     :: s2hkl(:)         ! Magnitude of the surface normal
  real*8,            allocatable     :: vhkl(:,:,:)      ! Vectors that define the surface plane
  real*8,            allocatable     :: vertex(:,:)      ! Coordinates of vertex
  real*8,            allocatable     :: vmid2v(:,:)      ! Vectors of vertices from mid point of face
  real*8                             :: cp(3)            ! Cross product vector of 2 surface normals
  real*8                             :: crossp           ! Cross product norm
  real*8                             :: dotp             ! Dot product
  real*8                             :: bulkenergy       ! Bulk energy per unit cell
  real*8                             :: energy           ! Energy of system
  real*8                             :: etmp             ! Temporary scalar for edge energy
  real*8                             :: lenergy          ! Energy due to edges
  real*8                             :: senergy          ! Energy due to surfaces
  real*8                             :: dvert            ! Distance between vertices
  real*8                             :: pintersect(3,3)  ! Workspace for points along axes
  real*8                             :: rv(3,3)          ! Lattice vectors
  real*8                             :: pmat(3,3)        ! Matrix used for vertex solve
  real*8                             :: rhs(3)           ! Right-hand side vector used for vertex solve
  real*8                             :: dscale           ! Size scaling factor
  real*8                             :: snorm            ! Norm of surface normal pre-scaling
  real*8                             :: tolerance        ! Criterion for planes not to meet
  real*8                             :: volume           ! Volume of polyhedron
  real*8                             :: unitarea         ! Unit conversion factor for surface area
  real*8                             :: wrk(6)           ! Workspace for matrix inversion
  real*8                             :: vnorm1           ! Norm of a vector
  real*8                             :: vnormk           ! Norm of a vector
!###################################################################################################
!  Initialise quantities
!###################################################################################################
  i3 = 3
  nfaces = 0
  nhkl = 0
  nvertex = 0
  linputrv = .false.
  dscale = 1.0d0
  tolerance = 1.0d-6
  unitarea = 1.0d0
  bulkenergy = 0.0d0
!
  maxvertex = maxhkl*(maxhkl-1)*(maxhkl-2)/2
!
  allocate(aface(maxhkl))
  allocate(ehkl(maxhkl))
  allocate(ihkl(3,maxhkl))
  allocate(lhkl(maxhkl*(maxhkl+1)/2))
  allocate(shkl(3,maxhkl))
  allocate(snhkl(3,maxhkl))
  allocate(s2hkl(maxhkl))
  allocate(vhkl(3,2,maxhkl))
  allocate(npavertex(maxvertex))
  allocate(nfaceptr(maxhkl))
  allocate(nfacevertex(maxhkl))
  allocate(ifacevertex(maxvertex,maxhkl))
  allocate(nhklvertex(maxhkl))
  allocate(ihklvertex(maxvertex,maxhkl))
  allocate(ivertex(maxpav,maxvertex))
  allocate(vertex(maxpav,maxvertex))
  allocate(lkeepvertex(maxvertex))
!
  lhkl(1:maxhkl*(maxhkl+1)/2) = 0.0d0
!###################################################################################################
!  Read input file
!###################################################################################################
!
!  Loop over option words
!
  leof = .false.
  do while (.not.leof)
!
!  Read line
!
    read(5,'(a)',end=10,err=10) line
    if (index(line,'vect').eq.1) then
!
!  Read lattice vectors for system
!
      read(5,*) (rv(j,1),j=1,3)
      read(5,*) (rv(j,2),j=1,3)
      read(5,*) (rv(j,3),j=1,3)
      linputrv = .true.
    elseif (index(line,'size').eq.1) then
!
!  Read scaling factor for distances
!
      read(5,*) dscale
      dscale = abs(dscale)
    elseif (index(line,'bulk').eq.1) then
!
!  Read bulk energy per unit cell
!
      read(5,*) bulkenergy
    elseif (index(line,'usurf').eq.1) then
      read(5,'(a)',end=10,err=10) line
      if (index(line,'jm2').ne.0) unitarea = 1.0d0/16.021917
    elseif (index(line,'surf').eq.1) then
!
!  Read surface info
!
      nhkl = nhkl + 1
      if (nhkl.gt.maxhkl) then
        write(6,'(/,'' Error : Number of planes exceeds maxhkl! '',/)')
        stop
      endif
      read(5,*) (ihkl(j,nhkl),j=1,3),ehkl(nhkl)
      if (ehkl(nhkl).lt.0.0d0) then
        write(6,'(/,'' Error : Surface energy must be positive! '',/)')
        stop
      elseif (ehkl(nhkl).eq.0.0d0) then
        write(6,'(/,'' Error : Surface energy cannot be zero! '',/)')
        stop
      endif
    elseif (index(line,'edge').eq.1) then
!
!  Read edge info
!
      read(5,*) (itmp(j,1),j=1,3),(itmp(k,2),k=1,3),etmp
!
!  Find corresponding edges
!
      iv1 = 0
      iv2 = 0
      lfound1 = .false.
      lfound2 = .false.
      i = 0
      do while (i.lt.nhkl.and.(.not.lfound1.or..not.lfound2))
        i = i + 1
        if (.not.lfound1) then
          lfound1 = (itmp(1,1).eq.ihkl(1,i).and. &
                     itmp(2,1).eq.ihkl(2,i).and. &
                     itmp(3,1).eq.ihkl(3,i))
          if (lfound1) iv1 = i
        endif
        if (.not.lfound2) then
          lfound2 = (itmp(1,2).eq.ihkl(1,i).and. &
                     itmp(2,2).eq.ihkl(2,i).and. &
                     itmp(3,2).eq.ihkl(3,i))
          if (lfound2) iv2 = i
        endif
      enddo
      if (.not.lfound1.or..not.lfound2) then
        write(6,'(/,'' Error : Surface for edge is not yet defined! '',/)')
        stop
      endif
      if (iv1.ge.iv2) then
        ind = iv1*(iv1 - 1)/2 + iv2
      else
        ind = iv2*(iv2 - 1)/2 + iv1
      endif
      lhkl(ind) = etmp
    endif
  enddo
!
!  End of input file has been reached
!
10 continue
!###################################################################################################
!  Initial output
!###################################################################################################
  write(6,'(/,'' Energy of Morphology : '',/)')
  write(6,'('' Input lattice vectors : '',/)')
  write(6,'(3f12.6)') (rv(j,1),j=1,3)
  write(6,'(3f12.6)') (rv(j,2),j=1,3)
  write(6,'(3f12.6)') (rv(j,3),j=1,3)
  write(6,'(/,'' Number of surface planes = '',i4,/)') nhkl
  write(6,'('' Surface planes : '',/)')
  write(6,'('' No. :    h   k   l   :     Energy per area'')')
  do i = 1,nhkl
    ehkl(i) = ehkl(i)*unitarea
    write(6,'(i4,'' : '',3i4,''   :'',f12.6)') i,(ihkl(j,i),j=1,3),ehkl(i)
  enddo
  write(6,'(/,'' Distance scale factor = '',f12.6,/)') dscale
!###################################################################################################
!  Post-read checks
!###################################################################################################
!
!  Have cell vectors been read in?
!
  if (.not.linputrv) then
    write(6,'(/,'' Error : No lattice vectors have been read in! '',/)')
    stop
  endif
!
!  Is the number of surface planes sufficient to define a solid?
!
  if (nhkl.lt.6) then
    write(6,'(/,'' Error : Insufficient planes to define a polyhedron! '',/)')
    stop
  endif
!###################################################################################################
!  Convert bulk energy per cell to bulk energy per volume
!###################################################################################################
  volume = rv(1,1)*(rv(2,2)*rv(3,3) - rv(3,2)*rv(2,3)) + &
           rv(2,1)*(rv(1,3)*rv(3,2) - rv(3,3)*rv(1,2)) + &
           rv(3,1)*(rv(1,2)*rv(2,3) - rv(2,2)*rv(1,3))
  volume = abs(volume)
  if (volume.lt.tolerance) then
    write(6,'(/,'' Error : Volume of unit cell is zero! '',/)')
    stop
  endif
  bulkenergy = bulkenergy/volume
!###################################################################################################
!  Set up vectors for surfaces
!###################################################################################################
  do i = 1,nhkl
!
!  Is this a special case - check number of zeros
!
    nzero = 0
    do j = 1,3
      if (ihkl(j,i).eq.0) nzero = nzero + 1
    enddo
!
!  If there are 3 zeros then this is wrong!
!
    if (nzero.eq.3) then
      write(6,'(/,'' Error : Invalid Miller indices! '',/)')
      stop
    elseif (nzero.eq.0) then
!
!  General case
!
      do j = 1,3
        pintersect(1:3,j) = rv(1:3,j)/dble(ihkl(j,i))
      enddo
      vhkl(1:3,1,i) = pintersect(1:3,2) - pintersect(1:3,1)
      vhkl(1:3,2,i) = pintersect(1:3,3) - pintersect(1:3,1)
    elseif (nzero.eq.2) then
!
!  Parallel to 2 axes
!
      if (ihkl(1,i).ne.0) then
        if (ihkl(1,i).gt.0) then
          vhkl(1:3,1,i) = rv(1:3,2)
          vhkl(1:3,2,i) = rv(1:3,3)
        else
          vhkl(1:3,1,i) = rv(1:3,3)
          vhkl(1:3,2,i) = rv(1:3,2)
        endif
      elseif (ihkl(2,i).ne.0) then
        if (ihkl(2,i).gt.0) then
          vhkl(1:3,1,i) = rv(1:3,3)
          vhkl(1:3,2,i) = rv(1:3,1)
        else
          vhkl(1:3,1,i) = rv(1:3,1)
          vhkl(1:3,2,i) = rv(1:3,3)
        endif
      elseif (ihkl(3,i).ne.0) then
        if (ihkl(3,i).gt.0) then
          vhkl(1:3,1,i) = rv(1:3,1)
          vhkl(1:3,2,i) = rv(1:3,2)
        else
          vhkl(1:3,1,i) = rv(1:3,2)
          vhkl(1:3,2,i) = rv(1:3,1)
        endif
      endif
    elseif (nzero.eq.1) then
!
!  Parallel to 1 axis
!
      if (ihkl(1,i).eq.0) then
        if (ihkl(2,i)*ihkl(3,i).gt.0) then
          pintersect(1:3,1) = rv(1:3,2)/dble(ihkl(2,i))
          pintersect(1:3,2) = rv(1:3,3)/dble(ihkl(3,i))
          vhkl(1:3,1,i) = pintersect(1:3,2) - pintersect(1:3,1)
          vhkl(1:3,2,i) = rv(1:3,1)
        else
          pintersect(1:3,1) = rv(1:3,2)/dble(ihkl(2,i))
          pintersect(1:3,2) = rv(1:3,3)/dble(ihkl(3,i))
          vhkl(1:3,1,i) = rv(1:3,1)
          vhkl(1:3,2,i) = pintersect(1:3,2) - pintersect(1:3,1)
        endif
      elseif (ihkl(2,i).eq.0) then
        if (ihkl(1,i)*ihkl(3,i).gt.0) then
          pintersect(1:3,1) = rv(1:3,1)/dble(ihkl(1,i))
          pintersect(1:3,2) = rv(1:3,3)/dble(ihkl(3,i))
          vhkl(1:3,1,i) = rv(1:3,2)
          vhkl(1:3,2,i) = pintersect(1:3,2) - pintersect(1:3,1)
        else
          pintersect(1:3,1) = rv(1:3,1)/dble(ihkl(1,i))
          pintersect(1:3,2) = rv(1:3,3)/dble(ihkl(3,i))
          vhkl(1:3,1,i) = pintersect(1:3,2) - pintersect(1:3,1)
          vhkl(1:3,2,i) = rv(1:3,2)
        endif
      elseif (ihkl(3,i).eq.0) then
        if (ihkl(1,i)*ihkl(2,i).gt.0) then
          pintersect(1:3,1) = rv(1:3,1)/dble(ihkl(1,i))
          pintersect(1:3,2) = rv(1:3,2)/dble(ihkl(2,i))
          vhkl(1:3,1,i) = pintersect(1:3,2) - pintersect(1:3,1)
          vhkl(1:3,2,i) = rv(1:3,3)
        else
          pintersect(1:3,1) = rv(1:3,1)/dble(ihkl(1,i))
          pintersect(1:3,2) = rv(1:3,2)/dble(ihkl(2,i))
          vhkl(1:3,1,i) = rv(1:3,3)
          vhkl(1:3,2,i) = pintersect(1:3,2) - pintersect(1:3,1)
        endif
      endif
    endif
  enddo
!###################################################################################################
!  Set up surface normals
!###################################################################################################
  do i = 1,nhkl
!
!  Take cross product of vectors in the plane of the surface
!
    shkl(1,i) = vhkl(2,1,i)*vhkl(3,2,i) - vhkl(3,1,i)*vhkl(2,2,i)
    shkl(2,i) = vhkl(3,1,i)*vhkl(1,2,i) - vhkl(1,1,i)*vhkl(3,2,i)
    shkl(3,i) = vhkl(1,1,i)*vhkl(2,2,i) - vhkl(2,1,i)*vhkl(1,2,i)
!
!  Normalise vector
!
    snorm = shkl(1,i)**2 + shkl(2,i)**2 + shkl(3,i)**2
    snorm = sqrt(snorm)
    snhkl(1,i) = shkl(1,i)/snorm
    snhkl(2,i) = shkl(2,i)/snorm
    snhkl(3,i) = shkl(3,i)/snorm
!
!  Scale vector according to crystal size factor : dscale * Ehkl
!
    shkl(1,i) = snhkl(1,i)*dscale*ehkl(i)/unitarea
    shkl(2,i) = snhkl(2,i)*dscale*ehkl(i)/unitarea
    shkl(3,i) = snhkl(3,i)*dscale*ehkl(i)/unitarea
!
!  Compute squared magnitude
!
    s2hkl(i) = sqrt(shkl(1,i)**2 + shkl(2,i)**2 + shkl(3,i)**2)
  enddo
  write(6,'('' Surface normals: '',/)')
  write(6,'('' No. :     Surface plane       :      Surface normal vector'')')
  write(6,'(''     :       i    j    k       :        x       y       z   '')')
  do i = 1,nhkl
    write(6,'(i4,'' :   '',3i5,''       :    '',3f8.3)') i,(ihkl(j,i),j=1,3),(shkl(k,i),k=1,3)
  enddo
  write(6,'(/)')
!###################################################################################################
!  Compute initial vertices based on the intersection of 3 planes
!###################################################################################################
!
!  Loop over triplets of planes
!
  do i = 3,nhkl
    do j = 2,i-1
!
!  Compute the cross product of the surface normals of first 2 
!
      cp(1) = shkl(2,i)*shkl(3,j) - shkl(3,i)*shkl(2,j)
      cp(2) = shkl(3,i)*shkl(1,j) - shkl(1,i)*shkl(3,j)
      cp(3) = shkl(1,i)*shkl(2,j) - shkl(2,i)*shkl(1,j)
      crossp = cp(1)**2 + cp(2)**2 + cp(3)**2
      if (abs(crossp).gt.tolerance) then
!
!  If cross product is non-zero then the planes will intersect in a line
!
        do k = 1,j-1
!
!  Compute the dot product of the cross product for i and j with k
!
          dotp = cp(1)*shkl(1,k) + cp(2)*shkl(2,k) + cp(3)*shkl(3,k)
          if (abs(dotp).gt.tolerance) then
!
!  This is a valid set of planes for a vertex so store info about vertex
!
            nvertex = nvertex + 1
            if (nvertex.gt.maxvertex) then
              write(6,'(/,'' Error : Number of vertices exceeds maxvertex! '',/)')
              stop
            endif
            npavertex(nvertex) = 3
            ivertex(1,nvertex) = i
            ivertex(2,nvertex) = j
            ivertex(3,nvertex) = k
!
!  Form matrix of equations of plane (surface normal)
!
            pmat(1:3,1) = shkl(1:3,i)
            pmat(1:3,2) = shkl(1:3,j)
            pmat(1:3,3) = shkl(1:3,k)
!
!  Form right hand side - dot product of surface normal with coordinates of point in plane
!  NB: For the present special case this is simplified since all points are relative to the origin
!
            rhs(1) = (shkl(1,i)**2 + shkl(2,i)**2 + shkl(3,i)**2)
            rhs(2) = (shkl(1,j)**2 + shkl(2,j)**2 + shkl(3,j)**2)
            rhs(3) = (shkl(1,k)**2 + shkl(2,k)**2 + shkl(3,k)**2)
!
!  Solve for coefficients in x, y and z
!
            call matrix_inversion(pmat,i3,i3,wrk,info)
            if (info.ne.0) then
              write(6,'(/,'' Error : Vertex solve has failed! '',i4,/)') info
              stop
            endif
!
!  Assign coordinates of vertex
!
            vertex(1:3,nvertex) = 0.0d0
            do l = 1,3
              vertex(1,nvertex) = vertex(1,nvertex) + pmat(l,1)*rhs(l)
              vertex(2,nvertex) = vertex(2,nvertex) + pmat(l,2)*rhs(l)
              vertex(3,nvertex) = vertex(3,nvertex) + pmat(l,3)*rhs(l)
            enddo
!
!  Before we accept the vertex, does the line from the vertex to the origin cross another plane?
!
            lvertexok = .true.
            l = 0
            do while (l.lt.nhkl.and.lvertexok)
              l = l + 1
!
!  Exclude current planes
!
              if (i.ne.l.and.j.ne.l.and.k.ne.l) then
!
!  Compute the dot product of the vertex position with the surface normal
!
                dotp = snhkl(1,l)*vertex(1,nvertex) + &
                       snhkl(2,l)*vertex(2,nvertex) + &
                       snhkl(3,l)*vertex(3,nvertex)
!
!  If the dot product is greater than the length of the surface normal then exclude vertex
!
                lvertexok = ((dotp-s2hkl(l)).le.tolerance)
              endif
            enddo
!
!  If vertex is not OK then remove
!
            if (.not.lvertexok) nvertex = nvertex - 1
          endif
        enddo
      endif
    enddo
  enddo
!###################################################################################################
!  Find reduced set of vertices where more than 3 planes meet
!###################################################################################################
  lkeepvertex(1:nvertex) = .true.
  do i = 1,nvertex-1
!
!  Search over other vertices looking for a match
!
    do j = i+1,nvertex
      if (lkeepvertex(j)) then
!
!  Compute distance between vertices
!
        dvert = (vertex(1,j) - vertex(1,i))**2 + &
                (vertex(2,j) - vertex(2,i))**2 + &
                (vertex(3,j) - vertex(3,i))**2
        if (dvert.lt.tolerance) then
!
!  Merge vertices
!
          lkeepvertex(j) = .false.
          do k = 1,npavertex(j)
!
!  Is this plane already listed at the vertex?
!
            iv = ivertex(k,j)
            il = 0
            l = 0
            do while (l.lt.npavertex(i).and.il.eq.0)
              l = l + 1
              if (ivertex(l,i).eq.iv) il = l
            enddo
            if (il.eq.0) then
!
!  Plane was not found and so needed to be added
!
              npavertex(i) = npavertex(i) + 1
              if (npavertex(i).gt.maxpav) then
                write(6,'(/,'' Error : Number of planes at a vertex has exceeded maxpav! '',/)')
                stop
              endif
              ivertex(npavertex(i),i) = iv
            endif
          enddo
        endif
      endif
    enddo
  enddo
!
!  Loop over vertices reducing the set
!
  il = nvertex
  nvertex = 0
  do i = 1,il
    if (lkeepvertex(i)) then
      nvertex = nvertex + 1
      npavertex(nvertex) = npavertex(i)
      ivertex(1:npavertex(i),nvertex) = ivertex(1:npavertex(i),i)
      vertex(1:npavertex(i),nvertex)  = vertex(1:npavertex(i),i)
    endif
  enddo
!
  write(6,'(/,'' Number of vertices = '',i4,/)') nvertex
  write(6,'('' Vertices: '',/)')
  write(6,'('' No. :   No. of planes   :    Coordinates of vertices   :  List of planes '')')
  write(6,'(''     :                   :      x       y       z       :  '')')
  do i = 1,nvertex
    write(6,'(i4,'' :     '',i5,''         :  '',3f8.3,''    : '',8i4)') &
      i,npavertex(i),(vertex(k,i),k=1,3),(ivertex(k,i),k=1,npavertex(i))
  enddo
  write(6,'(/)')
!###################################################################################################
!  Compute number of faces
!###################################################################################################
!
!  Loop over vertices to associate them with planes
!
  nhklvertex(1:nhkl) = 0
  do i = 1,nvertex
    do j = 1,npavertex(i)
      nhklvertex(ivertex(j,i)) = nhklvertex(ivertex(j,i)) + 1
      ihklvertex(nhklvertex(ivertex(j,i)),ivertex(j,i)) = i
    enddo
  enddo
!
!  Find number of planes that appear in the morphology
!
  nfaces = 0
  do i = 1,nhkl
    if (nhklvertex(i).gt.0) then
      nfaces = nfaces + 1
      nfaceptr(nfaces) = i
      nfacevertex(nfaces) = nhklvertex(i)
      ifacevertex(1:nhklvertex(i),nfaces) = ihklvertex(1:nhklvertex(i),i)
    endif
  enddo
  write(6,'(/,'' Number of visible surface planes = '',i4,/)') nfaces
  write(6,'('' Visible surface planes: '',/)')
  write(6,'('' No. :   No. of plane    :    Number of vertices '')')
  write(6,'(''     :                   :     '')')
  do i = 1,nfaces
    write(6,'(i4,'' :     '',i5,''         :     '',i4)') &
      i,nfaceptr(i),nfacevertex(i)
  enddo
  write(6,'(/)')
!###################################################################################################
!  Compute volume and surface area
!###################################################################################################
  allocate(vmid2v(3,maxvertex))
  allocate(dotf(maxvertex))
  allocate(ifvordered(maxvertex,nfaces))
  allocate(lfvdone(maxvertex))
  allocate(lface(maxvertex,nfaces))
  allocate(nadjfaceptr(maxvertex,nfaces))
!
  volume = 0.0d0
!
!  Loop over faces
!
  do i = 1,nfaces
!
!  Compute vectors from mid point to vertices
!
    do j = 1,nfacevertex(i)
      iv = ifacevertex(j,i)
      vmid2v(1,j) = vertex(1,iv) - shkl(1,nfaceptr(i))
      vmid2v(2,j) = vertex(2,iv) - shkl(2,nfaceptr(i))
      vmid2v(3,j) = vertex(3,iv) - shkl(3,nfaceptr(i))
    enddo
    iv = 1
    ifvordered(1,i) = 1
    lfvdone(1) = .true.
    lfvdone(2:nfacevertex(i)) = .false.
    do j = 2,nfacevertex(i)
!
!  Loop over vertices finding closest to current reference
!
      vnorm1 = vmid2v(1,iv)*vmid2v(1,iv) + &
               vmid2v(2,iv)*vmid2v(2,iv) + &
               vmid2v(3,iv)*vmid2v(3,iv)
      vnorm1 = sqrt(vnorm1)
      do k = 1,nfacevertex(i)
        dotf(k) = vmid2v(1,k)*vmid2v(1,iv) + &
                  vmid2v(2,k)*vmid2v(2,iv) + &
                  vmid2v(3,k)*vmid2v(3,iv)
        vnormk = vmid2v(1,k)*vmid2v(1,k) + &
                 vmid2v(2,k)*vmid2v(2,k) + &
                 vmid2v(3,k)*vmid2v(3,k)
        vnormk = sqrt(vnormk)
        dotf(k) = dotf(k)/(vnormk*vnorm1)
      enddo
      iv = -1
      dotp = -2.0d0
      do k = 1,nfacevertex(i)
        if (.not.lfvdone(k)) then
          if (dotf(k).gt.dotp) then
            iv = k
            dotp = dotf(k)
          endif
        endif
      enddo
      lfvdone(iv) = .true.
      ifvordered(j,i) = iv
    enddo
!
!  Use pairs of vertices in order to compute the area of the surface, volume and line energy
!
    aface(i) = 0.0d0
    lface(1:nfacevertex(i),i) = 0.0d0
    do j = 1,nfacevertex(i)
      iv1 = ifvordered(j,i)
      if (j.eq.nfacevertex(i)) then
        iv2 = ifvordered(1,i)
      else
        iv2 = ifvordered(j+1,i)
      endif
      vnorm1 = (vmid2v(2,iv1)*vmid2v(3,iv2) - vmid2v(3,iv1)*vmid2v(2,iv2))**2 + &
               (vmid2v(1,iv1)*vmid2v(3,iv2) - vmid2v(3,iv1)*vmid2v(1,iv2))**2 + &
               (vmid2v(2,iv1)*vmid2v(1,iv2) - vmid2v(1,iv1)*vmid2v(2,iv2))**2
!
!  Area
!
      aface(i) = aface(i) + 0.5d0*sqrt(vnorm1)
!
!  Volume
!
      volume = volume + 0.5d0*sqrt(vnorm1)*s2hkl(nfaceptr(i))
!
!  Length
!
      lface(j,i) = (vmid2v(1,iv2) - vmid2v(1,iv1))**2 + &
                   (vmid2v(2,iv2) - vmid2v(2,iv1))**2 + &
                   (vmid2v(3,iv2) - vmid2v(3,iv1))**2
      lface(j,i) = sqrt(lface(j,i))
!
!  Find common plane between vertices for edge (excluding the current face)
!
      ivv1 = ifacevertex(iv1,i)
      ivv2 = ifacevertex(iv2,i)
      lfound1 = .false.
      k = 0
      do while (k.lt.npavertex(ivv1).and..not.lfound1) 
        k = k + 1
        if (ivertex(k,ivv1).ne.nfaceptr(i)) then
          l = 0
          do while (l.lt.npavertex(ivv2).and..not.lfound1) 
            l = l + 1
            if (ivertex(l,ivv2).ne.nfaceptr(i)) then
              lfound1 = (ivertex(k,ivv1).eq.ivertex(l,ivv2))
            endif
          enddo
        endif
      enddo
      if (.not.lfound1) then
        write(6,'(/,'' Error : Valid pair of faces not found for edge! '',/)')
        stop
      endif
      nadjfaceptr(j,i) = ivertex(k,ivv1)
    enddo
!
  enddo
!
  volume = volume/3.0d0
!###################################################################################################
!  Compute energetics
!###################################################################################################
!###################################################################################################
!  Output volume, surface area and energies
!###################################################################################################
  energy = 0.0d0
  write(6,'('' Area and energies of surface planes: '',/)')
  write(6,'('' No. :    h   k   l   :     Area     :      Energy '')')
  do i = 1,nfaces
    iv = nfaceptr(i)
    write(6,'(i4,'' : '',3i4,''   :'',f12.6,''  : '',f12.6)') &
      i,(ihkl(j,iv),j=1,3),aface(i),ehkl(iv)*aface(i)
    energy = energy + ehkl(iv)*aface(i)
  enddo
!
  senergy = energy
!
  write(6,'(/,'' Edge lengths and energies of surface planes: '',/)')
  write(6,'('' No. :    h   k   l   :    h   k   l   :     Length   :      Energy '')')
  do i = 1,nfaces
    iv = nfaceptr(i)
    do j = 1,nfacevertex(i)
      iv2 = nadjfaceptr(j,i)
      if (iv1.ge.iv2) then
        ind = iv1*(iv1-1)/2 + iv2
      else
        ind = iv2*(iv2-1)/2 + iv1
      endif
      write(6,'(i4,'' : '',3i4,''   : '',3i4,''   :'',f12.6,''  : '',f12.6)') &
        i,(ihkl(k,iv),k=1,3),(ihkl(l,iv2),l=1,3),lface(j,i),lhkl(ind)*lface(j,i)
!
!  Edge energy is multiplied by half as each edge will be counted as part of two faces
!
      energy = energy + 0.5d0*lhkl(ind)*lface(j,i)
    enddo
  enddo
!
  lenergy = energy - senergy
!
  write(6,'(/,'' Volume of polyhedron        = '',f12.6,/)') volume
  write(6,'(/,'' Surface energy contribution = '',f12.6)') senergy
  write(6,'(/,'' Line    energy contribution = '',f12.6)') lenergy
  write(6,'(/,'' Bulk    energy contribution = '',f12.6)') volume*bulkenergy
  energy = energy + volume*bulkenergy
  write(6,'(/,'' Total energy of polyhedron  = '',f12.6,/)') energy

end program emorph
