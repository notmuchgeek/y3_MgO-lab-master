  subroutine setspatial(lusefrac)
!
!  Sets up the spatial decomposition of the atoms into boxes
!
!   1/03 Created
!   1/03 Modified so that buffer region is explicitly stored
!   1/03 Fix added to finding of images in the box
!   5/03 Coordinates placed into module for generality
!   5/03 Modified so that coordinate check uses ge/lt
!   5/03 spmin for periodic directions shift to a small negative amount
!  10/03 Size of buffer region generalised for benefit of Brenner pots
!   3/04 Increment for memory allocation increased to save CPU
!   9/04 Buffer for Brenner estimated in more sophisticated way
!  11/04 x/y/zdivision inverted and multiplied for efficiency
!   4/07 Coordinates now generated with mod'ing to ensure that atoms
!        are within the central cell.
!   5/07 lusefrac logical added to control whether coordinates are
!        fractional or Cartesian based
!   6/07 Structure of arrays for storing distribution changed to 1-D
!   6/07 Total number of atoms over all spatial partions added
!   6/07 Problem found for case where angles aren't 90 degrees since
!        when the cell is strained to change orientation of cell. 
!        Spatial decomposition turned off for this case if cell is
!        being varied.
!   6/07 Parallel check added to ensure that there are enough cells 
!        with atoms to make spatial decomposition work.
!   4/08 Option to add input domain size added 
!   4/08 xvec1cell replaced by xvec2cell etc
!   4/08 Buffer size reduced to zero in non-periodic directions
!   3/09 -1 no longer subtracted from nearestx/y/z
!   3/09 Check on cutoff versus supercell size removed
!   3/09 Anisotropic decomposition added
!   7/09 cutoffmax now passed via general module
!   3/10 Verbose output format corrected for lower dimensional cases
!   7/11 nspcell2atptr added to point from atom to cell index
!   6/15 ncellsearch safety margin increased from 1 to 2
!   6/15 Tolerance on cell angles lowered from 10^-3 to 10^-6
!   9/15 Default domain size now set to 55% of rcut since smaller
!        values than typical cutoff sizes are often more efficient
!  10/17 Modified to ensure that Cartesian coordinates are returned to 
!        inside the box
!   1/18 Arrays added to keep track of shifts due to moving into the box
!   7/19 Modified for general cell angles
!   7/19 Correction for lower dimensionalities to fractions in box
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!  10/19 Trap added for nearest values going to zero
!
!  On entry :
!
!  rcut          = maximum potential cut-off radius
!
!  On exit :
!
!  lspatialok    = if .true. then the cell can be decomposed
!                  into boxes
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
!  Julian Gale, CIC, Curtin University, October 2019
!
  use brennerdata, only : bR2, bTR2
  use g_constants, only : degtorad
  use control,     only : keyword, lbrenner, lmodco
  use current
  use general,     only : rcut=>cutoffmax
  use iochannels
  use parallel,    only : ioproc, nprocs
  use spatial
  use symmetry,    only : lra, lsymopt
  implicit none
!
!  Passed variables
!
  logical,     intent(in)  :: lusefrac
!
!  Local variables
!
  integer(i4)              :: i
  integer(i4)              :: icx
  integer(i4)              :: icy
  integer(i4)              :: icz
  integer(i4)              :: ii
  integer(i4)              :: ind
  integer(i4)              :: maxxy
  integer(i4)              :: maxx
  integer(i4)              :: n
  integer(i4)              :: n2bufferplus1x
  integer(i4)              :: n2bufferplus1y
  integer(i4)              :: n2bufferplus1z
  integer(i4)              :: nbuffersizex
  integer(i4)              :: nbuffersizey
  integer(i4)              :: nbuffersizez
  integer(i4)              :: nearestx
  integer(i4)              :: nearesty
  integer(i4)              :: nearestz
  integer(i4)              :: np
  integer(i4)              :: npc
  integer(i4)              :: nx
  integer(i4)              :: ny
  integer(i4)              :: nz
  integer(i4)              :: non90
  logical                  :: linboxx
  logical                  :: linboxy
  logical                  :: linboxz
  logical                  :: lmodcosave
  logical                  :: lortho
  logical                  :: loutbox
  real(dp)                 :: bcutoffmax
  real(dp)                 :: rcutloc
  real(dp)                 :: rcutlocx
  real(dp)                 :: rcutlocy
  real(dp)                 :: rcutlocz
  real(dp)                 :: rshortx
  real(dp)                 :: rshorty
  real(dp)                 :: rshortz
  real(dp)                 :: small
  real(dp)                 :: xdivision
  real(dp)                 :: ydivision
  real(dp)                 :: zdivision
  real(dp)                 :: xi
  real(dp)                 :: yi
  real(dp)                 :: zi
  real(dp)                 :: xf
  real(dp)                 :: yf
  real(dp)                 :: zf
  real(dp)                 :: xfi
  real(dp)                 :: yfi
  real(dp)                 :: zfi
!
  small = 0.000001_dp
!
!  Decide whether cell is suitable for division into cells
!
  lspatialok = .true.
  if (lra) then
    lortho = .true.
  else
    non90 = 0
    if (ndim.eq.3) then
      if (abs(alpha - 90.0_dp).gt.1.0d-6) non90 = non90 + 1
      if (abs(beta  - 90.0_dp).gt.1.0d-6) non90 = non90 + 1
      if (abs(gamma - 90.0_dp).gt.1.0d-6) non90 = non90 + 1
    elseif (ndim.eq.2) then
      if (abs(alpha - 90.0_dp).gt.1.0d-6) non90 = non90 + 1
    endif
!
!  Set flag as to whether orthogonal shortcuts can be used
!
    lortho = (non90.eq.0) 
  endif
!
!  Check that cut-off is greater than zero
!
  if (rcut.lt.1.0d-10) lspatialok = .false.
!
!  If cell is not compatible return
!
  if (.not.lspatialok) return
!
  if (lusefrac) then
!
!  Generate atomic coordinates in box coordinate arrays - here
!  we have to use the mod function to avoid problems.
!
    if (ndim.eq.3) then
      do i = 1,numat
        xi = mod(xfrac(i)+100.0_dp,1.0_dp)
        yi = mod(yfrac(i)+100.0_dp,1.0_dp)
        zi = mod(zfrac(i)+100.0_dp,1.0_dp)
        icosxsp(i) = nint(xi - xfrac(i))
        icosysp(i) = nint(yi - yfrac(i))
        icoszsp(i) = nint(zi - zfrac(i))
        xinbox(i) = xi*r1x + yi*r2x + zi*r3x
        yinbox(i) = xi*r1y + yi*r2y + zi*r3y
        zinbox(i) = xi*r1z + yi*r2z + zi*r3z
      enddo
    elseif (ndim.eq.2) then
      do i = 1,numat
        xi = mod(xfrac(i)+100.0_dp,1.0_dp)
        yi = mod(yfrac(i)+100.0_dp,1.0_dp)
        icosxsp(i) = nint(xi - xfrac(i))
        icosysp(i) = nint(yi - yfrac(i))
        icoszsp(i) = 0
        xinbox(i) = xi*r1x + yi*r2x
        yinbox(i) = xi*r1y + yi*r2y
        zinbox(i) = zfrac(i)
      enddo
    elseif (ndim.eq.1) then
      do i = 1,numat
        xi = mod(xfrac(i)+100.0_dp,1.0_dp)
        icosxsp(i) = nint(xi - xfrac(i))
        icosysp(i) = 0
        icoszsp(i) = 0
        xinbox(i) = xi*r1x
        yinbox(i) = yfrac(i)
        zinbox(i) = zfrac(i)
      enddo
    else
      do i = 1,numat
        xinbox(i) = xfrac(i)
        yinbox(i) = yfrac(i)
        zinbox(i) = zfrac(i)
        icosxsp(i) = 0
        icosysp(i) = 0
        icoszsp(i) = 0
      enddo
    endif
  else
!
!  Use Cartesian coordinates but ensure that they are placed back into the cell
!
!  Save status of lmodco and force to be true for use here
!
    lmodcosave = lmodco
    lmodco = .true.
!
!  Convert to fractional and then back to Cartesian
!
    if (ndim.eq.3) then
      do i = 1,numat
        call cart2frac(ndim,xclat(i),yclat(i),zclat(i),rv,xf,yf,zf,icx,icy,icz)
        xinbox(i) = xf*r1x + yf*r2x + zf*r3x
        yinbox(i) = xf*r1y + yf*r2y + zf*r3y
        zinbox(i) = xf*r1z + yf*r2z + zf*r3z
        icosxsp(i) = icx
        icosysp(i) = icy
        icoszsp(i) = icz
      enddo
    elseif (ndim.eq.2) then
      do i = 1,numat
        call cart2frac(ndim,xclat(i),yclat(i),zclat(i),rv,xf,yf,zf,icx,icy,icz)
        xinbox(i) = xf*r1x + yf*r2x
        yinbox(i) = xf*r1y + yf*r2y
        zinbox(i) = zf
        icosxsp(i) = icx
        icosysp(i) = icy
        icoszsp(i) = icz
      enddo
    elseif (ndim.eq.1) then
      do i = 1,numat
        call cart2frac(ndim,xclat(i),yclat(i),zclat(i),rv,xf,yf,zf,icx,icy,icz)
        xinbox(i) = xf*r1x
        yinbox(i) = yf
        zinbox(i) = zf
        icosxsp(i) = icx
        icosysp(i) = icy
        icoszsp(i) = icz
      enddo
    else
      do i = 1,numat
        xinbox(i) = xclat(i)
        yinbox(i) = yclat(i)
        zinbox(i) = zclat(i)
        icosxsp(i) = 0
        icosysp(i) = 0
        icoszsp(i) = 0
      enddo
    endif
!
!  Reset lmodco
!
    lmodco = lmodcosave
  endif
!
!  Find cell lengths in each direction 
!
  if (ndim.eq.3) then
    spmin(1:3) = - small
    spcell(1) = a
    spcell(2) = b
    spcell(3) = c
  elseif (ndim.eq.2) then
    spmin(1:2) = - small
    spcell(1) = a
    spcell(2) = b
    spmin(3)  = zinbox(1)
    spcell(3) = zinbox(1)
    do i = 1,numat
      spmin(3)  = min(spmin(3),zinbox(i))
      spcell(3) = max(spcell(3),zinbox(i))
    enddo
    spcell(3) = spcell(3) - spmin(3) + small
  elseif (ndim.eq.1) then
    spmin(1)  = - small
    spcell(1) = a
    spmin(2)  = yinbox(1)
    spcell(2) = yinbox(1)
    spmin(3)  = zinbox(1)
    spcell(3) = zinbox(1)
    do i = 1,numat
      spmin(2)  = min(spmin(2),yinbox(i))
      spcell(2) = max(spcell(2),yinbox(i))
      spmin(3)  = min(spmin(3),zinbox(i))
      spcell(3) = max(spcell(3),zinbox(i))
    enddo
    spcell(2) = spcell(2) - spmin(2) + small
    spcell(3) = spcell(3) - spmin(3) + small
  else
    spmin(1)  = xinbox(1)
    spcell(1) = xinbox(1)
    spmin(2)  = yinbox(1)
    spcell(2) = yinbox(1)
    spmin(3)  = zinbox(1)
    spcell(3) = zinbox(1)
    do i = 1,numat
      spmin(1)  = min(spmin(1),xinbox(i))
      spcell(1) = max(spcell(1),xinbox(i))
      spmin(2)  = min(spmin(2),yinbox(i))
      spcell(2) = max(spcell(2),yinbox(i))
      spmin(3)  = min(spmin(3),zinbox(i))
      spcell(3) = max(spcell(3),zinbox(i))
    enddo
    spcell(1) = spcell(1) - spmin(1) + small
    spcell(2) = spcell(2) - spmin(2) + small
    spcell(3) = spcell(3) - spmin(3) + small
  endif
!
!  Set fractional coordinates in the box
!
  if (lusefrac) then
    if (ndim.eq.3) then
      do i = 1,numat
        xi = mod(xfrac(i)+100.0_dp,1.0_dp)
        yi = mod(yfrac(i)+100.0_dp,1.0_dp)
        zi = mod(zfrac(i)+100.0_dp,1.0_dp)
        xfinbox(i) = xi
        yfinbox(i) = yi
        zfinbox(i) = zi
      enddo
    elseif (ndim.eq.2) then
      do i = 1,numat
        xi = mod(xfrac(i)+100.0_dp,1.0_dp)
        yi = mod(yfrac(i)+100.0_dp,1.0_dp)
        xfinbox(i) = xi
        yfinbox(i) = yi
        zfinbox(i) = (zfrac(i) - spmin(3))/spcell(3)
      enddo
    elseif (ndim.eq.1) then
      do i = 1,numat
        xi = mod(xfrac(i)+100.0_dp,1.0_dp)
        xfinbox(i) = xi
        yfinbox(i) = (yfrac(i) - spmin(2))/spcell(2)
        zfinbox(i) = (zfrac(i) - spmin(3))/spcell(3)
      enddo
    else
      do i = 1,numat
        xfinbox(i) = (xfrac(i) - spmin(1))/spcell(1)
        yfinbox(i) = (yfrac(i) - spmin(2))/spcell(2)
        zfinbox(i) = (zfrac(i) - spmin(3))/spcell(3)
      enddo
    endif
  else
!
!  Use Cartesian coordinates but ensure that they are placed back into the cell
!
!  Save status of lmodco and force to be true for use here
!
    lmodcosave = lmodco
    lmodco = .true.
!
!  Convert to fractional and then back to Cartesian
!
    if (ndim.eq.3) then
      do i = 1,numat
        call cart2frac(ndim,xclat(i),yclat(i),zclat(i),rv,xf,yf,zf,icx,icy,icz)
        xfinbox(i) = xf
        yfinbox(i) = yf
        zfinbox(i) = zf
      enddo
    elseif (ndim.eq.2) then
      do i = 1,numat
        call cart2frac(ndim,xclat(i),yclat(i),zclat(i),rv,xf,yf,zf,icx,icy,icz)
        xfinbox(i) = xf
        yfinbox(i) = yf
        zfinbox(i) = (zclat(i) - spmin(3))/spcell(3)
      enddo
    elseif (ndim.eq.1) then
      do i = 1,numat
        call cart2frac(ndim,xclat(i),yclat(i),zclat(i),rv,xf,yf,zf,icx,icy,icz)
        xfinbox(i) = xf
        yfinbox(i) = (yclat(i) - spmin(2))/spcell(2)
        zfinbox(i) = (zclat(i) - spmin(3))/spcell(3)
      enddo
    else
      do i = 1,numat
        xfinbox(i) = (xclat(i) - spmin(1))/spcell(1)
        yfinbox(i) = (yclat(i) - spmin(2))/spcell(2)
        zfinbox(i) = (zclat(i) - spmin(3))/spcell(3)
      enddo
    endif
!
!  Reset lmodco
!
    lmodco = lmodcosave
  endif
!
!  Set local desired domain size based either on cutoff or input value
!
  if (lrcspatial_anisotropic) then
    if (rcspatialx.gt.0.0_dp) then
      rcutlocx = min(rcut,rcspatialx)
    else
      rcutlocx = 0.55_dp*rcut
    endif
    if (rcspatialy.gt.0.0_dp) then
      rcutlocy = min(rcut,rcspatialy)
    else
      rcutlocy = 0.55_dp*rcut
    endif
    if (rcspatialz.gt.0.0_dp) then
      rcutlocz = min(rcut,rcspatialz)
    else
      rcutlocz = 0.55_dp*rcut
    endif
  else
    if (rcspatial.gt.0.0_dp) then
      rcutloc = min(rcut,rcspatial)
      rcutlocx = rcutloc
      rcutlocy = rcutloc
      rcutlocz = rcutloc
    else
      rcutloc = 0.55_dp*rcut
      rcutlocx = 0.55_dp*rcut
      rcutlocy = 0.55_dp*rcut
      rcutlocz = 0.55_dp*rcut
    endif
  endif
!
!  Now find value that will divide cell lengths to give an integer number
!
  nearestx = spcell(1)/rcutlocx
  nearesty = spcell(2)/rcutlocy
  nearestz = spcell(3)/rcutlocz
!
!  Prevent nearest value being zero
!
  nearestx = max(nearestx,1_i4)
  nearesty = max(nearesty,1_i4)
  nearestz = max(nearestz,1_i4)
!
  rnearestx = spcell(1)/dble(nearestx)
  rnearesty = spcell(2)/dble(nearesty)
  rnearestz = spcell(3)/dble(nearestz)
  if (lortho) then
! 
!  Safety margin increased from 1 to 2 as this give problems when atoms are near the cell edge
!
    ncellsearch(1) = 2 + rcut/rnearestx
    ncellsearch(2) = 2 + rcut/rnearesty
    ncellsearch(3) = 2 + rcut/rnearestz
  else
!
!  For non-orthogonal cells find the project of each cell length onto the the other to determine buffer sizes
!
    if (ndim.eq.3) then
      rshortx = rnearestx*min(sin(beta*degtorad),sin(gamma*degtorad))
      rshorty = rnearesty*min(sin(alpha*degtorad),sin(gamma*degtorad))
      rshortz = rnearestz*min(sin(alpha*degtorad),sin(beta*degtorad))
      ncellsearch(1) = 2 + rcut/rshortx
      ncellsearch(2) = 2 + rcut/rshorty
      ncellsearch(3) = 2 + rcut/rshortz
    elseif (ndim.eq.2) then
      ncellsearch(1) = 2 + rcut/(rnearestx*sin(alpha*degtorad))
      ncellsearch(2) = 2 + rcut/(rnearesty*sin(alpha*degtorad))
      ncellsearch(3) = 2 + rcut/rnearestz
    endif
  endif
!
!  Set buffer region size
!
  if (ndim.eq.3) then
    nbufferx = ncellsearch(1)
    nbuffery = ncellsearch(2)
    nbufferz = ncellsearch(3)
  elseif (ndim.eq.2) then
    nbufferx = ncellsearch(1)
    nbuffery = ncellsearch(2)
    nbufferz = 0
  elseif (ndim.eq.1) then
    nbufferx = ncellsearch(1)
    nbuffery = 0
    nbufferz = 0
  elseif (ndim.eq.0) then
    nbufferx = 0
    nbuffery = 0
    nbufferz = 0
  endif
  if (lbrenner.and.nprocs.gt.1) then
    bcutoffmax = max(bR2(1),bR2(2),bR2(3),bTR2(1),bTR2(2),bTR2(3))
    nbuffersizex = 3*bcutoffmax/rcutlocx
    nbuffersizey = 3*bcutoffmax/rcutlocy
    nbuffersizez = 3*bcutoffmax/rcutlocz
    if (ndim.eq.3) then
      nbufferx = nbufferx + nbuffersizex
      nbuffery = nbuffery + nbuffersizey
      nbufferz = nbufferz + nbuffersizez
    elseif (ndim.eq.2) then
      nbufferx = nbufferx + nbuffersizex
      nbuffery = nbuffery + nbuffersizey
    elseif (ndim.eq.1) then
      nbufferx = nbufferx + nbuffersizex
    endif
  endif
!
!  Find numbers of cells in each direction - subtract a small amount for safety
!  Add 2*nbuffer to number of cells in each direction to allow for buffer on each side
!
  nspcell(1) = 2*nbufferx + (spcell(1) - 1.0d-6)/rnearestx + 1
  nspcell(2) = 2*nbuffery + (spcell(2) - 1.0d-6)/rnearesty + 1
  nspcell(3) = 2*nbufferz + (spcell(3) - 1.0d-6)/rnearestz + 1
!
!  Check that minimum number of cells in any direction is at least 2*nbuffer + 1
!
  n2bufferplus1x = 2*nbufferx + 1
  n2bufferplus1y = 2*nbuffery + 1
  n2bufferplus1z = 2*nbufferz + 1
!
  nspcell(1) = max(nspcell(1),n2bufferplus1x)
  nspcell(2) = max(nspcell(2),n2bufferplus1y)
  nspcell(3) = max(nspcell(3),n2bufferplus1z)
!
!  Check that it is viable to use the spatial option - must be more than n2bufferplus1 cells in periodic directions
!
  if (ndim.eq.3) then
    lspatialok = (nspcell(1).gt.n2bufferplus1x.or.nspcell(2).gt.n2bufferplus1y.or.nspcell(3).gt.n2bufferplus1z) 
  elseif (ndim.eq.2) then
    lspatialok = (nspcell(1).gt.n2bufferplus1x.or.nspcell(2).gt.n2bufferplus1y) 
  elseif (ndim.eq.1) then
    lspatialok = (nspcell(1).gt.n2bufferplus1x)
  endif
!
!  Parallel check that there is enough work to do
!
  if (nprocs.gt.1) then
    if ((nspcell(1)-2*nbufferx)*(nspcell(2)-2*nbuffery)*(nspcell(3)-2*nbufferz).lt.nprocs) lspatialok = .false.
  endif
  if (.not.lspatialok) return
!
!  Place atoms within box
!
  do i = 1,numat
    linboxx = (xfinbox(i).ge.0.0_dp.and.xfinbox(i).lt.1.0_dp)
    linboxy = (yfinbox(i).ge.0.0_dp.and.yfinbox(i).lt.1.0_dp)
    linboxz = (zfinbox(i).ge.0.0_dp.and.zfinbox(i).lt.1.0_dp)
    loutbox = (.not.linboxx.or..not.linboxy.or..not.linboxz) 
    if (loutbox.and.ndim.eq.0) then
      call outerror('could not find image in box in setspatial',0_i4)
      call stopnow('setspatial')
    endif
    if (loutbox) then
!
!  Find image in the box
!
      ii = 0
      do while (loutbox.and.ii.lt.iimax2)
        ii = ii + 1
        xfi = xfinbox(i) + dble(ivec2cell(1,ii))
        yfi = yfinbox(i) + dble(ivec2cell(2,ii))
        zfi = zfinbox(i) + dble(ivec2cell(3,ii))
        linboxx = (xfi.ge.0.0_dp.and.xfi.lt.1.0_dp)
        linboxy = (yfi.ge.0.0_dp.and.yfi.lt.1.0_dp)
        linboxz = (zfi.ge.0.0_dp.and.zfi.lt.1.0_dp)
        loutbox = (.not.linboxx.or..not.linboxy.or..not.linboxz) 
      enddo
      if (loutbox) then
        call outerror('could not find image in box in setspatial',0_i4)
        call stopnow('setspatial')
      else
        xi = xinbox(i) + xvec2cell(ii)
        yi = yinbox(i) + yvec2cell(ii)
        zi = zinbox(i) + zvec2cell(ii)
!
        xinbox(i) = xi
        yinbox(i) = yi
        zinbox(i) = zi
!
        xfinbox(i) = xfi
        yfinbox(i) = yfi
        zfinbox(i) = zfi
!
!  Correct cell shift indices for addition of cell vectors
!
        icosxsp(i) = icosxsp(i) + ivec2cell(1,ii)
        icosysp(i) = icosysp(i) + ivec2cell(2,ii)
        icoszsp(i) = icoszsp(i) + ivec2cell(3,ii)
! 
        if (lsymopt) then
          if (nrela2f(nrelf2a(i)).eq.i) then
            xalat(i) = xi
            yalat(i) = yi
            zalat(i) = zi
          endif
        endif
      endif
    endif
  enddo
!
!  Calculate the total number of cells and check that arrays are allocated to this size
!
  nspcelltot = nspcell(1)*nspcell(2)*nspcell(3)
  if (nspcelltot.gt.maxspcell) then
    maxspcell = nspcelltot
    call changemaxspcell
  endif
!
!  Initialise atom counters
!
  nspcellat(1:nspcelltot) = 0
!
!  Find which box each atom belongs in
!
!  Have to loop over cell images in order to fill buffer regions too
!
  xdivision = dble(nearestx)
  ydivision = dble(nearesty)
  zdivision = dble(nearestz)
!
  maxxy = nspcell(1)*nspcell(2)
  maxx  = nspcell(1)
!  
!  First pass to find numbers of atoms per partition
!
  do i = 1,numat
    do ii = 1,iimax2
      nx = (xfinbox(i) + dble(ivec2cell(1,ii)))*xdivision + nbufferx + 1
      if (nx.ge.1.and.nx.le.nspcell(1)) then
        ny = (yfinbox(i) + dble(ivec2cell(2,ii)))*ydivision + nbuffery + 1
        if (ny.ge.1.and.ny.le.nspcell(2)) then
          nz = (zfinbox(i) + dble(ivec2cell(3,ii)))*zdivision + nbufferz + 1
          if (nz.ge.1.and.nz.le.nspcell(3)) then
            ind = (nz-1)*maxxy + (ny-1)*maxx + nx
            nspcellat(ind) = nspcellat(ind) + 1
          endif
        endif
      endif
    enddo
  enddo
!  
!  Construct pointer to start of each partition and count total number of atoms
! 
  ind = 0
  nspcellattot = 0
  do i = 1,nspcelltot
    nspcellat1ptr(i) = ind
    ind = ind + nspcellat(i)
  enddo
  nspcellattot = ind
!
!  Check dimension of arrays for nspcellattot
!
  if (nspcellattot.gt.maxnspcellattot) then
    maxnspcellattot = nspcellattot + 10
    call changemaxnspcellattot
  endif
!
!  Second pass to populate with data
!
  nspcellat(1:nspcelltot) = 0
  do i = 1,numat
    do ii = 1,iimax2
      nx = (xfinbox(i) + dble(ivec2cell(1,ii)))*xdivision + nbufferx + 1
      if (nx.ge.1.and.nx.le.nspcell(1)) then
        ny = (yfinbox(i) + dble(ivec2cell(2,ii)))*ydivision + nbuffery + 1
        if (ny.ge.1.and.ny.le.nspcell(2)) then
          nz = (zfinbox(i) + dble(ivec2cell(3,ii)))*zdivision + nbufferz + 1
          if (nz.ge.1.and.nz.le.nspcell(3)) then
            ind = (nz-1)*maxxy + (ny-1)*maxx + nx
            nspcellat(ind) = nspcellat(ind) + 1
            nspcell2atptr(i) = ind
            nspcellatptr(nspcellat1ptr(ind)+nspcellat(ind)) = i
            nspcellatptrcell(nspcellat1ptr(ind)+nspcellat(ind)) = ii
          endif
        endif
      endif
    enddo
  enddo
!
!  Debugging information
!
  if (ioproc.and.index(keyword,'debu').ne.0) then
    write(ioout,'(/,''  Spatial decomposition data : General potentials : '',/)')
    write(ioout,'(''  Cutoff = '',f8.4,'' Angstroms'',/)') rcut
    if (lrcspatial_anisotropic) then
      write(ioout,'(''  Size input  X = '',f8.4,'' Angstroms'')') rcutlocx
      write(ioout,'(''              Y = '',f8.4,'' Angstroms'')') rcutlocy
      write(ioout,'(''              Z = '',f8.4,'' Angstroms'',/)') rcutlocz
    else
      write(ioout,'(''  Size input    = '',f8.4,'' Angstroms'',/)') rcutloc
    endif
    if (ndim.ge.1) then
      write(ioout,'(''  Size actual x = '',f8.4,'' Angstroms'')') rnearestx
    endif
    if (ndim.ge.2) then
      write(ioout,'(''  Size actual y = '',f8.4,'' Angstroms'')') rnearesty
    endif
    if (ndim.eq.3) then
      write(ioout,'(''  Size actual z = '',f8.4,'' Angstroms'')') rnearestz
    endif
    write(ioout,'(/,''  Spatial decomposition cells (containing atoms) :'',/)')
    write(ioout,'(6x,''Ncells'',4x,''Cell min'',4x,''Cell length'')')
    write(ioout,'(''  x : '',i6,2(3x,f9.4))') nspcell(1),spmin(1),spcell(1)
    write(ioout,'(''  y : '',i6,2(3x,f9.4))') nspcell(2),spmin(2),spcell(2)
    write(ioout,'(''  z : '',i6,2(3x,f9.4))') nspcell(3),spmin(3),spcell(3)
    write(ioout,'(/)')
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    write(ioout,'('' Cell   No. of atoms    Atom No.      x            y            z'')')
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    do i = 1,nspcelltot
      ind = nspcellat1ptr(i)
      if (nspcellat(i).gt.0) then
        np = nspcellatptr(ind+1)
        npc = nspcellatptrcell(ind+1)
        xi = xinbox(np) + xvec2cell(npc)
        yi = yinbox(np) + yvec2cell(npc)
        zi = zinbox(np) + zvec2cell(npc)
        write(ioout,'(i6,4x,i8,4x,i8,3(1x,f12.4))') i,nspcellat(i),np,xi,yi,zi
        do n = 2,nspcellat(i)
          np = nspcellatptr(ind+n)
          npc = nspcellatptrcell(ind+n)
          xi = xinbox(np) + xvec2cell(npc)
          yi = yinbox(np) + yvec2cell(npc)
          zi = zinbox(np) + zvec2cell(npc)
          write(ioout,'(22x,i8,3(1x,f12.4))') np,xi,yi,zi
        enddo
        write(ioout,'(''--------------------------------------------------------------------------------'')')
      endif
    enddo
    write(ioout,'(/)')
  endif
!
  return
  end
