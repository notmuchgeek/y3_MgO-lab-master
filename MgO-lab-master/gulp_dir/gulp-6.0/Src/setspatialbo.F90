  subroutine setspatialbo(lusefrac)
!
!  Sets up the spatial decomposition of the atoms into boxes.
!  Version for specific bond order decomposition.
!
!   9/04 Created from setspatial
!  10/04 References to brenner potential removed
!   4/07 Coordinates now mod'ed to ensure they are in the
!        central cell.
!   5/07 lusefrac logical added to control whether coordinates are
!        fractional or Cartesian based
!   6/07 Structure of arrays for storing distribution changed to 1-D
!   6/07 Total number of atoms over all spatial partions added
!   6/07 Error messages modified to indicate setspatialbo
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
!   6/15 ncellsearch safety margin increased from 1 to 2
!   6/15 Tolerance on cell angles lowered from 10^-3 to 10^-6
!   9/15 Default domain size now set to 55% of rcut since smaller
!        values than typical cutoff sizes are often more efficient
!  10/17 Modified to ensure that Cartesian coordinates are returned to 
!        inside the box
!   7/19 Modified to allow for general cell angles
!   7/19 Correction for lower dimensionalities to fractions in box
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!  10/19 Trap added for nearest values going to zero
!
!  On entry :
!
!  rcut          = maximum potential cut-off radius
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
  use bondorderdata, only : nbopot
  use g_constants,   only : degtorad
  use control,       only : keyword, lmodco
  use current
  use general,       only : rcut=>cutoffmaxbo
  use iochannels
  use parallel,      only : ioproc, nprocs
  use spatialbo
  use symmetry,      only : lra, lsymopt
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
  logical                  :: linboxbox
  logical                  :: linboxboy
  logical                  :: linboxboz
  logical                  :: lmodcosave
  logical                  :: lortho
  logical                  :: loutbox
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
  lspatialBOok = .true.
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
  if (rcut.lt.1.0d-10) lspatialBOok = .false.
!
!  If cell is not compatible return
!
  if (.not.lspatialBOok) return
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
        xinboxbo(i) = xi*r1x + yi*r2x + zi*r3x
        yinboxbo(i) = xi*r1y + yi*r2y + zi*r3y
        zinboxbo(i) = xi*r1z + yi*r2z + zi*r3z
      enddo
    elseif (ndim.eq.2) then
      do i = 1,numat
        xi = mod(xfrac(i)+100.0_dp,1.0_dp)
        yi = mod(yfrac(i)+100.0_dp,1.0_dp)
        xinboxbo(i) = xi*r1x + yi*r2x
        yinboxbo(i) = xi*r1y + yi*r2y
        zinboxbo(i) = zfrac(i)
      enddo
    elseif (ndim.eq.1) then
      do i = 1,numat
        xi = mod(xfrac(i)+100.0_dp,1.0_dp)
        xinboxbo(i) = xi*r1x
        yinboxbo(i) = yfrac(i)
        zinboxbo(i) = zfrac(i)
      enddo
    else
      do i = 1,numat
        xinboxbo(i) = xfrac(i)
        yinboxbo(i) = yfrac(i)
        zinboxbo(i) = zfrac(i)
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
        xinboxbo(i) = xf*r1x + yf*r2x + zf*r3x
        yinboxbo(i) = xf*r1y + yf*r2y + zf*r3y
        zinboxbo(i) = xf*r1z + yf*r2z + zf*r3z
      enddo
    elseif (ndim.eq.2) then
      do i = 1,numat
        call cart2frac(ndim,xclat(i),yclat(i),zclat(i),rv,xf,yf,zf,icx,icy,icz)
        xinboxbo(i) = xf*r1x + yf*r2x
        yinboxbo(i) = xf*r1y + yf*r2y
        zinboxbo(i) = zf
      enddo
    elseif (ndim.eq.1) then
      do i = 1,numat
        call cart2frac(ndim,xclat(i),yclat(i),zclat(i),rv,xf,yf,zf,icx,icy,icz)
        xinboxbo(i) = xf*r1x
        yinboxbo(i) = yf
        zinboxbo(i) = zf
      enddo
    else
      do i = 1,numat
        xinboxbo(i) = xclat(i)
        yinboxbo(i) = yclat(i)
        zinboxbo(i) = zclat(i)
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
    spminbo(1:3) = - small
    if (abs(alpha - 90.0_dp).gt.1.0d-3) then
      spcellbo(1) = a
      spcellbo(2) = b
      spcellbo(3) = c*sin(alpha*degtorad)
    elseif (abs(beta - 90.0_dp).gt.1.0d-3) then
      spcellbo(1) = a
      spcellbo(2) = b
      spcellbo(3) = c*sin(beta*degtorad)
    elseif (abs(gamma - 90.0_dp).gt.1.0d-3) then
      spcellbo(1) = a
      spcellbo(2) = b*sin(gamma*degtorad)
      spcellbo(3) = c
    else
      spcellbo(1) = a
      spcellbo(2) = b
      spcellbo(3) = c
    endif
  elseif (ndim.eq.2) then
    spminbo(1:2) = - small
    if (abs(alpha - 90.0_dp).gt.1.0d-3) then
      spcellbo(1) = a
      spcellbo(2) = b*sin(alpha*degtorad)
    else
      spcellbo(1) = a
      spcellbo(2) = b
    endif
    spminbo(3)  = zinboxbo(1)
    spcellbo(3) = zinboxbo(1)
    do i = 1,numat
      spminbo(3)  = min(spminbo(3),zinboxbo(i))
      spcellbo(3) = max(spcellbo(3),zinboxbo(i))
    enddo
    spcellbo(3) = spcellbo(3) - spminbo(3) + small
  elseif (ndim.eq.1) then
    spminbo(1)  = - small
    spcellbo(1) = a
    spminbo(2)  = yinboxbo(1)
    spcellbo(2) = yinboxbo(1)
    spminbo(3)  = zinboxbo(1)
    spcellbo(3) = zinboxbo(1)
    do i = 1,numat
      spminbo(2)  = min(spminbo(2),yinboxbo(i))
      spcellbo(2) = max(spcellbo(2),yinboxbo(i))
      spminbo(3)  = min(spminbo(3),zinboxbo(i))
      spcellbo(3) = max(spcellbo(3),zinboxbo(i))
    enddo
    spcellbo(2) = spcellbo(2) - spminbo(2) + small
    spcellbo(3) = spcellbo(3) - spminbo(3) + small
  else
    spminbo(1)  = xinboxbo(1)
    spcellbo(1) = xinboxbo(1)
    spminbo(2)  = yinboxbo(1)
    spcellbo(2) = yinboxbo(1)
    spminbo(3)  = zinboxbo(1)
    spcellbo(3) = zinboxbo(1)
    do i = 1,numat
      spminbo(1)  = min(spminbo(1),xinboxbo(i))
      spcellbo(1) = max(spcellbo(1),xinboxbo(i))
      spminbo(2)  = min(spminbo(2),yinboxbo(i))
      spcellbo(2) = max(spcellbo(2),yinboxbo(i))
      spminbo(3)  = min(spminbo(3),zinboxbo(i))
      spcellbo(3) = max(spcellbo(3),zinboxbo(i))
    enddo
    spcellbo(1) = spcellbo(1) - spminbo(1) + small
    spcellbo(2) = spcellbo(2) - spminbo(2) + small
    spcellbo(3) = spcellbo(3) - spminbo(3) + small
  endif
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
        xfinboxbo(i) = xi
        yfinboxbo(i) = yi
        zfinboxbo(i) = zi
      enddo
    elseif (ndim.eq.2) then
      do i = 1,numat
        xi = mod(xfrac(i)+100.0_dp,1.0_dp)
        yi = mod(yfrac(i)+100.0_dp,1.0_dp)
        xfinboxbo(i) = xi
        yfinboxbo(i) = yi
        zfinboxbo(i) = (zfrac(i) - spminbo(3))/spcellbo(3)
      enddo
    elseif (ndim.eq.1) then
      do i = 1,numat
        xi = mod(xfrac(i)+100.0_dp,1.0_dp)
        xfinboxbo(i) = xi
        yfinboxbo(i) = (yfrac(i) - spminbo(2))/spcellbo(2)
        zfinboxbo(i) = (zfrac(i) - spminbo(3))/spcellbo(3)
      enddo
    else
      do i = 1,numat
        xfinboxbo(i) = (xfrac(i) - spminbo(1))/spcellbo(1)
        yfinboxbo(i) = (yfrac(i) - spminbo(2))/spcellbo(2)
        zfinboxbo(i) = (zfrac(i) - spminbo(3))/spcellbo(3)
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
        xfinboxbo(i) = xf
        yfinboxbo(i) = yf
        zfinboxbo(i) = zf
      enddo
    elseif (ndim.eq.2) then
      do i = 1,numat
        call cart2frac(ndim,xclat(i),yclat(i),zclat(i),rv,xf,yf,zf,icx,icy,icz)
        xfinboxbo(i) = xf
        yfinboxbo(i) = yf
        zfinboxbo(i) = (zclat(i) - spminbo(3))/spcellbo(3)
      enddo
    elseif (ndim.eq.1) then
      do i = 1,numat
        call cart2frac(ndim,xclat(i),yclat(i),zclat(i),rv,xf,yf,zf,icx,icy,icz)
        xfinboxbo(i) = xf
        yfinboxbo(i) = (yclat(i) - spminbo(2))/spcellbo(2)
        zfinboxbo(i) = (zclat(i) - spminbo(3))/spcellbo(3)
      enddo
    else
      do i = 1,numat
        xfinboxbo(i) = (xclat(i) - spminbo(1))/spcellbo(1)
        yfinboxbo(i) = (yclat(i) - spminbo(2))/spcellbo(2)
        zfinboxbo(i) = (zclat(i) - spminbo(3))/spcellbo(3)
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
  if (lrcspatialBO_anisotropic) then
    if (rcspatialbox.gt.0.0_dp) then
      rcutlocx = min(rcut,rcspatialbox)
    else
      rcutlocx = 0.55_dp*rcut
    endif
    if (rcspatialboy.gt.0.0_dp) then
      rcutlocy = min(rcut,rcspatialboy)
    else
      rcutlocy = 0.55_dp*rcut
    endif
    if (rcspatialboz.gt.0.0_dp) then
      rcutlocz = min(rcut,rcspatialboz)
    else
      rcutlocz = 0.55_dp*rcut
    endif
  else
    if (rcspatialbo.gt.0.0_dp) then
      rcutloc = min(rcut,rcspatialbo)
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
  nearestx = spcellbo(1)/rcutlocx
  nearesty = spcellbo(2)/rcutlocy
  nearestz = spcellbo(3)/rcutlocz
!
!  Prevent nearest value being zero
!
  nearestx = max(nearestx,1_i4)
  nearesty = max(nearesty,1_i4)
  nearestz = max(nearestz,1_i4)
!
  rnearestxbo = spcellbo(1)/dble(nearestx)
  rnearestybo = spcellbo(2)/dble(nearesty)
  rnearestzbo = spcellbo(3)/dble(nearestz)
  if (lortho) then
! 
!  Safety margin increased from 1 to 2 as this give problems when atoms are near the cell edge
!
    ncellsearchbo(1) = 2 + rcut/rnearestxbo
    ncellsearchbo(2) = 2 + rcut/rnearestybo
    ncellsearchbo(3) = 2 + rcut/rnearestzbo
  else
!
!  For non-orthogonal cells find the project of each cell length onto the the other to determine buffer sizes
!
    if (ndim.eq.3) then
      rshortx = rnearestxbo*min(sin(beta*degtorad),sin(gamma*degtorad))
      rshorty = rnearestybo*min(sin(alpha*degtorad),sin(gamma*degtorad))
      rshortz = rnearestzbo*min(sin(alpha*degtorad),sin(beta*degtorad))
      ncellsearchbo(1) = 2 + rcut/rshortx
      ncellsearchbo(2) = 2 + rcut/rshorty
      ncellsearchbo(3) = 2 + rcut/rshortz
    elseif (ndim.eq.2) then
      ncellsearchbo(1) = 2 + rcut/(rnearestxbo*sin(alpha*degtorad))
      ncellsearchbo(2) = 2 + rcut/(rnearestybo*sin(alpha*degtorad))
      ncellsearchbo(3) = 2 + rcut/rnearestzbo
    endif
  endif
!       
!  Set buffer region size
!       
  if (nbopot.gt.0.and.nprocs.gt.1) then 
    nbuffersizex = 3.0_dp*rcut/rcutlocx
    nbuffersizey = 3.0_dp*rcut/rcutlocy
    nbuffersizez = 3.0_dp*rcut/rcutlocz
    if (ndim.eq.3) then
      nbufferxbo = max(nbuffersizex,ncellsearchbo(1))
      nbufferybo = max(nbuffersizey,ncellsearchbo(2))
      nbufferzbo = max(nbuffersizez,ncellsearchbo(3))
    elseif (ndim.eq.2) then
      nbufferxbo = max(nbuffersizex,ncellsearchbo(1))
      nbufferybo = max(nbuffersizey,ncellsearchbo(2))
      nbufferzbo = 0
    elseif (ndim.eq.1) then
      nbufferxbo = max(nbuffersizex,ncellsearchbo(1))
      nbufferybo = 0
      nbufferzbo = 0
    else
      nbufferxbo = 0
      nbufferybo = 0
      nbufferzbo = 0
    endif
  else  
    if (ndim.eq.3) then
      nbufferxbo = ncellsearchbo(1)
      nbufferybo = ncellsearchbo(2)
      nbufferzbo = ncellsearchbo(3)
    elseif (ndim.eq.2) then
      nbufferxbo = ncellsearchbo(1)
      nbufferybo = ncellsearchbo(2)
      nbufferzbo = 0
    elseif (ndim.eq.1) then
      nbufferxbo = ncellsearchbo(1)
      nbufferybo = 0
      nbufferzbo = 0
    elseif (ndim.eq.0) then
      nbufferxbo = 0
      nbufferybo = 0
      nbufferzbo = 0
    endif
  endif
!
!  Find numbers of cells in each direction - subtract a small amount for safety
!  Add 2*nbuffer to number of cells in each direction to allow for buffer on each side
!
  nspcellbo(1) = 2*nbufferxbo + (spcellbo(1) - 1.0d-6)/rnearestxbo + 1
  nspcellbo(2) = 2*nbufferybo + (spcellbo(2) - 1.0d-6)/rnearestybo + 1
  nspcellbo(3) = 2*nbufferzbo + (spcellbo(3) - 1.0d-6)/rnearestzbo + 1
!
!  Check that minimum number of cells in any direction is at least 2*nbuffer + 1
!
  n2bufferplus1x = 2*nbufferxbo + 1
  n2bufferplus1y = 2*nbufferybo + 1
  n2bufferplus1z = 2*nbufferzbo + 1
!
  nspcellbo(1) = max(nspcellbo(1),n2bufferplus1x)
  nspcellbo(2) = max(nspcellbo(2),n2bufferplus1y)
  nspcellbo(3) = max(nspcellbo(3),n2bufferplus1z)
!
!  Check that it is viable to use the spatial option - must be more than
!  n2bufferplus1 cells in periodic directions.
!
  if (ndim.eq.3) then
    lspatialBOok = (nspcellbo(1).gt.n2bufferplus1x.or.nspcellbo(2).gt.n2bufferplus1y.or.nspcellbo(3).gt.n2bufferplus1z) 
  elseif (ndim.eq.2) then
    lspatialBOok = (nspcellbo(1).gt.n2bufferplus1x.or.nspcellbo(2).gt.n2bufferplus1y) 
  elseif (ndim.eq.1) then
    lspatialBOok = (nspcellbo(1).gt.n2bufferplus1x)
  endif
!
!  Parallel check that there is enough work to do
!
  if (nprocs.gt.1) then
    if ((nspcellbo(1)-2*nbufferxbo)*(nspcellbo(2)-2*nbufferybo)*(nspcellbo(3)-2*nbufferzbo).lt.nprocs) lspatialBOok = .false.
  endif
  if (.not.lspatialBOok) return
!
!  Place atoms within box
!
  do i = 1,numat
    linboxbox = (xfinboxbo(i).ge.0.0_dp.and.xfinboxbo(i).lt.1.0_dp)
    linboxboy = (yfinboxbo(i).ge.0.0_dp.and.yfinboxbo(i).lt.1.0_dp)
    linboxboz = (zfinboxbo(i).ge.0.0_dp.and.zfinboxbo(i).lt.1.0_dp)
    loutbox = (.not.linboxbox.or..not.linboxboy.or..not.linboxboz) 
    if (loutbox.and.ndim.eq.0) then
      call outerror('could not find image in box in setspatialbo',0_i4)
      call stopnow('setspatial')
    endif
    if (loutbox) then
!
!  Find image in the box
!
      ii = 0
      do while (loutbox.and.ii.lt.iimax2) 
        ii = ii + 1
        xfi = xfinboxbo(i) + dble(ivec2cell(1,ii))
        yfi = yfinboxbo(i) + dble(ivec2cell(2,ii))
        zfi = zfinboxbo(i) + dble(ivec2cell(3,ii))
        linboxbox = (xfi.ge.0.0_dp.and.xfi.lt.1.0_dp)
        linboxboy = (yfi.ge.0.0_dp.and.yfi.lt.1.0_dp)
        linboxboz = (zfi.ge.0.0_dp.and.zfi.lt.1.0_dp)
        loutbox = (.not.linboxbox.or..not.linboxboy.or..not.linboxboz) 
      enddo
      if (loutbox) then
        call outerror('could not find image in box in setspatialbo',0_i4)
        call stopnow('setspatial')
      else
        xi = xinboxbo(i) + xvec2cell(ii)
        yi = yinboxbo(i) + yvec2cell(ii)
        zi = zinboxbo(i) + zvec2cell(ii)
!
        xinboxbo(i) = xi
        yinboxbo(i) = yi
        zinboxbo(i) = zi
!
        xfinboxbo(i) = xfi
        yfinboxbo(i) = yfi
        zfinboxbo(i) = zfi
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
  nspcelltotbo = nspcellbo(1)*nspcellbo(2)*nspcellbo(3)
  if (nspcelltotbo.gt.maxspcellbo) then
    maxspcellbo = nspcelltotbo
    call changemaxspcellbo
  endif
!
!  Initialise atom counters
!
  nspcellatbo(1:nspcelltotbo) = 0
!
!  Find which box each atom belongs in
!
!  Have to loop over cell images in order to fill buffer regions too
!
  xdivision = dble(nearestx)
  ydivision = dble(nearesty)
  zdivision = dble(nearestz)
!
  maxxy = nspcellbo(1)*nspcellbo(2)
  maxx  = nspcellbo(1)
!
!  First pass to find numbers of atoms per partition 
!
  do i = 1,numat
    do ii = 1,iimax2
      nx = (xfinboxbo(i) + dble(ivec2cell(1,ii)))*xdivision + nbufferxbo + 1
      if (nx.ge.1.and.nx.le.nspcellbo(1)) then
        ny = (yfinboxbo(i) + dble(ivec2cell(2,ii)))*ydivision + nbufferybo + 1
        if (ny.ge.1.and.ny.le.nspcellbo(2)) then
          nz = (zfinboxbo(i) + dble(ivec2cell(3,ii)))*zdivision + nbufferzbo + 1
          if (nz.ge.1.and.nz.le.nspcellbo(3)) then
            ind = (nz-1)*maxxy + (ny-1)*maxx + nx
            nspcellatbo(ind) = nspcellatbo(ind) + 1
          endif
        endif
      endif
    enddo
  enddo
!
!  Construct pointer to start of each partition
!
  ind = 0
  do i = 1,nspcelltotbo
    nspcellat1ptrbo(i) = ind
    ind = ind + nspcellatbo(i)
  enddo
  nspcellattotbo = ind
!
!  Check dimension of arrays for nspcellattot
!
  if (nspcellattotbo.gt.maxnspcellattotbo) then
    maxnspcellattotbo = nspcellattotbo
    call changemaxnspcellattotbo
  endif
!
!  Second pass to populate with data
!
  nspcellatbo(1:nspcelltotbo) = 0
  do i = 1,numat
    do ii = 1,iimax2
      nx = (xfinboxbo(i) + dble(ivec2cell(1,ii)))*xdivision + nbufferxbo + 1
      if (nx.ge.1.and.nx.le.nspcellbo(1)) then
        ny = (yfinboxbo(i) + dble(ivec2cell(2,ii)))*ydivision + nbufferybo + 1
        if (ny.ge.1.and.ny.le.nspcellbo(2)) then
          nz = (zfinboxbo(i) + dble(ivec2cell(3,ii)))*zdivision + nbufferzbo + 1
          if (nz.ge.1.and.nz.le.nspcellbo(3)) then
            ind = (nz-1)*maxxy + (ny-1)*maxx + nx
            nspcellatbo(ind) = nspcellatbo(ind) + 1
            nspcellatptrbo(nspcellat1ptrbo(ind)+nspcellatbo(ind)) = i
            nspcellatptrcellbo(nspcellat1ptrbo(ind)+nspcellatbo(ind)) = ii
          endif
        endif
      endif
    enddo
  enddo
!
!  Debugging information
!
  if (ioproc.and.index(keyword,'debu').ne.0) then
    write(ioout,'(/,''  Spatial decomposition data : Bond order potentials :'',/)')
    write(ioout,'(''  Cutoff = '',f8.4,'' Angstroms'',/)') rcut
    if (lrcspatialBO_anisotropic) then
      write(ioout,'(''  Size input  X = '',f8.4,'' Angstroms'')') rcutlocx
      write(ioout,'(''              Y = '',f8.4,'' Angstroms'')') rcutlocy
      write(ioout,'(''              Z = '',f8.4,'' Angstroms'',/)') rcutlocz
    else
      write(ioout,'(''  Size input    = '',f8.4,'' Angstroms'',/)') rcutloc
    endif
    if (ndim.ge.1) then
      write(ioout,'(''  Size actual x = '',f8.4,'' Angstroms'')') rnearestxbo
    endif
    if (ndim.ge.2) then
      write(ioout,'(''  Size actual y = '',f8.4,'' Angstroms'')') rnearestybo
    endif
    if (ndim.eq.3) then
      write(ioout,'(''  Size actual z = '',f8.4,'' Angstroms'')') rnearestzbo
    endif
    write(ioout,'(/,''  Spatial decomposition cells (containing atoms) :'',/)')
    write(ioout,'(6x,''Ncells'',4x,''Cell min'',4x,''Cell length'')')
    write(ioout,'(''  x : '',i6,2(3x,f9.4))') nspcellbo(1),spminbo(1),spcellbo(1)
    write(ioout,'(''  y : '',i6,2(3x,f9.4))') nspcellbo(2),spminbo(2),spcellbo(2)
    write(ioout,'(''  z : '',i6,2(3x,f9.4))') nspcellbo(3),spminbo(3),spcellbo(3)
    write(ioout,'(/)')
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    write(ioout,'('' Cell   No. of atoms    Atom No.      x            y            z'')')
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    do i = 1,nspcelltotbo
      ind = nspcellat1ptrbo(i)
      if (nspcellatbo(i).gt.0) then
        np = nspcellatptrbo(ind+1)
        npc = nspcellatptrcellbo(ind+1)
        xi = xinboxbo(np) + xvec2cell(npc)
        yi = yinboxbo(np) + yvec2cell(npc)
        zi = zinboxbo(np) + zvec2cell(npc)
        write(ioout,'(i6,4x,i8,4x,i8,3(1x,f12.4))') i,nspcellatbo(i),np,xi,yi,zi
        do n = 2,nspcellatbo(i)
          np = nspcellatptrbo(ind+n)
          npc = nspcellatptrcellbo(ind+n)
          xi = xinboxbo(np) + xvec2cell(npc)
          yi = yinboxbo(np) + yvec2cell(npc)
          zi = zinboxbo(np) + zvec2cell(npc)
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
