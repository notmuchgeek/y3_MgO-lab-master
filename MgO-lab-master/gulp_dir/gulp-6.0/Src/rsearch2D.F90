  subroutine rsearch2D(xcrd,ycrd,zcrd,lmolok,lcspair,i,j,nmi,ixj,iyj,izj,nor,nmolonly,lself,cut2)
!
!  Subroutine that searches for valid distances - called from real
!  space subroutines.
!
!  On input :
!
!  xcrd    = x-coordinate of i-j vector for central unit cell
!  ycrd    = y-coordinate of i-j vector for central unit cell
!  zcrd    = z-coordinate of i-j vector for central unit cell
!  lmolok  = whether it is possible for atoms to be in molecule
!  lcspair = is this a possible core-shell pair
!  i       = number of atom i
!  j       = number of atom j
!  nmi     = molecule number for atom i
!  ixj     = molecule cell index in x direction
!  iyj     = molecule cell index in y direction
!  izj     = molecule cell index in z direction
!  cut2    = cut-off distance squared for this pair of atoms
!
!  On exit :
!
!  nor     = number of valid vectors between atoms
!  dist    = distance between atoms i and j for each vector
!  xtmp    = x component of valid i-j vector
!  ytmp    = y component of valid i-j vector
!  ztmp    = z component of valid i-j vector
!  lbonded = logical for bonding/connectivity
!  l2bonds = logical for bonding/connectivity
!  l3bonds = logical for bonding/connectivity
!  lptrmol = logical for bonding/connectivity
!  nmolonly= vector within molecule if non-zero
!  lself   = if .true. then a self term was found and must be handled
!
!
!  l111    = if .true. all interactions lie within unit cell and immediately
!            adjacent cells only => use short cut
!
!   2/01 Created from rsearch3D
!   6/01 Setting of lsamemol altered
!   1/02 Handling of negative cell parameters added
!   1/02 Minimum image option incorporated
!  12/02 Nadd increased for more extreme angles
!   1/03 Wolf modifications made
!   3/03 Call to rfind replaces general case
!   9/03 rfind used in lra case as well
!  11/04 sqrt cut2e removed
!   1/05 rp removed as an argument since it is no longer needed
!   1/05 Storage of cell indices for valid vectors added
!  10/06 Default value of small lowered
!   3/07 Bonding types added
!  11/07 Modifications for reaxFF Coulomb term added
!   1/08 lreaxFFqreal removed
!   1/08 Calls to bonded3c modified
!   2/09 Argument removed from changemaxdis call
!   3/09 small replaced by global value smallself from general module
!   3/12 Position of call to getbonds moved as it is needed for rfind2D
!   2/18 Trace added
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
!  Julian Gale, CIC, Curtin University, Feb 2018
!
  use control,       only : lminimage, lwolf
  use current
  use general
  use kspace
  use molecule
  use realvectors
  use shells
  use symmetry,      only : lra
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  use two
  implicit none
!
!  Passed variables
!
  integer(i4)        :: i   
  integer(i4)        :: ixj 
  integer(i4)        :: iyj  
  integer(i4)        :: izj
  integer(i4)        :: j  
  integer(i4)        :: nmi
  integer(i4)        :: nmolonly
  integer(i4)        :: nor
  logical            :: lcspair
  logical            :: lmolok
  logical            :: lself
  real(dp)           :: cut2 
!
!  Local variables
!
  integer(i4)        :: ii
  integer(i4)        :: iim
  integer(i4)        :: jj
  integer(i4)        :: jjm
  logical            :: l111
  logical            :: lsamemol
  real(dp)           :: cmax
  real(dp)           :: cut2e
  real(dp)           :: cut2s
  real(dp)           :: r2
  real(dp)           :: r2min
  real(dp)           :: xcd
  real(dp)           :: xcd2
  real(dp)           :: ycd     
  real(dp)           :: xcdi
  real(dp)           :: ycdi
  real(dp)           :: xcrd   
  real(dp)           :: ycrd  
  real(dp)           :: zcrd 
  real(dp)           :: zcrd2
  real(dp)           :: xmin
  real(dp)           :: ymin
#ifdef TRACE
  call trace_in('rsearch2D')
#endif
!
!  Initialise return variables
!
  nor = 0
  nmolonly = 0
  lself = .false.
!     
!  For minimum image we must exclude i=j case since there is no unique nearest image
!     
  if (lminimage) then
    if (i.eq.j) then
      lself = .true.
#ifdef TRACE
      call trace_out('rsearch2D')
#endif
      return   
    endif 
  endif
!***************************
!  Set up local variables  *
!***************************
!
!  Set up cutoffs
!
  if (lwolf) then
    cut2e = cutw*cutw
  else
    cut2e = rmx2
  endif
  cut2s = cuts*cuts
!
!  Decide whether all interactions lie within unit cell and first
!  neighbours - saves time for large systems
!
  cmax = max(cut2,cut2e)
  l111 = .true.
  if (a*a.lt.cmax.or.b*b.lt.cmax) l111 = .false.
  if (alpha.lt.80.0_dp.or.alpha.gt.100.0_dp) l111 = .false.
!********************
!  Looping section  *
!********************
  zcrd2 = zcrd*zcrd
  if (lminimage) then
!***********************
!  Find minimum image  *
!***********************
    if (lra) then
      iim = dnint(xcrd/r1x)
      jjm = dnint(ycrd/r2y)
      xmin = xcrd - r1x*iim
      ymin = ycrd - r2y*jjm
      r2min = xmin*xmin + ymin*ymin + zcrd2
    else
      r2min = 1.0d10
      xcdi = xcrd - 2.0_dp*r1x
      ycdi = ycrd - 2.0_dp*r1y
      do ii = -1,1
        xcdi = xcdi + r1x  
        ycdi = ycdi + r1y
        xcd = xcdi - 2.0_dp*r2x
        ycd = ycdi - 2.0_dp*r2y
        do jj = -1,1
          xcd = xcd + r2x
          ycd = ycd + r2y
          r2 = xcd*xcd + ycd*ycd + zcrd2
          if (r2.lt.r2min) then
            r2min = r2
            xmin = xcd
            ymin = ycd
            iim = ii
            jjm = jj
          endif
        enddo
      enddo 
    endif   
!***************************
!  Molecule - check index  *
!***************************
    if (lmolok.and.(r2min.gt.cut2s.or..not.lcspair)) then
      call bonded3(lbonded(nor+1),l2bonds(nor+1),l3bonds(nor+1),nbotype(nor+1),nbotype2(nor+1),i,j,iim,jjm,0_i4)
      lsamemol = (lbonded(nor+1).or.l2bonds(nor+1).or.l3bonds(nor+1))
      if (.not.lsamemol) then
        call samemol(lsamemol,nmi,iim,jjm,0_i4,ixj,iyj,izj)
      endif
      lptrmol(nor+1) = lsamemol
      if (lsamemol) then
        if (r2min.gt.cut2e) nmolonly = 1
      else  
        lbonded(nor+1) = .false.
        l2bonds(nor+1) = .false.
        l3bonds(nor+1) = .false.
      endif
    else
      lptrmol(nor+1)  = .false.
      lbonded(nor+1)  = .false.
      l2bonds(nor+1)  = .false.
      l3bonds(nor+1)  = .false. 
      nbotype(nor+1)  = 0
      nbotype2(nor+1) = 0
    endif
!***************************
!  Check cut-off distance  *
!***************************
    if (r2min.lt.smallself) then
      lself = .true.
    elseif (r2min.le.cut2.or.lptrmol(nor+1)) then
!           
!  Store vector
!         
      nor = 1
      dist(nor) = r2min
      xtmp(nor) = xmin
      ytmp(nor) = ymin
      ztmp(nor) = zcrd
      cellindex(1,nor) = iim
      cellindex(2,nor) = jjm
      cellindex(3,nor) = 0
    endif
!********************
!  Find all images  *
!********************
  elseif (lra) then
!
!  Set up bonding
!
    if (lmolok) call getbonds(i,j)
!
!  Right angled cell
!
    if (l111) then
!
      xcd = xcrd - 2.0_dp*r1x
      do ii = -1,1
        xcd = xcd + r1x
        xcd2 = xcd*xcd
        ycd = ycrd - 2.0_dp*r2y
        do jj = -1,1
          ycd = ycd + r2y
          r2 = xcd2 + ycd*ycd + zcrd2
          if (nor.ge.maxdis) then
            maxdis = nor + 200
            call changemaxdis
          endif
!
!  Molecule - check index
!
          if (lmolok.and.(r2.gt.cut2s.or..not.lcspair)) then
            call bonded3c(lbonded(nor+1),l2bonds(nor+1),l3bonds(nor+1),nbotype(nor+1),nbotype2(nor+1),ii,jj,0_i4)
            lsamemol = (lbonded(nor+1).or.l2bonds(nor+1).or.l3bonds(nor+1))
            if (.not.lsamemol) then
              call samemol(lsamemol,nmi,ii,jj,0_i4,ixj,iyj,izj)
            endif
            lptrmol(nor+1) = lsamemol
            if (lsamemol) then
              if (r2.gt.cut2e) nmolonly = nor + 1
            else
              lbonded(nor+1) = .false.
              l2bonds(nor+1) = .false.
              l3bonds(nor+1) = .false.
            endif
          else
            lptrmol(nor+1)  = .false.
            lbonded(nor+1)  = .false.
            l2bonds(nor+1)  = .false.
            l3bonds(nor+1)  = .false.
            nbotype(nor+1)  = 0
            nbotype2(nor+1) = 0
          endif
          if (r2.lt.smallself) then
!
!  Self term - set flag to indicate it was found
!
            lself = .true.
          elseif (r2.le.cut2.or.lptrmol(nor+1)) then
!
!  Store vector
!
            nor = nor + 1
            dist(nor) = r2
            xtmp(nor) = xcd
            ytmp(nor) = ycd
            ztmp(nor) = zcrd
            cellindex(1,nor) = ii
            cellindex(2,nor) = jj
            cellindex(3,nor) = 0
          endif
        enddo
      enddo
    else
      call rfind2D(xcrd,ycrd,zcrd,cut2,cut2e,cut2s,lcspair,lmolok,i,j,nmi,ixj,iyj,izj,nmolonly,lself,0_i4,nor)
    endif
  else
!
!  Set up bonding
!
    if (lmolok) call getbonds(i,j)
!
!  General cell
!
    if (l111) then
!
      xcdi = xcrd - 2.0_dp*r1x
      ycdi = ycrd - 2.0_dp*r1y
!
!  Loop over unit cells to find interatomic distances
!
      do ii = - 1,1
        xcdi = xcdi + r1x
        ycdi = ycdi + r1y
!
        xcd = xcdi - 2.0_dp*r2x
        ycd = ycdi - 2.0_dp*r2y
        do jj = - 1,1
          xcd = xcd + r2x
          ycd = ycd + r2y
          r2 = xcd*xcd + ycd*ycd + zcrd2
          if (nor.ge.maxdis) then
            maxdis = nor + 200
            call changemaxdis
          endif
!
!  Molecule - check index
!
          if (lmolok.and.(r2.gt.cut2s.or..not.lcspair)) then
            call bonded3c(lbonded(nor+1),l2bonds(nor+1),l3bonds(nor+1),nbotype(nor+1),nbotype2(nor+1),ii,jj,0_i4)
            lsamemol = (lbonded(nor+1).or.l2bonds(nor+1).or.l3bonds(nor+1))
            if (.not.lsamemol) then
              call samemol(lsamemol,nmi,ii,jj,0_i4,ixj,iyj,izj)
            endif
            lptrmol(nor+1) = lsamemol
            if (lsamemol) then
              if (r2.gt.cut2e) nmolonly = nor + 1
            else
              lbonded(nor+1) = .false.
              l2bonds(nor+1) = .false.
              l3bonds(nor+1) = .false.
            endif
          else
            lptrmol(nor+1)  = .false.
            lbonded(nor+1)  = .false.
            l2bonds(nor+1)  = .false.
            l3bonds(nor+1)  = .false.
            nbotype(nor+1)  = 0
            nbotype2(nor+1) = 0
          endif
          if (r2.lt.smallself) then
!
!  Self term - set flag to indicate it was found
!
            lself = .true.
          elseif (r2.le.cut2.or.lptrmol(nor+1)) then
!
!  Store vector
!
            nor = nor + 1
            dist(nor) = r2
            xtmp(nor) = xcd
            ytmp(nor) = ycd
            ztmp(nor) = zcrd
            cellindex(1,nor) = ii
            cellindex(2,nor) = jj
            cellindex(3,nor) = 0
          endif
        enddo
      enddo
    else
      call rfind2D(xcrd,ycrd,zcrd,cut2,cut2e,cut2s,lcspair,lmolok,i,j,nmi,ixj,iyj,izj,nmolonly,lself,0_i4,nor)
    endif
  endif
#ifdef TRACE
  call trace_out('rsearch2D')
#endif
!
  return
  end
