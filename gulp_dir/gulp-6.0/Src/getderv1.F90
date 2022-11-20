  subroutine getderv1(n,xc,gc,lgrad2,lrotalg)
!
!  Collects the derivatives into a linear array, including
!  some symmetry and constraint handling
!
!  NB: lrotalg not supported for rigid molecules
!
!   8/97 Created from funct
!  10/97 If lrotalg then calculate asymmetric unit derivs
!        by rotation and summation - this is because K
!        point sampling doesn't lead to all derivatives
!        being equal in a free energy gradient calculation
!   5/03 Region 3 modifications added
!   6/03 No of strains corrected to nstrains from 6 in constraint
!        section
!  11/06 Nudging of gradients option added
!   6/09 Charge as a coordinate option added
!   7/13 Bug in handling of constrained breathing shell derivatives fixed
!   3/17 Modifications made to allow for new variable order in iopt
!   2/18 Trace added
!   3/19 iopt replaced by ioptindex and iopttype
!   9/19 Rigid molecules added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!  10/19 Constrained derivatives added for rigid molecules
!  11/19 Rigid molecule handling for periodic systems added
!  12/19 2D translation derivatives for rigid molecules corrected
!  12/19 Correction to referencing of internal derivatives for rigid molecule case
!   3/20 Rigid molecule modifications added
!   5/20 xfdrv, yfdrv, zfdrv added
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
!  Julian Gale, CIC, Curtin University, May 2020
!
  use configurations, only : lbsmat
  use control,        only : lrigid
  use current
  use derivatives
  use molecule
  use optimisation
  use symmetry
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  use xcgc,           only : lnudgegc
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)    :: n
  logical,     intent(in)    :: lgrad2
  logical,     intent(in)    :: lrotalg
  real(dp),    intent(out)   :: gc(*)
  real(dp),    intent(in)    :: xc(*)
!
!  Local variables
!
  integer(i4)                :: i
  integer(i4)                :: indf
  integer(i4)                :: indv
  integer(i4)                :: k
  integer(i4)                :: kk
  integer(i4)                :: j
  integer(i4)                :: neq
  integer(i4)                :: nj
  integer(i4)                :: nr
  integer(i4)                :: nv
  logical                    :: lfound
  real(dp)                   :: cellderv(6)
  real(dp)                   :: g1(3)
  real(dp)                   :: g2(3)
  real(dp)                   :: rneq
  real(dp)                   :: rotinv(3,3)
  real(dp)                   :: rvp(3,3)
  real(dp)                   :: xdv
  real(dp)                   :: ydv
  real(dp)                   :: zdv
#ifdef TRACE
  call trace_in('getderv1')
#endif
!
!  For symmetry optimisations correct internal first derivatives
!
  if (lsymopt.and.(.not.lsymderv.or.(lgrad2.and..not.lsymderv2))) then
!
!  Check for lrotalg and rigid molecules
!
    if (lrigid.and.lrotalg) then
      call outerror('rigid molecules not supported with lrotalg',0_i4)
      call stopnow('getderv1')
    endif
!
!  Generate non-primitive cell
!
    do i = 1,3
      rvp(1,i) = rv(1,i)
      rvp(2,i) = rv(2,i)
      rvp(3,i) = rv(3,i)
    enddo
    if (ncbl.gt.1) call uncentre(rvp)
!
!  Transform cartesian derivatives to internal coords before symmetry reduction so that rotation can be performed
!
    do kk = 1,numat
      xdv = xdrv(kk)
      ydv = ydrv(kk)
      zdv = zdrv(kk)
      xdrv(kk) = xdv*rvp(1,1) + ydv*rvp(2,1) + zdv*rvp(3,1)
      ydrv(kk) = xdv*rvp(1,2) + ydv*rvp(2,2) + zdv*rvp(3,2)
      zdrv(kk) = xdv*rvp(1,3) + ydv*rvp(2,3) + zdv*rvp(3,3)
    enddo
    if (lrigid) then
!
!  Rigid molecule derivatives - translation only
!
      do kk = 1,nmol
        xdv = molTdrv(1,kk)
        ydv = molTdrv(2,kk)
        zdv = molTdrv(3,kk)
        molTdrv(1,kk) = xdv*rvp(1,1) + ydv*rvp(2,1) + zdv*rvp(3,1)
        molTdrv(2,kk) = xdv*rvp(1,2) + ydv*rvp(2,2) + zdv*rvp(3,2)
        molTdrv(3,kk) = xdv*rvp(1,3) + ydv*rvp(2,3) + zdv*rvp(3,3)
      enddo
    endif
!
!  Save full cell derivatives for use in rigidmoleculeprop
!
    xfdrv(1:numat) = xdrv(1:numat)
    yfdrv(1:numat) = ydrv(1:numat)
    zfdrv(1:numat) = zdrv(1:numat)
!
    if (lrotalg) then
!*****************************
!  Rotate and sum algorithm  *
!*****************************
      do i = 1,nasym
        nr = nrela2f(i)
        neq = neqv(i)
        xdrv(i) = xdrv(nr)
        ydrv(i) = ydrv(nr)
        zdrv(i) = zdrv(nr)
        if (neq.gt.1) then
          do j = 2,neq
            nr = nr + 1
            g1(1) = xdrv(nr)
            g1(2) = ydrv(nr)
            g1(3) = zdrv(nr)
            do k = 1,3
              rotinv(1,k) = rop(1,k,nrotop(nr))
              rotinv(2,k) = rop(2,k,nrotop(nr))
              rotinv(3,k) = rop(3,k,nrotop(nr))
            enddo
            g2(1) = rotinv(1,1)*g1(1) + rotinv(2,1)*g1(2) + rotinv(3,1)*g1(3)
            g2(2) = rotinv(1,2)*g1(1) + rotinv(2,2)*g1(2) + rotinv(3,2)*g1(3)
            g2(3) = rotinv(1,3)*g1(1) + rotinv(2,3)*g1(2) + rotinv(3,3)*g1(3)
            xdrv(i) = xdrv(i) + g2(1)
            ydrv(i) = ydrv(i) + g2(2)
            zdrv(i) = zdrv(i) + g2(3)
          enddo
        endif
      enddo
    else
!***************************
!  Conventional algorithm  *
!***************************
      do i = 1,nasym
        nr = nrela2f(i)
        rneq = dble(neqv(i))
        xdrv(i) = rneq*xdrv(nr)
        ydrv(i) = rneq*ydrv(nr)
        zdrv(i) = rneq*zdrv(nr)
        if (lbsmat(i+nsft)) raderv(i) = rneq*raderv(nr)
      enddo
      if (lrigid) then
!
!  Rigid molecule derivatives
!
        do i = 1,nmolasym
          nr = nmola2f(i)
          rneq = dble(neqv(i))
          molTdrv(1:3,i) = rneq*molTdrv(1:3,nr)
          molQdrv(1:3,i) = rneq*molQdrv(1:3,nr)
        enddo
      endif
    endif
  elseif (ndim.eq.3) then
!
!  Generate non-primitive cell
!
    do i = 1,3
      rvp(1,i) = rv(1,i)
      rvp(2,i) = rv(2,i)
      rvp(3,i) = rv(3,i)
    enddo
    if (ncbl.gt.1) call uncentre(rvp)
!
!  Transform cartesian derivatives to internal coords
!
    do kk = 1,nasym
      xdv = xdrv(kk)
      ydv = ydrv(kk)
      zdv = zdrv(kk)
      xdrv(kk) = xdv*rvp(1,1) + ydv*rvp(2,1) + zdv*rvp(3,1)
      ydrv(kk) = xdv*rvp(1,2) + ydv*rvp(2,2) + zdv*rvp(3,2)
      zdrv(kk) = xdv*rvp(1,3) + ydv*rvp(2,3) + zdv*rvp(3,3)
    enddo
    if (lrigid) then
!
!  Rigid molecule derivatives
!
      do kk = 1,nmol
        xdv = molTdrv(1,kk)
        ydv = molTdrv(2,kk)
        zdv = molTdrv(3,kk)
        molTdrv(1,kk) = xdv*rvp(1,1) + ydv*rvp(2,1) + zdv*rvp(3,1)
        molTdrv(2,kk) = xdv*rvp(1,2) + ydv*rvp(2,2) + zdv*rvp(3,2)
        molTdrv(3,kk) = xdv*rvp(1,3) + ydv*rvp(2,3) + zdv*rvp(3,3)
      enddo
    endif
  elseif (ndim.eq.2) then
!
!  Transform cartesian derivatives to internal coords for x and y
!
    do kk = 1,nasym
      xdv = xdrv(kk)
      ydv = ydrv(kk)
      xdrv(kk) = xdv*rv(1,1) + ydv*rv(2,1)
      ydrv(kk) = xdv*rv(1,2) + ydv*rv(2,2)
    enddo
    if (lrigid) then
!
!  Rigid molecule derivatives
!
      do kk = 1,nmol
        xdv = molTdrv(1,kk)
        ydv = molTdrv(2,kk)
        molTdrv(1,kk) = xdv*rv(1,1) + ydv*rv(2,1)
        molTdrv(2,kk) = xdv*rv(1,2) + ydv*rv(2,2)
      enddo
    endif
  elseif (ndim.eq.1) then
!
!  Transform cartesian derivatives to internal coords for x only
!
    do kk = 1,nasym
      xdv = xdrv(kk)
      xdrv(kk) = xdv*rv(1,1)
    enddo
    if (lrigid) then
!
!  Rigid molecule derivatives
!
      do kk = 1,nmol
        xdv = molTdrv(1,kk)
        molTdrv(1,kk) = xdv*rv(1,1)
      enddo
    endif
  endif
!*********************************
!  Collect internal derivatives  *
!*********************************
  do i = ninternalmin,ninternalmax
    nj = ioptindex(i)
    if (iopttype(i).eq.iopt_xf) then
      gc(i) = xdrv(nasymnomolptr(nj))
    elseif (iopttype(i).eq.iopt_yf) then
      gc(i) = ydrv(nasymnomolptr(nj))
    elseif (iopttype(i).eq.iopt_zf) then
      gc(i) = zdrv(nasymnomolptr(nj))
    elseif (iopttype(i).eq.iopt_radius) then
      gc(i) = raderv(nasymnomolptr(nj))
    elseif (iopttype(i).eq.iopt_xcom) then
      gc(i) = molTdrv(1,nj)
    elseif (iopttype(i).eq.iopt_ycom) then
      gc(i) = molTdrv(2,nj)
    elseif (iopttype(i).eq.iopt_zcom) then
      gc(i) = molTdrv(3,nj)
    elseif (iopttype(i).eq.iopt_xqtn) then
      gc(i) = molQdrv(1,nj)
    elseif (iopttype(i).eq.iopt_yqtn) then
      gc(i) = molQdrv(2,nj)
    elseif (iopttype(i).eq.iopt_zqtn) then
      gc(i) = molQdrv(3,nj)
    endif
  enddo
!****************************
!  Constrained derivatives  *
!****************************
  if (ncon.gt.0) then
    do i = 1,ncon
      indf = ncfixtyp(i)
      indv = ncvartyp(i)
      if (indf.eq.iopt_radius) then
!
!  Radial derivatives
!
        nv = ncfixind(i)
        lfound = .false.
        j = ninternalmax - nbsm
        do while (.not.lfound.and.j.le.ninternalmax)
          j = j + 1
          if (indv.eq.iopttype(j).and.ncvarind(i).eq.ioptindex(j)) lfound = .true.
        enddo
        gc(j) = gc(j) + raderv(nv)*conco(i)
      elseif (indf.eq.iopt_xf) then
!
!  Internal x derivatives
!
        nv = ncfixind(i)
        lfound = .false.
        j = ninternalmin - 1
        do while (.not.lfound.and.j.le.ninternalmax)
          j = j + 1
          if (indv.eq.iopttype(j).and.ncvarind(i).eq.ioptindex(j)) lfound = .true.
        enddo
        gc(j) = gc(j) + xdrv(nv)*conco(i)
      elseif (indf.eq.iopt_yf) then
!
!  Internal y derivatives
!
        nv = ncfixind(i)
        lfound = .false.
        j = ninternalmin - 1
        do while (.not.lfound.and.j.le.ninternalmax)
          j = j + 1
          if (indv.eq.iopttype(j).and.ncvarind(i).eq.ioptindex(j)) lfound = .true.
        enddo
        gc(j) = gc(j) + ydrv(nv)*conco(i)
      elseif (indf.eq.iopt_zf) then
!
!  Internal z derivatives
!
        nv = ncfixind(i)
        lfound = .false.
        j = ninternalmin - 1
        do while (.not.lfound.and.j.le.ninternalmax)
          j = j + 1
          if (indv.eq.iopttype(j).and.ncvarind(i).eq.ioptindex(j)) lfound = .true.
        enddo
        gc(j) = gc(j) + zdrv(nv)*conco(i)
      elseif (indf.eq.iopt_strain.or.indf.eq.iopt_cell) then
!
!  Strain derivatives
!
        nv = ncfixind(i)
        j = ncvarind(i)
        if (indv.eq.iopt_strain.or.indv.eq.iopt_cell) then
          strderv(j) = strderv(j) + strderv(nv)*conco(i)
        endif
      elseif (indf.eq.iopt_xcom) then
!
!  Internal x com derivatives
!
        nv = ncfixind(i)
        lfound = .false.
        j = ninternalmin - 1
        do while (.not.lfound.and.j.le.ninternalmax)
          j = j + 1
          if (indv.eq.iopttype(j).and.ncvarind(i).eq.ioptindex(j)) lfound = .true.
        enddo
        gc(j) = gc(j) + molTdrv(1,nv)*conco(i)
      elseif (indf.eq.iopt_ycom) then
!
!  Internal y com derivatives
!
        nv = ncfixind(i)
        lfound = .false.
        j = ninternalmin - 1
        do while (.not.lfound.and.j.le.ninternalmax)
          j = j + 1
          if (indv.eq.iopttype(j).and.ncvarind(i).eq.ioptindex(j)) lfound = .true.
        enddo
        gc(j) = gc(j) + molTdrv(2,nv)*conco(i)
      elseif (indf.eq.iopt_zcom) then
!
!  Internal z com derivatives
!
        nv = ncfixind(i)
        lfound = .false.
        j = ninternalmin - 1
        do while (.not.lfound.and.j.le.ninternalmax)
          j = j + 1
          if (indv.eq.iopttype(j).and.ncvarind(i).eq.ioptindex(j)) lfound = .true.
        enddo
        gc(j) = gc(j) + molTdrv(3,nv)*conco(i)
      elseif (indf.eq.iopt_xqtn) then
!
!  Internal x qtn derivatives
!
        nv = ncfixind(i)
        lfound = .false.
        j = ninternalmin - 1
        do while (.not.lfound.and.j.le.ninternalmax)
          j = j + 1
          if (indv.eq.iopttype(j).and.ncvarind(i).eq.ioptindex(j)) lfound = .true.
        enddo
        gc(j) = gc(j) + molQdrv(1,nv)*conco(i)
      elseif (indf.eq.iopt_yqtn) then
!
!  Internal y qtn derivatives
!
        nv = ncfixind(i)
        lfound = .false.
        j = ninternalmin - 1
        do while (.not.lfound.and.j.le.ninternalmax)
          j = j + 1
          if (indv.eq.iopttype(j).and.ncvarind(i).eq.ioptindex(j)) lfound = .true.
        enddo
        gc(j) = gc(j) + molQdrv(2,nv)*conco(i)
      elseif (indf.eq.iopt_zqtn) then
!
!  Internal z qtn derivatives
!
        nv = ncfixind(i)
        lfound = .false.
        j = ninternalmin - 1
        do while (.not.lfound.and.j.le.ninternalmax)
          j = j + 1
          if (indv.eq.iopttype(j).and.ncvarind(i).eq.ioptindex(j)) lfound = .true.
        enddo
        gc(j) = gc(j) + molQdrv(3,nv)*conco(i)
      endif
    enddo
  endif
!*******************************
!  Collect strain derivatives  *
!*******************************
  if (ncell.gt.0) then
    if (loptcellpar) then
      call celldrv(strderv,cellderv) 
      do i = ncellmin,ncellmax
        gc(i) = cellderv(ioptindex(i))
      enddo
    else
      do i = ncellmin,ncellmax
        gc(i) = strderv(ioptindex(i))
      enddo
    endif
  endif
!*************************
!  Nudging of gradients  *
!*************************
  if (lnudgegc) then
    call nudge(n,xc,gc)
  endif
#ifdef TRACE
  call trace_out('getderv1')
#endif
!
  return
  end
