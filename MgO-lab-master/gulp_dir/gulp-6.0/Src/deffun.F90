  subroutine deffun(iflag,n,xc,fc,gc)
!
!  Supplies the function and derivatives for defects
!
!   6/95 Modified to allow for additive defect constraints
!   1/18 Trace added
!   3/19 x0 removed
!   3/19 Change of idopt to idoptindex and idopttype
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
!  Julian Gale, CIC, Curtin University, March 2019
!
  use control
  use current
  use defects
  use derivatives
  use general
  use optimisation
#ifdef TRACE
  use trace,       only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(inout) :: iflag
  integer(i4), intent(in)    :: n
  real(dp),    intent(out)   :: fc
  real(dp),    intent(in)    :: xc(*)
  real(dp),    intent(out)   :: gc(*)
!
!  Local variables
!
  integer(i4)                :: i
  integer(i4)                :: ind
  integer(i4)                :: indf
  integer(i4)                :: indv
  integer(i4)                :: j
  integer(i4)                :: n3r1
  integer(i4)                :: neq
  integer(i4)                :: nf
  integer(i4)                :: nr
  integer(i4)                :: nv
  logical                    :: lfound
  logical                    :: lgrad1
  logical                    :: lgrad2
  real(dp)                   :: csft(3)
  real(dp)                   :: g_cpu_time
  real(dp)                   :: t1
  real(dp)                   :: t2
  real(dp),             save :: tdmax = 0.0_dp
  real(dp)                   :: x0nf
#ifdef TRACE
  call trace_in('deffun')
#endif
!
  t1 = g_cpu_time()
  lgrad1 = (iflag.ge.1)
  lgrad2 = (iflag.ge.2)
!
!  First substitute parameters into place
!
  call defvartocfg(n,xc)
!
!  Now apply constraints
!
  if (ndcon.gt.0) then
    csft(1) = xdc
    csft(2) = ydc
    csft(3) = zdc
    do i = 1,ndcon
      nf = ncdfixind(i)
!
!  Symmetry adapt index
!
      if (ldsym) nf = ndsptr(nf)
!
      if (ncdfixtyp(i).eq.idopt_dx) then
        xdefe(nf) = csft(1)
      elseif (ncdfixtyp(i).eq.idopt_dy) then
        ydefe(nf) = csft(2)
      elseif (ncdfixtyp(i).eq.idopt_dz) then
        zdefe(nf) = csft(3)
      endif
    enddo
    do i = 1,ndcon
      nf = ncdfixind(i)
!
!  Symmetry adapt index
!
      if (ldsym) nf = ndsptr(nf)
!
      if (ncdfixtyp(i).eq.idopt_dx) then
        x0nf = xdefe(nf)
      elseif (ncdfixtyp(i).eq.idopt_dy) then
        x0nf = ydefe(nf)
      elseif (ncdfixtyp(i).eq.idopt_dz) then
        x0nf = zdefe(nf)
      endif
!
      nv = ncdvarind(i)
!
!  Symmetry adapt index
!
      if (ldsym) nv = ndsptr(nv)
!
      if (ncdvartyp(i).eq.idopt_dx) then
        x0nf = (xdefe(nv) - csft(1))*dconco(i) + x0nf
      elseif (ncdvartyp(i).eq.idopt_dy) then
        x0nf = (ydefe(nv) - csft(2))*dconco(i) + x0nf
      elseif (ncdvartyp(i).eq.idopt_dz) then
        x0nf = (zdefe(nv) - csft(3))*dconco(i) + x0nf
      endif
!
      if (ncdfixtyp(i).eq.idopt_dx) then
        xdefe(nf) = x0nf
      elseif (ncdfixtyp(i).eq.idopt_dy) then
        ydefe(nf) = x0nf
      elseif (ncdfixtyp(i).eq.idopt_dz) then
        zdefe(nf) = x0nf
      endif
    enddo
  endif
!
!  Apply symmetry if necessary
!
  if (ldsym) then
    call defequ
  endif
!
  lfirst = .true.
!********************************************
!  Evaluate function and first derivatives  *
!********************************************
  call defener(fc,lgrad1,lgrad2)
  if (lgrad1) then
!***************************************************************
!  Collect derivatives from xdrv,ydrv,zdrv and raderv into gc  *
!***************************************************************
    if (ldsym.and.((lgrad2.and..not.ld2sym).or..not.ld1sym)) then
!
!  If lgrad2, gradients were generated for full cell so they
!  must now be symmetrised
!
      do i = 1,ndasym
        nr = ndsptr(i)
        neq = ndeqv(i)
        xdrv(i) = neq*xdrv(nr)
        ydrv(i) = neq*ydrv(nr)
        zdrv(i) = neq*zdrv(nr)
        if (ldefbsmat(nr)) raderv(i) = neq*raderv(nr)
      enddo
    endif
    if (ldsym) then
      n3r1 = 3*ndasym
    else
      n3r1 = 3*nreg1
    endif
    do i = 1,n
      ind = idoptindex(i)
      if (idopttype(i).eq.idopt_dradius) then
!
!  Radial derviative
!
        gc(i) = raderv(ind)
      elseif (idopttype(i).eq.idopt_dx) then
!
!  Cartesian x derivative
!
        gc(i) = xdrv(ind)
      elseif (idopttype(i).eq.idopt_dy) then
!
!  Cartesian y derivative
!
        gc(i) = ydrv(ind)
      elseif (idopttype(i).eq.idopt_dz) then
!
!  Cartesian z derivative
!
        gc(i) = zdrv(ind)
      endif
    enddo
!****************************
!  Constrained derivatives  *
!****************************
    if (ndcon.gt.0) then
      do i = 1,ndcon
        indf = ncdfixind(i)
        indv = ncdvarind(i)
        if (ncdfixtyp(i).eq.idopt_dradius) then
!
!  Radial derivatives
!
          lfound = .false.
          j = 0
          do while (.not.lfound.and.j.lt.n)
            j = j + 1
            if (indv.eq.idoptindex(j).and.ncdvartyp(i).eq.idopttype(j)) lfound = .true.
          enddo
          gc(j) = gc(j) + raderv(indf)*dconco(i)
        elseif (ncdfixtyp(i).eq.idopt_dx) then
!
!  Cartesian x derivatives
!
          lfound = .false.
          j = 0
          do while (.not.lfound.and.j.lt.n)
            j = j + 1
            if (indv.eq.idoptindex(j).and.ncdvartyp(i).eq.idopttype(j)) lfound = .true.
          enddo
          gc(j) = gc(j) + xdrv(indf)*dconco(i)
        elseif (ncdfixtyp(i).eq.idopt_dy) then
!
!  Cartesian y derivatives
!
          lfound = .false.
          j = 0
          do while (.not.lfound.and.j.lt.n)
            j = j + 1
            if (indv.eq.idoptindex(j).and.ncdvartyp(i).eq.idopttype(j)) lfound = .true.
          enddo
          gc(j) = gc(j) + ydrv(indf)*dconco(i)
        elseif (ncdfixtyp(i).eq.idopt_dz) then
!
!  Cartesian z derivatives
!
          lfound = .false.
          j = 0
          do while (.not.lfound.and.j.lt.n)
            j = j + 1
            if (indv.eq.idoptindex(j).and.ncdvartyp(i).eq.idopttype(j)) lfound = .true.
          enddo
          gc(j) = gc(j) + zdrv(indf)*dconco(i)
        endif
      enddo
    endif
  endif
  t2 = g_cpu_time()
  if (t2-t1.gt.tdmax) tdmax = t2 - t1
  if (timmax.gt.0.0_dp) then
    if ((timmax-t2).lt.tdmax.and..not.lrelax) iflag = - 1
  endif
#ifdef TRACE
  call trace_out('deffun')
#endif
!
  return
  end
