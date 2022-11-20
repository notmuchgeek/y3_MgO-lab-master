  subroutine symgenkpt(kptin,maxkptout,nkptout,kptout)
!
!  Generate symmetry equivalent k points from a given k point
!
!  On entry:
!
!  kptin(3) = fractional coordinates of original k point
!  maxkptout = dimension of kptout array
!
!  On exit:
!
!  nkptout  = number of symmetry images of k point
!  kptout(3,nkptout) = fractional coordinates of k points
!
!  11/14 Created from symgenatom
!   2/18 Trace added
!
  use current,        only : ncblp
  use symmetry
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4),   intent(out) :: nkptout
  integer(i4),   intent(in)  :: maxkptout
  real(dp),      intent(in)  :: kptin(3)
  real(dp),      intent(out) :: kptout(3,maxkptout)
!
!  Local variables
!
  integer(i4)                :: mv
  integer(i4)                :: n
  integer(i4)                :: nin
  logical                    :: lfound
  real(dp),             save :: thresh = 0.00001_dp
  real(dp)                   :: x(3)
  real(dp)                   :: v(3)
#ifdef TRACE
  call trace_in('symgenkpt')
#endif
!
  nkptout = 0
!***********************************************************************
!  Loop over symmetry operators to generate symmetry related k points  *
!***********************************************************************
  do mv = 1,ngop
!
!  Roto-translation
!
    x(1) = 0.0_dp
    x(2) = 0.0_dp
    x(3) = 0.0_dp
    v(1) = vitp(1,mv)
    v(2) = vitp(2,mv)
    v(3) = vitp(3,mv)
    call GULP_mxmb(ropp(1,1,mv),1_i4,3_i4,kptin,1_i4,1_i4,v,1_i4,1_i4,3_i4,3_i4,1_i4)
!
!  Transform coordinates to primitive cell
!
    if (ngop.gt.1) then
      x(1) = v(1)
      x(2) = v(2)
      x(3) = v(3)
    else
      call GULP_mxmb(w1p(ncblp,1,1),7_i4,21_i4,v,1_i4,1_i4,x,1_i4,1_i4,3_i4,3_i4,1_i4)
    endif
!
!  Place fractional coordinates in range 0-1
!
    x(1) = x(1) + 3.0_dp
    nin  = x(1)
    x(1) = x(1) - nin
    x(2) = x(2) + 3.0_dp
    nin  = x(2)
    x(2) = x(2) - nin
    x(3) = x(3) + 3.0_dp
    nin  = x(3)
    x(3) = x(3) - nin
!
!  Compare k point with those previously generated
!
    lfound = .false.
    n = 0
    do while (n.lt.nkptout.and..not.lfound)
      n = n + 1
      v(1) = kptout(1,n) - x(1)
      v(2) = kptout(2,n) - x(2)
      v(3) = kptout(3,n) - x(3)
      v(1) = v(1) - nint(v(1))
      v(2) = v(2) - nint(v(2))
      v(3) = v(3) - nint(v(3))
      lfound = ((abs(v(1))+abs(v(2))+abs(v(3))).le.thresh)
    enddo 
!
!  Check array dimensions
!
    if (.not.lfound) then
      nkptout = nkptout + 1
      if (nkptout.le.maxkptout) then
        kptout(1,nkptout) = x(1)
        kptout(2,nkptout) = x(2)
        kptout(3,nkptout) = x(3)
      endif
    endif
  enddo 
#ifdef TRACE
  call trace_out('symgenkpt')
#endif
!
  return
  end
