  subroutine outkptfile
!
!  Generate symmetry equivalent k points from a given k point
!  and output this to a file
!
!   5/15 Created from symgenkpt
!   6/17 Module files renamed to gulp_files
!   6/18 Trap for parallel I/O added
!   8/19 Correction to the syntax of open
!
  use current,      only : ncblp, ncf
  use gulp_files,   only : lkpt, kptfile
  use ksample
  use parallel,     only : ioproc
  use symmetry
  implicit none
!
!  Local variables
!
  integer(i4)                :: i
  integer(i4)                :: iout
  integer(i4)                :: j
  integer(i4)                :: k
  integer(i4)                :: mv
  integer(i4)                :: n
  integer(i4)                :: nin
  integer(i4)                :: nk
  integer(i4)                :: nki
  logical                    :: lfound
  real(dp),             save :: thresh = 0.00001_dp
  real(dp)                   :: kp(3)
  real(dp)                   :: kpall(3,48)
  real(dp)                   :: x(3)
  real(dp)                   :: v(3)
!
!  If not lkpt then return
!
  if (.not.lkpt) return
!
!  Ensure that only a single node writes this file
!
  if (.not.ioproc) return
!
!  Open file
!
  iout = 24
!
!  If name has been given then open file
!
  if (kptfile(1:1).ne.' ') then
    if (ncf.eq.1) then
      open(iout,file=kptfile,status='unknown')
    else
      open(iout,file=kptfile,status='unknown',access='sequential',position='append')
    endif
  endif
!
!  Write out header
!
  write(iout,'(''#  Configuration = '',i6)') ncf
!
!  Loop over k points
!
  nk = 0
  do i = 1,nkpt
    if (nkptcfg(i).eq.ncf) then
      nk = nk + 1
      kp(1) = xkpt(i)
      kp(2) = ykpt(i)
      kp(3) = zkpt(i)
!***********************************************************************
!  Loop over symmetry operators to generate symmetry related k points  *
!***********************************************************************
      nki = 0
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
        call GULP_mxmb(ropp(1,1,mv),1_i4,3_i4,kp,1_i4,1_i4,v,1_i4,1_i4,3_i4,3_i4,1_i4)
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
!  Compare k point with those previously generated from this k point
!
        lfound = .false.
        n = 0
        do while (n.lt.nki.and..not.lfound)
          n = n + 1
          v(1) = kpall(1,n) - x(1)
          v(2) = kpall(2,n) - x(2)
          v(3) = kpall(3,n) - x(3)
          v(1) = v(1) - nint(v(1))
          v(2) = v(2) - nint(v(2))
          v(3) = v(3) - nint(v(3))
          lfound = ((abs(v(1))+abs(v(2))+abs(v(3))).le.thresh)
        enddo 
!
!  No duplicate found, so add k point to the list
!
        if (.not.lfound) then
          nki = nki + 1
          kpall(1,nki) = x(1)
          kpall(2,nki) = x(2)
          kpall(3,nki) = x(3)
        endif
      enddo 
!
!  Output k points to file
!
      write(iout,'(i8,'' Wedge k point  : '',8x,3(1x,f9.6))') nk,(kp(j),j=1,3)
      do k = 1,nki
        if (k.eq.1) then
          write(iout,'(8x,'' Image k points : '',i8,3(1x,f9.6))') k,(kpall(j,k),j=1,3)
        else
          write(iout,'(26x,i8,3(1x,f9.6))') k,(kpall(j,k),j=1,3)
        endif
      enddo
!
!  End of condition on being the correct configuration
!
    endif
  enddo
!
!  Close file
!
  close(iout)
!
  return
  end
