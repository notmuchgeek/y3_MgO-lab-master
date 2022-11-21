  subroutine testspacegroup(ngroup,lvalid)
!
!  Tests whether a space group is valid for a given set of atoms
!
!   8/13 Created from equpos/symgenatom
!   2/18 Trace added
!
!  On entry:
!
!  ngroup = Space group number to test
!
!  On exit:
!
!  lvalid = if true, then this would be a valid space group for the set of atoms
!
  use current
  use symmetry
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)   :: ngroup
  logical,     intent(out)  :: lvalid
!
!  Local variables
!
  integer(i4)               :: mv
  integer(i4)               :: i
  integer(i4)               :: n
  integer(i4)               :: nin
  logical                   :: lfound
  real(dp),            save :: thresh = 0.00001_dp
  real(dp)                  :: x(3)
  real(dp)                  :: xx(3)
  real(dp)                  :: v(3)
#ifdef TRACE
  call trace_in('testspacegroup')
#endif
!
!  Initialise return logical
!
  lvalid = .true.
!
!  Loop over symmetry operators to test whether all of them are valid for the current structure
!
!  NB: Start from second operator since the first operator is the identity and so bound to be valid
!
  mv = 1
  do while (mv.lt.ngo.and.lvalid)
    mv = mv + 1
!
!  Loop over atoms to check that there is an equivalent atom under each operator
!
    i = 0
    do while (i.lt.numat.and.lvalid)
      i = i + 1
!
!  Set fractional coordinates of atom
!
      xx(1) = xfrac(i)
      xx(2) = yfrac(i)
      xx(3) = zfrac(i)
!
!  Roto-translation
!
      x(1) = 0.0_dp
      x(2) = 0.0_dp
      x(3) = 0.0_dp
      v(1) = vit(1,mv)
      v(2) = vit(2,mv)
      v(3) = vit(3,mv)
      call GULP_mxmb(rop(1,1,mv),1_i4,3_i4,xx,1_i4,1_i4,v,1_i4,1_i4,3_i4,3_i4,1_i4)
!
!  Transform coordinates to primitive cell
!
      if (ngocfg(ncf).gt.1) then
        x(1) = v(1)
        x(2) = v(2)
        x(3) = v(3)
      else
        call GULP_mxmb(w1(ncbl,1,1),7_i4,21_i4,v,1_i4,1_i4,x,1_i4,1_i4,3_i4,3_i4,1_i4)
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
!  Compare atom with other atoms to see if there is a match
!
      lfound = .false.
      n = 0
      do while (n.lt.numat.and..not.lfound)
        n = n + 1
        if (nat(n).eq.nat(i)) then
          v(1) = xfrac(n) - x(1)
          v(2) = yfrac(n) - x(2)
          v(3) = zfrac(n) - x(3)
          v(1) = v(1) - nint(v(1))
          v(2) = v(2) - nint(v(2))
          v(3) = v(3) - nint(v(3))
          lfound = ((abs(v(1))+abs(v(2))+abs(v(3))).lt.thresh) 
        endif
!
!  End of loop over atoms looking for a match
!
      enddo 
!
!  If no match was found then operator is not valid
!
      if (.not.lfound) lvalid = .false.
!
!  End of loop over atoms in unit cell
!
    enddo
!
!  End of loop over symmetry operators
!
  enddo 
#ifdef TRACE
  call trace_out('testspacegroup')
#endif
!
  return
  end 
