  subroutine specialqtn(loptm,lcheck)
!
!  Check which rigid molecule quaternion variables are allowed to move due to special positions
!  Displace each variable in turn and see if the site multiplicity is changed. 
!  This may not be elegant but it works for all space groups!
!
!  NB: At the moment this doesn't allow for combinations of quaternions with constraints
!
!  loptm  = array of logical flags according to whether a
!           parameter can be varied or not. Passed from setcfg.
!  lcheck = logical variable. If true then setting of ltmp is
!           tested for correctness rather than set.
!
!   2/20 Created from specialcom
!   3/20 Molecule number added to setquaternion arguments
!   3/20 Magnitude of the delta for testing increased otherwise changes can be too small
!   7/20 Change to use of quaternions in light of change to symmetry algorithm
!   7/20 Centre of mass is now set from molcom in centring is present due to symmetry
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
!  Julian Gale, CIC, Curtin University, July 2020
!
  use configurations
  use control
  use current
  use iochannels
  use optimisation
  use molecule
  use parallel
  use symmetry
  use times
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  logical, intent(inout) :: loptm(6,nmolasym)
  logical, intent(in)    :: lcheck
!
!  Local variables
!
  character(len=6)       :: molno
  integer(i4)            :: icosxi
  integer(i4)            :: icosyi
  integer(i4)            :: icoszi
  integer(i4)            :: ii
  integer(i4)            :: indm
  integer(i4)            :: ix
  integer(i4)            :: iy
  integer(i4)            :: iz
  integer(i4)            :: j
  integer(i4)            :: m
  integer(i4)            :: mv
  integer(i4)            :: na
  integer(i4)            :: neqvna
  integer(i4)            :: nf
  integer(i4)            :: nfinish
  integer(i4)            :: nfirst
  integer(i4)            :: ngen
  integer(i4)            :: nin
  integer(i4)            :: nm
  integer(i4)            :: nmaf
  integer(i4)            :: nmf
  logical                :: lnocon
  logical                :: ltmploc(3)
  real(dp)               :: drQ(3,3,3)      ! First derivatives of the quaternion matrix - not used in this routine
  real(dp)               :: drQ2(3,3,3,3)   ! Second derivatives of the quaternion matrix - not used in this routine
  real(dp)               :: rQ(3,3)         ! Quaternion matrix
  real(dp)               :: g_cpu_time
  real(dp)               :: mQtmp(3)        ! Local copy of quaternions to be modified when testing
  real(dp)               :: ra(3)           ! Vector from COM to atom after quaternion rotation
  real(dp)               :: thresh
  real(dp)               :: time1
  real(dp)               :: time2
  real(dp)               :: v(3)
  real(dp)               :: xcom
  real(dp)               :: ycom
  real(dp)               :: zcom
  real(dp)               :: xcrd
  real(dp)               :: ycrd
  real(dp)               :: zcrd
  real(dp)               :: x(3)
  real(dp)               :: xorig
  real(dp)               :: xstor(48)
  real(dp)               :: xx(3)
  real(dp)               :: yorig
  real(dp)               :: ystor(48)
  real(dp)               :: zorig
  real(dp)               :: zstor(48)
#ifdef TRACE
  call trace_in('specialqtn')
#endif
!
  time1 = g_cpu_time()
  thresh = 1.0d-4
!
!  Setup local variables
!
  lnocon = (index(keyword,'noco').ne.0)
!********************************************************
!  Test for each asymmetric unit atom in each molecule  *
!********************************************************
  do nm = 1,nmolasym
    nmaf = nmola2f(nm)
!
!  Convert fractional coordinates for COM to Cartesian
!
    xcom = molcom(1,nmaf)*r1x + molcom(2,nmaf)*r2x + molcom(3,nmaf)*r3x
    ycom = molcom(1,nmaf)*r1y + molcom(2,nmaf)*r2y + molcom(3,nmaf)*r3y
    zcom = molcom(1,nmaf)*r1z + molcom(2,nmaf)*r2z + molcom(3,nmaf)*r3z
!
    ltmploc(1:3) = .true.
!*****************************************
!  First pass - individual quaternions  *
!*****************************************
    iiloop: do ii = 1,3
!
!  Set quaternion for testing
!
      mQtmp(1:3) = molQa(1:3,nm)
      mQtmp(ii) = mQtmp(ii) + 0.01_dp
!
!  Set up quaternion matrix
!
      call setquaternion(nmaf,mQtmp,rQ,drQ,drQ2,.false.,.false.)
!
      mloop: do m = 1,nmolasymno(nm)
        na = nmolasymptr(m,nm)
        nmf = natmol(nrela2f(na))
        nf = natinmol(nrela2f(na))
!
!  Find cell image for atom
!
        indm = nmolind(nrela2f(na))
        call mindtoijk(indm,ix,iy,iz)
!
!  Apply quaternions to molecule orientation
!
        do j = 1,3
          ra(j) = rQ(j,1)*molQxyz(1,nf,nm) + rQ(j,2)*molQxyz(2,nf,nm) + rQ(j,3)*molQxyz(3,nf,nm)
        enddo
!
        xcrd = xcom + ra(1) - dble(ix)*rv(1,1) &
                            - dble(iy)*rv(1,2) &
                            - dble(iz)*rv(1,3)
        ycrd = ycom + ra(2) - dble(ix)*rv(2,1) &
                            - dble(iy)*rv(2,2) &
                            - dble(iz)*rv(2,3)
        zcrd = zcom + ra(3) - dble(ix)*rv(3,1) &
                            - dble(iy)*rv(3,2) &
                            - dble(iz)*rv(3,3)
!
!  Convert molecule atom coordinates back to fractional
!
        call cart2frac(ndim,xcrd,ycrd,zcrd,rv,xorig,yorig,zorig,icosxi,icosyi,icoszi)
!
        nfirst = nrela2f(na)
        neqvna = neqv(na)
!
        xx(1) = xorig
        xx(2) = yorig
        xx(3) = zorig
!
!  First symmetry operator
!
        x(1) = 0.0_dp
        x(2) = 0.0_dp
        x(3) = 0.0_dp
        v(1) = xx(1)
        v(2) = xx(2)
        v(3) = xx(3)
!
!  Transform coordinates to primitive cell
!
        call GULP_mxmb(w1(ncbl,1,1),7_i4,21_i4,v,1_i4,1_i4,x,1_i4,1_i4,3_i4,3_i4,1_i4)
!
!  Place fractional coordinates in range 0-1
!
        x(1) = x(1) + 3.0_dp
        nin = x(1)
        xstor(1) = x(1) - nin
        x(2) = x(2) + 3.0_dp
        nin = x(2)
        ystor(1) = x(2) - nin
        x(3) = x(3) + 3.0_dp
        nin = x(3)
        zstor(1) = x(3) - nin
!
!  Loop over symmetry operators
!
        ngen = 1
        ngoloop: do mv = 2,ngo
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
          call GULP_mxmb(w1(ncbl,1,1),7_i4,21_i4,v,1_i4,1_i4,x,1_i4,1_i4,3_i4,3_i4,1_i4)
!
!  Place fractional coordinates in range 0-1
!
          x(1) = x(1) + 3.0_dp
          nin = x(1)
          x(1) = x(1) - nin
          x(2) = x(2) + 3.0_dp
          nin = x(2)
          x(2) = x(2) - nin
          x(3) = x(3) + 3.0_dp
          nin = x(3)
          x(3) = x(3) - nin
!
!  Compare atom with previously generated equivalent
!
          nfinish = nfirst + ngen - 1
          do j = 1,ngen
            v(1) = xstor(j) - x(1)
            v(2) = ystor(j) - x(2)
            v(3) = zstor(j) - x(3)
            v(1) = v(1) - nint(v(1))
            v(2) = v(2) - nint(v(2))
            v(3) = v(3) - nint(v(3))
            if ((abs(v(1))+abs(v(2))+abs(v(3))).lt.thresh) cycle ngoloop
          enddo
!
!  Atom is not equivalent to any previous atom.
!  If new unique atom has been generated by the same space group
!  operator as for the original coordinates then there has been
!  no change in multiplicity.
!
          ngen = ngen + 1
          if (ngen.gt.neqvna) then
            ltmploc(ii) = .false.
            cycle iiloop
          elseif (nrotop(nfinish+1).ne.mv) then
            ltmploc(ii) = .false.
            cycle iiloop
          endif
          xstor(ngen) = x(1)
          ystor(ngen) = x(2)
          zstor(ngen) = x(3)
        enddo ngoloop
!
!  End of loop over atoms in molecule
!
      enddo mloop
!
!  End of loop over quaternions
!
    enddo iiloop
!
    if (lcheck) then
!
!  Compare flags against input values if in checking mode
!
      if (loptm(4,nm).and..not.ltmploc(1)) then
        call itow(molno,nm,6_i4)
        call outerror('Badly defined flags for '//'molecule '//molno,0_i4)
        call stopnow('specialqtn')
      endif
      if (loptm(5,nm).and..not.ltmploc(2)) then
        call itow(molno,nm,6_i4)
        call outerror('Badly defined flags for '//'molecule '//molno,0_i4)
        call stopnow('specialqtn')
      endif
      if (loptm(6,nm).and..not.ltmploc(3)) then
        call itow(molno,nm,6_i4)
        call outerror('Badly defined flags for '//'molecule '//molno,0_i4)
        call stopnow('specialqtn')
      endif
    else
!
!  Transfer flags to main array - only do this if the flag has initially been 
!  set to true though so that user flags can be validated.
!
      if (loptm(4,nm)) loptm(4,nm) = ltmploc(1)
      if (loptm(5,nm)) loptm(5,nm) = ltmploc(2)
      if (loptm(6,nm)) loptm(6,nm) = ltmploc(3)
    endif
!
!  End of loop over molecules
!
  enddo
  time2 = g_cpu_time()
  tsym = tsym + time2 - time1
#ifdef TRACE
  call trace_out('specialqtn')
#endif
!
  return
  end
