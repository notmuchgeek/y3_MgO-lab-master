  subroutine equpos(lprint,lfirst)
!
!  Generate symmetry equivalent atoms from asymmetric unit
!
!   3/07 Calls to mxmb renamed
!   9/15 Added number of cells for lmodco increased from 10 to 1000
!   2/18 Trace added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!
!  nasym  = Number of "irreducible" atoms
!
!  lprint = controls whether or not to output equivalence error message
!  lfirst = indicates whether this is the first call for a structure - if
!           so some arrays are initialised for use on subsequent calls
!
!  This version only expands the coordinates to the full fractional
!  set, while the rest of GULP handles the conversion to cartesian.
!
  use configurations
  use current
  use times
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  logical, intent(in)  :: lprint
  logical, intent(in)  :: lfirst
!
!  Local variables
!
  integer(i4)          :: na
  integer(i4)          :: naf
  real(dp)             :: g_cpu_time
  real(dp)             :: time1
  real(dp)             :: time2
  real(dp)             :: xx(3)
#ifdef TRACE
  call trace_in('equpos')
#endif
!
  time1 = g_cpu_time()
  if (nasym.gt.maxat) then
    maxat = nasym + 10
    call changemaxat
  endif
  if (lfirst) then
!
!  Calculate equivalent positions
!
    naf = 0
    do na = 1,nasym
      xx(1) = xafrac(na)
      xx(2) = yafrac(na)
      xx(3) = zafrac(na)
      call symgenatom(na,naf,xx,lprint)
    enddo
    numat = naf
  else
!
!  Update existing atoms based on symmetry
!
    call symupdateatom
  endif
!
  time2 = g_cpu_time()
  tsym = tsym + time2 - time1
#ifdef TRACE
  call trace_out('equpos')
#endif
!
  return
  end
!
!  Subroutines called from equpos :
!
  subroutine symgenatom(natom,naf,xx,lprint)
!
!  Generate symmetry equivalent atoms from a single given atom
!
!  natom   = Number of atom to which symmetry operators are to be applied
!  naf     = Total number of atoms in unit cell so far
!  xx(3)   = fractional coordinates of atom
!  iatn    = "irreducible" atoms atomic number
!  neqv    = number of symmetry equivalent atoms for each asymmetric unit one
!  nrelf2a = pointer from full to asymmetric atom
!  nrela2f = pointer from asymmetric to full atom
!  lprint  = controls whether or not to output equivalence error message
!
!  11/06 Modified so that testing of operators is only performed on 
!        call with lfirst = .true. and subsequently saved information is
!        used to generate a consistent set of images
!   2/18 Trace added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!
  use configurations
  use control
  use current
  use element,       only : maxele
  use general,       only : nwarn
  use iochannels
  use parallel
  use symmetry
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4)       :: natom
  integer(i4)       :: naf
  logical           :: lprint
  real(dp)          :: xx(3)
!
!  Local variables
!
  character(len=4)  :: atno1
  character(len=4)  :: atno2
  integer(i4)       :: mv
  integer(i4)       :: n
  integer(i4)       :: naf0
  integer(i4)       :: nin
  integer(i4), save :: ncfold = 0_i4
  real(dp),    save :: thresh = 0.00001_dp
  real(dp)          :: x(3)
  real(dp)          :: xoffset(3)
  real(dp)          :: v(3)
#ifdef TRACE
  call trace_in('symgenatom')
#endif
!
!  Calculate equivalent positions
!
  naf0 = naf
  nrela2f(natom) = naf + 1
!*******************************************************************************
!  Loop over symmetry operators to generate symmetry related atom information  *
!*******************************************************************************
  mvloop: do mv = 1,ngo
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
!  Check array dimensions
!
    if (naf.ge.maxat) then
      maxat = naf + 10
      call changemaxat
    endif
!
    xoffset(1) = x(1)
    xoffset(2) = x(2)
    xoffset(3) = x(3)
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
    xoffset(1) = x(1) - xoffset(1)
    xoffset(2) = x(2) - xoffset(2)
    xoffset(3) = x(3) - xoffset(3)
!
!  Compare atom with those previously generated
!
    nloop: do n = 1,naf
      v(1) = xfrac(n) - x(1)
      v(2) = yfrac(n) - x(2)
      v(3) = zfrac(n) - x(3)
      v(1) = v(1) - nint(v(1))
      v(2) = v(2) - nint(v(2))
      v(3) = v(3) - nint(v(3))
      if ((abs(v(1))+abs(v(2))+abs(v(3))).ge.thresh) cycle nloop
!
!  Changed to allow for partial occupancies at same coordinates
!
      if (abs(iatn(natom)-nat(n)).gt.0.or.((occua(natom)+occuf(n)).le.1.000001d0.and. &
        nrelf2a(n).ne.natom).or.abs(iatn(natom)-nat(n)).eq.maxele) cycle nloop
      if (mv.ne.1) cycle mvloop
      if (lprint.and.ncfold.ne.ncf) then
        nwarn = nwarn + 1
        if (ioproc) then
          call itow(atno1,natom,4_i4)
          call itow(atno2,n,4_i4)
          call outwarning('The irreducible atom '//atno1//' is equivalent to atom '//atno2,0_i4)
        endif
      endif
      return
    enddo nloop
    naf = naf + 1
    nat(naf) = iatn(natom)
    nftype(naf) = natype(natom)
    nrelf2a(naf) = natom
    nrotop(naf) = mv
    xfrac(naf) = x(1)
    yfrac(naf) = x(2)
    zfrac(naf) = x(3)
    xfracimage(naf) = xoffset(1)
    yfracimage(naf) = xoffset(2)
    zfracimage(naf) = xoffset(3)
    c6f(naf) = c6a(natom)
    qf(naf) = qa(natom)
    occuf(naf) = occua(natom)
    radf(naf) = rada(natom)
    oxf(naf) = oxa(natom)
    cnf(naf) = cna(natom)
    neqv(natom) = naf - naf0
  enddo mvloop
#ifdef TRACE
  call trace_out('symgenatom')
#endif
!
  return
  end

  subroutine symupdateatom
!
!  Update symmetry equivalent atoms from those in the asymmetric unit
!
!  11/06 Created from symgenatom so that a consistent image would be
!        generated between calls of equpos for the same structure.
!  11/07 Unused variables removed
!   2/18 Trace added
!
  use configurations
  use control
  use current
  use symmetry
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Local variables
!
  integer(i4)       :: mv
  integer(i4)       :: naf
  integer(i4)       :: natom
  real(dp)          :: x(3)
  real(dp)          :: xx(3)
  real(dp)          :: v(3)
#ifdef TRACE
  call trace_in('symupdateatom')
#endif
!*******************************************************************
!  Update coordinates based on previously saved list of operators  *
!*******************************************************************
!
!  Loop over full list of atoms
!
  do naf = 1,numat
!
!  Update properties from asymmetric unit atom
!
    natom = nrelf2a(naf)
    mv = nrotop(naf)
    c6f(naf) = c6a(natom)
    qf(naf) = qa(natom)
    radf(naf) = rada(natom)
!
!  Set coordinates of symmetry related atom
!
    xx(1) = xafrac(natom)
    xx(2) = yafrac(natom)
    xx(3) = zafrac(natom)
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
!  Shift to correct cell image of atom and place in return array
!
    xfrac(naf) = x(1) + xfracimage(naf)
    yfrac(naf) = x(2) + yfracimage(naf)
    zfrac(naf) = x(3) + zfracimage(naf)
!
!  Make sure all fractional coords are between 0 and 1 if required
!
    if (lmodco) then
      xfrac(naf) = mod(xfrac(naf)+1000.0_dp,1.0_dp)
      yfrac(naf) = mod(yfrac(naf)+1000.0_dp,1.0_dp)
      zfrac(naf) = mod(zfrac(naf)+1000.0_dp,1.0_dp)
    endif
  enddo
#ifdef TRACE
  call trace_out('symupdateatom')
#endif
!
  return
  end
!
  subroutine symupdatemol
!
!  Update symmetry equivalent molecules from those in the asymmetric unit
!
!  10/19 Created from symupdateatom
!   3/20 Quaternion update added
!   3/20 Symmetry adaption of quaternion update added for lmolstdframe case
!   5/20 Symmetry update of quaternions no longer dependent on frame
!   7/20 Symmetry adaption of molecular quaternions changed
!        to allow for non-orthogonal cell case
!   7/20 w1 included in matrix for symmetry adaption of quaternions
!
  use configurations
  use control
  use current
  use molecule
  use symmetry
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Local variables
!
  integer(i4)          :: i
  integer(i4)          :: ifail
  integer(i4)          :: j
  integer(i4)          :: k
  integer(i4)          :: mv
  integer(i4)          :: nm
  integer(i4)          :: nma
  real(dp)             :: rvinv(3,3)
  real(dp)             :: Rcs(3,3)
  real(dp)             :: Rmat(3,3)
  real(dp)             :: v(3)
  real(dp)             :: wrk(6)
  real(dp)             :: x(3)
  real(dp)             :: xx(3)
#ifdef TRACE
  call trace_in('symupdatemol')
#endif
!*******************************************************************
!  Update coordinates based on previously saved list of operators  *
!*******************************************************************
!
!  Set up cell matrix: For centred cells this is a product of the 
!  primitive cell and centring transformation matrix w1
!
  if (ngocfg(ncf).gt.1) then
    rvinv(1:3,1:3) = rv(1:3,1:3)
  else
    do i = 1,3
      do j = 1,3
        rvinv(j,i) = 0.0_dp
        do k = 1,3
          rvinv(j,i) = rvinv(j,i) + w1(ncbl,j,k)*rv(k,i)
        enddo
      enddo
    enddo
  endif
  Rmat(1:3,1:3) = rvinv(1:3,1:3)
!
!  Invert cell matrix 
!
  call matrix_inversion(rvinv,3_i4,3_i4,wrk,ifail)
!
!  Loop over full list of molecules
!
  do nm = 1,nmol
!
!  Find first atom in molecule to get symmetry operators
!
    mv = nrotop(nmollist(nmolptr(nm)+1))
!
!  Set coordinates of symmetry related molecule
!
    xx(1) = molcoma(1,nmolf2a(nm))
    xx(2) = molcoma(2,nmolf2a(nm))
    xx(3) = molcoma(3,nmolf2a(nm))
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
    molcom(1,nm) = x(1)
    molcom(2,nm) = x(2)
    molcom(3,nm) = x(3)
!
!  Make sure all fractional coords are between 0 and 1 if required
!
    if (lmodco) then
      molcom(1,nm) = mod(molcom(1,nm)+1000.0_dp,1.0_dp)
      molcom(2,nm) = mod(molcom(2,nm)+1000.0_dp,1.0_dp)
      molcom(3,nm) = mod(molcom(3,nm)+1000.0_dp,1.0_dp)
    endif
!
!  Update quaternions - set based on symmetry-related molecule in the asymmetry unit
!
    nma = nmolf2a(nm)
    molQ(1:3,nm) = molQa(1:3,nma)
!
!  Form Rc*Rs*Rc-1
!
    do i = 1,3
      do j = 1,3
        Rcs(j,i) = 0.0_dp
        do k = 1,3
          Rcs(j,i) = Rcs(j,i) + Rmat(j,k)*rop(k,i,mv)
        enddo
      enddo
    enddo
!
    do i = 1,3
      do j = 1,3
        molQsym(j,i,nm) = 0.0_dp
        do k = 1,3
          molQsym(j,i,nm) = molQsym(j,i,nm) + Rcs(j,k)*rvinv(k,i)
        enddo
      enddo
    enddo
  enddo
#ifdef TRACE
  call trace_out('symupdatemol')
#endif
!
  return
  end
