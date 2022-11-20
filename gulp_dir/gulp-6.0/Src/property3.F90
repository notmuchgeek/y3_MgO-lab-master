  subroutine property3(lprint)
!
!  Calculate the properties that depend on the second derivatives. 3-D
!
!  Calculate high freq dielectric constant first so that
!  on return the inverse second derivative matrix is
!  still available in dervi for use by defect calculations.
!  If optical relaxation then the order must be reversed!
!
!  Modified for breathing shell model
!
!  nbs   = no. of full cell atoms with breathing shells
!  nbss  = no. of full cell shells with breathing shells
!  nbsptr= pointer to atoms with breathing shells
!
!   8/98 Breathing shell pointer generation moved to subroutine
!   8/98 Calculation of refractive index added
!   5/99 Calculation of Young's moduli and Poisson's ratios added
!  10/01 Pressure correction for C44,55,66 corrected by a half
!   5/02 All three conventions for bulk modulus now present
!   5/02 Shear modulus added
!   8/02 S and P wave velocities added for polycrystalline average
!   8/02 Made implicit none
!   9/02 Output of elastic compliances added
!   9/02 conve now defined using parameters
!  10/02 Problems with negative shear moduli in sqrt handled
!   2/03 All mechanical properties now internally in GPa
!   2/03 Labels for piezoelectric strain/stress switched to be correct
!   2/03 Matrix inversion accelerated
!   3/03 Partial occupancies with multiple ions on the same site handled
!   6/03 Calculation of compressibility corrected for switch to GPa
!   4/04 Dimensions of tmp now set to a minimum of 12 to avoid error
!  10/04 Eispack calls replaced by lapack
!   5/06 Mass now taken from species masses
!   5/07 Partial occupancy data moved to module
!  12/07 Unused variables removed
!   1/09 Integer datatypes all explicitly declared
!   2/09 Modified to accommodate new version of FoX and gulp_cml
!   4/09 Calling of changemaxd2 modified to handle case where no internal
!        second derivatives are required.
!   6/09 Extra space added ahead of "Mechanical properties"
!   1/10 Young's moduli and Poisson's ratios store in module arrays
!   1/10 Incorrect factor of 10 in velocities corrected
!  10/11 Stress tensor output moved to property from optout
!   5/12 Stresses moved to separate routine
!   8/13 More of the properties moved from being local variables to being
!        in the module, properties
!   9/13 Calculation of Raman susceptibilities for atoms added
!  10/13 Renamed from property to property3 
!   3/14 Calls to matinv renamed to matrix_inversion for benefit of ChemShell
!   9/16 Modified to allow for parallel second derivatives
!   9/16 Matrix inversion moved into subroutine
!   9/16 Shell pointers modified
!   2/17 Blocksize added to call to matrix_inversion_library
!   5/17 nshell and nshellonnode added to subroutine call
!   6/17 Parallelisation for shell only case added
!   7/17 Use of qshell corrected to qshellnode for diconh
!   8/17 Parallel handling of qD added
!   1/18 Trace added
!   5/18 Output of elastic constant eigenvalues added
!  11/18 Output of elastic constant eigenvectors added
!   7/19 Keyword added to control output of elastic eigen properties
!   8/19 Separate flag introduced for piezoelectric stress so that it
!        can be disable if compliances fail
!   4/20 Rigid molecule modifications added
!   5/20 Change to use of rotations for rigid molecules 
!   5/20 Removal of rotations from internal contribution when inertia is too small
!        or there are negligible angular forces
!   6/20 Parallel modifications made for rigid molecules
!   6/20 nmolcore changes added
!   7/20 tmp2 saved from elastic constants and re-used for piezoelectrics
!   9/20 lpolarisation added
!   9/20 tmp2 no longer allocated if linternal is false
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
!  Julian Gale, CIC, Curtin University, September 2020
!
  use g_constants
  use control
  use current
  use derivatives
  use element
  use gulp_cml,       only : lcml
  use gulp_cml_props, only : gulp_cml_output_derivs, gulp_cml_output_dielectric
  use general
  use iochannels
  use molecule
  use parallel
  use partial
  use properties
  use shells
  use species,        only : massspec
  use symmetry
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
#ifdef MPI
  include 'mpif.h'
#endif
!
!  Local variables
!
  integer(i4), dimension(:),   allocatable       :: nmOK
  integer(i4), dimension(:,:), allocatable       :: nmOKptr
  real(dp),    dimension(:,:), allocatable       :: d3
#ifdef MPI
  real(dp),    dimension(:,:), allocatable       :: d3nomol
  real(dp),    dimension(:,:), allocatable       :: d3nomols
#endif
  real(dp),    dimension(:),   allocatable       :: qfo
  real(dp),    dimension(:),   allocatable       :: qshell
  real(dp),    dimension(:),   allocatable       :: qshellnode
  real(dp),    dimension(:),   allocatable       :: tmp
  real(dp),    dimension(:,:), allocatable       :: tmp2
  real(dp),    dimension(:,:), allocatable       :: tmp2s
  integer(i4)                                    :: i
  integer(i4)                                    :: ifail
  integer(i4)                                    :: ii
  integer(i4)                                    :: is
  integer(i4)                                    :: ix
  integer(i4)                                    :: iy
  integer(i4)                                    :: iz
  integer(i4)                                    :: ixs
  integer(i4)                                    :: iys
  integer(i4)                                    :: izs
  integer(i4)                                    :: indif
  integer(i4)                                    :: indil
  integer(i4)                                    :: indjf
  integer(i4)                                    :: indjl
  integer(i4)                                    :: indkf
  integer(i4)                                    :: indkl
  integer(i4)                                    :: indk
  integer(i4)                                    :: indl
  integer(i4)                                    :: iptr
  integer(i4)                                    :: j
  integer(i4)                                    :: jj
#ifdef MPI
  integer(i4)                                    :: jloc
#endif
  integer(i4)                                    :: jptr
  integer(i4)                                    :: jx
  integer(i4)                                    :: jy
  integer(i4)                                    :: jz
  integer(i4)                                    :: jxs
  integer(i4)                                    :: jys
  integer(i4)                                    :: jzs
  integer(i4)                                    :: k
  integer(i4)                                    :: ki
  integer(i4)                                    :: kk
  integer(i4)                                    :: kx
  integer(i4)                                    :: l
  integer(i4)                                    :: mind
  integer(i4)                                    :: mindi
  integer(i4)                                    :: mindj
  integer(i4)                                    :: mindk
  integer(i4)                                    :: mint
  integer(i4)                                    :: mintu
  integer(i4)                                    :: msvar
  integer(i4)                                    :: ninv
  integer(i4)                                    :: ninvs
  integer(i4)                                    :: n3f
  integer(i4)                                    :: n3fu
  integer(i4)                                    :: natomsonnodem1
  integer(i4)                                    :: nbsc
  integer(i4)                                    :: nelconeig
  integer(i4)                                    :: ni
  integer(i4)                                    :: nj
  integer(i4)                                    :: nmolfree
  integer(i4)                                    :: nqfo
  integer(i4)                                    :: numatnomolfree
  integer(i4)                                    :: numatnomolfreeonnode
  integer(i4)                                    :: status
  logical                                        :: lcompliances
  logical                                        :: ldlch
  logical                                        :: ldlcs
  logical                                        :: lelastic
  logical                                        :: lelconeigok
  logical                                        :: linternal
  logical                                        :: lmoduli
  logical                                        :: lpiezo
  logical                                        :: lpiezs
  logical                                        :: lprint
  real(dp)                                       :: cfactor
  real(dp)                                       :: conve
  real(dp)                                       :: g_cpu_time
  real(dp)                                       :: dlcfactor
  real(dp)                                       :: ditmp(3,3)
  real(dp)                                       :: elconeig(6)
  real(dp)                                       :: elconvec(6,6)
  real(dp)                                       :: eltmp(6,6)
  real(dp)                                       :: eltmp2(6,6)
  real(dp)                                       :: epv
  real(dp)                                       :: piefactor
  real(dp)                                       :: polar_xyz(3)
  real(dp)                                       :: qak
  real(dp),    dimension(:,:), allocatable       :: qD
  real(dp),    dimension(:,:), allocatable       :: sum2
  real(dp)                                       :: r3(3)
  real(dp)                                       :: r33(3,3)
  real(dp)                                       :: r43
  real(dp)                                       :: rvol
  real(dp)                                       :: rmolfct
  real(dp)                                       :: sum
  real(dp)                                       :: t1p
  real(dp)                                       :: t2p
  real(dp)                                       :: totmass
  real(dp)                                       :: vol
  real(dp)                                       :: volume
  real(dp)                                       :: w1l(9)
  real(dp)                                       :: zero
#ifdef MPI
  integer(i4)                                    :: icount
  integer(i4)                                    :: iloc
  integer(i4)                                    :: ndi
  integer(i4)                                    :: ndof
  integer(i4)                                    :: node
  real(dp),      dimension(:,:), allocatable     :: dtmp
!
!  Local variables in Scalapack/Blacs/MPI integer precision
!
  integer                                       :: nsize
  integer                                       :: ntag
  integer                                       :: rnode
  integer                                       :: snode
  integer                                       :: Request
  integer                                       :: MPIerror
  integer(i4),  dimension(:), allocatable       :: pdig2l
  integer(i4),  dimension(:), allocatable       :: pdinode
  integer,      dimension(:), allocatable       :: StatMPI       ! Array for status from MPI
  logical                                       :: ld2loc
  logical                                       :: ldiloc
#endif
!
  data cfactor/6.241460893d-3/
#ifdef TRACE
  call trace_in('property3')
#endif
!
  t1p = g_cpu_time()
!
  lcompliances = .true.
  lelastic = .true.
  ldlch = .true.
  ldlcs = .true.
  lmoduli = .true.
  lpiezo = .true.
  lpiezs = .true.
  conve = evtoj*1.0d20
  zero = 1.0d-12
  vol = volume(rv)
!
  if (lrigid) then
    n3f = 3*numatnomol + 6*nmol
    mint = n3f - 3
!
!  Decide whether an atom or molecule will be fixed for internal response matrix
!
    if (numatnomol.gt.0) then
      numatnomolfree = numatnomol - 1
      nmolfree = nmol
    else
      numatnomolfree = 0
      nmolfree = nmol - 1
    endif
!
!  Set up rotational derivatives
!
    call rigidmoleculeprop
  else
    n3f = 3*numat
    mint = n3f - 3
  endif
!
!  Check second derivative memory
!
  if (nprocs.gt.1) then
    n3fu = 3_i4*natomsonnode
    mintu = n3fu
    if (procid.eq.atom2node(numat)) then
      mintu = mintu - 3
      natomsonnodem1 = natomsonnode - 1
    else
      natomsonnodem1 = natomsonnode
    endif
    if (lrigid) then
      if (numatnomol.gt.0) then
        if (procid.eq.atom2node(numatnomolptr(numatnomol))) then
          numatnomolfreeonnode = numatnomolonnode - 1
        else
          numatnomolfreeonnode = numatnomolonnode
        endif
      else
        numatnomolfreeonnode = numatnomolonnode
      endif
    endif
  else
    n3fu = n3f
    mintu = mint
    natomsonnodem1 = natomsonnode - 1
    if (lrigid) then
      numatnomolfreeonnode = numatnomolfree
    endif
  endif
  if (n3f+6.gt.maxd2.or.n3fu+6.gt.maxd2u) then
    maxd2  = max(n3f  + 6,maxd2)
    maxd2u = max(n3fu + 6,maxd2u)
    call changemaxd2
  endif
#ifdef MPI
  if (lrigid.and.nprocs.gt.1) then
!
!  Compute distribution of atoms and molecules for rigid molecule case
!
!  => Block cyclic distribution with numatnomolfree followed by nmolfree in block multiples of 3
!
    if (lnorotate) then
      ndof = numatnomolfree + nmolfree
    else
      ndof = numatnomolfree + nmolfree + nmol
    endif
    allocate(pdinode(ndof),stat=status)
    if (status/=0) call outofmemory('property3','pdinode')
    allocate(pdig2l(ndof),stat=status)
    if (status/=0) call outofmemory('property3','pdig2l')
!
    node = 0
    ndi = 0
    icount = 0
    do i = 1,ndof
      icount = icount + 1
      pdinode(i) = node
      if (node.eq.procid) then
        ndi = ndi + 1
        pdig2l(i) = ndi
      else
        pdig2l(i) = 0
      endif
      if (icount.eq.nblocksize) then
        icount = 0
        node = node + 1
        if (node.eq.nprocs) node = 0
      endif
    enddo
  endif
#endif
!
!  Allocate local memory
!
  allocate(qfo(numat),stat=status)
  if (status/=0) call outofmemory('property3','qfo')
  allocate(qshell(numat),stat=status)
  if (status/=0) call outofmemory('property3','qshell')
  if (nprocs.gt.1) then
    allocate(qshellnode(natomsonnode),stat=status)
    if (status/=0) call outofmemory('property3','qshellnode')
  endif
  allocate(tmp(max(8*numat,12)),stat=status)
  if (status/=0) call outofmemory('property3','tmp')
  if (lpocc) then
    allocate(d3(maxd2,6),stat=status)
    if (status/=0) call outofmemory('property3','d3')
  endif
  if (lrigid) then
    allocate(nmOK(nmol),stat=status)
    if (status/=0) call outofmemory('property3','nmOK')
    allocate(nmOKptr(3,nmol),stat=status)
    if (status/=0) call outofmemory('property3','nmOKptr')
  endif
#ifdef MPI
  if (lrigid.and.nprocs.gt.1) then
    mind = 3*numatnomolfree
    allocate(d3nomol(mind,6),stat=status)
    if (status/=0) call outofmemory('property3','d3nomol')
    allocate(d3nomols(mind,6),stat=status)
    if (status/=0) call outofmemory('property3','d3nomols')
  endif
#endif
!
!  Set logical as to whether internal contribution to elastic constants are needed or wanted
!
  linternal = (index(keyword,'noin').eq.0)
  if (mint+nbs.eq.0) linternal = .false.
!
  if (ldlcs.or.ldlch.or.lpiezo.or.lpiezs.or.lpolarisation) then
    dlcfactor = 4.0_dp*pi/vol
    piefactor = dlcfactor*conve
    dlcfactor = angstoev*dlcfactor
  endif
!
!  Calculate density
!
  totmass = 0.0_dp
  do i = 1,nasym
    ni = nspecptr(i)
    totmass = totmass + massspec(ni)*occua(i)*dble(neqv(i))
  enddo
  if (vol.gt.1d-12) then
    rmolfct = avogadro*1.0d-23
    density = (10.0_dp*totmass)/(vol*rmolfct)
  endif
!
!  Set up charge times occupancy array for later use
!
  if (lrigid) then
    nqfo = numatnomolfree + nmolfree
    do i = 1,numatnomolfree
      ii = numatnomolptr(i)
      qfo(i) = qf(ii)*occuf(ii)
    enddo
    do i = 1,nmolfree
      qfo(numatnomolfree+i) = 0.0_dp
      do j = 1,nmolcore(i)
        k = nmollist(nmolptr(i)+j)
        qfo(numatnomolfree+i) = qfo(numatnomolfree+i) + qf(k)*occuf(k)
      enddo
    enddo
  else
    nqfo = ncsfoc - 1
    do i = 1,numat
      qfo(i) = qf(i)*occuf(i)
    enddo
  endif
!**************************************
!  High frequency dielectric constant *
!**************************************
  if (.not.lshello) then
    if (ldlch.and.nshell.gt.0) then
      if (lraman) then
        allocate(qD(3*nsfoc,3),stat=status)
        if (status/=0) call outofmemory('property3','qD')
        if (nprocs.gt.1) then
          allocate(sum2(3*nsfoc,3),stat=status)
          if (status/=0) call outofmemory('property3','sum2')
        endif
      endif
!
!  Collect shell second derivative terms
!
      msvar = 3*nshell
      if (nprocs.gt.1) then
        do i = 1,nshell
          ni = nshptr(i)
          qshell(i) = qf(ni)*occuf(ni)
        enddo
        do i = 1,nshellonnode
          is = nshonnodeptr(i)
          ni = node2atom(is)
          qshellnode(i) = qf(ni)*occuf(ni)
          ix = 3*(is-1) + 1
          iy = ix + 1
          iz = ix + 2
          ixs = 3*(i-1) + 1
          iys = ixs + 1
          izs = ixs + 2
          do j = 1,nshell
            nj = nshptr(j)
            jx = 3*(nj-1) + 1
            jy = jx + 1
            jz = jx + 2
            jxs = 3*(j-1) + 1
            jys = jxs + 1
            jzs = jxs + 2
            dervi(jxs,ixs) = derv2(jx,ix)
            dervi(jys,ixs) = derv2(jy,ix)
            dervi(jzs,ixs) = derv2(jz,ix)
            dervi(jxs,iys) = derv2(jx,iy)
            dervi(jys,iys) = derv2(jy,iy)
            dervi(jzs,iys) = derv2(jz,iy)
            dervi(jxs,izs) = derv2(jx,iz)
            dervi(jys,izs) = derv2(jy,iz)
            dervi(jzs,izs) = derv2(jz,iz)
          enddo
        enddo
      else
        do i = 1,nshell
          ni = nshptr(i)
          qshell(i) = qf(ni)*occuf(ni)
          ix = 3*(ni-1) + 1
          iy = ix + 1
          iz = ix + 2
          ixs = 3*(i-1) + 1
          iys = ixs + 1
          izs = ixs + 2
          do j = 1,nshell
            nj = nshptr(j)
            jx = 3*(nj-1) + 1
            jy = jx + 1
            jz = jx + 2
            jxs = 3*(j-1) + 1
            jys = jxs + 1
            jzs = jxs + 2
            dervi(jxs,ixs) = derv2(jx,ix)
            dervi(jys,ixs) = derv2(jy,ix)
            dervi(jzs,ixs) = derv2(jz,ix)
            dervi(jxs,iys) = derv2(jx,iy)
            dervi(jys,iys) = derv2(jy,iy)
            dervi(jzs,iys) = derv2(jz,iy)
            dervi(jxs,izs) = derv2(jx,iz)
            dervi(jys,izs) = derv2(jy,iz)
            dervi(jzs,izs) = derv2(jz,iz)
          enddo
        enddo
      endif
      if (nbss.gt.0) then
!
!  Collect radial second derivatives
!
        if (lpocc) then
          do i = 1,nshell
            ni = nshptr(i)
            do j = 1,nshell
              nj = nshptr(j)
              dervi(msvar+j,msvar+i) = derv2(n3f+nj,n3f+ni)
            enddo
            do j = 1,nshell
              nj = nshptr(j)
              jx = 3*(nj-1) + 1
              jy = jx + 1
              jz = jx + 2
              jxs = 3*(j-1) + 1
              jys = jxs + 1
              jzs = jxs + 2
              dervi(jxs,msvar+i) = derv2(jx,n3f+ni)
              dervi(jys,msvar+i) = derv2(jy,n3f+ni)
              dervi(jzs,msvar+i) = derv2(jz,n3f+ni)
              dervi(msvar+i,jxs) = derv2(n3f+ni,jx)
              dervi(msvar+i,jys) = derv2(n3f+ni,jy)
              dervi(msvar+i,jzs) = derv2(n3f+ni,jz)
            enddo
          enddo
        else
          nbsc = nbs - nbss
          do i = 1,nbss
            iptr = n3f + nbsptr(nbsc+i)
            do j = 1,nbss
              jptr = n3f + nbsptr(nbsc+j)
              dervi(msvar+j,msvar+i) = derv2(jptr,iptr)
            enddo
            do j = 1,nshell
              nj = nshptr(j)
              jx = 3*(nj-1) + 1
              jy = jx + 1
              jz = jx + 2
              jxs = 3*(j-1) + 1
              jys = jxs + 1
              jzs = jxs + 2
              dervi(jxs,msvar+i) = derv2(jx,iptr)
              dervi(jys,msvar+i) = derv2(jy,iptr)
              dervi(jzs,msvar+i) = derv2(jz,iptr)
            enddo
          enddo
        endif
      endif
!
!  Compress derivatives for partial occupancies
!
      if (lpocc) then
        call compressd1(qshell,0_i4,nsfoc,nshell,iocshptr)
        call compressd2(dervi,maxd2,0_i4,nsfoc,nbsfoc,nshell,iocshptr,ibocshptr)
        ninvs = 3*nsfoc + nbsfoc
      else
        ninvs = msvar + nbss
      endif
!
!  Invert second derivative matrix
!
      ifail = 1
!
!  Matrix inversion
!
      call matrix_inversion_shells(ninvs,1_i4,maxd2,dervi,nshell,nshellonnode,ifail)
!
      if (ifail.ne.0) then
        nwarn = nwarn + 1
        call outwarning('Properties cannot be calculated - matrix is singular',0_i4)
        goto 999
      endif
!
!  Multiply inverse second derivatives by charge vectors:
!
!  First multiply by one charge vector and store in dervi for use in Raman calculation, if needed
!
      if (lpocc) then
        do k = 1,nsfoc
          qak = qshell(k)
          indk = 3*(k-1)
          do indl = 1,3*nsfoc
            dervi(indl,indk+1) = dervi(indl,indk+1)*qak
            dervi(indl,indk+2) = dervi(indl,indk+2)*qak
            dervi(indl,indk+3) = dervi(indl,indk+3)*qak
          enddo
        enddo
      else
        if (nprocs.gt.1) then
          do k = 1,nshellonnode
            qak = qshellnode(k)
            indk = 3*(k-1)
            do indl = 1,3*nshell
              dervi(indl,indk+1) = dervi(indl,indk+1)*qak
              dervi(indl,indk+2) = dervi(indl,indk+2)*qak
              dervi(indl,indk+3) = dervi(indl,indk+3)*qak
            enddo
          enddo
        else
          do k = 1,nshell
            qak = qshell(k)
            indk = 3*(k-1)
            do indl = 1,3*nshell
              dervi(indl,indk+1) = dervi(indl,indk+1)*qak
              dervi(indl,indk+2) = dervi(indl,indk+2)*qak
              dervi(indl,indk+3) = dervi(indl,indk+3)*qak
            enddo
          enddo
        endif
      endif
!
!  Second multiply by second charge vector and compute dielectric tensor
!
      if (lpocc) then
        do i = 1,3
          do j = 1,3
            sum = 0.0_dp
            do k = 1,nsfoc
              indk = 3*(k-1) + i
              do l = 1,nsfoc
                indl = 3*(l-1) + j
                sum = sum + dervi(indl,indk)*qshell(l)
              enddo
            enddo
            diconh(j,i) = dlcfactor*sum
          enddo
          diconh(i,i) = diconh(i,i) + 1.0_dp
        enddo
      else
        if (nprocs.gt.1) then
          do i = 1,3
            do j = 1,3
              sum = 0.0_dp
              do k = 1,nshellonnode
                indk = 3*(k-1) + i
                do l = 1,nshell
                  indl = 3*(l-1) + j
                  sum = sum + dervi(indl,indk)*qshell(l)
                enddo
              enddo
              diconh(j,i) = dlcfactor*sum
            enddo
          enddo
          call sumall(diconh,ditmp,9_i4,"property3","ditmp")
          diconh(1:3,1:3) = ditmp(1:3,1:3)
          do i = 1,3
            diconh(i,i) = diconh(i,i) + 1.0_dp
          enddo
        else
          do i = 1,3
            do j = 1,3
              sum = 0.0_dp
              do k = 1,nshell
                indk = 3*(k-1) + i
                do l = 1,nshell
                  indl = 3*(l-1) + j
                  sum = sum + dervi(indl,indk)*qshell(l)
                enddo
              enddo
              diconh(j,i) = dlcfactor*sum
            enddo
            diconh(i,i) = diconh(i,i) + 1.0_dp
          enddo
        endif
      endif
!
!  Calculate refractive indices
!
      do i = 1,3
        do j = 1,3
          r33(j,i) = diconh(j,i)
        enddo
      enddo
      call dsyev('N','U',3_i4,r33,3_i4,r3,w1l,9_i4,ifail)
      hfrefind(1) = sqrt(dabs(r3(1)))
      hfrefind(2) = sqrt(dabs(r3(2)))
      hfrefind(3) = sqrt(dabs(r3(3)))
      hfrefind(1) = sign(hfrefind(1),r3(1))
      hfrefind(2) = sign(hfrefind(2),r3(2))
      hfrefind(3) = sign(hfrefind(3),r3(3))
!
!  Option to compute atomic Raman susceptibility tensors
!
      if (lraman) then
!
!  Compact dervi down into qD
!
        if (lpocc) then
          qD(1:3*nsfoc,1:3) = 0.0_dp
          do k = 1,nsfoc
            indk = 3*(k-1)
            do i = 1,3
              do indl = 1,3*nsfoc
                qD(indl,i) = qD(indl,i) + dervi(indl,indk+i)
              enddo
            enddo
          enddo
        else
          if (nprocs.gt.1) then
            qD(1:3*nshell,1:3) = 0.0_dp
            do k = 1,nshellonnode
              indk = 3*(k-1)
              do i = 1,3
                do indl = 1,3*nshell
                  qD(indl,i) = qD(indl,i) + dervi(indl,indk+i)
                enddo
              enddo
            enddo
            call sumall(qD,sum2,9_i4*nshell,"property3","qD")
            qD(1:3*nshell,1:3) = sum2(1:3*nshell,1:3)
          else
            qD(1:3*nshell,1:3) = 0.0_dp
            do k = 1,nshell
              indk = 3*(k-1)
              do i = 1,3
                do indl = 1,3*nshell
                  qD(indl,i) = qD(indl,i) + dervi(indl,indk+i)
                enddo
              enddo
            enddo
          endif
        endif
!
!  Call routine to compute the third derivatives and generate the susceptibility tensors
!
        call raman3(qD,3*nsfoc)
!
        if (nprocs.gt.1) then
          deallocate(sum2,stat=status)
          if (status/=0) call deallocate_error('property3','sum2')
        endif
        deallocate(qD,stat=status)
        if (status/=0) call deallocate_error('property3','qD')
      endif
    else
      do i = 1,3
        do j = 1,3
          diconh(j,i) = 0.0_dp
        enddo
        diconh(i,i) = 1.0_dp
        hfrefind(i) = 1.0_dp
      enddo
    endif
  endif
!
!  Set up constants and invert second derivative matrix
!
  if (lelastic.or.ldlcs.or.lpiezo.or.lpiezs) then
    if (linternal) then
      if (lpocc) then
!
!  Copy derv2 into dummy array as matinv overwrites it
!
        do i = 1,n3fu
          do j = 1,n3f
            dervi(j,i) = derv2(j,i)
          enddo
        enddo
        do i = 1,6
          do j = 1,n3fu
            d3(j,i) = derv3(j,i)
          enddo
        enddo
!
!  Radial component
!
        if (nbs.gt.0) then
          do i = 1,numat
            do j = 1,numat
              dervi(n3f+j,n3f+i) = derv2(n3f+j,n3f+i)
            enddo
            do j = 1,mint+3
              dervi(j,n3f+i) = derv2(j,n3f+i)
              dervi(n3f+i,j) = derv2(n3f+i,j)
            enddo
          enddo
          do i = 1,6
            do j = 1,numat
              d3(n3f+j,i) = derv3(n3f+j,i)
            enddo
          enddo
        endif
!
!  Compress partial occupancy second derivatives to full sites
!
        call compressd1(qfo,ncfoc,nsfoc,numat,iocptr)
        call compressd2(dervi,maxd2,ncfoc,nsfoc,nbfoc,numat,iocptr,ibocptr)
        call compressd3(d3,maxd2,ncfoc,nsfoc,nbfoc,numat,iocptr,ibocptr)
!
        ninv = 3*(ncsfoc - 1) + nbfoc
!
!  Now remove last atom from list to avoid singular second derivatives
!
        if (nbfoc.gt.0) then
          indk = 3*(ncsfoc - 1)
          do i = 1,3*ncsfoc+nbfoc
            do j = 1,nbfoc
              dervi(indk+j,i) = dervi(indk+3+j,i)
            enddo
          enddo
          do i = 1,6
            do j = 1,nbfoc
              d3(indk+j,i) = d3(indk+3+j,i)
            enddo
          enddo
          do i = 1,nbfoc
            do j = 1,ninv
              dervi(j,indk+i) = dervi(j,indk+3+i)
            enddo
          enddo
          do i = 1,nbfoc
            do j = 1,nbfoc
              dervi(indk+j,indk+i) = dervi(indk+3+j,indk+3+i)
            enddo
          enddo
        endif
      else
!#################################################################################
!  Build second derivative matrix in dervi for inversion as this is overwritten  #
!#################################################################################
        if (lrigid) then
#ifdef MPI
          if (nprocs.gt.1) then
!-----------------------------------
!  Rigid molecule case - parallel  |
!-----------------------------------
            allocate(dtmp(3*numatnomolfree,3),stat=status)
            if (status/=0) call outofmemory('property3','dtmp')
            allocate(StatMPI(MPI_Status_Size),stat=status)
            if (status/=0) call outofmemory('property3','StatMPI')
!
!  Copy individual atoms
!
            do i = 1,numatnomolfree
!
!  Is derv2 local to this node?
!
              snode = atom2node(numatnomolptr(i))
              if (snode.eq.procid) then
                ld2loc = .true.
                indif = 3*(atom2local(numatnomolptr(i)) - 1)
              else
                ld2loc = .false.
              endif
!
!  Is dervi local to this node?
!
              rnode = pdinode(i)
              if (rnode.eq.procid) then
                ldiloc = .true.
                iloc = pdig2l(i)
                indil = 3*(iloc-1)
              else
                ldiloc = .false.
              endif
!
!  If derv2 and dervi are both local then copy
!
              if (ld2loc.and.ldiloc) then
                do j = 1,numatnomolfree
                  indjl = 3*(j-1)
                  indjf = 3*(numatnomolptr(j)-1)
                  do ix = 1,3
                    do jx = 1,3
                      dervi(indjl+jx,indil+ix) = derv2(indjf+jx,indif+ix)
                    enddo
                  enddo
                enddo
              elseif (ldiloc) then
!
!  If dervi is local then receive
!
                nsize = 9*numatnomolfree
                ntag = i
                call MPI_IRecv(dtmp,nsize,MPI_double_precision,snode,ntag,MPI_Comm_World,Request,MPIerror)
              elseif (ld2loc) then
!
!  If derv2 is local then send
!
                do j = 1,numatnomolfree
                  indjl = 3*(j-1)
                  indjf = 3*(numatnomolptr(j)-1)
                  do ix = 1,3
                    do jx = 1,3
                      dtmp(indjl+jx,ix) = derv2(indjf+jx,indif+ix)
                    enddo
                  enddo
                enddo
                nsize = 9*numatnomolfree
                ntag = i
                call MPI_ISend(dtmp,nsize,MPI_double_precision,snode,ntag,MPI_Comm_World,Request,MPIerror)
              endif
!
!  Wait for communication to finish
!
              if ((ld2loc.and..not.ldiloc).or.(ldiloc.and..not.ld2loc)) then
                call MPI_WaitAll(1,Request,StatMPI,MPIerror)
                if (ldiloc) then
!
!  After communication place data in dervi
!
                  do j = 1,numatnomolfree
                    indjl = 3*(j-1)
                    do ix = 1,3
                      do jx = 1,3
                        dervi(indjl+jx,indil+ix) = dtmp(indjl+jx,ix)
                      enddo
                    enddo
                  enddo
                endif
              endif
            enddo
!
            deallocate(StatMPI,stat=status)
            if (status/=0) call deallocate_error('property3','StatMPI')
            deallocate(dtmp,stat=status)
            if (status/=0) call deallocate_error('property3','dtmp')
!
            mind = 3*numatnomolfree
!
!  Copy breathing shells
!
! DEBUG - breathing shells not handled here for rigid molecules!
!
!  Copy molecule terms
!
            do i = 1,numatnomolfreeonnode
              indil = 3*(i-1)
              indif = 3*(numatnomolonnodeptr(i)-1)
              do j = 1,nmolfree
                if (lnorotate) then
                  mindj = mind + 3*(j-1)
                else
                  mindj = mind + 6*(j-1)
                endif
!
!  Copy rigid molecule - atom terms
!
                do ix = 1,3
                  do jx = 1,3
                    dervi(mindj+jx,indil+ix) = molTCdrv(indif+ix,jx,j)
                  enddo
                enddo
                if (.not.lnorotate) then
                  do ix = 1,3
                    do jx = 1,3
                      dervi(mindj+3+jx,indil+ix) = molQCdrv(indif+ix,jx,j)
                    enddo
                  enddo
                endif
              enddo
            enddo
!
            do i = 1,nmolfree
              if (lnorotate) then
                iloc = numatnomolfree + i
              else
                iloc = numatnomolfree + 2*(i-1) + 1
              endif
              if (pdinode(iloc).eq.procid) then
                if (lnorotate) then
                  mindi = 3*(pdig2l(iloc)-1)
                else
                  mindi = 3*(pdig2l(iloc)-1)
                endif
!
!  Copy rigid molecule - atom terms for translation
!
                do j = 1,numatnomolfree
                  indjl = 3*(j-1)
                  indjf = 3*(numatnomolptr(j)-1)
                  do ix = 1,3
                    do jx = 1,3
                      dervi(indjl+jx,mindi+ix) = molTCdrv(indjf+jx,ix,i)
                    enddo
                  enddo
                enddo
!
!  Copy rigid molecule - rigid molecule terms
!
                do j = 1,nmolfree
                  if (lnorotate) then
                    mindj = mind + 3*(j-1)
                  else
                    mindj = mind + 6*(j-1)
                  endif
                  do ix = 1,3
                    do jx = 1,3
                      dervi(mindj+jx,mindi+ix) = molTTdrv(jx,ix,j,i)
                    enddo
                  enddo
                  if (.not.lnorotate) then
                    do ix = 1,3
                      do jx = 1,3
                        dervi(mindj+3+jx,mindi+ix) = molQTdrv(jx,ix,j,i)
                      enddo
                    enddo
                  endif
                enddo
              endif
              if (.not.lnorotate.and.pdinode(iloc+1).eq.procid) then
                mindi = 3*(pdig2l(iloc+1)-1)
!
!  Copy rigid molecule - atom terms for rotation
!
                do j = 1,numatnomolfree
                  indjl = 3*(j-1)
                  indjf = 3*(numatnomolptr(j)-1)
                  do ix = 1,3
                    do jx = 1,3
                      dervi(indjl+jx,mindi+ix) = molQCdrv(indjf+jx,ix,i)
                    enddo
                  enddo
                enddo
!
!  Copy rigid molecule - rigid molecule terms
!
                do j = 1,nmolfree
                  mindj = mind + 6*(j-1)
                  do ix = 1,3
                    do jx = 1,3
                      dervi(mindj+jx,mindi+ix) = molQTdrv(ix,jx,i,j)
                      dervi(mindj+3+jx,mindi+ix) = molQQdrv(jx,ix,j,i)
                    enddo
                  enddo
                enddo
              endif
            enddo
!
!  If a molecule has been fixed then add quaternions for final molecule
!  NB: Final row for quaternions will be in mindi+1 to mindi+3 as there is no translation
!
            if (nmol.gt.nmolfree.and..not.lnorotate) then
              iloc = numatnomolfree + 2*nmolfree + 1
              if (pdinode(iloc).eq.procid) then
                mindi = 3*(pdig2l(iloc)-1)
!
!  Copy rigid molecule - atom terms
!
                do j = 1,numatnomolfree
                  indjl = 3*(j-1)
                  indjf = 3*(numatnomolptr(j)-1)
                  do ix = 1,3
                    do jx = 1,3
                      dervi(indjl+jx,mindi+ix) = molQCdrv(indjf+jx,ix,nmol)
                    enddo
                  enddo
                enddo
!
!  Copy rigid molecule - rigid molecule terms
!
                do j = 1,nmolfree
                  mindj = mind + 6*(j-1)
                  do ix = 1,3
                    do jx = 1,3
                      dervi(mindj+jx,mindi+ix) = molQTdrv(ix,jx,nmol,j)
                      dervi(mindj+3+jx,mindi+ix) = molQQdrv(jx,ix,j,nmol)
                    enddo
                  enddo
                enddo
!
!  Add self terms for last molecule - only molQQdrv as it has no translation
!
                mindj = mind + 6*(nmol-1)
                do ix = 1,3
                  do jx = 1,3
                    dervi(mindj+jx,mindi+ix) = molQQdrv(jx,ix,nmol,nmol)
                  enddo
                enddo
              endif
!
              mindj = 3*(iloc-1)
!
!  Copy rigid molecule - atom terms
!
              do i = 1,numatnomolfreeonnode
                indil = 3*(i-1)
                indif = 3*(numatnomolptr(i)-1)
                do ix = 1,3
                  do jx = 1,3
                    dervi(mindj+jx,indil+ix) = molQCdrv(indif+ix,jx,nmol)
                  enddo
                enddo
              enddo
!
!  Copy rigid molecule - rigid molecule terms
!
              do i = 1,nmolfree
                iloc = numatnomolfree + 2*(i-1) + 1
                if (pdinode(iloc).eq.procid) then
                  mindi = 3*(pdig2l(iloc)-1)
                  do ix = 1,3
                    do jx = 1,3
                      dervi(mindj+jx,mindi+ix) = molQTdrv(jx,ix,nmol,i)
                    enddo
                  enddo
                endif
                if (pdinode(iloc+1).eq.procid) then
                  mindi = 3*(pdig2l(iloc+1)-1)
                  do ix = 1,3
                    do jx = 1,3
                      dervi(mindj+jx,mindi+ix) = molQQdrv(jx,ix,nmol,i)
                    enddo
                  enddo
                endif
              enddo
            endif
          else
#endif
!---------------------------------
!  Rigid molecule case - serial  |
!---------------------------------
!
!  Copy individual atoms
!
            do i = 1,numatnomolfree
              indil = 3*(i-1)
              indif = 3*(numatnomolptr(i)-1)
              do j = 1,numatnomolfree
                indjl = 3*(j-1)
                indjf = 3*(numatnomolptr(j)-1)
                do ix = 1,3
                  do jx = 1,3
                    dervi(indjl+jx,indil+ix) = derv2(indjf+jx,indif+ix)
                  enddo
                enddo
              enddo
            enddo
            mind = 3*numatnomolfree
!
!  Copy breathing shells
!
! DEBUG - breathing shells not handled here for rigid molecules!
!
!  Copy molecule terms
!
            do i = 1,nmolfree
              if (lnorotate) then
                mindi = mind + 3*(i-1)
              else
                mindi = mind + 6*(i-1)
              endif
!
!  Copy rigid molecule - atom terms
!
              do j = 1,numatnomolfree
                indjl = 3*(j-1)
                indjf = 3*(numatnomolptr(j)-1)
                do ix = 1,3
                  do jx = 1,3
                    dervi(indjl+jx,mindi+ix) = molTCdrv(indjf+jx,ix,i)
                    dervi(mindi+ix,indjl+jx) = molTCdrv(indjf+jx,ix,i)
                  enddo
                enddo
                if (.not.lnorotate) then
                  do ix = 1,3
                    do jx = 1,3
                      dervi(indjl+jx,mindi+3+ix) = molQCdrv(indjf+jx,ix,i)
                      dervi(mindi+3+ix,indjl+jx) = molQCdrv(indjf+jx,ix,i)
                    enddo
                  enddo
                endif
              enddo
!
!  Copy rigid molecule - rigid molecule terms
!
              do j = 1,nmolfree
                if (lnorotate) then
                  mindj = mind + 3*(j-1)
                else
                  mindj = mind + 6*(j-1)
                endif
                do ix = 1,3
                  do jx = 1,3
                    dervi(mindj+jx,mindi+ix) = molTTdrv(jx,ix,j,i)
                  enddo
                enddo
                if (.not.lnorotate) then
                  do ix = 1,3
                    do jx = 1,3
                      dervi(mindj+3+jx,mindi+ix) = molQTdrv(jx,ix,j,i)
                      dervi(mindj+jx,mindi+3+ix) = molQTdrv(ix,jx,i,j)
                      dervi(mindj+3+jx,mindi+3+ix) = molQQdrv(jx,ix,j,i)
                    enddo
                  enddo
                endif
              enddo
            enddo
!
!  If a molecule has been fixed then add quaternions for final molecule
!  NB: Final row for quaternions will be in mindi+1 to mindi+3 as there is no translation
!
            if (nmol.gt.nmolfree.and..not.lnorotate) then
              mindi = mind + 6*(nmol-1)
!
!  Copy rigid molecule - atom terms
!
              do j = 1,numatnomolfree
                indjl = 3*(j-1)
                indjf = 3*(numatnomolptr(j)-1)
                do ix = 1,3
                  do jx = 1,3
                    dervi(indjl+jx,mindi+ix) = molQCdrv(indjf+jx,ix,nmol)
                    dervi(mindi+ix,indjl+jx) = molQCdrv(indjf+jx,ix,nmol)
                  enddo
                enddo
              enddo
!
!  Copy rigid molecule - rigid molecule terms
!
              do j = 1,nmolfree
                mindj = mind + 6*(j-1)
                do ix = 1,3
                  do jx = 1,3
                    dervi(mindj+jx,mindi+ix) = molQTdrv(ix,jx,nmol,j)
                    dervi(mindi+ix,mindj+jx) = molQTdrv(ix,jx,nmol,j)
                    dervi(mindj+3+jx,mindi+ix) = molQQdrv(jx,ix,j,nmol)
                    dervi(mindi+ix,mindj+3+jx) = molQQdrv(jx,ix,j,nmol)
                  enddo
                enddo
              enddo
!
!  Add self terms for last molecule - only molQQdrv as it has no translation
!
              do ix = 1,3
                do jx = 1,3
                  dervi(mindi+jx,mindi+ix) = molQQdrv(jx,ix,nmol,nmol)
                enddo
              enddo
            endif
#ifdef MPI
          endif
#endif
!--------------------------------------
!  Set size of matrix to diagonalise  |
!--------------------------------------
          if (lnorotate) then
            ninv = mind + 3*nmolfree
          else
            ninv = mind + 3*nmolfree + 3*nmol
          endif
!
          if (.not.lnorotate) then
!----------------------------
!  Handle linear molecules  |
!----------------------------
            do i = 1,nmol
              nmOK(i) = 0
              do ix = 1,3
!
!  Are moments of inertia zero or there are no angular forces?
!
                if (abs(molI(ix,i)).gt.too_small.and.abs(molQQdrv(ix,ix,i,i)).gt.0.001_dp) then
                  nmOK(i) = nmOK(i) + 1
                  nmOKptr(nmOK(i),i) = ix
                else
#ifdef MPI
!
!  If so then zero second derivatives and place 1 on the diagonal
!
                  if (nprocs.gt.1) then
                    if (i.gt.nmolfree) then
                      iloc = numatnomolfree + 2*(i-1) + 1
                    else
                      iloc = numatnomolfree + 2*(i-1) + 2
                    endif
                    if (i.gt.nmolfree) then
                      mindj = mind + 6*(i-1) + ix
                    else
                      mindj = mind + 6*(i-1) + 3 + ix
                    endif
                    dervi(mindj,1:ndi) = 0.0_dp
                    if (pdinode(iloc).eq.procid) then
                      mindi = 3*(pdig2l(iloc)-1)
                      dervi(1:ninv,mindi+ix) = 0.0_dp
                      dervi(mindj,mindi+ix) = 1.0_dp
                    endif
                  else
#endif
                    if (i.gt.nmolfree) then
                      mindi = mind + 6*(i-1) + ix
                    else
                      mindi = mind + 6*(i-1) + 3 + ix
                    endif
                    dervi(1:ninv,mindi) = 0.0_dp
                    dervi(mindi,1:ninv) = 0.0_dp
                    dervi(mindi,mindi) = 1.0_dp
#ifdef MPI
                  endif
#endif
                endif
              enddo
            enddo
          endif
        else
!---------------------------------------
!  Old algorithm - no rigid molecules  |
!---------------------------------------
          do i = 1,n3fu
            do j = 1,n3f
              dervi(j,i) = derv2(j,i)
            enddo
          enddo
!
!  Radial component
!
          if (nbs.gt.0) then
            do i = 1,nbs
              iptr = n3f + nbsptr(i)
              do j = 1,nbs
                jptr = n3f + nbsptr(j)
                dervi(mint+j,mint+i) = derv2(jptr,iptr)
              enddo
              do j = 1,mint
                dervi(j,mint+i) = derv2(j,iptr)
              enddo
            enddo
          endif
          ninv = mint + nbs
        endif
      endif
!
!  Invert internal derivative matrix:
!
!  Ignore second derivatives of last atom to prevent a singularity - corresponds
!  to removing the 3 translational degrees of freedom of the lattice
!
      ifail = 1
      call matrix_inversion_library(ninv,1_i4,maxd2,3_i4*nblocksize,dervi,0_i4,ifail)
!
      if (ifail.gt.0) then
        nwarn = nwarn + 1
        call outwarning('Properties cannot be calculated - matrix is singular',0_i4)
        goto 999
      endif
    endif
  endif
!
!  Allocate tmp2 workspace for elastic terms also used for piezoelectric constants
!
  if ((lelastic.or.lpiezo).and.linternal) then
    if (lpocc) then
      allocate(tmp2(ninv,6),stat=status)
      if (status/=0) call outofmemory('property3','tmp2')
    elseif (nbs.gt.0) then
      allocate(tmp2(mint+nbs,6),stat=status)
      if (status/=0) call outofmemory('property3','tmp2')
    else
      allocate(tmp2(ninv,6),stat=status)
      if (status/=0) call outofmemory('property3','tmp2')
    endif
  endif
!*********************
!  Elastic constants *
!*********************
  if (lelastic) then
    elcon(1:6,1:6) = sderv2(1:6,1:6)
!
!  Correct for pressure term in the second derivatives
!
!  First remove enthalpy terms.
!  Second apply non-zero pressure corrections
!
    if (abs(press).gt.0.0_dp) then
      epv = press*vol*cfactor
      elcon(1,1) = elcon(1,1) - epv
      elcon(2,2) = elcon(2,2) - epv
      elcon(3,3) = elcon(3,3) - epv
      elcon(2,1) = elcon(2,1) - epv
      elcon(3,1) = elcon(3,1) - epv
      elcon(3,2) = elcon(3,2) - epv
      elcon(1,2) = elcon(1,2) - epv
      elcon(1,3) = elcon(1,3) - epv
      elcon(2,3) = elcon(2,3) - epv
      elcon(4,4) = elcon(4,4) - 0.5_dp*epv
      elcon(5,5) = elcon(5,5) - 0.5_dp*epv
      elcon(6,6) = elcon(6,6) - 0.5_dp*epv
! non-zero pressure correction
! infinitesimal (small) strain correction at finite hydostatic pressure
! Barron & Klein Proc. Phys. Soc.vol 85 1965 Eqn.5.5
! 1/2P(2DijDkl-DilDjk-DikDjl      
! Cijkl=1/V(d2U/dEijdEkl)+P/2*(2DijDkl-DilDjk-DikDjl) pressure correction
      elcon(2,1) = elcon(2,1) + epv
      elcon(3,1) = elcon(3,1) + epv
      elcon(3,2) = elcon(3,2) + epv
      elcon(1,2) = elcon(1,2) + epv
      elcon(1,3) = elcon(1,3) + epv
      elcon(2,3) = elcon(2,3) + epv
      elcon(4,4) = elcon(4,4) - 0.5_dp*epv
      elcon(5,5) = elcon(5,5) - 0.5_dp*epv
      elcon(6,6) = elcon(6,6) - 0.5_dp*epv
    endif
!
    if (linternal) then
      if (lpocc) then
        do i = 1,6
          do j = 1,ninv
            tmp2(j,i) = 0.0_dp
            do k = 1,ninv
              tmp2(j,i) = tmp2(j,i) + dervi(k,j)*d3(k,i)
            enddo
          enddo
          do j = 1,6
            do k = 1,ninv
              elcon(j,i) = elcon(j,i) - d3(k,j)*tmp2(k,i)
            enddo
          enddo
        enddo
      elseif (nbs.gt.0) then
        do i = 1,6
          do j = 1,mint+nbs
            tmp2(j,i) = 0.0_dp
            do k = 1,mint+nbs
              if (k.gt.mint) then
                ki = n3f + nbsptr(k-mint)
              else
                ki = k
              endif
              tmp2(j,i) = tmp2(j,i) + dervi(k,j)*derv3(ki,i)
            enddo
          enddo
          do j = 1,6
            do k = 1,mint+nbs
              if (k.gt.mint) then
                ki = n3f + nbsptr(k-mint)
              else
                ki = k
              endif
              elcon(j,i) = elcon(j,i) - derv3(ki,j)*tmp2(k,i)
            enddo
          enddo
        enddo
      else
        if (nprocs.gt.1) then
!-------------
!  Parallel  |
!-------------
          allocate(tmp2s(ninv,6),stat=status)
          if (status/=0) call outofmemory('property3','tmp2s')
!
#ifdef MPI
          if (lrigid) then
!
!  Create a copy of derv3 with only numatnomolfree component
!
            d3nomol(1:mind,1:6) = 0.0_dp
            do i = 1,numatnomolfree
              snode = atom2node(numatnomolptr(i))
              if (snode.eq.procid) then
                indif = 3*(atom2local(numatnomolptr(i)) - 1)
                indil = 3*(i-1)
                do ix = 1,3
                  d3nomol(indil+ix,1:6) = derv3(indif+ix,1:6)
                enddo
              endif
            enddo
!
            call sumall(d3nomol,d3nomols,6_i4*mind,"property3","d3nomol")
!
            tmp2(1:ninv,1:6) = 0.0_dp
            do i = 1,6
              do j = 1,numatnomolfree
                jloc = pdig2l(j)
                if (jloc.gt.0) then
                  indjl = 3*(jloc-1)
                  indjf = 3*(j-1)
                  do jx = 1,3
                    do k = 1,ninv
                      tmp2(k,i) = tmp2(k,i) + dervi(k,indjl+jx)*d3nomols(indjf+jx,i)
                    enddo
                  enddo
                endif
              enddo
              if (lnorotate) then
                do j = 1,nmolfree
                  jloc = pdig2l(numatnomolfree+j)
                  if (jloc.gt.0) then
                    mindj = 3*(jloc-1)
                    do jx = 1,3
                      do k = 1,ninv
                        tmp2(k,i) = tmp2(k,i) + dervi(k,mindj+jx)*molTSdrv(i,jx,j)
                      enddo
                    enddo
                  endif
                enddo
              else
                do j = 1,nmolfree
                  jloc = pdig2l(numatnomolfree+2*(j-1)+1)
                  if (jloc.gt.0) then
                    mindj = 3*(jloc-1)
                    do jx = 1,3
                      do k = 1,ninv
                        tmp2(k,i) = tmp2(k,i) + dervi(k,mindj+jx)*molTSdrv(i,jx,j)
                      enddo
                    enddo
                  endif
                  jloc = pdig2l(numatnomolfree+2*(j-1)+2)
                  if (jloc.gt.0) then
                    mindj = 3*(jloc-1)
                    do jj = 1,nmOK(j)
                      jx = nmOKptr(jj,j)
                      do k = 1,ninv
                        tmp2(k,i) = tmp2(k,i) + dervi(k,mindj+jx)*molQSdrv(i,jx,j)
                      enddo
                    enddo
                  endif
                enddo
              endif
              if (nmol.gt.nmolfree.and..not.lnorotate) then
                jloc = pdig2l(numatnomolfree+nmolfree+nmol)
                if (jloc.gt.0) then
                  mindj = 3*(jloc-1)
                  do jj = 1,nmOK(nmol)
                    jx = nmOKptr(jj,nmol)
                    do k = 1,ninv
                      tmp2(k,i) = tmp2(k,i) + dervi(k,mindj+jx)*molQSdrv(i,jx,nmol)
                    enddo
                  enddo
                endif
              endif
            enddo
!
            call sumall(tmp2,tmp2s,6_i4*ninv,"property3","tmp2")
!
            eltmp(1:6,1:6) = 0.0_dp
!
            do i = 1,6
              do j = 1,6
                do k = 1,numatnomolfree
                  if (pdinode(k).eq.procid) then
                    indkl = 3*(k-1)
                    do kx = 1,3
                      eltmp(j,i) = eltmp(j,i) - tmp2s(indkl+kx,i)*d3nomols(indkl+kx,j)
                    enddo
                  endif
                enddo
                if (lnorotate) then
                  do k = procid+1,nmolfree,nprocs
                    mindk = mind + 3*(k-1)
                    do kx = 1,3
                      eltmp(j,i) = eltmp(j,i) - tmp2s(mindk+kx,i)*molTSdrv(j,kx,k)
                    enddo
                  enddo
                else
                  do k = procid+1,nmolfree,nprocs
                    mindk = mind + 6*(k-1)
                    do kx = 1,3
                      eltmp(j,i) = eltmp(j,i) - tmp2s(mindk+kx,i)*molTSdrv(j,kx,k)
                    enddo
                    do kk = 1,nmOK(k)
                      kx = nmOKptr(kk,k)
                      eltmp(j,i) = eltmp(j,i) - tmp2s(mindk+3+kx,i)*molQSdrv(j,kx,k)
                    enddo
                  enddo
                endif
                if (nmol.gt.nmolfree.and..not.lnorotate.and.ioproc) then
                  mindk = mind + 6*(nmol-1)
                  do kk = 1,nmOK(nmol)
                    kx = nmOKptr(kk,nmol)
                    eltmp(j,i) = eltmp(j,i) - tmp2s(mindk+kx,i)*molQSdrv(j,kx,nmol)
                  enddo
                endif
              enddo
            enddo
!
            call sumall(eltmp,eltmp2,6_i4*6_i4,"property3","eltmp")
!
            do i = 1,6
              do j = 1,6
                elcon(j,i) = elcon(j,i) + eltmp2(j,i)
              enddo
            enddo
          else
#endif
            tmp2(1:mint,1:6) = 0.0_dp
!
            do i = 1,6
              do j = 1,mintu
                do k = 1,mint
                  tmp2(k,i) = tmp2(k,i) + dervi(k,j)*derv3(j,i)
                enddo
              enddo
            enddo
!
            call sumall(tmp2,tmp2s,6_i4*mint,"property3","tmp2")
!
            eltmp(1:6,1:6) = 0.0_dp
!
            do i = 1,6
              do j = 1,6
                indk = 0
                do k = 1,natomsonnodem1
                  ki = 3*(node2atom(k) - 1)
                  do ii = 1,3
                    ki = ki + 1
                    indk = indk + 1
                    eltmp(j,i) = eltmp(j,i) - derv3(indk,j)*tmp2s(ki,i)
                  enddo
                enddo
              enddo
            enddo
!
            call sumall(eltmp,eltmp2,6_i4*6_i4,"property3","eltmp")
!
            do i = 1,6
              do j = 1,6
                elcon(j,i) = elcon(j,i) + eltmp2(j,i)
              enddo
            enddo
#ifdef MPI
          endif
#endif
!
          deallocate(tmp2s,stat=status)
          if (status/=0) call deallocate_error('property3','tmp2s')
        else
!-----------
!  Serial  |
!-----------
          if (lrigid) then
            do i = 1,6
              tmp2(1:ninv,i) = 0.0_dp
              do j = 1,numatnomolfree
                indjl = 3*(j-1)
                indjf = 3*(numatnomolptr(j)-1)
                do jx = 1,3
                  do k = 1,ninv
                    tmp2(k,i) = tmp2(k,i) + dervi(k,indjl+jx)*derv3(indjf+jx,i)
                  enddo
                enddo
              enddo
              if (lnorotate) then
                do j = 1,nmolfree
                  mindj = mind + 3*(j-1)
                  do jx = 1,3
                    do k = 1,ninv
                      tmp2(k,i) = tmp2(k,i) + dervi(k,mindj+jx)*molTSdrv(i,jx,j)
                    enddo
                  enddo
                enddo
              else
                do j = 1,nmolfree
                  mindj = mind + 6*(j-1)
                  do jx = 1,3
                    do k = 1,ninv
                      tmp2(k,i) = tmp2(k,i) + dervi(k,mindj+jx)*molTSdrv(i,jx,j)
                    enddo
                  enddo
                  do jj = 1,nmOK(j)
                    jx = nmOKptr(jj,j)
                    do k = 1,ninv
                      tmp2(k,i) = tmp2(k,i) + dervi(k,mindj+3+jx)*molQSdrv(i,jx,j)
                    enddo
                  enddo
                enddo
              endif
              if (nmol.gt.nmolfree.and..not.lnorotate) then
                mindj = mind + 6*(nmol-1)
                do jj = 1,nmOK(nmol)
                  jx = nmOKptr(jj,nmol)
                  do k = 1,ninv
                    tmp2(k,i) = tmp2(k,i) + dervi(k,mindj+jx)*molQSdrv(i,jx,nmol)
                  enddo
                enddo
              endif
              do j = 1,6
                do k = 1,numatnomolfree
                  indkl = 3*(k-1)
                  indkf = 3*(numatnomolptr(k)-1)
                  do kx = 1,3
                    elcon(j,i) = elcon(j,i) - tmp2(indkl+kx,i)*derv3(indkf+kx,j)
                  enddo
                enddo
                if (lnorotate) then
                  do k = 1,nmolfree
                    mindk = mind + 3*(k-1)
                    do kx = 1,3
                      elcon(j,i) = elcon(j,i) - tmp2(mindk+kx,i)*molTSdrv(j,kx,k)
                    enddo
                  enddo
                else
                  do k = 1,nmolfree
                    mindk = mind + 6*(k-1)
                    do kx = 1,3
                      elcon(j,i) = elcon(j,i) - tmp2(mindk+kx,i)*molTSdrv(j,kx,k)
                    enddo
                    do kk = 1,nmOK(k)
                      kx = nmOKptr(kk,k)
                      elcon(j,i) = elcon(j,i) - tmp2(mindk+3+kx,i)*molQSdrv(j,kx,k)
                    enddo
                  enddo
                endif
                if (nmol.gt.nmolfree.and..not.lnorotate) then
                  mindk = mind + 6*(nmol-1)
                  do kk = 1,nmOK(nmol)
                    kx = nmOKptr(kk,nmol)
                    elcon(j,i) = elcon(j,i) - tmp2(mindk+kx,i)*molQSdrv(j,kx,nmol)
                  enddo
                endif
              enddo
            enddo
          else
            do i = 1,6
              do j = 1,mint
                tmp2(j,i) = 0.0_dp
                do k = 1,mint
                  tmp2(j,i) = tmp2(j,i) + dervi(k,j)*derv3(k,i)
                enddo
              enddo
              do j = 1,6
                do k = 1,mint
                  elcon(j,i) = elcon(j,i) - derv3(k,j)*tmp2(k,i)
                enddo
              enddo
            enddo
          endif
        endif
      endif
      rvol = 10.0_dp*conve/vol
      do i = 1,6
        do j = 1,6
          elcon(j,i) = rvol*elcon(j,i)
        enddo
      enddo
    else
      rvol = 10.0_dp*conve/vol
      do i = 1,6
        do j = 1,6
          elcon(j,i) = rvol*elcon(j,i)
        enddo
      enddo
    endif
!
!  Compute eigenvalues of elastic constant tensor
!
    elconvec(1:6,1:6) = elcon(1:6,1:6)
    call dsyev('V','U',6_i4,elconvec,6_i4,elconeig,eltmp2,36_i4,ifail)
    lelconeigok = (ifail.eq.0)
    nelconeig = 0
    do i = 1,6
      if (elconeig(i).lt.1.0d-12) nelconeig = nelconeig + 1
    enddo
  endif
  if (lpiezs.or.lmoduli) then
!**********************
!  Compliance tensor  *
!**********************
    do i = 1,6
      do j = 1,6
        compliances(j,i) = elcon(j,i)
      enddo
    enddo
    ifail = 1
    call matrix_inversion(compliances,6_i4,6_i4,tmp,ifail)
    if (ifail.gt.0) then
      lcompliances = .false.
      lpiezs = .false.
    endif
  endif
!****************************
!  Piezoelectric constants  *
!****************************
  if (lpiezo) then
    if (linternal) then
!
!  Piezoelectric stress constants
!
      do i = 1,6
        do j = 1,3
          sum = 0.0_dp
          do k = 1,nqfo
            indk = 3*(k-1) + j
            sum = sum + tmp2(indk,i)*qfo(k)
          enddo
          piezo(i,j) = - piefactor*sum
        enddo
      enddo
      if (nprocs.gt.1) then
        call sumall(piezo,eltmp,6_i4*3_i4,"property3","piezo")
        piezo(1:6,1:3) = eltmp(1:6,1:3)
      endif
    else
!
!  If not linternal then piezo must be zero
!
      piezo(1:6,1:3) = 0.0_dp
    endif
    if (lpiezs) then
!
!  Piezoelectric strain constants
!
      do i = 1,3
        do j = 1,6
          piezs(j,i) = 0.0_dp
          do k = 1,6
            piezs(j,i) = piezs(j,i) + compliances(k,j)*piezo(k,i)
          enddo
        enddo
      enddo
    else
      piezs(1:6,1:3) = 0.0_dp
    endif
  endif
  if (lmoduli) then
!
!  Bulk modulus - calculate all 3 possible definitions
!
    bulkmod_hill = 0.0_dp
    bulkmod_reuss = 0.0_dp
    bulkmod_voigt = 0.0_dp
    do i = 1,3
      bulkmod_voigt = bulkmod_voigt + elcon(1,i) + elcon(2,i) + elcon(3,i)
    enddo
    bulkmod_voigt = bulkmod_voigt/9.0_dp
    if (lcompliances) then
      do i = 1,3
        bulkmod_reuss = bulkmod_reuss + compliances(1,i) + compliances(2,i) + compliances(3,i)
      enddo
      if (abs(bulkmod_reuss).gt.1.0d-8) then
        bulkmod_reuss = 1.0_dp/bulkmod_reuss
      else
        bulkmod_reuss = 0.0_dp
      endif
    else
      bulkmod_reuss = 0.0_dp
    endif
    bulkmod_hill = 0.5_dp*(bulkmod_reuss + bulkmod_voigt)
!
!  Shear modulus - calculate all 3 possible definitions
!
    shearmod_voigt = elcon(1,1) + elcon(2,2) + elcon(3,3) &
      + 3.0_dp*(elcon(4,4) + elcon(5,5) + elcon(6,6)) &
      - elcon(1,2) - elcon(1,3) - elcon(2,3)
    shearmod_voigt = shearmod_voigt/15.0_dp
    if (lcompliances) then
      shearmod_reuss = 4.0_dp*(compliances(1,1) + compliances(2,2) + compliances(3,3)) &
                     - 4.0_dp*(compliances(1,2) + compliances(1,3) + compliances(2,3)) &
                     + 3.0_dp*(compliances(4,4) + compliances(5,5) + compliances(6,6))
      if (abs(shearmod_reuss).gt.1.0d-8) then
        shearmod_reuss = 15.0_dp/shearmod_reuss
      else
        shearmod_reuss = 0.0_dp
      endif
    else
      shearmod_reuss = 0.0_dp
    endif
    shearmod_hill = 0.5_dp*(shearmod_reuss + shearmod_voigt)
!
!  Return the requested moduli values for fitting
!
    if (index(keyword,'voi').ne.0) then
      bulkmod = bulkmod_voigt
      shearmod = shearmod_voigt
    elseif (index(keyword,'hill').ne.0) then
      bulkmod = bulkmod_hill
      shearmod = shearmod_hill
    else
      bulkmod = bulkmod_reuss
      shearmod = shearmod_reuss
    endif
!
!  Acoustic wave velocities - polycrystalline average
!
    r43 = 4.0_dp/3.0_dp
    vs_hill = sign(sqrt(abs(shearmod_hill)/density),shearmod_hill)
    vs_reuss = sign(sqrt(abs(shearmod_reuss)/density),shearmod_reuss)
    vs_voigt = sign(sqrt(abs(shearmod_voigt)/density),shearmod_voigt)
    vp_hill = sqrt(abs(bulkmod_hill + r43*shearmod_hill)/density)
    vp_reuss = sqrt(abs(bulkmod_reuss + r43*shearmod_reuss)/density)
    vp_voigt = sqrt(abs(bulkmod_voigt + r43*shearmod_voigt)/density)
!
!  Young's Moduli
!
    if (abs(compliances(1,1)).gt.1.0d-8) then
      ym(1) = 1.0_dp/compliances(1,1)
    else
      ym(1) = 0.0_dp
    endif
    if (abs(compliances(2,2)).gt.1.0d-8) then
      ym(2) = 1.0_dp/compliances(2,2)
    else
      ym(2) = 0.0_dp
    endif
    if (abs(compliances(3,3)).gt.1.0d-8) then
      ym(3) = 1.0_dp/compliances(3,3)
    else
      ym(3) = 0.0_dp
    endif
!
!  Poisson's ratios
!
    poissonratio(1) = - compliances(2,1)*ym(2)
    poissonratio(2) = - compliances(3,1)*ym(3)
    poissonratio(3) = - compliances(3,2)*ym(3)
  endif
!******************************
!  Static dielectric constant *
!******************************
  if (ldlcs) then
    if (lpocc) then
      do i = 1,3
        do j = 1,3
          sum = 0.0_dp
          do k = 1,nqfo
            qak = qfo(k)
            indk = 3*(k-1) + i
            do l = 1,nqfo
              indl = 3*(l-1) + j
              sum = sum + dervi(indl,indk)*qak*qfo(l)
            enddo
          enddo
          dicons(j,i) = dlcfactor*sum
        enddo
      enddo
    else
#ifdef MPI
      if (nprocs.gt.1) then
        if (lrigid) then
          do i = 1,3
            do j = 1,3
              sum = 0.0_dp
              do k = 1,nqfo
                if (pdinode(k).eq.procid) then
                  ki = pdig2l(k)
                  qak = qfo(k)
                  indk = 3*(ki-1) + i
                  do l = 1,nqfo
                    indl = 3*(l-1) + j
                    sum = sum + dervi(indl,indk)*qak*qfo(l)
                  enddo
                endif
              enddo
              dicons(j,i) = dlcfactor*sum
            enddo
          enddo
        else
          do i = 1,3
            do j = 1,3
              sum = 0.0_dp
              do k = 1,natomsonnodem1
                ki = node2atom(k)
                qak = qfo(ki)
                indk = 3*(k-1) + i
                do l = 1,nqfo
                  indl = 3*(l-1) + j
                  sum = sum + dervi(indl,indk)*qak*qfo(l)
                enddo
              enddo
              dicons(j,i) = dlcfactor*sum
            enddo
          enddo
        endif
      else
#endif
        do i = 1,3
          do j = 1,3
            sum = 0.0_dp
            do k = 1,nqfo
              qak = qfo(k)
              indk = 3*(k-1) + i
              do l = 1,nqfo
                indl = 3*(l-1) + j
                sum = sum + dervi(indl,indk)*qak*qfo(l)
              enddo
            enddo
            dicons(j,i) = dlcfactor*sum
          enddo
        enddo
#ifdef MPI
      endif
#endif
    endif
    if (nprocs.gt.1) then
      call sumall(dicons,ditmp,3_i4*3_i4,"property3","ditmp")
      dicons(1:3,1:3) = ditmp(1:3,1:3)
    endif
    do i = 1,3
      dicons(i,i) = dicons(i,i) + 1.0_dp
    enddo
!
!  Calculate refractive indices
!
    do i = 1,3
      do j = 1,3
        r33(j,i) = dicons(j,i)
      enddo
    enddo
    call dsyev('N','U',3_i4,r33,3_i4,r3,w1l,9_i4,ifail)
    srefind(1) = sqrt(abs(r3(1)))
    srefind(2) = sqrt(abs(r3(2)))
    srefind(3) = sqrt(abs(r3(3)))
    srefind(1) = sign(srefind(1),r3(1))
    srefind(2) = sign(srefind(2),r3(2))
    srefind(3) = sign(srefind(3),r3(3))
  endif
!**************************************
!  High frequency dielectric constant *
!**************************************
  if (lshello) then
    if (ldlch.and.nshell.gt.0) then
!
!  Collect shell second derivative terms
!
      msvar = 3*nshell
      if (nprocs.gt.1) then
        do i = 1,nshell
          ni = nshptr(i)
          qshell(i) = qf(ni)*occuf(ni)
        enddo
        do i = 1,nshellonnode
          is = nshonnodeptr(i)
          ni = node2atom(is)
          qshellnode(i) = qf(ni)*occuf(ni)
          ix = 3*(is-1) + 1
          iy = ix + 1
          iz = ix + 2
          ixs = 3*(i-1) + 1
          iys = ixs + 1
          izs = ixs + 2
          do j = 1,nshell
            nj = nshptr(j)
            jx = 3*(nj-1) + 1
            jy = jx + 1
            jz = jx + 2
            jxs = 3*(j-1) + 1
            jys = jxs + 1
            jzs = jxs + 2
            dervi(jxs,ixs) = derv2(jx,ix)
            dervi(jys,ixs) = derv2(jy,ix)
            dervi(jzs,ixs) = derv2(jz,ix)
            dervi(jxs,iys) = derv2(jx,iy)
            dervi(jys,iys) = derv2(jy,iy)
            dervi(jzs,iys) = derv2(jz,iy)
            dervi(jxs,izs) = derv2(jx,iz)
            dervi(jys,izs) = derv2(jy,iz)
            dervi(jzs,izs) = derv2(jz,iz)
          enddo
        enddo
      else
        do i = 1,nshell
          ni = nshptr(i)
          qshell(i) = qf(ni)*occuf(ni)
          ix = 3*(ni-1) + 1
          iy = ix + 1
          iz = ix + 2
          ixs = 3*(i-1) + 1
          iys = ixs + 1
          izs = ixs + 2
          do j = 1,i
            nj = nshptr(j)
            jx = 3*(nj-1) + 1
            jy = jx + 1
            jz = jx + 2
            jxs = 3*(j-1) + 1
            jys = jxs + 1
            jzs = jxs + 2
            dervi(jxs,ixs) = derv2(jx,ix)
            dervi(jys,ixs) = derv2(jy,ix)
            dervi(jzs,ixs) = derv2(jz,ix)
            dervi(jxs,iys) = derv2(jx,iy)
            dervi(jys,iys) = derv2(jy,iy)
            dervi(jzs,iys) = derv2(jz,iy)
            dervi(jxs,izs) = derv2(jx,iz)
            dervi(jys,izs) = derv2(jy,iz)
            dervi(jzs,izs) = derv2(jz,iz)
          enddo
        enddo
      endif
      if (nbss.gt.0) then
!
!  Collect radial second derivatives
!
        nbsc = nbs - nbss
        do i = 1,nbss
          iptr = n3f + nbsptr(nbsc+i)
          do j = 1,i
            jptr = n3f + nbsptr(nbsc+j)
            dervi(msvar+j,msvar+i) = derv2(jptr,iptr)
          enddo
          do j = 1,nshell
            nj = nshptr(j)
            jx = 3*(nj-1) + 1
            jy = jx + 1
            jz = jx + 2
            jxs = 3*(j-1) + 1
            jys = jxs + 1
            jzs = jxs + 2
            dervi(jxs,msvar+i) = derv2(jx,iptr)
            dervi(jys,msvar+i) = derv2(jy,iptr)
            dervi(jzs,msvar+i) = derv2(jz,iptr)
          enddo
        enddo
      endif
!
!  Compress derivatives for partial occupancies
!
      if (lpocc) then
        call compressd1(qshell,0_i4,nsfoc,nshell,iocshptr)
        call compressd2(dervi,maxd2,0_i4,nsfoc,nbsfoc,nshell,iocshptr,ibocshptr)
        ninvs = 3*nsfoc + nbsfoc
      else
        ninvs = msvar + nbss
      endif
!
!  Invert second derivative matrix
!
      ifail = 1
      call matrix_inversion_shells(ninvs,1_i4,maxd2,dervi,nshell,nshellonnode,ifail)
!
      if (ifail.gt.0) then
        nwarn = nwarn + 1
        call outwarning('Properties cannot be calculated - matrix is singular',0_i4)
        goto 999
      endif
!
!  Multiply inverse second derivatives by charge vectors
!
      if (nprocs.gt.1) then
        do i = 1,3
          do j = 1,3
            sum = 0.0_dp
            do k = 1,nshellonnode
              qak = qshellnode(k)
              indk = 3*(k-1) + i
              do l = 1,nsfoc
                indl = 3*(l-1) + j
                sum = sum + dervi(indl,indk)*qak*qshell(l)
              enddo
            enddo
            diconh(j,i) = dlcfactor*sum
          enddo
        enddo
!
!  Global sum
!
        call sumall(diconh,ditmp,3_i4*3_i4,"property3","ditmp")
        diconh(1:3,1:3) = ditmp(1:3,1:3)
!
!  Add one to diagonal terms
!
        do i = 1,3
          diconh(i,i) = diconh(i,i) + 1.0_dp
        enddo
      else
        do i = 1,3
          do j = 1,3
            sum = 0.0_dp
            do k = 1,nsfoc
              qak = qshell(k)
              indk = 3*(k-1) + i
              do l = 1,nsfoc
                indl = 3*(l-1) + j
                sum = sum + dervi(indl,indk)*qak*qshell(l)
              enddo
            enddo
            diconh(j,i) = dlcfactor*sum
          enddo
          diconh(i,i) = diconh(i,i) + 1.0_dp
        enddo
      endif
!
!  Calculate refractive indices
!
      do i = 1,3
        do j = 1,3
          r33(j,i) = diconh(j,i)
        enddo
      enddo
      call dsyev('N','U',3_i4,r33,3_i4,r3,w1l,9_i4,ifail)
      hfrefind(1) = sqrt(abs(r3(1)))
      hfrefind(2) = sqrt(abs(r3(2)))
      hfrefind(3) = sqrt(abs(r3(3)))
      hfrefind(1) = sign(hfrefind(1),r3(1))
      hfrefind(2) = sign(hfrefind(2),r3(2))
      hfrefind(3) = sign(hfrefind(3),r3(3))
    else
      do i = 1,3
        do j = 1,3
          diconh(j,i) = 0.0_dp
        enddo
        diconh(i,i) = 1.0_dp
        hfrefind(i) = 1.0_dp
      enddo
    endif
  endif
!*****************
!  Polarisation  *
!*****************
  if (lpolarisation) then
    polar_xyz(1:3) = 0.0_dp
!
!  Sum over charge times position
!
    do i = 1,numat
      polar_xyz(1) = polar_xyz(1) + qf(i)*xclat(i)
      polar_xyz(2) = polar_xyz(2) + qf(i)*yclat(i)
      polar_xyz(3) = polar_xyz(3) + qf(i)*zclat(i)
    enddo
!
!  Multiply by conversion factor for 4*pi/V and units to C/m^2
!
    polar_xyz(1:3) = piefactor*polar_xyz(1:3)
  endif
!**********************
!  Output properties  *
!**********************
  if (lprint.and.ioproc) then
    if (lelastic.and.abs(elcon(1,1)).gt.zero) then
      write(ioout,'(/)')
      if (index(keyword,'oldu').ne.0) then
        write(ioout,'(''  Elastic Constant Matrix: (Units=10**11 Dyne/cm**2= 10 GPa)'',/)')
      else
        write(ioout,'(''  Elastic Constant Matrix: (Units=GPa)'',/)')
      endif
      write(ioout,'(''-------------------------------------------------------------------------------'')')
      write(ioout,'(''  Indices      1         2         3         4         5         6    '')')
      write(ioout,'(''-------------------------------------------------------------------------------'')')
      if (index(keyword,'oldu').ne.0) then
        do i = 1,6
          write(ioout,'(7x,i1,2x,6f10.5)')i,(0.1_dp*elcon(j,i),j=1,6)
        enddo
      else
        do i = 1,6
          write(ioout,'(7x,i1,2x,6f10.4)')i,(elcon(j,i),j=1,6)
        enddo
      endif
      write(ioout,'(''-------------------------------------------------------------------------------'')')
      if (lceigen.and.lelconeigok) then
        write(ioout,'(/,''  Elastic Constant Tensor Eigenvalues and Eigenvectors: '',/)')
        write(ioout,'(''-------------------------------------------------------------------------------'')')
        write(ioout,'('' Eigenvalue      1          2          3          4          5          6    '')')
        write(ioout,'(''   (GPa)        (xx)       (yy)       (zz)       (yz)       (xz)       (xy)'')')
        write(ioout,'(''-------------------------------------------------------------------------------'')')
        do i = 1,6
          write(ioout,'(f10.5,1x,6f11.6)') elconeig(i),(elconvec(j,i),j=1,6)
        enddo
        write(ioout,'(''-------------------------------------------------------------------------------'')')
        if (nelconeig.ne.0) then
          write(ioout,'(''  No. of non-positive eigenvalues = '',i1,'' => system is elastically unstable'')') nelconeig
          write(ioout,'(''-------------------------------------------------------------------------------'')')
        endif
      endif
    endif
    if (lcompliances) then
      write(ioout,'(/)')
      write(ioout,'(''  Elastic Compliance Matrix: (Units=1/GPa)'',/)')
      write(ioout,'(''-------------------------------------------------------------------------------'')')
      write(ioout,'(''  Indices      1         2         3         4         5         6    '')')
      write(ioout,'(''-------------------------------------------------------------------------------'')')
      do i = 1,6
        write(ioout,'(7x,i1,2x,6f10.6)')i,(compliances(j,i),j=1,6)
      enddo
      write(ioout,'(''-------------------------------------------------------------------------------'')')
    endif
    if (lelastic) then
      write(ioout,'(/,''  Mechanical properties :'',/)')
      write(ioout,'(''-------------------------------------------------------------------------------'')')
      write(ioout,'(''  Convention : '',19x,''Reuss'',9x,''Voigt'',9x,''Hill'')')
      write(ioout,'(''-------------------------------------------------------------------------------'')')
      write(ioout,'(''  Bulk  Modulus (GPa)     = '',3(f13.5,1x))') &
        bulkmod_reuss,bulkmod_voigt,bulkmod_hill
      write(ioout,'(''  Shear Modulus (GPa)     = '',3(f13.5,1x))') &
        shearmod_reuss,shearmod_voigt,shearmod_hill
      write(ioout,'(''-------------------------------------------------------------------------------'')')
      write(ioout,'(''  Velocity S-wave (km/s)  = '',3(f13.5,1x))') vs_reuss,vs_voigt,vs_hill
      write(ioout,'(''  Velocity P-wave (km/s)  = '',3(f13.5,1x))') vp_reuss,vp_voigt,vp_hill
      write(ioout,'(''-------------------------------------------------------------------------------'')')
      if (abs(bulkmod_reuss).gt.1.0d-8) then
        write(ioout,'(''  Compressibility (1/GPa) = '',f13.8)') 1.0_dp/bulkmod_reuss
      endif
      write(ioout,'(''-------------------------------------------------------------------------------'')')
      write(ioout,'(''  Stress axis :'',21x,''x'',13x,''y'',13x,''z'')')
      write(ioout,'(''-------------------------------------------------------------------------------'')')
      write(ioout,'(''  Youngs Moduli (GPa)     = '',3(f13.5,1x))') (ym(i),i=1,3)
      write(ioout,'(''-------------------------------------------------------------------------------'')')
      write(ioout,'(''  Poissons Ratio (x)      = '',14x,2(f13.5,1x))')-compliances(2,1)*ym(2),-compliances(3,1)*ym(3)
      write(ioout,'(''  Poissons Ratio (y)      = '',f13.5,15x,f13.5)')-compliances(2,1)*ym(1),-compliances(3,2)*ym(3)
      write(ioout,'(''  Poissons Ratio (z)      = '',2(f13.5,1x))')-compliances(3,1)*ym(1),-compliances(3,2)*ym(2)
      write(ioout,'(''-------------------------------------------------------------------------------'')')
      write(ioout,'(/)')
      if (lcml) then
        call gulp_cml_output_derivs(elcon, compliances, bulkmod_reuss, bulkmod_voigt, bulkmod_hill,  &
             shearmod_reuss, shearmod_voigt, shearmod_hill, vs_reuss, vs_voigt, vs_hill,  &
             vp_reuss, vp_voigt, vp_hill)
      endif
    endif
    if (lpiezo) then
      write(ioout,'(''  Piezoelectric Strain Matrix: (Units=C/m**2)'',/)')
      write(ioout,'(''-------------------------------------------------------------------------------'')')
      write(ioout,'(''  Indices      1         2         3         4         5         6    '')')
      write(ioout,'(''-------------------------------------------------------------------------------'')')
      write(ioout,'(7x,''x'',2x,6f10.5)')(piezo(i,1),i=1,6)
      write(ioout,'(7x,''y'',2x,6f10.5)')(piezo(i,2),i=1,6)
      write(ioout,'(7x,''z'',2x,6f10.5)')(piezo(i,3),i=1,6)
      write(ioout,'(''-------------------------------------------------------------------------------'')')
      write(ioout,'(/)')
    endif
    if (lpiezs) then
      write(ioout,'(''  Piezoelectric Stress Matrix: (Units=10**-11 C/N)'',/)')
      write(ioout,'(''-------------------------------------------------------------------------------'')')
      write(ioout,'(''  Indices      1         2         3         4         5         6    '')')
      write(ioout,'(''-------------------------------------------------------------------------------'')')
      write(ioout,'(7x,''x'',2x,6f10.5)')(100.0_dp*piezs(i,1),i=1,6)
      write(ioout,'(7x,''y'',2x,6f10.5)')(100.0_dp*piezs(i,2),i=1,6)
      write(ioout,'(7x,''z'',2x,6f10.5)')(100.0_dp*piezs(i,3),i=1,6)
      write(ioout,'(''-------------------------------------------------------------------------------'')')
    endif
    if (ldlcs) then
      write(ioout,'(/)')
      write(ioout,'(''  Static dielectric constant tensor : '',/)')
      write(ioout,'(''-------------------------------------------------------------------------------'')')
      write(ioout,'(10x,4x,''x'',9x,''y'',9x,''z'')')
      write(ioout,'(''-------------------------------------------------------------------------------'')')
      write(ioout,'(7x,''x'',2x,3f10.5)')(dicons(1,i),i=1,3)
      write(ioout,'(7x,''y'',2x,3f10.5)')(dicons(2,i),i=1,3)
      write(ioout,'(7x,''z'',2x,3f10.5)')(dicons(3,i),i=1,3)
      write(ioout,'(''--------------------------------------------------------------------------------'',/)')
    endif
    if (ldlch.and.(nshell.gt.0)) then
      write(ioout,'(''  High frequency dielectric constant tensor : '',/)')
      write(ioout,'(''-------------------------------------------------------------------------------'')')
      write(ioout,'(10x,4x,''x'',9x,''y'',9x,''z'')')
      write(ioout,'(''-------------------------------------------------------------------------------'')')
      write(ioout,'(7x,''x'',2x,3f10.5)')(diconh(1,i),i=1,3)
      write(ioout,'(7x,''y'',2x,3f10.5)')(diconh(2,i),i=1,3)
      write(ioout,'(7x,''z'',2x,3f10.5)')(diconh(3,i),i=1,3)
      write(ioout,'(''-------------------------------------------------------------------------------'',/)')
    endif
    if (ldlcs) then
      write(ioout,'(''  Static refractive indices : '',/)')
      write(ioout,'(''-------------------------------------------------------------------------------'')')
      write(ioout,'(4x,''1 = '',f10.5,6x,''2 = '',f10.5,6x,''3 = '',f10.5)')(srefind(i),i=1,3)
      write(ioout,'(''-------------------------------------------------------------------------------'',/)')
    endif
    if (ldlch.and.(nshell.gt.0)) then
      write(ioout,'(''  High frequency refractive indices : '',/)')
      write(ioout,'(''-------------------------------------------------------------------------------'')')
      write(ioout,'(4x,''1 = '',f10.5,6x,''2 = '',f10.5,6x,''3 = '',f10.5)')(hfrefind(i),i=1,3)
      write(ioout,'(''-------------------------------------------------------------------------------'',/)')
    endif
    if (lpolarisation) then
      write(ioout,'(''  Polarisation: (Units=C/m**2)'',/)')
      write(ioout,'(''-------------------------------------------------------------------------------'')')
      write(ioout,'(4x,''x = '',f10.5,6x,''y = '',f10.5,6x,''z = '',f10.5)')(polar_xyz(i),i=1,3)
      write(ioout,'(''-------------------------------------------------------------------------------'',/)')
    endif
    if (lcml) then
      if (ldlch.and.(nshell.gt.0).and.(ldlcs).and.(lpiezo)) then
          call gulp_cml_output_dielectric(sdct=dicons, srind=srefind, &
                 & hfct=diconh, hrind=hfrefind, piezo=piezo, piezs=piezs)
      elseif (ldlch.and.(nshell.gt.0).and.(ldlcs)) then
          call gulp_cml_output_dielectric(sdct=dicons, srind=srefind, &
                 & hfct=diconh, hrind=hfrefind)
      elseif ( (ldlcs).and.(lpiezo) ) then
          call gulp_cml_output_dielectric(sdct=dicons, srind=srefind, piezo=piezo, piezs=piezs)
      elseif (ldlcs) then
          call gulp_cml_output_dielectric(sdct=dicons, srind=srefind)
      endif
    endif
  endif
!
!  Flush the output buffer
!
  call gflush(ioout)
!***************
!  Exit tasks  *
!***************
999 continue
!
!  Timings
!
  t2p = g_cpu_time()
  tprop = t2p - t1p + tprop
!
!  Deallocate memory
!
#ifdef MPI
  if (lrigid.and.nprocs.gt.1) then
    deallocate(d3nomols,stat=status)
    if (status/=0) call deallocate_error('property3','d3nomols')
    deallocate(d3nomol,stat=status)
    if (status/=0) call deallocate_error('property3','d3nomol')
  endif
#endif
  if (allocated(tmp2)) then
    deallocate(tmp2,stat=status)
    if (status/=0) call deallocate_error('property3','tmp2')
  endif
  if (lrigid) then
    deallocate(nmOKptr,stat=status)
    if (status/=0) call deallocate_error('property3','nmOKptr')
    deallocate(nmOK,stat=status)
    if (status/=0) call deallocate_error('property3','nmOK')
  endif
  if (lpocc) then
    deallocate(d3,stat=status)
    if (status/=0) call deallocate_error('property3','d3')
  endif
  deallocate(tmp,stat=status)
  if (status/=0) call deallocate_error('property3','tmp')
  if (nprocs.gt.1) then
    deallocate(qshellnode,stat=status)
    if (status/=0) call deallocate_error('property3','qshellnode')
  endif
  deallocate(qshell,stat=status)
  if (status/=0) call deallocate_error('property3','qshell')
  deallocate(qfo,stat=status)
  if (status/=0) call deallocate_error('property3','qfo')
#ifdef MPI
  if (lrigid.and.nprocs.gt.1) then
    deallocate(pdig2l,stat=status)
    if (status/=0) call deallocate_error('property3','pdig2l')
    deallocate(pdinode,stat=status)
    if (status/=0) call deallocate_error('property3','pdinode')
  endif
#endif
#ifdef TRACE
  call trace_out('property3')
#endif
!
  return
  end
