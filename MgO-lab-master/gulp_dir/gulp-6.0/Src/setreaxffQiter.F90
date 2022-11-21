  subroutine setreaxffQiter(nboatom,nboatomptr,nboatomRptr,nbos,nbosptr,qreaxFF,eself,loutputQ)
!
!  Subroutine for performing electronegativity equilisation calcns
!  in order to determine reaxFF charges. Iterative version.
!
!   3/08 Created based setreaxffQ
!   4/08 Spatial decomposition added
!   4/08 Qtot now properly initialised
!   4/08 qreaxFF returned to full atom array at end
!   6/08 Shell structure correction added
!   6/08 Option to fix charges in ReaxFF added
!   6/08 Shell structure correction added along with iterative solution
!   6/08 Option to fix charges in ReaxFF added
!   6/08 Dimension of problem being passed to bcgsolve corrected to nfree+1
!   7/08 Bugs in fixed charge handling fixed
!   7/08 Bug in spatial decomposition version fixed
!   7/08 Mixed fixed charge/variable charge case enabled for spatial
!  11/09 Modified for parallel execution of standard algorithm
!  11/09 Check on whether charges need to be reset corrected
!   4/10 Name of routine corrected in memory calls
!   8/10 Spatial decomposition algorithm corrected to use full D matrix
!   8/10 Parallelisation for spatial decomposition enabled with iterative charges.
!   1/11 Charge extrapolation added as an option
!   1/11 Tolerance increased for convergence in bcgsolve
!   1/11 Pointer to diagonal elements in listD added to speed up bcgsolve
!   2/11 Bug in setting of z matrix for spatial decomposition fixed - for mixed
!        case the contribution to z was being double counted as the algorithm is
!        not lower half triangular.
!   8/11 Modified to handle the case where a non-reaxff atom type is inserted into
!        a structure without an explicit frozen charge.
!   9/11 Overwriting of ni by nati fixed
!   9/11 Parallelisation of spatial decomposition corrected
!  10/11 Site energy terms added
!  11/11 Tolerance for convergence in bcgsolve made smaller by 10 to improve finite
!        difference second derivatives.
!  12/11 Dimension of nbos corrected to numat from nboatom
!  12/11 lreaxFFqfix and reaxFFgamma referenced by species
!  12/11 Shell structure correction changed to match that in setreaxffQ
!  11/13 qtot removed as an argument from bcgsolve as it wasn't being used.
!   2/14 Trap for iterative charge solve failing corrected
!   4/14 Fixed parameters for iterative charge solve replaced by input values
!   7/15 External potential added
!   8/15 Argument added to flag whether charges should be output
!   6/17 Handling of nfree changed to point to nboatoms
!   6/17 nboatomptr added as input argument
!   9/17 Correction to referencing of extpotcfg in parallel 
!   4/18 Parallel I/O corrected
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!   9/19 ReaxFF charge parameters now given by species
!   3/20 Location of angstoev changed to current
!   8/20 if statement for not cell buffer moved to avoid unnecessary calculation
!  12/20 Atom number format increased to allow for larger values
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
!  Julian Gale, CIC, Curtin University, December 2020
!
  use configurations, only : extpotcfg
  use control
  use current
  use element
  use energies,       only : siteenergy
  use iochannels
  use numbers,        only : third
  use parallel
  use reallocate
  use realvectors,    only : dist
  use reaxFFdata,     only : reaxFFcutoffQ, reaxFFqfix, reaxFFchi, reaxFFmu, reaxFFgamma
  use reaxFFdata,     only : reaxFFqdamp, nreaxFFqiter, reaxFFqconverged, reaxFFqconverged1
  use reaxFFdata,     only : reaxFFshell, lreaxFFqfix
  use spatialbo
  use times
  implicit none
#ifdef MPI
  include 'mpif.h'
#endif
!
!  Passed variables
!
  integer(i4), intent(in)                      :: nboatom
  integer(i4), intent(in)                      :: nboatomptr(numat)
  integer(i4), intent(in)                      :: nboatomRptr(nboatom)
  integer(i4), intent(in)                      :: nbos(numat)
  integer(i4), intent(in)                      :: nbosptr(numat)
  logical,     intent(in)                      :: loutputQ
  real(dp),    intent(out)                     :: qreaxFF(nboatom+1)
  real(dp),    intent(out)                     :: eself
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ic
  integer(i4)                                  :: ierror
  integer(i4)                                  :: ifree
  integer(i4)                                  :: ii
  integer(i4)                                  :: imx
  integer(i4)                                  :: imy
  integer(i4)                                  :: imz
  integer(i4)                                  :: ind1
  integer(i4)                                  :: ind2
  integer(i4)                                  :: indn
  integer(i4)                                  :: iptr
  integer(i4)                                  :: iter
  integer(i4)                                  :: ixyz
  integer(i4)                                  :: ix
  integer(i4)                                  :: iy
  integer(i4)                                  :: iz
  integer(i4), dimension(:,:), pointer, save   :: listD => null()
  integer(i4), dimension(:), allocatable       :: listDrowRptr
  integer(i4), dimension(:), allocatable       :: nDdiag
  integer(i4), dimension(:), allocatable       :: numD
  integer(i4), dimension(:), allocatable       :: numDold
  integer(i4)                                  :: j
  integer(i4)                                  :: jc
  integer(i4)                                  :: jfree
  integer(i4)                                  :: jj
  integer(i4),                          save   :: maxD = 1
  integer(i4)                                  :: maxx
  integer(i4)                                  :: maxxy
  integer(i4)                                  :: n
  integer(i4)                                  :: n1i
  integer(i4)                                  :: n1j
  integer(i4)                                  :: nn
  integer(i4)                                  :: nD
  integer(i4)                                  :: nfree
  integer(i4)                                  :: nfreeloc
  integer(i4), dimension(:), allocatable       :: nfreeptr
  integer(i4), dimension(:), allocatable       :: nfreelocptr
  integer(i4), dimension(:), allocatable       :: nfreelocnptr
  integer(i4), dimension(:), allocatable       :: nfreelocRptr
  integer(i4), dimension(:), allocatable       :: nfreeRptr
  integer(i4)                                  :: nfi
  integer(i4)                                  :: nfii
  integer(i4)                                  :: nfj
  integer(i4)                                  :: ni
  integer(i4)                                  :: nj
  integer(i4)                                  :: niter
  integer(i4)                                  :: nitermax
  integer(i4)                                  :: nmolonly
  integer(i4)                                  :: nor
  integer(i4)                                  :: nspeci
  integer(i4)                                  :: nspecj
  integer(i4)                                  :: nsplower(3)
  integer(i4)                                  :: nspupper(3)
  integer(i4), dimension(:), allocatable       :: nwrk
  integer(i4)                                  :: numDcurrent
  integer(i4)                                  :: status
#ifdef MPI
  integer(i4)                                  :: numDtmp
  integer                                      :: MPIerror
  integer                                      :: ntag
  integer                                      :: nnode
  integer                                      :: ntmp
  integer                                      :: Request
  integer(i4), dimension(:), allocatable       :: StatMPI       ! Array for status from MPI
#endif
  logical                                      :: lfixQi
  logical                                      :: lfixQj
  logical                                      :: literative
  logical                                      :: lmixed
  logical                                      :: lself
  real(dp),    dimension(:,:),   pointer, save :: D => null()
  real(dp),    dimension(:,:), allocatable     :: Dsave
#ifdef MPI
  real(dp),    dimension(:),   allocatable     :: dtmp
#endif
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: cut2
  real(dp)                                     :: damp
  real(dp)                                     :: err
  real(dp)                                     :: esite
  real(dp)                                     :: gam
  real(dp)                                     :: gammai
  real(dp)                                     :: gammaij
  real(dp)                                     :: gammaj
  real(dp)                                     :: qdiff
  real(dp)                                     :: qdiff1
  real(dp)                                     :: qfree
  real(dp),    dimension(:),   allocatable     :: qreaxFFfree
  real(dp),    dimension(:),   allocatable     :: qreaxFFsave
  real(dp)                                     :: qi
  real(dp)                                     :: qsum
  real(dp)                                     :: qtot
  real(dp)                                     :: rij
  real(dp)                                     :: rij2
  real(dp)                                     :: tp
  real(dp)                                     :: tsuml
  real(dp)                                     :: dtpdr
  real(dp)                                     :: d2tpdr2
  real(dp)                                     :: xi
  real(dp)                                     :: yi
  real(dp)                                     :: zi
  real(dp)                                     :: xji
  real(dp)                                     :: yji
  real(dp)                                     :: zji
  real(dp),    dimension(:),   allocatable     :: z
  real(dp),    dimension(:),   allocatable     :: ztmp
! 
!  Allocate pointer array for atoms to act on  
! 
  allocate(nfreeptr(nboatom),stat=status)      
  if (status/=0) call outofmemory('setreaxffQiter','nfreeptr')
  allocate(nfreeRptr(nboatom),stat=status)      
  if (status/=0) call outofmemory('setreaxffQiter','nfreeRptr')
  allocate(nfreelocptr(nboatom),stat=status)      
  if (status/=0) call outofmemory('setreaxffQiter','nfreelocptr')
  allocate(nfreelocnptr(nboatom),stat=status)      
  if (status/=0) call outofmemory('setreaxffQiter','nfreelocnptr')
  allocate(nfreelocRptr(nboatom),stat=status)      
  if (status/=0) call outofmemory('setreaxffQiter','nfreelocRptr')
! 
!  Find out whether an iterative solution is needed due to shell structure or whether there are any fixed charges
! 
  literative = .false.
  qsum = 0.0_dp
  qfree = 0.0_dp
  nfree = 0
  nfreeRptr(1:nboatom) = 0
  do i = 1,nboatom
    ii = nboatomRptr(i)
    if (ii.gt.0) then
      nspeci = nbosptr(ii)
      if (abs(reaxFFshell(1,nspeci)).gt.1.0d-12) then
        literative = .true.
      endif
      if (lreaxFFqfix(nspeci)) then
        qsum = qsum - reaxFFqfix(nspeci)*occuf(ii)
        qreaxFF(i) = reaxFFqfix(nspeci)
      elseif (nbos(ii).gt.0) then
        nfree = nfree + 1
        nfreeptr(nfree) = i
        nfreeRptr(i) = nfree
        qfree = qfree + qreaxFF(i)
      endif
    endif
  enddo
  if (nprocs.gt.1) then
    nfreelocnptr(1:nboatom) = 0_i4
    nfreelocRptr(1:nboatom) = 0_i4
    if (lspatialboOK) then
!
!  Loop over all local spatial cells except buffer regions
!
      nfreeloc = 0
      do ixyz = 1,ncellpernodebo
        if (.not.lbuffercellbo(ixyz)) then
          ind1 = ncellnodeptrbo(ixyz)
!
!  Get number of atoms in this cell
!
          ni = nspcellatbo(ind1)
          n1i = nspcellat1ptrbo(ind1)
!
!  Outer loop over atoms within this cell
!
          do ii = 1,ni
            i = nspcellatptrbo(n1i+ii)
            nspeci = nbosptr(i)
            ifree = nfreeRptr(i)
            lfixQi = lreaxFFqfix(nspeci)
!
!  Does i have a reaxFF species?
!
            if (nbos(i).gt.0.and..not.lfixQi) then
!
!  End check on whether i has a reaxFF species
!
              nfreeloc = nfreeloc + 1
              nfreelocptr(nfreeloc) = ifree
              nfreelocnptr(ifree) = procid
              nfreelocRptr(ifree) = nfreeloc
            endif
!
!  End loop over atom i
!
          enddo
!
!  End checks on whether cell is required
!
        endif
!
!  End loop over cells on node
!
      enddo
    else
!
!  For parallel case find number of local free atoms
!
      nfreeloc = 0
      do i = procid+1,nboatom,nprocs
        ii = nboatomRptr(i)
        if (ii.gt.0) then
          nspeci = nbosptr(ii)
          if (.not.lreaxFFqfix(nspeci)) then
            nfreeloc = nfreeloc + 1
            nfreelocptr(nfreeloc) = nfreeRptr(i)
            nfreelocnptr(nfreeRptr(i)) = procid
            nfreelocRptr(nfreeRptr(i)) = nfreeloc
          endif
        endif
      enddo
    endif
!
!  Globalise nfreelocnptr
!
    allocate(nwrk(nboatom),stat=status)      
    if (status/=0) call outofmemory('setreaxffQiter','nwrk')
!
    call isumall(nfreelocnptr,nwrk,nboatom,"setreaxffQiter","nfreelocnptr")
    nfreelocnptr(1:nboatom) = nwrk(1:nboatom)
!
    deallocate(nwrk,stat=status)
    if (status/=0) call deallocate_error('setreaxffQiter','nwrk')
!
  else
    nfreeloc = nfree
    do i = 1,nfree
      nfreelocptr(i) = i
      nfreelocnptr(i) = 0_i4
      nfreelocRptr(i) = i
    enddo
  endif
!
!  Decide on total EEM fragment charge - for cluster
!  there is no constraint so use sum of initial charges.
!  For periodic system, charge must be equal to the
!  negative sum of the non-EEM ion charges.
!
  qtot = totalcharge + qsum
!
!  If initial charges don't sum to the right amount then correct
!
  if (abs(qfree-qtot).gt.1.0d-5) then
    qfree = qtot/dble(nfree)
    do i = 1,nfree
      ii = nfreeptr(i)
      qreaxFF(ii) = qfree
    enddo
  endif
!
!  Allocate local memory that depends on nboatom
!
  n = nfree + 1
  allocate(nDdiag(nfreeloc),stat=status)
  if (status/=0) call outofmemory('setreaxffQiter','nDdiag')
  allocate(numD(nfreeloc),stat=status)
  if (status/=0) call outofmemory('setreaxffQiter','numD')
  allocate(z(n),stat=status)
  if (status/=0) call outofmemory('setreaxffQiter','z')
  allocate(qreaxFFfree(nfree+1),stat=status)
  if (status/=0) call outofmemory('setreaxffQiter','qreaxFFfree')
  allocate(qreaxFFsave(nboatom+1),stat=status)
  if (status/=0) call outofmemory('setreaxffQiter','qreaxFFsave')
!
  maxD = max(maxD,min(20,nfreeloc))
  call realloc(listD,maxD,nfreeloc,ierror)
  if (ierror.ne.0) call outofmemory('setreaxffQiter','listD')
  call realloc(D,maxD,nfreeloc,ierror)
  if (ierror.ne.0) call outofmemory('setreaxffQiter','D')
!  
!  Form right hand vector
! 
  z(1:n) = 0.0_dp
  if (nprocs.gt.1) then
    if (ioproc) then
      do i = 1,nfree
        ii = nboatomRptr(nfreeptr(i))
        nspeci = nbosptr(ii)
        z(i) = z(i) - reaxFFchi(nspeci) - extpotcfg(nsft+nrelf2a(ii))
      enddo
      z(n) = qtot
    else
      z(n) = qtot
    endif
  else
    do i = 1,nfree
      ii = nboatomRptr(nfreeptr(i))
      nspeci = nbosptr(ii)
      z(i) = z(i) - reaxFFchi(nspeci) - extpotcfg(nsft+nrelf2a(ii))
    enddo
    z(n) = qtot
  endif
!*******************************************
!  Generate electrostatic potential terms  *
!*******************************************
!
!  Initialise sparsity related counters
!
  numD(1:nfreeloc) = 0
!
!  Set cutoffs
!
  cut2 = reaxFFcutoffQ**2
!
!  Loop over first atom
!
  if (lspatialboOK) then
!************************************
!  Spatial decomposition algorithm  *
!************************************
    maxxy = nspcellbo(1)*nspcellbo(2)
    maxx  = nspcellbo(1)
!
    allocate(listDrowRptr(nfree),stat=status)
    if (status/=0) call outofmemory('setreaxffQiter','listDrowRptr')
!
    listDrowRptr(1:nfree) = 0
!
!  Set diagonal to be first element in D
!
    do i = 1,nfreeloc
      numD(i) = 1
      nDdiag(i) = 1
      listD(1,i) = nfreelocptr(i)
      D(1,i) = 0.0_dp
    enddo
    if (nprocs.gt.1) then
!
!  Loop over all local spatial cells except buffer regions
!
      nfi = 0
      do ixyz = 1,ncellpernodebo
        if (.not.lbuffercellbo(ixyz)) then
          ind1 = ncellnodeptrbo(ixyz)
          ind2 = ind1 - 1
          iz = ind2/maxxy
          ind2 = ind2 - maxxy*iz
          iy = ind2/maxx
          ix = ind2 - maxx*iy + 1
          iy = iy + 1
          iz = iz + 1
!
!  Set cell search bounds
!  
          nspupper(1) = min(ix+ncellsearchbo(1),nspcellbo(1))
          nspupper(2) = min(iy+ncellsearchbo(2),nspcellbo(2))
          nspupper(3) = min(iz+ncellsearchbo(3),nspcellbo(3))
          nsplower(1) = max(ix-ncellsearchbo(1),1)
          nsplower(2) = max(iy-ncellsearchbo(2),1)
          nsplower(3) = max(iz-ncellsearchbo(3),1)
!
!  Get number of atoms in this cell
!
          ni = nspcellatbo(ind1)
          n1i = nspcellat1ptrbo(ind1)
!
!  Outer loop over atoms within this cell
!
          do ii = 1,ni
            i = nspcellatptrbo(n1i+ii)
            nspeci = nbosptr(i)
            ifree = nfreeRptr(nboatomptr(i))
            lfixQi = lreaxFFqfix(nspeci)
!
!  Does i have a reaxFF species?
!
            if (nbos(i).gt.0.or.lfixQi) then
              gammai = reaxFFgamma(nspeci)
              ic = nspcellatptrcellbo(n1i+ii)
!
!  Set coordinates of atom i
!
              xi = xinboxbo(i) + xvec2cell(ic)
              yi = yinboxbo(i) + yvec2cell(ic)
              zi = zinboxbo(i) + zvec2cell(ic)
!
!  Set up row reverse pointer for i
!
              if (.not.lfixQi) then
                nfi = nfi + 1
                do iptr = 1,numD(nfi)
                  listDrowRptr(listD(iptr,nfi)) = iptr
                enddo
              endif
!
!  Loop over neighbouring cells
!         
              do imz = nsplower(3),nspupper(3)
                do imy = nsplower(2),nspupper(2)
                  do imx = nsplower(1),nspupper(1)
                    indn = (imz-1)*maxxy + (imy-1)*maxx + imx
!  
!  Loop over atoms within neighbouring cells
!               
                    nj = nspcellatbo(indn)
                    n1j = nspcellat1ptrbo(indn)
                    jjloopp: do jj = 1,nj
                      j = nspcellatptrbo(n1j+jj)
!
                      nspecj = nbosptr(j)
                      jfree = nfreeRptr(nboatomptr(j))
                      lfixQj = lreaxFFqfix(nspecj)
!
!  Does j have a reaxFF species?
!
                      if (nbos(j).gt.0.or.lfixQj) then
!
!  Skip loop for pairs of fixed atoms
!
                        if (lfixQi.and.lfixQj) then
                          cycle jjloopp
                        elseif (.not.lfixQi.and..not.lfixQj) then
!
!  Only increment matrix element for pairs of free atoms
!
                          lmixed = .false.
                        else
                          lmixed = .true.
                        endif
!
                        gammaj = reaxFFgamma(nspecj)
                        jc = nspcellatptrcellbo(n1j+jj)
!  
!  Exclude self-interaction
!                 
                        if (i.ne.j.or.(ind1.ne.indn)) then
!  
!  Set coordinate differences and calculate square of distance
!                   
                          xji = xvec2cell(jc) + xinboxbo(j) - xi
                          yji = yvec2cell(jc) + yinboxbo(j) - yi
                          zji = zvec2cell(jc) + zinboxbo(j) - zi
                          rij2 = xji*xji + yji*yji + zji*zji
                          if (rij2.lt.cut2) then
                            rij = sqrt(rij2)
!  
!  Compute Coulomb shielding parameters
!             
                            gammaij  = sqrt(gammai*gammaj)
                            gammaij  = 1.0_dp/(gammaij**3)
!
!  Compute taper function(s)
!
                            call p7reaxFFqtaper(rij,reaxFFcutoffQ,tp,dtpdr,d2tpdr2,.false.,.false.)
!
!  Pure real space tapered lreaxFF case
!
                            gam = tp/(rij**3 + gammaij)**third
                            if (lmixed) then
!
!  For mixed case add to Z array
!
                              if (lfixQi) then
                                z(jfree) = z(jfree) - gam*reaxFFqfix(nspeci)*angstoev
                              endif
                            else
!
!  Check whether this element exists in sparse array and if not add
!
                              if (listDrowRptr(jfree).eq.0) then
                                numD(nfi) = numD(nfi) + 1
                                if (numD(nfi).gt.maxD) then
                                  maxD = numD(nfi) + 10
                                  call realloc(listD,maxD,nfreeloc,ierror)
                                  if (ierror.ne.0) call outofmemory('setreaxffQiter','listD')
                                  call realloc(D,maxD,nfreeloc,ierror)
                                  if (ierror.ne.0) call outofmemory('setreaxffQiter','D')
                                endif
                                listD(numD(nfi),nfi) = jfree
                                listDrowRptr(jfree) = numD(nfi)
                                D(numD(nfi),nfi) = 0.0_dp
                              endif
!
!  Add to sparse array
!
                              D(listDrowRptr(jfree),nfi) = D(listDrowRptr(jfree),nfi) + gam
                            endif
                          endif
                        endif
                      endif
!
!  End loop over j
!
                    enddo jjloopp
!               
!  End loops over neighbouring cells
!  
                  enddo
                enddo
              enddo
!
!  Reinitialise row reverse pointer for i
!
              if (.not.lfixQi) then
                do iptr = 1,numD(nfi)
                  listDrowRptr(listD(iptr,nfi)) = 0
                enddo
              endif
!
!  End check on whether i has a reaxFF species
!
            endif
!
!  End loop over atom i
!
          enddo
!
!  End checks on whether cell is required
!
        endif
!
!  End loop over cells on node
!
      enddo
! 
!  Globalization of z
! 
      tsuml = g_cpu_time()
      allocate(ztmp(n))
      call sumall(z,ztmp,nfree,"setreaxffQiter","z") 
      z(1:nfree) = ztmp(1:nfree)
      deallocate(ztmp)
      tsum = tsum + g_cpu_time() - tsuml
    else
!
!  Loop over all local spatial cells except buffer regions
!
      do ixyz = 1,ncellpernodebo
        if (.not.lbuffercellbo(ixyz)) then
          ind1 = ncellnodeptrbo(ixyz)
          ind2 = ind1 - 1
          iz = ind2/maxxy
          ind2 = ind2 - maxxy*iz
          iy = ind2/maxx
          ix = ind2 - maxx*iy + 1
          iy = iy + 1
          iz = iz + 1
!
!  Set cell search bounds
!  
          nspupper(1) = min(ix+ncellsearchbo(1),nspcellbo(1))
          nspupper(2) = min(iy+ncellsearchbo(2),nspcellbo(2))
          nspupper(3) = min(iz+ncellsearchbo(3),nspcellbo(3))
          nsplower(1) = max(ix-ncellsearchbo(1),1)
          nsplower(2) = max(iy-ncellsearchbo(2),1)
          nsplower(3) = max(iz-ncellsearchbo(3),1)
!
!  Get number of atoms in this cell
!
          ni = nspcellatbo(ind1)
          n1i = nspcellat1ptrbo(ind1)
!
!  Outer loop over atoms within this cell
!
          do ii = 1,ni
            i = nspcellatptrbo(n1i+ii)
            nspeci = nbosptr(i)
            ifree = nfreeRptr(nboatomptr(i))
            lfixQi = lreaxFFqfix(nspeci)
!
!  Does i have a reaxFF species?
!
            if (nbos(i).gt.0.or.lfixQi) then
              gammai = reaxFFgamma(nspeci)
              ic = nspcellatptrcellbo(n1i+ii)
!
!  Set coordinates of atom i
!
              xi = xinboxbo(i) + xvec2cell(ic)
              yi = yinboxbo(i) + yvec2cell(ic)
              zi = zinboxbo(i) + zvec2cell(ic)
!
!  Set up row reverse pointer for i
!
              if (.not.lfixQi) then
                do iptr = 1,numD(ifree)
                  listDrowRptr(listD(iptr,ifree)) = iptr
                enddo
              endif
!
!  Loop over neighbouring cells
!         
              do imz = nsplower(3),nspupper(3)
                do imy = nsplower(2),nspupper(2)
                  do imx = nsplower(1),nspupper(1)
                    indn = (imz-1)*maxxy + (imy-1)*maxx + imx
!  
!  Loop over atoms within neighbouring cells
!               
                    nj = nspcellatbo(indn)
                    n1j = nspcellat1ptrbo(indn)
                    jjloop: do jj = 1,nj
                      j = nspcellatptrbo(n1j+jj)
!
                      nspecj = nbosptr(j)
                      jfree = nfreeRptr(nboatomptr(j))
                      lfixQj = lreaxFFqfix(nspecj)
!
!  Does j have a reaxFF species?
!
                      if (nbos(j).gt.0.or.lfixQj) then
!
!  Skip loop for pairs of fixed atoms
!
                        if (lfixQi.and.lfixQj) then
                          cycle jjloop
                        elseif (.not.lfixQi.and..not.lfixQj) then
!
!  Only increment matrix element for pairs of free atoms
!
                          lmixed = .false.
                        else
                          lmixed = .true.
                        endif
!
                        gammaj = reaxFFgamma(nspecj)
                        jc = nspcellatptrcellbo(n1j+jj)
!  
!  Exclude self-interaction
!                 
                        if (i.ne.j.or.(ind1.ne.indn)) then
!  
!  Set coordinate differences and calculate square of distance
!                   
                          xji = xvec2cell(jc) + xinboxbo(j) - xi
                          yji = yvec2cell(jc) + yinboxbo(j) - yi
                          zji = zvec2cell(jc) + zinboxbo(j) - zi
                          rij2 = xji*xji + yji*yji + zji*zji
                          if (rij2.lt.cut2) then
                            rij = sqrt(rij2)
!  
!  Compute Coulomb shielding parameters
!             
                            gammaij  = sqrt(gammai*gammaj)
                            gammaij  = 1.0_dp/(gammaij**3)
!
!  Compute taper function(s)
!
                            call p7reaxFFqtaper(rij,reaxFFcutoffQ,tp,dtpdr,d2tpdr2,.false.,.false.)
!
!  Pure real space tapered lreaxFF case
!
                            gam = tp/(rij**3 + gammaij)**third
                            if (lmixed) then
!
!  For mixed case add to Z array
!
                              if (lfixQi) then
                                z(jfree) = z(jfree) - gam*reaxFFqfix(nspeci)*angstoev
                              endif
                            else
!
!  Check whether this element exists in sparse array and if not add
!
                              if (listDrowRptr(jfree).eq.0) then
                                numD(ifree) = numD(ifree) + 1
                                if (ifree.eq.jfree) then
                                  nDdiag(ifree) = numD(ifree)
                                endif
                                if (numD(ifree).gt.maxD) then
                                  maxD = numD(ifree) + 10
                                  call realloc(listD,maxD,nfree,ierror)
                                  if (ierror.ne.0) call outofmemory('setreaxffQiter','listD')
                                  call realloc(D,maxD,nfree,ierror)
                                  if (ierror.ne.0) call outofmemory('setreaxffQiter','D')
                                endif
                                listD(numD(ifree),ifree) = jfree
                                listDrowRptr(jfree) = numD(ifree)
                                D(numD(ifree),ifree) = 0.0_dp
                              endif
!
!  Add to sparse array
!
                              D(listDrowRptr(jfree),ifree) = D(listDrowRptr(jfree),ifree) + gam
                            endif
                          endif
                        endif
                      endif
!
!  End loop over j
!
                    enddo jjloop
!               
!  End loops over neighbouring cells
!  
                  enddo
                enddo
              enddo
!
!  Reinitialise row reverse pointer for i
!
              if (.not.lfixQi) then
                do iptr = 1,numD(ifree)
                  listDrowRptr(listD(iptr,ifree)) = 0
                enddo
              endif
!
!  End check on whether i has a reaxFF species
!
            endif
!
!  End loop over atom i
!
          enddo
!
!  End checks on whether cell is required
!
        endif
!
!  End loop over cells on node
!
      enddo
    endif
!     
!  Scale matrix elements by conversion factor
!       
    do i = 1,nfreeloc
      do j = 1,numD(i)
        D(j,i) = D(j,i)*angstoev
      enddo
    enddo
!
!  Deallocate spatial specific arrays
!
    deallocate(listDrowRptr,stat=status)
    if (status/=0) call deallocate_error('setreaxffQiter','listDrowRptr')
  else
!***************************
!  Standard N^2 algorithm  *
!***************************
    if (nprocs.gt.1) then
!
!  Parallel version
!
      nfi = 0
      do ii = procid+1,nboatom,nprocs
        i = nboatomRptr(ii)
        if (i.gt.0) then
          nspeci = nbosptr(i)
          lfixQi = lreaxFFqfix(nspeci)
          if (.not.lfixQi) then
            nfi = nfi + 1
            nfii = nfreelocptr(nfi)
          endif
!
!  Does i have a reaxFF species?
!
          if (nbos(i).gt.0.or.lfixQi) then
            gammai = reaxFFgamma(nspeci)
            if (.not.lfixQi) then
!
!  Add self-term to sparse matrix structure as first matrix element in row
!
              numD(nfi) = numD(nfi) + 1
              nDdiag(nfi) = 1
              if (numD(nfi).gt.maxD) then
                maxD = numD(nfi) + 10
                call realloc(listD,maxD,nfreeloc,ierror)
                if (ierror.ne.0) call outofmemory('setreaxffQiter','listD')
                call realloc(D,maxD,nfreeloc,ierror)
                if (ierror.ne.0) call outofmemory('setreaxffQiter','D')
              endif
              listD(numD(nfi),nfi) = nfii
              D(numD(nfi),nfi) = 0.0_dp
            endif
!
!  If i is fixed then skip j atoms
!
            if (.not.lfixQi) then
!
!  Loop over second atom
!
            nfj = 0
            pjloop: do jj = 1,nboatom
              j = nboatomRptr(jj)
              if (j.gt.0) then
                nspecj = nbosptr(j)
                lfixQj = lreaxFFqfix(nspecj)
                if (.not.lfixQj) nfj = nfj + 1
!
!  Does j have a reaxFF species?
!
                if (nbos(j).gt.0.or.lfixQj) then
!
!  Skip loop for pairs of fixed atoms
!
                  if (lfixQi.and.lfixQj) then
                    cycle pjloop
                  elseif (.not.lfixQi.and..not.lfixQj) then
!
!  Only increment matrix element for pairs of free atoms
!               
                    lmixed = .false.
                  else
                    lmixed = .true.
                  endif
!             
                  gammaj = reaxFFgamma(nspecj)
!
!  Compute basic interatomic vector
!
                  xji = xclat(j) - xclat(i)
                  yji = yclat(j) - yclat(i)
                  zji = zclat(j) - zclat(i)
!
!  Find valid vectors
!
                  nor = 0
                  if (ndim.eq.3) then
                    call rsearch3D(xji,yji,zji,.false.,.false.,i,j,0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2)
                  elseif (ndim.eq.2) then
                    call rsearch2D(xji,yji,zji,.false.,.false.,i,j,0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2)
                  elseif (ndim.eq.1) then
                    call rsearch1D(xji,yji,zji,.false.,.false.,i,j,0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2)
                  elseif (ndim.eq.0) then
                    rij2 = xji*xji + yji*yji + zji*zji
                    if (i.ne.j) then
                      if (rij2.lt.cut2) then
                        nor = nor + 1
                      endif
                    else
                      nor = 0
                    endif
                  endif
                  if (nor.gt.0.and.ii.ne.jj.and..not.lmixed) then
!
!  Add to sparse matrix structure
!
                    numD(nfi) = numD(nfi) + 1
                    if (numD(nfi).gt.maxD) then
                      maxD = numD(nfi) + 10
                      call realloc(listD,maxD,nfreeloc,ierror)
                      if (ierror.ne.0) call outofmemory('setreaxffQiter','listD')
                      call realloc(D,maxD,nfreeloc,ierror)
                      if (ierror.ne.0) call outofmemory('setreaxffQiter','D')
                    endif
                    listD(numD(nfi),nfi) = nfj
                    D(numD(nfi),nfi) = 0.0_dp
                  endif
!
!  Set pointer to element of sparse D where matrix elements should be added
!
                  if (ii.eq.jj) then
                    numDcurrent = 1
                  else
                    numDcurrent = numD(nfi)
                  endif
!
!  If no distances then cycle
!
                  if (nor.eq.0) cycle pjloop
!  
!  Compute Coulomb shielding parameters
!             
                  gammaij  = sqrt(gammai*gammaj)
                  gammaij  = 1.0_dp/(gammaij**3)
!
!  Loop over valid distances and calculate contributions
!
                  if (ndim.gt.0) then
                    do nn = 1,nor
                      rij2 = dist(nn)
                      rij = sqrt(rij2)
!
!  Compute taper function
!
                      call p7reaxFFqtaper(rij,reaxFFcutoffQ,tp,dtpdr,d2tpdr2,.false.,.false.)
!                     
!  Pure real space tapered lreaxFF case
!                     
                      gam = tp/(rij2*rij + gammaij)**third
                      if (lmixed) then
                        if (.not.lfixQi) then
                          z(nfii) = z(nfii) - gam*reaxFFqfix(nspecj)*angstoev
                        endif
                      else
                        D(numDcurrent,nfi) = D(numDcurrent,nfi) + gam
                      endif
                    enddo
                  else
                    rij = sqrt(rij2)
!
!  Compute taper function
!
                    call p7reaxFFqtaper(rij,reaxFFcutoffQ,tp,dtpdr,d2tpdr2,.false.,.false.)
!                     
!  Pure real space tapered lreaxFF case
!                     
                    gam = tp/(rij2*rij + gammaij)**third
                    if (lmixed) then
                      if (.not.lfixQi) then
                        z(nfii) = z(nfii) - gam*reaxFFqfix(nspecj)*angstoev
                      endif
                    else
                      D(numDcurrent,nfi) = D(numDcurrent,nfi) + gam
                    endif
                  endif
                endif
              endif
!
!  End of loop over j
!
            enddo pjloop
!
!  End of fixQi exclusion
!
            endif
          endif
        endif
        if (.not.lfixQi) then
!
!  Scale matrix elements by conversion factor 
!
          do j = 1,numD(nfi)
            D(j,nfi) = D(j,nfi)*angstoev
          enddo
        endif
!
!  End of loop over i
!
      enddo
! 
!  Globalization of z
! 
      tsuml = g_cpu_time()
      allocate(ztmp(n))
      call sumall(z,ztmp,nfree,"setreaxffQiter","z") 
      z(1:nfree) = ztmp(1:nfree)
      deallocate(ztmp)
      tsum = tsum + g_cpu_time() - tsuml
    else
!
!  Serial version
!
      nfi = 0
      do ii = 1,nboatom
        i = nboatomRptr(ii)
        if (i.gt.0) then
          nspeci = nbosptr(i)
          lfixQi = lreaxFFqfix(nspeci)
          if (.not.lfixQi.and.nbos(i).gt.0) nfi = nfi + 1
!
!  Does i have a reaxFF species?
!
          if (nbos(i).gt.0.or.lfixQi) then
            gammai = reaxFFgamma(nspeci)
!
!  Loop over second atom
!
            if (lfixQi) then
              nfj = nfi
            else
              nfj = nfi - 1
            endif
            jloop: do jj = ii,nboatom
              j = nboatomRptr(jj)
              if (j.gt.0) then
                nspecj = nbosptr(j)
                lfixQj = lreaxFFqfix(nspecj)
                if (.not.lfixQj.and.nbos(j).gt.0) nfj = nfj + 1
!
!  Does j have a reaxFF species?
!
                if (nbos(j).gt.0.or.lfixQj) then
!
!  Skip loop for pairs of fixed atoms
!
                  if (lfixQi.and.lfixQj) then
                    cycle jloop
                  elseif (.not.lfixQi.and..not.lfixQj) then
!
!  Only increment matrix element for pairs of free atoms
!               
                    lmixed = .false.
                  else
                    lmixed = .true.
                  endif
!             
                  gammaj = reaxFFgamma(nspecj)
!
!  Compute basic interatomic vector
!
                  xji = xclat(j) - xclat(i)
                  yji = yclat(j) - yclat(i)
                  zji = zclat(j) - zclat(i)
!
!  Find valid vectors
!
                  nor = 0
                  if (ndim.eq.3) then
                    call rsearch3D(xji,yji,zji,.false.,.false.,i,j,0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2)
                  elseif (ndim.eq.2) then
                    call rsearch2D(xji,yji,zji,.false.,.false.,i,j,0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2)
                  elseif (ndim.eq.1) then
                    call rsearch1D(xji,yji,zji,.false.,.false.,i,j,0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2)
                  elseif (ndim.eq.0) then
                    rij2 = xji*xji + yji*yji + zji*zji
                    if (i.ne.j) then
                      if (rij2.lt.cut2) then
                        nor = nor + 1
                      endif
                    else
                      nor = 0
                    endif
                  endif
                  if ((nor.gt.0.or.ii.eq.jj).and..not.lmixed) then
!
!  Add to sparse matrix structure
!
                    numD(nfi) = numD(nfi) + 1
                    if (nfi.eq.nfj) then
                      nDdiag(nfi) = numD(nfi)
                    endif
                    if (numD(nfi).gt.maxD) then
                      maxD = numD(nfi) + 10
                      call realloc(listD,maxD,nfree,ierror)
                      if (ierror.ne.0) call outofmemory('setreaxffQiter','listD')
                      call realloc(D,maxD,nfree,ierror)
                      if (ierror.ne.0) call outofmemory('setreaxffQiter','D')
                    endif
                    listD(numD(nfi),nfi) = nfj
                    D(numD(nfi),nfi) = 0.0_dp
                  endif
!
!  If no distances then cycle
!
                  if (nor.eq.0) cycle jloop
!  
!  Compute Coulomb shielding parameters
!             
                  gammaij  = sqrt(gammai*gammaj)
                  gammaij  = 1.0_dp/(gammaij**3)
!
!  Loop over valid distances and calculate contributions
!
                  if (ndim.gt.0) then
                    do nn = 1,nor
                      rij2 = dist(nn)
                      rij = sqrt(rij2)
!
!  Compute taper function
!
                      call p7reaxFFqtaper(rij,reaxFFcutoffQ,tp,dtpdr,d2tpdr2,.false.,.false.)
!                     
!  Pure real space tapered lreaxFF case
!                     
                      gam = tp/(rij2*rij + gammaij)**third
                      if (lmixed) then
                        if (lfixQi) then
                          z(nfj) = z(nfj) - gam*reaxFFqfix(nspeci)*angstoev
                        else
                          z(nfi) = z(nfi) - gam*reaxFFqfix(nspecj)*angstoev
                        endif
                      else
                        D(numD(nfi),nfi) = D(numD(nfi),nfi) + gam
                      endif
                    enddo
                  else
                    rij = sqrt(rij2)
!
!  Compute taper function
!
                    call p7reaxFFqtaper(rij,reaxFFcutoffQ,tp,dtpdr,d2tpdr2,.false.,.false.)
!                     
!  Pure real space tapered lreaxFF case
!                     
                    gam = tp/(rij2*rij + gammaij)**third
                    if (lmixed) then
                      if (lfixQi) then
                        z(nfj) = z(nfj) - gam*reaxFFqfix(nspeci)*angstoev
                      else
                        z(nfi) = z(nfi) - gam*reaxFFqfix(nspecj)*angstoev
                      endif
                    else
                      D(numD(nfi),nfi) = D(numD(nfi),nfi) + gam
                    endif
                  endif
                endif
              endif
!
!  End of loop over j
!
            enddo jloop
          endif
        endif
        if (.not.lfixQi.and.nbos(i).gt.0) then
!
!  Scale matrix elements by conversion factor 
!
          do j = 1,numD(nfi)
            D(j,nfi) = D(j,nfi)*angstoev
          enddo
        endif
!
!  End of loop over i
!
      enddo
!
!  Now build other side of the matrix
!
      allocate(numDold(nfree),stat=status)
      if (status/=0) call outofmemory('setreaxffQiter','numDold')
      numDold(1:nfree) = numD(1:nfree)
      do i = 1,nfree
        do jj = 1,numDold(i)
          j = listD(jj,i)
          if (j.ne.i) then
            numD(j) = numD(j) + 1
          endif
        enddo
      enddo
      nD = 0
      do i = 1,nfree
        nD = max(nD,numD(i))
      enddo
      if (maxD.lt.nD) then
        maxD = nD
        call realloc(listD,maxD,nfree,ierror)
        if (ierror.ne.0) call outofmemory('setreaxffQiter','listD')
        call realloc(D,maxD,nfree,ierror)
        if (ierror.ne.0) call outofmemory('setreaxffQiter','D')
      endif
      numD(1:nfree) = numDold(1:nfree)
      do i = 1,nfree
        do jj = 1,numDold(i)
          j = listD(jj,i)
          if (j.ne.i) then
            numD(j) = numD(j) + 1
            listD(numD(j),j) = i
            D(numD(j),j) = D(jj,i)
          endif
        enddo
      enddo
      deallocate(numDold,stat=status)
      if (status/=0) call deallocate_error('setreaxffQiter','numDold')
    endif
  endif
!
  if (literative) then
!
!  Create a copy of D
!
    allocate(Dsave(maxD,nfreeloc),stat=status)
    if (status/=0) call outofmemory('setreaxffQiter','Dsave')
    do i = 1,nfreeloc
      do j = 1,numD(i)
        Dsave(j,i) = D(j,i)
      enddo
    enddo
    nitermax = nreaxFFqiter
  else
    nitermax = 1
  endif
!****************************
!  Start of iterative loop  *
!****************************
  do niter = 1,nitermax
!  
!  If this is not the first time then copy saved dpacked back
!   
    if (niter.gt.1) then
      do i = 1,nfreeloc
        do j = 1,numD(i)
          D(j,i) = Dsave(j,i)
        enddo
      enddo
    endif
!  
!  Save qreaxFF for later comparison
!   
    qreaxFFsave(1:nboatom) = qreaxFF(1:nboatom)
!********************************
!  Form matrix of coefficients  *
!********************************
    do i = 1,nfreeloc
      ii = nboatomRptr(nfreeptr(nfreelocptr(i)))
      nspeci = nbosptr(ii)
      D(1,i) = D(1,i) + 2.0_dp*reaxFFmu(nspeci)*occuf(ii)
!  
!  Shell structure option
!     
      if (abs(reaxFFshell(1,nspeci)).gt.1.0d-12) then
        if (qreaxFF(ii).gt.reaxFFshell(2,nspeci)) then
          D(1,i) = D(1,i) + 4.0_dp*reaxFFshell(1,nspeci)*occuf(ii)*(qreaxFF(ii) - reaxFFshell(2,nspeci))**3
        endif
      endif
    enddo
!
!  Debugging output
!
    if (index(keyword,'debu').ne.0) then
      if (ioproc) then
        write(ioout,'(/,''  ReaxFF Charge Matrix :'',/)')
        call gflush(ioout)
      endif
#ifdef MPI
      if (lioproconly) then
!
!  Allocate temporary workspace for communication
!
        ntmp = nfree + 1
        ntag = 1
        allocate(dtmp(ntmp),stat=status)
        if (status/=0) call outofmemory('setreaxffQiter','dtmp')
        allocate(StatMPI(MPI_Status_Size),stat=status)
        if (status/=0) call outofmemory('setreaxffQiter','StatMPI')
!
        do i = 1,nfree
          ii = nfreelocRptr(i)
          if (ioproc.or.ii.ne.0_i4) then
            if (ioproc.and.ii.ne.0_i4) then
!
!  Data is local ioproc
!
              write(ioout,'(1x,f9.5,4x,9(1x,f9.5))') z(i),(D(j,ii),j=1,numD(ii))
            else
!
!  Post receive
!
              if (ioproc) then
                nnode = nfreelocnptr(i)
                call MPI_IRecv(dtmp,ntmp,MPI_double_precision,nnode, &
                               ntag,MPI_Comm_World,Request,MPIerror)
              endif
!
!  Pass data to ioproc for writing
!
              if (ii.ne.0_i4) then
                dtmp(1:ntmp) = 0.0_dp
                dtmp(1) = dble(numD(ii))
                do j = 1,numD(ii)
                  dtmp(j+1) = D(j,ii)
                enddo
!
!  Post send
!
                call MPI_ISend(dtmp,ntmp,MPI_double_precision,0, &
                               ntag,MPI_Comm_World,Request,MPIerror)
              endif
              call MPI_WaitAll(1,Request,StatMPI,MPIerror)
              if (ioproc) then
!
!  Write on I/O node
!
                numDtmp = nint(dtmp(1))
                write(ioout,'(1x,f9.5,4x,9(1x,f9.5))') z(i),(dtmp(1+j),j=1,numDtmp)
              endif
            endif
          endif
        enddo
!
        deallocate(StatMPI,stat=status)
        if (status/=0) call deallocate_error('setreaxffQiter','StatMPI')
        deallocate(dtmp,stat=status)
        if (status/=0) call deallocate_error('setreaxffQiter','dtmp')
      else
#endif
        call mpbarrier
        do i = 1,nfree
          ii = nfreelocRptr(i)
          if (ii.gt.0_i4) then
            write(ioout,'(1x,f9.5,4x,9(1x,f9.5))') z(i),(D(j,ii),j=1,numD(ii))
          endif
          call mpbarrier
        enddo
#ifdef MPI
      endif
#endif
    endif
    if (nfree.gt.0) then
!******************
!  Invert matrix  *
!******************
!
!  Collect free charges
!
      do i = 1,nfree
        ii = nfreeptr(i)
        qreaxFFfree(i) = qreaxFF(ii)
      enddo
      qreaxFFfree(nfree+1) = 0.0_dp
!
!  Solve using iterative route
!
      call bcgsolve(nfree+1,nfreeloc,nfreelocptr,D,maxD,numD,nDdiag,listD,z,qreaxFFfree,qitertol,nqitermax,iter,err)
!
      if (ioproc) then
        if (index(keyword,'verb').ne.0) then
          write(ioout,'('' Number of iterations / error in bcgsolve = '',i4,1x,f16.14)') iter,err
          call gflush(ioout)
        endif
      endif
!
!  Was iterative solution successful?
!
      if (iter.ge.nqitermax) then
        call outerror('charge solution failed in setreaxffQiter',0_i4)
        call stopnow('setreaxffQiter')
      endif
!
!  Return free charges
!
      do i = 1,nfree
        ii = nfreeptr(i)
        qreaxFF(ii) = qreaxFFfree(i)
      enddo
    endif
    if (literative) then
!
!  Check for convergence
!
      qdiff  = 0.0_dp
      qdiff1 = 0.0_dp
      do i = 1,nboatom
        qdiff  = qdiff + abs(qreaxFF(i) - qreaxFFsave(i))
        qdiff1 = max(qdiff1,abs(qreaxFF(i) - qreaxFFsave(i)))
      enddo
      qdiff = qdiff/dble(nboatom)
      if (index(keyword,'debu').ne.0.and.ioproc) then
        write(ioout,'(''  ** Cycle : '',i4,'' Qdiffs : '',2f10.8)') niter,qdiff,qdiff1
      endif
      if (qdiff.lt.reaxFFqconverged.and.qdiff1.lt.reaxFFqconverged1) exit
!
!  If not converged then damp charges for next iteration except for first iteration
!
      if (niter.gt.1) then
        damp = reaxFFqdamp
        do i = 1,nboatom
          qreaxFF(i) = damp*qreaxFF(i) + (1.0_dp - damp)*qreaxFFsave(i)
        enddo
      endif
      if (index(keyword,'debu').ne.0.and.ioproc) then
        write(ioout,'(//,''  Charges for ReaxFF during iteration :'',/)')
        write(ioout,'(''--------------------------------------------------------------------------------'')')
        write(ioout,'(''    Atom no.            Atomic No.             Charge'')')
        write(ioout,'(''--------------------------------------------------------------------------------'')')
        do i = 1,nboatom
          ii = nboatomRptr(i)
          if (ii.gt.0) then
            write(ioout,'(6x,i4,18x,i2,16x,f10.7)') i,nat(ii),qreaxFF(i)
          endif
        enddo
        write(ioout,'(''--------------------------------------------------------------------------------'',/)')
        call gflush(ioout)
      endif
    endif
!**************************
!  End of iterative loop  *
!**************************
  enddo
!**************************
!  Calculate self energy  *
!**************************
  eself = 0.0_dp
  do i = 1,nfreeloc
    ii = nfreeptr(nfreelocptr(i))
    qi = qreaxFF(ii)
    ii = nboatomRptr(ii)
    nspeci = nbosptr(ii)
    esite = qi*occuf(ii)*(reaxFFchi(nspeci)+extpotcfg(nsft+nrelf2a(ii))+qi*reaxFFmu(nspeci))
    eself = eself + esite
    siteenergy(ii) = siteenergy(ii) + esite
!
!  Shell structure correction
!   
    if (abs(reaxFFshell(1,nspeci)).gt.1.0d-12) then
      if (qi.gt.reaxFFshell(2,nspeci)) then
        esite = occuf(ii)*reaxFFshell(1,nspeci)*(qi - reaxFFshell(2,nspeci))**4
        eself = eself + esite
        siteenergy(ii) = siteenergy(ii) + esite
      endif
    endif
  enddo
!*******************
!  Output results  *
!*******************
  if ((loutputQ.or.index(keyword,'debu').ne.0).and.ioproc) then
    write(ioout,'(//,''  Final charges from ReaxFF :'',/)')
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    write(ioout,'(''    Atom no.            Atomic No.             Charge'')')
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    do i = 1,nboatom
      ii = nboatomRptr(i)
      if (ii.gt.0) then
        write(ioout,'(4x,i6,18x,i2,16x,f10.7)') i,nat(ii),qreaxFF(i)
      endif
    enddo
    write(ioout,'(''--------------------------------------------------------------------------------'',/)')
    call gflush(ioout)
  endif
!
!  Free local memory 
!
  if (literative) then
    deallocate(Dsave,stat=status)
    if (status/=0) call deallocate_error('setreaxffQiter','Dsave')
  endif
  deallocate(qreaxFFsave,stat=status)
  if (status/=0) call deallocate_error('setreaxffQiter','qreaxFFsave')
  deallocate(qreaxFFfree,stat=status)
  if (status/=0) call deallocate_error('setreaxffQiter','qreaxFFfree')
  deallocate(z,stat=status)
  if (status/=0) call deallocate_error('setreaxffQiter','z')
  deallocate(numD,stat=status)
  if (status/=0) call deallocate_error('setreaxffQiter','numD')
  deallocate(nDdiag,stat=status)
  if (status/=0) call deallocate_error('setreaxffQiter','nDdiag')
  deallocate(nfreelocRptr,stat=status)
  if (status/=0) call deallocate_error('setreaxffQiter','nfreelocRptr')
  deallocate(nfreelocptr,stat=status)
  if (status/=0) call deallocate_error('setreaxffQiter','nfreelocptr')
  deallocate(nfreeRptr,stat=status)
  if (status/=0) call deallocate_error('setreaxffQiter','nfreeRptr')
  deallocate(nfreeptr,stat=status)
  if (status/=0) call deallocate_error('setreaxffQiter','nfreeptr')
!
  return
  end
