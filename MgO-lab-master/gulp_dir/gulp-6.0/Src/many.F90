  subroutine many(emany,esregion12,esregion2,eattach,lgrad1,lgrad2)
!
!  Subroutine for calculating the many-body energy from the
!  Sutton-Chen potential(s). Requires real space routine to
!  be called first to evaluate the repulsive contribution
!  and the rho parameters.
!
!  On entry the array scrho must contain the density at each atomic site.
!
!  lfreeze = if .true. only do derivatives for atoms with a degree of freedom
!
!   4/97 Created from reale
!   7/97 Adapted for density and functional options
!   6/00 nadd default values reduced to save CPU time as they were
!        over cautious and the many-body densities decay rapidly
!        with distance.
!   2/01 Modifications for general dimensionality added
!   4/01 Calculation of region 2 energies for surfaces added
!   9/01 Marvin compatibility option for SE added
!  11/01 Attachment energy added
!   8/02 Surface energy calculation algorithm changed
!   3/03 Speed up added to second derivative distance checking
!  11/03 ndennat/ndentyp replaced
!  11/03 Alloy changes added
!  10/04 Symmetrisation of derivatives moved to subroutine
!   7/05 Constant scaling added to sqrt and power functionals
!   7/05 Use of EAM species pointers introduced
!   9/05 rhoderv called to get density derivatives
!   9/05 Strain derivatives corrected
!   9/05 Call to eamfnderv used to replace EAM function derivative code
!  10/05 Call to eamfnderv added for energy
!   4/06 Modified to handle species specific densities
!   3/07 Printing of EAM densities and energies added as an option
!   5/07 QM/MM schemes added
!   5/07 Call to rhoderv modified by adding rpot
!  12/07 Call to rfindeither corrected
!  12/07 Unused variables removed
!  12/07 Call to rfindeither modified now that lincludeself is of intent in
!  11/08 Call to rhoderv updated to include the density
!  11/08 Call to rhoderv modified to include x,y,z Cartesian components
!  11/08 rho arrays changed to 2-D to benefit MEAM
!  11/08 call to rhoderv replaced by calls to meamrho/eamrho according to MEAM vs EAM
!  12/08 rho switched back to 1-D array with condensed components
!  12/08 scrho, and defect versions, converted to 2-D array for MEAM
!  12/08 MEAM calculation of density added as option
!   1/09 Total density now printed for MEAM case
!   1/09 Derivatives modified to accommodate MEAM : deriv -> deriv(3)
!   2/09 Derivatives for MEAM added
!   2/09 Expanding and contracting of maxdis arrays removed
!   3/09 Skipping when density is zero modified to ensure no terms are missed
!   3/09 Sign of return arguments from EAM/MEAM function routines corrected for
!   4/09 MEAM screening function derivatives added
!   5/09 MEAM third derivative arguments removed
!   9/09 Strain derivatives added to rstrd instead of rstrdnr for EAM case
!  11/09 Region derivatives added
!   4/12 Explicit virial calculation removed as no longer needed
!   5/12 Atomic stresses added
!   1/14 vectorpair saved to avoid memory leak
!   8/14 MEAM screening made species specific
!   8/14 Taper range passed to eamrho/meamrho
!   8/14 Pair potential derivatives added here so that screening can be included
!   8/14 All MEAM species (i,j,k) now passed to meanscreen
!   8/14 neamspec passed psibaskes
!   2/15 Cycling of loops for zero density removed for lMEAMden case
!   2/15 lvalidij set to be true by baskes potential
!   8/15 npartial corrected to npartialik/npartialjk in two blocks
!   9/16 Second derivatives involving i-j-k triples modified to improve speed
!   9/16 Correction to indexing of twopot in 2nd deriv cross terms
!   9/16 eamXcutfactor added
!   1/17 Unused variables (nff, nfi, nfj) removed
!   6/17 dt3/dt4 referencing restricted to MEAM code
!   2/18 Trace added
!   9/18 Modified for changes due to lstraincell algorithm
!   9/18 Strain module introduced
!  10/18 Second derivatives of baskes energy added
!  11/18 Baskes strain derivatives added to radial strain for MEAM
!  11/18 Non-radial arrays removed since they are no longer needed
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!  12/19 Rigid molecule modifications added
!   7/20 scrho now passed to meamfnderv as frho
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
  use configurations, only : nregions, nregionno, nregiontype, QMMMmode
  use control
  use current
  use derivatives
  use eam
  use general
  use iochannels,     only : ioout
  use m_strain,       only : real1strterm
  use mdlogic
  use numbers,        only : third, sixth
  use optimisation
  use progress,       only : lduring_opt
  use realvectors,    only : dist, dist2, xtmp, ytmp, ztmp, xtmp2, ytmp2, ztmp2
  use sutton
  use symmetry
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  use two
  use vectors,        only : vector_pair
  implicit none
!
!  Passed variables
!
  real(dp),    intent(inout)                   :: eattach
  real(dp),    intent(inout)                   :: emany
  real(dp),    intent(inout)                   :: esregion12
  real(dp),    intent(inout)                   :: esregion2
  logical,     intent(in)                      :: lgrad1
  logical,     intent(in)                      :: lgrad2
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ii
  integer(i4)                                  :: ind
  integer(i4)                                  :: indij
  integer(i4)                                  :: indik
  integer(i4)                                  :: indjk
  integer(i4)                                  :: is
  integer(i4)                                  :: ix
  integer(i4)                                  :: iy
  integer(i4)                                  :: iz
  integer(i4)                                  :: ixc
  integer(i4)                                  :: iyc
  integer(i4)                                  :: izc
  integer(i4)                                  :: j
  integer(i4)                                  :: jj
  integer(i4)                                  :: js
  integer(i4)                                  :: jx
  integer(i4)                                  :: jy
  integer(i4)                                  :: jz
  integer(i4)                                  :: jxc
  integer(i4)                                  :: jyc
  integer(i4)                                  :: jzc
  integer(i4)                                  :: k
  integer(i4)                                  :: km
  integer(i4)                                  :: kl
  integer(i4)                                  :: l
  integer(i4)                                  :: lvec
  integer(i4)                                  :: lxc
  integer(i4)                                  :: lyc
  integer(i4)                                  :: lzc
  integer(i4)                                  :: mp
  integer(i4)                                  :: n
  integer(i4)                                  :: nat1
  integer(i4)                                  :: nat2
  integer(i4)                                  :: nati
  integer(i4)                                  :: natj
  integer(i4)                                  :: natk
  integer(i4)                                  :: neamfnspec2
  integer(i4)                                  :: neamspeci
  integer(i4)                                  :: neamspecj
  integer(i4)                                  :: neamspeck
  integer(i4)                                  :: neamspecl
  integer(i4)                                  :: nmolonly
  integer(i4)                                  :: nor
  integer(i4)                                  :: nor2
  integer(i4)                                  :: np1
  integer(i4)                                  :: npartial
  integer(i4)                                  :: npartialik
  integer(i4)                                  :: npartialjk
  integer(i4)                                  :: npot
  integer(i4)                                  :: npotik
  integer(i4)                                  :: npotjk
  integer(i4)                                  :: npots
  integer(i4)                                  :: ns
  integer(i4)                                  :: nregioni
  integer(i4)                                  :: nregionj
  integer(i4)                                  :: nregionl
  integer(i4)                                  :: nregiontypi
  integer(i4)                                  :: nregiontypj
  integer(i4)                                  :: ntyp1  
  integer(i4)                                  :: ntyp2
  integer(i4)                                  :: ntypi
  integer(i4)                                  :: ntypj    
  integer(i4)                                  :: ntypk
  integer(i4), dimension(:), allocatable       :: npotikptr
  integer(i4), dimension(:), allocatable       :: npotjkptr
  integer(i4), dimension(:), allocatable       :: npotl
  integer(i4)                                  :: nvec
  integer(i4)                                  :: nvec0
  integer(i4)                                  :: status
  logical                                      :: lanyvalidik
  logical                                      :: lanyvalidjk
  logical,                                save :: lfirstcall = .true.
  logical                                      :: lnonzeroSij
  logical                                      :: lnonzeroSik
  logical                                      :: lnonzeroSjk
  logical                                      :: lopi
  logical                                      :: lopj
  logical                                      :: lopl
  logical                                      :: lpartial
  logical                                      :: lQMMMok
  logical                                      :: lreg2ok
  logical                                      :: lself 
  logical                                      :: lsg1 
  logical                                      :: lsg2
  logical                                      :: lvalidij
  logical                                      :: lvalidik
  logical                                      :: lvalidjk
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: cut2
  real(dp)                                     :: cut2Xscale
  real(dp)                                     :: cut2k
  real(dp)                                     :: cut2p
  real(dp)                                     :: cut2r
  real(dp)                                     :: cut2rk
  real(dp)                                     :: deltarho(maxmeamcomponent)
  real(dp)                                     :: deriv(3)
  real(dp)                                     :: derivs(6)
  real(dp)                                     :: deriv2(6)
  real(dp)                                     :: deriv2s(21)
  real(dp)                                     :: deriv2m(6,3)
  real(dp)                                     :: drhoij(3,maxmeamcomponent)
  real(dp)                                     :: drhoijs(6,maxmeamcomponent)
  real(dp)                                     :: drhoij2(6,maxmeamcomponent)
  real(dp)                                     :: drhoij2s(21,maxmeamcomponent)
  real(dp)                                     :: drhoij2m(6,3,maxmeamcomponent)
  real(dp)                                     :: drhoik(3,maxmeamcomponent)
  real(dp)                                     :: drhoiks(6,maxmeamcomponent)
  real(dp)                                     :: drhoik2(6,maxmeamcomponent)
  real(dp)                                     :: drhoik2s(21,maxmeamcomponent)
  real(dp)                                     :: drhoik2m(6,3,maxmeamcomponent)
  real(dp)                                     :: drhoji(3,maxmeamcomponent)
  real(dp)                                     :: drhojis(6,maxmeamcomponent)
  real(dp)                                     :: drhoji2(6,maxmeamcomponent)
  real(dp)                                     :: drhoji2s(21,maxmeamcomponent)
  real(dp)                                     :: drhoji2m(6,3,maxmeamcomponent)
  real(dp)                                     :: drhojk(3,maxmeamcomponent)
  real(dp)                                     :: drhojks(6,maxmeamcomponent)
  real(dp)                                     :: drhojk2(6,maxmeamcomponent)
  real(dp)                                     :: drhojk2s(21,maxmeamcomponent)
  real(dp)                                     :: drhojk2m(6,3,maxmeamcomponent)
  real(dp)                                     :: drhoki(3,maxmeamcomponent)
  real(dp)                                     :: drhokis(6,maxmeamcomponent)
  real(dp)                                     :: drhoki2(6,maxmeamcomponent)
  real(dp)                                     :: drhoki2s(21,maxmeamcomponent)
  real(dp)                                     :: drhoki2m(6,3,maxmeamcomponent)
  real(dp)                                     :: drhokj(3,maxmeamcomponent)
  real(dp)                                     :: drhokjs(6,maxmeamcomponent)
  real(dp)                                     :: drhokj2(6,maxmeamcomponent)
  real(dp)                                     :: drhokj2s(21,maxmeamcomponent)
  real(dp)                                     :: drhokj2m(6,3,maxmeamcomponent)
  real(dp)                                     :: drhototij(3)
  real(dp)                                     :: drhototijs(6)
  real(dp)                                     :: drhototij2(6)
  real(dp)                                     :: drhototij2s(21)
  real(dp)                                     :: drhototij2m(6,3)
  real(dp)                                     :: drhototij3(10)
  real(dp)                                     :: drhototik(3)
  real(dp)                                     :: drhototiksum(3)
  real(dp)                                     :: drhototiks(6)
  real(dp)                                     :: drhototikssum(6)
  real(dp)                                     :: drhototik2(6)
  real(dp)                                     :: drhototik2s(21)
  real(dp)                                     :: drhototik2m(6,3)
  real(dp)                                     :: drhototik3(10)
  real(dp)                                     :: drhototji(3)
  real(dp)                                     :: drhototjis(6)
  real(dp)                                     :: drhototji2(6)
  real(dp)                                     :: drhototji2s(21)
  real(dp)                                     :: drhototji2m(6,3)
  real(dp)                                     :: drhototji3(10)
  real(dp)                                     :: drhototjk(3)
  real(dp)                                     :: drhototjksum(3)
  real(dp)                                     :: drhototjks(6)
  real(dp)                                     :: drhototjkssum(6)
  real(dp)                                     :: drhototjk2(6)
  real(dp)                                     :: drhototjk2s(21)
  real(dp)                                     :: drhototjk2m(6,3)
  real(dp)                                     :: drhototjk3(10)
  real(dp)                                     :: drhototki(3)
  real(dp)                                     :: drhototkis(6)
  real(dp)                                     :: drhototki2(6)
  real(dp)                                     :: drhototki2s(21)
  real(dp)                                     :: drhototki2m(6,3)
  real(dp)                                     :: drhototki3(10)
  real(dp)                                     :: drhototkj(3)
  real(dp)                                     :: drhototkjs(6)
  real(dp)                                     :: drhototkj2(6)
  real(dp)                                     :: drhototkj2s(21)
  real(dp)                                     :: drhototkj2m(6,3)
  real(dp)                                     :: drhototkj3(10)
  real(dp)                                     :: drhototijk2(3,3)
  real(dp)                                     :: drhototijk2sum(3,3)
  real(dp)                                     :: drhototjik2(3,3)
  real(dp)                                     :: drhototjik2sum(3,3)
  real(dp)                                     :: drhototkij2(3,3)
  real(dp)                                     :: drhototijk2s(6,6)
  real(dp)                                     :: drhototijk2ssum(6,6)
  real(dp)                                     :: drhototjik2s(6,6)
  real(dp)                                     :: drhototjik2ssum(6,6)
  real(dp)                                     :: drhototkij2s(6,6)
  real(dp)                                     :: drhototijk2m(6,3)
  real(dp)                                     :: drhototijk2msum(6,3)
  real(dp)                                     :: drhototjik2m(6,3)
  real(dp)                                     :: drhototjik2msum(6,3)
  real(dp)                                     :: drhototkij2m(6,3)
  real(dp)                                     :: dr2ds(6)
  real(dp)                                     :: d2r2dx2(6)
  real(dp)                                     :: d2r2ds2(6,6)
  real(dp)                                     :: d2r2dsdx(6,3)
  real(dp)                                     :: dt1
  real(dp)                                     :: dt1l
  real(dp)                                     :: dt2
  real(dp)                                     :: dt2l
  real(dp)                                     :: dt3
  real(dp)                                     :: dt3l
  real(dp)                                     :: dt4
  real(dp)                                     :: dt4l
  real(dp)                                     :: ebas
  real(dp)                                     :: d1bas
  real(dp)                                     :: d2bas
  real(dp)                                     :: eeam
  real(dp)                                     :: elcom(6,6)
  real(dp)                                     :: emanytrm
  real(dp)                                     :: frho(maxmeamcomponent)
  real(dp)                                     :: oci      
  real(dp)                                     :: ocj  
  real(dp)                                     :: ock
  real(dp)                                     :: ofct
  real(dp)                                     :: ofctijk
  real(dp)                                     :: oij
  real(dp)                                     :: r
  real(dp)                                     :: r2
  real(dp)                                     :: rcut2
  real(dp)                                     :: rcutfactor
  real(dp)                                     :: rhoi
  real(dp)                                     :: rhoj
  real(dp)                                     :: rhok
  real(dp)                                     :: rhoij(maxmeamcomponent)
  real(dp)                                     :: rhoji(maxmeamcomponent)
  real(dp)                                     :: rhoik(maxmeamcomponent)
  real(dp)                                     :: rhoki(maxmeamcomponent)
  real(dp)                                     :: rhojk(maxmeamcomponent)
  real(dp)                                     :: rhokj(maxmeamcomponent)
  real(dp)                                     :: rhorr  
  real(dp)                                     :: rik
  real(dp)                                     :: rjk
  real(dp)                                     :: rik2
  real(dp)                                     :: rjk2
  real(dp)                                     :: rk
  real(dp)                                     :: rp
  real(dp)                                     :: rpik
  real(dp)                                     :: rpjk
  real(dp)                                     :: rpijk
  real(dp)                                     :: rscrhoi
  real(dp)                                     :: rscrhoi3
  real(dp)                                     :: rscrhoi5
  real(dp)                                     :: rscrhoj
  real(dp)                                     :: rscrhoj3
  real(dp)                                     :: rscrhoj5
  real(dp)                                     :: rscrhok
  real(dp)                                     :: rscrhok3
  real(dp)                                     :: rscrhok5
  real(dp)                                     :: scmax
  real(dp)                                     :: Sij
  real(dp)                                     :: Sik
  real(dp)                                     :: Sjk
  real(dp)                                     :: Silj
  real(dp)                                     :: Silk
  real(dp)                                     :: Sjlk
  real(dp)                                     :: dSiljdr(3)
  real(dp)                                     :: dSilkdr(3)
  real(dp)                                     :: dSjlkdr(3)
  real(dp)                                     :: time1
  real(dp)                                     :: time2  
  real(dp)                                     :: rstrl(6)
  real(dp)                                     :: xal 
  real(dp)                                     :: yal    
  real(dp)                                     :: zal
  real(dp)                                     :: xcd 
  real(dp)                                     :: ycd    
  real(dp)                                     :: zcd
  real(dp)                                     :: xcd1
  real(dp)                                     :: ycd1   
  real(dp)                                     :: zcd1
  real(dp)                                     :: xcd2
  real(dp)                                     :: ycd2   
  real(dp)                                     :: zcd2
  real(dp)                                     :: xcrd 
  real(dp)                                     :: ycrd    
  real(dp)                                     :: zcrd
  real(dp)                                     :: xcrdik
  real(dp)                                     :: ycrdik
  real(dp)                                     :: zcrdik
  real(dp)                                     :: xcrdjk
  real(dp)                                     :: ycrdjk
  real(dp)                                     :: zcrdjk
  real(dp)                                     :: xil0
  real(dp)                                     :: yil0
  real(dp)                                     :: zil0
  real(dp)                                     :: xjl0
  real(dp)                                     :: yjl0
  real(dp)                                     :: zjl0
  real(dp)                                     :: xkl0
  real(dp)                                     :: ykl0
  real(dp)                                     :: zkl0
  type(screening_atoms)                        :: partial
  type(screening_atoms)                        :: partialik
  type(screening_atoms)                        :: partialjk
  type(vector_pair),                      save :: vectorpair
#ifdef TRACE
  call trace_in('many')
#endif
!
  time1 = g_cpu_time()
!
!  If this is the first call then initialise vectorpair
!
  if (lfirstcall) then
    lfirstcall = .false.
    call changemaxvectorpair(vectorpair,0_i4)
  endif
!
!  For screened MEAM, set scale factor that determines searching cutoff based on which axis of the ellipse is largest.
!  Note: the cutoff is applied to the mid point of the i-j vector and so this is guaranteed to find all distances.
!
  rcutfactor = 1.0_dp
  if (lanyMEAMscreen) then
    neamfnspec2 = neamfnspec*(neamfnspec+1)/2
    do i = 1,neamfnspec
      do j = 1,neamfnspec2
        if (lMEAMscreen(j,i)) then
          rcutfactor = max(1.0_dp,0.25_dp*(1.0_dp + meam_Cmax(j,i)))
        endif
      enddo
    enddo
  endif
!
!  Scale density
!
  call eamscalescrho(1_i4)
  if (lPrintEAM) then
!
!  Opening banner for energy decomposition
!
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    write(ioout,'(''  EAM : Atom No.                Density                 Atom energy (eV) '')')
    write(ioout,'(''--------------------------------------------------------------------------------'')')
  endif
!
!  Energy calculation
!
  emany = 0.0_dp
  lreg2ok = (lseok.and.nregions(ncf).gt.1)
  do i = 1,numat
    neamspeci = neamfnspecptr(i)
    rhoi = scrho(1,i)
    nregioni = nregionno(nsft+nrelf2a(i))
    nregiontypi = nregiontype(nregioni,ncf)
!     
!  QM/MM handling : i is a QM atom => exclude
!     
    lQMMMok = .true.
    if (QMMMmode(ncf).gt.0) then
      if (nregiontypi.eq.1) lQMMMok = .false.
    endif
    if (neamspeci.gt.0.and.rhoi.gt.1.0d-12.and.lQMMMok) then
      if (lMEAMfn) then
        frho(1:maxmeamcomponent) = scrho(1:maxmeamcomponent,i)
        call meamfnderv(neamfn,neamspeci,frho,rhoi,eeam,rscrhoi,rscrhoi3,rscrhoi5,.false.,.false.,.false.)
      else
        call eamfnderv(neamfn,neamspeci,scrho(1,i),eeam,rscrhoi,rscrhoi3,rscrhoi5,.false.,.false.,.false.)
      endif
      emanytrm = occuf(i)*eeam
      if (lreg2ok.and.nregionno(nsft+nrelf2a(i)).gt.1) then
        esregion2 = esregion2 + emanytrm
      else
        emany = emany + emanytrm
      endif
      if (lPrintEAM) then
        write(ioout,'(7x,i8,4x,f24.10,4x,f24.12)') i,rhoi,emanytrm
      endif
      eattach = eattach + emanytrm
      if (lMEAMfn) then
        deltarho(1:maxmeamcomponent) = scrho(1:maxmeamcomponent,i) - scrho12(1:maxmeamcomponent,i)
        call meamfnderv(neamfn,neamspeci,deltarho,rhoi,eeam,rscrhoi,rscrhoi3,rscrhoi5,.false.,.false.,.false.)
      else
        rhorr = rhoi - scrho12(1,i)
        call eamfnderv(neamfn,neamspeci,rhorr,eeam,rscrhoi,rscrhoi3,rscrhoi5,.false.,.false.,.false.)
      endif
      emanytrm = occuf(i)*eeam
      eattach = eattach - emanytrm
    endif
  enddo
  if (lPrintEAM) then
!
!  Closing banner for energy decomposition
!
    write(ioout,'(''--------------------------------------------------------------------------------'')')
  endif
!
!  If no forces are needed then we don't need to do loops over atoms, so just return
!
  if (.not.lgrad1) goto 1000
!
  if (lduring_opt) then
    cut2Xscale = eamXcutfactor**2
  else
    cut2Xscale = 1.0_dp
  endif
!
!  From here on we can assume that lgrad1 = .true.
!
!  Local variables
!
  lsg1 = lstr
  lsg2 = (lgrad2.and.lstr)
  if (lsg1) then                               
    do i = 1,nstrains
      rstrl(i) = 0.0_dp
    enddo
    if (lsg2) then                             
      do i = 1,nstrains                               
        do j = 1,nstrains                             
          elcom(j,i) = 0.0_dp                  
        enddo                                  
      enddo                                    
    endif
  endif
!***************************
!  Set up local variables  *
!***************************
!
!  Find maximum cut-off radius
!
  scmax = 0.0_dp
  do i = 1,npote
    if (nptype(i).eq.19.or.nptype(i).eq.45.or.nptype(i).eq.55) then
      if (rpot(i).gt.scmax) scmax = rpot(i)
    endif
  enddo
!
!  Allocate local memory
!
  allocate(npotl(npote),stat=status)
  if (status/=0) call outofmemory('many','npotl')
  if (lgrad2) then
    allocate(npotikptr(npote),stat=status)
    if (status/=0) call outofmemory('many','npotikptr')
    allocate(npotjkptr(npote),stat=status)
    if (status/=0) call outofmemory('many','npotjkptr')
  endif
!
!  Set up cutoffs
!
  cut2p = cutp*cutp
  if (lnoreal) goto 999
!
!  Outer loop over sites
!
  ixc = - 2
  iyc = - 1
  izc = 0
  iloop: do i = 1,numat
!
!  Handle ixc-izc before rho check otherwise
!  values go wrong
!
    lopi = (.not.lfreeze.or.lopf(nrelf2a(i)))
    if (lopi) then
      ixc = ixc + 3
      iyc = iyc + 3
      izc = izc + 3
    endif
!
    xal = xclat(i)
    yal = yclat(i)
    zal = zclat(i)
    nati = nat(i)
    ntypi = nftype(i)
    nregioni = nregionno(nsft+nrelf2a(i))
    nregiontypi = nregiontype(nregioni,ncf)
    oci = occuf(i)
!     
!  Find EAM species for i
!  
    neamspeci = neamfnspecptr(i)
!
!  Evaluate functional derivatives
!
    if (lMEAMfn) then
      frho(1:maxmeamcomponent) = scrho(1:maxmeamcomponent,i)
      call meamfnderv(neamfn,neamspeci,frho,rhoi,eeam,rscrhoi,rscrhoi3,rscrhoi5,.true.,lgrad2,.false.)
    else
      call eamfnderv(neamfn,neamspeci,scrho(1,i),eeam,rscrhoi,rscrhoi3,rscrhoi5,.true.,lgrad2,.false.)
      rhoi = scrho(1,i)
    endif
!
!  Start of second atom loop
!
    jxc = - 2
    jyc = - 1
    jzc = 0
    jloop: do j = 1,i
!
!  Freezing flag - need to do this before checking rho values otherwise jxc-jzc go wrong
!
      lopj = (.not.lfreeze.or.lopf(nrelf2a(j)))
      if (.not.lopi.and..not.lopj) cycle jloop
      if (lopj) then
        jxc = jxc + 3
        jyc = jyc + 3
        jzc = jzc + 3
      endif
!
      natj = nat(j)
      ntypj = nftype(j)
      nregionj = nregionno(nsft+nrelf2a(j))
      nregiontypj = nregiontype(nregionj,ncf)
!     
!  QM/MM handling : i & j are both QM atoms => no forces to compute
!     
      if (QMMMmode(ncf).gt.0) then
        if (nregiontypi.eq.1.and.nregiontypj.eq.1) cycle jloop
      endif
!     
!  Find EAM species for j
!  
      neamspecj = neamfnspecptr(j)
!
!  Find index for i-j
!
      if (neamspeci.gt.neamspecj) then
        indij = neamspeci*(neamspeci - 1)/2 + neamspecj
      else
        indij = neamspecj*(neamspecj - 1)/2 + neamspeci
      endif
!
!  Evaluate functional derivatives
!
      if (lMEAMfn) then
        frho(1:maxmeamcomponent) = scrho(1:maxmeamcomponent,j)
        call meamfnderv(neamfn,neamspecj,frho,rhoj,eeam,rscrhoj,rscrhoj3,rscrhoj5,.true.,lgrad2,.false.)
      else
        call eamfnderv(neamfn,neamspecj,scrho(1,j),eeam,rscrhoj,rscrhoj3,rscrhoj5,.true.,lgrad2,.false.)
        rhoj = scrho(1,j)
      endif
!
!  If there is no density at either centre and this is not a second derivative run then cycle
!
      if (.not.lMEAMden.and.rhoi.eq.0.0_dp.and.rhoj.eq.0.0_dp.and..not.lgrad2) cycle jloop
!
!  Arrange nat1 to be lower than nat2. If they are equal
!  then arrange that ntyp1 is lower than ntyp2
!
      if (nati.eq.natj) then
        nat1 = nati
        nat2 = natj
        if (ntypi.lt.ntypj) then
          ntyp1 = ntypi
          ntyp2 = ntypj
        else
          ntyp1 = ntypj
          ntyp2 = ntypi
        endif
      elseif (nati.lt.natj) then
        nat1 = nati
        nat2 = nat(j)
        ntyp1 = ntypi
        ntyp2 = nftype(j)
      else
        nat1 = nat(j)
        nat2 = nati
        ntyp1 = nftype(j)
        ntyp2 = ntypi
      endif
      xcrd = xclat(j) - xal
      ycrd = yclat(j) - yal
      zcrd = zclat(j) - zal
      ocj = occuf(j)
!
!  Set flags such that if one atom is not being optimised place 
!  3 x 3 second derivative matrix in the on-diagonal block
!
      if (lopi.and.lopj) then
        ix = ixc
        iy = iyc
        iz = izc
        jx = jxc
        jy = jyc
        jz = jzc
      elseif (lopi) then
        ix = ixc
        iy = iyc
        iz = izc
        jx = ixc
        jy = iyc
        jz = izc
      elseif (lopj) then
        ix = jxc
        iy = jyc
        iz = jzc
        jx = jxc
        jy = jyc
        jz = jzc
      endif
      ofct = oci*ocj
!
!  Locate potential number
!  Check whether potential requires specific types
!  Calculate sum of all dispersion terms for pair
!
      npots = 0
      rp = 0.0_dp
      do n = 1,npote
        if (nptype(n).eq.19.or.nptype(n).eq.45.or.nptype(n).eq.55) then
          if (nat1.eq.nspec1(n).and.nat2.eq.nspec2(n)) then
            if ((ntyp1.eq.nptyp1(n).or.nptyp1(n).eq.0).and.(ntyp2.eq.nptyp2(n).or.nptyp2(n).eq.0)) then
              npots = npots + 1
              npotl(npots) = n
              if (rpot(n).gt.rp) rp = rpot(n)
            endif
          endif
        endif
      enddo
!
!  If no valid potentials then there is no need
!  to continue with this pair, unless this is
!  a second derivative calculation, in which
!  case there may be a contribution from triangles
!  of interactions.
!
      if (npots.eq.0.and..not.lgrad2) cycle jloop
      if (lgrad2) then
!
!  Need to make cut-off equal to double the maximum
!  to ensure all triangles are included
!
        rp = 2.0_dp*scmax
        cut2r = rp*rp
        if (cut2r.gt.4.0_dp*cut2p) cut2r = cut2p
      else
        cut2r = rp*rp
        if (cut2r.gt.cut2p) cut2r = cut2p
      endif
      cut2 = cut2r
      rp = sqrt(cut2)
!*********************************
!  Find valid vectors for i - j  *
!*********************************
      call rfind(xcrd,ycrd,zcrd,cut2,0.0_dp,0.0_dp,.false.,.false.,i,j,0_i4,0_i4,0_i4,0_i4,nmolonly,lself,0_i4,nor)
!
!  Loop over valid vectors
!
      do ii = 1,nor
!***************************************************
!  Calculate many-body contribution in real space  *
!***************************************************
        deriv(1:3) = 0.0_dp
        if (lsg1) then
          derivs(1:nstrains) = 0.0_dp
        endif
        if (lgrad2) then
          deriv2(1:6) = 0.0_dp
          if (lstr) then
            deriv2s(1:nstrains2) = 0.0_dp
            deriv2m(1:nstrains,1:3) = 0.0_dp
          endif
        endif
        r2 = dist(ii)
        r = sqrt(r2)
        xcd = xtmp(ii)
        ycd = ytmp(ii)
        zcd = ztmp(ii)
!***************************************
!  Valid many-body potentials for i-j  *
!***************************************
        if (lMEAM) then
          rhoij(1:maxmeamcomponent) = 0.0_dp
          rhoji(1:maxmeamcomponent) = 0.0_dp
          drhoij(1:3,1:maxmeamcomponent) = 0.0_dp
          drhoji(1:3,1:maxmeamcomponent) = 0.0_dp
          if (lstr) then
            drhoijs(1:nstrains,1:maxmeamcomponent) = 0.0_dp
            drhojis(1:nstrains,1:maxmeamcomponent) = 0.0_dp
          endif
          if (lgrad2) then
            drhoij2(1:6,1:maxmeamcomponent) = 0.0_dp
            drhoji2(1:6,1:maxmeamcomponent) = 0.0_dp
            if (lstr) then
              drhoij2s(1:nstrains2,1:maxmeamcomponent) = 0.0_dp
              drhoji2s(1:nstrains2,1:maxmeamcomponent) = 0.0_dp
              drhoij2m(1:nstrains,1:3,1:maxmeamcomponent) = 0.0_dp
              drhoji2m(1:nstrains,1:3,1:maxmeamcomponent) = 0.0_dp
            endif
          endif
        else
          rhoij(1) = 0.0_dp
          rhoji(1) = 0.0_dp
        endif
        drhototij(1:3) = 0.0_dp
        drhototji(1:3) = 0.0_dp
        drhototijs(1:nstrains) = 0.0_dp
        drhototjis(1:nstrains) = 0.0_dp
        if (lgrad2) then
          drhototij2(1:6) = 0.0_dp
          drhototji2(1:6) = 0.0_dp
          if (lstr) then
            drhototij2s(1:nstrains2) = 0.0_dp
            drhototji2s(1:nstrains2) = 0.0_dp
            drhototij2m(1:nstrains,1:3) = 0.0_dp
            drhototji2m(1:nstrains,1:3) = 0.0_dp
          endif
        endif
        lvalidij = .false.
        if (npots.gt.0) then
          if (lMEAMden) then
            do mp = 1,npots
              npot = npotl(mp)
              if (nptype(npot).eq.19) then
                if (r.gt.rpot2(npot).and.r.le.rpot(npot).and.r.le.rpmax) then
                  lvalidij = .true.
!**********************************
!  Calculate density derivatives  *
!**********************************
                  call meamrho(nati,ntypi,natj,ntypj,r,rpot(npot),xcd,ycd,zcd,rhoij,rhoji,drhoij,drhoji, &
                               drhoijs,drhojis,drhoij2,drhoji2,drhoij2s,drhoji2s,drhoij2m,drhoji2m, &
                               1.0_dp,1.0_dp,.true.,lstr,.true.,lgrad2,twopot(1,npot))
                endif
              endif
            enddo
!*******************************
!  Compute screening function  *
!*******************************
!
!  Initialise screening function to 1
!
            lnonzeroSij = .true.
            Sij = 1.0_dp
!
            if (lanyMEAMscreen) then
!
!  Set cutoffs for possible screening atoms
!
              rcut2 = rcutfactor*r2
!
!  Loop over atoms to search for images that may contribute to the screening
!
              l = 0
              npartial = 0
              lxc = - 2
              lyc = - 1
              lzc =   0
              do while (l.lt.numat.and.lnonzeroSij)
                l = l + 1
                neamspecl = neamfnspecptr(l)
                lopl = (.not.lfreeze.or.lopf(nrelf2a(l)))
                if (lopl) then
                  lxc = lxc + 3
                  lyc = lyc + 3
                  lzc = lzc + 3
                endif
                if (lMEAMscreen(indij,neamspecl)) then
!
!  Set basic vectors between atoms
!
                  xil0 = xclat(l) - xal
                  yil0 = yclat(l) - yal
                  zil0 = zclat(l) - zal
                  xjl0 = xil0 - xcd
                  yjl0 = yil0 - ycd
                  zjl0 = zil0 - zcd
!
!  Find images within cutoffs of both atoms - excluding self images
!
                  nvec0 = 0
                  call rfindmid(xil0,yil0,zil0,xjl0,yjl0,zjl0,rcut2,.false.,nvec0,nvec,vectorpair)
!
!  Loop over results of search
!
                  lvec = 0
                  do while (lvec.lt.nvec.and.lnonzeroSij)
                    lvec = lvec + 1
!
!  Compute screening function
!
                    call meamscreen(neamspeci,neamspecj,neamspecl,r2,vectorpair%distance_pair1(lvec), &
                                    vectorpair%distance_pair2(lvec),Silj,dSiljdr,lpartial,.true.)
!
!  If screening function contribution is 0, then no need to continue for this pair
!
                    if (Silj.eq.0.0_dp) then
                      lnonzeroSij = .false.
                      Sij = 0.0_dp
                    else
!
!  Multiply total screening product
!
                      Sij = Sij*Silj
                      if (lpartial) then
!
!  If this atom has a screening factor between 0 and 1, we need to keep track of it since it will generate non-zero derivatives
!
                        npartial = npartial + 1
                        if (npartial.gt.partial%sa_maxdim) then
                          call changemaxsa(partial,npartial)
                        endif
                        partial%sa_atom(npartial) = l
                        partial%sa_kxc(npartial) = lxc
                        partial%sa_xik(npartial) = vectorpair%xvector_pair1(lvec)
                        partial%sa_yik(npartial) = vectorpair%yvector_pair1(lvec)
                        partial%sa_zik(npartial) = vectorpair%zvector_pair1(lvec)
                        partial%sa_xjk(npartial) = vectorpair%xvector_pair2(lvec)
                        partial%sa_yjk(npartial) = vectorpair%yvector_pair2(lvec)
                        partial%sa_zjk(npartial) = vectorpair%zvector_pair2(lvec)
                        partial%sa_Sikj(npartial) = Silj
                        partial%sa_dSikjdr(1:3,npartial)   = dSiljdr(1:3)
                      endif
                    endif
!
!  End loop over images of possible screening atoms
!
                  enddo
                endif
!
!  End loop over possible screening atoms
!
              enddo
!
!  End of screening function
!
            endif
!
!  Only do remainder of work if the screening factor is non-zero
!
            if (lnonzeroSij) then
              do mp = 1,npots
                npot = npotl(mp)
!
!  Pair potential contribution
!
                if (nptype(npot).eq.45.or.nptype(npot).eq.55) then
                  if (r.gt.rpot2(npot).and.r.le.rpot(npot)) then
                    lvalidij = .true.
                    rk = 1.0_dp/r
                    call baskes(nati,ntypi,neamspeci,natj,ntypj,neamspecj,npot,r,rk,1.0_dp,ebas,d1bas,d2bas,.true.,lgrad2)
                    ebas  = ebas*ofct
                    d1bas = rk*d1bas*ofct
!
                    deriv(1) = deriv(1) + d1bas*xcd*Sij
                    deriv(2) = deriv(2) + d1bas*ycd*Sij
                    deriv(3) = deriv(3) + d1bas*zcd*Sij
!
                    if (lstr.or.lgrad2) then
                      call real1strterm(ndim,xcd,ycd,zcd,0.0_dp,0.0_dp,0.0_dp,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,lgrad2)
                    endif
!
                    if (lgrad2) then
                      d2bas = d2bas*ofct
                      d2bas = rk*rk*(d2bas - d1bas)
!
                      deriv2(1) = deriv2(1) + d2bas*Sij*d2r2dx2(1) + d1bas*Sij
                      deriv2(2) = deriv2(2) + d2bas*Sij*d2r2dx2(6)
                      deriv2(3) = deriv2(3) + d2bas*Sij*d2r2dx2(5)
                      deriv2(4) = deriv2(4) + d2bas*Sij*d2r2dx2(2) + d1bas*Sij
                      deriv2(5) = deriv2(5) + d2bas*Sij*d2r2dx2(4)
                      deriv2(6) = deriv2(6) + d2bas*Sij*d2r2dx2(3) + d1bas*Sij
                    endif
                    if (lstr) then
                      do is = 1,nstrains
                        ns = nstrptr(is)
                        derivs(is) = derivs(is) + d1bas*Sij*dr2ds(ns)
                      enddo
                      if (lgrad2) then
!
!  Strain-strain derivatives
!
                        ind = 0
                        do is = 1,nstrains
                          kl = nstrptr(is)
                          do js = 1,is
                            km = nstrptr(js)
                            ind = ind + 1
                            deriv2s(ind) = deriv2s(ind) + d2bas*Sij*dr2ds(kl)*dr2ds(km)
                            deriv2s(ind) = deriv2s(ind) + d1bas*Sij*d2r2ds2(km,kl)
                          enddo
                        enddo
!
!  Mixed strain-Cartesian derivatives
!
                        deriv2m(1:nstrains,1:3) = deriv2m(1:nstrains,1:3) + d1bas*Sij*d2r2dsdx(1:nstrains,1:3)
                        deriv2m(1:nstrains,1) = deriv2m(1:nstrains,1) + xcd*d2bas*Sij*dr2ds(1:nstrains)
                        deriv2m(1:nstrains,2) = deriv2m(1:nstrains,2) + ycd*d2bas*Sij*dr2ds(1:nstrains)
                        deriv2m(1:nstrains,3) = deriv2m(1:nstrains,3) + zcd*d2bas*Sij*dr2ds(1:nstrains)
                      endif
                    endif
!
                    call psiscreenderv(i,j,npartial,partial,xcd,ycd,zcd,ebas,Sij,lstr)
                  endif
                endif
!
!  Density contribution
!
                if (nptype(npot).eq.19) then
                  if (lanyMEAMscreen) then
                    if (npartial.gt.0) then
!
!  Compute derivative contributions of the screening function w.r.t. total density of atom i
!
                      call meamtotalrhoscreenderv(neamspeci,npartial,partial,xcd,ycd,zcd,scrho(1,i), &
                                                  rhoi,rscrhoi,rhoij,Sij,lsg1)
!
                      do np1 = 1,npartial
                        l = partial%sa_atom(np1)
                        lopl = (.not.lfreeze.or.lopf(nrelf2a(l)))
                        nregionl = nregionno(nsft+nrelf2a(l))
!
!  i-j contribution
!
                        if (lopi) then
                          xdrv(i) = xdrv(i) - partial%sa_drhototij(1,np1)*ofct
                          ydrv(i) = ydrv(i) - partial%sa_drhototij(2,np1)*ofct
                          zdrv(i) = zdrv(i) - partial%sa_drhototij(3,np1)*ofct
                        endif
                        if (lopj) then
                          xdrv(j) = xdrv(j) + partial%sa_drhototij(1,np1)*ofct
                          ydrv(j) = ydrv(j) + partial%sa_drhototij(2,np1)*ofct
                          zdrv(j) = zdrv(j) + partial%sa_drhototij(3,np1)*ofct
                        endif
                        if (nregioni.ne.nregionj) then
                          xregdrv(nregioni) = xregdrv(nregioni) - partial%sa_drhototij(1,np1)*ofct
                          yregdrv(nregioni) = yregdrv(nregioni) - partial%sa_drhototij(2,np1)*ofct
                          zregdrv(nregioni) = zregdrv(nregioni) - partial%sa_drhototij(3,np1)*ofct
                          xregdrv(nregionj) = xregdrv(nregionj) + partial%sa_drhototij(1,np1)*ofct
                          yregdrv(nregionj) = yregdrv(nregionj) + partial%sa_drhototij(2,np1)*ofct
                          zregdrv(nregionj) = zregdrv(nregionj) + partial%sa_drhototij(3,np1)*ofct
                        endif
!
!  i-l contribution
!
                        if (lopi) then
                          xdrv(i) = xdrv(i) - partial%sa_drhototik(1,np1)*ofct
                          ydrv(i) = ydrv(i) - partial%sa_drhototik(2,np1)*ofct
                          zdrv(i) = zdrv(i) - partial%sa_drhototik(3,np1)*ofct
                        endif
                        if (lopl) then
                          xdrv(l) = xdrv(l) + partial%sa_drhototik(1,np1)*ofct
                          ydrv(l) = ydrv(l) + partial%sa_drhototik(2,np1)*ofct
                          zdrv(l) = zdrv(l) + partial%sa_drhototik(3,np1)*ofct
                        endif
                        if (nregioni.ne.nregionl) then
                          xregdrv(nregioni) = xregdrv(nregioni) - partial%sa_drhototik(1,np1)*ofct
                          yregdrv(nregioni) = yregdrv(nregioni) - partial%sa_drhototik(2,np1)*ofct
                          zregdrv(nregioni) = zregdrv(nregioni) - partial%sa_drhototik(3,np1)*ofct
                          xregdrv(nregionl) = xregdrv(nregionl) + partial%sa_drhototik(1,np1)*ofct
                          yregdrv(nregionl) = yregdrv(nregionl) + partial%sa_drhototik(2,np1)*ofct
                          zregdrv(nregionl) = zregdrv(nregionl) + partial%sa_drhototik(3,np1)*ofct
                        endif
!
!  j-l contribution
!
                        if (lopj) then
                          xdrv(j) = xdrv(j) - partial%sa_drhototjk(1,np1)*ofct
                          ydrv(j) = ydrv(j) - partial%sa_drhototjk(2,np1)*ofct
                          zdrv(j) = zdrv(j) - partial%sa_drhototjk(3,np1)*ofct
                        endif
                        if (lopl) then
                          xdrv(l) = xdrv(l) + partial%sa_drhototjk(1,np1)*ofct
                          ydrv(l) = ydrv(l) + partial%sa_drhototjk(2,np1)*ofct
                          zdrv(l) = zdrv(l) + partial%sa_drhototjk(3,np1)*ofct
                        endif
                        if (nregionj.ne.nregionl) then
                          xregdrv(nregionj) = xregdrv(nregionj) - partial%sa_drhototjk(1,np1)*ofct
                          yregdrv(nregionj) = yregdrv(nregionj) - partial%sa_drhototjk(2,np1)*ofct
                          zregdrv(nregionj) = zregdrv(nregionj) - partial%sa_drhototjk(3,np1)*ofct
                          xregdrv(nregionl) = xregdrv(nregionl) + partial%sa_drhototjk(1,np1)*ofct
                          yregdrv(nregionl) = yregdrv(nregionl) + partial%sa_drhototjk(2,np1)*ofct
                          zregdrv(nregionl) = zregdrv(nregionl) + partial%sa_drhototjk(3,np1)*ofct
                        endif
                        if (lstr) then
                          do kl = 1,nstrains
                            rstrl(kl) = rstrl(kl) + partial%sa_drhotots(kl,np1)*ofct
                          enddo
                          if (latomicstress) then
                            if (.not.lfreeze) then
                              do kl = 1,nstrains
                                atomicstress(kl,i) = atomicstress(kl,i) + sixth*partial%sa_drhotots(kl,np1)*ofct
                                atomicstress(kl,j) = atomicstress(kl,j) + sixth*partial%sa_drhotots(kl,np1)*ofct
                                atomicstress(kl,k) = atomicstress(kl,k) + sixth*partial%sa_drhotots(kl,np1)*ofct
                              enddo
                            else
                              do kl = 1,nstrains
                                atomicstress(kl,i) = atomicstress(kl,i) + third*partial%sa_drhotots(kl,np1)*ofct
                                atomicstress(kl,j) = atomicstress(kl,j) + third*partial%sa_drhotots(kl,np1)*ofct
                                atomicstress(kl,k) = atomicstress(kl,k) + third*partial%sa_drhotots(kl,np1)*ofct
                              enddo
                            endif
                          endif
                        endif
                      enddo
!
!  Compute derivative contributions of the screening function w.r.t. total density of atom j
!
                      call meamtotalrhoscreenderv(neamspecj,npartial,partial,xcd,ycd,zcd,scrho(1,j), &
                                                  rhoj,rscrhoj,rhoji,Sij,lsg1)
!
                      do np1 = 1,npartial
                        l = partial%sa_atom(np1)
                        lopl = (.not.lfreeze.or.lopf(nrelf2a(l)))
                        nregionl = nregionno(nsft+nrelf2a(l))
!
!  i-j contribution
!
                        if (lopi) then
                          xdrv(i) = xdrv(i) - partial%sa_drhototij(1,np1)*ofct
                          ydrv(i) = ydrv(i) - partial%sa_drhototij(2,np1)*ofct
                          zdrv(i) = zdrv(i) - partial%sa_drhototij(3,np1)*ofct
                        endif
                        if (lopj) then
                          xdrv(j) = xdrv(j) + partial%sa_drhototij(1,np1)*ofct
                          ydrv(j) = ydrv(j) + partial%sa_drhototij(2,np1)*ofct
                          zdrv(j) = zdrv(j) + partial%sa_drhototij(3,np1)*ofct
                        endif
                        if (nregioni.ne.nregionj) then
                          xregdrv(nregioni) = xregdrv(nregioni) - partial%sa_drhototij(1,np1)*ofct
                          yregdrv(nregioni) = yregdrv(nregioni) - partial%sa_drhototij(2,np1)*ofct
                          zregdrv(nregioni) = zregdrv(nregioni) - partial%sa_drhototij(3,np1)*ofct
                          xregdrv(nregionj) = xregdrv(nregionj) + partial%sa_drhototij(1,np1)*ofct
                          yregdrv(nregionj) = yregdrv(nregionj) + partial%sa_drhototij(2,np1)*ofct
                          zregdrv(nregionj) = zregdrv(nregionj) + partial%sa_drhototij(3,np1)*ofct
                        endif
!
!  i-l contribution
!
                        if (lopi) then
                          xdrv(i) = xdrv(i) - partial%sa_drhototik(1,np1)*ofct
                          ydrv(i) = ydrv(i) - partial%sa_drhototik(2,np1)*ofct
                          zdrv(i) = zdrv(i) - partial%sa_drhototik(3,np1)*ofct
                        endif
                        if (lopl) then
                          xdrv(l) = xdrv(l) + partial%sa_drhototik(1,np1)*ofct
                          ydrv(l) = ydrv(l) + partial%sa_drhototik(2,np1)*ofct
                          zdrv(l) = zdrv(l) + partial%sa_drhototik(3,np1)*ofct
                        endif
                        if (nregioni.ne.nregionl) then
                          xregdrv(nregioni) = xregdrv(nregioni) - partial%sa_drhototik(1,np1)*ofct
                          yregdrv(nregioni) = yregdrv(nregioni) - partial%sa_drhototik(2,np1)*ofct
                          zregdrv(nregioni) = zregdrv(nregioni) - partial%sa_drhototik(3,np1)*ofct
                          xregdrv(nregionl) = xregdrv(nregionl) + partial%sa_drhototik(1,np1)*ofct
                          yregdrv(nregionl) = yregdrv(nregionl) + partial%sa_drhototik(2,np1)*ofct
                          zregdrv(nregionl) = zregdrv(nregionl) + partial%sa_drhototik(3,np1)*ofct
                        endif
!
!  j-l contribution
!
                        if (lopj) then
                          xdrv(j) = xdrv(j) - partial%sa_drhototjk(1,np1)*ofct
                          ydrv(j) = ydrv(j) - partial%sa_drhototjk(2,np1)*ofct
                          zdrv(j) = zdrv(j) - partial%sa_drhototjk(3,np1)*ofct
                        endif
                        if (lopl) then
                          xdrv(l) = xdrv(l) + partial%sa_drhototjk(1,np1)*ofct
                          ydrv(l) = ydrv(l) + partial%sa_drhototjk(2,np1)*ofct
                          zdrv(l) = zdrv(l) + partial%sa_drhototjk(3,np1)*ofct
                        endif
                        if (nregionj.ne.nregionl) then
                          xregdrv(nregionj) = xregdrv(nregionj) - partial%sa_drhototjk(1,np1)*ofct
                          yregdrv(nregionj) = yregdrv(nregionj) - partial%sa_drhototjk(2,np1)*ofct
                          zregdrv(nregionj) = zregdrv(nregionj) - partial%sa_drhototjk(3,np1)*ofct
                          xregdrv(nregionl) = xregdrv(nregionl) + partial%sa_drhototjk(1,np1)*ofct
                          yregdrv(nregionl) = yregdrv(nregionl) + partial%sa_drhototjk(2,np1)*ofct
                          zregdrv(nregionl) = zregdrv(nregionl) + partial%sa_drhototjk(3,np1)*ofct
                        endif
                        if (lstr) then
                          do kl = 1,nstrains
                            rstrl(kl) = rstrl(kl) + partial%sa_drhotots(kl,np1)*ofct
                          enddo
                          if (latomicstress) then
                            if (.not.lfreeze) then
                              do kl = 1,nstrains
                                atomicstress(kl,i) = atomicstress(kl,i) + sixth*partial%sa_drhotots(kl,np1)*ofct
                                atomicstress(kl,j) = atomicstress(kl,j) + sixth*partial%sa_drhotots(kl,np1)*ofct
                                atomicstress(kl,k) = atomicstress(kl,k) + sixth*partial%sa_drhotots(kl,np1)*ofct
                              enddo
                            else
                              do kl = 1,nstrains
                                atomicstress(kl,i) = atomicstress(kl,i) + third*partial%sa_drhotots(kl,np1)*ofct
                                atomicstress(kl,j) = atomicstress(kl,j) + third*partial%sa_drhotots(kl,np1)*ofct
                                atomicstress(kl,k) = atomicstress(kl,k) + third*partial%sa_drhotots(kl,np1)*ofct
                              enddo
                            endif
                          endif
                        endif
                      enddo
                    endif
                  endif
!
!  Scale density and derivatives by screening factor
!
                  drhoij(1:3,1:maxmeamcomponent) = Sij*drhoij(1:3,1:maxmeamcomponent)
                  drhoji(1:3,1:maxmeamcomponent) = Sij*drhoji(1:3,1:maxmeamcomponent)
                  if (lsg1) then
                    drhoijs(1:6,1:maxmeamcomponent) = Sij*drhoijs(1:6,1:maxmeamcomponent)
                    drhojis(1:6,1:maxmeamcomponent) = Sij*drhojis(1:6,1:maxmeamcomponent)
                  endif
                  if (lgrad2) then
                    drhoij2(1:6,1:maxmeamcomponent) = Sij*drhoij2(1:6,1:maxmeamcomponent)
                    drhoji2(1:6,1:maxmeamcomponent) = Sij*drhoji2(1:6,1:maxmeamcomponent)
                    if (lsg2) then
                      drhoij2s(1:21,1:maxmeamcomponent) = Sij*drhoij2s(1:21,1:maxmeamcomponent)
                      drhoji2s(1:21,1:maxmeamcomponent) = Sij*drhoji2s(1:21,1:maxmeamcomponent)
                      drhoij2m(1:6,1:3,1:maxmeamcomponent) = Sij*drhoij2m(1:6,1:3,1:maxmeamcomponent)
                      drhoji2m(1:6,1:3,1:maxmeamcomponent) = Sij*drhoji2m(1:6,1:3,1:maxmeamcomponent)
                    endif
                  endif
!
!  Compute total derivatives of MEAM density
!
                  call meamtotalrhoderv(neamspeci,scrho(1,i),rhoi,drhoij,drhototij,drhoijs,drhototijs, &
                                        drhoij2,drhototij2,drhoij2s,drhototij2s,drhoij2m,drhototij2m, &
                                        lstr,lgrad2)
                  call meamtotalrhoderv(neamspecj,scrho(1,j),rhoj,drhoji,drhototji,drhojis,drhototjis, &
                                        drhoji2,drhototji2,drhoji2s,drhototji2s,drhoji2m,drhototji2m, &
                                        lstr,lgrad2)
                endif
              enddo
            endif
          else
            do mp = 1,npots
              npot = npotl(mp)
              if (nptype(npot).eq.19) then
                if (r.gt.rpot2(npot).and.r.le.rpot(npot)) then
                  lvalidij = .true.
                  call eamrho(nati,ntypi,natj,ntypj,r,rpot(npot),xcd,ycd,zcd,rhoij,rhoji,drhototij,drhototji, &
                              drhototijs,drhototjis,drhototij2,drhototji2,drhototij2s,drhototji2s, &
                              drhototij2m,drhototji2m,drhototij3,drhototji3,1.0_dp,1.0_dp,lsg1,.true.,lgrad2,.false., &
                              twopot(1,npot))
                endif
              endif
            enddo
          endif
!
!  Combine derivative terms
!
          if (QMMMmode(ncf).gt.0) then
            if (nregiontypi.ne.1) then
              deriv(1:3) = deriv(1:3) + rscrhoi*drhototij(1:3)*ofct
              if (lstr) then
                derivs(1:nstrains) = derivs(1:nstrains) + rscrhoi*drhototijs(1:nstrains)*ofct
              endif
            endif
            if (nregiontypj.ne.1) then
              deriv(1:3) = deriv(1:3) + rscrhoj*drhototji(1:3)*ofct
              if (lstr) then
                derivs(1:nstrains) = derivs(1:nstrains) + rscrhoj*drhototjis(1:nstrains)*ofct
              endif
            endif
          else
            deriv(1:3) = deriv(1:3) + (rscrhoi*drhototij(1:3) + rscrhoj*drhototji(1:3))*ofct
            if (lstr) then
              derivs(1:nstrains) = derivs(1:nstrains) + (rscrhoi*drhototijs(1:nstrains) + rscrhoj*drhototjis(1:nstrains))*ofct
            endif
          endif
          if (lgrad2) then
            ind = 0
            do kl = 1,3
              do km = kl,3
                ind = ind + 1
                deriv2(ind) = deriv2(ind) + (rscrhoi*drhototij2(ind) + rscrhoj*drhototji2(ind))*ofct
                deriv2(ind) = deriv2(ind) + ocj*rscrhoi3*drhototij(kl)*drhototij(km)*ofct
                deriv2(ind) = deriv2(ind) + oci*rscrhoj3*drhototji(kl)*drhototji(km)*ofct
              enddo
            enddo
            if (lstr) then
              deriv2s(1:nstrains2) = deriv2s(1:nstrains2) + (rscrhoi*drhototij2s(1:nstrains2) &
                                                          + rscrhoj*drhototji2s(1:nstrains2))*ofct
              deriv2m(1:nstrains,1:3) = deriv2m(1:nstrains,1:3) + (rscrhoi*drhototij2m(1:nstrains,1:3) &
                                                                + rscrhoj*drhototji2m(1:nstrains,1:3))*ofct
              ind = 0
              do kl = 1,nstrains
                do km = 1,kl
                  ind = ind + 1
                  deriv2s(ind) = deriv2s(ind) + ocj*rscrhoi3*ofct*drhototijs(kl)*drhototijs(km)
                  deriv2s(ind) = deriv2s(ind) + oci*rscrhoj3*ofct*drhototjis(kl)*drhototjis(km)
                enddo
                deriv2m(kl,1) = deriv2m(kl,1) + ocj*rscrhoi3*ofct*drhototijs(kl)*drhototij(1)
                deriv2m(kl,2) = deriv2m(kl,2) + ocj*rscrhoi3*ofct*drhototijs(kl)*drhototij(2)
                deriv2m(kl,3) = deriv2m(kl,3) + ocj*rscrhoi3*ofct*drhototijs(kl)*drhototij(3)
                deriv2m(kl,1) = deriv2m(kl,1) + oci*rscrhoj3*ofct*drhototjis(kl)*drhototji(1)
                deriv2m(kl,2) = deriv2m(kl,2) + oci*rscrhoj3*ofct*drhototjis(kl)*drhototji(2)
                deriv2m(kl,3) = deriv2m(kl,3) + oci*rscrhoj3*ofct*drhototjis(kl)*drhototji(3)
              enddo
            endif
          endif
        endif
        if (lvalidij) then
!
!  Only need to do internal derivatives if i not equals j
!
          if (i.ne.j) then
!******************************
!  Internal first derivatives *
!******************************
            if (lopi) then
              xdrv(i) = xdrv(i) - deriv(1)
              ydrv(i) = ydrv(i) - deriv(2)
              zdrv(i) = zdrv(i) - deriv(3)
            endif
            if (lopj) then
              xdrv(j) = xdrv(j) + deriv(1)
              ydrv(j) = ydrv(j) + deriv(2)
              zdrv(j) = zdrv(j) + deriv(3)
            endif
            if (nregioni.ne.nregionj) then
              xregdrv(nregioni) = xregdrv(nregioni) - deriv(1)
              yregdrv(nregioni) = yregdrv(nregioni) - deriv(2)
              zregdrv(nregioni) = zregdrv(nregioni) - deriv(3)
              xregdrv(nregionj) = xregdrv(nregionj) + deriv(1)
              yregdrv(nregionj) = yregdrv(nregionj) + deriv(2)
              zregdrv(nregionj) = zregdrv(nregionj) + deriv(3)
            endif
!********************************
!  Internal second derivatives  *
!********************************
            if (lgrad2) then
              derv2(jx,ix) = derv2(jx,ix) - deriv2(1)
              derv2(jy,ix) = derv2(jy,ix) - deriv2(2)
              derv2(jz,ix) = derv2(jz,ix) - deriv2(3)
              derv2(jx,iy) = derv2(jx,iy) - deriv2(2)
              derv2(jy,iy) = derv2(jy,iy) - deriv2(4)
              derv2(jz,iy) = derv2(jz,iy) - deriv2(5)
              derv2(jx,iz) = derv2(jx,iz) - deriv2(3)
              derv2(jy,iz) = derv2(jy,iz) - deriv2(5)
              derv2(jz,iz) = derv2(jz,iz) - deriv2(6)
            endif
!*****************************
!  Mixed strain derivatives  *
!*****************************
            if (lsg2) then
              do kl = 1,nstrains
                derv3(ix,kl) = derv3(ix,kl) - deriv2m(kl,1)
                derv3(iy,kl) = derv3(iy,kl) - deriv2m(kl,2)
                derv3(iz,kl) = derv3(iz,kl) - deriv2m(kl,3)
                derv3(jx,kl) = derv3(jx,kl) + deriv2m(kl,1)
                derv3(jy,kl) = derv3(jy,kl) + deriv2m(kl,2)
                derv3(jz,kl) = derv3(jz,kl) + deriv2m(kl,3)
              enddo
            endif
          endif
!*****************************
!  Strain first derivatives  *
!*****************************
          if (lsg1.or.lgrad2) then
            if (i.eq.j) then
              oij = 1.0_dp
            else
              oij = 2.0_dp
            endif
            do kl = 1,nstrains
              rstrl(kl)  = rstrl(kl)  + derivs(kl)*oij
            enddo
            if (latomicstress) then
              if (.not.lfreeze) then
                do kl = 1,nstrains
                  atomicstress(kl,i) = atomicstress(kl,i) + 0.25_dp*derivs(kl)*oij
                  atomicstress(kl,j) = atomicstress(kl,j) + 0.25_dp*derivs(kl)*oij
                enddo
              else
                do kl = 1,nstrains
                  atomicstress(kl,i) = atomicstress(kl,i) + 0.5_dp*derivs(kl)*oij
                  atomicstress(kl,j) = atomicstress(kl,j) + 0.5_dp*derivs(kl)*oij
                enddo
              endif
            endif
!******************************
!  Strain second derivatives  *
!******************************
            if (lsg2) then
              do kl = 1,nstrains
                do km = 1,nstrains
                  if (km.ge.kl) then
                    ind = km*(km - 1)/2 + kl
                  else
                    ind = kl*(kl - 1)/2 + km
                  endif
                  elcom(km,kl) = elcom(km,kl) + deriv2s(ind)*oij
                enddo
              enddo
            endif
          endif
        endif
        if (lgrad2) then
!******************************************************************
!  Start of third atom loop - only needed for second derivatives  *
!******************************************************************
          kloop: do k = 1,numat
            natk = nat(k)
            ntypk = nftype(k)
            xcrdik = xclat(k) - xal
            ycrdik = yclat(k) - yal
            zcrdik = zclat(k) - zal
            xcrdjk = xcrdik  - xcd
            ycrdjk = ycrdik - ycd
            zcrdjk = zcrdik - zcd
            ock = occuf(k)
            ofctijk = oci*ocj*ock
!
!  Check whether there are any potentials between i-k or j-k
!
            npotik = 0
            npotjk = 0
            rpik = 0.0_dp
            rpjk = 0.0_dp
            do n = 1,npote
              if (nptype(n).eq.19) then
                lvalidik = .false.
                if (nati.eq.nspec1(n).and.natk.eq.nspec2(n)) then
                  if (ntypi.eq.nptyp1(n).or.nptyp1(n).eq.0) then
                    if (ntypk.eq.nptyp2(n).or.nptyp2(n).eq.0) lvalidik = .true.
                  endif
                elseif (nati.eq.nspec2(n).and.natk.eq.nspec1(n)) then
                  if (ntypi.eq.nptyp2(n).or.nptyp2(n).eq.0) then
                    if (ntypk.eq.nptyp1(n).or.nptyp1(n).eq.0) lvalidik = .true.
                  endif
                endif
                if (lvalidik) then
                  npotik = npotik + 1
                  npotikptr(npotik) = n
                  if (rpot(n).gt.rpik) rpik = rpot(n)
                endif
!
                lvalidjk = .false.
                if (natj.eq.nspec1(n).and.natk.eq.nspec2(n)) then
                  if (ntypj.eq.nptyp1(n).or.nptyp1(n).eq.0) then
                    if (ntypk.eq.nptyp2(n).or.nptyp2(n).eq.0) lvalidjk = .true.
                  endif
                elseif (natj.eq.nspec2(n).and.natk.eq.nspec1(n)) then
                  if (ntypj.eq.nptyp2(n).or.nptyp2(n).eq.0) then
                    if (ntypk.eq.nptyp1(n).or.nptyp1(n).eq.0) lvalidjk = .true.
                  endif
                endif
                if (lvalidjk) then
                  npotjk = npotjk + 1
                  npotjkptr(npotjk) = n
                  if (rpot(n).gt.rpjk) rpjk = rpot(n)
                endif
              endif
            enddo
            rpijk = max(rpik,rpjk)
!
!  If no valid potentials for i-k or j-k then skip
!
            if ((npotik+npotjk).eq.0) cycle kloop
            cut2rk = rpijk*rpijk
            if (cut2rk.gt.cut2p) cut2rk = cut2p
            cut2k = cut2Xscale*cut2rk
!     
!  Find EAM species for k
!  
            neamspeck = neamfnspecptr(k)
!
!  Find index for i-k
!
            if (neamspeci.gt.neamspeck) then
              indik = neamspeci*(neamspeci - 1)/2 + neamspeck
            else
              indik = neamspeck*(neamspeck - 1)/2 + neamspeci
            endif
!
!  Find index for j-k
!
            if (neamspecj.gt.neamspeck) then
              indjk = neamspecj*(neamspecj - 1)/2 + neamspeck
            else
              indjk = neamspeck*(neamspeck - 1)/2 + neamspecj
            endif
!
!  Evaluate functional derivatives
!
            if (lMEAMfn) then
              frho(1:maxmeamcomponent) = scrho(1:maxmeamcomponent,k)
              call meamfnderv(neamfn,neamspeck,frho,rhok,eeam,rscrhok,rscrhok3,rscrhok5,.true.,lgrad2,.false.)
            else
              call eamfnderv(neamfn,neamspeck,scrho(1,k),eeam,rscrhok,rscrhok3,rscrhok5,.true.,lgrad2,.false.)
              rhok = scrho(1,k)
            endif
!
!  If no rho then skip
!
            if (abs(rhoi+rhoj+rhok).eq.0.0_dp) cycle kloop
!***********************
!  General cell loops  *
!***********************
            call rfindeither(xcrdik,ycrdik,zcrdik,xcrdjk,ycrdjk,zcrdjk,cut2k,cut2k,.false.,nor,nor2)
!
            drhototiksum(1:3) = 0.0_dp
            drhototjksum(1:3) = 0.0_dp
            drhototijk2sum(1:3,1:3) = 0.0_dp
            drhototjik2sum(1:3,1:3) = 0.0_dp
            if (lstr) then
              drhototikssum(1:nstrains) = 0.0_dp
              drhototjkssum(1:nstrains) = 0.0_dp
              drhototijk2ssum(1:nstrains,1:nstrains) = 0.0_dp
              drhototjik2ssum(1:nstrains,1:nstrains) = 0.0_dp
              drhototijk2msum(1:nstrains,1:3) = 0.0_dp
              drhototjik2msum(1:nstrains,1:3) = 0.0_dp
            endif
!
            do jj = 1,nor2
              rik2 = dist(nor+jj)
              rjk2 = dist2(nor+jj)
              xcd1 = xtmp(nor+jj)
              ycd1 = ytmp(nor+jj)
              zcd1 = ztmp(nor+jj)
              xcd2 = xtmp2(nor+jj)
              ycd2 = ytmp2(nor+jj)
              zcd2 = ztmp2(nor+jj)
!************************************************************
!  Calculate triangular contribution to second derivatives  *
!************************************************************
              lanyvalidik = .false.
              lanyvalidjk = .false.
              if (rik2.le.cut2k.and.npotik.gt.0) then
!*********************
!  i-k contribution  *
!*********************
!
!  Zero terms
!
                if (lMEAM) then
                  rhoik(1:maxmeamcomponent) = 0.0_dp
                  rhoki(1:maxmeamcomponent) = 0.0_dp
                  drhoik(1:3,1:maxmeamcomponent) = 0.0_dp
                  drhoki(1:3,1:maxmeamcomponent) = 0.0_dp
                  if (lstr) then
                    drhoiks(1:nstrains,1:maxmeamcomponent) = 0.0_dp
                    drhokis(1:nstrains,1:maxmeamcomponent) = 0.0_dp
                  endif
                else
                  rhoik(1) = 0.0_dp
                  rhoki(1) = 0.0_dp
                endif
                drhototik(1:3) = 0.0_dp
                drhototki(1:3) = 0.0_dp
                drhototijk2(1:3,1:3) = 0.0_dp
                if (lstr) then
                  drhototiks(1:nstrains) = 0.0_dp
                  drhototkis(1:nstrains) = 0.0_dp
                  drhototijk2s(1:nstrains,1:nstrains) = 0.0_dp
                  drhototijk2m(1:nstrains,1:3) = 0.0_dp
                endif
!
                rik = sqrt(rik2)
!
!  Loop over potentials to find many-body ones
!
                do n = 1,npotik
                  mp = npotikptr(n)
                  lvalidik = .false.
                  if (rik.gt.rpot2(mp).and.rik.le.rpot(mp)) then
                    lvalidik = .true.
                    if (lvalidik) then
!**********************************
!  Calculate density derivatives  *
!**********************************
                      if (lMEAMden) then
                        call meamrho(nati,ntypi,natk,ntypk,rik,rpot(mp),xcd1,ycd1,zcd1,rhoik,rhoki,drhoik,drhoki, &
                                     drhoiks,drhokis,drhoik2,drhoki2,drhoik2s,drhoki2s,drhoik2m,drhoki2m, &
                                     1.0_dp,1.0_dp,.true.,lstr,.true.,.false.,twopot(1,mp))
!*******************************
!  Compute screening function  *
!*******************************
!
!  Initialise screening function to 1
!
                        lnonzeroSik = .true.
                        Sik = 1.0_dp
!
                        if (lanyMEAMscreen) then
!
!  Set cutoffs for possible screening atoms
!
                          rcut2 = rcutfactor*rik2
!
!  Loop over atoms to search for images that may contribute to the screening
!
                          l = 0
                          npartialik = 0
                          do while (l.lt.numat.and.lnonzeroSik)
                            l = l + 1
                            neamspecl = neamfnspecptr(l)
                            if (lMEAMscreen(indik,neamspecl)) then
!
!  Set basic vectors between atoms
!
                              xil0 = xclat(l) - xal
                              yil0 = yclat(l) - yal
                              zil0 = zclat(l) - zal
                              xkl0 = xil0 - xcd1
                              ykl0 = yil0 - ycd1
                              zkl0 = zil0 - zcd1
!
!  Find images within cutoffs of both atoms - excluding self images
!
                              nvec0 = 0
                              call rfindmid(xil0,yil0,zil0,xkl0,ykl0,zkl0,rcut2,.false.,nvec0,nvec,vectorpair)
!
!  Loop over results of search
!
                              lvec = 0
                              do while (lvec.lt.nvec.and.lnonzeroSij)
                                lvec = lvec + 1
!
!  Compute screening function
!
                                call meamscreen(neamspeci,neamspeck,neamspecl,rik2,vectorpair%distance_pair1(lvec), &
                                                vectorpair%distance_pair2(lvec),Silk,dSilkdr,lpartial,.true.)
!
!  If screening function contribution is 0, then no need to continue for this pair
!
                                if (Silk.eq.0.0_dp) then
                                  lnonzeroSik = .false.
                                  Sik = 0.0_dp
                                else
!
!  Multiply total screening product
!
                                  Sik = Sik*Silk
                                  if (lpartial) then
!
!  If this atom has a screening factor between 0 and 1, we need to keep track of it since it will generate non-zero derivatives
!
                                    npartialik = npartialik + 1
                                    if (npartialik.gt.partialik%sa_maxdim) then
                                      call changemaxsa(partialik,npartialik)
                                    endif
                                    partialik%sa_atom(npartialik) = l
                                    partialik%sa_xik(npartialik) = vectorpair%xvector_pair1(lvec)
                                    partialik%sa_yik(npartialik) = vectorpair%yvector_pair1(lvec)
                                    partialik%sa_zik(npartialik) = vectorpair%zvector_pair1(lvec)
                                    partialik%sa_xjk(npartialik) = vectorpair%xvector_pair2(lvec)
                                    partialik%sa_yjk(npartialik) = vectorpair%yvector_pair2(lvec)
                                    partialik%sa_zjk(npartialik) = vectorpair%zvector_pair2(lvec)
                                    partialik%sa_Sikj(npartialik) = Silk
                                    partialik%sa_dSikjdr(1:3,npartialik) = dSilkdr(1:3)
                                  endif
                                endif
!
!  End loop over images of possible screening atoms
!
                              enddo
                            endif
!
!  End loop over possible screening atoms
!
                          enddo
!
!  End of screening function
!
                        endif
!
!  Only do remainder of work if the screening factor is non-zero
!
                        if (lnonzeroSik) then
!
!  Scale density and derivatives by screening factor
!
                          drhoik(1:3,1:maxmeamcomponent) = Sik*drhoik(1:3,1:maxmeamcomponent)
                          drhoki(1:3,1:maxmeamcomponent) = Sik*drhoki(1:3,1:maxmeamcomponent)
                          if (lsg1) then
                            drhoiks(1:6,1:maxmeamcomponent) = Sik*drhoiks(1:6,1:maxmeamcomponent)
                            drhokis(1:6,1:maxmeamcomponent) = Sik*drhokis(1:6,1:maxmeamcomponent)
                          endif
!
                          call meamtotalrhoderv(neamspeci,scrho(1,i),rhoi,drhoik,drhototik,drhoiks,drhototiks, &
                                                drhoik2,drhototik2,drhoik2s,drhototik2s,drhoik2m,drhototik2m, &
                                                lstr,.false.)
                          call meamtotalrhoderv(neamspeck,scrho(1,k),rhok,drhoki,drhototki,drhokis,drhototkis, &
                                                drhoki2,drhototki2,drhoki2s,drhototki2s,drhoki2m,drhototki2m, &
                                                lstr,.false.)
                          call meamtotalrhocrossderv(neamspeci,scrho(1,i),rhoi,drhoij,drhototij,drhoik,drhototik, &
                                                     drhoijs,drhoiks,drhototijk2,drhototijk2s,drhototijk2m,lstr)
                        endif
                      else
                        call eamrho(nati,ntypi,natk,ntypk,rik,rpot(mp),xcd1,ycd1,zcd1,rhoik,rhoki,drhototik,drhototki, &
                                    drhototiks,drhototkis,drhototik2,drhototki2,drhototik2s,drhototki2s, &
                                    drhototik2m,drhototki2m,drhototik3,drhototki3,1.0_dp,1.0_dp,lsg1,.true.,.false.,.false., &
                                    twopot(1,mp))
                      endif
                      lanyvalidik = .true.
                    endif
                  endif
                enddo
              endif
              if (rjk2.le.cut2.and.npotjk.gt.0) then
!*********************
!  j-k contribution  *
!*********************
!
!  Zero terms
!
                if (lMEAM) then
                  rhojk(1:maxmeamcomponent) = 0.0_dp
                  rhokj(1:maxmeamcomponent) = 0.0_dp
                  drhojk(1:3,1:maxmeamcomponent) = 0.0_dp
                  drhokj(1:3,1:maxmeamcomponent) = 0.0_dp
                  if (lstr) then
                    drhojks(1:nstrains,1:maxmeamcomponent) = 0.0_dp
                    drhokjs(1:nstrains,1:maxmeamcomponent) = 0.0_dp
                  endif
                else
                  rhojk(1) = 0.0_dp
                  rhokj(1) = 0.0_dp
                endif
                drhototjk(1:3) = 0.0_dp
                drhototkj(1:3) = 0.0_dp
                drhototjik2(1:3,1:3) = 0.0_dp
                if (lstr) then
                  drhototjks(1:nstrains) = 0.0_dp
                  drhototkjs(1:nstrains) = 0.0_dp
                  drhototjik2s(1:nstrains,1:nstrains) = 0.0_dp
                  drhototjik2m(1:nstrains,1:3) = 0.0_dp
                endif
!
                rjk = sqrt(rjk2)
!
!  Loop over potentials to find many-body ones
!
                do n = 1,npotjk
                  mp = npotjkptr(n)
                  lvalidjk = .false.
                  if (rjk.gt.rpot2(mp).and.rjk.le.rpot(mp)) then
                    lvalidjk = .true.
                    if (lvalidjk) then
!**********************************
!  Calculate density derivatives  *
!**********************************
                      if (lMEAMden) then
                        call meamrho(natj,ntypj,natk,ntypk,rjk,rpot(mp),xcd2,ycd2,zcd2,rhojk,rhokj,drhojk,drhokj, &
                                     drhojks,drhokjs,drhojk2,drhokj2,drhojk2s,drhokj2s,drhojk2m,drhokj2m, &
                                     1.0_dp,1.0_dp,.true.,lstr,.true.,.false.,twopot(1,mp))
!*******************************
!  Compute screening function  *
!*******************************
!
!  Initialise screening function to 1
!
                        lnonzeroSjk = .true.
                        Sjk = 1.0_dp
!
                        if (lanyMEAMscreen) then
!
!  Set cutoffs for possible screening atoms
!
                          rcut2 = rcutfactor*rjk2
!
!  Loop over atoms to search for images that may contribute to the screening
!
                          l = 0
                          npartialjk = 0
                          do while (l.lt.numat.and.lnonzeroSjk)
                            l = l + 1
                            neamspecl = neamfnspecptr(l)
                            if (lMEAMscreen(indjk,neamspecl)) then
!
!  Set basic vectors between atoms
!
                              xjl0 = xclat(l) - xcd + xal
                              yjl0 = yclat(l) - ycd + yal
                              zjl0 = zclat(l) - zcd + zal
                              xkl0 = xjl0 - xcd2
                              ykl0 = yjl0 - ycd2
                              zkl0 = zjl0 - zcd2
!
!  Find images within cutoffs of both atoms - excluding self images
!
                              nvec0 = 0
                              call rfindmid(xjl0,yjl0,zjl0,xkl0,ykl0,zkl0,rcut2,.false.,nvec0,nvec,vectorpair)
!
!  Loop over results of search
!
                              lvec = 0
                              do while (lvec.lt.nvec.and.lnonzeroSjk)
                                lvec = lvec + 1
!
!  Compute screening function
!
                                call meamscreen(neamspecj,neamspeck,neamspecl,rjk2,vectorpair%distance_pair1(lvec), &
                                                vectorpair%distance_pair2(lvec),Sjlk,dSjlkdr,lpartial,.true.)
!
!  If screening function contribution is 0, then no need to continue for this pair
!
                                if (Sjlk.eq.0.0_dp) then
                                  lnonzeroSjk = .false.
                                  Sjk = 0.0_dp
                                else
!
!  Multiply total screening product
!
                                  Sjk = Sjk*Sjlk
                                  if (lpartial) then
!
!  If this atom has a screening factor between 0 and 1, we need to keep track of it since it will generate non-zero derivatives
!
                                    npartialjk = npartialjk + 1
                                    if (npartialjk.gt.partialjk%sa_maxdim) then
                                      call changemaxsa(partialjk,npartialjk)
                                    endif
                                    partialjk%sa_atom(npartialjk) = l
                                    partialjk%sa_xik(npartialjk) = vectorpair%xvector_pair1(lvec)
                                    partialjk%sa_yik(npartialjk) = vectorpair%yvector_pair1(lvec)
                                    partialjk%sa_zik(npartialjk) = vectorpair%zvector_pair1(lvec)
                                    partialjk%sa_xjk(npartialjk) = vectorpair%xvector_pair2(lvec)
                                    partialjk%sa_yjk(npartialjk) = vectorpair%yvector_pair2(lvec)
                                    partialjk%sa_zjk(npartialjk) = vectorpair%zvector_pair2(lvec)
                                    partialjk%sa_Sikj(npartialjk) = Sjlk
                                    partialjk%sa_dSikjdr(1:3,npartialjk) = dSjlkdr(1:3)
                                  endif
                                endif
!
!  End loop over images of possible screening atoms
!
                              enddo
                            endif
!
!  End loop over possible screening atoms
!
                          enddo
!
!  End of screening function
!
                        endif
!
!  Only do remainder of work if the screening factor is non-zero
!
                        if (lnonzeroSjk) then
!
!  Scale density and derivatives by screening factor
!
                          drhojk(1:3,1:maxmeamcomponent) = Sjk*drhojk(1:3,1:maxmeamcomponent)
                          drhokj(1:3,1:maxmeamcomponent) = Sjk*drhokj(1:3,1:maxmeamcomponent)
                          if (lsg1) then
                            drhojks(1:6,1:maxmeamcomponent) = Sjk*drhojks(1:6,1:maxmeamcomponent)
                            drhokjs(1:6,1:maxmeamcomponent) = Sjk*drhokjs(1:6,1:maxmeamcomponent)
                          endif
                          call meamtotalrhoderv(neamspecj,scrho(1,j),rhoj,drhojk,drhototjk,drhojks,drhototjks, &
                                                drhojk2,drhototjk2,drhojk2s,drhototjk2s,drhojk2m,drhototjk2m, &
                                                lstr,.false.)
                          call meamtotalrhoderv(neamspeck,scrho(1,k),rhok,drhokj,drhototkj,drhokjs,drhototkjs, &
                                                drhokj2,drhototkj2,drhokj2s,drhototkj2s,drhokj2m,drhototkj2m, &
                                                lstr,.false.)
                          call meamtotalrhocrossderv(neamspecj,scrho(1,j),rhoj,drhoji,drhototji,drhojk,drhototjk, &
                                                     drhojis,drhojks,drhototjik2,drhototjik2s,drhototjik2m,lstr)
                        endif
                      else
                        call eamrho(natj,ntypj,natk,ntypk,rjk,rpot(mp),xcd2,ycd2,zcd2,rhojk,rhokj,drhototjk,drhototkj, &
                                    drhototjks,drhototkjs,drhototjk2,drhototkj2,drhototjk2s,drhototkj2s, &
                                    drhototjk2m,drhototkj2m,drhototjk3,drhototkj3,1.0_dp,1.0_dp,lsg1,.true.,.false.,.false., &
                                    twopot(1,mp))
                      endif
                      lanyvalidjk = .true.
                    endif
                  endif
                enddo
              endif
!
!  Cross term derivative for k-i / k-j
!
              if (lanyvalidik.and.lanyvalidjk.and.lMEAMden) then
                drhototkij2(1:3,1:3) = 0.0_dp
                if (lstr) then
                  drhototkij2s(1:nstrains,1:nstrains) = 0.0_dp
                  drhototkij2m(1:nstrains,1:3) = 0.0_dp
                endif
                call meamtotalrhocrossderv(neamspeck,scrho(1,k),rhok,drhoki,drhototki,drhokj,drhototkj,drhokis,drhokjs, &
                                           drhototkij2,drhototkij2s,drhototkij2m,lstr)
              endif
!**********************
!  Add terms to sums  *
!**********************
              if (lanyvalidik) then
                drhototiksum(1:3) = drhototiksum(1:3) + drhototik(1:3)
                drhototijk2sum(1:3,1:3) = drhototijk2sum(1:3,1:3) + drhototijk2(1:3,1:3)
                if (lsg2) then
                  drhototikssum(1:nstrains) = drhototikssum(1:nstrains) + drhototiks(1:nstrains)
                  drhototijk2ssum(1:nstrains,1:nstrains) = drhototijk2ssum(1:nstrains,1:nstrains) + &
                    drhototijk2s(1:nstrains,1:nstrains)
                  drhototijk2msum(1:nstrains,1:3) = drhototijk2msum(1:nstrains,1:3) + &
                    drhototijk2m(1:nstrains,1:3)
                endif
              endif
              if (lanyvalidjk) then
                drhototjksum(1:3) = drhototjksum(1:3) + drhototjk(1:3)
                drhototjik2sum(1:3,1:3) = drhototjik2sum(1:3,1:3) + drhototjik2(1:3,1:3)
                if (lsg2) then
                  drhototjkssum(1:nstrains) = drhototjkssum(1:nstrains) + drhototjks(1:nstrains)
                  drhototjik2ssum(1:nstrains,1:nstrains) = drhototjik2ssum(1:nstrains,1:nstrains) + &
                    drhototjik2s(1:nstrains,1:nstrains)
                  drhototjik2msum(1:nstrains,1:3) = drhototjik2msum(1:nstrains,1:3) + &
                    drhototjik2m(1:nstrains,1:3)
                endif
              endif
!
!  Only need to do internal derivatives if i not equals j
!
              if (i.ne.j) then
!****************************************************************************************************************
!  Calculate second derivatives for i-k/j-k terms - terms that have to be summed for each distance combination  *
!****************************************************************************************************************
!
!  i-k/j-k
!
                if (lanyvalidik.and.lanyvalidjk) then
                  if (lMEAM) then
                    dt1 = rscrhok3*ofctijk
                    dt2 = rscrhok*ofctijk
                    derv2(jx,ix) = derv2(jx,ix) + dt1*drhototki(1)*drhototkj(1) + dt2*drhototkij2(1,1)
                    derv2(jy,ix) = derv2(jy,ix) + dt1*drhototki(1)*drhototkj(2) + dt2*drhototkij2(1,2)
                    derv2(jz,ix) = derv2(jz,ix) + dt1*drhototki(1)*drhototkj(3) + dt2*drhototkij2(1,3)
                    derv2(jx,iy) = derv2(jx,iy) + dt1*drhototki(2)*drhototkj(1) + dt2*drhototkij2(2,1)
                    derv2(jy,iy) = derv2(jy,iy) + dt1*drhototki(2)*drhototkj(2) + dt2*drhototkij2(2,2)
                    derv2(jz,iy) = derv2(jz,iy) + dt1*drhototki(2)*drhototkj(3) + dt2*drhototkij2(2,3)
                    derv2(jx,iz) = derv2(jx,iz) + dt1*drhototki(3)*drhototkj(1) + dt2*drhototkij2(3,1)
                    derv2(jy,iz) = derv2(jy,iz) + dt1*drhototki(3)*drhototkj(2) + dt2*drhototkij2(3,2)
                    derv2(jz,iz) = derv2(jz,iz) + dt1*drhototki(3)*drhototkj(3) + dt2*drhototkij2(3,3)
                  else
                    dt1 = rscrhok3*ofctijk
                    derv2(jx,ix) = derv2(jx,ix) + dt1*drhototki(1)*drhototkj(1)
                    derv2(jy,ix) = derv2(jy,ix) + dt1*drhototki(1)*drhototkj(2)
                    derv2(jz,ix) = derv2(jz,ix) + dt1*drhototki(1)*drhototkj(3)
                    derv2(jx,iy) = derv2(jx,iy) + dt1*drhototki(2)*drhototkj(1)
                    derv2(jy,iy) = derv2(jy,iy) + dt1*drhototki(2)*drhototkj(2)
                    derv2(jz,iy) = derv2(jz,iy) + dt1*drhototki(2)*drhototkj(3)
                    derv2(jx,iz) = derv2(jx,iz) + dt1*drhototki(3)*drhototkj(1)
                    derv2(jy,iz) = derv2(jy,iz) + dt1*drhototki(3)*drhototkj(2)
                    derv2(jz,iz) = derv2(jz,iz) + dt1*drhototki(3)*drhototkj(3)
                  endif
                endif
              endif
!*************************************************
!  End of second derivative loop over distances  *
!*************************************************
            enddo
!
!  Only need to do internal derivatives if i not equals j
!
            if (i.ne.j) then
!*****************************************************************************************
!  Calculate second derivatives for i-k/j-k terms that can be summed over all distances  *
!*****************************************************************************************
!
!  i-k
!
              if (lMEAM) then
                dt1 = rscrhoi3*ofctijk
                dt2 = rscrhoi*ofctijk
                derv2(jx,ix) = derv2(jx,ix) - dt1*drhototiksum(1)*drhototij(1) - dt2*drhototijk2sum(1,1)
                derv2(jy,ix) = derv2(jy,ix) - dt1*drhototiksum(1)*drhototij(2) - dt2*drhototijk2sum(2,1)
                derv2(jz,ix) = derv2(jz,ix) - dt1*drhototiksum(1)*drhototij(3) - dt2*drhototijk2sum(3,1)
                derv2(jx,iy) = derv2(jx,iy) - dt1*drhototiksum(2)*drhototij(1) - dt2*drhototijk2sum(1,2)
                derv2(jy,iy) = derv2(jy,iy) - dt1*drhototiksum(2)*drhototij(2) - dt2*drhototijk2sum(2,2)
                derv2(jz,iy) = derv2(jz,iy) - dt1*drhototiksum(2)*drhototij(3) - dt2*drhototijk2sum(3,2)
                derv2(jx,iz) = derv2(jx,iz) - dt1*drhototiksum(3)*drhototij(1) - dt2*drhototijk2sum(1,3)
                derv2(jy,iz) = derv2(jy,iz) - dt1*drhototiksum(3)*drhototij(2) - dt2*drhototijk2sum(2,3)
                derv2(jz,iz) = derv2(jz,iz) - dt1*drhototiksum(3)*drhototij(3) - dt2*drhototijk2sum(3,3)
              else
                dt1 = rscrhoi3*ofctijk
                derv2(jx,ix) = derv2(jx,ix) - dt1*drhototiksum(1)*drhototij(1)
                derv2(jy,ix) = derv2(jy,ix) - dt1*drhototiksum(1)*drhototij(2)
                derv2(jz,ix) = derv2(jz,ix) - dt1*drhototiksum(1)*drhototij(3)
                derv2(jx,iy) = derv2(jx,iy) - dt1*drhototiksum(2)*drhototij(1)
                derv2(jy,iy) = derv2(jy,iy) - dt1*drhototiksum(2)*drhototij(2)
                derv2(jz,iy) = derv2(jz,iy) - dt1*drhototiksum(2)*drhototij(3)
                derv2(jx,iz) = derv2(jx,iz) - dt1*drhototiksum(3)*drhototij(1)
                derv2(jy,iz) = derv2(jy,iz) - dt1*drhototiksum(3)*drhototij(2)
                derv2(jz,iz) = derv2(jz,iz) - dt1*drhototiksum(3)*drhototij(3)
              endif
!
!  j-k
!
              if (lMEAM) then
                dt1 = rscrhoj3*ofctijk
                dt2 = rscrhoj*ofctijk
                derv2(jx,ix) = derv2(jx,ix) + dt1*drhototjksum(1)*drhototji(1) + dt2*drhototjik2sum(1,1)
                derv2(jy,ix) = derv2(jy,ix) + dt1*drhototjksum(2)*drhototji(1) + dt2*drhototjik2sum(1,2)
                derv2(jz,ix) = derv2(jz,ix) + dt1*drhototjksum(3)*drhototji(1) + dt2*drhototjik2sum(1,3)
                derv2(jx,iy) = derv2(jx,iy) + dt1*drhototjksum(1)*drhototji(2) + dt2*drhototjik2sum(2,1)
                derv2(jy,iy) = derv2(jy,iy) + dt1*drhototjksum(2)*drhototji(2) + dt2*drhototjik2sum(2,2)
                derv2(jz,iy) = derv2(jz,iy) + dt1*drhototjksum(3)*drhototji(2) + dt2*drhototjik2sum(2,3)
                derv2(jx,iz) = derv2(jx,iz) + dt1*drhototjksum(1)*drhototji(3) + dt2*drhototjik2sum(3,1)
                derv2(jy,iz) = derv2(jy,iz) + dt1*drhototjksum(2)*drhototji(3) + dt2*drhototjik2sum(3,2)
                derv2(jz,iz) = derv2(jz,iz) + dt1*drhototjksum(3)*drhototji(3) + dt2*drhototjik2sum(3,3)
              else
                dt1 = rscrhoj3*ofctijk
                derv2(jx,ix) = derv2(jx,ix) + dt1*drhototjksum(1)*drhototji(1)
                derv2(jy,ix) = derv2(jy,ix) + dt1*drhototjksum(2)*drhototji(1)
                derv2(jz,ix) = derv2(jz,ix) + dt1*drhototjksum(3)*drhototji(1)
                derv2(jx,iy) = derv2(jx,iy) + dt1*drhototjksum(1)*drhototji(2)
                derv2(jy,iy) = derv2(jy,iy) + dt1*drhototjksum(2)*drhototji(2)
                derv2(jz,iy) = derv2(jz,iy) + dt1*drhototjksum(3)*drhototji(2)
                derv2(jx,iz) = derv2(jx,iz) + dt1*drhototjksum(1)*drhototji(3)
                derv2(jy,iz) = derv2(jy,iz) + dt1*drhototjksum(2)*drhototji(3)
                derv2(jz,iz) = derv2(jz,iz) + dt1*drhototjksum(3)*drhototji(3)
              endif
            endif
            if (lsg2) then
!
!  Mixed internal-strain derivatives
!
              if (lMEAM) then
                dt1 = rscrhoi3*ofctijk
                dt2 = rscrhoj3*ofctijk
                dt3 = rscrhoi*ofctijk
                dt4 = rscrhoj*ofctijk
                do kl = 1,nstrains
                  derv3(ix,kl) = derv3(ix,kl) - dt1*drhototikssum(kl)*drhototij(1) - dt2*drhototjkssum(kl)*drhototji(1)
                  derv3(iy,kl) = derv3(iy,kl) - dt1*drhototikssum(kl)*drhototij(2) - dt2*drhototjkssum(kl)*drhototji(2)
                  derv3(iz,kl) = derv3(iz,kl) - dt1*drhototikssum(kl)*drhototij(3) - dt2*drhototjkssum(kl)*drhototji(3)
                  derv3(jx,kl) = derv3(jx,kl) + dt1*drhototikssum(kl)*drhototij(1) + dt2*drhototjkssum(kl)*drhototji(1)
                  derv3(jy,kl) = derv3(jy,kl) + dt1*drhototikssum(kl)*drhototij(2) + dt2*drhototjkssum(kl)*drhototji(2)
                  derv3(jz,kl) = derv3(jz,kl) + dt1*drhototikssum(kl)*drhototij(3) + dt2*drhototjkssum(kl)*drhototji(3)
!
                  derv3(ix,kl) = derv3(ix,kl) - dt3*drhototijk2msum(kl,1) - dt4*drhototjik2msum(kl,1)
                  derv3(iy,kl) = derv3(iy,kl) - dt3*drhototijk2msum(kl,2) - dt4*drhototjik2msum(kl,2)
                  derv3(iz,kl) = derv3(iz,kl) - dt3*drhototijk2msum(kl,3) - dt4*drhototjik2msum(kl,3)
                  derv3(jx,kl) = derv3(jx,kl) + dt3*drhototijk2msum(kl,1) + dt4*drhototjik2msum(kl,1)
                  derv3(jy,kl) = derv3(jy,kl) + dt3*drhototijk2msum(kl,2) + dt4*drhototjik2msum(kl,2)
                  derv3(jz,kl) = derv3(jz,kl) + dt3*drhototijk2msum(kl,3) + dt4*drhototjik2msum(kl,3)
                enddo
              else
                dt1 = rscrhoi3*ofctijk
                dt2 = rscrhoj3*ofctijk
                do kl = 1,nstrains
                  derv3(ix,kl) = derv3(ix,kl) - dt1*drhototikssum(kl)*drhototij(1) - dt2*drhototjkssum(kl)*drhototji(1)
                  derv3(iy,kl) = derv3(iy,kl) - dt1*drhototikssum(kl)*drhototij(2) - dt2*drhototjkssum(kl)*drhototji(2)
                  derv3(iz,kl) = derv3(iz,kl) - dt1*drhototikssum(kl)*drhototij(3) - dt2*drhototjkssum(kl)*drhototji(3)
                  derv3(jx,kl) = derv3(jx,kl) + dt1*drhototikssum(kl)*drhototij(1) + dt2*drhototjkssum(kl)*drhototji(1)
                  derv3(jy,kl) = derv3(jy,kl) + dt1*drhototikssum(kl)*drhototij(2) + dt2*drhototjkssum(kl)*drhototji(2)
                  derv3(jz,kl) = derv3(jz,kl) + dt1*drhototikssum(kl)*drhototij(3) + dt2*drhototjkssum(kl)*drhototji(3)
                enddo
              endif
!
!  Strain-strain second derivatives
!
              dt1l = oij*dt1
              if (lMEAM) then
                dt3l = oij*dt3
                do kl = 1,nstrains
                  do km = 1,nstrains
                    elcom(km,kl) = elcom(km,kl) + dt1l*drhototikssum(kl)*drhototijs(km)
                    elcom(km,kl) = elcom(km,kl) + dt3l*drhototijk2ssum(km,kl)
                  enddo
                enddo
              else
                do kl = 1,nstrains
                  do km = 1,nstrains
                    elcom(km,kl) = elcom(km,kl) + dt1l*drhototikssum(kl)*drhototijs(km)
                  enddo
                enddo
              endif
              dt2l = oij*dt2
              if (lMEAM) then
                dt4l = oij*dt4
                do kl = 1,nstrains
                  do km = 1,nstrains
                    elcom(km,kl) = elcom(km,kl) + dt2l*drhototjkssum(kl)*drhototjis(km)
                    elcom(km,kl) = elcom(km,kl) + dt4l*drhototjik2ssum(km,kl)
                  enddo
                enddo
              else
                do kl = 1,nstrains
                  do km = 1,nstrains
                    elcom(km,kl) = elcom(km,kl) + dt2l*drhototjkssum(kl)*drhototjis(km)
                  enddo
                enddo
              endif
            endif
!***************************
!  End of third atom loop  *
!***************************
          enddo kloop
        endif
!**************************************
!  End of valid distance i-j section  *
!**************************************
      enddo
    enddo jloop
  enddo iloop
!  
!  Double counting correction
!  
  if (.not.lfreeze) then
    if (lsg1) then
      do i = 1,nstrains
        rstrd(i) = rstrd(i) + 0.5_dp*rstrl(i)
      enddo
    endif
    if (lsg2) then
      do i = 1,nstrains
        do j = 1,i
          sderv2(i,j) = sderv2(i,j) + 0.5_dp*elcom(i,j)
        enddo
      enddo
    endif
  else
    if (lsg1) then
      do i = 1,nstrains
        rstrd(i) = rstrd(i) + rstrl(i)
      enddo
    endif
    if (lsg2) then
      do i = 1,nstrains
        do j = 1,i
          sderv2(i,j) = sderv2(i,j) + elcom(i,j)
        enddo
      enddo
    endif
  endif
!
!  End of real space part - perform general tasks
!
999 continue
!
!  Free local memory
!
  if (lgrad2) then
    deallocate(npotjkptr,stat=status)
    if (status/=0) call deallocate_error('many','npotjkptr')
    deallocate(npotikptr,stat=status)
    if (status/=0) call deallocate_error('many','npotikptr')
  endif
  deallocate(npotl,stat=status)
  if (status/=0) call deallocate_error('many','npotl')
!
!  Exit point
!
1000 continue
!
!  Unscale density
!
  call eamscalescrho(-1_i4)
!
!  Timing
!
  time2 = g_cpu_time()
  tmany = tmany + time2 - time1
#ifdef TRACE
  call trace_out('many')
#endif
!
  return
  end
