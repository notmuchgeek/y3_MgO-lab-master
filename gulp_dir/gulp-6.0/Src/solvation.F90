  subroutine solvation(ecosmo,lgrad1,lgrad2)
!
!  Subroutine for calculating the solvation energy according to the COSMO model.
!
!  11/04 Use of sasparticle data structures added
!  12/04 Matrix cosmoBq introduced in energy calc
!  12/04 Computation of charges on SAS moved to new routine
!   1/05 Setting of deltaq removed since this is now done
!        in setqonsas
!   1/05 New COSMIC energy implemented in which deltaq is
!        multiplied by segment weight
!  12/08 Migrated to version 3.5 and converted to f90 format
!   4/17 Changes for parallelisation with distributed memory made
!   4/17 Energy calculation parallelised
!   1/18 Trace added
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
!  Julian Gale, CIC, Curtin University, January 2018
!
  use g_constants
  use cosmic
  use current
  use maths,       only : lcosmo2D
  use parallel,    only : node2pts, nptsonnode
  use times,       only : tcosmo, tcosmoderv
#ifdef TRACE
  use trace,       only : trace_in, trace_out
#endif
  implicit none
!
!  Passed arguments
!
  real(dp),    intent(out) :: ecosmo
  logical,     intent(in)  :: lgrad1
  logical,     intent(in)  :: lgrad2
!
!  Local variables
!
  integer(i4)              :: ipts
  integer(i4)              :: iptsloc
  integer(i4)              :: maxval
  real(dp)                 :: fact
  real(dp)                 :: t1
  real(dp)                 :: t2
  real(dp)                 :: t3
!
!  Functions
!
  real(dp)                 :: g_cpu_time
#ifdef TRACE
  call trace_in('solvation')
#endif
!
!  Initialise timer for this routine
!     
  t1 = g_cpu_time()                                                                              
!
!  Set conversion factor x dielectric factor
!
  fact = 0.5_dp*autoev*autoangs*cosmofneps
  maxnptsonnode = max(maxnptsonnode,nptsonnode)
!
!  Check that dimensions for two-dimensional arrays are up to date
!
  if (lcosmo2D) then
    if (maxnpts.gt.maxcosmoA) then
      maxcosmoA = maxnpts
      maxcosmoAu = max(maxcosmoAu,maxnptsonnode)
      call changemaxcosmoA
    elseif (maxnptsonnode.gt.maxcosmoAu) then
      maxcosmoAu = maxnptsonnode
      maxcosmoA = max(maxcosmoA,maxnpts)
      call changemaxcosmoA
    endif
  else
    maxval = maxnpts*(maxnpts+1_i4)/2_i4
    if (maxval.gt.maxcosmoA) then
      maxcosmoA = maxval
      maxcosmoAu = 1_i4
      call changemaxcosmoA
    endif
  endif
!
!  Initialise Coulomb matrix element terms
!
  call initqmatrix
  call initqmatrixc
  call initktrm
!
!  Generate segment weighting factors
!
  call setsegweight
!
!  Generate SAS - SAS matrix - cosmoA
!
  if (lcosmo2D) then
    call setcosmoamat2D
  else
    call setcosmoamat
  endif
!
!  Generate SAS - atom matrix - cosmoB
!
  call setcosmobmat
!
!  Calculate charges on SAS
!
  if (lcosmo2D) then
    call setqonsas2D(lgrad2)
  else
    call setqonsas(lgrad2)
  endif
!
!  Calculate energy
!
  ecosmo = 0.0_dp
  do iptsloc = 1,nptsonnode
    ipts = node2pts(iptsloc)
    ecosmo = ecosmo + (qonsas(ipts) - deltaq*segweight(ipts))*cosmoBq(ipts)
  enddo
  ecosmo = fact*ecosmo
!
  t2 = g_cpu_time()                                                                              
!
!  Calculate gradients
!
  if (lgrad1.or.lgrad2) then
    if (lcosmo2D) then
      call cosmoderv2D(lgrad2)
    else
      call cosmoderv(lgrad2)
    endif
  endif
!
!  Finalise timer for this routine
!     
  t3 = g_cpu_time()                                                                              
  tcosmo = tcosmo + t2 - t1
  tcosmoderv = tcosmoderv + t3 - t2
#ifdef TRACE
  call trace_out('solvation')
#endif
!
  return
  end
