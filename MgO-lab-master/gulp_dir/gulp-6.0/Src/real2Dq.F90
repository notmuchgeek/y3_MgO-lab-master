  subroutine real2Dq(ereal,erecip,qtot,lgrad1,lgrad2)
!
!  Calculates the real space contribution to the interaction of an ion with it's own 2-D images.
!
!   6/01 Created from reale
!   1/03 Wolf sum modifcations made
!   1/05 rp no longer sqrt'd and passed to rsearch routines
!   3/07 Printing of twobody energies added as an option
!   5/07 Argument list for twobody call modified
!  11/07 Modifications for reaxFF Coulomb term added
!  12/07 Unused variables removed
!   1/08 lreaxFFqreal removed
!   1/09 Integer datatypes all explicitly declared
!   2/09 Expanding and contracting of maxdis arrays removed
!   3/09 lorder12 added to twobody argument list
!   4/12 Explicit virial calculation removed as no longer needed
!   7/13 Call to twobody modified to include core-shell occupancy factor
!  12/14 rtrm1 changed from scalar to array
!   2/18 Trace added
!   8/18 Call to twostrterms introduced for setting of rpd
!   9/18 Strain module introduced
!   1/19 Finite strain second derivatives finished
!   8/19 Corrections to call to twobody
!  11/19 Rigid molecule modifications added
!   4/20 derv3c and d2r2dsdc added for benefit of rigid molecules
!   4/20 derv3c changes reversed as they are no longer required
!   6/20 Second derivatives restricted to ioproc to avoid duplication
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
!  Julian Gale, CIC, Curtin University, June 2020
!
  use g_constants
  use control
  use current
  use derivatives
  use eam,            only : maxmeamcomponent
  use element
  use general,        only : cutw
  use kspace
  use m_strain,       only : twostrterms
  use mdlogic
  use molecule
  use optimisation
  use parallel
  use realvectors
  use shells
  use symmetry
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  logical,     intent(in)    :: lgrad1
  logical,     intent(in)    :: lgrad2
  real(dp)                   :: ereal
  real(dp)                   :: erecip
  real(dp)                   :: qtot
!
!  Local variables
!
  integer(i4)                :: k
  integer(i4)                :: kk
  integer(i4)                :: kl
  integer(i4)                :: ks
  integer(i4)                :: kt
  integer(i4)                :: nmolonly
  integer(i4)                :: nor
  logical                    :: lmdl
  logical                    :: lself
  logical                    :: lsg1
  logical                    :: lsg2
  real(dp)                   :: g_cpu_time
  real(dp)                   :: cut2
  real(dp)                   :: cut2q
  real(dp)                   :: cut2s
  real(dp)                   :: eatom
  real(dp)                   :: ec6
  real(dp)                   :: factor
  real(dp)                   :: fct
  real(dp)                   :: ofct
  real(dp)                   :: ospfct
  real(dp)                   :: sctrm1(maxmeamcomponent)
  real(dp)                   :: sctrm2(maxmeamcomponent)
  real(dp)                   :: time1
  real(dp)                   :: time2
#ifdef TRACE
  call trace_in('real2Dq')
#endif
!
  time1 = g_cpu_time()
!
!  Local variables
!
  lmdl = lmd
  if (.not.lgrad1) lmdl = .false.
!
  lsg1 = (lgrad1.and.lstr)
  lsg2 = (lgrad2.and.lstr)
!
!  Set up cutoffs
!
  cut2s = cuts*cuts
  if (lwolf) then        
    cut2q = cutw*cutw
  else
    cut2q = rmx2
  endif
  cut2 = cut2q
!
  if (lnoreal) then
#ifdef TRACE   
    call trace_out('real2Dq')
#endif
    return
  endif
!***************************************************************
!  Atomistic and real space electrostatic component of energy  *
!***************************************************************
  ofct = 0.5_dp
  ospfct = ofct
  fct = ofct*angstoev
  factor = qtot*qtot*fct
!***********************
!  Find valid vectors  *
!***********************
  call rsearch2D(0.0_dp,0.0_dp,0.0_dp,.false.,.false.,1_i4,1_i4, &
                 0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2)
!
!  Self term
!
  erecip = erecip + factor*tweatpi
!
  if (nor.eq.0) goto 1110
!
!  Sqrt distances
!
  do k = 1,nor
    dist(k) = sqrt(dist(k))
  enddo
!*************************************************
!  Calculate twobody contribution in real space  *
!*************************************************
  call twobody(eatom,ereal,ec6,lgrad1,lgrad2,.false.,nor,1_i4,0_i4,0_i4,0.0_dp,cut2q,cut2s, &
               0_i4,factor,ofct,ospfct,0.0_dp,sctrm1,sctrm2,qtot,qtot,.false.,.true.,.false., &
               .false.,.false.)
!
!  Change sign of ereal
!
  ereal = - ereal
!
!  Generate products for derivatives
!
  if (lmdl.or.lsg1.or.lgrad2) then
    call twostrterms(ndim,maxdis,nor,xtmp,ytmp,ztmp,0.0_dp,0.0_dp,0.0_dp,dr2ds,rpd,d2r2dsdx,d2r2ds2,lgrad2)
  endif
!***********************
!  Strain derivatives  *
!***********************
!
!  First derivatives 
!
  if (lsg1.or.lgrad2) then
    do kl = 1,nstrains
      ks = nstrptr(kl)
      do k = 1,nor
        rstrd(kl) = rstrd(kl) - deriv(k)*dr2ds(k,ks)
      enddo
    enddo
!
!  Second derivatives
!
    if (lsg2.and.ioproc) then
      do kk = 1,nstrains
        ks = nstrptr(kk)
        do kl = 1,nstrains
          kt = nstrptr(kl)
          do k = 1,nor
            sderv2(kl,kk) = sderv2(kl,kk) - deriv2(k)*dr2ds(k,kt)*dr2ds(k,ks)
            sderv2(kl,kk) = sderv2(kl,kk) - deriv(k)*d2r2ds2(k,kt,ks)
          enddo
        enddo
      enddo
    endif
  endif
1110 continue
!
!  End of real space part - perform general tasks
!
!  Timing
!
  time2 = g_cpu_time()
  tatom = tatom + time2 - time1
#ifdef TRACE   
  call trace_out('real2Dq')
#endif
!
  return
  end
