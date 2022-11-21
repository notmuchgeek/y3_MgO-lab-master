  subroutine recip2Dq(erecip,qtot,lgrad1,lgrad2)
!
!  Calculates the energy of interaction of an ion with a
!  periodic array in 2-D of itself. Assumes that it is 
!  being called after recip2D so that setup phase is
!  unnecessary.
!
!   6/01 Created from recip2D
!  11/02 Parallel modifications made
!   6/07 Bug in parallel execution corrected
!  12/07 Unused variables removed
!   6/12 Handling of esum corrected 
!  12/16 derfc changed to g_derfc for benefit of ChemShell
!   2/18 Trace added
!   9/18 Handling of lstraincell algorithm added
!  11/18 Finite strain flag introduced instead of lstraincell
!   1/19 Handling of norecip case corrected to avoid memory issue
!   1/19 Error in setting of strderv fixed when lgrad2 is false
!   5/19 Modified so that strfin additions to sderv2 are no longer needed
!   6/20 Corrected for changes in distribution of k vectors
!   7/20 Separate routine for sumall with 1 argument added
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
  use g_constants
  use control
  use current
  use derivatives
  use kspace
  use m_strain,       only : gstrterms, strainddetds, straindet, straind2detds2
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
  logical,    intent(in)                        :: lgrad1
  logical,    intent(in)                        :: lgrad2
  real(dp),   intent(inout)                     :: erecip
  real(dp),   intent(in)                        :: qtot
!
!  Local variables
!
  integer(i4)                                   :: i
  integer(i4)                                   :: iv
  integer(i4)                                   :: j
  integer(i4)                                   :: kk
  integer(i4)                                   :: ll
  integer(i4)                                   :: nlocalkvec
  integer(i4)                                   :: nremainder
  integer(i4)                                   :: status
  logical                                       :: lsg1
  real(dp)                                      :: cosq
  real(dp)                                      :: g_cpu_time
  real(dp)                                      :: csinq
  real(dp)                                      :: darg1
  real(dp)                                      :: g_derfc
  real(dp)                                      :: derfc1
  real(dp)                                      :: dexp1
  real(dp)                                      :: eltrm
  real(dp)                                      :: eltrm1
  real(dp)                                      :: eltrm2
  real(dp)                                      :: erecipl
  real(dp)                                      :: esum
  real(dp)                                      :: gseta
  real(dp)                                      :: kexperfc
  real(dp)                                      :: kvec
  real(dp)                                      :: qfct
  real(dp)                                      :: rkvec
  real(dp)                                      :: strdrvl(3)
  real(dp)                                      :: strm1
  real(dp),   dimension(:),   allocatable       :: sum
  real(dp)                                      :: time0
  real(dp)                                      :: time1
  real(dp)                                      :: tsum0
  real(dp)                                      :: twoqv
!
  if (lnorecip) return
!
#ifdef TRACE
  call trace_in('recip2Dq')
#endif
  time0 = g_cpu_time()
  erecipl = 0.0_dp
  lsg1 = (lstr.and.lgrad1)
!
!  Distribute kvec loops 
!
  if (lgrad2.and.nprocs.gt.1) then
    nlocalkvec = nkvec
  else
    nlocalkvec = (nkvec/nprocs)
    nremainder = nkvec - nlocalkvec*nprocs
    if (procid.lt.nremainder) nlocalkvec = nlocalkvec + 1
  endif
!
!  Allocate local memory
!
  allocate(sum(max(1_i4,nstrains)),stat=status)
  if (status/=0) call outofmemory('recip2Dq','sum')
!
!  Define constants
!
  rpieta = 1.0_dp/sqrt(pi*eta)
  rhseta = 0.5_dp/seta
!
!  Build products of K vector components
!
  if (lsg1.or.lgrad2) then
    call gstrterms(ndim,maxkvec,nlocalkvec,xrk,yrk,zrk,dg2ds,d2g2dx2,d2g2ds2,lgrad2)
  endif
  qfct = 0.5_dp*qtot*qtot*angstoev
!
!  First term - K vector independent
!
  twoqv = qfct*vol4pi
  if (lgrad2.and.nprocs.gt.1) then
    erecipl = erecipl + twoqv*rpieta
  elseif (ioproc) then
    erecipl = erecipl + twoqv*rpieta
  endif
!
!  Second term - K vector dependent
!
  csinq = 0.0_dp
  if (lgrad1) then
    if (lsg1) then
      strdrvl(1:3) = 0.0_dp
    endif
    do iv = 1,nlocalkvec
      cosq = qfct*ktrm(iv)
      kvec = kmod(iv)
      darg1 = kvec*rhseta
      dexp1 = exp(-(darg1)**2)
      derfc1 = g_derfc(darg1)
      kexperfc = 2.0_dp*derfc1
!
!  Energy
!
      csinq = csinq + cosq*kexperfc
!
      if (lsg1) then
!
!  Strain first derivatives
!
        rkvec = 1.0_dp/kvec
        strm1 = rkvec*(-rkvec*kexperfc - 2.0*rpieta*dexp1)
        do kk = 1,nstrains
          strdrvl(kk) = strdrvl(kk) - strm1*cosq*dg2ds(iv,kk)
        enddo
        if (lgrad2.and.ioproc) then
          eltrm1 = - 2.0_dp*strm1*cosq*rkvec*rkvec
          eltrm2 = cosq*rkvec*rkvec*(2.0*derfc1*(rkvec*rkvec) + &
            rpieta*(2.0_dp*dexp1*(rkvec+0.5_dp*kvec/eta)))
          eltrm = (eltrm1 + eltrm2)
          do kk = 1,nstrains
            do ll = 1,nstrains
              sderv2(kk,ll) = sderv2(kk,ll) - eltrm*dg2ds(iv,kk)*dg2ds(iv,ll) - strm1*cosq*d2g2ds2(iv,ll,kk)
            enddo
          enddo
        endif
      endif
    enddo
  else
    do iv = 1,nlocalkvec
      kvec = kmod(iv)
      gseta = kvec*rhseta
      kexperfc = 2.0_dp*g_derfc(gseta)
      csinq = csinq + ktrm(iv)*kexperfc
    enddo
    csinq = csinq*qfct
  endif
!
!  Lattice energy
!
  erecipl = erecipl - csinq
!****************
!  Global sums  *
!****************
  if (.not.lgrad2.or.nprocs.eq.1) then
    call sumone(erecipl,esum,"recip2Dq","erecipl")
    erecipl = esum
  endif
!
  erecip = erecip + erecipl
!
  if (.not.lgrad2) then
    tsum0 = g_cpu_time()
    if (lsg1) then
      call sumall(strdrvl,sum,nstrains,"recip2Dq","strderv")
      do i = 1,nstrains
        strdrvl(i) = sum(i)
        strderv(i) = strderv(i) + sum(i)
      enddo
    endif
    tsum = tsum + g_cpu_time() - tsum0
  else
    do i = 1,nstrains
      strderv(i) = strderv(i) + strdrvl(i)
    enddo
  endif
!****************************************************
!  Complete strain derivatives in reciprocal space  *
!****************************************************
  if (lsg1) then
    esum = erecipl
    if (lgrad2.and.ioproc) then
!
!  Area corrections to strain second derivatives
!
      if (lfinitestrain) then
        do i = 1,3
          do j = 1,3
            sderv2(j,i) = sderv2(j,i) - strdrvl(j)*strainddetds(i)*straindet &
                                      - strdrvl(i)*strainddetds(j)*straindet &
                                      + 2.0_dp*esum*strainddetds(j)*strainddetds(i)*straindet**2 &
                                      - esum*straind2detds2(j,i)*straindet
          enddo
        enddo
      else
        do i = 1,2
          do j = 1,2
            sderv2(j,i) = sderv2(j,i) - strdrvl(j) - strdrvl(i)
            sderv2(j,i) = sderv2(j,i) + esum
          enddo
          sderv2(3,i) = sderv2(3,i) - strdrvl(3)
          sderv2(i,3) = sderv2(i,3) - strdrvl(3)
        enddo
      endif
    endif
    if (lfinitestrain) then
      strderv(1) = strderv(1) - strainddetds(1)*straindet*esum
      strderv(2) = strderv(2) - strainddetds(2)*straindet*esum
      strderv(3) = strderv(3) - strainddetds(3)*straindet*esum
    else
      strderv(1) = strderv(1) - esum
      strderv(2) = strderv(2) - esum
    endif
  endif
!
!  Free local memory
!
  deallocate(sum,stat=status)
  if (status/=0) call deallocate_error('recip2Dq','sum')
!
!  Timing
!
  time1 = g_cpu_time()
  tion = tion + time1 - time0
#ifdef TRACE
  call trace_out('recip2Dq')
#endif
!
  return
  end
