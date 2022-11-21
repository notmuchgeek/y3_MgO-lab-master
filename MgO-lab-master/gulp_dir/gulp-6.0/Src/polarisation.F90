  subroutine polarisation(epolar,esregion12,esregion2,eattach,lsymnotused,lgrad1,lgrad2)
!
!  Main subroutine for calculating the point-ion polarisability
!  contribution to the lattice energy.
!
!   5/00 Created from energy.f
!   2/01 Logic for call of derivatives routines corrected
!   4/01 Region 2 self energy calculation for surfaces added
!  11/01 Attachment energy added
!   8/02 Surface energy calculation algorithm changed
!  11/02 Parallel modifications made
!   6/09 Site energies added
!   2/18 Trace added
!   4/19 Second derivatives added
!   7/19 Corrections to symmetry handling
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!  10/19 Langevin damping of dipoles added
!   2/20 Correction to Langevin damped formalism
!   7/20 Separate routine for sumall with 1 argument added
!
!  Note : Needs correction for polarisation energy between
!         surface region 1 and surface region 2.
!
!  If lsymnotused is true then vx/vy/vz were generated without use
!  of symmetry and derivatives should be generated likewise.
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
  use configurations, only : nregions, nregionno
  use control,        only : lseok
  use g_constants
  use current
  use energies,       only : siteenergy
  use parallel
  use polarise
  use potentialxyz
  use symmetry
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  logical,                      intent(in)    :: lgrad1
  logical,                      intent(in)    :: lgrad2
  logical,                      intent(in)    :: lsymnotused
  real(dp),                     intent(inout) :: epolar
  real(dp),                     intent(inout) :: esregion12
  real(dp),                     intent(inout) :: esregion2
  real(dp),                     intent(inout) :: eattach
!
!  Local variables
!
  integer(i4)                                 :: i
  integer(i4)                                 :: status
  real(dp)                                    :: eattachl
  real(dp)                                    :: epolar2
  real(dp)                                    :: epolarsum
  real(dp)                                    :: Es
  real(dp)                                    :: g_cpu_time
  real(dp)                                    :: Ep
  real(dp)                                    :: dEpdEs
  real(dp)                                    :: d2EpdEs2
  real(dp)                                    :: tsum0
  real(dp)                                    :: vxi
  real(dp)                                    :: vyi
  real(dp)                                    :: vzi
  real(dp), dimension(:),   allocatable       :: sum
#ifdef TRACE
  call trace_in('polarisation')
#endif
  if (nprocs.gt.1) then
!************************************
!  Perform global sums on vx,vy,vz  *
!************************************
    allocate(sum(numat),stat=status)
    if (status/=0) call outofmemory('polarisation','sum')
!
    tsum0 = g_cpu_time()
!
    call sumall(vx,sum,numat,"polarisation","vx")
    vx(1:numat) = sum(1:numat)
!
    call sumall(vy,sum,numat,"polarisation","vy")
    vy(1:numat) = sum(1:numat)
!
    call sumall(vz,sum,numat,"polarisation","vz")
    vz(1:numat) = sum(1:numat)
!
    call sumall(vx12,sum,numat,"polarisation","vx12")
    vx12(1:numat) = sum(1:numat)
!
    call sumall(vy12,sum,numat,"polarisation","vy12")
    vy12(1:numat) = sum(1:numat)
!
    call sumall(vz12,sum,numat,"polarisation","vz12")
    vz12(1:numat) = sum(1:numat)
!
    tsum = tsum + g_cpu_time() - tsum0
!
    deallocate(sum,stat=status)
    if (status/=0) call deallocate_error('polarisation','sum')
  endif
!*********************
!  Calculate energy  *
!*********************
  if (lsymopt.and.lsymderv.and..not.lsymnotused) then
    if (lpollangevin) then
      do i = procid+1,nasym,nprocs
        if (dpolar(i).gt.1.0d-12) then
          vxi = vx(i)
          vyi = vy(i)
          vzi = vz(i)
          Es = sqrt(vxi*vxi + vyi*vyi + vzi*vzi)
          call langevinpol(Es,dpolarmax(i),dpolar(i),Ep,dEpdEs,d2EpdEs2,.false.,.false.)
          epolar = epolar - Ep*occua(i)*dble(neqv(i))
          siteenergy(i) = siteenergy(i) - 0.5_dp*Ep*occua(i)*dble(neqv(i))/angstoev
        endif
      enddo
    else
      do i = procid+1,nasym,nprocs
        if (dpolar(i).gt.1.0d-12) then
          vxi = vx(i)
          vyi = vy(i)
          vzi = vz(i)
          epolar = epolar - dpolar(i)*(vxi*vxi + vyi*vyi + vzi*vzi)*occua(i)*dble(neqv(i))
          siteenergy(i) = siteenergy(i) - 0.5_dp*dpolar(i)*(vxi*vxi + vyi*vyi + vzi*vzi)*occua(i)*dble(neqv(i))/angstoev
        endif
      enddo
    endif
  else
    if (lpollangevin) then
      do i = procid+1,numat,nprocs
        if (dpolar(nrelf2a(i)).gt.1.0d-12) then
          vxi = vx(i)
          vyi = vy(i)
          vzi = vz(i)
          Es = sqrt(vxi*vxi + vyi*vyi + vzi*vzi)
          call langevinpol(Es,dpolarmax(nrelf2a(i)),dpolar(nrelf2a(i)),Ep,dEpdEs,d2EpdEs2,.false.,.false.)
          epolar = epolar - Ep*occuf(i)
          siteenergy(i) = siteenergy(i) - 0.5_dp*Ep*occuf(i)/angstoev
        endif
      enddo
    else
      do i = procid+1,numat,nprocs
        if (dpolar(nrelf2a(i)).gt.1.0d-12) then
          vxi = vx(i)
          vyi = vy(i)
          vzi = vz(i)
          epolar = epolar - dpolar(nrelf2a(i))*(vxi*vxi + vyi*vyi + vzi*vzi)*occuf(i)
          siteenergy(i) = siteenergy(i) - 0.5_dp*dpolar(nrelf2a(i))*(vxi*vxi + vyi*vyi + vzi*vzi)*occuf(i)/angstoev
        endif
      enddo
    endif
  endif
  epolar = 0.5_dp*epolar/angstoev
  if (lseok.and.nregions(ncf).gt.1) then
!
!  If there are multiple regions then symmetry isn't being used
!
    epolar2 = 0.0_dp
    if (lpollangevin) then
      do i = procid+1,numat,nprocs
        if (nregionno(nsft+i).gt.1.and.dpolar(i).gt.1.0d-12) then
          vxi = vx(i)
          vyi = vy(i)
          vzi = vz(i)
          Es = sqrt(vxi*vxi + vyi*vyi + vzi*vzi)
          call langevinpol(Es,dpolarmax(i),dpolar(i),Ep,dEpdEs,d2EpdEs2,.false.,.false.)
          epolar2 = epolar2 - Ep*occuf(i)
        endif
      enddo
    else
      do i = procid+1,numat,nprocs
        if (nregionno(nsft+i).gt.1.and.dpolar(i).gt.1.0d-12) then
          vxi = vx(i)
          vyi = vy(i)
          vzi = vz(i)
          epolar2 = epolar2 - dpolar(i)*(vxi*vxi + vyi*vyi + vzi*vzi)*occuf(i)
        endif
      enddo
    endif
    esregion2 = esregion2 + 0.5_dp*epolar2/angstoev
!
!  In new algorithm region2 self contribution must be excluded from epolar
!
    epolar = epolar - 0.5_dp*epolar2/angstoev
  endif
  eattachl = 0.0_dp
  if (lpollangevin) then
    do i = procid+1,numat,nprocs
      if (dpolar(i).gt.1.0d-12) then
        vxi = vx12(i)
        vyi = vy12(i)
        vzi = vz12(i)
        Es = sqrt(vxi*vxi + vyi*vyi + vzi*vzi)
        call langevinpol(Es,dpolarmax(i),dpolar(i),Ep,dEpdEs,d2EpdEs2,.false.,.false.)
        eattachl = eattachl - Ep*occuf(i)
      endif
    enddo
  else
    do i = procid+1,numat,nprocs
      if (dpolar(i).gt.1.0d-12) then
        eattachl = eattachl - dpolar(i)*(vx12(i)*vx12(i) + vy12(i)*vy12(i) + vz12(i)*vz12(i))*occuf(i)
      endif
    enddo
  endif
  eattachl = 0.5_dp*eattachl/angstoev
  eattach = eattach + eattachl
!
  if (nprocs.gt.1) then
!
!  If in parallel check that sum of polarisation energies are zero
!
    call sumone(epolar,epolarsum,"polarisation","epolar")
  else
    epolarsum = epolar
  endif
!
!  If derivatives are required and polarisation is significant then continue further, else return
!
  if (abs(epolarsum).lt.1.0d-12.or..not.lgrad1) then
#ifdef TRACE
    call trace_out('polarisation')
#endif
    return
  endif
!*******************************
!  Electrostatic contribution  *
!*******************************
  if (ndim.gt.0) then
    if (lsymopt.and.lsymderv.and..not.lsymnotused) then
      call pirealrecips
    else
      if (nprocs.gt.1) then
        call pirealrecipd(lgrad2)
      else
        call pirealrecip(lgrad2)
      endif
    endif
  else
    if (nprocs.gt.1) then
      call pireal0dd(lgrad2)
    else
      call pireal0d(lgrad2)
    endif
  endif
#ifdef TRACE
  call trace_out('polarisation')
#endif
!
  return
  end
