  subroutine epot(lfrst,nsite,vsite,xsite,ysite,zsite,lgrad1,vx,vy,vz,lgrad2,efg,lshell)
!
!  Routine calculates the electrostatic potential at a set of general sites.
!
!  nsite = number of sites for potential calculation
!  vsite = calculated potential on output
!  xsite = x cartesian coordinate of site
!  ysite = y cartesian coordinate of site
!  zsite = z cartesian coordinate of site
!  lshell= if .true. then subtract core-shell terms
!
!   7/95 Electric field and EFG added
!  11/95 Bug fixed - now works with "noenergy" keyword / call to kindex
!        was needed
!   6/00 Dimensions of efg switch and msites argument removed
!   2/01 Modifications for 2-D periodic systems added
!   8/01 Modified to use call to qmatrixelement for terms. Now
!        handles the 0-D case as well.
!   1/03 Modifications made for Wolf sum
!   4/03 Sign of EFG changed since this second derivative block
!        refers to an on diagonal block, rather than off diagonal
!   7/03 Parallel modifications added
!  11/07 Use of the noelectro keyword now trapped so that quantities are
!        initialised and then routine exits.
!   3/09 small replaced by global value smallself from general module
!  10/09 Modified to allow for situation where nouterloop is too large to compute
!   5/12 Intent added for passed variables
!   2/18 Trace added
!  12/19 Algorithm for reciprocal space in 3D changed to increase the speed for large systems
!   7/20 Calls to qmatrixelement changed to qmatrixelemente which is a version that doesn't
!        duplicate the 3D reciprocal space code here
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
  use datatypes,     only : i4_limit
  use general,       only : smallself
  use kspace
  use parallel
  use shells
  use symmetry
  use times
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4),               intent(in)       :: nsite
  logical,                   intent(in)       :: lfrst
  logical,                   intent(in)       :: lgrad1
  logical,                   intent(in)       :: lgrad2
  logical,                   intent(in)       :: lshell
  real(dp),                  intent(out)      :: efg(6,*)
  real(dp),                  intent(out)      :: vsite(*)
  real(dp),                  intent(out)      :: vx(*)
  real(dp),                  intent(out)      :: vy(*)
  real(dp),                  intent(out)      :: vz(*)
  real(dp),                  intent(in)       :: xsite(*)
  real(dp),                  intent(in)       :: ysite(*)
  real(dp),                  intent(in)       :: zsite(*)
!
!  Local variables
!
  integer(i4)                                 :: i
  integer(i4)                                 :: ii
  integer(i4)                                 :: iv
  integer(i4)                                 :: j
  integer(i4)                                 :: nout
  integer(i4)                                 :: nouterloop
  integer(i4)                                 :: status
  real(dp)                                    :: g_cpu_time
  real(dp)                                    :: argk
  real(dp)                                    :: cosi
  real(dp)                                    :: sini
  real(dp)                                    :: qlj
  real(dp)                                    :: qme
  real(dp)                                    :: dqme(3)
  real(dp)                                    :: d2qme(6)
  real(dp)                                    :: small2
  real(dp),   dimension(:), allocatable       :: sum
  real(dp)                                    :: tsum0
  real(dp)                                    :: xd
  real(dp)                                    :: yd
  real(dp)                                    :: zd
#ifdef TRACE
  call trace_in('epot')
#endif
!
  do i = 1,nsite
    vsite(i) = 0.0_dp
  enddo
  if (lgrad1) then
    do i = 1,nsite
      vx(i) = 0.0_dp
      vy(i) = 0.0_dp
      vz(i) = 0.0_dp
    enddo
    if (lgrad2) then
      do i = 1,nsite
        do j = 1,6
          efg(j,i) = 0.0_dp
        enddo
      enddo
    endif
  endif
!
!  If electrostatics have been turned off then there is no point in proceeding.
!
  if (.not.lDoElectrostatics) then
#ifdef TRACE
    call trace_out('epot')
#endif
    return
  endif
!
  if (lshell) then
    small2 = cuts*cuts
  else
    small2 = smallself
  endif
  if (lewald.and.ndim.gt.1) then
    call kindex
  endif
!
!  Initialise terms
!
  call initqmatrix
  if (lfrst.and.lewald.and..not.lnorecip) then
    call initktrm
  endif
!***************************
!  Reciprocal space in 3D  *
!***************************
  if (lewald.and.ndim.eq.3.and..not.lnorecip) then
!
!  Start loop over cluster atom - unit cell atom pairs
!
    do iv = 1,nkvec
      csin(iv) = 0.0_dp
      sine(iv) = 0.0_dp
    enddo
    do j = procid+1,numat,nprocs
      qlj = qf(j)*occuf(j)
      do iv = 1,nkvec
        argk = xrk(iv)*xclat(j) + yrk(iv)*yclat(j) + zrk(iv)*zclat(j)
        csin(iv) = csin(iv) + qlj*cos(argk)
        sine(iv) = sine(iv) + qlj*sin(argk)
      enddo
    enddo
!
!  Global sums
!
    if (nprocs.gt.1) then
      tsum0 = g_cpu_time()
      allocate(sum(nkvec),stat=status)
      if (status/=0) call outofmemory('epot','sum')
!
      call sumall(csin,sum,nkvec,"psumall","csin")
      do i = 1,nkvec
        csin(i) = sum(i)
      enddo
!
      call sumall(sine,sum,nkvec,"psumall","sine")
      do i = 1,nkvec
        sine(i) = sum(i)
      enddo
!
      deallocate(sum,stat=status)
      if (status/=0) call deallocate_error('epot','sum')
      tsum = tsum + g_cpu_time() - tsum0
    endif
!
!  Add reciprocal space part of site potential
!
    do i = procid+1,nsite,nprocs
      do iv = 1,nkvec
        argk = xrk(iv)*xsite(i) + yrk(iv)*ysite(i) + zrk(iv)*zsite(i)
        cosi = cos(argk)
        sini = sin(argk)
        vsite(i) = vsite(i) + ktrm(iv)*(csin(iv)*cosi + sine(iv)*sini)
        if (lgrad1) then
          argk = ktrm(iv)*(-csin(iv)*sini + sine(iv)*cosi)
          vx(i) = vx(i) + argk*xrk(iv)
          vy(i) = vy(i) + argk*yrk(iv)
          vz(i) = vz(i) + argk*zrk(iv)
          if (lgrad2) then
            argk = - ktrm(iv)*(csin(iv)*cosi + sine(iv)*sini)
            efg(1,i) = efg(1,i) + argk*xrk(iv)*xrk(iv)
            efg(2,i) = efg(2,i) + argk*xrk(iv)*yrk(iv)
            efg(3,i) = efg(3,i) + argk*yrk(iv)*yrk(iv)
            efg(4,i) = efg(4,i) + argk*xrk(iv)*zrk(iv)
            efg(5,i) = efg(5,i) + argk*yrk(iv)*zrk(iv)
            efg(6,i) = efg(6,i) + argk*zrk(iv)*zrk(iv)
          endif
        endif
      enddo
    enddo
  endif
!*****************************
!  Loop over pairs of sites  *
!*****************************
  if (numat.gt.i4_limit) then
!
!  If nsite*numat is greater than integer precision will allow change algorithm
!
    do i = 1,nsite
      do j = procid+1,numat,nprocs
        xd = xclat(j) - xsite(i)
        yd = yclat(j) - ysite(i)
        zd = zclat(j) - zsite(i)
        qlj = qf(j)*occuf(j)
        call qmatrixelemente(xd,yd,zd,small2,lgrad1,lgrad2,qme,dqme,d2qme)
        vsite(i) = vsite(i) + qme*qlj
        if (lgrad1) then
          vx(i) = vx(i) + dqme(1)*qlj
          vy(i) = vy(i) + dqme(2)*qlj
          vz(i) = vz(i) + dqme(3)*qlj
          if (lgrad2) then
            efg(1,i) = efg(1,i) - d2qme(1)*qlj
            efg(2,i) = efg(2,i) - d2qme(2)*qlj
            efg(3,i) = efg(3,i) - d2qme(3)*qlj
            efg(4,i) = efg(4,i) - d2qme(4)*qlj
            efg(5,i) = efg(5,i) - d2qme(5)*qlj
            efg(6,i) = efg(6,i) - d2qme(6)*qlj
          endif
        endif
      enddo
    enddo
  else
    nouterloop = nsite*numat
    do nout = procid+1,nouterloop,nprocs
      i = ((nout - 1)/numat) + 1
      j = nout - (i - 1)*numat
      xd = xclat(j) - xsite(i)
      yd = yclat(j) - ysite(i)
      zd = zclat(j) - zsite(i)
      qlj = qf(j)*occuf(j)
      call qmatrixelemente(xd,yd,zd,small2,lgrad1,lgrad2,qme,dqme,d2qme)
      vsite(i) = vsite(i) + qme*qlj
      if (lgrad1) then
        vx(i) = vx(i) + dqme(1)*qlj
        vy(i) = vy(i) + dqme(2)*qlj
        vz(i) = vz(i) + dqme(3)*qlj
        if (lgrad2) then
          efg(1,i) = efg(1,i) - d2qme(1)*qlj
          efg(2,i) = efg(2,i) - d2qme(2)*qlj
          efg(3,i) = efg(3,i) - d2qme(3)*qlj
          efg(4,i) = efg(4,i) - d2qme(4)*qlj
          efg(5,i) = efg(5,i) - d2qme(5)*qlj
          efg(6,i) = efg(6,i) - d2qme(6)*qlj
        endif
      endif
    enddo
  endif
!
!  Global sums
!
  if (nprocs.gt.1) then
    tsum0 = g_cpu_time()
    if (lgrad2) then
      allocate(sum(6*nsite),stat=status)
    else
      allocate(sum(nsite),stat=status)
    endif
    if (status/=0) call outofmemory('epot','sum')
    call sumall(vsite,sum,nsite,"psumall","vsite")
    do i = 1,nsite
      vsite(i) = sum(i)
    enddo
    if (lgrad1) then
      call sumall(vx,sum,nsite,"psumall","vx")
      do i = 1,nsite
        vx(i) = sum(i)
      enddo
      call sumall(vy,sum,nsite,"psumall","vy")
      do i = 1,nsite
        vy(i) = sum(i)
      enddo
      call sumall(vz,sum,nsite,"psumall","vz")
      do i = 1,nsite
        vz(i) = sum(i)
      enddo
      if (lgrad2) then
        call sumall(efg,sum,6_i4*nsite,"psumall","efg")
        ii = 0
        do i = 1,nsite
          do j = 1,6
            ii = ii + 1
            efg(j,i) = sum(ii)
          enddo
        enddo
      endif
    endif
    deallocate(sum,stat=status)
    if (status/=0) call deallocate_error('epot','sum')
    tsum = tsum + g_cpu_time() - tsum0
  endif
!
!  Convert units to eV
!
  do i = 1,nsite
    vsite(i) = vsite(i)*angstoev
  enddo
  if (lgrad1) then
    do i = 1,nsite
      vx(i) = vx(i)*angstoev
      vy(i) = vy(i)*angstoev
      vz(i) = vz(i)*angstoev
    enddo
    if (lgrad2) then
      do i = 1,nsite
        do j = 1,6
          efg(j,i) = efg(j,i)*angstoev
        enddo
      enddo
    endif
  endif
#ifdef TRACE
  call trace_out('epot')
#endif
!
  return
  end
