  subroutine potcut(np)
!
!  Calculates energy and/or gradient shifts to smooth
!  potentials at cutoff boundary.
!
!   2/04 Boundary correction for qerfc potential added
!  11/04 Sqrt pi constant used
!   4/05 Mods for cosh-spring potential added
!  11/06 qoverr2 potential added
!   3/07 Check on whether potential should have an energy/gradient 
!        shift added
!   8/07 Possibility of non-integer powers in L-J potentials added
!   1/09 Arrays apot, bpot, cpot & dpot replaced by new 2-D array twopot
!   1/09 leshift/lgshift logicals introduced to replace use of tpot(3/4,) for energy
!        and gradient shifts
!   3/14 derfc changed to g_derfc for benefit of ChemShell
!   2/15 MM3buck added
!   2/18 Trace added
!   2/18 Slater potential added
!   8/18 Energy shift added for Lennard-Jones epsilon-sigma
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
!  Copyright Curtin University 2018
!
!  Julian Gale, CIC, Curtin University, August 2018
!
  use g_constants,    only : sqrtpi
  use control
  use element
  use general,        only : nwarn
  use iochannels
  use numbers
  use parallel
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  use two
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in) :: np
!
!  Local variables
!
  integer(i4)             :: mpt
  integer(i4)             :: npt
  integer(i4)             :: nptyp
  logical                 :: lOKpot
  real(dp)                :: apt
  real(dp)                :: bpt
  real(dp)                :: br
  real(dp)                :: cpt
  real(dp)                :: g_derfc
  real(dp)                :: dpt
  real(dp)                :: deriv
  real(dp)                :: eatom
  real(dp)                :: r6
  real(dp)                :: r12
  real(dp)                :: rd
  real(dp)                :: rd2
  real(dp)                :: rk
  real(dp)                :: rk2
  real(dp)                :: rmpt
  real(dp)                :: rnpt
  real(dp)                :: rp
  real(dp)                :: sig
  real(dp)                :: trm1
  real(dp)                :: trm2
  real(dp)                :: trme
#ifdef TRACE
  call trace_in('potcut')
#endif
!
  rp = rpot(np)
  rk = 1.0_dp/rp
  rk2 = rk*rk
  nptyp = nptype(np)
  lOKpot = .false.
  if (nptyp.eq.1) then
! 
!  Buckingham potential
!
    apt = twopot(1,np)
    bpt = 1.0_dp/twopot(2,np)
    cpt = twopot(3,np)
    r6 = rk2**3
    trm1 = apt*exp(-bpt*rp)
    eatom = trm1 - cpt*r6
    deriv = 6.0_dp*cpt*r6*rk - trm1*bpt
    lOKpot = .true.
  elseif (nptyp.eq.2) then
!
!  Lennard-Jones potential
!
    apt = twopot(1,np)
    bpt = twopot(2,np)
    rmpt = tpot(1,np)
    rnpt = tpot(2,np)
    r6 = rk**rnpt
    r12 = rk**rmpt
    eatom = (apt*r12 - bpt*r6)
    deriv = (rnpt*bpt*r6 - rmpt*apt*r12)*rk
    lOKpot = .true.
  elseif (nptyp.eq.3) then
!
!  Morse potential - no coulomb offset
!
    apt = twopot(1,np)
    bpt = twopot(2,np)
    cpt = twopot(3,np)
    trme = exp(-bpt*(rp-cpt))
    trm1 = 1.0_dp - trme
    trm2 = 2.0_dp*apt*bpt*trme
    deriv = trm2*trm1
    trm1 = trm1*trm1
    eatom = apt*(trm1 - 1.0_dp)
    lOKpot = .true.
  elseif (nptyp.eq.4) then
!
!  Morse potential - coulomb offset
!
    nwarn = nwarn + 1
    if (ioproc) then
      call outwarning('Energy/gradient shifts not allowed for Morse-q pot',0_i4)
    endif
  elseif (nptyp.eq.5) then
!
!  Harmonic potential - no coulomb offset
!
    apt = twopot(1,np)
    bpt = twopot(2,np)
    cpt = twopot(3,np)
    dpt = twopot(4,np)
    rd = rp - bpt
    rd2 = rd*rd
    eatom = 0.5_dp*apt*rd2
    deriv = apt*rd
    if (cpt.ne.0.0_dp) then
      eatom = eatom + sixth*cpt*rd2*rd
      deriv = deriv + 0.5_dp*rd2*cpt
    endif
    if (dpt.ne.0.0_dp) then
      eatom = eatom + 0.25_dp*sixth*dpt*rd2*rd2
      deriv = deriv + sixth*rd2*dpt
    endif
    lOKpot = .true.
  elseif (nptyp.eq.6) then
!
!  Harmonic potential - coulomb offset
!
    nwarn = nwarn + 1
    if (ioproc) then
      call outwarning('Energy/gradient shifts not allowed for Harmonic-q pot',0_i4)
    endif
  elseif (nptyp.eq.7) then
!
!  General potential
!
    apt = twopot(1,np)
    bpt = twopot(2,np)
    if (bpt.ne.0.0_dp) bpt = 1.0_dp/bpt
    cpt = twopot(3,np)
    mpt = nint(tpot(1,np))
    npt = nint(tpot(2,np))
    r6 = rk**npt
    r12 = rk**mpt
    trm1 = apt*exp(-bpt*rp)
    eatom = trm1*r12 - cpt*r6
    deriv = npt*cpt*r6*rk - trm1*bpt*r12
    deriv = deriv - trm1*mpt*rk*r12
    lOKpot = .true.
  elseif (nptyp.eq.8) then
!
!  Spring potential
!
    nwarn = nwarn + 1
    if (ioproc) then
      call outwarning('Energy/gradient shifts not allowed for Spring potential',0_i4)
    endif
  elseif (nptyp.eq.9) then
!
!  Coulomb potential
!
    nwarn = nwarn + 1
    if (ioproc) then
      call outwarning('Energy/gradient shifts not allowed for Coulomb potential',0_i4)
    endif
  elseif (nptyp.eq.12.or.nptyp.eq.13.or.nptyp.eq.59) then
!
!  Lennard-Jones potential - epsilon/sigma form
!
    apt = twopot(1,np)*twopot(3,np)
    bpt = twopot(1,np)*twopot(4,np)
    sig = twopot(2,np)
    rmpt = tpot(1,np)
    rnpt = tpot(2,np)
    r6 = bpt*(sig*rk)**rnpt
    r12 = apt*(sig*rk)**rmpt
    eatom = (r12 - r6)
    r6 = rnpt*r6
    r12 = rmpt*r12
    deriv = (r6-r12)*rk2
    lOKpot = .true.
  elseif (nptyp.eq.16) then
!
!  Inverse gaussian
!
    apt = twopot(1,np)
    bpt = twopot(2,np)
    cpt = twopot(3,np)
    trm1 = bpt*(rp-cpt)
    trme = apt*exp(-trm1*(rp-cpt))
    eatom = - trme
    deriv = 2.0_dp*trm1*trme
    lOKpot = .true.
  elseif (nptyp.eq.18) then
!
!  Damped dispersion potential
!
    nwarn = nwarn + 1
    if (ioproc) then
      call outwarning('Energy/gradient shifts not allowed for damped dispersion',0_i4)
    endif
  elseif (nptyp.eq.24) then
!
!  Qerfc potential
!
    apt = 1.0_dp/twopot(1,np)
    eatom = g_derfc(rp*apt)*rk
    trm1 = 2.0_dp*apt*exp(-rp*rp*apt*apt)/sqrtpi
    deriv = - (eatom + trm1)*rk
    lOKpot = .true.
  elseif (nptyp.eq.28) then
!
!  Squared harmonic potential - no coulomb offset
!
    apt = twopot(1,np)
    bpt = twopot(2,np)
    rd = rp*rp - bpt*bpt
    rd2 = rd*rd
    eatom = 0.25_dp*apt*rd2
    deriv = apt*rd*rp
    lOKpot = .true.
  elseif (nptyp.eq.29) then
!
!  Squared harmonic potential - coulomb offset
!
    nwarn = nwarn + 1
    if (ioproc) then
      call outwarning('Energy/gradient shifts not allowed for Harmonic-q pot',0_i4)
    endif
  elseif (nptyp.eq.33) then
!
!  Cosh-spring potential
!
    nwarn = nwarn + 1
    if (ioproc) then
      call outwarning('Energy/gradient shifts not allowed for Cosh-spring potential',0_i4)
    endif
  elseif (nptyp.eq.36) then
!
!  QoverR2 potential
!
    eatom = rk2
    deriv = - 2.0_dp*eatom*rk
    lOKpot = .true.
  elseif (nptyp.eq.57) then
! 
!  MM3 Buckingham potential
!
    apt = twopot(1,np)*twopot(4,np)
    bpt = twopot(2,np)/twopot(5,np)
    cpt = twopot(3,np)*twopot(4,np)*twopot(5,np)**6
    r6 = rk2**3
    trm1 = apt*exp(-bpt*rp)
    eatom = trm1 - cpt*r6
    deriv = 6.0_dp*cpt*r6*rk - trm1*bpt
    lOKpot = .true.
  elseif (nptyp.eq.61) then
!
!  Slater potential
!
    apt = twopot(1,np)
    bpt = twopot(2,np)
    br  = bpt*rp
    trm1 = apt*exp(-br)
    trm2 = (br*br/3.0_dp + br + 1.0_dp)
    eatom = trm1*trm2
    deriv = - trm1*bpt*br*(br + 1.0_dp)
    lOKpot = .true.
  endif
!
!  Energy shift
!
  if (leshift(np).and.lOKpot) then
    eshift(np) = - eatom
  endif
!
!  Gradient shift
!
  if (lgshift(np).and.lOKpot) then
    eshift(np) = eshift(np) + deriv*rp
    gshift(np) = - deriv
  endif
#ifdef TRACE
  call trace_out('potcut')
#endif
!
  return
  end
