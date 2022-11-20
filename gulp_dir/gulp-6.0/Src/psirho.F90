  subroutine psirho(nati,ntypi,neamspeci,natj,ntypj,neamspecj,npot,r,rhoi,rhoj,drhoidr,drhojdr, &
                    d2rhoidr2,d2rhojdr2,lgrad1,lgrad2)
!
!  Subroutine for calculating rho(R) for use in MEAM twobody potential
!
!   8/14 Created from psibaskes
!   3/15 Trap for floating point error due to small densities added
!   2/18 Trace added
!  11/18 Bug in second derivatives fixed
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
!  Julian Gale, CIC, Curtin University, November 2018
!
  use eam
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  use two
  implicit none
!
!  Passed variables
!
  integer(i4),    intent(in)      :: nati
  integer(i4),    intent(in)      :: natj
  integer(i4),    intent(in)      :: neamspeci
  integer(i4),    intent(in)      :: neamspecj
  integer(i4),    intent(in)      :: ntypi
  integer(i4),    intent(in)      :: ntypj
  integer(i4),    intent(in)      :: npot
  logical,        intent(in)      :: lgrad1
  logical,        intent(in)      :: lgrad2
  real(dp),       intent(in)      :: r         ! Distance
  real(dp),       intent(out)     :: rhoi      ! Density for reference structure at i for this distance
  real(dp),       intent(out)     :: rhoj      ! Density for reference structure at j for this distance
  real(dp),       intent(out)     :: drhoidr   ! First derivative of density at i 
  real(dp),       intent(out)     :: drhojdr   ! First derivative of density at j 
  real(dp),       intent(out)     :: d2rhoidr2 ! Second derivative of density at i
  real(dp),       intent(out)     :: d2rhojdr2 ! Second derivative of density at j
!
!  Local variables
!
  integer(i4)                     :: i5
  integer(i4)                     :: i6
  integer(i4)                     :: i7
  integer(i4)                     :: i8
  integer(i4)                     :: l               ! Looping index over higher order MEAM components
  integer(i4)                     :: mi              ! Pointer to meam functional for i
  integer(i4)                     :: mj              ! Pointer to meam functional for j
  logical                         :: lfound
  logical                         :: lmatch
  real(dp)                        :: betaoverr0i
  real(dp)                        :: betaoverr0j
  real(dp)                        :: ctrm
  real(dp)                        :: exptrm
  real(dp)                        :: gamma
  real(dp)                        :: dgammadr
  real(dp)                        :: d2gammadr2
  real(dp)                        :: drhoi2dr
  real(dp)                        :: drhoj2dr
  real(dp)                        :: drhoi20dr
  real(dp)                        :: drhoj20dr
  real(dp)                        :: drhoi20sdr
  real(dp)                        :: drhoj20sdr
  real(dp)                        :: drhoi20_2NNdr
  real(dp)                        :: drhoj20_2NNdr
  real(dp)                        :: drhoi2sdr
  real(dp)                        :: drhoj2sdr
  real(dp)                        :: d2rhoi2dr2
  real(dp)                        :: d2rhoj2dr2
  real(dp)                        :: d2rhoi2sdr2
  real(dp)                        :: d2rhoj2sdr2
  real(dp)                        :: d2rhoi20dr2
  real(dp)                        :: d2rhoj20dr2
  real(dp)                        :: d2rhoi20sdr2
  real(dp)                        :: d2rhoj20sdr2
  real(dp)                        :: d2rhoi20_2NNdr2
  real(dp)                        :: d2rhoj20_2NNdr2
  real(dp)                        :: gofrho2
  real(dp)                        :: dgofrho2dr
  real(dp)                        :: d2gofrho2dr2
  real(dp)                        :: rhoi2
  real(dp)                        :: rhoj2
  real(dp)                        :: rhoi2l
  real(dp)                        :: rhoj2l
  real(dp)                        :: rhoi2s
  real(dp)                        :: rhoj2s
  real(dp)                        :: rhoi20
  real(dp)                        :: rhoj20
  real(dp)                        :: rhoi20s
  real(dp)                        :: rhoj20s
  real(dp)                        :: rhoi20_2NN
  real(dp)                        :: rhoj20_2NN
  real(dp)                        :: rhoi2total
  real(dp)                        :: rhoj2total
  real(dp)                        :: rtrm1
  real(dp)                        :: rtrm2
  real(dp)                        :: scai
  real(dp)                        :: scaj
  real(dp)                        :: Z2Si
  real(dp)                        :: Z2Sj
#ifdef TRACE
  call trace_in('psirho')
#endif
!
!  Check potential types
!
  if (nptype(npot).ne.45.and.nptype(npot).ne.55) then
    rhoi = 0.0_dp
    rhoj = 0.0_dp
    if (lgrad1) then
      drhoidr = 0.0_dp
      drhojdr = 0.0_dp
      if (lgrad2) then
        d2rhoidr2 = 0.0_dp
        d2rhojdr2 = 0.0_dp
      endif
    endif
#ifdef TRACE
    call trace_out('psirho')
#endif
    return
  endif
!
!  Check order of potential
!
  if (lmatch(nati,ntypi,nspec1(npot),nptyp1(npot),.true.).and. &
      lmatch(natj,ntypj,nspec2(npot),nptyp2(npot),.true.)) then
    i5 = 5
    i7 = 7
    i6 = 6
    i8 = 8
  else
    i5 = 7
    i7 = 5
    i6 = 8
    i8 = 6
  endif
!
!  2NN terms
!
  Z2Si = tpot(i5,npot)  ! Z2 x S => number of second neighbours multiplied by screening factor for i
  Z2Sj = tpot(i7,npot)  ! Z2 x S => number of second neighbours multiplied by screening factor for j
  scai = tpot(i6,npot)  ! Scale factor to convert first neighbour distance to second neighbour one for i
  scaj = tpot(i8,npot)  ! Scale factor to convert first neighbour distance to second neighbour one for j
!
!  Search for matching species for i
!
  mi = 0
  lfound = .false.
  do while (.not.lfound.and.mi.lt.neamspec) 
    mi = mi + 1
    lfound = ((nati.eq.neamnat(mi).and.(ntypi.eq.neamtyp(mi).or.neamtyp(mi).eq.0)).and. &
              (neamnat2(mi).eq.0.or.(neamnat2(mi).eq.nati.and.(ntypi.eq.neamtyp2(mi).or.neamtyp2(mi).eq.0))))
  enddo
  if (mi.eq.0) then
    call outerror('no MEAM functional available for atom in Baskes potential',0_i4)
    call stopnow('psirho')
  endif
!
!  Search for matching species for j
!
  mj = 0
  lfound = .false.
  do while (.not.lfound.and.mj.lt.neamspec) 
    mj = mj + 1
    lfound = ((natj.eq.neamnat(mj).and.(ntypj.eq.neamtyp(mj).or.neamtyp(mj).eq.0)).and. &
              (neamnat2(mj).eq.0.or.(neamnat2(mj).eq.natj.and.(ntypj.eq.neamtyp2(mj).or.neamtyp2(mj).eq.0))))
  enddo
  if (mj.eq.0) then
    call outerror('no MEAM functional available for atom in Baskes potential',0_i4)
    call stopnow('psirho')
  endif
!**********************************************************
!  Compute the reference densities for each lattice type  *
!**********************************************************
  if (ipot(2,npot).eq.1.or.ipot(2,npot).eq.2.or.ipot(2,npot).eq.3.or.ipot(2,npot).eq.4.or. &
      ipot(2,npot).eq.8.or.ipot(2,npot).eq.9) then
!------------------------------------------
!  FCC / BCC / NaCl / HCP / Diamond / ZnS |
!------------------------------------------
!
!  MEAM 0th order density - first neighbours
!
    betaoverr0i = 2.0_dp*denpar(2,1,1,mj)/denpar(3,1,1,mj)
    betaoverr0j = 2.0_dp*denpar(2,1,1,mi)/denpar(3,1,1,mi)
    rhoi20 = (denpar(1,1,1,mj)**2)*meamlatticecoeff(1,ipot(2,npot))*exp(-betaoverr0i*(r-denpar(3,1,1,mj)))
    rhoj20 = (denpar(1,1,1,mi)**2)*meamlatticecoeff(1,ipot(2,npot))*exp(-betaoverr0j*(r-denpar(3,1,1,mi)))
    if (lgrad1) then
      drhoi20dr = - rhoi20*betaoverr0i
      drhoj20dr = - rhoj20*betaoverr0j
      if (lgrad2) then
        d2rhoi20dr2 = rhoi20*betaoverr0i**2
        d2rhoj20dr2 = rhoj20*betaoverr0j**2
      endif
    endif
!
!  MEAM 0th order density - second neighbours
!
    if (Z2Si.gt.0.0_dp) then
      betaoverr0i = 2.0_dp*denpar(2,1,1,mi)/denpar(3,1,1,mi)
      rhoi20_2NN = Z2Si*Z2Si*(denpar(1,1,1,mi)**2)*exp(-betaoverr0i*(scai*r-denpar(3,1,1,mi)))
      if (lgrad1) then
        drhoi20_2NNdr = - rhoi20_2NN*betaoverr0i*scai
        if (lgrad2) then
          d2rhoi20_2NNdr2 = rhoi20_2NN*(betaoverr0i*scai)**2
        endif
      endif
    else
      rhoi20_2NN = 0.0_dp
      if (lgrad1) then
        drhoi20_2NNdr = 0.0_dp
        if (lgrad2) then
          d2rhoi20_2NNdr2 = 0.0_dp
        endif
      endif
    endif
    if (Z2Sj.gt.0.0_dp) then
      betaoverr0j = 2.0_dp*denpar(2,1,1,mj)/denpar(3,1,1,mj)
      rhoj20_2NN = Z2Sj*Z2Sj*(denpar(1,1,1,mj)**2)*exp(-betaoverr0j*(scaj*r-denpar(3,1,1,mj)))
      if (lgrad1) then
        drhoj20_2NNdr = - rhoj20_2NN*betaoverr0j*scaj
        if (lgrad2) then
          d2rhoj20_2NNdr2 = rhoj20_2NN*(betaoverr0j*scaj)**2
        endif
      endif
    else
      rhoj20_2NN = 0.0_dp
      if (lgrad1) then
        drhoj20_2NNdr = 0.0_dp
        if (lgrad2) then
          d2rhoj20_2NNdr2 = 0.0_dp
        endif
      endif
    endif
!
!  MEAM higher order densities
!
    rhoi2 = 0.0_dp
    drhoi2dr = 0.0_dp
    d2rhoi2dr2 = 0.0_dp
    rhoj2 = 0.0_dp
    drhoj2dr = 0.0_dp
    d2rhoj2dr2 = 0.0_dp
!
!  Loop over MEAM higher orders
!
    do l = 2,neammeamorder(1,mj)
      betaoverr0i = 2.0_dp*denpar(2,l,1,mj)/denpar(3,l,1,mj)
      rhoi2l = (denpar(1,l,1,mj)**2)*meamlatticecoeff(l,ipot(2,npot))* &
                eamfnmeamcoeff(l,neamspeci)*exp(-betaoverr0i*(r-denpar(3,l,1,mj)))
      rhoi2 = rhoi2 + rhoi2l
      if (lgrad1) then
        drhoi2dr = drhoi2dr - rhoi2l*betaoverr0i
        if (lgrad2) then
          d2rhoi2dr2 = d2rhoi2dr2 + rhoi2l*betaoverr0i**2
        endif
      endif
    enddo
    do l = 2,neammeamorder(1,mi)
      betaoverr0j = 2.0_dp*denpar(2,l,1,mi)/denpar(3,l,1,mi)
      rhoj2l = (denpar(1,l,1,mi)**2)*meamlatticecoeff(l,ipot(2,npot))* &
                eamfnmeamcoeff(l,neamspecj)*exp(-betaoverr0j*(r-denpar(3,l,1,mi)))
      rhoj2 = rhoj2 + rhoj2l
      if (lgrad1) then
        drhoj2dr = drhoj2dr - rhoj2l*betaoverr0j
        if (lgrad2) then
          d2rhoj2dr2 = d2rhoj2dr2 + rhoj2l*betaoverr0j**2
        endif
      endif
    enddo
  elseif (ipot(2,npot).eq.5) then
!----------
!  Dimer  |
!----------
!
!  MEAM 0th order density - first neighbours
!
    betaoverr0i = 2.0_dp*denpar(2,1,1,mj)/denpar(3,1,1,mj)
    betaoverr0j = 2.0_dp*denpar(2,1,1,mi)/denpar(3,1,1,mi)
    rhoi20 = (denpar(1,1,1,mj)**2)*meamlatticecoeff(1,ipot(2,npot))*exp(-betaoverr0i*(r-denpar(3,1,1,mj)))
    rhoj20 = (denpar(1,1,1,mi)**2)*meamlatticecoeff(1,ipot(2,npot))*exp(-betaoverr0j*(r-denpar(3,1,1,mi)))
    if (lgrad1) then
      drhoi20dr = - rhoi20*betaoverr0i
      drhoj20dr = - rhoj20*betaoverr0j
      if (lgrad2) then
        d2rhoi20dr2 = rhoi20*betaoverr0i**2
        d2rhoj20dr2 = rhoj20*betaoverr0j**2
      endif
    endif
!
!  MEAM 0th order density - for a dimer there are no second neighbours
!
    rhoi20_2NN = 0.0_dp
    if (lgrad1) then
      drhoi20_2NNdr = 0.0_dp
      if (lgrad2) then
        d2rhoi20_2NNdr2 = 0.0_dp
      endif
    endif
    rhoj20_2NN = 0.0_dp
    if (lgrad1) then
      drhoj20_2NNdr = 0.0_dp
      if (lgrad2) then
        d2rhoj20_2NNdr2 = 0.0_dp
      endif
    endif
!
!  MEAM higher order densities
!
    rhoi2 = 0.0_dp
    drhoi2dr = 0.0_dp
    d2rhoi2dr2 = 0.0_dp
    rhoj2 = 0.0_dp
    drhoj2dr = 0.0_dp
    d2rhoj2dr2 = 0.0_dp
!
!  Loop over MEAM higher orders
!
    do l = 2,neammeamorder(1,mj)
      betaoverr0i = 2.0_dp*denpar(2,l,1,mj)/denpar(3,l,1,mj)
      rhoi2l = (denpar(1,l,1,mj)**2)*meamlatticecoeff(l,ipot(2,npot))* &
                eamfnmeamcoeff(l,neamspeci)*exp(-betaoverr0i*(r-denpar(3,l,1,mj)))
      rhoi2 = rhoi2 + rhoi2l
      if (lgrad1) then
        drhoi2dr = drhoi2dr - rhoi2l*betaoverr0i
        if (lgrad2) then
          d2rhoi2dr2 = d2rhoi2dr2 + rhoi2l*betaoverr0i**2
        endif
      endif
    enddo
    do l = 2,neammeamorder(1,mi)
      betaoverr0j = 2.0_dp*denpar(2,l,1,mi)/denpar(3,l,1,mi)
      rhoj2l = (denpar(1,l,1,mi)**2)*meamlatticecoeff(l,ipot(2,npot))* &
                eamfnmeamcoeff(l,neamspecj)*exp(-betaoverr0j*(r-denpar(3,l,1,mi)))
      rhoj2 = rhoj2 + rhoj2l
      if (lgrad1) then
        drhoj2dr = drhoj2dr - rhoj2l*betaoverr0j
        if (lgrad2) then
          d2rhoj2dr2 = d2rhoj2dr2 + rhoj2l*betaoverr0j**2
        endif
      endif
    enddo
  elseif ((ipot(2,npot).eq.6.and.lorder12(npot)).or.(ipot(2,npot).eq.7.and..not.lorder12(npot))) then
!-----------
!  L12AB3  |
!-----------
!
!  MEAM 0th order density - first neighbours : i
!
    betaoverr0j = 2.0_dp*denpar(2,1,1,mj)/denpar(3,1,1,mj)
    rhoi20 = 144.0_dp*(denpar(1,1,1,mj)**2)*exp(-betaoverr0j*(r-denpar(3,1,1,mj)))
    if (lgrad1) then
      drhoi20dr = - rhoi20*betaoverr0j
      if (lgrad2) then
        d2rhoi20dr2 = rhoi20*betaoverr0j**2
      endif
    endif
!
!  MEAM 0th order density - first neighbours : j. For j this shell is a mix of i and j
!
    betaoverr0i = denpar(2,1,1,mi)/denpar(3,1,1,mi)
    betaoverr0j = denpar(2,1,1,mj)/denpar(3,1,1,mj)
    rtrm1 = 4.0_dp*denpar(1,1,1,mi)*exp(-betaoverr0i*(r-denpar(3,1,1,mi)))
    rtrm2 = 8.0_dp*denpar(1,1,1,mj)*exp(-betaoverr0j*(r-denpar(3,1,1,mj)))
    rhoj20 = (rtrm1 + rtrm2)**2
    if (lgrad1) then
      drhoj20dr = - 2.0_dp*(rtrm1 + rtrm2)*(betaoverr0i*rtrm1 + betaoverr0j*rtrm2)
      if (lgrad2) then
        d2rhoj20dr2 = 2.0_dp*(betaoverr0i*rtrm1 + betaoverr0j*rtrm2)**2 + &
                      2.0_dp*(rtrm1 + rtrm2)*(rtrm1*betaoverr0i**2 + rtrm2*betaoverr0j**2)
      endif
    endif
!
!  MEAM 0th order density - second neighbours
!
    if (Z2Si.gt.0.0_dp) then
      betaoverr0i = 2.0_dp*denpar(2,1,1,mi)/denpar(3,1,1,mi)
      rhoi20_2NN = ((Z2Si*denpar(1,1,1,mi))**2)*exp(-betaoverr0i*(scai*r-denpar(3,1,1,mi)))
      if (lgrad1) then
        drhoi20_2NNdr = - rhoi20_2NN*betaoverr0i*scai
        if (lgrad2) then
          d2rhoi20_2NNdr2 = rhoi20_2NN*(betaoverr0i*scai)**2
        endif
      endif
    else
      rhoi20_2NN = 0.0_dp
      if (lgrad1) then
        drhoi20_2NNdr = 0.0_dp
        if (lgrad2) then
          d2rhoi20_2NNdr2 = 0.0_dp
        endif
      endif
    endif
    if (Z2Sj.gt.0.0_dp) then
      betaoverr0j = 2.0_dp*denpar(2,1,1,mj)/denpar(3,1,1,mj)
      rhoj20_2NN = ((Z2Sj*denpar(1,1,1,mj))**2)*exp(-betaoverr0j*(scaj*r-denpar(3,1,1,mj)))
      if (lgrad1) then
        drhoj20_2NNdr = - rhoj20_2NN*betaoverr0j*scaj
        if (lgrad2) then
          d2rhoj20_2NNdr2 = rhoj20_2NN*(betaoverr0j*scaj)**2
        endif
      endif
    else
      rhoj20_2NN = 0.0_dp
      if (lgrad1) then
        drhoj20_2NNdr = 0.0_dp
        if (lgrad2) then
          d2rhoj20_2NNdr2 = 0.0_dp
        endif
      endif
    endif
!
!  MEAM higher order densities
!
    rhoi2 = 0.0_dp
    drhoi2dr = 0.0_dp
    d2rhoi2dr2 = 0.0_dp
    rhoj2 = 0.0_dp
    drhoj2dr = 0.0_dp
    d2rhoj2dr2 = 0.0_dp
!
!  Loop over MEAM higher orders: For this structure only l = 3 for j is non-zero
!
    l = 3
    betaoverr0i = denpar(2,l,1,mj)/denpar(3,l,1,mj)
    betaoverr0j = denpar(2,l,1,mi)/denpar(3,l,1,mi)
    rtrm1 = denpar(1,l,1,mj)*exp(-betaoverr0i*(r-denpar(3,l,1,mj)))
    rtrm2 = denpar(1,l,1,mi)*exp(-betaoverr0j*(r-denpar(3,l,1,mi)))
    ctrm  = 8.0_dp*eamfnmeamcoeff(l,neamspecj)/3.0_dp
    rhoj2 = ctrm*(rtrm1 - rtrm2)**2
    if (lgrad1) then
      drhoj2dr = - 2.0_dp*ctrm*(rtrm1 - rtrm2)*(betaoverr0i*rtrm1 - betaoverr0j*rtrm2)
      if (lgrad2) then
        d2rhoj2dr2 = 2.0_dp*ctrm*((betaoverr0i*rtrm1 - betaoverr0j*rtrm2)**2 + &
                                  (rtrm1 - rtrm2)*(rtrm1*betaoverr0i**2 - rtrm2*betaoverr0j**2))
      endif
    endif
  elseif ((ipot(2,npot).eq.7.and.lorder12(npot)).or.(ipot(2,npot).eq.6.and..not.lorder12(npot))) then
!-----------
!  L12A3B  |
!-----------
!
!  MEAM 0th order density - first neighbours : For i this shell is a mix of i and j
!
    betaoverr0i = denpar(2,1,1,mi)/denpar(3,1,1,mi)
    betaoverr0j = denpar(2,1,1,mj)/denpar(3,1,1,mj)
    rtrm1 = 4.0_dp*denpar(1,1,1,mj)*exp(-betaoverr0j*(r-denpar(3,1,1,mj)))
    rtrm2 = 8.0_dp*denpar(1,1,1,mi)*exp(-betaoverr0i*(r-denpar(3,1,1,mi)))
    rhoi20 = (rtrm1 + rtrm2)**2
    if (lgrad1) then
      drhoi20dr = - 2.0_dp*(rtrm1 + rtrm2)*(betaoverr0j*rtrm1 + betaoverr0i*rtrm2)
      if (lgrad2) then
        d2rhoi20dr2 = 2.0_dp*(betaoverr0j*rtrm1 + betaoverr0i*rtrm2)**2 + &
                      2.0_dp*(rtrm1 + rtrm2)*(rtrm1*betaoverr0j**2 + rtrm2*betaoverr0i**2)
      endif
    endif
!
    betaoverr0i = 2.0_dp*denpar(2,1,1,mi)/denpar(3,1,1,mi)
    rhoj20 = 144.0_dp*(denpar(1,1,1,mi)**2)*exp(-betaoverr0i*(r-denpar(3,1,1,mi)))
    if (lgrad1) then
      drhoj20dr = - rhoj20*betaoverr0i
      if (lgrad2) then
        d2rhoj20dr2 = rhoj20*betaoverr0i**2
      endif
    endif
!
!  MEAM 0th order density - second neighbours
!
    if (Z2Si.gt.0.0_dp) then
      betaoverr0i = denpar(2,1,1,mj)/denpar(3,1,1,mj)
      rhoi20_2NN = Z2Si*denpar(1,1,1,mj)*exp(-betaoverr0i*(scai*r-denpar(3,1,1,mj)))
      if (lgrad1) then
        drhoi20_2NNdr = - rhoi20_2NN*betaoverr0i*scai
        if (lgrad2) then
          d2rhoi20_2NNdr2 = rhoi20_2NN*(betaoverr0i*scai)**2
        endif
      endif
    else
      rhoi20_2NN = 0.0_dp
      if (lgrad1) then
        drhoi20_2NNdr = 0.0_dp
        if (lgrad2) then
          d2rhoi20_2NNdr2 = 0.0_dp
        endif
      endif
    endif
    if (Z2Sj.gt.0.0_dp) then
      betaoverr0j = 2.0_dp*denpar(2,1,1,mj)/denpar(3,1,1,mj)
      rhoj20_2NN = Z2Sj*Z2Sj*(denpar(1,1,1,mj)**2)*exp(-betaoverr0j*(scaj*r-denpar(3,1,1,mj)))
      if (lgrad1) then
        drhoj20_2NNdr = - rhoj20_2NN*betaoverr0j*scaj
        if (lgrad2) then
          d2rhoj20_2NNdr2 = rhoj20_2NN*(betaoverr0j*scaj)**2
        endif
      endif
    else
      rhoj20_2NN = 0.0_dp
      if (lgrad1) then
        drhoj20_2NNdr = 0.0_dp
        if (lgrad2) then
          d2rhoj20_2NNdr2 = 0.0_dp
        endif
      endif
    endif
!
!  MEAM higher order densities
!
    rhoi2 = 0.0_dp
    drhoi2dr = 0.0_dp
    d2rhoi2dr2 = 0.0_dp
    rhoj2 = 0.0_dp
    drhoj2dr = 0.0_dp
    d2rhoj2dr2 = 0.0_dp
!
!  Loop over MEAM higher orders: For this structure only l = 3 for i is non-zero
!
    l = 3
    betaoverr0i = denpar(2,l,1,mj)/denpar(3,l,1,mj)
    betaoverr0j = denpar(2,l,1,mi)/denpar(3,l,1,mi)
    rtrm1 = denpar(1,l,1,mj)*exp(-betaoverr0i*(r-denpar(3,l,1,mj)))
    rtrm2 = denpar(1,l,1,mi)*exp(-betaoverr0j*(r-denpar(3,l,1,mi)))
    ctrm  = 8.0_dp*eamfnmeamcoeff(l,neamspeci)/3.0_dp
    rhoi2 = ctrm*(rtrm1 - rtrm2)**2
    if (lgrad1) then
      drhoi2dr = - 2.0_dp*ctrm*(rtrm1 - rtrm2)*(betaoverr0i*rtrm1 - betaoverr0j*rtrm2)
      if (lgrad2) then
        d2rhoi2dr2 = 2.0_dp*ctrm*((betaoverr0i*rtrm1 - betaoverr0j*rtrm2)**2 + &
                                  (rtrm1 - rtrm2)*(rtrm1*betaoverr0i**2 - rtrm2*betaoverr0j**2))
      endif
    endif
  endif
!
!  Compute the total density from the components
!
  if (nmeamcombotype.eq.1) then
!-----------------------------
!  MEAM 1NN combination rule |
!-----------------------------
!
!  Add order = 0 density squared to total
!
    if (Z2Si.gt.0.0_dp) then
      rhoi2s = rhoi2 + rhoi20
      rhoi2total = (sqrt(rhoi2s) + sqrt(rhoi20_2NN))**2
      if (lgrad1) then
        drhoi2sdr = drhoi2dr + drhoi20dr
        drhoidr = sqrt(rhoi2total)*(drhoi2sdr/sqrt(rhoi2s) + drhoi20_2NNdr/sqrt(rhoi20_2NN))
        if (lgrad2) then
          d2rhoi2sdr2 = d2rhoi2dr2 + d2rhoi20dr2
          d2rhoidr2 = 0.5_dp*(drhoi2sdr/sqrt(rhoi2s) + drhoi20_2NNdr/sqrt(rhoi20_2NN))**2 + &
                     sqrt(rhoi2total)*(d2rhoi2sdr2/sqrt(rhoi2s) + d2rhoi20_2NNdr2/sqrt(rhoi20_2NN) &
                     - 0.5_dp*((drhoi2sdr**2)/(rhoi2s**1.5_dp) + (drhoi20_2NNdr**2)/(rhoi20_2NN**1.5_dp)))
        endif
      endif
    else
      rhoi2total = rhoi2 + rhoi20
      if (lgrad1) then
        drhoidr = drhoi2dr + drhoi20dr
        if (lgrad2) then
          d2rhoidr2 = d2rhoi2dr2 + d2rhoi20dr2
        endif
      endif
    endif
!
    if (Z2Sj.gt.0.0_dp) then
      rhoj2s = rhoj2 + rhoj20
      rhoj2total = (sqrt(rhoj2s) + sqrt(rhoj20_2NN))**2
      if (lgrad1) then
        drhoj2sdr = drhoj2dr + drhoj20dr
        drhojdr = sqrt(rhoj2total)*(drhoj2sdr/sqrt(rhoj2s) + drhoj20_2NNdr/sqrt(rhoj20_2NN))
        if (lgrad2) then
          d2rhoj2sdr2 = d2rhoj2dr2 + d2rhoj20dr2
          d2rhojdr2 = 0.5_dp*(drhoj2sdr/sqrt(rhoj2s) + drhoj20_2NNdr/sqrt(rhoj20_2NN))**2 + &
                     sqrt(rhoj2total)*(d2rhoj2sdr2/sqrt(rhoj2s) + d2rhoj20_2NNdr2/sqrt(rhoj20_2NN) &
                     - 0.5_dp*((drhoj2sdr**2)/(rhoj2s**1.5_dp) + (drhoj20_2NNdr**2)/(rhoj20_2NN**1.5_dp)))
        endif
      endif
    else
      rhoj2total = rhoj2 + rhoj20
      if (lgrad1) then
        drhojdr = drhoj2dr + drhoj20dr
        if (lgrad2) then
          d2rhojdr2 = d2rhoj2dr2 + d2rhoj20dr2
        endif
      endif
    endif
  else
!-----------------------------
!  MEAM 2NN combination rule |
!-----------------------------
!
!  If rhoi20 is less than a threshold then need to trap
!
    if (rhoi20.lt.1.0d-12) then
!
!  Here we assume that if the spherical contribution to rho is close to zero, 
!  that the higher order harmonics will be even closer to zero and therefore that
!  G(gamma) is tending to 1. Error in this limit is not likely to matter given 
!  that prefactor is very small.
!
      rhoi2total = rhoi20
      if (lgrad1) then
        drhoidr = drhoi20dr
        if (lgrad2) then
          d2rhoidr2 = d2rhoi20dr2
        endif
      endif
      rhoi2total = rhoi20
    else
!
!  General case
!
      if (Z2Si.gt.0.0_dp) then
        rhoi20s = rhoi20
        rhoi20 = (sqrt(rhoi20s) + sqrt(rhoi20_2NN))**2
        if (lgrad1) then
          drhoi20sdr = drhoi20dr
          drhoi20dr = sqrt(rhoi20)*(drhoi20sdr/sqrt(rhoi20s) + drhoi20_2NNdr/sqrt(rhoi20_2NN))
          if (lgrad2) then
            d2rhoi20sdr2 = d2rhoi20dr2
            d2rhoi20dr2 = 0.5_dp*(drhoi20sdr/sqrt(rhoi20s) + drhoi20_2NNdr/sqrt(rhoi20_2NN))**2 + &
                         sqrt(rhoi20)*(d2rhoi20sdr2/sqrt(rhoi20s) + d2rhoi20_2NNdr2/sqrt(rhoi20_2NN) &
                         - 0.5_dp*((drhoi20sdr**2)/(rhoi20s**1.5_dp) + (drhoi20_2NNdr**2)/(rhoi20_2NN**1.5_dp)))
          endif
        endif
      endif
      gamma    = rhoi2/rhoi20
      exptrm   = exp(-gamma)
      gofrho2  = 2.0_dp/(1.0_dp + exptrm)
      if (gofrho2.lt.1.0d-15) then
        rhoi2total = 0.0_dp
        if (lgrad1) then
          drhoidr = 0.0_dp
          if (lgrad2) then
            d2rhoidr2 = 0.0_dp
          endif
        endif
      else
        rhoi2total = rhoi20*gofrho2**2
        if (lgrad1) then
          dgammadr = (drhoi2dr/rhoi20) - (rhoi2/rhoi20**2)*drhoi20dr
          dgofrho2dr = (2.0_dp*exptrm/(1.0_dp + exptrm)**2)*dgammadr
          drhoidr = drhoi20dr*gofrho2**2 + 2.0_dp*gofrho2*rhoi20*dgofrho2dr
          if (lgrad2) then
            d2gammadr2 = (d2rhoi2dr2/rhoi20) - 2.0_dp*(drhoi2dr/rhoi20**2)*drhoi20dr + &
                         (2.0_dp*rhoi2/rhoi20**3)*drhoi20dr*drhoi20dr - (rhoi2/rhoi20**2)*d2rhoi20dr2
            d2gofrho2dr2 = (2.0_dp*exptrm/(1.0_dp + exptrm)**2)*d2gammadr2 - &
                           (2.0_dp*exptrm/(1.0_dp + exptrm)**2)*dgammadr*dgammadr + &
                           (4.0_dp*exptrm*exptrm/(1.0_dp + exptrm)**3)*dgammadr*dgammadr
            d2rhoidr2 = d2rhoi20dr2*gofrho2**2 + 4.0_dp*drhoi20dr*gofrho2*dgofrho2dr + &
                       2.0_dp*rhoi20*(dgofrho2dr*dgofrho2dr + gofrho2*d2gofrho2dr2)
          endif
        endif
      endif
    endif
!
!  If rhoj20 is less than a threshold then need to trap
!
    if (rhoj20.lt.1.0d-12) then
!
!  Here we assume that if the spherical contribution to rho is close to zero,
!  that the higher order harmonics will be even closer to zero and therefore that
!  G(gamma) is tending to 1. Error in this limit is not likely to matter given
!  that prefactor is very small.
!
      rhoj2total = rhoj20
      if (lgrad1) then
        drhojdr = drhoj20dr
        if (lgrad2) then
          d2rhojdr2 = d2rhoj20dr2
        endif
      endif
      rhoj2total = rhoj20
    else
!
!  General case
!
      if (Z2Sj.gt.0.0_dp) then
        rhoj20s = rhoj20
        rhoj20 = (sqrt(rhoj20s) + sqrt(rhoj20_2NN))**2
        if (lgrad1) then
          drhoj20sdr = drhoj20dr
          drhoj20dr = sqrt(rhoj20)*(drhoj20sdr/sqrt(rhoj20s) + drhoj20_2NNdr/sqrt(rhoj20_2NN))
          if (lgrad2) then
            d2rhoj20sdr2 = d2rhoj20dr2
            d2rhoj20dr2 = 0.5_dp*(drhoj20sdr/sqrt(rhoj20s) + drhoj20_2NNdr/sqrt(rhoj20_2NN))**2 + &
                         sqrt(rhoj20)*(d2rhoj20sdr2/sqrt(rhoj20s) + d2rhoj20_2NNdr2/sqrt(rhoj20_2NN) &
                         - 0.5_dp*((drhoj20sdr**2)/(rhoj20s**1.5_dp) + (drhoj20_2NNdr**2)/(rhoj20_2NN**1.5_dp)))
          endif
        endif
      endif
      gamma    = rhoj2/rhoj20
      exptrm   = exp(-gamma)
      gofrho2  = 2.0_dp/(1.0_dp + exptrm)
      if (gofrho2.lt.1.0d-15) then
        rhoj2total = 0.0_dp
        if (lgrad1) then
          drhojdr = 0.0_dp
          if (lgrad2) then
            d2rhojdr2 = 0.0_dp
          endif
        endif
      else
        rhoj2total = rhoj20*gofrho2**2
        if (lgrad1) then
          dgammadr = (drhoj2dr/rhoj20) - (rhoj2/rhoj20**2)*drhoj20dr
          dgofrho2dr = (2.0_dp*exptrm/(1.0_dp + exptrm)**2)*dgammadr
          drhojdr = drhoj20dr*gofrho2**2 + 2.0_dp*gofrho2*rhoj20*dgofrho2dr
          if (lgrad2) then
            d2gammadr2 = (d2rhoj2dr2/rhoj20) - 2.0_dp*(drhoj2dr/rhoj20**2)*drhoj20dr + &
                         (2.0_dp*rhoj2/rhoj20**3)*drhoj20dr*drhoj20dr - (rhoj2/rhoj20**2)*d2rhoj20dr2
            d2gofrho2dr2 = (2.0_dp*exptrm/(1.0_dp + exptrm)**2)*d2gammadr2 - &
                           (2.0_dp*exptrm/(1.0_dp + exptrm)**2)*dgammadr*dgammadr + &
                           (4.0_dp*exptrm*exptrm/(1.0_dp + exptrm)**3)*dgammadr*dgammadr
            d2rhojdr2 = d2rhoj20dr2*gofrho2**2 + 4.0_dp*drhoj20dr*gofrho2*dgofrho2dr + &
                       2.0_dp*rhoj20*(dgofrho2dr*dgofrho2dr + gofrho2*d2gofrho2dr2)
          endif
        endif
      endif
    endif
  endif
!
!  Square root to obtain total density
!
  rhoi = sqrt(rhoi2total)
  rhoj = sqrt(rhoj2total)
#ifdef TRACE
  call trace_out('psirho')
#endif
!
  return
  end
