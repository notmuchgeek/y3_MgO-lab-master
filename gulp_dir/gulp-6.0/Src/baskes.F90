  subroutine baskes(nati,ntypi,neamspeci,natj,ntypj,neamspecj,npot,r,rk,psfct, &
                    ebas,d1bas,d2bas,lgrad1,lgrad2)
!
!  Subroutine for calculating psi(R) in MEAM twobody potential
!
!  On return deriv and deriv2 contain the complete first and
!  second derivative terms, while derive and derive2 contain 
!  the derivatives for the charge only terms.
!
!   8/14 Created from psibaskes/density/many routines
!   9/14 Alloys added
!   1/18 Trace added
!  10/18 Second derivative terms added
!  11/18 Tapering added to baskes potential
!   1/21 Checking of convergence of iterative screening to prevent
!        possible underflow and save work
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
!  Copyright Curtin University 2021
!
!  Julian Gale, CIC, Curtin University, January 2021
!
  use eam
#ifdef TRACE
  use trace,      only : trace_in, trace_out
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
  real(dp),       intent(in)      :: r      ! Distance
  real(dp),       intent(in)      :: rk     ! Inverse distance
  real(dp),       intent(in)      :: psfct  ! Scale factor
  real(dp),       intent(out)     :: ebas   ! Energy of potential
  real(dp),       intent(out)     :: d1bas  ! First derivative of potential divided by r
  real(dp),       intent(out)     :: d2bas  ! Second derivative of potential 
!
!  Local variables
!
  integer(i4)                     :: i5
  integer(i4)                     :: i6
  integer(i4)                     :: i7
  integer(i4)                     :: i8
  integer(i4)                     :: n
  integer(i4)                     :: nselfi
  integer(i4)                     :: nselfj
  logical                         :: ldoi
  logical                         :: ldoj
  logical                         :: lmatch
  logical                         :: lswap
  real(dp)                        :: apt
  real(dp)                        :: cpti
  real(dp)                        :: cptj
  real(dp)                        :: d1bas0
  real(dp)                        :: d1basi
  real(dp)                        :: d1basij
  real(dp)                        :: d1basj
  real(dp)                        :: d1bast0
  real(dp)                        :: d1basti
  real(dp)                        :: d1bastj
  real(dp)                        :: d2bas0
  real(dp)                        :: d2basi
  real(dp)                        :: d2basij
  real(dp)                        :: d2basj
  real(dp)                        :: d2bast0
  real(dp)                        :: d2basti
  real(dp)                        :: d2bastj
  real(dp)                        :: dtpfn
  real(dp)                        :: d2tpfn
  real(dp)                        :: d3tpfn
  real(dp)                        :: ebas0
  real(dp)                        :: ebasi
  real(dp)                        :: ebasij
  real(dp)                        :: ebasj
  real(dp)                        :: ebasold
  real(dp)                        :: ebast0
  real(dp)                        :: ebasti
  real(dp)                        :: ebastj
  real(dp)                        :: oneZ
  real(dp)                        :: s
  real(dp)                        :: s0
  real(dp)                        :: sn
  real(dp)                        :: sr
  real(dp)                        :: srk
  real(dp)                        :: tpfn
  real(dp)                        :: wi
  real(dp)                        :: wj
  real(dp)                        :: Z2SoverZ1
#ifdef TRACE
  call trace_in('baskes')
#endif
!
  ebas = 0.0_dp
  if (lgrad1) then
    d1bas = 0.0_dp
    if (lgrad2) then
      d2bas = 0.0_dp
    endif
  endif
!
!  Find weight of functionals taking care with regard to the order
!
  if (lmatch(nati,ntypi,nspec1(npot),nptyp1(npot),.true.).and. &
      lmatch(natj,ntypj,nspec2(npot),nptyp2(npot),.true.)) then
    wi = twopot(6,npot)
    wj = 1.0_dp - wi
    nselfi = ipot(3,npot)
    nselfj = ipot(4,npot)
    i5 = 5
    i7 = 7
    i6 = 6
    i8 = 8
    lswap = .false.
  else
    wj = twopot(6,npot)
    wi = 1.0_dp - wj
    nselfi = ipot(4,npot)
    nselfj = ipot(3,npot)
    i5 = 7
    i7 = 5
    i6 = 8
    i8 = 6
    lswap = .true.
  endif
!
  oneZ = - psfct/twopot(5,npot) ! - (1/Z)
!
  if (nselfi.ne.nselfj.and.ipot(1,npot).ne.5) then
!*************************
!  Mixed term for alloy  *
!*************************
!-----------------------------
!  Compute potential for i-j |
!-----------------------------
    if (lswap) then
      call psibaskes(natj,ntypj,neamspecj,nati,ntypi,neamspeci,npot,r,rk,psfct, &
                     ebas,d1bas,d2bas,lgrad1,lgrad2)
    else
      call psibaskes(nati,ntypi,neamspeci,natj,ntypj,neamspecj,npot,r,rk,psfct, &
                     ebas,d1bas,d2bas,lgrad1,lgrad2)
    endif
!------------------------------------
!  Compute potential for i in alloy |
!------------------------------------
    ebasti = 0.0_dp
    if (lgrad1) then
      d1basti = 0.0_dp
      if (lgrad2) then
        d2basti = 0.0_dp
      endif
    endif
    s0 = tpot(i6,npot)
    s  = tpot(i6,nselfi)
    sr = s0*r
    srk = 1.0_dp/sr
    call psibaskes(nati,ntypi,neamspeci,nati,ntypi,neamspeci,nselfi,sr,srk,psfct, &
                   ebasi,d1basi,d2basi,lgrad1,lgrad2)
    ebasti = ebasti + ebasi
    if (lgrad1) then
      d1basti = d1basti + s0*d1basi
      if (lgrad2) then
        d2basti = d2basti + s0*d2basi
      endif
    endif
!
!  2NN formalism, if requested
!
    if (ipot(1,npot).gt.0) then
      Z2SoverZ1 = tpot(5,nselfi)/twopot(5,nselfi)
      do n = 1,ipot(1,npot)
        apt = (-Z2SoverZ1)**n
!
!  Scale distances
!
        sn = s0*s**n
        sr = sn*r
        srk = 1.0_dp/sr
!
        call psibaskes(nati,ntypi,neamspeci,nati,ntypi,neamspeci,nselfi,sr,srk,psfct, &
                       ebasi,d1basi,d2basi,lgrad1,lgrad2)
!
        ebasti = ebasti + apt*ebasi
        if (lgrad1) then
          d1basti = d1basti + apt*d1basi*sn
          if (lgrad2) then
            d2basti = d2basti + apt*d2basi*sn*sn
          endif
        endif
!
      enddo
    endif
!------------------------------------
!  Compute potential for j in alloy |
!------------------------------------
    ebastj = 0.0_dp
    if (lgrad1) then
      d1bastj = 0.0_dp
      if (lgrad2) then
        d2bastj = 0.0_dp
      endif
    endif
    s0 = tpot(i8,npot)
    s  = tpot(i8,nselfj)
    sr = s0*r
    srk = 1.0_dp/sr
    call psibaskes(natj,ntypj,neamspecj,natj,ntypj,neamspecj,nselfj,sr,srk,psfct, &
                   ebasj,d1basj,d2basj,lgrad1,lgrad2)
    ebastj = ebastj + ebasj
    if (lgrad1) then
      d1bastj = d1bastj + s0*d1basj
      if (lgrad2) then
        d2bastj = d2bastj + s0*d2basj
      endif
    endif
!
!  2NN formalism, if requested
!
    if (ipot(1,npot).gt.0) then
      Z2SoverZ1 = tpot(5,nselfj)/twopot(5,nselfj)
      do n = 1,ipot(1,npot)
        apt = (-Z2SoverZ1)**n
!
!  Scale distances
!
        sn = s0*s**n
        sr = sn*r
        srk = 1.0_dp/sr
!
        call psibaskes(natj,ntypj,neamspecj,natj,ntypj,neamspecj,nselfj,sr,srk,psfct, &
                       ebasj,d1basj,d2basj,lgrad1,lgrad2)
!
        ebastj = ebastj + apt*ebasj
        if (lgrad1) then
          d1bastj = d1bastj + apt*d1basj*sn
          if (lgrad2) then
            d2bastj = d2bastj + apt*d2basj*sn*sn
          endif
        endif
!
      enddo
    endif
!
!  For certain lattice types we need a genuine self term too
!
    ebast0 = 0.0_dp
    if (lgrad1) then
      d1bast0 = 0.0_dp
      if (lgrad2) then
        d2bast0 = 0.0_dp
      endif
    endif
!
!  Work out logic of lattice type and self term choice
!
    ldoi = .false.
    ldoj = .false.
    if (ipot(2,npot).eq.6) then
!
!  L12AB3 => extra term is for B
!
      if (lswap.and.lorder12(npot)) then
        ldoi = .true.
      elseif (.not.lswap.and..not.lorder12(npot)) then
        ldoi = .true.
      else
        ldoj = .true.
      endif
    elseif (ipot(2,npot).eq.7) then
!
!  L12A3B => extra term is for A
!
      if (lswap.and.lorder12(npot)) then
        ldoj = .true.
      elseif (.not.lswap.and..not.lorder12(npot)) then
        ldoj = .true.
      else
        ldoi = .true.
      endif
    endif
!
    if (ldoj) then
      s  = tpot(i6,nselfj)
      call psibaskes(natj,ntypj,neamspecj,natj,ntypj,neamspecj,nselfj,r,rk,psfct, &
                     ebas0,d1bas0,d2bas0,lgrad1,lgrad2)
      ebast0 = ebast0 + ebas0
      if (lgrad1) then
        d1bast0 = d1bast0 + d1bas0
        if (lgrad2) then
          d2bast0 = d2bast0 + d2bas0
        endif
      endif
!
!  2NN formalism, if requested
!
      if (ipot(1,npot).gt.0) then
        Z2SoverZ1 = tpot(5,nselfj)/twopot(5,nselfj)
        do n = 1,ipot(1,npot)
          apt = (-Z2SoverZ1)**n
!
!  Scale distances
!
          sn = s**n
          sr = sn*r
          srk = 1.0_dp/sr
!
          call psibaskes(natj,ntypj,neamspecj,natj,ntypj,neamspecj,nselfj,sr,srk,psfct, &
                         ebas0,d1bas0,d2bas0,lgrad1,lgrad2)
!
          ebast0 = ebast0 + apt*ebas0
          if (lgrad1) then
            d1bast0 = d1bast0 + apt*d1bas0*sn
            if (lgrad2) then
              d2bast0 = d2bast0 + apt*d2bas0*sn*sn
            endif
          endif
        enddo
      endif
    elseif (ldoi) then
      s  = tpot(i6,nselfi)
      call psibaskes(nati,ntypi,neamspeci,nati,ntypi,neamspeci,nselfi,r,rk,psfct, &
                     ebas0,d1bas0,d2bas0,lgrad1,lgrad2)
      ebast0 = ebast0 + ebas0
      if (lgrad1) then
        d1bast0 = d1bast0 + d1bas0
        if (lgrad2) then
          d2bast0 = d2bast0 + d2bas0
        endif
      endif
!
!  2NN formalism, if requested
!
      if (ipot(1,npot).gt.0) then
        Z2SoverZ1 = tpot(5,nselfi)/twopot(5,nselfi)
        do n = 1,ipot(1,npot)
          apt = (-Z2SoverZ1)**n
!
!  Scale distances
!
          sn = s**n
          sr = sn*r
          srk = 1.0_dp/sr
!
          call psibaskes(nati,ntypi,neamspeci,nati,ntypi,neamspeci,nselfi,sr,srk,psfct, &
                         ebas0,d1bas0,d2bas0,lgrad1,lgrad2)
!
          ebast0 = ebast0 + apt*ebas0
          if (lgrad1) then
            d1bast0 = d1bast0 + apt*d1bas0*sn
            if (lgrad2) then
              d2bast0 = d2bast0 + apt*d2bas0*sn*sn
            endif
          endif
        enddo
      endif
    endif
!------------------------------------------
!  Combine potential components together  |
!------------------------------------------
    if (ipot(2,npot).eq.1.or.ipot(2,npot).eq.2.or.ipot(2,npot).eq.3) then
!
!  FCC / BCC / NaCl
!
      cpti = oneZ*wi*tpot(i5,npot)
      cptj = oneZ*wj*tpot(i7,npot)
      ebas = ebas + cpti*ebasti + cptj*ebastj
      if (lgrad1) then
        d1bas = d1bas + cpti*d1basti + cptj*d1bastj
        if (lgrad2) then
          d2bas = d2bas + cpti*d2basti + cptj*d2bastj
        endif
      endif
    elseif (ipot(2,npot).eq.6.or.ipot(2,npot).eq.7) then
!
!  L12AB3 / L12A3B
!
      cpti = oneZ*wi*tpot(i5,npot)
      cptj = oneZ*wj*tpot(i7,npot)
      ebas = 2.0_dp*(ebas + cpti*ebasti + cptj*ebastj) - ebast0
      if (lgrad1) then
        d1bas = 2.0_dp*(d1bas + cpti*d1basti + cptj*d1bastj) - d1bast0
        if (lgrad2) then
          d2bas = 2.0_dp*(d2bas + cpti*d2basti + cptj*d2bastj) - d2bast0
        endif
      endif
    elseif (ipot(2,npot).eq.8) then
!
!  ZnS
!
      cpti = oneZ*wi*tpot(i5,npot)
      cptj = oneZ*wj*tpot(i7,npot)
      ebas = ebas + cpti*ebasti + cptj*ebastj
      if (lgrad1) then
        d1bas = d1bas + cpti*d1basti + cptj*d1bastj
        if (lgrad2) then
          d2bas = d2bas + cpti*d2basti + cptj*d2bastj
        endif
      endif
    endif
  else
!**************
!  Self term  *
!**************
!-----------------------------
!  Compute potential for i-j |
!-----------------------------
    call psibaskes(nati,ntypi,neamspeci,natj,ntypj,neamspecj,npot,r,rk,psfct, &
                   ebas,d1bas,d2bas,lgrad1,lgrad2)
!
!  2NN formalism, if requested
!
    if (ipot(1,npot).gt.0) then
      Z2SoverZ1 = tpot(5,npot)/twopot(5,npot)
      s = tpot(6,npot)
      ebasold = ebas
      do n = 1,ipot(1,npot)
        apt = (-Z2SoverZ1)**n
!
!  Scale distances
!
        sn = s**n
        sr = sn*r
        srk = 1.0_dp/sr
!
        call psibaskes(nati,ntypi,neamspeci,natj,ntypj,neamspecj,npot,sr,srk,psfct, &
                       ebasij,d1basij,d2basij,lgrad1,lgrad2)
!
        ebas = ebas + apt*ebasij
        if (lgrad1) then
          d1bas = d1bas + apt*d1basij*sn
          if (lgrad2) then
            d2bas = d2bas + apt*d2basij*sn*sn
          endif
        endif
!
!  Convergence check
!
        if (abs(ebas-ebasold).lt.1.0d-15) exit
        ebasold = ebas
      enddo
    endif
  endif
!********************
!  Taper potential  *
!********************
  if (r.gt.tapermin.and.r.le.tapermax) then
!
!  Polynomial, cosine, MDF or exponential taper
!
    if (tapertype.eq.1) then
      call p5taper(r,tapermin,tapermax,tpfn,dtpfn,d2tpfn,d3tpfn,lgrad1,lgrad2,.false.)
    elseif (tapertype.eq.4) then
      call etaper(r,tapermin,tapermax,tpfn,dtpfn,d2tpfn,d3tpfn,lgrad1,lgrad2,.false.)
    elseif (tapertype.eq.5) then
      call mdftaper(r,tapermin,tapermax,tpfn,dtpfn,d2tpfn,d3tpfn,lgrad1,lgrad2,.false.)
    elseif (tapertype.eq.6) then
      call p7taper(r,tapermax,tpfn,dtpfn,d2tpfn,d3tpfn,lgrad1,lgrad2,.false.)
    else
      call ctaper(r,tapermin,tapermax,tpfn,dtpfn,d2tpfn,d3tpfn,lgrad1,lgrad2,.false.)
    endif
!
    ebas0 = ebas
    d1bas0 = d1bas
    d2bas0 = d2bas
!
    ebas = tpfn*ebas0
    if (lgrad1) then
      d1bas = tpfn*d1bas0 + dtpfn*ebas0
      if (lgrad2) then
        d2bas = tpfn*d2bas0 + 2.0_dp*dtpfn*d1bas0 + d2tpfn*ebas0
      endif
    endif
  endif
#ifdef TRACE
  call trace_out('baskes')
#endif
!
  return
  end
