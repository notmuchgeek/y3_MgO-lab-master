  subroutine projd3dVd(mcvloc,d3rs,d3is,d3s,i,ix,iy,iz,j,jx,jy,jz,lcorei,lcorej, &
                       rmassi,rmassj,derv2,dervi,eigr,eigi,maxeigc,grueneisen, &
                       d33s,d33rs,d33is,lmany,nmanyk,d34s,d34rs,d34is,nforkl)
!
!  Projects the 3rd derivatives on to the phonon modes to compute volume derivatives
!  Used in the calculation of Grueneisen parameters. Distributed memory version.
!
!   2/18 Created from projd3dV
!
!  NB: Sign of imaginary eigenvectors changed since the sign is inverted
!      relative to those in the free energy derivative routine.
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
!  Julian Gale, CIC, Curtin University, Feb 2018
!
  use current
  use datatypes
  use derivatives,   only : maxd2
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4)         :: i
  integer(i4)         :: ix
  integer(i4)         :: iy
  integer(i4)         :: iz
  integer(i4)         :: j
  integer(i4)         :: jx
  integer(i4)         :: jy
  integer(i4)         :: jz
  integer(i4)         :: maxeigc
  integer(i4)         :: mcvloc
  integer(i4)         :: nmanyk
  integer(i4)         :: nforkl
  logical             :: lcorei
  logical             :: lcorej
  logical             :: lmany
  real(dp)            :: d3d(3,3)
  real(dp)            :: d3is(3,3,6)
  real(dp)            :: d3rs(3,3,6)
  real(dp)            :: d3s(3,3,6)
  real(dp)            :: d33is(108,*)
  real(dp)            :: d33rs(108,*)
  real(dp)            :: d33s(108,*)
  real(dp)            :: d34is(54,*)
  real(dp)            :: d34rs(54,*)
  real(dp)            :: d34s(54,*)
  real(dp)            :: derv2(maxd2,*)
  real(dp)            :: dervi(maxd2,*)
  real(dp)            :: eigr(maxeigc,*)
  real(dp)            :: eigi(maxeigc,*)
  real(dp)            :: grueneisen(*)
  real(dp)            :: rmassi
  real(dp)            :: rmassj
!
!  Local variables
!
  integer(i4)         :: ii
  integer(i4)         :: k
  integer(i4)         :: kk
  real(dp)            :: d1ij11
  real(dp)            :: d1ij21
  real(dp)            :: d1ij31
  real(dp)            :: d1ij12
  real(dp)            :: d1ij22
  real(dp)            :: d1ij32
  real(dp)            :: d1ij13
  real(dp)            :: d1ij23
  real(dp)            :: d1ij33
  real(dp)            :: d2ij11
  real(dp)            :: d2ij21
  real(dp)            :: d2ij31
  real(dp)            :: d2ij12
  real(dp)            :: d2ij22
  real(dp)            :: d2ij32
  real(dp)            :: d2ij13
  real(dp)            :: d2ij23
  real(dp)            :: d2ij33
  real(dp)            :: dii11
  real(dp)            :: dii21
  real(dp)            :: dii31
  real(dp)            :: dii22
  real(dp)            :: dii32
  real(dp)            :: dii33
  real(dp)            :: diix
  real(dp)            :: diiy
  real(dp)            :: diiz
  real(dp)            :: dijx
  real(dp)            :: dijy
  real(dp)            :: dijz
  real(dp)            :: djj11
  real(dp)            :: djj21
  real(dp)            :: djj31
  real(dp)            :: djj22
  real(dp)            :: djj32
  real(dp)            :: djj33
  real(dp)            :: drix
  real(dp)            :: driy
  real(dp)            :: driz
  real(dp)            :: drjx
  real(dp)            :: drjy
  real(dp)            :: drjz
  real(dp)            :: dt11
  real(dp)            :: dt21
  real(dp)            :: dt31
  real(dp)            :: dt12
  real(dp)            :: dt22
  real(dp)            :: dt32
  real(dp)            :: dt13
  real(dp)            :: dt23
  real(dp)            :: dt33
  real(dp)            :: dt11i
  real(dp)            :: dt21i
  real(dp)            :: dt31i
  real(dp)            :: dt12i
  real(dp)            :: dt22i
  real(dp)            :: dt32i
  real(dp)            :: dt13i
  real(dp)            :: dt23i
  real(dp)            :: dt33i
  real(dp)            :: dt11r
  real(dp)            :: dt21r
  real(dp)            :: dt31r
  real(dp)            :: dt12r
  real(dp)            :: dt22r
  real(dp)            :: dt32r
  real(dp)            :: dt13r
  real(dp)            :: dt23r
  real(dp)            :: dt33r
  real(dp)            :: dw2ds12
  real(dp)            :: dw2ds34
  real(dp)            :: dw2dsii
  real(dp)            :: dw2dsij
  real(dp)            :: dw2dsjj
  real(dp)            :: eiix
  real(dp)            :: eiiy
  real(dp)            :: eiiz
  real(dp)            :: eijx
  real(dp)            :: eijy
  real(dp)            :: eijz
  real(dp)            :: erix
  real(dp)            :: eriy
  real(dp)            :: eriz
  real(dp)            :: erjx
  real(dp)            :: erjy
  real(dp)            :: erjz
  real(dp)            :: rmii
  real(dp)            :: rmij
  real(dp)            :: rmjj
#ifdef TRACE
  call trace_in('projd3dVd')
#endif
!
!*********************************************
!  Project contributions on to phonon modes  *
!*********************************************
  rmij = rmassi*rmassj
  rmii = rmassi*rmassi
  rmjj = rmassj*rmassj
!
  if (i.ne.j) then
!$$$$$$$$$$$$$$$$$$$$$$$$$$$
!  i and j different case  $
!$$$$$$$$$$$$$$$$$$$$$$$$$$$
    if (lcorei.and.lcorej) then
!****************
!  Core - Core  *
!****************
      do k = 1,mcvloc
        erix = eigr(ix,k)
        eriy = eigr(iy,k)
        eriz = eigr(iz,k)
        erjx = eigr(jx,k)
        erjy = eigr(jy,k)
        erjz = eigr(jz,k)
!
!  Change sign of eigenvectors to match those in free energy case
!
        eiix = - eigi(ix,k)
        eiiy = - eigi(iy,k)
        eiiz = - eigi(iz,k)
        eijx = - eigi(jx,k)
        eijy = - eigi(jy,k)
        eijz = - eigi(jz,k)
!%%%%%%%%%%%%%%%%%
!  Off diagonal  %
!%%%%%%%%%%%%%%%%%
!
!  Real-real-real  +  Imag-real-imag
!
        d1ij11 = erix*erjx + eiix*eijx
        d1ij21 = eriy*erjx + eiiy*eijx
        d1ij31 = eriz*erjx + eiiz*eijx
        d1ij12 = erix*erjy + eiix*eijy
        d1ij22 = eriy*erjy + eiiy*eijy
        d1ij32 = eriz*erjy + eiiz*eijy
        d1ij13 = erix*erjz + eiix*eijz
        d1ij23 = eriy*erjz + eiiy*eijz
        d1ij33 = eriz*erjz + eiiz*eijz
!
!  Imag-imag-real - real*imag*imag
!
        d2ij11 = eiix*erjx - erix*eijx
        d2ij21 = eiiy*erjx - eriy*eijx
        d2ij31 = eiiz*erjx - eriz*eijx
        d2ij12 = eiix*erjy - erix*eijy
        d2ij22 = eiiy*erjy - eriy*eijy
        d2ij32 = eiiz*erjy - eriz*eijy
        d2ij13 = eiix*erjz - erix*eijz
        d2ij23 = eiiy*erjz - eriy*eijz
        d2ij33 = eiiz*erjz - eriz*eijz
!%%%%%%%%%%%%%%%%
!  On diagonal  %
!%%%%%%%%%%%%%%%%
!
!  Real-real-real + Imag-real-imag
!
        dii11 = erix*erix + eiix*eiix
        dii21 = erix*eriy + eiix*eiiy
        dii31 = erix*eriz + eiix*eiiz
        dii22 = eriy*eriy + eiiy*eiiy
        dii32 = eriy*eriz + eiiy*eiiz
        dii33 = eriz*eriz + eiiz*eiiz
        djj11 = erjx*erjx + eijx*eijx
        djj21 = erjx*erjy + eijx*eijy
        djj31 = erjx*erjz + eijx*eijz
        djj22 = erjy*erjy + eijy*eijy
        djj32 = erjy*erjz + eijy*eijz
        djj33 = erjz*erjz + eijz*eijz
!
!  Convert to free energy derivatives
!
        do kk = 1,3
          dw2ds12 = d3rs(1,1,kk)*d1ij11 + d3rs(2,1,kk)*d1ij21 +  &
                    d3rs(3,1,kk)*d1ij31 + d3rs(1,2,kk)*d1ij12 +  &
                    d3rs(2,2,kk)*d1ij22 + d3rs(3,2,kk)*d1ij32 +  &
                    d3rs(1,3,kk)*d1ij13 + d3rs(2,3,kk)*d1ij23 +  &
                    d3rs(3,3,kk)*d1ij33
          dw2ds34 = d3is(1,1,kk)*d2ij11 + d3is(2,1,kk)*d2ij21 +  &
                    d3is(3,1,kk)*d2ij31 + d3is(1,2,kk)*d2ij12 +  &
                    d3is(2,2,kk)*d2ij22 + d3is(3,2,kk)*d2ij32 +  &
                    d3is(1,3,kk)*d2ij13 + d3is(2,3,kk)*d2ij23 +  &
                    d3is(3,3,kk)*d2ij33
          dw2dsij = 2.0_dp*(dw2ds12 + dw2ds34)*rmij
          dw2dsii = d3s(1,1,kk)*dii11 + (d3s(2,1,kk) + d3s(1,2,kk))*dii21 +  &
                    (d3s(3,1,kk) + d3s(1,3,kk))*dii31 + d3s(2,2,kk)*dii22 +  &
                    (d3s(3,2,kk) + d3s(2,3,kk))*dii32 + d3s(3,3,kk)*dii33
          dw2dsjj = d3s(1,1,kk)*djj11 + (d3s(2,1,kk) + d3s(1,2,kk))*djj21 +  &
                    (d3s(3,1,kk) + d3s(1,3,kk))*djj31 + d3s(2,2,kk)*djj22 +  &
                    (d3s(3,2,kk) + d3s(2,3,kk))*djj32 + d3s(3,3,kk)*djj33
          grueneisen(k) = grueneisen(k) + (dw2dsij - dw2dsii*rmii - dw2dsjj*rmjj)
        enddo
!
        if (lmany.and.nmanyk.gt.0) then
!*****************
!  Off diagonal  *
!*****************
!
!  Real
!
          dt11r = 2.0_dp*d1ij11*rmij
          dt21r = 2.0_dp*d1ij21*rmij
          dt31r = 2.0_dp*d1ij31*rmij
          dt12r = 2.0_dp*d1ij12*rmij
          dt22r = 2.0_dp*d1ij22*rmij
          dt32r = 2.0_dp*d1ij32*rmij
          dt13r = 2.0_dp*d1ij13*rmij
          dt23r = 2.0_dp*d1ij23*rmij
          dt33r = 2.0_dp*d1ij33*rmij
!
!  Imaginary
!
          dt11i = 2.0_dp*d2ij11*rmij
          dt21i = 2.0_dp*d2ij21*rmij
          dt31i = 2.0_dp*d2ij31*rmij
          dt12i = 2.0_dp*d2ij12*rmij
          dt22i = 2.0_dp*d2ij22*rmij
          dt32i = 2.0_dp*d2ij32*rmij
          dt13i = 2.0_dp*d2ij13*rmij
          dt23i = 2.0_dp*d2ij23*rmij
          dt33i = 2.0_dp*d2ij33*rmij
!
!  Real - on-diagonal
!
          dt11 =  - (dii11*rmii + djj11*rmjj)
          dt21 =  - (dii21*rmii + djj21*rmjj)
          dt31 =  - (dii31*rmii + djj31*rmjj)
          dt12 =  - (dii21*rmii + djj21*rmjj)
          dt22 =  - (dii22*rmii + djj22*rmjj)
          dt32 =  - (dii32*rmii + djj32*rmjj)
          dt13 =  - (dii31*rmii + djj31*rmjj)
          dt23 =  - (dii32*rmii + djj32*rmjj)
          dt33 =  - (dii33*rmii + djj33*rmjj)
!
          call projthbk3dV(nmanyk,d33s,d33rs,d33is, &
            dt11r,dt21r,dt31r,dt12r,dt22r,dt32r,dt13r, &
            dt23r,dt33r,dt11i,dt21i,dt31i,dt12i,dt22i, &
            dt32i,dt13i,dt23i,dt33i,dt11,dt21,dt31,dt12, &
            dt22,dt32,dt13,dt23,dt33,grueneisen(k))
          if (nforkl.gt.0) then
            call projfork3dV(nforkl,d34s,d34rs,d34is, &
              dt11r,dt21r,dt31r,dt12r,dt22r,dt32r,dt13r, &
              dt23r,dt33r,dt11i,dt21i,dt31i,dt12i,dt22i, &
              dt32i,dt13i,dt23i,dt33i,dt11,dt21,dt31,dt12, &
              dt22,dt32,dt13,dt23,dt33,grueneisen(k))
          endif
        endif
      enddo
    elseif (.not.lcorei.and..not.lcorej) then
!******************
!  Shell - Shell  *
!******************
      do k = 1,mcvloc
        drix = derv2(ix,k)
        driy = derv2(iy,k)
        driz = derv2(iz,k)
        drjx = derv2(jx,k)
        drjy = derv2(jy,k)
        drjz = derv2(jz,k)
        diix = dervi(ix,k)
        diiy = dervi(iy,k)
        diiz = dervi(iz,k)
        dijx = dervi(jx,k)
        dijy = dervi(jy,k)
        dijz = dervi(jz,k)
!%%%%%%%%%%%%%%%%%%%
!  Off - diagonal  %
!%%%%%%%%%%%%%%%%%%%
!
!  Real-real-real + Imag-real-imag
!
        d1ij11 = drix*drjx + diix*dijx
        d1ij21 = driy*drjx + diiy*dijx
        d1ij31 = driz*drjx + diiz*dijx
        d1ij12 = drix*drjy + diix*dijy
        d1ij22 = driy*drjy + diiy*dijy
        d1ij32 = driz*drjy + diiz*dijy
        d1ij13 = drix*drjz + diix*dijz
        d1ij23 = driy*drjz + diiy*dijz
        d1ij33 = driz*drjz + diiz*dijz
!
!  Imag-imag-real - real*imag*imag
!
        d2ij11 = diix*drjx - drix*dijx
        d2ij21 = diiy*drjx - driy*dijx
        d2ij31 = diiz*drjx - driz*dijx
        d2ij12 = diix*drjy - drix*dijy
        d2ij22 = diiy*drjy - driy*dijy
        d2ij32 = diiz*drjy - driz*dijy
        d2ij13 = diix*drjz - drix*dijz
        d2ij23 = diiy*drjz - driy*dijz
        d2ij33 = diiz*drjz - driz*dijz
!%%%%%%%%%%%%%%%%
!  On diagonal  %
!%%%%%%%%%%%%%%%%
!
!  Real-real-real + Imag-real-imag
!
        dii11 = drix*drix + diix*diix
        dii21 = drix*driy + diix*diiy
        dii31 = drix*driz + diix*diiz
        dii22 = driy*driy + diiy*diiy
        dii32 = driy*driz + diiy*diiz
        dii33 = driz*driz + diiz*diiz
        djj11 = drjx*drjx + dijx*dijx
        djj21 = drjx*drjy + dijx*dijy
        djj31 = drjx*drjz + dijx*dijz
        djj22 = drjy*drjy + dijy*dijy
        djj32 = drjy*drjz + dijy*dijz
        djj33 = drjz*drjz + dijz*dijz
!
!  Convert to free energy derivatives
!
        do kk = 1,3
          dw2ds12 = d3rs(1,1,kk)*d1ij11 + d3rs(2,1,kk)*d1ij21 +  &
                    d3rs(3,1,kk)*d1ij31 + d3rs(1,2,kk)*d1ij12 +  &
                    d3rs(2,2,kk)*d1ij22 + d3rs(3,2,kk)*d1ij32 +  &
                    d3rs(1,3,kk)*d1ij13 + d3rs(2,3,kk)*d1ij23 +  &
                    d3rs(3,3,kk)*d1ij33
          dw2ds34 = d3is(1,1,kk)*d2ij11 + d3is(2,1,kk)*d2ij21 +  &
                    d3is(3,1,kk)*d2ij31 + d3is(1,2,kk)*d2ij12 +  &
                    d3is(2,2,kk)*d2ij22 + d3is(3,2,kk)*d2ij32 +  &
                    d3is(1,3,kk)*d2ij13 + d3is(2,3,kk)*d2ij23 +  &
                    d3is(3,3,kk)*d2ij33
          dw2dsij = 2.0_dp*(dw2ds12 + dw2ds34)
          dw2dsii = d3s(1,1,kk)*dii11 + (d3s(2,1,kk) + d3s(1,2,kk))*dii21 +  &
                    (d3s(3,1,kk) + d3s(1,3,kk))*dii31 + d3s(2,2,kk)*dii22 +  &
                    (d3s(3,2,kk) + d3s(2,3,kk))*dii32 + d3s(3,3,kk)*dii33
          dw2dsjj = d3s(1,1,kk)*djj11 + (d3s(2,1,kk) + d3s(1,2,kk))*djj21 +  &
                    (d3s(3,1,kk) + d3s(1,3,kk))*djj31 + d3s(2,2,kk)*djj22 +  &
                    (d3s(3,2,kk) + d3s(2,3,kk))*djj32 + d3s(3,3,kk)*djj33
          grueneisen(k) = grueneisen(k) + (dw2dsij - dw2dsii - dw2dsjj)
        enddo
!
        if (lmany.and.nmanyk.gt.0) then
!
!  Real - off-diagonal
!
          dt11r = 2.0_dp*d1ij11
          dt21r = 2.0_dp*d1ij21
          dt31r = 2.0_dp*d1ij31
          dt12r = 2.0_dp*d1ij12
          dt22r = 2.0_dp*d1ij22
          dt32r = 2.0_dp*d1ij32
          dt13r = 2.0_dp*d1ij13
          dt23r = 2.0_dp*d1ij23
          dt33r = 2.0_dp*d1ij33
!
!  Imaginary
!
          dt11i = 2.0_dp*d2ij11
          dt21i = 2.0_dp*d2ij21
          dt31i = 2.0_dp*d2ij31
          dt12i = 2.0_dp*d2ij12
          dt22i = 2.0_dp*d2ij22
          dt32i = 2.0_dp*d2ij32
          dt13i = 2.0_dp*d2ij13
          dt23i = 2.0_dp*d2ij23
          dt33i = 2.0_dp*d2ij33
!
!  Real - on-diagonal
!
          dt11 = - dii11 - djj11
          dt21 = - dii21 - djj21
          dt31 = - dii31 - djj31
          dt12 = - dii21 - djj21
          dt22 = - dii22 - djj22
          dt32 = - dii32 - djj32
          dt13 = - dii31 - djj31
          dt23 = - dii32 - djj32
          dt33 = - dii33 - djj33
!
          call projthbk3dV(nmanyk,d33s,d33rs,d33is, &
            dt11r,dt21r,dt31r,dt12r,dt22r,dt32r,dt13r, &
            dt23r,dt33r,dt11i,dt21i,dt31i,dt12i,dt22i, &
            dt32i,dt13i,dt23i,dt33i,dt11,dt21,dt31,dt12, &
            dt22,dt32,dt13,dt23,dt33,grueneisen(k))
          if (nforkl.gt.0) then
            call projfork3dV(nforkl,d34s,d34rs,d34is, &
              dt11r,dt21r,dt31r,dt12r,dt22r,dt32r,dt13r, &
              dt23r,dt33r,dt11i,dt21i,dt31i,dt12i,dt22i, &
              dt32i,dt13i,dt23i,dt33i,dt11,dt21,dt31,dt12, &
              dt22,dt32,dt13,dt23,dt33,grueneisen(k))
          endif
        endif
      enddo
    elseif (lcorei.and..not.lcorej) then
!*****************
!  Core - Shell  *
!*****************
!
!  As shells are sorted to after cores, this combination should never happen!
!
      call outerror('This combination should never be reached!',0_i4)
      call stopnow('projd3dVd')
    elseif (.not.lcorei.and.lcorej) then
!*****************
!  Shell - Core  *
!*****************
      do k = 1,mcvloc
        erjx = eigr(jx,k)
        erjy = eigr(jy,k)
        erjz = eigr(jz,k)
!
!  Change sign of eigenvectors to match those in free energy case
!
        eijx = - eigi(jx,k)
        eijy = - eigi(jy,k)
        eijz = - eigi(jz,k)
!
        drix = derv2(ix,k)
        driy = derv2(iy,k)
        driz = derv2(iz,k)
        diix = dervi(ix,k)
        diiy = dervi(iy,k)
        diiz = dervi(iz,k)
!%%%%%%%%%%%%%%%%%%%%%%%%%
!  e.(D3).Psn component  %
!%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  Real-real-real + Imag-real-imag
!
        d1ij11 = drix*erjx + diix*eijx
        d1ij21 = driy*erjx + diiy*eijx
        d1ij31 = driz*erjx + diiz*eijx
        d1ij12 = drix*erjy + diix*eijy
        d1ij22 = driy*erjy + diiy*eijy
        d1ij32 = driz*erjy + diiz*eijy
        d1ij13 = drix*erjz + diix*eijz
        d1ij23 = driy*erjz + diiy*eijz
        d1ij33 = driz*erjz + diiz*eijz
!
!  Imag-imag-real - real*imag*imag
!
        d2ij11 = diix*erjx - drix*eijx
        d2ij21 = diiy*erjx - driy*eijx
        d2ij31 = diiz*erjx - driz*eijx
        d2ij12 = diix*erjy - drix*eijy
        d2ij22 = diiy*erjy - driy*eijy
        d2ij32 = diiz*erjy - driz*eijy
        d2ij13 = diix*erjz - drix*eijz
        d2ij23 = diiy*erjz - driy*eijz
        d2ij33 = diiz*erjz - driz*eijz
!%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  Pns.(D3).Psn component  %
!%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  Real-real-real + Imag-real-imag
!
        dii11 = drix*drix + diix*diix
        dii21 = drix*driy + diix*diiy
        dii31 = drix*driz + diix*diiz
        dii22 = driy*driy + diiy*diiy
        dii32 = driy*driz + diiy*diiz
        dii33 = driz*driz + diiz*diiz
        djj11 = erjx*erjx + eijx*eijx
        djj21 = erjx*erjy + eijx*eijy
        djj31 = erjx*erjz + eijx*eijz
        djj22 = erjy*erjy + eijy*eijy
        djj32 = erjy*erjz + eijy*eijz
        djj33 = erjz*erjz + eijz*eijz
!
!  Convert to free energy derivatives
!
        do kk = 1,3
          dw2ds12 = d3rs(1,1,kk)*d1ij11 + d3rs(2,1,kk)*d1ij21 +  &
                    d3rs(3,1,kk)*d1ij31 + d3rs(1,2,kk)*d1ij12 +  &
                    d3rs(2,2,kk)*d1ij22 + d3rs(3,2,kk)*d1ij32 +  &
                    d3rs(1,3,kk)*d1ij13 + d3rs(2,3,kk)*d1ij23 +  &
                    d3rs(3,3,kk)*d1ij33
          dw2ds34 = d3is(1,1,kk)*d2ij11 + d3is(2,1,kk)*d2ij21 +  &
                    d3is(3,1,kk)*d2ij31 + d3is(1,2,kk)*d2ij12 +  &
                    d3is(2,2,kk)*d2ij22 + d3is(3,2,kk)*d2ij32 +  &
                    d3is(1,3,kk)*d2ij13 + d3is(2,3,kk)*d2ij23 +  &
                    d3is(3,3,kk)*d2ij33
          dw2dsij = 2.0_dp*(dw2ds12 + dw2ds34)*rmassj
          dw2dsii = d3s(1,1,kk)*dii11 + (d3s(2,1,kk) + d3s(1,2,kk))*dii21 +  &
                    (d3s(3,1,kk) + d3s(1,3,kk))*dii31 + d3s(2,2,kk)*dii22 +  &
                    (d3s(3,2,kk) + d3s(2,3,kk))*dii32 + d3s(3,3,kk)*dii33
          dw2dsjj = d3s(1,1,kk)*djj11 + (d3s(2,1,kk) + d3s(1,2,kk))*djj21 +  &
                    (d3s(3,1,kk) + d3s(1,3,kk))*djj31 + d3s(2,2,kk)*djj22 +  &
                    (d3s(3,2,kk) + d3s(2,3,kk))*djj32 + d3s(3,3,kk)*djj33
          grueneisen(k) = grueneisen(k) - (dw2dsij + dw2dsii + dw2dsjj*rmjj)
        enddo
!
        if (lmany.and.nmanyk.gt.0) then
!
!  Real - off-diagonal
!
          dt11r = - 2.0_dp*d1ij11*rmassj
          dt21r = - 2.0_dp*d1ij21*rmassj
          dt31r = - 2.0_dp*d1ij31*rmassj
          dt12r = - 2.0_dp*d1ij12*rmassj
          dt22r = - 2.0_dp*d1ij22*rmassj
          dt32r = - 2.0_dp*d1ij32*rmassj
          dt13r = - 2.0_dp*d1ij13*rmassj
          dt23r = - 2.0_dp*d1ij23*rmassj
          dt33r = - 2.0_dp*d1ij33*rmassj
!
!  Imaginary
!
          dt11i = - 2.0_dp*d2ij11*rmassj
          dt21i = - 2.0_dp*d2ij21*rmassj
          dt31i = - 2.0_dp*d2ij31*rmassj
          dt12i = - 2.0_dp*d2ij12*rmassj
          dt22i = - 2.0_dp*d2ij22*rmassj
          dt32i = - 2.0_dp*d2ij32*rmassj
          dt13i = - 2.0_dp*d2ij13*rmassj
          dt23i = - 2.0_dp*d2ij23*rmassj
          dt33i = - 2.0_dp*d2ij33*rmassj
!
!  Real - on-diagonal
!
          dt11 = - (dii11 + djj11*rmjj)
          dt21 = - (dii21 + djj21*rmjj)
          dt31 = - (dii31 + djj31*rmjj)
          dt12 = - (dii21 + djj21*rmjj)
          dt22 = - (dii22 + djj22*rmjj)
          dt32 = - (dii32 + djj32*rmjj)
          dt13 = - (dii31 + djj31*rmjj)
          dt23 = - (dii32 + djj32*rmjj)
          dt33 = - (dii33 + djj33*rmjj)
!
          call projthbk3dV(nmanyk,d33s,d33rs,d33is, &
            dt11r,dt21r,dt31r,dt12r,dt22r,dt32r,dt13r, &
            dt23r,dt33r,dt11i,dt21i,dt31i,dt12i,dt22i, &
            dt32i,dt13i,dt23i,dt33i,dt11,dt21,dt31,dt12, &
            dt22,dt32,dt13,dt23,dt33,grueneisen(k))
          if (nforkl.gt.0) then
            call projfork3dV(nforkl,d34s,d34rs,d34is, &
              dt11r,dt21r,dt31r,dt12r,dt22r,dt32r,dt13r, &
              dt23r,dt33r,dt11i,dt21i,dt31i,dt12i,dt22i, &
              dt32i,dt13i,dt23i,dt33i,dt11,dt21,dt31,dt12, &
              dt22,dt32,dt13,dt23,dt33,grueneisen(k))
          endif
        endif
      enddo
    endif
  elseif (i.eq.j) then
!$$$$$$$$$$$$$$$$$
!  i  =  j case  $
!$$$$$$$$$$$$$$$$$
    if (lcorei.and.lcorej) then
!****************
!  Core - Core  *
!****************
      do k = 1,mcvloc
        erix = eigr(ix,k)
        eriy = eigr(iy,k)
        eriz = eigr(iz,k)
!
!  Change sign of eigenvectors to match those in free energy case
!
        eiix = - eigi(ix,k)
        eiiy = - eigi(iy,k)
        eiiz = - eigi(iz,k)
!%%%%%%%%%%%%%%%%
!  On diagonal  %
!%%%%%%%%%%%%%%%%
!
!  Real-real-real + Imag-real-imag
!
        dii11 = erix*erix + eiix*eiix
        dii21 = erix*eriy + eiix*eiiy
        dii31 = erix*eriz + eiix*eiiz
        dii22 = eriy*eriy + eiiy*eiiy
        dii32 = eriy*eriz + eiiy*eiiz
        dii33 = eriz*eriz + eiiz*eiiz
!
!  Multiply by mass factor
!
        dii11 = dii11*rmii
        dii21 = dii21*rmii
        dii31 = dii31*rmii
        dii22 = dii22*rmii
        dii32 = dii32*rmii
        dii33 = dii33*rmii
!
!  Convert to free energy derivatives
!
        do kk = 1,3
!
!  Take difference of phased on-diagonal element and pure real term
!
          do ii = 1,3
            d3d(1,ii) = d3rs(1,ii,kk) - d3s(1,ii,kk)
            d3d(2,ii) = d3rs(2,ii,kk) - d3s(2,ii,kk)
            d3d(3,ii) = d3rs(3,ii,kk) - d3s(3,ii,kk)
          enddo
          dw2dsii = d3d(1,1)*dii11 + 2.0_dp*d3d(2,1)*dii21 +  &
                    2.0_dp*d3d(3,1)*dii31 + d3d(2,2)*dii22 +  &
                    2.0_dp*d3d(3,2)*dii32 + d3d(3,3)*dii33
          grueneisen(k) = grueneisen(k) + dw2dsii
        enddo
        if (lmany.and.nmanyk.gt.0) then
!
!  Three-body free energy derivatives of diagonal elements
!
          call projthbk3ddV(nmanyk,d33s,d33rs,dii11,dii21,dii31,dii22,dii32,dii33,grueneisen(k))
          if (nforkl.gt.0) then
            call projfork3ddV(nforkl,d34s,d34rs,dii11,dii21,dii31,dii22,dii32,dii33,grueneisen(k))
          endif
        endif
      enddo
    elseif (.not.lcorei.and..not.lcorej) then
!******************
!  Shell - Shell  *
!******************
      do k = 1,mcvloc
        drix = derv2(ix,k)
        driy = derv2(iy,k)
        driz = derv2(iz,k)
        diix = dervi(ix,k)
        diiy = dervi(iy,k)
        diiz = dervi(iz,k)
!%%%%%%%%%%%%%%%%
!  On diagonal  %
!%%%%%%%%%%%%%%%%
!
!  Real-real-real + Imag-real-imag
!
        dii11 = drix*drix + diix*diix
        dii21 = drix*driy + diix*diiy
        dii31 = drix*driz + diix*diiz
        dii22 = driy*driy + diiy*diiy
        dii32 = driy*driz + diiy*diiz
        dii33 = driz*driz + diiz*diiz
!
!  Convert to free energy derivatives
!
        do kk = 1,3
!
!  Take difference of phased on-diagonal element and pure real term
!
          do ii = 1,3
            d3d(1,ii) = d3rs(1,ii,kk) - d3s(1,ii,kk)
            d3d(2,ii) = d3rs(2,ii,kk) - d3s(2,ii,kk)
            d3d(3,ii) = d3rs(3,ii,kk) - d3s(3,ii,kk)
          enddo
          dw2dsii = d3d(1,1)*dii11 + 2.0_dp*d3d(2,1)*dii21 +  &
                    2.0_dp*d3d(3,1)*dii31 + d3d(2,2)*dii22 +  &
                    2.0_dp*d3d(3,2)*dii32 + d3d(3,3)*dii33
          grueneisen(k) = grueneisen(k) + dw2dsii
        enddo
!
        if (lmany.and.nmanyk.gt.0) then
!
!  Three-body free energy derivatives of diagonal elements
!
          call projthbk3ddV(nmanyk,d33s,d33rs,dii11,dii21,dii31,dii22,dii32,dii33,grueneisen(k))
          if (nforkl.gt.0) then
            call projfork3ddV(nforkl,d34s,d34rs,dii11,dii21,dii31,dii22,dii32,dii33,grueneisen(k))
          endif
        endif
      enddo
    endif
  endif
#ifdef TRACE
  call trace_out('projd3dVd')
#endif
!
  return
  end
