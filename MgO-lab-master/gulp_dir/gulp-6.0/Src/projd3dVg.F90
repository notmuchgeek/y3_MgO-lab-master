  subroutine projd3dVg(mcv,d3s,i,ix,iy,iz,j,jx,jy,jz,lcorei,lcorej, &
                       rmassi,rmassj,derv2,eigr,maxd2,grueneisen,d33s, &
                       lmany,nmanyk,d34s,nforkl)
!
!  Projects the 3rd derivatives on to the phonon modes to compute volume derivatives
!  Used in the calculation of Grueneisen parameters
!
!   1/18 Created from projd3a
!   2/18 Trace added
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
  integer(i4)         :: maxd2
  integer(i4)         :: mcv
  integer(i4)         :: nmanyk
  integer(i4)         :: nforkl
  logical             :: lcorei
  logical             :: lcorej
  logical             :: lmany
  real(dp)            :: d3s(3,3,3)
  real(dp)            :: d33s(108,*)
  real(dp)            :: d34s(54,*)
  real(dp)            :: derv2(maxd2,*)
  real(dp)            :: eigr(maxd2,*)
  real(dp)            :: grueneisen(*)
  real(dp)            :: rmassi
  real(dp)            :: rmassj
!
!  Local variables
!
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
  real(dp)            :: dii11
  real(dp)            :: dii21
  real(dp)            :: dii31
  real(dp)            :: dii22
  real(dp)            :: dii32
  real(dp)            :: dii33
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
  real(dp)            :: dw2dsii
  real(dp)            :: dw2dsij
  real(dp)            :: dw2dsjj
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
  call trace_in('projd3dVg')
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
      do k = 1,mcv
        erix = eigr(ix,k)
        eriy = eigr(iy,k)
        eriz = eigr(iz,k)
        erjx = eigr(jx,k)
        erjy = eigr(jy,k)
        erjz = eigr(jz,k)
!%%%%%%%%%%%%%%%%%
!  Off diagonal  %
!%%%%%%%%%%%%%%%%%
!
!  Real-real-real
!
        d1ij11 = erix*erjx
        d1ij21 = eriy*erjx
        d1ij31 = eriz*erjx
        d1ij12 = erix*erjy
        d1ij22 = eriy*erjy
        d1ij32 = eriz*erjy
        d1ij13 = erix*erjz
        d1ij23 = eriy*erjz
        d1ij33 = eriz*erjz
!%%%%%%%%%%%%%%%%
!  On diagonal  %
!%%%%%%%%%%%%%%%%
!
!  Real-real-real
!
        dii11 = erix*erix
        dii21 = erix*eriy
        dii31 = erix*eriz
        dii22 = eriy*eriy
        dii32 = eriy*eriz
        dii33 = eriz*eriz
        djj11 = erjx*erjx
        djj21 = erjx*erjy
        djj31 = erjx*erjz
        djj22 = erjy*erjy
        djj32 = erjy*erjz
        djj33 = erjz*erjz
!
!  Convert to free energy derivatives
!
        do kk = 1,3
          dw2ds12 = d3s(1,1,kk)*d1ij11 + d3s(2,1,kk)*d1ij21 +  &
                    d3s(3,1,kk)*d1ij31 + d3s(1,2,kk)*d1ij12 +  &
                    d3s(2,2,kk)*d1ij22 + d3s(3,2,kk)*d1ij32 +  &
                    d3s(1,3,kk)*d1ij13 + d3s(2,3,kk)*d1ij23 +  &
                    d3s(3,3,kk)*d1ij33
          dw2dsij = 2.0_dp*dw2ds12*rmij
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
          call projthbk3dVg(nmanyk,d33s, &
            dt11r,dt21r,dt31r,dt12r,dt22r,dt32r,dt13r, &
            dt23r,dt33r,dt11,dt21,dt31,dt12,dt22,dt32, &
            dt13,dt23,dt33,grueneisen(k))
          if (nforkl.gt.0) then
            call projfork3dVg(nforkl,d34s, &
              dt11r,dt21r,dt31r,dt12r,dt22r,dt32r,dt13r, &
              dt23r,dt33r,dt11,dt21,dt31,dt12,dt22,dt32, &
              dt13,dt23,dt33,grueneisen(k))
          endif
        endif
      enddo
    elseif (.not.lcorei.and..not.lcorej) then
!******************
!  Shell - Shell  *
!******************
      do k = 1,mcv
        drix = derv2(k,ix)
        driy = derv2(k,iy)
        driz = derv2(k,iz)
        drjx = derv2(k,jx)
        drjy = derv2(k,jy)
        drjz = derv2(k,jz)
!%%%%%%%%%%%%%%%%%%%
!  Off - diagonal  %
!%%%%%%%%%%%%%%%%%%%
!
!  Real-real-real
!
        d1ij11 = drix*drjx
        d1ij21 = driy*drjx
        d1ij31 = driz*drjx
        d1ij12 = drix*drjy
        d1ij22 = driy*drjy
        d1ij32 = driz*drjy
        d1ij13 = drix*drjz
        d1ij23 = driy*drjz
        d1ij33 = driz*drjz
!%%%%%%%%%%%%%%%%
!  On diagonal  %
!%%%%%%%%%%%%%%%%
!
!  Real-real-real
!
        dii11 = drix*drix
        dii21 = drix*driy
        dii31 = drix*driz
        dii22 = driy*driy
        dii32 = driy*driz
        dii33 = driz*driz
        djj11 = drjx*drjx
        djj21 = drjx*drjy
        djj31 = drjx*drjz
        djj22 = drjy*drjy
        djj32 = drjy*drjz
        djj33 = drjz*drjz
!
!  Convert to free energy derivatives
!
        do kk = 1,3
          dw2ds12 = d3s(1,1,kk)*d1ij11 + d3s(2,1,kk)*d1ij21 +  &
                    d3s(3,1,kk)*d1ij31 + d3s(1,2,kk)*d1ij12 +  &
                    d3s(2,2,kk)*d1ij22 + d3s(3,2,kk)*d1ij32 +  &
                    d3s(1,3,kk)*d1ij13 + d3s(2,3,kk)*d1ij23 +  &
                    d3s(3,3,kk)*d1ij33
          dw2dsij = 2.0_dp*dw2ds12
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
          call projthbk3dVg(nmanyk,d33s, &
            dt11r,dt21r,dt31r,dt12r,dt22r,dt32r,dt13r, &
            dt23r,dt33r,dt11,dt21,dt31,dt12,dt22,dt32, &
            dt13,dt23,dt33,grueneisen(k))
          if (nforkl.gt.0) then
            call projfork3dVg(nforkl,d34s, &
              dt11r,dt21r,dt31r,dt12r,dt22r,dt32r,dt13r, &
              dt23r,dt33r,dt11,dt21,dt31,dt12,dt22,dt32, &
              dt13,dt23,dt33,grueneisen(k))
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
      call stopnow('projd3dVg')
    elseif (.not.lcorei.and.lcorej) then
!*****************
!  Shell - Core  *
!*****************
      do k = 1,mcv
        erjx = eigr(jx,k)
        erjy = eigr(jy,k)
        erjz = eigr(jz,k)
        drix = derv2(k,ix)
        driy = derv2(k,iy)
        driz = derv2(k,iz)
!%%%%%%%%%%%%%%%%%%%%%%%%%
!  e.(D3).Psn component  %
!%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  Real-real-real
!
        d1ij11 = drix*erjx
        d1ij21 = driy*erjx
        d1ij31 = driz*erjx
        d1ij12 = drix*erjy
        d1ij22 = driy*erjy
        d1ij32 = driz*erjy
        d1ij13 = drix*erjz
        d1ij23 = driy*erjz
        d1ij33 = driz*erjz
!%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  Pns.(D3).Psn component  %
!%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  Real-real-real + Imag-real-imag
!
        dii11 = drix*drix
        dii21 = drix*driy
        dii31 = drix*driz
        dii22 = driy*driy
        dii32 = driy*driz
        dii33 = driz*driz
        djj11 = erjx*erjx
        djj21 = erjx*erjy
        djj31 = erjx*erjz
        djj22 = erjy*erjy
        djj32 = erjy*erjz
        djj33 = erjz*erjz
!
!  Convert to free energy derivatives
!
        do kk = 1,3
          dw2ds12 = d3s(1,1,kk)*d1ij11 + d3s(2,1,kk)*d1ij21 +  &
                    d3s(3,1,kk)*d1ij31 + d3s(1,2,kk)*d1ij12 +  &
                    d3s(2,2,kk)*d1ij22 + d3s(3,2,kk)*d1ij32 +  &
                    d3s(1,3,kk)*d1ij13 + d3s(2,3,kk)*d1ij23 +  &
                    d3s(3,3,kk)*d1ij33
          dw2dsij = 2.0_dp*dw2ds12*rmassj
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
          call projthbk3dVg(nmanyk,d33s, &
            dt11r,dt21r,dt31r,dt12r,dt22r,dt32r,dt13r, &
            dt23r,dt33r,dt11,dt21,dt31,dt12,dt22,dt32, &
            dt13,dt23,dt33,grueneisen(k))
          if (nforkl.gt.0) then
            call projfork3dVg(nforkl,d34s, &
              dt11r,dt21r,dt31r,dt12r,dt22r,dt32r,dt13r, &
              dt23r,dt33r,dt11,dt21,dt31,dt12,dt22,dt32, &
              dt13,dt23,dt33,grueneisen(k))
          endif
        endif
      enddo
    endif
  endif
#ifdef TRACE
  call trace_out('projd3dVg')
#endif
!
  return
  end
