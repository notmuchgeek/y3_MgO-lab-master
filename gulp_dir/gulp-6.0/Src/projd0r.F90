  subroutine projd0r(d3,i,ix,iy,iz,j,jx,jy,jz,lcorei,lcorej,d2dx,d2dy,d2dz,qD,maxqD, &
                     lmany,nmanyk,nptrmanyk,d33,nforkl,nptrfork,nptrforl,d34)
!
!  Computes the 3rd derivatives with respect to atom positions for Raman susceptibilities
!
!   9/13 Created from projd0
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
  use datatypes
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4)       :: i
  integer(i4)       :: ix
  integer(i4)       :: iy
  integer(i4)       :: iz
  integer(i4)       :: j
  integer(i4)       :: jx
  integer(i4)       :: jy
  integer(i4)       :: jz
  integer(i4)       :: maxqD
  integer(i4)       :: nforkl
  integer(i4)       :: nmanyk
  integer(i4)       :: nptrfork(*)
  integer(i4)       :: nptrforl(*)
  integer(i4)       :: nptrmanyk(*)
  logical           :: lcorei
  logical           :: lcorej
  logical           :: lmany
  real(dp)          :: d3(3,3,3)
  real(dp)          :: d33(54,*)
  real(dp)          :: d34(27,*)
  real(dp)          :: qD(maxqD,3)
  real(dp)          :: d2dx(3,3)
  real(dp)          :: d2dy(3,3)
  real(dp)          :: d2dz(3,3)
!
!  Local variables
!
  integer(i4)       :: ia
  integer(i4)       :: ib
  real(dp)          :: d1ij11
  real(dp)          :: d1ij21
  real(dp)          :: d1ij31
  real(dp)          :: d1ij12
  real(dp)          :: d1ij22
  real(dp)          :: d1ij32
  real(dp)          :: d1ij13
  real(dp)          :: d1ij23
  real(dp)          :: d1ij33
  real(dp)          :: dii11
  real(dp)          :: dii21
  real(dp)          :: dii31
  real(dp)          :: dii22
  real(dp)          :: dii32
  real(dp)          :: dii33
  real(dp)          :: djj11
  real(dp)          :: djj21
  real(dp)          :: djj31
  real(dp)          :: djj22
  real(dp)          :: djj32
  real(dp)          :: djj33
  real(dp)          :: davg
  real(dp)          :: drixa
  real(dp)          :: driya
  real(dp)          :: driza
  real(dp)          :: drixb
  real(dp)          :: driyb
  real(dp)          :: drizb
  real(dp)          :: drjxa
  real(dp)          :: drjya
  real(dp)          :: drjza
  real(dp)          :: drjxb
  real(dp)          :: drjyb
  real(dp)          :: drjzb
  real(dp)          :: dt11
  real(dp)          :: dt21
  real(dp)          :: dt31
  real(dp)          :: dt12
  real(dp)          :: dt22
  real(dp)          :: dt32
  real(dp)          :: dt13
  real(dp)          :: dt23
  real(dp)          :: dt33
  real(dp)          :: dw2ds12x
  real(dp)          :: dw2ds12y
  real(dp)          :: dw2ds12z
  real(dp)          :: dw2dsiix
  real(dp)          :: dw2dsiiy
  real(dp)          :: dw2dsiiz
  real(dp)          :: dw2dsijx
  real(dp)          :: dw2dsijy
  real(dp)          :: dw2dsijz
  real(dp)          :: dw2dsjjx
  real(dp)          :: dw2dsjjy
  real(dp)          :: dw2dsjjz
!
  d2dx(1:3,1:3) = 0.0_dp
  d2dy(1:3,1:3) = 0.0_dp
  d2dz(1:3,1:3) = 0.0_dp
!
!  If i = j or both are cores then skip
!
  if ((lcorei.and.lcorej).or.i.eq.j) return
  if (.not.lcorei.and..not.lcorej) then
#ifdef TRACE
  call trace_in('prod0r')
#endif
!******************
!  Shell - Shell  *
!******************
!
!  Loop over Cartesian components
!
    do ia = 1,3
      do ib = 1,3
        drixa = qD(ix,ia)
        driya = qD(iy,ia)
        driza = qD(iz,ia)
        drixb = qD(ix,ib)
        driyb = qD(iy,ib)
        drizb = qD(iz,ib)
        drjxa = qD(jx,ia)
        drjya = qD(jy,ia)
        drjza = qD(jz,ia)
        drjxb = qD(jx,ib)
        drjyb = qD(jy,ib)
        drjzb = qD(jz,ib)
!%%%%%%%%%%%%%%%%%%%
!  Off - diagonal  %
!%%%%%%%%%%%%%%%%%%%
        d1ij11 = drixa*drjxb
        d1ij21 = driya*drjxb
        d1ij31 = driza*drjxb
        d1ij12 = drixa*drjyb
        d1ij22 = driya*drjyb
        d1ij32 = driza*drjyb
        d1ij13 = drixa*drjzb
        d1ij23 = driya*drjzb
        d1ij33 = driza*drjzb
!%%%%%%%%%%%%%%%%
!  On diagonal  %
!%%%%%%%%%%%%%%%%
        dii11 = drixa*drixb
        dii21 = drixa*driyb
        dii31 = drixa*drizb
        dii22 = driya*driyb
        dii32 = driya*drizb
        dii33 = driza*drizb
        djj11 = drjxb*drjxa
        djj21 = drjxb*drjya
        djj31 = drjxb*drjza
        djj22 = drjyb*drjya
        djj32 = drjyb*drjza
        djj33 = drjzb*drjza
!
!  Convert to Raman susceptibility derivatives
!
        dw2ds12x = d3(1,1,1)*d1ij11 + d3(2,1,1)*d1ij21 +  &
                   d3(3,1,1)*d1ij31 + d3(1,2,1)*d1ij12 +  &
                   d3(2,2,1)*d1ij22 + d3(3,2,1)*d1ij32 +  &
                   d3(1,3,1)*d1ij13 + d3(2,3,1)*d1ij23 +  &
                   d3(3,3,1)*d1ij33
        dw2ds12y = d3(1,1,2)*d1ij11 + d3(2,1,2)*d1ij21 +  &
                   d3(3,1,2)*d1ij31 + d3(1,2,2)*d1ij12 +  &
                   d3(2,2,2)*d1ij22 + d3(3,2,2)*d1ij32 +  &
                   d3(1,3,2)*d1ij13 + d3(2,3,2)*d1ij23 +  &
                   d3(3,3,2)*d1ij33
        dw2ds12z = d3(1,1,3)*d1ij11 + d3(2,1,3)*d1ij21 +  &
                   d3(3,1,3)*d1ij31 + d3(1,2,3)*d1ij12 +  &
                   d3(2,2,3)*d1ij22 + d3(3,2,3)*d1ij32 +  &
                   d3(1,3,3)*d1ij13 + d3(2,3,3)*d1ij23 +  &
                   d3(3,3,3)*d1ij33
!
        dw2dsijx = dw2ds12x
        dw2dsijy = dw2ds12y
        dw2dsijz = dw2ds12z
!
        dw2dsiix = d3(1,1,1)*dii11 + (d3(2,1,1) + d3(1,2,1))*dii21 +  &
                  (d3(3,1,1) + d3(1,3,1))*dii31 + d3(2,2,1)*dii22 +  &
                  (d3(3,2,1) + d3(2,3,1))*dii32 + d3(3,3,1)*dii33
        dw2dsiiy = d3(1,1,2)*dii11 + (d3(2,1,2) + d3(1,2,2))*dii21 +  &
                  (d3(3,1,2) + d3(1,3,2))*dii31 + d3(2,2,2)*dii22 +  &
                  (d3(3,2,2) + d3(2,3,2))*dii32 + d3(3,3,2)*dii33
        dw2dsiiz = d3(1,1,3)*dii11 + (d3(2,1,3) + d3(1,2,3))*dii21 +  &
                  (d3(3,1,3) + d3(1,3,3))*dii31 + d3(2,2,3)*dii22 +  &
                  (d3(3,2,3) + d3(2,3,3))*dii32 + d3(3,3,3)*dii33
        dw2dsjjx = d3(1,1,1)*djj11 + (d3(2,1,1) + d3(1,2,1))*djj21 +  &
                  (d3(3,1,1) + d3(1,3,1))*djj31 + d3(2,2,1)*djj22 +  &
                  (d3(3,2,1) + d3(2,3,1))*djj32 + d3(3,3,1)*djj33
        dw2dsjjy = d3(1,1,2)*djj11 + (d3(2,1,2) + d3(1,2,2))*djj21 +  &
                  (d3(3,1,2) + d3(1,3,2))*djj31 + d3(2,2,2)*djj22 +  &
                  (d3(3,2,2) + d3(2,3,2))*djj32 + d3(3,3,2)*djj33
        dw2dsjjz = d3(1,1,3)*djj11 + (d3(2,1,3) + d3(1,2,3))*djj21 +  &
                  (d3(3,1,3) + d3(1,3,3))*djj31 + d3(2,2,3)*djj22 +  &
                  (d3(3,2,3) + d3(2,3,3))*djj32 + d3(3,3,3)*djj33
!
        d2dx(ib,ia) = (dw2dsijx - 0.5*(dw2dsiix + dw2dsjjx))
        d2dy(ib,ia) = (dw2dsijy - 0.5*(dw2dsiiy + dw2dsjjy))
        d2dz(ib,ia) = (dw2dsijz - 0.5*(dw2dsiiz + dw2dsjjz))
!
        if (lmany.and.nmanyk.gt.0) then
          dt11 = - (d1ij11 - 0.5_dp*(dii11 + djj11))
          dt21 = - (d1ij21 - 0.5_dp*(dii21 + djj21))
          dt31 = - (d1ij31 - 0.5_dp*(dii31 + djj31))
          dt12 = - (d1ij12 - 0.5_dp*(dii21 + djj21))
          dt22 = - (d1ij22 - 0.5_dp*(dii22 + djj22))
          dt32 = - (d1ij32 - 0.5_dp*(dii32 + djj32))
          dt13 = - (d1ij13 - 0.5_dp*(dii31 + djj31))
          dt23 = - (d1ij23 - 0.5_dp*(dii32 + djj32))
          dt33 = - (d1ij33 - 0.5_dp*(dii33 + djj33))
          call projthbk0r(i,j,ia,ib,nmanyk,nptrmanyk,d33,dt11,dt21,dt31,dt12,dt22,dt32,dt13,dt23,dt33)
          if (nforkl.gt.0) then
            call projfork0r(ia,ib,nforkl,nptrfork,nptrforl,d34,dt11,dt21,dt31,dt12,dt22,dt32,dt13,dt23,dt33)
          endif
        endif
      enddo
    enddo
  elseif (.not.lcorei.and.lcorej) then
!*****************
!  Shell - core  *
!*****************
!
!  Loop over Cartesian components
!
    do ia = 1,3
      do ib = 1,3
        drixa = qD(ix,ia)
        driya = qD(iy,ia)
        driza = qD(iz,ia)
        drixb = qD(ix,ib)
        driyb = qD(iy,ib)
        drizb = qD(iz,ib)
!%%%%%%%%%%%%%%%%
!  On diagonal  %
!%%%%%%%%%%%%%%%%
        d1ij11 = drixa*drixb
        d1ij21 = driya*drixb
        d1ij31 = driza*drixb
        d1ij12 = drixa*driyb
        d1ij22 = driya*driyb
        d1ij32 = driza*driyb
        d1ij13 = drixa*drizb
        d1ij23 = driya*drizb
        d1ij33 = driza*drizb
!
!  Convert to Raman susceptibility derivatives
!
        dw2dsijx = d3(1,1,1)*d1ij11 + d3(2,1,1)*d1ij21 +  &
                   d3(3,1,1)*d1ij31 + d3(1,2,1)*d1ij12 +  &
                   d3(2,2,1)*d1ij22 + d3(3,2,1)*d1ij32 +  &
                   d3(1,3,1)*d1ij13 + d3(2,3,1)*d1ij23 +  &
                   d3(3,3,1)*d1ij33
        dw2dsijy = d3(1,1,2)*d1ij11 + d3(2,1,2)*d1ij21 +  &
                   d3(3,1,2)*d1ij31 + d3(1,2,2)*d1ij12 +  &
                   d3(2,2,2)*d1ij22 + d3(3,2,2)*d1ij32 +  &
                   d3(1,3,2)*d1ij13 + d3(2,3,2)*d1ij23 +  &
                   d3(3,3,2)*d1ij33
        dw2dsijz = d3(1,1,3)*d1ij11 + d3(2,1,3)*d1ij21 +  &
                   d3(3,1,3)*d1ij31 + d3(1,2,3)*d1ij12 +  &
                   d3(2,2,3)*d1ij22 + d3(3,2,3)*d1ij32 +  &
                   d3(1,3,3)*d1ij13 + d3(2,3,3)*d1ij23 +  &
                   d3(3,3,3)*d1ij33
!
!  Negative sign because this the derivative of the on-diagonal sum
!
        d2dx(ib,ia) = - dw2dsijx
        d2dy(ib,ia) = - dw2dsijy
        d2dz(ib,ia) = - dw2dsijz
!
        if (lmany.and.nmanyk.gt.0) then
          dt11 = d1ij11
          dt21 = d1ij21
          dt31 = d1ij31
          dt12 = d1ij12
          dt22 = d1ij22
          dt32 = d1ij32
          dt13 = d1ij13
          dt23 = d1ij23
          dt33 = d1ij33
          call projthbk0r(i,j,ia,ib,nmanyk,nptrmanyk,d33,dt11,dt21,dt31,dt12,dt22,dt32,dt13,dt23,dt33)
          if (nforkl.gt.0) then
            call projfork0r(ia,ib,nforkl,nptrfork,nptrforl,d34,dt11,dt21,dt31,dt12,dt22,dt32,dt13,dt23,dt33)
          endif
        endif
      enddo
    enddo
  endif
!
!  Average off diagonals to restore symmetry. This compensates for the asymmetric atom loops
!  combined with the many-body terms. 
!
  do ia = 2,3
    do ib = 1,ia
      davg = 0.5_dp*(d2dx(ia,ib) + d2dx(ib,ia))
      d2dx(ia,ib) = davg
      d2dx(ib,ia) = davg
      davg = 0.5_dp*(d2dy(ia,ib) + d2dy(ib,ia))
      d2dy(ia,ib) = davg
      d2dy(ib,ia) = davg
      davg = 0.5_dp*(d2dz(ia,ib) + d2dz(ib,ia))
      d2dz(ia,ib) = davg
      d2dz(ib,ia) = davg
    enddo
  enddo
#ifdef TRACE
  call trace_out('prod0r')
#endif
!
  return
  end
