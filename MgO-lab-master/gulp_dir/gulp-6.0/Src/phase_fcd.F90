  subroutine phase_fcd(xk,yk,zk)
!
!  Calculates the phased dynamical matrix for a k point from
!  the second derivatives stored for supercells.
!  Distributed memory parallel version.
!
!  Note that at present the reciprocal space and 1-D real space
!  Coulomb sums are computed using the standard algorithm.
!
!   4/17 Created from phase_fc
!   2/18 Trace added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
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
!  Copyright Curtin University 2019
!
!  Julian Gale, CIC, Curtin University, September 2019
!
  use configurations, only : lbsmat
  use g_constants
  use control
  use current
  use datatypes
  use derivatives
  use element
  use kspace
  use molecule
  use parallel
  use shells
  use symmetry
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  real(dp),    intent(in)                      :: xk
  real(dp),    intent(in)                      :: yk
  real(dp),    intent(in)                      :: zk
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ib
  integer(i4)                                  :: iblast
  integer(i4)                                  :: ii
  integer(i4)                                  :: iis
  integer(i4)                                  :: iisg
  integer(i4)                                  :: iloc
  integer(i4)                                  :: ix
  integer(i4)                                  :: iy
  integer(i4)                                  :: iz
  integer(i4)                                  :: ixg
  integer(i4)                                  :: iyg
  integer(i4)                                  :: izg
  integer(i4)                                  :: j
  integer(i4)                                  :: jj
  integer(i4)                                  :: jjs
  integer(i4)                                  :: jx
  integer(i4)                                  :: jy
  integer(i4)                                  :: jz
  integer(i4)                                  :: kk
  integer(i4)                                  :: maxlim
  integer(i4)                                  :: maxlimloc
  integer(i4)                                  :: mint
  integer(i4)                                  :: mintloc
  integer(i4)                                  :: nbsi
  integer(i4)                                  :: nbsig
  integer(i4)                                  :: nbsj
  integer(i4)                                  :: ncind
  logical                                      :: lbsi
  logical                                      :: lbsj
  real(dp)                                     :: cosk
  real(dp)                                     :: oneij
  real(dp)                                     :: sink
  real(dp)                                     :: kvf(3,3)
  real(dp)                                     :: rvx
  real(dp)                                     :: rvy
  real(dp)                                     :: rvz
  real(dp)                                     :: xi
  real(dp)                                     :: yi
  real(dp)                                     :: zi
  real(dp)                                     :: xj
  real(dp)                                     :: yj
  real(dp)                                     :: zj
  real(dp)                                     :: xkv
  real(dp)                                     :: ykv
  real(dp)                                     :: zkv
#ifdef TRACE
  call trace_in('phase_fcd')
#endif
!****************************
!  Zero second derivatives  *
!****************************
  mint = 3*numat
  mintloc = 3*natomsonnode
  maxlim = mint
  maxlimloc = mintloc
  if (nbsmat.gt.0) then
    maxlim = maxlim + numat
    maxlimloc = maxlimloc + natomsonnode
  endif
  if (maxlimloc.gt.maxd2u) then
    maxd2u = maxlimloc
    call changemaxd2
  endif
  if (maxlim.gt.maxd2) then
    maxd2 = maxlim
    call changemaxd2
  endif
  do i = 1,maxlimloc
    do j = 1,maxlim
      derv2(j,i) = 0.0_dp
      dervi(j,i) = 0.0_dp
    enddo
  enddo
!
!  If group velocities are to be computed then zero k derivatives of dynamical matrix
!
  if (lgroupvelocity) then
    derv2dk(1:3,1:maxlim,1:maxlimloc) = 0.0_dpc
  endif
!
!  Select appropriate K vectors
!
  if (lkfull.and.ndim.eq.3) then
    call kvector3Df(kvf)
  else
    kvf(1:3,1:3) = kv(1:3,1:3)
  endif
!***************************
!  Calculate phase factor  *
!***************************
  if (ndim.eq.3) then
    xkv = xk*kvf(1,1) + yk*kvf(1,2) + zk*kvf(1,3)
    ykv = xk*kvf(2,1) + yk*kvf(2,2) + zk*kvf(2,3)
    zkv = xk*kvf(3,1) + yk*kvf(3,2) + zk*kvf(3,3)
  elseif (ndim.eq.2) then
    xkv = xk*kvf(1,1) + yk*kvf(1,2)
    ykv = xk*kvf(2,1) + yk*kvf(2,2)
    zkv = 0.0_dp
  elseif (ndim.eq.1) then
    xkv = xk*kvf(1,1)
    ykv = 0.0_dp
    zkv = 0.0_dp
  endif
!***************************************************************************
!  Reciprocal space component - for now has to be calculated analytically  *
!***************************************************************************
  if (lewald.and.ndim.gt.1) then
    call kindex
    if (ndim.eq.3) then
      call recip3Dpd(xkv,ykv,zkv)
    elseif (ndim.eq.2) then
      call recip2Dpd(xkv,ykv)
    endif
  endif
!*************************************************
!  Extra real space component for 1-D case only  *
!*************************************************
  if (ndim.eq.1) call real1Dpd(xkv)
!********************
!  Loop over cells  *
!********************
  ncind = 0
  do ii = -nd2cell(1),nd2cell(1)
    do jj = -nd2cell(2),nd2cell(2)
      do kk = -nd2cell(3),nd2cell(3)
        ncind = ncind + 1
        rvx = dble(ii)*rv(1,1) + dble(jj)*rv(1,2) + dble(kk)*rv(1,3)
        rvy = dble(ii)*rv(2,1) + dble(jj)*rv(2,2) + dble(kk)*rv(2,3)
        rvz = dble(ii)*rv(3,1) + dble(jj)*rv(3,2) + dble(kk)*rv(3,3)
!
!  Loop over first atom
!
        ix = -2
        iy = -1
        iz =  0
        nbsi = 0
        iblast = 0
        do iloc = 1,natomsonnode
          i = node2atom(iloc)
          xi = xclat(i)
          yi = yclat(i)
          zi = zclat(i)
          ix = ix + 3
          iy = iy + 3
          iz = iz + 3
          ixg = 3*(i-1) + 1
          iyg = ixg + 1
          izg = ixg + 2
          if (lbsmat(nsft+nrelf2a(i))) then
            lbsi = .true.
            nbsi = nbsi + 1
            iis = 3*natomsonnode + nbsi
!
!  Compute nbsi for global equivalent
!
            do ib = iblast+1,i
              if (lbsmat(nsft+nrelf2a(ib))) then
                nbsig = nbsig + 1
              endif
            enddo
            iblast = i
            iisg = 3*numat + nbsig
          else
            lbsi = .false.
          endif
!
          if (lbsi) then
!
!  Radial components
!
            derv2(iisg,iis) = derv2(iisg,iis) + d2cell(iisg,iis,ncind)
!
            derv2(ixg,iis) = derv2(ixg,iis) + d2cell(ixg,iis,ncind)
            derv2(iyg,iis) = derv2(iyg,iis) + d2cell(iyg,iis,ncind)
            derv2(izg,iis) = derv2(izg,iis) + d2cell(izg,iis,ncind)
!
            derv2(iisg,ix) = derv2(iisg,ix) + d2cell(ixg,iis,ncind)
            derv2(iisg,iy) = derv2(iisg,iy) + d2cell(iyg,iis,ncind)
            derv2(iisg,iz) = derv2(iisg,iz) + d2cell(izg,iis,ncind)
          endif
!
!  Loop over second atom
!
          jx = -2
          jy = -1
          jz =  0
          nbsj = 0
          do j = 1,numat
            xj = xclat(j) + rvx
            yj = yclat(j) + rvy
            zj = zclat(j) + rvz
            jx = jx + 3
            jy = jy + 3
            jz = jz + 3
            if (lbsmat(nsft+nrelf2a(j))) then
              lbsj = .true.
              nbsj = nbsj + 1
              jjs = 3*numat + nbsj
            else
              lbsj = .false.
            endif
!
!  Compute pair-wise phase factor
!
            if (i.eq.j) then
              oneij = 1.0_dp
            else
              oneij = 0.0_dp
            endif
            cosk = xkv*(xj-xi) + ykv*(yj-yi) + zkv*(zj-zi)
            sink = sin(cosk)
            cosk = cos(cosk) - oneij
!
!  Add phased contribution to the main arrays
!
            derv2(jx,ix) = derv2(jx,ix) + d2cell(jx,ix,ncind)*cosk
            derv2(jy,ix) = derv2(jy,ix) + d2cell(jy,ix,ncind)*cosk
            derv2(jz,ix) = derv2(jz,ix) + d2cell(jz,ix,ncind)*cosk
            derv2(jx,iy) = derv2(jx,iy) + d2cell(jx,iy,ncind)*cosk
            derv2(jy,iy) = derv2(jy,iy) + d2cell(jy,iy,ncind)*cosk
            derv2(jz,iy) = derv2(jz,iy) + d2cell(jz,iy,ncind)*cosk
            derv2(jx,iz) = derv2(jx,iz) + d2cell(jx,iz,ncind)*cosk
            derv2(jy,iz) = derv2(jy,iz) + d2cell(jy,iz,ncind)*cosk
            derv2(jz,iz) = derv2(jz,iz) + d2cell(jz,iz,ncind)*cosk
!
            dervi(jx,ix) = dervi(jx,ix) + d2cell(jx,ix,ncind)*sink
            dervi(jy,ix) = dervi(jy,ix) + d2cell(jy,ix,ncind)*sink
            dervi(jz,ix) = dervi(jz,ix) + d2cell(jz,ix,ncind)*sink
            dervi(jx,iy) = dervi(jx,iy) + d2cell(jx,iy,ncind)*sink
            dervi(jy,iy) = dervi(jy,iy) + d2cell(jy,iy,ncind)*sink
            dervi(jz,iy) = dervi(jz,iy) + d2cell(jz,iy,ncind)*sink
            dervi(jx,iz) = dervi(jx,iz) + d2cell(jx,iz,ncind)*sink
            dervi(jy,iz) = dervi(jy,iz) + d2cell(jy,iz,ncind)*sink
            dervi(jz,iz) = dervi(jz,iz) + d2cell(jz,iz,ncind)*sink
!
!  Radial components
!
            if (lbsi) then
              derv2(jx,iis) = derv2(jx,iis) + d2cell(jx,iis,ncind)*cosk
              derv2(jy,iis) = derv2(jy,iis) + d2cell(jy,iis,ncind)*cosk
              derv2(jz,iis) = derv2(jz,iis) + d2cell(jz,iis,ncind)*cosk
!
              dervi(jx,iis) = dervi(jx,iis) + d2cell(jx,iis,ncind)*sink
              dervi(jy,iis) = dervi(jy,iis) + d2cell(jy,iis,ncind)*sink
              dervi(jz,iis) = dervi(jz,iis) + d2cell(jz,iis,ncind)*sink
            endif
            if (lbsj) then
              derv2(jjs,ix) = derv2(jjs,ix) + d2cell(jjs,ix,ncind)*cosk
              derv2(jjs,iy) = derv2(jjs,iy) + d2cell(jjs,iy,ncind)*cosk
              derv2(jjs,iz) = derv2(jjs,iz) + d2cell(jjs,iz,ncind)*cosk
!
              dervi(jjs,ix) = dervi(jjs,ix) + d2cell(jjs,ix,ncind)*sink
              dervi(jjs,iy) = dervi(jjs,iy) + d2cell(jjs,iy,ncind)*sink
              dervi(jjs,iz) = dervi(jjs,iz) + d2cell(jjs,iz,ncind)*sink
!
              if (lbsi) then
                derv2(jjs,iis) = derv2(jjs,iis) + d2cell(jjs,iis,ncind)*cosk
                dervi(jjs,iis) = dervi(jjs,iis) + d2cell(jjs,iis,ncind)*sink
              endif
            endif
!
!  End of loop over second atom
!
          enddo
!
!  End of loop over first atom
!
        enddo
!
!  End of loop over cells
!
      enddo
    enddo
  enddo
#ifdef TRACE
  call trace_out('phase_fcd')
#endif
!
  return
  end
