  subroutine BOselfp(xkv,ykv,zkv)
!
!  Calculates the phonon contribution from the self energy 
!  for the bond order charges.
!
!   2/15 Created from BOself
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
!  Julian Gale, CIC, Curtin University, February 2018
!
  use datatypes
  use bondorderdata
  use control,        only : lDoQDeriv1, lDoQDeriv2
  use current
  use derivatives,    only : dqdxyz, d2qdxyz2, nqatoms, nqatomptr
  use derivatives,    only : derv2, dervi, qatomxyz
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
!
!  Passed variables
!
  real(dp),   intent(in)                         :: xkv
  real(dp),   intent(in)                         :: ykv
  real(dp),   intent(in)                         :: zkv
!
!  Local variables
!
  integer(i4)                                    :: i
  integer(i4)                                    :: ix
  integer(i4)                                    :: iy
  integer(i4)                                    :: iz
  integer(i4)                                    :: j
  integer(i4)                                    :: jx
  integer(i4)                                    :: jy
  integer(i4)                                    :: jz
  integer(i4)                                    :: k
  integer(i4)                                    :: kx
  integer(i4)                                    :: ky
  integer(i4)                                    :: kz
  integer(i4)                                    :: m
  integer(i4)                                    :: mm
  integer(i4)                                    :: mn
  integer(i4)                                    :: mx
  integer(i4)                                    :: my
  integer(i4)                                    :: mz
  integer(i4)                                    :: mmxx
  integer(i4)                                    :: mmxy
  integer(i4)                                    :: mmxz
  integer(i4)                                    :: mmyy
  integer(i4)                                    :: mmyz
  integer(i4)                                    :: mmzz
  integer(i4)                                    :: nati
  integer(i4)                                    :: ntypi
  logical                                        :: lDoTerm
  real(dp)                                       :: g_cpu_time
  real(dp)                                       :: dedq
  real(dp)                                       :: d2edq2
  real(dp)                                       :: dneqv
  real(dp)                                       :: cosij
  real(dp)                                       :: cosjk
  real(dp)                                       :: sinij
  real(dp)                                       :: sinjk
  real(dp)                                       :: esum
  real(dp)                                       :: oci
  real(dp)                                       :: one
  real(dp)                                       :: qi
  real(dp)                                       :: rho
  real(dp)                                       :: t1
  real(dp)                                       :: t2
  real(dp)                                       :: xjk
  real(dp)                                       :: yjk
  real(dp)                                       :: zjk
#ifdef TRACE
  call trace_in('boselfp')
#endif
!
  t1 = g_cpu_time()
!
  do i = 1,numat
    nati = nat(i)
    ntypi = nftype(i)
    oci = occuf(i)
    dneqv = oci*oci
    qi = qf(i)
    ix = 3*(i - 1) + 1
    iy = ix + 1
    iz = iy + 1
    do m = 1,nboQ0
      esum = 0.0_dp
      if (nati.eq.nBOspecQ0(m).and.(ntypi.eq.nBOtypQ0(m).or.nBOtypQ0(m).eq.0)) then
        lDoTerm = .false.
        rho = BOq0rho(m)
        if (BOq0ref(m).gt.0.0_dp.and.qi.gt.BOq0ref(m)) then
!
!  Cation case
!
          lDoTerm = .true.
          esum = dneqv*BOq0pot(m)*exp(-rho/(qi - BOq0ref(m)))
          dedq = esum*rho/(qi - BOq0ref(m))**2
          d2edq2 = dedq*(rho/(qi - BOq0ref(m)) - 2.0_dp)/(qi - BOq0ref(m))
        elseif (BOq0ref(m).lt.0.0_dp.and.qi.lt.BOq0ref(m)) then
!
!  Anion case
!
          lDoTerm = .true.
          esum = dneqv*BOq0pot(m)*exp(rho/(qi - BOq0ref(m)))
          if (lDoQDeriv1) then
            dedq = - esum*rho/(qi - BOq0ref(m))**2
            d2edq2 = - dedq*(rho/(qi - BOq0ref(m)) + 2.0_dp)/(qi - BOq0ref(m))
          endif
        endif
!
        if (lDoTerm) then
          if (lDoQDeriv2) then
            do mm = 1,nqatoms(i)
              j = nqatomptr(mm,i)
              jx = 3*(j - 1) + 1
              jy = jx + 1
              jz = jy + 1
!
              mx = 3*(mm - 1) + 1
              my = mx + 1 
              mz = my + 1
!   
              mmxx = mx*(mx + 1)/2 
              mmyy = my*(my + 1)/2 
              mmzz = mz*(mz + 1)/2 
              mmxy = mmyy - 1
              mmxz = mmzz - 2
              mmyz = mmzz - 1
!
!  Calculate phase factor 
!
              if (i.eq.j) then
                one = 1.0_dp
              else
                one = 0.0_dp
              endif
!
              cosij = xkv*qatomxyz(1,mm,i) + ykv*qatomxyz(2,mm,i) + zkv*qatomxyz(3,mm,i)
              sinij = sin(cosij)
              cosij = cos(cosij) - one
!
!  Second internal derivatives : dE/dQ x d2Q/da.db
!
              derv2(jx,ix) = derv2(jx,ix) + dedq*d2qdxyz2(mmxx,i)*cosij
              derv2(jy,ix) = derv2(jy,ix) + dedq*d2qdxyz2(mmxy,i)*cosij
              derv2(jz,ix) = derv2(jz,ix) + dedq*d2qdxyz2(mmxz,i)*cosij
              derv2(jx,iy) = derv2(jx,iy) + dedq*d2qdxyz2(mmxy,i)*cosij
              derv2(jy,iy) = derv2(jy,iy) + dedq*d2qdxyz2(mmyy,i)*cosij
              derv2(jz,iy) = derv2(jz,iy) + dedq*d2qdxyz2(mmyz,i)*cosij
              derv2(jx,iz) = derv2(jx,iz) + dedq*d2qdxyz2(mmxz,i)*cosij
              derv2(jy,iz) = derv2(jy,iz) + dedq*d2qdxyz2(mmyz,i)*cosij
              derv2(jz,iz) = derv2(jz,iz) + dedq*d2qdxyz2(mmzz,i)*cosij
!
              dervi(jx,ix) = dervi(jx,ix) + dedq*d2qdxyz2(mmxx,i)*sinij
              dervi(jy,ix) = dervi(jy,ix) + dedq*d2qdxyz2(mmxy,i)*sinij
              dervi(jz,ix) = dervi(jz,ix) + dedq*d2qdxyz2(mmxz,i)*sinij
              dervi(jx,iy) = dervi(jx,iy) + dedq*d2qdxyz2(mmxy,i)*sinij
              dervi(jy,iy) = dervi(jy,iy) + dedq*d2qdxyz2(mmyy,i)*sinij
              dervi(jz,iy) = dervi(jz,iy) + dedq*d2qdxyz2(mmyz,i)*sinij
              dervi(jx,iz) = dervi(jx,iz) + dedq*d2qdxyz2(mmxz,i)*sinij
              dervi(jy,iz) = dervi(jy,iz) + dedq*d2qdxyz2(mmyz,i)*sinij
              dervi(jz,iz) = dervi(jz,iz) + dedq*d2qdxyz2(mmzz,i)*sinij
!
              if (i.ne.j) then
                derv2(ix,jx) = derv2(ix,jx) + dedq*d2qdxyz2(mmxx,i)*cosij
                derv2(iy,jx) = derv2(iy,jx) + dedq*d2qdxyz2(mmxy,i)*cosij
                derv2(iz,jx) = derv2(iz,jx) + dedq*d2qdxyz2(mmxz,i)*cosij
                derv2(ix,jy) = derv2(ix,jy) + dedq*d2qdxyz2(mmxy,i)*cosij
                derv2(iy,jy) = derv2(iy,jy) + dedq*d2qdxyz2(mmyy,i)*cosij
                derv2(iz,jy) = derv2(iz,jy) + dedq*d2qdxyz2(mmyz,i)*cosij
                derv2(ix,jz) = derv2(ix,jz) + dedq*d2qdxyz2(mmxz,i)*cosij
                derv2(iy,jz) = derv2(iy,jz) + dedq*d2qdxyz2(mmyz,i)*cosij
                derv2(iz,jz) = derv2(iz,jz) + dedq*d2qdxyz2(mmzz,i)*cosij
!
                dervi(ix,jx) = dervi(ix,jx) - dedq*d2qdxyz2(mmxx,i)*sinij
                dervi(iy,jx) = dervi(iy,jx) - dedq*d2qdxyz2(mmxy,i)*sinij
                dervi(iz,jx) = dervi(iz,jx) - dedq*d2qdxyz2(mmxz,i)*sinij
                dervi(ix,jy) = dervi(ix,jy) - dedq*d2qdxyz2(mmxy,i)*sinij
                dervi(iy,jy) = dervi(iy,jy) - dedq*d2qdxyz2(mmyy,i)*sinij
                dervi(iz,jy) = dervi(iz,jy) - dedq*d2qdxyz2(mmyz,i)*sinij
                dervi(ix,jz) = dervi(ix,jz) - dedq*d2qdxyz2(mmxz,i)*sinij
                dervi(iy,jz) = dervi(iy,jz) - dedq*d2qdxyz2(mmyz,i)*sinij
                dervi(iz,jz) = dervi(iz,jz) - dedq*d2qdxyz2(mmzz,i)*sinij
              endif
!
!  Second internal derivatives : d2E/dQ2 x dQ/da x dQ/db : i - j
!
              derv2(jx,ix) = derv2(jx,ix) + d2edq2*dqdxyz(ix,i)*dqdxyz(jx,i)*cosij
              derv2(jy,ix) = derv2(jy,ix) + d2edq2*dqdxyz(ix,i)*dqdxyz(jy,i)*cosij
              derv2(jz,ix) = derv2(jz,ix) + d2edq2*dqdxyz(ix,i)*dqdxyz(jz,i)*cosij
              derv2(jx,iy) = derv2(jx,iy) + d2edq2*dqdxyz(iy,i)*dqdxyz(jx,i)*cosij
              derv2(jy,iy) = derv2(jy,iy) + d2edq2*dqdxyz(iy,i)*dqdxyz(jy,i)*cosij
              derv2(jz,iy) = derv2(jz,iy) + d2edq2*dqdxyz(iy,i)*dqdxyz(jz,i)*cosij
              derv2(jx,iz) = derv2(jx,iz) + d2edq2*dqdxyz(iz,i)*dqdxyz(jx,i)*cosij
              derv2(jy,iz) = derv2(jy,iz) + d2edq2*dqdxyz(iz,i)*dqdxyz(jy,i)*cosij
              derv2(jz,iz) = derv2(jz,iz) + d2edq2*dqdxyz(iz,i)*dqdxyz(jz,i)*cosij
!
              dervi(jx,ix) = dervi(jx,ix) + d2edq2*dqdxyz(ix,i)*dqdxyz(jx,i)*sinij
              dervi(jy,ix) = dervi(jy,ix) + d2edq2*dqdxyz(ix,i)*dqdxyz(jy,i)*sinij
              dervi(jz,ix) = dervi(jz,ix) + d2edq2*dqdxyz(ix,i)*dqdxyz(jz,i)*sinij
              dervi(jx,iy) = dervi(jx,iy) + d2edq2*dqdxyz(iy,i)*dqdxyz(jx,i)*sinij
              dervi(jy,iy) = dervi(jy,iy) + d2edq2*dqdxyz(iy,i)*dqdxyz(jy,i)*sinij
              dervi(jz,iy) = dervi(jz,iy) + d2edq2*dqdxyz(iy,i)*dqdxyz(jz,i)*sinij
              dervi(jx,iz) = dervi(jx,iz) + d2edq2*dqdxyz(iz,i)*dqdxyz(jx,i)*sinij
              dervi(jy,iz) = dervi(jy,iz) + d2edq2*dqdxyz(iz,i)*dqdxyz(jy,i)*sinij
              dervi(jz,iz) = dervi(jz,iz) + d2edq2*dqdxyz(iz,i)*dqdxyz(jz,i)*sinij
!
              if (i.ne.j) then
                derv2(ix,jx) = derv2(ix,jx) + d2edq2*dqdxyz(jx,i)*dqdxyz(ix,i)*cosij
                derv2(iy,jx) = derv2(iy,jx) + d2edq2*dqdxyz(jx,i)*dqdxyz(iy,i)*cosij
                derv2(iz,jx) = derv2(iz,jx) + d2edq2*dqdxyz(jx,i)*dqdxyz(iz,i)*cosij
                derv2(ix,jy) = derv2(ix,jy) + d2edq2*dqdxyz(jy,i)*dqdxyz(ix,i)*cosij
                derv2(iy,jy) = derv2(iy,jy) + d2edq2*dqdxyz(jy,i)*dqdxyz(iy,i)*cosij
                derv2(iz,jy) = derv2(iz,jy) + d2edq2*dqdxyz(jy,i)*dqdxyz(iz,i)*cosij
                derv2(ix,jz) = derv2(ix,jz) + d2edq2*dqdxyz(jz,i)*dqdxyz(ix,i)*cosij
                derv2(iy,jz) = derv2(iy,jz) + d2edq2*dqdxyz(jz,i)*dqdxyz(iy,i)*cosij
                derv2(iz,jz) = derv2(iz,jz) + d2edq2*dqdxyz(jz,i)*dqdxyz(iz,i)*cosij
!
                dervi(ix,jx) = dervi(ix,jx) - d2edq2*dqdxyz(jx,i)*dqdxyz(ix,i)*sinij
                dervi(iy,jx) = dervi(iy,jx) - d2edq2*dqdxyz(jx,i)*dqdxyz(iy,i)*sinij
                dervi(iz,jx) = dervi(iz,jx) - d2edq2*dqdxyz(jx,i)*dqdxyz(iz,i)*sinij
                dervi(ix,jy) = dervi(ix,jy) - d2edq2*dqdxyz(jy,i)*dqdxyz(ix,i)*sinij
                dervi(iy,jy) = dervi(iy,jy) - d2edq2*dqdxyz(jy,i)*dqdxyz(iy,i)*sinij
                dervi(iz,jy) = dervi(iz,jy) - d2edq2*dqdxyz(jy,i)*dqdxyz(iz,i)*sinij
                dervi(ix,jz) = dervi(ix,jz) - d2edq2*dqdxyz(jz,i)*dqdxyz(ix,i)*sinij
                dervi(iy,jz) = dervi(iy,jz) - d2edq2*dqdxyz(jz,i)*dqdxyz(iy,i)*sinij
                dervi(iz,jz) = dervi(iz,jz) - d2edq2*dqdxyz(jz,i)*dqdxyz(iz,i)*sinij
              endif
!
              do mn = 1,mm-1
                k = nqatomptr(mn,i)
                kx = 3*(k - 1) + 1
                ky = kx + 1
                kz = ky + 1
!
!  Calculate phase factor
!
                if (j.eq.k) then
                  one = 1.0_dp
                else
                  one = 0.0_dp
                endif
!
                xjk = qatomxyz(1,mn,i) - qatomxyz(1,mm,i)
                yjk = qatomxyz(2,mn,i) - qatomxyz(2,mm,i)
                zjk = qatomxyz(3,mn,i) - qatomxyz(3,mm,i)
!
                cosjk = xkv*xjk + ykv*yjk + zkv*zjk
                sinjk = sin(cosjk)
                cosjk = cos(cosjk) - one
!
!  Second internal derivatives : d2E/dQ2 x dQ/da x dQ/db : k - j
!
                derv2(kx,jx) = derv2(kx,jx) + d2edq2*dqdxyz(jx,i)*dqdxyz(kx,i)*cosjk
                derv2(ky,jx) = derv2(ky,jx) + d2edq2*dqdxyz(jx,i)*dqdxyz(ky,i)*cosjk
                derv2(kz,jx) = derv2(kz,jx) + d2edq2*dqdxyz(jx,i)*dqdxyz(kz,i)*cosjk
                derv2(kx,jy) = derv2(kx,jy) + d2edq2*dqdxyz(jy,i)*dqdxyz(kx,i)*cosjk
                derv2(ky,jy) = derv2(ky,jy) + d2edq2*dqdxyz(jy,i)*dqdxyz(ky,i)*cosjk
                derv2(kz,jy) = derv2(kz,jy) + d2edq2*dqdxyz(jy,i)*dqdxyz(kz,i)*cosjk
                derv2(kx,jz) = derv2(kx,jz) + d2edq2*dqdxyz(jz,i)*dqdxyz(kx,i)*cosjk
                derv2(ky,jz) = derv2(ky,jz) + d2edq2*dqdxyz(jz,i)*dqdxyz(ky,i)*cosjk
                derv2(kz,jz) = derv2(kz,jz) + d2edq2*dqdxyz(jz,i)*dqdxyz(kz,i)*cosjk
!
                dervi(kx,jx) = dervi(kx,jx) + d2edq2*dqdxyz(jx,i)*dqdxyz(kx,i)*sinjk
                dervi(ky,jx) = dervi(ky,jx) + d2edq2*dqdxyz(jx,i)*dqdxyz(ky,i)*sinjk
                dervi(kz,jx) = dervi(kz,jx) + d2edq2*dqdxyz(jx,i)*dqdxyz(kz,i)*sinjk
                dervi(kx,jy) = dervi(kx,jy) + d2edq2*dqdxyz(jy,i)*dqdxyz(kx,i)*sinjk
                dervi(ky,jy) = dervi(ky,jy) + d2edq2*dqdxyz(jy,i)*dqdxyz(ky,i)*sinjk
                dervi(kz,jy) = dervi(kz,jy) + d2edq2*dqdxyz(jy,i)*dqdxyz(kz,i)*sinjk
                dervi(kx,jz) = dervi(kx,jz) + d2edq2*dqdxyz(jz,i)*dqdxyz(kx,i)*sinjk
                dervi(ky,jz) = dervi(ky,jz) + d2edq2*dqdxyz(jz,i)*dqdxyz(ky,i)*sinjk
                dervi(kz,jz) = dervi(kz,jz) + d2edq2*dqdxyz(jz,i)*dqdxyz(kz,i)*sinjk
!
                if (j.ne.k) then
                  derv2(jx,kx) = derv2(jx,kx) + d2edq2*dqdxyz(kx,i)*dqdxyz(jx,i)*cosjk
                  derv2(jy,kx) = derv2(jy,kx) + d2edq2*dqdxyz(kx,i)*dqdxyz(jy,i)*cosjk
                  derv2(jz,kx) = derv2(jz,kx) + d2edq2*dqdxyz(kx,i)*dqdxyz(jz,i)*cosjk
                  derv2(jx,ky) = derv2(jx,ky) + d2edq2*dqdxyz(ky,i)*dqdxyz(jx,i)*cosjk
                  derv2(jy,ky) = derv2(jy,ky) + d2edq2*dqdxyz(ky,i)*dqdxyz(jy,i)*cosjk
                  derv2(jz,ky) = derv2(jz,ky) + d2edq2*dqdxyz(ky,i)*dqdxyz(jz,i)*cosjk
                  derv2(jx,kz) = derv2(jx,kz) + d2edq2*dqdxyz(kz,i)*dqdxyz(jx,i)*cosjk
                  derv2(jy,kz) = derv2(jy,kz) + d2edq2*dqdxyz(kz,i)*dqdxyz(jy,i)*cosjk
                  derv2(jz,kz) = derv2(jz,kz) + d2edq2*dqdxyz(kz,i)*dqdxyz(jz,i)*cosjk
!
                  dervi(jx,kx) = dervi(jx,kx) - d2edq2*dqdxyz(kx,i)*dqdxyz(jx,i)*sinjk
                  dervi(jy,kx) = dervi(jy,kx) - d2edq2*dqdxyz(kx,i)*dqdxyz(jy,i)*sinjk
                  dervi(jz,kx) = dervi(jz,kx) - d2edq2*dqdxyz(kx,i)*dqdxyz(jz,i)*sinjk
                  dervi(jx,ky) = dervi(jx,ky) - d2edq2*dqdxyz(ky,i)*dqdxyz(jx,i)*sinjk
                  dervi(jy,ky) = dervi(jy,ky) - d2edq2*dqdxyz(ky,i)*dqdxyz(jy,i)*sinjk
                  dervi(jz,ky) = dervi(jz,ky) - d2edq2*dqdxyz(ky,i)*dqdxyz(jz,i)*sinjk
                  dervi(jx,kz) = dervi(jx,kz) - d2edq2*dqdxyz(kz,i)*dqdxyz(jx,i)*sinjk
                  dervi(jy,kz) = dervi(jy,kz) - d2edq2*dqdxyz(kz,i)*dqdxyz(jy,i)*sinjk
                  dervi(jz,kz) = dervi(jz,kz) - d2edq2*dqdxyz(kz,i)*dqdxyz(jz,i)*sinjk
                endif
              enddo
            enddo
          endif
        endif
      endif
    enddo
  enddo
!
  t2 = g_cpu_time()
  tbondorder = tbondorder + t2 - t1
#ifdef TRACE
  call trace_out('boselfp')
#endif
!
  return
  end
