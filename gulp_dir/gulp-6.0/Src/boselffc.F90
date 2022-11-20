  subroutine BOselffc
!
!  Calculates the self energy for the bond order charges.
!  Unphased phonon version.
!
!  11/14 Created from BOself
!   2/15 Cell indices now handled
!   7/17 lDoQDeriv flag handling improved
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
  use bondorderdata
  use control,        only : lDoQDeriv2
  use current
  use derivatives,    only : dqdxyz, d2qdxyz2, nqatoms, nqatomptr, nqatomcell
  use derivatives,    only : d2cell, nd2central, nd2cell, nd2cellptr
  use parallel,       only : procid, nprocs
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Local variables
!
  integer(i4)                                    :: i
  integer(i4)                                    :: ij1
  integer(i4)                                    :: ij2
  integer(i4)                                    :: ij3
  integer(i4)                                    :: jk1
  integer(i4)                                    :: jk2
  integer(i4)                                    :: jk3
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
  integer(i4)                                    :: ncindijm
  integer(i4)                                    :: ncindijp
  integer(i4)                                    :: ncindjkm
  integer(i4)                                    :: ncindjkp
  integer(i4)                                    :: ntypi
  logical                                        :: lDoTerm
  real(dp)                                       :: g_cpu_time
  real(dp)                                       :: dedq
  real(dp)                                       :: d2edq2
  real(dp)                                       :: dneqv
  real(dp)                                       :: esum
  real(dp)                                       :: oci
  real(dp)                                       :: qi
  real(dp)                                       :: rho
  real(dp)                                       :: t1
  real(dp)                                       :: t2
!
  if (.not.lDoQDeriv2) return
#ifdef TRACE
  call trace_in('boselffc')
#endif
!
  t1 = g_cpu_time()
!*************************************
!  Calculate self-energy of charges  *
!*************************************
  do i = 1+procid,numat,nprocs
    nati = nat(i)
    ntypi = nftype(i)
    oci = occuf(i)
    dneqv = oci*oci
    qi = qf(i)
    ix = 3*(i - 1) + 1
    iy = ix + 1
    iz = iy + 1
!       
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
          dedq = - esum*rho/(qi - BOq0ref(m))**2
          d2edq2 = - dedq*(rho/(qi - BOq0ref(m)) + 2.0_dp)/(qi - BOq0ref(m))
        endif
!
        if (lDoTerm) then
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
            ij1 = nqatomcell(1,mm,i)
            ij2 = nqatomcell(2,mm,i)
            ij3 = nqatomcell(3,mm,i)
!
            if (abs(ij1).gt.nd2cell(1).or. &
                abs(ij2).gt.nd2cell(2).or. &
                abs(ij3).gt.nd2cell(3)) then
!
!  Find cell index - if outside user range then assign to central cell
!
              ncindijp = nd2central
              ncindijm = nd2central
            else
!
!  Compute index
!
              ncindijp = nd2cellptr(nd2cell(1)+1+ij1,nd2cell(2)+1+ij2,nd2cell(3)+1+ij3)
              ncindijm = nd2cellptr(nd2cell(1)+1-ij1,nd2cell(2)+1-ij2,nd2cell(3)+1-ij3)
            endif
!
!  Second internal derivatives : dE/dQ x d2Q/da.db
!
            if (i.gt.j) then
              d2cell(jx,ix,ncindijp) = d2cell(jx,ix,ncindijp) + dedq*d2qdxyz2(mmxx,i)
              d2cell(jy,ix,ncindijp) = d2cell(jy,ix,ncindijp) + dedq*d2qdxyz2(mmxy,i)
              d2cell(jz,ix,ncindijp) = d2cell(jz,ix,ncindijp) + dedq*d2qdxyz2(mmxz,i)
              d2cell(jx,iy,ncindijp) = d2cell(jx,iy,ncindijp) + dedq*d2qdxyz2(mmxy,i)
              d2cell(jy,iy,ncindijp) = d2cell(jy,iy,ncindijp) + dedq*d2qdxyz2(mmyy,i)
              d2cell(jz,iy,ncindijp) = d2cell(jz,iy,ncindijp) + dedq*d2qdxyz2(mmyz,i)
              d2cell(jx,iz,ncindijp) = d2cell(jx,iz,ncindijp) + dedq*d2qdxyz2(mmxz,i)
              d2cell(jy,iz,ncindijp) = d2cell(jy,iz,ncindijp) + dedq*d2qdxyz2(mmyz,i)
              d2cell(jz,iz,ncindijp) = d2cell(jz,iz,ncindijp) + dedq*d2qdxyz2(mmzz,i)
            elseif (i.lt.j) then
              d2cell(ix,jx,ncindijm) = d2cell(ix,jx,ncindijm) + dedq*d2qdxyz2(mmxx,i)
              d2cell(iy,jx,ncindijm) = d2cell(iy,jx,ncindijm) + dedq*d2qdxyz2(mmxy,i)
              d2cell(iz,jx,ncindijm) = d2cell(iz,jx,ncindijm) + dedq*d2qdxyz2(mmxz,i)
              d2cell(ix,jy,ncindijm) = d2cell(ix,jy,ncindijm) + dedq*d2qdxyz2(mmxy,i)
              d2cell(iy,jy,ncindijm) = d2cell(iy,jy,ncindijm) + dedq*d2qdxyz2(mmyy,i)
              d2cell(iz,jy,ncindijm) = d2cell(iz,jy,ncindijm) + dedq*d2qdxyz2(mmyz,i)
              d2cell(ix,jz,ncindijm) = d2cell(ix,jz,ncindijm) + dedq*d2qdxyz2(mmxz,i)
              d2cell(iy,jz,ncindijm) = d2cell(iy,jz,ncindijm) + dedq*d2qdxyz2(mmyz,i)
              d2cell(iz,jz,ncindijm) = d2cell(iz,jz,ncindijm) + dedq*d2qdxyz2(mmzz,i)
            endif
!
!  Second internal derivatives : d2E/dQ2 x dQ/da x dQ/db : i - j
!
            if (i.gt.j) then
              d2cell(jx,ix,ncindijp) = d2cell(jx,ix,ncindijp) + d2edq2*dqdxyz(ix,i)*dqdxyz(jx,i)
              d2cell(jy,ix,ncindijp) = d2cell(jy,ix,ncindijp) + d2edq2*dqdxyz(ix,i)*dqdxyz(jy,i)
              d2cell(jz,ix,ncindijp) = d2cell(jz,ix,ncindijp) + d2edq2*dqdxyz(ix,i)*dqdxyz(jz,i)
              d2cell(jx,iy,ncindijp) = d2cell(jx,iy,ncindijp) + d2edq2*dqdxyz(iy,i)*dqdxyz(jx,i)
              d2cell(jy,iy,ncindijp) = d2cell(jy,iy,ncindijp) + d2edq2*dqdxyz(iy,i)*dqdxyz(jy,i)
              d2cell(jz,iy,ncindijp) = d2cell(jz,iy,ncindijp) + d2edq2*dqdxyz(iy,i)*dqdxyz(jz,i)
              d2cell(jx,iz,ncindijp) = d2cell(jx,iz,ncindijp) + d2edq2*dqdxyz(iz,i)*dqdxyz(jx,i)
              d2cell(jy,iz,ncindijp) = d2cell(jy,iz,ncindijp) + d2edq2*dqdxyz(iz,i)*dqdxyz(jy,i)
              d2cell(jz,iz,ncindijp) = d2cell(jz,iz,ncindijp) + d2edq2*dqdxyz(iz,i)*dqdxyz(jz,i)
            elseif (i.lt.j) then
              d2cell(ix,jx,ncindijm) = d2cell(ix,jx,ncindijm) + d2edq2*dqdxyz(jx,i)*dqdxyz(ix,i)
              d2cell(iy,jx,ncindijm) = d2cell(iy,jx,ncindijm) + d2edq2*dqdxyz(jx,i)*dqdxyz(iy,i)
              d2cell(iz,jx,ncindijm) = d2cell(iz,jx,ncindijm) + d2edq2*dqdxyz(jx,i)*dqdxyz(iz,i)
              d2cell(ix,jy,ncindijm) = d2cell(ix,jy,ncindijm) + d2edq2*dqdxyz(jy,i)*dqdxyz(ix,i)
              d2cell(iy,jy,ncindijm) = d2cell(iy,jy,ncindijm) + d2edq2*dqdxyz(jy,i)*dqdxyz(iy,i)
              d2cell(iz,jy,ncindijm) = d2cell(iz,jy,ncindijm) + d2edq2*dqdxyz(jy,i)*dqdxyz(iz,i)
              d2cell(ix,jz,ncindijm) = d2cell(ix,jz,ncindijm) + d2edq2*dqdxyz(jz,i)*dqdxyz(ix,i)
              d2cell(iy,jz,ncindijm) = d2cell(iy,jz,ncindijm) + d2edq2*dqdxyz(jz,i)*dqdxyz(iy,i)
              d2cell(iz,jz,ncindijm) = d2cell(iz,jz,ncindijm) + d2edq2*dqdxyz(jz,i)*dqdxyz(iz,i)
            endif
!
            do mn = 1,mm-1
              k = nqatomptr(mn,i)
              kx = 3*(k - 1) + 1
              ky = kx + 1
              kz = ky + 1
!
              jk1 = nqatomcell(1,mn,i) - nqatomcell(1,mm,i)
              jk2 = nqatomcell(2,mn,i) - nqatomcell(2,mm,i)
              jk3 = nqatomcell(3,mn,i) - nqatomcell(3,mm,i)
!
              if (abs(jk1).gt.nd2cell(1).or. &
                  abs(jk2).gt.nd2cell(2).or. &
                  abs(jk3).gt.nd2cell(3)) then
!
!  Find cell index - if outside user range then assign to central cell
!
                ncindjkp = nd2central
                ncindjkm = nd2central
              else
!
!  Compute index
!
                ncindjkp = nd2cellptr(nd2cell(1)+1+jk1,nd2cell(2)+1+jk2,nd2cell(3)+1+jk3)
                ncindjkm = nd2cellptr(nd2cell(1)+1-jk1,nd2cell(2)+1-jk2,nd2cell(3)+1-jk3)
              endif
!
!  Second internal derivatives : d2E/dQ2 x dQ/da x dQ/db : k - j
!
              if (j.gt.k) then
                d2cell(kx,jx,ncindjkp) = d2cell(kx,jx,ncindjkp) + d2edq2*dqdxyz(jx,i)*dqdxyz(kx,i)
                d2cell(ky,jx,ncindjkp) = d2cell(ky,jx,ncindjkp) + d2edq2*dqdxyz(jx,i)*dqdxyz(ky,i)
                d2cell(kz,jx,ncindjkp) = d2cell(kz,jx,ncindjkp) + d2edq2*dqdxyz(jx,i)*dqdxyz(kz,i)
                d2cell(kx,jy,ncindjkp) = d2cell(kx,jy,ncindjkp) + d2edq2*dqdxyz(jy,i)*dqdxyz(kx,i)
                d2cell(ky,jy,ncindjkp) = d2cell(ky,jy,ncindjkp) + d2edq2*dqdxyz(jy,i)*dqdxyz(ky,i)
                d2cell(kz,jy,ncindjkp) = d2cell(kz,jy,ncindjkp) + d2edq2*dqdxyz(jy,i)*dqdxyz(kz,i)
                d2cell(kx,jz,ncindjkp) = d2cell(kx,jz,ncindjkp) + d2edq2*dqdxyz(jz,i)*dqdxyz(kx,i)
                d2cell(ky,jz,ncindjkp) = d2cell(ky,jz,ncindjkp) + d2edq2*dqdxyz(jz,i)*dqdxyz(ky,i)
                d2cell(kz,jz,ncindjkp) = d2cell(kz,jz,ncindjkp) + d2edq2*dqdxyz(jz,i)*dqdxyz(kz,i)
              elseif (j.lt.k) then
                d2cell(jx,kx,ncindjkm) = d2cell(jx,kx,ncindjkm) + d2edq2*dqdxyz(kx,i)*dqdxyz(jx,i)
                d2cell(jy,kx,ncindjkm) = d2cell(jy,kx,ncindjkm) + d2edq2*dqdxyz(kx,i)*dqdxyz(jy,i)
                d2cell(jz,kx,ncindjkm) = d2cell(jz,kx,ncindjkm) + d2edq2*dqdxyz(kx,i)*dqdxyz(jz,i)
                d2cell(jx,ky,ncindjkm) = d2cell(jx,ky,ncindjkm) + d2edq2*dqdxyz(ky,i)*dqdxyz(jx,i)
                d2cell(jy,ky,ncindjkm) = d2cell(jy,ky,ncindjkm) + d2edq2*dqdxyz(ky,i)*dqdxyz(jy,i)
                d2cell(jz,ky,ncindjkm) = d2cell(jz,ky,ncindjkm) + d2edq2*dqdxyz(ky,i)*dqdxyz(jz,i)
                d2cell(jx,kz,ncindjkm) = d2cell(jx,kz,ncindjkm) + d2edq2*dqdxyz(kz,i)*dqdxyz(jx,i)
                d2cell(jy,kz,ncindjkm) = d2cell(jy,kz,ncindjkm) + d2edq2*dqdxyz(kz,i)*dqdxyz(jy,i)
                d2cell(jz,kz,ncindjkm) = d2cell(jz,kz,ncindjkm) + d2edq2*dqdxyz(kz,i)*dqdxyz(jz,i)
              endif
            enddo
          enddo
        endif
      endif
    enddo
  enddo
!
  t2 = g_cpu_time()
  tbondorder = tbondorder + t2 - t1
#ifdef TRACE
  call trace_out('boselffc')
#endif
!
  return
  end
