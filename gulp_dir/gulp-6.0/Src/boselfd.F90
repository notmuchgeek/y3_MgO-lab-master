  subroutine BOselfd(eboQself,lgrad1,lgrad2,lPhonon)
!
!  Calculates the self energy for the bond order charges.
!
!  Distributed memory parallel version.
!
!  On entry :
!
!  lgrad1       = if .true. calculate first derivatives
!  lgrad2       = if .true. calculate second derivatives
!  lPhonon      = if .true. then exclude strain derivatives
!
!  On exit :
!
!  eboQself     = the self energy of the bond order charges
!
!   7/12 Created from boself
!   9/16 cputime renamed to g_cpu_time
!   1/17 Parallelisation added
!   2/17 Modifications added for distributed storage of 
!        q arrays
!   7/17 Correction to addressing of nqatoms made
!   7/17 Adding of esum to eregion2region moved to match boself
!   7/17 Atomic stress added
!   7/17 Site energy added
!   7/17 lDoQDeriv2 calls removed for now
!   2/18 Trace added
!  10/18 Handling of terms for strfin removed since these are no longer needed
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
  use datatypes
  use bondorderdata
  use configurations, only : nregions, nregionno
  use control,        only : lDoQDeriv1, lDoQDeriv2, latomicstress
  use current
  use derivatives,    only : dqdxyz, d2qdxyz2, d2qdxyzs, dqds, d2qds2, nqatoms, nqatomptr
  use derivatives,    only : xdrv, ydrv, zdrv, rstrd, derv2, derv3, sderv2
  use derivatives,    only : xregdrv, yregdrv, zregdrv, atomicstress
  use energies,       only : esregion2, eregion2region, siteenergy
  use optimisation,   only : lfreeze, lopf
  use parallel,       only : natomsonnode, node2atom
  use symmetry,       only : lstr
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  real(dp),    intent(out)                       :: eboQself
  logical,     intent(in)                        :: lgrad1
  logical,     intent(in)                        :: lgrad2
  logical,     intent(in)                        :: lPhonon
!
!  Local variables
!
  integer(i4)                                    :: i
  integer(i4)                                    :: iloc
  integer(i4)                                    :: ind
  integer(i4)                                    :: ix
  integer(i4)                                    :: iy
  integer(i4)                                    :: iz
  integer(i4)                                    :: j
  integer(i4)                                    :: jx
  integer(i4)                                    :: jy
  integer(i4)                                    :: jz
  integer(i4)                                    :: kk
  integer(i4)                                    :: kl
  integer(i4)                                    :: m
  integer(i4)                                    :: mm
  integer(i4)                                    :: mx
  integer(i4)                                    :: my
  integer(i4)                                    :: mz
  integer(i4)                                    :: mmxx
  integer(i4)                                    :: mmxy
  integer(i4)                                    :: mmxz
  integer(i4)                                    :: mmyy
  integer(i4)                                    :: mmyz
  integer(i4)                                    :: mmzz
  integer(i4)                                    :: n
  integer(i4)                                    :: nati
  integer(i4)                                    :: nregioni
  integer(i4)                                    :: nregionj
  integer(i4)                                    :: ntypi
  logical                                        :: lDoTerm
  logical                                        :: lopi
  logical                                        :: lopj
  logical                                        :: lreg2pair
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
#ifdef TRACE
  call trace_in('boselfd')
#endif
!
  t1 = g_cpu_time()
!
!  Initialise Bond Order charge self energy
!
  eboQself = 0.0_dp
!*************************************
!  Calculate self-energy of charges  *
!  Outer loop over atoms on node     *
!*************************************
  do iloc = 1,natomsonnode
    i = node2atom(iloc)
    nati = nat(i)
    ntypi = nftype(i)
    nregioni = nregionno(nsft+nrelf2a(i))
    oci = occuf(i)
    lopi = (.not.lfreeze.or.lopf(nrelf2a(i)))
    dneqv = oci*oci
    qi = qf(i)
    ix = 3*(iloc - 1) + 1
    iy = ix + 1
    iz = iy + 1
!       
!  Set region 2 pair flag
!       
    lreg2pair = .false.
    if (nregions(ncf).ge.2) then
      lreg2pair = (nregionno(nsft+nrelf2a(i)).gt.1)
    endif
    do m = 1,nboQ0
      esum = 0.0_dp
      if (nati.eq.nBOspecQ0(m).and.(ntypi.eq.nBOtypQ0(m).or.nBOtypQ0(m).eq.0)) then
        lDoTerm = .false.
        rho = BOq0rho(m)
        if (BOq0ref(m).gt.0.0d0.and.qi.gt.BOq0ref(m)) then
!
!  Cation case
!
          lDoTerm = .true.
          esum = dneqv*BOq0pot(m)*exp(-rho/(qi - BOq0ref(m)))
          if (lgrad1) then
            dedq = esum*rho/(qi - BOq0ref(m))**2
            if (lgrad2) then
              d2edq2 = dedq*(rho/(qi - BOq0ref(m)) - 2.0_dp)/(qi - BOq0ref(m))
            endif
          endif
        elseif (BOq0ref(m).lt.0.0d0.and.qi.lt.BOq0ref(m)) then
!
!  Anion case
!
          lDoTerm = .true.
          esum = dneqv*BOq0pot(m)*exp(rho/(qi - BOq0ref(m)))
          if (lgrad1.and.lDoQDeriv1) then
            dedq = - esum*rho/(qi - BOq0ref(m))**2
            if (lgrad2) then
              d2edq2 = - dedq*(rho/(qi - BOq0ref(m)) + 2.0_dp)/(qi - BOq0ref(m))
            endif
          endif
        endif
        eregion2region(nregioni,nregioni) = eregion2region(nregioni,nregioni) + esum
!
        siteenergy(i) = siteenergy(i) + esum
!
!  If no derivatives are needed don't do terms
!
        if (.not.lgrad1) lDoTerm = .false.
!
        if (lDoTerm) then
          if (lstr.and..not.lPhonon) then
!
!  Strain first derivatives
!
            do kl = 1,nstrains
              rstrd(kl) = rstrd(kl) + dedq*dqds(kl,iloc)
            enddo
            if (latomicstress) then
              do kl = 1,nstrains
                atomicstress(kl,i) = atomicstress(kl,i) + dedq*dqds(kl,iloc)
              enddo
            endif
            if (lgrad2) then
!
!  Strain-strain second derivatives : d2E/dQ2 x dQ/ds1 x dQ/ds2
!
              do kk = 1,nstrains
                do kl = 1,nstrains
                  sderv2(kl,kk) = sderv2(kl,kk) + d2edq2*dqds(kk,iloc)*dqds(kl,iloc)
                enddo
              enddo
!
!  Strain-strain second derivatives : dE/dQ x d2Q/ds1.ds2
!
              ind = 0
              do kk = 1,nstrains
                do kl = 1,kk-1
                  ind = ind + 1
                  sderv2(kl,kk) = sderv2(kl,kk) + dedq*d2qds2(ind,iloc)
                  sderv2(kk,kl) = sderv2(kk,kl) + dedq*d2qds2(ind,iloc)
                enddo
                ind = ind + 1
                sderv2(kk,kk) = sderv2(kk,kk) + dedq*d2qds2(ind,iloc)
              enddo
            endif
          endif
          do n = 1,nqatoms(iloc)
            j = nqatomptr(n,iloc)
            jx = 3*(j - 1) + 1
            jy = jx + 1
            jz = jy + 1
            nregionj = nregionno(nsft+nrelf2a(j))
            lopj = (.not.lfreeze.or.lopf(nrelf2a(j)))
!
!  Internal first derivatives
!
            if (lopi) then
              xdrv(i) = xdrv(i) - dedq*dqdxyz(jx,iloc)
              ydrv(i) = ydrv(i) - dedq*dqdxyz(jy,iloc)
              zdrv(i) = zdrv(i) - dedq*dqdxyz(jz,iloc)
            endif
            xregdrv(nregioni) = xregdrv(nregioni) - dedq*dqdxyz(jx,iloc)
            yregdrv(nregioni) = yregdrv(nregioni) - dedq*dqdxyz(jy,iloc)
            zregdrv(nregioni) = zregdrv(nregioni) - dedq*dqdxyz(jz,iloc)
            if (nrela2f(nrelf2a(j)).eq.j) then
              if (lopj) then
                xdrv(j) = xdrv(j) + dedq*dqdxyz(jx,iloc)
                ydrv(j) = ydrv(j) + dedq*dqdxyz(jy,iloc)
                zdrv(j) = zdrv(j) + dedq*dqdxyz(jz,iloc)
              endif
              xregdrv(nregionj) = xregdrv(nregionj) + dedq*dqdxyz(jx,iloc)
              yregdrv(nregionj) = yregdrv(nregionj) + dedq*dqdxyz(jy,iloc)
              zregdrv(nregionj) = zregdrv(nregionj) + dedq*dqdxyz(jz,iloc)
            endif
          enddo
          if (lgrad2.and.lDoQDeriv2) then
            do mm = 1,nqatoms(iloc)
              j = nqatomptr(mm,iloc)
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
!  Second internal derivatives : dE/dQ x d2Q/da.db
!
              derv2(jx,ix) = derv2(jx,ix) + dedq*d2qdxyz2(mmxx,iloc)
              derv2(jy,ix) = derv2(jy,ix) + dedq*d2qdxyz2(mmxy,iloc)
              derv2(jz,ix) = derv2(jz,ix) + dedq*d2qdxyz2(mmxz,iloc)
              derv2(jx,iy) = derv2(jx,iy) + dedq*d2qdxyz2(mmxy,iloc)
              derv2(jy,iy) = derv2(jy,iy) + dedq*d2qdxyz2(mmyy,iloc)
              derv2(jz,iy) = derv2(jz,iy) + dedq*d2qdxyz2(mmyz,iloc)
              derv2(jx,iz) = derv2(jx,iz) + dedq*d2qdxyz2(mmxz,iloc)
              derv2(jy,iz) = derv2(jy,iz) + dedq*d2qdxyz2(mmyz,iloc)
              derv2(jz,iz) = derv2(jz,iz) + dedq*d2qdxyz2(mmzz,iloc)
!
!  Second internal derivatives : d2E/dQ2 x dQ/da x dQ/db : i - j
!
              derv2(jx,ix) = derv2(jx,ix) + d2edq2*dqdxyz(ix,iloc)*dqdxyz(jx,iloc)
              derv2(jy,ix) = derv2(jy,ix) + d2edq2*dqdxyz(ix,iloc)*dqdxyz(jy,iloc)
              derv2(jz,ix) = derv2(jz,ix) + d2edq2*dqdxyz(ix,iloc)*dqdxyz(jz,iloc)
              derv2(jx,iy) = derv2(jx,iy) + d2edq2*dqdxyz(iy,iloc)*dqdxyz(jx,iloc)
              derv2(jy,iy) = derv2(jy,iy) + d2edq2*dqdxyz(iy,iloc)*dqdxyz(jy,iloc)
              derv2(jz,iy) = derv2(jz,iy) + d2edq2*dqdxyz(iy,iloc)*dqdxyz(jz,iloc)
              derv2(jx,iz) = derv2(jx,iz) + d2edq2*dqdxyz(iz,iloc)*dqdxyz(jx,iloc)
              derv2(jy,iz) = derv2(jy,iz) + d2edq2*dqdxyz(iz,iloc)*dqdxyz(jy,iloc)
              derv2(jz,iz) = derv2(jz,iz) + d2edq2*dqdxyz(iz,iloc)*dqdxyz(jz,iloc)
!
              do kk = 1,nstrains
!                   
!  Strain-internal second derivatives : dE/dQ x d2Q/ds.da
!               
                derv3(ix,kk) = derv3(ix,kk) - dedq*d2qdxyzs(kk,mx,iloc)
                derv3(iy,kk) = derv3(iy,kk) - dedq*d2qdxyzs(kk,my,iloc)
                derv3(iz,kk) = derv3(iz,kk) - dedq*d2qdxyzs(kk,mz,iloc)
!               
!  Strain-internal second derivatives : d2E/dQ2 x dQ/ds x dQ/da
!  
                derv3(ix,kk) = derv3(ix,kk) - d2edq2*dqdxyz(jx,iloc)*dqds(kk,iloc)
                derv3(iy,kk) = derv3(iy,kk) - d2edq2*dqdxyz(jy,iloc)*dqds(kk,iloc)
                derv3(iz,kk) = derv3(iz,kk) - d2edq2*dqdxyz(jz,iloc)*dqds(kk,iloc)
              enddo
            enddo
          endif
        endif
      endif
      if (lreg2pair) then
        esregion2 = esregion2 + esum
      else
        eboQself = eboQself + esum
      endif
    enddo
  enddo
  if (lgrad2.and.lDoQDeriv2) then
!*******************************************************
!  Calculate self-energy of charges                    *
!  Loop over all atoms for extra 2nd derivative terms  *
!*******************************************************
    call outerror('parallel second derivatives for variable Q not done',0_i4) 
    call stopnow('boselfd') 
  endif
!
  t2 = g_cpu_time()
  tbondorder = tbondorder + t2 - t1
#ifdef TRACE
  call trace_out('boselfd')
#endif
!
  return
  end
