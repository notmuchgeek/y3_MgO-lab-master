  subroutine real1Dpd(xkv)
!
!  This subroutine calculates the extra contributions to the
!  second derivatives of the electrostatic energy needed for
!  a phonon calculation.
!  Distributed memory version.
!
!   7/12 Created from real1Dp
!   9/16 cputime renamed to g_cpu_time
!   9/16 constants renamed to g_constants
!  12/16 derv2dk added back
!   8/17 On-diagonal blocks corrected by excluding double counting
!   8/17 Energies now handled as per serial to ensure matching
!        convergence
!   2/18 Trace added
!   3/20 Location of angstoev changed to current
!   3/20 Centre of mass positions passed to hfunc/emfunc 
!   5/20 Rigid molecule modifications added
!   5/20 Phasing based on centre of mass added
!   5/20 Phasing of group velocities modified for rigid molecules
!   6/20 Corrections to self terms to handle +acell and -acell
!   7/20 Separate routine for sumall with 1 argument added
!  11/20 Handling of ereal = 0 added to convergence check & check for non-zero charge
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
!  No claim is made that this program is free from errors and
!  no liability will be accepted for any loss or damage that
!  may result. The user is responsible for checking the validity
!  of their results.
!
!  Copyright Curtin University 2020
!
!  Julian Gale, CIC, Curtin University, November 2020
!
  use g_constants,    only : too_small
  use control,        only : lDoQDeriv2, lgroupvelocity, lrigid, keyword
  use current
  use derivatives
  use element,        only : maxele
  use general,        only : nmaxcells, nemorder, smallself
  use iochannels,     only : ioout
  use kspace,         only : accf1D
  use molecule
  use parallel,       only : natomsonnode, node2atom, nprocs, ioproc
  use qmedata,        only : maxloop
  use shells,         only : cuts
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed arguments
!
  real(dp),                    intent(in)    :: xkv
!
!  Local variables
!
  integer(i4)                                :: i,j 
  integer(i4)                                :: iloc
  integer(i4)                                :: ix,iy,iz 
  integer(i4)                                :: ixf,iyf,izf
  integer(i4)                                :: jx,jy,jz 
  integer(i4)                                :: m
  integer(i4)                                :: nati
  integer(i4)                                :: natj
  integer(i4)                                :: nmi
  integer(i4)                                :: nmj
  logical                                    :: lconverged
  logical                                    :: lcspair
  logical                                    :: lgamma
  logical                                    :: lqnz
  complex(dpc)                               :: cdk
  real(dp)                                   :: accf
  real(dp)                                   :: acell
  real(dp)                                   :: cosk
  real(dp)                                   :: cos2k
  real(dp)                                   :: g_cpu_time
  real(dp)                                   :: cut2s
  real(dp)                                   :: d0
  real(dp)                                   :: d0i
  real(dp)                                   :: d0j
  real(dp)                                   :: d0term
  real(dp)                                   :: d1
  real(dp)                                   :: d1i
  real(dp)                                   :: d1ix
  real(dp)                                   :: d1iy
  real(dp)                                   :: d1iz
  real(dp)                                   :: d1j
  real(dp)                                   :: d1jx
  real(dp)                                   :: d1jy
  real(dp)                                   :: d1jz
  real(dp)                                   :: d2
  real(dp)                                   :: d2i
  real(dp)                                   :: d2i2
  real(dp)                                   :: d2ij
  real(dp)                                   :: d2j2
  real(dp)                                   :: dh1(3)
  real(dp)                                   :: dh2(3)
  real(dp)                                   :: d2h1(6)
  real(dp)                                   :: d2h2(6)
  real(dp)                                   :: d3h1(10)
  real(dp)                                   :: ediff
  real(dp)                                   :: elast
  real(dp)                                   :: ereal
  real(dp)                                   :: esum
  real(dp)                                   :: esumem
  real(dp)                                   :: esumh
  real(dp)                                   :: e1
  real(dp)                                   :: e2
  real(dp)                                   :: h1
  real(dp)                                   :: h2
  real(dp)                                   :: lna
  real(dp)                                   :: oci     
  real(dp)                                   :: ocj 
  real(dp)                                   :: qi  
  real(dp)                                   :: qj
  real(dp)                                   :: qii
  real(dp)                                   :: qij
  real(dp)                                   :: r
  real(dp)                                   :: rcut
  real(dp)                                   :: rr
  real(dp)                                   :: sink
  real(dp)                                   :: t1, t2
  real(dp)                                   :: u
  real(dp)                                   :: x
  real(dp)                                   :: y
  real(dp)                                   :: z
  real(dp)                                   :: xcom
  real(dp)                                   :: xcomi
#ifdef TRACE
  call trace_in('real1Dpd')
#endif
!
!  Check for non-zero charge
!
  lqnz = .false.
  do i = 1,numat
    lqnz = (abs(qf(i)).gt.1.0d-10)
    if (lqnz) exit
  enddo
  if (.not.lqnz) return
!
  t1 = g_cpu_time()
!********************************************************
!  Calculate Coulomb sum converged to desired accuracy  *
!********************************************************
!
!  Loop over number of cells in sum
!
  accf = 10.0**(-accf1D)
  lna = log(a)
  cut2s = cuts*cuts
  m = - 1
  lconverged = .false.
  lgamma = (abs(xkv).lt.too_small)
!
  if (ioproc.and.index(keyword,'verb').ne.0) then
    write(ioout,'(/,'' Convergence of 1-D electrostatic summation:'',/)')
  endif
!
  elast = 0.0_dp
  ereal = 0.0_dp
  esum = 0.0_dp
  do while (m.lt.nmaxcells.and..not.lconverged) 
    m = m + 1
!
!  Direct sum component over neutral cells
!
    acell = dble(m)*a
    ix = -2
    iy = -1
    iz =  0
    do iloc = 1,natomsonnode
      i = node2atom(iloc)
      nati = nat(i)
      oci = occuf(i)
      qi = qf(i)*oci
!
!  Set COM coordinates
!
      nmi = natmol(i)
      if (lrigid.and.nmi.gt.0) then
        xcomi = molxyz(1,natinmol(i),nmi)
      else
        xcomi = 0.0_dp
      endif
!
      ix = ix + 3
      iy = iy + 3
      iz = iz + 3
      jx = -2
      jy = -1
      jz =  0
      do j = 1,numat
        jx = jx + 3
        jy = jy + 3
        jz = jz + 3
!
!  Exclude self term which is handled separately
!
        if (j.ne.i) then
          natj = nat(j)
          ocj = occuf(j)
          qj = qf(j)*ocj
          lcspair = (abs(nati-natj).eq.maxele.or.(oci+ocj).lt.1.0001d0)
          if (lcspair) then
            rcut = cut2s
          else
            rcut = smallself
          endif
!
!  Set COM coordinates
!
          nmj = natmol(j)
          if (lrigid.and.nmj.gt.0) then
            xcom = molxyz(1,natinmol(j),nmj) - xcomi
          else
            xcom = - xcomi
          endif
!
          x = acell + xclat(j) - xclat(i)
          y = yclat(j) - yclat(i)
          z = zclat(j) - zclat(i)
          r = x*x + y*y + z*z
          if (r.gt.rcut) then
            if (lphasecom) then
              cosk = xkv*(x - xcom)
            else
              cosk = xkv*x
            endif
            sink = sin(cosk)
            cosk = cos(cosk)
            r = sqrt(r)
            rr = 1.0_dp/r
            d0 = qi*qj*rr
            esum = esum + 0.5_dp*d0
            d1 = d0*rr*rr*angstoev
            d2 = - 3.0_dp*d1*rr*rr
            d1i = d1*sink
            d2i = d2*sink
            d1 = d1*cosk
            d2 = d2*cosk
!
!  Real terms
!
            derv2(jx,ix) = derv2(jx,ix) + d2*x*x
            derv2(jy,ix) = derv2(jy,ix) + d2*y*x
            derv2(jz,ix) = derv2(jz,ix) + d2*z*x
            derv2(jx,iy) = derv2(jx,iy) + d2*x*y
            derv2(jy,iy) = derv2(jy,iy) + d2*y*y
            derv2(jz,iy) = derv2(jz,iy) + d2*z*y
            derv2(jx,iz) = derv2(jx,iz) + d2*x*z
            derv2(jy,iz) = derv2(jy,iz) + d2*y*z
            derv2(jz,iz) = derv2(jz,iz) + d2*z*z
            derv2(jx,ix) = derv2(jx,ix) + d1
            derv2(jy,iy) = derv2(jy,iy) + d1
            derv2(jz,iz) = derv2(jz,iz) + d1
!
!  Imaginary terms
!
            dervi(jx,ix) = dervi(jx,ix) + d2i*x*x
            dervi(jy,ix) = dervi(jy,ix) + d2i*y*x
            dervi(jz,ix) = dervi(jz,ix) + d2i*z*x
            dervi(jx,iy) = dervi(jx,iy) + d2i*x*y
            dervi(jy,iy) = dervi(jy,iy) + d2i*y*y
            dervi(jz,iy) = dervi(jz,iy) + d2i*z*y
            dervi(jx,iz) = dervi(jx,iz) + d2i*x*z
            dervi(jy,iz) = dervi(jy,iz) + d2i*y*z
            dervi(jz,iz) = dervi(jz,iz) + d2i*z*z
            dervi(jx,ix) = dervi(jx,ix) + d1i
            dervi(jy,iy) = dervi(jy,iy) + d1i
            dervi(jz,iz) = dervi(jz,iz) + d1i
!
!  Group velocities
!
            if (lgroupvelocity) then
              if (lphasecom) then
                cdk = dcmplx(d2*(x-xcom),d2i*(x-xcom))*dcmplx(0.0_dp,1.0_dp)
              else
                cdk = dcmplx(d2*x,d2i*x)*dcmplx(0.0_dp,1.0_dp)
              endif
!
              derv2dk(1,jx,ix) = derv2dk(1,jx,ix) + cdk*x*x
              derv2dk(1,jy,ix) = derv2dk(1,jy,ix) + cdk*y*x
              derv2dk(1,jz,ix) = derv2dk(1,jz,ix) + cdk*z*x
              derv2dk(1,jx,iy) = derv2dk(1,jx,iy) + cdk*x*y
              derv2dk(1,jy,iy) = derv2dk(1,jy,iy) + cdk*y*y
              derv2dk(1,jz,iy) = derv2dk(1,jz,iy) + cdk*z*y
              derv2dk(1,jx,iz) = derv2dk(1,jx,iz) + cdk*x*z
              derv2dk(1,jy,iz) = derv2dk(1,jy,iz) + cdk*y*z
              derv2dk(1,jz,iz) = derv2dk(1,jz,iz) + cdk*z*z
!
              if (lphasecom) then
                cdk = dcmplx(d1*(x-xcom),d1i*(x-xcom))*dcmplx(0.0_dp,1.0_dp)
              else
                cdk = dcmplx(d1*x,d1i*x)*dcmplx(0.0_dp,1.0_dp)
              endif
!
              derv2dk(1,jx,ix) = derv2dk(1,jx,ix) + cdk
              derv2dk(1,jy,iy) = derv2dk(1,jy,iy) + cdk
              derv2dk(1,jz,iz) = derv2dk(1,jz,iz) + cdk
            endif
!**********************************
!  Variable charge contributions  *
!**********************************
            if (lDoQDeriv2) then
              d0i  = qj*rr*angstoev
              d0j  = qi*rr*angstoev
              d1i  = - d0i*rr*rr
              d1j  = - d0j*rr*rr
              d2i2 = 0.0_dp
              d2ij = rr*angstoev
              d2j2 = 0.0_dp
              d1ix = d1i*x
              d1iy = d1i*y
              d1iz = d1i*z
              d1jx = d1j*x
              d1jy = d1j*y
              d1jz = d1j*z
              call d2chargep(i,j,1_i4,ix,iy,iz,jx,jy,jz,xkv,0.0_dp,0.0_dp,d0i,d0j,d1ix,d1iy,d1iz, &
                             d1jx,d1jy,d1jz,d2i2,d2ij,d2j2,0.0_dp,0.0_dp,0.0_dp,.true.)
            endif
          endif
          if (m.gt.0) then
            x = - acell + xclat(j) - xclat(i)
            r = x*x + y*y + z*z
            if (r.gt.rcut) then
              if (lphasecom) then
                cosk = xkv*(x - xcom)
              else
                cosk = xkv*x
              endif
              sink = sin(cosk)
              cosk = cos(cosk)
              r = sqrt(r)
              rr = 1.0_dp/r
              d0 = qi*qj*rr
              esum = esum + 0.5_dp*d0
              d1 = d0*rr*rr*angstoev
              d2 = - 3.0_dp*d1*rr*rr
              d1i = d1*sink
              d2i = d2*sink
              d1 = d1*cosk
              d2 = d2*cosk
!
!  Real terms
!
              derv2(jx,ix) = derv2(jx,ix) + d2*x*x
              derv2(jy,ix) = derv2(jy,ix) + d2*y*x
              derv2(jz,ix) = derv2(jz,ix) + d2*z*x
              derv2(jx,iy) = derv2(jx,iy) + d2*x*y
              derv2(jy,iy) = derv2(jy,iy) + d2*y*y
              derv2(jz,iy) = derv2(jz,iy) + d2*z*y
              derv2(jx,iz) = derv2(jx,iz) + d2*x*z
              derv2(jy,iz) = derv2(jy,iz) + d2*y*z
              derv2(jz,iz) = derv2(jz,iz) + d2*z*z
              derv2(jx,ix) = derv2(jx,ix) + d1
              derv2(jy,iy) = derv2(jy,iy) + d1
              derv2(jz,iz) = derv2(jz,iz) + d1
!
!  Imaginary terms
!
              dervi(jx,ix) = dervi(jx,ix) + d2i*x*x
              dervi(jy,ix) = dervi(jy,ix) + d2i*y*x
              dervi(jz,ix) = dervi(jz,ix) + d2i*z*x
              dervi(jx,iy) = dervi(jx,iy) + d2i*x*y
              dervi(jy,iy) = dervi(jy,iy) + d2i*y*y
              dervi(jz,iy) = dervi(jz,iy) + d2i*z*y
              dervi(jx,iz) = dervi(jx,iz) + d2i*x*z
              dervi(jy,iz) = dervi(jy,iz) + d2i*y*z
              dervi(jz,iz) = dervi(jz,iz) + d2i*z*z
              dervi(jx,ix) = dervi(jx,ix) + d1i
              dervi(jy,iy) = dervi(jy,iy) + d1i
              dervi(jz,iz) = dervi(jz,iz) + d1i
!
!  Group velocities
!
              if (lgroupvelocity) then
                if (lphasecom) then
                  cdk = dcmplx(d2*(x-xcom),d2i*(x-xcom))*dcmplx(0.0_dp,1.0_dp)
                else
                  cdk = dcmplx(d2*x,d2i*x)*dcmplx(0.0_dp,1.0_dp)
                endif
!
                derv2dk(1,jx,ix) = derv2dk(1,jx,ix) + cdk*x*x
                derv2dk(1,jy,ix) = derv2dk(1,jy,ix) + cdk*y*x
                derv2dk(1,jz,ix) = derv2dk(1,jz,ix) + cdk*z*x
                derv2dk(1,jx,iy) = derv2dk(1,jx,iy) + cdk*x*y
                derv2dk(1,jy,iy) = derv2dk(1,jy,iy) + cdk*y*y
                derv2dk(1,jz,iy) = derv2dk(1,jz,iy) + cdk*z*y
                derv2dk(1,jx,iz) = derv2dk(1,jx,iz) + cdk*x*z
                derv2dk(1,jy,iz) = derv2dk(1,jy,iz) + cdk*y*z
                derv2dk(1,jz,iz) = derv2dk(1,jz,iz) + cdk*z*z
!
                if (lphasecom) then
                  cdk = dcmplx(d1*(x-xcom),d1i*(x-xcom))*dcmplx(0.0_dp,1.0_dp)
                else
                  cdk = dcmplx(d1*x,d1i*x)*dcmplx(0.0_dp,1.0_dp)
                endif
!
                derv2dk(1,jx,ix) = derv2dk(1,jx,ix) + cdk
                derv2dk(1,jy,iy) = derv2dk(1,jy,iy) + cdk
                derv2dk(1,jz,iz) = derv2dk(1,jz,iz) + cdk
              endif
!**********************************
!  Variable charge contributions  *
!**********************************
              if (lDoQDeriv2) then
                d0i  = qj*rr*angstoev
                d0j  = qi*rr*angstoev
                d1i  = - d0i*rr*rr
                d1j  = - d0j*rr*rr
                d2i2 = 0.0_dp
                d2ij = rr*angstoev
                d2j2 = 0.0_dp
                d1ix = d1i*x 
                d1iy = d1i*y
                d1iz = d1i*z 
                d1jx = d1j*x 
                d1jy = d1j*y 
                d1jz = d1j*z 
                call d2chargep(i,j,1_i4,ix,iy,iz,jx,jy,jz,xkv,0.0_dp,0.0_dp,d0i,d0j,d1ix,d1iy,d1iz, &
                               d1jx,d1jy,d1jz,d2i2,d2ij,d2j2,0.0_dp,0.0_dp,0.0_dp,.true.)                   
              endif
            endif
          endif
        endif
      enddo
    enddo
!
!  Self interactions
!
    ix = -2
    iy = -1
    iz =  0
    do iloc = 1,natomsonnode
      i = node2atom(iloc)
      oci = occuf(i)
      qi = qf(i)*oci
      ix = ix + 3
      iy = iy + 3
      iz = iz + 3
      ixf = 3*(i-1) + 1
      iyf = ixf + 1
      izf = iyf + 1
      if (acell.gt.smallself) then
        rr = 1.0_dp/acell
        d0 = qi*qi*rr
        esum = esum + d0
! 
!  Phase factors
!
        cosk = xkv*acell
        sink = sin(cosk)
        cosk = cos(cosk) - 1.0_dp
        cos2k = 2.0_dp*cosk
!
        d1 = d0*rr*rr*angstoev
        d2 = - 3.0_dp*d1*rr*rr
!
!  Real terms
!
        derv2(ixf,ix) = derv2(ixf,ix) + d2*cos2k*acell*acell
        derv2(ixf,ix) = derv2(ixf,ix) + d1*cos2k
        derv2(iyf,iy) = derv2(iyf,iy) + d1*cos2k
        derv2(izf,iz) = derv2(izf,iz) + d1*cos2k
!
!  Group velocities
!
        if (lgroupvelocity) then
!
!  +acell and -acell
!
          cdk = dcmplx(d2*cosk*acell,d2*sink*acell)*dcmplx(0.0_dp,1.0_dp) + &
                dcmplx(-d2*cosk*acell,d2*sink*acell)*dcmplx(0.0_dp,1.0_dp)
!
          derv2dk(1,ixf,ix) = derv2dk(1,ixf,ix) + cdk*acell*acell
!
!
!  +acell and -acell
!
          cdk = dcmplx(d1*cosk*acell,d1*sink*acell)*dcmplx(0.0_dp,1.0_dp) + &
                dcmplx(-d1*cosk*acell,d1*sink*acell)*dcmplx(0.0_dp,1.0_dp)
!
          derv2dk(1,ixf,ix) = derv2dk(1,ixf,ix) + cdk
          derv2dk(1,iyf,iy) = derv2dk(1,iyf,iy) + cdk
          derv2dk(1,izf,iz) = derv2dk(1,izf,iz) + cdk
        endif
!*******************************************************
!  Variable charge contribution to second derivatives  *
!*******************************************************
        if (lDoQDeriv2) then
          d0i  = qi*rr*angstoev
          d1i  = - d0i*rr*rr
          d2i2 = 0.0_dp
          d2ij = rr*angstoev
          d2j2 = 0.0_dp
          d1ix = d1i*acell
          call d2chargep(i,i,1_i4,ix,iy,iz,ix,iy,iz,xkv,0.0_dp,0.0_dp,d0i,d0i,d1ix,0.0_dp,0.0_dp, &
                         d1ix,0.0_dp,0.0_dp,d2i2,d2ij,d2j2,0.0_dp,0.0_dp,0.0_dp,.true.)                   
        endif
      endif
    enddo
!
!  Neutralising terms
!
!  Background
!
!  and
!
!  Euler-MacLaurin component
!
    esumh = 0.0_dp
    esumem = 0.0_dp
    if (lgamma.and.m.gt.0) then
      u = (dble(m)+0.5_dp)*a
      do iloc = 1,natomsonnode
        i = node2atom(iloc)
        oci = occuf(i)
        qi = qf(i)*oci
        do j = 1,numat
          ocj = occuf(j)
          qj = qf(j)*ocj
          qij = qi*qj
          x = xclat(j) - xclat(i)
          y = yclat(j) - yclat(i)
          z = zclat(j) - zclat(i)
          call hfunc(u,+x,1.0_dp,y,z,h1,dh1,d2h1,d3h1,.false.,.false.,.false.)
          call hfunc(u,-x,-1.0_dp,y,z,h2,dh2,d2h2,d3h1,.false.,.false.,.false.)
          esumh = esumh - 0.5_dp*qij*(h1 + h2 - 2.0_dp*lna)/a
          call emfunc(nemorder,u,+x,1.0_dp,y,z,a,e1,dh1,d2h1,d3h1,.false.,.false.,.false.)
          call emfunc(nemorder,u,-x,-1.0_dp,y,z,a,e2,dh2,d2h2,d3h1,.false.,.false.,.false.)
          esumem = esumem + 0.5_dp*qij*(e1 + e2)
        enddo
      enddo
    endif
!
!  Sum up terms
!
    if (nprocs.gt.1) then
      call sumone(esum+esumh+esumem,ereal,"real1Dpd","ereal")
    else
      ereal = esum + esumh + esumem
    endif
!
!  Compare to energy with previous number of cells
!  and check for convergence to the required
!  accuracy.
!
    if (abs(ereal).lt.accf.and.m.gt.0) lconverged = .true.
    if (.not.lconverged) then
      if (abs(ereal).gt.1.0d-8) then
        ediff = abs((ereal - elast)/ereal)
      else
        ediff = abs(ereal - elast)
      endif
      lconverged = (ediff.lt.accf)
    endif
    elast = ereal
    if (ioproc.and.index(keyword,'verb').ne.0) then
      write(ioout,'('' Mcell = '',i4,''  Ereal(1-D) = '',f15.6,'' eV'')') m,ereal
    endif
  enddo
  ereal = ereal*angstoev
!
!  Save number of cells needed
!
  maxloop(1) = m
  if (lgamma.and.m.gt.0) then
!
!  Derivatives of Euler-MacLaurin terms since these are not cumulative
!
    u = (dble(m)+0.5_dp)*a
    ix = -2
    iy = -1
    iz =  0
    do iloc = 1,natomsonnode
      i = node2atom(iloc)
      oci = occuf(i)  
      qi = qf(i)*oci
!
!  Set COM coordinates
!
      nmi = natmol(i)
      if (lrigid.and.nmi.gt.0) then
        xcomi = molxyz(1,natinmol(i),nmi)
      else
        xcomi = 0.0_dp
      endif
      ix = ix + 3
      iy = iy + 3
      iz = iz + 3
      jx = -2
      jy = -1
      jz =  0
      do j = 1,numat
        jx = jx + 3
        jy = jy + 3
        jz = jz + 3  
!
!  Exclude self interactions since these are handled separately
!
        if (j.ne.i) then
          ocj = occuf(j)
          qj = qf(j)*ocj
          qij = qi*qj
!
!  Set COM coordinates
!
          nmj = natmol(j)
          if (lrigid.and.nmj.gt.0) then
            xcom = molxyz(1,natinmol(j),nmj) - xcomi
          else
            xcom = - xcomi
          endif
          x = xclat(j) - xclat(i)
          y = yclat(j) - yclat(i)
          z = zclat(j) - zclat(i)
!
!  Phase factor
!
          if (lphasecom) then
            cosk = xkv*(x - xcom)
          else
            cosk = xkv*x
          endif
          sink = sin(cosk)
          cosk = cos(cosk)
!
!  H term
!
          call hfunc(u,+x,1.0_dp,y,z,h1,dh1,d2h1,d3h1,.true.,.true.,.false.)
          call hfunc(u,-x,-1.0_dp,y,z,h2,dh2,d2h2,d3h1,.true.,.true.,.false.)
          d2 = qij*angstoev/a
          d2i = d2*sink
          d2 = d2*cosk
!
!  Real terms
!
          derv2(jx,ix) = derv2(jx,ix) + d2*(d2h1(1)+d2h2(1))
          derv2(jy,ix) = derv2(jy,ix) + d2*(d2h1(6)+d2h2(6))
          derv2(jz,ix) = derv2(jz,ix) + d2*(d2h1(5)+d2h2(5))
          derv2(jx,iy) = derv2(jx,iy) + d2*(d2h1(6)+d2h2(6))
          derv2(jy,iy) = derv2(jy,iy) + d2*(d2h1(2)+d2h2(2))
          derv2(jz,iy) = derv2(jz,iy) + d2*(d2h1(4)+d2h2(4))
          derv2(jx,iz) = derv2(jx,iz) + d2*(d2h1(5)+d2h2(5))
          derv2(jy,iz) = derv2(jy,iz) + d2*(d2h1(4)+d2h2(4))
          derv2(jz,iz) = derv2(jz,iz) + d2*(d2h1(3)+d2h2(3))
!
!  Imaginary terms
!
          dervi(jx,ix) = dervi(jx,ix) + d2i*(d2h1(1)+d2h2(1))
          dervi(jy,ix) = dervi(jy,ix) + d2i*(d2h1(6)+d2h2(6))
          dervi(jz,ix) = dervi(jz,ix) + d2i*(d2h1(5)+d2h2(5))
          dervi(jx,iy) = dervi(jx,iy) + d2i*(d2h1(6)+d2h2(6))
          dervi(jy,iy) = dervi(jy,iy) + d2i*(d2h1(2)+d2h2(2))
          dervi(jz,iy) = dervi(jz,iy) + d2i*(d2h1(4)+d2h2(4))
          dervi(jx,iz) = dervi(jx,iz) + d2i*(d2h1(5)+d2h2(5))
          dervi(jy,iz) = dervi(jy,iz) + d2i*(d2h1(4)+d2h2(4))
          dervi(jz,iz) = dervi(jz,iz) + d2i*(d2h1(3)+d2h2(3))
!
!  Group velocities
!
          if (lgroupvelocity) then
            if (lphasecom) then
              cdk = dcmplx(d2*(x-xcom),d2i*(x-xcom))*dcmplx(0.0_dp,1.0_dp)
            else
              cdk = dcmplx(d2*x,d2i*x)*dcmplx(0.0_dp,1.0_dp)
            endif
!
            derv2dk(1,jx,ix) = derv2dk(1,jx,ix) + cdk*(d2h1(1)+d2h2(1))
            derv2dk(1,jy,ix) = derv2dk(1,jy,ix) + cdk*(d2h1(6)+d2h2(6))
            derv2dk(1,jz,ix) = derv2dk(1,jz,ix) + cdk*(d2h1(5)+d2h2(5))
            derv2dk(1,jx,iy) = derv2dk(1,jx,iy) + cdk*(d2h1(6)+d2h2(6))
            derv2dk(1,jy,iy) = derv2dk(1,jy,iy) + cdk*(d2h1(2)+d2h2(2))
            derv2dk(1,jz,iy) = derv2dk(1,jz,iy) + cdk*(d2h1(4)+d2h2(4))
            derv2dk(1,jx,iz) = derv2dk(1,jx,iz) + cdk*(d2h1(5)+d2h2(5))
            derv2dk(1,jy,iz) = derv2dk(1,jy,iz) + cdk*(d2h1(4)+d2h2(4))
            derv2dk(1,jz,iz) = derv2dk(1,jz,iz) + cdk*(d2h1(3)+d2h2(3))
          endif
!
          if (lDoQDeriv2) then
!
!  Accumulate terms that contribute to the variable charge second derivatives
!                 
            d0term = - angstoev*(h1 + h2 - 2.0_dp*lna)/a
            d1ix = - qj*angstoev*(dh1(1) + dh2(1))/a
            d1iy = - qj*angstoev*(dh1(2) + dh2(2))/a
            d1iz = - qj*angstoev*(dh1(3) + dh2(3))/a
            d1jx = - qi*angstoev*(dh1(1) + dh2(1))/a
            d1jy = - qi*angstoev*(dh1(2) + dh2(2))/a
            d1jz = - qi*angstoev*(dh1(3) + dh2(3))/a
          endif
!
!  E-M term
!
          call emfunc(nemorder,u,+x,1.0_dp,y,z,a,e1,dh1,d2h1,d3h1,.true.,.true.,.false.)
          call emfunc(nemorder,u,-x,-1.0_dp,y,z,a,e2,dh2,d2h2,d3h1,.true.,.true.,.false.)
          d2 = qij*angstoev
          d2i = d2*sink
          d2 = d2*cosk
!
!  Real terms
!
          derv2(jx,ix) = derv2(jx,ix) - d2*(d2h1(1)+d2h2(1))
          derv2(jy,ix) = derv2(jy,ix) - d2*(d2h1(6)+d2h2(6))
          derv2(jz,ix) = derv2(jz,ix) - d2*(d2h1(5)+d2h2(5))
          derv2(jx,iy) = derv2(jx,iy) - d2*(d2h1(6)+d2h2(6))
          derv2(jy,iy) = derv2(jy,iy) - d2*(d2h1(2)+d2h2(2))
          derv2(jz,iy) = derv2(jz,iy) - d2*(d2h1(4)+d2h2(4))
          derv2(jx,iz) = derv2(jx,iz) - d2*(d2h1(5)+d2h2(5))
          derv2(jy,iz) = derv2(jy,iz) - d2*(d2h1(4)+d2h2(4))
          derv2(jz,iz) = derv2(jz,iz) - d2*(d2h1(3)+d2h2(3))
!
!  Imaginary terms
!
          dervi(jx,ix) = dervi(jx,ix) - d2i*(d2h1(1)+d2h2(1))
          dervi(jy,ix) = dervi(jy,ix) - d2i*(d2h1(6)+d2h2(6))
          dervi(jz,ix) = dervi(jz,ix) - d2i*(d2h1(5)+d2h2(5))
          dervi(jx,iy) = dervi(jx,iy) - d2i*(d2h1(6)+d2h2(6))
          dervi(jy,iy) = dervi(jy,iy) - d2i*(d2h1(2)+d2h2(2))
          dervi(jz,iy) = dervi(jz,iy) - d2i*(d2h1(4)+d2h2(4))
          dervi(jx,iz) = dervi(jx,iz) - d2i*(d2h1(5)+d2h2(5))
          dervi(jy,iz) = dervi(jy,iz) - d2i*(d2h1(4)+d2h2(4))
          dervi(jz,iz) = dervi(jz,iz) - d2i*(d2h1(3)+d2h2(3))
!
!  Group velocities
!
          if (lgroupvelocity) then
            if (lphasecom) then
              cdk = dcmplx(d2*(x-xcom),d2i*(x-xcom))*dcmplx(0.0_dp,1.0_dp)
            else
              cdk = dcmplx(d2*x,d2i*x)*dcmplx(0.0_dp,1.0_dp)
            endif
!
            derv2dk(1,jx,ix) = derv2dk(1,jx,ix) - cdk*(d2h1(1)+d2h2(1))
            derv2dk(1,jy,ix) = derv2dk(1,jy,ix) - cdk*(d2h1(6)+d2h2(6))
            derv2dk(1,jz,ix) = derv2dk(1,jz,ix) - cdk*(d2h1(5)+d2h2(5))
            derv2dk(1,jx,iy) = derv2dk(1,jx,iy) - cdk*(d2h1(6)+d2h2(6))
            derv2dk(1,jy,iy) = derv2dk(1,jy,iy) - cdk*(d2h1(2)+d2h2(2))
            derv2dk(1,jz,iy) = derv2dk(1,jz,iy) - cdk*(d2h1(4)+d2h2(4))
            derv2dk(1,jx,iz) = derv2dk(1,jx,iz) - cdk*(d2h1(5)+d2h2(5))
            derv2dk(1,jy,iz) = derv2dk(1,jy,iz) - cdk*(d2h1(4)+d2h2(4))
            derv2dk(1,jz,iz) = derv2dk(1,jz,iz) - cdk*(d2h1(3)+d2h2(3))
          endif
!*******************************************************
!  Variable charge contribution to second derivatives  *
!*******************************************************
          if (lDoQDeriv2) then
            d0term = d0term + (e1 + e2)*angstoev
            d0i  = qj*d0term
            d0j  = qi*d0term
            d2i2 = 0.0_dp
            d2ij = d0term
            d2j2 = 0.0_dp
            d1ix = d1ix + qj*angstoev*(dh1(1) + dh2(1))
            d1iy = d1iy + qj*angstoev*(dh1(2) + dh2(2))
            d1iz = d1iz + qj*angstoev*(dh1(3) + dh2(3))
            d1jx = d1jx + qi*angstoev*(dh1(1) + dh2(1))
            d1jy = d1jy + qi*angstoev*(dh1(2) + dh2(2))
            d1jz = d1jz + qi*angstoev*(dh1(3) + dh2(3))
            call d2chargep(i,j,1_i4,ix,iy,iz,jx,jy,jz,xkv,0.0_dp,0.0_dp,d0i,d0j,d1ix,d1iy,d1iz, &
                           d1jx,d1jy,d1jz,d2i2,d2ij,d2j2,0.0_dp,0.0_dp,0.0_dp,.true.)
          endif
        endif
      enddo
    enddo
!
!  Self-interaction terms
!
    ix = -2
    iy = -1
    iz =  0
    do iloc = 1,natomsonnode
      i = node2atom(iloc)
      oci = occuf(i)
      qi = qf(i)*oci
      ix = ix + 3
      iy = iy + 3
      iz = iz + 3
      ixf = 3*(i-1) + 1
      iyf = ixf + 1
      izf = iyf + 1
      qii = qi*qi
!
!  Phase factors
! 
      cosk = xkv*acell
      sink = sin(cosk)
      cosk = cos(cosk) - 1.0_dp
      cos2k = 2.0_dp*cosk
!
!  H term
!
      call hfunc(u,0.0_dp,1.0_dp,0.0_dp,0.0_dp,h1,dh1,d2h1,d3h1,.true.,.true.,.false.)
      d2 = qii*angstoev/a
!
!  Real terms
!
      derv2(ixf,ix) = derv2(ixf,ix) + d2*d2h1(1)*cos2k
      derv2(iyf,ix) = derv2(iyf,ix) + d2*d2h1(6)*cos2k
      derv2(izf,ix) = derv2(izf,ix) + d2*d2h1(5)*cos2k
      derv2(ixf,iy) = derv2(ixf,iy) + d2*d2h1(6)*cos2k
      derv2(iyf,iy) = derv2(iyf,iy) + d2*d2h1(2)*cos2k
      derv2(izf,iy) = derv2(izf,iy) + d2*d2h1(4)*cos2k
      derv2(ixf,iz) = derv2(ixf,iz) + d2*d2h1(5)*cos2k
      derv2(iyf,iz) = derv2(iyf,iz) + d2*d2h1(4)*cos2k
      derv2(izf,iz) = derv2(izf,iz) + d2*d2h1(3)*cos2k
!
!  Group velocities
!
      if (lgroupvelocity) then
!
!  +acell and -acell
!
        cdk = dcmplx(d2*cosk*acell,d2*sink*acell)*dcmplx(0.0_dp,1.0_dp) + &
              dcmplx(-d2*cosk*acell,d2*sink*acell)*dcmplx(0.0_dp,1.0_dp)
!
        derv2dk(1,ixf,ix) = derv2dk(1,ixf,ix) + cdk*d2h1(1)
        derv2dk(1,iyf,ix) = derv2dk(1,iyf,ix) + cdk*d2h1(6)
        derv2dk(1,izf,ix) = derv2dk(1,izf,ix) + cdk*d2h1(5)
        derv2dk(1,ixf,iy) = derv2dk(1,ixf,iy) + cdk*d2h1(6)
        derv2dk(1,iyf,iy) = derv2dk(1,iyf,iy) + cdk*d2h1(2)
        derv2dk(1,izf,iy) = derv2dk(1,izf,iy) + cdk*d2h1(4)
        derv2dk(1,ixf,iz) = derv2dk(1,ixf,iz) + cdk*d2h1(5)
        derv2dk(1,iyf,iz) = derv2dk(1,iyf,iz) + cdk*d2h1(4)
        derv2dk(1,izf,iz) = derv2dk(1,izf,iz) + cdk*d2h1(3)
      endif
!
      if (lDoQDeriv2) then
!
!  Accumulate terms that contribute to the variable charge second derivatives
!             
        d0term = - angstoev*(h1 - lna)/a
        d1ix = - qi*angstoev*dh1(1)/a
        d1iy = - qi*angstoev*dh1(2)/a
        d1iz = - qi*angstoev*dh1(3)/a
      endif
!
      call emfunc(nemorder,u,0.0_dp,1.0_dp,0.0_dp,0.0_dp,a,e1,dh1,d2h1,d3h1,.true.,.true.,.false.)
      d2 = qii*angstoev
!
!  Real terms
!
      derv2(ixf,ix) = derv2(ixf,ix) - d2*d2h1(1)*cos2k
      derv2(iyf,ix) = derv2(iyf,ix) - d2*d2h1(6)*cos2k
      derv2(izf,ix) = derv2(izf,ix) - d2*d2h1(5)*cos2k
      derv2(ixf,iy) = derv2(ixf,iy) - d2*d2h1(6)*cos2k
      derv2(iyf,iy) = derv2(iyf,iy) - d2*d2h1(2)*cos2k
      derv2(izf,iy) = derv2(izf,iy) - d2*d2h1(4)*cos2k
      derv2(ixf,iz) = derv2(ixf,iz) - d2*d2h1(5)*cos2k
      derv2(iyf,iz) = derv2(iyf,iz) - d2*d2h1(4)*cos2k
      derv2(izf,iz) = derv2(izf,iz) - d2*d2h1(3)*cos2k
!
!  Group velocities
!
      if (lgroupvelocity) then
!
!  +acell and -acell
!
        cdk = dcmplx(d2*cosk*acell,d2*sink*acell)*dcmplx(0.0_dp,1.0_dp) + &
              dcmplx(-d2*cosk*acell,d2*sink*acell)*dcmplx(0.0_dp,1.0_dp)
!
        derv2dk(1,ixf,ix) = derv2dk(1,ixf,ix) - cdk*d2h1(1)
        derv2dk(1,iyf,ix) = derv2dk(1,iyf,ix) - cdk*d2h1(6)
        derv2dk(1,izf,ix) = derv2dk(1,izf,ix) - cdk*d2h1(5)
        derv2dk(1,ixf,iy) = derv2dk(1,ixf,iy) - cdk*d2h1(6)
        derv2dk(1,iyf,iy) = derv2dk(1,iyf,iy) - cdk*d2h1(2)
        derv2dk(1,izf,iy) = derv2dk(1,izf,iy) - cdk*d2h1(4)
        derv2dk(1,ixf,iz) = derv2dk(1,ixf,iz) - cdk*d2h1(5)
        derv2dk(1,iyf,iz) = derv2dk(1,iyf,iz) - cdk*d2h1(4)
        derv2dk(1,izf,iz) = derv2dk(1,izf,iz) - cdk*d2h1(3)
      endif
!*******************************************************
!  Variable charge contribution to second derivatives  *
!*******************************************************
      if (lDoQDeriv2) then
        d0term = d0term + e1*angstoev
        d0i  = qi*d0term
        d2i2 = 0.0_dp
        d2ij = d0term
        d2j2 = 0.0_dp
        d1ix = d1ix + qi*angstoev*dh1(1)
        d1iy = d1iy + qi*angstoev*dh1(2)
        d1iz = d1iz + qi*angstoev*dh1(3)
        call d2chargep(i,i,1_i4,ix,iy,iz,ix,iy,iz,xkv,0.0_dp,0.0_dp,d0i,d0i,d1ix,0.0_dp,0.0_dp, &
                       d1ix,0.0_dp,0.0_dp,d2i2,d2ij,d2j2,0.0_dp,0.0_dp,0.0_dp,.true.)                   
      endif
    enddo
  endif
!
  t2 = g_cpu_time()
  tatom = tatom + t2 - t1
#ifdef TRACE
  call trace_out('real1Dpd')
#endif
!
  return
  end
