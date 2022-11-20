  subroutine cosmoamatdaddd(ipts,jpts,f0,f1,f2,qsij,qsipj,sumAinvipts,sumAinvjpts,ldqneeded,lcosmicd2, &
                            lgrad2,dcosmoA,dcosmoAA,d2qfct,dsegweight,d2segweight,nearsas,nearsasrptr, &
                            dcosmoA2)
!
!  Subroutine adds contributions from a term in the cosmo A matrix to the first and second 
!  derivative matrices, as appropriate.
!  Distributed memory version for parallel second derivatives.
!
!  On input :
!
!  ipts     = first segment pointer, associated with atom i
!  jpts     = second segment pointer, associated with atom j
!  f0      = A matrix term
!  f1(3)   = Cartesian first derivatives of A matrix term with respect to rij
!  f2(6)   = Cartesian second derivatives of A matrix term with respect to rij
!
!   4/17 Created from cosmoamatdadd
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
  use configurations, only : nregionno
  use cosmic
  use current,        only : nsft, nrelf2a
  use derivatives
  use parallel 
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4)                                   :: ipts
  integer(i4)                                   :: jpts
  integer(i4)                                   :: nearsas
  integer(i4)                                   :: nearsasrptr(*)
  logical                                       :: lcosmicd2
  logical                                       :: ldqneeded
  logical                                       :: lgrad2
  real(dp)                                      :: d2qfct
  real(dp)                                      :: dcosmoA(3,npts,*)
  real(dp)                                      :: dcosmoA2(3,npts,*)
  real(dp)                                      :: dcosmoAA(3,nearsas,*)
  real(dp)                                      :: dsegweight(3,maxnearseg,*)
  real(dp)                                      :: d2segweight(3,3,maxnearseg,maxnearseg,*)
  real(dp)                                      :: f0
  real(dp)                                      :: f1(3)
  real(dp)                                      :: f2(6)
  real(dp)                                      :: qsij
  real(dp)                                      :: qsipj
  real(dp)                                      :: sumAinvipts
  real(dp)                                      :: sumAinvjpts
!
!  Local variables
!
  integer(i4)                                   :: i
  integer(i4)                                   :: ia
  integer(i4)                                   :: iloc
  integer(i4)                                   :: in
  integer(i4)                                   :: ix
  integer(i4)                                   :: ixl
  integer(i4)                                   :: iy
  integer(i4)                                   :: iyl
  integer(i4)                                   :: iz
  integer(i4)                                   :: izl
  integer(i4)                                   :: j
  integer(i4)                                   :: ja
  integer(i4)                                   :: jloc
  integer(i4)                                   :: jn
  integer(i4)                                   :: jx
  integer(i4)                                   :: jxl
  integer(i4)                                   :: jy
  integer(i4)                                   :: jyl
  integer(i4)                                   :: jz
  integer(i4)                                   :: jzl
  integer(i4)                                   :: k
  integer(i4)                                   :: kk
  integer(i4)                                   :: kloc
  integer(i4)                                   :: kn
  integer(i4)                                   :: kx
  integer(i4)                                   :: kxl
  integer(i4)                                   :: ky
  integer(i4)                                   :: kyl
  integer(i4)                                   :: kz
  integer(i4)                                   :: kzl
  integer(i4)                                   :: l
  integer(i4)                                   :: ll
  integer(i4)                                   :: lloc
  integer(i4)                                   :: ln
  integer(i4)                                   :: lx
  integer(i4)                                   :: lxl
  integer(i4)                                   :: ly
  integer(i4)                                   :: lyl
  integer(i4)                                   :: lz
  integer(i4)                                   :: lzl
  integer(i4)                                   :: m
  integer(i4)                                   :: mm
  integer(i4)                                   :: mn
  integer(i4)                                   :: nregioni
  integer(i4)                                   :: nregionj
  integer(i4)                                   :: nregionm
  real(dp)                                      :: f0trm
  real(dp)                                      :: qsi
  real(dp)                                      :: qsj
  real(dp)                                      :: qsijw
  real(dp)                                      :: qsipjw
  real(dp)                                      :: smtrm
  real(dp)                                      :: smtrmi
  real(dp)                                      :: smtrmj
  real(dp)                                      :: smtrmp
  real(dp)                                      :: smtrmpi
  real(dp)                                      :: smtrmpj
  real(dp)                                      :: swi
  real(dp)                                      :: swj
  real(dp)                                      :: trm
  real(dp)                                      :: trmi
  real(dp)                                      :: trmj
#ifdef TRACE
  call trace_in('cosmoamatdaddd')
#endif
!
!  Set local constants
!
  i = cosmoatomptr(ipts)
  ia = nrelf2a(i)
  in = nearsasrptr(i)
  qsi = qonsas(ipts)
  swi = segweight(ipts)
  nregioni = nregionno(nsft+ia)
!
  j = cosmoatomptr(jpts)
  ja = nrelf2a(j)
  jn = nearsasrptr(j)
  qsj = qonsas(jpts)
  swj = segweight(jpts)
  nregionj = nregionno(nsft+ja)
!
  qsijw = qsij*swi*swj
  qsipjw = qsipj*swi*swj
!**************************************
!  First derivatives of the A matrix  *
!**************************************
  if (i.ne.j) then
    xdrv(i) = xdrv(i) - qsijw*f1(1)
    ydrv(i) = ydrv(i) - qsijw*f1(2)
    zdrv(i) = zdrv(i) - qsijw*f1(3)
    xdrv(j) = xdrv(j) + qsijw*f1(1)
    ydrv(j) = ydrv(j) + qsijw*f1(2)
    zdrv(j) = zdrv(j) + qsijw*f1(3)
    if (nregioni.ne.nregionj) then
      xregdrv(nregioni) = xregdrv(nregioni) - qsijw*f1(1)
      yregdrv(nregioni) = yregdrv(nregioni) - qsijw*f1(2)
      zregdrv(nregioni) = zregdrv(nregioni) - qsijw*f1(3)
      xregdrv(nregionj) = xregdrv(nregionj) + qsijw*f1(1)
      yregdrv(nregionj) = yregdrv(nregionj) + qsijw*f1(2)
      zregdrv(nregionj) = zregdrv(nregionj) + qsijw*f1(3)
    endif
  endif
!
!  Derivatives with respect to smoothing of segments
!
  if (lsegsmooth) then
    if (ipts.eq.jpts) then
      smtrm = 0.5_dp*qsij*f0
    else
      smtrm = qsij*f0
    endif
    smtrmi = smtrm*swi
    smtrmj = smtrm*swj
    do mm = 1,nnearseg(ipts)
      m = nnearsegptr(mm,ipts)
      xdrv(i) = xdrv(i) - smtrmj*dsegweight(1,mm,ipts)
      ydrv(i) = ydrv(i) - smtrmj*dsegweight(2,mm,ipts)
      zdrv(i) = zdrv(i) - smtrmj*dsegweight(3,mm,ipts)
      xdrv(m) = xdrv(m) + smtrmj*dsegweight(1,mm,ipts)
      ydrv(m) = ydrv(m) + smtrmj*dsegweight(2,mm,ipts)
      zdrv(m) = zdrv(m) + smtrmj*dsegweight(3,mm,ipts)
      nregionm = nregionno(nsft+nrelf2a(m))
      if (nregioni.ne.nregionm) then
        xregdrv(nregioni) = xregdrv(nregioni) - smtrmj*dsegweight(1,mm,ipts)
        yregdrv(nregioni) = yregdrv(nregioni) - smtrmj*dsegweight(2,mm,ipts)
        zregdrv(nregioni) = zregdrv(nregioni) - smtrmj*dsegweight(3,mm,ipts)
        xregdrv(nregionm) = xregdrv(nregionm) + smtrmj*dsegweight(1,mm,ipts)
        yregdrv(nregionm) = yregdrv(nregionm) + smtrmj*dsegweight(2,mm,ipts)
        zregdrv(nregionm) = zregdrv(nregionm) + smtrmj*dsegweight(3,mm,ipts)
      endif
    enddo
    do mm = 1,nnearseg(jpts)
      m = nnearsegptr(mm,jpts)
      xdrv(j) = xdrv(j) - smtrmi*dsegweight(1,mm,jpts)
      ydrv(j) = ydrv(j) - smtrmi*dsegweight(2,mm,jpts)
      zdrv(j) = zdrv(j) - smtrmi*dsegweight(3,mm,jpts)
      xdrv(m) = xdrv(m) + smtrmi*dsegweight(1,mm,jpts)
      ydrv(m) = ydrv(m) + smtrmi*dsegweight(2,mm,jpts)
      zdrv(m) = zdrv(m) + smtrmi*dsegweight(3,mm,jpts)
      nregionm = nregionno(nsft+nrelf2a(m))
      if (nregionj.ne.nregionm) then
        xregdrv(nregionj) = xregdrv(nregionj) - smtrmi*dsegweight(1,mm,jpts)
        yregdrv(nregionj) = yregdrv(nregionj) - smtrmi*dsegweight(2,mm,jpts)
        zregdrv(nregionj) = zregdrv(nregionj) - smtrmi*dsegweight(3,mm,jpts)
        xregdrv(nregionm) = xregdrv(nregionm) + smtrmi*dsegweight(1,mm,jpts)
        yregdrv(nregionm) = yregdrv(nregionm) + smtrmi*dsegweight(2,mm,jpts)
        zregdrv(nregionm) = zregdrv(nregionm) + smtrmi*dsegweight(3,mm,jpts)
      endif
    enddo
  endif
!
  if (ldqneeded) then
!
!  Save first derivatives of A matrix term
!
    if (i.ne.j) then
      dcosmoA(1,ipts,jn) = dcosmoA(1,ipts,jn) - f1(1)*qsj*swi*swj
      dcosmoA(2,ipts,jn) = dcosmoA(2,ipts,jn) - f1(2)*qsj*swi*swj
      dcosmoA(3,ipts,jn) = dcosmoA(3,ipts,jn) - f1(3)*qsj*swi*swj
      dcosmoA(1,jpts,in) = dcosmoA(1,jpts,in) + f1(1)*qsi*swi*swj
      dcosmoA(2,jpts,in) = dcosmoA(2,jpts,in) + f1(2)*qsi*swi*swj
      dcosmoA(3,jpts,in) = dcosmoA(3,jpts,in) + f1(3)*qsi*swi*swj
    endif
!
!  Smoothing contribution to dcosmoA
!
    if (lsegsmooth) then
      if (ipts.eq.jpts) then
        smtrm = 0.5_dp*f0
      else
        smtrm = f0
      endif
      smtrmi = smtrm*qsi*swj
      smtrmj = smtrm*qsj*swj
      do mm = 1,nnearseg(ipts)
        m = nnearsegptr(mm,ipts)
        mn = nearsasrptr(m)
        dcosmoA(1,ipts,mn) = dcosmoA(1,ipts,mn) - smtrmj*dsegweight(1,mm,ipts)
        dcosmoA(2,ipts,mn) = dcosmoA(2,ipts,mn) - smtrmj*dsegweight(2,mm,ipts)
        dcosmoA(3,ipts,mn) = dcosmoA(3,ipts,mn) - smtrmj*dsegweight(3,mm,ipts)
        dcosmoA2(1,jpts,mn) = dcosmoA2(1,jpts,mn) - smtrmi*dsegweight(1,mm,ipts)
        dcosmoA2(2,jpts,mn) = dcosmoA2(2,jpts,mn) - smtrmi*dsegweight(2,mm,ipts)
        dcosmoA2(3,jpts,mn) = dcosmoA2(3,jpts,mn) - smtrmi*dsegweight(3,mm,ipts)
        dcosmoA2(1,jpts,in) = dcosmoA2(1,jpts,in) + smtrmi*dsegweight(1,mm,ipts)
        dcosmoA2(2,jpts,in) = dcosmoA2(2,jpts,in) + smtrmi*dsegweight(2,mm,ipts)
        dcosmoA2(3,jpts,in) = dcosmoA2(3,jpts,in) + smtrmi*dsegweight(3,mm,ipts)
      enddo
      smtrmi = smtrm*qsi*swi
      smtrmj = smtrm*qsj*swi
      do mm = 1,nnearseg(jpts)
        m = nnearsegptr(mm,jpts)
        mn = nearsasrptr(m)
        dcosmoA(1,jpts,mn) = dcosmoA(1,jpts,mn) - smtrmi*dsegweight(1,mm,jpts)
        dcosmoA(2,jpts,mn) = dcosmoA(2,jpts,mn) - smtrmi*dsegweight(2,mm,jpts)
        dcosmoA(3,jpts,mn) = dcosmoA(3,jpts,mn) - smtrmi*dsegweight(3,mm,jpts)
        dcosmoA2(1,ipts,mn) = dcosmoA2(1,ipts,mn) - smtrmj*dsegweight(1,mm,jpts)
        dcosmoA2(2,ipts,mn) = dcosmoA2(2,ipts,mn) - smtrmj*dsegweight(2,mm,jpts)
        dcosmoA2(3,ipts,mn) = dcosmoA2(3,ipts,mn) - smtrmj*dsegweight(3,mm,jpts)
        dcosmoA2(1,ipts,jn) = dcosmoA2(1,ipts,jn) + smtrmj*dsegweight(1,mm,jpts)
        dcosmoA2(2,ipts,jn) = dcosmoA2(2,ipts,jn) + smtrmj*dsegweight(2,mm,jpts)
        dcosmoA2(3,ipts,jn) = dcosmoA2(3,ipts,jn) + smtrmj*dsegweight(3,mm,jpts)
      enddo
    endif
  endif
  if (lgrad2) then
!***********************************
!  Second derivatives of A matrix  *
!***********************************
    ix = 3*(i - 1) + 1
    iy = ix + 1
    iz = ix + 2
    iloc = atom2local(i)
    if (iloc.gt.0) then
      ixl = 3*(iloc - 1) + 1
      iyl = ixl + 1
      izl = ixl + 2
    endif
    jx = 3*(j - 1) + 1
    jy = jx + 1
    jz = jx + 2
    jloc = atom2local(j)
    if (jloc.gt.0) then
      jxl = 3*(jloc - 1) + 1
      jyl = jxl + 1
      jzl = jxl + 2
    endif
!
    if (iloc.gt.0) then
      derv2(jx,ixl) = derv2(jx,ixl) - qsijw*f2(1)
      derv2(jy,ixl) = derv2(jy,ixl) - qsijw*f2(2)
      derv2(jz,ixl) = derv2(jz,ixl) - qsijw*f2(4)
      derv2(jx,iyl) = derv2(jx,iyl) - qsijw*f2(2)
      derv2(jy,iyl) = derv2(jy,iyl) - qsijw*f2(3)
      derv2(jz,iyl) = derv2(jz,iyl) - qsijw*f2(5)
      derv2(jx,izl) = derv2(jx,izl) - qsijw*f2(4)
      derv2(jy,izl) = derv2(jy,izl) - qsijw*f2(5)
      derv2(jz,izl) = derv2(jz,izl) - qsijw*f2(6)
    endif
    if (jloc.gt.0) then
      derv2(ix,jxl) = derv2(ix,jxl) - qsijw*f2(1)
      derv2(iy,jxl) = derv2(iy,jxl) - qsijw*f2(2)
      derv2(iz,jxl) = derv2(iz,jxl) - qsijw*f2(4)
      derv2(ix,jyl) = derv2(ix,jyl) - qsijw*f2(2)
      derv2(iy,jyl) = derv2(iy,jyl) - qsijw*f2(3)
      derv2(iz,jyl) = derv2(iz,jyl) - qsijw*f2(5)
      derv2(ix,jzl) = derv2(ix,jzl) - qsijw*f2(4)
      derv2(iy,jzl) = derv2(iy,jzl) - qsijw*f2(5)
      derv2(iz,jzl) = derv2(iz,jzl) - qsijw*f2(6)
    endif
    if (lcosmicd2) then
      if (iloc.gt.0) then
        derv2(jx,ixl) = derv2(jx,ixl) - d2qfct*qsipjw*f2(1)
        derv2(jy,ixl) = derv2(jy,ixl) - d2qfct*qsipjw*f2(2)
        derv2(jz,ixl) = derv2(jz,ixl) - d2qfct*qsipjw*f2(4)
        derv2(jx,iyl) = derv2(jx,iyl) - d2qfct*qsipjw*f2(2)
        derv2(jy,iyl) = derv2(jy,iyl) - d2qfct*qsipjw*f2(3)
        derv2(jz,iyl) = derv2(jz,iyl) - d2qfct*qsipjw*f2(5)
        derv2(jx,izl) = derv2(jx,izl) - d2qfct*qsipjw*f2(4)
        derv2(jy,izl) = derv2(jy,izl) - d2qfct*qsipjw*f2(5)
        derv2(jz,izl) = derv2(jz,izl) - d2qfct*qsipjw*f2(6)
      endif
      if (jloc.gt.0) then
        derv2(ix,jxl) = derv2(ix,jxl) - d2qfct*qsipjw*f2(1)
        derv2(iy,jxl) = derv2(iy,jxl) - d2qfct*qsipjw*f2(2)
        derv2(iz,jxl) = derv2(iz,jxl) - d2qfct*qsipjw*f2(4)
        derv2(ix,jyl) = derv2(ix,jyl) - d2qfct*qsipjw*f2(2)
        derv2(iy,jyl) = derv2(iy,jyl) - d2qfct*qsipjw*f2(3)
        derv2(iz,jyl) = derv2(iz,jyl) - d2qfct*qsipjw*f2(5)
        derv2(ix,jzl) = derv2(ix,jzl) - d2qfct*qsipjw*f2(4)
        derv2(iy,jzl) = derv2(iy,jzl) - d2qfct*qsipjw*f2(5)
        derv2(iz,jzl) = derv2(iz,jzl) - d2qfct*qsipjw*f2(6)
      endif
!
      dcosmoAA(1,jn,ipts) = dcosmoAA(1,jn,ipts) - f1(1)*sumAinvjpts*swi*swj
      dcosmoAA(2,jn,ipts) = dcosmoAA(2,jn,ipts) - f1(2)*sumAinvjpts*swi*swj
      dcosmoAA(3,jn,ipts) = dcosmoAA(3,jn,ipts) - f1(3)*sumAinvjpts*swi*swj
      dcosmoAA(1,in,jpts) = dcosmoAA(1,in,jpts) + f1(1)*sumAinvipts*swi*swj
      dcosmoAA(2,in,jpts) = dcosmoAA(2,in,jpts) + f1(2)*sumAinvipts*swi*swj
      dcosmoAA(3,in,jpts) = dcosmoAA(3,in,jpts) + f1(3)*sumAinvipts*swi*swj
!
      dcosmoAA(1,jn,jpts) = dcosmoAA(1,jn,jpts) - f1(1)*sumAinvipts*swi*swj
      dcosmoAA(2,jn,jpts) = dcosmoAA(2,jn,jpts) - f1(2)*sumAinvipts*swi*swj
      dcosmoAA(3,jn,jpts) = dcosmoAA(3,jn,jpts) - f1(3)*sumAinvipts*swi*swj
      dcosmoAA(1,in,ipts) = dcosmoAA(1,in,ipts) + f1(1)*sumAinvjpts*swi*swj
      dcosmoAA(2,in,ipts) = dcosmoAA(2,in,ipts) + f1(2)*sumAinvjpts*swi*swj
      dcosmoAA(3,in,ipts) = dcosmoAA(3,in,ipts) + f1(3)*sumAinvjpts*swi*swj
    endif
    if (lsegsmooth) then
!
!  Derivatives due to segment smoothing
!
      if (ipts.eq.jpts) then
        smtrm = 0.5_dp*qsij
        smtrmp = 0.5_dp*qsipj
        f0trm = 0.5_dp*f0
      else
        smtrm = qsij
        smtrmp = qsipj
        f0trm = f0
      endif
      smtrmi = smtrm*swi
      smtrmj = smtrm*swj
      smtrmpi = smtrmp*swi*d2qfct
      smtrmpj = smtrmp*swj*d2qfct
      do kk = 1,nnearseg(ipts)
        k = nnearsegptr(kk,ipts)
        kn = nearsasrptr(k)
!
        if (lcosmicd2) then
          dcosmoAA(1,in,ipts) = dcosmoAA(1,in,ipts) + f0trm*sumAinvjpts*swj*dsegweight(1,kk,ipts)
          dcosmoAA(2,in,ipts) = dcosmoAA(2,in,ipts) + f0trm*sumAinvjpts*swj*dsegweight(2,kk,ipts)
          dcosmoAA(3,in,ipts) = dcosmoAA(3,in,ipts) + f0trm*sumAinvjpts*swj*dsegweight(3,kk,ipts)
          dcosmoAA(1,kn,ipts) = dcosmoAA(1,kn,ipts) - f0trm*sumAinvjpts*swj*dsegweight(1,kk,ipts)
          dcosmoAA(2,kn,ipts) = dcosmoAA(2,kn,ipts) - f0trm*sumAinvjpts*swj*dsegweight(2,kk,ipts)
          dcosmoAA(3,kn,ipts) = dcosmoAA(3,kn,ipts) - f0trm*sumAinvjpts*swj*dsegweight(3,kk,ipts)
          dcosmoAA(1,in,jpts) = dcosmoAA(1,in,jpts) + f0trm*sumAinvipts*swj*dsegweight(1,kk,ipts)
          dcosmoAA(2,in,jpts) = dcosmoAA(2,in,jpts) + f0trm*sumAinvipts*swj*dsegweight(2,kk,ipts)
          dcosmoAA(3,in,jpts) = dcosmoAA(3,in,jpts) + f0trm*sumAinvipts*swj*dsegweight(3,kk,ipts)
          dcosmoAA(1,kn,jpts) = dcosmoAA(1,kn,jpts) - f0trm*sumAinvipts*swj*dsegweight(1,kk,ipts)
          dcosmoAA(2,kn,jpts) = dcosmoAA(2,kn,jpts) - f0trm*sumAinvipts*swj*dsegweight(2,kk,ipts)
          dcosmoAA(3,kn,jpts) = dcosmoAA(3,kn,jpts) - f0trm*sumAinvipts*swj*dsegweight(3,kk,ipts)
        endif
!
        kx = 3*(k - 1) + 1
        ky = kx + 1
        kz = kx + 2
        kloc = atom2local(k)
        if (kloc.gt.0) then
          kxl = 3*(kloc - 1) + 1
          kyl = kxl + 1
          kzl = kxl + 2
        endif
!
        if (i.ne.j) then
          if (iloc.gt.0) then
            derv2(jx,ixl) = derv2(jx,ixl) - smtrmj*f1(1)*dsegweight(1,kk,ipts)
            derv2(jy,ixl) = derv2(jy,ixl) - smtrmj*f1(1)*dsegweight(2,kk,ipts)
            derv2(jz,ixl) = derv2(jz,ixl) - smtrmj*f1(1)*dsegweight(3,kk,ipts)
            derv2(jx,iyl) = derv2(jx,iyl) - smtrmj*f1(2)*dsegweight(1,kk,ipts)
            derv2(jy,iyl) = derv2(jy,iyl) - smtrmj*f1(2)*dsegweight(2,kk,ipts)
            derv2(jz,iyl) = derv2(jz,iyl) - smtrmj*f1(2)*dsegweight(3,kk,ipts)
            derv2(jx,izl) = derv2(jx,izl) - smtrmj*f1(3)*dsegweight(1,kk,ipts)
            derv2(jy,izl) = derv2(jy,izl) - smtrmj*f1(3)*dsegweight(2,kk,ipts)
            derv2(jz,izl) = derv2(jz,izl) - smtrmj*f1(3)*dsegweight(3,kk,ipts)
            if (lcosmicd2) then
              derv2(jx,ixl) = derv2(jx,ixl) - smtrmpj*f1(1)*dsegweight(1,kk,ipts)
              derv2(jy,ixl) = derv2(jy,ixl) - smtrmpj*f1(1)*dsegweight(2,kk,ipts)
              derv2(jz,ixl) = derv2(jz,ixl) - smtrmpj*f1(1)*dsegweight(3,kk,ipts)
              derv2(jx,iyl) = derv2(jx,iyl) - smtrmpj*f1(2)*dsegweight(1,kk,ipts)
              derv2(jy,iyl) = derv2(jy,iyl) - smtrmpj*f1(2)*dsegweight(2,kk,ipts)
              derv2(jz,iyl) = derv2(jz,iyl) - smtrmpj*f1(2)*dsegweight(3,kk,ipts)
              derv2(jx,izl) = derv2(jx,izl) - smtrmpj*f1(3)*dsegweight(1,kk,ipts)
              derv2(jy,izl) = derv2(jy,izl) - smtrmpj*f1(3)*dsegweight(2,kk,ipts)
              derv2(jz,izl) = derv2(jz,izl) - smtrmpj*f1(3)*dsegweight(3,kk,ipts)
            endif
          endif
          if (jloc.gt.0) then
            derv2(ix,jxl) = derv2(ix,jxl) - smtrmj*f1(1)*dsegweight(1,kk,ipts)
            derv2(iy,jxl) = derv2(iy,jxl) - smtrmj*f1(1)*dsegweight(2,kk,ipts)
            derv2(iz,jxl) = derv2(iz,jxl) - smtrmj*f1(1)*dsegweight(3,kk,ipts)
            derv2(ix,jyl) = derv2(ix,jyl) - smtrmj*f1(2)*dsegweight(1,kk,ipts)
            derv2(iy,jyl) = derv2(iy,jyl) - smtrmj*f1(2)*dsegweight(2,kk,ipts)
            derv2(iz,jyl) = derv2(iz,jyl) - smtrmj*f1(2)*dsegweight(3,kk,ipts)
            derv2(ix,jzl) = derv2(ix,jzl) - smtrmj*f1(3)*dsegweight(1,kk,ipts)
            derv2(iy,jzl) = derv2(iy,jzl) - smtrmj*f1(3)*dsegweight(2,kk,ipts)
            derv2(iz,jzl) = derv2(iz,jzl) - smtrmj*f1(3)*dsegweight(3,kk,ipts)
            if (lcosmicd2) then
              derv2(ix,jxl) = derv2(ix,jxl) - smtrmpj*f1(1)*dsegweight(1,kk,ipts)
              derv2(iy,jxl) = derv2(iy,jxl) - smtrmpj*f1(1)*dsegweight(2,kk,ipts)
              derv2(iz,jxl) = derv2(iz,jxl) - smtrmpj*f1(1)*dsegweight(3,kk,ipts)
              derv2(ix,jyl) = derv2(ix,jyl) - smtrmpj*f1(2)*dsegweight(1,kk,ipts)
              derv2(iy,jyl) = derv2(iy,jyl) - smtrmpj*f1(2)*dsegweight(2,kk,ipts)
              derv2(iz,jyl) = derv2(iz,jyl) - smtrmpj*f1(2)*dsegweight(3,kk,ipts)
              derv2(ix,jzl) = derv2(ix,jzl) - smtrmpj*f1(3)*dsegweight(1,kk,ipts)
              derv2(iy,jzl) = derv2(iy,jzl) - smtrmpj*f1(3)*dsegweight(2,kk,ipts)
              derv2(iz,jzl) = derv2(iz,jzl) - smtrmpj*f1(3)*dsegweight(3,kk,ipts)
            endif
          endif
        endif
        if (iloc.gt.0) then
          derv2(kx,ixl) = derv2(kx,ixl) - smtrmj*(f1(1)*dsegweight(1,kk,ipts) - f0*d2segweight(1,1,kk,kk,ipts))
          derv2(ky,ixl) = derv2(ky,ixl) - smtrmj*(f1(1)*dsegweight(2,kk,ipts) - f0*d2segweight(2,1,kk,kk,ipts))
          derv2(kz,ixl) = derv2(kz,ixl) - smtrmj*(f1(1)*dsegweight(3,kk,ipts) - f0*d2segweight(3,1,kk,kk,ipts))
          derv2(kx,iyl) = derv2(kx,iyl) - smtrmj*(f1(2)*dsegweight(1,kk,ipts) - f0*d2segweight(1,2,kk,kk,ipts))
          derv2(ky,iyl) = derv2(ky,iyl) - smtrmj*(f1(2)*dsegweight(2,kk,ipts) - f0*d2segweight(2,2,kk,kk,ipts))
          derv2(kz,iyl) = derv2(kz,iyl) - smtrmj*(f1(2)*dsegweight(3,kk,ipts) - f0*d2segweight(3,2,kk,kk,ipts))
          derv2(kx,izl) = derv2(kx,izl) - smtrmj*(f1(3)*dsegweight(1,kk,ipts) - f0*d2segweight(1,3,kk,kk,ipts))
          derv2(ky,izl) = derv2(ky,izl) - smtrmj*(f1(3)*dsegweight(2,kk,ipts) - f0*d2segweight(2,3,kk,kk,ipts))
          derv2(kz,izl) = derv2(kz,izl) - smtrmj*(f1(3)*dsegweight(3,kk,ipts) - f0*d2segweight(3,3,kk,kk,ipts))
          if (lcosmicd2) then
            derv2(kx,ixl) = derv2(kx,ixl) + smtrmpj*f0*d2segweight(1,1,kk,kk,ipts)
            derv2(ky,ixl) = derv2(ky,ixl) + smtrmpj*f0*d2segweight(2,1,kk,kk,ipts)
            derv2(kz,ixl) = derv2(kz,ixl) + smtrmpj*f0*d2segweight(3,1,kk,kk,ipts)
            derv2(kx,iyl) = derv2(kx,iyl) + smtrmpj*f0*d2segweight(1,2,kk,kk,ipts)
            derv2(ky,iyl) = derv2(ky,iyl) + smtrmpj*f0*d2segweight(2,2,kk,kk,ipts)
            derv2(kz,iyl) = derv2(kz,iyl) + smtrmpj*f0*d2segweight(3,2,kk,kk,ipts)
            derv2(kx,izl) = derv2(kx,izl) + smtrmpj*f0*d2segweight(1,3,kk,kk,ipts)
            derv2(ky,izl) = derv2(ky,izl) + smtrmpj*f0*d2segweight(2,3,kk,kk,ipts)
            derv2(kz,izl) = derv2(kz,izl) + smtrmpj*f0*d2segweight(3,3,kk,kk,ipts)
!
            derv2(kx,ix) = derv2(kx,ix) - smtrmpj*f1(1)*dsegweight(1,kk,ipts)
            derv2(ky,ix) = derv2(ky,ix) - smtrmpj*f1(1)*dsegweight(2,kk,ipts)
            derv2(kz,ix) = derv2(kz,ix) - smtrmpj*f1(1)*dsegweight(3,kk,ipts)
            derv2(kx,iy) = derv2(kx,iy) - smtrmpj*f1(2)*dsegweight(1,kk,ipts)
            derv2(ky,iy) = derv2(ky,iy) - smtrmpj*f1(2)*dsegweight(2,kk,ipts)
            derv2(kz,iy) = derv2(kz,iy) - smtrmpj*f1(2)*dsegweight(3,kk,ipts)
            derv2(kx,iz) = derv2(kx,iz) - smtrmpj*f1(3)*dsegweight(1,kk,ipts)
            derv2(ky,iz) = derv2(ky,iz) - smtrmpj*f1(3)*dsegweight(2,kk,ipts)
            derv2(kz,iz) = derv2(kz,iz) - smtrmpj*f1(3)*dsegweight(3,kk,ipts)
          endif
        endif
        if (kloc.gt.0) then
          derv2(ix,kxl) = derv2(ix,kxl) - smtrmj*(f1(1)*dsegweight(1,kk,ipts) - f0*d2segweight(1,1,kk,kk,ipts))
          derv2(iy,kxl) = derv2(iy,kxl) - smtrmj*(f1(2)*dsegweight(1,kk,ipts) - f0*d2segweight(1,2,kk,kk,ipts))
          derv2(iz,kxl) = derv2(iz,kxl) - smtrmj*(f1(3)*dsegweight(1,kk,ipts) - f0*d2segweight(1,3,kk,kk,ipts))
          derv2(ix,kyl) = derv2(ix,kyl) - smtrmj*(f1(1)*dsegweight(2,kk,ipts) - f0*d2segweight(2,1,kk,kk,ipts))
          derv2(iy,kyl) = derv2(iy,kyl) - smtrmj*(f1(2)*dsegweight(2,kk,ipts) - f0*d2segweight(2,2,kk,kk,ipts))
          derv2(iz,kyl) = derv2(iz,kyl) - smtrmj*(f1(3)*dsegweight(2,kk,ipts) - f0*d2segweight(2,3,kk,kk,ipts))
          derv2(ix,kzl) = derv2(ix,kzl) - smtrmj*(f1(1)*dsegweight(3,kk,ipts) - f0*d2segweight(3,1,kk,kk,ipts))
          derv2(iy,kzl) = derv2(iy,kzl) - smtrmj*(f1(2)*dsegweight(3,kk,ipts) - f0*d2segweight(3,2,kk,kk,ipts))
          derv2(iz,kzl) = derv2(iz,kzl) - smtrmj*(f1(3)*dsegweight(3,kk,ipts) - f0*d2segweight(3,3,kk,kk,ipts))
          if (lcosmicd2) then
            derv2(ix,kxl) = derv2(ix,kxl) + smtrmpj*f0*d2segweight(1,1,kk,kk,ipts)
            derv2(iy,kxl) = derv2(iy,kxl) + smtrmpj*f0*d2segweight(1,2,kk,kk,ipts)
            derv2(iz,kxl) = derv2(iz,kxl) + smtrmpj*f0*d2segweight(1,3,kk,kk,ipts)
            derv2(ix,kyl) = derv2(ix,kyl) + smtrmpj*f0*d2segweight(2,1,kk,kk,ipts)
            derv2(iy,kyl) = derv2(iy,kyl) + smtrmpj*f0*d2segweight(2,2,kk,kk,ipts)
            derv2(iz,kyl) = derv2(iz,kyl) + smtrmpj*f0*d2segweight(2,3,kk,kk,ipts)
            derv2(ix,kzl) = derv2(ix,kzl) + smtrmpj*f0*d2segweight(3,1,kk,kk,ipts)
            derv2(iy,kzl) = derv2(iy,kzl) + smtrmpj*f0*d2segweight(3,2,kk,kk,ipts)
            derv2(iz,kzl) = derv2(iz,kzl) + smtrmpj*f0*d2segweight(3,3,kk,kk,ipts)
!
            derv2(ix,kxl) = derv2(ix,kxl) - smtrmpj*f1(1)*dsegweight(1,kk,ipts)
            derv2(iy,kxl) = derv2(iy,kxl) - smtrmpj*f1(2)*dsegweight(1,kk,ipts)
            derv2(iz,kxl) = derv2(iz,kxl) - smtrmpj*f1(3)*dsegweight(1,kk,ipts)
            derv2(ix,kyl) = derv2(ix,kyl) - smtrmpj*f1(1)*dsegweight(2,kk,ipts)
            derv2(iy,kyl) = derv2(iy,kyl) - smtrmpj*f1(2)*dsegweight(2,kk,ipts)
            derv2(iz,kyl) = derv2(iz,kyl) - smtrmpj*f1(3)*dsegweight(2,kk,ipts)
            derv2(ix,kzl) = derv2(ix,kzl) - smtrmpj*f1(1)*dsegweight(3,kk,ipts)
            derv2(iy,kzl) = derv2(iy,kzl) - smtrmpj*f1(2)*dsegweight(3,kk,ipts)
            derv2(iz,kzl) = derv2(iz,kzl) - smtrmpj*f1(3)*dsegweight(3,kk,ipts)
          endif
        endif
        if (jloc.gt.0) then
          derv2(kx,jxl) = derv2(kx,jxl) + smtrmj*f1(1)*dsegweight(1,kk,ipts)
          derv2(ky,jxl) = derv2(ky,jxl) + smtrmj*f1(1)*dsegweight(2,kk,ipts)
          derv2(kz,jxl) = derv2(kz,jxl) + smtrmj*f1(1)*dsegweight(3,kk,ipts)
          derv2(kx,jyl) = derv2(kx,jyl) + smtrmj*f1(2)*dsegweight(1,kk,ipts)
          derv2(ky,jyl) = derv2(ky,jyl) + smtrmj*f1(2)*dsegweight(2,kk,ipts)
          derv2(kz,jyl) = derv2(kz,jyl) + smtrmj*f1(2)*dsegweight(3,kk,ipts)
          derv2(kx,jzl) = derv2(kx,jzl) + smtrmj*f1(3)*dsegweight(1,kk,ipts)
          derv2(ky,jzl) = derv2(ky,jzl) + smtrmj*f1(3)*dsegweight(2,kk,ipts)
          derv2(kz,jzl) = derv2(kz,jzl) + smtrmj*f1(3)*dsegweight(3,kk,ipts)
          if (lcosmicd2) then
            derv2(kx,jxl) = derv2(kx,jxl) + smtrmpj*f1(1)*dsegweight(1,kk,ipts)
            derv2(ky,jxl) = derv2(ky,jxl) + smtrmpj*f1(1)*dsegweight(2,kk,ipts)
            derv2(kz,jxl) = derv2(kz,jxl) + smtrmpj*f1(1)*dsegweight(3,kk,ipts)
            derv2(kx,jyl) = derv2(kx,jyl) + smtrmpj*f1(2)*dsegweight(1,kk,ipts)
            derv2(ky,jyl) = derv2(ky,jyl) + smtrmpj*f1(2)*dsegweight(2,kk,ipts)
            derv2(kz,jyl) = derv2(kz,jyl) + smtrmpj*f1(2)*dsegweight(3,kk,ipts)
            derv2(kx,jzl) = derv2(kx,jzl) + smtrmpj*f1(3)*dsegweight(1,kk,ipts)
            derv2(ky,jzl) = derv2(ky,jzl) + smtrmpj*f1(3)*dsegweight(2,kk,ipts)
            derv2(kz,jzl) = derv2(kz,jzl) + smtrmpj*f1(3)*dsegweight(3,kk,ipts)
          endif
        endif
        if (kloc.gt.0) then
          derv2(jx,kxl) = derv2(jx,kxl) + smtrmj*f1(1)*dsegweight(1,kk,ipts)
          derv2(jy,kxl) = derv2(jy,kxl) + smtrmj*f1(2)*dsegweight(1,kk,ipts)
          derv2(jz,kxl) = derv2(jz,kxl) + smtrmj*f1(3)*dsegweight(1,kk,ipts)
          derv2(jx,kyl) = derv2(jx,kyl) + smtrmj*f1(1)*dsegweight(2,kk,ipts)
          derv2(jy,kyl) = derv2(jy,kyl) + smtrmj*f1(2)*dsegweight(2,kk,ipts)
          derv2(jz,kyl) = derv2(jz,kyl) + smtrmj*f1(3)*dsegweight(2,kk,ipts)
          derv2(jx,kzl) = derv2(jx,kzl) + smtrmj*f1(1)*dsegweight(3,kk,ipts)
          derv2(jy,kzl) = derv2(jy,kzl) + smtrmj*f1(2)*dsegweight(3,kk,ipts)
          derv2(jz,kzl) = derv2(jz,kzl) + smtrmj*f1(3)*dsegweight(3,kk,ipts)
          if (lcosmicd2) then
            derv2(jx,kxl) = derv2(jx,kxl) + smtrmpj*f1(1)*dsegweight(1,kk,ipts)
            derv2(jy,kxl) = derv2(jy,kxl) + smtrmpj*f1(2)*dsegweight(1,kk,ipts)
            derv2(jz,kxl) = derv2(jz,kxl) + smtrmpj*f1(3)*dsegweight(1,kk,ipts)
            derv2(jx,kyl) = derv2(jx,kyl) + smtrmpj*f1(1)*dsegweight(2,kk,ipts)
            derv2(jy,kyl) = derv2(jy,kyl) + smtrmpj*f1(2)*dsegweight(2,kk,ipts)
            derv2(jz,kyl) = derv2(jz,kyl) + smtrmpj*f1(3)*dsegweight(2,kk,ipts)
            derv2(jx,kzl) = derv2(jx,kzl) + smtrmpj*f1(1)*dsegweight(3,kk,ipts)
            derv2(jy,kzl) = derv2(jy,kzl) + smtrmpj*f1(2)*dsegweight(3,kk,ipts)
            derv2(jz,kzl) = derv2(jz,kzl) + smtrmpj*f1(3)*dsegweight(3,kk,ipts)
          endif
        endif
      enddo
      do ll = 1,nnearseg(jpts)
        l = nnearsegptr(ll,jpts)
        ln = nearsasrptr(l)
!
        if (lcosmicd2) then
          dcosmoAA(1,jn,jpts) = dcosmoAA(1,jn,jpts) + f0trm*sumAinvipts*swi*dsegweight(1,ll,jpts)
          dcosmoAA(2,jn,jpts) = dcosmoAA(2,jn,jpts) + f0trm*sumAinvipts*swi*dsegweight(2,ll,jpts)
          dcosmoAA(3,jn,jpts) = dcosmoAA(3,jn,jpts) + f0trm*sumAinvipts*swi*dsegweight(3,ll,jpts)
          dcosmoAA(1,ln,jpts) = dcosmoAA(1,ln,jpts) - f0trm*sumAinvipts*swi*dsegweight(1,ll,jpts)
          dcosmoAA(2,ln,jpts) = dcosmoAA(2,ln,jpts) - f0trm*sumAinvipts*swi*dsegweight(2,ll,jpts)
          dcosmoAA(3,ln,jpts) = dcosmoAA(3,ln,jpts) - f0trm*sumAinvipts*swi*dsegweight(3,ll,jpts)
          dcosmoAA(1,jn,ipts) = dcosmoAA(1,jn,ipts) + f0trm*sumAinvjpts*swi*dsegweight(1,ll,jpts)
          dcosmoAA(2,jn,ipts) = dcosmoAA(2,jn,ipts) + f0trm*sumAinvjpts*swi*dsegweight(2,ll,jpts)
          dcosmoAA(3,jn,ipts) = dcosmoAA(3,jn,ipts) + f0trm*sumAinvjpts*swi*dsegweight(3,ll,jpts)
          dcosmoAA(1,ln,ipts) = dcosmoAA(1,ln,ipts) - f0trm*sumAinvjpts*swi*dsegweight(1,ll,jpts)
          dcosmoAA(2,ln,ipts) = dcosmoAA(2,ln,ipts) - f0trm*sumAinvjpts*swi*dsegweight(2,ll,jpts)
          dcosmoAA(3,ln,ipts) = dcosmoAA(3,ln,ipts) - f0trm*sumAinvjpts*swi*dsegweight(3,ll,jpts)
        endif
!
        lx = 3*(l - 1) + 1
        ly = lx + 1
        lz = lx + 2
        lloc = atom2local(l)
        if (lloc.gt.0) then
          lxl = 3*(lloc - 1) + 1
          lyl = lxl + 1
          lzl = lxl + 2
        endif
!
        if (i.ne.j) then
          if (iloc.gt.0) then
            derv2(jx,ixl) = derv2(jx,ixl) + smtrmi*f1(1)*dsegweight(1,ll,jpts)
            derv2(jy,ixl) = derv2(jy,ixl) + smtrmi*f1(2)*dsegweight(1,ll,jpts)
            derv2(jz,ixl) = derv2(jz,ixl) + smtrmi*f1(3)*dsegweight(1,ll,jpts)
            derv2(jx,iyl) = derv2(jx,iyl) + smtrmi*f1(1)*dsegweight(2,ll,jpts)
            derv2(jy,iyl) = derv2(jy,iyl) + smtrmi*f1(2)*dsegweight(2,ll,jpts)
            derv2(jz,iyl) = derv2(jz,iyl) + smtrmi*f1(3)*dsegweight(2,ll,jpts)
            derv2(jx,izl) = derv2(jx,izl) + smtrmi*f1(1)*dsegweight(3,ll,jpts)
            derv2(jy,izl) = derv2(jy,izl) + smtrmi*f1(2)*dsegweight(3,ll,jpts)
            derv2(jz,izl) = derv2(jz,izl) + smtrmi*f1(3)*dsegweight(3,ll,jpts)
            if (lcosmicd2) then
              derv2(jx,ixl) = derv2(jx,ixl) + smtrmpi*f1(1)*dsegweight(1,ll,jpts)
              derv2(jy,ixl) = derv2(jy,ixl) + smtrmpi*f1(2)*dsegweight(1,ll,jpts)
              derv2(jz,ixl) = derv2(jz,ixl) + smtrmpi*f1(3)*dsegweight(1,ll,jpts)
              derv2(jx,iyl) = derv2(jx,iyl) + smtrmpi*f1(1)*dsegweight(2,ll,jpts)
              derv2(jy,iyl) = derv2(jy,iyl) + smtrmpi*f1(2)*dsegweight(2,ll,jpts)
              derv2(jz,iyl) = derv2(jz,iyl) + smtrmpi*f1(3)*dsegweight(2,ll,jpts)
              derv2(jx,izl) = derv2(jx,izl) + smtrmpi*f1(1)*dsegweight(3,ll,jpts)
              derv2(jy,izl) = derv2(jy,izl) + smtrmpi*f1(2)*dsegweight(3,ll,jpts)
              derv2(jz,izl) = derv2(jz,izl) + smtrmpi*f1(3)*dsegweight(3,ll,jpts)
            endif
          endif
          if (jloc.gt.0) then
            derv2(ix,jxl) = derv2(ix,jxl) + smtrmi*f1(1)*dsegweight(1,ll,jpts)
            derv2(iy,jxl) = derv2(iy,jxl) + smtrmi*f1(2)*dsegweight(1,ll,jpts)
            derv2(iz,jxl) = derv2(iz,jxl) + smtrmi*f1(3)*dsegweight(1,ll,jpts)
            derv2(ix,jyl) = derv2(ix,jyl) + smtrmi*f1(1)*dsegweight(2,ll,jpts)
            derv2(iy,jyl) = derv2(iy,jyl) + smtrmi*f1(2)*dsegweight(2,ll,jpts)
            derv2(iz,jyl) = derv2(iz,jyl) + smtrmi*f1(3)*dsegweight(2,ll,jpts)
            derv2(ix,jzl) = derv2(ix,jzl) + smtrmi*f1(1)*dsegweight(3,ll,jpts)
            derv2(iy,jzl) = derv2(iy,jzl) + smtrmi*f1(2)*dsegweight(3,ll,jpts)
            derv2(iz,jzl) = derv2(iz,jzl) + smtrmi*f1(3)*dsegweight(3,ll,jpts)
            if (lcosmicd2) then
              derv2(ix,jxl) = derv2(ix,jxl) + smtrmpi*f1(1)*dsegweight(1,ll,jpts)
              derv2(iy,jxl) = derv2(iy,jxl) + smtrmpi*f1(2)*dsegweight(1,ll,jpts)
              derv2(iz,jxl) = derv2(iz,jxl) + smtrmpi*f1(3)*dsegweight(1,ll,jpts)
              derv2(ix,jyl) = derv2(ix,jyl) + smtrmpi*f1(1)*dsegweight(2,ll,jpts)
              derv2(iy,jyl) = derv2(iy,jyl) + smtrmpi*f1(2)*dsegweight(2,ll,jpts)
              derv2(iz,jyl) = derv2(iz,jyl) + smtrmpi*f1(3)*dsegweight(2,ll,jpts)
              derv2(ix,jzl) = derv2(ix,jzl) + smtrmpi*f1(1)*dsegweight(3,ll,jpts)
              derv2(iy,jzl) = derv2(iy,jzl) + smtrmpi*f1(2)*dsegweight(3,ll,jpts)
              derv2(iz,jzl) = derv2(iz,jzl) + smtrmpi*f1(3)*dsegweight(3,ll,jpts)
            endif
          endif
        endif
        if (iloc.gt.0) then
          derv2(lx,ixl) = derv2(lx,ixl) - smtrmi*f1(1)*dsegweight(1,ll,jpts)
          derv2(ly,ixl) = derv2(ly,ixl) - smtrmi*f1(1)*dsegweight(2,ll,jpts)
          derv2(lz,ixl) = derv2(lz,ixl) - smtrmi*f1(1)*dsegweight(3,ll,jpts)
          derv2(lx,iyl) = derv2(lx,iyl) - smtrmi*f1(2)*dsegweight(1,ll,jpts)
          derv2(ly,iyl) = derv2(ly,iyl) - smtrmi*f1(2)*dsegweight(2,ll,jpts)
          derv2(lz,iyl) = derv2(lz,iyl) - smtrmi*f1(2)*dsegweight(3,ll,jpts)
          derv2(lx,izl) = derv2(lx,izl) - smtrmi*f1(3)*dsegweight(1,ll,jpts)
          derv2(ly,izl) = derv2(ly,izl) - smtrmi*f1(3)*dsegweight(2,ll,jpts)
          derv2(lz,izl) = derv2(lz,izl) - smtrmi*f1(3)*dsegweight(3,ll,jpts)
          if (lcosmicd2) then
            derv2(lx,ixl) = derv2(lx,ixl) - smtrmpi*f1(1)*dsegweight(1,ll,jpts)
            derv2(ly,ixl) = derv2(ly,ixl) - smtrmpi*f1(1)*dsegweight(2,ll,jpts)
            derv2(lz,ixl) = derv2(lz,ixl) - smtrmpi*f1(1)*dsegweight(3,ll,jpts)
            derv2(lx,iyl) = derv2(lx,iyl) - smtrmpi*f1(2)*dsegweight(1,ll,jpts)
            derv2(ly,iyl) = derv2(ly,iyl) - smtrmpi*f1(2)*dsegweight(2,ll,jpts)
            derv2(lz,iyl) = derv2(lz,iyl) - smtrmpi*f1(2)*dsegweight(3,ll,jpts)
            derv2(lx,izl) = derv2(lx,izl) - smtrmpi*f1(3)*dsegweight(1,ll,jpts)
            derv2(ly,izl) = derv2(ly,izl) - smtrmpi*f1(3)*dsegweight(2,ll,jpts)
            derv2(lz,izl) = derv2(lz,izl) - smtrmpi*f1(3)*dsegweight(3,ll,jpts)
          endif
        endif
        if (lloc.gt.0) then
          derv2(ix,lxl) = derv2(ix,lxl) - smtrmi*f1(1)*dsegweight(1,ll,jpts)
          derv2(iy,lxl) = derv2(iy,lxl) - smtrmi*f1(2)*dsegweight(1,ll,jpts)
          derv2(iz,lxl) = derv2(iz,lxl) - smtrmi*f1(3)*dsegweight(1,ll,jpts)
          derv2(ix,lyl) = derv2(ix,lyl) - smtrmi*f1(1)*dsegweight(2,ll,jpts)
          derv2(iy,lyl) = derv2(iy,lyl) - smtrmi*f1(2)*dsegweight(2,ll,jpts)
          derv2(iz,lyl) = derv2(iz,lyl) - smtrmi*f1(3)*dsegweight(2,ll,jpts)
          derv2(ix,lzl) = derv2(ix,lzl) - smtrmi*f1(1)*dsegweight(3,ll,jpts)
          derv2(iy,lzl) = derv2(iy,lzl) - smtrmi*f1(2)*dsegweight(3,ll,jpts)
          derv2(iz,lzl) = derv2(iz,lzl) - smtrmi*f1(3)*dsegweight(3,ll,jpts)
          if (lcosmicd2) then
            derv2(ix,lxl) = derv2(ix,lxl) - smtrmpi*f1(1)*dsegweight(1,ll,jpts)
            derv2(iy,lxl) = derv2(iy,lxl) - smtrmpi*f1(2)*dsegweight(1,ll,jpts)
            derv2(iz,lxl) = derv2(iz,lxl) - smtrmpi*f1(3)*dsegweight(1,ll,jpts)
            derv2(ix,lyl) = derv2(ix,lyl) - smtrmpi*f1(1)*dsegweight(2,ll,jpts)
            derv2(iy,lyl) = derv2(iy,lyl) - smtrmpi*f1(2)*dsegweight(2,ll,jpts)
            derv2(iz,lyl) = derv2(iz,lyl) - smtrmpi*f1(3)*dsegweight(2,ll,jpts)
            derv2(ix,lzl) = derv2(ix,lzl) - smtrmpi*f1(1)*dsegweight(3,ll,jpts)
            derv2(iy,lzl) = derv2(iy,lzl) - smtrmpi*f1(2)*dsegweight(3,ll,jpts)
            derv2(iz,lzl) = derv2(iz,lzl) - smtrmpi*f1(3)*dsegweight(3,ll,jpts)
          endif
        endif
        if (jloc.gt.0) then
          derv2(lx,jxl) = derv2(lx,jxl) + smtrmi*(f1(1)*dsegweight(1,ll,jpts) + f0*d2segweight(1,1,ll,ll,jpts))
          derv2(ly,jxl) = derv2(ly,jxl) + smtrmi*(f1(1)*dsegweight(2,ll,jpts) + f0*d2segweight(2,1,ll,ll,jpts))
          derv2(lz,jxl) = derv2(lz,jxl) + smtrmi*(f1(1)*dsegweight(3,ll,jpts) + f0*d2segweight(3,1,ll,ll,jpts))
          derv2(lx,jyl) = derv2(lx,jyl) + smtrmi*(f1(2)*dsegweight(1,ll,jpts) + f0*d2segweight(1,2,ll,ll,jpts))
          derv2(ly,jyl) = derv2(ly,jyl) + smtrmi*(f1(2)*dsegweight(2,ll,jpts) + f0*d2segweight(2,2,ll,ll,jpts))
          derv2(lz,jyl) = derv2(lz,jyl) + smtrmi*(f1(2)*dsegweight(3,ll,jpts) + f0*d2segweight(3,2,ll,ll,jpts))
          derv2(lx,jzl) = derv2(lx,jzl) + smtrmi*(f1(3)*dsegweight(1,ll,jpts) + f0*d2segweight(1,3,ll,ll,jpts))
          derv2(ly,jzl) = derv2(ly,jzl) + smtrmi*(f1(3)*dsegweight(2,ll,jpts) + f0*d2segweight(2,3,ll,ll,jpts))
          derv2(lz,jzl) = derv2(lz,jzl) + smtrmi*(f1(3)*dsegweight(3,ll,jpts) + f0*d2segweight(3,3,ll,ll,jpts))
          if (lcosmicd2) then
            derv2(lx,jxl) = derv2(lx,jxl) + smtrmpi*f0*d2segweight(1,1,ll,ll,jpts)
            derv2(ly,jxl) = derv2(ly,jxl) + smtrmpi*f0*d2segweight(2,1,ll,ll,jpts)
            derv2(lz,jxl) = derv2(lz,jxl) + smtrmpi*f0*d2segweight(3,1,ll,ll,jpts)
            derv2(lx,jyl) = derv2(lx,jyl) + smtrmpi*f0*d2segweight(1,2,ll,ll,jpts)
            derv2(ly,jyl) = derv2(ly,jyl) + smtrmpi*f0*d2segweight(2,2,ll,ll,jpts)
            derv2(lz,jyl) = derv2(lz,jyl) + smtrmpi*f0*d2segweight(3,2,ll,ll,jpts)
            derv2(lx,jzl) = derv2(lx,jzl) + smtrmpi*f0*d2segweight(1,3,ll,ll,jpts)
            derv2(ly,jzl) = derv2(ly,jzl) + smtrmpi*f0*d2segweight(2,3,ll,ll,jpts)
            derv2(lz,jzl) = derv2(lz,jzl) + smtrmpi*f0*d2segweight(3,3,ll,ll,jpts)
!
            derv2(lx,jxl) = derv2(lx,jxl) + smtrmpi*f1(1)*dsegweight(1,ll,jpts)
            derv2(ly,jxl) = derv2(ly,jxl) + smtrmpi*f1(1)*dsegweight(2,ll,jpts)
            derv2(lz,jxl) = derv2(lz,jxl) + smtrmpi*f1(1)*dsegweight(3,ll,jpts)
            derv2(lx,jyl) = derv2(lx,jyl) + smtrmpi*f1(2)*dsegweight(1,ll,jpts)
            derv2(ly,jyl) = derv2(ly,jyl) + smtrmpi*f1(2)*dsegweight(2,ll,jpts)
            derv2(lz,jyl) = derv2(lz,jyl) + smtrmpi*f1(2)*dsegweight(3,ll,jpts)
            derv2(lx,jzl) = derv2(lx,jzl) + smtrmpi*f1(3)*dsegweight(1,ll,jpts)
            derv2(ly,jzl) = derv2(ly,jzl) + smtrmpi*f1(3)*dsegweight(2,ll,jpts)
            derv2(lz,jzl) = derv2(lz,jzl) + smtrmpi*f1(3)*dsegweight(3,ll,jpts)
          endif
        endif
        if (lloc.gt.0) then
          derv2(jx,lxl) = derv2(jx,lxl) + smtrmi*(f1(1)*dsegweight(1,ll,jpts) + f0*d2segweight(1,1,ll,ll,jpts))
          derv2(jy,lxl) = derv2(jy,lxl) + smtrmi*(f1(2)*dsegweight(1,ll,jpts) + f0*d2segweight(1,2,ll,ll,jpts))
          derv2(jz,lxl) = derv2(jz,lxl) + smtrmi*(f1(3)*dsegweight(1,ll,jpts) + f0*d2segweight(1,3,ll,ll,jpts))
          derv2(jx,lyl) = derv2(jx,lyl) + smtrmi*(f1(1)*dsegweight(2,ll,jpts) + f0*d2segweight(2,1,ll,ll,jpts))
          derv2(jy,lyl) = derv2(jy,lyl) + smtrmi*(f1(2)*dsegweight(2,ll,jpts) + f0*d2segweight(2,2,ll,ll,jpts))
          derv2(jz,lyl) = derv2(jz,lyl) + smtrmi*(f1(3)*dsegweight(2,ll,jpts) + f0*d2segweight(2,3,ll,ll,jpts))
          derv2(jx,lzl) = derv2(jx,lzl) + smtrmi*(f1(1)*dsegweight(3,ll,jpts) + f0*d2segweight(3,1,ll,ll,jpts))
          derv2(jy,lzl) = derv2(jy,lzl) + smtrmi*(f1(2)*dsegweight(3,ll,jpts) + f0*d2segweight(3,2,ll,ll,jpts))
          derv2(jz,lzl) = derv2(jz,lzl) + smtrmi*(f1(3)*dsegweight(3,ll,jpts) + f0*d2segweight(3,3,ll,ll,jpts))
          if (lcosmicd2) then
            derv2(jx,lxl) = derv2(jx,lxl) + smtrmpi*f0*d2segweight(1,1,ll,ll,jpts)
            derv2(jy,lxl) = derv2(jy,lxl) + smtrmpi*f0*d2segweight(1,2,ll,ll,jpts)
            derv2(jz,lxl) = derv2(jz,lxl) + smtrmpi*f0*d2segweight(1,3,ll,ll,jpts)
            derv2(jx,lyl) = derv2(jx,lyl) + smtrmpi*f0*d2segweight(2,1,ll,ll,jpts)
            derv2(jy,lyl) = derv2(jy,lyl) + smtrmpi*f0*d2segweight(2,2,ll,ll,jpts)
            derv2(jz,lyl) = derv2(jz,lyl) + smtrmpi*f0*d2segweight(2,3,ll,ll,jpts)
            derv2(jx,lzl) = derv2(jx,lzl) + smtrmpi*f0*d2segweight(3,1,ll,ll,jpts)
            derv2(jy,lzl) = derv2(jy,lzl) + smtrmpi*f0*d2segweight(3,2,ll,ll,jpts)
            derv2(jz,lzl) = derv2(jz,lzl) + smtrmpi*f0*d2segweight(3,3,ll,ll,jpts)
!
            derv2(jx,lxl) = derv2(jx,lxl) + smtrmpi*f1(1)*dsegweight(1,ll,jpts)
            derv2(jy,lxl) = derv2(jy,lxl) + smtrmpi*f1(2)*dsegweight(1,ll,jpts)
            derv2(jz,lxl) = derv2(jz,lxl) + smtrmpi*f1(3)*dsegweight(1,ll,jpts)
            derv2(jx,lyl) = derv2(jx,lyl) + smtrmpi*f1(1)*dsegweight(2,ll,jpts)
            derv2(jy,lyl) = derv2(jy,lyl) + smtrmpi*f1(2)*dsegweight(2,ll,jpts)
            derv2(jz,lyl) = derv2(jz,lyl) + smtrmpi*f1(3)*dsegweight(2,ll,jpts)
            derv2(jx,lzl) = derv2(jx,lzl) + smtrmpi*f1(1)*dsegweight(3,ll,jpts)
            derv2(jy,lzl) = derv2(jy,lzl) + smtrmpi*f1(2)*dsegweight(3,ll,jpts)
            derv2(jz,lzl) = derv2(jz,lzl) + smtrmpi*f1(3)*dsegweight(3,ll,jpts)
          endif
        endif
      enddo
!
!  Combined loops over weighting atoms - one for i and one for j
!
      if (ipts.eq.jpts) then
        smtrm = 0.5_dp*qsij*f0
        smtrmp = 0.5_dp*qsipj*f0*d2qfct
      else
        smtrm = qsij*f0
        smtrmp = qsipj*f0*d2qfct
      endif
      if (lcosmicd2) then
        trm = smtrm + smtrmp
      else
        trm = smtrm 
      endif
      do kk = 1,nnearseg(ipts)
        k = nnearsegptr(kk,ipts)
        kx = 3*(k - 1) + 1
        ky = kx + 1
        kz = kx + 2
        kloc = atom2local(k)
        if (kloc.gt.0) then
          kxl = 3*(kloc - 1) + 1
          kyl = kxl + 1
          kzl = kxl + 2
        endif
        do ll = 1,nnearseg(jpts)
          l = nnearsegptr(ll,jpts)
          lx = 3*(l - 1) + 1
          ly = lx + 1
          lz = lx + 2
          lloc = atom2local(l)
          if (lloc.gt.0) then
            lxl = 3*(lloc - 1) + 1
            lyl = lxl + 1
            lzl = lxl + 2
          endif
          if (i.ne.j) then
            if (iloc.gt.0) then
              derv2(jx,ixl) = derv2(jx,ixl) + trm*dsegweight(1,ll,jpts)*dsegweight(1,kk,ipts)
              derv2(jy,ixl) = derv2(jy,ixl) + trm*dsegweight(1,ll,jpts)*dsegweight(2,kk,ipts)
              derv2(jz,ixl) = derv2(jz,ixl) + trm*dsegweight(1,ll,jpts)*dsegweight(3,kk,ipts)
              derv2(jx,iyl) = derv2(jx,iyl) + trm*dsegweight(2,ll,jpts)*dsegweight(1,kk,ipts)
              derv2(jy,iyl) = derv2(jy,iyl) + trm*dsegweight(2,ll,jpts)*dsegweight(2,kk,ipts)
              derv2(jz,iyl) = derv2(jz,iyl) + trm*dsegweight(2,ll,jpts)*dsegweight(3,kk,ipts)
              derv2(jx,izl) = derv2(jx,izl) + trm*dsegweight(3,ll,jpts)*dsegweight(1,kk,ipts)
              derv2(jy,izl) = derv2(jy,izl) + trm*dsegweight(3,ll,jpts)*dsegweight(2,kk,ipts)
              derv2(jz,izl) = derv2(jz,izl) + trm*dsegweight(3,ll,jpts)*dsegweight(3,kk,ipts)
            endif
            if (jloc.gt.0) then
              derv2(ix,jxl) = derv2(ix,jxl) + trm*dsegweight(1,ll,jpts)*dsegweight(1,kk,ipts)
              derv2(iy,jxl) = derv2(iy,jxl) + trm*dsegweight(1,ll,jpts)*dsegweight(2,kk,ipts)
              derv2(iz,jxl) = derv2(iz,jxl) + trm*dsegweight(1,ll,jpts)*dsegweight(3,kk,ipts)
              derv2(ix,jyl) = derv2(ix,jyl) + trm*dsegweight(2,ll,jpts)*dsegweight(1,kk,ipts)
              derv2(iy,jyl) = derv2(iy,jyl) + trm*dsegweight(2,ll,jpts)*dsegweight(2,kk,ipts)
              derv2(iz,jyl) = derv2(iz,jyl) + trm*dsegweight(2,ll,jpts)*dsegweight(3,kk,ipts)
              derv2(ix,jzl) = derv2(ix,jzl) + trm*dsegweight(3,ll,jpts)*dsegweight(1,kk,ipts)
              derv2(iy,jzl) = derv2(iy,jzl) + trm*dsegweight(3,ll,jpts)*dsegweight(2,kk,ipts)
              derv2(iz,jzl) = derv2(iz,jzl) + trm*dsegweight(3,ll,jpts)*dsegweight(3,kk,ipts)
            endif
          endif
          if (iloc.gt.0) then
            derv2(lx,ixl) = derv2(lx,ixl) - trm*dsegweight(1,ll,jpts)*dsegweight(1,kk,ipts)
            derv2(ly,ixl) = derv2(ly,ixl) - trm*dsegweight(2,ll,jpts)*dsegweight(1,kk,ipts)
            derv2(lz,ixl) = derv2(lz,ixl) - trm*dsegweight(3,ll,jpts)*dsegweight(1,kk,ipts)
            derv2(lx,iyl) = derv2(lx,iyl) - trm*dsegweight(1,ll,jpts)*dsegweight(2,kk,ipts)
            derv2(ly,iyl) = derv2(ly,iyl) - trm*dsegweight(2,ll,jpts)*dsegweight(2,kk,ipts)
            derv2(lz,iyl) = derv2(lz,iyl) - trm*dsegweight(3,ll,jpts)*dsegweight(2,kk,ipts)
            derv2(lx,izl) = derv2(lx,izl) - trm*dsegweight(1,ll,jpts)*dsegweight(3,kk,ipts)
            derv2(ly,izl) = derv2(ly,izl) - trm*dsegweight(2,ll,jpts)*dsegweight(3,kk,ipts)
            derv2(lz,izl) = derv2(lz,izl) - trm*dsegweight(3,ll,jpts)*dsegweight(3,kk,ipts)
          endif
          if (lloc.gt.0) then
            derv2(ix,lxl) = derv2(ix,lxl) - trm*dsegweight(1,ll,jpts)*dsegweight(1,kk,ipts)
            derv2(iy,lxl) = derv2(iy,lxl) - trm*dsegweight(1,ll,jpts)*dsegweight(2,kk,ipts)
            derv2(iz,lxl) = derv2(iz,lxl) - trm*dsegweight(1,ll,jpts)*dsegweight(3,kk,ipts)
            derv2(ix,lyl) = derv2(ix,lyl) - trm*dsegweight(2,ll,jpts)*dsegweight(1,kk,ipts)
            derv2(iy,lyl) = derv2(iy,lyl) - trm*dsegweight(2,ll,jpts)*dsegweight(2,kk,ipts)
            derv2(iz,lyl) = derv2(iz,lyl) - trm*dsegweight(2,ll,jpts)*dsegweight(3,kk,ipts)
            derv2(ix,lzl) = derv2(ix,lzl) - trm*dsegweight(3,ll,jpts)*dsegweight(1,kk,ipts)
            derv2(iy,lzl) = derv2(iy,lzl) - trm*dsegweight(3,ll,jpts)*dsegweight(2,kk,ipts)
            derv2(iz,lzl) = derv2(iz,lzl) - trm*dsegweight(3,ll,jpts)*dsegweight(3,kk,ipts)
          endif
          if (jloc.gt.0) then
            derv2(kx,jxl) = derv2(kx,jxl) - trm*dsegweight(1,ll,jpts)*dsegweight(1,kk,ipts)
            derv2(ky,jxl) = derv2(ky,jxl) - trm*dsegweight(1,ll,jpts)*dsegweight(2,kk,ipts)
            derv2(kz,jxl) = derv2(kz,jxl) - trm*dsegweight(1,ll,jpts)*dsegweight(3,kk,ipts)
            derv2(kx,jyl) = derv2(kx,jyl) - trm*dsegweight(2,ll,jpts)*dsegweight(1,kk,ipts)
            derv2(ky,jyl) = derv2(ky,jyl) - trm*dsegweight(2,ll,jpts)*dsegweight(2,kk,ipts)
            derv2(kz,jyl) = derv2(kz,jyl) - trm*dsegweight(2,ll,jpts)*dsegweight(3,kk,ipts)
            derv2(kx,jzl) = derv2(kx,jzl) - trm*dsegweight(3,ll,jpts)*dsegweight(1,kk,ipts)
            derv2(ky,jzl) = derv2(ky,jzl) - trm*dsegweight(3,ll,jpts)*dsegweight(2,kk,ipts)
            derv2(kz,jzl) = derv2(kz,jzl) - trm*dsegweight(3,ll,jpts)*dsegweight(3,kk,ipts)
          endif
          if (kloc.gt.0) then
            derv2(jx,kxl) = derv2(jx,kxl) - trm*dsegweight(1,ll,jpts)*dsegweight(1,kk,ipts)
            derv2(jy,kxl) = derv2(jy,kxl) - trm*dsegweight(2,ll,jpts)*dsegweight(1,kk,ipts)
            derv2(jz,kxl) = derv2(jz,kxl) - trm*dsegweight(3,ll,jpts)*dsegweight(1,kk,ipts)
            derv2(jx,kyl) = derv2(jx,kyl) - trm*dsegweight(1,ll,jpts)*dsegweight(2,kk,ipts)
            derv2(jy,kyl) = derv2(jy,kyl) - trm*dsegweight(2,ll,jpts)*dsegweight(2,kk,ipts)
            derv2(jz,kyl) = derv2(jz,kyl) - trm*dsegweight(3,ll,jpts)*dsegweight(2,kk,ipts)
            derv2(jx,kzl) = derv2(jx,kzl) - trm*dsegweight(1,ll,jpts)*dsegweight(3,kk,ipts)
            derv2(jy,kzl) = derv2(jy,kzl) - trm*dsegweight(2,ll,jpts)*dsegweight(3,kk,ipts)
            derv2(jz,kzl) = derv2(jz,kzl) - trm*dsegweight(3,ll,jpts)*dsegweight(3,kk,ipts)
          endif
          if (lloc.gt.0) then
            derv2(kx,lxl) = derv2(kx,lxl) + trm*dsegweight(1,ll,jpts)*dsegweight(1,kk,ipts)
            derv2(ky,lxl) = derv2(ky,lxl) + trm*dsegweight(1,ll,jpts)*dsegweight(2,kk,ipts)
            derv2(kz,lxl) = derv2(kz,lxl) + trm*dsegweight(1,ll,jpts)*dsegweight(3,kk,ipts)
            derv2(kx,lyl) = derv2(kx,lyl) + trm*dsegweight(2,ll,jpts)*dsegweight(1,kk,ipts)
            derv2(ky,lyl) = derv2(ky,lyl) + trm*dsegweight(2,ll,jpts)*dsegweight(2,kk,ipts)
            derv2(kz,lyl) = derv2(kz,lyl) + trm*dsegweight(2,ll,jpts)*dsegweight(3,kk,ipts)
            derv2(kx,lzl) = derv2(kx,lzl) + trm*dsegweight(3,ll,jpts)*dsegweight(1,kk,ipts)
            derv2(ky,lzl) = derv2(ky,lzl) + trm*dsegweight(3,ll,jpts)*dsegweight(2,kk,ipts)
            derv2(kz,lzl) = derv2(kz,lzl) + trm*dsegweight(3,ll,jpts)*dsegweight(3,kk,ipts)
          endif
          if (kloc.gt.0) then
            derv2(lx,kxl) = derv2(lx,kxl) + trm*dsegweight(1,ll,jpts)*dsegweight(1,kk,ipts)
            derv2(ly,kxl) = derv2(ly,kxl) + trm*dsegweight(2,ll,jpts)*dsegweight(1,kk,ipts)
            derv2(lz,kxl) = derv2(lz,kxl) + trm*dsegweight(3,ll,jpts)*dsegweight(1,kk,ipts)
            derv2(lx,kyl) = derv2(lx,kyl) + trm*dsegweight(1,ll,jpts)*dsegweight(2,kk,ipts)
            derv2(ly,kyl) = derv2(ly,kyl) + trm*dsegweight(2,ll,jpts)*dsegweight(2,kk,ipts)
            derv2(lz,kyl) = derv2(lz,kyl) + trm*dsegweight(3,ll,jpts)*dsegweight(2,kk,ipts)
            derv2(lx,kzl) = derv2(lx,kzl) + trm*dsegweight(1,ll,jpts)*dsegweight(3,kk,ipts)
            derv2(ly,kzl) = derv2(ly,kzl) + trm*dsegweight(2,ll,jpts)*dsegweight(3,kk,ipts)
            derv2(lz,kzl) = derv2(lz,kzl) + trm*dsegweight(3,ll,jpts)*dsegweight(3,kk,ipts)
          endif
        enddo
      enddo
!
!  Combined loops over weighting atoms - both on the same centre
!
      smtrmi = smtrm*swi
      smtrmj = smtrm*swj
      smtrmpi = smtrmp*swi
      smtrmpj = smtrmp*swj
      if (lcosmicd2) then
        trmi = smtrmi + smtrmpi
        trmj = smtrmj + smtrmpj
      else
        trmi = smtrmi
        trmj = smtrmj
      endif
      do kk = 2,nnearseg(ipts)
        k = nnearsegptr(kk,ipts)
        kx = 3*(k - 1) + 1
        ky = kx + 1
        kz = kx + 2
        kloc = atom2local(k)
        if (kloc.gt.0) then
          kxl = 3*(kloc - 1) + 1
          kyl = kxl + 1
          kzl = kxl + 2
        endif
        do ll = 1,kk-1
          l = nnearsegptr(ll,ipts)
          lx = 3*(l - 1) + 1
          ly = lx + 1
          lz = lx + 2
          lloc = atom2local(l)
          if (lloc.gt.0) then
            lxl = 3*(lloc - 1) + 1
            lyl = lxl + 1
            lzl = lxl + 2
          endif
          if (iloc.gt.0) then
            derv2(kx,ixl) = derv2(kx,ixl) - trmj*d2segweight(1,1,ll,kk,ipts)
            derv2(ky,ixl) = derv2(ky,ixl) - trmj*d2segweight(1,2,ll,kk,ipts)
            derv2(kz,ixl) = derv2(kz,ixl) - trmj*d2segweight(1,3,ll,kk,ipts)
            derv2(kx,iyl) = derv2(kx,iyl) - trmj*d2segweight(2,1,ll,kk,ipts)
            derv2(ky,iyl) = derv2(ky,iyl) - trmj*d2segweight(2,2,ll,kk,ipts)
            derv2(kz,iyl) = derv2(kz,iyl) - trmj*d2segweight(2,3,ll,kk,ipts)
            derv2(kx,izl) = derv2(kx,izl) - trmj*d2segweight(3,1,ll,kk,ipts)
            derv2(ky,izl) = derv2(ky,izl) - trmj*d2segweight(3,2,ll,kk,ipts)
            derv2(kz,izl) = derv2(kz,izl) - trmj*d2segweight(3,3,ll,kk,ipts)
          endif
          if (kloc.gt.0) then
            derv2(ix,kxl) = derv2(ix,kxl) - trmj*d2segweight(1,1,ll,kk,ipts)
            derv2(iy,kxl) = derv2(iy,kxl) - trmj*d2segweight(2,1,ll,kk,ipts)
            derv2(iz,kxl) = derv2(iz,kxl) - trmj*d2segweight(3,1,ll,kk,ipts)
            derv2(ix,kyl) = derv2(ix,kyl) - trmj*d2segweight(1,2,ll,kk,ipts)
            derv2(iy,kyl) = derv2(iy,kyl) - trmj*d2segweight(2,2,ll,kk,ipts)
            derv2(iz,kyl) = derv2(iz,kyl) - trmj*d2segweight(3,2,ll,kk,ipts)
            derv2(ix,kzl) = derv2(ix,kzl) - trmj*d2segweight(1,3,ll,kk,ipts)
            derv2(iy,kzl) = derv2(iy,kzl) - trmj*d2segweight(2,3,ll,kk,ipts)
            derv2(iz,kzl) = derv2(iz,kzl) - trmj*d2segweight(3,3,ll,kk,ipts)
          endif
          if (iloc.gt.0) then
            derv2(lx,ixl) = derv2(lx,ixl) - trmj*d2segweight(1,1,ll,kk,ipts)
            derv2(ly,ixl) = derv2(ly,ixl) - trmj*d2segweight(2,1,ll,kk,ipts)
            derv2(lz,ixl) = derv2(lz,ixl) - trmj*d2segweight(3,1,ll,kk,ipts)
            derv2(lx,iyl) = derv2(lx,iyl) - trmj*d2segweight(1,2,ll,kk,ipts)
            derv2(ly,iyl) = derv2(ly,iyl) - trmj*d2segweight(2,2,ll,kk,ipts)
            derv2(lz,iyl) = derv2(lz,iyl) - trmj*d2segweight(3,2,ll,kk,ipts)
            derv2(lx,izl) = derv2(lx,izl) - trmj*d2segweight(1,3,ll,kk,ipts)
            derv2(ly,izl) = derv2(ly,izl) - trmj*d2segweight(2,3,ll,kk,ipts)
            derv2(lz,izl) = derv2(lz,izl) - trmj*d2segweight(3,3,ll,kk,ipts)
          endif
          if (lloc.gt.0) then
            derv2(ix,lxl) = derv2(ix,lxl) - trmj*d2segweight(1,1,ll,kk,ipts)
            derv2(iy,lxl) = derv2(iy,lxl) - trmj*d2segweight(1,2,ll,kk,ipts)
            derv2(iz,lxl) = derv2(iz,lxl) - trmj*d2segweight(1,3,ll,kk,ipts)
            derv2(ix,lyl) = derv2(ix,lyl) - trmj*d2segweight(2,1,ll,kk,ipts)
            derv2(iy,lyl) = derv2(iy,lyl) - trmj*d2segweight(2,2,ll,kk,ipts)
            derv2(iz,lyl) = derv2(iz,lyl) - trmj*d2segweight(2,3,ll,kk,ipts)
            derv2(ix,lzl) = derv2(ix,lzl) - trmj*d2segweight(3,1,ll,kk,ipts)
            derv2(iy,lzl) = derv2(iy,lzl) - trmj*d2segweight(3,2,ll,kk,ipts)
            derv2(iz,lzl) = derv2(iz,lzl) - trmj*d2segweight(3,3,ll,kk,ipts)
          endif
          if (kloc.gt.0) then
            derv2(lx,kxl) = derv2(lx,kxl) + trmj*d2segweight(1,1,ll,kk,ipts)
            derv2(ly,kxl) = derv2(ly,kxl) + trmj*d2segweight(2,1,ll,kk,ipts)
            derv2(lz,kxl) = derv2(lz,kxl) + trmj*d2segweight(3,1,ll,kk,ipts)
            derv2(lx,kyl) = derv2(lx,kyl) + trmj*d2segweight(1,2,ll,kk,ipts)
            derv2(ly,kyl) = derv2(ly,kyl) + trmj*d2segweight(2,2,ll,kk,ipts)
            derv2(lz,kyl) = derv2(lz,kyl) + trmj*d2segweight(3,2,ll,kk,ipts)
            derv2(lx,kzl) = derv2(lx,kzl) + trmj*d2segweight(1,3,ll,kk,ipts)
            derv2(ly,kzl) = derv2(ly,kzl) + trmj*d2segweight(2,3,ll,kk,ipts)
            derv2(lz,kzl) = derv2(lz,kzl) + trmj*d2segweight(3,3,ll,kk,ipts)
          endif
          if (lloc.gt.0) then
            derv2(kx,lxl) = derv2(kx,lxl) + trmj*d2segweight(1,1,ll,kk,ipts)
            derv2(ky,lxl) = derv2(ky,lxl) + trmj*d2segweight(1,2,ll,kk,ipts)
            derv2(kz,lxl) = derv2(kz,lxl) + trmj*d2segweight(1,3,ll,kk,ipts)
            derv2(kx,lyl) = derv2(kx,lyl) + trmj*d2segweight(2,1,ll,kk,ipts)
            derv2(ky,lyl) = derv2(ky,lyl) + trmj*d2segweight(2,2,ll,kk,ipts)
            derv2(kz,lyl) = derv2(kz,lyl) + trmj*d2segweight(2,3,ll,kk,ipts)
            derv2(kx,lzl) = derv2(kx,lzl) + trmj*d2segweight(3,1,ll,kk,ipts)
            derv2(ky,lzl) = derv2(ky,lzl) + trmj*d2segweight(3,2,ll,kk,ipts)
            derv2(kz,lzl) = derv2(kz,lzl) + trmj*d2segweight(3,3,ll,kk,ipts)
          endif
        enddo
      enddo
!
      do kk = 2,nnearseg(jpts)
        k = nnearsegptr(kk,jpts)
        kx = 3*(k - 1) + 1
        ky = kx + 1
        kz = kx + 2
        kloc = atom2local(k)
        if (kloc.gt.0) then
          kxl = 3*(kloc - 1) + 1
          kyl = kxl + 1
          kzl = kxl + 2
        endif
        do ll = 1,kk-1
          l = nnearsegptr(ll,jpts)
          lx = 3*(l - 1) + 1
          ly = lx + 1
          lz = lx + 2
          lloc = atom2local(l)
          if (lloc.gt.0) then
            lxl = 3*(lloc - 1) + 1
            lyl = lxl + 1
            lzl = lxl + 2
          endif
          if (jloc.gt.0) then
            derv2(kx,jxl) = derv2(kx,jxl) - trmi*d2segweight(1,1,ll,kk,jpts)
            derv2(ky,jxl) = derv2(ky,jxl) - trmi*d2segweight(1,2,ll,kk,jpts)
            derv2(kz,jxl) = derv2(kz,jxl) - trmi*d2segweight(1,3,ll,kk,jpts)
            derv2(kx,jyl) = derv2(kx,jyl) - trmi*d2segweight(2,1,ll,kk,jpts)
            derv2(ky,jyl) = derv2(ky,jyl) - trmi*d2segweight(2,2,ll,kk,jpts)
            derv2(kz,jyl) = derv2(kz,jyl) - trmi*d2segweight(2,3,ll,kk,jpts)
            derv2(kx,jzl) = derv2(kx,jzl) - trmi*d2segweight(3,1,ll,kk,jpts)
            derv2(ky,jzl) = derv2(ky,jzl) - trmi*d2segweight(3,2,ll,kk,jpts)
            derv2(kz,jzl) = derv2(kz,jzl) - trmi*d2segweight(3,3,ll,kk,jpts)
          endif
          if (kloc.gt.0) then
            derv2(jx,kxl) = derv2(jx,kxl) - trmi*d2segweight(1,1,ll,kk,jpts)
            derv2(jy,kxl) = derv2(jy,kxl) - trmi*d2segweight(2,1,ll,kk,jpts)
            derv2(jz,kxl) = derv2(jz,kxl) - trmi*d2segweight(3,1,ll,kk,jpts)
            derv2(jx,kyl) = derv2(jx,kyl) - trmi*d2segweight(1,2,ll,kk,jpts)
            derv2(jy,kyl) = derv2(jy,kyl) - trmi*d2segweight(2,2,ll,kk,jpts)
            derv2(jz,kyl) = derv2(jz,kyl) - trmi*d2segweight(3,2,ll,kk,jpts)
            derv2(jx,kzl) = derv2(jx,kzl) - trmi*d2segweight(1,3,ll,kk,jpts)
            derv2(jy,kzl) = derv2(jy,kzl) - trmi*d2segweight(2,3,ll,kk,jpts)
            derv2(jz,kzl) = derv2(jz,kzl) - trmi*d2segweight(3,3,ll,kk,jpts)
          endif
          if (jloc.gt.0) then
            derv2(lx,jxl) = derv2(lx,jxl) - trmi*d2segweight(1,1,ll,kk,jpts)
            derv2(ly,jxl) = derv2(ly,jxl) - trmi*d2segweight(2,1,ll,kk,jpts)
            derv2(lz,jxl) = derv2(lz,jxl) - trmi*d2segweight(3,1,ll,kk,jpts)
            derv2(lx,jyl) = derv2(lx,jyl) - trmi*d2segweight(1,2,ll,kk,jpts)
            derv2(ly,jyl) = derv2(ly,jyl) - trmi*d2segweight(2,2,ll,kk,jpts)
            derv2(lz,jyl) = derv2(lz,jyl) - trmi*d2segweight(3,2,ll,kk,jpts)
            derv2(lx,jzl) = derv2(lx,jzl) - trmi*d2segweight(1,3,ll,kk,jpts)
            derv2(ly,jzl) = derv2(ly,jzl) - trmi*d2segweight(2,3,ll,kk,jpts)
            derv2(lz,jzl) = derv2(lz,jzl) - trmi*d2segweight(3,3,ll,kk,jpts)
          endif
          if (lloc.gt.0) then
            derv2(jx,lxl) = derv2(jx,lxl) - trmi*d2segweight(1,1,ll,kk,jpts)
            derv2(jy,lxl) = derv2(jy,lxl) - trmi*d2segweight(1,2,ll,kk,jpts)
            derv2(jz,lxl) = derv2(jz,lxl) - trmi*d2segweight(1,3,ll,kk,jpts)
            derv2(jx,lyl) = derv2(jx,lyl) - trmi*d2segweight(2,1,ll,kk,jpts)
            derv2(jy,lyl) = derv2(jy,lyl) - trmi*d2segweight(2,2,ll,kk,jpts)
            derv2(jz,lyl) = derv2(jz,lyl) - trmi*d2segweight(2,3,ll,kk,jpts)
            derv2(jx,lzl) = derv2(jx,lzl) - trmi*d2segweight(3,1,ll,kk,jpts)
            derv2(jy,lzl) = derv2(jy,lzl) - trmi*d2segweight(3,2,ll,kk,jpts)
            derv2(jz,lzl) = derv2(jz,lzl) - trmi*d2segweight(3,3,ll,kk,jpts)
          endif
          if (kloc.gt.0) then
            derv2(lx,kxl) = derv2(lx,kxl) + trmi*d2segweight(1,1,ll,kk,jpts)
            derv2(ly,kxl) = derv2(ly,kxl) + trmi*d2segweight(2,1,ll,kk,jpts)
            derv2(lz,kxl) = derv2(lz,kxl) + trmi*d2segweight(3,1,ll,kk,jpts)
            derv2(lx,kyl) = derv2(lx,kyl) + trmi*d2segweight(1,2,ll,kk,jpts)
            derv2(ly,kyl) = derv2(ly,kyl) + trmi*d2segweight(2,2,ll,kk,jpts)
            derv2(lz,kyl) = derv2(lz,kyl) + trmi*d2segweight(3,2,ll,kk,jpts)
            derv2(lx,kzl) = derv2(lx,kzl) + trmi*d2segweight(1,3,ll,kk,jpts)
            derv2(ly,kzl) = derv2(ly,kzl) + trmi*d2segweight(2,3,ll,kk,jpts)
            derv2(lz,kzl) = derv2(lz,kzl) + trmi*d2segweight(3,3,ll,kk,jpts)
          endif
          if (lloc.gt.0) then
            derv2(kx,lxl) = derv2(kx,lxl) + trmi*d2segweight(1,1,ll,kk,jpts)
            derv2(ky,lxl) = derv2(ky,lxl) + trmi*d2segweight(1,2,ll,kk,jpts)
            derv2(kz,lxl) = derv2(kz,lxl) + trmi*d2segweight(1,3,ll,kk,jpts)
            derv2(kx,lyl) = derv2(kx,lyl) + trmi*d2segweight(2,1,ll,kk,jpts)
            derv2(ky,lyl) = derv2(ky,lyl) + trmi*d2segweight(2,2,ll,kk,jpts)
            derv2(kz,lyl) = derv2(kz,lyl) + trmi*d2segweight(2,3,ll,kk,jpts)
            derv2(kx,lzl) = derv2(kx,lzl) + trmi*d2segweight(3,1,ll,kk,jpts)
            derv2(ky,lzl) = derv2(ky,lzl) + trmi*d2segweight(3,2,ll,kk,jpts)
            derv2(kz,lzl) = derv2(kz,lzl) + trmi*d2segweight(3,3,ll,kk,jpts)
          endif
        enddo
      enddo
    endif
  endif
#ifdef TRACE
  call trace_out('cosmoamatdaddd')
#endif
!
  return
  end
