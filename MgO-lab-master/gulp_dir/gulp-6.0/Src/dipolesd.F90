  subroutine dipolesd(edipole,lgrad1,lgrad2)
!
!  Subroutine for dipole correction energy - symmetry adapted version
!
!   7/20 Re-created from dipole3D after rigid molecule modifications
!   8/20 Corrected for the effect of the dielectric constant option
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
!  Copyright Curtin University 2020
!
!  Julian Gale, CIC, Curtin University, August 2020
!
  use g_constants
  use control,        only : lddipole, lrigid
  use current
  use derivatives
  use m_strain,       only : strainddetds, straindet, straind2detds2, cartstrterm
  use mdlogic
  use molecule
  use optimisation
  use parallel
  use symmetry
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  logical,     intent(in)    :: lgrad1
  logical,     intent(in)    :: lgrad2
  real(dp),    intent(out)   :: edipole
!
!  Local variables
!
  integer(i4)                :: i
  integer(i4)                :: ii
  integer(i4)                :: indm
  integer(i4)                :: indmj
  integer(i4)                :: is
  integer(i4)                :: ix
  integer(i4)                :: iy
  integer(i4)                :: iz
  integer(i4)                :: ixf
  integer(i4)                :: iyf
  integer(i4)                :: izf
  integer(i4)                :: ixfo
  integer(i4)                :: iyfo
  integer(i4)                :: izfo
  integer(i4)                :: j
  integer(i4)                :: jj
  integer(i4)                :: js
  integer(i4)                :: jx
  integer(i4)                :: jy
  integer(i4)                :: jz
  integer(i4)                :: jxc
  integer(i4)                :: jyc
  integer(i4)                :: jzc
  integer(i4)                :: kk
  integer(i4)                :: nmi
  integer(i4)                :: nmj
  integer(i4)                :: ns
  integer(i4)                :: nsi
  integer(i4)                :: nsj
  logical                    :: lopi
  logical                    :: lopj
  real(dp)                   :: const
  real(dp)                   :: dip2
  real(dp)                   :: dtrm1
  real(dp)                   :: dtrm2
  real(dp)                   :: dx
  real(dp)                   :: dy
  real(dp)                   :: dz
  real(dp)                   :: drxyzds(6,3)
  real(dp)                   :: drjxyzds(6,3)
  real(dp)                   :: d2rxyzdsdx(6,3,3)
  real(dp)                   :: d2rjxyzdsdx(6,3,3)
  real(dp)                   :: d2rxyzds2(6,6,3)
  real(dp)                   :: d2rjxyzds2(6,6,3)
  real(dp)                   :: q
  real(dp)                   :: qi
  real(dp)                   :: qj
  real(dp)                   :: rstrdloc(6)
  real(dp)                   :: rtmp(6)
  real(dp)                   :: rvol
  real(dp)                   :: sumd3(6,3)
  real(dp)                   :: vol
  real(dp)                   :: volume
  real(dp)                   :: xcom
  real(dp)                   :: ycom
  real(dp)                   :: zcom
  real(dp)                   :: xcomj
  real(dp)                   :: ycomj
  real(dp)                   :: zcomj
  real(dp)                   :: xi
  real(dp)                   :: yi
  real(dp)                   :: zi
  real(dp)                   :: xj
  real(dp)                   :: yj
  real(dp)                   :: zj
#ifdef TRACE
  call trace_in('dipolesd')
#endif
!
!  Check dimensionality
!
  if (ndim.eq.2) then
    call outerror('dipole correction not available for surface',0_i4)
    call stopnow('dipolesd')
  elseif (ndim.eq.1) then
    call outerror('dipole correction not available for polymer',0_i4)
    call stopnow('dipolesd')
  endif
!
  vol = volume(rv)
  rvol = 1.0_dp/vol
  const = 2.0_dp*pi*angstoev*rvol/(2.0_dp*dielectric + 1.0_dp)
!
!  Calculate charge and dipole moment per unit cell
!
  q = 0.0_dp
  dx = 0.0_dp
  dy = 0.0_dp
  dz = 0.0_dp
  if (lddipole) then
!
!  Use dipole relative to initial coordinates
!
    do i = 1,numat
      qi = qf(i)*occuf(i)
      q = q + qi
      dx = dx + qi*(xalat(i) - xinitial(i))
      dy = dy + qi*(yalat(i) - yinitial(i))
      dz = dz + qi*(zalat(i) - zinitial(i))
    enddo
  else
    do i = 1,numat
      qi = qf(i)*occuf(i)
      q = q + qi
      nmi = natmol(i)
      if (lmol.and.nmi.ne.0) then
!
!  If molecules are present that make sure that images
!  belong to the same molecule are used to ensure a
!  consistent dipole moment
!
        indm = nmolind(i)
        call mindtoijk(indm,ii,jj,kk)
        xi = xclat(i) + ii*rv(1,1) + jj*rv(2,1) + kk*rv(3,1)
        yi = yclat(i) + ii*rv(1,2) + jj*rv(2,2) + kk*rv(3,2)
        zi = zclat(i) + ii*rv(1,3) + jj*rv(2,3) + kk*rv(3,3)
        dx = dx + qi*xi
        dy = dy + qi*yi
        dz = dz + qi*zi
      else
        dx = dx + qi*xclat(i)
        dy = dy + qi*yclat(i)
        dz = dz + qi*zclat(i)
      endif
    enddo
  endif
  if (abs(q).gt.1.0d-8) then
    call outwarning('Dipole correction is origin dependent as charge is not zero',0_i4)
  endif
  dip2 = dx*dx + dy*dy + dz*dz
!
!  If dipole is zero then we are finished
!
  if (dip2.eq.0.0_dp) then
    edipole = 0.0_dp
#ifdef TRACE
    call trace_out('dipole3D')
#endif
    return
  endif
!
!  Calculate energy
!
  edipole = const*dip2
!
!  If no derivatives are needed then we are finished
!
  if (.not.lgrad1) then
#ifdef TRACE
    call trace_out('dipole3D')
#endif
    return
  endif
!
!  Internal derivatives
!
  dtrm1 = 2.0_dp*const
  if (lgrad2) then
    ix =  - 2
    iy =  - 1
    iz = 0
    ixf = 1
    iyf = 2
    izf = 3
    ixfo = 1
    iyfo = 2
    izfo = 3
    do i = 1,nasym
      qi = qa(i)*dble(neqv(i))*occua(i)
      xdrv(i) = xdrv(i) + dtrm1*dx*qi
      ydrv(i) = ydrv(i) + dtrm1*dy*qi
      zdrv(i) = zdrv(i) + dtrm1*dz*qi
!
      lopi = (.not.lfreeze.or.lopf(i))
      if (lopi) then
        ix = ix + 3
        iy = iy + 3
        iz = iz + 3
        ixf = ixfo
        iyf = iyfo
        izf = izfo
        ixfo = ixfo + 3*neqv(i)
        iyfo = iyfo + 3*neqv(i)
        izfo = izfo + 3*neqv(i)
        jxc =  - 2
        jyc =  - 1
        jzc = 0
        do j = 1,numat
          lopj = (.not.lfreeze.or.lopf(nrelf2a(j)))
          if (lopj) then
            jxc = jxc + 3
            jyc = jyc + 3
            jzc = jzc + 3
            jx = jxc
            jy = jyc
            jz = jzc
          else
            jx = ixf
            jy = iyf
            jz = izf
          endif
          dtrm2 = dtrm1*qi*qf(j)
          derv2(jx,ix) = derv2(jx,ix) + dtrm2
          derv2(jy,iy) = derv2(jy,iy) + dtrm2
          derv2(jz,iz) = derv2(jz,iz) + dtrm2
        enddo
      endif
    enddo
  else
    do i = procid + 1,nasym,nprocs
      qi = qa(i)*dble(neqv(i))*occua(i)
      xdrv(i) = xdrv(i) + dtrm1*dx*qi
      ydrv(i) = ydrv(i) + dtrm1*dy*qi
      zdrv(i) = zdrv(i) + dtrm1*dz*qi
    enddo
  endif
!
!  Strain derivatives
!
  if (lstr) then
!
!  Initialise local strain array
!
    rstrdloc(1:nstrains) = 0.0_dp
    if (lgrad2) then
      sumd3(1:nstrains,1:3) = 0.0_dp
    endif
!
!  Coordinate strain contribution
!
    ix = - 2
    iy = - 1
    iz =   0
    do i = procid+1,numat,nprocs
      qi = qf(i)*occuf(i)
      if (lddipole) then
        xi = xalat(i) - xinitial(i)
        yi = yalat(i) - yinitial(i)
        zi = zalat(i) - zinitial(i)
        nmi = natmol(i)
        if (nmi.ne.0.and.lrigid) then
          xcom = molxyz(1,natinmol(i),nmi)
          ycom = molxyz(2,natinmol(i),nmi)
          zcom = molxyz(3,natinmol(i),nmi)
        else
          xcom = 0.0_dp
          ycom = 0.0_dp
          zcom = 0.0_dp
        endif
      else
        nmi = natmol(i)
        if (nmi.ne.0) then
!
!  If molecules are present that make sure that images belong to the same molecule are used to 
!  ensure a consistent dipole moment
!
          indm = nmolind(i)
          call mindtoijk(indm,ii,jj,kk)
!
          if (lrigid) then
            xcom = molxyz(1,natinmol(i),nmi)
            ycom = molxyz(2,natinmol(i),nmi)
            zcom = molxyz(3,natinmol(i),nmi)
          else
            xcom = 0.0_dp
            ycom = 0.0_dp
            zcom = 0.0_dp
          endif
!
          xi = xclat(i) + ii*rv(1,1) + jj*rv(2,1) + kk*rv(3,1)
          yi = yclat(i) + ii*rv(1,2) + jj*rv(2,2) + kk*rv(3,2)
          zi = zclat(i) + ii*rv(1,3) + jj*rv(2,3) + kk*rv(3,3)
        else
          xi = xclat(i)
          yi = yclat(i)
          zi = zclat(i)
          xcom = 0.0_dp
          ycom = 0.0_dp
          zcom = 0.0_dp
        endif
      endif
!
      call cartstrterm(ndim,xi,yi,zi,xcom,ycom,zcom,drxyzds,d2rxyzdsdx,d2rxyzds2,lgrad2)
!
      do is = 1,nstrains
        ns = nstrptr(is)
        rstrdloc(is) = rstrdloc(is) + qi*dtrm1*(dx*drxyzds(ns,1) + &
                                                dy*drxyzds(ns,2) + &
                                                dz*drxyzds(ns,3))
      enddo
      if (lgrad2) then
        do is = 1,nstrains
          ns = nstrptr(is)
          sumd3(is,1) = sumd3(is,1) + qi*drxyzds(ns,1)
          sumd3(is,2) = sumd3(is,2) + qi*drxyzds(ns,2)
          sumd3(is,3) = sumd3(is,3) + qi*drxyzds(ns,3)
        enddo
      endif
!
      lopi = (.not.lfreeze.or.lopf(nrelf2a(i)))
      if (lgrad2.and.lopi.and.nrela2f(nrelf2a(i)).eq.i) then
        ix = ix + 3
        iy = iy + 3
        iz = iz + 3
!
!  Mixed internal derivatives
!
        dtrm2 = dtrm1*qi*dble(neqv(nrelf2a(i)))
        if (lfinitestrain) then
          do j = 1,nstrains
            nsj = nstrptr(j)
            derv3(ix,j) = derv3(ix,j) - dtrm2*dx*strainddetds(j)*straindet 
            derv3(iy,j) = derv3(iy,j) - dtrm2*dy*strainddetds(j)*straindet 
            derv3(iz,j) = derv3(iz,j) - dtrm2*dz*strainddetds(j)*straindet 
          enddo
        else
          derv3(ix,1) = derv3(ix,1) - dtrm2*dx
          derv3(iy,1) = derv3(iy,1) - dtrm2*dy
          derv3(iz,1) = derv3(iz,1) - dtrm2*dz
!
          derv3(ix,2) = derv3(ix,2) - dtrm2*dx
          derv3(iy,2) = derv3(iy,2) - dtrm2*dy
          derv3(iz,2) = derv3(iz,2) - dtrm2*dz
!
          derv3(ix,3) = derv3(ix,3) - dtrm2*dx
          derv3(iy,3) = derv3(iy,3) - dtrm2*dy
          derv3(iz,3) = derv3(iz,3) - dtrm2*dz
        endif
        do is = 1,nstrains
          nsi = nstrptr(is)
          derv3(ix,is) = derv3(ix,is) + dtrm2*(dx*d2rxyzdsdx(nsi,1,1) + &
                                               dy*d2rxyzdsdx(nsi,1,2) + &
                                               dz*d2rxyzdsdx(nsi,1,3))
          derv3(iy,is) = derv3(iy,is) + dtrm2*(dx*d2rxyzdsdx(nsi,2,1) + &
                                               dy*d2rxyzdsdx(nsi,2,2) + &
                                               dz*d2rxyzdsdx(nsi,2,3))
          derv3(iz,is) = derv3(iz,is) + dtrm2*(dx*d2rxyzdsdx(nsi,3,1) + &
                                               dy*d2rxyzdsdx(nsi,3,2) + &
                                               dz*d2rxyzdsdx(nsi,3,3))
        enddo
      endif
!
      if (lgrad2) then
        do is = 1,nstrains
          nsi = nstrptr(is)
          do js = 1,nstrains
            nsj = nstrptr(js)
            sderv2(js,is) = sderv2(js,is) + dtrm1*qi*(dx*d2rxyzds2(nsj,nsi,1) + &
                                                      dy*d2rxyzds2(nsj,nsi,2) + &
                                                      dz*d2rxyzds2(nsj,nsi,3))
          enddo
        enddo
!
        do j = 1,numat
          qj = qf(j)*occuf(j)
          if (lddipole) then
            xj = xalat(j) - xinitial(j)
            yj = yalat(j) - yinitial(j)
            zj = zalat(j) - zinitial(j)
            nmj = natmol(j)
            if (nmj.ne.0.and.lrigid) then
              xcom = molxyz(1,natinmol(j),nmj)
              ycom = molxyz(2,natinmol(j),nmj)
              zcom = molxyz(3,natinmol(j),nmj)
            else
              xcom = 0.0_dp
              ycom = 0.0_dp
              zcom = 0.0_dp
            endif
          else
            nmj = natmol(j)
            if (nmj.ne.0) then
!
!  If molecules are present that make sure that images belong to the same molecule are used to ensure a consistent dipole moment
!
              indmj = nmolind(j)
              call mindtoijk(indmj,ii,jj,kk)
!
              if (lrigid) then
                xcom = molxyz(1,natinmol(j),nmj)
                ycom = molxyz(2,natinmol(j),nmj)
                zcom = molxyz(3,natinmol(j),nmj)
              else
                xcom = 0.0_dp
                ycom = 0.0_dp
                zcom = 0.0_dp
              endif
!
              xj = xclat(j) + ii*rv(1,1) + jj*rv(2,1) + kk*rv(3,1)
              yj = yclat(j) + ii*rv(1,2) + jj*rv(2,2) + kk*rv(3,2)
              zj = zclat(j) + ii*rv(1,3) + jj*rv(2,3) + kk*rv(3,3)
            else
              xj = xclat(j)
              yj = yclat(j)
              zj = zclat(j)
              xcom = 0.0_dp
              ycom = 0.0_dp
              zcom = 0.0_dp
            endif
          endif
!
          call cartstrterm(ndim,xj,yj,zj,xcomj,ycomj,zcomj,drjxyzds,d2rjxyzdsdx,d2rjxyzds2,.false.)
!
          do is = 1,nstrains
            nsi = nstrptr(is)
            do js = 1,nstrains
              nsj = nstrptr(js)
              sderv2(js,is) = sderv2(js,is) + dtrm1*qi*qj*(drjxyzds(nsj,1)*drxyzds(nsi,1) + &
                                                           drjxyzds(nsj,2)*drxyzds(nsi,2) + &
                                                           drjxyzds(nsj,3)*drxyzds(nsi,3))
            enddo
          enddo
        enddo
      endif
    enddo
!
    if (lgrad2) then
      ix = - 2
      iy = - 1
      iz =   0
!
!  Cummulative derv3 term for i-j
!
      do i = procid+1,nasym,nprocs
        qi = qa(i)*occua(i)*dble(neqv(i))
        lopi = (.not.lfreeze.or.lopf(i))
        if (lopi) then
          ix = ix + 3
          iy = iy + 3
          iz = iz + 3
          do is = 1,nstrains
            derv3(ix,is) = derv3(ix,is) + dtrm1*qi*sumd3(is,1)
            derv3(iy,is) = derv3(iy,is) + dtrm1*qi*sumd3(is,2)
            derv3(iz,is) = derv3(iz,is) + dtrm1*qi*sumd3(is,3)
          enddo
        endif
      enddo
!
!  Volume corrections to strain second derivatives
!
      if (lfinitestrain) then
        do i = 1,nstrains
          do j = 1,nstrains
            sderv2(j,i) = sderv2(j,i) - rstrdloc(j)*strainddetds(i)*straindet &
                                      - rstrdloc(i)*strainddetds(j)*straindet &
                                      + 2.0_dp*edipole*strainddetds(j)*strainddetds(i)*straindet**2 &
                                      - edipole*straind2detds2(j,i)*straindet
          enddo
        enddo
      else
        do i = 1,3
          do j = 1,3
            sderv2(j,i) = sderv2(j,i) - rstrdloc(j) - rstrdloc(i)
            sderv2(j,i) = sderv2(j,i) + edipole
          enddo
!
          do j = 4,6
            sderv2(j,i) = sderv2(j,i) - rstrdloc(j)
            sderv2(i,j) = sderv2(i,j) - rstrdloc(j)
          enddo
        enddo
      endif
    endif
!
!  Globalise strain first derivatives
!
    if (nprocs.gt.1) then
      call sumall(rstrdloc,rtmp,nstrains,"dipole3D","rstrdloc")
      rstrdloc(1:nstrains) = rtmp(1:nstrains)
    endif
!
!  Volume strain contribution
!
    if (lfinitestrain) then
      do i = 1,nstrains
        rstrdloc(i) = rstrdloc(i) - edipole*strainddetds(i)*straindet
      enddo
    else
      rstrdloc(1) = rstrdloc(1) - edipole
      rstrdloc(2) = rstrdloc(2) - edipole
      rstrdloc(3) = rstrdloc(3) - edipole
    endif
!
!  Add local strain array to main array
!
    strderv(1:nstrains) = strderv(1:nstrains) + rstrdloc(1:nstrains)
  endif
#ifdef TRACE
  call trace_out('dipolesd')
#endif
!
  return
  end
