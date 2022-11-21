  subroutine thermalconductivity_af(mtv,derv2,eigr,Sij,freq,nphonatc,ncfoc,nphonatptr,maxd2,fscale)
!
!  Compute the thermal conductivity in a quasiharmonic supercell approximation
!  according to the method of Allen and Feldman, PRB, 48, 12581 (1993)
!
!   6/12 Created 
!   1/13 JL modified
!   1/13 Loop over Cartesian degrees of freedom added
!   3/13 Correction of constants made
!   6/13 Calculation of thermal conductivity added
!   6/13 Changed to give mode thermal conductivities
!   7/13 Option to have fixed, rather than scaled, broadening
!   8/13 Temperature steps added
!   8/13 Cutoff frequency added for Allen-Feldman contribution
!   8/13 Calculation of propagating contribution added
!   9/13 Modified to allow for separate p-wave velocity
!  12/15 thermalconductivity renamed to thermalconductivity_af
!   1/18 Modified for ghostcell case
!   1/18 Trace added
!   3/19 Multiple temperature ramps added
!   6/20 mcv now changed for mtv for rigid molecules
!   6/20 Rigid molecule modifications added
!   6/20 Check that nearest distances are not ambiguous added
!   1/21 nfreqmin now used as minimum mode number in double loop over phonons
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
!  Copyright Curtin University 2021
!
!  Julian Gale, CIC, Curtin University, January 2021
!
  use configurations, only : nsuperghost
  use g_constants,    only : pi, speedl, avogadro, evtoj, boltz, planck
  use control,        only : lbroaden_scale, lghost, lrigid
  use current
  use general,        only : bfactor, Lor_tol
  use iochannels
  use molecule
  use parallel
  use properties,     only : vs_reuss, vp_reuss
  use thermalcond
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                      :: nphonatptr(*)
  integer(i4), intent(in)                      :: maxd2
  integer(i4), intent(in)                      :: mtv
  integer(i4), intent(in)                      :: ncfoc
  integer(i4), intent(in)                      :: nphonatc
  real(dp),    intent(inout)                   :: derv2(maxd2,*)
  real(dp),    intent(in)                      :: eigr(maxd2,*)
  real(dp),    intent(out)                     :: Sij(mtv,mtv)
  real(dp),    intent(in)                      :: freq(*)
  real(dp),    intent(in)                      :: fscale
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ig
  integer(i4)                                  :: ii
  integer(i4)                                  :: im
  integer(i4)                                  :: ind
  integer(i4)                                  :: ind0
  integer(i4)                                  :: iv
  integer(i4)                                  :: ix
  integer(i4)                                  :: iy
  integer(i4)                                  :: iz
  integer(i4)                                  :: ixyz
  integer(i4)                                  :: j
  integer(i4)                                  :: jg
  integer(i4)                                  :: jj
  integer(i4)                                  :: jm
  integer(i4)                                  :: jx
  integer(i4)                                  :: jy
  integer(i4)                                  :: jz
  integer(i4)                                  :: mindi
  integer(i4)                                  :: mindj
  integer(i4)                                  :: ncfoc2
  integer(i4)                                  :: nfreqmin
  integer(i4)                                  :: nghostcell
  integer(i4)                                  :: nt
  integer(i4)                                  :: nt0
  integer(i4)                                  :: ntr
  integer(i4)                                  :: status
  logical                                      :: lbelow
  logical                                      :: lpropagating
  logical                                      :: lvecok
  logical,        allocatable,            save :: ldone(:)
  real(dp)                                     :: B_pr
  real(dp)                                     :: Bconvert
  real(dp)                                     :: broad
  real(dp)                                     :: constant
  real(dp)                                     :: cmfact
  real(dp)                                     :: cv_i
  real(dp),       allocatable,            save :: Di(:)
  real(dp)                                     :: Di_loc
  real(dp)                                     :: D_pr_max
  real(dp)                                     :: dwavg
  real(dp)                                     :: dwij
  real(dp)                                     :: dxyz
  real(dp)                                     :: expfreq
  real(dp)                                     :: fcut
  real(dp)                                     :: f_pr_max
  real(dp),       allocatable,            save :: freqinv(:)
  real(dp)                                     :: kappa_af
  real(dp)                                     :: kappa_pr
  real(dp)                                     :: kappafct
  real(dp)                                     :: rij
  real(dp)                                     :: rhalf(3)
  real(dp)                                     :: rv2(3,3)
  real(dp)                                     :: tem
  real(dp)                                     :: xd
  real(dp)                                     :: yd
  real(dp)                                     :: zd
  real(dp)                                     :: vol
  real(dp)                                     :: volume
  real(dp)                                     :: vprod
  real(dp),       allocatable,            save :: Vij(:,:)
  real(dp)                                     :: xfreq
  real(dp)                                     :: v_s
  real(dp)                                     :: v_p
#ifdef TRACE
  call trace_in('thermalconductivity')
#endif
!
!  Allocate local array ldone to avoid duplicate multiplies in case of partial occupancy
!
  if (lrigid) then
    ncfoc2 = numatnomol*(numatnomol+1)/2 + numatnomol*nmol + nmol*(nmol+1)/2
  else
    ncfoc2 = ncfoc*(ncfoc+1)/2
  endif
  allocate(freqinv(mtv),stat=status)
  if (status/=0) call outofmemory('thermalconductivity_af','freqinv')
  allocate(ldone(ncfoc2),stat=status)
  if (status/=0) call outofmemory('thermalconductivity_af','ldone')
  allocate(Di(mtv),stat=status)
  if (status/=0) call outofmemory('thermalconductivity_af','Di')
  allocate(Vij(mtv,mtv),stat=status)
  if (status/=0) call outofmemory('thermalconductivity_af','Vij')
!
!  Create inverse frequency factors while trapping translations and imaginary modes
!
  nfreqmin = 0
  do i = 1,mtv
    if (nfreqmin.eq.0.and.freq(i).gt.0.01_dp) nfreqmin = i
    if (freq(i).gt.0.0_dp) then
      freqinv(i) = 1.0_dp/sqrt(freq(i))
    else
      freqinv(i) = 0.0_dp
    endif
  enddo
!
!  Find mean level spacing
!
  dwavg = 0.0_dp
  do i = 1,mtv-1
    if (freq(i).gt.0.0_dp) then
      dwavg = freq(i+1) - freq(i) + dwavg
    elseif (freq(i+1).gt.0.0_dp) then
      dwavg = freq(i+1) + dwavg
    endif
  enddo
  dwavg = dwavg/(mtv - 1)
!
!  Set broadening factor
!
  if (lbroaden_scale) then
    broad = bfactor*dwavg
  else
    broad = bfactor
  endif
!
!  Set constant that includes conversion factor to m^2/s
!
  constant = ((1.0d-17*evtoj*avogadro)**0.5_dp)*(fscale**3)
!
  constant = pi*constant/48.0_dp    ! 1/3 convoluted with 1/2 squared from A7 and 1/2 squared from A8
!
!  Initialise thermal conductivities for each mode
!
  Di(1:mtv) = 0.0_dp
!
!  Find number of ghost cells
!
  if (lghost) then
    nghostcell = nsuperghost(1,ncf)*nsuperghost(2,ncf)*nsuperghost(3,ncf)
  else
    nghostcell = 1
  endif
!
!  Create half cell vectors
!
  do ix = 1,3
    rhalf(ix) = 0.0_dp
    do jx = 1,3
      rhalf(ix) = rhalf(ix) + rv(jx,ix)**2
    enddo
    rhalf(ix) = sqrt(rhalf(ix))
    if (rhalf(ix).gt.1.0d-6) then
      do jx = 1,3
        rv2(jx,ix) = rv(jx,ix)/rhalf(ix)
      enddo
    else
      do jx = 1,3
        rv2(jx,ix) = 0.0_dp
      enddo
    endif
    rhalf(ix) = 0.5_dp*rhalf(ix)
  enddo
!
!  Loop over Cartesian directions
!
  do ixyz = 1,3
!
!  Initialise Vij
!
    Vij(1:mtv,1:mtv) = 0.0_dp
!
!  Initialise ldone
!
    ldone(1:ncfoc2) = .false.
!
!  Scale dynamical matrix elements by minimum image nearest distance between sites
!
    if (lrigid) then
!
!  Non-molecule atoms - non-molecule atoms
!
      do i = 1,numatnomol
        ii = numatnomolptr(i)
        ix = 3*i - 2
        iy = ix + 1
        iz = ix + 2
        do j = 1,i-1
          jj = numatnomolptr(j)
          ind = i*(i-1)/2 + j
          if (.not.ldone(ind)) then
            jx = 3*j - 2
            jy = jx + 1
            jz = jx + 2
!
!  Compute initial vector
!
            xd = xclat(jj) - xclat(ii)
            yd = yclat(jj) - yclat(ii)
            zd = zclat(jj) - zclat(ii)
!
!  Find minimum distance between images
!
            call nearestr(ndim,xd,yd,zd,rv,rij)
!
!  Check that vector is not half a cell vector and therefore ambiguous
!
            lvecok = .true.
            iv = 0
            do while (lvecok.and.iv.lt.3)
              iv = iv + 1
              vprod = xd*rv2(1,iv) + yd*rv2(2,iv) + zd*rv2(3,iv)
              if (abs(abs(vprod)-rhalf(iv)).lt.1.0d-2) lvecok = .false.
            enddo
!
            if (lvecok) then
              if (ixyz.eq.1) then
                dxyz = xd
              elseif (ixyz.eq.2) then
                dxyz = yd
              else
                dxyz = zd
              endif
!  
              Vij(jx,ix) = derv2(jx,ix)*dxyz
              Vij(jy,ix) = derv2(jy,ix)*dxyz
              Vij(jz,ix) = derv2(jz,ix)*dxyz
              Vij(jx,iy) = derv2(jx,iy)*dxyz
              Vij(jy,iy) = derv2(jy,iy)*dxyz
              Vij(jz,iy) = derv2(jz,iy)*dxyz
              Vij(jx,iz) = derv2(jx,iz)*dxyz
              Vij(jy,iz) = derv2(jy,iz)*dxyz
              Vij(jz,iz) = derv2(jz,iz)*dxyz
!  
              Vij(ix,jx) = - derv2(ix,jx)*dxyz
              Vij(iy,jx) = - derv2(iy,jx)*dxyz
              Vij(iz,jx) = - derv2(iz,jx)*dxyz
              Vij(ix,jy) = - derv2(ix,jy)*dxyz
              Vij(iy,jy) = - derv2(iy,jy)*dxyz
              Vij(iz,jy) = - derv2(iz,jy)*dxyz
              Vij(ix,jz) = - derv2(ix,jz)*dxyz
              Vij(iy,jz) = - derv2(iy,jz)*dxyz
              Vij(iz,jz) = - derv2(iz,jz)*dxyz
            endif
            ldone(ind) = .true.
          endif
        enddo
!
!  Self term is zero
!
        Vij(ix,ix) = 0.0_dp
        Vij(iy,ix) = 0.0_dp
        Vij(iz,ix) = 0.0_dp
        Vij(ix,iy) = 0.0_dp
        Vij(iy,iy) = 0.0_dp
        Vij(iz,iy) = 0.0_dp
        Vij(ix,iz) = 0.0_dp
        Vij(iy,iz) = 0.0_dp
        Vij(iz,iz) = 0.0_dp
      enddo
!
!  Molecule - molecule
!
      ind0 = numatnomol*(numatnomol+1)/2 + numatnomol*nmol
      do im = 1,nmol
        mindi = 3*(numatnomol + 2*(im-1))
        i = nmollist(nmolptr(im)+1)
        do jm = 1,im-1
          mindj = 3*(numatnomol + 2*(jm-1))
          ind = ind0 + im*(im-1)/2 + jm
          if (.not.ldone(ind)) then
            j = nmollist(nmolptr(jm)+1)
!
!  Compute vector between centres of mass
!
            xd = xclat(j) - xclat(i) - (molxyz(1,1,jm) - molxyz(1,1,im))
            yd = yclat(j) - yclat(i) - (molxyz(2,1,jm) - molxyz(2,1,im))
            zd = zclat(j) - zclat(i) - (molxyz(3,1,jm) - molxyz(3,1,im))
!
!  Find minimum distance between images
!
            call nearestr(ndim,xd,yd,zd,rv,rij)
!
!  Check that vector is not half a cell vector and therefore ambiguous
!
            lvecok = .true.
            iv = 0
            do while (lvecok.and.iv.lt.3)
              iv = iv + 1
              vprod = xd*rv2(1,iv) + yd*rv2(2,iv) + zd*rv2(3,iv)
              if (abs(abs(vprod)-rhalf(iv)).lt.1.0d-2) lvecok = .false.
            enddo
!
            if (lvecok) then
              if (ixyz.eq.1) then
                dxyz = xd
              elseif (ixyz.eq.2) then
                dxyz = yd
              else
                dxyz = zd
              endif
!  
              do ix = 1,6
                do jx = 1,6
                  Vij(mindj+jx,mindi+ix) = derv2(mindj+jx,mindi+ix)*dxyz
                  Vij(mindi+ix,mindj+jx) = - derv2(mindi+ix,mindj+jx)*dxyz
                enddo
              enddo
            endif
            ldone(ind) = .true.
          endif
        enddo
!
!  Self term is zero
!
        do ix = 1,6
          do jx = 1,6
            Vij(mindi+jx,mindi+ix) = 0.0_dp
          enddo
        enddo
      enddo
!
!  Non-molecule atoms - non-molecule atoms
!
      ind0 = numatnomol*(numatnomol+1)/2
      do i = 1,numatnomol
        ii = numatnomolptr(i)
        ix = 3*i - 2
        iy = ix + 1
        iz = ix + 2
        do jm = 1,nmol
          mindj = 3*(numatnomol + 2*(jm-1))
          ind = ind0 + nmol*(i-1) + jm
          if (.not.ldone(ind)) then
            j = nmollist(nmolptr(jm)+1)
!
!  Compute initial vector
!
            xd = xclat(j) - xclat(ii) - molxyz(1,1,jm)
            yd = yclat(j) - yclat(ii) - molxyz(2,1,jm)
            zd = zclat(j) - zclat(ii) - molxyz(3,1,jm)
!
!  Find minimum distance between images
!
            call nearestr(ndim,xd,yd,zd,rv,rij)
!
!  Check that vector is not half a cell vector and therefore ambiguous
!
            lvecok = .true.
            iv = 0
            do while (lvecok.and.iv.lt.3)
              iv = iv + 1
              vprod = xd*rv2(1,iv) + yd*rv2(2,iv) + zd*rv2(3,iv)
              if (abs(abs(vprod)-rhalf(iv)).lt.1.0d-2) lvecok = .false.
            enddo
!
            if (lvecok) then
              if (ixyz.eq.1) then
                dxyz = xd
              elseif (ixyz.eq.2) then
                dxyz = yd
              else
                dxyz = zd
              endif
!  
              do jx = 1,6
                Vij(mindj+jx,ix) = derv2(mindj+jx,ix)*dxyz
                Vij(mindj+jx,iy) = derv2(mindj+jx,iy)*dxyz
                Vij(mindj+jx,iz) = derv2(mindj+jx,iz)*dxyz
                Vij(ix,mindj+jx) = - derv2(ix,mindj+jx)*dxyz
                Vij(iy,mindj+jx) = - derv2(iy,mindj+jx)*dxyz
                Vij(iz,mindj+jx) = - derv2(iz,mindj+jx)*dxyz
              enddo
            endif
            ldone(ind) = .true.
          endif
        enddo
      enddo
    else
      do i = 1,nphonatc,nghostcell
        ii = nphonatptr(i)
        ig = (ii - 1)/nghostcell + 1
        ix = 3*ig - 2
        iy = ix + 1
        iz = ix + 2
        do j = 1,i-1,nghostcell
          jj = nphonatptr(j)
          jg = (jj - 1)/nghostcell + 1
          ind = ig*(ig-1)/2 + jg
          if (.not.ldone(ind)) then
            jx = 3*jg - 2
            jy = jx + 1
            jz = jx + 2
!
!  Compute initial vector
!
            xd = xclat(jj) - xclat(ii)
            yd = yclat(jj) - yclat(ii)
            zd = zclat(jj) - zclat(ii)
!
!  Find minimum distance between images
!
            call nearestr(ndim,xd,yd,zd,rv,rij)
!
!  Check that vector is not half a cell vector and therefore ambiguous
!
            lvecok = .true.
            iv = 0
            do while (lvecok.and.iv.lt.3)
              iv = iv + 1
              vprod = xd*rv2(1,iv) + yd*rv2(2,iv) + zd*rv2(3,iv)
              if (abs(abs(vprod)-rhalf(iv)).lt.1.0d-2) lvecok = .false.
            enddo
!
            if (lvecok) then
              if (ixyz.eq.1) then
                dxyz = xd
              elseif (ixyz.eq.2) then
                dxyz = yd
              else
                dxyz = zd
              endif
!  
              Vij(jx,ix) = derv2(jx,ix)*dxyz
              Vij(jy,ix) = derv2(jy,ix)*dxyz
              Vij(jz,ix) = derv2(jz,ix)*dxyz
              Vij(jx,iy) = derv2(jx,iy)*dxyz
              Vij(jy,iy) = derv2(jy,iy)*dxyz
              Vij(jz,iy) = derv2(jz,iy)*dxyz
              Vij(jx,iz) = derv2(jx,iz)*dxyz
              Vij(jy,iz) = derv2(jy,iz)*dxyz
              Vij(jz,iz) = derv2(jz,iz)*dxyz
!  
              Vij(ix,jx) = - derv2(ix,jx)*dxyz
              Vij(iy,jx) = - derv2(iy,jx)*dxyz
              Vij(iz,jx) = - derv2(iz,jx)*dxyz
              Vij(ix,jy) = - derv2(ix,jy)*dxyz
              Vij(iy,jy) = - derv2(iy,jy)*dxyz
              Vij(iz,jy) = - derv2(iz,jy)*dxyz
              Vij(ix,jz) = - derv2(ix,jz)*dxyz
              Vij(iy,jz) = - derv2(iy,jz)*dxyz
              Vij(iz,jz) = - derv2(iz,jz)*dxyz
            endif
            ldone(ind) = .true.
          endif
        enddo
!
!  Self term is zero
!
        Vij(ix,ix) = 0.0_dp
        Vij(iy,ix) = 0.0_dp
        Vij(iz,ix) = 0.0_dp
        Vij(ix,iy) = 0.0_dp
        Vij(iy,iy) = 0.0_dp
        Vij(iz,iy) = 0.0_dp
        Vij(ix,iz) = 0.0_dp
        Vij(iy,iz) = 0.0_dp
        Vij(iz,iz) = 0.0_dp
      enddo
    endif
!
!  Multiply eigenvectors by distance weighted dynamical matrix from both sides
!
    call dgemm('N','N',mtv,mtv,mtv,1.0_dp,Vij,mtv,eigr,maxd2,0.0_dp,Sij,mtv)
    call dgemm('T','N',mtv,mtv,mtv,1.0_dp,eigr,maxd2,Sij,mtv,0.0_dp,Vij,mtv)
!
!  Scale by constants and frequency factors to get to Sij
!
    do i = 1,mtv
      do j = 1,mtv
        Sij(j,i) = Vij(j,i)*freqinv(i)*freqinv(j)*(freq(i) + freq(j))
      enddo
    enddo
!
!  Compute Di values (factors of pi have been cancelled)
!
    do i = nfreqmin,mtv
      Di_loc = 0.0_dp
!
!  Sum over coupling with mode j weighted by Lorentzian factor
!
      do j = nfreqmin,mtv
        dwij = (1.0/pi)*broad/( (freq(j) - freq(i))**2 + broad**2 )
        if (dwij.gt.Lor_tol) then
          Di_loc = Di_loc + dwij*Sij(j,i)**2
        endif
      enddo
!
!  Scale by constants and inverse frequency squared
!
      Di(i) = Di(i) + Di_loc*constant/(freq(i)**2)
    enddo
!
!  End loop over Cartesian degrees of freedom
!
  enddo
!
!  Set cutoff for Allen-Feldman thermal conductivity
!
  fcut = omega_af(ncf)
!
!  Set wave velocities
!
  if (v_s_cfg(ncf).gt.0.0_dp) then
    v_s = v_s_cfg(ncf)
  else
    v_s = vs_reuss
  endif
  if (v_p_cfg(ncf).gt.0.0_dp) then
    v_p = v_p_cfg(ncf)
  else
    v_p = vp_reuss
  endif
!
!  Do we need to do propagating contribution?
!
  lpropagating = .false.
  if (fcut.gt.0.0_dp) then
    lpropagating = .true.
    B_pr = B_pr_cfg(ncf)
    if (B_pr.eq.0.0_dp) then
!
!  If fcut > 0 and B hasn't been input then find maximum mode diffusivity below the cut off in order to estimate B
!
      D_pr_max = 0.0_dp
      f_pr_max = 0.0_dp
      lbelow = .true.
      i = 0
      do while (lbelow.and.i.lt.mtv)
        i = i + 1
        lbelow = (freq(i).lt.fcut)
        if (lbelow) then
          if (Di(i).gt.D_pr_max) then
            D_pr_max = Di(i)
            f_pr_max = freq(i)
          endif
        endif
      enddo
      B_pr = 3.0_dp*D_pr_max*(f_pr_max**2)/(100.0_dp*v_s**2)
!
!  If D_pr_max = 0 then we can't do propagating contribution
!
      if (f_pr_max.eq.0.0_dp) lpropagating = .false.
    endif
  endif  
  if (lpropagating) then
!
!  If propagating contribution is to be computed then find B value in s/m^2
!
    Bconvert = (2.0_dp*pi*speedl*1.0d-9)**2
!
!  Now compute thermal conductivity contribution from propagation
!
    kappa_pr = 4.0_dp*pi*boltz*B_pr*fcut*(speedl**3)*1.0d-7*(2.0_dp/v_s + 1.0_dp/v_p)/3.0_dp
  endif
!*********************************************************************
!  Calculation of thermal conductivity as a function of temperature  *
!*********************************************************************
  if (temperature.lt.1.0d-6.and.ntemperatureramp.eq.0) then
!
!  Check that temperature is not zero before computing the thermal conductivity
!
    call outerror('thermal conductivity cannot be calculated as temperature is too low',0_i4)
    call stopnow('thermalconductivity_af')
  endif
  if (ioproc) then
    write(ioout,'(/,''  Thermal conductivity: '',/)')
    write(ioout,'(''  Lorentzian broadening factor = '',f14.6,'' cm-1'')') broad
    write(ioout,'(''                               = '',f14.6,'' meV'')') broad/8.0655_dp
    write(ioout,'(''  Frequency cutoff for AF term = '',f14.6,'' cm-1'')') fcut
    if (lpropagating) then
      write(ioout,'(''                               = '',f14.6,'' meV'')') fcut/8.0655_dp
      write(ioout,'(''  Transverse   speed of sound  = '',f14.6,'' m/s'')') 1000.0_dp*v_s
      write(ioout,'(''  Longitudinal speed of sound  = '',f14.6,'' m/s'')') 1000.0_dp*v_p
      write(ioout,'(''  Estimated B for propagation  = '',f14.6,'' s/km**2'')') B_pr*1.0d6
      write(ioout,'(''                               = '',f14.6,'' 10**14 rads**2/s'',/)') B_pr*Bconvert
    else
      write(ioout,'(''                               = '',f14.6,'' meV'',/)') fcut/8.0655_dp
    endif
  endif
!
!  Divide by volume and convert units to SI
!
  vol = volume(rv)
  kappafct = 1.0d30/vol
  if (ntemperatureramp.eq.0) then
    cmfact = planck*speedl/(boltz*temperature)
!
!  Output banner for thermal conductivity
!
    if (ioproc) then
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'('' Mode    : Frequency              Mode diffusivity        Thermal conductivity  '')')
      write(ioout,'(''         :   (cm-1)                  (cm**2/s)                  (W/m.K)         '')')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
    endif
!
!  Compute Di values (factors of pi have been cancelled)
!
    kappa_af = 0.0_dp
    do i = nfreqmin,mtv
      if (freq(i).gt.fcut) then
        xfreq = freq(i)*cmfact
        expfreq = exp(xfreq)
        cv_i = boltz*xfreq*xfreq*expfreq/(expfreq - 1.0_dp)**2
        write(ioout,'(i6,2x,f12.4,6x,f22.10,3x,f22.10)') i,freq(i),Di(i)*1.0d4,cv_i*kappafct*Di(i)
        kappa_af = kappa_af + cv_i*kappafct*Di(i)
      endif
    enddo
!
!  Close output
!
    if (ioproc) then
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''  Thermal conductivity (Allen-Feldman) = '',f12.6,'' W/(m.K) '')') kappa_af
      if (lpropagating) then
        write(ioout,'(''  Thermal conductivity (propagation)   = '',f12.6,'' W/(m.K) '')') kappa_pr
        write(ioout,'(''--------------------------------------------------------------------------------'')')
        write(ioout,'(''  Thermal conductivity (total)         = '',f12.6,'' W/(m.K) '')') kappa_af + kappa_pr
      endif
      write(ioout,'(''--------------------------------------------------------------------------------'',/)')
    endif
  else
    if (ioproc) then
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''  Thermal conductivity :'')')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      if (lpropagating) then
        write(ioout,'(''  Temperature               Allen-Feldman     Propagation         Total    '')')
        write(ioout,'(''      (K)                      (W/m/K)          (W/m/K)          (W/m/K)   '')')
      else
        write(ioout,'(''  Temperature               Allen-Feldman     '')')
        write(ioout,'(''      (K)                      (W/m/K)        '')')
      endif
      write(ioout,'(''--------------------------------------------------------------------------------'')')
    endif
!
!  Loop over temperatures
!
    do ntr = 1,ntemperatureramp
!
!  For first ramp need to start from initial temperature
!  For subsequent ramps this would be the same as the end of the previous ramp
!
      if (ntr.eq.1) then
        nt0 = 0
      else
        nt0 = 1
      endif
      do nt = nt0,ntemperaturestep(ntr)
        tem = temperaturestart(ntr) + dble(nt)*temperaturestep(ntr)
        if (tem.gt.1.0d-6) then
          cmfact = planck*speedl/(boltz*tem)
!
!  Compute Di values (factors of pi have been cancelled)
!
          kappa_af = 0.0_dp
          do i = nfreqmin,mtv
            if (freq(i).gt.fcut) then
              xfreq = freq(i)*cmfact
              expfreq = exp(xfreq)
              cv_i = boltz*xfreq*xfreq*expfreq/(expfreq - 1.0_dp)**2
              kappa_af = kappa_af + cv_i*kappafct*Di(i)
            endif
          enddo
          if (ioproc) then
            if (lpropagating) then
              write(ioout,'(2x,f10.3,11x,3(5x,f12.6))') tem,kappa_af,kappa_pr,kappa_af+kappa_pr
            else
              write(ioout,'(2x,f10.3,16x,f12.6)') tem,kappa_af
            endif
          endif
        endif
      enddo
    enddo
    if (ioproc) then
      write(ioout,'(''--------------------------------------------------------------------------------'',/)')
    endif
  endif
!
!  Deallocate local memory
!
  deallocate(Vij,stat=status)
  if (status/=0) call deallocate_error('thermalconductivity_af','Vij')
  deallocate(Di,stat=status)
  if (status/=0) call deallocate_error('thermalconductivity_af','Di')
  deallocate(ldone,stat=status)
  if (status/=0) call deallocate_error('thermalconductivity_af','ldone')
  deallocate(freqinv,stat=status)
  if (status/=0) call deallocate_error('thermalconductivity_af','freqinv')
#ifdef TRACE
  call trace_out('thermalconductivity')
#endif
!
  return
  end
