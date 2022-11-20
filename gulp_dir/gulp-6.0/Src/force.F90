  subroutine force(eforce,lgrad1,lgrad2)
!
!  Subroutine for calculating the correction to the energy
!  due to an external force.
!
!   8/02 Created
!  10/02 Weighting by neqv added
!   2/04 lgrad1 flag added
!   2/04 Force delay added for MD
!   6/04 End of force time added for MD
!   7/04 Sign of the force corrected for the TD part
!   8/04 Bug in parallel use corrected - forces only applied
!        once over all nodes.
!   8/04 Trap added for zero atom case
!   8/04 Addition of external force to force moved here from 
!        initdervs.f
!  11/04 Pi accessed from module
!   6/09 Site energy added
!  11/09 Region derivatives added
!  11/11 Region-region energy contributions stored
!   4/12 xvir, yvir and zvir removed
!   2/18 Trace added
!  12/18 Shear force added
!   1/19 Site energies corrected for shear force
!   2/19 Setting of shear force origin moved to force
!   2/19 Freezing added
!   7/20 Rigid molecule modifications added
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
!  Julian Gale, CIC, Curtin University, July 2020
!
  use configurations, only : forcecfg, tdforcecfg, ltdforcecfg, nregionno
  use configurations, only : shearforcecfg, lshearforcecfg
  use configurations, only : shearforcedircfg, shearforcenormcfg
  use g_constants,    only : pi
  use control,        only : lrigid
  use current
  use derivatives,    only : xdrv, ydrv, zdrv, xregdrv, yregdrv, zregdrv, rstrd
  use derivatives,    only : derv3, sderv2
  use energies,       only : siteenergy, eregion2region
  use general,        only : timesofar
  use m_strain,       only : cartstrterm
  use mdlogic,        only : lmd
  use moldyn,         only : tmdforcestart,tmdforcestop
  use molecule
  use optimisation,   only : nfixatom, lopf, lfreeze
  use parallel,       only : nprocs, procid
  use symmetry,       only : lstr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  logical,  intent(in)  :: lgrad1
  logical,  intent(in)  :: lgrad2
  real(dp), intent(out) :: eforce
!
!  Local variables
!
  integer(i4)           :: i
  integer(i4)           :: is
  integer(i4)           :: ix
  integer(i4)           :: iy
  integer(i4)           :: iz
  integer(i4)           :: js
  integer(i4)           :: nmi
  integer(i4)           :: nregioni
  integer(i4)           :: nsi
  integer(i4)           :: nsj
  logical               :: lopi
  real(dp)              :: distn
  real(dp)              :: dx
  real(dp)              :: dy
  real(dp)              :: dz
  real(dp)              :: drxyzds(6,3)
  real(dp)              :: d2rxyzdsdx(6,3,3)
  real(dp)              :: d2rxyzds2(6,6,3)
  real(dp)              :: esum
  real(dp)              :: fs(6)
  real(dp)              :: fxi
  real(dp)              :: fyi
  real(dp)              :: fzi
  real(dp)              :: sforce
  real(dp)              :: sumf(3)
  real(dp)              :: tdf
  real(dp)              :: tmp(3)
  real(dp)              :: tsf
  real(dp)              :: twopi
  real(dp)              :: vx
  real(dp)              :: vy
  real(dp)              :: vz
  real(dp)              :: xcom
  real(dp)              :: ycom
  real(dp)              :: zcom
  real(dp)              :: xi
  real(dp)              :: yi
  real(dp)              :: zi
#ifdef TRACE
  call trace_in('force')
#endif
!
!  Initialise integral of force x distance
!
  eforce = 0.0_dp
!
!  Force delay for MD
!
  if (lmd) then
    if (timesofar.lt.tmdforcestart(ncf)) then
#ifdef TRACE
      call trace_out('force')
#endif
      return
    endif
    if (tmdforcestop(ncf).gt.0.0_dp) then
      if (timesofar.gt.tmdforcestop(ncf)) then
#ifdef TRACE
        call trace_out('force')
#endif
        return
      endif
    endif
  endif
!
!  If number of atoms is zero then return as there is nothing to do
!
  if (nasym.eq.0) then
#ifdef TRACE
    call trace_out('force')
#endif
    return
  endif
!
!  Loop over asymmetric unit performing integral
!
  do i = 1+procid,nasym,nprocs
    lopi = (.not.lfreeze.or.lopf(i))
    if (lopi) then
      esum = dble(neqv(i))*(forcecfg(1,nsft+i)*(xalat(i) - xinitial(i)) + &
                            forcecfg(2,nsft+i)*(yalat(i) - yinitial(i)) + &
                            forcecfg(3,nsft+i)*(zalat(i) - zinitial(i)))
      eforce = eforce - esum
      nregioni = nregionno(nsft+i)
      eregion2region(nregioni,nregioni) = eregion2region(nregioni,nregioni) - esum
      siteenergy(i) = siteenergy(i) - esum
    endif
  enddo
  if (lgrad1) then
    do i = 1+procid,nasym,nprocs
      lopi = (.not.lfreeze.or.lopf(i))
      if (lopi) then
        xdrv(i) = xdrv(i) - forcecfg(1,nsft+i)
        ydrv(i) = ydrv(i) - forcecfg(2,nsft+i)
        zdrv(i) = zdrv(i) - forcecfg(3,nsft+i)
!
        nregioni = nregionno(nsft+i)
        xregdrv(nregioni) = xregdrv(nregioni) - forcecfg(1,nsft+i)
        yregdrv(nregioni) = yregdrv(nregioni) - forcecfg(2,nsft+i)
        zregdrv(nregioni) = zregdrv(nregioni) - forcecfg(3,nsft+i)
      endif
    enddo
  endif
!
!  Time-dependent forces
!
  if (lmd.and.lgrad1) then
    twopi = 2.0_dp*pi
    tsf = timesofar - tmdforcestart(ncf)
    do i = 1+procid,nasym,nprocs
      lopi = (.not.lfreeze.or.lopf(i))
      if (lopi) then
        nregioni = nregionno(nsft+i)
        if (ltdforcecfg(1,nsft+i)) then
          tdf = - tdforcecfg(1,1,nsft+i)*cos(twopi*(tsf*tdforcecfg(2,1,nsft+i) + tdforcecfg(3,1,nsft+i)))
          xdrv(i) = xdrv(i) + tdf
          xregdrv(nregioni) = xregdrv(nregioni) + tdf
        endif
        if (ltdforcecfg(2,nsft+i)) then
          tdf = - tdforcecfg(1,2,nsft+i)*cos(twopi*(tsf*tdforcecfg(2,2,nsft+i) + tdforcecfg(3,2,nsft+i)))
          ydrv(i) = ydrv(i) + tdf
          yregdrv(nregioni) = yregdrv(nregioni) + tdf
        endif
        if (ltdforcecfg(3,nsft+i)) then
          tdf = - tdforcecfg(1,3,nsft+i)*cos(twopi*(tsf*tdforcecfg(2,3,nsft+i) + tdforcecfg(3,3,nsft+i)))
          zdrv(i) = zdrv(i) + tdf
          zregdrv(nregioni) = zregdrv(nregioni) + tdf
        endif
      endif
    enddo
  endif
!
!  Shear force
!
  if (lshearforcecfg(ncf)) then
!
!  Set up origin of shear force
!
    if (lshearforcecfg(ncf)) then
      if (nfixatom.gt.0.and.nfixatom.lt.numat) then
        shearforceorigin(1) = xclat(nfixatom)
        shearforceorigin(2) = yclat(nfixatom)
        shearforceorigin(3) = zclat(nfixatom)
      else
        shearforceorigin(1) = xclat(1)
        shearforceorigin(2) = yclat(1)
        shearforceorigin(3) = zclat(1)
      endif
    endif
!
    sumf(1:3) = 0.0_dp
!
    ix = - 2
    iy = - 1
    iz =   0
!
    do i = 1+procid,nasym,nprocs
!
!  Find vector to the shear origin for the original coordinates
!
      vx = xinitial(i) - shearforceorigin(1)
      vy = yinitial(i) - shearforceorigin(2)
      vz = zinitial(i) - shearforceorigin(3)
!
!  Compute normal distance from origin shear plane for atom from original structure
!
      distn = vx*shearforcenormcfg(1,ncf) + &
              vy*shearforcenormcfg(2,ncf) + &
              vz*shearforcenormcfg(3,ncf)
!
!  Scale force constant by normal distance
!
      sforce = distn*shearforcecfg(ncf)*dble(neqv(i))
!
      dx = xalat(i) - xinitial(i)
      dy = yalat(i) - yinitial(i)
      dz = zalat(i) - zinitial(i)
!
      esum = (sforce*shearforcedircfg(1,ncf)*dx + &
              sforce*shearforcedircfg(2,ncf)*dy + &
              sforce*shearforcedircfg(3,ncf)*dz)
      eforce = eforce - esum
      siteenergy(i) = siteenergy(i) - esum
!
!  Apply force to atom
!
      fxi = sforce*shearforcedircfg(1,ncf)
      fyi = sforce*shearforcedircfg(2,ncf)
      fzi = sforce*shearforcedircfg(3,ncf)
!
      sumf(1) = sumf(1) - fxi
      sumf(2) = sumf(2) - fyi
      sumf(3) = sumf(3) - fzi
!
      if (lgrad1) then
        xdrv(i) = xdrv(i) - fxi 
        ydrv(i) = ydrv(i) - fyi
        zdrv(i) = zdrv(i) - fzi
        if (lstr) then
!
!  Compute strain derivatives
!
          xi = xalat(i)
          yi = yalat(i)
          zi = zalat(i)
!
!  Molecule handling
!
          if (lrigid) then
            nmi = natmol(nrela2f(i))
!
!  Set COM coordinates
!
            if (nmi.gt.0) then
              xcom = molxyz(1,natinmol(nrela2f(i)),nmi)
              ycom = molxyz(2,natinmol(nrela2f(i)),nmi)
              zcom = molxyz(3,natinmol(nrela2f(i)),nmi)
            endif
          else
            xcom = 0.0_dp
            ycom = 0.0_dp
            zcom = 0.0_dp
          endif
!
          call cartstrterm(ndim,xi,yi,zi,xcom,ycom,zcom,drxyzds,d2rxyzdsdx,d2rxyzds2,lgrad2)
!
          fs(1:nstrains) = 0.0_dp
!
          do is = 1,nstrains
            nsi = nstrptr(is)
            fs(is) = fs(is) + sforce*shearforcedircfg(1,ncf)*drxyzds(nsi,1)
            fs(is) = fs(is) + sforce*shearforcedircfg(2,ncf)*drxyzds(nsi,2)
            fs(is) = fs(is) + sforce*shearforcedircfg(3,ncf)*drxyzds(nsi,3)
          enddo
!
          rstrd(1:nstrains) = rstrd(1:nstrains) - fs(1:nstrains)
!
          if (lgrad2) then
            lopi = (.not.lfreeze.or.lopf(i))
            if (lopi) then
              ix = ix + 3
              iy = iy + 3
              iz = iz + 3
!
              do is = 1,nstrains
                nsi = nstrptr(is)
                derv3(ix,is) = derv3(ix,is) - sforce*shearforcedircfg(1,ncf)*d2rxyzdsdx(nsi,1,1)
                derv3(ix,is) = derv3(ix,is) - sforce*shearforcedircfg(2,ncf)*d2rxyzdsdx(nsi,1,2)
                derv3(ix,is) = derv3(ix,is) - sforce*shearforcedircfg(3,ncf)*d2rxyzdsdx(nsi,1,3)
!
                derv3(iy,is) = derv3(iy,is) - sforce*shearforcedircfg(1,ncf)*d2rxyzdsdx(nsi,2,1)
                derv3(iy,is) = derv3(iy,is) - sforce*shearforcedircfg(2,ncf)*d2rxyzdsdx(nsi,2,2)
                derv3(iy,is) = derv3(iy,is) - sforce*shearforcedircfg(3,ncf)*d2rxyzdsdx(nsi,2,3)
!
                derv3(iz,is) = derv3(iz,is) - sforce*shearforcedircfg(1,ncf)*d2rxyzdsdx(nsi,3,1)
                derv3(iz,is) = derv3(iz,is) - sforce*shearforcedircfg(2,ncf)*d2rxyzdsdx(nsi,3,2)
                derv3(iz,is) = derv3(iz,is) - sforce*shearforcedircfg(3,ncf)*d2rxyzdsdx(nsi,3,3)
              enddo
            endif
!
            do is = 1,nstrains
              nsi = nstrptr(is)
              do js = 1,nstrains
                nsj = nstrptr(js)
                sderv2(js,is) = sderv2(js,is) - sforce*shearforcedircfg(1,ncf)*d2rxyzds2(nsj,nsi,1)
                sderv2(js,is) = sderv2(js,is) - sforce*shearforcedircfg(2,ncf)*d2rxyzds2(nsj,nsi,2)
                sderv2(js,is) = sderv2(js,is) - sforce*shearforcedircfg(3,ncf)*d2rxyzds2(nsj,nsi,3)
              enddo
            enddo
          endif
        endif
      endif
    enddo
!
!  Global sum of energy/forces if needed
!
    if (nprocs.gt.1) then
      call sumall(sumf,tmp,3_i4,"force","sumf")
      sumf(1:3) = tmp(1:3)
    endif
!
    sumf(1:3) = sumf(1:3)/dble(numat)
!
!  Shift forces by sum to restore translational invariance
!
    ix = - 2
    iy = - 1
    iz =   0
!
    do i = 1+procid,nasym,nprocs
      dx = xalat(i) - xinitial(i)
      dy = yalat(i) - yinitial(i)
      dz = zalat(i) - zinitial(i)
!
      esum = (sumf(1)*dx + sumf(2)*dy + sumf(3)*dz)
      eforce = eforce - esum
      siteenergy(i) = siteenergy(i) - esum
!
      if (lgrad1) then
        xdrv(i) = xdrv(i) - sumf(1)
        ydrv(i) = ydrv(i) - sumf(2)
        zdrv(i) = zdrv(i) - sumf(3)
!
        if (lstr) then
!
!  Compute strain derivatives
!
          xi = xalat(i)
          yi = yalat(i)
          zi = zalat(i)
!
!  Molecule handling
!
          if (lrigid) then
            nmi = natmol(nrela2f(i))
!
!  Set COM coordinates
!
            if (nmi.gt.0) then
              xcom = molxyz(1,natinmol(nrela2f(i)),nmi)
              ycom = molxyz(2,natinmol(nrela2f(i)),nmi)
              zcom = molxyz(3,natinmol(nrela2f(i)),nmi)
            endif
          else
            xcom = 0.0_dp
            ycom = 0.0_dp
            zcom = 0.0_dp
          endif
!
          call cartstrterm(ndim,xi,yi,zi,xcom,ycom,zcom,drxyzds,d2rxyzdsdx,d2rxyzds2,lgrad2)
!
          fs(1:nstrains) = 0.0_dp
!
          do is = 1,nstrains
            nsi = nstrptr(is)
            fs(is) = fs(is) + sumf(1)*drxyzds(nsi,1)
            fs(is) = fs(is) + sumf(2)*drxyzds(nsi,2)
            fs(is) = fs(is) + sumf(3)*drxyzds(nsi,3)
          enddo
!
          rstrd(1:nstrains) = rstrd(1:nstrains) - fs(1:nstrains)
!
          if (lgrad2) then
            lopi = (.not.lfreeze.or.lopf(i))
            if (lopi) then
              ix = ix + 3
              iy = iy + 3
              iz = iz + 3
!
              do is = 1,nstrains
                nsi = nstrptr(is)
                derv3(ix,is) = derv3(ix,is) - sumf(1)*d2rxyzdsdx(nsi,1,1)
                derv3(ix,is) = derv3(ix,is) - sumf(2)*d2rxyzdsdx(nsi,1,2)
                derv3(ix,is) = derv3(ix,is) - sumf(3)*d2rxyzdsdx(nsi,1,3)
!
                derv3(iy,is) = derv3(iy,is) - sumf(1)*d2rxyzdsdx(nsi,2,1)
                derv3(iy,is) = derv3(iy,is) - sumf(2)*d2rxyzdsdx(nsi,2,2)
                derv3(iy,is) = derv3(iy,is) - sumf(3)*d2rxyzdsdx(nsi,2,3)
!
                derv3(iz,is) = derv3(iz,is) - sumf(1)*d2rxyzdsdx(nsi,3,1)
                derv3(iz,is) = derv3(iz,is) - sumf(2)*d2rxyzdsdx(nsi,3,2)
                derv3(iz,is) = derv3(iz,is) - sumf(3)*d2rxyzdsdx(nsi,3,3)
              enddo
            endif
!
            do is = 1,nstrains
              nsi = nstrptr(is)
              do js = 1,nstrains
                nsj = nstrptr(js)
                sderv2(js,is) = sderv2(js,is) - sumf(1)*d2rxyzds2(nsj,nsi,1)
                sderv2(js,is) = sderv2(js,is) - sumf(2)*d2rxyzds2(nsj,nsi,2)
                sderv2(js,is) = sderv2(js,is) - sumf(3)*d2rxyzds2(nsj,nsi,3)
              enddo
            enddo
          endif
        endif
      endif
    enddo
  endif
#ifdef TRACE
  call trace_out('force')
#endif
!
  return
  end
