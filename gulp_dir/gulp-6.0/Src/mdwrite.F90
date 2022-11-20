  subroutine mdwrite(fc,lprod,ncstep)
!
!  Writes out MD dumpfile on channel 31 plus the normal dumpfile.
!  The MD dumpfile contains information about all sampled 
!  configurations so far, whereas the normal dumpfile contains
!  only information about the current structure and any running
!  averages which are required to allow a restart to be performed.
!
!  lprod  = .true. => in production phase of run
!         = .false.=> in equilibration phase of run
!  ncstep = number of current step
!
!  if (lprod) then only write out ordinary dumpfile so that
!  equilibration can be restarted.
!
!  Format of MD dumpfile is:
!
!  numat
!  
!  then for each configuration:
!
!  current_time, kinetic energy, potential energy, temperature
!  all the absolute x coordinates
!  all the absolute y coordinates
!  all the absolute z coordinates
!  all the velocities in the x direction in Angs/ps
!  all the velocities in the y direction in Angs/ps
!  all the velocities in the z direction in Angs/ps
!  if (conp) then
!    all the cell components
!    all the cell component velocities
!  endif
!
!   2/97 Modifications added for GC/shell model MD (JRH)
!   4/98 Writing out of numat every cycle stopped 
!   5/98 XYZ format added as an option
!   7/98 History file format added as an option
!   8/02 ASCII option added
!   1/04 Option to write during equilibration added
!   4/04 Output of pressure file added
!  10/04 al/bl/cl -> aloc/bloc/cloc
!   6/05 Format of XYZ atom labels corrected
!   5/06 Mass now taken from species values
!   5/07 Call to date_and_time corrected
!   2/08 Missing unit conversion for ascii trajectory velocities added
!   7/08 Writing to arc file moved to subroutine addframe2arc
!   9/11 Forces added to MD trajectory file
!  10/11 Site energies added to MD trajectory file for thermal conductivity
!   3/12 Output of DCD file added
!   8/17 Flushing of buffers added
!  10/17 Output of dipoles added
!   2/18 Trace added
!   4/18 Sequential XSF file output added
!   6/18 Parallel I/O trap added
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
!  Julian Gale, CIC, Curtin University, June 2018
!
  use control
  use g_constants,    only : degtorad
  use current
  use derivatives,    only : strderv, xdrv, ydrv, zdrv
  use dump
  use element
  use energies,       only : siteenergy
  use gulp_files
  use general
  use moldyn
  use molecule,       only : natmol, nmolind, lmol
  use parallel,       only : ioproc
  use shells
  use species,        only : massspec
  use symmetry
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  use velocities
  implicit none
!
!  Passed variables
!
  real(dp)                     :: fc
  logical                      :: lprod
  integer(i4)                  :: ncstep
!
!  Local variables
!
  character(len=4)             :: lab
  character(len=8)             :: lab8
  integer(i4)                  :: i
  integer(i4)                  :: ii
  integer(i4)                  :: indm
  integer(i4)                  :: ix
  integer(i4)                  :: iy
  integer(i4)                  :: iz
  integer(i4)                  :: inat
  integer(i4)                  :: itype
  integer(i4)                  :: j
  integer(i4)                  :: keytrj
  logical                      :: lshel
  real(dp)                     :: ctmp(6)
  real(dp)                     :: dx
  real(dp)                     :: dy
  real(dp)                     :: dz
  real(dp)                     :: pconv
  real(dp)                     :: qi
  real(dp)                     :: rnsteps
  real(dp)                     :: tfct
  real(dp)                     :: volume
!
!  Arrays for DCD file - must be real*4
!
  real*4,         allocatable, save :: xloc(:)
  real*4,         allocatable, save :: yloc(:)
  real*4,         allocatable, save :: zloc(:)
!
!  Ensure that only a single node writes this file
!
  if (.not.ioproc) return
#ifdef TRACE
  call trace_in('mdwrite')
#endif
!
  tfct = 1.0_dp/tstep(ncf)
!
!  Ordinary dumpfile
!
  if (idump.gt.0) call dumpdur(idump,1_i4)
  if (.not.lprod.and..not.ltrjequil) then
#ifdef TRACE
    call trace_out('mdwrite')
#endif
    return
  endif
!
!  MD dumpfile on channel 31
!
  if (ltrj) then
    if (ltrjascii) then
      write(31,'(''#  Time/KE/E/T'')')
      write(31,'(4(g25.15,1x))') timesofar,ekin,fc,temperature
      write(31,'(''#  Coordinates'')')
      do i = 1,numat
        write(31,'(3(g25.15,1x))') xalat(i),yalat(i),zalat(i)
      enddo
      write(31,'(''#  Velocities'')')
      do i = 1,numat
        write(31,'(3(g25.15,1x))') velx(i)*tfct,vely(i)*tfct,velz(i)*tfct
      enddo
      write(31,'(''#  Derivatives '')')
      do i = 1,numat
        write(31,'(3(g25.15,1x))') xdrv(i),ydrv(i),zdrv(i)
      enddo
      write(31,'(''#  Site energies '')')
      do i = 1,numat
        write(31,'(g25.15)') siteenergy(i)
      enddo
      if (lconp) then
        write(31,'(''#  Cell'')')
        do i = 1,3
          write(31,'(3(g25.15,1x))') (xcell(3*(i-1)+j),j=1,3)
        enddo
        write(31,'(''#  Cell strain velocities'')')
        do i = 1,2
          write(31,'(3(g25.15,1x))') (velc(3*(i-1)+j)*tfct,j=1,3)
        enddo
      endif
    else
      write(31) timesofar,ekin,fc,temperature
      write(31) (xalat(i),i=1,numat)
      write(31) (yalat(i),i=1,numat)
      write(31) (zalat(i),i=1,numat)
      write(31) (velx(i)*tfct,i=1,numat)
      write(31) (vely(i)*tfct,i=1,numat)
      write(31) (velz(i)*tfct,i=1,numat)
      write(31) (xdrv(i),i=1,numat)
      write(31) (ydrv(i),i=1,numat)
      write(31) (zdrv(i),i=1,numat)
      write(31) (siteenergy(i),i=1,numat)
      if (lconp) then
        write(31) (xcell(i),i=1,9)
        write(31) (velc(i)*tfct,i=1,nstrains)
      endif
    endif
    call gflush(31_i4)
  endif
  if (larc.and.lmovie) then
    call addframe2arc(32_i4,fc,.true.,ncstep)
    call gflush(32_i4)
  endif
!
!  XYZ file output
!
  if (lxyz.and.lxyzmovie) then
    write(33,'(i8)') numat - nshell
    write(33,'(''SCF Done '',f24.8)') fc
    do i = 1,numat
      inat = nat(i)
      if (inat.le.maxele) then
        call label(inat,0_i4,lab)
        write(33,'(a4,1x,3f20.9)') lab,xalat(i),yalat(i),zalat(i)
      endif
    enddo
    call gflush(33_i4)
  endif
!
!  History file output
!
  if (lhis) then
    keytrj = 1
    write(34,'(''timestep'',4i10,f12.6)') ncstep,numat,keytrj,ndim,tstep(ncf)
    if (ndim.eq.3) then
      do i = 1,3
        write(34,'(3g12.4)')(rv(j,i),j=1,3)
      enddo
    endif
    do i = 1,numat
      inat = nat(i)
      lshel = .false.
      if (inat.gt.maxele) then
        inat = inat - maxele
        lshel = .true.
      endif
      itype = natype(i)
      call label(inat,itype,lab)
      lab8 = lab
      if (lshel) then
        ii = index(lab8,' ')
        lab8(ii:ii+1) = '_s'
      endif
      write(34,'(a8,i10,2f12.6)') lab8,i,massspec(nspecptr(i)),qa(i)
      write(34,'(1p,3e12.4)') xalat(i),yalat(i),zalat(i)
      write(34,'(1p,3e12.4)') velx(i)*tfct,vely(i)*tfct,velz(i)*tfct
    enddo
    call gflush(34_i4)
  endif
!
!  Pressure file output
!
  if (lpre.and.ndim.eq.3) then
    pconv = 3.0_dp*pfct/volume(rv)
    rnsteps = pconv/dble(naverpt)
    ctmp(1:6) = - strderv(1:6)
    call mdkestrain(ctmp)
    write(35,'(''Instant: '',i10,6f20.8)') naverpt,(pconv*ctmp(j),j=1,6)
    write(35,'(''Average: '',i10,6f20.8)') naverpt,(rnsteps*sumstress(j),j=1,6)
    call gflush(35_i4)
  endif
!
!  DCD file output
!
  if (ldcd) then
    if (ndim.eq.3) then
      write(36) a,cos(gamma*degtorad),b,cos(beta*degtorad),cos(alpha*degtorad),c
    endif
!
    allocate(xloc(numat))
    allocate(yloc(numat))
    allocate(zloc(numat))
!
    do i = 1,numat
      xloc(i) = real(xalat(i),4)
      yloc(i) = real(yalat(i),4)
      zloc(i) = real(zalat(i),4)
    enddo
!
    write(36) xloc
    write(36) yloc
    write(36) zloc
!
    deallocate(zloc)
    deallocate(yloc)
    deallocate(xloc)
!
    call gflush(36_i4)
  endif
!
!  Dipole file output
!
  if (ldip) then
    dx = 0.0_dp
    dy = 0.0_dp
    dz = 0.0_dp
    if (lddipole) then
!
!  Use dipole relative to initial coordinates
!
      do i = 1,numat
        qi = qf(i)*occuf(i)
        dx = dx + qi*(xalat(i) - xinitial(i))
        dy = dy + qi*(yalat(i) - yinitial(i))
        dz = dz + qi*(zalat(i) - zinitial(i))
      enddo
    else
      do i = 1,numat
        qi = qf(i)*occuf(i)
        if (lmol.and.natmol(i).ne.0) then
!
!  If molecules are present that make sure that images belong to the same molecule 
!  are used to ensure a consistent dipole moment
!
          indm = nmolind(i)
          call mindtoijk(indm,ix,iy,iz)
          dx = dx + qi*(xclat(i) + ix*rv(1,1) + iy*rv(2,1) + iz*rv(3,1))
          dy = dy + qi*(yclat(i) + ix*rv(1,2) + iy*rv(2,2) + iz*rv(3,2))
          dz = dz + qi*(zclat(i) + ix*rv(1,3) + iy*rv(2,3) + iz*rv(3,3))
        else
          dx = dx + qi*xclat(i)
          dy = dy + qi*yclat(i)
          dz = dz + qi*zclat(i)
        endif
      enddo
    endif
    write(37,'(f12.4,3f12.6)') timesofar,4.8_dp*dx,4.8_dp*dy,4.8_dp*dz
    call gflush(37_i4)
  endif
!
!  XSF sequential file output
!
  if (lxsf.and.lxsfseq) then
    nxsfseq = nxsfseq + 1
    call outxsfseq(38_i4,nxsfseq,.true.,fc)
  endif
#ifdef TRACE
  call trace_out('mdwrite')
#endif
!
  return
  end
