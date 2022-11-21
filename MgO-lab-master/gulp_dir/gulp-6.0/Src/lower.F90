  subroutine lower(mtv,mtvrptr,freq,nphonatc,nphonatptr,ncfoc,iocptr,maxd2,eigr)
!
!  Lowers the symmetry of a system according to the imaginary eigenvalue eigenvectors.
!
!   3/02 Created from peigen/peigeng
!   7/02 Modified to allow for region 1 phonons only
!   7/02 Probable bug in assignment of shift atoms in xcfg corrected -
!        displacements were being applied to i not j
!   4/04 Lowering of shell position added
!   3/07 Gauss renamed to GULP_gauss
!   8/15 If lower has been applied then set flag to indicate this
!  12/16 mcvrptr added to arguments for the parallel case
!   1/18 Modified to allow for the present of ghost cells
!   2/18 Trace added
!   3/18 Tolerance on imaginary frequencies now can be set from input files
!  11/18 Try to ensure consistent choice of direction
!  11/18 Option added to switch direction of lowering
!   5/20 Modified for rigid molecules
!   5/20 mcv changed to mtv
!   6/20 nmolcore changes added
!   8/20 nghostcell removed as it is not used
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
  use configurations, only : xcfg, ycfg, zcfg, llowered
  use control
  use current
  use general,        only : lowerscale, frqtol, lowersign
  use iochannels
  use molecule
  use parallel
  use shells,         only : ncsptr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4),      intent(in)   :: iocptr(*)
  integer(i4),      intent(in)   :: maxd2
  integer(i4),      intent(in)   :: mtv            ! Total number of modes = mcv + mmv
  integer(i4),      intent(in)   :: mtvrptr(mtv)
  integer(i4),      intent(in)   :: ncfoc
  integer(i4),      intent(in)   :: nphonatc
  integer(i4),      intent(in)   :: nphonatptr(*)
  real(dp),         intent(in)   :: freq(*)
  real(dp),         intent(in)   :: eigr(maxd2,*)
!
!  Local variables
!
  integer(i4)                    :: i
  integer(i4)                    :: ind
  integer(i4)                    :: ixlarge
  integer(i4)                    :: iylarge
  integer(i4)                    :: izlarge
  integer(i4)                    :: j
  integer(i4)                    :: jj
  integer(i4)                    :: jjs
  integer(i4)                    :: k
  integer(i4)                    :: l
  integer(i4)                    :: n
  integer(i4)                    :: nm
  integer(i4)                    :: nimag
  integer(i4)                    :: nmode
  integer(i4)                    :: status
  real(dp)                       :: crd(3)
  real(dp)                       :: rmat(3,4)
  real(dp)                       :: rotm(3,3)
  real(dp)                       :: scale
  real(dp),    allocatable, save :: sum1(:)
  real(dp),    allocatable, save :: sum2(:)
  real(dp)                       :: theta
  real(dp)                       :: tmp(3)
  real(dp),    allocatable, save :: xshiftT(:)   ! Shifts for translation
  real(dp),    allocatable, save :: yshiftT(:)
  real(dp),    allocatable, save :: zshiftT(:)
  real(dp),    allocatable, save :: xshiftR(:)   ! Shifts for rotation
  real(dp),    allocatable, save :: yshiftR(:)
  real(dp),    allocatable, save :: zshiftR(:)
#ifdef TRACE
  call trace_in('lower')
#endif
!*****************************************************
!  Lowering of symmetry according to imaginary modes  *
!******************************************************
  if (index(keyword,'lowe').ne.0) then
!
!  Allocate arrays to hold coordinate shifts
!
    allocate(xshiftT(numat),stat=status)
    if (status/=0) call outofmemory('lower','xshiftT')
    allocate(yshiftT(numat),stat=status)
    if (status/=0) call outofmemory('lower','yshiftT')
    allocate(zshiftT(numat),stat=status)
    if (status/=0) call outofmemory('lower','zshiftT')
!
    if (lrigid) then
      allocate(xshiftR(numat),stat=status)
      if (status/=0) call outofmemory('lower','xshiftR')
      allocate(yshiftR(numat),stat=status)
      if (status/=0) call outofmemory('lower','yshiftR')
      allocate(zshiftR(numat),stat=status)
      if (status/=0) call outofmemory('lower','zshiftR')
    endif
!
!  Zero shift arrays
!
    xshiftT(1:numat) = 0.0_dp
    yshiftT(1:numat) = 0.0_dp
    zshiftT(1:numat) = 0.0_dp
!
    if (lrigid) then
      xshiftR(1:numat) = 0.0_dp
      yshiftR(1:numat) = 0.0_dp
      zshiftR(1:numat) = 0.0_dp
    endif
!
    nimag = 1
    do while (freq(nimag).lt.-frqtol.and.nimag.le.mtv)
      nmode = mtvrptr(nimag)
      if (nmode.ne.0) then
        ind = 0
        scale = lowerscale*abs(freq(nimag))
        if (lrigid) then
!
!  Rigid molecule algorithm - no partial occupancy
!
          do i = 1,numatnomol
            jj = numatnomolptr(i)
            rmat(1,4) = scale*eigr(ind+1,nmode)
            rmat(2,4) = scale*eigr(ind+2,nmode)
            rmat(3,4) = scale*eigr(ind+3,nmode)
            if (ndim.eq.3) then
              rmat(1:3,1:3) = rv(1:3,1:3)
              call GULP_gauss(3_i4,3_i4,1_i4,rmat)
            elseif (ndim.eq.2) then
              rmat(1:2,1:2) = rv(1:2,1:2)
              call GULP_gauss(2_i4,3_i4,1_i4,rmat)
            elseif (ndim.eq.1) then
              rmat(1,4) = rv(1,1)
            endif
!
!  Move atom
!
            xshiftT(jj) = xshiftT(jj) + rmat(1,4)
            yshiftT(jj) = yshiftT(jj) + rmat(2,4)
            zshiftT(jj) = zshiftT(jj) + rmat(3,4)
!
!  Move associated shell if there is one
!
            jjs = ncsptr(jj)
            if (jjs.gt.0) then
              xshiftT(jjs) = xshiftT(jjs) + rmat(1,4)
              yshiftT(jjs) = yshiftT(jjs) + rmat(2,4)
              zshiftT(jjs) = zshiftT(jjs) + rmat(3,4)
            endif
            ind = ind + 3
          enddo
!
!  Loop over rigid molecules apply shifts
!
          do nm = 1,nmol
            do n = 1,nmolcore(nm)
              jj = nmollist(nmolptr(nm)+n)
!
!  Translational shift for molecule
!
              rmat(1,4) = scale*eigr(ind+1,nmode)
              rmat(2,4) = scale*eigr(ind+2,nmode)
              rmat(3,4) = scale*eigr(ind+3,nmode)
!
!  Convert Cartesian shift back to fractional
!
              if (ndim.eq.3) then
                rmat(1:3,1:3) = rv(1:3,1:3)
                call GULP_gauss(3_i4,3_i4,1_i4,rmat)
              elseif (ndim.eq.2) then
                rmat(1:2,1:2) = rv(1:2,1:2)
                call GULP_gauss(2_i4,3_i4,1_i4,rmat)
              elseif (ndim.eq.1) then
                rmat(1,4) = rv(1,1)
              endif
!
!  Move atoms
!
              xshiftT(jj) = xshiftT(jj) + rmat(1,4)
              yshiftT(jj) = yshiftT(jj) + rmat(2,4)
              zshiftT(jj) = zshiftT(jj) + rmat(3,4)
!
!  Move associated shell if there is one
!
              jjs = ncsptr(jj)
              if (jjs.gt.0) then
                xshiftT(jjs) = xshiftT(jjs) + rmat(1,4)
                yshiftT(jjs) = yshiftT(jjs) + rmat(2,4)
                zshiftT(jjs) = zshiftT(jjs) + rmat(3,4)
              endif
!
!  Rotational shift - must be cumulative over all rotations
!
              crd(1) = molxyz(1,n,nm) + xshiftR(jj)
              crd(2) = molxyz(2,n,nm) + yshiftR(jj)
              crd(3) = molxyz(3,n,nm) + zshiftR(jj)
              do l = 1,3    ! Loop over axes
                theta = scale*eigr(ind+3+l,nmode)
                call getrotationmatrix(molaxes(1,l,nm),theta,rotm)
                tmp(1:3) = 0.0_dp
                do k = 1,3
                  tmp(1) = tmp(1) + rotm(1,k)*crd(k)
                  tmp(2) = tmp(2) + rotm(2,k)*crd(k)
                  tmp(3) = tmp(3) + rotm(3,k)*crd(k)
                enddo
                crd(1:3) = tmp(1:3)
              enddo
              xshiftR(jj) = crd(1) - molxyz(1,n,nm)
              yshiftR(jj) = crd(2) - molxyz(2,n,nm)
              zshiftR(jj) = crd(3) - molxyz(3,n,nm)
            enddo
            ind = ind + 6
          enddo
        else
          do i = 1,ncfoc
            rmat(1,4) = scale*eigr(ind+1,nmode)
            rmat(2,4) = scale*eigr(ind+2,nmode)
            rmat(3,4) = scale*eigr(ind+3,nmode)
            if (ndim.eq.3) then
              do j = 1,3
                rmat(1,j) = rv(1,j)
                rmat(2,j) = rv(2,j)
                rmat(3,j) = rv(3,j)
              enddo
              call GULP_gauss(3_i4,3_i4,1_i4,rmat)
            elseif (ndim.eq.2) then
              do j = 1,2
                rmat(1,j) = rv(1,j)
                rmat(2,j) = rv(2,j)
              enddo
              call GULP_gauss(2_i4,3_i4,1_i4,rmat)
            elseif (ndim.eq.1) then
              rmat(1,4) = rv(1,1)
            endif
            do j = 1,nphonatc
              if (iocptr(j).eq.i) then
                jj = nphonatptr(j)
                xshiftT(jj) = xshiftT(jj) + rmat(1,4)
                yshiftT(jj) = yshiftT(jj) + rmat(2,4)
                zshiftT(jj) = zshiftT(jj) + rmat(3,4)
!
!  Move associated shells too
!
                jjs = ncsptr(jj)
                if (jjs.gt.0) then
                  xshiftT(jjs) = xshiftT(jjs) + rmat(1,4)
                  yshiftT(jjs) = yshiftT(jjs) + rmat(2,4)
                  zshiftT(jjs) = zshiftT(jjs) + rmat(3,4)
                endif
              endif
            enddo
            ind = ind + 3
          enddo
        endif
!
!  End of condition on mode being local to this processor
!
      endif
      nimag = nimag + 1
    enddo
    nimag = nimag - 1
    if (nimag.gt.0) then
!
!  For rigid molecules apply the rotational component of the imaginary modes and add to translational shifts
!
      if (lrigid) then
        do nm = 1,nmol
          do n = 1,nmolcore(nm)
            jj = nmollist(nmolptr(nm)+n)
            rmat(1,4) = xshiftR(jj)
            rmat(2,4) = yshiftR(jj)
            rmat(3,4) = zshiftR(jj)
!
!  Convert Cartesian shift back to fractional
!
            if (ndim.eq.3) then
              rmat(1:3,1:3) = rv(1:3,1:3)
              call GULP_gauss(3_i4,3_i4,1_i4,rmat)
            elseif (ndim.eq.2) then
              rmat(1:2,1:2) = rv(1:2,1:2)
              call GULP_gauss(2_i4,3_i4,1_i4,rmat)
            elseif (ndim.eq.1) then
              rmat(1,4) = rv(1,1)
            endif
!
!  Move atoms
!
            xshiftT(jj) = xshiftT(jj) + rmat(1,4)
            yshiftT(jj) = yshiftT(jj) + rmat(2,4)
            zshiftT(jj) = zshiftT(jj) + rmat(3,4)
!
!  Move associated shell if there is one
!
            jjs = ncsptr(jj)
            if (jjs.gt.0) then
              xshiftT(jjs) = xshiftT(jjs) + rmat(1,4)
              yshiftT(jjs) = yshiftT(jjs) + rmat(2,4)
              zshiftT(jjs) = zshiftT(jjs) + rmat(3,4)
            endif
          enddo
        enddo
      endif
!
!  Only apply shifts if there were imaginary modes
!
      if (nprocs.gt.1) then
!
!  Allocate workspace arrays for communication
!
        allocate(sum1(3*numat),stat=status)
        if (status/=0) call outofmemory('lower','sum1')
        allocate(sum2(3*numat),stat=status)
        if (status/=0) call outofmemory('lower','sum2')
!
!  Transfer shifts to single array
!
        do i = 1,numat
          sum1(i) = xshiftT(i)
        enddo
        ind = numat
        do i = 1,numat
          sum1(ind+i) = yshiftT(i)
        enddo
        ind = 2*numat
        do i = 1,numat
          sum1(ind+i) = zshiftT(i)
        enddo
!
!  Globalise data
!
        call sumall(sum1,sum2,3_i4*numat,"lower","sum1")
!
!  Transfer global data back to shift arrays
!
        do i = 1,numat
          xshiftT(i) = sum2(i)
        enddo
        ind = numat
        do i = 1,numat
          yshiftT(i) = sum2(ind+i)
        enddo
        ind = 2*numat
        do i = 1,numat
          zshiftT(i) = sum2(ind+i)
        enddo
!
!  Free workspace arrays for communication
!
        deallocate(sum2,stat=status)
        if (status/=0) call deallocate_error('lower','sum2')
        deallocate(sum1,stat=status)
        if (status/=0) call deallocate_error('lower','sum1')
      endif
!
!  Try to ensure that a consistent choice of direction is made since this can be
!  subject to numerical noise and the arbitrary sign of the phonon eigenvectors
!
      ixlarge = 0
      i = 0
      do while (ixlarge.eq.0.and.i.lt.numat)
        i = i + 1
        if (abs(xshiftT(i)).gt.1.0d-3) then
          ixlarge = i
        endif
      enddo
      if (ixlarge.ne.0) then
        scale = sign(1.0_dp,xshiftT(ixlarge))
      else
        iylarge = 0
        i = 0
        do while (iylarge.eq.0.and.i.lt.numat)
          i = i + 1
          if (abs(yshiftT(i)).gt.1.0d-3) then
            iylarge = i
          endif
        enddo
        if (iylarge.ne.0) then
          scale = sign(1.0_dp,yshiftT(iylarge))
        else
          izlarge = 0
          i = 0
          do while (izlarge.eq.0.and.i.lt.numat)
            i = i + 1
            if (abs(zshiftT(i)).gt.1.0d-3) then
              izlarge = i
            endif
          enddo
          if (izlarge.ne.0) then
            scale = sign(1.0_dp,zshiftT(izlarge))
          else
            scale = 1.0_dp
          endif
        endif
      endif
!
!  Multiply scale by lowersign to change direction if requested
!
      scale = scale*lowersign
!
!  Add shifts to configuration arrays
!
      do i = 1,numat
        xcfg(nsft+i) = xcfg(nsft+i) + xshiftT(i)*scale
        ycfg(nsft+i) = ycfg(nsft+i) + yshiftT(i)*scale
        zcfg(nsft+i) = zcfg(nsft+i) + zshiftT(i)*scale
      enddo
    endif
!
!  Set flag to indicate whether lowering has occured
!
    if (nimag.ge.1) then
      llowered(ncf) = .true.
    endif
!
!  Output
!
    if (ioproc) then
      if (nimag.gt.1) then
        write(ioout,'(''  Symmetry lowered according to '',i3,'' imaginary mode eigenvectors'',/)') nimag
      elseif (nimag.eq.1) then
        write(ioout,'(''  Symmetry lowered according to one imaginary mode eigenvector'',/)')
      else
        write(ioout,'(''  No imaginary modes present - current symmetry is correct'',/)')
      endif
    endif
!
!  Free workspace arrays
!
    if (lrigid) then
      deallocate(zshiftR,stat=status)
      if (status/=0) call deallocate_error('lower','zshiftR')
      deallocate(yshiftR,stat=status)
      if (status/=0) call deallocate_error('lower','yshiftR')
      deallocate(xshiftR,stat=status)
      if (status/=0) call deallocate_error('lower','xshiftR')
    endif
!
    deallocate(zshiftT,stat=status)
    if (status/=0) call deallocate_error('lower','zshiftT')
    deallocate(yshiftT,stat=status)
    if (status/=0) call deallocate_error('lower','yshiftT')
    deallocate(xshiftT,stat=status)
    if (status/=0) call deallocate_error('lower','xshiftT')
  endif
#ifdef TRACE
  call trace_out('lower')
#endif
!
  return
  end
