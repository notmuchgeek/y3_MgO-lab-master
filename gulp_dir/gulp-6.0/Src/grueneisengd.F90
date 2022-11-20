  subroutine gruneisengd(nk,mcv,mcvloc,mcvptr,msv,ncfoc,eigr,maxeigr,derv2,maxd2, &
                         fscale,lprint)
!
!  Calculates the Grueneisen parameters for the modes at this k point.
!  Gamma point version. Distributed memory parallel version.
!
!   2/18 Created from gruneiseng
!   4/20 Mass arrays now use fmass and rfmass
!
!  On entry :
!
!  nk          = k point for storing group velocities
!  mcv         = no. of modes
!  mcvloc      = no. of modes on local node
!  mcvptr      = pointer from local to global mode
!  msv         = no. of shell coordinates
!  ncfoc       = number of condensed core sites
!  eigr        = eigenvectors of dynamical matrix
!  maxeigr     = left-hand dimension of eigr
!  derv2       = contains matrices needed for shells
!  maxd2       = left-hand dimension of derv2
!  fscale      = unit conversion from (eV/Ang**2/amu)^(1/2) to cm^-1
!  lprint      = if true then output Grueneisen parameters
!
!  On exit :
!
!  grueneisen = mode Grueneisen parameters (unit less)
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
!  Julian Gale, CIC, Curtin University, April 2020
!
#ifdef MPI
  use current,     only : rv
#endif
  use element
  use frequencies
  use iochannels
#ifdef MPI
  use numbers,     only : third
#endif
  use parallel
#ifdef TRACE
  use trace,       only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4),                    intent(in)    :: maxd2
  integer(i4),                    intent(in)    :: maxeigr
  integer(i4),                    intent(in)    :: mcv
  integer(i4),                    intent(in)    :: mcvloc
  integer(i4),                    intent(in)    :: mcvptr(*)
  integer(i4),                    intent(in)    :: msv
  integer(i4),                    intent(in)    :: nk
  integer(i4),                    intent(in)    :: ncfoc
  logical,                        intent(in)    :: lprint
  real(dp),                       intent(inout) :: derv2(maxd2,*)
  real(dp),                       intent(in)    :: eigr(maxeigr,*)
  real(dp),                       intent(in)    :: fscale
#ifdef MPI
!
!  Local variables
!
  integer                                       :: idesd(9)
  integer                                       :: idesr(9)
  integer                                       :: ifails
  integer                                       :: ld
  integer                                       :: nb
  integer                                       :: ncs
!
  integer(i4)                                   :: i
  integer(i4)                                   :: indi
  integer(i4)                                   :: ix
  integer(i4)                                   :: iy
  integer(i4)                                   :: iz
  integer(i4)                                   :: j
  integer(i4)                                   :: m
  integer(i4)                                   :: status
  real(dp)                                      :: rtmp
  real(dp)                                      :: vol
  real(dp)                                      :: volume
  real(dp),     dimension(:,:), allocatable     :: wrk2
#endif
#ifdef TRACE
  call trace_in('grueneisengd')
#endif
#ifdef MPI
!
!  Calculate shell projection matrices
!
!  Msc => already stored in derv2(mcv+j,i)/dervi(mcv+j,i)
!  Pns => store in derv2(i,mcv+j)/dervi(i,mcv+j)
!      =  Enc.Mcs, where Enc = eigenvector for mode n
!
  if (msv.gt.0) then
!
!  Scale Msc by mass factor for use in derivatives
!
    ix = - 2
    iy = - 1
    iz =   0
    do i = 1,ncoreonnode
      indi = 3*(node2atom(i)-1)
      ix = ix + 3
      iy = iy + 3
      iz = iz + 3
      do j = 1,msv
        derv2(mcv+j,ix) = rfmass(indi+1)*derv2(mcv+j,ix)
        derv2(mcv+j,iy) = rfmass(indi+2)*derv2(mcv+j,iy)
        derv2(mcv+j,iz) = rfmass(indi+3)*derv2(mcv+j,iz)
      enddo
    enddo
!
!  Set up Blacs descriptors for matrices
!
    nb = nblocksize
    ifails = 0
    ncs = mcv + msv
    ld = maxeigr
    call descinit( idesr, ncs, ncs, 3*nb, 3*nb, 0, 0, iBlacsContext, ld, ifails )
!
    if (ifails.ne.0) then
      call outerror('initialisation in descinit failed',0_i4)
      call stopnow('grueneisengd')
    endif
!
    ld = maxd2
    call descinit( idesd, ncs, ncs, 3*nb, 3*nb, 0, 0, iBlacsContext, ld, ifails )
!
    if (ifails.ne.0) then
      call outerror('initialisation in descinit failed',0_i4)
      call stopnow('grueneisengd')
    endif
!
!  Generate Pns from Msc and eigenvectors
!  NB position in derv2 of final matrix is flipped relative to old position
!
    allocate(wrk2(maxd2,mcvloc),stat=status)
    if (status/=0) call outofmemory('grueneisengd','wrk2')
!
    call pdgemm('N','N',msv,mcv,mcv,1.0d0,derv2,mcv+1,1,idesd,eigr,1,1,idesr,0.0d0,wrk2,1,1,idesd)
!
!  Copy wrk2 back to derv2
!
    do i = 1,mcvloc
      do j = 1,msv
        derv2(mcv+j,i) = wrk2(j,i)
      enddo
    enddo
!
    deallocate(wrk2,stat=status)
    if (status/=0) call deallocate_error('grueneisengd','wrk2')
  endif
!********************************
!  Evaluate phonon derivatives  *
!********************************
  call realrecip3d3dVgd(nk,mcv,mcvloc,mcvptr,eigr,maxeigr,grueneisen(1,nk))
!
!  Compute volume
!
  vol = volume(rv)
!
!  Convert units
!
  do m = 1,mcv
    if (abs(freq(m,nk)).gt.1.0d-3) then
      rtmp = 0.5_dp*third*fscale*fscale/freq(m,nk)**2
      grueneisen(m,nk) = rtmp*grueneisen(m,nk)
    else
      grueneisen(m,nk) = 0.0_dp
    endif
  enddo
  if (lprint.and.ioproc) then
!
!  Output Grueneisen parameters
!
    write(ioout,'(/,''  Grueneisen parameters : '',/)')
    write(ioout,'(''  Mode '',4x,'' Frequency '',10x,''G'')')
    write(ioout,'(11x,''  (cm-1)   '',6x,''(unitless)'',/)')
    do m = 1,mcv
      if (freq(m,nk).gt.0.0_dp) then
        write(ioout,'(i6,1x,f14.6,3x,f14.6)') m,freq(m,nk),grueneisen(m,nk)
      elseif (abs(freq(m,nk)).gt.1.0d-2) then
        write(ioout,'(i6,1x,f14.6,'' i '',f14.6)') m,freq(m,nk),grueneisen(m,nk)
      else
        write(ioout,'(i6,1x,f14.6,3x,f14.6)') m,freq(m,nk),grueneisen(m,nk)
      endif
    enddo
    write(ioout,'(/)')
  endif
#else
  call outerror('grueneisengd called when not compiled with MPI',0_i4)
  call stopnow('grueneisengd')
#endif
#ifdef TRACE
  call trace_out('grueneisengd')
#endif
!
  return
  end
