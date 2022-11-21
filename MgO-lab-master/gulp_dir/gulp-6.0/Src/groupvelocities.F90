  subroutine groupvelocities(nk,mtv,eigr,eigi,maxd2,fscale,lprint)
!
!  Calculates the group velocities for the modes at this k point.
!
!   8/14 Created from oscillatorstrength
!  10/14 Restriction on core-only calculations removed
!  12/16 Modified to reduce memory for ceig
!   2/18 Trace added
!   4/20 Mass arrays now use fmass and rfmass
!   5/20 mcv changed to mtv and ncfoc removed from arguments
!
!  On entry :
!
!  nk          = k point for storing group velocities
!  mtv         = no. of modes
!  eigr        = eigenvectors of dynamical matrix - real part
!  eigi        = eigenvectors of dynamical matrix - imaginary part
!  maxd2       = left-hand dimension of eigr
!  fscale      = unit conversion from (eV/Ang**2/amu)^(1/2) to cm^-1
!  lprint      = if true then output group velocities
!
!  On exit :
!
!  groupvelocity = group velocities in cm^-1.Ang
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
!  Julian Gale, CIC, Curtin University, May 2020
!
  use derivatives,    only : derv2dk
  use element
  use frequencies
  use iochannels
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4),      intent(in)    :: maxd2
  integer(i4),      intent(in)    :: mtv
  integer(i4),      intent(in)    :: nk
  logical,          intent(in)    :: lprint
  real(dp),         intent(in)    :: eigr(maxd2,mtv)
  real(dp),         intent(in)    :: eigi(maxd2,mtv)
  real(dp),         intent(in)    :: fscale
!
!  Local variables
!
  integer(i4)                     :: i
  integer(i4)                     :: ii
  integer(i4)                     :: j
  integer(i4)                     :: m
  integer(i4)                     :: status
  complex(dpc)                    :: cmassij
  complex(dpc), allocatable, save :: ceig(:)
  real(dp)                        :: rtmp
#ifdef TRACE
  call trace_in('groupvelocities')
#endif
!
!  Mass weight dynamical matrix derivatives w.r.t. k
!
  do i = 1,mtv
    do j = 1,mtv
      cmassij = dcmplx(rfmass(i)*rfmass(j),0.0_dp)
      derv2dk(1:3,j,i) = cmassij*derv2dk(1:3,j,i)
    enddo
  enddo
!
!  Allocate workspace array
!
  allocate(ceig(mtv),stat=status)
  if (status/=0) call outofmemory('groupvelocities','ceig')
!
!  Loop over modes and project according to eigenvectors
!
  do m = 1,mtv
    groupvelocity(1:3,m,nk) = 0.0_dp
!
!  Create complex eigenvectors
!
    do j = 1,mtv
      ceig(j) = dcmplx(eigr(j,m),eigi(j,m))
    enddo
    do i = 1,mtv
      do j = 1,mtv
        do ii = 1,3
          groupvelocity(ii,m,nk) = groupvelocity(ii,m,nk) + real(conjg(ceig(j))*derv2dk(ii,j,i)*ceig(i))
        enddo
      enddo
    enddo
  enddo
!
!  Free workspace array
!
  deallocate(ceig,stat=status)
  if (status/=0) call deallocate_error('groupvelocities','ceig')
!
!  Divide by conversion factors and 2 x frequency to convert derivatives
!  of frequency squared to those of frequency w.r.t. k
!
  do m = 1,mtv
    if (abs(freq(m,nk)).gt.1.0d-3) then
      rtmp = 0.5_dp*fscale*fscale/freq(m,nk)
      do ii = 1,3
        groupvelocity(ii,m,nk) = rtmp*groupvelocity(ii,m,nk)
      enddo
    else
      do ii = 1,3
        groupvelocity(ii,m,nk) = 0.0_dp
      enddo
    endif
  enddo
  if (lprint) then
!
!  Output group velocities
!
    write(ioout,'('' Group velocities (cm-1.Ang) : '',/)')
    write(ioout,'('' Mode '',8x,''X'',13x,''Y'',13x,''Z''/)')
    do m = 1,mtv
      if (freq(m,nk).gt.0.0_dp) then
        write(ioout,'(i6,1x,3f14.6)') m,(groupvelocity(ii,m,nk),ii=1,3)
      else
        write(ioout,'(i6,1x,3f14.6,'' i'')') m,(groupvelocity(ii,m,nk),ii=1,3)
      endif
    enddo
    write(ioout,'(/)')
  endif
#ifdef TRACE
  call trace_out('groupvelocities')
#endif
!
  return
  end
