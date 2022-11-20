  subroutine gruneisen(nk,mcv,msv,ncfoc,eigr,eigi,dervr,dervi,maxd2,iocptr,fscale,lprint)
!
!  Calculates the Grueneisen parameters for the modes at this k point.
!
!   1/18 Created from groupvelocities and fenergy3
!   2/18 Trace added
!   4/20 Mass arrays now use fmass and rfmass
!
!  On entry :
!
!  nk          = k point for storing group velocities
!  mcv         = no. of modes
!  msv         = no. of shell coordinates
!  ncfoc       = number of condensed core sites
!  iocptr      = pointer to fully occupied sites
!  eigr        = eigenvectors of dynamical matrix - real part
!  eigi        = eigenvectors of dynamical matrix - imaginary part
!  dervr       = contains matrices needed for shells - real part
!  dervi       = contains matrices needed for shells - imaginary part
!  maxd2       = left-hand dimension of eigr
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
  use current,     only : rv
  use element
  use frequencies
  use iochannels
  use numbers,     only : third
#ifdef TRACE
  use trace,       only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4),      intent(in)    :: maxd2
  integer(i4),      intent(in)    :: mcv
  integer(i4),      intent(in)    :: msv
  integer(i4),      intent(in)    :: nk
  integer(i4),      intent(in)    :: ncfoc
  integer(i4),      intent(in)    :: iocptr(*)
  logical,          intent(in)    :: lprint
  real(dp),         intent(inout) :: dervr(maxd2,*)
  real(dp),         intent(inout) :: dervi(maxd2,*)
  real(dp),         intent(in)    :: eigr(maxd2,mcv)
  real(dp),         intent(in)    :: eigi(maxd2,mcv)
  real(dp),         intent(in)    :: fscale
!
!  Local variables
!
  integer(i4)                     :: i
  integer(i4)                     :: j
  integer(i4)                     :: m
  real(dp)                        :: rtmp
  real(dp)                        :: vol
  real(dp)                        :: volume
#ifdef TRACE
  call trace_in('grueneisen')
#endif
!
!  Calculate shell projection matrices
!
!  Msc => already stored in dervr(mcv+j,i)/dervi(mcv+j,i)
!  Pns => store in dervr(i,mcv+j)/dervi(i,mcv+j)
!      =  Enc.Mcs, where Enc = eigenvector for mode n
!
  if (msv.gt.0) then
!
!  Scale Msc by mass factor for use in derivatives
!
    do i = 1,3*ncfoc
      do j = 1,msv
        dervr(mcv+j,i) = rfmass(i)*dervr(mcv+j,i)
        dervi(mcv+j,i) = rfmass(i)*dervi(mcv+j,i)
      enddo
    enddo
!
!  Generate Pns from Msc and eigenvectors
!
    call setfeshtmat(mcv,msv,dervr,dervi,eigr,eigi,maxd2)
  endif
!********************************
!  Evaluate phonon derivatives  *
!********************************
  call realrecip3d3dV(nk,mcv,eigr,eigi,maxd2,iocptr,grueneisen(1,nk))
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
  if (lprint) then
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
#ifdef TRACE
  call trace_out('grueneisen')
#endif
!
  return
  end
