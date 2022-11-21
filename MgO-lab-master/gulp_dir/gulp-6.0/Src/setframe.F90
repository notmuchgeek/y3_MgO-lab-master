  subroutine setframe
!
!  This routine rotates the coordinates of a finite system to be
!  in the frame of reference for the moment of inertia tensor.
!
!   3/18 Created from part of setsas
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
!  Julian Gale, CIC, Curtin University, March 2018
!
  use configurations,   only : xcfg, ycfg, zcfg
  use control,          only : ldebug
  use current
  use iochannels
  use parallel,         only : ioproc
  use shells,           only : ncore, ncoptr
#ifdef TRACE
  use trace,            only : trace_in, trace_out
#endif
  implicit none
!
!  Local variables
!
  integer(i4)                                    :: i
  integer(i4)                                    :: j
  integer(i4)                                    :: ierror
  integer(i4)                                    :: ix
  integer(i4)                                    :: nc
  real(dp)                                       :: eig(3)
  real(dp)                                       :: eigen(3,3)
  real(dp)                                       :: fv1(9)
  real(dp)                                       :: moments(3,3)
  real(dp)                                       :: originx
  real(dp)                                       :: originy
  real(dp)                                       :: originz
  real(dp)                                       :: xc, yc, zc
  real(dp)                                       :: xo, yo, zo
  real(dp)                                       :: xyz(3)
!
!  If this is not a molecule then return
!
  if (ndim.ne.0) return
!
#ifdef TRACE
  call trace_in('setframe')
#endif
!
!  Find centre of system
!
  originx = 0.0_dp
  originy = 0.0_dp
  originz = 0.0_dp
  do nc = 1,ncore
    i = ncoptr(nc)
    originx = originx + xclat(i)
    originy = originy + yclat(i)
    originz = originz + zclat(i)
  enddo
  originx = originx/dble(ncore)
  originy = originy/dble(ncore)
  originz = originz/dble(ncore)
!
!  Calculate moments about centre
!
  moments(1:3,1:3) = 0.0_dp
  do nc = 1,ncore
    i = ncoptr(nc)
    xc = xclat(i) - originx
    yc = yclat(i) - originy
    zc = zclat(i) - originz
    moments(1,1) = moments(1,1) + xc*xc
    moments(2,1) = moments(2,1) + yc*xc
    moments(3,1) = moments(3,1) + zc*xc
    moments(2,2) = moments(2,2) + yc*yc
    moments(3,2) = moments(3,2) + zc*yc
    moments(3,3) = moments(3,3) + zc*zc
  enddo
  moments(1,2) = moments(2,1)
  moments(1,3) = moments(3,1)
  moments(2,3) = moments(3,2)
!
!  Find eigensystem of moments
!
  eigen(1:3,1:3) = moments(1:3,1:3)
  call dsyev('V','U',3_i4,eigen,3_i4,eig,fv1,9_i4,ierror)
!
!  Debugging output
!
  if (ldebug.and.ioproc) then
    write(ioout,'(/)')
    write(ioout,'(''  Frame transformation matrix and eigenvalues :'',/)')
    do i = 1,3
      write(ioout,'(''  No. '',i1,'' E = '',f12.6,'' TM = '',3(f12.6,1x))') i,eig(i),(eigen(j,i),j=1,3)
    enddo
    write(ioout,'(/)')
  endif
!
!  Place the coordinates of all atoms in the reference frame
!
  do i = 1,numat
    xc = xclat(i)
    yc = yclat(i)
    zc = zclat(i)
!
!  Shift relative to the origin
!
    xo = xc - originx
    yo = yc - originy
    zo = zc - originz
!
!  Calculate dot product of vector to the origin with moment eigenvectors
!
    do ix = 1,3
      xyz(ix) = xo*eigen(1,ix) + yo*eigen(2,ix) + zo*eigen(3,ix)
    enddo
!
!  Place back in main arrays
!
    xclat(i) = xyz(1)
    yclat(i) = xyz(2)
    zclat(i) = xyz(3)
!
    xafrac(i) = xclat(i)
    yafrac(i) = yclat(i)
    zafrac(i) = zclat(i)
!
    xcfg(nsft+i) = xclat(i)
    ycfg(nsft+i) = yclat(i)
    zcfg(nsft+i) = zclat(i)
!
!  End of loop over atoms
!
  enddo 
#ifdef TRACE
  call trace_out('setframe')
#endif
!
  return
  end
