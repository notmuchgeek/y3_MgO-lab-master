  subroutine ramanstrength(mcv,mcvloc,mcvptr,nphonatc,nphonatptr,ncfoc,iocptr,eigr,maxd2,ramstrength)
!
!  Calculates the Raman strengths for the modes. Gamma point version.
!
!  10/13 Created by separating from oscillatorstrengthg
!  12/16 Modified for distributed memory parallelism
!   1/18 Modified for ghostcell case
!   2/18 Trace added
!
!  On entry :
!
!  mcv         = no. of modes
!  mcvloc      = no. of modes on local node
!  mcvptr      = pointer from local to global mode
!  nphonatc    = total number of cores in relevant regions
!  nphonatptr  = pointer from reduce to full atom set
!  ncfoc       = number of condensed core sites
!  iocptr      = pointer that connects sites to condensed sites
!  eigr        = eigenvectors of dynamical matrix
!  maxd2       = left-hand dimension of eigr
!
!  On exit : 
!
!  ramstrength = Raman susceptibility tensors, if lraman = .true.
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
!  Julian Gale, CIC, Curtin University, Feb 2018
!
  use configurations, only : nsuperghost
  use control,        only : lraman, lghost
  use current
  use element
  use iochannels
  use parallel,       only : nprocs
  use properties,     only : ramanasus
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4),                   intent(in)    :: iocptr(*)
  integer(i4),                   intent(in)    :: maxd2
  integer(i4),                   intent(in)    :: mcv
  integer(i4),                   intent(in)    :: mcvloc
  integer(i4),                   intent(in)    :: mcvptr(mcvloc)
  integer(i4),                   intent(in)    :: nphonatc
  integer(i4),                   intent(in)    :: nphonatptr(*)
  integer(i4),                   intent(in)    :: ncfoc
  real(dp),                      intent(in)    :: eigr(maxd2,mcv)
  real(dp),                      intent(out)   :: ramstrength(3,3,mcv)
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: iocj
  integer(i4)                                  :: ix
  integer(i4)                                  :: iy
  integer(i4)                                  :: iz
  integer(i4)                                  :: j
  integer(i4)                                  :: m
  integer(i4)                                  :: mm
  integer(i4)                                  :: nghostcell
  integer(i4)                                  :: status
  real(dp),     dimension(:,:,:), allocatable  :: sum3
#ifdef TRACE
  call trace_in('ramanstrength')
#endif
!
  if (lraman) then
!**********************************
!  Oscillator strengths per mode  *
!**********************************
!
!  Zero arrays
!
    ramstrength(1:3,1:3,1:mcv) = 0.0_dp
!
!  Find number of ghost cells
!
    if (lghost) then
      nghostcell = nsuperghost(1,ncf)*nsuperghost(2,ncf)*nsuperghost(3,ncf)
    else
      nghostcell = 1
    endif
!
!  Loop over modes
!
    do m = 1,mcvloc
      mm = mcvptr(m)
!
!  Loop over full sites
!
      do i = 1,ncfoc
        ix = 3*(i-1) + 1
        iy = ix + 1
        iz = ix + 2
!
!  Find all cores associated with full site
!
        do j = 1,nphonatc,nghostcell
          iocj = (iocptr(j) - 1)/nghostcell + 1 
          if (iocj.eq.i) then
!
!  Multiple inverse mass weighted eigenvectors by Raman susceptibility
!
            ramstrength(1:3,1:3,mm) = ramstrength(1:3,1:3,mm) + ramanasus(1:3,1:3,1,j)*eigr(ix,m) &
                                                              + ramanasus(1:3,1:3,2,j)*eigr(iy,m) &
                                                              + ramanasus(1:3,1:3,3,j)*eigr(iz,m) 
          endif
        enddo
      enddo
    enddo
    if (nprocs.gt.1) then
!
!  Globalise oscillator strengths
!
      allocate(sum3(3,3,mcv),stat=status)
      if (status/=0) call outofmemory('ramanstrength','sum3')
!
      call sumall(ramstrength,sum3,9_i4*mcv,"ramanstrength","ramstrength")
      ramstrength(1:3,1:3,1:mcv) = sum3(1:3,1:3,1:mcv)
!
      deallocate(sum3,stat=status)
      if (status/=0) call deallocate_error('ramanstrength','sum3')
    endif
  else
!
!  Loop over modes
!
    ramstrength(1:3,1:3,1:mcv) = 0.0_dp
  endif
#ifdef TRACE
  call trace_out('ramanstrength')
#endif
!
  return
  end
