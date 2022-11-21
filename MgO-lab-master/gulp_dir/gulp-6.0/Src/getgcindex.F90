  subroutine getgcindex(i,igc,ix,iy,iz)
!
!  Convert atom index in full cell to equivalent index in the original ghost cell
!
!   9/20 Created
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
!  Julian Gale, CIC, Curtin University, September 2020
!
  use configurations
  use current
  use shells,         only : ncore
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4),                     intent(in)  :: i
  integer(i4),                     intent(out) :: igc
  integer(i4),                     intent(out) :: ix
  integer(i4),                     intent(out) :: iy
  integer(i4),                     intent(out) :: iz
!
!  Local variables
!
  integer(i4)                                  :: image
  integer(i4)                                  :: ngatoms
  integer(i4)                                  :: ngcells
  integer(i4)                                  :: nx
  integer(i4)                                  :: ny
  integer(i4)                                  :: nz
#ifdef TRACE
  call trace_in('getgcindex')
#endif
!
  nx = nsuperghost(1,ncf)
  ny = nsuperghost(1,ncf)
  nz = nsuperghost(1,ncf)
!
!  Set total number of ghostcells
!
  ngcells = nx*ny*nz
  ngatoms = ncore/ngcells
!
!  Divide atom by number of ghostcells to find atom in ghost cell
!
  igc = (i-1)/ngcells + 1
!
!  Find image number
!
  image = i - (igc-1)*ngatoms
!
!  Convert image into cell indices
!
  if (ordersuper(ncf).eq.2) then
!
!  Order of generator loops = Z -> Y -> X
!
    iz = (image-1)/nx*ny + 1
    image = image - (iz-1)*nx*ny
    iy = (image-1)/nx + 1
    ix = image - (iy-1)*nx
  else
!
!  Order of generator loops = X -> Y -> Z
!
    ix = (image-1)/ny*nz + 1
    image = image - (ix-1)*ny*nz
    iy = (image-1)/nz + 1
    iz = image - (iy-1)*nz
  endif
#ifdef TRACE
  call trace_out('getgcindex')
#endif
!
  return
  end
