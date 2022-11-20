  subroutine checkquaternions(lgeometryOK,lgeoreset)
!
!  Checks quaternions and resets if required
!
!   7/20 Created
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
  use control,        only : keyword
  use current
  use iochannels
  use molecule
  use optimisation,   only : iopt_xqtn, iopt_yqtn, iopt_zqtn
  use parallel,       only : ioproc
  use xcgc,           only : xc
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  logical,                       intent(in)    :: lgeoreset   ! If true then reset xc
  logical,                       intent(out)   :: lgeometryOK ! Quaternions were OK
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ind
  integer(i4)                                  :: nm
  real(dp)                                     :: e02
  real(dp)                                     :: e1
  real(dp)                                     :: e2
  real(dp)                                     :: e3
#ifdef TRACE
  call trace_in('checkquaternions')
#endif
!***************************************************
!  Loop over rigid molecules to check quaternions  *
!***************************************************
  lgeometryOK = .true.
  do nm = 1,nmol
    e1 = molQa(1,nm)
    e2 = molQa(2,nm)
    e3 = molQa(3,nm)
    e02 = 1.0_dp - e1**2 - e2**2 - e3**2
    if (e02.lt.0.01_dp) then
!
!  First quaternion component is close to going complex => reset
!
      lgeometryOK = .false.
      if (lgeoreset) then
        call resetrigidmol(nm)
        if (ioproc.and.index(keyword,'verb').ne.0) then
          write(ioout,'(/,''  Quaternions for molecule '',i4,'' reset as q0 squared = '',f12.6,/)') nm,e02
        endif
      endif
    endif
  enddo
  if (lgeoreset.and..not.lgeometryOK) then
!***************************************************
!  Loop over variables to reset quaternion values  *
!***************************************************
    do i = 1,nvar
      ind = ioptindex(i)
      if (iopttype(i).eq.iopt_xqtn) then
!
!  x molecule quaternion
!
        xc(i) = molQa(1,ind)
      elseif (iopttype(i).eq.iopt_yqtn) then
!
!  y molecule quaternion
!
        xc(i) = molQa(2,ind)
      elseif (iopttype(i).eq.iopt_zqtn) then
!
!  z molecule quaternion
!
        xc(i) = molQa(3,ind)
      endif
    enddo
  endif
#ifdef TRACE
  call trace_out('checkquaternions')
#endif
!
  return
  end
!
  subroutine resetquaternions
!
!  Resets quaternions 
!
!   7/20 Created from checkquaternions
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
  use current
  use molecule
  use optimisation,   only : iopt_xqtn, iopt_yqtn, iopt_zqtn
  use xcgc,           only : xc
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ind
  integer(i4)                                  :: nm
#ifdef TRACE
  call trace_in('resetquaternions')
#endif
!***************************************************
!  Loop over rigid molecules to check quaternions  *
!***************************************************
  do nm = 1,nmol
    call resetrigidmol(nm)
  enddo
!****************************************************
!  Loop over variables to reset quaternions values  *
!****************************************************
  do i = 1,nvar
    ind = ioptindex(i)
    if (iopttype(i).eq.iopt_xqtn) then
!
!  x molecule quaternion
!
      xc(i) = molQa(1,ind)
    elseif (iopttype(i).eq.iopt_yqtn) then
!
!  y molecule quaternion
!
      xc(i) = molQa(2,ind)
    elseif (iopttype(i).eq.iopt_zqtn) then
!
!  z molecule quaternion
!
      xc(i) = molQa(3,ind)
    endif
  enddo
#ifdef TRACE
  call trace_out('resetquaternions')
#endif
!
  return
  end
