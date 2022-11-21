  subroutine outstresses(lprint)
!
!  Calculate the stresses and output them
!
!   5/12 Created from property
!  12/13 Calculation of stresses moved to strfin
!   7/17 Condition added to check that this is a 3-D system
!        otherwise the stresses are not yet defined
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
!  Copyright Curtin University 2017
!
!  Julian Gale, CIC, Curtin University, July 2017
!
  use control,        only : lstressout
  use current,        only : ndim
  use derivatives,    only : stresses
  use iochannels
  use parallel,       only : ioproc
  implicit none
!
!  Passed variables
!
  logical, intent(in)   :: lprint
!
  if (.not.lstressout.or.ndim.ne.3) return
!*************************
!  Stress tensor output  *
!*************************
  if (lprint.and.ioproc) then
    write(ioout,'(''  Final stress tensor components (GPa):'',/)')
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    write(ioout,'(''       xx '',1x,f14.6,''    yz '',1x,f14.6)') stresses(1),stresses(4)
    write(ioout,'(''       yy '',1x,f14.6,''    xz '',1x,f14.6)') stresses(2),stresses(5)
    write(ioout,'(''       zz '',1x,f14.6,''    xy '',1x,f14.6)') stresses(3),stresses(6)
    write(ioout,'(''--------------------------------------------------------------------------------'',/)')
  endif
!
!  Flush the output buffer
!
  call gflush(ioout)
!
  return
  end
