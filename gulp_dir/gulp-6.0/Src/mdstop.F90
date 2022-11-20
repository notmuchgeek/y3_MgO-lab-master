  subroutine mdstop
!
!  This routine is called when a job runs out of time
!
!   9/16 lclose argument no longer passed to finish
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
!  Copyright Curtin University 2016
!
  use iochannels
  use parallel
  implicit none
!
  if (ioproc) then
    write(ioout,'(/,''  **** Insufficent time for another cycle of molecular dynamics ****'',/)')
  endif
  call finish
!
  stop
  end
