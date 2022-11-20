  subroutine outofmemory(routine,arrayname)
!
!  Handles an error due to memory exceeded
!
!   9/16 lclose argument no longer passed to finish
!
!  Julian Gale, Curtin University, September 2016
!
  use iochannels
  use parallel
!
  implicit none
  character(len=*) :: routine, arrayname
!
!  Write out error message
!
  call outerror('memory allocation failed',0_i4)
  write(ioout,'(''  Calling routine = '',a)') routine
  write(ioout,'(''  Array name      = '',a,/)') arrayname
!
!  Perform final tasks before stopping
!
  call finish
!
  end

  subroutine outofsize(routine,arrayname)
!
!  Handles an error due to memory allocation failure due to integer precision
!
!   9/16 lclose argument no longer passed to finish
!
!  Julian Gale, Curtin University, September 2016
!
  use iochannels
  use parallel
!
  implicit none
  character(len=*) :: routine, arrayname
!
!  Write out error message
!
  call outerror('memory allocation failed due to integer data type limits',0_i4)
  write(ioout,'(''  Calling routine = '',a)') routine
  write(ioout,'(''  Array name      = '',a,/)') arrayname
!
!  Perform final tasks before stopping
!
  call finish
!
  end

  subroutine deallocate_error(routine,arrayname)
!
!  Handles an error at deallocation time
!
!   9/16 lclose argument no longer passed to finish
!
!  Victor Milman, Accelrys, September 2016
!
  use iochannels
  use parallel
!
  implicit none
  character(len=*) :: routine, arrayname
!
!  Write out error message
!
  call outerror('memory deallocation failed',0_i4)
  write(ioout,'(''  Calling routine = '',a)') routine
  write(ioout,'(''  Array name      = '',a,/)') arrayname
!
!  Perform final tasks before stopping
!
  call finish
!
  end
