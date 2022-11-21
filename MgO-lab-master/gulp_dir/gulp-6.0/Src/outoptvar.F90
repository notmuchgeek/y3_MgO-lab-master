  subroutine outoptvar(xc)
!
!  Output optimisation variables
!
!   3/20 Created
!   5/20 Parallel modifications made
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
  use current
  use iochannels
  use optimisation
  use parallel,      only : ioproc
  implicit none
!
!  Passed variables
!
  real(dp)         :: xc(*)
!
!  Local variables
!
  integer(i4)      :: i
  integer(i4)      :: ind
!
  if (ioproc) then
!***********************
!  Table of variables  *
!***********************
    write(ioout,'(/,''  Optimisation variables and initial values:'',/)')
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    write(ioout,'(''     Variable No.       Variable type             Index        Initial value    '')')
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    do i = 1,nvar
      ind = ioptindex(i)
      if (iopttype(i).eq.iopt_cell) then
        write(ioout,'(4x,i10,10x,''Cell parameter      '',6x,i8,f16.6)') i,ind,xc(i)
      elseif (iopttype(i).eq.iopt_strain) then
        write(ioout,'(4x,i10,10x,''Cell strain         '',6x,i8,f16.6)') i,ind,xc(i)
      elseif (iopttype(i).eq.iopt_xf) then
        write(ioout,'(4x,i10,10x,''X coordinate        '',6x,i8,f16.6)') i,ind,xc(i)
      elseif (iopttype(i).eq.iopt_yf) then
        write(ioout,'(4x,i10,10x,''Y coordinate        '',6x,i8,f16.6)') i,ind,xc(i)
      elseif (iopttype(i).eq.iopt_zf) then
        write(ioout,'(4x,i10,10x,''Z coordinate        '',6x,i8,f16.6)') i,ind,xc(i)
      elseif (iopttype(i).eq.iopt_radius) then
        write(ioout,'(4x,i10,10x,''Radius              '',6x,i8,f16.6)') i,ind,xc(i)
      elseif (iopttype(i).eq.iopt_xcom) then
        write(ioout,'(4x,i10,10x,''X centre of mass    '',6x,i8,f16.6)') i,ind,xc(i)
      elseif (iopttype(i).eq.iopt_ycom) then
        write(ioout,'(4x,i10,10x,''Y centre of mass    '',6x,i8,f16.6)') i,ind,xc(i)
      elseif (iopttype(i).eq.iopt_zcom) then
        write(ioout,'(4x,i10,10x,''Z centre of mass    '',6x,i8,f16.6)') i,ind,xc(i)
      elseif (iopttype(i).eq.iopt_xqtn) then
        write(ioout,'(4x,i10,10x,''Quaternion about X  '',6x,i8,f16.6)') i,ind,xc(i)
      elseif (iopttype(i).eq.iopt_yqtn) then
        write(ioout,'(4x,i10,10x,''Quaternion about Y  '',6x,i8,f16.6)') i,ind,xc(i)
      elseif (iopttype(i).eq.iopt_zqtn) then
        write(ioout,'(4x,i10,10x,''Quaternion about Z  '',6x,i8,f16.6)') i,ind,xc(i)
      endif
    enddo
    write(ioout,'(''--------------------------------------------------------------------------------'',/)')
  endif
!
  return
  end
