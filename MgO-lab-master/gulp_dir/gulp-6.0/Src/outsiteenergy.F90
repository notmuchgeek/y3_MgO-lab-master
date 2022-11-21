  subroutine outsiteenergy(fc)
!
!  Output site energies
!
!   7/20 Created from optout
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
!  Copyright Curtin University, 2020
!
!  Julian Gale, CIC, Curtin University, July 2020
!
  use control,       only : lsiteenergy
  use current
  use energies,      only : siteenergy, epv
  use iochannels
  use parallel,      only : ioproc
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  real(dp),                      intent(in)    :: fc
!
!  Local variables
!
  integer(i4)                                  :: i
  real(dp)                                     :: sum
#ifdef TRACE
  call trace_in('outsiteenergy')
#endif
!
!  Site energies
!
  if (ioproc.and.lsiteenergy) then
    write(ioout,'(/,''  Site energies: '',/)')
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    write(ioout,'(''  Atom No.                Atom energy (eV) '')')
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    sum = 0.0_dp
    do i = 1,numat
      write(ioout,'(i10,4x,f32.8)') i,siteenergy(i)
      sum = sum + siteenergy(i)
    enddo
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    write(ioout,'(''  Sum of site energies  = '',f20.8,'' eV'')') sum
    write(ioout,'(''  Total internal energy = '',f20.8,'' eV'')') fc-epv
    write(ioout,'(''--------------------------------------------------------------------------------'',/)')
  endif
#ifdef TRACE
  call trace_out('outsiteenergy')
#endif
!
  return
  end
