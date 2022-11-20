  subroutine gulp_banner
!
!  Outputs banner for GULP
!
!   4/01 Created from GULP main routine for tidiness
!  10/02 Version incremented to 1.4.1 to signify SE changes
!  10/02 Date of modifcation added to banner
!   9/06 Number incremented due to major torsion modification
!   3/07 Number incremented due to switch to f90 format
!   6/09 Renamed to gulp_banner for benefit of Chemshell
!   8/11 Updated to MS studio 6.0 and GULP 4.0.
!  11/13 Updated to MS studio 7.0 and GULP 4.2.
!   8/15 Updated to MS studio 9.0 and GULP 4.4.
!   9/16 Updated to MS studio 10.0 and GULP 4.6.
!   4/17 Changed to GULP 5.0
!  12/17 Changed to GULP 5.1
!   1/18 Trace added
!   8/18 Changed to GULP 5.2
!   8/19 Changed to GULP 5.3
!   7/20 Changed to GULP 6.0
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
!  Copyright Curtin University 2021
!
!  Julian Gale, Curtin University, February 2021
!
  use iochannels
  use parallel
#ifdef TRACE
  use trace,      only : trace_in, trace_out
#endif
  implicit none
#ifdef TRACE
  call trace_in('banner')
#endif
!******************
!  Output header  *
!******************
  if (ioproc) then
#ifdef ACCELRYS
    write(ioout,'(''********************************************************************************'')')
    write(ioout,'(''*                       GENERAL UTILITY LATTICE PROGRAM  v6.0                  *'')')
    write(ioout,'(''*                                 Julian Gale                                  *'')')
    write(ioout,'(''*                       Curtin Institute for Computation                       *'')')
    write(ioout,'(''*                    School of Molecular and Life Sciences                     *'')')
    write(ioout,'(''*                    Curtin University, Western Australia                      *'')')
    write(ioout,'(''********************************************************************************'')')
    write(ioout,'(''*                     BIOVIA Materials Studio 2021 Release                     *'')')
    write(ioout,'(''********************************************************************************'')')
#else
    write(ioout,'(''********************************************************************************'')')
    write(ioout,'(''*                       GENERAL UTILITY LATTICE PROGRAM                        *'')')
    write(ioout,'(''*                                 Julian Gale                                  *'')')
    write(ioout,'(''*                       Curtin Institute for Computation                       *'')')
    write(ioout,'(''*                    School of Molecular and Life Sciences                     *'')')
    write(ioout,'(''*                    Curtin University, Western Australia                      *'')')
    write(ioout,'(''********************************************************************************'')')
    write(ioout,'(''* Version = 6.0.0 * Last modified =   8th February 2021                        *'')')
    write(ioout,'(''********************************************************************************'')')
#endif
  endif
#ifdef TRACE
  call trace_out('banner')
#endif
!
  return
  end
