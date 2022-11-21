

! Dummy routines to link against a library that 
! support ChemShell builds

  function opc_root()
  implicit none
  logical                   :: opc_root

  opc_root=.true.

  return 
  end function opc_root

!----------------------------------------------------------------      
  subroutine GetGulpFileNames(infile, outfile, len )
  use datatypes
  implicit none

  character(*),intent(out) :: infile
  character(*),intent(out) :: outfile
  integer(i4),intent(in) :: len

  infile = "stdout"
  outfile = "stdin"

  return
  end subroutine GetGulpFileNames
!----------------------------------------------------------------      

#ifdef OLDCS
  subroutine GetShellForces(shell_force,iret,nsh)
  use datatypes
  implicit none

  integer(i4),intent(out) :: iret
  integer(i4),intent(out) :: nsh
  real(dp),dimension(:,:),intent(out) :: shell_force

  iret = 0_i4
  nsh = 0_i4
  shell_force = 0.0_dp

  return
  end subroutine GetShellForces
!----------------------------------------------------------------

  subroutine ExportGulpCharges(ngq, gq)
  use datatypes
  implicit none

  integer(i4),intent(in) :: ngq
  real(dp),dimension(:),intent(in) :: gq

  return
  end subroutine ExportGulpCharges
!----------------------------------------------------------------
#endif
