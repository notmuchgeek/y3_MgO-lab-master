  subroutine outeem
!
!  Output electronegativity information
!
!   5/18 Created from outspec
!   6/18 e0range added
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
!  Julian Gale, CIC, Curtin University, June 2018
!
  use control,       only : leem
  use element
  use eemdata
  use iochannels
  use parallel,      only : ioproc
  use species
  implicit none
!
!  Local variables
!
  character(len=7) :: wqmax
  character(len=7) :: wqmin
  integer(i4)      :: iele
  integer(i4)      :: nqr
  logical          :: lfound
  logical          :: lout
!
!  If no EEM then return
!
  if (.not.leem.or..not.ioproc) return
  lout = .false.
!
!  Loop over elements looking for those that are present
!
  do iele = 1,maxele
    lfound = any(natspec(1:nspec).eq.iele)
    if (lfound.and.nqrange(iele,neemtype).gt.0) then
!
!  If this is the first time then write out banner
!
      if (.not.lout) then
        lout = .true.
        if (neemtype.eq.2) then
          write(ioout,'(/,''  Charge equilibration parameters for QEq: '',/)')
          write(ioout,'(''--------------------------------------------------------------------------------'')')
          write(ioout,'(''  Element   Chi        Mu         Q0       E0      Radius      Qmin       Qmax  '')')
          write(ioout,'(''            (eV)       (eV)      (au)     (eV)     (Ang)       (au)       (au)  '')')
          write(ioout,'(''--------------------------------------------------------------------------------'')')
        elseif (neemtype.eq.3) then
          write(ioout,'(/,''  Charge equilibration parameters for Streitz and Mintmire: '',/)')
          write(ioout,'(''--------------------------------------------------------------------------------'')')
          write(ioout,'(''  Element   Chi        Mu       Zeta      Znuc     Q0       E0     Qmin    Qmax '')')
          write(ioout,'(''            (eV)      (eV)     (1/Ang)    (au)    (au)     (eV)    (au)    (au) '')')
          write(ioout,'(''--------------------------------------------------------------------------------'')')
        elseif (neemtype.eq.4) then
          write(ioout,'(/,''  Charge equilibration parameters for Pacha: '',/)')
          write(ioout,'(''--------------------------------------------------------------------------------'')')
          write(ioout,'(''  Element   Chi         Mu         Q0        E0         Qmin        Qmax  '')')
          write(ioout,'(''            (eV)       (eV)       (au)      (eV)        (au)        (au)  '')')
          write(ioout,'(''--------------------------------------------------------------------------------'')')
        else
          write(ioout,'(/,''  Charge equilibration parameters for EEM: '',/)')
          write(ioout,'(''--------------------------------------------------------------------------------'')')
          write(ioout,'(''  Element   Chi         Mu         Q0        E0         Qmin        Qmax '')')
          write(ioout,'(''            (eV)       (eV)       (au)      (eV)        (au)        (au) '')')
          write(ioout,'(''--------------------------------------------------------------------------------'')')
        endif
      endif
!
!  Loop over ranges
!
      do nqr = 1,nqrange(iele,neemtype)
!
!  Set range parameters
!
        wqmax = '   -   '
        wqmin = '   -   '
        if (nqrangetype(nqr,iele,neemtype).eq.3) then
          write(wqmin,'(f7.3)') qrangemin(nqr,iele,neemtype)
          write(wqmax,'(f7.3)') qrangemax(nqr,iele,neemtype)
        elseif (nqrangetype(nqr,iele,neemtype).eq.2) then
          write(wqmax,'(f7.3)') qrangemax(nqr,iele,neemtype)
        elseif (nqrangetype(nqr,iele,neemtype).eq.1) then
          write(wqmin,'(f7.3)') qrangemin(nqr,iele,neemtype)
        endif
!
!  Write out parameters
!
        if (neemtype.eq.2) then
          write(ioout,'(3x,a2,1x,f11.5,1x,f11.5,1x,f9.4,1x,f8.4,1x,f7.3,1x,2(4x,a7))') &
            atsym(iele),chirange(nqr,iele,neemtype),murange(nqr,iele,neemtype),q0range(nqr,iele,neemtype), &
            e0range(nqr,iele,neemtype),radrange(nqr,iele,neemtype),wqmin,wqmax
        elseif (neemtype.eq.3) then
          write(ioout,'(3x,a2,1x,f11.5,1x,f11.5,2(1x,f8.4),1x,f8.4,1x,f7.3,2(1x,a7))') &
            atsym(iele),chirange(nqr,iele,neemtype),murange(nqr,iele,neemtype),zetarange(nqr,iele,neemtype), &
            znucrange(nqr,iele,neemtype),q0range(nqr,iele,neemtype),e0range(nqr,iele,neemtype),wqmin,wqmax
        elseif (neemtype.eq.4) then
          write(ioout,'(3x,a2,1x,f11.5,1x,f11.5,1x,f9.4,1x,f8.3,2(5x,a7))') &
            atsym(iele),chirange(nqr,iele,neemtype),murange(nqr,iele,neemtype),q0range(nqr,iele,neemtype), &
            e0range(nqr,iele,neemtype),wqmin,wqmax
        else
          write(ioout,'(3x,a2,1x,f11.5,1x,f11.5,1x,f9.4,1x,f8.3,2(5x,a7))') &
            atsym(iele),chirange(nqr,iele,neemtype),murange(nqr,iele,neemtype),q0range(nqr,iele,neemtype), &
            e0range(nqr,iele,neemtype),wqmin,wqmax
        endif
      enddo
    endif
  enddo
  if (lout) then
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    write(ioout,'(/)')
  endif
!
  return
  end
