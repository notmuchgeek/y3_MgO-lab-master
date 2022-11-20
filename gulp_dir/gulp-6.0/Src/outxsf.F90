  subroutine outxsf(iout)
!
!  Write out XCrySDen XSF files 
!
!   4/18 Created from outxtl
!   6/18 Parallel I/O trap added
!   1/19 maxwordlength changes added
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
!  Copyright Curtin University 2019
!
!  Julian Gale, CIC, Curtin University, January 2019
!
  use control
  use configurations
  use current
  use element
  use gulp_files
  use gulp_lengths
  use general
  use iochannels
  use parallel,         only : ioproc
  use shells
  use symmetry
  implicit none
!
!  Passed variables
!
  integer(i4)              :: iout
!
!  Local variables
!
  character(len=1)             :: numbers(10)
  character(len=maxwordlength) :: xsffile2
  character(len=maxwordlength) :: xsffilel
  integer(i4)                  :: i
  integer(i4)                  :: ii
  integer(i4)                  :: incf
  integer(i4)                  :: ind
  integer(i4)                  :: j
  integer(i4)                  :: jj
  integer(i4)                  :: nc
!
  data numbers/'0','1','2','3','4','5','6','7','8','9'/
!
!  Ensure that only a single node writes this file
!
  if (.not.ioproc) return
!
  do nc = 1,ncfg
    ncf = nc
!***************************
!  Initialisation of file  *
!***************************
    xsffile2 = ' '
    xsffilel = xsffile
!
!  If xsf file name has been given then open file
!
    if (xsffilel(1:1).eq.' ') then
      xsffilel = 'gulp.xsf'
    endif
    ind = index(xsffilel,'.xsf')
    if (ind.eq.0) then
      call endstring(xsffilel,len(xsffilel),ind)
      xsffilel(ind:ind+3) = '.xsf'
    endif
    xsffile2(1:ind-1) = xsffilel(1:ind-1)
    if (ncf.ge.100) then
      xsffile2(ind+4:ind+7) = xsffilel(ind:ind+3)
      xsffile2(ind:ind) = '_'
      incf = ncf
      ii = incf/10
      incf = incf - ii*10
      jj = incf/10
      incf = incf - jj*10
      xsffile2(ind+1:ind+1) = numbers(ii+1)
      xsffile2(ind+2:ind+2) = numbers(jj+1)
      xsffile2(ind+3:ind+3) = numbers(incf+1)
    elseif (ncf.ge.10) then
      xsffile2(ind+3:ind+6) = xsffilel(ind:ind+3)
      xsffile2(ind:ind) = '_'
      incf = ncf
      ii = incf/10
      incf = incf - ii*10
      xsffile2(ind+1:ind+1) = numbers(ii+1)
      xsffile2(ind+2:ind+2) = numbers(incf+1)
    elseif (ncfg.gt.1) then
      xsffile2(ind+2:ind+5) = xsffilel(ind:ind+3)
      xsffile2(ind:ind) = '_'
      xsffile2(ind+1:ind+1) = numbers(ncf+1)
    else
      xsffile2(ind:ind+3) = xsffilel(ind:ind+3)
    endif
    open(iout,file=xsffile2,form='formatted',status='unknown')
    if (ncfg.gt.1) then
      write(ioout,'(''  XSF file written for configuration '',i4,'' as '',a)') ncf,trim(xsffile2)
    else
      write(ioout,'(''  XSF file written as '',a)') trim(xsffile2)
    endif
!**********************************************************************
!  Write out energy if available
!**********************************************************************
    if (abs(energycfg(nc)).gt.1.0d-8) then
      if (ndim.eq.3) then
!
!  Correct energy to primitive value
!
        if (index(hmssg(1,1),'R').ne.0) then
          write(iout,'(''# total energy = '',f15.8,'' eV'',/)') energycfg(nc)/3.0_dp
        elseif (index(hmssg(1,1),'F').ne.0) then
          write(iout,'(''# total energy = '',f15.8,'' eV'',/)') energycfg(nc)/4.0_dp
        elseif (index(hmssg(1,1),'A').ne.0) then
          write(iout,'(''# total energy = '',f15.8,'' eV'',/)') energycfg(nc)/2.0_dp
        elseif (index(hmssg(1,1),'B').ne.0) then
          write(iout,'(''# total energy = '',f15.8,'' eV'',/)') energycfg(nc)/2.0_dp
        elseif (index(hmssg(1,1),'C').ne.0) then
          write(iout,'(''# total energy = '',f15.8,'' eV'',/)') energycfg(nc)/2.0_dp
        elseif (index(hmssg(1,1),'I').ne.0) then
          write(iout,'(''# total energy = '',f15.8,'' eV'',/)') energycfg(nc)/2.0_dp
        else
          write(iout,'(''# total energy = '',f15.8,'' eV'',/)') energycfg(nc)
        endif
      else
        write(iout,'(''# total energy = '',f15.8,'' eV'',/)') energycfg(nc)
      endif
    endif
    call setup(.false.)
!**********************************************************************
!  Write out keyword for dimensionality (and cell) for periodic systems
!**********************************************************************
    if (ndimen(ncf).eq.3) then
      write(iout,'(''CRYSTAL'',/)')
      write(iout,'(''PRIMVEC'')')
      write(iout,'(3f16.10)')(rvcfg(j,1,ncf),j=1,3)
      write(iout,'(3f16.10)')(rvcfg(j,2,ncf),j=1,3)
      write(iout,'(3f16.10)')(rvcfg(j,3,ncf),j=1,3)
    else
      write(iout,'(''MOLECULE'',/)')
    endif
!**********************************************************************
!  Write out coordinates
!**********************************************************************
    if (ndimen(ncf).eq.3) then
      write(iout,'(''PRIMCOORD'')')
      write(iout,'(i6,'' 1'')') ncore
      do i = 1,ncore
        ii = ncoptr(i)
        if (lxsfsym) then
          write(iout,'(a2,1x,3f12.6)') atsym(nat(ii)),xclat(ii),yclat(ii),zclat(ii)
        else
          write(iout,'(i3,3f12.6)') nat(ii),xclat(ii),yclat(ii),zclat(ii)
        endif
      enddo
    else
      write(iout,'(''ATOMS'')')
      do i = 1,ncore
        ii = ncoptr(i)
        if (lxsfsym) then
          write(iout,'(a2,1x,3f12.6)') atsym(nat(ii)),xclat(ii),yclat(ii),zclat(ii)
        else
          write(iout,'(i3,3f12.6)') nat(ii),xclat(ii),yclat(ii),zclat(ii)
        endif
      enddo
    endif
    close(iout)
!******************************
!  End of configuration loop  *
!******************************
  enddo
!
  return
  end
!
  subroutine outxsfseq(iout,nseq,lincforce,fc)
!
!  Write out XCrySDen XSF file - sequential writring of file during a run
!
!   4/18 Created from outxsf
!   6/18 Parallel I/O trap added
!   1/19 maxwordlength changes added
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
!  Copyright Curtin University 2019
!
!  Julian Gale, CIC, Curtin University, June 2019
!
  use control
  use configurations
  use current
  use derivatives,      only : xdrv, ydrv, zdrv
  use element
  use gulp_files
  use gulp_lengths
  use general
  use iochannels
  use parallel,         only : ioproc
  use shells
  use symmetry
  implicit none
!
!  Passed variables
!
  integer(i4),      intent(in)   :: iout        ! Output channel number to use
  integer(i4),      intent(in)   :: nseq        ! Sequence number for file
  logical,          intent(in)   :: lincforce   ! If true print forces after coordinates
  real(dp),         intent(in)   :: fc          ! Energy for structure
!
!  Local variables
!
  character(len=6)               :: seqnum
  character(len=maxwordlength)   :: xsffilel
  integer(i4)                    :: i
  integer(i4)                    :: ii
  integer(i4)                    :: ind
  integer(i4)                    :: j
!
!  Ensure that only a single node writes this file
!
  if (.not.ioproc) return
!***************************
!  Initialisation of file  *
!***************************
  xsffilel = ' '
  xsffilel(1:maxwordlength) = xsffile(1:maxwordlength)
!
!  Set file name as root + number (<999,999) + .xsf
!
  if (nseq.lt.999999_i4) then
    write(seqnum,'(i6)') nseq
  else
    call outerror('sequence number too large for XSF file name',0_i4)
    call stopnow('outxsfseq')
  endif
  do i = 1,6
    if (seqnum(i:i).eq.' ') seqnum(i:i) = '0'
  enddo
  ind = index(xsffilel,'.xsf')
  if (ind.eq.0) then
    ind = index(xsffilel,' ')
  endif
  xsffilel(ind:ind+5) = seqnum
  ind = ind + 6
  xsffilel(ind:ind+3) = '.xsf'
!
!  Open file
!
  open(iout,file=xsffilel,form='formatted',status='unknown')
!**********************************************************************
!  Write out energy 
!**********************************************************************
  write(iout,'(''# total energy = '',f15.8,'' eV'',/)') fc
!**********************************************************************
!  Write out keyword for dimensionality (and cell) for periodic systems
!**********************************************************************
  if (ndimen(ncf).eq.3) then
    write(iout,'(''CRYSTAL'',/)')
    write(iout,'(''PRIMVEC'')')
    write(iout,'(3f16.10)')(rv(j,1),j=1,3)
    write(iout,'(3f16.10)')(rv(j,2),j=1,3)
    write(iout,'(3f16.10)')(rv(j,3),j=1,3)
  else
    write(iout,'(''MOLECULE'',/)')
  endif
!**********************************************************************
!  Write out coordinates
!**********************************************************************
  if (ndim.eq.3) then
    write(iout,'(''PRIMCOORD'')')
    write(iout,'(i6,'' 1'')') ncore
    do i = 1,ncore
      ii = ncoptr(i)
      if (lincforce) then
        if (lxsfsym) then
          write(iout,'(a2,1x,6f12.6)') atsym(nat(ii)),xclat(ii),yclat(ii),zclat(ii), &
                                            xdrv(ii),ydrv(ii),zdrv(ii)
        else
          write(iout,'(i3,6f12.6)') nat(ii),xclat(ii),yclat(ii),zclat(ii), &
                                            xdrv(ii),ydrv(ii),zdrv(ii)
        endif
      else
        if (lxsfsym) then
          write(iout,'(a2,1x,3f12.6)') atsym(nat(ii)),xclat(ii),yclat(ii),zclat(ii)
        else
          write(iout,'(i3,3f12.6)') nat(ii),xclat(ii),yclat(ii),zclat(ii)
        endif
      endif
    enddo
  else
    write(iout,'(''ATOMS'')')
    do i = 1,ncore
      ii = ncoptr(i)
      if (lincforce) then
        if (lxsfsym) then
          write(iout,'(a2,1x,6f12.6)') atsym(nat(ii)),xclat(ii),yclat(ii),zclat(ii), &
                                            xdrv(ii),ydrv(ii),zdrv(ii)
        else
          write(iout,'(i3,6f12.6)') nat(ii),xclat(ii),yclat(ii),zclat(ii), &
                                            xdrv(ii),ydrv(ii),zdrv(ii)
        endif
      else
        if (lxsfsym) then
          write(iout,'(a2,1x,3f12.6)') atsym(nat(ii)),xclat(ii),yclat(ii),zclat(ii)
        else
          write(iout,'(i3,3f12.6)') nat(ii),xclat(ii),yclat(ii),zclat(ii)
        endif
      endif
    enddo
  endif
!
!  Close file
!
  close(iout)
!
  return
  end
