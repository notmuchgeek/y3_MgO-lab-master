  subroutine phonword(nru,word,lwordok,iline,line,l55,l1000,ncurr)
!
!  Processes input for phonon related words
!
!  nru = fortran channel for reading input
!
!  12/00 Generalised for 1/2-D
!   6/01 Initialisation of line added for benefit of some compilers
!   3/02 Input of Born approach direction added
!   3/02 Input for frequency dependent properties added
!   6/05 outwarning called for warnings
!   2/06 Setting of omegadirtype added
!   3/06 Omega damping factor added
!   8/06 nru passed to linepro
!   1/08 Code style updated
!  12/08 Module input renamed to gulpinput
!   3/09 lkptdispersion flag added
!   6/09 lshrinkset(ncurr) & ldispersionset(ncurr) added for the 
!        PDF functionality (ers29)
!   6/13 Lorentzian drop tolerance added
!   7/13 lbroaden_scale added
!   8/13 omega_af added
!   9/13 v_p_cfg added
!  10/13 lbroaden_gauss added
!   4/14 lbornkin added
!  10/14 fc_supercell option added
!   1/15 Error in memory allocation for projected DOS fixed
!   2/18 Alamode options added
!   2/18 nxks, nyks, nzks converted to a single array
!   1/19 maxwordlength changes added
!   3/19 Multiple temperature ramps added
!   3/19 Correction to temperature unit handling for temperature steps
!   2/20 lksorigin added
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
!  Julian Gale, CIC, Curtin University, February 2020
!
  use configurations
  use g_constants,    only : pi, speedl
  use control
  use dispersion
  use general,        only : bfactor, nwarn, Lor_tol
  use gulpinput
  use iochannels
  use ksample
  use m_alamode,      only : cutoff_fc2, cutoff_fc3, nalaprocs
  use m_alamode,      only : disp_fc2, disp_fc3, Bohr_radius, ala_mode
  use m_pdfneutron,   only : lshrinkset, ldispersionset
  use parallel
  use phonout
  use projectdos
  use thermalcond,    only : omega_af, lomega_af_in, B_pr_cfg, v_s_cfg, v_p_cfg
  implicit none
!
!  Passed variables
!
  character(len=maxwordlength) :: word
  character(len=maxlinelength) :: line
  integer(i4)                  :: iline
  integer(i4)                  :: ncurr
  integer(i4)                  :: nru
  logical                      :: lfound
  logical                      :: lwordok
  logical                      :: l55
  logical                      :: l1000
!
!  Local variables
!
  integer(i4)                  :: i
  integer(i4)                  :: ii
  integer(i4)                  :: ifptr
  integer(i4)                  :: inat
  integer(i4)                  :: ind
  integer(i4)                  :: itype
  integer(i4)                  :: j
  integer(i4)                  :: jj
  integer(i4)                  :: kk
  integer(i4)                  :: nksin
  integer(i4)                  :: nn
  integer(i4)                  :: npdb
  integer(i4), save            :: nproj = 0
  integer(i4), save            :: nproji = 0
  integer(i4)                  :: nprojiold
  integer(i4)                  :: npt
  integer(i4)                  :: nread
  logical                      :: lmultiT
  logical                      :: lnone
  logical                      :: lTinC
  logical                      :: lTinF
  real(dp)                     :: tem
  real(dp)                     :: temstp
  real(dp)                     :: units
  real(dp)                     :: units2
!
  if (index(word,'fc_s').eq.1) goto 600
  if (index(word,'odir').eq.1) goto 610
  if (index(word,'omega ').eq.1) goto 620
  if (index(word,'gamma_d').eq.1) goto 630
  if (index(word,'gamma_a').eq.1) goto 635
  if (index(word,'kpoi').eq.1) goto 640
  if (index(word,'disp').eq.1) goto 650
  if (index(word,'shri').eq.1) goto 660
  if (index(word,'box ').eq.1) goto 670
  if (index(word,'temp').eq.1) goto 680
  if (index(word,'eige').eq.1) goto 690
  if (index(word,'lowe').eq.1) goto 700
  if (index(word,'proj').eq.1) goto 710
  if (index(word,'omega_d').eq.1) goto 720
  if (index(word,'broa').eq.1) goto 730
  if (index(word,'lore').eq.1) goto 740
  if (index(word,'omega_a').eq.1) goto 750
  if (index(word,'rdir').eq.1) goto 760
  if (index(word,'ala_c').eq.1) goto 770
  if (index(word,'ala_d').eq.1) goto 780
  if (index(word,'ala_s').eq.1) goto 790
  if (index(word,'ala_p').eq.1) goto 800
  if (index(word,'ala_m').eq.1) goto 810
  return
!************************
!  FC_supercell option  *
!************************
600 if (ndimen(ncurr).eq.2) then
    if (nfloat.ge.2) then
      ii = nint(abs(floats(1)))
      jj = nint(abs(floats(2)))
      kk = 0
    else
      line = '  '
      read(nru,'(a)')line
      iline = iline + 1
      call linepro(nru,line,iline)
      if (nfloat.ge.2) then
        ii = nint(abs(floats(1)))
        jj = nint(abs(floats(2)))
        kk = 0
      else
        call outerror('Insufficient data for fc_supercell input',iline)
        call stopnow('phonword')
      endif
    endif
  elseif (ndimen(ncurr).eq.1) then
    if (nfloat.ge.1) then
      ii = nint(abs(floats(1)))
      jj = 0
      kk = 0
    else
      line = '  '
      read(nru,'(a)') line
      iline = iline + 1
      call linepro(nru,line,iline)
      if (nfloat.ge.1) then
        ii = nint(abs(floats(1)))
        jj = 0
        kk = 0
      else
        call outerror('Insufficient data for fc_supercell input',iline)
        call stopnow('phonword')
      endif
    endif
  else
    if (nfloat.ge.3) then
      ii = nint(abs(floats(1)))
      jj = nint(abs(floats(2)))
      kk = nint(abs(floats(3)))
    else
      line = '  '
      read(nru,'(a)') line
      iline = iline + 1
      call linepro(nru,line,iline)
      if (nfloat.ge.3) then
        ii = nint(abs(floats(1)))
        jj = nint(abs(floats(2)))
        kk = nint(abs(floats(3)))
      else
        call outerror('Insufficient data for fc_supercell input',iline)
        call stopnow('phonword')
      endif
    endif
  endif
  nd2cellcfg(1,ncurr) = ii
  nd2cellcfg(2,ncurr) = jj
  nd2cellcfg(3,ncurr) = kk
  lwordok = .true.
  return
!************************************************************
!  Directions for frequency dependent property observation  *
!************************************************************
610 omegadirtype(ncurr) = 1
  if (nword.gt.1) then
    do i = 2,nword
      if (index(words(i),'frac').eq.1) omegadirtype(ncurr) = 2
    enddo
  endif
  if (nfloat.ge.6) then
    omegadir(1:6,ncurr) = floats(1:6)
  else
    line = '  '                
    read(nru,'(a)') line       
    iline = iline + 1          
    call linepro(nru,line,iline)
    if (nfloat.ge.6) then
      omegadir(1:6,ncurr) = floats(1:6)
    else
      call outerror('Insufficient data for odirection option',iline)
      call stopnow('phonword') 
    endif
  endif
  lwordok = .true.
  return
!************************************************
!  Frequencies for phonon property calculation  *
!************************************************
620 if (nfloat.ge.1) then
    omega(ncurr) = floats(1)
    if (nfloat.ge.3) then
!
!  Omega range specified
!
      omegastep(ncurr)=floats(2)
      nomegastep(ncurr)=nint(floats(3))
      if ((omega(ncurr)+omegastep(ncurr)*nomegastep(ncurr)).lt.0.0_dp) then
        call outerror('Omega step is such that omega will go negative',iline)
        call stopnow('phonword')
      endif
    else
      nomegastep(ncurr) = 0
      omegastep(ncurr) = 0.0_dp
    endif
  else
    line = '  '
    read(nru,'(a)') line
    iline = iline + 1
    call linepro(nru,line,iline)
    if (nfloat.ge.1) then
      omega(ncurr) = floats(1)
      if (nfloat.ge.3) then
!
!  Omega range specified
!
        omegastep(ncurr) = floats(2)
        nomegastep(ncurr) = nint(floats(3))
        if ((omega(ncurr)+omegastep(ncurr)*nomegastep(ncurr)).lt.0.0_dp) then
          call outerror('Omega step is such that omega will go negative',iline)
          call stopnow('phonword')
        endif
      else
        nomegastep(ncurr) = 0
        omegastep(ncurr) = 0.0_dp
      endif
    else
      call outerror('Omega is missing from input',iline)
      call stopnow('phonword')
    endif
  endif
!
!  Set flag for freqeuncy-dependent property calculation
!
  lomega(ncurr) = .true.
  lwordok = .true.
  return
!*****************************
!  Gamma approach direction  *
!*****************************
630 if (nfloat.lt.3) then
    call outerror('Insufficient data in gamma approach direction',iline)
    call stopnow('phonword')
  endif
  lbornkin(ncurr) = .true.
  bornk(1,ncurr) = floats(1)
  bornk(2,ncurr) = floats(2)
  bornk(3,ncurr) = floats(3)
  lwordok = .true.
  return
!***********************************
!  Gamma angular integration steps *
!***********************************
635 if (nfloat.eq.0) then
    call outerror('Insufficient data in gamma integration steps',iline)
    call stopnow('phonword')
  endif
  nbornstep(ncurr) = nint(abs(floats(1)))
  lwordok = .true.
  return
!******************************************
!  K points input for phonon calculation  *
!******************************************
640 if (nfloat.gt.0) then
    nread = nint(floats(1))
  else
    nread = 1
  endif
  if (ndimen(ncurr).eq.0) then
    call outerror('Kpoints option not valid for a cluster',iline)
    call stopnow('phonword')
  endif
  norigkpt(ncurr) = norigkpt(ncurr) + nread
  if (nkpt+nread.gt.maxkpt) then
    maxkpt = nkpt + nread
    call changemaxkpt
  endif
  do i = 1,nread
645 line = '  '
    read(nru,'(a)')line
    iline = iline + 1
    if (index(line,'#').eq.1) goto 645
    call linepro(nru,line,iline)
    nkpt = nkpt + 1
    nkptcfg(nkpt) = ncurr
    lkptdispersion(nkpt) = .false.
    if (nfloat.ge.4) then
      xkpt(nkpt) = floats(1)
      ykpt(nkpt) = floats(2)
      zkpt(nkpt) = floats(3)
      wkpt(nkpt) = floats(4)
    elseif (nfloat.eq.3) then
      xkpt(nkpt) = floats(1)
      ykpt(nkpt) = floats(2)
      zkpt(nkpt) = floats(3)
      wkpt(nkpt) = 1.0_dp
    elseif (nfloat.eq.2) then
      xkpt(nkpt) = floats(1)
      ykpt(nkpt) = floats(2)
      zkpt(nkpt) = 0.0_dp
      wkpt(nkpt) = 1.0_dp
    elseif (nfloat.eq.1) then
      xkpt(nkpt) =floats(1)
      ykpt(nkpt) = 0.0_dp
      zkpt(nkpt) = 0.0_dp
      wkpt(nkpt) = 1.0_dp
    else
      call outerror('Insufficient input in kpoints option',iline)
      call stopnow('phonword')
    endif
  enddo
  lwordok = .true.
  return
!*****************************
!  Phonon dispersion curves  *
!*****************************
650 if (nfloat.gt.0) then
    nread = nint(floats(1))
    if (nfloat.ge.2) ndispres = nint(floats(2))
  else
    nread = 1
  endif
  if (ndimen(ncurr).eq.0) then
    call outerror('Dispersion option not valid for a cluster',iline)
    call stopnow('phonword')
  endif
  ldispersionset(ncurr) = .true.
  do i = 1,nread
    ndline = ndline + 1
    ndstart(ndline) = ndpoint + 1
655 line = '  '
    read(nru,'(a)')line
    iline = iline + 1
    if (index(line,'#').eq.1) goto 655
    call linepro(nru,line,iline)
    if (ndimen(ncurr).eq.3) then
      npt = nfloat/3
      ind = 0
      do j = 1,npt
        ndpoint = ndpoint + 1
        if (ndpoint.gt.maxkpt) then
          maxkpt = ndpoint + 10
          call changemaxkpt
        endif
        xdisp(ndpoint) = floats(ind+1)
        ydisp(ndpoint) = floats(ind+2)
        zdisp(ndpoint) = floats(ind+3)
        ind = ind + 3
      enddo
    elseif (ndimen(ncurr).eq.2) then
      npt = nfloat/2
      ind = 0
      do j = 1,npt
        ndpoint = ndpoint + 1
        if (ndpoint.gt.maxkpt) then
          maxkpt = ndpoint + 10
          call changemaxkpt
        endif
        xdisp(ndpoint) = floats(ind+1)
        ydisp(ndpoint) = floats(ind+2)
        zdisp(ndpoint) = 0.0_dp
        ind = ind + 2
      enddo
    elseif (ndimen(ncurr).eq.1) then
      npt = nfloat
      ind = 0
      do j = 1,npt
        ndpoint = ndpoint + 1
        if (ndpoint.gt.maxkpt) then
          maxkpt = ndpoint + 10
          call changemaxkpt
        endif
        xdisp(ndpoint) = floats(ind+1)
        ydisp(ndpoint) = 0.0_dp
        zdisp(ndpoint) = 0.0_dp
        ind = ind + 1
      enddo
    endif
    ndend(ndline) = ndpoint
    ndispcfg(ndline) = ncurr
    if (ndpoint.gt.maxkpt) then
      maxkpt = ndpoint + 10
      call changemaxkpt
    endif
!
!  If there is a "to" statement at the end then continue reading lines
!
    if (nword.eq.nfloat) goto 655
    if ((ndend(ndline)-ndstart(ndline)).lt.1) then
      call outerror('Dispersion input has only one point',iline)
      call stopnow('phonword')
    endif
  enddo
  lwordok = .true.
  return
!*****************************************
!  Shrinking factors for Brillouin zone  *
!*****************************************
660 if (ndimen(ncurr).eq.0) then
    call outerror('Shrink option not valid for a cluster',iline)
    call stopnow('phonword')
  endif
  if (nword.gt.1) then
    do i = 2,nword
      if (index(words(i),'orig').eq.1) lksorigin(ncurr) = .true.
    enddo
  endif
  lshrinkset(ncurr) = .true.
  if (nfloat.gt.0) then
    if (nfloat.eq.1) then
      nksin = nint(floats(1))
      nks(1,ncurr) = nksin
      nks(2,ncurr) = nksin
      nks(3,ncurr) = nksin
    elseif (nfloat.eq.2) then
      nksin = nint(floats(1))
      nks(1,ncurr) = nksin
      nksin = nint(floats(2))
      nks(2,ncurr) = nksin
      nks(3,ncurr) = nksin
    else
      nks(1,ncurr) = nint(floats(1))
      nks(2,ncurr) = nint(floats(2))
      nks(3,ncurr) = nint(floats(3))
    endif
  else
665 line  =  '  '
    read(nru,'(a)') line
    iline = iline + 1
    if (index(line,'#').eq.1) goto 665
    call linepro(nru,line,iline)
    if (nfloat.eq.1) then
      nksin = nint(floats(1))
      nks(1,ncurr) = nksin
      nks(2,ncurr) = nksin
      nks(3,ncurr) = nksin
    elseif (nfloat.eq.2) then
      nksin = nint(floats(1))
      nks(1,ncurr) = nksin
      nks(2,ncurr) = nksin
      nks(3,ncurr) = nint(floats(2))
    elseif (nfloat.ge.3) then
      nks(1,ncurr) = nint(floats(1))
      nks(2,ncurr) = nint(floats(2))
      nks(3,ncurr) = nint(floats(3))
    else
      call outerror('No shrinking factors given in shrink option',iline)
      call stopnow('phonword')
    endif
  endif
  nks(1,ncurr) = abs(nks(1,ncurr))
  nks(2,ncurr) = abs(nks(2,ncurr))
  nks(3,ncurr) = abs(nks(3,ncurr))
!
!  Handle 1-D and 2-D case
!
  if (ndimen(ncurr).eq.2) then
    nks(3,ncurr) = 1
  elseif (ndimen(ncurr).eq.1) then
    nks(2,ncurr) = 1
    nks(3,ncurr) = 1
  endif
  lwordok = .true.
  return
!*************************************
!  Box dimensions for phonon output  *
!*************************************
670 if (nword.eq.1) then
675 line = '  '
    read(nru,'(a)') line
    iline = iline + 1
    if (index(line,'#').eq.1) goto 675
    call linepro(nru,line,iline)
  endif
  ifptr = 1
  do i = 1,nword
    if (index(words(i),'disp').eq.1) then
      j = i + 1
      lfound = .false.
      do while (j.le.nword.and..not.lfound)
        if (index(words(j),'size').eq.1) then
          flbox = floats(ifptr)
          ifptr = ifptr + 1
          lfound = .true.
        elseif (index(words(j),'numb').eq.1) then
          nlbox = nint(floats(ifptr))
          ifptr = ifptr + 1
          lfound = .true.
        endif
        j = j + 1
      enddo
      if (.not.lfound) then
        nwarn = nwarn + 1
        call outwarning('type missing in box specification',0_i4)
      endif
    elseif (index(words(i),'dens').eq.1) then
      j = i + 1
      lfound = .false.
      do while (j.le.nword.and..not.lfound)
        if (index(words(j),'size').eq.1) then
          fbox = floats(ifptr)
          ifptr = ifptr + 1
          lfound = .true.
        elseif (index(words(j),'numb').eq.1) then
          nbox = nint(floats(ifptr))
          ifptr = ifptr + 1
          lfound = .true.
        endif
        j = j + 1
      enddo
      if (.not.lfound) then
        nwarn = nwarn + 1
        call outwarning('type missing in box specification',0_i4)
      endif
    endif
  enddo
  lwordok = .true.
  return
!*********************************
!  Temperature of configuration  *
!*********************************
680 continue
  lmultiT = .false.
  lTinC = .false.
  lTinF = .false.
  if (nword.gt.1) then
    do i = 2,nword
      if (index(words(i),'mult').eq.1) lmultiT = .true.
      if (index(words(i),'c').eq.1) lTinC = .true.
      if (index(words(i),'f').eq.1) lTinC = .true.
    enddo
  endif
  if (nfloat.ge.1) then
    tem = floats(1)
    if (nfloat.ge.3) then
!
!  T range specified
!
      ntempramp(ncurr) = 1
      temstp = floats(2)
      ntempstp(1,ncurr) = nint(floats(3))
      if (nfloat.ge.4) then
        ntempstpstart(1,ncurr) = nint(floats(4))
      endif
      if ((tem+temstp*ntempstp(1,ncurr)).lt.0.0_dp) then
        call outerror('Temperature step is such that T will go negative',iline)
        call stopnow('phonword')
      endif
    else
      temstp = 0.0_dp
    endif
!
!  Convert units
!
    if (lTinC) then
      tem = tem + 273.0_dp
    elseif (lTinF) then
      tem = 273.0_dp + 5.0_dp*(tem - 32.0_dp)/9.0_dp
      temstp = 5.0_dp*(temstp - 32.0_dp)/9.0_dp
    endif
  else
    line = '  '
    read(nru,'(a)') line
    iline = iline + 1
    call linepro(nru,line,iline)
    if (nfloat.ge.1) then
      tem = floats(1)
      if (nfloat.ge.3) then
!
!  T range specified
!
        ntempramp(ncurr) = 1
        temstp = floats(2)
        ntempstp(1,ncurr) = nint(floats(3))
        if (nfloat.ge.4) then
          ntempstpstart(1,ncurr) = nint(floats(4))
        endif
        if ((tem+temstp*ntempstp(1,ncurr)).lt.0.0_dp) then
          call outerror('Temperature step is such that T will go negative',iline)
          call stopnow('phonword')
        endif
      else
        temstp = 0.0_dp
      endif
!
!  Convert units
!
      if (lTinC) then
        tem = tem + 273.0_dp
      elseif (lTinF) then
        tem = 273.0_dp + 5.0_dp*(tem - 32.0_dp)/9.0_dp
        temstp = 5.0_dp*(temstp - 32.0_dp)/9.0_dp
      endif
    else
      call outerror('Temperature is missing from input',iline)
      call stopnow('phonword')
    endif
  endif
  tempcfg(ncurr) = tem
  tempstp(1,ncurr) = temstp
!
  if (lmultiT) then
!
!  Multiple temperature ramps
!
685 line = ' '
    read(nru,'(a)',end=688) line
    iline = iline + 1
    call linepro(nru,line,iline)
!
!  Check whether line is empty
!
    if ((nword+nfloat).eq.0) goto 685
!
!  If there are words on the line then this is a new option
!
    if (nword.gt.0) then
      l55 = .true.
      goto 688
    endif
!
!  Increment running temperature for start of this ramp
!
    tem = tem + temstp*ntempstp(ntempramp(ncurr),ncurr)
!
    if (nfloat.ge.2) then
!
!  Add a new T range
!
      ntempramp(ncurr) = ntempramp(ncurr) + 1
      if (ntempramp(ncurr).gt.maxtempramp) then
        maxtempramp = ntempramp(ncurr) + 1
        call changemaxtempramp
      endif
      temstp = floats(1)
      ntempstp(ntempramp(ncurr),ncurr) = nint(floats(2))
      if (nfloat.ge.3) then
        ntempstpstart(ntempramp(ncurr),ncurr) = nint(floats(3))
      else
        ntempstpstart(ntempramp(ncurr),ncurr) = 0
      endif
!
!  Unit handling
!
      if (lTinF) then
        temstp = 5.0_dp*(temstp - 32.0_dp)/9.0_dp
      endif
!
      if ((tem+temstp*ntempstp(ntempramp(ncurr),ncurr)).lt.0.0_dp) then
        call outerror('Temperature step is such that T will go negative',iline)
        call stopnow('phonword')
      endif
      if (nword.ge.2) then
        call stolc(words(2),maxword)
        if (index(words(2),'c ').eq.1) then
          tem = tem + 273.0_dp
          temstp = temstp + 273.0_dp
        elseif (index(words(2),'f ').eq.1) then
          tem = 273.0_dp + 5.0_dp*(tem - 32.0_dp)/9.0_dp
          temstp = 273.0_dp + 5.0_dp*(temstp - 32.0_dp)/9.0_dp
        endif
      endif
      tempstp(ntempramp(ncurr),ncurr) = temstp
    else
      call outerror('Insufficient input for temperature ramp',iline)
      call stopnow('phonword')
    endif
    goto 685
688 if (.not.l55) l1000 = .true.
  endif
  lwordok = .true.
  return
!**********************
!  Eigenvector range  *
!**********************
690 if (nfloat.eq.0) then
    call outerror('No eigenvector range specified in option',iline)
    call stopnow('phonword')
  endif
  neiglow(ncurr) = abs(nint(floats(1)))
  if (nfloat.ge.2) then
    neighigh(ncurr) = abs(nint(floats(2)))
    if (neiglow(ncurr).gt.neighigh(ncurr)) then
!
!  Switch values round as high is lower than low
!
      nn = neiglow(ncurr)
      neiglow(ncurr) = neighigh(ncurr)
      neighigh(ncurr) = nn
    endif
  endif
  lwordok = .true.
  return
!***************************************
!  Lowest phonon mode for free energy  *
!***************************************
700 if (nfloat.gt.0) then
    minmodecfg(ncurr) = nint(floats(1))
    if (nfloat.gt.1) maxmodecfg(ncurr) = nint(floats(2))
  endif
  if (minmodecfg(ncurr).le.0) minmodecfg(ncurr) = 1
  lwordok = .true.
  return
!************************************************
!  Project phonon density of states onto atoms  *
!************************************************
710 if (nfloat.ge.1) then
    nread = nint(floats(1))
  else
    read(nru,*) nread
    iline = iline + 1
  endif
  nprojcfg(ncfg) = nprojcfg(ncfg) + nread
  leigen = .true.
!
!  Set logical according to whether bulk or defect projection
!
  if (nword.gt.1.and.index(words(2),'def').eq.1) then
    npdb = 2
    nprojdef(ncfg) = nprojdef(ncfg) + nread
  else
    npdb = 1
  endif
  do i = 1,nread
    line = '  '
    read(nru,'(a)') line
    iline = iline + 1
    call linepro(nru,line,iline)
    nprojiold = nproji
    nproj = nproj + 1
    if (nproj.gt.maxproj) then
      maxproj = nproj
      call changemaxproj
    endif
    nprojdb(nproj) = npdb
!
!  Check bounds
!
    if (nproji+nfloat+nword.gt.maxproji) then
      maxproji = nproji + nfloat + nword + 10
      call changemaxproji
    endif
!
!  Do explicit atom numbers
!
    do j = 1,nfloat
      nproji = nproji + 1
      nprojnat(nproji) = nint(floats(j))
      nprojtyp(nproji) = 100
      nprojptr(nproji) = nproj
    enddo
!
!  Do species
!
    do j = 1,nword
      call ltont(words(j),inat,itype)
      nproji = nproji + 1
      nprojnat(nproji) = inat
      nprojtyp(nproji) = itype
      nprojptr(nproji) = nproj
    enddo
    nprojit(nproj) = nproji - nprojiold
  enddo
  lwordok = .true.
  return
!*******************************************
!  Broadening factor for Omega properties  *
!*******************************************
720 if (nfloat.ge.1) then
    omegadamping(ncurr) = floats(1)
  else
    read(nru,*) omegadamping(ncurr)
    iline = iline + 1
  endif
  omegadamping(ncurr) = max(omegadamping(ncurr),1.0d-6)
  lwordok = .true.
  return
!******************************
!  Broadening factor for DOS  *
!******************************
730 continue
  lbroaden_gauss = .false.
  lbroaden_scale = .false.
  if (nword.gt.1) then
    lbroaden_gauss = (index(words(2),'gau').eq.1) 
    lbroaden_scale = (index(words(2),'sca').eq.1) 
  endif
  if (nfloat.ge.1) then
    bfactor = abs(floats(1))
  else
    call outerror('Broadening factor is missing from input',iline)
    call stopnow('phonword')
  endif
  lbroad = .true.
  lwordok = .true.
  return
!*******************************************************
!  Lorentzian drop tolerance for thermal conductivity  *
!*******************************************************
740 if (nfloat.ge.1) then
    Lor_tol = abs(floats(1))
  else
    read(nru,*) Lor_tol
    iline = iline + 1
  endif
  lwordok = .true.
  return
!*****************************************************************************
!  Minimum frequency for Allen-Feldman contribution to thermal conductivity  *
!*****************************************************************************
750 continue
  if (nfloat.ge.1) then
    lomega_af_in(ncurr) = .true.
    units  = 1.0_dp
    units2 = 1.0d-6
    if (nword.gt.1) then
      if (index(words(2),'rad').eq.1) then
        units = 1.0d12/(2.0_dp*pi*speedl)
        units2 = 1.0d18/(2.0_dp*pi*speedl)**2
      endif
    endif
    if (nfloat.ge.4) then
      omega_af(ncurr) = floats(1)*units
      B_pr_cfg(ncurr) = abs(floats(2))*units2
      v_s_cfg(ncurr) = abs(floats(3))/1000.0_dp
      v_p_cfg(ncurr) = abs(floats(4))/1000.0_dp
    elseif (nfloat.ge.3) then
      omega_af(ncurr) = floats(1)*units
      B_pr_cfg(ncurr) = abs(floats(2))*units2
      v_s_cfg(ncurr) = abs(floats(3))/1000.0_dp
      v_p_cfg(ncurr) = v_s_cfg(ncurr)
    elseif (nfloat.ge.2) then
      omega_af(ncurr) = floats(1)*units
      B_pr_cfg(ncurr) = abs(floats(2))*units2
    else
      omega_af(ncurr) = floats(1)*units
    endif
  else
    call outerror('Omega_af value is missing from input',iline)
    call stopnow('phonword')
  endif
  lwordok = .true.
  return
!***********************************
!  Directions for Raman intensity  *
!***********************************
760 ramandirtype(ncurr) = 1
  if (nword.gt.1) then
    do i = 2,nword
      if (index(words(i),'frac').eq.1) ramandirtype(ncurr) = 2
    enddo
    if (ndimen(ncurr).eq.0) then
      call outerror('Fractional rdirection specified for molecule',iline)
      call stopnow('phonword')
    endif
  endif
  if (nfloat.ge.6) then
    ramandir(1:6,ncurr) = floats(1:6)
  else
    line = '  '
    read(nru,'(a)') line
    iline = iline + 1
    call linepro(nru,line,iline)
    if (nfloat.ge.6) then
      ramandir(1:6,ncurr) = floats(1:6)
    else
      call outerror('Insufficient data for rdirection option',iline)
      call stopnow('phonword')
    endif
  endif
!
!  Check that directions are normalisable
!
  if (abs(ramandir(1,ncurr))+abs(ramandir(2,ncurr))+abs(ramandir(3,ncurr)).eq.0.0_dp) then
    call outerror('In rdirection is equal to zero',iline)
    call stopnow('phonword')
  endif
  if (abs(ramandir(4,ncurr))+abs(ramandir(5,ncurr))+abs(ramandir(6,ncurr)).eq.0.0_dp) then
    call outerror('Out rdirection is equal to zero',iline)
    call stopnow('phonword')
  endif
  lwordok = .true.
  return
!****************************************
!  Cutoffs for Alamode force constants  *
!****************************************
770 units = 1.0_dp
  lnone = .false.
  if (nword.ge.2) then
    call stolc(words(2),maxword)
    if (index(words(2),'none').ne.0) then
      cutoff_fc2 = -1.0_dp
      cutoff_fc3 = -1.0_dp
      lnone = .true.
    elseif (index(words(2),'au').ne.0) then
      units = Bohr_radius
    endif
  endif
  if (.not.lnone) then
    if (nfloat.ge.2) then
      cutoff_fc2 = floats(1)*units
      cutoff_fc3 = floats(2)*units
    elseif (nfloat.ge.1) then
      cutoff_fc2 = floats(1)*units
      cutoff_fc3 = floats(1)*units
    else
      call outerror('cutoff for Alamode force constants is missing',iline)
      call stopnow('phonword')
    endif
  endif
  lwordok = .true.
  return
!**********************************************
!  Displacements for Alamode force constants  *
!**********************************************
780 units = 1.0_dp
  if (nword.ge.2) then
    call stolc(words(2),maxword)
    if (index(words(2),'au').ne.0) then
      units = Bohr_radius
    endif
  endif
  if (nfloat.ge.2) then
    disp_fc2 = abs(floats(1))*units
    disp_fc3 = abs(floats(2))*units
  elseif (nfloat.ge.1) then
    disp_fc2 = abs(floats(1))*units
    disp_fc3 = abs(floats(1))*units
  else
    call outerror('displacement for Alamode force constants is missing',iline)
    call stopnow('phonword')
  endif
!
!  Check magnitude of displacements
!
  if (disp_fc2.gt.0.1_dp.or.disp_fc3.gt.0.1_dp) then
    nwarn = nwarn + 1
    call outwarning('displacements for Alamode seem too large',0_i4)
  endif
  if (disp_fc2.lt.1.0d-5.or.disp_fc3.lt.1.0d-5) then
    nwarn = nwarn + 1
    call outwarning('displacements for Alamode seem too small',0_i4)
  endif
  lwordok = .true.
  return
!****************************************************
!  Shrinking factors for Brillouin zone in Alamode  *
!****************************************************
790 if (ndimen(ncurr).eq.0) then
    call outerror('Shrink option not valid for a cluster',iline)
    call stopnow('phonword')
  endif
  if (nfloat.gt.0) then
    if (nfloat.eq.1) then
      nksin = nint(floats(1))
      nksala(1,ncurr) = nksin
      nksala(2,ncurr) = nksin
      nksala(3,ncurr) = nksin
    elseif (nfloat.eq.2) then
      nksin = nint(floats(1))
      nksala(1,ncurr) = nksin
      nksin = nint(floats(2))
      nksala(2,ncurr) = nksin
      nksala(3,ncurr) = nksin
    else
      nksala(1,ncurr) = nint(floats(1))
      nksala(2,ncurr) = nint(floats(2))
      nksala(3,ncurr) = nint(floats(3))
    endif
  else
795 line  =  '  '
    read(nru,'(a)') line
    iline = iline + 1
    if (index(line,'#').eq.1) goto 795
    call linepro(nru,line,iline)
    if (nfloat.eq.1) then
      nksin = nint(floats(1))
      nksala(1,ncurr) = nksin
      nksala(2,ncurr) = nksin
      nksala(3,ncurr) = nksin
    elseif (nfloat.eq.2) then
      nksin = nint(floats(1))
      nksala(1,ncurr) = nksin
      nksala(2,ncurr) = nksin
      nksala(3,ncurr) = nint(floats(2))
    elseif (nfloat.ge.3) then
      nksala(1,ncurr) = nint(floats(1))
      nksala(2,ncurr) = nint(floats(2))
      nksala(3,ncurr) = nint(floats(3))
    else
      call outerror('No shrinking factors given in shrink option',iline)
      call stopnow('phonword')
    endif
  endif
  nksala(1,ncurr) = abs(nksala(1,ncurr))
  nksala(2,ncurr) = abs(nksala(2,ncurr))
  nksala(3,ncurr) = abs(nksala(3,ncurr))
!
!  Handle 1-D and 2-D case
!
  if (ndimen(ncurr).eq.2) then
    nksala(3,ncurr) = 1
  elseif (ndimen(ncurr).eq.1) then
    nksala(2,ncurr) = 1
    nksala(3,ncurr) = 1
  endif
  lwordok = .true.
  return
!*************************************
!  Number of processors for Alamode  *
!*************************************
800 continue
  if (nfloat.ge.1) then
    nalaprocs = abs(nint(floats(1)))
  else
    call outerror('number of processors is missing for Alamode',iline)
    call stopnow('phonword')
  endif
  lwordok = .true.
  return
!*****************************
!  Mode control for Alamode  *
!*****************************
810 continue
  if (nword.gt.1) then
    ala_mode = words(2)
    if (index(ala_mode,'all').ne.1.and. &
        index(ala_mode,'pat').ne.1.and. &
        index(ala_mode,'dis').ne.1.and. &
        index(ala_mode,'fit').ne.1.and. &
        index(ala_mode,'bte').ne.1) then
      call outerror('mode type is unknown for Alamode',iline)
      call stopnow('phonword')
    endif
  else
    call outerror('mode type is missing for Alamode',iline)
    call stopnow('phonword')
  endif
  lwordok = .true.
  return
!
  end
