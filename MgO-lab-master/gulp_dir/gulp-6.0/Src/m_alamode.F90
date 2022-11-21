  module m_alamode
    use datatypes
    use gulp_lengths
#ifdef NAG
!For NAG compiler only:
    use f90_unix_env
#endif
!
!  Variables for use with Alamode
!
!  ala_mode      => string that describes the mode of the call to these routines
!  alamodedir    => path of Alamode main directory
!  ioalamode     => I/O channel for Alamode files
!  lalamode_fc2  => If true compute second order force constants
!  lalamode_fc3  => If true compute third  order force constants
!  cutoff_fc2    => cutoff for second order force constants (Ang). NB Negative means none
!  cutoff_fc3    => cutoff for second order force constants (Ang). NB Negative means none
!  disp_fc2      => displacement magnitude for second order force constants (Ang)
!  disp_fc3      => displacement magnitude for third  order force constants (Ang)
!  ioalamode     => channel number for I/O of files for Alamode
!  iodisp        => channel number for displacement file
!  ioforce       => channel number for force file
!
    character(len=maxlinelength),               save :: ala_mode
    character(len=maxlinelength),               save :: alamodedir
    character(len=maxlinelength),               save :: almexe
    character(len=maxlinelength),               save :: anphonexe
    character(len=5),                           save :: pairstring
    integer(i4),                                save :: ioalamode = 91_i4
    integer(i4),                                save :: iodisp = 92_i4
    integer(i4),                                save :: ioforce = 93_i4
    integer(i4),                                save :: iotmp = 94_i4
    integer(i4), dimension(:,:),   allocatable, save :: idisplace
    integer(i4),                                save :: nalaprocs = 0_i4
    logical,                                    save :: lalamode_fc2 = .false.
    logical,                                    save :: lalamode_fc3 = .true.
    real(dp),                                   save :: cutoff_fc2 = 8.0_dp
    real(dp),                                   save :: cutoff_fc3 = 8.0_dp
    real(dp),                                   save :: disp_fc2 = 0.005_dp
    real(dp),                                   save :: disp_fc3 = 0.01_dp
    real(dp),                                   save :: Bohr_radius = 0.52917721067  ! These conversion factors are from Alamode
    real(dp),                                   save :: Rydberg_to_eV = 13.60569253  ! These conversion factors are from Alamode

  contains

  subroutine gen_fc_alamode
!
!  Computes force constants for Alamode
!
!   2/18 Created
!   2/18 Shell relaxation added
!   3/18 Saving of coordinates added for reinitialisation
!   3/18 Testing of loptsuccess added for shells
!   1/19 maxwordlength changes added
!   9/19 Changes for modes added
!  12/20 Changes for gfortran v10
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
!  Julian Gale, CIC, Curtin University, December 2020
!
  use configurations, only : names
  use control
  use current
  use derivatives,    only : xdrv, ydrv, zdrv
  use element,        only : atsym, maxele
  use iochannels,     only : ioout
  use optimisation,   only : loptsuccess
  use parallel,       only : ioproc
#ifdef MPI
  use parallel,       only : nprocs, procid, MPI_comm_GULP
#endif
  use shells,         only : ncore, ncoptr, ncsptr, nshell
  use symmetry,       only : lsymopt
  use times,          only : talamode_fc
  implicit none
#ifdef MPI
  include 'mpif.h'
#endif
!
!  Local variables
!
  character(len=maxwordlength)                     :: alminput
  character(len=maxwordlength)                     :: almoutput
  character(len=maxwordlength)                     :: dispfile_c
  character(len=maxwordlength)                     :: dispfile_h
  character(len=maxwordlength)                     :: forcefile_c
  character(len=maxwordlength)                     :: forcefile_h
  character(len=maxlinelength)                     :: line
  character(len=maxwordlength)                     :: prefix
  character(len=maxlinelength)                     :: runalm
  integer(i4)                                      :: i
#ifdef MPI
  integer(i4)                                      :: ier
  integer(i4)                                      :: n1(1)
#endif
  integer(i4)                                      :: ii
  integer(i4)                                      :: j
  integer(i4)                                      :: nelement
  integer(i4), dimension(:),     allocatable, save :: nelementno
  integer(i4), dimension(:),     allocatable, save :: nelementrno
  integer(i4)                                      :: ni
  integer(i4)                                      :: ndisp
  integer(i4)                                      :: ndisp2
  integer(i4)                                      :: ndisp3
  integer(i4)                                      :: nr
  integer(i4)                                      :: status
  logical                                          :: lmode_all
  logical                                          :: lmode_displace
  logical                                          :: lmode_fit2
  logical                                          :: lmode_fit3
  logical                                          :: lmode_pattern
  logical                                          :: leof
  logical                                          :: lfound
  logical                                          :: lshellfail    ! Used to track whether shell optimisation failed
  logical                                          :: lshellopt     ! If true then use shell optimisation
  real(dp)                                         :: etot
  real(dp)                                         :: fconv
  real(dp)                                         :: g_cpu_time
  real(dp)                                         :: t1
  real(dp)                                         :: t2
  real(dp),    dimension(:,:),   allocatable, save :: dispxyz
  real(dp),    dimension(:,:),   allocatable, save :: savexyz
!
  t1 = g_cpu_time()
!
  lshellopt = (index(keyword,'shop').ne.0.and.nshell.gt.0)
!
!  Find what mode we are running in
!
  lmode_all = (index(ala_mode,'all').eq.1)
  lmode_displace = (index(ala_mode,'dis').eq.1)
  lmode_fit2     = (index(ala_mode,'fit_h').eq.1)
  lmode_fit3     = (index(ala_mode,'fit_c').eq.1)
  lmode_pattern  = (index(ala_mode,'pat').eq.1)
!
!  Set conversion factor for GULP derivatives to forces in Ry/au
!
  fconv = - Bohr_radius/Rydberg_to_eV
!
!  Allocate memory
!
  allocate(nelementno(ncore),stat=status)
  if (status/=0) call outofmemory('gen_fc_alamode','nelementno')
  allocate(nelementrno(maxele),stat=status)
  if (status/=0) call outofmemory('gen_fc_alamode','nelementrno')
  allocate(idisplace(4,3*ncore),stat=status)
  if (status/=0) call outofmemory('gen_fc_alamode','idisplace')
  allocate(dispxyz(3,ncore),stat=status)
  if (status/=0) call outofmemory('gen_fc_alamode','dispxyz')
  allocate(savexyz(3,numat),stat=status)
  if (status/=0) call outofmemory('gen_fc_alamode','savexyz')
!###############################################################################
!  Setup for running alamode externally if mode = all
!###############################################################################
  if (lmode_all) then
!
!  Set environment variable for Alamode path
!
    alamodedir=' '
    call getenv('ALAMODE_DIR',alamodedir)
!
!  Ensure directory ends in / if not blank
!
    call endstring(alamodedir,len(alamodedir),i)
    if (i.gt.1) then
      if (alamodedir(i-1:i-1).ne.'/') then
        alamodedir(i:i) = '/'
      endif
    endif
!
!  Create strings with exe names in
!
    almexe = trim(alamodedir) // 'alm/alm'
! 
!  Check that alm exists
!
    inquire(file=trim(almexe),exist=lfound)
!
!  If executable doesn't exist then stop
!
    if (.not.lfound) then
      call outerror('Executable alm for Alamode not found',0_i4)
      call stopnow('gen_fc_alamode')
    endif
  endif
!###############################################################################
!  General setup 
!###############################################################################
!
!  Find information for configuration needed for input
!
  nelementrno(1:maxele) = 0
  nelement = 0
  do i = 1,ncore
    ni = nat(ncoptr(i))
    if (nelementrno(ni).eq.0) then
      nelement = nelement + 1
      nelementno(nelement) = ni
      nelementrno(ni) = nelement
    endif
  enddo
!
!  Set prefix for runs based on configuration name
!
  prefix = 'gulp'
  if (names(ncf).ne.' ') then
    prefix = names(ncf)
  endif
!
!  Create name for alm input file
!
  alminput = prefix
  i = min(len_trim(prefix),maxwordlength-7_i4)
  alminput(i+1:i+7) = '_alm.in'
!
!  Create file names for displacements and forces
!
  dispfile_h = trim(prefix)//'_disp_harm.dat'
  forcefile_h = trim(prefix)//'_force_harm.dat'
  dispfile_c = trim(prefix)//'_disp_cubic.dat'
  forcefile_c = trim(prefix)//'_force_cubic.dat'
!###############################################################################
!  Output information about Alamode call
!###############################################################################
  if (ioproc) then
    write(ioout,'(/,''--------------------------------------------------------------------------------'')')
    write(ioout,'(''  Alamode Calculation : '')')
    write(ioout,'(''--------------------------------------------------------------------------------'',/)')
    write(ioout,'(''  Mode of calculation = '',a)') trim(ala_mode)
    if (lmode_all) then
      write(ioout,'(''  Alamode directory = '',a)') trim(alamodedir)
    endif
    if (lmode_all.or.lmode_displace) then
      if (lalamode_fc2) then
        write(ioout,'(/,''  Harmonic force constants: '')')
        write(ioout,'(''    Displacement = '',f8.6,'' Angstrom '')') disp_fc2
        if (cutoff_fc2.gt.0.0_dp) then
          write(ioout,'(''    Cutoff       = '',f8.6,'' Angstrom '')') cutoff_fc2
          if (cutoff_fc2.gt.0.5_dp*a.or.cutoff_fc2.gt.0.5_dp*b.or.cutoff_fc2.gt.0.5_dp*c) then
            write(ioout,'(''    WARNING: Cell is too small for this cutoff! '')') 
          endif
        else
          write(ioout,'(''    Cutoff       = None '')') 
        endif
        write(ioout,'(''    Displacements to be written to '',a)') trim(dispfile_h)
        write(ioout,'(''    Forces        to be written to '',a)') trim(forcefile_h)
      endif
      if (lalamode_fc3) then
        write(ioout,'(/,''  Cubic force constants: '')')
        write(ioout,'(''    Displacement = '',f8.6,'' Angstrom '')') disp_fc3
        if (cutoff_fc3.gt.0.0_dp) then
          write(ioout,'(''    Cutoff       = '',f8.6,'' Angstrom '')') cutoff_fc3
          if (cutoff_fc3.gt.0.5_dp*a.or.cutoff_fc3.gt.0.5_dp*b.or.cutoff_fc3.gt.0.5_dp*c) then
            write(ioout,'(''    WARNING: Cell is too small for this cutoff! '')') 
          endif
        else
          write(ioout,'(''    Cutoff       = None '')') 
        endif
        write(ioout,'(''    Displacements to be written to '',a)') trim(dispfile_c)
        write(ioout,'(''    Forces        to be written to '',a)') trim(forcefile_c)
      endif
    endif
    write(ioout,'(/)')
  endif
!
!  Save initial coordinates of atoms
!
  do i = 1,numat
    savexyz(1,i) = xclat(i)
    savexyz(2,i) = yclat(i)
    savexyz(3,i) = zclat(i)
  enddo
!###############################################################################
!  Create input required for displacement patterns
!###############################################################################
  if (lmode_all.or.lmode_pattern) then
    if (ioproc) then
!
!  Open input file for alm input
!
      open(ioalamode,file=alminput,status='unknown',form='formatted')
!
!  General block
!
      write(ioalamode,'(''&general'')')
      write(ioalamode,'(''  PREFIX = '',a)') trim(prefix)
      write(ioalamode,'(''  MODE = suggest'')')
      write(ioalamode,'(''  NAT = '',i6)') ncore
      write(ioalamode,'(''  NKD = '',i6)') nelement
      write(ioalamode,'(''  KD = '',20(a2,1x))') (atsym(nelementno(i)),i=1,nelement)
      write(ioalamode,'(''/'',/)')
!
!  Interaction block
!
      write(ioalamode,'(''&interaction'')')
      if (lalamode_fc3) then
        write(ioalamode,'(''  NORDER = 2 '')')
      else
        write(ioalamode,'(''  NORDER = 1 '')')
      endif
      write(ioalamode,'(''/'',/)')
!
!  Cell block
!
      write(ioalamode,'(''&cell'')')
      write(ioalamode,'(2x,f12.6)') 1.0_dp/Bohr_radius
      write(ioalamode,'(2x,3f12.6)') (rv(i,1),i=1,3)
      write(ioalamode,'(2x,3f12.6)') (rv(i,2),i=1,3)
      write(ioalamode,'(2x,3f12.6)') (rv(i,3),i=1,3)
      write(ioalamode,'(''/'',/)')
!
!  Cutoff block
!
      write(ioalamode,'(''&cutoff'')')
      do i = 1,nelement
        do j = 1,i
          pairstring = trim(atsym(nelementno(i)))//'-'//trim(atsym(nelementno(j)))
          if (cutoff_fc2.lt.0.0_dp.and.cutoff_fc3.lt.0.0_dp) then
            write(ioalamode,'(2x,a5,'' None None '')') pairstring
          elseif (cutoff_fc2.lt.0.0_dp) then
            write(ioalamode,'(2x,a5,'' None '',f8.4)') pairstring,cutoff_fc3/Bohr_radius
          elseif (cutoff_fc3.lt.0.0_dp) then
            write(ioalamode,'(2x,a5,1x,f8.4,'' None '')') pairstring,cutoff_fc2/Bohr_radius
          else
            write(ioalamode,'(2x,a5,1x,f8.4,1x,f8.4)') pairstring,cutoff_fc2/Bohr_radius,cutoff_fc3/Bohr_radius
          endif
        enddo
      enddo
      write(ioalamode,'(''/'',/)')
!
!  Position block
!
      write(ioalamode,'(''&position'')')
      do i = 1,ncore
        ii = ncoptr(i)
        ni = nelementrno(nat(ii))
        write(ioalamode,'(2x,i4,2x,3(f12.9,1x))') ni,xfrac(ii),yfrac(ii),zfrac(ii)
      enddo
      write(ioalamode,'(''/'',/)')
!
!  Close input file for alm input
!
      close(ioalamode)
!
      if (lmode_pattern) then
        write(ioout,'(''  Alamode input for pattern written as '',a)') trim(alminput)
      endif
    endif
  endif
!###############################################################################
!  Call alamode to set up the displacement patterns if mode = all
!###############################################################################
  if (lmode_all) then
!
!  Set up command line to run alm
!
    runalm = trim(almexe)//' '//trim(alminput)//' > alm_pattern.log'
    if (ioproc) then
      write(ioout,'(''  Computing displacement patterns..'')')
    endif
!
!  Run alm to get displacement patterns
!
    call execute_command_line(runalm,wait=.true.,exitstat=status)
    if (status/=0) then
      call outerror('attempt to run alm for displacements failed',0_i4)
      call stopnow('gen_fc_alamode')
    endif
  endif
  if (lmode_all.or.lmode_displace.or.lmode_fit2.or.lmode_fit3) then
    if (lalamode_fc2) then
!###############################################################################
!  Harmonic displacements
!###############################################################################
      lshellfail = .false.
      if (ioproc) then
!
!  Create name for alm output file with displacement patterns
!
        almoutput = trim(prefix)//'.pattern_HARMONIC'
!
!  Open file for alm displacement patterns
!
        open(ioalamode,file=almoutput,status='unknown',form='formatted')
        if (lmode_all.or.lmode_displace) then
!
!  Open output file for displacements
!
          open(iodisp,file=dispfile_h,status='unknown',form='formatted')
!
!  Open output file for forces
!
          open(ioforce,file=forcefile_h,status='unknown',form='formatted')
        endif
!
!  Read in header from file with Basis
!
        read(ioalamode,'(a)') line
!
        if (lmode_all.or.lmode_displace) then
          if (lshellopt) then
            write(ioout,'(''  Start of harmonic displacements with shell optimisation..'')')
          else
            write(ioout,'(''  Start of harmonic displacements..'')')
          endif
        endif
!******************************************************
!  Loop over displacements until the end of the file  *
!******************************************************
        leof = .false.
        ndisp2 = 0
        do while (.not.leof)
!
!  Read in displacement patterns 
!
          read(ioalamode,'(a)',end=10,err=10) line
          i = index(line,':')
          if (i.eq.0) then
            call outerror('format error in displacement file for Alamode',0_i4)
            call stopnow('gen_fc_alamode')
          endif
!
!  Increment number of displacements
!
          ndisp2 = ndisp2 + 1
!
!  Read number of atom shifts within displacement
!
          read(line(i+1:len_trim(line)),*) ni
!
!  Loop over number of displacements
!
          do i = 1,ni
            read(ioalamode,*,end=10,err=10) (idisplace(j,i),j=1,4)
          enddo
        enddo
!**************************************
!  End of harmonic displacement loop  *
!**************************************
  10    continue
!
!  Rewind file to the start
!
        rewind(ioalamode)
!
!  Read in header from file with Basis
!
        read(ioalamode,'(a)') line
      endif
#ifdef MPI
      if (nprocs.gt.1) then
!
!  Broadcast displacement patterns
!
        if (procid.eq.0) n1(1) = ndisp2
        call MPI_bcast(n1,1,MPI_integer,0,MPI_comm_GULP,ier)
        if (procid.ne.0) ndisp2 = n1(1)
      endif
#endif
      if (lmode_all.or.lmode_displace) then
!***********************************************
!  Loop over displacements and compute forces  *
!***********************************************
        do ndisp = 1,ndisp2
          if (ioproc) then
!
!  Read in displacement patterns 
!
            read(ioalamode,'(a)') line
            i = index(line,':')
!
!  Read number of atom shifts within displacement
!
            read(line(i+1:len_trim(line)),*) ni
!
!  Loop over number of displacements
!
            do i = 1,ni
              read(ioalamode,*) (idisplace(j,i),j=1,4)
            enddo
!
!  Output displacement information to indicate progress
!
            if (numat.le.99_i4) then
              write(ioout,'(''    Doing Displacement = '',i6,'' of '',i6,'' for atoms '',10i3)') ndisp,ndisp2, &
                    (idisplace(1,j),j=1,min(ni,10_i4))
            elseif (numat.le.999_i4) then
              write(ioout,'(''    Doing Displacement = '',i6,'' of '',i6,'' for atoms '',8i4)') ndisp,ndisp2, &
                    (idisplace(1,j),j=1,min(ni,8_i4))
            elseif (numat.le.9999_i4) then
              write(ioout,'(''    Doing Displacement = '',i6,'' of '',i6,'' for atoms '',6i5)') ndisp,ndisp2, &
                    (idisplace(1,j),j=1,min(ni,6_i4))
            elseif (numat.le.99999_i4) then
              write(ioout,'(''    Doing Displacement = '',i6,'' of '',i6,'' for atoms '',5i6)') ndisp,ndisp2, &
                    (idisplace(1,j),j=1,min(ni,5_i4))
            else
              write(ioout,'(''    Doing Displacement = '',i6,'' of '',i6)') ndisp,ndisp2
            endif
          endif
#ifdef MPI
          if (nprocs.gt.1) then
!
!  Broadcast displacement patterns
!
            if (procid.eq.0) n1(1) = ni
            call MPI_bcast(n1,1,MPI_integer,0,MPI_comm_GULP,ier)
            if (procid.ne.0) ni = n1(1)
            call MPI_bcast(idisplace,4*ni,MPI_integer,0,MPI_comm_GULP,ier)
          endif
#endif
!
!  Add displacements to coordinates
!
          dispxyz(1:3,1:ncore) = 0.0_dp
          do i = 1,ni
            ii = ncoptr(idisplace(1,i))
            xclat(ii) = xclat(ii) + dble(idisplace(2,i))*disp_fc2
            yclat(ii) = yclat(ii) + dble(idisplace(3,i))*disp_fc2
            zclat(ii) = zclat(ii) + dble(idisplace(4,i))*disp_fc2
!
            ii = idisplace(1,i)
            dispxyz(1,ii) = dispxyz(1,ii) + dble(idisplace(2,i))*disp_fc2/Bohr_radius
            dispxyz(2,ii) = dispxyz(2,ii) + dble(idisplace(3,i))*disp_fc2/Bohr_radius
            dispxyz(3,ii) = dispxyz(3,ii) + dble(idisplace(4,i))*disp_fc2/Bohr_radius
          enddo
          if (lshellopt) then
            do i = 1,ni
              ii = ncoptr(idisplace(1,i))
!
!  If this core has a shell then displace shell as well to avoid large separations
!
              if (ncsptr(ii).gt.0) then
                j = ncsptr(ii)
                xclat(j) = xclat(j) + dble(idisplace(2,i))*disp_fc2
                yclat(j) = yclat(j) + dble(idisplace(3,i))*disp_fc2
                zclat(j) = zclat(j) + dble(idisplace(4,i))*disp_fc2
              endif
            enddo
!
!  Relax the shells before calculating the force constants
!
            call shellopt
!
!  Check whether optimisation was successful
!
            if (.not.loptsuccess) lshellfail = .true.
          endif
!
!  Handle symmetry 
!
          if (lsymopt) then
            do i = 1,nasym
              nr = nrela2f(i)
              xalat(i) = xclat(nr)
              yalat(i) = yclat(nr)
              zalat(i) = zclat(nr)
            enddo
          else
            do i = 1,numat
              xalat(i) = xclat(i)
              yalat(i) = yclat(i)
              zalat(i) = zclat(i)
            enddo
          endif
          if (ioproc) then
!
!  Output to displacement file
!
            do i = 1,ncore
              write(iodisp,'(3f21.14)') (dispxyz(j,i),j=1,3)
            enddo
          endif
!
!  Compute forces for displacements
!
          call energy(etot,.true.,.false.)
          if (ioproc) then
!
!  Output forces to file
!
            do i = 1,ncore
              ii = ncoptr(i)
              write(ioforce,'(3e19.11)') xdrv(ii)*fconv,ydrv(ii)*fconv,zdrv(ii)*fconv
            enddo
          endif
!
!  Restore coordinates to initial values
!
          do i = 1,numat
            xclat(i) = savexyz(1,i)
            yclat(i) = savexyz(2,i)
            zclat(i) = savexyz(3,i)
          enddo
!
!  Handle symmetry 
!
          if (lsymopt) then
            do i = 1,nasym
              nr = nrela2f(i)
              xalat(i) = xclat(nr)
              yalat(i) = yclat(nr)
              zalat(i) = zclat(nr)
            enddo
          else
            do i = 1,numat
              xalat(i) = xclat(i)
              yalat(i) = yclat(i)
              zalat(i) = zclat(i)
            enddo
          endif
        enddo
!*************************************************
!  End of harmonic displacement loop for forces  *
!*************************************************
        if (ioproc) then
!
!  Close input file for alm output
!
          close(ioalamode)
!
!  Close files for coordinates and forces
!
          close(iodisp)
          close(ioforce)
          if (lshellfail) then
            write(ioout,'(''  Cubic displacements finished'')')
            write(ioout,'(''  WARNING: Not all shell optimisations were successful'',/)')
          else
            write(ioout,'(''  Cubic displacements finished'',/)')
          endif
        endif
      endif
      if (ioproc) then
!
!  Close input file for alm output
!
        close(ioalamode)
      endif
    endif
    if (lalamode_fc3) then
!###############################################################################
!  Anharmonic displacements
!###############################################################################
      lshellfail = .false.
      if (ioproc) then
!
!  Create name for alm output file with displacement patterns
!
        almoutput = trim(prefix)//'.pattern_ANHARM3'
!
!  Open file for alm displacement patterns
!
        open(ioalamode,file=almoutput,status='unknown',form='formatted')
        if (lmode_all.or.lmode_displace) then
!
!  Open output file for displacements
!
          open(iodisp,file=dispfile_c,status='unknown',form='formatted')
!
!  Open output file for forces
!
          open(ioforce,file=forcefile_c,status='unknown',form='formatted')
        endif
!
!  Read in header from file with Basis
!
        read(ioalamode,'(a)') line
!
        if (lmode_all.or.lmode_displace) then
          if (lshellopt) then
            write(ioout,'(''  Start of cubic displacements with shell optimisation..'')')
          else
            write(ioout,'(''  Start of cubic displacements..'')')
          endif
        endif
!******************************************************
!  Loop over displacements until the end of the file  *
!******************************************************
        leof = .false.
        ndisp3 = 0
        do while (.not.leof)
!
!  Read in displacement patterns 
!
          read(ioalamode,'(a)',end=20,err=20) line
          i = index(line,':')
          if (i.eq.0) then
            call outerror('format error in displacement file for Alamode',0_i4)
            call stopnow('gen_fc_alamode')
          endif
!
!  Increment number of displacements
!
          ndisp3 = ndisp3 + 1 
!
!  Read number of atom shifts within displacement
!
          read(line(i+1:len_trim(line)),*) ni
!
!  Loop over number of displacements
!
          do i = 1,ni
            read(ioalamode,*,end=20,err=20) (idisplace(j,i),j=1,4)
          enddo
        enddo
!****************************************
!  End of anharmonic displacement loop  *
!****************************************
  20    continue
!
!  Rewind file
!
        rewind(ioalamode)
!
!  Read in header from file with Basis
!
        read(ioalamode,'(a)') line
      endif
#ifdef MPI
      if (nprocs.gt.1) then
!
!  Broadcast displacement patterns
!
        if (procid.eq.0) n1(1) = ndisp3
        call MPI_bcast(n1,1,MPI_integer,0,MPI_comm_GULP,ier)
        if (procid.ne.0) ndisp3 = n1(1)
      endif
#endif
!***********************************************
!  Loop over displacements and compute forces  *
!***********************************************
      if (lmode_all.or.lmode_displace) then
        do ndisp = 1,ndisp3
          if (ioproc) then
!
!  Read in displacement patterns 
!
            read(ioalamode,'(a)') line
            i = index(line,':')
!
!  Read number of atom shifts within displacement
!
            read(line(i+1:len_trim(line)),*) ni
!
!  Loop over number of displacements
!
            do i = 1,ni
              read(ioalamode,*) (idisplace(j,i),j=1,4)
            enddo
!
!  Output displacement information to indicate progress
!
            if (numat.le.99_i4) then
              write(ioout,'(''    Doing Displacement = '',i6,'' of '',i6,'' for atoms '',10i3)') ndisp,ndisp3, &
                    (idisplace(1,j),j=1,min(ni,10_i4))
            elseif (numat.le.999_i4) then
              write(ioout,'(''    Doing Displacement = '',i6,'' of '',i6,'' for atoms '',8i4)') ndisp,ndisp3, &
                    (idisplace(1,j),j=1,min(ni,8_i4))
            elseif (numat.le.9999_i4) then
              write(ioout,'(''    Doing Displacement = '',i6,'' of '',i6,'' for atoms '',6i5)') ndisp,ndisp3, &
                    (idisplace(1,j),j=1,min(ni,6_i4))
            elseif (numat.le.99999_i4) then
              write(ioout,'(''    Doing Displacement = '',i6,'' of '',i6,'' for atoms '',5i6)') ndisp,ndisp3, &
                    (idisplace(1,j),j=1,min(ni,5_i4))
            else
              write(ioout,'(''    Doing Displacement = '',i6,'' of '',i6)') ndisp,ndisp3
            endif
          endif
#ifdef MPI
          if (nprocs.gt.1) then
!
!  Broadcast displacement patterns
!
            if (procid.eq.0) n1(1) = ni
            call MPI_bcast(n1,1,MPI_integer,0,MPI_comm_GULP,ier)
            if (procid.ne.0) ni = n1(1)
            call MPI_bcast(idisplace,4*ni,MPI_integer,0,MPI_comm_GULP,ier)
          endif
#endif
!
!  Add displacements to coordinates
!
          dispxyz(1:3,1:ncore) = 0.0_dp
          do i = 1,ni
            ii = ncoptr(idisplace(1,i))
            xclat(ii) = xclat(ii) + dble(idisplace(2,i))*disp_fc3
            yclat(ii) = yclat(ii) + dble(idisplace(3,i))*disp_fc3
            zclat(ii) = zclat(ii) + dble(idisplace(4,i))*disp_fc3
!
            ii = idisplace(1,i)
            dispxyz(1,ii) = dispxyz(1,ii) + dble(idisplace(2,i))*disp_fc3/Bohr_radius
            dispxyz(2,ii) = dispxyz(2,ii) + dble(idisplace(3,i))*disp_fc3/Bohr_radius
            dispxyz(3,ii) = dispxyz(3,ii) + dble(idisplace(4,i))*disp_fc3/Bohr_radius
          enddo
          if (lshellopt) then
            do i = 1,ni
              ii = ncoptr(idisplace(1,i))
!
!  If this core has a shell then displace shell as well to avoid large separations
!
              if (ncsptr(ii).gt.0) then
                j = ncsptr(ii)
                xclat(j) = xclat(j) + dble(idisplace(2,i))*disp_fc2
                yclat(j) = yclat(j) + dble(idisplace(3,i))*disp_fc2
                zclat(j) = zclat(j) + dble(idisplace(4,i))*disp_fc2
              endif
            enddo
!
!  Relax the shells before calculating the force constants
!
            call shellopt
!
!  Check whether optimisation was successful
!
            if (.not.loptsuccess) lshellfail = .true.
          endif
!
!  Handle symmetry 
!
          if (lsymopt) then
            do i = 1,nasym
              nr = nrela2f(i)
              xalat(i) = xclat(nr)
              yalat(i) = yclat(nr)
              zalat(i) = zclat(nr)
            enddo
          else
            do i = 1,numat
              xalat(i) = xclat(i)
              yalat(i) = yclat(i)
              zalat(i) = zclat(i)
            enddo
          endif
          if (ioproc) then
!
!  Output to displacement file
!
            do i = 1,ncore
              write(iodisp,'(3f21.14)') (dispxyz(j,i),j=1,3)
            enddo
          endif
!
!  Compute forces for displacements
!
          call energy(etot,.true.,.false.)
          if (ioproc) then
!
!  Output forces to file
!
            do i = 1,ncore
              ii = ncoptr(i)
              write(ioforce,'(3e19.11)') xdrv(ii)*fconv,ydrv(ii)*fconv,zdrv(ii)*fconv
            enddo
          endif
!
!  Restore coordinates to initial values
!
          do i = 1,numat
            xclat(i) = savexyz(1,i) 
            yclat(i) = savexyz(2,i) 
            zclat(i) = savexyz(3,i) 
          enddo
!
!  Handle symmetry 
!
          if (lsymopt) then
            do i = 1,nasym
              nr = nrela2f(i)
              xalat(i) = xclat(nr)
              yalat(i) = yclat(nr)
              zalat(i) = zclat(nr)
            enddo
          else
            do i = 1,numat
              xalat(i) = xclat(i)
              yalat(i) = yclat(i)
              zalat(i) = zclat(i)
            enddo
          endif
        enddo
      endif
!***************************************************
!  End of anharmonic displacement loop for forces  *
!***************************************************
      if (ioproc) then
!
!  Close input file for alm output
!
        close(ioalamode)
        if (lmode_all.or.lmode_displace) then
!
!  Close files for coordinates and forces
!
          close(iodisp)
          close(ioforce)
          if (lshellfail) then
            write(ioout,'(''  Cubic displacements finished'')')
            write(ioout,'(''  WARNING: Not all shell optimisations were successful'',/)')
          else
            write(ioout,'(''  Cubic displacements finished'',/)')
          endif
        endif
      endif
    endif
    if (lmode_all.or.lmode_displace) then
!
!  If this is a shell optimisation then we need to call shellopt reset the fractional coordinates
!
      if (lshellopt) then
        call shellopt
      endif
    endif
  endif
!
  if (lmode_all.or.lmode_fit2) then
    if (lalamode_fc2) then
!###############################################################################
!  Create input required for fitting harmonic force constants
!###############################################################################
      if (ioproc) then
!
!  Open input file for alm input
!
        open(ioalamode,file=alminput,status='unknown',form='formatted')
!
!  General block
!
        write(ioalamode,'(''&general'')')
        write(ioalamode,'(''  PREFIX = '',a)') trim(prefix)//'_harm'
        write(ioalamode,'(''  MODE = fitting'')')
        write(ioalamode,'(''  NAT = '',i6)') ncore
        write(ioalamode,'(''  NKD = '',i6)') nelement
        write(ioalamode,'(''  KD = '',20(a2,1x))') (atsym(nelementno(i)),i=1,nelement)
        write(ioalamode,'(''/'',/)')
!
!  Fitting block
!
        write(ioalamode,'(''&fitting'')')
        write(ioalamode,'(''  NDATA = '',i6)') ndisp2
        write(ioalamode,'(''  DFILE = '',a)') dispfile_h
        write(ioalamode,'(''  FFILE = '',a)') forcefile_h
        write(ioalamode,'(''/'',/)')
!
!  Interaction block
!
        write(ioalamode,'(''&interaction'')')
        write(ioalamode,'(''  NORDER = 1 '')')
        write(ioalamode,'(''/'',/)')
!
!  Cell block
!
        write(ioalamode,'(''&cell'')')
        write(ioalamode,'(2x,f12.6)') 1.0_dp/Bohr_radius
        write(ioalamode,'(2x,3f12.6)') (rv(i,1),i=1,3)
        write(ioalamode,'(2x,3f12.6)') (rv(i,2),i=1,3)
        write(ioalamode,'(2x,3f12.6)') (rv(i,3),i=1,3)
        write(ioalamode,'(''/'',/)')
!
!  Cutoff block
!
        write(ioalamode,'(''&cutoff'')')
        do i = 1,nelement
          do j = 1,i
            pairstring = trim(atsym(nelementno(i)))//'-'//trim(atsym(nelementno(j)))
            if (cutoff_fc2.lt.0.0_dp.and.cutoff_fc3.lt.0.0_dp) then
              write(ioalamode,'(2x,a5,'' None None '')') pairstring
            elseif (cutoff_fc2.lt.0.0_dp) then
              write(ioalamode,'(2x,a5,'' None '',f8.4)') pairstring,cutoff_fc3/Bohr_radius
            elseif (cutoff_fc3.lt.0.0_dp) then
              write(ioalamode,'(2x,a5,1x,f8.4,'' None '')') pairstring,cutoff_fc2/Bohr_radius
            else
              write(ioalamode,'(2x,a5,1x,f8.4,1x,f8.4)') pairstring,cutoff_fc2/Bohr_radius,cutoff_fc3/Bohr_radius
            endif
          enddo
        enddo
        write(ioalamode,'(''/'',/)')
!
!  Position block
!
        write(ioalamode,'(''&position'')')
        do i = 1,ncore
          ii = ncoptr(i)
          ni = nelementrno(nat(ii))
          write(ioalamode,'(2x,i4,2x,3(f12.9,1x))') ni,xfrac(ii),yfrac(ii),zfrac(ii)
        enddo
        write(ioalamode,'(''/'',/)')
!
!  Close input file for alm input
!
        close(ioalamode)
        if (lmode_all) then
!
!  Set up command line to run alm
!
          runalm = trim(almexe)//' '//trim(alminput)//' >> alm_fit.log'
!
!  Run alm to get harmonic force constants
!
          call execute_command_line(runalm,wait=.true.,exitstat=status)
          if (status/=0) then
            call outerror('attempt to run alm for force constants failed',0_i4)
            call stopnow('gen_fc_alamode')
          endif
          write(ioout,'(''  Harmonic force constants fitted'')')
        else
          write(ioout,'(''  Alamode input for harmonic fit written as '',a)') trim(alminput)
        endif
      endif
    endif
  endif
  if (lmode_all.or.lmode_fit3) then
    if (lalamode_fc3) then
!###############################################################################
!  Create input required for fitting cubic force constants
!###############################################################################
      if (ioproc) then
!
!  Open input file for alm input
!
        open(ioalamode,file=alminput,status='unknown',form='formatted')
!
!  General block
!
        write(ioalamode,'(''&general'')')
        write(ioalamode,'(''  PREFIX = '',a)') trim(prefix)//'_cubic'
        write(ioalamode,'(''  MODE = fitting'')')
        write(ioalamode,'(''  NAT = '',i6)') ncore
        write(ioalamode,'(''  NKD = '',i6)') nelement
        write(ioalamode,'(''  KD = '',20(a2,1x))') (atsym(nelementno(i)),i=1,nelement)
        write(ioalamode,'(''/'',/)')
!
!  Fitting block
!
        write(ioalamode,'(''&fitting'')')
        write(ioalamode,'(''  NDATA = '',i6)') ndisp3
        write(ioalamode,'(''  DFILE = '',a)') dispfile_c
        write(ioalamode,'(''  FFILE = '',a)') forcefile_c
        write(ioalamode,'(''/'',/)')
!
!  Interaction block
!
        write(ioalamode,'(''&interaction'')')
        write(ioalamode,'(''  NORDER = 2 '')')
        write(ioalamode,'(''/'',/)')
!
!  Cell block
!
        write(ioalamode,'(''&cell'')')
        write(ioalamode,'(2x,f12.6)') 1.0_dp/Bohr_radius
        write(ioalamode,'(2x,3f12.6)') (rv(i,1),i=1,3)
        write(ioalamode,'(2x,3f12.6)') (rv(i,2),i=1,3)
        write(ioalamode,'(2x,3f12.6)') (rv(i,3),i=1,3)
        write(ioalamode,'(''/'',/)')
!
!  Cutoff block
!
        write(ioalamode,'(''&cutoff'')')
        do i = 1,nelement
          do j = 1,i
            pairstring = trim(atsym(nelementno(i)))//'-'//trim(atsym(nelementno(j)))
            if (cutoff_fc2.lt.0.0_dp.and.cutoff_fc3.lt.0.0_dp) then
              write(ioalamode,'(2x,a5,'' None None '')') pairstring
            elseif (cutoff_fc2.lt.0.0_dp) then
              write(ioalamode,'(2x,a5,'' None '',f8.4)') pairstring,cutoff_fc3/Bohr_radius
            elseif (cutoff_fc3.lt.0.0_dp) then
              write(ioalamode,'(2x,a5,1x,f8.4,'' None '')') pairstring,cutoff_fc2/Bohr_radius
            else
              write(ioalamode,'(2x,a5,1x,f8.4,1x,f8.4)') pairstring,cutoff_fc2/Bohr_radius,cutoff_fc3/Bohr_radius
            endif
          enddo
        enddo
        write(ioalamode,'(''/'',/)')
!
!  Position block
!
        write(ioalamode,'(''&position'')')
        do i = 1,ncore
          ii = ncoptr(i)
          ni = nelementrno(nat(ii))
          write(ioalamode,'(2x,i4,2x,3(f12.9,1x))') ni,xfrac(ii),yfrac(ii),zfrac(ii)
        enddo
        write(ioalamode,'(''/'',/)')
!
!  Close input file for alm input
!
        close(ioalamode)
!
        if (lmode_all) then
!
!  Set up command line to run alm
!
          runalm = trim(almexe)//' '//trim(alminput)//' >> alm_fit.log'
!
!  Run alm to get cubic force constants
!
          call execute_command_line(runalm,wait=.true.,exitstat=status)
          if (status/=0) then
            call outerror('attempt to run alm for force constants failed',0_i4)
            call stopnow('gen_fc_alamode')
          endif
          write(ioout,'(''  Cubic force constants fitted'')')
        else
          write(ioout,'(''  Alamode input for cubic fit written as '',a)') trim(alminput)
        endif
      endif
    endif
  endif
!
!  Deallocate memory
!
  deallocate(dispxyz,stat=status)
  if (status/=0) call deallocate_error('gen_fc_alamode','dispxyz')
  deallocate(idisplace,stat=status)
  if (status/=0) call deallocate_error('gen_fc_alamode','idisplace')
  deallocate(nelementrno,stat=status)
  if (status/=0) call deallocate_error('gen_fc_alamode','nelementrno')
  deallocate(nelementno,stat=status)
  if (status/=0) call deallocate_error('gen_fc_alamode','nelementno')
!
  t2 = g_cpu_time()
  talamode_fc = talamode_fc + t2 - t1
!
  return
  end subroutine gen_fc_alamode

  subroutine anharmonic_alamode
!
!  Computes anharmonic properties from Alamode:
!  - Thermal conductivity
!  - Phonon lifetimes
!
!   2/18 Created
!   8/18 Additional options passed through and defaults set
!   1/19 maxwordlength changes added
!   9/19 Use of intochar replace by string write
!   9/19 Changes for mode added
!  11/19 Check for version number added to handle CLASSICAL option
!   7/20 Version number now obtained alm_fit.log as name has changed
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
!  Julian Gale, CIC, Curtin University, July 2020
!
  use configurations, only : names, rvpcfg, ncorepcfg
  use configurations, only : ntempstp, tempstp, tempcfg, ntempramp
  use control
  use current
  use element,        only : atsym, atmass, maxele
  use iochannels,     only : ioout
  use ksample,        only : nksala
  use parallel,       only : ioproc, nprocs
  use shells,         only : ncore, ncoptr
  use times,          only : talamode_tc
  implicit none
#ifdef MPI
  include 'mpif.h'
#endif
!
!  Local variables
!
  character(len=maxwordlength)                     :: anphoninput
  character(len=maxwordlength)                     :: anphonoutput
  character(len=maxwordlength)                     :: fcsfile
  character(len=maxwordlength)                     :: line
  character(len=maxwordlength)                     :: prefix
  character(len=5)                                 :: procs
  character(len=maxlinelength)                     :: runanphon
  integer(i4)                                      :: i
  integer(i4)                                      :: ind1
  integer(i4)                                      :: ind2
  integer(i4)                                      :: j
  integer(i4)                                      :: nelement
  integer(i4), dimension(:),     allocatable, save :: nelementno
  integer(i4), dimension(:),     allocatable, save :: nelementrno
  integer(i4)                                      :: ni
  integer(i4)                                      :: nprim
  integer(i4)                                      :: status
  logical                                          :: lerror
  logical                                          :: lfound
  logical                                          :: lmode_all
  logical                                          :: lmode_bte
  logical                                          :: lversion0
  real(dp)                                         :: g_cpu_time
  real(dp)                                         :: tmax
  real(dp)                                         :: tmin
  real(dp)                                         :: t1
  real(dp)                                         :: t2
!
  t1 = g_cpu_time()
!
!  Set mode flags
!
  lmode_all = (index(ala_mode,'all').eq.1)
  lmode_bte = (index(ala_mode,'bte').eq.1)
!
  if (lmode_all.or.lmode_bte) then
!
    if (ioproc) then
!
!  Allocate memory
!
      allocate(nelementno(ncore),stat=status)
      if (status/=0) call outofmemory('anharmonic_alamode','nelementno')
      allocate(nelementrno(maxele),stat=status)
      if (status/=0) call outofmemory('anharmonic_alamode','nelementrno')
!
      if (lmode_all) then
!
!  Set environment variable for Alamode path
!
        alamodedir=' '
        call getenv('ALAMODE_DIR',alamodedir)
!
!  Ensure directory ends in / if not blank
!
        call endstring(alamodedir,len(alamodedir),i)
        if (i.gt.1) then
          if (alamodedir(i-1:i-1).ne.'/') then
            alamodedir(i:i) = '/'
          endif
        endif
!
!  Create strings with exe names in
!
        anphonexe = trim(alamodedir) // 'anphon/anphon'
! 
!  Check that anphon exists
!
        inquire(file=trim(anphonexe),exist=lfound)
!
!  If executable doesn't exist then stop
!
        if (.not.lfound) then
          call outerror('Executable anphon for Alamode not found',0_i4)
          call stopnow('anharmonic_alamode')
        endif
      endif
!
!  Find information for configuration needed for input
!
      nelementrno(1:maxele) = 0
      nelement = 0
      do i = 1,ncore
        ni = nat(ncoptr(i))
        if (nelementrno(ni).eq.0) then
          nelement = nelement + 1
          nelementno(nelement) = ni
          nelementrno(ni) = nelement
        endif
      enddo
!
!  Set prefix for runs based on configuration name
!
      prefix = 'gulp'
      if (names(ncf).ne.' ') then
        prefix = names(ncf)
      endif
!
!  Create name for anphon input file
!
      anphoninput = prefix
      i = min(len_trim(prefix),maxwordlength-10_i4)
      anphoninput(i+1:i+10) = '_anphon.in'
!
!  Create name for anphon output file with thermal conductivity
!
      anphonoutput = trim(prefix)//'_rta.kl'
!
!  Create file names for force constants
!
      fcsfile = trim(prefix)//'_cubic.xml'
!
!  Find out which version of Alamode we are using
!
      runanphon = 'grep Ver alm_fit.log > almgrep.tmp'
!
!  Perform command to get version number
!
      call execute_command_line(runanphon,wait=.true.,exitstat=status)
      if (status/=0) then
        call outerror('failed to find Alamode version number',0_i4)
        call stopnow('anharmonic_alamode')
      endif
      lerror = .true.
      lversion0 = .false.
      open(iotmp,file='almgrep.tmp',status='old',form='formatted',err=5)
      read(iotmp,'(a)',err=15,end=15) line
      lversion0 = (index(line,'0.9').ne.0)
  15  continue
      close(iotmp,status='delete')
   5  continue
!
!  Find out how many atoms there should be in the primitive cell according to alm
!
      runanphon = 'grep Primitive alm_fit.log | tail -n 1 > almgrep.tmp'
!
!  Run anphon to get thermal conductivity
!
      call execute_command_line(runanphon,wait=.true.,exitstat=status)
      if (status/=0) then
        call outerror('failed to find number of atoms in the primitive cell',0_i4)
        call stopnow('anharmonic_alamode')
      endif
      lerror = .true.
      nprim = 0
      open(iotmp,file='almgrep.tmp',status='old',form='formatted',err=10)
      read(iotmp,'(a)',err=20,end=20) line
      ind1 = index(line,'Primitive cell contains') + 23_i4
      ind2 = index(line,'atoms')
      if (ind2.gt.ind1) then
        read(line(ind1:ind2),*,err=20) nprim
      endif
      lerror = .false.
  20  continue
      close(iotmp,status='delete')
  10  continue
!###############################################################################
!  Output information about Alamode anharmonic property files
!###############################################################################
      if (lalamode_fc3) then
        write(ioout,'(/,''  Anharmonic phonon calculation: '')')
        write(ioout,'(''    Primitive cell vectors (Angstrom): '')')
        write(ioout,'(''      '',3(f10.5,1x))') (rvpcfg(j,1,ncf),j=1,3)
        write(ioout,'(''      '',3(f10.5,1x))') (rvpcfg(j,2,ncf),j=1,3)
        write(ioout,'(''      '',3(f10.5,1x))') (rvpcfg(j,3,ncf),j=1,3)
        if (lerror) then
          write(ioout,'(''    WARNING: Unable to check number of atoms in the primitive cell! '')')
        else
          if (nprim.eq.ncorepcfg(ncf)) then
            write(ioout,'(''    Number of atoms in the primitive cell = '',i6)') nprim
          else
            write(ioout,'(''    WARNING: Number of primitive atoms does not match Alamode: '',i6, &
                        & '' vs '',i6)') nprim,ncorepcfg(ncf)
          endif
        endif
        write(ioout,'(''    Shrinking factors = '',3(i4,1x))') nksala(1,ncf),nksala(2,ncf),nksala(3,ncf)
        write(ioout,'(''    Force constants to be read from       '',a)') trim(fcsfile)
        write(ioout,'(''    Thermal conductivity to be written to '',a)') trim(anphonoutput)
      endif
      write(ioout,'(/)')
!###############################################################################
!  Create input required for thermal conductivity
!###############################################################################
!
!  Open input file for anphon input
!
      open(ioalamode,file=anphoninput,status='unknown',form='formatted')
!
!  General block
!
      write(ioalamode,'(''&general'')')
      write(ioalamode,'(''  PREFIX = '',a)') trim(prefix)//'_rta'
      write(ioalamode,'(''  MODE = RTA'')')
      write(ioalamode,'(''  NKD = '',i6)') nelement
      write(ioalamode,'(''  KD = '',10(a2,1x))') (atsym(nelementno(i)),i=1,nelement)
      write(ioalamode,'(''  MASS = '',10(f8.4,1x))') (atmass(nelementno(i)),i=1,nelement)
      if (ntempramp(ncf).gt.0) then
        tmin = tempcfg(ncf)
        tmax = tempcfg(ncf) + dble(ntempstp(1,ncf))*tempstp(1,ncf)
        if (tmax.gt.tmin) then
          write(ioalamode,'(''  TMIN = '',f12.4)') tmin
          write(ioalamode,'(''  TMAX = '',f12.4)') tmax
          write(ioalamode,'(''  DT =   '',f12.4)') tempstp(1,ncf)
        else
          write(ioalamode,'(''  TMIN = '',f12.4)') tmax
          write(ioalamode,'(''  TMAX = '',f12.4)') tmin
          write(ioalamode,'(''  DT =   '',f12.4)') tempstp(1,ncf)
        endif
      else
        tmin = tempcfg(ncf)
        write(ioalamode,'(''  TMIN = '',f12.4)') tmin
        write(ioalamode,'(''  TMAX = '',f12.4)') tmin
        write(ioalamode,'(''  DT =   '',f12.4)') 10.0_dp
      endif
      if (.not.lversion0) then
        write(ioalamode,'(''  CLASSICAL = 0 '')')
      endif
      write(ioalamode,'(''  FCSXML = '',a)') fcsfile
      write(ioalamode,'(''  RESTART = 0 '')')
      write(ioalamode,'(''/'',/)')
!
!  Cell block
!
      write(ioalamode,'(''&cell'')')
      write(ioalamode,'(2x,f12.6)') 1.0_dp/Bohr_radius
      write(ioalamode,'(2x,3f12.6)') (rvpcfg(i,1,ncf),i=1,3)
      write(ioalamode,'(2x,3f12.6)') (rvpcfg(i,2,ncf),i=1,3)
      write(ioalamode,'(2x,3f12.6)') (rvpcfg(i,3,ncf),i=1,3)
      write(ioalamode,'(''/'',/)')
!
!  K points
!
      write(ioalamode,'(''&kpoint'')')
      write(ioalamode,'(2x,''2'')')  ! Mode indicator - uniform k point grid
      if (nksala(1,ncf)*nksala(2,ncf)*nksala(3,ncf).eq.0) then
        write(ioalamode,'(2x,3(i4,1x))') 1_i4,1_i4,1_i4
      else
        write(ioalamode,'(2x,3(i4,1x))') nksala(1,ncf),nksala(2,ncf),nksala(3,ncf)
      endif
      write(ioalamode,'(''/'',/)')
!
!  Close input file for anphon input
!
      close(ioalamode)
!
      if (lmode_all) then
!
!  Set up command line to run anphon
!
        if (nalaprocs.ne.0) then
          if (nalaprocs.gt.1) then
            write(procs,'(i5)') nalaprocs
            runanphon = 'mpirun -np '//trim(adjustl(procs))//' '//trim(anphonexe)//' '//trim(anphoninput)//' >> anphon.log'
          else
            runanphon = trim(anphonexe)//' '//trim(anphoninput)//' >> anphon.log'
          endif
        else
          if (nprocs.gt.1) then
            write(procs,'(i5)') nalaprocs
            runanphon = 'mpirun -np '//trim(adjustl(procs))//' '//trim(anphonexe)//' '//trim(anphoninput)//' >> anphon.log'
          else
            runanphon = trim(anphonexe)//' '//trim(anphoninput)//' >> anphon.log'
          endif
        endif
!
!  Run anphon to get thermal conductivity
!
        call execute_command_line(runanphon,wait=.true.,exitstat=status)
        if (status/=0) then
          call outerror('attempt to run anphon for thermal conductivity failed',0_i4)
          call stopnow('anharmonic_alamode')
        endif
      else
        write(ioout,'(''  Alamode input for BTE written as '',a)') trim(anphoninput)
      endif
!
!  Deallocate memory
!
      deallocate(nelementrno,stat=status)
      if (status/=0) call deallocate_error('anharmonic_alamode','nelementrno')
      deallocate(nelementno,stat=status)
      if (status/=0) call deallocate_error('anharmonic_alamode','nelementno')
    endif
#ifdef MPI
    if (nprocs.gt.1) then
!
!  Get all processors to wait here for completion of thermal conductivity
!
      call mpbarrier
    endif
#endif
!###############################################################################
!  Post-processing
!###############################################################################
  endif
!
  t2 = g_cpu_time()
  talamode_tc = talamode_tc + t2 - t1
!
  return
  end subroutine anharmonic_alamode

  end module m_alamode
