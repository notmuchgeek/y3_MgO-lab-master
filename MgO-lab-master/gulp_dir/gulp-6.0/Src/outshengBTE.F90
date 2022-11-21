  subroutine outshengBTEcontrol(iout,lusefc)
!
!  Outputs file CONTROL for ShengBTE.
!
!   3/14 Created
!   4/14 Corrections made based on debug testing
!   2/15 Modified to allow for supercell force constants
!   3/15 Supercell corrected to 2*nd2cell+1
!   4/15 Modified to use ghost supercell information
!   4/15 Corrections made to output format
!   4/15 Minimum temperature set to 1 K
!   7/15 Corrected so that shells are excluded from output
!  12/16 Modified so that only ioproc executes this routine
!   2/18 nxks, nyks, nzks converted to a single array
!   1/19 Changes for maxlinelength made
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!
!  Assumptions:
!
!  - cores are sorted to come before shells
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
!  Julian Gale, CIC, Curtin University, September 2019
!
  use configurations, only : nsuperghost
  use control
  use current 
  use derivatives,    only : nd2cell
  use element,        only : maxele
  use gulp_lengths
  use ksample
  use m_pdfneutron,   only : lshrinkset
  use parallel,       only : ioproc
  use properties,     only : dicons
  use shells,         only : ncore
  use species,        only : nspec, natspec
  implicit none
!
!  Passed variables
!
  integer(i4),                         intent(in)   :: iout             ! Output channel for file
  logical,                             intent(in)   :: lusefc           ! If true then set scell to be fc_supercell
!
!  Local variables
!
  character(len=5)                                  :: lab
  integer(i4)                                       :: i
  integer(i4)                                       :: ii
  integer(i4)                                       :: ig
  integer(i4)                                       :: j
  integer(i4)                                       :: nghostatom
  integer(i4)                                       :: nghostcell
  integer(i4)                                       :: nghostcore
  integer(i4)                                       :: nghostspec
  integer(i4),      dimension(:), allocatable, save :: nghostspecptr
  integer(i4)                                       :: status
  real(dp)                                          :: rvtmp(3,3)
  real(dp)                                          :: xf
  real(dp)                                          :: yf
  real(dp)                                          :: zf
!
!  Only ioproc needs to execute this routine
!
  if (.not.ioproc) return
!
!  Allocate memory
!
  allocate(nghostspecptr(nspec),stat=status)
  if (status/=0) call outofmemory('outshengBTEcontrol','nghostspecptr')
!
!  Setup up pointer to ghost species
!
  nghostspec = 0
  do i = 1,nspec
    if (natspec(i).le.maxele) then
      nghostspec = nghostspec + 1
      nghostspecptr(nghostspec) = i
    endif
  enddo
!-------------------------------------------------------------------------------
!  Setup tasks for ghost cells
!-------------------------------------------------------------------------------
!
!  Find number of ghost cells
!
  nghostcell = nsuperghost(1,ncf)*nsuperghost(2,ncf)*nsuperghost(3,ncf)
!
!  Check the number of atoms is divisable by the number of ghost cells
!
  if (mod(numat,nghostcell).ne.0) then
    call outerror('number of atoms is not compatible with ghost_supercell',0_i4)
    call stopnow('outshengBTEcontrol')
  endif
!
  nghostatom = numat/nghostcell
  nghostcore = ncore/nghostcell
!
!  Generate original cell prior to supercell
!
  rvtmp(1:3,1) = rv(1:3,1)/dble(nsuperghost(1,ncf))
  rvtmp(1:3,2) = rv(1:3,2)/dble(nsuperghost(2,ncf))
  rvtmp(1:3,3) = rv(1:3,3)/dble(nsuperghost(3,ncf))
!-------------------------------------------------------------------------------
!  Open file
!-------------------------------------------------------------------------------
  open(iout,file='CONTROL',status='unknown',form='formatted')
!-------------------------------------------------------------------------------
!  Allocations block
!-------------------------------------------------------------------------------
  write(iout,'(''&allocations'')')
  write(iout,'(8x,''nelements='',i6)',advance='no') nghostspec
  write(iout,'('','')')
  write(iout,'(8x,''natoms='',i6)',advance='no') nghostcore
  write(iout,'('','')')
  if (lshrinkset(ncf)) then
    write(iout,'(8x,''ngrid(:)='',3i4)') nks(1,ncf),nks(2,ncf),nks(3,ncf)
  else
    write(iout,'(8x,''ngrid(:)=1 1 1'')') 
  endif
  write(iout,'(''&end'')')
!-------------------------------------------------------------------------------
!  Crystal block
!-------------------------------------------------------------------------------
  write(iout,'(''&crystal'')')
!
!  Output lattice factor as conversion from Angstroms to nm
!
  write(iout,'(8x,''lfactor=0.1'')',advance='no') 
  write(iout,'('','')')
!
!  Lattice vectors
!
  write(iout,'(8x,''lattvec(:,1)='',3(f12.6,1x))',advance='no') (rvtmp(j,1),j=1,3)
  write(iout,'('','')')
  write(iout,'(8x,''lattvec(:,2)='',3(f12.6,1x))',advance='no') (rvtmp(j,2),j=1,3)
  write(iout,'('','')')
  write(iout,'(8x,''lattvec(:,3)='',3(f12.6,1x))',advance='no') (rvtmp(j,3),j=1,3)
  write(iout,'('','')')
!
!  Elements
!
  write(iout,'(8x,''elements='')',advance='no') 
  do ii = 1,nghostspec
    i = nghostspecptr(ii)
    write(iout,'(''"'')',advance='no') 
    call label(natspec(i),0_i4,lab)
    write(iout,'(a)',advance='no') trim(lab)
    if (ii.eq.nghostspec) then
      write(iout,'(''"'')') 
    else
      write(iout,'(''" '')',advance='no') 
    endif
  enddo
  write(iout,'(8x,''types='')',advance='no')
  do ig = 1,nghostcore
    i = (ig - 1)*nghostcell + 1
    write(iout,'(i3,1x)',advance='no') nspecptr(nrelf2a(i))
  enddo
  write(iout,'('','')')
!
!  Atomic fractional coordinates
!
  do ig = 1,nghostcore
    i = (ig - 1)*nghostcell + 1
    xf = xfrac(i)*dble(nsuperghost(1,ncf))
    yf = yfrac(i)*dble(nsuperghost(2,ncf))
    zf = zfrac(i)*dble(nsuperghost(3,ncf))
    write(iout,'(8x,''positions(:,'',i3,'')='',3(f12.6,1x))',advance='no') ig,xf,yf,zf
    write(iout,'('','')')
  enddo
!
!  Static dielectric constant tensor
!
  write(iout,'(8x,''epsilon(:,1)='',3(f12.6,1x))',advance='no') (dicons(j,1),j=1,3)
  write(iout,'('','')')
  write(iout,'(8x,''epsilon(:,2)='',3(f12.6,1x))',advance='no') (dicons(j,2),j=1,3)
  write(iout,'('','')')
  write(iout,'(8x,''epsilon(:,3)='',3(f12.6,1x))',advance='no') (dicons(j,3),j=1,3)
  write(iout,'('','')')
!
!  Born effective charges
!
  do ig = 1,nghostcore
    i = (ig - 1)*nghostcell + 1
    write(iout,'(8x,''born(:,1,'',i3,'')='',3(f12.6,1x))',advance='no') ig,(bornq(j,1,i),j=1,3)
    write(iout,'('','')')
    write(iout,'(8x,''born(:,2,'',i3,'')='',3(f12.6,1x))',advance='no') ig,(bornq(j,2,i),j=1,3)
    write(iout,'('','')')
    write(iout,'(8x,''born(:,3,'',i3,'')='',3(f12.6,1x))',advance='no') ig,(bornq(j,3,i),j=1,3)
    write(iout,'('','')')
  enddo
  if (lusefc) then
    write(iout,'(8x,''scell(:)='',3i3)') (2*nd2cell(j)+1_i4,j=1,3)
  else
    write(iout,'(8x,''scell(:)='',3i3)') (nsuperghost(j,ncf),j=1,3)
  endif

  write(iout,'(''&end'')')
!-------------------------------------------------------------------------------
!  Parameters block
!-------------------------------------------------------------------------------
  write(iout,'(''&parameters'')')
!
!  Temperature
!
  write(iout,'(''        T='',f12.6)') max(temperature,1.0_dp)

  write(iout,'(''&end'')')
!-------------------------------------------------------------------------------
!  Flags block
!-------------------------------------------------------------------------------
  write(iout,'(''&flags'')')
!
!  Non-analytic correction
!
  if (index(keyword,'nono').eq.0) then
    write(iout,'(''        nonanalytic=.TRUE.'')')
  else
    write(iout,'(''        nonanalytic=.FALSE.'')')
  endif

  write(iout,'(''&end'')')
!
!  Close file
!
  close(iout)
!
!  Deallocate memory
!
  deallocate(nghostspecptr,stat=status)
  if (status/=0) call deallocate_error('outshengBTEcontrol','nghostspecptr')

  return
  end
!
  subroutine outshengBTEfc2(iout)
!
!  Outputs file FORCE_CONSTANTS_2ND for ShengBTE.
!
!   3/14 Created
!   7/15 Corrected for shell model so that only core-core component
!        is output
!   6/17 Parallel modifications added
!
!  Assumptions:
!
!  - cores are sorted to come before shells
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
!  Copyright Curtin University 2017
!
!  Julian Gale, CIC, Curtin University, June 2017
!
  use current
  use derivatives,    only : derv2
  use parallel
  use shells,         only : ncore, ncoptr
  implicit none
!
!  Passed variables
!
  integer(i4),                    intent(in)   :: iout
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ix
  integer(i4)                                  :: iy
  integer(i4)                                  :: iz
  integer(i4)                                  :: j
  integer(i4)                                  :: jx
  integer(i4)                                  :: jy
  integer(i4)                                  :: jz
!
!  Compute second order force constants (corrected for shells if necessary)
!
  if (nprocs.gt.1) then
    call secondorderfcd(.true.)
  else
    call secondorderfc(.true.)
  endif
!
!  Open file for the first time
!
  if (ioproc) then
    open(iout,file='FORCE_CONSTANTS_2ND',status='unknown',form='formatted')
!
!  Write out number of atoms
!
    write(iout,'(i6)') ncore
    close(iout)
  endif
  call mpbarrier
!
!  Loop over second derivative matrix force constant pairs
!
  ix = - 2
  iy = - 1
  iz =   0
  do i = 1,ncore
    if (atom2local(ncoptr(i)).gt.0) then
!
!  Re-open file
!
      open(iout,file='FORCE_CONSTANTS_2ND',position='append',status='old',form='formatted')
!
      ix = ix + 3
      iy = iy + 3
      iz = iz + 3
      jx = - 2
      jy = - 1
      jz =   0
      do j = 1,ncore
        jx = jx + 3
        jy = jy + 3
        jz = jz + 3
        write(iout,'(i6,1x,i6)') i,j
        write(iout,'(3f18.12)') derv2(jx,ix),derv2(jy,ix),derv2(jz,ix)
        write(iout,'(3f18.12)') derv2(jx,iy),derv2(jy,iy),derv2(jz,iy)
        write(iout,'(3f18.12)') derv2(jx,iz),derv2(jy,iz),derv2(jz,iz)
      enddo
      close(iout)
    endif
    call mpbarrier
  enddo
!
  return
  end
!
  subroutine outshengBTEfc3(iout,nfc3)
!
!  Outputs file FORCE_CONSTANTS_3RD for ShengBTE.
!
!   3/14 Created
!   6/17 Parallel modification made
!   6/17 Trap to prevent writing when nfc3=0 added
!   1/19 Changes for maxlinelength made
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
  use current
  use gulp_lengths
  use parallel
  implicit none
!
!  Passed variables
!
  integer(i4),                    intent(in)   :: iout      ! Channel number to read force constants from
  integer(i4),                    intent(in)   :: nfc3      ! Number of three-body force constants
!
!  Local variables
!
  character(len=maxlinelength)                 :: line
  integer(i4)                                  :: iout2
  logical                                      :: leof
!
!  Only ioproc needs to execute this routine
!
  if (.not.ioproc.or.nfc3.eq.0) return
!
!  Open file
!
  iout2 = iout + 1
  open(iout2,file='FORCE_CONSTANTS_3RD',status='unknown',form='formatted')
!
!  Write out number of third derivative terms
!
  write(iout2,'(i6)') nfc3
!
!  Copy remaining lines from temporary file to final one
!
  rewind(iout)
  leof = .false.
  do while (.not.leof)
    read(iout,'(a)',end=10,err=10) line
    write(iout2,'(a)') trim(line)
  enddo
10 continue
!
!  Close file
!
  close(iout2)

  return
  end
!
  subroutine outshengBTEfc2s(iout)
!
!  Outputs file FORCE_CONSTANTS_2ND for ShengBTE.
!  Supercell version.
!
!   3/15 Created from outshengBTEfc2
!   7/15 Corrected for shell model so that only core-core component
!        is output
!   6/17 Parallel modifications added
!
!  Assumptions:
!
!  - cores are sorted to come before shells
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
!  Copyright Curtin University 2017
!
!  Julian Gale, CIC, Curtin University, June 2017
!
  use current
  use derivatives,    only : d2cell, nd2cell, nd2cellptr, derv2
  use parallel
  use shells,         only : ncore, ncoptr
  implicit none
!
!  Passed variables
!
  integer(i4),                    intent(in)   :: iout
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: id
  integer(i4)                                  :: id2
  integer(i4)                                  :: ii
  integer(i4)                                  :: ii2
  integer(i4)                                  :: ij
  integer(i4)                                  :: ij2
  integer(i4)                                  :: is
  integer(i4)                                  :: ix
  integer(i4)                                  :: iy
  integer(i4)                                  :: iz
  integer(i4)                                  :: j
  integer(i4)                                  :: jd
  integer(i4)                                  :: jd2
  integer(i4)                                  :: ji
  integer(i4)                                  :: ji2
  integer(i4)                                  :: jj
  integer(i4)                                  :: jj2
  integer(i4)                                  :: js
  integer(i4)                                  :: jx
  integer(i4)                                  :: jy
  integer(i4)                                  :: jz
  integer(i4)                                  :: kd
  integer(i4)                                  :: kd2
  integer(i4)                                  :: ki
  integer(i4)                                  :: ki2
  integer(i4)                                  :: kj
  integer(i4)                                  :: kj2
  integer(i4)                                  :: kk
  integer(i4)                                  :: ncind
  real(dp)                                     :: d2i(3,3)
!
!  Open file
!
  if (ioproc) then
    open(iout,file='FORCE_CONSTANTS_2ND',status='unknown',form='formatted')
!
!  Write out number of atoms
!
    write(iout,'(i6)') ncore*(2*nd2cell(1)+1)*(2*nd2cell(2)+1)*(2*nd2cell(3)+1)
    close(iout)
  endif
  call mpbarrier
!
!  Loop over first atom
!
  ix = - 2
  iy = - 1
  iz =   0
  is =   0
  do i = 1,ncore
    ix = ix + 3
    iy = iy + 3
    iz = iz + 3
    if (atom2local(ncoptr(i)).gt.0) then
!
!  Re-open file
!
      open(iout,file='FORCE_CONSTANTS_2ND',position='append',status='old',form='formatted')
!
!  Find diagonal block less the d2cell self contributions
!
      d2i(1,1) = derv2(ix,ix)
      d2i(2,1) = derv2(iy,ix)
      d2i(3,1) = derv2(iz,ix)
      d2i(1,2) = derv2(ix,iy)
      d2i(2,2) = derv2(iy,iy)
      d2i(3,2) = derv2(iz,iy)
      d2i(1,3) = derv2(ix,iz)
      d2i(2,3) = derv2(iy,iz)
      d2i(3,3) = derv2(iz,iz)
      ncind = 0
      do ii = -nd2cell(1),nd2cell(1)
        do jj = -nd2cell(1),nd2cell(1)
          do kk = -nd2cell(1),nd2cell(1)
            ncind = ncind + 1
            d2i(1,1) = d2i(1,1) - d2cell(ix,ix,ncind)
            d2i(2,1) = d2i(2,1) - d2cell(iy,ix,ncind)
            d2i(3,1) = d2i(3,1) - d2cell(iz,ix,ncind)
            d2i(1,2) = d2i(1,2) - d2cell(ix,iy,ncind)
            d2i(2,2) = d2i(2,2) - d2cell(iy,iy,ncind)
            d2i(3,2) = d2i(3,2) - d2cell(iz,iy,ncind)
            d2i(1,3) = d2i(1,3) - d2cell(ix,iz,ncind)
            d2i(2,3) = d2i(2,3) - d2cell(iy,iz,ncind)
            d2i(3,3) = d2i(3,3) - d2cell(iz,iz,ncind)
          enddo
        enddo
      enddo
!
      do ii = 1,2*nd2cell(1)+1
        if (ii.gt.nd2cell(1)+1) then
          ii2 = ii - 2*nd2cell(1) - 2
        else
          ii2 = ii - 1 
        endif
        do ji = 1,2*nd2cell(2)+1
          if (ji.gt.nd2cell(2)+1) then
            ji2 = ji - 2*nd2cell(2) - 2
          else
            ji2 = ji - 1
          endif
          do ki = 1,2*nd2cell(3)+1
            if (ki.gt.nd2cell(3)+1) then
              ki2 = ki - 2*nd2cell(3) - 2
            else
              ki2 = ki - 1
            endif
            is = is + 1
!
!  Loop over second atom
!
            jx = - 2
            jy = - 1
            jz =   0
            js =   0
            do j = 1,ncore
              jx = jx + 3
              jy = jy + 3
              jz = jz + 3
              do ij = 1,2*nd2cell(1)+1
                if (ij.gt.nd2cell(1)+1) then
                  ij2 = ij - 2*nd2cell(1) - 2
                else
                  ij2 = ij - 1
                endif
                do jj = 1,2*nd2cell(2)+1
                  if (jj.gt.nd2cell(2)+1) then
                    jj2 = jj - 2*nd2cell(2) - 2
                  else
                    jj2 = jj - 1
                  endif
                  do kj = 1,2*nd2cell(3)+1
                    if (kj.gt.nd2cell(3)+1) then
                      kj2 = kj - 2*nd2cell(3) - 2
                    else
                      kj2 = kj - 1
                    endif
                    js = js + 1
!
!  Compute cell index difference
!
                    id = ij2 - ii2
                    jd = jj2 - ji2
                    kd = kj2 - ki2
!
!  Wrap cell index difference back to allowed range
!
                    if (id.gt.nd2cell(1)) then
                      id2 = id - 2*nd2cell(1) - 1
                    elseif (id.lt.-nd2cell(1)) then
                      id2 = id + 2*nd2cell(1) + 1
                    else
                      id2 = id 
                    endif
                    if (jd.gt.nd2cell(2)) then
                      jd2 = jd - 2*nd2cell(2) - 1
                    elseif (jd.lt.-nd2cell(2)) then
                      jd2 = jd + 2*nd2cell(2) + 1
                    else
                      jd2 = jd 
                    endif
                    if (kd.gt.nd2cell(3)) then
                      kd2 = kd - 2*nd2cell(3) - 1
                    elseif (kd.lt.-nd2cell(3)) then
                      kd2 = kd + 2*nd2cell(3) + 1
                    else
                      kd2 = kd 
                    endif
                    if (is.eq.js) then
!
!  Special case - diagonal block
!
                      write(iout,'(i6,1x,i6)') is,js
                      write(iout,'(3f18.12)') d2i(1,1),d2i(2,1),d2i(3,1)
                      write(iout,'(3f18.12)') d2i(1,2),d2i(2,2),d2i(3,2)
                      write(iout,'(3f18.12)') d2i(1,3),d2i(2,3),d2i(3,3)
                    else
!
!  General case - off diagonal
!
                      ncind = nd2cellptr(nd2cell(1)+1+id2,nd2cell(2)+1+jd2,nd2cell(3)+1+kd2)
!
                      write(iout,'(i6,1x,i6)') is,js
                      write(iout,'(3f18.12)') d2cell(jx,ix,ncind),d2cell(jy,ix,ncind),d2cell(jz,ix,ncind)
                      write(iout,'(3f18.12)') d2cell(jx,iy,ncind),d2cell(jy,iy,ncind),d2cell(jz,iy,ncind)
                      write(iout,'(3f18.12)') d2cell(jx,iz,ncind),d2cell(jy,iz,ncind),d2cell(jz,iz,ncind)
                    endif
                  enddo
                enddo
              enddo
            enddo
!
!  End of loops over js
!
          enddo
        enddo
      enddo
      close(iout)
    endif
    call mpbarrier
  enddo
!
  return
  end
