  subroutine outinertia(iout,lappend)
!
!  Write out file containing information on momentum of 
!  inertia for molecules.
!
!  lappend = if .true. then structure should be appended
!            to file (for optimisation sequence)
!
!  12/14 Created from outxyz
!   1/15 Output eigenvectors and restrict to region 1
!   6/18 Trap added for parallel I/O
!   1/19 maxwordlength changes added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
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
  use configurations
  use g_constants,     only : avogadro
  use control
  use current
  use element
  use gulp_files
  use gulp_lengths
  use iochannels
  use molecule
  use parallel,        only : ioproc
  use species,         only : massspec
  implicit none
!
!  Passed variables
!
  integer(i4),         intent(in) :: iout
  logical,             intent(in) :: lappend
!
!  Local variables
!
  character(len=1)                :: numbers(10)
  character(len=maxwordlength)    :: inertfile2
  character(len=maxwordlength)    :: inertfilel
  integer(i4)                     :: i
  integer(i4)                     :: ifail
  integer(i4)                     :: ii
  integer(i4)                     :: incf
  integer(i4)                     :: ind
  integer(i4)                     :: j
  integer(i4)                     :: jj
  integer(i4)                     :: kk
  integer(i4)                     :: n
  real(dp)                        :: dx
  real(dp)                        :: dy
  real(dp)                        :: dz
  real(dp)                        :: inertia(3,3)
  real(dp)                        :: mi
  real(dp)                        :: totalmass
  real(dp)                        :: trace(3)
  real(dp)                        :: wrk(9)
  real(dp)                        :: xcom
  real(dp)                        :: ycom
  real(dp)                        :: zcom
  real(dp)                        :: xi
  real(dp)                        :: yi
  real(dp)                        :: zi
!
  data numbers/'0','1','2','3','4','5','6','7','8','9'/
!
!  Ensure that only a single node writes this file
!
  if (.not.ioproc) return
!***************************
!  Initialisation of file  *
!***************************
  if (.not.lappend) then
    inertfile2 = ' '
    inertfilel = inertfile
!
!  If inertia file name has been given then open file
!
    if (inertfilel(1:1).eq.' ') then
      inertfilel = 'gulp.ine'
    endif
    ind = index(inertfilel,'.ine')
    if (ind.eq.0) then
      call endstring(inertfilel,len(inertfilel),ind)
      inertfilel(ind:ind+3) = '.ine'
    endif
    inertfile2(1:ind-1) = inertfilel(1:ind-1)
    if (ncf.ge.100) then
      inertfile2(ind+4:ind+7) = inertfilel(ind:ind+3)
      inertfile2(ind:ind) = '_'
      incf = ncf
      ii = incf/10
      incf = incf - ii*10
      jj = incf/10
      incf = incf - jj*10
      inertfile2(ind+1:ind+1) = numbers(ii+1)
      inertfile2(ind+2:ind+2) = numbers(jj+1)
      inertfile2(ind+3:ind+3) = numbers(incf+1)
    elseif (ncf.ge.10) then
      inertfile2(ind+3:ind+6) = inertfilel(ind:ind+3)
      inertfile2(ind:ind) = '_'
      incf = ncf
      ii = incf/10
      incf = incf - ii*10
      inertfile2(ind+1:ind+1) = numbers(ii+1)
      inertfile2(ind+2:ind+2) = numbers(incf+1)
    elseif (ncfg.gt.1) then
      inertfile2(ind+2:ind+5) = inertfilel(ind:ind+3)
      inertfile2(ind:ind) = '_'
      inertfile2(ind+1:ind+1) = numbers(ncf+1)
    else
      inertfile2(ind:ind+3) = inertfilel(ind:ind+3)
    endif
    open(iout,file=inertfile2,status='unknown')
!
!  Write header
!
    write(iout,'(''# '')')
    write(iout,'(''#  Moment of inertia for molecules:'')')
    write(iout,'(''# '')')
    write(iout,'(''#  Mol No / Value  /         Eigenvector components '')')
  endif
!***********************
!  Configuration dump  *
!***********************
  do n = 1,nmol
!
!  Restrict to region 1
!
    if (nregionno(nsft+nrelf2a(nmollist(nmolptr(n)+1))).eq.1) then
!
!  Find centre of mass and total mass
!
      xcom = 0.0_dp
      ycom = 0.0_dp
      zcom = 0.0_dp
      totalmass = 0.0_dp
      do j = 1,nmolatom(n)
        i = nmollist(nmolptr(n)+j)
        mi = occuf(i)*massspec(nspecptr(i))
        call mindtoijk(nmolind(i),ii,jj,kk)
        xi = xclat(i) + dble(ii)*rv(1,1) + dble(jj)*rv(2,1) + dble(kk)*rv(3,1)
        yi = yclat(i) + dble(ii)*rv(1,2) + dble(jj)*rv(2,2) + dble(kk)*rv(3,2)
        zi = zclat(i) + dble(ii)*rv(1,3) + dble(jj)*rv(2,3) + dble(kk)*rv(3,3)
        xcom = xcom + xi*mi
        ycom = ycom + yi*mi
        zcom = zcom + zi*mi
        totalmass = totalmass + mi
      enddo
!
!  Initialise tensor
!
      inertia(1:3,1:3) = 0.0_dp
!
      xcom = xcom/totalmass
      ycom = ycom/totalmass
      zcom = zcom/totalmass
!
!  Calculate tensor
!
      do j = 1,nmolatom(n)
        i = nmollist(nmolptr(n)+j)
        mi = occuf(i)*massspec(nspecptr(i))
        call mindtoijk(nmolind(i),ii,jj,kk)
        xi = xclat(i) + ii*rv(1,1) + jj*rv(2,1) + kk*rv(3,1)
        yi = yclat(i) + ii*rv(1,2) + jj*rv(2,2) + kk*rv(3,2)
        zi = zclat(i) + ii*rv(1,3) + jj*rv(2,3) + kk*rv(3,3)
        dx = xi - xcom
        dy = yi - ycom
        dz = zi - zcom
        inertia(1,1) = inertia(1,1) + mi*dy*dy + mi*dz*dz
        inertia(2,1) = inertia(2,1) - mi*dx*dy
        inertia(3,1) = inertia(3,1) - mi*dx*dz
        inertia(1,2) = inertia(1,2) - mi*dy*dx
        inertia(2,2) = inertia(2,2) + mi*dx*dx + mi*dz*dz
        inertia(3,2) = inertia(3,2) - mi*dy*dz
        inertia(1,3) = inertia(1,3) - mi*dz*dx
        inertia(2,3) = inertia(2,3) - mi*dz*dy
        inertia(3,3) = inertia(3,3) + mi*dx*dx + mi*dy*dy
      enddo
!
!  Convert units to 10**-46 x kg.m**2
!
      inertia(1:3,1:3) = inertia(1:3,1:3)*1.0d23/avogadro
!
!  Calculate diagonalised tensor components
!
      call dsyev('V','U',3_i4,inertia,3_i4,trace,wrk,9_i4,ifail)
!
!  Write inertia information to file
!
      write(iout,'(i6,4(1x,f12.6))') n,trace(1),(inertia(j,1),j=1,3)
      write(iout,'(6x,4(1x,f12.6))')   trace(2),(inertia(j,2),j=1,3)
      write(iout,'(6x,4(1x,f12.6))')   trace(3),(inertia(j,3),j=1,3)
    endif
  enddo
!
  return
  end
