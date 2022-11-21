  subroutine outgaopt
!
!  Output final configurations from genetic algorithm optimisation
!
!  11/93 First created
!   6/95 Modified for additive constraints
!  12/07 Parallel error in setting nga fixed
!   3/17 Modifications made to allow for new variable order in iopt
!   1/18 jbase for fconf and xconf separated
!   8/18 Adding 1 to strains 1-3 removed
!   3/19 x0 removed
!   3/19 iopt replaced by ioptindex and iopttype
!   3/19 Constraint handling now moved to subroutine
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
!  Julian Gale, CIC, Curtin University, March 2019
!
  use configurations
  use control,             only : lstraincell
  use current
  use gaconf
  use genetic
  use iochannels
  use optimisation
  use parallel
  implicit none
!
!  Local variables
!
  character(len=5)                                  :: namecell(6)
  character(len=1)                                  :: namecrd(3)
  integer(i4)                                       :: i
  integer(i4)                                       :: ip
  integer(i4)                                       :: j
  integer(i4)                                       :: jbase
  integer(i4)                                       :: jbasef
  integer(i4)                                       :: jmax
  integer(i4)                                       :: ng
  integer(i4)                                       :: nga
  integer(i4)                                       :: ngroup
  real(dp)                                          :: cellp(6,4)
  real(dp)                                          :: rvp(3,3)
!
  data namecell/'a    ','b    ','c    ','alpha','beta','gamma'/
  data namecrd/'x','y','z'/
!***************************
!  Decide type of results  *
!***************************
  if (ioproc) then
    write(ioout,'(/)')
    if (ngabest.eq.0) then
      write(ioout,'(''  Final configurations from genetic algoritthm optimisation :'',/)')
      nga = ngacfg
    else
      write(ioout,'(''  Best configurations from genetic algorithm optimisation :'',/)')
      nga = ngabest
    endif
  else
    if (ngabest.eq.0) then
      nga = ngacfg
    else
      nga = ngabest
    endif
  endif
!
!  How many lots of 4?
!
  ngroup = ((nga-1)/4) + 1
!***********************
!  Table of variables  *
!***********************
  do ng = 1,ngroup
    if (ng.eq.ngroup) then
      jmax = nga - 4*(ngroup-1)
    else
      jmax = 4
    endif
    jbase = 4*(ng - 1)
    jbasef = jbase
    if (ngabest.gt.0) jbasef = jbasef + mgacfg
    if (ioproc) then
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(25x,4(5x,i3,5x))')(4*(ng-1)+j,j=1,jmax)
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''  Sum of squares =  '',5x,4(1x,f12.5))') (fconf(jbasef+j),j=1,jmax)
      write(ioout,'(''--------------------------------------------------------------------------------'')')
    endif
    do i = 1,jmax
      do j = 1,3
        rv(1,j) = rvcfg(1,j,ncf)
        rv(2,j) = rvcfg(2,j,ncf)
        rv(3,j) = rvcfg(3,j,ncf)
      enddo
      if (lstraincell) then
        strain(1:nstrains) = straincfg(1:nstrains,ncf)
      else
        strain(1:nstrains) = 0.0_dp
      endif
!****************************************
!  Transfer variables to configuration  *
!****************************************
      call vartocfg(nvar,xconf(1,i+jbase))
!**********************
!  Apply constraints  *
!**********************
      if (ncon.gt.0) then
        call applyconstraints
      endif
      if (ncell.gt.0) then
!******************
!  Apply strains  *
!******************
        call strain3D(strain,rv)
!
!  Generate non-primitive cell
!
        do j = 1,3
          rvp(1,j) = rv(1,j)
          rvp(2,j) = rv(2,j)
          rvp(3,j) = rv(3,j)
        enddo
        if (ncbl.gt.1) call uncentre(rvp)
        call uncell3D(rvp,cellp(1,i),cellp(2,i),cellp(3,i),cellp(4,i),cellp(5,i),cellp(6,i))
      endif
    enddo
    if (ioproc) then
      do i = 1,nvar
        ip = ioptindex(i)
        if (iopttype(i).eq.iopt_cell) then
!
!  Cell parameter
!
          write(ioout,'(i3,2x,''Unit cell'',4x,a5,2x,4(1x,f12.5))') i,namecell(ip),(cellp(ip,j),j=1,jmax)
        elseif (iopttype(i).eq.iopt_strain) then
!
!  Cell strain
!
          write(ioout,'(i3,2x,''Unit strain'',2x,a5,2x,4(1x,f12.5))') i,namecell(ip),(cellp(ip,j),j=1,jmax)
        elseif (iopttype(i).eq.iopt_xf) then
!
!  x fractional coordinate
!
          write(ioout,'(i3,2x,''Fractional'',2x,i3,1x,a1,2x,4(1x,f12.5))') i,ip,namecrd(1),(xconf(i,j+jbase),j=1,jmax)
        elseif (iopttype(i).eq.iopt_yf) then
!
!  y fractional coordinate
!
          write(ioout,'(i3,2x,''Fractional'',2x,i3,1x,a1,2x,4(1x,f12.5))') i,ip,namecrd(2),(xconf(i,j+jbase),j=1,jmax)
        elseif (iopttype(i).eq.iopt_zf) then
!
!  z fractional coordinate
!
          write(ioout,'(i3,2x,''Fractional'',2x,i3,1x,a1,2x,4(1x,f12.5))') i,ip,namecrd(3),(xconf(i,j+jbase),j=1,jmax)
        elseif (iopttype(i).eq.iopt_radius) then
!
!  Radius
!
          write(ioout,'(i3,2x,''Radius    '',2x,i3,4x,4(1x,f12.5))') i,ip,(xconf(i,j+jbase),j=1,jmax)
        endif
      enddo
      write(ioout,'(''--------------------------------------------------------------------------------'',/)')
    endif
  enddo
!
  return
  end
