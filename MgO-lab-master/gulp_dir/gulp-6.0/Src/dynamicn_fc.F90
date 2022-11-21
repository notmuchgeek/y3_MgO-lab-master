  subroutine dynamicn_fc
!
!  Subroutine for calculating the dynamical matrix by finite differences
!
!  11/14 Created from dynamicn and dynamic_fc
!  12/14 Changed so that displaced atom doesn't change main 
!        arrays.
!  12/17 FastFD algorithm added
!   2/18 Trace added
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
!  Julian Gale, CIC, Curtin University, Feb 2018
!
  use bondorderdata, only : nbopot
  use cellmultipole
  use configurations
  use g_constants
  use control
  use current
  use derivatives
  use four
  use general,       only : phondiff
  use iochannels
  use kspace
  use m_bondorder,   only : setbondorderneigh, unsetbondorderneigh
  use m_edip,        only : setedipneigh, unsetedipneigh
  use m_three
  use molecule
  use optimisation
  use parallel
  use polarise
  use potentialxyz
  use shells,        only : nbsptr
  use six
  use sutton
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  use two
  implicit none
!
!  Local variables
!
  integer(i4)                                        :: i
  integer(i4)                                        :: ix
  integer(i4)                                        :: iy
  integer(i4)                                        :: iz
  integer(i4)                                        :: m
  integer(i4)                                        :: matom
  integer(i4)                                        :: maxlim
  integer(i4)                                        :: mcrd
  integer(i4)                                        :: mint
  integer(i4)                                        :: mm
  integer(i4)                                        :: n
  integer(i4)                                        :: nbsi
  integer(i4)                                        :: status
  logical                                            :: lforward
  logical                                            :: lradial
  real(dp)                                           :: rstep
  real(dp)                                           :: step
  real(dp),    dimension(:,:,:), allocatable,   save :: d1cell
  real(dp)                                           :: vmatom(4)
#ifdef TRACE
  call trace_in('dynamicn_fc')
#endif
!****************************
!  Zero second derivatives  *
!****************************
  mint  = 3*numat
  maxlim = mint
  if (nbsmat.gt.0) maxlim = maxlim + nbsmat
!
  if (maxlim.gt.maxd2u) then
    maxd2u = maxlim
    call changemaxd2
  endif
  if (maxlim.gt.maxd2) then
    maxd2 = maxlim
    call changemaxd2
  endif
  d2cell(1:maxlim,1:maxlim,1:nd2cells) = 0.0_dp
!
!  Allocate memory
!
  allocate(d1cell(4,numat,nd2cells),stat=status)
  if (status/=0) call outofmemory('dynamicn_fc','d1cell')
!
!  Initialise energy setup
!
  if (lfastfd) then
    if (nbopot.gt.0) then
      call setbondorderneigh
    endif
    if (lEDIP) then
      call setedipneigh
    endif
  endif
!
  rstep = 0.5_dp/phondiff
!************************************************
!  Loop over coordinates for finite difference  *
!************************************************
  do m = 1,2*maxlim
!
!  Find the degree of freedom and whether this is the forward or backward step
!
    if (m.gt.maxlim) then
      lforward = .true.
      mm = m - maxlim
    else
      lforward = .false.
      mm = m
    endif
    if (mm.gt.mint) then
      lradial = .true.
      matom = nbsptr(mm - mint)
      mcrd  = 4
    else
      lradial = .false.
      matom = (mm-1)/3 + 1
      mcrd  = mm - 3*(matom - 1)
    endif
!
!  Select step
!
    if (lforward) then
      step = phondiff
    else
      step = - phondiff
    endif
!
!  Initialise variables for matom
!
    vmatom(1) = xclat(matom)
    vmatom(2) = yclat(matom)
    vmatom(3) = zclat(matom)
    vmatom(4) = radf(matom)
!
!  Alter coordinate
!
    if (lradial) then
      vmatom(4) = vmatom(4) + step
    else
      if (mcrd.eq.1) then
        vmatom(1) = vmatom(1) + step
      elseif (mcrd.eq.2) then
        vmatom(2) = vmatom(2) + step
      elseif (mcrd.eq.3) then
        vmatom(3) = vmatom(3) + step
      endif
    endif
!
!  Evaluate function and first derivatives
!
    call forcebycell(numat,d1cell,matom,vmatom)
!
!  Add gradients to dynamical array
!
    if (lforward) then
      do n = 1,nd2cells
        ix = -2
        iy = -1
        iz =  0
        do i = 1,numat
          ix = ix + 3
          iy = iy + 3
          iz = iz + 3
          d2cell(ix,mm,n) = d2cell(ix,mm,n) + rstep*d1cell(1,i,n)
          d2cell(iy,mm,n) = d2cell(iy,mm,n) + rstep*d1cell(2,i,n)
          d2cell(iz,mm,n) = d2cell(iz,mm,n) + rstep*d1cell(3,i,n)
        enddo
      enddo
      if (nbsmat.gt.0) then
        do n = 1,nd2cells
          do nbsi = 1,nbsmat
            i = nbsptr(nbsi)
            d2cell(mint+nbsi,mm,n) = d2cell(mint+nbsi,mm,n) + rstep*d1cell(4,i,n)
          enddo
        enddo
      endif
    else
      do n = 1,nd2cells
        ix = -2
        iy = -1
        iz =  0
        do i = 1,numat
          ix = ix + 3
          iy = iy + 3
          iz = iz + 3
          d2cell(ix,mm,n) = d2cell(ix,mm,n) - rstep*d1cell(1,i,n)
          d2cell(iy,mm,n) = d2cell(iy,mm,n) - rstep*d1cell(2,i,n)
          d2cell(iz,mm,n) = d2cell(iz,mm,n) - rstep*d1cell(3,i,n)
        enddo
      enddo
      if (nbsmat.gt.0) then
        do n = 1,nd2cells
          do nbsi = 1,nbsmat
            i = nbsptr(nbsi)
            d2cell(mint+nbsi,mm,n) = d2cell(mint+nbsi,mm,n) - rstep*d1cell(4,i,n)
          enddo
        enddo
      endif
    endif
!*************************************
!  End loop over finite differences  *
!*************************************
  enddo
!
!  Finalise energy terms
!
  if (lfastfd) then
    if (lEDIP) then
      call unsetedipneigh
    endif
    if (nbopot.gt.0) then
      call unsetbondorderneigh
    endif
  endif
!
!  Zero on diagonal blocks for central cell
!
  ix = -2
  iy = -1
  iz =  0
  do i = 1,numat
    ix = ix + 3
    iy = iy + 3
    iz = iz + 3
    d2cell(ix,ix,nd2central) = 0.0_dp
    d2cell(iy,ix,nd2central) = 0.0_dp
    d2cell(iz,ix,nd2central) = 0.0_dp
    d2cell(ix,iy,nd2central) = 0.0_dp
    d2cell(iy,iy,nd2central) = 0.0_dp
    d2cell(iz,iy,nd2central) = 0.0_dp
    d2cell(ix,iz,nd2central) = 0.0_dp
    d2cell(iy,iz,nd2central) = 0.0_dp
    d2cell(iz,iz,nd2central) = 0.0_dp
  enddo
!
!  Free local memory
!
  deallocate(d1cell,stat=status)
  if (status/=0) call deallocate_error('dynamicn_fc','d1cell')
#ifdef TRACE
  call trace_out('dynamicn_fc')
#endif
!
  return
  end
