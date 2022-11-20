  subroutine move2a1
!
!  Moves ions from region 2a to region 1 at the end of a defect
!  run using the relaxed positions.
!
!  12/07 Unused variables removed
!   3/13 Use of scratch files for defect information removed
!   3/13 xdis, ydis, zdis used directly from arrays rather than
!        read from channel 48 since arrays shouldn't be modified
!        since defopt
!   3/13 After move then reg1 radius is now changed to reflect 
!        the new size.
!   3/13 Displacements not applied in region 2 as in previous 
!        version since this causes problems with restarting. 
!   5/17 Calling of changemaxr1at changed to ensure that there is
!        nreg1 + 1 words in arrays.
!   2/18 Trace added
!   3/19 Change of idopt to idoptindex and idopttype
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
  use control
  use current
  use defects
  use optimisation
  use parallel
  use region2a
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ind
  integer(i4), dimension(:), allocatable       :: itmp
  integer(i4)                                  :: n2last
  integer(i4)                                  :: ndind
  integer(i4)                                  :: nreg1o
  integer(i4)                                  :: status
  real(dp)                                     :: r2
  real(dp)                                     :: r2last
  real(dp)                                     :: r2max
  real(dp)                                     :: small
  real(dp)                                     :: xa
  real(dp)                                     :: ya
  real(dp)                                     :: za
  real(dp)                                     :: xi
  real(dp)                                     :: yi
  real(dp)                                     :: zi
#ifdef TRACE
  call trace_in('move2a1')
#endif
!
!  Local variables
!
  r2max = reg2a1(ncf)
  r2max = r2max*r2max
  small = 1.0d-6
  nreg1o = nreg1
!*******************************
!  Move region 2a to region 1  *
!*******************************
  r2last = reg1(ncf)*reg1(ncf)
  n2last = 0
  i = 0
  do while (r2last.le.r2max.and.i.lt.nreg2) 
    i = i + 1
    nreg1 = nreg1 + 1
    if (nreg1.ge.maxr1at) then
      maxr1at = nreg1 + 20
      call changemaxr1at
    endif
!
!  Transfer data between arrays
!
    xa = xr2a(i)
    ya = yr2a(i)
    za = zr2a(i)
    qdefe(nreg1) = qr2a(i)
    occdefe(nreg1) = or2a(i)
    radefe(nreg1) = rr2a(i)
    natdefe(nreg1) = nr2a(i)
    ntypdefe(nreg1) = ntr2a(i)
    ndefmol(nreg1) = nmr2a(i)
    if (ndefmol(nreg1).gt.0) then
      ndefind(nreg1) = nmir2a(i)
    else
      ndefind(nreg1) = 0
    endif
    ldefbsmat(nreg1) = ldbr2a(i)
    xdefe(nreg1) = xa
    ydefe(nreg1) = ya
    zdefe(nreg1) = za
    ndrel(nreg1) = nreg1
    ndrelop(nreg1) = 1
    npsite(nreg1) = nps(i)
    nreldef(nreg1) = nps(i)
!
!  Check distance to see if new shell has started
!
    xi = xa - xdc
    yi = ya - ydc
    zi = za - zdc
    r2 = xi*xi + yi*yi + zi*zi
    if ((r2 - r2last).gt.small) then
      n2last = nreg1
      r2last = r2
    endif
  enddo
  if (r2last.gt.r2max) then
!
!  Exited because of radius being reached
!
    nreg1 = n2last - 1
    reg1(ncf) = reg2a1(ncf)
  endif
!
!  Add optimisation flags
!
  do i = nreg1o + 1,nreg1
    ind = 3*(i - 1)
    idoptindex(nvar+1) = i
    idopttype(nvar+1) = idopt_dx
    idoptindex(nvar+2) = i
    idopttype(nvar+2) = idopt_dy
    idoptindex(nvar+3) = i
    idopttype(nvar+3) = idopt_dz
    nvar = nvar + 3
  enddo
!
!  Prior to calling sort1 we need to change the format of idoptindex
!
  allocate(itmp(4*nreg1),stat=status)
  if (status/=0) call outofmemory('move2a1','itmp')
!
  itmp(1:nvar) = idoptindex(1:nvar)
  do i = 1,4*nreg1
    idoptindex(i) = 0
  enddo 
  do i = 1,nvar
    if (idopttype(i).eq.idopt_dx) then
      idoptindex(3*itmp(i)-2) = 1
    elseif (idopttype(i).eq.idopt_dy) then
      idoptindex(3*itmp(i)-1) = 1
    elseif (idopttype(i).eq.idopt_dz) then
      idoptindex(3*itmp(i)) = 1
    elseif (idopttype(i).eq.idopt_dradius) then
      idoptindex(3*nreg1+itmp(i)) = 1
    endif
  enddo
!
!  Sort data to correct order
!
  call sort1(.false.,.false.)
!
!  Return idopt format to previous one
!
  itmp(1:4*nreg1) = idoptindex(1:4*nreg1)
  nvar = 0
  do i = 1,nreg1
    if (itmp(3*i-2).eq.1) then
      nvar = nvar + 1
      idoptindex(nvar) = i
      idopttype(nvar) = idopt_dx
    endif
    if (itmp(3*i-1).eq.1) then
      nvar = nvar + 1
      idoptindex(nvar) = i
      idopttype(nvar) = idopt_dy
    endif
    if (itmp(3*i).eq.1) then
      nvar = nvar + 1
      idoptindex(nvar) = i
      idopttype(nvar) = idopt_dz
    endif
  enddo
  do i = 1,nreg1
    if (itmp(3*nreg1+i).eq.1) then
      nvar = nvar + 1
      idoptindex(nvar) = i
      idopttype(nvar) = idopt_dradius
    endif
  enddo
  deallocate(itmp,stat=status)
  if (status/=0) call deallocate_error('move2a1','itmp')
!
!  Find start of configuration in total arrays
!
  ndind = 0
  do i = 1,ncf-1
    ndind = ndind + nreg1cfg(i)
  enddo
!
!  Check dimension of configuration arrays
!
  if (nreg1tot+nreg1-nreg1o.gt.maxtotr1at) then
    maxtotr1at = nreg1tot + nreg1 - nreg1o
    call changemaxtotr1at
  endif
!
!  Make space for expanded region 1
!
  do i = nreg1tot,ndind+nreg1o+1,-1
    natdefecfg(i+nreg1-nreg1o) = natdefecfg(i)
    ntypdefecfg(i+nreg1-nreg1o) = ntypdefecfg(i)
    nreldefcfg(i+nreg1-nreg1o) = nreldefcfg(i)
    xdefecfg(i+nreg1-nreg1o) = xdefecfg(i)
    ydefecfg(i+nreg1-nreg1o) = ydefecfg(i)
    zdefecfg(i+nreg1-nreg1o) = zdefecfg(i)
    qdefecfg(i+nreg1-nreg1o) = qdefecfg(i)
    occdefecfg(i+nreg1-nreg1o) = occdefecfg(i)
    radefecfg(i+nreg1-nreg1o) = radefecfg(i)
    ndefmolcfg(i+nreg1-nreg1o) = ndefmolcfg(i) 
    ndefindcfg(i+nreg1-nreg1o) = ndefindcfg(i)
    ldefbsmatcfg(i+nreg1-nreg1o) = ldefbsmatcfg(i)
    ldqmatomcfg(i+nreg1-nreg1o) = ldqmatomcfg(i)
    idoptcfg(1,i+nreg1-nreg1o) = idoptcfg(1,i)
    idoptcfg(2,i+nreg1-nreg1o) = idoptcfg(2,i)
    idoptcfg(3,i+nreg1-nreg1o) = idoptcfg(3,i)
  enddo
!
!  Insert new information
!
  do i = 1,nreg1
    natdefecfg(ndind+i) = natdefe(i)
    ntypdefecfg(ndind+i) = ntypdefe(i)
    xdefecfg(ndind+i) = xdefe(i)
    ydefecfg(ndind+i) = ydefe(i)
    zdefecfg(ndind+i) = zdefe(i)
    qdefecfg(ndind+i) = qdefe(i)
    occdefecfg(ndind+i) = occdefe(i)
    radefecfg(ndind+i) = radefe(i)
    ndefmolcfg(ndind+i) = ndefmol(i)
    ndefindcfg(ndind+i) = ndefind(i)
    ldefbsmatcfg(ndind+i) = ldefbsmat(i)
    ldqmatomcfg(ndind+i) = ldqmatom(i)
    if (i.gt.nreg1o) then
      idoptcfg(1,ndind+i) = 1
      idoptcfg(2,ndind+i) = 1
      idoptcfg(3,ndind+i) = 1
      nreldefcfg(ndind+i) = nreldef(i)
    else
      nreldefcfg(ndind+i) = npsite(i)
    endif
  enddo
!
!  Set final quantities for new region 1 size
!
  nreg1cfg(ncf) = nreg1
  nreg1tot = nreg1tot + nreg1 - nreg1o
  reg1(ncf) = min(reg2a1(ncf),reg2(ncf))
!
!  Set move radius to zero as move has been performed
!
  reg2a1(ncf) = 0.0_dp
!
!  Final task having changed the region size is to adjust the deflist pointer in ndptr.
!  To do this we have to pretend to create the new perfect region 1 otherwise the 
!  vacancy list will be incorrect as the atom numbers change.
!
#ifdef TRACE
  call trace_out('move2a1')
#endif
!
  return
  end
