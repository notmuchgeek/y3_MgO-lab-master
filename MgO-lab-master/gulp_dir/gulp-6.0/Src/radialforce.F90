  subroutine radialforce(eradial,lgrad1,lgrad2)
!
!  Subroutine for calculating the energy and derivatives
!  due to a radial force.
!
!   3/07 Created
!   6/09 Site energy and virials added
!  11/09 Region derivatives added
!   4/12 xvir, yvir and zvir removed
!   2/17 Second derivatives added
!   2/17 Occupancy factor added
!   2/18 Trace added
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
  use configurations, only : nregionno
  use current
  use derivatives,    only : xdrv, ydrv, zdrv
  use derivatives,    only : xregdrv, yregdrv, zregdrv
  use derivatives,    only : derv2d
  use energies,       only : siteenergy
  use optimisation,   only : lopf, lfreeze
  use parallel,       only : nprocs, procid, natomsonnode, node2atom
  use radial
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  logical,  intent(in)  :: lgrad1
  logical,  intent(in)  :: lgrad2
  real(dp), intent(out) :: eradial
!
!  Local variables
!
  integer(i4)           :: i
  integer(i4)           :: iloc
  integer(i4)           :: ix
  integer(i4)           :: iy
  integer(i4)           :: iz
  integer(i4)           :: nregioni
  logical               :: lopi
  real(dp)              :: oci
  real(dp)              :: k
  real(dp)              :: x
  real(dp)              :: y
  real(dp)              :: z
!
!  Initialise energy
!
  eradial = 0.0_dp
!
!  If number of atoms is zero then return as there is nothing to do
!
  if (numat.eq.0) return
#ifdef TRACE
  call trace_in('radialforce')
#endif
!
!  Loop over atoms calculating radial force
!
  k = radialKcfg(ncf)
  x = radialXYZcfg(1,ncf)
  y = radialXYZcfg(2,ncf)
  z = radialXYZcfg(3,ncf)
  do i = 1+procid,numat,nprocs
    oci = occuf(i)
    eradial = eradial + k*oci*(xclat(i) - x)**2
    eradial = eradial + k*oci*(yclat(i) - y)**2
    eradial = eradial + k*oci*(zclat(i) - z)**2
    siteenergy(i) = siteenergy(i) + 0.5_dp*k*oci*((xclat(i) - x)**2 + (yclat(i) - y)**2 + (zclat(i) - z)**2)
  enddo
  if (lgrad1) then
    do i = 1+procid,nasym,nprocs
      oci = occua(i)
      xdrv(i) = xdrv(i) + k*oci*(xclat(i) - x)
      ydrv(i) = ydrv(i) + k*oci*(yclat(i) - y)
      zdrv(i) = zdrv(i) + k*oci*(zclat(i) - z)
!
      nregioni = nregionno(nsft+i)
      xregdrv(nregioni) = xregdrv(nregioni) + k*oci*(xclat(i) - x)
      yregdrv(nregioni) = yregdrv(nregioni) + k*oci*(yclat(i) - y)
      zregdrv(nregioni) = zregdrv(nregioni) + k*oci*(zclat(i) - z)
    enddo
  endif
  if (lgrad2) then
    ix =  -2
    iy =  -1
    iz =   0
    if (nprocs.gt.1) then
      do iloc = 1,natomsonnode
        i = node2atom(iloc)
!
!  Set attributes of i
!
        oci = occuf(i)
        lopi = (.not.lfreeze.or.lopf(nrelf2a(i)))
        if (lopi.or..not.lfreeze) then
          ix = ix + 3 
          iy = iy + 3 
          iz = iz + 3 
          derv2d(ix) = derv2d(ix) + k*oci
          derv2d(iy) = derv2d(iy) + k*oci
          derv2d(iz) = derv2d(iz) + k*oci
        endif
      enddo
    else
      do i = 1,nasym
!
!  Set attributes of i
!
        oci = occua(i)
        lopi = (.not.lfreeze.or.lopf(nrelf2a(i)))
        if (lopi.or..not.lfreeze) then
          ix = ix + 3
          iy = iy + 3
          iz = iz + 3
          derv2d(ix) = derv2d(ix) + k*oci
          derv2d(iy) = derv2d(iy) + k*oci
          derv2d(iz) = derv2d(iz) + k*oci
        endif
      enddo
    endif
  endif
!
!  Multiply by factor of a half
!
  eradial = 0.5_dp*eradial
#ifdef TRACE
  call trace_out('radialforce')
#endif
!
  return
  end
