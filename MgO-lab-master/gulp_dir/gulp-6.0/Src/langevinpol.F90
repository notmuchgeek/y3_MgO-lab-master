  subroutine langevinpol(Es,rmus,dpols,Ep,dEpdEs,d2EpdEs2,lgrad1,lgrad2)
!
!  Subroutine for calculating the Langevin-damped local dipole polarisation energy
!
!   2/20 Created
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
  use datatypes
  implicit none
!
!  Passed variables
!
  logical,                      intent(in)    :: lgrad1          ! If true then compute first derivative
  logical,                      intent(in)    :: lgrad2          ! If true then compute second derivative
  real(dp),                     intent(in)    :: Es              ! Modulus of the electric field 
  real(dp),                     intent(in)    :: rmus            ! Saturation dipole
  real(dp),                     intent(in)    :: dpols           ! Dipolar polarisability
  real(dp),                     intent(out)   :: Ep              ! Dipolar polarisation energy with Langevin damping
  real(dp),                     intent(out)   :: dEpdEs          ! First derivative of Ep with respect to modulus of electric field 
  real(dp),                     intent(out)   :: d2EpdEs2        ! Second derivative of Ep with respect to modulus of electric field
!
!  Local variables
!
  real(dp)                                    :: const
  real(dp)                                    :: dxdEs
  real(dp)                                    :: ex
  real(dp)                                    :: exm
  real(dp)                                    :: rex
  real(dp)                                    :: drexdx
  real(dp)                                    :: dsinhxdx
  real(dp)                                    :: sinhx
  real(dp)                                    :: rsinhx
  real(dp)                                    :: x
!
!  Trap small maximum dipole
!
  if (rmus.lt.1.0d-12) then
    Ep = 0.0_dp
    dEpdEs = 0.0_dp
    d2EpdEs2 = 0.0_dp
    return
  endif
!
!  Set the argument for the Langevin function, x
!
  dxdEs = 3.0_dp*dpols/rmus
  x = dxdEs*Es
!
!  Initialise return values of derivatives
!
  dEpdEs = 0.0_dp
  d2EpdEs2 = 0.0_dp
!
!  Trap zero polarisability
!
  if (dpols.lt.1.0d-8) return
!
!  Set up constant
!
  const = (rmus**2)/(3.0_dp*dpols)
!
  if (x.lt.1.0d-3) then
!
!  Handle small values of x using series expansion of sinh(x)/x
!
    sinhx = 1.0_dp + (x**2)/6.0_dp + (x**4)/120.0_dp
    Ep = 2.0_dp*const*log(sinhx)
    if (lgrad1) then
      rsinhx = 1.0_dp/sinhx
      dsinhxdx = x*(1.0_dp/3.0_dp + (x**2)/30.0_dp)
      dEpdEs = const*dxdEs*rsinhx*dsinhxdx/Es
      if (lgrad2) then
        d2EpdEs2 = (const*dxdEs*dxdEs*rsinhx*((1.0_dp/3.0_dp + (x**2)/10.0_dp) - rsinhx*dsinhxdx) - dEpdEs)/Es**2
      endif
    endif
  else
    ex = exp(x)
    exm = 1.0_dp/ex
    sinhx = 0.5_dp*(ex - exm)/x
    Ep = 2.0_dp*const*log(sinhx)
    if (lgrad1) then
      rex = (ex + exm)/(ex - exm)
      dEpdEs = const*dxdEs*(rex - 1.0_dp/x)/Es
      if (lgrad2) then
        drexdx = 1.0_dp - rex**2
        d2EpdEs2 = (const*dxdEs*dxdEs*(drexdx + 1.0_dp/x**2) - dEpdEs)/Es**2
      endif
    endif
  endif
!
  return
  end
