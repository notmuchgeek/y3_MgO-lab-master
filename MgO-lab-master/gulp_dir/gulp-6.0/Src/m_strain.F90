module m_strain

  use datatypes

  implicit none

  real(dp),                               save :: straindet               ! Inverse determinant of strain matrix
  real(dp),                               save :: strainddetds(6)         ! First derivatives of determinant of strain matrix
  real(dp),                               save :: straind2detds2(6,6)     ! Second derivatives of determinant of strain matrix
  real(dp),                               save :: strainfull(3,3)         ! Strain matrix in full form with unit matrix added to diagonal
  real(dp),                               save :: straininv(3,3)          ! Inverse strain matrix before division by determinant
  real(dp),                               save :: dstraininvds(6,3,3)     ! First derivatives of elements of straininv
  real(dp),                               save :: d2straininvds2(6,6,3,3) ! Second derivatives of elements of straininv
  real(dp),                               save :: straininverse(3,3)      ! Complete inversion of the strain matrix

  public :: twostrterms, gstrterms, realstrterms, real1strterm, vecprestrain, cartstrterm, gxyzterms

CONTAINS

  subroutine getinversestrain(strainin)
!
!  Computes the inverse of the strain matrix
!
!   9/18 Created
!   9/18 Extra setup added to avoid repeating in twostrgterms
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
!  Julian Gale, CIC, Curtin University, September 2018
!
  use datatypes
  use current
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  real(dp),   intent(in)    :: strainin(6)
!
!  Local variables
!
  integer(i4)               :: ifail
  real(dp)                  :: e(6)
  real(dp)                  :: wrk(6)
#ifdef TRACE
  call trace_in('getinversestrain')
#endif
  if (ndim.eq.3) then
!
!  Build the 3-D strain matrix
!
    strainfull(1,1) = 1.0_dp + strainin(1)
    strainfull(2,2) = 1.0_dp + strainin(2)
    strainfull(3,3) = 1.0_dp + strainin(3)
    strainfull(3,2) = 0.5_dp*strainin(4)
    strainfull(2,3) = 0.5_dp*strainin(4)
    strainfull(3,1) = 0.5_dp*strainin(5)
    strainfull(1,3) = 0.5_dp*strainin(5)
    strainfull(2,1) = 0.5_dp*strainin(6)
    strainfull(1,2) = 0.5_dp*strainin(6)
!
!  Compute the inverse of the strain matrix
!
    straininverse(1:3,1:3) = strainfull(1:3,1:3)
    call matrix_inversion(straininverse,3_i4,3_i4,wrk,ifail)
!
    e(1) = strainfull(1,1)
    e(2) = strainfull(2,2)
    e(3) = strainfull(3,3)
    e(4) = strainfull(3,2)
    e(5) = strainfull(3,1)
    e(6) = strainfull(2,1)
!
!  Compute determinant of inverse matrix
!
    straindet = e(1)*(e(2)*e(3) - e(4)**2) - &
                e(6)*(e(3)*e(6) - e(4)*e(5)) + &
                e(5)*(e(4)*e(6) - e(2)*e(5))
!
    straindet = 1.0_dp/straindet
!
!  Compute first strain derivatives of determinant
!
    strainddetds(1) = e(2)*e(3) - e(4)**2
    strainddetds(2) = e(1)*e(3) - e(5)**2
    strainddetds(3) = e(1)*e(2) - e(6)**2
    strainddetds(4) = e(5)*e(6) - e(1)*e(4)
    strainddetds(5) = e(4)*e(6) - e(2)*e(5)
    strainddetds(6) = e(4)*e(5) - e(3)*e(6)
!
!  Compute second strain derivatives of determinant
!
    straind2detds2(1:6,1:6) = 0.0_dp
!
    straind2detds2(2,1) = e(3)
    straind2detds2(3,1) = e(2)
    straind2detds2(4,1) = - e(4)
!
    straind2detds2(1,2) = e(3)
    straind2detds2(3,2) = e(1)
    straind2detds2(5,2) = - e(5)
!
    straind2detds2(1,3) = e(2)
    straind2detds2(2,3) = e(1)
    straind2detds2(6,3) = - e(6)
!
    straind2detds2(1,4) = - e(4)
    straind2detds2(4,4) = - 0.5_dp*e(1)
    straind2detds2(5,4) = 0.5_dp*e(6)
    straind2detds2(6,4) = 0.5_dp*e(5)
!
    straind2detds2(2,5) = - e(5)
    straind2detds2(5,5) = - 0.5_dp*e(2)
    straind2detds2(4,5) = 0.5_dp*e(6)
    straind2detds2(6,5) = 0.5_dp*e(4)
!
    straind2detds2(3,6) = - e(6)
    straind2detds2(6,6) = - 0.5_dp*e(3)
    straind2detds2(4,6) = 0.5_dp*e(5)
    straind2detds2(5,6) = 0.5_dp*e(4)
  elseif (ndim.eq.2) then
!
!  Build the 2-D strain matrix
!
    strainfull(1,1) = 1.0_dp + strainin(1)
    strainfull(2,2) = 1.0_dp + strainin(2)
    strainfull(2,1) = 0.5_dp*strainin(3)
    strainfull(1,2) = 0.5_dp*strainin(3)
!
!  Compute the inverse of the strain matrix
!
    straininverse(1:2,1:2) = strainfull(1:2,1:2)
    call matrix_inversion(straininverse,3_i4,2_i4,wrk,ifail)
!
!  Compute determinant of inverse matrix
!
    straindet = 1.0_dp/(strainfull(1,1)*strainfull(2,2) - strainfull(2,1)**2)
!
!  Compute strain derivatives of determinant
!
    strainddetds(1) = strainfull(2,2)
    strainddetds(2) = strainfull(1,1)
    strainddetds(3) = - strainfull(2,1)
!
!  Compute second strain derivatives of determinant
!
    straind2detds2(1:6,1:6) = 0.0_dp
    straind2detds2(2,1) = 1.0_dp
    straind2detds2(1,2) = 1.0_dp
    straind2detds2(3,3) = - 0.5_dp
  elseif (ndim.eq.1) then
!
!  Compute the inverse 1-D strain matrix
!
    strainfull(1,1) = 1.0_dp + strainin(1)
    straininverse(1,1) = 1.0_dp/strainfull(1,1)
    straindet = strainfull(1,1)
    strainddetds(1) = 1.0_dp
    straind2detds2(1,1) = 0.0_dp
  endif
#ifdef TRACE
  call trace_out('getinversestrain')
#endif
!
  return
  end

  subroutine twostrterms(ndim,maxvec,nvec,xvec,yvec,zvec,xcom,ycom,zcom,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,lgrad2)
!
!  Subroutine for calculating two-body strain terms
!
!  On entry :
!
!   ndim              = number of dimensions
!   nvec              = number of vectors to be processed
!   maxvec            = maximum number of vectors (lower dimension of d2r2dx2, etc)
!   xvec, yvec, zvec  = Cartesian components of interatomic vectors
!   xcom, ycom, zcom  = Cartesian components of intramolecular vectors to the centre of mass
!   lgrad2            = if true then compute d2r2ds2 and d2r2dsdx
!
!  On exit :
!
!   d2r2dx2           = second Cartesian derivative component
!   dr2ds             = 1/2 x first derivative of r^2 with respect to strain (can be the same as d2r2dx2)
!   d2r2ds2           = 1/2 x second derivative of r^2 with respect to two strains
!   d2r2dsdx          = 1/2 x second derivative of r^2 with respect to strains and a Cartesian component - including cell strain
!  
!   8/18 Created from threestrterms
!   9/18 Second derivatives added to allow for finite strain
!  11/18 Finite strain flag introduced instead of lstraincell
!   5/19 Second derivatives with respect to strain now complete - no need for strfin
!   7/19 Factor of two correction to the finite strain second derivatives with 
!        respect to strain
!   4/20 Rigid molecule derivatives added
!   4/20 d2r2dsdc added
!   4/20 derv3c changes reversed as they are no longer required
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
!  Julian Gale, CIC, Curtin University, April 2020
!
  use datatypes
  use control,       only : lrigid
  use current,       only : strain
  use derivatives,   only : lfinitestrain
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)  :: ndim
  integer(i4), intent(in)  :: nvec
  integer(i4), intent(in)  :: maxvec
  logical,     intent(in)  :: lgrad2
  real(dp),    intent(out) :: d2r2dx2(maxvec,6)
  real(dp),    intent(out) :: dr2ds(maxvec,6)
  real(dp),    intent(out) :: d2r2ds2(maxvec,6,6)
  real(dp),    intent(out) :: d2r2dsdx(maxvec,6,3)
  real(dp),    intent(in)  :: xvec(nvec)
  real(dp),    intent(in)  :: yvec(nvec)
  real(dp),    intent(in)  :: zvec(nvec)
  real(dp),    intent(in)  :: xcom
  real(dp),    intent(in)  :: ycom
  real(dp),    intent(in)  :: zcom
!
!  Local variables
!
  integer(i4)              :: i
  integer(i4)              :: j
  integer(i4)              :: k
  integer(i4)              :: l
  real(dp)                 :: rpd(6)
  real(dp)                 :: tmp(3)
  real(dp)                 :: xyzloc(3)
#ifdef TRACE
  call trace_in('twostrterms')
#endif
!
!  Set up Cartesian second derivatives - this is the same for all dimensionalities & also needed for strains
!
  if (lrigid) then
    do k = 1,nvec
      d2r2dx2(k,1) = xvec(k)*(xvec(k) - xcom)
      d2r2dx2(k,2) = yvec(k)*(yvec(k) - ycom)
      d2r2dx2(k,3) = zvec(k)*(zvec(k) - zcom)
      d2r2dx2(k,4) = 0.5_dp*(yvec(k)*(zvec(k) - zcom) + zvec(k)*(yvec(k) - ycom))
      d2r2dx2(k,5) = 0.5_dp*(xvec(k)*(zvec(k) - zcom) + zvec(k)*(xvec(k) - xcom))
      d2r2dx2(k,6) = 0.5_dp*(xvec(k)*(yvec(k) - ycom) + yvec(k)*(xvec(k) - xcom))
    enddo
  else
    do k = 1,nvec
      d2r2dx2(k,1) = xvec(k)*xvec(k)
      d2r2dx2(k,2) = yvec(k)*yvec(k)
      d2r2dx2(k,3) = zvec(k)*zvec(k)
      d2r2dx2(k,4) = yvec(k)*zvec(k)
      d2r2dx2(k,5) = xvec(k)*zvec(k)
      d2r2dx2(k,6) = xvec(k)*yvec(k)
    enddo
  endif
!
  if (ndim.eq.3) then
    if (lfinitestrain) then
      do k = 1,nvec
!
!  Create vector scaled by inverse strain matrix
!
        xyzloc(1) = straininverse(1,1)*(xvec(k) - xcom) + &
                    straininverse(2,1)*(yvec(k) - ycom) + &
                    straininverse(3,1)*(zvec(k) - zcom)
        xyzloc(2) = straininverse(1,2)*(xvec(k) - xcom) + &
                    straininverse(2,2)*(yvec(k) - ycom) + &
                    straininverse(3,2)*(zvec(k) - zcom)
        xyzloc(3) = straininverse(1,3)*(xvec(k) - xcom) + &
                    straininverse(2,3)*(yvec(k) - ycom) + &
                    straininverse(3,3)*(zvec(k) - zcom)
!
        rpd(1) = xyzloc(1)*xyzloc(1)
        rpd(2) = xyzloc(2)*xyzloc(2)
        rpd(3) = xyzloc(3)*xyzloc(3)
        rpd(4) = xyzloc(2)*xyzloc(3)
        rpd(5) = xyzloc(1)*xyzloc(3)
        rpd(6) = xyzloc(1)*xyzloc(2)
!
        dr2ds(k,1) = (1.0_dp + strain(1))*rpd(1) + &
                      0.5_dp*(strain(6)*rpd(6) + strain(5)*rpd(5)) + xyzloc(1)*xcom
        dr2ds(k,2) = (1.0_dp + strain(2))*rpd(2) + &
                      0.5_dp*(strain(6)*rpd(6) + strain(4)*rpd(4)) + xyzloc(2)*ycom
        dr2ds(k,3) = (1.0_dp + strain(3))*rpd(3) + &
                      0.5_dp*(strain(5)*rpd(5) + strain(4)*rpd(4)) + xyzloc(3)*zcom
        dr2ds(k,4) = (1.0_dp + 0.5_dp*(strain(2) + strain(3)))*rpd(4) + &
                      0.25_dp*(strain(4)*(rpd(2) + rpd(3)) + &
                               strain(5)*rpd(6) + strain(6)*rpd(5)) + &
                      0.5_dp*(xyzloc(2)*zcom + xyzloc(3)*ycom)
        dr2ds(k,5) = (1.0_dp + 0.5_dp*(strain(1) + strain(3)))*rpd(5) + &
                      0.25_dp*(strain(5)*(rpd(1) + rpd(3)) + &
                               strain(4)*rpd(6) + strain(6)*rpd(4)) + &
                      0.5_dp*(xyzloc(1)*zcom + xyzloc(3)*xcom)
        dr2ds(k,6) = (1.0_dp + 0.5_dp*(strain(1) + strain(2)))*rpd(6) + &
                      0.25_dp*(strain(6)*(rpd(1) + rpd(2)) + &
                               strain(4)*rpd(5) + strain(5)*rpd(4)) + &
                      0.5_dp*(xyzloc(1)*ycom + xyzloc(2)*xcom)
        if (lgrad2) then
!
!  Initialise second derivatives
!
          d2r2ds2(k,1:6,1:6) = 0.0_dp
!
!  Second derivatives with respect to 2 strains
!
          d2r2ds2(k,1,1) = rpd(1)
          d2r2ds2(k,5,1) = 0.5_dp*rpd(5)
          d2r2ds2(k,6,1) = 0.5_dp*rpd(6)
!
          d2r2ds2(k,2,2) = rpd(2)
          d2r2ds2(k,4,2) = 0.5_dp*rpd(4)
          d2r2ds2(k,6,2) = 0.5_dp*rpd(6)
!
          d2r2ds2(k,3,3) = rpd(3)
          d2r2ds2(k,4,3) = 0.5_dp*rpd(4)
          d2r2ds2(k,5,3) = 0.5_dp*rpd(5)
!
          d2r2ds2(k,4,4) = 0.25_dp*(rpd(2) + rpd(3))
          d2r2ds2(k,5,4) = 0.25_dp*rpd(6)
          d2r2ds2(k,6,4) = 0.25_dp*rpd(5)
!
          d2r2ds2(k,5,5) = 0.25_dp*(rpd(1) + rpd(3))
          d2r2ds2(k,6,5) = 0.25_dp*rpd(4)
!
          d2r2ds2(k,6,6) = 0.25_dp*(rpd(1) + rpd(2))
!
          do i = 2,6
            do j = 1,i-1
              d2r2ds2(k,j,i) = d2r2ds2(k,i,j)
            enddo
          enddo
!
!  Second derivatives with respect to 1 strain and a Cartesian coordinate
!
          d2r2dsdx(k,1,1) = 2.0_dp*xyzloc(1)*(1.0_dp + strain(1)) + &
                            0.5_dp*(strain(6)*xyzloc(2) + strain(5)*xyzloc(3)) + xcom
          d2r2dsdx(k,1,2) = 0.5_dp*strain(6)*xyzloc(1)
          d2r2dsdx(k,1,3) = 0.5_dp*strain(5)*xyzloc(1)
!
          d2r2dsdx(k,2,1) = 0.5_dp*strain(6)*xyzloc(2)
          d2r2dsdx(k,2,2) = 2.0_dp*xyzloc(2)*(1.0_dp + strain(2)) + &
                            0.5_dp*(strain(6)*xyzloc(1) + strain(4)*xyzloc(3)) + ycom
          d2r2dsdx(k,2,3) = 0.5_dp*strain(4)*xyzloc(2)
!
          d2r2dsdx(k,3,1) = 0.5_dp*strain(5)*xyzloc(3)
          d2r2dsdx(k,3,2) = 0.5_dp*strain(4)*xyzloc(3)
          d2r2dsdx(k,3,3) = 2.0_dp*xyzloc(3)*(1.0_dp + strain(3)) + &
                            0.5_dp*(strain(5)*xyzloc(1) + strain(4)*xyzloc(2)) + zcom
!
          d2r2dsdx(k,4,1) = 0.25_dp*(strain(5)*xyzloc(2) + strain(6)*xyzloc(3))
          d2r2dsdx(k,4,2) = (1.0_dp + 0.5_dp*(strain(2) + strain(3)))*xyzloc(3) + &
                             0.25_dp*(2.0_dp*strain(4)*xyzloc(2) + strain(5)*xyzloc(1)) + 0.5_dp*zcom
          d2r2dsdx(k,4,3) = (1.0_dp + 0.5_dp*(strain(2) + strain(3)))*xyzloc(2) + &
                             0.25_dp*(2.0_dp*strain(4)*xyzloc(3) + strain(6)*xyzloc(1)) + 0.5_dp*ycom
!
          d2r2dsdx(k,5,1) = (1.0_dp + 0.5_dp*(strain(1) + strain(3)))*xyzloc(3) + &
                             0.25_dp*(2.0_dp*strain(5)*xyzloc(1) + strain(4)*xyzloc(2)) + 0.5_dp*zcom
          d2r2dsdx(k,5,2) = 0.25_dp*(strain(4)*xyzloc(1) + strain(6)*xyzloc(3))
          d2r2dsdx(k,5,3) = (1.0_dp + 0.5_dp*(strain(1) + strain(3)))*xyzloc(1) + &
                             0.25_dp*(2.0_dp*strain(5)*xyzloc(3) + strain(6)*xyzloc(2)) + 0.5_dp*xcom
!
          d2r2dsdx(k,6,1) = (1.0_dp + 0.5_dp*(strain(1) + strain(2)))*xyzloc(2) + &
                             0.25_dp*(2.0_dp*strain(6)*xyzloc(1) + strain(4)*xyzloc(3)) + 0.5_dp*ycom
          d2r2dsdx(k,6,2) = (1.0_dp + 0.5_dp*(strain(1) + strain(2)))*xyzloc(1) + &
                             0.25_dp*(2.0_dp*strain(6)*xyzloc(2) + strain(5)*xyzloc(3)) + 0.5_dp*xcom
          d2r2dsdx(k,6,3) = 0.25_dp*(strain(4)*xyzloc(1) + strain(5)*xyzloc(2))
!
!  Transform mixed derivatives by inverse strain matrix
!
          do i = 1,6
            tmp(1:3) = 0.0_dp
            do j = 1,3
              do l = 1,3
                tmp(j) = tmp(j) + d2r2dsdx(k,i,l)*straininverse(l,j)
              enddo
            enddo
            d2r2dsdx(k,i,1:3) = tmp(1:3)
          enddo
        endif
      enddo
    else
      dr2ds(1:nvec,1:6) = d2r2dx2(1:nvec,1:6)
      if (lgrad2) then
        do k = 1,nvec
          rpd(1) = xvec(k)*xvec(k)
          rpd(2) = yvec(k)*yvec(k)
          rpd(3) = zvec(k)*zvec(k)
          rpd(4) = yvec(k)*zvec(k)
          rpd(5) = xvec(k)*zvec(k)
          rpd(6) = xvec(k)*yvec(k)
!
!  Initialise second derivatives
!
          d2r2ds2(k,1:6,1:6) = 0.0_dp
          d2r2dsdx(k,1:6,1:3) = 0.0_dp
!
!  Second derivatives with respect to 2 strains
!
          d2r2ds2(k,1,1) = (xvec(k) - xcom)*((xvec(k) - xcom) + xvec(k))
          d2r2ds2(k,5,1) = 0.5_dp*((zvec(k) - zcom)*(xvec(k) - xcom) + xvec(k)*(zvec(k) - zcom))
          d2r2ds2(k,6,1) = 0.5_dp*((yvec(k) - ycom)*(xvec(k) - xcom) + xvec(k)*(yvec(k) - ycom))
!
          d2r2ds2(k,2,2) = (yvec(k) - ycom)*((yvec(k) - ycom) + yvec(k))
          d2r2ds2(k,4,2) = 0.5_dp*((zvec(k) - zcom)*(yvec(k) - ycom) + yvec(k)*(zvec(k) - zcom))
          d2r2ds2(k,6,2) = 0.5_dp*((xvec(k) - xcom)*(yvec(k) - ycom) + yvec(k)*(xvec(k) - xcom))
!
          d2r2ds2(k,3,3) = (zvec(k) - zcom)*((zvec(k) - zcom) + zvec(k))
          d2r2ds2(k,4,3) = 0.5_dp*((yvec(k) - ycom)*(zvec(k) - zcom) + zvec(k)*(yvec(k) - ycom))
          d2r2ds2(k,5,3) = 0.5_dp*((xvec(k) - xcom)*(zvec(k) - zcom) + zvec(k)*(xvec(k) - xcom))
!
          d2r2ds2(k,4,4) = 0.25_dp*((yvec(k) - ycom)**2 + (zvec(k) - zcom)**2 + &
                                    yvec(k)*(yvec(k) - ycom) + zvec(k)*(zvec(k) - zcom))
          d2r2ds2(k,5,4) = 0.25_dp*((xvec(k) - xcom)*(yvec(k) - ycom) + yvec(k)*(xvec(k) - xcom))
          d2r2ds2(k,6,4) = 0.25_dp*((xvec(k) - xcom)*(zvec(k) - zcom) + zvec(k)*(xvec(k) - xcom))
!
          d2r2ds2(k,5,5) = 0.25_dp*((xvec(k) - xcom)**2 + (zvec(k) - zcom)**2 + &
                                    xvec(k)*(xvec(k) - xcom) + zvec(k)*(zvec(k) - zcom))
          d2r2ds2(k,6,5) = 0.25_dp*((yvec(k) - ycom)*(zvec(k) - zcom) + zvec(k)*(yvec(k) - ycom))
!
          d2r2ds2(k,6,6) = 0.25_dp*((xvec(k) - xcom)**2 + (yvec(k) - ycom)**2 + &
                                    xvec(k)*(xvec(k) - xcom) + yvec(k)*(yvec(k) - ycom))
!
          do i = 2,6
            do j = 1,i-1
              d2r2ds2(k,j,i) = d2r2ds2(k,i,j)
            enddo
          enddo
!
!  Second derivatives with respect to 1 strain and a Cartesian coordinate for fractionals
!
          d2r2dsdx(k,1,1) = xvec(k) + (xvec(k) - xcom)
          d2r2dsdx(k,2,2) = yvec(k) + (yvec(k) - ycom)
          d2r2dsdx(k,3,3) = zvec(k) + (zvec(k) - zcom)
!
          d2r2dsdx(k,4,2) = 0.5_dp*(zvec(k) + (zvec(k) - zcom))
          d2r2dsdx(k,4,3) = 0.5_dp*(yvec(k) + (yvec(k) - ycom))
!
          d2r2dsdx(k,5,1) = 0.5_dp*(zvec(k) + (zvec(k) - zcom))
          d2r2dsdx(k,5,3) = 0.5_dp*(xvec(k) + (xvec(k) - xcom))
!
          d2r2dsdx(k,6,1) = 0.5_dp*(yvec(k) + (yvec(k) - ycom))
          d2r2dsdx(k,6,2) = 0.5_dp*(xvec(k) + (xvec(k) - xcom))
        enddo
      endif
    endif
  elseif (ndim.eq.2) then
    if (lfinitestrain) then
      do k = 1,nvec
!
!  Create vector scaled by inverse strain matrix
!
        xyzloc(1) = straininverse(1,1)*(xvec(k) - xcom) + &
                    straininverse(2,1)*(yvec(k) - ycom)
        xyzloc(2) = straininverse(1,2)*(xvec(k) - xcom) + &
                    straininverse(2,2)*(yvec(k) - ycom)
!
        rpd(1) = xyzloc(1)*xyzloc(1)
        rpd(2) = xyzloc(2)*xyzloc(2)
        rpd(3) = xyzloc(1)*xyzloc(2)
!
        dr2ds(k,1) = (1.0_dp + strain(1))*rpd(1) + &
                      0.5_dp*strain(3)*rpd(3) + xyzloc(1)*xcom
        dr2ds(k,2) = (1.0_dp + strain(2))*rpd(2) + &
                      0.5_dp*strain(3)*rpd(3) + xyzloc(2)*ycom
        dr2ds(k,6) = (1.0_dp + 0.5_dp*(strain(1) + strain(2)))*rpd(3) + &
                      0.25_dp*strain(3)*(rpd(1) + rpd(2)) + &
                      0.5_dp*(xyzloc(1)*ycom + xyzloc(2)*xcom)
        dr2ds(k,3:5) = d2r2dx2(k,3:5)
        if (lgrad2) then
!
!  Initialise second derivatives
!
          d2r2ds2(k,1:6,1:6) = 0.0_dp
          d2r2dsdx(k,1:6,1:3) = 0.0_dp
!
!  Second derivatives with respect to 2 strains
!
          d2r2ds2(k,1,1) = rpd(1)
          d2r2ds2(k,6,1) = 0.5_dp*rpd(3)
          d2r2ds2(k,1,6) = 0.5_dp*rpd(3)
!
          d2r2ds2(k,2,2) = rpd(2)
          d2r2ds2(k,6,2) = 0.5_dp*rpd(3)
          d2r2ds2(k,2,6) = 0.5_dp*rpd(3)
!
          d2r2ds2(k,6,6) = 0.25_dp*(rpd(1) + rpd(2))
!
!  Second derivatives with respect to 1 strain and a Cartesian coordinate for fractionals
!
          d2r2dsdx(k,1,1) = 2.0_dp*xyzloc(1)*(1.0_dp + strain(1)) + &
                            0.5_dp*strain(3)*xyzloc(2) + xcom
          d2r2dsdx(k,1,2) = 0.5_dp*strain(3)*xyzloc(1)
!
          d2r2dsdx(k,2,1) = 0.5_dp*strain(3)*xyzloc(2)
          d2r2dsdx(k,2,2) = 2.0_dp*xyzloc(2)*(1.0_dp + strain(2)) + &
                            0.5_dp*strain(3)*xyzloc(1) + ycom
!
          d2r2dsdx(k,6,1) = (1.0_dp + 0.5_dp*(strain(1) + strain(2)))*xyzloc(2) + &
                             0.5_dp*strain(3)*xyzloc(1) + 0.5_dp*ycom
          d2r2dsdx(k,6,2) = (1.0_dp + 0.5_dp*(strain(1) + strain(2)))*xyzloc(1) + &
                             0.5_dp*strain(3)*xyzloc(2) + 0.5_dp*xcom
!
!  Transform mixed derivatives by inverse strain matrix
!
          do i = 1,6
            tmp(1:2) = 0.0_dp
            do j = 1,2
              do l = 1,2
                tmp(j) = tmp(j) + d2r2dsdx(k,i,l)*straininverse(l,j)
              enddo
            enddo
            d2r2dsdx(k,i,1:2) = tmp(1:2)
          enddo
        endif
      enddo
    else
      dr2ds(1:nvec,1:6) = d2r2dx2(1:nvec,1:6)
      if (lgrad2) then
        do k = 1,nvec
          rpd(1) = xvec(k)*xvec(k)
          rpd(2) = yvec(k)*yvec(k)
          rpd(3) = xvec(k)*yvec(k)
!
!  Initialise second derivatives
!
          d2r2ds2(k,1:6,1:6) = 0.0_dp
          d2r2dsdx(k,1:6,1:3) = 0.0_dp
!
!  Second derivatives with respect to 2 strains
!
          d2r2ds2(k,1,1) = (xvec(k) - xcom)*((xvec(k) - xcom) + xvec(k))
          d2r2ds2(k,6,1) = 0.5_dp*((yvec(k) - ycom)*(xvec(k) - xcom) + xvec(k)*(yvec(k) - ycom))
          d2r2ds2(k,1,6) = 0.5_dp*((yvec(k) - ycom)*(xvec(k) - xcom) + xvec(k)*(yvec(k) - ycom))
!
          d2r2ds2(k,2,2) = (yvec(k) - ycom)*((yvec(k) - ycom) + yvec(k))
          d2r2ds2(k,6,2) = 0.5_dp*((xvec(k) - xcom)*(yvec(k) - ycom) + yvec(k)*(xvec(k) - xcom))
          d2r2ds2(k,2,6) = 0.5_dp*((xvec(k) - xcom)*(yvec(k) - ycom) + yvec(k)*(xvec(k) - xcom))
!
          d2r2ds2(k,6,6) = 0.25_dp*((xvec(k) - xcom)**2 + (yvec(k) - ycom)**2 + &
                                    xvec(k)*(xvec(k) - xcom) + yvec(k)*(yvec(k) - ycom))
!
!  Second derivatives with respect to 1 strain and a Cartesian coordinate for fractionals
!
          d2r2dsdx(k,1,1) = xvec(k) + (xvec(k) - xcom)
          d2r2dsdx(k,2,2) = yvec(k) + (yvec(k) - ycom)
          d2r2dsdx(k,6,1) = 0.5_dp*(yvec(k) + (yvec(k) - ycom))
          d2r2dsdx(k,6,2) = 0.5_dp*(xvec(k) + (xvec(k) - xcom))
        enddo
      endif
    endif
  elseif (ndim.eq.1) then
    if (lfinitestrain) then
      do k = 1,nvec
        xyzloc(1) = (xvec(k) - xcom)*straininverse(1,1)
        rpd(1) = xyzloc(1)*xyzloc(1)
        dr2ds(k,1) = (1.0_dp + strain(1))*rpd(1) + xyzloc(1)*xcom
        if (lgrad2) then
          d2r2ds2(k,1,1) = 0.0_dp
          d2r2dsdx(k,1,1:3) = 0.0_dp
!
!  Second derivatives with respect to 2 strains
!
          d2r2ds2(k,1,1) = rpd(1)
!
!  Second derivatives with respect to 1 strain and a Cartesian coordinate for fractionals
!
          d2r2dsdx(k,1,1) = 2.0_dp*xyzloc(1)*(1.0_dp + strain(1)) + xcom
!
!  Transform mixed derivatives by inverse strain matrix
!
          d2r2dsdx(k,1,1) = d2r2dsdx(k,1,1)*straininverse(1,1)
        endif
      enddo
    else
      dr2ds(1:nvec,1) = d2r2dx2(1:nvec,1)
      if (lgrad2) then
        do k = 1,nvec
          rpd(1) = xvec(k)*xvec(k)
!
!  Second derivatives with respect to 2 strains
!
          d2r2ds2(k,1,1) = (xvec(k) - xcom)*((xvec(k) - xcom) + xvec(k))
!
!  Second derivatives with respect to 1 strain and a Cartesian coordinate for fractionals
!
          d2r2dsdx(k,1,1) = xvec(k) + (xvec(k) - xcom)
          d2r2dsdx(k,1,2) = 0.0_dp
          d2r2dsdx(k,1,3) = 0.0_dp
        enddo
      endif
    endif
  endif
!
!  Reset Cartesian second derivatives without centre of mass subtraction
!
  if (lrigid.and.lgrad2) then
    do k = 1,nvec
      d2r2dx2(k,1) = xvec(k)*xvec(k)
      d2r2dx2(k,2) = yvec(k)*yvec(k)
      d2r2dx2(k,3) = zvec(k)*zvec(k)
      d2r2dx2(k,4) = yvec(k)*zvec(k)
      d2r2dx2(k,5) = xvec(k)*zvec(k)
      d2r2dx2(k,6) = xvec(k)*yvec(k)
    enddo
  endif
#ifdef TRACE
  call trace_out('twostrterms')
#endif
!
  return
  end
!
  subroutine gstrterms(ndim,maxgvec,ngvec,xrk,yrk,zrk,dg2ds,d2g2dx2,d2g2ds2,lgrad2)
!
!  Subroutine for calculating strain terms in reciprocal space
!
!  On entry :
!
!   ndim              = number of dimensions
!   ngvec             = number of reciprocal space vectors to be processed
!   maxgvec           = maximum number of reciprocal space vectors (lower dimension of arrays)
!   xrk, yrk, zrk     = Cartesian components of reciprocal space vectors
!   lgrad2            = if true then compute d2g2ds2 and d2g2dsdx
!
!  On exit :
!
!   d2g2dx2           = second Cartesian derivative components in reciprocal space
!   dg2ds             = 1/2 x first derivative of G^2 with respect to strain (can be the same as d2g2dx2)
!   d2g2ds2           = 1/2 x second derivative of G^2 with respect to two strains
!   d2g2dsdx          = 1/2 x second derivative of G^2 with respect to strains and a Cartesian component
!  
!   9/18 Created from twostrterms
!  11/18 Finite strain flag introduced instead of lstraincell
!   5/19 Second derivatives of g2 changed so that strfin modification is no longer needed
!   8/19 Correction to 1-D case
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
!  Julian Gale, CIC, Curtin University, August 2019
!
  use datatypes
  use derivatives,   only : lfinitestrain
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)  :: ndim
  integer(i4), intent(in)  :: ngvec
  integer(i4), intent(in)  :: maxgvec
  logical,     intent(in)  :: lgrad2
  real(dp),    intent(out) :: d2g2dx2(maxgvec,6)
  real(dp),    intent(out) :: dg2ds(maxgvec,6)
  real(dp),    intent(out) :: d2g2ds2(maxgvec,6,6)
  real(dp),    intent(in)  :: xrk(ngvec)
  real(dp),    intent(in)  :: yrk(ngvec)
  real(dp),    intent(in)  :: zrk(ngvec)
!
!  Local variables
!
  integer(i4)              :: i
  integer(i4)              :: j
  integer(i4)              :: k
  integer(i4)              :: l
  integer(i4)              :: m
  real(dp)                 :: rpd(6)
  real(dp)                 :: gloc(3)
  real(dp)                 :: Imat(3,3)
  real(dp)                 :: dIds(6,3,3)
  real(dp)                 :: d2Ids2(6,6,3,3)
#ifdef TRACE
  call trace_in('gstrterms')
#endif
!
  if (ndim.eq.3) then
!
!  Set up strain products for 3-D case
!
    do k = 1,ngvec
      d2g2dx2(k,1) = xrk(k)*xrk(k)
      d2g2dx2(k,2) = yrk(k)*yrk(k)
      d2g2dx2(k,3) = zrk(k)*zrk(k)
      d2g2dx2(k,4) = yrk(k)*zrk(k)
      d2g2dx2(k,5) = xrk(k)*zrk(k)
      d2g2dx2(k,6) = xrk(k)*yrk(k)
    enddo
    if (lfinitestrain) then
!
!  Compute inverse matrix
!
      straininv(1,1) = strainfull(2,2)*strainfull(3,3) - strainfull(3,2)**2
      straininv(2,1) = strainfull(3,2)*strainfull(3,1) - strainfull(2,1)*strainfull(3,3)
      straininv(3,1) = strainfull(3,2)*strainfull(2,1) - strainfull(3,1)*strainfull(2,2)
      straininv(2,2) = strainfull(1,1)*strainfull(3,3) - strainfull(3,1)**2
      straininv(3,2) = strainfull(3,1)*strainfull(2,1) - strainfull(3,2)*strainfull(1,1)
      straininv(3,3) = strainfull(1,1)*strainfull(2,2) - strainfull(2,1)**2
      straininv(1,2) = straininv(2,1)
      straininv(1,3) = straininv(3,1)
      straininv(2,3) = straininv(3,2)
!
!  Compute unique derivatives of inverse matrix with respect to strains
!  LHS = strain index; RHS = matrix element index
!
      dstraininvds(1:6,1:3,1:3) = 0.0_dp
!
      dstraininvds(2,1,1) = strainfull(3,3)
      dstraininvds(3,1,1) = strainfull(2,2)
      dstraininvds(4,1,1) = - strainfull(3,2)
!
      dstraininvds(1,2,2) = strainfull(3,3)
      dstraininvds(3,2,2) = strainfull(1,1)
      dstraininvds(5,2,2) = - strainfull(3,1)
!
      dstraininvds(1,3,3) = strainfull(2,2)
      dstraininvds(2,3,3) = strainfull(1,1)
      dstraininvds(6,3,3) = - strainfull(2,1)
!
      dstraininvds(1,3,2) = - strainfull(3,2)
      dstraininvds(4,3,2) = - 0.5_dp*strainfull(1,1)
      dstraininvds(5,3,2) = 0.5_dp*strainfull(2,1)
      dstraininvds(6,3,2) = 0.5_dp*strainfull(3,1)
!
      dstraininvds(2,3,1) = - strainfull(3,1)
      dstraininvds(4,3,1) = 0.5_dp*strainfull(2,1)
      dstraininvds(5,3,1) = - 0.5_dp*strainfull(2,2)
      dstraininvds(6,3,1) = 0.5_dp*strainfull(3,2)
!
      dstraininvds(3,2,1) = - strainfull(2,1)
      dstraininvds(4,2,1) = 0.5_dp*strainfull(3,1)
      dstraininvds(5,2,1) = 0.5_dp*strainfull(3,2)
      dstraininvds(6,2,1) = - 0.5_dp*strainfull(3,3)
!
      dstraininvds(3,1,2) = dstraininvds(3,2,1)
      dstraininvds(4,1,2) = dstraininvds(4,2,1)
      dstraininvds(5,1,2) = dstraininvds(5,2,1)
      dstraininvds(6,1,2) = dstraininvds(6,2,1)
!
      dstraininvds(2,1,3) = dstraininvds(2,3,1)
      dstraininvds(4,1,3) = dstraininvds(4,3,1)
      dstraininvds(5,1,3) = dstraininvds(5,3,1)
      dstraininvds(6,1,3) = dstraininvds(6,3,1)
!
      dstraininvds(1,2,3) = dstraininvds(1,3,2)
      dstraininvds(4,2,3) = dstraininvds(4,3,2)
      dstraininvds(5,2,3) = dstraininvds(5,3,2)
      dstraininvds(6,2,3) = dstraininvds(6,3,2)
!
!  Create overall strain derivatives of inverse matrix
!
      do i = 1,3
        do j = 1,3
          Imat(j,i) = straininv(j,i)*straindet
          do l = 1,6
            dIds(l,j,i) = straindet*(dstraininvds(l,j,i) - straindet*strainddetds(l)*straininv(j,i))
          enddo
        enddo
      enddo
      if (lgrad2) then
!
!  Compute second derivatives of inverse strain matrix
!
        d2straininvds2(1:6,1:6,1:3,1:3) = 0.0_dp
!
        d2straininvds2(3,2,1,1) = 1.0_dp
        d2straininvds2(2,3,1,1) = 1.0_dp
        d2straininvds2(4,4,1,1) = - 0.5_dp
!
        d2straininvds2(3,1,2,2) = 1.0_dp
        d2straininvds2(1,3,2,2) = 1.0_dp
        d2straininvds2(5,5,2,2) = - 0.5_dp
!
        d2straininvds2(2,1,3,3) = 1.0_dp
        d2straininvds2(1,2,3,3) = 1.0_dp
        d2straininvds2(6,6,3,3) = - 0.5_dp
!
        d2straininvds2(4,1,3,2) = - 0.5_dp
        d2straininvds2(1,4,3,2) = - 0.5_dp
        d2straininvds2(6,5,3,2) = 0.25_dp
        d2straininvds2(5,6,3,2) = 0.25_dp
!
        d2straininvds2(4,1,2,3) = - 0.5_dp
        d2straininvds2(1,4,2,3) = - 0.5_dp
        d2straininvds2(6,5,2,3) = 0.25_dp
        d2straininvds2(5,6,2,3) = 0.25_dp
!
        d2straininvds2(5,2,3,1) = - 0.5_dp
        d2straininvds2(6,4,3,1) = 0.25_dp
        d2straininvds2(2,5,3,1) = - 0.5_dp
        d2straininvds2(4,6,3,1) = 0.25_dp
!
        d2straininvds2(5,2,1,3) = - 0.5_dp
        d2straininvds2(6,4,1,3) = 0.25_dp
        d2straininvds2(2,5,1,3) = - 0.5_dp
        d2straininvds2(4,6,1,3) = 0.25_dp
!
        d2straininvds2(6,3,2,1) = - 0.5_dp
        d2straininvds2(5,4,2,1) = 0.25_dp
        d2straininvds2(4,5,2,1) = 0.25_dp
        d2straininvds2(3,6,2,1) = - 0.5_dp
!
        d2straininvds2(6,3,1,2) = - 0.5_dp
        d2straininvds2(5,4,1,2) = 0.25_dp
        d2straininvds2(4,5,1,2) = 0.25_dp
        d2straininvds2(3,6,1,2) = - 0.5_dp
!
!  Create overall strain second derivatives of inverse matrix
!
        do i = 1,3
          do j = 1,3
            do l = 1,6
              do m = 1,6
                d2Ids2(m,l,j,i) = straindet*(d2straininvds2(m,l,j,i) &
                                  + 2.0_dp*straindet*straindet*strainddetds(m)*strainddetds(l)*straininv(j,i) &
                                  - straindet*straind2detds2(m,l)*straininv(j,i) &
                                  - straindet*strainddetds(l)*dstraininvds(m,j,i) &
                                  - straindet*strainddetds(m)*dstraininvds(l,j,i))
              enddo
            enddo
          enddo
        enddo
      endif
!
      do k = 1,ngvec
!
!  Create vector scaled by inverse strain matrix
!
        gloc(1) = strainfull(1,1)*xrk(k) + &
                  strainfull(2,1)*yrk(k) + &
                  strainfull(3,1)*zrk(k)
        gloc(2) = strainfull(1,2)*xrk(k) + &
                  strainfull(2,2)*yrk(k) + &
                  strainfull(3,2)*zrk(k)
        gloc(3) = strainfull(1,3)*xrk(k) + &
                  strainfull(2,3)*yrk(k) + &
                  strainfull(3,3)*zrk(k)
!
        rpd(1) = gloc(1)*gloc(1)
        rpd(2) = gloc(2)*gloc(2)
        rpd(3) = gloc(3)*gloc(3)
        rpd(4) = gloc(2)*gloc(3)
        rpd(5) = gloc(1)*gloc(3)
        rpd(6) = gloc(1)*gloc(2)
!
!  Compute strain first derivatives of G^2
!
        dg2ds(k,1:6) = 0.0_dp
        do i = 1,3
          do j = 1,6
            dg2ds(k,j) = dg2ds(k,j) + (dIds(j,1,i)*Imat(i,1)*rpd(1) + &
                                       dIds(j,2,i)*Imat(i,2)*rpd(2) + &
                                       dIds(j,3,i)*Imat(i,3)*rpd(3) + &
                                       dIds(j,3,i)*Imat(i,2)*rpd(4) + &
                                       dIds(j,2,i)*Imat(i,3)*rpd(4) + &
                                       dIds(j,3,i)*Imat(i,1)*rpd(5) + &
                                       dIds(j,1,i)*Imat(i,3)*rpd(5) + &
                                       dIds(j,2,i)*Imat(i,1)*rpd(6) + &
                                       dIds(j,1,i)*Imat(i,2)*rpd(6))
          enddo
        enddo
        if (lgrad2) then
!
!  Compute second derivatives at finite strain - factor of half included by
!  only summing over half of the terms due to symmetry
!
          d2g2ds2(k,1:6,1:6) = 0.0_dp
          do i = 1,3
            do j = 1,6
              do l = 1,6
                d2g2ds2(k,l,j) = d2g2ds2(k,l,j) + (d2Ids2(l,j,1,i)*Imat(i,1)*rpd(1) + &
                                                   d2Ids2(l,j,2,i)*Imat(i,2)*rpd(2) + &
                                                   d2Ids2(l,j,3,i)*Imat(i,3)*rpd(3) + &
                                                   d2Ids2(l,j,3,i)*Imat(i,2)*rpd(4) + &
                                                   d2Ids2(l,j,2,i)*Imat(i,3)*rpd(4) + &
                                                   d2Ids2(l,j,3,i)*Imat(i,1)*rpd(5) + &
                                                   d2Ids2(l,j,1,i)*Imat(i,3)*rpd(5) + &
                                                   d2Ids2(l,j,2,i)*Imat(i,1)*rpd(6) + &
                                                   d2Ids2(l,j,1,i)*Imat(i,2)*rpd(6))
                d2g2ds2(k,l,j) = d2g2ds2(k,l,j) + (dIds(j,1,i)*dIds(l,i,1)*rpd(1) + &
                                                   dIds(j,2,i)*dIds(l,i,2)*rpd(2) + &
                                                   dIds(j,3,i)*dIds(l,i,3)*rpd(3) + &
                                                   dIds(j,3,i)*dIds(l,i,2)*rpd(4) + &
                                                   dIds(j,2,i)*dIds(l,i,3)*rpd(4) + &
                                                   dIds(j,3,i)*dIds(l,i,1)*rpd(5) + &
                                                   dIds(j,1,i)*dIds(l,i,3)*rpd(5) + &
                                                   dIds(j,2,i)*dIds(l,i,1)*rpd(6) + &
                                                   dIds(j,1,i)*dIds(l,i,2)*rpd(6))
              enddo
            enddo
          enddo
        endif
      enddo
    else
      dg2ds(1:ngvec,1:6) = - d2g2dx2(1:ngvec,1:6)
      if (lgrad2) then
        d2g2ds2(1:ngvec,1:6,1:6) = 0.0_dp
!
        do k = 1,ngvec
          rpd(1) = xrk(k)*xrk(k)
          rpd(2) = yrk(k)*yrk(k)
          rpd(3) = zrk(k)*zrk(k)
          rpd(4) = yrk(k)*zrk(k)
          rpd(5) = xrk(k)*zrk(k)
          rpd(6) = xrk(k)*yrk(k)
!
          d2g2ds2(k,1,1) = 2.0_dp*rpd(1)
          d2g2ds2(k,2,2) = 2.0_dp*rpd(2)
          d2g2ds2(k,3,3) = 2.0_dp*rpd(3)
!
          d2g2ds2(k,5,1) = rpd(5)
          d2g2ds2(k,6,1) = rpd(6)
          d2g2ds2(k,4,2) = rpd(4)
          d2g2ds2(k,6,2) = rpd(6)
          d2g2ds2(k,4,3) = rpd(4)
          d2g2ds2(k,5,3) = rpd(5)
!
          d2g2ds2(k,1,5) = rpd(5)
          d2g2ds2(k,1,6) = rpd(6)
          d2g2ds2(k,2,4) = rpd(4)
          d2g2ds2(k,2,6) = rpd(6)
          d2g2ds2(k,3,4) = rpd(4)
          d2g2ds2(k,3,5) = rpd(5)
!
          d2g2ds2(k,4,4) = 0.5_dp*(rpd(2) + rpd(3))
          d2g2ds2(k,5,5) = 0.5_dp*(rpd(1) + rpd(3))
          d2g2ds2(k,6,6) = 0.5_dp*(rpd(1) + rpd(2))
!
          d2g2ds2(k,5,4) = 0.5_dp*rpd(6)
          d2g2ds2(k,6,4) = 0.5_dp*rpd(5)
          d2g2ds2(k,6,5) = 0.5_dp*rpd(4)
          d2g2ds2(k,4,5) = 0.5_dp*rpd(6)
          d2g2ds2(k,4,6) = 0.5_dp*rpd(5)
          d2g2ds2(k,5,6) = 0.5_dp*rpd(4)
!
          d2g2ds2(k,4,5) = 0.5_dp*rpd(6)
          d2g2ds2(k,4,6) = 0.5_dp*rpd(5)
          d2g2ds2(k,5,6) = 0.5_dp*rpd(4)
          d2g2ds2(k,5,4) = 0.5_dp*rpd(6)
          d2g2ds2(k,6,4) = 0.5_dp*rpd(5)
          d2g2ds2(k,6,5) = 0.5_dp*rpd(4)
        enddo
      endif
    endif
  elseif (ndim.eq.2) then
!
!  Set up strain products for 2-D case
!
    do k = 1,ngvec
      d2g2dx2(k,1) = xrk(k)*xrk(k)
      d2g2dx2(k,2) = yrk(k)*yrk(k)
      d2g2dx2(k,3) = xrk(k)*yrk(k)
    enddo
    if (lfinitestrain) then
!
!  Compute inverse matrix
!
      straininv(1,1) = strainfull(2,2)
      straininv(2,1) = - strainfull(2,1)
      straininv(1,2) = - strainfull(2,1)
      straininv(2,2) = strainfull(1,1)
!
!  Compute unique derivatives of inverse matrix with respect to strains
!  LHS = strain index; RHS = matrix element index
!
      dstraininvds(1:3,1:2,1:2) = 0.0_dp
!
      dstraininvds(2,1,1) = 1.0_dp
      dstraininvds(1,2,2) = 1.0_dp
      dstraininvds(3,2,1) = - 0.5_dp
      dstraininvds(3,1,2) = - 0.5_dp
!
!  Create overall strain derivatives of inverse matrix
!
      do i = 1,2
        do j = 1,2
          Imat(j,i) = straininv(j,i)*straindet
          do l = 1,3
            dIds(l,j,i) = straindet*(dstraininvds(l,j,i) - straindet*strainddetds(l)*straininv(j,i))
          enddo
        enddo
      enddo
      if (lgrad2) then
!
!  Create overall strain second derivatives of inverse matrix
!
        do i = 1,2
          do j = 1,2
            do l = 1,3
              do m = 1,3
                d2Ids2(m,l,j,i) = straindet*( &
                                  + 2.0_dp*straindet*straindet*strainddetds(m)*strainddetds(l)*straininv(j,i) &
                                  - straindet*straind2detds2(m,l)*straininv(j,i) &
                                  - straindet*strainddetds(l)*dstraininvds(m,j,i) &
                                  - straindet*strainddetds(m)*dstraininvds(l,j,i))
              enddo
            enddo
          enddo
        enddo
      endif
!
      do k = 1,ngvec
!
!  Create vector scaled by inverse strain matrix
!
        gloc(1) = strainfull(1,1)*xrk(k) + &
                  strainfull(2,1)*yrk(k)
        gloc(2) = strainfull(1,2)*xrk(k) + &
                  strainfull(2,2)*yrk(k)
!
        rpd(1) = gloc(1)*gloc(1)
        rpd(2) = gloc(2)*gloc(2)
        rpd(3) = gloc(1)*gloc(2)
!
!  Compute strain first derivatives of G^2
!
        dg2ds(k,1:3) = 0.0_dp
        do i = 1,2
          do j = 1,3
            dg2ds(k,j) = dg2ds(k,j) + (dIds(j,1,i)*Imat(i,1)*rpd(1) + &
                                       dIds(j,2,i)*Imat(i,2)*rpd(2) + &
                                       dIds(j,2,i)*Imat(i,1)*rpd(3) + &
                                       dIds(j,1,i)*Imat(i,2)*rpd(3))
          enddo
        enddo
        if (lgrad2) then
!
!  Compute second derivatives at finite strain - factor of half included by
!  only summing over half of the terms due to symmetry
!
          d2g2ds2(k,1:3,1:3) = 0.0_dp
          do i = 1,2
            do j = 1,3
              do l = 1,3
                d2g2ds2(k,l,j) = d2g2ds2(k,l,j) + (d2Ids2(l,j,1,i)*Imat(i,1)*rpd(1) + &
                                                   d2Ids2(l,j,2,i)*Imat(i,2)*rpd(2) + &
                                                   d2Ids2(l,j,2,i)*Imat(i,1)*rpd(3) + &
                                                   d2Ids2(l,j,1,i)*Imat(i,2)*rpd(3))
                d2g2ds2(k,l,j) = d2g2ds2(k,l,j) + (dIds(j,1,i)*dIds(l,i,1)*rpd(1) + &
                                                   dIds(j,2,i)*dIds(l,i,2)*rpd(2) + &
                                                   dIds(j,2,i)*dIds(l,i,1)*rpd(3) + &
                                                   dIds(j,1,i)*dIds(l,i,2)*rpd(3))
              enddo
            enddo
          enddo
        endif
      enddo
    else
      dg2ds(1:ngvec,1:3) = - d2g2dx2(1:ngvec,1:3)
      if (lgrad2) then
        d2g2ds2(1:ngvec,1:3,1:3) = 0.0_dp
!
        do k = 1,ngvec
          rpd(1) = xrk(k)*xrk(k)
          rpd(2) = yrk(k)*yrk(k)
          rpd(3) = xrk(k)*yrk(k)
!
          d2g2ds2(k,1,1) = 2.0_dp*rpd(1)
          d2g2ds2(k,2,2) = 2.0_dp*rpd(2)
          d2g2ds2(k,3,3) = 0.5_dp*(rpd(1) + rpd(2))
!
          d2g2ds2(k,3,1) = rpd(3)
          d2g2ds2(k,3,2) = rpd(3)
          d2g2ds2(k,1,3) = rpd(3)
          d2g2ds2(k,2,3) = rpd(3)
        enddo
      endif
    endif
  elseif (ndim.eq.1) then
!
!  Set up strain products for 1-D case
!
    do k = 1,ngvec
      d2g2dx2(k,1) = xrk(k)*xrk(k)
    enddo
    if (lfinitestrain) then
      do k = 1,ngvec
        gloc(1) = xrk(k)*strainfull(1,1)
        rpd(1) = gloc(1)*gloc(1)
        dg2ds(k,1) = - rpd(1)*straindet**3
        if (lgrad2) then
!
!  Second derivatives
!
          d2g2ds2(k,1,1) = 3.0_dp*rpd(1)*straindet**4
        endif
      enddo
    else
      dg2ds(1:ngvec,1) = - d2g2dx2(1:ngvec,1)
      if (lgrad2) then
        do k = 1,ngvec
          d2g2ds2(k,1,1) = 2.0_dp*xrk(k)*xrk(k)
        enddo
      endif
    endif
  endif
#ifdef TRACE
  call trace_out('gstrterms')
#endif
!
  return
  end

  subroutine realstrterms(ndim,maxvec,nvec,xvec,yvec,zvec,xcom,ycom,zcom,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,lgrad2)
!
!  Subroutine for calculating real space two-body strain terms
!
!  NB: Order of array dimensions is swapped relative to twostrterms
!
!  On entry :
!
!   ndim              = number of dimensions
!   nvec              = number of vectors to be processed
!   maxvec            = maximum number of vectors (lower dimension of d2r2dx2, etc)
!   xvec, yvec, zvec  = Cartesian components of interatomic vectors
!   xcom, ycom, zcom  = Cartesian components of intramolecular vectors to the centre of mass
!   lgrad2            = if true then compute d2r2ds2 and d2r2dsdx
!
!  On exit :
!
!   d2r2dx2           = second Cartesian derivative component
!   dr2ds             = 1/2 x first derivative of r^2 with respect to strain (can be the same as d2r2dx2)
!   d2r2ds2           = 1/2 x second derivative of r^2 with respect to two strains
!   d2r2dsdx          = 1/2 x second derivative of r^2 with respect to strains and a Cartesian component
!  
!   9/18 Created from twostrterms
!  11/18 Finite strain flag introduced instead of lstraincell
!   5/19 Second derivatives with respect to strain now complete - no need for strfin
!   4/20 Centre of mass added for rigid molecules
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
!  Julian Gale, CIC, Curtin University, April 2020
!
  use datatypes
  use control,       only : lrigid
  use current,       only : strain
  use derivatives,   only : lfinitestrain
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)  :: ndim
  integer(i4), intent(in)  :: nvec
  integer(i4), intent(in)  :: maxvec
  logical,     intent(in)  :: lgrad2
  real(dp),    intent(out) :: d2r2dx2(6,maxvec)
  real(dp),    intent(out) :: dr2ds(6,maxvec)
  real(dp),    intent(out) :: d2r2ds2(6,6,maxvec)
  real(dp),    intent(out) :: d2r2dsdx(6,3,maxvec)
  real(dp),    intent(in)  :: xvec(nvec)
  real(dp),    intent(in)  :: yvec(nvec)
  real(dp),    intent(in)  :: zvec(nvec)
  real(dp),    intent(in)  :: xcom(nvec)
  real(dp),    intent(in)  :: ycom(nvec)
  real(dp),    intent(in)  :: zcom(nvec)
!
!  Local variables
!
  integer(i4)              :: i
  integer(i4)              :: j
  integer(i4)              :: k
  integer(i4)              :: l
  real(dp)                 :: rpd(6)
  real(dp)                 :: tmp(3)
  real(dp)                 :: xyzloc(3)
#ifdef TRACE
  call trace_in('realstrterms')
#endif
!
!  Set up Cartesian second derivatives - this is the same for all dimensionalities and needed for strains
!
  if (lrigid) then
    do k = 1,nvec
      d2r2dx2(1,k) = xvec(k)*(xvec(k) - xcom(k))
      d2r2dx2(2,k) = yvec(k)*(yvec(k) - ycom(k))
      d2r2dx2(3,k) = zvec(k)*(zvec(k) - zcom(k))
      d2r2dx2(4,k) = 0.5_dp*(yvec(k)*(zvec(k) - zcom(k)) + zvec(k)*(yvec(k) - ycom(k)))
      d2r2dx2(5,k) = 0.5_dp*(xvec(k)*(zvec(k) - zcom(k)) + zvec(k)*(xvec(k) - xcom(k)))
      d2r2dx2(6,k) = 0.5_dp*(xvec(k)*(yvec(k) - ycom(k)) + yvec(k)*(xvec(k) - xcom(k)))
    enddo
  else
    do k = 1,nvec
      d2r2dx2(1,k) = xvec(k)*xvec(k)
      d2r2dx2(2,k) = yvec(k)*yvec(k)
      d2r2dx2(3,k) = zvec(k)*zvec(k)
      d2r2dx2(4,k) = yvec(k)*zvec(k)
      d2r2dx2(5,k) = xvec(k)*zvec(k)
      d2r2dx2(6,k) = xvec(k)*yvec(k)
    enddo
  endif
!
  if (ndim.eq.3) then
    if (lfinitestrain) then
      do k = 1,nvec
!
!  Create vector scaled by inverse strain matrix
!
        xyzloc(1) = straininverse(1,1)*(xvec(k) - xcom(k)) + &
                    straininverse(2,1)*(yvec(k) - ycom(k)) + &
                    straininverse(3,1)*(zvec(k) - zcom(k))
        xyzloc(2) = straininverse(1,2)*(xvec(k) - xcom(k)) + &
                    straininverse(2,2)*(yvec(k) - ycom(k)) + &
                    straininverse(3,2)*(zvec(k) - zcom(k))
        xyzloc(3) = straininverse(1,3)*(xvec(k) - xcom(k)) + &
                    straininverse(2,3)*(yvec(k) - ycom(k)) + &
                    straininverse(3,3)*(zvec(k) - zcom(k))
!
        rpd(1) = xyzloc(1)*xyzloc(1)
        rpd(2) = xyzloc(2)*xyzloc(2)
        rpd(3) = xyzloc(3)*xyzloc(3)
        rpd(4) = xyzloc(2)*xyzloc(3)
        rpd(5) = xyzloc(1)*xyzloc(3)
        rpd(6) = xyzloc(1)*xyzloc(2)
!
        dr2ds(1,k) = (1.0_dp + strain(1))*rpd(1) + &
                      0.5_dp*(strain(6)*rpd(6) + strain(5)*rpd(5)) + xyzloc(1)*xcom(k)
        dr2ds(2,k) = (1.0_dp + strain(2))*rpd(2) + &
                      0.5_dp*(strain(6)*rpd(6) + strain(4)*rpd(4)) + xyzloc(2)*ycom(k)
        dr2ds(3,k) = (1.0_dp + strain(3))*rpd(3) + &
                      0.5_dp*(strain(5)*rpd(5) + strain(4)*rpd(4)) + xyzloc(3)*zcom(k)
        dr2ds(4,k) = (1.0_dp + 0.5_dp*(strain(2) + strain(3)))*rpd(4) + &
                      0.25_dp*(strain(4)*(rpd(2) + rpd(3)) + &
                               strain(5)*rpd(6) + strain(6)*rpd(5)) + &
                      0.5_dp*(xyzloc(2)*zcom(k) + xyzloc(3)*ycom(k))
        dr2ds(5,k) = (1.0_dp + 0.5_dp*(strain(1) + strain(3)))*rpd(5) + &
                      0.25_dp*(strain(5)*(rpd(1) + rpd(3)) + &
                               strain(4)*rpd(6) + strain(6)*rpd(4)) + &
                      0.5_dp*(xyzloc(1)*zcom(k) + xyzloc(3)*xcom(k))
        dr2ds(6,k) = (1.0_dp + 0.5_dp*(strain(1) + strain(2)))*rpd(6) + &
                      0.25_dp*(strain(6)*(rpd(1) + rpd(2)) + &
                               strain(4)*rpd(5) + strain(5)*rpd(4)) + &
                      0.5_dp*(xyzloc(1)*ycom(k) + xyzloc(2)*xcom(k))
        if (lgrad2) then
!
!  Initialise second derivatives
!
          d2r2ds2(1:6,1:6,k) = 0.0_dp
!
!  Second derivatives with respect to 2 strains
!
          d2r2ds2(1,1,k) = rpd(1)
          d2r2ds2(5,1,k) = 0.5_dp*rpd(5)
          d2r2ds2(6,1,k) = 0.5_dp*rpd(6)
!
          d2r2ds2(2,2,k) = rpd(2)
          d2r2ds2(4,2,k) = 0.5_dp*rpd(4)
          d2r2ds2(6,2,k) = 0.5_dp*rpd(6)
!
          d2r2ds2(3,3,k) = rpd(3)
          d2r2ds2(4,3,k) = 0.5_dp*rpd(4)
          d2r2ds2(5,3,k) = 0.5_dp*rpd(5)
!
          d2r2ds2(4,4,k) = 0.25_dp*(rpd(2) + rpd(3))
          d2r2ds2(5,4,k) = 0.25_dp*rpd(6)
          d2r2ds2(6,4,k) = 0.25_dp*rpd(5)
!
          d2r2ds2(5,5,k) = 0.25_dp*(rpd(1) + rpd(3))
          d2r2ds2(6,5,k) = 0.25_dp*rpd(4)
!
          d2r2ds2(6,6,k) = 0.25_dp*(rpd(1) + rpd(2))
!
          do i = 2,6
            do j = 1,i-1
              d2r2ds2(j,i,k) = d2r2ds2(i,j,k)
            enddo
          enddo
!
!  Second derivatives with respect to 1 strain and a Cartesian coordinate
!
          d2r2dsdx(1,1,k) = 2.0_dp*xyzloc(1)*(1.0_dp + strain(1)) + &
                            0.5_dp*(strain(6)*xyzloc(2) + strain(5)*xyzloc(3)) + xcom(k)
          d2r2dsdx(1,2,k) = 0.5_dp*strain(6)*xyzloc(1)
          d2r2dsdx(1,3,k) = 0.5_dp*strain(5)*xyzloc(1)
!
          d2r2dsdx(2,1,k) = 0.5_dp*strain(6)*xyzloc(2)
          d2r2dsdx(2,2,k) = 2.0_dp*xyzloc(2)*(1.0_dp + strain(2)) + &
                            0.5_dp*(strain(6)*xyzloc(1) + strain(4)*xyzloc(3)) + ycom(k)
          d2r2dsdx(2,3,k) = 0.5_dp*strain(4)*xyzloc(2)
!
          d2r2dsdx(3,1,k) = 0.5_dp*strain(5)*xyzloc(3)
          d2r2dsdx(3,2,k) = 0.5_dp*strain(4)*xyzloc(3)
          d2r2dsdx(3,3,k) = 2.0_dp*xyzloc(3)*(1.0_dp + strain(3)) + &
                            0.5_dp*(strain(5)*xyzloc(1) + strain(4)*xyzloc(2)) + zcom(k)
!
          d2r2dsdx(4,1,k) = 0.25_dp*(strain(5)*xyzloc(2) + strain(6)*xyzloc(3))
          d2r2dsdx(4,2,k) = (1.0_dp + 0.5_dp*(strain(2) + strain(3)))*xyzloc(3) + &
                             0.25_dp*(2.0_dp*strain(4)*xyzloc(2) + strain(5)*xyzloc(1)) + 0.5_dp*zcom(k)
          d2r2dsdx(4,3,k) = (1.0_dp + 0.5_dp*(strain(2) + strain(3)))*xyzloc(2) + &
                             0.25_dp*(2.0_dp*strain(4)*xyzloc(3) + strain(6)*xyzloc(1)) + 0.5_dp*ycom(k)
!
          d2r2dsdx(5,1,k) = (1.0_dp + 0.5_dp*(strain(1) + strain(3)))*xyzloc(3) + &
                             0.25_dp*(2.0_dp*strain(5)*xyzloc(1) + strain(4)*xyzloc(2)) + 0.5_dp*zcom(k)
          d2r2dsdx(5,2,k) = 0.25_dp*(strain(4)*xyzloc(1) + strain(6)*xyzloc(3))
          d2r2dsdx(5,3,k) = (1.0_dp + 0.5_dp*(strain(1) + strain(3)))*xyzloc(1) + &
                             0.25_dp*(2.0_dp*strain(5)*xyzloc(3) + strain(6)*xyzloc(2)) + 0.5_dp*xcom(k)
!
          d2r2dsdx(6,1,k) = (1.0_dp + 0.5_dp*(strain(1) + strain(2)))*xyzloc(2) + &
                             0.25_dp*(2.0_dp*strain(6)*xyzloc(1) + strain(4)*xyzloc(3)) + 0.5_dp*ycom(k)
          d2r2dsdx(6,2,k) = (1.0_dp + 0.5_dp*(strain(1) + strain(2)))*xyzloc(1) + &
                             0.25_dp*(2.0_dp*strain(6)*xyzloc(2) + strain(5)*xyzloc(3)) + 0.5_dp*xcom(k)
          d2r2dsdx(6,3,k) = 0.25_dp*(strain(4)*xyzloc(1) + strain(5)*xyzloc(2))
!
!  Transform mixed derivatives by inverse strain matrix
!
          do i = 1,6
            tmp(1:3) = 0.0_dp
            do j = 1,3
              do l = 1,3
                tmp(j) = tmp(j) + d2r2dsdx(i,l,k)*straininverse(l,j)
              enddo
            enddo
            d2r2dsdx(i,1:3,k) = tmp(1:3)
          enddo
        endif
      enddo
    else
      dr2ds(1:6,1:nvec) = d2r2dx2(1:6,1:nvec)
      if (lgrad2) then
        do k = 1,nvec
          rpd(1) = xvec(k)*xvec(k)
          rpd(2) = yvec(k)*yvec(k)
          rpd(3) = zvec(k)*zvec(k)
          rpd(4) = yvec(k)*zvec(k)
          rpd(5) = xvec(k)*zvec(k)
          rpd(6) = xvec(k)*yvec(k)
!
!  Initialise second derivatives
!
          d2r2ds2(1:6,1:6,k) = 0.0_dp
          d2r2dsdx(1:6,1:3,k) = 0.0_dp
!
!  Second derivatives with respect to 2 strains
!
          d2r2ds2(1,1,k) = (xvec(k) - xcom(k))*((xvec(k) - xcom(k)) + xvec(k))
          d2r2ds2(5,1,k) = 0.5_dp*((zvec(k) - zcom(k))*(xvec(k) - xcom(k)) + xvec(k)*(zvec(k) - zcom(k)))
          d2r2ds2(6,1,k) = 0.5_dp*((yvec(k) - ycom(k))*(xvec(k) - xcom(k)) + xvec(k)*(yvec(k) - ycom(k)))
!
          d2r2ds2(2,2,k) = (yvec(k) - ycom(k))*((yvec(k) - ycom(k)) + yvec(k))
          d2r2ds2(4,2,k) = 0.5_dp*((zvec(k) - zcom(k))*(yvec(k) - ycom(k)) + yvec(k)*(zvec(k) - zcom(k)))
          d2r2ds2(6,2,k) = 0.5_dp*((xvec(k) - xcom(k))*(yvec(k) - ycom(k)) + yvec(k)*(xvec(k) - xcom(k)))
!
          d2r2ds2(3,3,k) = (zvec(k) - zcom(k))*((zvec(k) - zcom(k)) + zvec(k))
          d2r2ds2(4,3,k) = 0.5_dp*((yvec(k) - ycom(k))*(zvec(k) - zcom(k)) + zvec(k)*(yvec(k) - ycom(k)))
          d2r2ds2(5,3,k) = 0.5_dp*((xvec(k) - xcom(k))*(zvec(k) - zcom(k)) + zvec(k)*(xvec(k) - xcom(k)))
!
          d2r2ds2(4,4,k) = 0.25_dp*((yvec(k) - ycom(k))**2 + (zvec(k) - zcom(k))**2 + &
                                    yvec(k)*(yvec(k) - ycom(k)) + zvec(k)*(zvec(k) - zcom(k)))
          d2r2ds2(5,4,k) = 0.25_dp*((xvec(k) - xcom(k))*(yvec(k) - ycom(k)) + yvec(k)*(xvec(k) - xcom(k)))
          d2r2ds2(6,4,k) = 0.25_dp*((xvec(k) - xcom(k))*(zvec(k) - zcom(k)) + zvec(k)*(xvec(k) - xcom(k)))
!
          d2r2ds2(5,5,k) = 0.25_dp*((xvec(k) - xcom(k))**2 + (zvec(k) - zcom(k))**2 + &
                                    xvec(k)*(xvec(k) - xcom(k)) + zvec(k)*(zvec(k) - zcom(k)))
          d2r2ds2(6,5,k) = 0.25_dp*((yvec(k) - ycom(k))*(zvec(k) - zcom(k)) + zvec(k)*(yvec(k) - ycom(k)))
!
          d2r2ds2(6,6,k) = 0.25_dp*((xvec(k) - xcom(k))**2 + (yvec(k) - ycom(k))**2 + &
                                    xvec(k)*(xvec(k) - xcom(k)) + yvec(k)*(yvec(k) - ycom(k)))
!
          do i = 2,6
            do j = 1,i-1
              d2r2ds2(j,i,k) = d2r2ds2(i,j,k)
            enddo
          enddo
!
!  Second derivatives with respect to 1 strain and a Cartesian coordinate
!
          d2r2dsdx(1,1,k) = xvec(k) + (xvec(k) - xcom(k))
          d2r2dsdx(2,2,k) = yvec(k) + (yvec(k) - ycom(k))
          d2r2dsdx(3,3,k) = zvec(k) + (zvec(k) - zcom(k))
!
          d2r2dsdx(4,2,k) = 0.5_dp*(zvec(k) + (zvec(k) - zcom(k)))
          d2r2dsdx(4,3,k) = 0.5_dp*(yvec(k) + (yvec(k) - ycom(k)))
!
          d2r2dsdx(5,1,k) = 0.5_dp*(zvec(k) + (zvec(k) - zcom(k)))
          d2r2dsdx(5,3,k) = 0.5_dp*(xvec(k) + (xvec(k) - xcom(k)))
!
          d2r2dsdx(6,1,k) = 0.5_dp*(yvec(k) + (yvec(k) - ycom(k)))
          d2r2dsdx(6,2,k) = 0.5_dp*(xvec(k) + (xvec(k) - xcom(k)))
        enddo
      endif
    endif
  elseif (ndim.eq.2) then
    if (lfinitestrain) then
      do k = 1,nvec
!
!  Create vector scaled by inverse strain matrix
!
        xyzloc(1) = straininverse(1,1)*(xvec(k) - xcom(k)) + &
                    straininverse(2,1)*(yvec(k) - ycom(k))
        xyzloc(2) = straininverse(1,2)*(xvec(k) - xcom(k)) + &
                    straininverse(2,2)*(yvec(k) - ycom(k))
!
        rpd(1) = xyzloc(1)*xyzloc(1)
        rpd(2) = xyzloc(2)*xyzloc(2)
        rpd(3) = xyzloc(1)*xyzloc(2)
!
        dr2ds(1,k) = (1.0_dp + strain(1))*rpd(1) + &
                      0.5_dp*strain(3)*rpd(3) + xyzloc(1)*xcom(k)
        dr2ds(2,k) = (1.0_dp + strain(2))*rpd(2) + &
                      0.5_dp*strain(3)*rpd(3) + xyzloc(2)*ycom(k)
        dr2ds(6,k) = (1.0_dp + 0.5_dp*(strain(1) + strain(2)))*rpd(3) + &
                      0.25_dp*strain(3)*(rpd(1) + rpd(2)) + &
                      0.5_dp*(xyzloc(1)*ycom(k) + xyzloc(2)*xcom(k))
        dr2ds(3:5,k) = d2r2dx2(3:5,k)
        if (lgrad2) then
!
!  Initialise second derivatives
!
          d2r2ds2(1:6,1:6,k) = 0.0_dp
          d2r2dsdx(1:6,1:3,k) = 0.0_dp
!
!  Second derivatives with respect to 2 strains
!
          d2r2ds2(1,1,k) = rpd(1)
          d2r2ds2(6,1,k) = 0.5_dp*rpd(3)
          d2r2ds2(1,6,k) = 0.5_dp*rpd(3)
!
          d2r2ds2(2,2,k) = rpd(2)
          d2r2ds2(6,2,k) = 0.5_dp*rpd(3)
          d2r2ds2(2,6,k) = 0.5_dp*rpd(3)
!
          d2r2ds2(6,6,k) = 0.25_dp*(rpd(1) + rpd(2))
!
!  Second derivatives with respect to 1 strain and a Cartesian coordinate
!
          d2r2dsdx(1,1,k) = 2.0_dp*xyzloc(1)*(1.0_dp + strain(1)) + &
                            0.5_dp*strain(3)*xyzloc(2) + xcom(k)
          d2r2dsdx(1,2,k) = 0.5_dp*strain(3)*xyzloc(1)
!
          d2r2dsdx(2,1,k) = 0.5_dp*strain(3)*xyzloc(2)
          d2r2dsdx(2,2,k) = 2.0_dp*xyzloc(2)*(1.0_dp + strain(2)) + &
                            0.5_dp*strain(3)*xyzloc(1) + ycom(k)
!
          d2r2dsdx(6,1,k) = (1.0_dp + 0.5_dp*(strain(1) + strain(2)))*xyzloc(2) + &
                             0.5_dp*strain(3)*xyzloc(1) + 0.5_dp*ycom(k)
          d2r2dsdx(6,2,k) = (1.0_dp + 0.5_dp*(strain(1) + strain(2)))*xyzloc(1) + &
                             0.5_dp*strain(3)*xyzloc(2) + 0.5_dp*xcom(k)
!
!  Transform mixed derivatives by inverse strain matrix
!
          do i = 1,6
            tmp(1:2) = 0.0_dp
            do j = 1,2
              do l = 1,2
                tmp(j) = tmp(j) + d2r2dsdx(i,l,k)*straininverse(l,j)
              enddo
            enddo
            d2r2dsdx(i,1:2,k) = tmp(1:2)
          enddo
        endif
      enddo
    else
      dr2ds(1:6,1:nvec) = d2r2dx2(1:6,1:nvec)
      if (lgrad2) then
        do k = 1,nvec
          rpd(1) = xvec(k)*xvec(k)
          rpd(2) = yvec(k)*yvec(k)
          rpd(3) = xvec(k)*yvec(k)
!
!  Initialise second derivatives
!
          d2r2ds2(1:6,1:6,k) = 0.0_dp
          d2r2dsdx(1:6,1:3,k) = 0.0_dp
!
!  Second derivatives with respect to 2 strains
!
          d2r2ds2(1,1,k) = (xvec(k) - xcom(k))*((xvec(k) - xcom(k)) + xvec(k))
          d2r2ds2(6,1,k) = 0.5_dp*((yvec(k) - ycom(k))*(xvec(k) - xcom(k)) + xvec(k)*(yvec(k) - ycom(k)))
          d2r2ds2(1,6,k) = 0.5_dp*((yvec(k) - ycom(k))*(xvec(k) - xcom(k)) + xvec(k)*(yvec(k) - ycom(k)))
!
          d2r2ds2(2,2,k) = (yvec(k) - ycom(k))*((yvec(k) - ycom(k)) + yvec(k))
          d2r2ds2(6,2,k) = 0.5_dp*((xvec(k) - xcom(k))*(yvec(k) - ycom(k)) + yvec(k)*(xvec(k) - xcom(k)))
          d2r2ds2(2,6,k) = 0.5_dp*((xvec(k) - xcom(k))*(yvec(k) - ycom(k)) + yvec(k)*(xvec(k) - xcom(k)))
!
          d2r2ds2(6,6,k) = 0.25_dp*((xvec(k) - xcom(k))**2 + (yvec(k) - ycom(k))**2 + &
                                    xvec(k)*(xvec(k) - xcom(k)) + yvec(k)*(yvec(k) - ycom(k)))
!
!  Second derivatives with respect to 1 strain and a Cartesian coordinate
!
          d2r2dsdx(1,1,k) = xvec(k) + (xvec(k) - xcom(k))
          d2r2dsdx(2,2,k) = yvec(k) + (yvec(k) - ycom(k))
          d2r2dsdx(6,1,k) = 0.5_dp*(yvec(k) + (yvec(k) - ycom(k)))
          d2r2dsdx(6,2,k) = 0.5_dp*(xvec(k) + (xvec(k) - xcom(k)))
        enddo
      endif
    endif
  elseif (ndim.eq.1) then
    if (lfinitestrain) then
      do k = 1,nvec
        xyzloc(1) = (xvec(k) - xcom(k))*straininverse(1,1)
        rpd(1) = xyzloc(1)*xyzloc(1)
        dr2ds(1,k) = (1.0_dp + strain(1))*rpd(1) + xyzloc(1)*xcom(k)
        if (lgrad2) then
          d2r2ds2(1,1,k) = 0.0_dp
          d2r2dsdx(1,1:3,k) = 0.0_dp
!
!  Second derivatives with respect to 2 strains
!
          d2r2ds2(1,1,k) = rpd(1)
!
!  Second derivatives with respect to 1 strain and a Cartesian coordinate
!
          d2r2dsdx(1,1,k) = 2.0_dp*xyzloc(1)*(1.0_dp + strain(1)) + xcom(k)
!
!  Transform mixed derivatives by inverse strain matrix
!
          d2r2dsdx(1,1,k) = d2r2dsdx(1,1,k)*straininverse(1,1)
        endif
      enddo
    else
      dr2ds(1,1:nvec) = d2r2dx2(1,1:nvec)
      if (lgrad2) then
        do k = 1,nvec
          rpd(1) = xvec(k)*xvec(k)
!
!  Second derivatives with respect to 2 strains
!
          d2r2ds2(1,1,k) = (xvec(k) - xcom(k))*((xvec(k) - xcom(k)) + xvec(k))
!
!  Second derivatives with respect to 1 strain and a Cartesian coordinate
!
          d2r2dsdx(1,1,k) = xvec(k) + (xvec(k) - xcom(k))
          d2r2dsdx(1,2,k) = 0.0_dp
          d2r2dsdx(1,3,k) = 0.0_dp
        enddo
      endif
    endif
  endif
#ifdef TRACE
  call trace_out('realstrterms')
#endif
!
  return
  end

  subroutine real1strterm(ndim,xvec,yvec,zvec,xcom,ycom,zcom,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,lgrad2)
!
!  Subroutine for calculating real space two-body strain terms
!  for a single distance.
!
!  On entry :
!
!   ndim              = number of dimensions
!   xvec, yvec, zvec  = Cartesian components of interatomic vectors
!   xcom, ycom, zcom  = Cartesian components of intramolecular vectors to the centre of molecule
!   lgrad2            = if true then compute d2r2ds2 and d2r2dsdx
!
!  On exit :
!
!   d2r2dx2           = second Cartesian derivative component
!   dr2ds             = 1/2 x first derivative of r^2 with respect to strain (can be the same as d2r2dx2)
!   d2r2ds2           = 1/2 x second derivative of r^2 with respect to two strains
!   d2r2dsdx          = 1/2 x second derivative of r^2 with respect to strains and a Cartesian component
!  
!   9/18 Created from realstrterms
!  11/18 Finite strain flag introduced instead of lstraincell
!   5/19 Second derivatives with respect to strain now complete - no need for strfin
!   4/20 Rigid molecule derivatives added
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
!  Julian Gale, CIC, Curtin University, April 2020
!
  use datatypes
  use control,       only : lrigid
  use current,       only : strain
  use derivatives,   only : lfinitestrain
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)  :: ndim
  logical,     intent(in)  :: lgrad2
  real(dp),    intent(out) :: d2r2dx2(6)
  real(dp),    intent(out) :: dr2ds(6)
  real(dp),    intent(out) :: d2r2ds2(6,6)
  real(dp),    intent(out) :: d2r2dsdx(6,3)
  real(dp),    intent(in)  :: xvec
  real(dp),    intent(in)  :: yvec
  real(dp),    intent(in)  :: zvec
  real(dp),    intent(in)  :: xcom
  real(dp),    intent(in)  :: ycom
  real(dp),    intent(in)  :: zcom
!
!  Local variables
!
  integer(i4)              :: i
  integer(i4)              :: j
  integer(i4)              :: l
  real(dp)                 :: rpd(6)
  real(dp)                 :: tmp(3)
  real(dp)                 :: xyzloc(3)
#ifdef TRACE
  call trace_in('real1strterm')
#endif
!
!  Set up Cartesian second derivatives - this is the same for all dimensionalities & also needed for strains
!
  if (lrigid) then
    d2r2dx2(1) = xvec*(xvec - xcom)
    d2r2dx2(2) = yvec*(yvec - ycom)
    d2r2dx2(3) = zvec*(zvec - zcom)
    d2r2dx2(4) = 0.5_dp*(yvec*(zvec - zcom) + zvec*(yvec - ycom))
    d2r2dx2(5) = 0.5_dp*(xvec*(zvec - zcom) + zvec*(xvec - xcom))
    d2r2dx2(6) = 0.5_dp*(xvec*(yvec - ycom) + yvec*(xvec - xcom))
  else
    d2r2dx2(1) = xvec*xvec
    d2r2dx2(2) = yvec*yvec
    d2r2dx2(3) = zvec*zvec
    d2r2dx2(4) = yvec*zvec
    d2r2dx2(5) = xvec*zvec
    d2r2dx2(6) = xvec*yvec
  endif
!
  if (ndim.eq.3) then
    if (lfinitestrain) then
!
!  Create vector scaled by inverse strain matrix
!
      xyzloc(1) = straininverse(1,1)*(xvec - xcom) + &
                  straininverse(2,1)*(yvec - ycom) + &
                  straininverse(3,1)*(zvec - zcom)
      xyzloc(2) = straininverse(1,2)*(xvec - xcom) + &
                  straininverse(2,2)*(yvec - ycom) + &
                  straininverse(3,2)*(zvec - zcom)
      xyzloc(3) = straininverse(1,3)*(xvec - xcom) + &
                  straininverse(2,3)*(yvec - ycom) + &
                  straininverse(3,3)*(zvec - zcom)
!
      rpd(1) = xyzloc(1)*xyzloc(1)
      rpd(2) = xyzloc(2)*xyzloc(2)
      rpd(3) = xyzloc(3)*xyzloc(3)
      rpd(4) = xyzloc(2)*xyzloc(3)
      rpd(5) = xyzloc(1)*xyzloc(3)
      rpd(6) = xyzloc(1)*xyzloc(2)
!
      dr2ds(1) = (1.0_dp + strain(1))*rpd(1) + &
                  0.5_dp*(strain(6)*rpd(6) + strain(5)*rpd(5)) + xyzloc(1)*xcom
      dr2ds(2) = (1.0_dp + strain(2))*rpd(2) + &
                  0.5_dp*(strain(6)*rpd(6) + strain(4)*rpd(4)) + xyzloc(2)*ycom
      dr2ds(3) = (1.0_dp + strain(3))*rpd(3) + &
                  0.5_dp*(strain(5)*rpd(5) + strain(4)*rpd(4)) + xyzloc(3)*zcom
      dr2ds(4) = (1.0_dp + 0.5_dp*(strain(2) + strain(3)))*rpd(4) + &
                  0.25_dp*(strain(4)*(rpd(2) + rpd(3)) + &
                           strain(5)*rpd(6) + strain(6)*rpd(5)) + &
                      0.5_dp*(xyzloc(2)*zcom + xyzloc(3)*ycom)
      dr2ds(5) = (1.0_dp + 0.5_dp*(strain(1) + strain(3)))*rpd(5) + &
                  0.25_dp*(strain(5)*(rpd(1) + rpd(3)) + &
                           strain(4)*rpd(6) + strain(6)*rpd(4)) + &
                      0.5_dp*(xyzloc(1)*zcom + xyzloc(3)*xcom)
      dr2ds(6) = (1.0_dp + 0.5_dp*(strain(1) + strain(2)))*rpd(6) + &
                  0.25_dp*(strain(6)*(rpd(1) + rpd(2)) + &
                           strain(4)*rpd(5) + strain(5)*rpd(4)) + &
                      0.5_dp*(xyzloc(1)*ycom + xyzloc(2)*xcom)
      if (lgrad2) then
!
!  Initialise second derivatives
!
        d2r2ds2(1:6,1:6) = 0.0_dp
!
!  Second derivatives with respect to 2 strains
!
        d2r2ds2(1,1) = rpd(1)
        d2r2ds2(5,1) = 0.5_dp*rpd(5)
        d2r2ds2(6,1) = 0.5_dp*rpd(6)
!
        d2r2ds2(2,2) = rpd(2)
        d2r2ds2(4,2) = 0.5_dp*rpd(4)
        d2r2ds2(6,2) = 0.5_dp*rpd(6)
!
        d2r2ds2(3,3) = rpd(3)
        d2r2ds2(4,3) = 0.5_dp*rpd(4)
        d2r2ds2(5,3) = 0.5_dp*rpd(5)
!
        d2r2ds2(4,4) = 0.25_dp*(rpd(2) + rpd(3))
        d2r2ds2(5,4) = 0.25_dp*rpd(6)
        d2r2ds2(6,4) = 0.25_dp*rpd(5)
!
        d2r2ds2(5,5) = 0.25_dp*(rpd(1) + rpd(3))
        d2r2ds2(6,5) = 0.25_dp*rpd(4)
!
        d2r2ds2(6,6) = 0.25_dp*(rpd(1) + rpd(2))
!
        do i = 2,6
          do j = 1,i-1
            d2r2ds2(j,i) = d2r2ds2(i,j)
          enddo
        enddo
!
!  Second derivatives with respect to 1 strain and a Cartesian coordinate
!
        d2r2dsdx(1,1) = 2.0_dp*xyzloc(1)*(1.0_dp + strain(1)) + &
                        0.5_dp*(strain(6)*xyzloc(2) + strain(5)*xyzloc(3)) + xcom
        d2r2dsdx(1,2) = 0.5_dp*strain(6)*xyzloc(1)
        d2r2dsdx(1,3) = 0.5_dp*strain(5)*xyzloc(1)
!
        d2r2dsdx(2,1) = 0.5_dp*strain(6)*xyzloc(2)
        d2r2dsdx(2,2) = 2.0_dp*xyzloc(2)*(1.0_dp + strain(2)) + &
                        0.5_dp*(strain(6)*xyzloc(1) + strain(4)*xyzloc(3)) + ycom
        d2r2dsdx(2,3) = 0.5_dp*strain(4)*xyzloc(2)
!
        d2r2dsdx(3,1) = 0.5_dp*strain(5)*xyzloc(3)
        d2r2dsdx(3,2) = 0.5_dp*strain(4)*xyzloc(3)
        d2r2dsdx(3,3) = 2.0_dp*xyzloc(3)*(1.0_dp + strain(3)) + &
                        0.5_dp*(strain(5)*xyzloc(1) + strain(4)*xyzloc(2)) + zcom
!
        d2r2dsdx(4,1) = 0.25_dp*(strain(5)*xyzloc(2) + strain(6)*xyzloc(3))
        d2r2dsdx(4,2) = (1.0_dp + 0.5_dp*(strain(2) + strain(3)))*xyzloc(3) + &
                         0.25_dp*(2.0_dp*strain(4)*xyzloc(2) + strain(5)*xyzloc(1)) + 0.5_dp*zcom
        d2r2dsdx(4,3) = (1.0_dp + 0.5_dp*(strain(2) + strain(3)))*xyzloc(2) + &
                         0.25_dp*(2.0_dp*strain(4)*xyzloc(3) + strain(6)*xyzloc(1)) + 0.5_dp*ycom
!
        d2r2dsdx(5,1) = (1.0_dp + 0.5_dp*(strain(1) + strain(3)))*xyzloc(3) + &
                         0.25_dp*(2.0_dp*strain(5)*xyzloc(1) + strain(4)*xyzloc(2)) + 0.5_dp*zcom
        d2r2dsdx(5,2) = 0.25_dp*(strain(4)*xyzloc(1) + strain(6)*xyzloc(3))
        d2r2dsdx(5,3) = (1.0_dp + 0.5_dp*(strain(1) + strain(3)))*xyzloc(1) + &
                         0.25_dp*(2.0_dp*strain(5)*xyzloc(3) + strain(6)*xyzloc(2)) + 0.5_dp*xcom
!
        d2r2dsdx(6,1) = (1.0_dp + 0.5_dp*(strain(1) + strain(2)))*xyzloc(2) + &
                         0.25_dp*(2.0_dp*strain(6)*xyzloc(1) + strain(4)*xyzloc(3)) + 0.5_dp*ycom
        d2r2dsdx(6,2) = (1.0_dp + 0.5_dp*(strain(1) + strain(2)))*xyzloc(1) + &
                         0.25_dp*(2.0_dp*strain(6)*xyzloc(2) + strain(5)*xyzloc(3)) + 0.5_dp*xcom
        d2r2dsdx(6,3) = 0.25_dp*(strain(4)*xyzloc(1) + strain(5)*xyzloc(2))
!
!  Transform mixed derivatives by inverse strain matrix
!
        do i = 1,6
          tmp(1:3) = 0.0_dp
          do j = 1,3
            do l = 1,3
              tmp(j) = tmp(j) + d2r2dsdx(i,l)*straininverse(l,j)
            enddo
          enddo
          d2r2dsdx(i,1:3) = tmp(1:3)
        enddo
      endif
    else
      dr2ds(1:6) = d2r2dx2(1:6)
      if (lgrad2) then
!
!  Initialise second derivatives
!
        d2r2ds2(1:6,1:6) = 0.0_dp
        d2r2dsdx(1:6,1:3) = 0.0_dp
!
!  Second derivatives with respect to 2 strains
!
        d2r2ds2(1,1) = (xvec - xcom)*((xvec - xcom) + xvec)
        d2r2ds2(5,1) = 0.5_dp*((zvec - zcom)*(xvec - xcom) + xvec*(zvec - zcom))
        d2r2ds2(6,1) = 0.5_dp*((yvec - ycom)*(xvec - xcom) + xvec*(yvec - ycom))
!
        d2r2ds2(2,2) = (yvec - ycom)*((yvec - ycom) + yvec)
        d2r2ds2(4,2) = 0.5_dp*((zvec - zcom)*(yvec - ycom) + yvec*(zvec - zcom))
        d2r2ds2(6,2) = 0.5_dp*((xvec - xcom)*(yvec - ycom) + yvec*(xvec - xcom))
!
        d2r2ds2(3,3) = (zvec - zcom)*((zvec - zcom) + zvec)
        d2r2ds2(4,3) = 0.5_dp*((yvec - ycom)*(zvec - zcom) + zvec*(yvec - ycom))
        d2r2ds2(5,3) = 0.5_dp*((xvec - xcom)*(zvec - zcom) + zvec*(xvec - xcom))
!
        d2r2ds2(4,4) = 0.25_dp*((yvec - ycom)**2 + (zvec - zcom)**2 + &
                                yvec*(yvec - ycom) + zvec*(zvec - zcom))
        d2r2ds2(5,4) = 0.25_dp*((xvec - xcom)*(yvec - ycom) + yvec*(xvec - xcom))
        d2r2ds2(6,4) = 0.25_dp*((xvec - xcom)*(zvec - zcom) + zvec*(xvec - xcom))
!
        d2r2ds2(5,5) = 0.25_dp*((xvec - xcom)**2 + (zvec - zcom)**2 + &
                                xvec*(xvec - xcom) + zvec*(zvec - zcom))
        d2r2ds2(6,5) = 0.25_dp*((yvec - ycom)*(zvec - zcom) + zvec*(yvec - ycom))
!
        d2r2ds2(6,6) = 0.25_dp*((xvec - xcom)**2 + (yvec - ycom)**2 + &
                                xvec*(xvec - xcom) + yvec*(yvec - ycom))
!
        do i = 2,6
          do j = 1,i-1
            d2r2ds2(j,i) = d2r2ds2(i,j)
          enddo
        enddo
!
!  Second derivatives with respect to 1 strain and a Cartesian coordinate
!
        d2r2dsdx(1,1) = xvec + (xvec - xcom)
        d2r2dsdx(2,2) = yvec + (yvec - ycom)
        d2r2dsdx(3,3) = zvec + (zvec - zcom)
!
        d2r2dsdx(4,2) = 0.5_dp*(zvec + (zvec - zcom))
        d2r2dsdx(4,3) = 0.5_dp*(yvec + (yvec - ycom))
!
        d2r2dsdx(5,1) = 0.5_dp*(zvec + (zvec - zcom))
        d2r2dsdx(5,3) = 0.5_dp*(xvec + (xvec - xcom))
!
        d2r2dsdx(6,1) = 0.5_dp*(yvec + (yvec - ycom))
        d2r2dsdx(6,2) = 0.5_dp*(xvec + (xvec - xcom))
      endif
    endif
  elseif (ndim.eq.2) then
    if (lfinitestrain) then
!
!  Create vector scaled by inverse strain matrix
!
      xyzloc(1) = straininverse(1,1)*(xvec - xcom) + &
                  straininverse(2,1)*(yvec - ycom)
      xyzloc(2) = straininverse(1,2)*(xvec - xcom) + &
                  straininverse(2,2)*(yvec - ycom)
!
      rpd(1) = xyzloc(1)*xyzloc(1)
      rpd(2) = xyzloc(2)*xyzloc(2)
      rpd(3) = xyzloc(1)*xyzloc(2)
!
      dr2ds(1) = (1.0_dp + strain(1))*rpd(1) + &
                  0.5_dp*strain(3)*rpd(3) + xyzloc(1)*xcom
      dr2ds(2) = (1.0_dp + strain(2))*rpd(2) + &
                  0.5_dp*strain(3)*rpd(3) + xyzloc(2)*ycom
      dr2ds(6) = (1.0_dp + 0.5_dp*(strain(1) + strain(2)))*rpd(3) + &
                  0.25_dp*strain(3)*(rpd(1) + rpd(2)) + &
                  0.5_dp*(xyzloc(1)*ycom + xyzloc(2)*xcom)
      dr2ds(3:5) = d2r2dx2(3:5)
      if (lgrad2) then
!
!  Initialise second derivatives
!
        d2r2ds2(1:6,1:6) = 0.0_dp
        d2r2dsdx(1:6,1:3) = 0.0_dp
!
!  Second derivatives with respect to 2 strains
!
        d2r2ds2(1,1) = rpd(1)
        d2r2ds2(6,1) = 0.5_dp*rpd(3)
        d2r2ds2(1,6) = 0.5_dp*rpd(3)
!
        d2r2ds2(2,2) = rpd(2)
        d2r2ds2(6,2) = 0.5_dp*rpd(3)
        d2r2ds2(2,6) = 0.5_dp*rpd(3)
!
        d2r2ds2(6,6) = 0.25_dp*(rpd(1) + rpd(2))
!
!  Second derivatives with respect to 1 strain and a Cartesian coordinate
!
        d2r2dsdx(1,1) = 2.0_dp*xyzloc(1)*(1.0_dp + strain(1)) + &
                        0.5_dp*strain(3)*xyzloc(2) + xcom
        d2r2dsdx(1,2) = 0.5_dp*strain(3)*xyzloc(1)
!
        d2r2dsdx(2,1) = 0.5_dp*strain(3)*xyzloc(2)
        d2r2dsdx(2,2) = 2.0_dp*xyzloc(2)*(1.0_dp + strain(2)) + &
                        0.5_dp*strain(3)*xyzloc(1) + ycom
!
        d2r2dsdx(6,1) = (1.0_dp + 0.5_dp*(strain(1) + strain(2)))*xyzloc(2) + &
                         0.5_dp*strain(3)*xyzloc(1) + 0.5_dp*ycom
        d2r2dsdx(6,2) = (1.0_dp + 0.5_dp*(strain(1) + strain(2)))*xyzloc(1) + &
                         0.5_dp*strain(3)*xyzloc(2) + 0.5_dp*xcom
!
!  Transform mixed derivatives by inverse strain matrix
!
        do i = 1,6
          tmp(1:2) = 0.0_dp
          do j = 1,2
            do l = 1,2
              tmp(j) = tmp(j) + d2r2dsdx(i,l)*straininverse(l,j)
            enddo
          enddo
          d2r2dsdx(i,1:2) = tmp(1:2)
        enddo
      endif
    else
      dr2ds(1:6) = d2r2dx2(1:6)
      if (lgrad2) then
!
!  Initialise second derivatives
!
        d2r2ds2(1:6,1:6) = 0.0_dp
        d2r2dsdx(1:6,1:3) = 0.0_dp
!
!  Second derivatives with respect to 2 strains
!
        d2r2ds2(1,1) = (xvec - xcom)*((xvec - xcom) + xvec)
        d2r2ds2(6,1) = 0.5_dp*((yvec - ycom)*(xvec - xcom) + xvec*(yvec - ycom))
        d2r2ds2(1,6) = 0.5_dp*((yvec - ycom)*(xvec - xcom) + xvec*(yvec - ycom))
!
        d2r2ds2(2,2) = (yvec - ycom)*((yvec - ycom) + yvec)
        d2r2ds2(6,2) = 0.5_dp*((xvec - xcom)*(yvec - ycom) + yvec*(xvec - xcom))
        d2r2ds2(2,6) = 0.5_dp*((xvec - xcom)*(yvec - ycom) + yvec*(xvec - xcom))
!
        d2r2ds2(6,6) = 0.25_dp*((xvec - xcom)**2 + (yvec - ycom)**2 + &
                                xvec*(xvec - xcom) + yvec*(yvec - ycom))
!
!  Second derivatives with respect to 1 strain and a Cartesian coordinate
!
        d2r2dsdx(1,1) = xvec + (xvec - xcom)
        d2r2dsdx(2,2) = yvec + (yvec - ycom)
        d2r2dsdx(6,1) = 0.5_dp*(yvec + (yvec - ycom))
        d2r2dsdx(6,2) = 0.5_dp*(xvec + (xvec - xcom))
      endif
    endif
  elseif (ndim.eq.1) then
    if (lfinitestrain) then
      xyzloc(1) = (xvec - xcom)*straininverse(1,1)
      rpd(1) = xyzloc(1)*xyzloc(1)
      dr2ds(1) = (1.0_dp + strain(1))*rpd(1) + xyzloc(1)*xcom
      if (lgrad2) then
        d2r2ds2(1,1) = 0.0_dp
        d2r2dsdx(1,1:3) = 0.0_dp
!
!  Second derivatives with respect to 2 strains
!
        d2r2ds2(1,1) = rpd(1)
!
!  Second derivatives with respect to 1 strain and a Cartesian coordinate
!
        d2r2dsdx(1,1) = 2.0_dp*xyzloc(1)*(1.0_dp + strain(1)) + xcom
!
!  Transform mixed derivatives by inverse strain matrix
!
        d2r2dsdx(1,1) = d2r2dsdx(1,1)*straininverse(1,1)
      endif
    else
      dr2ds(1) = d2r2dx2(1)
      if (lgrad2) then
!
!  Second derivatives with respect to 2 strains
!
        d2r2ds2(1,1) = (xvec - xcom)*((xvec - xcom) + xvec)
!
!  Second derivatives with respect to 1 strain and a Cartesian coordinate
!
        d2r2dsdx(1,1) = xvec + (xvec - xcom)
        d2r2dsdx(1,2) = 0.0_dp
        d2r2dsdx(1,3) = 0.0_dp
      endif
    endif
  endif
!
!  Reset Cartesian second derivatives without centre of mass subtraction
!
  if (lrigid.and.lgrad2) then
    d2r2dx2(1) = xvec*xvec
    d2r2dx2(2) = yvec*yvec
    d2r2dx2(3) = zvec*zvec
    d2r2dx2(4) = yvec*zvec
    d2r2dx2(5) = xvec*zvec
    d2r2dx2(6) = xvec*yvec
  endif
#ifdef TRACE
  call trace_out('real1strterm')
#endif
!
  return
  end

  subroutine vecprestrain(ndim,xvec,yvec,zvec)
!
!  Convert a vector to its form prior to strain being applied.
!
!  On entry :
!
!   ndim              = number of dimensions
!   xvec, yvec, zvec  = Cartesian components of interatomic vectors
!
!  On exit :
!
!   xvec, yvec, zvec  = Cartesian components of interatomic vectors pre-strain
!  
!  10/18 Created from real1strterm
!  11/18 Finite strain flag introduced instead of lstraincell
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
!  Julian Gale, CIC, Curtin University, November 2018
!
  use datatypes
  use derivatives,   only : lfinitestrain
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)     :: ndim
  real(dp),    intent(inout)  :: xvec
  real(dp),    intent(inout)  :: yvec
  real(dp),    intent(inout)  :: zvec
!
!  Local variables
!
  real(dp)                    :: xtmp
  real(dp)                    :: ytmp
  real(dp)                    :: ztmp
#ifdef TRACE
  call trace_in('vecprestrain')
#endif
!
  if (.not.lfinitestrain) return
!
  if (ndim.eq.3) then
    xtmp = xvec
    ytmp = yvec
    ztmp = zvec
!
!  Create vector scaled by inverse strain matrix
!
    xvec = straininverse(1,1)*xtmp + &
           straininverse(2,1)*ytmp + &
           straininverse(3,1)*ztmp
    yvec = straininverse(1,2)*xtmp + &
           straininverse(2,2)*ytmp + &
           straininverse(3,2)*ztmp
    zvec = straininverse(1,3)*xtmp + &
           straininverse(2,3)*ytmp + &
           straininverse(3,3)*ztmp
  elseif (ndim.eq.2) then
    xtmp = xvec
    ytmp = yvec
!
!  Create vector scaled by inverse strain matrix
!
    xvec = straininverse(1,1)*xtmp + &
           straininverse(2,1)*ytmp
    yvec = straininverse(1,2)*xtmp + &
           straininverse(2,2)*ytmp
  elseif (ndim.eq.1) then
    xvec = xvec*straininverse(1,1)
  endif
#ifdef TRACE
  call trace_out('vecprestrain')
#endif
!
  return
  end

  subroutine cartstrterm(ndim,xvec,yvec,zvec,xcom,ycom,zcom,dxyzds,d2xyzdsdx,d2xyzds2,lgrad2)
!
!  Subroutine for calculating real space two-body strain terms
!  for Cartesian components
!
!  On entry :
!
!   ndim              = number of dimensions
!   xvec, yvec, zvec  = Cartesian components of interatomic vectors
!   xcom, ycom, zcom  = Cartesian components of vectors to centre of mass of rigid molecule
!   lgrad2            = if true then compute d2xyzdsdx and d2xyzds2
!
!  On exit :
!
!   dxyzds            = first derivative of x,y,z with respect to strain
!   d2xyzdsdx         = second derivative of x,y,z with respect to strains and a Cartesian component for fractionals (if lgrad2 is true)
!   d2xyzds2          = second derivative of x,y,z with respect to two strains (if lgrad2 is true)
!  
!   1/19 Created from real1strterm
!   5/19 d2xyzds2 added
!   7/19 Correction to finite strain second derivatives
!  11/19 Rigid molecules modifications added
!   3/20 Corrections to 2-D strain terms
!   4/20 d2xyzdsdc added
!   4/20 derv3c changes reversed as they are no longer required
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
!  Julian Gale, CIC, Curtin University, April 2020
!
  use datatypes
  use control,       only : lrigid
  use derivatives,   only : lfinitestrain
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)  :: ndim
  logical,     intent(in)  :: lgrad2
  real(dp),    intent(out) :: dxyzds(6,3)
  real(dp),    intent(out) :: d2xyzdsdx(6,3,3)
  real(dp),    intent(out) :: d2xyzds2(6,6,3)
  real(dp),    intent(in)  :: xvec
  real(dp),    intent(in)  :: yvec
  real(dp),    intent(in)  :: zvec
  real(dp),    intent(in)  :: xcom
  real(dp),    intent(in)  :: ycom
  real(dp),    intent(in)  :: zcom
!
!  Local variables
!
  real(dp)                 :: xlp
  real(dp)                 :: ylp
  real(dp)                 :: zlp
#ifdef TRACE
  call trace_in('cartstrterm')
#endif
!
!  For straincell algorithm create pre-strain values of Cartesian coordinates
!
  if (lrigid) then
    xlp = xvec - xcom
    ylp = yvec - ycom
    zlp = zvec - zcom
  else
    xlp = xvec
    ylp = yvec
    zlp = zvec
  endif
!
  if (lfinitestrain) then
    call vecprestrain(ndim,xlp,ylp,zlp)
  endif
!
!  First derivatives with respect to strain of x, y and z
!
  if (ndim.eq.3) then
    dxyzds(1:6,1:3) = 0.0_dp
    dxyzds(1,1) = xlp
    dxyzds(5,1) = 0.5_dp*zlp
    dxyzds(6,1) = 0.5_dp*ylp
    dxyzds(2,2) = ylp
    dxyzds(4,2) = 0.5_dp*zlp
    dxyzds(6,2) = 0.5_dp*xlp
    dxyzds(3,3) = zlp
    dxyzds(4,3) = 0.5_dp*ylp
    dxyzds(5,3) = 0.5_dp*xlp
  elseif (ndim.eq.2) then
    dxyzds(1:3,1:3) = 0.0_dp
    dxyzds(1,1) = xlp
    dxyzds(3,1) = 0.5_dp*ylp
    dxyzds(2,2) = ylp
    dxyzds(3,2) = 0.5_dp*xlp
  elseif (ndim.eq.1) then
    dxyzds(1,1:3) = 0.0_dp
    dxyzds(1,1) = xlp
  endif
  if (lgrad2) then
    if (lfinitestrain) then
!
!  For finite strains the second derivatives with respect to two strains are zero
!
      if (ndim.eq.3) then
        d2xyzds2(1:6,1:6,1:3) = 0.0_dp
        d2xyzdsdx(1:6,1:3,1:3) = 0.0_dp
!
        d2xyzdsdx(1,1,1) = straininverse(1,1)
        d2xyzdsdx(1,1,2) = straininverse(2,1)
        d2xyzdsdx(1,1,3) = straininverse(3,1)
        d2xyzdsdx(5,1,1) = 0.5_dp*straininverse(1,3)
        d2xyzdsdx(5,1,2) = 0.5_dp*straininverse(2,3)
        d2xyzdsdx(5,1,3) = 0.5_dp*straininverse(3,3)
        d2xyzdsdx(6,1,1) = 0.5_dp*straininverse(1,2)
        d2xyzdsdx(6,1,2) = 0.5_dp*straininverse(2,2)
        d2xyzdsdx(6,1,3) = 0.5_dp*straininverse(3,2)
        d2xyzdsdx(2,2,1) = straininverse(1,2)
        d2xyzdsdx(2,2,2) = straininverse(2,2)
        d2xyzdsdx(2,2,3) = straininverse(3,2)
        d2xyzdsdx(4,2,1) = 0.5_dp*straininverse(1,3)
        d2xyzdsdx(4,2,2) = 0.5_dp*straininverse(2,3)
        d2xyzdsdx(4,2,3) = 0.5_dp*straininverse(3,3)
        d2xyzdsdx(6,2,1) = 0.5_dp*straininverse(1,1)
        d2xyzdsdx(6,2,2) = 0.5_dp*straininverse(2,1)
        d2xyzdsdx(6,2,3) = 0.5_dp*straininverse(3,1)
        d2xyzdsdx(3,3,1) = straininverse(1,3)
        d2xyzdsdx(3,3,2) = straininverse(2,3)
        d2xyzdsdx(3,3,3) = straininverse(3,3)
        d2xyzdsdx(4,3,1) = 0.5_dp*straininverse(1,2)
        d2xyzdsdx(4,3,2) = 0.5_dp*straininverse(2,2)
        d2xyzdsdx(4,3,3) = 0.5_dp*straininverse(3,2)
        d2xyzdsdx(5,3,1) = 0.5_dp*straininverse(1,1)
        d2xyzdsdx(5,3,2) = 0.5_dp*straininverse(2,1)
        d2xyzdsdx(5,3,3) = 0.5_dp*straininverse(3,1)
      elseif (ndim.eq.2) then
        d2xyzds2(1:3,1:6,1:3) = 0.0_dp
        d2xyzdsdx(1:3,1:3,1:3) = 0.0_dp
!
        d2xyzdsdx(1,1,1) = straininverse(1,1)
        d2xyzdsdx(1,1,2) = straininverse(2,1)
        d2xyzdsdx(3,1,1) = 0.5_dp*straininverse(1,2)
        d2xyzdsdx(3,1,2) = 0.5_dp*straininverse(2,2)
        d2xyzdsdx(2,2,1) = straininverse(1,2)
        d2xyzdsdx(2,2,2) = straininverse(2,2)
        d2xyzdsdx(3,2,1) = 0.5_dp*straininverse(1,1)
        d2xyzdsdx(3,2,2) = 0.5_dp*straininverse(2,1)
      elseif (ndim.eq.1) then
        d2xyzds2(1,1:6,1:3) = 0.0_dp
        d2xyzdsdx(1,1:3,1:3) = 0.0_dp
!
        d2xyzdsdx(1,1,1) = straininverse(1,1)
      endif
    else
!
!  Second derivatives with respect to Cartesian component and strain of x, y and z
!
      if (ndim.eq.3) then
        d2xyzdsdx(1:6,1:3,1:3) = 0.0_dp
        d2xyzdsdx(1,1,1) = 1.0_dp
        d2xyzdsdx(5,3,1) = 0.5_dp
        d2xyzdsdx(6,2,1) = 0.5_dp
        d2xyzdsdx(2,2,2) = 1.0_dp
        d2xyzdsdx(4,3,2) = 0.5_dp
        d2xyzdsdx(6,1,2) = 0.5_dp
        d2xyzdsdx(3,3,3) = 1.0_dp
        d2xyzdsdx(4,2,3) = 0.5_dp
        d2xyzdsdx(5,1,3) = 0.5_dp
      elseif (ndim.eq.2) then
        d2xyzdsdx(1:3,1:3,1:3) = 0.0_dp
        d2xyzdsdx(1,1,1) = 1.0_dp
        d2xyzdsdx(3,2,1) = 0.5_dp
        d2xyzdsdx(2,2,2) = 1.0_dp
        d2xyzdsdx(3,1,2) = 0.5_dp
      elseif (ndim.eq.1) then
        d2xyzdsdx(1,1:3,1:3) = 0.0_dp
        d2xyzdsdx(1,1,1) = 1.0_dp
      endif
!
!  Second derivatives with respect to two strains of x, y and z
!
      if (ndim.eq.3) then
        d2xyzds2(1:6,1:6,1:3) = 0.0_dp
!
        d2xyzds2(1,1,1) = xlp
        d2xyzds2(5,1,1) = 0.5_dp*zlp
        d2xyzds2(6,1,1) = 0.5_dp*ylp
!
        d2xyzds2(3,5,1) = 0.5_dp*zlp
        d2xyzds2(4,5,1) = 0.25_dp*ylp
        d2xyzds2(5,5,1) = 0.25_dp*xlp
!
        d2xyzds2(2,6,1) = 0.5_dp*ylp
        d2xyzds2(4,6,1) = 0.25_dp*zlp
        d2xyzds2(6,6,1) = 0.25_dp*xlp
!
        d2xyzds2(2,2,2) = ylp
        d2xyzds2(4,2,2) = 0.5_dp*zlp
        d2xyzds2(6,2,2) = 0.5_dp*xlp
!
        d2xyzds2(3,4,2) = 0.5_dp*zlp
        d2xyzds2(4,4,2) = 0.25_dp*ylp
        d2xyzds2(5,4,2) = 0.25_dp*xlp
!
        d2xyzds2(1,6,2) = 0.5_dp*xlp
        d2xyzds2(5,6,2) = 0.25_dp*zlp
        d2xyzds2(6,6,2) = 0.25_dp*ylp
!
        d2xyzds2(3,3,3) = zlp
        d2xyzds2(4,3,3) = 0.5_dp*ylp
        d2xyzds2(5,3,3) = 0.5_dp*xlp
!
        d2xyzds2(2,4,3) = 0.5_dp*ylp
        d2xyzds2(4,4,3) = 0.25_dp*zlp
        d2xyzds2(6,4,3) = 0.25_dp*xlp
!
        d2xyzds2(1,5,3) = 0.5_dp*xlp
        d2xyzds2(5,5,3) = 0.25_dp*zlp
        d2xyzds2(6,5,3) = 0.25_dp*ylp
      elseif (ndim.eq.2) then
        d2xyzds2(1:3,1:3,1:3) = 0.0_dp
!
        d2xyzds2(1,1,1) = xlp
        d2xyzds2(1,3,1) = 0.5_dp*ylp
!
        d2xyzds2(3,2,1) = 0.5_dp*ylp
        d2xyzds2(3,3,1) = 0.25_dp*xlp
!
        d2xyzds2(2,2,2) = ylp
        d2xyzds2(2,3,2) = 0.5_dp*xlp
!
        d2xyzds2(3,1,2) = 0.5_dp*xlp
        d2xyzds2(3,3,2) = 0.25_dp*ylp
      elseif (ndim.eq.1) then
        d2xyzds2(1,1,1:3) = 0.0_dp
!
        d2xyzds2(1,1,1) = xlp
      endif
    endif
  endif
#ifdef TRACE
  call trace_out('cartstrterm')
#endif
!
  return
  end

  subroutine gxyzterms(ndim,maxgvec,ngvec,xrk,yrk,zrk,dgds,d2gds2,lgrad2)
!
!  Subroutine for calculating strain derivatives terms in reciprocal space
!
!  On entry :
!
!   ndim              = number of dimensions
!   ngvec             = number of reciprocal space vectors to be processed
!   maxgvec           = maximum number of reciprocal space vectors (lower dimension of arrays)
!   xrk, yrk, zrk     = Cartesian components of reciprocal space vectors
!   lgrad2            = if true then compute d2gds2
!
!  On exit :
!
!   dgds              = first derivative of G components with respect to strain
!   d2g2ds2           = second derivative of G components with respect to two strains
!  
!   5/19 Created from gstrterms
!   7/19 Corrected for finite strain case
!   2/20 Strain second derivatives corrected for non-finite strain case
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
  use derivatives,   only : lfinitestrain
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)  :: ndim
  integer(i4), intent(in)  :: ngvec
  integer(i4), intent(in)  :: maxgvec
  logical,     intent(in)  :: lgrad2
  real(dp),    intent(out) :: dgds(maxgvec,3,6)
  real(dp),    intent(out) :: d2gds2(maxgvec,3,6,6)
  real(dp),    intent(in)  :: xrk(ngvec)
  real(dp),    intent(in)  :: yrk(ngvec)
  real(dp),    intent(in)  :: zrk(ngvec)
!
!  Local variables
!
  integer(i4)              :: i
  integer(i4)              :: j
  integer(i4)              :: k
  integer(i4)              :: l
  integer(i4)              :: m
  integer(i4)              :: ns1
  integer(i4)              :: ns2
  real(dp)                 :: gloc(3)
  real(dp)                 :: dIds(6,3,3)
  real(dp)                 :: d2Ids2(6,6,3,3)
#ifdef TRACE
  call trace_in('gxyzterms')
#endif
!
  if (ndim.eq.3) then
!
!  Initialise derivatives
!
    dgds(1:ngvec,1:3,1:6) = 0.0_dp
    if (lgrad2) then
      d2gds2(1:ngvec,1:3,1:6,1:6) = 0.0_dp
    endif
!
    if (lfinitestrain) then
!
!  Compute inverse matrix
!
      straininv(1,1) = strainfull(2,2)*strainfull(3,3) - strainfull(3,2)**2
      straininv(2,1) = strainfull(3,2)*strainfull(3,1) - strainfull(2,1)*strainfull(3,3)
      straininv(3,1) = strainfull(3,2)*strainfull(2,1) - strainfull(3,1)*strainfull(2,2)
      straininv(2,2) = strainfull(1,1)*strainfull(3,3) - strainfull(3,1)**2
      straininv(3,2) = strainfull(3,1)*strainfull(2,1) - strainfull(3,2)*strainfull(1,1)
      straininv(3,3) = strainfull(1,1)*strainfull(2,2) - strainfull(2,1)**2
      straininv(1,2) = straininv(2,1)
      straininv(1,3) = straininv(3,1)
      straininv(2,3) = straininv(3,2)
!
!  Compute unique derivatives of inverse matrix with respect to strains
!  LHS = strain index; RHS = matrix element index
!
      dstraininvds(1:6,1:3,1:3) = 0.0_dp
!
      dstraininvds(2,1,1) = strainfull(3,3)
      dstraininvds(3,1,1) = strainfull(2,2)
      dstraininvds(4,1,1) = - strainfull(3,2)
!
      dstraininvds(1,2,2) = strainfull(3,3)
      dstraininvds(3,2,2) = strainfull(1,1)
      dstraininvds(5,2,2) = - strainfull(3,1)
!
      dstraininvds(1,3,3) = strainfull(2,2)
      dstraininvds(2,3,3) = strainfull(1,1)
      dstraininvds(6,3,3) = - strainfull(2,1)
!
      dstraininvds(1,3,2) = - strainfull(3,2)
      dstraininvds(4,3,2) = - 0.5_dp*strainfull(1,1)
      dstraininvds(5,3,2) = 0.5_dp*strainfull(2,1)
      dstraininvds(6,3,2) = 0.5_dp*strainfull(3,1)
!
      dstraininvds(2,3,1) = - strainfull(3,1)
      dstraininvds(4,3,1) = 0.5_dp*strainfull(2,1)
      dstraininvds(5,3,1) = - 0.5_dp*strainfull(2,2)
      dstraininvds(6,3,1) = 0.5_dp*strainfull(3,2)
!
      dstraininvds(3,2,1) = - strainfull(2,1)
      dstraininvds(4,2,1) = 0.5_dp*strainfull(3,1)
      dstraininvds(5,2,1) = 0.5_dp*strainfull(3,2)
      dstraininvds(6,2,1) = - 0.5_dp*strainfull(3,3)
!
      dstraininvds(3,1,2) = dstraininvds(3,2,1)
      dstraininvds(4,1,2) = dstraininvds(4,2,1)
      dstraininvds(5,1,2) = dstraininvds(5,2,1)
      dstraininvds(6,1,2) = dstraininvds(6,2,1)
!
      dstraininvds(2,1,3) = dstraininvds(2,3,1)
      dstraininvds(4,1,3) = dstraininvds(4,3,1)
      dstraininvds(5,1,3) = dstraininvds(5,3,1)
      dstraininvds(6,1,3) = dstraininvds(6,3,1)
!
      dstraininvds(1,2,3) = dstraininvds(1,3,2)
      dstraininvds(4,2,3) = dstraininvds(4,3,2)
      dstraininvds(5,2,3) = dstraininvds(5,3,2)
      dstraininvds(6,2,3) = dstraininvds(6,3,2)
!
!  Create overall strain derivatives of inverse matrix
!
      do i = 1,3
        do j = 1,3
          do k = 1,6
            dIds(k,j,i) = straindet*(dstraininvds(k,j,i) - straindet*strainddetds(k)*straininv(j,i))
          enddo
        enddo
      enddo
      if (lgrad2) then
!
!  Compute second derivatives of inverse strain matrix
!
        d2straininvds2(1:6,1:6,1:3,1:3) = 0.0_dp
!
        d2straininvds2(3,2,1,1) = 1.0_dp
        d2straininvds2(2,3,1,1) = 1.0_dp
        d2straininvds2(4,4,1,1) = - 0.5_dp
!
        d2straininvds2(3,1,2,2) = 1.0_dp
        d2straininvds2(1,3,2,2) = 1.0_dp
        d2straininvds2(5,5,2,2) = - 0.5_dp
!
        d2straininvds2(2,1,3,3) = 1.0_dp
        d2straininvds2(1,2,3,3) = 1.0_dp
        d2straininvds2(6,6,3,3) = - 0.5_dp
!
        d2straininvds2(4,1,3,2) = - 0.5_dp
        d2straininvds2(1,4,3,2) = - 0.5_dp
        d2straininvds2(6,5,3,2) = 0.25_dp
        d2straininvds2(5,6,3,2) = 0.25_dp
!
        d2straininvds2(4,1,2,3) = - 0.5_dp
        d2straininvds2(1,4,2,3) = - 0.5_dp
        d2straininvds2(6,5,2,3) = 0.25_dp
        d2straininvds2(5,6,2,3) = 0.25_dp
!
        d2straininvds2(5,2,3,1) = - 0.5_dp
        d2straininvds2(6,4,3,1) = 0.25_dp
        d2straininvds2(2,5,3,1) = - 0.5_dp
        d2straininvds2(4,6,3,1) = 0.25_dp
!
        d2straininvds2(5,2,1,3) = - 0.5_dp
        d2straininvds2(6,4,1,3) = 0.25_dp
        d2straininvds2(2,5,1,3) = - 0.5_dp
        d2straininvds2(4,6,1,3) = 0.25_dp
!
        d2straininvds2(6,3,2,1) = - 0.5_dp
        d2straininvds2(5,4,2,1) = 0.25_dp
        d2straininvds2(4,5,2,1) = 0.25_dp
        d2straininvds2(3,6,2,1) = - 0.5_dp
!
        d2straininvds2(6,3,1,2) = - 0.5_dp
        d2straininvds2(5,4,1,2) = 0.25_dp
        d2straininvds2(4,5,1,2) = 0.25_dp
        d2straininvds2(3,6,1,2) = - 0.5_dp
!
!  Create overall strain second derivatives of inverse matrix
!
        do i = 1,3
          do j = 1,3
            do l = 1,6
              do m = 1,6
                d2Ids2(m,l,j,i) = straindet*(d2straininvds2(m,l,j,i) &
                                  + 2.0_dp*straindet*straindet*strainddetds(m)*strainddetds(l)*straininv(j,i) &
                                  - straindet*straind2detds2(m,l)*straininv(j,i) &
                                  - straindet*strainddetds(l)*dstraininvds(m,j,i) &
                                  - straindet*strainddetds(m)*dstraininvds(l,j,i))
              enddo
            enddo
          enddo
        enddo
      endif
!
      do k = 1,ngvec
!
!  Create vector scaled by inverse strain matrix
!
        gloc(1) = strainfull(1,1)*xrk(k) + &
                  strainfull(2,1)*yrk(k) + &
                  strainfull(3,1)*zrk(k)
        gloc(2) = strainfull(1,2)*xrk(k) + &
                  strainfull(2,2)*yrk(k) + &
                  strainfull(3,2)*zrk(k)
        gloc(3) = strainfull(1,3)*xrk(k) + &
                  strainfull(2,3)*yrk(k) + &
                  strainfull(3,3)*zrk(k)
!
        do i = 1,3
          do j = 1,6
            dgds(k,1,j) = dgds(k,1,j) + dIds(j,1,i)*gloc(i)
            dgds(k,2,j) = dgds(k,2,j) + dIds(j,2,i)*gloc(i)
            dgds(k,3,j) = dgds(k,3,j) + dIds(j,3,i)*gloc(i)
          enddo
        enddo
      enddo
    else
      do k = 1,ngvec
        dgds(k,1,1) = - xrk(k)
        dgds(k,1,5) = - 0.5_dp*zrk(k)
        dgds(k,1,6) = - 0.5_dp*yrk(k)
        dgds(k,2,2) = - yrk(k)
        dgds(k,2,4) = - 0.5_dp*zrk(k)
        dgds(k,2,6) = - 0.5_dp*xrk(k)
        dgds(k,3,3) = - zrk(k)
        dgds(k,3,4) = - 0.5_dp*yrk(k)
        dgds(k,3,5) = - 0.5_dp*xrk(k)
      enddo
    endif
    if (lgrad2) then
      if (lfinitestrain) then
        do k = 1,ngvec
!
!  Create vector scaled by inverse strain matrix
!
          gloc(1) = strainfull(1,1)*xrk(k) + &
                    strainfull(2,1)*yrk(k) + &
                    strainfull(3,1)*zrk(k)
          gloc(2) = strainfull(1,2)*xrk(k) + &
                    strainfull(2,2)*yrk(k) + &
                    strainfull(3,2)*zrk(k)
          gloc(3) = strainfull(1,3)*xrk(k) + &
                    strainfull(2,3)*yrk(k) + &
                    strainfull(3,3)*zrk(k)
!
          do i = 1,3
            do ns1 = 1,6
              do ns2 = 1,6
                d2gds2(k,1,ns2,ns1) = d2gds2(k,1,ns2,ns1) + d2Ids2(ns2,ns1,1,i)*gloc(i)
                d2gds2(k,2,ns2,ns1) = d2gds2(k,2,ns2,ns1) + d2Ids2(ns2,ns1,2,i)*gloc(i)
                d2gds2(k,3,ns2,ns1) = d2gds2(k,3,ns2,ns1) + d2Ids2(ns2,ns1,3,i)*gloc(i)
              enddo
            enddo
          enddo
        enddo
      else
        do k = 1,ngvec
          d2gds2(k,1,1,1) = xrk(k)
          d2gds2(k,1,5,1) = 0.5_dp*zrk(k)
          d2gds2(k,1,6,1) = 0.5_dp*yrk(k)
!
          d2gds2(k,1,3,5) = 0.5_dp*zrk(k)
          d2gds2(k,1,4,5) = 0.25_dp*yrk(k)
          d2gds2(k,1,5,5) = 0.25_dp*xrk(k)
!
          d2gds2(k,1,2,6) = 0.5_dp*yrk(k)
          d2gds2(k,1,4,6) = 0.25_dp*zrk(k)
          d2gds2(k,1,6,6) = 0.25_dp*xrk(k)
!
          d2gds2(k,2,2,2) = yrk(k)
          d2gds2(k,2,4,2) = 0.5_dp*zrk(k)
          d2gds2(k,2,6,2) = 0.5_dp*xrk(k)
!
          d2gds2(k,2,3,4) = 0.5_dp*zrk(k)
          d2gds2(k,2,4,4) = 0.25_dp*yrk(k)
          d2gds2(k,2,5,4) = 0.25_dp*xrk(k)
!
          d2gds2(k,2,1,6) = 0.5_dp*xrk(k)
          d2gds2(k,2,5,6) = 0.25_dp*zrk(k)
          d2gds2(k,2,6,6) = 0.25_dp*yrk(k)
!
          d2gds2(k,3,3,3) = zrk(k)
          d2gds2(k,3,4,3) = 0.5_dp*yrk(k)
          d2gds2(k,3,5,3) = 0.5_dp*xrk(k)
!
          d2gds2(k,3,2,4) = 0.5_dp*yrk(k)
          d2gds2(k,3,4,4) = 0.25_dp*zrk(k)
          d2gds2(k,3,6,4) = 0.25_dp*xrk(k)
!
          d2gds2(k,3,1,5) = 0.5_dp*xrk(k)
          d2gds2(k,3,5,5) = 0.25_dp*zrk(k)
          d2gds2(k,3,6,5) = 0.25_dp*yrk(k)
        enddo
      endif
    endif
  elseif (ndim.eq.2) then
!
!  Initialise derivatives
!
    dgds(1:ngvec,1:3,1:3) = 0.0_dp
    if (lgrad2) then
      d2gds2(1:ngvec,1:3,1:3,1:3) = 0.0_dp
    endif
!
    if (lfinitestrain) then
!
!  Compute inverse matrix
!
      straininv(1,1) = strainfull(2,2)
      straininv(2,1) = - strainfull(2,1)
      straininv(1,2) = - strainfull(2,1)
      straininv(2,2) = strainfull(1,1)
!
!  Compute unique derivatives of inverse matrix with respect to strains
!  LHS = strain index; RHS = matrix element index
!
      dstraininvds(1:3,1:2,1:2) = 0.0_dp
!
      dstraininvds(2,1,1) = 1.0_dp
      dstraininvds(1,2,2) = 1.0_dp
      dstraininvds(3,2,1) = - 0.5_dp
      dstraininvds(3,1,2) = - 0.5_dp
!
!  Create overall strain derivatives of inverse matrix
!
      do i = 1,2
        do j = 1,2
          do k = 1,3
            dIds(k,j,i) = straindet*(dstraininvds(k,j,i) - straindet*strainddetds(k)*straininv(j,i))
          enddo
        enddo
      enddo
      if (lgrad2) then
!
!  Create overall strain second derivatives of inverse matrix
!
        do i = 1,2
          do j = 1,2
            do l = 1,3
              do m = 1,3
                d2Ids2(m,l,j,i) = straindet*( &
                                  + 2.0_dp*straindet*straindet*strainddetds(m)*strainddetds(l)*straininv(j,i) &
                                  - straindet*straind2detds2(m,l)*straininv(j,i) &
                                  - straindet*strainddetds(l)*dstraininvds(m,j,i) &
                                  - straindet*strainddetds(m)*dstraininvds(l,j,i))
              enddo
            enddo
          enddo
        enddo
      endif
    endif
!
    if (lfinitestrain) then
      do k = 1,ngvec
!
!  Create vector scaled by inverse strain matrix
!
        gloc(1) = strainfull(1,1)*xrk(k) + &
                  strainfull(2,1)*yrk(k)
        gloc(2) = strainfull(1,2)*xrk(k) + &
                  strainfull(2,2)*yrk(k)
!
        do i = 1,2
          do j = 1,3
            dgds(k,1,j) = dgds(k,1,j) + dIds(j,1,i)*gloc(i)
            dgds(k,2,j) = dgds(k,2,j) + dIds(j,2,i)*gloc(i)
          enddo
        enddo
      enddo
    else
      do k = 1,ngvec
        dgds(k,1,1) = - xrk(k)
        dgds(k,1,3) = - 0.5_dp*yrk(k)
        dgds(k,2,2) = - yrk(k)
        dgds(k,2,3) = - 0.5_dp*xrk(k)
      enddo
    endif
    if (lgrad2) then
      if (lfinitestrain) then
        do k = 1,ngvec
!
!  Create vector scaled by inverse strain matrix
!
          gloc(1) = strainfull(1,1)*xrk(k) + &
                    strainfull(2,1)*yrk(k)
          gloc(2) = strainfull(1,2)*xrk(k) + &
                    strainfull(2,2)*yrk(k)
!
          do i = 1,2
            do ns1 = 1,3
              do ns2 = 1,3
                d2gds2(k,1,ns2,ns1) = d2gds2(k,1,ns2,ns1) + d2Ids2(ns2,ns1,1,i)*gloc(i)
                d2gds2(k,2,ns2,ns1) = d2gds2(k,2,ns2,ns1) + d2Ids2(ns2,ns1,2,i)*gloc(i)
              enddo
            enddo
          enddo
        enddo
      else
        do k = 1,ngvec
          d2gds2(k,1,1,1) = xrk(k)
          d2gds2(k,1,3,1) = 0.5_dp*yrk(k)
!
          d2gds2(k,1,2,3) = 0.5_dp*yrk(k)
          d2gds2(k,1,3,3) = 0.25_dp*xrk(k)
!
          d2gds2(k,2,2,2) = yrk(k)
          d2gds2(k,2,3,2) = 0.5_dp*xrk(k)
!
          d2gds2(k,2,1,3) = 0.5_dp*xrk(k)
          d2gds2(k,2,3,3) = 0.25_dp*yrk(k)
        enddo
      endif
    endif
  elseif (ndim.eq.1) then
!
!  Initialise derivatives
!
    dgds(1:ngvec,1:3,1) = 0.0_dp
    if (lgrad2) then
      d2gds2(1:ngvec,1:3,1,1) = 0.0_dp
    endif
    if (lfinitestrain) then
!
!  Compute inverse matrix
!
      straininv(1,1) = strainfull(2,2)
!
!  Create overall strain derivatives of inverse matrix
!
      dIds(1,1,1) = straindet*straindet*strainddetds(1)*straininv(1,1)
    endif
!
    if (lfinitestrain) then
      do k = 1,ngvec
        gloc(1) = xrk(k)*strainfull(1,1)
!
!  Length derivative
!
        dgds(k,1,1) = dgds(k,1,1) + dIds(1,1,1)*gloc(1)
      enddo
    else
      do k = 1,ngvec
        dgds(k,1,1) = - xrk(k)
      enddo
    endif
    if (lgrad2) then
      do k = 1,ngvec
        d2gds2(k,1,1,1) = - dgds(k,1,1)
      enddo
    endif
  endif
#ifdef TRACE
  call trace_out('gxyzterms')
#endif
!
  return
  end

end module m_strain
