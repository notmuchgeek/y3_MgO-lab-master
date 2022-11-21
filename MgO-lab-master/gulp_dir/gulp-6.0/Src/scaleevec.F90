  subroutine scaleevec(mcv,mcvloc,eigr,maxd2,scale)
!
!  Scales the eigenvectors by a mode dependant constant
!
!   9/97 Created 
!   4/17 Changes made for parallel version
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
  use datatypes
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4)        :: maxd2
  integer(i4)        :: mcv
  integer(i4)        :: mcvloc
  real(dp)           :: eigr(maxd2,*)
  real(dp)           :: scale(*)
!
!  Local variables
!
  integer(i4)        :: i
  integer(i4)        :: j
  real(dp)           :: rscale
#ifdef TRACE
  call trace_in('scaleevec')
#endif
!
!  Scale by mode
!
  do i = 1,mcvloc
    rscale = scale(i)
    do j = 1,mcv
      eigr(j,i) = rscale*eigr(j,i)
    enddo
  enddo
#ifdef TRACE
  call trace_out('scaleevec')
#endif
!
  return
  end
!
  subroutine scaleevecc(mcv,mcvloc,eigc,maxd2,scale)
!
!  Scales the complex eigenvectors by a mode dependant constant
!
!   4/17 Created from scaleevec
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
  use datatypes
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4)        :: maxd2
  integer(i4)        :: mcv
  integer(i4)        :: mcvloc
  complex(dpc)       :: eigc(maxd2,*)
  real(dp)           :: scale(*)
!
!  Local variables
!
  integer(i4)        :: i
  integer(i4)        :: j
  complex(dpc)       :: cscale
#ifdef TRACE
  call trace_in('scaleevecc')
#endif
!
!  Scale by mode
!
  do i = 1,mcvloc
    cscale = dcmplx(scale(i),0.0_dp)
    do j = 1,mcv
      eigc(j,i) = cscale*eigc(j,i)
    enddo
  enddo
#ifdef TRACE
  call trace_out('scaleevecc')
#endif
!
  return
  end
