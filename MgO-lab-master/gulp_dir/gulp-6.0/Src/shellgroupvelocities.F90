  subroutine shellgroupvelocities(mcv,msv,derv2,dervi,maxd2)
!
!  Corrects the group velocities for shell contributions
!
!   8/14 Created from groupvelocities
!   2/18 Trace added
!   6/20 Corrections to second term
!
!  On entry :
!
!  mcv         = no. of core coordinates
!  msv         = no. of shell coordinates
!  derv2       = inverse shell-shell matrix multiplied by shell-core part - real part
!  dervi       = inverse shell-shell matrix multiplied by shell-core part - imaginary part
!  maxd2       = left-hand dimension of eigr
!
!  On exit :
!
!  Group velocities are reduced to core-core only part
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
!  Julian Gale, CIC, Curtin University, June 2020
!
  use derivatives,   only : derv2dk
  use element
  use frequencies
  use iochannels
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4),      intent(in)    :: maxd2
  integer(i4),      intent(in)    :: mcv
  integer(i4),      intent(in)    :: msv
  real(dp),         intent(in)    :: derv2(maxd2,mcv)
  real(dp),         intent(in)    :: dervi(maxd2,mcv)
!
!  Local variables
!
  integer(i4)                     :: i
  integer(i4)                     :: j
  integer(i4)                     :: k
  integer(i4)                     :: l
  complex(dpc)                    :: cjl
  complex(dpc)                    :: cki
  complex(dpc)                    :: ckj
  complex(dpc)                    :: dk(3)
#ifdef TRACE
  call trace_in('shellgroupvelocities')
#endif
!
!  Add shell corrections to derv2dk
!
  do i = 1,mcv
    do j = 1,mcv
      do k = 1,msv
!
!  Make complex terms from separate real and imaginary parts
!
        cki = dcmplx(derv2(mcv+k,i),dervi(mcv+k,i))
        do l = 1,msv
!
!  Make complex terms from separate real and imaginary parts
!
          cjl = dcmplx(derv2(mcv+l,j),-dervi(mcv+l,j))
!
!  First term  - derivatives of inverse of shell-shell matrix
!
          derv2dk(1:3,j,i) = derv2dk(1:3,j,i) + cjl*derv2dk(1:3,mcv+l,mcv+k)*cki
        enddo
!
!  Second terms - derivatives of core-shell matrix
!
        ckj = dcmplx(derv2(mcv+k,j),dervi(mcv+k,j))
        dk(1:3) = derv2dk(1:3,j,mcv+k)*cki + conjg(ckj)*derv2dk(1:3,mcv+k,i)
        derv2dk(1:3,j,i) = derv2dk(1:3,j,i) - dk(1:3)
      enddo
    enddo
  enddo
#ifdef TRACE
  call trace_out('shellgroupvelocities')
#endif
!
  return
  end
