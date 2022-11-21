  function epsilon(x)
  use datatypes
  real(dp) :: epsilon
  real(dp) :: x
!
!     estimate unit roundoff in quantities of size x.
!
  real(dp) :: a,b,c,eps
!
!**** warning:  rumor has it that some machines (including the IBM PC?)
!****   have compilers that defeat the following scheme.
!     this program should function properly on all systems
!     satisfying the following two assumptions,
!        1.  the base used in representing floating point
!            numbers is not a power of three.
!        2.  the quantity  a  in statement 10 is represented to 
!            the accuracy used in floating point variables
!            that are stored in memory.
!     the statement number 10 and the go to 10 are intended to
!     force optimizing compilers to generate code satisfying 
!     assumption 2.
!     under these assumptions, it should be true that,
!            a  is not exactly equal to four-thirds,
!            b  has a zero for its last bit or digit,
!            c  is not exactly equal to one,
!            eps  measures the separation of 1.0 from
!                 the next larger floating point number.
!     the developers of eispack would appreciate being informed
!     about any systems where these assumptions do not hold.
!
!     this version dated 4/6/83.
!
  a = 4.0_dp/3.0_dp
10 b = a - 1.0_dp
  c = b + b + b
  eps = dabs(c-1.0_dp)
  if (eps .eq. 0.0_dp) go to 10
  epsilon = eps*dabs(x)
  return
  end
