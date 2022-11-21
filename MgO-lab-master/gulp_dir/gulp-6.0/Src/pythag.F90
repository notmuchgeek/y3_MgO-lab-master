  function e_pythag(a,b)
!
!  finds sqrt(a**2+b**2) without overflow or destructive underflow
!
  use datatypes
  real(dp)  :: e_pythag
  real(dp)  :: a
  real(dp)  :: b
!
  real(dp)  :: p
  real(dp)  :: r
  real(dp)  :: s
  real(dp)  :: t
  real(dp)  :: u
!
  p = dmax1(dabs(a),dabs(b))
  if (p .eq. 0.0_dp) go to 20
  r = (dmin1(dabs(a),dabs(b))/p)**2
10 continue
     t = 4.0_dp + r
     if (t .eq. 4.0_dp) go to 20
     s = r/t
     u = 1.0_dp + 2.0_dp*s
     p = u*p
     r = (s/u)**2 * r
  go to 10
20 e_pythag = p
  return
  end
