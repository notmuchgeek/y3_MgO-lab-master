!*****************************************************
!  Serial matrix inversion using standalone routine  *
!*****************************************************
  subroutine matrix_inversion(a,ia,n,wrk,ifail)
!
!  Matrix inverter
!
!  On entry :
!
!  a     = matrix to be inverted
!  ia    = lower dimension of a
!  n     = actual size of matrix to be inverted
!  wrk   = workspace array of length 2*n
!
!  On exit :
!
!  a     = inverse matrix
!  ifail = 0, if OK
!
!   3/14 Renamed from matinv for benefit of ChemShell
!   2/18 Trace added
!
  implicit none
!
!  Passed variables
!
  integer*4    :: ia
  integer*4    :: ifail
  integer*4    :: n
  real*8       :: a(ia,*)
  real*8       :: wrk(*)
!
!  Local variables
!
  integer*4    :: i
  integer*4    :: iarg1
  integer*4    :: iarg2
  integer*4    :: i1
  integer*4    :: ins
  integer*4    :: j
  integer*4    :: jj1
  integer*4    :: k
  integer*4    :: kk
  integer*4    :: l
  integer*4    :: l1
  integer*4    :: m
  integer*4    :: nupper
  real*8       :: acc
  real*8       :: aloc
  real*8       :: t
!
  acc = 1.0d-8
  ifail = 0
  do j = 1,n
    if (j.ne.1) then
      iarg1 = n-j+1
      iarg2 = j-1
      call mxm1(a(j,1),ia,a(1,j),ia,wrk,ia,iarg1,iarg2)
      nupper = n - j + 1
      do l = 1,nupper
        a(j+l-1,j) = a(j+l-1,j) + wrk(l)
      enddo
    endif
    t = abs(a(j,j))
    k = j
    if (j.ne.n) then
      do i = j+1,n
        if (abs(a(i,j)).gt.t) then
          t = abs(a(i,j))
          k = i
        endif
      enddo
    endif
    wrk(j+n) = dble(k)
    if (t.le.acc) then
      ifail = 1
      return
    endif
    if (k.ne.j) then
      do m = 1,n
        t = a(j,m)
        a(j,m) = a(k,m)
        a(k,m) = t
      enddo
    endif
    a(j,j) = 1.0d0/a(j,j)
    if (j.ne.n) then
      if (j.ne.1) then
        iarg1 = n - j
        iarg2 = j - 1
        call mxm2(a(1,j+1),ia,a(j,1),ia,wrk,iarg1,iarg2)
        nupper = n - j
        do l1 = 1,nupper
          a(j,j+l1) = a(j,j+l1) + wrk(l1)
        enddo
      endif
      t = - a(j,j)
      nupper = n - (j+1)
      do i1 = j+1,n
        a(j,i1) = t*a(j,i1)
      enddo
    endif
  enddo
!
!  Use cminv method to solve for a**-1
!
  do k = 2,n
    nupper = k - 1
    do m = 1,nupper
      wrk(m) = 0.0d0
    enddo
    do j = 1,k-1
      aloc = a(k,j)
      do m = 1,j
        wrk(m) = wrk(m)-aloc*a(j,m)
      enddo
    enddo
    aloc = a(k,k)
    nupper = k - 1
    do m = 1,nupper
      a(k,m) = wrk(m)*aloc
    enddo
  enddo
!
!  Now back substitution
!
  k = n
  do kk = 2,n
    k = k - 1
    jj1 = kk - 1
    call mxm2(a(k+1,1),ia,a(k,k+1),ia,wrk,n,jj1)
    do l = 1,k
      wrk(l) = wrk(l)+a(k,l)
    enddo
    do j = 1,n
      a(k,j) = wrk(j)
    enddo
  enddo
!
!  Multiply solution by inverse of permutation matrix
!
  k = n + 1
  do i = 1,n
    k = k-1
    ins = int(wrk(k+n))
    if (ins.ne.k) then
      do j = 1,n
        wrk(j) = a(j,ins)
        a(j,ins) = a(j,k)
        a(j,k) = wrk(j)
      enddo
    endif
  enddo
  return
  end
!******************************
!  Ancillary matrix routines  *
!******************************
  subroutine mxm1(rarr,ir,sarr,jr,tarr,kr,id1,id2)
!
!  Matrix multiplier
!
  implicit none
!
!  Passed variables
!
  integer*4    :: id1
  integer*4    :: id2
  integer*4    :: ir
  integer*4    :: jr
  integer*4    :: kr
  real*8       :: rarr(*)
  real*8       :: sarr(*)
  real*8       :: tarr(*)
!
!  Local variables
!
  integer*4    :: i
  integer*4    :: ia
  integer*4    :: j
  real*8       :: sum
!
  do i = 1,id1
    ia = i - ir
    sum = 0.0d0
    do j = 1,id2
      ia = ia + ir
      sum = sum + rarr(ia)*sarr(j)
    enddo
    tarr(i) = sum
  enddo
  return
  end
!
  subroutine mxm2(rarr,ic,sarr,jc,tarr,id1,id2)
!
!  Matrix multiplier
!
  implicit none
!
!  Passed variables
!
  integer*4    :: id1
  integer*4    :: id2
  integer*4    :: ic
  integer*4    :: jc
  real*8       :: rarr(*)
  real*8       :: sarr(*)
  real*8       :: tarr(*)
!
!  Local variables
!
  integer*4    :: i
  integer*4    :: ia
  integer*4    :: ira
  integer*4    :: j
  integer*4    :: ja
  real*8       :: sum
!
  ira = 1 - ic
  do i = 1,id1
    ira = ira + ic
    ia = ira - 1
    ja = 1 - jc
    sum = 0.0d0
    do j = 1,id2
      ja = ja + jc
      sum = sum + rarr(ia+j)*sarr(ja)
    enddo
    tarr(i) = sum
  enddo
  return
  end
