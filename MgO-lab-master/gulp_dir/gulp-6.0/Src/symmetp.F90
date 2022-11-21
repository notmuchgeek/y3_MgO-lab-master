  subroutine symmetp
!***********************************************************************
!
!     Adapted from routine in CRYSTAL88 : QCPE
!
!     Generation of the reciprocal space group for a 
!     three-dimensional lattice. Adapted from symmet.
!
!  12/14 Setting of hmssg removed as this messes up the real space
!        symmetry and is not needed
!
!***********************************************************************
  use control
  use current
  use iochannels
  use numbers
  use parallel
  use symmetry
  use times
#ifdef TRACE
  use trace,       only : trace_in, trace_out
#endif
  implicit none
!
!  Local variables
!
  character(len=16)             :: gronam2(17)
  character(len=1)              :: hmssg2(16)
  integer(i4)                   :: i
  integer(i4)                   :: id
  integer(i4)                   :: iflags2
  integer(i4)                   :: ifso2
  integer(i4)                   :: ii
  integer(i4)                   :: il
  integer(i4)                   :: isc(3,3)
  integer(i4)                   :: in(3,4)
  integer(i4)                   :: indgro2(17)
  integer(i4),             save :: iq(3,4,5)
  integer(i4)                   :: is(3)
  integer(i4),             save :: it(3,48)
  integer(i4)                   :: j
  integer(i4)                   :: j1
  integer(i4)                   :: j2
  integer(i4)                   :: j3
  integer(i4)                   :: j4
  integer(i4)                   :: jl
  integer(i4)                   :: js(3)
  integer(i4)                   :: k
  integer(i4)                   :: k1
  integer(i4)                   :: k2
  integer(i4)                   :: k3
  integer(i4)                   :: k4
  integer(i4)                   :: kk
  integer(i4)                   :: kl
  integer(i4)                   :: l
  integer(i4)                   :: lgt(48)
  integer(i4)                   :: lt
  integer(i4)                   :: m
  integer(i4)                   :: mcc
  integer(i4)                   :: mv
  integer(i4)                   :: n
  integer(i4),             save :: ncfold = 0
  integer(i4)                   :: nd
  integer(i4)                   :: ne
  integer(i4)                   :: net(3)
  integer(i4)                   :: ngt(48,48)
  integer(i4)                   :: nrt
  integer(i4)                   :: nspcg2
  real(dp)                      :: g_cpu_time
  real(dp)                      :: r24
  real(dp)                      :: time1
  real(dp)                      :: time2
!
  equivalence (j2,in(2,2)),(j3,in(2,3)),(k2,in(3,2)),(k3,in(3,3)), &
              (j1,in(2,1)),(k1,in(3,1)),(j4,in(2,4)),(k4,in(3,4))
!
  data (gronam2(i),i=1,17)/'P M -3          ', &
    'P N -3          ', 'F M -3          ', 'F D -3          ', &
    'I M -3          ', 'P A -3          ', 'I A -3          ', &
    'P M -3 M        ', 'P N -3 N        ', 'P M -3 N        ', &
    'P N -3 M        ', 'F M -3 M        ', 'F M -3 C        ', &
    'F D -3 M        ', 'F D -3 C        ', 'I M -3 M        ', &
    'I A -3 D        '/
  data (indgro2(i),i=1,17)/200,201,202,203,204,205,206,221,222,223,224,225,226,227,228,229,230/
#ifdef TRACE
  call trace_in('symmetp')
#endif
!
!  Copy centring matrices into w/w1 arrays
!
  ii = 0
  do i = 1,3
    do j = 1,3
      do k = 1,7
        ii = ii+1
        wp(k,j,i) = wdat(ii)
      enddo
    enddo
  enddo
  do i = 1,3
    do j = 1,3
      do k = 1,7
        ii = ii + 1
        w1p(k,j,i) = wdat(ii)
      enddo
    enddo
  enddo
!
  time1 = g_cpu_time()
!
!  Set local variables
!
  iflags2 = 0
  ifso2 = 0
  nspcg2 = ipatptr(ipatgrp(nspcg(ncf)))
!  nspcg2 = nspcgp(ncf)
  hmssg2 = ' '
  do j = 1,9
    hmssg2(j) = patter(ipatgrp(nspcg(ncf)))(j:j)
  enddo
!  if (iperm(ncf).gt.1) then
!    do j = 1,16
!      hmssg2(j) = gronam(nspcg2)(j:j)
!    enddo
!  else
!    do j = 1,16
!      hmssg2(j) = hmssg(j,ncf)
!    enddo
!  endif
  if (ncf.ne.ncfold) lalterp = .false.
  ncfold = ncf
!
  ii = 0
  do i = 1,5
    do j = 1,4
      do k = 1,3
        ii = ii + 1
        iq(k,j,i) = iqdd(ii)
      enddo
    enddo
  enddo
  n = 1
  k = 1
  nccsp = 2
  id = 0
  il = 0
  lt = 0
  nrt = 0
  mcc = 0
  ncsp = 0
  ngop = 0
  do i = 1,3
    js(i) = 0
    in(i,1) = 0
    in(i,2) = 0
    in(i,3) = 0
    in(i,4) = 0
    is(i) = 1
    net(i) = 1
    do j = 1,48
      it(i,j) = 0
    enddo
  enddo
  do j = 1,48
    lgt(j) = 0
  enddo
!
!  Division by 3 of W matrix for Rhombohedral cell
!
  do i = 1,3
    do j = 1,3
      wp(7,j,i) = third*wp(7,j,i)
    enddo
  enddo
!
!  Determine Bravais lattice and crystal family
!
  if (iflags2.eq.0) goto 577
  goto 270   
!
!  Space group number input case
!
577 if (nspcg2.eq.0) nspcg2 = 1
  if (nspcg2.gt.232.or.nspcg2.le.0) then
    call outerror('space group number outside range',0_i4)
    call stopnow('symmetp')
  endif
  do i = 1,16
    hmssg2(i) = gronam(nspcg2)(i:i)
  enddo
!
!  Determine space group number
!
270 if (nspcg2.le.1) then
    do i = 1,232
      do j = 1,16
        if (hmssg2(j).ne.gronam(i)(j:j)) goto 345
      enddo
      nspcg2 = i
      nspcgp(ncf) = i
      goto 347
345   continue
    enddo
!
!  Check alternative symbols
!
    do i = 1,17
      do j = 1,16
        if (hmssg2(j).ne.gronam2(i)(j:j)) goto 346
      enddo
      nspcg2 = indgro2(i)
      nspcgp(ncf) = nspcg2
      do j = 1,16
        hmssg2(j) = gronam(nspcg2)(j:j)
      enddo
      goto 347
346   continue
    enddo
!
    call settings(hmssg2,nspcg2)
    if (nspcg2.gt.0) then
      lalterp = .true.
      nspcgp(ncf) = nspcg2
      goto 347
    endif
    if (ioproc) then
      write(ioout,'(/,''  **** Invalid space group symbol supplied : '',16a1,'' ****'',/)')(hmssg2(i),i = 1,16)
    endif
    call stopnow('symmetp')
347 continue
  endif
!
!  Start processing symbol
!
  kk = 0
275 kk = kk + 1
  if (kk.gt.14) then
    if (ioproc) then
      write(ioout,'('' **** Too many blanks before space group symbol ****'')')
    endif
    call stopnow('symmetp')
  endif
  if (hmssg2(kk).eq.hbr(8)) goto 275
  do i = kk,16
    hmssg2(i-kk+1) = hmssg2(i)
  enddo
  do i = 18-kk,16
    hmssg2(i) = hbr(8)
  enddo
!
!  Remove excess blanks from symbol
!
  do i = 3,15
    if (hmssg2(i-1).eq.hbr(8).and.hmssg2(i).eq.hbr(8).and.hmssg2(i+1).ne.hbr(8)) then
      do j = i,15
        hmssg2(j) = hmssg2(j+1)
      enddo
    endif
  enddo
  do i = 1,7
    if (hmssg2(1).eq.hbr(i)) ncblp = i
  enddo
  if (ncblp.ne.7.and.ifhr(ncf).ne.0) ifhr(ncf) = 0
!
!  If B centred then the W1 matrix must be corrected
!
  if (ncblp.eq.3) then
    w1p(3,1,1) = 1.0_dp
    w1p(3,2,1) = 0.0_dp
    w1p(3,3,1) = 1.0_dp
    w1p(3,1,2) = 0.0_dp
    w1p(3,2,2) = 1.0_dp
    w1p(3,3,2) = 0.0_dp
    w1p(3,1,3) = -1.0_dp
    w1p(3,2,3) = 0.0_dp
    w1p(3,3,3) = 1.0_dp
  endif
!
!  Standardisation of the symbol
!
  do i = 3,16
    if (hmssg2(i).eq.hbr(8)) then
      k = 0
      n = n + 1
    else
      if (hmssg2(i).eq.hs(11).or.hmssg2(i).eq.hs(14)) nccsp = 5
      do j = 1,15
        if (hmssg2(i).eq.hs(j)) in(n,k) = j - 8
      enddo
    endif
    k = k + 1
  enddo
  if (in(1,1).eq.4.or.in(1,1).eq.-1.and.in(1,2).eq.4) nccsp = 4
  if (j1.eq.3) nccsp = 6
  if (j1.eq.0.and.in(1,1).eq.1.or.in(1,1).eq.-1.and.in(1,2).eq.1) nccsp = 1
  if (nccsp.ne.2.or.j1.ne.0) then
    if (nccsp.eq.2.and.in(1,1).ne.1.and.j1.ne.1.and.k1.ne.1) nccsp = 3
  else
    in(2,1) = in(1,1)
    in(1,1) = 0
    in(2,2) = in(1,2)
    in(1,2) = 0
    in(2,3) = in(1,3)
    in(1,3) = 0
    in(2,4) = in(1,4)
    in(1,4) = 0
    in(1,1) = 1
    k1 = 1
  endif
!
!  Determination of code no. ncpgp for the point group
!
  do i = 1,3
    if (in(i,1).ne.0) is(i) = in(i,1)
    if (is(i).eq.-1)  is(i) = -in(i,2)
    if (in(i,1).lt.-1) is(i) = 7
    if (in(i,1).le.-2) then
      in(i,4) = in(i,1)
      in(i,1) = 1
    endif
    if (in(i,2).eq.7) then
      in(i,4) = in(i,3)
      in(i,3) = 7
      in(i,2) = 0
    endif
  enddo
  if (k1.eq.2.and.(j1.eq.2.or.j1.eq.3)) nrt = 1
  m = 100*is(1) + 10*is(2) + is(3) + in(1,3) + j3 + k3
  do i = 1,45
    if (m.eq.iss(i)) ncpgp = i
  enddo
!
!  Trap case if ncpgp is zero
!
  if (ncpgp.eq.0) then
    call outerror('Problem with space group symbol - ncpgp is zero',0_i4)
    call stopnow('symmetp')
  endif
!
!  Fix code no. ncsp for centrosymmetry
!
  if (ncpgp.le.31) then
    ncsp = 1
    nd = ncpgp
  else
    nd = iss(ncpgp+14)
  endif
!
!  Rotation parts
!
  do i = 1,216
    if (lf(i).eq.1) id = id + 1
    if (id-nd.gt.0) then
      goto 420
    elseif (id-nd.eq.0) then
      m = abs(lf(i))
      if (nccsp.eq.5) m = m + 24
      ngop = ngop + 1
      do k = 1,3
        do l = 1,3
          ropp(l,k,ngop) = ige(l,k,m)*sign(1_i4,lf(i))
        enddo
      enddo
    endif
  enddo
!
!  Multiplication table
!
420 continue
  do i = 1,ngop
    do j = 1,3
      do k = 1,3
        if (ncsp.eq.0) ropp(k,j,ngop+i) = - ropp(k,j,i)
      enddo
    enddo
    do 480 j = 1,ngop
      do k = 1,3
        do l = 1,3
          isc(l,k) = 0
        enddo
      enddo
      do k = 1,3
        do l = 1,3
          do m = 1,3
            isc(k,l) = isc(k,l) + ropp(k,m,i)*ropp(m,l,j)
          enddo
        enddo
      enddo
      do 470 k = 1,ngop
        do l = 1,3
          do m = 1,3
            if (isc(m,l).ne.ropp(m,l,k)) goto 470
          enddo
        enddo
        ngt(i,j) = k
        if (ncsp.eq.0) then
          ngt(ngop+i,j) = k + ngop
          ngt(i,ngop+j) = k + ngop
          ngt(ngop+i,ngop+j) = k
        endif
        goto 480
470   continue
480 continue
  enddo
  if (ncsp.eq.0) ngop = ngop + ngop
!
!  Select translation parts for generators
!
  if (ncpgp.eq.30.or.ncpgp.eq.31) il = 1
  if (nccsp.eq.3) mcc = 400*ncblp + in(1,4)*49 + j4*7 + k4 + 399
  if (nccsp.ge.4) mcc = 1000*nd + 100*ncblp + 25*((k4+7)/2) + 5*in(1,4) + j4 + 42 + 4*in(1,2) + j2*2
  if (nccsp.eq.5) mcc = 3000 + 100*ncpgp + in(1,2) + j1 + k1
  if (ncpgp.eq.11) mcc = in(1,2)*4 + j2*2 + k2
  do i = 1,3
    l = in(i,4) + 8
    if (l.lt.7) then
      il = il + 1
      net(il) = mst1(i,l)
      if (nccsp.gt.3) net(il) = mst2(i,l)
    endif
  enddo
  if (ncpgp.eq.30.and.k4.eq.-2) net(il) = 13
  if (il.lt.2) then
    do i = 1,3
      if (in(i,1).eq.2) then
        il = il + 1
        net(il) = 1 + i*in(i,2)
        if (nccsp.gt.3.and.in(i,1).ne.-1) net(il) = 1 + in(i,2)
        if (nccsp.eq.6) net(il) = 1 + 5*in(i,2)
      endif
    enddo
    if (in(1,1).ge.3.and.in(1,1).le.6.and.in(1,2).ne.0) then
      il = il + 1
      net(il) = 13 + (10*in(1,2))/in(1,1)-(2*in(1,2))/in(1,1)
    endif
  endif
  n = net(3)
  if (ncpgp.eq.31) n = 1 + 7*(in(1,2)/2)+in(1,2)/3+in(1,2)*(10*in(1,2)*in(1,2)-42*in(1,2)+43+mod(ncblp,5_i4))
  do i = 1,3
    id = lgn(i,ncpgp)
    if (id.eq.0) goto 530
    m = net(i)
    do j = 1,3
      it(j,id) = ket(j,m) - nrt*(i/2)*(2/i)*(1-2*(nccsp/6))*ket(j,n)
      if (nrt.eq.0.and.nccsp.lt.4.and.ncblp.ne.5) js(j) = ket(j,m)/2 + js(j)
    enddo
  enddo
!
!  Complete translation parts
!
530 do i = 1,3
    if (lgn(i,ncpgp).ne.0) then
      lt = lt + 1
      lgt(lt) = lgn(i,ncpgp)
    endif
  enddo
550 if (lt.eq.ngop) goto 588
  n = 0
  do i = 1,lt
    do j = 1,lt
      il = lgt(i)
      jl = lgt(j)
      kl = ngt(il,jl)
      do m = 1,ngop
        if (lgt(m).eq.kl) goto 590
      enddo
      do k = 1,3
        it(k,kl) = it(k,kl) + ropp(k,1,il)*it(1,jl) + ropp(k,2,il)*it(2,jl) + ropp(k,3,il)*it(3,jl)
        it(k,kl) = it(k,kl) + it(k,il) + 48
        it(k,kl) = mod(it(k,kl),24_i4)
      enddo
      n = n + 1
      lgt(lt+n) = kl
590   continue
    enddo
  enddo
  lt = lt + n
  goto 550
!****************************
!  Selection of the origin  *
!****************************
588 if (abs(ishorg(1,nspcg2))+abs(ishorg(2,nspcg2))+abs(ishorg(3,nspcg2)).gt.0) then
    if (ifso2.eq.1) goto 620
  else
    if (ifso2.eq.0) goto 620
  endif
  if (ifso2.gt.1) then
    js(1) = js(1) + ivso(1,ncf)
    js(2) = js(2) + ivso(2,ncf)
    js(3) = js(3) + ivso(3,ncf)
  elseif (nspcg2.gt.0) then
    js(1) = js(1) + ishorg(1,nspcg2)
    js(2) = js(2) + ishorg(2,nspcg2)
    js(3) = js(3) + ishorg(3,nspcg2)
  endif
620 n = 0
  do i = 1,78
    if (mcc.eq.lsh(1,i)) n = lsh(2,i)
  enddo
  if (n.ne.0) then
    js(1) = js(1) + ket(1,n)
    js(2) = js(2) + ket(2,n)
    js(3) = js(3) + ket(3,n)
  endif
  do i = 1,ngop
    do j = 1,3
      it(j,i) = it(j,i) + ropp(j,1,i)*js(1) + ropp(j,2,i)*js(2) + ropp(j,3,i)*js(3) - js(j) + 48
      it(j,i) = mod(it(j,i),24_i4)
    enddo
    if (ncblp.ne.1.and.ncblp.ne.7) then
      nd = 72
      do l = 1,4
        ne = 0
        do k = 1,3
          is(k) = it(k,i) + iq(k,l,ncblp-1)
          ne = ne + mod(is(k),24_i4)
        enddo
        if (nd.gt.ne) then
          nd = ne
          it(1,i) = mod(is(1),24_i4)
          it(2,i) = mod(is(2),24_i4)
          it(3,i) = mod(is(3),24_i4)
        endif
      enddo
    endif
  enddo
!
!  Divide it by 24
!
  r24 = 1.0_dp/24.0_dp
  do mv = 1,ngop
    vitp(1,mv) = it(1,mv)*r24
    vitp(2,mv) = it(2,mv)*r24
    vitp(3,mv) = it(3,mv)*r24
  enddo
!
!  Debugging print statements!
!
  if (index(keyword,'opera').ne.0.and.ioproc) then
    write(ioout,'(/,''  Symmetry operators :'')')
    do i = 1,ngop
      write(ioout,'(/,''  Operator no.  =  '',i2,/)') i
      write(ioout,'(3f5.1,f6.2)') ropp(1,1,i),ropp(2,1,i),ropp(3,1,i),vitp(1,i)
      write(ioout,'(3f5.1,f6.2)') ropp(1,2,i),ropp(2,2,i),ropp(3,2,i),vitp(2,i)
      write(ioout,'(3f5.1,f6.2)') ropp(1,3,i),ropp(2,3,i),ropp(3,3,i),vitp(3,i)
    enddo
  endif
  time2 = g_cpu_time()
  tsym = tsym + time2 - time1
#ifdef TRACE
  call trace_out('symmetp')
#endif
!
  return
  end
