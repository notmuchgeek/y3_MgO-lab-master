  subroutine setdcosmoamat1p(ldqneeded,dcosmoA,dsegweight,nvec,nmid,nearsas,nearsasptr, &
                             nearsasrptr,lnearsas,lanynonunit,xvec,yvec,zvec,dcosmoA2)
!
!  Subroutine calculates the first derivatives of the COSMO A matrix in parallel
!
!   4/17 Created from setdcosmoamat
!   8/19 Correction to one call to cosmoamatdadd1p
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
  use cosmic
  use g_constants
  use control
  use current
  use derivatives
  use iochannels
  use optimisation
  use parallel
  use reallocate
  use wolfcosmo,      only : lPureCoulomb0D
  implicit none
!
!  Passed variables
!
  integer(i4),       intent(in)                 :: nearsas
  integer(i4),       intent(in)                 :: nearsasptr(*)
  integer(i4),       intent(in)                 :: nearsasrptr(*)
  integer(i4),       intent(in)                 :: nmid
  integer(i4),       intent(in)                 :: nvec
  logical,           intent(in)                 :: lanynonunit(*)
  logical,           intent(in)                 :: ldqneeded
  logical,           intent(in)                 :: lnearsas(*)
  real(dp),          intent(inout)              :: dcosmoA(3,nptsonnode,*)
  real(dp),          intent(inout)              :: dcosmoA2(3,nptsonnode,*)
  real(dp),          intent(inout)              :: dsegweight(3,maxnearseg,*)
  real(dp),          intent(in)                 :: xvec(*)
  real(dp),          intent(in)                 :: yvec(*)
  real(dp),          intent(in)                 :: zvec(*)
!
!  Local variables
!
  integer(i4)                                   :: i
  integer(i4)                                   :: ia
  integer(i4)                                   :: ierror
  integer(i4)                                   :: ii
  integer(i4)                                   :: in
  integer(i4)                                   :: ipts
  integer(i4)                                   :: iptsloc
  integer(i4)                                   :: j
  integer(i4)                                   :: j1
  integer(i4)                                   :: j2
  integer(i4)                                   :: ja
  integer(i4)                                   :: jn
  integer(i4)                                   :: jpts
  integer(i4)                                   :: k
  integer(i4)                                   :: l
  integer(i4)                                   :: maxnpwt2
  integer(i4)                                   :: maxPperS
  integer(i4)                                   :: n
  integer(i4)                                   :: nn
  integer(i4)                                   :: nari
  integer(i4)                                   :: narj
  integer(i4)                                   :: nearsas2
  integer(i4)                                   :: npwt1
  integer(i4)                                   :: npwt2
  integer(i4)                                   :: nregioni
  integer(i4)                                   :: nregionj
  integer(i4)                                   :: nregionn
  integer(i4)                                   :: nsetfi
  integer(i4)                                   :: nsetfj
  integer(i4)                                   :: status
  logical,      dimension(:), allocatable, save :: ladrv
  logical                                       :: lopi
  logical                                       :: lopj
  logical                                       :: lperiodic
  logical                                       :: lrweight
  real(dp)                                      :: aij
  real(dp)                                      :: aij2
  real(dp)                                      :: d1rdist
  real(dp)                                      :: ddist
  real(dp)                                      :: dist
  real(dp)                                      :: dists
  real(dp)                                      :: dists2
  real(dp)                                      :: dqme(3)
  real(dp)                                      :: d2qme(6)
  real(dp)                                      :: d2wdr2
  real(dp)                                      :: dwdr
  real(dp)                                      :: dwtm1
  real(dp)                                      :: dwtm1swiswj
  real(dp)                                      :: dwt
  real(dp)                                      :: f0
  real(dp)                                      :: f1(3)
  real(dp)                                      :: f2(6)
  real(dp)                                      :: fact
  real(dp)                                      :: fdiag
  real(dp)                                      :: ff
  real(dp)                                      :: ffa
  real(dp)                                      :: pw1
  real(dp)                                      :: pw2
  real(dp)                                      :: qme
  real(dp)                                      :: qmeself
  real(dp)                                      :: qsi
  real(dp)                                      :: qsi2
  real(dp)                                      :: qsij
  real(dp)                                      :: qsijdwtm1swiswj
  real(dp)                                      :: qsipj
  real(dp)                                      :: qsipjdwtm1swiswj
  real(dp)                                      :: qsj
  real(dp)                                      :: ri
  real(dp)                                      :: rj
  real(dp)                                      :: rdist
  real(dp)                                      :: rdists
  real(dp)                                      :: rfact
  real(dp)                                      :: rnari
  real(dp)                                      :: rtrm
  real(dp)                                      :: rtrm2
  real(dp)                                      :: sumAinvipts
  real(dp)                                      :: sumAinvjpts
  real(dp)                                      :: swi
  real(dp)                                      :: swj
  real(dp)                                      :: x1, x2, x3
  real(dp)                                      :: xa, ya, za
  real(dp)                                      :: xc, yc, zc
  real(dp)                                      :: xi, yi, zi
  real(dp)                                      :: xij, yij, zij
  real(dp)                                      :: xj, yj, zj
  real(dp)                                      :: xji, yji, zji
  real(dp)                                      :: y1, y2, y3
  real(dp)                                      :: w1
  real(dp)                                      :: w2
  real(dp), dimension(:,:),   pointer,     save :: adrvi => null()
  real(dp), dimension(:,:),   pointer,     save :: adrvj => null()
  real(dp), dimension(:,:,:), pointer,     save :: a2drv => null()
  real(dp), dimension(:,:,:), pointer,     save :: dwti => null()
  real(dp), dimension(:,:,:), pointer,     save :: dwtj => null()
  real(dp), dimension(:,:,:,:), pointer,   save :: d2wti => null()
  real(dp), dimension(:,:,:,:), pointer,   save :: d2wtj => null()
  real(dp), dimension(:),     allocatable, save :: totalwt
!
!  Allocate local memory
!
  allocate(totalwt(npts),stat=status)
  if (status/=0) call outofmemory('setdcosmoamat','totalwt')
!
!  Set local constants
!
  fdiag = 1.05_dp*sqrt(nppa+0.0_dp)
  fact  = 0.5_dp*autoev*autoangs*cosmofneps
  lperiodic = (ndim.gt.0)
  nearsas2 = nearsas*(nearsas + 1)/2
!
!  Calculate self term
!
  if (lperiodic.and.lsegsmooth) then
    call qmatrixelementc(0.0_dp,0.0_dp,0.0_dp,0.0_dp,.true.,.false.,.false.,qmeself,dqme,d2qme)
  endif
!
!  Find maximum number of points per segment, maxPperS
!
  maxPperS = 0
  do ipts = 1,npts
    nari = nar(ipts)
    maxPperS = max(maxPperS,nari)
  enddo
!
!  Allocate arrays to store derivatives between i and other atoms
!
!  dwti   = derivative of weight w.r.t. position of i
!  dwtj   = derivative of weight w.r.t. position of j
!  d2wti  = second derivative of weight w.r.t. position of i
!  d2wtj  = second derivative of weight w.r.t. position of j
!  adrvi  = derivative of matrix element Aij w.r.t. i
!  adrvj  = derivative of matrix element Aij w.r.t. j
!  a2drv  = second derivative of matrix element Aij 
!
  call realloc(dwti,3_i4,maxnpwt,maxPperS,ierror)
  if (ierror.ne.0) call outofmemory('cosmoderv','dwti')
  call realloc(dwtj,3_i4,maxnpwt,maxPperS,ierror)
  if (ierror.ne.0) call outofmemory('cosmoderv','dwtj')
  call realloc(adrvi,3_i4,nearsas,ierror)
  if (ierror.ne.0) call outofmemory('cosmoderv','adrvi')
  call realloc(adrvj,3_i4,nearsas,ierror)
  if (ierror.ne.0) call outofmemory('cosmoderv','adrvj')
  allocate(ladrv(nearsas),stat=status)
  if (status/=0) call outofmemory('setdcosmoamat','ladrv')
  maxnpwt2 = 1
  call realloc(a2drv,1_i4,1_i4,1_i4,ierror)
  if (ierror.ne.0) call outofmemory('cosmoderv','a2drv')
  call realloc(d2wti,3_i4,3_i4,maxnpwt2,maxPperS,ierror)
  if (ierror.ne.0) call outofmemory('cosmoderv','d2wti')
  call realloc(d2wtj,3_i4,3_i4,maxnpwt2,maxPperS,ierror)
  if (ierror.ne.0) call outofmemory('cosmoderv','d2wtj')
!
!  Calculate the sum of the weighting factors for each segment
!
  do ipts = 1,npts
    i = cosmoatomptr(ipts)
    nari = nar(ipts)
    nsetfi = nsetf(ipts)
    rnari = 0.0_dp
    do k = 1,nari
      j1 = nset(k+nsetfi-1)
      w1 = cosmowt(j1,i)*cosmopwt(k+nsetfi-1)
      rnari = rnari + w1
    enddo
    totalwt(ipts) = 1.0_dp/rnari
  enddo
!
!  First initialisation of variables
!
  do n = 1,nearsas
    adrvi(1:3,n) = 0.0_dp
    adrvj(1:3,n) = 0.0_dp
  enddo
  ladrv(1:nearsas) = .false.
!********************************
!  Derivatives of the A matrix  *
!********************************
  do iptsloc = 1,nptsonnode
    ipts = node2pts(iptsloc)
    i = cosmoatomptr(ipts)
    ia = nrelf2a(i)
    in = nearsasrptr(i)
    npwt1 = npwt(ipts)
    nregioni = nregionno(nsft+ia)
    lopi = (.not.lfreeze.or.lopf(i).or.lnearsas(i))
    ri = atsrad(i) 
    swi = segweight(ipts)
    nari = nar(ipts)
    nsetfi = nsetf(ipts)
    qsi = qonsas(ipts)
    qsi2 = qsi*qsi*fact
    sumAinvipts = 0.0_dp
    qsipj = 0.0_dp
!
!  Get derivatives of weighting factors
!
    call dtotcosmowt(ipts,totalwt(ipts),dwti,d2wti,.false.)
    if (lanynonunit(i)) then
!
!  Self-term - derivatives come from the weighting factor
!
      do k = 1,nari
        j1 = nset(k+nsetfi-1)
!
!  Calculate derivatives of all weighting factors
!
        w1 = cosmowt(j1,i)
        pw1 = cosmopwt(k+nsetfi-1)
        w1 = w1*pw1*totalwt(ipts)
        call cosmoaiiderv(ipts,in,k,k,nearsasrptr,npwt1,maxnpwt2,fdiag,w1,w1,dwti,d2wti,adrvi,a2drv,ladrv,.false.)
!
        x1 = sphere2(1,j1)
        x2 = sphere2(2,j1)
        x3 = sphere2(3,j1)
        do l = 1,k-1
          j2 = nset(l+nsetfi-1)
          w2 = cosmowt(j2,i)
          pw2 = cosmopwt(l+nsetfi-1)
          w2 = w2*pw2*totalwt(ipts)
          dist = (x1-sphere2(1,j2))**2 + (x2-sphere2(2,j2))**2 + (x3-sphere2(3,j2))**2
          dist = 2.0_dp/sqrt(dist)
          call cosmoaiiderv(ipts,in,k,l,nearsasrptr,npwt1,maxnpwt2,dist,w1,w2,dwti,d2wti,adrvi,a2drv,ladrv,.false.)
        enddo
      enddo
      rfact = 1.0_dp/ri
      do n = 1,nearsas
        adrvi(1:3,n) = rfact*adrvi(1:3,n)
      enddo
!
!  Add derivatives connected with self-term on to main arrays
!
      do nn = 1,nearsas
        if (ladrv(nn)) then
          n = nearsasptr(nn)
          if (n .ne. i) then
            xdrv(i) = xdrv(i) - qsi2*adrvi(1,nn)
            ydrv(i) = ydrv(i) - qsi2*adrvi(2,nn)
            zdrv(i) = zdrv(i) - qsi2*adrvi(3,nn)
            xdrv(n) = xdrv(n) + qsi2*adrvi(1,nn)
            ydrv(n) = ydrv(n) + qsi2*adrvi(2,nn)
            zdrv(n) = zdrv(n) + qsi2*adrvi(3,nn)
            nregionn = nregionno(nsft+nrelf2a(n))
            if (nregioni.ne.nregionn) then
              xregdrv(nregioni) = xregdrv(nregioni) - qsi2*adrvi(1,nn)
              yregdrv(nregioni) = yregdrv(nregioni) - qsi2*adrvi(2,nn)
              zregdrv(nregioni) = zregdrv(nregioni) - qsi2*adrvi(3,nn)
              xregdrv(nregionn) = xregdrv(nregionn) + qsi2*adrvi(1,nn)
              yregdrv(nregionn) = yregdrv(nregionn) + qsi2*adrvi(2,nn)
              zregdrv(nregionn) = zregdrv(nregionn) + qsi2*adrvi(3,nn)
            endif
          endif
          if (ldqneeded) then
!  
!  Save first derivatives of A matrix term
!     
            dcosmoA(1,iptsloc,nn) = dcosmoA(1,iptsloc,nn) - qsi*adrvi(1,nn)
            dcosmoA(2,iptsloc,nn) = dcosmoA(2,iptsloc,nn) - qsi*adrvi(2,nn)
            dcosmoA(3,iptsloc,nn) = dcosmoA(3,iptsloc,nn) - qsi*adrvi(3,nn)
          endif
        endif
      enddo
!
!  Rezero derivatives
!
      call cosmoaiizero(ipts,in,nearsasrptr,npwt1,adrvi,a2drv,ladrv,.false.)
    endif
!
!  Self-term derivatives
!
    if (lperiodic.and.lsegsmooth) then
      qsi2 = 2.0_dp*qsi2
      qsipj = 2.0_dp*qsipj
      f1(1:3) = 0.0_dp
      f2(1:6) = 0.0_dp
      call cosmoamatdadd1p(ipts,ipts,qmeself,f1,f2,qsi2,qsipj,sumAinvipts,sumAinvipts,ldqneeded, &
        dcosmoA,dsegweight,nearsas,nearsasrptr,dcosmoA2)
    endif
!
!  Set constants for first SAS point
!
    xi = xclat(i)
    yi = yclat(i)
    zi = zclat(i)
    xa = spxyz(1,ipts) 
    ya = spxyz(2,ipts) 
    za = spxyz(3,ipts) 
!
!  Calculation terms for A-matrix between SAS points
!
    do jpts = 1,npts
      j = cosmoatomptr(jpts)
      ja = nrelf2a(j)
      npwt2 = npwt(jpts)
      jn = nearsasrptr(j)
      nregionj = nregionno(nsft+ja)
      lopj = (.not.lfreeze.or.lopf(j).or.lnearsas(j))
      if (lopi.or.lopj) then
        swj = segweight(jpts)
        narj = nar(jpts)
        nsetfj = nsetf(jpts)
        qsj = qonsas(jpts)
        qsij = qsi*qsj*fact
        sumAinvjpts = 0.0_dp
        qsipj = 0.0_dp
        xji = spxyz(1,jpts) - xa
        yji = spxyz(2,jpts) - ya
        zji = spxyz(3,jpts) - za
!
!  Get derivatives of weighting factors
!
        call dtotcosmowt(jpts,totalwt(jpts),dwtj,d2wtj,.false.)
!
!  Loop over cell images
!
        do ii = 1,nvec
          if (ii.ne.nmid.or.jpts.ne.ipts) then
            xij = xji + xvec(ii)
            yij = yji + yvec(ii)
            zij = zji + zvec(ii)
            dists2 = xij*xij + yij*yij + zij*zij
!
            if (dists2 .lt. cosmormax2) then
              if (dists2 .gt. cosmormax2s) then
                lrweight = .true.
                dists = sqrt(dists2)
                call sasswitch2(cosmormax-dists,dwt,.true.,.false.,dwdr,d2wdr2)
!
!  Switch sign of dwdr since the negative of the distance was passed into switch
!
                dwdr = - dwdr
!
                rdists = 1.0_dp/dists
                dwtm1 = 1.0_dp - dwt
                dwdr  = rdists*dwdr
                d2wdr2 = rdists*rdists*(d2wdr2 - dwdr)
!
                if (lPureCoulomb0D) then
                  f0    = dwt*rdists
                  rtrm = rdists*(dwdr - dwt*rdists*rdists)
                  f1(1) = rtrm*xij
                  f1(2) = rtrm*yij
                  f1(3) = rtrm*zij
                  rtrm2 = rdists*(d2wdr2 + rdists*rdists*(3.0_dp*dwt*rdists*rdists - 2.0_dp*dwdr))
                  f2(1) = rtrm2*xij*xij + rtrm
                  f2(2) = rtrm2*xij*yij
                  f2(3) = rtrm2*yij*yij + rtrm
                  f2(4) = rtrm2*xij*zij
                  f2(5) = rtrm2*yij*zij
                  f2(6) = rtrm2*zij*zij + rtrm
                else
                  call qmatrixelementc(xij,yij,zij,0.0_dp,.false.,.true.,.false.,aij2,dqme,d2qme)
                  f0 = dwt*aij2
                  f1(1) = - dwt*dqme(1) + dwdr*aij2*xij
                  f1(2) = - dwt*dqme(2) + dwdr*aij2*yij
                  f1(3) = - dwt*dqme(3) + dwdr*aij2*zij
                  f2(1) = - dwt*d2qme(1) + d2wdr2*aij2*xij*xij - 2.0_dp*dwdr*xij*dqme(1) + dwdr*aij2
                  f2(2) = - dwt*d2qme(2) + d2wdr2*aij2*xij*yij - dwdr*(xij*dqme(2) + yij*dqme(1))
                  f2(3) = - dwt*d2qme(3) + d2wdr2*aij2*yij*yij - 2.0_dp*dwdr*yij*dqme(2) + dwdr*aij2
                  f2(4) = - dwt*d2qme(4) + d2wdr2*aij2*xij*zij - dwdr*(xij*dqme(3) + zij*dqme(1))
                  f2(5) = - dwt*d2qme(5) + d2wdr2*aij2*yij*zij - dwdr*(zij*dqme(2) + yij*dqme(3))
                  f2(6) = - dwt*d2qme(6) + d2wdr2*aij2*zij*zij - 2.0_dp*dwdr*zij*dqme(3) + dwdr*aij2
                endif
                call cosmoamatdadd1p(ipts,jpts,f0,f1,f2,qsij,qsipj,sumAinvipts,sumAinvjpts,ldqneeded, &
                  dcosmoA,dsegweight,nearsas,nearsasrptr,dcosmoA2)
              else
                lrweight = .false.
                dwdr = 0.0_dp
                d2wdr2 = 0.0_dp
                dwt = 0.0_dp
                dwtm1 = 1.0_dp
              endif
              xj = xclat(j) - xi + xvec(ii)
              yj = yclat(j) - yi + yvec(ii)
              zj = zclat(j) - zi + zvec(ii)
              rj = atsrad(j)
              do k = 1,nari
                j1 = nset(k+nsetfi-1)
                xc = sphere2(1,j1)*ri
                yc = sphere2(2,j1)*ri
                zc = sphere2(3,j1)*ri
                w1 = cosmowt(j1,i)
                pw1 = cosmopwt(k+nsetfi-1)
                w1 = w1*pw1*totalwt(ipts)
                if (i .ne. j .or. ii.ne.nmid) then
                  x1 = xc*cosmotm(1,1,i) + yc*cosmotm(2,1,i) + zc*cosmotm(3,1,i) - xj
                  x2 = xc*cosmotm(1,2,i) + yc*cosmotm(2,2,i) + zc*cosmotm(3,2,i) - yj
                  x3 = xc*cosmotm(1,3,i) + yc*cosmotm(2,3,i) + zc*cosmotm(3,3,i) - zj
                  do l = 1,narj
                    j2 = nset(l+nsetfj-1)
                    w2 = cosmowt(j2,j)
                    pw2 = cosmopwt(l+nsetfj-1)
                    w2 = w2*pw2*totalwt(jpts)
                    xc = sphere2(1,j2)*rj
                    yc = sphere2(2,j2)*rj
                    zc = sphere2(3,j2)*rj
                    y1 = xc*cosmotm(1,1,j) + yc*cosmotm(2,1,j) + zc*cosmotm(3,1,j) - x1
                    y2 = xc*cosmotm(1,2,j) + yc*cosmotm(2,2,j) + zc*cosmotm(3,2,j) - x2
                    y3 = xc*cosmotm(1,3,j) + yc*cosmotm(2,3,j) + zc*cosmotm(3,3,j) - x3
                    dist = y1*y1 + y2*y2 + y3*y3
                    rdist = 1.0_dp/sqrt(dist)
                    d1rdist = - rdist*rdist*rdist
                    call cosmoaijderv(ipts,in,k,jpts,jn,l,nearsasrptr,npwt1,npwt2,maxnpwt2, &
                      y1,y2,y3,rdist,d1rdist,w1,w2,dwti,dwtj,d2wti,d2wtj,adrvi,adrvj, &
                      a2drv,ladrv,.false.)
!
!  Derivatives with respect to smoothing
!
                    ddist = w1*w2*d1rdist
                    ff = dwtm1*ddist      
                    f0    = dwtm1*w1*w2*rdist
                    f1(1) = ff*y1 
                    f1(2) = ff*y2
                    f1(3) = ff*y3
                    if (lrweight) then
                      ffa = - w1*w2*rdist*dwdr
                      f1(1) = f1(1) + ffa*xij
                      f1(2) = f1(2) + ffa*yij
                      f1(3) = f1(3) + ffa*zij
                    endif
                    rtrm  = - ddist*dwtm1
                    rtrm2 = - 3.0_dp*rdist*rdist*ff
                    f2(1) = rtrm2*y1*y1 - rtrm
                    f2(2) = rtrm2*y1*y2  
                    f2(3) = rtrm2*y2*y2 - rtrm
                    f2(4) = rtrm2*y1*y3       
                    f2(5) = rtrm2*y2*y3     
                    f2(6) = rtrm2*y3*y3 - rtrm
                    if (lrweight) then
                      rtrm2 = - ddist*dwdr
                      f2(1) = f2(1) + rtrm2*(y1*xij + xij*y1)
                      f2(2) = f2(2) + rtrm2*(y2*xij + yij*y1)
                      f2(3) = f2(3) + rtrm2*(y2*yij + yij*y2)
                      f2(4) = f2(4) + rtrm2*(y3*xij + zij*y1)
                      f2(5) = f2(5) + rtrm2*(y3*yij + zij*y2)
                      f2(6) = f2(6) + rtrm2*(y3*zij + zij*y3)
!
                      rtrm2 = - w1*w2*rdist*d2wdr2
                      rtrm  = - w1*w2*rdist*dwdr
                      f2(1) = f2(1) + rtrm2*xij*xij + rtrm
                      f2(2) = f2(2) + rtrm2*xij*yij
                      f2(3) = f2(3) + rtrm2*yij*yij + rtrm
                      f2(4) = f2(4) + rtrm2*xij*zij
                      f2(5) = f2(5) + rtrm2*yij*zij
                      f2(6) = f2(6) + rtrm2*zij*zij + rtrm
                    endif
                    call cosmoamatdadd1p(ipts,jpts,f0,f1,f2,qsij,qsipj,sumAinvipts,sumAinvjpts,ldqneeded, &
                      dcosmoA,dsegweight,nearsas,nearsasrptr,dcosmoA2)
                  enddo
                else
                  do l = 1,narj
                    j2 = nset(l+nsetfj-1)
                    w2 = cosmowt(j2,j)
                    pw2 = cosmopwt(l+nsetfj-1)
                    w2 = w2*pw2*totalwt(jpts)
                    y1 = sphere2(1,j2)*rj - xc
                    y2 = sphere2(2,j2)*rj - yc
                    y3 = sphere2(3,j2)*rj - zc
                    dist = y1*y1 + y2*y2 + y3*y3                                                                           
                    if (dist.gt.1.0d-8) then
                      rdist = 1.0_dp/sqrt(dist)
                      d1rdist = 0.0_dp
                      call cosmoaijderv(ipts,in,k,jpts,jn,l,nearsasrptr,npwt1,npwt2,maxnpwt2, &
                        y1,y2,y3,rdist,d1rdist,w1,w2,dwti,dwtj,d2wti,d2wtj,adrvi,adrvj, &
                        a2drv,ladrv,.false.)
!
!  Derivatives with respect to smoothing
!
                      f0 = dwtm1*w1*w2*rdist
                      f1(1:3) = 0.0_dp
                      if (lrweight) then
                        ffa = - w1*w2*rdist*dwdr
                        f1(1) = f1(1) + ffa*xij
                        f1(2) = f1(2) + ffa*yij
                        f1(3) = f1(3) + ffa*zij
                      endif
                      f2(1:6) = 0.0_dp
                      if (lrweight) then
                        rtrm2 = - w1*w2*rdist*d2wdr2
                        rtrm  = - w1*w2*rdist*dwdr
                        f2(1) = f2(1) + rtrm2*xij*xij + rtrm
                        f2(2) = f2(2) + rtrm2*xij*yij
                        f2(3) = f2(3) + rtrm2*yij*yij + rtrm
                        f2(4) = f2(4) + rtrm2*xij*zij
                        f2(5) = f2(5) + rtrm2*yij*zij
                        f2(6) = f2(6) + rtrm2*zij*zij + rtrm
                      endif
                      call cosmoamatdadd1p(ipts,jpts,f0,f1,f2,qsij,qsipj,sumAinvipts,sumAinvjpts,ldqneeded, &
                        dcosmoA,dsegweight,nearsas,nearsasrptr,dcosmoA2)
                    endif
                  enddo
                endif
              enddo
!
!  Add contributions of weighting factors to overall derivatives
!                     
              dwtm1swiswj = dwtm1*swi*swj
              qsijdwtm1swiswj = qsij*dwtm1swiswj
              qsipjdwtm1swiswj = qsipj*dwtm1swiswj
              do nn = 1,nearsas
                if (ladrv(nn)) then
                  n = nearsasptr(nn) 
                  nregionn = nregionno(nsft+nrelf2a(n))
                  if (n .ne. i) then
                    xdrv(i) = xdrv(i) - qsijdwtm1swiswj*adrvi(1,nn)
                    ydrv(i) = ydrv(i) - qsijdwtm1swiswj*adrvi(2,nn)
                    zdrv(i) = zdrv(i) - qsijdwtm1swiswj*adrvi(3,nn)
                    xdrv(n) = xdrv(n) + qsijdwtm1swiswj*adrvi(1,nn)
                    ydrv(n) = ydrv(n) + qsijdwtm1swiswj*adrvi(2,nn)
                    zdrv(n) = zdrv(n) + qsijdwtm1swiswj*adrvi(3,nn)
                    if (nregioni.ne.nregionn) then
                      xregdrv(nregioni) = xregdrv(nregioni) - qsijdwtm1swiswj*adrvi(1,nn)
                      yregdrv(nregioni) = yregdrv(nregioni) - qsijdwtm1swiswj*adrvi(2,nn)
                      zregdrv(nregioni) = zregdrv(nregioni) - qsijdwtm1swiswj*adrvi(3,nn)
                      xregdrv(nregionn) = xregdrv(nregionn) + qsijdwtm1swiswj*adrvi(1,nn)
                      yregdrv(nregionn) = yregdrv(nregionn) + qsijdwtm1swiswj*adrvi(2,nn)
                      zregdrv(nregionn) = zregdrv(nregionn) + qsijdwtm1swiswj*adrvi(3,nn)
                    endif
                  endif
!
!  Save first derivatives of A matrix term
!
                  if (ldqneeded) then
                    dcosmoA(1,iptsloc,nn) = dcosmoA(1,iptsloc,nn) - qsj*adrvi(1,nn)*dwtm1swiswj
                    dcosmoA(2,iptsloc,nn) = dcosmoA(2,iptsloc,nn) - qsj*adrvi(2,nn)*dwtm1swiswj
                    dcosmoA(3,iptsloc,nn) = dcosmoA(3,iptsloc,nn) - qsj*adrvi(3,nn)*dwtm1swiswj
                  endif
                  if (n.ne.j) then
                    xdrv(j) = xdrv(j) - qsijdwtm1swiswj*adrvj(1,nn)
                    ydrv(j) = ydrv(j) - qsijdwtm1swiswj*adrvj(2,nn)
                    zdrv(j) = zdrv(j) - qsijdwtm1swiswj*adrvj(3,nn)
                    xdrv(n) = xdrv(n) + qsijdwtm1swiswj*adrvj(1,nn)
                    ydrv(n) = ydrv(n) + qsijdwtm1swiswj*adrvj(2,nn)
                    zdrv(n) = zdrv(n) + qsijdwtm1swiswj*adrvj(3,nn)
                    if (nregionj.ne.nregionn) then
                      xregdrv(nregionj) = xregdrv(nregionj) - qsijdwtm1swiswj*adrvj(1,nn)
                      yregdrv(nregionj) = yregdrv(nregionj) - qsijdwtm1swiswj*adrvj(2,nn)
                      zregdrv(nregionj) = zregdrv(nregionj) - qsijdwtm1swiswj*adrvj(3,nn)
                      xregdrv(nregionn) = xregdrv(nregionn) + qsijdwtm1swiswj*adrvj(1,nn)
                      yregdrv(nregionn) = yregdrv(nregionn) + qsijdwtm1swiswj*adrvj(2,nn)
                      zregdrv(nregionn) = zregdrv(nregionn) + qsijdwtm1swiswj*adrvj(3,nn)
                    endif
                  endif
!
!  Save first derivatives of A matrix term
!
                  if (ldqneeded) then
                    dcosmoA(1,iptsloc,jn) = dcosmoA(1,iptsloc,jn) + qsj*adrvj(1,nn)*dwtm1swiswj
                    dcosmoA(2,iptsloc,jn) = dcosmoA(2,iptsloc,jn) + qsj*adrvj(2,nn)*dwtm1swiswj
                    dcosmoA(3,iptsloc,jn) = dcosmoA(3,iptsloc,jn) + qsj*adrvj(3,nn)*dwtm1swiswj
                    dcosmoA(1,iptsloc,nn) = dcosmoA(1,iptsloc,nn) - qsj*adrvj(1,nn)*dwtm1swiswj
                    dcosmoA(2,iptsloc,nn) = dcosmoA(2,iptsloc,nn) - qsj*adrvj(2,nn)*dwtm1swiswj
                    dcosmoA(3,iptsloc,nn) = dcosmoA(3,iptsloc,nn) - qsj*adrvj(3,nn)*dwtm1swiswj
                  endif
                endif
              enddo
!
!  For periodic case, subtract 1/r term to avoid double counting with the term from the long range sum.
!
              if (lperiodic.and.dists2.gt.1.0d-15) then
                call qmatrixelementc(xij,yij,zij,0.0_dp,.false.,.true.,.false.,aij,dqme,d2qme)
                f0    = - aij
                f1(1:3) = dqme(1:3)
                f2(1:6) = d2qme(1:6)
                call cosmoamatdadd1p(ipts,jpts,f0,f1,f2,qsij,qsipj,sumAinvipts,sumAinvjpts,ldqneeded, &
                  dcosmoA,dsegweight,nearsas,nearsasrptr,dcosmoA2)
              endif
!
!  Rezero derivatives
!
              call cosmoaijzero(ipts,in,jpts,jn,nearsasrptr,npwt1,npwt2,adrvi,adrvj,a2drv,ladrv,.false.)
            else
!
!  Coulomb element only needs adding for a cluster since this is handled for periodic systems via qmatrixelementc.
!
              if (.not.lperiodic.and.dists2.gt.1.0d-8) then
                if (lPureCoulomb0D) then
                  aij = 1.0_dp/sqrt(dists2)
                  f0    = aij
                  rtrm  = aij**3
                  f1(1) = - rtrm*xij
                  f1(2) = - rtrm*yij
                  f1(3) = - rtrm*zij
                  f2(1) = rtrm*(3.0_dp*xij*xij*aij**2 - 1.0_dp)
                  f2(2) = rtrm*3.0_dp*xij*yij*aij**2
                  f2(3) = rtrm*(3.0_dp*yij*yij*aij**2 - 1.0_dp)
                  f2(4) = rtrm*3.0_dp*xij*zij*aij**2
                  f2(5) = rtrm*3.0_dp*yij*zij*aij**2
                  f2(6) = rtrm*(3.0_dp*zij*zij*aij**2 - 1.0_dp)
                else
                  call qmatrixelementc(xji,yji,zji,0.0_dp,.false.,.true.,.false.,qme,dqme,d2qme)
                  f0 = qme
                  f1(1:3) = - dqme(1:3)
                  f2(1:6) = - d2qme(1:6)
                endif
                call cosmoamatdadd1p(ipts,jpts,f0,f1,f2,qsij,qsipj,sumAinvipts,sumAinvjpts,ldqneeded, &
                  dcosmoA,dsegweight,nearsas,nearsasrptr,dcosmoA2)
              endif
            endif
          endif
!
!  End loop over cell images
!
        enddo
!
!  Long range contribution to derivatives for periodic systems
!
        if (lperiodic.and.(ipts.ne.jpts)) then
          call qmatrixelementc(xji,yji,zji,0.0_dp,.true.,.true.,.false.,qme,dqme,d2qme)
          f0 = qme
          f1(1:3) = - dqme(1:3)
          f2(1:6) = - d2qme(1:6)
          call cosmoamatdadd1p(ipts,jpts,f0,f1,f2,qsij,qsipj,sumAinvipts,sumAinvjpts,ldqneeded, &
            dcosmoA,dsegweight,nearsas,nearsasrptr,dcosmoA2)
        endif
!
!  End if for lopi/lopj
!
      endif
!
!  End loops over segment points
!
    enddo
  enddo
!
!  Deallocate array for local derivative storage
!
  call realloc(d2wtj,0_i4,0_i4,0_i4,0_i4,ierror)
  call realloc(d2wti,0_i4,0_i4,0_i4,0_i4,ierror)
  call realloc(a2drv,0_i4,0_i4,0_i4,ierror)
  deallocate(ladrv,stat=status)
  if (status/=0) call deallocate_error('setdcosmoamat','ladrv')
  call realloc(adrvj,0_i4,0_i4,ierror)
  call realloc(adrvi,0_i4,0_i4,ierror)
  call realloc(dwtj,0_i4,0_i4,0_i4,ierror)
  call realloc(dwti,0_i4,0_i4,0_i4,ierror)
!
!  Free local memory
!
  deallocate(totalwt,stat=status)
  if (status/=0) call deallocate_error('setdcosmoamat','totalwt')
!
  return
  end
