  subroutine eckart(mcv,mcvloc,mcvptr,maxd2,derv2,maxeigr,eigr)
!
!  Performs the Eckart purification of the second derivatives for clusters
!  Distributed second derivative version. 
!
!  12/16 Created from deffreq
!   8/17 Parallel options wrapped with ifdef 
!   2/18 Trace added
!   5/20 Correction to calculation of moment of inertia tensor
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
!  Julian Gale, CIC, Curtin University, May 2020
!
  use g_constants
  use control
  use current
  use element
  use gulp_files
  use frequencies
  use general
  use iochannels
  use parallel
  use phononatoms
  use species,        only : massspec
  use shells
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4),  intent(in)                          :: mcv
  integer(i4),  intent(in)                          :: mcvloc
  integer(i4),  intent(in)                          :: mcvptr(mcvloc)
  integer(i4),  intent(in)                          :: maxd2
  integer(i4),  intent(in)                          :: maxeigr
  real(dp),     intent(inout)                       :: derv2(maxd2,*)
  real(dp),     intent(out)                         :: eigr(maxeigr,*)
!
!  Local variables
!
  integer(i4)                                       :: i
  integer(i4)                                       :: ifail
  integer(i4)                                       :: j
  integer(i4)                                       :: k
  integer(i4)                                       :: nrottran
  integer(i4)                                       :: status
#ifdef MPI
  integer                                           :: idesd(9)
  integer                                           :: idese(9)
  integer                                           :: idest(9)
  integer                                           :: idesw(9)
  integer                                           :: ifails
  integer                                           :: ld
  integer                                           :: nb
  integer                                           :: ncs
#endif
!
  real(dp)                                          :: dinert(3,3)
  real(dp)                                          :: dx
  real(dp)                                          :: dy
  real(dp)                                          :: dz
  real(dp)                                          :: mi
  real(dp)                                          :: totalmass
  real(dp),     dimension(:,:),   allocatable       :: tr_vectors
  real(dp),     dimension(:),     allocatable       :: w1
  real(dp),     dimension(:,:),   allocatable       :: ww1
  real(dp)                                          :: xcom
  real(dp)                                          :: ycom
  real(dp)                                          :: zcom
  real(dp)                                          :: xi
  real(dp)                                          :: yi
  real(dp)                                          :: zi
!
!  Check the system is finite
!
  if (ndim.ne.0) return
#ifdef TRACE
  call trace_in('eckart')
#endif
!
!  Find centre of mass
!
  xcom = 0.0_dp
  ycom = 0.0_dp
  zcom = 0.0_dp
  totalmass = 0.0_dp
  do i = 1,ncore
    mi = occuf(i)*massspec(nspecptr(i))
    totalmass = totalmass + mi
    xcom = xcom + xclat(i)*mi
    ycom = ycom + yclat(i)*mi
    zcom = zcom + zclat(i)*mi
  enddo
  xcom = xcom/totalmass
  ycom = ycom/totalmass
  zcom = zcom/totalmass
!
!  Calculate inertia tensor
!
  dinert(1:3,1:3) = 0.0_dp
  do i = 1,ncore
    mi = occuf(i)*massspec(nspecptr(i))
    dx = xclat(i) - xcom
    dy = yclat(i) - ycom
    dz = zclat(i) - zcom
    dinert(1,1) = dinert(1,1) + mi*(dy*dy + dz*dz)
    dinert(2,1) = dinert(2,1) - mi*dy*dx
    dinert(3,1) = dinert(3,1) - mi*dz*dx
    dinert(1,2) = dinert(1,2) - mi*dx*dy
    dinert(2,2) = dinert(2,2) + mi*(dx*dx + dz*dz)
    dinert(3,2) = dinert(3,2) - mi*dz*dy
    dinert(1,3) = dinert(1,3) - mi*dx*dz
    dinert(2,3) = dinert(2,3) - mi*dy*dz
    dinert(3,3) = dinert(3,3) + mi*(dx*dx + dy*dy)
  enddo
!
  allocate(w1(3*mcv),stat=status)
  if (status/=0) call outofmemory('eckart','w1')
!
  call dsyev('V','U',3_i4,dinert,3_i4,freq,w1,3_i4*mcv,ifail)
!
  deallocate(w1,stat=status)
  if (status/=0) call deallocate_error('eckart','w1')
!
!  Find whether molecule is linear or not by looking at principal components of inertia tensor
!
  if (freq(1,1).lt.1.0d-6.or.ncore.eq.2) then
    nrottran = 5_i4
  else
    nrottran = 6_i4
  endif
!
!  Allocate workspace for translational and rotational vectors
!
  allocate(tr_vectors(mcv,6),stat=status)
  if (status/=0) call outofmemory('eckart','tr_vectors')
!
!  Create translation vectors
!
  do i = 1,ncore
    tr_vectors(3*(i-1)+1,1) = 1.0_dp
    tr_vectors(3*(i-1)+2,1) = 0.0_dp
    tr_vectors(3*(i-1)+3,1) = 0.0_dp
    tr_vectors(3*(i-1)+1,2) = 0.0_dp
    tr_vectors(3*(i-1)+2,2) = 1.0_dp
    tr_vectors(3*(i-1)+3,2) = 0.0_dp
    tr_vectors(3*(i-1)+1,3) = 0.0_dp
    tr_vectors(3*(i-1)+2,3) = 0.0_dp
    tr_vectors(3*(i-1)+3,3) = 1.0_dp
  enddo
!
!  Create rotation vectors
!
  if (nrottran.eq.6) then
    do i = 1,ncore
      xi = (xclat(i) - xcom)
      yi = (yclat(i) - ycom)
      zi = (zclat(i) - zcom)
!
      tr_vectors(3*(i-1)+1,4) = yi*dinert(3,1) - zi*dinert(2,1)
      tr_vectors(3*(i-1)+2,4) = zi*dinert(1,1) - xi*dinert(3,1)
      tr_vectors(3*(i-1)+3,4) = xi*dinert(2,1) - yi*dinert(1,1)
      tr_vectors(3*(i-1)+1,5) = yi*dinert(3,2) - zi*dinert(2,2)
      tr_vectors(3*(i-1)+2,5) = zi*dinert(1,2) - xi*dinert(3,2)
      tr_vectors(3*(i-1)+3,5) = xi*dinert(2,2) - yi*dinert(1,2)
      tr_vectors(3*(i-1)+1,6) = yi*dinert(3,3) - zi*dinert(2,3)
      tr_vectors(3*(i-1)+2,6) = zi*dinert(1,3) - xi*dinert(3,3)
      tr_vectors(3*(i-1)+3,6) = xi*dinert(2,3) - yi*dinert(1,3)
    enddo
  else
    do i = 1,ncore
      xi = (xclat(i) - xcom)
      yi = (yclat(i) - ycom)
      zi = (zclat(i) - zcom)
!
      tr_vectors(3*(i-1)+1,4) = yi*dinert(3,2) - zi*dinert(2,2)
      tr_vectors(3*(i-1)+2,4) = zi*dinert(1,2) - xi*dinert(3,2)
      tr_vectors(3*(i-1)+3,4) = xi*dinert(2,2) - yi*dinert(1,2)
      tr_vectors(3*(i-1)+1,5) = yi*dinert(3,3) - zi*dinert(2,3)
      tr_vectors(3*(i-1)+2,5) = zi*dinert(1,3) - xi*dinert(3,3)
      tr_vectors(3*(i-1)+3,5) = xi*dinert(2,3) - yi*dinert(1,3)
    enddo
  endif
!
!  Normalise
!
  do i = 1,nrottran
    mi = 0.0_dp
    do j = 1,mcv
      mi = mi + tr_vectors(j,i)**2
    enddo
    mi = 1.0_dp/sqrt(mi)
    do j = 1,mcv
      tr_vectors(j,i) = tr_vectors(j,i)*mi
    enddo
  enddo
!
!  G-S orthogonalisation of rotations to translations
!
  do i = 4,nrottran
    do j = 1,i-1
      mi = 0.0_dp
      do k = 1,mcv
        mi = mi + tr_vectors(k,i)*tr_vectors(k,j)
      enddo
      do k = 1,mcv
        tr_vectors(k,i) = tr_vectors(k,i) - mi*tr_vectors(k,j)
      enddo
    enddo
!
!  Renormalise
!
    mi = 0.0_dp
    do j = 1,mcv
      mi = mi + tr_vectors(j,i)**2
    enddo
    mi = 1.0_dp/sqrt(mi)
    do j = 1,mcv
      tr_vectors(j,i) = tr_vectors(j,i)*mi
    enddo
  enddo
!
  allocate(ww1(mcv,mcvloc),stat=status)
  if (status/=0) call outofmemory('eckart','ww1')
!
  eigr(1:mcv,1:mcvloc) = derv2(1:mcv,1:mcvloc)
!
!  Eckart transformation to remove the rotational and translation modes : (1-P)*H*(1-P)
!
  ww1(1:mcv,1:mcvloc) = 0.0_dp
  do i = 1,mcvloc
    ww1(mcvptr(i),i) = 1.0_dp
  enddo
#ifdef MPI
  if (nprocs.gt.1) then
!
!  Set up Blacs descriptors for matrices
!
    nb = nblocksize
    ifails = 0
    ncs = mcv
    ld = maxd2
    call descinit( idesd, ncs, ncs, 3*nb, 3*nb, 0, 0, iBlacsContext, ld, ifails )
!
    if (ifails.ne.0) then
      call outerror('initialisation in descinit failed',0_i4)
      call stopnow('eckart')
    endif
!
    ld = maxeigr
    call descinit( idese, ncs, ncs, 3*nb, 3*nb, 0, 0, iBlacsContext, ld, ifails )
!
    if (ifails.ne.0) then
      call outerror('initialisation in descinit failed',0_i4)
      call stopnow('eckart')
    endif
!
    ld = mcv
    call descinit( idesw, ncs, ncs, 3*nb, 3*nb, 0, 0, iBlacsContext, ld, ifails )
!
    if (ifails.ne.0) then
      call outerror('initialisation in descinit failed',0_i4)
      call stopnow('eckart')
    endif
!
    ld = mcv
    call descinit( idest, ncs, 6, 3*nb, 6, 0, 0, iBlacsContext, ld, ifails )
!
    if (ifails.ne.0) then
      call outerror('initialisation in descinit failed',0_i4)
      call stopnow('eckart')
    endif
!
    call pdgemm('n','t',mcv,mcv,nrottran,-1.0d0,tr_vectors,1,1,idest,tr_vectors,1,1,idest,1.0d0,ww1,1,1,idesw) ! Creates 1-P where P = v.v^t
    call pdgemm('n','n',mcv,mcv,mcv,1.0_dp,ww1,1,1,idesw,derv2,1,1,idesd,0.0_dp,eigr,1,1,idese)   ! (1-P)*H
    call pdgemm('n','n',mcv,mcv,mcv,1.0_dp,eigr,1,1,idese,ww1,1,1,idesw,0.0_dp,derv2,1,1,idesd)   ! [(1-P)*H]*(1-P)
  else
#endif
    call dgemm('n','t',mcv,mcv,nrottran,-1.0_dp,tr_vectors,mcv,tr_vectors,mcv,1.0_dp,ww1,mcv)   ! Creates 1-P where P = v.v^t
    call dgemm('n','n',mcv,mcv,mcv,1.0_dp,ww1,mcv,derv2,maxd2,0.0_dp,eigr,maxeigr)             ! (1-P)*H
    call dgemm('n','n',mcv,mcv,mcv,1.0_dp,eigr,maxeigr,ww1,mcv,0.0_dp,derv2,maxd2)             ! [(1-P)*H]*(1-P)
#ifdef MPI
  endif
#endif
!
!  Deallocate memory for purification
!
  deallocate(ww1,stat=status)
  if (status/=0) call deallocate_error('eckart','ww1')
  deallocate(tr_vectors,stat=status)
  if (status/=0) call deallocate_error('eckart','tr_vectors')
#ifdef TRACE
  call trace_out('eckart')
#endif
!
  return
  end
