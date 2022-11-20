  subroutine sec0
!
!  Generate symmetry adapted cluster second derivative matrix.
!
!  Freezing now added.
!
!   8/95 Modifications added to avoid using transformation
!        matrix unless constraints are present as this is
!        much slower.
!  10/04 Style updated for new version
!   3/17 fix_atom option added
!   2/18 Trace added
!   4/18 Bug corrected for case where nfixatom is zero due to flags
!   3/19 iopt replaced by ioptindex and iopttype
!  12/19 Correction to referencing of internal derivatives for rigid molecule case
!   3/20 Rigid molecule modifications added
!   5/20 Corrections to molQTdrv handling
!   5/20 Trapping of molecules with no on-diagonal molQQdrv added
!   5/20 ltmat algorithm added for rigid molecules
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
  use control
  use current
  use derivatives
  use iochannels
  use molecule
  use optimisation
  use parallel
  use times
  use transform
#ifdef TRACE
  use trace,        only : trace_in, trace_out
#endif
  implicit none
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ii
  integer(i4)                                  :: indj
  integer(i4)                                  :: indk
  integer(i4)                                  :: ix
  integer(i4)                                  :: j
  integer(i4)                                  :: jj
  integer(i4)                                  :: jx
  integer(i4)                                  :: k
  integer(i4)                                  :: kk
  integer(i4)                                  :: n3f
  integer(i4)                                  :: n3m
  integer(i4)                                  :: ni
  integer(i4)                                  :: nint
  integer(i4)                                  :: nintmin
  integer(i4)                                  :: nj
  integer(i4)                                  :: status
  logical                                      :: lsdebug
  logical                                      :: ltmat
  logical                                      :: lmolQi
  logical                                      :: lmolTi
  logical                                      :: lmolQj
  logical                                      :: lmolTj
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: t1
  real(dp)                                     :: t2
  real(dp), dimension(:), allocatable          :: tmp2
  real(dp)                                     :: tr1
#ifdef TRACE
  call trace_in('sec0')
#endif
!
  t1 = g_cpu_time()
  nint = nvar
  lsdebug = (index(keyword,'derv2').ne.0)
  n3f = 3*numat
  n3m = 3*nmol
  if (nbsm.gt.0) then
    n3f = n3f + numat
  endif
!
!  Work out whether full tmat multiplication is needed - use faster method to handle P1 fully relaxed case
!
  if (ncon.eq.0) then
    ltmat = .false.
  else
    ltmat = .true.
  endif
!
!  Print out full second derivatives
!
  if (lsdebug.and.ioproc) then
    write(ioout,'(/,''  Second Derivative Matrix  :'',/)')
    do i = 1,n3f
      write(ioout,'(2x,9(f9.4))')(derv2(i,j),j=1,(3*numat+nbsmat))
    enddo
    write(ioout,'(/)')
    if (lrigid) then
      write(ioout,'(/,''  Second Derivative Matrix : Rigid molecule translations :'',/)')
      do i = 1,nmol
        do j = 1,nmol
          write(ioout,'(2x,2i5,'' x'',1x,3(f9.4))') i,j,(molTTdrv(1,jj,i,j),jj=1,3)
          write(ioout,'(2x,2i5,'' y'',1x,3(f9.4))') i,j,(molTTdrv(2,jj,i,j),jj=1,3)
          write(ioout,'(2x,2i5,'' z'',1x,3(f9.4))') i,j,(molTTdrv(3,jj,i,j),jj=1,3)
        enddo
      enddo
      write(ioout,'(/,''  Second Derivative Matrix : Rigid molecule quaternions :'',/)')
      do i = 1,nmol
        do j = 1,nmol
          write(ioout,'(2x,2i5,'' x'',1x,3(f9.4))') i,j,(molQQdrv(1,jj,i,j),jj=1,3)
          write(ioout,'(2x,2i5,'' y'',1x,3(f9.4))') i,j,(molQQdrv(2,jj,i,j),jj=1,3)
          write(ioout,'(2x,2i5,'' z'',1x,3(f9.4))') i,j,(molQQdrv(3,jj,i,j),jj=1,3)
        enddo
      enddo
      write(ioout,'(/,''  Second Derivative Matrix : Rigid molecule quaternion-translation :'',/)')
      do i = 1,nmol
        do j = 1,nmol
          write(ioout,'(2x,2i5,'' x'',1x,3(f9.4))') i,j,(molQTdrv(1,jj,i,j),jj=1,3)
          write(ioout,'(2x,2i5,'' y'',1x,3(f9.4))') i,j,(molQTdrv(2,jj,i,j),jj=1,3)
          write(ioout,'(2x,2i5,'' z'',1x,3(f9.4))') i,j,(molQTdrv(3,jj,i,j),jj=1,3)
        enddo
      enddo
      write(ioout,'(/)')
    endif
  endif
!********************************************
!  Rigid molecules - check for linear case  *
!********************************************
  if (lrigid) then
    do i = 1,nmol
      do ix = 1,3
!
!  If molecule is linear or there are no rotational forces that set the diagonal element to be a dummy value
!  so that the Hessian can be inverted
!
        if (abs(molQQdrv(ix,ix,i,i)).lt.0.001_dp) then
          do jx = 1,3
            molQQdrv(jx,ix,i,i) = 0.0_dp
            molQQdrv(ix,jx,i,i) = 0.0_dp
          enddo
          molQQdrv(ix,ix,i,i) = 1.0_dp
        endif
      enddo
    enddo
  endif
!******************************************
!  Internal derivatives :                 *
!  Symmetry transform second derivatives  *
!******************************************
  if (ltmat) then
    if ((2*ninternalatm).le.maxd2.and..not.lrigid) then
      do i = 1,ninternalatm
        do j = 1,n3f
          tr1 = 0.0_dp
          do k = 1,n3f
            tr1 = tr1 + derv2(k,j)*tmat(k,i)
          enddo
          tmat(j,ninternalatm+i) = tr1
        enddo
      enddo
      do i = 1,ninternalatm
        do j = 1,ninternalatm
          derv2(j,i) = 0.0_dp
          do k = 1,n3f
            derv2(j,i) = derv2(j,i) + tmat(k,j)*tmat(k,ninternalatm+i)
          enddo
        enddo
      enddo
    else
      allocate(tmp2(n3f),stat=status)
      if (status/=0) call outofmemory('sec0','tmp2')
      do i = 1,ninternalatm
        do j = 1,n3f
          tmp2(j) = 0.0_dp
          do k = 1,n3f
            tmp2(j) = tmp2(j) + derv2(j,k)*tmat(k,i)
          enddo
        enddo
        do j = 1,ninternalatm
          derv2(j,i) = 0.0_dp
          do k = 1,n3f
            derv2(j,i) = derv2(j,i) + tmp2(k)*tmat(k,j)
          enddo
        enddo
      enddo
      if (lrigid) then
!------------------------------------------
!  Rigid molecule - rigid molecule terms  |
!------------------------------------------
!
!  Translation - translation
!
        do i = 1,ninternalmolT
          ii = ninternalatm + i
          tmp2(1:n3m) = 0.0_dp
          do j = 1,nmol
            indj = 3*(j-1)
            do jj = 1,3
              do k = 1,nmol
                indk = 3*(k-1)
                do kk = 1,3
                  tmp2(indj+jj) = tmp2(indj+jj) + molTTdrv(jj,kk,j,k)*tmatT(indk+kk,i)
                enddo
              enddo
            enddo
          enddo
!
          do j = 1,ninternalmolT
            jj = ninternalatm + j
            derv2(jj,ii) = 0.0_dp
            do k = 1,nmol
              indk = 3*(k-1)
              do kk = 1,3
                derv2(jj,ii) = derv2(jj,ii) + tmp2(indk+kk)*tmatT(indk+kk,j)
              enddo
            enddo
          enddo
        enddo
!
!  Quaternion - translation
!
        do i = 1,ninternalmolT
          ii = ninternalatm + i
          tmp2(1:n3m) = 0.0_dp
          do j = 1,nmol
            indj = 3*(j-1)
            do jj = 1,3
              do k = 1,nmol
                indk = 3*(k-1)
                do kk = 1,3
                  tmp2(indj+jj) = tmp2(indj+jj) + molQTdrv(jj,kk,j,k)*tmatT(indk+kk,i)
                enddo
              enddo
            enddo
          enddo
!
          do j = 1,ninternalmolQ
            jj = ninternalatm + ninternalmolT + j
            derv2(ii,jj) = 0.0_dp
            derv2(jj,ii) = 0.0_dp
            do k = 1,nmol
              indk = 3*(k-1)
              do kk = 1,3
                derv2(ii,jj) = derv2(ii,jj) + tmp2(indk+kk)*tmatQ(indk+kk,j)
                derv2(jj,ii) = derv2(jj,ii) + tmp2(indk+kk)*tmatQ(indk+kk,j)
              enddo
            enddo
          enddo
        enddo
!
!  Quaternion - quaternion
!
        do i = 1,ninternalmolQ
          ii = ninternalatm + ninternalmolT + i
          tmp2(1:n3m) = 0.0_dp
          do j = 1,nmol
            indj = 3*(j-1)
            do jj = 1,3
              do k = 1,nmol
                indk = 3*(k-1)
                do kk = 1,3
                  tmp2(indj+jj) = tmp2(indj+jj) + molQQdrv(jj,kk,j,k)*tmatQ(indk+kk,i)
                enddo
              enddo
            enddo
          enddo
!
          do j = 1,ninternalmolQ
            jj = ninternalatm + ninternalmolT + j
            derv2(jj,ii) = 0.0_dp
            do k = 1,nmol
              indk = 3*(k-1)
              do kk = 1,3
                derv2(jj,ii) = derv2(jj,ii) + tmp2(indk+kk)*tmatQ(indk+kk,j)
              enddo
            enddo
          enddo
        enddo
!--------------------------------
!  Rigid molecule - atom terms  |
!--------------------------------
        do i = 1,ninternalatm
          tmp2(1:n3m) = 0.0_dp
          do j = 1,nmol
            indj = 3*(j-1)
            do jj = 1,3
              do k = 1,n3f
                tmp2(indj+jj) = tmp2(indj+jj) + molTCdrv(k,jj,j)*tmat(k,i)
              enddo
            enddo
          enddo
          do j = 1,ninternalmolT
            jj = ninternalatm + j
            derv2(jj,i) = 0.0_dp
            derv2(i,jj) = 0.0_dp
            do k = 1,nmol
              indk = 3*(k-1)
              do kk = 1,3
                derv2(jj,i) = derv2(jj,i) + tmp2(indk+kk)*tmatT(indk+kk,j)
                derv2(i,jj) = derv2(i,jj) + tmp2(indk+kk)*tmatT(indk+kk,j)
              enddo
            enddo
          enddo
        enddo
!
        do i = 1,ninternalatm
          tmp2(1:n3m) = 0.0_dp
          do j = 1,nmol
            indj = 3*(j-1)
            do jj = 1,3
              do k = 1,n3f
                tmp2(indj+jj) = tmp2(indj+jj) + molQCdrv(k,jj,j)*tmat(k,i)
              enddo
            enddo
          enddo
          do j = 1,ninternalmolQ
            jj = ninternalatm + ninternalmolT + j
            derv2(jj,i) = 0.0_dp
            derv2(i,jj) = 0.0_dp
            do k = 1,nmol
              indk = 3*(k-1)
              do kk = 1,3
                derv2(jj,i) = derv2(jj,i) + tmp2(indk+kk)*tmatQ(indk+kk,j)
                derv2(i,jj) = derv2(i,jj) + tmp2(indk+kk)*tmatQ(indk+kk,j)
              enddo
            enddo
          enddo
        enddo
      endif
!
      deallocate(tmp2,stat=status)
      if (status/=0) call deallocate_error('sec0','tmp2')
    endif
  elseif (nfixatom.ne.0.and.nint.eq.(n3f-3).and..not.lrigid) then
!
!  All atoms are variable bar one & so only need to copy from fixed atom onwards
!
    nintmin = 3*(nfixatom-1)
    do i = nintmin+1,nint
      do j = nintmin+1,nint
        derv2(j,i) = derv2(j+3,i+3)
      enddo
    enddo
    do i = nintmin+1,nint
      do j = 1,nintmin
        derv2(j,i) = derv2(j,i+3)
      enddo
    enddo
    do i = 1,nintmin
      do j = nintmin+1,nint
        derv2(j,i) = derv2(j+3,i)
      enddo
    enddo
  else
    do i = 1,nint
      lmolQi = .false.
      lmolTi = .false.
      ni = ioptindex(i)
      if (iopttype(i).eq.iopt_xf) then
        ii = 3*nasymnomolptr(ni) - 2
      elseif (iopttype(i).eq.iopt_yf) then
        ii = 3*nasymnomolptr(ni) - 1
      elseif (iopttype(i).eq.iopt_zf) then
        ii = 3*nasymnomolptr(ni)
      elseif (iopttype(i).eq.iopt_radius) then
        ii = 3*numat + nasymnomolptr(ni)
      elseif (iopttype(i).eq.iopt_xcom) then
        lmolTi = .true.
        ii = 1
      elseif (iopttype(i).eq.iopt_ycom) then
        lmolTi = .true.
        ii = 2
      elseif (iopttype(i).eq.iopt_zcom) then
        lmolTi = .true.
        ii = 3
      elseif (iopttype(i).eq.iopt_xqtn) then
        lmolQi = .true.
        ii = 1
      elseif (iopttype(i).eq.iopt_yqtn) then
        lmolQi = .true.
        ii = 2
      elseif (iopttype(i).eq.iopt_zqtn) then
        lmolQi = .true.
        ii = 3
      endif
      do j = 1,nint
        lmolQj = .false.
        lmolTj = .false.
        nj = ioptindex(j)
        if (iopttype(j).eq.iopt_xf) then
          jj = 3*nasymnomolptr(nj) - 2
        elseif (iopttype(j).eq.iopt_yf) then
          jj = 3*nasymnomolptr(nj) - 1
        elseif (iopttype(j).eq.iopt_zf) then
          jj = 3*nasymnomolptr(nj)
        elseif (iopttype(j).eq.iopt_radius) then
          jj = 3*numat + nasymnomolptr(nj)
        elseif (iopttype(j).eq.iopt_xcom) then
          lmolTj = .true.
          jj = 1
        elseif (iopttype(j).eq.iopt_ycom) then
          lmolTj = .true.
          jj = 2
        elseif (iopttype(j).eq.iopt_zcom) then
          lmolTj = .true.
          jj = 3
        elseif (iopttype(j).eq.iopt_xqtn) then
          lmolQj = .true.
          jj = 1
        elseif (iopttype(j).eq.iopt_yqtn) then
          lmolQj = .true.
          jj = 2
        elseif (iopttype(j).eq.iopt_zqtn) then
          lmolQj = .true.
          jj = 3
        endif
        if (lmolQi) then
          if (lmolQj) then
            derv2(j,i) = molQQdrv(jj,ii,nj,ni)
          elseif (lmolTj) then
            derv2(j,i) = molQTdrv(ii,jj,ni,nj)
          else
            derv2(j,i) = molQCdrv(jj,ii,ni)
          endif
        elseif (lmolTi) then
          if (lmolQj) then
            derv2(j,i) = molQTdrv(jj,ii,nj,ni)
          elseif (lmolTj) then
            derv2(j,i) = molTTdrv(jj,ii,nj,ni)
          else
            derv2(j,i) = molTCdrv(jj,ii,ni)
          endif
        else
          if (lmolQj) then
            derv2(j,i) = molQCdrv(ii,jj,nj)
          elseif (lmolTj) then
            derv2(j,i) = molTCdrv(ii,jj,nj)
          else
            derv2(j,i) = derv2(jj,ii)
          endif
        endif
      enddo
    enddo
  endif
  if (lsdebug.and.ioproc) then
    write(ioout,'(/,''  Symmetrised Second Derivative Matrix  :'',/)')
    do i = 1,nvar
      write(ioout,'(2x,10(f14.4))')(derv2(i,j),j=1,i)
    enddo
    write(ioout,'(/)')
  endif
!
  t2 = g_cpu_time()
  thes = thes + t2 - t1
#ifdef TRACE
  call trace_out('sec0')
#endif
!
  return
  end
