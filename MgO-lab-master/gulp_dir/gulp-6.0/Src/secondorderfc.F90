  subroutine secondorderfc(lshelleliminate)
!
!  Calculates the second order force constant matrix for the cores
!
!  The block diagonal elements of derv2 must be stored, as these are
!  not recalculated in dynamic as the summing of off diagonals is no
!  longer applicable.
!
!   4/15 Created from phonon
!   7/15 Argument added to control whether shell contribution is
!        explicitly eliminated or not
!   7/15 It is assumed that derv3 is already set to contain the diagonal
!        blocks due to a prior call to phonon
!   7/15 Modified so that on-diagonal blocks are correctly generated for
!        the current geometry
!   9/16 nphonatptr no longer allocated or deallocated since it is in a module
!   2/18 Trace added
!   5/19 Finite difference flag split for first and second derivatives
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
!  Julian Gale, CIC, Curtin University, May 2019
!
  use configurations
  use g_constants
  use control
  use current
  use derivatives
  use element
  use gulp_files
  use frequencies
  use general
  use iochannels
  use parallel
  use partial
  use phononatoms
  use properties
  use shells
  use times
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  logical,                           intent(in)     :: lshelleliminate
!
!  Local variables
!
  integer(i4)                                       :: i
  integer(i4)                                       :: ifail
  integer(i4)                                       :: indi
  integer(i4)                                       :: indii
  integer(i4)                                       :: indj
  integer(i4)                                       :: indjj
  integer(i4),  dimension(:),     allocatable       :: ipivot
  integer(i4)                                       :: j
  integer(i4)                                       :: kk
  integer(i4)                                       :: l
  integer(i4)                                       :: maxlim
  integer(i4)                                       :: mcv
  integer(i4)                                       :: mint
  integer(i4)                                       :: msv
  integer(i4)                                       :: status
  logical                                           :: lnoanald2loc
  real(dp),     dimension(:),     allocatable       :: eigr
  real(dp)                                          :: fc
  real(dp)                                          :: wr
  real(dp),     dimension(:),     allocatable       :: wtmp
#ifdef TRACE
  call trace_in('secondorderfc')
#endif
!
!  Set logicals
!
  lnoanald2loc = (lnoanald2)
!************************************************
!  Find number of atoms for phonon calculation  *
!************************************************
!
!  Set up points for atoms involve in phonon calculations
!
  call setphonptr
!
!  Allocate local pointer arrays
!
  lpocc = (nsfoc+ncfoc.ne.nphonat)
!
!  Calculate a few constants to do with the size of the problem
!
  mint = 3*nphonat
  maxlim = mint
  if (nbsmat.gt.0) maxlim = maxlim + nphonat
  msv = 3*nsfoc + nbfoc
  mcv = 3*ncfoc
!
!  Check that maxd2 is greater than or equal to mcv
!
  if (maxd2.lt.mcv) then
    maxd2 = mcv
    call changemaxd2
  endif
!************************************************
!  Allocate second derivative workspace memory  *
!************************************************
  if (lshelleliminate.and.msv.gt.0) then
    allocate(eigr(msv*msv),stat=status)
  else
    allocate(eigr(1),stat=status)
  endif
  if (status/=0) call outofmemory('secondorderfc','eigr')
!
!  Evaluate second derivatives
!
  fc = 0.0_dp
  call energy(fc,.true.,.true.)
!
  if (.not.lnoanald2loc.and..not.lfinitediff2) then
!
!  Store diagonal blocks in derv3 to avoid recalculation
!
    do i = 1,nphonat
      indi = 3*(nphonatptr(i)-1)
      derv3(indi+1,1) = derv2(indi+1,indi+1)
      derv3(indi+2,1) = derv2(indi+2,indi+1)
      derv3(indi+3,1) = derv2(indi+3,indi+1)
      derv3(indi+1,2) = derv2(indi+1,indi+2)
      derv3(indi+2,2) = derv2(indi+2,indi+2)
      derv3(indi+3,2) = derv2(indi+3,indi+2)
      derv3(indi+1,3) = derv2(indi+1,indi+3)
      derv3(indi+2,3) = derv2(indi+2,indi+3)
      derv3(indi+3,3) = derv2(indi+3,indi+3)
    enddo
  endif
!
!  Generate second derivatives
!
  if (lnoanald2loc.or.lfinitediff2) then
    call dynamicn
  else
    call dynamic(0.0_dp,0.0_dp,0.0_dp)
  endif
!
  if (.not.lnoanald2loc.and..not.lfinitediff2) then
!
!  Include diagonal blocks, stored in derv3
!
    do i = 1,nphonat
      indi = 3*(nphonatptr(i) - 1)
      derv2(indi+1,indi+1) = derv2(indi+1,indi+1) + derv3(indi+1,1)
      derv2(indi+2,indi+1) = derv2(indi+2,indi+1) + derv3(indi+2,1)
      derv2(indi+3,indi+1) = derv2(indi+3,indi+1) + derv3(indi+3,1)
      derv2(indi+1,indi+2) = derv2(indi+1,indi+2) + derv3(indi+1,2)
      derv2(indi+2,indi+2) = derv2(indi+2,indi+2) + derv3(indi+2,2)
      derv2(indi+3,indi+2) = derv2(indi+3,indi+2) + derv3(indi+3,2)
      derv2(indi+1,indi+3) = derv2(indi+1,indi+3) + derv3(indi+1,3)
      derv2(indi+2,indi+3) = derv2(indi+2,indi+3) + derv3(indi+2,3)
      derv2(indi+3,indi+3) = derv2(indi+3,indi+3) + derv3(indi+3,3)
    enddo
  endif
!*****************************************************************
!  Compress second derivatives according to partial occupancies  *
!*****************************************************************
  if (lpocc) then
    ncsfoc = ncfoc + nsfoc
    call compressd2(derv2,maxd2,ncfoc,nsfoc,nbfoc,numat,iocptr,ibocptr)
  elseif (numat.ne.nphonat) then
!*********************************************************************
!  Compress full second derivatives down to region 1 only if needed  *
!*********************************************************************
    do i = 1,nphonat
      indi  = 3*(i-1)
      indii = 3*(nphonatptr(i)-1)
      do j = 1,nphonat
        indj  = 3*(j-1)
        indjj = 3*(nphonatptr(j)-1)
        derv2(indj+1,indi+1) = derv2(indjj+1,indii+1)
        derv2(indj+2,indi+1) = derv2(indjj+2,indii+1)
        derv2(indj+3,indi+1) = derv2(indjj+3,indii+1)
        derv2(indj+1,indi+2) = derv2(indjj+1,indii+2)
        derv2(indj+2,indi+2) = derv2(indjj+2,indii+2)
        derv2(indj+3,indi+2) = derv2(indjj+3,indii+2)
        derv2(indj+1,indi+3) = derv2(indjj+1,indii+3)
        derv2(indj+2,indi+3) = derv2(indjj+2,indii+3)
        derv2(indj+3,indi+3) = derv2(indjj+3,indii+3)
      enddo
    enddo
  endif
!
!  Output the uncompressed second derivatives for debugging
!
  if (index(keyword,'force').ne.0.and.index(keyword,'debu').ne.0.and.ioproc) then
    write(ioout,'(/,''  Uncompressed Force Constant matrix :'',/)')
    do i = 1,mcv+msv
      write(ioout,'(12f11.6)')(derv2(j,i),j=1,mcv+msv)
    enddo
  endif
!**********************************
!  Eliminate shell contributions  *
!**********************************
  if (lshelleliminate.and.msv.gt.0) then
!**************************
!  Real Matrix Inversion  *
!**************************
    ifail = 0                   
!     
!  Allocate workspace for inversion
!     
    allocate(ipivot(msv),stat=status)
    if (status/=0) call outofmemory('secondorderfc','ipivot')                  
    allocate(wtmp(3*msv),stat=status)
    if (status/=0) call outofmemory('secondorderfc','wtmp')
!
!  Transfer data to packed storage
!    
    kk = 0
    do i = 1,msv
      do j = 1,i
        kk = kk + 1
        eigr(kk) = derv2(mcv+j,mcv+i)
      enddo
    enddo                                    
!         
!  Factorise matrix
!  
    call dsptrf('U',msv,eigr,ipivot,ifail)
    if (ifail.eq.0) then
!
!  Form inverse
!
      call dsptri('U',msv,eigr,ipivot,wtmp,ifail)
!
!  Transfer data back
!
      kk = 0
      do i = 1,msv
        do j = 1,i
          kk = kk + 1
          derv2(mcv+j,mcv+i) = eigr(kk)
          derv2(mcv+i,mcv+j) = eigr(kk)
        enddo
      enddo
    endif
!
!  Free workspace  
!
    deallocate(wtmp,stat=status)
    if (status/=0) call deallocate_error('secondorderfc','wtmp')
    deallocate(ipivot,stat=status)
    if (status/=0) call deallocate_error('secondorderfc','ipivot')  
!
    if (ifail.ne.0) then
      call outerror('inversion of shell 2nd derivatives failed',0_i4)
      call stopnow('secondorderfc')
    endif
!***********************************************
!  Corrected second derivatives = R - T*S-1*T  *
!***********************************************
!
!  First pass : S-1*T stored in diagonally opposite copy of T
!
    do i = 1,mcv
      do j = 1,msv
        wr = 0.0_dp
        do l = 1,msv
          wr = wr + derv2(j+mcv,l+mcv)*derv2(i,l+mcv)
        enddo
        derv2(mcv+j,i) = wr
      enddo
    enddo
!
!  Second pass : T*(S-1*T)
!
    do i = 1,mcv
      do j = 1,mcv
        wr = 0.0_dp
        do l = 1,msv
          wr = wr - derv2(j,l+mcv)*derv2(mcv+l,i)
        enddo
        derv2(j,i) = derv2(j,i) + wr
      enddo
    enddo
  endif
!****************************
!  End of shell correction  *
!****************************
!
!  If debugging print out dynamical matrix
!
  if (index(keyword,'force').ne.0.and.ioproc) then
    write(ioout,'(/,''  Real Force Constant Matrix :'',/)')
    do i = 1,mcv
      write(ioout,'(12f11.6)')(derv2(j,i),j=1,mcv)
    enddo
  endif
!**************************************************
!  Deallocate second derivative workspace memory  *
!**************************************************
  deallocate(eigr,stat=status)
  if (status/=0) call deallocate_error('secondorderfc','eigr')
#ifdef TRACE
  call trace_out('secondorderfc')
#endif
!
  return
  end
