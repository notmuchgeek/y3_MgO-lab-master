  subroutine getshellmatrices(Dall,maxDall,Dsc,maxDsc)
!
!  Calculates the matrices needed for handling shells in a phonon/force
!  constant calculation. The two matrices needed are:
!
!  Dsc    - second derivatives of shells to cores (and transpose, Dcs)
!  Dss^-1 - inverse of second derivatives of shells to shells
!
!  On entry:
!
!  Dall   = matrix of all second derivatives between cores and shells
!
!  On exit:
!
!  Dsc    = (Dss^-1)*Dsc
!
!  Assumptions:
!
!  - cores are sorted to come before shells
!  - there are no partial occupancies
!  - this is a single region calculation
!
!   5/15 Created from phonon
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
  use configurations
  use control
  use current
  use general
  use shells
#ifdef TRACE
  use trace,       only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4),                          intent(in)  :: maxDall
  integer(i4),                          intent(in)  :: maxDsc
  real(dp),                             intent(in)  :: Dall(maxDall,*)
  real(dp),                             intent(out) :: Dsc(maxDsc,*)
!
!  Local variables
!
  integer(i4)                                       :: i
  integer(i4)                                       :: ifail
  integer(i4),  dimension(:),     allocatable       :: ipivot
  integer(i4)                                       :: j
  integer(i4)                                       :: kk
  integer(i4)                                       :: l
  integer(i4)                                       :: mcv
  integer(i4)                                       :: msv
  integer(i4)                                       :: status
  real(dp),     dimension(:),     allocatable       :: D2tmp
  real(dp),     dimension(:),     allocatable       :: wtmp
!
!  If there are no shells then exit
!
  if (nshell.eq.0) return
#ifdef TRACE
  call trace_in('getshellmatrices')
#endif
!
!  Calculate a few constants to do with the size of the problem
!
  msv = 3*nshell
  mcv = 3*ncore
!************************************************
!  Allocate second derivative workspace memory  *
!************************************************
  allocate(D2tmp(msv*(msv+1)/2),stat=status)
  if (status/=0) call outofmemory('getshellmatrices','D2tmp')
!********************************
!  Invert shell - shell matrix  *
!********************************
  ifail = 0
!     
!  Allocate workspace for inversion
!     
  allocate(ipivot(msv),stat=status)
  if (status/=0) call outofmemory('getshellmatrices','ipivot')                  
  allocate(wtmp(3*msv),stat=status)
  if (status/=0) call outofmemory('getshellmatrices','wtmp')
!
!  Transfer data to packed storage
!    
  kk = 0
  do i = 1,msv
    do j = 1,i
      kk = kk + 1
      D2tmp(kk) = Dall(mcv+j,mcv+i)
    enddo
  enddo                                    
!         
!  Factorise matrix
!  
  call dsptrf('U',msv,D2tmp,ipivot,ifail)
  if (ifail.eq.0) then
!
!  Form inverse
!
    call dsptri('U',msv,D2tmp,ipivot,wtmp,ifail)
  endif
!
!  Free workspace  
!
  deallocate(wtmp,stat=status)
  if (status/=0) call deallocate_error('getshellmatrices','wtmp')
  deallocate(ipivot,stat=status)
  if (status/=0) call deallocate_error('getshellmatrices','ipivot')  
!
  if (ifail.ne.0) then
    call outerror('inversion of shell 2nd derivatives failed',0_i4)
    call stopnow('getshellmatrices')
  endif
!**************************************************************
!  Multiply inverted shell-shell matrix by core-shell matrix  *
!**************************************************************
  Dsc(1:msv,1:mcv) = 0.0_dp
!
  do i = 1,mcv
    kk = 0
    do j = 1,msv
!
!  Off diagonal
!
      do l = 1,j-1
        kk = kk + 1
        Dsc(j,i) = Dsc(j,i) + D2tmp(kk)*Dall(l+mcv,i)
        Dsc(l,i) = Dsc(l,i) + D2tmp(kk)*Dall(j+mcv,i)
      enddo
!
!  On diagonal
!
      kk = kk + 1
      Dsc(j,i) = Dsc(j,i) + D2tmp(kk)*Dall(j+mcv,i)
    enddo
  enddo
!**************************************************
!  Deallocate second derivative workspace memory  *
!**************************************************
  deallocate(D2tmp,stat=status)
  if (status/=0) call deallocate_error('getshellmatrices','D2tmp')
#ifdef TRACE
  call trace_out('getshellmatrices')
#endif
!
  return
  end
