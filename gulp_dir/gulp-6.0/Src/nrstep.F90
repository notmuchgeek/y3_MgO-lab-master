  subroutine nrstep(pvect,hessian,maxhess,lhess2D,gc,nvar)
!
!  Generates the Newton-Raphson step / search direction from the
!  hessian in lower-half triangular form and the gradient vector.
!
!  10/04 Rewritten to use blas call
!   1/09 Integer datatypes all explicitly declared
!   2/17 nmin removed from arguments 
!   2/17 maxhess & lhess2D added as arguments
!   3/17 Parallelisation added
!   4/17 Modified so that parallel algorithm is only used for 2D
!        Hessian case since fitting currently uses replicated 
!        data for the Hessian
!   5/17 maxvar replaced by nvar since the latter is used to set
!        the dimensions of gc
!   1/18 Trace added
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
!  Julian Gale, CIC, Curtin University, January 2018
!
  use datatypes
  use parallel 
#ifdef TRACE
  use trace,       only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)  :: nvar
  integer(i4), intent(in)  :: maxhess
  logical,     intent(in)  :: lhess2D
  real(dp),    intent(in)  :: gc(*)
  real(dp),    intent(in)  :: hessian(maxhess,*)
  real(dp),    intent(out) :: pvect(*)
#ifdef MPI
!
!  Local variables in Scalapack/Blacs integer precision
!
  integer                  :: ifails
  integer                  :: idesg(9)
  integer                  :: idesh(9)
  integer                  :: ldm
  integer                  :: nb
  integer                  :: np
  integer                  :: nv
#endif
#ifdef TRACE
  call trace_in('nrstep')
#endif
#ifdef MPI
!
  if (nprocs.gt.1.and.lhess2D) then
!**************************************************
!  Call pblas routine for matrix-vector multiply  *
!**************************************************
!
!  Set local block size
!
    nb = nblocksizevar
    np = nprocs
    nv = nvar
    ifails = 0
!
!  Set up Blacs descriptors
!
    ldm = maxhess
    call descinit( idesh, nv, nv, nb, nb, 0, 0, iBlacsContext, ldm, ifails )
!
    if (ifails.ne.0) then
      call outerror('initialisation in descinit failed - 1',0_i4)
      call stopnow('nrstep')
    endif
!
    ldm = nvar
    call descinit( idesg, nv, 1, nb, nb, 0, 0, iBlacsContext, ldm, ifails )
!
    if (ifails.ne.0) then
      call outerror('initialisation in descinit failed - 2',0_i4)
      call stopnow('nrstep')
    endif
!
!  Call pblas 
!
    call pdgemv('N',nvar,nvar,1.0d0,hessian,1,1,idesh,gc,1,1,idesg,1,0.0d0,pvect,1,1,idesg,1)
!
!  Broadcase the result from node 0
!
    call sendall(pvect,nvar,0_i4,"nrstep","pvect")
  else
#endif
!*************************************************
!  Call blas routine for matrix-vector multiply  *
!*************************************************
    if (lhess2D) then
      call dgemv('N',nvar,nvar,1.0_dp,hessian,maxhess,gc,1_i4,0.0_dp,pvect,1_i4)
    else
      call dspmv('U',nvar,1.0_dp,hessian(1,1),gc,1_i4,0.0_dp,pvect,1_i4)
    endif
#ifdef MPI
  endif
#endif
#ifdef TRACE
  call trace_out('nrstep')
#endif
!
  return
  end
