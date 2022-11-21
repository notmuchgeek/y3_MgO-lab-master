  subroutine applyconstraints
!
!  Apply structural constraints 
!
!   3/19 Created from setup as part of changes to iopt
!   5/20 Modified to handle rigid molecules based on xctox0
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
  use configurations
  use current
  use datatypes
  use molecule
  use optimisation
  use parallel
  use reallocate
  use shells
  use species
  use symmetry
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ii
  integer(i4)                                  :: iopt_ix
  integer(i4)                                  :: ix
  integer(i4)                                  :: j
  integer(i4)                                  :: k
  integer(i4), dimension(:,:), allocatable     :: ncount
  integer(i4)                                  :: status
  logical                                      :: ldothis
  logical                                      :: lfound
  real(dp)                                     :: diff
  real(dp)                                     :: xcv
  real(dp),    dimension(:),   allocatable     :: xcon        ! Array for holding constrained values
  real(dp)                                     :: xj
  real(dp)                                     :: xk
#ifdef TRACE
  call trace_in('applyconstraints')
#endif
!
!  Apply constraints
!
  if (ncon.gt.0) then
    allocate(xcon(ncon),stat=status)
    if (status/=0) call outofmemory('applyconstraints','xcon')
!
    do i = 1,ncon
      xcon(i) = 0.0_dp
    enddo
!
    do i = 1,ncon
      if (ncvartyp(i).eq.iopt_cell) then
        if (ncvarind(i).eq.1) then
          xcv = a
        elseif (ncvarind(i).eq.2) then
          xcv = b
        elseif (ncvarind(i).eq.3) then
          xcv = c
        elseif (ncvarind(i).eq.4) then
          xcv = alpha
        elseif (ncvarind(i).eq.5) then
          xcv = beta
        elseif (ncvarind(i).eq.6) then
          xcv = gamma
        endif
      elseif (ncvartyp(i).eq.iopt_strain) then
        xcv = strain(ncvarind(i))
      elseif (ncvartyp(i).eq.iopt_radius) then
        xcv = rada(ncvarind(i))
      elseif (ncvartyp(i).eq.iopt_xf) then
        xcv = xafrac(ncvarind(i))
      elseif (ncvartyp(i).eq.iopt_yf) then
        xcv = yafrac(ncvarind(i))
      elseif (ncvartyp(i).eq.iopt_zf) then
        xcv = zafrac(ncvarind(i))
      elseif (ncvartyp(i).eq.iopt_xcom) then
        xcv = molcoma(1,ncvarind(i))
      elseif (ncvartyp(i).eq.iopt_ycom) then
        xcv = molcoma(2,ncvarind(i))
      elseif (ncvartyp(i).eq.iopt_zcom) then
        xcv = molcoma(3,ncvarind(i))
      elseif (ncvartyp(i).eq.iopt_xqtn) then
        xcv = molQa(1,ncvarind(i))
      elseif (ncvartyp(i).eq.iopt_yqtn) then
        xcv = molQa(2,ncvarind(i))
      elseif (ncvartyp(i).eq.iopt_zqtn) then
        xcv = molQa(3,ncvarind(i))
      endif
      xcon(i) = xcon(i) + xcv*conco(i) + conadd(i)
    enddo
!
    do i = 1,ncon
      if (ncfixtyp(i).eq.iopt_cell) then
        if (ncfixind(i).eq.1) then
          a = xcon(i)
        elseif (ncfixind(i).eq.2) then
          b = xcon(i)
        elseif (ncfixind(i).eq.3) then
          c = xcon(i)
        elseif (ncfixind(i).eq.4) then
          alpha = xcon(i)
        elseif (ncfixind(i).eq.5) then
          beta = xcon(i)
        elseif (ncfixind(i).eq.6) then
          gamma = xcon(i)
        endif
      elseif (ncfixtyp(i).eq.iopt_strain) then
        strain(ncfixind(i)) = xcon(i)
      elseif (ncfixtyp(i).eq.iopt_radius) then
        rada(ncfixind(i)) = xcon(i)
      elseif (ncfixtyp(i).eq.iopt_xf) then
        xafrac(ncfixind(i)) = xcon(i)
      elseif (ncfixtyp(i).eq.iopt_yf) then
        yafrac(ncfixind(i)) = xcon(i)
      elseif (ncfixtyp(i).eq.iopt_zf) then
        zafrac(ncfixind(i)) = xcon(i)
      elseif (ncfixtyp(i).eq.iopt_xcom) then
        molcoma(1,ncfixind(i)) = xcon(i)
      elseif (ncfixtyp(i).eq.iopt_ycom) then
        molcoma(2,ncfixind(i)) = xcon(i)
      elseif (ncfixtyp(i).eq.iopt_zcom) then
        molcoma(3,ncfixind(i)) = xcon(i)
      elseif (ncfixtyp(i).eq.iopt_xqtn) then
        molQa(1,ncfixind(i)) = xcon(i)
      elseif (ncfixtyp(i).eq.iopt_yqtn) then
        molQa(2,ncfixind(i)) = xcon(i)
      elseif (ncfixtyp(i).eq.iopt_zqtn) then
        molQa(3,ncfixind(i)) = xcon(i)
      endif
    enddo
!
    deallocate(xcon,stat=status)
    if (status/=0) call deallocate_error('applyconstraints','xcon')
!
!  Handle additive constraints for fractional coordinates - take nearest pair of images
!
    if (ndim.gt.0) then
      allocate(ncount(3,nasym),stat=status)
      if (status/=0) call outofmemory('applyconstraints','ncount')
!
      ncount(1:3,1:nasym) = 0
!
      do i = 1,ncon
        ii = ncfixind(i)
        if (ncfixtyp(i).eq.iopt_xf) then
          ncount(1,ii) = ncount(1,ii) + 1
        elseif (ncfixtyp(i).eq.iopt_yf) then
          ncount(2,ii) = ncount(2,ii) + 1
        elseif (ncfixtyp(i).eq.iopt_zf) then
          ncount(3,ii) = ncount(3,ii) + 1
        endif
      enddo
      do i = 1,nasym
        do ix = 1,3
          if (ix.eq.1) then
            iopt_ix = iopt_xf
          elseif (ix.eq.2) then
            iopt_ix = iopt_yf
          elseif (ix.eq.3) then
            iopt_ix = iopt_zf
          endif
!
!  Select only those coordinates which are fractional
!
          if (ndim.eq.3) then
            ldothis = .true.
          elseif (ndim.eq.2) then
            ldothis = (ix.le.2)
          elseif (ndim.eq.1) then
            ldothis = (ix.eq.1)
          endif
          if (ncount(ix,i).ge.2.and.ldothis) then
            lfound = .false.
            j = 0
            do while (.not.lfound.and.j.lt.ncon-1)
              j = j + 1
              if (ncfixind(j).eq.i.and.ncfixtyp(j).eq.iopt_ix) then
                k = j
                do while (.not.lfound.and.k.lt.ncon) 
                  k = k + 1
                  lfound = (ncfixind(k).eq.i.and.ncfixtyp(k).eq.iopt_ix)
                enddo
              endif
            enddo
            if (lfound) then
              if (ncvartyp(j).eq.iopt_xf) then
                xj = xafrac(ncvarind(j))
              elseif (ncvartyp(j).eq.iopt_yf) then
                xj = yafrac(ncvarind(j))
              elseif (ncvartyp(j).eq.iopt_zf) then
                xj = zafrac(ncvarind(j))
              endif
              if (ncvartyp(k).eq.iopt_xf) then
                xk = xafrac(ncvarind(k))
              elseif (ncvartyp(k).eq.iopt_yf) then
                xk = yafrac(ncvarind(k))
              elseif (ncvartyp(k).eq.iopt_zf) then
                xk = zafrac(ncvarind(k))
              endif
              diff = abs(xj-xk)
              if (diff.gt.0.5_dp) then
                if (iopt_ix.eq.iopt_xf) then
                  xafrac(i) = xafrac(i) + 0.5_dp
                  xafrac(i) = mod(xafrac(i),1.0_dp)
                elseif (iopt_ix.eq.iopt_yf) then
                  yafrac(i) = yafrac(i) + 0.5_dp
                  yafrac(i) = mod(yafrac(i),1.0_dp)
                elseif (iopt_ix.eq.iopt_zf) then
                  zafrac(i) = zafrac(i) + 0.5_dp
                  zafrac(i) = mod(zafrac(i),1.0_dp)
                endif
              endif
            endif
          endif
        enddo
      enddo
!
      deallocate(ncount,stat=status)
      if (status/=0) call deallocate_error('applyconstraints','ncount')
    endif
  endif
#ifdef TRACE
  call trace_out('applyconstraints')
#endif
!
  return
  end
