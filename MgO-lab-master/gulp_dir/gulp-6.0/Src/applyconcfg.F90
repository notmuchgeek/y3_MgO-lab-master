  subroutine applyconcfg
!
!  Apply constraints to configuration arrays to re-symmetrise.
!  Currently needed for simul fitting.
!
!  12/17 Created from setup
!   1/18 Trace added
!   3/19 Constraint arrays changed to have index and type
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
!  Julian Gale, CIC, Curtin University, March 2019
!
  use control
  use configurations
  use current
  use datatypes
  use optimisation
  use reallocate
  use shells
  use species
  use symmetry
#ifdef TRACE
  use trace,      only : trace_in, trace_out
#endif
  implicit none
!
!  Local variables
!
  integer(i4), dimension(:,:), allocatable     :: ncount
  integer(i4)                                  :: i
  integer(i4)                                  :: ii
  integer(i4)                                  :: iopt_ix
  integer(i4)                                  :: ix
  integer(i4)                                  :: j
  integer(i4)                                  :: k
  integer(i4)                                  :: nf
  integer(i4)                                  :: nv
  integer(i4)                                  :: nvj
  integer(i4)                                  :: nvk
  integer(i4)                                  :: status
  logical                                      :: lfound
  real(dp)                                     :: diff
  real(dp)                                     :: vcrd
  real(dp)                                     :: vcrdj
  real(dp)                                     :: vcrdk
#ifdef TRACE
  call trace_in('applyconcfg')
#endif
!
!  Find first atom shift
!
  nsft = 0
  if (ncf.gt.1) then
    do i = 1,ncf-1
      nsft = nsft + nascfg(i)
    enddo
  endif
!
!  Set dimensionality
!
  ndim = ndimen(ncf)
!
!  Set number of strains according to the dimensionality and pointer
!
  if (ndim.eq.3) then
    nstrains = 6
  elseif (ndim.eq.2) then
    nstrains = 3
  elseif (ndim.eq.1) then
    nstrains = 1
  else
    nstrains = 0
  endif
!
!  Transfer stored configuration data into working arrays
!
  nasym = nascfg(ncf)
!
!  Set up constraint pointers
!
  if (ncf.eq.ncfg) then
    if (ncontot.gt.0) ncon = ncontot + 1 - n1con(ncf)
  else
    ncon = n1con(ncf+1) - n1con(ncf)
  endif
  if (ncon.gt.maxcon) then
    maxcon = ncon
    call changemaxcon
  endif
  ncfst = n1con(ncf) - 1
!
!  Apply constraints
!
  if (ncon.gt.0) then
    do j = 1,ncon
      ncfixind(j) = ncfixindcfg(ncfst+j)
      ncvarind(j) = ncvarindcfg(ncfst+j)
      conco(j) = concocfg(ncfst+j)
      conadd(j) = conaddcfg(ncfst+j)
    enddo
    do j = 1,ncon
      nf = ncfixind(j)
      if (ncfixtyp(j).eq.iopt_xf) then
        xcfg(nsft+nf) = 0.0_dp
      elseif (ncfixtyp(j).eq.iopt_yf) then
        ycfg(nsft+nf) = 0.0_dp
      elseif (ncfixtyp(j).eq.iopt_zf) then
        zcfg(nsft+nf) = 0.0_dp
      endif
    enddo
    do j = 1,ncon
      nv = ncvarind(j)
      if (ncvartyp(j).eq.iopt_xf) then
        vcrd = xcfg(nsft+nv)
      elseif (ncvartyp(j).eq.iopt_yf) then
        vcrd = ycfg(nsft+nv)
      elseif (ncvartyp(j).eq.iopt_zf) then
        vcrd = zcfg(nsft+nv)
      endif
      nf = ncfixind(j)
      if (ncfixtyp(j).eq.iopt_xf) then
        xcfg(nsft+nf) = vcrd*conco(j) + conadd(j) + xcfg(nsft+nf)
      elseif (ncfixtyp(j).eq.iopt_yf) then
        ycfg(nsft+nf) = vcrd*conco(j) + conadd(j) + ycfg(nsft+nf)
      elseif (ncfixtyp(j).eq.iopt_zf) then
        zcfg(nsft+nf) = vcrd*conco(j) + conadd(j) + zcfg(nsft+nf)
      endif
    enddo
!
!  Handle additive constraints for fractional coordinates - take nearest pair of images
!
    if (ndim.eq.3) then
      allocate(ncount(3,nasym),stat=status)
      if (status/=0) call outofmemory('applyconcfg','ncount')
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
!
      do i = 1,nasym
        do ix = 1,3
          if (ncount(ix,i).ge.2) then
            if (ix.eq.1) then
              iopt_ix = iopt_xf
            elseif (ix.eq.2) then
              iopt_ix = iopt_yf
            elseif (ix.eq.3) then
              iopt_ix = iopt_zf
            endif
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
              nvj = ncvarind(j)
              if (ncvartyp(j).eq.iopt_xf) then
                vcrdj = xcfg(nsft+nvj)
              elseif (ncvartyp(j).eq.iopt_yf) then
                vcrdj = ycfg(nsft+nvj)
              elseif (ncvartyp(j).eq.iopt_zf) then
                vcrdj = zcfg(nsft+nvj)
              endif
              nvk = ncvarind(k)
              if (ncvartyp(k).eq.iopt_xf) then
                vcrdk = xcfg(nsft+nvk)
              elseif (ncvartyp(k).eq.iopt_yf) then
                vcrdk = ycfg(nsft+nvk)
              elseif (ncvartyp(k).eq.iopt_zf) then
                vcrdk = zcfg(nsft+nvk)
              endif
              diff = abs(vcrdk - vcrdj)
              if (diff.gt.0.5_dp) then
                nf = ncfixind(j)
                if (ncfixtyp(j).eq.iopt_xf) then
                  xcfg(nsft+nf) = xcfg(nsft+nf) + 0.5_dp
                  xcfg(nsft+nf) = mod(xcfg(nsft+nf),1.0_dp)
                elseif (ncfixtyp(j).eq.iopt_yf) then
                  ycfg(nsft+nf) = ycfg(nsft+nf) + 0.5_dp
                  ycfg(nsft+nf) = mod(ycfg(nsft+nf),1.0_dp)
                elseif (ncfixtyp(j).eq.iopt_zf) then
                  zcfg(nsft+nf) = zcfg(nsft+nf) + 0.5_dp
                  zcfg(nsft+nf) = mod(zcfg(nsft+nf),1.0_dp)
                endif
              endif
            endif
          endif
        enddo
      enddo
      deallocate(ncount,stat=status)
      if (status/=0) call deallocate_error('applyconcfg','ncount')
    endif
  endif
#ifdef TRACE
  call trace_out('applyconcfg')
#endif
!
  return
  end
