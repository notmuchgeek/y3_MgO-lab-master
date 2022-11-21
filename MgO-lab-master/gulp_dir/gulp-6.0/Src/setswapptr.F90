  subroutine setswapptr(nswap2try)
!
!  Sets up pointer to atoms that can be swapped
!
!   1/09 Created from setmoveptr
!   5/16 Modified for multiple swap pairs
!   4/19 Call to lmatchany modified to include wildcard argument
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
!  Julian Gale, CIC, Curtin University, April 2019
!
  use current
  use montecarlo
  use molecule,     only : natmol
  use reallocate
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)  :: nswap2try
!
!  Local variables
!
  integer(i4)              :: i, ierror
  integer(i4)              :: status
  logical,            save :: lfirstime = .true.
  logical                  :: lmatchany
  logical                  :: lmatched
!
!  Initialise memory on first call
!
  if (lfirstime) then
    maxswapable = numat
    allocate(nptrswapable(maxswapable,maxmcswaps),stat=status)
    if (status/=0) call outofmemory('setswapptr','nptrswapable')
    lfirstime = .false.
  endif
!
!  Setup pointer to atoms that can be swapped 
!
  nswapable(nswap2try) = 0
  do i = 1,numat
!
!  Exclude atoms in molecules
!
    if (natmol(i).eq.0) then
!
      if (lmcswapany(nswap2try)) then
        nswapable(nswap2try) = nswapable(nswap2try) + 1
        if (nswapable(nswap2try).gt.maxswapable) then
          maxswapable = nswapable(nswap2try) + 10
          call realloc(nptrswapable,maxswapable,maxmcswaps,ierror)
          if (ierror.ne.0) call outofmemory('setswapptr','nptrswapable')
        endif
        nptrswapable(nswapable(nswap2try),nswap2try) = i
      else
        lmatched = lmatchany(nat(i),nftype(i),nmcswapspec(nswap2try), &
                     nmcswapnat(1,nswap2try),nmcswaptype(1,nswap2try),.false.)
        if (lmatched) then
          nswapable(nswap2try) = nswapable(nswap2try) + 1
          if (nswapable(nswap2try).gt.maxswapable) then
            maxswapable = nswapable(nswap2try) + 10
            call realloc(nptrswapable,maxswapable,maxmcswaps,ierror)
            if (ierror.ne.0) call outofmemory('setswapptr','nptrswapable')
          endif
          nptrswapable(nswapable(nswap2try),nswap2try) = i
        endif
      endif
    endif
  enddo
!
  return
  end
