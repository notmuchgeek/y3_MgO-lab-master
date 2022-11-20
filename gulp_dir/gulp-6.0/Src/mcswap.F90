  subroutine mcswap(mode,nswap2try)
!
!  MC routine for swapping of atoms. Approach taken is to swap
!  coordinates rather than attributes since this is simpler.
!
!  mode = if mode = 1, choose atoms to apply swap to
!         if mode = 2, then create new trial swap
!         if mode = 3, then undo previous swap
!
!   1/09 Created from mcmove
!   4/16 Variable number of pairs to swap added
!   5/16 Modified for multiple possible swaps
!   5/16 nswap2try added to arguments
!   5/16 Corrections to swap algorithm pointers
!   2/18 Trace added
!   2/19 x0 removed
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
!  Julian Gale, CIC, Curtin University, February 2019
!
  use current
  use general
  use genetic,       only : iseed
  use montecarlo
  use parallel
  use reallocate
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
!
!  Passed variables
!
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                  :: mode
  integer(i4), intent(in)                  :: nswap2try    ! Which if the nmcswaps to use
!
!  Local variables
!
  integer(i4), dimension(:), pointer, save :: nptrswap1 => null()
  integer(i4), dimension(:), pointer, save :: nptrswap2 => null()
  integer(i4),                        save :: nswap1
  integer(i4),                        save :: nswap2
  integer(i4),                        save :: nswap2l
!
  integer(i4), dimension(:), pointer, save :: nptrswapabletot => null()
  integer(i4), dimension(:), pointer, save :: nptrswapable2 => null()
  integer(i4)                              :: i
  integer(i4)                              :: ii
  integer(i4)                              :: nat1
  integer(i4)                              :: ns
  integer(i4)                              :: nswapable2
  integer(i4)                              :: nswapabletot
  integer(i4)                              :: ntype1
  integer(i4)                              :: status
  logical,                            save :: lfirstcall = .true.
  logical                                  :: lmatch
  real(dp)                                 :: randnum
  real(dp)                                 :: GULP_random
  real(dp)                                 :: x1
  real(dp)                                 :: y1
  real(dp)                                 :: z1
  real(dp)                                 :: x2
  real(dp)                                 :: y2
  real(dp)                                 :: z2
#ifdef TRACE
  call trace_in('mcswap')
#endif
!
!  Initialisation of swap pointer arrays
!
  if (lfirstcall) then
    lfirstcall = .false.
    allocate(nptrswap1(nmcswappair(nswap2try)),stat=status)
    if (status/=0) call outofmemory('mcswap','nptrswap1')
    allocate(nptrswap2(nmcswappair(nswap2try)),stat=status)
    if (status/=0) call outofmemory('mcswap','nptrswap2')
  endif
!
  if (mode.eq.3) then
!****************************
!  Mode 3 : Undo last swap  *
!****************************
    do ns = 1,nmcswappair(nswap2try)
      x1 = xafrac(nptrswap1(ns))
      y1 = yafrac(nptrswap1(ns))
      z1 = zafrac(nptrswap1(ns))  
      x2 = xafrac(nptrswap2(ns))
      y2 = yafrac(nptrswap2(ns))
      z2 = zafrac(nptrswap2(ns))  
!
      xafrac(nptrswap1(ns)) = x2
      yafrac(nptrswap1(ns)) = y2
      zafrac(nptrswap1(ns)) = z2
      xafrac(nptrswap2(ns)) = x1
      yafrac(nptrswap2(ns)) = y1
      zafrac(nptrswap2(ns)) = z1
    enddo
  elseif (mode.eq.1) then
!**********************
!  Mode 1 : New swap  *
!**********************
!
!  Allocate array to store details of remaining swapable atoms
!
    allocate(nptrswapable2(nswapable(nswap2try)),stat=status)
    if (status/=0) call outofmemory('mcswap','nptrswapable2')
    allocate(nptrswapabletot(nswapable(nswap2try)),stat=status)
    if (status/=0) call outofmemory('mcswap','nptrswapabletot')
!
!  Build local list of atoms to swap from which atoms will be removed if already swapped
!
    nswapabletot = nswapable(nswap2try)
    do i = 1,nswapable(nswap2try)
      nptrswapabletot(i) = i
    enddo
!
!  Loop over number of pairs to swap
!
    ntrialatom = 0
    do ns = 1,nmcswappair(nswap2try)
!
!  Choose first atom to swap
!
      randnum = GULP_random(iseed,1_i4)
      nswap1 = nswapabletot*randnum + 1_i4
      if (nswap1.gt.nswapabletot) nswap1 = nswapabletot
!
!  Translate first atom from remaining swapable list to overall pointer
!
      nptrswap1(ns) = nptrswapable(nptrswapabletot(nswap1),nswap2try)
!
      nat1 = nat(nptrswap1(ns))
      ntype1 = nftype(nptrswap1(ns))
!
!  Remove atom 1 from local swapable list
!
      do i = nswap1+1,nswapable(nswap2try)
        nptrswapabletot(i-1) = nptrswapabletot(i)
      enddo
      nswapabletot = nswapabletot - 1
!
!  Build list of valid atom 2s to swap with atom 1
!
      nswapable2 = 0
      do i = 1,nswapabletot
        ii = nptrswapable(nptrswapabletot(i),nswap2try)
        if (.not.lmatch(nat1,ntype1,nat(ii),nftype(ii),.false.)) then
          nswapable2 = nswapable2 + 1
          nptrswapable2(nswapable2) = i
        endif
      enddo
!
!  Check that there are some swaps possible
!
      if (nswapable2.eq.0) then
        call outerror('swap requested but no meaningful swap possible',0_i4)
        call stopnow('mcswap')
      endif
!
!  Choose second atom to swap
!
      randnum = GULP_random(iseed,1_i4)
      nswap2l = nswapable2*randnum + 1_i4
      if (nswap2l.gt.nswapable2) nswap2l = nswapable2
!
!  Translate second atom from remaining swapable list to overall pointer
!
      nswap2 = nptrswapable2(nswap2l)
      nptrswap2(ns) = nptrswapable(nptrswapabletot(nswap2),nswap2try)
!
!  Remove atom 2 from local swapable list
!
      do i = nswap2+1,nswapable(nswap2try)
        nptrswapabletot(i-1) = nptrswapabletot(i)
      enddo
      nswapabletot = nswapabletot - 1
!
!  Copy pointers to main arrays
!
      nptrtrialatom(ntrialatom+1) = nptrswap1(ns)
      nptrtrialatom(ntrialatom+2) = nptrswap2(ns)
      ntrialatom = ntrialatom + 2
    enddo
!
!  Deallocate array used to store details of second atom choice
!
    deallocate(nptrswapabletot,stat=status)
    if (status/=0) call deallocate_error('mcswap','nptrswapabletot')
    deallocate(nptrswapable2,stat=status)
    if (status/=0) call deallocate_error('mcswap','nptrswapable2')
  elseif (mode.eq.2) then
!************************
!  Mode 2 : Apply swap  *
!************************
!
!  Apply swap to configuration array
!
    do ns = 1,nmcswappair(nswap2try)
      x1 = xafrac(nptrswap1(ns))
      y1 = yafrac(nptrswap1(ns))
      z1 = zafrac(nptrswap1(ns))  
      x2 = xafrac(nptrswap2(ns))
      y2 = yafrac(nptrswap2(ns))
      z2 = zafrac(nptrswap2(ns))  
!
      xafrac(nptrswap1(ns)) = x2
      yafrac(nptrswap1(ns)) = y2
      zafrac(nptrswap1(ns)) = z2
      xafrac(nptrswap2(ns)) = x1
      yafrac(nptrswap2(ns)) = y1
      zafrac(nptrswap2(ns)) = z1
    enddo
  endif
#ifdef TRACE
  call trace_out('mcswap')
#endif
!
  return
  end
