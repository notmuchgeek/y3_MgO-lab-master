  subroutine init2bodydefaults(npot)
!
!  Initialises the default values of variables for 2 body potentials
!
!  11/06 Created from changemaxpot
!   7/07 n2botype(2,*) defaults to 0
!   1/09 Arrays apot, bpot, cpot & dpot replaced by new 2-D array twopot
!   1/09 leshift/lgshift logicals introduced to replace use of tpot(3/4,) for energy
!        and gradient shifts
!   4/09 Size of first dimension of twopot array increased to 6
!   4/09 Size of first dimension of twopot array increased to 7
!   3/10 Arrays that flag rule generated potentials added
!   9/10 Checks added that npot is between 1 and maxpot
!   8/14 Size of first dimension of twopot array increased to 12
!   8/14 ipot added
!   8/14 lorder12 added
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
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  use two
  implicit none
!
  integer(i4), intent(in) :: npot
!
  integer(i4)             :: j
#ifdef TRACE
  call trace_in('init2bodydefaults')
#endif
!
!  Check that npot passed in is sensible
!
  if (npot.le.maxpot.and.npot.ge.1) then
!
!  Initialise defaults 
!
    lcombine(npot) = .false.
    lgenerated2(npot) = .false.
    leshift(npot) = .false.
    lgshift(npot) = .false.
    lorder12(npot) = .true.
    ncombipower(npot) = 0
    ipot(1,npot) = 0
    ipot(2,npot) = 1
    ipot(3,npot) = 0
    ipot(4,npot) = 0
    mmexc(npot) = 0
    n2botype(1,npot) = 0
    n2botype(2,npot) = 0
    twopot(1:12,npot) = 0.0_dp
    rpot(npot) = 0.0_dp
    rpot2(npot) = 0.0_dp
    scale14(npot) = 1.0_dp
    taperpot(npot) = 0.0_dp
    tapergrad(npot) = 0.0_dp
    eshift(npot) = 0.0_dp
    gshift(npot) = 0.0_dp
    do j = 1,12
      tpot(j,npot) = 0.0_dp
    enddo
  endif
#ifdef TRACE
  call trace_out('init2bodydefaults')
#endif
!
  return
  end
