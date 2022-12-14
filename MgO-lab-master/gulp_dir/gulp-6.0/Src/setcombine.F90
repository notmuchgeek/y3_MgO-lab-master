  subroutine setcombine
!
!  Uses combination rules and epsilon and sigma values to
!  set up lennard jones potential coefficients.
!
!   4/98 ESFF Lennard-Jones combination rules added
!   3/04 Geometric combination rules added
!   1/09 Arrays apot, bpot, cpot & dpot replaced by new 2-D array twopot
!   7/13 Product combination rule added
!   2/15 Modified to allow for lmm3se flag
!   2/15 Routine renamed from setlenn to setcombine to reflect increased generality
!   2/15 potcut called to ensure that any energy/gradient shifts are updated.
!   2/15 Combination rules for MM3buck added
!   3/15 Modifications for bond type specific combination of epsilon/sigma added
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
!  Copyright Curtin University 2015
!
!  Julian Gale, CIC, Curtin University, March 2015
!
  use two
  implicit none
!
!  Local variables
!
  integer(i4)           :: i
  integer(i4)           :: j
  integer(i4)           :: nat1
  integer(i4)           :: nat2
  integer(i4)           :: ntyp1
  integer(i4)           :: ntyp2
  logical               :: lfound
  logical               :: lfound1
  logical               :: lfound2
  real(dp)              :: a1
  real(dp)              :: a2
  real(dp)              :: b1
  real(dp)              :: b2
  real(dp)              :: eps1
  real(dp)              :: eps2
  real(dp)              :: sig1
  real(dp)              :: sig13
  real(dp)              :: sig2
  real(dp)              :: sig23
  real(dp)              :: trm1
!
  do i = 1,npote
    if (lcombine(i)) then
      if (nptype(i).eq.2) then
!
!  A/B combination rules
!
        nat1 = nspec1(i)
        nat2 = nspec2(i)
        ntyp1 = nptyp1(i)
        ntyp2 = nptyp2(i)
        a1 = 0.0_dp
        a2 = 0.0_dp
        b1 = 0.0_dp
        b2 = 0.0_dp
        lfound = .false.
        j = 0
        do while (j.lt.natab.and..not.lfound)
          j = j + 1
          if (nat1.eq.nattab(j).and.(ntyp1.eq.ntypab(j).or.ntypab(j).eq.0)) then
            lfound = .true.
            a1 = atoma(j)
            b1 = atomb(j)
          endif
        enddo
        lfound = .false.
        j = 0
        do while (j.lt.natab.and..not.lfound)
          j = j + 1
          if (nat2.eq.nattab(j).and.(ntyp2.eq.ntypab(j).or.ntypab(j).eq.0)) then
            lfound = .true.
            a2 = atoma(j)
            b2 = atomb(j)
          endif
        enddo
        twopot(1,i) = sqrt(a1*a2)
        twopot(2,i) = sqrt(b1*b2)
      elseif (nptype(i).eq.21) then
!
!  ESFF combination rules
!
        nat1 = nspec1(i)
        nat2 = nspec2(i)
        ntyp1 = nptyp1(i)
        ntyp2 = nptyp2(i)
        sig1 = 0.0_dp
        sig2 = 0.0_dp
        eps1 = 0.0_dp
        eps2 = 0.0_dp
        lfound = .false.
        j = 0
        sig1 = 0.0_dp
        eps1 = 0.0_dp
        do while (j.lt.nseps.and..not.lfound)
          j = j + 1
          if (.not.lmm3se(j)) then
            if (nat1.eq.natse(j).and.(ntyp1.eq.ntypse(j).or.ntypse(j).eq.0)) then
              lfound = .true.
              sig1 = sigma(j)
              eps1 = epsilon(j)
            endif
          endif
        enddo
        lfound = .false.
        j = 0
        sig2 = 0.0_dp
        eps2 = 0.0_dp
        do while (j.lt.nseps.and..not.lfound)
          j = j + 1
          if (.not.lmm3se(j)) then
            if (nat2.eq.natse(j).and.(ntyp2.eq.ntypse(j).or.ntypse(j).eq.0)) then
              lfound = .true.
              sig2 = sigma(j)
              eps2 = epsilon(j)
            endif
          endif
        enddo
        eps1 = sqrt(eps1)
        eps2 = sqrt(eps2)
        sig1 = sig1*sig1*sig1
        sig2 = sig2*sig2*sig2
        a1 = eps1*sig1*sig1
        b1 = eps1*sig1
        a2 = eps2*sig2*sig2
        b2 = eps2*sig2
        twopot(1,i) = (a1*b2 + a2*b1)
        twopot(2,i) = 3.0_dp*b1*b2
      elseif (nptype(i).eq.57) then
!
!  MM3 buckingham combination rules
!
        nat1 = nspec1(i)
        nat2 = nspec2(i)
        ntyp1 = nptyp1(i)
        ntyp2 = nptyp2(i)
        sig1 = 0.0_dp
        sig2 = 0.0_dp
        eps1 = 0.0_dp
        eps2 = 0.0_dp
        lfound = .false.
        j = 0
        sig1 = 0.0_dp
        eps1 = 0.0_dp
        do while (j.lt.nseps.and..not.lfound) 
          j = j + 1
          if (lmm3se(j)) then
            if (nat1.eq.natse(j).and.(ntyp1.eq.ntypse(j).or.ntypse(j).eq.0)) then
              lfound = .true.
              sig1 = sigma(j)
              eps1 = epsilon(j)
            endif
          endif
        enddo
        lfound = .false.
        j = 0
        sig2 = 0.0_dp
        eps2 = 0.0_dp
        do while (j.lt.nseps.and..not.lfound)
          j = j + 1
          if (lmm3se(j)) then
            if (nat2.eq.natse(j).and.(ntyp2.eq.ntypse(j).or.ntypse(j).eq.0)) then
              lfound = .true.
              sig2 = sigma(j)
              eps2 = epsilon(j)
            endif
          endif
        enddo
!
!  Combine and assign to epsilon and rv
!
        trm1 = sig1 + sig2
        twopot(4,i) = sqrt(eps1*eps2)
        twopot(5,i) = 0.5_dp*trm1
      else
!
!  Sigma and epsilon combination rules
!
        nat1 = nspec1(i)
        nat2 = nspec2(i)
        ntyp1 = nptyp1(i)
        ntyp2 = nptyp2(i)
        sig1 = 0.0_dp
        sig2 = 0.0_dp
        eps1 = 0.0_dp
        eps2 = 0.0_dp
        lfound1 = .false.
        j = 0
        sig1 = 0.0_dp
        eps1 = 0.0_dp
        do while (j.lt.nseps.and..not.lfound1) 
          j = j + 1
          if (.not.lmm3se(j)) then
            if (nat1.eq.natse(j).and.(ntyp1.eq.ntypse(j).or.ntypse(j).eq.0)) then
              if (mmexcse(j).eq.mmexc(i).or.mmexcse(j).eq.0) then
                lfound1 = .true.
                sig1 = sigma(j)
                eps1 = epsilon(j)
              endif
            endif
          endif
        enddo
        lfound2 = .false.
        j = 0
        sig2 = 0.0_dp
        eps2 = 0.0_dp
        do while (j.lt.nseps.and..not.lfound2)
          j = j + 1
          if (.not.lmm3se(j)) then
            if (nat2.eq.natse(j).and.(ntyp2.eq.ntypse(j).or.ntypse(j).eq.0)) then
              if (mmexcse(j).eq.mmexc(i).or.mmexcse(j).eq.0) then
                lfound2 = .true.
                sig2 = sigma(j)
                eps2 = epsilon(j)
              endif
            endif
          endif
        enddo
        if (lfound1.and.lfound2) then
!
!  Combine and assign to apot and bpot
!
          if (ncombipower(i).eq.6) then
            sig13 = sig1**3
            sig23 = sig2**3
            trm1 = sig13*sig13 + sig23*sig23
            twopot(1,i) = 2.0_dp*sqrt(eps1*eps2)*sig13*sig23/trm1
            twopot(2,i) = (0.5_dp*trm1)**(1.0_dp/6.0_dp)
          elseif (ncombipower(i).eq.2) then
            trm1 = sig1 + sig2
            twopot(1,i) = sqrt(eps1*eps2)
            twopot(2,i) = 0.5_dp*trm1
          elseif (ncombipower(i).eq.1) then
            twopot(1,i) = sqrt(eps1*eps2)
            twopot(2,i) = sqrt(sig1*sig2)
          endif
        else
!
!  One or more species were not found so zero parameters
!
          twopot(1,i) = 0.0_dp
          twopot(2,i) = 0.0_dp
        endif
      endif
    endif
  enddo
!
!  Correct energy/gradient shifts
!
  do i = 1,npote
    if (leshift(i)) call potcut(i)
  enddo
!
  return
  end
