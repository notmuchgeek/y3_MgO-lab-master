  subroutine kinthewedge(kpt)
!
!  Finds the symmetry related image for a k point in the asymmetric wedge.
!
!  11/14 Created from setkpt
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
  use control
  use current
  use ksample
  use symmetry
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  real(dp),    intent(inout) :: kpt(3)
!
!  Local variables
!
  integer(i4)                :: iprimgp(232)
  integer(i4)                :: nkimage
  integer(i4)                :: npgrp
  integer(i4)                :: nwedge
  logical                    :: lfound
  real(dp)                   :: kimage(3,48)
!
!  Primitive
!
  data iprimgp/2*1, &
               2*2,1,2*2,2*1,2*2,1,2*2,1, & !  Monoclinic
               4*3,2*2,3*1,10*3,7*2,5*1,16*3,6*2,6*1, & !  Orthorhombic
               4*4,2*1,4,1,4*4,2*1,8*5,2*1,8*5,4*1,8*5,4*1,16*5,4*1, & !  Tetragonal
               3*6,13,6,13,8,7,8,7,8,7,14,7,8,7,8,2*14,2*8,2*7,2*14, & !  Trigonal
               9*9,18*10, & !  Hexagonal
               11,2*13,11,13,2*11,3*13,11,13,2*12,3*14,2*12,14,12, & !  Cubic
               2*14,12,2*14,4*12,6*14, & 
               2*1/ !  Additional groups!
!
!  Check that this is 3-D
!
  if (ndim.ne.3) return
#ifdef TRACE
  call trace_in('kinthewedge')
#endif
!
  if (nspcg(ncf).gt.0) then
    npgrp = iprimgp(nspcg(ncf))
  else
    npgrp = 1
  endif
!
!  Generate all images of this k point
!
  call symgenkpt(kpt,48_i4,nkimage,kimage)
!
!  Loop over images looking for the one that is in the wedge
!
  lfound = .false.
  nwedge = 0
  do while (nwedge.lt.nkimage.and..not.lfound)
    nwedge = nwedge + 1
    kpt(1) = kimage(1,nwedge)
    kpt(2) = kimage(2,nwedge)
    kpt(3) = kimage(3,nwedge)
!
    if (npgrp.eq.12.or.npgrp.eq.13.or.npgrp.eq.14) then
      if (npgrp.eq.12) then
        if (kpt(1).le.0.5_dp.and.kpt(2).le.0.5_dp.and.kpt(3).le.0.5_dp) then
          lfound = (kpt(1).ge.kpt(2).and.kpt(2).ge.kpt(3)) 
        else
          lfound = .false.
        endif
      else
        lfound = (kpt(1).ge.kpt(2).and.kpt(2).ge.kpt(3)) 
      endif
    elseif (npgrp.eq.11) then
      if (kpt(1).le.0.5_dp.and.kpt(2).le.0.5_dp.and.kpt(3).le.0.5_dp) then
        if (kpt(2).le.kpt(1)) then
          lfound = (kpt(2).ge.kpt(3)) 
        else
          lfound = (kpt(1).ge.kpt(3)) 
        endif
      else
        lfound = .false.
      endif
    elseif (npgrp.eq.10) then
      if (kpt(1).le.0.5_dp) then
        lfound = (kpt(1).ge.kpt(2).and.kpt(3).le.0.5_dp)
      else
        lfound = (kpt(2).le.(1.0_dp - kpt(1)).and.kpt(3).le.0.5_dp)
      endif
    elseif (npgrp.eq.9) then
      lfound = (kpt(1).le.0.5_dp.and.kpt(3).le.0.5_dp)
    elseif (npgrp.eq.8) then
      if (kpt(1).le.0.5_dp) then
        lfound = (kpt(1).ge.kpt(2))
      else
        if (kpt(2).eq.(1.0_dp - kpt(1))) then
          lfound = (kpt(3).le.0.5_dp)
        else
          lfound = (kpt(2).le.(1.0_dp - kpt(1)))
        endif
      endif
    elseif (npgrp.eq.7) then
      if (kpt(1).le.0.5_dp) then
        if (kpt(1).eq.kpt(2)) then
          lfound = (kpt(3).le.0.5_dp)
        else
          lfound = (kpt(1).ge.kpt(2))
        endif
      else
        lfound = (kpt(2).le.(1.0_dp - kpt(1)))
      endif
    elseif (npgrp.eq.6) then
      lfound = (kpt(1).le.0.5_dp)
    elseif (npgrp.eq.5) then
      if (kpt(1).le.0.5_dp.and.kpt(2).le.0.5_dp.and.kpt(3).le.0.5_dp) then
        lfound = (kpt(2).le.kpt(1))
      else
        lfound = .false. 
      endif
    elseif (npgrp.eq.4) then
      lfound = (kpt(1).le.0.5_dp.and.kpt(2).le.0.5_dp.and.kpt(3).le.0.5_dp) 
    elseif (npgrp.eq.3) then
      lfound = (kpt(1).le.0.5_dp.and.kpt(2).le.0.5_dp.and.kpt(3).le.0.5_dp)
    elseif (npgrp.eq.2) then
      lfound = (kpt(1).le.0.5_dp.and.kpt(2).le.0.5_dp)
    elseif (npgrp.eq.1) then
      lfound = (kpt(1).le.0.5_dp)
    endif
  enddo
!
!  If no image was found then stop with an error
!
  if (.not.lfound) then
    call outerror('no k point image found in the asymmetric wedge',0_i4)
    call stopnow('kinthewedge')
  endif
#ifdef TRACE
  call trace_out('kinthewedge')
#endif
!
  return
  end
