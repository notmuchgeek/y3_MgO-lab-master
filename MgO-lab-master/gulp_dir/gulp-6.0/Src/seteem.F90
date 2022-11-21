  subroutine seteem
!
!  Set up electronegativity parameters if required
!
!   5/18 Created
!   3/20 angstoev changed to inverse_angstroms_to_ev
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
!  Julian Gale, CIC, Curtin University, March 2020
!
  use control,       only : leem, keyword
  use element
  use eemdata
  use g_constants,   only : inverse_angstroms_to_ev
  use iochannels
  use species
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Local variables
!
  integer(i4)      :: i
  integer(i4)      :: j
  integer(i4)      :: iele
  logical          :: lfound
  logical          :: lglobal
!
!  If no EEM then return
!
  if (.not.leem) return
#ifdef TRACE
  call trace_in('seteem')
#endif
!
!  Substitute parameters for oldeem keyword:
!
!  Note that there are parameter sets for hydrogen:
!     nat = 1 => H+
!     nat = 2 => H-
!
  if (.not.lqeq.and..not.lSandM.and.index(keyword,'oldeem').ne.0) then
    chi(14) = 3.478_dp
    rmu(14) = 6.408_dp
    rmu(8) = 9.466_dp
    chi(1) = 3.398_dp
    rmu(1) = 20.818_dp
    chi(2) = 4.706_dp
    rmu(2) = 8.899_dp
  endif
!
!  Find type of EEM for this run
!
  if (lqeq) then
    neemtype = 2
  elseif (lSandM) then
    neemtype = 3
  elseif (lpacha) then
    neemtype = 4
  elseif (leem) then
    neemtype = 1
  endif
!
!  Initialise lelementOK
!
  lelementOK(1:maxele) = .false.
!
!  Loop over elements looking for those that are present
!
  do iele = 1,maxele
    lfound = any(natspec(1:nspec).eq.iele)
    if (lfound.and.nqrange(iele,neemtype).eq.0) then
!##########################################################################################
!  Parameters needed but none specified by the user and so apply defaults where possible  #
!##########################################################################################
      nqrange(iele,neemtype) = 1
      if (nqrange(iele,neemtype).gt.maxqrange) then
        maxqrange = 1
        call changemaxqrange
      endif
!
!  Default to the same parameters for all charge values
!
      nqrangetype(1,iele,neemtype) = 0
!
!  Set method specific default parameters
!
      if (neemtype.eq.2) then
!
!  QEq
!
        chirange(1,iele,neemtype) = qeqchi(iele)
        murange(1,iele,neemtype) = qeqmu(iele)
        q0range(1,iele,neemtype) = qeqq0(iele)
        radrange(1,iele,neemtype) = qeqrad(iele)
      elseif (neemtype.eq.3) then
!
!  Streitz and Mintmire
!
        chirange(1,iele,neemtype) = smchi(iele)
        murange(1,iele,neemtype) = smmu(iele)
        q0range(1,iele,neemtype) = smq0(iele)
        zetarange(1,iele,neemtype) = smzeta(iele)
        znucrange(1,iele,neemtype) = smZnuc(iele)
      elseif (neemtype.eq.4) then
!
!  Pacha
!
        chirange(1,iele,neemtype) = chi_pacha(iele)
        murange(1,iele,neemtype) = 0.5_dp*inverse_angstroms_to_ev/rad_pacha(i)
        q0range(1,iele,neemtype) = q0_pacha(iele)
      else
!
!  EEM
!
        chirange(1,iele,neemtype) = chi(iele)
        murange(1,iele,neemtype) = rmu(iele)
        q0range(1,iele,neemtype) = q0(iele)
      endif
      if (abs(chirange(1,iele,neemtype)).gt.1.0d-6.or.abs(murange(1,iele,neemtype)).gt.1.0d-6) then
        lelementOK(iele) = .true.
      endif
    elseif (lfound) then
      if (nqrange(iele,neemtype).gt.1) lmultiqrange = .true.
      lelementOK(iele) = .true.
!########################################################
!  Check parameters to ensure that ranges are sensible  #
!########################################################
!
!  Is there a case where qmin > qmax?
!
      do i = 1,nqrange(iele,neemtype)
        if (nqrangetype(i,iele,neemtype).eq.3) then
          if (qrangemin(i,iele,neemtype).gt.qrangemax(i,iele,neemtype)) then
            call outerror('charge parameters with qmin greater than qmax',0_i4)
            call stopnow('seteem')
          endif
        endif
      enddo
      if (nqrange(iele,neemtype).gt.1) then
!
!  Is there more than one range when global parameters are specified?
!
        lglobal = any(nqrangetype(1:nqrange(iele,neemtype),iele,neemtype).eq.0_i4)
        if (lglobal) then
          call outerror('charge parameters with overlapping q ranges',0_i4)
          call stopnow('seteem')
        endif
!
!  Do any pair of parameters have overlapping ranges?
!
        do i = 1,nqrange(iele,neemtype)
          if (nqrangetype(i,iele,neemtype).eq.1) then
!
!  Minimum range
!
            do j = 1,nqrange(iele,neemtype)
              if (i.ne.j) then
                if (nqrangetype(j,iele,neemtype).eq.1) then
                  call outerror('only one set of charge parameters can have qmin only',0_i4)
                  call stopnow('seteem')
                elseif (nqrangetype(j,iele,neemtype).eq.2) then
                  if (qrangemin(i,iele,neemtype).lt.qrangemax(j,iele,neemtype)) then
                    call outerror('charge parameters with overlapping q ranges',0_i4)
                    call stopnow('seteem')
                  endif
                elseif (nqrangetype(j,iele,neemtype).eq.3) then
                  if (qrangemin(i,iele,neemtype).lt.qrangemax(j,iele,neemtype)) then
                    call outerror('charge parameters with overlapping q ranges',0_i4)
                    call stopnow('seteem')
                  endif
                endif
              endif
            enddo
          elseif (nqrangetype(i,iele,neemtype).eq.2) then
!
!  Maximum range
!
            do j = 1,nqrange(iele,neemtype)
              if (i.ne.j) then
                if (nqrangetype(j,iele,neemtype).eq.2) then
                  call outerror('only one set of charge parameters can have qmax only',0_i4)
                  call stopnow('seteem')
                elseif (nqrangetype(j,iele,neemtype).eq.1) then
                  if (qrangemax(i,iele,neemtype).gt.qrangemin(j,iele,neemtype)) then
                    call outerror('charge parameters with overlapping q ranges',0_i4)
                    call stopnow('seteem')
                  endif
                elseif (nqrangetype(j,iele,neemtype).eq.3) then
                  if (qrangemax(i,iele,neemtype).gt.qrangemin(j,iele,neemtype)) then
                    call outerror('charge parameters with overlapping q ranges',0_i4)
                    call stopnow('seteem')
                  endif
                endif
              endif
            enddo
          elseif (nqrangetype(i,iele,neemtype).eq.3) then
!         
!  Maximum and minimum range
!     
            do j = 1,nqrange(iele,neemtype)
              if (i.ne.j) then
                if (nqrangetype(j,iele,neemtype).eq.1) then
                  if (qrangemin(j,iele,neemtype).lt.qrangemax(i,iele,neemtype)) then
                    call outerror('charge parameters with overlapping q ranges',0_i4)
                    call stopnow('seteem')
                  endif
                elseif (nqrangetype(j,iele,neemtype).eq.2) then
                  if (qrangemax(j,iele,neemtype).gt.qrangemin(i,iele,neemtype)) then
                    call outerror('charge parameters with overlapping q ranges',0_i4)
                    call stopnow('seteem')
                  endif
                elseif (nqrangetype(j,iele,neemtype).eq.3) then
                  if (qrangemin(i,iele,neemtype).gt.qrangemin(j,iele,neemtype)) then
                    if (qrangemax(j,iele,neemtype).gt.qrangemin(i,iele,neemtype)) then
                      call outerror('charge parameters with overlapping q ranges',0_i4)
                      call stopnow('seteem')
                    endif
                  else
                    if (qrangemax(i,iele,neemtype).gt.qrangemin(j,iele,neemtype)) then
                      call outerror('charge parameters with overlapping q ranges',0_i4)
                      call stopnow('seteem')
                    endif
                  endif
                endif
              endif
            enddo
          endif
        enddo
      endif
    endif
  enddo
#ifdef TRACE
  call trace_out('seteem')
#endif
!
  return
  end
