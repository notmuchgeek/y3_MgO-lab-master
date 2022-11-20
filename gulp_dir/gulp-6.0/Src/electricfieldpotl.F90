  subroutine electricfieldpotl(Vfield)
!
!  Subroutine for calculating the potential due to an external
!  electric field. Note that this can should only be done in non-periodic
!  directions and that the potential is defined as the integral of the
!  field/per charge relative to the initial geometry.
!
!   8/11 Created from electricfield.f90
!  12/12 Time-dependent field added
!  12/12 Delay and end added for time-dependent field
!  12/12 Modified to allow for multiple time-dependent fields
!   7/17 Global sum of vfield added
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
  use g_constants,    only : pi
  use current
  use field
  use general,        only : timesofar
  use mdlogic,        only : lmd
  use moldyn,         only : tmdfieldstart, tmdfieldstop
  use parallel,       only : nprocs, procid
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  real(dp),                      intent(out)   :: Vfield(*)
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: status
  real(dp)                                     :: fieldx
  real(dp)                                     :: fieldy
  real(dp)                                     :: fieldz
  real(dp)                                     :: fnorm
  real(dp)                                     :: tdf
  real(dp)                                     :: tsf
  real(dp)                                     :: twopi
  real(dp),    dimension(:), allocatable       :: tmp
!
  Vfield(1:numat) = 0.0_dp
!
!  If number of atoms is zero then return as there is nothing to do
!
  if (numat.eq.0) return
#ifdef TRACE
  call trace_in('electricfieldpotl')
#endif
!
!  Static electric field
!
  if (lfieldcfg(ncf)) then
!
!  For 3-D convert the field direction from cell vector to Cartesian
!
    if (ndim.eq.3) then
      fieldx = fielddirectioncfg(1,ncf)*rv(1,1) + fielddirectioncfg(2,ncf)*rv(1,2) + fielddirectioncfg(3,ncf)*rv(1,3)
      fieldy = fielddirectioncfg(1,ncf)*rv(2,1) + fielddirectioncfg(2,ncf)*rv(2,2) + fielddirectioncfg(3,ncf)*rv(2,3)
      fieldz = fielddirectioncfg(1,ncf)*rv(3,1) + fielddirectioncfg(2,ncf)*rv(3,2) + fielddirectioncfg(3,ncf)*rv(3,3)
      fnorm = fieldx**2 + fieldy**2 + fieldz**2
      fnorm = 1.0_dp/sqrt(fnorm)
      fieldx = fieldcfg(ncf)*fieldx*fnorm
      fieldy = fieldcfg(ncf)*fieldy*fnorm
      fieldz = fieldcfg(ncf)*fieldz*fnorm
    else
!
!  Find norm of field direction and scale components
!
      fnorm = fielddirectioncfg(1,ncf)**2 + fielddirectioncfg(2,ncf)**2 + fielddirectioncfg(3,ncf)**2
      fnorm = 1.0_dp/sqrt(fnorm)
      fieldx = fieldcfg(ncf)*fielddirectioncfg(1,ncf)*fnorm
      fieldy = fieldcfg(ncf)*fielddirectioncfg(2,ncf)*fnorm
      fieldz = fieldcfg(ncf)*fielddirectioncfg(3,ncf)*fnorm
    endif
  else
    fieldx = 0.0_dp
    fieldy = 0.0_dp
    fieldz = 0.0_dp
  endif
!
!  Time-dependent electric field
!
  if (lmd.and.(ntdfieldcfg(ncf).gt.0)) then
    tsf = timesofar - tmdfieldstart(ncf)
    if (tsf.gt.0.0_dp.and.(tsf.lt.tmdfieldstop(ncf).or.tmdfieldstop(ncf).eq.0.0_dp)) then
      twopi = 2.0_dp*pi
      do i = 1,ntdfieldcfg(ncf)
        tdf = td_fieldcfg(1,i,ncf)*cos(twopi*(tsf*td_fieldcfg(2,i,ncf) + td_fieldcfg(3,i,ncf)))
        fnorm = td_fielddirectioncfg(1,i,ncf)**2 + td_fielddirectioncfg(2,i,ncf)**2 + td_fielddirectioncfg(3,i,ncf)**2
        fnorm = 1.0_dp/sqrt(fnorm)
        fieldx = fieldx + tdf*td_fielddirectioncfg(1,i,ncf)*fnorm
        fieldy = fieldy + tdf*td_fielddirectioncfg(2,i,ncf)*fnorm
        fieldz = fieldz + tdf*td_fielddirectioncfg(3,i,ncf)*fnorm
      enddo
    endif
  endif
!
!  Loop over asymmetric unit performing integral
!
  do i = 1+procid,numat,nprocs
    Vfield(i) = Vfield(i) + (fieldx*xalat(i) + fieldy*yalat(i) + fieldz*zalat(i))
  enddo
!
  if (nprocs.gt.1) then
    allocate(tmp(numat),stat=status)
    if (status/=0) call outofmemory('electricfieldpotl','tmp')
!
    call sumall(Vfield,tmp,numat,"electricfieldpotl","Vfield")
    Vfield(1:numat) = tmp(1:numat)
!
    deallocate(tmp,stat=status)
    if (status/=0) call deallocate_error('electricfieldpotl','tmp')
  endif
#ifdef TRACE
  call trace_out('electricfieldpotl')
#endif
!
  return
  end
