  subroutine mcstorecfg
!
!  Stores final configuration back in main structure arrays.
!  Uses the assumption that there is only one configuration
!  present as per much of the rest of the Monte Carlo code.
!
!   2/18 Trace added
!   3/19 iopt changed to ioptindex and iopttype
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
  use configurations
  use current
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
#ifdef TRACE
  call trace_in('mcstorecfg')
#endif
!
!  Check dimensions
!
  if (nasym.gt.maxatot) then
    maxatot = nasym + 10
    call changemaxatot
  endif
  if (nvar.gt.maxvar) then
    maxvar = nvar
    call changemaxvar
  endif
!
!  Copy back relevant current data into configuration arrays
!
  nascfg(ncf) = nasym
  natcfg(1:nasym) = iatn(1:nasym)
  ntypcfg(1:nasym) = natype(1:nasym)
  cncfg(1:nasym) = cna(1:nasym)
  occucfg(1:nasym) = occua(1:nasym)
  oxcfg(1:nasym) = oxa(1:nasym)
  qlcfg(1:nasym) = qa(1:nasym)
  radcfg(1:nasym) = rada(1:nasym)
  xcfg(1:nasym) = xafrac(1:nasym)
  ycfg(1:nasym) = yafrac(1:nasym)
  zcfg(1:nasym) = zafrac(1:nasym)
  ioptindexcfg(1:nvar) = ioptindex(1:nvar)
  iopttypecfg(1:nvar) = iopttype(1:nvar)
#ifdef TRACE
  call trace_out('mcstorecfg')
#endif
!
  return
  end
