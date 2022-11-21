  subroutine cosmoamatdadd1p(ipts,jpts,f0,f1,f2,qsij,qsipj,sumAinvipts,sumAinvjpts,ldqneeded, &
                             dcosmoA,dsegweight,nearsas,nearsasrptr,dcosmoA2)
!
!  Subroutine adds contributions from a term in the cosmo A matrix to the first 
!  derivatives, as appropriate.
!  Designed to work in parallel
!
!  On input :
!
!  ipts     = first segment pointer, associated with atom i
!  jpts     = second segment pointer, associated with atom j
!  f0      = A matrix term
!  f1(3)   = Cartesian first derivatives of A matrix term with respect to rij
!  f2(6)   = Cartesian second derivatives of A matrix term with respect to rij
!
!   4/17 Created from cosmoamatdadd
!   2/18 Trace added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
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
!  Julian Gale, CIC, Curtin University, September 2019
!
  use configurations, only : nregionno
  use cosmic
  use current,        only : nsft, nrelf2a
  use derivatives
  use parallel
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4)                                   :: ipts
  integer(i4)                                   :: iptsloc
  integer(i4)                                   :: jpts
  integer(i4)                                   :: nearsas
  integer(i4)                                   :: nearsasrptr(*)
  logical                                       :: ldqneeded
  real(dp)                                      :: dcosmoA(3,nptsonnode,*)
  real(dp)                                      :: dcosmoA2(3,nptsonnode,*)
  real(dp)                                      :: dsegweight(3,maxnearseg,*)
  real(dp)                                      :: f0
  real(dp)                                      :: f1(3)
  real(dp)                                      :: f2(6)
  real(dp)                                      :: qsij
  real(dp)                                      :: qsipj
  real(dp)                                      :: sumAinvipts
  real(dp)                                      :: sumAinvjpts
!
!  Local variables
!
  integer(i4)                                   :: i
  integer(i4)                                   :: ia
  integer(i4)                                   :: in
  integer(i4)                                   :: j
  integer(i4)                                   :: ja
  integer(i4)                                   :: jn
  integer(i4)                                   :: m
  integer(i4)                                   :: mm
  integer(i4)                                   :: mn
  integer(i4)                                   :: nregioni
  integer(i4)                                   :: nregionj
  integer(i4)                                   :: nregionm
  real(dp)                                      :: qsi
  real(dp)                                      :: qsj
  real(dp)                                      :: qsijw
  real(dp)                                      :: qsipjw
  real(dp)                                      :: smtrm
  real(dp)                                      :: smtrmi
  real(dp)                                      :: smtrmj
  real(dp)                                      :: swi
  real(dp)                                      :: swj
#ifdef TRACE
  call trace_in('cosmoamatdadd1p')
#endif
!
!  Set local constants
!
  iptsloc = npts2local(ipts)
  i = cosmoatomptr(ipts)
  ia = nrelf2a(i)
  in = nearsasrptr(i)
  qsi = qonsas(ipts)
  swi = segweight(ipts)
  nregioni = nregionno(nsft+ia)
!
  j = cosmoatomptr(jpts)
  ja = nrelf2a(j)
  jn = nearsasrptr(j)
  qsj = qonsas(jpts)
  swj = segweight(jpts)
  nregionj = nregionno(nsft+ja)
!
  qsijw = qsij*swi*swj
  qsipjw = qsipj*swi*swj
!**************************************
!  First derivatives of the A matrix  *
!**************************************
  if (i.ne.j) then
    xdrv(i) = xdrv(i) - qsijw*f1(1)
    ydrv(i) = ydrv(i) - qsijw*f1(2)
    zdrv(i) = zdrv(i) - qsijw*f1(3)
    xdrv(j) = xdrv(j) + qsijw*f1(1)
    ydrv(j) = ydrv(j) + qsijw*f1(2)
    zdrv(j) = zdrv(j) + qsijw*f1(3)
    if (nregioni.ne.nregionj) then
      xregdrv(nregioni) = xregdrv(nregioni) - qsijw*f1(1)
      yregdrv(nregioni) = yregdrv(nregioni) - qsijw*f1(2)
      zregdrv(nregioni) = zregdrv(nregioni) - qsijw*f1(3)
      xregdrv(nregionj) = xregdrv(nregionj) + qsijw*f1(1)
      yregdrv(nregionj) = yregdrv(nregionj) + qsijw*f1(2)
      zregdrv(nregionj) = zregdrv(nregionj) + qsijw*f1(3)
    endif
  endif
!
!  Derivatives with respect to smoothing of segments
!
  if (lsegsmooth) then
    if (ipts.eq.jpts) then
      smtrm = 0.5_dp*qsij*f0
    else
      smtrm = qsij*f0
    endif
    smtrmi = smtrm*swi
    smtrmj = smtrm*swj
    do mm = 1,nnearseg(ipts)
      m = nnearsegptr(mm,ipts)
      xdrv(i) = xdrv(i) - smtrmj*dsegweight(1,mm,ipts)
      ydrv(i) = ydrv(i) - smtrmj*dsegweight(2,mm,ipts)
      zdrv(i) = zdrv(i) - smtrmj*dsegweight(3,mm,ipts)
      xdrv(m) = xdrv(m) + smtrmj*dsegweight(1,mm,ipts)
      ydrv(m) = ydrv(m) + smtrmj*dsegweight(2,mm,ipts)
      zdrv(m) = zdrv(m) + smtrmj*dsegweight(3,mm,ipts)
      nregionm = nregionno(nsft+nrelf2a(m))
      if (nregioni.ne.nregionm) then
        xregdrv(nregioni) = xregdrv(nregioni) - smtrmj*dsegweight(1,mm,ipts)
        yregdrv(nregioni) = yregdrv(nregioni) - smtrmj*dsegweight(2,mm,ipts)
        zregdrv(nregioni) = zregdrv(nregioni) - smtrmj*dsegweight(3,mm,ipts)
        xregdrv(nregionm) = xregdrv(nregionm) + smtrmj*dsegweight(1,mm,ipts)
        yregdrv(nregionm) = yregdrv(nregionm) + smtrmj*dsegweight(2,mm,ipts)
        zregdrv(nregionm) = zregdrv(nregionm) + smtrmj*dsegweight(3,mm,ipts)
      endif
    enddo
    do mm = 1,nnearseg(jpts)
      m = nnearsegptr(mm,jpts)
      xdrv(j) = xdrv(j) - smtrmi*dsegweight(1,mm,jpts)
      ydrv(j) = ydrv(j) - smtrmi*dsegweight(2,mm,jpts)
      zdrv(j) = zdrv(j) - smtrmi*dsegweight(3,mm,jpts)
      xdrv(m) = xdrv(m) + smtrmi*dsegweight(1,mm,jpts)
      ydrv(m) = ydrv(m) + smtrmi*dsegweight(2,mm,jpts)
      zdrv(m) = zdrv(m) + smtrmi*dsegweight(3,mm,jpts)
      nregionm = nregionno(nsft+nrelf2a(m))
      if (nregionj.ne.nregionm) then
        xregdrv(nregionj) = xregdrv(nregionj) - smtrmi*dsegweight(1,mm,jpts)
        yregdrv(nregionj) = yregdrv(nregionj) - smtrmi*dsegweight(2,mm,jpts)
        zregdrv(nregionj) = zregdrv(nregionj) - smtrmi*dsegweight(3,mm,jpts)
        xregdrv(nregionm) = xregdrv(nregionm) + smtrmi*dsegweight(1,mm,jpts)
        yregdrv(nregionm) = yregdrv(nregionm) + smtrmi*dsegweight(2,mm,jpts)
        zregdrv(nregionm) = zregdrv(nregionm) + smtrmi*dsegweight(3,mm,jpts)
      endif
    enddo
  endif
!
  if (ldqneeded) then
!
!  Save first derivatives of A matrix term
!
    if (i.ne.j) then
      dcosmoA(1,iptsloc,jn) = dcosmoA(1,iptsloc,jn) - f1(1)*qsj*swi*swj
      dcosmoA(2,iptsloc,jn) = dcosmoA(2,iptsloc,jn) - f1(2)*qsj*swi*swj
      dcosmoA(3,iptsloc,jn) = dcosmoA(3,iptsloc,jn) - f1(3)*qsj*swi*swj
    endif
!
!  Smoothing contribution to dcosmoA
!
    if (lsegsmooth) then
      if (ipts.eq.jpts) then
        smtrm = 0.5_dp*f0
      else
        smtrm = f0
      endif
      smtrmi = smtrm*qsi*swj
      smtrmj = smtrm*qsj*swj
      do mm = 1,nnearseg(ipts)
        m = nnearsegptr(mm,ipts)
        mn = nearsasrptr(m)
        dcosmoA(1,iptsloc,mn) = dcosmoA(1,iptsloc,mn) - smtrmj*dsegweight(1,mm,ipts)
        dcosmoA(2,iptsloc,mn) = dcosmoA(2,iptsloc,mn) - smtrmj*dsegweight(2,mm,ipts)
        dcosmoA(3,iptsloc,mn) = dcosmoA(3,iptsloc,mn) - smtrmj*dsegweight(3,mm,ipts)
      enddo
      smtrmi = smtrm*qsi*swi
      smtrmj = smtrm*qsj*swi
      do mm = 1,nnearseg(jpts)
        m = nnearsegptr(mm,jpts)
        mn = nearsasrptr(m)
        dcosmoA2(1,iptsloc,mn) = dcosmoA2(1,iptsloc,mn) - smtrmj*dsegweight(1,mm,jpts)
        dcosmoA2(2,iptsloc,mn) = dcosmoA2(2,iptsloc,mn) - smtrmj*dsegweight(2,mm,jpts)
        dcosmoA2(3,iptsloc,mn) = dcosmoA2(3,iptsloc,mn) - smtrmj*dsegweight(3,mm,jpts)
        dcosmoA2(1,iptsloc,jn) = dcosmoA2(1,iptsloc,jn) + smtrmj*dsegweight(1,mm,jpts)
        dcosmoA2(2,iptsloc,jn) = dcosmoA2(2,iptsloc,jn) + smtrmj*dsegweight(2,mm,jpts)
        dcosmoA2(3,iptsloc,jn) = dcosmoA2(3,iptsloc,jn) + smtrmj*dsegweight(3,mm,jpts)
      enddo
    endif
  endif
#ifdef TRACE
  call trace_out('cosmoamatdadd1p')
#endif
!
  return
  end
