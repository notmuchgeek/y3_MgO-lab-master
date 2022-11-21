  subroutine special(loptx,lopty,loptz,lcheck)
!
!  Check which atomic variables are special positions
!  Displace each coordinate in turn and see if the
!  site multiplicity is changed. This may not be elegant
!  but it works for all space groups!
!
!  ltmp   = array of logical flags according to whether a
!           parameter can be varied or not. Passed from setcfg.
!  lcheck = logical variable. If true then setting of ltmp is
!           tested for correctness rather than set.
!
!   6/95 Automatic creation of constraints for partial occupancy
!        sites added.
!   3/03 Error in incrementing of maxcontot fixed
!   4/04 Checks on nrotop modified to avoid out of bounds reference
!  12/05 Modified to handle case where general space group operators
!        have been input
!  11/06 lcheck mode added & ftow replaced with itow
!  11/06 Setting of indx/y/z corrected and warning messages standardised
!   3/07 Calls to mxmb renamed to GULP_mxmb
!   2/09 Modified so that only constraints allowed by the flags are added
!   3/09 Position of indx, indy and indz calculation returned to old one
!        since move led to errors in flags.
!   2/18 Trace added
!   3/19 Constraint arrays changed to have index and type
!   3/19 ltmp changed to loptx, lopty, loptz
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!   2/20 Modifications added to skip atoms in rigid molecules
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
!  Julian Gale, CIC, Curtin University, February 2020
!
  use configurations
  use control
  use current
  use general,       only : nwarn
  use iochannels
  use molecule,      only : natinmol
  use optimisation
  use parallel
  use symmetry
  use times
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  logical, intent(inout) :: loptx(*)
  logical, intent(inout) :: lopty(*)
  logical, intent(inout) :: loptz(*)
  logical, intent(in)    :: lcheck
!
!  Local variables
!
  character(len=3)       :: conno
  character(len=6)       :: atomno
  integer(i4)            :: ii
  integer(i4)            :: indi
  integer(i4)            :: j
  integer(i4)            :: k
  integer(i4)            :: mv
  integer(i4)            :: na
  integer(i4)            :: ncv
  integer(i4)            :: neqvna
  integer(i4)            :: nfinish
  integer(i4)            :: nfirst
  integer(i4)            :: nfk
  integer(i4)            :: ngen
  integer(i4)            :: nin
  integer(i4)            :: nposs
  integer(i4)            :: nvk
  integer(i4)            :: nvk1
  integer(i4)            :: nvk2
  integer(i4)            :: nxy1
  integer(i4)            :: nxz1
  integer(i4)            :: nyz1
  logical                :: lfound
  logical                :: lfoundxy
  logical                :: lfoundxz
  logical                :: lfoundyz
  logical                :: lnocon
  logical                :: ltmploc(3)
  logical                :: lxposs
  logical                :: lyposs
  logical                :: lzposs
  real(dp)               :: coeff
  real(dp)               :: g_cpu_time
  real(dp)               :: diff
  real(dp)               :: rxp
  real(dp)               :: ryp
  real(dp)               :: rzp
  real(dp)               :: thresh
  real(dp)               :: time1
  real(dp)               :: time2
  real(dp)               :: v(3)
  real(dp)               :: x(3)
  real(dp)               :: xorig
  real(dp)               :: xposs(10)
  real(dp)               :: xstor(48)
  real(dp)               :: xx(3)
  real(dp)               :: yorig
  real(dp)               :: yposs(10)
  real(dp)               :: ystor(48)
  real(dp)               :: zorig
  real(dp)               :: zposs(10)
  real(dp)               :: zstor(48)
#ifdef TRACE
  call trace_in('special')
#endif
!
  time1 = g_cpu_time()
  thresh = 1.0d-4
!
!  Setup local variables
!
  lnocon = (index(keyword,'noco').ne.0)
!********************************************
!  Setup combination shifts of coordinates  *
!********************************************
!
!  Space group input case
!
  if (nspcg(ncf).gt.1) then
    if (nccs.eq.4) then
      nposs = 2
      xposs(1) = 0.0001_dp
      yposs(1) = 0.0001_dp
      zposs(1) = 0.0_dp
      xposs(2) = 0.0001_dp
      yposs(2) = -0.0001_dp
      zposs(2) = 0.0_dp
    elseif (nccs.eq.5) then
      nposs = 4
      xposs(1) = 0.0001_dp
      yposs(1) = 0.0001_dp
      zposs(1) = 0.0_dp
      xposs(2) = 0.0001_dp
      yposs(2) = -0.0001_dp
      zposs(2) = 0.0_dp
      xposs(3) = 0.0001_dp
      yposs(3) = 0.0002_dp
      zposs(3) = 0.0_dp
      xposs(4) = 0.0002_dp
      yposs(4) = 0.0001_dp
      zposs(4) = 0.0_dp
    elseif (nccs.eq.6) then
      if (nspcg(ncf).lt.207) then
        nposs = 4
        xposs(1) = 0.0001_dp
        yposs(1) = 0.0001_dp
        zposs(1) = 0.0001_dp
        xposs(2) = 0.0001_dp
        yposs(2) = 0.0001_dp
        zposs(2) = -0.0001_dp
        xposs(3) = 0.0001_dp
        yposs(3) = -0.0001_dp
        zposs(3) = 0.0001_dp
        xposs(4) = -0.0001_dp
        yposs(4) = 0.0001_dp
        zposs(4) = 0.0001_dp
      else
        nposs = 10
        xposs(1) = 0.0001_dp
        yposs(1) = 0.0001_dp
        zposs(1) = 0.0001_dp
        xposs(2) = 0.0001_dp
        yposs(2) = 0.0001_dp
        zposs(2) = -0.0001_dp
        xposs(3) = 0.0001_dp
        yposs(3) = -0.0001_dp
        zposs(3) = 0.0001_dp
        xposs(4) = -0.0001_dp
        yposs(4) = 0.0001_dp
        zposs(4) = 0.0001_dp
        xposs(5) = 0.0001_dp
        yposs(5) = 0.0001_dp
        zposs(5) = 0.0_dp
        xposs(6) = 0.0001_dp
        yposs(6) = -0.0001_dp
        zposs(6) = 0.0_dp
        xposs(7) = 0.0001_dp
        yposs(7) = 0.0_dp
        zposs(7) = 0.0001_dp
        xposs(8) = 0.0001_dp
        yposs(8) = 0.0_dp
        zposs(8) = -0.0001_dp
        xposs(9) = 0.0_dp
        yposs(9) = 0.0001_dp
        zposs(9) = 0.0001_dp
        xposs(10) = 0.0_dp
        yposs(10) = 0.0001_dp
        zposs(10) = -0.0001_dp
      endif
    else
      nposs = 0
    endif
  elseif (ngocfg(ncf).gt.1) then
!
!  General operator input case
!
    nposs = 10
    xposs(1) = 0.0001_dp
    yposs(1) = 0.0001_dp
    zposs(1) = 0.0001_dp
    xposs(2) = 0.0001_dp
    yposs(2) = 0.0001_dp
    zposs(2) = -0.0001_dp
    xposs(3) = 0.0001_dp
    yposs(3) = -0.0001_dp
    zposs(3) = 0.0001_dp
    xposs(4) = -0.0001_dp
    yposs(4) = 0.0001_dp
    zposs(4) = 0.0001_dp
    xposs(5) = 0.0001_dp
    yposs(5) = 0.0001_dp
    zposs(5) = 0.0_dp
    xposs(6) = 0.0001_dp
    yposs(6) = -0.0001_dp
    zposs(6) = 0.0_dp
    xposs(7) = 0.0001_dp
    yposs(7) = 0.0_dp
    zposs(7) = 0.0001_dp
    xposs(8) = 0.0001_dp
    yposs(8) = 0.0_dp
    zposs(8) = -0.0001_dp
    xposs(9) = 0.0_dp
    yposs(9) = 0.0001_dp
    zposs(9) = 0.0001_dp
    xposs(10) = 0.0_dp
    yposs(10) = 0.0001_dp
    zposs(10) = -0.0001_dp
  else
!
!  No symmetry case
!
    nposs = 0
  endif
!**************************************
!  Test for each atom and coordinate  *
!**************************************
  indi = 0
  naloop: do na = 1,nasym
!
!  Check for rigid molecules
!
    if (lrigid) then
      if (natinmol(nrela2f(na)).gt.0) then
        loptx(indi+1) = .false.
        lopty(indi+1) = .false.
        loptz(indi+1) = .false.
        cycle naloop
      endif
    endif
!
    xorig = xafrac(na)
    yorig = yafrac(na)
    zorig = zafrac(na)
    nfirst = nrela2f(na)
    neqvna = neqv(na)
!***************************************
!  First pass - individual x, y and z  *
!***************************************
    iiloop: do ii = 1,3
      xx(1) = xorig
      xx(2) = yorig
      xx(3) = zorig
      xx(ii) = xx(ii) + 0.0001_dp
!
!  First symmetry operator
!
      x(1) = 0.0_dp
      x(2) = 0.0_dp
      x(3) = 0.0_dp
      v(1) = xx(1)
      v(2) = xx(2)
      v(3) = xx(3)
!
!  Transform coordinates to primitive cell
!
      call GULP_mxmb(w1(ncbl,1,1),7_i4,21_i4,v,1_i4,1_i4,x,1_i4,1_i4,3_i4,3_i4,1_i4)
!
!  Place fractional coordinates in range 0-1
!
      x(1) = x(1) + 3.0_dp
      nin = x(1)
      xstor(1) = x(1) - nin
      x(2) = x(2) + 3.0_dp
      nin = x(2)
      ystor(1) = x(2) - nin
      x(3) = x(3) + 3.0_dp
      nin = x(3)
      zstor(1) = x(3) - nin
!
!  Loop over symmetry operators
!
      ngen = 1
      ngoloop: do mv = 2,ngo
!
!  Roto-translation
!
        x(1) = 0.0_dp
        x(2) = 0.0_dp
        x(3) = 0.0_dp
        v(1) = vit(1,mv)
        v(2) = vit(2,mv)
        v(3) = vit(3,mv)
        call GULP_mxmb(rop(1,1,mv),1_i4,3_i4,xx,1_i4,1_i4,v,1_i4,1_i4,3_i4,3_i4,1_i4)
!
!  Transform coordinates to primitive cell
!
        call GULP_mxmb(w1(ncbl,1,1),7_i4,21_i4,v,1_i4,1_i4,x,1_i4,1_i4,3_i4,3_i4,1_i4)
!
!  Place fractional coordinates in range 0-1
!
        x(1) = x(1) + 3.0_dp
        nin = x(1)
        x(1) = x(1) - nin
        x(2) = x(2) + 3.0_dp
        nin = x(2)
        x(2) = x(2) - nin
        x(3) = x(3) + 3.0_dp
        nin = x(3)
        x(3) = x(3) - nin
!
!  Compare atom with previously generated equivalent
!
        nfinish = nfirst + ngen - 1
        do j = 1,ngen
          v(1) = xstor(j) - x(1)
          v(2) = ystor(j) - x(2)
          v(3) = zstor(j) - x(3)
          v(1) = v(1) - nint(v(1))
          v(2) = v(2) - nint(v(2))
          v(3) = v(3) - nint(v(3))
          if ((abs(v(1))+abs(v(2))+abs(v(3))).lt.thresh) cycle ngoloop
        enddo
!
!  Atom is not equivalent to any previous atom.
!  If new unique atom has been generated by the same space group
!  operator as for the original coordinates then there has been
!  no change in multiplicity.
!
        ngen = ngen + 1
        if (ngen.gt.neqvna) then
          ltmploc(ii) = .false.
          cycle iiloop
        elseif (nrotop(nfinish+1).ne.mv) then
          ltmploc(ii) = .false.
          cycle iiloop
        endif
        xstor(ngen) = x(1)
        ystor(ngen) = x(2)
        zstor(ngen) = x(3)
      enddo ngoloop
      ltmploc(ii) = .true.
    enddo iiloop
    indi = indi + 1
!*********************************************
!  Second pass - combinations of x, y and z  *
!*********************************************
    if (nposs.gt.0) then
      npossloop: do ii = 1,nposs
        xx(1) = xorig + xposs(ii)
        xx(2) = yorig + yposs(ii)
        xx(3) = zorig + zposs(ii)
!
!  Only test combinations where shifted coordinates haven't
!  already been flagged for optimisation
!
        lxposs = (abs(xposs(ii)).gt.1.0d-6)
        lyposs = (abs(yposs(ii)).gt.1.0d-6)
        lzposs = (abs(zposs(ii)).gt.1.0d-6)
        if (ltmploc(1).and.lxposs) cycle npossloop
        if (ltmploc(2).and.lyposs) cycle npossloop
        if (ltmploc(3).and.lzposs) cycle npossloop
!****************************
!  First symmetry operator  *
!****************************
        x(1) = 0.0_dp
        x(2) = 0.0_dp
        x(3) = 0.0_dp
        v(1) = xx(1)
        v(2) = xx(2)
        v(3) = xx(3)
!
!  Transform coordinates to primitive cell
!
        call GULP_mxmb(w1(ncbl,1,1),7_i4,21_i4,v,1_i4,1_i4,x,1_i4,1_i4,3_i4,3_i4,1_i4)
!
!  Place fractional coordinates in range 0-1
!
        x(1) = x(1) + 3.0_dp
        nin = x(1)
        xstor(1) = x(1) - nin
        x(2) = x(2) + 3.0_dp
        nin = x(2)
        ystor(1) = x(2) - nin
        x(3) = x(3) + 3.0_dp
        nin = x(3)
        zstor(1) = x(3) - nin
!*********************************
!  Loop over symmetry operators  *
!*********************************
        ngen = 1
        mvloop: do mv = 2,ngo
!
!  Roto-translation
!
          x(1) = 0.0_dp
          x(2) = 0.0_dp
          x(3) = 0.0_dp
          v(1) = vit(1,mv)
          v(2) = vit(2,mv)
          v(3) = vit(3,mv)
          call GULP_mxmb(rop(1,1,mv),1_i4,3_i4,xx,1_i4,1_i4,v,1_i4,1_i4,3_i4,3_i4,1_i4)
!
!  Transform coordinates to primitive cell
!
          call GULP_mxmb(w1(ncbl,1,1),7_i4,21_i4,v,1_i4,1_i4,x,1_i4,1_i4,3_i4,3_i4,1_i4)
!
!  Place fractional coordinates in range 0-1
!
          x(1) = x(1) + 3.0_dp
          nin = x(1)
          x(1) = x(1) - nin
          x(2) = x(2) + 3.0_dp
          nin = x(2)
          x(2) = x(2) - nin
          x(3) = x(3) + 3.0_dp
          nin = x(3)
          x(3) = x(3) - nin
!
!  Compare atom with previously generated equivalent
!
          nfinish = nfirst + ngen - 1
          do j = 1,ngen
            v(1) = xstor(j) - x(1)
            v(2) = ystor(j) - x(2)
            v(3) = zstor(j) - x(3)
            v(1) = v(1) - nint(v(1))
            v(2) = v(2) - nint(v(2))
            v(3) = v(3) - nint(v(3))
            if ((abs(v(1))+abs(v(2))+abs(v(3))).lt.thresh) cycle mvloop
          enddo
!
!  Atom is not equivalent to any previous atom.
!  If new unique atom has been generated by the same space group
!  operator as for the original coordinates then there has been
!  no change in multiplicity.
!
          ngen = ngen + 1
          if (ngen.gt.neqvna) then
            cycle npossloop
          elseif (nrotop(nfinish+1).ne.mv) then
            cycle npossloop
          endif
          xstor(ngen) = x(1)
          ystor(ngen) = x(2)
          zstor(ngen) = x(3)
        enddo mvloop
!**************************************************
!  Valid combination - set flags and constraints  *
!**************************************************
        if (lnocon) then
          if (lxposs) ltmploc(1) = .true.
          if (lyposs) ltmploc(2) = .true.
          if (lzposs) ltmploc(3) = .true.
        else
          if (lxposs.and.lyposs.and.lzposs) then
!**********************
!  x y z combination  *
!**********************
            rxp = 1.0_dp/xposs(ii)
            ryp = 1.0_dp/yposs(ii)
            rzp = 1.0_dp/zposs(ii)
!
!  Check to see if any constraints already exist for these variables
!
            lfound = .false.
            lfoundxy = .false.
            lfoundyz = .false.
            lfoundxz = .false.
            do k = 1,ncon
              nvk = ncvarindcfg(ncfst+k)
              nfk = ncfixindcfg(ncfst+k)
              if (ncvartypcfg(ncfst+k).eq.iopt_xf.and.ncfixtypcfg(ncfst+k).eq.iopt_yf.and. &
                  nvk.eq.indi.and.nfk.eq.indi) then
                lfoundxy = .true.
                nxy1 = k
              elseif (ncvartypcfg(ncfst+k).eq.iopt_yf.and.ncfixtypcfg(ncfst+k).eq.iopt_xf.and. &
                  nvk.eq.indi.and.nfk.eq.indi) then
                lfoundxy = .true.
                nxy1 = k
              elseif (ncvartypcfg(ncfst+k).eq.iopt_yf.and.ncfixtypcfg(ncfst+k).eq.iopt_zf.and. &
                  nvk.eq.indi.and.nfk.eq.indi) then
                lfoundyz = .true.
                nyz1 = k
              elseif (ncvartypcfg(ncfst+k).eq.iopt_zf.and.ncfixtypcfg(ncfst+k).eq.iopt_yf.and. &
                  nvk.eq.indi.and.nfk.eq.indi) then
                lfoundyz = .true.
                nyz1 = k
              elseif (ncvartypcfg(ncfst+k).eq.iopt_xf.and.ncfixtypcfg(ncfst+k).eq.iopt_zf.and. &
                  nvk.eq.indi.and.nfk.eq.indi) then
                lfoundxz = .true.
                nxz1 = k
              elseif (ncvartypcfg(ncfst+k).eq.iopt_zf.and.ncfixtypcfg(ncfst+k).eq.iopt_xf.and. &
                  nvk.eq.indi.and.nfk.eq.indi) then
                lfoundxz = .true.
                nxz1 = k
              endif
            enddo
            if (lfoundxy.and.lfoundyz.and.lfoundxz) then
!
!  All pairs constrained
!
              lfound = .true.
              if (ioproc) then
                write(ioout,'(''  **** Advice - redundant constraint supplied for atom '',i4,'' ****'')')na
              endif
!
!  If constraints have been found for both pairs
!
            elseif (lfoundxy.and.lfoundxz) then
              nvk1 = ncvarind(nxy1)
              nvk2 = ncvarind(nxz1)
              if (nvk1.ne.nvk2) then
                call itow(atomno,na,6_i4)
                call outerror('Badly defined constraints for '//'atom '//atomno,0_i4)
                call stopnow('special')
              endif
              lfound = .true.
              if (ncvartyp(nxy1).eq.iopt_xf) then
                ltmploc(1) = .true.
              elseif (ncvartyp(nxy1).eq.iopt_yf) then
                ltmploc(2) = .true.
              elseif (ncvartyp(nxy1).eq.iopt_zf) then
                ltmploc(3) = .true.
              endif
            elseif (lfoundxy.and.lfoundyz) then
              nvk1 = ncvarind(nxy1)
              nvk2 = ncvarind(nyz1)
              if (nvk1.ne.nvk2) then
                call itow(atomno,na,6_i4)
                call outerror('Badly defined constraints for '//'atom '//atomno,0_i4)
                call stopnow('special')
              endif
              lfound = .true.
              if (ncvartyp(nxy1).eq.iopt_xf) then
                ltmploc(1) = .true.
              elseif (ncvartyp(nxy1).eq.iopt_yf) then
                ltmploc(2) = .true.
              elseif (ncvartyp(nxy1).eq.iopt_zf) then
                ltmploc(3) = .true.
              endif
            elseif (lfoundxz.and.lfoundyz) then
              nvk1 = ncvarind(nxz1)
              nvk2 = ncvarind(nyz1)
              if (nvk1.ne.nvk2) then
                call itow(atomno,na,6_i4)
                call outerror('Badly defined constraints for '//'atom '//atomno,0_i4)
                call stopnow('special')
              endif
              lfound = .true.
              if (ncvartyp(nxz1).eq.iopt_xf) then
                ltmploc(1) = .true.
              elseif (ncvartyp(nxz1).eq.iopt_yf) then
                ltmploc(2) = .true.
              elseif (ncvartyp(nxz1).eq.iopt_zf) then
                ltmploc(3) = .true.
              endif
!
!  If constraints have been found for one pair
!
            elseif (lfoundxy.and.loptx(indi).and.lopty(indi)) then
              ncv = ncvarind(nxy1)
              if (ncv.eq.indi.and.ncvartyp(nxy1).eq.iopt_xf) then
                if (ncontot.ge.maxcontot) then
                  maxcontot = ncontot + 10
                  call changemaxcontot
                endif
                if (ncf.lt.ncfg) then
                  do k = ncontot,n1con(ncf+1),-1
                    ncvarind(k+1) = ncvarind(k)
                    ncvartyp(k+1) = ncvartyp(k)
                    ncfixind(k+1) = ncfixind(k)
                    ncfixtyp(k+1) = ncfixtyp(k)
                    conco(k+1) = conco(k)
                    nconcfg(k+1) = nconcfg(k)
                    conadd(k+1) = conadd(k)
                  enddo
                  do k = ncf+1,ncfg
                    n1con(k) = n1con(k) + 1
                  enddo
                endif
                ncontot = ncontot + 1
                ltmploc(1) = .true.
                ncon = ncon + 1
                ncvarindcfg(ncfst+ncon) = indi
                ncvartypcfg(ncfst+ncon) = iopt_xf
                ncfixindcfg(ncfst+ncon) = indi
                ncfixtypcfg(ncfst+ncon) = iopt_zf
                concocfg(ncfst+ncon) = zposs(ii)*rxp
                conaddcfg(ncfst+ncon) = zorig - xorig*concocfg(ncfst+ncon)
                nconcfg(ncfst+ncon) = ncf
                lfound = .true.
              else
                if (ncontot.ge.maxcontot) then
                  maxcontot = ncontot + 10
                  call changemaxcontot
                endif
                if (ncf.lt.ncfg) then
                  do k = ncontot,n1con(ncf+1),-1
                    ncvarindcfg(k+1) = ncvarindcfg(k)
                    ncvartypcfg(k+1) = ncvartypcfg(k)
                    ncfixindcfg(k+1) = ncfixindcfg(k)
                    ncfixtypcfg(k+1) = ncfixtypcfg(k)
                    concocfg(k+1) = concocfg(k)
                    conaddcfg(k+1) = conaddcfg(k)
                    nconcfg(k+1) = nconcfg(k)
                  enddo
                  do k = ncf+1,ncfg
                    n1con(k) = n1con(k) + 1
                  enddo
                endif
                ncontot = ncontot + 1
                ltmploc(2) = .true.
                ncon = ncon + 1
                ncvarindcfg(ncfst+ncon) = indi
                ncvartypcfg(ncfst+ncon) = iopt_yf
                ncfixindcfg(ncfst+ncon) = indi
                ncfixtypcfg(ncfst+ncon) = iopt_zf
                concocfg(ncfst+ncon) = zposs(ii)*ryp
                conaddcfg(ncfst+ncon) = zorig - yorig*concocfg(ncfst+ncon)
                nconcfg(ncfst+ncon) = ncf
                lfound = .true.
              endif
            elseif (lfoundyz.and.lopty(indi).and.loptz(indi)) then
              ncv = ncvarind(nyz1)
              if (ncv.eq.indi.and.ncvartyp(nyz1).eq.iopt_yf) then
                if (ncontot.ge.maxcontot) then
                  maxcontot = ncontot + 10
                  call changemaxcontot
                endif
                if (ncf.lt.ncfg) then
                  do k = ncontot,n1con(ncf+1),-1
                    ncvarindcfg(k+1) = ncvarindcfg(k)
                    ncvartypcfg(k+1) = ncvartypcfg(k)
                    ncfixindcfg(k+1) = ncfixindcfg(k)
                    ncfixtypcfg(k+1) = ncfixtypcfg(k)
                    concocfg(k+1) = concocfg(k)
                    conaddcfg(k+1) = conaddcfg(k)
                    nconcfg(k+1) = nconcfg(k)
                  enddo
                  do k = ncf+1,ncfg
                    n1con(k) = n1con(k) + 1
                  enddo
                endif
                ncontot = ncontot + 1
                ltmploc(2) = .true.
                ncon = ncon + 1
                ncvarindcfg(ncfst+ncon) = indi
                ncvartypcfg(ncfst+ncon) = iopt_yf
                ncfixindcfg(ncfst+ncon) = indi
                ncfixtypcfg(ncfst+ncon) = iopt_xf
                concocfg(ncfst+ncon) = xposs(ii)*ryp
                conaddcfg(ncfst+ncon) = xorig - yorig*concocfg(ncfst+ncon)
                nconcfg(ncfst+ncon) = ncf
                lfound = .true.
              else
                if (ncontot.ge.maxcontot) then
                  maxcontot = ncontot + 10
                  call changemaxcontot
                endif
                if (ncf.lt.ncfg) then
                  do k = ncontot,n1con(ncf+1),-1
                    ncvarindcfg(k+1) = ncvarindcfg(k)
                    ncvartypcfg(k+1) = ncvartypcfg(k)
                    ncfixindcfg(k+1) = ncfixindcfg(k)
                    ncfixtypcfg(k+1) = ncfixtypcfg(k)
                    concocfg(k+1) = concocfg(k)
                    conaddcfg(k+1) = conaddcfg(k)
                    nconcfg(k+1) = nconcfg(k)
                  enddo
                  do k = ncf+1,ncfg
                    n1con(k) = n1con(k) + 1
                  enddo
                endif
                ncontot = ncontot + 1
                ltmploc(3) = .true.
                ncon = ncon + 1
                ncvarindcfg(ncfst+ncon) = indi
                ncvartypcfg(ncfst+ncon) = iopt_zf
                ncfixindcfg(ncfst+ncon) = indi
                ncfixtypcfg(ncfst+ncon) = iopt_xf
                concocfg(ncfst+ncon) = xposs(ii)*rzp
                conaddcfg(ncfst+ncon) = xorig - zorig*conco(ncfst+ncon)
                nconcfg(ncfst+ncon) = ncf
                lfound = .true.
              endif
            elseif (lfoundxz.and.loptx(indi).and.loptz(indi)) then
              ncv = ncvarind(nxz1)
              if (ncv.eq.indi.and.ncvartyp(nxz1).eq.iopt_xf) then
                if (ncontot.ge.maxcontot) then
                  maxcontot = ncontot + 10
                  call changemaxcontot
                endif
                if (ncf.lt.ncfg) then
                  do k = ncontot,n1con(ncf+1),-1
                    ncvarindcfg(k+1) = ncvarindcfg(k)
                    ncvartypcfg(k+1) = ncvartypcfg(k)
                    ncfixindcfg(k+1) = ncfixindcfg(k)
                    ncfixtypcfg(k+1) = ncfixtypcfg(k)
                    concocfg(k+1) = concocfg(k)
                    conaddcfg(k+1) = conaddcfg(k)
                    nconcfg(k+1) = nconcfg(k)
                  enddo
                  do k = ncf+1,ncfg
                    n1con(k) = n1con(k) + 1
                  enddo
                endif
                ncontot = ncontot + 1
                ltmploc(1) = .true.
                ncon = ncon + 1
                ncvarindcfg(ncfst+ncon) = indi
                ncvartypcfg(ncfst+ncon) = iopt_xf
                ncfixindcfg(ncfst+ncon) = indi
                ncfixtypcfg(ncfst+ncon) = iopt_yf
                concocfg(ncfst+ncon) = yposs(ii)*rxp
                conaddcfg(ncfst+ncon) = yorig - xorig*concocfg(ncfst+ncon)
                nconcfg(ncfst+ncon) = ncf
                lfound = .true.
              else
                if (ncontot.ge.maxcontot) then
                  maxcontot = ncontot + 10
                  call changemaxcontot
                endif
                if (ncf.lt.ncfg) then
                  do k = ncontot,n1con(ncf+1),-1
                    ncvarindcfg(k+1) = ncvarindcfg(k)
                    ncvartypcfg(k+1) = ncvartypcfg(k)
                    ncfixindcfg(k+1) = ncfixindcfg(k)
                    ncfixtypcfg(k+1) = ncfixtypcfg(k)
                    concocfg(k+1) = concocfg(k)
                    conaddcfg(k+1) = conaddcfg(k)
                    nconcfg(k+1) = nconcfg(k)
                  enddo
                  do k = ncf+1,ncfg
                    n1con(k) = n1con(k) + 1
                  enddo
                endif
                ncontot = ncontot + 1
                ltmploc(3) = .true.
                ncon = ncon + 1
                ncvarindcfg(ncfst+ncon) = indi
                ncvartypcfg(ncfst+ncon) = iopt_zf
                ncfixindcfg(ncfst+ncon) = indi
                ncfixtypcfg(ncfst+ncon) = iopt_yf
                concocfg(ncfst+ncon) = yposs(ii)*rzp
                conaddcfg(ncfst+ncon) = yorig - zorig*concocfg(ncfst+ncon)
                nconcfg(ncfst+ncon) = ncf
                lfound = .true.
              endif
            endif
!
!  No constraints set already
!
            if (.not.lfound) then
!
!  Move constraints along to make space for new one
!
              if (ncontot+2.gt.maxcontot) then
                maxcontot = ncontot + 10
                call changemaxcontot
              endif
              if (ncf.lt.ncfg) then
                do k = ncontot,n1con(ncf+1),-1
                  ncvarindcfg(k+2) = ncvarindcfg(k)
                  ncvartypcfg(k+2) = ncvartypcfg(k)
                  ncfixindcfg(k+2) = ncfixindcfg(k)
                  ncfixtypcfg(k+2) = ncfixtypcfg(k)
                  concocfg(k+2) = concocfg(k)
                  conaddcfg(k+2) = conaddcfg(k)
                  nconcfg(k+2) = nconcfg(k)
                enddo
                do k = ncf+1,ncfg
                  n1con(k) = n1con(k) + 2
                enddo
              endif
              ncontot = ncontot + 2
              ltmploc(1) = .true.
              ncon = ncon + 1
              ncvarindcfg(ncfst+ncon) = indi
              ncvartypcfg(ncfst+ncon) = iopt_xf
              ncfixindcfg(ncfst+ncon) = indi
              ncfixtypcfg(ncfst+ncon) = iopt_yf
              concocfg(ncfst+ncon) = yposs(ii)*rxp
              conaddcfg(ncfst+ncon) = yorig - xorig*concocfg(ncfst+ncon)
              nconcfg(ncfst+ncon) = ncf
              ncon = ncon+1
              ncvarindcfg(ncfst+ncon) = indi
              ncvartypcfg(ncfst+ncon) = iopt_xf
              ncfixindcfg(ncfst+ncon) = indi
              ncfixtypcfg(ncfst+ncon) = iopt_zf
              concocfg(ncfst+ncon) = zposs(ii)*rxp
              conaddcfg(ncfst+ncon) = zorig - xorig*concocfg(ncfst+ncon)
              nconcfg(ncfst+ncon) = ncf
            endif
          elseif (lxposs.and.lyposs) then
!********************
!  x y combination  *
!********************
            rxp = 1.0_dp/xposs(ii)
            ryp = 1.0_dp/yposs(ii)
!
!  Check to see if any constraints already exist for these variables
!
            lfound = .false.
            do k = 1,ncon
              nvk = ncvarindcfg(ncfst+k)
              nfk = ncfixindcfg(ncfst+k)
              if (ncvartypcfg(ncfst+k).eq.iopt_xf.and.ncfixtypcfg(ncfst+k).eq.iopt_yf.and. &
                  nvk.eq.indi.and.nfk.eq.indi) then
                lfound = .true.
                ltmploc(1) = .true.
                coeff = yposs(ii)*rxp
                diff = abs(coeff-concocfg(ncfst+k))
                if (diff.gt.1.0d-4) then
                  nwarn = nwarn + 1
                  if (ioproc) then
                    call itow(conno,k,3_i4)
                    call outerror('Coefficient for constraint '//conno//' is being reset ',0_i4)
                  endif
                endif
                concocfg(ncfst+k) = coeff
                coeff = yorig - xorig*coeff
                diff = abs(coeff-conaddcfg(ncfst+k))
                if (diff.gt.1.0d-4) then
                  nwarn = nwarn + 1
                  if (ioproc) then
                    call itow(conno,k,3_i4)
                    call outerror('Coefficient for constraint '//conno//' is being reset ',0_i4)
                  endif
                endif
                conaddcfg(ncfst+k) = coeff
              elseif (ncvartypcfg(ncfst+k).eq.iopt_yf.and.ncfixtypcfg(ncfst+k).eq.iopt_xf.and. &
                  nvk.eq.indi.and.nfk.eq.indi) then
                lfound = .true.
                ltmploc(2) = .true.
                coeff = xposs(ii)*ryp
                diff = abs(coeff-concocfg(ncfst+k))
                if (diff.gt.1.0d-4) then
                  nwarn = nwarn + 1
                  if (ioproc) then
                    call itow(conno,k,3_i4)
                    call outerror('Coefficient for constraint '//conno//' is being reset ',0_i4)
                  endif
                endif
                concocfg(ncfst+k) = coeff
                coeff = xorig - yorig*coeff
                diff = abs(coeff-conaddcfg(ncfst+k))
                if (diff.gt.1.0d-4) then
                  nwarn = nwarn + 1
                  if (ioproc) then
                    call itow(conno,k,3_i4)
                    call outerror('Coefficient for constraint '//conno//' is being reset ',0_i4)
                  endif
                endif
                conaddcfg(ncfst+k) = coeff
              endif
            enddo
            if (.not.lfound) then
!
!  Move constraints along to make space for new one
!
              if (ncontot.ge.maxcontot) then
                maxcontot = ncontot + 10
                call changemaxcontot
              endif
              if (ncf.lt.ncfg) then
                do k = ncontot,n1con(ncf+1),-1
                  ncvarindcfg(k+1) = ncvarindcfg(k)
                  ncvartypcfg(k+1) = ncvartypcfg(k)
                  ncfixindcfg(k+1) = ncfixindcfg(k)
                  ncfixtypcfg(k+1) = ncfixtypcfg(k)
                  concocfg(k+1) = concocfg(k)
                  conaddcfg(k+1) = conaddcfg(k)
                  nconcfg(k+1) = nconcfg(k)
                enddo
                do k = ncf+1,ncfg
                  n1con(k) = n1con(k) + 1
                enddo
              endif
              ncontot = ncontot + 1
              ltmploc(1) = .true.
              ncon = ncon + 1
              ncvarindcfg(ncfst+ncon) = indi
              ncvartypcfg(ncfst+ncon) = iopt_xf
              ncfixindcfg(ncfst+ncon) = indi
              ncfixtypcfg(ncfst+ncon) = iopt_yf
              concocfg(ncfst+ncon) = yposs(ii)*rxp
              conaddcfg(ncfst+ncon) = yorig - xorig*concocfg(ncfst+ncon)
              nconcfg(ncfst+ncon) = ncf
            endif
          elseif (lyposs.and.lzposs) then
!********************
!  y z combination  *
!********************
            ryp = 1.0_dp/yposs(ii)
            rzp = 1.0_dp/zposs(ii)
!
!  Check to see if any constraints already exist for these variables
!
            lfound = .false.
            do k = 1,ncon
              nvk = ncvarindcfg(ncfst+k)
              nfk = ncfixindcfg(ncfst+k)
              if (ncvartypcfg(ncfst+k).eq.iopt_yf.and.ncfixtypcfg(ncfst+k).eq.iopt_zf.and. &
                  nvk.eq.indi.and.nfk.eq.indi) then
                lfound = .true.
                ltmploc(2) = .true.
                coeff = zposs(ii)*ryp
                diff = abs(coeff-concocfg(ncfst+k))
                if (diff.gt.1.0d-4) then
                  nwarn = nwarn + 1
                  if (ioproc) then
                    call itow(conno,k,3_i4)
                    call outerror('Coefficient for constraint '//conno//' is being reset ',0_i4)
                  endif
                endif
                concocfg(ncfst+k) = coeff
                coeff = zorig - yorig*coeff
                diff = abs(coeff-conaddcfg(ncfst+k))
                if (diff.gt.1.0d-4) then
                  nwarn = nwarn + 1
                  if (ioproc) then
                    call itow(conno,k,3_i4)
                    call outerror('Coefficient for constraint '//conno//' is being reset ',0_i4)
                  endif
                endif
                conaddcfg(ncfst+k) = coeff
              elseif (ncvartypcfg(ncfst+k).eq.iopt_zf.and.ncfixtypcfg(ncfst+k).eq.iopt_yf.and. &
                  nvk.eq.indi.and.nfk.eq.indi) then
                lfound = .true.
                ltmploc(3) = .true.
                coeff = yposs(ii)*rzp
                diff = abs(coeff-concocfg(ncfst+k))
                if (diff.gt.1.0d-4) then
                  nwarn = nwarn + 1
                  if (ioproc) then
                    call itow(conno,k,3_i4)
                    call outerror('Coefficient for constraint '//conno//' is being reset ',0_i4)
                  endif
                endif
                concocfg(ncfst+k) = coeff
                coeff = yorig - zorig*coeff
                diff = abs(coeff-conaddcfg(ncfst+k))
                if (diff.gt.1.0d-4) then
                  nwarn = nwarn + 1
                  if (ioproc) then
                    call itow(conno,k,3_i4)
                    call outerror('Coefficient for constraint '//conno//' is being reset ',0_i4)
                  endif
                endif
                conaddcfg(ncfst+k) = coeff
              endif
            enddo
            if (.not.lfound) then
!
!  Move constraints along to make space for new one
!
              if (ncontot.ge.maxcontot) then
                maxcontot = ncontot + 10
                call changemaxcontot
              endif
              if (ncf.lt.ncfg) then
                do k = ncontot,n1con(ncf+1),-1
                  ncvarindcfg(k+1) = ncvarindcfg(k)
                  ncvartypcfg(k+1) = ncvartypcfg(k)
                  ncfixindcfg(k+1) = ncfixindcfg(k)
                  ncfixtypcfg(k+1) = ncfixtypcfg(k)
                  concocfg(k+1) = concocfg(k)
                  conaddcfg(k+1) = conaddcfg(k)
                  nconcfg(k+1) = nconcfg(k)
                enddo
                do k = ncf+1,ncfg
                  n1con(k) = n1con(k) + 1
                enddo
              endif
              ncontot = ncontot + 1
              ltmploc(2) = .true.
              ncon = ncon + 1
              ncvarindcfg(ncfst+ncon) = indi
              ncvartypcfg(ncfst+ncon) = iopt_yf
              ncfixindcfg(ncfst+ncon) = indi
              ncfixtypcfg(ncfst+ncon) = iopt_zf
              concocfg(ncfst+ncon) = zposs(ii)*ryp
              conaddcfg(ncfst+ncon) = zorig - yorig*concocfg(ncfst+ncon)
              nconcfg(ncfst+ncon) = ncf
            endif
          elseif (lxposs.and.lzposs) then
!********************
!  x z combination  *
!********************
            rxp = 1.0_dp/xposs(ii)
            rzp = 1.0_dp/zposs(ii)
!
!  Check to see if any constraints already exist for these variables
!
            lfound = .false.
            do k = 1,ncon
              nvk = ncvarindcfg(ncfst+k)
              nfk = ncfixindcfg(ncfst+k)
              if (ncvartypcfg(ncfst+k).eq.iopt_xf.and.ncfixtypcfg(ncfst+k).eq.iopt_zf.and. &
                  nvk.eq.indi.and.nfk.eq.indi) then
                lfound = .true.
                ltmploc(1) = .true.
                coeff = zposs(ii)*rxp
                diff = abs(coeff-concocfg(ncfst+k))
                if (diff.gt.1.0d-4) then
                  nwarn = nwarn + 1
                  if (ioproc) then
                    call itow(conno,k,3_i4)
                    call outerror('Coefficient for constraint '//conno//' is being reset ',0_i4)
                  endif
                endif
                concocfg(ncfst+k) = coeff
                coeff = zorig - xorig*coeff
                diff = abs(coeff-conaddcfg(ncfst+k))
                if (diff.gt.1.0d-4) then
                  nwarn = nwarn + 1
                  if (ioproc) then
                    call itow(conno,k,3_i4)
                    call outerror('Coefficient for constraint '//conno//' is being reset ',0_i4)
                  endif
                endif
                conaddcfg(ncfst+k) = coeff
              elseif (ncvartypcfg(ncfst+k).eq.iopt_zf.and.ncfixtypcfg(ncfst+k).eq.iopt_xf.and. &
                  nvk.eq.indi.and.nfk.eq.indi) then
                lfound = .true.
                ltmploc(3) = .true.
                coeff = xposs(ii)*rzp
                diff = abs(coeff-concocfg(ncfst+k))
                if (diff.gt.1.0d-4) then
                  nwarn = nwarn + 1
                  if (ioproc) then
                    call itow(conno,k,3_i4)
                    call outerror('Coefficient for constraint '//conno//' is being reset ',0_i4)
                  endif
                endif
                concocfg(ncfst+k) = coeff
                coeff = xorig - zorig*coeff
                diff = abs(coeff-conaddcfg(ncfst+k))
                if (diff.gt.1.0d-4) then
                  nwarn = nwarn + 1
                  if (ioproc) then
                    call itow(conno,k,3_i4)
                    call outerror('Coefficient for constraint '//conno//' is being reset ',0_i4)
                  endif
                endif
                conaddcfg(ncfst+k) = coeff
              endif
            enddo
            if (.not.lfound) then
!
!  Move constraints along to make space for new one
!
              if (ncontot.ge.maxcontot) then
                maxcontot = ncontot + 10
                call changemaxcontot
              endif
              if (ncf.lt.ncfg) then
                do k = ncontot,n1con(ncf+1),-1
                  ncvarindcfg(k+1) = ncvarindcfg(k)
                  ncvartypcfg(k+1) = ncvartypcfg(k)
                  ncfixindcfg(k+1) = ncfixindcfg(k)
                  ncfixtypcfg(k+1) = ncfixtypcfg(k)
                  concocfg(k+1) = concocfg(k)
                  conaddcfg(k+1) = conaddcfg(k)
                  nconcfg(k+1) = nconcfg(k)
                enddo
                do k = ncf+1,ncfg
                  n1con(k) = n1con(k) + 1
                enddo
              endif
              ncontot = ncontot + 1
              ltmploc(1) = .true.
              ncon = ncon + 1
              ncvarindcfg(ncfst+ncon) = indi
              ncvartypcfg(ncfst+ncon) = iopt_xf
              ncfixindcfg(ncfst+ncon) = indi
              ncfixtypcfg(ncfst+ncon) = iopt_zf
              concocfg(ncfst+ncon) = zposs(ii)*rxp
              conaddcfg(ncfst+ncon) = zorig - xorig*concocfg(ncfst+ncon)
              nconcfg(ncfst+ncon) = ncf
            endif
          endif
        endif
      enddo npossloop
    endif
    if (lcheck) then
!
!  Compare flags against input values if in checking mode
!
      if (loptx(indi).and..not.ltmploc(1)) then
        call itow(atomno,na,6_i4)
        call outerror('Badly defined flags for '//'atom '//atomno,0_i4)
        call stopnow('special')
      endif
      if (lopty(indi).and..not.ltmploc(2)) then
        call itow(atomno,na,6_i4)
        call outerror('Badly defined flags for '//'atom '//atomno,0_i4)
        call stopnow('special')
      endif
      if (loptz(indi).and..not.ltmploc(3)) then
        call itow(atomno,na,6_i4)
        call outerror('Badly defined flags for '//'atom '//atomno,0_i4)
        call stopnow('special')
      endif
    else
!
!  Transfer flags to main array - only do this if the flag has initially been 
!  set to true though so that user flags can be validated.
!
      if (loptx(indi)) loptx(indi) = ltmploc(1)
      if (lopty(indi)) lopty(indi) = ltmploc(2)
      if (loptz(indi)) loptz(indi) = ltmploc(3)
    endif
!
!  End of loop over atoms
!
  enddo naloop
  time2 = g_cpu_time()
  tsym = tsym + time2 - time1
#ifdef TRACE
  call trace_out('special')
#endif
!
  return
  end
