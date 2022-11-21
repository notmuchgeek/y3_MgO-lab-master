  subroutine setregiontrans(loptx,lopty,loptz,loptm)
!
!  This routine creates constraints that are needed to perform
!  the rigid translation of regions.
!
!  lopt = array of logical flags according to whether a
!         parameter can be varied or not. Passed from
!         setcfg.
!
!   5/03 Created
!   9/03 Modified to use logical pointer as to whether region
!        is to be rigid or not.
!   3/19 Constraint arrays changed to have index and type
!   3/19 ltmp changed for individual arrays for x, y and z
!   4/19 Indices and types corrected
!   7/20 Rigid molecule modifications added including check for
!        molecules being split over multiple regions
!   7/20 Symmetry handling added for rigid molecule modifications
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
!  Julian Gale, CIC, Curtin University, July 2020
!
  use configurations
  use control
  use current
  use molecule
  use optimisation
  implicit none
!
!  Passed variables
!
  logical,     intent(inout) :: loptm(6,*)
  logical,     intent(inout) :: loptx(*)
  logical,     intent(inout) :: lopty(*)
  logical,     intent(inout) :: loptz(*)
!
!  Local variables
!
  integer(i4)                :: i
  integer(i4)                :: ifirst
  integer(i4)                :: ii
  integer(i4)                :: k
  integer(i4)                :: natomthisregion
  integer(i4)                :: ndirection
  integer(i4)                :: nnewcon
  integer(i4)                :: nr
  logical                    :: lfirstatom
  logical                    :: lfirstisatom
!
!  For rigid molecules check that they are not split over multiple regions since this will cause problems
!
  if (lrigid) then
    do ii = 1,nmolasym
      i = nmolasymptr(1,ii)
      nr = nregionno(nsft+i)
      do k = 2,nmolasymno(ii)
        i = nmolasymptr(k,ii)
        if (nregionno(nsft+i).ne.nr) then
          call outerror('rigid molecules cannot be split over multiple regions',0_i4)
          call stopnow('setregiontrans')
        endif
      enddo
    enddo
  endif
!**********************
!  Loop over regions  *
!**********************
  do nr = 1,nregions(ncf)
    if (lregionrigid(nr,ncf)) then
!
!  Count number of directions of motion
!
      ndirection = 0
      if (lopfreg(3*(nr-1)+1,ncf)) ndirection = ndirection + 1
      if (lopfreg(3*(nr-1)+2,ncf)) ndirection = ndirection + 2
      if (lopfreg(3*(nr-1)+3,ncf)) ndirection = ndirection + 3
      if (ndirection.gt.0) then
!
!  Find atoms for this region and count
!
        natomthisregion = 0
        if (lrigid) then
          do ii = 1,nasymnomol
            i = nasymnomolptr(ii)
            if (nregionno(nsft+i).eq.nr) then
              natomthisregion = natomthisregion + 1
            endif
          enddo
          do ii = 1,nmolasym
            i = nmolasymptr(1,ii)
            if (nregionno(nsft+i).eq.nr) then
              natomthisregion = natomthisregion + 1
            endif
          enddo
        else
          do i = 1,nasym
            if (nregionno(nsft+i).eq.nr) then
              natomthisregion = natomthisregion + 1
            endif
          enddo
        endif
!
!  Check that there is space in constraint arrays
!
        nnewcon = ndirection*(natomthisregion - 1)
        if (ncontot+nnewcon.ge.maxcontot) then
          maxcontot = ncontot + nnewcon
          call changemaxcontot
        endif
!********************************************
!  Move data to make space for constraints  *
!********************************************
        if (ncf.lt.ncfg) then
          do k = ncontot,n1con(ncf+1),-1
            ncvarindcfg(k+nnewcon) = ncvarindcfg(k)
            ncvartypcfg(k+nnewcon) = ncvartypcfg(k)
            ncfixindcfg(k+nnewcon) = ncfixindcfg(k)
            ncfixtypcfg(k+nnewcon) = ncfixtypcfg(k)
            concocfg(k+nnewcon) = concocfg(k)
            nconcfg(k+nnewcon) = nconcfg(k)
            conaddcfg(k+nnewcon) = conaddcfg(k)
          enddo
          do k = ncf+1,ncfg
            n1con(k) = n1con(k) + nnewcon
          enddo
        endif
!
!  Find atoms/molecules for this region
!
        lfirstatom = .true.
        if (lrigid) then
!--------------------
!  Rigid molecules  |
!--------------------
!
!  First loop over atoms not in molecules
!
          do ii = 1,nasymnomol
            i = nasymnomolptr(ii)
            if (nregionno(nsft+i).eq.nr) then
!
!  Is this the first atom for the region?
!
              if (lfirstatom) then
                ifirst = i
                lfirstatom = .false.
                lfirstisatom = .true.
              endif
              if (i.ne.ifirst) then
!*********************************************
!  Add constraint to atom if not first atom  *
!*********************************************
!
!  X direction
!
                if (lopfreg(3*(nr-1)+1,ncf)) then
                  ncon = ncon + 1
                  ncontot = ncontot + 1
                  loptx(i) = .false.
                  ncvarindcfg(ncon) = ifirst
                  ncvartypcfg(ncon) = iopt_xf
                  ncfixindcfg(ncon) = i
                  ncfixtypcfg(ncon) = iopt_xf
                  concocfg(ncon) = 1.0_dp
                  conaddcfg(ncon) = xcfg(nsft+i) - xcfg(nsft+ifirst)
                  nconcfg(ncon) = ncf
                endif
!
!  Y direction
!
                if (lopfreg(3*(nr-1)+2,ncf)) then
                  ncon = ncon + 1
                  ncontot = ncontot + 1
                  lopty(i) = .false.
                  ncvarindcfg(ncon) = ifirst
                  ncvartypcfg(ncon) = iopt_yf
                  ncfixindcfg(ncon) = i
                  ncfixtypcfg(ncon) = iopt_yf
                  concocfg(ncon) = 1.0_dp
                  conaddcfg(ncon) = ycfg(nsft+i) - ycfg(nsft+ifirst)
                  nconcfg(ncon) = ncf
                endif
!
!  Z direction
!
                if (lopfreg(3*(nr-1)+3,ncf)) then
                  ncon = ncon + 1
                  ncontot = ncontot + 1
                  loptz(i) = .false.
                  ncvarindcfg(ncon) = ifirst
                  ncvartypcfg(ncon) = iopt_zf
                  ncfixindcfg(ncon) = i
                  ncfixtypcfg(ncon) = iopt_zf
                  concocfg(ncon) = 1.0_dp
                  conaddcfg(ncon) = zcfg(nsft+i) - zcfg(nsft+ifirst)
                  nconcfg(ncon) = ncf
                endif
              endif
            endif
          enddo
!
!  Second loop over rigid molecules
!
          do ii = 1,nmolasym
            i = nmolasymptr(1,ii)
            if (nregionno(nsft+i).eq.nr) then
!
!  Is this the first molecule for the region?
!
              if (lfirstatom) then
                ifirst = ii
                lfirstatom = .false.
                lfirstisatom = .false.
              endif
              if (ii.ne.ifirst.or.lfirstisatom) then
!*********************************************
!  Add constraint to atom if not first atom  *
!*********************************************
!
!  X direction
!
                if (lopfreg(3*(nr-1)+1,ncf)) then
                  ncon = ncon + 1
                  ncontot = ncontot + 1
                  loptm(1,ii) = .false.
                  ncvarindcfg(ncon) = ifirst
                  if (lfirstisatom) then
                    ncvartypcfg(ncon) = iopt_xf
                  else
                    ncvartypcfg(ncon) = iopt_xcom
                  endif
                  ncfixindcfg(ncon) = ii
                  ncfixtypcfg(ncon) = iopt_xcom
                  concocfg(ncon) = 1.0_dp
                  if (lfirstisatom) then
                    conaddcfg(ncon) = molcom(1,ii) - xcfg(nsft+ifirst)
                  else
                    conaddcfg(ncon) = molcom(1,ii) - molcom(1,ifirst)
                  endif
                  nconcfg(ncon) = ncf
                endif
!
!  Y direction
!
                if (lopfreg(3*(nr-1)+2,ncf)) then
                  ncon = ncon + 1
                  ncontot = ncontot + 1
                  loptm(2,ii) = .false.
                  ncvarindcfg(ncon) = ifirst
                  if (lfirstisatom) then
                    ncvartypcfg(ncon) = iopt_yf
                  else
                    ncvartypcfg(ncon) = iopt_ycom
                  endif
                  ncfixindcfg(ncon) = ii
                  ncfixtypcfg(ncon) = iopt_ycom
                  concocfg(ncon) = 1.0_dp
                  if (lfirstisatom) then
                    conaddcfg(ncon) = molcom(2,ii) - ycfg(nsft+ifirst)
                  else
                    conaddcfg(ncon) = molcom(2,ii) - molcom(2,ifirst)
                  endif
                  nconcfg(ncon) = ncf
                endif
!
!  Z direction
!
                if (lopfreg(3*(nr-1)+3,ncf)) then
                  ncon = ncon + 1
                  ncontot = ncontot + 1
                  loptm(3,ii) = .false.
                  ncvarindcfg(ncon) = ifirst
                  if (lfirstisatom) then
                    ncvartypcfg(ncon) = iopt_zf
                  else
                    ncvartypcfg(ncon) = iopt_zcom
                  endif
                  ncfixindcfg(ncon) = ii
                  ncfixtypcfg(ncon) = iopt_zcom
                  concocfg(ncon) = 1.0_dp
                  if (lfirstisatom) then
                    conaddcfg(ncon) = molcom(3,ii) - zcfg(nsft+ifirst)
                  else
                    conaddcfg(ncon) = molcom(3,ii) - molcom(3,ifirst)
                  endif
                  nconcfg(ncon) = ncf
                endif
              endif
            endif
          enddo
        else
!-----------------------
!  No rigid molecules  |
!-----------------------
          do i = 1,nasym
            if (nregionno(nsft+i).eq.nr) then
!
!  Is this the first atom for the region?
!
              if (lfirstatom) then
                ifirst = i
                lfirstatom = .false.
              endif
              if (i.ne.ifirst) then
!*********************************************
!  Add constraint to atom if not first atom  *
!*********************************************
!
!  X direction
!
                if (lopfreg(3*(nr-1)+1,ncf)) then
                  ncon = ncon + 1
                  ncontot = ncontot + 1
                  loptx(i) = .false.
                  ncvarindcfg(ncon) = ifirst
                  ncvartypcfg(ncon) = iopt_xf
                  ncfixindcfg(ncon) = i
                  ncfixtypcfg(ncon) = iopt_xf
                  concocfg(ncon) = 1.0_dp
                  conaddcfg(ncon) = xcfg(nsft+i) - xcfg(nsft+ifirst)
                  nconcfg(ncon) = ncf
                endif
!
!  Y direction
!
                if (lopfreg(3*(nr-1)+2,ncf)) then
                  ncon = ncon + 1
                  ncontot = ncontot + 1
                  lopty(i) = .false.
                  ncvarindcfg(ncon) = ifirst
                  ncvartypcfg(ncon) = iopt_yf
                  ncfixindcfg(ncon) = i
                  ncfixtypcfg(ncon) = iopt_yf
                  concocfg(ncon) = 1.0_dp
                  conaddcfg(ncon) = ycfg(nsft+i) - ycfg(nsft+ifirst)
                  nconcfg(ncon) = ncf
                endif
!
!  Z direction
!
                if (lopfreg(3*(nr-1)+3,ncf)) then
                  ncon = ncon + 1
                  ncontot = ncontot + 1
                  loptz(i) = .false.
                  ncvarindcfg(ncon) = ifirst
                  ncvartypcfg(ncon) = iopt_zf
                  ncfixindcfg(ncon) = i
                  ncfixtypcfg(ncon) = iopt_zf
                  concocfg(ncon) = 1.0_dp
                  conaddcfg(ncon) = zcfg(nsft+i) - zcfg(nsft+ifirst)
                  nconcfg(ncon) = ncf
                endif
              endif
            endif
          enddo
        endif
      endif
    endif
  enddo
!
  return
  end
