  subroutine gabcfg(mc,mcfg)
!
!  Subroutine for obtaining possible crystal structure
!  Structures to be optimised are stored in xbest(,)
!
!  11/06 lfirst argument added to equpos call
!   6/12 File handling changed to allow for spaces in names
!   2/18 Trace added
!   8/18 Adding 1 to strains 1-3 removed
!   2/19 x0 removed
!   3/19 iopt replaced by ioptindex and iopttype
!   3/19 Constraint handling now moved to subroutine
!
!  Scott Woodley, R.I.G.B., June 1997
!  Julian Gale, CIC, Curtin University, March 2019
!
  use configurations
  use control
  use current
  use dump
  use gulp_files
  use gaconf,        only : xbest, ithbest
  use genetic
  use iochannels
  use optimisation
  use parallel
  use symmetry
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4)                            :: mc
  integer(i4)                            :: mcfg
!
!  Local variables
!
  integer(i4)                            :: i
  integer(i4)                            :: i0
  integer(i4)                            :: iend
  integer(i4)                            :: mmc
  integer(i4)                            :: nc
  integer(i4),                      save :: nend
  logical                                :: lstuck
  real(dp),                         save :: x, y, z
#ifdef TRACE
  call trace_in('gabcfg')
#endif
!
  lstuck = (index(keyword,'hhhh').ne.0)
!
  if (ioproc) then
    write(ioout,'(''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'')')
    if (.not.lstuck) then
      write(ioout,'(''################################################################################'')')
    endif
    write(ioout,*)'         ',mc-1,' structures optimised out of ',mcfg
    if (.not.lstuck) then
      write(ioout,'(''################################################################################'')')
      write(ioout,'(''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'')')
    endif
  endif

  if ((lrelax.or.lcomp).and.mcfg.ne.1) then
    call outerror('problem in gabcfg',0_i4)
    call stopnow('gabcfg')
  endif
!
!  Select appropriate element
!
  if (lgaconjg) then
    mmc = ithbest(mc)
  else
    mmc = mc
  endif
!
!  Set initial atom (which is not randomised) at (0,0,0)?
!
  if (mc.eq.1) then
    x = xafrac(1)
    y = yafrac(1)
    z = zafrac(1)
  else
    xafrac(1) = x
    yafrac(1) = y
    zafrac(1) = z
  endif
!
!  Substitute parameters into place
!
  do i = 1,nvar
    if (iopttypecfg(i+nfst).eq.iopt_cell) then
      if (ioptindexcfg(i+nfst).eq.1) then
        a = xbest(i,mmc)
      elseif (ioptindexcfg(i+nfst).eq.2) then
        b = xbest(i,mmc)
      elseif (ioptindexcfg(i+nfst).eq.3) then
        c = xbest(i,mmc)
      elseif (ioptindexcfg(i+nfst).eq.4) then
        alpha = xbest(i,mmc)
      elseif (ioptindexcfg(i+nfst).eq.5) then
        beta = xbest(i,mmc)
      elseif (ioptindexcfg(i+nfst).eq.6) then
        gamma = xbest(i,mmc)
      endif
    elseif (iopttypecfg(i+nfst).eq.iopt_strain) then
      strain(ioptindexcfg(i+nfst)) = xbest(i,mmc)
    elseif (iopttypecfg(i+nfst).eq.iopt_xf) then
      xafrac(ioptindexcfg(i+nfst)) = xbest(i,mmc)
    elseif (iopttypecfg(i+nfst).eq.iopt_yf) then
      yafrac(ioptindexcfg(i+nfst)) = xbest(i,mmc)
    elseif (iopttypecfg(i+nfst).eq.iopt_zf) then
      zafrac(ioptindexcfg(i+nfst)) = xbest(i,mmc)
    elseif (iopttypecfg(i+nfst).eq.iopt_radius) then
      rada(ioptindexcfg(i+nfst)) = xbest(i,mmc)
    endif
  enddo
!
!  Change output filename for insight package
!
  if (.not.lgaconjg) then
    if (lxtl) then
      if (ncf.eq.1.and.mc.eq.1) then
        nend = index(xtlfile,'.xtl')
        if (nend.gt.1) then
          if (nend.gt.20) then
            nend = 20
            xtlfile(20:30) = '.xtl       '
          endif
          xtlfile(nend-1:nend-1) = char(48+ncf)
        elseif (index(names(ncf),' ').ne.1) then
          xtlfile(1:4) = names(ncf)(1:4)
          xtlfile(5:8) = '.xtl'
        else
          xtlfile(1:3) = 'str'
          xtlfile(4:4) = char(48+ncf)
          xtlfile(5:8) = '.xtl'
        endif
      elseif (mc.eq.1) then
        if (nend.gt.0) then
          xtlfile(nend-1:nend-1) = char(48+ncf)
        elseif (index(names(ncf),' ').ne.1) then
          xtlfile(1:4) = names(ncf)(1:4)
          xtlfile(5:8) = '.xtl'
        else
          xtlfile(1:3) = 'str'
          xtlfile(4:4) = char(48+ncf)
          xtlfile(5:8) = '.xtl'
        endif
      endif
      iend = index(xtlfile,'.xtl')
      if (mc.ne.1) iend = iend - 2
      if (iend.gt.25) iend = iend - 5
      if (mc.gt.30) then
        if (ioproc) write(ioout,*)'Too many structures to optimise!'
        call stopnow('gabcfg')
      elseif (mc.gt.20) then
        xtlfile(iend:iend) = '2'
        nc = mc - 20
      elseif (mc.gt.10)  then
        xtlfile(iend:iend) = '1'
        nc = mc - 10
      else
        xtlfile(iend:iend) = '0'
        nc = mc
      endif
      iend = iend + 1
      i0 = ichar('0')
      i0 = i0 - 1
      xtlfile(iend:iend) = char(nc+i0)
      xtlfile(iend+1:iend+4) = '.xtl'
    endif
!
!  Similarly for (individual restart) dump files
!
    if (idump.eq.12) then
      call endstring(dfile,len(dfile),iend)
      if (mc.ne.1) then
        iend = iend - 2
      endif
      if (iend.gt.25) iend = iend - 5
      if (mc.gt.30) then
        if (ioproc) write(ioout,*)'Too many structures to optimise!'
        call stopnow('gabcfg')
      elseif (mc.gt.20) then
        dfile(iend:iend) = '2'
        nc = mc - 20
      elseif (mc.gt.10)  then
        dfile(iend:iend) = '1'
        nc = mc - 10
      else
        dfile(iend:iend) = '0'
        nc = mc
      endif
      iend = iend + 1
      i0 = ichar('0')
      i0 = i0 - 1
      dfile(iend:iend)=char(nc+i0)
    endif
  endif
!**********************
!  Apply constraints  *
!**********************
  if (ncon.gt.0) then
    call applyconstraints
  endif
  if (ncell.gt.0) then
!
!  Apply strains
!
    do i = 1,3
      rv(1,i) = rvcfg(1,i,ncf)
      rv(2,i) = rvcfg(2,i,ncf)
      rv(3,i) = rvcfg(3,i,ncf)
    enddo
    call strain3D(strain,rv)
    call uncell3D(rv,a,b,c,alpha,beta,gamma)
    if (.not.lrelax) then
      do i = 1,3
        rvcfg(1,i,ncf) = rv(1,i)
        rvcfg(2,i,ncf) = rv(2,i)
        rvcfg(3,i,ncf) = rv(3,i)
      enddo
    endif
  endif
!
  if (.not.lrelax) then
    do i = 1,nasym
      xafrac(i) = dmod(xafrac(i)+1.0_dp,1.0_dp)
      yafrac(i) = dmod(xafrac(i)+1.0_dp,1.0_dp)
      zafrac(i) = dmod(xafrac(i)+1.0_dp,1.0_dp)
    enddo
    if (lsymopt) then
      call equpos(.true.,.false.)
    else
      do i = 1,numat
        xfrac(i) = xafrac(i)
        yfrac(i) = yafrac(i)
        zfrac(i) = zafrac(i)
      enddo
    endif
  endif
#ifdef TRACE
  call trace_out('gabcfg')
#endif
!
  return
  end
