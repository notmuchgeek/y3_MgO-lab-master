  subroutine shmin(fc,spring,bspring,gtol,loutput,niter,ifail)
!
!  Conjugate gradients optimisation routine for shell positions
!
!  10/13 Created from qtmin
!  12/14 nbs, nbss and nbsptr moved to shells modul
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
!  Julian Gale, CIC, Curtin University, Feb 2018
!
  use control
  use current,     only : xalat, yalat, zalat, radf, nbsmat
  use derivatives, only : xdrv, ydrv, zdrv, raderv
  use iochannels
  use m_conjgr,    only : conjgr
  use moldyn,      only : lfix
  use parallel
  use shells,      only : nshell, nshptr, moptit, nbsptr
#ifdef TRACE
  use trace,       only : trace_in, trace_out
#endif

  implicit none
!
!  Passed variables
!
  integer(i4),          intent(out)         :: ifail
  integer(i4),          intent(out)         :: niter
  logical,              intent(in)          :: loutput
  real(dp),             intent(out)         :: fc
  real(dp),             intent(in)          :: gtol
  real(dp),             intent(in)          :: spring(*)
  real(dp),             intent(in)          :: bspring(*)
!
!  Local variables
!
  integer(i4)                               :: i
  integer(i4)                               :: ibp
  integer(i4)                               :: isp
  integer(i4)                               :: n
  integer(i4)                               :: nfree
  integer(i4)                               :: status
  logical                                   :: okf
  real(dp)                                  :: dxmax
  real(dp)                                  :: gnorm
  real(dp)                                  :: cgcntr(0:20)
  real(dp), dimension(:), allocatable, save :: aux
  real(dp), dimension(:), allocatable, save :: gc
  real(dp), dimension(:), allocatable, save :: xc
#ifdef TRACE
  call trace_in('shmin')
#endif
!
!  Initialisation
!
  okf = .false.
  cgcntr(0:20) = 0.0_dp
  dxmax = 0.1_dp
!
!  Allocate workspace
!
  allocate(aux(8*nshell),stat=status)
  if (status/=0) call outofmemory('shmin','aux')
  allocate(gc(4*nshell),stat=status)
  if (status/=0) call outofmemory('shmin','gc')
  allocate(xc(4*nshell),stat=status)
  if (status/=0) call outofmemory('shmin','xc')
!
!  Place coordinates in xc
!
  n = 0
  do i = 1,nshell
    isp = nshptr(i)
    if (.not.lfix(isp)) then
      xc(n+1) = xalat(isp)
      xc(n+2) = yalat(isp)
      xc(n+3) = zalat(isp)
      n = n + 3
    endif
  enddo
  if (nbsmat.ne.0) then
    do i = 1,nbsmat
      ibp = nbsptr(i)
      if (.not.lfix(ibp)) then
        n = n + 1
        xc(n) = radf(ibp)
      endif
    enddo
  endif
!
!  Output heading
!
  if (loutput) then
    write(ioout,'(/,''  SHmin: minimisation of shell positions: '',/)')
  endif
!
!  Start loop over calls to conjugate gradients until converged
!
  niter = 0
  do while (.not.okf.and.niter.lt.moptit)
    niter = niter + 1
!
!  Move positions from xc to coordinate arrays
!
    nfree = 0
    do i = 1,nshell
      isp = nshptr(i)
      if (.not.lfix(isp)) then
        xalat(isp) = xc(nfree+1)
        yalat(isp) = xc(nfree+2)
        zalat(isp) = xc(nfree+3)
        nfree = nfree + 3
      endif
    enddo
    if (nbsmat.ne.0) then
      do i = 1,nbsmat
        ibp = nbsptr(i)
        if (.not.lfix(ibp)) then
          nfree = nfree + 1
          radf(ibp) = xc(nfree)
        endif
      enddo
    endif
!
!  Call function
!
    call mdfunct(0_i4,fc,.false.,niter.eq.1,.true.)
!
!  Move gradients into gc force array
!
    gnorm = 0.0_dp
    nfree = 0
    do i = 1,nshell
      isp = nshptr(i)
      if (.not.lfix(isp)) then
        gc(nfree+1) = - xdrv(isp)*spring(i)
        gc(nfree+2) = - ydrv(isp)*spring(i)
        gc(nfree+3) = - zdrv(isp)*spring(i)
        gnorm = gnorm + xdrv(isp)*xdrv(isp) + ydrv(isp)*ydrv(isp) + zdrv(isp)*zdrv(isp)
        nfree = nfree + 3
      endif
    enddo
    if (nbsmat.ne.0) then
      do i = 1,nbsmat
        ibp = nbsptr(i)
        if (.not.lfix(ibp)) then
          nfree = nfree + 1
          gc(nfree) = - raderv(ibp)*bspring(i)
          gnorm = gnorm + raderv(ibp)*raderv(ibp)
        endif
      enddo
    endif
    gnorm = gnorm/dble(n)
    gnorm = sqrt(gnorm)
!
!  Print latest info
!
    if (loutput) then
      write(ioout,'(''  SHmin: '',i5,'' Energy = '',f18.8,'' Gnorm = '',g24.12)') niter,fc,gnorm
    endif
!
!  Call minimiser
!
    call conjgr( n, xc, gc, dxmax, gtol, cgcntr, aux )
!
!  Test for convergence
!
    okf = (int(cgcntr(0)) .eq. 0)  
  enddo
!
!  Finalisation
!
  if (okf) then
    ifail = 0
  elseif (niter.ge.moptit) then
    ifail = 1
  endif
!
!  Deallocate workspace
!
  deallocate(xc,stat=status)
  if (status/=0) call deallocate_error('shmin','xc')
  deallocate(gc,stat=status)
  if (status/=0) call deallocate_error('shmin','gc')
  deallocate(aux,stat=status)
  if (status/=0) call deallocate_error('shmin','aux')
#ifdef TRACE
  call trace_out('shmin')
#endif
  return
  end
