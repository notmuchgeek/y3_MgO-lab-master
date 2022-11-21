  subroutine setqonsas2D(lgrad2)
!
!  Subroutine calculates the charges on the solvent accessible surface
!  2D matrix for cosmoA version
!
!  Note that the iterative algorithm for finding the charges on the 
!  solvent accessible surface can only be used for COSMO and not for
!  COSMIC since while the charges could be computed the derivatives
!  rely on the presence of the inverse of A being present in cosmoA
!
!   4/17 Created from setqonsas
!   4/17 qonsas option added
!   4/17 Parallelisation added
!   6/17 Module files renamed to gulp_files
!   6/18 Extra flag passed to dbcgsolve
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
!  Julian Gale, CIC, Curtin University, June 2018
!
  use configurations, only : totalchargecfg
  use control,        only : keyword, literativeQ
  use cosmic
  use current
  use element,        only : nqitermax, qitertol
  use gulp_files,     only : lsas, sasfile
  use iochannels
  use maths,          only : qsolver
  use parallel
  implicit none
!
!  Passed variables
!
  logical,      intent(in)                    :: lgrad2
!
!  Local variables
!
  integer(i4)                                 :: i
  integer(i4)                                 :: info
  integer(i4)                                 :: ipts
  integer(i4)                                 :: iter
  integer(i4)                                 :: j
  integer(i4)                                 :: lwrk
  integer(i4)                                 :: nblock
  integer(i4),                           save :: nptslast = 0
  integer(i4)                                 :: status
  integer(i4),              allocatable, save :: nlocptr(:)
  real(dp)                                    :: err
  real(dp)                                    :: sumfinalq
  real(dp)                                    :: dqstot(3)
  real(dp),                 allocatable, save :: wrk(:)
!
!  Functions
!
  integer(i4)                                 :: ilaenv
!
#ifdef MPI
  integer                                     :: idesA(9)
  integer                                     :: idesB(9)
  integer                                     :: ifails
  integer,     dimension(:), allocatable      :: ipivot
  integer                                     :: ldm
  integer                                     :: nb
  integer                                     :: np
  integer                                     :: nv
#endif
!
  if (.not.lgrad2.and..not.lcosmic.and.literativeQ) then
!************************
!  Iterative algorithm  *
!************************
!
!  On first call initialise charges on sas to zero
!
    if (npts.ne.nptslast) then
      qonsas(1:npts) = 0.0_dp
      nptslast = npts
    endif
    if (qsolver.eq.3.or.(nprocs.gt.1.and.qsolver.eq.1)) then
!
!  Solve using biconjugate gradients - now default in parallel for iterative approach
!  since itpack is not parallelised. 
!
      allocate(nlocptr(npts))
      do i = 1,npts
        nlocptr(i) = i
      enddo
      call dbcgsolve(1_i4,npts,nptsonnode,node2pts,cosmoA,maxcosmoA,cosmoBq,qonsas,qitertol,nqitermax,iter,err)
      deallocate(nlocptr)
      do i = 1,npts
        qonsas(i) = - qonsas(i)
      enddo
    elseif (qsolver.eq.2) then
#ifdef MPI
      if (nprocs.gt.1) then
!**************************
!  Solve using scalapack  *
!**************************
!
!  Set local block size
!
        nb = nblocksizesas
        np = nprocs
        nv = npts
        ifails = 0
!
!  Set up Blacs descriptors
!
        ldm = maxcosmoA
        call descinit( idesA, nv, nv, nb, nb, 0, 0, iBlacsContext, ldm, ifails )
!
        if (ifails.ne.0) then
          call outerror('initialisation in descinit failed',0_i4)
          call stopnow('setqonsas2D')
        endif
!
        ldm = maxnpts
        call descinit( idesB, nv, 1, nb, nb, 0, 0, iBlacsContext, ldm, ifails )
!
        if (ifails.ne.0) then
          call outerror('initialisation in descinit failed',0_i4)
          call stopnow('setqonsas2D')
        endif
!
        allocate(ipivot(4*npts),stat=status)
        if (status/=0) call outofmemory('setqonsas2D','ipivot')
!
        qonsas(1:npts) = - cosmoBq(1:npts)
        call pdgesv(npts,1,cosmoA,1,1,idesA,ipivot,qonsas,1,1,idesB,info)
!
        if (info.ne.0) then
          call outerror('pdgesv failed in setqonsas2D',0_i4)
          call stopnow('setqonsas2D')
        endif
!
!  Broadcase the result from node 0
!
        call sendall(qonsas,npts,0_i4,"setqonsas2D","qonsas")
!
        deallocate(ipivot,stat=status)
        if (status/=0) call deallocate_error('setqonsas2D','ipivot')
      else
#endif
!***********************
!  Solve using lapack  *
!***********************
        qonsas(1:npts) = - cosmoBq(1:npts)
!
!  Compute optimal workspace
!
        nblock = ilaenv(1,'DSYTRF','U',npts,-1,-1,-1)
        lwrk = nblock*npts
        allocate(wrk(lwrk),stat=status)
        if (status/=0) call outofmemory('setqonsas2D','wrk')
        allocate(nlocptr(npts),stat=status)
        if (status/=0) call outofmemory('setqonsas2D','nlocptr')
!
        call dsysv('U',npts,1_i4,cosmoA,maxcosmoA,nlocptr,qonsas,npts,wrk,lwrk,info)
!
        deallocate(wrk,stat=status)
        if (status/=0) call deallocate_error('setqonsas2D','wrk')
        deallocate(nlocptr,stat=status)
        if (status/=0) call deallocate_error('setqonsas2D','nlocptr')
#ifdef MPI
      endif
#endif
    else
!
!  Solve using itpack
!
      call iterinv(npts,cosmoA,maxcosmoA,cosmoBq,qonsas,info)
    endif
  else
!*********************
!  Direct algorithm  *
!*********************
#ifdef MPI
    if (nprocs.gt.1) then
!
!  Invert A matrix
!
      call setcosmoainv2D
!
!  Set local block size
!
      nb = nblocksizesas
      np = nprocs
      nv = npts
      ifails = 0
!
!  Set up Blacs descriptors
!
      ldm = maxcosmoA
      call descinit( idesA, nv, nv, nb, nb, 0, 0, iBlacsContext, ldm, ifails )
!
      if (ifails.ne.0) then
        call outerror('initialisation in descinit failed',0_i4)
        call stopnow('setqonsas2D')
      endif
!
      ldm = maxnpts
      call descinit( idesB, nv, 1, nb, nb, 0, 0, iBlacsContext, ldm, ifails )
!
      if (ifails.ne.0) then
        call outerror('initialisation in descinit failed',0_i4)
        call stopnow('setqonsas2D')
      endif
!
!  Calculate charges on SAS
!
      call pdgemv('N',npts,npts,-1.0d0,cosmoA,1,1,idesA,cosmoBq,1,1,idesB,1,0.0d0,qonsas,1,1,idesB,1)
!
!  Broadcase the result from node 0
!
      call sendall(qonsas,npts,0_i4,"setqonsas2D","qonsas")
    else
#endif
!
!  Invert A matrix
!
      call setcosmoainv2D
!
!  Calculate charges on SAS
!
      call dsymv('U',npts,-1.0d0,cosmoA,maxcosmoA,cosmoBq,1,0.0d0,qonsas,1)
#ifdef MPI
    endif
#endif
  endif
!
!  Calculate sum of SAS charges  
!
  qonsastot = 0.0_dp
  do ipts = 1,npts
    qonsastot = qonsastot + qonsas(ipts)
  enddo
!**********************
!  COSMIC correction  *
!**********************
  if (lcosmic.and.totsegweight.gt.1.0d-12) then
    deltaq = (qonsastot+totalchargecfg(ncf)-qonsastarget(ncf))/totsegweight
  else 
    deltaq = 0.0_dp
  endif
!*********************************
!  Print out charges on surface  *
!*********************************
  if (index(keyword,'qsas').ne.0.and.ioproc) then   
!
!  Calculate dipole for output
!
    dqstot(1:3) = 0.0_dp
    do ipts = 1,npts
      dqstot(1) = dqstot(1) + (qonsas(ipts) - deltaq*segweight(ipts))*spxyz(1,ipts)
      dqstot(2) = dqstot(2) + (qonsas(ipts) - deltaq*segweight(ipts))*spxyz(2,ipts)
      dqstot(3) = dqstot(3) + (qonsas(ipts) - deltaq*segweight(ipts))*spxyz(3,ipts)
    enddo
!
    write(ioout,'(/,''  Total charge on Solvent Accessible Surface = '',f12.6)') qonsastot
    if (ndim.eq.0) then
      write(ioout,'(''  Dipole in  X on Solvent Accessible Surface = '',f12.6)') dqstot(1)
      write(ioout,'(''  Dipole in  Y on Solvent Accessible Surface = '',f12.6)') dqstot(2)
      write(ioout,'(''  Dipole in  Z on Solvent Accessible Surface = '',f12.6)') dqstot(3)
    elseif (ndim.eq.1) then
      write(ioout,'(''  Dipole in  Y on Solvent Accessible Surface = '',f12.6)') dqstot(2)
      write(ioout,'(''  Dipole in  Z on Solvent Accessible Surface = '',f12.6)') dqstot(3)
    elseif (ndim.eq.2) then
      write(ioout,'(''  Dipole in  Z on Solvent Accessible Surface = '',f12.6)') dqstot(3)
    endif
    write(ioout,'(/,''  Charges on Solvent Accessible Surface :'',/)')
    write(ioout,'(''-------------------------------------------------------------------------------'')')
    write(ioout,'(''  Segment    x (Angs)     y (Angs)     z (Angs)       Charge      Weight'')')
    write(ioout,'(''-------------------------------------------------------------------------------'')')
    sumfinalq = 0.0_dp
    do ipts = 1,npts
      write(ioout,'(i8,2x,3(1x,f12.5),2(2x,f12.8))') ipts,spxyz(1,ipts),spxyz(2,ipts),spxyz(3,ipts),  &
        (qonsas(ipts) - deltaq*segweight(ipts)),segweight(ipts)
      sumfinalq = sumfinalq + qonsas(ipts) - deltaq*segweight(ipts)
    enddo
    write(ioout,'(''-------------------------------------------------------------------------------'')')
    write(ioout,'(''    Total '',41x,f12.8)') sumfinalq
    write(ioout,'(''-------------------------------------------------------------------------------'',/)')
  endif
!********************
!  Output SAS file  *
!********************
  if (lsas.and.ioproc) then
!
!  Open file 
!
    if (sasfile(1:1).ne.' ') then
      open(22,file=sasfile,status='unknown')
    else
      open(22,status='unknown')
    endif
!     
!  Output number of points
!
    write(22,'(i8)') npts
!     
!  Output dimensionality
!
    write(22,'(i8)') ndim
!     
!  Output vectors
!
    if (ndim.gt.0) then
      do i = 1,ndim
        write(22,'(3(f14.6,1x))') (rv(j,i),j=1,ndim)
      enddo
    endif
!     
!  Output point information
!     
    do ipts = 1,npts
      write(22,'(i8,1x,i6,1x,3(1x,f12.6),1x,f14.8)') &
        ipts,cosmoatomptr(ipts),spxyz(1,ipts),spxyz(2,ipts),spxyz(3,ipts), &
        (qonsas(ipts) - deltaq*segweight(ipts))
    enddo
!
!  Close file
!
    close(22)
  endif
!
  return
  end
