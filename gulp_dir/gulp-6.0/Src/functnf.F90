  subroutine functnf(iflag,n,xc,fc,gc,lnonfitcall)
!
!  Supplies the function and derivatives where the second
!  derivatives are determined by finite differences. This
!  routine differs from functn in that it determines the 
!  full second derivative matrices rather than just the 
!  Hessian matrix elements for the optimisation variables.
!
!   5/09 Created from funct
!   6/10 Structure was not being reset after final finite difference
!        point and before calculation of gradients for return. Call to
!        x0 routines added to correct this.
!   9/13 Outderv now called by all processors
!   7/17 Modified to allow for parallel second derivative calls
!   9/17 Modified to handle case where loptcellpar has been used during
!        optimisation and then strain second derivatives are needed
!   1/18 Trace added
!   6/18 Clean out of unused variables
!   8/18 Modified due to changes in lstraincell algorithm
!   8/18 Adding 1 to strains 1-3 removed
!   9/18 Handling of strain reseting corrected
!   7/19 Corrected for second strain derivatives when not lstraincellprop
!   8/19 Correction to initialisation of x0str
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
!  Julian Gale, CIC, Curtin University, August 2019
!
  use configurations, only : lbsmat, rvcfg
  use control
  use current,        only : numat, nstrains, xclat, yclat, zclat, ndim, radf, nbsmat, strain, rv, ncf
  use derivatives,    only : xdrv, ydrv, zdrv, strderv, derv2, derv3, sderv2, raderv
  use gulp_files
  use general
  use parallel
  use symmetry,       only : lstr, lsymderv, lsymderv2
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4),            intent(inout) :: iflag            ! Flag to indicate level of derivatives
  integer(i4),            intent(in)    :: n                ! Number of variables
  logical,                intent(in)    :: lnonfitcall      ! Flag to indicate whether this is a call from fitfun or not
  real(dp),               intent(out)   :: fc               ! Energy (on return)
  real(dp),               intent(in)    :: xc(*)            ! Structure in linear array
  real(dp),               intent(out)   :: gc(*)            ! Gradients in linear array (on return)
!
!  Local variables
!
  integer(i4)                           :: i                ! Looping variable
  integer(i4)                           :: iloc             ! Index of atom on local node
  integer(i4)                           :: indi             ! Index for Cartesian coordinate of i in derv2
  integer(i4)                           :: indif            ! Index for Cartesian coordinate of i in derv2
  integer(i4)                           :: ind3             ! 3*numat - used in breathing shell indexing in derv2
  integer(i4)                           :: indj             ! Index for Cartesian coordinate of i in derv2
  integer(i4)                           :: ix               ! Looping variable over Cartesian components for i
  integer(i4)                           :: j                ! Looping variable
  integer(i4)                           :: jx               ! Looping variable over Cartesian components for j
  integer(i4)                           :: status           ! Error flag for memory allocation
  logical                               :: lbreathing       ! Logical to indicate whether there are any breathing shells
  logical                               :: lgeometryOK
  logical                               :: lgrad1           ! Flag to indicate whether first derivatives are to be calculated
  logical                               :: lgrad2           ! Flag to indicate whether second derivatives are to be calculated
  logical                               :: lsymdervsave     ! Saves value of lsymderv
  logical                               :: lsymderv2save    ! Saves value of lsymderv2
  real(dp)                              :: g_cpu_time
  real(dp)                              :: dfindiffcorr     ! Term used in diagonal strain correction (1 - delta**2)
  real(dp)                              :: rfindiffc        ! 1/(2 x finite difference value for Cartesians)
  real(dp)                              :: rfindiffs        ! 1/(2 x finite difference value for strains)
  real(dp),                        save :: tdmax = 0.0_dp
  real(dp)                              :: t1
  real(dp)                              :: t2
  real(dp)                              :: xci              ! Stores the original value of the variable being finite differenced
  real(dp), dimension(:),   allocatable :: gbb              ! Workspace for backward radial derivatives
  real(dp), dimension(:),   allocatable :: gbf              ! Workspace for forward radial derivatives
  real(dp), dimension(:,:), allocatable :: gcb              ! Workspace for backward Cartesian derivatives
  real(dp), dimension(:,:), allocatable :: gcf              ! Workspace for forward Cartesian derivatives
  real(dp)                              :: gsb(6)           ! Workspace for backward strain derivatives
  real(dp)                              :: gsf(6)           ! Workspace for forward strain derivatives
  real(dp)                              :: x0save(6)        ! Strain workspace
  real(dp)                              :: x0str(6)         ! Strain workspace
  real(dp)                              :: rvsave(3,3)      ! Save initial cell
#ifdef TRACE
  call trace_in('functnf')
#endif
!
  t1 = g_cpu_time()
  lgrad1 = (iflag.ge.1)
  lgrad2 = (iflag.ge.2)
  lbreathing = (nbsmat.gt.0) 
!
!  Store values of lsymderv(2) since these must be turned off during finite differencing as it breaks symmetry
!
  lsymdervsave = lsymderv
  lsymderv2save = lsymderv2
  lsymderv = .false.
  lsymderv2 = .false.
!************************************************************
!  Convert optimisation variables to linear structure array *
!************************************************************
  if (lnonfitcall) call xctox0(n,xc,lgeometryOK)
!************************************************************
!  Convert linear structure array to main structure arrays  *
!************************************************************
  if (lx0centroid) then
    call x0tostrcentroid
  else
    call x0tostr
  endif
!
!  Save initial cell
!
  if (ndim.gt.0) then
    if (lstraincell) then
      rvsave(1:3,1:3) = rvcfg(1:3,1:3,ncf)
    else
      rvsave(1:3,1:3) = rv(1:3,1:3)
    endif
    if (lstraincell.and..not.lstraincellprop) then
      x0str(1:nstrains) = strain(1:nstrains)
      x0save(1:nstrains) = x0str(1:nstrains)
!
!  Temporarily set x0str to zero for working here
!
      x0str(1:nstrains) = 0.0_dp
!
!  If properties are to be determined with an applied strain, but not at finite strain
!  then applied the initial strains before finite differences
!
      call x0strainf(x0str,lgeometryOK)
      rvsave(1:3,1:3) = rv(1:3,1:3)
    endif
  endif
!
  lfirst = .true.
  if (lgrad2) then
!
!  Allocate local arrays
!
    allocate(gcb(3,numat),stat=status)
    if (status/=0) call outofmemory('functnf','gcb')
    allocate(gcf(3,numat),stat=status)
    if (status/=0) call outofmemory('functnf','gcf')
    if (lbreathing) then
      allocate(gbb(numat),stat=status)
      if (status/=0) call outofmemory('functnf','gbb')
      allocate(gbf(numat),stat=status)
      if (status/=0) call outofmemory('functnf','gbf')
    endif
!
!  Initialise second derivatives
!
    call initdervs(lgrad1,lgrad2)
!*********************************************************************************************
!  Compute second derivatives by finite differences with respect to all Cartesian components *
!*********************************************************************************************
    rfindiffc = 0.5_dp/findiffc
    rfindiffs = 0.5_dp/findiffs
    dfindiffcorr = 0.25_dp*findiffs**2/(1.0_dp - 0.25_dp*findiffs**2)
    ind3 = 3*numat
    do i = 1,numat
      do ix = 1,3
!
!  Save coordinate
!
        if (ix.eq.1) then
          xci = xclat(i)
        elseif (ix.eq.2) then
          xci = yclat(i)
        else
          xci = zclat(i)
        endif
!
!  Forward step
!
        if (ix.eq.1) then
          xclat(i) = xci + findiffc
        elseif (ix.eq.2) then
          yclat(i) = xci + findiffc
        else
          zclat(i) = xci + findiffc
        endif
        call energy(fc,lgrad1,.false.)
        do j = 1,numat
          gcf(1,j) = xdrv(j)
          gcf(2,j) = ydrv(j)
          gcf(3,j) = zdrv(j)
        enddo
        if (lstr) then
          call strfin(.false.)
          gsf(1:nstrains) = strderv(1:nstrains)
        endif
        if (lbreathing) then
          gbf(1:numat) = raderv(1:numat)
        endif
!
!  Backward step
!
        if (ix.eq.1) then
          xclat(i) = xci - findiffc
        elseif (ix.eq.2) then
          yclat(i) = xci - findiffc
        else
          zclat(i) = xci - findiffc
        endif
        call energy(fc,lgrad1,.false.)
        do j = 1,numat
          gcb(1,j) = xdrv(j)
          gcb(2,j) = ydrv(j)
          gcb(3,j) = zdrv(j)
        enddo
        if (lstr) then
          call strfin(.false.)
          gsb(1:nstrains) = strderv(1:nstrains)
        endif
        if (lbreathing) then
          gbb(1:numat) = raderv(1:numat)
        endif
!
!  Find node that has atom i locally
!
        iloc = atom2local(i)
        if (iloc.gt.0) then
          indi = 3*(iloc-1)
!
!  Calculate Cartesian-Cartesian second derivatives
!
          do j = 1,numat
            indj = 3*(j-1)
            do jx = 1,3
              derv2(indj+jx,indi+ix) = rfindiffc*(gcf(jx,j) - gcb(jx,j))
            enddo
          enddo
          if (lstr) then
!
!  Calculate mixed Cartesian-strain second derivatives
!
            do jx = 1,nstrains
              derv3(indi+ix,jx) = rfindiffc*(gsf(jx) - gsb(jx))
            enddo
          endif
          if (lbreathing) then
!
!  Calculate mixed Cartesian-radial second derivatives
!
            do j = 1,numat
              derv2(ind3+j,indi+ix) = rfindiffc*(gbf(j) - gbb(j))
            enddo
          endif
        endif
!
!  Restore coordinate
!
        if (ix.eq.1) then
          xclat(i) = xci
        elseif (ix.eq.2) then
          yclat(i) = xci
        else
          zclat(i) = xci
        endif
      enddo
    enddo
    if (lstr) then
!******************************************************************************************
!  Compute second derivatives by finite differences with respect to all strain components *
!******************************************************************************************
!
!  Set initial strain values
!
      if (lstraincell.and.lstraincellprop) then
        x0str(1:nstrains) = strain(1:nstrains)
      else
        x0str(1:nstrains) = 0.0_dp
      endif
!
      do ix = 1,nstrains
!
!  Set cell
!
        rv(1:3,1:3) = rvsave(1:3,1:3)
!
!  Save coordinate
!
        xci = x0str(ix)
!
!  Forward step
!
        x0str(ix) = xci + findiffs
        call x0strainf(x0str,lgeometryOK)
        if (lx0centroid) then
          call x0tostrcentroid
        else
          call x0tostr
        endif
        call energy(fc,lgrad1,.false.)
        call strfin(.false.)
        gsf(1:nstrains) = strderv(1:nstrains)
        if (lbreathing) then
          gbf(1:numat) = raderv(1:numat)
        endif
!
!  Reset cell
!
        rv(1:3,1:3) = rvsave(1:3,1:3)
!
!  Backward step
!
        x0str(ix) = xci - findiffs
        call x0strainf(x0str,lgeometryOK)
        if (lx0centroid) then
          call x0tostrcentroid
        else
          call x0tostr
        endif
        call energy(fc,lgrad1,.false.)
        call strfin(.false.)
        gsb(1:nstrains) = strderv(1:nstrains)
        if (lbreathing) then
          gbb(1:numat) = raderv(1:numat)
        endif
!
!  Calculate strain-strain second derivatives
!
        do jx = 1,nstrains
          sderv2(jx,ix) = rfindiffs*(gsf(jx) - gsb(jx))
        enddo
        if (lbreathing) then
!
!  Calculate mixed strain-radial second derivatives
!
          do j = 1,numat
            derv3(ind3+j,ix) = rfindiffs*(gbf(j) - gbb(j))
          enddo
        endif
!
!  Reset cell
!
        rv(1:3,1:3) = rvsave(1:3,1:3)
!
!  Restore coordinate
!
        x0str(ix) = xci
        call x0strainf(x0str,lgeometryOK)
        if (lx0centroid) then
          call x0tostrcentroid
        else
          call x0tostr
        endif
      enddo
    endif
    if (lbreathing) then
!**********************************************************************************************
!  Compute second derivatives by finite differences with respect to all breathing shell radii *
!**********************************************************************************************
      do i = 1,numat
        if (lbsmat(i)) then
!
!  Save radius
!
          xci = radf(i)
!
!  Forward step
!
          radf(i) = xci + findiffc
          call energy(fc,lgrad1,.false.)
          do j = 1,numat
            gcf(1,j) = xdrv(j)
            gcf(2,j) = ydrv(j)
            gcf(3,j) = zdrv(j)
          enddo
          gbf(1:numat) = raderv(1:numat)
!
!  Backward step
!
          radf(i) = xci - findiffc
          call energy(fc,lgrad1,.false.)
          do j = 1,numat
            gcb(1,j) = xdrv(j)
            gcb(2,j) = ydrv(j)
            gcb(3,j) = zdrv(j)
          enddo
          gbb(1:numat) = raderv(1:numat)
!
!  Calculate mixed Cartesian-radial second derivatives
!
          do j = 1,numat
            indj = 3*(j-1)
            derv2(indj+1,ind3+i) = rfindiffc*(gcf(1,j) - gcb(1,j))
            derv2(indj+2,ind3+i) = rfindiffc*(gcf(2,j) - gcb(2,j))
            derv2(indj+3,ind3+i) = rfindiffc*(gcf(3,j) - gcb(3,j))
          enddo
!
!  Calculate radial-radial second derivatives
!
          do j = 1,numat
            derv2(ind3+j,ind3+i) = rfindiffc*(gbf(j) - gbb(j))
          enddo
!
!  Restore coordinate
!
          radf(i) = xci
        endif
      enddo
    endif
!****************************************************************************
!  Average duplicate elements in matrices to ensure that they are symmetric *
!****************************************************************************
    if (nprocs.eq.1) then
      do i = 2,numat
        indi = 3*(i-1)
        do j = 1,i-1
          indj = 3*(j-1)
          do ix = 1,3
            do jx = 1,3
              derv2(indj+jx,indi+ix) = 0.5_dp*(derv2(indj+jx,indi+ix) + derv2(indi+ix,indj+jx))
              derv2(indi+ix,indj+jx) = derv2(indj+jx,indi+ix)
            enddo
          enddo
        enddo
      enddo
    endif
    if (lstr) then
      do ix = 2,nstrains
        do jx = 1,ix-1
          sderv2(jx,ix) = 0.5_dp*(sderv2(jx,ix) + sderv2(ix,jx))
          sderv2(ix,jx) = sderv2(jx,ix)
        enddo
      enddo
    endif
    if (lbreathing) then
      do i = 1,numat
        if (lbsmat(i)) then
          do jx = 1,3*numat
            derv2(jx,ind3+i) = 0.5_dp*(derv2(jx,ind3+i) + derv2(ind3+i,jx))
            derv2(ind3+i,jx) = derv2(jx,ind3+i)
          enddo
        endif
      enddo
      do i = 2,numat
        if (lbsmat(i)) then
          do j = 1,i-1
            if (lbsmat(j)) then
              derv2(ind3+j,ind3+i) = 0.5_dp*(derv2(ind3+j,ind3+i) + derv2(ind3+i,ind3+j))
              derv2(ind3+i,ind3+j) = derv2(ind3+j,ind3+i)
            endif
          enddo
        endif
      enddo
    endif
!***********************************************************************************
!  Ensure that translational invariance is obeyed by setting the diagonal elements *
!***********************************************************************************
    do i = 1,numat
!
!  Find node that has atom i locally
!
      iloc = atom2local(i)
      if (iloc.gt.0) then
        indi = 3*(iloc-1)
        indif = 3*(i-1)
        do ix = 1,3
          do jx = 1,3
            derv2(indif+jx,indi+ix) = 0.0_dp
          enddo
        enddo
        do j = 1,i-1
          indj = 3*(j-1)
          do ix = 1,3
            do jx = 1,3
              derv2(indif+jx,indi+ix) = derv2(indif+jx,indi+ix) - derv2(indj+jx,indi+ix)
            enddo
          enddo
        enddo
        do j = i+1,numat
          indj = 3*(j-1)
          do ix = 1,3
            do jx = 1,3
              derv2(indif+jx,indi+ix) = derv2(indif+jx,indi+ix) - derv2(indj+jx,indi+ix)
            enddo
          enddo
        enddo
      endif
    enddo
!
!  Free local arrays
!
    if (lbreathing) then
      deallocate(gbf,stat=status)
      if (status/=0) call deallocate_error('functnf','gbf')
      deallocate(gbb,stat=status)
      if (status/=0) call deallocate_error('functnf','gbb')
    endif
    deallocate(gcf,stat=status)
    if (status/=0) call deallocate_error('functnf','gcf')
    deallocate(gcb,stat=status)
    if (status/=0) call deallocate_error('functnf','gcb')
  endif
!***********************************************************************************************
!  Compute the first derivatives and associated quantities according to normal funct procedure *
!***********************************************************************************************
!
!  Restore original strains
!
  if (ndim.gt.0) then
    if (lstraincell.and..not.lstraincellprop) then
      x0str(1:nstrains) = x0save(1:nstrains)
    endif
  endif
!
!  Reset structure
!
  if (lx0centroid) then
    call x0tostrcentroid
  else
    call x0tostr
  endif
!
!  Evaluate function and first derivatives
!
  call energy(fc,lgrad1,.false.)
!
!  For surface, get surface energy
!
  if (lseok) call surfaceenergy(fc)
!
!  Complete strain derivatives
!
  if (lstr) call strfin(.false.)
!
!  Output second derivatives 
!
  if (lgrad2) call outderv
!
!  Output energy and derivatives to a file
!
  if (ldrv.and.ioproc) call outdrv(fc,lgrad1,lgrad2)
!
!  Option to write out a .frc file for QMPOT
!
  if (lfrc.and.ioproc) call outfrc(fc,lgrad1,.false.)
!
!  First derivative handling
!
  if (lgrad1.and.lnonfitcall) then
    call getderv1(n,xc,gc,lgrad2,.false.)
  endif
!
!  Restore values of lsymderv(2) before call to obtain first derivatives
!
  lsymderv = lsymdervsave
  lsymderv2 = lsymderv2save
!
  t2 = g_cpu_time()
  if (t2-t1.gt.tdmax) tdmax = t2 - t1
  if (timmax.gt.0.0_dp) then
    if ((timmax-t2).lt.tdmax.and..not.lrelax) iflag = -1
  endif
  if (.not.lgeometryOK) iflag = -2
#ifdef TRACE
  call trace_out('functnf')
#endif
!
  return
  end
