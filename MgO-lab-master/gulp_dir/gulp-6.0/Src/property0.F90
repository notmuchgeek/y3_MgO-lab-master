  subroutine property0(lprint)
!
!  Calculate the properties of clusters
!
!  10/05 Created
!   5/06 Mass now set using species values
!   8/06 Missing occupancy factor in inertia calculation added
!  12/07 Unused variables removed
!   1/09 Integer datatypes all explicitly declared
!   2/09 Modified to accommodate new version of FoX and gulp_cml
!  12/10 Wrong spelling of principal corrected
!   7/13 Calculation of translation partition function and free 
!        energy added
!   8/13 Trap for numat = 1 added for rotation
!   9/13 Calculation of Raman susceptibilities added
!   3/14 n3f properly defined
!  12/14 nbs, nbss and nbsptr moved to shells modul
!   1/16 Separate rotational and translational free energies printed
!   9/16 Modified to allow for parallel second derivatives
!   9/16 Matrix inversion moved into subroutine
!   2/17 Blocksize added to call to matrix_inversion_library
!   8/17 Parallel handling of qD added
!   1/18 Trace added
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
!  Julian Gale, CIC, Curtin University, January 2018
!
  use g_constants,      only : avogadro, evtoj, boltz, kjmtoev
  use control,          only : lraman
  use current
  use derivatives,      only : derv2, dervi, maxd2
  use element
  use general,          only : nwarn
  use gulp_cml,         only : lcml
  use gulp_cml_props,   only : gulp_cml_output_inertia
  use iochannels
  use parallel
  use partial
  use shells,           only : ncore, nshell, nshptr, nbs, nbss, nbsptr
  use species,          only : massspec
  use symmetry,         only : symnocfg
  use times
#ifdef TRACE
  use trace,            only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  logical,                            intent(in) :: lprint
!
!  Local variables
!
  integer(i4)                                    :: i
  integer(i4)                                    :: ifail
  integer(i4)                                    :: indk
  integer(i4)                                    :: indl
  integer(i4)                                    :: iptr
  integer(i4)                                    :: ix
  integer(i4)                                    :: iy
  integer(i4)                                    :: iz
  integer(i4)                                    :: ixs
  integer(i4)                                    :: iys
  integer(i4)                                    :: izs
  integer(i4)                                    :: j
  integer(i4)                                    :: jptr
  integer(i4)                                    :: jx
  integer(i4)                                    :: jy
  integer(i4)                                    :: jz
  integer(i4)                                    :: jxs
  integer(i4)                                    :: jys
  integer(i4)                                    :: jzs
  integer(i4)                                    :: k
  integer(i4)                                    :: msvar
  integer(i4)                                    :: n
  integer(i4)                                    :: n3f
  integer(i4)                                    :: nbsc
  integer(i4)                                    :: ni
  integer(i4)                                    :: nj
  integer(i4)                                    :: status
  real(dp)                                       :: g_cpu_time
  real(dp)                                       :: dx
  real(dp)                                       :: dy
  real(dp)                                       :: dz
  real(dp)                                       :: Gn
  real(dp)                                       :: Gr
  real(dp)                                       :: Gt
  real(dp)                                       :: inertia(3,3)
  real(dp)                                       :: inertia_product
  real(dp)                                       :: lnZr
  real(dp)                                       :: lnZt
  real(dp)                                       :: mi
  real(dp)                                       :: press_local
  real(dp)                                       :: qak
  real(dp),    dimension(:,:), allocatable       :: qD
  real(dp),    dimension(:),   allocatable       :: qshell
  real(dp)                                       :: r3(3)
  real(dp)                                       :: r33(3,3)
  real(dp),    dimension(:,:), allocatable       :: sum2
  real(dp)                                       :: t1p
  real(dp)                                       :: t2p
  real(dp)                                       :: totalmass
  real(dp)                                       :: w1l(9)
  real(dp)                                       :: xcom
  real(dp)                                       :: ycom
  real(dp)                                       :: zcom
#ifdef TRACE
  call trace_in('property0')
#endif
!
  t1p = g_cpu_time()
!
  n3f = 3*numat
!
!  If pressure is zero (the default) then set to 1 atm
!
  if (press.eq.0.0_dp) then
    press_local = 1.0_dp
  else
    press_local = press/101.325d-6
  endif
!**************************
!  Raman susceptibilites  *
!**************************
  if (lraman.and.nshell.gt.0) then
!
!  Allocate local memory
!
    allocate(qshell(numat),stat=status)
    if (status/=0) call outofmemory('property0','qshell')
    allocate(qD(3*nsfoc,3),stat=status)
    if (status/=0) call outofmemory('property0','qD')
    if (nprocs.gt.1) then
      allocate(sum2(3*nsfoc,3),stat=status)
      if (status/=0) call outofmemory('property0','sum2')
    endif
!
!  Collect shell second derivative terms
!
    msvar = 3*nshell
    do i = 1,nshell
      ni = nshptr(i)
      qshell(i) = qf(ni)*occuf(ni)
      ix = 3*(ni-1) + 1
      iy = ix + 1
      iz = ix + 2
      ixs = 3*(i-1) + 1
      iys = ixs + 1
      izs = ixs + 2
      do j = 1,nshell
        nj = nshptr(j)
        jx = 3*(nj-1) + 1
        jy = jx + 1
        jz = jx + 2
        jxs = 3*(j-1) + 1
        jys = jxs + 1
        jzs = jxs + 2
        dervi(jxs,ixs) = derv2(jx,ix)
        dervi(jys,ixs) = derv2(jy,ix)
        dervi(jzs,ixs) = derv2(jz,ix)
        dervi(jxs,iys) = derv2(jx,iy)
        dervi(jys,iys) = derv2(jy,iy)
        dervi(jzs,iys) = derv2(jz,iy)
        dervi(jxs,izs) = derv2(jx,iz)
        dervi(jys,izs) = derv2(jy,iz)
        dervi(jzs,izs) = derv2(jz,iz)
      enddo
    enddo
    if (nbss.gt.0) then
!
!  Collect radial second derivatives
!
      if (lpocc) then
        do i = 1,nshell
          ni = nshptr(i)
          do j = 1,nshell
            nj = nshptr(j)
            dervi(msvar+j,msvar+i) = derv2(n3f+nj,n3f+ni)
          enddo
          do j = 1,nshell
            nj = nshptr(j)
            jx = 3*(nj-1) + 1
            jy = jx + 1
            jz = jx + 2
            jxs = 3*(j-1) + 1
            jys = jxs + 1
            jzs = jxs + 2
            dervi(jxs,msvar+i) = derv2(jx,n3f+ni)
            dervi(jys,msvar+i) = derv2(jy,n3f+ni)
            dervi(jzs,msvar+i) = derv2(jz,n3f+ni)
            dervi(msvar+i,jxs) = derv2(n3f+ni,jx)
            dervi(msvar+i,jys) = derv2(n3f+ni,jy)
            dervi(msvar+i,jzs) = derv2(n3f+ni,jz)
          enddo
        enddo
      else
        nbsc = nbs - nbss
        do i = 1,nbss
          iptr = n3f + nbsptr(nbsc+i)
          do j = 1,nbss
            jptr = n3f + nbsptr(nbsc+j)
            dervi(msvar+j,msvar+i) = derv2(jptr,iptr)
          enddo
          do j = 1,nshell
            nj = nshptr(j)
            jx = 3*(nj-1) + 1
            jy = jx + 1
            jz = jx + 2
            jxs = 3*(j-1) + 1
            jys = jxs + 1
            jzs = jxs + 2
            dervi(jxs,msvar+i) = derv2(jx,iptr)
            dervi(jys,msvar+i) = derv2(jy,iptr)
            dervi(jzs,msvar+i) = derv2(jz,iptr)
          enddo
        enddo
      endif
    endif
!
!  Compress derivatives for partial occupancies
!
    if (lpocc) then
      call compressd1(qshell,0_i4,nsfoc,nshell,iocshptr)
      call compressd2(dervi,maxd2,0_i4,nsfoc,nbsfoc,nshell,iocshptr,ibocshptr)
      n = 3*nsfoc + nbsfoc
    else
      n = msvar + nbss
    endif
!
!  Invert second derivative matrix
!
    ifail = 1
    call matrix_inversion_library(n,1_i4,maxd2,3_i4*nblocksize,dervi,0_i4,ifail)
!
    if (ifail.ne.0) then
      nwarn = nwarn + 1
      call outwarning('Properties cannot be calculated - matrix is singular',0_i4)
      goto 999
    endif
!
!  Multiply inverse second derivatives by charge vectors:
!
!  First multiply by one charge vector and store in dervi for use in Raman calculation, if needed
!
    do k = 1,nsfoc
      qak = qshell(k)
      indk = 3*(k-1)
      do indl = 1,3*nsfoc
        dervi(indl,indk+1) = dervi(indl,indk+1)*qak
      enddo
      do indl = 1,3*nsfoc
        dervi(indl,indk+2) = dervi(indl,indk+2)*qak
      enddo
      do indl = 1,3*nsfoc
        dervi(indl,indk+3) = dervi(indl,indk+3)*qak
      enddo
    enddo
!
!  Compact dervi down into qD
!
    if (lpocc) then
      qD(1:3*nsfoc,1:3) = 0.0_dp
      do k = 1,nsfoc
        indk = 3*(k-1)
        do i = 1,3
          do indl = 1,3*nsfoc
            qD(indl,i) = qD(indl,i) + dervi(indl,indk+i)
          enddo
        enddo
      enddo
    else
      if (nprocs.gt.1) then
        qD(1:3*nshell,1:3) = 0.0_dp
        do k = 1,nshellonnode
          indk = 3*(k-1)
          do i = 1,3
            do indl = 1,3*nshell
              qD(indl,i) = qD(indl,i) + dervi(indl,indk+i)
            enddo
          enddo
        enddo
        call sumall(qD,sum2,9_i4*nshell,"property0","qD")
        qD(1:3*nshell,1:3) = sum2(1:3*nshell,1:3)
      else
        qD(1:3*nshell,1:3) = 0.0_dp
        do k = 1,nshell
          indk = 3*(k-1)
          do i = 1,3
            do indl = 1,3*nshell
              qD(indl,i) = qD(indl,i) + dervi(indl,indk+i)
            enddo
          enddo
        enddo
      endif
    endif
!
!  Call routine to compute the third derivatives and generate the susceptibility tensors
!
    call raman0(qD,3*nsfoc)
!
!  Deallocate memory
!
    if (nprocs.gt.1) then
      deallocate(sum2,stat=status)
      if (status/=0) call deallocate_error('property0','sum2')
    endif
    deallocate(qD,stat=status)
    if (status/=0) call deallocate_error('property0','qD')
    deallocate(qshell,stat=status)
    if (status/=0) call deallocate_error('property0','qshell')
  endif
999 continue
!*****************************
!  Moment of inertia tensor  *
!*****************************
!
!  Find centre of mass and total mass
!
  xcom = 0.0_dp
  ycom = 0.0_dp
  zcom = 0.0_dp
  totalmass = 0.0_dp
  do i = 1,ncore
    mi = occuf(i)*massspec(nspecptr(i))
    xcom = xcom + xclat(i)*mi
    ycom = ycom + yclat(i)*mi
    zcom = zcom + zclat(i)*mi
    totalmass = totalmass + mi
  enddo
  if (numat.gt.1) then
!
!  Initialise tensor
!
    inertia(1:3,1:3) = 0.0_dp
!
    xcom = xcom/totalmass
    ycom = ycom/totalmass
    zcom = zcom/totalmass
!
!  Calculate tensor
!
    do i = 1,numat
      mi = occuf(i)*massspec(nspecptr(i))
      dx = xclat(i) - xcom
      dy = yclat(i) - ycom
      dz = zclat(i) - zcom
      inertia(1,1) = inertia(1,1) + mi*dy*dy + mi*dz*dz
      inertia(2,1) = inertia(2,1) - mi*dx*dy
      inertia(3,1) = inertia(3,1) - mi*dx*dz
      inertia(1,2) = inertia(1,2) - mi*dy*dx
      inertia(2,2) = inertia(2,2) + mi*dx*dx + mi*dz*dz
      inertia(3,2) = inertia(3,2) - mi*dy*dz
      inertia(1,3) = inertia(1,3) - mi*dz*dx
      inertia(2,3) = inertia(2,3) - mi*dz*dy
      inertia(3,3) = inertia(3,3) + mi*dx*dx + mi*dy*dy
    enddo
!
!  Convert units to 10**-46 x kg.m**2
!
    inertia(1:3,1:3) = inertia(1:3,1:3)*1.0d23/avogadro
!
!  Calculate diagonalised tensor components
!
    do i = 1,3
      do j = 1,3
        r33(j,i) = inertia(j,i)
      enddo
    enddo
    call dsyev('N','U',3_i4,r33,3_i4,r3,w1l,9_i4,ifail)
!
!  Rotational partition function
!
!  If T > 0 then compute thermodynamics too
!
    if (temperature.gt.1.0d-3) then
      if (r3(1).lt.1.0d-12) then
!
!  Linear molecule
!
        inertia_product = r3(3)/16.60538921_dp
        lnZr = log(inertia_product) + log(temperature) - log(dble(symnocfg(ncf))) + 1.418_dp
      else
!
!  Non-linear molecule
!
        inertia_product = r3(1)*r3(2)*r3(3)/(16.60538921_dp)**3
        lnZr = 0.5_dp*log(inertia_product) + 1.5_dp*log(temperature) - log(dble(symnocfg(ncf))) + 1.5_dp*1.418_dp
      endif
      Gr = - boltz*temperature*lnZr/evtoj
    endif
  endif
!
!  Translational partition function
!
!  If T > 0 then compute thermodynamics too
!
  if (temperature.gt.1.0d-3) then
    lnZt = 1.5_dp*log(totalmass) + 2.5_dp*log(temperature) - log(press_local) + 51.104_dp
    Gt = - boltz*temperature*lnZt/evtoj
    Gn =   boltz*temperature*log(avogadro)/evtoj
  endif
!**********************
!  Output properties  *
!**********************
  if (lprint.and.ioproc) then
    if (numat.gt.1) then
      write(ioout,'(/)')
      write(ioout,'(''  Moment of inertia tensor (10^-46 x kgm^2): '',/)')
      write(ioout,'(''-------------------------------------------------------------------------------'')')
      write(ioout,'(10x,6x,''x'',11x,''y'',11x,''z'',11x,''Principal axis system'')')
      write(ioout,'(''-------------------------------------------------------------------------------'')')
      write(ioout,'(7x,''x'',2x,3f12.2,2x,f18.6)')(inertia(1,i),i=1,3),r3(1)
      write(ioout,'(7x,''y'',2x,3f12.2,2x,f18.6)')(inertia(2,i),i=1,3),r3(2)
      write(ioout,'(7x,''z'',2x,3f12.2,2x,f18.6)')(inertia(3,i),i=1,3),r3(3)
      write(ioout,'(''-------------------------------------------------------------------------------'')')
      write(ioout,'(''  Centre of mass (Ang) = '',3(f12.6,1x))') xcom,ycom,zcom
      write(ioout,'(''-------------------------------------------------------------------------------'',/)')
    endif
    if (temperature.gt.1.0d-3) then
      write(ioout,'(''-------------------------------------------------------------------------------'')')
      if (numat.gt.1) then
        write(ioout,'(''  Thermodynamics from rotation and translation at T = '',f12.4,'' K  '')') temperature
        write(ioout,'(''                                              and p = '',f12.4,'' atm  '')') press_local
      else
        write(ioout,'(''  Thermodynamics from translation at T = '',f12.4,'' K  '')') temperature
        write(ioout,'(''                                 and p = '',f12.4,'' atm  '')') press_local
      endif
      write(ioout,'(''-------------------------------------------------------------------------------'')')
      if (numat.gt.1) then
        write(ioout,'(''  Rotational    partition function = '',f18.6)') lnZr
        write(ioout,'(''  Translational partition function = '',f18.6)') lnZt
        write(ioout,'(''  Rotational free energy           = '',f18.6,'' eV'')') Gr
        write(ioout,'(''                                   = '',f18.6,'' kJ/mol'')') Gr/kjmtoev
        write(ioout,'(''  Translational free energy        = '',f18.6,'' eV'')') Gt + Gn
        write(ioout,'(''                                   = '',f18.6,'' kJ/mol'')') (Gt + Gn)/kjmtoev
        write(ioout,'(''  Total free energy                = '',f18.6,'' eV'')') Gt + Gr + Gn
        write(ioout,'(''                                   = '',f18.6,'' kJ/mol'')') (Gt + Gr + Gn)/kjmtoev
      else
        write(ioout,'(''  Translational partition function = '',f18.6)') lnZt
        write(ioout,'(''  Total free energy                = '',f18.6,'' eV'')') Gt + Gn
        write(ioout,'(''                                   = '',f18.6,'' kJ/mol'')') (Gt + Gn)/kjmtoev
      endif
      write(ioout,'(''-------------------------------------------------------------------------------'',/)')
    endif
    if (lcml) call gulp_cml_output_inertia(inertia)
    call gflush(ioout)
  endif
!***************
!  Exit tasks  *
!***************
!
!  Timings
!
  t2p = g_cpu_time()
  tprop = t2p - t1p + tprop
#ifdef TRACE
  call trace_out('property0')
#endif
!
  return
  end
