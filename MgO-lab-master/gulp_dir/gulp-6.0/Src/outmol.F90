  subroutine outmol
!
!  Outputs details of molecules
!
!   3/95 Dimensionality of molecule added to output
!  10/01 Format of print changed
!   6/07 Format of output revised slightly
!   4/08 Format for atom number changed to i5 to handle larger systems
!   6/09 Modified to use new atom in molecule pointers
!   4/20 Debug printing added for rigid molecules
!   6/20 Output of multipole moments for rigid molecules added
!   6/20 nmolcore changes added
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
!  Julian Gale, CIC, Curtin University, June 2020
!
  use control,        only : lrigid, ldebug
  use current
  use element
  use iochannels
  use molecule
  implicit none
!
!  Local variables
!
  character(len=1), dimension(:), allocatable       :: char
  character(len=5), dimension(:), allocatable       :: cha
  character(len=2)                                  :: dchar(6)
  character(len=5)                                  :: lab
  character(len=5)                                  :: lab1
  character(len=5)                                  :: lab2
  integer(i4)                                       :: i
  integer(i4)                                       :: j
  integer(i4)                                       :: jj
  integer(i4)                                       :: k
  integer(i4)                                       :: n
  integer(i4)                                       :: n1
  integer(i4)                                       :: n2
  integer(i4)                                       :: na
  integer(i4)                                       :: nend
  integer(i4)                                       :: nhigh
  integer(i4)                                       :: ninc
  integer(i4)                                       :: nline
  integer(i4)                                       :: nlow
  integer(i4)                                       :: nm
  integer(i4)                                       :: ntyp
  integer(i4)                                       :: ntvar1
  integer(i4)                                       :: ntvar2
  integer(i4)                                       :: nvar1
  integer(i4)                                       :: nvar2
  integer(i4)                                       :: status
  real(dp)                                          :: qi
  real(dp)                                          :: r2
  real(dp)                                          :: rmdipole(3)
  real(dp)                                          :: rmquad(6)
  real(dp)                                          :: xi
  real(dp)                                          :: yi
  real(dp)                                          :: zi
!
  data dchar/'x ','y ','z ','xy','xz','yz'/
!
  write(ioout,'(/,''  Molecule list generated from bond lengths :'',/)')
  if (nnobo.gt.0) then
    write(ioout,'(''  Bonds between the following atom types will be excluded : '',/)')
    do i = 1,nnobo
      nvar1 = nobond(i)/1000
      nvar2 = nobond(i) - 1000*nvar1
      if (nvar1.gt.maxele) then
        n1 = nvar1 - maxele
      else
        n1 = nvar1
      endif
      if (nvar2.gt.maxele) then
        n2 = nvar2 - maxele
      else
        n2 = nvar2
      endif
      ntvar1 = nobotyp(i)/1000
      ntvar2 = nobotyp(i) - 1000*ntvar1
      call label(n1,ntvar1,lab1)
      call label(n2,ntvar2,lab2)
      if (nvar1.gt.maxele.and.nvar2.gt.maxele) then
        write(ioout,'(3x,a5,1x,''Shell'',2x,''-'',1x,a5,1x,''Shell'')') lab1,lab2
      elseif (nvar1.gt.maxele.and.nvar2.le.maxele) then
        write(ioout,'(3x,a5,1x,''Shell'',2x,''-'',1x,a5,1x,''Core'')') lab1,lab2
      elseif (nvar1.le.maxele.and.nvar2.gt.maxele) then
        write(ioout,'(3x,a5,1x,''Core'',3x,''-'',1x,a5,1x,''Shell'')') lab1,lab2
      else
        write(ioout,'(3x,a5,1x,''Core'',3x,''-'',1x,a5,1x,''Core'')') lab1,lab2
      endif
    enddo
    write(ioout,'(/)')
  endif
!
!  Output molecule details
!
  write(ioout,'(''  Total number of molecules = '',i5,/)') nmol
  if (nmol.gt.0) then
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    write(ioout,'(''Molecule No./:  Atoms'')')
    write(ioout,'(''Periodicity  :  '')')
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    allocate(cha(numat),stat=status)
    if (status/=0) call outofmemory('outmol','cha')
    allocate(char(numat),stat=status)
    if (status/=0) call outofmemory('outmol','char')
    do i = 1,nmol
!
!  Loop over atoms in molecule number i
!
      ninc = 0
      do jj = 1,nmolatom(i)
        ninc = ninc + 1
        j = nmollist(nmolptr(i)+jj)
        na = nat(j)
        ntyp = nftype(j)
        call label(na,ntyp,lab)
        if (na.gt.maxele) then
          char(ninc) = 's'
          na = na - maxele
        else
          char(ninc) = 'c'
        endif
        cha(ninc) = lab
      enddo
      nline = (ninc-1)/5 + 1
      nend = min(ninc,5)
      if (moldim(i).eq.1.or.moldim(i).eq.2) then
        write(ioout,'(1x,i5,2x,i1,1x,a2,1x,'':'',5(1x,a5,a1,1x,i5))') &
          i,moldim(i),dchar(3*(moldim(i)-1)+moldimi(i)),(cha(j),char(j),nmollist(nmolptr(i)+j),j=1,nend)
      else
        write(ioout,'(1x,i5,2x,i1,1x,''   :'',5(1x,a5,a1,1x,i5))') i,moldim(i),(cha(j),char(j),nmollist(nmolptr(i)+j),j=1,nend)
      endif
      if (nline.gt.2) then
        do k = 1,nline - 2
          nlow = 5*k + 1
          nhigh = 5*(k + 1)
          write(ioout,'(13x,'':'',5(1x,a5,a1,1x,i5))') (cha(j),char(j),nmollist(nmolptr(i)+j),j=nlow,nhigh)
        enddo
      endif
      if (nline.ge.2) then
        nlow = 5*(nline - 1) + 1
        write(ioout,'(13x,'':'',5(1x,a5,a1,1x,i5))') (cha(j),char(j),nmollist(nmolptr(i)+j),j=nlow,ninc)
      endif
    enddo
    write(ioout,'(''--------------------------------------------------------------------------------'',/)')
    deallocate(char,stat=status)
    if (status/=0) call deallocate_error('outmol','char')
    deallocate(cha,stat=status)
    if (status/=0) call deallocate_error('outmol','cha')
!
!  Rigid molecules
!
    if (lrigid.and.ldebug) then
      write(ioout,'(''  Rigid molecule information: '',/)')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      do nm = 1,nmol
        write(ioout,'(''  Molecule = '',i6,'' : Quaternions = '',3f12.6)') nm,(molQ(j,nm),j=1,3)
        rmdipole(1:3) = 0.0_dp
        rmquad(1:6) = 0.0_dp
        do n = 1,nmolcore(nm)
          i = nmollist(nmolptr(nm)+n)
          write(ioout,'(''      Atom = '',i6,'' : Ref. coord   = '',3f12.6)') i,(molQxyz(j,n,nm),j=1,3)
!
          qi = qf(i)
          xi = molQxyz(1,n,nm)
          yi = molQxyz(2,n,nm)
          zi = molQxyz(3,n,nm)
          r2 = xi*xi + yi*yi + zi*zi
!
          rmdipole(1) = rmdipole(1) + qi*xi
          rmdipole(2) = rmdipole(2) + qi*yi
          rmdipole(3) = rmdipole(3) + qi*zi
!
          rmquad(1) = rmquad(1) + qi*(3.0_dp*xi*xi - r2)
          rmquad(2) = rmquad(2) + qi*(3.0_dp*xi*yi)
          rmquad(3) = rmquad(3) + qi*(3.0_dp*xi*zi)
          rmquad(4) = rmquad(4) + qi*(3.0_dp*yi*yi - r2)
          rmquad(5) = rmquad(5) + qi*(3.0_dp*yi*zi)
          rmquad(6) = rmquad(6) + qi*(3.0_dp*zi*zi - r2)
        enddo
        rmquad(1:6) = 0.5_dp*rmquad(1:6)
!
        write(ioout,'(''      Dipole     moment : '','' x  = '',f12.6,'' y  = '',f12.6,'' z  = '',f12.6)') rmdipole(1:3)
        write(ioout,'(''      Quadrupole moment : '','' xx = '',f12.6,'' xy = '',f12.6,'' xz = '',f12.6)') &
                      rmquad(1),rmquad(2),rmquad(3)
        write(ioout,'(''      Quadrupole moment : '','' yx = '',f12.6,'' yy = '',f12.6,'' yz = '',f12.6)') &
                      rmquad(2),rmquad(4),rmquad(5)
        write(ioout,'(''      Quadrupole moment : '','' zx = '',f12.6,'' zy = '',f12.6,'' zz = '',f12.6)') &
                      rmquad(3),rmquad(5),rmquad(6)
      enddo
      write(ioout,'(''--------------------------------------------------------------------------------'',/)')
    endif
  endif
!
  return
  end
