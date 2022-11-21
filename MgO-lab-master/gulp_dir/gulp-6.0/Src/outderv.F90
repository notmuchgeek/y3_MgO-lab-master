  subroutine outderv
!
!  Outputs second derivatives for debugging
!
!   5/02 Created from strfin
!   5/09 Output banners modified and space added after strain second derivs
!   7/13 Modified to allow for the fact that the second derivative matrix
!        may not have been symmetrised prior to this call
!   4/17 Further calls to mpbarrier and gflush added
!   3/18 Parallel I/O corrected
!   4/19 Setting of matom changed to allow for freezing
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
!  Julian Gale, CIC, Curtin University, April 2019
!
  use control
  use current
  use derivatives
  use iochannels
  use optimisation,    only : lopf, lfreeze
  use parallel
  use symmetry
  implicit none
#ifdef MPI
  include 'mpif.h'
#endif
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ii
  integer(i4)                                  :: ind
  integer(i4)                                  :: j
  integer(i4)                                  :: maxlim
  integer(i4)                                  :: matom
#ifdef MPI
  integer(i4)                                  :: iloc
  integer(i4)                                  :: status
  integer                                      :: MPIerror
  integer                                      :: ntag
  integer                                      :: nnode
  integer                                      :: ntmp
  integer                                      :: Request
  integer,     dimension(:),   allocatable     :: StatMPI       ! Array for status from MPI
  real(dp),    dimension(:,:), allocatable     :: dtmp
#endif
!
!  Find number of unfrozen atoms
!
  if (lsymderv2) then
    if (lfreeze) then
      matom = 0
      do i = 1,nasym
        if (lopf(i)) matom = matom + 1
      enddo
    else
      matom = nasym
    endif
  else
    if (lfreeze) then
      matom = 0
      do i = 1,nasym
        if (lopf(i)) matom = matom + neqv(i)
      enddo
    else
      matom = numat
    endif
  endif
!
  maxlim = 3*matom
  if (nbsm.gt.0) maxlim = maxlim + matom
!
!  Write out a space of any of the options are to be output
!
  if (ioproc) then
    if (index(keyword,'derv2').ne.0.or.index(keyword,'derv3').ne.0) then
      write(ioout,'(/)')
    endif
  endif
!********************
!  Strain - strain  *
!********************
  if (ioproc.and.index(keyword,'derv2').ne.0.and.lstr) then
    write(ioout,'(''  Strain-strain second derivative matrix : (eV)'',/)')
    do i = 1,nstrains
      write(ioout,'(6f12.5)')(sderv2(j,i),j=1,nstrains)
    enddo
    write(ioout,'(/)')
    call gflush(ioout)
  endif
  if (nprocs.gt.1) then
!**********************
!  Internal - strain  *
!**********************
    if (index(keyword,'derv3').ne.0.and.lstr) then
      call mpbarrier
      if (ioproc) then
        write(ioout,'(''  Mixed strain-internal second derivative matrix : '',''(eV/Angstrom)'',/)')
        call gflush(ioout)
      endif
#ifdef MPI
      if (lioproconly) then
!
!  Allocate temporary workspace for communication
!
        ntmp = 3*nstrains
        ntag = 1
        allocate(dtmp(nstrains,3_i4),stat=status)
        if (status/=0) call outofmemory('outderv','dtmp')
        allocate(StatMPI(MPI_Status_Size),stat=status)
        if (status/=0) call outofmemory('outderv','StatMPI')
!
        do i = 1,numat
          iloc = atom2local(i)
          if (atom2node(i).ne.0_i4) then
!
!  Post receive
!
            if (ioproc) then
              nnode = atom2node(i)
              call MPI_IRecv(dtmp,ntmp,MPI_double_precision,nnode, &
                             ntag,MPI_Comm_World,Request,MPIerror)
            endif
!
!  Pass data to ioproc for writing
!
            if (iloc.gt.0) then
              ind = 3*(atom2local(i)-1)
              do ii = 1,3
                do j = 1,nstrains
                  dtmp(j,ii) = derv3(ind+ii,j)
                enddo
              enddo
!
!  Post send
!
              call MPI_ISend(dtmp,ntmp,MPI_double_precision,0, &
                             ntag,MPI_Comm_World,Request,MPIerror)
            endif
            if (ioproc.or.iloc.gt.0) then
              call MPI_WaitAll(1,Request,StatMPI,MPIerror)
            endif
            if (ioproc) then
!
!  Write on I/O node
!
              do ii = 1,3
                write(ioout,'(9f12.5)')(dtmp(j,ii),j=1,nstrains)
              enddo
            endif
          else
            if (iloc.gt.0) then
              ind = 3*(atom2local(i)-1)
              do ii = 1,3
                write(ioout,'(9f12.5)')(derv3(ind+ii,j),j=1,nstrains)
              enddo
            endif
          endif
        enddo
!
        deallocate(StatMPI,stat=status)
        if (status/=0) call deallocate_error('outderv','StatMPI')
        deallocate(dtmp,stat=status)
        if (status/=0) call deallocate_error('outderv','dtmp')
      else
#endif
        call mpbarrier
        do i = 1,numat
          if (procid.eq.atom2node(i)) then
            ind = 3*(atom2local(i)-1)
            do ii = 1,3
              write(ioout,'(9f12.5)')(derv3(ind+ii,j),j=1,nstrains)
            enddo
            call gflush(ioout)
          endif
          call mpbarrier
        enddo
#ifdef MPI
      endif
      call mpbarrier
#endif
      if (ioproc) then
        write(ioout,'(/)')
        call gflush(ioout)
      endif
    endif
!************************
!  Internal - internal  *
!************************
    if (index(keyword,'derv2').ne.0) then
      call mpbarrier
      if (ioproc) then
        write(ioout,'(''  Internal-internal second derivative matrix : '',''(eV/Angstrom**2)'',/)')
        call gflush(ioout)
      endif
#ifdef MPI
      if (lioproconly) then
!
!  Allocate temporary workspace for communication
!
        ntmp = 3*maxlim
        ntag = 1
        allocate(dtmp(maxlim,3_i4),stat=status)
        if (status/=0) call outofmemory('outderv','dtmp')
        allocate(StatMPI(MPI_Status_Size),stat=status)
        if (status/=0) call outofmemory('outderv','StatMPI')
!
        do i = 1,numat
          iloc = atom2local(i)
          if (atom2node(i).ne.0_i4) then
!
!  Post receive
!
            if (ioproc) then
              nnode = atom2node(i)
              call MPI_IRecv(dtmp,ntmp,MPI_double_precision,nnode, &
                             ntag,MPI_Comm_World,Request,MPIerror)
            endif
!
!  Pass data to ioproc for writing
!
            if (iloc.gt.0) then
              ind = 3*(atom2local(i)-1)
              do ii = 1,3
                do j = 1,maxlim
                  dtmp(j,ii) = derv2(j,ind+ii)
                enddo
              enddo
!
!  Post send
!
              call MPI_ISend(dtmp,ntmp,MPI_double_precision,0, &
                             ntag,MPI_Comm_World,Request,MPIerror)
            endif
            if (ioproc.or.iloc.gt.0) then
              call MPI_WaitAll(1,Request,StatMPI,MPIerror)
            endif
            if (ioproc) then
!
!  Write on I/O node
!
              do ii = 1,3
                write(ioout,'(9f12.5)')(dtmp(j,ii),j=1,maxlim)
              enddo
            endif
          else
            if (iloc.gt.0) then
              ind = 3*(atom2local(i)-1)
              do ii = 1,3
                write(ioout,'(9f12.5)')(derv2(j,ind+ii),j=1,maxlim)
              enddo
            endif
          endif
        enddo
!
        deallocate(StatMPI,stat=status)
        if (status/=0) call deallocate_error('outderv','StatMPI')
        deallocate(dtmp,stat=status)
        if (status/=0) call deallocate_error('outderv','dtmp')
      else
#endif
        call mpbarrier
        do i = 1,numat
          if (procid.eq.atom2node(i)) then
            ind = 3*(atom2local(i)-1)
            do ii = 1,3
              write(ioout,'(9f12.5)')(derv2(j,ind+ii),j=1,maxlim)
            enddo
            call gflush(ioout)
          endif
          call mpbarrier
        enddo
        call mpbarrier
#ifdef MPI
      endif
#endif
      if (ioproc) then
        write(ioout,'(/)')
        call gflush(ioout)
      endif
    endif
  else
!**********************
!  Internal - strain  *
!**********************
    if (index(keyword,'derv3').ne.0.and.lstr) then
      write(ioout,'(''  Mixed strain-internal second derivative matrix : '',''(eV/Angstrom)'',/)')
      do i = 1,maxlim
        write(ioout,'(9f12.5)')(derv3(i,j),j=1,nstrains)
      enddo
      write(ioout,'(/)')
    endif
!************************
!  Internal - internal  *
!************************
    if (index(keyword,'derv2').ne.0) then
      write(ioout,'(''  Internal-internal second derivative matrix : '',''(eV/Angstrom**2)'',/)')
      do i = 1,maxlim
        write(ioout,'(9f12.5)')(derv2(j,i),j=1,maxlim)
      enddo
      write(ioout,'(/)')
    endif
  endif
!
  return
  end
