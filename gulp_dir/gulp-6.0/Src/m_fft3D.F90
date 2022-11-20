module m_fft3d
!
!  This module contains the variables and subroutines associated 
!  with the 3-D FFT used in the SPME routines. Everything is based
!  on the use of FFTW to do the FFT.
!
!  Julian Gale, Curtin University, June 2016
!
  use datatypes
  use, intrinsic :: iso_c_binding
  implicit none
#ifdef FFTW3
#ifdef MPI
  include 'fftw3-mpi.f03'
#else
  include 'fftw3.f03'
#endif
#endif
  integer(C_INTPTR_T)                    :: nallocloc
  integer(C_INTPTR_T)                    :: ncfft1
  integer(C_INTPTR_T)                    :: ncfft2
  integer(C_INTPTR_T)                    :: ncfft3
  integer(C_INTPTR_T)                    :: nclfft3
  integer(C_INTPTR_T)                    :: nclfft3_offset
  type(C_PTR)                            :: ccq
  type(C_PTR)                            :: fftplanforward
  type(C_PTR)                            :: fftplanbackward
  complex(C_DOUBLE_COMPLEX), pointer     :: cq(:,:,:)
#ifdef FFTW3
  logical                                :: lfftw3 = .true.
#else
  logical                                :: lfftw3 = .false.
#endif

CONTAINS

  subroutine initfft3D(K1,K2,K3,K3loc,K3start)
#ifdef FFTW3
#ifdef MPI
  use parallel, only : MPI_comm_GULP, nprocs
#endif
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)  :: K1
  integer(i4), intent(in)  :: K2
  integer(i4), intent(in)  :: K3
  integer(i4), intent(out) :: K3loc
  integer(i4), intent(out) :: K3start
!
!  Local variables
!
  integer       :: ifront
  integer       :: iback
  integer       :: iguess
  integer       :: error
!
!  Initialise FFTW
!
  ifront = -1
  iback = 1
  error = 0
  iguess = 0
!
!  Copy dimensions of FFT to FFTW variables
!
  ncfft1 = K1
  ncfft2 = K2
  ncfft3 = K3
#ifdef FFTW3
#ifdef MPI
  if (nprocs.gt.1) then
!
!  Parallel FFT
!
!  Initialise FFTW with MPI
!
    call fftw_mpi_init()
!
!  Get local problem sizes
!
    nallocloc = fftw_mpi_local_size_3d(ncfft3, ncfft2, ncfft1, MPI_comm_GULP, &
                                       nclfft3, nclfft3_offset)
!
!  Create array and Fortran pointer
!
    ccq = fftw_alloc_complex(nallocloc)
    call c_f_pointer(ccq, cq, [ncfft1,ncfft2,nclfft3])
!
    fftplanforward  = fftw_mpi_plan_dft_3d(ncfft3,ncfft2,ncfft1,cq,cq,MPI_comm_GULP,ifront,iguess)
    fftplanbackward = fftw_mpi_plan_dft_3d(ncfft3,ncfft2,ncfft1,cq,cq,MPI_comm_GULP,iback,iguess)
  else
#endif
!
!  Create array and Fortran pointer
!
    nallocloc = ncfft1*ncfft2*ncfft3
    ccq = fftw_alloc_complex(nallocloc)
    call c_f_pointer(ccq, cq, [ncfft1,ncfft2,ncfft3])
    nclfft3 = ncfft3
    nclfft3_offset = 0
!
    fftplanforward  = fftw_plan_dft_3d(K3,K2,K1,cq,cq,ifront,iguess)
    fftplanbackward = fftw_plan_dft_3d(K3,K2,K1,cq,cq,iback,iguess)
#ifdef MPI
  end if
!
!  Set parallel distribution in third direction
!
  K3loc = nclfft3
  K3start = nclfft3_offset
#else
  K3loc = K3
  K3start = 0_i4
#endif
#else
  K3loc = K3
  K3start = 0_i4
#endif
  return
  end subroutine
!
  subroutine fft3Dforward
!
!  Forward 3D-FFT for SPME using FFTW
!
#ifdef FFTW3
#ifdef MPI
  use parallel, only : nprocs
#endif
#endif
  implicit none
#ifdef FFTW3
#ifdef MPI
  if (nprocs.gt.1) then
    call fftw_mpi_execute_dft(fftplanforward,cq,cq)
  else
#endif
    call fftw_execute_dft(fftplanforward,cq,cq)
#ifdef MPI
  endif
#endif
!
#endif
  return
  end subroutine
!
  subroutine fft3Dbackward
!
!  Forward 3D-FFT for SPME using FFTW
!
#ifdef FFTW3
#ifdef MPI
  use parallel, only : nprocs
#endif
#endif
  implicit none
!
!     perform a single 3-D backward transform using FFTW
!
#ifdef FFTW3
#ifdef MPI
  if (nprocs.gt.1) then
    call fftw_mpi_execute_dft(fftplanbackward,cq,cq)
  else
#endif
    call fftw_execute_dft(fftplanbackward,cq,cq)
#ifdef MPI
  endif
#endif
#endif
!
  return
  end subroutine

end module m_fft3d
