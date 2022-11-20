  subroutine phoncopy1(derv2,dervi,cmat,maxd2,mci,mcj,maxd2c,msv)
  use datatypes
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
  integer(i4)  :: maxd2
  integer(i4)  :: maxd2c
  integer(i4)  :: mci
  integer(i4)  :: mcj
  integer(i4)  :: msv
  real(dp)     :: derv2(maxd2,*),dervi(maxd2,*)
  complex(dpc) :: cmat(maxd2c,*)
!
  integer(i4)  :: i
  integer(i4)  :: j
#ifdef TRACE
  call trace_in('phoncopy1')
#endif
!
  do i = 1,msv
    do j = 1,msv
      cmat(j,i) = dcmplx(derv2(mci+j,mcj+i),dervi(mci+j,mcj+i))
    enddo
  enddo
#ifdef TRACE
  call trace_out('phoncopy1')
#endif
  return
  end
!
  subroutine phoncopy1r(derv2,dmat,maxd2,mcv,msv)
  use datatypes
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
  integer(i4) :: maxd2
  integer(i4) :: mcv
  integer(i4) :: msv
  real(dp)    :: derv2(maxd2,*),dmat(msv,*)
!
  integer(i4) :: i
  integer(i4) :: j
#ifdef TRACE
  call trace_in('phoncopy1r')
#endif
!
  do i = 1,msv
    do j = 1,msv
      dmat(j,i) = derv2(mcv+j,mcv+i)
    enddo
  enddo
#ifdef TRACE
  call trace_out('phoncopy1r')
#endif
  return
  end
