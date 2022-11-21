  subroutine setdcosmobmat1p(ldqneeded,dcosmoB,dsegweight)
!
!  Subroutine calculates the first derivatives of the COSMO B matrix in parallel
!
!   4/17 Created from setdcosmobmat
!
!  It is assumed that lgrad1 = .true. on entry otherwise routine
!  would not have been called. If lgrad = .true. then second
!  derivatives will be calculated as well.
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
!  Copyright Curtin University 2017
!
!  Julian Gale, CIC, Curtin University, April 2017
!
  use cosmic
  use g_constants
  use control
  use current
  use derivatives
  use optimisation
  use parallel
  use reallocate
  implicit none
!
!  Passed variables
!
  logical,        intent(in)                    :: ldqneeded
  real(dp),       intent(inout)                 :: dcosmoB(3,numat,*)
  real(dp),       intent(in)                    :: dsegweight(3,maxnearseg,*)
!
!  Local variables
!
  integer(i4)                                   :: i
  integer(i4)                                   :: ipts
  integer(i4)                                   :: iptsloc
  integer(i4)                                   :: j
  integer(i4)                                   :: k
  integer(i4)                                   :: kk
  integer(i4)                                   :: n
  real(dp)                                      :: dqme(3)
  real(dp)                                      :: d2qme(6)
  real(dp)                                      :: fact
  real(dp)                                      :: ff
  real(dp)                                      :: ffs
  real(dp)                                      :: qj
  real(dp)                                      :: qme
  real(dp)                                      :: qsk
  real(dp)                                      :: sumAinvipts
  real(dp)                                      :: swi
  real(dp)                                      :: xd, yd, zd
!
!  Set local constants
!
  fact  = 0.5_dp*autoev*autoangs*cosmofneps
!****************************
!  Derivatives of B matrix  *
!****************************
  do iptsloc = 1,nptsonnode
    ipts = node2pts(iptsloc)
    i = cosmoatomptr(ipts)
    swi = segweight(ipts)
    qsk = 2.0_dp*qonsas(ipts) 
    sumAinvipts = 0.0_dp
    do n = 1,nsasparticles
      j = nsasparticleptr(n)
      qj = qsasparticles(n)
      xd = spxyz(1,ipts) - xclat(j)                                                                         
      yd = spxyz(2,ipts) - yclat(j)                                                                   
      zd = spxyz(3,ipts) - zclat(j)
      call qmatrixelement(xd,yd,zd,0.0_dp,.true.,.false.,qme,dqme,d2qme)
      ff = (qsk - deltaq*swi)*fact*qj
      ffs = swi*ff
      if (i.ne.j) then
!
!  First derivatives of B matrix due to B matrix element
!
        xdrv(i) = xdrv(i) - ffs*dqme(1)
        ydrv(i) = ydrv(i) - ffs*dqme(2)
        zdrv(i) = zdrv(i) - ffs*dqme(3)
        xdrv(j) = xdrv(j) + ffs*dqme(1)
        ydrv(j) = ydrv(j) + ffs*dqme(2)
        zdrv(j) = zdrv(j) + ffs*dqme(3)
      endif
      if (lsegsmooth) then
!
!  First derivatives of B matrix due to smoothing term
!
        do kk = 1,nnearseg(ipts)
          k = nnearsegptr(kk,ipts)
          xdrv(i) = xdrv(i) - ff*qme*dsegweight(1,kk,ipts)
          ydrv(i) = ydrv(i) - ff*qme*dsegweight(2,kk,ipts)
          zdrv(i) = zdrv(i) - ff*qme*dsegweight(3,kk,ipts)
          xdrv(k) = xdrv(k) + ff*qme*dsegweight(1,kk,ipts)
          ydrv(k) = ydrv(k) + ff*qme*dsegweight(2,kk,ipts)
          zdrv(k) = zdrv(k) + ff*qme*dsegweight(3,kk,ipts)
        enddo
!
!  Smoothing contribution to dcosmoB
!
        if (ldqneeded) then
          do kk = 1,nnearseg(ipts)
            k = nnearsegptr(kk,ipts)
            dcosmoB(1,k,iptsloc) = dcosmoB(1,k,iptsloc) + qme*qj*dsegweight(1,kk,ipts)
            dcosmoB(2,k,iptsloc) = dcosmoB(2,k,iptsloc) + qme*qj*dsegweight(2,kk,ipts)
            dcosmoB(3,k,iptsloc) = dcosmoB(3,k,iptsloc) + qme*qj*dsegweight(3,kk,ipts)
            dcosmoB(1,i,iptsloc) = dcosmoB(1,i,iptsloc) - qme*qj*dsegweight(1,kk,ipts)
            dcosmoB(2,i,iptsloc) = dcosmoB(2,i,iptsloc) - qme*qj*dsegweight(2,kk,ipts)
            dcosmoB(3,i,iptsloc) = dcosmoB(3,i,iptsloc) - qme*qj*dsegweight(3,kk,ipts)
          enddo
        endif
      endif
!
      if (i.ne.j) then
        if (ldqneeded) then
!
!  Save first derivative terms for use in second derivatives
!
          dcosmoB(1,j,iptsloc) = dcosmoB(1,j,iptsloc) + swi*dqme(1)*qj
          dcosmoB(2,j,iptsloc) = dcosmoB(2,j,iptsloc) + swi*dqme(2)*qj
          dcosmoB(3,j,iptsloc) = dcosmoB(3,j,iptsloc) + swi*dqme(3)*qj
          dcosmoB(1,i,iptsloc) = dcosmoB(1,i,iptsloc) - swi*dqme(1)*qj
          dcosmoB(2,i,iptsloc) = dcosmoB(2,i,iptsloc) - swi*dqme(2)*qj
          dcosmoB(3,i,iptsloc) = dcosmoB(3,i,iptsloc) - swi*dqme(3)*qj
        endif
      endif
    enddo
  enddo
!
  return
  end
