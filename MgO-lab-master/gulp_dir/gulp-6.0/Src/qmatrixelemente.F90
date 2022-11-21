  subroutine qmatrixelemente(xji,yji,zji,cuts2,lgrad1,lgrad2,qme,dqme,d2qme)
!
!  Routine calculates the long range sum matrix element
!  between two points whose coordinates are supplied.
!  For 3-D this invokes the Ewald sum, while for 2-D it
!  uses the Parry sum. The reciprocal space terms are
!  assumed to have been initialised by a call to a
!  previous routine.
!
!  Version for calling from epot when reciprocal space terms in 3D have
!  been computed by an alternative algorithm that is faster.
!
!   7/20 Created from qmatrixelement
!
!  On entry:
!
!    xji    = difference in X coordinates of two points
!    yji    = difference in Y coordinates of two points
!    zji    = difference in Z coordinates of two points
!    cuts2  = core-shell cutoff, if applicable
!    lgrad1 = if .true. then calculate the first derivatives
!    lgrad2 = if .true. then calculate the second derivatives
!
!  On exit:
!
!    qme    = Coulomb matrix element between points (in Angs**-1)
!    dqme   = first Cartesian derivatives of qme (dQ/d(alpha) in
!             Angs**-2) if lgrad1 is .true.
!    d2qme  = second Cartesian derivatives of qme (d2Q/d(alpha)d(beta)
!             in Angs**-3) if lgrad2 is .true.
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
!  Copyright Curtin Univerisity 2020
!
!  Julian Gale, CIC, Curtin University, July 2020
!
  use control
  use current
  use general,       only : cutw, etaw, rkw, selfwolf, selfwolfforce, selfwolfenergy
  use kspace
  use qmedata
  use shells
  use symmetry,      only : lra
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  real(dp), intent(in)    :: xji
  real(dp), intent(in)    :: yji
  real(dp), intent(in)    :: zji
  real(dp), intent(in)    :: cuts2
  real(dp), intent(out)   :: qme
  real(dp), intent(out)   :: dqme(3)
  real(dp), intent(out)   :: d2qme(6)
  logical,  intent(in)    :: lgrad1
  logical,  intent(in)    :: lgrad2
!
!  Local variables
!
  integer(i4)             :: ii
  integer(i4)             :: iim
  integer(i4)             :: iimn
  integer(i4)             :: iimx
  integer(i4)             :: jj
  integer(i4)             :: jjm
  integer(i4)             :: jjmn
  integer(i4)             :: jjmx
  integer(i4)             :: kk
  integer(i4)             :: kkm
  integer(i4)             :: kkmn
  integer(i4)             :: kkmx
  integer(i4)             :: iv
  integer(i4)             :: ml1
  integer(i4)             :: ml2
  integer(i4)             :: ml3
  real(dp)                :: accf2
  real(dp)                :: arg
  real(dp)                :: argtest
  real(dp)                :: cosa
  real(dp)                :: darg1
  real(dp)                :: darg2
  real(dp)                :: derfez
  real(dp)                :: derfc1
  real(dp)                :: derfc2
  real(dp)                :: dexp1
  real(dp)                :: dexp2
  real(dp)                :: dexp3
  real(dp)                :: dexp4
  real(dp)                :: dexpz
  real(dp)                :: dtrm
  real(dp)                :: dtrm2
  real(dp)                :: dtrm2w
  real(dp)                :: dtrmw
  real(dp)                :: etaloc
  real(dp)                :: etaz
  real(dp)                :: etaz2
  real(dp)                :: Gmax
  real(dp)                :: kexperfc
  real(dp)                :: kvec
  real(dp)                :: r2
  real(dp)                :: rexp
  real(dp)                :: rexpw
  real(dp)                :: rkw2
  real(dp)                :: rl
  real(dp)                :: rrl
  real(dp)                :: rrl2
  real(dp)                :: rr2
  real(dp)                :: rxi, ryi, rzi
  real(dp)                :: rxj, ryj, rzj
  real(dp)                :: rxk, ryk, rzk
  real(dp)                :: setaloc
  real(dp)                :: sina
  real(dp)                :: smallestG
  real(dp)                :: trm
  real(dp)                :: xcd, ycd, zcd
  real(dp)                :: xij
  real(dp)                :: yji2
  real(dp)                :: zji2
  real(dp)                :: ztrm
  logical                 :: lreciploc
!
!  Functions
!
  real(dp)    :: g_derf
  real(dp)    :: g_derfc
#ifdef TRACE
  call trace_in('qmatrixelemente')
#endif
!
!  Zero matrix element / derivatives
!
  qme = 0.0_dp
  if (lgrad1) then
    dqme(1:3) = 0.0_dp
    if (lgrad2) then
      d2qme(1:6) = 0.0_dp
    endif
  endif
!
!  Set up local variables
!
  lreciploc = (.not.lnorecip.and.ndim.gt.1.and..not.lwolf)
  if (lwolf) then
    etaloc = etaw*etaw
    setaloc = etaw
    if (lwolforiginal) then
      if (lgrad1) then
        dtrmw = selfwolfforce
        if (lgrad2) then
          rkw2 = rkw*rkw
          rexpw = tweatpi*exp(-etaloc*cutw*cutw)
          dtrm2w = (rexpw*((3.0_dp*rkw2)+2.0_dp*etaloc) + 3.0_dp*selfwolf*rkw2)*rkw2
        endif
      endif
    elseif (lwolffennell) then
      if (lgrad1) then
        dtrmw = selfwolfforce
        if (lgrad2) then
          dtrm2w = 0.0_dp
        endif
      endif
    endif
  elseif (lewald) then
    etaloc = eta
    setaloc = seta
  endif
!*********************************************
!  Reciprocal space contribution in 2D only  *
!*********************************************
  if (lreciploc) then
    if (ndim.eq.2) then
!
!  2-D sum in reciprocal space
!
      accf2 = accf*accf
      argtest = sqrt(3.0_dp+0.5_dp*accf2) - sqrt(3.0_dp)
      smallestG = min(abs(kv(1,1)),abs(kv(2,2)))
      etaz = seta*zji
      etaz2 = etaz * etaz
      derfez = g_derf(etaz)
      dexpz  = exp(-etaz2)
      qme = qme - vol4pi*(zji*derfez + dexpz*rpieta)
      if (lgrad1) then
        dqme(3) = dqme(3) + vol4pi*derfez
        if (lgrad2) then
          d2qme(6) = d2qme(6) + vol4pi*tweatpi*dexpz
        endif
      endif
!
!  Find local kvector cut-off
!
      if (abs(etaz).gt.argtest) then
        Gmax = abs(accf2/zji)
      else
        Gmax = sqrt(4.0_dp*eta*(accf2-etaz2))
      endif
      if (Gmax.ge.smallestG) then
        iv = 1
        kvec = kmod(iv)
        do while (iv.le.nkvec.and.kvec.le.Gmax)
          arg = xrk(iv)*xji + yrk(iv)*yji
          cosa = cos(arg)*ktrm(iv)
          kvec = kmod(iv)
          dexp1 = exp(kvec*zji)
          dexp2 = 1.0_dp/dexp1
          darg1 = kvec*rhseta + etaz
          darg2 = kvec*rhseta - etaz
          derfc1 = g_derfc(darg1)
          derfc2 = g_derfc(darg2)
          kexperfc = dexp1*derfc1 + dexp2*derfc2
          qme = qme + cosa*kexperfc
          if (lgrad1) then
            sina = sin(arg)*ktrm(iv)
            dexp3 = exp(-darg1**2)
            dexp4 = exp(-darg2**2)
            ztrm  = (kvec*(dexp1*derfc1-dexp2*derfc2) - tweatpi*(dexp1*dexp3-dexp2*dexp4))
            dqme(1) = dqme(1) + sina*kexperfc*xrk(iv)
            dqme(2) = dqme(2) + sina*kexperfc*yrk(iv)
            dqme(3) = dqme(3) - cosa*ztrm
            if (lgrad2) then
              d2qme(1) = d2qme(1) + cosa*xrk(iv)*xrk(iv)*kexperfc
              d2qme(2) = d2qme(2) + cosa*xrk(iv)*yrk(iv)*kexperfc
              d2qme(3) = d2qme(3) + cosa*yrk(iv)*yrk(iv)*kexperfc
              d2qme(4) = d2qme(4) + xrk(iv)*sina*ztrm
              d2qme(5) = d2qme(5) + yrk(iv)*sina*ztrm
              d2qme(6) = d2qme(6) - cosa*(kvec*kvec*kexperfc - &
                2.0_dp*tweatpi*(kvec*(dexp1*dexp3+dexp2*dexp4) - &
                seta*(darg1*dexp1*dexp3+darg2*dexp2*dexp4)))                                       
            endif
          endif
          iv = iv + 1
          kvec = kmod(iv)
        enddo
      endif
    endif
  endif
!*************************
!  Real space summation  *
!*************************
  if (.not.lnoreal) then
    if (ndim.gt.1.or.lwolf) then
      if (lminimage) then
!****************************
!  Minimum image algorithm  *
!****************************
!
!  Set up looping extents
!
        if (ndim.eq.3) then
          iimn = -1
          iimx =  1
          jjmn = -1
          jjmx =  1
          kkmn = -1
          kkmx =  1
        elseif (ndim.eq.2) then
          iimn = -1
          iimx =  1
          jjmn = -1
          jjmx =  1
          kkmn =  0
          kkmx =  0
        elseif (ndim.eq.1) then
          iimn = -1
          iimx =  1
          jjmn =  0
          jjmx =  0
          kkmn =  0
          kkmx =  0
        endif
!
!  Find minimum image
!
        if (lra) then
          rr2 = 1.0d10
          rxi = xji + (iimn-1)*r1x
          do ii = iimn,iimx
            rxi = rxi + r1x
            ryi = yji + (jjmn-1)*r2y
            do jj = jjmn,jjmx
              ryi = ryi + r2y
              rzi = zji + (kkmn-1)*r3z
              do kk = kkmn,kkmx
                rzi = rzi + r3z
                r2 = rxi*rxi + ryi*ryi + rzi*rzi
                if (r2.lt.rr2) then
                  rr2 = r2
                  rxk = rxi
                  ryk = ryi
                  rzk = rzi
                  iim = ii
                  jjm = jj
                  kkm = kk
                endif
              enddo
            enddo
          enddo
        else
          rr2 = 1.0d10
          rxi = xji + (iimn-1)*r1x
          ryi = yji + (iimn-1)*r1y
          rzi = zji + (iimn-1)*r1z
          do ii = iimn,iimx
            rxi = rxi + r1x
            ryi = ryi + r1y
            rzi = rzi + r1z
            rxj = rxi + (jjmn-1)*r2x
            ryj = ryi + (jjmn-1)*r2y
            rzj = rzi + (jjmn-1)*r2z
            do jj = jjmn,jjmx
              rxj = rxj + r2x
              ryj = ryj + r2y
              rzj = rzj + r2z
              xcd = rxj + (kkmn-1)*r3x
              ycd = ryj + (kkmn-1)*r3y
              zcd = rzj + (kkmn-1)*r3z
              do kk = kkmn,kkmx
                xcd = xcd + r3x
                ycd = ycd + r3y
                zcd = zcd + r3z
                r2 = xcd*xcd + ycd*ycd + zcd*zcd
                if (r2.lt.rr2) then
                  rr2 = r2
                  rxk = xcd
                  ryk = ycd
                  rzk = zcd
                  iim = ii
                  jjm = jj
                  kkm = kk
                endif
              enddo
            enddo
          enddo
        endif
!
!
!  Exclude distances outside maximum cutoff
!
        if (rr2.le.rmax2) then
!
!  Trap self term
!
          if (rr2.lt.1.0d-15) then
            qme = qme - tweatpi
            if (lwolf) then
              qme = qme - selfwolf - selfwolfenergy*cutw
            endif
            if (lgrad2) then
              trm = 2.0_dp*tweatpi*etaloc/3.0_dp
              d2qme(1) = d2qme(1) - trm
              d2qme(3) = d2qme(3) - trm
              d2qme(6) = d2qme(6) - trm
            endif
          else
            rl = sqrt(rr2)
            rrl = 1.0_dp/rl
            if (rr2.lt.cuts2) then
!
!  Core-shell interaction
!
              qme = qme - rrl
              if (lgrad1) then
                dtrm = - rrl*rrl*rrl
                if (lgrad2) then
                  dtrm2 = - 3.0_dp*rrl*rrl*rrl
                endif
              endif
            else
              dtrm = 0.0_dp
              dtrm2 = 0.0_dp
            endif
            trm = g_derfc(setaloc*rl)*rrl
            qme = qme + trm
            if (lwolf) then
              qme = qme - selfwolf + selfwolfenergy*(rl - cutw)
            endif
            if (lgrad1) then
              rexp = tweatpi*exp(-etaloc*rr2)
              rrl2 = rrl*rrl
              dtrm = dtrm + (trm + rexp)*rrl2
              if (lwolforiginal.or.lwolffennell) then
                dtrm = dtrm - dtrmw
              endif
              dqme(1) = dqme(1) + dtrm*rxk
              dqme(2) = dqme(2) + dtrm*ryk
              dqme(3) = dqme(3) + dtrm*rzk
              if (lgrad2) then
                dtrm2 = dtrm2 + rexp*((3.0_dp*rrl2)+2.0_dp*etaloc) + 3.0_dp*trm*rrl2
                dtrm2 = dtrm2*rrl2
                if (lwolforiginal) then
                  dtrm2 = dtrm2 - dtrm2w
                endif
                d2qme(1) = d2qme(1) - dtrm2*rxk*rxk
                d2qme(2) = d2qme(2) - dtrm2*rxk*ryk
                d2qme(3) = d2qme(3) - dtrm2*ryk*ryk
                d2qme(4) = d2qme(4) - dtrm2*rxk*rzk
                d2qme(5) = d2qme(5) - dtrm2*ryk*rzk
                d2qme(6) = d2qme(6) - dtrm2*rzk*rzk
                d2qme(1) = d2qme(1) + dtrm
                d2qme(3) = d2qme(3) + dtrm
                d2qme(6) = d2qme(6) + dtrm
              endif
            endif
          endif
        endif
      else
!******************************
!  Periodic 2-D and 3-D case  *
!******************************
        ml1 = maxloop(1)
        ml2 = maxloop(2)
        ml3 = maxloop(3)
!
!  Loop over cell vectors
!
        rxi = xji - (ml1+1)*r1x
        ryi = yji - (ml1+1)*r1y
        rzi = zji - (ml1+1)*r1z
        do ii = -ml1,ml1
          rxi = rxi + r1x
          ryi = ryi + r1y
          rzi = rzi + r1z
          rxj = rxi - (ml2+1)*r2x
          ryj = ryi - (ml2+1)*r2y
          rzj = rzi - (ml2+1)*r2z
          do jj = -ml2,ml2
            rxj = rxj + r2x
            ryj = ryj + r2y
            rzj = rzj + r2z
            rxk = rxj - (ml3+1)*r3x
            ryk = ryj - (ml3+1)*r3y
            rzk = rzj - (ml3+1)*r3z
            do kk = -ml3,ml3
              rxk = rxk + r3x
              ryk = ryk + r3y
              rzk = rzk + r3z
!
!  Calculate distance squared
!
              rr2 = rxk*rxk + ryk*ryk + rzk*rzk
!
!  Exclude distances outside maximum cutoff
!
              if (rr2.le.rmax2) then
!
!  Trap self term
!
                if (rr2.lt.1.0d-15) then
                  qme = qme - tweatpi
                  if (lwolf) then
                    qme = qme - selfwolf - selfwolfenergy*cutw
                  endif
                  if (lgrad2) then
                    trm = 2.0_dp*tweatpi*etaloc/3.0_dp
                    d2qme(1) = d2qme(1) - trm
                    d2qme(3) = d2qme(3) - trm
                    d2qme(6) = d2qme(6) - trm
                  endif
                else
                  rl = sqrt(rr2)
                  rrl = 1.0_dp/rl
                  if (rr2.lt.cuts2) then
!
!  Core-shell interaction
!
                    qme = qme - rrl
                    if (lgrad1) then
                      dtrm = - rrl*rrl*rrl
                      if (lgrad2) then
                        dtrm2 = - 3.0_dp*rrl*rrl*rrl
                      endif
                    endif
                  else
                    dtrm = 0.0_dp
                    dtrm2 = 0.0_dp
                  endif
                  trm = g_derfc(setaloc*rl)*rrl
                  qme = qme + trm
                  if (lwolf) then
                    qme = qme - selfwolf + selfwolfenergy*(rl - cutw)
                  endif
                  if (lgrad1) then
                    rexp = tweatpi*exp(-etaloc*rr2)
                    rrl2 = rrl*rrl
                    dtrm = dtrm + (trm + rexp)*rrl2
                    if (lwolforiginal.or.lwolffennell) then
                      dtrm = dtrm - dtrmw
                    endif
                    dqme(1) = dqme(1) + dtrm*rxk
                    dqme(2) = dqme(2) + dtrm*ryk
                    dqme(3) = dqme(3) + dtrm*rzk
                    if (lgrad2) then
                      dtrm2 = dtrm2 + rexp*((3.0_dp*rrl2)+2.0_dp*etaloc) + 3.0_dp*trm*rrl2
                      dtrm2 = dtrm2*rrl2
                      if (lwolforiginal) then
                        dtrm2 = dtrm2 - dtrm2w
                      endif
                      d2qme(1) = d2qme(1) - dtrm2*rxk*rxk
                      d2qme(2) = d2qme(2) - dtrm2*rxk*ryk
                      d2qme(3) = d2qme(3) - dtrm2*ryk*ryk
                      d2qme(4) = d2qme(4) - dtrm2*rxk*rzk
                      d2qme(5) = d2qme(5) - dtrm2*ryk*rzk
                      d2qme(6) = d2qme(6) - dtrm2*rzk*rzk
                      d2qme(1) = d2qme(1) + dtrm
                      d2qme(3) = d2qme(3) + dtrm
                      d2qme(6) = d2qme(6) + dtrm
                    endif
                  endif
                endif
              endif
!
!  End of loops over lattice vectors
!
            enddo
          enddo
        enddo
      endif
    else
!*************************
!  Cluster and 1-D case  *
!*************************
      if (ndim.eq.1) then
        ml1 = maxloop(1)
        call qmatrix1D(xji,yji,zji,lgrad1,lgrad2,qme,dqme,d2qme)
      else
        ml1 = 0
        r1x = 0.0_dp
      endif
      xij = xji - (ml1+1)*r1x
      yji2 = yji*yji
      zji2 = zji*zji
      do ii = -ml1,ml1
        xij = xij + r1x
        rr2 = xij*xij + yji2 + zji2
!
!  Exclude core-shell interaction
!
        if (rr2.ge.cuts2.and.rr2.ge.1.0d-10) then
          rl = sqrt(rr2)
          rrl = 1.0_dp/rl
          qme = qme + rrl
          if (lgrad1) then
            rrl2 = rrl*rrl
            dtrm = rrl*rrl2
            dqme(1) = dqme(1) + dtrm*xij
            dqme(2) = dqme(2) + dtrm*yji
            dqme(3) = dqme(3) + dtrm*zji
            if (lgrad2) then
              dtrm2 = 3.0_dp*rrl*rrl2*rrl2
              d2qme(1) = d2qme(1) - dtrm2*xij*xij
              d2qme(2) = d2qme(2) - dtrm2*xij*yji
              d2qme(3) = d2qme(3) - dtrm2*yji*yji
              d2qme(4) = d2qme(4) - dtrm2*xij*zji
              d2qme(5) = d2qme(5) - dtrm2*yji*zji
              d2qme(6) = d2qme(6) - dtrm2*zji*zji
              d2qme(1) = d2qme(1) + dtrm
              d2qme(3) = d2qme(3) + dtrm
              d2qme(6) = d2qme(6) + dtrm
            endif
          endif
        endif
      enddo
    endif
  endif
#ifdef TRACE
  call trace_out('qmatrixelemente')
#endif
!
  return
  end
