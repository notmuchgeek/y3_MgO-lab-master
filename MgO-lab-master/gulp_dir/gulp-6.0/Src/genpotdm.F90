  subroutine genpotdm(potl,maxd2,minuschi,imode)
!
!  Routine calculates the electrostatic potential at the
!  atoms and stores the result in a matrix form - only 
!  used with electronegativity equalisation. 
!
!  Distributed memory version. 
!
!   1/17 Created from genpot
!   7/17 Core-shell pair checks propagated to distance check
!   7/17 Global sum added for minuschi in S & M case
!   1/18 Trace added
!   5/18 Multiple qranges added
!   8/19 Trap for neemrptr being zero added
!
!  On entry:
!
!    imode = 1 => calculate A
!          = 2 => calculate A + (dA/dq).q
!    maxd2 = lower dimension of potl
!
!  On exit:
!
!    potl  = array A (+ (dA/dq).q) needed in electronegativity
!            equalisation
!    minuschi = modified in Streitz-Mintmire scheme
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
  use g_constants
  use control
  use current
  use eemdata
  use element
  use kspace
  use parallel
  use shells
  use symmetry
#ifdef TRACE
  use trace,       only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4),                 intent(in)    :: imode
  integer(i4),                 intent(in)    :: maxd2
  real(dp),                    intent(inout) :: potl(maxd2,*)
  real(dp),                    intent(inout) :: minuschi(*)
!
!  Local variables
!
  integer(i4)                                :: i
  integer(i4)                                :: ii
  integer(i4)                                :: iloc
  integer(i4)                                :: j
  integer(i4)                                :: jj
  integer(i4)                                :: kk
  integer(i4)                                :: maxloop2(3)
  integer(i4)                                :: ni
  integer(i4)                                :: nj
  integer(i4)                                :: npqni
  integer(i4)                                :: npqnj
  integer(i4)                                :: nqr
  integer(i4)                                :: nqrj
  integer(i4)                                :: status
  real(dp)                                   :: cuts2
  real(dp)                                   :: dgam
  real(dp)                                   :: dgamifj
  real(dp)                                   :: dgamjfi
  real(dp)                                   :: d2gamr2
  real(dp)                                   :: d2gamifj
  real(dp)                                   :: d2gamjfi
  real(dp)                                   :: dzetai
  real(dp)                                   :: dzetaj
  real(dp)                                   :: d2zetaii
  real(dp)                                   :: d2zetaij
  real(dp)                                   :: d2zetajj
  real(dp)                                   :: d2zetari
  real(dp)                                   :: d2zetarj
  real(dp)                                   :: gam
  real(dp)                                   :: gamifj
  real(dp)                                   :: gamjfi
  real(dp),    dimension(:), allocatable     :: minuschiloc
  real(dp)                                   :: qlii
  real(dp)                                   :: qljj
  real(dp)                                   :: rconv
  real(dp)                                   :: rconv2
  real(dp)                                   :: rl
  real(dp)                                   :: rqeq2
  real(dp)                                   :: rr2
  real(dp)                                   :: rtrmi1
  real(dp)                                   :: rtrmi3
  real(dp)                                   :: rtrmij
  real(dp)                                   :: rtrmj1
  real(dp)                                   :: rtrmj3
  real(dp)                                   :: rv2
  real(dp)                                   :: rx
  real(dp)                                   :: ry
  real(dp)                                   :: rz
  real(dp)                                   :: rxi
  real(dp)                                   :: ryi
  real(dp)                                   :: rzi
  real(dp)                                   :: rxj
  real(dp)                                   :: ryj
  real(dp)                                   :: rzj
  real(dp)                                   :: rxk
  real(dp)                                   :: ryk
  real(dp)                                   :: rzk
  real(dp),    dimension(:), allocatable     :: tmp
  real(dp)                                   :: trmi
  real(dp)                                   :: vij
  real(dp)                                   :: dvij(3)
  real(dp)                                   :: d2vij(6)
  real(dp)                                   :: xci
  real(dp)                                   :: yci
  real(dp)                                   :: zci
  real(dp)                                   :: xd
  real(dp)                                   :: yd
  real(dp)                                   :: zd
  real(dp)                                   :: zetai
  real(dp)                                   :: zetaj
  real(dp)                                   :: znuci
  real(dp)                                   :: znucj
  logical                                    :: lhi
  logical                                    :: lhj
#ifdef TRACE
  call trace_in('genpotdm')
#endif
!
!  Zero matrices
!
  do i = 1,natomsonnode
    do j = 1,numat
      potl(j,i) = 0.0_dp
    enddo
  enddo
!
!  If electrostatics have been turned off then there is no point in proceeding.
!
  if (.not.lDoElectrostatics) then
#ifdef TRACE
    call trace_out('genpotdm')
#endif
    return
  endif
!
!  Allocate workspace arrays
!
  if (lSandM) then
    allocate(minuschiloc(numat),stat=status)
    if (status/=0) call outofmemory('genpotdm','minuschiloc')
    allocate(tmp(numat),stat=status)
    if (status/=0) call outofmemory('genpotdm','tmp')
!
    minuschiloc(1:numat) = 0.0_dp
  endif
!
!  Local variables
!
  rconv = 1.0_dp/autoangs
  rconv2 = rconv*rconv
  cuts2 = cuts*cuts
  rqeq2 = rqeq*rqeq
  if (ndim.gt.0) then
!******************
!  Periodic case  *
!******************
!  
!  Initialise terms   
!     
    if (lewald.and.ndim.gt.1) then
      call kindex
      call initktrm
    endif
    call initqmatrix   
!***********************************
!  Calculate Coulomb contribution  *
!***********************************
    do iloc = 1,natomsonnode
      i = node2atom(iloc)
      xci = xclat(i)
      yci = yclat(i)
      zci = zclat(i)
      do j = 1,numat
        xd = xclat(j) - xci
        yd = yclat(j) - yci
        zd = zclat(j) - zci
        vij = 0.0_dp
!
!  If atoms are not a core-shell pair then pass cuts2 as 0
!
        if (j.eq.ncsptr(i)) then
          call qmatrixelement(xd,yd,zd,cuts2,.false.,.false.,vij,dvij,d2vij)
        else
          call qmatrixelement(xd,yd,zd,0.0_dp,.false.,.false.,vij,dvij,d2vij)
        endif
!
!  Evaluate potential due to reciprocal lattice summation.
!
        potl(j,iloc) = potl(j,iloc) + vij
      enddo
    enddo
    if (.not.lqeq.and..not.lSandM) goto 135
!**************************************************
!  QEq/SandM/ReaxFF corrections to Coulomb terms  *
!**************************************************
!
!  Estimate extent of summations over lattice vectors
!
    if (ndim.eq.3) then
      do i = 1,3
        rv2 = rv(1,i)**2 + rv(2,i)**2 + rv(3,i)**2
        rv2 = sqrt(rv2)
        maxloop2(i) = (rqeq/rv2) + 2
      enddo
    elseif (ndim.eq.2) then
      do i = 1,2
        rv2 = rv(1,i)**2 + rv(2,i)**2
        rv2 = sqrt(rv2)
        maxloop2(i) = (rqeq/rv2) + 2
      enddo
      maxloop2(3) = 0
    elseif (ndim.eq.1) then
      rv2 = rv(1,1)
      maxloop2(1) = (rqeq/rv2) + 1
      maxloop2(2) = 0
      maxloop2(3) = 0
    endif
!
!  Start loop over local atoms
!
    do iloc = 1,natomsonnode
      i = node2atom(iloc)
      ni = nat(i)
      qlii = qf(i)
      xci = xclat(i)
      yci = yclat(i)
      zci = zclat(i)
!
!  If QEq work out principal quantum number
!
      if (lqeq) then
        if (lmultiqrange.and.neemrptr(i).ne.0) then
          nqr = nqrnow(neemrptr(i))
        else
          nqr = 1
        endif
        if (ni.le.2) then
          npqni = 1
        elseif (ni.le.10) then
          npqni = 2
        elseif (ni.le.18) then
          npqni = 3
        elseif (ni.le.36) then
          npqni = 4
        elseif (ni.le.54) then
          npqni = 5
        elseif (ni.le.86) then
          npqni = 6
        else
          npqni = 7
        endif
        zetai = 0.5_dp*qeqlambda*(2*npqni+1)/radrange(nqr,ni,neemtype)
        lhi = (ni.eq.1)
        if (lhi) then
!
!  Special case for hydrogen
!
          zetai = zetai + qlii*rconv
        endif
      elseif (lSandM) then
        if (lmultiqrange.and.neemrptr(i).ne.0) then
          nqr = nqrnow(neemrptr(i))
        else
          nqr = 1
        endif
        zetai = zetarange(nqr,ni,neemtype)
        znuci = znucrange(nqr,ni,neemtype)
      endif
      do j = 1,numat
        qljj = qf(j)
        rx = xclat(j) - xci
        ry = yclat(j) - yci
        rz = zclat(j) - zci
        nj = nat(j)
!
!  If QEq work out principal quantum number
!
        if (lqeq) then
          if (lmultiqrange.and.neemrptr(j).ne.0) then
            nqrj = nqrnow(neemrptr(j))
          else
            nqrj = 1
          endif
          if (nj.le.2) then
            npqnj = 1
          elseif (nj.le.10) then
            npqnj = 2
          elseif (nj.le.18) then
            npqnj = 3
          elseif (nj.le.36) then
            npqnj = 4
          elseif (nj.le.54) then
            npqnj = 5
          elseif (nj.le.86) then
            npqnj = 6
          else
            npqnj = 7
          endif
          zetaj = 0.5_dp*qeqlambda*(2*npqnj+1)/radrange(nqrj,nj,neemtype)
          lhj = (nj.eq.1)
          if (lhj) then
!
!  Special case for hydrogen
!
            zetaj = zetaj + qljj*rconv
          endif
        elseif (lSandM) then
          if (lmultiqrange.and.neemrptr(j).ne.0) then
            nqrj = nqrnow(neemrptr(j))
          else
            nqrj = 1
          endif
          zetaj = zetarange(nqrj,nj,neemtype)
          znucj = znucrange(nqrj,nj,neemtype)
        endif
!
!  Loop over cell vectors
!
        rxi = rx - (maxloop2(1) + 1)*r1x
        ryi = ry - (maxloop2(1) + 1)*r1y
        rzi = rz - (maxloop2(1) + 1)*r1z
        do ii = - maxloop2(1),maxloop2(1)
          rxi = rxi + r1x
          ryi = ryi + r1y
          rzi = rzi + r1z
          rxj = rxi - (maxloop2(2) + 1)*r2x
          ryj = ryi - (maxloop2(2) + 1)*r2y
          rzj = rzi - (maxloop2(2) + 1)*r2z
          do jj = - maxloop2(2),maxloop2(2)
            rxj = rxj + r2x
            ryj = ryj + r2y
            rzj = rzj + r2z
            rxk = rxj - (maxloop2(3) + 1)*r3x
            ryk = ryj - (maxloop2(3) + 1)*r3y
            rzk = rzj - (maxloop2(3) + 1)*r3z
            do 120 kk = - maxloop2(3),maxloop2(3)
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
              if (rr2.gt.rqeq2) goto 120
!
!  Trap self term
!
              if (rr2.lt.cuts2.and.(j.eq.ncsptr(i))) goto 120
!
!  Trap zero term
!
              if (i.eq.j.and.rr2.lt.1.0d-6) goto 120
              rl = sqrt(rr2)
              if (lqeq) then
!
!  Calculate Coulomb interaction according to QEq scheme and
!  subtract 1/r term from Ewald sum.
!
                call gammas(npqni,npqnj,zetai,zetaj,rl,gam,dgam,dzetai,dzetaj,d2zetaii,d2zetaij, &
                            d2zetajj,d2zetari,d2zetarj,d2gamr2)
                potl(j,iloc) = potl(j,iloc) - 1.0_dp/rl + gam
                if (index(keyword,'oldqeq').eq.0) then
                  if (lhi) then
                    potl(j,iloc) = potl(j,iloc) + qlii*dzetai*rconv
                  endif
                  if (imode.eq.2) then
!
!  Add (dA/dq).q term to A
!
                    if (lhi) then
                      rtrmi1 = dzetai*rconv
                      rtrmi3 = qlii*d2zetaii*rconv2
                      potl(i,iloc) = potl(i,iloc) + qljj*(2.0_dp*rtrmi1+rtrmi3)
                    endif
                    if (lhj.and.i.ne.j) then
                      rtrmj1 = dzetaj*rconv
                      rtrmj3 = qljj*d2zetajj*rconv2
                      potl(j,iloc) = potl(j,iloc) + qljj*rtrmj1
                    endif
!
!  Term if both atoms are hydrogen
!
                    if (lhi.and.lhj) then
                      rtrmij = qlii*qljj*d2zetaij*rconv2
                      potl(j,iloc) = potl(j,iloc) + rtrmij
                    endif
                  endif
                endif
              elseif (lSandM) then
                call gammasm(zetai,zetaj,rl,gam,dgam,d2gamr2,gamifj,gamjfi,dgamifj,dgamjfi,d2gamifj,d2gamjfi)
                potl(j,iloc) = potl(j,iloc) + gam
                minuschiloc(i) = minuschiloc(i) - angstoev*znucj*(gamjfi - gam)
              endif
!
!  End of loops over lattice vectors
!
120           continue
          enddo
        enddo
!
!  End of loops over atom pairs
!
      enddo
    enddo
  else
!*****************
!  Cluster case  *
!*****************
!
!  Start loop over cluster atom - unit cell atom pairs
!
    do iloc = 1,natomsonnode
      i = node2atom(iloc)
      qlii = qf(i)
      xci = xclat(i)
      yci = yclat(i)
      zci = zclat(i)
      ni = nat(i)
!
!  If QEq work out principal quantum number
!
      if (lqeq) then
        if (lmultiqrange.and.neemrptr(i).ne.0) then
          nqr = nqrnow(neemrptr(i))
        else
          nqr = 1
        endif
        if (ni.le.2) then
          npqni = 1
        elseif (ni.le.10) then
          npqni = 2
        elseif (ni.le.18) then
          npqni = 3
        elseif (ni.le.36) then
          npqni = 4
        elseif (ni.le.54) then
          npqni = 5
        elseif (ni.le.86) then
          npqni = 6
        else
          npqni = 7
        endif
        zetai = 0.5_dp*qeqlambda*(2*npqni+1)/radrange(nqr,ni,neemtype)
        lhi = (ni.eq.1)
        if (lhi) then
!
!  Special case for hydrogen
!
          zetai = zetai + qlii*rconv
        endif
      elseif (lSandM) then
        if (lmultiqrange.and.neemrptr(i).ne.0) then
          nqr = nqrnow(neemrptr(i))
        else
          nqr = 1
        endif
        zetai = zetarange(nqr,ni,neemtype)
        znuci = znucrange(nqr,ni,neemtype)
      endif
      jloop: do j = 1,numat
!
!  Trap zero term
!
        if (i.eq.j) cycle jloop
!
        qljj = qf(j)
        nj = nat(j)
!
!  If QEq work out principal quantum number
!
        if (lqeq) then
          if (lmultiqrange.and.neemrptr(j).ne.0) then
            nqrj = nqrnow(neemrptr(j))
          else
            nqrj = 1
          endif
          if (nj.le.2) then
            npqnj = 1
          elseif (nj.le.10) then
            npqnj = 2
          elseif (nj.le.18) then
            npqnj = 3
          elseif (nj.le.36) then
            npqnj = 4
          elseif (nj.le.54) then
            npqnj = 5
          elseif (nj.le.86) then
            npqnj = 6
          else
            npqnj = 7
          endif
          zetaj = 0.5_dp*qeqlambda*(2*npqnj+1)/radrange(nqrj,nj,neemtype)
          lhj = (nj.eq.1)
          if (lhj) then
!
!  Special case for hydrogen
!
            zetaj = zetaj + qljj*rconv
          endif
        elseif (lSandM) then
          if (lmultiqrange.and.neemrptr(j).ne.0) then
            nqrj = nqrnow(neemrptr(j))
          else
            nqrj = 1
          endif
          zetaj = zetarange(nqrj,nj,neemtype)
          znucj = znucrange(nqrj,nj,neemtype)
        endif
!
!  Find relative vector between atoms
!
        rx = xclat(j) - xci
        ry = yclat(j) - yci
        rz = zclat(j) - zci
        rr2 = rx*rx + ry*ry + rz*rz
!
!  Exclude core-shell interaction
!
        if (rr2.lt.cuts2.and.(j.eq.ncsptr(i))) cycle jloop
        rl = sqrt(rr2)
        if (lqeq.and.rl.lt.rqeq) then
!***************
!  QEq scheme  *
!***************
          call gammas(npqni,npqnj,zetai,zetaj,rl,gam,dgam,dzetai,dzetaj, &
                      d2zetaii,d2zetaij,d2zetajj,d2zetari,d2zetarj,d2gamr2)
          potl(j,iloc) = potl(j,iloc) + gam
          if (index(keyword,'oldqeq').eq.0) then
            if (lhi) then
              potl(j,iloc) = potl(j,iloc) + qlii*dzetai*rconv
            endif
            if (imode.eq.2) then
!
!  Add (dA/dq).q term to A
!
              if (lhi) then
                rtrmi1 = dzetai*rconv
                rtrmi3 = qlii*d2zetaii*rconv2
                potl(i,iloc) = potl(i,iloc) + qljj*(2.0_dp*rtrmi1 + rtrmi3)
              endif
              if (lhj) then
                rtrmj1 = dzetaj*rconv
                rtrmj3 = qljj*d2zetajj*rconv2
                potl(j,iloc) = potl(j,iloc) + qljj*rtrmj1
              endif
!
!  Term if both atoms are hydrogen
!
              if (lhi.and.lhj) then
                rtrmij = qlii*qljj*d2zetaij*rconv2
                potl(j,iloc) = potl(j,iloc) + rtrmij
              endif
!
            endif
          endif
        elseif (lSandM.and.rl.lt.rqeq) then
!*******************
!  S and M scheme  *
!*******************
          call gammasm(zetai,zetaj,rl,gam,dgam,d2gamr2,gamifj,gamjfi,dgamifj,dgamjfi,d2gamifj,d2gamjfi)
          potl(j,iloc) = potl(j,iloc) + gam + 1.0_dp/rl
          minuschiloc(i) = minuschiloc(i) - angstoev*znucj*(gamjfi - gam)
        else
!***************
!  EEM scheme  *
!***************
          trmi = 1.0_dp/rl
          potl(j,iloc) = trmi + potl(j,iloc)
        endif
      enddo jloop
    enddo
  endif
!
!  Global sum of minuschi
!
  if (lSandM) then
    call sumall(minuschiloc,tmp,numat,"genpotdm","minuschi")
    do i = 1,numat
      minuschi(i) = minuschi(i) + tmp(i)
    enddo
  endif
!
135 continue
!
!  Convert units to eV
!
  do i = 1,natomsonnode
    do j = 1,numat
      potl(j,i) = potl(j,i)*angstoev
    enddo
  enddo
!
!  Free workspace memory
!
  if (lSandM) then
    deallocate(tmp,stat=status)
    if (status/=0) call deallocate_error('genpotdm','tmp')
    deallocate(minuschiloc,stat=status)
    if (status/=0) call deallocate_error('genpotdm','minuschiloc')
  endif
#ifdef TRACE
  call trace_out('genpotdm')
#endif
!
  return
  end
