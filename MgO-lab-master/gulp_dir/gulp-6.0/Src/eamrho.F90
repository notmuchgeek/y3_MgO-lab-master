  subroutine eamrho(nati,ntypi,natj,ntypj,r,rmax,x,y,z,rhoij,rhoji,drhoij,drhoji, &
                    drhoijs,drhojis,drhoij2,drhoji2,drhoij2s,drhoji2s,drhoij2m,drhoji2m, &
                    drhoij3,drhoji3,oci,ocj,lstr,lgrad1,lgrad2,lgrad3,rtaper)
!
!  Calculates the density function and its derivatives in the Embedded Atom Method. 
!  Assumes that derivatives have been pre-initialised by the calling routine.
!
!  On entry :
!
!  nati     = atomic number of i
!  ntypi    = type number of i
!  natj     = atomic number of j
!  ntypj    = type number of j
!  r        = distance from i to j (Angstroms)
!  rmax     = cutoff distance from i to j (Angstroms)
!  x        = Cartesian component of r in x direction (Angstroms)
!  y        = Cartesian component of r in y direction (Angstroms)
!  z        = Cartesian component of r in z direction (Angstroms)
!  oci      = occupancy of site i
!  ocj      = occupancy of site j
!  rtaper   = range of taper function
!  lgrad2   = if .true. calculate first derivatives 
!  lgrad2   = if .true. and lgrad1, calculate second derivatives 
!  lgrad3   = if .true. and lgrad2, calculate third derivatives
!  lstr     = if .true. calculate strain derivatives 
!
!  On exit :
!
!  rhoij    = density of j at i (if MEAM then this is the contributions for each order)
!  rhoji    = density of i at j (if MEAM then this is the contributions for each order)
!  drhoij   = 1st derivative of j density at i w.r.t. Cartesian directions (if lgrad1 = .true.)
!  drhoji   = 1st derivative of i density at j w.r.t. Cartesian directions (if lgrad1 = .true.)
!  drhoijs  = 1st derivative of j density at i w.r.t. strain (if lgrad1 = .true.)
!  drhojis  = 1st derivative of i density at j w.r.t. strain (if lgrad1 = .true.)
!  drhoij2  = 2nd derivative of j density at i w.r.t. Cartesian directions (if lgrad2 = .true.)
!  drhoji2  = 2nd derivative of i density at j w.r.t. Cartesian directions (if lgrad2 = .true.)
!  drhoij2s = 2nd derivative of j density at i w.r.t. strain (if lgrad2 & lstr = .true.)
!  drhoji2s = 2nd derivative of i density at j w.r.t. strain (if lgrad2 & lstr = .true.)
!  drhoij2m = 2nd derivative of j density at i w.r.t. mixed Cartesian direction/strain (if lgrad2 & lstr = .true.)
!  drhoji2m = 2nd derivative of i density at j w.r.t. mixed Cartesian direction/strain (if lgrad2 & lstr = .true.)
!  drhoij3  = 3rd derivative of j density at i w.r.t. Cartesian directions (if lgrad3 = .true.)
!  drhoji3  = 3rd derivative of i density at j w.r.t. Cartesian directions (if lgrad3 = .true.)
!
!  10/99 Cubic density function added
!  11/03 ndennat/ndentyp replaced
!   7/05 Style updated
!   9/05 Voter form of density added
!  11/05 Tapering of density added
!   2/06 Quadratic and quartic densities added
!   3/06 Power law EAM densities truncated after r0
!   3/06 Modified to allow for density component number
!   4/06 Species specific density added
!   3/07 Glue density added
!   5/07 eVoter EAM density added
!  11/07 Mei-Davenport density added
!  11/07 MDF tapering of density added
!  10/08 MEAM modifications added - denpar increased in dimension
!  11/08 Calculation of density added to this routine so that all 
!        terms are collected in a single place
!  11/08 Baskes form of exponential density added
!  11/08 Modified so that rho arguments are arrays rather than scalars to 
!        handle the needs of MEAM
!  11/08 rho arrays made 2-D to accommodate MEAM
!  11/08 rho*sum variables removed since terms for each order have to be kept separate
!  11/08 Created from rhoderv
!   1/09 VBO density added
!   1/09 x, y, z now passed in for consistency with meamrho
!   1/09 Density derivative arguments are now arrays to accommodate MEAM
!   2/09 Strain derivatives added
!   7/09 Use of drhojiloc/drhoijloc when not defined corrected
!  10/11 Fractional power density added
!  10/12 Intent of rhoij/rhoji changed to inout
!   7/13 Cubic spline EAM density added
!   8/14 Taper range passed as argument
!   8/14 Specific taper added for EAM density
!  12/14 eamalloy parameters added
!   2/18 Trace added
!   9/18 Modified for changes due to lstraincell algorithm
!   9/18 Strain module introduced
!  12/19 Rigid molecule modifications added
!   3/20 Morse squared added
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
!  Julian Gale, CIC, Curtin University, March 2020
!
  use current,        only : ndim, nstrptr, nstrains
  use eam
  use m_strain,       only : real1strterm
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)    :: nati
  integer(i4), intent(in)    :: natj
  integer(i4), intent(in)    :: ntypi
  integer(i4), intent(in)    :: ntypj
  logical,     intent(in)    :: lgrad1
  logical,     intent(in)    :: lgrad2
  logical,     intent(in)    :: lgrad3
  logical,     intent(in)    :: lstr
  real(dp),    intent(inout) :: rhoij
  real(dp),    intent(inout) :: rhoji
  real(dp),    intent(inout) :: drhoij(3)
  real(dp),    intent(inout) :: drhoijs(6)
  real(dp),    intent(inout) :: drhoij2(6)
  real(dp),    intent(inout) :: drhoij2s(21)
  real(dp),    intent(inout) :: drhoij2m(6,3)
  real(dp),    intent(inout) :: drhoij3(10)
  real(dp),    intent(inout) :: drhoji(3)
  real(dp),    intent(inout) :: drhojis(6)
  real(dp),    intent(inout) :: drhoji2(6)
  real(dp),    intent(inout) :: drhoji2s(21)
  real(dp),    intent(inout) :: drhoji2m(6,3)
  real(dp),    intent(inout) :: drhoji3(10)
  real(dp),    intent(in)    :: oci
  real(dp),    intent(in)    :: ocj
  real(dp),    intent(in)    :: r
  real(dp),    intent(in)    :: rmax
  real(dp),    intent(in)    :: rtaper
  real(dp),    intent(in)    :: x
  real(dp),    intent(in)    :: y
  real(dp),    intent(in)    :: z
!
!  Local variables
!
  integer(i4)                :: ind
  integer(i4)                :: is1
  integer(i4)                :: is2
  integer(i4)                :: j
  integer(i4)                :: l
  integer(i4)                :: mm
  integer(i4)                :: npt
  integer(i4)                :: ns1
  integer(i4)                :: ns2
  real(dp)                   :: apt
  real(dp)                   :: alloy
  real(dp)                   :: bpt
  real(dp)                   :: cpt
  real(dp)                   :: dpt
  real(dp)                   :: ept
  real(dp)                   :: dr
  real(dp)                   :: dr2
  real(dp)                   :: dr3
  real(dp)                   :: etrm
  real(dp)                   :: exptrm
  real(dp)                   :: r2
  real(dp)                   :: rd
  real(dp)                   :: rd2
  real(dp)                   :: dbdx
  real(dp)                   :: d2bdx2
  real(dp)                   :: d3bdx3
  real(dp)                   :: drhoijloc
  real(dp)                   :: drhoij2loc
  real(dp)                   :: drhoij3loc
  real(dp)                   :: drhojiloc
  real(dp)                   :: drhoji2loc
  real(dp)                   :: drhoji3loc
  real(dp)                   :: drhoijlocsum
  real(dp)                   :: drhoij2locsum
  real(dp)                   :: drhoij3locsum
  real(dp)                   :: drhojilocsum
  real(dp)                   :: drhoji2locsum
  real(dp)                   :: drhoji3locsum
  real(dp)                   :: detpfn
  real(dp)                   :: d2etpfn
  real(dp)                   :: d3etpfn
  real(dp)                   :: dxdr
  real(dp)                   :: d2xdr2
  real(dp)                   :: d3xdr3
  real(dp)                   :: dr2ds(6)
  real(dp)                   :: d2r2dx2(3,3)
  real(dp)                   :: d2r2ds2(6,6)
  real(dp)                   :: d2r2dsdx(6,3)
  real(dp)                   :: etpfn
  real(dp)                   :: r0
  real(dp)                   :: rhoijloc
  real(dp)                   :: rhojiloc
  real(dp)                   :: rhoijlocsum
  real(dp)                   :: rhojilocsum
  real(dp)                   :: rk
  real(dp)                   :: rk2
  real(dp)                   :: rmin
  real(dp)                   :: rn1
  real(dp)                   :: rpt
  real(dp)                   :: rr0
  real(dp)                   :: rr12
  real(dp)                   :: rx
  real(dp)                   :: rxm1
  real(dp)                   :: term0
  real(dp)                   :: term1
  real(dp)                   :: term2
  real(dp)                   :: term3
  real(dp)                   :: tpfn
  real(dp)                   :: dtpfn
  real(dp)                   :: d2tpfn
  real(dp)                   :: d3tpfn
  real(dp)                   :: trm1
  real(dp)                   :: trm2
  real(dp)                   :: trm3
  real(dp)                   :: trm4
  real(dp)                   :: dtrm2
  real(dp)                   :: dtrm3
  real(dp)                   :: d2trm2
  real(dp)                   :: d2trm3
  real(dp)                   :: d3trm2
  real(dp)                   :: d3trm3
#ifdef TRACE
  call trace_in('eamrho')
#endif
!
!  Calculate local variables
!
  r2 = r*r
  rk = 1.0_dp/r
  rk2 = rk*rk
  if (lstr) then
    call real1strterm(ndim,x,y,z,0.0_dp,0.0_dp,0.0_dp,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,lgrad2)
  endif
!
!  Set up general taper values
!
  if (rtaper.gt.1.0d-12) then
!
!  EAM taper
!
    call eamtaper(r,rmax-rtaper,rmax,tpfn,dtpfn,d2tpfn,d3tpfn,lgrad1,lgrad2,lgrad3)
!
!  Convert taper terms to allow for 1/r factors
!
    if (lgrad1) then
      dtpfn = rk*dtpfn
      if (lgrad2) then
        d2tpfn = rk2*d2tpfn
        if (lgrad3) then
          d3tpfn = rk2*(rk*d3tpfn - 3.0_dp*(d2tpfn - rk2*dtpfn))
        endif
        d2tpfn = d2tpfn - rk2*dtpfn
      endif
    endif
  endif
!
  do mm = 1,neamspec
!
!  Terms for i->j
!
    if (nati.eq.neamnat(mm).and.(ntypi.eq.neamtyp(mm).or.neamtyp(mm).eq.0)) then
      if (neamnat2(mm).eq.0.or.(neamnat2(mm).eq.natj.and.(ntypj.eq.neamtyp2(mm).or.neamtyp2(mm).eq.0))) then
        rhojilocsum = 0.0_dp
        drhojilocsum = 0.0_dp
        drhoji2locsum = 0.0_dp
        drhoji3locsum = 0.0_dp
        alloy = eamalloy(1,mm)
        do j = 1,ndenfncomp(mm)
!
!  Density functional forms
!
          npt = nint(denpar(6,1,j,mm))
          if (ndenfn(j,mm).eq.1) then
!
!  Power law
!
            trm1 = oci*denpar(1,1,j,mm)*(rk**npt)
            rhojiloc = trm1
            if (lgrad1) then
              trm1 = - dble(npt)*trm1*rk2
              drhojiloc = trm1
              if (lgrad2) then
                trm2 = - trm1*dble(npt+2)*rk2
                drhoji2loc = trm2
                if (lgrad3) then
                  drhoji3loc = - trm2*dble(npt+4)*rk2
                endif
              endif
            endif
          elseif (ndenfn(j,mm).eq.2) then
!
!  Exponential
!
            apt = denpar(1,1,j,mm)
            bpt = denpar(2,1,j,mm)
            cpt = denpar(3,1,j,mm)
            trm1 = apt*oci*exp(-bpt*(r-cpt))
            if (npt.ne.0) then
              rn1 = r**(npt-1)
              rhojiloc = trm1*rn1*r
              if (lgrad1) then
                drhojiloc = trm1*rn1*(npt*rk-bpt)
                if (lgrad2) then
                  drhoji2loc = trm1*rn1*rk*(npt*(npt-1)*rk2 - bpt*dble(2*npt-1)*rk + bpt*bpt)
                  if (lgrad3) then
                    drhoji3loc = trm1*rn1*rk2*(dble(npt*(npt-2)*(npt-4))*rk2*rk &
                      - 3.0_dp*bpt*dble(npt*(npt-3)+1)*rk2 + 3.0_dp*bpt*bpt*dble(npt-1)*rk - bpt*bpt*bpt)
                  endif
                endif
              endif
            else
              rhojiloc = trm1
              if (lgrad1) then
                drhojiloc = - bpt*trm1*rk
                if (lgrad2) then
                  drhoji2loc = bpt*trm1*rk2*(bpt + rk)
                  if (lgrad3) then
                    drhoji3loc = - bpt*trm1*rk2*rk*(bpt*(bpt+rk) - rk2)
                  endif
                endif
              endif
            endif
          elseif (ndenfn(j,mm).eq.3) then
!
!  Gaussian
!
            apt = denpar(1,1,j,mm)
            bpt = denpar(2,1,j,mm)
            cpt = denpar(3,1,j,mm)
            rd = (r - cpt)
            rd2 = rd*rd
            trm1 = apt*oci*exp(-bpt*rd2)
            if (npt.ne.0) then
              rn1 = r**(npt-1)
              rhojiloc = trm1*rn1*r
              if (lgrad1) then
                drhojiloc = trm1*rn1*(npt*rk - 2.0_dp*bpt*rd)
                if (lgrad2) then
                  drhoji2loc = trm1*rn1*rk*(npt*(npt-2)*rk2- &
                    2.0_dp*bpt*rd*dble(2*npt-1)*rk+2.0_dp*bpt*(2.0_dp*bpt*rd2-1._dp))
                  if (lgrad3) then
                    drhoji3loc = trm1*rn1*rk2*(dble(npt*(npt-2)*(npt-4))*rk2*rk-6.0_dp*bpt*rd*dble(npt &
                      *npt-3*npt+1)*rk2-6.0_dp*bpt*dble(npt-1)*rk*(1.0_dp-2.0_dp*bpt*rd2)+4.0_dp*bpt*bpt*rd*( &
                      3.0_dp-2.0_dp*bpt*rd2))
                  endif
                endif
              endif
            else
              rhojiloc = trm1
              if (lgrad1) then
                drhojiloc = - 2.0_dp*bpt*rd*trm1*rk
                if (lgrad2) then
                  drhoji2loc = 2.0_dp*trm1*rk2*bpt*(rd*rk-1.0_dp+4.0_dp*bpt*rd2)
                  if (lgrad3) then
                    drhoji3loc = - 4.0_dp*trm1*bpt*rk2*(bpt*rd*rk+rk2)* &
                      (rd*rk-1.0_dp+4.0_dp*bpt*rd2)+2.0_dp*trm1*rk2*bpt*( &
                      -rd*rk2+8.0_dp*bpt*rd)*rk
                  endif
                endif
              endif
            endif
          elseif (ndenfn(j,mm).eq.4) then
!
!  Cubic
!
            apt = denpar(1,1,j,mm)
            bpt = denpar(2,1,j,mm)
            if (r.lt.bpt) then
              rd = (r - bpt)
              trm1 = 3.0_dp*apt*rk
              rhojiloc = apt*rd*rd*rd
              if (lgrad1) then
                drhojiloc = trm1*rd*rd
                if (lgrad2) then
                  trm1 = trm1*rk
                  drhoji2loc = trm1*rd*(2.0_dp - rd*rk)
                  if (lgrad3) then
                    trm1 = trm1*rk
                    drhoji3loc = trm1*(2.0_dp - 6.0_dp*rd*rk + 3.0_dp*rd*rd*rk*rk)
                  endif
                endif
              endif
            else
              rhojiloc = 0.0_dp
              if (lgrad1) then
                drhojiloc = 0.0_dp
                if (lgrad2) then
                  drhoji2loc = 0.0_dp
                  if (lgrad3) then
                    drhoji3loc = 0.0_dp
                  endif
                endif
              endif
            endif
          elseif (ndenfn(j,mm).eq.5) then
!
!  Voter-Chen
!
            apt = denpar(1,1,j,mm)
            bpt = denpar(2,1,j,mm)
            etrm = exp(-bpt*r)
            cpt = 2.0_dp**9
            trm1 = apt*etrm*(1.0_dp + cpt*etrm)
            trm2 = apt*etrm*(1.0_dp + 2.0_dp*cpt*etrm)
            rhojiloc = r2*r2*r2*trm1
            if (lgrad1) then
              drhojiloc = r2*r2*(6.0_dp*trm1 - bpt*r*trm2)
              if (lgrad2) then
                trm3 = apt*etrm*(1.0_dp + 4.0_dp*cpt*etrm)
                drhoji2loc = r2*(24.0_dp*trm1 - 11.0_dp*bpt*r*trm2 + bpt*bpt*r2*trm3)
                if (lgrad3) then
                  trm4 = apt*etrm*(1.0_dp + 8.0_dp*cpt*etrm)
                  drhoji3loc = 48.0_dp*trm1 - bpt*r*(57.0_dp*trm2 - bpt*r*(15.0_dp*trm3 - bpt*r*trm4))
                endif
              endif
            endif
          elseif (ndenfn(j,mm).eq.6) then
!
!  Quadratic
!
            apt = denpar(1,1,j,mm)
            bpt = denpar(2,1,j,mm)
            if (r.lt.bpt) then
              rd = (r - bpt)
              trm1 = 2.0_dp*apt*rk
              rhojiloc = apt*rd*rd
              if (lgrad1) then
                drhojiloc = trm1*rd
                if (lgrad2) then
                  trm1 = trm1*rk
                  drhoji2loc = trm1*(1.0_dp - rd*rk)
                  if (lgrad3) then
                    trm1 = 3.0_dp*trm1*rk*rk
                    drhoji3loc = - trm1*(1.0_dp - rd*rk)
                  endif
                endif
              endif
            else
              rhojiloc = 0.0_dp
              if (lgrad1) then
                drhojiloc = 0.0_dp
                if (lgrad2) then
                  drhoji2loc = 0.0_dp
                  if (lgrad3) then
                    drhoji3loc = 0.0_dp
                  endif
                endif
              endif
            endif
          elseif (ndenfn(j,mm).eq.7) then
!
!  Quartic
!
            apt = denpar(1,1,j,mm)
            bpt = denpar(2,1,j,mm)
            if (r.lt.bpt) then
              rd = (r - bpt)
              trm1 = 4.0_dp*apt*rk*rd
              rhojiloc = apt*rd*rd*rd*rd
              if (lgrad1) then
                drhojiloc = trm1*rd*rd
                if (lgrad2) then
                  trm1 = trm1*rk
                  drhoji2loc = trm1*rd*(3.0_dp - rd*rk)
                  if (lgrad3) then
                    trm1 = trm1*rk
                    drhoji3loc = 3.0_dp*trm1*(2.0_dp - 3.0_dp*rd*rk + rd*rd*rk*rk)
                  endif
                endif
              endif
            else
              rhojiloc = 0.0_dp
              if (lgrad1) then
                drhojiloc = 0.0_dp
                if (lgrad2) then
                  drhoji2loc = 0.0_dp
                  if (lgrad3) then
                    drhoji3loc = 0.0_dp
                  endif
                endif
              endif
            endif
          elseif (ndenfn(j,mm).eq.8) then
!
!  Glue
!
            if (r.lt.denpar(1,1,j,mm)) then
              dr = (r - denpar(1,1,j,mm))
              dr2 = dr*dr
              dr3 = dr2*dr
              trm1 = 3.0_dp*denpar(4,1,j,mm)*dr2 + 2.0_dp*denpar(5,1,j,mm)*dr + denpar(6,1,j,mm)
              rhojiloc = denpar(4,1,j,mm)*dr3 + denpar(5,1,j,mm)*dr2 + denpar(6,1,j,mm)*dr + denpar(7,1,j,mm)
              if (lgrad1) then
                drhojiloc = trm1*rk
                if (lgrad2) then
                  trm2 = 6.0_dp*denpar(4,1,j,mm)*dr + 2.0_dp*denpar(5,1,j,mm)
                  drhoji2loc = rk2*(trm2 - rk*trm1)
                  if (lgrad3) then
                    trm3 = 6.0_dp*denpar(4,1,j,mm)
                    drhoji3loc = rk2*rk*(trm3 - 3.0_dp*rk*trm2 + 3.0_dp*rk2*trm1)
                  endif
                endif
              endif
            elseif (r.lt.denpar(2,1,j,mm)) then
              dr = (r - denpar(1,1,j,mm))
              dr2 = dr*dr
              dr3 = dr2*dr
              trm1 = 3.0_dp*denpar(8,1,j,mm)*dr2 + 2.0_dp*denpar(9,1,j,mm)*dr + denpar(10,1,j,mm)
              rhojiloc = denpar(8,1,j,mm)*dr3 + denpar(9,1,j,mm)*dr2 + denpar(10,1,j,mm)*dr + denpar(11,1,j,mm)
              if (lgrad1) then
                drhojiloc = trm1*rk
                if (lgrad2) then
                  trm2 = 6.0_dp*denpar(8,1,j,mm)*dr + 2.0_dp*denpar(9,1,j,mm)
                  drhoji2loc = rk2*(trm2 - rk*trm1)
                  if (lgrad3) then
                    trm3 = 6.0_dp*denpar(8,1,j,mm)
                    drhoji3loc = rk2*rk*(trm3 - 3.0_dp*rk*trm2 + 3.0_dp*rk2*trm1)
                  endif
                endif
              endif
            elseif (r.lt.denpar(3,1,j,mm)) then
              dr = (r - denpar(3,1,j,mm))
              dr2 = dr*dr
              dr3 = dr2*dr
              trm1 = 3.0_dp*denpar(12,1,j,mm)*dr2 + 2.0_dp*denpar(13,1,j,mm)*dr + denpar(14,1,j,mm)
              rhojiloc = denpar(12,1,j,mm)*dr3 + denpar(13,1,j,mm)*dr2 + denpar(14,1,j,mm)*dr + denpar(15,1,j,mm)
              if (lgrad1) then
                drhojiloc = trm1*rk
                if (lgrad2) then
                  trm2 = 6.0_dp*denpar(12,1,j,mm)*dr + 2.0_dp*denpar(13,1,j,mm)
                  drhoji2loc = rk2*(trm2 - rk*trm1)
                  if (lgrad3) then
                    trm3 = 6.0_dp*denpar(12,1,j,mm)
                    drhoji3loc = rk2*rk*(trm3 - 3.0_dp*rk*trm2 + 3.0_dp*rk2*trm1)
                  endif
                endif
              endif
            else
              rhojiloc = 0.0_dp
              if (lgrad1) then
                drhojiloc = 0.0_dp
                if (lgrad2) then
                  drhoji2loc = 0.0_dp
                  if (lgrad3) then
                    drhoji3loc = 0.0_dp
                  endif
                endif
              endif
            endif
          elseif (ndenfn(j,mm).eq.9) then
!
!  eVoter-Chen
!    
            call etaper(r,0.0_dp,rmax,etpfn,detpfn,d2etpfn,d3etpfn,lgrad1,lgrad2,lgrad3)
            detpfn = rk*detpfn
            d2etpfn = rk2*d2etpfn
            d3etpfn = rk2*(rk*d3etpfn - 3.0_dp*(d2etpfn - rk2*detpfn))
            d2etpfn = d2etpfn - rk2*detpfn
!
            apt = denpar(1,1,j,mm)
            bpt = denpar(2,1,j,mm)
            etrm = exp(-bpt*r)
            cpt = 2.0_dp**9
            trm1 = apt*etrm*(1.0_dp + cpt*etrm)
            trm2 = apt*etrm*(1.0_dp + 2.0_dp*cpt*etrm)
            rhojiloc = r2*r2*r2*trm1
            if (lgrad1) then
              drhojiloc = r2*r2*(6.0_dp*trm1 - bpt*r*trm2)
              if (lgrad2) then
                trm3 = apt*etrm*(1.0_dp + 4.0_dp*cpt*etrm)
                drhoji2loc = r2*(24.0_dp*trm1 - 11.0_dp*bpt*r*trm2 + bpt*bpt*r2*trm3)
                if (lgrad3) then
                  trm4 = apt*etrm*(1.0_dp + 8.0_dp*cpt*etrm)
                  drhoji3loc = 48.0_dp*trm1 - bpt*r*(57.0_dp*trm2 - bpt*r*(15.0_dp*trm3 - bpt*r*trm4))
                  drhoji3loc = drhoji3loc*etpfn + 3.0_dp*(drhoji2loc*detpfn + drhojiloc*d2etpfn) + rhojiloc*d3etpfn
                endif
                drhoji2loc = drhoji2loc*etpfn + 2.0_dp*drhojiloc*detpfn + rhojiloc*d2etpfn
              endif
              drhojiloc = drhojiloc*etpfn + rhojiloc*detpfn
            endif
          elseif (ndenfn(j,mm).eq.10) then
!
!  Mei-Davenport
!   
            rr0 = denpar(7,1,j,mm)
            rr12 = 1.0_dp/12.0_dp
            trm1 = rk*rr0
            rr0 = 1.0_dp/denpar(7,1,j,mm)
            rhojiloc = 0.0_dp
            if (lgrad1) then
              drhojiloc = 0.0_dp
              if (lgrad2) then
                drhoji2loc = 0.0_dp
                if (lgrad3) then
                  drhoji3loc = 0.0_dp
                endif
              endif
            endif
            do l = 0,5
              trm2 = rr0*rr0
              rhojiloc = rhojiloc + rr12*denpar(l+1,1,j,mm)*(trm1)**(l)
              if (lgrad1) then
                drhojiloc = drhojiloc - rr12*denpar(l+1,1,j,mm)*trm2*dble(l)*(trm1)**(l+2)
                if (lgrad2) then
                  trm2 = trm2*trm2
                  drhoji2loc = drhoji2loc + rr12*denpar(l+1,1,j,mm)*trm2*dble(l*(l+2))*(trm1)**(l+4)
                  if (lgrad3) then
                    trm2 = trm2*rr0*rr0
                    drhoji3loc = drhoji3loc - rr12*denpar(l+1,1,j,mm)*trm2*dble(l*(l+2)*(l+4))*(trm1)**(l+6)
                  endif
                endif
              endif
            enddo
          elseif (ndenfn(j,mm).eq.12) then
!
!  Baskes
!
            apt = denpar(1,1,j,mm)
            bpt = denpar(2,1,j,mm)/denpar(3,1,j,mm)
            cpt = denpar(3,1,j,mm)
            trm1 = apt*oci*exp(-bpt*(r-cpt))
            rhojiloc = trm1
            if (lgrad1) then
              drhojiloc = - bpt*trm1*rk
              if (lgrad2) then
                drhoji2loc = bpt*trm1*rk2*(bpt + rk)
                if (lgrad3) then
                  drhoji3loc = - bpt*trm1*rk2*rk*(bpt*(bpt+rk) - rk2)
                endif
              endif
            endif
          elseif (ndenfn(j,mm).eq.14) then
!
!  Fractional power law
!
            rpt = denpar(2,1,j,mm)
            trm1 = oci*denpar(1,1,j,mm)*(rk**rpt)
            rhojiloc = trm1
            if (lgrad1) then
              trm1 = - rpt*trm1*rk2
              drhojiloc = trm1
              if (lgrad2) then
                trm2 = - trm1*(rpt+2.0_dp)*rk2
                drhoji2loc = trm2
                if (lgrad3) then
                  drhoji3loc = - trm2*(rpt+4.0_dp)*rk2
                endif
              endif
            endif
          elseif (ndenfn(j,mm).eq.13) then
!
!  VBO
!
            apt = denpar(1,1,j,mm)  ! c
            bpt = denpar(2,1,j,mm)  ! sigma
            cpt = denpar(3,1,j,mm)  ! gamma
            dpt = denpar(4,1,j,mm)  ! r0
            ept = denpar(5,1,j,mm)  ! delta
!
            rx = sqrt(r/ept)
            rxm1 = 1.0_dp/(1.0_dp - rx)
            exptrm = exp(-cpt*rxm1)
!  Multiply c by sigma and normalisation factor
            apt = apt*bpt*exp(cpt/(1.0_dp - sqrt(dpt/ept)))
!
            rhojiloc = apt*exptrm
            if (lgrad1) then
              dbdx = - rhojiloc*cpt*rxm1*rxm1
              dxdr = 0.5_dp/(ept*rx)
              trm1 = rk*dbdx*dxdr
              drhojiloc = trm1
              if (lgrad2) then
                d2bdx2 = - dbdx*rxm1*(cpt*rxm1 - 2.0_dp)
                d2xdr2 = - dxdr*dxdr/rx
                trm2 = d2bdx2*dxdr*dxdr + dbdx*d2xdr2
                drhoji2loc = rk2*(trm2 - trm1)
                if (lgrad3) then
                  d3bdx3 = - d2bdx2*cpt*rxm1*rxm1 - dbdx*(4.0_dp*cpt*rxm1 - 6.0_dp)*rxm1*rxm1
                  d3xdr3 = - 3.0_dp*d2xdr2*dxdr/rx
                  trm3 = d3bdx3*dxdr*dxdr*dxdr + 3.0_dp*d2bdx2*d2xdr2*dxdr + dbdx*d3xdr3
                  drhoji3loc = rk2*rk*(trm3 - 3.0_dp*rk*(trm2 - trm1))
                endif
              endif
            endif
          elseif (ndenfn(j,mm).eq.15) then
!
!  Spline (cubic)
!
            apt  = denpar(1,1,j,mm)
            bpt  = denpar(2,1,j,mm)
            cpt  = denpar(3,1,j,mm)
            dpt  = denpar(4,1,j,mm)
            rmin = denpar(5,1,j,mm)
            r0   = denpar(6,1,j,mm)
            if (r.lt.r0.and.r.ge.rmin) then
              rd = (r - r0)
              rhojiloc = apt*rd*rd*rd + bpt*rd*rd + cpt*rd + dpt
              if (lgrad1) then
                trm3 = 3.0_dp*apt*rk
                trm2 = 2.0_dp*bpt*rk
                trm1 = cpt*rk
                drhojiloc = trm3*rd*rd + trm2*rd + trm1
                if (lgrad2) then
                  trm3 = trm3*rk
                  trm2 = trm2*rk
                  trm1 = trm1*rk*rk
                  drhoji2loc = trm3*rd*(2.0_dp - rd*rk) + trm2*(1.0_dp - rd*rk) - trm1
                  if (lgrad3) then
                    trm3 = trm3*rk
                    trm2 = 3.0_dp*trm2*rk*rk
                    trm1 = trm1*rk*rk
                    drhoji3loc = trm3*(2.0_dp - 6.0_dp*rd*rk + 3.0_dp*rd*rd*rk*rk) - trm2*(1.0_dp - rd*rk) + 3.0_dp*trm1
                  endif
                endif
              endif
            else
              rhojiloc = 0.0_dp
              if (lgrad1) then
                drhojiloc = 0.0_dp
                if (lgrad2) then
                  drhoji2loc = 0.0_dp
                  if (lgrad3) then
                    drhoji3loc = 0.0_dp
                  endif
                endif
              endif
            endif
          elseif (ndenfn(j,mm).eq.16) then
!
!  Morse squared
!
            apt = denpar(1,1,j,mm)   ! De
            bpt = denpar(2,1,j,mm)   ! a
            cpt = denpar(3,1,j,mm)   ! r0
            trm1 = exp(-bpt*(r-cpt))
            trm2 = (1 - trm1)
            trm3 = (trm2*trm2 - 1.0_dp)
            trm4 = oci*apt*apt*trm3*trm3
            rhojiloc = trm4
            if (lgrad1) then
              dtrm2 = bpt*trm1
              dtrm3 = 2.0_dp*trm2*dtrm2
              drhojiloc = 2.0_dp*oci*apt*apt*trm3*dtrm3*rk
              if (lgrad2) then
                d2trm2 = - bpt*bpt*trm1
                d2trm3 = 2.0_dp*(trm2*d2trm2 + dtrm2*dtrm2)
                drhoji2loc = 2.0_dp*oci*apt*apt*((trm3*d2trm3 + dtrm3*dtrm3)*rk2 - trm3*dtrm3*rk2*rk)
                if (lgrad3) then
                  d3trm2 = bpt*bpt*bpt*trm1
                  d3trm3 = 2.0_dp*(3.0_dp*dtrm2*d2trm2 + trm2*d3trm2)
                  drhoji3loc = 2.0_dp*oci*apt*apt*((3.0_dp*dtrm3*d2trm3 + trm3*d3trm3)*rk2*rk &
                               - 3.0_dp*(trm3*d2trm3 + dtrm3*dtrm3)*rk2*rk2 &
                               + 3.0_dp*trm3*dtrm3*rk2*rk2*rk)
                endif
              endif
            endif
          endif
!
          rhojilocsum  = rhojilocsum  + rhojiloc*alloy
          if (lgrad1) then
            drhojilocsum = drhojilocsum + drhojiloc*alloy
            if (lgrad2) then
              drhoji2locsum = drhoji2locsum + drhoji2loc*alloy
              if (lgrad3) then
                drhoji3locsum = drhoji3locsum + drhoji3loc*alloy
              endif
            endif
          endif
        enddo
!
!  Add contributions 
!
        if (rtaper.gt.1.0d-12) then
!
!  EAM taper
!
          term0 = rhojilocsum*tpfn
          if (lgrad1) then
            term1 = drhojilocsum*tpfn + rhojilocsum*dtpfn
            if (lgrad2) then
              term2 = drhoji2locsum*tpfn + 2.0_dp*drhojilocsum*dtpfn + rhojilocsum*d2tpfn
              if (lgrad3) then
                term3 = drhoji3locsum*tpfn + 3.0_dp*(drhoji2locsum*dtpfn + drhojilocsum*d2tpfn) + rhojilocsum*d3tpfn
              endif
            endif
          endif
        else
          term0 = rhojilocsum
          if (lgrad1) then
            term1 = drhojilocsum
            if (lgrad2) then
              term2 = drhoji2locsum
              if (lgrad3) then
                term3 = drhoji3locsum
              endif
            endif
          endif
        endif
!
!  Add derivatives to arrays for ji case
!
        rhoji = rhoji + term0
        if (lgrad1) then
          drhoji(1) = drhoji(1) + term1*x
          drhoji(2) = drhoji(2) + term1*y
          drhoji(3) = drhoji(3) + term1*z
          if (lstr) then
            do is1 = 1,nstrains
              ns1 = nstrptr(is1)
              drhojis(is1) = drhojis(is1) + term1*dr2ds(ns1)
            enddo
          endif
          if (lgrad2) then
            drhoji2(1) = drhoji2(1) + term2*x*x + term1
            drhoji2(2) = drhoji2(2) + term2*x*y
            drhoji2(3) = drhoji2(3) + term2*x*z
            drhoji2(4) = drhoji2(4) + term2*y*y + term1
            drhoji2(5) = drhoji2(5) + term2*y*z
            drhoji2(6) = drhoji2(6) + term2*z*z + term1
            if (lstr) then
              ind = 0
              do is1 = 1,nstrains
                ns1 = nstrptr(is1)
                do is2 = 1,is1
                  ns2 = nstrptr(is2)
                  ind = ind + 1
                  drhoji2s(ind) = drhoji2s(ind) + term2*dr2ds(ns1)*dr2ds(ns2) + term1*d2r2ds2(ns2,ns1)
                enddo
                drhoji2m(is1,1) = drhoji2m(is1,1) + term2*dr2ds(ns1)*x + term1*d2r2dsdx(ns1,1)
                drhoji2m(is1,2) = drhoji2m(is1,2) + term2*dr2ds(ns1)*y + term1*d2r2dsdx(ns1,2)
                drhoji2m(is1,3) = drhoji2m(is1,3) + term2*dr2ds(ns1)*z + term1*d2r2dsdx(ns1,3)
              enddo
            endif
            if (lgrad3) then
              drhoji3(1)  = drhoji3(1)  + term3*x*x*x
              drhoji3(2)  = drhoji3(2)  + term3*x*x*y
              drhoji3(3)  = drhoji3(3)  + term3*x*x*z
              drhoji3(4)  = drhoji3(4)  + term3*x*y*y
              drhoji3(5)  = drhoji3(5)  + term3*x*y*z
              drhoji3(6)  = drhoji3(6)  + term3*x*z*z
              drhoji3(7)  = drhoji3(7)  + term3*y*y*y
              drhoji3(8)  = drhoji3(8)  + term3*y*y*z
              drhoji3(9)  = drhoji3(9)  + term3*y*z*z
              drhoji3(10) = drhoji3(10) + term3*z*z*z
            endif
          endif
        endif
      endif
    endif
!
!  Terms for j->i
!
    if (natj.eq.neamnat(mm).and.(ntypj.eq.neamtyp(mm).or.neamtyp(mm).eq.0)) then
      if (neamnat2(mm).eq.0.or.(neamnat2(mm).eq.nati.and.(ntypi.eq.neamtyp2(mm).or.neamtyp2(mm).eq.0))) then
        rhoijlocsum = 0.0_dp
        drhoijlocsum = 0.0_dp
        drhoij2locsum = 0.0_dp
        drhoij3locsum = 0.0_dp
        alloy = eamalloy(1,mm)
        do j = 1,ndenfncomp(mm)
!
!  Density functional forms
!
          npt = nint(denpar(6,1,j,mm))
          if (ndenfn(j,mm).eq.1) then
!
!  Power-law
!
            trm1 = oci*denpar(1,1,j,mm)*(rk**npt)
            rhoijloc = trm1
            if (lgrad1) then
              trm1 = - dble(npt)*trm1*rk2
              drhoijloc = trm1
              if (lgrad2) then
                trm2 = - trm1*dble(npt+2)*rk2
                drhoij2loc = trm2
                if (lgrad3) then
                  drhoij3loc = - trm2*dble(npt+4)*rk2
                endif
              endif
            endif
          elseif (ndenfn(j,mm).eq.2) then
!
!  Exponential
!
            apt = denpar(1,1,j,mm)
            bpt = denpar(2,1,j,mm)
            cpt = denpar(3,1,j,mm)
            trm1 = apt*ocj*exp(-bpt*(r-cpt))
            if (npt.ne.0) then
              rn1 = r**(npt-1)
              rhoijloc = trm1*rn1*r
              if (lgrad1) then
                drhoijloc = trm1*rn1*(npt*rk-bpt)
                if (lgrad2) then
                  drhoij2loc = trm1*rn1*rk*(npt*(npt-1)*rk2-bpt*dble(2*npt-1)*rk+bpt*bpt)
                  if (lgrad3) then
                    drhoij3loc = trm1*rn1*rk2*(dble(npt*(npt-2)*(npt-4))*rk2*rk-3.0_dp*bpt*dble(npt*(npt-3)+1) &
                      *rk2+3.0_dp*bpt*bpt*dble(npt-1)*rk-bpt*bpt*bpt)
                  endif
                endif
              endif
            else
              rhoijloc = trm1
              if (lgrad1) then
                drhoijloc = - bpt*trm1*rk
                if (lgrad2) then
                  drhoij2loc = bpt*trm1*rk2*(bpt+rk)
                  if (lgrad3) then
                    drhoij3loc = - bpt*trm1*rk2*rk*(bpt*(bpt+rk)-rk2)
                  endif
                endif
              endif
            endif
          elseif (ndenfn(j,mm).eq.3) then
!
!  Gaussian
!
            apt = denpar(1,1,j,mm)
            bpt = denpar(2,1,j,mm)
            cpt = denpar(3,1,j,mm)
            rd = (r - cpt)
            rd2 = rd*rd
            trm1 = apt*ocj*exp(-bpt*rd2)
            if (npt.ne.0) then
              rn1 = r**(npt-1)
              rhoijloc = trm1*rn1*r
              if (lgrad1) then
                drhoijloc = trm1*rn1*(npt*rk-2.0_dp*bpt*rd)
                if (lgrad2) then
                  drhoij2loc = trm1*rn1*rk*(npt*(npt-2)*rk2-2.0_dp*bpt*rd*dble(2*npt-1)*rk+2.0_dp*bpt*( &
                    2.0_dp*bpt*rd2-1.0_dp))
                  if (lgrad3) then
                    drhoij3loc = trm1*rn1*rk2*(dble(npt*(npt-2)*(npt-4))*rk2*rk-6.0_dp*bpt*rd*dble(npt &
                      *npt-3*npt+1)*rk2-6.0_dp*bpt*dble(npt-1)*rk*(1.0_dp-2.0_dp*bpt*rd2)+4.0_dp*bpt*bpt*rd*( &
                      3.0_dp-2.0_dp*bpt*rd2))
                  endif
                endif
              endif
            else
              rhoijloc = trm1
              if (lgrad1) then
                drhoijloc = - 2.0_dp*bpt*rd*trm1*rk
                if (lgrad2) then
                  drhoij2loc = 2.0_dp*trm1*rk2*bpt*(rd*rk-1.0_dp+4.0_dp*bpt*rd2)
                  if (lgrad3) then
                    drhoij3loc = - 4.0_dp*trm1*bpt*rk2*(bpt*rd*rk+rk2)* &
                      (rd*rk-1.0_dp+4.0_dp*bpt*rd2)+2.0_dp*trm1*rk2*bpt*(-rd*rk2+8.0_dp*bpt*rd)*rk
                  endif
                endif
              endif
            endif
          elseif (ndenfn(j,mm).eq.4) then
!
!  Cubic
!
            apt = denpar(1,1,j,mm)
            bpt = denpar(2,1,j,mm)
            if (r.lt.bpt) then
              rd = (r - bpt)
              trm1 = 3.0_dp*apt*rk
              rhoijloc = apt*rd*rd*rd
              if (lgrad1) then
                drhoijloc = trm1*rd*rd
                if (lgrad2) then
                  trm1 = trm1*rk
                  drhoij2loc = trm1*rd*(2.0_dp - rd*rk)
                  if (lgrad3) then
                    trm1 = trm1*rk
                    drhoij3loc = trm1*(2.0_dp - 6.0_dp*rd*rk + 3.0_dp*rd*rd*rk*rk)
                  endif
                endif
              endif
            else
              rhoijloc = 0.0_dp
              if (lgrad1) then
                drhoijloc = 0.0_dp
                if (lgrad2) then
                  drhoij2loc = 0.0_dp
                  if (lgrad3) then
                    drhoij3loc = 0.0_dp
                  endif
                endif
              endif
            endif
          elseif (ndenfn(j,mm).eq.5) then
!
!  Voter-Chen
!
            apt = denpar(1,1,j,mm)
            bpt = denpar(2,1,j,mm)
            etrm = exp(-bpt*r)
            cpt = 2.0_dp**9
            trm1 = apt*etrm*(1.0_dp + cpt*etrm)
            trm2 = apt*etrm*(1.0_dp + 2.0_dp*cpt*etrm)
            rhoijloc = r2*r2*r2*trm1
            if (lgrad1) then
              drhoijloc = r2*r2*(6.0_dp*trm1 - bpt*r*trm2)
              if (lgrad2) then
                trm3 = apt*etrm*(1.0_dp + 4.0_dp*cpt*etrm)
                drhoij2loc = r2*(24.0_dp*trm1 - 11.0_dp*bpt*r*trm2 + bpt*bpt*r2*trm3)
                if (lgrad3) then
                  trm4 = apt*etrm*(1.0_dp + 8.0_dp*cpt*etrm)
                  drhoij3loc = 48.0_dp*trm1 - bpt*r*(57.0_dp*trm2 - bpt*r*(15.0_dp*trm3 - bpt*r*trm4))
                endif
              endif
            endif
          elseif (ndenfn(j,mm).eq.6) then
!
!  Quadratic
!
            apt = denpar(1,1,j,mm)
            bpt = denpar(2,1,j,mm)
            if (r.lt.bpt) then
              rd = (r - bpt)
              trm1 = 2.0_dp*apt*rk
              rhoijloc = apt*rd*rd
              if (lgrad1) then
                drhoijloc = trm1*rd
                if (lgrad2) then
                  trm1 = trm1*rk
                  drhoij2loc = trm1*(1.0_dp - rd*rk)
                  if (lgrad3) then
                    trm1 = 3.0_dp*trm1*rk*rk
                    drhoij3loc = - trm1*(1.0_dp - rd*rk)
                  endif
                endif
              endif
            else
              rhoijloc = 0.0_dp
              if (lgrad1) then
                drhoijloc = 0.0_dp
                if (lgrad2) then
                  drhoij2loc = 0.0_dp
                  if (lgrad3) then
                    drhoij3loc = 0.0_dp
                  endif
                endif
              endif
            endif
          elseif (ndenfn(j,mm).eq.7) then
!
!  Quartic
!
            apt = denpar(1,1,j,mm)
            bpt = denpar(2,1,j,mm)
            if (r.lt.bpt) then
              rd = (r - bpt)
              trm1 = 4.0_dp*apt*rk*rd
              rhoijloc = apt*rd*rd*rd*rd
              if (lgrad1) then
                drhoijloc = trm1*rd*rd
                if (lgrad2) then
                  trm1 = trm1*rk
                  drhoij2loc = trm1*rd*(3.0_dp - rd*rk)
                  if (lgrad3) then
                    trm1 = trm1*rk
                    drhoij3loc = 3.0_dp*trm1*(2.0_dp - 3.0_dp*rd*rk + rd*rd*rk*rk)
                  endif
                endif
              endif
            else
              rhoijloc = 0.0_dp
              if (lgrad1) then
                drhoijloc = 0.0_dp
                if (lgrad2) then
                  drhoij2loc = 0.0_dp
                  if (lgrad3) then
                    drhoij3loc = 0.0_dp
                  endif
                endif
              endif
            endif
          elseif (ndenfn(j,mm).eq.8) then
!
!  Glue
!
            if (r.lt.denpar(1,1,j,mm)) then
              dr = (r - denpar(1,1,j,mm))
              dr2 = dr*dr
              dr3 = dr2*dr
              trm1 = 3.0_dp*denpar(4,1,j,mm)*dr2 + 2.0_dp*denpar(5,1,j,mm)*dr + denpar(6,1,j,mm)
              rhoijloc = denpar(4,1,j,mm)*dr3 + denpar(5,1,j,mm)*dr2 + denpar(6,1,j,mm)*dr + denpar(7,1,j,mm)
              if (lgrad1) then
                drhoijloc = trm1*rk
                if (lgrad2) then
                  trm2 = 6.0_dp*denpar(4,1,j,mm)*dr + 2.0_dp*denpar(5,1,j,mm)
                  drhoij2loc = rk2*(trm2 - rk*trm1)
                  if (lgrad3) then
                    trm3 = 6.0_dp*denpar(4,1,j,mm)
                    drhoij3loc = rk2*rk*(trm3 - 3.0_dp*rk*trm2 + 3.0_dp*rk2*trm1)
                  endif
                endif
              endif
            elseif (r.lt.denpar(2,1,j,mm)) then
              dr = (r - denpar(1,1,j,mm))
              dr2 = dr*dr
              dr3 = dr2*dr
              trm1 = 3.0_dp*denpar(8,1,j,mm)*dr2 + 2.0_dp*denpar(9,1,j,mm)*dr + denpar(10,1,j,mm)
              rhoijloc = denpar(8,1,j,mm)*dr3 + denpar(9,1,j,mm)*dr2 + denpar(10,1,j,mm)*dr + denpar(11,1,j,mm)
              if (lgrad1) then
                drhoijloc = trm1*rk
                if (lgrad2) then
                  trm2 = 6.0_dp*denpar(8,1,j,mm)*dr + 2.0_dp*denpar(9,1,j,mm)
                  drhoij2loc = rk2*(trm2 - rk*trm1)
                  if (lgrad3) then
                    trm3 = 6.0_dp*denpar(8,1,j,mm)
                    drhoij3loc = rk2*rk*(trm3 - 3.0_dp*rk*trm2 + 3.0_dp*rk2*trm1)
                  endif
                endif
              endif
            elseif (r.lt.denpar(3,1,j,mm)) then
              dr = (r - denpar(3,1,j,mm))
              dr2 = dr*dr
              dr3 = dr2*dr
              trm1 = 3.0_dp*denpar(12,1,j,mm)*dr2 + 2.0_dp*denpar(13,1,j,mm)*dr + denpar(14,1,j,mm)
              rhoijloc = denpar(12,1,j,mm)*dr3 + denpar(13,1,j,mm)*dr2 + denpar(14,1,j,mm)*dr + denpar(15,1,j,mm)
              if (lgrad1) then
                drhoijloc = trm1*rk
                if (lgrad2) then
                  trm2 = 6.0_dp*denpar(12,1,j,mm)*dr + 2.0_dp*denpar(13,1,j,mm)
                  drhoij2loc = rk2*(trm2 - rk*trm1)
                  if (lgrad3) then
                    trm3 = 6.0_dp*denpar(12,1,j,mm)
                    drhoij3loc = rk2*rk*(trm3 - 3.0_dp*rk*trm2 + 3.0_dp*rk2*trm1)
                  endif
                endif
              endif
            else
              rhoijloc = 0.0_dp
              if (lgrad1) then
                drhoijloc = 0.0_dp
                if (lgrad2) then
                  drhoij2loc = 0.0_dp
                  if (lgrad3) then
                    drhoij3loc = 0.0_dp
                  endif
                endif
              endif
            endif
          elseif (ndenfn(j,mm).eq.9) then
!
!  eVoter-Chen
!
            call etaper(r,0.0_dp,rmax,etpfn,detpfn,d2etpfn,d3etpfn,.true.,lgrad2,lgrad3)
            detpfn = rk*detpfn
            d2etpfn = rk2*d2etpfn
            d3etpfn = rk2*(rk*d3etpfn - 3.0_dp*(d2etpfn - rk2*detpfn))
            d2etpfn = d2etpfn - rk2*detpfn
!
            apt = denpar(1,1,j,mm)
            bpt = denpar(2,1,j,mm)
            etrm = exp(-bpt*r)
            cpt = 2.0_dp**9
            trm1 = apt*etrm*(1.0_dp + cpt*etrm)
            trm2 = apt*etrm*(1.0_dp + 2.0_dp*cpt*etrm)
            rhoijloc = r2*r2*r2*trm1
            if (lgrad1) then
              drhoijloc = r2*r2*(6.0_dp*trm1 - bpt*r*trm2)
              if (lgrad2) then
                trm3 = apt*etrm*(1.0_dp + 4.0_dp*cpt*etrm)
                drhoij2loc = r2*(24.0_dp*trm1 - 11.0_dp*bpt*r*trm2 + bpt*bpt*r2*trm3)
                if (lgrad3) then
                  trm4 = apt*etrm*(1.0_dp + 8.0_dp*cpt*etrm)
                  drhoij3loc = 48.0_dp*trm1 - bpt*r*(57.0_dp*trm2 - bpt*r*(15.0_dp*trm3 - bpt*r*trm4))
                  drhoij3loc = drhoij3loc*etpfn + 3.0_dp*(drhoij2loc*detpfn + drhoijloc*d2etpfn) + rhoijloc*d3etpfn
                endif
                drhoij2loc = drhoij2loc*etpfn + 2.0_dp*drhoijloc*detpfn + rhoijloc*d2etpfn
              endif
              drhoijloc = drhoijloc*etpfn + rhoijloc*detpfn
            endif
          elseif (ndenfn(j,mm).eq.10) then
!
!  Mei-Davenport
!   
            rr0 = denpar(7,1,j,mm)
            rr12 = 1.0_dp/12.0_dp
            trm1 = rk*rr0
            rr0 = 1.0_dp/denpar(7,1,j,mm)
            rhoijloc = 0.0_dp
            if (lgrad1) then
              drhoijloc = 0.0_dp
              if (lgrad2) then
                drhoij2loc = 0.0_dp
                if (lgrad3) then
                  drhoij3loc = 0.0_dp
                endif
              endif
            endif
            do l = 0,5
              trm2 = rr0*rr0
              rhoijloc = rhoijloc + rr12*denpar(l+1,1,j,mm)*(trm1)**(l)
              if (lgrad1) then
                drhoijloc = drhoijloc - rr12*denpar(l+1,1,j,mm)*trm2*dble(l)*(trm1)**(l+2)
                if (lgrad2) then
                  trm2 = trm2*trm2
                  drhoij2loc = drhoij2loc + rr12*denpar(l+1,1,j,mm)*trm2*dble(l*(l+2))*(trm1)**(l+4)
                  if (lgrad3) then
                    trm2 = trm2*rr0*rr0
                    drhoij3loc = drhoij3loc - rr12*denpar(l+1,1,j,mm)*trm2*dble(l*(l+2)*(l+4))*(trm1)**(l+6)
                  endif
                endif
              endif
            enddo
          elseif (ndenfn(j,mm).eq.12) then
!
!  Baskes
!
            apt = denpar(1,1,j,mm)
            bpt = denpar(2,1,j,mm)/denpar(3,1,j,mm)
            cpt = denpar(3,1,j,mm)
            trm1 = apt*oci*exp(-bpt*(r-cpt))
            rhoijloc = trm1
            if (lgrad1) then
              drhoijloc = - bpt*trm1*rk
              if (lgrad2) then
                drhoij2loc = bpt*trm1*rk2*(bpt + rk)
                if (lgrad3) then
                  drhoij3loc = - bpt*trm1*rk2*rk*(bpt*(bpt+rk) - rk2)
                endif
              endif
            endif
          elseif (ndenfn(j,mm).eq.14) then
!
!  Fractional power law
!
            rpt = denpar(2,1,j,mm)
            trm1 = oci*denpar(1,1,j,mm)*(rk**rpt)
            rhoijloc = trm1
            if (lgrad1) then
              trm1 = - rpt*trm1*rk2
              drhoijloc = trm1
              if (lgrad2) then
                trm2 = - trm1*(rpt+2.0_dp)*rk2
                drhoij2loc = trm2
                if (lgrad3) then
                  drhoij3loc = - trm2*(rpt+4.0_dp)*rk2
                endif
              endif
            endif
          elseif (ndenfn(j,mm).eq.13) then
!
!  VBO
!
            apt = denpar(1,1,j,mm)  ! c
            bpt = denpar(2,1,j,mm)  ! sigma
            cpt = denpar(3,1,j,mm)  ! gamma
            dpt = denpar(4,1,j,mm)  ! r0
            ept = denpar(5,1,j,mm)  ! delta
!
            rx = sqrt(r/ept)
            rxm1 = 1.0_dp/(1.0_dp - rx)
            exptrm = exp(-cpt*rxm1)
!  Multiply c by sigma and normalisation factor
            apt = apt*bpt*exp(cpt/(1.0_dp - sqrt(dpt/ept)))
!
            rhoijloc = apt*exptrm
            if (lgrad1) then
              dbdx = - rhoijloc*cpt*rxm1*rxm1
              dxdr = 0.5_dp/(ept*rx)
              trm1 = rk*dbdx*dxdr
              drhoijloc = trm1
              if (lgrad2) then
                d2bdx2 = - dbdx*rxm1*(cpt*rxm1 - 2.0_dp)
                d2xdr2 = - dxdr*dxdr/rx
                trm2 = d2bdx2*dxdr*dxdr + dbdx*d2xdr2
                drhoij2loc = rk2*(trm2 - trm1)
                if (lgrad3) then
                  d3bdx3 = - d2bdx2*cpt*rxm1*rxm1 - dbdx*(4.0_dp*cpt*rxm1 - 6.0_dp)*rxm1*rxm1
                  d3xdr3 = - 3.0_dp*d2xdr2*dxdr/rx
                  trm3 = d3bdx3*dxdr*dxdr*dxdr + 3.0_dp*d2bdx2*d2xdr2*dxdr + dbdx*d3xdr3
                  drhoij3loc = rk2*rk*(trm3 - 3.0_dp*rk*(trm2 - trm1))
                endif
              endif
            endif
          elseif (ndenfn(j,mm).eq.15) then
!
!  Spline (cubic)
!
            apt  = denpar(1,1,j,mm)
            bpt  = denpar(2,1,j,mm)
            cpt  = denpar(3,1,j,mm)
            dpt  = denpar(4,1,j,mm)
            rmin = denpar(5,1,j,mm)
            r0   = denpar(6,1,j,mm)
            if (r.lt.r0.and.r.ge.rmin) then
              rd = (r - r0)
              rhoijloc = apt*rd*rd*rd + bpt*rd*rd + cpt*rd + dpt
              if (lgrad1) then
                trm3 = 3.0_dp*apt*rk
                trm2 = 2.0_dp*bpt*rk
                trm1 = cpt*rk
                drhoijloc = trm3*rd*rd + trm2*rd + trm1
                if (lgrad2) then
                  trm3 = trm3*rk
                  trm2 = trm2*rk
                  trm1 = trm1*rk*rk
                  drhoij2loc = trm3*rd*(2.0_dp - rd*rk) + trm2*(1.0_dp - rd*rk) - trm1
                  if (lgrad3) then
                    trm3 = trm3*rk
                    trm2 = 3.0_dp*trm2*rk*rk
                    trm1 = trm1*rk*rk
                    drhoij3loc = trm3*(2.0_dp - 6.0_dp*rd*rk + 3.0_dp*rd*rd*rk*rk) - trm2*(1.0_dp - rd*rk) + 3.0_dp*trm1
                  endif
                endif
              endif
            else
              rhoijloc = 0.0_dp
              if (lgrad1) then
                drhoijloc = 0.0_dp
                if (lgrad2) then
                  drhoij2loc = 0.0_dp
                  if (lgrad3) then
                    drhoij3loc = 0.0_dp
                  endif
                endif
              endif
            endif
          elseif (ndenfn(j,mm).eq.16) then
!
!  Morse squared
!
            apt = denpar(1,1,j,mm)   ! De
            bpt = denpar(2,1,j,mm)   ! a
            cpt = denpar(3,1,j,mm)   ! r0
            trm1 = exp(-bpt*(r-cpt))
            trm2 = (1 - trm1)
            trm3 = (trm2*trm2 - 1.0_dp)
            trm4 = oci*apt*apt*trm3*trm3
            rhoijloc = trm4
            if (lgrad1) then
              dtrm2 = bpt*trm1
              dtrm3 = 2.0_dp*trm2*dtrm2
              drhoijloc = 2.0_dp*oci*apt*apt*trm3*dtrm3*rk
              if (lgrad2) then
                d2trm2 = - bpt*bpt*trm1
                d2trm3 = 2.0_dp*(trm2*d2trm2 + dtrm2*dtrm2)
                drhoij2loc = 2.0_dp*oci*apt*apt*((trm3*d2trm3 + dtrm3*dtrm3)*rk2 - trm3*dtrm3*rk2*rk)
                if (lgrad3) then
                  d3trm2 = bpt*bpt*bpt*trm1
                  d3trm3 = 2.0_dp*(3.0_dp*dtrm2*d2trm2 + trm2*d3trm2)
                  drhoij3loc = 2.0_dp*oci*apt*apt*((3.0_dp*dtrm3*d2trm3 + trm3*d3trm3)*rk2*rk &
                               - 3.0_dp*(trm3*d2trm3 + dtrm3*dtrm3)*rk2*rk2 &
                               + 3.0_dp*trm3*dtrm3*rk2*rk2*rk)
                endif
              endif
            endif
          endif
          rhoijlocsum  = rhoijlocsum  + rhoijloc*alloy
          if (lgrad1) then
            drhoijlocsum = drhoijlocsum + drhoijloc*alloy
            if (lgrad2) then
              drhoij2locsum = drhoij2locsum + drhoij2loc*alloy
              if (lgrad3) then
                drhoij3locsum = drhoij3locsum + drhoij3loc*alloy
              endif
            endif
          endif
        enddo
!               
!  Add contributions 
!    
        if (rtaper.gt.1.0d-12) then
!
!  EAM taper
!
          term0 = rhoijlocsum*tpfn
          if (lgrad1) then
            term1 = drhoijlocsum*tpfn + rhoijlocsum*dtpfn
            if (lgrad2) then
              term2 = drhoij2locsum*tpfn + 2.0_dp*drhoijlocsum*dtpfn + rhoijlocsum*d2tpfn
              if (lgrad3) then
                term3 = drhoij3locsum*tpfn + 3.0_dp*(drhoij2locsum*dtpfn + drhoijlocsum*d2tpfn) + rhoijlocsum*d3tpfn
              endif
            endif
          endif
        else
          term0 = rhoijlocsum
          if (lgrad1) then
            term1 = drhoijlocsum
            if (lgrad2) then
              term2 = drhoij2locsum
              if (lgrad3) then 
                term3 = drhoij3locsum
              endif
            endif
          endif
        endif
!
!  Add derivatives to arrays for ji case
!
        rhoij = rhoij + term0
        if (lgrad1) then
          drhoij(1) = drhoij(1) + term1*x
          drhoij(2) = drhoij(2) + term1*y
          drhoij(3) = drhoij(3) + term1*z
          if (lstr) then
            do is1 = 1,nstrains
              ns1 = nstrptr(is1)
              drhoijs(is1) = drhoijs(is1) + term1*dr2ds(ns1)
            enddo
          endif
          if (lgrad2) then
            drhoij2(1) = drhoij2(1) + term2*x*x + term1
            drhoij2(2) = drhoij2(2) + term2*x*y
            drhoij2(3) = drhoij2(3) + term2*x*z
            drhoij2(4) = drhoij2(4) + term2*y*y + term1
            drhoij2(5) = drhoij2(5) + term2*y*z
            drhoij2(6) = drhoij2(6) + term2*z*z + term1
            if (lstr) then
              ind = 0
              do is1 = 1,nstrains
                ns1 = nstrptr(is1)
                do is2 = 1,is1
                  ns2 = nstrptr(is2)
                  ind = ind + 1
                  drhoij2s(ind) = drhoij2s(ind) + term2*dr2ds(ns1)*dr2ds(ns2) + term1*d2r2ds2(ns2,ns1)
                enddo
                drhoij2m(is1,1) = drhoij2m(is1,1) + term2*dr2ds(ns1)*x + term1*d2r2dsdx(ns1,1)
                drhoij2m(is1,2) = drhoij2m(is1,2) + term2*dr2ds(ns1)*y + term1*d2r2dsdx(ns1,2)
                drhoij2m(is1,3) = drhoij2m(is1,3) + term2*dr2ds(ns1)*z + term1*d2r2dsdx(ns1,3)
              enddo
            endif
            if (lgrad3) then
              drhoij3(1)  = drhoij3(1)  + term3*x*x*x
              drhoij3(2)  = drhoij3(2)  + term3*x*x*y
              drhoij3(3)  = drhoij3(3)  + term3*x*x*z
              drhoij3(4)  = drhoij3(4)  + term3*x*y*y
              drhoij3(5)  = drhoij3(5)  + term3*x*y*z
              drhoij3(6)  = drhoij3(6)  + term3*x*z*z
              drhoij3(7)  = drhoij3(7)  + term3*y*y*y
              drhoij3(8)  = drhoij3(8)  + term3*y*y*z
              drhoij3(9)  = drhoij3(9)  + term3*y*z*z
              drhoij3(10) = drhoij3(10) + term3*z*z*z
            endif
          endif
        endif
      endif
    endif
  enddo
#ifdef TRACE
  call trace_out('eamrho')
#endif
!
  return
  end
