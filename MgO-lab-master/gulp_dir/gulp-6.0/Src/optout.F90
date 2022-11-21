  subroutine optout(ifail,nc,fc,gc,oldcell)
!
!  Output final point of optimisation
!
!  12/00 call to rlist added before geometry measurements for safety
!   2/01 extra call to property added if phonons are calculated prior
!        to a defect calculation
!   4/02 units now output with final energy
!   6/03 XML output added
!   9/03 Option to use cell parameters as optimisation variables added
!   3/04 Setting of gradient pointers used in nptr to improve efficiency
!   5/06 Mass now uses species values
!   9/06 Output of cell vectors after opt added
!   3/07 Call to bond changed to GULP_bond
!   4/07 Clean to parallel I/O
!   4/08 Print out of maximum derivative added
!   4/08 Output of ReaxFF charges added
!   2/09 xml removed
!   6/09 Module name changed from three to m_three
!  12/09 Option to output the forces on a region added
!  12/10 Hiding of shells added as an option
!   1/11 Direction labels corrected for fractional internal derivatives
!   1/11 Output modified for force minimisation case
!   3/11 Stresses added as an output option
!   3/11 Stresses moved to end of output and modified to be output regardless
!        of whether this is a constant volume or pressure calculation.
!  10/11 Output of site energies added as a debugging option
!  10/11 Output of stress tensor moved to property
!   1/12 Format statement for charges modified to allow for more atoms
!  10/12 Call to angle and torsion modified to allow for new arguments being passed
!  11/12 pv removed from argument list since it is unused
!  12/12 site_energy keyword added
!   3/14 Call to angle changed to getangles & torsion to gettorsions
!   8/15 Output of ReaxFF charges removed since this duplicates output
!   3/17 Modifications made to allow for new variable order in iopt
!   1/18 Trace added
!   6/18 Output of final strain added for strain cell option
!  11/18 Handling of cell strain option added for non-primitive cells
!   3/19 iopt replaced by ioptindex and iopttype
!   9/19 Rigid molecules added
!  10/19 Rigid molecule modifications made
!  12/19 ifail = 6 handling added for trap on extreme function change
!   7/20 Format for final gnorm output adjusted for enthalpy case
!   7/20 Site energies moved to a separate routine
!   7/20 Units of final energy now reflect the keywords kcal and kjmol
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
!  Copyright Curtin University, 2020
!
!  Julian Gale, CIC, Curtin University, July 2020
!
  use configurations
  use g_constants
  use control
  use current
  use derivatives
  use element
  use four
  use iochannels
  use m_three
  use molecule,      only : nmolasym
  use optimisation
  use parallel
  use species,       only : massspec
  use symmetry
  use terse,         only : ltersederivs, lterseoutcoords, lterseoutcell
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4),                   intent(inout) :: ifail
  integer(i4),                   intent(in)    :: nc
  real(dp),                      intent(in)    :: fc
  real(dp),                      intent(in)    :: gc(*)
  real(dp),                      intent(inout) :: oldcell(*)
!
!  Local variables
!
  character(len=10)                            :: cotype(3)
  character(len=2)                             :: cstype
  character(len=5)                             :: lab
  character(len=5)                             :: namecell(6)
  character(len=8)                             :: unitword
  integer(i4)                                  :: i
  integer(i4)                                  :: i1max
  integer(i4)                                  :: i2max
  integer(i4)                                  :: i2min
  integer(i4)                                  :: ii
  integer(i4)                                  :: inat
  integer(i4)                                  :: itype
  integer(i4)                                  :: j
  integer(i4)                                  :: nangtot
  integer(i4)                                  :: ni
  integer(i4)                                  :: nphitot
  real(dp)                                     :: ad
  real(dp)                                     :: adiff
  real(dp)                                     :: alpo
  real(dp)                                     :: ao
  real(dp)                                     :: ara
  real(dp)                                     :: arao
  real(dp)                                     :: area
  real(dp)                                     :: bd
  real(dp)                                     :: beto
  real(dp)                                     :: bo
  real(dp)                                     :: cd
  real(dp)                                     :: co
  real(dp)                                     :: cellnew(6)
  real(dp)                                     :: gamo
  real(dp)                                     :: padiff
  real(dp)                                     :: pdiff(6)
  real(dp)                                     :: pvdiff
  real(dp)                                     :: rd
  real(dp)                                     :: rdmax
  real(dp)                                     :: rmolfct
  real(dp)                                     :: rv2(3,3)
  real(dp)                                     :: rv3(3,3)
  real(dp)                                     :: sdiff(6)
  real(dp)                                     :: totmass
  real(dp)                                     :: units
  real(dp)                                     :: vdiff
  real(dp)                                     :: vol
  real(dp)                                     :: volf
  real(dp)                                     :: volo
  real(dp)                                     :: volo1
  real(dp)                                     :: volume
  real(dp)                                     :: xci
  real(dp)                                     :: yci
  real(dp)                                     :: zci
  real(dp)                                     :: xd
  real(dp)                                     :: yd
  real(dp)                                     :: zd
  real(dp)                                     :: xdmax
  real(dp)                                     :: ydmax
  real(dp)                                     :: zdmax
  real(dp)                                     :: xsi
  real(dp)                                     :: ysi
  real(dp)                                     :: zsi
!
  data namecell/'a    ','b    ','c    ','alpha','beta','gamma'/
#ifdef TRACE
  call trace_in('optout')
#endif
!
!  Set units for energy output
!
  if (lkcal) then
    units = 1.0_dp/kcaltoev
    unitword = 'kcal/mol'
  elseif (lkjmol) then
    units = 1.0_dp/kjmtoev
    unitword = 'kJ/mol  '
  else
    units = 1.0_dp
    unitword = 'eV      '
  endif
!
  rmolfct = avogadro*1.0d-23
  if (ndim.eq.3) then
    cotype(1) = 'Fractional'
    cotype(2) = 'Fractional'
    cotype(3) = 'Fractional'
  elseif (ndim.eq.2) then
    cotype(1) = 'Fractional'
    cotype(2) = 'Fractional'
    cotype(3) = 'Cartesian '
  elseif (ndim.eq.1) then
    cotype(1) = 'Fractional'
    cotype(2) = 'Cartesian '
    cotype(3) = 'Cartesian '
  else
    cotype(1) = 'Cartesian '
    cotype(2) = 'Cartesian '
    cotype(3) = 'Cartesian '
  endif
!
!  Check gradient norm
!
  gnorm = 0.0_dp
  do i = 1,nvar
    gnorm = gnorm + gc(i)*gc(i)
  enddo
  if (nvar.gt.0) gnorm = sqrt(gnorm)/nvar
!
!  Has the calculation been successful?
!
  if (ioproc) then
    write(ioout,'(/)')
  endif
  if (ifail.eq.3.and.gnorm.lt.0.001_dp) ifail = 0
  if (ifail.lt.0) then
    if (ioproc) then
      write(ioout,'(''  **** CPU limit has been exceeded - restart optimisation ****'',/)')
    endif
  elseif (ifail.eq.1) then
    if (ioproc) then
      write(ioout,'(''  **** Too many failed attempts to optimise ****'')')
    endif
  elseif (ifail.eq.2) then
    if (ioproc) then
      write(ioout,'(''  **** Maximum number of function calls has been reached ****'',/)')
    endif
  elseif (ifail.eq.3) then
    if (ioproc) then
      write(ioout,'(''  **** Conditions for a minimum have not been satisfied. However ****'')')
      write(ioout,'(''  **** no lower point can be found - treat results with caution  ****'')')
      write(ioout,'(''  **** unless gradient norm is small (less than 0.1)             ****'',/)')
    endif
  elseif (ifail.eq.5) then
#ifdef TRACE
    call trace_out('optout')
#endif
    return
  elseif (ifail.eq.6) then
    if (ioproc) then
      write(ioout,'(''  **** Extreme change of function was trapped ****'',/)')
    endif
  elseif (ifail.eq.0) then
    if (ioproc) then
      write(ioout,'(''  **** Optimisation achieved ****'',/)')
    endif
  elseif (ifail.ne.4) then
    if (ioproc) then
      write(ioout,'(''  **** Unexpected termination of optimisation ****'',/)')
    endif
  endif
  if (ifail.ne.4) then
!
!  Output final energy breakdown
!
    if (ioproc) then
      if (lforcemin) then
        write(ioout,'(/,''  Final Gnorm       = '',f16.8)') gnorm
      else
        if (lfree) then
          write(ioout,'(/,''  Final free energy = '',f16.8,1x,a8)') fc*units,unitword
          write(ioout,'(''  Final Gnorm       = '',f16.8)') gnorm
        elseif (abs(press).gt.0.0_dp) then
          write(ioout,'(/,''  Final enthalpy = '',f16.8,1x,a8)') fc*units,unitword
          write(ioout,'(''  Final Gnorm    = '',f16.8)') gnorm
        else
          write(ioout,'(/,''  Final energy = '',f16.8,1x,a8)') fc*units,unitword
          write(ioout,'(''  Final Gnorm  = '',f16.8)') gnorm
        endif
      endif
      call outener
    endif
!
!  Output final geometry and derivatives
!
    if (.not.lterseoutcoords) then
      call outstructure
    endif
  endif
  xd = 0.0_dp
  yd = 0.0_dp
  zd = 0.0_dp
  ad = 0.0_dp
  bd = 0.0_dp
  cd = 0.0_dp
  do i = ncellmin,ncellmax
    if (ioptindex(i).eq.1) xd = gc(i)
    if (ioptindex(i).eq.2) yd = gc(i)
    if (ioptindex(i).eq.3) zd = gc(i)
    if (ioptindex(i).eq.4) ad = gc(i)
    if (ioptindex(i).eq.5) bd = gc(i)
    if (ioptindex(i).eq.6) cd = gc(i)
  enddo
  if (.not.lconv.and..not.ltersederivs) then
    if (ndim.eq.3) then
      call uncell3D(rv,a,b,c,alpha,beta,gamma)
      if (ioproc) then
        if (.not.lterseoutcell) then
          write(ioout,'(''  Final Cartesian lattice vectors (Angstroms) :'',/)')
          do i = 1,3
            write(ioout,'(4x,3f12.6)')(rv(j,i),j=1,3)
          enddo
          if (lstraincell) then
            write(ioout,'(/,''  Final strains :'',/)')
            write(ioout,'(4x,6f12.6)')(straincfg(j,ncf),j=1,6)
          endif
          write(ioout,'(/)')
        endif
        if (loptcellpar) then
          write(ioout,'(''  Final cell parameters and derivatives :'',/)')
          write(ioout,'(''--------------------------------------------------------------------------------'')')
          if (lkcal) then
            write(ioout,'(''       a      '',f14.6,'' Angstrom     dE/a    '',1x,f12.6,'' kcal/mol/Ang'')') a,xd*units
            write(ioout,'(''       b      '',f14.6,'' Angstrom     dE/b    '',1x,f12.6,'' kcal/mol/Ang'')') b,yd*units
            write(ioout,'(''       c      '',f14.6,'' Angstrom     dE/c    '',1x,f12.6,'' kcal/mol/Ang'')') c,zd*units
            write(ioout,'(''       alpha  '',f14.6,'' Degrees      dE/alpha'',1x,f12.6,'' kcal/mol/Deg'')') alpha,ad*units
            write(ioout,'(''       beta   '',f14.6,'' Degrees      dE/beta '',1x,f12.6,'' kcal/mol/Deg'')') beta,bd*units
            write(ioout,'(''       gamma  '',f14.6,'' Degrees      dE/gamma'',1x,f12.6,'' kcal/mol/Deg'')') gamma,cd*units
          elseif (lkjmol) then
            write(ioout,'(''       a      '',f14.6,'' Angstrom     dE/a    '',1x,f12.6,'' kJ/mol/Ang'')') a,xd*units
            write(ioout,'(''       b      '',f14.6,'' Angstrom     dE/b    '',1x,f12.6,'' kJ/mol/Ang'')') b,yd*units
            write(ioout,'(''       c      '',f14.6,'' Angstrom     dE/c    '',1x,f12.6,'' kJ/mol/Ang'')') c,zd*units
            write(ioout,'(''       alpha  '',f14.6,'' Degrees      dE/alpha'',1x,f12.6,'' kJ/mol/Deg'')') alpha,ad*units
            write(ioout,'(''       beta   '',f14.6,'' Degrees      dE/beta '',1x,f12.6,'' kJ/mol/Deg'')') beta,bd*units
            write(ioout,'(''       gamma  '',f14.6,'' Degrees      dE/gamma'',1x,f12.6,'' kJ/mol/Deg'')') gamma,cd*units
          else
            write(ioout,'(''       a      '',f14.6,'' Angstrom     dE/a    '',1x,f12.6,'' eV/Ang'')') a,xd
            write(ioout,'(''       b      '',f14.6,'' Angstrom     dE/b    '',1x,f12.6,'' eV/Ang'')') b,yd
            write(ioout,'(''       c      '',f14.6,'' Angstrom     dE/c    '',1x,f12.6,'' eV/Ang'')') c,zd
            write(ioout,'(''       alpha  '',f14.6,'' Degrees      dE/alpha'',1x,f12.6,'' eV/Deg'')') alpha,ad
            write(ioout,'(''       beta   '',f14.6,'' Degrees      dE/beta '',1x,f12.6,'' eV/Deg'')') beta,bd
            write(ioout,'(''       gamma  '',f14.6,'' Degrees      dE/gamma'',1x,f12.6,'' eV/Deg'')') gamma,cd
          endif
          write(ioout,'(''--------------------------------------------------------------------------------'',/)')
        else
          write(ioout,'(''  Final cell parameters and derivatives :'',/)')
          write(ioout,'(''--------------------------------------------------------------------------------'')')
          if (lkcal) then
            write(ioout,'(''       a      '',f14.6,'' Angstrom     dE/de1(xx)'',1x,f12.6,'' kcal/mol/strain'')') a,xd*units
            write(ioout,'(''       b      '',f14.6,'' Angstrom     dE/de2(yy)'',1x,f12.6,'' kcal/mol/strain'')') b,yd*units
            write(ioout,'(''       c      '',f14.6,'' Angstrom     dE/de3(zz)'',1x,f12.6,'' kcal/mol/strain'')') c,zd*units
            write(ioout,'(''       alpha  '',f14.6,'' Degrees      dE/de4(yz)'',1x,f12.6,'' kcal/mol/strain'')') alpha,ad*units
            write(ioout,'(''       beta   '',f14.6,'' Degrees      dE/de5(xz)'',1x,f12.6,'' kcal/mol/strain'')') beta,bd*units
            write(ioout,'(''       gamma  '',f14.6,'' Degrees      dE/de6(xy)'',1x,f12.6,'' kcal/mol/strain'')') gamma,cd*units
          elseif (lkjmol) then
            write(ioout,'(''       a      '',f14.6,'' Angstrom     dE/de1(xx)'',1x,f12.6,'' kJ/mol/strain'')') a,xd*units
            write(ioout,'(''       b      '',f14.6,'' Angstrom     dE/de2(yy)'',1x,f12.6,'' kJ/mol/strain'')') b,yd*units
            write(ioout,'(''       c      '',f14.6,'' Angstrom     dE/de3(zz)'',1x,f12.6,'' kJ/mol/strain'')') c,zd*units
            write(ioout,'(''       alpha  '',f14.6,'' Degrees      dE/de4(yz)'',1x,f12.6,'' kJ/mol/strain'')') alpha,ad*units
            write(ioout,'(''       beta   '',f14.6,'' Degrees      dE/de5(xz)'',1x,f12.6,'' kJ/mol/strain'')') beta,bd*units
            write(ioout,'(''       gamma  '',f14.6,'' Degrees      dE/de6(xy)'',1x,f12.6,'' kJ/mol/strain'')') gamma,cd*units
          else
            write(ioout,'(''       a      '',f14.6,'' Angstrom     dE/de1(xx)'',1x,f12.6,'' eV/strain'')') a,xd
            write(ioout,'(''       b      '',f14.6,'' Angstrom     dE/de2(yy)'',1x,f12.6,'' eV/strain'')') b,yd
            write(ioout,'(''       c      '',f14.6,'' Angstrom     dE/de3(zz)'',1x,f12.6,'' eV/strain'')') c,zd
            write(ioout,'(''       alpha  '',f14.6,'' Degrees      dE/de4(yz)'',1x,f12.6,'' eV/strain'')') alpha,ad
            write(ioout,'(''       beta   '',f14.6,'' Degrees      dE/de5(xz)'',1x,f12.6,'' eV/strain'')') beta,bd
            write(ioout,'(''       gamma  '',f14.6,'' Degrees      dE/de6(xy)'',1x,f12.6,'' eV/strain'')') gamma,cd
          endif
          write(ioout,'(''--------------------------------------------------------------------------------'',/)')
        endif
      endif
    elseif (ndim.eq.2) then
      if (.not.lterseoutcell.and.ioproc) then
        write(ioout,'(''  Final Cartesian surface vectors (Angstroms) :'',/)')
        do i = 1,2
          write(ioout,'(4x,3f12.6)')(rv(j,i),j=1,3)
        enddo
        if (lstraincell) then
          write(ioout,'(/,''  Final strains :'',/)')
          write(ioout,'(4x,3f12.6)')(straincfg(j,ncf),j=1,3)
        endif
        write(ioout,'(/)')
      endif
      call uncell2D(rv,a,b,alpha)
      if (ioproc) then
        write(ioout,'(''  Final surface cell parameters and derivatives :'',/)')
        write(ioout,'(''--------------------------------------------------------------------------------'')')
        if (lkcal) then
          write(ioout,'(''       a      '',f14.6,'' Angstrom     dE/de1(xx)'',1x,f12.6,'' kcal/mol/strain'')') a,xd*units
          write(ioout,'(''       b      '',f14.6,'' Angstrom     dE/de2(yy)'',1x,f12.6,'' kcal/mol/strain'')') b,yd*units
          write(ioout,'(''       alpha  '',f14.6,'' Degrees      dE/de3(xy)'',1x,f12.6,'' kcal/mol/strain'')') alpha,zd*units
        elseif (lkjmol) then
          write(ioout,'(''       a      '',f14.6,'' Angstrom     dE/de1(xx)'',1x,f12.6,'' kJ/mol/strain'')') a,xd*units
          write(ioout,'(''       b      '',f14.6,'' Angstrom     dE/de2(yy)'',1x,f12.6,'' kJ/mol/strain'')') b,yd*units
          write(ioout,'(''       alpha  '',f14.6,'' Degrees      dE/de3(xy)'',1x,f12.6,'' kJ/mol/strain'')') alpha,zd*units
        else
          write(ioout,'(''       a      '',f14.6,'' Angstrom     dE/de1(xx)'',1x,f12.6,'' eV/strain'')') a,xd
          write(ioout,'(''       b      '',f14.6,'' Angstrom     dE/de2(yy)'',1x,f12.6,'' eV/strain'')') b,yd
          write(ioout,'(''       alpha  '',f14.6,'' Degrees      dE/de3(xy)'',1x,f12.6,'' eV/strain'')') alpha,zd
        endif
        write(ioout,'(''--------------------------------------------------------------------------------'',/)')
      endif
    elseif (ndim.eq.1) then
      if (.not.lterseoutcell.and.ioproc) then
        write(ioout,'(''  Final Cartesian polymer vector (Angstroms) :'',/)')
        write(ioout,'(4x,3f12.6)')(rv(j,1),j=1,3)
        if (lstraincell) then
          write(ioout,'(/,''  Final strain :'',/)')
          write(ioout,'(4x,1f12.6)') straincfg(1,ncf)
        endif
        write(ioout,'(/)')
      endif
      call uncell1D(rv,a)
      if (ioproc) then
        write(ioout,'(''  Final polymer cell parameter and derivative :'',/)')
        write(ioout,'(''--------------------------------------------------------------------------------'')')
        if (lkcal) then
          write(ioout,'(''       a      '',f14.6,'' Angstrom     dE/de1(xx)'',1x,f12.6,'' kcal/mol/strain'')') a,xd*units
        elseif (lkjmol) then
          write(ioout,'(''       a      '',f14.6,'' Angstrom     dE/de1(xx)'',1x,f12.6,'' kJ/mol/strain'')') a,xd*units
        else
          write(ioout,'(''       a      '',f14.6,'' Angstrom     dE/de1(xx)'',1x,f12.6,'' eV/strain'')') a,xd
        endif
        write(ioout,'(''--------------------------------------------------------------------------------'',/)')
      endif
    endif
!
!  If in compare mode then store primitive cell parameters
!
    if (lcomp) then
      if (ndim.eq.3) then
        cellnew(1) = a
        cellnew(2) = b
        cellnew(3) = c
        cellnew(4) = alpha
        cellnew(5) = beta
        cellnew(6) = gamma
      elseif (ndim.eq.2) then
        cellnew(1) = a
        cellnew(2) = b
        cellnew(3) = alpha
      elseif (ndim.eq.1) then
        cellnew(1) = a
      endif
    endif
    if (ndim.eq.3) then
!
!  Calculate primitive cell volume
!
      vol = volume(rv)
      if (ioproc) then
        write(ioout,'(''  Primitive cell volume = '',f20.6,'' Angs**3'')') vol
      endif
!
!  Calculate density
!
      totmass = 0.0_dp
      do i = 1,nasym
        ni = nspecptr(i)
        totmass = totmass + massspec(ni)*occua(i)*dble(neqv(i))
      enddo
      if (vol.gt.1d-12) then
        density = (10.0_dp*totmass)/(vol*rmolfct)
        if (ioproc) then
          write(ioout,'(/,''  Density of cell = '',f13.6,'' g/cm**3'')') density
        endif
      endif
!
!  Recentre
!
      if (ncbl.gt.1) then
        do i = 1,3
          rv2(1,i) = rvcfg(1,i,nc)
          rv2(2,i) = rvcfg(2,i,nc)
          rv2(3,i) = rvcfg(3,i,nc)
        enddo
        if (lstraincell) then
          call strain3D(strain,rv2)
        endif
        call uncentre(rv2)
        call uncell3D(rv2,a,b,c,alpha,beta,gamma)
        if (ioproc) then
          write(ioout,'(/,''  Non-primitive lattice parameters :'',/)')
          write(ioout,'(''  a    = '',f12.6,''  b   = '',f12.6,''  c    = '',f12.6)') a,b,c
          write(ioout,'(''  alpha= '',f12.6,''  beta= '',f12.6,''  gamma= '',f12.6)') alpha,beta,gamma
        endif
      endif
!
!  Calculate full cell volume
!
      volf = volume(rv)*dble(icentfct(ncbl))
      if (ioproc) then
        write(ioout,'(/,''  Non-primitive cell volume = '',f20.6,'' Angs**3'')') volf
        write(ioout,'(/)')
      endif
    endif
    if (ndim.eq.2) then
!
!  Calculate surface cell area
!
      ara = area(rv)
      if (ioproc) then
        write(ioout,'(''  Surface cell area = '',f13.6,'' Angs**2'',/)') ara
      endif
    endif
  endif
!
  if (.not.ltersederivs) then
    if (ioproc) then
      if (ndim.eq.3) then
        write(ioout,'(''  Final internal derivatives :'',/)')
        write(ioout,'(''--------------------------------------------------------------------------------'')')
        write(ioout,'(''   No.  Atomic          a             b             c           Radius'')')
        if (lkcal) then
          write(ioout,'(''        Label        (kcal/mol)    (kcal/mol)    (kcal/mol)  (kcal/m/Angs)'')')
        elseif (lkjmol) then
          write(ioout,'(''        Label         (kJ/mol)      (kJ/mol)      (kJ/mol)    (kJ/m/Angs)'')')
        else
          write(ioout,'(''        Label          (eV)          (eV)          (eV)        (eV/Angs)'')')
        endif
      elseif (ndim.eq.2) then
        write(ioout,'(''  Final internal/Cartesian derivatives :'',/)')
        write(ioout,'(''--------------------------------------------------------------------------------'')')
        write(ioout,'(''   No.  Atomic          a             b             z           Radius'')')
        if (lkcal) then
          write(ioout,'(''        Label       (kcal/mol)    (kcal/mol)  (kcal/m/Angs)  (kcal/m/Angs)'')')
        elseif (lkjmol) then
          write(ioout,'(''        Label        (kJ/mol)      (kJ/mol)    (kJ/m/Angs)    (kj/m/Angs)'')')
        else
          write(ioout,'(''        Label          (eV)          (eV)       (eV/Angs)      (eV/Angs)'')')
        endif
      elseif (ndim.eq.1) then
        write(ioout,'(''  Final internal/Cartesian derivatives :'',/)')
        write(ioout,'(''--------------------------------------------------------------------------------'')')
        write(ioout,'(''   No.  Atomic          a             y             z           Radius'')')
        if (lkcal) then
          write(ioout,'(''        Label       (kcal/mol)  (kcal/m/Angs) (kcal/m/Angs)  (kcal/m/Angs)'')')
        elseif (lkjmol) then
          write(ioout,'(''        Label        (kJ/mol)    (kJ/m/Angs)   (kJ/m/Angs)    (kJ/m/Angs)'')')
        else
          write(ioout,'(''        Label          (eV)       (eV/Angs)     (eV/Angs)      (eV/Angs)'')')
        endif
      else
        write(ioout,'(/,''  Final Cartesian derivatives :'',/)')
        write(ioout,'(''--------------------------------------------------------------------------------'')')
        write(ioout,'(''   No.  Atomic          x             y             z           Radius'')')
        if (lkcal) then
          write(ioout,'(''        Label     (kcal/m/Angs) (kcal/m/Angs) (kcal/m/Angs)  (kcal/m/Angs)'')')
        elseif (lkjmol) then
          write(ioout,'(''        Label      (kJ/m/Angs)   (kJ/m/Angs)   (kJ/m/Angs)    (kJ/m/Angs)'')')
        else
          write(ioout,'(''        Label       (eV/Angs)     (eV/Angs)     (eV/Angs)      (eV/Angs)'')')
        endif
      endif
      write(ioout,'(''--------------------------------------------------------------------------------'')')
    endif
  endif
!
  xdmax = 0.0_dp
  ydmax = 0.0_dp
  zdmax = 0.0_dp
  rdmax = 0.0_dp
!
  do ii = 1,nasymnomol
    i = nasymnomolptr(ii)
    inat = iatn(i)
    itype = natype(i)
!
!  Hide shells?
!
    if (inat.le.maxele.or..not.lhideshells) then
      call label(inat,itype,lab)
      if (lbsmat(i+nsft)) then
        cstype = 'bc'
        if (inat.gt.maxele) cstype = 'bs'
      else
        cstype = 'c '
        if (inat.gt.maxele) cstype = 's '
      endif
      xd = 0.0_dp
      yd = 0.0_dp
      zd = 0.0_dp
      rd = 0.0_dp
!
!  Search through variables to find those relevant to atom i
!
      do j = ninternalmin,ninternalmax
        if (ioptindex(j).eq.ii) then
          if (iopttype(j).eq.iopt_xf) then
            xd = gc(j)*units
          elseif (iopttype(j).eq.iopt_yf) then
            yd = gc(j)*units
          elseif (iopttype(j).eq.iopt_zf) then
            zd = gc(j)*units
          elseif (iopttype(j).eq.iopt_radius) then
            rd = gc(j)*units
          endif
        endif
      enddo
!
      xdmax = max(xdmax,abs(xd))
      ydmax = max(ydmax,abs(yd))
      zdmax = max(zdmax,abs(zd))
      rdmax = max(rdmax,abs(rd))
      if (.not.ltersederivs) then
        if (ioproc) then
          write(ioout,'(i7,1x,a5,1x,a2,4f14.6)') i,lab,cstype,xd,yd,zd,rd
        endif
      endif
    endif
  enddo
  if (ioproc) then
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    write(ioout,'(''  Maximum abs   '',4f14.6)') xdmax,ydmax,zdmax,rdmax
    write(ioout,'(''--------------------------------------------------------------------------------'',/)')
  endif
!
!  Output forces on rigid molecules
!
  if (ioproc.and.lrigid.and..not.ltersederivs.and.nmolasym.gt.0) then
    write(ioout,'(''  Final derivatives for centre of mass of rigid molecules: '',/)')
    if (ndim.eq.3) then
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''   Molecule             a             b             c       '')')
      if (lkcal) then
        write(ioout,'(''                     (kcal/mol)    (kcal/mol)    (kcal/mol)     '')')
      elseif (lkjmol) then
        write(ioout,'(''                      (kJ/mol)      (kJ/mol)      (kJ/mol)     '')')
      else
        write(ioout,'(''                       (eV)          (eV)          (eV)     '')')
      endif
    elseif (ndim.eq.2) then
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''   Molecule             a             b             z       '')')
      if (lkcal) then
        write(ioout,'(''                     (kcal/mol)    (kcal/mol) (kcal/m/Angs)   '')')
      elseif (lkjmol) then
        write(ioout,'(''                      (kJ/mol)      (kJ/mol)   (kJ/m/Angs)   '')')
      else
        write(ioout,'(''                       (eV)          (eV)       (eV/Angs)   '')')
      endif
    elseif (ndim.eq.1) then
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''   Molecule             a             y             z       '')')
      if (lkcal) then
        write(ioout,'(''                     (kcal/mol) (kcal/m/Angs) (kcal/m/Angs)   '')')
      elseif (lkjmol) then
        write(ioout,'(''                      (kJ/mol)   (kJ/m/Angs)   (kJ/m/Angs)   '')')
      else
        write(ioout,'(''                       (eV)       (eV/Angs)     (eV/Angs)   '')')
      endif
    else
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''   Molecule             x             y             z       '')')
      if (lkcal) then
        write(ioout,'(''                  (kcal/m/Angs) (kcal/m/Angs) (kcal/m/Angs)    '')')
      elseif (lkjmol) then
        write(ioout,'(''                   (kJ/m/Angs)   (kJ/m/Angs)   (kJ/m/Angs)    '')')
      else
        write(ioout,'(''                    (eV/Angs)     (eV/Angs)     (eV/Angs)    '')')
      endif
    endif
    write(ioout,'(''--------------------------------------------------------------------------------'')')
!
    do i = 1,nmolasym
      xd = 0.0_dp
      yd = 0.0_dp
      zd = 0.0_dp
!
!  Search through variables to find those relevant to atom i
!
      do j = ninternalmin,ninternalmax
        if (ioptindex(j).eq.i) then
          if (iopttype(j).eq.iopt_xcom) then
            xd = gc(j)*units
          elseif (iopttype(j).eq.iopt_ycom) then
            yd = gc(j)*units
          elseif (iopttype(j).eq.iopt_zcom) then
            zd = gc(j)*units
          endif
        endif
      enddo
      write(ioout,'(i7,9x,3f14.6)') i,xd,yd,zd
    enddo
    write(ioout,'(''--------------------------------------------------------------------------------'',/)')
!
    write(ioout,'(''  Final derivatives for quaternions of rigid molecules: '',/)')
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    write(ioout,'(''   Molecule             Q1            Q2            Q3      '')')
    if (lkcal) then
      write(ioout,'(''                     (kcal/mol)    (kcal/mol)    (kcal/mol) '')')
    elseif (lkjmol) then
      write(ioout,'(''                      (kJ/mol)      (kJ/mol)      (kJ/mol)  '')')
    else
      write(ioout,'(''                       (eV)          (eV)          (eV)     '')')
    endif
    write(ioout,'(''--------------------------------------------------------------------------------'')')
!
    do i = 1,nmolasym
      xd = 0.0_dp
      yd = 0.0_dp
      zd = 0.0_dp
!
!  Search through variables to find those relevant to atom i
!
      do j = ninternalmin,ninternalmax
        if (ioptindex(j).eq.i) then
          if (iopttype(j).eq.iopt_xqtn) then
            xd = gc(j)*units
          elseif (iopttype(j).eq.iopt_yqtn) then
            yd = gc(j)*units
          elseif (iopttype(j).eq.iopt_zqtn) then
            zd = gc(j)*units
          endif
        endif
      enddo
      write(ioout,'(i7,9x,3f14.6)') i,xd,yd,zd
    enddo
    write(ioout,'(''--------------------------------------------------------------------------------'')')
  endif
!
!  Option to output forces on a region
!
  if (ioproc.and.lregionforce) then
    write(ioout,'(''  Forces acting on regions : '',/)')
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    write(ioout,'(''  Region No             x               y               z      '')')
    if (lkcal) then
      write(ioout,'(''                  (kcal/m/Angs)   (kcal/m/Angs)   (kcal/m/Angs)  '')')
    elseif (lkjmol) then
      write(ioout,'(''                   (kJ/m/Angs)     (kJ/m/Angs)     (kJ/m/Angs)  '')')
    else
      write(ioout,'(''                    (eV/Angs)       (eV/Angs)       (eV/Angs)  '')')
    endif
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    do i = 1,nregions(ncf)
      write(ioout,'(2x,i8,5x,3f16.6)') i,xregdrv(i)*units,yregdrv(i)*units,zregdrv(i)*units
    enddo
    write(ioout,'(''--------------------------------------------------------------------------------'',/)')
  endif
!
!  Option to compare initial and final structures
!
  if (lcomp.and.ioproc) then
    write(ioout,'(''  Comparison of initial and final structures : '',/)')
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    write(ioout,'(''  Parameter   Initial value   Final value   Difference    Units      Percent'')')
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    if (.not.lconv.and.ndim.gt.0) then
      if (ndim.eq.3) then
        cellnew(1) = a
        cellnew(2) = b
        cellnew(3) = c
        cellnew(4) = alpha
        cellnew(5) = beta
        cellnew(6) = gamma
!
!  Calculate old cell volume
!
        ao = oldcell(1)
        bo = oldcell(2)
        co = oldcell(3)
        alpo = oldcell(4)
        beto = oldcell(5)
        gamo = oldcell(6)
        call cell3D(rv3,ao,bo,co,alpo,beto,gamo)
        volo1 = volume(rv3)
        call uncentre(rv3)
        call uncell3D(rv3,ao,bo,co,alpo,beto,gamo)
        oldcell(1) = ao
        oldcell(2) = bo
        oldcell(3) = co
        oldcell(4) = alpo
        oldcell(5) = beto
        oldcell(6) = gamo
        volo = volume(rv3)
        vol = vol*volo/volo1
        vdiff = vol - volo
        pvdiff = 100.0_dp*vdiff/volo
        write(ioout,'(4x,''Volume'',5x,f12.6,2x,f12.6,1x,f12.6,4x,''Angs**3  '',f8.2)')  &
          volo,vol,vdiff,pvdiff
        i1max = 3
        i2min = 4
        i2max = 6
      elseif (ndim.eq.2) then
        cellnew(1) = a
        cellnew(2) = b
        cellnew(3) = alpha
!
!  Calculate old cell area
!
        ao = oldcell(1)
        bo = oldcell(2)
        alpo = oldcell(3)
        call cell2D(rv3,ao,bo,alpo)
        arao = area(rv3)
        adiff = ara - arao
        padiff = 100.0_dp*adiff/ara
        write(ioout,'(4x,''Area  '',5x,f12.6,2x,f12.6,1x,f12.6,4x,''Angs**2  '',f8.2)') &
          arao,ara,adiff,padiff
        i1max = 2
        i2min = 3
        i2max = 3
      elseif (ndim.eq.1) then
        cellnew(1) = a
        i1max = 1
      endif
!
!  Cell parameters first
!
      do i = 1,nstrains
        sdiff(i) = cellnew(i) - oldcell(i)
        pdiff(i) = 100.0_dp*sdiff(i)/oldcell(i)
      enddo
      do i = 1,i1max
        write(ioout,'(4x,a5,6x,f12.6,2x,f12.6,3x,f10.6,4x,''Angstroms'',f8.2)') &
          namecell(i),oldcell(i),cellnew(i),sdiff(i),pdiff(i)
      enddo
      if (ndim.gt.1) then
        do i = i2min,i2max
          write(ioout,'(4x,a5,6x,f12.6,2x,f12.6,3x,f10.6,4x,''Degrees'',2x,f8.2)') &
            namecell(i),oldcell(i),cellnew(i),sdiff(i),pdiff(i)
        enddo
      endif
    endif
!
!  Fractional or cartesian coordinates
!
    do i = 1,nasym
!
!  Hide shells?
!
      if (iatn(i).le.maxele.or..not.lhideshells) then
        xci = xcfg(i+nsft)
        yci = ycfg(i+nsft)
        zci = zcfg(i+nsft)
        xsi = xstore(i)
        ysi = ystore(i)
        zsi = zstore(i)
        sdiff(1) = xci - xsi
        sdiff(2) = yci - ysi
        sdiff(3) = zci - zsi
        if (abs(xsi).gt.0.0_dp) then
          if (ndim.ge.1) then
!
!  Check if coordinate has moved across boundary
!
            if ((xci/xsi).lt.0.0_dp) then
              xsi = xsi + 1.0_dp
            endif
            sdiff(1) = sdiff(1) + 2.0_dp
            sdiff(1) = dmod(sdiff(1),1.0_dp)
            if (sdiff(1).gt.0.5_dp) then
              sdiff(1) = - (sdiff(1) - 1.0_dp)
              if (xsi.lt.0.01_dp) xsi = xsi + 1.0_dp
            endif
          endif
          pdiff(1) = 100.0_dp*sdiff(1)/xsi
        else
          pdiff(1) = 0.0_dp
        endif
        if (abs(ysi).gt.0.0_dp) then
          if (ndim.ge.2) then
!
!  Check if coordinate has moved across boundary
!
            if ((yci/ysi).lt.0.0_dp) then
              ysi = ysi + 1.0_dp
            endif
            sdiff(2) = sdiff(2) + 2.0_dp
            sdiff(2) = dmod(sdiff(2),1.0_dp)
            if (sdiff(2).gt.0.5_dp) then
              sdiff(2) = - (sdiff(2) - 1.0_dp)
              if (ysi.lt.0.01_dp) ysi = ysi + 1.0_dp
            endif
          endif
          pdiff(2) = 100.0_dp*sdiff(2)/ysi
        else
          pdiff(2) = 0.0_dp
        endif
        if (abs(zsi).gt.0.0_dp) then
          if (ndim.eq.3) then
!
!  Check if coordinate has moved across boundary
!
            if ((zci/zsi).lt.0.0_dp) then
              zsi = zsi + 1.0_dp
            endif
            sdiff(3) = sdiff(3) + 2.0_dp
            sdiff(3) = dmod(sdiff(3),1.0_dp)
            if (sdiff(3).gt.0.5_dp) then
              sdiff(3) = - (sdiff(3) - 1.0_dp)
              if (zsi.lt.0.01_dp) zsi=zsi + 1.0_dp
            endif
          endif
          pdiff(3) = 100.0_dp*sdiff(3)/zsi
        else
          pdiff(3) = 0.0_dp
        endif
        write(ioout,'(1x,i6,1x,''x'',6x,f12.6,2x,f12.6,3x,f10.6,4x,a10,f7.2)') &
          i,xstore(i),xci,sdiff(1),cotype(1),pdiff(1)
        write(ioout,'(1x,i6,1x,''y'',6x,f12.6,2x,f12.6,3x,f10.6,4x,a10,f7.2)') &
          i,ystore(i),yci,sdiff(2),cotype(2),pdiff(2)
        write(ioout,'(1x,i6,1x,''z'',6x,f12.6,2x,f12.6,3x,f10.6,4x,a10,f7.2)') &
          i,zstore(i),zci,sdiff(3),cotype(3),pdiff(3)
      endif
    enddo
    write(ioout,'(''--------------------------------------------------------------------------------'',/)')
  endif
!
!  Option to print out final angles for three-body terms
!
  if (lopt) then
    call rlist
    if (ioproc) then
      if (lbond) call GULP_bond
      if (ldist) call distance
      if (nthb.gt.0.and.langle) call getangles(ioout,0_i4,nangtot)
      if (nfor.gt.0.and.ltors) call gettorsions(ioout,0_i4,nphitot)
    endif
  endif
#ifdef TRACE
  call trace_out('optout')
#endif
!
  return
  end
