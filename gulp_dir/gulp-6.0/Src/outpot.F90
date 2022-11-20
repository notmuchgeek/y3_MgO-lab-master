  subroutine outpot
!
!  Outputs interatomic potentials 
!
!  Potential type for onebody:
!
!   Only type is ashift
!
!  Potential type is given by nptype for twobody:
!
!   1 = buckingham
!   2 = lennard-jones
!   3 = morse
!   4 = morse-coulomb subtracted
!   5 = harmonic
!   6 = harmonic-coulomb subtracted
!   7 = general potential
!   8 = spring (core-shell only)
!   9 = coulomb
!  10 = buckingham 4 range
!  11 = spline
!  12 = lennard-jones - epsilon/sigma form : sigma = r at E=0
!  13 = lennard-jones - epsilon/sigma form : sigma = r at minimum
!  14 = BSM - breathing shell model
!  15 = Stillinger-Weber 2 body potential
!  16 = Inverse gaussian potential
!  17 = BSM - exponential
!  18 = Damped dispersion
!  19 = Many-body potential for metals
!  20 = Rydberg / Rose-Smith-Guinea-Ferrante potential
!  21 = Lennard-Jones with epsilon/sigma input : ESFF combination rules
!  22 = qtaper (short range Coulomb taper)
!  23 = polynomial 
!  24 = qerfc (Coulomb with erfc)
!  25 = CovExp (Covalent-Expotential form)
!  26 = Fermi-Dirac form
!  27 = Lennard-Jones buffered
!  28 = Squared harmonic potential
!  29 = Squared harmonic with Coulomb correction
!  30 = Tsuneyuki Coulomb correction potential
!  31 = BSM - single exponential breathing shell model
!  32 = Stillinger-Weber - Jiang & Brown modification
!  33 = Cosh spring potential
!  34 = EAM potential shift
!  35 = Poly harmonic 
!  36 = qoverr2
!  37 = force constant
!  38 = sr_glue
!  39 = morse with etaper
!  40 = morse-coulomb subtracted with etaper
!  41 = Mei-Davenport twobody
!  42 = erferfc potential
!  43 = reperfc potential
!  44 = erf potential
!  45 = Baskes twobody potential
!  46 = VBO_twobody potential
!  47 = exppowers
!  48 = Grimme_C6
!  49 = cfm_harmonic
!  50 = cfm_gaussian
!  51 = cfm_power
!  52 = cfm_fermi
!  53 = gcoulomb potential
!  54 = Becke_Johnson_C6
!  55 = Baskes twobody potential with a4 term
!  56 = ZBL potential
!  57 = MM3 Buckingham
!  58 = MM3 stretch
!  59 = lennard-jones - epsilon/sigma form : fixed coefficients
!  60 = Buffered Lennard-Jones
!  61 = Slater
!  62 = Constant J2 potential
!  63 = Stillinger-Weber 2 body potential - general form
!  64 = Screened Coulomb potential
!
!  Potential type is given by nthrty for three-body:
!
!   1 = conventional harmonic (k2/k3/k4)
!   2 = exponentially decaying harmonic
!   3 = Axilrod-Teller
!   4 = Exponential three-body form
!   5 = Stillinger-Weber form
!   6 = Bcross (bond-bond cross term)
!   7 = Urey-Bradley
!   8 = exponentially decaying - Vessal form
!   9 = cosine harmonic (k2/k3/k4)
!  10 = Murrell-Mottram
!  11 = BAcross - theta
!  12 = Linear-three
!  13 = Bcoscross
!  14 = Stillinger-Weber - Jiang & Brown modification
!  15 = Hydrogen-bond
!  16 = Equatorial
!  17 = UFF3
!  18 = BAcoscross
!  19 = 3Coulomb
!  20 = exp2
!  21 = g3Coulomb
!  22 = BAGcross
!  23 = BALcross
!  24 = MM3 angle bending
!  25 = Stillinger-Weber with Garofalini form
!  26 = j3
!  27 = PPP 3-body
!
!  Potential type is given by nforty for four-body:
!
!   1 = standard torsion
!   2 = Ryckaert-Bellemanns
!   3 = Out of plane
!   4 = ESFF torsion
!   5 = Harmonic
!   6 = Exponentially decaying standard
!   7 = Exponentially decaying ESFF
!   8 = Tapered standard
!   9 = Tapered ESFF
!  10 = Angle cross term
!  11 = Inversion
!  12 = Inversion squared
!  13 = UFF4
!  14 = Angle-angle cross potential
!  15 = UFFoop
!  16 = Cos angle - cos angle cross potential
!  17 = Torsion - cosine angle cross term
!
!  Called from gulp and fit
!
!   1/95 Intra/intermolecular option added for 3 and 4 body terms
!   1/95 K3 added for three-body harmonic potential
!   2/95 K3 and K4 added for harmonic potential
!   2/95 phi0 added for standard torsional potential 
!   2/95 Exponential and SW 3-body forms added
!   2/95 Bonded specification for 3 and 4 body potentials added
!   3/95 Bcross potential added
!   3/95 SW2 potential added
!   4/95 AB combination rules added for Lennard-Jones potential
!   3/96 Urey-Bradley potential added
!   4/96 Exponentially decaying Vessal form added
!   5/96 Inverse gaussian potential added
!   6/96 General theta0 added to SW3
!   2/97 BSM exponential and damped dispersion added
!   3/97 Outofplane potential added
!   4/97 Many-body potential added
!   7/97 EAM densities added
!   7/97 EAM functionals added
!   7/97 Rose-SGF potential added
!   3/98 Cosine based harmonic three body potential added
!   4/98 ESFF form of Lennard-Jones combination rules added
!   5/98 qtaper potential added
!   6/98 Murrell-Mottram potential added
!   8/98 ESFF torsional potential added
!  10/98 Bond-Angle cross potential added
!   1/99 Modifications to allow for 1-4 only potentials added
!   8/99 Linear-threebody term added
!  10/99 Cubic EAM density added
!   7/00 CovExp added
!   2/01 Fermi-Dirac form added
!   5/01 Output of 3-body potentials modified to allow for
!        minimum cut-offs
!   5/02 Scaling of shifts added
!   5/02 Brenner potential added
!   7/02 Out of plane K4 added
!   8/02 L-J buffered added
!  10/02 Torharm potential added
!  10/02 ReaxFF forcefield added
!   4/03 Exponentially decaying torsion added
!   4/03 Tapered torsion potentials added
!   7/03 Checking removed to separate routine
!   7/03 Bond order potentials added
!  11/03 ndennat/ndentyp replaced
!  11/03 Output of EAM alloy factors added
!   3/04 Squared harmonic potential added
!   7/04 Output for bond order charge potential added
!   9/04 Output for bond order charge self energy added
!  11/04 Output for torangle potential added
!  11/04 Output for six-body potentials added
!   2/05 Rho added for boselfenergy
!   4/05 Cosh spring potential added
!   6/05 Error in output of coefficients for Ryckaert-Bellemanns potl fixed
!   7/05 Output of constant terms for EAM functionals sqrt & power added
!   7/05 Output for Brenner 1990 potential added
!   9/05 Voter form of density for EAM added
!  10/05 Modified to handle numerical EAM functional
!  10/05 Hydrogen-bond potential added
!  10/05 Modified to handle Voter style tapering of Morse
!  10/05 Inversion outofplane potential added
!  11/05 EAM potential shift added
!  12/05 Equatorial ESFF three-body term added
!   2/06 Quadratic and quartic densities added
!   3/06 Modified to allow for density component number
!   3/06 Poly harmonic potential added
!   4/06 Species specific density added
!   4/06 Error in output for taper option fixed
!   5/06 Modified to handle introduction of neamfnspec
!   6/06 Squared inversion added
!  11/06 qoverr2 potential added
!   1/07 UFF3 potential added
!   1/07 UFF4 potential added
!   1/07 force constant potential added
!   2/07 Output of bond type for 4 & 6 body potentials added
!   2/07 End of CML output moved inside conditional statement
!   2/07 Johnson and Glue EAM functionals added
!   2/07 Glue density added
!   3/07 sr_glue potential added
!   4/07 Cyclic/exocyclic label now output for torsions
!   4/07 Units of linear-threebody corrected to eV in output
!   5/07 Exponential taper added
!   5/07 Morse with etaper added
!   5/07 eVoter EAM density added
!   5/07 Foiles EAM functional added
!   7/07 Plane potential added
!   8/07 Format statement for maximum range of potentials changed
!  10/07 Angle-angle cross potential added
!  11/07 Mei-Davenport potential added
!  12/07 Error in type setting for output of sixbody potential labels fixed
!   3/08 erferfc, erfpot and reperfc potentials added
!   4/08 Minimum cutoffs added for out of plane
!   5/08 UFF oop potential added
!   5/08 Output of bond type controls added for twobody pots
!  10/08 Output of taper function name corrected for MDF case
!  10/08 MEAM modifications added - denpar increased in dimension
!  11/08 bacoscross form added
!  11/08 xcosangleangle potential added
!  11/08 torcosangle potential added
!  11/08 Output format statement for bond-bond cross potential K & standard
!        torsion K modified from g10.5 to g10.4
!  11/08 Baskes form of exponential density added
!  11/08 Baskes form of functional added
!  11/08 Output modified to accommodate MEAM densities
!  12/08 eammeamcoeff array changed from density to function related
!  12/08 rho0 added to Baskes functional
!   1/09 Output for Murrell-Mottram potential now includes c0 - c10
!   1/09 Arrays apot, bpot, cpot & dpot replaced by new 2-D array twopot
!   1/09 leshift/lgshift logicals introduced to replace use of tpot(3/4,) for energy
!        and gradient shifts
!   1/09 Baskes twobody potential added
!   1/09 VBO_twobody potential added
!   1/09 VBO density added 
!   1/09 VBO functional added 
!   2/09 Modified to accommodate new version of FoX and gulp_cml
!   3/09 3coulomb potential added
!   4/09 MEAM type output added
!   4/09 Output of MEAM screening parameters added
!   6/09 Module name changed from three to m_three
!   7/09 exp2 potential added
!   7/09 exppowers potential added
!   8/09 Grimme_C6 potential added
!   9/09 Belashchenko EAM functional added
!   9/09 Maximum order of polynomial potential increased to 8
!   1/10 One-body potentials added
!   2/10 Central force model potentials added
!   3/10 Format of output statements adjusted to accommodate
!        uo2_tiwary library
!   5/10 gcoulomb potential added
!   6/10 Bond order header wrapped with if statements for theta
!        option by C. Fisher
!   6/10 Modified to handle optional second atom type for bond-order
!        attractive and repulsive terms.
!   9/10 EDIP force field added
!  10/10 Output of EDIP parameters added
!  10/10 EDIP linear threebody modifications added
!   3/11 namepot moved to module
!   9/11 Twobody potential output format changed for 4th term to g9.3
!  10/11 Fractional power density added
!   1/12 Output for mmexc = 4 case added for missing non-polynomial case
!   1/12 Output for mmexc = 5 case added
!   6/12 Numerical EAM density label added for output
!  10/12 Modified to allow for OpenKIM potentials
!   2/13 bagcross potential added
!   5/13 Custom bond type added
!   6/13 balcross potential added
!   6/13 Igarashi EAM functional added
!   7/13 Cubic spline EAM density added
!   7/13 Improper torsion sub-option added
!   8/13 Spline EAM functional added
!   1/14 EDIP2p added
!   3/14 Fractional bond orders added
!   7/14 lEDIP3mod added
!   8/14 MEAM screening parameters made species specific
!   8/14 Baskes potential output modified
!   8/14 MEAM screening made dependent on trio of atom types
!   9/14 Baskes potential output tied
!   9/14 ZBL potential added
!   1/15 k3 and k4 terms added for cosine form of threebody potential
!   1/15 Modified Stillinger-Weber 3-body added
!   1/15 Bug in output for bond-order attractive potential fixed
!   2/15 MM3buck added
!   2/15 Taper output corrected for MDF
!   3/15 MM3angle potential added
!   3/15 Twobody output now in separate routines
!   3/15 namepot renamed to name2pot
!   3/15 Three and four body potentials now moved to subroutine
!   8/15 Garofalini form of sw3 added
!   9/15 BOdcoeff replaced by extended BOccoeff array
!   9/15 Kumagai variant of Tersoff added
!   3/16 Bond order coordination potential added
!   7/16 KIM descriptors now output in debug mode
!   1/18 Modified for MMP potential
!   3/18 nmeamcombotype = 3 added
!   6/18 Grammar correction to output
!   8/18 7th order polynominal taper added
!   8/18 Modified for version 2 of OpenKIM
!   9/18 Modified for multiple OpenKIM models
!   1/19 EDIP sigma/pi labels corrected
!   3/19 Smoothing of ReaxFF output added
!   8/19 Unused data removed
!  11/19 ppp3body added
!   3/20 Morse squared added
!  11/20 Tersoff reorganised
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
!  Julian Gale, CIC, Curtin University, November 2020
!
  use bondorderdata
  use brennerdata
  use configurations
  use g_constants
  use control
  use current
  use eam
  use EDIPdata
  use element,        only : maxele
  use four
  use iochannels
  use kim_models
  use m_three
  use molecule
  use one
  use parallel
  use plane
  use reaxFFdata
  use shells
  use shifts
  use six
  use species
  use terse,          only : ltersepotentials
  use two
  implicit none
!
!  Local variables
!
  integer(i4),      dimension(:), allocatable  :: iptyp
  character(len=5), dimension(:), allocatable  :: cfg_string
  character(len=1)                             :: cstyp1
  character(len=1)                             :: cstyp2
  character(len=1)                             :: cstyp3
  character(len=1)                             :: cstyp4
  character(len=1)                             :: cstyp5
  character(len=1)                             :: cstyp6
  character(len=5)                             :: lab1
  character(len=5)                             :: lab2
  character(len=5)                             :: lab3
  character(len=5)                             :: lab4
  character(len=5)                             :: lab5
  character(len=5)                             :: lab6
  character(len=9)                             :: botyword(11)
  character(len=12)                            :: denfntyp
  integer(i4)                                  :: i
  integer(i4)                                  :: icfg
  integer(i4)                                  :: ii
  integer(i4)                                  :: ind
  integer(i4)                                  :: j
  integer(i4)                                  :: k
  integer(i4)                                  :: nat1
  integer(i4)                                  :: nat2
  integer(i4)                                  :: nat3
  integer(i4)                                  :: nat4
  integer(i4)                                  :: nat5
  integer(i4)                                  :: nat6
  integer(i4)                                  :: nbotyptr
  integer(i4)                                  :: nm
  integer(i4)                                  :: nmpt(3)
  integer(i4)                                  :: ntype1
  integer(i4)                                  :: ntype2
  integer(i4)                                  :: ntype3
  integer(i4)                                  :: ntype4
  integer(i4)                                  :: ntype5
  integer(i4)                                  :: ntype6
  integer(i4)                                  :: order
  integer(i4)                                  :: status
  logical                                      :: lout
  real(dp)                                     :: dcutp
!
  data botyword/'         ','single   ','double   ','triple   ', &
                'quadruple','resonant ','amide    ','custom   ', &
                'half     ','quarter  ','third    '/
!
  dcutp = 50.0_dp
!
!  Beyond this point everything is I/O
!
  if (.not.ioproc.or.ltersepotentials) return
!
  if (none.gt.0.or.npote.gt.0.or.lbrenner.or.lreaxFF.or.lEDIP.or.lkim_model) then
!**********************************
!  Output interatomic potentials  *
!**********************************
    if (cutp.ne.dcutp) then
      write(ioout,'(/,''  Maximum range for interatomic potentials = '',f16.6,'' Angstroms'')') cutp
    endif
    if ((tapermax - tapermin).gt.0.01_dp) then
      if (tapertype.eq.1) then
        write(ioout,'(/,''  Taper potentials to zero over '',f8.4,'' Angstroms using polynomial'')') (tapermax - tapermin)
      elseif (tapertype.eq.3) then
        write(ioout,'(/,''  Taper potentials to zero over '',f8.4,'' Angstroms using Voter m = '',f8.4)') tapermax,taperm
      elseif (tapertype.eq.4) then
        write(ioout,'(/,''  Taper potentials to zero over '',f8.4,'' Angstroms using exponential'')') tapermax
      elseif (tapertype.eq.5) then
        write(ioout,'(/,''  Taper potentials to zero over '',f8.4,'' Angstroms using M-D-F'')') (tapermax - tapermin)
      elseif (tapertype.eq.6) then
        write(ioout,'(/,''  Taper potentials to zero over '',f8.4,'' Angstroms using 7th order ReaxFF polynominal'')') tapermax
      else
        write(ioout,'(/,''  Taper potentials to zero over '',f8.4,'' Angstroms using cosine'')') (tapermax - tapermin)
      endif
    endif
    if (lc6) then
      if (lc6one) then
        write(ioout,'(/,''  C6 terms to be calculated in real/reciprocal space by one-centre decomposition'')')
      else
        write(ioout,'(/,''  C6 terms to be calculated in real and reciprocal space '')')
      endif
    endif
    if (lbrenner) then
      if (nbrennertype.eq.3) then
        if (lbrennersplineh) then
          write(ioout,'(/,''  Brenner et al 2002 potential to be used with explicit splines for P'',/)')
        else
          write(ioout,'(/,''  Brenner et al 2002 potential to be used with pretabulated spline for P'',/)')
        endif
      elseif (nbrennertype.eq.1) then
        write(ioout,'(/,''  Brenner 1990 and Dyson-Smith 1996 potential to be used '',/)')
      endif
    endif
    if (lreaxFF) then
      write(ioout,'(/,''  ReaxFF forcefield to be used'',/)')
      write(ioout,'(''  ReaxFF Coulomb cutoff = '',f8.4,'' Ang'')') reaxFFcutoffQ
      write(ioout,'(''  ReaxFF VDW     cutoff = '',f8.4,'' Ang'')') reaxFFcutoffVDW
      write(ioout,'(''  ReaxFF H-bond  cutoff = '',f8.4,'' Ang'',/)') reaxFFrhtol
      write(ioout,'(''  ReaxFF pairwise bond order       threshold = '',f16.14)') reaxFFtol
      write(ioout,'(''  ReaxFF angle/torsion bond order  threshold = '',f16.14)') reaxFFatol
      write(ioout,'(''  ReaxFF bond order double product threshold = '',f16.14)') reaxFFatol2
      write(ioout,'(''  ReaxFF bond order triple product threshold = '',f16.14)') reaxFFatol3
      write(ioout,'(''  ReaxFF hydrogen-bond bond order  threshold = '',f16.14,/)') reaxFFhtol
      if (lreaxFFlpsmooth) then
        write(ioout,'(''  ReaxFF lone pair energy smoothed with exponent = '',f10.4)') reaxFFlpsmooth
      endif
      if (lreaxFFovsmooth) then
        write(ioout,'(''  ReaxFF overcoordination smoothed with exponent = '',f10.4)') reaxFFlpsmooth
      endif
      if (lreaxFFlpsmooth.or.lreaxFFovsmooth) then
        write(ioout,'(/)')
      endif
    endif
!
!  OpenKIM
!
    if (lkim_model) then
      write(ioout,'(/,''  OpenKIM model :'',/)')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'('' Model : name / configurations successfully coupled with'')')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      allocate(cfg_string(ncfg))
!
      do nm = 1,nkimmodel
        write(ioout,'(i6,'' : '',a72)') nm,kim_model_name(nm)(1:72)
        ii = 0
        do icfg = 1,ncfg
          if (lkim_model_cfg_OK(icfg,nm)) then
            ii = ii + 1
            call itow(cfg_string(ii),icfg,5_i4)
          endif
        enddo
        write(ioout,'(6x,'' : '',12a5)') (cfg_string(j),j=1,ii)
      enddo
!
      deallocate(cfg_string)
      write(ioout,'(''--------------------------------------------------------------------------------'')')
    endif
    if (lEDIP) then
      write(ioout,'(/,''  EDIP forcefield to be used'')')
      if (nEDIPspec.gt.0) then
!
!  Coordination number
!
        write(ioout,'(/,''  EDIP Coordination Number :'',/)')
        write(ioout,'(''--------------------------------------------------------------------------------'')')
        write(ioout,'('' Atom Types   Sigma    alpha      f_low      f_high     p_low      p_high       '')') 
        write(ioout,'(''               or                 (Ang)      (Ang)      (Ang)      (Ang)        '')') 
        write(ioout,'(''  1     2     Pi       Zdih       Zrep       c0                                 '')')
        write(ioout,'(''                                             (Ang)                              '')') 
        write(ioout,'(''--------------------------------------------------------------------------------'')')
        ind = 0
        do i = 1,nEDIPspec
          nat1 = natEDIPspec(i)
          ntype1 = ntypEDIPspec(i)
          call label(nat1,ntype1,lab1)
          cstyp1 = 'c'
          if (nat1.gt.maxele) cstyp1 = 's'
          do j = 1,i
            ind = ind + 1
            if (lEDIPpairOK(ind)) then
              nat2 = natEDIPspec(j)
              ntype2 = ntypEDIPspec(j)
              call label(nat2,ntype2,lab2)
              cstyp2 = 'c'
              if (nat2.gt.maxele) cstyp2 = 's'
              if (lEDIPpi(ind)) then
                write(ioout,'(2(a5,a1,1x),''pi'',4x,f10.5,4(1x,f10.5))') &
                      lab1,cstyp1,lab2,cstyp2,EDIPalpha(ind),EDIPflow(ind),EDIPfhigh(ind),EDIPplow(ind),EDIPphigh(ind)
                write(ioout,'(20x,f10.5,2(1x,f10.5))') &
                      EDIPZdih(ind),EDIPZrep(ind),EDIPc0(ind)
              else
                write(ioout,'(2(a5,a1,1x),''sigma'',1x,f10.5,2(1x,f10.5))') &
                      lab1,cstyp1,lab2,cstyp2,EDIPalpha(ind),EDIPflow(ind),EDIPfhigh(ind)
              endif
            endif
          enddo
        enddo
        write(ioout,'(''--------------------------------------------------------------------------------'')')
!
!  Twobody energy
!
        write(ioout,'(/,''  EDIP Twobody Energy :'',/)')
        write(ioout,'(''--------------------------------------------------------------------------------'')')
        write(ioout,'(''Atom  Types   Epsilon     p     B       Beta      Sigma       a       a_prime   '')')
        write(ioout,'(''  1     2      (eV)           (Ang)     (Ang)     (Ang)     (Ang)     (Ang)     '')')
        write(ioout,'(''--------------------------------------------------------------------------------'')')
        ind = 0
        do i = 1,nEDIPspec
          nat1 = natEDIPspec(i)
          ntype1 = ntypEDIPspec(i)
          call label(nat1,ntype1,lab1)
          cstyp1 = 'c'
          if (nat1.gt.maxele) cstyp1 = 's'
          do j = 1,i
            ind = ind + 1
            if (lEDIPpairOK(ind)) then
              nat2 = natEDIPspec(j)
              ntype2 = ntypEDIPspec(j)
              call label(nat2,ntype2,lab2)
              cstyp2 = 'c'
              if (nat2.gt.maxele) cstyp2 = 's'
              write(ioout,'(2(a5,a1,1x),f10.5,f5.3,5(1x,f9.4))') &
                    lab1,cstyp1,lab2,cstyp2,EDIP2epsilon(ind),EDIP2p(ind),EDIP2B(ind),EDIP2beta(ind), &
                    EDIP2sigma(ind),EDIP2a(ind),EDIP2aprime(ind)
            endif
          enddo
        enddo
        write(ioout,'(''--------------------------------------------------------------------------------'')')
!
!  Threebody energy - original form
!
        lout = .true.
        do i = 1,nEDIPspec
          nat1 = natEDIPspec(i)
          ntype1 = ntypEDIPspec(i)
          call label(nat1,ntype1,lab1)
          cstyp1 = 'c'
          if (nat1.gt.maxele) cstyp1 = 's'
          ind = 0
          do j = 1,nEDIPspec
            nat2 = natEDIPspec(j)
            ntype2 = ntypEDIPspec(j)
            call label(nat2,ntype2,lab2)
            cstyp2 = 'c'
            if (nat2.gt.maxele) cstyp2 = 's'
            do k = 1,j
              ind = ind + 1
              if (lEDIPtriadOK(ind,i).and.lEDIP3orig(ind,i)) then
                if (lout) then
                  write(ioout,'(/,''  EDIP Threebody Energy (Original form) :'',/)')
                  write(ioout,'(''--------------------------------------------------------------------------------'')')
                  write(ioout,'('' Atom  Types               Lambda0/          Mu/                 Eta/           '')')
                  write(ioout,'(''  1     2     3            Gamma             Gamma_prime         q              '')')
                  write(ioout,'(''--------------------------------------------------------------------------------'')')
                  lout = .false.
                endif
                nat3 = natEDIPspec(k)
                ntype3 = ntypEDIPspec(k)
                call label(nat3,ntype3,lab3)
                cstyp3 = 'c'
                if (nat3.gt.maxele) cstyp3 = 's'
                write(ioout,'(3(a5,a1,1x),f15.6,2(2x,f15.6))') &
                      lab1,cstyp1,lab2,cstyp2,lab3,cstyp3,EDIP3lambda0(ind,i),EDIP3lambdap(ind,i),EDIP3Z0(ind,i)
                write(ioout,'(21x,f15.6,2(2x,f15.6))') EDIP3gamma0(ind,i),EDIP3gammap(ind,i),EDIP3q(ind,i)
              endif
            enddo
          enddo
        enddo
        if (.not.lout) then
          write(ioout,'(''--------------------------------------------------------------------------------'')')
        endif
!
!  Threebody energy - new form, unmodified
!
        lout = .true.
        do i = 1,nEDIPspec
          nat1 = natEDIPspec(i)
          ntype1 = ntypEDIPspec(i)
          call label(nat1,ntype1,lab1)
          cstyp1 = 'c'
          if (nat1.gt.maxele) cstyp1 = 's'
          ind = 0
          do j = 1,nEDIPspec
            nat2 = natEDIPspec(j)
            ntype2 = ntypEDIPspec(j)
            call label(nat2,ntype2,lab2)
            cstyp2 = 'c'
            if (nat2.gt.maxele) cstyp2 = 's'
            do k = 1,j
              ind = ind + 1
              if (lEDIPtriadOK(ind,i).and..not.lEDIP3orig(ind,i).and..not.lEDIP3mod(ind,i)) then
                if (lout) then
                  write(ioout,'(/,''  EDIP Threebody Energy (New form) :'',/)')
                  write(ioout,'(''--------------------------------------------------------------------------------'')')
                  write(ioout,'('' Atom  Types               Lambda0/          Lambda_prime/       Z0/            '')')
                  write(ioout,'(''  1     2     3            Gamma             Gamma_prime         q              '')')
                  write(ioout,'(''--------------------------------------------------------------------------------'')')
                  lout = .false.
                endif
                nat3 = natEDIPspec(k)
                ntype3 = ntypEDIPspec(k)
                call label(nat3,ntype3,lab3)
                cstyp3 = 'c'
                if (nat3.gt.maxele) cstyp3 = 's'
                write(ioout,'(3(a5,a1,1x),f15.6,2(2x,f15.6))') &
                      lab1,cstyp1,lab2,cstyp2,lab3,cstyp3,EDIP3lambda0(ind,i),EDIP3lambdap(ind,i),EDIP3Z0(ind,i)
                write(ioout,'(21x,f15.6,2(2x,f15.6))') EDIP3gamma0(ind,i),EDIP3gammap(ind,i),EDIP3q(ind,i)
              endif
            enddo
          enddo
        enddo
        if (.not.lout) then
          write(ioout,'(''--------------------------------------------------------------------------------'')')
        endif
!
!  Threebody energy - new form, modified
!
        lout = .true.
        do i = 1,nEDIPspec
          nat1 = natEDIPspec(i)
          ntype1 = ntypEDIPspec(i)
          call label(nat1,ntype1,lab1)
          cstyp1 = 'c'
          if (nat1.gt.maxele) cstyp1 = 's'
          ind = 0
          do j = 1,nEDIPspec
            nat2 = natEDIPspec(j)
            ntype2 = ntypEDIPspec(j)
            call label(nat2,ntype2,lab2)
            cstyp2 = 'c'
            if (nat2.gt.maxele) cstyp2 = 's'
            do k = 1,j
              ind = ind + 1
              if (lEDIPtriadOK(ind,i).and.lEDIP3mod(ind,i)) then
                if (lout) then
                  write(ioout,'(/,''  EDIP Threebody Energy (New form, modified):'',/)')
                  write(ioout,'(''--------------------------------------------------------------------------------'')')
                  write(ioout,'('' Atom  Types               Lambda0/          Lambda_prime/       Z0/            '')')
                  write(ioout,'(''  1     2     3            Gamma             Gamma_prime         q              '')')
                  write(ioout,'(''                                             K_q2                               '')')
                  write(ioout,'(''--------------------------------------------------------------------------------'')')
                  lout = .false.
                endif
                nat3 = natEDIPspec(k)
                ntype3 = ntypEDIPspec(k)
                call label(nat3,ntype3,lab3)
                cstyp3 = 'c'
                if (nat3.gt.maxele) cstyp3 = 's'
                write(ioout,'(3(a5,a1,1x),f15.6,2(2x,f15.6))') &
                      lab1,cstyp1,lab2,cstyp2,lab3,cstyp3,EDIP3lambda0(ind,i),EDIP3lambdap(ind,i),EDIP3Z0(ind,i)
                write(ioout,'(21x,f15.6,2(2x,f15.6))') EDIP3gamma0(ind,i),EDIP3gammap(ind,i),EDIP3q(ind,i)
                write(ioout,'(36x,2x,f15.6)') EDIP3kq2(ind,i)
              endif
            enddo
          enddo
        enddo
        if (.not.lout) then
          write(ioout,'(''--------------------------------------------------------------------------------'')')
        endif
      endif
    endif
!
!  One-body
!
    if (none.gt.0) then
      write(ioout,'(/,''  One-body self-energy potentials :'',/)')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''  Atom Type    Potential         Energy (eV) '')')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      do i = 1,none
        nat1 = nspec11(i)
        ntype1 = nptyp11(i)
        call label(nat1,ntype1,lab1)
        cstyp1 = 'c'
        if (nat1.gt.maxele) cstyp1 = 's'
        write(ioout,'(2x,a5,1x,a1,6x,''Ashift'',4x,f20.8)') lab1,cstyp1,onepot(i)
      enddo
      write(ioout,'(''--------------------------------------------------------------------------------'')')
    endif
    if (npote.gt.0) then
!
!  Two-body
!
      call outtwo
    endif
  endif
!******************
!  EAM densities  *
!******************
  if (neamspec.gt.0) then
    lout = .true.
    do i = 1,neamspec
      nat1 = neamnat(i)
      ntype1 = neamtyp(i)
      call label(nat1,ntype1,lab1)
      cstyp1 = 'c'
      if (nat1.gt.maxele) cstyp1 = 's'
      nat2 = neamnat2(i)
      if (nat2.gt.0) then
        ntype2 = neamtyp2(i)
        call label(nat2,ntype2,lab2)
        cstyp2 = 'c'
        if (nat2.gt.maxele) cstyp1 = 's'
      endif
      do j = 1,ndenfncomp(i)
        if (neammeamorder(j,i).eq.1) then
          if (lout) then
            write(ioout,'(/,''  Embedded Atom Model Densities :'',/)')
            write(ioout,'(''--------------------------------------------------------------------------------'')')
            write(ioout,'(''  Atom(s)       Functional Form          A          B          C         n      '')')
            write(ioout,'(''--------------------------------------------------------------------------------'')')
            lout = .false.
          endif
          if (ndenfn(j,i).gt.0) then
            if (ndenfn(j,i).eq.1) then
              denfntyp = 'Power Law   '
            elseif (ndenfn(j,i).eq.2) then
              denfntyp = 'Exponential '
            elseif (ndenfn(j,i).eq.3) then
              denfntyp = 'Gaussian    '
            elseif (ndenfn(j,i).eq.4) then
              denfntyp = 'Cubic       '
            elseif (ndenfn(j,i).eq.5) then
              denfntyp = 'Voter-Chen  '
            elseif (ndenfn(j,i).eq.6) then
              denfntyp = 'Quadratic   '
            elseif (ndenfn(j,i).eq.7) then
              denfntyp = 'Quartic     '
            elseif (ndenfn(j,i).eq.8) then
              denfntyp = 'Glue        '
            elseif (ndenfn(j,i).eq.9) then
              denfntyp = 'eVoter-Chen '
            elseif (ndenfn(j,i).eq.10) then
              denfntyp = 'Mei-Davenprt'
            elseif (ndenfn(j,i).eq.11) then
              denfntyp = 'Numerical'
            elseif (ndenfn(j,i).eq.12) then
              denfntyp = 'Baskes      '
            elseif (ndenfn(j,i).eq.13) then
              denfntyp = 'VBO         '
            elseif (ndenfn(j,i).eq.14) then
              denfntyp = 'Fractl power'
            elseif (ndenfn(j,i).eq.15) then
              denfntyp = 'Spline      '
            elseif (ndenfn(j,i).eq.16) then
              denfntyp = 'Morse2      '
            endif
            if (ndenfn(j,i).eq.8) then
              if (nat2.gt.0) then
                write(ioout,'(2x,a4,1x,a1,1x,a4,1x,a1,1x,a12,4x,g14.4,2(1x,g10.4))') &
                  lab1(1:4),cstyp1,lab2(1:4),cstyp2,denfntyp,(denpar(k,1,j,i),k=1,3)
              else
                write(ioout,'(2x,a4,1x,a1,8x,a12,4x,g14.4,2(1x,g10.4))') &
                  lab1(1:4),cstyp1,denfntyp,(denpar(k,1,j,i),k=1,3)
              endif
              write(ioout,'(16x,4(f12.6,1x))') (denpar(k,1,j,i),k=4,7)
              write(ioout,'(16x,4(f12.6,1x))') (denpar(k,1,j,i),k=8,11)
              write(ioout,'(16x,4(f12.6,1x))') (denpar(k,1,j,i),k=12,15)
            elseif (ndenfn(j,i).eq.10) then
              if (nat2.gt.0) then
                write(ioout,'(2x,a4,1x,a1,1x,a4,1x,a1,1x,a12,6x,3(f11.6,1x))') &
                  lab1(1:4),cstyp1,lab2(1:4),cstyp2,denfntyp,(denpar(k,1,j,i),k=1,3)
              else
                write(ioout,'(2x,a4,1x,a1,8x,a12,6x,3(f11.6,1x))') &
                  lab1(1:4),cstyp1,denfntyp,(denpar(k,1,j,i),k=1,3)
              endif
              write(ioout,'(34x,3(f11.6,1x),f8.4)') (denpar(k,1,j,i),k=4,7)
            elseif (ndenfn(j,i).eq.13) then
              if (nat2.gt.0) then
                write(ioout,'(2x,a4,1x,a1,1x,a4,1x,a1,1x,a12,6x,3(f11.6,1x))') &
                  lab1(1:4),cstyp1,lab2(1:4),cstyp2,denfntyp,(denpar(k,1,j,i),k=1,3)
              else
                write(ioout,'(2x,a4,1x,a1,8x,a12,6x,3(f11.6,1x))') &
                  lab1(1:4),cstyp1,denfntyp,(denpar(k,1,j,i),k=1,3)
              endif
              write(ioout,'(34x,2(f11.6,1x))') (denpar(k,1,j,i),k=4,5)
            elseif (ndenfn(j,i).eq.15) then
              if (nat2.gt.0) then
                write(ioout,'(2x,a4,1x,a1,1x,a4,1x,a1,1x,a12,6x,3(f11.6,1x))') &
                  lab1(1:4),cstyp1,lab2(1:4),cstyp2,denfntyp,(denpar(k,1,j,i),k=1,3)
              else
                write(ioout,'(2x,a4,1x,a1,8x,a12,6x,3(f11.6,1x))') &
                  lab1(1:4),cstyp1,denfntyp,(denpar(k,1,j,i),k=1,3)
              endif
              write(ioout,'(34x,3(f11.6,1x))') (denpar(k,1,j,i),k=4,6)
            else
              if (nat2.gt.0) then
                write(ioout,'(2x,a4,1x,a1,1x,a4,1x,a1,1x,a12,4x,g14.4,2(1x,g10.4),4x,i2)') &
                  lab1(1:4),cstyp1,lab2(1:4),cstyp2,denfntyp,(denpar(k,1,j,i),k=1,3),nint(denpar(6,1,j,i))
              else
                write(ioout,'(2x,a4,1x,a1,8x,a12,4x,g14.4,2(1x,g10.4),4x,i2)') &
                  lab1(1:4),cstyp1,denfntyp,(denpar(k,1,j,i),k=1,3),nint(denpar(6,1,j,i))
              endif
            endif
          endif
        endif
      enddo
    enddo
    if (.not.lout) then
      write(ioout,'(''--------------------------------------------------------------------------------'',/)')
    endif
!*******************
!  MEAM densities  *
!*******************
    lout = .true.
    do i = 1,neamspec
      nat1 = neamnat(i)
      ntype1 = neamtyp(i)
      call label(nat1,ntype1,lab1)
      cstyp1 = 'c'
      if (nat1.gt.maxele) cstyp1 = 's'
      nat2 = neamnat2(i)
      if (nat2.gt.0) then
        ntype2 = neamtyp2(i)
        call label(nat2,ntype2,lab2)
        cstyp2 = 'c'
        if (nat2.gt.maxele) cstyp1 = 's'
      endif
      do j = 1,ndenfncomp(i)
        if (neammeamorder(j,i).gt.1) then
          if (lout) then
            write(ioout,'(/,''  Modified Embedded Atom Model Densities :'',/)')
            write(ioout,'(''--------------------------------------------------------------------------------'')')
            write(ioout,'(''  Atom(s)       Functional/Order        A          B          C         n      '')')
            write(ioout,'(''--------------------------------------------------------------------------------'')')
            lout = .false.
          endif
          if (ndenfn(j,i).gt.0) then
            if (ndenfn(j,i).eq.1) then
              denfntyp = 'Power Law   '
            elseif (ndenfn(j,i).eq.2) then
              denfntyp = 'Exponential '
            elseif (ndenfn(j,i).eq.3) then
              denfntyp = 'Gaussian    '
            elseif (ndenfn(j,i).eq.4) then
              denfntyp = 'Cubic       '
            elseif (ndenfn(j,i).eq.5) then
              denfntyp = 'Voter-Chen  '
            elseif (ndenfn(j,i).eq.6) then
              denfntyp = 'Quadratic   '
            elseif (ndenfn(j,i).eq.7) then
              denfntyp = 'Quartic     '
            elseif (ndenfn(j,i).eq.8) then
              denfntyp = 'Glue        '
            elseif (ndenfn(j,i).eq.9) then
              denfntyp = 'eVoter-Chen '
            elseif (ndenfn(j,i).eq.10) then
              denfntyp = 'Mei-Davenprt'
            elseif (ndenfn(j,i).eq.12) then
              denfntyp = 'Baskes      '
            elseif (ndenfn(j,i).eq.14) then
              denfntyp = 'Fractl power'
            elseif (ndenfn(j,i).eq.15) then
              denfntyp = 'Spline      '
            elseif (ndenfn(j,i).eq.16) then
              denfntyp = 'Morse2      '
            endif
            if (ndenfn(j,i).eq.8) then
              if (nat2.gt.0) then
                write(ioout,'(2x,a4,1x,a1,1x,a4,1x,a1,1x,a12,1x,''0'',2x,g14.4,2(1x,g10.4))') &
                  lab1(1:4),cstyp1,lab2(1:4),cstyp2,denfntyp,(denpar(k,1,j,i),k=1,3)
              else
                write(ioout,'(2x,a4,1x,a1,8x,a12,1x,''0'',2x,g14.4,2(1x,g10.4))') &
                  lab1(1:4),cstyp1,denfntyp,(denpar(k,1,j,i),k=1,3)
              endif
              write(ioout,'(16x,4(f12.6,1x))') (denpar(k,1,j,i),k=4,7)
              write(ioout,'(16x,4(f12.6,1x))') (denpar(k,1,j,i),k=8,11)
              write(ioout,'(16x,4(f12.6,1x))') (denpar(k,1,j,i),k=12,15)
            elseif (ndenfn(j,i).eq.10) then
              if (nat2.gt.0) then
                write(ioout,'(2x,a4,1x,a1,1x,a4,1x,a1,1x,a12,1x,''0'',4x,3(f11.6,1x))') &
                  lab1(1:4),cstyp1,lab2(1:4),cstyp2,denfntyp,(denpar(k,1,j,i),k=1,3)
              else
                write(ioout,'(2x,a4,1x,a1,8x,a12,1x,''0'',4x,3(f11.6,1x))') &
                  lab1(1:4),cstyp1,denfntyp,(denpar(k,1,j,i),k=1,3)
              endif
              write(ioout,'(34x,3(f11.6,1x),f8.4)') (denpar(k,1,j,i),k=4,7)
            elseif (ndenfn(j,i).eq.15) then
              if (nat2.gt.0) then
                write(ioout,'(2x,a4,1x,a1,1x,a4,1x,a1,1x,a12,1x,''0'',4x,3(f11.6,1x))') &
                  lab1(1:4),cstyp1,lab2(1:4),cstyp2,denfntyp,(denpar(k,1,j,i),k=1,3)
              else
                write(ioout,'(2x,a4,1x,a1,8x,a12,1x,''0'',4x,3(f11.6,1x))') &
                  lab1(1:4),cstyp1,denfntyp,(denpar(k,1,j,i),k=1,3)
              endif
              write(ioout,'(34x,3(f11.6,1x))') (denpar(k,1,j,i),k=4,6)
            else
              if (nat2.gt.0) then
                write(ioout,'(2x,a4,1x,a1,1x,a4,1x,a1,1x,a12,1x,''0'',2x,g14.4,2(1x,g10.4),4x,i2)') &
                  lab1(1:4),cstyp1,lab2(1:4),cstyp2,denfntyp,(denpar(k,1,j,i),k=1,3),nint(denpar(6,1,j,i))
              else
                write(ioout,'(2x,a4,1x,a1,8x,a12,1x,''0'',2x,g14.4,2(1x,g10.4),4x,i2)') &
                  lab1(1:4),cstyp1,denfntyp,(denpar(k,1,j,i),k=1,3),nint(denpar(6,1,j,i))
              endif
            endif
            do order = 2,neammeamorder(j,i)
              if (ndenfn(j,i).eq.8) then
                write(ioout,'(29x,i1,2x,g14.4,2(1x,g10.4))') order-1,(denpar(k,order,j,i),k=1,3)
                write(ioout,'(16x,4(f12.6,1x))') (denpar(k,order,j,i),k=4,7)
                write(ioout,'(16x,4(f12.6,1x))') (denpar(k,order,j,i),k=8,11)
                write(ioout,'(16x,4(f12.6,1x))') (denpar(k,order,j,i),k=12,15)
              elseif (ndenfn(j,i).eq.10) then
                write(ioout,'(29x,i1,4x,3(f11.6,1x))') order-1,(denpar(k,order,j,i),k=1,3)
                write(ioout,'(34x,3(f11.6,1x),f8.4)') (denpar(k,order,j,i),k=4,7)
              else
                write(ioout,'(29x,i1,2x,g14.4,2(1x,g10.4),4x,i2)') order-1,(denpar(k,order,j,i),k=1,3),nint(denpar(6,order,j,i))
              endif
            enddo
          endif
        endif
      enddo
    enddo
    if (.not.lout) then
      write(ioout,'(''--------------------------------------------------------------------------------'',/)')
    endif
!**********************
!  EAM alloy factors  *
!**********************
    write(ioout,'(/,''  Embedded Atom Model alloy parameters :'',/)')
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    write(ioout,'(''  Atom(s)                          Scale factor '')')
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    do i = 1,neamspec
      nat1 = neamnat(i)
      ntype1 = neamtyp(i)
      call label(nat1,ntype1,lab1)
      cstyp1 = 'c'
      if (nat1.gt.maxele) cstyp1 = 's'   
      nat2 = neamnat2(i)
      if (nat2.gt.0) then
        ntype2 = neamtyp2(i)
        call label(nat2,ntype2,lab2)
        cstyp2 = 'c'
        if (nat2.gt.maxele) cstyp1 = 's'
        write(ioout,'(2x,a5,1x,a1,4x,a5,1x,a1,4x,f21.6)') lab1,cstyp1,lab2,cstyp2,eamalloy(1,i)
      else
        write(ioout,'(2x,a5,1x,a1,15x,f21.6)') lab1,cstyp1,eamalloy(1,i)
      endif
    enddo
    write(ioout,'(''--------------------------------------------------------------------------------'',/)')
!********************
!  EAM functionals  *
!********************
    if (neamfn.eq.1) then
      write(ioout,'(''  Embedded Atom Model functional = Square Root'',/)')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''  Atom                     A     '')')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      do i = 1,neamfnspec
        nat1 = neamfnnat(i)
        ntype1 = neamfntyp(i)
        call label(nat1,ntype1,lab1)
        cstyp1 = 'c'
        if (nat1.gt.maxele) cstyp1 = 's'
        write(ioout,'(2x,a5,1x,a1,9x,f14.6)') lab1,cstyp1,eamfnpar(1,i)
        if (neamfnmeamorder(i).gt.1) then
          if (nmeamcombotype.eq.1.and.nmeamrhotype.eq.1) then
            write(ioout,'(10x,''MEAM type = sum : 21 term'')')
          elseif (nmeamcombotype.eq.1.and.nmeamrhotype.eq.2) then
            write(ioout,'(10x,''MEAM type = sum : 24 term'')')
          elseif (nmeamcombotype.eq.2.and.nmeamrhotype.eq.1) then
            write(ioout,'(10x,''MEAM type = exponential : 21 term'')')
          elseif (nmeamcombotype.eq.3.and.nmeamrhotype.eq.1) then
            write(ioout,'(10x,''MEAM type = exphalf : 21 term'')')
          elseif (nmeamcombotype.eq.3.and.nmeamrhotype.eq.2) then
            write(ioout,'(10x,''MEAM type = exphalf : 24 term'')')
          else
            write(ioout,'(10x,''MEAM type = exponential : 24 term'')')
          endif
          write(ioout,'(10x,''MEAM order/coefficients = '',i2,4(1x,f8.4))') &
            neamfnmeamorder(i),(eamfnmeamcoeff(j,i),j=1,neamfnmeamorder(i))
        endif
      enddo
      write(ioout,'(''--------------------------------------------------------------------------------'',/)')
    elseif (neamfn.eq.2) then
      write(ioout,'(''  Embedded Atom Model functional = Inverse Power Law of '',i2,/)') neampower
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''  Atom                     A     '')')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      do i = 1,neamfnspec
        nat1 = neamfnnat(i)
        ntype1 = neamfntyp(i)
        call label(nat1,ntype1,lab1)
        cstyp1 = 'c'
        if (nat1.gt.maxele) cstyp1 = 's'
        write(ioout,'(2x,a5,1x,a1,9x,f14.6)') lab1,cstyp1,eamfnpar(1,i)
        if (neamfnmeamorder(i).gt.1) then
          if (nmeamcombotype.eq.1.and.nmeamrhotype.eq.1) then
            write(ioout,'(10x,''MEAM type = sum : 21 term'')')
          elseif (nmeamcombotype.eq.1.and.nmeamrhotype.eq.2) then
            write(ioout,'(10x,''MEAM type = sum : 24 term'')')
          elseif (nmeamcombotype.eq.2.and.nmeamrhotype.eq.1) then
            write(ioout,'(10x,''MEAM type = exponential : 21 term'')')
          elseif (nmeamcombotype.eq.3.and.nmeamrhotype.eq.1) then
            write(ioout,'(10x,''MEAM type = exphalf : 21 term'')')
          elseif (nmeamcombotype.eq.3.and.nmeamrhotype.eq.2) then
            write(ioout,'(10x,''MEAM type = exphalf : 24 term'')')
          else
            write(ioout,'(10x,''MEAM type = exponential : 24 term'')')
          endif
          write(ioout,'(10x,''MEAM order/coefficients = '',i2,4(1x,f8.4))') &
            neamfnmeamorder(i),(eamfnmeamcoeff(j,i),j=1,neamfnmeamorder(i))
        endif
      enddo
      write(ioout,'(''--------------------------------------------------------------------------------'',/)')
    elseif (neamfn.eq.3) then
      write(ioout,'(''  Embedded Atom Model functional = Banerjea/Smith : Power = '',i2,/)') neampower
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''  Atom                     F0         F1        rho0               '')')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      do i = 1,neamfnspec
        nat1 = neamfnnat(i)
        ntype1 = neamfntyp(i)
        call label(nat1,ntype1,lab1)
        cstyp1 = 'c'
        if (nat1.gt.maxele) cstyp1 = 's'
        write(ioout,'(2x,a5,1x,a1,9x,g14.4,2(1x,g10.4))') lab1,cstyp1,(eamfnpar(j,i),j=1,3)
        if (neamfnmeamorder(i).gt.1) then
          if (nmeamcombotype.eq.1.and.nmeamrhotype.eq.1) then
            write(ioout,'(10x,''MEAM type = sum : 21 term'')')
          elseif (nmeamcombotype.eq.1.and.nmeamrhotype.eq.2) then
            write(ioout,'(10x,''MEAM type = sum : 24 term'')')
          elseif (nmeamcombotype.eq.2.and.nmeamrhotype.eq.1) then
            write(ioout,'(10x,''MEAM type = exponential : 21 term'')')
          elseif (nmeamcombotype.eq.3.and.nmeamrhotype.eq.1) then
            write(ioout,'(10x,''MEAM type = exphalf : 21 term'')')
          elseif (nmeamcombotype.eq.3.and.nmeamrhotype.eq.2) then
            write(ioout,'(10x,''MEAM type = exphalf : 24 term'')')
          else
            write(ioout,'(10x,''MEAM type = exponential : 24 term'')')
          endif
          write(ioout,'(10x,''MEAM order/coefficients = '',i2,4(1x,f8.4))') &
            neamfnmeamorder(i),(eamfnmeamcoeff(j,i),j=1,neamfnmeamorder(i))
        endif
      enddo
      write(ioout,'(''--------------------------------------------------------------------------------'',/)')
    elseif (neamfn.eq.4) then
      write(ioout,'(''  Embedded Atom Model functional = Numerical '',/)') 
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''  Atom    Filename '')')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      do i = 1,neamfnspec
        nat1 = neamfnnat(i)
        ntype1 = neamfntyp(i)
        call label(nat1,ntype1,lab1)
        cstyp1 = 'c'
        if (nat1.gt.maxele) cstyp1 = 's'
        write(ioout,'(2x,a5,1x,a1,1x,a70)') lab1,cstyp1,eamfnfile(i)(1:70)
        if (neamfnmeamorder(i).gt.1) then
          if (nmeamcombotype.eq.1.and.nmeamrhotype.eq.1) then
            write(ioout,'(10x,''MEAM type = sum : 21 term'')')
          elseif (nmeamcombotype.eq.1.and.nmeamrhotype.eq.2) then
            write(ioout,'(10x,''MEAM type = sum : 24 term'')')
          elseif (nmeamcombotype.eq.2.and.nmeamrhotype.eq.1) then
            write(ioout,'(10x,''MEAM type = exponential : 21 term'')')
          elseif (nmeamcombotype.eq.3.and.nmeamrhotype.eq.1) then
            write(ioout,'(10x,''MEAM type = exphalf : 21 term'')')
          elseif (nmeamcombotype.eq.3.and.nmeamrhotype.eq.2) then
            write(ioout,'(10x,''MEAM type = exphalf : 24 term'')')
          else
            write(ioout,'(10x,''MEAM type = exponential : 24 term'')')
          endif
          write(ioout,'(10x,''MEAM order/coefficients = '',i2,4(1x,f8.4))') &
            neamfnmeamorder(i),(eamfnmeamcoeff(j,i),j=1,neamfnmeamorder(i))
        endif
      enddo
      write(ioout,'(''--------------------------------------------------------------------------------'',/)')
    elseif (neamfn.eq.5) then
      write(ioout,'(''  Embedded Atom Model functional = Johnson '',/)') 
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''  Atom                  F0         F1        rho0       Alpha    Beta    Gamma  '')')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      do i = 1,neamfnspec
        nat1 = neamfnnat(i)
        ntype1 = neamfntyp(i)
        call label(nat1,ntype1,lab1)
        cstyp1 = 'c'
        if (nat1.gt.maxele) cstyp1 = 's'
        write(ioout,'(2x,a5,1x,a1,9x,g14.4,2(1x,g10.4),3(1x,f7.5))') lab1,cstyp1,(eamfnpar(j,i),j=1,6)
        if (neamfnmeamorder(i).gt.1) then
          if (nmeamcombotype.eq.1.and.nmeamrhotype.eq.1) then
            write(ioout,'(10x,''MEAM type = sum : 21 term'')')
          elseif (nmeamcombotype.eq.1.and.nmeamrhotype.eq.2) then
            write(ioout,'(10x,''MEAM type = sum : 24 term'')')
          elseif (nmeamcombotype.eq.2.and.nmeamrhotype.eq.1) then
            write(ioout,'(10x,''MEAM type = exponential : 21 term'')')
          elseif (nmeamcombotype.eq.3.and.nmeamrhotype.eq.1) then
            write(ioout,'(10x,''MEAM type = exphalf : 21 term'')')
          elseif (nmeamcombotype.eq.3.and.nmeamrhotype.eq.2) then
            write(ioout,'(10x,''MEAM type = exphalf : 24 term'')')
          else
            write(ioout,'(10x,''MEAM type = exponential : 24 term'')')
          endif
          write(ioout,'(10x,''MEAM order/coefficients = '',i2,4(1x,f8.4))') &
            neamfnmeamorder(i),(eamfnmeamcoeff(j,i),j=1,neamfnmeamorder(i))
        endif
      enddo
      write(ioout,'(''--------------------------------------------------------------------------------'',/)')
    elseif (neamfn.eq.6) then
      write(ioout,'(''  Embedded Atom Model functional = Glue '',/)') 
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''  Atom                     rho1         rho2        '')')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      do i = 1,neamfnspec
        nat1 = neamfnnat(i)
        ntype1 = neamfntyp(i)
        call label(nat1,ntype1,lab1)
        cstyp1 = 'c'
        if (nat1.gt.maxele) cstyp1 = 's'
        write(ioout,'(2x,a5,1x,a1,12x,2(1x,f12.6))') lab1,cstyp1,(eamfnpar(j,i),j=1,2)
        write(ioout,'(2x,5(f12.8,1x))') (eamfnpar(j,i),j=3,7)
        write(ioout,'(2x,5(f12.8,1x))') (eamfnpar(j,i),j=8,12)
        write(ioout,'(2x,4(f12.8,1x))') (eamfnpar(j,i),j=13,16)
        if (neamfnmeamorder(i).gt.1) then
          if (nmeamcombotype.eq.1.and.nmeamrhotype.eq.1) then
            write(ioout,'(10x,''MEAM type = sum : 21 term'')')
          elseif (nmeamcombotype.eq.1.and.nmeamrhotype.eq.2) then
            write(ioout,'(10x,''MEAM type = sum : 24 term'')')
          elseif (nmeamcombotype.eq.2.and.nmeamrhotype.eq.1) then
            write(ioout,'(10x,''MEAM type = exponential : 21 term'')')
          elseif (nmeamcombotype.eq.3.and.nmeamrhotype.eq.1) then
            write(ioout,'(10x,''MEAM type = exphalf : 21 term'')')
          elseif (nmeamcombotype.eq.3.and.nmeamrhotype.eq.2) then
            write(ioout,'(10x,''MEAM type = exphalf : 24 term'')')
          else
            write(ioout,'(10x,''MEAM type = exponential : 24 term'')')
          endif
          write(ioout,'(10x,''MEAM order/coefficients = '',i2,4(1x,f8.4))') &
            neamfnmeamorder(i),(eamfnmeamcoeff(j,i),j=1,neamfnmeamorder(i))
        endif
      enddo
      write(ioout,'(''--------------------------------------------------------------------------------'',/)')
    elseif (neamfn.eq.7) then
      write(ioout,'(''  Embedded Atom Model functional = Foiles '',/)') 
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''  Atom                  F0         F1        F2       F3                   '')')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      do i = 1,neamfnspec
        nat1 = neamfnnat(i)
        ntype1 = neamfntyp(i)
        call label(nat1,ntype1,lab1)
        cstyp1 = 'c'
        if (nat1.gt.maxele) cstyp1 = 's'
        write(ioout,'(2x,a5,1x,a1,9x,g14.4,2(1x,g10.4),1x,f7.5)') lab1,cstyp1,(eamfnpar(j,i),j=1,4)
        if (neamfnmeamorder(i).gt.1) then
          if (nmeamcombotype.eq.1.and.nmeamrhotype.eq.1) then
            write(ioout,'(10x,''MEAM type = sum : 21 term'')')
          elseif (nmeamcombotype.eq.1.and.nmeamrhotype.eq.2) then
            write(ioout,'(10x,''MEAM type = sum : 24 term'')')
          elseif (nmeamcombotype.eq.2.and.nmeamrhotype.eq.1) then
            write(ioout,'(10x,''MEAM type = exponential : 21 term'')')
          elseif (nmeamcombotype.eq.3.and.nmeamrhotype.eq.1) then
            write(ioout,'(10x,''MEAM type = exphalf : 21 term'')')
          elseif (nmeamcombotype.eq.3.and.nmeamrhotype.eq.2) then
            write(ioout,'(10x,''MEAM type = exphalf : 24 term'')')
          else
            write(ioout,'(10x,''MEAM type = exponential : 24 term'')')
          endif
          write(ioout,'(10x,''MEAM order/coefficients = '',i2,4(1x,f8.4))') &
            neamfnmeamorder(i),(eamfnmeamcoeff(j,i),j=1,neamfnmeamorder(i))
        endif
      enddo
      write(ioout,'(''--------------------------------------------------------------------------------'',/)')
    elseif (neamfn.eq.8) then
      write(ioout,'(''  Embedded Atom Model functional = Mei-Davenport '',/)') 
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''  Atom                 Ec          Alpha      Beta     Gamma     Delta       '')')
      write(ioout,'(''                      (eV)          phi0       s1       s2        s3         '')')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      do i = 1,neamfnspec
        nat1 = neamfnnat(i)
        ntype1 = neamfntyp(i)
        call label(nat1,ntype1,lab1)
        cstyp1 = 'c'
        if (nat1.gt.maxele) cstyp1 = 's'
        write(ioout,'(2x,a5,1x,a1,6x,f14.6,2x,4(1x,f9.5))') lab1,cstyp1,(eamfnpar(j,i),j=1,5)
        write(ioout,'(31x,4(1x,f9.5))') (eamfnpar(j,i),j=6,9)
        if (neamfnmeamorder(i).gt.1) then
          if (nmeamcombotype.eq.1.and.nmeamrhotype.eq.1) then
            write(ioout,'(10x,''MEAM type = sum : 21 term'')')
          elseif (nmeamcombotype.eq.1.and.nmeamrhotype.eq.2) then
            write(ioout,'(10x,''MEAM type = sum : 24 term'')')
          elseif (nmeamcombotype.eq.2.and.nmeamrhotype.eq.1) then
            write(ioout,'(10x,''MEAM type = exponential : 21 term'')')
          elseif (nmeamcombotype.eq.3.and.nmeamrhotype.eq.1) then
            write(ioout,'(10x,''MEAM type = exphalf : 21 term'')')
          elseif (nmeamcombotype.eq.3.and.nmeamrhotype.eq.2) then
            write(ioout,'(10x,''MEAM type = exphalf : 24 term'')')
          else
            write(ioout,'(10x,''MEAM type = exponential : 24 term'')')
          endif
          write(ioout,'(10x,''MEAM order/coefficients = '',i2,4(1x,f8.4))') &
            neamfnmeamorder(i),(eamfnmeamcoeff(j,i),j=1,neamfnmeamorder(i))
        endif
      enddo
      write(ioout,'(''--------------------------------------------------------------------------------'',/)')
    elseif (neamfn.eq.9) then
      write(ioout,'(''  Embedded Atom Model functional = Baskes'',/)')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''  Atom                     Ec          A            rho0'')')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      do i = 1,neamfnspec
        nat1 = neamfnnat(i)
        ntype1 = neamfntyp(i)
        call label(nat1,ntype1,lab1)
        cstyp1 = 'c'
        if (nat1.gt.maxele) cstyp1 = 's'
        write(ioout,'(2x,a5,1x,a1,9x,f14.6,2(1x,g14.6))') lab1,cstyp1,eamfnpar(1,i),eamfnpar(2,i),eamfnpar(3,i)
        if (neamfnmeamorder(i).gt.1) then
          if (nmeamcombotype.eq.1.and.nmeamrhotype.eq.1) then
            write(ioout,'(10x,''MEAM type = sum : 21 term'')')
          elseif (nmeamcombotype.eq.1.and.nmeamrhotype.eq.2) then
            write(ioout,'(10x,''MEAM type = sum : 24 term'')')
          elseif (nmeamcombotype.eq.2.and.nmeamrhotype.eq.1) then
            write(ioout,'(10x,''MEAM type = exponential : 21 term'')')
          elseif (nmeamcombotype.eq.3.and.nmeamrhotype.eq.1) then
            write(ioout,'(10x,''MEAM type = exphalf : 21 term'')')
          elseif (nmeamcombotype.eq.3.and.nmeamrhotype.eq.2) then
            write(ioout,'(10x,''MEAM type = exphalf : 24 term'')')
          else
            write(ioout,'(10x,''MEAM type = exponential : 24 term'')')
          endif
          write(ioout,'(10x,''MEAM order/coefficients = '',i2,4(1x,f8.4))') &
            neamfnmeamorder(i),(eamfnmeamcoeff(j,i),j=1,neamfnmeamorder(i))
        endif
      enddo
      write(ioout,'(''--------------------------------------------------------------------------------'',/)')
    elseif (neamfn.eq.10) then
      write(ioout,'(''  Embedded Atom Model functional = VBO'',/)')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''  Atom                     A         rn'')')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      do i = 1,neamfnspec
        nat1 = neamfnnat(i)
        ntype1 = neamfntyp(i)
        call label(nat1,ntype1,lab1)
        cstyp1 = 'c'
        if (nat1.gt.maxele) cstyp1 = 's'
        write(ioout,'(2x,a5,1x,a1,9x,f14.6,1x,g10.4)') lab1,cstyp1,eamfnpar(1,i),eamfnpar(2,i)
        if (neamfnmeamorder(i).gt.1) then
          if (nmeamcombotype.eq.1.and.nmeamrhotype.eq.1) then
            write(ioout,'(10x,''MEAM type = sum : 21 term'')')
          elseif (nmeamcombotype.eq.1.and.nmeamrhotype.eq.2) then
            write(ioout,'(10x,''MEAM type = sum : 24 term'')')
          elseif (nmeamcombotype.eq.2.and.nmeamrhotype.eq.1) then
            write(ioout,'(10x,''MEAM type = exponential : 21 term'')')
          elseif (nmeamcombotype.eq.3.and.nmeamrhotype.eq.1) then
            write(ioout,'(10x,''MEAM type = exphalf : 21 term'')')
          elseif (nmeamcombotype.eq.3.and.nmeamrhotype.eq.2) then
            write(ioout,'(10x,''MEAM type = exphalf : 24 term'')')
          else
            write(ioout,'(10x,''MEAM type = exponential : 24 term'')')
          endif
          write(ioout,'(10x,''MEAM order/coefficients = '',i2,4(1x,f8.4))') &
            neamfnmeamorder(i),(eamfnmeamcoeff(j,i),j=1,neamfnmeamorder(i))
        endif
      enddo
      write(ioout,'(''--------------------------------------------------------------------------------'',/)')
    elseif (neamfn.eq.11) then
      write(ioout,'(''  Embedded Atom Model functional = Belashchenko '',/)')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''  Atom      Rho1         Rho2          Rho3        '')')
      write(ioout,'(''            Rho4         Rho5          Rho6         Rho7 '')')
      write(ioout,'(''            a1           m             n    '')')
      write(ioout,'(''            c1           c2            c3           c4   '')')
      write(ioout,'(''            c5           c6            c7           c8   '')')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      do i = 1,neamfnspec
        nat1 = neamfnnat(i)
        ntype1 = neamfntyp(i)
        call label(nat1,ntype1,lab1)
        cstyp1 = 'c'
        if (nat1.gt.maxele) cstyp1 = 's'
        write(ioout,'(2x,a5,1x,a1,3(1x,f12.6))') lab1,cstyp1,(eamfnpar(j,i),j=1,3)
        write(ioout,'(9x,4(1x,f12.6))') (eamfnpar(j,i),j=4,7)
        write(ioout,'(10x,3(f12.8,1x))') eamfnpar(8,i),eamfnpar(32,i),eamfnpar(33,i)
        write(ioout,'(10x,4(f12.8,1x))') (eamfnpar(j,i),j=24,27)
        write(ioout,'(10x,4(f12.8,1x))') (eamfnpar(j,i),j=28,31)
      enddo
      write(ioout,'(''--------------------------------------------------------------------------------'',/)')
    elseif (neamfn.eq.12) then
      write(ioout,'(''  Embedded Atom Model functional = Igarashi'',/)')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''  Atom                     A         B'')')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      do i = 1,neamfnspec
        nat1 = neamfnnat(i)
        ntype1 = neamfntyp(i)
        call label(nat1,ntype1,lab1)
        cstyp1 = 'c'
        if (nat1.gt.maxele) cstyp1 = 's'
        write(ioout,'(2x,a5,1x,a1,9x,f14.6,1x,g14.6)') lab1,cstyp1,eamfnpar(1,i),eamfnpar(2,i)
        if (neamfnmeamorder(i).gt.1) then
          if (nmeamcombotype.eq.1.and.nmeamrhotype.eq.1) then
            write(ioout,'(10x,''MEAM type = sum : 21 term'')')
          elseif (nmeamcombotype.eq.1.and.nmeamrhotype.eq.2) then
            write(ioout,'(10x,''MEAM type = sum : 24 term'')')
          elseif (nmeamcombotype.eq.2.and.nmeamrhotype.eq.1) then
            write(ioout,'(10x,''MEAM type = exponential : 21 term'')')
          elseif (nmeamcombotype.eq.3.and.nmeamrhotype.eq.1) then
            write(ioout,'(10x,''MEAM type = exphalf : 21 term'')')
          elseif (nmeamcombotype.eq.3.and.nmeamrhotype.eq.2) then
            write(ioout,'(10x,''MEAM type = exphalf : 24 term'')')
          else
            write(ioout,'(10x,''MEAM type = exponential : 24 term'')')
          endif
          write(ioout,'(10x,''MEAM order/coefficients = '',i2,4(1x,f8.4))') &
            neamfnmeamorder(i),(eamfnmeamcoeff(j,i),j=1,neamfnmeamorder(i))
        endif
      enddo
      write(ioout,'(''--------------------------------------------------------------------------------'',/)')
    elseif (neamfn.eq.13) then
      write(ioout,'(''  Embedded Atom Model functional = Spline '',/)')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''  Atom           A          B          C          D          rho0     rho_max '')')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      do i = 1,neamfnspec
        nat1 = neamfnnat(i)
        ntype1 = neamfntyp(i)
        call label(nat1,ntype1,lab1)
        cstyp1 = 'c'
        if (nat1.gt.maxele) cstyp1 = 's'
        write(ioout,'(2x,a5,1x,a1,5x,4(1x,g10.4),2(1x,f9.4))') lab1,cstyp1,(eamfnpar(j,i),j=1,6)
        if (neamfnmeamorder(i).gt.1) then
          if (nmeamcombotype.eq.1.and.nmeamrhotype.eq.1) then
            write(ioout,'(10x,''MEAM type = sum : 21 term'')')
          elseif (nmeamcombotype.eq.1.and.nmeamrhotype.eq.2) then
            write(ioout,'(10x,''MEAM type = sum : 24 term'')')
          elseif (nmeamcombotype.eq.2.and.nmeamrhotype.eq.1) then
            write(ioout,'(10x,''MEAM type = exponential : 21 term'')')
          elseif (nmeamcombotype.eq.3.and.nmeamrhotype.eq.1) then
            write(ioout,'(10x,''MEAM type = exphalf : 21 term'')')
          elseif (nmeamcombotype.eq.3.and.nmeamrhotype.eq.2) then
            write(ioout,'(10x,''MEAM type = exphalf : 24 term'')')
          else
            write(ioout,'(10x,''MEAM type = exponential : 24 term'')')
          endif
          write(ioout,'(10x,''MEAM order/coefficients = '',i2,4(1x,f8.4))') &
            neamfnmeamorder(i),(eamfnmeamcoeff(j,i),j=1,neamfnmeamorder(i))
        endif
      enddo
      write(ioout,'(''--------------------------------------------------------------------------------'',/)')
    endif
!**********************
!  Cross term cutoff  *
!**********************
    if (eamXcutfactor.ne.1.0_dp) then
      write(ioout,'(/,''  Manybody second derivative cross term cutoffs to be scaled by = '',f10.8)') eamXcutfactor
    endif
!*******************
!  MEAM screening  *
!*******************
    lout = .true.
    do i = 1,neamfnspec
      ind = 0
      do j = 1,neamfnspec
        do k = 1,j
          ind = ind + 1
          if (lMEAMscreen(ind,i)) then
            if (lout) then
              write(ioout,'(''  Modified Embedded Atom Model screening parameters '',/)')
              write(ioout,'(''--------------------------------------------------------------------------------'')')
              write(ioout,'(''  Atom 1    Atom 2    Atom 3            C_min       C_max '')')
              write(ioout,'(''--------------------------------------------------------------------------------'')')
              lout = .false. 
            endif
            nat1 = neamfnnat(i)
            ntype1 = neamfntyp(i)
            call label(nat1,ntype1,lab1)
            cstyp1 = 'c'
            if (nat1.gt.maxele) cstyp1 = 's'
!
            nat2 = neamfnnat(j)
            ntype2 = neamfntyp(j)
            call label(nat2,ntype2,lab2)
            cstyp2 = 'c'
            if (nat2.gt.maxele) cstyp2 = 's'
!
            nat3 = neamfnnat(k)
            ntype3 = neamfntyp(k)
            call label(nat3,ntype3,lab3)
            cstyp3 = 'c'
            if (nat3.gt.maxele) cstyp3 = 's'
!
            write(ioout,'(2x,3(a5,1x,a1,3x),3x,2f12.6)') lab1,cstyp1,lab2,cstyp2,lab3,cstyp3, &
                  meam_Cmin(ind,i),meam_Cmax(ind,i)
          endif
        enddo
      enddo
    enddo
    if (.not.lout) then
      write(ioout,'(''--------------------------------------------------------------------------------'',/)')
    endif
  endif
!**************************
!  Shift operator output  *
!**************************
  if (nshift.gt.0) then
    write(ioout,'(/,''  Energy shifts for configurations :'',/)')
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    write(ioout,'(''      Configuration           Energy (eV)        Scale Factor'')')
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    do i=1,ncfg
      write(ioout,'(11x,i3,15x,g12.6,10x,f8.3)')i,shift(nshcfg(i)),shscalecfg(i)
    enddo
    write(ioout,'(''--------------------------------------------------------------------------------'',/)')
  endif
!***********************************
!  Output three-bodied potentials  *
!***********************************
  if (nthb.gt.0) then
!
!  Three-body
!
    call outthree
  endif
!**********************************
!  Output four-bodied potentials  *
!**********************************
  if (nfor.gt.0) then
!
!  Four-body
!
    call outfour
  endif
!**********************************
!  Output six-bodied potentials  *
!**********************************
  if (nsix.gt.0) then
!
!  Count number of general, intra- and inter-molecular potentials
!
    nmpt(1) = 0
    nmpt(2) = 0
    nmpt(3) = 0
    allocate(iptyp(nsix),stat=status)
    if (status/=0) call outofmemory('outpot','iptyp')
    if (lmol) then
      do i = 1,nsix
        if (lsintra(i).and..not.lsinter(i)) then
          nmpt(2) = nmpt(2) + 1
          iptyp(i) = 2
        elseif (lsinter(i).and..not.lsintra(i)) then
          nmpt(3) = nmpt(3) + 1
          iptyp(i) = 3
        else
          nmpt(1) = nmpt(1) + 1
          iptyp(i) = 1
        endif
      enddo
    else
      nmpt(1) = nsix
      do i = 1,nsix
        iptyp(i) = 1
      enddo
    endif
!
!  Standard six-body potentials
!
    do k = 1,3
      if (nmpt(k).gt.0) then
        if (k.eq.1) then
          write(ioout,'(/,''  General Six-body potentials :'')')
        elseif (k.eq.2) then
          write(ioout,'(/,''  Intramolecular Six-body potentials :'')')
        elseif (k.eq.3) then
          write(ioout,'(/,''  Intermolecular Six-body potentials :'')')
        endif
        lout = .true.
        do i = 1,nsix
          if (nsixty(i).eq.1.and.iptyp(i).eq.k) then
            if (lout) then
              lout = .false.
              write(ioout,'(/,''  Cross - out of plane :'',/)')
              write(ioout,'(''--------------------------------------------------------------------------------'')')
              write(ioout,'(''          Atom Types              Force cst                       Cutoffs       '')')
              write(ioout,'(''           1       2                 (eV)                           1-2         '')')
              write(ioout,'(''   3       4       5       6                                1-3  1-4  2-5  2-6  '')')
              write(ioout,'(''--------------------------------------------------------------------------------'')')
            endif
            nat1 = nsspec1(i)
            nat2 = nsspec2(i)
            nat3 = nsspec3(i)
            nat4 = nsspec4(i)
            nat5 = nsspec5(i)
            nat6 = nsspec6(i)
            ntype1 = nsptyp1(i)
            ntype2 = nsptyp2(i)
            ntype3 = nsptyp3(i)
            ntype4 = nsptyp4(i)
            ntype5 = nsptyp5(i)
            ntype6 = nsptyp6(i)
            call label(nat1,ntype1,lab1)
            call label(nat2,ntype2,lab2)
            call label(nat3,ntype3,lab3)
            call label(nat4,ntype4,lab4)
            call label(nat5,ntype5,lab5)
            call label(nat6,ntype6,lab6)
            cstyp1 = 'c'
            cstyp2 = 'c'
            cstyp3 = 'c'
            cstyp4 = 'c'
            cstyp5 = 'c'
            cstyp6 = 'c'
            if (nat1.gt.maxele) cstyp1 = 's'
            if (nat2.gt.maxele) cstyp2 = 's'
            if (nat3.gt.maxele) cstyp3 = 's'
            if (nat4.gt.maxele) cstyp4 = 's'
            if (nat5.gt.maxele) cstyp5 = 's'
            if (nat6.gt.maxele) cstyp6 = 's'
            if (mmsexc(i).eq.1) then
              nbotyptr = n6botype(1,i) + 1
              write(ioout,'(8x,2(a5,1x,a1,1x),9x,g12.5,15x,''Bonded'',1x,a9)') &
                lab1,cstyp1,lab2,cstyp2,sixk(i),botyword(nbotyptr)
              write(ioout,'(4(a5,1x,a1,1x))') &
                lab3,cstyp3,lab4,cstyp4,lab5,cstyp5,lab6,cstyp6
            else
              write(ioout,'(8x,2(a5,1x,a1,1x),9x,g12.5,21x,f5.2)') &
                lab1,cstyp1,lab2,cstyp2,sixk(i),six1(i)
              write(ioout,'(4(a5,1x,a1,1x),26x,4f5.2)') &
                lab3,cstyp3,lab4,cstyp4,lab5,cstyp5,lab6,cstyp6,six2(i),six3(i),six4(i),six5(i)
            endif
          endif
        enddo
        if (.not.lout) then
          write(ioout,'(''--------------------------------------------------------------------------------'')')
        endif
      endif
    enddo
    deallocate(iptyp,stat=status)
    if (status/=0) call deallocate_error('outpot','iptyp')
  endif
!***********************************
!  Bond-order one-body potentials  *
!***********************************
  if (nboO.gt.0) then
    lout = .true.
    do i = 1,nboO
      nat1 = nBOspec0(i)
      ntype1 = nBOtyp0(i)
      call label(nat1,ntype1,lab1)
      cstyp1 = 'c'
      if (nat1.gt.maxele) cstyp1 = 's'
!
      if (lout) then
        lout = .false.
        write(ioout,'(/,''  Bond-order one-body potentials :'',/)')
        write(ioout,'(''--------------------------------------------------------------------------------'')')
        write(ioout,'(''Atom  Type            AlphaA               AlphaB          nA       nB          '')')
        write(ioout,'(''--------------------------------------------------------------------------------'')')
      endif
      write(ioout,'(a5,a1,1x,4x,2(f20.10,1x),1x,2(1x,f8.4))') lab1,cstyp1,BOecoeffR(i),BOecoeffA(i),BOncoeffR(i),BOncoeffA(i)
      if (.not.lout) then
        write(ioout,'(''--------------------------------------------------------------------------------'')')
      endif
    enddo
  endif
!***********************************
!  Bond-order two-body potentials  *
!***********************************
  if (nbopot.gt.0) then
    lout = .true.
    do i = 1,nbopot
      nat1 = nBOspec1(i)
      nat2 = nBOspec2(i)
      ntype1 = nBOtyp1(i)
      ntype2 = nBOtyp2(i)
      call label(nat1,ntype1,lab1)
      call label(nat2,ntype2,lab2)
      cstyp1 = 'c'
      cstyp2 = 'c'
      if (nat1.gt.maxele) cstyp1 = 's'
      if (nat2.gt.maxele) cstyp2 = 's'
!
      if (lout) then
        lout = .false.
        write(ioout,'(/,''  Bond-order two-body potentials :'',/)')
        write(ioout,'(''--------------------------------------------------------------------------------'')')
        write(ioout,'(''Atom  Types                 A                 B                  Cutoffs(Ang)   '')')
        write(ioout,'(''  1     2                 ZetaA             ZetaB                Taper    Max   '')')
        write(ioout,'(''--------------------------------------------------------------------------------'')')
      endif
      write(ioout,'(2(a5,a1,1x),4x,2(f20.10,1x),1x,2(1x,f8.4))') lab1,cstyp1,lab2,cstyp2,BOacoeff(i), &
        BObcoeff(i),rBOmin(i),rBOmax(i)
      write(ioout,'(21x,2(g20.10,1x))') BOzacoeff(i),BOzbcoeff(i)
      if (.not.lout) then
        write(ioout,'(''--------------------------------------------------------------------------------'')')
      endif
    enddo
  endif
!
  if (nboR.gt.0) then
    lout = .true.
    do i = 1,nboR
      nat1 = nBOspecR1(i)
      ntype1 = nBOtypR1(i)
      call label(nat1,ntype1,lab1)
      cstyp1 = 'c'
      if (nat1.gt.maxele) cstyp1 = 's'
      nat2 = nBOspecR2(i)
      ntype2 = nBOtypR2(i)
      call label(nat2,ntype2,lab2)
      cstyp2 = 'c'
      if (nat2.gt.maxele) cstyp2 = 's'
!
      if (lout) then
        lout = .false.
        write(ioout,'(/,''  Bond-order repulsive :'',/)')
        write(ioout,'(''--------------------------------------------------------------------------------'')')
        write(ioout,'(''Atom    Type         Omega            Lambda          M                           '')')
        if (nBOtypeR(i).eq.2) then
          write(ioout,'('' 1/2     1/2           C              D                H                        '')')
        elseif (nBOtypeR(i).eq.3) then
          write(ioout,'('' 1/2     1/2           C1             C2               H                        '')')
          write(ioout,'(''                       C3             C4               C5                       '')')
        elseif (nBOtypeR(i).eq.4) then
          write(ioout,'('' 1/2     1/2           C0             C1               C2                       '')')
          write(ioout,'(''                       C3             C4                                        '')')
        else
          write(ioout,'('' 1/2     1/2                                                                    '')')
        endif
        write(ioout,'(''--------------------------------------------------------------------------------'')')
      endif
      if (nBOtypeR(i).eq.1) then
        write(ioout,'(a5,a1,6x,1x,2(f16.8,1x),5x,i3)') lab1,cstyp1,BOocoeffR(i),BOlcoeffR(i), &
          nint(BOmcoeffR(i))
        write(ioout,'(a5,a1)') lab2,cstyp2
      else
        if (nBOtypeR(i).eq.3) then
          write(ioout,'(a5,a1,2x,''Kuma-'',1x,2(f16.8,1x),4x,i3)') lab1,cstyp1,BOocoeffR(i), &
            BOlcoeffR(i),nint(BOmcoeffR(i))
          write(ioout,'(a5,a1,2x,''-gai '',1x,3(f16.8,1x))') lab2,cstyp2,BOccoeffR(1,i),BOccoeffR(2,i), &
            BOhcoeffR(i)
          write(ioout,'(14x,3(f16.8,1x))') BOccoeffR(3,i),BOccoeffR(4,i),BOccoeffR(5,i)
        elseif (nBOtypeR(i).eq.4) then
          write(ioout,'(a5,a1,2x,''MMP  '',1x,3(f16.8,1x),4x,i3)') lab1,cstyp1,BOecoeffR(i),BOncoeffR(i), &
            BOlcoeffR(i),nint(BOmcoeffR(i))
          write(ioout,'(a5,a1,8x,3(f16.8,1x))') lab2,cstyp2,BOccoeffR(1,i),BOccoeffR(2,i), &
            BOccoeffR(3,i)
          write(ioout,'(14x,2(f16.8,1x))') BOccoeffR(4,i),BOccoeffR(5,i)
        else
          write(ioout,'(a5,a1,2x,''Theta'',1x,2(f16.8,1x),4x,i3)') lab1,cstyp1,BOocoeffR(i), &
            BOlcoeffR(i),nint(BOmcoeffR(i))
          write(ioout,'(a5,a1,8x,3(f16.8,1x))') lab2,cstyp2,BOccoeffR(1,i),BOccoeffR(2,i),BOhcoeffR(i)
        endif
      endif
      if (.not.lout) then
        write(ioout,'(''--------------------------------------------------------------------------------'')')
      endif
    enddo
  endif
!
  if (nboA.gt.0) then
    lout = .true.
    do i = 1,nboA
      nat1 = nBOspecA1(i)      
      ntype1 = nBOtypA1(i)      
      call label(nat1,ntype1,lab1)      
      cstyp1 = 'c'      
      if (nat1.gt.maxele) cstyp1 = 's'      
      nat2 = nBOspecA2(i)
      ntype2 = nBOtypA2(i)
      call label(nat2,ntype2,lab2)
      cstyp2 = 'c'
      if (nat2.gt.maxele) cstyp2 = 's'
!
      if (lout) then
        lout = .false.
        write(ioout,'(/,''  Bond-order attractive :'',/)')
        write(ioout,'(''--------------------------------------------------------------------------------'')')
        write(ioout,'(''Atom    Type         Omega            Lambda          M                           '')')
        if (nBOtypeA(i).eq.2) then
          write(ioout,'('' 1/2     1/2           C              D                H                        '')')
        elseif (nBOtypeA(i).eq.3) then
          write(ioout,'('' 1/2     1/2           C1             C2               H                        '')')
          write(ioout,'(''                       C3             C4               C5                       '')')
        elseif (nBOtypeA(i).eq.4) then
          write(ioout,'('' 1/2     1/2           C0             C1               C2                       '')')
          write(ioout,'(''                       C3             C4                                        '')')
        else
          write(ioout,'('' 1/2     1/2                                                                    '')')
        endif
        write(ioout,'(''--------------------------------------------------------------------------------'')')
      endif
      if (nBOtypeA(i).eq.1) then
        write(ioout,'(a5,a1,6x,1x,2(f16.8,1x),5x,i3)') lab1,cstyp1,BOocoeffA(i),BOlcoeffA(i), &
          nint(BOmcoeffA(i))
        write(ioout,'(a5,a1)') lab2,cstyp2
      else
        if (nBOtypeA(i).eq.3) then
          write(ioout,'(a5,a1,2x,''Kuma-'',1x,2(f16.8,1x),4x,i3)') lab1,cstyp1,BOocoeffA(i), &
            BOlcoeffA(i),nint(BOmcoeffA(i))
          write(ioout,'(a5,a1,2x,''-gai '',1x,3(f16.8,1x))') lab2,cstyp2,BOccoeffA(1,i),BOccoeffA(2,i), &
            BOhcoeffA(i)
          write(ioout,'(14x,3(f16.8,1x))') BOccoeffA(3,i),BOccoeffA(4,i),BOccoeffA(5,i)
        elseif (nBOtypeA(i).eq.4) then
          write(ioout,'(a5,a1,2x,''MMP  '',1x,2(f16.8,1x),4x,i3)') lab1,cstyp1,BOocoeffA(i), &
            BOlcoeffA(i),nint(BOmcoeffA(i))
          write(ioout,'(a5,a1,8x,3(f16.8,1x))') lab2,cstyp2,BOccoeffA(1,i),BOccoeffA(2,i), &
            BOccoeffA(3,i)
          write(ioout,'(14x,2(f16.8,1x))') BOccoeffA(4,i),BOccoeffA(5,i)
        else
          write(ioout,'(a5,a1,2x,''Theta'',1x,2(f16.8,1x),4x,i3)') lab1,cstyp1,BOocoeffA(i), &
            BOlcoeffA(i),nint(BOmcoeffA(i))
          write(ioout,'(a5,a1,8x,3(f16.8,1x))') lab2,cstyp2,BOccoeffA(1,i),BOccoeffA(2,i),BOhcoeffA(i)
        endif
      endif
      if (.not.lout) then
        write(ioout,'(''--------------------------------------------------------------------------------'')')
      endif
    enddo
  endif
!*********************************
!  Bond-order charge potentials  *
!*********************************
  if (nboQ.gt.0) then
    lout = .true.
    do i = 1,nboQ
      nat1 = nBOspecQ1(i)
      nat2 = nBOspecQ2(i)
      ntype1 = nBOtypQ1(i)
      ntype2 = nBOtypQ2(i)
      call label(nat1,ntype1,lab1)
      call label(nat2,ntype2,lab2)
      cstyp1 = 'c'
      cstyp2 = 'c'
      if (nat1.gt.maxele) cstyp1 = 's'
      if (nat2.gt.maxele) cstyp2 = 's'
!
      if (lout) then
        lout = .false.
        write(ioout,'(/,''  Bond-order charge potentials :'',/)')
        write(ioout,'(''--------------------------------------------------------------------------------'')')
        write(ioout,'(''Atom  Types                         q0           Taper           Cutoffs(Ang)   '')')
        write(ioout,'(''  1     2                                        type            Taper    Max   '')')
        write(ioout,'(''--------------------------------------------------------------------------------'')')
      endif
      if (nBOtaperQ(i).eq.2) then
        write(ioout,'(2(a5,a1,1x),15x,f15.6,5x,''sine  '',6x,2(1x,f8.4))') lab1,cstyp1,lab2,cstyp2,BOq0(i), &
          rBOminQ(i),rBOmaxQ(i)
      else
        write(ioout,'(2(a5,a1,1x),15x,f15.6,5x,''cosine'',6x,2(1x,f8.4))') lab1,cstyp1,lab2,cstyp2,BOq0(i), &
          rBOminQ(i),rBOmaxQ(i)
      endif
      if (.not.lout) then
        write(ioout,'(''--------------------------------------------------------------------------------'')')
      endif
    enddo
  endif
!**********************************
!  Bond-order charge self-energy  *
!**********************************
  if (nboQ0.gt.0) then
    lout = .true.
    do i = 1,nboQ0
      nat1 = nBOspecQ0(i)
      ntype1 = nBOtypQ0(i)
      call label(nat1,ntype1,lab1)
      cstyp1 = 'c'
      if (nat1.gt.maxele) cstyp1 = 's'
!
      if (lout) then
        lout = .false.
        write(ioout,'(/,''  Bond-order charge self energy :'',/)')
        write(ioout,'(''--------------------------------------------------------------------------------'')')
        write(ioout,'(''Atom  Type           e0               rho              q0                       '')')
        write(ioout,'(''--------------------------------------------------------------------------------'')')
      endif
      write(ioout,'(a5,1x,a1,3(6x,f15.6))') lab1,cstyp1,BOq0pot(i),BOq0rho(i),BOq0ref(i)
      if (.not.lout) then
        write(ioout,'(''--------------------------------------------------------------------------------'')')
      endif
    enddo
  endif
!***********************************
!  Bond-order coordination number  *
!***********************************
  if (nboZ.gt.0) then
    lout = .true.
    do i = 1,nboZ
      nat1 = nBOspecZ(i)      
      ntype1 = nBOtypZ(i)      
      call label(nat1,ntype1,lab1)      
      cstyp1 = 'c'      
      if (nat1.gt.maxele) cstyp1 = 's'      
!
      if (lout) then
        lout = .false.
        write(ioout,'(/,''  Bond-order coordination :'',/)')
        write(ioout,'(''--------------------------------------------------------------------------------'')')
        write(ioout,'(''Atom              C1              C2              Z0              E0            '')')
        write(ioout,'(''--------------------------------------------------------------------------------'')')
      endif
      write(ioout,'(a5,a1,4x,4(f16.8,1x),4x,i3)') lab1,cstyp1,BOccoeffZ(1,i),BOccoeffZ(2,i),BOzcoeffZ(i), &
          BOecoeffZ(i)
      if (.not.lout) then
        write(ioout,'(''--------------------------------------------------------------------------------'')')
      endif
    enddo
  endif
!*********************
!  Plane potentials  *
!*********************
  if (nplanepot.gt.0) then
    lout = .true.
    do i = 1,nplanepot
      nat1 = natplanepot(i)
      ntype1 = ntypplanepot(i)
      call label(nat1,ntype1,lab1)
      cstyp1 = 'c'
      if (nat1.gt.maxele) cstyp1 = 's'
!
      if (lout) then
        lout = .false.
        write(ioout,'(/,''  Plane potentials (L-J):'',/)')
        write(ioout,'(''--------------------------------------------------------------------------------'')')
        write(ioout,'(''Atom  Type    m  n           z              A              B         rmin  rmax '')')
        write(ioout,'(''                           (Ang)       (eV/Ang**m)     (eV/Ang**n)   (Ang) (Ang)'')')
        write(ioout,'(''--------------------------------------------------------------------------------'')')
      endif
      write(ioout,'(a5,1x,a1,5x,2i3,3(1x,f15.6),2(1x,f6.3))') &
        lab1,cstyp1,nplanepotpower(1,i),nplanepotpower(2,i),(planepot(j,i),j=1,3),planepotrmin(i),planepotrmax(i)
      if (.not.lout) then
        write(ioout,'(''--------------------------------------------------------------------------------'')')
      endif
    enddo
  endif
!
  return
  end
