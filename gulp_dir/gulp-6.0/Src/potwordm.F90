  subroutine potwordm(iin,word,lwordok,iline,line,l55,l1000,llibrary,lfflags)
!
!  Processes potential input for many-body potential related
!  information.
!
!  iin = input fortran channel
!
!  ndenfn  = pointer to density functional form
!            1 => power law (Sutton-Chen)
!            2 => exponential
!            3 => Gaussian
!            4 => cubic 
!            5 => Voter
!            6 => Cubic
!            7 => Quartic
!            8 => Glue
!            9 => eVoter
!           10 => Mei-Davenport
!           11 => numerical
!           12 => Baskes
!           13 => VBO
!           14 => fpower (Fractional power law)
!           15 => Spline (cubic)
!           16 => Morse squared
!  denpar  = parameters (up to sixteen) for density functional form
!  neamfn  = pointer to EAM/MEAM density functional form
!            1 => square root
!            2 => general inverse power of n
!            3 => Banerjea and Smith
!            4 => numerical
!            5 => Johnson
!            6 => Glue
!            7 => Foiles
!            8 => Mei-Davenport
!            9 => Baskes
!           10 => VBO
!           11 => Belashchenko
!           12 => Igarashi
!           13 => Spline
!  neampower = power of n for functional
!  neamspec  = number of functional species
!  neamnat   = atom symbol
!  neamtyp   = atom type 
!  eamfnpar  = parameters (up to three) for functional form
!  eamfnfile = name of file for numerical functional
!  eamfnmeamcoeff = MEAM coefficients for functional
!  nmeamcombotype = indicates type of combining rule used for MEAM density components
!                   1 => sqrt of sum of densities squared
!                   2 => gamma exponential formalism
!  nmeamrhotype = indicates whether the 21 or 24 component MEAM density is to be used:
!                   1 => 21 component
!                   2 => 24 component (default)
!
!   7/97 Created from potword2
!   7/97 Optional EAM density functionals added
!   4/98 Assignment of species type for library call changed
!        so that potential file form is accepted
!  10/98 Codes for fitting variables simplified
!  10/99 Cubic density form added
!   6/01 Initialisation of line added for benefit of some compilers
!  10/02 ReaxFF mods added
!  11/03 eam_alloy option added
!  11/03 input changed to make nden & neamspec species order the same
!   7/05 Option of scaling parameter for square root functional form added
!   9/05 Voter form of EAM density added
!  10/05 Option to read numerical embedding function added
!   2/06 Quadratic and quartic densities added
!   3/06 Modified to allow for density component number
!   4/06 Species specific density added
!   4/06 Flag added to indicate whether any target species 
!        for EAM density has been given
!   5/06 Modified to accommodate separation of neamspec and neamfnspec
!   8/06 ltype01 added to call to okspec
!   8/06 iin passed to linepro
!   8/06 Separate flag added for fitting flags
!   2/07 Johnson & Glue EAM functional added
!   2/07 Glue density added
!   5/07 eVoter added
!   5/07 Foiles EAM functional added
!   7/07 reaxFF option word removed from this routine
!  11/07 Mei-Davenport density added
!  11/07 Mei-Davenport functional added
!   2/08 Option to read density from a file added
!   4/08 Error in the order of flag assignment for Mei-Davenport density fixed
!  10/08 MEAM modifications added - denpar increased in dimension
!  11/08 nfvar3 added to handle pointer to MEAM order
!  11/08 Baskes form of exponential density added
!  11/08 Baskes form of functional added
!  12/08 Module input renamed to gulpinput
!  12/08 eammeamcoeff array changed from density to function related
!  12/08 rho0 added to Baskes form of functional
!   1/09 VBO density added
!   1/09 Energy unit conversion removed from densities as inappropriate
!   1/09 VBO functional added
!   3/09 lTargetEAM removed
!   4/09 neamfnmeamtype added
!   4/09 neamfnmeamcombotype added
!   4/09 meam_screening option added
!   5/09 Analytic second derivatives disabled for screened MEAM
!   5/09 Analytic third derivatives disabled for MEAM in all forms
!   7/09 Handling of species labels modified for library case so that
!        library symbol is translated into actual symbol.
!   7/09 Saving of symbols added
!   9/09 Belashchenko EAM functional added
!   9/09 Initialization of ndenfncomp added when neamspec is increased
!  10/11 Fractional power density type added
!   6/13 Igarashi functional added
!   7/13 Cubic spline EAM density added
!   8/13 Cubic spline EAM functional added
!   3/14 Unused arguments removed
!   8/14 MEAM screening made species specific
!   8/14 neamfnmeamtype replaced by global variable nmeamrhotype 
!   8/14 neamfnmeamcombotype replaced by global variable nmeamcombotype 
!   8/14 Trio specific MEAM screening parameters allowed for
!   8/14 A and Ec separated in Baskes MEAM functional
!   8/15 Correction to lnoanald3 logic made for two options
!   9/16 eamXcutfactor added
!   1/19 Worsy changed to use word length parameter
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
  use g_constants
  use control
  use eam
  use element,     only : maxele
  use fitting
  use gulpinput
  use gulp_lengths
  use parallel
  use species
  use two,         only : eamXcutfactor
  implicit none
!
!  Passed variables
!
  character(len=maxwordlength) :: word
  character(len=maxlinelength) :: line
  integer(i4)                  :: iin
  integer(i4)                  :: iline
  logical                      :: l55
  logical                      :: l1000
  logical                      :: lfflags
  logical                      :: llibrary
  logical                      :: lwordok
!
!  Local variables
!
  character(len=5)             :: sym1
  character(len=5)             :: sym2
  character(len=5)             :: sym3
  integer(i4)                  :: i
  integer(i4)                  :: ilp1
  integer(i4)                  :: ilp2
  integer(i4)                  :: ilp3
  integer(i4)                  :: ind
  integer(i4)                  :: itype1
  integer(i4)                  :: itype2
  integer(i4)                  :: itype3
  integer(i4)                  :: iword
  integer(i4)                  :: j
  integer(i4)                  :: meamorder
  integer(i4)                  :: meamtype
  integer(i4)                  :: n1
  integer(i4)                  :: n2
  integer(i4)                  :: n3
  integer(i4)                  :: n4
  integer(i4)                  :: n5
  integer(i4)                  :: n6
  integer(i4)                  :: n7
  integer(i4)                  :: n8
  integer(i4)                  :: n9
  integer(i4)                  :: nbeg
  integer(i4)                  :: ndenfnloc
  integer(i4)                  :: neamfnspecloc
  integer(i4)                  :: neamfnspecloc2
  integer(i4)                  :: neamfnspecloc3
  integer(i4)                  :: neamspecloc
  integer(i4)                  :: nrho
  integer(i4)                  :: nvar1
  integer(i4)                  :: nvar2
  integer(i4)                  :: nvar3
  integer(i4)                  :: order
  logical                      :: lfound
  logical                      :: lok1
  logical                      :: lok2
  logical                      :: lok3
  logical                      :: lself
  logical                      :: lsymbol
  logical                      :: ltype01
  logical                      :: ltype02
  logical                      :: ltype03
  real(dp)                     :: denpar6loc
  real(dp)                     :: rdrho
  real(dp)                     :: units
!
!  Initialise local variables
!
  if (index(word,'eam_d').eq.1) goto 100
  if (index(word,'eam_f').eq.1) goto 110
  if (index(word,'eam_a').eq.1) goto 120
  if (index(word,'meam_d').eq.1) goto 130
  if (index(word,'meam_f').eq.1) goto 140
  if (index(word,'meam_s').eq.1) goto 150
  if (index(word,'meam_r').eq.1) goto 160
  if (index(word,'cutm').eq.1) goto 170
  return
!*********************************
!  Parameters for EAM densities  *
!*********************************
100 units = 1.0_dp
  ndenfnloc = 1
  if (nword.gt.1) then
    i = 2
    do while (i.le.nword)
      if (index(words(i),'exp').eq.1) then
        ndenfnloc = 2
      elseif (index(words(i),'gau').eq.1) then
        ndenfnloc = 3
      elseif (index(words(i),'pow').eq.1) then
        ndenfnloc = 1
      elseif (index(words(i),'cub').eq.1) then
        ndenfnloc = 4
      elseif (index(words(i),'vot').eq.1) then
        ndenfnloc = 5
      elseif (index(words(i),'quad').eq.1) then
        ndenfnloc = 6
      elseif (index(words(i),'quar').eq.1) then
        ndenfnloc = 7
      elseif (index(words(i),'glue').eq.1) then
        ndenfnloc = 8
      elseif (index(words(i),'evo').eq.1) then
        ndenfnloc = 9
      elseif (index(words(i),'mei').eq.1) then
        ndenfnloc = 10
      elseif (index(words(i),'num').eq.1) then
        ndenfnloc = 11
      elseif (index(words(i),'bas').eq.1) then
        ndenfnloc = 12
      elseif (index(words(i),'vbo').eq.1) then
        ndenfnloc = 13
      elseif (index(words(i),'fpo').eq.1) then
        ndenfnloc = 14
      elseif (index(words(i),'spl').eq.1) then
        ndenfnloc = 15
      elseif (index(words(i),'mor').eq.1) then
        ndenfnloc = 16
      endif
      i = i + 1
    enddo
  endif
  if (nfloat.gt.0) then
    denpar6loc = nint(floats(1))
  elseif (ndenfnloc.eq.1) then
    denpar6loc = 6
  else
    denpar6loc = 0
  endif
105 line = '  '
  read(iin,'(a)',end=108) line
  iline = iline + 1
  call linepro(iin,line,iline)
!
!  Check whether input is symbol or option
!
  if ((nword+nfloat).eq.0) goto 105
  if (nword.gt.0) then
    word = words(1)
    call worsy(word,lsymbol,.true.)
    if (.not.lsymbol) then
      l55 = .true.
      goto 108
    endif
  endif
!
!  Symbols used in input
!
  nvar2 = 0
  itype2 = 0
  ltype01 = .false.
  sym1 = ' '
  sym2 = ' '
  if (nword.gt.0) then
    if (llibrary) then
      call okspec(lok1,words(1),ilp1,.true.,ltype01)
      if (.not.lok1) then
        goto 105
      endif
      if (ilp1.gt.0.and.ilp1.le.nspec) then
        nvar1 = natspec(ilp1)
        if (ltype01) then
          itype1 = 0
        else
          itype1 = ntypspec(ilp1)
        endif
      elseif (ilp1.eq.-1) then
        nvar1 = maxele
        itype1 = 0
      endif
      if (ltype01) itype1 = 0
    else
      call ltont(words(1),nvar1,itype1)
    endif
    sym1 = words(1)(1:5)
    if (nword.ge.4) then
      if (index(words(2),'she').eq.1) then
        nvar1 = nvar1 + maxele
      endif
      call ltont(words(3),nvar2,itype2)
      if (index(words(4),'she').eq.1) then
        nvar2 = nvar2 + maxele
      endif
      sym2 = words(3)(1:5)
    elseif (nword.eq.3) then
      if (index(words(2),'she').eq.1) then
        nvar1 = nvar1 + maxele
      elseif (index(words(2),'cor').eq.0) then
        call ltont(words(2),nvar2,itype2)
        if (index(words(3),'she').eq.1) then
          nvar2 = nvar2 + maxele
        endif
        sym2 = words(2)(1:5)
      else
        call ltont(words(3),nvar2,itype2)
        sym2 = words(3)(1:5)
      endif
    elseif (nword.eq.2) then
      if (index(words(2),'she').eq.1) then
        nvar1 = nvar1 + maxele
      elseif (index(words(2),'cor').eq.0) then
        call ltont(words(2),nvar2,itype2)
        sym2 = words(2)(1:5)
      endif
    endif
    nbeg = 0
  else
!
!  Numeric input
!
    nvar1 = int(floats(1))
    if (nvar1.gt.100) nvar1 = nvar1 - 100 + maxele
    itype1 = 0
    nbeg = 1
    nfloat = nfloat - 1
    call label(nvar1,itype1,sym1)
  endif
!
!  Check to see whether this atomic species is already in the list and if not add it
!
  i = 0
  lfound = .false.       
  do while (i.lt.neamspec.and..not.lfound)
    i = i + 1
    lfound = (neamnat(i).eq.nvar1.and.neamtyp(i).eq.itype1.and.neamnat2(i).eq.nvar2.and.neamtyp2(i).eq.itype2)
  enddo    
  if (lfound) then
    neamspecloc = i
  else
    neamspec = neamspec + 1
    neamspecloc = neamspec
    if (neamspec.gt.maxeamspec) then
      maxeamspec = neamspec + 20
      call changemaxeamspec
    endif     
!
!  When adding a new species initialize related variables
!
    ndenfncomp(neamspec) = 0
!
    neamnat(neamspec) = nvar1
    neamtyp(neamspec) = itype1
    neamnat2(neamspec) = nvar2
    neamtyp2(neamspec) = itype2
    symboleamspec(1,neamspec) = sym1
    symboleamspec(2,neamspec) = sym2
  endif
!
!  Copy first line info from first density
!
  ndenfncomp(neamspecloc) = ndenfncomp(neamspecloc) + 1
  if (ndenfncomp(neamspecloc).gt.maxeamden) then
    maxeamden = ndenfncomp(neamspecloc) + 3
    call changemaxeamden
  endif
  ndenfn(ndenfncomp(neamspecloc),neamspecloc) = ndenfnloc
  denpar(6,1,ndenfncomp(neamspecloc),neamspecloc) = denpar6loc
!
!  Fitting flags
!
  if (lfit.and.lfflags) then
    if (ndenfnloc.eq.2.or.ndenfnloc.eq.3.or.ndenfnloc.eq.12.or.ndenfnloc.eq.16) then
      n1 = int(floats(nfloat-2))
      n2 = int(floats(nfloat-1))
      n3 = int(floats(nfloat))
      nfloat = nfloat - 3
      if (n1.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamspecloc
        nfvar(nfit) = 1
        nfvar2(nfit) = ndenfncomp(neamspecloc)
        nfvar3(nfit) = 1
      endif
      if (n2.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamspecloc
        nfvar(nfit) = 2
        nfvar2(nfit) = ndenfncomp(neamspecloc)
        nfvar3(nfit) = 1
      endif
      if (n3.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamspecloc
        nfvar(nfit) = 3
        nfvar2(nfit) = ndenfncomp(neamspecloc)
        nfvar3(nfit) = 1
      endif
    elseif (ndenfnloc.eq.1) then
      n1 = int(floats(nfloat))
      nfloat = nfloat - 1
      if (n1.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamspecloc
        nfvar(nfit) = 1
        nfvar2(nfit) = ndenfncomp(neamspecloc)
        nfvar3(nfit) = 1
      endif
    elseif (ndenfnloc.eq.4.or.ndenfnloc.eq.5.or.ndenfnloc.eq.6.or.ndenfnloc.eq.7.or.ndenfnloc.eq.9.or.ndenfnloc.eq.14) then
      n1 = int(floats(nfloat-1))
      n2 = int(floats(nfloat))
      nfloat = nfloat - 2
      if (n1.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamspecloc
        nfvar(nfit) = 1
        nfvar2(nfit) = ndenfncomp(neamspecloc)
        nfvar3(nfit) = 1
      endif
      if (n2.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamspecloc
        nfvar(nfit) = 2
        nfvar2(nfit) = ndenfncomp(neamspecloc)
        nfvar3(nfit) = 1
      endif
    elseif (ndenfnloc.eq.10) then
      do n1 = 1,7
        n2 = int(floats(nfloat+n1-7))
        if (n2.eq.1) then 
          nfit = nfit + 1
          if (nfit.ge.maxfit) then
            maxfit = nfit + 10
            call changemaxfit
          endif
          nftyp(nfit) = 5
          nfpot(nfit) = neamspecloc
          if (n1.le.3) then
            nfvar(nfit) = n1
          else
            nfvar(nfit) = n1 + 11
          endif
          nfvar2(nfit) = ndenfncomp(neamspecloc)
          nfvar3(nfit) = 1
        endif
      enddo
      nfloat = nfloat - 7
    elseif (ndenfnloc.eq.13) then
      do n1 = 1,5
        n2 = int(floats(nfloat+n1-5))
        if (n2.eq.1) then 
          nfit = nfit + 1
          if (nfit.ge.maxfit) then
            maxfit = nfit + 10
            call changemaxfit
          endif
          nftyp(nfit) = 5
          nfpot(nfit) = neamspecloc
          if (n1.le.3) then
            nfvar(nfit) = n1
          else
            nfvar(nfit) = n1 + 11
          endif
          nfvar2(nfit) = ndenfncomp(neamspecloc)
          nfvar3(nfit) = 1
        endif
      enddo
      nfloat = nfloat - 5
    elseif (ndenfnloc.eq.15) then
      do n1 = 1,4
        n2 = int(floats(nfloat+n1-4))
        if (n2.eq.1) then
          nfit = nfit + 1
          if (nfit.ge.maxfit) then
            maxfit = nfit + 10
            call changemaxfit
          endif
          nftyp(nfit) = 5
          nfpot(nfit) = neamspecloc
          if (n1.le.3) then
            nfvar(nfit) = n1
          else
            nfvar(nfit) = n1 + 11
          endif
          nfvar2(nfit) = ndenfncomp(neamspecloc)
          nfvar3(nfit) = 1
        endif
      enddo
      nfloat = nfloat - 4
    endif
  endif
  if (ndenfnloc.eq.1) then
!
!  If number of floats is greater than 1, assume that fitting flags have been left on line
!
    if (nfloat.gt.1) nfloat = nfloat - 1
!
!  Assign coefficients 
!
    if (nfloat.ge.1) then
      denpar(1,1,ndenfncomp(neamspecloc),neamspecloc) = floats(1+nbeg)*units
    else
      call outerror('Incorrect coefficient input for EAM density',iline)
      call stopnow('potwordm')
    endif
  elseif (ndenfnloc.eq.4.or.ndenfnloc.eq.5.or.ndenfnloc.eq.6.or.ndenfnloc.eq.7.or.ndenfnloc.eq.9.or.ndenfnloc.eq.14) then
!
!  If number of floats is greater than 2, assume that fitting flags have been left on line
!
    if (nfloat.gt.2) nfloat = nfloat - 2
!
!  Assign coefficients 
!
    if (nfloat.ge.2) then
      denpar(1,1,ndenfncomp(neamspecloc),neamspecloc) = floats(1+nbeg)*units
      denpar(2,1,ndenfncomp(neamspecloc),neamspecloc) = floats(2+nbeg)
    elseif (nfloat.eq.1) then
      denpar(1,1,ndenfncomp(neamspecloc),neamspecloc) = floats(1+nbeg)*units
      denpar(2,1,ndenfncomp(neamspecloc),neamspecloc) = 0.0_dp
    else
      call outerror('Incorrect coefficient input for EAM density',iline)
      call stopnow('potwordm')
    endif
  elseif (ndenfnloc.eq.10) then
!
!  If number of floats is greater than 7, assume that fitting flags have been left on line
!
    if (nfloat.gt.7) nfloat = nfloat - 7
!
!  Assign coefficients 
!
    if (nfloat.ge.7) then
      denpar(1,1,ndenfncomp(neamspecloc),neamspecloc) = floats(1+nbeg)*units
      denpar(2,1,ndenfncomp(neamspecloc),neamspecloc) = floats(2+nbeg)*units
      denpar(3,1,ndenfncomp(neamspecloc),neamspecloc) = floats(3+nbeg)*units
      denpar(4,1,ndenfncomp(neamspecloc),neamspecloc) = floats(4+nbeg)*units
      denpar(5,1,ndenfncomp(neamspecloc),neamspecloc) = floats(5+nbeg)*units
      denpar(6,1,ndenfncomp(neamspecloc),neamspecloc) = floats(6+nbeg)*units
      denpar(7,1,ndenfncomp(neamspecloc),neamspecloc) = floats(7+nbeg)
    else
      call outerror('Incorrect coefficient input for EAM density',iline)
      call stopnow('potwordm')
    endif
    if (abs(denpar(7,1,ndenfncomp(neamspecloc),neamspecloc)).lt.1.0d-12) then
      call outerror('r0 is too close to zero for EAM density',iline)
      call stopnow('potwordm')
    endif
  elseif (ndenfnloc.eq.13) then
!
!  VBO density
!
!  If number of floats is greater than 5, assume that fitting flags have been left on line
!
    if (nfloat.gt.5) nfloat = nfloat - 5
!
!  Assign coefficients 
!
    if (nfloat.ge.5) then
      denpar(1,1,ndenfncomp(neamspecloc),neamspecloc) = floats(1+nbeg)*units
      denpar(2,1,ndenfncomp(neamspecloc),neamspecloc) = floats(2+nbeg)
      denpar(3,1,ndenfncomp(neamspecloc),neamspecloc) = floats(3+nbeg)
      denpar(4,1,ndenfncomp(neamspecloc),neamspecloc) = abs(floats(4+nbeg))
      denpar(5,1,ndenfncomp(neamspecloc),neamspecloc) = abs(floats(5+nbeg))
    else
      call outerror('Incorrect coefficient input for EAM density',iline)
      call stopnow('potwordm')
    endif
    if (abs(denpar(4,1,ndenfncomp(neamspecloc),neamspecloc)).lt.1.0d-12) then
      call outerror('r0 is too close to zero for EAM density',iline)
      call stopnow('potwordm')
    endif
    if (abs(denpar(5,1,ndenfncomp(neamspecloc),neamspecloc)).lt.1.0d-12) then
      call outerror('delta is too close to zero for EAM density',iline)
      call stopnow('potwordm')
    endif
  elseif (ndenfnloc.eq.15) then
!
!  Spline density
!
!  If number of floats is greater than 6, assume that fitting flags have been left on line
!
    if (nfloat.gt.6) nfloat = nfloat - 4
!
!  Assign coefficients
!
    if (nfloat.ge.6) then
      denpar(1,1,ndenfncomp(neamspecloc),neamspecloc) = floats(1+nbeg)*units
      denpar(2,1,ndenfncomp(neamspecloc),neamspecloc) = floats(2+nbeg)*units
      denpar(3,1,ndenfncomp(neamspecloc),neamspecloc) = floats(3+nbeg)*units
      denpar(4,1,ndenfncomp(neamspecloc),neamspecloc) = floats(4+nbeg)*units
      denpar(5,1,ndenfncomp(neamspecloc),neamspecloc) = abs(floats(5+nbeg))
      denpar(6,1,ndenfncomp(neamspecloc),neamspecloc) = abs(floats(6+nbeg))
    else
      call outerror('Incorrect coefficient input for EAM density',iline)
      call stopnow('potwordm')
    endif
    if (abs(denpar(5,1,ndenfncomp(neamspecloc),neamspecloc)).lt.1.0d-12) then
      call outerror('rmin is too close to zero for EAM density',iline)
      call stopnow('potwordm')
    endif
    if (abs(denpar(6,1,ndenfncomp(neamspecloc),neamspecloc)).lt.1.0d-12) then
      call outerror('r0 is too close to zero for EAM density',iline)
      call stopnow('potwordm')
    endif
  elseif (ndenfnloc.eq.8) then
!
!  Assign distance cutoffs : d / r_b / r_m
!
    if (nfloat.ge.3) then
      denpar(1,1,ndenfncomp(neamspecloc),neamspecloc) = floats(1+nbeg)
      denpar(2,1,ndenfncomp(neamspecloc),neamspecloc) = floats(2+nbeg)
      denpar(3,1,ndenfncomp(neamspecloc),neamspecloc) = floats(3+nbeg)
    else
      call outerror('Incorrect coefficient input for EAM density',iline)
      call stopnow('potwordm')
    endif
!
!  Read second line
!
    line = '  '
    read(iin,'(a)',end=108) line
    iline = iline + 1
    call linepro(iin,line,iline)
!
!  Assign coefficients 
!
    if (nfloat.ge.4) then
      denpar(4,1,ndenfncomp(neamspecloc),neamspecloc) = floats(1+nbeg)*units
      denpar(5,1,ndenfncomp(neamspecloc),neamspecloc) = floats(2+nbeg)*units
      denpar(6,1,ndenfncomp(neamspecloc),neamspecloc) = floats(3+nbeg)*units
      denpar(7,1,ndenfncomp(neamspecloc),neamspecloc) = floats(4+nbeg)*units
    else
      call outerror('Incorrect coefficient input for EAM density',iline)
      call stopnow('potwordm')
    endif
!
!  Read third line
!
    line = '  '
    read(iin,'(a)',end=108) line
    iline = iline + 1
    call linepro(iin,line,iline)
!
!  Assign coefficients 
!
    if (nfloat.ge.4) then
      denpar(8,1,ndenfncomp(neamspecloc),neamspecloc)  = floats(1+nbeg)*units
      denpar(9,1,ndenfncomp(neamspecloc),neamspecloc)  = floats(2+nbeg)*units
      denpar(10,1,ndenfncomp(neamspecloc),neamspecloc) = floats(3+nbeg)*units
      denpar(11,1,ndenfncomp(neamspecloc),neamspecloc) = floats(4+nbeg)*units
    else
      call outerror('Incorrect coefficient input for EAM density',iline)
      call stopnow('potwordm')
    endif
!
!  Read fourth line
!
    line = '  '
    read(iin,'(a)',end=108) line
    iline = iline + 1
    call linepro(iin,line,iline)
!
!  Assign coefficients 
!
    if (nfloat.ge.4) then
      denpar(12,1,ndenfncomp(neamspecloc),neamspecloc) = floats(1+nbeg)*units
      denpar(13,1,ndenfncomp(neamspecloc),neamspecloc) = floats(2+nbeg)*units
      denpar(14,1,ndenfncomp(neamspecloc),neamspecloc) = floats(3+nbeg)*units
      denpar(15,1,ndenfncomp(neamspecloc),neamspecloc) = floats(4+nbeg)*units
    else
      call outerror('Incorrect coefficient input for EAM density',iline)
      call stopnow('potwordm')
    endif
  elseif (ndenfnloc.eq.11) then
!********************************
!  Numerical density from file  *
!********************************
!
!  Assign file name 
!
!    eamdenfile(neamspecloc) = words(iword)
!
!  Check that file exists by opening and reading
!
!    open(9,file=eamfnfile(neamfnspecloc),status='old',form='formatted',err=119)
!    read(9,'(a)') line
!    read(9,'(a)') line
!    read(9,'(i5,e24.16)') neamfnnumeric(neamfnspecloc),eamfnnumericdrho(neamfnspecloc)
  else
!
!  If number of floats is greater than 3, assume that fitting flags have been left on line
!
    if (nfloat.gt.3) nfloat = nfloat - 3
!
!  Assign coefficients 
!
    if (nfloat.ge.3) then
      denpar(1,1,ndenfncomp(neamspecloc),neamspecloc) = floats(1+nbeg)*units
      denpar(2,1,ndenfncomp(neamspecloc),neamspecloc) = floats(2+nbeg)
      denpar(3,1,ndenfncomp(neamspecloc),neamspecloc) = floats(3+nbeg)
    elseif (nfloat.eq.2) then
      denpar(1,1,ndenfncomp(neamspecloc),neamspecloc) = floats(1+nbeg)*units
      denpar(2,1,ndenfncomp(neamspecloc),neamspecloc) = floats(2+nbeg)
      denpar(3,1,ndenfncomp(neamspecloc),neamspecloc) = 0.0_dp
    else
      call outerror('Incorrect coefficient input for EAM density',iline)
      call stopnow('potwordm')
    endif
  endif
  goto 105
108 if (.not.l55) l1000 = .true.
  lwordok = .true.
  return
!***********************************
!  Parameters for EAM functionals  *
!***********************************
110 units = 1.0_dp
  if (nword.gt.1) then
    i = 2
    do while (i.le.nword)
      if (index(words(i),'sq').eq.1) then
        neamfn = 1
        neampower = 2
      elseif (index(words(i),'po').eq.1) then
        neamfn = 2
      elseif (index(words(i),'ban').eq.1) then
        neamfn = 3
      elseif (index(words(i),'nu').eq.1) then
        neamfn = 4
      elseif (index(words(i),'jo').eq.1) then
        neamfn = 5
      elseif (index(words(i),'gl').eq.1) then
        neamfn = 6
      elseif (index(words(i),'fo').eq.1) then
        neamfn = 7
      elseif (index(words(i),'me').eq.1) then
        neamfn = 8
      elseif (index(words(i),'bas').eq.1) then
        neamfn = 9
      elseif (index(words(i),'vb').eq.1) then
        neamfn = 10
      elseif (index(words(i),'be').eq.1) then
        neamfn = 11
      elseif (index(words(i),'ig').eq.1) then
        neamfn = 12
      elseif (index(words(i),'sp').eq.1) then
        neamfn = 13
      elseif (index(words(i),'au').eq.1) then
        units = autoev
      elseif (index(words(i),'kc').eq.1) then
        units = kcaltoev
      elseif (index(words(i),'kj').eq.1) then
        units = kjmtoev
      endif
      i = i + 1
    enddo
  endif
  if (nfloat.gt.0.and.neamfn.gt.1) then
    neampower = nint(floats(1))
  endif
  if (neamfn.gt.1.and.neampower.eq.0) then
    call outerror('inverse power of zero in EAM functional',iline)
    call stopnow('potwordm')
  endif
!
!  Read in species parameters for functional
!
115 line = '  '
  read(iin,'(a)',end=118) line
  iline = iline + 1
  call linepro(iin,line,iline)
!
!  Check whether input is symbol or option
!
  if ((nword+nfloat).eq.0) goto 115
  if (nword.gt.0) then
    word = words(1)
    call worsy(word,lsymbol,.true.)
    if (.not.lsymbol) then
      l55 = .true.
      goto 118
    endif
  endif
!
!  Symbol used in input
!
  ltype01 = .false.
  sym1 = ' '
  if (nword.gt.0) then
    if (llibrary) then
      call okspec(lok1,words(1),ilp1,.true.,ltype01)
      if (.not.lok1) then
        goto 115
      endif
      if (ilp1.gt.0.and.ilp1.le.nspec) then
        nvar1 = natspec(ilp1)
        if (ltype01) then
          itype1 = 0
        else
          itype1 = ntypspec(ilp1)
        endif
      elseif (ilp1.eq.-1) then
        nvar1 = maxele
        itype1 = 0
      endif
    else
      call ltont(words(1),nvar1,itype1)
    endif
    sym1 = words(1)(1:5)
    if (nword.ge.2) then
      if (index(words(2),'she').eq.1) then
        nvar1 = nvar1 + maxele
      endif
    endif
    nbeg = 0
  else
!
!  Numeric input
!
    nvar1 = int(floats(1))
    if (nvar1.gt.100) nvar1 = nvar1 - 100 + maxele
    itype1 = 0
    nbeg = 1
    nfloat = nfloat - 1
    call label(nvar1,itype1,sym1)
  endif
!
!  Check to see whether this atomic species is already in the list and if not add it
!
  i = 0
  lfound = .false.       
  do while (i.lt.neamfnspec.and..not.lfound)
    i = i + 1
    lfound = (neamfnnat(i).eq.nvar1.and.neamfntyp(i).eq.itype1)
  enddo    
  if (lfound) then
    neamfnspecloc = i
  else
    neamfnspec = neamfnspec + 1
    neamfnspecloc = neamfnspec
    if (neamfnspec.gt.maxeamfnspec) then
      maxeamfnspec = neamfnspec + 20
      call changemaxeamfnspec
    endif     
    neamfnnat(neamfnspec) = nvar1
    neamfntyp(neamfnspec) = itype1
    symboleamfnspec(neamfnspec) = sym1
  endif
!
  lwordok = .true.
  if (neamfn.eq.3) then
!***********************
!  Banerjea and Smith  *
!***********************
!
!  Fitting flags
!
    if (lfit.and.lfflags) then
      n1 = int(floats(nfloat-2))
      n2 = int(floats(nfloat-1))
      n3 = int(floats(nfloat))
      nfloat = nfloat - 3
      if (n1.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 4
      endif
      if (n2.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 5
      endif
      if (n3.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 6
      endif
    endif
!
!  If number of floats is greater than 3, assume that fitting flags have been left on line
!
    if (nfloat.gt.3) nfloat = nfloat - 3
!
!  Assign coefficients 
!
    if (nfloat.ge.3) then
      eamfnpar(1,neamfnspecloc) = floats(1+nbeg)
      eamfnpar(2,neamfnspecloc) = floats(2+nbeg)
      eamfnpar(3,neamfnspecloc) = floats(3+nbeg)
    elseif (nfloat.eq.2) then
      eamfnpar(1,neamfnspecloc) = floats(1+nbeg)
      eamfnpar(2,neamfnspecloc) = floats(2+nbeg)
      eamfnpar(3,neamfnspecloc) = 0.0_dp
    else
      call outerror('Incorrect coefficient input for EAM functional',iline)
      call stopnow('potwordm')
    endif
    if (abs(eamfnpar(3,neamfnspecloc)).lt.1.0d-12) then
      call outerror('equilibrium density too close to zero in EAM functional',iline)
      call stopnow('potwordm')
    endif
  elseif (neamfn.eq.4) then
!*****************************
!  Numerical EAM functional  *
!*****************************
!
!  Is second word core or shel? If so skip
!
    iword = 2
    if (index(words(2),'cor').eq.1.or.index(words(2),'she').eq.1) iword = 3
    if (iword.gt.nword) then
      call outerror('filename missing for numerical EAM',iline)
      call stopnow('potwordm')
    endif
!
!  Assign file name 
!
    eamfnfile(neamfnspecloc) = words(iword)
!
!  Check that file exists by opening and reading
!
    open(9,file=eamfnfile(neamfnspecloc),status='old',form='formatted',err=119)
    read(9,'(a)') line
    read(9,'(a)') line
    read(9,'(i5,e24.16)') neamfnnumeric(neamfnspecloc),eamfnnumericdrho(neamfnspecloc)
!
!  Check the size of the array for numeric function
!
    if (neamfnnumeric(neamfnspecloc).gt.maxneamfnnumeric) then
      maxneamfnnumeric = neamfnnumeric(neamfnspecloc)
      call changemaxneamfnnumeric
    endif
!
!  Read embedding function data
!
    read(9,'(5e24.16)') (eamfnnumeric(j,neamfnspecloc),j=1,neamfnnumeric(neamfnspecloc))
!
!  Close file having read all relevant information 
!
    close(9)
!
!  Interpolation setup
!
    nrho = neamfnnumeric(neamfnspecloc)
    eamfnnumeric1(1,neamfnspecloc) = eamfnnumeric(2,neamfnspecloc) - eamfnnumeric(1,neamfnspecloc)
    eamfnnumeric1(2,neamfnspecloc) = 0.5_dp*(eamfnnumeric(3,neamfnspecloc) - eamfnnumeric(1,neamfnspecloc))
    eamfnnumeric1(nrho-1,neamfnspecloc) = &
      0.5_dp*(eamfnnumeric(nrho,neamfnspecloc) - eamfnnumeric(nrho-2,neamfnspecloc))
    eamfnnumeric1(nrho,neamfnspecloc) = eamfnnumeric(nrho,neamfnspecloc) - eamfnnumeric(nrho-1,neamfnspecloc)
    do j = 3,nrho-2
      eamfnnumeric1(j,neamfnspecloc) = ((eamfnnumeric(j-2,neamfnspecloc) - eamfnnumeric(j+2,neamfnspecloc)) + &
        8.0_dp*(eamfnnumeric(j+1,neamfnspecloc) - eamfnnumeric(j-1,neamfnspecloc)))/12.0_dp
    enddo 
    do j = 1,nrho-1
      eamfnnumeric2(j,neamfnspecloc) = 3.0_dp*(eamfnnumeric(j+1,neamfnspecloc) - eamfnnumeric(j,neamfnspecloc)) &
        - 2.0_dp*eamfnnumeric1(j,neamfnspecloc) - eamfnnumeric1(j+1,neamfnspecloc)
      eamfnnumeric3(j,neamfnspecloc) = eamfnnumeric1(j,neamfnspecloc) + eamfnnumeric1(j+1,neamfnspecloc) &
        - 2.0_dp*(eamfnnumeric(j+1,neamfnspecloc) - eamfnnumeric(j,neamfnspecloc))
    enddo
    eamfnnumeric2(nrho,neamfnspecloc) = 0.0_dp
    eamfnnumeric3(nrho,neamfnspecloc) = 0.0_dp
    rdrho = 1.0_dp/eamfnnumericdrho(neamfnspecloc)
!
!  Force interpolation set up
!
    do j = 1,nrho
      eamfnnumeric4(j,neamfnspecloc) = eamfnnumeric1(j,neamfnspecloc)*rdrho
      eamfnnumeric5(j,neamfnspecloc) = 2.0_dp*eamfnnumeric2(j,neamfnspecloc)*rdrho
      eamfnnumeric6(j,neamfnspecloc) = 3.0_dp*eamfnnumeric3(j,neamfnspecloc)*rdrho
    enddo
!
!  Second derivative interpolation set up
!
    do j = 1,nrho
      eamfnnumeric7(j,neamfnspecloc) = eamfnnumeric5(j,neamfnspecloc)*rdrho
      eamfnnumeric8(j,neamfnspecloc) = 2.0_dp*eamfnnumeric6(j,neamfnspecloc)*rdrho
    enddo
  elseif (neamfn.eq.5) then
!************
!  Johnson  *
!************
!
!  Fitting flags
!
    if (lfit.and.lfflags) then
      n1 = int(floats(nfloat-5))
      n2 = int(floats(nfloat-4))
      n3 = int(floats(nfloat-3))
      n4 = int(floats(nfloat-2))
      n5 = int(floats(nfloat-1))
      n6 = int(floats(nfloat))
      nfloat = nfloat - 6
      if (n1.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 4
      endif
      if (n2.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 5
      endif
      if (n3.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 6
      endif
      if (n4.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 9
      endif
      if (n5.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 10
      endif
      if (n6.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 11
      endif
    endif
!
!  If number of floats is greater than 6, assume that fitting flags have been left on line
!
    if (nfloat.gt.6) nfloat = nfloat - 6
!
!  Assign coefficients 
!
    if (nfloat.ge.6) then
      eamfnpar(1,neamfnspecloc) = floats(1+nbeg)
      eamfnpar(2,neamfnspecloc) = floats(2+nbeg)
      eamfnpar(3,neamfnspecloc) = floats(3+nbeg)
      eamfnpar(4,neamfnspecloc) = floats(4+nbeg)
      eamfnpar(5,neamfnspecloc) = floats(5+nbeg)
      eamfnpar(6,neamfnspecloc) = floats(6+nbeg)
    elseif (nfloat.eq.5) then
      eamfnpar(1,neamfnspecloc) = floats(1+nbeg)
      eamfnpar(2,neamfnspecloc) = 0.0_dp
      eamfnpar(3,neamfnspecloc) = floats(2+nbeg)
      eamfnpar(4,neamfnspecloc) = floats(3+nbeg)
      eamfnpar(5,neamfnspecloc) = floats(4+nbeg)
      eamfnpar(6,neamfnspecloc) = floats(5+nbeg)
    else
      call outerror('Incorrect coefficient input for EAM functional',iline)
      call stopnow('potwordm')
    endif
    if (abs(eamfnpar(3,neamfnspecloc)).lt.1.0d-12) then
      call outerror('equilibrium density too close to zero in EAM functional',iline)
      call stopnow('potwordm')
    endif
    if (abs(eamfnpar(5,neamfnspecloc)).lt.1.0d-12) then
      call outerror('denominator for density power is zero in EAM functional',iline)
      call stopnow('potwordm')
    endif
  elseif (neamfn.eq.6) then
!*********
!  Glue  *
!*********
!
!  Assign density cutoffs
!
    if (nfloat.ge.2) then
      eamfnpar(1,neamfnspecloc) = floats(1+nbeg)
      eamfnpar(2,neamfnspecloc) = floats(2+nbeg)
    else
      call outerror('Incorrect coefficient input for EAM functional',iline)
      call stopnow('potwordm')
    endif
!
!  Read second line
!
    line = '  '
    read(iin,'(a)',end=118) line
    iline = iline + 1
    call linepro(iin,line,iline)
!
!  Assign coefficients 
!
    if (nfloat.ge.5) then
      eamfnpar(3,neamfnspecloc) = floats(1+nbeg)
      eamfnpar(4,neamfnspecloc) = floats(2+nbeg)
      eamfnpar(5,neamfnspecloc) = floats(3+nbeg)
      eamfnpar(6,neamfnspecloc) = floats(4+nbeg)
      eamfnpar(7,neamfnspecloc) = floats(5+nbeg)
    else
      call outerror('Incorrect coefficient input for EAM functional',iline)
      call stopnow('potwordm')
    endif
!
!  Read third line
!
    line = '  '
    read(iin,'(a)',end=118) line
    iline = iline + 1
    call linepro(iin,line,iline)
!
!  Assign coefficients 
!
    if (nfloat.ge.5) then
      eamfnpar(8,neamfnspecloc) = floats(1+nbeg)
      eamfnpar(9,neamfnspecloc) = floats(2+nbeg)
      eamfnpar(10,neamfnspecloc) = floats(3+nbeg)
      eamfnpar(11,neamfnspecloc) = floats(4+nbeg)
      eamfnpar(12,neamfnspecloc) = floats(5+nbeg)
    else
      call outerror('Incorrect coefficient input for EAM functional',iline)
      call stopnow('potwordm')
    endif
!
!  Read fourth line
!
    line = '  '
    read(iin,'(a)',end=118) line
    iline = iline + 1
    call linepro(iin,line,iline)
!
!  Assign coefficients 
!
    if (nfloat.ge.4) then
      eamfnpar(13,neamfnspecloc) = floats(1+nbeg)
      eamfnpar(14,neamfnspecloc) = floats(2+nbeg)
      eamfnpar(15,neamfnspecloc) = floats(3+nbeg)
      eamfnpar(16,neamfnspecloc) = floats(4+nbeg)
    else
      call outerror('Incorrect coefficient input for EAM functional',iline)
      call stopnow('potwordm')
    endif
  elseif (neamfn.eq.7) then
!***********
!  Foiles  *
!***********
!
!  Fitting flags
!
    if (lfit.and.lfflags) then
      n1 = int(floats(nfloat-3))
      n2 = int(floats(nfloat-2))
      n3 = int(floats(nfloat-1))
      n4 = int(floats(nfloat))
      nfloat = nfloat - 4
      if (n1.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 4
      endif
      if (n2.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 5
      endif
      if (n3.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 6
      endif
      if (n4.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 9
      endif
    endif
!
!  If number of floats is greater than 4, assume that fitting flags have been left on line
!
    if (nfloat.gt.4) nfloat = nfloat - 4
!
!  Assign coefficients 
!
    if (nfloat.ge.4) then
      eamfnpar(1,neamfnspecloc) = floats(1+nbeg)
      eamfnpar(2,neamfnspecloc) = floats(2+nbeg)
      eamfnpar(3,neamfnspecloc) = floats(3+nbeg)
      eamfnpar(4,neamfnspecloc) = floats(4+nbeg)
    else
      call outerror('Incorrect coefficient input for EAM functional',iline)
      call stopnow('potwordm')
    endif
  elseif (neamfn.eq.8) then
!******************
!  Mei-Davenport  *
!******************
!
!  Fitting flags
!
    if (lfit.and.lfflags) then
      n1 = int(floats(nfloat-8))
      n2 = int(floats(nfloat-7))
      n3 = int(floats(nfloat-6))
      n4 = int(floats(nfloat-5))
      n5 = int(floats(nfloat-4))
      n6 = int(floats(nfloat-3))
      n7 = int(floats(nfloat-2))
      n8 = int(floats(nfloat-1))
      n9 = int(floats(nfloat))
      nfloat = nfloat - 9
      if (n1.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 4
      endif
      if (n2.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 5
      endif
      if (n3.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 6
      endif
      if (n4.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 9
      endif
      if (n5.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 10
      endif
      if (n6.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 11
      endif
      if (n7.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 12
      endif
      if (n8.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 13
      endif
      if (n9.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 14
      endif
    endif
!
!  If number of floats is greater than 8, assume that fitting flags have been left on line
!
    if (nfloat.gt.9) nfloat = nfloat - 9
!
!  Assign coefficients 
!
    if (nfloat.ge.8) then
      eamfnpar(1,neamfnspecloc) = floats(1+nbeg)*units
      eamfnpar(2,neamfnspecloc) = floats(2+nbeg)
      eamfnpar(3,neamfnspecloc) = floats(3+nbeg)
      eamfnpar(4,neamfnspecloc) = floats(4+nbeg)
      eamfnpar(5,neamfnspecloc) = floats(5+nbeg)
      eamfnpar(6,neamfnspecloc) = floats(6+nbeg)*units
      eamfnpar(7,neamfnspecloc) = floats(7+nbeg)
      eamfnpar(8,neamfnspecloc) = floats(8+nbeg)
      eamfnpar(9,neamfnspecloc) = floats(9+nbeg)
    else
      call outerror('Incorrect coefficient input for EAM functional',iline)
      call stopnow('potwordm')
    endif
  elseif (neamfn.eq.9) then
!***********
!  Baskes  *
!***********
!
!  Fitting flags
!
    if (lfit.and.lfflags) then
      n1 = int(floats(nfloat-2))
      n2 = int(floats(nfloat-1))
      n3 = int(floats(nfloat))
      nfloat = nfloat - 3
      if (n1.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 4
      endif
      if (n2.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 5
      endif
      if (n3.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 6
      endif
    endif
!
!  If number of floats is greater than 3, assume that fitting flags have been left on line
!
    if (nfloat.gt.3) nfloat = nfloat - 3
!
!  Assign coefficients 
!
    if (nfloat.ge.3) then
      eamfnpar(1,neamfnspecloc) = floats(1+nbeg)
      eamfnpar(2,neamfnspecloc) = floats(2+nbeg)
      eamfnpar(3,neamfnspecloc) = floats(3+nbeg)
    else
      call outerror('Incorrect coefficient input for EAM functional',iline)
      call stopnow('potwordm')
    endif
    if (abs(eamfnpar(3,neamfnspecloc)).lt.1.0d-12) then
      call outerror('equilibrium density too close to zero in EAM functional',iline)
      call stopnow('potwordm')
    endif
  elseif (neamfn.eq.10) then
!********
!  VBO  *
!********
!
!  Fitting flags
!
    if (lfit.and.lfflags) then
      n1 = int(floats(nfloat-1))
      n2 = int(floats(nfloat))
      nfloat = nfloat - 2
      if (n1.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 4
      endif
      if (n2.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 5
      endif
    endif
!
!  If number of floats is greater than 2, assume that fitting flags have been left on line
!
    if (nfloat.gt.2) nfloat = nfloat - 2
!
!  Assign coefficients 
!
    if (nfloat.ge.2) then
      eamfnpar(1,neamfnspecloc) = floats(1+nbeg)
      eamfnpar(2,neamfnspecloc) = abs(floats(2+nbeg))
    else
      call outerror('Incorrect coefficient input for EAM functional',iline)
      call stopnow('potwordm')
    endif
    if (abs(eamfnpar(2,neamfnspecloc)).gt.1.0_dp) then
      call outerror('exponent in EAM functional should be less than 1',iline)
      call stopnow('potwordm')
    endif
  elseif (neamfn.eq.11) then
!*****************
!  Belashchenko  *
!*****************
!
!  eamfnpar  1- 7  => rho1 - rho7
!  eamfnpar  8-17  => a1 - a8
!  eamfnpar 18-23  => b1 - b8
!  eamfnpar 24-31  => c1 - c8
!  eamfnpar 32-33  => m / n
!
!  Assign density cutoffs
!
    if (nfloat.ge.3) then
      eamfnpar(1,neamfnspecloc) = floats(1+nbeg)
      eamfnpar(2,neamfnspecloc) = floats(2+nbeg)
      eamfnpar(3,neamfnspecloc) = floats(3+nbeg)
    else
      call outerror('Incorrect coefficient input for EAM functional',iline)
      call stopnow('potwordm')
    endif
!
!  Read second line
!
    line = '  '
    read(iin,'(a)',end=118) line
    iline = iline + 1
    call linepro(iin,line,iline)
!
!  Assign coefficients 
!
    if (nfloat.ge.4) then
      eamfnpar(4,neamfnspecloc) = floats(1+nbeg)
      eamfnpar(5,neamfnspecloc) = floats(2+nbeg)
      eamfnpar(6,neamfnspecloc) = floats(3+nbeg)
      eamfnpar(7,neamfnspecloc) = floats(4+nbeg)
    else
      call outerror('Incorrect coefficient input for EAM functional',iline)
      call stopnow('potwordm')
    endif
!
!  Read third line
!
    line = '  '
    read(iin,'(a)',end=118) line
    iline = iline + 1
    call linepro(iin,line,iline)
!
!  Assign coefficients
!
    if (nfloat.ge.3) then
      eamfnpar(8,neamfnspecloc) = floats(1+nbeg)
      eamfnpar(32,neamfnspecloc) = floats(2+nbeg)
      eamfnpar(33,neamfnspecloc) = floats(3+nbeg)
    else
      call outerror('Incorrect coefficient input for EAM functional',iline)
      call stopnow('potwordm')
    endif
!
!  Read fourth line
!
    line = '  '
    read(iin,'(a)',end=118) line
    iline = iline + 1
    call linepro(iin,line,iline)
!
!  Assign coefficients 
!
    if (nfloat.ge.4) then
      eamfnpar(24,neamfnspecloc) = floats(1+nbeg)
      eamfnpar(25,neamfnspecloc) = floats(2+nbeg)
      eamfnpar(26,neamfnspecloc) = floats(3+nbeg)
      eamfnpar(27,neamfnspecloc) = floats(4+nbeg)
    else
      call outerror('Incorrect coefficient input for EAM functional',iline)
      call stopnow('potwordm')
    endif
!
!  Read fifth line
!
    line = '  '
    read(iin,'(a)',end=118) line
    iline = iline + 1
    call linepro(iin,line,iline)
!
!  Assign coefficients 
!
    if (nfloat.ge.4) then
      eamfnpar(28,neamfnspecloc) = floats(1+nbeg)
      eamfnpar(29,neamfnspecloc) = floats(2+nbeg)
      eamfnpar(30,neamfnspecloc) = floats(3+nbeg)
      eamfnpar(31,neamfnspecloc) = floats(4+nbeg)
    else
      call outerror('Incorrect coefficient input for EAM functional',iline)
      call stopnow('potwordm')
    endif
!
!  Compute dependent parameters
!
    call setbelashchenko(neamfnspecloc)
  elseif (neamfn.eq.12) then
!*************
!  Igarashi  *
!*************
!
!  Fitting flags
!
    if (lfit.and.lfflags) then
      n1 = int(floats(nfloat-1))
      n2 = int(floats(nfloat))
      nfloat = nfloat - 2
      if (n1.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 4
      endif
      if (n2.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 5
      endif
    endif
!
!  If number of floats is greater than 2, assume that fitting flags have been left on line
!
    if (nfloat.gt.2) nfloat = nfloat - 2
!
!  Assign coefficients
!
    if (nfloat.ge.2) then
      eamfnpar(1,neamfnspecloc) = floats(1+nbeg)
      eamfnpar(2,neamfnspecloc) = floats(2+nbeg)
    else
      call outerror('Incorrect coefficient input for EAM functional',iline)
      call stopnow('potwordm')
    endif
  elseif (neamfn.eq.13) then
!***********
!  Spline  *
!***********
!
!  Fitting flags
!
    if (lfit.and.lfflags) then
      n1 = int(floats(nfloat-3))
      n2 = int(floats(nfloat-2))
      n3 = int(floats(nfloat-1))
      n4 = int(floats(nfloat))
      nfloat = nfloat - 4
      if (n1.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 4
      endif
      if (n2.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 5
      endif
      if (n3.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 6
      endif
      if (n4.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 9
      endif
    endif
!
!  If number of floats is greater than 6, assume that fitting flags have been left on line
!
    if (nfloat.gt.6) nfloat = nfloat - 4
!
!  Assign coefficients 
!
    if (nfloat.ge.6) then
      eamfnpar(1,neamfnspecloc) = floats(1+nbeg)
      eamfnpar(2,neamfnspecloc) = floats(2+nbeg)
      eamfnpar(3,neamfnspecloc) = floats(3+nbeg)
      eamfnpar(4,neamfnspecloc) = floats(4+nbeg)
      eamfnpar(5,neamfnspecloc) = abs(floats(5+nbeg))
      eamfnpar(6,neamfnspecloc) = abs(floats(6+nbeg))
    else
      call outerror('Incorrect coefficient input for EAM functional',iline)
      call stopnow('potwordm')
    endif
    if (eamfnpar(5,neamfnspecloc).gt.eamfnpar(6,neamfnspecloc)) then
      call outerror('rho0 is greater than rho_max in EAM functional',iline)
      call stopnow('potwordm')
    endif
  else
!***************************
!  Default EAM functional  *
!***************************
!
!  Fitting flags
!
    if (lfit.and.lfflags) then
      n1 = int(floats(nfloat))
      nfloat = nfloat - 1
      if (n1.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 4
      endif
    endif
!
!  If number of floats is greater than 1, assume that fitting flags have been left on line
!
    if (nfloat.gt.1) nfloat = nfloat - 1
!
!  Assign coefficients 
!
    if (nfloat.ge.1) then
      eamfnpar(1,neamfnspecloc) = floats(1+nbeg)
    else
      eamfnpar(1,neamfnspecloc) = 1.0_dp
    endif
  endif
  goto 115
118 if (.not.l55) l1000 = .true.
  return
119 call outerror('file for numerical EAM cannot be opened',iline)
  call stopnow('potwordm')
!********************************************
!  Parameters for EAM alloy transformations *
!********************************************
120 continue
!
!  Read in species parameters for functional
!
125 line = '  '
  read(iin,'(a)',end=128) line
  iline = iline + 1
  call linepro(iin,line,iline)
!
!  Check whether input is symbol or option
!
  if ((nword+nfloat).eq.0) goto 125
  if (nword.gt.0) then
    word = words(1)
    call worsy(word,lsymbol,.true.)
    if (.not.lsymbol) then
      l55 = .true.
      goto 128
    endif
  endif
!
!  Find atomic species
!
  ltype01 = .false.
  sym1 = ' '
  if (nword.gt.0) then
    if (llibrary) then
      call okspec(lok1,words(1),ilp1,.true.,ltype01)
      if (.not.lok1) then
        goto 125
      endif
    endif
    call ltont(words(1),nvar1,itype1)
    sym1 = words(1)(1:5)
    if (ltype01) itype1 = 0
    if (nword.ge.2) then
      if (index(words(2),'she').eq.1) then
        nvar1 = nvar1 + maxele
      endif
    endif
    nbeg = 0
  else
!
!  Numeric input
!
    nvar1 = int(floats(1))
    if (nvar1.gt.100) nvar1 = nvar1 - 100 + maxele
    itype1 = 0
    nbeg = 1
    nfloat = nfloat - 1
    call label(nvar1,itype1,sym1)
  endif
!
!  Check to see whether this atomic species is already in the list and if not add it
!
  i = 0 
  lfound = .false.
  do while (i.lt.neamspec.and..not.lfound)
    i = i + 1
    lfound = (neamnat(i).eq.nvar1.and.neamtyp(i).eq.itype1)
  enddo
  if (lfound) then
    neamspecloc = i
  else
    neamspec = neamspec + 1
    neamspecloc = neamspec
    if (neamspec.gt.maxeamspec) then
      maxeamspec = neamspec + 20
      call changemaxeamspec
    endif
!
!  When adding a new species initialize related variables
!
    ndenfncomp(neamspec) = 0
!
    neamnat(neamspec) = nvar1
    neamtyp(neamspec) = itype1
    symboleamspec(1,neamspec) = sym1
  endif
!
!  Fitting flags
!
  if (lfit.and.lfflags) then
    n1 = int(floats(nfloat))
    nfloat = nfloat - 1
    if (n1.eq.1) then
      nfit = nfit + 1
      if (nfit.ge.maxfit) then
        maxfit = nfit + 10
        call changemaxfit
      endif
      nftyp(nfit) = 5
      nfpot(nfit) = neamspecloc
      nfvar(nfit) = 7
    endif
  endif
!
!  If number of floats is greater than 2, assume that fitting flags have been left on line
!
  if (nfloat.gt.1) nfloat = nfloat - 1
!
!  Assign coefficients 
!
  if (nfloat.ge.1) then
    eamalloy(1,neamspecloc) = floats(1+nbeg)
    eamalloy(2,neamspecloc) = 0.0_dp
  else
    call outerror('Incorrect coefficient input for EAM alloy option',iline)
    call stopnow('potwordm')
  endif
  goto 125
128 if (.not.l55) l1000 = .true.
  lwordok = .true.
  return
!**********************************
!  Parameters for MEAM densities  *
!**********************************
130 units = 1.0_dp
  ndenfnloc = 1
  if (nword.gt.1) then
    i = 2
    do while (i.le.nword)
      if (index(words(i),'exp').eq.1) then
        ndenfnloc = 2
      elseif (index(words(i),'gau').eq.1) then
        ndenfnloc = 3
      elseif (index(words(i),'pow').eq.1) then
        ndenfnloc = 1
      elseif (index(words(i),'cub').eq.1) then
        ndenfnloc = 4
      elseif (index(words(i),'vot').eq.1) then
        ndenfnloc = 5
      elseif (index(words(i),'quad').eq.1) then
        ndenfnloc = 6
      elseif (index(words(i),'quar').eq.1) then
        ndenfnloc = 7
      elseif (index(words(i),'glue').eq.1) then
        ndenfnloc = 8
      elseif (index(words(i),'evo').eq.1) then
        ndenfnloc = 9
      elseif (index(words(i),'mei').eq.1) then
        ndenfnloc = 10
      elseif (index(words(i),'num').eq.1) then
        ndenfnloc = 11
      elseif (index(words(i),'bas').eq.1) then
        ndenfnloc = 12
      elseif (index(words(i),'vbo').eq.1) then
        ndenfnloc = 13
      elseif (index(words(i),'fpo').eq.1) then
        ndenfnloc = 14
      elseif (index(words(i),'spl').eq.1) then
        ndenfnloc = 15
      endif
      i = i + 1
    enddo
  endif
  if (nfloat.gt.0) then
!
!  Set MEAM order
!
    meamorder = nint(floats(1)) + 1
!
!  Check that order is in allowed range
!
    if (meamorder.lt.1.or.meamorder.gt.maxmeamorder) then
      call outerror('Order of MEAM density is invalid',iline)
      call stopnow('potwordm')
    endif
  else
!
!  Set default MEAM order
!
    meamorder = 4
  endif
!
!  Set global flag as to whether this is really a MEAM density or not
!
  lMEAMden = .true.
  if (lMEAMden.and.lMEAMfn) lMEAM = .true.
!
!  If this is MEAM then disable third derivatives
!
  lnoanald3 = .true.
!
  if (nfloat.gt.1) then
    denpar6loc = nint(floats(2))
  elseif (ndenfnloc.eq.1) then
    denpar6loc = 6
  else
    denpar6loc = 0
  endif
135 line = '  '
  read(iin,'(a)',end=138) line
  iline = iline + 1
  call linepro(iin,line,iline)
!
!  Check whether input is symbol or option
!
  if ((nword+nfloat).eq.0) goto 135
  if (nword.gt.0) then
    word = words(1)
    call worsy(word,lsymbol,.true.)
    if (.not.lsymbol) then
      l55 = .true.
      goto 138
    endif
  endif
!
!  Symbols used in input
!
  nvar2 = 0
  itype2 = 0
  ltype01 = .false.
  sym1 = ' '
  sym2 = ' '
  if (nword.gt.0) then
    sym1 = words(1)(1:5)
    if (llibrary) then
      call okspec(lok1,words(1),ilp1,.true.,ltype01)
      if (.not.lok1) then
!
!  Read remaining potential lines before trying next potential
!
        do order = 2,meamorder
          read(iin,'(a)',end=138) line
          iline = iline + 1
        enddo
        goto 135
      endif
    endif
    call ltont(words(1),nvar1,itype1)
    if (ltype01) itype1 = 0
    if (nword.ge.4) then
      if (index(words(2),'she').eq.1) then
        nvar1 = nvar1 + maxele
      endif
      call ltont(words(3),nvar2,itype2)
      sym2 = words(3)(1:5)
      if (index(words(4),'she').eq.1) then
        nvar2 = nvar2 + maxele
      endif
    elseif (nword.eq.3) then
      if (index(words(2),'she').eq.1) then
        nvar1 = nvar1 + maxele
      elseif (index(words(2),'cor').eq.0) then
        call ltont(words(2),nvar2,itype2)
        sym2 = words(2)(1:5)
        if (index(words(3),'she').eq.1) then
          nvar2 = nvar2 + maxele
        endif
      else
        call ltont(words(3),nvar2,itype2)
      endif
    elseif (nword.eq.2) then
      if (index(words(2),'she').eq.1) then
        nvar1 = nvar1 + maxele
      elseif (index(words(2),'cor').eq.0) then
        call ltont(words(2),nvar2,itype2)
        sym2 = words(2)(1:5)
      endif
    endif
    nbeg = 0
  else
!
!  Numeric input
!
    nvar1 = int(floats(1))
    if (nvar1.gt.100) nvar1 = nvar1 - 100 + maxele
    itype1 = 0
    nbeg = 1
    nfloat = nfloat - 1
    call label(nvar1,itype1,sym1)
  endif
!
!  Check to see whether this atomic species is already in the list and if not add it
!
  i = 0
  lfound = .false.       
  do while (i.lt.neamspec.and..not.lfound)
    i = i + 1
    lfound = (neamnat(i).eq.nvar1.and.neamtyp(i).eq.itype1.and.neamnat2(i).eq.nvar2.and.neamtyp2(i).eq.itype2)
  enddo    
  if (lfound) then
    neamspecloc = i
  else
    neamspec = neamspec + 1
    neamspecloc = neamspec
    if (neamspec.gt.maxeamspec) then
      maxeamspec = neamspec + 20
      call changemaxeamspec
    endif     
    neamnat(neamspec) = nvar1
    neamtyp(neamspec) = itype1
    neamnat2(neamspec) = nvar2
    neamtyp2(neamspec) = itype2
    symboleamspec(1,neamspec) = sym1
    symboleamspec(2,neamspec) = sym2
  endif
!
!  Copy first line info from first density
!
  ndenfncomp(neamspecloc) = ndenfncomp(neamspecloc) + 1
  if (ndenfncomp(neamspecloc).gt.maxeamden) then
    maxeamden = ndenfncomp(neamspecloc) + 3
    call changemaxeamden
  endif
  ndenfn(ndenfncomp(neamspecloc),neamspecloc) = ndenfnloc
  neammeamorder(ndenfncomp(neamspecloc),neamspecloc) = meamorder
!
!  Set species flag as to whether this is really MEAM or not
!
  if (meamorder.gt.1) lmeamspec(neamspecloc) = .true.
!
!  Loop over orders of MEAM density
!
  do order = 1,meamorder
!
!  If not the first line then read new line
!
    if (order.gt.1) then
136   line = '  '
      read(iin,'(a)',end=138) line
      iline = iline + 1
      call linepro(iin,line,iline)
      if ((nword+nfloat).eq.0) goto 136
    endif
!
!  Copy first line info from first density that are specific to each order
!
    denpar(6,order,ndenfncomp(neamspecloc),neamspecloc) = denpar6loc
!
!  Fitting flags
!
    if (lfit.and.lfflags) then
      if (ndenfnloc.eq.2.or.ndenfnloc.eq.3.or.ndenfnloc.eq.12) then
        n1 = int(floats(nfloat-2))
        n2 = int(floats(nfloat-1))
        n3 = int(floats(nfloat))
        nfloat = nfloat - 3
        if (n1.eq.1) then
          nfit = nfit + 1
          if (nfit.ge.maxfit) then
            maxfit = nfit + 10
            call changemaxfit
          endif
          nftyp(nfit) = 5
          nfpot(nfit) = neamspecloc
          nfvar(nfit) = 1
          nfvar2(nfit) = ndenfncomp(neamspecloc)
          nfvar3(nfit) = order
        endif
        if (n2.eq.1) then
          nfit = nfit + 1
          if (nfit.ge.maxfit) then
            maxfit = nfit + 10
            call changemaxfit
          endif
          nftyp(nfit) = 5
          nfpot(nfit) = neamspecloc
          nfvar(nfit) = 2
          nfvar2(nfit) = ndenfncomp(neamspecloc)
          nfvar3(nfit) = order
        endif
        if (n3.eq.1) then
          nfit = nfit + 1
          if (nfit.ge.maxfit) then
            maxfit = nfit + 10
            call changemaxfit
          endif
          nftyp(nfit) = 5
          nfpot(nfit) = neamspecloc
          nfvar(nfit) = 3
          nfvar2(nfit) = ndenfncomp(neamspecloc)
          nfvar3(nfit) = order
        endif
      elseif (ndenfnloc.eq.1) then
        n1 = int(floats(nfloat))
        nfloat = nfloat - 1
        if (n1.eq.1) then
          nfit = nfit + 1
          if (nfit.ge.maxfit) then
            maxfit = nfit + 10
            call changemaxfit
          endif
          nftyp(nfit) = 5
          nfpot(nfit) = neamspecloc
          nfvar(nfit) = 1
          nfvar2(nfit) = ndenfncomp(neamspecloc)
          nfvar3(nfit) = order
        endif
      elseif (ndenfnloc.eq.4.or.ndenfnloc.eq.5.or.ndenfnloc.eq.6.or.ndenfnloc.eq.7.or.ndenfnloc.eq.9.or.ndenfnloc.eq.14) then
        n1 = int(floats(nfloat-1))
        n2 = int(floats(nfloat))
        nfloat = nfloat - 2
        if (n1.eq.1) then
          nfit = nfit + 1
          if (nfit.ge.maxfit) then
            maxfit = nfit + 10
            call changemaxfit
          endif
          nftyp(nfit) = 5
          nfpot(nfit) = neamspecloc
          nfvar(nfit) = 1
          nfvar2(nfit) = ndenfncomp(neamspecloc)
          nfvar3(nfit) = order
        endif
        if (n2.eq.1) then
          nfit = nfit + 1
          if (nfit.ge.maxfit) then
            maxfit = nfit + 10
            call changemaxfit
          endif
          nftyp(nfit) = 5
          nfpot(nfit) = neamspecloc
          nfvar(nfit) = 2
          nfvar2(nfit) = ndenfncomp(neamspecloc)
          nfvar3(nfit) = order
        endif
      elseif (ndenfnloc.eq.10) then
        do n1 = 1,7
          n2 = int(floats(nfloat+n1-7))
          if (n2.eq.1) then 
            nfit = nfit + 1
            if (nfit.ge.maxfit) then
              maxfit = nfit + 10
              call changemaxfit
            endif
            nftyp(nfit) = 5
            nfpot(nfit) = neamspecloc
            if (n1.le.3) then
              nfvar(nfit) = n1
            else
              nfvar(nfit) = n1 + 11
            endif
            nfvar2(nfit) = ndenfncomp(neamspecloc)
            nfvar3(nfit) = order
          endif
        enddo
        nfloat = nfloat - 8
      elseif (ndenfnloc.eq.13) then
        do n1 = 1,5
          n2 = int(floats(nfloat+n1-5))
          if (n2.eq.1) then
            nfit = nfit + 1
            if (nfit.ge.maxfit) then
              maxfit = nfit + 10
              call changemaxfit
            endif
            nftyp(nfit) = 5
            nfpot(nfit) = neamspecloc
            if (n1.le.3) then
              nfvar(nfit) = n1
            else
              nfvar(nfit) = n1 + 11
            endif
            nfvar2(nfit) = ndenfncomp(neamspecloc)
            nfvar3(nfit) = 1
          endif
        enddo
        nfloat = nfloat - 5
      elseif (ndenfnloc.eq.15) then
        do n1 = 1,4
          n2 = int(floats(nfloat+n1-4))
          if (n2.eq.1) then
            nfit = nfit + 1
            if (nfit.ge.maxfit) then
              maxfit = nfit + 10
              call changemaxfit
            endif
            nftyp(nfit) = 5
            nfpot(nfit) = neamspecloc
            if (n1.le.3) then
              nfvar(nfit) = n1
            else
              nfvar(nfit) = n1 + 11
            endif
            nfvar2(nfit) = ndenfncomp(neamspecloc)
            nfvar3(nfit) = 1
          endif
        enddo
        nfloat = nfloat - 4
      endif
    endif
    if (ndenfnloc.eq.1) then
!
!  If number of floats is greater than 2, assume that fitting flags have been left on line
!
      if (nfloat.gt.2) nfloat = nfloat - 2
!
!  Assign coefficients 
!
      if (nfloat.ge.1) then
        denpar(1,order,ndenfncomp(neamspecloc),neamspecloc) = floats(1+nbeg)*units
      else
        call outerror('Incorrect coefficient input for MEAM density',iline)
        call stopnow('potwordm')
      endif
    elseif (ndenfnloc.eq.4.or.ndenfnloc.eq.5.or.ndenfnloc.eq.6.or.ndenfnloc.eq.7.or.ndenfnloc.eq.9.or.ndenfnloc.eq.14) then
!
!  If number of floats is greater than 3, assume that fitting flags have been left on line
!
      if (nfloat.gt.3) nfloat = nfloat - 3
!
!  Assign coefficients 
!
      if (nfloat.ge.2) then
        denpar(1,order,ndenfncomp(neamspecloc),neamspecloc) = floats(1+nbeg)*units
        denpar(2,order,ndenfncomp(neamspecloc),neamspecloc) = floats(2+nbeg)
      elseif (nfloat.eq.1) then
        denpar(1,order,ndenfncomp(neamspecloc),neamspecloc) = floats(1+nbeg)*units
        denpar(2,order,ndenfncomp(neamspecloc),neamspecloc) = 0.0_dp
      else
        call outerror('Incorrect coefficient input for MEAM density',iline)
        call stopnow('potwordm')
      endif
    elseif (ndenfnloc.eq.10) then
!
!  If number of floats is greater than 8, assume that fitting flags have been left on line
!
      if (nfloat.gt.8) nfloat = nfloat - 8
!
!  Assign coefficients 
!
      if (nfloat.ge.7) then
        denpar(1,order,ndenfncomp(neamspecloc),neamspecloc) = floats(1+nbeg)*units
        denpar(2,order,ndenfncomp(neamspecloc),neamspecloc) = floats(2+nbeg)*units
        denpar(3,order,ndenfncomp(neamspecloc),neamspecloc) = floats(3+nbeg)*units
        denpar(4,order,ndenfncomp(neamspecloc),neamspecloc) = floats(4+nbeg)*units
        denpar(5,order,ndenfncomp(neamspecloc),neamspecloc) = floats(5+nbeg)*units
        denpar(6,order,ndenfncomp(neamspecloc),neamspecloc) = floats(6+nbeg)*units
        denpar(7,order,ndenfncomp(neamspecloc),neamspecloc) = floats(7+nbeg)
      else
        call outerror('Incorrect coefficient input for MEAM density',iline)
        call stopnow('potwordm')
      endif
      if (abs(denpar(7,order,ndenfncomp(neamspecloc),neamspecloc)).lt.1.0d-12) then
        call outerror('r0 is too close to zero for MEAM density',iline)
        call stopnow('potwordm')
      endif
    elseif (ndenfnloc.eq.13) then
!
!  VBO density
!
!  If number of floats is greater than 5, assume that fitting flags have been left on line
!
      if (nfloat.gt.5) nfloat = nfloat - 5
!
!  Assign coefficients
!
      if (nfloat.ge.5) then
        denpar(1,order,ndenfncomp(neamspecloc),neamspecloc) = floats(1+nbeg)*units
        denpar(2,order,ndenfncomp(neamspecloc),neamspecloc) = floats(2+nbeg)
        denpar(3,order,ndenfncomp(neamspecloc),neamspecloc) = floats(3+nbeg)
        denpar(4,order,ndenfncomp(neamspecloc),neamspecloc) = abs(floats(4+nbeg))
        denpar(5,order,ndenfncomp(neamspecloc),neamspecloc) = abs(floats(5+nbeg))
      else
        call outerror('Incorrect coefficient input for MEAM density',iline)
        call stopnow('potwordm')
      endif
      if (abs(denpar(4,order,ndenfncomp(neamspecloc),neamspecloc)).lt.1.0d-12) then
        call outerror('r0 is too close to zero for MEAM density',iline)
        call stopnow('potwordm')
      endif
      if (abs(denpar(5,order,ndenfncomp(neamspecloc),neamspecloc)).lt.1.0d-12) then
        call outerror('delta is too close to zero for MEAM density',iline)
        call stopnow('potwordm')
      endif
    elseif (ndenfnloc.eq.15) then
!
!  Spline density
!
!  If number of floats is greater than 6, assume that fitting flags have been left on line
!
      if (nfloat.gt.6) nfloat = nfloat - 4
!
!  Assign coefficients
!
      if (nfloat.ge.6) then
        denpar(1,order,ndenfncomp(neamspecloc),neamspecloc) = floats(1+nbeg)*units
        denpar(2,order,ndenfncomp(neamspecloc),neamspecloc) = floats(2+nbeg)*units
        denpar(3,order,ndenfncomp(neamspecloc),neamspecloc) = floats(3+nbeg)*units
        denpar(4,order,ndenfncomp(neamspecloc),neamspecloc) = floats(4+nbeg)*units
        denpar(5,order,ndenfncomp(neamspecloc),neamspecloc) = abs(floats(5+nbeg))
        denpar(6,order,ndenfncomp(neamspecloc),neamspecloc) = abs(floats(6+nbeg))
      else
        call outerror('Incorrect coefficient input for MEAM density',iline)
        call stopnow('potwordm')
      endif
      if (abs(denpar(5,order,ndenfncomp(neamspecloc),neamspecloc)).lt.1.0d-12) then
        call outerror('rmin is too close to zero for MEAM density',iline)
        call stopnow('potwordm')
      endif
      if (abs(denpar(6,order,ndenfncomp(neamspecloc),neamspecloc)).lt.1.0d-12) then
        call outerror('r0 is too close to zero for MEAM density',iline)
        call stopnow('potwordm')
      endif
    elseif (ndenfnloc.eq.8) then
!
!  Assign distance cutoffs : d / r_b / r_m
!
      if (nfloat.ge.3) then
        denpar(1,order,ndenfncomp(neamspecloc),neamspecloc) = floats(1+nbeg)
        denpar(2,order,ndenfncomp(neamspecloc),neamspecloc) = floats(2+nbeg)
        denpar(3,order,ndenfncomp(neamspecloc),neamspecloc) = floats(3+nbeg)
      else
        call outerror('Incorrect coefficient input for MEAM density',iline)
        call stopnow('potwordm')
      endif
!
!  Read second line
!
      line = '  '
      read(iin,'(a)',end=138) line
      iline = iline + 1
      call linepro(iin,line,iline)
!
!  Assign coefficients 
!
      if (nfloat.ge.4) then
        denpar(4,order,ndenfncomp(neamspecloc),neamspecloc) = floats(1+nbeg)*units
        denpar(5,order,ndenfncomp(neamspecloc),neamspecloc) = floats(2+nbeg)*units
        denpar(6,order,ndenfncomp(neamspecloc),neamspecloc) = floats(3+nbeg)*units
        denpar(7,order,ndenfncomp(neamspecloc),neamspecloc) = floats(4+nbeg)*units
      else
        call outerror('Incorrect coefficient input for MEAM density',iline)
        call stopnow('potwordm')
      endif
!
!  Read third line
!
      line = '  '
      read(iin,'(a)',end=138) line
      iline = iline + 1
      call linepro(iin,line,iline)
!
!  Assign coefficients 
!
      if (nfloat.ge.4) then
        denpar(8,order,ndenfncomp(neamspecloc),neamspecloc)  = floats(1+nbeg)*units
        denpar(9,order,ndenfncomp(neamspecloc),neamspecloc)  = floats(2+nbeg)*units
        denpar(10,order,ndenfncomp(neamspecloc),neamspecloc) = floats(3+nbeg)*units
        denpar(11,order,ndenfncomp(neamspecloc),neamspecloc) = floats(4+nbeg)*units
      else
        call outerror('Incorrect coefficient input for MEAM density',iline)
        call stopnow('potwordm')
      endif
!
!  Read fourth line
!
      line = '  '
      read(iin,'(a)',end=138) line
      iline = iline + 1
      call linepro(iin,line,iline)
!
!  Assign coefficients 
!
      if (nfloat.ge.4) then
        denpar(12,order,ndenfncomp(neamspecloc),neamspecloc) = floats(1+nbeg)*units
        denpar(13,order,ndenfncomp(neamspecloc),neamspecloc) = floats(2+nbeg)*units
        denpar(14,order,ndenfncomp(neamspecloc),neamspecloc) = floats(3+nbeg)*units
        denpar(15,order,ndenfncomp(neamspecloc),neamspecloc) = floats(4+nbeg)*units
      else
        call outerror('Incorrect coefficient input for MEAM density',iline)
        call stopnow('potwordm')
      endif
    else
!
!  If number of floats is greater than 4, assume that fitting flags have been left on line
!
      if (nfloat.gt.4) nfloat = nfloat - 4
!
!  Assign coefficients 
!
      if (nfloat.ge.3) then
        denpar(1,order,ndenfncomp(neamspecloc),neamspecloc) = floats(1+nbeg)*units
        denpar(2,order,ndenfncomp(neamspecloc),neamspecloc) = floats(2+nbeg)
        denpar(3,order,ndenfncomp(neamspecloc),neamspecloc) = floats(3+nbeg)
      elseif (nfloat.eq.2) then
        denpar(1,order,ndenfncomp(neamspecloc),neamspecloc) = floats(1+nbeg)*units
        denpar(2,order,ndenfncomp(neamspecloc),neamspecloc) = floats(2+nbeg)
        denpar(3,order,ndenfncomp(neamspecloc),neamspecloc) = 0.0_dp
      else
        call outerror('Incorrect coefficient input for MEAM density',iline)
        call stopnow('potwordm')
      endif
    endif
!
!  End of loop over MEAM order
!
  enddo
!
  goto 135
138 if (.not.l55) l1000 = .true.
  lwordok = .true.
  return
!************************************
!  Parameters for MEAM functionals  *
!************************************
140 units = 1.0_dp
  meamtype = 2
  if (nword.gt.1) then
    i = 2
    do while (i.le.nword)
      if (index(words(i),'sq').eq.1) then
        neamfn = 1
        neampower = 2
      elseif (index(words(i),'po').eq.1) then
        neamfn = 2
      elseif (index(words(i),'ban').eq.1) then
        neamfn = 3
      elseif (index(words(i),'nu').eq.1) then
        neamfn = 4
      elseif (index(words(i),'jo').eq.1) then
        neamfn = 5
      elseif (index(words(i),'gl').eq.1) then
        neamfn = 6
      elseif (index(words(i),'fo').eq.1) then
        neamfn = 7
      elseif (index(words(i),'me').eq.1) then
        neamfn = 8
      elseif (index(words(i),'bas').eq.1) then
        neamfn = 9
      elseif (index(words(i),'vb').eq.1) then
        neamfn = 10
      elseif (index(words(i),'sp').eq.1) then
        neamfn = 13
      elseif (index(words(i),'au').eq.1) then
        units = autoev
      elseif (index(words(i),'kc').eq.1) then
        units = kcaltoev
      elseif (index(words(i),'kj').eq.1) then
        units = kjmtoev
      endif
      i = i + 1
    enddo
  endif
  if (nfloat.gt.0) then
!
!  Set MEAM order
!
    meamorder = nint(floats(1)) + 1
!
!  Check that order is in allowed range
!
    if (meamorder.lt.1.or.meamorder.gt.maxmeamorder) then
      call outerror('Order of MEAM functional is invalid',iline)
      call stopnow('potwordm')
    endif
  else
!
!  Set default MEAM order
!
    meamorder = 4
  endif
  if (nfloat.gt.1.and.neamfn.gt.1) then
    neampower = nint(floats(2))
  endif
  if (neamfn.gt.1.and.neampower.eq.0) then
    call outerror('inverse power of zero in MEAM functional',iline)
    call stopnow('potwordm')
  endif
!
!  Read in species parameters for functional
!
145 line = '  '
  read(iin,'(a)',end=148) line
  iline = iline + 1
  call linepro(iin,line,iline)
!
!  Check whether input is symbol or option
!
  if ((nword+nfloat).eq.0) goto 145
  if (nword.gt.0) then
    word = words(1)
    call worsy(word,lsymbol,.true.)
    if (.not.lsymbol) then
      l55 = .true.
      goto 148
    endif
  endif
!
!  Symbol used in input
!
  ltype01 = .false.
  sym1 = ' '
  if (nword.gt.0) then
    sym1 = words(1)(1:5)
    if (llibrary) then
      call okspec(lok1,words(1),ilp1,.true.,ltype01)
      if (.not.lok1) then
!
!  Read coefficient line if needed before checking next symbol
!
        line = '  '
        read(iin,'(a)',end=148) line
        iline = iline + 1
        goto 145
      endif
      if (ilp1.gt.0.and.ilp1.le.nspec) then
        nvar1 = natspec(ilp1)
        if (ltype01) then
          itype1 = 0
        else
          itype1 = ntypspec(ilp1)
        endif
      elseif (ilp1.eq.-1) then
        nvar1 = maxele
        itype1 = 0
      endif
    else
      call ltont(words(1),nvar1,itype1)
    endif
    if (nword.ge.2) then
      if (index(words(2),'she').eq.1) then
        nvar1 = nvar1 + maxele
      endif
    endif
    nbeg = 0
  else
!
!  Numeric input
!
    nvar1 = int(floats(1))
    if (nvar1.gt.100) nvar1 = nvar1 - 100 + maxele
    itype1 = 0
    nbeg = 1
    nfloat = nfloat - 1
    call label(nvar1,itype1,sym1)
  endif
!
!  Check to see whether this atomic species is already in the list and if not add it
!
  i = 0
  lfound = .false.       
  do while (i.lt.neamfnspec.and..not.lfound)
    i = i + 1
    lfound = (neamfnnat(i).eq.nvar1.and.neamfntyp(i).eq.itype1)
  enddo    
  if (lfound) then
    neamfnspecloc = i
  else
    neamfnspec = neamfnspec + 1
    neamfnspecloc = neamfnspec
    if (neamfnspec.gt.maxeamfnspec) then
      maxeamfnspec = neamfnspec + 20
      call changemaxeamfnspec
    endif     
    neamfnnat(neamfnspec) = nvar1
    neamfntyp(neamfnspec) = itype1
    symboleamfnspec(neamfnspec) = sym1
  endif
!
!  Set MEAM order
!
  neamfnmeamorder(neamfnspec) = meamorder
!
!  Set global flag as to whether this is really a MEAM density or not
!
  lMEAMfn = .true.
  if (lMEAMden.and.lMEAMfn) lMEAM = .true.
!
!  If this is MEAM then disable third derivatives
!
  lnoanald3 = .true.
!
  lwordok = .true.
  if (neamfn.eq.3) then
!***********************
!  Banerjea and Smith  *
!***********************
!
!  Fitting flags
!
    if (lfit.and.lfflags) then
      n1 = int(floats(nfloat-2))
      n2 = int(floats(nfloat-1))
      n3 = int(floats(nfloat))
      nfloat = nfloat - 3
      if (n1.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 4
      endif
      if (n2.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 5
      endif
      if (n3.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 6
      endif
    endif
!
!  If number of floats is greater than 3, assume that fitting flags have been left on line
!
    if (nfloat.gt.3) nfloat = nfloat - 3
!
!  Assign coefficients 
!
    if (nfloat.ge.3) then
      eamfnpar(1,neamfnspecloc) = floats(1+nbeg)
      eamfnpar(2,neamfnspecloc) = floats(2+nbeg)
      eamfnpar(3,neamfnspecloc) = floats(3+nbeg)
    elseif (nfloat.eq.2) then
      eamfnpar(1,neamfnspecloc) = floats(1+nbeg)
      eamfnpar(2,neamfnspecloc) = floats(2+nbeg)
      eamfnpar(3,neamfnspecloc) = 0.0_dp
    else
      call outerror('Incorrect coefficient input for MEAM functional',iline)
      call stopnow('potwordm')
    endif
    if (abs(eamfnpar(3,neamfnspecloc)).lt.1.0d-12) then
      call outerror('equilibrium density too close to zero in MEAM functional',iline)
      call stopnow('potwordm')
    endif
  elseif (neamfn.eq.4) then
!******************************
!  Numerical MEAM functional  *
!******************************
!
!  Is second word core or shel? If so skip
!
    iword = 2
    if (index(words(2),'cor').eq.1.or.index(words(2),'she').eq.1) iword = 3
    if (iword.gt.nword) then
      call outerror('filename missing for numerical MEAM',iline)
      call stopnow('potwordm')
    endif
!
!  Assign file name 
!
    eamfnfile(neamfnspecloc) = words(iword)
!
!  Check that file exists by opening and reading
!
    open(9,file=eamfnfile(neamfnspecloc),status='old',form='formatted',err=149)
    read(9,'(a)') line
    read(9,'(a)') line
    read(9,'(i5,e24.16)') neamfnnumeric(neamfnspecloc),eamfnnumericdrho(neamfnspecloc)
!
!  Check the size of the array for numeric function
!
    if (neamfnnumeric(neamfnspecloc).gt.maxneamfnnumeric) then
      maxneamfnnumeric = neamfnnumeric(neamfnspecloc)
      call changemaxneamfnnumeric
    endif
!
!  Read embedding function data
!
    read(9,'(5e24.16)') (eamfnnumeric(j,neamfnspecloc),j=1,neamfnnumeric(neamfnspecloc))
!
!  Close file having read all relevant information 
!
    close(9)
!
!  Interpolation setup
!
    nrho = neamfnnumeric(neamfnspecloc)
    eamfnnumeric1(1,neamfnspecloc) = eamfnnumeric(2,neamfnspecloc) - eamfnnumeric(1,neamfnspecloc)
    eamfnnumeric1(2,neamfnspecloc) = 0.5_dp*(eamfnnumeric(3,neamfnspecloc) - eamfnnumeric(1,neamfnspecloc))
    eamfnnumeric1(nrho-1,neamfnspecloc) = &
      0.5_dp*(eamfnnumeric(nrho,neamfnspecloc) - eamfnnumeric(nrho-2,neamfnspecloc))
    eamfnnumeric1(nrho,neamfnspecloc) = eamfnnumeric(nrho,neamfnspecloc) - eamfnnumeric(nrho-1,neamfnspecloc)
    do j = 3,nrho-2
      eamfnnumeric1(j,neamfnspecloc) = ((eamfnnumeric(j-2,neamfnspecloc) - eamfnnumeric(j+2,neamfnspecloc)) + &
        8.0_dp*(eamfnnumeric(j+1,neamfnspecloc) - eamfnnumeric(j-1,neamfnspecloc)))/12.0_dp
    enddo 
    do j = 1,nrho-1
      eamfnnumeric2(j,neamfnspecloc) = 3.0_dp*(eamfnnumeric(j+1,neamfnspecloc) - eamfnnumeric(j,neamfnspecloc)) &
        - 2.0_dp*eamfnnumeric1(j,neamfnspecloc) - eamfnnumeric1(j+1,neamfnspecloc)
      eamfnnumeric3(j,neamfnspecloc) = eamfnnumeric1(j,neamfnspecloc) + eamfnnumeric1(j+1,neamfnspecloc) &
        - 2.0_dp*(eamfnnumeric(j+1,neamfnspecloc) - eamfnnumeric(j,neamfnspecloc))
    enddo
    eamfnnumeric2(nrho,neamfnspecloc) = 0.0_dp
    eamfnnumeric3(nrho,neamfnspecloc) = 0.0_dp
    rdrho = 1.0_dp/eamfnnumericdrho(neamfnspecloc)
!
!  Force interpolation set up
!
    do j = 1,nrho
      eamfnnumeric4(j,neamfnspecloc) = eamfnnumeric1(j,neamfnspecloc)*rdrho
      eamfnnumeric5(j,neamfnspecloc) = 2.0_dp*eamfnnumeric2(j,neamfnspecloc)*rdrho
      eamfnnumeric6(j,neamfnspecloc) = 3.0_dp*eamfnnumeric3(j,neamfnspecloc)*rdrho
    enddo
!
!  Second derivative interpolation set up
!
    do j = 1,nrho
      eamfnnumeric7(j,neamfnspecloc) = eamfnnumeric5(j,neamfnspecloc)*rdrho
      eamfnnumeric8(j,neamfnspecloc) = 2.0_dp*eamfnnumeric6(j,neamfnspecloc)*rdrho
    enddo
  elseif (neamfn.eq.5) then
!************
!  Johnson  *
!************
!
!  Fitting flags
!
    if (lfit.and.lfflags) then
      n1 = int(floats(nfloat-5))
      n2 = int(floats(nfloat-4))
      n3 = int(floats(nfloat-3))
      n4 = int(floats(nfloat-2))
      n5 = int(floats(nfloat-1))
      n6 = int(floats(nfloat))
      nfloat = nfloat - 6
      if (n1.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 4
      endif
      if (n2.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 5
      endif
      if (n3.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 6
      endif
      if (n4.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 9
      endif
      if (n5.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 10
      endif
      if (n6.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 11
      endif
    endif
!
!  If number of floats is greater than 6, assume that fitting flags have been left on line
!
    if (nfloat.gt.6) nfloat = nfloat - 6
!
!  Assign coefficients 
!
    if (nfloat.ge.6) then
      eamfnpar(1,neamfnspecloc) = floats(1+nbeg)
      eamfnpar(2,neamfnspecloc) = floats(2+nbeg)
      eamfnpar(3,neamfnspecloc) = floats(3+nbeg)
      eamfnpar(4,neamfnspecloc) = floats(4+nbeg)
      eamfnpar(5,neamfnspecloc) = floats(5+nbeg)
      eamfnpar(6,neamfnspecloc) = floats(6+nbeg)
    elseif (nfloat.eq.5) then
      eamfnpar(1,neamfnspecloc) = floats(1+nbeg)
      eamfnpar(2,neamfnspecloc) = 0.0_dp
      eamfnpar(3,neamfnspecloc) = floats(2+nbeg)
      eamfnpar(4,neamfnspecloc) = floats(3+nbeg)
      eamfnpar(5,neamfnspecloc) = floats(4+nbeg)
      eamfnpar(6,neamfnspecloc) = floats(5+nbeg)
    else
      call outerror('Incorrect coefficient input for MEAM functional',iline)
      call stopnow('potwordm')
    endif
    if (abs(eamfnpar(3,neamfnspecloc)).lt.1.0d-12) then
      call outerror('equilibrium density too close to zero in MEAM functional',iline)
      call stopnow('potwordm')
    endif
    if (abs(eamfnpar(5,neamfnspecloc)).lt.1.0d-12) then
      call outerror('denominator for density power is zero in MEAM functional',iline)
      call stopnow('potwordm')
    endif
  elseif (neamfn.eq.6) then
!*********
!  Glue  *
!*********
!
!  Assign density cutoffs
!
    if (nfloat.ge.2) then
      eamfnpar(1,neamfnspecloc) = floats(1+nbeg)
      eamfnpar(2,neamfnspecloc) = floats(2+nbeg)
    else
      call outerror('Incorrect coefficient input for MEAM functional',iline)
      call stopnow('potwordm')
    endif
!
!  Read second line
!
    line = '  '
    read(iin,'(a)',end=148) line
    iline = iline + 1
    call linepro(iin,line,iline)
!
!  Assign coefficients 
!
    if (nfloat.ge.5) then
      eamfnpar(3,neamfnspecloc) = floats(1+nbeg)
      eamfnpar(4,neamfnspecloc) = floats(2+nbeg)
      eamfnpar(5,neamfnspecloc) = floats(3+nbeg)
      eamfnpar(6,neamfnspecloc) = floats(4+nbeg)
      eamfnpar(7,neamfnspecloc) = floats(5+nbeg)
    else
      call outerror('Incorrect coefficient input for MEAM functional',iline)
      call stopnow('potwordm')
    endif
!
!  Read third line
!
    line = '  '
    read(iin,'(a)',end=148) line
    iline = iline + 1
    call linepro(iin,line,iline)
!
!  Assign coefficients 
!
    if (nfloat.ge.5) then
      eamfnpar(8,neamfnspecloc) = floats(1+nbeg)
      eamfnpar(9,neamfnspecloc) = floats(2+nbeg)
      eamfnpar(10,neamfnspecloc) = floats(3+nbeg)
      eamfnpar(11,neamfnspecloc) = floats(4+nbeg)
      eamfnpar(12,neamfnspecloc) = floats(5+nbeg)
    else
      call outerror('Incorrect coefficient input for MEAM functional',iline)
      call stopnow('potwordm')
    endif
!
!  Read fourth line
!
    line = '  '
    read(iin,'(a)',end=148) line
    iline = iline + 1
    call linepro(iin,line,iline)
!
!  Assign coefficients 
!
    if (nfloat.ge.4) then
      eamfnpar(13,neamfnspecloc) = floats(1+nbeg)
      eamfnpar(14,neamfnspecloc) = floats(2+nbeg)
      eamfnpar(15,neamfnspecloc) = floats(3+nbeg)
      eamfnpar(16,neamfnspecloc) = floats(4+nbeg)
    else
      call outerror('Incorrect coefficient input for MEAM functional',iline)
      call stopnow('potwordm')
    endif
  elseif (neamfn.eq.7) then
!***********
!  Foiles  *
!***********
!
!  Fitting flags
!
    if (lfit.and.lfflags) then
      n1 = int(floats(nfloat-3))
      n2 = int(floats(nfloat-2))
      n3 = int(floats(nfloat-1))
      n4 = int(floats(nfloat))
      nfloat = nfloat - 4
      if (n1.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 4
      endif
      if (n2.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 5
      endif
      if (n3.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 6
      endif
      if (n4.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 9
      endif
    endif
!
!  If number of floats is greater than 4, assume that fitting flags have been left on line
!
    if (nfloat.gt.4) nfloat = nfloat - 4
!
!  Assign coefficients 
!
    if (nfloat.ge.4) then
      eamfnpar(1,neamfnspecloc) = floats(1+nbeg)
      eamfnpar(2,neamfnspecloc) = floats(2+nbeg)
      eamfnpar(3,neamfnspecloc) = floats(3+nbeg)
      eamfnpar(4,neamfnspecloc) = floats(4+nbeg)
    else
      call outerror('Incorrect coefficient input for MEAM functional',iline)
      call stopnow('potwordm')
    endif
  elseif (neamfn.eq.8) then
!******************
!  Mei-Davenport  *
!******************
!
!  Fitting flags
!
    if (lfit.and.lfflags) then
      n1 = int(floats(nfloat-8))
      n2 = int(floats(nfloat-7))
      n3 = int(floats(nfloat-6))
      n4 = int(floats(nfloat-5))
      n5 = int(floats(nfloat-4))
      n6 = int(floats(nfloat-3))
      n7 = int(floats(nfloat-2))
      n8 = int(floats(nfloat-1))
      n9 = int(floats(nfloat))
      nfloat = nfloat - 9
      if (n1.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 4
      endif
      if (n2.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 5
      endif
      if (n3.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 6
      endif
      if (n4.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 9
      endif
      if (n5.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 10
      endif
      if (n6.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 11
      endif
      if (n7.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 12
      endif
      if (n8.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 13
      endif
      if (n9.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 14
      endif
    endif
!
!  If number of floats is greater than 8, assume that fitting flags have been left on line
!
    if (nfloat.gt.9) nfloat = nfloat - 9
!
!  Assign coefficients 
!
    if (nfloat.ge.8) then
      eamfnpar(1,neamfnspecloc) = floats(1+nbeg)*units
      eamfnpar(2,neamfnspecloc) = floats(2+nbeg)
      eamfnpar(3,neamfnspecloc) = floats(3+nbeg)
      eamfnpar(4,neamfnspecloc) = floats(4+nbeg)
      eamfnpar(5,neamfnspecloc) = floats(5+nbeg)
      eamfnpar(6,neamfnspecloc) = floats(6+nbeg)*units
      eamfnpar(7,neamfnspecloc) = floats(7+nbeg)
      eamfnpar(8,neamfnspecloc) = floats(8+nbeg)
      eamfnpar(9,neamfnspecloc) = floats(9+nbeg)
    else
      call outerror('Incorrect coefficient input for MEAM functional',iline)
      call stopnow('potwordm')
    endif
  elseif (neamfn.eq.9) then
!***********
!  Baskes  *
!***********
!
!  Fitting flags
!
    if (lfit.and.lfflags) then
      n1 = int(floats(nfloat-2))
      n2 = int(floats(nfloat-1))
      n3 = int(floats(nfloat))
      nfloat = nfloat - 3
      if (n1.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 4
      endif
      if (n2.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 5
      endif
      if (n3.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 6
      endif
    endif
!
!  If number of floats is greater than 3, assume that fitting flags have been left on line
!
    if (nfloat.gt.3) nfloat = nfloat - 3
!
!  Assign coefficients 
!
    if (nfloat.ge.3) then
      eamfnpar(1,neamfnspecloc) = floats(1+nbeg)
      eamfnpar(2,neamfnspecloc) = floats(2+nbeg)
      eamfnpar(3,neamfnspecloc) = floats(3+nbeg)
    else
      call outerror('Incorrect coefficient input for MEAM functional',iline)
      call stopnow('potwordm')
    endif
    if (abs(eamfnpar(3,neamfnspecloc)).lt.1.0d-12) then
      call outerror('equilibrium density too close to zero in MEAM functional',iline)
      call stopnow('potwordm')
    endif
  elseif (neamfn.eq.10) then
!********
!  VBO  *
!********
!
!  Fitting flags
!
    if (lfit.and.lfflags) then
      n1 = int(floats(nfloat-1))
      n2 = int(floats(nfloat))
      nfloat = nfloat - 2
      if (n1.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 4
      endif
      if (n2.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 5
      endif
    endif
!
!  If number of floats is greater than 2, assume that fitting flags have been left on line
!
    if (nfloat.gt.2) nfloat = nfloat - 2
!
!  Assign coefficients 
!
    if (nfloat.ge.2) then
      eamfnpar(1,neamfnspecloc) = floats(1+nbeg)
      eamfnpar(2,neamfnspecloc) = abs(floats(2+nbeg))
    else
      call outerror('Incorrect coefficient input for MEAM functional',iline)
      call stopnow('potwordm')
    endif
    if (abs(eamfnpar(2,neamfnspecloc)).gt.1.0_dp) then
      call outerror('exponent in MEAM functional should be less than 1',iline)
      call stopnow('potwordm')
    endif
  elseif (neamfn.eq.13) then
!***********
!  Spline  *
!***********
!
!  Fitting flags
!
    if (lfit.and.lfflags) then
      n1 = int(floats(nfloat-3))
      n2 = int(floats(nfloat-2))
      n3 = int(floats(nfloat-1))
      n4 = int(floats(nfloat))
      nfloat = nfloat - 4
      if (n1.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 4
      endif
      if (n2.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 5
      endif
      if (n3.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 6
      endif
      if (n4.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 9
      endif
    endif
!
!  If number of floats is greater than 6, assume that fitting flags have been left on line
!
    if (nfloat.gt.6) nfloat = nfloat - 4
!
!  Assign coefficients 
!
    if (nfloat.ge.6) then
      eamfnpar(1,neamfnspecloc) = floats(1+nbeg)
      eamfnpar(2,neamfnspecloc) = floats(2+nbeg)
      eamfnpar(3,neamfnspecloc) = floats(3+nbeg)
      eamfnpar(4,neamfnspecloc) = floats(4+nbeg)
      eamfnpar(5,neamfnspecloc) = abs(floats(5+nbeg))
      eamfnpar(6,neamfnspecloc) = abs(floats(6+nbeg))
    else
      call outerror('Incorrect coefficient input for MEAM functional',iline)
      call stopnow('potwordm')
    endif
    if (eamfnpar(5,neamfnspecloc).gt.eamfnpar(6,neamfnspecloc)) then
      call outerror('rho0 is greater than rho_max in MEAM functional',iline)
      call stopnow('potwordm')
    endif
  else
!****************************
!  Default MEAM functional  *
!****************************
!
!  Fitting flags
!
    if (lfit.and.lfflags) then
      n1 = int(floats(nfloat))
      nfloat = nfloat - 1
      if (n1.eq.1) then
        nfit = nfit + 1
        if (nfit.ge.maxfit) then
          maxfit = nfit + 10
          call changemaxfit
        endif
        nftyp(nfit) = 5
        nfpot(nfit) = neamfnspecloc
        nfvar(nfit) = 4
      endif
    endif
!
!  If number of floats is greater than 1, assume that fitting flags have been left on line
!
    if (nfloat.gt.1) nfloat = nfloat - 1
!
!  Assign coefficients 
!
    if (nfloat.ge.1) then
      eamfnpar(1,neamfnspecloc) = floats(1+nbeg)
    else
      eamfnpar(1,neamfnspecloc) = 1.0_dp
    endif
  endif
!
!  Read coefficient line
!
  line = ' '
  read(iin,'(a)',end=148) line
  iline = iline + 1
  call linepro(iin,line,iline)
!
!  Fitting flags for coefficients
!
  if (lfit.and.lfflags) then
    if (nfloat.gt.meamorder) then
      do order = 1,meamorder
        n1 = int(floats(nfloat+order-meamorder))
        if (n1.eq.1) then
          nfit = nfit + 1
          if (nfit.ge.maxfit) then
            maxfit = nfit + 10
            call changemaxfit
          endif
          nftyp(nfit) = 5
          nfpot(nfit) = neamfnspecloc
          nfvar(nfit) = 19
          nfvar2(nfit) = order
        endif
      enddo
      nfloat = nfloat - meamorder
    else
      call outerror('Incorrect coefficient input for MEAM functional',iline)
      call stopnow('potwordm')
    endif
  endif
!
!  MEAM functional coefficients
!
  if (nfloat.ge.meamorder) then
    do order = 1,meamorder
      eamfnmeamcoeff(order,neamfnspecloc) = floats(order)
    enddo
  else
    call outerror('Incorrect coefficient input for MEAM functional',iline)
    call stopnow('potwordm')
  endif
!
  goto 145
!
148 if (.not.l55) l1000 = .true.
  return
149 call outerror('file for numerical MEAM cannot be opened',iline)
  call stopnow('potwordm')
!**********************************
!  Parameters for MEAM screening  *
!**********************************
150 units = 1.0_dp
!
!  Read in species parameters for functional
!
155 line = '  '
  read(iin,'(a)',end=158) line
  iline = iline + 1
  call linepro(iin,line,iline)
!
!  Check whether input is symbol or option
!
  if ((nword+nfloat).eq.0) goto 155
  if (nword.gt.0) then
    word = words(1)
    call worsy(word,lsymbol,.true.)
    if (.not.lsymbol) then
      l55 = .true.
      goto 158
    endif
  endif
!
!  Symbol used in input
!
  ltype01 = .false.
  sym1 = ' '
  lself = .false.
  if (nword.gt.0) then
    if (nword.lt.3) then
!
!  Only single species specified => Self term
!
      lself = .true.
      if (llibrary) then
        call okspec(lok1,words(1),ilp1,.true.,ltype01)
        if (.not.lok1) then
          goto 155
        endif
        if (ilp1.gt.0.and.ilp1.le.nspec) then
          nvar1 = natspec(ilp1)
          if (ltype01) then
            itype1 = 0
          else
            itype1 = ntypspec(ilp1)
          endif
        elseif (ilp1.eq.-1) then
          nvar1 = maxele
          itype1 = 0
        endif
      else
        call ltont(words(1),nvar1,itype1)
      endif
      sym1 = words(1)(1:5)
      if (nword.ge.2) then
        if (index(words(2),'she').eq.1) then
          nvar1 = nvar1 + maxele
        endif
      endif
    elseif (nword.eq.6) then
!
!  All three species specified with core/shell attributes
!
      lself = .false.
      if (llibrary) then
        call okspec(lok1,words(1),ilp1,.true.,ltype01)
        if (.not.lok1) then
          goto 155
        endif
        if (ilp1.gt.0.and.ilp1.le.nspec) then
          nvar1 = natspec(ilp1)
          if (ltype01) then
            itype1 = 0
          else
            itype1 = ntypspec(ilp1)
          endif
        elseif (ilp1.eq.-1) then
          nvar1 = maxele
          itype1 = 0
        endif
!
        call okspec(lok2,words(3),ilp2,.true.,ltype02)
        if (.not.lok2) then
          goto 155
        endif
        if (ilp2.gt.0.and.ilp2.le.nspec) then
          nvar2 = natspec(ilp2)
          if (ltype02) then
            itype2 = 0
          else
            itype2 = ntypspec(ilp2)
          endif
        elseif (ilp2.eq.-1) then
          nvar2 = maxele
          itype2 = 0
        endif
!
        call okspec(lok3,words(5),ilp3,.true.,ltype03)
        if (.not.lok3) then
          goto 155
        endif
        if (ilp3.gt.0.and.ilp3.le.nspec) then
          nvar3 = natspec(ilp3)
          if (ltype03) then
            itype3 = 0
          else
            itype3 = ntypspec(ilp3)
          endif
        elseif (ilp3.eq.-1) then
          nvar3 = maxele
          itype3 = 0
        endif
      else
        call ltont(words(1),nvar1,itype1)
        call ltont(words(3),nvar2,itype2)
        call ltont(words(5),nvar3,itype3)
      endif
      sym1 = words(1)(1:5)
      sym2 = words(3)(1:5)
      sym3 = words(5)(1:5)
      if (index(words(2),'she').eq.1) then
        nvar1 = nvar1 + maxele
      endif
      if (index(words(4),'she').eq.1) then
        nvar2 = nvar2 + maxele
      endif
      if (index(words(6),'she').eq.1) then
        nvar3 = nvar3 + maxele
      endif
    elseif (nword.eq.3) then
!
!  All three species specified without core/shell attributes
!
      lself = .false.
      if (llibrary) then
        call okspec(lok1,words(1),ilp1,.true.,ltype01)
        if (.not.lok1) then
          goto 155
        endif
        if (ilp1.gt.0.and.ilp1.le.nspec) then
          nvar1 = natspec(ilp1)
          if (ltype01) then
            itype1 = 0
          else
            itype1 = ntypspec(ilp1)
          endif
        elseif (ilp1.eq.-1) then
          nvar1 = maxele
          itype1 = 0
        endif
!
        call okspec(lok2,words(2),ilp2,.true.,ltype02)
        if (.not.lok2) then
          goto 155
        endif
        if (ilp2.gt.0.and.ilp2.le.nspec) then
          nvar2 = natspec(ilp2)
          if (ltype02) then
            itype2 = 0
          else
            itype2 = ntypspec(ilp2)
          endif
        elseif (ilp2.eq.-1) then
          nvar2 = maxele
          itype2 = 0
        endif
        call okspec(lok3,words(3),ilp3,.true.,ltype03)
        if (.not.lok3) then
          goto 155
        endif
        if (ilp3.gt.0.and.ilp3.le.nspec) then
          nvar3 = natspec(ilp3)
          if (ltype03) then
            itype3 = 0
          else
            itype3 = ntypspec(ilp3)
          endif
        elseif (ilp3.eq.-1) then
          nvar3 = maxele
          itype3 = 0
        endif
      else
        call ltont(words(1),nvar1,itype1)
        call ltont(words(2),nvar2,itype2)
        call ltont(words(3),nvar3,itype3)
      endif
      sym1 = words(1)(1:5)
      sym2 = words(2)(1:5)
      sym3 = words(3)(1:5)
    endif
    nbeg = 0
  else
!
!  Numeric input
!
    nvar1 = int(floats(1))
    if (nvar1.gt.100) nvar1 = nvar1 - 100 + maxele
    itype1 = 0
    nfloat = nfloat - 1
    call label(nvar1,itype1,sym1)
!
    nvar2 = int(floats(2))
    if (nvar2.gt.100) nvar2 = nvar2 - 100 + maxele
    itype2 = 0
    nfloat = nfloat - 1
    call label(nvar2,itype2,sym2)
!
    nvar3 = int(floats(3))
    if (nvar3.gt.100) nvar3 = nvar3 - 100 + maxele
    itype3 = 0
    nfloat = nfloat - 1
    call label(nvar3,itype3,sym3)
!
    nbeg = 3
  endif
!
!  Check to see whether the atomic species are already in the list and if not add them
!
  i = 0
  lfound = .false.       
  do while (i.lt.neamfnspec.and..not.lfound)
    i = i + 1
    lfound = (neamfnnat(i).eq.nvar1.and.neamfntyp(i).eq.itype1)
  enddo    
  if (lfound) then
    neamfnspecloc = i
  else
    neamfnspec = neamfnspec + 1
    neamfnspecloc = neamfnspec
    if (neamfnspec.gt.maxeamfnspec) then
      maxeamfnspec = neamfnspec + 20
      call changemaxeamfnspec
    endif     
    neamfnnat(neamfnspec) = nvar1
    neamfntyp(neamfnspec) = itype1
    symboleamfnspec(neamfnspec) = sym1
  endif
!
  if (lself) then
!
!  Copy species information from self
!
    neamfnspecloc2 = neamfnspecloc
    neamfnspecloc3 = neamfnspecloc
  else
!
!  Check other species
!
    i = 0
    lfound = .false.
    do while (i.lt.neamfnspec.and..not.lfound)
      i = i + 1
      lfound = (neamfnnat(i).eq.nvar2.and.neamfntyp(i).eq.itype2)
    enddo
    if (lfound) then
      neamfnspecloc2 = i
    else
      neamfnspec = neamfnspec + 1
      neamfnspecloc2 = neamfnspec
      if (neamfnspec.gt.maxeamfnspec) then
        maxeamfnspec = neamfnspec + 20
        call changemaxeamfnspec
      endif
      neamfnnat(neamfnspec) = nvar2
      neamfntyp(neamfnspec) = itype2
      symboleamfnspec(neamfnspec) = sym2
    endif
!
    i = 0
    lfound = .false.
    do while (i.lt.neamfnspec.and..not.lfound)
      i = i + 1
      lfound = (neamfnnat(i).eq.nvar3.and.neamfntyp(i).eq.itype3)
    enddo
    if (lfound) then
      neamfnspecloc3 = i
    else
      neamfnspec = neamfnspec + 1
      neamfnspecloc3 = neamfnspec
      if (neamfnspec.gt.maxeamfnspec) then
        maxeamfnspec = neamfnspec + 20
        call changemaxeamfnspec
      endif
      neamfnnat(neamfnspec) = nvar3
      neamfntyp(neamfnspec) = itype3
      symboleamfnspec(neamfnspec) = sym3
    endif
  endif
!
!  Work out index for 2-3 location of parameters
!
  if (neamfnspecloc3.gt.neamfnspecloc2) then
    ind = neamfnspecloc3*(neamfnspecloc3 - 1)/2 + neamfnspecloc2
  else
    ind = neamfnspecloc2*(neamfnspecloc2 - 1)/2 + neamfnspecloc3
  endif
!
  if (nfloat.ge.2) then
    lMEAMscreen(ind,neamfnspecloc) = .true.
    meam_Cmin(ind,neamfnspecloc) = abs(floats(1))
    meam_Cmax(ind,neamfnspecloc) = abs(floats(2))
    if (meam_Cmin(ind,neamfnspecloc).gt.meam_Cmax(ind,neamfnspecloc)) then
      call outerror('Cmin is greater than Cmax for MEAM screening',iline)
      call stopnow('potwordm')
    elseif ((meam_Cmax(ind,neamfnspecloc)-meam_Cmin(ind,neamfnspecloc)).lt.1.0d-12) then
      call outerror('Cmin is too close to Cmax for MEAM screening',iline)
      call stopnow('potwordm')
    endif
!
!  Turn off analytic second derivatives
!
    lanyMEAMscreen = .true.
    lnoanald2 = .true.
  else
    call outerror('Missing parameters for MEAM screening',iline)
    call stopnow('potwordm')
  endif
  goto 155
158 if (.not.l55) l1000 = .true.
  lwordok = .true.
  return
!*************************
!  Type of MEAM density  *
!*************************
160 units = 1.0_dp
  if (nword.ge.2) then
    do i = 2,nword
      if (index(words(i),'t21').eq.1) then
        nmeamrhotype = 1
      elseif (index(words(i),'t24').eq.1) then
        nmeamrhotype = 2
      elseif (index(words(i),'sum').eq.1) then
        nmeamcombotype = 1
      elseif (index(words(i),'exp').eq.1) then
        nmeamcombotype = 2
      else
        call outerror('Unknown MEAM rho type sub-option',iline)
        call stopnow('potwordm')
      endif
    enddo
  endif
  lwordok = .true.
  return
!********************************************
!  Manybody cross term cutoff scale factor  *
!********************************************
170 continue
  if (nfloat.gt.0) then
    eamXcutfactor = abs(floats(1))
  else
    line = '  '
    read(iin,'(a)') line
    iline = iline + 1
    call linepro(iin,line,iline)
    eamXcutfactor = abs(floats(1))
  endif
  if (eamXcutfactor.gt.1.0_dp) then
    call outwarning('Manybody cross term scale factor greater than one',iline)
  endif
  lwordok = .true.
  return
!
  end
