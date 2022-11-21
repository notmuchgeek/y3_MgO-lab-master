  subroutine outtwo
!
!  Outputs twobody interatomic potentials in new format
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
!  59 = Lennard-Jones with epsilon/sigma input : c1 = c2 = 4
!  60 = Buffered Lennard-Jones
!  61 = Slater
!  62 = Constant J2 potential
!  63 = Stillinger-Weber 2 body potential - general form
!  64 = Screened Coulomb potential 
!
!  Called from outpot
!
!   3/15 Created from outpot
!   3/15 x14 sub-option added
!   3/15 namepot renamed to name2pot
!   4/16 Handling of case where fitlabel is set to Unknown type added
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
!  Copyright Curtin University 2016
!
!  Julian Gale, CIC, Curtin University, April 2016
!
  use configurations
  use g_constants
  use control
  use current
  use element,        only : maxele
  use fitting,        only : nfitlabel, fitlabel, fitlabelunits
  use general,        only : nwarn
  use iochannels
  use molecule
  use one
  use parallel
  use plane
  use potentialnames, only : name2pot
  use shells
  use shifts
  use species
  use two
  implicit none
!
!  Local variables
!
  integer(i4),      dimension(:), allocatable  :: iptyp
  integer(i4),      dimension(:), allocatable  :: ntemp
  character(len=1)                             :: cstyp1
  character(len=1)                             :: cstyp2
  character(len=3)                             :: botyword2(3)
  character(len=5)                             :: lab1
  character(len=5)                             :: lab2
  character(len=9)                             :: botyword(11)
  integer(i4)                                  :: i
  integer(i4)                                  :: il
  integer(i4)                                  :: k
  integer(i4)                                  :: mpt
  integer(i4)                                  :: nat1
  integer(i4)                                  :: nat2
  integer(i4)                                  :: nmpt(3)
  integer(i4)                                  :: np
  integer(i4)                                  :: npi
  integer(i4)                                  :: npt
  integer(i4)                                  :: ntmp
  integer(i4)                                  :: ntype1
  integer(i4)                                  :: ntype2
  integer(i4)                                  :: status
  logical                                      :: lfirstline
  logical                                      :: lwarnpn
  real(dp)                                     :: value
!
  data botyword/'         ','single   ','double   ','triple   ', &
                'quadruple','resonant ','amide    ','custom   ', &
                'half     ','quarter  ','third    '/
  data botyword2/'   ','cyc','exo'/
!
  if (npote.gt.0) then
!
!  Sort potentials for those which are inter- and intra- molecular
!
    nmpt(2) = 0
    nmpt(3) = 0
    allocate(iptyp(npote),stat=status)
    if (status/=0) call outofmemory('outtwo','iptyp')
    if (lmol) then
      nmpt(1) = 0
      do i = 1,npote
        if (lintra(i).and..not.linter(i)) then
          nmpt(2) = nmpt(2) + 1
          iptyp(i) = 2
        elseif (linter(i).and..not.lintra(i)) then
          nmpt(3) = nmpt(3) + 1
          iptyp(i) = 3
        else
          nmpt(1) = nmpt(1) + 1
          iptyp(i) = 1
        endif
      enddo
    else
      nmpt(1) = npote
      do i = 1,npote
        iptyp(i) = 1
      enddo
    endif
!**********************************
!  Output interatomic potentials  *
!**********************************
    lwarnpn = .false.
    do k = 1,3
      if (nmpt(k).gt.0) then
        if (k.eq.1) then
          write(ioout,'(/,''  General interatomic potentials :'',/)')
        elseif (k.eq.2) then
          write(ioout,'(/,''  Intramolecular potentials :'',/)')
        else
          write(ioout,'(/,''  Intermolecular potentials :'',/)')
        endif
        write(ioout,'(''--------------------------------------------------------------------------------'')')
        write(ioout,'(''Atom  Types   Potential      Parameter       Value         Units   Cutoffs(Ang)'')')
        write(ioout,'(''  1     2                                                            Min /  Max '')')
        write(ioout,'(''--------------------------------------------------------------------------------'')')
        do i = 1,npote
          if (iptyp(i).eq.k) then
            np = nptype(i)
            nat1 = nspec1(i)
            nat2 = nspec2(i)
            ntype1 = nptyp1(i)
            ntype2 = nptyp2(i)
            call label(nat1,ntype1,lab1)
            call label(nat2,ntype2,lab2)
            cstyp1 = 'c'
            cstyp2 = 'c'
            if (nat1.gt.maxele) cstyp1 = 's'
            if (nat2.gt.maxele) cstyp2 = 's'
!
            lfirstline = .true.
            if (nfitlabel(np,2).gt.0) then
              do il = 1,nfitlabel(np,2)
                if (il.le.9) then
                  value = twopot(il,i)
                elseif (il.ge.10) then
                  value = tpot(il-9,i)
                endif
!
!  Check that this isn't a blank unused term for this potential or unknown type
!
                if (fitlabel(il,np,2)(1:1).ne.' '.and.fitlabel(il,np,2)(1:7).ne.'Unknown') then
                  if (lfirstline) then
!
!  First line has labels and cut-offs
!
                    lfirstline = .false.
                    if (mmexc(i).eq.0) then
                      write(ioout,'(2(a5,a1,1x),a13,1x,a15,g15.8,a10,f5.3,1x,f6.3)') &
                        lab1,cstyp1,lab2,cstyp2,name2pot(np),fitlabel(il,np,2),value,fitlabelunits(il,np,2),rpot2(i),rpot(i)
                    elseif (mmexc(i).eq.1) then
                      write(ioout,'(2(a5,a1,1x),a13,1x,a15,g15.8,a10,f5.3,1x,''1 Bond'')') &
                        lab1,cstyp1,lab2,cstyp2,name2pot(np),fitlabel(il,np,2),value,fitlabelunits(il,np,2),rpot2(i)
                      if (k.eq.3) lwarnpn = .true.
                    elseif (mmexc(i).eq.2) then
                      write(ioout,'(2(a5,a1,1x),a13,1x,a15,g15.8,a10,''2Bond'',1x,f6.3)') &
                        lab1,cstyp1,lab2,cstyp2,name2pot(np),fitlabel(il,np,2),value,fitlabelunits(il,np,2),rpot(i)
                    elseif (mmexc(i).eq.3) then
                      write(ioout,'(2(a5,a1,1x),a13,1x,a15,g15.8,a10,''3Bond'',1x,f6.3)') &
                        lab1,cstyp1,lab2,cstyp2,name2pot(np),fitlabel(il,np,2),value,fitlabelunits(il,np,2),rpot(i)
                    elseif (mmexc(i).eq.4) then
                      write(ioout,'(2(a5,a1,1x),a13,1x,a15,g15.8,a10,''1-4only'',f5.2)') &
                        lab1,cstyp1,lab2,cstyp2,name2pot(np),fitlabel(il,np,2),value,fitlabelunits(il,np,2),rpot(i)
                    elseif (mmexc(i).eq.5) then
                      write(ioout,'(2(a5,a1,1x),a13,1x,a15,g15.8,a10,''1-4gtr '',f5.2)') &
                        lab1,cstyp1,lab2,cstyp2,name2pot(np),fitlabel(il,np,2),value,fitlabelunits(il,np,2),rpot(i)
                    elseif (mmexc(i).eq.6) then
                      write(ioout,'(2(a5,a1,1x),a13,1x,a15,g15.8,a10,''4Bond'',1x,f6.3)') &
                        lab1,cstyp1,lab2,cstyp2,name2pot(np),fitlabel(il,np,2),value,fitlabelunits(il,np,2),rpot(i)
                    endif
                  else
!
!  Subsequent lines have just potential term
!
                    write(ioout,'(28x,a15,g15.8,a10)') &
                      fitlabel(il,np,2),value,fitlabelunits(il,np,2)
                  endif
                endif
              enddo
              if (np.eq.7.or.np.eq.2.or.np.eq.12.or.np.eq.13.or.np.eq.21) then
                mpt = int(tpot(1,i))
                npt = int(tpot(2,i))
                write(ioout,'(28x,''Exponent m'',5x,i5,10x,''None'')') mpt
                write(ioout,'(28x,''Exponent n'',5x,i5,10x,''None'')') npt
              elseif (np.eq.56) then
                write(ioout,'(28x,''Exponent n'',5x,g15.8,''None'')') twopot(10,i)
              endif
            else
!
!  Special case of potential with no parameters
!
              if (mmexc(i).eq.0) then
                write(ioout,'(2(a5,a1,1x),a13,41x,f5.3,1x,f6.3)') &
                  lab1,cstyp1,lab2,cstyp2,name2pot(np),rpot2(i),rpot(i)
              elseif (mmexc(i).eq.1) then
                write(ioout,'(2(a5,a1,1x),a13,41x,f5.3,1x,''1 Bond'')') &
                  lab1,cstyp1,lab2,cstyp2,name2pot(np),rpot2(i)
                if (k.eq.3) lwarnpn = .true.
              elseif (mmexc(i).eq.2) then
                write(ioout,'(2(a5,a1,1x),a13,41x,''2Bond'',1x,f6.3)') &
                  lab1,cstyp1,lab2,cstyp2,name2pot(np),rpot(i)
              elseif (mmexc(i).eq.3) then
                write(ioout,'(2(a5,a1,1x),a13,41x,''3Bond'',1x,f6.3)') &
                  lab1,cstyp1,lab2,cstyp2,name2pot(np),rpot(i)
              elseif (mmexc(i).eq.4) then
                write(ioout,'(2(a5,a1,1x),a13,41x,''1-4only'',f5.2)') &
                  lab1,cstyp1,lab2,cstyp2,name2pot(np),rpot(i)
              elseif (mmexc(i).eq.5) then
                write(ioout,'(2(a5,a1,1x),a13,41x,''1-4gtr '',f5.2)') &
                  lab1,cstyp1,lab2,cstyp2,name2pot(np),rpot(i)
              elseif (mmexc(i).eq.6) then
                write(ioout,'(2(a5,a1,1x),a13,41x,''4Bond'',1x,f6.3)') &
                  lab1,cstyp1,lab2,cstyp2,name2pot(np),rpot(i)
              endif
            endif
!
!  Bond types
!
            if (n2botype(1,i).ne.0.or.n2botype(2,i).ne.0) then
              write(ioout,'(68x,a9,a3)') botyword(n2botype(1,i)+1),botyword2(n2botype(2,i)+1)
            endif
            write(ioout,'(''--------------------------------------------------------------------------------'')')
          endif
        enddo
      endif
    enddo
    deallocate(iptyp,stat=status)
    if (status/=0) call deallocate_error('outtwo','iptyp')
    if (lwarnpn) then
!
!  Intermolecular potential has been specified as being bonded only
!
      nwarn = nwarn + 1
      call outwarning('Intermolecular potential is specified as bonded potential',0_i4)
    endif
  endif
!**************************************************
!  Output potentials with energy/gradient shifts  *
!**************************************************
  ntmp = 0
  allocate(ntemp(npote),stat=status)
  if (status/=0) call outofmemory('outtwo','ntemp')
  do i = 1,npote
    if (leshift(i).and.(nptype(i).ne.10.and.nptype(i).ne.23).and.nptype(i).lt.100) then
      ntmp = ntmp + 1
      ntemp(ntmp) = i
    endif
  enddo
  if (ntmp.gt.0) then
    write(ioout,'(/,''  Potentials with energy / gradient boundary corrections : '',/)')
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    write(ioout,'(''    Potential       Atom1   Atom2    Type of correction'')')
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    do i = 1,ntmp
      npi = ntemp(i)
      np = nptype(npi)
      nat1 = nspec1(npi)
      nat2 = nspec2(npi)
      ntype1 = nptyp1(npi)
      ntype2 = nptyp2(npi)
      call label(nat1,ntype1,lab1)
      call label(nat2,ntype2,lab2)
      cstyp1 = 'c'
      cstyp2 = 'c'
      if (nat1.gt.maxele) cstyp1 = 's'
      if (nat2.gt.maxele) cstyp2 = 's'
      if (np.eq.7) then
        mpt = int(tpot(1,npi))
        npt = int(tpot(2,npi))
      endif
      if (tpot(5,npi).ne.0.0_dp.and.np.ne.10) then
        if (np.eq.7) then
          write(ioout,'(4x,a7,2i3,3x,a5,1x,a1,1x,a5,1x,a1,2x,''Voter'')') &
            name2pot(7),mpt,npt,lab1,cstyp1,lab2,cstyp2
        else
          write(ioout,'(4x,a13,3x,a5,1x,a1,1x,a5,1x,a1,2x,''Voter'')') &
            name2pot(np),lab1,cstyp1,lab2,cstyp2
        endif
      elseif (leshift(npi).and.lgshift(npi).and.np.ne.10) then
        if (np.eq.7) then
          write(ioout,'(4x,a7,2i3,3x,a5,1x,a1,1x,a5,1x,a1,2x,''Gradient'')') &
            name2pot(7),mpt,npt,lab1,cstyp1,lab2,cstyp2
        else
          write(ioout,'(4x,a13,3x,a5,1x,a1,1x,a5,1x,a1,2x,''Gradient'')') &
            name2pot(np),lab1,cstyp1,lab2,cstyp2
        endif
      else
        if (np.eq.7) then
          write(ioout,'(4x,a7,2i3,3x,a5,1x,a1,1x,a5,1x,a1,2x,''Energy'')') &
            name2pot(7),mpt,npt,lab1,cstyp1,lab2,cstyp2
        else
          write(ioout,'(4x,a13,3x,a5,1x,a1,1x,a5,1x,a1,2x,''Energy'')') &
            name2pot(np),lab1,cstyp1,lab2,cstyp2
        endif
      endif
    enddo
    write(ioout,'(''--------------------------------------------------------------------------------'')')
  endif
  deallocate(ntemp,stat=status)
  if (status/=0) call deallocate_error('outtwo','ntemp')
!
  return
  end
