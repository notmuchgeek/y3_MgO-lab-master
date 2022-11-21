  subroutine outfour
!
!  Outputs fourbody interatomic potentials in new format
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
!  Called from outpot
!
!   3/15 Created from outthree
!   4/16 Handling of case where fitlabel is set to Unknown type added
!   4/16 Error in addressing of botyword2 corrected
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
  use four
  use iochannels
  use molecule
  use one
  use parallel
  use plane
  use potentialnames, only : name4pot
  use species
  implicit none
!
!  Local variables
!
  integer(i4),      dimension(:), allocatable  :: iptyp
  character(len=1)                             :: sgs(2)
  character(len=3)                             :: botyword2(3)
  character(len=4)                             :: cstyp1
  character(len=4)                             :: cstyp2
  character(len=4)                             :: cstyp3
  character(len=4)                             :: cstyp4
  character(len=5)                             :: lab1
  character(len=5)                             :: lab2
  character(len=5)                             :: lab3
  character(len=5)                             :: lab4
  character(len=9)                             :: botyword(11)
  integer(i4)                                  :: i
  integer(i4)                                  :: il
  integer(i4)                                  :: ind
  integer(i4)                                  :: j
  integer(i4)                                  :: k
  integer(i4)                                  :: nat1
  integer(i4)                                  :: nat2
  integer(i4)                                  :: nat3
  integer(i4)                                  :: nat4
  integer(i4)                                  :: nbotyptr
  integer(i4)                                  :: nbotyptr2
  integer(i4)                                  :: nmpt(3)
  integer(i4)                                  :: np
  integer(i4)                                  :: ntype1
  integer(i4)                                  :: ntype2
  integer(i4)                                  :: ntype3
  integer(i4)                                  :: ntype4
  integer(i4)                                  :: status
  logical                                      :: lparamline
  real(dp)                                     :: value
!
  data (sgs(j),j=1,2)/'+','-'/
  data botyword/'         ','single   ','double   ','triple   ', &
                'quadruple','resonant ','amide    ','custom   ', &
                'half     ','quarter  ','third    '/
  data botyword2/'   ','cyc','exo'/
!
  if (nfor.gt.0) then
!
!  Count number of general, intra- and inter-molecular potentials
!
    nmpt(1) = 0
    nmpt(2) = 0
    nmpt(3) = 0
    allocate(iptyp(nfor),stat=status)
    if (status/=0) call outofmemory('outfour','iptyp')
    if (lmol) then
      nmpt(1) = 0
      do i = 1,nfor
        if (lfintra(i).and..not.lfinter(i)) then
          nmpt(2) = nmpt(2) + 1
          iptyp(i) = 2
        elseif (lfinter(i).and..not.lfintra(i)) then
          nmpt(3) = nmpt(3) + 1
          iptyp(i) = 3
        else
          nmpt(1) = nmpt(1) + 1
          iptyp(i) = 1
        endif
      enddo
    else
      nmpt(1) = nfor
      do i = 1,nfor
        iptyp(i) = 1
      enddo
    endif
!********************************
!  Output four-body potentials  *
!********************************
    do k = 1,3
      if (nmpt(k).gt.0) then
        if (k.eq.1) then
          write(ioout,'(/,''  General Four-body potentials :'',/)')
        elseif (k.eq.2) then
          write(ioout,'(/,''  Intramolecular Four-body potentials :'',/)')
        elseif (k.eq.3) then
          write(ioout,'(/,''  Intermolecular Four-body potentials :'',/)')
        endif
        write(ioout,'(''--------------------------------------------------------------------------------'')')
        write(ioout,'(''  Atoms       Potential      Parameter       Value         Units   Cutoffs(Ang)'')')
        write(ioout,'(''  1 - 4                                                             Min /  Max '')')
        write(ioout,'(''--------------------------------------------------------------------------------'')')
        do i = 1,nfor
          if (iptyp(i).eq.k) then
            np = nforty(i)
            nat1 = nfspec1(i)
            nat2 = nfspec2(i)
            nat3 = nfspec3(i)
            nat4 = nfspec4(i)
            ntype1 = nfptyp1(i)
            ntype2 = nfptyp2(i)
            ntype3 = nfptyp3(i)
            ntype4 = nfptyp4(i)
            call label(nat1,ntype1,lab1)
            call label(nat2,ntype2,lab2)
            call label(nat3,ntype3,lab3)
            call label(nat4,ntype4,lab4)
            cstyp1 = 'core'
            cstyp2 = 'core'
            cstyp3 = 'core'
            cstyp4 = 'core'
            if (nat1.gt.maxele) cstyp1 = 'shel'
            if (nat2.gt.maxele) cstyp2 = 'shel'
            if (nat3.gt.maxele) cstyp3 = 'shel'
            if (nat4.gt.maxele) cstyp4 = 'shel'
            if (npfor(i).gt.0) then
              ind = 1
            else
              ind = 2
            endif
            nbotyptr  = n4botype(1,i) + 1
            nbotyptr2 = n4botype(2,i)
!
            lparamline = .true.
            if (nfitlabel(np,4).gt.0) then
              do il = 1,max(4_i4,nfitlabel(np,4))
                lparamline = (il.le.nfitlabel(np,4))
                if (fitlabel(il,np,4)(1:1).eq.' '.or.fitlabel(il,np,4)(1:7).eq.'Unknown') lparamline = .false.
                if (lparamline) then
                  if (il.eq.1) then
                    value = fork(i)
                  elseif (il.ge.2) then
                    value = forpoly(il-1,i)
                  endif
                endif
!
                if (lparamline) then
                  if (il.eq.1) then
!
!  First line 
!
                    if (mmfexc(i).eq.0) then
                      write(ioout,'(a5,1x,a4,4x,a13,1x,a15,g15.8,a10,f5.3,1x,f6.3)') &
                        lab1,cstyp1,name4pot(np),fitlabel(il,np,4),value,fitlabelunits(il,np,4),for1min(i),for1(i)
                    elseif (mmfexc(i).eq.1) then
                      write(ioout,'(a5,1x,a4,4x,a13,1x,a15,g15.8,a10,1x,''Bonded'')') &
                        lab1,cstyp1,name4pot(np),fitlabel(il,np,4),value,fitlabelunits(il,np,4)
                    elseif (mmfexc(i).eq.2) then
                      write(ioout,'(a5,1x,a4,4x,a13,1x,a15,g15.8,a10,1x,''Improper'')') &
                        lab1,cstyp1,name4pot(np),fitlabel(il,np,4),value,fitlabelunits(il,np,4)
                    endif
                  elseif (il.eq.2) then
!
!  Second line
!
                    if (npfor(i).ne.0) then
!
!  Potential has n value to output
!
                      if (mmfexc(i).eq.0) then
                        write(ioout,'(a5,1x,a4,'' Phase = '',a1,1x,i2,5x,a15,g15.8,a10,f5.3,1x,f6.3)') &
                          lab2,cstyp2,sgs(ind),abs(npfor(i)),fitlabel(il,np,4),value,fitlabelunits(il,np,4),for2min(i),for2(i)
                      elseif (mmfexc(i).eq.1) then
                        write(ioout,'(a5,1x,a4,'' Phase = '',a1,1x,i2,5x,a15,g15.8,a10,1x,a9)') &
                          lab2,cstyp2,sgs(ind),abs(npfor(i)),fitlabel(il,np,4),value,fitlabelunits(il,np,4),botyword(nbotyptr)
                      elseif (mmfexc(i).eq.2) then
                        write(ioout,'(a5,1x,a4,'' Phase = '',a1,1x,i2,5x,a15,g15.8,a10,1x,a9)') &
                          lab2,cstyp2,sgs(ind),abs(npfor(i)),fitlabel(il,np,4),value,fitlabelunits(il,np,4),botyword(nbotyptr)
                      endif
                    else
!
!  No n value
!
                      if (mmfexc(i).eq.0) then
                        write(ioout,'(a5,1x,a4,18x,a15,g15.8,a10,f5.3,1x,f6.3)') &
                          lab2,cstyp2,fitlabel(il,np,4),value,fitlabelunits(il,np,4),for2min(i),for2(i)
                      elseif (mmfexc(i).eq.1) then
                        write(ioout,'(a5,1x,a4,18x,a15,g15.8,a10,1x,a9)') &
                          lab2,cstyp2,fitlabel(il,np,4),value,fitlabelunits(il,np,4),botyword(nbotyptr)
                      elseif (mmfexc(i).eq.2) then
                        write(ioout,'(a5,1x,a4,18x,a15,g15.8,a10,1x,a9)') &
                          lab2,cstyp2,fitlabel(il,np,4),value,fitlabelunits(il,np,4),botyword(nbotyptr)
                      endif
                    endif
                  elseif (il.eq.3) then
!
!  Third line
!
                    if (mmfexc(i).eq.0) then
                      write(ioout,'(a5,1x,a4,18x,a15,g15.8,a10,f5.3,1x,f6.3)') &
                        lab3,cstyp3,fitlabel(il,np,4),value,fitlabelunits(il,np,4),for3min(i),for3(i)
                    elseif (mmfexc(i).eq.1) then
                      write(ioout,'(a5,1x,a4,18x,a15,g15.8,a10,1x,a3)') &
                        lab3,cstyp3,fitlabel(il,np,4),value,fitlabelunits(il,np,4),botyword2(nbotyptr2)
                    elseif (mmfexc(i).eq.2) then
                      write(ioout,'(a5,1x,a4,18x,a15,g15.8,a10,1x,a3)') &
                        lab3,cstyp3,fitlabel(il,np,4),value,fitlabelunits(il,np,4),botyword2(nbotyptr2)
                    endif
                  elseif (il.eq.4) then
!
!  Fourth line
!
                    if (mmfexc(i).eq.0) then
                      write(ioout,'(a5,1x,a4,18x,a15,g15.8,a10,f5.3,1x,f6.3)') &
                        lab4,cstyp4,fitlabel(il,np,4),value,fitlabelunits(il,np,4),for4min(i),for4(i)
                    else
                      write(ioout,'(a5,1x,a4,18x,a15,g15.8,a10)') &
                        lab4,cstyp4,fitlabel(il,np,4),value,fitlabelunits(il,np,4)
                    endif
                  else
!
!  Subsequent lines have just potential term
!
                    write(ioout,'(28x,a15,g15.8,a10)') &
                      fitlabel(il,np,4),value,fitlabelunits(il,np,4)
                  endif
                else
!
!  Parameter free lines
!
                  if (il.eq.1) then
!
!  First line 
!
                    if (mmfexc(i).eq.0) then
                      write(ioout,'(a5,1x,a4,4x,a13,41x,f5.3,1x,f6.3)') &
                        lab1,cstyp1,name4pot(np),for1min(i),for1(i)
                    elseif (mmfexc(i).eq.1) then
                      write(ioout,'(a5,1x,a4,4x,a13,42x,''Bonded'')') &
                        lab1,cstyp1,name4pot(np)
                    elseif (mmfexc(i).eq.2) then
                      write(ioout,'(a5,1x,a4,4x,a13,42x,''Improper'')') &
                        lab1,cstyp1,name4pot(np)
                    endif
                  elseif (il.eq.2) then
!
!  Second line
!
                    if (npfor(i).ne.0) then
!
!  Potential has n value to output
!
                      if (mmfexc(i).eq.0) then
                        write(ioout,'(a5,1x,a4,'' Phase = '',a1,1x,i2,45x,f5.3,1x,f6.3)') &
                          lab2,cstyp2,sgs(ind),abs(npfor(i)),for2min(i),for2(i)
                      elseif (mmfexc(i).eq.1) then
                        write(ioout,'(a5,1x,a4,'' Phase = '',a1,1x,i2,46x,a9)') &
                          lab2,cstyp2,sgs(ind),abs(npfor(i)),botyword(nbotyptr)
                      elseif (mmfexc(i).eq.2) then
                        write(ioout,'(a5,1x,a4,'' Phase = '',a1,1x,i2,46x,a9)') &
                          lab2,cstyp2,sgs(ind),abs(npfor(i)),botyword(nbotyptr)
                      endif
                    else
!
!  No n value
!
                      if (mmfexc(i).eq.0) then
                        write(ioout,'(a5,1x,a4,58x,f5.3,1x,f6.3)') &
                          lab2,cstyp2,for2min(i),for2(i)
                      elseif (mmfexc(i).eq.1) then
                        write(ioout,'(a5,1x,a4,59x,a9)') &
                          lab2,cstyp2,botyword(nbotyptr)
                      elseif (mmfexc(i).eq.2) then
                        write(ioout,'(a5,1x,a4,59x,a9)') &
                          lab2,cstyp2,botyword(nbotyptr)
                      endif
                    endif
                  elseif (il.eq.3) then
!
!  Third line
!
                    if (mmfexc(i).eq.0) then
                      write(ioout,'(a5,1x,a4,58x,f5.3,1x,f6.3)') &
                        lab3,cstyp3,for3min(i),for3(i)
                    elseif (mmfexc(i).eq.1) then
                      write(ioout,'(a5,1x,a4,59x,a3)') & 
                        lab3,cstyp3,botyword2(nbotyptr2)
                    elseif (mmfexc(i).eq.2) then
                      write(ioout,'(a5,1x,a4,59x,a3)') & 
                        lab3,cstyp3,botyword2(nbotyptr2)
                    endif
                  elseif (il.eq.4) then
!
!  Fourth line
!
                    if (mmfexc(i).eq.0) then
                      write(ioout,'(a5,1x,a4,58x,f5.3,1x,f6.3)') &
                        lab4,cstyp4,for4min(i),for4(i)
                    elseif (mmfexc(i).eq.1) then
                      write(ioout,'(a5,1x,a4)') & 
                        lab4,cstyp4
                    elseif (mmfexc(i).eq.2) then
                      write(ioout,'(a5,1x,a4)') & 
                        lab4,cstyp4
                    endif
                  endif
                endif
              enddo
            endif
            write(ioout,'(''--------------------------------------------------------------------------------'')')
          endif
        enddo
      endif
    enddo
    deallocate(iptyp,stat=status)
    if (status/=0) call deallocate_error('outfour','iptyp')
  endif
!
  return
  end
