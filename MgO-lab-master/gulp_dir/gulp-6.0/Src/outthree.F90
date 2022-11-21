  subroutine outthree
!
!  Outputs threebody interatomic potentials in new format
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
!  25 = SW3 with Garofalini form
!  26 = j3
!  27 = ppp3body
!
!  Called from outpot
!
!   3/15 Created from outtwo
!   8/15 Garofalini form of sw3 added
!   4/16 Handling of case where fitlabel is set to Unknown type added
!  11/19 ppp3body added
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
!  Julian Gale, CIC, Curtin University, November 2019
!
  use configurations
  use g_constants
  use control
  use current
  use element,        only : maxele
  use fitting,        only : nfitlabel, fitlabel, fitlabelunits
  use iochannels
  use m_three
  use molecule
  use one
  use parallel
  use plane
  use potentialnames, only : name3pot
  use shells
  use species
  implicit none
!
!  Local variables
!
  integer(i4),      dimension(:), allocatable  :: iptyp
  character(len=4)                             :: cstyp1
  character(len=4)                             :: cstyp2
  character(len=4)                             :: cstyp3
  character(len=5)                             :: lab1
  character(len=5)                             :: lab2
  character(len=5)                             :: lab3
  integer(i4)                                  :: i
  integer(i4)                                  :: il
  integer(i4)                                  :: k
  integer(i4)                                  :: nat1
  integer(i4)                                  :: nat2
  integer(i4)                                  :: nat3
  integer(i4)                                  :: nmpt(3)
  integer(i4)                                  :: np
  integer(i4)                                  :: ntype1
  integer(i4)                                  :: ntype2
  integer(i4)                                  :: ntype3
  integer(i4)                                  :: status
  logical                                      :: lparamline
  real(dp)                                     :: value
!
  if (nthb.gt.0) then
!
!  Count number of general, intra- and inter-molecular potentials
!
    nmpt(1) = 0
    nmpt(2) = 0
    nmpt(3) = 0
    allocate(iptyp(nthb),stat=status)
    if (status/=0) call outofmemory('outthree','iptyp')
    if (lmol) then
      nmpt(1) = 0
      do i = 1,nthb
        if (ltintra(i).and..not.ltinter(i)) then
          nmpt(2) = nmpt(2) + 1
          iptyp(i) = 2
        elseif (ltinter(i).and..not.ltintra(i)) then
          nmpt(3) = nmpt(3) + 1
          iptyp(i) = 3
        else
          nmpt(1) = nmpt(1) + 1
          iptyp(i) = 1
        endif
      enddo
    else
      nmpt(1) = nthb
      do i = 1,nthb
        iptyp(i) = 1
      enddo
    endif
!*********************************
!  Output three-body potentials  *
!*********************************
    do k = 1,3
      if (nmpt(k).gt.0) then
        if (k.eq.1) then
          write(ioout,'(/,''  General Three-body potentials :'',/)')
        elseif (k.eq.2) then
          write(ioout,'(/,''  Intramolecular Three-body potentials :'',/)')
        elseif (k.eq.3) then
          write(ioout,'(/,''  Intermolecular Three-body potentials :'',/)')
        endif
        write(ioout,'(''--------------------------------------------------------------------------------'')')
        write(ioout,'(''  Atoms       Potential      Parameter       Value         Units   Cutoffs(Ang)'')')
        write(ioout,'(''  1 / 2 / 3                                                          Min /  Max '')')
        write(ioout,'(''--------------------------------------------------------------------------------'')')
        do i = 1,nthb
          if (iptyp(i).eq.k) then
            np = nthrty(i)
            nat1 = ntspec1(i)
            nat2 = ntspec2(i)
            nat3 = ntspec3(i)
            ntype1 = ntptyp1(i)
            ntype2 = ntptyp2(i)
            ntype3 = ntptyp3(i)
            call label(nat1,ntype1,lab1)
            call label(nat2,ntype2,lab2)
            call label(nat3,ntype3,lab3)
            cstyp1 = 'core'
            cstyp2 = 'core'
            cstyp3 = 'core'
            if (nat1.gt.maxele) cstyp1 = 'shel'
            if (nat2.gt.maxele) cstyp2 = 'shel'
            if (nat3.gt.maxele) cstyp3 = 'shel'
!
            lparamline = .true.
            if (nfitlabel(np,3).gt.0) then
              do il = 1,max(3_i4,nfitlabel(np,3))
                lparamline = (il.le.nfitlabel(np,3))
                if (fitlabel(il,np,3)(1:1).eq.' '.or.fitlabel(il,np,3)(1:7).eq.'Unknown') lparamline = .false.
                if (lparamline) then
                  if (il.eq.1) then
                    value = thbk(i)
                  elseif (il.eq.2) then
                    value = theta(i)
                  elseif (il.eq.3) then
                    value = thrho1(i)
                  elseif (il.eq.4) then
                    value = thrho2(i)
                  elseif (il.eq.5) then
                    value = thrho3(i)
                  elseif (il.ge.5) then
                    value = threepoly(il-5,i)
                  endif
                endif
!
                if (lparamline) then
                  if (il.eq.1) then
!
!  First line 
!
                    if (mmtexc(i).eq.0) then
                      write(ioout,'(a5,1x,a4,4x,a13,1x,a15,g15.8,a10,f5.3,1x,f6.3)') &
                        lab1,cstyp1,name3pot(np),fitlabel(il,np,3),value,fitlabelunits(il,np,3),thr1min(i),thr1(i)
                    elseif (mmtexc(i).eq.1) then
                      write(ioout,'(a5,1x,a4,4x,a13,1x,a15,g15.8,a10,f5.3,1x,''1 Bond'')') &
                        lab1,cstyp1,name3pot(np),fitlabel(il,np,3),value,fitlabelunits(il,np,3),thr1min(i)
                    elseif (mmtexc(i).eq.2) then
                      write(ioout,'(a5,1x,a4,4x,a13,1x,a15,g15.8,a10,''2Bond'',1x,f6.3)') &
                        lab1,cstyp1,name3pot(np),fitlabel(il,np,3),value,fitlabelunits(il,np,3),thr1(i)
                    elseif (mmtexc(i).eq.3) then
                      write(ioout,'(a5,1x,a4,4x,a13,1x,a15,g15.8,a10,''3Bond'',1x,f6.3)') &
                        lab1,cstyp1,name3pot(np),fitlabel(il,np,3),value,fitlabelunits(il,np,3),thr1(i)
                    elseif (mmtexc(i).eq.4) then
                      write(ioout,'(a5,1x,a4,4x,a13,1x,a15,g15.8,a10,''1-4only'',f5.2)') &
                        lab1,cstyp1,name3pot(np),fitlabel(il,np,3),value,fitlabelunits(il,np,3),thr1(i)
                    elseif (mmtexc(i).eq.5) then
                      write(ioout,'(a5,1x,a4,4x,a13,1x,a15,g15.8,a10,''1-4gtr '',f5.2)') &
                        lab1,cstyp1,name3pot(np),fitlabel(il,np,3),value,fitlabelunits(il,np,3),thr1(i)
                    elseif (mmtexc(i).eq.6) then
                      write(ioout,'(a5,1x,a4,4x,a13,1x,a15,g15.8,a10,''4Bond'',1x,f6.3)') &
                        lab1,cstyp1,name3pot(np),fitlabel(il,np,3),value,fitlabelunits(il,np,3),thr1(i)
                    endif
                  elseif (il.eq.2) then
!
!  Second line
!
                    if (mmtexc(i).eq.0) then
                      write(ioout,'(a5,1x,a4,18x,a15,g15.8,a10,f5.3,1x,f6.3)') &
                        lab2,cstyp2,fitlabel(il,np,3),value,fitlabelunits(il,np,3),thr2min(i),thr2(i)
                    elseif (mmtexc(i).eq.1) then
                      write(ioout,'(a5,1x,a4,18x,a15,g15.8,a10,f5.3,1x,''1 Bond'')') &
                        lab2,cstyp2,fitlabel(il,np,3),value,fitlabelunits(il,np,3),thr2min(i)
                    elseif (mmtexc(i).eq.2) then
                      write(ioout,'(a5,1x,a4,18x,a15,g15.8,a10,''2Bond'',1x,f6.3)') &
                        lab2,cstyp2,fitlabel(il,np,3),value,fitlabelunits(il,np,3),thr2(i)
                    elseif (mmtexc(i).eq.3) then
                      write(ioout,'(a5,1x,a4,18x,a15,g15.8,a10,''3Bond'',1x,f6.3)') &
                        lab2,cstyp2,fitlabel(il,np,3),value,fitlabelunits(il,np,3),thr2(i)
                    elseif (mmtexc(i).eq.4) then
                      write(ioout,'(a5,1x,a4,18x,a15,g15.8,a10,''1-4only'',f5.2)') &
                        lab2,cstyp2,fitlabel(il,np,3),value,fitlabelunits(il,np,3),thr2(i)
                    elseif (mmtexc(i).eq.5) then
                      write(ioout,'(a5,1x,a4,18x,a15,g15.8,a10,''1-4gtr '',f5.2)') &
                        lab2,cstyp2,fitlabel(il,np,3),value,fitlabelunits(il,np,3),thr2(i)
                    elseif (mmtexc(i).eq.6) then
                      write(ioout,'(a5,1x,a4,18x,a15,g15.8,a10,''4Bond'',1x,f6.3)') &
                        lab2,cstyp2,fitlabel(il,np,3),value,fitlabelunits(il,np,3),thr2(i)
                    endif
                  elseif (il.eq.3) then
!
!  Third line
!
                    if (mmtexc(i).eq.0) then
                      write(ioout,'(a5,1x,a4,18x,a15,g15.8,a10,f5.3,1x,f6.3)') &
                        lab3,cstyp3,fitlabel(il,np,3),value,fitlabelunits(il,np,3),thr3min(i),thr3(i)
                    elseif (mmtexc(i).eq.1) then
                      write(ioout,'(a5,1x,a4,18x,a15,g15.8,a10,f5.3,1x,''1 Bond'')') &
                        lab3,cstyp3,fitlabel(il,np,3),value,fitlabelunits(il,np,3),thr3min(i)
                    elseif (mmtexc(i).eq.2) then
                      write(ioout,'(a5,1x,a4,18x,a15,g15.8,a10,''2Bond'',1x,f6.3)') &
                        lab3,cstyp3,fitlabel(il,np,3),value,fitlabelunits(il,np,3),thr3(i)
                    elseif (mmtexc(i).eq.3) then
                      write(ioout,'(a5,1x,a4,18x,a15,g15.8,a10,''3Bond'',1x,f6.3)') &
                        lab3,cstyp3,fitlabel(il,np,3),value,fitlabelunits(il,np,3),thr3(i)
                    elseif (mmtexc(i).eq.4) then
                      write(ioout,'(a5,1x,a4,18x,a15,g15.8,a10,''1-4only'',f5.2)') &
                        lab3,cstyp3,fitlabel(il,np,3),value,fitlabelunits(il,np,3),thr3(i)
                    elseif (mmtexc(i).eq.5) then
                      write(ioout,'(a5,1x,a4,18x,a15,g15.8,a10,''1-4gtr '',f5.2)') &
                        lab3,cstyp3,fitlabel(il,np,3),value,fitlabelunits(il,np,3),thr3(i)
                    elseif (mmtexc(i).eq.6) then
                      write(ioout,'(a5,1x,a4,18x,a15,g15.8,a10,''4Bond'',1x,f6.3)') &
                        lab3,cstyp3,fitlabel(il,np,3),value,fitlabelunits(il,np,3),thr3(i)
                    endif
                  else
!
!  Subsequent lines have just potential term
!
                    write(ioout,'(28x,a15,g15.8,a10)') &
                      fitlabel(il,np,3),value,fitlabelunits(il,np,3)
                  endif
                else
!
!  Parameter free lines
!
                  if (il.eq.1) then
!
!  First line 
!
                    if (mmtexc(i).eq.0) then
                      write(ioout,'(a5,1x,a4,4x,a13,41x,f5.3,1x,f6.3)') &
                        lab1,cstyp1,name3pot(np),thr1min(i),thr1(i)
                    elseif (mmtexc(i).eq.1) then
                      write(ioout,'(a5,1x,a4,4x,a13,41x,f5.3,1x,''1 Bond'')') &
                        lab1,cstyp1,name3pot(np),thr1min(i)
                    elseif (mmtexc(i).eq.2) then
                      write(ioout,'(a5,1x,a4,4x,a13,41x,''2Bond'',1x,f6.3)') &
                        lab1,cstyp1,name3pot(np),thr1(i)
                    elseif (mmtexc(i).eq.3) then
                      write(ioout,'(a5,1x,a4,4x,a13,41x,''3Bond'',1x,f6.3)') &
                        lab1,cstyp1,name3pot(np),thr1(i)
                    elseif (mmtexc(i).eq.4) then
                      write(ioout,'(a5,1x,a4,4x,a13,41x,''1-4only'',f5.2)') &
                        lab1,cstyp1,name3pot(np),thr1(i)
                    elseif (mmtexc(i).eq.5) then
                      write(ioout,'(a5,1x,a4,4x,a13,41x,''1-4gtr '',f5.2)') &
                        lab1,cstyp1,name3pot(np),thr1(i)
                    elseif (mmtexc(i).eq.6) then
                      write(ioout,'(a5,1x,a4,4x,a13,41x,''4Bond'',1x,f6.3)') &
                        lab1,cstyp1,name3pot(np),thr1(i)
                    endif
                  elseif (il.eq.2) then
!
!  Second line
!
                    if (mmtexc(i).eq.0) then
                      write(ioout,'(a5,1x,a4,58x,f5.3,1x,f6.3)') &
                        lab2,cstyp2,thr2min(i),thr2(i)
                    elseif (mmtexc(i).eq.1) then
                      write(ioout,'(a5,1x,a4,58x,f5.3)') &
                        lab2,cstyp2,thr2min(i)
                    elseif (mmtexc(i).eq.2) then
                      write(ioout,'(a5,1x,a4,64x,f6.3)') &
                        lab2,cstyp2,thr2(i)
                    elseif (mmtexc(i).eq.3) then
                      write(ioout,'(a5,1x,a4,64x,f6.3)') &
                        lab2,cstyp2,thr2(i)
                    elseif (mmtexc(i).eq.4) then
                      write(ioout,'(a5,1x,a4,65x,f5.2)') &
                        lab2,cstyp2,thr2(i)
                    elseif (mmtexc(i).eq.5) then
                      write(ioout,'(a5,1x,a4,65x,f5.2)') &
                        lab2,cstyp2,thr2(i)
                    elseif (mmtexc(i).eq.6) then
                      write(ioout,'(a5,1x,a4,64x,f6.3)') &
                        lab2,cstyp2,thr2(i)
                    endif
                  elseif (il.eq.3) then
!
!  Third line
!
                    if (mmtexc(i).eq.0) then
                      write(ioout,'(a5,1x,a4,58x,f5.3,1x,f6.3)') &
                        lab3,cstyp3,thr3min(i),thr3(i)
                    elseif (mmtexc(i).eq.1) then
                      write(ioout,'(a5,1x,a4,58x,f5.3)') & 
                        lab3,cstyp3,thr3min(i)
                    elseif (mmtexc(i).eq.2) then
                      write(ioout,'(a5,1x,a4,64x,f6.3)') &
                        lab3,cstyp3,thr3(i)
                    elseif (mmtexc(i).eq.3) then
                      write(ioout,'(a5,1x,a4,64x,f6.3)') &
                        lab3,cstyp3,thr3(i)
                    elseif (mmtexc(i).eq.4) then
                      write(ioout,'(a5,1x,a4,65x,f5.2)') &
                        lab3,cstyp3,thr3(i)
                    elseif (mmtexc(i).eq.5) then
                      write(ioout,'(a5,1x,a4,65x,f5.2)') &
                        lab3,cstyp3,thr3(i)
                    elseif (mmtexc(i).eq.6) then
                      write(ioout,'(a5,1x,a4,64x,f6.3)') &
                        lab3,cstyp3,thr3(i)
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
    if (status/=0) call deallocate_error('outthree','iptyp')
  endif
!
  return
  end
