  subroutine setatomdistribution(mode)
!
!  Selects which parallel distribution of atoms should be used
!
!  mode = 'a' or 'A' => atoms
!  mode = 'v' or 'V' => variables
!
!   2/17 Created from setatomnoded2
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
!  Julian Gale, CIC, Curtin University, February 2017
!
  use current
  use parallel
  implicit none
!
!  Passed variables
!
  character(len=*),    intent(in)  :: mode
!
!  Local variables
!
  logical                          :: latomcentric
!
  latomcentric = .true.
  if (index(mode,'V').eq.1.or.index(mode,'v').eq.1) latomcentric = .false.
!
  if (latomcentric) then
    natomsonnode = natomsonnodea
  else
    natomsonnode = natomsonnodev
  endif
!
!  Update maxatloc
!
  if (natomsonnode.gt.maxatloc) then
    maxatloc = natomsonnode 
    call changemaxatloc
  endif
  if (latomcentric) then
    natomsonnode = natomsonnodea
    node2atom(1:natomsonnode) = node2atoma(1:natomsonnode)
    atom2local(1:numat) = atom2locala(1:numat)
    atom2node(1:numat) = atom2nodea(1:numat)
  else
    natomsonnode = natomsonnodev
    node2atom(1:natomsonnode) = node2atomv(1:natomsonnode)
    atom2local(1:numat) = atom2localv(1:numat)
    atom2node(1:numat) = atom2nodev(1:numat)
  endif
!
  return
  end
