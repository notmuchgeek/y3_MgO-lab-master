  subroutine worsy(word,lsymbol,lpotcall)
!
!  Word OR SYmbol - if word is a valid atomic symbol then return
!  with lsymbol = .true. else return with lsymbol = .false.
!  This routine is used during input to decide when the end of
!  a list of atom specifications has been reached.
!  If lpotcall is .true. then this is a call relating to a potential
!  and therefore wildcards are allowed. If called elsewhere, wildcards
!  are rejected.
!
!  10/02 Style updated
!  10/02 Wildcard atom allowed for
!   1/05 Dummy atom X now allowed for structural input
!   9/13 Algorithm for checking letters and numbers corrected
!   2/15 Trap for il > 2 added
!   1/19 Worsy changed to use word length parameter
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
!  Julian Gale, CIC, Curtin University, January 2019
!
  use element
  use gulp_lengths
  use library
  use species
  implicit none
!
!  Passed variables
!
  character(len=maxwordlength), intent(inout) :: word
  logical,                      intent(in)    :: lpotcall
  logical,                      intent(out)   :: lsymbol
!
!  Local variables
!
  character(len=2)                            :: lab
  character(len=maxwordlength)                :: wordin
  integer(i4)                                 :: i
  integer(i4)                                 :: i0
  integer(i4)                                 :: ial
  integer(i4)                                 :: iau
  integer(i4)                                 :: ic
  integer(i4)                                 :: il
  integer(i4)                                 :: ilen
  integer(i4)                                 :: ino
  integer(i4)                                 :: ishift
  integer(i4)                                 :: maxsymbol
  logical                                     :: lblank
!
  i0  = ichar('0')
  ial = ichar('a')
  iau = ichar('A')
  ishift = iau - ial
  lsymbol = .false.
!***********************************
!  Check first for library symbol  *
!***********************************
  do i = 1,nspec
    if (libspec(i).eq.word(1:16)) lsymbol = .true.
  enddo
  if (lsymbol) return
!
!  Save word as input for return
!
  wordin = word
!*****************
!  Check format  *
!*****************
!
!  Check length
!
  lblank = .false.
  i = 0
  do while (.not.lblank.and.i.lt.maxwordlength)
    i = i + 1
    if (word(i:i).eq.' ') lblank = .true.
    if (word(i:i).eq.'_') lblank = .true.
  enddo
  ilen = i - 1
  if (ilen.eq.0) return
  if (ilen.gt.5) return
!************************
!  Conventional symbol  *
!************************
!
!  Check number of characters and numbers and order
!
!  - Make sure there are no more than 2 letters in a row at the start
!
  lblank = .false.     ! Note that here lblank is being used as a flag to detect an invalid symbol
  il = 0
  ino = 0
  do i = 1,ilen
    ic = ichar(word(i:i))
    if ((ic-i0).lt.10.and.(ic-i0).ge.0) then
      ino = ino + 1
    elseif ((ic-ial).lt.26.and.(ic-ial).ge.0) then
      il = il + 1
      if (i.eq.3.and.il.eq.3) lblank = .true.
      if (il.eq.1) then
        word(i:i) = char(ic+ishift)
      endif
    elseif ((ic-iau).lt.26.and.(ic-iau).ge.0) then
      il = il + 1
      if (i.eq.3.and.il.eq.3) lblank = .true.
      if (il.eq.2) then
        word(i:i) = char(ic-ishift)
      endif
    endif
  enddo
!
!  Need to make sure that il doesn't exceed 2
!
  il = min(il,2_i4)
  if (lblank) then
    word = wordin
    return
  endif
!
!  Finally check symbol against known list
!
!  If symbol is X then this is allowed as a wildcard if lpotcall is true
!  Now X is allowed for structure too as this can be useful
!
  lab(1:2) = '  '
  lab(1:il) = word(1:il)
  if (lpotcall) then
    maxsymbol = maxele
  else
    maxsymbol = maxele 
! This is the old line if X is not allowed in the structure
!        maxsymbol = maxele - 1
!
  endif
  i = 0
  do while (.not.lsymbol.and.i.lt.maxsymbol)
    i = i + 1
    lsymbol = (lab.eq.atsym(i))
  enddo
!
!  If this isn't a symbol then reset word to original form
!
  if (.not.lsymbol) then
    word = wordin
  endif
!
  return
  end
