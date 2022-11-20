  subroutine firstpass
!
!  Performs first pass of input file before proper reading in and processing
!  begins. The tasks performed during this first pass are as follows:
!
!    (1) Check for any "keyword" option lines as it is important that all
!        keywords present are known about before other options are read.
!    (2) Write out input to channel 4 for second pass.
!
!  10/02 Start option allowed for
!  11/06 Stop command now calls gulpfinish for clean exit
!  12/08 Module input renamed to gulpinput
!   6/09 Name of getline changed to gulp_getline
!  12/09 Check added to prevent adding duplicate keywords to the keyword line
!   2/17 Check for include directive
!   1/18 Trace added
!   2/18 Trap for missing include files added
!   1/19 maxwordlength changes added
!   3/19 Modified to avoid string length overflow
!   8/19 Modified to avoid compiler warning
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
  use control
  use gulpinput
  use gulp_lengths
  use iochannels,  only : ioin
  use parallel
#ifdef TRACE
  use trace,       only : trace_in, trace_out
#endif
  implicit none
!
!  Local variables
!
  character(len=maxwordlength) :: word
  character(len=maxlinelength) :: error
  character(len=maxlinelength) :: incfile
  character(len=maxlinelength) :: line
  character(len=maxlinelength) :: linesave
  integer(i4)                  :: i
  integer(i4)                  :: iblank
  integer(i4)                  :: iline
  integer(i4)                  :: ioinc
  integer(i4)                  :: nptr
  logical                      :: lend
  logical                      :: lendinclude
#ifdef TRACE
  call trace_in('firstpass')
#endif
!
  nptr = 1
  iline = 0
!**********************************
!  Loop over lines of input file  *
!**********************************
!
10 continue
!
!  Fetch next line
!
  call gulp_getline(ioin,line,lend)
!
!  Check for end of input file
!
  if (lend) goto 100
!
!  Save line before processing
!
  linesave = line
!
!  Process line
!
  call linepronc(line,iline)
!
!  Put the words into lower case
!
  do i = 1,nword
    call stolc(words(i),maxword)
  enddo
  if (index(words(1),'help').ne.0) then
!
!  Check for help
!
    call help
    goto 10
  elseif (index(words(1),'inc').ne.0.and.nword.gt.1) then
!
!  Find name of include file as full length string
!
    call lineprown(linesave,2_i4,incfile,iline)
!
!  Open include file
!
    ioinc = 8_i4
    open(ioinc,file=incfile,status='old',form='formatted',err=20)
    lendinclude = .false.
    do while (.not.lendinclude)
!
!  Read lines from include file
!
      call gulp_getline(ioinc,line,lendinclude)
!
!  Write line
!
      if (.not.lendinclude) then
        write(4,'(a)') line
        iline = iline + 1
      endif
    enddo
    close(ioinc)
  else
!
!  Write line
!
    write(4,'(a)') line
    iline = iline + 1
  endif
!
!  Check for comment line
!
  if (index(words(1),'#').eq.1) goto 10
!
!  Check for blank line
!
  if ((nfloat+nword).eq.0) goto 10
!
!  Check for stop
!
  if (index(words(1),'stop').eq.1) call gulpfinish
!
!  Check for run
!
  if (index(words(1),'run').eq.1) goto 100
!
!  Check for start
!
  if (index(words(1),'star').eq.1) goto 100
!
!  Is this a keyword line?
!
  if (index(words(1),'key').eq.1) goto 200
!
!  End of processing of this line -> go back to read
!
  goto 10
!
!  End of reading
!
100 continue
#ifdef TRACE
  call trace_out('firstpass')
#endif
  return
!*********************
!  Get the keywords  *
!*********************
!
!  Keyword option line - add to keyword line
!
200 do i = 2,nword
    word = words(i)
!
!  Check whether word is already in keyword line - if so, don't duplicate
!
    iblank = index(word,' ')
    if (iblank.eq.0) iblank = 21
    if (index(keyword,word(1:iblank-1)).eq.0) then
      if (nptr+iblank.gt.400) then
        call outerror('keyword line has exceeded maximum length',0_i4)
        stop
      endif
      keyword(nptr:nptr+iblank-2) = word(1:iblank-1)
      nptr = nptr + iblank
    endif
  enddo
  goto 10
20 error = 'failed to open include file '//trim(incfile(1:maxlinelength-28))
  call outerror(trim(error),0_i4)
  call stopnow('firstpass')
!
  end
