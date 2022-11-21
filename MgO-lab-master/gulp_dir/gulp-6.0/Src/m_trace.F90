  module trace
!
!  Debug trace option for GULP to track the progress through the code
!
!   1/19 maxwordlength changes added
!
!  Julian Gale, CIC, Curtin University, January 2019
!
    use datatypes
    use gulp_lengths
    use parallel
    character(len=maxwordlength),         save :: tracefile
    integer(i4),                          save :: iotrace = 3_i4  ! Trace file channel
    integer(i4),                          save :: maxtracelevel = 40
    integer(i4),                          save :: ntracelevel = 0
    logical,                              save :: ltrace = .false.

  contains

    subroutine init_trace
      character(len=20),                  save :: procstring
!   
!  Called to initialise the trace
!
      if (nprocs.eq.1) then
        tracefile = 'gulp_trace'
      else
        write(procstring,'(i20)') procid
        procstring = adjustl(procstring)
        tracefile = 'gulp_trace_'//trim(procstring)
      endif
      open(iotrace,file=tracefile,form='formatted',status='unknown',err=999)
!
!  Successful open of trace file
!
      return
!
  999 call outerror('failed to open trace file',0_i4)
      call stopnow('init_trace')

    end subroutine init_trace

    subroutine trace_in(routine)
!
!  Called when entering a routine
!
      character(len=*),  intent(in) :: routine
!
      character(len=maxlinelength)  :: line
      character(len=maxwordlength)  :: routine_local
      integer(i4)                   :: nindent
!
!  Copy routine name to local string
!
      routine_local = ' '
      routine_local = trim(routine)
!
!  Set the indent level
!
      nindent = min(ntracelevel,maxtracelevel)
!
!  Prefix routine name
!
      line = ' '
      if (nindent.eq.0) then
        line(nindent+1:nindent+4) = "TOP "
      else
        line(nindent+1:nindent+4) = " -> "
      endif
!
!  Place the routine name in the output string
!
      line(nindent+5:nindent+85) = routine_local(1:80)
!
!  Write out line
!
      write(iotrace,'(a)') trim(line)
      call gflush(iotrace)
!
!  Increment the level
!
      ntracelevel = ntracelevel + 1

    end subroutine trace_in

    subroutine trace_out(routine)
!   
!  Called when exiting a routine
!
      character(len=*),  intent(in) :: routine
!
      character(len=maxlinelength)  :: line
      character(len=maxwordlength)  :: routine_local
      character(len=60)             :: errorstring
      integer(i4)                   :: nindent
!
!  Copy routine name to local string
!
      routine_local = ' '
      routine_local = trim(routine)
!
!  Decrease the level
!
      ntracelevel = ntracelevel - 1
!
!  Has the level gone too low?
!
      if (ntracelevel.lt.0) then
        errorstring = 'trace level has gone below zero in '//trim(routine_local)
        call outerror(errorstring,0_i4)
        call stopnow('trace_out')
      endif
!
!  Set the indent level
!
      nindent = min(ntracelevel,maxtracelevel)
!
!  Prefix routine name
!
      line = ' '
      if (nindent.eq.0) then
        line(nindent+1:nindent+4) = "END "
      else
        line(nindent+1:nindent+4) = " <- "
      endif
!
!  Place the routine name in the output string
!
      line(nindent+5:nindent+85) = routine_local(1:80)
!
!  Write out line
!
      write(iotrace,'(a)') trim(line)
      call gflush(iotrace)

    end subroutine trace_out

  end module trace
