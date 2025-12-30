module fileop_module
  ! This module contains implementations of low-level file operations.
  ! At present, only stat and glob are implemented.

  use, intrinsic :: iso_c_binding, only: C_INT64_T, C_DOUBLE, C_INTPTR_T, &
       C_PTR, C_NULL_CHAR, C_CHAR, C_NULL_PTR, C_ASSOCIATED
  implicit none
  private

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type, public :: stat
     ! This class is a wrapper around the underlying system calls
     ! lstat(2), stat(2) and the C type "struct stat".  It provides a
     ! Fortran-accessible interface to both.  To use it:
     !
     !   type(stat), target :: s
     !   s%stat("filename")
     !
     ! The s%lstat does an lstat instead, which stats the specified
     ! filename, even if it is a link.  The s%stat will follow the
     ! link.
     !
     ! To detect whether the stat succeeded, you can do one of these:
     !   1. Check for a .true. return from s%stat or s%lstat
     !   2. Check for s%ierr==0 .and. s%initialized=.true.
     ! 
     ! The convenience variable s%initialized is set to .true. after
     ! stat or lstat, even if the stat or lstat failed.  
     ! 
     ! All of the member variables in type stat are public but are
     ! meant to be read, not written.  Only the member functions
     ! (stat, lstat, stat_lstat_file) should modify this class's
     ! member variables.  Read on for meanings of each variable.

     ! These fields have the same meanings as in stat(2), with the
     ! "st_" part of the name removed:
     integer(kind=C_INT64_T) :: dev=-1,ino=-1,mode=0,nlink=0,uid=-1,gid=-1
     integer(kind=C_INT64_T) :: rdev=-1,size=0,blksize=-1,blocks=0
     real(kind=C_DOUBLE) :: atime=-1,mtime=-1,ctime=-1

     ! The justperm is the permission bits of the mode:
     integer(kind=C_INT64_T) :: justperm=0

     ! The ierr is used to check what error happened.  It is set to
     ! the value of errno from the call to the C stat or lstat routine.
     integer(kind=C_INT64_T) :: ierr=0 ! 0 = no error

     ! These are the return values of the corresponding POSIX macros
     ! (S_ISDIR=>isdir, S_ISREG=>isreg, etc.) for checking the type
     ! of the file:
     logical :: isreg=.false., isdir=.false., ischr=.false., isblk=.false.
     logical :: isfifo=.false., islnk=.false., issock=.false.

     ! lstat_flag - was this from an lstat, or a stat?  This is .false. if
     ! the stat followed links, and .true. if it did not.
     logical :: lstat_flag=.false.

     ! The "initialized" flag tells whether this was initialized via a
     ! stat_file or lstat_file call.  .True. means it was.  Even if
     ! stat or lstat fail, the initialized=.true.
     logical :: initialized=.false.
   contains
     procedure stat_nonempty   ! test -s
     procedure same_file => stat_same_file
     procedure clear => stat_clear
     procedure stat => stat_file
     procedure lstat => lstat_file
     procedure print => stat_print
     procedure, private :: stat_lstat_file
  end type stat

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  public :: init_glob
  type, public :: glob
     ! Globbing support, a wrapper around the C glob() and globfree()
     ! functions.  Globbing is expansion of "/dir/path/*.txt" to a
     ! list of files in "/dir/path" whose names end in ".txt".  See
     ! the glob(3) manpage for details.  

     ! Example usage:

     ! type(glob) :: globby
     ! integer :: ierr
     ! call init_glob(globby,"/dir/path/*.txt",ierr)
     ! if(ierr/=0) then
     !    write(0,*) 'could not find /dir/path/*.txt'
     !    stop 9
     ! endif
     ! do while(globby%next())   ! note the initial next() before accessing data
     !  7 format('Found .txt file: "',A,'"')
     !    print 7,globby%match
     ! enddo
     ! call globby%free()
     type(C_PTR), private :: pglob
     character(len=:), pointer, private :: buffer=>NULL()
     logical, private :: quiet, reverse

     ! NOTE: The public variables are intended to be read-only.
     ! Changing them will modify the internal state of the class,
     ! possibly in inconsistent ways.  Use the subroutines, in the
     ! "contains" section, to modify this object.

     integer(kind=C_INT64_T) :: idx=-999, ierr=-999, len=0
     character(len=:), pointer :: match=>NULL()

     ! match - a matching filename.  Do NOT use trim.
     ! len - the number of matches
     ! idx - the index of this match, in the list of matching filenames
     ! ierr - last error status
     
   contains
     procedure rewind => glob_rewind ! reset the match iterator

     procedure set_order => glob_set_reverse ! sets the reverse flag
     procedure get_order => glob_get_reverse ! gets the reverse flag
     ! reverse - .true. = glob%rewind will jump to the last element,
     !           .false. = glob%rewind will jump to the first element

     ! NOTE about rewind: if reverse=.true. in init_glob, then this will
     ! jump to the last index.  If reverse=.false. or if reverse was not
     ! given, then it will jump to the first index.

     procedure jump => glob_jump ! jump to specific match index
     procedure next => glob_next ! go to next match index, return false at end
     procedure prev => glob_prev ! go to previous match, return false at first

     procedure free => glob_free ! free used resources
  end type glob

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Glob error return statuses.  NOTE: These must match the F_
  ! equivalents in c_glob.c
  !   GLOB_NOSPACE - cannot allocate memory
  !   GLOB_ABORTED - globbing aborted due to I/O error, and GLOB_ERR
  !     (err=.true.) was specified.
  !   GLOB_NOMATCH - glob had no match
  ! These last two should not be possible:
  !   GLOB_UNKNOWN_ERROR - something went wrong, but reason was unknown
  integer(kind=C_INT64_T), public, parameter :: &
       GLOB_NOSPACE=1, GLOB_ABORTED=2, GLOB_NOMATCH=4, &
       GLOB_UNKNOWN_ERROR=8, GLOB_INVALID_INDEX=16

  ! GLOB_ERR, etc.: Flags used internally to pass options to
  ! c_glob_start:
  integer(kind=C_INT64_T), parameter :: &
       GLOB_ERR=1, GLOB_MARK=2, GLOB_NOSORT=4, GLOB_DOOFFS=8, &
       GLOB_PERIOD=16, GLOB_BRACE=32, GLOB_TILDE=64, &
       GLOB_TILDE_CHECK=128, GLOB_ONLYDIR=256, GLOB_QUIET=512, &
       GLOB_REVERSE_ITERATE=1024

  integer, parameter :: GLOB_BUFFER_PADDING = 300
  ! Number of bytes to pad string buffers on the right, to allow for
  ! longer strings in future glob matches.

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  interface 
     subroutine c_glob_start(pattern,flags,pglob,idx,offs,len,ierr) bind(C)
       use iso_c_binding
       implicit none
       type(C_PTR) :: pglob
       integer(kind=C_INT64_T), intent(in) :: flags
       integer(kind=C_INT64_T), intent(out) :: idx,offs,len,ierr
       character(kind=C_CHAR), intent(in) :: pattern(*)
     end subroutine c_glob_start
  end interface

  interface 
     subroutine c_glob_len(pglob,idx,strlen) bind(C)
       use iso_c_binding
       implicit none
       type(C_PTR), value :: pglob
       integer(kind=C_INT64_T), intent(in) :: idx
       integer(kind=C_INT64_T), intent(out) :: strlen
     end subroutine c_glob_len
  end interface

  interface 
     subroutine c_glob_fill(pglob,idx,strlen,buf) bind(C)
       use iso_c_binding
       implicit none
       type(C_PTR), value :: pglob
       integer(kind=C_INT64_T), intent(in) :: idx,strlen
       character(len=1), intent(out) :: buf(*)
     end subroutine c_glob_fill
  end interface

  interface
     subroutine c_glob_free(pglob) bind(C)
       use iso_c_binding
       implicit none
       type(C_PTR) :: pglob
     end subroutine c_glob_free
  end interface
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! ####################################################################
  ! ## STAT CLASS IMPLEMENTATION #######################################
  ! ####################################################################

  logical function stat_same_file(s,t)
    ! Are these two stats for the same filesystem object?  Returns
    ! false if either stat is invalid, or if either stat failed.
    class(stat) :: s,t
    if(s%ierr==0 .and. t%ierr==0 .and. s%initialized .and. t%initialized) then
       stat_same_file = s%dev==t%dev .and. s%ino==t%ino
    else
       ! One or both stats failed.  Assume different files.
       stat_same_file = .false.
    endif
  end function stat_same_file

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine stat_clear(s)
    ! Resets all values to invalid data, and sets initialized=.false.
    class(stat) :: s
    s%dev=-1
    s%ino=-1
    s%mode=0
    s%nlink=0
    s%uid=-1
    s%gid=-1
    s%rdev=-1
    s%size=0
    s%blksize=-1
    s%blocks=0
    s%atime=-1
    s%mtime=-1
    s%ctime=-1
    s%justperm=0
    s%ierr=0
    s%isreg=.false.
    s%isdir=.false.
    s%ischr=.false.
    s%isblk=.false.
    s%isfifo=.false.
    s%islnk=.false.
    s%issock=.false.
    s%initialized=.false.
  end subroutine stat_clear

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine stat_print(s,unit,filename)
    use datetime_module, only: grads_date, grads_date_len
    ! Prints a textual representation of the stat buffer to the
    ! specified output unit, or to unit 6 (stdout) if unspecified.
    ! The optional filename argument is printed in the first line of
    ! output (Stat "filename") if given.
    class(stat) :: s
    integer, intent(in), optional :: unit
    character*(*), intent(in), optional :: filename
    character(len=grads_date_len) :: atime,ctime,mtime
    integer :: iunit

    if(present(unit)) then
       iunit=unit
    else
       iunit=6
    endif

10  format('Stat "',A,'"')
20  format('Perm 0',O0,' times: a=',F0.1,' m=',F0.1,' c=',F0.1)
25  format('Perm 0',O0,' times: a="',A,'" m="',A,'" c="',A,'"')
30  format('Dev=',I0,' Ino=',I0,' Nlink=',I0,' Uid=',I0,' Gid=',I0)
40  format('Size=',I0,' Blocks=',I0,' Block Size=',I0)
50  format(A)
55  format(' ',A)

    call grads_date(atime,s%atime)
    call grads_date(mtime,s%mtime)
    call grads_date(ctime,s%ctime)

    if(present(filename)) then
       write(iunit,10) trim(filename)
    endif
    write(iunit,20) s%justperm,s%atime,s%mtime,s%ctime
    write(iunit,25) s%justperm,atime,mtime,ctime
    write(iunit,30) s%dev,s%ino,s%nlink,s%uid,s%gid
    write(iunit,40) s%size,s%blocks,s%blksize

    write(iunit,50,advance='no') 'Flags: ['
    if(s%isreg)  write(iunit,55,advance='no') 'REG'
    if(s%isdir)  write(iunit,55,advance='no') 'DIR'
    if(s%ischr)  write(iunit,55,advance='no') 'CHR'
    if(s%isblk)  write(iunit,55,advance='no') 'BLK'
    if(s%isfifo) write(iunit,55,advance='no') 'FIFO'
    if(s%islnk)  write(iunit,55,advance='no') 'LNK'
    if(s%issock) write(iunit,55,advance='no') 'SOCK'
    write(iunit,55) ']'
  end subroutine stat_print

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  logical function stat_nonempty(s)
    ! Returns .true. if the statted file is non-empty and exists, or
    ! .false. otherwise.  Returns .false. if the object is
    ! uninitialized.
    class(stat) :: s
    stat_nonempty = s%size>0 .and. s%isreg .and. s%initialized &
         .and. s%ierr==0
  end function stat_nonempty

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  logical function stat_file(s,filename,quiet)
    ! Requests a stat(2) call that follows symbolic links.  The actual
    ! work is done by the private stat_lstat_file routine.  Set
    ! quiet=.true. to disable warning messages.
    class(stat) :: s
    character*(*), intent(in) :: filename
    logical, intent(in), optional :: quiet
    logical :: pquiet
    pquiet=.false.
    if(present(quiet)) pquiet=quiet
    stat_file=s%stat_lstat_file(filename,len(filename),.false.,pquiet)
  end function stat_file

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  logical function lstat_file(s,filename,quiet)
    ! Requests a stat(2) call that does NOT follow symbolic links.
    ! The actual work is done by the private stat_lstat_file routine.
    ! Set quiet=.true. to disable warning messages.
    class(stat) :: s
    character*(*), intent(in) :: filename
    logical, intent(in), optional :: quiet
    logical :: pquiet
    pquiet=.false.
    if(present(quiet)) pquiet=quiet
    lstat_file=s%stat_lstat_file(filename,len(filename),.true.,pquiet)
  end function lstat_file

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  logical function stat_lstat_file(s,filename,flen,lstat,quiet)
    ! Calls the c_stat function in c_stat.c, which calls the stat or
    ! lstat routines.  Set quiet=.true. to disable c_stat's warning
    ! messages.  
    class(stat) :: s
    character*(*), intent(in) :: filename
    integer, intent(in) :: flen
    logical, intent(in) :: lstat,quiet
    ! locals
    character(len=flen+1) :: cfilename
    integer(kind=C_INT64_T) :: iflags,ilstat,iquiet

    s%lstat_flag=lstat
    s%initialized=.true.
    ilstat=0
    if(lstat) ilstat=1
    iquiet=0
    if(quiet) iquiet=1
    s%ierr=-999
    iflags=0
    
    cfilename(1:flen)=filename(1:flen)
    cfilename(flen+1:flen+1)=char(0)

    call c_stat(cfilename,ilstat,s%ierr,iquiet,&
         s%dev, s%ino, s%mode, s%justperm, s%nlink, s%uid, s%gid, s%rdev, &
         s%size, s%blksize, s%blocks, s%atime, s%mtime, s%ctime, iflags)

    stat_lstat_file = (s%ierr==0)

    s%isreg  = ( 0 /= iand(iflags,1) )
    s%isdir  = ( 0 /= iand(iflags,2) )
    s%ischr  = ( 0 /= iand(iflags,4) )
    s%isblk  = ( 0 /= iand(iflags,8) )
    s%isfifo = ( 0 /= iand(iflags,16) )
    s%islnk  = ( 0 /= iand(iflags,32) )
    s%issock = ( 0 /= iand(iflags,64) )
  end function stat_lstat_file

  ! ####################################################################
  ! ## GLOB CLASS IMPLEMENTATION #######################################
  ! ####################################################################

  subroutine init_glob(G,pattern,ierr,offs,&
                       opt_err,opt_mark,opt_nosort,&
                       gnu_period,gnu_brace,gnu_tilde,gnu_tilde_check,&
                       gnu_onlydir,&
                       quiet,reverse)
    class(glob) :: G
    character*(*) :: pattern
    integer, intent(out) :: ierr

    ! Optional arguments:
    logical, intent(in), optional :: opt_err,opt_mark,opt_nosort ! glob options
    integer, intent(in), optional :: offs ! guess for number of matches
    logical, intent(in), optional :: quiet, reverse ! implementation flags

    ! GNU glob options, which may not be supported on lower-quality
    ! (ie.: commercial) unix implementations:
    logical, intent(in), optional :: gnu_period,gnu_brace,gnu_tilde,&
         gnu_tilde_check,gnu_onlydir

    ! Locals
    integer(kind=C_INT64_T) :: c_offs,c_flags
    logical :: first

    c_flags=0
    c_offs=0
    
    ! Add POSIX glob options to c_flags
    if(present(opt_err)) then
       if(opt_err)          c_flags=ior(c_flags,GLOB_ERR)
    endif
    if(present(opt_mark)) then
       if(opt_mark)         c_flags=ior(c_flags,GLOB_MARK)
    endif
    if(present(opt_nosort)) then
       if(opt_nosort)       c_flags=ior(c_flags,GLOB_NOSORT)
    endif

    ! If requested, add POSIX GLOB_DOOFFS option, and set c_offs.
    if(present(offs)) then
       c_offs=max(1,offs)
       c_flags=ior(c_flags,GLOB_DOOFFS)
    endif

    ! Add GNU glob options to c_flags
    if(present(gnu_period)) then
       if(gnu_period)       c_flags=ior(c_flags,GLOB_PERIOD)
    endif
    if(present(gnu_brace)) then
       if(gnu_brace)        c_flags=ior(c_flags,GLOB_BRACE)
    endif
    if(present(gnu_tilde)) then
       if(gnu_tilde)        c_flags=ior(c_flags,GLOB_TILDE)
    endif
    if(present(gnu_tilde_check)) then
       if(gnu_tilde_check)  c_flags=ior(c_flags,GLOB_TILDE_CHECK)
    endif
    if(present(gnu_onlydir)) then
       if(gnu_onlydir)      c_flags=ior(c_flags,GLOB_ONLYDIR)
    endif

    ! Add flags specific to this implementation
    G%reverse=.false.
    if(present(reverse)) then
       if(reverse)          G%reverse=.true.
    endif
    G%quiet=.false.
    if(present(quiet)) then
       if(quiet) then
          c_flags=ior(c_flags,GLOB_QUIET)
          G%quiet=.true.
       endif
    endif

    call c_glob_start(pattern//C_NULL_CHAR,c_flags,G%pglob,G%idx,c_offs,&
                      G%len,G%ierr)
    ierr=G%ierr
  end subroutine init_glob

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine glob_set_reverse(G,reverse)
    class(glob) :: G
    logical, intent(in) :: reverse
    G%reverse=reverse
    ! nothing else to do here yet...
  end subroutine glob_set_reverse

  logical function glob_get_reverse(G)
    class(glob) :: G
    glob_get_reverse=G%reverse
    ! nothing else to do here yet...
  end function glob_get_reverse

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine glob_jump(G,idx)
    class(glob) :: G
    integer(kind=C_INT64_T), intent(in) :: idx
    ! Locals
    integer(kind=C_INT64_T) :: strlen
    integer :: ierr
    logical :: need_alloc

    if(idx<1 .or. idx>G%len) then
       if(.not.G%quiet) then
202       format('ERROR: invalid index ',I0,' specified to glob%jump.  ',&
                 'Must be in [',I0,'...',I0,'].')
          write(0,202) idx,1,G%len
       endif
       G%ierr=GLOB_INVALID_INDEX
    endif

    ! Set the current index.
    G%idx=idx
    nullify(G%match)

    ! Get the length of the match:
    call c_glob_len(G%pglob,idx,strlen)
    if(strlen<1) then
       if(.not.G%quiet) then
303       format('ERROR: Invalid index ',I0,' specified to glob%jump:',&
               'c_glob_len returned invalid strlen=',I0,'.')
          write(0,303) idx,strlen
       endif
       G%ierr=GLOB_INVALID_INDEX
       nullify(G%match)
       return
    endif

    ! Make sure we have room for it in the buffer:
    need_alloc=.false.
    if(.not.associated(G%buffer)) then
       need_alloc=.true.
    elseif(len(G%buffer)>strlen) then
       need_alloc=.true.
       deallocate(G%buffer)
       nullify(G%buffer)
    endif

    if(need_alloc) then
       allocate(character(len=strlen+GLOB_BUFFER_PADDING) :: G%buffer,&
            stat=ierr)
       if(ierr/=0) then
          if(.not.G%quiet) then
404          format('ERROR: cannot allocate ',I0,' bytes for a glob match.')
             write(0,404) strlen+GLOB_BUFFER_PADDING
          endif
          G%ierr=GLOB_NOSPACE
          return
       endif
    endif

    ! Copy the data to the buffer
    call c_glob_fill(G%pglob,idx,strlen,G%buffer)

    ! Point the match to the used portion of the buffer:
    G%match => G%buffer(1:strlen)
  end subroutine glob_jump

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Aliases of glob_jump for common cases

  logical function glob_next(G)
    class(glob) :: G
    if(G%idx==G%len) then
       glob_next=.false.
    else
       call G%jump(G%idx+1)
       glob_next=.true.
    endif
  end function glob_next

  logical function glob_prev(G)
    class(glob) :: G
    if(G%idx==1) then
       glob_prev=.false.
    elseif(G%idx==0) then
       ! Special case: just initialized.  Have no index yet.
       call G%jump(G%len)
       glob_prev=.true.
    else
       call G%jump(G%idx-1)
       glob_prev=.true.
    endif
  end function glob_prev

  subroutine glob_rewind(G)
    class(glob) :: G
    logical :: discard
    G%idx=0
    nullify(G%match)
    if(G%reverse) then
       discard=G%prev()
    else
       discard=G%next()
    endif
  end subroutine glob_rewind

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine glob_free(G)
    class(glob) :: G
    if(C_ASSOCIATED(G%pglob)) then
       call c_glob_free(G%pglob)
       G%pglob=C_NULL_PTR
    endif
    if(associated(G%buffer)) then
       deallocate(G%buffer)
    endif
    nullify(G%buffer, G%match)
    G%idx=-999
    G%ierr=-999
    G%len=0
  end subroutine glob_free

end module fileop_module
