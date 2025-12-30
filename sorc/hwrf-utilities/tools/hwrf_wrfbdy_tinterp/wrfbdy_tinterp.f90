module module_wrfbdy_tinterp
  use netcdf
  implicit none

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! External C functions for time manipulation from etime.c
  interface
     subroutine abort() bind(c)
     end subroutine abort
  end interface

  interface
     subroutine c_etimeinit(ierr)
       integer(kind=4) :: ierr
     end subroutine c_etimeinit
  end interface

  interface
     subroutine c_etime(ymdh,epochtime)
       integer(kind=4) :: ymdh
       integer(kind=8) :: epochtime
     end subroutine c_etime
  end interface

  interface
     subroutine c_wrf2epoch(wrftime,epochtime,ierr)
       integer(kind=4) :: ierr
       integer(kind=8) :: epochtime
       character(len=19) :: wrftime
     end subroutine c_wrf2epoch
  end interface

  interface
     subroutine c_epoch2wrf(wrftime,epochtime,ierr)
       integer(kind=4) :: ierr
       integer(kind=8) :: epochtime
       character(len=19) :: wrftime
     end subroutine c_epoch2wrf
  end interface

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type wrfbdy_var_info
     character(len=:), pointer :: name=>NULL()
     integer :: ndims=-1,natts=-1,xtype=-1
     integer :: time_dim=-1 !! if >0, index in dims_in() of the time dimension
     integer(kind=8) :: nelem=-1
     integer, pointer :: dims_in(:)=>NULL(), dims_out(:)=>NULL()
     integer, pointer :: count(:)=>NULL(), start(:)=>NULL()

     integer(kind=8), pointer :: i8in(:)=>NULL()
     integer(kind=4), pointer :: i4in(:)=>NULL()
     real(kind=8), pointer :: r8in(:)=>NULL()
     real(kind=4), pointer :: r4in(:)=>NULL()
     integer(kind=8), pointer :: i8out(:)=>NULL()
     integer(kind=4), pointer :: i4out(:)=>NULL()
     real(kind=8), pointer :: r8out(:)=>NULL()
     real(kind=4), pointer :: r4out(:)=>NULL()
  end type wrfbdy_var_info

  type wrfbdy_tinterp
     character(len=19), pointer :: times_in(:)=>NULL(), times_out(:)=>NULL()
     character(len=19) :: new_time0
     integer(kind=8), pointer :: epochs_in(:)=>NULL(), epochs_out(:)=>NULL()
     integer, pointer :: varids_in(:)=>NULL(), varids_out(:)=>NULL()
     type(wrfbdy_var_info), pointer :: vars(:)=>NULL()

     integer :: ncid_in=-1, ncid_out=-1, inf_id_in=-1, inf_id_out=-1
     integer :: ndims=-1, nvars=-1, in_format=-1, ntimes_out=-1, natts=-1
     integer :: itime_in=0, itime_out=0, ntimes_in=-1
     real :: dt=-1
     integer(kind=8) :: bdystep=-1

     integer, pointer :: dims_in(:)=>NULL(), dims_out(:)=>NULL()
     integer, pointer :: dimlen(:)=>NULL()
     character(len=:), pointer :: infile=>NULL()
     character(len=:), pointer :: outfile=>NULL()
   contains
     ! Public functions
     procedure initial_io
     procedure interp
     procedure close

     ! Protected (implementation) functiosn
     procedure nfcall
     procedure fail
     procedure interp_bdy_time_loop
     procedure generate_output_times
     procedure dimlen_in
     procedure bdyopen
     procedure get_info
     procedure define_file
     procedure interp_bdy_time
     procedure find_var
     procedure read_time
     procedure copy_apply_time
     procedure copy_dims
     procedure create_vars
     procedure create_var
     procedure copy_global_attr
     procedure copy_var_attr
  end type wrfbdy_tinterp
  
contains !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Public functions

  integer function initial_io(this,infile,outfile,time0)
    implicit none
    class(wrfbdy_tinterp) :: this
    character(len=*), intent(in) :: infile,outfile
    character(len=19), intent(in) :: time0

    initial_io=0
    initial_io=this%bdyopen(infile,outfile)
    if(initial_io/=0) return

    initial_io=this%get_info()
    if(initial_io/=0) return

    this%new_time0=time0
    initial_io=this%define_file()
    if(initial_io/=0) return

    initial_io=this%nfcall('finish defining the NetCDF file',this%outfile,&
         nf90_enddef(this%ncid_out))
  end function initial_io

  ! --------------------------------------------------------------------

  integer function interp(this,time0,bdystep,ntimes_out)
    implicit none
    class(wrfbdy_tinterp) :: this
    character(len=19), intent(in) :: time0
    integer(kind=8),intent(in) :: bdystep
    integer, intent(in) :: ntimes_out

    interp=0
    interp=this%generate_output_times(time0,bdystep,ntimes_out)
    if(interp/=0) return

    interp=this%interp_bdy_time_loop()
    if(interp/=0) return
  end function interp

  ! --------------------------------------------------------------------

  integer function close(this)
    implicit none
    class(wrfbdy_tinterp) :: this
    close=0
    close=this%nfcall('close input file',this%infile,nf90_close(this%ncid_in))
    if(close/=0) return
    close=this%nfcall('close output file',this%outfile,nf90_close(this%ncid_out))
  end function close

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Protected (implementation) functions

  integer function nfcall(this,what,file,status,varid,varname,&
                          dimid,dimname,attnum,attname)
    implicit none
    class(wrfbdy_tinterp) :: this
    character(len=*), intent(in) :: what,file
    integer, intent(in) :: status
    character(len=*), intent(in), optional :: varname,dimname,attname
    integer, intent(in), optional :: varid,dimid,attnum

    integer :: ierr
    character(len = 80)  :: strerror
    character(len=len(what)+len(file)+100) :: failbuf

    nfcall=status
    if(status==0) return

    
10  format(A,': NetCDF error #',I0,'.  Operation: ',A)
20  format('NetCDF says: "',A,'"')
30  format('Variable name "',A,'"')
35  format('Variable ID #',I0)
40  format('Dimension name "',A,'"')
45  format('Dimension ID #',I0)
50  format('Attribute name "',A,'"')
55  format('Attribute number #',I0)
    print 10,trim(file),status,trim(what)
    print 20,trim(nf90_strerror(status))
    if(present(varname)) print 30,trim(varname)
    if(present(varid)) print 35,varid
    if(present(dimname)) print 40,trim(dimname)
    if(present(dimid)) print 45,dimid
    if(present(attname)) print 50,attname
    if(present(attnum)) print 55,attnum
    write(failbuf,10) trim(file),status,trim(what)
    call this%fail(trim(failbuf))
  end function nfcall

  ! --------------------------------------------------------------------

  subroutine fail(this,why)
    implicit none
    class(wrfbdy_tinterp) :: this
    character(len=*), intent(in) :: why
    write(0,'("ABORT: ",A)') why
    !call abort
    stop 9
    return
  end subroutine fail

  ! --------------------------------------------------------------------

  integer function interp_bdy_time_loop(this)
    implicit none
    class(wrfbdy_tinterp) :: this
    integer :: itime
    interp_bdy_time_loop=0
    do itime=1,this%ntimes_out
       interp_bdy_time_loop=this%interp_bdy_time(this%epochs_out(itime),&
            this%epochs_out(itime)+this%bdystep)
       if(interp_bdy_time_loop/=0) return
    enddo
  end function interp_bdy_time_loop

  ! --------------------------------------------------------------------

  integer function generate_output_times(this,time0,bdystep,ntimes_out)
    implicit none
    class(wrfbdy_tinterp) :: this
    character(len=19), intent(in) :: time0
    integer(kind=8),intent(in) :: bdystep
    integer, intent(in) :: ntimes_out
    character*180 :: message
    integer :: itime
    integer(kind=8) :: epochtime

    this%bdystep=bdystep
    this%ntimes_out=ntimes_out
    allocate(this%times_out(ntimes_out))
    allocate(this%epochs_out(ntimes_out))

    if(time0/=this%new_time0) then
10     format('New time0 sent to interp ("',A19,'") does not match initial_io ("',A19,'")')
       write(message,10) time0,this%new_time0
       call this%fail(trim(message))
       generate_output_times=-1
       return
    endif

    generate_output_times=0
    call c_wrf2epoch(time0,epochtime,generate_output_times)
    if(generate_output_times/=0) return
    
    do itime=1,ntimes_out
       this%epochs_out(itime)=epochtime
       call c_epoch2wrf(this%times_out(itime),epochtime,generate_output_times)
       if(generate_output_times/=0) return
       epochtime=epochtime+bdystep
    enddo
  end function generate_output_times

  ! --------------------------------------------------------------------

  integer function dimlen_in(this,dimid_in)
    implicit none
    class(wrfbdy_tinterp) :: this
    integer, intent(in) :: dimid_in
    integer :: dimidx

    do dimidx=1,this%ndims
       if(dimid_in==this%dims_in(dimidx)) then
          dimlen_in=this%dimlen(dimidx)
          return
       endif
    enddo
    dimlen_in=-1
  end function dimlen_in

  ! --------------------------------------------------------------------

  integer function bdyopen(this,infile,outfile)
    implicit none
    character(len=*), intent(in) :: infile,outfile
    class(wrfbdy_tinterp) :: this
    integer :: ierr

    allocate(character(len=len(infile)) :: this%infile)
    this%infile=infile

    bdyopen=this%nfcall('open for read',infile,nf90_open( &
         infile,NF90_NOWRITE,this%ncid_in))
    if(bdyopen/=0) return

    bdyopen=this%nfcall('inquire input file info',infile,this%get_info())
    if(bdyopen/=0) return

    allocate(character(len=len(outfile)) :: this%outfile)
    this%outfile=outfile
    
    bdyopen=this%nfcall('create new file',outfile,nf90_create( &
         outfile,this%in_format,this%ncid_out,0))
    if(bdyopen/=0) return
  end function bdyopen

  ! --------------------------------------------------------------------

  integer function get_info(this)
    implicit none
    class(wrfbdy_tinterp) :: this
    integer :: timesid

    get_info=this%nfcall('inquire global info',this%infile,nf90_inquire( &
       this%ncid_in,this%ndims,this%nvars,this%natts,this%inf_id_in,&
       this%in_format))
    if(get_info/=0) return

    allocate(this%dims_in(this%ndims), this%dims_out(this%ndims), &
         this%dimlen(this%ndims))    
  end function get_info

  ! --------------------------------------------------------------------

  integer function define_file(this)
    implicit none
    class(wrfbdy_tinterp) :: this

    define_file=this%copy_dims()
    if(define_file/=0) return

    define_file=this%create_vars()
    if(define_file/=0) return

    define_file=this%copy_global_attr()
    !if(define_file/=0) return

  end function define_file

  ! --------------------------------------------------------------------

  integer function interp_bdy_time(this,epoch,next_epoch)
    implicit none
    class(wrfbdy_tinterp) :: this
    integer(kind=8) :: epoch, next_epoch
    character(len=19) :: nowtime,nexttime ! output values
    character(len=19) :: nowtime_in,nexttime_in ! inputs, ignored
    integer :: itime
    real :: tendfact

    tendfact=-999
    itime=1
    interp_bdy_time=0
    if(epoch>=this%epochs_in(1)) then
       ! Time is at or after simulation start.  Is it between two bdy
       ! times?
       timeloop: do itime=1,this%ntimes_in-1
          if(epoch >= this%epochs_in(itime) .and. &
               epoch <= this%epochs_in(itime+1)) then
             ! Found the bounding bdy times, so exit loop to set tendfact.
             exit timeloop
          endif
       enddo timeloop
    endif
5   format('Time ',I0,' epoch=',I0,' relative to ',I0,' with dt=',F0.3,' tendfact=',F0.3)
    tendfact = real(epoch-this%epochs_in(itime)) !/ this%dt (docs were wrong)
    print 5,itime,epoch,this%epochs_in(itime),this%dt,tendfact

    if(this%itime_in/=itime) then
       ! We need to read in the data for this time.
       interp_bdy_time=this%read_time(itime,nowtime_in,nexttime_in)
       if(interp_bdy_time/=0) return
10     format('Input:  ',A,' => ',A)
       print 10,nowtime_in,nexttime_in
    endif

    call c_epoch2wrf(nowtime,epoch,interp_bdy_time)
    if(interp_bdy_time/=0) return

    call c_epoch2wrf(nexttime,next_epoch,interp_bdy_time)
    if(interp_bdy_time/=0) return

20  format('Output: ',A,' => ',A,' tendfact=',F0.7)

    print 20,nowtime,nexttime,tendfact

    this%itime_out=this%itime_out+1
    interp_bdy_time=this%copy_apply_time(tendfact,nowtime,nexttime,this%itime_out)
  end function interp_bdy_time

  ! --------------------------------------------------------------------

  subroutine find_var(this,ivar, bdyvar, tendvar)
    implicit none
    class(wrfbdy_tinterp) :: this
    integer, intent(in) :: ivar
    integer, intent(inout) :: bdyvar,tendvar

    character(len=:), pointer :: iname,bname,tname
    integer :: leni,i
    
    iname=>this%vars(ivar)%name
    leni=len(iname)
    bdyvar=-1
    tendvar=-1

10  format('find_var(',I0,'="',A,'") look for bdy and tend...')
    !write(0,10) ivar,iname

    if(iname(leni-3:leni-2)=='_B') then
       ! This is a boundary variable
       allocate(character(len=leni+1) :: tname)
       tname(1:leni-2)=iname(1:leni-2)
       tname(leni-1:leni-1)='T'
       tname(leni:leni+1)=iname(leni-1:leni)
20     format('find_var(',I0,'="',A,'") = is bdy var.  Look for tend var "',A,'"...')
       !write(0,20) ivar,iname,tname
       bdyvar=ivar
       do i=1,this%nvars
          if(this%vars(i)%name == tname) then
             tendvar=i
25           format('find_var(',I0,'="',A,'") found tend var "',A,'" at ',I0,'.')
             !print 25,ivar,iname,tname,tendvar
             deallocate(tname)
             return
          endif
       enddo
    elseif(iname(leni-4:leni-2)=='_BT') then
       ! This is a boundary tendency variable
       allocate(character(len=leni) :: bname)
       bname(1:leni-3)=iname(1:leni-3)
       bname(leni-2:leni-1)=iname(leni-2:leni-1)
       tendvar=ivar
30     format('find_var(',I0,'="',A,'") = is tend var.  Look for bdy var "',A,'"...')
       !write(0,20) ivar,iname,bname
       do i=1,this%nvars
          if(this%vars(i)%name == bname) then
             bdyvar=i
35           format('find_var(',I0,'="',A,'") found bdy var "',A,'" at ',I0,'.')
             !print 35,ivar,iname,bname,bdyvar
             deallocate(bname)
             return
          endif
       enddo
       deallocate(bname)
    endif
  end subroutine find_var

  ! --------------------------------------------------------------------

  integer function read_time(this,itime,nowtime,nexttime)
    implicit none
    class(wrfbdy_tinterp) :: this
    integer, intent(in) :: itime
    type(wrfbdy_var_info), pointer :: varinfo
    integer :: varid_in, ivar
    character(len=19),intent(inout) :: nowtime,nexttime
    character(len=300) :: message

25  format('Read var ',I0,' id ',I0,' type ',I0,' name "',A,'" size ',I0)

    read_time=0
    varloop: do ivar=1,this%nvars
       varinfo=>this%vars(ivar)
       varid_in=this%varids_in(ivar)
       if(this%vars(ivar)%time_dim>0) then
          this%vars(ivar)%start(this%vars(ivar)%time_dim)=itime
       endif
       vartype: select case(varinfo%xtype)
       case(NF90_INT)
          if(.not.associated(varinfo%i4in)) then
             allocate(varinfo%i4in(varinfo%nelem))
             allocate(varinfo%i4out(varinfo%nelem))
          endif
          print 25,ivar,this%varids_in(ivar),varinfo%xtype,varinfo%name,varinfo%nelem
          read_time=this%nfcall('read int data',this%infile,&
               nf90_get_var(this%ncid_in,varid_in,varinfo%i4in,varinfo%start,&
               varinfo%count),varname=varinfo%name)
          if(read_time/=0) return
       case(NF90_FLOAT)
          if(.not.associated(varinfo%r4in)) then
             allocate(varinfo%r4in(varinfo%nelem))
             allocate(varinfo%r4out(varinfo%nelem))
          endif
          print 25,ivar,this%varids_in(ivar),varinfo%xtype,varinfo%name,varinfo%nelem
          read_time=this%nfcall('read float data',this%infile,&
               nf90_get_var(this%ncid_in,varid_in,varinfo%r4in,varinfo%start,&
               varinfo%count),varname=varinfo%name)
          if(read_time/=0) return
       case(NF90_DOUBLE)
          if(.not.associated(varinfo%r8in)) then
             allocate(varinfo%r8in(varinfo%nelem))
             allocate(varinfo%r8out(varinfo%nelem))
          endif 
          print 25,ivar,this%varids_in(ivar),varinfo%xtype,varinfo%name,varinfo%nelem
          read_time=this%nfcall('read double data',this%infile,&
               nf90_get_var(this%ncid_in,varid_in,varinfo%r8in,varinfo%start,&
               varinfo%count),varname=varinfo%name)
          if(read_time/=0) return
       case(NF90_CHAR)
          varname: select case(varinfo%name)
          case('md___thisbdytimee_x_t_d_o_m_a_i_n_m_e_t_a_data_')
             print 25,ivar,this%varids_in(ivar),varinfo%xtype,varinfo%name,varinfo%nelem
             read_time=this%nfcall('read thisbdytime',this%infile,&
                  nf90_get_var(this%ncid_in,varid_in,nowtime,varinfo%start,&
                  varinfo%count),varname=varinfo%name)
             if(read_time/=0) return
          case('md___nextbdytimee_x_t_d_o_m_a_i_n_m_e_t_a_data_')
             print 25,ivar,this%varids_in(ivar),varinfo%xtype,varinfo%name,varinfo%nelem
             read_time=this%nfcall('read nextbdytime',this%infile,&
                  nf90_get_var(this%ncid_in,varid_in,nexttime,varinfo%start,&
                  varinfo%count),varname=varinfo%name)
             if(read_time/=0) return
          end select varname
          cycle varloop
       case default
          write(message,77) varinfo%name
          call this%fail(trim(message))
          read_time=-1
          return
77        format(A,": do not know what to do with this variable's datatype")
       end select vartype
    end do varloop
    this%itime_in=itime
  end function read_time

  ! --------------------------------------------------------------------

  integer function copy_apply_time(this,tendfact,nowtime,nexttime,itime)
    implicit none
    class(wrfbdy_tinterp) :: this
    type(wrfbdy_var_info), pointer :: varinfo
    real, intent(in) :: tendfact
    character(len=19), intent(in) :: nowtime, nexttime
    integer, intent(in) :: itime

    integer :: bdyvar,tendvar, varid_in, varid_out, btype,ttype, idim
    logical :: done(this%nvars)
    integer :: starts(2), counts(2) ! for writing character data
    integer :: ivar,ndims,dims(NF90_MAX_VAR_DIMS),&
         stride(NF90_MAX_VAR_DIMS),dimidx
    character*255 :: message

    copy_apply_time=0

    ! Initialize all variable flags to "I have not copied this one yet."
    done=.false.

!10  format('Read var ',I0,' id ',I0,' type ',I0,' name "',A,'" size ',I0)
10  format('Bdy ',A,' tend ',A,' type ',I0)

    ! First pass: find data for bdy+tend variable pairs and integer variables.
    varloop1: do ivar=1,this%nvars

       varinfo=>this%vars(ivar)
       varid_in=this%varids_in(ivar)
       varid_out=this%varids_out(ivar)

       if(done(ivar)) cycle varloop1

       call this%find_var(ivar,bdyvar,tendvar)
       bdytend: if(bdyvar>=1 .and. tendvar>=1) then
          ! This is a variable with a boundary and tendency
          btype=this%vars(bdyvar)%xtype
          ttype=this%vars(tendvar)%xtype
          whichinterp: select case(ttype)
          case(NF90_FLOAT)
             print 10,this%vars(bdyvar)%name,this%vars(tendvar)%name,ttype
             this%vars(bdyvar)%r4out=this%vars(bdyvar)%r4in+this%vars(tendvar)%r4in*tendfact
             this%vars(tendvar)%r4out=this%vars(tendvar)%r4in
             done(tendvar)=.true.
             done(bdyvar)=.true.
             cycle varloop1
          case(NF90_DOUBLE)
             print 10,this%vars(bdyvar)%name,this%vars(tendvar)%name,ttype
             this%vars(bdyvar)%r8out=this%vars(bdyvar)%r8in+this%vars(tendvar)%r8in*tendfact
             this%vars(tendvar)%r8out=this%vars(tendvar)%r8in
             done(tendvar)=.true.
             done(bdyvar)=.true.
             cycle varloop1
          end select whichinterp
       endif bdytend
       select case(varinfo%xtype)
       case(NF90_CHAR)
          ! Skip these.  We handle them outside the loop.
          cycle varloop1
       case(NF90_INT)
          ! Copy ints without interpolation
          varinfo%i4out=varinfo%i4in
          done(ivar)=.true.
          cycle varloop1
       end select
    end do varloop1

    ! Second pass: find data not yet found, and copy data to new file.
    varloop2: do ivar=1,this%nvars

       varinfo=>this%vars(ivar)
       varid_in=this%varids_in(ivar)
       varid_out=this%varids_out(ivar)

       vardone: if(done(ivar)) then
          ! Set start, counts, ndims for writing variable.
          if(this%vars(ivar)%time_dim>0) then
             this%vars(ivar)%start(this%vars(ivar)%time_dim)=this%itime_out
          endif

          selectwrite: select case(this%vars(ivar)%xtype)
             case(NF90_FLOAT)
                copy_apply_time=this%nfcall('write float data',this%outfile,&
                     nf90_put_var(this%ncid_out,this%varids_out(ivar),&
                     varinfo%r4out,varinfo%start,varinfo%count),&
                     varname=varinfo%name)
                if(copy_apply_time/=0) return
                cycle varloop2
             case(NF90_DOUBLE)
                copy_apply_time=this%nfcall('write double data',this%outfile,&
                     nf90_put_var(this%ncid_out,this%varids_out(ivar),&
                     varinfo%r8out,varinfo%start,varinfo%count),&
                     varname=varinfo%name)
                if(copy_apply_time/=0) return
                cycle varloop2
             case(NF90_INT)
                copy_apply_time=this%nfcall('write int data',this%outfile,&
                     nf90_put_var(this%ncid_out,this%varids_out(ivar),&
                     varinfo%i4out,varinfo%start,varinfo%count),&
                     varname=varinfo%name)
                if(copy_apply_time/=0) return
                cycle varloop2
             end select selectwrite

       end if vardone
       charname: select case(this%vars(ivar)%name)
       case('Times')
          starts=(/1,this%itime_out/)
          counts=(/19,1/)
          copy_apply_time=this%nfcall('write Times',this%outfile,&
               nf90_put_var(this%ncid_out,this%varids_out(ivar),&
               nowtime,starts,counts),varname='Times')
          done(ivar)=.true.
          if(copy_apply_time/=0) return
       case('md___thisbdytimee_x_t_d_o_m_a_i_n_m_e_t_a_data_')
          starts=(/1,this%itime_out/)
          counts=(/19,1/)
          copy_apply_time=this%nfcall('write thisbdytime',this%outfile,&
               nf90_put_var(this%ncid_out,this%varids_out(ivar),&
               nowtime,starts,counts),varname=this%vars(ivar)%name)
          done(ivar)=.true.
          if(copy_apply_time/=0) return
       case('md___nextbdytimee_x_t_d_o_m_a_i_n_m_e_t_a_data_')
          starts=(/1,this%itime_out/)
          counts=(/19,1/)
          copy_apply_time=this%nfcall('write nextbdytime',this%outfile,&
               nf90_put_var(this%ncid_out,this%varids_out(ivar),&
               nexttime,starts,counts),varname=this%vars(ivar)%name)
          done(ivar)=.true.
          if(copy_apply_time/=0) return
       end select charname
       if(.not.done(ivar)) then
38        format(A,': have no way of processing this variable.  Giving up.')
          write(message,38) trim(varinfo%name)
          call this%fail(trim(message))
          copy_apply_time=-1
          return
       endif
    end do varloop2
  end function copy_apply_time

  ! --------------------------------------------------------------------

  integer function copy_dims(this)
    implicit none
    class(wrfbdy_tinterp) :: this
    character(len=NF90_MAX_NAME) :: namebuf
    integer :: dimbuf(this%ndims), ndims, ierr, idim, dimlen
    integer :: include_parents
    character(len=255) :: message

    copy_dims=0
    include_parents=0
    copy_dims=this%nfcall('inquire dimension list',this%infile,nf90_inq_dimids( &
         this%ncid_in,ndims,dimbuf,include_parents))
    if(copy_dims/=0) return

    if(ndims/=this%ndims) then
90     format(A,": ndims=",I0," from nf90_inquire does not match ndims=",I0," from nf90_inq_dimids.")
       write(message,90) this%infile,this%ndims,ndims
       call this%fail(trim(message))
       copy_dims=1
       return
    endif

    this%dims_in=dimbuf(1:ndims)
    !write(0,*) 'copy_dims says file dims: ',ndims,this%dims_in

    do idim=1,ndims
       copy_dims=this%nfcall('inquire dimension info',this%infile,&
          nf90_inquire_dimension(this%ncid_in,dimbuf(idim),namebuf,dimlen),&
          dimid=dimbuf(idim))
       if(copy_dims/=0) return

       this%dimlen(idim)=dimlen

       if(dimbuf(idim)==this%inf_id_in) then
          copy_dims=this%nfcall('create unlimited dimension',this%outfile,&
               nf90_def_dim(this%ncid_out,trim(namebuf),NF90_UNLIMITED,this%dims_out(idim)),&
               dimname=trim(namebuf))
          this%inf_id_out=this%dims_out(idim)
          this%ntimes_in=dimlen
       else
          copy_dims=this%nfcall('create dimension',this%outfile,&
               nf90_def_dim(this%ncid_out,trim(namebuf),dimlen,this%dims_out(idim)),&
               dimname=trim(namebuf))
       endif
       if(copy_dims/=0) return
    enddo
  end function copy_dims

  ! --------------------------------------------------------------------

  integer function create_vars(this)
    implicit none
    class(wrfbdy_tinterp) :: this
    integer ivar, nvars, nvars2, starts(2),counts(2), itime
    character*100 :: message
    character*19 :: timebuf

    ! Get the number of variables.
    create_vars=this%nfcall('get variable count',this%infile,&
         nf90_inquire(this%ncid_in, nvariables = nvars))
    if(create_vars/=0) return
    if(nvars<1) then
30     format(A,': number of variables from nf90_inquire (',I0,') is implausibly low.')
       write(message,30) this%infile,nvars
       call this%fail(trim(message))
       create_vars=-1
       return
    endif
    this%nvars=nvars

    allocate(this%varids_in(nvars),this%varids_out(nvars))
    allocate(this%vars(nvars))
    nvars2=nvars
    create_vars=this%nfcall('get variable list',this%infile,&
         nf90_inq_varids(this%ncid_in,nvars2,this%varids_in))
    if(create_vars/=0) return
    if(nvars/=nvars2) then
40     format(A,': number of variables from nf90_inquire (',I0,') does not match the value from nf90_inq_varids (',I0,')')
       write(message,30) this%infile,nvars,nvars2
       call this%fail(trim(message))
       create_vars=-1
       return
    endif

    do ivar=1,nvars
       create_vars=create_var(this,ivar)
       if(create_vars/=0) return
       create_vars=copy_var_attr(this,this%vars(ivar)%name,&
            this%varids_in(ivar),this%varids_out(ivar),this%vars(ivar)%natts)
       if(create_vars/=0) return
       if(this%vars(ivar)%name == 'Times') then
          ! Need to read times.
          starts=(/1,1/)
          counts=(/19,this%ntimes_in/)
          allocate(this%times_in(this%ntimes_in),this%epochs_in(this%ntimes_in))
          create_vars=this%nfcall('read Times array',this%infile,&
               nf90_get_var(this%ncid_in,this%varids_in(ivar),this%times_in,&
               starts,counts),varname='Times',varid=this%varids_in(ivar))
          if(create_vars/=0) return
       endif
    enddo

    if(associated(this%times_in)) then
       do itime=1,this%ntimes_in
          call c_wrf2epoch(this%times_in(itime),this%epochs_in(itime),create_vars)
          if(create_vars/=0) return
          call c_epoch2wrf(timebuf,this%epochs_in(itime),create_vars)
          if(create_vars/=0) return
          if(timebuf/=this%times_in(itime)) then
55           format('Times(',I3,') = "',A19,'" = ',I0,' converts to "',A19,'" by strptime=>strftime.  Is the punctuation wrong?  Is the date impossible?  Aborting.')
             write(message,55) itime,this%times_in(itime),this%epochs_in(itime),timebuf
             call this%fail(message)
             create_vars=-1
             return
          endif
50        format('Times(',I3,') = "',A19,'" = ',I0,' seconds after 1970 begin.')
          print 50,itime,this%times_in(itime),this%epochs_in(itime)
       enddo
    else
       call this%fail('Cannot find Times array in the input file.')
       create_vars=-1
       return
    endif
  end function create_vars

  ! --------------------------------------------------------------------

  integer function create_var(file,ivar)
    implicit none
    class(wrfbdy_tinterp) :: file
    type(wrfbdy_var_info),pointer :: varinfo
    integer, intent(in) :: ivar
    character(len=NF90_MAX_NAME) :: namebuf
    integer :: dimbuf(NF90_MAX_VAR_DIMS),vardim,filedim,varid_in,nname,ndims
    integer(kind=8) :: nelem
    character(len=500) :: message

    varinfo=>file%vars(ivar)
    varid_in=file%varids_in(ivar)

    create_var=file%nfcall('inquire variable information',file%infile,&
         nf90_inquire_variable(file%ncid_in,varid_in,namebuf,varinfo%xtype, &
         varinfo%ndims,dimbuf,varinfo%natts),varid=varid_in)
    if(create_var/=0) return
    ndims=varinfo%ndims ! to shorten the below code

    ! Store name:
    nname=len_trim(namebuf)
    allocate(character(len=nname) :: varinfo%name)
    varinfo%name(1:nname)=namebuf(1:nname)

    ! Store dim ids from input and output files:
    allocate(varinfo%dims_in(ndims),varinfo%dims_out(ndims))
    allocate(varinfo%start(ndims),varinfo%count(ndims))
    nelem=1
    varinfo%dims_in(1:ndims)=dimbuf(1:ndims)
    !print *,'dims in file:',file%ndims,file%dims_in
    !print *,'dims in var:',size(varinfo%dims_in),varinfo%dims_in
    outer: do vardim=1,varinfo%ndims
       varinfo%start(vardim)=1
       varinfo%count(vardim)=1
       do filedim=1,file%ndims
          if(file%dims_in(filedim)==varinfo%dims_in(vardim)) then
             varinfo%dims_out(vardim)=file%dims_out(filedim)
             if(file%dims_in(filedim) /= file%inf_id_in) then
                nelem = nelem * file%dimlen(filedim)
                varinfo%count(vardim)=file%dimlen(filedim)
             else
                varinfo%time_dim=vardim
             endif
             cycle outer
          endif
       enddo
       ! Should never get here.  Variable has a dimension that was not
       ! defined in the file.
19     format(A,': ',A,': dimension ',I0,' is dimid=',I0, &
            ' in input file, but has no matching dimension from nf90_inquire.')
       write(message,19) file%infile,varinfo%name,vardim,varinfo%dims_in(vardim)
       call file%fail(trim(message))
       create_var=-1
       return
    enddo outer
    varinfo%nelem=nelem

    create_var=file%nfcall('create variable',file%outfile,&
       nf90_def_var(file%ncid_out,varinfo%name,varinfo%xtype,&
       varinfo%dims_out,file%varids_out(ivar)),varname=varinfo%name)
  end function create_var

  ! --------------------------------------------------------------------

  integer function copy_global_attr(this)
    implicit none
    class(wrfbdy_tinterp) :: this
    integer :: iatt
    character(len=NF90_MAX_NAME) :: attname
    logical found_dt
    character*255 :: message
    found_dt=.false.
    copy_global_attr=0
    do iatt=1,this%natts
       copy_global_attr=this%nfcall('inquire global attribute name',this%infile,&
            nf90_inq_attname(this%ncid_in,NF90_GLOBAL,iatt,attname),&
            attnum=iatt)
       if(copy_global_attr/=0) return

       if(trim(attname)=='SIMULATION_START_DATE' .or. &
            trim(attname)=='START_DATE') then
          ! Special case.  Replace with new start date.
          copy_global_attr=this%nfcall('write new start date/time',this%outfile,&
               nf90_put_att(this%ncid_out,NF90_GLOBAL,trim(attname),this%new_time0),&
               attname=trim(attname))
          if(copy_global_attr/=0) return
       else
          ! Not so special.  Copy attribute.
          copy_global_attr=this%nfcall('copy global attribute',this%outfile,&
               nf90_copy_att(this%ncid_in,NF90_GLOBAL,trim(attname),this%ncid_out,NF90_GLOBAL),&
               attname=trim(attname))
          if(copy_global_attr/=0) return
       endif

       if(trim(attname)=='DT') then
          ! Need value of DT for time stepping
          found_dt=.true.
          this%dt=-1.
          copy_global_attr=this%nfcall('read DT',this%infile,&
               nf90_get_att(this%ncid_in,NF90_GLOBAL,'DT',this%dt),&
               attname='DT')
          if(this%dt<1.0) then
20           format(A,': time step must be 1 second or greater (it is ',F0.3,').')
             write(message,20) this%infile,this%dt
             call this%fail(trim(message))
             copy_global_attr=-1
             return
          endif
       endif
    end do
    if(.not.found_dt) then
10     format(A,': could not find global attribute "DT".')
       write(message,10) this%infile
       call this%fail(trim(message))
       copy_global_attr=-1
       return
    endif
  end function copy_global_attr

  ! --------------------------------------------------------------------
  
  integer function copy_var_attr(this,varname,varid_in,varid_out,natts)
    implicit none
    class(wrfbdy_tinterp) :: this
    integer,intent(in) :: varid_in,varid_out,natts
    character(len=*), intent(in) :: varname
    integer :: iatt
    character(len=NF90_MAX_NAME) :: attname

    copy_var_attr=0
    do iatt=1,natts
       copy_var_attr=this%nfcall('inquire var attribute name',this%infile,&
            nf90_inq_attname(this%ncid_in,varid_in,iatt,attname),&
            attnum=iatt,varname=varname,varid=varid_in)
       if(copy_var_attr/=0) return
       
       copy_var_attr=this%nfcall('copy var attribute',this%outfile,&
            nf90_copy_att(this%ncid_in,varid_in,trim(attname),this%ncid_out,varid_out),&
            varname=varname,attname=trim(attname),attnum=iatt)
       if(copy_var_attr/=0) return
    end do

  end function copy_var_attr

end module module_wrfbdy_tinterp

program prog_wrfbdy_tinterp
  use module_wrfbdy_tinterp
  implicit none
  type(wrfbdy_tinterp) :: interper
  character(len=19) :: times(2),timedup
  character(len=120) :: charbuf
  integer(kind=8) :: epoch(2), bdystep
  integer :: ierr,narg,nout,arglen, ntimes_out,iarg

  narg=command_argument_count()
  if(narg/=3) call usage('wrong number of arguments')

  ! Parse start time:
  do iarg=1,2
     call get_command_argument(iarg,length=arglen,status=ierr)
     if(ierr/=0) call iusage(iarg,'absent or empty')
     if(arglen/=19) call iusage (iarg,'must be 19 characters')
     call get_command_argument(iarg,times(iarg),status=ierr)
     if(ierr/=0) call iusage(iarg,'absent or empty')
     call c_wrf2epoch(times(iarg),epoch(iarg),ierr)
     if(ierr/=0) call iusage(iarg,'is unparseable')
     timedup=' '
     call c_epoch2wrf(timedup,epoch(iarg),ierr)
     if(ierr/=0) call iusage(iarg,'cannot convert date+time unix epoch back to WRF time')
     if(times(iarg)/=timedup) then
10      format('Date+time input "',A,'" is the wrong format; should be "',A,'"')
        write(charbuf,10) times(iarg),timedup
     endif
  enddo

  ! Parse bdystep
  charbuf=' '
  call get_command_argument(3,charbuf,status=ierr)
  if(ierr/=0) call usage('error reading third argument')
  read(charbuf,'(I120)',iostat=ierr) bdystep
  if(ierr/=0) call usage('error parsing third argument (bdystep)')
  if(bdystep<1) call usage('bdystep must be 1 or larger')

  ! Get ntimes_out
  ntimes_out=(epoch(2)-epoch(1)+bdystep/10)/bdystep
  if(ntimes_out<1) call usage('end time must be at least one bdystep after start time')
  if(abs(epoch(1)+ntimes_out*bdystep - epoch(2))>bdystep/10) then
     call usage('end time is not on or near a boundary step')
  endif

20 format('Interp to ',A19,' => ',A19,' step ',I0,' seconds total ',I0,' times.')
  print 20,times(1),times(2),bdystep,ntimes_out
  
  ! Open files and define structure:
  ierr=interper%initial_io('./wrfbdy_d01.in','./wrfbdy_d01.out',times(1))
  if(ierr/=0) then
     write(0,'(A)') 'Error in initial I/O.  See above for details.'
     stop 7
  endif

  ! Interpolate boundary times from source to target:
  ierr=interper%interp(times(1),bdystep,ntimes_out)
  if(ierr/=0) then
     write(0,'(A)') 'Error during boundary interpolation.  See above for details.'
     stop 7
  endif

  ! Close output and input files, flushing contents to disk:
  ierr=interper%close()
  if(ierr/=0) then
     write(0,'(A)') 'Error closing files.  See above for details.'
     stop 7
  endif

contains !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine iusage(iarg,why)
    integer, intent(in) :: iarg
    character(len=*), intent(in) :: why
    character(len=len_trim(why)+50) :: buf
10  format('error in argument ',I0,': ',A)
    write(buf,10) iarg,trim(why)
    call usage(trim(buf))
  end subroutine iusage

  subroutine usage(why)
    character(len=*), intent(in), optional :: why

    write(0,'(A)') 'wrfbdy_tinterp version 0.1'
    write(0,'(A)') 'Syntax:'
    write(0,'(A)') '    wrfbdy_tinterp start-time end-time bdystep'
    write(0,'(A)') 'Synopsis: '
    write(0,'(A)') '    Interpolates NetCDF WRF-NMM boundary data temporally from file'
    write(0,'(A)') '    ./wrfbdy_d01.in to ./wrfbdy_d01.out'
    write(0,'(A)') 'Arguments:'
    write(0,'(A)') '    start-time = "2018-12-31 17:18:19"'
    write(0,'(A)') '      end-time = "2018-12-31-23:18:19"'
    write(0,'(A)') '   bdystep (s) = "3600"'
    if(present(why)) then
       write(0,"('ABORTING DUE TO INCORRECT ARGUMENTS: ',A)") why
       stop 2
    else
       stop 1
    endif
end subroutine usage
end program prog_wrfbdy_tinterp
