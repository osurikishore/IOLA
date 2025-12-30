module ens_prob_module
  use fileop_module, only: stat
  use datetime_module, only: grads_date, grads_date_len, epochnow
  use projection_module, only: proj_latlon, init_proj_latlon, projection
  use decomp_module, only: decomp, init_decomp
  use wrfdiag_module, only: wrfdiag_file, init_wrfdiag_file
  use vardata_module, only: vardata_real, init_vardata_real, &
       full_init_vardata_real, vardata_int, full_init_vardata_int
  use wrfdiag_module, only: init_wrfdiag_file, wrfdiag_var, &
       wrfdiag_file, init_wrfdiag_var
  use interp_module, only: bilinear_real, pd_init_bilinear_real, &
       OP_MAX, OP_MIN
  use tcf_module, only: atcf,read_atcf
  use distance_module, only: greatarc

  implicit none

  integer, public, parameter :: maxensmax=99, minensmin=0
  integer, public, parameter :: maxdommax=21, mindommin=1
  integer, public, parameter :: max_output_times=9999
  integer, public, parameter :: max_thresholds=9999

  integer, public, parameter :: filenamelen = 2000
  integer, public, parameter :: varnamelen = 40
  real, public, parameter :: dgepsilon=1e-5, dgepsilon2=dgepsilon*dgepsilon
  real, public, parameter :: mepsilon=11
  real, public, parameter :: grads_badval=-999.0
  integer, public, parameter :: grads_first_unit=1000

  public :: outputter, member, ensemble
  public :: init_outputter, init_member, init_ensemble
  public :: insert_ens_domain, dist_weight

  type outputter
     real, pointer :: output_times(:) ! forecast hours at which to output
     integer :: notimes ! number of output times (length of output_times)
     real :: tg_freq ! how long (hours) to accum. and output tornado vars
     real, pointer :: wthresh(:),rthresh(:)
     integer :: nwthresh=0,nrthresh=0
   contains
     procedure free => free_outputter

     ! Time loop:
     procedure start => start_nothing         ! called before time loop
     procedure step_time => step_time_outputter ! called once per out time
     procedure finish => do_nothing       ! called after time loop

     ! Data types:
     procedure trackprob => output_nothing
     procedure rainswath => output_nothing
     procedure windprob => output_nothing_thresh
     procedure rainprob => output_nothing_thresh
     procedure udhelfrac => output_nothing_thresh
     procedure udhel => output_nothing_levels
     procedure maxw => output_nothing_tornado
     procedure minw => output_nothing_tornado
     procedure maxu10 => output_nothing_tornado
     procedure member_rain => output_nothing_member
     procedure member_wind => output_nothing_member
     procedure member_prob => output_nothing_member
  end type outputter

  type member
     class(proj_latlon), pointer :: proj
     class(decomp), pointer :: dec
     type(vardata_real) :: wind, rain, prob, udh03, udh25, minw, maxw, maxu10
     type(wrfdiag_file), pointer :: files(:) => NULL()
     character(len=filenamelen) :: indiag_format, inatcf_format
     real :: lasthour=-9e9
     integer :: mindom=-999,maxdom=-999,trackdom=-999,lastitime=-999,ensid=-999
     real(kind=8) :: atime=-9e9
     logical :: local=.false., valid=.false.
     real, pointer :: lats(:),lons(:),track_lat(:),track_lon(:),ftimes(:)
     real :: clon=0.0
     integer :: ntimes=0, nx=0, ny=0
   contains
     !procedure init => init_member
     procedure check_input => check_input_member
     procedure open => open_member
     procedure read_tracks => read_tracks_member
     procedure center_longitude => center_longitude_member
     procedure bound => bound_member
     procedure init_arrays => init_arrays_member
     procedure track_probability => track_probability_member
     procedure accum_prob
     procedure loop_prep => loop_prep_member
     procedure open_var => open_var_member
     procedure step_time => step_time_member
     procedure set_grid
  end type member

  type ensemble
     logical, pointer :: ensflag(:)=>NULL(), sent(:)=>NULL()
     real, pointer :: ftimes(:)
     integer :: enscount=0,lastitime=-999
     type(proj_latlon) :: proj
     type(decomp) :: dec
     type(member), pointer :: members(:)=>NULL()
     type(vardata_real) :: rainaccum, rainall
     type(vardata_real) :: hold, extrema, prob
     type(vardata_int) :: raincount
     real(kind=8) :: atime=-9e9, clon=0.0
     integer :: ensmin=-999,ensmax=-999,mindom=-999,maxdom=-999,trackdom=-999
     integer :: ntimes=0, nx=0,ny=0, most_probable_member=-999
     real :: lat1=0., lon1=0., latres=0., lonres=0., last_tg_reset=-9e9
     real, pointer :: lats(:), lons(:)
     logical :: reinit=.false.
   contains
     procedure init => init_ensemble
     procedure check_input => check_input_ensemble
     procedure open => open_ensemble
     procedure read_tracks => read_tracks_ensemble
     procedure center_longitude => center_longitude_ensemble
     procedure bound => bound_ensemble
     procedure init_arrays => init_arrays_ensemble
     procedure track_probability => track_probability_ensemble
     procedure combine_track_prob
     procedure most_probable_track
     !procedure tg_reset
     procedure loop_prep => loop_prep_ensemble
     procedure step_time => step_time_ensemble

     !procedure output_track
     procedure output_swath
     procedure output_tornado

     procedure run

  end type ensemble

  ! --------------------------------------------------------------------

contains

  ! --------------------------------------------------------------------

  subroutine init_outputter(this,otimes,notimes,wthresh,&
       nwthresh,rthresh,nrthresh,tg_freq)
    class(outputter) :: this
    real, intent(in) :: otimes(notimes),rthresh(nrthresh),wthresh(nwthresh)
    integer, intent(in) :: notimes,nrthresh,nwthresh
    real, intent(in) :: tg_freq
    ! Local
    integer :: itime
30  format('init outputter with notimes=',I0,' nrthresh=',I0,' nwthresh=',I0)
    print 30,notimes,nrthresh,nwthresh
    allocate(this%output_times(notimes),this%rthresh(nrthresh),&
         this%wthresh(nwthresh))
    this%output_times=otimes
    this%wthresh=wthresh
    this%rthresh=rthresh
    this%tg_freq=tg_freq
    this%nrthresh=nrthresh
    this%nwthresh=nwthresh
    this%notimes=notimes
  end subroutine init_outputter

  subroutine free_outputter(this)
    class(outputter) :: this
    if(associated(this%output_times)) then
       deallocate(this%output_times)
       nullify(this%output_times)
    endif
    if(associated(this%wthresh)) then
       deallocate(this%wthresh)
       nullify(this%wthresh)
    endif
    if(associated(this%rthresh)) then
       deallocate(this%rthresh)
       nullify(this%rthresh)
    endif
  end subroutine free_outputter

  subroutine do_nothing(this)
    class(outputter) :: this
  end subroutine do_nothing

  subroutine step_time_outputter(this,itime,rtime)
    class(outputter) :: this
    integer, intent(in) :: itime
    real, intent(in) :: rtime
  end subroutine step_time_outputter

  subroutine start_nothing(this,lats,lons,atime,nx,ny,fcstlen)
    class(outputter) :: this
    integer, intent(in) :: nx, ny
    real, intent(in) :: lats(ny), lons(nx)
    real(kind=8), intent(in) :: atime, fcstlen
  end subroutine start_nothing

  subroutine output_nothing(this,data,itime)
    class(outputter) :: this
    class(vardata_real) :: data
    integer, intent(in) :: itime
  end subroutine output_nothing

  subroutine output_nothing_tornado(this,data,itgtime)
    class(outputter) :: this
    class(vardata_real) :: data
    integer, intent(in) :: itgtime
  end subroutine output_nothing_tornado

  subroutine output_nothing_thresh(this,data,itime,ithresh,rthresh)
    class(outputter) :: this
    class(vardata_real) :: data
    integer, intent(in) :: itime
    integer, intent(in) :: ithresh
    real, intent(in) :: rthresh
  end subroutine output_nothing_thresh

  subroutine output_nothing_levels(this,data,itgtime,lev1,lev2)
    class(outputter) :: this
    class(vardata_real) :: data
    integer, intent(in) :: itgtime
    real, intent(in) :: lev1,lev2
  end subroutine output_nothing_levels

  subroutine output_nothing_member(this,data,itime,imember)
    class(outputter) :: this
    class(vardata_real) :: data
    integer, intent(in) :: itime
    integer, intent(in) :: imember
  end subroutine output_nothing_member

  subroutine output_nothing_thresh_levels(this,data,itime,ithresh,rthresh,lev1,lev2)
    class(outputter) :: this
    class(vardata_real) :: data
    integer, intent(in) :: itime
    integer, intent(in) :: ithresh
    real, intent(in) :: rthresh,lev1,lev2
  end subroutine output_nothing_thresh_levels

  ! --------------------------------------------------------------------

  subroutine init_member(this,ensid,indiag_format,inatcf_format,mindom,maxdom,trackdom)
    class(member) :: this
    character(len=filenamelen),intent(in) :: indiag_format,inatcf_format
    integer,intent(in) :: mindom,maxdom,ensid,trackdom

    allocate(this%files(mindom:maxdom))
    this%local=.false.
    this%valid=.false.
    this%trackdom=trackdom
    this%atime=-9e9
    this%lasthour=-9e9
    this%ensid=ensid
    this%mindom=mindom
    this%maxdom=maxdom
    this%indiag_format=indiag_format
    this%inatcf_format=inatcf_format
  end subroutine init_member

  subroutine init_ensemble(this, indiag_format, inatcf_format, ensmin,ensmax, &
       mindom,maxdom, trackdom, latres,lonres)
    class(ensemble) :: this
    integer, intent(in) :: ensmin,ensmax, mindom,maxdom, trackdom
    character(len=filenamelen) :: indiag_format, inatcf_format
    real, intent(in) :: latres,lonres
    integer :: ensid
    this%latres=latres
    this%lonres=lonres
    this%ensmin=ensmin
    this%ensmax=ensmax
    this%mindom=mindom
    this%maxdom=maxdom
    this%trackdom=trackdom
    allocate(this%ensflag(ensmin:ensmax))
    allocate(this%members(ensmin:ensmax))
    do ensid=ensmin,ensmax
       call init_member(this%members(ensid), ensid, indiag_format, &
            inatcf_format, mindom,maxdom, trackdom)
    enddo
  end subroutine init_ensemble

  ! --------------------------------------------------------------------

  logical function check_input_member(this,min_wrfdiag_sizes,minage)
    ! Determines if all inputs for the specified ensemble member are present.
    class(member) :: this
    integer(kind=8) :: min_wrfdiag_sizes(this%maxdom-this%mindom+1)
    real(kind=8), intent(in) :: minage
    ! locals
    type(stat) :: s
    integer(kind=8) :: minsize
    real(kind=8) :: age
    logical :: check_ens
    integer :: domid
    character(len=filenamelen) :: indiag,inatcf
10  format('ENS ',I0,A)
20  format('ENS ',I0,' ATCF file ',A,': ',A)
25  format('ENS ',I0,' ATCF file is too new (',F0.3,'<',F0.3,'): ',A)
30  format('ENS ',I0,' domain ',I0,' wrfdiag file ',A,': ',A)
35  format('ENS ',I0,' domain ',I0,' wrfdiag file  is too new (',F0.3,&
         '<',F0.3,'): ',A)
40  format('ENS ',I0,' domain ',I0,' wrfdiag file is too small (',I0,'<',I0,&
         ' bytes): ',A)

    print 10,this%ensid,": Do I have this member's data?"

    check_ens=.false.
    this%valid=.false.

    call insert_ens_domain(this%inatcf_format,inatcf,this%ensid,0)
    if( .not.s%stat(trim(inatcf))) then
       print 20,this%ensid,'does not exist',trim(inatcf)
       goto 100 ! error handler
    elseif(.not.s%isreg) then
       print 20,this%ensid,'is not a file',trim(inatcf)
       goto 100 ! error handler
    else
       age=epochnow()-s%mtime
       if(age < minage) then
          print 25,this%ensid,age,minage,trim(inatcf)
          goto 100 ! error handler
       else
          print 20,this%ensid,'is present',trim(inatcf)
       endif
    endif

    do domid=this%mindom,this%maxdom
       call insert_ens_domain(this%indiag_format,indiag,this%ensid,domid)
       minsize=min_wrfdiag_sizes(domid-this%mindom+1)
       if(.not.s%stat(trim(indiag))) then
          print 30,this%ensid,domid,'does not exist',trim(indiag)
          goto 100 ! error handler
       elseif(.not.s%isreg) then
          print 30,this%ensid,domid,'is not a regular file',trim(indiag)
          goto 100 ! error handler
       elseif(s%size < minsize) then
          print 40,this%ensid,domid,s%size,minsize,trim(indiag)
          goto 100 ! error handler
       else
          age=epochnow()-s%mtime
          if(age < minage) then
             print 35,this%ensid,domid,age,minage,trim(indiag)
             goto 100 ! error handler
         else
             print 30,this%ensid,domid,'is present, large enough, and old enough',trim(indiag)
          endif
       endif
    enddo

    this%valid=.true.
    check_input_member=.true.
    print 10,this%ensid,": Member is usable."
    return

100 continue ! Error handling section (why no try...catch in fortran?)
    print 10,this%ensid,": Member is not usable."
    this%valid=.false.
    check_input_member=.false.
    return
  end function check_input_member

  logical function check_input_ensemble(this,min_wrfdiag_sizes,minage)
    class(ensemble) :: this
    integer(kind=8) :: min_wrfdiag_sizes(this%maxdom-this%mindom+1)
    real(kind=8), intent(in) :: minage
    ! locals
    integer :: ensid

    this%ensflag=.false.
    this%enscount=0
10  format(A)
    print 10,'Check ensemble input.'

    do ensid=this%ensmin,this%ensmax
       this%ensflag(ensid) = this%members(ensid)%&
            check_input(min_wrfdiag_sizes,minage)
       if(this%ensflag(ensid)) this%enscount=this%enscount+1
    end do
    if(this%enscount>0) then
       write(0,222) this%enscount,this%ensmax-this%ensmin+1
222    format('Have ',I0,' of ',I0,' members.')
       check_input_ensemble=.true.
    else
       write(0,220) this%ensmax-this%ensmin+1
       write(6,220) this%ensmax-this%ensmin+1
       check_input_ensemble=.false.
220    format('Have NONE of ',I0,' members.  Aborting.')
    endif
  end function check_input_ensemble

  logical function have_members(this)
    class(ensemble) :: this
    have_members=this%enscount>0
  end function have_members

  ! --------------------------------------------------------------------

  subroutine open_member(this)
    class(member) :: this
    character(len=filenamelen) :: indiag
    integer :: domid,maxntimes
    logical :: have_atime

    have_atime=.false.
    maxntimes=0
    do domid=this%mindom,this%maxdom
       call insert_ens_domain(this%indiag_format,indiag,this%ensid,domid)
33     format('ENS ',I0,': open file "',A,'"')
       print 33,this%ensid,trim(indiag)
       call init_wrfdiag_file(this%files(domid),indiag)
       if(have_atime) then
          if(abs(this%atime-this%files(domid)%epochtimes(1))>1) then
41           format(A)
             print 41,'ERROR: Analysis times do not match between various wrfdiag files.'
          endif
       else
          this%atime=this%files(domid)%epochtimes(1)
          have_atime=.true.
       endif
       maxntimes=max(maxntimes,this%files(domid)%ntimes)
22     format('maxntimes is now ',I0,' from dom ',I0,' ntimes=',I0)
       !print 22,maxntimes,domid,this%files(domid)%ntimes
    end do
    this%ntimes=maxntimes
  end subroutine open_member

  subroutine open_ensemble(this)
    class(ensemble) :: this
    integer :: ensid,maxntimes
    logical :: have_atime
    real(kind=8) :: atime
    have_atime=.false.
    maxntimes=0
    do ensid=this%ensmin,this%ensmax
10     format('Open ens ',I0,'?')
       !print 10,ensid
       if(.not.this%ensflag(ensid)) cycle
20     format('Call this%members(',I0,')%open()')
       !print 20,ensid
       call this%members(ensid)%open()
       if(have_atime) then
          if(abs(this%members(ensid)%atime-atime)>1) then
41           format(A)
             print 41,'ERROR: Analysis times do not match between various wrfdiag files.'
          endif
       else
          atime=this%members(ensid)%atime
       endif
       maxntimes=max(this%members(ensid)%ntimes,maxntimes)
    enddo
    this%atime=atime
    this%ntimes=maxntimes
  end subroutine open_ensemble

  ! --------------------------------------------------------------------
  
  subroutine read_tracks_member(this)
    class(member) :: this
    ! Determines the last forecast hour seen in an ATCF file.
    type(atcf), target :: datcf(500)
    class(wrfdiag_file), pointer :: file
    character(len=filenamelen) :: inatcf
    integer :: used, itime, xmid,ymid
    real :: lasthour, tlat,tlon
    used=0
    lasthour=0
403 format('Read ATCFUNIX file "',A,'"...')

    ! FIRST: read the atcf track
    call insert_ens_domain(this%inatcf_format,inatcf,this%ensid,0)
    print 403,trim(inatcf)
    call read_atcf(inatcf,datcf,500,used)
    if(used>0) then
       lasthour=datcf(used)%fcsthour
404    format('Last hour was ',F0.1)
       print 404,lasthour
    else
405    format(A)
       print 405,'No data read in from file!!!'
       print 405,'Will assume last hour is 0.'
    endif
    this%lasthour=lasthour

    ! SECOND: read the 126hr track
    file=>this%files(this%trackdom)
    allocate(this%track_lat(this%ntimes))
    allocate(this%track_lon(this%ntimes))
    allocate(this%ftimes(this%ntimes))
    do itime=1,this%ntimes
       tlat=file%track_lat(itime)
       tlon=file%track_lon(itime)
       if(tlat<-998. .and. tlon<-998.) then
          call file%make_latlonrot(itime)
          xmid=ceiling(file%de_ij%ide/2.0)
          ymid=ceiling(file%de_ij%jde/2.0)
          tlon=file%lons%rdata(xmid,ymid,1)
          tlat=file%lats%rdata(xmid,ymid,1)
          call file%free_latlonrot()
          file%track_lat(itime)=tlat
          file%track_lon(itime)=tlon
38        format('ENS ',I0,' itime=',I0,' no lat/lon.  Use dom center lat=',F0.3,' lon=',F0.3)
          print 38,this%ensid,itime,tlat,tlon
       else
          !print 39,this%ensid,itime,tlat,tlon
39        format('ENS ',I0,' itime=',I0,' lat/lon from internal track lat=',F0.3,' lon=',F0.3)
       end if
       this%track_lat(itime)=tlat
       this%track_lon(itime)=tlon
       this%ftimes(itime)=file%epochtimes(itime)-file%epochtimes(1)
    enddo
  end subroutine read_tracks_member

  subroutine center_longitude_member(this,clon)
    ! Ensures all longitudes are within +/- 180 degrees of clon.
    ! Internally records clon for later calculations, to ensure those
    ! are within the same range.
    class(member) :: this
    real, intent(in) :: clon
    integer :: itime
    this%clon=clon
    do itime=1,this%ntimes
       this%track_lon(itime)=mod(3600.0+180.0+this%track_lon(itime)-clon,&
            360.0)-180.0+clon
    end do
  end subroutine center_longitude_member

  subroutine read_tracks_ensemble(this)
    class(ensemble) :: this
    integer :: ensid
    logical :: have_clon
    real :: clon
    real(kind=8) :: lat_accum(this%ntimes), lon_accum(this%ntimes)
    integer :: count
    have_clon=.false.
    do ensid=this%ensmin,this%ensmax
       if(.not.this%ensflag(ensid)) cycle
       call this%members(ensid)%read_tracks()
       if(.not.have_clon) then
          clon=this%members(ensid)%track_lon(1)
          this%clon=clon
          have_clon=.true.
          allocate(this%ftimes(this%ntimes))
          this%ftimes=this%members(ensid)%ftimes
       endif
    enddo

    call this%center_longitude(clon)
  end subroutine read_tracks_ensemble

  subroutine center_longitude_ensemble(this,clon)
    ! Ensures all longitudes are within +/- 180 degrees of clon.
    ! Internally records clon for later calculations, to ensure those
    ! are within the same range.
    class(ensemble) :: this
    real, intent(in) :: clon
    integer ensid
    this%clon=clon
    do ensid=this%ensmin,this%ensmax
       if(.not.this%ensflag(ensid)) cycle
       call this%members(ensid)%center_longitude(clon)
    end do
  end subroutine center_longitude_ensemble

  ! --------------------------------------------------------------------

  subroutine bound_member(this,blon1,blat1,blon2,blat2,blonmid,blatmid,midcount)
    class(member) :: this
    real, intent(inout) :: blon1,blat1,blon2,blat2,blonmid,blatmid
    integer, intent(inout) :: midcount
    real :: elat1,elatmid,elat2, elon1,elonmid,elon2 ! ensemble member grid
    real :: clon
    integer :: domid

    clon=this%clon
    do domid=this%mindom,this%maxdom
       call this%files(domid)%latlon_bounds(elat1,elatmid,elat2, &
            elon1,elonmid,elon2)
       if (clon<-998.) then
          clon=elonmid
       endif
       elon1=mod(3600.0+180.0+elon1-clon,360.0)-180.0+clon
       elon2=mod(3600.0+180.0+elon2-clon,360.0)-180.0+clon
       elonmid=mod(3600.0+180.0+elonmid-clon,360.0)-180.0+clon
       
       !filecount=filecount+1
       blon1=min(blon1,elon1)
       blat1=min(blat1,elat1)
       blon2=max(blon2,elon2)
       blat2=max(blat2,elat2)
       blonmid=blonmid+elonmid
       blatmid=blatmid+elatmid
       midcount=midcount+1

5      format('ENS ',I0,' DOM ',I0)
10     format('Bound lat=[',F0.3,',',F0.3,'] lon=[',F0.3,',',F0.3,']')
20     format('Accumulated mid lat=',F0.3,' lon=',F0.3,' cenlon=',F0.3)
       !print 5,this%ensid,domid
       !print 10,blat1,blat2,blon1,blon2
       !print 20,blatmid,blonmid,clon
    enddo

  end subroutine bound_member

  subroutine bound_ensemble(this)
    class(ensemble) :: this
    real :: blon1,blat1,blon2,blat2,blonmid,blatmid
    real :: elat1,elatmid,elat2, elon1,elonmid,elon2 ! ensemble member grid
    real :: clon
    integer :: ensid,i,j,midcount

    blon1=9e9
    blat1=9e9
    blon2=-9e9
    blat2=-9e9
    blonmid=0.
    blatmid=0.
    clon=this%clon
    midcount=0

    if(this%enscount<1) then
       write(0,*) 'ERROR: enscount<1.  No members.  Abort.'
       write(6,*) 'ERROR: enscount<1.  No members.  Abort.'
    endif

    do ensid=this%ensmin,this%ensmax
       !print *,'Check bound for ens ',ensid,'???'
       if(.not.this%ensflag(ensid)) cycle
       !print *,'Check bound for ens ',ensid,'...'
       call this%members(ensid)%bound(blon1,blat1,blon2,blat2,blonmid,&
            blatmid,midcount)
    enddo

    !print *,'lon diff',blon2-blon1
    !print *,'lat diff',blat2-blat1
    !print *,'res lat,lon',this%latres,this%lonres
    this%nx=ceiling( (blon2-blon1)/this%lonres )
    this%ny=ceiling( (blat2-blat1)/this%latres )
    blonmid=(blon2+blon1)/2.
    blatmid=(blat2+blat1)/2.
    this%lat1=blatmid-this%ny/2.0*this%latres
    this%lon1=blonmid-this%nx/2.0*this%lonres
    
    allocate(this%lats(this%ny))
    allocate(this%lons(this%nx))

    do i=1,this%nx
       this%lons(i)=this%lon1+i*this%lonres
    enddo
    do j=1,this%ny
       this%lats(j)=this%lat1+j*this%latres
    enddo

    call init_proj_latlon(this%proj,this%lat1,this%lon1,this%latres,this%lonres,this%nx,this%ny,.true.)
    call init_decomp(this%dec, 1,this%nx,1,this%ny,1,1, 1,this%nx,1,this%ny,1,1)
22  format('Call set_grid for ensid=',I0,'...',I0,' nx=',I0,' ny=',I0)
    !print 22,this%ensmin,this%ensmax,this%nx,this%ny
    do ensid=this%ensmin,this%ensmax
       if(.not.this%ensflag(ensid)) cycle
       call this%members(ensid)%set_grid(this%proj,this%dec,&
            this%nx,this%ny,this%lats,this%lons)
    end do

10  format('Bound lat=[',F0.3,',',F0.3,'] lon=[',F0.3,',',F0.3,']')
20  format('Mid lat=',F0.3,' lon=',F0.3,' cenlon=',F0.3,' nx=',I0,' ny=',I0)
30  format('Grid lat1=',F0.3,' lon1=',F0.3,' latres=',F0.5,' lonres=',F0.5)
    print 10,blat1,blat2,blon1,blon2
    print 20,blatmid,blonmid,clon,this%nx,this%ny
    print 30,this%lat1,this%lon1,this%latres,this%lonres

  end subroutine bound_ensemble

  subroutine set_grid(this,proj,dec,nx,ny,lats,lons)
    class(member) :: this
    class(proj_latlon), target :: proj
    class(decomp), target :: dec
    integer, intent(in) :: nx,ny
    real,intent(in) :: lats(ny), lons(nx)

    allocate(this%lats(ny), this%lons(nx))

    this%proj=>proj
    this%dec=>dec
    this%nx=nx
    this%ny=ny

    this%lons=lons
    this%lats=lats

23  format('ENS ',I0,' nx=',I0,' ny=',I0)
    !print 23,this%ensid,this%nx,this%ny
  end subroutine set_grid

  ! --------------------------------------------------------------------

  subroutine init_iarray(iarray,val,proj,dec)
    type(vardata_int),intent(inout) :: iarray
    class(projection) :: proj
    class(decomp) :: dec
    integer, intent(in) :: val
    integer :: i,j

    call full_init_vardata_int(iarray,proj,dec,.true.,.true.,.false.,.false.)
    !$OMP PARALLEL DO PRIVATE(i,j)
    do j=1,dec%jde
       do i=1,dec%ide
          iarray%mask(i,j,1)=.false.
          iarray%idata(i,j,1)=val
       enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine init_iarray

  subroutine init_rarray(rarray,val,proj,dec)
    class(vardata_real),intent(inout) :: rarray
    class(projection) :: proj
    class(decomp) :: dec
    real, intent(in) :: val
    integer :: i,j

    call full_init_vardata_real(rarray,proj,dec,.true.,.true.,.false.,.false.)
    !$OMP PARALLEL DO PRIVATE(i,j)
    do j=1,dec%jde
       do i=1,dec%ide
          rarray%mask(i,j,1)=.false.
          rarray%rdata(i,j,1)=val
       enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine init_rarray

  subroutine init_arrays_member(this)
    class(member) :: this
    call init_rarray(this%wind,0.0,this%proj,this%dec)
    call init_rarray(this%rain,0.0,this%proj,this%dec)
    call init_rarray(this%prob,0.0,this%proj,this%dec)

    call init_rarray(this%udh03,0.0,this%proj,this%dec)
    call init_rarray(this%udh25,0.0,this%proj,this%dec)
    call init_rarray(this%minw,0.0,this%proj,this%dec)
    call init_rarray(this%maxw,0.0,this%proj,this%dec)
    call init_rarray(this%maxu10,0.0,this%proj,this%dec)
  end subroutine init_arrays_member

  subroutine init_arrays_ensemble(this)
    class(ensemble) :: this
    integer :: ensid
    call init_rarray(this%rainall,0.0,this%proj,this%dec)
    call init_rarray(this%rainaccum,0.0,this%proj,this%dec)
    call init_iarray(this%raincount,0,this%proj,this%dec)
    call init_rarray(this%hold,0.0,this%proj,this%dec)
    call init_rarray(this%prob,0.0,this%proj,this%dec)
    call init_rarray(this%extrema,0.0,this%proj,this%dec)

    ! Initialize arrays in each ensemble member:
    do ensid=this%ensmin,this%ensmax
       if(.not.this%ensflag(ensid)) cycle
       call this%members(ensid)%init_arrays()
    end do
  end subroutine init_arrays_ensemble

  ! --------------------------------------------------------------------

  subroutine track_probability_member(this, itime, reset, last, dt2, &
       stdev0, sd_per_hour, maxdist, cendist, minprob)
    class(member) :: this
    logical, intent(in) :: reset
    integer, intent(in) :: itime
    real, intent(in) :: stdev0, sd_per_hour, maxdist, cendist, minprob
    real, intent(inout) :: dt2
    logical, intent(inout) :: last
    !locals
    integer :: i,j,nmask
    real :: tlat1,tlon1,tlat2,tlon2,dt1, stdev1,stdev2
    if(reset) then
       !$OMP PARALLEL DO PRIVATE(i,j)
       do j=1,this%ny
          do i=1,this%nx
             this%prob%mask(i,j,1)=.false.
             this%prob%rdata(i,j,1)=.false.
          enddo
       enddo
       !$OMP END PARALLEL DO
    endif

    tlat1=this%track_lat(itime)
    tlon1=this%track_lon(itime)
    tlat2=this%track_lat(itime+1)
    tlon2=this%track_lon(itime+1)
    dt1=this%ftimes(itime)/3600.0
    dt2=this%ftimes(itime+1)/3600.
    stdev1 = stdev0 + dt1*sd_per_hour
    stdev2 = stdev0 + dt2*sd_per_hour

37  format('track prob ENS ',I0,' cendist=',F0.3)
38  format('At hour ',F0.1,'-',F0.1', stdev1=',F0.1,' and stdev2=',F0.1)
39  format('   lat1=',F0.2,' lon1=',F0.2,' lat2=',F0.2,' lon2=',F0.2)
    print 37,this%ensid,cendist
    print 38,dt1,dt2,stdev1,stdev2
    print 39,tlat1,tlon1,tlat2,tlon2
    
    call dist_weight(this%nx,this%ny, this%lats,this%lons, &
         this%prob%rdata(:,:,1), this%prob%mask(:,:,1), &
         tlat1,tlon1,tlat2,tlon2, maxdist, cendist,stdev1, &
         cendist,stdev2, minprob)

    nmask=0
    !$OMP PARALLEL DO PRIVATE(i,j) REDUCTION(+:nmask)
    do j=1,this%ny
       do i=1,this%nx
          if(this%prob%mask(i,j,1)) then
             nmask=nmask+1
          endif
       enddo
    enddo
    !$OMP END PARALLEL DO

    last = dt2>this%lasthour-0.05
    if(last) then
104    format('Last ATCF hour was ',F0.1, &
            ' and current end of interval is ',F0.1, &
            '.  Exiting loop for ensid=',I0,'.')
       print 104,this%lasthour,dt2,this%ensid
    endif
  end subroutine track_probability_member
  
  subroutine track_probability_ensemble(this, out,&
       stdev0, sd_per_hour, maxdist, cendist, minprob)
    class(ensemble) :: this
    class(outputter) :: out
    class(member), pointer :: mem
    real, intent(in) :: stdev0, sd_per_hour, maxdist, cendist, minprob
    ! locals
    integer :: ensid,i,j,itrack, notime
    logical :: last, sent(this%ntimes), done(this%ensmin:this%ensmax)
    real :: dt2

    done=.false.
    timeloop: do itrack=1,this%ntimes-1
       ensloop: do ensid=this%ensmin,this%ensmax
          if(.not.this%ensflag(ensid)) cycle ensloop
          if(done(ensid)) cycle ensloop

200       format('Update prob for member ',I0,' at itime=',I0)
          print 200,ensid,itrack
          call this%members(ensid)%track_probability(&
               itrack, itrack==1, done(ensid), dt2, &
               stdev0, sd_per_hour, maxdist, cendist, minprob)
       end do ensloop
       
       output: do notime=1,out%notimes
          if(sent(notime)) then
403          format('Skipping hour ',F0.1,' because I already wrote its output.')
             print 403,out%output_times(notime)
          elseif(dt2+0.05>=out%output_times(notime)) then
404          format('Write data for forecast hour ',F0.1,'.')
             print 404,dt2
             call this%combine_track_prob(nint(dt2),notime,out)
             sent(notime)=.true.
          else
405          format('Hour ',F0.1,' is too early for output time ',F0.1,'.')
             print 405,dt2,out%output_times(notime)
          endif
       enddo output
    enddo timeloop
  end subroutine track_probability_ensemble

  subroutine accum_prob(this,prob)
    class(member) :: this
    class(vardata_real) :: prob
    integer :: i,j
    !$OMP PARALLEL DO PRIVATE(i,j)
    probaccum: do j=1,this%ny
       do i=1,this%nx
          if(this%prob%mask(i,j,1)) then
             if(prob%mask(i,j,1)) then
                prob%rdata(i,j,1) = prob%rdata(i,j,1)+this%prob%rdata(i,j,1)
             else
                prob%mask(i,j,1) = .true.
                prob%rdata(i,j,1) = this%prob%rdata(i,j,1)
             endif
          endif
       enddo
    enddo probaccum
    !$OMP END PARALLEL DO
  end subroutine accum_prob

  subroutine combine_track_prob(this,ifhr,notime,out)
    class(ensemble) :: this
    class(outputter) :: out
    integer, intent(in) :: ifhr
    integer, intent(in) :: notime
    !locals
    integer :: ensid, i, j

    !$OMP PARALLEL DO PRIVATE(i,j)
    do j=1,this%ny
       do i=1,this%nx
          this%prob%mask(i,j,1)=.false.
          this%prob%rdata(i,j,1)=.false.
       enddo
    enddo
    !$OMP END PARALLEL DO

    ensloop: do ensid=this%ensmin,this%ensmax
       if(.not.this%ensflag(ensid)) cycle ensloop
       call out%member_prob(this%prob,notime,ensid)
       !write_pgm pgmone_format
100    format('Accumulate probability from member ',I0)
       call this%members(ensid)%accum_prob(this%prob)
    end do ensloop

    !$OMP PARALLEL DO PRIVATE(i,j)
    probdivide: do j=1,this%ny
       do i=1,this%nx
          if(this%prob%mask(i,j,1)) then
             this%prob%rdata(i,j,1) = this%prob%rdata(i,j,1) &
                  / this%enscount
          endif
       enddo
    enddo probdivide
    !$OMP END PARALLEL DO

    print *,'Write track prob for otime ',notime
    call out%trackprob(this%prob,notime)
    ! write_pgm pgmany_format
    ! write_grads probany itime ivar=1
  end subroutine combine_track_prob

  ! --------------------------------------------------------------------
  
  subroutine most_probable_track(this,ensprob_file, cendist, stdev0, sd_per_hour)
    class(ensemble) :: this
    character(len=filenamelen), intent(in), optional :: ensprob_file
    real, intent(in) :: cendist, stdev0, sd_per_hour
    !locals
    integer :: ensid, ensid1, ensid2, itime, maxntimes, ensmin,ensmax
    real :: hour1, hour2, stdevdist, dist, lasthr1, lasthr2
    logical :: first1,first2
    integer :: ens_winner
    real :: winner_prob, prob

    real :: Pdist(this%ensmin:this%ensmax)
    real :: Ptime(this%ensmin:this%ensmax)
    integer :: Ndist(this%ensmin:this%ensmax)

    maxntimes=this%ntimes
    ensmin=this%ensmin
    ensmax=this%ensmax

    do ensid=this%ensmin,this%ensmax
       if(.not.this%ensflag(ensid)) cycle
90     format('File for ensemble member ',I2,' has last time of ',F0.1)
       print 90,ensid,this%members(ensid)%lasthour
    enddo

    have_multiple: if(this%enscount>1) then
       Ndist=0
       Pdist=0
       Ptime=0
       !$OMP PARALLEL DO PRIVATE(itime,hour1,hour2,ensid1,ensid2, &
       !$OMP                     dist,prob,stdevdist,first1,lasthr1,lasthr2) &
       !$OMP          REDUCTION(+:Ptime) REDUCTION(+:Ndist) REDUCTION(+:Pdist)
       timeloop: do itime=1,maxntimes
33        format('ens range=',I0,'..',I0)
          !print 33,ensmin,ensmax
          first1=.true.
          ens1loop: do ensid1=ensmin,ensmax
             if(.not.this%ensflag(ensid1)) cycle ens1loop
             hour1 = this%members(ensid1)%ftimes(itime)/3600.0
             lasthr1 = this%members(ensid1)%lasthour
             if(first1) then
                first1=.false.
100             format('ENS: Processing time ',F0.1)
                print 100, hour1
             else
1005            format('ENS: Still processing time ',F0.1)
                !print 100,hour1
             endif
             ens2loop: do ensid2=ensmin,ensid1-1
                if(.not.this%ensflag(ensid2)) cycle ens2loop
                hour2 = this%members(ensid2)%ftimes(itime)/3600.0
                lasthr2 = this%members(ensid2)%lasthour
                if(abs(hour1-hour2)>0.05) then
                   write(0,*) 'Mismatch in times between files.  All files must have the same output frequency.'
26                 format('Time index ',I0,' ensid ',I0,' has ',F0.1,' ensid ',I0,' has ',F0.1)
                   write(0,26) itime,ensid1,hour1,ensid2,hour2
                   stop 31
                endif
                if(lasthr1-0.05 < hour1) then
                   if(lasthr2-0.05 < hour2) then
101                   format('ENS ',I2,' vs. ',I2,': both ended: hours ',F0.1,' and ',F0.1)
                      print 101,ensid1,ensid2,hour1,hour2
                      Ptime(ensid1)=Ptime(ensid1)+1
                      Ptime(ensid2)=Ptime(ensid2)+1
                   else
                      ! Mismatch in ending hour, so nothing added to dist or time
102                   format('ENS ',I2,' vs. ',I2,': first ended at hour ',F0.1,' and second has not ended (end ',F0.1,')')
                      print 102,ensid1,ensid2,hour1,lasthr2
                   endif
                elseif(lasthr1-0.05 < hour2) then
                   ! Mismatch in ending hour, so nothing added to dist or time
103                format('ENS ',I2,' vs. ',I2,': first still going (end ',F0.1,') but second ended at hour ',F0.1)
                   print 103,ensid1,ensid2,lasthr1,hour2
                else
                   Ptime(ensid1)=Ptime(ensid1)+1
                   Ptime(ensid2)=Ptime(ensid2)+1
                   dist=greatarc(this%members(ensid1)%track_lat(itime),&
                        this%members(ensid1)%track_lon(itime), &
                        this%members(ensid2)%track_lat(itime), &
                        this%members(ensid2)%track_lon(itime))
3308               format('ENS ',I2,' lat=',F0.3,' lon=',F0.3,' vs. ',I2,' lat=',F0.3,' lon=',F0.3)
                   print 3308,ensid1,this%members(ensid1)%track_lat(itime),this%members(ensid1)%track_lon(itime),&
                        ensid2,this%members(ensid2)%track_lat(itime),this%members(ensid2)%track_lon(itime)
                   Ndist(ensid1)=Ndist(ensid1)+1
                   Ndist(ensid2)=Ndist(ensid2)+1
                   if(dist<cendist) then
                      prob=1.0
104                   format('ENS ',I2,' vs. ',I2,': dist=',F0.1,' within cendist=',F0.1)
                      print 104,ensid1,ensid2,dist,cendist
                   else
                      stdevdist = stdev0 + hour1 * sd_per_hour
                      !prob=exp(- 2*( ( (dist-cendist) / stdevdist ) ** 2 ) )
                      prob=point_prob(dist,cendist,stdevdist)
105                   format('ENS ',I2,' vs. ',I2,': dist=',F0.1,' outside cendist, stdev=',F0.1,' prob=',F0.7)
                      print 105,ensid1,ensid2,dist,stdevdist,prob
                   endif
                   Pdist(ensid1)=Pdist(ensid1)+prob
                   Pdist(ensid2)=Pdist(ensid2)+prob
                endif
             end do ens2loop
          enddo ens1loop
       enddo timeloop
       !$OMP END PARALLEL DO
    endif have_multiple

    if(present(ensprob_file)) then
       open(unit=200,file=trim(ensprob_file),form='FORMATTED')
200    format(F10.7,' probability for ensemble member ',I2)
    endif

    ens_winner=-999
    winner_prob=-9e9
    who_wins: do ensid=ensmin,ensmax
       if(.not.this%ensflag(ensid)) cycle who_wins
       if(this%enscount>1) then
          prob=Pdist(ensid)/max(1,Ndist(ensid))*Ptime(ensid)/max(1,this%enscount-1)/maxntimes
       else
          prob=1.0
       endif
       print 41,ensid,Pdist(ensid),Ndist(ensid),Ptime(ensid),prob
       if(present(ensprob_file)) then
          write(200,200) prob,ensid
       endif
       if(prob>winner_prob) then
          print 42,ensid
          ens_winner=ensid
          winner_prob=prob
       endif
    enddo who_wins
    if(ens_winner==-999) then
       ! Cannot reach this line.
       write(0,*) 'Internal error: no ensemble member had the highest probability track.  (!?)'
       stop 23
    endif
    print 44,ens_winner
    print 45,winner_prob
    if(present(ensprob_file)) then
       close(200)
    endif
41  format('Ensemble member ',I0,' Pdist=',F0.7,' Ndist=',I0,' Ptime=',F0.7,' P=',F0.7)
42  format('And we have a new winner!! Ensemble member ',I0,'!!!')
44  format('Ensemble member ',I0,' won the "Most Probable Member Award"')
45  format('with a probability of ',F0.7)

    this%most_probable_member=ens_winner
  end subroutine most_probable_track

  ! --------------------------------------------------------------------

  subroutine loop_prep_member(this,lastitime)
    class(member) :: this
    integer, intent(in) :: lastitime
    this%lastitime=lastitime

    call reset_rgrid(this%rain,.false.)
    call reset_rgrid(this%wind,.false.)
    call reset_rgrid(this%prob,.false.)
    call reset_rgrid(this%udh03,.false.)
    call reset_rgrid(this%udh25,.false.)
    call reset_rgrid(this%minw,.false.)
    call reset_rgrid(this%maxw,.false.)
    call reset_rgrid(this%maxu10,.false.)
  end subroutine loop_prep_member

  subroutine loop_prep_ensemble(this)
    class(ensemble) :: this
    integer :: ensid, itime
    real :: lastfhr !, smlasthour
    integer :: lastitime, maxntimes, trackdom, swathmem

    maxntimes=this%ntimes
    trackdom=this%trackdom
    swathmem=this%most_probable_member
    this%last_tg_reset=-9e9

    ! Find the last time index that has track information for the
    ! "most probable member".  That is the member for which we are
    ! generating the rainfall swath.
    lastitime=maxntimes
    !smlasthour=this%members(swathmem)%ftimes(1)/3600.0
    find_lastitime: do itime=1,maxntimes
       lastfhr=this%members(swathmem)%ftimes(itime)/3600.0
       if(lastfhr+0.05>=this%members(swathmem)%lasthour) then
          lastitime=itime
          !print 25,itime,lastfhr,swathmem
25        format('Last time for mean rain swath is time #',I0,'=',F0.1, &
               ' from member ',I0)
          exit find_lastitime
       endif
    enddo find_lastitime

    this%lastitime=lastitime

    do ensid=this%ensmin,this%ensmax
       if(.not.this%ensflag(ensid)) cycle
       call this%members(ensid)%loop_prep(lastitime)
    end do

    this%reinit=.false.
  end subroutine loop_prep_ensemble

  ! --------------------------------------------------------------------

  subroutine open_var_member(this,itime,rwork,domid,varname,masktype)
    class(member) :: this
    integer, intent(in) :: itime, masktype, domid
    class(wrfdiag_var) :: rwork
    character*(*), intent(in) :: varname
    !locals
    integer :: i,j

    call rwork%open(this%files(domid),varname,itime)
    call rwork%alloc_mask()

    if(masktype==1) then
       !$OMP PARALLEL DO PRIVATE(i,j)
       do j=rwork%dc%jps,rwork%dc%jpe
          do i=rwork%dc%ips,rwork%dc%ipe
             rwork%mask(i,j,1) = rwork%rdata(i,j,1)>1e-5
          enddo
       enddo
       !$OMP END PARALLEL DO
    elseif(masktype==-1) then
       !$OMP PARALLEL DO PRIVATE(i,j)
       do j=rwork%dc%jps,rwork%dc%jpe
          do i=rwork%dc%ips,rwork%dc%ipe
             rwork%mask(i,j,1) = rwork%rdata(i,j,1)<1e-5
          enddo
       enddo
       !$OMP END PARALLEL DO
    else ! masktype==0
       !$OMP PARALLEL DO PRIVATE(i,j)
       do j=rwork%dc%jps,rwork%dc%jpe
          do i=rwork%dc%ips,rwork%dc%ipe
             rwork%mask(i,j,1) = abs(rwork%rdata(i,j,1))>1e-5
          enddo
       enddo
       !$OMP END PARALLEL DO
    endif
  end subroutine open_var_member

  ! --------------------------------------------------------------------

  subroutine reset_rgrid(var,maskfill)
    class(vardata_real) :: var
    integer :: i,j
    logical, intent(in) :: maskfill
    !$OMP PARALLEL DO PRIVATE(i,j)
    do j=var%dc%jps,var%dc%jpe
       do i=var%dc%ips,var%dc%ipe
          var%mask(i,j,1)=maskfill
          var%rdata(i,j,1)=0.
       enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine reset_rgrid

  subroutine reset_igrid(var,maskfill)
    class(vardata_int) :: var
    logical, intent(in) :: maskfill
    integer :: i,j
    !$OMP PARALLEL DO PRIVATE(i,j)
    do j=var%dc%jps,var%dc%jpe
       do i=var%dc%ips,var%dc%ipe
          var%mask(i,j,1)=maskfill
          var%idata(i,j,1)=0
       enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine reset_igrid

  ! --------------------------------------------------------------------

  subroutine step_time_member(this,itime,rwork,reinit,rainswath,&
       hold,rainaccum,raincount,rainall)
    class(member) :: this
    integer, intent(in) :: itime
    class(wrfdiag_var) :: rwork
    class(vardata_real) :: hold,rainaccum,rainall
    class(vardata_int) :: raincount
    logical, intent(in) :: rainswath
    logical, intent(inout) :: reinit
    ! locals
    integer :: domid,j,i,totalpoints
    real :: totaldiff,raindiff,maxdiff
    type(bilinear_real) :: interp

    if(rainswath) then
       ! Create a copy of the rain swath so we will know how much it
       ! is increased during this update.
35     format('Duplicate ensemble swath for member ',I0)
       !print 35,this%ensid
       !$OMP PARALLEL DO PRIVATE(i,j)
       do j=1,this%ny
          do i=1,this%nx
             hold%mask(i,j,1)=this%rain%mask(i,j,1)
             if(hold%mask(i,j,1)) then
                hold%rdata(i,j,1)=this%rain%rdata(i,j,1)
             else
                hold%rdata(i,j,1)=0
             endif
          enddo
       enddo
       !$OMP END PARALLEL DO
    endif

    do domid=this%mindom,this%maxdom
       ! Open the wind swath variable to get projection information:
       call this%open_var(itime,rwork,domid,'WINDSQ_SWATH',1)

       ! Create an interpolator that will interpolate from the native
       ! grid to the lat-lon grid:
       call pd_init_bilinear_real(interp,this%proj,this%dec,&
            rwork%pj,rwork%dc)

       ! Prepare interpolation arrays:
       call interp%prep()

       ! Interpolate the wind swath to the target grid:
       call interp%scalar(this%wind,rwork,op=OP_MAX)

       ! And free the arrays that stored the native grid mask and data:
       call rwork%free()

       ! Now interpolate other variables...

       call this%open_var(itime,rwork,domid,'PRECIP_SWATH',1)
       call interp%scalar(this%rain,rwork,op=OP_MAX)
       call rwork%free()

       call this%open_var(itime,rwork,domid,'TG_MAX_UPDHEL_03',1)
       call interp%scalar(this%udh03,rwork,op=OP_MAX)
       call rwork%free()

       call this%open_var(itime,rwork,domid,'TG_MAX_UPDHEL_25',1)
       call interp%scalar(this%udh25,rwork,op=OP_MAX)
       call rwork%free()

       call this%open_var(itime,rwork,domid,'TG_MIN_WWIND',-1)
       call interp%scalar(this%minw,rwork,op=OP_MIN)
       call rwork%free()

       call this%open_var(itime,rwork,domid,'TG_MAX_WWIND',1)
       call interp%scalar(this%maxw,rwork,op=OP_MAX)
       call rwork%free()

       call this%open_var(itime,rwork,domid,'TG_MAX_M10WIND',1)
       call interp%scalar(this%maxu10,rwork,op=OP_MAX)
       call rwork%free()

       ! Free the (huge) arrays associated with the interpolator:
       call interp%free()
    enddo

    if(rainswath) then
       !print *,'UPDATE RAINALL'
       totaldiff=0
       totalpoints=0
       maxdiff=0
       !$OMP PARALLEL DO PRIVATE(i,j,raindiff) &
       !$OMP   REDUCTION(+:totaldiff) REDUCTION(+:totalpoints) &
       !$OMP   REDUCTION(max:maxdiff)
       do j=1,this%ny
          do i=1,this%nx
             raindiff = this%rain%rdata(i,j,1) - hold%rdata(i,j,1)
             if(raindiff>0) then
                rainaccum%mask(i,j,1)  = .true.
                rainaccum%rdata(i,j,1) = rainaccum%rdata(i,j,1) + raindiff
                raincount%mask(i,j,1)  = .true.
                raincount%idata(i,j,1) = raincount%idata(i,j,1) + 1
                totaldiff              = totaldiff+raindiff
                totalpoints            = totalpoints+1
                maxdiff                = max(maxdiff,raindiff)
             elseif(raindiff<-1e-5) then
                write(0,13012) raindiff,i,j
13012           format('ERROR: negative accumulation ',F0.5,' from i=',I0,' j=',I0)
                stop 26
             endif
          enddo
       enddo
       !$OMP END PARALLEL DO

12013  format('ENS ',I0,' itime=',I0,': total accumulation is ',F0.5,&
              ' from ',I0,' points with max = ',F0.5,'.')
       print 12013,this%ensid,itime,totaldiff,totalpoints,maxdiff
    endif
  end subroutine step_time_member

  subroutine step_time_ensemble(this,itime,rainswath,&
       okdist0,okdist_per_hour)
    class(ensemble) :: this
    integer, intent(in) :: itime
    type(wrfdiag_var) :: rwork
    logical, intent(in) :: rainswath
    real, intent(in) :: okdist0,okdist_per_hour
    ! Locals
    logical :: reinit
    real :: swathlat, swathlon, rfhr, okdist, dist, raindiff, &
            totaldiff, maxdiff
    integer :: ifhr, swathmem, ensid, nswath, totalpoints, i,j
    logical :: membswath(this%ensmin:this%ensmax), xrainswath

    swathmem=this%most_probable_member

30  format('Process time ',I0)
    print 30,itime
    
314 format('The swathmem=',I0,' and trackdom=',I0)
    print 314,swathmem,this%trackdom
    
    ! Get the forecast hour at this time.  We get it from the
    ! swathmem:
    xrainswath=itime<this%lastitime.and.rainswath
    rfhr=this%members(swathmem)%ftimes(itime)/3600.0
    ifhr=nint(rfhr)
    
    if(rainswath) then
       okdist=okdist0+okdist_per_hour*rfhr
       if(okdist<1e3) then
315       format('ERROR: okdist<1km.  That is inadvisable due to numerical error.  okdist=',F0.3,'=',F0.3,'+',F0.3,'*',F0.3)
          print 315,okdist,okdist0,okdist_per_hour,rfhr
          write(0,315) okdist,okdist0,okdist_per_hour,rfhr
          stop 23
       else
          print 316,okdist,okdist0,okdist_per_hour,rfhr
316       format('Using okdist=',F0.3,'=',F0.3,'+',F0.3,'*',F0.3)
       endif

       swathlat=this%members(swathmem)%track_lat(itime)
       swathlon=this%members(swathmem)%track_lon(itime)
       nswath=0
       do ensid=this%ensmin,this%ensmax
          if(.not.this%ensflag(ensid)) cycle
4048      format('ens ',I0,' at lat=',F0.3,' lon=',F0.3, &
               ' vs. swath ens ',I0,' at lat=',F0.3,' lon=',F0.3)
          print 4048,ensid,this%members(ensid)%track_lat(itime), &
               this%members(ensid)%track_lon(itime), &
               swathmem,swathlat,swathlon
          dist=greatarc(swathlat,swathlon,&
               this%members(ensid)%track_lat(itime), &
               this%members(ensid)%track_lon(itime))
          membswath(ensid) = dist<=okdist
          nswath=nswath+1
32        format('Member ',I0,' is ',F0.1,' meters from most probable member ',I0)
          !print 32,ensid,dist,swathmem
          if(membswath(ensid)) then
33           format('Include member ',I0,' in mean swath at itime=',I0,'.')
             print 33,ensid,itime
          else
             print 34,ensid,itime
34           format('Do NOT include member ',I0,' in mean swath at itime=',I0,'.')
          endif
       end do

       call reset_rgrid(this%rainaccum,.false.)
       call reset_igrid(this%raincount,.false.)
    endif

    reinit = this%reinit
    do ensid=this%ensmin,this%ensmax
       if(.not.this%ensflag(ensid)) cycle
       call this%members(ensid)%step_time(itime,rwork,reinit, &
            membswath(ensid),this%hold,this%rainaccum, &
            this%raincount,this%rainall)
       this%reinit=.true.
    end do
    
    this%reinit=reinit

    if(rainswath) then    
       ! Update the mean rainfall swath from the rainaccum and
       ! raincount arrays
       maxdiff=0
       totaldiff=0
       totalpoints=0
       !$OMP PARALLEL DO PRIVATE(i,j,raindiff) &
       !$OMP REDUCTION(+:totaldiff) REDUCTION(+:totalpoints) &
       !$OMP REDUCTION(max:maxdiff)
       do j=1,this%ny
          do i=1,this%nx
             if(this%raincount%mask(i,j,1) .and. this%rainaccum%mask(i,j,1)) then
                if(this%raincount%idata(i,j,1)>0 .and. this%rainaccum%rdata(i,j,1)>0) then
                   this%rainall%mask(i,j,1)=.true.
                   raindiff=this%rainaccum%rdata(i,j,1)/this%raincount%idata(i,j,1)
                   totaldiff=totaldiff+raindiff
                   totalpoints=totalpoints+1
                   this%rainall%rdata(i,j,1) = this%rainall%rdata(i,j,1) + raindiff
                   maxdiff=max(maxdiff,raindiff)
                endif
             endif
          enddo
       enddo
       !$OMP END PARALLEL DO
95     format('Updated rain on ',I0,' points with total diff = ',F0.5, &
            'and max diff = ',F0.5)
       print 95,totalpoints,totaldiff,maxdiff
    endif

    if(reinit) then
       call rwork%free()
    end if
  end subroutine step_time_ensemble

  ! --------------------------------------------------------------------

  subroutine output_swath(this,out,itime)
    class(ensemble) :: this
    class(outputter) :: out
    integer, intent(in) :: itime
    ! Locals
    logical :: outflag, rainswath
    integer :: chosen_otime, otime, nthresh, ithresh, itype, i,j, swathmem, ifhr
    real, pointer :: thresh(:)
    class(vardata_real), pointer :: indata
    real :: cutoff,rfhr
    integer :: ensid

    if(.not.associated(this%sent)) then
       allocate(this%sent(out%notimes))
       this%sent=.false.
    endif

    swathmem=this%most_probable_member

30  format('Swath output time ',I0)
    print 30,itime
    
314 format('The swathmem=',I0,' and trackdom=',I0)
    print 314,swathmem,this%trackdom
    
    ! Get the forecast hour at this time.  We get it from the
    ! swathmem:
    rainswath=itime<this%lastitime
    rfhr=this%members(swathmem)%ftimes(itime)/3600.0
    ifhr=nint(rfhr)

100 format(A)
    print 100,'DECIDE WHETHER TO OUTPUT'
    outflag=.false.
    chosen_otime=-999
    do otime=1,out%notimes
       if(this%sent(otime)) then
403       format('Skipping hour ',F0.1,' because I already wrote its output.')
          !print 403,out%output_times(otime)
       elseif(rfhr+0.05>=out%output_times(otime)) then
          if(chosen_otime<0) chosen_otime=otime
          print 404,rfhr
          outflag=.true.
          this%sent(otime)=.true.
404       format('Write data for forecast hour ',F0.1,'.')
       else
          print 405,rfhr,out%output_times(otime)
405       format('Hour ',F0.1,' is too early for output time ',F0.1,'.')
       endif
    enddo
    if(.not.outflag) then
58     format('No output to do at time ',F0.1)
       print 58,rfhr
       return
    endif

    call out%rainswath(this%rainall,chosen_otime)

    ! LOOP OVER ALL THRESHOLDS
    typeloop: do itype=1,2
       if(itype==1) then
          nthresh=out%nwthresh
          thresh=>out%wthresh
       else
          nthresh=out%nrthresh
          thresh=>out%rthresh
       endif
       wthreshloop: do ithresh=1,nthresh
          if(thresh(ithresh)<-998) then
49           format('Invalid threshold type ',I0,' #',I0,': ',F0.5,'<-998')
             print 49,itype,ithresh,thresh(ithresh)
             write(0,49) itype,ithresh,thresh(ithresh) 
             stop 27
          endif
          if(itype==1) then
             cutoff=thresh(ithresh)**2
          else
             cutoff=thresh(ithresh)
          endif
50        format('Threshold type ',I0,' number ',I0,' is ',F0.1,': process wthreshold')
          print 50, itype,ithresh,cutoff
          ! STEP 2: CLEAR COMBINED PROBABILITY ARRAY
          !$OMP PARALLEL DO PRIVATE(i,j)
          do j=1,this%ny
             do i=1,this%nx
                this%prob%mask(i,j,1)=.false.
                this%prob%rdata(i,j,1)=0.0
             enddo
          enddo
          !$OMP END PARALLEL DO
          
          ! STEP 3: MAKE PROBABILITY FOR EACH ENSEMBLE MEMBER
          do ensid=this%ensmin,this%ensmax
             if(.not.this%ensflag(ensid)) cycle
             if(itype==1) then
                indata=>this%members(ensid)%wind
             else
                indata=>this%members(ensid)%rain
             endif
             !$OMP PARALLEL DO PRIVATE(i,j)
             do j=1,this%ny
                do i=1,this%nx
                   if(indata%mask(i,j,1)) then
                      if(indata%rdata(i,j,1)>cutoff) then
                         this%prob%mask(i,j,1)=.true.
                         this%prob%rdata(i,j,1)=this%prob%rdata(i,j,1)+1
                      endif
                   endif
                enddo
             enddo
             !$OMP END PARALLEL DO
          enddo

          ! STEP 4: DIVIDE COMBINED PROBABILITY BY ENSCOUNT
          !$OMP PARALLEL DO PRIVATE(i,j)
          do j=1,this%ny
             do i=1,this%nx
                this%prob%rdata(i,j,1)=this%prob%rdata(i,j,1)/this%enscount
             enddo
          enddo
          !$OMP END PARALLEL DO

          if(itype==1) then
             call out%windprob(this%prob,chosen_otime,ithresh,thresh(ithresh))
          else
             call out%rainprob(this%prob,chosen_otime,ithresh,thresh(ithresh))
          endif
       enddo wthreshloop
    enddo typeloop
  end subroutine output_swath

  ! --------------------------------------------------------------------

  subroutine output_tornado(this,out,itgtime,rtime)
    class(ensemble) :: this
    class(outputter) :: out
    integer, intent(in) :: itgtime
    real, intent(in) :: rtime
    ! Locals
    class(vardata_real), pointer :: indata
    integer :: ivar,ensid,op,i,j

30  format('Outputting tornado data for time #',I0,' = ',F0.1,'hrs')
    print 30,itgtime,rtime
    
    do ivar=1,5
       call reset_rgrid(this%extrema,.false.)
       do ensid=this%ensmin,this%ensmax
          if(.not.this%ensflag(ensid)) cycle
          op=OP_MAX
          select case(ivar)
          case(1) ; indata=>this%members(ensid)%udh03
          case(2) ; indata=>this%members(ensid)%udh25
          case(3) ; indata=>this%members(ensid)%maxw
          case(4) ; indata=>this%members(ensid)%maxu10
          case(5)
             indata=>this%members(ensid)%minw
             op=OP_MIN
          end select

          if(op==OP_MAX) then
             !$OMP PARALLEL DO PRIVATE(i,j)
             do j=1,this%ny
                do i=1,this%nx
                   if(indata%mask(i,j,1)) then
                      this%extrema%mask(i,j,1)=.true.
                      this%extrema%rdata(i,j,1)=max(this%extrema%rdata(i,j,1),&
                           indata%rdata(i,j,1))
                   endif
                end do
             end do
             !$OMP END PARALLEL DO
          else
             !$OMP PARALLEL DO PRIVATE(i,j)
             do j=1,this%ny
                do i=1,this%nx
                   if(indata%mask(i,j,1)) then
                      this%extrema%mask(i,j,1)=.true.
                      this%extrema%rdata(i,j,1)=min(this%extrema%rdata(i,j,1),&
                           indata%rdata(i,j,1))
                   endif
                end do
             end do
             !$OMP END PARALLEL DO
          endif

          select case(ivar)
          case(1) ; call out%udhel(this%extrema,itgtime,0.,3000.)
          case(2) ; call out%udhel(this%extrema,itgtime,2000.,5000.)
          case(3) ; call out%maxw(this%extrema,itgtime)
          case(4) ; call out%maxu10(this%extrema,itgtime)
          case(5) ; call out%minw(this%extrema,itgtime)
          end select

          call reset_rgrid(indata,.false.)
       enddo
    enddo
  end subroutine output_tornado

  ! --------------------------------------------------------------------

  subroutine run(this,out, &
       stdev0, sd_per_hour, maxdist, cendist, minprob,&
       ensprob_file, rainswath, okdist0, okdist_per_hour)
    real, intent(in) :: stdev0, sd_per_hour, maxdist, cendist, minprob, &
         okdist0, okdist_per_hour

    class(ensemble) :: this
    class(outputter) :: out
    character(len=filenamelen), intent(in), optional :: ensprob_file
    logical, intent(in), optional :: rainswath
    ! locals
    integer :: itime,itgtime
    real :: last_tg_time, rtime
    logical :: lrainswath, is_tg_time

    lrainswath=.true.
    if(present(rainswath)) lrainswath=rainswath

    ! Open wrfdiag files for valid times:
    call this%open()
    
    ! Read tracks: both ATCF and wrfdiag hourly internal track:
    call this%read_tracks()

    ! Get the bounds for the lat-lon grid:
    call this%bound()

    ! Initialize work arrays:
    call this%init_arrays()

    ! Prepare the outputter:
    call out%start(this%lats,this%lons,this%atime,this%nx,this%ny,&
         (this%ftimes(this%ntimes)-this%ftimes(1))/3600.0_8)

    ! Determine track probability:
    call this%track_probability(out,&
         stdev0, sd_per_hour, maxdist, cendist, minprob)

    ! Find the "most probable track:"
    call this%most_probable_track(ensprob_file, &
         cendist, stdev0, sd_per_hour)

    ! Prepare for the time loop:
    call this%loop_prep()

    last_tg_time=0
    itgtime=0
    ! Loop over all times
    do itime=1,this%ntimes
       rtime=this%ftimes(itime)/3600.0
       is_tg_time = rtime-out%tg_freq>last_tg_time-0.05

       ! Inform the outputter that the time is increasing:
       call out%step_time(itime,rtime)

       ! Read the next time's data and update accumulator arrays:
       call this%step_time(itime,lrainswath,okdist0,okdist_per_hour)

       ! Output swath and per-member information:
       call this%output_swath(out,itime)
       
       if(is_tg_time) then
          itgtime=itgtime+1
          call this%output_tornado(out,itgtime,rtime)
          last_tg_time=rtime
       endif
    end do

    ! Flush output:
    call out%finish()
  end subroutine run

    ! ----------------------------------------------------------------------------

    subroutine insert_ens_domain(inpath,outpath,ensid,domain,ifhr,ithresh,&
         var,name,subname)
      character(len=filenamelen), intent(in) :: inpath
      character(len=filenamelen), intent(out) :: outpath
      character*(*), intent(in), optional :: var, name, subname
      integer, intent(in), optional :: ifhr, ithresh
      integer :: i,o,n,varlen
      integer, intent(in), optional :: ensid,domain
      integer :: not_initialized
38    format('Insert ens, domain, etc. into "',A,'"')
      !print 38,trim(inpath)
      n=len_trim(inpath)
      i=1
      o=1
      do while(i<=n)
         if(inpath(i:i)=='<') then
            ! <thresh> = 8 chars
            ! <domain> = 8 chars
            ! <ensid> = 7 chars
            ! <ifhr> = 6 chars
            ! <var> = 5 chars
            if(present(var) .and. i+4<=n) then
               if(inpath(i:(4+i))=='<var>') then
                  varlen=len_trim(var)
                  outpath(o:(varlen+o-1))=var(1:varlen)
                  i=i+5
                  o=o+varlen
                  cycle
               endif
            endif
            if(i+5<=n) then
               if(present(ifhr) .and. inpath(i:(5+i))=='<ifhr>') then
                  if(ifhr<10) then
                     write(outpath(o:(2+o)),'("00",I1)') ifhr
                  elseif(ifhr<100) then
                     write(outpath(o:(2+o)),'("0",I2)') ifhr
                  else
                     write(outpath(o:(2+o)),'(I3)') ifhr
                  endif
                  i=i+6
                  o=o+3
                  cycle
               elseif(present(name) .and. inpath(i:(5+i))=='<name>') then
                  varlen=len_trim(name)
                  outpath(o:(varlen+o-1))=name(1:varlen)
                  i=i+6
                  o=o+varlen
               endif
            endif
            if(i+6<=n) then
               if(present(ensid) .and. inpath(i:(6+i))=='<ensid>') then
                  if(ensid<10) then
                     write(outpath(o:(1+o)),'("0",I1)') ensid
                  else
                     write(outpath(o:(1+o)),'(I2)') ensid
                  endif
                  i=i+7
                  o=o+2
                  cycle
               endif
            endif
            if(i+7<=n) then
               if(present(ensid) .and. inpath(i:(7+i))=='<ensid3>') then
                  if(ensid<10) then
                     write(outpath(o:(2+o)),'("00",I1)') ensid
                  elseif(ensid<100) then
                     write(outpath(o:(2+o)),'("0",I2)') ensid
                  else
                     write(outpath(o:(2+o)),'(I3)') ensid
                  endif
                  i=i+8
                  o=o+3
                  cycle
               elseif(present(ensid) .and. inpath(i:(7+i))=='<ensid2>') then
                  if(ensid<10) then
                     write(outpath(o:(1+o)),'("0",I1)') ensid
                  else
                     write(outpath(o:(1+o)),'(I2)') ensid
                  endif
                  i=i+8
                  o=o+2
                  cycle
               elseif(present(domain) .and. inpath(i:(7+i))=='<domain>') then
                  if(domain<10) then
                     write(outpath(o:(1+o)),'("0",I1)') domain
                  else
                     write(outpath(o:(1+o)),'(I2)') domain
                  endif
                  i=i+8
                  o=o+2
                  cycle
               elseif(present(ithresh) .and. inpath(i:(7+i))=='<thresh>') then
                  if(ithresh<10) then
                     write(outpath(o:(2+o)),'("00",I1)') ithresh
                  elseif(ithresh<100) then
                     write(outpath(o:(2+o)),'("0",I2)') ithresh
                  else
                     write(outpath(o:(2+o)),'(I3)') ithresh
                  endif
                  i=i+8
                  o=o+3
                  cycle
               endif
            elseif(present(subname)) then
               if(i+8<=n .and. inpath(i:(i+8))=='<subname>') then
                  varlen=len_trim(subname)
                  outpath(o:(varlen+o-1))=subname(1:varlen)
                  i=i+8
                  o=o+varlen
               endif
            endif
31          format('Unexpected "<" found in "',A,'" near "',A,'"')
            write(0,31) trim(inpath), trim(inpath(i:min(i+7,n)))
            print *, not_initialized
            stop 9
         else
            outpath(o:o)=inpath(i:i)
            o=o+1
            i=i+1
         end if
      enddo
      outpath(o:)=' '
39    format('   ... result="',A,'"')
      !print 39,trim(outpath)
    end subroutine insert_ens_domain

    ! ----------------------------------------------------------------------------

    real function point_prob(dist,cendist,stdevdist)
      real, intent(in) :: dist,cendist,stdevdist

      real :: x1,x2

      x1=(dist-cendist)/stdevdist
      x2=(dist+cendist)/stdevdist

      point_prob=min(1.0,max(0.0,(erfc(-x2/sqrt(2.))-erfc(-x1/sqrt(2.)))/2))
    end function point_prob

    subroutine dist_weight(nx,ny, lats,lons, out,mask, l1lat,l1lon, &
         l2lat,l2lon, maxdist, cendist1, stdevdist1, &
         cendist2, stdevdist2, minprob)
      use constants_module, only: Rearth,pi
      integer :: i,j
      integer, intent(in) :: nx,ny
      real, intent(in) :: lats(ny), lons(nx)
      real, intent(inout) :: out(nx,ny)
      logical, intent(out) :: mask(nx,ny)
      real, intent(in) :: l1lat,l1lon, l2lat,l2lon, maxdist, minprob
      real, intent(in) :: cendist1, stdevdist1, cendist2, stdevdist2

      real :: dx1,dx2,dy1,dy2, xlen,ylen, len2, xmid,ymid, t, dist
      real :: px,py, lx1,ly1,lx2,ly2, l1dist, xlenlen, ylenlen
      real :: cendist, stdevdist, prob
      real, parameter :: s2pi=2.50662827463 ! sqrt(2*pi)
      integer :: istart,istop
      real :: ilonmax,ilonmin,norm1,norm2

10    format('Get dist from P=(',F0.1,',',F0.1,') to L1=(',F0.1,',',F0.1, &
           ') L2=(',F0.1,',',F0.1,').')

11    format('Point is near L1 with dist=',F0.2)
12    format('Point is near L2 with dist=',F0.2)
13    format('Point is not near ends.  The dist=',F0.2)

20    format('L1=L2 so get dist to L1: dist=',F0.2)
40    format('Too far from segment.  Use maxdist: dist=',F0.2)
50    format('Within maxdist.  Keep dist: dist=',F0.2)

      lx1 = l1lon
      ly1 = l1lat

      lx2 = mod(3600.0 + l2lon - l1lon + 180, 360.0) + l1lon - 180
      ly2 = l2lat

      xlen=lx2-lx1
      ylen=ly2-ly1
      len2=xlen*xlen+ylen*ylen
      xlenlen=xlen/len2
      ylenlen=ylen/len2

      norm1=erf(cendist1/stdevdist1)
      norm2=erf(cendist2/stdevdist2)

      ilonmin=(min(lx1,lx2)-lons(1))/(lons(2)-lons(1))+1
      ilonmax=(max(lx1,lx2)-lons(1))/(lons(2)-lons(1))+1

55    format('dist_weight: lx1=',F0.3,' lx2=',F0.3,' ilonmin=',F0.3,&
             ' ilonmax=',F0.3,' lons(1)=',F0.3,' norm1=',F0.5,' norm2=',F0.5)
      print 55,lx1,lx2,ilonmin,ilonmax,lons(1),norm1,norm2

      if(len2<dgepsilon2) then
         ! The "line segment" is a point.  Get the distance to that point.
         cendist=max(cendist1,cendist2)
         stdevdist=max(stdevdist1,stdevdist2)
         !$OMP PARALLEL DO PRIVATE(i,j,px,py,dist,prob,istart,istop)
         do j=1,ny

            dist=abs(maxdist*1.1*cosd(lats(j))/Rearth*180/pi / (lons(2)-lons(1)))
            istart=min(nx,max(1, floor(ilonmin-dist) ))
            istop=min(nx,max(1, floor(ilonmax+dist) ))
            
            do i=1,nx ! istart,istop
               px  = mod(3600.0 +  lons(i) - l1lon + 180, 360.0) + l1lon - 180
               py  = lats(j)
               !print 10,px,py,lx1,ly1,lx2,ly2
               dist=greatarc(py,px,ly1,lx1)
               if(dist>maxdist) then
                  out(i,j)=max(out(i,j),0.0)
               else
                  !prob=exp(- 2*( ( (dist-cendist) / stdevdist ) ** 2 ) ) * &
                  !     erf(cendist/stdevdist)
                  prob=point_prob(dist,cendist,stdevdist)
                  if(prob<minprob) then
                     out(i,j)=max(out(i,j),0.0)
                  else
                     out(i,j)=max(out(i,j),prob)
                     mask(i,j)=.true.
                  endif
               endif
               !print 20,dist
            enddo
         enddo
         !$OMP END PARALLEL DO
         return
      endif

      !$OMP PARALLEL DO PRIVATE(i,j,px,py,dx1,dy1,t,dist,cendist,stdevdist,dx2,&
      !$OMP                     dy2,prob,istart,istop)
      do j=1,ny

         dist=abs(maxdist*1.1*cosd(lats(j))/Rearth*180/pi / (lons(2)-lons(1)))
         istart=min(nx,max(1, floor(ilonmin-dist) ))
         istop=min(nx,max(1, floor(ilonmax+dist) ))

606      format('At j=',I4,', dist=',F0.3,' i=',I0,'...',I0)
         !print 606,j,dist,istart,istop

         do i=1,nx ! istart,istop
            px  = mod(3600.0 + lons(i) - l1lon + 180, 360.0) + l1lon - 180
            py  = lats(j)

            !print 10,px,py,lx1,ly1,lx2,ly2
            dx1=px-lx1
            dy1=py-ly1

            ! Parameterize the line going through the segment as l1+t*(l2-l1):
            t = dx1*xlenlen + dy1*ylenlen

            ! Are we past the end of the line segment?
            if(t<0) then
               dist=greatarc(py,px,ly1,lx1)
               !print 11,dist
               cendist=cendist1
               stdevdist=stdevdist1
            else
               dx2=px-lx2
               dy2=py-ly2
               if(t>1) then
                  dist=greatarc(py,px,ly2,lx2)
                  !print 12,dist
                  cendist=cendist2
                  stdevdist=stdevdist2
               else
                  ! Point P is closest to a point on the segment that is
                  ! betwen point 1 and point 2.  Project point P to nearest
                  ! point on the line segment, and get the distance between
                  ! the two.
                  dist=greatarc(py,px, &
                       ly1 + t * ylen, &
                       lx1 + t * xlen)
                  cendist = cendist1 + t*(cendist2-cendist1)
                  stdevdist = stdevdist1 + t*(stdevdist2-stdevdist1)
                  !print 13,dist
               endif
            endif

            if(dist>maxdist) then
               ! The line is too far from the maximum distance
               dist=maxdist
               out(i,j)=max(out(i,j),0.0)
               !print 40,dist
               cycle
            else
               !print 50,dist
               !prob=exp(- 2*( ( (dist-cendist) / stdevdist ) ** 2 ) ) * &
               !     erf(cendist/stdevdist)
               prob=point_prob(dist,cendist,stdevdist)
               if(prob<minprob) then
                  out(i,j)=max(out(i,j),0.0)
               else
                  out(i,j)=max(out(i,j),prob)
                  mask(i,j)=.true.
               endif
            endif
         enddo
      enddo
      !$OMP END PARALLEL DO
    end subroutine dist_weight


end module ens_prob_module
