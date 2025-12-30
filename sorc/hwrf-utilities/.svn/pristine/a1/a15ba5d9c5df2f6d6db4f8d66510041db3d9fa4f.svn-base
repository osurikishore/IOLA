module swathgen_module
  use, intrinsic :: iso_c_binding, only: C_DOUBLE, C_INT64_T
  use sysutil_module, only: fail, warn
  use wrfdiag_module!, only: wrfdiag_file, init_wrfdiag_file, wrfdiag_var
  use vardata_module!, only: vardata_real, init_vardata_real
  use projection_module !, only: proj_nmme, proj_latlon, init_proj_latlon
  use decomp_module !, only: decomp
  use interp_module !, only: bilinear_real, init_bilinear_real

  implicit none

  private
  public :: swathgen, init_swathgen

  integer , parameter :: ix_wind=1, ix_precip=2, ix_suswind=3, &
       ix_cuprecip=4, ix_maxw=5, ix_minw=6, ix_updhel25=7,     &
       ix_updhel03=8, num_vars=8

  integer, parameter :: interp_op(num_vars) = (/ &
       OP_MAX, OP_MAX, OP_MAX, &
       OP_MAX, OP_MAX, OP_MIN, OP_MAX, &
       OP_MAX/)

  logical, parameter :: swath_clear(num_vars) = &
       (/ .false., .false., .false., .false., & ! do not clear swath vars
          .true., .true., .true., .true. /)     ! clear hourly vars

  !    '123456789012345678901234567890'
  character(len=30), parameter :: swath_wrfname(num_vars) = (/  &
       'WINDSQ_SWATH                  ', & !  1 = ix_wind
       'PRECIP_SWATH                  ', & !  2 = ix_precip
       'SUSWIND_SWATH                 ', & !  3 = ix_suswind
       'CUPRECIP_SWATH                ', & !  4 = ix_cuprecip
       'TG_MAX_WWIND                  ', & !  5 = ix_maxw
       'TG_MIN_WWIND                  ', & !  6 = ix_minw
       'TG_MAX_UPDHEL_25              ', & !  7 = ix_updhel25
       'TG_MAX_UPDHEL_03              ' /) !  8 = ix_updhel03

  real, parameter :: mps2knots=3600./1852. ! exact
  real, parameter :: meters2inches=100./2.54 ! exact

  real, parameter :: swath_to_crazy_us_units(num_vars) = (/ &
       mps2knots                       , & !  1 = wind: m/s => knots
       meters2inches                   , & !  2 = precip: m => inches
       mps2knots                       , & !  3 = sus.wind: m/s => knots
       meters2inches                   , & !  4 = cu.precip: m => inches
       mps2knots                       , & !  5 = max w: m/s => knots
       mps2knots                       , & !  6 = min w: m/s => knots
       1.                              , & !  7 = updraft hel (no conv.)
       1.                               /) !  8 = updraft hel (no conv.)

  type swathgen
     type(wrfdiag_file), pointer :: files(:) => NULL()
     type(bilinear_real), pointer :: interp(:) => NULL()

     ! ALL vars for wrfdiag input:
     type(wrfdiag_var), pointer :: filevars(:,:)=>NULL()

     ! ALL vars for GRIB2 output:
     type(vardata_real), pointer :: datavars(:)=>NULL()

     ! Lat and lon vars for swath:
     type(vardata_real) :: lats,lons

     ! Projection information:
     type(proj_latlon) :: proj
     type(decomp) :: dec

     ! Track information:
     real, pointer :: track_lat(:)=>NULL(), track_lon(:)=>NULL(), &
          track_wind(:)=>NULL()

     real :: file_dtime, nauxphr
     integer :: maxnfile=0, nfiles=0, itime=0, ntimes=0, nx=0, ny=0
     logical :: have_proj=.false., inited=.false., did_first_interp=.false.
     logical :: in_crazy_us_units=.false., have_track=.false.

     ! GRIB write information:
     logical :: write_grib=.false.
     integer :: grib_unit=0
     integer(kind=8) :: grib_pos=0
   contains
     procedure :: add_file
     procedure, private :: first_file
     procedure :: open_grib2
     procedure :: write_grib2
     procedure :: free => free_swathgen
     procedure :: setup_proj
     procedure :: interp_time
     procedure :: write
     procedure :: to_crazy_us_units
     procedure :: discard_track
     procedure :: get_track
  end type swathgen

contains

  subroutine to_crazy_us_units(this)
    ! Convert winds to knots and precip to inches
    class(swathgen), intent(inout) :: this
    integer :: i,j, ivar

    if(this%in_crazy_us_units) return
    do ivar=1,num_vars
       !$OMP PARALLEL DO PRIVATE(i,j)
       do  j = 1 , this%ny
          do  i = 1 , this%nx
             this%datavars(ivar)%rdata(i,j,1) = &
                  this%datavars(ivar)%rdata(i,j,1) * &
                  swath_to_crazy_us_units(ivar)
          enddo
       enddo
    !$OMP END PARALLEL DO
    enddo
    this%in_crazy_us_units=.true.
  end subroutine to_crazy_us_units

  subroutine discard_track(this)
    class(swathgen), intent(inout) :: this
    ! Do not deallocate track: it is allocated in another object.
    nullify(this%track_lat,this%track_wind,this%track_lon)
  end subroutine discard_track

  subroutine get_track(this)
    class(swathgen), intent(inout) :: this
    integer :: i
    if(this%have_track) return
    call this%discard_track()
    fileloop: do i=1,this%nfiles
       if(this%files(i)%have_track) then
          if(size(this%files(i)%track_lat)>=this%ntimes) then
             this%track_lat=>this%files(i)%track_lat
             this%track_lon=>this%files(i)%track_lon
             this%track_wind=>this%files(i)%track_wind
             this%have_track=.true.
             exit fileloop
          endif
       endif
    enddo fileloop

    if(.not. this%have_track) then
       call fail('swathgen: in get_track: no input files have valid internal HWRF tracker information.  Did you forget to include domain 3?  Or did you forget to set vortex_tracker=6 or 7 in your namelist.input?')
    endif

  end subroutine get_track

  ! --------------------------------------------------------------------

  subroutine open_grib2(this,filename,iostat)
    class(swathgen), intent(inout) :: this
    integer, intent(inout) :: iostat
    character(len=*), intent(in) :: filename

    integer :: unit

    iostat=-999

    if(this%write_grib) then
       write(0,*) 'FATAL ERROR: TRIED TO OPEN TWO OUTPUT GRIB2 FILES FROM THE SAME SWATH GENERATOR'
       return
    endif
    open(newunit=unit,iostat=iostat,access='STREAM',file=filename,&
         status='UNKNOWN',form='UNFORMATTED')
303 format(A,': ',A)
404 format(A,': ',A,I0)
    if(iostat/=0) then
       write(0,404) filename,'cannot open for write, status=',iostat
       return
    else
       print 303, filename,'write GRIB2 swath here'
    endif

    this%grib_unit=unit
    this%grib_pos=1
    this%write_grib=.true.

33  format(A,'open for write (truncate) unit=',I0,' pos=',I0)
    print 33, filename,this%grib_unit,this%grib_pos
  end subroutine open_grib2

  ! --------------------------------------------------------------------

  subroutine write_grib2(this, ierr)
    class(swathgen), intent(inout) :: this
    integer, intent(inout) :: ierr

    character(len=:), allocatable :: cgrib(:), writeme

    integer :: lengribs(num_vars)
    integer(kind=8) :: lenwrite, writepos(num_vars), writesize, pos
    type(vardata_real), pointer :: v
    real(kind=C_DOUBLE) :: atime, ftime, window
    integer :: lcgrib, nx, ny, lengrib, ivar, fitime, i, j
    real :: rmiss, lat0, lon0, dlat, dlon
    logical :: lon_inner

    if(.not.this%write_grib) return

    atime=this%files(1)%epochtimes(1)
    fitime=this%filevars(1,1)%itime
    ftime=this%files(1)%epochtimes(fitime)
    window=this%files(1)%epochtimes(2)-atime
    nx=this%nx
    ny=this%ny
    lon0=this%proj%lon1
    lat0=this%proj%lat1
    dlon=this%proj%dlon
    dlat=this%proj%dlat
    lon_inner=this%proj%lon_inner
    rmiss=-9999.0
    lcgrib=nx*ny*3+1000

    allocate(character(len=lcgrib) :: cgrib(num_vars))

    lenwrite=0
    ierr=0
    !$OMP PARALLEL DO PRIVATE(i,j,ivar,lengrib,v) &
    !$OMP&  REDUCTION(max:ierr) REDUCTION(+:lenwrite) DEFAULT(SHARED)
    do ivar=1,num_vars
       v=>this%datavars(ivar)
       ierr=0
       call impl_make_grib2(atime,ftime,window,ivar,nx,ny,&
            v%mask(:,:,1),v%rdata(:,:,1),rmiss,lat0,lon0, &
            dlat,dlon, lon_inner, cgrib(ivar),lcgrib,lengrib, ierr)
       if(ierr/=0) then
20        format(A,I0)
          write(0,20) 'Could not encode grib record: ierr=',ierr
          lengribs(ivar)=0
       else
          lengribs(ivar)=lengrib
       endif
    enddo
    !$OMP END PARALLEL DO

    if(ierr/=0) goto 2000

    if(this%grib_pos<1) then
       write(0,*) 'ASSERTION FAILURE: GRIB_POS<1 (FILE NOT OPEN?)'
       stop 23
    endif

    writepos(1)=0 ! write position relative to current file pointer
    do ivar=2,num_vars
       if(lengribs(ivar-1)>0) then
          writepos(ivar)=writepos(ivar-1)+lengribs(ivar-1)+8
       else
          writepos(ivar)=writepos(ivar-1)
       endif
    enddo
    if(lengribs(num_vars)>0) then
       writesize=writepos(num_vars)+lengribs(num_vars)+8
    else
       writesize=writepos(num_vars)
    endif
    allocate(character(len=writesize) :: writeme)
    
    !$OMP PARALLEL DO PRIVATE(lengrib,i,ivar) DEFAULT(SHARED)
    do ivar=1,num_vars
       !write(0,*) 'copy var',ivar
       lengrib=lengribs(ivar)
       if(lengrib>0) then
          i=writepos(ivar)+1

          ! Store fortran record size in big endian
          writeme(i:i)=char(iand(255,ishft(lengrib,-24)))
          writeme(i+1:i+1)=char(iand(255,ishft(lengrib,-16)))
          writeme(i+2:i+2)=char(iand(255,ishft(lengrib,-8)))
          writeme(i+3:i+3)=char(iand(255,lengrib))

          ! Store grib record
          writeme(i+4:i+4+lengrib-1)=cgrib(ivar)(1:lengrib)

          ! Store fortran record size again (big endian)
          writeme(i+4+lengrib:i+4+lengrib+3)=writeme(i:i+3)
       endif
    enddo
    !$OMP END PARALLEL DO

19  format('Write ',I0,' bytes to grib file at pos=',I0,'.')
    write(0,19) len(writeme),this%grib_pos
    write(this%grib_unit,pos=this%grib_pos) writeme
29  format('Wrote ',I0,' bytes to grib file at pos=',I0,'.  Hurrah!')
    write(0,29) len(writeme),this%grib_pos
    this%grib_pos=this%grib_pos+len(writeme)

2000 continue ! cleanup section

    if(allocated(writeme))   deallocate(writeme)
    if(allocated(cgrib))     deallocate(cgrib)
  end subroutine write_grib2

  subroutine impl_make_grib2(atime, ftime, window, ivar, nx, ny, &
       mask, data, rmiss, lat0, lon0, dlat, dlon, lon_inner, &
       cgrib, lcgrib, lengrib, ierr)
    use datetime_module
    implicit none

    interface
      subroutine addfield(cgrib,lcgrib,ipdsnum,ipdstmpl,ipdstmplen, &
                          coordlist,numcoord,idrsnum,idrstmpl,      &
                          idrstmplen,fld,ngrdpts,ibmap,bmap,ierr)
        implicit none
        character(len=1),intent(inout) :: cgrib(lcgrib)
        integer,intent(in) :: ipdsnum,ipdstmpl(*)
        integer,intent(in) :: idrsnum,numcoord,ipdstmplen,idrstmplen
        integer,intent(in) :: lcgrib,ngrdpts,ibmap
        real,intent(in) :: coordlist(numcoord)
        real,target,intent(in) :: fld(ngrdpts)
        integer,intent(out) :: ierr
        integer,intent(inout) :: idrstmpl(*)
        logical*1,intent(in) :: bmap(ngrdpts)
      end subroutine addfield
    end interface
    real(kind=C_DOUBLE), intent(in) :: atime,ftime, window
    logical, intent(in) :: mask(nx,ny), lon_inner
    real, intent(in) :: lat0, lon0, dlat, dlon, rmiss
    real, intent(inout) :: data(nx,ny)
    integer, intent(inout) :: ierr,lengrib
    integer,intent(in) :: lcgrib, ivar, nx,ny
    character(len=lcgrib), intent(inout) :: cgrib

    real, allocatable :: flipdata(:,:)
    logical(kind=1), allocatable :: flipmask(:,:)
    integer, allocatable :: ideflist(:)

    real(kind=C_DOUBLE) :: stime, faminutes, saminutes
    integer(kind=C_INT64_T) :: ayear, amon, aday, ahr, amin
    integer(kind=C_INT64_T) :: fyear, fmon, fday, fhr, fmin
    integer(kind=C_INT64_T) :: syear, smon, sday, shr, smin
    integer(kind=C_INT64_T) :: fahr, famin, sahr, samin
    real :: rdummy(2)
    integer :: listsec0(2),listsec1(20),ipdst(200),ipdsnum
    integer :: i,j,igds(200),igdst(200),ilat0,ilon0,idlat, idlon
    integer :: iscan, ilat1,ilon1, idrsnum,idrs(200)
    integer :: idefnum,ibmap,nmask
    real(kind=8) :: sum

    allocate(flipmask(nx,ny))
    allocate(flipdata(nx,ny))
    allocate(ideflist(nx*ny+100))
 
    stime=ftime-window
    call epoch_split_ymdhm(atime,ayear, amon, aday, ahr, amin)
    call epoch_split_ymdhm(ftime,fyear, fmon, fday, fhr, fmin)
    call epoch_split_ymdhm(stime,syear, smon, sday, shr, smin)

    faminutes=nint((ftime-atime)/60.0)
    fahr=sign(abs(faminutes)/60.0,faminutes)
    famin=sign(mod(abs(faminutes),60.0),faminutes)

    saminutes=nint((stime-atime)/60.0)
    sahr=sign(abs(saminutes)/60.0,saminutes)
    samin=sign(mod(abs(saminutes),60.0),saminutes)

10  format('stime=',F0.1,' fcst: n=',F0.1,' h:mm=',I0,':',I2,' strt: n=',F0.1,&
           ' h:mm=',I0,':',I2)
    !print 10,stime,faminutes,fahr,famin,saminutes,sahr,samin

    ilat0=nint(lat0*1000000)
    ilon0=nint(mod(3600.0+lon0,360.0)*1000000)

    if(lon_inner) then
       ilat1=nint((lat0+(ny-1)*dlat)*1000000)
       ilon1=nint(mod(3600.0+lon0+(nx-1)*dlon,360.0)*1000000)
    else
       ilat1=nint((lat0+(nx-1)*dlat)*1000000)
       ilon1=nint(mod(3600.0+lon0+(ny-1)*dlon,360.0)*1000000)
    endif

    idlon=abs(nint(dlon*1000000))
    idlat=abs(nint(dlat*1000000))

    iscan=0
    if(dlon<0) iscan=iand(iscan,128)
    if(dlat>0) iscan=iand(iscan,64)
    if(.not.lon_inner) iscan=iand(iscan,32)
    
    nmask=0
    sum=0
    do j=1,ny
       do i=1,nx
          flipmask(i,j)=mask(i,ny-j+1)
          if(flipmask(i,j)) then
             nmask=nmask+1
             flipdata(i,j)=data(i,ny-j+1)
             sum=sum+flipdata(i,j)
          else
             flipdata(i,j)=rmiss
          endif
       enddo
    enddo

20  format('var ',A,': mean=',F0.3,', lon=',I,' ... ',I0,' by ',I0,', lat=',I0,' ... ',I0,' by ',I0,', iscan=',I0,', nmask=',I0)
    if(nmask>0) then
       print 20,trim(swath_wrfname(ivar)),sum/nmask,ilon0,ilon1,idlon, ilat0,ilat1,idlat, iscan, nmask
    else
       print 20,trim(swath_wrfname(ivar)),0.,ilon0,ilon1,idlon, ilat0,ilat1,idlat, iscan, nmask
    endif
    
    listsec0=0
    select case(ivar)
    case(ix_wind, ix_suswind, ix_updhel25, ix_updhel03)
       listsec0(1)=0         ! Discipline: meterology
    end select
    listsec0(2)=2         ! GRIB version 2

    listsec1=0
    listsec1(1)=7         ! Id of orginating centre (Common Code Table C-1)
    listsec1(2)=0         ! Id of orginating sub-centre (local
                          ! table)/Table C of ON388
    listsec1(3)=8         ! GRIB Master Tables Version Number (Code Table 1.0)
    listsec1(4)=1         ! GRIB Local Tables Version Number (Code Table 1.1)
    listsec1(5)=1         ! Significance of Reference Time (Code Table 1.2)
    listsec1(6)=ayear     ! Reference Time - Year (4 digits)
    listsec1(7)=amon      ! Reference Time - Month
    listsec1(8)=aday      ! End day
    listsec1(9)=ahr       ! End hour
    listsec1(10)=amin     ! Reference Time - Minute
    listsec1(11)=0        ! Reference Time - Second
    listsec1(12)=0        ! Production status of data (Code Table 1.3)
    listsec1(13)=1        ! Type of processed data (Code Table 1.4)

    ipdst=0
    ipdst(3)=2            ! Type of generating process = forecast
    ipdst(4)=255          ! Analysis or forecast generating process
                          ! identified (255=unknown/missing)
    ! 4 & 5 are generating process; HWRF doesn't have one
    ! 6 & 7 are related to observational data; unused
    ipdst(8)=0            ! Significance of time range: minutes

    ! Table 4.8 = Average, Accumulation and/or Extreme values or other
    ! Statistically-processed values at a horizontal level or in a
    ! horizontal layer in a continuous or non-continuous time interval
    ipdsnum=8

    ipdst(21)=0           ! End second
    ipdst(22)=1           ! Number of time ranges
    ipdst(23)=0           ! Number of missing values
    ipdst(26)=0           ! Unit of time range (minutes)
    !write(0,*) 'end',(endhr*60+endmin)-(starthr*60+startmin) ! Time range
    ipdst(28)=0           ! Unit of time range again (minutes)

    select case(ivar)
    case(ix_maxw,ix_minw)
       ipdst(1)=2             ! Category: momentum
       ipdst(2)=9             ! Vertical velocity (geometric) m/s
       ipdst(8)=0             ! Reference time units: minutes
       ipdst(9)=0             ! Reference time is analysis time
       ipdst(10)=108          ! 1st Level type: pressure relative to ground
       ipdst(11)=-4           ! 1st Level scale
       ipdst(12)=0            ! 1st Level value
       ipdst(13)=108          ! 2nd Level type: pressure relative to ground
       ipdst(14)=-4           ! 2nd Level scale
       ipdst(15)=4            ! 2nd Level value
       ipdst(16)=fyear        ! End year
       ipdst(17)=fmon         ! End month
       ipdst(18)=fday         ! End Time - Day
       ipdst(19)=fhr          ! End Time - Hour
       ipdst(20)=fmin         ! End minute
       ipdst(21)=0            ! End second
       if(ivar==ix_maxw) then
          ipdst(24)=2            ! Statistical method: maximum
       else
          ipdst(24)=3            ! Statistical method: minimum
       endif
       ipdst(25)=2            ! Successive times processed have same
                              ! start time of forecast, forecast time
                              ! is incremented.
       ipdst(26)=0            ! Time range units (hours)
       ipdst(27)=nint(faminutes) ! Time range
    case(ix_precip,ix_cuprecip)
       ipdst(1)=1             ! Category: moisture
       if(ivar==ix_precip) then
          ipdst(2)=7          ! Precipitation rate
       else
          ipdst(2)=37         ! Convective precipitation rate
       endif
       ipdst(8)=0             ! Reference time units: minutes
       ipdst(9)=0             ! Reference time is analysis time
       ipdst(10)=1            ! First level type: ground
       ipdst(13)=8            ! Second level type: top of atmosphere
       ipdst(16)=fyear        ! End year
       ipdst(17)=fmon         ! End month
       ipdst(18)=fday         ! End Time - Day
       ipdst(19)=fhr          ! End Time - Hour
       ipdst(20)=fmin         ! End minute
       ipdst(21)=0            ! End second
       ipdst(24)=1            ! Statistical method: accumulation
       ipdst(25)=2            ! Successive times processed have same
                              ! start time of forecast, forecast time
                              ! is incremented.
       ipdst(26)=0            ! Time range units (hours)
       ipdst(27)=nint(faminutes) ! Time range
    case(ix_wind,ix_suswind)
       ipdst(1)=2             ! Category: momentum
       if(ivar==ix_suswind) then
          ipdst(2)=1          ! Wind Speed
       else
          ipdst(2)=22         ! Wind Gust
       endif
       ipdst(8)=0             ! Reference time units: minutes
       ipdst(9)=0             ! Reference time is analysis time
       ipdst(10)=103          ! First level type: specified height above ground
       ipdst(11)=0            ! First level scale
       ipdst(12)=10           ! First level value
       ipdst(13)=103          ! Second level type: specified height above ground
       ipdst(14)=0            ! Second level scale
       ipdst(15)=10           ! Second level value
       ipdst(16)=fyear        ! End year
       ipdst(17)=fmon         ! End month
       ipdst(18)=fday         ! End Time - Day
       ipdst(19)=fhr          ! End Time - Hour
       ipdst(20)=fmin         ! End minute
       ipdst(21)=0            ! End second
       if(ivar==ix_suswind) then
          ipdst(22)=2
       endif
       ipdst(24)=2            ! Statistical process: maximum
       ipdst(25)=2            ! Successive times processed have same
                              ! start time of forecast, forecast time
                              ! is incremented.
       ipdst(26)=0            ! Time range units (hours)
       ipdst(27)=nint(faminutes) ! Time range
       if(ivar==ix_suswind) then
          !ipdst(28)=0         ! Unit = minutes
          !ipdst(29)=1         ! one minute
          ipdst(30)=3         ! Minimum
          ipdst(31)=2         ! Successive times processed have same
                              ! start time of forecast, forecast time
                              ! is incremented.
          ipdst(32)=0         ! Unit = minutes
          ipdst(33)=1         ! one minute
       endif
    case(ix_updhel25, ix_updhel03)
       ipdst(1)=7             ! Category: thermodynamic stability indices
       ipdst(2)=15            ! Updraft helicity
       ipdst(8)=0             ! Reference time units (minutes)
       ipdst(9)=nint(saminutes) ! Reference time
       ipdst(10)=103          ! First level type: specified height above ground
       ipdst(11)=-3           ! First level scale
       ipdst(13)=103          ! Second level type: specified height above ground
       ipdst(14)=-3           ! Second level scale

       if(ivar==ix_updhel03) then
          ipdst(12)=0         ! First level value
          ipdst(15)=3         ! Second level value
       else
          ipdst(12)=2         ! First level value
          ipdst(15)=5         ! Second level value
       endif
       ipdst(16)=syear        ! End year
       ipdst(17)=smon         ! End month
       ipdst(18)=sday         ! End Time - Day
       ipdst(19)=shr          ! End Time - Hour
       ipdst(20)=smin         ! End minute
       ipdst(24)=2            ! Statistical process: maximum
       ipdst(25)=2            ! Analysis time kept fixed, sliding
                              ! forecast lead time
       ipdst(26)=0            ! Time range units (hours)
       ipdst(27)=nint(window/60) ! Time range
    end select

    igds=0
    igds(1)=0             ! Grid definition is specified in template
    igds(2)=nx*ny         ! Number of grid points
    igds(3)=0             ! Regular grid
    igds(4)=0             ! Regular grid
    igds(5)=0             ! Lat-lon grid

    igdst=0
    ! igdst 1:7 deal with shape of earth
    if(lon_inner) then
       igdst(8)=nx           ! Number of points in parallel
       igdst(9)=ny           ! Number of points in meridian
    else
       igdst(8)=ny           ! Number of points in parallel
       igdst(9)=nx           ! Number of points in meridian
    endif
    ! igdst(10:11) are to change angle of projection
    igdst(12)=ilat1       ! Last latitude in microdegrees
    igdst(13)=ilon0       ! First longitude in microdegres
    igdst(14)=48          ! 48 = I & J direction increments given
    igdst(15)=ilat0       ! First latitude in microdegrees
    igdst(16)=ilon1       ! Last longitude in microdegres
    igdst(17)=idlon       ! I direction increment
    igdst(18)=idlat       ! J direction increment
    igdst(19)=iscan       ! Scanning mode flags

    idrs=0
    if(nmask<1) then
       idrsnum=0             ! Simple packing
    else
       idrsnum=3             ! Complex encoding
       idrs(6)=1             ! General group split
       idrs(7)=1             ! Missing value management
       call mkieee(rmiss,idrs(8),1) ! Missing value
       idrs(17)=2            ! Order of secondary differences (1 or 2)
    endif
    select case(ivar)
    case(ix_wind,ix_suswind,ix_maxw,ix_minw)
       idrs(3)=3          ! Decimal scaling factor = 1000
    case(ix_precip,ix_cuprecip)
       idrs(3)=5          ! Decimal scaling factor = 100000
    case(ix_updhel03,ix_updhel25)
       idrs(3)=2          ! Decimal scaling factor = 100
    end select

    call gribcreate(cgrib,lcgrib,listsec0,listsec1,ierr)
    if(ierr/=0) then
       write(0,300) 'error creating new GRIB2 record',ierr
       goto 2000
    endif

    call addgrid(cgrib,lcgrib,igds,igdst,200,ideflist,idefnum,ierr)
    if(ierr/=0) then
       write(0,300) 'error adding grid to GRIB2 record',ierr
       goto 2000
    endif

    ibmap=0
    call addfield(cgrib,lcgrib,ipdsnum,ipdst,200,&
                  rdummy,0,idrsnum,idrs,200, &
                  flipdata,nx*ny,ibmap,flipmask,ierr)
    if(ierr/=0) then
       write(0,300) 'error adding field to GRIB2 record',ierr
       goto 2000
    endif

    call gribend(cgrib,lcgrib,lengrib,ierr)
    if(ierr/=0) then
       write(0,300) 'error ending GRIB2 record',ierr
       goto 2000
    endif 

300 format(A,': ',I0)

    2000 continue ! cleanup

    deallocate(flipmask)
    deallocate(flipdata)
    deallocate(ideflist)
  end subroutine impl_make_grib2
  
  subroutine write(this,outprefix,track,ntrack,grads_byteswap)
    use tcf_module, only: atcf
    use sysutil_module, only: get_unit
    use datetime_module, only: grads_date, grads_date_len
    ! This subroutine writes out the various swath output files.  This
    ! is translated from Bob Tuleya's original source code.  It now
    ! uses the internal HWRF tracker, and the swathgen output swath
    ! data.  The external tracker's ATCF track is only used to decide
    ! the cutoff hour, beyond which to produce no swath.
    class(swathgen), intent(inout) :: this
    type(vardata_real), pointer :: wind, suswind, precip
    character*(*), intent(in) :: outprefix
    type(atcf), intent(in) :: track(ntrack)
    integer, intent(in) :: ntrack
    logical, intent(in) :: grads_byteswap
    integer :: iwindhrly,iwind10m,irainfall,isuswind10m, i,j, iswathdat,igrads
    integer :: lastitime,itime,ihour, nx,ny, ifile, ivar
    real :: lat,lon,tmns, hour, lasthour
    character(len=grads_date_len) :: gr2date
    
    if(ntrack<1) return ! no ATCF, so no storm
    lasthour=track(ntrack)%fcsthour*this%nauxphr
    if(lasthour<1e-5) return ! no forecast, so no swath

    call this%get_track() ! find hourly track from internal HWRF tracker

    ! Aliases for convenience:
    nx=this%nx
    ny=this%ny
    wind => this%datavars(ix_wind)
    precip => this%datavars(ix_precip)
    suswind => this%datavars(ix_suswind)

    ! Open all output files:
    iwindhrly=get_unit()
    open(iwindhrly,file=trim(outprefix)//'wind10hrly.ascii',form='FORMATTED',status='UNKNOWN')

    iwind10m=get_unit()
    open(iwind10m,file=trim(outprefix)//'wind10m.ascii',form='FORMATTED',status='UNKNOWN')

    isuswind10m=get_unit()
    open(isuswind10m,file=trim(outprefix)//'suswind10m.ascii',form='FORMATTED',status='UNKNOWN')

    irainfall=get_unit()
    open(irainfall,file=trim(outprefix)//'rainfall.ascii',form='FORMATTED',status='UNKNOWN')

    iswathdat=get_unit()
    open(iswathdat,file=trim(outprefix)//'swath.dat',form='unformatted', &
         ACCESS='DIRECT',RECL=(this%nx*this%ny*4),STATUS='UNKNOWN')

    igrads=get_unit()
    open(igrads,file=trim(outprefix)//'swath.ctl',form='FORMATTED', &
         status='UNKNOWN')

    ! Fill with badval outside the mask:
    do ivar=1,num_vars
       !$OMP PARALLEL DO PRIVATE(i,j)
       badvalloop: do j=1,this%ny
          do i=1,this%nx
             if(.not.this%datavars(ivar)%mask(i,j,1)) then
                this%datavars(ivar)%rdata(i,j,1)=-999.0
             endif
          enddo
       enddo badvalloop
    enddo
    lastitime=this%ntimes

    itimeloop: do itime=2,this%ntimes
       hour=(itime-1) / this%nauxphr
       if(hour>lasthour) then
          lastitime=itime-1
          exit itimeloop
       endif
    enddo itimeloop
    
    ! Original comment from Bob Tuleya:
    !     now write data out to text file for NHC
    !     make the ascii files for Michelle's color plots
3333 format(3f9.2)
3332 format(5f9.2,2i5)
2900 format(1x,'HOUR: ', f5.1,5x, 'MAX SURF WIND (KNOTS): ', f5.1)
    ! Write the lats & lons to the wind10m and rainfall files for the
    ! color plots.  Note no OpenMP here since it is a write loop.

    tmns = -99.0

    !$OMP PARALLEL DO PRIVATE(i,j,ifile,itime,hour) DEFAULT(SHARED) SCHEDULE(STATIC)
    ifile_loop: do ifile=1,9
303    format('TXT/BIN output file #',I0,'/',I0)
       write(0,303) ifile,9
       select case(ifile)
       case(1)
          ! Print the max wind for forecast hours, and find the last index
          ! with a forecast time that is at or before the last ATCF time.
          windloop2: do itime=2,this%ntimes
             hour=(itime-1) / this%nauxphr
             if(hour>lasthour) then
                exit windloop2
             endif
             write(iwindhrly,2900) hour, this%track_wind(itime)
          enddo windloop2
       case(2)
          write(iwind10m,3332) this%proj%lon1, this%proj%lon1+this%proj%dlon*(nx-1), &
               this%proj%lat1, this%proj%lat1+this%proj%dlat*(ny-1), &
               (this%proj%dlon+this%proj%dlat)/2., this%nx,this%ny
          do  j = 1 , ny
             do  i = 1 , nx
                write(iwind10m,3333) this%lats%rdata(i,j,1),this%lons%rdata(i,j,1), &
                     wind%rdata(i,j,1)
             enddo
          enddo
          write(iwind10m,3333)tmns,tmns,tmns
          write(iwind10m,3335)nint(lasthour)
          do i = 2 , lastitime
             write(iwind10m,3334) this%track_lat(i),this%track_lon(i)
          enddo
       case(3)
          write(isuswind10m,3332) this%proj%lon1, this%proj%lon1+this%proj%dlon*(nx-1), &
               this%proj%lat1, this%proj%lat1+this%proj%dlat*(ny-1), &
               (this%proj%dlon+this%proj%dlat)/2., this%nx,this%ny
          do  j = 1 , ny
             do  i = 1 , nx
                write(isuswind10m,3333) this%lats%rdata(i,j,1),this%lons%rdata(i,j,1), &
                     suswind%rdata(i,j,1)
             enddo
          enddo
          write(isuswind10m,3333)tmns,tmns,tmns
          write(isuswind10m,3335)nint(lasthour)
          do i = 2 , lastitime
             write(isuswind10m,3334) this%track_lat(i),this%track_lon(i)
          enddo
       case(4)
          write(irainfall,3332) this%proj%lon1, this%proj%lon1+this%proj%dlon*(nx-1), &
               this%proj%lat1, this%proj%lat1+this%proj%dlat*(ny-1), &
               (this%proj%dlon+this%proj%dlat)/2., this%nx,this%ny
          do  j = 1 , ny
             do  i = 1 , nx
                write(irainfall,3333) this%lats%rdata(i,j,1),this%lons%rdata(i,j,1), &
                     precip%rdata(i,j,1)
             enddo
          enddo
          write(irainfall,3333)tmns,tmns,tmns
          write(irainfall,3335)nint(lasthour)
          do i = 2 , lastitime
             write(irainfall,3334) this%track_lat(i),this%track_lon(i)
          enddo
       case(5)
          ! Output the swath.dat binary file used for plotting swath data
          write(iswathdat,rec=1) wind%rdata
          write(iswathdat,rec=2) precip%rdata
          write(iswathdat,rec=3) suswind%rdata
       case(6)
          ! Output a GrADS control file:
          write(igrads,*)'DSET ^'//trim(outprefix)//'swath.dat'
          write(igrads,*)'TITLE 10m max wind and rain'
          write(igrads,*)'UNDEF -999'
          if(grads_byteswap) then
             write(igrads,*)'OPTIONS byteswapped'
          endif
          write(igrads,9200)this%nx,this%proj%lon1,this%proj%dlon
9200      format(' XDEF ',I5,' LINEAR ',F7.2,1x,F6.3)
          write(igrads,9201)this%ny,this%proj%lat1,this%proj%dlat
9201      format(' YDEF ',I5,' LINEAR ',F7.2,1x,F6.3)
          write(igrads,*)'ZDEF 1 LEVELS 0'
          call grads_date(gr2date,this%files(1)%epochtimes(1))
9202      format(' ',A,' ',A,' ',A)
          write(igrads,9202)'tdef 1 linear',gr2date,'1mon'
          write(igrads,*)'VARS 3'
          if(this%in_crazy_us_units) then
             write(igrads,*)'gustwind10 0 99 10-m max gust wind (kt)'
             write(igrads,*)'rain       0 99 accumulated precipitation(IN)'
             write(igrads,*)'suswind10  0 99 10-m max sustained wind (kt)'
          else
             write(igrads,*)'gustwind10 0 99 10-m max gust wind (m/s)'
             write(igrads,*)'rain       0 99 accumulated precipitation(m)'
             write(igrads,*)'suswind10  0 99 10-m max sustained wind (m/s)'
          endif
          write(igrads,*)'ENDVARS'
       case(7)
          call make_pgm('debug-wind.pgm',this%nx,this%ny,wind%rdata,wind%mask,0.,4.)
       case(8)
          call make_pgm('debug-suswind.pgm',this%nx,this%ny,suswind%rdata,suswind%mask,0.,4.)
       case(9)
          call make_pgm('debug-precip.pgm',this%nx,this%ny,precip%rdata,precip%mask,0.,200.)
       end select
    enddo ifile_loop
    ! NOT !$OMP END PARALLEL DO

3335 format(i5)
3334 format(2f8.2)

    close(iwindhrly)
    close(iwind10m)
    close(isuswind10m)
    close(irainfall)
    close(iswathdat)
    close(igrads)

  end subroutine write

  subroutine make_pgm(filename,nx,ny,var,mask,sub,scale)
    implicit none
    character*(*), intent(in) :: filename
    integer(kind=1) :: bytes(nx,ny)
    real, intent(in) :: var(nx,ny), sub,scale
    integer, intent(in) :: nx,ny
    logical, intent(in) :: mask(nx,ny)
    integer :: i,j,unit
    ! NOTE: Don't parallelize this since it is already within an OpenMP section.
    do j=1,ny
       do i=1,nx
          if(mask(i,j)) then
             bytes(i,j) = max(0,min(199,nint((var(i,j)-sub) * scale))) + 50
          else
             bytes(i,j) = 0
          endif
       enddo
    enddo

33  format('P5',/,I0,' ',I0,' 255')

    open(newunit=unit,file=trim(filename),status='UNKNOWN',form='FORMATTED')
    write(unit,33) nx,ny
    close(unit)
    open(newunit=unit,file=trim(filename),status='OLD',access='STREAM',position='APPEND')
    write(unit) bytes
    close(unit)
  end subroutine make_pgm

  subroutine free_swathgen(this)
    class(swathgen), intent(inout) :: this
    integer :: i,j,iostat
    call this%discard_track()
    if(associated(this%files)) then
       do i=1,min(this%nfiles,size(this%files))
          call this%files(i)%close()
       enddo
       deallocate(this%files)
       nullify(this%files)
    endif
    if(associated(this%interp)) then
       do i=1,min(this%nfiles,size(this%interp))
          call this%interp(i)%free()
       enddo
       deallocate(this%interp)
       nullify(this%interp)
    endif
    if(associated(this%filevars)) then
       do j=1,size(this%filevars,2)
          do i=1,size(this%filevars,1)
             call this%filevars(i,j)%free()
          enddo
       enddo
       deallocate(this%filevars)
       nullify(this%filevars)
    endif
    if(this%have_proj) then
       do i=1,size(this%datavars)
          call this%datavars(i)%free()
       enddo
       deallocate(this%datavars)
       nullify(this%datavars)
       call this%proj%free()
       call this%dec%free()
       call this%lats%free()
       call this%lons%free()
    endif
    this%maxnfile=0
    this%nfiles=0
    this%itime=0
    this%ntimes=0
    this%nx=0
    this%ny=0
    this%have_proj=.false.
    this%inited=.false.
    this%did_first_interp=.false.
    this%have_track=.false.
    if(this%write_grib .and. this%grib_unit/=0) then
       close(this%grib_unit)
    endif
    this%write_grib=.false.
    this%grib_unit=0
    this%grib_pos=1
  end subroutine free_swathgen

  subroutine init_swathgen(this,maxnfile)
    integer, intent(in) :: maxnfile
    class(swathgen), intent(inout) :: this
    integer :: i,ios

    allocate(this%files(maxnfile))
    allocate(this%datavars(num_vars))
    allocate(this%filevars(maxnfile,num_vars))
    do i=1,maxnfile
       call init_wrfdiag_file(this%files(i))
    enddo
    this%nfiles=0
    this%maxnfile=maxnfile
    this%in_crazy_us_units=.false.
    call init_vardata_real(this%lats)
    call init_vardata_real(this%lons)
  end subroutine init_swathgen

  subroutine add_file(this,filename)
    class(swathgen), intent(inout) :: this
    character*(*), intent(in) :: filename
    class(wrfdiag_file), pointer :: f
    integer :: ifile
    character*255 message
    ifile=this%nfiles+1

    if(ifile>this%maxnfile) then
3011   format(A,': too many files.  Could not add this one.  Increase maxnfile in call to init_swathgen.  You provided maxnfile=',I0) ! ,/,' -- PROGRAM IS ABORTING DUE TO INSUFFICIENT MAXNFILE DIMENSION IN CALL TO INIT_SWATHGEN --')
       write(message,3011) trim(filename),this%maxnfile
       call fail(message)
    endif

    this%in_crazy_us_units=.false.
    this%nfiles=ifile
    f=>this%files(ifile)
    call init_wrfdiag_file(f,trim(filename))
    !call f%make_latlonrot(2) ! 2 = first non-analysis time
    if(ifile==1) call this%first_file(f)
    if(f%ntimes /= this%ntimes) then
       write(message,30) trim(f%path),f%ntimes,trim(this%files(1)%path),this%files(1)%ntimes
       call fail(message)
30     format(A,': ntimes mismatch.  File has ',I0,' times, while first file (',A,') has ',I0,'.',/,'  -- PROGRAM IS ABORTING DUE TO MISMATCH IN TIME COUNT --')
    endif
  end subroutine add_file

  subroutine first_file(this,f)
    ! Called when the first file is read in.  Only argument is the
    ! first file's wrfdiag_file object.  
    class(swathgen), intent(inout) :: this
    class(wrfdiag_file), pointer, intent(inout) :: f
    this%ntimes = f%ntimes
    this%nauxphr = 3600./(f%epochtimes(2)-f%epochtimes(1))
  end subroutine first_file

  subroutine setup_proj(this, latres, lonres, latpad1, latpad2, lonpad1, lonpad2, first_time, last_time)
    class(swathgen), intent(inout) :: this
    class(wrfdiag_file), pointer :: f
    integer, optional, intent(in) :: first_time, last_time
    real, intent(in) :: latres,lonres, lonpad1,lonpad2, latpad1,latpad2
    integer :: ifile, nlat,nlon, ilast
    real :: lat1f,lon1f,lat2f,lon2f,latmidf,lonmidf
    real :: lat1,lon1,lat2,lon2,latmid,lonmid
    real(kind=C_DOUBLE) :: hourdiff, hourdiff0
    character*256 :: message

    print '(A)','Setup swath projection...'
    if(this%nfiles<1) then
       call fail('In swathgen, no files provided before call to setup_proj.')
    endif

    lat1=1e19 ; lat2=-1e19
    lon1=1e19 ; lon2=-1e19
    do ifile=1,this%nfiles
       f=>this%files(ifile)
       call f%latlon_bounds(lat1f,latmidf,lat2f, lon1f,lonmidf,lon2f, &
            first_time=first_time, last_time=last_time)
       if(lat1f<lat1) lat1=lat1f
       if(lon1f<lon1) lon1=lon1f
       if(lat2f>lat2) lat2=lat2f
       if(lon2f>lon2) lon2=lon2f
       hourdiff=(f%epochtimes(2)-f%epochtimes(1))/3600.
       if(ifile==1) then
          hourdiff0=hourdiff
       elseif(hourdiff0/=hourdiff) then
          write(message,3301) hourdiff0,hourdiff
3301      format('Error: mismatch in output frequencies between files (',F0.3,' hours vs. ',F0.3,' hours).  All files must have the same output frequency.')
       endif
    enddo
    if(hourdiff0<1/3600.) then
       call fail('Likely an error: wrfdiag output frequency is less than a second.  Aborting.')
    endif
    this%nauxphr=1./hourdiff0

    lat1 = latres *  floor (lat1/latres) - latpad1
    lon1 = lonres *  floor (lon1/lonres) - lonpad1
    lat2 = latres * ceiling(lat2/latres) + latpad2
    lon2 = lonres * ceiling(lon2/lonres) + lonpad2
    nlat = floor((lat2-lat1-latres/2)/latres)+1
    nlon = floor(mod(3600+lon2-lon1-lonres/2,360.)/lonres)+1
    this%ny=nlat
    this%nx=nlon

    call init_proj_latlon(this%proj,lat1,lon1,latres,lonres,nlon,nlat, &
         lon_inner=.true.,cyclic_lon=.false.,north_polar=.false.,      &
         south_polar=.false.)

    call init_decomp(this%dec,          &
                     1,nlon,1,nlat,1,1, &
                     1,nlon,1,nlat,1,1  )

104 format('Generate projection with nx=',I0,' ny=',I0,' lat1=',F0.3,' lon1=',F0.3,' dlat=',F0.3,' dlon=',F0.3)
    print 104, nlon,nlat,lat1,lon1,latres,lonres

    this%in_crazy_us_units=.false.
    this%have_proj=.true.
  end subroutine setup_proj

  subroutine interp_time(this,itime)
    class(swathgen), target, intent(inout) :: this
    integer, intent(in) :: itime
    integer :: ifile, i,j,ivar
    class(decomp), pointer :: fd, sd ! file and swath dimensions
    class(wrfdiag_file), pointer :: f
    class(wrfdiag_var), pointer :: w1,s1,p1 ! gust, sust. wind, and precip input
    real :: w,p,s
    logical :: m
    character*255 message

    if(this%nfiles<1) then
       call fail('In swathgen, no files provided before call to interp_next_time.')
    endif

    if(itime>this%ntimes .or. itime<2) then
30     format('Error: interp_time called with invalid itime.  The itime argument must be between 2 and ntimes (',I0,'), inclusive.  You provided ',I0,'.')
       write(message,30) this%ntimes,itime
       call fail(message)
    endif
    
    this%itime=itime
    sd=>this%dec
    f=>this%files(1)

    ! Do the initial stuff at the first time:
    if(.not. this%did_first_interp) then
       print '(A)','Do first interp initialization...'
       if(.not. this%have_proj) &
            call this%setup_proj(0.1,0.1,0.5,0.5,0.5,0.5)

       ! Allocate data and initialize to 0:
       do ivar=1,num_vars
          ! Allocate space for data:
          call init_vardata_real(this%datavars(ivar),this%proj,this%dec,&
               mask=.false.,alloc=.true.)

          ! Fill new arrays with 0 using appropriate thread:
          !$OMP PARALLEL DO PRIVATE(i,j)
          do j=sd%jps,sd%jpe
             do i=sd%ips,sd%ipe
                this%datavars(ivar)%rdata(i,j,1)=0
             enddo
          enddo

          ! Allocate mask and fill with "no data here" for all locations:
          call this%datavars(ivar)%alloc_mask(fill=.false.)
       enddo

       call init_vardata_real(this%lats,this%proj,this%dec,mask=.false.,alloc=.true.)
       call init_vardata_real(this%lons,this%proj,this%dec,mask=.false.,alloc=.true.)

       ! Get lats & lons:
       call this%proj%projinfo(1,this%nx,1,this%ny,1,1, &
                               1,this%nx,1,this%ny,1,1, &
                               1,this%nx,1,this%ny,1,1, &
                               lat=this%lats%rdata,     &
                               lon=this%lons%rdata      )

       ! Allocate the mask and set it to "no data here" for all locations:
       do ivar=1,num_vars
          call this%datavars(ivar)%alloc_mask(fill=.false.)
       enddo

       ! Allocate the arrays used to store per-file precip and wind:
       allocate(this%interp(this%nfiles))
       allocate(this%filevars(this%nfiles,num_vars))
       do ifile=1,this%nfiles
          do ivar=1,num_vars
             call init_wrfdiag_var(this%filevars(ifile,ivar))
          enddo
          call init_bilinear_real(this%interp(ifile))
       enddo
       this%did_first_interp=.true.
    else
       do ivar=1,num_vars
          if(swath_clear(ivar)) then
             write(0,*) 'clear ',swath_wrfname(ivar)
             !$OMP PARALLEL DO PRIVATE(i,j)
             do j=1,this%ny
                do i=1,this%nx
                   this%datavars(ivar)%mask(i,j,1)=.false.
                   this%datavars(ivar)%rdata(i,j,1)=0
                enddo
             enddo
          endif
       enddo
    endif

    ! Process all files in order received:
    fileloop: do ifile=1,this%nfiles
       f=>this%files(ifile)
       fd=>f%de_ij

333    format('Projection for file ',I0,' time ',I0)
       print 333, ifile,itime
       call f%eproj(itime)%print(unit=6)

       ! Make aliases for this file precip, wind and sus wind for convenience:
       p1=>this%filevars(ifile,ix_precip)
       w1=>this%filevars(ifile,ix_wind)
       s1=>this%filevars(ifile,ix_suswind)

       ! Read the variables, and add a mask if there isn't one already:
       write(0,*) 'read data for ifile=',ifile
       read_loop: do ivar=1,num_vars
          ! Read the next time's variables for this file:
          call init_wrfdiag_var(this%filevars(ifile,ivar),f,&
               trim(swath_wrfname(ivar)),itime)
          ! Allocate the mask if it isn't there already:
          call this%filevars(ifile,ivar)%alloc_mask(fill=.false.)
       enddo read_loop
       write(0,*) 'interpolate...'

       ! Initialize the interpolator:
       call this%interp(ifile)%free()
       call init_bilinear_real(this%interp(ifile),this%datavars(ix_wind),&
            this%filevars(ifile,ix_wind))
       call this%interp(ifile)%prep()

       ! Replace the square of the wind with the wind, and update
       ! the bitmask:
       !$OMP PARALLEL DO PRIVATE(i,j,w,p,m,s)
       jfixit: do j=fd%jps,fd%jpe
          ifixit: do i=fd%ips,fd%ipe
             w=max(0.,w1%rdata(i,j,1)) ! discard negative windsq
             s=max(0.,s1%rdata(i,j,1)) ! discard negative sustained wind mag
             p=max(0.,p1%rdata(i,j,1)) ! discard negative precip
             w=sqrt(w)                 ! wind is sqrt(windsq)

             ! Store back corrected values and non-squared winds:
             w1%rdata(i,j,1) = w
             s1%rdata(i,j,1) = s
             p1%rdata(i,j,1) = p

             ! Calculate the mask.
             m=(w>0 .or. p>0 .or. s>0) 
             w1%mask(i,j,1) = m
             s1%mask(i,j,1) = m
             p1%mask(i,j,1) = m
             ! In the mask generation, note the intentional bitwise
             ! zero comparison, with no epsilon.  This is because
             ! the only bitwise zeros in windsq should be places
             ! where there is no swath.  Near-zero wind is possible
             ! in areas where wind is very small due to dynamical
             ! reasons.  Wind that is bitwise zero is extremely
             ! unlikely, except where there is no swath data.
          enddo ifixit
       enddo jfixit

       ! Store masks for other vars:
       mask_other_vars: do ivar=1,num_vars
          if(ivar==ix_wind) cycle mask_other_vars
          if(ivar==ix_precip) cycle mask_other_vars
          if(ivar==ix_suswind) cycle mask_other_vars
          jcopymask: do j=fd%jps,fd%jpe
             icopymask: do i=fd%ips,fd%ipe
                this%filevars(ifile,ivar)%mask(i,j,1)=w1%mask(i,j,1)
             end do icopymask
          end do jcopymask
       end do mask_other_vars

       ! Interpolate wind and precip, storing the maximum of the
       ! current and prior data:
       interp_loop: do ivar=1,num_vars
          call this%interp(ifile)%scalar(this%datavars(ivar),&
               this%filevars(ifile,ivar), op=interp_op(ivar))
       enddo interp_loop
    enddo fileloop
  end subroutine interp_time
  
end module swathgen_module
