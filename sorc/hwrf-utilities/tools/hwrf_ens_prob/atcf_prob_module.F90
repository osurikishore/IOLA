module atcf_prob_module
  use, intrinsic :: iso_c_binding, only: C_DOUBLE
  use fileop_module, only: glob, init_glob
  use vardata_module, only: vardata_real, full_init_vardata_real
  use projection_module, only: proj_latlon, init_proj_latlon
  use decomp_module, only: decomp, init_decomp
  use ens_prob_module, only: insert_ens_domain, dist_weight, filenamelen, &
       max_output_times
  use ens_output_module, only: write_pgm
  implicit none

  private
  public :: member_atcf, model_atcf, multimodel

  integer, parameter, public :: maxmodel=10 ! maximum number of models

  integer, parameter, public :: ERR_ATCF_NO_VALID_LINES = 100, &
       ERR_ATCF_FIRST_TIME_NOT_0FHR=200, &
       ERR_ATCF_MULTIPLE_TRACKS=300, &
       ERR_ATCF_ATIME_MISMATCH=400

  integer, parameter, public :: minensmin=0, maxensmax=999 ! max members per ensemble
  ! NOTE: Different ens range than ens_prob_module for multi-model
  ! ensemble to allow bigger ensembles: 999 members instead of 99.

  integer, parameter, public :: namelen=4 !ATCF name length

  type member_atcf
     integer :: ntimes=0, alloced=0, last_time=0
     real, pointer :: lats(:)=>NULL(), lons(:)=>NULL(), times(:)=>NULL()
   contains
     procedure free => free_member
     procedure read_atcf => read_atcf_member
     procedure update_prob => update_prob_member
  end type member_atcf

  type model_atcf
     private
     integer :: ensmin=1,ensmax=0
     real :: stdev0=30e3, stdev_per_hr=120e3, cendist=10e3, maxdist=240e3
     logical, pointer :: ensflag(:)=>NULL()
     character(len=filenamelen) :: atcf_glob
     character(len=namelen) :: name,subname

     real, public :: weight=1.0
     type(member_atcf), pointer, public :: atcf(:)=>NULL()
     integer, public :: enscount=0
   contains
     procedure, public :: print => print_model
     procedure, public :: set_values
     procedure, public :: set_weight => set_weight_model
     procedure, public :: free => free_model
     procedure, public :: read_atcf => read_atcf_model
     procedure, public :: update_prob => update_prob_model
  end type model_atcf

  type multimodel
     real :: out_times(max_output_times), minprob
     logical :: sent(max_output_times)
     real :: ftime, time_step, last_ftime, latres, lonres
     real(kind=C_DOUBLE) :: atime ! analysis time
     character(len=filenamelen) :: &
          out_pgm_format, out_dat_format, out_ctl_format
     integer :: weight_method=1,nx=0,ny=0,notimes=0,itime=0
     logical :: modelflag(maxmodel)
     type(proj_latlon) :: proj
     type(decomp) :: dec
     type(vardata_real), pointer :: prob=>NULL(), work=>NULL()
     real, pointer :: lats(:)=>NULL(), lons(:)=>NULL()
     type(model_atcf) :: mod(maxmodel)
   contains
     procedure free => free_multimodel
     procedure read_namelist
     procedure read_atcf => read_atcf_multimodel
     procedure next_time
     procedure set_weight => set_weight_multimodel
     procedure time_for_output
     procedure update_prob => update_prob_multimodel
  end type multimodel

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine free_member(m)
    class(member_atcf) :: m
    if(associated(m%lats)) then
       deallocate(m%lats)
       nullify(m%lats)
    endif
    if(associated(m%lons)) then
       deallocate(m%lons)
       nullify(m%lons)
    endif
    if(associated(m%times)) then
       deallocate(m%times)
       nullify(m%times)
    endif
    m%ntimes=0
    m%alloced=0
  end subroutine free_member

  subroutine free_model(m)
    class(model_atcf) :: m
    integer :: ensid
    if(associated(m%atcf)) then
       do ensid=m%ensmin,m%ensmax
          call m%atcf(ensid)%free()
       enddo
       deallocate(m%atcf)
       nullify(m%atcf)
    endif

    if(associated(m%ensflag)) then
       deallocate(m%ensflag)
       nullify(m%ensflag)
    endif
  end subroutine free_model

  subroutine free_multimodel(m)
    class(multimodel) :: m
    integer :: i
    if(associated(m%prob)) then
       call m%prob%free()
       deallocate(m%prob)
       nullify(m%prob)
    endif
    if(associated(m%work)) then
       call m%work%free()
       deallocate(m%work)
       nullify(m%work)
    endif

    if(associated(m%lats)) then
       deallocate(m%lats)
       nullify(m%lats)
    endif

    if(associated(m%lons)) then
       deallocate(m%lons)
       nullify(m%lons)
    endif

    do i=1,maxmodel
       call m%mod(i)%free()
    enddo
  end subroutine free_multimodel

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine set_values(m,ensmin,ensmax,stdev0,stdev_per_hr, &
                        atcf_glob,name,subname,cendist,maxdist)
    class(model_atcf) :: m
    integer, intent(in) :: ensmin,ensmax
    real, intent(in) :: stdev0,stdev_per_hr,cendist,maxdist
    character(len=filenamelen) :: atcf_glob
    character(len=namelen) :: name,subname

    if(associated(m%atcf) .or. associated(m%ensflag)) then
       call m%free()
    endif

    if(m%ensmax>m%ensmin) then
       ! Should never get here.
303    format('ERROR: Invalid range ensid=[',I0,'..',I0,'] sent to set_values.')
       write(0,303) ensmin,ensmax
       write(6,303) ensmin,ensmax
       stop 3
    endif

    m%ensmin=ensmin
    m%ensmax=ensmax
    m%stdev0=stdev0
    m%stdev_per_hr=stdev_per_hr
    m%atcf_glob=atcf_glob
    m%name=name
    m%subname=subname
    m%cendist=cendist
    m%maxdist=maxdist

    allocate(m%atcf(ensmin:ensmax))
    allocate(m%ensflag(ensmin:ensmax))
    m%ensflag=.false.
  end subroutine set_values

  subroutine read_namelist(mm,filename,outname,ierr)
    class(multimodel) :: mm
    character*(*), intent(in) :: filename, outname
    integer, intent(inout) :: ierr
    ! Local copies of mm%stuff:
    character(len=filenamelen) :: atcf_glob(maxmodel), &
         out_pgm_format, out_dat_format, out_ctl_format
    integer :: ensmin(maxmodel), ensmax(maxmodel), weight_method
    real :: stdev0(maxmodel), stdev120(maxmodel), out_times(max_output_times)
    real :: time_step_hrs,latres,lonres, cendist(maxmodel), maxdist(maxmodel)
    real :: last_fhr
    character(len=namelen) :: name(maxmodel), subname(maxmodel)
    ! Work variables
    integer :: i, modelcount
    real :: stdev_per_hr
    logical :: good

    namelist/input/ atcf_glob, name, subname, ensmax, last_fhr, time_step_hrs
    namelist/proc/ stdev0, stdev120, weight_method, latres,lonres, &
                   cendist,maxdist
    namelist/output/ out_pgm_format, out_times

    ! Initialize namelist values to reasonable defaults.
    latres=0.1
    lonres=0.1
    out_pgm_format='multimodel_track_prob_f<ifhr>.pgm'
    out_dat_format='multimodel_track_prob_f<ifhr>.dat'
    out_ctl_format='multimodel_track_prob_f<ifhr>.ctl'
    ensmin=1
    ensmax=0
    out_times=-999.
    out_times(1:4)=(/ 72.0, 78.0, 120.0, 126.0 /)
    last_fhr=126.
    weight_method=1
    do i=1,maxmodel
       atcf_glob(i)='#'
       name(i)='#'
       subname(i)='#'
       stdev0(i)=30e3
       stdev120(i)=120e3
    end do
    time_step_hrs=6.0

    ! Read the namelist
333 format(A,': error reading namelist &',A,': iostat=',I0)
    open(unit=50,file=trim(filename),form='FORMATTED',status='OLD')

    read(50,nml=proc,iostat=ierr)
    if(ierr/=0) then
       print 333, trim(filename), 'proc', ierr
       return
    endif

    rewind(50)
    read(50,nml=input,iostat=ierr)
    if(ierr/=0) then
       print 333, trim(filename), 'input', ierr
       return
    endif

    rewind(50)
    read(50,nml=output,iostat=ierr)
    if(ierr/=0) then
       print 333, trim(filename), 'input', ierr
       return
    endif

    close(50)

    ! Check the weight method
20  format(A)
    if(weight_method==1) then
       print 20,'Weighting method = 1 = members are equal.  Weight = 1.0.'
    elseif(weight_method==2) then
       print 20,'Weighting method = 2 = models are equal, members are not.'
       print 20,'                                number of members present for model A'
       print 20,' weight for member of model A = -------------------------------------'
       print 20,'                                  number of members for all models'
    else
       print 25,weight_method
       write(0,25) weight_method
       stop 2
25     format('ERROR: Invalid weighting method ',I0,' must be 1 or 2.')
    endif

    ! Figure out which models exist.
    modelcount=0
    do i=1,maxmodel
30     format('Model #',I0,': ',A)
       good=.true.
       print 30,i,'checking namelist parameters...'
301    format('  name="',A,'" glob="',A,'" ensmax=',I0)
       print 301,trim(name(i)),trim(atcf_glob(i)),ensmax(i)
       if(trim(name(i))=='#') then
          print 30,i,'model name is unspecified'
          good=.false.
       endif
       if(trim(atcf_glob(i))=='#') then
          print 30,i,'atcf_glob is unspecified'
          good=.false.
       endif
       if(ensmax(i)<ensmin(i) .or. ensmin(i)<minensmin .or. &
            ensmax(i)>maxensmax) then
          print 31,i,ensmin(i),ensmax(i)
          good=.false.
31        format('Model #',I0,': invalid ensemble number bounds ensid=',&
               I0,'...',I0,'.')
       endif
       if(.not.good) then
          print 33,i
33        format('Model #',I0,': will not process this model.')
          mm%modelflag(i)=.false.
          cycle
       endif
       if(trim(subname(i))=='#') then
          subname(i)=name(i)
37        format('Model #',I0,' has no submodel name.  ',&
               'Will use model name "',A,'" as submodel name.')
          print 37,i,trim(name(i))
       endif
38     format('Model #',I0,' name="',A,'" subname="',A,'" ensid=',I0,&
            '...',I0,' (',I0,' members)')
       print 38,i,trim(name(i)),trim(subname(i)),ensmin(i),ensmax(i),&
            ensmax(i)-ensmin(i)+1
       modelcount=modelcount+1
       mm%modelflag(i)=.true.
    end do

    ! Copy inputs to mm object and subobjects.
    do i=1,maxmodel
       if(mm%modelflag(i)) then
          print 50,i
50        format('Model #',I0,': initializing from namelist values.')
          call mm%mod(i)%set_values(ensmin(i),ensmax(i),stdev0(i),&
               (stdev120(i)-stdev0(i))/120.0,atcf_glob(i),&
               name(i),subname(i),cendist(i),maxdist(i))
       else
          print 100,i
100       format('Model #',I0,': not processing this.')
       endif
    enddo
    mm%out_pgm_format=out_pgm_format
    mm%out_dat_format=out_dat_format
    mm%out_ctl_format=out_ctl_format
    mm%out_times=out_times
    mm%weight_method=weight_method
    mm%latres=max(0.001,latres)
    mm%lonres=max(0.001,lonres)
    mm%sent=.false.
    mm%time_step=time_step_hrs
    mm%last_ftime=last_fhr

    ! Write namelist back out. 
    open(unit=50,file=trim(outname),form='FORMATTED')
    write(50,nml=proc)
    write(50,nml=input)
    write(50,nml=output)
    close(50)

    ! List models to process
    do i=1,maxmodel
       if(mm%modelflag(i)) then
200       format('Process model #',I0,'.  Description:')
          print 200,i
          call mm%mod(i)%print(6)
       endif
    enddo
  end subroutine read_namelist

  subroutine print_model(m,unit)
    class(model_atcf), intent(in) :: m
    integer, intent(in) :: unit
    print 200,m%name,m%subname,trim(m%atcf_glob)
    print 201,m%enscount,m%weight
    print 202,m%stdev0,m%stdev_per_hr,m%cendist,m%maxdist
200 format('Model name=',A4,' subname=',A4,' glob=',A)
201 format('      enscount=',I0,' weight=',F0.5)
202 format('      stdev0=',F0.1,' stdev_per_hr=',F0.1,' cendist=',F0.1, &
           ' maxdist=',F0.1)
  end subroutine print_model

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine set_weight_model(m,weight)
    class(model_atcf) :: m
    real, intent(in) :: weight
    m%weight=weight
  end subroutine set_weight_model

  subroutine set_weight_multimodel(mm,weight)
    class(multimodel) :: mm
    real, intent(in) :: weight
    integer :: i
    do i=1,maxmodel
       if(.not.mm%modelflag(i)) cycle
       call mm%mod(i)%set_weight(weight)
    end do
  end subroutine set_weight_multimodel
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_atcf_member(m,filename,max_atcf_records,ierr,&
       reflon,have_reflon,minlat,maxlat,minlon,maxlon,atime)
    use tcf_module, only : atcf, read_atcf
    character*(*), intent(in) :: filename
    real, intent(inout) :: reflon,minlat,maxlat,minlon,maxlon
    real(kind=C_DOUBLE) :: atime
    logical, intent(inout) :: have_reflon
    class(member_atcf) :: m
    type(atcf) :: atcfs(max_atcf_records)
    integer, intent(in) :: max_atcf_records
    integer, intent(out) :: ierr

    integer :: used, i, iread

20  format('ERROR: ',A,': ',A)

    call read_atcf(trim(filename),atcfs,max_atcf_records,used, iostat=iread)
    if(iread/=0) then
       print 20,trim(filename),'i/o error or file not found'
       goto 100 ! error handler
    endif

    if(used<1) then
       ierr=ERR_ATCF_NO_VALID_LINES
       print 20,trim(filename),'no valid ATCF track found in file.'
       goto 100 ! error handler
    elseif(atcfs(1)%fcsthour>1e-5) then
       ierr=ERR_ATCF_FIRST_TIME_NOT_0FHR
       print 20,trim(filename),'first ATCF time is not 0 hours.  '&
            //'Assuming file is corrupted.'
       goto 100 ! error handler
    else
       do i=2,used
          if(atcfs(i)%fcsthour < atcfs(i-1)%fcsthour) then
30           format('ERROR: ',A,': multiple tracks detected in file.  '&
                  //'Assuming file is corrupted.')
             print 30,trim(filename)
             ierr=ERR_ATCF_MULTIPLE_TRACKS
             goto 100 ! error handler
          endif
       enddo
    endif

    call m%free()
    ierr=0
    allocate(m%lats(used), m%lons(used), m%times(used))
    m%ntimes=used
    m%alloced=used
    do i=1,used
       m%lats(i)=atcfs(i)%lat
       if(i==1) then
          if(.not.have_reflon) then
             reflon=atcfs(i)%lon
             atime=atcfs(i)%epoch
             have_reflon=.true.
          else
             if(abs(atime-atcfs(i)%epoch)>30) then
                ! atime varies by more than 30 seconds
                ierr=ERR_ATCF_ATIME_MISMATCH
37              format(A,': analysis time does not match other files.  Discarding this file.')
                print 37,trim(filename)
                goto 100 ! error handler
             endif
          endif
       endif
       m%lons(i)=mod(3600+180+atcfs(i)%lon-reflon,360.)-180+reflon
       m%times(i)=atcfs(i)%fcsthour

       minlon=min(minlon,m%lons(i))
       maxlon=max(maxlon,m%lons(i))
       minlat=min(minlat,m%lats(i))
       maxlat=max(maxlat,m%lats(i))
    enddo
    m%last_time=1
    return

100 continue ! Error handler
    call m%free()
  end subroutine read_atcf_member

  subroutine read_atcf_model(m,max_atcf_records,ierr,&
       reflon,have_reflon,minlat,maxlat,minlon,maxlon,atime)
    use fileop_module, only: glob, init_glob
    class(model_atcf) :: m
    real(kind=C_DOUBLE) :: atime
    real, intent(inout) :: reflon,minlat,maxlat,minlon,maxlon
    logical, intent(inout) :: have_reflon
    integer, intent(out) :: ierr
    integer, intent(in) :: max_atcf_records
    !Locals
    integer :: ensid, ierr1, so_far,gerr
    character(len=filenamelen) :: filename
    type(glob) :: search

    ensid=m%ensmin
    ierr=0
    m%enscount=0
28  format(A,': ATCF found here.  Will use ensid=',I0,'.')
29  format(A,': ATCF file appears invalid: ierr1=',I0)
30  format(A,': skipping; ran out of ensemble entries for this model.')
    call init_glob(search,trim(m%atcf_glob),gerr)
    if(gerr/=0) then
20     format(A,': nothing found: ',I0)
       print 20,trim(m%atcf_glob),gerr
       ierr=gerr
       return
    endif
    do while(search.next())
       if(ensid<=m%ensmax) then
          call m%atcf(ensid)%read_atcf(search%match,max_atcf_records,ierr1,&
               reflon,have_reflon,minlat,maxlat,minlon,maxlon,atime)
          if(ierr1==0) then
             print 28,search%match,ensid
             m%ensflag(ensid)=.true.
             m%enscount=m%enscount+1
             ensid=ensid+1
          else
             print 29,search%match,ierr1
          endif
       else
          print 30,search%match
       endif
    enddo

    so_far=ensid
    do ensid=so_far+1,m%ensmax
       m%ensflag(ensid)=.false.
    enddo
    
    if(m%enscount==0) then
38     format(A)
       ierr=1
    endif
  end subroutine read_atcf_model

  subroutine read_atcf_multimodel(mm,max_atcf_records,ierr)
    class(multimodel) :: mm
    integer, intent(in) :: max_atcf_records
    integer, intent(inout) :: ierr

    integer :: i, j, nx, ny, ensaccum
    real :: reflon,minlat,maxlat,minlon,maxlon
    logical :: have_reflon
    real :: midlon,midlat,lat1,lon1 ,latres,lonres
    real(kind=C_DOUBLE) :: atime

    latres=mm%latres
    lonres=mm%lonres

    ensaccum=0
    atime=-9e9
    reflon=-9e9
    minlat=9e9
    maxlat=-9e9
    minlon=9e9
    maxlon=-9e9
    have_reflon=.false.

    do i=1,maxmodel
       ierr=-999
       if(.not.mm%modelflag(i)) cycle
       call mm%mod(i)%read_atcf(max_atcf_records,ierr,&
            reflon,have_reflon,minlat,maxlat,minlon,maxlon,atime)
       if(ierr/=0) then
38        format('ERROR: Model #',I0,' (model "',A,'" submodel "',A,&
               '"): no valid data: ierr=',I0)
          print 38, i, trim(mm%mod(i)%name), trim(mm%mod(i)%subname), ierr
          write(0,38) i, trim(mm%mod(i)%name), trim(mm%mod(i)%subname), ierr

          mm%modelflag(i)=.false.
          call mm%mod(i)%free()
          cycle
       endif
39     format('Model #',I0,' (model "',A,'" submodel "',A,&
              '"): have ',I0,' members with valid data.')
       print 39,i,trim(mm%mod(i)%name),trim(mm%mod(i)%subname),&
            mm%mod(i)%enscount
       ensaccum=ensaccum+mm%mod(i)%enscount
       if(mm%weight_method==1) then
          call mm%set_weight(1.0)
       endif
    end do
    if(mm%weight_method/=1) then
       do i=1,maxmodel
          call mm%mod(i)%set_weight(real(mm%mod(i)%enscount)/ensaccum)
       enddo
    endif

    if(.not.have_reflon) then
40     format(A)
       write(0,40) 'ERROR: no ATCF files found.  Giving up.'
       write(6,40) 'ERROR: no ATCF files found.  Giving up.'
       ierr=40
       return
    endif

41  format('Result: reflon=',F0.2,' minlon=',F0.2,' maxlon=',F0.2)
42  format('        atime=',F0.1', minlat=',F0.2,' maxlat=',F0.2)
    print 41,reflon,minlon,maxlon
    print 42,atime,minlat,maxlat

    ! Allocate the grid
    mm%atime=atime
    midlon=(minlon+maxlon)/2
    midlat=(minlat+maxlat)/2
    ny=ceiling( (maxlat-minlat)/latres/2+1 )*2
    nx=ceiling( (maxlon-minlon)/lonres/2+1 )*2
    if(nx<2 .or. ny<2) then
45     format('ERROR: Bad grid dims: nx=',I0,' ny=',I0)
       write(0,45) nx,ny
       write(6,45) nx,ny
       ierr=30
       return
    endif
    mm%ny=ny
    mm%nx=nx
    lat1=midlat-ny/2*latres
    lon1=midlon-nx/2*lonres
54  format('Decided lat1=',F0.2,' lon1=',F0.2,' nx=',I0,' ny=',I0)
    print 54,lat1,lon1,nx,ny
    call init_proj_latlon(mm%proj,lat1,lon1,latres,lonres,nx,ny,.true.)
    call init_decomp(mm%dec,  1,nx,1,ny,1,1,  1,nx,1,ny,1,1)
    allocate(mm%prob,mm%work)
    call full_init_vardata_real(mm%prob,mm%proj,mm%dec,.true.,.true.,&
         .false.,.false.)
    call full_init_vardata_real(mm%work,mm%proj,mm%dec,.true.,.true.,&
         .false.,.false.)

    !$OMP PARALLEL DO PRIVATE(i,j)
    do j=1,ny
       do i=1,nx
          mm%prob%rdata(i,j,1)=0
          mm%prob%mask(i,j,1)=.false.
       enddo
    enddo

    !$OMP PARALLEL DO PRIVATE(i,j)
    do j=1,ny
       do i=1,nx
          mm%work%rdata(i,j,1)=0
          mm%work%mask(i,j,1)=.false.
       enddo
    enddo

    allocate(mm%lats(ny), mm%lons(nx))
    do j=1,ny
       mm%lats(j)=lat1+(j-1)*latres
    enddo
    do i=1,nx
       mm%lons(i)=lon1+(i-1)*lonres
    enddo

407 format('Grid lat=[',F0.1,'..',F0.1,'] lon=[',F0.1,'..',F0.1,']')
    print 407,mm%lats(1),mm%lats(ny),mm%lons(1),mm%lons(nx)
408 format('     reflon=',F0.1,' latres=',F0.5,' lonres=',F0.5)
    print 408,reflon,latres,lonres

    ierr=0
  end subroutine read_atcf_multimodel

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine update_prob_member(m,prob,when,lats,lons,maxdist,cendist, &
       stdev0, stdev_per_hr, minprob)
    use ens_prob_module, only: dist_weight
    class(member_atcf) :: m
    class(vardata_real) :: prob
    real, intent(in) :: when
    real, intent(in) :: maxdist,cendist,stdev0,stdev_per_hr,minprob
    real, pointer, intent(in) :: lons(:)
    real, pointer, intent(in) :: lats(:)
    ! Locals
    real :: tnext,stdev1,stdev2
    integer :: i,j,inext,iend
    
    do ! infinite loop, terminated by the IF block inside it
       inext=m%last_time+1
       iend=inext+1
       m%last_time=inext
       if(m%ntimes==1) then
          ! Special case: track has one entry.  Use distance to point
          ! instead of distance to line segment.
          iend=inext
       endif
       if(iend>m%ntimes) then
10        format('Time index #',I0,' hit maximum ',I0,'.  Done with this track.')
          print 10,inext-1,m%ntimes
          m%last_time=m%ntimes
          return
       elseif(when+1e-5<m%times(inext)) then
          print 20,when,m%times(inext)
          m%last_time=inext-1
20        format('Time ',F0.1,'hrs too early for timespan end ',F0.1,'.')
          return
       endif

       stdev1=stdev0+m%times(inext)*stdev_per_hr
       stdev2=stdev0+m%times(iend)*stdev_per_hr

30     format('Call dist_weight: nx=',I0,' ny=',I0,' lat1=',F0.1,' lon1=',F0.1)
31     format('                  lat2=',F0.1,' lon2=',F0.1,' maxdist=',F0.1)
32     format('                  cendist=',F0.1,' stdev1=',F0.1', stdev2=',F0.1)
33     format('                  minprob=',F0.1,' for times ',I0,'..',I0)

       print 30,prob%dc%ide-prob%dc%ids+1, prob%dc%jde-prob%dc%jds+1,&
            m%lats(inext),m%lons(inext)
       print 31,m%lats(iend),m%lons(iend),maxdist
       print 32,cendist,stdev1,stdev2
       print 33,minprob,inext,iend

       call dist_weight(prob%dc%ide-prob%dc%ids+1, prob%dc%jde-prob%dc%jds+1, &
            lats,lons, prob%rdata(:,:,1), prob%mask(:,:,1), &
            m%lats(inext),m%lons(inext), m%lats(iend),m%lons(iend), &
            maxdist, cendist, stdev1, cendist, stdev2, minprob)

    enddo

  end subroutine update_prob_member

  subroutine update_prob_model(m,prob,when,lats,lons,minprob)
    class(model_atcf) :: m
    class(vardata_real) :: prob
    real, intent(in) :: when
    real, intent(in) :: minprob
    real, pointer, intent(in) :: lons(:)
    real, pointer, intent(in) :: lats(:)
    ! Locals
    integer :: ensid
47  format('Update probability for model name=',A,' subname=',A,' member ',&
           I0,'/',I0)
    do ensid=m%ensmin,m%ensmax
       if(.not.m%ensflag(ensid)) cycle
       print 47,trim(m%name),trim(m%subname),ensid,m%enscount
       call m%atcf(ensid)%update_prob(prob,when,lats,lons,m%maxdist,m%cendist,&
                                      m%stdev0,m%stdev_per_hr,minprob)
    enddo
  end subroutine update_prob_model

  subroutine update_prob_multimodel(mm,when)
    class(multimodel) :: mm
    integer :: i,j,m
    real, intent(in) :: when

    ! !$OMP PARALLEL DO PRIVATE(i,j)
    ! do j=1,mm%ny
    !    do i=1,mm%nx
    !       mm%prob%mask(i,j,1)=.false.
    !       mm%prob%rdata(i,j,1)=0.
    !       mm%work%mask(i,j,1)=.false.
    !       mm%work%rdata(i,j,1)=0.
    !    enddo
    ! enddo
    ! !$OMP END PARALLEL DO

    model: do m=1,maxmodel
       if(.not.mm%modelflag(m)) cycle
       call mm%mod(m)%update_prob(mm%prob,when,mm%lats,mm%lons,mm%minprob)
    enddo model
  end subroutine update_prob_multimodel

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  logical function time_for_output(mm,ftime)
    class(multimodel) :: mm
    real, intent(in) :: ftime
    integer :: otime

    do otime=1,mm%notimes
       if(.not.mm%sent(otime) .and. ftime+.05>mm%out_times(otime)) then
          time_for_output=.true.
          return
       endif
    enddo

    time_for_output=.false.
  end function time_for_output

  logical function next_time(mm,ftime)
    class(multimodel) :: mm
    integer :: itime
    real, intent(inout) :: ftime
    character(len=filenamelen) :: pgm

    next_time=.true.

    ftime=mm%ftime
    if(mm%time_for_output(ftime)) return

    do while(mm%ftime+mm%time_step<mm%last_ftime+.05)
300    format('next_time: time step +',F0.3,' ftime now ',F0.3,'.')
       mm%ftime=mm%ftime+mm%time_step
       ftime=mm%ftime
       print 300,mm%time_step,ftime
       call mm%update_prob(mm%ftime)
       if(len_trim(mm%out_pgm_format)>1) then
          call insert_ens_domain(mm%out_pgm_format,pgm,ifhr=nint(ftime))
          call write_pgm(pgm,mm%nx,mm%ny,mm%prob%mask,mm%prob%rdata,mult=200.)
       endif
       if(mm%time_for_output(mm%ftime)) return
    enddo

400 format('next_time: terminating at ftime=',F0.1)
    print 400,ftime
    next_time=.false. ! No more times to process
  end function next_time

end module atcf_prob_module
