module ens_output_module
  use ens_prob_module, only: outputter, init_outputter, insert_ens_domain, &
       filenamelen
  use vardata_module, only: vardata_real
  implicit none
  private
  public :: grads_pgm, init_grads_pgm, write_pgm

  integer, parameter :: ntgvars=5 ! number of tornado genesis vars

  type, extends(outputter) :: grads_pgm
     character(len=filenamelen) :: &
          pgmone_format, pgmany_format, windpgm_format, windprob_format,&
          rainpgm_format, rainprob_format, rainswath_format, &
          grads_swath_data, grads_swath_ctl, grads_tornado_data, &
          grads_tornado_ctl, tornado_format
     real :: undef=-999.0
     integer :: grads_first_unit=3000
   contains
     procedure free
     
     procedure start
     procedure step_time
     procedure finish
     
     procedure trackprob
     procedure rainswath
     procedure windprob
     procedure rainprob
     !procedure udhelfrac
     procedure udhel
     procedure maxw
     procedure minw
     procedure maxu10
     procedure member_rain
     procedure member_wind
     procedure member_prob

     procedure, private :: write_grads
  end type grads_pgm
  
contains

  subroutine free(this)
    class(grads_pgm) :: this
    ! Nothing to free yet...
  end subroutine free

  subroutine finish(this)
    class(grads_pgm) :: this
    integer :: ixunit,iunit
    do ixunit=0,this%notimes
       iunit=this%grads_first_unit+ixunit
       print *,'close unit ',iunit
       close(iunit)
    end do
  end subroutine finish

  subroutine start(this,lats,lons,atime,nx,ny,fcstlen)
    use datetime_module, only: grads_date_len, grads_date

    integer, intent(in) :: nx, ny
    real, intent(in) :: lats(ny), lons(nx)
    real(kind=8), intent(in) :: atime,fcstlen

    class(grads_pgm) :: this

    ! Locals

    character(len=grads_date_len) :: stratime
    character(len=filenamelen) :: data_filename, ctl_filename
    real, allocatable :: junk(:,:)
    integer :: i,j,itime,ctlu,ifhr,nvars, grads_unit, ithresh

    ! Magic starts here

    nvars=2+this%nrthresh+this%nwthresh
    call grads_date(stratime,atime)

    allocate(junk(nx,ny))
    !$OMP PARALLEL DO PRIVATE(i,j)
    do j=1,ny
       do i=1,nx
          junk(i,j)=this%undef
       enddo
    enddo
    !$OMP END PARALLEL DO

    do itime=0,this%notimes
       grads_unit=this%grads_first_unit+itime
       if(itime==0) then ! tornado genesis file
          ifhr=nint(this%output_times(1))
          call insert_ens_domain(this%grads_tornado_data,data_filename,0,0,ifhr=ifhr)
       else
          ifhr=nint(this%output_times(itime))
          call insert_ens_domain(this%grads_swath_data,data_filename,0,0,ifhr=ifhr)
       endif
48     format('Unit ',I0,' recl=',I0,' open data file="',A,'"')
       print 48,grads_unit,nx*ny,trim(data_filename)
       open(grads_unit,access='DIRECT',recl=nx*ny,file=trim(data_filename))
       do i=1,nvars
          write(grads_unit,rec=i) junk
       enddo
       if(itime==0) then ! tonardo genesis control file
          call insert_ens_domain(this%grads_tornado_ctl,ctl_filename,0,0,ifhr=ifhr)
       else
          call insert_ens_domain(this%grads_swath_ctl,ctl_filename,0,0,ifhr=ifhr)
       endif
       ctlu=1000+grads_unit
33     format('Unit ',I0,': open CTL file "',A,'"')
       print 33,ctlu,trim(ctl_filename)
       open(ctlu,form='FORMATTED',file=trim(ctl_filename))
       write(ctlu,10) trim(data_filename)
       write(ctlu,20) ifhr,this%undef
       write(ctlu,30) nx,lons(1),(lons(nx)-lons(1))/(nx-1)
       write(ctlu,40) ny,lats(1),(lats(ny)-lats(1))/(ny-1)
       write(ctlu,50) ' '
       whichfile: if(itime==0) then ! tornado genesis control file
          write(ctlu,60) floor((fcstlen+1e-5)/this%tg_freq),trim(stratime)
          write(ctlu,70) ntgvars
          write(ctlu,75) 'udhel03','Max. 0-3km Updraft Helicity (m2/s2)'
          write(ctlu,75) 'udhel25','Max. 2-5km Updraft Helicity (m2/s2)'
          write(ctlu,75) 'minw','Minimum Vertical Motion Below 400mbar (m/s)'
          write(ctlu,75) 'maxw','Maximum Vertical Motion Below 400mbar (m/s)'
          write(ctlu,75) 'maxu10','Maximum 10m Altitude Wind (m/s)'
       else
          write(ctlu,60) 1,trim(stratime)
          write(ctlu,70) nvars
          write(ctlu,80) ' '
          write(ctlu,81) ' '
          do i=1,this%nrthresh
             ithresh=nint(this%rthresh(i)*100)
             if(ithresh>=100) then
                write(ctlu,93) 'rain',ithresh,'cm','rain accumulation',ithresh,'cm'
             elseif(ithresh>=10) then
                write(ctlu,92) 'rain',ithresh,'cm','rain accumulation',ithresh,'cm'
             else
                write(ctlu,91) 'rain',ithresh,'cm','rain accumulation',ithresh,'cm'
             endif
          enddo
          do i=1,this%nwthresh
             ithresh=nint(this%wthresh(i))
             if(ithresh>=10) then
                write(ctlu,95) 'wind',ithresh,'mps','wind',ithresh,'mps'
             else
                write(ctlu,94) 'wind',ithresh,'mps','wind',ithresh,'mps'
             endif
          enddo
       endif whichfile
       write(ctlu,100) ' '
       close(ctlu)
    enddo

10  format("DSET  ^",A)
20  format("TITLE Swath Data Fcst Hour ",I0,/,"UNDEF ",F0.1)
30  format("XDEF ",I0," LINEAR ",F0.5,"  ",F0.5)
40  format("YDEF ",I0," LINEAR ",F0.5,"  ",F0.5)
50  format("ZDEF 1 LEVELS 1000",A)
60  format("TDEF ",I0," LINEAR ",A," 1hr")
70  format("VARS ",I0)
75  format(A10,      "       0 99 ",A)
80  format("track_prob       0 99 track probability",A)
81  format("rain_swath       0 99 mean rain swath",A)
    !       rain_prob_000cm  0 99 probability of more rain than 000cm"
91  format(A,"_prob_00",I1,A,"  0 99 probability of ",A," greater than   ",I1,A)
92  format(A,"_prob_0", I2,A,"  0 99 probability of ",A," greater than  ", I2,A)
93  format(A,"_prob_",  I3,A,"  0 99 probability of ",A," greater than ",  I3,A)
94  format(A,"_prob_0", I1,A,"  0 99 probability of ",A," greater than  ", I1,A)
95  format(A,"_prob_",  I2,A,"  0 99 probability of ",A," greater than ",  I2,A)
100 format('ENDVARS',A)

    deallocate(junk)

  end subroutine start

  ! ----------------------------------------------------------------------------

  subroutine step_time(this,itime,rtime)
    class(grads_pgm) :: this
    integer, intent(in) :: itime
    real, intent(in) :: rtime
  end subroutine step_time

  ! ----------------------------------------------------------------------------

  subroutine init_grads_pgm(this,&
       otimes,notimes,wthresh,nwthresh,rthresh,nrthresh,tg_freq,&
       pgmone_format, pgmany_format, windpgm_format, windprob_format,&
       rainpgm_format, rainprob_format, rainswath_format, &
       grads_swath_data, grads_swath_ctl, grads_tornado_data, &
       grads_tornado_ctl,tornado_format)
    character(len=filenamelen), intent(in) :: &
         pgmone_format, pgmany_format, windpgm_format, windprob_format,&
         rainpgm_format, rainprob_format, rainswath_format, &
         grads_swath_data, grads_swath_ctl, grads_tornado_data, &
         grads_tornado_ctl, tornado_format
    class(grads_pgm) :: this
    real, intent(in) :: otimes(notimes),rthresh(nrthresh),wthresh(nwthresh)
    integer, intent(in) :: notimes,nrthresh,nwthresh
    real, intent(in) :: tg_freq
    call init_outputter(this,otimes,notimes,wthresh,nwthresh,&
         rthresh,nrthresh,tg_freq)
    this%pgmone_format=pgmone_format
    this%pgmany_format=pgmany_format
    this%windpgm_format=windpgm_format
    this%windprob_format=windprob_format
    this%rainpgm_format=rainpgm_format
    this%rainprob_format=rainprob_format
    this%rainswath_format=rainswath_format
    this%grads_swath_data=grads_swath_data
    this%grads_swath_ctl=grads_swath_ctl
    this%grads_tornado_data=grads_tornado_data
    this%grads_tornado_ctl=grads_tornado_ctl
    this%tornado_format=tornado_format
  end subroutine init_grads_pgm

  ! ----------------------------------------------------------------------------

  subroutine trackprob(this,data,itime)
    class(grads_pgm) :: this
    class(vardata_real) :: data
    integer, intent(in) :: itime
    character(len=filenamelen) :: pgm
    integer :: ifhr, nx,ny
    ifhr=nint(this%output_times(itime))
    nx=data%dc%ide-data%dc%ids+1
    ny=data%dc%jde-data%dc%jds+1
    if(len_trim(this%pgmany_format)>1) then
       call insert_ens_domain(this%pgmany_format,pgm,99,0,ifhr)
       call write_pgm(trim(pgm),nx,ny,data%mask(:,:,1),data%rdata(:,:,1))
    endif
    call this%write_grads(data,itime,ivar=1)
  end subroutine trackprob

  ! ----------------------------------------------------------------------------

  subroutine rainswath(this,data,itime)
    class(grads_pgm) :: this
    class(vardata_real) :: data
    integer, intent(in) :: itime
    character(len=filenamelen) :: pgm
    integer :: ifhr, nx,ny
    ifhr=nint(this%output_times(itime))
    nx=data%dc%ide-data%dc%ids+1
    ny=data%dc%jde-data%dc%jds+1
    if(len_trim(this%rainswath_format)>1) then
59     format('Write swath for hour ',I0,' to "',A,'".')
       call insert_ens_domain(this%rainswath_format,pgm,99,0,ifhr)
       print 59,ifhr,trim(pgm)
       call write_pgm(trim(pgm),nx,ny,data%mask(:,:,1), &
            data%rdata(:,:,1),scale=.false.,mult=100.)
    endif
    call this%write_grads(data,itime,ivar=2)
  end subroutine rainswath

  ! ----------------------------------------------------------------------------

  subroutine windprob(this,data,itime,ithresh,rthresh)
    class(grads_pgm) :: this
    class(vardata_real) :: data
    integer, intent(in) :: itime
    integer, intent(in) :: ithresh  ! threshold index
    real, intent(in) :: rthresh ! threshold
    integer :: irthresh ! integer representation of threshold for filename
    integer :: ifhr, nx,ny
    character(len=filenamelen) :: pgm
    ifhr=nint(this%output_times(itime))
    nx=data%dc%ide-data%dc%ids+1
    ny=data%dc%jde-data%dc%jds+1
    irthresh=nint(rthresh)

    ! STEP 5: WRITE DATA
60  format('Write data for wind prob threshold #',I0,'to "',A,'".')
    if(len_trim(this%windprob_format)>1) then
       call insert_ens_domain(this%windprob_format,pgm,99,0,&
            nint(this%output_times(itime)),ithresh=irthresh)
       print 60,ithresh,trim(pgm)
       call write_pgm(trim(pgm),nx,ny,data%mask(:,:,1), &
            data%rdata(:,:,1),scale=.false.,mult=100.)
    endif
    call this%write_grads(data,itime,iwthresh=ithresh)
  end subroutine windprob

  ! ----------------------------------------------------------------------------

  subroutine rainprob(this,data,itime,ithresh,rthresh)
    class(grads_pgm) :: this
    class(vardata_real) :: data
    integer, intent(in) :: itime
    integer, intent(in) :: ithresh  ! threshold index
    real, intent(in) :: rthresh ! threshold
    character(len=filenamelen) :: pgm
    integer :: irthresh ! integer representation of threshold for filename
    integer :: ifhr, nx,ny
    ifhr=nint(this%output_times(itime))
    nx=data%dc%ide-data%dc%ids+1
    ny=data%dc%jde-data%dc%jds+1
    irthresh=nint(rthresh*100)

    ! STEP 5: WRITE DATA
60  format('Write data for rain prob threshold #',I0,'to "',A,'".')
    if(len_trim(this%rainprob_format)>1) then
       call insert_ens_domain(this%rainprob_format,pgm,99,0,&
            nint(this%output_times(itime)),ithresh=irthresh)
       print 60,ithresh,trim(pgm)
       call write_pgm(trim(pgm),nx,ny,data%mask(:,:,1), &
            data%rdata(:,:,1),scale=.false.,mult=100.)
    endif
    call this%write_grads(data,itime,irthresh=ithresh)
  end subroutine rainprob

  ! ----------------------------------------------------------------------------

  subroutine udhel(this,data,itgtime,lev1,lev2)
    class(grads_pgm) :: this
    class(vardata_real) :: data
    integer, intent(in) :: itgtime
    real, intent(in) :: lev1,lev2
    character(len=filenamelen) :: pgm
    integer :: ifhr, nx,ny
    ifhr=nint(itgtime*this%tg_freq)
    nx=data%dc%ide-data%dc%ids+1
    ny=data%dc%jde-data%dc%jds+1

    if(lev1>1990.) then
       if(len_trim(this%tornado_format)>1) then
          call insert_ens_domain(this%tornado_format,pgm,99,0,&
               ifhr,var='udhel25')
          call write_pgm(trim(pgm),nx,ny,data%mask(:,:,1), &
               data%rdata(:,:,1),scale=.false.,mult=100.)
       endif
       call this%write_grads(data,itgtime,tgvar=2)
    else
       if(len_trim(this%tornado_format)>1) then
          call insert_ens_domain(this%tornado_format,pgm,99,0,&
               ifhr,var='udhel03')
          call write_pgm(trim(pgm),nx,ny,data%mask(:,:,1), &
               data%rdata(:,:,1),scale=.false.,mult=100.)
       endif
       call this%write_grads(data,itgtime,tgvar=1)
    endif
  end subroutine udhel

  ! ----------------------------------------------------------------------------

  subroutine minw(this,data,itgtime)
    class(grads_pgm) :: this
    class(vardata_real) :: data
    integer, intent(in) :: itgtime
    character(len=filenamelen) :: pgm
    integer :: ifhr, nx,ny
    ifhr=nint(itgtime*this%tg_freq)
    nx=data%dc%ide-data%dc%ids+1
    ny=data%dc%jde-data%dc%jds+1
    if(len_trim(this%tornado_format)>1) then
       call insert_ens_domain(this%tornado_format,pgm,99,0,ifhr,var='minw')
       call write_pgm(trim(pgm),nx,ny,data%mask(:,:,1),data%rdata(:,:,1))
    endif
    call this%write_grads(data,itgtime,tgvar=3)
  end subroutine minw

  ! ----------------------------------------------------------------------------

  subroutine maxw(this,data,itgtime)
    class(grads_pgm) :: this
    class(vardata_real) :: data
    integer, intent(in) :: itgtime
    character(len=filenamelen) :: pgm
    integer :: ifhr, nx,ny
    ifhr=nint(itgtime*this%tg_freq)
    nx=data%dc%ide-data%dc%ids+1
    ny=data%dc%jde-data%dc%jds+1
    if(len_trim(this%tornado_format)>1) then
       call insert_ens_domain(this%tornado_format,pgm,99,0,ifhr,var='maxw')
       call write_pgm(trim(pgm),nx,ny,data%mask(:,:,1),data%rdata(:,:,1))
    endif
    call this%write_grads(data,itgtime,tgvar=4)
  end subroutine maxw

  ! ----------------------------------------------------------------------------

  subroutine maxu10(this,data,itgtime)
    class(grads_pgm) :: this
    class(vardata_real) :: data
    integer, intent(in) :: itgtime
    character(len=filenamelen) :: pgm
    integer :: ifhr, nx,ny
    ifhr=nint(itgtime*this%tg_freq)
    nx=data%dc%ide-data%dc%ids+1
    ny=data%dc%jde-data%dc%jds+1
    if(len_trim(this%tornado_format)>1) then
       call insert_ens_domain(this%tornado_format,pgm,99,0,ifhr,var='maxu10')
       call write_pgm(trim(pgm),nx,ny,data%mask(:,:,1),data%rdata(:,:,1))
    endif
    call this%write_grads(data,itgtime,tgvar=5)
  end subroutine maxu10

  ! ----------------------------------------------------------------------------

  subroutine member_prob(this,data,itime,imember)
    class(grads_pgm) :: this
    class(vardata_real) :: data
    integer, intent(in) :: itime
    integer, intent(in) :: imember
    character(len=filenamelen) :: pgm
    integer :: ifhr, nx,ny
    ifhr=nint(this%output_times(itime))
    nx=data%dc%ide-data%dc%ids+1
    ny=data%dc%jde-data%dc%jds+1

    if(len_trim(this%windpgm_format)>1) then
       call insert_ens_domain(this%pgmone_format,pgm,imember,0,ifhr)
16     format('Write track prob of member ',I0,' to "',A,'".')
       print 16,imember,trim(pgm)
       call write_pgm(trim(pgm),nx,ny,data%mask(:,:,1), &
            data%rdata(:,:,1),scale=.false.,takesqrt=.true.)
    endif
  end subroutine member_prob

  ! ----------------------------------------------------------------------------

  subroutine member_wind(this,data,itime,imember)
    class(grads_pgm) :: this
    class(vardata_real) :: data
    integer, intent(in) :: itime
    integer, intent(in) :: imember
    character(len=filenamelen) :: pgm
    integer :: ifhr, nx,ny
    ifhr=nint(this%output_times(itime))
    nx=data%dc%ide-data%dc%ids+1
    ny=data%dc%jde-data%dc%jds+1

    if(len_trim(this%windpgm_format)>1) then
       call insert_ens_domain(this%windpgm_format,pgm,imember,0,ifhr)
16     format('Write wind to "',A,'".')
       print 16,trim(pgm)
       call write_pgm(trim(pgm),nx,ny,data%mask(:,:,1), &
            data%rdata(:,:,1),scale=.false.,takesqrt=.true.)
    endif
  end subroutine member_wind

  ! ----------------------------------------------------------------------------

  subroutine member_rain(this,data,itime,imember)
    class(grads_pgm) :: this
    class(vardata_real) :: data
    integer, intent(in) :: itime
    integer, intent(in) :: imember
    character(len=filenamelen) :: pgm
    integer :: ifhr, nx,ny
    ifhr=nint(this%output_times(itime))
    nx=data%dc%ide-data%dc%ids+1
    ny=data%dc%jde-data%dc%jds+1

    if(len_trim(this%rainpgm_format)>1) then
       call insert_ens_domain(this%rainpgm_format,pgm,imember,0,ifhr)
16     format('Write rain to "',A,'".')
       print 16,trim(pgm)
       call write_pgm(trim(pgm),nx,ny,data%mask(:,:,1), &
            data%rdata(:,:,1),scale=.false.,takesqrt=.true.)
    endif
  end subroutine member_rain

  ! ----------------------------------------------------------------------------

  subroutine write_pgm(filename,nx,ny,mask,data,scale,mult,takesqrt)
    integer, intent(in) :: nx,ny
    logical, intent(in), optional :: scale
    real, intent(in) :: data(nx,ny)
    logical, intent(in) :: mask(nx,ny)
    character*(*), intent(in) :: filename
    real, intent(in), optional :: mult
    logical, intent(in), optional :: takesqrt

    logical :: fscale
    integer(kind=1) :: bytes(nx,ny)
    integer :: i,j, jj
    real :: rmax,rmin,D

    rmax=-9e20
    rmin=9e20

    fscale=.true.
    if(present(scale)) fscale=scale

    if(fscale) then
       !$OMP PARALLEL DO PRIVATE(i,j,D)
       minmax: do j=1,ny
          do i=1,nx
             if(mask(i,j)) then
                D=data(i,j)

                if(present(takesqrt)) then
                   if(takesqrt) then
                      D=sqrt(D)
                   endif
                endif

                if(present(mult)) then
                   D=D*mult
                endif

                if(D<rmin) rmin=D
                if(D>rmax) rmax=D
             endif
          enddo
       enddo minmax
       !$OMP END PARALLEL DO
    endif

    !$OMP PARALLEL DO PRIVATE(i,j,jj,D)
    byteloop: do j=1,ny
       do i=1,nx
          jj=ny+1-j
          D=data(i,j)

          if(present(takesqrt)) then
             if(takesqrt) then
                D=sqrt(D)
             endif
          endif

          if(present(mult)) then
             D=D*mult
          endif

          if(mask(i,j)) then
             if(fscale) then
                bytes(i,jj)=min(200,max(0,floor((D-rmin)/(rmax-rmin)*200.0)))+55
             else
                bytes(i,jj)=min(200,max(0,nint(D)))+55
             endif
          else
             bytes(i,jj) = 0
          endif
       enddo
    enddo byteloop
    !$OMP END PARALLEL DO

    open(unit=300,file=filename,form='FORMATTED')
77  format('P5',/,I0,' ',I0,' 255')
    write(300,77) nx,ny
    close(300)
    open(unit=300,file=filename,access='STREAM',position='APPEND', &
         status='OLD')
    write(300) bytes
    close(300)
  end subroutine write_pgm

  ! ----------------------------------------------------------------------------

  subroutine write_grads(this,data,itime,iwthresh,irthresh,ivar,tgvar)
    ! NOTE: This routine uses grads_unit, num_output_times, nx, ny and
    ! nvars from the program level.
    class(grads_pgm) :: this
    class(vardata_real), intent(inout) :: data
    integer, intent(in) :: itime
    integer, intent(in), optional :: ivar,iwthresh,irthresh,tgvar

    ! locals
    integer :: iloc,grads_unit
    
    !$OMP MASTER
    whichfile: if (present(tgvar)) then
       ! Write to tornado genesis file, which is always at
       ! grads_first_unit
       iloc = (itime-1)*5 + tgvar
       write(this%grads_first_unit,rec=iloc) data%rdata(:,:,1)
    else
       ! Write to non-tornado file
       iloc=1
       if(present(ivar)) then
20        format('Var #',I0,' is at rec=',I0)
          iloc=ivar
          print 20,ivar,iloc
       elseif(present(iwthresh)) then
          iloc=2+this%nrthresh+iwthresh
          print 30,iwthresh,iloc
30        format('Wind thresh #',I0,' is at rec=',I0)
       elseif(present(irthresh)) then
          iloc=2+irthresh
          print 40,irthresh,iloc
40        format('Rain thresh #',I0,' is at rec=',I0)
       endif
       grads_unit = this%grads_first_unit + itime
       write(grads_unit,rec=iloc) data%rdata(:,:,1)
    endif whichfile
    !$OMP END MASTER
  end subroutine write_grads

end module ens_output_module
