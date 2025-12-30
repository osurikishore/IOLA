program hwrf_ens_prob
  use ens_output_module, only: grads_pgm, init_grads_pgm
  use ens_prob_module, only: ensemble, init_ensemble, mindommin,maxdommax, &
       max_output_times, max_thresholds, filenamelen

  implicit none
  
  character(len=filenamelen) :: &
       pgmone_format, pgmany_format, windpgm_format, windprob_format,&
       rainpgm_format, rainprob_format, rainswath_format, &
       grads_data_file, grads_ctl_file, grads_tornado_data, &
       grads_tornado_ctl, tornado_format

  ! FIXME: ADD THESE TWO:
  character(len=filenamelen) :: ensprob_file, accumpgm_format

  real :: undef=-999

  character(len=filenamelen) :: indiag_format, inatcf_format
  integer(kind=8) :: min_wrfdiag_sizes(maxdommax-mindommin+1)
  real :: lonres, latres, stdev0, stdev120
  real(kind=8) :: min_input_age
  real :: sd_per_hour, cendist, minprob, maxdist, tg_freq
  real :: output_times(max_output_times)
  real :: wind_thresholds(max_thresholds)
  real :: rain_thresholds(max_thresholds), okdist0, okdist_per_hour
  integer :: trackdom,dommin,dommax,ensmin,ensmax
  integer :: num_wind_thresh,num_rain_thresh,num_output_times

  type(ensemble) :: ens
  type(grads_pgm) :: out

  namelist /proc/ lonres,latres, dommin,dommax, trackdom, &
       ensmin,ensmax, stdev0,stdev120, cendist,maxdist,minprob, &
       okdist0,okdist_per_hour, wind_thresholds,rain_thresholds
  namelist /input/ indiag_format, inatcf_format, min_wrfdiag_sizes, &
       min_input_age
  namelist /output/ pgmone_format,pgmany_format,windpgm_format,&
       windprob_format,rainpgm_format,rainprob_format, &
       rainswath_format,ensprob_file,grads_data_file,grads_ctl_file, &
       output_times,grads_tornado_data,grads_tornado_ctl, tornado_format, &
       accumpgm_format,tg_freq

  ! FIXME: ADD THIS: accumpgm_format

  call read_namelist()

  call init_grads_pgm(out,&
       output_times,num_output_times,wind_thresholds,num_wind_thresh,&
       rain_thresholds,num_rain_thresh,tg_freq,&
       pgmone_format, pgmany_format, windpgm_format, windprob_format,&
       rainpgm_format, rainprob_format, rainswath_format, &
       grads_data_file, grads_ctl_file, grads_tornado_data, &
       grads_tornado_ctl,tornado_format)

  call init_ensemble(ens,&
       indiag_format, inatcf_format, ensmin,ensmax, dommin,dommax, &
       trackdom, latres,lonres)

  if(.not.ens%check_input(min_wrfdiag_sizes,min_input_age)) then
10   format(A)
     print 10,'HAVE NO INPUT!! Aborting.'
     write(0,10) 'HAVE NO INPUT!! Aborting.'
     stop 22
  endif

  call ens%run(out, &
       stdev0, sd_per_hour, maxdist, cendist, minprob,&
       ensprob_file, .true., okdist0, okdist_per_hour)

contains
  subroutine read_namelist()
    integer :: i,t, o
    lonres=0.03
    latres=0.03
    dommin=2
    dommax=3
    tg_freq=3.0
    indiag_format='*'
    inatcf_format='*'
    !pgmone_format  ='probens_m<ensid>_f<ifhr>.pgm'
    pgmone_format  ='*'
    !pgmany_format  ='probability_f<ifhr>.pgm'
    pgmany_format  ='*'
    !accumpgm_format='rainaccum_f<ifhr>.pgm' ! rain accumulator array for mean
    accumpgm_format='*'
    !windpgm_format ='windens_m<ensid>_f<ifhr>.pgm'
    windpgm_format ='*'
    !windprob_format='windprob_<thresh>mps_f<ifhr>.pgm'
    windprob_format='*'
    !rainpgm_format ='rainens_m<ensid>_f<ifhr>.pgm'
    rainpgm_format ='*'
    !rainprob_format='rainprob_<thresh>cm_f<ifhr>.pgm'
    rainprob_format='*'
    rainswath_format='*'
    ensprob_file='ensprob.lst'

    !tornado_format='tg_<var>_f<ifhr>.pgm'
    tornado_format='*'

    grads_tornado_data='hwrf_hourly.dat'
    grads_tornado_ctl='hwrf_hourly.ctl'
    grads_data_file='hwrf_ens_prob_f<ifhr>.dat'
    grads_ctl_file='hwrf_ens_prob_f<ifhr>.ctl'
    
    stdev0=10e3
    stdev120=130e3;
    cendist=20e3
    maxdist=1110e3  ! ~10 degrees
    minprob=0.05
    ensmin=1
    ensmax=3
    
    num_wind_thresh=0
    wind_thresholds=-999
    !wind_thresholds(1:3)=(/ 34, 50, 65 /) * 1852 / 3600
    
    num_rain_thresh=0
    rain_thresholds=-999
    !rain_thresholds(1:6)=(/ 2, 4, 6, 10, 15, 20 /) * 25.4e-3
    
    okdist0=30e3
    okdist_per_hour=300 ! meters
    
    num_output_times=0
    output_times=-999
    trackdom=3

    open(unit=50,file='ens_prob.nml',form='formatted')
    read(50,nml=proc)
    rewind(50)
    read(50,nml=input)
    rewind(50)
    read(50,nml=output)
    close(50)

    open(unit=51,file='ens_prob.nml.out',form='formatted')
    write(51,nml=proc)
    write(51,nml=input)
    write(51,nml=output)
    close(51)

    
    if(len_trim(indiag_format)<2) then
       write(0,*) 'ERROR: Provide &input indiag_format in namelist.'
       stop 2
    elseif(len_trim(inatcf_format)<2) then
       write(0,*) 'ERROR: Provide &input inatcf_format in namelist.'
       stop 2
    endif

    do o=1,max_output_times
       if(output_times(o)>=0) then
12        format('Output time ',I0,' is ',F0.3)
          print 12,o,output_times(o)
       else
          num_output_times=o-1
          print 13,num_output_times
13        format('There are ',I0,' output times.')
          exit
       endif
    enddo

    if(num_output_times<1) then
       write(0,*) 'ERROR: Provide at least one output time.'
       stop 2
    endif

    do t=1,max_thresholds
       if(wind_thresholds(t)>=0) then
14        format('Wind threshold ',I0,' is ',F0.3)
          print 14,o,wind_thresholds(t)
       else
          num_wind_thresh=t-1
          print 15,num_wind_thresh
15        format('There are ',I0,' wind thresholds.')
          exit
       endif
    enddo

    if(num_wind_thresh<1) then
       write(0,*) 'ERROR: Provide at least one wind threshold.'
       stop 2
    endif

    do t=1,max_thresholds
       if(rain_thresholds(t)>=0) then
16        format('Rain threshold ',I0,' is ',F0.3)
          print 16,o,rain_thresholds(t)
       else
          num_rain_thresh=t-1
          print 17,num_rain_thresh
17        format('There are ',I0,' rain thresholds.')
          exit
       endif
    enddo

    if(num_rain_thresh<1) then
       write(0,*) 'ERROR: Provide at least one rain threshold.'
       stop 2
    endif

    sd_per_hour=(stdev120-stdev0)/120.0
  end subroutine read_namelist
  
end program hwrf_ens_prob
