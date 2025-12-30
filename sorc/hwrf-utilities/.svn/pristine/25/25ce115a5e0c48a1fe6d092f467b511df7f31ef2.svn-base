module namelist

  !=======================================================================

  ! Define associated modules and subroutines

  !-----------------------------------------------------------------------

  use kinds

  !-----------------------------------------------------------------------

  use mpi_interface

  !-----------------------------------------------------------------------

  implicit none

  !-----------------------------------------------------------------------
  
  ! Define global variables

  character(len=500)                             :: ensemble_mean_filename             = 'NOT USED'
  character(len=500)                             :: ensemble_variance_filename         = 'NOT USED'
  character(len=500)                             :: ensemble_recenter_filename         = 'NOT USED'
  character(len=500)                             :: ensemble_member_filename_list      = 'NOT USED'
  character(len=500)                             :: ensemble_analysis_filename_list    = 'NOT USED'
  character(len=500)                             :: ensemble_recenter_geogrid_filename = 'NOT USED'
  character(len=500)                             :: ensemble_geogrid_filename          = 'NOT USED'
  character(len=500)                             :: ensemble_recenter_interp_filename  = 'NOT USED'
  character(len=500)                             :: ens_infilename                     = 'NOT USED'
  character(len=500)                             :: ens_outfilename                    = 'NOT USED'
  logical                                        :: is_ensemble_recenter               = .false.
  logical                                        :: is_ensemble_mean                   = .false.
  logical                                        :: is_ensemble_variance               = .false.
  logical                                        :: is_ensemble_inflate                = .false.
  logical                                        :: is_ensemble_update                 = .false.
  logical                                        :: debug                              = .false.
  logical                                        :: update_clouds                      = .false.
  logical                                        :: recenter_compute                   = .false.
  real(r_kind)                                   :: inflate_coeff                      = 0.0
  real(r_kind)                                   :: barnes_weights_threshold           = 1.0e-5
  real(r_kind)                                   :: barnes_distance_threshold          = 3.0e6
  integer                                        :: barnes_nneighbors                  = 10
  integer                                        :: barnes_npasses                     = 3
  integer                                        :: blndgrid_nhalo                     = 8
  integer                                        :: blndgrid_nsmooth                   = 1
  integer                                        :: nanals                             = 0
  namelist /share/                debug, is_ensemble_mean,               &
       & is_ensemble_variance, is_ensemble_recenter,                     &
       & is_ensemble_inflate, is_ensemble_update,                        &
       & ensemble_member_filename_list
  namelist /ensemble_mean_io/     ensemble_mean_filename
  namelist /ensemble_variance_io/ ensemble_variance_filename
  namelist /ensemble_recenter_io/ ensemble_recenter_filename,            &
       & ensemble_recenter_geogrid_filename, ensemble_geogrid_filename,  &
       & ensemble_recenter_interp_filename,                              &
       & barnes_nneighbors, barnes_npasses, barnes_weights_threshold,    &
       & barnes_distance_threshold, blndgrid_nhalo, blndgrid_nsmooth,    &
       & update_clouds, recenter_compute
  namelist /ensemble_inflate_io/  ensemble_analysis_filename_list,       &
       & inflate_coeff
  namelist /ensemble_update_io/   ens_infilename, ens_outfilename

  !---------------------------------------------------------------------

contains

  !=====================================================================

  ! namelistparams.f90:

  !---------------------------------------------------------------------

  subroutine namelistparams()

    ! Define variables computed within routine

    character(len=500)                                                   :: dummy
    logical                                                              :: is_it_there
    integer                                                              :: unit_nml

    ! Define counting variables

    integer                                                              :: i, j, k

    !===================================================================

    ! Define local variables

    unit_nml    = 9
    is_it_there = .false.
    inquire(file='wrf-ensemble.input',exist = is_it_there)

    ! Check local variable and proceed accordingly

    if(is_it_there) then

       ! Define local variables

       open(file   = 'wrf-ensemble.input',                               &
            unit   = unit_nml        ,                                   &
            status = 'old'         ,                                     &
            form   = 'formatted'     ,                                   &
            action = 'read'        ,                                     &
            access = 'sequential'  )
       read(unit_nml,NML = share)
       if(is_ensemble_mean)     read(unit_nml,NML=ensemble_mean_io)
       if(is_ensemble_variance) read(unit_nml,NML=ensemble_variance_io)
       if(is_ensemble_recenter) read(unit_nml,NML=ensemble_recenter_io)
       if(is_ensemble_inflate)  read(unit_nml,NML=ensemble_inflate_io)
       if(is_ensemble_update)   read(unit_nml,NML=ensemble_update_io)
       close(unit_nml)

    end if ! if(is_it_there)

    ! Check local variable and proceed accordingly

    if(.not. is_it_there) then 

       ! Print message to user

       if(mpi_procid .eq. mpi_masternode) write(6,500)

       ! Enable the root task to catch up from I/O and calculations

       call mpi_barrier(mpi_comm_world,mpi_ierror)

       ! Terminate all MPI tasks

       call mpi_interface_terminate()

    end if ! if(.not. is_it_there)

    ! Initialize local variables

    nanals = 0

    ! Define local variables

    open(99,file=trim(adjustl(ensemble_member_filename_list)),form=      &
         & 'formatted')
1000 read(99,*,end=1001) dummy
    nanals = nanals + 1
    goto 1000
1001 continue 
    close(99)

    !===================================================================

    ! Check local variable and proceed accordingly

    if(mpi_procid .eq. mpi_masternode) then

       ! Print message to user
    
       write(6,*) '&SHARE'
       write(6,*) 'DEBUG                              = ', debug
       write(6,*) 'IS_ENSEMBLE_MEAN                   = ',               &
            & is_ensemble_mean
       write(6,*) 'IS_ENSEMBLE_VARIANCE               = ',               &
            & is_ensemble_variance
       write(6,*) 'IS_ENSEMBLE_RECENTER               = ',               &
            & is_ensemble_recenter
       write(6,*)  'IS_ENSEMBLE_INFLATE               = ',               &
            & is_ensemble_inflate
       write(6,*)  'IS_ENSEMBLE_UPDATE                = ',               &
            & is_ensemble_update
       write(6,*) 'ENSEMBLE_MEMBER_FILENAME_LIST      = ',               &
            & trim(adjustl(ensemble_member_filename_list))
       write(6,*) '/'
       if(is_ensemble_mean) then
          write(6,*) '&ENSEMBLE_MEAN_IO'
          write(6,*) 'ENSEMBLE_MEAN_FILENAME             = ',            &
               & trim(adjustl(ensemble_mean_filename))
          write(6,*) '/'
       end if ! if(is_ensemble_mean)
       if(is_ensemble_variance) then
          write(6,*) '&ENSEMBLE_VARIANCE_IO'
          write(6,*) 'ENSEMBLE_VARIANCE_FILENAME         = ',            &
               & trim(adjustl(ensemble_variance_filename))
          write(6,*) '/'
       end if ! if(is_ensemble_variance)
       if(is_ensemble_recenter) then
          write(6,*) '&ENSEMBLE_RECENTER_IO'
          write(6,*) 'ENSEMBLE_RECENTER_FILENAME         = ',            &
               & trim(adjustl(ensemble_recenter_filename))
          write(6,*) 'ENSEMBLE_RECENTER_GEOGRID_FILENAME = ',            &
               & trim(adjustl(ensemble_recenter_geogrid_filename))
          write(6,*) 'ENSEMBLE_GEOGRID_FILENAME          = ',            &
               & trim(adjustl(ensemble_geogrid_filename))
          write(6,*) 'RECENTER_COMPUTE                   = ',            &
               & recenter_compute
          write(6,*) 'BARNES_NNEIGHBORS                  = ',            &
               & barnes_nneighbors
          write(6,*) 'BARNES_NPASSES                     = ',            &
               & barnes_npasses
          write(6,*) 'BARNES_WEIGHTS_THRESHOLD           = ',            &
               & barnes_weights_threshold
          write(6,*) 'BARNES_DISTANCE_THRESHOLD          = ',            &
               & barnes_distance_threshold
          write(6,*) 'BLNDGRID_NHALO                     = ',            &
               & blndgrid_nhalo
          write(6,*) 'BLNDGRID_NSMOOTH                   = ',            &
               & blndgrid_nsmooth
          write(6,*) 'update_clouds                      = ',            &
               & update_clouds
          write(6,*) '/'
       end if ! if(is_ensemble_recenter)
       if(is_ensemble_inflate) then
          write(6,*) '&ENSEMBLE_INFLATE_IO'      
          write(6,*) 'ENSEMBLE_ANALYSIS_FILENAME_LIST    = ',            &
               & trim(adjustl(ensemble_analysis_filename_list))
          write(6,*) 'INFLATE_COEFF                      = ',            &
               & inflate_coeff
          write(6,*) '/' 
       end if ! if(is_ensemble_inflate)
       if(is_ensemble_update) then
          write(6,*) '&ENSEMBLE_UPDATE_IO'      
          write(6,*) 'ENS_INFILENAME                     = ',            &
               & trim(adjustl(ens_infilename))
          write(6,*) 'ENS_OUTFILENAME                    = ',            &
               & trim(adjustl(ens_outfilename))
          write(6,*) '/' 
       end if ! if(is_ensemble_update)
       write(6,*) ' '

    end if ! if(mpi_procid .eq. mpi_masternode)

    ! Enable the root task to catch up from I/O and calculations

    call mpi_barrier(mpi_comm_world,mpi_ierror)

    !===================================================================

    ! Define format statements

500 format('NAMELISTPARAMS: wrf-ensemble.input not found in the ',       &
         & 'current working directory. ABORTING!!!!')

    !===================================================================

  end subroutine namelistparams

  !=====================================================================
  
end module namelist
