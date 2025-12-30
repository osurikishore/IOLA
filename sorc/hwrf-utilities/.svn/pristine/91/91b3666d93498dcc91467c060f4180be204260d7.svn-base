module ensemble_recenter_interface

  !=======================================================================

  ! Define associated modules and subroutines

  !-----------------------------------------------------------------------

  use constants
  use kinds

  !-----------------------------------------------------------------------

  use common_interface
  use interpolation_interface
  use mpi_interface
  use namelist
  use netcdf

  !-----------------------------------------------------------------------

  implicit none

  !-----------------------------------------------------------------------

contains

  !=======================================================================

  ! ensemble_recenter.f90:

  !-----------------------------------------------------------------------

  subroutine ensemble_recenter()

    !=====================================================================

    ! Define local variables

    call common_analysis_file_initialize()

    ! Compute local variables

    call ensemble_recenter_wrfnmm()

    ! Deallocate memory for local variables

    call common_analysis_file_cleanup()

    !=====================================================================

  end subroutine ensemble_recenter

  !=======================================================================

  ! ensemble_recenter_wrfnmm.f90:

  !-----------------------------------------------------------------------

  subroutine ensemble_recenter_wrfnmm()

    ! Define variables computed within routine

    type(wrfnmm_variables)                                               :: ens_member
    type(wrfnmm_variables)                                               :: ens_mean
    type(wrfnmm_prognostics)                                             :: ens_recenter
    type(wrfnmm_prognostics)                                             :: ens_recenter_interp

    !=====================================================================

    ! Define local variables

    ensemble_mean_filename      = trim(adjustl(analysis_file(1)))
    call common_wrfnmm_ensemble_mean_initialize(ens_mean)
    ens_member%west_east        = ens_mean%west_east
    ens_member%south_north      = ens_mean%south_north
    ens_member%bottom_top       = ens_mean%bottom_top
    ens_member%bottom_top_stag  = ens_mean%bottom_top_stag
    ens_member%soil_layers_stag = ens_mean%soil_layers_stag
    if(recenter_compute)then
       call common_wrfnmm_ensemble_recenter_initialize(ens_recenter,          &
           & ens_mean,ens_recenter_interp)
    else
       call common_wrfnmm_ensemble_recenter_interp_initialize(ens_recenter_interp)
    end if
    call common_wrfnmm_ensemble_member_initialize(ens_member)
    call common_wrfnmm_ensemble_mean_compute(ens_member,ens_mean)

    ! Compute local variables

    if(recenter_compute)then
       call ensemble_recenter_wrfnmm_compute(ens_member,ens_mean,             &
            ens_recenter,ens_recenter_interp)
    else
       call ensemble_recenter_wrfnmm_replace(ens_member,ens_mean,             &
            ens_recenter_interp)
    end if       

    ! Write variables to external file

    call ensemble_recenter_wrfnmm_write(ens_member)

    ! Deallocate memory for local variables

    call common_wrfnmm_ensemble_recenter_cleanup(ens_recenter,             &
         & ens_recenter_interp)
    call common_wrfnmm_ensemble_member_cleanup(ens_member)
    call common_wrfnmm_ensemble_member_cleanup(ens_mean)

    !=====================================================================

  end subroutine ensemble_recenter_wrfnmm

  !=======================================================================

  ! ensemble_recenter_wrfnmm_write.f90:

  !-----------------------------------------------------------------------

  subroutine ensemble_recenter_wrfnmm_write(grid)

    ! Define variables passed to routine

    type(wrfnmm_variables)                                               :: grid

    !=====================================================================

    ! Check local variable and proceed accordingly

    if(mpi_procid .ne. mpi_masternode .and. mpi_procid .le. nanals) then

       ! Define local variables

       ncstatus = nf90_open(path=trim(analysis_file(mpi_procid)),          &
            & mode=nf90_write,ncid=ncfileid)
       ncstatus = nf90_inq_varid(ncfileid,'PD',ncvarid)
       ncstatus = nf90_put_var(ncfileid,ncvarid,grid%pd)
       ncstatus = nf90_inq_varid(ncfileid,'T',ncvarid)
       ncstatus = nf90_put_var(ncfileid,ncvarid,grid%t)
       ncstatus = nf90_inq_varid(ncfileid,'Q',ncvarid)
       ncstatus = nf90_put_var(ncfileid,ncvarid,grid%q)
       ncstatus = nf90_inq_varid(ncfileid,'U',ncvarid)
       ncstatus = nf90_put_var(ncfileid,ncvarid,grid%u)
       ncstatus = nf90_inq_varid(ncfileid,'V',ncvarid)
       ncstatus = nf90_put_var(ncfileid,ncvarid,grid%v)
       ncstatus = nf90_inq_varid(ncfileid,'CWM',ncvarid)
       ncstatus = nf90_put_var(ncfileid,ncvarid,grid%cwm)
       ncstatus = nf90_close(ncfileid)

    end if ! if(mpi_procid .ne. mpi_masternode .and. mpi_procid
           ! .le. nanals)

   ! Enable the root task to catch up from I/O and calculations

    call mpi_barrier(mpi_comm_world,mpi_ierror)

    !=====================================================================

  end subroutine ensemble_recenter_wrfnmm_write

  !=======================================================================
  
  ! ensemble_recenter_wrfnmm_compute.f90:

  !-----------------------------------------------------------------------

  subroutine ensemble_recenter_wrfnmm_compute(ens_member,ens_mean,         &
       & ens_recenter,ens_recenter_interp)

    ! Define variables passed to routine

    type(wrfnmm_variables)                                               :: ens_member
    type(wrfnmm_variables)                                               :: ens_mean
    type(wrfnmm_prognostics)                                             :: ens_recenter
    type(wrfnmm_prognostics)                                             :: ens_recenter_interp

    ! Define variables computed within routine
    
    type(interpgrid)                                                     :: srcgrid
    type(interpgrid)                                                     :: dstgrid

    ! Define counting variables

    integer                                                              :: i, j, k

    !=====================================================================

    ! Initialize local variables

    call init_constants_derived()
    call init_constants(.true.)

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)
   
    call mpi_bcast(ens_mean%u,(ens_mean%west_east*ens_mean%south_north*   &
         & ens_mean%bottom_top),mpi_real,mpi_masternode,mpi_comm_world,   &
         & mpi_ierror) 
    call mpi_bcast(ens_mean%v,(ens_mean%west_east*ens_mean%south_north*   &
         & ens_mean%bottom_top),mpi_real,mpi_masternode,mpi_comm_world,   &
         & mpi_ierror) 
    call mpi_bcast(ens_mean%t,(ens_mean%west_east*ens_mean%south_north*   &
         & ens_mean%bottom_top),mpi_real,mpi_masternode,mpi_comm_world,   &
         & mpi_ierror) 
    call mpi_bcast(ens_mean%q,(ens_mean%west_east*ens_mean%south_north*   &
         & ens_mean%bottom_top),mpi_real,mpi_masternode,mpi_comm_world,   &
         & mpi_ierror) 
    call mpi_bcast(ens_mean%cwm,(ens_mean%west_east*ens_mean%south_north* &
         & ens_mean%bottom_top),mpi_real,mpi_masternode,mpi_comm_world,   &
         & mpi_ierror) 
    call mpi_bcast(ens_mean%pd,(ens_mean%west_east*ens_mean%south_north), &
         & mpi_real,mpi_masternode,mpi_comm_world,mpi_ierror)

    ! Define local variables

    srcgrid%ncoords    = ens_recenter%west_east*ens_recenter%south_north
    dstgrid%ncoords    = ens_recenter_interp%west_east*                   &
         & ens_recenter_interp%south_north
    srcgrid%neighbors  = barnes_nneighbors
    dstgrid%neighbors  = barnes_nneighbors
    srcgrid%npasses    = barnes_npasses
    dstgrid%npasses    = barnes_npasses
    srcgrid%distthresh = barnes_distance_threshold
    dstgrid%distthresh = barnes_distance_threshold
    srcgrid%xdim       = ens_recenter%west_east
    srcgrid%ydim       = ens_recenter%south_north
    dstgrid%xdim       = ens_recenter_interp%west_east
    dstgrid%ydim       = ens_recenter_interp%south_north
    call interpolation_initialize_grid(srcgrid)
    call interpolation_initialize_grid(dstgrid)
    call interpolation_initialize_task_balance(srcgrid)
    call interpolation_initialize_task_balance(dstgrid)
    srcgrid%xlong      = reshape(ens_recenter%xlong_m,                    &
         & (/srcgrid%ncoords/))*deg2rad
    dstgrid%xlong      = reshape(ens_recenter_interp%xlong_m,             &
         & (/dstgrid%ncoords/))*deg2rad
    srcgrid%xlat       = reshape(ens_recenter%xlat_m,                     &
         & (/srcgrid%ncoords/))*deg2rad
    dstgrid%xlat       = reshape(ens_recenter_interp%xlat_m,              &
         & (/dstgrid%ncoords/))*deg2rad

    ! Compute local variables

    call interpolation_define_kdtree_mpi(srcgrid)
    call interpolation_define_kdtree_mpi(dstgrid)
    call interpolation_define_weights_mpi(srcgrid,dstgrid)

    ! Allocate memory for local variables
       
    if(.not. allocated(srcgrid_var))                                      &
         & allocate(srcgrid_var(srcgrid%ncoords))
    if(.not. allocated(dstgrid_var))                                      &
         & allocate(dstgrid_var(dstgrid%ncoords))
    if(.not. allocated(blndgrid_a))                                       &
         & allocate(blndgrid_a(dstgrid%ncoords))
    if(.not. allocated(blndgrid_b))                                       &
         & allocate(blndgrid_b(dstgrid%ncoords))

    ! Define local variables

    srcgrid_var = reshape(ens_recenter%pd,(/srcgrid%ncoords/))

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)
    
    call mpi_bcast(srcgrid_var,srcgrid%ncoords,mpi_real,                  &
         & mpi_masternode,mpi_comm_world,mpi_ierror)

    ! Compute local variables

    call interpolation_barnes_analysis_mpi(srcgrid,dstgrid)

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)
    
    call mpi_bcast(dstgrid_var,dstgrid%ncoords,mpi_real,                  &
         & mpi_masternode,mpi_comm_world,mpi_ierror)

    ! Define local variables

    where(dstgrid%weights(:,1,1) .eq. 0) dstgrid_var = -1.e30

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)
    
    call mpi_bcast(dstgrid_var,dstgrid%ncoords,mpi_real,                  &
         & mpi_masternode,mpi_comm_world,mpi_ierror)

    ! Define local variables

    blndgrid_a = dstgrid_var

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)
    
    call mpi_bcast(blndgrid_a,dstgrid%ncoords,mpi_real,                   &
         & mpi_masternode,mpi_comm_world,mpi_ierror)

    ! Define local variables

    where(dstgrid%weights(:,1,1) .eq. 0)                                  &
         & blndgrid_a = reshape(ens_mean%pd,(/dstgrid%ncoords/))

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)
    
    call mpi_bcast(blndgrid_a,dstgrid%ncoords,mpi_real,                   &
         & mpi_masternode,mpi_comm_world,mpi_ierror)

    ! Compute local variables

    call interpolation_blend_analysis(dstgrid)

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)
    
    call mpi_bcast(dstgrid_var,dstgrid%ncoords,mpi_real,                  &
         & mpi_masternode,mpi_comm_world,mpi_ierror)

    ! Define local variables

    ens_recenter_interp%pd = reshape(dstgrid_var,                         &
         & (/ens_recenter_interp%west_east,                               &
         & ens_recenter_interp%south_north/))

    ! Loop through local variable

    do k = 1, ens_mean%bottom_top

       ! Define local variables

       srcgrid_var = reshape(ens_recenter%t(:,:,k),(/srcgrid%ncoords/))

       ! Broadcast all necessary variables to slave (compute) nodes
       ! (tasks)
    
       call mpi_bcast(srcgrid_var,srcgrid%ncoords,mpi_real,               &
            & mpi_masternode,mpi_comm_world,mpi_ierror)

       ! Compute local variables

       call interpolation_barnes_analysis_mpi(srcgrid,dstgrid)

       ! Broadcast all necessary variables to slave (compute) nodes
       ! (tasks)
    
       call mpi_bcast(dstgrid_var,dstgrid%ncoords,mpi_real,               &
            & mpi_masternode,mpi_comm_world,mpi_ierror)

       ! Define local variables

       blndgrid_a = dstgrid_var

       ! Broadcast all necessary variables to slave (compute) nodes
       ! (tasks)
    
       call mpi_bcast(blndgrid_a,dstgrid%ncoords,mpi_real,                &
            & mpi_masternode,mpi_comm_world,mpi_ierror)

       ! Define local variables

       where(dstgrid%weights(:,1,1) .eq. 0)                               &
            & blndgrid_a = reshape(ens_mean%t(:,:,k),                     &
            & (/dstgrid%ncoords/))

       ! Broadcast all necessary variables to slave (compute) nodes
       ! (tasks)
    
       call mpi_bcast(blndgrid_a,dstgrid%ncoords,mpi_real,                &
            & mpi_masternode,mpi_comm_world,mpi_ierror)

       ! Compute local variables

       call interpolation_blend_analysis(dstgrid)

       ! Broadcast all necessary variables to slave (compute) nodes
       ! (tasks)
    
       call mpi_bcast(dstgrid_var,dstgrid%ncoords,mpi_real,               &
            & mpi_masternode,mpi_comm_world,mpi_ierror)

       ! Define local variables

       ens_recenter_interp%t(:,:,k) = reshape(dstgrid_var,                &
            & (/ens_recenter_interp%west_east,                            &
            & ens_recenter_interp%south_north/))    

       ! Define local variables

       srcgrid_var = reshape(ens_recenter%q(:,:,k),(/srcgrid%ncoords/))

       ! Broadcast all necessary variables to slave (compute) nodes
       ! (tasks)
    
       call mpi_bcast(srcgrid_var,srcgrid%ncoords,mpi_real,               &
            & mpi_masternode,mpi_comm_world,mpi_ierror)

       ! Compute local variables

       call interpolation_barnes_analysis_mpi(srcgrid,dstgrid)

       ! Broadcast all necessary variables to slave (compute) nodes
       ! (tasks)
    
       call mpi_bcast(dstgrid_var,dstgrid%ncoords,mpi_real,               &
            & mpi_masternode,mpi_comm_world,mpi_ierror)

       ! Define local variables

       blndgrid_a = dstgrid_var

       ! Broadcast all necessary variables to slave (compute) nodes
       ! (tasks)
    
       call mpi_bcast(blndgrid_a,dstgrid%ncoords,mpi_real,                &
            & mpi_masternode,mpi_comm_world,mpi_ierror)

       ! Define local variables

       where(dstgrid%weights(:,1,1) .eq. 0)                               &
            & blndgrid_a = reshape(ens_mean%q(:,:,k),                     &
            & (/dstgrid%ncoords/))

       ! Broadcast all necessary variables to slave (compute) nodes
       ! (tasks)
    
       call mpi_bcast(blndgrid_a,dstgrid%ncoords,mpi_real,                &
            & mpi_masternode,mpi_comm_world,mpi_ierror)

       ! Compute local variables

       call interpolation_blend_analysis(dstgrid)

       ! Broadcast all necessary variables to slave (compute) nodes
       ! (tasks)
    
       call mpi_bcast(dstgrid_var,dstgrid%ncoords,mpi_real,               &
            & mpi_masternode,mpi_comm_world,mpi_ierror)

       ! Define local variables

       ens_recenter_interp%q(:,:,k) = reshape(dstgrid_var,                &
            & (/ens_recenter_interp%west_east,                            &
            & ens_recenter_interp%south_north/))   

       ! Define local variables

       srcgrid_var = reshape(ens_recenter%cwm(:,:,k),                     &
            & (/srcgrid%ncoords/))

       ! Broadcast all necessary variables to slave (compute) nodes
       ! (tasks)
    
       call mpi_bcast(srcgrid_var,srcgrid%ncoords,mpi_real,               &
            & mpi_masternode,mpi_comm_world,mpi_ierror)

       ! Compute local variables

       call interpolation_barnes_analysis_mpi(srcgrid,dstgrid)

       ! Broadcast all necessary variables to slave (compute) nodes
       ! (tasks)
    
       call mpi_bcast(dstgrid_var,dstgrid%ncoords,mpi_real,               &
            & mpi_masternode,mpi_comm_world,mpi_ierror)

       ! Define local variables

       blndgrid_a = dstgrid_var

       ! Broadcast all necessary variables to slave (compute) nodes
       ! (tasks)
    
       call mpi_bcast(blndgrid_a,dstgrid%ncoords,mpi_real,                &
            & mpi_masternode,mpi_comm_world,mpi_ierror)

       ! Define local variables

       where(dstgrid%weights(:,1,1) .eq. 0)                               &
            & blndgrid_a = reshape(ens_mean%cwm(:,:,k),                   &
            & (/dstgrid%ncoords/))

       ! Broadcast all necessary variables to slave (compute) nodes
       ! (tasks)
    
       call mpi_bcast(blndgrid_a,dstgrid%ncoords,mpi_real,                &
            & mpi_masternode,mpi_comm_world,mpi_ierror)

       ! Compute local variables

       call interpolation_blend_analysis(dstgrid)

       ! Broadcast all necessary variables to slave (compute) nodes
       ! (tasks)
    
       call mpi_bcast(dstgrid_var,dstgrid%ncoords,mpi_real,               &
            & mpi_masternode,mpi_comm_world,mpi_ierror)

       ! Define local variables

       ens_recenter_interp%cwm(:,:,k) = reshape(dstgrid_var,              &
            & (/ens_recenter_interp%west_east,                            &
            & ens_recenter_interp%south_north/))   

    end do !  do k = 1, ens_mean%bottom_top

    ! Deallocate memory for local variables
       
    if(allocated(srcgrid_var)) deallocate(srcgrid_var)
    if(allocated(dstgrid_var)) deallocate(dstgrid_var)
    call interpolation_cleanup_grid(srcgrid)
    call interpolation_cleanup_grid(dstgrid)
    call interpolation_cleanup_task_balance(srcgrid)
    call interpolation_cleanup_task_balance(dstgrid)

    ! Define local variables

    srcgrid%ncoords    = ens_recenter%west_east*ens_recenter%south_north
    dstgrid%ncoords    = ens_recenter_interp%west_east*                   &
         & ens_recenter_interp%south_north
    srcgrid%neighbors  = barnes_nneighbors
    dstgrid%neighbors  = barnes_nneighbors
    srcgrid%npasses    = barnes_npasses
    dstgrid%npasses    = barnes_npasses
    srcgrid%distthresh = barnes_distance_threshold
    dstgrid%distthresh = barnes_distance_threshold
    srcgrid%xdim       = ens_recenter%west_east
    srcgrid%ydim       = ens_recenter%south_north
    dstgrid%xdim       = ens_recenter_interp%west_east
    dstgrid%ydim       = ens_recenter_interp%south_north
    call interpolation_initialize_grid(srcgrid)
    call interpolation_initialize_grid(dstgrid)
    call interpolation_initialize_task_balance(srcgrid)
    call interpolation_initialize_task_balance(dstgrid)
    srcgrid%xlong      = reshape(ens_recenter%xlong_v,                    &
         & (/srcgrid%ncoords/))*deg2rad
    dstgrid%xlong      = reshape(ens_recenter_interp%xlong_v,             &
         & (/dstgrid%ncoords/))*deg2rad
    srcgrid%xlat       = reshape(ens_recenter%xlat_v,                     &
         & (/srcgrid%ncoords/))*deg2rad
    dstgrid%xlat       = reshape(ens_recenter_interp%xlat_v,              &
         & (/dstgrid%ncoords/))*deg2rad

    ! Compute local variables

    call interpolation_define_kdtree_mpi(srcgrid)
    call interpolation_define_kdtree_mpi(dstgrid)
    call interpolation_define_weights_mpi(srcgrid,dstgrid)

    ! Allocate memory for local variables
       
    if(.not. allocated(srcgrid_var))                                      &
         & allocate(srcgrid_var(srcgrid%ncoords))
    if(.not. allocated(dstgrid_var))                                      &
         & allocate(dstgrid_var(dstgrid%ncoords))

    ! Loop through local variable

    do k = 1, ens_mean%bottom_top

       ! Define local variables

       srcgrid_var = reshape(ens_recenter%u(:,:,k),(/srcgrid%ncoords/))

       ! Broadcast all necessary variables to slave (compute) nodes
       ! (tasks)
    
       call mpi_bcast(srcgrid_var,srcgrid%ncoords,mpi_real,               &
            & mpi_masternode,mpi_comm_world,mpi_ierror)

       ! Compute local variables

       call interpolation_barnes_analysis_mpi(srcgrid,dstgrid)

       ! Broadcast all necessary variables to slave (compute) nodes
       ! (tasks)
    
       call mpi_bcast(dstgrid_var,dstgrid%ncoords,mpi_real,               &
            & mpi_masternode,mpi_comm_world,mpi_ierror)

       ! Define local variables

       blndgrid_a = dstgrid_var

       ! Broadcast all necessary variables to slave (compute) nodes
       ! (tasks)
    
       call mpi_bcast(blndgrid_a,dstgrid%ncoords,mpi_real,                &
            & mpi_masternode,mpi_comm_world,mpi_ierror)

       ! Define local variables

       where(dstgrid%weights(:,1,1) .eq. 0)                               &
            & blndgrid_a = reshape(ens_mean%u(:,:,k),                     &
            & (/dstgrid%ncoords/))

       ! Broadcast all necessary variables to slave (compute) nodes
       ! (tasks)
    
       call mpi_bcast(blndgrid_a,dstgrid%ncoords,mpi_real,                &
            & mpi_masternode,mpi_comm_world,mpi_ierror)

       ! Compute local variables

       call interpolation_blend_analysis(dstgrid)

       ! Broadcast all necessary variables to slave (compute) nodes
       ! (tasks)
    
       call mpi_bcast(dstgrid_var,dstgrid%ncoords,mpi_real,               &
            & mpi_masternode,mpi_comm_world,mpi_ierror)

       ! Define local variables

       ens_recenter_interp%u(:,:,k) = reshape(dstgrid_var,                &
            & (/ens_recenter_interp%west_east,                            &
            & ens_recenter_interp%south_north/))   

       ! Define local variables

       srcgrid_var = reshape(ens_recenter%v(:,:,k),(/srcgrid%ncoords/))

       ! Broadcast all necessary variables to slave (compute) nodes
       ! (tasks)
    
       call mpi_bcast(srcgrid_var,srcgrid%ncoords,mpi_real,               &
            & mpi_masternode,mpi_comm_world,mpi_ierror)

       ! Compute local variables

       call interpolation_barnes_analysis_mpi(srcgrid,dstgrid)

       ! Broadcast all necessary variables to slave (compute) nodes
       ! (tasks)
    
       call mpi_bcast(dstgrid_var,dstgrid%ncoords,mpi_real,               &
            & mpi_masternode,mpi_comm_world,mpi_ierror)

       ! Define local variables

       blndgrid_a = dstgrid_var

       ! Broadcast all necessary variables to slave (compute) nodes
       ! (tasks)
    
       call mpi_bcast(blndgrid_a,dstgrid%ncoords,mpi_real,                &
            & mpi_masternode,mpi_comm_world,mpi_ierror)

       ! Define local variables

       where(dstgrid%weights(:,1,1) .eq. 0)                               &
            & blndgrid_a = reshape(ens_mean%v(:,:,k),                     &
            & (/dstgrid%ncoords/))

       ! Broadcast all necessary variables to slave (compute) nodes
       ! (tasks)
    
       call mpi_bcast(blndgrid_a,dstgrid%ncoords,mpi_real,                &
            & mpi_masternode,mpi_comm_world,mpi_ierror)

       ! Compute local variables

       call interpolation_blend_analysis(dstgrid)

       ! Broadcast all necessary variables to slave (compute) nodes
       ! (tasks)
    
       call mpi_bcast(dstgrid_var,dstgrid%ncoords,mpi_real,               &
            & mpi_masternode,mpi_comm_world,mpi_ierror)

       ! Define local variables

       ens_recenter_interp%v(:,:,k) = reshape(dstgrid_var,                &
            & (/ens_recenter_interp%west_east,                            &
            & ens_recenter_interp%south_north/))   

    end do ! do k = 1, ens_mean%bottom_top

    ! Check local variable and proceed accordingly

    if(mpi_procid .ne. mpi_masternode .and. mpi_procid .le. nanals) then

       ! Compute local variables

       ens_member%pd  = (ens_member%pd  - ens_mean%pd)  +                 &
            & ens_recenter_interp%pd
       ens_member%u   = (ens_member%u   - ens_mean%u)   +                 &
            & ens_recenter_interp%u
       ens_member%v   = (ens_member%v   - ens_mean%v)   +                 &
            & ens_recenter_interp%v
       ens_member%t   = (ens_member%t   - ens_mean%t)   +                 &
            & ens_recenter_interp%t
       ens_member%q   = (ens_member%q   - ens_mean%q)   +                 &
            & ens_recenter_interp%q
       ens_member%cwm = (ens_member%cwm - ens_mean%cwm) +                 &
            & ens_recenter_interp%cwm

    end if ! if(mpi_procid .ne. mpi_masternode .and. mpi_procid
           ! .le. nanals)

    ! Enable the root task to catch up from I/O and calculations

    call mpi_barrier(mpi_comm_world,mpi_ierror)

    ! Deallocate memory for local variables
       
    if(allocated(srcgrid_var)) deallocate(srcgrid_var)
    if(allocated(dstgrid_var)) deallocate(dstgrid_var)
    call interpolation_cleanup_grid(srcgrid)
    call interpolation_cleanup_grid(dstgrid)
    call interpolation_cleanup_task_balance(srcgrid)
    call interpolation_cleanup_task_balance(dstgrid)

    !=====================================================================

  end subroutine ensemble_recenter_wrfnmm_compute

  !=======================================================================

  ! ensemble_recenter_wrfnmm_compute.f90:

  !-----------------------------------------------------------------------

  subroutine ensemble_recenter_wrfnmm_replace(ens_member,ens_mean,         &
       & ens_recenter_interp)

    ! Define variables passed to routine

    type(wrfnmm_variables)                                               :: ens_member
    type(wrfnmm_variables)                                               :: ens_mean
    type(wrfnmm_prognostics)                                             :: ens_recenter_interp

    integer                                                              :: im, jm, km

    !=====================================================================

    ! Initialize local variables

    call init_constants_derived()
    call init_constants(.true.)

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)

    im=ens_mean%west_east
    jm=ens_mean%south_north
    km=ens_mean%bottom_top
    call mpi_bcast(ens_mean%u,(im*jm*km),   &
         & mpi_real,mpi_masternode,mpi_comm_world,   &
         & mpi_ierror)
    call mpi_bcast(ens_mean%v,(im*jm*km),   &
         & mpi_real,mpi_masternode,mpi_comm_world,   &
         & mpi_ierror)
    call mpi_bcast(ens_mean%t,(im*jm*km),   &
         & mpi_real,mpi_masternode,mpi_comm_world,   &
         & mpi_ierror)
    call mpi_bcast(ens_mean%q,(im*jm*km),   &
         & mpi_real,mpi_masternode,mpi_comm_world,   &
         & mpi_ierror)
    call mpi_bcast(ens_mean%cwm,(im*jm*km), &
         & mpi_real,mpi_masternode,mpi_comm_world,   &
         & mpi_ierror)
    call mpi_bcast(ens_mean%pd,(im*jm), &
         & mpi_real,mpi_masternode,mpi_comm_world,mpi_ierror)


    im=ens_recenter_interp%west_east
    jm=ens_recenter_interp%south_north
    km=ens_recenter_interp%bottom_top
    call mpi_bcast(ens_recenter_interp%u,(im*jm*km),   &
         & mpi_real,mpi_masternode,mpi_comm_world,   &
         & mpi_ierror)
    call mpi_bcast(ens_recenter_interp%v,(im*jm*km),   &
         & mpi_real,mpi_masternode,mpi_comm_world,   &
         & mpi_ierror)
    call mpi_bcast(ens_recenter_interp%t,(im*jm*km),   &
         & mpi_real,mpi_masternode,mpi_comm_world,   &
         & mpi_ierror)
    call mpi_bcast(ens_recenter_interp%q,(im*jm*km),   &
         & mpi_real,mpi_masternode,mpi_comm_world,   &
         & mpi_ierror)
    call mpi_bcast(ens_recenter_interp%cwm,(im*jm*km), &
         & mpi_real,mpi_masternode,mpi_comm_world,   &
         & mpi_ierror)
    call mpi_bcast(ens_recenter_interp%pd,(im*jm), &
         & mpi_real,mpi_masternode,mpi_comm_world,mpi_ierror)


    if(mpi_procid .ne. mpi_masternode .and. mpi_procid .le. nanals) then

       ! Compute local variables

       ens_member%pd  = (ens_member%pd  - ens_mean%pd)  +                 &
            & ens_recenter_interp%pd
       ens_member%u   = (ens_member%u   - ens_mean%u)   +                 &
            & ens_recenter_interp%u
       ens_member%v   = (ens_member%v   - ens_mean%v)   +                 &
            & ens_recenter_interp%v
       ens_member%t   = (ens_member%t   - ens_mean%t)   +                 &
            & ens_recenter_interp%t
       ens_member%q   = (ens_member%q   - ens_mean%q)   +                 &
            & ens_recenter_interp%q
       if(update_clouds)then
          ens_member%cwm = (ens_member%cwm - ens_mean%cwm) +              &
               & ens_recenter_interp%cwm
       end if

    end if ! if(mpi_procid .ne. mpi_masternode .and. mpi_procid
           ! .le. nanals)

    ! Enable the root task to catch up from I/O and calculations

    call mpi_barrier(mpi_comm_world,mpi_ierror)

    !=====================================================================

  end subroutine ensemble_recenter_wrfnmm_replace

  !=======================================================================

end module ensemble_recenter_interface
