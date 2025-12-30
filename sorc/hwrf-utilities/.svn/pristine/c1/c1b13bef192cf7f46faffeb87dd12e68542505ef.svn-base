module common_interface

  !=======================================================================

  ! Define associated modules and subroutines

  !-----------------------------------------------------------------------

  use constants
  use kinds

  !-----------------------------------------------------------------------

  use mpi_interface
  use namelist
  use netcdf

  !-----------------------------------------------------------------------

  implicit none

  !-----------------------------------------------------------------------

  ! Define all data and structure types for routine; these variables
  ! are variables required by the subroutines within this module

  type wrfnmm_variables
     real(r_kind),                 dimension(:,:,:),   allocatable :: pint
     real(r_kind),                 dimension(:,:,:),   allocatable :: t
     real(r_kind),                 dimension(:,:,:),   allocatable :: q
     real(r_kind),                 dimension(:,:,:),   allocatable :: u
     real(r_kind),                 dimension(:,:,:),   allocatable :: v
     real(r_kind),                 dimension(:,:,:),   allocatable :: cwm
     real(r_kind),                 dimension(:,:,:),   allocatable :: f_ice
     real(r_kind),                 dimension(:,:,:),   allocatable :: f_rain
     real(r_kind),                 dimension(:,:,:),   allocatable :: f_rimef
     real(r_kind),                 dimension(:,:,:),   allocatable :: smc
     real(r_kind),                 dimension(:,:,:),   allocatable :: stc
     real(r_kind),                 dimension(:,:),     allocatable :: pd
     real(r_kind),                 dimension(:,:),     allocatable :: tsk
     real(r_kind),                 dimension(:,:),     allocatable :: sst
     real(r_kind),                 dimension(:,:),     allocatable :: sno
     real(r_kind),                 dimension(:,:),     allocatable :: sice
     real(r_kind),                 dimension(:,:),     allocatable :: u10
     real(r_kind),                 dimension(:,:),     allocatable :: v10
     real(r_kind),                 dimension(:,:),     allocatable :: xlong_m
     real(r_kind),                 dimension(:,:),     allocatable :: xlat_m
     real(r_kind),                 dimension(:,:),     allocatable :: xlong_v
     real(r_kind),                 dimension(:,:),     allocatable :: xlat_v
     integer                                                       :: west_east
     integer                                                       :: south_north
     integer                                                       :: bottom_top
     integer                                                       :: bottom_top_stag
     integer                                                       :: soil_layers_stag
  end type wrfnmm_variables

  type wrfnmm_prognostics
     logical                                                       :: is_exist_pint
     logical                                                       :: is_exist_t
     logical                                                       :: is_exist_q
     logical                                                       :: is_exist_u
     logical                                                       :: is_exist_v
     logical                                                       :: is_exist_cwm
     logical                                                       :: is_exist_pd
     real(r_kind),                 dimension(:,:,:),   allocatable :: pint
     real(r_kind),                 dimension(:,:,:),   allocatable :: t
     real(r_kind),                 dimension(:,:,:),   allocatable :: q
     real(r_kind),                 dimension(:,:,:),   allocatable :: u
     real(r_kind),                 dimension(:,:,:),   allocatable :: v
     real(r_kind),                 dimension(:,:,:),   allocatable :: cwm
     real(r_kind),                 dimension(:,:),     allocatable :: pd
     real(r_kind),                 dimension(:,:),     allocatable :: xlong_m
     real(r_kind),                 dimension(:,:),     allocatable :: xlat_m
     real(r_kind),                 dimension(:,:),     allocatable :: xlong_v
     real(r_kind),                 dimension(:,:),     allocatable :: xlat_v
     integer                                                       :: west_east
     integer                                                       :: south_north
     integer                                                       :: bottom_top
     integer                                                       :: bottom_top_stag
  end type wrfnmm_prognostics

  ! Define global variables

  character(len=500),              dimension(:),       allocatable :: analysis_file
  integer                                                          :: ncfileid
  integer                                                          :: ncvarid
  integer                                                          :: ncdimid
  integer                                                          :: ncstatus

  !-----------------------------------------------------------------------

contains

  !=======================================================================

  ! common_wrfnmm_ensemble_mean_compute.f90:

  !-----------------------------------------------------------------------

  subroutine common_wrfnmm_ensemble_mean_compute(ens_member,ens_mean)

    ! Define variables computed within routine

    type(wrfnmm_variables)                                               :: ens_member
    type(wrfnmm_variables)                                               :: ens_mean

    !=====================================================================

    ! Compute local variables

    call mpi_reduce(ens_member%pint,ens_mean%pint,(ens_mean%west_east*     &
         & ens_mean%south_north*ens_mean%bottom_top_stag),mpi_real,        &
         & mpi_sum,mpi_masternode,mpi_comm_world,mpi_ierror)
    call mpi_reduce(ens_member%t,ens_mean%t,(ens_mean%west_east*           &
         & ens_mean%south_north*ens_mean%bottom_top),mpi_real,mpi_sum,     &
         & mpi_masternode,mpi_comm_world,mpi_ierror)
    call mpi_reduce(ens_member%q,ens_mean%q,(ens_mean%west_east*           &
         & ens_mean%south_north*ens_mean%bottom_top),mpi_real,mpi_sum,     &
         & mpi_masternode,mpi_comm_world,mpi_ierror)
    call mpi_reduce(ens_member%u,ens_mean%u,(ens_mean%west_east*           &
         & ens_mean%south_north*ens_mean%bottom_top),mpi_real,mpi_sum,     &
         & mpi_masternode,mpi_comm_world,mpi_ierror)
    call mpi_reduce(ens_member%v,ens_mean%v,(ens_mean%west_east*           &
         & ens_mean%south_north*ens_mean%bottom_top),mpi_real,mpi_sum,     &
         & mpi_masternode,mpi_comm_world,mpi_ierror)
    call mpi_reduce(ens_member%cwm,ens_mean%cwm,(ens_mean%west_east*       &
         & ens_mean%south_north*ens_mean%bottom_top),mpi_real,mpi_sum,     &
         & mpi_masternode,mpi_comm_world,mpi_ierror)
    call mpi_reduce(ens_member%f_ice,ens_mean%f_ice,(ens_mean%west_east*   &
         & ens_mean%south_north*ens_mean%bottom_top),mpi_real,mpi_sum,     &
         & mpi_masternode,mpi_comm_world,mpi_ierror)
    call mpi_reduce(ens_member%f_rain,ens_mean%f_rain,(ens_mean%west_east* &
         & ens_mean%south_north*ens_mean%bottom_top),mpi_real,mpi_sum,     &
         & mpi_masternode,mpi_comm_world,mpi_ierror)
    call mpi_reduce(ens_member%f_rimef,ens_mean%f_rimef,                   &
         & (ens_mean%west_east*ens_mean%south_north*ens_mean%bottom_top),  &
         & mpi_real,mpi_sum,mpi_masternode,mpi_comm_world,mpi_ierror)
    call mpi_reduce(ens_member%smc,ens_mean%smc,(ens_mean%west_east*       &
         & ens_mean%south_north*ens_mean%soil_layers_stag),mpi_real,       &
         & mpi_sum,mpi_masternode,mpi_comm_world,mpi_ierror)
    call mpi_reduce(ens_member%stc,ens_mean%stc,(ens_mean%west_east*       &
         & ens_mean%south_north*ens_mean%soil_layers_stag),mpi_real,       &
         & mpi_sum,mpi_masternode,mpi_comm_world,mpi_ierror)
    call mpi_reduce(ens_member%pd,ens_mean%pd,(ens_mean%west_east*         &
         & ens_mean%south_north),mpi_real,mpi_sum,mpi_masternode,          &
         & mpi_comm_world,mpi_ierror)
    call mpi_reduce(ens_member%tsk,ens_mean%tsk,(ens_mean%west_east*       &
         & ens_mean%south_north),mpi_real,mpi_sum,mpi_masternode,          &
         & mpi_comm_world,mpi_ierror)
    call mpi_reduce(ens_member%sst,ens_mean%sst,(ens_mean%west_east*       &
         & ens_mean%south_north),mpi_real,mpi_sum,mpi_masternode,          &
         & mpi_comm_world,mpi_ierror)
    call mpi_reduce(ens_member%sno,ens_mean%sno,(ens_mean%west_east*       &
         & ens_mean%south_north),mpi_real,mpi_sum,mpi_masternode,          &
         & mpi_comm_world,mpi_ierror)
    call mpi_reduce(ens_member%sice,ens_mean%sice,(ens_mean%west_east*     &
         & ens_mean%south_north),mpi_real,mpi_sum,mpi_masternode,          &
         & mpi_comm_world,mpi_ierror)
    call mpi_reduce(ens_member%u10,ens_mean%u10,(ens_mean%west_east*       &
         & ens_mean%south_north),mpi_real,mpi_sum,mpi_masternode,          &
         & mpi_comm_world,mpi_ierror)
    call mpi_reduce(ens_member%v10,ens_mean%v10,(ens_mean%west_east*       &
         & ens_mean%south_north),mpi_real,mpi_sum,mpi_masternode,          &
         & mpi_comm_world,mpi_ierror)
    ens_mean%pint    = ens_mean%pint   /real(nanals)
    ens_mean%t       = ens_mean%t      /real(nanals)
    ens_mean%q       = ens_mean%q      /real(nanals)
    ens_mean%u       = ens_mean%u      /real(nanals)
    ens_mean%v       = ens_mean%v      /real(nanals)
    ens_mean%cwm     = ens_mean%cwm    /real(nanals)
    ens_mean%f_ice   = ens_mean%f_ice  /real(nanals)
    ens_mean%f_rain  = ens_mean%f_rain /real(nanals)
    ens_mean%f_rimef = ens_mean%f_rimef/real(nanals)
    ens_mean%smc     = ens_mean%smc    /real(nanals)
    ens_mean%stc     = ens_mean%stc    /real(nanals)
    ens_mean%pd      = ens_mean%pd     /real(nanals)
    ens_mean%tsk     = ens_mean%tsk    /real(nanals)
    ens_mean%sst     = ens_mean%sst    /real(nanals)
    ens_mean%sno     = ens_mean%sno    /real(nanals)
    ens_mean%sice    = ens_mean%sice   /real(nanals)
    ens_mean%u10     = ens_mean%u10    /real(nanals)
    ens_mean%v10     = ens_mean%v10    /real(nanals)

    ! Check local variables and proceed accordingly

    where(ens_mean%q       .lt. 0.0) ens_mean%q       = 0.0
    where(ens_mean%cwm     .lt. 0.0) ens_mean%cwm     = 0.0
    where(ens_mean%f_ice   .lt. 0.0) ens_mean%f_ice   = 0.0
    where(ens_mean%f_rain  .lt. 0.0) ens_mean%f_rain  = 0.0
    where(ens_mean%f_rimef .lt. 0.0) ens_mean%f_rimef = 0.0
    where(ens_mean%sice    .lt. 0.0) ens_mean%sice    = 0.0
    where(ens_mean%sice    .gt. 1.0) ens_mean%sice    = 1.0

    !=====================================================================

  end subroutine common_wrfnmm_ensemble_mean_compute

  !=======================================================================

  ! common_wrfnmm_ensemble_variance_compute.f90:

  !-----------------------------------------------------------------------

  subroutine common_wrfnmm_ensemble_variance_compute(ens_member,ens_mean,  &
       & ens_vari)

    ! Define variables computed within routine

    type(wrfnmm_variables)                                               :: ens_member
    type(wrfnmm_variables)                                               :: ens_mean
    type(wrfnmm_variables)                                               :: ens_vari

    ! Define counting variables

    integer                                                              :: i, j, k

    !=====================================================================

    ! Compute local variables

    call common_wrfnmm_ensemble_mean_compute(ens_member,ens_mean)

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)
    
    call mpi_bcast(ens_mean%pint,(ens_mean%west_east*                      &
         & ens_mean%south_north*ens_mean%bottom_top_stag),mpi_real,        &
         & mpi_masternode,mpi_comm_world,mpi_ierror)
    call mpi_bcast(ens_mean%t,(ens_mean%west_east*ens_mean%south_north*    &
         & ens_mean%bottom_top),mpi_real,mpi_masternode,mpi_comm_world,    &
         & mpi_ierror)
    call mpi_bcast(ens_mean%q,(ens_mean%west_east*ens_mean%south_north*    &
         & ens_mean%bottom_top),mpi_real,mpi_masternode,mpi_comm_world,    &
         & mpi_ierror)
    call mpi_bcast(ens_mean%u,(ens_mean%west_east*ens_mean%south_north*    &
         & ens_mean%bottom_top),mpi_real,mpi_masternode,mpi_comm_world,    &
         & mpi_ierror)
    call mpi_bcast(ens_mean%v,(ens_mean%west_east*ens_mean%south_north*    &
         & ens_mean%bottom_top),mpi_real,mpi_masternode,mpi_comm_world,    &
         & mpi_ierror)
    call mpi_bcast(ens_mean%cwm,(ens_mean%west_east*ens_mean%south_north*  &
         & ens_mean%bottom_top),mpi_real,mpi_masternode,mpi_comm_world,    &
         & mpi_ierror)
    call mpi_bcast(ens_mean%f_ice,(ens_mean%west_east*                     &
         & ens_mean%south_north*ens_mean%bottom_top),mpi_real,             &
         & mpi_masternode,mpi_comm_world,mpi_ierror)
    call mpi_bcast(ens_mean%f_rain,(ens_mean%west_east*                    &
         & ens_mean%south_north*ens_mean%bottom_top),mpi_real,             &
         & mpi_masternode,mpi_comm_world,mpi_ierror)
    call mpi_bcast(ens_mean%f_rimef,(ens_mean%west_east*                   &
         & ens_mean%south_north*ens_mean%bottom_top),mpi_real,             &
         & mpi_masternode,mpi_comm_world,mpi_ierror)
    call mpi_bcast(ens_mean%smc,(ens_mean%west_east*                       &
         & ens_mean%south_north),mpi_real,mpi_masternode,mpi_comm_world,   &
         & mpi_ierror)
    call mpi_bcast(ens_mean%stc,(ens_mean%west_east*                       &
         & ens_mean%south_north),mpi_real,mpi_masternode,mpi_comm_world,   &
         & mpi_ierror)
    call mpi_bcast(ens_mean%pd,(ens_mean%west_east*                        &
         & ens_mean%south_north),mpi_real,mpi_masternode,mpi_comm_world,   &
         & mpi_ierror)
    call mpi_bcast(ens_mean%tsk,(ens_mean%west_east*                       &
         & ens_mean%south_north),mpi_real,mpi_masternode,mpi_comm_world,   &
         & mpi_ierror)
    call mpi_bcast(ens_mean%sst,(ens_mean%west_east*                       &
         & ens_mean%south_north),mpi_real,mpi_masternode,mpi_comm_world,   &
         & mpi_ierror)
    call mpi_bcast(ens_mean%sno,(ens_mean%west_east*                       &
         & ens_mean%south_north),mpi_real,mpi_masternode,mpi_comm_world,   &
         & mpi_ierror)
    call mpi_bcast(ens_mean%sice,(ens_mean%west_east*                      &
         & ens_mean%south_north),mpi_real,mpi_masternode,mpi_comm_world,   &
         & mpi_ierror)
    call mpi_bcast(ens_mean%u10,(ens_mean%west_east*                       &
         & ens_mean%south_north),mpi_real,mpi_masternode,mpi_comm_world,   &
         & mpi_ierror)
    call mpi_bcast(ens_mean%v10,(ens_mean%west_east*                       &
         & ens_mean%south_north),mpi_real,mpi_masternode,mpi_comm_world,   &
         & mpi_ierror)

    ! Check local variable and proceed accordingly

    if(mpi_procid .ne. mpi_masternode .and. mpi_procid .le. nanals) then

       ! Compute local variables

       ens_member%pint    = (ens_member%pint    - ens_mean%pint)**2.0
       ens_member%t       = (ens_member%t       - ens_mean%t)**2.0
       ens_member%q       = (ens_member%q       - ens_mean%q)**2.0
       ens_member%u       = (ens_member%u       - ens_mean%u)**2.0
       ens_member%v       = (ens_member%v       - ens_mean%v)**2.0
       ens_member%cwm     = (ens_member%cwm     - ens_mean%cwm)**2.0
       ens_member%f_ice   = (ens_member%f_ice   - ens_mean%f_ice)**2.0
       ens_member%f_rain  = (ens_member%f_rain  - ens_mean%f_rain)**2.0
       ens_member%f_rimef = (ens_member%f_rimef - ens_mean%f_rimef)**2.0
       ens_member%smc     = (ens_member%smc     - ens_mean%smc)**2.0
       ens_member%stc     = (ens_member%stc     - ens_mean%stc)**2.0
       ens_member%pd      = (ens_member%pd      - ens_mean%pd)**2.0
       ens_member%tsk     = (ens_member%tsk     - ens_mean%tsk)**2.0
       ens_member%sst     = (ens_member%sst     - ens_mean%sst)**2.0
       ens_member%sno     = (ens_member%sno     - ens_mean%sno)**2.0
       ens_member%sice    = (ens_member%sice    - ens_mean%sice)**2.0
       ens_member%u10     = (ens_member%u10     - ens_mean%u10)**2.0
       ens_member%v10     = (ens_member%v10     - ens_mean%v10)**2.0

    end if ! if(mpi_procid .ne. mpi_masternode .and. mpi_procid
           ! .le. nanals)

    ! Compute local variables

    call mpi_reduce(ens_member%pint,ens_vari%pint,(ens_vari%west_east*     &
         & ens_vari%south_north*ens_vari%bottom_top_stag),mpi_real,        &
         & mpi_sum,mpi_masternode,mpi_comm_world,mpi_ierror)
    call mpi_reduce(ens_member%t,ens_vari%t,(ens_vari%west_east*           &
         & ens_vari%south_north*ens_vari%bottom_top),mpi_real,mpi_sum,     &
         & mpi_masternode,mpi_comm_world,mpi_ierror)
    call mpi_reduce(ens_member%q,ens_vari%q,(ens_vari%west_east*           &
         & ens_vari%south_north*ens_vari%bottom_top),mpi_real,mpi_sum,     &
         & mpi_masternode,mpi_comm_world,mpi_ierror)
    call mpi_reduce(ens_member%u,ens_vari%u,(ens_vari%west_east*           &
         & ens_vari%south_north*ens_vari%bottom_top),mpi_real,mpi_sum,     &
         & mpi_masternode,mpi_comm_world,mpi_ierror)
    call mpi_reduce(ens_member%v,ens_vari%v,(ens_vari%west_east*           &
         & ens_vari%south_north*ens_vari%bottom_top),mpi_real,mpi_sum,     &
         & mpi_masternode,mpi_comm_world,mpi_ierror)
    call mpi_reduce(ens_member%cwm,ens_vari%cwm,(ens_vari%west_east*       &
         & ens_vari%south_north*ens_vari%bottom_top),mpi_real,mpi_sum,     &
         & mpi_masternode,mpi_comm_world,mpi_ierror)
    call mpi_reduce(ens_member%f_ice,ens_vari%f_ice,(ens_vari%west_east*   &
         & ens_vari%south_north*ens_vari%bottom_top),mpi_real,mpi_sum,     &
         & mpi_masternode,mpi_comm_world,mpi_ierror)
    call mpi_reduce(ens_member%f_rain,ens_vari%f_rain,(ens_vari%west_east* &
         & ens_vari%south_north*ens_vari%bottom_top),mpi_real,mpi_sum,     &
         & mpi_masternode,mpi_comm_world,mpi_ierror)
    call mpi_reduce(ens_member%f_rimef,ens_vari%f_rimef,                   &
         & (ens_vari%west_east*ens_vari%south_north*ens_vari%bottom_top),  &
         & mpi_real,mpi_sum,mpi_masternode,mpi_comm_world,mpi_ierror)
    call mpi_reduce(ens_member%smc,ens_vari%smc,(ens_vari%west_east*       &
         & ens_vari%south_north*ens_vari%soil_layers_stag),mpi_real,       &
         & mpi_sum,mpi_masternode,mpi_comm_world,mpi_ierror)
    call mpi_reduce(ens_member%stc,ens_vari%stc,(ens_vari%west_east*       &
         & ens_vari%south_north*ens_vari%soil_layers_stag),mpi_real,       &
         & mpi_sum,mpi_masternode,mpi_comm_world,mpi_ierror)
    call mpi_reduce(ens_member%pd,ens_vari%pd,(ens_vari%west_east*         &
         & ens_vari%south_north),mpi_real,mpi_sum,mpi_masternode,          &
         & mpi_comm_world,mpi_ierror)
    call mpi_reduce(ens_member%tsk,ens_vari%tsk,(ens_vari%west_east*       &
         & ens_vari%south_north),mpi_real,mpi_sum,mpi_masternode,          &
         & mpi_comm_world,mpi_ierror)
    call mpi_reduce(ens_member%sst,ens_vari%sst,(ens_vari%west_east*       &
         & ens_vari%south_north),mpi_real,mpi_sum,mpi_masternode,          &
         & mpi_comm_world,mpi_ierror)
    call mpi_reduce(ens_member%sno,ens_vari%sno,(ens_vari%west_east*       &
         & ens_vari%south_north),mpi_real,mpi_sum,mpi_masternode,          &
         & mpi_comm_world,mpi_ierror)
    call mpi_reduce(ens_member%sice,ens_vari%sice,(ens_vari%west_east*     &
         & ens_vari%south_north),mpi_real,mpi_sum,mpi_masternode,          &
         & mpi_comm_world,mpi_ierror)
    call mpi_reduce(ens_member%u10,ens_vari%u10,(ens_vari%west_east*       &
         & ens_vari%south_north),mpi_real,mpi_sum,mpi_masternode,          &
         & mpi_comm_world,mpi_ierror)
    call mpi_reduce(ens_member%v10,ens_vari%v10,(ens_vari%west_east*       &
         & ens_vari%south_north),mpi_real,mpi_sum,mpi_masternode,          &
         & mpi_comm_world,mpi_ierror)
    ens_vari%pint    = sqrt(ens_vari%pint   /real(nanals))
    ens_vari%t       = sqrt(ens_vari%t      /real(nanals))
    ens_vari%q       = sqrt(ens_vari%q      /real(nanals))
    ens_vari%u       = sqrt(ens_vari%u      /real(nanals))
    ens_vari%v       = sqrt(ens_vari%v      /real(nanals))
    ens_vari%cwm     = sqrt(ens_vari%cwm    /real(nanals))
    ens_vari%f_ice   = sqrt(ens_vari%f_ice  /real(nanals))
    ens_vari%f_rain  = sqrt(ens_vari%f_rain /real(nanals))
    ens_vari%f_rimef = sqrt(ens_vari%f_rimef/real(nanals))
    ens_vari%smc     = sqrt(ens_vari%smc    /real(nanals))
    ens_vari%stc     = sqrt(ens_vari%stc    /real(nanals))
    ens_vari%pd      = sqrt(ens_vari%pd     /real(nanals))
    ens_vari%tsk     = sqrt(ens_vari%tsk    /real(nanals))
    ens_vari%sst     = sqrt(ens_vari%sst    /real(nanals))
    ens_vari%sno     = sqrt(ens_vari%sno    /real(nanals))
    ens_vari%sice    = sqrt(ens_vari%sice   /real(nanals))
    ens_vari%u10     = sqrt(ens_vari%u10    /real(nanals))
    ens_vari%v10     = sqrt(ens_vari%v10    /real(nanals))

    !=====================================================================

  end subroutine common_wrfnmm_ensemble_variance_compute

  !=======================================================================

  ! common_wrfnmm_ensemble_recenter_initialize.f90:

  !-----------------------------------------------------------------------

  subroutine common_wrfnmm_ensemble_recenter_initialize(grid,grid_ens,     &
       & grid_interp)

    ! Define variables passed to routine

    type(wrfnmm_prognostics)                                             :: grid
    type(wrfnmm_variables)                                               :: grid_ens
    type(wrfnmm_prognostics)                                             :: grid_interp

    !=====================================================================

    ! Initialize local variables

    call init_constants_derived()
    call init_constants(.true.)

    ! Define local variables

    ncstatus = nf90_open(path=trim(ensemble_recenter_filename),            &
         & mode=nf90_nowrite,ncid=ncfileid)
    ncstatus = nf90_inq_dimid(ncfileid,'west_east',ncdimid)
    ncstatus = nf90_inquire_dimension(ncfileid,ncdimid,                    &
         & len=grid%west_east)
    ncstatus = nf90_inq_dimid(ncfileid,'south_north',ncdimid)
    ncstatus = nf90_inquire_dimension(ncfileid,ncdimid,                    &
         & len=grid%south_north)
    ncstatus = nf90_inq_dimid(ncfileid,'bottom_top',ncdimid)
    ncstatus = nf90_inquire_dimension(ncfileid,ncdimid,                    &
         & len=grid%bottom_top)
    ncstatus = nf90_close(ncfileid)

    ! Allocate memory for local variables

    if(.not. allocated(grid%t))                                            &
         & allocate(grid%t(grid%west_east,grid%south_north,                &
         & grid%bottom_top))
    if(.not. allocated(grid%q))                                            &
         & allocate(grid%q(grid%west_east,grid%south_north,                &
         & grid%bottom_top))
    if(.not. allocated(grid%u))                                            &
         & allocate(grid%u(grid%west_east,grid%south_north,                &
         & grid%bottom_top))
    if(.not. allocated(grid%v))                                            &
         & allocate(grid%v(grid%west_east,grid%south_north,                &
         & grid%bottom_top))
    if(.not. allocated(grid%cwm))                                          &
         & allocate(grid%cwm(grid%west_east,grid%south_north,              &
         & grid%bottom_top))
    if(.not. allocated(grid%pd))                                           &
         & allocate(grid%pd(grid%west_east,grid%south_north))
    if(.not. allocated(grid%xlong_m))                                      &
         & allocate(grid%xlong_m(grid%west_east,grid%south_north))
    if(.not. allocated(grid%xlat_m))                                       &
         & allocate(grid%xlat_m(grid%west_east,grid%south_north))
    if(.not. allocated(grid%xlong_v))                                      &
         & allocate(grid%xlong_v(grid%west_east,grid%south_north))
    if(.not. allocated(grid%xlat_v))                                       &
         & allocate(grid%xlat_v(grid%west_east,grid%south_north))

    ! Define local variables
    
    ncstatus = nf90_open(path=trim(ensemble_recenter_filename),            &
         & mode=nf90_nowrite,ncid=ncfileid)
    ncstatus = nf90_inq_varid(ncfileid,'T',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,grid%t)
    ncstatus = nf90_inq_varid(ncfileid,'Q',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,grid%q)
    ncstatus = nf90_inq_varid(ncfileid,'U',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,grid%u)
    ncstatus = nf90_inq_varid(ncfileid,'V',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,grid%v)
    ncstatus = nf90_inq_varid(ncfileid,'CWM',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,grid%cwm)
    ncstatus = nf90_inq_varid(ncfileid,'PD',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,grid%pd)
    ncstatus = nf90_close(ncfileid)
    ncstatus =                                                             &
         & nf90_open(path=trim(ensemble_recenter_geogrid_filename),        &
         & mode=nf90_nowrite,ncid=ncfileid)
    ncstatus = nf90_inq_varid(ncfileid,'XLONG_M',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,grid%xlong_m)
    ncstatus = nf90_inq_varid(ncfileid,'XLAT_M',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,grid%xlat_m)
    ncstatus = nf90_inq_varid(ncfileid,'XLONG_V',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,grid%xlong_v)
    ncstatus = nf90_inq_varid(ncfileid,'XLAT_V',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,grid%xlat_v)
    ncstatus = nf90_close(ncfileid)
    ncstatus = nf90_open(path=trim(ensemble_geogrid_filename),mode=        &
         & nf90_nowrite,ncid=ncfileid)
    ncstatus = nf90_inq_varid(ncfileid,'XLONG_M',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,grid_ens%xlong_m)
    ncstatus = nf90_inq_varid(ncfileid,'XLAT_M',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,grid_ens%xlat_m)
    ncstatus = nf90_inq_varid(ncfileid,'XLONG_V',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,grid_ens%xlong_v)
    ncstatus = nf90_inq_varid(ncfileid,'XLAT_V',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,grid_ens%xlat_v)
    ncstatus = nf90_close(ncfileid)

    ! Define local variables

    grid_interp%west_east   = grid_ens%west_east
    grid_interp%south_north = grid_ens%south_north
    grid_interp%bottom_top  = grid_ens%bottom_top
    
    ! Allocate memory for local variables

    if(.not. allocated(grid_interp%t))                                     &
         & allocate(grid_interp%t(grid_interp%west_east,                   &
         & grid_interp%south_north,grid_interp%bottom_top))
    if(.not. allocated(grid_interp%q))                                     &
         & allocate(grid_interp%q(grid_interp%west_east,                   &
         & grid_interp%south_north,grid_interp%bottom_top))
    if(.not. allocated(grid_interp%u))                                     &
         & allocate(grid_interp%u(grid_interp%west_east,                   &
         & grid_interp%south_north,grid_interp%bottom_top))
    if(.not. allocated(grid_interp%v))                                     &
         & allocate(grid_interp%v(grid_interp%west_east,                   &
         & grid_interp%south_north,grid_interp%bottom_top))
    if(.not. allocated(grid_interp%cwm))                                   &
         & allocate(grid_interp%cwm(grid_interp%west_east,                 &
         & grid_interp%south_north,grid_interp%bottom_top))
    if(.not. allocated(grid_interp%pd))                                    &
         & allocate(grid_interp%pd(grid_interp%west_east,                  &
         & grid_interp%south_north))
    if(.not. allocated(grid_interp%xlong_m))                               &
         & allocate(grid_interp%xlong_m(grid_interp%west_east,             &
         & grid_interp%south_north))
    if(.not. allocated(grid_interp%xlat_m))                                &
         & allocate(grid_interp%xlat_m(grid_interp%west_east,              &
         & grid_interp%south_north))
    if(.not. allocated(grid_interp%xlong_v))                               &
         & allocate(grid_interp%xlong_v(grid_interp%west_east,             &
         & grid_interp%south_north))
    if(.not. allocated(grid_interp%xlat_v))                                &
         & allocate(grid_interp%xlat_v(grid_interp%west_east,              &
         & grid_interp%south_north))

    ! Define local variables

    grid_interp%xlong_m = grid_ens%xlong_m
    grid_interp%xlat_m  = grid_ens%xlat_m
    grid_interp%xlong_v = grid_ens%xlong_v
    grid_interp%xlat_v  = grid_ens%xlat_v

    !=====================================================================

  end subroutine common_wrfnmm_ensemble_recenter_initialize

  !=======================================================================

  ! common_wrfnmm_ensemble_recenter_initialize.f90:

  !-----------------------------------------------------------------------

  subroutine common_wrfnmm_ensemble_recenter_interp_initialize(grid)

    ! Define variables passed to routine

    type(wrfnmm_prognostics)                                             :: grid

    !=====================================================================

    ! Initialize local variables

    call init_constants_derived()
    call init_constants(.true.)

    ! Define local variables

    ncstatus = nf90_open(path=trim(ensemble_recenter_interp_filename),     &
         & mode=nf90_nowrite,ncid=ncfileid)
    ncstatus = nf90_inq_dimid(ncfileid,'west_east',ncdimid)
    ncstatus = nf90_inquire_dimension(ncfileid,ncdimid,                    &
         & len=grid%west_east)
    ncstatus = nf90_inq_dimid(ncfileid,'south_north',ncdimid)
    ncstatus = nf90_inquire_dimension(ncfileid,ncdimid,                    &
         & len=grid%south_north)
    ncstatus = nf90_inq_dimid(ncfileid,'bottom_top',ncdimid)
    ncstatus = nf90_inquire_dimension(ncfileid,ncdimid,                    &
         & len=grid%bottom_top)

    ! Allocate memory for local variables

    if(.not. allocated(grid%t))                                            &
         & allocate(grid%t(grid%west_east,grid%south_north,                &
         & grid%bottom_top))
    if(.not. allocated(grid%q))                                            &
         & allocate(grid%q(grid%west_east,grid%south_north,                &
         & grid%bottom_top))
    if(.not. allocated(grid%u))                                            &
         & allocate(grid%u(grid%west_east,grid%south_north,                &
         & grid%bottom_top))
    if(.not. allocated(grid%v))                                            &
         & allocate(grid%v(grid%west_east,grid%south_north,                &
         & grid%bottom_top))
    if(.not. allocated(grid%cwm))                                          &
         & allocate(grid%cwm(grid%west_east,grid%south_north,              &
         & grid%bottom_top))
    if(.not. allocated(grid%pd))                                           &
         & allocate(grid%pd(grid%west_east,grid%south_north))

    ! Define local variables

    ncstatus = nf90_inq_varid(ncfileid,'T',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,grid%t)
    ncstatus = nf90_inq_varid(ncfileid,'Q',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,grid%q)
    ncstatus = nf90_inq_varid(ncfileid,'U',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,grid%u)
    ncstatus = nf90_inq_varid(ncfileid,'V',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,grid%v)
    ncstatus = nf90_inq_varid(ncfileid,'CWM',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,grid%cwm)
    ncstatus = nf90_inq_varid(ncfileid,'PD',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,grid%pd)
    ncstatus = nf90_close(ncfileid)

    !=====================================================================

  end subroutine common_wrfnmm_ensemble_recenter_interp_initialize

  !=======================================================================

  ! common_wrfnmm_ensemble_variance_initialize.f90:

  !-----------------------------------------------------------------------

  subroutine common_wrfnmm_ensemble_variance_initialize(grid)

    ! Define variables passed to routine

    type(wrfnmm_variables)                                               :: grid
    
    !=====================================================================

    ! Define local variables

    ncstatus = nf90_open(path=trim(ensemble_variance_filename),            &
         & mode=nf90_nowrite,ncid=ncfileid)
    ncstatus = nf90_inq_dimid(ncfileid,'west_east',ncdimid)
    ncstatus = nf90_inquire_dimension(ncfileid,ncdimid,                    &
         & len=grid%west_east)
    ncstatus = nf90_inq_dimid(ncfileid,'south_north',ncdimid)
    ncstatus = nf90_inquire_dimension(ncfileid,ncdimid,                    &
         & len=grid%south_north)
    ncstatus = nf90_inq_dimid(ncfileid,'bottom_top',ncdimid)
    ncstatus = nf90_inquire_dimension(ncfileid,ncdimid,                    &
         & len=grid%bottom_top)
    ncstatus = nf90_inq_dimid(ncfileid,'bottom_top_stag',ncdimid)
    ncstatus = nf90_inquire_dimension(ncfileid,ncdimid,                    &
         & len=grid%bottom_top_stag)
    ncstatus = nf90_inq_dimid(ncfileid,'soil_layers_stag',ncdimid)
    ncstatus = nf90_inquire_dimension(ncfileid,ncdimid,                    &
         & len=grid%soil_layers_stag)
    ncstatus = nf90_close(ncfileid)

    ! Allocate memory for local variables

    if(.not. allocated(grid%pint))                                         &
         & allocate(grid%pint(grid%west_east,grid%south_north,             &
         & grid%bottom_top_stag))
    if(.not. allocated(grid%t))                                            &
         & allocate(grid%t(grid%west_east,grid%south_north,                &
         & grid%bottom_top))
    if(.not. allocated(grid%q))                                            &
         & allocate(grid%q(grid%west_east,grid%south_north,                &
         & grid%bottom_top))
    if(.not. allocated(grid%u))                                            &
         & allocate(grid%u(grid%west_east,grid%south_north,                &
         & grid%bottom_top))
    if(.not. allocated(grid%v))                                            &
         & allocate(grid%v(grid%west_east,grid%south_north,                &
         & grid%bottom_top))
    if(.not. allocated(grid%cwm))                                          &
         & allocate(grid%cwm(grid%west_east,grid%south_north,              &
         & grid%bottom_top))
    if(.not. allocated(grid%f_ice))                                        &
         & allocate(grid%f_ice(grid%west_east,grid%south_north,            &
         & grid%bottom_top))
    if(.not. allocated(grid%f_rain))                                       &
          & allocate(grid%f_rain(grid%west_east,grid%south_north,           &
         & grid%bottom_top))
    if(.not. allocated(grid%f_rimef))                                      &
         & allocate(grid%f_rimef(grid%west_east,grid%south_north,          &
         & grid%bottom_top))
    if(.not. allocated(grid%smc))                                          &
         & allocate(grid%smc(grid%west_east,grid%south_north,              &
         & grid%soil_layers_stag))
    if(.not. allocated(grid%stc))                                          &
         & allocate(grid%stc(grid%west_east,grid%south_north,              &
         & grid%soil_layers_stag))
    if(.not. allocated(grid%pd))                                           &
         & allocate(grid%pd(grid%west_east,grid%south_north))
    if(.not. allocated(grid%tsk))                                          &
         & allocate(grid%tsk(grid%west_east,grid%south_north))
    if(.not. allocated(grid%sst))                                          &
         & allocate(grid%sst(grid%west_east,grid%south_north))
    if(.not. allocated(grid%sno))                                          &
         & allocate(grid%sno(grid%west_east,grid%south_north))
    if(.not. allocated(grid%sice))                                         &
         & allocate(grid%sice(grid%west_east,grid%south_north))
    if(.not. allocated(grid%u10))                                          &
         & allocate(grid%u10(grid%west_east,grid%south_north))
    if(.not. allocated(grid%v10))                                          &
         & allocate(grid%v10(grid%west_east,grid%south_north))    
    if(.not. allocated(grid%xlong_m))                                      &
         & allocate(grid%xlong_m(grid%west_east,grid%south_north))   
    if(.not. allocated(grid%xlat_m))                                       &
         & allocate(grid%xlat_m(grid%west_east,grid%south_north))       
    if(.not. allocated(grid%xlong_v))                                      &
         & allocate(grid%xlong_v(grid%west_east,grid%south_north))   
    if(.not. allocated(grid%xlat_v))                                       &
         & allocate(grid%xlat_v(grid%west_east,grid%south_north)) 

    ! Initialize local variables

    grid%pint    = 0.0
    grid%t       = 0.0
    grid%q       = 0.0
    grid%u       = 0.0
    grid%v       = 0.0
    grid%cwm     = 0.0
    grid%f_ice   = 0.0
    grid%f_rain  = 0.0
    grid%f_rimef = 0.0
    grid%smc     = 0.0    
    grid%stc     = 0.0
    grid%pd      = 0.0
    grid%tsk     = 0.0
    grid%sst     = 0.0
    grid%sno     = 0.0
    grid%sice    = 0.0
    grid%u10     = 0.0
    grid%v10     = 0.0

    !=====================================================================

  end subroutine common_wrfnmm_ensemble_variance_initialize

  !=======================================================================

  ! common_wrfnmm_ensemble_mean_initialize.f90:

  !-----------------------------------------------------------------------

  subroutine common_wrfnmm_ensemble_mean_initialize(grid)

    ! Define variables passed to routine

    type(wrfnmm_variables)                                               :: grid
    
    !=====================================================================

    ! Define local variables

    ncstatus = nf90_open(path=trim(ensemble_mean_filename),                &
         & mode=nf90_nowrite,ncid=ncfileid)
    ncstatus = nf90_inq_dimid(ncfileid,'west_east',ncdimid)
    ncstatus = nf90_inquire_dimension(ncfileid,ncdimid,                    &
         & len=grid%west_east)
    ncstatus = nf90_inq_dimid(ncfileid,'south_north',ncdimid)
    ncstatus = nf90_inquire_dimension(ncfileid,ncdimid,                    &
         & len=grid%south_north)
    ncstatus = nf90_inq_dimid(ncfileid,'bottom_top',ncdimid)
    ncstatus = nf90_inquire_dimension(ncfileid,ncdimid,                    &
         & len=grid%bottom_top)
    ncstatus = nf90_inq_dimid(ncfileid,'bottom_top_stag',ncdimid)
    ncstatus = nf90_inquire_dimension(ncfileid,ncdimid,                    &
         & len=grid%bottom_top_stag)
    ncstatus = nf90_inq_dimid(ncfileid,'soil_layers_stag',ncdimid)
    ncstatus = nf90_inquire_dimension(ncfileid,ncdimid,                    &
         & len=grid%soil_layers_stag)
    ncstatus = nf90_close(ncfileid)

    ! Allocate memory for local variables

    if(.not. allocated(grid%pint))                                         &
         & allocate(grid%pint(grid%west_east,grid%south_north,             &
         & grid%bottom_top_stag))
    if(.not. allocated(grid%t))                                            &
         & allocate(grid%t(grid%west_east,grid%south_north,                &
         & grid%bottom_top))
    if(.not. allocated(grid%q))                                            &
         & allocate(grid%q(grid%west_east,grid%south_north,                &
         & grid%bottom_top))
    if(.not. allocated(grid%u))                                            &
         & allocate(grid%u(grid%west_east,grid%south_north,                &
         & grid%bottom_top))
    if(.not. allocated(grid%v))                                            &
         & allocate(grid%v(grid%west_east,grid%south_north,                &
         & grid%bottom_top))
    if(.not. allocated(grid%cwm))                                          &
         & allocate(grid%cwm(grid%west_east,grid%south_north,              &
         & grid%bottom_top))
    if(.not. allocated(grid%f_ice))                                        &
         & allocate(grid%f_ice(grid%west_east,grid%south_north,            &
         & grid%bottom_top))
    if(.not. allocated(grid%f_rain))                                       &
          & allocate(grid%f_rain(grid%west_east,grid%south_north,           &
         & grid%bottom_top))
    if(.not. allocated(grid%f_rimef))                                      &
         & allocate(grid%f_rimef(grid%west_east,grid%south_north,          &
         & grid%bottom_top))
    if(.not. allocated(grid%smc))                                          &
         & allocate(grid%smc(grid%west_east,grid%south_north,              &
         & grid%soil_layers_stag))
    if(.not. allocated(grid%stc))                                          &
         & allocate(grid%stc(grid%west_east,grid%south_north,              &
         & grid%soil_layers_stag))
    if(.not. allocated(grid%pd))                                           &
         & allocate(grid%pd(grid%west_east,grid%south_north))
    if(.not. allocated(grid%tsk))                                          &
         & allocate(grid%tsk(grid%west_east,grid%south_north))
    if(.not. allocated(grid%sst))                                          &
         & allocate(grid%sst(grid%west_east,grid%south_north))
    if(.not. allocated(grid%sno))                                          &
         & allocate(grid%sno(grid%west_east,grid%south_north))
    if(.not. allocated(grid%sice))                                         &
         & allocate(grid%sice(grid%west_east,grid%south_north))
    if(.not. allocated(grid%u10))                                          &
         & allocate(grid%u10(grid%west_east,grid%south_north))
    if(.not. allocated(grid%v10))                                          &
         & allocate(grid%v10(grid%west_east,grid%south_north))    
    if(.not. allocated(grid%xlong_m))                                      &
         & allocate(grid%xlong_m(grid%west_east,grid%south_north))   
    if(.not. allocated(grid%xlat_m))                                       &
         & allocate(grid%xlat_m(grid%west_east,grid%south_north))       
    if(.not. allocated(grid%xlong_v))                                      &
         & allocate(grid%xlong_v(grid%west_east,grid%south_north))   
    if(.not. allocated(grid%xlat_v))                                       &
         & allocate(grid%xlat_v(grid%west_east,grid%south_north)) 

    ! Initialize local variables

    grid%pint    = 0.0
    grid%t       = 0.0
    grid%q       = 0.0
    grid%u       = 0.0
    grid%v       = 0.0
    grid%cwm     = 0.0
    grid%f_ice   = 0.0
    grid%f_rain  = 0.0
    grid%f_rimef = 0.0
    grid%smc     = 0.0    
    grid%stc     = 0.0
    grid%pd      = 0.0
    grid%tsk     = 0.0
    grid%sst     = 0.0
    grid%sno     = 0.0
    grid%sice    = 0.0
    grid%u10     = 0.0
    grid%v10     = 0.0

    !=====================================================================

  end subroutine common_wrfnmm_ensemble_mean_initialize

  !=======================================================================

  ! common_wrfnmm_ensemble_update_initialize.f90:

  !-----------------------------------------------------------------------

  subroutine common_wrfnmm_ensemble_update_initialize(grid)

    ! Define variables passed to routine

    type(wrfnmm_prognostics)                                             :: grid
    
    !=====================================================================

    ! Define local variables

    ncstatus = nf90_open(path=trim(ens_infilename),mode=nf90_nowrite,      &
         & ncid=ncfileid)
    ncstatus = nf90_inq_dimid(ncfileid,'west_east',ncdimid)
    ncstatus = nf90_inquire_dimension(ncfileid,ncdimid,                    &
         & len=grid%west_east)
    ncstatus = nf90_inq_dimid(ncfileid,'south_north',ncdimid)
    ncstatus = nf90_inquire_dimension(ncfileid,ncdimid,                    &
         & len=grid%south_north)
    ncstatus = nf90_inq_dimid(ncfileid,'bottom_top',ncdimid)
    ncstatus = nf90_inquire_dimension(ncfileid,ncdimid,                    &
         & len=grid%bottom_top)
    ncstatus = nf90_inq_dimid(ncfileid,'bottom_top_stag',ncdimid)
    ncstatus = nf90_inquire_dimension(ncfileid,ncdimid,                    &
         & len=grid%bottom_top_stag)
    ncstatus = nf90_close(ncfileid)

    ! Allocate memory for local variables

    if(.not. allocated(grid%pint))                                         &
         & allocate(grid%pint(grid%west_east,grid%south_north,             &
         & grid%bottom_top_stag))
    if(.not. allocated(grid%t))                                            &
         & allocate(grid%t(grid%west_east,grid%south_north,                &
         & grid%bottom_top))
    if(.not. allocated(grid%q))                                            &
         & allocate(grid%q(grid%west_east,grid%south_north,                &
         & grid%bottom_top))
    if(.not. allocated(grid%u))                                            &
         & allocate(grid%u(grid%west_east,grid%south_north,                &
         & grid%bottom_top))
    if(.not. allocated(grid%v))                                            &
         & allocate(grid%v(grid%west_east,grid%south_north,                &
         & grid%bottom_top))
    if(.not. allocated(grid%cwm))                                          &
         & allocate(grid%cwm(grid%west_east,grid%south_north,              &
         & grid%bottom_top)) 
    if(.not. allocated(grid%pd))                                           &
         & allocate(grid%pd(grid%west_east,grid%south_north))

    ! Define local variables

    ncstatus = nf90_open(path=trim(ens_infilename),mode=nf90_nowrite,      &
         & ncid=ncfileid)
    ncstatus = nf90_inq_varid(ncfileid,'PINT',ncvarid)
    if(ncstatus .eq. 0) grid%is_exist_pint = .true.
    ncstatus = nf90_get_var(ncfileid,ncvarid,grid%pint)
    ncstatus = nf90_inq_varid(ncfileid,'T',ncvarid)
    if(ncstatus .eq. 0) grid%is_exist_t    = .true.
    ncstatus = nf90_get_var(ncfileid,ncvarid,grid%t)
    ncstatus = nf90_inq_varid(ncfileid,'Q',ncvarid)
    if(ncstatus .eq. 0) grid%is_exist_q    = .true.
    ncstatus = nf90_get_var(ncfileid,ncvarid,grid%q)
    ncstatus = nf90_inq_varid(ncfileid,'U',ncvarid)
    if(ncstatus .eq. 0) grid%is_exist_u    = .true.
    ncstatus = nf90_get_var(ncfileid,ncvarid,grid%u)
    ncstatus = nf90_inq_varid(ncfileid,'V',ncvarid)
    if(ncstatus .eq. 0) grid%is_exist_v    = .true.
    ncstatus = nf90_get_var(ncfileid,ncvarid,grid%v)
    ncstatus = nf90_inq_varid(ncfileid,'CWM',ncvarid)
    if(ncstatus .eq. 0) grid%is_exist_cwm  = .true.
    ncstatus = nf90_get_var(ncfileid,ncvarid,grid%cwm)
    ncstatus = nf90_inq_varid(ncfileid,'PD',ncvarid)
    if(ncstatus .eq. 0) grid%is_exist_pd   = .true.
    ncstatus = nf90_get_var(ncfileid,ncvarid,grid%pd)
    ncstatus = nf90_close(ncfileid)

    !=====================================================================

  end subroutine common_wrfnmm_ensemble_update_initialize

  !=======================================================================

  ! common_wrfnmm_ensemble_member_initialize.f90:

  !-----------------------------------------------------------------------

  subroutine common_wrfnmm_ensemble_member_initialize(grid)

    ! Define variables passed to routine

    type(wrfnmm_variables)                                               :: grid
    
    !=====================================================================

    ! Allocate memory for local variables

    if(.not. allocated(grid%pint))                                         &
         & allocate(grid%pint(grid%west_east,grid%south_north,             &
         & grid%bottom_top_stag))
    if(.not. allocated(grid%t))                                            &
         & allocate(grid%t(grid%west_east,grid%south_north,                &
         & grid%bottom_top))
    if(.not. allocated(grid%q))                                            &
         & allocate(grid%q(grid%west_east,grid%south_north,                &
         & grid%bottom_top))
    if(.not. allocated(grid%u))                                            &
         & allocate(grid%u(grid%west_east,grid%south_north,                &
         & grid%bottom_top))
    if(.not. allocated(grid%v))                                            &
         & allocate(grid%v(grid%west_east,grid%south_north,                &
         & grid%bottom_top))
    if(.not. allocated(grid%cwm))                                          &
         & allocate(grid%cwm(grid%west_east,grid%south_north,              &
         & grid%bottom_top))
    if(.not. allocated(grid%f_ice))                                        &
         & allocate(grid%f_ice(grid%west_east,grid%south_north,            &
         & grid%bottom_top))
    if(.not. allocated(grid%f_rain))                                       &
         & allocate(grid%f_rain(grid%west_east,grid%south_north,           &
         & grid%bottom_top))
    if(.not. allocated(grid%f_rimef))                                      &
         & allocate(grid%f_rimef(grid%west_east,grid%south_north,          &
         & grid%bottom_top))
    if(.not. allocated(grid%smc))                                          &
         & allocate(grid%smc(grid%west_east,grid%south_north,              &
         & grid%soil_layers_stag))
    if(.not. allocated(grid%stc))                                          &
         & allocate(grid%stc(grid%west_east,grid%south_north,              &
         & grid%soil_layers_stag))
    if(.not. allocated(grid%pd))                                           &
         & allocate(grid%pd(grid%west_east,grid%south_north))
    if(.not. allocated(grid%tsk))                                          &
         & allocate(grid%tsk(grid%west_east,grid%south_north))
    if(.not. allocated(grid%sst))                                          &
         & allocate(grid%sst(grid%west_east,grid%south_north))
    if(.not. allocated(grid%sno))                                          &
         & allocate(grid%sno(grid%west_east,grid%south_north))
    if(.not. allocated(grid%sice))                                         &
         & allocate(grid%sice(grid%west_east,grid%south_north))
    if(.not. allocated(grid%u10))                                          &
         & allocate(grid%u10(grid%west_east,grid%south_north))
    if(.not. allocated(grid%v10))                                          &
         & allocate(grid%v10(grid%west_east,grid%south_north))    

    ! Initialize local variables

    grid%pint    = 0.0
    grid%t       = 0.0
    grid%q       = 0.0
    grid%u       = 0.0
    grid%v       = 0.0
    grid%cwm     = 0.0
    grid%f_ice   = 0.0
    grid%f_rain  = 0.0
    grid%f_rimef = 0.0
    grid%smc     = 0.0    
    grid%stc     = 0.0
    grid%pd      = 0.0
    grid%tsk     = 0.0
    grid%sst     = 0.0
    grid%sno     = 0.0
    grid%sice    = 0.0
    grid%u10     = 0.0
    grid%v10     = 0.0

    ! Check local variable and proceed accordingly

    if(mpi_nprocs - 1 .lt. nanals) then
       if(mpi_procid .eq. mpi_masternode) write (*,*)  &
          'ERROR: mpi_nprocs-1 is smaller than ensemble size.', &
          'mpi_nprocs has to be at least 1+ensemble size'
       call mpi_abort(mpi_comm_world,mpi_ierror)
    end if

    if(mpi_procid .ne. mpi_masternode .and. mpi_procid .le. nanals) then

       ! Define local variables

       ncstatus = nf90_open(path=trim(analysis_file(mpi_procid)),          &
            & mode=nf90_nowrite,ncid=ncfileid)
       ncstatus = nf90_inq_varid(ncfileid,'PINT',ncvarid)
       ncstatus = nf90_get_var(ncfileid,ncvarid,grid%pint)
       ncstatus = nf90_inq_varid(ncfileid,'T',ncvarid)
       ncstatus = nf90_get_var(ncfileid,ncvarid,grid%t)
       ncstatus = nf90_inq_varid(ncfileid,'Q',ncvarid)
       ncstatus = nf90_get_var(ncfileid,ncvarid,grid%q)
       ncstatus = nf90_inq_varid(ncfileid,'U',ncvarid)
       ncstatus = nf90_get_var(ncfileid,ncvarid,grid%u)
       ncstatus = nf90_inq_varid(ncfileid,'V',ncvarid)
       ncstatus = nf90_get_var(ncfileid,ncvarid,grid%v)
       ncstatus = nf90_inq_varid(ncfileid,'CWM',ncvarid)
       ncstatus = nf90_get_var(ncfileid,ncvarid,grid%cwm)
       ncstatus = nf90_inq_varid(ncfileid,'F_ICE',ncvarid)
       ncstatus = nf90_get_var(ncfileid,ncvarid,grid%f_ice)
       ncstatus = nf90_inq_varid(ncfileid,'F_RAIN',ncvarid)
       ncstatus = nf90_get_var(ncfileid,ncvarid,grid%f_rain)
       ncstatus = nf90_inq_varid(ncfileid,'F_RIMEF',ncvarid)
       ncstatus = nf90_get_var(ncfileid,ncvarid,grid%f_rimef)
       ncstatus = nf90_inq_varid(ncfileid,'SMC',ncvarid)
       ncstatus = nf90_get_var(ncfileid,ncvarid,grid%smc)
       ncstatus = nf90_inq_varid(ncfileid,'STC',ncvarid)
       ncstatus = nf90_get_var(ncfileid,ncvarid,grid%stc)
       ncstatus = nf90_inq_varid(ncfileid,'PD',ncvarid)
       ncstatus = nf90_get_var(ncfileid,ncvarid,grid%pd)
       ncstatus = nf90_inq_varid(ncfileid,'TSK',ncvarid)
       ncstatus = nf90_get_var(ncfileid,ncvarid,grid%tsk)
       ncstatus = nf90_inq_varid(ncfileid,'SST',ncvarid)
       ncstatus = nf90_get_var(ncfileid,ncvarid,grid%sst)
       ncstatus = nf90_inq_varid(ncfileid,'SNO',ncvarid)
       ncstatus = nf90_get_var(ncfileid,ncvarid,grid%sno)
       ncstatus = nf90_inq_varid(ncfileid,'SICE',ncvarid)
       ncstatus = nf90_get_var(ncfileid,ncvarid,grid%sice)
       ncstatus = nf90_inq_varid(ncfileid,'U10',ncvarid)
       ncstatus = nf90_get_var(ncfileid,ncvarid,grid%u10)
       ncstatus = nf90_inq_varid(ncfileid,'V10',ncvarid)
       ncstatus = nf90_get_var(ncfileid,ncvarid,grid%v10)
       ncstatus = nf90_close(ncfileid)

    end if ! if(mpi_procid .ne. mpi_masternode)

    ! Enable the root task to catch up from I/O and calculations

    call mpi_barrier(mpi_comm_world,mpi_ierror)

    !=====================================================================

  end subroutine common_wrfnmm_ensemble_member_initialize

  !=======================================================================

  ! common_analysis_file_initialize.f90:

  !-----------------------------------------------------------------------

  subroutine common_analysis_file_initialize()

    ! Define counting variables

    integer                                                              :: i, j, k

    !=====================================================================

    ! Deallocate memory for local variables

    call common_analysis_file_cleanup()

    ! Allocate memory for local variable

    if(.not. allocated(analysis_file)) allocate(analysis_file(nanals))

    ! Define local variables

    open(99,file=ensemble_member_filename_list,form='formatted')

    ! Loop through local variable

    do k = 1, nanals

       ! Define local variables

       read(99,*) analysis_file(k)

    end do ! do k = 1, nanals

    ! Define local variables

    close(99)

    !=====================================================================

  end subroutine common_analysis_file_initialize

  !=======================================================================
  
  ! common_wrfnmm_ensemble_recenter_cleanup.f90:

  !-----------------------------------------------------------------------

  subroutine common_wrfnmm_ensemble_recenter_cleanup(grid,grid_interp)

    ! Define variables passed to routine

    type(wrfnmm_prognostics)                                             :: grid
    type(wrfnmm_prognostics)                                             :: grid_interp

    !=====================================================================

    ! Deallocate memory for local variables

    if(allocated(grid%pint))           deallocate(grid%pint)
    if(allocated(grid%t))              deallocate(grid%t)
    if(allocated(grid%q))              deallocate(grid%q)
    if(allocated(grid%u))              deallocate(grid%u)
    if(allocated(grid%v))              deallocate(grid%v)
    if(allocated(grid%cwm))            deallocate(grid%cwm)
    if(allocated(grid%xlong_m))        deallocate(grid%xlong_m)
    if(allocated(grid%xlat_m))         deallocate(grid%xlat_m)
    if(allocated(grid%xlong_v))        deallocate(grid%xlong_v)
    if(allocated(grid%xlat_v))         deallocate(grid%xlat_v)
    if(allocated(grid_interp%pint))    deallocate(grid_interp%pint)
    if(allocated(grid_interp%t))       deallocate(grid_interp%t)
    if(allocated(grid_interp%q))       deallocate(grid_interp%q)
    if(allocated(grid_interp%u))       deallocate(grid_interp%u)
    if(allocated(grid_interp%v))       deallocate(grid_interp%v)
    if(allocated(grid_interp%cwm))     deallocate(grid_interp%cwm)
    if(allocated(grid_interp%xlong_m)) deallocate(grid_interp%xlong_m)
    if(allocated(grid_interp%xlat_m))  deallocate(grid_interp%xlat_m)
    if(allocated(grid_interp%xlong_v)) deallocate(grid_interp%xlong_v)
    if(allocated(grid_interp%xlat_v))  deallocate(grid_interp%xlat_v)

    !=====================================================================

  end subroutine common_wrfnmm_ensemble_recenter_cleanup

  !=======================================================================
  
  ! common_wrfnmm_ensemble_member_cleanup.f90:

  !-----------------------------------------------------------------------

  subroutine common_wrfnmm_ensemble_member_cleanup(grid)

    ! Define variables passed to routine

    type(wrfnmm_variables)                                               :: grid
    
    !=====================================================================

    ! Deallocate memory for local variables

    if(allocated(grid%pint))    deallocate(grid%pint)
    if(allocated(grid%t))       deallocate(grid%t)
    if(allocated(grid%q))       deallocate(grid%q)
    if(allocated(grid%u))       deallocate(grid%u)
    if(allocated(grid%v))       deallocate(grid%v)
    if(allocated(grid%cwm))     deallocate(grid%cwm)
    if(allocated(grid%f_ice))   deallocate(grid%f_ice)
    if(allocated(grid%f_rain))  deallocate(grid%f_rain)
    if(allocated(grid%f_rimef)) deallocate(grid%f_rimef)
    if(allocated(grid%smc))     deallocate(grid%smc)
    if(allocated(grid%stc))     deallocate(grid%stc)
    if(allocated(grid%pd))      deallocate(grid%pd)
    if(allocated(grid%tsk))     deallocate(grid%tsk)
    if(allocated(grid%sst))     deallocate(grid%sst)
    if(allocated(grid%sno))     deallocate(grid%sno)
    if(allocated(grid%sice))    deallocate(grid%sice)
    if(allocated(grid%u10))     deallocate(grid%u10)
    if(allocated(grid%v10))     deallocate(grid%v10)
    if(allocated(grid%xlong_m)) deallocate(grid%xlong_m)
    if(allocated(grid%xlat_m))  deallocate(grid%xlat_m)
    if(allocated(grid%xlong_v)) deallocate(grid%xlong_v)
    if(allocated(grid%xlat_v))  deallocate(grid%xlat_v)

    !=====================================================================

  end subroutine common_wrfnmm_ensemble_member_cleanup

  !=======================================================================

  ! common_analysis_file_cleanup.f90:

  !-----------------------------------------------------------------------

  subroutine common_analysis_file_cleanup()

    !=====================================================================

    ! Deallocate memory for local variables

    if(allocated(analysis_file)) deallocate(analysis_file)
    
    !=====================================================================

  end subroutine common_analysis_file_cleanup

  !=======================================================================

end module common_interface
