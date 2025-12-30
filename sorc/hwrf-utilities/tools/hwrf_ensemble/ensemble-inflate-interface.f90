module ensemble_inflate_interface

  !=======================================================================

  ! Define associated modules and subroutines

  !-----------------------------------------------------------------------

  use constants
  use kinds

  !-----------------------------------------------------------------------

  use common_interface
  use mpi_interface
  use namelist
  use netcdf

  !-----------------------------------------------------------------------

  implicit none

  !-----------------------------------------------------------------------

contains

  !=======================================================================

  ! ensemble_inflate.f90:

  !-----------------------------------------------------------------------

  subroutine ensemble_inflate()

    !=====================================================================

    ! Define local variables

    call common_analysis_file_initialize()

    ! Compute local variables

    call ensemble_inflate_wrfnmm()

    ! Deallocate memory for local variables

    call common_analysis_file_cleanup()

    !=====================================================================

  end subroutine ensemble_inflate

  !=======================================================================

  ! ensemble_inflate_wrfnmm.f90:

  !-----------------------------------------------------------------------

  subroutine ensemble_inflate_wrfnmm()

    ! Define variables computed within routine

    type(wrfnmm_variables)                                               :: ens_member

    !=====================================================================

    ! Compute local variables

    call ensemble_inflate_wrfnmm_compute(ens_member)

    ! Write variables to external file

    call ensemble_inflate_wrfnmm_write(ens_member)

    ! Deallocate memory for local variables

    call common_wrfnmm_ensemble_member_cleanup(ens_member)

    !=====================================================================

  end subroutine ensemble_inflate_wrfnmm

  !=======================================================================

  ! ensemble_inflate_wrfnmm_write.f90:

  !-----------------------------------------------------------------------

  subroutine ensemble_inflate_wrfnmm_write(grid)

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

  end subroutine ensemble_inflate_wrfnmm_write

  !=======================================================================
  
  ! ensemble_inflate_wrfnmm_compute.f90:

  !-----------------------------------------------------------------------

  subroutine ensemble_inflate_wrfnmm_compute(ens_member)

    ! Define variables passed to routine

    type(wrfnmm_variables)                                               :: ens_member

    ! Define variables computed within routine

    type(wrfnmm_variables)                                               :: ens_mean

    !=====================================================================

    ! Initialize local variables

    call init_constants_derived()
    call init_constants(.true.)

    ! Define local variables

    ensemble_mean_filename      = trim(adjustl(analysis_file(1)))
    call common_wrfnmm_ensemble_mean_initialize(ens_mean)
    ens_member%west_east        = ens_mean%west_east
    ens_member%south_north      = ens_mean%south_north
    ens_member%bottom_top       = ens_mean%bottom_top
    ens_member%bottom_top_stag  = ens_mean%bottom_top_stag
    ens_member%soil_layers_stag = ens_mean%soil_layers_stag
    call common_wrfnmm_ensemble_member_initialize(ens_member)
    
    ! Compute local variables

    call common_wrfnmm_ensemble_mean_compute(ens_member,ens_mean)

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

    ! Deallocate memory for local variables

    call common_wrfnmm_ensemble_member_cleanup(ens_member)
    call common_analysis_file_cleanup()

    ! Define local variables

    ensemble_member_filename_list = ensemble_analysis_filename_list
    call common_analysis_file_initialize()
    ens_member%west_east          = ens_mean%west_east
    ens_member%south_north        = ens_mean%south_north
    ens_member%bottom_top         = ens_mean%bottom_top
    ens_member%bottom_top_stag    = ens_mean%bottom_top_stag
    ens_member%soil_layers_stag   = ens_mean%soil_layers_stag
    call common_wrfnmm_ensemble_member_initialize(ens_member)

    ! Check local variable and proceed accordingly

    if(mpi_procid .ne. mpi_masternode .and. mpi_procid .le. nanals) then

       ! Compute local variables

       ens_member%pd  = ens_member%pd  - inflate_coeff*(ens_mean%pd  -    &
            & ens_member%pd) 
       ens_member%u   = ens_member%u   - inflate_coeff*(ens_mean%u   -    &
            & ens_member%u) 
       ens_member%v   = ens_member%v   - inflate_coeff*(ens_mean%v   -    &
            & ens_member%v)
       ens_member%t   = ens_member%t   - inflate_coeff*(ens_mean%t   -    &
            & ens_member%t)
       ens_member%q   = ens_member%q   - inflate_coeff*(ens_mean%q   -    &
            & ens_member%q)
       ens_member%cwm = ens_member%cwm - inflate_coeff*(ens_mean%cwm -    &
            & ens_member%cwm)

    end if ! if(mpi_procid .ne. mpi_masternode .and. mpi_procid
           ! .le. nanals)

    ! Enable the root task to catch up from I/O and calculations

    call mpi_barrier(mpi_comm_world,mpi_ierror)

    ! Deallocate memory for local variables

    call common_analysis_file_cleanup()
    call common_wrfnmm_ensemble_member_cleanup(ens_mean)

    !=====================================================================

  end subroutine ensemble_inflate_wrfnmm_compute

  !=======================================================================

end module ensemble_inflate_interface
