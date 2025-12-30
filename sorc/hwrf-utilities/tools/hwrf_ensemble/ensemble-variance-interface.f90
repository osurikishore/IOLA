module ensemble_variance_interface

  !=======================================================================

  ! Define associated modules and subroutines

  !-----------------------------------------------------------------------

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

  ! ensemble_variance.f90:

  !-----------------------------------------------------------------------

  subroutine ensemble_variance()

    !=====================================================================

    ! Define local variables

    call common_analysis_file_initialize()

    ! Compute local variables

    call ensemble_variance_wrfnmm()

    ! Deallocate memory for local variables

    call common_analysis_file_cleanup()

    !=====================================================================

  end subroutine ensemble_variance

  !=======================================================================

  ! ensemble_variance_wrfnmm.f90:

  !-----------------------------------------------------------------------

  subroutine ensemble_variance_wrfnmm()

    ! Define variables computed within routine

    type(wrfnmm_variables)                                               :: ens_member
    type(wrfnmm_variables)                                               :: ens_mean
    type(wrfnmm_variables)                                               :: ens_vari

    ! Define counting variables

    integer                                                              :: i, j, k

    !=====================================================================

    ! Define local variables

    ensemble_mean_filename      = ensemble_variance_filename
    call common_wrfnmm_ensemble_mean_initialize(ens_mean)
    ens_member%west_east        = ens_mean%west_east
    ens_member%south_north      = ens_mean%south_north
    ens_member%bottom_top       = ens_mean%bottom_top
    ens_member%bottom_top_stag  = ens_mean%bottom_top_stag
    ens_member%soil_layers_stag = ens_mean%soil_layers_stag
    call common_wrfnmm_ensemble_member_initialize(ens_member)
    call common_wrfnmm_ensemble_mean_initialize(ens_vari)
    call common_wrfnmm_ensemble_variance_compute(ens_member,ens_mean,      &
         & ens_vari)

    ! Check local variable and proceed accordingly

    if(mpi_procid .eq. mpi_masternode) then

       ! Write variables to external file

       call ensemble_variance_wrfnmm_write(ens_vari)

    end if ! if(mpi_procid .eq. mpi_masternode)

    ! Enable the root task to catch up from I/O and calculations

    call mpi_barrier(mpi_comm_world,mpi_ierror)

    ! Deallocate memory for local variables

    call common_wrfnmm_ensemble_member_cleanup(ens_member)
    call common_wrfnmm_ensemble_member_cleanup(ens_mean)
    call common_wrfnmm_ensemble_member_cleanup(ens_vari)

    ! Enable the root task to catch up from I/O and calculations

    call mpi_barrier(mpi_comm_world,mpi_ierror)

    !=====================================================================

  end subroutine ensemble_variance_wrfnmm

  !=======================================================================

  ! ensemble_variance_wrfnmm_write.f90:

  !-----------------------------------------------------------------------

  subroutine ensemble_variance_wrfnmm_write(grid)

    ! Define variables passed to routine

    type(wrfnmm_variables)                                               :: grid

    !=====================================================================

    ! Define local variables

    ncstatus = nf90_open(path=trim(ensemble_variance_filename),            &
         & mode=nf90_write,ncid=ncfileid)
    ncstatus = nf90_inq_varid(ncfileid,'PINT',ncvarid)
    ncstatus = nf90_put_var(ncfileid,ncvarid,grid%pint)
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
    ncstatus = nf90_inq_varid(ncfileid,'F_ICE',ncvarid)
    ncstatus = nf90_put_var(ncfileid,ncvarid,grid%f_ice)
    ncstatus = nf90_inq_varid(ncfileid,'F_RAIN',ncvarid)
    ncstatus = nf90_put_var(ncfileid,ncvarid,grid%f_rain)
    ncstatus = nf90_inq_varid(ncfileid,'F_RIMEF',ncvarid)
    ncstatus = nf90_put_var(ncfileid,ncvarid,grid%f_rimef)
    ncstatus = nf90_inq_varid(ncfileid,'SMC',ncvarid)
    ncstatus = nf90_put_var(ncfileid,ncvarid,grid%smc)
    ncstatus = nf90_inq_varid(ncfileid,'STC',ncvarid)
    ncstatus = nf90_put_var(ncfileid,ncvarid,grid%stc)
    ncstatus = nf90_inq_varid(ncfileid,'PD',ncvarid)
    ncstatus = nf90_put_var(ncfileid,ncvarid,grid%pd)
    ncstatus = nf90_inq_varid(ncfileid,'TSK',ncvarid)
    ncstatus = nf90_put_var(ncfileid,ncvarid,grid%tsk)
    ncstatus = nf90_inq_varid(ncfileid,'SST',ncvarid)
    ncstatus = nf90_put_var(ncfileid,ncvarid,grid%sst)
    ncstatus = nf90_inq_varid(ncfileid,'SNO',ncvarid)
    ncstatus = nf90_put_var(ncfileid,ncvarid,grid%sno)
    ncstatus = nf90_inq_varid(ncfileid,'SICE',ncvarid)
    ncstatus = nf90_put_var(ncfileid,ncvarid,grid%sice)
    ncstatus = nf90_inq_varid(ncfileid,'U10',ncvarid)
    ncstatus = nf90_put_var(ncfileid,ncvarid,grid%u10)
    ncstatus = nf90_inq_varid(ncfileid,'V10',ncvarid)
    ncstatus = nf90_put_var(ncfileid,ncvarid,grid%v10)
    ncstatus = nf90_close(ncfileid)

    !=====================================================================

  end subroutine ensemble_variance_wrfnmm_write

  !=======================================================================

end module ensemble_variance_interface
