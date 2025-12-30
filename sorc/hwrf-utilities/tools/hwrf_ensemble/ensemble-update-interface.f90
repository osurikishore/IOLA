module ensemble_update_interface

  !=======================================================================

  ! Define associated modules and subroutines

  !-----------------------------------------------------------------------

  use kinds

  !-----------------------------------------------------------------------

  use common_interface
  use namelist
  use netcdf

  !-----------------------------------------------------------------------

  implicit none

  !-----------------------------------------------------------------------

contains

  !=======================================================================

  ! ensemble_update.f90:

  !-----------------------------------------------------------------------

  subroutine ensemble_update()

    !=====================================================================

    ! Compute local variables

    call ensemble_update_wrfnmm()

    !=====================================================================

  end subroutine ensemble_update

  !=======================================================================

  ! ensemble_update_wrfnmm.f90:

  !-----------------------------------------------------------------------

  subroutine ensemble_update_wrfnmm()

    ! Define variables computed within routine

    type(wrfnmm_prognostics)                                             :: ens_member

    !=====================================================================

    ! Define local variables

    call common_wrfnmm_ensemble_update_initialize(ens_member)

    ! Write external file

    call ensemble_update_wrfnmm_write(ens_member)

    ! Deallocate memory for local variables

    call common_wrfnmm_ensemble_recenter_cleanup(ens_member,ens_member)

    !=====================================================================

  end subroutine ensemble_update_wrfnmm

  !=======================================================================

  ! ensemble_update_wrfnmm_write.f90:

  !-----------------------------------------------------------------------

  subroutine ensemble_update_wrfnmm_write(grid)

    ! Define variables passed to routine

    type(wrfnmm_prognostics)                                             :: grid

    !=====================================================================

    ! Define local variables

    ncstatus = nf90_open(path=trim(ens_outfilename),mode=nf90_write,       &
         & ncid=ncfileid)
    ncstatus = nf90_inq_varid(ncfileid,'PINT',ncvarid)
    if(grid%is_exist_pint) ncstatus = nf90_put_var(ncfileid,ncvarid,       &
         & grid%pint)
    ncstatus = nf90_inq_varid(ncfileid,'T',ncvarid)
    if(grid%is_exist_t)   ncstatus = nf90_put_var(ncfileid,ncvarid,grid%t)
    ncstatus = nf90_inq_varid(ncfileid,'Q',ncvarid)
    if(grid%is_exist_q)   ncstatus = nf90_put_var(ncfileid,ncvarid,grid%q)
    ncstatus = nf90_inq_varid(ncfileid,'U',ncvarid)
    if(grid%is_exist_u)   ncstatus = nf90_put_var(ncfileid,ncvarid,grid%u)
    ncstatus = nf90_inq_varid(ncfileid,'V',ncvarid)
    if(grid%is_exist_v)   ncstatus = nf90_put_var(ncfileid,ncvarid,grid%v)
    ncstatus = nf90_inq_varid(ncfileid,'CWM',ncvarid)
    if(grid%is_exist_cwm) ncstatus = nf90_put_var(ncfileid,ncvarid,        &
         & grid%cwm)
    ncstatus = nf90_inq_varid(ncfileid,'PD',ncvarid)
    if(grid%is_exist_pd) ncstatus = nf90_put_var(ncfileid,ncvarid,         &
         & grid%pd)
    ncstatus = nf90_close(ncfileid)

    !=====================================================================

  end subroutine ensemble_update_wrfnmm_write

  !=======================================================================

end module ensemble_update_interface
