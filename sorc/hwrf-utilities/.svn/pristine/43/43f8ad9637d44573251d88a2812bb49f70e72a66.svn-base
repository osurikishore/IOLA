module wrf_ensemble_interface

  !=======================================================================

  ! Define associated modules and subroutines

  !-----------------------------------------------------------------------

  use kinds

  !-----------------------------------------------------------------------

  use ensemble_inflate_interface
  use ensemble_mean_interface
  use ensemble_recenter_interface
  use ensemble_variance_interface
  use ensemble_update_interface
  use namelist

  !-----------------------------------------------------------------------

  implicit none

  !-----------------------------------------------------------------------

contains

  !=======================================================================

  ! wrf_ensemble.f90:

  !-----------------------------------------------------------------------

  subroutine wrf_ensemble()

    !=====================================================================

    ! Define local variables

    call namelistparams()

    ! Check local variable and proceed accordingly

    if(is_ensemble_mean)     call ensemble_mean()
    if(is_ensemble_variance) call ensemble_variance()
    if(is_ensemble_recenter) call ensemble_recenter()
    if(is_ensemble_inflate)  call ensemble_inflate()
    if(is_ensemble_update)   call ensemble_update()

    !=====================================================================

  end subroutine wrf_ensemble

  !=======================================================================

end module wrf_ensemble_interface
