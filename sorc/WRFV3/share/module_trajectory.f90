   
   module module_trajectory

   use module_driver_constants,  only : max_domains
   CONTAINS

   subroutine trajectory_init( grid, config_flags, &
                               ims,ime, jms,jme, kms,kme )

   use module_domain
   use module_configure,         only : grid_config_rec_type




   integer, intent(in)      :: ims,ime, jms,jme, kms,kme
   type(domain), intent(inout)            :: grid
   type(grid_config_rec_type), intent(in) :: config_flags

   end subroutine trajectory_init

   subroutine trajectory_driver( grid )

   use module_domain




   type(domain), intent(in) :: grid

   end subroutine trajectory_driver


   end module module_trajectory
