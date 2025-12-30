module bndrycnd_interface

  !=======================================================================

  !$$$ MODULE DOCUMENTATION BLOCK
  
  ! da-update-analysis :: bndrycnd_interface
  ! Copyright (C) 2018 Henry R. Winterbottom

  ! Email: henry.winterbottom@noaa.gov

  ! This program is free software: you can redistribute it and/or
  ! modify it under the terms of the GNU General Public License as
  ! published by the Free Software Foundation, either version 3 of the
  ! License, or (at your option) any later version.

  ! This program is distributed in the hope that it will be useful,
  ! but WITHOUT ANY WARRANTY; without even the implied warranty of
  ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  ! General Public License for more details.

  ! You should have received a copy of the GNU General Public License
  ! along with this program.  If not, see
  ! <http://www.gnu.org/licenses/>.

  !=======================================================================

  ! Define associated modules and subroutines

  use fileio_interface
  use namelist_interface
  use variable_interface

  ! Define interfaces and attributes for module routines
  
  implicit none
  private
  public :: bndrycnd_update
  interface bndrycnd_update
     module procedure hwrf_update
  end interface bndrycnd_update

  !-----------------------------------------------------------------------

contains

  !=======================================================================

  ! SUBROUTINE:

  ! compute_tendencies_hwrf.f90

  ! DESCRIPTION:

  ! This subroutine computes the updated boundary conditions and
  ! tendencies along each of the edges of the HWRF domain.

  ! INPUT VARIABLES:

  ! * hwrf_bdy_grid; a FORTRAN hwrf_bdy_grid_struct variable
  !   containing (at minimum) the analysis variable and the boundary
  !   update frequency.

  ! * hwrf_json_bdy; a FORTRAN hwrf_json_bdy_struct variable
  !   containing the user specifications regarding the boundary
  !   condition and tendency updates.

  ! OUTPUT VARIABLES:

  ! * hwrf_bdy_grid; a FORTRAN hwrf_bdy_struct variable containing the
  !   updated state variable and tendencies for the HWRF forecast
  !   domain.

  !-----------------------------------------------------------------------

  subroutine compute_tendencies_hwrf(hwrf_bdy_grid,hwrf_json_bdy)

    ! Define variables passed to routine

    type(hwrf_bdy_grid_struct)                                          :: hwrf_bdy_grid
    type(hwrf_json_bdy_struct)                                          :: hwrf_json_bdy

    !=====================================================================

    ! Check local variable and proceed accordingly

    if(hwrf_json_bdy%rezero) then

       ! Define local variables

       hwrf_bdy_grid%bxs  = 0.0
       hwrf_bdy_grid%bxe  = 0.0
       hwrf_bdy_grid%btxs = 0.0
       hwrf_bdy_grid%btxe = 0.0
       hwrf_bdy_grid%bys  = 0.0
       hwrf_bdy_grid%bye  = 0.0
       hwrf_bdy_grid%btys = 0.0
       hwrf_bdy_grid%btye = 0.0

    end if ! if(hwrf_json_bdy%rezero)

    ! Check local variable and proceed accordingly

    if(.not. hwrf_json_bdy%rezero) then

       ! Compute local variables

       call west_bdy_hwrf(hwrf_bdy_grid,hwrf_json_bdy)
       call east_bdy_hwrf(hwrf_bdy_grid,hwrf_json_bdy)
       call south_bdy_hwrf(hwrf_bdy_grid,hwrf_json_bdy)
       call north_bdy_hwrf(hwrf_bdy_grid,hwrf_json_bdy)

    end if ! if(.not. hwrf_json_bdy%rezero)

    !=====================================================================

  end subroutine compute_tendencies_hwrf

  !=======================================================================

  ! SUBROUTINE:

  ! east_bdy_hwrf.f90

  ! DESCRIPTION:

  ! This subroutine computes the HWRF forecast domain eastern-boundary
  ! tendencies and defines the state variable such that it is
  ! consistent with the user-specified initial conditions.

  ! INPUT VARIABLES:

  ! * hwrf_bdy_grid; a FORTRAN hwrf_bdy_grid_struct variable
  !   containing (at minimum) the analysis variable and the boundary
  !   update frequency.

  ! * hwrf_json_bdy; a FORTRAN hwrf_json_bdy_struct variable
  !   containing the user specifications regarding the boundary
  !   condition and tendency updates.

  ! OUTPUT VARIABLES:

  ! * hwrf_bdy_grid; a FORTRAN hwrf_bdy_struct variable containing the
  !   updated state variable and tendencies (bxe and btxe,
  !   respectively) for the HWRF forecast domain eastern boundary.

  !-----------------------------------------------------------------------

  subroutine east_bdy_hwrf(hwrf_bdy_grid,hwrf_json_bdy)

    ! Define variables passed to routine

    type(hwrf_bdy_grid_struct)                                          :: hwrf_bdy_grid
    type(hwrf_json_bdy_struct)                                          :: hwrf_json_bdy

    ! Define variables computed within routine

    logical                                                             :: write_diff
    real(r_kind)                                                        :: diff
    real(r_kind)                                                        :: var_ba
    real(r_kind)                                                        :: var_bt

    ! Define counting variables

    integer                                                             :: i, j, k

    !=====================================================================

    ! Loop through local variable

    do k = 1, hwrf_bdy_grid%nb

       ! Loop through local variable

       do j = 1, hwrf_bdy_grid%nz

          ! Loop through local variable

          do i = 1, hwrf_bdy_grid%ny

             ! Define local variables

             write_diff = .false.
             diff       = 0.0

             ! Check local variable and proceed accordingly
             
             if(hwrf_json_bdy%gridtype .eq. 1) then

                ! Check local variable and proceed accordingly
                
                if(mod(i,2) .eq. 1) then

                   ! Compute local variables

                   write_diff                                  = .true.
                   diff                                        =           &
                        & hwrf_bdy_grid%analvar(hwrf_bdy_grid%nx,i,j) -    &
                        & hwrf_bdy_grid%bxe(hwrf_bdy_grid%nx,i,j,1,k) 
                   var_ba                                      =           &
                        & hwrf_bdy_grid%analvar(hwrf_bdy_grid%nx,i,j)
                   var_bt                                      =           &
                        & hwrf_bdy_grid%bxe(hwrf_bdy_grid%nx,i,j,1,k) +    &
                        & hwrf_bdy_grid%bdyfrq*                            &
                        & hwrf_bdy_grid%btxe(hwrf_bdy_grid%nx,i,j,1,k)
                   hwrf_bdy_grid%btxe(hwrf_bdy_grid%nx,i,j,k,1) =          &
                        & (var_bt - var_ba)/hwrf_bdy_grid%bdyfrq
                   hwrf_bdy_grid%bxe(hwrf_bdy_grid%nx,i,j,k,1)  = var_ba

                end if ! if(mod(i,2) .eq. 1) 

             end if ! if(hwrf_json_bdy%gridtype .eq. 1)

             ! Check local variable and proceed accordingly
             
             if((hwrf_json_bdy%gridtype .eq. 2) .or.                       &
                  & (hwrf_json_bdy%gridtype .eq. 3)) then

                ! Check local variable and proceed accordingly
                
                if(mod(i,2) .eq. 0) then

                   ! Compute local variables
                
                   write_diff                                   = .true.
                   diff                                         =          &
                        & hwrf_bdy_grid%analvar(hwrf_bdy_grid%nx,i,j) -    &
                        & hwrf_bdy_grid%bxe(hwrf_bdy_grid%nx,i,j,1,k)
                   var_ba                                       =          &
                        & hwrf_bdy_grid%analvar(hwrf_bdy_grid%nx,i,j)
                   var_bt                                       =          &
                        & hwrf_bdy_grid%bxe(hwrf_bdy_grid%nx,i,j,1,k) +    &
                        & hwrf_bdy_grid%bdyfrq*                            &
                        & hwrf_bdy_grid%btxe(hwrf_bdy_grid%nx,i,j,1,k)
                   hwrf_bdy_grid%btxe(hwrf_bdy_grid%nx,i,j,k,1) =          &
                        & (var_bt - var_ba)/hwrf_bdy_grid%bdyfrq
                   hwrf_bdy_grid%bxe(hwrf_bdy_grid%nx,i,j,k,1)  = var_ba

                endif ! if(mod(i,2) .eq. 1)

             end if ! if((hwrf_json_bdy%gridtype .eq. 2)
                    ! .or. (hwrf_json_bdy%gridtype .eq. 3))

             ! Check local variable and proceed accordingly

             if(debug .and. write_diff) write(6,500)                       &
                  & trim(adjustl(hwrf_json_bdy%variable_name)), i, j, k,   &
                  & diff

          end do ! do i = 1, hwrf_bdy_grid%ny

       end do ! do j = 1, hwrf_bdy_grid%nz

    end do ! do k = 1, hwrf_bdy_grid%nb

    ! Define local variables

500 format('HWRF ', a,' EAST BOUNDARY DIFFERENCE (i,j,k): ',i4,1x,i4,1x,   &
         & i4,1x,f13.5)

    !=====================================================================

  end subroutine east_bdy_hwrf

  !=======================================================================

  ! SUBROUTINE:

  ! hwrf_update.f90

  ! DESCRIPTION:

  ! This subroutine is the driver routine to compute and update HWRF
  ! boundary conditions and tendencies in accordance with the user
  ! specifications.

  ! INPUT VARIABLES:

  ! * hwrf_json_bdy; a FORTRAN hwrf_json_bdy_struct variable
  !   containing the user specifications regarding the boundary
  !   condition and tendency updates.

  !-----------------------------------------------------------------------

  subroutine hwrf_update(hwrf_json_bdy)

    ! Define variables passed to routine

    type(hwrf_json_bdy_struct)                                          :: hwrf_json_bdy

    ! Define variables computed within routine

    type(hwrf_bdy_grid_struct)                                          :: hwrf_bdy_grid
    type(netcdf_var_struct)                                             :: netcdf_var
    type(nmm_json_var_struct)                                           :: json_var

    !=====================================================================

    ! Define local variables

    json_var%variable_name = hwrf_json_bdy%variable_name
    json_var%gridtype      = hwrf_json_bdy%gridtype
    json_var%levtype       = hwrf_json_bdy%levtype
    call fileio_interface_read(hwrf_an_filename,json_var,netcdf_var)
    call fileio_interface_read(hwrf_bdy_filename,json_var,hwrf_bdy_grid)
    hwrf_bdy_grid%analvar  = netcdf_var%var(:,:,:,1)

    ! Deallocate memory for local variables

    call variable_interface_cleanup_struct(netcdf_var)

    ! Compute local variables

    call compute_tendencies_hwrf(hwrf_bdy_grid,hwrf_json_bdy)

    ! Define local variables

    call fileio_interface_write(hwrf_bdy_filename,hwrf_json_bdy,           &
         & hwrf_bdy_grid)

    ! Deallocate memory for local variables

    call variable_interface_cleanup_struct(hwrf_bdy_grid)

    !=====================================================================

  end subroutine hwrf_update

  !=======================================================================

  ! SUBROUTINE:

  ! north_bdy_hwrf.f90

  ! DESCRIPTION:

  ! This subroutine computes the HWRF forecast domain
  ! northern-boundary tendencies and defines the state variable such
  ! that it is consistent with the user-specified initial conditions.

  ! INPUT VARIABLES:

  ! * hwrf_bdy_grid; a FORTRAN hwrf_bdy_grid_struct variable
  !   containing (at minimum) the analysis variable and the boundary
  !   update frequency.

  ! * hwrf_json_bdy; a FORTRAN hwrf_json_bdy_struct variable
  !   containing the user specifications regarding the boundary
  !   condition and tendency updates.

  ! OUTPUT VARIABLES:

  ! * hwrf_bdy_grid; a FORTRAN hwrf_bdy_struct variable containing the
  !   updated state variable and tendencies (bye and btye,
  !   respectively) for the HWRF forecast domain northern boundary.

  !-----------------------------------------------------------------------

  subroutine north_bdy_hwrf(hwrf_bdy_grid,hwrf_json_bdy)

    ! Define variables passed to routine

    type(hwrf_bdy_grid_struct)                                          :: hwrf_bdy_grid
    type(hwrf_json_bdy_struct)                                          :: hwrf_json_bdy

    ! Define variables computed within routine

    real(r_kind)                                                        :: diff
    real(r_kind)                                                        :: var_ba
    real(r_kind)                                                        :: var_bt

    ! Define counting variables

    integer                                                             :: i, j, k

    !=====================================================================

    ! Loop through local variable

    do k = 1, hwrf_bdy_grid%nb

       ! Loop through local variable

       do j = 1, hwrf_bdy_grid%nz

          ! Loop through local variable

          do i = 1, hwrf_bdy_grid%nx

             ! Compute local variables

             diff                                        =                 &
                  & hwrf_bdy_grid%analvar(i,hwrf_bdy_grid%ny,j) -          &
                  & hwrf_bdy_grid%bye(i,hwrf_bdy_grid%ny,j,1,k)
             var_ba                                      =                 &
                  & hwrf_bdy_grid%analvar(i,hwrf_bdy_grid%ny,j)
             var_bt                                      =                 &
                  & hwrf_bdy_grid%bye(i,hwrf_bdy_grid%ny,j,1,k) +          &
                  & hwrf_bdy_grid%bdyfrq*                                  &
                  & hwrf_bdy_grid%btye(i,hwrf_bdy_grid%ny,j,1,k)
             hwrf_bdy_grid%btye(i,hwrf_bdy_grid%ny,j,k,1) =                &
                  & (var_bt - var_ba)/hwrf_bdy_grid%bdyfrq
             hwrf_bdy_grid%bye(i,hwrf_bdy_grid%ny,j,k,1)  = var_ba

             ! Define local variables

             if(debug) write(6,500)                                        &
                  & trim(adjustl(hwrf_json_bdy%variable_name)), i, j, k,   &
                  & diff

          end do ! do i = 1, hwrf_bdy_grid%nx

       end do ! do j = 1, hwrf_bdy_grid%nz

    end do ! do k = 1, hwrf_bdy_grid%nb

    ! Define local variables

500 format('HWRF ', a,' NORTH BOUNDARY DIFFERENCE (i,j,k): ',i4,1x,i4,1x,  &
         & i4,1x,f13.5)

    !=====================================================================

  end subroutine north_bdy_hwrf

  !=======================================================================

  ! SUBROUTINE:

  ! south_bdy_hwrf.f90

  ! DESCRIPTION:

  ! This subroutine computes the HWRF forecast domain
  ! southern-boundary tendencies and defines the state variable such
  ! that it is consistent with the user-specified initial conditions.

  ! INPUT VARIABLES:

  ! * hwrf_bdy_grid; a FORTRAN hwrf_bdy_grid_struct variable
  !   containing (at minimum) the analysis variable and the boundary
  !   update frequency.

  ! * hwrf_json_bdy; a FORTRAN hwrf_json_bdy_struct variable
  !   containing the user specifications regarding the boundary
  !   condition and tendency updates.

  ! OUTPUT VARIABLES:

  ! * hwrf_bdy_grid; a FORTRAN hwrf_bdy_struct variable containing the
  !   updated state variable and tendencies (bys and btys,
  !   respectively) for the HWRF forecast domain southern boundary.

  !-----------------------------------------------------------------------

  subroutine south_bdy_hwrf(hwrf_bdy_grid,hwrf_json_bdy)

    ! Define variables passed to routine

    type(hwrf_bdy_grid_struct)                                          :: hwrf_bdy_grid
    type(hwrf_json_bdy_struct)                                          :: hwrf_json_bdy

    ! Define variables computed within routine

    real(r_kind)                                                        :: diff
    real(r_kind)                                                        :: var_ba
    real(r_kind)                                                        :: var_bt

    ! Define counting variables

    integer                                                             :: i, j, k

    !=====================================================================

    ! Loop through local variable

    do k = 1, hwrf_bdy_grid%nb

       ! Loop through local variable

       do j = 1, hwrf_bdy_grid%nz

          ! Loop through local variable
          
          do i = 1, hwrf_bdy_grid%nx

             ! Compute local variables

             diff                          = hwrf_bdy_grid%analvar(i,1,j)  &
                  & - hwrf_bdy_grid%bys(i,1,j,1,k)
             var_ba                        =                               &
                  & hwrf_bdy_grid%analvar(i,1,j)
             var_bt                        =                               &
                  & hwrf_bdy_grid%bys(i,1,j,1,k)                           &
                  & + hwrf_bdy_grid%bdyfrq*                                &
                  & hwrf_bdy_grid%btys(i,1,j,1,k)
             hwrf_bdy_grid%btys(i,1,j,k,1) = (var_bt - var_ba)/            &
                  & hwrf_bdy_grid%bdyfrq
             hwrf_bdy_grid%bys(i,1,j,k,1)  = var_ba

             ! Define local variables

             if(debug) write(6,500)                                        &
                  & trim(adjustl(hwrf_json_bdy%variable_name)), i, j, k,   &
                  & diff

          end do ! do i = 1, hwrf_bdy_grid%nx

       end do ! do j = 1, hwrf_bdy_grid%nz

    end do ! do k = 1, hwrf_bdy_grid%nb

    ! Define local variables

500 format('HWRF ', a,' SOUTH BOUNDARY DIFFERENCE (i,j,k): ',i4,1x,i4,1x,  &
         & i4,1x,f13.5)

    !=====================================================================

  end subroutine south_bdy_hwrf

  !=======================================================================

  ! SUBROUTINE:

  ! west_bdy_hwrf.f90

  ! DESCRIPTION:

  ! This subroutine computes the HWRF forecast domain western-boundary
  ! tendencies and defines the state variable such that it is
  ! consistent with the user-specified initial conditions.

  ! INPUT VARIABLES:

  ! * hwrf_bdy_grid; a FORTRAN hwrf_bdy_grid_struct variable
  !   containing (at minimum) the analysis variable and the boundary
  !   update frequency.

  ! * hwrf_json_bdy; a FORTRAN hwrf_json_bdy_struct variable
  !   containing the user specifications regarding the boundary
  !   condition and tendency updates.

  ! OUTPUT VARIABLES:

  ! * hwrf_bdy_grid; a FORTRAN hwrf_bdy_struct variable containing the
  !   updated state variable and tendencies (bxs and btxs,
  !   respectively) for the HWRF forecast domain northern boundary.

  !-----------------------------------------------------------------------

  subroutine west_bdy_hwrf(hwrf_bdy_grid,hwrf_json_bdy)

    ! Define variables passed to routine

    type(hwrf_bdy_grid_struct)                                          :: hwrf_bdy_grid
    type(hwrf_json_bdy_struct)                                          :: hwrf_json_bdy

    ! Define variables computed within routine

    logical                                                             :: write_diff
    real(r_kind)                                                        :: diff
    real(r_kind)                                                        :: var_ba
    real(r_kind)                                                        :: var_bt

    ! Define counting variables

    integer                                                             :: i, j, k

    !=====================================================================

    ! Loop through local variable

    do k = 1, hwrf_bdy_grid%nb

       ! Loop through local variable

       do j = 1, hwrf_bdy_grid%nz

          ! Loop through local variable
          
          do i = 1, hwrf_bdy_grid%ny

             ! Define local variables

             write_diff = .false.
             diff       = 0.0
             
             ! Check local variable and proceed accordingly
             
             if(hwrf_json_bdy%gridtype .eq. 1) then

                ! Check local variable and proceed accordingly
                
                if(mod(i,2) .eq. 1) then

                   ! Compute local variables
                
                   write_diff                                  = .true.
                   diff                                        =           &
                        & hwrf_bdy_grid%analvar(1,i,j) -                   &
                        & hwrf_bdy_grid%bxs(1,i,j,1,k)
                   var_ba                                      =           &
                        & hwrf_bdy_grid%analvar(1,i,j)
                   var_bt                                      =           &
                     & hwrf_bdy_grid%bxs(1,i,j,1,k) +                      &
                     & hwrf_bdy_grid%bdyfrq*hwrf_bdy_grid%btxs(1,i,j,1,k)
                   hwrf_bdy_grid%btxs(1,i,j,k,1)               =           &
                        & (var_bt - var_ba)/hwrf_bdy_grid%bdyfrq
                   hwrf_bdy_grid%bxs(1,i,j,k,1)                = var_ba

                endif ! if(mod(i,2) .eq. 1)

             end if ! if(hwrf_json_bdy%gridtype .eq. 1)

             ! Check local variable and proceed accordingly
             
             if((hwrf_json_bdy%gridtype .eq. 2) .or.                       &
                  & (hwrf_json_bdy%gridtype .eq. 3)) then

                ! Check local variable and proceed accordingly
                
                if(mod(i,2) .eq. 0) then

                   ! Compute local variables

                   write_diff                                  = .true.                
                   diff                                        =           &
                        & hwrf_bdy_grid%analvar(1,i,j) -                   &
                        & hwrf_bdy_grid%bxs(1,i,j,1,k)
                   var_ba                                      =           &
                        & hwrf_bdy_grid%analvar(1,i,j)
                   var_bt                                      =           &
                        & hwrf_bdy_grid%bxs(1,i,j,1,k) +                   &
                        & hwrf_bdy_grid%bdyfrq*                            &
                        & hwrf_bdy_grid%btxs(1,i,j,1,k)
                   hwrf_bdy_grid%btxs(1,i,j,k,1)               =           &
                        & (var_bt - var_ba)/hwrf_bdy_grid%bdyfrq
                   hwrf_bdy_grid%bxs(1,i,j,k,1)                = var_ba

                endif ! if(mod(i,2) .eq. 1)

             end if ! if((hwrf_json_bdy%gridtype .eq. 2)
                    ! .or. (hwrf_json_bdy%gridtype .eq. 3) then

             ! Check local variable and proceed accordingly

             if(debug .and. write_diff) write(6,500)                       &
                  & trim(adjustl(hwrf_json_bdy%variable_name)), i, j, k,   &
                  & diff
             
          end do ! do i = 1, hwrf_bdy_grid%ny

       end do ! do j = 1, hwrf_bdy_grid%nz

    end do ! do k = 1, hwrf_bdy_grid%nb

    ! Define local variables

500 format('HWRF ', a,' WEST BOUNDARY DIFFERENCE (i,j,k): ',i4,1x,i4,1x,   &
         & i4,1x,f13.5)

    !=====================================================================

  end subroutine west_bdy_hwrf

  !=======================================================================

end module bndrycnd_interface
