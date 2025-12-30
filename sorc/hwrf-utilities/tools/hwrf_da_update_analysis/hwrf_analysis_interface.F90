module hwrf_analysis_interface

  !=======================================================================

  !$$$ MODULE DOCUMENTATION BLOCK
  
  ! da-update-analysis :: hwrf_analysis_interface
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

  use analysis_interface
  use bndrycnd_interface
  use constants_interface
  use fileio_interface
  use interpolation_interface
  use namelist_interface
  use nmm_interpolation_interface
  use variable_interface

  ! Define interfaces and attributes for module routines
  
  implicit none
  private
  public :: hwrf_analysis

  !-----------------------------------------------------------------------

contains

  !=======================================================================

  ! SUBROUTINE:

  ! hwrf_analysis.f90

  ! DESCRIPTION:

  ! This is the driver routine for HWRF analysis variables actions to
  ! be accomplished in accordance with the user specifications.

  !-----------------------------------------------------------------------

  subroutine hwrf_analysis()

    !=====================================================================

    ! Compute local variables

    if(is_bdy_update) call hwrf_bdycnd_update()
    if(is_da)         call hwrf_da_analysis()
    if(is_ghst)       call hwrf_ghost_analysis()
    if(is_nest)       call hwrf_nest_analysis()
    if(is_nest_cycle) call hwrf_cycle_nest_analysis()

    !=====================================================================

  end subroutine hwrf_analysis

  !=======================================================================

  ! SUBROUTINE:

  ! hwrf_bdycnd_update.f90

  ! DESCRIPTION:

  ! This is the driver routine to computed updated HWRF boundary
  ! conditions in accordance with the user specifications.

  !-----------------------------------------------------------------------

  subroutine hwrf_bdycnd_update()

    ! Define variables computed within routine

    type(hwrf_json_bdy_struct), dimension(:),               allocatable :: hwrf_json_bdy
    integer                                                             :: nvars

    ! Define counting variables

    integer                                                             :: i

    !=====================================================================

    ! Define local variables

    call fileio_interface_read(hwrf_bdy_json_vtable,nvars,hwrf_json_bdy)

    ! Loop through local variable

    do i = 1, size(hwrf_json_bdy)

       ! Compute local variables

       call bndrycnd_update(hwrf_json_bdy(i))

    end do ! do i = 1, size(hwrf_json_bdy)
 
    ! Deallocate memory for local variables

    if(allocated(hwrf_json_bdy)) deallocate(hwrf_json_bdy)

    !=====================================================================

  end subroutine hwrf_bdycnd_update

  !=======================================================================

  ! SUBROUTINE:

  ! hwrf_cycle_nest_analysis.f90

  ! DESCRIPTION:

  ! This subroutine is the driver routine to cycle HWRF nested grids
  ! from (subsequent) data-assimilation analysis cycles.

  ! NOTE:

  ! This routine is not yet supported within the operational HWRF
  ! workflow system; use with caution!!

  !-----------------------------------------------------------------------

  subroutine hwrf_cycle_nest_analysis()

    ! Define variables computed within routine

    type(nmm_json_var_struct),  dimension(:),               allocatable :: nmm_json_in
    type(slint_struct)                                                  :: slint(3)
    type(merge_struct)                                                  :: merge
    type(grid_struct)                                                   :: grid_an
    type(grid_struct)                                                   :: grid_fg
    type(hwrf_vert_grid_struct)                                         :: hwrf_vert_grid
    type(netcdf_var_struct)                                             :: netcdf_an
    type(netcdf_var_struct)                                             :: netcdf_fg
    type(nmm_json_grid_struct)                                          :: nmm_json_grid_an
    type(nmm_json_grid_struct)                                          :: nmm_json_grid_fg
    type(nmm_json_var_struct)                                           :: nmm_json_out
    logical                                                             :: mgrdn
    logical                                                             :: ugrdn
    logical                                                             :: vgrdn
    integer                                                             :: gridtype
    integer                                                             :: intrptype
    integer                                                             :: nvars

    ! Define counting variables

    integer                                                             :: i, j

    !=====================================================================

    ! Define local variables

    nmm_in_filename  = hwrf_fg_filename
    nmm_out_filename = hwrf_an_filename
    call fileio_interface_read(hwrf_grid_json_vtable,nmm_json_grid_fg,    &
         & nmm_json_grid_an) 
    call nmm_setup_nest_interpolation(nmm_json_grid_fg,nmm_json_grid_an,  &
         & slint)
    call fileio_interface_read(hwrf_var_json_vtable,nvars,nmm_json_in)
    if(is_slint) hwrf_vert_grid%ncoords = slint(1)%dst_ncoords

    ! Loop through local variable

    do i = 1, nvars

       ! Check local variable and proceed accordingly

       if(trim(adjustl(nmm_json_in(i)%variable_name)) .eq. 'PINT') then

          ! Define local variables

          gridtype                   = nmm_json_in(i)%gridtype
          intrptype                  = nmm_json_in(i)%intrptype
          call fileio_interface_read(hwrf_fg_filename,nmm_json_in(i),     &
               & netcdf_fg)
          nmm_json_out               = nmm_json_in(i)
          nmm_json_out%variable_name =                                    &
               & nmm_json_in(i)%remap_variable_name
          call fileio_interface_read(hwrf_an_filename,nmm_json_out,       &
               & netcdf_an)
          mgrdn                      = .false.
          ugrdn                      = .false.
          vgrdn                      = .false.
          if(gridtype .eq. 1) mgrdn = .true.
          if(gridtype .eq. 2) ugrdn = .true.
          if(gridtype .eq. 3) vgrdn = .true.
          grid_fg%nz                 = netcdf_fg%nz
          grid_an%nz                 = netcdf_an%nz
          call nmm_setup_grid(hwrf_fg_filename,nmm_json_grid_fg,grid_fg,  &
               & mgrd=mgrdn,ugrd=ugrdn,vgrd=vgrdn)
          call nmm_setup_grid(hwrf_an_filename,nmm_json_grid_an,grid_an,  &
               & mgrd=mgrdn,ugrd=ugrdn,vgrd=vgrdn)
          hwrf_vert_grid%nx          = grid_an%nx
          hwrf_vert_grid%ny          = grid_an%ny
          hwrf_vert_grid%nz          = (grid_an%nz - 1)
          call variable_interface_setup_struct(hwrf_vert_grid)

          ! Loop through local variable

          do j = 1, grid_fg%nz

             ! Define local variables

             grid_fg%var(:,j)   = reshape(netcdf_fg%var(:,:,j,1),         &
                  & shape(grid_fg%var(:,j)))
             grid_an%var(:,j) = reshape(netcdf_an%var(:,:,j,1),           &
                  & shape(grid_an%var(:,j)))
          
             ! Check local variable and proceed accordingly

             if(is_slint) then

                ! Compute local variables

                call nmm_nest_interpolation(nmm_json_in(i),               &
                     & slint(gridtype),grid_an%var(:,j),grid_fg%var(:,j))

             end if ! if(is_slint)

             ! Define local variables

             where(grid_an%var(:,j) .eq. spval) grid_an%var(:,j)          &
                  & = reshape(netcdf_an%var(:,:,j,1),                     &
                  & shape(grid_an%var(:,j)))

          end do ! do j = 1, grid_fg%nz

          ! Define local variables

          hwrf_vert_grid%nest_pint = grid_an%var

          ! Deallocate memory for local variables

          call variable_interface_cleanup_struct(netcdf_fg)
          call variable_interface_cleanup_struct(netcdf_an)
          call variable_interface_cleanup_struct(grid_fg)
          call variable_interface_cleanup_struct(grid_an)

       end if ! if(trim(adjustl(nmm_json_in(i)%variable_name))
              ! .eq. 'PINT')

    end do ! do i = 1, nvars

    ! Loop through local variable

    do i = 1, nvars

       ! Check local variable and proceed accordingly

       if(trim(adjustl(nmm_json_in(i)%variable_name)) .eq. 'PD') then

          ! Define local variables

          gridtype                   = nmm_json_in(i)%gridtype
          intrptype                  = nmm_json_in(i)%intrptype
          call fileio_interface_read(hwrf_fg_filename,nmm_json_in(i),     &
               & netcdf_fg)
          nmm_json_out               = nmm_json_in(i)
          nmm_json_out%variable_name =                                    &
               & nmm_json_in(i)%remap_variable_name
          call fileio_interface_read(hwrf_an_filename,nmm_json_out,       &
               & netcdf_an)
          mgrdn                      = .false.
          ugrdn                      = .false.
          vgrdn                      = .false.
          if(gridtype .eq. 1) mgrdn = .true.
          if(gridtype .eq. 2) ugrdn = .true.
          if(gridtype .eq. 3) vgrdn = .true.
          grid_an%nz                 = netcdf_fg%nz
          grid_fg%nz                 = netcdf_an%nz
          call nmm_setup_grid(hwrf_fg_filename,nmm_json_grid_fg,grid_fg,  &
               & mgrd=mgrdn,ugrd=ugrdn,vgrd=vgrdn)
          call nmm_setup_grid(hwrf_an_filename,nmm_json_grid_an,grid_an,  &
               & mgrd=mgrdn,ugrd=ugrdn,vgrd=vgrdn)
          merge%nx                   = netcdf_an%nx
          merge%ny                   = netcdf_an%ny
          merge%nz                   = netcdf_an%nz
          call variable_interface_setup_struct(merge)
          merge%npts_relax           = npts_relax
          merge%lat_an               = grid_an%lat
          merge%lon_an               = grid_an%lon
          merge%lat_fg               = grid_fg%lat
          merge%lon_fg               = grid_fg%lon
          merge%dx                   = grid_fg%dx
          merge%dy                   = grid_fg%dy

          ! Loop through local variable

          do j = 1, grid_fg%nz

             ! Define local variables

             grid_fg%var(:,j)  = reshape(netcdf_fg%var(:,:,j,1),          &
                  & shape(grid_fg%var(:,j)))
             grid_an%var(:,j)  = reshape(netcdf_an%var(:,:,j,1),          &
                  & shape(grid_an%var(:,j)))
             merge%var_fg(:,j) = grid_an%var(:,j)
          
             ! Check local variable and proceed accordingly

             if(is_slint) then

                ! Compute local variables

                call nmm_nest_interpolation(nmm_json_in(i),               &
                     & slint(gridtype),grid_fg%var(:,j),grid_an%var(:,j))

             end if ! if(is_slint)

             ! Define local variables

             where(grid_an%var(:,j) .eq. spval) grid_an%var(:,j)          &
                  & = reshape(netcdf_an%var(:,:,j,1),                     &
                  & shape(grid_an%var(:,j)))
             merge%var_an(:,j) = grid_an%var(:,j)

          end do ! do j = 1, grid_nest%nz

          ! Compute local variables

          call blend_interface_cycling(merge)
          
          ! Define local variables

          hwrf_vert_grid%pd = merge%var(:,1)
          
          ! Deallocate memory for local variables

          call variable_interface_cleanup_struct(grid_fg)
          call variable_interface_cleanup_struct(grid_an)
          call variable_interface_cleanup_struct(merge)
          call variable_interface_cleanup_struct(netcdf_fg)
          call variable_interface_cleanup_struct(netcdf_an)
          
       end if ! if(trim(adjustl(nmm_json_in(i)%variable_name))
              ! .eq. 'PD')

    end do ! do i = 1, nvars

    ! Compute local variables

    call nmm_vertical_levels(hwrf_an_filename,hwrf_vert_grid)

    ! Define local variables

    call fileio_interface_write(hwrf_an_filename,hwrf_vert_grid)

    ! Loop through local variable

    do i = 1, nvars

       ! Check local variable and proceed accordingly
       
       if((trim(adjustl(nmm_json_in(i)%variable_name)) .ne. 'PD') .and.   &
            & (trim(adjustl(nmm_json_in(i)%variable_name)) .ne. 'PINT'))  &
            & then

          ! Define local variables
          
          gridtype                   = nmm_json_in(i)%gridtype
          intrptype                  = nmm_json_in(i)%intrptype
          call fileio_interface_read(hwrf_fg_filename,nmm_json_in(i),     &
               & netcdf_fg)
          nmm_json_out               = nmm_json_in(i)
          nmm_json_out%variable_name = nmm_json_in(i)%remap_variable_name
          call fileio_interface_read(hwrf_an_filename,nmm_json_out,       &
               & netcdf_an)
          mgrdn                      = .false.
          ugrdn                      = .false.
          vgrdn                      = .false.
          if(gridtype .eq. 1) mgrdn = .true.
          if(gridtype .eq. 2) ugrdn = .true.
          if(gridtype .eq. 3) vgrdn = .true.
          grid_fg%nz                 = netcdf_fg%nz
          grid_an%nz                 = netcdf_an%nz
          call nmm_setup_grid(hwrf_fg_filename,nmm_json_grid_fg,grid_fg,  &
               & mgrd=mgrdn,ugrd=ugrdn,vgrd=vgrdn)
          call nmm_setup_grid(hwrf_an_filename,nmm_json_grid_an,grid_an,  &
               & mgrd=mgrdn,ugrd=ugrdn,vgrd=vgrdn)
          merge%nx                   = netcdf_an%nx
          merge%ny                   = netcdf_an%ny
          merge%nz                   = netcdf_an%nz
          call variable_interface_setup_struct(merge)
          merge%npts_blend           = npts_blend
          merge%npts_relax           = npts_relax
          merge%lat_an               = grid_an%lat
          merge%lon_an               = grid_an%lon
          merge%lat_fg               = grid_fg%lat
          merge%lon_fg               = grid_fg%lon
          merge%dx                   = grid_fg%dx
          merge%dy                   = grid_fg%dy

          ! Loop through local variable
          
          do j = 1, grid_fg%nz
             
             ! Define local variables
             
             grid_fg%var(:,j)  = reshape(netcdf_fg%var(:,:,j,1),          &
                  & shape(grid_fg%var(:,j)))
             grid_an%var(:,j)  = reshape(netcdf_an%var(:,:,j,1),          &
                  & shape(grid_an%var(:,j)))
             merge%var_fg(:,j) = grid_fg%var(:,j)
          
             ! Check local variable and proceed accordingly
             
             if(is_slint) then
                
                ! Compute local variables
                
                call nmm_nest_interpolation(nmm_json_in(i),               &
                     & slint(gridtype),grid_fg%var(:,j),                  &
                     & grid_an%var(:,j))

             end if ! if(is_slint)

             ! Define local variables

             where(grid_an%var(:,j) .eq. spval) grid_an%var(:,j)          &
                  & = reshape(netcdf_an%var(:,:,j,1),                     &
                  & shape(grid_an%var(:,j)))
             merge%var_an(:,j) = grid_an%var(:,j)
             
          end do ! do j = 1, grid_fg%nz
          
          ! Compute local variables
          
          call blend_interface_cycling(merge)

          ! Check local variable and proceed accordingly

          if(nmm_json_in(i)%levtype .eq. 3) then

             ! Define local variables

             hwrf_vert_grid%var(:,1:merge%nz) = merge%var

             ! Compute local variables

             call nmm_vertical_interpolation(nmm_json_in(i),              &
                  & hwrf_vert_grid)

             ! Define local variables

             merge%var = hwrf_vert_grid%var(:,1:merge%nz)

          end if ! if(nmm_json_in(i)%levtype .eq. 3)

          ! Loop through local variable

          do j = 1, grid_an%nz

             ! Define local variables
             
             netcdf_an%var(:,:,j,1) = reshape(merge%var(:,j),             &
                  & shape(netcdf_an%var(:,:,j,1)))

          end do ! do j = 1, grid_an%nz

          ! Define local variables
          
          call fileio_interface_write(hwrf_an_filename,nmm_json_out,      &
               & netcdf_an)

          ! Deallocate memory for local variables
          
          call variable_interface_cleanup_struct(grid_an)
          call variable_interface_cleanup_struct(grid_fg)
          call variable_interface_cleanup_struct(merge)
          call variable_interface_cleanup_struct(netcdf_an)
          call variable_interface_cleanup_struct(netcdf_fg)

       end if ! if((trim(adjustl(nmm_json_in(i)%variable_name))
              ! .ne. 'PD')
              ! .and. (trim(adjustl(nmm_json_in(i)%variable_name))
              ! .ne. 'PINT'))

    end do ! do i = 1, nvars

    ! Deallocate memory for local variables

    call interpolation_interface_cleanup(slint)

    !=====================================================================

  end subroutine hwrf_cycle_nest_analysis

  !=======================================================================

  ! SUBROUTINE:

  ! hwrf_da_analysis.f90

  ! DESCRIPTION:

  ! This is the driver routine to blend HWRF analysis increment
  ! variables, computed by a data-assimilation (DA) system, relative
  ! to the location of tropical cyclones (TC) within the analysis
  ! fields.

  ! NOTE:

  ! The the Grid-point Statistical Interpolation (GSI) algorithm is
  ! the only DA system supported by this software at the moment.

  !-----------------------------------------------------------------------

  subroutine hwrf_da_analysis()

    ! Define variables computed within routine

    type(nmm_json_var_struct),  dimension(:),               allocatable :: nmm_json_in
    type(nmm_json_grid_struct)                                          :: nmm_json_grid_in
    type(nmm_json_grid_struct)                                          :: nmm_json_grid_out
    type(nmm_json_var_struct)                                           :: nmm_json_out
    type(netcdf_var_struct)                                             :: netcdf_fg
    type(netcdf_var_struct)                                             :: netcdf_an
    type(grid_struct)                                                   :: grid_an
    type(grid_struct)                                                   :: grid_fg
    type(grid_struct)                                                   :: grid_ic
    logical                                                             :: mgrdn
    logical                                                             :: ugrdn
    logical                                                             :: vgrdn
    integer                                                             :: nvars

    ! Define counting variables

    integer                                                             :: i, j

    !=====================================================================

    ! Define local variables

    nmm_in_filename  = hwrf_an_filename
    nmm_out_filename = hwrf_fg_filename
    call fileio_interface_read(hwrf_grid_json_vtable,nmm_json_grid_in,    &
         & nmm_json_grid_out)    
    call fileio_interface_read(hwrf_var_json_vtable,nvars,nmm_json_in)

    ! Loop through local variable

    do i = 1, nvars

       ! Define local variables

       call fileio_interface_read(hwrf_an_filename,nmm_json_in(i),        &
            & netcdf_an)
       call fileio_interface_read(hwrf_fg_filename,nmm_json_in(i),        &
            & netcdf_fg)
       nmm_json_out               = nmm_json_in(i)
       nmm_json_out%variable_name = nmm_json_in(i)%remap_variable_name
       grid_an%nx                 = netcdf_an%nx
       grid_an%ny                 = netcdf_an%ny
       grid_an%nz                 = netcdf_an%nz
       call variable_interface_setup_struct(grid_an)
       grid_fg%nx                 = netcdf_fg%nx
       grid_fg%ny                 = netcdf_fg%ny
       grid_fg%nz                 = netcdf_fg%nz
       call variable_interface_setup_struct(grid_fg)
       
       ! Loop through local variable

       do j = 1, grid_fg%nz

          ! Define local variables

          grid_fg%var(:,j) = reshape(netcdf_fg%var(:,:,j,1),              &
               & shape(grid_fg%var(:,j)))
          grid_an%var(:,j) = reshape(netcdf_an%var(:,:,j,1),              &
               & shape(grid_an%var(:,j)))

       end do ! do j = 1, grid_fg%nz

       ! Define local variables

       mgrdn      = .false.
       ugrdn      = .false.
       vgrdn      = .false.
       if(nmm_json_in(i)%gridtype .eq. 1) mgrdn = .true.
       if(nmm_json_in(i)%gridtype .eq. 2) ugrdn = .true.
       if(nmm_json_in(i)%gridtype .eq. 3) vgrdn = .true.
       grid_ic%nz = grid_fg%nz
       call nmm_setup_grid(hwrf_an_filename,nmm_json_grid_in,grid_ic,     &
            & mgrd=mgrdn,ugrd=ugrdn,vgrd=vgrdn)

       ! Compute local variables

       call analysis_increment(grid_fg,grid_an,grid_ic)
       if(is_wnblnd)    call blend_interface_wnd(grid_ic)
       grid_an%var = grid_fg%var + grid_ic%var

       ! Check local variable and proceed accordingly

       if(is_blend_env) then

          ! Define local variables
          
          grid_an%lat = grid_ic%lat
          grid_an%lon = grid_ic%lon

          ! Compute local variables

          call blend_interface_tc_environment(grid_fg,grid_an)

       end if ! if(is_blend_env)

       ! Loop through local variable

       do j = 1, grid_an%nz

          ! Define local variables

          netcdf_an%var(:,:,j,1) = reshape(grid_an%var(:,j),              &
               & shape(netcdf_an%var(:,:,j,1)))

       end do ! do j = 1, grid_an%nz

       ! Define local variables

       call fileio_interface_write(hwrf_an_filename,nmm_json_out,         &
            & netcdf_an)
       
       ! Deallocate memory for local variables

       call variable_interface_cleanup_struct(netcdf_fg)
       call variable_interface_cleanup_struct(netcdf_an)
       call variable_interface_cleanup_struct(grid_an)
       call variable_interface_cleanup_struct(grid_fg)
       call variable_interface_cleanup_struct(grid_ic)
       
    end do ! do i = 1, nvars

    ! Deallocate memory for local variables

    if(allocated(nmm_json_in)) deallocate(nmm_json_in)

    !=====================================================================

  end subroutine hwrf_da_analysis

  !=======================================================================

  ! SUBROUTINE:

  ! hwrf_ghost_analysis.f90

  ! DESCRIPTION:

  ! This subroutine replaces Hurricane WRF (HWRF) forecast input file
  ! variables using those from HWRF data-assimilation (e.g.,
  ! ghost-domain) variables.

  !-----------------------------------------------------------------------

  subroutine hwrf_ghost_analysis()

    ! Define variables computed within routine

    type(nmm_json_var_struct),  dimension(:),               allocatable :: nmm_json_in
    type(slint_struct)                                                  :: slint(3)
    type(nmm_json_grid_struct)                                          :: nmm_json_grid_in
    type(nmm_json_grid_struct)                                          :: nmm_json_grid_out
    type(nmm_json_var_struct)                                           :: nmm_json_out
    type(netcdf_var_struct)                                             :: netcdf_ghst
    type(netcdf_var_struct)                                             :: netcdf_inp
    real(r_kind),               dimension(:),               allocatable :: workgrid
    integer                                                             :: gridtype
    integer                                                             :: nvars

    ! Define counting variables

    integer                                                             :: i, j

    !=====================================================================

    ! Define local variables

    nmm_in_filename  = hwrf_fg_filename
    nmm_out_filename = hwrf_an_filename
    call fileio_interface_read(hwrf_grid_json_vtable,nmm_json_grid_in,    &
         & nmm_json_grid_out)    
    call fileio_interface_read(hwrf_var_json_vtable,nvars,nmm_json_in)
    call nmm_setup_interpolation(nmm_json_grid_in,nmm_json_grid_out,      &
         & slint)

    ! Loop through local variable

    do i = 1, nvars

       ! Define local variables

       gridtype                    = nmm_json_in(i)%gridtype
       call fileio_interface_read(hwrf_fg_filename,nmm_json_in(i),        &
            & netcdf_ghst)
       nmm_json_out                = nmm_json_in(i)
       nmm_json_out%variable_name  = nmm_json_in(i)%remap_variable_name
       call fileio_interface_read(hwrf_an_filename,nmm_json_out,          &
            & netcdf_inp)

       ! Check local variable and proceed accordingly

       if(is_slint) then

          ! Allocate memory for local variables

          if(.not. allocated(workgrid))                                   &
               & allocate(workgrid(slint(gridtype)%src_ncoords))

       end if ! if(is_slint)

       ! Loop through local variable

       do j = 1, netcdf_ghst%nz

          ! Define local variables

          workgrid = reshape(netcdf_ghst%var(:,:,j,1),shape(workgrid))

          ! Check local variable and proceed accordingly

          if(is_slint) then

             ! Compute local variables

             call nmm_interpolation(nmm_json_in(i),slint(gridtype),       &
                  & workgrid)

             ! Define local variables

             netcdf_inp%var(:,:,j,1) = reshape(slint(gridtype)%var,       &
                  & (/netcdf_inp%nx,netcdf_inp%ny/))

          end if ! if(is_slint)
          
       end do ! do j = 1, netcdf_ghst%nz

       ! Define local variables

       call fileio_interface_write(hwrf_an_filename,nmm_json_out,         &
            & netcdf_inp)

       ! Deallocate memory for local variables

       if(allocated(workgrid)) deallocate(workgrid)
       call variable_interface_cleanup_struct(netcdf_ghst)
       call variable_interface_cleanup_struct(netcdf_inp)

    end do ! do i = 1, nvars

    ! Deallocate memory for local variables

    if(allocated(nmm_json_in)) deallocate(nmm_json_in)
    call interpolation_interface_cleanup(slint)

    !=====================================================================

  end subroutine hwrf_ghost_analysis

  !=======================================================================

  ! SUBROUTINE:

  ! hwrf_nest_analysis.f90

  ! DESCRIPTION:

  ! This is the driver routine to blend HWRF analysis variables,
  ! within a parent domain nest, with the same analysis variables
  ! within the respective parent domain.

  !-----------------------------------------------------------------------

  subroutine hwrf_nest_analysis()

    ! Define variables computed within routine

    type(nmm_json_var_struct),  dimension(:),               allocatable :: nmm_json_in
    type(tcv_struct),           dimension(:),               allocatable :: tcv
    type(hwrf_sfc_grid_struct)                                          :: hwrf_sfc_grid
    type(slint_struct)                                                  :: slint(3)
    type(blend_struct)                                                  :: blend
    type(grid_struct)                                                   :: grid_nest
    type(grid_struct)                                                   :: grid_parent
    type(hwrf_vert_grid_struct)                                         :: hwrf_vert_grid
    type(netcdf_var_struct)                                             :: netcdf_nest
    type(netcdf_var_struct)                                             :: netcdf_parent
    type(nmm_json_grid_struct)                                          :: nmm_json_grid_in
    type(nmm_json_grid_struct)                                          :: nmm_json_grid_out
    type(nmm_json_var_struct)                                           :: nmm_json_out
    logical                                                             :: mgrdn
    logical                                                             :: ugrdn
    logical                                                             :: vgrdn
    real(r_kind),               dimension(:,:),             allocatable :: nest_incr
    integer                                                             :: gridtype
    integer                                                             :: intrptype
    integer                                                             :: nvars
  
    ! Define counting variables

    integer                                                             :: i, j

    !=====================================================================
    
    ! Define local variables

    call fileio_interface_read(tcv_filename,tcv)
    nmm_in_filename  = hwrf_nest_filename
    nmm_out_filename = hwrf_parent_filename
    call fileio_interface_read(hwrf_grid_json_vtable,nmm_json_grid_in,    &
         & nmm_json_grid_out)    
    call nmm_setup_nest_interpolation(nmm_json_grid_in,nmm_json_grid_out, &
         & slint,tcv(1))
    call fileio_interface_read(hwrf_var_json_vtable,nvars,nmm_json_in)
    if(is_slint) hwrf_vert_grid%ncoords = slint(1)%dst_ncoords

    ! Loop through local variable

    do i = 1, nvars

       ! Check local variable and proceed accordingly

       if(trim(adjustl(nmm_json_in(i)%variable_name)) .eq. 'PINT') then

          ! Define local variables

          gridtype                   = nmm_json_in(i)%gridtype
          intrptype                  = nmm_json_in(i)%intrptype
          call fileio_interface_read(hwrf_nest_filename,nmm_json_in(i),   &
               & netcdf_nest)
          nmm_json_out               = nmm_json_in(i)
          nmm_json_out%variable_name =                                    &
               & nmm_json_in(i)%remap_variable_name
          call fileio_interface_read(hwrf_parent_filename,nmm_json_out,   &
               & netcdf_parent)
          mgrdn                      = .false.
          ugrdn                      = .false.
          vgrdn                      = .false.
          if(gridtype .eq. 1) mgrdn = .true.
          if(gridtype .eq. 2) ugrdn = .true.
          if(gridtype .eq. 3) vgrdn = .true.
          grid_nest%clat             = tcv(1)%lat
          grid_nest%clon             = tcv(1)%lon
          grid_nest%nz               = netcdf_nest%nz
          grid_parent%clat           = tcv(1)%lat
          grid_parent%clon           = tcv(1)%lon
          grid_parent%nz             = netcdf_parent%nz
          call nmm_setup_grid(hwrf_nest_filename,nmm_json_grid_in,        &
               & grid_nest,mgrd=mgrdn,ugrd=ugrdn,vgrd=vgrdn)
          call nmm_setup_grid(hwrf_parent_filename,nmm_json_grid_out,     &
               & grid_parent,mgrd=mgrdn,ugrd=ugrdn,vgrd=vgrdn)
          hwrf_vert_grid%nx          = grid_parent%nx
          hwrf_vert_grid%ny          = grid_parent%ny
          hwrf_vert_grid%nz          = (grid_parent%nz - 1)
          call variable_interface_setup_struct(hwrf_vert_grid)

          ! Loop through local variable

          do j = 1, grid_nest%nz

             ! Define local variables

             grid_nest%var(:,j)   = reshape(netcdf_nest%var(:,:,j,1),     &
                  & shape(grid_nest%var(:,j)))
             grid_parent%var(:,j) = reshape(netcdf_parent%var(:,:,j,1),   &
                  & shape(grid_parent%var(:,j)))
          
             ! Check local variable and proceed accordingly

             if(is_slint) then

                ! Compute local variables

                call nmm_nest_interpolation(nmm_json_in(i),               &
                     & slint(gridtype),grid_nest%var(:,j),                &
                     & grid_parent%var(:,j))

             end if ! if(is_slint)

             ! Define local variables

             where(grid_parent%var(:,j) .eq. spval) grid_parent%var(:,j)  &
                  & = reshape(netcdf_parent%var(:,:,j,1),                 &
                  & shape(grid_parent%var(:,j)))

          end do ! do j = 1, grid_nest%nz

          ! Define local variables

          hwrf_vert_grid%nest_pint = grid_parent%var

          ! Deallocate memory for local variables

          call variable_interface_cleanup_struct(netcdf_nest)
          call variable_interface_cleanup_struct(netcdf_parent)
          call variable_interface_cleanup_struct(grid_nest)
          call variable_interface_cleanup_struct(grid_parent)

       end if ! if(trim(adjustl(nmm_json_in(i)%variable_name))
              ! .eq. 'PINT')

    end do ! do i = 1, nvars

    ! Define local variables

    hwrf_sfc_grid%nx = grid_parent%nx
    hwrf_sfc_grid%ny = grid_parent%ny
    call nmm_surface(nmm_out_filename,hwrf_sfc_grid)

    ! Loop through local variable

    do i = 1, nvars

       ! Check local variable and proceed accordingly

       if(trim(adjustl(nmm_json_in(i)%variable_name)) .eq. 'PD') then

          ! Define local variables

          gridtype                   = nmm_json_in(i)%gridtype
          intrptype                  = nmm_json_in(i)%intrptype
          call fileio_interface_read(hwrf_nest_filename,nmm_json_in(i),   &
               & netcdf_nest)
          nmm_json_out               = nmm_json_in(i)
          nmm_json_out%variable_name =                                    &
               & nmm_json_in(i)%remap_variable_name
          call fileio_interface_read(hwrf_parent_filename,nmm_json_out,   &
               & netcdf_parent)
          mgrdn                      = .false.
          ugrdn                      = .false.
          vgrdn                      = .false.
          if(gridtype .eq. 1) mgrdn = .true.
          if(gridtype .eq. 2) ugrdn = .true.
          if(gridtype .eq. 3) vgrdn = .true.
          grid_nest%clat             = tcv(1)%lat
          grid_nest%clon             = tcv(1)%lon
          grid_nest%nz               = netcdf_nest%nz
          grid_parent%nz             = netcdf_parent%nz
          call nmm_setup_grid(hwrf_nest_filename,nmm_json_grid_in,        &
               & grid_nest,mgrd=mgrdn,ugrd=ugrdn,vgrd=vgrdn)
          call nmm_setup_grid(hwrf_parent_filename,nmm_json_grid_out,     &
               & grid_parent,mgrd=mgrdn,ugrd=ugrdn,vgrd=vgrdn)
          blend%nx                   = netcdf_parent%nx
          blend%ny                   = netcdf_parent%ny
          blend%nz                   = netcdf_parent%nz
          call variable_interface_setup_struct(blend)
          blend%npts_blend           = npts_blend
          blend%npts_relax           = npts_relax
          blend%lat                  = grid_parent%lat
          blend%lon                  = grid_parent%lon  
          blend%clat                 = tcv(1)%lat
          blend%clon                 = tcv(1)%lon

          ! Loop through local variable

          do j = 1, grid_nest%nz

             ! Define local variables

             grid_nest%var(:,j)   = reshape(netcdf_nest%var(:,:,j,1),     &
                  & shape(grid_nest%var(:,j)))
             grid_parent%var(:,j) = reshape(netcdf_parent%var(:,:,j,1),   &
                  & shape(grid_parent%var(:,j)))
             blend%var_fixed(:,j) = grid_parent%var(:,j)
          
             ! Check local variable and proceed accordingly

             if(is_slint) then

                ! Compute local variables

                call nmm_nest_interpolation(nmm_json_in(i),               &
                     & slint(gridtype),grid_nest%var(:,j),                &
                     & grid_parent%var(:,j))

             end if ! if(is_slint)

             ! Define local variables

             where(grid_parent%var(:,j) .eq. spval) grid_parent%var(:,j)  &
                  & = reshape(netcdf_parent%var(:,:,j,1),                 &
                  & shape(grid_parent%var(:,j)))
             blend%var_blend(:,j) = grid_parent%var(:,j)

          end do ! do j = 1, grid_nest%nz

          ! Compute local variables

          call blend_interface_radii(grid_nest,blend)
          call blend_interface_nest(blend)

          ! Define local variables

          hwrf_vert_grid%pd = blend%var_blend(:,1)
         
          ! Deallocate memory for local variables

          call variable_interface_cleanup_struct(netcdf_nest)
          call variable_interface_cleanup_struct(netcdf_parent)
          call variable_interface_cleanup_struct(grid_nest)
          call variable_interface_cleanup_struct(grid_parent)
          call variable_interface_cleanup_struct(blend)

       end if ! if(trim(adjustl(nmm_json_in(i)%variable_name))
              ! .eq. 'PD')

    end do ! do i = 1, nvars

    ! Compute local variables

    call nmm_vertical_levels(hwrf_parent_filename,hwrf_vert_grid)

    ! Define local variables

    call fileio_interface_write(hwrf_parent_filename,hwrf_vert_grid)

    ! Loop through local variable

    do i = 1, nvars
    
       ! Check local variable and proceed accordingly
       
       if((trim(adjustl(nmm_json_in(i)%variable_name)) .ne. 'PD') .and.   &
            & (trim(adjustl(nmm_json_in(i)%variable_name)) .ne. 'PINT'))  &
            & then

          ! Define local variables
          
          gridtype                   = nmm_json_in(i)%gridtype
          intrptype                  = nmm_json_in(i)%intrptype
          call fileio_interface_read(hwrf_nest_filename,nmm_json_in(i),   &
               & netcdf_nest)
          nmm_json_out               = nmm_json_in(i)
          nmm_json_out%variable_name = nmm_json_in(i)%remap_variable_name
          call fileio_interface_read(hwrf_parent_filename,                &
               & nmm_json_out,netcdf_parent)
          mgrdn                      = .false.
          ugrdn                      = .false.
          vgrdn                      = .false.
          if(gridtype .eq. 1) mgrdn = .true.
          if(gridtype .eq. 2) ugrdn = .true.
          if(gridtype .eq. 3) vgrdn = .true.
          grid_nest%clat             = tcv(1)%lat
          grid_nest%clon             = tcv(1)%lon
          grid_nest%nz               = netcdf_nest%nz
          grid_parent%nz             = netcdf_parent%nz
          call nmm_setup_grid(hwrf_nest_filename,nmm_json_grid_in,        &
               & grid_nest,mgrd=mgrdn,ugrd=ugrdn,vgrd=vgrdn)
          call nmm_setup_grid(hwrf_parent_filename,nmm_json_grid_out,     &
               & grid_parent,mgrd=mgrdn,ugrd=ugrdn,vgrd=vgrdn)
          blend%nx                   = netcdf_parent%nx
          blend%ny                   = netcdf_parent%ny
          blend%nz                   = netcdf_parent%nz
          call variable_interface_setup_struct(blend)
          blend%npts_blend           = npts_blend
          blend%npts_relax           = npts_relax
          blend%lat                  = grid_parent%lat
          blend%lon                  = grid_parent%lon
          blend%clat                 = tcv(1)%lat
          blend%clon                 = tcv(1)%lon

          ! Loop through local variable
          
          do j = 1, grid_nest%nz
             
             ! Define local variables
             
             grid_nest%var(:,j)   = reshape(netcdf_nest%var(:,:,j,1),     &
                  & shape(grid_nest%var(:,j)))
             grid_parent%var(:,j) =                                       &
                  & reshape(netcdf_parent%var(:,:,j,1),                   &
                  & shape(grid_parent%var(:,j)))
             blend%var_fixed(:,j) = grid_parent%var(:,j)
          
             ! Check local variable and proceed accordingly
             
             if(is_slint) then
                
                ! Compute local variables

                call nmm_nest_interpolation(nmm_json_in(i),               &
                     & slint(gridtype),grid_nest%var(:,j),                &
                     & grid_parent%var(:,j))

             end if ! if(is_slint)
             
             ! Define local variables

             blend%var_blend(:,j) = grid_parent%var(:,j)

          end do ! do j = 1, grid_nest%nz

          ! Check local variable and proceed accordingly

          if(nmm_json_in(i)%levtype .eq. 3) then

             ! Define local variables

             hwrf_vert_grid%var(:,1:blend%nz) = blend%var_blend

             ! Compute local variables

             call nmm_vertical_interpolation(nmm_json_in(i),              &
                  & hwrf_vert_grid)

             ! Define local variables

             blend%var_blend = hwrf_vert_grid%var(:,1:blend%nz)

          end if ! if(nmm_json_in(i)%levtype .eq. 3)

          ! Compute local variables
          
          call blend_interface_radii(grid_nest,blend)
          call blend_interface_nest(blend)
          
          ! Loop through local variable

          do j = 1, grid_nest%nz

             ! Define local variables

             netcdf_parent%var(:,:,j,1) = reshape(blend%var_blend(:,j),   &
                  & shape(netcdf_parent%var(:,:,j,1)))

          end do ! do j = 1, grid_nest%nz

          ! Define local variables
          
          call fileio_interface_write(hwrf_parent_filename,               &
               & nmm_json_out,netcdf_parent)
          
          ! Deallocate memory for local variables
          
          call variable_interface_cleanup_struct(netcdf_nest)
          call variable_interface_cleanup_struct(netcdf_parent)
          call variable_interface_cleanup_struct(grid_nest)
          call variable_interface_cleanup_struct(grid_parent)
          call variable_interface_cleanup_struct(blend)

       end if ! if((trim(adjustl(nmm_json_in(i)%variable_name))
              ! .ne. 'PD')
              ! .and. (trim(adjustl(nmm_json_in(i)%variable_name))
              ! .ne. 'PINT'))

    end do ! do i = 1, nvars

    ! Deallocate memory for local variables

    if(allocated(nmm_json_in)) deallocate(nmm_json_in)
    if(allocated(tcv))         deallocate(tcv)
    call variable_interface_cleanup_struct(hwrf_vert_grid)
    call variable_interface_cleanup_struct(hwrf_sfc_grid)
    call interpolation_interface_cleanup(slint)

    !=====================================================================

  end subroutine hwrf_nest_analysis

  !=======================================================================

  ! SUBROUTINE:

  ! state_merge.f90

  ! DESCRIPTION:

  ! This subroutine merges (e.g., blends) HWRF state variables defined
  ! on identical grid projections.

  ! INPUT VARIABLES:

  ! * json; a FORTRAN nmm_json_grid_struct variable.

  ! * gridtype; a FORTRAN integer value specifying the grid type;
  !   current options are:

  !   1: a mass (e.g., thermodynamic) variable.

  !   2: a zonal coordinate-defined (i.e., u-wind) variable.

  !   3: a meridional coordinate-defined (i.e., v-wind) variable.

  ! * netcdf_in; a FORTRAN netcdf_var_struct containing the variable
  !   array to be blended with the contents of the netcdf_out variable
  !   array (see below).

  ! * netcdf_out; a FORTRAN netcdf_var_struct variable containing the
  !   variable array within which to blend the netcdf_in variable
  !   array.

  ! OUTPUT VARIABLES:

  ! * netcdf_out; a FORTRAN netcdf_var_struct variable containing the
  !   blended netcdf_in and the (input) netcdf_out variable arrays.

  !-----------------------------------------------------------------------

  subroutine state_merge(json,gridtype,netcdf_in,netcdf_out)

    ! Define variables passed to routine

    type(nmm_json_grid_struct)                                          :: json
    type(netcdf_var_struct)                                             :: netcdf_in
    type(netcdf_var_struct)                                             :: netcdf_out
    integer                                                             :: gridtype

    ! Define variables computed within routine

    type(netcdf_grid_struct)                                            :: netcdf_grid
    type(grid_struct)                                                   :: in_grid
    type(grid_struct)                                                   :: out_grid

    ! Define counting variables

    integer                                                             :: i

    !=====================================================================

    ! Define local variables

    call fileio_interface_read(hwrf_fg_filename,json,netcdf_grid)
    in_grid%nz  = netcdf_out%nz
    out_grid%nz = netcdf_out%nz

    ! Check local variable and proceed accordingly
    
    if(gridtype .eq. 1) then

       ! Define local variables

       in_grid%nx  = netcdf_grid%nx
       in_grid%ny  = netcdf_grid%ny
       out_grid%nx = netcdf_grid%nx
       out_grid%ny = netcdf_grid%ny
       
    end if ! if(gridtype .eq. 1)

    ! Check local variable and proceed accordingly
    
    if(gridtype .eq. 2) then

       ! Define local variables

       in_grid%nx  = netcdf_grid%nxs
       in_grid%ny  = netcdf_grid%ny
       out_grid%nx = netcdf_grid%nxs
       out_grid%ny = netcdf_grid%ny
       
    end if ! if(gridtype .eq. 2)

    ! Check local variable and proceed accordingly
    
    if(gridtype .eq. 3) then

       ! Define local variables

       in_grid%nx  = netcdf_grid%nx
       in_grid%ny  = netcdf_grid%nys
       out_grid%nx = netcdf_grid%nx
       out_grid%ny = netcdf_grid%nys
       
    end if ! if(gridtype .eq. 3)

    ! Define local variables

    call variable_interface_setup_struct(in_grid)
    call variable_interface_setup_struct(out_grid)

    ! Check local variable and proceed accordingly

    if(gridtype .eq. 1) then

       ! Define local variables

       in_grid%lat = reshape(netcdf_grid%m_lat,shape(in_grid%lat))
       if(json%m_lat_radians) in_grid%lat = in_grid%lat*rad2deg
       in_grid%lon = reshape(netcdf_grid%m_lon,shape(in_grid%lon))
       if(json%m_lon_radians) in_grid%lon = in_grid%lon*rad2deg

    end if ! if(gridtype .eq. 1)

    ! Check local variable and proceed accordingly

    if(gridtype .eq. 2) then

       ! Define local variables

       in_grid%lat = reshape(netcdf_grid%u_lat,shape(in_grid%lat))
       if(json%u_lat_radians) in_grid%lat = in_grid%lat*rad2deg
       in_grid%lon = reshape(netcdf_grid%u_lon,shape(in_grid%lon))
       if(json%u_lon_radians) in_grid%lon = in_grid%lon*rad2deg

    end if ! if(gridtype .eq. 2)

    ! Check local variable and proceed accordingly

    if(gridtype .eq. 3) then

       ! Define local variables

       in_grid%lat = reshape(netcdf_grid%v_lat,shape(in_grid%lat))
       if(json%v_lat_radians) in_grid%lat = in_grid%lat*rad2deg
       in_grid%lon = reshape(netcdf_grid%v_lon,shape(in_grid%lon))
       if(json%v_lon_radians) in_grid%lon = in_grid%lon*rad2deg

    end if ! if(gridtype .eq. 3)

    ! Define local variables

    out_grid%lat = in_grid%lat
    out_grid%lon = in_grid%lon

    ! Loop through local variable

    do i = 1, netcdf_out%nz

       ! Define local variables

       in_grid%var(:,i)  = reshape(netcdf_in%var(:,:,i,1),                 &
            & shape(in_grid%var(:,i)))
       out_grid%var(:,i) = reshape(netcdf_out%var(:,:,i,1),                &
            & shape(out_grid%var(:,i)))

    end do ! do i = 1, netcdf_out%nz

    ! Compute local variables

    call blend_interface_boundary(in_grid,out_grid)

    ! Loop through local variable

    do i = 1, netcdf_out%nz

       ! Define local variables
       
       netcdf_out%var(:,:,i,1) = reshape(out_grid%var(:,i),               &
            & shape(netcdf_out%var(:,:,i,1)))

    end do ! do i = 1, netcdf_out%nz

    ! Deallocate memory for local variables

    call variable_interface_cleanup_struct(in_grid)
    call variable_interface_cleanup_struct(out_grid)

    !=====================================================================

  end subroutine state_merge

  !=======================================================================

end module hwrf_analysis_interface
