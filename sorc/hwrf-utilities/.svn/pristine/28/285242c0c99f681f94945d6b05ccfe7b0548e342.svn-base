module nmm_interpolation_interface

  !=======================================================================

  !$$$ MODULE DOCUMENTATION BLOCK
  
  ! da-update-analysis :: nmm_interpolation_interface
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

  use constants_interface
  use fileio_interface
  use grid_methods_interface
  use interpolation_interface
  use math_methods_interface
  use namelist_interface
  use netcdf_interface
  use variable_interface

  ! Define interfaces and attributes for module routines
  
  implicit none
  private
  public :: nmm_in_filename
  public :: nmm_interpolation
  public :: nmm_nest_interpolation
  public :: nmm_out_filename
  public :: nmm_setup_grid
  public :: nmm_setup_interpolation
  public :: nmm_setup_nest_interpolation
  public :: nmm_surface
  public :: nmm_vertical_interpolation
  public :: nmm_vertical_levels
  interface nmm_interpolation
     module procedure slint_interpolation
  end interface nmm_interpolation
  interface nmm_nest_interpolation
     module procedure slint_nest_interpolation
  end interface nmm_nest_interpolation
  interface nmm_surface
     module procedure hwrf_surface
  end interface nmm_surface
  interface nmm_vertical_interpolation
     module procedure hwrf_vertical_interpolation
  end interface nmm_vertical_interpolation
  interface nmm_vertical_levels
     module procedure hwrf_vertical_levels
  end interface nmm_vertical_levels

  ! Define local variables

  character(len=500)                                                    :: nmm_in_filename
  character(len=500)                                                    :: nmm_out_filename

  !-----------------------------------------------------------------------

contains

  !=======================================================================

  ! SUBROUTINE:

  ! hwrf_surface.f90

  ! DESCRIPTION:
  
  ! This subroutine defines the Hurricane WRF (HWRF) surface (e.g.,
  ! topography and sea-surface temperature) grid using the attributes
  ! within the user specified file.

  ! INPUT VARIABLES:

  ! * filename; a FORTRAN character string specifying the HWRF file
  !   name containing the attributes necessary to compute the vertical
  !   coordinate grid.

  ! * sfc_grid; a FORTRAN hwrf_sfc_grid_struct variable; it is assumed
  !   that all arrays have been allocated (e.g., ncoords is defined).

  ! OUTPUT VARIABLES:

  ! * sfc_grid; a FORTRAN hwrf_vert_sfc_struct variable containing the
  !   topography (fis) and the sea-surface temperature (sst). 

  !-----------------------------------------------------------------------

  subroutine hwrf_surface(filename,sfc_grid)

    ! Define variables passed to routine

    type(hwrf_sfc_grid_struct)                                          :: sfc_grid
    character(len=500)                                                  :: filename

    ! Define variables computed within routine

    type(netcdf_var_struct)                                             :: netcdf_var
    character(len=100)                                                  :: varname

    !=====================================================================

    ! Define local variables

    netcdf_var%nx    = sfc_grid%nx
    netcdf_var%ny    = sfc_grid%ny
    netcdf_var%nz    = 1
    netcdf_var%nt    = 1
    call variable_interface_setup_struct(netcdf_var)
    sfc_grid%ncoords = (netcdf_var%nx*netcdf_var%ny)
    call variable_interface_setup_struct(sfc_grid)
    varname          = 'FIS'
    call netcdf_interface_getvar(filename,varname,netcdf_var%var)
    sfc_grid%fis     = reshape(netcdf_var%var,shape(sfc_grid%fis))
    varname          = 'SST'
    call netcdf_interface_getvar(filename,varname,netcdf_var%var)
    sfc_grid%sst     = reshape(netcdf_var%var,shape(sfc_grid%sst))
    varname          = 'SEAICE'
    call netcdf_interface_getvar(filename,varname,netcdf_var%var)
    sfc_grid%seaice  = reshape(netcdf_var%var,shape(sfc_grid%seaice))

    ! Deallocate memory for local variables

    call variable_interface_cleanup_struct(netcdf_var)

    !=====================================================================

  end subroutine hwrf_surface

  !=======================================================================

  ! SUBROUTINE:

  ! hwrf_vertical_interpolation.f90

  ! DESCRIPTION:

  ! This subroutine interpolates a variable, defined on the Hurricane
  ! WRF (HWRF) vertical grid from a source grid profile to a
  ! destination grid profile in accordance with the interpolation
  ! specifications of the user.

  ! INPUT VARIABLES:

  ! * json; a FORTRAN nmm_json_grid_struct variable containing the
  !   respective grid attributes.

  ! * hwrf_vert_grid; a FORTRAN hwrf_vert_grid_struct variable
  !   containing the HWRF vertical grid attributes for both the source
  !   and destination grids as well as the source variable to be
  !   interpolated.

  ! OUTPUT VARIABLES:

  ! * hwrf_vert_grid; a FORTRAN hwrf_vert_grid_struct variable
  !   containing the interpolated HWRF source variable.

  !-----------------------------------------------------------------------

  subroutine hwrf_vertical_interpolation(json,vert_grid)

    ! Define variables passed to routine

    type(hwrf_vert_grid_struct)                                         :: vert_grid
    type(nmm_json_var_struct)                                           :: json

    !=====================================================================

    ! Compute local variables

    call interpolation_interface_vertical(json,vert_grid)   

    !=====================================================================

  end subroutine hwrf_vertical_interpolation

  !=======================================================================

  ! SUBROUTINE:

  ! hwrf_vertical_levels.f90

  ! DESCRIPTION:

  ! This subroutine defines the Hurricane WRF (HWRF) vertical
  ! coordinate grid using the attributes within the user specified
  ! file.

  ! INPUT VARIABLES:

  ! * filename; a FORTRAN character string specifying the HWRF file
  !   name containing the attributes necessary to compute the vertical
  !   coordinate grid.

  ! * vert_grid; a FORTRAN hwrf_vert_grid_struct variable; it is
  !   assumed that all arrays have been allocated (e.g., ncoords and
  !   nz attributes are defined), the dry-hydrostatic surface pressure
  !   (e.g., pd attribute), and nest-interpolated interface levels
  !   (e.g., nest_pint, if remapping from a nested region) have been
  !   defined prior to calling this routine.

  ! OUTPUT VARIABLES:

  ! * vert_grid; a FORTRAN hwrf_vert_grid_struct variable containing
  !   the pressure levels (plev), the interface pressure levels
  !   (pint), the vertical grid attributes (eta1, eta2, aeta1, aeta2,
  !   pdtop, and pt), the nest-interpolated pressure levels (e.g.,
  !   nest_plev), and the surface pressure (psfc).

  !-----------------------------------------------------------------------

  subroutine hwrf_vertical_levels(filename,vert_grid)

    ! Define variables passed to routine

    type(hwrf_vert_grid_struct)                                         :: vert_grid
    character(len=500)                                                  :: filename

    ! Define variables computed within routine

    type(netcdf_var_struct)                                             :: netcdf_var
    character(len=100)                                                  :: varname

    ! Define counting variables

    integer                                                             :: i

    !=====================================================================

    ! Loop through local variable

    do i = 1, vert_grid%nz

       ! Compute local variables

       vert_grid%nest_plev(:,i) = 0.5*(vert_grid%nest_pint(:,i) +          &
            & vert_grid%nest_pint(:,(i + 1)))

    end do ! do i = 1, vert_grid%nz

    ! Define local variables

    netcdf_var%nx   = vert_grid%nz
    netcdf_var%ny   = 1
    netcdf_var%nz   = 1    
    netcdf_var%nt   = 1
    call variable_interface_setup_struct(netcdf_var)
    varname         = 'AETA1'
    call netcdf_interface_getvar(filename,varname,netcdf_var%var)
    vert_grid%aeta1 = netcdf_var%var(:,1,1,1)
    varname         = 'AETA2'
    call netcdf_interface_getvar(filename,varname,netcdf_var%var)
    vert_grid%aeta2 = netcdf_var%var(:,1,1,1)

    ! Deallocate memory for local variables

    call variable_interface_cleanup_struct(netcdf_var)

    ! Define local variables

    netcdf_var%nx  = (vert_grid%nz + 1)
    netcdf_var%ny  = 1
    netcdf_var%nz  = 1    
    netcdf_var%nt  = 1
    call variable_interface_setup_struct(netcdf_var)
    varname        = 'ETA1'
    call netcdf_interface_getvar(filename,varname,netcdf_var%var)
    vert_grid%eta1 = netcdf_var%var(:,1,1,1)
    varname        = 'ETA2'
    call netcdf_interface_getvar(filename,varname,netcdf_var%var)
    vert_grid%eta2 = netcdf_var%var(:,1,1,1)

    ! Deallocate memory for local variables

    call variable_interface_cleanup_struct(netcdf_var)

    ! Define local variables

    netcdf_var%nx   = 1
    netcdf_var%ny   = 1
    netcdf_var%nz   = 1    
    netcdf_var%nt   = 1
    call variable_interface_setup_struct(netcdf_var)
    varname         = 'PT'
    call netcdf_interface_getvar(filename,varname,netcdf_var%var)
    vert_grid%pt    = netcdf_var%var(1,1,1,1)
    varname         = 'PDTOP'
    call netcdf_interface_getvar(filename,varname,netcdf_var%var)
    vert_grid%pdtop = netcdf_var%var(1,1,1,1)

    ! Deallocate memory for local variables

    call variable_interface_cleanup_struct(netcdf_var)

    ! Compute local variables

    vert_grid%psfc = (vert_grid%pd + vert_grid%pdtop + vert_grid%pt)

    ! Loop through local variable

    do i = 1, vert_grid%nz

       ! Compute local variables

       vert_grid%plev(:,i) = (vert_grid%aeta1(i)*vert_grid%pdtop) +        &
            & vert_grid%aeta2(i)*(vert_grid%psfc(:) - vert_grid%pdtop -    &
            & vert_grid%pt) + vert_grid%pt

    end do ! do i = 1, vert_grid%nz

    ! Loop through local variable

    do i = 1, (vert_grid%nz + 1)

       ! Compute local variables

       vert_grid%pint(:,i) = (vert_grid%eta1(i)*vert_grid%pdtop) +         &
            & vert_grid%eta2(i)*(vert_grid%psfc(:) - vert_grid%pdtop -     &
            & vert_grid%pt) + vert_grid%pt

    end do ! do i = 1, (vert_grid%nz + 1)

    !=====================================================================

  end subroutine hwrf_vertical_levels

  !=======================================================================

  ! SUBROUTINE:

  ! nmm_setup_grid.f90

  ! DESCRIPTION:

  ! This subroutine defines the grid projection in accordance with the
  ! user grid-type (e.g., mass or momentum component).

  ! INPUT VARIABLES:

  ! * nmm_filename; a FORTRAN character string specifying the
  !   full-path to the NMM netcdf file containing the grid projection
  !   variables.

  ! * json; a FORTRAN nmm_json_grid_struct variable containing the
  !   respective grid attributes.

  ! * grid; a FORTRAN grid_struct variable.

  ! OPTIONAL INPUT VARIABLES:

  ! * mgrd; a FORTRAN logical variable; a value of .true. specifies
  !   that the user wishes to define the geographical locations within
  !   the grid_struct variable (e.g., grid) along the mass variable
  !   geographical coordinate values.

  ! * ugrd; a FORTRAN logical variable; a value of .true. specifies
  !   that the user wishes to define the geographical locations within
  !   the grid_struct variable (e.g., grid) along the zonal-wind
  !   variable geographical coordinate values.

  ! * vgrd; a FORTRAN logical variable; a value of .true. specifies
  !   that the user wishes to define the geographical locations within
  !   the grid_struct variable (e.g., grid) along the meridional-wind
  !   variable geographical coordinate values.

  ! OUTPUT VARIABLES:

  ! * grid; a FORTRAN grid_struct variable containing the geographical
  !   coordinate variables in accordance with the user specified
  !   grid-type.

  !-----------------------------------------------------------------------

  subroutine nmm_setup_grid(nmm_filename,json,grid,mgrd,ugrd,vgrd)

    ! Define variables passed to routine

    type(nmm_json_grid_struct)                                          :: json
    type(grid_struct)                                                   :: grid
    character(len=500)                                                  :: nmm_filename
    logical, optional, intent(in)                                       :: mgrd
    logical, optional, intent(in)                                       :: ugrd
    logical, optional, intent(in)                                       :: vgrd

    ! Define variables computed within routine

    type(grid_struct)                                                   :: dst_grid
    type(grid_struct)                                                   :: src_grid
    type(netcdf_grid_struct)                                            :: netcdf
    type(slint_struct)                                                  :: slint
    logical                                                             :: mgrdn
    logical                                                             :: ugrdn
    logical                                                             :: vgrdn

    !=====================================================================

    ! Check local variable and proceed accordingly

    if(.not. present(mgrd)) then

       ! Define local variables

       mgrdn = .false.

    else   ! if(.not. present(mgrd))

       ! Define local variables

       mgrdn = mgrd

    end if ! if(.not. present(mgrd))

    ! Check local variable and proceed accordingly

    if(.not. present(ugrd)) then

       ! Define local variables

       ugrdn = .false.

    else   ! if(.not. present(ugrd))

       ! Define local variables

       ugrdn = ugrd

    end if ! if(.not. present(ugrd))

    ! Check local variable and proceed accordingly

    if(.not. present(vgrd)) then

       ! Define local variables

       vgrdn = .false.

    else   ! if(.not. present(vgrd))

       ! Define local variables

       vgrdn = vgrd

    end if ! if(.not. present(vgrd))

    ! Define local variables

    call fileio_interface_read(nmm_filename,json,netcdf)
    if(json%m_lat_radians) netcdf%m_lat = netcdf%m_lat*rad2deg
    if(json%m_lon_radians) netcdf%m_lon = netcdf%m_lon*rad2deg
    if(json%u_lat_radians) netcdf%u_lat = netcdf%u_lat*rad2deg
    if(json%u_lon_radians) netcdf%u_lon = netcdf%u_lon*rad2deg
    if(json%v_lat_radians) netcdf%v_lat = netcdf%v_lat*rad2deg
    if(json%v_lon_radians) netcdf%v_lon = netcdf%v_lon*rad2deg
    if(grid%nz .eq. 0)     grid%nz = netcdf%nz

    ! Check local variable and proceed accordingly

    if(mgrdn) then

       ! Define local variables

       grid%nx            = netcdf%nx
       grid%ny            = netcdf%ny
       call variable_interface_setup_struct(grid)
       grid%lat           = reshape(netcdf%m_lat,shape(grid%lat))
       grid%lon           = reshape(netcdf%m_lon,shape(grid%lon))

    end if ! if(mgrdn)

    ! Check local variable and proceed accordingly

    if(ugrdn) then

       ! Define local variables

       grid%nx          = netcdf%nxs
       grid%ny          = netcdf%ny
       call variable_interface_setup_struct(grid)
       grid%lat         = reshape(netcdf%u_lat,shape(grid%lat))
       grid%lon         = reshape(netcdf%u_lon,shape(grid%lon))

    end if ! if(ugrdn)

    ! Check local variable and proceed accordingly

    if(vgrdn) then

       ! Define local variables

       grid%nx          = netcdf%nx
       grid%ny          = netcdf%nys
       call variable_interface_setup_struct(grid)
       grid%lat         = reshape(netcdf%v_lat,shape(grid%lat))
       grid%lon         = reshape(netcdf%v_lon,shape(grid%lon))

    end if ! if(vgrdn)

    ! Compute local variables

    call grid_methods_gcdist(grid,netcdf%dx,netcdf%dy)
    call math_methods_polarcoords(grid,grid%clon,grid%clat)

    ! Deallocate memory for local variables

    call variable_interface_cleanup_struct(netcdf) 
    call variable_interface_cleanup_struct(slint)

    !=====================================================================

  end subroutine nmm_setup_grid

  !=======================================================================

  ! SUBROUTINE:

  ! nmm_setup_interpolation.f90

  ! DESCRIPTION:

  ! This subroutine initializes the respective remapping method (e.g.,
  ! SLINT, ESMF, etc.) FORTRAN structures for the various Arakawa-type
  ! grid projections.

  ! INPUT VARIABLES:

  ! * json_in; a FORTRAN nmm_json_grid_struct variable containing the
  !   source (e.g., input) grid attributes.

  ! * json_out; a FORTRAN nmm_json_grid_struct variable containing the
  !   destination (e.g., output) grid attributes.

  ! * slint; a FORTRAN slint_struct variable array of dimension 3.

  ! OUTPUT VARIABLES:

  ! * slint; a FORTRAN slint_struct variable array, of dimension 3,
  !   containing the remapping attributes for the various Arakawa-type
  !   grid projections.

  !-----------------------------------------------------------------------

  subroutine nmm_setup_interpolation(json_in,json_out,slint)

    ! Define variables passed to routine

    type(slint_struct)                                                  :: slint(3)
    type(nmm_json_grid_struct)                                          :: json_in
    type(nmm_json_grid_struct)                                          :: json_out

    ! Define variables computed within routine

    type(grid_struct)                                                   :: dst_grid
    type(grid_struct)                                                   :: src_grid

    !=====================================================================

    ! Define local variables

    call nmm_setup_grid(nmm_in_filename,json_in,src_grid,mgrd=.true.)
    call nmm_setup_grid(nmm_out_filename,json_out,dst_grid,mgrd=.true.)

    ! Check local variable and proceed accordingly

    if(is_slint) then

       ! Define local variables

       slint(1)%ncoords = dst_grid%ncoords
       call variable_interface_setup_struct(slint(1))

       ! Compute local variables

       call interpolation_interface_init(src_grid,dst_grid,slint(1))

    end if ! if(is_slint)

    ! Deallocate memory for local variables

    call variable_interface_cleanup_struct(src_grid)
    call variable_interface_cleanup_struct(dst_grid)

    ! Define local variables

    call nmm_setup_grid(nmm_in_filename,json_in,src_grid,ugrd=.true.)
    call nmm_setup_grid(nmm_out_filename,json_out,dst_grid,ugrd=.true.)

    ! Check local variable and proceed accordingly

    if(is_slint) then

       ! Define local variables

       slint(2)%ncoords = dst_grid%ncoords
       call variable_interface_setup_struct(slint(2))

       ! Compute local variables

       call interpolation_interface_init(src_grid,dst_grid,slint(2))

    end if ! if(is_slint)

    ! Deallocate memory for local variables

    call variable_interface_cleanup_struct(src_grid)
    call variable_interface_cleanup_struct(dst_grid)

    ! Define local variables

    call nmm_setup_grid(nmm_in_filename,json_in,src_grid,vgrd=.true.)
    call nmm_setup_grid(nmm_out_filename,json_out,dst_grid,vgrd=.true.)

    ! Check local variable and proceed accordingly

    if(is_slint) then

       ! Define local variables

       slint(3)%ncoords = dst_grid%ncoords
       call variable_interface_setup_struct(slint(3))

       ! Compute local variables

       call interpolation_interface_init(src_grid,dst_grid,slint(3))

    end if ! if(is_slint)

    ! Deallocate memory for local variables

    call variable_interface_cleanup_struct(src_grid)
    call variable_interface_cleanup_struct(dst_grid)

    !=====================================================================

  end subroutine nmm_setup_interpolation

  !=======================================================================

  ! SUBROUTINE:

  ! nmm_setup_nest_interpolation.f90

  ! DESCRIPTION:

  ! This subroutine initializes the respective remapping method (e.g.,
  ! SLINT, ESMF, etc.) FORTRAN structures for the various Arakawa-type
  ! grid projections to be used to replace co-located nest variable
  ! locations within their respective parent grid(s).

  ! INPUT VARIABLES:

  ! * json_in; a FORTRAN nmm_json_grid_struct variable containing the
  !   source (e.g., input) grid attributes.

  ! * json_out; a FORTRAN nmm_json_grid_struct variable containing the
  !   destination (e.g., output) grid attributes.

  ! * slint; a FORTRAN slint_struct variable array of dimension 3.

  ! OUTPUT VARIABLES:

  ! * slint; a FORTRAN slint_struct variable array, of dimension 3,
  !   containing the remapping attributes for the various Arakawa-type
  !   grid projections.

  !-----------------------------------------------------------------------

  subroutine nmm_setup_nest_interpolation(json_in,json_out,slint,tcv)

    ! Define variables passed to routine

    type(slint_struct)                                                  :: slint(3)
    type(nmm_json_grid_struct)                                          :: json_in
    type(nmm_json_grid_struct)                                          :: json_out
    type(tcv_struct), optional, intent(in)                              :: tcv

    ! Define variables computed within routine

    type(grid_struct)                                                   :: dst_grid
    type(grid_struct)                                                   :: dst_sub_grid
    type(grid_struct)                                                   :: src_grid
    type(kdtree_struct)                                                 :: kdtree
    integer,                    dimension(:),               allocatable :: idx
    real(r_kind)                                                        :: lats(4)
    real(r_kind)                                                        :: lons(4)
    real(r_kind)                                                        :: src_mn_lat
    real(r_kind)                                                        :: src_mn_lon
    real(r_kind)                                                        :: src_mx_lat
    real(r_kind)                                                        :: src_mx_lon
    integer                                                             :: dst_sub_cnt

    ! Define counting variables

    integer                                                             :: i

    !=====================================================================
    
    ! Check local variable and proceed accordingly

    if(present(tcv)) then

       ! Define local variables

       src_grid%clat = tcv%lat
       src_grid%clon = tcv%lon

    end if ! if(present(tcv))
    
    ! Define local variables

    call nmm_setup_grid(nmm_in_filename,json_in,src_grid,mgrd=.true.)
    call nmm_setup_grid(nmm_out_filename,json_out,dst_grid,mgrd=.true.)

    ! Check local variable and proceed accordingly

    if(present(tcv)) then

       ! Define local variables

       kdtree%ncoords = 1
       kdtree%nn      = dst_grid%ncoords
       call variable_interface_setup_struct(kdtree)
       kdtree%r2      = maxval(src_grid%radius)*maxval(src_grid%radius)
       kdtree%nalloc  = dst_grid%ncoords
    
       ! Allocate memory for local variables

       if(.not. allocated(idx)) allocate(idx(dst_grid%ncoords))

       ! Compute local variables

       call math_methods_kdtree_r2(src_grid,dst_grid,kdtree)

       ! Define local variables
       
       dst_sub_grid%ncoords = kdtree%nfound
       
    else  ! if(present(tcv))
    
       ! Define local variables

       src_mn_lat           = minval(src_grid%lat)
       src_mn_lon           = minval(src_grid%lon)
       src_mx_lat           = maxval(src_grid%lat)
       src_mx_lon           = maxval(src_grid%lon)
       dst_sub_grid%ncoords = count((dst_grid%lat .ge. src_mn_lat) .and.   &
            & (dst_grid%lat .le. src_mx_lat) .and. (dst_grid%lon .ge.      &
            & src_mn_lon) .and. (dst_grid%lon .le. src_mx_lon))

    end if ! if(present(tcv))

    ! Define local variables
    
    dst_sub_grid%nz = dst_grid%nz
    call variable_interface_setup_struct(dst_sub_grid)
    dst_sub_cnt     = 0

    ! Check local variable and proceed accordingly

    if(present(tcv)) then

       ! Loop through local variable

       do i = 1, kdtree%nfound

          ! Define local variables

          dst_sub_cnt                   = dst_sub_cnt + 1
          dst_sub_grid%lat(dst_sub_cnt) = dst_grid%lat(kdtree%idx(1,i))
          dst_sub_grid%lon(dst_sub_cnt) = dst_grid%lon(kdtree%idx(1,i))
          dst_sub_grid%idx(dst_sub_cnt) = kdtree%idx(1,i)
          
       end do ! do i = 1, kdtree%nfound

    else  ! if(present(tcv))
       
       ! Loop through local variable

       do i = 1, dst_grid%ncoords

          ! Check local variable and proceed accordingly
 
          if((dst_grid%lat(i) .ge. src_mn_lat) .and. (dst_grid%lat(i)      &
               & .le. src_mx_lat) .and. (dst_grid%lon(i) .ge. src_mn_lon)  &
               & .and. (dst_grid%lon(i) .le. src_mx_lon)) then
    
             ! Define local variables

             dst_sub_cnt                   = dst_sub_cnt + 1
             dst_sub_grid%lat(dst_sub_cnt) = dst_grid%lat(i)
             dst_sub_grid%lon(dst_sub_cnt) = dst_grid%lon(i)
             dst_sub_grid%idx(dst_sub_cnt) = i

          end if ! if((dst_grid%lat(i) .ge. src_mn_lat)
                 ! .and. (dst_grid%lat(i) .le. src_mx_lat)
                 ! .and. (dst_grid%lon(i) .ge. src_mn_lon)
                 ! .and. (dst_grid%lon(i) .le. src_mx_lon))

       end do ! do i = 1, dst_grid%ncoords

    end if ! if(present(tcv))
       
    ! Check local variable and proceed accordingly

    if(is_slint) then

       ! Define local variables

       slint(1)%ncoords     = dst_sub_grid%ncoords
       slint(1)%dst_ncoords = dst_grid%ncoords
       slint(1)%src_ncoords = src_grid%ncoords
       call variable_interface_setup_struct(slint(1))
       slint(1)%idx         = dst_sub_grid%idx
       
       ! Compute local variables

       call interpolation_interface_init(src_grid,dst_sub_grid,slint(1))
       
    end if ! if(is_slint)

    ! Deallocate memory for local variables

    call variable_interface_cleanup_struct(kdtree)
    call variable_interface_cleanup_struct(src_grid)
    call variable_interface_cleanup_struct(dst_grid)
    call variable_interface_cleanup_struct(dst_sub_grid)
    
    ! Define local variables

    call nmm_setup_grid(nmm_in_filename,json_in,src_grid,ugrd=.true.)
    call nmm_setup_grid(nmm_out_filename,json_out,dst_grid,ugrd=.true.)

    ! Check local variable and proceed accordingly

    if(present(tcv)) then

       ! Define local variables

       kdtree%ncoords = 1
       kdtree%nn      = dst_grid%ncoords
       call variable_interface_setup_struct(kdtree)
       kdtree%r2      = maxval(src_grid%radius)*maxval(src_grid%radius)
       kdtree%nalloc  = dst_grid%ncoords
    
       ! Allocate memory for local variables

       if(.not. allocated(idx)) allocate(idx(dst_grid%ncoords))

       ! Compute local variables

       call math_methods_kdtree_r2(src_grid,dst_grid,kdtree)

       ! Define local variables
       
       dst_sub_grid%ncoords = kdtree%nfound
       
    else  ! if(present(tcv))
    
       ! Define local variables

       src_mn_lat = minval(src_grid%lat)
       src_mn_lon = minval(src_grid%lon)
       src_mx_lat = maxval(src_grid%lat)
       src_mx_lon = maxval(src_grid%lon)

    end if ! if(present(tcv))

    ! Define local variables
    
    dst_sub_grid%nz = dst_grid%nz
    call variable_interface_setup_struct(dst_sub_grid)
    dst_sub_cnt     = 0

    ! Check local variable and proceed accordingly

    if(present(tcv)) then

       ! Loop through local variable

       do i = 1, kdtree%nfound

          ! Define local variables

          dst_sub_cnt                   = dst_sub_cnt + 1
          dst_sub_grid%lat(dst_sub_cnt) = dst_grid%lat(kdtree%idx(1,i))
          dst_sub_grid%lon(dst_sub_cnt) = dst_grid%lon(kdtree%idx(1,i))
          dst_sub_grid%idx(dst_sub_cnt) = kdtree%idx(1,i)
          
       end do ! do i = 1, kdtree%nfound

    else  ! if(present(tcv))
       
       ! Loop through local variable

       do i = 1, dst_grid%ncoords

          ! Check local variable and proceed accordingly
 
          if((dst_grid%lat(i) .ge. src_mn_lat) .and. (dst_grid%lat(i)      &
               & .le. src_mx_lat) .and. (dst_grid%lon(i) .ge. src_mn_lon)  &
               & .and. (dst_grid%lon(i) .le. src_mx_lon)) then
    
             ! Define local variables

             dst_sub_cnt                   = dst_sub_cnt + 1
             dst_sub_grid%lat(dst_sub_cnt) = dst_grid%lat(i)
             dst_sub_grid%lon(dst_sub_cnt) = dst_grid%lon(i)
             dst_sub_grid%idx(dst_sub_cnt) = i

          end if ! if((dst_grid%lat(i) .ge. src_mn_lat)
                 ! .and. (dst_grid%lat(i) .le. src_mx_lat)
                 ! .and. (dst_grid%lon(i) .ge. src_mn_lon)
                 ! .and. (dst_grid%lon(i) .le. src_mx_lon))

       end do ! do i = 1, dst_grid%ncoords

    end if ! if(present(tcv))

    ! Check local variable and proceed accordingly

    if(is_slint) then

       ! Define local variables

       slint(2)%ncoords     = dst_sub_grid%ncoords
       slint(2)%dst_ncoords = dst_grid%ncoords
       slint(2)%src_ncoords = src_grid%ncoords
       call variable_interface_setup_struct(slint(2))
       slint(2)%idx         = dst_sub_grid%idx

       ! Compute local variables

       call interpolation_interface_init(src_grid,dst_sub_grid,slint(2))

    end if ! if(is_slint)

    ! Deallocate memory for local variables

    call variable_interface_cleanup_struct(kdtree)
    call variable_interface_cleanup_struct(src_grid)
    call variable_interface_cleanup_struct(dst_grid)
    call variable_interface_cleanup_struct(dst_sub_grid)

    ! Define local variables

    call nmm_setup_grid(nmm_in_filename,json_in,src_grid,vgrd=.true.)
    call nmm_setup_grid(nmm_out_filename,json_out,dst_grid,vgrd=.true.)

    ! Check local variable and proceed accordingly

    if(present(tcv)) then

       ! Define local variables

       kdtree%ncoords = 1
       kdtree%nn      = dst_grid%ncoords
       call variable_interface_setup_struct(kdtree)
       kdtree%r2      = maxval(src_grid%radius)*maxval(src_grid%radius)
       kdtree%nalloc  = dst_grid%ncoords
    
       ! Allocate memory for local variables

       if(.not. allocated(idx)) allocate(idx(dst_grid%ncoords))

       ! Compute local variables

       call math_methods_kdtree_r2(src_grid,dst_grid,kdtree)

       ! Define local variables
       
       dst_sub_grid%ncoords = kdtree%nfound
       
    else  ! if(present(tcv))
    
       ! Define local variables

       src_mn_lat = minval(src_grid%lat)
       src_mn_lon = minval(src_grid%lon)
       src_mx_lat = maxval(src_grid%lat)
       src_mx_lon = maxval(src_grid%lon)

    end if ! if(present(tcv))

    ! Define local variables
    
    dst_sub_grid%nz = dst_grid%nz
    call variable_interface_setup_struct(dst_sub_grid)
    dst_sub_cnt     = 0

    ! Check local variable and proceed accordingly

    if(present(tcv)) then

       ! Loop through local variable

       do i = 1, kdtree%nfound

          ! Define local variables

          dst_sub_cnt                   = dst_sub_cnt + 1
          dst_sub_grid%lat(dst_sub_cnt) = dst_grid%lat(kdtree%idx(1,i))
          dst_sub_grid%lon(dst_sub_cnt) = dst_grid%lon(kdtree%idx(1,i))
          dst_sub_grid%idx(dst_sub_cnt) = kdtree%idx(1,i)
          
       end do ! do i = 1, kdtree%nfound

    else  ! if(present(tcv))
       
       ! Loop through local variable

       do i = 1, dst_grid%ncoords

          ! Check local variable and proceed accordingly
 
          if((dst_grid%lat(i) .ge. src_mn_lat) .and. (dst_grid%lat(i)      &
               & .le. src_mx_lat) .and. (dst_grid%lon(i) .ge. src_mn_lon)  &
               & .and. (dst_grid%lon(i) .le. src_mx_lon)) then
    
             ! Define local variables

             dst_sub_cnt                   = dst_sub_cnt + 1
             dst_sub_grid%lat(dst_sub_cnt) = dst_grid%lat(i)
             dst_sub_grid%lon(dst_sub_cnt) = dst_grid%lon(i)
             dst_sub_grid%idx(dst_sub_cnt) = i

          end if ! if((dst_grid%lat(i) .ge. src_mn_lat)
                 ! .and. (dst_grid%lat(i) .le. src_mx_lat)
                 ! .and. (dst_grid%lon(i) .ge. src_mn_lon)
                 ! .and. (dst_grid%lon(i) .le. src_mx_lon))

       end do ! do i = 1, dst_grid%ncoords

    end if ! if(present(tcv))

    ! Check local variable and proceed accordingly

    if(is_slint) then

       ! Define local variables

       slint(3)%ncoords     = dst_sub_grid%ncoords
       slint(3)%dst_ncoords = dst_grid%ncoords
       slint(3)%src_ncoords = src_grid%ncoords
       call variable_interface_setup_struct(slint(3))
       slint(3)%idx         = dst_sub_grid%idx

       ! Compute local variables

       call interpolation_interface_init(src_grid,dst_sub_grid,slint(3))
       
    end if ! if(is_slint)

    ! Deallocate memory for local variables

    call variable_interface_cleanup_struct(kdtree)
    call variable_interface_cleanup_struct(src_grid)
    call variable_interface_cleanup_struct(dst_grid)
    call variable_interface_cleanup_struct(dst_sub_grid)
    
    !=====================================================================
    
  end subroutine nmm_setup_nest_interpolation

  !=======================================================================
  
  ! SUBROUTINE:

  ! slint_interpolation.f90

  ! DESCRIPTION:

  ! This subroutine remaps a user specified variable using the SLINT
  ! remapping bi-linear and nearest-neighbor methods.

  ! INPUT VARIABLES:

  ! * json; a FORTRAN nmm_json_var_struct containing the respective
  !   variable, to be remapped, attributes.

  ! * slint; a FORTRAN slint_struct variable containing the remapping
  !   attributes for the respective variable type.

  ! * var; a FORTRAN 4-byte real-valued variable array of dimension
  !   'src_ncoords' which contains the source (e.g., input) grid
  !   variable to be remapped.

  ! OUTPUT VARIABLES:

  ! * slint; a FORTRAN slint_struct variable containing the remapped
  !   variable within the struture 'var' array.

  !-----------------------------------------------------------------------

  subroutine slint_interpolation(json,slint,var)

    ! Define variables passed to routine

    type(nmm_json_var_struct)                                           :: json
    type(slint_struct)                                                  :: slint
    real(r_kind)                                                        :: var(slint%src_ncoords)

    !=====================================================================

    ! Check local variable and proceed accordingly

    if(json%intrptype .eq. 1) then

       ! Compute local variables

       call interpolation_interface_nnghbr(slint,var)

    end if ! if(json%intrptype .eq. 1)

    ! Check local variable and proceed accordingly

    if(json%intrptype .eq. 2) then

       ! Compute local variables

       call interpolation_interface_bilinear(slint,var)

    end if ! if(json%intrptype .eq. 2)

    !=====================================================================

  end subroutine slint_interpolation

  !=======================================================================

  ! SUBROUTINE:

  ! slint_nest_interpolation.f90

  ! DESCRIPTION:

  ! This subroutine remaps a user specified variable using the SLINT
  ! remapping bi-linear and nearest-neighbor methods.

  ! INPUT VARIABLES:

  ! * json; a FORTRAN nmm_json_var_struct containing the respective
  !   variable, to be remapped, attributes.

  ! * slint; a FORTRAN slint_struct variable containing the remapping
  !   attributes for the respective variable type.

  ! * var_nest; a FORTRAN 4-byte real-valued variable array of
  !   dimension 'src_ncoords' which contains the source (e.g., input
  !   and/or nest) grid variable to be remapped.

  ! * var_parent; a FORTRAN 4-byte real-valued variable array of
  !   dimension 'dst_ncoords' which contains the destination (e.g.,
  !   output and/or parent) grid variable within which to remap the
  !   nest variable.

  ! OUTPUT VARIABLES:

  ! * slint; a FORTRAN slint_struct variable containing the remapped
  !   variable within the struture 'var' array.

  !-----------------------------------------------------------------------

  subroutine slint_nest_interpolation(json,slint,var_nest,var_parent)

    ! Define variables passed to routine

    type(nmm_json_var_struct)                                           :: json
    type(slint_struct)                                                  :: slint
    real(r_kind)                                                        :: var_nest(slint%src_ncoords)
    real(r_kind)                                                        :: var_parent(slint%dst_ncoords)

    ! Define variables computed within routine

    integer                                                             :: dst_idx

    ! Define counting variables

    integer                                                             :: i

    !=====================================================================

    ! Define local variables

    var_parent = spval

    ! Compute local variables

    call slint_interpolation(json,slint,var_nest)

    ! Loop through local variable

    do i = 1, slint%ncoords
       
       ! Define local variables

       dst_idx             = slint%idx(i)
       var_parent(dst_idx) = slint%var(i)

    end do ! do i = 1, slint%ncoords

    !=====================================================================

  end subroutine slint_nest_interpolation

  !=======================================================================

end module nmm_interpolation_interface
