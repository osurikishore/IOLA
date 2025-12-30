module json_interface

  !=======================================================================

  !$$$ MODULE DOCUMENTATION BLOCK
  
  ! da-update-analysis :: json_interface
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

  use fson
  use fson_value_m
  use kinds_interface
  use variable_interface

  ! Define interfaces and attributes for module routines
  
  implicit none
  private
  public :: json_interface_nrecs
  public :: json_interface_read
  interface json_interface_read
     module procedure read_hwrf_bdy
     module procedure read_nmm_grid
     module procedure read_nmm_var
  end interface json_interface_read

  !-----------------------------------------------------------------------

contains

  !=======================================================================

  ! SUBROUTINE:

  ! json_interface_nrecs.f90

  ! DESCRIPTION:

  ! This subroutine defines the number of records within the external
  ! JSON formatted file specified by the user using the FSON API
  ! utilities.

  ! INPUT VARIABLES:

  ! * filename; a FORTRAN character string specifying the full-path to
  !   a JSON formatted file.

  !-----------------------------------------------------------------------

  subroutine json_interface_nrecs(filename,json_size)

    ! Define variables passed to routine

    character(len=500)                                                  :: filename
    integer                                                             :: json_size

    ! Define variables computed within routine
    
    type(fson_value),                                           pointer :: json_file

    !=====================================================================

    ! Define local variables

    json_file => fson_parse(trim(adjustl(filename)))
    json_size = fson_value_count(json_file)

    !=====================================================================

  end subroutine json_interface_nrecs

  !=======================================================================

  ! SUBROUTINE:

  ! read_hwrf_bdy.f90

  ! DESCRIPTION:

  ! This subroutine parses a JSON formatted file and defines elements
  ! of the hwrf_json_bdy_struct variable necessary to parse the JSON
  ! formatted variable table (vtable) specfied by the user.

  ! INPUT VARIABLES:

  ! * filename; a FORTRAN character string specifying the full-path to
  !   a JSON formatted file.

  ! * json_size; a FORTRAN integer to define the size of the
  !   hwrf_json_bdy_struct variable array (below).

  ! * json; a FORTRAN hwrf_json_bdy_struct variable.

  ! OUTPUT VARIABLES:

  ! * json_size; a FORTRAN integer defining the size of the
  !   hwrf_json_bdy_struct variable array.

  ! * json; a FORTRAN hwrf_json_bdy_struct variable containing the
  !   contents of the user JSON formatted file.

  !-----------------------------------------------------------------------

  subroutine read_hwrf_bdy(filename,json_size,json)

    ! Define variables passed to routine

    integer                                                             :: json_size
    type(hwrf_json_bdy_struct)                                          :: json(json_size)
    character(len=500)                                                  :: filename

    ! Define variables computed within routine

    type(fson_value),                                           pointer :: json_file  
    type(fson_value),                                           pointer :: json_item

    ! Define counting variables

    integer                                                             :: i

    !=====================================================================

    ! Define local variables

    json_file => fson_parse(trim(adjustl(filename)))

    ! Loop through local variable

    do i = 1, json_size

       ! Define local variables

       json_item => fson_value_get(json_file,i)
       call fson_get(json_item,'variable_name',json(i)%variable_name)
       call fson_get(json_item,'gridtype',json(i)%gridtype)
       call fson_get(json_item,'levtype',json(i)%levtype)
       call fson_get(json_item,'rezero',json(i)%rezero)

    end do ! do i = 1, json_size

    ! Define local variables

    call fson_destroy(json_file)

    !=====================================================================

  end subroutine read_hwrf_bdy
  
  !=======================================================================

  ! SUBROUTINE:

  ! read_nmm_grid.f90

  ! DESCRIPTION:

  ! This subroutine parses a JSON formatted file and defines elements
  ! of the nmm_json_grid_struct variable necessary to parse the JSON
  ! formatted variable table (vtable) specfied by the user.

  ! INPUT VARIABLES:

  ! * filename; a FORTRAN character string specifying the full-path to
  !   a JSON formatted file.

  ! * json_in; a FORTRAN nmm_json_grid_struct variable.

  ! * json_out; a FORTRAN nmm_json_grid_struct variable.

  ! OUTPUT VARIABLES:

  ! * json_in; a FORTRAN nmm_json_grid_struct variable containing the
  !   contents of the input-grid from the user-specified JSON
  !   formatted file.

  ! * json_out; a FORTRAN nmm_json_grid_struct variable containing the
  !   contents of the output-grid from the user-specified JSON
  !   formatted file.

  !-----------------------------------------------------------------------

  subroutine read_nmm_grid(filename,json_in,json_out)

    ! Define variables passed to routine

    type(nmm_json_grid_struct)                                          :: json_in
    type(nmm_json_grid_struct)                                          :: json_out
    character(len=500)                                                  :: filename

    ! Define variables computed within routine

    type(fson_value),                                           pointer :: json_data
    type(fson_value),                                           pointer :: json_item
    type(nmm_json_grid_struct)                                          :: json
    character(len=100)                                                  :: grid_name
    character(len=1)                                                    :: grid_id
    integer                                                             :: json_size

    ! Define counting variables

    integer                                                             :: i

    !=====================================================================

    ! Define local variables

    call json_interface_nrecs(filename,json_size)
    json_data => fson_parse(trim(adjustl(filename)))

    ! Loop through local variable
    
    do i = 1, json_size

       ! Define local variables

       json_item => fson_value_get(json_data,i)
       call fson_get(json_item,'grid_name',grid_name)
       call fson_get(json_item,'arakawa_type',json%arakawa_type)
       call fson_get(json_item,'lmsk_name',json%lmsk_name)
       call fson_get(json_item,'m_lat_name',json%m_lat_name)
       call fson_get(json_item,'m_lat_radians',json%m_lat_radians)
       call fson_get(json_item,'m_lon_name',json%m_lon_name)
       call fson_get(json_item,'m_lon_radians',json%m_lon_radians)
       call fson_get(json_item,'topo_name',json%topo_name)
       call fson_get(json_item,'u_lat_name',json%u_lat_name)
       call fson_get(json_item,'u_lat_radians',json%u_lat_radians)
       call fson_get(json_item,'u_lon_name',json%u_lon_name)
       call fson_get(json_item,'u_lon_radians',json%u_lon_radians)
       call fson_get(json_item,'v_lat_name',json%v_lat_name)
       call fson_get(json_item,'v_lat_radians',json%v_lat_radians)
       call fson_get(json_item,'v_lon_name',json%v_lon_name)
       call fson_get(json_item,'v_lon_radians',json%v_lon_radians)
       write(grid_id,'(a1)') trim(adjustl(grid_name))
       if((grid_id .eq. 'i') .or. (grid_id .eq. 'I')) json_in  = json
       if((grid_id .eq. 'o') .or. (grid_id .eq. 'O')) json_out = json
       
    end do ! do i = 1, json_size

    ! Define local variables

    call fson_destroy(json_data)

    !=====================================================================

  end subroutine read_nmm_grid

  !=======================================================================

  ! SUBROUTINE:

  ! read_nmm_var.f90

  ! DESCRIPTION:

  ! This subroutine parses a JSON formatted file and defines elements
  ! of the nmm_json_var_struct variable necessary to parse the JSON
  ! formatted variable table (vtable) specfied by the user.

  ! INPUT VARIABLES:

  ! * filename; a FORTRAN character string specifying the full-path to
  !   a JSON formatted file.

  ! * json_size; a FORTRAN integer to define the size of the
  !   nmm_json_var_struct variable array (below).

  ! * json; a FORTRAN nmm_json_var_struct variable.

  ! OUTPUT VARIABLES:

  ! * json_size; a FORTRAN integer defining the size of the
  !   nmm_json_struct variable array.

  ! * json; a FORTRAN nmm_json_var_struct variable containing the
  !   contents of the user JSON formatted file.

  !-----------------------------------------------------------------------

  subroutine read_nmm_var(filename,json_size,json)

    ! Define variables passed to routine

    integer                                                             :: json_size
    type(nmm_json_var_struct)                                           :: json(json_size)
    character(len=500)                                                  :: filename

    ! Define variables computed within routine

    type(fson_value),                                           pointer :: json_file  
    type(fson_value),                                           pointer :: json_item
    
    ! Define counting variables

    integer                                                             :: i

    !=====================================================================

    ! Define local variables

    json_file => fson_parse(trim(adjustl(filename)))

    ! Loop through local variable

    do i = 1, json_size

       ! Define local variables

       json_item => fson_value_get(json_file,i)
       call fson_get(json_item,'variable_name',json(i)%variable_name)
       call fson_get(json_item,'remap_variable_name',                      &
            & json(i)%remap_variable_name)
       call fson_get(json_item,'levtype',json(i)%levtype)
       call fson_get(json_item,'intrptype',json(i)%intrptype)
       call fson_get(json_item,'gridtype',json(i)%gridtype)
       call fson_get(json_item,'z_staggered',json(i)%z_staggered)
       call fson_get(json_item,'replace',json(i)%replace)
       
    end do ! do i = 1, json_size

    ! Define local variables

    call fson_destroy(json_file)

    !=====================================================================

  end subroutine read_nmm_var
  
  !=======================================================================

end module json_interface
