module fileio_interface

  !=======================================================================

  !$$$ MODULE DOCUMENTATION BLOCK
  
  ! da-update-analysis :: fileio_interface
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

  ! Define associate modules and subroutines

  use constants_interface
  use json_interface
  use kinds_interface
  use namelist_interface
  use netcdf_interface
  use time_methods_interface
  use variable_interface
  
  ! Define interfaces and attributes for module routines
  
  implicit none
  private
  public :: fileio_interface_read
  public :: fileio_interface_write
  interface fileio_interface_read
     module procedure hwrf_json_bdy_read
     module procedure hwrf_netcdf_bdy_read
     module procedure nmm_json_grid_read
     module procedure nmm_json_var_read
     module procedure nmm_netcdf_grid_read 
     module procedure nmm_netcdf_var_read
     module procedure tcv_read
  end interface fileio_interface_read
  interface fileio_interface_read_init
     module procedure hwrf_netcdf_bdy_init
     module procedure nmm_json_read_init
     module procedure nmm_netcdf_read_grid_init
     module procedure nmm_netcdf_read_var_init
  end interface fileio_interface_read_init
  interface fileio_interface_write
     module procedure hwrf_netcdf_vert_grid_write
     module procedure hwrf_netcdf_bdy_write
     module procedure nmm_netcdf_var_write
  end interface fileio_interface_write

  !-----------------------------------------------------------------------

contains

  !=======================================================================

  ! SUBROUTINE:

  ! hwrf_json_bdy_read.f90

  ! DESCRIPTION:

  ! This subroutine reads the contents of a user specified JSON
  ! formatted file and fills the json variable array.

  ! INPUT VARIABLES:

  ! * filename; a FORTRAN character string specifying the full-path to
  !   the JSON file.

  ! * nvars; a FORTRAN integer to contain the total number of
  !   variables in the JSON file.

  ! * json; a FORTRAN hwrf_json_bdy_struct variable array.

  ! OUTPUT VARIABLES:

  ! * nvars; a FORTRAN integer specifying the total number of variable
  !   records within the JSON file.

  ! * json; a FORTRAN hwrf_json_bdy_struct variable containing the
  !   contents of the user-specified JSON formatted file.

  !-----------------------------------------------------------------------

  subroutine hwrf_json_bdy_read(filename,nvars,json)

    ! Define variables passed to routine

    type(hwrf_json_bdy_struct), dimension(:),               allocatable :: json
    character(len=500)                                                  :: filename
    integer                                                             :: nvars

    !=====================================================================

    ! Define local variables

    call json_interface_nrecs(filename,nvars)

    ! Allocate memory for local variables

    if(.not. allocated(json)) allocate(json(nvars))

    ! Define local variables

    call json_interface_read(filename,nvars,json)

    !=====================================================================

  end subroutine hwrf_json_bdy_read

  !=======================================================================

  ! SUBROUTINE:

  ! hwrf_netcdf_vert_grid_write.f90

  ! DESCRIPTION:
  
  ! This subroutine writes the Hurricane WRF (HWRF) prognostic
  ! vertical level variables to a user specified netcdf file.

  ! INPUT VARIABLES:

  ! * filename; a FORTRAN character string specifying the full-path to
  !   the netcdf file.

  ! * vert_grid; a FORTRAN hwrf_vert_grid_struct variable.

  !-----------------------------------------------------------------------

  subroutine hwrf_netcdf_vert_grid_write(filename,vert_grid)

    ! Define variables passed to routine

    type(hwrf_vert_grid_struct)                                         :: vert_grid
    character(len=500)                                                  :: filename

    ! Define variables computed within routine

    character(len=100)                                                  :: varname
    real(r_kind),               dimension(:,:,:,:),         allocatable :: workgrid

    ! Define counting variables

    integer                                                             :: i

    !=====================================================================

    ! Allocate memory for local variables

    if(.not. allocated(workgrid))                                          &
         & allocate(workgrid(vert_grid%nx,vert_grid%ny,1,1))

    ! Define local variables

    workgrid(:,:,1,1) = reshape(vert_grid%pd,shape(workgrid(:,:,1,1)))
    varname           = 'PD'
    call netcdf_interface_putvar(filename,varname,workgrid)

    ! Deallocate memory for local variables

    if(allocated(workgrid)) deallocate(workgrid)

    ! Allocate memory for local variables

    if(.not. allocated(workgrid))                                          &
         & allocate(workgrid(vert_grid%nx,vert_grid%ny,                    &
         & (vert_grid%nz + 1),1))

    ! Loop through local variable

    do i = 1, (vert_grid%nz + 1)

       ! Define local variables

       workgrid(:,:,i,1) = reshape(vert_grid%pint(:,i),                    &
            & shape(workgrid(:,:,i,1)))

    end do ! do i = 1, (vert_grid%nz + 1)

    ! Define local variables

    varname = 'PINT'
    call netcdf_interface_putvar(filename,varname,workgrid)

    ! Deallocate memory for local variables

    if(allocated(workgrid)) deallocate(workgrid)

    !=====================================================================

  end subroutine hwrf_netcdf_vert_grid_write

  !=======================================================================

  ! SUBROUTINE:

  ! hwrf_netcdf_bdy_init.f90

  ! DESCRIPTION:

  ! This subroutine initializes the netcdf_struct variable arrays in
  ! preparation for netcdf file reading.
  
  ! INPUT VARIABLES:

  ! * filename; a FORTRAN character string specifying the full-path to
  !   the netcdf file.

  ! * json; a FORTRAN nmm_json_var_struct variable.

  ! * netcdf; a FORTRAN hwrf_bdy_grid_struct variable.

  ! OUTPUT VARIABLES:

  ! * netcdf; 

  !-----------------------------------------------------------------------

  subroutine hwrf_netcdf_bdy_init(filename,json,netcdf)

    ! Define variables passed to routine

    type(hwrf_bdy_grid_struct)                                          :: netcdf
    type(nmm_json_var_struct)                                           :: json
    character(len=500)                                                  :: filename

    ! Define variables computed within routine

    character(len=100)                                                  :: dimname

    !=====================================================================

    ! Define local variables

    dimname = 'bdy_width'
    call netcdf_interface_getdim(filename,dimname,netcdf%nb)
    dimname = 'Time'
    call netcdf_interface_getdim(filename,dimname,netcdf%nt)
    dimname = 'west_east'
    call netcdf_interface_getdim(filename,dimname,netcdf%nx)
    dimname = 'south_north'
    call netcdf_interface_getdim(filename,dimname,netcdf%ny)

    ! Check local variable and proceed accordingly

    if(json%levtype .eq. 1) then

       ! Define local variables

       netcdf%nz = 1

    end if ! if(json%levtype .eq. 1)

    ! Check local variable and proceed accordingly

    if(json%levtype .eq. 3) then

       ! Define local variables

       dimname = 'bottom_top'
       call netcdf_interface_getdim(filename,dimname,netcdf%nz)

    end if ! if(json%levtype .eq. 3)

    ! Define local variables

    call netcdf_interface_close()
    call variable_interface_setup_struct(netcdf)

    !=====================================================================

  end subroutine hwrf_netcdf_bdy_init

  !=======================================================================

  ! SUBROUTINE:

  ! hwrf_netcdf_bdy_read.f90

  ! DESCRIPTION:

  ! This subroutine reads a user specified variable from the user
  ! specified netcdf file.
  
  ! INPUT VARIABLES:

  ! * filename; a FORTRAN character string specifying the full-path to
  !   the netcdf file.

  ! * json; a FORTRAN nmm_json_var_struct variable.

  ! * hwrf_bdy_grid; a FORTRAN hwrf_bdy_grid_struct variable.

  ! OUTPUT VARIABLES:

  ! * hwrf_bdy_grid; a FORTRAN hwrf_bdy_grid_struct variable
  !   containing the user requested boundary-condition variable
  !   collected from the user specified file.

  !-----------------------------------------------------------------------

  subroutine hwrf_netcdf_bdy_read(filename,json,hwrf_bdy_grid)

    ! Define variables passed to routine

    type(hwrf_bdy_grid_struct)                                          :: hwrf_bdy_grid
    type(nmm_json_var_struct)                                           :: json
    character(len=500)                                                  :: filename

    ! Define variables computed within routine

    type(netcdf_var_struct)                                             :: netcdf
    character(len=19),          dimension(:),               allocatable :: times
    character(len=100)                                                  :: varname

    !=====================================================================

    ! Define local variables

    call fileio_interface_read_init(filename,json,hwrf_bdy_grid)

    ! Allocate memory for local variables

    if(.not. allocated(times)) allocate(times(hwrf_bdy_grid%nt))

    ! Define local variables

    varname = 'Times'
    call netcdf_interface_getvar(filename,varname,19,times)

    ! Compute local variables

    call time_methods_interface_dseconds(times(2),times(1),                &
         & hwrf_bdy_grid%bdyfrq)

    ! Deallocate memory for local variables

    if(allocated(times)) deallocate(times)

    ! Define local variables

    call variable_interface_setup_struct(hwrf_bdy_grid)
    netcdf%nx = hwrf_bdy_grid%ny

    ! Check local variable and proceed accordingly

    if(json%levtype .eq. 1) then

       ! Define local variables

       netcdf%nt = 1
       netcdf%ny = hwrf_bdy_grid%nb
       netcdf%nz = hwrf_bdy_grid%nt

    end if ! if(json%levtype .eq. 1)

    ! Check local variable and proceed accordingly

    if(json%levtype .eq. 3) then

       ! Define local variables

       netcdf%nt = hwrf_bdy_grid%nt
       netcdf%ny = hwrf_bdy_grid%nz
       netcdf%nz = hwrf_bdy_grid%nb

    end if ! if(json%levtype .eq. 3)

    ! Define local variables

    call variable_interface_setup_struct(netcdf)

    ! Check local variable and proceed accordingly

    if(json%levtype .eq. 1) then    

       ! Define local variables
       
       write(varname,500) trim(adjustl(json%variable_name))
       call netcdf_interface_getvar(filename,varname,netcdf%var)
       hwrf_bdy_grid%bxs(1,:,:,:,1)                 =                      &
            & netcdf%var(1:netcdf%nx,1:netcdf%ny,1:netcdf%nt,1)
       write(varname,501) trim(adjustl(json%variable_name))
       call netcdf_interface_getvar(filename,varname,netcdf%var)
       hwrf_bdy_grid%bxe(hwrf_bdy_grid%nx,:,:,:,1)  =                      &
            & netcdf%var(1:netcdf%nx,1:netcdf%ny,1:netcdf%nt,1)
       write(varname,502) trim(adjustl(json%variable_name))
       call netcdf_interface_getvar(filename,varname,netcdf%var)
       hwrf_bdy_grid%btxs(1,:,:,:,1)                =                      &
            & netcdf%var(1:netcdf%nx,1:netcdf%ny,1:netcdf%nt,1)
       write(varname,503) trim(adjustl(json%variable_name))
       call netcdf_interface_getvar(filename,varname,netcdf%var)
       hwrf_bdy_grid%btxe(hwrf_bdy_grid%nz,:,:,:,1) =                      &
            & netcdf%var(1:netcdf%nx,1:netcdf%ny,1:netcdf%nt,1)

    end if ! if(json%levtype .eq. 1)
    
    ! Check local variable and proceed accordingly

    if(json%levtype .eq. 3) then    

       ! Define local variables
       
       write(varname,500) trim(adjustl(json%variable_name))
       call netcdf_interface_getvar(filename,varname,netcdf%var)
       hwrf_bdy_grid%bxs(1,:,:,:,1)                 =                      &
            & netcdf%var(1:netcdf%nx,1:netcdf%ny,1:netcdf%nz,1)
       write(varname,501) trim(adjustl(json%variable_name))
       call netcdf_interface_getvar(filename,varname,netcdf%var)
       hwrf_bdy_grid%bxe(hwrf_bdy_grid%nx,:,:,:,1)  =                      &
            & netcdf%var(1:netcdf%nx,1:netcdf%ny,1:netcdf%nz,1)
       write(varname,502) trim(adjustl(json%variable_name))
       call netcdf_interface_getvar(filename,varname,netcdf%var)
       hwrf_bdy_grid%btxs(1,:,:,:,1)                =                      &
            & netcdf%var(1:netcdf%nx,1:netcdf%ny,1:netcdf%nz,1)
       write(varname,503) trim(adjustl(json%variable_name))
       call netcdf_interface_getvar(filename,varname,netcdf%var)
       hwrf_bdy_grid%btxe(hwrf_bdy_grid%nz,:,:,:,1) =                      &
            & netcdf%var(1:netcdf%nx,1:netcdf%ny,1:netcdf%nz,1)

    end if ! if(json%levtype .eq. 3)

    ! Deallocate memory for local variables

    call variable_interface_cleanup_struct(netcdf)

    ! Define local variables

    netcdf%nt = hwrf_bdy_grid%nt
    netcdf%nx = hwrf_bdy_grid%nx

    ! Check local variable and proceed accordingly

    if(json%levtype .eq. 1) then

       ! Define local variables

       netcdf%nt = 1
       netcdf%ny = hwrf_bdy_grid%nb
       netcdf%nz = hwrf_bdy_grid%nt

    end if ! if(json%levtype .eq. 1)

    ! Check local variable and proceed accordingly

    if(json%levtype .eq. 3) then

       ! Define local variables

       netcdf%nt = hwrf_bdy_grid%nt
       netcdf%ny = hwrf_bdy_grid%nz
       netcdf%nz = hwrf_bdy_grid%nb

    end if ! if(json%levtype .eq. 3)

    ! Check local variable and proceed accordingly

    if(json%levtype .eq. 1) then
       
       ! Define local variables

       call variable_interface_setup_struct(netcdf)
       write(varname,504) trim(adjustl(json%variable_name))
       call netcdf_interface_getvar(filename,varname,netcdf%var)
       hwrf_bdy_grid%bys(:,1,:,:,1)                 =                      &
            & netcdf%var(1:netcdf%nx,1:netcdf%ny,1:netcdf%nt,1)
       write(varname,505) trim(adjustl(json%variable_name))
       call netcdf_interface_getvar(filename,varname,netcdf%var)
       hwrf_bdy_grid%bye(:,hwrf_bdy_grid%ny,:,:,1)  =                      &
            & netcdf%var(1:netcdf%nx,1:netcdf%ny,1:netcdf%nt,1)
       write(varname,506) trim(adjustl(json%variable_name))
       call netcdf_interface_getvar(filename,varname,netcdf%var)
       hwrf_bdy_grid%btys(:,1,:,:,1)                =                      &
            & netcdf%var(1:netcdf%nx,1:netcdf%ny,1:netcdf%nt,1)
       write(varname,507) trim(adjustl(json%variable_name))
       call netcdf_interface_getvar(filename,varname,netcdf%var)
       hwrf_bdy_grid%btye(:,hwrf_bdy_grid%ny,:,:,1) =                      &
            & netcdf%var(1:netcdf%nx,1:netcdf%ny,1:netcdf%nt,1)

    end if ! if(json%levtype .eq. 1)

    ! Check local variable and proceed accordingly

    if(json%levtype .eq. 3) then
       
       ! Define local variables

       call variable_interface_setup_struct(netcdf)
       write(varname,504) trim(adjustl(json%variable_name))
       call netcdf_interface_getvar(filename,varname,netcdf%var)
       hwrf_bdy_grid%bys(:,1,:,:,1)                 =                      &
            & netcdf%var(1:netcdf%nx,1:netcdf%ny,1:netcdf%nz,1)
       write(varname,505) trim(adjustl(json%variable_name))
       call netcdf_interface_getvar(filename,varname,netcdf%var)
       hwrf_bdy_grid%bye(:,hwrf_bdy_grid%ny,:,:,1)  =                      &
            & netcdf%var(1:netcdf%nx,1:netcdf%ny,1:netcdf%nz,1)
       write(varname,506) trim(adjustl(json%variable_name))
       call netcdf_interface_getvar(filename,varname,netcdf%var)
       hwrf_bdy_grid%btys(:,1,:,:,1)                =                      &
            & netcdf%var(1:netcdf%nx,1:netcdf%ny,1:netcdf%nz,1)
       write(varname,507) trim(adjustl(json%variable_name))
       call netcdf_interface_getvar(filename,varname,netcdf%var)
       hwrf_bdy_grid%btye(:,hwrf_bdy_grid%ny,:,:,1) =                      &
            & netcdf%var(1:netcdf%nx,1:netcdf%ny,1:netcdf%nz,1)

    end if ! if(json%levtype .eq. 3)

    ! Deallocate memory for local variables

    call variable_interface_cleanup_struct(netcdf)

    ! Define local variables

500 format(a,'_BXS')
501 format(a,'_BXE')
502 format(a,'_BTXS')
503 format(a,'_BTXE')
504 format(a,'_BYS')
505 format(a,'_BYE')
506 format(a,'_BTYS')
507 format(a,'_BTYE')

    !=====================================================================

  end subroutine hwrf_netcdf_bdy_read

  !=======================================================================

  ! SUBROUTINE:

  ! hwrf_netcdf_bdy_write.f90

  ! DESCRIPTION:
  
  ! This subroutine writes the updated HWRF boundary condition
  ! variables to the user specified netcdf file.

  ! INPUT VARIABLES:

  ! * filename; a FORTRAN character string specifying the full-path to
  !   the netcdf file.

  ! * hwrf_bdy_grid; a FORTRAN hwrf_bdy_grid_struct variable.

  !-----------------------------------------------------------------------

  subroutine hwrf_netcdf_bdy_write(filename,hwrf_json_bdy,hwrf_bdy_grid)

    ! Define variables passed to routine

    type(hwrf_bdy_grid_struct)                                          :: hwrf_bdy_grid
    type(hwrf_json_bdy_struct)                                          :: hwrf_json_bdy
    character(len=500)                                                  :: filename

    ! Define variables computed within routine

    type(hwrf_bdy_grid_struct)                                          :: hwrf_bdy_grid_local
    type(netcdf_var_struct)                                             :: netcdf
    type(nmm_json_var_struct)                                           :: json
    character(len=100)                                                  :: varname
    real(r_kind),               dimension(:,:,:,:),         allocatable :: workgrid

    ! Define counting variables

    integer                                                             :: i

    !=====================================================================

    ! Define local variables

    json%levtype = hwrf_json_bdy%levtype
    call fileio_interface_read_init(filename,json,hwrf_bdy_grid_local)
    call variable_interface_setup_struct(hwrf_bdy_grid_local)
    netcdf%nx    = hwrf_bdy_grid_local%ny

    ! Check local variable and proceed accordingly

    if(json%levtype .eq. 1) then

       ! Define local variables

       netcdf%nt = 1
       netcdf%ny = hwrf_bdy_grid_local%nb
       netcdf%nz = hwrf_bdy_grid_local%nt

       ! Allocate memory for local variables

       if(.not. allocated(workgrid))                                       &
            & allocate(workgrid(hwrf_bdy_grid_local%ny,1,                  &
            & hwrf_bdy_grid_local%nb,hwrf_bdy_grid_local%nt))

    end if ! if(json%levtype .eq. 1)

    ! Check local variable and proceed accordingly

    if(json%levtype .eq. 3) then

       ! Define local variables

       netcdf%nt = hwrf_bdy_grid_local%nt
       netcdf%ny = hwrf_bdy_grid_local%nz
       netcdf%nz = hwrf_bdy_grid_local%nb

       ! Allocate memory for local variables

       if(.not. allocated(workgrid))                                       &
            & allocate(workgrid(hwrf_bdy_grid_local%ny,                    &
            & hwrf_bdy_grid_local%nz,hwrf_bdy_grid_local%nb,               &
            & hwrf_bdy_grid_local%nt))

    end if ! if(json%levtype .eq. 3)    

    ! Define local variables

    call variable_interface_setup_struct(netcdf)
    write(varname,500) trim(adjustl(hwrf_json_bdy%variable_name))
    call netcdf_interface_getvar(filename,varname,netcdf%var)

    ! Check local variable and proceed accordingly

    if(json%levtype .eq. 1) then

       ! Define local variables

       hwrf_bdy_grid_local%bxs(1,:,1,:,:) = netcdf%var(:,:,:,1)
       hwrf_bdy_grid_local%bxs(1,:,1,:,1) = hwrf_bdy_grid%bxs(1,:,1,:,1)
       workgrid(:,1,:,:)                  =                                &
            & hwrf_bdy_grid_local%bxs(1,:,1,:,:)

    end if ! if(json%levtype .eq. 1)

    ! Check local variable and proceed accordingly

    if(json%levtype .eq. 3) then

       ! Define local variables

       hwrf_bdy_grid_local%bxs(1,:,:,:,:) = netcdf%var
       hwrf_bdy_grid_local%bxs(1,:,:,:,1) = hwrf_bdy_grid%bxs(1,:,:,:,1)
       workgrid                           =                                &
            & hwrf_bdy_grid_local%bxs(1,:,:,:,:)

    end if ! if(json%levtype .eq. 3)

    ! Define local variables

    call netcdf_interface_putvar(filename,varname,workgrid)
    write(varname,501) trim(adjustl(hwrf_json_bdy%variable_name))
    call netcdf_interface_getvar(filename,varname,netcdf%var)

    ! Check local variable and proceed accordingly

    if(json%levtype .eq. 1) then

       ! Define local variables

       hwrf_bdy_grid_local%bxe(hwrf_bdy_grid_local%nx,:,1,:,:)  =          &
            & netcdf%var(:,:,:,1)
       hwrf_bdy_grid_local%bxe(hwrf_bdy_grid_local%nx,:,1,:,1)  =          &
            & hwrf_bdy_grid%bxe(hwrf_bdy_grid%nx,:,1,:,1)
       workgrid(:,1,:,:)                                        =          &
            & hwrf_bdy_grid_local%bxe(hwrf_bdy_grid%nx,:,1,:,:)

    end if ! if(json%levtype .eq. 1)

    ! Check local variable and proceed accordingly

    if(json%levtype .eq. 3) then

       ! Define local variables

       hwrf_bdy_grid_local%bxe(hwrf_bdy_grid_local%nx,:,:,:,:)  =          &
            & netcdf%var
       hwrf_bdy_grid_local%bxe(hwrf_bdy_grid_local%nx,:,:,:,1)  =          &
            & hwrf_bdy_grid%bxe(hwrf_bdy_grid%nx,:,:,:,1)
       workgrid                                                 =          &
            & hwrf_bdy_grid_local%bxe(hwrf_bdy_grid%nx,:,:,:,:)

    end if ! if(json%levtype .eq. 3)

    ! Define local variables

    call netcdf_interface_putvar(filename,varname,workgrid)
    write(varname,502) trim(adjustl(hwrf_json_bdy%variable_name))
    call netcdf_interface_getvar(filename,varname,netcdf%var)

    ! Check local variable and proceed accordingly

    if(json%levtype .eq. 1) then

       ! Define local variables

       hwrf_bdy_grid_local%btxs(1,:,1,:,:) = netcdf%var(:,:,:,1)
       hwrf_bdy_grid_local%btxs(1,:,1,:,1) =                               &
            & hwrf_bdy_grid%btxs(1,:,1,:,1)
       workgrid(:,1,:,:)                   =                               &
            & hwrf_bdy_grid_local%btxs(1,:,1,:,:)

    end if ! if(json%levtype .eq. 3)

    ! Check local variable and proceed accordingly

    if(json%levtype .eq. 3) then

       ! Define local variables

       hwrf_bdy_grid_local%btxs(1,:,:,:,:) = netcdf%var
       hwrf_bdy_grid_local%btxs(1,:,:,:,1) =                               &
            & hwrf_bdy_grid%btxs(1,:,:,:,1)
       workgrid                            =                               &
            & hwrf_bdy_grid_local%btxs(1,:,:,:,:)

    end if ! if(json%levtype .eq. 3)

    ! Define local variables

    call netcdf_interface_putvar(filename,varname,workgrid)
    write(varname,503) trim(adjustl(hwrf_json_bdy%variable_name))
    call netcdf_interface_getvar(filename,varname,netcdf%var)

    ! Check local variable and proceed accordingly

    if(json%levtype .eq. 1) then

       ! Define local variables

       hwrf_bdy_grid_local%btxe(hwrf_bdy_grid_local%nx,:,1,:,:) =          &
            & netcdf%var(:,:,:,1)
       hwrf_bdy_grid_local%btxe(hwrf_bdy_grid_local%nx,:,1,:,1) =          &
            & hwrf_bdy_grid%btxe(hwrf_bdy_grid%nx,:,1,:,1)
       workgrid(:,1,:,:)                                        =          &
            & hwrf_bdy_grid_local%btxe(hwrf_bdy_grid%nx,:,1,:,:)

    end if ! if(json%levtype .eq. 1)

    ! Check local variable and proceed accordingly

    if(json%levtype .eq. 3) then

       ! Define local variables

       hwrf_bdy_grid_local%btxe(hwrf_bdy_grid_local%nx,:,:,:,:) =          &
            & netcdf%var
       hwrf_bdy_grid_local%btxe(hwrf_bdy_grid_local%nx,:,:,:,1) =          &
            & hwrf_bdy_grid%btxe(hwrf_bdy_grid%nx,:,:,:,1)
       workgrid                                                 =          &
            & hwrf_bdy_grid_local%btxe(hwrf_bdy_grid%nx,:,:,:,:)

    end if ! if(json%levtype .eq. 3)

    ! Define local variables

    call netcdf_interface_putvar(filename,varname,workgrid)

    ! Deallocate memory for local variables

    if(allocated(workgrid)) deallocate(workgrid)
    call variable_interface_cleanup_struct(netcdf)

    ! Define local variables

    netcdf%nx = hwrf_bdy_grid_local%nx

    ! Check local variable and proceed accordingly

    if(json%levtype .eq. 1) then

       ! Define local variables

       netcdf%nt = 1
       netcdf%ny = hwrf_bdy_grid_local%nb
       netcdf%nz = hwrf_bdy_grid_local%nt

       ! Allocate memory for local variables

       if(.not. allocated(workgrid))                                       &
            & allocate(workgrid(hwrf_bdy_grid_local%nx,1,                  &
            & hwrf_bdy_grid_local%nb,hwrf_bdy_grid_local%nt))

    end if ! if(json%levtype .eq. 1)

    ! Check local variable and proceed accordingly

    if(json%levtype .eq. 3) then

       ! Define local variables

       netcdf%nt = hwrf_bdy_grid_local%nt
       netcdf%ny = hwrf_bdy_grid_local%nz
       netcdf%nz = hwrf_bdy_grid_local%nb

       ! Allocate memory for local variables

       if(.not. allocated(workgrid))                                       &
            & allocate(workgrid(hwrf_bdy_grid_local%nx,                    &
            & hwrf_bdy_grid_local%nz,hwrf_bdy_grid_local%nb,               &
            & hwrf_bdy_grid_local%nt))

    end if ! if(json%levtype .eq. 3)    

    ! Define local variables

    call variable_interface_setup_struct(netcdf)
    write(varname,504) trim(adjustl(hwrf_json_bdy%variable_name))
    call netcdf_interface_getvar(filename,varname,netcdf%var)

    ! Check local variable and proceed accordingly

    if(json%levtype .eq. 1) then

       ! Define local variables

       hwrf_bdy_grid_local%bys(:,1,1,:,:) = netcdf%var(:,:,:,1)
       hwrf_bdy_grid_local%bys(:,1,1,:,1) = hwrf_bdy_grid%bys(:,1,1,:,1)
       workgrid(:,1,:,:)                  =                                &
            & hwrf_bdy_grid_local%bys(:,1,1,:,:)

    end if ! if(json%levtype .eq. 1)

    ! Check local variable and proceed accordingly

    if(json%levtype .eq. 3) then

       ! Define local variables

       hwrf_bdy_grid_local%bys(:,1,:,:,:) = netcdf%var
       hwrf_bdy_grid_local%bys(:,1,:,:,1) = hwrf_bdy_grid%bys(:,1,:,:,1)
       workgrid                           =                                &
            & hwrf_bdy_grid_local%bys(:,1,:,:,:)

    end if ! if(json%levtype .eq. 3)

    ! Define local variables

    call netcdf_interface_putvar(filename,varname,workgrid)
    write(varname,505) trim(adjustl(hwrf_json_bdy%variable_name))
    call netcdf_interface_getvar(filename,varname,netcdf%var)

    ! Check local variable and proceed accordingly

    if(json%levtype .eq. 1) then

       ! Define local variables

       hwrf_bdy_grid_local%bye(:,hwrf_bdy_grid_local%ny,1,:,:)  =          &
            & netcdf%var(:,:,:,1)
       hwrf_bdy_grid_local%bye(:,hwrf_bdy_grid_local%ny,1,:,1)  =          &
            & hwrf_bdy_grid%bye(:,hwrf_bdy_grid%ny,1,:,1)
       workgrid(:,1,:,:)                                        =          &
            & hwrf_bdy_grid_local%bye(:,hwrf_bdy_grid%ny,1,:,:)

    end if ! if(json%levtype .eq. 1)

    ! Check local variable and proceed accordingly

    if(json%levtype .eq. 3) then

       ! Define local variables

       hwrf_bdy_grid_local%bye(:,hwrf_bdy_grid_local%ny,:,:,:)  =          &
            & netcdf%var
       hwrf_bdy_grid_local%bye(:,hwrf_bdy_grid_local%ny,:,:,1)  =          &
            & hwrf_bdy_grid%bye(:,hwrf_bdy_grid%ny,:,:,1)
       workgrid                                                 =          &
            & hwrf_bdy_grid_local%bye(:,hwrf_bdy_grid%ny,:,:,:)

    end if ! if(json%levtype .eq. 3)

    ! Define local variables

    call netcdf_interface_putvar(filename,varname,workgrid)
    write(varname,506) trim(adjustl(hwrf_json_bdy%variable_name))
    call netcdf_interface_getvar(filename,varname,netcdf%var)

    ! Check local variable and proceed accordingly

    if(json%levtype .eq. 1) then

       ! Define local variables

       hwrf_bdy_grid_local%btys(:,1,1,:,:) = netcdf%var(:,:,:,1)
       hwrf_bdy_grid_local%btys(:,1,1,:,1) = hwrf_bdy_grid%btys(:,1,1,:,1)
       workgrid(:,1,:,:)                  =                                &
            & hwrf_bdy_grid_local%btys(:,1,1,:,:)

    end if ! if(json%levtype .eq. 1)

    ! Check local variable and proceed accordingly

    if(json%levtype .eq. 3) then

       ! Define local variables

       hwrf_bdy_grid_local%btys(:,1,:,:,:) = netcdf%var
       hwrf_bdy_grid_local%btys(:,1,:,:,1) = hwrf_bdy_grid%btys(:,1,:,:,1)
       workgrid                           =                                &
            & hwrf_bdy_grid_local%btys(:,1,:,:,:)

    end if ! if(json%levtype .eq. 3)

    ! Define local variables

    call netcdf_interface_putvar(filename,varname,workgrid)
    write(varname,507) trim(adjustl(hwrf_json_bdy%variable_name))
    call netcdf_interface_getvar(filename,varname,netcdf%var)

    ! Check local variable and proceed accordingly

    if(json%levtype .eq. 1) then

       ! Define local variables

       hwrf_bdy_grid_local%btye(:,hwrf_bdy_grid_local%ny,1,:,:)  =         &
            & netcdf%var(:,:,:,1)
       hwrf_bdy_grid_local%btye(:,hwrf_bdy_grid_local%ny,1,:,1)  =         &
            & hwrf_bdy_grid%btye(:,hwrf_bdy_grid%ny,1,:,1)
       workgrid(:,1,:,:)                                        =          &
            & hwrf_bdy_grid_local%btye(:,hwrf_bdy_grid%ny,1,:,:)

    end if ! if(json%levtype .eq. 1)

    ! Check local variable and proceed accordingly

    if(json%levtype .eq. 3) then

       ! Define local variables

       hwrf_bdy_grid_local%btye(:,hwrf_bdy_grid_local%ny,:,:,:) =          &
            & netcdf%var
       hwrf_bdy_grid_local%btye(:,hwrf_bdy_grid_local%ny,:,:,1) =          &
            & hwrf_bdy_grid%btye(:,hwrf_bdy_grid%ny,:,:,1)
       workgrid                                                 =          &
            & hwrf_bdy_grid_local%btye(:,hwrf_bdy_grid%ny,:,:,:)

    end if ! if(json%levtype .eq. 3)

    ! Define local variables

    call netcdf_interface_putvar(filename,varname,workgrid)

    ! Deallocate memory for local variables

    call variable_interface_cleanup_struct(hwrf_bdy_grid_local)
    call variable_interface_cleanup_struct(netcdf)

    ! Define local variables

500 format(a,'_BXS')
501 format(a,'_BXE')
502 format(a,'_BTXS')
503 format(a,'_BTXE')
504 format(a,'_BYS')
505 format(a,'_BYE')
506 format(a,'_BTYS')
507 format(a,'_BTYE')

    !=====================================================================

  end subroutine hwrf_netcdf_bdy_write

  !=======================================================================

  ! SUBROUTINE:

  ! nmm_json_grid_read.f90

  ! DESCRIPTION:

  ! This subroutine reads the contents of a user specified JSON
  ! formatted file and fills the json variable array.

  ! INPUT VARIABLES:

  ! * filename; a FORTRAN character string specifying the full-path to
  !   the JSON file.

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

  subroutine nmm_json_grid_read(filename,json_in,json_out)

    ! Define variables passed to routine

    type(nmm_json_grid_struct)                                          :: json_in
    type(nmm_json_grid_struct)                                          :: json_out
    character(len=500)                                                  :: filename

    !=====================================================================

    ! Define local variables

    call json_interface_read(filename,json_in,json_out)

    !=====================================================================

  end subroutine nmm_json_grid_read

  !=======================================================================

  ! SUBROUTINE:

  ! nmm_json_var_read.f90

  ! DESCRIPTION:

  ! This subroutine reads the contents of a user specified JSON
  ! formatted file and fills the json variable array.

  ! INPUT VARIABLES:

  ! * filename; a FORTRAN character string specifying the full-path to
  !   the JSON file.

  ! * nvars; a FORTRAN integer to contain the total number of
  !   variables in the JSON file.

  ! * json; a FORTRAN nmm_json_var_struct variable array.

  ! OUTPUT VARIABLES:

  ! * nvars; a FORTRAN integer specifying the total number of variable
  !   records within the JSON file.

  ! * json; a FORTRAN nmm_json_var_struct variable array containing
  !   the respective variable record attributes.

  !-----------------------------------------------------------------------

  subroutine nmm_json_var_read(filename,nvars,json)

    ! Define variables passed to routine

    type(nmm_json_var_struct),  dimension(:),               allocatable :: json
    character(len=500)                                                  :: filename
    integer                                                             :: nvars

    !=====================================================================

    ! Define local variables

    call fileio_interface_read_init(filename,nvars,json) 
    call json_interface_read(filename,nvars,json)

    !=====================================================================

  end subroutine nmm_json_var_read

  !=======================================================================

  ! SUBROUTINE:

  ! nmm_json_read_init.f90

  ! DESCRIPTION:

  ! This subroutine initializes the FORTRAN nmm_json_struct variable
  ! in accordance with the variables provided in the user specified
  ! JSON variable table filename.

  ! INPUT VARIABLES:

  ! * filename; a FORTRAN character string specifying the full-path to
  !   the JSON file.

  ! * nvars; a FORTRAN integer to contain the total number of
  !   variables in the JSON file.

  ! * json; a FORTRAN nmm_json_struct variable array.

  ! OUTPUT VARIABLES:

  ! * nvars; a FORTRAN integer specifying the total number of variable
  !   records within the JSON file.

  ! * json; a FORTRAN nmm_json_struct variable array containing the
  !   respective variable record attributes.

  !-----------------------------------------------------------------------

  subroutine nmm_json_read_init(filename,nvars,json)

    ! Define variables passed to routine

    type(nmm_json_var_struct),  dimension(:),               allocatable :: json
    character(len=500)                                                  :: filename
    integer                                                             :: nvars

    !=====================================================================
    
    ! Define local variables

    call json_interface_nrecs(filename,nvars)

    ! Allocate memory for local variables

    if(.not. allocated(json)) allocate(json(nvars))

    !=====================================================================

  end subroutine nmm_json_read_init

  !=======================================================================

  ! SUBROUTINE:

  ! nmm_netcdf_grid_read.f90

  ! DESCRIPTION:

  ! This subroutine reads a user specified variable and the respective
  ! grid projection (e.g., geographical location) attributes from the
  ! user specified file.

  ! INPUT VARIABLES:

  ! * filename; a FORTRAN character string specifying the full-path to
  !   the netcdf file.

  ! * json; a FORTRAN nmm_json_grid_struct variable.

  ! * netcdf; a FORTRAN netcdf_grid_struct variable.

  ! OUTPUT VARIABLES:

  ! * netcdf; a FORTRAN netcdf_grid_struct variable containing the
  !   user requested variable and the grid projection (e.g.,
  !   geographical location) attributes collected from the user
  !   specified file.

  !-----------------------------------------------------------------------

  subroutine nmm_netcdf_grid_read(filename,json,netcdf)

    ! Define variables passed to routine

    type(netcdf_grid_struct)                                            :: netcdf
    type(nmm_json_grid_struct)                                          :: json
    character(len=500)                                                  :: filename

    ! Define variables computed within routine

    character(len=100)                                                  :: attrname
    character(len=100)                                                  :: varname
    real(r_kind)                                                        :: dx
    real(r_kind)                                                        :: dy

    !=====================================================================

    ! Define local variables

    call fileio_interface_read_init(filename,json,netcdf)
    call variable_interface_setup_struct(netcdf)
    varname  = json%lmsk_name
    call netcdf_interface_getvar(filename,varname,netcdf%landmask)
    varname  = json%m_lat_name
    call netcdf_interface_getvar(filename,varname,netcdf%m_lat)
    varname  = json%m_lon_name
    call netcdf_interface_getvar(filename,varname,netcdf%m_lon)
    varname  = json%topo_name
    call netcdf_interface_getvar(filename,varname,netcdf%topo)
    varname  = json%u_lat_name
    call netcdf_interface_getvar(filename,varname,netcdf%u_lat)
    varname  = json%u_lon_name
    call netcdf_interface_getvar(filename,varname,netcdf%u_lon)
    varname  = json%v_lat_name
    call netcdf_interface_getvar(filename,varname,netcdf%v_lat)
    varname  = json%v_lon_name
    call netcdf_interface_getvar(filename,varname,netcdf%v_lon)
    attrname = 'DX'
    call netcdf_interface_getattr(filename,attrname,netcdf%dx)
    attrname = 'DY'
    call netcdf_interface_getattr(filename,attrname,netcdf%dy)

    !=====================================================================

  end subroutine nmm_netcdf_grid_read

  !=======================================================================

  ! SUBROUTINE:

  ! nmm_netcdf_var_read.f90

  ! DESCRIPTION:

  ! This subroutine reads a user specified variable from the user
  ! specified netcdf file.

  ! INPUT VARIABLES:

  ! * filename; a FORTRAN character string specifying the full-path to
  !   the netcdf file.

  ! * json; a FORTRAN nmm_json_struct variable.

  ! * netcdf; a FORTRAN netcdf_struct variable.

  ! OUTPUT VARIABLES:

  ! * netcdf; a FORTRAN netcdf_struct variable containing the user
  !   requested variable collected from the user specified file.

  !-----------------------------------------------------------------------

  subroutine nmm_netcdf_var_read(filename,json,netcdf)

    ! Define variables passed to routine

    type(netcdf_var_struct)                                             :: netcdf
    type(nmm_json_var_struct)                                           :: json
    character(len=500)                                                  :: filename

    ! Define variables computed within routine

    character(len=100)                                                  :: attrname
    character(len=100)                                                  :: varname
    
    !=====================================================================

    ! Define local variables

    call fileio_interface_read_init(filename,json,netcdf)
    call variable_interface_setup_struct(netcdf)
    varname  = json%variable_name
    call netcdf_interface_getvar(filename,varname,netcdf%var)
    attrname = 'CEN_LAT'
    call netcdf_interface_getattr(filename,attrname,netcdf%clat)
    attrname = 'CEN_LON'
    call netcdf_interface_getattr(filename,attrname,netcdf%clon)

    !=====================================================================

  end subroutine nmm_netcdf_var_read

  !=======================================================================

  ! SUBROUTINE:

  ! nmm_netcdf_var_write.f90

  ! DESCRIPTION:

  ! This subroutine writes a user specified variable to the user
  ! specified netcdf file.

  ! INPUT VARIABLES:

  ! * filename; a FORTRAN character string specifying the full-path to
  !   the netcdf file.

  ! * json; a FORTRAN nmm_json_struct variable.

  ! * netcdf; a FORTRAN netcdf_struct variable.

  ! OUTPUT VARIABLES:

  ! * netcdf; a FORTRAN netcdf_struct variable containing the user
  !   requested variable collected from the user specified file.

  !-----------------------------------------------------------------------

  subroutine nmm_netcdf_var_write(filename,json,netcdf)

    ! Define variables passed to routine

    type(netcdf_var_struct)                                             :: netcdf
    type(nmm_json_var_struct)                                           :: json
    character(len=500)                                                  :: filename

    ! Define variables computed within routine

    character(len=100)                                                  :: varname

    !=====================================================================

    ! Define local variables

    varname = json%variable_name
    call netcdf_interface_putvar(filename,varname,netcdf%var)

    !=====================================================================

  end subroutine nmm_netcdf_var_write

  !=======================================================================

  ! SUBROUTINE:

  ! nmm_netcdf_read_grid_init.f90

  ! DESCRIPTION:

  ! This subroutine initializes the netcdf_struct variable arrays in
  ! preparation for netcdf file reading.

  ! INPUT VARIABLES:

  ! * filename; a FORTRAN character string specifying the full-path to
  !   the netcdf file.

  ! * json; a FORTRAN nmm_json_grid_struct variable containing the
  !   grid attributes.

  ! * netcdf; a FORTRAN netcdf_grid_struct variable.

  ! OUTPUT VARIABLES:

  ! * netcdf; a FORTRAN netcdf_grid_struct variable with appropriate
  !   grid dimension variables defined and structure arrays allocated.

  !-----------------------------------------------------------------------

  subroutine nmm_netcdf_read_grid_init(filename,json,netcdf)

    ! Define variables passed to routine

    type(netcdf_grid_struct)                                            :: netcdf
    type(nmm_json_grid_struct)                                          :: json
    character(len=500)                                                  :: filename

    ! Define variables computed within routine

    character(len=100)                                                  :: dimname

    !=====================================================================

    ! Define local variables

    dimname = 'west_east'
    call netcdf_interface_getdim(filename,dimname,netcdf%nx)
    dimname = 'south_north'
    call netcdf_interface_getdim(filename,dimname,netcdf%ny)
    dimname = 'bottom_top'
    call netcdf_interface_getdim(filename,dimname,netcdf%nz)
    call netcdf_interface_close()

    ! Check local variable and proceed accordingly

    if((trim(adjustl(json%arakawa_type)) .eq. 'B') .or.                    &
         & (trim(adjustl(json%arakawa_type)) .eq. 'b')) then

       ! Define local variables

       netcdf%nxs = netcdf%nx
       netcdf%nys = netcdf%ny

    end if ! if((trim(adjustl(json%arakawa_type)) .eq. 'B')
           ! .or. (trim(adjustl(json%arakawa_type)) .eq. 'b'))

    ! Check local variable and proceed accordingly

    if((trim(adjustl(json%arakawa_type)) .eq. 'C') .or.                    &
         & (trim(adjustl(json%arakawa_type)) .eq. 'c')) then

       ! Define local variables

       netcdf%nxs = (netcdf%nx + 1)
       netcdf%nys = (netcdf%ny + 1)

    end if ! if((trim(adjustl(json%arakawa_type)) .eq. 'C')
           ! .or. (trim(adjustl(json%arakawa_type)) .eq. 'c'))

    !=====================================================================

  end subroutine nmm_netcdf_read_grid_init

  !=======================================================================

  ! SUBROUTINE:

  ! nmm_netcdf_read_var_init.f90

  ! DESCRIPTION:

  ! This subroutine initializes the netcdf_struct variable arrays in
  ! preparation for netcdf file reading.

  ! INPUT VARIABLES:

  ! * filename; a FORTRAN character string specifying the full-path to
  !   the netcdf file.

  ! * json; a FORTRAN nmm_json_struct variable.

  ! * netcdf; a FORTRAN netcdf_var_struct variable.

  ! OUTPUT VARIABLES:

  ! * netcdf; a FORTRAN netcdf_var_struct variable with appropriate
  !   grid dimension variables defined and structure arrays allocated.

  !-----------------------------------------------------------------------

  subroutine nmm_netcdf_read_var_init(filename,json,netcdf)

    ! Define variables passed to routine

    type(netcdf_var_struct)                                             :: netcdf
    type(nmm_json_var_struct)                                           :: json
    character(len=500)                                                  :: filename

    ! Define variables computed within routine

    character(len=100)                                                  :: dimname

    !=====================================================================

    ! Define local variables

    dimname = 'Time'
    call netcdf_interface_getdim(filename,dimname,netcdf%nt)
    dimname = 'west_east'
    call netcdf_interface_getdim(filename,dimname,netcdf%nx)
    dimname = 'south_north'
    call netcdf_interface_getdim(filename,dimname,netcdf%ny)
    
    ! Check local variable and proceed accordingly

    if(json%levtype .eq. 1) then

       ! Define local variables

       netcdf%nz = 1

    end if ! if(json%levtype .eq. 3)

    ! Check local variable and proceed accordingly

    if(json%levtype .eq. 2) then

       ! Define local variables

       dimname = 'soil_layers_stag'
       call netcdf_interface_getdim(filename,dimname,netcdf%nz)

    end if ! if(json%levtype .eq. 2)

    ! Check local variable and proceed accordingly

    if(json%levtype .eq. 3) then

       ! Check local variable and proceed accordingly

       if(json%z_staggered) then

          ! Define local variables

          dimname = 'bottom_top_stag'
          
       else   ! if(json%z_staggered) 

          ! Define local variables

          dimname = 'bottom_top'

       end if ! if(json%z_staggered) 

       ! Define local variables

       call netcdf_interface_getdim(filename,dimname,netcdf%nz)

    end if ! if(json%levtype .eq. 3)

    ! Define local variables

    call netcdf_interface_close()

    !=====================================================================

  end subroutine nmm_netcdf_read_var_init

  !=======================================================================

  ! SUBROUTINE:

  ! tcv_read.f90

  ! DESCRIPTION:

  ! This subroutine ingests an external file containing tropical
  ! cyclone attributes (e.g., TC-vitals).

  ! INPUT VARIABLES:

  ! * filename; a FORTRAN character string specifying the path to the
  !   TCV file.

  ! * tcv; a FORTRAN tcv_struct variable.

  ! OUTPUT VARIABLES:

  ! * tcv; a FORTRAN tcv_struct variable now containing the tropical
  !   cyclone attributes retrieved from the user specified file.

  !-----------------------------------------------------------------------

  subroutine tcv_read(filename,tcv)

    ! Define variables passed to routine

    type(tcv_struct),           dimension(:),               allocatable :: tcv
    character(len=500)                                                  :: filename

    ! Define variables computed within routine

    character(len=1)                                                    :: dummy
    real(r_kind)                                                        :: lat_scale
    real(r_kind)                                                        :: lon_scale

    ! Define counting variables

    integer                                                             :: i

    !=====================================================================

    ! Define local variables

    ntcs = 0
    open(99,file=trim(adjustl(filename)),form='formatted')
1000 read(99,*,end=1001) dummy
    ntcs = ntcs + 1
    goto 1000
1001 continue
    close(99)

    ! Allocate memory for local variables

    if(.not. allocated(tcv)) allocate(tcv(ntcs))

    ! Define local variables

    open(99,file=trim(adjustl(filename)),form='formatted')

    ! Loop through local variable

    do i = 1, ntcs

       ! Define local variables

       read(99,500) tcv(i)%center, tcv(i)%id, tcv(i)%name, tcv(i)%century, &
            & tcv(i)%yymmdd, tcv(i)%hhmm, tcv(i)%lati, tcv(i)%latns,       &
            & tcv(i)%loni, tcv(i)%lonew, tcv(i)%stdir, tcv(i)%stspd,       &
            & tcv(i)%pcen, tcv(i)%penv, tcv(i)%penvrad, tcv(i)%vmax,       &
            & tcv(i)%vmaxrad, tcv(i)%r15ne, tcv(i)%r15se, tcv(i)%r15sw,    &
            & tcv(i)%r15nw, tcv(i)%depth
       lat_scale = 1.0/10.0
       lon_scale = 1.0/10.0

       ! Check local variable and proceed accordingly

       if(tcv(i)%latns .eq. 'S') lat_scale = -1.0*lat_scale
       if(tcv(i)%lonew .eq. 'W') lon_scale = -1.0*lon_scale

       ! Compute local variables

       tcv(i)%lat = real(tcv(i)%lati)*lat_scale
       tcv(i)%lon = real(tcv(i)%loni)*lon_scale

    end do ! do i = 1, ntcs

    ! Define local variables

    close(99)
500 format(a4,1x,a3,1x,a9,1x,i2,i6,1x,i4,1x,i3,a1,1x,i4,a1,1x,i3,1x,i3,    &
         & 3(1x,i4),1x,i2,1x,i3,1x,4(i4,1x),a1)

    !=====================================================================

  end subroutine tcv_read

  !=======================================================================

end module fileio_interface
