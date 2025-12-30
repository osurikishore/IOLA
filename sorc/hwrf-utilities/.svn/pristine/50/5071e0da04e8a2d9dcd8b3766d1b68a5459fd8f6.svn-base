module grid_methods_interface

  !=======================================================================

  !$$$ MODULE DOCUMENTATION BLOCK
  
  ! da-update-analysis :: grid_methods_interface
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
  use kinds_interface
  use namelist_interface
  use variable_interface

  ! Define interfaces and attributes for module routines
  
  implicit none
  private
  public :: grid_methods_gcdist
  public :: grid_methods_gcgeo

  !-----------------------------------------------------------------------

contains

  !=======================================================================
  
  ! SUBROUTINE:

  ! grid_methods_gcdist.f90

  ! DESCRIPTION:

  ! This subroutine estimates the respective grid grid-cell resolution
  ! using the grid spacing (specified by the user) and the estimated
  ! longitude and latitude grid spacings (assuming a spherical Earth)
  ! at the equator).

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN grid_struct variable.

  ! * dx; a FORTRAN 4-byte real value specifying the zonal (e.g., x-)
  !   coordinate resolution; units are degrees.

  ! * dy; a FORTRAN 4-byte real value specifying the meridional (e.g.,
  !   y-) coordinate resolution; units are degrees.

  ! OUTPUT VARIABLES:

  ! * grid; a FORTRAN grid_struct variable where the attributes dx and
  !   dy have been define accordingly.

  !-----------------------------------------------------------------------

  subroutine grid_methods_gcdist(grid,dx,dy)

    ! Define variables passed to routine

    type(grid_struct)                                                   :: grid
    real(r_kind)                                                        :: dx
    real(r_kind)                                                        :: dy

    ! Define counting variables

    integer                                                             :: i

    !=====================================================================

    ! Define local variables

    call constants_interface_init()
    grid%dy = dlat_equator*dy
    
    ! Loop through local variable

    do i = 1, (grid%nx*grid%ny)

       ! Compute local variables

       grid%dx(i) = dlon_equator*cos(grid%lat(i)*deg2rad)*dx

    end do ! do i = 1, (grid%nx*grid%ny)

    !=====================================================================

  end subroutine grid_methods_gcdist

  !=======================================================================

  ! SUBROUTINE:

  ! grid_methods_gcgeo.f90

  ! DESCRIPTION:

  ! This subroutine computes the longitude and latitude (e.g.,
  ! geographical) coordinates, assuming a great-circle distance, of an
  ! advected position provided an initial geographical position
  ! (defined by the longitude and latitude coordinates) and a distance
  ! and heading.

  ! REFERENCES:

  ! http://www.edwilliams.org/avform.htm

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN grid_struct variable containing the initial
  !   geographical position (lon and lat) and the heading (head;
  !   degrees) and distance (dist; meters); all geographical position
  !   units are assumed to be degrees.

  ! OUTPUT VARIABLES:

  ! * grid; a FORTRAN grid_struct variable containing the newly
  !   computed geographical position (lon and lat).

  !-----------------------------------------------------------------------

  subroutine grid_methods_gcgeo(grid)

    ! Define variables passed to routine

    type(grid_struct)                                                   :: grid

    ! Define variables computed within routine

    real(r_kind)                                                        :: lon
    real(r_kind)                                                        :: lat
    real(r_kind)                                                        :: dlon
    real(r_kind)                                                        :: head
    real(r_kind)                                                        :: dist

    !=====================================================================

    ! Define local variables

    lat  = grid%gclat*deg2rad
    lon  = grid%gclon*deg2rad
    head = grid%gchead*deg2rad
    dist = (grid%gcdist/rearth_equator)

    ! Compute local variables

    lat  = asin(sin(lat)*cos(dist) + cos(lat)*sin(dist)*cos(head))
    dlon = atan2(sin(head)*sin(dist)*cos(lat),cos(dist) - sin(lat)*        &
         & sin(lat))
    lon  = mod(lon - dlon + pi, 2.0*pi) - pi

    ! Define local variables

    grid%gclon = lon*rad2deg
    grid%gclat = lat*rad2deg
    lon        = grid%gclon
       
    !=====================================================================

  end subroutine grid_methods_gcgeo

  !=======================================================================

end module grid_methods_interface
