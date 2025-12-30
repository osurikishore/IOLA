module interpolation_interface

  !=======================================================================

  !$$$ MODULE DOCUMENTATION BLOCK
  
  ! da-update-analysis :: interpolation_interface
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
  use math_methods_interface
  use namelist_interface
  use slint
  use variable_interface

  ! Define interfaces and attributes for module routines
  
  implicit none
  private
  public :: interpolation_interface_bilinear
  public :: interpolation_interface_cleanup
  public :: interpolation_interface_init
  public :: interpolation_interface_nest 
  public :: interpolation_interface_nnghbr
  public :: interpolation_interface_remap
  public :: interpolation_interface_vertical
  interface interpolation_interface_bilinear
     module procedure recenter_bilinear
     module procedure slint_bilinear
  end interface interpolation_interface_bilinear
  interface interpolation_interface_cleanup
     module procedure slint_cleanup
  end interface interpolation_interface_cleanup
  interface interpolation_interface_init
     module procedure init_recenter
     module procedure init_slint
  end interface interpolation_interface_init
  interface interpolation_interface_nest
     module procedure slint_nest
  end interface interpolation_interface_nest
  interface interpolation_interface_nnghbr
     module procedure slint_nnghbr
  end interface interpolation_interface_nnghbr
  interface interpolation_interface_vertical
     module procedure hwrf_vertical
     module procedure presremap_vertical
  end interface interpolation_interface_vertical

  !-----------------------------------------------------------------------

contains

  !=======================================================================

  ! SUBROUTINE:

  ! hwrf_vertical.f90

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

  subroutine hwrf_vertical(json,hwrf_vert_grid)

    ! Define variables passed to routine

    type(nmm_json_var_struct)                                           :: json
    type(hwrf_vert_grid_struct)                                         :: hwrf_vert_grid

    ! Define variables computed within routine

    type(vertgrid_struct)                                               :: vertgrid

    ! Define counting variables

    integer                                                             :: i

    !=====================================================================

    ! Define local variables

    vertgrid%nz = hwrf_vert_grid%nz
    if(json%z_staggered) vertgrid%nz = vertgrid%nz + 1
    call variable_interface_setup_struct(vertgrid)
    
    ! Loop through local variable

    do i = 1, hwrf_vert_grid%ncoords

       ! Check local variable and proceed accordingly

       if(json%z_staggered) then
          
          ! Define local variables

          vertgrid%dst_prs = hwrf_vert_grid%pint(i,:)
          vertgrid%src_prs = hwrf_vert_grid%nest_pint(i,:)

       else   ! if(json%z_staggered)

          ! Define local variables
          
          vertgrid%dst_prs = hwrf_vert_grid%plev(i,:)
          vertgrid%src_prs = hwrf_vert_grid%nest_plev(i,:)

       end if ! if(json%z_staggered)

       ! Define local variables
       
       vertgrid%src_var(1:vertgrid%nz) =                                   &
            & hwrf_vert_grid%var(i,1:vertgrid%nz)
       
       ! Compute local variables

       if(is_llp) call linear_log_pressure(vertgrid)

       ! Define local variables

       hwrf_vert_grid%var(i,1:vertgrid%nz) =                               &
            & vertgrid%dst_var(1:vertgrid%nz)

    end do ! do i = 1, hwrf_vert_grid%ncoords
       
    ! Deallocate memory for local variables
    
    call variable_interface_cleanup_struct(vertgrid)

    !=====================================================================

  end subroutine hwrf_vertical

  !=======================================================================

  ! SUBROUTINE:

  ! init_recenter.f90

  ! DESCRIPTION:

  ! This subroutine initializes and defines the grid variable used to
  ! remap a variable field from a Cartesian (e.g., longitude/latitude
  ! grid) to a coordinate system recentered relative to a user
  ! specified location.

  ! INPUT VARIABLES:

  ! * remap_grid; a FORTRAN recenter_struct variable.

  ! OUTPUT VARIABLES:

  ! * remap_grid; a FORTRAN recenter_struct variable containing the
  !   radii and angles for the coordinate projection recentered
  !   relative to the user specified position.

  !-----------------------------------------------------------------------  

  subroutine init_recenter(remap_grid)

    ! Define variables passed to routine

    type(recenter_struct)                                               :: remap_grid

    ! Define variables computed within routine

    type(grid_struct)                                                   :: grid
    real(r_kind)                                                        :: latmin
    real(r_kind)                                                        :: latmax
    real(r_kind)                                                        :: lonmin
    real(r_kind)                                                        :: lonmax

    ! Define counting variables

    integer                                                             :: i, j

    !=====================================================================

    ! Define local variables

    latmin             = remap_grid%clat - remap_grid%area
    latmax             = remap_grid%clat + remap_grid%area
    lonmin             = remap_grid%clon - remap_grid%area
    lonmax             = remap_grid%clon + remap_grid%area
    remap_grid%nx      = ((lonmax - lonmin)/remap_grid%darea) + 1
    remap_grid%ny      = ((latmax - latmin)/remap_grid%darea) + 1
    call variable_interface_setup_struct(remap_grid)
    remap_grid%ncoords = 0

    ! Loop through local variable

    do j = 1, remap_grid%ny

       ! Loop through local variable

       do i = 1, remap_grid%nx
          
          ! Define local variables

          remap_grid%ncoords                 = remap_grid%ncoords + 1
          remap_grid%lon(remap_grid%ncoords) = (lonmin +                   &
               & (i-1)*remap_grid%darea)
          remap_grid%lat(remap_grid%ncoords) = (latmin +                   &
               & (j-1)*remap_grid%darea)

       end do ! do i = 1, remap_grid%nx

    end do ! do j = 1, remap_grid%ny

    ! Define local variables

    grid%nx      = remap_grid%nx
    grid%ny      = remap_grid%ny
    grid%ncoords = (grid%nx*grid%ny)
    call variable_interface_setup_struct(grid)
    grid%lon     = remap_grid%lon
    grid%lat     = remap_grid%lat

    ! Compute local variables

    call math_methods_polarcoords(grid,remap_grid%clon,remap_grid%clat)

    ! Define local variables

    remap_grid%radius = grid%radius
    remap_grid%angle  = grid%angle

    ! Deallocate memory for local variables

    call variable_interface_cleanup_struct(grid)

    !=====================================================================

  end subroutine init_recenter

  !=======================================================================

  ! SUBROUTINE: 

  ! init_slint.f90

  ! DESCRIPTION: 

  ! This subroutine initializes and computes the remapping
  ! coefficients and nearest neighbor locations using the geographical
  ! locations passed using the source (src) and destination (dst)
  ! grid_struct FORTRAN structures; all geographical locations are
  ! assumed to be units of degrees.

  ! INPUT VARIABLES:

  ! * src_grid; a FORTRAN grid_struct variable containing the source
  !   grid information; all geographical location units are degrees.

  ! * dst_grid; a FORTRAN grid_struct variable containing the
  !   destination grid information; all geographical location units
  !   are degrees.

  ! * remap_grid; a FORTRAN slint_struct variable.

  ! OUTPUT VARIABLES:

  ! * remap_grid; a FORTRAN slint_struct containing the remapping
  !   coefficients and nearest neighbor indices determined from the
  !   source and destination grid_struct FORTRAN structures.

  !-----------------------------------------------------------------------

  subroutine init_slint(src_grid,dst_grid,remap_grid)

    ! Define variables passed to routine

    type(grid_struct)                                                   :: src_grid
    type(grid_struct)                                                   :: dst_grid
    type(slint_struct)                                                  :: remap_grid

    ! Define variables computed within routine

    type(grid)                                                          :: slint_grid
    real(r_kind),               dimension(:,:),             allocatable :: grid1
    real(r_kind),               dimension(:,:),             allocatable :: grid2

    !=====================================================================

    ! Allocate memory for local variables

    if(.not. allocated(grid1))                                             &
         & allocate(grid1(src_grid%ncoords,2))
    if(.not. allocated(grid2))                                             &
         & allocate(grid2(dst_grid%ncoords,2))

    ! Define local variables

    grid1(:,1) = src_grid%lat*deg2rad
    grid1(:,2) = src_grid%lon*deg2rad
    grid2(:,1) = dst_grid%lat*deg2rad
    grid2(:,2) = dst_grid%lon*deg2rad

    ! Compute local variables
    
    call slint_init_compute(grid1,src_grid%ncoords,grid2,                  &
         & dst_grid%ncoords,slint_grid)

    ! Define local variables

    remap_grid%coeffs      = slint_grid%coeffs
    remap_grid%nn          = slint_grid%nn
    remap_grid%src_ncoords = src_grid%ncoords

    ! Deallocate memory for local variables

    if(allocated(grid2)) deallocate(grid2)
    if(allocated(grid1)) deallocate(grid1)

    !=====================================================================

  end subroutine init_slint

  !=======================================================================

  ! SUBROUTINE:

  ! interpolation_interface_remap.f90

  ! DESCRIPTION:

  ! This subroutine remaps a variable grid, which as been relative to
  ! a geographical location to another projection; this is namely the
  ! case when a variable field has been analyzed in a space other than
  ! real (or Cartesian) space (e.g., spectral space or polar
  ! coordinates).

  ! INPUT VARIABLES:

  ! * src_grid; a FORTRAN grid_struct variable containing the source
  !   variable projection (e.g., geographical locations) and variable
  !   grid to be remapped.

  ! * remap_grid; a FORTRAN recenter_struct variable containing the
  !   geographical locations for remapping.

  ! OUTPUT VARIABLES:

  ! * remap_grid; a FORTRAN recenter_struct variable containing the
  !   remapped variable values (within the remap_grid variable 'var'
  !   attribute).

  !-----------------------------------------------------------------------

  subroutine interpolation_interface_remap(src_grid,remap_grid)

    ! Define variables passed to routine

    type(grid_struct)                                                   :: src_grid
    type(recenter_struct)                                               :: remap_grid

    ! Define variables computed within routine

    type(grid_struct)                                                   :: dst_grid
    type(slint_struct)                                                  :: slint

    ! Define counting variables

    integer                                                             :: i

    !=====================================================================

    ! Define local variables

    dst_grid%ncoords = (remap_grid%nx*remap_grid%ny)
    dst_grid%nz      = src_grid%nz
    call variable_interface_setup_struct(dst_grid)
    dst_grid%lat     = remap_grid%lat
    dst_grid%lon     = remap_grid%lon

    ! Check local variable and proceed accordingly

    if(is_slint) then

       ! Define local variables

       slint%ncoords = dst_grid%ncoords
       call variable_interface_setup_struct(slint)
       call interpolation_interface_init(src_grid,dst_grid,slint)

       ! Loop through local variable

       do i = 1, remap_grid%nz

          ! Compute local variables

          call interpolation_interface_bilinear(slint,src_grid%var(:,i))

          ! Define local variables

          remap_grid%var(:,i) = slint%var

       end do ! do i = 1, remap_grid%nz

       ! Deallocate memory for local variables

       call variable_interface_cleanup_struct(slint)

    end if ! if(is_slint)

    ! Deallocate memory for local variables

    call variable_interface_cleanup_struct(dst_grid)

    !=====================================================================

  end subroutine interpolation_interface_remap

  !=======================================================================

  ! SUBROUTINE:

  ! interpolation_interface_spline.f90

  ! DESCRIPTION:

  ! This subroutine interpolates, using cubic splines, to find the
  ! value of a variable at a location specified by the user; this
  ! subroutine implements the sort2 subroutine of SLATEC and the
  ! spline_cubic_set and spline_cubic_val subroutines of SPLINE.

  ! REFERENCES:

  ! de Boor, C., 1978: A practical guide to splines. Springer-Verlag,
  ! 346pp.

  ! https://people.sc.fsu.edu/~jburkardt/f_src/spline/spline.html

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN spline_struct containing the location (xa) and
  !   variable (ya) arrays as well as the interpolation location (x).

  ! OUTPUT VARIABLES:

  ! * grid; a FORTRAN interp_spline_struct containing the interpolated
  !   value (y).

  !-----------------------------------------------------------------------

  subroutine interpolation_interface_spline(grid)

    ! Define variables passed to routine

    type(spline_struct)                                                 :: grid

    ! Define variables computed within routine

    type(spline_struct)                                                 :: gridl
    real(r_double),             dimension(:),               allocatable :: xa
    real(r_double),             dimension(:),               allocatable :: ya
    real(r_double),             dimension(:),               allocatable :: y2a
    real(r_double)                                                      :: x
    real(r_double)                                                      :: y
    real(r_double)                                                      :: yp
    real(r_double)                                                      :: ypp
    real(r_double)                                                      :: yp1
    real(r_double)                                                      :: ypn
    integer                                                             :: n

    ! Define counting variables

    integer                                                             :: i, j

    !=====================================================================

    ! Define local variables

    n       = count(grid%xa .ne. spval)
    x       = dble(grid%x)
    y       = spval
    gridl%n = n
    call variable_interface_setup_struct(gridl)

    ! Check local variable and proceed accordingly

    if(n .ge. 2) then

       ! Allocate memory for local variables
       
       if(.not. allocated(xa))  allocate(xa(n))
       if(.not. allocated(ya))  allocate(ya(n))
       if(.not. allocated(y2a)) allocate(y2a(n))

       ! Define local variables

       j = 0

       ! Loop through local variable

       do i = 1, size(grid%xa)

          ! Check local variable and proceed accordingly

          if(grid%xa(i) .ne. spval .and. grid%ya(i) .ne. spval) then

             ! Define local variables
             
             j     = j + 1
             xa(j) = dble(grid%xa(i))
             ya(j) = dble(grid%ya(i))

          end if ! if(grid%xa(i) .ne. spval .and. grid%ya(i)
                 ! .ne. spval)

       end do ! do i = 1, size(grid%xa)

       ! Define local variables

       gridl%xa = real(xa)
       gridl%ya = real(ya)
       call math_methods_unique_array(gridl)
       call math_methods_sort_array(gridl,.true.,.false.)
       n        = gridl%n

       ! Deallocate memory for local variables

       if(allocated(xa))  deallocate(xa)
       if(allocated(ya))  deallocate(ya)       
       if(allocated(y2a)) deallocate(y2a)

       ! Allocate memory for local variables
       
       if(.not. allocated(xa))  allocate(xa(n))
       if(.not. allocated(ya))  allocate(ya(n))
       if(.not. allocated(y2a)) allocate(y2a(n))
       
       ! Define local variables

       xa(1:n) = dble(gridl%xa(1:n))
       ya(1:n) = dble(gridl%ya(1:n))
       yp1     = dble(0.0)
       ypn     = dble(0.0)

       ! Compute local variables

       call spline_cubic_set(n,xa(1:n),ya(1:n),3,yp1,3,ypn,y2a)

       ! Check local variable and proceed accordingly

       if(minval(xa(1:n)) .le. x .and. maxval(xa(1:n)) .ge. x) then

          ! Compute local variables

          call spline_cubic_val(n,xa(1:n),ya(1:n),y2a(1:n),x,y,yp,ypp)

       end if ! if(minval(xa(1:n)) .le. x .and. maxval(xa(1:n))
              ! .ge. x)

       ! Deallocate memory for local variables

       if(allocated(xa))  deallocate(xa)
       if(allocated(ya))  deallocate(ya)       
       if(allocated(y2a)) deallocate(y2a)

    end if ! if(n .ge. 2)

    ! Deallocate memory for local variables

    call variable_interface_cleanup_struct(gridl)

    ! Define local variables

    grid%y = real(y)

    !=====================================================================

  end subroutine interpolation_interface_spline

  !=======================================================================

  ! SUBROUTINE:

  ! linear_log_pressure.f90

  ! DESCRIPTION:

  ! This subroutine interpolates from a source grid vertical pressure
  ! coordinate to a destination grid vertical pressure coordinate
  ! using linear-log pressure (llp) interpolation.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN vertgrid_struct variable containing the source
  !   and destination pressure profiles as well as the source grid
  !   variable to be interpolated.

  ! OPTIONAL INPUT VARIABLES:

  ! * flag; a FORTRAN logical variable specifying whether to flag
  !   interpolated values (with spval) if the level of interpolation
  !   is above or below the source profile pressure level; default
  !   value is .false..

  ! OUTPUT VARIABLES:

  ! * grid; a FORTRAN vertgrid_struct variable where the 'var'
  !   variable has been updated in accordance with the interpolation !
  !   to destination grid pressure profile.

  !-----------------------------------------------------------------------

  subroutine linear_log_pressure(grid,flag)

    ! Define variables passed to routine

    type(vertgrid_struct)                                               :: grid
    logical, optional, intent(in)                                       :: flag

    ! Define variables computed within routine

    logical                                                             :: flagn
    integer                                                             :: idxb
    integer                                                             :: idxt

    ! Define counting variables

    integer                                                             :: i, j

    !=====================================================================

    ! Check local variable and proceed accordingly

    if(.not. present(flag)) then
       
       ! Define local variables
       
       flagn = .false.

    else   ! if(.not. present(flag))

       ! Define local variables

       flagn = flag

    end if ! if(.not. present(flag))

    ! Loop through local variable

    do j = 1, grid%nz

       ! Define local variables

       idxb = 0
       idxt = 0

       ! Loop through local variable

       do i = 2, grid%nz

          ! Check local variable and proceed accordingly

          if((grid%dst_prs(j) .le. grid%src_prs(i-1)) .and.                &
               & (grid%dst_prs(j) .gt. grid%src_prs(i))) then

             ! Define local variables

             idxb = i - 1
             idxt = i
             goto 1000

          end if ! if((grid%dst_prs(j) .le. grid%src_prs(i-1))
                 ! .and. (grid%dst_prs(j) .gt. grid%src_prs(i)))

       end do ! do i = 2, grid%nz

       ! Define local variables

1000   continue

       ! Check local variable and proceed accordingly

       if(flagn) then

          ! Check local variable and proceed accordingly

          if((idxb .eq. 0) .or. (idxt .eq. 0)) then

             ! Define local variables

             grid%dst_var(j) = spval

          end if ! if((idxb .eq. 0) .or. (idxt .eq. 0))

       end if ! if(flagn)

       ! Check local variable and proceed accordingly

       if(.not. flagn) then

          ! Check local variable and proceed accordingly

          if((idxb .eq. 0) .or. (idxt .eq. 0)) then

             ! Define local variables

             idxb = 1
             idxt = grid%nz
          
          end if ! if((idxb .eq. 0) .or. (idxt .eq. 0))

          ! Compute local variables

          grid%dst_var(j) = grid%src_var(idxb) + (grid%src_var(idxt) -     &
               & grid%src_var(idxb))*((log(grid%dst_prs(j)) -              &
               & log(grid%src_prs(idxb)))/(log(grid%src_prs(idxt)) -       &
               & log(grid%src_prs(idxb))))

       end if ! if(.not. flagn)

    end do !  do j = 1, grid%nz

    !=====================================================================

  end subroutine linear_log_pressure

  !=======================================================================

  ! SUBROUTINE:

  ! presremap_vertical.f90

  ! DESCRIPTION:

  ! This subroutine remaps a variable array from a source grid
  ! pressure coordinate to a destination grid pressure coordinate
  ! using (optional, logarithmic) cubic-splines.

  ! INPUT VARIABLES:

  ! * presremap; a FORTRAN presremap_struct variable containing the
  !   source variable pressure coordinate and variable values
  !   (src_pres and src_var attributes, respectively) and the
  !   destination pressure coordinate (dst_pres).

  ! OUTPUT VARIABLES:

  ! * presremap; a FORTRAN presremap_struct variable containing the
  !   source grid variable remapped to the destination pressure
  !   coordinate values.

  !-----------------------------------------------------------------------

  subroutine presremap_vertical(presremap)

    ! Define variables passed to routine

    type(presremap_struct)                                              :: presremap

    ! Define variables computed within routine

    type(spline_struct)                                                 :: spline

    ! Define counting variables

    integer                                                             :: i, j

    !=====================================================================

    ! Define local variables
    
    spline%n            = presremap%src_nz
    call variable_interface_setup_struct(spline)
    presremap%remap_var = presremap%dst_var

    ! Loop through local variable

    do i = 1, presremap%src_ncoords

       ! Define local variables

       spline%xa = presremap%src_pres(i,:)
       spline%ya = presremap%src_var(i,:)

       ! Check local variable and proceed accordingly

       if(is_llp) spline%xa = log(spline%xa)

       ! Loop through local variable

       do j = 1, presremap%dst_nz

          ! Define local variables

          spline%x = presremap%dst_pres(i,j)

          ! Check local variable and proceed accordingly

          if(is_llp) spline%x = log(spline%x)

          ! Compute local variables

          call interpolation_interface_spline(spline)

          ! Check local variable and proceed accordingly

          if(spline%y .ne. spval) then

             ! Define local variables

             presremap%remap_var(i,j) = spline%y

          end if ! if(spline%y .eq. spval)

       end do ! do j = 1, presremap%dst_ncoords

    end do ! do i = 1, presremap%src_ncoords

    ! Deallocate memory for local variables
    
    call variable_interface_cleanup_struct(spline)

    !=====================================================================

  end subroutine presremap_vertical

  !=======================================================================

  ! SUBROUTINE:

  ! recenter_bilinear.f90

  ! DESCRIPTION:

  ! This subroutine interpolates a Cartesian grid to grid where the
  ! center is defined by the reference longitude and latitude
  ! coordinate specified by the user; the original grid (i.e., grid)
  ! may be of any allowed projection; the recentered grid (i.e.,
  ! remap_grid) is defined by computing a longitude and latitude grid,
  ! centered at the user specified reference coordinates, and
  ! remapping via bi-linear interpolation to original (e.g.,
  ! Cartesian) grid to the new grid.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN grid_struct variable containing the Cartesian
  !   grid longitude and latitude values.

  ! * var; a FORTRAN 4-byte real array containing the variable values
  !   to be interpolated.

  ! * remap_grid; a FORTRAN recenter_struct variable containing the
  !   reference longitude and latitude coordinates for the new grid
  !   remapping.

  ! OUTPUT VARIABLES:

  ! * remap_grid; a FORTRAN recenter_struct variable containing the
  !   user specified variable field remapped to the grid projection
  !   centered at the reference location specified by the user.

  !-----------------------------------------------------------------------

  subroutine recenter_bilinear(grid,remap_grid)

    ! Define variables passed to routine

    type(grid_struct)                                                   :: grid
    type(recenter_struct)                                               :: remap_grid

    ! Define variables computed within routine

    type(slint_struct)                                                  :: slint
    type(grid_struct)                                                   :: src_grid
    type(grid_struct)                                                   :: dst_grid

    ! Define counting variables

    integer                                                             :: i

    !=====================================================================
    
    ! Define local variables

    call interpolation_interface_init(remap_grid)
    src_grid%ncoords = (grid%nx*grid%ny)
    src_grid%nz      = grid%nz
    dst_grid%ncoords = (remap_grid%nx*remap_grid%ny)
    dst_grid%nz      = grid%nz
    call variable_interface_setup_struct(dst_grid)
    call variable_interface_setup_struct(src_grid)
    src_grid%lon     = grid%lon
    src_grid%lat     = grid%lat
    dst_grid%lon     = remap_grid%lon
    dst_grid%lat     = remap_grid%lat
    src_grid%var     = grid%var

    ! Check local variable and proceed accordingly

    if(is_slint) then

       ! Define local variables
       
       slint%ncoords = (remap_grid%nx*remap_grid%ny)
       call variable_interface_setup_struct(slint)
       call interpolation_interface_init(src_grid,dst_grid,slint)

       ! Loop through local variable

       do i = 1, grid%nz

          ! Compute local variables

          call interpolation_interface_bilinear(slint,src_grid%var(:,i))

          ! Define local variables

          remap_grid%var(:,i) = slint%var

       end do ! do i = 1, grid%nz

    end if ! if(is_slint)

    ! Deallocate memory for local variables

    call variable_interface_cleanup_struct(src_grid)
    call variable_interface_cleanup_struct(dst_grid)
    call variable_interface_cleanup_struct(slint)

    !=====================================================================

  end subroutine recenter_bilinear

  !=======================================================================

  ! SUBROUTINE: 

  ! slint_bilinear.f90

  ! DESCRIPTION: 

  ! This subroutine interpolates a variable to a destination grid
  ! using bi-linear interpolation means via the SLINT algorithm.

  ! INPUT VARIABLES:

  ! * remap_grid; a FORTRAN slint_struct variable containing the
  !   remapping (from source to destination) coefficients and nearest
  !   neighbors.

  ! * var; a FORTRAN 4-byte real value containing the gridded variable
  !   on the source grid.

  ! OUTPUT VARIABLES:

  ! * remap_grid; a FORTRAN slint_struct variable containing the
  !   gridded source variable interpolated to the destination grid.

  !-----------------------------------------------------------------------

  subroutine slint_bilinear(remap_grid,var)

    ! Define variables passed to routine

    type(slint_struct)                                                  :: remap_grid
    real(r_kind)                                                        :: var(remap_grid%src_ncoords)

    ! Define variables computed within routine

    real(r_double)                                                      :: c(3)
    real(r_double)                                                      :: v(3)

    ! Define counting variables

    integer                                                             :: i

    !=====================================================================

    ! Define local variables

    remap_grid%var = 0.0
    
    ! Check local variable and proceed accordingly

    if(minval(var) .ne. spval .and. maxval(var) .ne. spval) then

       ! Loop through local variable

       do i = 1, remap_grid%ncoords

          ! Define local variables

          c    = remap_grid%coeffs(:,i)
          v(1) = var(remap_grid%nn(1,i))
          v(2) = var(remap_grid%nn(2,i))
          v(3) = var(remap_grid%nn(3,i))

          ! Compute local variables

          remap_grid%var(i) = c(1)*v(1) + c(2)*v(2) + c(3)*v(3)
             
       end do ! do i = 1, remap_grid%ncoords

    else   ! if(minval(var) .ne. spval .and. maxval(var) .ne. spval)

       ! Define local variables

       remap_grid%var = spval

    end if ! if(minval(var) .ne. spval .and. maxval(var) .ne. spval)

    !=====================================================================

  end subroutine slint_bilinear

  !=======================================================================

  ! SUBROUTINE:

  ! slint_cleanup.f90

  ! DESCRIPTION:

  ! This subroutine deallocates all memory for all SLINT algorithm
  ! interpolation structures.

  ! INPUT VARIABLES:

  ! * slint; a FORTRAN slint_struct variable, of dimension 3,
  !   containing the mass and momentum (both zonal and meridional wind
  !   component grids) remapping grid attributes.

  !-----------------------------------------------------------------------

  subroutine slint_cleanup(slint)

    ! Define variables passed to routine

    type(slint_struct)                                                  :: slint(3)

    ! Define counting variables

    integer                                                             :: i

    !=====================================================================

    ! Loop through local variable

    do i = 1, 3

       ! Deallocate memory for local variables
       
       call variable_interface_cleanup_struct(slint(i))

    end do ! do i = 1, 3

    !=====================================================================

  end subroutine slint_cleanup

  !=======================================================================

  ! SUBROUTINE:

  ! slint_nest.f90

  ! DESCRIPTION:
 
  ! This subroutine interpolates a nested grid variable to a parent
  ! grid using the SLINT algorithm in accordance with the user
  ! specified interpolation type.

  ! INPUT VARIABLES:

  ! * grid_nest; a FORTRAN grid_stuct variable containing the nest
  !   grid attributes.

  ! * grid_parent; a FORTRAN grid_struct variable containing the
  !   parent grid attributes.

  ! * remap_grid; a FORTRAN slint_struct variable containing the
  !   remapping (from source to destination) coefficients and nearest
  !   neighbors.

  ! * intrptype; a FORTRAN integer valued variable indicating the
  !   interpolation/remapping type; the current options are:

  ! 1: SLINT nearest-neighbor interpolation.

  ! 2: SLINT bilinear interpolation.

  ! OUTPUT VARIABLES:

  ! * remap_grid; a FORTRAN slint_struct variable containing the
  !   gridded source variable interpolated to the destination grid. 

  !-----------------------------------------------------------------------

  subroutine slint_nest(grid_nest,grid_parent,remap_grid,intrptype)

    ! Define variables passed to routine

    type(grid_struct)                                                   :: grid_nest
    type(grid_struct)                                                   :: grid_parent
    type(slint_struct)                                                  :: remap_grid
    integer                                                             :: intrptype

    ! Define variables computed within routine

    type(grid_struct)                                                   :: dst_grid
    type(grid_struct)                                                   :: src_grid
    integer                                                             :: dst_idx

    ! Define counting variables

    integer                                                             :: i, j

    !=====================================================================

    ! Define local variables

    dst_grid%ncoords = grid_parent%ncoords
    dst_grid%nz      = grid_parent%nz
    call variable_interface_setup_struct(dst_grid)
    dst_grid%var     = spval
    src_grid%ncoords = grid_nest%ncoords
    src_grid%nz      = grid_nest%nz
    call variable_interface_setup_struct(src_grid)
    src_grid%var     = grid_nest%var

    ! Loop through local variable

    do i = 1, dst_grid%nz

       ! Check local variable and proceed accordingly

       if(intrptype .eq. 1) then

          ! Compute local variables
          
          call slint_nnghbr(remap_grid,src_grid%var(:,i))

       end if ! if(intrptype .eq. 1)

       ! Check local variable and proceed accordingly

       if(intrptype .eq. 2) then

          ! Compute local variables
          
          call slint_bilinear(remap_grid,src_grid%var(:,i))

       end if ! if(intrptype .eq. 2)

       ! Loop through local variable

       do j = 1, remap_grid%ncoords
       
          ! Define local variables

          dst_idx                 = remap_grid%idx(j)
          dst_grid%var(dst_idx,i) = remap_grid%var(j)

       end do ! do j = 1, remap_grid%ncoords
          
       ! Define local variables

       where(dst_grid%var(:,i) .eq. spval) dst_grid%var(:,i) =             &
            & grid_parent%var(:,i)
       grid_parent%var(:,i) = dst_grid%var(:,i)

    end do ! do i = 1, dst_grid%nz

    ! Deallocate memory for local variables

    call variable_interface_cleanup_struct(dst_grid)
    call variable_interface_cleanup_struct(src_grid)

    !=====================================================================

  end subroutine slint_nest

  !=======================================================================

  ! SUBROUTINE: 
  
  ! slint_nnghbr.f90

  ! DESCRIPTION:

  ! This subroutine interpolates a variable to a destination grid
  ! using nearest-neighbor interpolation means via the SLINT
  ! algorithm.

  ! INPUT VARIABLES:

  ! * remap_grid; a FORTRAN slint_struct variable containing the
  !   remapping (from source to destination) coefficients and nearest
  !   neighbors.

  ! * var; a FORTRAN 4-byte real value containing the gridded variable
  !   on the source grid.

  ! OUTPUT VARIABLES:

  ! * remap_grid; a FORTRAN slint_struct variable containing the
  !   gridded source variable interpolated to the destination grid.

  !-----------------------------------------------------------------------

  subroutine slint_nnghbr(remap_grid,var)

    ! Define variables passed to routine

    type(slint_struct)                                                  :: remap_grid
    real(r_kind)                                                        :: var(remap_grid%src_ncoords)

    ! Define variables computed within routine

    real(r_double)                                                      :: v(3)

    ! Define counting variables

    integer                                                             :: i

    !=====================================================================

    ! Define local variables

    remap_grid%var = 0.0

    ! Loop through local variable

    do i = 1, remap_grid%ncoords
       
       ! Define local variables
       
       v(1)              = var(remap_grid%nn(1,i))
       remap_grid%var(i) = v(1)
          
    end do ! do i = 1, remap_grid%ncoords

    !=====================================================================

  end subroutine slint_nnghbr

  !=======================================================================

end module interpolation_interface
