module math_methods_interface

  !=======================================================================

  !$$$ MODULE DOCUMENTATION BLOCK
  
  ! da-update-analysis :: math_methods_interface
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
  use grid_methods_interface
  use iso_c_binding
  use kdtree2_module
  use kinds_interface
  use variable_interface

  ! Define interfaces and attributes for module routines
  
  implicit none
#ifdef DARWIN
  include 'fftw3.f03'
#endif
#ifdef LINUX
  include 'fftw3.f'
#endif
  private
  public :: math_methods_9ptsmthr
  public :: math_methods_anglegrid
  public :: math_methods_fftw_dft
  public :: math_methods_kdtree_nn
  public :: math_methods_kdtree_r2
  public :: math_methods_polarcoords
  public :: math_methods_unique_array
  public :: math_methods_radialdist
  public :: math_methods_radialgrid
  public :: math_methods_sort_array
  public :: math_methods_stats
  public :: math_methods_wnd
  interface math_methods_fftw_dft
     module procedure fftw_dft_c2c_2d
  end interface math_methods_fftw_dft
  interface math_methods_unique_array
     module procedure unique_array_spline
  end interface math_methods_unique_array
  interface math_methods_wnd
     module procedure wvnbrdcmp2d
  end interface math_methods_wnd

  !-----------------------------------------------------------------------

contains

  !=======================================================================

  ! SUBROUTINE:

  ! fftw_dft_r2c_2d.f90

  ! DESCRIPTION:

  ! This subroutine computes the discrete Fast Fourier transform (FFT)
  ! for a 2-dimensional real-valued variable using the Fastest Fourier
  ! Transform in the West (FFTW) version 3 software.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN fft2d_struct variable containing, at minimum,
  !   the input/output array and the appropriately defined array
  !   dimensions; it is assumed this FORTRAN structure is initialized
  !   prior to calling this subroutine.

  ! * frwd; a FORTRAN logical variable specifying whether to computed
  !   the forward FFT.

  ! * bkwd; a FORTRAN logical variable specifying whether to computed
  !   the backward FFT.

  ! OUTPUT VARIABLES:

  ! * grid; a FORTRAN fft2d_struct variable containing the output
  !   variable (frwd = .true.) or the input variable (bkwd = .true.).

  !-----------------------------------------------------------------------

  subroutine fftw_dft_c2c_2d(grid,frwd,bkwd)

    ! Define variables passed to routine

    type(fft2d_struct)                                                  :: grid
    logical                                                             :: frwd
    logical                                                             :: bkwd

    ! Define variables computed within routine

    integer(i_long)                                                     :: fftw_frwd
    integer(i_long)                                                     :: fftw_bkwd

    !=====================================================================

    ! Check local variable and proceed accordingly

    if(frwd) then

       ! Define local variables

       call dfftw_plan_dft_2d(fftw_frwd,grid%nx,grid%ny,grid%in,grid%out,   &
            & FFTW_FORWARD,FFTW_ESTIMATE)

       ! Compute local variables

       call dfftw_execute(fftw_frwd)

    end if ! if(frwd)

    ! Check local variable and proceed accordingly

    if(bkwd) then

       ! Define local variables

       call dfftw_plan_dft_2d(fftw_bkwd,grid%nx,grid%ny,grid%out,grid%in,   &
            & FFTW_BACKWARD,FFTW_ESTIMATE)

       ! Compute local variables

       call dfftw_execute(fftw_bkwd)

       ! Define local variables

       grid%in = grid%in/dble(grid%nx*grid%ny)

    end if ! if(bkwd)

    !=====================================================================

  end subroutine fftw_dft_c2c_2d

  !=======================================================================

  ! SUBROUTINE: 

  ! init_stats.f90

  ! DESCRIPTION:

  ! This subroutine initializes a statgrid_struct variable.

  ! INPUT VARIABLES:

  ! * statgrid; an uninitialized FORTRAN statgrid_struct variable.

  ! OUTPUT VARIABLES:

  ! * statgrid; an initialized FORTRAN statgrid_struct variable.

  !-----------------------------------------------------------------------

  subroutine init_stats(statgrid)

    ! Define variables passed to routine

    type(statgrid_struct)                                               :: statgrid

    !=====================================================================

    ! Define local variables

    statgrid%varmin = spval
    statgrid%varmax = spval
    statgrid%mean   = spval
    statgrid%vari   = spval
    statgrid%nvals  = 0

    !=====================================================================

  end subroutine init_stats

  !=======================================================================

  ! SUBROUTINE: 

  ! math_methods_9ptsmthr.f90

  ! DESCRIPTION:

  ! This subroutine applies a 9-point smoother to a user specified
  ! variable grid; the 9-point smoother is applied in accordance
  ! within the range of index positions (e.g., the minx, maxx, miny,
  ! maxy attributes).

  ! INPUT VARIABLES:

  ! * smthr9pt; a FORTRAN smthr9pt_struct variable containing the
  !   variable to be smooth and the index positions.

  ! OUTPUT VARIABLES:

  ! * smthr9pt; a FORTRAN smthr9pt_struct variable containing the
  !   smoothed variable.

  !-----------------------------------------------------------------------

  subroutine math_methods_9ptsmthr(smthr9pt)

    ! Define variables passed to routine

    type(smthr9pt_struct)                                               :: smthr9pt

    ! Define variables computed within routine

    real(r_kind),               dimension(:,:),             allocatable :: ingrid
    real(r_kind),               dimension(:,:),             allocatable :: outgrid
    integer                                                             :: xmin
    integer                                                             :: xmax
    integer                                                             :: ymin
    integer                                                             :: ymax

    ! Define counting variables

    integer                                                             :: i, j

    !=====================================================================

    ! Allocate memory for local variables

    if(.not. allocated(ingrid))                                            &
         & allocate(ingrid(smthr9pt%nx,smthr9pt%ny))
    if(.not. allocated(outgrid))                                           &
         & allocate(outgrid(smthr9pt%nx,smthr9pt%ny))
    
    ! Define local variables

    ingrid  = reshape(smthr9pt%var,shape(ingrid))
    outgrid = ingrid

    ! Loop through local variable
    
    do j = smthr9pt%miny, smthr9pt%maxy
          
       ! Loop throught local variable
       
       do i = smthr9pt%minx, smthr9pt%maxx

          ! Check local variable and proceed accordingly

          if((ingrid(i-1,j-1) .ne. spval) .and. (ingrid(i,j-1) .ne.        &
               & spval) .and. (ingrid(i+1,j-1) .ne. spval) .and.           &
               & (ingrid(i-1,j) .ne. spval) .and. (ingrid(i,j) .ne.        &
               & spval) .and. (ingrid(i+1,j) .ne. spval) .and.             &
               & (ingrid(i-1,j+1) .ne. spval) .and. (ingrid(i,j+1) .ne.    &
               & spval) .and. (ingrid(i+1,j+1) .ne. spval)) then
             
             ! Compute local variables

             outgrid(i,j) = ((ingrid(i-1,j-1) + ingrid(i,j-1) +            &
                  & ingrid(i+1,j-1) + ingrid(i-1,j) + ingrid(i,j) +        &
                  & ingrid(i+1,j) + ingrid(i-1,j+1) + ingrid(i,j+1) +      &
                  & ingrid(i+1,j+1))/9.0)

          end if ! if((ingrid(i-1,j-1) .ne. spval)
                 ! .and. (ingrid(i,j-1) .ne. spval)
                 ! .and. (ingrid(i+1,j-1) .ne. spval)
                 ! .and. (ingrid(i-1,j) .ne. spval) .and. (ingrid(i,j)
                 ! .ne. spval) .and. (ingrid(i+1,j) .ne. spval)
                 ! .and. (ingrid(i-1,j+1) .ne. spval)
                 ! .and. (ingrid(i,j+1) .ne. spval)
                 ! .and. (ingrid(i+1,j+1) .ne. spval))

       end do ! do i = smthr9pt%minx, smthr9pt%maxx

    end do ! do j = smthr9pt%miny, smthr9pt%maxy

    ! Define local variables

    smthr9pt%var = reshape(outgrid,shape(smthr9pt%var))

    ! Deallocate memory for local variables

    if(allocated(ingrid))  deallocate(ingrid)
    if(allocated(outgrid)) deallocate(outgrid)

    !=====================================================================

  end subroutine math_methods_9ptsmthr

  !=======================================================================

  ! SUBROUTINE: 

  ! math_methods_anglegrid.f90

  ! DESCRIPTION:

  ! This subroutine computes a grid of angle values relative to a
  ! reference Earth-bound geographical location and a grid of
  ! longitude and latitude values assuming the equitorial radius of
  ! the Earth.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN grid_struct variable containing grid longitude
  !   and latitude values; units for longitude and latitude are
  !   degrees.

  ! * reflat; a FORTRAN 4-byte real value specifying the geographical
  !   latitude coordinate of the reference location; units are
  !   degrees.

  ! * reflon; a FORTRAN 4-byte real value specifying the geographical
  !   longitude coordinate of the reference location; units are
  !   degrees.

  ! OUTPUT VARIABLES:

  ! * grid; a FORTRAN grid_struct variable containing the computed
  !   angle relative to a reference Earth-bound geographical location
  !   assuming the equitorial radius of the Earth; units are degrees.

  !-----------------------------------------------------------------------

  subroutine math_methods_anglegrid(grid,reflon,reflat)

    ! Define variables passed to routine

    type(grid_struct)                                                   :: grid
    real(r_kind)                                                        :: reflon
    real(r_kind)                                                        :: reflat

    ! Define variables computed within routine

    real(r_kind)                                                        :: ecircum
    real(r_kind)                                                        :: dlon
    real(r_kind)                                                        :: dlat
    real(r_kind)                                                        :: hyp_dist
    real(r_kind)                                                        :: opp_dist
    real(r_kind)                                                        :: adj_dist
    real(r_kind)                                                        :: sin_value
    real(r_kind)                                                        :: cos_value
    real(r_kind)                                                        :: sin_angle
    real(r_kind)                                                        :: cos_angle
    real(r_kind)                                                        :: tmp_angle    

    ! Define counting variables

    integer                                                             :: i

    !=====================================================================

    ! Compute local variables

    ecircum = 2.0*pi*(rearth_equator/1000.0)

    ! Define local variables

    grid%angle = 0.0

    ! Loop through local variable

    do i = 1, (grid%nx*grid%ny)

       ! Compute local variables

       dlon = abs(reflon - grid%lon(i))
       dlat = abs(reflat - grid%lat(i))
       call math_methods_radialdist(reflon,reflat,grid%lon(i),              &
            & grid%lat(i),hyp_dist)

       ! Check local variable and proceed accordingly

       if(dlon .eq. 0.0 .and. dlat .gt. 0.0) then

          ! Define local variables

          if(reflat .gt. grid%lat(i)) grid%angle(i) = 180.0
          if(reflat .lt. grid%lat(i)) grid%angle(i) = 0.0

       else if(dlon .gt. 0.0 .and. dlat .eq. 0.0) then
                
          ! Define local variables

          if(reflon .gt. grid%lon(i)) grid%angle(i) = 270.0
          if(reflon .lt. grid%lon(i)) grid%angle(i) = 90.0
                
       else   ! if(dlon .eq. 0.0 .and. dlat .gt. 0.0)

          ! Compute local variables

          opp_dist  = dlat/360.0*ecircum
          sin_value = opp_dist/(hyp_dist/1000.0)
          call math_methods_radialdist(reflon,reflat,grid%lon(i),          &
               & reflat,adj_dist)
          cos_value = (adj_dist/1000.0)/(hyp_dist/1000.0)

          ! Define local variables

          if(sin_value .gt. 1.0) sin_value = 0.99999
          if(cos_value .gt. 1.0) cos_value = 0.99999

          ! Compute local variables

          sin_angle = asin(sin_value)/deg2rad
          cos_angle = acos(cos_value)/deg2rad
          tmp_angle = 0.5*(sin_angle + cos_angle)

          ! Check local variable and proceed accordingly

          if(reflat .le. grid%lat(i) .and. reflon .le. grid%lon(i)) then

             ! Compute local variables

             grid%angle(i) = 90.0 - tmp_angle
                   
          else if(reflat .gt. grid%lat(i) .and. reflon .le. grid%lon(i))   &
               & then

             ! Compute local variables

             grid%angle(i) = 90.0 + tmp_angle

          else if(reflat .ge. grid%lat(i) .and. reflon .ge. grid%lon(i))   &
               & then

             ! Compute local variables

             grid%angle(i) = 270.0 - tmp_angle
             
          else if(reflat .lt. grid%lat(i) .and. reflon .ge. grid%lon(i))   &
               & then

             ! Compute local variables

             grid%angle(i) = 270.0 + tmp_angle

          end if ! if(reflat .le. grid%lat(i) .and. reflon
                 ! .le. grid%lon(i))

       end if ! if(dlon .eq. 0.0 .and. dlat .gt. 0.0)

    end do ! i = 1, (grid%nx*grid%ny)

    !=====================================================================
    
  end subroutine math_methods_anglegrid

  !=======================================================================

  ! SUBROUTINE:

  ! math_methods_kdtree_nn.f90

  ! DESCRIPTION:

  ! This subroutine implements a KD-tree search algorithm to locate
  ! the N-nearest neighbor locations for a user specified geographical
  ! location.

  ! REFERENCES:

  ! Kennel, M. B., 2004: KDTREE2: Fortran 95 and C++ software to
  ! efficiently search for near neighbors in a multi-dimensional
  ! Euclidean space. 

  ! http://arxiv.org/PScache/phvsics/pdf/0408/0408067.pdf.

  ! INPUT VARIABLES:
  
  ! * src_grid; a FORTRAN grid_struct variable containing the
  !   Euclidean space geographical locations within which to find
  !   nearest-neighbors; geographical locations (e.g., latitudes and
  !   longitudes) are assumed to have units of degrees.

  ! * dst_grid; a FORTRAN grid_struct variable containing the
  !   geographical locations for which to find the N-nearest neighbor
  !   locations; geographical locations (e.g., latitudes and
  !   longitudes) are assumed to have units of degrees.

  ! * kdtree; a FORTRAN kdtree_struct variable containing (at minimum)
  !   the number of coordinate values (the ncoords attribute should be
  !   equal to the dst_grid variable ncoords attribute) and the number
  !   of nearest-neighbors to seek (the nn attribute should be less
  !   than or equal to the src_grid variable ncoords attribute).

  ! OUTPUT VARIABLES:

  ! * kdtree; a FORTRAN kdtree_struct variable containing the
  !   N-nearest neighbor R^2 distances and src_grid variable
  !   coordinate values (r2dist and idx attributes, respectively); the
  !   units for the r2dist attribute are meters squared.

  !-----------------------------------------------------------------------

  subroutine math_methods_kdtree_nn(src_grid,dst_grid,kdtree)

    ! Define variables passed to routine

    type(grid_struct)                                                   :: dst_grid
    type(grid_struct)                                                   :: src_grid
    type(kdtree_struct)                                                 :: kdtree

    ! Define variables computed within routine

    type(kdtree2), pointer                                              :: kdtree2
    type(kdtree2_result)                                                :: sresults(kdtree%nn)
    real(r_kind),               dimension(:,:),             allocatable :: src_grdloc
    real(r_kind)                                                        :: dst_grdloc(3)

    ! Define counting variables

    integer                                                             :: i

    !=====================================================================

    ! Allocate memory for local variables

    if(.not. allocated(src_grdloc))                                        &
         & allocate(src_grdloc(3,src_grid%ncoords))

    ! Loop through local variables

    do i = 1, src_grid%ncoords

       ! Compute local variables

       src_grdloc(1,i) = rearth_equator*cos(src_grid%lat(i)*deg2rad)*      &
            & cos(src_grid%lon(i)*deg2rad)
       src_grdloc(2,i) = rearth_equator*cos(src_grid%lat(i)*deg2rad)*      &
            & sin(src_grid%lon(i)*deg2rad)
       src_grdloc(3,i) = rearth_equator*sin(src_grid%lat(i)*deg2rad)

    end do ! do i = 1, src_grid%ncoords

    ! Compute local variables

    kdtree2 => kdtree2_create(src_grdloc,sort=.true.,rearrange=.true.)

    ! Loop through local variables

    do i = 1, dst_grid%ncoords

       ! Compute local variables

       dst_grdloc(1) = rearth_equator*cos(dst_grid%lat(i)*deg2rad)*        &
            & cos(dst_grid%lon(i)*deg2rad)
       dst_grdloc(2) = rearth_equator*cos(dst_grid%lat(i)*deg2rad)*        &
            & sin(dst_grid%lon(i)*deg2rad)
       dst_grdloc(3) = rearth_equator*sin(dst_grid%lat(i)*deg2rad)
       
       ! Define local variables

       call kdtree2_n_nearest(tp=kdtree2,qv=dst_grdloc,nn=kdtree%nn,       &
            & results=sresults)
       kdtree%r2dist(i,1:kdtree%nn) = sresults(1:kdtree%nn)%dis
       kdtree%idx(i,1:kdtree%nn)    = sresults(1:kdtree%nn)%idx

    end do ! do i = 1, dst_grid%ncoords

    ! Deallocate memory for local variables

    call kdtree2_destroy(kdtree2)
    if(allocated(src_grdloc)) deallocate(src_grdloc)

    !=====================================================================

  end subroutine math_methods_kdtree_nn

  !=======================================================================

  ! SUBROUTINE:

  ! math_methods_kdtree_r2.f90

  ! DESCRIPTION:

  ! This subroutine finds the N nearest-neighbors within a user
  ! specified radius of a given location.

  ! REFERENCES:

  ! Kennel, M. B., 2004: KDTREE2: Fortran 95 and C++ software to
  ! efficiently search for near neighbors in a multi-dimensional
  ! Euclidean space. 

  ! http://arxiv.org/PScache/phvsics/pdf/0408/0408067.pdf.

  ! INPUT VARIABLES:
  
  ! * dst_grid; a FORTRAN grid_struct variable containing the
  !   Euclidean space geographical locations within which to find
  !   nearest-neighbors; geographical locations (e.g., latitudes and
  !   longitudes) are assumed to have units of degrees.

  ! * src_grid; a FORTRAN grid_struct variable containing the
  !   geographical locations for which to find the N-nearest neighbor
  !   locations; geographical locations (e.g., latitudes and
  !   longitudes) are assumed to have units of degrees.

  ! * kdtree; a FORTRAN kdtree_struct variable containing (at minimum)
  !   the number of coordinate values (the ncoords attribute should be
  !   equal to the dst_grid variable ncoords attribute) and the radial
  !   distance within which to seek the nearest neighbors.

  ! OUTPUT VARIABLES:

  ! * kdtree; a FORTRAN kdtree_struct variable containing the
  !   N-nearest neighbors within the R^2 distance specified by the
  !   user.

  !-----------------------------------------------------------------------

  subroutine math_methods_kdtree_r2(src_grid,dst_grid,kdtree)

    ! Define variables passed to routine

    type(grid_struct)                                                   :: dst_grid
    type(grid_struct)                                                   :: src_grid
    type(kdtree_struct)                                                 :: kdtree

    ! Define variables computed within routine

    type(kdtree2), pointer                                              :: kdtree2
    type(kdtree2_result)                                                :: sresults(kdtree%nalloc)
    real(r_kind),               dimension(:,:),             allocatable :: src_grdloc
    real(r_kind)                                                        :: dst_grdloc(3)

    ! Define counting variables

    integer                                                             :: i

    !=====================================================================

    ! Allocate memory for local variables

    if(.not. allocated(src_grdloc))                                        &
         & allocate(src_grdloc(3,dst_grid%ncoords))

    ! Loop through local variables

    do i = 1, dst_grid%ncoords

       ! Compute local variables

       src_grdloc(1,i) = rearth_equator*cos(dst_grid%lat(i)*deg2rad)*      &
            & cos(dst_grid%lon(i)*deg2rad)
       src_grdloc(2,i) = rearth_equator*cos(dst_grid%lat(i)*deg2rad)*      &
            & sin(dst_grid%lon(i)*deg2rad)
       src_grdloc(3,i) = rearth_equator*sin(dst_grid%lat(i)*deg2rad)

    end do ! do i = 1, dst_grid%ncoords

    ! Compute local variables

    kdtree2 => kdtree2_create(src_grdloc,sort=.true.,rearrange=.true.)
    dst_grdloc(1) = rearth_equator*cos(src_grid%clat*deg2rad)*             &
         & cos(src_grid%clon*deg2rad)
    dst_grdloc(2) = rearth_equator*cos(src_grid%clat*deg2rad)*             &
         & sin(src_grid%clon*deg2rad)
    dst_grdloc(3) = rearth_equator*sin(src_grid%clat*deg2rad)
       
    ! Define local variables

    call kdtree2_r_nearest(tp=kdtree2,qv=dst_grdloc,r2=kdtree%r2,nfound=   &
         & kdtree%nfound,nalloc=kdtree%nalloc,results=sresults)
    kdtree%r2dist(1,1:kdtree%nfound) = sresults(1:kdtree%nfound)%dis
    kdtree%idx(1,1:kdtree%nfound)    = sresults(1:kdtree%nfound)%idx

    ! Deallocate memory for local variables

    call kdtree2_destroy(kdtree2)
    if(allocated(src_grdloc)) deallocate(src_grdloc)

    !=====================================================================

  end subroutine math_methods_kdtree_r2

  !=======================================================================

  ! SUBROUTINE:

  ! math_methods_polarcoords.f90

  ! DESCRIPTION:

  ! This subroutine will compute the polar coordinate values of radial
  ! distance and azimuthal angle relative to a reference geographical
  ! location.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN grid_struct variable containing grid longitude
  !   and latitude values; units for longitude and latitude are
  !   degrees.

  ! * reflat; a FORTRAN 4-byte real value specifying the geographical
  !   latitude coordinate of the reference location; units are
  !   degrees.

  ! * reflon; a FORTRAN 4-byte real value specifying the geographical
  !   longitude coordinate of the reference location; units are
  !   degrees.

  ! OUTPUT VARIABLES:

  ! * grid; a FORTRAN grid_struct variable containing the computed
  !   radial distance relative to a reference Earth-bound geographical
  !   location assuming the equitorial radius of the Earth (units are
  !   meters) and the angle relative to a reference Earth-bound
  !   geographical location assuming the equitorial radius of the
  !   Earth (units are degrees).

  !-----------------------------------------------------------------------

  subroutine math_methods_polarcoords(grid,reflon,reflat)

    ! Define variables passed to routine

    type(grid_struct)                                                   :: grid
    real(r_kind)                                                        :: reflon
    real(r_kind)                                                        :: reflat

    !=====================================================================

    ! Compute local variables

    call math_methods_radialgrid(grid,reflon,reflat)
    call math_methods_anglegrid(grid,reflon,reflat)

    !=====================================================================

  end subroutine math_methods_polarcoords

  !=======================================================================

  ! SUBROUTINE: 

  ! math_methods_radialdist.f90

  ! DESCRIPTION:

  ! This subroutine computes the radial distance, using the
  ! Earth-bound geographical coordinates, between two locations
  ! assuming the equitorial radius of the Earth.

  ! INPUT VARIABLES:

  ! * lon1; a FORTRAN 4-byte real value specifying the geographical
  !   longitude coordinate of location 1; units are degrees.

  ! * lat1; a FORTRAN 4-byte real value specifying the geographical
  !   latitude coordinate of location 1; units are degrees.

  ! * lon2; a FORTRAN 4-byte real value specifying the geographical
  !   longitude coordinate of location 2; units are degrees.

  ! * lat2; a FORTRAN 4-byte real value specifying the geographical
  !   latitude coordinate of location 2; units are degrees.

  ! OUTPUT VARIABLES:

  ! * dist; a FORTRAN 4-byte real value specifying the radial distance
  !   between two locations assuming the equitorial radius of the
  !   Earth; units are meters.

  !-----------------------------------------------------------------------

  subroutine math_methods_radialdist(lon1,lat1,lon2,lat2,dist)

    ! Define variables passed to routine

    real(r_kind)                                                        :: lon1
    real(r_kind)                                                        :: lat1
    real(r_kind)                                                        :: lon2
    real(r_kind)                                                        :: lat2
    real(r_kind)                                                        :: dist

    ! Define variables computed within routine

    real(r_double)                                                      :: lat1_deg
    real(r_double)                                                      :: lat2_deg
    real(r_double)                                                      :: lat1_rad
    real(r_double)                                                      :: lat2_rad
    real(r_double)                                                      :: lon1_deg
    real(r_double)                                                      :: lon2_deg
    real(r_double)                                                      :: lon1_rad
    real(r_double)                                                      :: lon2_rad
    real(r_double)                                                      :: x1
    real(r_double)                                                      :: y1
    real(r_double)                                                      :: z1
    real(r_double)                                                      :: x2
    real(r_double)                                                      :: y2
    real(r_double)                                                      :: z2
    real(r_kind)                                                        :: dr

    !=====================================================================

    ! Compute local variables

    lon1_deg = lon1
    lon1_deg = mod(lon1_deg,360.d0)
    if(lon1_deg .lt. 0.d0) then
       lon1_deg = lon1_deg + 360.d0
    end if ! if(lon1_deg .lt. 0.d0)
    lon2_deg = lon2
    lon2_deg = mod(lon2_deg,360.d0)
    if(lon2_deg .lt. 0.d0) then
       lon2_deg = lon2_deg + 360.d0
    end if ! if(lon2_deg .lt. 0.d0)
    if(lat1 .lt. lat2) then
       lon1_rad = lon1_deg*deg2rad
       lat1_rad = (90.d0 - lat1)*deg2rad
       lon2_rad = lon2_deg*deg2rad
       lat2_rad = (90.d0 - lat2)*deg2rad
    else if(lat1 .eq. lat2 .and. lon1_deg .le. lon2_deg) then
       lon1_rad = lon1_deg*deg2rad
       lat1_rad = (90.d0 - lat1)*deg2rad
       lon2_rad = lon2_deg*deg2rad
       lat2_rad = (90.d0 - lat2)*deg2rad
    else
       lon1_rad = lon2_deg*deg2rad
       lat1_rad = (90.d0 - lat2)*deg2rad
       lon2_rad = lon1_deg*deg2rad
       lat2_rad = (90.d0 - lat1)*deg2rad
    end if ! if(lat1 .lt. lat2)
    x1   = sin(lat1_rad)*cos(lon1_rad)
    y1   = sin(lat1_rad)*sin(lon1_rad)
    z1   = cos(lat1_rad)
    x2   = sin(lat2_rad)*cos(lon2_rad)
    y2   = sin(lat2_rad)*sin(lon2_rad)
    z2   = cos(lat2_rad)
    dr   = acos(min(1.d0,x1*x2 + y1*y2 + z1*z2))
    dist = dr*rearth_equator

    !=====================================================================

  end subroutine math_methods_radialdist

  !=======================================================================

  ! SUBROUTINE: 

  ! math_methods_radialgrid.f90

  ! DESCRIPTION:

  ! This subroutine computes a grid of radial distance values relative
  ! to a reference Earth-bound geographical location and a grid of
  ! longitude and latitude values assuming the equitorial radius of
  ! the Earth.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN grid_struct variable containing grid longitude
  !   and latitude values; units for longitude and latitude are
  !   degrees.

  ! * reflat; a FORTRAN 4-byte real value specifying the geographical
  !   latitude coordinate of the reference location; units are
  !   degrees.

  ! * reflon; a FORTRAN 4-byte real value specifying the geographical
  !   longitude coordinate of the reference location; units are
  !   degrees.

  ! OUTPUT VARIABLES:

  ! * grid; a FORTRAN grid_struct variable containing the computed
  !   radial distance relative to a reference Earth-bound geographical
  !   location assuming the equitorial radius of the Earth; units are
  !   meters.

  !-----------------------------------------------------------------------

  subroutine math_methods_radialgrid(grid,reflon,reflat)

    ! Define variables passed to routine

    type(grid_struct)                                                   :: grid
    real(r_kind)                                                        :: reflon
    real(r_kind)                                                        :: reflat

    ! Define counting variables

    integer                                                             :: i

    !=====================================================================

    ! Loop through local variable

    do i = 1, (grid%nx*grid%ny)

       ! Compute local variables

       call math_methods_radialdist(reflon,reflat,grid%lon(i),             &
            & grid%lat(i),grid%radius(i))

    end do ! i = 1, (grid%nx*grid%ny)

    !=====================================================================
    
  end subroutine math_methods_radialgrid

  !=======================================================================

  ! SUBROUTINE:

  ! math_methods_sort_array.f90

  ! DESCRIPTION:

  ! This subroutine implements the SLATEC ssort routine to sort a
  ! dependent array (ya) relative to a sorted independent array (xa);
  ! the arrays may be sorted in either the ascending or descending
  ! direction.

  ! REFERENCES:

  ! Singleton, R.C., 1969. Algorithm 347: an efficient algorithm for
  ! sorting with minimal storage [M1]. Communications of the ACM,
  ! 12(3), pp.185-186.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN spline_struct containing the independent (xa)
  !   and dependent (ya) variable arrays to be sorted.

  ! * ascend; a FORTRAN logical variable specifying whether the arrays
  !   are to be sorted in the ascending direction.

  ! * descend; a FORTRAN logical variable specifying whether the
  !   arrays are to be sorted in the descending direction.

  ! OUTPUT VARIABLES:

  ! * grid; a FORTRAN spline_struct containing the sorted independent
  !   (xa) and dependent (ya) variable arrays.

  !-----------------------------------------------------------------------

  subroutine math_methods_sort_array(grid,ascend,descend)

    ! Define variables passed to routine

    type(spline_struct)                                                 :: grid
    logical                                                             :: ascend
    logical                                                             :: descend

    ! Define variables computed within routine

    real(r_kind),               dimension(:),               allocatable :: xa
    real(r_kind),               dimension(:),               allocatable :: ya
    integer                                                             :: kflag
    integer                                                             :: n

    !=====================================================================

    ! Check local variable and proceed accordingly

    if(ascend)  kflag = 2
    if(descend) kflag = -2

    ! Define local variables

    n = grid%n

    ! Check local variable and proceed accordingly

    if(n .ge. 2) then

       ! Allocate memory for local variables
       
       if(.not. allocated(xa)) allocate(xa(n))
       if(.not. allocated(ya)) allocate(ya(n))

       ! Define local variables

       xa      = grid%xa
       ya      = grid%ya
       call ssort(xa,ya,n,kflag)
       grid%xa = xa
       grid%ya = ya

       ! Deallocate memory for local variables
       
       if(allocated(xa)) deallocate(xa)
       if(allocated(ya)) deallocate(ya)

    end if ! if(n .ge. 2)

    !=====================================================================

  end subroutine math_methods_sort_array

  !=======================================================================

  ! SUBROUTINE: 

  ! math_methods_stats.f90

  ! DESCRIPTION:

  ! This subroutine defines/computes the attributes of a variable
  ! array and returns the respective attibutes.

  ! INPUT VARIABLES:

  ! * vargrid; a FORTRAN vargrid_struct variable.

  ! OUTPUT VARIABLES:

  ! * statgrid; a FORTRAN statgrid_struct variable.

  !-----------------------------------------------------------------------

  subroutine math_methods_stats(vargrid,statgrid)

    ! Define variables passed to routine

    type(vargrid_struct)                                                :: vargrid
    type(statgrid_struct)                                               :: statgrid

    ! Define variables computed within routine

    real(r_kind)                                                        :: sum
    real(r_kind)                                                        :: sumsq
    real(r_kind)                                                        :: varmin
    real(r_kind)                                                        :: varmax
    integer                                                             :: count

    ! Define counting variables

    integer                                                             :: i

    !=====================================================================

    ! Define local variables

    call init_stats(statgrid)
    varmin = spval
    varmax = -spval
    sum    = 0.0
    count  = 0

    ! Loop through local variable

    do i = 1, vargrid%nvals

       ! Check local variable and proceed accordingly

       if(vargrid%var(i) .ne. spval) then

          ! Define local variables

          varmin = min(varmin,vargrid%var(i))
          varmax = max(varmax,vargrid%var(i))

          ! Compute local variables

          sum   = sum + vargrid%var(i)
          count = count + 1

       end if ! if(vargrid%var(i) .ne. spval)
       
    end do ! do i = 1, vargrid%nvals

    ! Check local variable and proceed accordingly

    if(abs(varmax) .eq. spval) varmax          = spval
    if(varmin .ne. spval)      statgrid%varmin = varmin
    if(varmax .ne. spval)      statgrid%varmax = varmax
    if(count .gt. 0)           statgrid%mean   = sum/real(count)
    if(statgrid%mean .ne. spval) then
    
       ! Define local variables

       sum   = 0.0
       sumsq = 0.0

       ! Loop through local variable

       do i = 1, vargrid%nvals

          ! Check local variable and proceed accordingly

          if(vargrid%var(i) .ne. spval) then

             ! Compute local variables

             sum   = sum   + (vargrid%var(i) - statgrid%mean)
             sumsq = sumsq + ((vargrid%var(i) - statgrid%mean)*           &
                  & (vargrid%var(i) - statgrid%mean))

          end if ! if(vargrid%var(i) .ne. spval)

       end do ! do i = 1, vargrid%nvals

       ! Check local variable and proceed accordingly

       if(count .gt. 1) then

          ! Compute local variables

          statgrid%vari  = (sumsq - (sum*sum)/count)/(count - 1)
          statgrid%stdev = sqrt(statgrid%vari)

       end if ! if(count .gt. 1)

    end if ! if(statgrid%mean .ne. spval)

    ! Define local variables

    statgrid%nvals = count

    !=====================================================================

  end subroutine math_methods_stats

  !=======================================================================

  ! SUBROUTINE:

  ! unique_array_spline.f90

  ! DESCRIPTION:

  ! This subroutine returns an array containing only unique values;
  ! this method is utilitized by spline routines which require
  ! absolute ascending/descending arrays/lists in order to perform
  ! their respective algorithms.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN spline_struct variable containing the values
  !   defined by the user for spline calculations;

  ! OUTPUT VARIABLES:

  ! * grid; a FORTRAN spline_struct variable which is a modified
  !   version of the input FORTRAN spline_struct variable assuming the
  !   FORTRAN spline_struct variable does not contain entirely unique
  !   values.

  !-----------------------------------------------------------------------

  subroutine unique_array_spline(grid)

    ! Define variables passed to routine

    type(spline_struct)                                                 :: grid

    ! Define variables computed within routine

    type(spline_struct)                                                 :: gridl
    real(r_kind)                                                        :: max_val
    real(r_kind)                                                        :: min_val
    integer                                                             :: nn

    ! Define counting variables

    integer                                                             :: i

    !=====================================================================

    ! Define local variables

    gridl%n  = grid%n
    call variable_interface_setup_struct(gridl)
    gridl%xa = spval
    gridl%ya = spval
    max_val  = maxval(grid%xa)
    min_val  = minval(grid%xa)
    nn       = 0

    ! Loop through local variable

    do while(min_val .lt. max_val)

       ! Define local variables

       nn           = nn + 1
       min_val      = minval(grid%xa,mask=(grid%xa .gt. min_val))
       gridl%xa(nn) = min_val

       ! Loop through local variable

       do i = 1, grid%n

          ! Check local variable and proceed accordingly

          if(gridl%xa(nn) .eq. grid%xa(i)) then

             ! Define local variables

             gridl%ya(nn) = grid%ya(i)
             goto 1000

          end if ! if(gridl%xa(nn) .eq. grid%xa(i))

       end do ! do i = 1, grid%n
       
       ! Define local variables
       
1000   continue
       
    end do ! do while(min_val .lt. max_val)

    ! Deallocate memory for local variables

    call variable_interface_cleanup_struct(grid)

    ! Define local variables

    grid%n  = count(gridl%xa .ne. spval)
    call variable_interface_setup_struct(grid)
    grid%xa = gridl%xa(1:grid%n)
    grid%ya = gridl%ya(1:grid%n)

    ! Deallocate memory for local variables

    call variable_interface_cleanup_struct(gridl)

    !=====================================================================

  end subroutine unique_array_spline

  !=======================================================================

  ! SUBROUTINE: 

  ! wvnbrdcmp2d.f90

  ! DESCRIPTION:

  ! This subroutine computes the Fast Fourier transform (FFT) of a
  ! 2-dimensional variable field and subsequently decomposes the
  ! respective variable field into it's wave number components; the
  ! total number of wave numbers is equal to (nx/2 + 1).

  ! INPUT VARIABLES:

  ! * wnd2d; a FORTRAN wnd2d_struct variable where the zonal- and
  !   meridional-dimensions (i.e., nx and ny, respectively) and the
  !   total wavenumber (nh) are defined in addition to the
  !   2-dimensional variable field to be decomposed (var).

  ! * remap_grid; a FORTRAN recenter_struct variable containing the
  !   radial and azimuthal (i.e., polar-projection) components of the
  !   analysis variable grid to be utilized by the FFT.

  ! * grid; a FORTRAN grid_struct variable.

  ! OUTPUT VARIABLES:

  ! * wnd2d; a FORTRAN wnd2d_struct variable containing the wavenumber
  !   decomposition (wnvar) of the variable field specified by the
  !   user (var).

  ! * grid; a FORTRAN grid_struct variable containing the grid
  !   projection information required to remap the respective
  !   wavenumbers to the original grid projection.

  !-----------------------------------------------------------------------

  subroutine wvnbrdcmp2d(wnd2d,remap_grid,grid)

    ! Define variables passed to routine

    type(wnd2d_struct)                                                  :: wnd2d
    type(recenter_struct)                                               :: remap_grid
    type(grid_struct)                                                   :: grid

    ! Define variables computed within routine

    type(fft2d_struct)                                                  :: fft2d
    type(vargrid_struct)                                                :: vargrid
    type(statgrid_struct)                                               :: statgrid
    complex(r_double),          dimension(:,:),             allocatable :: fftsave
    complex(r_double),          dimension(:,:),             allocatable :: fftvar
    real(r_kind),               dimension(:),               allocatable :: var
    real(r_kind)                                                        :: radius_min
    real(r_kind)                                                        :: radius_max
    real(r_kind)                                                        :: angle_min
    real(r_kind)                                                        :: angle_max
    integer                                                             :: ncoord

    ! Define counting variables

    integer                                                             :: i, j

    !=====================================================================

    ! Define local variables

    fft2d%nx      = wnd2d%nr
    fft2d%ny      = wnd2d%na
    call variable_interface_setup_struct(fft2d)
    vargrid%nx    = wnd2d%nx
    vargrid%ny    = wnd2d%ny
    vargrid%nvals = (vargrid%nx*vargrid%ny)
    call variable_interface_setup_struct(vargrid) 
    grid%nx       = fft2d%nx
    grid%ny       = fft2d%ny
    grid%ncoords  = (grid%nx*grid%ny)
    grid%nz       = remap_grid%nz
    call variable_interface_setup_struct(grid)

    ! Allocate memory for local variables

    if(.not. allocated(fftsave)) allocate(fftsave(fft2d%nx,fft2d%ny))
    if(.not. allocated(fftvar))  allocate(fftvar(fft2d%nx,fft2d%ny))
    if(.not. allocated(var))     allocate(var(wnd2d%nx*wnd2d%ny))

    ! Define local variables

    ncoord = 0

    ! Loop through local variable
    
    do j = 1, fft2d%ny

       ! Loop through local variable
       
       do i = 1, fft2d%nx

          ! Define local variables

          ncoord      = ncoord + 1
          vargrid%var = spval
          var         = spval
          radius_min  = 0.0 + (i-1)*wnd2d%dradius
          radius_max  = radius_min + wnd2d%dradius
          angle_min   = 0.0 + (j-1)*wnd2d%dangle
          angle_max   = angle_min + wnd2d%dangle
          grid%gcdist = radius_max
          grid%gchead = 0.5*(angle_max + angle_min) + 90.0
          grid%gclon  = remap_grid%clon
          grid%gclat  = remap_grid%clat
          where(remap_grid%radius .ge. radius_min .and.                   &
               & remap_grid%radius .lt. radius_max .and.                  &
               & remap_grid%angle .ge. angle_min .and.                    &
               & remap_grid%angle .lt. angle_max) var = wnd2d%var
          vargrid%var = var

          ! Compute local variables

          call grid_methods_gcgeo(grid)
          call math_methods_stats(vargrid,statgrid)
          if(statgrid%mean .eq. spval) statgrid%mean = 0.0

          ! Define local variables

          grid%lon(ncoord) = grid%gclon
          grid%lat(ncoord) = grid%gclat

          ! Define local variables

          fftvar(i,j) = cmplx(dble(statgrid%mean),dble(0.0),kind=r_double)

       end do ! do i = 1, fft2d%nx

    end do ! do j = 1, fft2d%ny

    ! Define local variables

    fft2d%in = fftvar

    ! Compute local variables

    call math_methods_fftw_dft(fft2d,.true.,.false.)

    ! Define local variables

    fftsave = fft2d%out

    ! Loop through local variable

    do i = 1, wnd2d%nh

       ! Define local variables

       fft2d%out      = cmplx(0.0,0.0,kind=r_double)
       fft2d%in       = cmplx(0.0,0.0,kind=r_double)
       fft2d%out(:,i) = fftsave(:,i)

       ! Compute local variables

       call math_methods_fftw_dft(fft2d,.false.,.true.)

       ! Define local variables

       wnd2d%wnvar(i,:,:) = real(fft2d%in)

    end do ! do i = 1, wnd2d%nh

    ! Deallocate memory for local variables

    if(allocated(fftsave)) deallocate(fftsave)
    if(allocated(fftvar))  deallocate(fftvar)
    if(allocated(var))     deallocate(var)
    call variable_interface_cleanup_struct(fft2d)
    call variable_interface_cleanup_struct(vargrid)

    !=====================================================================

  end subroutine wvnbrdcmp2d

  !=======================================================================

end module math_methods_interface
