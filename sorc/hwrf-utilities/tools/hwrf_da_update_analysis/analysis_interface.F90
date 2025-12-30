module analysis_interface

  !=======================================================================

  !$$$ MODULE DOCUMENTATION BLOCK
  
  ! da-update-analysis :: analysis_interface
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
  use f95_precision
  use fileio_interface
  use interpolation_interface
  use math_methods_interface
  use namelist_interface
  use variable_interface

  ! Define interfaces and attributes for module routines
  
  implicit none
  private
  public :: analysis_increment
  public :: blend_interface_boundary
  public :: blend_interface_cycling
  public :: blend_interface_nest
  public :: blend_interface_radii
  public :: blend_interface_tc_environment
  public :: blend_interface_wnd

  !-----------------------------------------------------------------------

contains

  !=======================================================================

  ! SUBROUTINE:

  ! analysis_increment.f90

  ! DESCRIPTION:

  ! This subroutine computes the increment applied to the first-guess
  ! (i.e., background) forecast variable, usually via a data
  ! assimilative methodology, to define the analysis (e.g., state)
  ! variable; this subroutine assumes that the first-guess and
  ! analysis variables, used to compute the increment values, are
  ! defined on identical grid projections.

  ! INPUT VARIABLES:

  ! * fg_grid; a FORTRAN grid_struct variable containing the
  !   first-guess (e.g., background) forecast variable.

  ! * an_grid; a FORTRAN grid_struct variable containing the analysis
  !   (i.e., the data-assimialtion updated state) variable.

  ! * ic_grid; a FORTRAN grid_struct variable to contain the
  !   increments (e.g., analysis - background) variable.

  ! OUTPUT VARIABLES:

  ! * ic_grid; a FORTRAN grid_struct variable containing the
  !   increments applied to the first-guess forecast variable to
  !   define the respective analysis variable as well as the
  !   first-guess grid projection.

  !-----------------------------------------------------------------------

  subroutine analysis_increment(fg_grid,an_grid,ic_grid)

    ! Define variables passed to routine

    type(grid_struct)                                                   :: an_grid
    type(grid_struct)                                                   :: fg_grid
    type(grid_struct)                                                   :: ic_grid

    ! Define counting variables

    integer                                                             :: i

    !=====================================================================

    ! Define local variables

    ic_grid%nx = fg_grid%nx
    ic_grid%ny = fg_grid%ny
    ic_grid%nz = fg_grid%nz
    call variable_interface_setup_struct(ic_grid)

    ! Loop through local variable

    do i = 1, ic_grid%nz

       ! Compute local variables

       ic_grid%var(:,i) = (an_grid%var(:,i) - fg_grid%var(:,i))

    end do ! do i = 1, ic_grid%nz

    !=====================================================================

  end subroutine analysis_increment

  !=======================================================================

  ! SUBROUTINE:

  ! blend_interface_boundary.f90

  ! DESCRIPTION:

  ! This subroutine blends two analysis variables; the projections for
  ! the input and output grids (in_grid and out_grid, respectively)
  ! are assumed identical; the user configuration defines the region
  ! within the input grid to merge (e.g., blend) the analysis (see
  ! namelist variable npts_blend) and the number of points within the
  ! merge region to relax to the output grid variable.
  
  ! INPUT VARIABLES:

  ! * in_grid; a FORTRAN grid_struct variable containing the input
  !   grid projection and variable to be relaxed to the output
  !   variable.

  ! * out_grid; a FORTRAN grid_struct variable containing the output
  !   grid projection and the variable to be relaxed to at the edges
  !   of the blending (e.g., merge) region.

  ! OUTPUT VARIABLES:

  ! * out_grid; a FORTRAN grid_struct variable containing the merged
  !   (e.g., blended) variable.

  !-----------------------------------------------------------------------

  subroutine blend_interface_boundary(in_grid,out_grid)

    ! Define variables passed to routine

    type(grid_struct)                                                   :: in_grid
    type(grid_struct)                                                   :: out_grid
    
    ! Define variables computed within routine

    type(blend_struct)                                                  :: blend
    type(merge_struct)                                                  :: merge
    type(slint_struct)                                                  :: slint
    real(r_kind),               dimension(:,:,:),           allocatable :: in_work
    real(r_kind),               dimension(:,:,:),           allocatable :: out_work
    integer                                                             :: dst_idx
    integer                                                             :: bmnx
    integer                                                             :: bmny
    integer                                                             :: bmxx
    integer                                                             :: bmxy
    integer                                                             :: rmnx
    integer                                                             :: rmny
    integer                                                             :: rmxx
    integer                                                             :: rmxy

    ! Define counting variables                  

    integer                                                             :: i, j, k

    !=====================================================================

    ! Allocate memory for local variables

    if(.not. allocated(in_work))                                           &
         & allocate(in_work(in_grid%nx,in_grid%ny,in_grid%nz))
    if(.not. allocated(out_work))                                          &
         & allocate(out_work(out_grid%nx,out_grid%ny,out_grid%nz))

    ! Define local variables

    bmnx         = min(npts_blend,in_grid%nx)
    bmny         = min(npts_blend,in_grid%ny)
    bmxx         = max((in_grid%nx - npts_blend),1)
    bmxy         = max((in_grid%ny - npts_blend),1)
    rmnx         = min(max(bmnx,(bmnx + npts_relax)),in_grid%nx)
    rmny         = min(max(bmny,(bmny + npts_relax)),in_grid%ny)
    rmxx         = max(min(bmxx,(bmxx - npts_relax)),1)
    rmxy         = max(min(bmxy,(bmxy - npts_relax)),1)
    merge%bmnx   = bmnx
    merge%bmxx   = bmxx
    merge%bmny   = bmny
    merge%bmxy   = bmxy
    merge%rmnx   = rmnx
    merge%rmxx   = rmxx
    merge%rmny   = rmny
    merge%rmxy   = rmxy
    merge%nx     = in_grid%nx
    merge%ny     = in_grid%ny
    merge%nz     = in_grid%nz
    call variable_interface_setup_struct(merge)
    merge%lat_an = out_grid%lat
    merge%lon_an = out_grid%lon
    merge%lat_fg = in_grid%lat
    merge%lon_fg = in_grid%lon

    ! Loop through local variable

    do k = 1, out_grid%nz

       ! Define local variables

       in_work(:,:,k)  = reshape(in_grid%var(:,k),shape(in_work(:,:,k)))   
       out_work(:,:,k) = reshape(out_grid%var(:,k),                        &
            & shape(out_work(:,:,k))) 

       ! Loop through local variable

       do j = rmny, rmxy

          ! Loop through local variable

          do i = rmnx, rmxx

             ! Define local variables

             in_work(i,j,k) = out_work(i,j,k)

          end do ! do i = rmnx, rmxx

       end do ! do j = rmny, rmxy

       ! Loop through local variable

       do j = rmny, rmxy

          ! Loop through local variable

          do i = bmnx, rmnx

             ! Define local variables

             in_work(i,j,k) = spval

          end do ! do i = bmnx, rmnx

          ! Loop through local variable

          do i = rmxx, bmxx

             ! Define local variables

             in_work(i,j,k) = spval

          end do ! do i = rmxx, bmxx

       end do ! do j = bmny, rmny

       ! Loop through local variable

       do i = rmnx, rmxx

          ! Loop through local variable

          do j = bmny, rmny

             ! Define local variables

             in_work(i,j,k) = spval

          end do ! do j = bmny, rmny

          ! Loop through local variable

          do j = rmxy, bmxy

             ! Define local variables

             in_work(i,j,k) = spval

          end do ! do j = rmxy, bmxy

       end do ! do i = rmnx, rmxx

       ! Define local variables

       in_grid%var(:,k) = reshape(in_work(:,:,k),shape(in_grid%var(:,k)))

    end do ! do k = 1, out_grid%nz

    ! Define local variables

    merge%npts_blend = npts_blend
    merge%npts_relax = npts_relax
    merge%var        = in_grid%var
    call define_mask_boundary(merge)
    merge%var_fg     = in_grid%var 
    merge%var_an     = out_grid%var
    blend%nx         = merge%nx
    blend%ny         = merge%ny
    blend%nz         = merge%nz
    call variable_interface_setup_struct(blend)
    blend%lat        = merge%lat_an
    blend%lon        = merge%lon_an
    blend%mask       = merge%mask
    blend%npts_blend = merge%npts_blend
    blend%npts_relax = merge%npts_relax

    ! Loop through local variable
    
    do j = 1, blend%nz

       ! Define local variables

       blend%var(:,j) = merge%var_fg(:,j)
       where(blend%var(:,j) .eq. spval) blend%var(:,j) =                   &
            & merge%var_an(:,j)

    end do ! do j = 1, blend%nz

    ! Define local variables

    blend%var_blend = blend%var

    ! Loop through local variable

    do k = 1, (blend%npts_blend + blend%npts_relax)

       ! Compute local variables

       if(is_slint) call grid_blend(blend,slint,real(k))

       ! Loop through local variable

       do j = 1, blend%nz

          ! Compute local variables

          if(is_slint) call interpolation_interface_bilinear(slint,        &
               & blend%var(:,j))

          ! Loop through local variable

          do i = 1, blend%ncoords

             ! Define local variables

             dst_idx = blend%idx(i)

             ! Check local variable and proceed accordingly

             if((dst_idx .gt. 0) .and. (dst_idx .le. (blend%nx*blend%ny))  &
                  & .and. (size(slint%var) .gt. 0)) then

                ! Check local variable and proceed accordingly

                if(is_slint) then

                   ! Define local variables

                   blend%var_blend(i,j) = slint%var(dst_idx)

                end if ! if(is_slint)

             end if ! if((dst_idx .gt. 0) .and. (dst_idx
                    ! .le. (blend%nx*blend%ny)) .and. (size(slint%var)
                    ! .gt. 0))

          end do ! do i = 1, blend%ncoords

       end do ! do j = 1, blend%nz

       ! Define local variables

       where(blend%mask .eq. real(k)) blend%mask = spval
       blend%idx = -999

       ! Deallocate memory for local variables

       call variable_interface_cleanup_struct(slint)

    end do ! do k = 1, (blend%npts_blend + blend%npts_relax)

    ! Define local variables

    out_grid%var = blend%var_blend

    ! Deallocate memory for local variables

    if(allocated(in_work))  deallocate(in_work)
    if(allocated(out_work)) deallocate(out_work)
    call variable_interface_cleanup_struct(merge)
    call variable_interface_cleanup_struct(blend)

    !=====================================================================

  end subroutine blend_interface_boundary

  !=======================================================================

  ! SUBROUTINE:

  ! blend_interface_cycling.f90

  ! DESCRIPTION:

  ! This subroutine merges two domains, of identical dimensions, using
  ! a relaxation blending mask for the overlapping regions.

  ! INPUT VARIABLES:

  ! * merge; a FORTRAN merge_struct variable; it is assumed that all
  !   arrays within the structure and the respective blending
  !   attributes have been specified before calling this routine;
  !   also, the lat_an, lon_an, lat_fg, lon_fg, var_an and var_fg
  !   attribute arrays have also been populated.

  ! OUTPUT VARIABLES:

  ! * merge; a FORTRAN merge_struct variable containing the merged
  !   grid variable.

  !-----------------------------------------------------------------------

  subroutine blend_interface_cycling(merge)

    ! Define variables passed to routine

    type(merge_struct)                                                  :: merge

    ! Define counting variables

    integer                                                             :: i

    !=====================================================================

    ! Define local variables
    
    merge%var = merge%var_an
    
    ! Compute local variables
    
    call define_overlap_mask(merge)

    ! Loop through local variable

    do i = 1, merge%nz
    
       ! Define local variables

       merge%var(:,i) = merge%var_fg(:,i)*merge%mask +             &
            & merge%var_an(:,i)*(1.0 - merge%mask)

    end do ! do i = 1, merge%nz
    
    !=====================================================================

  end subroutine blend_interface_cycling

  !=======================================================================

  ! SUBROUTINE:

  ! blend_interface_nest.f90

  ! DESCRIPTION:

  ! This subroutine blends a nested grid within the parent grid region
  ! using an interative interpolation method; this method interpolates
  ! from outer to inner; an optional 9-point smoothing algorithm may
  ! be applied to smooth any gradients which may arise from large
  ! discontinuities between nested/blended variable fields.

  ! INPUT VARIABLES:

  ! * blend; a FORTRAN blend_struct variable; it is assumed that all
  !   arrays within the structure and the respective blending
  !   attributes have been specified before calling this routine;
  !   also, the var_fixed and var_blend attribute arrays have also
  !   been populated.

  ! OUTPUT VARIABLES:

  ! * blend; a FORTRAN blend_struct variable containing the blended
  !   nest and parent grid variables. 

  !-----------------------------------------------------------------------

  subroutine blend_interface_nest(blend)

    ! Define variables passed to routine

    type(blend_struct)                                                  :: blend

    ! Define variables computed within routine

    type(grid_struct)                                                   :: grid
    real(r_kind),               dimension(:),               allocatable :: mask_relax
    real(r_kind),               dimension(:),               allocatable :: nest_incr
    
    ! Define counting variables

    integer                                                             :: i

    !=====================================================================

    ! Define local variables

    grid%nx  = blend%nx
    grid%ny  = blend%ny
    call variable_interface_setup_struct(grid)
    grid%lat = blend%lat
    grid%lon = blend%lon
    
    ! Compute local variables

    call math_methods_polarcoords(grid,blend%clon,blend%clat)

    ! Allocate memory for local variables

    if(.not. allocated(mask_relax)) allocate(mask_relax(blend%ncoords))
    if(.not. allocated(nest_incr))  allocate(nest_incr(blend%ncoords))

    ! Define local variables

    mask_relax = (grid%radius - blend%tc_nest_min_radius)/                 &
         & (blend%tc_nest_max_radius - blend%tc_nest_min_radius)
    where(mask_relax .le. 0.0) mask_relax = 0.0
    where(mask_relax .ge. 1.0) mask_relax = 1.0
    mask_relax = 1.0 - mask_relax
    
    ! Loop through local variable

    do i = 1, blend%nz

       ! Define local variables
       
       nest_incr            = blend%var_fixed(:,i)*(1.0 - mask_relax)      &
            & + blend%var_blend(:,i)*mask_relax
       blend%var_blend(:,i) = nest_incr
       
    end do ! do i = 1, blend%nz

    ! Deallocate memory for local variables

    if(allocated(mask_relax)) deallocate(mask_relax)
    if(allocated(nest_incr))  deallocate(nest_incr)
    call variable_interface_cleanup_struct(grid)
    
    !=====================================================================

  end subroutine blend_interface_nest
  
  !=======================================================================

  ! SUBROUTINE:

  ! blend_interface_radii.f90

  ! DESCRIPTION:

  ! This subroutine computes the nest to parent blending radii
  ! thresholds; the radii are determined by the minimum value(s) for
  ! the maximum nest radii and then scaled according to the namelist
  ! parameter tc_nest_ratio_blend.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN grid_struct variable containing the nested-grid
  !   and it's respective grid projection attributes.

  ! * blend; a FORTRAN blend_struct variable.

  ! OUTPUT VARIABLES:

  ! * blend; a FORTRAN blend_struct variable containing the threshold
  !   variable values to be used for blending the nest and parent
  !   domain co-located regions (see the tc_nest_max_radius and
  !   tc_nest_min_radius attributes).

  !-----------------------------------------------------------------------

  subroutine blend_interface_radii(grid,blend)

    ! Define variables passed to routine

    type(blend_struct)                                                  :: blend
    type(grid_struct)                                                   :: grid

    ! Define variables computed within routine

    real(r_kind),               dimension(:,:),             allocatable :: radius
    real(r_kind)                                                        :: maxradii

    ! Define counting variables

    integer                                                             :: i

    !=====================================================================

    ! Compute local variables

    call math_methods_polarcoords(grid,blend%clon,blend%clat)

    ! Allocate memory for local variables

    if(.not. allocated(radius)) allocate(radius(grid%nx,grid%ny))
    
    ! Define local variables

    radius   = reshape(grid%radius,shape(radius))
    maxradii = rearth_equator

    ! Loop through local variable

    do i = 1, grid%nx

       ! Define local variables

       maxradii = min(maxradii,maxval(radius(i,:)))

    end do ! do i = 1, grid%nx    
    
    ! Loop through local variable

    do i = 1, grid%ny

       ! Define local variables

       maxradii = min(maxradii,maxval(radius(:,i)))

    end do ! do i = 1, grid%ny

    ! Define local variables

    blend%tc_nest_max_radius = maxradii
    blend%tc_nest_min_radius = blend%tc_nest_max_radius*                   &
         & tc_nest_ratio_blend
    
    ! Deallocate memory for local variables

    if(allocated(radius)) deallocate(radius)
       
    !=====================================================================

  end subroutine blend_interface_radii

  !=======================================================================

  ! SUBROUTINE:

  ! blend_interface_tc_environment.f90

  ! DESCRIPTION:

  ! This subroutine blends the an analysis variable field with a
  ! first-guess variable (e.g., environmental) field.

  ! INPUT VARIABLES:

  ! * grid_fg; a FORTRAN grid_struct variable containing the
  !   first-guess variable (and it's respective grid projection
  !   attributes).

  ! * grid_an; a FORTRAN grid_struct variable containing the analysis
  !   variable (and it's respective grid projection attributes) which
  !   to blend in accordance with the user specifications.

  ! OUTPUT VARIABLES:

  ! * grid_an; a FORTRAN grid_struct variable containing the analysis
  !   variable blended, for each TC, in accordance with the user
  !   specifications.

  !-----------------------------------------------------------------------

  subroutine blend_interface_tc_environment(grid_fg,grid_an)

    ! Define variables passed to routine

    type(grid_struct)                                                   :: grid_an
    type(grid_struct)                                                   :: grid_fg

    ! Define variables computed within routine

    type(tcv_struct),           dimension(:),               allocatable :: tcv
    type(grid_struct)                                                   :: dst_grid
    type(grid_struct)                                                   :: dst_sub_grid
    type(grid_struct)                                                   :: src_grid

    ! Define counting variables

    integer                                                             :: i, j

    !=====================================================================

    ! Define local variables

    call fileio_interface_read(tcv_filename,tcv)
    dst_grid%nx  = grid_an%nx
    dst_grid%ny  = grid_an%ny
    dst_grid%nz  = grid_an%nz
    call variable_interface_setup_struct(dst_grid)
    dst_grid%lat = grid_an%lat
    dst_grid%lon = grid_an%lon
    dst_grid%var = grid_fg%var

    ! Check local variable and proceed accordingly

    if(ic_blend_opt .ne. 0) then

       ! Loop through local variable

       do i = 1, size(tcv)

          ! Compute local variables

          call math_methods_polarcoords(dst_grid,tcv(i)%lon,tcv(i)%lat)

          ! Loop through local variable

          do j = 1, grid_an%nz

             ! Define local variables

             where(dst_grid%radius .lt.                                    &
                  & ((1.0 - region_blend_ratio)*tc_large_scale_radius))    &
                  & dst_grid%var(:,j) = grid_an%var(:,j)
             where(dst_grid%radius .ge.                                    &
                  & ((1.0 - region_blend_ratio)*tc_large_scale_radius))    &
                  & dst_grid%var(:,j) = spval
             where(dst_grid%radius .gt. tc_large_scale_radius)             &
                  & dst_grid%var(:,j) = grid_fg%var(:,j) 

          end do ! do j = 1, grid_an%nz

          ! Compute local variables

          call environment_blend(dst_grid,grid_an,grid_fg)

          ! Define local variables

          grid_an%var = dst_grid%var

       end do ! do i = 1, size(tcv)

    end if ! if(ic_blend_opt .ne. 0)

    ! Deallocate memory for local variables

    if(allocated(tcv)) deallocate(tcv)
    call variable_interface_cleanup_struct(dst_grid)

    !=====================================================================

  end subroutine blend_interface_tc_environment

  !=======================================================================

  ! SUBROUTINE:

  ! blend_interface_wnd.f90

  ! DESCRIPTION:

  ! This subroutine blends the analysis increments surrounding
  ! tropical cyclone (TC) events in accordance with the user
  ! specifications; this subroutine currently supports the following
  ! methods (e.g., ic_blend_opt values):

  ! 1: compute the wavenumber spectra for the respective variable and
  !    blend only up to and including the maximum wavenumber specified
  !    by the user (wndcmp_mxwvn).

  ! 2: set all variable values within the user specified TC inner-core
  !    radius to zero.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN grid_struct variable containing the analysis
  !   increment variable (and it's respective grid projection
  !   attributes) which to blend in accordance with the user
  !   specifications.

  ! OUTPUT VARIABLES:

  ! * grid; a FORTRAN grid_struct variable containing the analysis
  !   increment variable blended, for each TC, in accordance with the
  !   user specifications.

  !-----------------------------------------------------------------------

  subroutine blend_interface_wnd(grid)

    ! Define variables passed to routine

    type(grid_struct)                                                   :: grid

    ! Define variables computed within routine

    type(tcv_struct),           dimension(:),               allocatable :: tcv
    type(grid_struct)                                                   :: dst_grid
    type(grid_struct)                                                   :: dst_sub_grid
    type(grid_struct)                                                   :: src_grid

    ! Define counting variables

    integer                                                             :: i, j

    !=====================================================================

    ! Define local variables

    call fileio_interface_read(tcv_filename,tcv)
    dst_grid%nx  = grid%nx
    dst_grid%ny  = grid%ny
    dst_grid%nz  = grid%nz
    call variable_interface_setup_struct(dst_grid)
    dst_grid%lat = grid%lat
    dst_grid%lon = grid%lon

    ! Check local variable and proceed accordingly

    if(ic_blend_opt .ne. 0) then

       ! Loop through local variable

       do i = 1, size(tcv)

          ! Compute local variables

          call math_methods_polarcoords(dst_grid,tcv(i)%lon,tcv(i)%lat)

          ! Define local variables
          
          grid%radius = dst_grid%radius

          ! Compute local variables
          
          if(ic_blend_opt .eq. 1) call wnd2d_tc_increment(grid,tcv(i),     &
               & dst_grid)
          if(ic_blend_opt .eq. 2) call zero_tc_increment(grid,dst_grid)
          call region_blend(dst_grid,grid)

       end do ! do i = 1, size(tcv)

    end if ! if(ic_blend_opt .ne. 0)

    ! Deallocate memory for local variables

    if(allocated(tcv)) deallocate(tcv)
    call variable_interface_cleanup_struct(dst_grid)

    !=====================================================================

  end subroutine blend_interface_wnd

  !=======================================================================

  ! SUBROUTINE:

  ! define_mask_boundary.f90

  ! DESCRIPTION:

  ! This subroutine defines the blending mask used to merge two
  ! identical dimension grid variables along the boundary edges.

  ! INPUT VARIABLES:

  ! * merge; a FORTRAN merge_struct variable; it is assumed that all
  !   arrays within the structure and the respective blending
  !   attributes have been specified before calling this routine.

  ! OUTPUT VARIABLES:

  ! * merge; a FORTRAN merge_struct variable containing the merging
  !   and blending mask regions.

  !-----------------------------------------------------------------------

  subroutine define_mask_boundary(merge)

    ! Define variables passed to routine

    type(merge_struct)                                                  :: merge

    ! Define variables computed within routine

    type(blend_struct)                                                  :: blend

    !=====================================================================

    ! Define local variables

    merge%mask       = 1.0
    where(merge%var(:,1) .eq. spval) merge%mask = -1.0
    blend%nx         = merge%nx
    blend%ny         = merge%ny
    blend%nz         = merge%nz
    blend%bmnx       = merge%bmnx
    blend%bmxx       = merge%bmxx
    blend%bmny       = merge%bmny
    blend%bmxy       = merge%bmxy
    blend%rmnx       = merge%rmnx
    blend%rmxx       = merge%rmxx
    blend%rmny       = merge%rmny
    blend%rmxy       = merge%rmxy
    call variable_interface_setup_struct(blend)
    blend%npts_relax = merge%npts_relax
    blend%npts_blend = merge%npts_blend
    blend%mask       = merge%mask
    call mask_relax_boundary(blend)
    merge%mask       = blend%mask
    merge%idx        = -999

    ! Deallocate memory for local variables

    call variable_interface_cleanup_struct(blend)

    !=====================================================================

  end subroutine define_mask_boundary
  
  !=======================================================================

  ! SUBROUTINE:

  ! define_mask_merge.f90

  ! DESCRIPTION:

  ! This subroutine defines the blending mask used to merge two
  ! identical dimension grid variables.

  ! INPUT VARIABLES:

  ! * merge; a FORTRAN merge_struct variable; it is assumed that all
  !   arrays within the structure and the respective blending
  !   attributes have been specified before calling this routine.

  ! OUTPUT VARIABLES:

  ! * merge; a FORTRAN merge_struct variable containing the merging
  !   and blending mask regions.

  !-----------------------------------------------------------------------

  subroutine define_mask_merge(merge)

    ! Define variables passed to routine

    type(merge_struct)                                                  :: merge

    ! Define variables computed within routine

    real(r_kind)                                                        :: mnlat_fg
    real(r_kind)                                                        :: mnlon_fg
    real(r_kind)                                                        :: mxlat_fg
    real(r_kind)                                                        :: mxlon_fg

    !=====================================================================

    ! Define local variables

    call init_mask_merge(merge)
    call mask_relax_merge(merge)
    merge%idx = -999

    !=====================================================================

  end subroutine define_mask_merge

  !=======================================================================

  ! SUBROUTINE:

  ! define_overlap_mask.f90

  ! DESCRIPTION

  ! This subroutine determines the overlap regions for a two
  ! respective domains (e.g., the first-guess denoted by 'fg' and the
  ! analysis denoted by 'an') and defines a relaxation mask of width
  ! 'npts_relax' as specified in the FORTRAN merge_struct variable.

  ! INPUT VARIABLES:

  ! * merge; a FORTRAN merge_struct variable; it is assumed that all
  !   arrays within the structure and the respective blending
  !   attributes have been specified before calling this routine.

  ! OUTPUT VARIABLES:
  
  ! * merge; a FORTRAN merge_struct variable containing the merging
  !   and blending mask regions.

  !-----------------------------------------------------------------------

  subroutine define_overlap_mask(merge)

    ! Define variables passed to routine

    type(merge_struct)                                                  :: merge

    ! Define variables computed within routine
  
    type(grid_struct)                                                   :: grid_an
    type(grid_struct)                                                   :: grid_fg
    type(kdtree_struct)                                                 :: kdtree
    real(r_kind),               dimension(:,:),             allocatable :: mask
    real(r_kind),               dimension(:,:),             allocatable :: mask_save
    real(r_kind)                                                        :: r2dist
    real(r_kind)                                                        :: maxdist
    integer                                                             :: icoord
    integer                                                             :: idx_strt
    integer                                                             :: idx_stop
    integer                                                             :: strt_idx
    integer                                                             :: stop_idx

    ! Define counting variables

    integer                                                             :: i, j, k

    !=====================================================================

    ! Define local variables

    grid_an%nx     = merge%nx
    grid_an%ny     = merge%ny
    call variable_interface_setup_struct(grid_an)
    grid_an%lat    = merge%lat_an 
    grid_an%lon    = merge%lon_an 
    grid_fg%nx     = merge%nx
    grid_fg%ny     = merge%ny
    call variable_interface_setup_struct(grid_fg)
    grid_fg%lat    = merge%lat_fg 
    grid_fg%lon    = merge%lon_fg 
    kdtree%ncoords = (merge%nx*merge%ny)
    kdtree%nn      = 1
    call variable_interface_setup_struct(kdtree)
    merge%mask     = spval
    
    ! Compute local variables

    call math_methods_kdtree_nn(grid_fg,grid_an,kdtree)

    ! Define local variables

    maxdist = -spval

    ! Loop through local variable

    do i = 1, merge%ncoords

       ! Compute local variables

       r2dist = sqrt(merge%dx(i)*merge%dx(i) + merge%dy(i)*merge%dy(i))

       ! Define local variables

       maxdist = max(maxdist,r2dist)

    end do ! do i = 1, merge%ncoords

    ! Define local variables

    where(sqrt(kdtree%r2dist(:,1)) .gt. maxdist) merge%mask = 1.0

    ! Allocate memory for local variables

    if(.not. allocated(mask))      allocate(mask(merge%nx,merge%ny))
    if(.not. allocated(mask_save)) allocate(mask_save(merge%nx,merge%ny))

    ! Define local variables
    
    mask      = reshape(merge%mask,shape(mask))
    mask_save = mask
    
    ! Loop through local variable

    do j = 1, merge%ny

       ! Define local variables

       idx_strt = imin
       idx_stop = imax

       ! Loop through local variable

       do i = 1, merge%nx

          ! Check local variable and proceed accordingly
          
          if(mask(i,j) .eq. spval) idx_strt = min(idx_strt,i)
          if(mask(i,j) .eq. 1.0)   idx_stop = max(idx_stop,i)
          
       end do ! do i = 1, merge%nx
       
       ! Check local variable and proceed accordingly
       
       if((idx_stop .ne. imax) .and. (idx_strt .ne. imin) .and. (idx_strt  &
            & .ne. 1) .and. (idx_stop .ne. merge%nx)) then
       
          ! Define local variables
          
          icoord   = 1          
          strt_idx = min((idx_strt + merge%npts_relax),merge%nx)
          stop_idx = idx_strt 
          icoord   = strt_idx
             
          ! Loop through local variable

          do while (icoord .ge. stop_idx)
 
             ! Define local variables

             mask_save(stop_idx:icoord,j) = 1.0 - (real(icoord) -          &
                  & real(stop_idx))/(real(strt_idx) - real(stop_idx))
             icoord                       = icoord - 1

          end do ! do while (icoord .le. stop_idx)
          
       end if ! if((idx_stop .ne. imax) .and. (idx_strt .ne. imin)
              ! .and. (idx_strt .ne. 1) .and. (idx_stop
              ! .ne. merge%nx))

    end do ! do j = 1, merge%ny

    ! Loop through local variable

    do j = 1, merge%ny
    
       ! Define local variables

       idx_strt = imax
       idx_stop = imin

       ! Loop through local variable

       do i = 1, merge%nx

          ! Check local variable and proceed accordingly
          
          if(mask(i,j) .eq. spval) idx_strt = max(idx_strt,i)
          if(mask(i,j) .eq. 1.0)   idx_stop = min(idx_stop,i)
          
       end do ! do i = 1, merge%nx

       ! Check local variable and proceed accordingly
       
       if((idx_stop .ne. imin) .and. (idx_strt .ne. imax) .and. (idx_strt  &
            & .ne. merge%nx) .and. (idx_stop .ne. 1)) then
          
          ! Define local variables
          
          icoord   = 1          
          strt_idx = max((idx_stop - merge%npts_relax),1)
          stop_idx = idx_stop
          icoord   = strt_idx
             
          ! Loop through local variable

          do while (icoord .le. stop_idx)
 
             ! Define local variables

             mask_save(icoord:stop_idx,j) = 1.0 - (real(icoord) -          &
                  & real(stop_idx))/(real(strt_idx) - real(stop_idx))
             icoord                       = icoord + 1

          end do ! do while (icoord .le. stop_idx)
          
       end if ! if((idx_stop .ne. imin) .and. (idx_strt .ne. imax)
              ! .and. (idx_strt .ne. merge%nx) .and. (idx_stop
              ! .ne. 1))
       
    end do ! do j = 1, merge%ny

    ! Loop through local variable

    do i = 1, merge%nx

       ! Define local variables

       idx_strt = imin
       idx_stop = imax

       ! Loop through local variable

       do j = 1, merge%ny

          ! Check local variable and proceed accordingly
          
          if(mask(i,j) .eq. spval) idx_strt = min(idx_strt,j)
          if(mask(i,j) .eq. 1.0)   idx_stop = max(idx_stop,j)
          
       end do ! do j = 1, merge%ny

       ! Check local variable and proceed accordingly

       if((idx_stop .ne. imax) .and. (idx_strt .ne. imin) .and. (idx_stop   &
            & .ne. merge%ny) .and. (idx_strt .ne. 1)) then
          
          ! Define local variables
          
          icoord   = 1
          strt_idx = min((idx_strt + merge%npts_relax),merge%ny)
          stop_idx = idx_strt
          icoord   = strt_idx
             
          ! Loop through local variable

          do while (icoord .ge. stop_idx)
 
             ! Define local variables
                
             mask_save(i,stop_idx:icoord) = 1.0 - (real(stop_idx) -         &
                  & real(icoord))/(real(stop_idx) - real(strt_idx))
             icoord                       = icoord - 1

          end do ! do while (icoord .le. stop_idx)

       end if ! if((idx_stop .ne. imax) .and. (idx_strt .ne. imin)
              ! .and. (idx_stop .ne. merge%ny) .and. (idx_strt
              ! .ne. 1))

    end do ! do i = 1, merge%nx

    ! Loop through local variable

    do i = 1, merge%nx

       ! Define local variables

       idx_strt = imax
       idx_stop = imin

       ! Loop through local variable

       do j = 1, merge%ny

          ! Check local variable and proceed accordingly
          
          if(mask(i,j) .eq. spval) idx_strt = max(idx_strt,j)
          if(mask(i,j) .eq. 1.0)   idx_stop = min(idx_stop,j)
          
       end do ! do j = 1, merge%ny

       ! Check local variable and proceed accordingly

       if((idx_stop .ne. imax) .and. (idx_strt .ne. imin) .and. (idx_stop   &
            & .ne. 1) .and. (idx_strt .ne. merge%ny)) then
          
          ! Define local variables
          
          icoord = 1
          strt_idx = max((idx_stop - merge%npts_relax),1)
          stop_idx = idx_stop 
          icoord   = strt_idx
             
          ! Loop through local variable

          do while (icoord .le. stop_idx)
             
             ! Define local variables
                
             mask_save(i,icoord:stop_idx) = 1.0 - (real(icoord) -           &
                     & real(stop_idx))/(real(strt_idx) - real(stop_idx))
             icoord                       = icoord + 1

          end do ! do while (icoord .le. stop_idx)

       end if ! if((idx_stop .ne. imax) .and. (idx_strt .ne. imin)
              ! .and. (idx_stop .ne. 1) .and. (idx_strt
              ! .ne. merge%ny))

    end do ! do i = 1, merge%nx

    ! Define local variables

    merge%mask = reshape(mask_save,shape(merge%mask))
    where(merge%mask .eq. spval) merge%mask = 0.0
    where(merge%mask .gt. 1.0)   merge%mask = 1.0
    where(merge%mask .lt. 0.0)   merge%mask = 0.0
    
    ! Deallocate memory for local variables

    if(allocated(mask))      deallocate(mask)
    if(allocated(mask_save)) deallocate(mask_save)
    call variable_interface_cleanup_struct(grid_an)
    call variable_interface_cleanup_struct(grid_fg)
    call variable_interface_cleanup_struct(kdtree)

    !=====================================================================
    
  end subroutine define_overlap_mask

 !=======================================================================

  ! SUBROUTINE:

  ! environment_blend.f90

  ! DESCRIPTION:

  ! This subroutine blends values to fill a missing data region
  ! (signified by a value of 'spval') within the input variable grid
  ! (in_grid) and fills the respective locations (determined via
  ! interpolation) using the analysis and first-guess variable arrays
  ! (within the respective an_grid and fg_grid FORTRAN grid_struct
  ! variables) to relax the regional values into the environmental
  ! values.

  ! INPUT VARIABLES:

  ! * dst_grid; a FORTRAN grid_struct variable containing a region to
  !   be populated via linear-interpolation.

  ! * grid_an; a FORTRAN grid_struct variable containing the analysis
  !   variable attributes.

  ! * grid_fg; a FORTRAN grid_struct variable containing the
  !   first-guess variable attributes.

  ! OUTPUT VARIABLES:

  ! * an_grid; a FORTRAN grid_struct variable containing analysis
  !   variable values where the input grid values are determined via
  !   linear-interpolation.

  !-----------------------------------------------------------------------

  subroutine environment_blend(dst_grid,grid_an,grid_fg)

    ! Define variables passed to routine

    type(grid_struct)                                                   :: dst_grid
    type(grid_struct)                                                   :: grid_an
    type(grid_struct)                                                   :: grid_fg

    ! Define variables computed within routine

    real(r_kind)                                                        :: mask_relax

    ! Define counting variables                  

    integer                                                             :: i, j

    !=====================================================================

    ! Loop through local variable

    do j = 1, grid_an%nz

       ! Define local variables

       where(dst_grid%radius .le. tc_ic_radius) dst_grid%var(:,j) =        &
            & grid_an%var(:,j)
       where(dst_grid%radius .ge. tc_large_scale_radius)                   &
            & dst_grid%var(:,j) = grid_fg%var(:,j)
       
       ! Loop through local variable

       do i = 1, grid_an%ncoords
       
          ! Check local variable and proceed accordingly

          if((dst_grid%radius(i) .gt. tc_ic_radius) .and.                  &
               & (dst_grid%radius(i) .lt. tc_large_scale_radius)) then

             ! Define local variables

             mask_relax        = (dst_grid%radius(i) - tc_ic_radius)/      &
                  & (tc_large_scale_radius - tc_ic_radius)
             dst_grid%var(i,j) = grid_an%var(i,j)*(1.0 - mask_relax) +     &
                  & grid_fg%var(i,j)*mask_relax

          end if ! if((dst_grid%radius(i) .gt. tc_ic_radius)
                 ! .and. (dst_grid%radius(i)
                 ! .lt. tc_large_scale_radius)

       end do ! do i = 1, grid_an%ncoords

    end do ! do j = 1, grid_an%nz

    !=====================================================================

  end subroutine environment_blend
  
  !=======================================================================

  ! SUBROUTINE:

  ! grid_blend.f90

  ! DESCRIPTION:

  ! This subroutine defines the remapping attributes required to blend
  ! two grids together.

  ! INPUT VARIABLES:

  ! * blend; a FORTRAN blend_struct variable; it is assumed that all
  !   arrays within the structure and the respective blending
  !   attributes have been specified before calling this routine.

  ! * slint; a FORTRAN slint_struct variable; only utilized if SLINT
  !   interpolation is specified by the user.

  ! * mask; a FORTRAN 4-byte real variable specifying the mask value
  !   for which to interpolate.

  ! OUTPUT VARIABLES:

  ! * slint; a FORTRAN slint_struct variable containing the SLINT
  !   defined remapping attributes if SLINT interpolation is specified
  !   by the user.

  !-----------------------------------------------------------------------

  subroutine grid_blend(blend,slint,mask)

    ! Define variables passed to routine

    type(blend_struct)                                                  :: blend
    type(slint_struct)                                                  :: slint
    real(r_kind)                                                        :: mask

    ! Define variables computed within routine

    type(grid_struct)                                                   :: dst_grid
    type(grid_struct)                                                   :: src_grid
    integer                                                             :: dst_ncoords
    integer                                                             :: src_ncoords

    ! Define counting variables

    integer                                                             :: i

    !=====================================================================

    ! Define local variables

    src_grid%ncoords = count(blend%mask .eq. spval)
    src_grid%nz      = 1
    call variable_interface_setup_struct(src_grid)
    dst_grid%ncoords = count(blend%mask .eq. mask)
    dst_grid%nz      = 1
    call variable_interface_setup_struct(dst_grid)
    src_ncoords      = 0
    dst_ncoords      = 0

    ! Allocate memory for local variables

    if(.not. allocated(blend%var))                                        &
         & allocate(blend%var(src_grid%ncoords,src_grid%nz))

    ! Loop through local variable

    do i = 1, blend%ncoords

       ! Check local variable and proceed accordingly

       if(blend%mask(i) .eq. spval) then

          ! Define local variables

          src_ncoords               = src_ncoords + 1
          src_grid%lat(src_ncoords) = blend%lat(i)
          src_grid%lon(src_ncoords) = blend%lon(i)
          blend%var(src_ncoords,:)  = blend%var_blend(i,:)

       end if ! if(blend%mask(i) .eq. spval)

       ! Check local variable and proceed accordingly

       if(blend%mask(i) .eq. mask) then

          ! Define local variables

          dst_ncoords               = dst_ncoords + 1
          dst_grid%lat(dst_ncoords) = blend%lat(i)
          dst_grid%lon(dst_ncoords) = blend%lon(i)
          blend%idx(i)              = dst_ncoords

       end if ! if(blend%mask(i) .eq. mask)

    end do ! do i = 1, blend%ncoords

    ! Check local variable and proceed accordingly

    if(is_slint) then

       ! Define local variables

       slint%ncoords = dst_ncoords
       call variable_interface_setup_struct(slint)

       ! Compute local variables
       
       call interpolation_interface_init(src_grid,dst_grid,slint)
       
    end if ! if(is_slint) 
    
    ! Define local variables
    
    blend%scoords = src_ncoords
    
    ! Deallocate memory for local variables
    
    call variable_interface_cleanup_struct(src_grid)
    call variable_interface_cleanup_struct(dst_grid)

    !=====================================================================

  end subroutine grid_blend

  !=======================================================================

  ! SUBROUTINE:

  ! init_mask_merge.f90

  ! DESCRIPTION:

  ! This subroutine initializes the blending mask used to merge two
  ! identical dimension grid variables.

  ! INPUT VARIABLES:

  ! * merge; a FORTRAN merge_struct variable; it is assumed that all
  !   arrays within the structure and the respective blending
  !   attributes have been specified before calling this routine.

  ! OUTPUT VARIABLES:

  ! * merge; a FORTRAN merge_struct variable containing the merging
  !   and blending mask regions.

  !-----------------------------------------------------------------------

  subroutine init_mask_merge(merge)

    ! Define variables passed to routine

    type(merge_struct)                                                  :: merge

    ! Define variables computed within routine

    real(r_kind)                                                        :: mnlat_fg
    real(r_kind)                                                        :: mnlon_fg
    real(r_kind)                                                        :: mxlat_fg
    real(r_kind)                                                        :: mxlon_fg

    !=====================================================================

    ! Define local variables

    mnlat_fg   = minval(merge%lat_fg)
    mnlon_fg   = minval(merge%lon_fg)
    mxlat_fg   = maxval(merge%lat_fg)
    mxlon_fg   = maxval(merge%lon_fg)
    merge%mask = 1.0
    where(merge%lat_an .lt. mnlat_fg) merge%mask = -1.0
    where(merge%lon_an .lt. mnlon_fg) merge%mask = -1.0
    where(merge%lat_an .gt. mxlat_fg) merge%mask = -1.0
    where(merge%lon_an .gt. mxlon_fg) merge%mask = -1.0

    !=====================================================================

  end subroutine init_mask_merge
  
  !=======================================================================
  
  ! SUBROUTINE:

  ! interpolate.f90

  ! DESCRIPTION:

  ! This subroutine remaps a user specified variable using the user
  ! specified interpolation remapping scheme bi-linear method.

  ! INPUT VARIABLES:

  ! * nlevs; a FORTRAN integer value specifying the total number of
  !   vertical levels.

  ! * src_grid; a FORTRAN grid_struct variable containing the
  !   geographical locations which do not contain missing values
  !   (e.g., spval) and the associated variable values.

  ! * dst_grid; a FORTRAN grid_struct variable containing the
  !   geographical locations of the destination grid.

  ! OUTPUT VARIABLES:

  ! * dst_grid; a FORTRAN grid_struct variable containing the
  !   interpolated variable values within the structure variable
  !   'var'.

  !-----------------------------------------------------------------------

  subroutine interpolate(nlevs,src_grid,dst_grid)
    
    ! Define variables passed to routine

    integer                                                             :: nlevs
    type(grid_struct)                                                   :: src_grid
    type(grid_struct)                                                   :: dst_grid

    ! Define variables computed within routine

    type(slint_struct)                                                  :: slint

    ! Define counting variables

    integer                                                             :: i

    !=====================================================================

    ! Check local variable and proceed accordingly

    if(is_slint) then

       ! Define local variables
       
       slint%ncoords     = dst_grid%ncoords
       slint%src_ncoords = src_grid%ncoords
       call variable_interface_setup_struct(slint)

       ! Compute local variables
    
       call interpolation_interface_init(src_grid,dst_grid,slint)

    end if ! if(is_slint)

    ! Loop through local variable

    do i = 1, nlevs

       ! Check local variable and proceed accordingly

       if(is_slint) then

          ! Compute local variables

          call interpolation_interface_bilinear(slint,src_grid%var(:,i))

          ! Define local variables

          dst_grid%var(:,i) = slint%var

       end if ! if(is_slint)

    end do ! do i = 1, nlevs

    ! Deallocate memory for local variables

    call variable_interface_cleanup_struct(slint)

    !=====================================================================

  end subroutine interpolate

  !=======================================================================

  ! SUBROUTINE:

  ! interpolate_subregion.f90

  ! DESCRIPTION:

  ! This subroutine interpolates a sub-region (dst_sub_grid) of the
  ! destination geographical projection (dst_grid).

  ! INPUT VARIABLES:

  ! * src_grid; a FORTRAN grid_struct variable containing geographical
  !   locations and values to remapped to the sub-region grid
  !   projection.

  ! * dst_sub_grid; a FORTRAN grid_struct variable containing the
  !   geographical locations within the destination grid projection
  !   for remapping of the source geographical projection (src_grid)
  !   values.

  ! * dst_grid; a FORTRAN grid_struct variable to contain the source
  !   grid values remapped within the destination grid sub-region.

  ! OUTPUT VARIABLES:

  ! * dst_grid; a FORTRAN grid_struct variable containing the source
  !   grid values remapped within the destination grid sub-region.

  !-----------------------------------------------------------------------

  subroutine interpolate_subregion(src_grid,dst_sub_grid,dst_grid)

    ! Define variables passed to routine

    type(grid_struct)                                                   :: dst_grid
    type(grid_struct)                                                   :: dst_sub_grid
    type(grid_struct)                                                   :: src_grid

    ! Define variables computed within routine

    integer                                                             :: dst_idx

    ! Define counting variables

    integer                                                             :: i, j

    !=====================================================================

    ! Compute local variables

    call interpolate(src_grid%nz,src_grid,dst_sub_grid)

    ! Loop through local variable

    do j = 1, dst_grid%nz

       ! Loop through local variable

       do i = 1, dst_sub_grid%ncoords

          ! Define local variables

          dst_idx                 = dst_sub_grid%idx(i)
          dst_grid%var(dst_idx,j) = dst_sub_grid%var(i,j)

       end do ! do i = 1, dst_sub_grid%ncoords

    end do ! do j = 1, dst_grid%nz

    !=====================================================================

  end subroutine interpolate_subregion

  !=======================================================================

  ! SUBROUTINE:

  ! mask_fill.f90

  ! DESCRIPTION:

  ! This subroutine determines the edges of the masked region (e.g.,
  ! where the local values are greater than 0.0) and fills all values
  ! within that array interval with the mask value.

  ! INPUT VARIABLES:

  ! * blend; a FORTRAN blend_struct variable with the initial mask
  !   region (e.g., all mask values with a non-zero variable value
  !   difference are set to 1.0).

  ! OUTPUT VARIABLES:

  ! * blend; a FORTRAN blend_struct variable with a mask region that
  !   has been filled within the respective interval.

  !-----------------------------------------------------------------------

  subroutine mask_fill(blend)

    ! Define variables passed to routine

    type(blend_struct)                                                  :: blend

    ! Define variables computed within routine

    real(r_kind),               dimension(:,:),             allocatable :: grid_mask
    real(r_kind)                                                        :: mxval
    integer                                                             :: mnxidx
    integer                                                             :: mxxidx
    integer                                                             :: mnyidx
    integer                                                             :: mxyidx

    ! Define counting variables

    integer                                                             :: i, j

    !=====================================================================

    ! Allocate memory for local variables

    if(.not. allocated(grid_mask)) allocate(grid_mask(blend%nx,blend%ny))

    ! Define local variables

    grid_mask = reshape(blend%mask,shape(grid_mask))

    ! Loop through local variable

    do j = 1, blend%ny

       ! Define local variables

       mxval = maxval(grid_mask(:,j))

       ! Check local variable and proceed accordingly

       if(mxval .gt. 0.0) then

          ! Define local variables

          mnxidx = blend%nx
          mxxidx = 1

          ! Loop through local variable

          do i = 1, blend%nx

             ! Check local variable and proceed accordingly

             if(grid_mask(i,j) .gt. 0.0) then

                ! Define local variables

                mnxidx = min(mnxidx,i)
                mxxidx = max(mxxidx,i)

             end if ! if(grid_mask(i,j) .gt. 0.0)

          end do ! do i = 1, blend%nx

          ! Define local variables

          grid_mask(mnxidx:mxxidx,j) = 1.0

       end if ! if(mxval .gt. 0.0)

    end do ! do j = 1, blend%ny

    ! Loop through local variable

    do i = 1, blend%nx

       ! Define local variables

       mxval = maxval(grid_mask(i,:))

       ! Check local variable and proceed accordingly

       if(mxval .gt. 0.0) then

          ! Define local variables

          mnyidx = blend%ny
          mxyidx = 1

          ! Loop through local variable

          do j = 1, blend%ny

             ! Check local variable and proceed accordingly

             if(grid_mask(i,j) .gt. 0.0) then

                ! Define local variables

                mnyidx = min(mnyidx,j)
                mxyidx = max(mxyidx,j)

             end if ! if(grid_mask(i,j) .gt. 0.0)

          end do ! do j = 1, blend%ny

          ! Define local variables

          grid_mask(i,mnyidx:mxyidx) = 1.0
          
       end if ! if(mxval .gt. 0.0)

    end do ! do i = 1, blend%nx

    ! Define local variables

    blend%mask = reshape(grid_mask,shape(blend%mask))

    ! Deallocate memory for local variables

    if(allocated(grid_mask)) deallocate(grid_mask)

    !=====================================================================

  end subroutine mask_fill

  !=======================================================================

  ! SUBROUTINE:

  ! mask_relax_boundary.f90

  ! DESCRIPTION:

  ! This subroutine defines the blending (i.e., relaxation) region for
  ! two variables of identical dimension and identical grid projectio;
  ! the size of the blending region is specified by the user and
  ! carried in the FORTRAN blend_struct variable 'npts_blend'
  ! attribute.

  ! INPUT VARIABLES:

  ! * blend; a FORTRAN blend_struct variable with a mask region that
  !   has been filled within the respective interval and contains the
  !   attribute 'npts_blend'.

  ! OUTPUT VARIABLES:

  ! * blend; a FORTRAN blend_struct variable within a mask region that
  !   now contains the relaxation region also; the relaxation region
  !   is specified by the 'spval' attribute while non-mask regions are
  !   denoted by a value of 0.0 within the mask.

  !-----------------------------------------------------------------------

  subroutine mask_relax_boundary(blend)

    ! Define variables passed to routine

    type(blend_struct)                                                  :: blend  

    ! Define variables computed within routine

    real(r_kind),               dimension(:,:),             allocatable :: grid_mask
    real(r_kind)                                                        :: grid_value
    real(r_kind)                                                        :: mxval
    integer                                                             :: mnxidx
    integer                                                             :: mxxidx
    integer                                                             :: mnyidx
    integer                                                             :: mxyidx
    integer                                                             :: mnx
    integer                                                             :: mxx
    integer                                                             :: mny
    integer                                                             :: mxy
    
    ! Define counting variables

    integer                                                             :: i, j, k

    !=====================================================================

    ! Allocate memory for local variables

    if(.not. allocated(grid_mask)) allocate(grid_mask(blend%nx,blend%ny))

    ! Define local variables

    grid_mask = reshape(blend%mask,shape(grid_mask))
    mnx       = imin
    mxx       = imax
    mny       = imin
    mxy       = imax
    
    ! Loop through local variable

    do j = blend%bmny, blend%bmxy

       ! Loop through local variable

       do i = blend%bmnx, blend%rmnx

          ! Define local variables

          grid_value     = (i - blend%bmnx) + 1
          grid_mask(i,j) = grid_value

       end do ! do i = blend%bmnx, blend%rmnx

       ! Loop through local variable

       do i = blend%rmxx, blend%bmxx

          ! Define local variables

          grid_value     = (blend%bmxx - i) + 1
          grid_mask(i,j) = grid_value

       end do ! do i = blend%bmnx, blend%rmnx

    end do ! do j = 1, blend%ny

    ! Loop through local variable

    do i = blend%bmnx, blend%bmxx

       ! Loop through local variable

       do j = blend%bmny, blend%rmny

          ! Define local variables

          grid_value     = (j - blend%bmny) + 1
          grid_mask(i,j) = grid_value

       end do ! do j = blend%bmny, blend%rmny

       ! Loop through local variable

       do j = blend%rmxy, blend%bmxy

          ! Define local variables

          grid_value     = (blend%bmxy - j) + 1
          grid_mask(i,j) = grid_value

       end do ! do j = blend%bmny, blend%rmny

    end do ! do i = blend%bmnx, blend%rmnx

    ! Define local variables

    where(grid_mask .le. 1.0) grid_mask = spval
    blend%mask = reshape(grid_mask,shape(blend%mask))

    ! Deallocate memory for local variables

    if(allocated(grid_mask)) deallocate(grid_mask)

    !=====================================================================

  end subroutine mask_relax_boundary
  
  !=======================================================================

  ! SUBROUTINE:

  ! mask_relax_merge.f90

  ! DESCRIPTION:

  ! This subroutine defines the blending (i.e., relaxation) region
  ! surrounding a first-guess grid relative to the analysis grid
  ! within which it is to be inserted (i.e., merged); the size of the
  ! blending region is specified by the user and carried in the
  ! FORTRAN merge_struct variable 'npts_blend' and 'npts_relax' attributes.

  ! INPUT VARIABLES:

  ! * merge; a FORTRAN merge_struct variable with a mask region that
  !   has been filled within the respective interval and contains the
  !   attributes 'npts_blend' and 'npts_relax'.

  ! OUTPUT VARIABLES:

  ! * merge; a FORTRAN merge_struct variable within a mask region that
  !   now contains the relaxation region also; the relaxation region
  !   is specified by the 'spval' attribute while non-mask regions are
  !   denoted by a value of 0.0 within the mask.

  !-----------------------------------------------------------------------

  subroutine mask_relax_merge(merge)

    ! Define variables passed to routine

    type(merge_struct)                                                  :: merge  

    ! Define variables computed within routine

    type(blend_struct)                                                  :: blend
    real(r_kind),               dimension(:,:),             allocatable :: grid_mask
    real(r_kind),               dimension(:,:),             allocatable :: workgrid
    
    ! Define counting variables

    integer                                                             :: i, j, k

    !=====================================================================

    ! Allocate memory for local variables

    if(.not. allocated(grid_mask)) allocate(grid_mask(merge%nx,merge%ny))
    if(.not. allocated(workgrid))  allocate(workgrid(merge%nx,merge%ny))

    ! Define local variables

    grid_mask        = reshape(merge%mask,shape(grid_mask))
    workgrid         = grid_mask
    where(workgrid .lt. 0.0) workgrid = spval
    blend%nx         = merge%nx
    blend%ny         = merge%ny
    blend%nz         = merge%nz
    call variable_interface_setup_struct(blend)
    blend%npts_relax = merge%npts_relax
    blend%npts_blend = merge%npts_blend
    blend%mask       = merge%mask
    call mask_relax_nest(blend)
    grid_mask        = reshape(blend%mask,(/merge%nx,merge%ny/))
  
    ! Loop through local variable

    do i = 1, merge%nx

       ! Check local variable and proceed accordingly

       if(maxval(workgrid(i,1:blend%mny)) .ne. spval) then

          ! Define local variables

          grid_mask(i,1:blend%mny) = spval

       end if ! if(maxval(workgrid(i,1:blend%mny)) .ne. spval)

       ! Check local variable and proceed accordingly

       if(maxval(workgrid(i,blend%mxy:blend%nx)) .ne. spval) then

          ! Define local variables

          grid_mask(i,blend%mxy:blend%nx) = spval

       end if ! if(maxval(workgrid(i,blend%mxy:blend%nx)) .ne. spval)

    end do ! do i = 1, merge%nx

    ! Loop through local variable

    do j = 1, merge%ny

       ! Check local variable and proceed accordingly

       if(maxval(workgrid(1:blend%mnx,j)) .ne. spval) then

          ! Define local variables

          grid_mask(1:blend%mnx,j) = spval

       end if ! if(maxval(workgrid(1:blend%mnx,j)) .ne. spval)

       ! Check local variable and proceed accordingly

       if(maxval(workgrid(blend%mxx:blend%nx,j)) .ne. spval) then

          ! Define local variables

          grid_mask(blend%mxx:blend%nx,j) = spval

       end if ! if(maxval(workgrid(blend%mxx:blend%nx,j)) .ne. spval)

    end do ! do j = 1, merge%ny

    ! Define local variables

    merge%mask = reshape(grid_mask,shape(merge%mask))

    ! Deallocate memory for local variables

    call variable_interface_cleanup_struct(blend)
    if(allocated(grid_mask)) deallocate(grid_mask)
    if(allocated(workgrid))  deallocate(workgrid)

    !=====================================================================

  end subroutine mask_relax_merge

  !=======================================================================

  ! SUBROUTINE:

  ! mask_relax_nest.f90

  ! DESCRIPTION:

  ! This subroutine defines the blending (i.e., relaxation) region
  ! surrounding the nest and it's location within it's respective
  ! parent domain; the size of the blending region is specified by the
  ! user and carried in the FORTRAN blend_struct variable 'npts_blend'
  ! attribute.

  ! INPUT VARIABLES:

  ! * blend; a FORTRAN blend_struct variable with a mask region that
  !   has been filled within the respective interval and contains the
  !   attribute 'npts_blend'.

  ! OUTPUT VARIABLES:

  ! * blend; a FORTRAN blend_struct variable within a mask region that
  !   now contains the relaxation region also; the relaxation region
  !   is specified by the 'spval' attribute while non-mask regions are
  !   denoted by a value of 0.0 within the mask.

  !-----------------------------------------------------------------------

  subroutine mask_relax_nest(blend)

    ! Define variables passed to routine

    type(blend_struct)                                                  :: blend  

    ! Define variables computed within routine

    real(r_kind),               dimension(:,:),             allocatable :: grid_mask
    real(r_kind)                                                        :: mxval
    integer                                                             :: mnxidx
    integer                                                             :: mxxidx
    integer                                                             :: mnyidx
    integer                                                             :: mxyidx
    integer                                                             :: mnx
    integer                                                             :: mxx
    integer                                                             :: mny
    integer                                                             :: mxy
    
    ! Define counting variables

    integer                                                             :: i, j, k

    !=====================================================================

    ! Allocate memory for local variables

    if(.not. allocated(grid_mask)) allocate(grid_mask(blend%nx,blend%ny))

    ! Define local variables

    grid_mask = reshape(blend%mask,shape(grid_mask))
    mnx       = imin
    mxx       = imax
    mny       = imin
    mxy       = imax

    ! Loop through local variable

    do j = 1, blend%ny

       ! Define local variables

       mxval = maxval(grid_mask(:,j))

       ! Check local variable and proceed accordingly

       if(mxval .gt. 0.0) then

          ! Define local variables

          mnxidx = blend%nx
          mxxidx = 1

          ! Loop through local variable

          do i = 1, blend%nx

             ! Check local variable and proceed accordingly

             if(grid_mask(i,j) .gt. 0.0) then

                ! Define local variables

                mnxidx = min(mnxidx,(i + blend%npts_relax))
                mxxidx = max(mxxidx,(i - blend%npts_relax))

             end if ! if(grid_mask(i,j) .gt. 0.0)

          end do ! do i = 1, blend%nx

          ! Loop through local variable

          do i = (blend%npts_blend + blend%npts_relax), 1, -1

             ! Define local variables

             mnx                     = max(1,(mnxidx - i))
             mxx                     = min(blend%nx,(mxxidx + i))
             grid_mask(mnx:mnxidx,j) = real(i)
             grid_mask(mxxidx:mxx,j) = real(i)

          end do ! do i = (blend%npts_blend + blend%npts_relax), 1, -1

       end if ! if(mxval .gt. 0.0)

    end do ! do j = 1, blend%ny

    ! Loop through local variable

    do i = 1, blend%nx

       ! Define local variables

       mxval = maxval(grid_mask(i,:))

       ! Check local variable and proceed accordingly

       if(mxval .gt. 0.0) then

          ! Define local variables

          mnyidx = blend%ny
          mxyidx = 1

          ! Loop through local variable

          do j = 1, blend%ny

             ! Check local variable and proceed accordingly

             if((grid_mask(i,j) .gt. 0.0)) then

                ! Define local variables

                mnyidx = min(mnyidx,(j + blend%npts_relax))
                mxyidx = max(mxyidx,(j - blend%npts_relax))

             end if ! if((grid_mask(i,j) .gt. 0.0))

          end do ! do j = 1, blend%ny

          ! Loop through local variable

          do j = (blend%npts_blend + blend%npts_relax), 1, -1

             ! Define local variables

             mny                     = max(1,(mnyidx - j))
             mxy                     = min(blend%ny,(mxyidx + j))
             grid_mask(i,mny:mnyidx) = real(j)
             grid_mask(i,mxyidx:mxy) = real(j)

          end do ! do j = (blend%npts_blend + blend%npts_relax), 1, -1

       end if ! if(mxval .gt. 0.0)

    end do ! do i = 1, blend%nx

    ! Define local variables

    blend%mnx  = mnx
    blend%mny  = mny
    blend%mxx  = mxx
    blend%mxy  = mxy
    where(grid_mask .le. 1.0) grid_mask = spval
    blend%mask = reshape(grid_mask,shape(blend%mask))

    ! Deallocate memory for local variables

    if(allocated(grid_mask)) deallocate(grid_mask)

    !=====================================================================

  end subroutine mask_relax_nest

  !=======================================================================

  ! SUBROUTINE:

  ! region_blend.f90

  ! DESCRIPTION:

  ! This subroutine blends values to fill a missing data region
  ! (signified by a value of 'spval') within the input variable grid
  ! (in_grid) and fills the respective locations (determined via
  ! interpolation) within the output grid (out_grid).

  ! INPUT VARIABLES:

  ! * in_grid; a FORTRAN grid_struct variable containing a region of
  !   missing data values to be populated via interpolation.

  ! * out_grid; a FORTRAN grid_struct variable containing the
  !   geographical locations for the destination grid.

  ! OUTPUT VARIABLES:

  ! * out_grid; a FORTRAN grid_struct variable containing variable
  !   values where the input grid once contained missing data value;
  !   these values are determined via interpolation.

  !-----------------------------------------------------------------------

  subroutine region_blend(in_grid,out_grid)

    ! Define variables passed to routine

    type(grid_struct)                                                   :: in_grid
    type(grid_struct)                                                   :: out_grid

    ! Define variables computed within routine

    real(r_kind)                                                        :: mask_relax

    ! Define counting variables                  

    integer                                                             :: i, j

    !=====================================================================
    
    ! Loop through local variable

    do i = 1, out_grid%ncoords

       ! Check local variable and proceed accordingly

       if(out_grid%radius(i) .le. tc_ic_radius) then
          
          ! Loop through local variable

          do j = 1, out_grid%nz
          
             ! Define local variables

             out_grid%var(i,j) = in_grid%var(i,j)

          end do ! do j = 1, out_grid%nz
             
       end if ! if(out_grid%radius(i) .le. tc_ic_radius)

       ! Check local variable and proceed accordingly

       if((out_grid%radius(i) .gt. tc_ic_radius) .and.                     &
            & (out_grid%radius(i) .lt. tc_env_radius)) then

          ! Compute local variables

          mask_relax = (in_grid%radius(i) - tc_ic_radius)/(tc_env_radius   &
               & - tc_ic_radius)

          ! Loop through local variable

          do j = 1, out_grid%nz          

             ! Compute local variables

             out_grid%var(i,j) = in_grid%var(i,j)*(1.0 - mask_relax) +     &
                  & out_grid%var(i,j)*mask_relax

          end do ! do j = 1, out_grid%nz

       end if ! if((out_grid%radius(i) .gt. tc_ic_radius)
              ! .and. (out_grid%radius(i) .lt. tc_env_radius))

    end do ! do i = 1, out_grid%ncoords
    
    !=====================================================================

  end subroutine region_blend

  !=======================================================================

  ! SUBROUTINE:

  ! wnd2d_tc_increment.f90

  ! DESCRIPTION:

  ! This subroutine computes the wave-number spectra for the
  ! respective variable using the Fast-Fourier Transform (FFT); the
  ! variable is first recentered relative to the geographical location
  ! of the TC; the 2-dimensional FFT is computed along each azimuth
  ! for all radii within the user specified radius of influence; the
  ! spectral components to be retained by the user are than projected
  ! back to a geographical map projection and returned via a FORTRAN
  ! grid_struct variable (dst_grid).

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN grid_struct variable containing the analysis
  !   increment and the corresponding grid-projection attributes.

  ! * tcv; a FORTRAN tcv_struct variable containing the geographical
  !   location for the respective TC about which to center the FFT.

  ! * dst_grid; a FORTRAN grid_struct variable to contain the
  !   wave-number spectra remapped to the geographical projection of
  !   the input grid (grid).

  ! OUTPUT VARIABLES:

  ! * dst_grid; a FORTRAN grid_struct variable to containing the
  !   wave-number spectra remapped to the geographical projection of
  !   the input grid (grid).

  !-----------------------------------------------------------------------

  subroutine wnd2d_tc_increment(grid,tcv,dst_grid)

    ! Define variables passed to routine

    type(grid_struct)                                                   :: dst_grid
    type(grid_struct)                                                   :: grid
    type(tcv_struct)                                                    :: tcv

    ! Define variables computed within routine

    type(grid_struct)                                                   :: dst_sub_grid
    type(grid_struct)                                                   :: src_grid
    type(recenter_struct)                                               :: recenter
    type(wnd2d_struct)                                                  :: wnd2d
    integer                                                             :: dst_sub_cnt

    ! Define counting variables

    integer                                                             :: i, j

    !=====================================================================

    ! Define local variables

    recenter%clon  = tcv%lon
    recenter%clat  = tcv%lat
    recenter%area  = wndcmp_area
    recenter%darea = wndcmp_darea
    recenter%nz    = grid%nz
    call interpolation_interface_init(recenter)
    
    ! Compute local variables

    call interpolation_interface_bilinear(grid,recenter)

    ! Define local variables

    wnd2d%nx        = recenter%nx
    wnd2d%ny        = recenter%ny
    wnd2d%maxradius = maxval(recenter%radius)
    wnd2d%dradius   = wndcmp_dradius
    wnd2d%dangle    = wndcmp_dangle
    call variable_interface_setup_struct(wnd2d)

    ! Loop through local variable

    do j = 1, recenter%nz

       ! Define local variables

       wnd2d%var = recenter%var(:,j)

       ! Compute local variables

       call math_methods_wnd(wnd2d,recenter,src_grid)

       ! Define local variables
       
       src_grid%var(:,j) = 0.0

       ! Loop through local variable

       do i = 1, min((wndcmp_mxwvn + 1),wnd2d%nh)

          ! Compute local variables
          
          src_grid%var(:,j) = src_grid%var(:,j) +                          &
               & reshape(wnd2d%wnvar(i,:,:),shape(src_grid%var(:,j)))

       end do ! do i = 1, min((wndcmp_mxwvn + 1),wnd2d%nh)

    end do ! do j = 1, recenter%nz

    ! Define local variables

    dst_sub_cnt          = count((grid%lat .ge. minval(src_grid%lat))      &
         & .and. (grid%lat .le. maxval(src_grid%lat)) .and. (grid%lon      &
         & .ge. minval(src_grid%lon)) .and. (grid%lon .le.                 &
         & maxval(src_grid%lon)))
    dst_sub_grid%ncoords = dst_sub_cnt
    dst_sub_grid%nz      = grid%nz
    call variable_interface_setup_struct(dst_sub_grid)
    
    ! Define local variables

    dst_sub_cnt = 0

    ! Loop through local variable

    do i = 1, grid%ncoords

       ! Check local variable and proceed accordingly

       if((grid%lat(i) .ge. minval(src_grid%lat)) .and. (grid%lat(i) .le.  &
            & maxval(src_grid%lat)) .and. (grid%lon(i) .ge.                &
            & minval(src_grid%lon)) .and. (grid%lon(i) .le.                &
            & maxval(src_grid%lon))) then
    
          ! Define local variables

          dst_sub_cnt                   = dst_sub_cnt + 1
          dst_sub_grid%lat(dst_sub_cnt) = grid%lat(i)
          dst_sub_grid%lon(dst_sub_cnt) = grid%lon(i)
          dst_sub_grid%idx(dst_sub_cnt) = i

       end if ! if((grid%lat(i) .ge. minval(src_grid%lat))
              ! .and. (grid%lat(i) .le. maxval(src_grid%lat))
              ! .and. (grid%lon(i) .ge. minval(src_grid%lon))
              ! .and. (grid%lon(i) .le. maxval(src_grid%lon)))

    end do ! do i = 1, grid%ncoords

    ! Compute local variables

    call interpolate_subregion(src_grid,dst_sub_grid,dst_grid)

    ! Deallocate memory for local variables
    
    call variable_interface_cleanup_struct(src_grid)
    call variable_interface_cleanup_struct(dst_sub_grid)

    !=====================================================================

  end subroutine wnd2d_tc_increment

  !=======================================================================

  ! SUBROUTINE: 

  ! zero_tc_increment.f90

  ! DESCRIPTION:

  ! This subroutine sets all analysis increment values within the user
  ! specified inner-core region of the tropical cyclone (TC) to zero.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN grid_struct variable containing the analysis
  !   increment and the corresponding grid-projection attributes.

  ! * dst_grid; a FORTRAN grid_struct variable to contain the updated
  !   analysis increments.

  ! OUTPUT VARIABLES:

  ! * dst_grid; a FORTRAN grid_struct variable to containing the
  !   updated analysis increments.

  !-----------------------------------------------------------------------

  subroutine zero_tc_increment(grid,dst_grid)

    ! Define variables passed to routine

    type(grid_struct)                                                   :: dst_grid
    type(grid_struct)                                                   :: grid

    ! Define counting variables

    integer                                                             :: i

    !=====================================================================

    ! Define local variables

    dst_grid%var = grid%var

    ! Loop through local variable

    do i = 1, dst_grid%nz
    
       ! Define local variables

       where(dst_grid%radius .le. tc_ic_radius) dst_grid%var(:,i) = 0.0

    end do ! do i = 1, dst_grid%nz

    !=====================================================================

  end subroutine zero_tc_increment

  !=======================================================================

end module analysis_interface
