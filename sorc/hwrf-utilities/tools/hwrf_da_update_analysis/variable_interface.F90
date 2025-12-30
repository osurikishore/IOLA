module variable_interface

  !=======================================================================

  !$$$ MODULE DOCUMENTATION BLOCK
  
  ! da-update-analysis :: variable_interface
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

  use kinds_interface
  use namelist_interface

  ! Define interfaces and attributes for module routines
  
  implicit none
  private
  public :: blend_struct
  public :: clpval
  public :: fft2d_struct
  public :: grid_struct
  public :: hwrf_bdy_grid_struct
  public :: hwrf_json_bdy_struct
  public :: hwrf_sfc_grid_struct
  public :: hwrf_vert_grid_struct
  public :: imax
  public :: imin
  public :: kdtree_struct
  public :: merge_struct
  public :: meteo_struct
  public :: moisture_struct
  public :: netcdf_grid_struct
  public :: netcdf_var_struct
  public :: nmm_json_grid_struct
  public :: nmm_json_var_struct
  public :: ntcs
  public :: presremap_struct
  public :: recenter_struct
  public :: smthr9pt_struct
  public :: spval
  public :: slint_struct
  public :: spline_struct
  public :: statgrid_struct
  public :: tcv_struct
  public :: timeinfo_struct
  public :: variable_interface_cleanup_struct
  public :: variable_interface_setup_struct
  public :: vargrid_struct
  public :: varinfo_struct
  public :: vertgrid_struct
  public :: wnd2d_struct
  interface variable_interface_cleanup_struct
     module procedure finalize_blend_struct
     module procedure finalize_fft2d_struct
     module procedure finalize_grid_struct
     module procedure finalize_hwrf_bdy_grid_struct
     module procedure finalize_hwrf_sfc_grid_struct
     module procedure finalize_hwrf_vert_grid_struct
     module procedure finalize_kdtree_struct
     module procedure finalize_merge_struct
     module procedure finalize_meteo_struct
     module procedure finalize_moisture_struct
     module procedure finalize_netcdf_grid_struct
     module procedure finalize_netcdf_var_struct
     module procedure finalize_presremap_struct
     module procedure finalize_recenter_struct
     module procedure finalize_slint_struct
     module procedure finalize_smthr9pt_struct
     module procedure finalize_spline_struct
     module procedure finalize_vargrid_struct
     module procedure finalize_varinfo_struct
     module procedure finalize_vertgrid_struct
     module procedure finalize_wnd2d_struct
  end interface variable_interface_cleanup_struct
  interface variable_interface_setup_struct
     module procedure initialize_blend_struct
     module procedure initialize_fft2d_struct
     module procedure initialize_grid_struct
     module procedure initialize_hwrf_bdy_grid_struct
     module procedure initialize_hwrf_sfc_grid_struct
     module procedure initialize_hwrf_vert_grid_struct
     module procedure initialize_kdtree_struct
     module procedure initialize_merge_struct
     module procedure initialize_meteo_struct
     module procedure initialize_moisture_struct
     module procedure initialize_netcdf_grid_struct
     module procedure initialize_netcdf_var_struct
     module procedure initialize_presremap_struct
     module procedure initialize_recenter_struct
     module procedure initialize_slint_struct
     module procedure initialize_smthr9pt_struct
     module procedure initialize_spline_struct
     module procedure initialize_vargrid_struct
     module procedure initialize_varinfo_struct
     module procedure initialize_vertgrid_struct
     module procedure initialize_wnd2d_struct
  end interface variable_interface_setup_struct
  
  ! Define local variables

  type blend_struct
     real(r_kind),              dimension(:,:),             allocatable :: var_fixed
     real(r_kind),              dimension(:,:),             allocatable :: var_blend
     real(r_kind),              dimension(:,:),             allocatable :: var
     real(r_kind),              dimension(:),               allocatable :: lat
     real(r_kind),              dimension(:),               allocatable :: lon
     real(r_kind),              dimension(:),               allocatable :: mask
     real(r_kind)                                                       :: clat
     real(r_kind)                                                       :: clon
     real(r_kind)                                                       :: tc_nest_max_radius
     real(r_kind)                                                       :: tc_nest_min_radius
     integer,                   dimension(:),               allocatable :: idx
     integer                                                            :: ncoords
     integer                                                            :: npts_blend
     integer                                                            :: npts_relax
     integer                                                            :: bmnx
     integer                                                            :: bmny
     integer                                                            :: bmxx
     integer                                                            :: bmxy
     integer                                                            :: mnx
     integer                                                            :: mny
     integer                                                            :: mxx
     integer                                                            :: mxy
     integer                                                            :: rmnx
     integer                                                            :: rmny
     integer                                                            :: rmxx
     integer                                                            :: rmxy
     integer                                                            :: nx
     integer                                                            :: ny
     integer                                                            :: nz
     integer                                                            :: scoords
  end type blend_struct           ! type blend_struct
  type fft2d_struct
     complex(r_double),         dimension(:,:),             allocatable :: out
     complex(r_double),         dimension(:,:),             allocatable :: in
     integer                                                            :: nx
     integer                                                            :: ny
  end type fft2d_struct           ! type fft2d_struct
  type grid_struct
     real(r_kind),              dimension(:,:),             allocatable :: var
     real(r_kind),              dimension(:),               allocatable :: angle
     real(r_kind),              dimension(:),               allocatable :: dx
     real(r_kind),              dimension(:),               allocatable :: dy
     real(r_kind),              dimension(:),               allocatable :: lat
     real(r_kind),              dimension(:),               allocatable :: lon
     real(r_kind),              dimension(:),               allocatable :: mapfac
     real(r_kind),              dimension(:),               allocatable :: radius
     real(r_kind)                                                       :: clat
     real(r_kind)                                                       :: clon
     real(r_kind)                                                       :: gchead
     real(r_kind)                                                       :: gcdist
     real(r_kind)                                                       :: gclat
     real(r_kind)                                                       :: gclon
     integer,                   dimension(:),               allocatable :: idx
     integer                                                            :: nx
     integer                                                            :: ny
     integer                                                            :: nz
     integer                                                            :: ncoords
  end type grid_struct            ! type grid_struct
  type hwrf_bdy_grid_struct
     real(r_kind),              dimension(:,:,:,:,:),       allocatable :: btxs
     real(r_kind),              dimension(:,:,:,:,:),       allocatable :: btxe
     real(r_kind),              dimension(:,:,:,:,:),       allocatable :: btye
     real(r_kind),              dimension(:,:,:,:,:),       allocatable :: btys
     real(r_kind),              dimension(:,:,:,:,:),       allocatable :: bxs
     real(r_kind),              dimension(:,:,:,:,:),       allocatable :: bxe
     real(r_kind),              dimension(:,:,:,:,:),       allocatable :: bye
     real(r_kind),              dimension(:,:,:,:,:),       allocatable :: bys
     real(r_kind),              dimension(:,:,:),           allocatable :: analvar
     real(r_kind)                                                       :: bdyfrq
     integer                                                            :: nb
     integer                                                            :: nt
     integer                                                            :: nx
     integer                                                            :: ny
     integer                                                            :: nz
  end type hwrf_bdy_grid_struct   ! type hwrf_bdy_grid_struct
  type hwrf_json_bdy_struct
     character(len=20)                                                  :: variable_name
     logical                                                            :: rezero
     integer                                                            :: gridtype
     integer                                                            :: levtype
  end type hwrf_json_bdy_struct   ! type hwrf_json_bdy_struct
  type hwrf_sfc_grid_struct
     real(r_kind),              dimension(:),               allocatable :: fis
     real(r_kind),              dimension(:),               allocatable :: sst
     real(r_kind),              dimension(:),               allocatable :: seaice
     integer                                                            :: nx
     integer                                                            :: ny
     integer                                                            :: ncoords     
  end type hwrf_sfc_grid_struct   ! type hwrf_sfc_grid_struct
  type hwrf_vert_grid_struct
     real(r_kind),              dimension(:,:),             allocatable :: pint
     real(r_kind),              dimension(:,:),             allocatable :: nest_pint
     real(r_kind),              dimension(:,:),             allocatable :: plev
     real(r_kind),              dimension(:,:),             allocatable :: nest_plev
     real(r_kind),              dimension(:,:),             allocatable :: var
     real(r_kind),              dimension(:),               allocatable :: eta1
     real(r_kind),              dimension(:),               allocatable :: eta2
     real(r_kind),              dimension(:),               allocatable :: aeta1
     real(r_kind),              dimension(:),               allocatable :: aeta2
     real(r_kind),              dimension(:),               allocatable :: pd
     real(r_kind),              dimension(:),               allocatable :: psfc
     real(r_kind)                                                       :: pdtop
     real(r_kind)                                                       :: pt
     integer                                                            :: nx
     integer                                                            :: ny
     integer                                                            :: nz
     integer                                                            :: ncoords     
  end type hwrf_vert_grid_struct  ! type hwrf_vert_grid_struct
  type kdtree_struct
     real(r_kind),              dimension(:,:),             allocatable :: r2dist
     real(r_kind)                                                       :: r2
     integer,                   dimension(:,:),             allocatable :: idx
     integer                                                            :: nalloc
     integer                                                            :: ncoords
     integer                                                            :: nfound
     integer                                                            :: nn
  end type kdtree_struct          ! type kdtree_struct  
  type merge_struct
     real(r_kind),              dimension(:,:),             allocatable :: var_an
     real(r_kind),              dimension(:,:),             allocatable :: var_fg
     real(r_kind),              dimension(:,:),             allocatable :: var
     real(r_kind),              dimension(:),               allocatable :: dx
     real(r_kind),              dimension(:),               allocatable :: dy
     real(r_kind),              dimension(:),               allocatable :: lat_an
     real(r_kind),              dimension(:),               allocatable :: lon_an
     real(r_kind),              dimension(:),               allocatable :: lat_fg
     real(r_kind),              dimension(:),               allocatable :: lon_fg
     real(r_kind),              dimension(:),               allocatable :: mask
     real(r_kind),              dimension(:),               allocatable :: radius
     integer,                   dimension(:),               allocatable :: idx
     integer                                                            :: bmnx
     integer                                                            :: bmny
     integer                                                            :: bmxx
     integer                                                            :: bmxy
     integer                                                            :: rmnx
     integer                                                            :: rmny
     integer                                                            :: rmxx
     integer                                                            :: rmxy
     integer                                                            :: ncoords
     integer                                                            :: npts_blend
     integer                                                            :: npts_relax
     integer                                                            :: nx
     integer                                                            :: ny
     integer                                                            :: nz
     integer                                                            :: scoords
  end type merge_struct           ! type merge_struct
  type meteo_struct
     real(r_kind),              dimension(:,:),             allocatable :: q
     real(r_kind),              dimension(:,:),             allocatable :: rh
     real(r_kind),              dimension(:,:),             allocatable :: t
     real(r_kind),              dimension(:,:),             allocatable :: u
     real(r_kind),              dimension(:,:),             allocatable :: v
     real(r_kind),              dimension(:,:),             allocatable :: virtt
     real(r_kind),              dimension(:,:),             allocatable :: wdir
     real(r_kind),              dimension(:,:),             allocatable :: wspd
     real(r_kind),              dimension(:,:),             allocatable :: wvmxrt
     integer                                                            :: ncoords
     integer                                                            :: nx
     integer                                                            :: ny
     integer                                                            :: nz     
  end type meteo_struct           ! type meteo_struct
  type moisture_struct
     real(r_kind),              dimension(:,:),             allocatable :: p
     real(r_kind),              dimension(:,:),             allocatable :: t     
     real(r_kind),              dimension(:,:),             allocatable :: q
     real(r_kind),              dimension(:,:),             allocatable :: rh
     real(r_kind),              dimension(:,:),             allocatable :: svp
     real(r_kind),              dimension(:,:),             allocatable :: swvmxrt
     real(r_kind),              dimension(:,:),             allocatable :: vp
     real(r_kind),              dimension(:,:),             allocatable :: wvmxrt
     integer                                                            :: ncoords
     integer                                                            :: nx
     integer                                                            :: ny
     integer                                                            :: nz 
  end type moisture_struct        ! type moisture_struct
  type netcdf_grid_struct
     real(r_kind),              dimension(:,:),             allocatable :: landmask
     real(r_kind),              dimension(:,:),             allocatable :: m_lat
     real(r_kind),              dimension(:,:),             allocatable :: m_lon
     real(r_kind),              dimension(:,:),             allocatable :: topo
     real(r_kind),              dimension(:,:),             allocatable :: u_lat
     real(r_kind),              dimension(:,:),             allocatable :: u_lon
     real(r_kind),              dimension(:,:),             allocatable :: v_lat
     real(r_kind),              dimension(:,:),             allocatable :: v_lon
     real(r_kind)                                                       :: dx
     real(r_kind)                                                       :: dy
     integer                                                            :: nx
     integer                                                            :: nxs
     integer                                                            :: ny
     integer                                                            :: nys
     integer                                                            :: nz
  end type netcdf_grid_struct     ! type netcdf_grid_struct
  type netcdf_var_struct
     character(len=20)                                                  :: name
     real(r_kind),              dimension(:,:,:,:),         allocatable :: pres
     real(r_kind),              dimension(:,:,:,:),         allocatable :: var
     real(r_kind),              dimension(:,:,:),           allocatable :: psfc
     real(r_kind)                                                       :: clat
     real(r_kind)                                                       :: clon
     integer                                                            :: filetype
     integer                                                            :: gridtype
     integer                                                            :: nt
     integer                                                            :: nx
     integer                                                            :: ny
     integer                                                            :: nz
  end type netcdf_var_struct      ! type netcdf_var_struct
  type nmm_json_grid_struct
     character(len=20)                                                  :: lmsk_name
     character(len=20)                                                  :: m_lat_name
     character(len=20)                                                  :: m_lon_name
     character(len=20)                                                  :: topo_name
     character(len=20)                                                  :: u_lat_name
     character(len=20)                                                  :: u_lon_name
     character(len=20)                                                  :: v_lat_name
     character(len=20)                                                  :: v_lon_name
     character(len=1)                                                   :: arakawa_type
     logical                                                            :: m_lat_radians
     logical                                                            :: m_lon_radians
     logical                                                            :: u_lat_radians
     logical                                                            :: u_lon_radians
     logical                                                            :: v_lat_radians
     logical                                                            :: v_lon_radians
  end type nmm_json_grid_struct   ! type nmm_json_grid_struct
  type nmm_json_var_struct
     character(len=20)                                                  :: remap_variable_name
     character(len=20)                                                  :: variable_name
     logical                                                            :: z_staggered
     logical                                                            :: replace
     integer                                                            :: gridtype
     integer                                                            :: intrptype
     integer                                                            :: levtype
  end type nmm_json_var_struct    ! type nmm_json_var_struct
  type presremap_struct
     logical                                                            :: dst_zstagger
     logical                                                            :: src_zstagger
     real(r_kind),              dimension(:,:),             allocatable :: dst_var  
     real(r_kind),              dimension(:,:),             allocatable :: dst_pres 
     real(r_kind),              dimension(:,:),             allocatable :: remap_var 
     real(r_kind),              dimension(:,:),             allocatable :: src_var     
     real(r_kind),              dimension(:,:),             allocatable :: src_pres
     integer                                                            :: dst_ncoords
     integer                                                            :: dst_nz
     integer                                                            :: src_ncoords
     integer                                                            :: src_nz
  end type presremap_struct       ! type presremap_struct
  type recenter_struct
     real(r_kind),              dimension(:,:),             allocatable :: var
     real(r_kind),              dimension(:),               allocatable :: lon
     real(r_kind),              dimension(:),               allocatable :: lat
     real(r_kind),              dimension(:),               allocatable :: radius
     real(r_kind),              dimension(:),               allocatable :: angle
     real(r_kind)                                                       :: area
     real(r_kind)                                                       :: darea
     real(r_kind)                                                       :: dradius
     real(r_kind)                                                       :: clon
     real(r_kind)                                                       :: clat
     integer                                                            :: nx
     integer                                                            :: ny
     integer                                                            :: nz
     integer                                                            :: ncoords
  end type recenter_struct        ! type recenter_struct
  type slint_struct
     real(r_double),            dimension(:,:),             allocatable :: coeffs
     real(r_kind),              dimension(:),               allocatable :: var
     real(r_kind),              dimension(:),               allocatable :: lon
     real(r_kind),              dimension(:),               allocatable :: lat
     integer,                   dimension(:,:),             allocatable :: nn
     integer,                   dimension(:),               allocatable :: idx
     integer                                                            :: ncoeffs
     integer                                                            :: ncoords
     integer                                                            :: dst_ncoords
     integer                                                            :: src_ncoords
  end type slint_struct           ! type slint_struct
  type smthr9pt_struct
     real(r_kind),              dimension(:),               allocatable :: filter
     real(r_kind),              dimension(:),               allocatable :: var
     real(r_kind)                                                       :: thrshratio
     integer                                                            :: cidx
     integer                                                            :: minx
     integer                                                            :: maxx
     integer                                                            :: miny
     integer                                                            :: maxy
     integer                                                            :: ncoords
     integer                                                            :: nx
     integer                                                            :: ny
  end type smthr9pt_struct        ! smthr9pt_struct 
  type spline_struct
     real(r_kind),              dimension(:),               allocatable :: xa
     real(r_kind),              dimension(:),               allocatable :: ya
     real(r_kind)                                                       :: x
     real(r_kind)                                                       :: y
     integer                                                            :: n
  end type spline_struct          ! type spline_struct
  type statgrid_struct
     real(r_kind)                                                       :: mean
     real(r_kind)                                                       :: vari
     real(r_kind)                                                       :: stdev
     real(r_kind)                                                       :: varmin
     real(r_kind)                                                       :: varmax
     integer                                                            :: nvals
  end type statgrid_struct        ! type statgrid_struct
  type tcv_struct
     character(len=9)                                                   :: name
     character(len=4)                                                   :: center
     character(len=3)                                                   :: id
     character(len=1)                                                   :: depth
     character(len=1)                                                   :: latns
     character(len=1)                                                   :: lonew
     real(r_kind)                                                       :: area_mnlat
     real(r_kind)                                                       :: area_mnlon
     real(r_kind)                                                       :: area_mxlat
     real(r_kind)                                                       :: area_mxlon
     real(r_kind)                                                       :: lat
     real(r_kind)                                                       :: lon
     integer                                                            :: century
     integer                                                            :: hhmm
     integer                                                            :: lati
     integer                                                            :: loni
     integer                                                            :: pcen
     integer                                                            :: penv
     integer                                                            :: penvrad
     integer                                                            :: r15ne
     integer                                                            :: r15se
     integer                                                            :: r15sw
     integer                                                            :: r15nw
     integer                                                            :: stdir
     integer                                                            :: stspd
     integer                                                            :: vmax
     integer                                                            :: vmaxrad
     integer                                                            :: yymmdd
  end type tcv_struct             ! type tcv_struct
  type timeinfo_struct
     character(len=19)                                                  :: timestamp
     real(r_double)                                                     :: jday
     integer                                                            :: dd
     integer                                                            :: hh
     integer                                                            :: mm
     integer                                                            :: nn
     integer                                                            :: ss
     integer                                                            :: yy
  end type timeinfo_struct        ! type timeinfo_struct
  type vargrid_struct
     real(r_kind),              dimension(:),               allocatable :: var
     real(r_kind)                                                       :: cutoff
     integer                                                            :: xmin
     integer                                                            :: xmax
     integer                                                            :: ymin
     integer                                                            :: ymax
     integer                                                            :: nx
     integer                                                            :: ny
     integer                                                            :: nvals
  end type vargrid_struct         ! type vargrid_struct
  type varinfo_struct
     character(len=500),        dimension(:,:,:),           allocatable :: varattrs
     character(len=25),         dimension(:),               allocatable :: varname
     character(len=10),         dimension(:),               allocatable :: dimname
     character(len=10),         dimension(:),               allocatable :: vartype
     integer,                   dimension(:,:),             allocatable :: vardimid
     integer,                   dimension(:),               allocatable :: dimid
     integer,                   dimension(:),               allocatable :: dimval
     integer,                   dimension(:),               allocatable :: varid
     integer,                   dimension(:),               allocatable :: varndims
     integer,                   dimension(:),               allocatable :: varnattrs
     integer                                                            :: nvars
     integer                                                            :: ndims
     integer                                                            :: nattrs
  end type varinfo_struct         ! type varinfo_struct
  type vertgrid_struct
     real(r_kind),              dimension(:),               allocatable :: dst_prs
     real(r_kind),              dimension(:),               allocatable :: dst_var
     real(r_kind),              dimension(:),               allocatable :: src_prs
     real(r_kind),              dimension(:),               allocatable :: src_var
     integer                                                            :: nz
  end type vertgrid_struct        ! type vertgrid_struct
  type wnd2d_struct
     real(r_kind),              dimension(:,:,:),           allocatable :: wnvar
     real(r_kind),              dimension(:),               allocatable :: var
     real(r_kind)                                                       :: maxradius
     real(r_kind)                                                       :: dangle
     real(r_kind)                                                       :: dradius
     integer                                                            :: nx
     integer                                                            :: ny
     integer                                                            :: nh
     integer                                                            :: na
     integer                                                            :: nr
  end type wnd2d_struct           ! type wnd2d_struct
  real(r_kind), parameter                                               :: clpval = tiny(0.0)
  real(r_kind), parameter                                               :: spval  = huge(0.0)
  integer,      parameter                                               :: imin   = huge(1)
  integer,      parameter                                               :: imax   = -huge(1)
  integer                                                               :: ntcs

  !-----------------------------------------------------------------------

contains

  !=======================================================================

  ! SUBROUTINE:

  ! finalize_blend_struct.f90

  ! DESCRIPTION:

  ! This subroutine deallocates memory for all arrays within the
  ! blend_struct FORTRAN structure.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN blend_struct variable.

  !-----------------------------------------------------------------------

  subroutine finalize_blend_struct(grid)

    ! Define variables passed to routine

    type(blend_struct)                                                  :: grid

    !=====================================================================

    ! Deallocate memory for local variables

    if(allocated(grid%var_blend)) deallocate(grid%var_blend)
    if(allocated(grid%var_fixed)) deallocate(grid%var_fixed)
    if(allocated(grid%var))       deallocate(grid%var)
    if(allocated(grid%lat))       deallocate(grid%lat)
    if(allocated(grid%lon))       deallocate(grid%lon)
    if(allocated(grid%mask))      deallocate(grid%mask)
    if(allocated(grid%idx))       deallocate(grid%idx)

    !=====================================================================

  end subroutine finalize_blend_struct

  !=======================================================================

  ! SUBROUTINE:

  ! finalize_fft2d_struct.f90

  ! DESCRIPTION:

  ! This subroutine deallocates memory for all arrays within the
  ! fft2d_struct FORTRAN structure.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN fft2d_struct variable.

  !-----------------------------------------------------------------------

  subroutine finalize_fft2d_struct(grid)

    ! Define variables passed to routine

    type(fft2d_struct)                                                  :: grid

    !=====================================================================

    ! Deallocate memory for local variables

    if(allocated(grid%out)) deallocate(grid%out)
    if(allocated(grid%in))  deallocate(grid%in)

    !=====================================================================

  end subroutine finalize_fft2d_struct

  !=======================================================================

  ! SUBROUTINE: 

  ! finalize_grid_struct.f90

  ! DESCRIPTION:

  ! This subroutine deallocates memory for all arrays within the
  ! grid_struct FORTRAN structure.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN grid_struct variable.

  !-----------------------------------------------------------------------

  subroutine finalize_grid_struct(grid)

    ! Define variables passed to routine

    type(grid_struct)                                                   :: grid

    !=====================================================================

    ! Deallocate memory for local variables

    if(allocated(grid%var))    deallocate(grid%var)
    if(allocated(grid%angle))  deallocate(grid%angle)
    if(allocated(grid%dx))     deallocate(grid%dx)
    if(allocated(grid%dy))     deallocate(grid%dy)
    if(allocated(grid%lon))    deallocate(grid%lon)
    if(allocated(grid%lat))    deallocate(grid%lat)
    if(allocated(grid%mapfac)) deallocate(grid%mapfac)
    if(allocated(grid%radius)) deallocate(grid%radius)
    if(allocated(grid%idx))    deallocate(grid%idx)

    !=====================================================================

  end subroutine finalize_grid_struct

  !=======================================================================

  ! SUBROUTINE:

  ! finalize_hwrf_bdy_grid_struct.f90

  ! DESCRIPTION:

  ! This subroutine deallocated memory for all arrays within the
  ! hwrf_bdy_grid_struct FORTRAN structure.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN hwrf_bdy_grid_struct variable.

  !-----------------------------------------------------------------------

  subroutine finalize_hwrf_bdy_grid_struct(grid)

    ! Define variables passed to routine

    type(hwrf_bdy_grid_struct)                                          :: grid

    !=====================================================================

    ! Deallocate memory for local variables

    if(allocated(grid%analvar)) deallocate(grid%analvar)
    if(allocated(grid%btxs))    deallocate(grid%btxs)
    if(allocated(grid%btxe))    deallocate(grid%btxe)
    if(allocated(grid%btys))    deallocate(grid%btys)
    if(allocated(grid%btye))    deallocate(grid%btye)
    if(allocated(grid%bxs))     deallocate(grid%bxs)
    if(allocated(grid%bxe))     deallocate(grid%bxe)
    if(allocated(grid%bys))     deallocate(grid%bys)
    if(allocated(grid%bye))     deallocate(grid%bye)

    !=====================================================================

  end subroutine finalize_hwrf_bdy_grid_struct

  !=======================================================================

  ! SUBROUTINE:

  ! finalize_hwrf_sfc_grid_struct.f90

  ! DESCRIPTION:

  ! This subroutine deallocated memory for all arrays within the
  ! hwrf_sfc_grid_struct FORTRAN structure.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN hwrf_sfc_grid_struct variable.

  !-----------------------------------------------------------------------

  subroutine finalize_hwrf_sfc_grid_struct(grid)

    ! Define variables passed to routine

    type(hwrf_sfc_grid_struct)                                          :: grid

    !=====================================================================

    ! Deallocate memory for local variables

    if(allocated(grid%fis))    deallocate(grid%fis)
    if(allocated(grid%sst))    deallocate(grid%sst)
    if(allocated(grid%seaice)) deallocate(grid%seaice)

    !=====================================================================

  end subroutine finalize_hwrf_sfc_grid_struct

  !=======================================================================

  ! SUBROUTINE:

  ! finalize_hwrf_vert_grid_struct.f90

  ! DESCRIPTION:

  ! This subroutine deallocated memory for all arrays within the
  ! hwrf_vert_grid_struct FORTRAN structure.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN hwrf_vert_grid_struct variable.

  !-----------------------------------------------------------------------

  subroutine finalize_hwrf_vert_grid_struct(grid)

    ! Define variables passed to routine

    type(hwrf_vert_grid_struct)                                         :: grid

    !=====================================================================

    ! Deallocate memory for local variables

    if(allocated(grid%pint))      deallocate(grid%pint)
    if(allocated(grid%nest_pint)) deallocate(grid%nest_pint)
    if(allocated(grid%plev))      deallocate(grid%plev)
    if(allocated(grid%nest_plev)) deallocate(grid%nest_plev)
    if(allocated(grid%var))       deallocate(grid%var)
    if(allocated(grid%eta1))      deallocate(grid%eta1)
    if(allocated(grid%eta2))      deallocate(grid%eta2)
    if(allocated(grid%aeta1))     deallocate(grid%aeta1)
    if(allocated(grid%aeta2))     deallocate(grid%aeta2)
    if(allocated(grid%pd))        deallocate(grid%pd)
    if(allocated(grid%psfc))      deallocate(grid%psfc)

    !=====================================================================

  end subroutine finalize_hwrf_vert_grid_struct
  
  !=======================================================================

  ! SUBROUTINE: 

  ! finalize_kdtree_struct.f90

  ! DESCRIPTION:

  ! This subroutine deallocates memory for all arrays within the
  ! kdtree_struct FORTRAN structure.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN kdtree_struct variable.

  !-----------------------------------------------------------------------

  subroutine finalize_kdtree_struct(grid)

    ! Define variables passed to routine

    type(kdtree_struct)                                                 :: grid

    !=====================================================================

    ! Deallocate memory for local variables

    if(allocated(grid%r2dist)) deallocate(grid%r2dist)
    if(allocated(grid%idx))    deallocate(grid%idx)

    !=====================================================================

  end subroutine finalize_kdtree_struct
  
  !=======================================================================

  ! SUBROUTINE:

  ! finalize_merge_struct.f90

  ! DESCRIPTION:

  ! This subroutine deallocates memory for all arrays within the
  ! merge_struct FORTRAN structure.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN merge_struct variable.

  !-----------------------------------------------------------------------

  subroutine finalize_merge_struct(grid)

    ! Define variables passed to routine

    type(merge_struct)                                                  :: grid

    !=====================================================================

    ! Deallocate memory for local variables

    if(allocated(grid%var_an)) deallocate(grid%var_an)
    if(allocated(grid%var_fg)) deallocate(grid%var_fg)
    if(allocated(grid%var))    deallocate(grid%var)
    if(allocated(grid%dx))     deallocate(grid%dx)
    if(allocated(grid%dy))     deallocate(grid%dy)
    if(allocated(grid%lat_an)) deallocate(grid%lat_an)
    if(allocated(grid%lon_an)) deallocate(grid%lon_an)
    if(allocated(grid%lat_fg)) deallocate(grid%lat_fg)
    if(allocated(grid%lon_fg)) deallocate(grid%lon_fg)
    if(allocated(grid%mask))   deallocate(grid%mask)
    if(allocated(grid%radius)) deallocate(grid%radius)
    if(allocated(grid%idx))    deallocate(grid%idx)

    !=====================================================================

  end subroutine finalize_merge_struct

  !=======================================================================

  ! SUBROUTINE:

  ! finalize_meteo_struct.f90

  ! DESCRIPTION:

  ! This subroutine deallocates memory for all arrays within the
  ! FORTRAN meteo_struct variable.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN meteo_struct variable.

  !-----------------------------------------------------------------------

  subroutine finalize_meteo_struct(grid)

    ! Define variables passed to routine

    type(meteo_struct)                                                  :: grid

    !=====================================================================

    ! Deallocate memory for local variables

    if(allocated(grid%q))      deallocate(grid%q)
    if(allocated(grid%rh))     deallocate(grid%rh)
    if(allocated(grid%t))      deallocate(grid%t)
    if(allocated(grid%u))      deallocate(grid%u)
    if(allocated(grid%v))      deallocate(grid%v)
    if(allocated(grid%virtt))  deallocate(grid%virtt)
    if(allocated(grid%wdir))   deallocate(grid%wdir)
    if(allocated(grid%wspd))   deallocate(grid%wspd)
    if(allocated(grid%wvmxrt)) deallocate(grid%wvmxrt)

    !=====================================================================

  end subroutine finalize_meteo_struct

  !=======================================================================

  ! SUBROUTINE:

  ! finalize_moisture_struct.f90

  ! DESCRIPTION:

  ! This subroutine deallocates memory for all arrays within the
  ! FORTRAN moisture_struct variable.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN moisture_struct variable.

  !-----------------------------------------------------------------------

  subroutine finalize_moisture_struct(grid)

    ! Define variables passed to routine

    type(moisture_struct)                                               :: grid

    !=====================================================================

    ! Deallocate memory for local variables

    if(allocated(grid%p))       deallocate(grid%p)
    if(allocated(grid%t))       deallocate(grid%t)
    if(allocated(grid%q))       deallocate(grid%q)
    if(allocated(grid%rh))      deallocate(grid%rh)
    if(allocated(grid%svp))     deallocate(grid%svp)
    if(allocated(grid%swvmxrt)) deallocate(grid%swvmxrt)
    if(allocated(grid%vp))      deallocate(grid%vp)
    if(allocated(grid%wvmxrt))  deallocate(grid%wvmxrt)

    !=====================================================================

  end subroutine finalize_moisture_struct

  !=======================================================================

  ! SUBROUTINE:

  ! finalize_netcdf_grid_struct.f90

  ! DESCRIPTION:

  ! This subroutine deallocated memory for all arrays within the
  ! netcdf_grid_struct FORTRAN structure.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN netcdf_grid_struct variable.

  !-----------------------------------------------------------------------

  subroutine finalize_netcdf_grid_struct(grid)

    ! Define variables passed to routine

    type(netcdf_grid_struct)                                            :: grid

    !=====================================================================

    ! Deallocate memory for local variables

    if(allocated(grid%landmask)) deallocate(grid%landmask)
    if(allocated(grid%m_lat))    deallocate(grid%m_lat)
    if(allocated(grid%m_lon))    deallocate(grid%m_lon)
    if(allocated(grid%topo))     deallocate(grid%topo)
    if(allocated(grid%u_lat))    deallocate(grid%u_lat)
    if(allocated(grid%u_lon))    deallocate(grid%u_lon)
    if(allocated(grid%v_lat))    deallocate(grid%v_lat)
    if(allocated(grid%v_lon))    deallocate(grid%v_lon)

    !=====================================================================

  end subroutine finalize_netcdf_grid_struct

  !=======================================================================

  ! SUBROUTINE:

  ! finalize_netcdf_var_struct.f90

  ! DESCRIPTION:

  ! This subroutine deallocated memory for all arrays within the
  ! netcdf_var_struct FORTRAN structure.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN netcdf_var_struct variable.

  !-----------------------------------------------------------------------

  subroutine finalize_netcdf_var_struct(grid)

    ! Define variables passed to routine

    type(netcdf_var_struct)                                             :: grid

    !=====================================================================

    ! Deallocate memory for local variables

    if(allocated(grid%pres)) deallocate(grid%pres)
    if(allocated(grid%var))  deallocate(grid%var)
    if(allocated(grid%psfc)) deallocate(grid%psfc)

    !=====================================================================

  end subroutine finalize_netcdf_var_struct

  !=======================================================================

  ! SUBROUTINE:

  ! finalize_presremap_struct.f90

  ! DESCRIPTION:

  ! This subroutine deallocated memory for all arrays within the
  ! presremap_struct FORTRAN structure.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN presremap_struct variable.

  !-----------------------------------------------------------------------

  subroutine finalize_presremap_struct(grid)

    ! Define variables passed to routine

    type(presremap_struct)                                              :: grid

    !=====================================================================

    ! Deallocate memory for local variables

    if(allocated(grid%dst_var))   deallocate(grid%dst_var)
    if(allocated(grid%dst_pres))  deallocate(grid%dst_pres)
    if(allocated(grid%remap_var)) deallocate(grid%remap_var)
    if(allocated(grid%src_var))   deallocate(grid%src_var)
    if(allocated(grid%src_pres))  deallocate(grid%src_pres)

    !=====================================================================

  end subroutine finalize_presremap_struct

  !=======================================================================

  ! SUBROUTINE: 

  ! finalize_recenter_struct.f90

  ! DESCRIPTION:

  ! This subroutine deallocates memory for all arrays within the
  ! recenter_struct FORTRAN structure.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN recenter_struct variable.

  !-----------------------------------------------------------------------

  subroutine finalize_recenter_struct(grid)

    ! Define variables passed to routine

    type(recenter_struct)                                               :: grid

    !=====================================================================

    ! Deallocate memory for local variables

    if(allocated(grid%lon))    deallocate(grid%lon)
    if(allocated(grid%lat))    deallocate(grid%lat)
    if(allocated(grid%radius)) deallocate(grid%radius)
    if(allocated(grid%angle))  deallocate(grid%angle)
    if(allocated(grid%var))    deallocate(grid%var)

    !=====================================================================

  end subroutine finalize_recenter_struct

  !=======================================================================

  ! SUBROUTINE: 

  ! finalize_slint_struct.f90

  ! DESCRIPTION:

  ! This subroutine deallocates memory for all arrays within the
  ! slint_struct FORTRAN structure.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN slint_struct variable.

  !-----------------------------------------------------------------------

  subroutine finalize_slint_struct(grid)

    ! Define variables passed to routine

    type(slint_struct)                                                  :: grid

    !=====================================================================

    ! Deallocate memory for local variables

    if(allocated(grid%coeffs)) deallocate(grid%coeffs)
    if(allocated(grid%nn))     deallocate(grid%nn)
    if(allocated(grid%var))    deallocate(grid%var)
    if(allocated(grid%lon))    deallocate(grid%lon)
    if(allocated(grid%lat))    deallocate(grid%lat)
    if(allocated(grid%idx))    deallocate(grid%idx)

    !=====================================================================

  end subroutine finalize_slint_struct

  !=======================================================================

  ! SUBROUTINE: 

  ! finalize_smthr9pt_struct.f90

  ! DESCRIPTION:

  ! This subroutine deallocates memory for all arrays within the
  ! smthr9pt_struct FORTRAN structure.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN smthr9pt_struct variable.

  !-----------------------------------------------------------------------

  subroutine finalize_smthr9pt_struct(grid)

    ! Define variables passed to routine

    type(smthr9pt_struct)                                               :: grid

    !=====================================================================

    ! Deallocate memory for local variables

    if(allocated(grid%filter)) deallocate(grid%filter)
    if(allocated(grid%var))    deallocate(grid%var)

    !=====================================================================
    
  end subroutine finalize_smthr9pt_struct
  
  !=======================================================================

  ! SUBROUTINE: 

  ! finalize_spline_struct.f90

  ! DESCRIPTION:

  ! This subroutine deallocates memory for all arrays within the
  ! spline_struct FORTRAN structure.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN spline_struct variable.

  !-----------------------------------------------------------------------

  subroutine finalize_spline_struct(grid)

    ! Define variables passed to routine

    type(spline_struct)                                                 :: grid

    !=====================================================================

    ! Deallocate memory for local variables

    if(allocated(grid%xa)) deallocate(grid%xa)
    if(allocated(grid%ya)) deallocate(grid%ya)

    !=====================================================================
    
  end subroutine finalize_spline_struct

  !=======================================================================

  ! SUBROUTINE:

  ! finalize_vargrid_struct.f90

  ! DESCRIPTION:

  ! This subroutine deallocates memory for all arrays within the
  ! vargrid_struct FORTRAN structure.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN vargrid_struct variable.

  !-----------------------------------------------------------------------

  subroutine finalize_vargrid_struct(grid)

    ! Define variables passed to routine

    type(vargrid_struct)                                                :: grid

    !=====================================================================

    ! Deallocate memory for local variables

    if(allocated(grid%var)) deallocate(grid%var)

    !=====================================================================

  end subroutine finalize_vargrid_struct

  !=======================================================================

  ! SUBROUTINE:

  ! finalize_varinfo_struct.f90

  ! DESCRIPTION:

  ! This subroutine deallocates memory for all arrays within the
  ! varinfo_struct FORTRAN structure.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN varinfo_struct variable.

  !-----------------------------------------------------------------------

  subroutine finalize_varinfo_struct(grid)

    ! Define variables passed to routine

    type(varinfo_struct)                                                :: grid

    !=====================================================================

    ! Deallocate memory for local variables

    if(allocated(grid%varattrs))  deallocate(grid%varattrs)
    if(allocated(grid%varname))   deallocate(grid%varname)
    if(allocated(grid%vartype))   deallocate(grid%vartype)
    if(allocated(grid%dimname))   deallocate(grid%dimname)
    if(allocated(grid%dimval))    deallocate(grid%dimval)
    if(allocated(grid%dimid))     deallocate(grid%dimid)
    if(allocated(grid%vardimid))  deallocate(grid%vardimid)
    if(allocated(grid%varid))     deallocate(grid%varid)
    if(allocated(grid%varndims))  deallocate(grid%varndims)
    if(allocated(grid%varnattrs)) deallocate(grid%varnattrs)

    !=====================================================================

  end subroutine finalize_varinfo_struct
  
  !=======================================================================

  ! SUBROUTINE:

  ! finalize_vertgrid_struct.f90

  ! DESCRIPTION:

  ! This subroutine deallocates memory for all arrays within the
  ! vertgrid_struct FORTRAN structure.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN vertgrid_struct variable.

  !-----------------------------------------------------------------------

  subroutine finalize_vertgrid_struct(grid)

    ! Define variables passed to routine

    type(vertgrid_struct)                                               :: grid

    !=====================================================================

    ! Deallocate memory for local variables

    if(allocated(grid%dst_prs)) deallocate(grid%dst_prs)
    if(allocated(grid%dst_var)) deallocate(grid%dst_var)
    if(allocated(grid%src_prs)) deallocate(grid%src_prs)
    if(allocated(grid%src_var)) deallocate(grid%src_var)

    !=====================================================================

  end subroutine finalize_vertgrid_struct
  
  !=======================================================================

  ! SUBROUTINE:

  ! finalize_wnd2d_struct.f90

  ! DESCRIPTION:

  ! This subroutine deallocates memory for all arrays within the
  ! wnd2d_struct FORTRAN structure.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN wnd2d_struct variable.

  !-----------------------------------------------------------------------

  subroutine finalize_wnd2d_struct(grid)

    ! Define variables passed to routine

    type(wnd2d_struct)                                                  :: grid

    !=====================================================================

    ! Deallocate memory for local variables

    if(allocated(grid%wnvar)) deallocate(grid%wnvar)
    if(allocated(grid%var))   deallocate(grid%var)

    !=====================================================================

  end subroutine finalize_wnd2d_struct

  !=======================================================================

  ! SUBROUTINE: 

  ! initialize_blend_struct.f90

  ! DESCRIPTION:

  ! This subroutine allocates memory for all arrays within the
  ! blend_struct FORTRAN structure.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN blend_struct variable.

  ! OUTPUT VARIABLES:

  ! * grid; a FORTRAN blend_struct variable where all arrays are
  !   allocated and initialized (when necessary).

  !-----------------------------------------------------------------------

  subroutine initialize_blend_struct(grid)

    ! Define variables passed to routine

    type(blend_struct)                                                  :: grid

    !=====================================================================

    ! Define local variables

    grid%ncoords = (grid%nx*grid%ny)

    ! Allocate memory for local variables

    if(.not. allocated(grid%var_blend))                                    &
         & allocate(grid%var_blend(grid%ncoords,grid%nz))
    if(.not. allocated(grid%var_fixed))                                    &
         & allocate(grid%var_fixed(grid%ncoords,grid%nz))
    if(.not. allocated(grid%var))                                          &
         & allocate(grid%var(grid%ncoords,grid%nz))
    if(.not. allocated(grid%lat))                                          &
         & allocate(grid%lat(grid%ncoords))
    if(.not. allocated(grid%lon))                                          &
         & allocate(grid%lon(grid%ncoords))
    if(.not. allocated(grid%mask))                                         &
         & allocate(grid%mask(grid%ncoords))
    if(.not. allocated(grid%idx))                                          &
         & allocate(grid%idx(grid%ncoords))

    !=====================================================================

  end subroutine initialize_blend_struct

  !=======================================================================

  ! SUBROUTINE: 

  ! initialize_fft2d_struct.f90

  ! DESCRIPTION:

  ! This subroutine allocates memory for all arrays within the
  ! fft2d_struct FORTRAN structure.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN fft2d_struct variable.

  ! OUTPUT VARIABLES:

  ! * grid; a FORTRAN fft2d_struct variable where all arrays are
  !   allocated and initialized (when necessary).

  !-----------------------------------------------------------------------

  subroutine initialize_fft2d_struct(grid)

    ! Define variables passed to routine

    type(fft2d_struct)                                                  :: grid

    !=====================================================================

    ! Allocate memory for local variables

    if(.not. allocated(grid%out)) allocate(grid%out(grid%nx,grid%ny))
    if(.not. allocated(grid%in))  allocate(grid%in(grid%nx,grid%ny))

    !=====================================================================

  end subroutine initialize_fft2d_struct

  !=======================================================================

  ! SUBROUTINE: 

  ! initialize_grid_struct.f90

  ! DESCRIPTION:

  ! This subroutine allocates memory for all arrays within the
  ! grid_struct FORTRAN structure.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN grid_struct variable.

  ! OUTPUT VARIABLES:

  ! * grid; a FORTRAN grid_struct variable where all arrays are
  !   allocated and initialized (when necessary).

  !-----------------------------------------------------------------------

  subroutine initialize_grid_struct(grid)

    ! Define variables passed to routine

    type(grid_struct)                                                   :: grid

    !=====================================================================

    ! Define local variables

    if(grid%ncoords .le. 0) grid%ncoords = (grid%nx*grid%ny)

    ! Allocate memory for local variables

    if(.not. allocated(grid%var))                                          &
         & allocate(grid%var(grid%ncoords,grid%nz))
    if(.not. allocated(grid%angle))                                        &
         & allocate(grid%angle(grid%ncoords))
    if(.not. allocated(grid%dx))                                           &
         & allocate(grid%dx(grid%ncoords))
    if(.not. allocated(grid%dy))                                           &
         & allocate(grid%dy(grid%ncoords))
    if(.not. allocated(grid%lon))                                          &
         & allocate(grid%lon(grid%ncoords))
    if(.not. allocated(grid%lat))                                          &
         & allocate(grid%lat(grid%ncoords))
    if(.not. allocated(grid%mapfac))                                       &
         & allocate(grid%mapfac(grid%ncoords))
    if(.not. allocated(grid%radius))                                       &
         & allocate(grid%radius(grid%ncoords))
    if(.not. allocated(grid%idx))                                          &
         & allocate(grid%idx(grid%ncoords))

    !=====================================================================

  end subroutine initialize_grid_struct

  !=======================================================================

  ! SUBROUTINE:

  ! initialize_hwrf_bdy_grid_struct.f90

  ! DESCRIPTION:

  ! This subroutine allocates (and initializes, when necessary) all
  ! arrays within the hwrf_bdy_grid_struct FORTRAN structure.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN hwrf_bdy_grid_struct variable.

  ! OUTPUT VARIABLES:

  ! * grid; a FORTRAN hwrf_bdy_grid_struct variable where all arrays
  !   are allocated and initialized (when necessary).

  !-----------------------------------------------------------------------

  subroutine initialize_hwrf_bdy_grid_struct(grid)

    ! Define variables passed to routine

    type(hwrf_bdy_grid_struct)                                          :: grid

    !=====================================================================

    ! Allocate memory for local variables

    if(.not. allocated(grid%analvar))                                      &
         & allocate(grid%analvar(grid%nx,grid%ny,grid%nz))
    if(.not. allocated(grid%btxs))                                         &
         & allocate(grid%btxs(grid%nx,grid%ny,grid%nz,grid%nb,grid%nt))
    if(.not. allocated(grid%btxe))                                         &
         & allocate(grid%btxe(grid%nx,grid%ny,grid%nz,grid%nb,grid%nt))
    if(.not. allocated(grid%btys))                                         &
         & allocate(grid%btys(grid%nx,grid%ny,grid%nz,grid%nb,grid%nt))
    if(.not. allocated(grid%btye))                                         &
         & allocate(grid%btye(grid%nx,grid%ny,grid%nz,grid%nb,grid%nt))
    if(.not. allocated(grid%bxs))                                          &
         & allocate(grid%bxs(grid%nx,grid%ny,grid%nz,grid%nb,grid%nt))
    if(.not. allocated(grid%bxe))                                          &
         & allocate(grid%bxe(grid%nx,grid%ny,grid%nz,grid%nb,grid%nt))
    if(.not. allocated(grid%bys))                                          &
         & allocate(grid%bys(grid%nx,grid%ny,grid%nz,grid%nb,grid%nt))
    if(.not. allocated(grid%bye))                                          &
         & allocate(grid%bye(grid%nx,grid%ny,grid%nz,grid%nb,grid%nt))

    ! Define local variables

    grid%btxs = 0.0
    grid%btxe = 0.0
    grid%btys = 0.0
    grid%btye = 0.0
    grid%bxs  = 0.0
    grid%bxe  = 0.0
    grid%bys  = 0.0
    grid%bye  = 0.0

    !=====================================================================

  end subroutine initialize_hwrf_bdy_grid_struct

  !=======================================================================

  ! SUBROUTINE:

  ! initialize_hwrf_sfc_grid_struct.f90

  ! DESCRIPTION:

  ! This subroutine allocates (and initializes, when necessary) all
  ! arrays within the hwrf_sfc_grid_struct FORTRAN structure.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN hwrf_sfc_grid_struct variable.

  ! OUTPUT VARIABLES:

  ! * grid; a FORTRAN hwrf_sfc_grid_struct variable where all arrays
  !   are allocated and initialized (when necessary).

  !-----------------------------------------------------------------------

  subroutine initialize_hwrf_sfc_grid_struct(grid)

    ! Define variables passed to routine

    type(hwrf_sfc_grid_struct)                                          :: grid

    !=====================================================================

    ! Allocate memory for local variables

    if(.not. allocated(grid%fis))    allocate(grid%fis(grid%ncoords))
    if(.not. allocated(grid%sst))    allocate(grid%sst(grid%ncoords))
    if(.not. allocated(grid%seaice)) allocate(grid%seaice(grid%ncoords))

    !=====================================================================

  end subroutine initialize_hwrf_sfc_grid_struct

  !=======================================================================

  ! SUBROUTINE:

  ! initialize_hwrf_vert_grid_struct.f90

  ! DESCRIPTION:

  ! This subroutine allocates (and initializes, when necessary) all
  ! arrays within the hwrf_vert_grid_struct FORTRAN structure.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN hwrf_vert_grid_struct variable.

  ! OUTPUT VARIABLES:

  ! * grid; a FORTRAN hwrf_vert_grid_struct variable where all arrays
  !   are allocated and initialized (when necessary).

  !-----------------------------------------------------------------------

  subroutine initialize_hwrf_vert_grid_struct(grid)

    ! Define variables passed to routine

    type(hwrf_vert_grid_struct)                                         :: grid

    !=====================================================================

    ! Allocate memory for local variables

    if(.not. allocated(grid%pint))                                         &
         & allocate(grid%pint(grid%ncoords,(grid%nz + 1)))
    if(.not. allocated(grid%nest_pint))                                    &
         & allocate(grid%nest_pint(grid%ncoords,(grid%nz + 1)))
    if(.not. allocated(grid%plev))                                         &
         & allocate(grid%plev(grid%ncoords,grid%nz))
    if(.not. allocated(grid%nest_plev))                                    &
         & allocate(grid%nest_plev(grid%ncoords,grid%nz))
    if(.not. allocated(grid%var))                                          &
         & allocate(grid%var(grid%ncoords,(grid%nz + 1)))
    if(.not. allocated(grid%eta1))                                         &
         & allocate(grid%eta1((grid%nz + 1)))
    if(.not. allocated(grid%eta2))                                         &
         & allocate(grid%eta2((grid%nz + 1)))
    if(.not. allocated(grid%aeta1))                                        &
         & allocate(grid%aeta1(grid%nz))
    if(.not. allocated(grid%aeta2))                                        &
         & allocate(grid%aeta2(grid%nz))
    if(.not. allocated(grid%pd))                                           &
         & allocate(grid%pd(grid%ncoords))
    if(.not. allocated(grid%psfc))                                         &
         & allocate(grid%psfc(grid%ncoords))

    !=====================================================================

  end subroutine initialize_hwrf_vert_grid_struct

  !=======================================================================

  ! SUBROUTINE: 

  ! initialize_kdtree_struct.f90

  ! DESCRIPTION:

  ! This subroutine allocates memory for all arrays within the
  ! kdtree_struct FORTRAN structure.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN kdtree_struct variable.

  ! OUTPUT VARIABLES:

  ! * grid; a FORTRAN kdtree_struct variable where all arrays are
  !   allocated and initialized (when necessary).

  !-----------------------------------------------------------------------

  subroutine initialize_kdtree_struct(grid)

    ! Define variables passed to routine

    type(kdtree_struct)                                                 :: grid

    !=====================================================================

    ! Allocate memory for local variables

    if(.not. allocated(grid%r2dist))                                       &
         & allocate(grid%r2dist(grid%ncoords,grid%nn))
    if(.not. allocated(grid%idx))                                          &
         & allocate(grid%idx(grid%ncoords,grid%nn))

    !=====================================================================

  end subroutine initialize_kdtree_struct  

  !=======================================================================

  ! SUBROUTINE: 

  ! initialize_merge_struct.f90

  ! DESCRIPTION:

  ! This subroutine allocates memory for all arrays within the
  ! merge_struct FORTRAN structure.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN merge_struct variable.

  ! OUTPUT VARIABLES:

  ! * grid; a FORTRAN merge_struct variable where all arrays are
  !   allocated and initialized (when necessary).

  !-----------------------------------------------------------------------

  subroutine initialize_merge_struct(grid)

    ! Define variables passed to routine

    type(merge_struct)                                                  :: grid

    !=====================================================================

    ! Define local variables

    grid%ncoords = (grid%nx*grid%ny)

    ! Allocate memory for local variables

    if(.not. allocated(grid%var_an))                                       &
         & allocate(grid%var_an(grid%ncoords,grid%nz))
    if(.not. allocated(grid%var_fg))                                       &
         & allocate(grid%var_fg(grid%ncoords,grid%nz))
    if(.not. allocated(grid%var))                                          &
         & allocate(grid%var(grid%ncoords,grid%nz))
    if(.not. allocated(grid%dx))                                           &
         & allocate(grid%dx(grid%ncoords))
    if(.not. allocated(grid%dy))                                           &
         & allocate(grid%dy(grid%ncoords))
    if(.not. allocated(grid%lat_an))                                       &
         & allocate(grid%lat_an(grid%ncoords))
    if(.not. allocated(grid%lon_an))                                       &
         & allocate(grid%lon_an(grid%ncoords))
    if(.not. allocated(grid%lat_fg))                                       &
         & allocate(grid%lat_fg(grid%ncoords))
    if(.not. allocated(grid%lon_fg))                                       &
         & allocate(grid%lon_fg(grid%ncoords))
    if(.not. allocated(grid%mask))                                         &
         & allocate(grid%mask(grid%ncoords))
    if(.not. allocated(grid%radius))                                       &
         & allocate(grid%radius(grid%ncoords))
    if(.not. allocated(grid%idx))                                          &
         & allocate(grid%idx(grid%ncoords))

    !=====================================================================

  end subroutine initialize_merge_struct

  !=======================================================================

  ! SUBROUTINE:

  ! initialize_meteo_struct.f90

  ! DESCRIPTION:

  ! This subroutine allocates memory and initializes (if applicable)
  ! all arrays within the FORTRAN meteo_struct variable.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN meteo_struct variable.
  
  ! OUTPUT VARIABLES:

  ! * grid; a FORTRAN meteo_struct variable where all arrays are
  !   allocated and initialized (when necessary).

  !-----------------------------------------------------------------------

  subroutine initialize_meteo_struct(grid)

    ! Define variables passed to routine

    type(meteo_struct)                                                  :: grid

    !=====================================================================

    ! Define local variables

    grid%ncoords = (grid%nx*grid%ny)
    
    ! Allocate memory for local variables

    if(.not. allocated(grid%q))                                            &
         & allocate(grid%q(grid%ncoords,grid%nz))
    if(.not. allocated(grid%rh))                                           &
         & allocate(grid%rh(grid%ncoords,grid%nz))
    if(.not. allocated(grid%t))                                            &
         & allocate(grid%t(grid%ncoords,grid%nz))
    if(.not. allocated(grid%u))                                            &
         & allocate(grid%u(grid%ncoords,grid%nz))
    if(.not. allocated(grid%v))                                            &
         & allocate(grid%v(grid%ncoords,grid%nz))
    if(.not. allocated(grid%virtt))                                        &
         & allocate(grid%virtt(grid%ncoords,grid%nz))
    if(.not. allocated(grid%wdir))                                         &
         & allocate(grid%wdir(grid%ncoords,grid%nz))
    if(.not. allocated(grid%wspd))                                         &
         & allocate(grid%wspd(grid%ncoords,grid%nz))
    if(.not. allocated(grid%wvmxrt))                                       &
         & allocate(grid%wvmxrt(grid%ncoords,grid%nz))

    ! Define local variables

    grid%q      = spval
    grid%rh     = spval
    grid%u      = spval
    grid%v      = spval
    grid%virtt  = spval
    grid%wdir   = spval
    grid%wspd   = spval
    grid%wvmxrt = spval

    !=====================================================================

  end subroutine initialize_meteo_struct

  !=======================================================================

  ! SUBROUTINE:

  ! initialize_moisture_struct.f90

  ! DESCRIPTION:

  ! This subroutine allocates memory and initializes (if applicable)
  ! all arrays within the FORTRAN moisture_struct variable.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN moisture_struct variable.
  
  ! OUTPUT VARIABLES:

  ! * grid; a FORTRAN moisture_struct variable where all arrays are
  !   allocated and initialized (when necessary).

  !-----------------------------------------------------------------------

  subroutine initialize_moisture_struct(grid)

    ! Define variables passed to routine

    type(moisture_struct)                                               :: grid

    !=====================================================================

    ! Define local variables

    grid%ncoords = (grid%nx*grid%ny)

    ! Allocate memory for local variables

    if(.not. allocated(grid%p))                                            &
         & allocate(grid%p(grid%ncoords,grid%nz))
    if(.not. allocated(grid%t))                                            &
         & allocate(grid%t(grid%ncoords,grid%nz))
    if(.not. allocated(grid%q))                                            &
         & allocate(grid%q(grid%ncoords,grid%nz))
    if(.not. allocated(grid%rh))                                           &
         & allocate(grid%rh(grid%ncoords,grid%nz))
    if(.not. allocated(grid%svp))                                          &
         & allocate(grid%svp(grid%ncoords,grid%nz))
    if(.not. allocated(grid%swvmxrt))                                      &
         & allocate(grid%swvmxrt(grid%ncoords,grid%nz))
    if(.not. allocated(grid%vp))                                           &
         & allocate(grid%vp(grid%ncoords,grid%nz))
    if(.not. allocated(grid%wvmxrt))                                       &
         & allocate(grid%wvmxrt(grid%ncoords,grid%nz))

    ! Define local variables

    grid%q       = spval
    grid%rh      = spval
    grid%vp      = spval
    grid%wvmxrt  = spval

    !=====================================================================

  end subroutine initialize_moisture_struct

  !=======================================================================

  ! SUBROUTINE:

  ! initialize_netcdf_grid_struct.f90

  ! DESCRIPTION:

  ! This subroutine allocates (and initializes, when necessary) all
  ! arrays within the netcdf_grid_struct FORTRAN structure.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN netcdf_grid_struct variable.

  ! OUTPUT VARIABLES:

  ! * grid; a FORTRAN netcdf_grid_struct variable where all arrays are
  !   allocated and initialized (when necessary).

  !-----------------------------------------------------------------------

  subroutine initialize_netcdf_grid_struct(grid)

    ! Define variables passed to routine

    type(netcdf_grid_struct)                                            :: grid

    !=====================================================================

    ! Check local variable and proceed accordingly

    if(is_hwrf) then

       ! Allocate memory for local variables

       if(.not. allocated(grid%landmask))                                  &
            & allocate(grid%landmask(grid%nx,grid%ny))
       if(.not. allocated(grid%m_lat))                                     &
            & allocate(grid%m_lat(grid%nx,grid%ny))
       if(.not. allocated(grid%m_lon))                                     &
            & allocate(grid%m_lon(grid%nx,grid%ny))
       if(.not. allocated(grid%topo))                                      &
            & allocate(grid%topo(grid%nx,grid%ny))
       if(.not. allocated(grid%u_lat))                                     &
            & allocate(grid%u_lat(grid%nxs,grid%ny))
       if(.not. allocated(grid%u_lon))                                     &
            & allocate(grid%u_lon(grid%nxs,grid%ny))
       if(.not. allocated(grid%v_lat))                                     &
            & allocate(grid%v_lat(grid%nx,grid%nys))
       if(.not. allocated(grid%v_lon))                                     &
            & allocate(grid%v_lon(grid%nx,grid%nys))

    end if ! if(is_hwrf)

    !=====================================================================

  end subroutine initialize_netcdf_grid_struct

  !=======================================================================

  ! SUBROUTINE:

  ! initialize_netcdf_var_struct.f90

  ! DESCRIPTION:

  ! This subroutine allocates (and initializes, when necessary) all
  ! arrays within the netcdf_var_struct FORTRAN structure.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN netcdf_var_struct variable.

  ! OUTPUT VARIABLES:

  ! * grid; a FORTRAN netcdf_var_struct variable where all arrays are
  !   allocated and initialized (when necessary).

  !-----------------------------------------------------------------------

  subroutine initialize_netcdf_var_struct(grid)

    ! Define variables passed to routine

    type(netcdf_var_struct)                                             :: grid

    !=====================================================================

    ! Allocate memory for local variables

    if(.not. allocated(grid%pres))                                         &
         & allocate(grid%pres(grid%nx,grid%ny,grid%nz,grid%nt))
    if(.not. allocated(grid%var))                                          &
         & allocate(grid%var(grid%nx,grid%ny,grid%nz,grid%nt))
    if(.not. allocated(grid%psfc))                                         &
         & allocate(grid%psfc(grid%nx,grid%ny,grid%nt))

    !=====================================================================

  end subroutine initialize_netcdf_var_struct

  !=======================================================================

  ! SUBROUTINE: 

  ! initialize_presremap_struct.f90

  ! DESCRIPTION:

  ! This subroutine allocates memory for all arrays within the
  ! presremap_struct FORTRAN structure.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN presremap_struct variable.

  ! OUTPUT VARIABLES:

  ! * grid; a FORTRAN presremap_struct variable where all arrays are
  !   allocated and initialized (when necessary).

  !-----------------------------------------------------------------------

  subroutine initialize_presremap_struct(grid)

    ! Define variables passed to routine

    type(presremap_struct)                                              :: grid

    !=====================================================================

    ! Allocate memory for local variables

    if(.not. allocated(grid%dst_var))                                      &
         & allocate(grid%dst_var(grid%dst_ncoords,grid%dst_nz))
    if(.not. allocated(grid%dst_pres))                                     &
         & allocate(grid%dst_pres(grid%dst_ncoords,grid%dst_nz))
    if(.not. allocated(grid%remap_var))                                    &
         & allocate(grid%remap_var(grid%src_ncoords,grid%dst_nz))
    if(.not. allocated(grid%src_var))                                      &
         & allocate(grid%src_var(grid%src_ncoords,grid%src_nz))
    if(.not. allocated(grid%src_pres))                                     &
         & allocate(grid%src_pres(grid%src_ncoords,grid%src_nz))

    !=====================================================================
    
  end subroutine initialize_presremap_struct

  !=======================================================================

  ! SUBROUTINE: 

  ! initialize_recenter_struct.f90

  ! DESCRIPTION:

  ! This subroutine allocates memory for all arrays within the
  ! recenter_struct FORTRAN structure.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN recenter_struct variable.

  ! OUTPUT VARIABLES:

  ! * grid; a FORTRAN recenter_struct variable where all arrays are
  !   allocated and initialized (when necessary).

  !-----------------------------------------------------------------------

  subroutine initialize_recenter_struct(grid)

    ! Define variables passed to routine

    type(recenter_struct)                                               :: grid

    !=====================================================================

    ! Define local variables

    grid%ncoords = (grid%nx*grid%ny)

    ! Allocate memory for local variables

    if(.not. allocated(grid%var))                                          &
         & allocate(grid%var(grid%ncoords,grid%nz))
    if(.not. allocated(grid%lon))                                          &
         & allocate(grid%lon(grid%ncoords))
    if(.not. allocated(grid%lat))                                          &
         & allocate(grid%lat(grid%ncoords))
    if(.not. allocated(grid%radius))                                       &
         & allocate(grid%radius(grid%ncoords))
    if(.not. allocated(grid%angle))                                        &
         & allocate(grid%angle(grid%ncoords))

    !=====================================================================
    
  end subroutine initialize_recenter_struct

  !=======================================================================

  ! SUBROUTINE: 

  ! initialize_slint_struct.f90

  ! DESCRIPTION:

  ! This subroutine allocates memory for all arrays within the
  ! slint_struct FORTRAN structure.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN slint_struct variable.

  ! OUTPUT VARIABLES:

  ! * grid; a FORTRAN slint_struct variable where all arrays are
  !   allocated and initialized (when necessary).

  !-----------------------------------------------------------------------

  subroutine initialize_slint_struct(grid)

    ! Define variables passed to routine

    type(slint_struct)                                                  :: grid

    !=====================================================================

    ! Define local variables

    grid%ncoeffs = 3

    ! Allocate memory for local variables

    if(.not. allocated(grid%coeffs))                                       &
         & allocate(grid%coeffs(grid%ncoeffs,grid%ncoords))
    if(.not. allocated(grid%nn))                                           &
         & allocate(grid%nn(grid%ncoeffs,grid%ncoords))
    if(.not. allocated(grid%var))                                          &
         & allocate(grid%var(grid%ncoords))
    if(.not. allocated(grid%lon))                                          &
         & allocate(grid%lon(grid%ncoords))
    if(.not. allocated(grid%lat))                                          &
         & allocate(grid%lat(grid%ncoords))
    if(.not. allocated(grid%idx))                                          &
         & allocate(grid%idx(grid%ncoords))

    !=====================================================================
    
  end subroutine initialize_slint_struct

  !=======================================================================

  ! SUBROUTINE: 
  
  ! initialize_smthr9pt_struct.f90

  ! DESCRIPTION:

  ! This subroutine allocates memory for all arrays within the
  ! smthr9pt_struct FORTRAN structure.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN smthr9pt_struct variable.

  ! OUTPUT VARIABLES:

  ! * grid; a FORTRAN smthr9pt_struct variable where all arrays are
  !   allocated and initialized (when necessary).

  !-----------------------------------------------------------------------

  subroutine initialize_smthr9pt_struct(grid)

    ! Define variables passed to routine

    type(smthr9pt_struct)                                               :: grid

    !=====================================================================

    ! Check local variable and proceed accordingly

    if(grid%ncoords .eq. 0) grid%ncoords = (grid%nx*grid%ny)

    ! Allocate memory for local variables

    if(.not. allocated(grid%filter)) allocate(grid%filter(grid%ncoords))
    if(.not. allocated(grid%var))    allocate(grid%var(grid%ncoords))

    !=====================================================================
    
  end subroutine initialize_smthr9pt_struct  

  !=======================================================================

  ! SUBROUTINE: 
  
  ! initialize_spline_struct.f90

  ! DESCRIPTION:

  ! This subroutine allocates memory for all arrays within the
  ! spline_struct FORTRAN structure.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN spline_struct variable.

  ! OUTPUT VARIABLES:

  ! * grid; a FORTRAN spline_struct variable where all arrays are
  !   allocated and initialized (when necessary).

  !-----------------------------------------------------------------------

  subroutine initialize_spline_struct(grid)

    ! Define variables passed to routine

    type(spline_struct)                                                 :: grid

    !=====================================================================

    ! Allocate memory for local variables

    if(.not. allocated(grid%xa)) allocate(grid%xa(grid%n))
    if(.not. allocated(grid%ya)) allocate(grid%ya(grid%n))

    !=====================================================================
    
  end subroutine initialize_spline_struct

  !=======================================================================

  ! SUBROUTINE: 

  ! initialize_vargrid_struct.f90

  ! DESCRIPTION:

  ! This subroutine allocates memory for all arrays within the
  ! vargrid_struct FORTRAN structure.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN vargrid_struct variable.

  ! OUTPUT VARIABLES:

  ! * grid; a FORTRAN vargrid_struct variable where all arrays are
  !   allocated and initialized (when necessary).

  !-----------------------------------------------------------------------

  subroutine initialize_vargrid_struct(grid)

    ! Define variables passed to routine

    type(vargrid_struct)                                                :: grid

    !=====================================================================

    ! Check local variable and proceed accordingly

    if(grid%nvals .eq. 0) then

       ! Define local variables

       grid%nvals = (grid%nx*grid%ny)

    end if ! if(grid%nvals .eq. 0)

    ! Allocate memory for local variables

    if(.not. allocated(grid%var)) allocate(grid%var(grid%nvals))
    
    !=====================================================================

  end subroutine initialize_vargrid_struct

  !=======================================================================

  ! SUBROUTINE: 

  ! initialize_varinfo_struct.f90

  ! DESCRIPTION:

  ! This subroutine allocates memory for all arrays within the
  ! varinfo_struct FORTRAN structure.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN varinfo_struct variable.

  ! OUTPUT VARIABLES:

  ! * grid; a FORTRAN varinfo_struct variable where all arrays are
  !   allocated and initialized (when necessary).

  !-----------------------------------------------------------------------

  subroutine initialize_varinfo_struct(grid)

    ! Define variables passed to routine

    type(varinfo_struct)                                                :: grid

    !=====================================================================

    ! Allocate memory for local variables

    if(.not. allocated(grid%varattrs))                                     &
         & allocate(grid%varattrs(grid%nvars,grid%nattrs,2))
    if(.not. allocated(grid%vardimid))                                     &
         & allocate(grid%vardimid(grid%nvars,grid%ndims))
    if(.not. allocated(grid%varname))                                      &
         & allocate(grid%varname(grid%nvars))
    if(.not. allocated(grid%vartype))                                      &
         & allocate(grid%vartype(grid%nvars))
    if(.not. allocated(grid%varndims))                                     &
         & allocate(grid%varndims(grid%nvars))
    if(.not. allocated(grid%varnattrs))                                    &
         & allocate(grid%varnattrs(grid%nvars))
    if(.not. allocated(grid%varid))                                        &
         & allocate(grid%varid(grid%nvars))
    if(.not. allocated(grid%dimval))                                       &
         & allocate(grid%dimval(grid%ndims))
    if(.not. allocated(grid%dimname))                                      &
         & allocate(grid%dimname(grid%ndims))
    if(.not. allocated(grid%dimid))                                        &
         & allocate(grid%dimid(grid%ndims))

    !=====================================================================

  end subroutine initialize_varinfo_struct  

  !=======================================================================

  ! SUBROUTINE: 

  ! initialize_vertgrid_struct.f90

  ! DESCRIPTION:

  ! This subroutine allocates memory for all arrays within the
  ! vertgrid_struct FORTRAN structure.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN vertgrid_struct variable.

  ! OUTPUT VARIABLES:

  ! * grid; a FORTRAN vertgrid_struct variable where all arrays are
  !   allocated and initialized (when necessary).

  !-----------------------------------------------------------------------

  subroutine initialize_vertgrid_struct(grid)

    ! Define variables passed to routine

    type(vertgrid_struct)                                               :: grid

    !=====================================================================

    ! Allocate memory for local variables

    if(.not. allocated(grid%dst_prs)) allocate(grid%dst_prs(grid%nz))
    if(.not. allocated(grid%dst_var)) allocate(grid%dst_var(grid%nz))    
    if(.not. allocated(grid%src_prs)) allocate(grid%src_prs(grid%nz))
    if(.not. allocated(grid%src_var)) allocate(grid%src_var(grid%nz))

    !=====================================================================

  end subroutine initialize_vertgrid_struct
    
  !=======================================================================

  ! SUBROUTINE: 

  ! initialize_wnd2d_struct.f90

  ! DESCRIPTION:

  ! This subroutine allocates memory for all arrays within the
  ! wnd2d_struct FORTRAN structure.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN wnd2d_struct variable.

  ! OUTPUT VARIABLES:

  ! * grid; a FORTRAN wnd2d_struct variable where all arrays are
  !   allocated and initialized (when necessary).

  !-----------------------------------------------------------------------

  subroutine initialize_wnd2d_struct(grid)

    ! Define variables passed to routine

    type(wnd2d_struct)                                                  :: grid

    !=====================================================================

    ! Define local variables

    grid%nr = grid%maxradius/grid%dradius
    grid%na = 360.0/grid%dangle
    grid%nh = ((grid%na/2) + 1)

    ! Allocate memory for local variables

    if(.not. allocated(grid%wnvar))                                        &
         & allocate(grid%wnvar(grid%nh,grid%nr,grid%na))
    if(.not. allocated(grid%var))                                          &
         & allocate(grid%var(grid%nx*grid%ny))

    !=====================================================================

  end subroutine initialize_wnd2d_struct

  !=======================================================================

end module variable_interface
