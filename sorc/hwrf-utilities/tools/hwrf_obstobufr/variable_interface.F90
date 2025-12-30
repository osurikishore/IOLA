module variable_interface

  !=======================================================================

  !$$$ MODULE DOCUMENTATION BLOCK
  
  ! obs_to_bufr :: variable_interface
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

  ! Define interfaces and attributes for module routines
  
  implicit none
  private
  public :: bufr_struct
  public :: error_struct
  public :: hsa_spval
  public :: hsa_struct
  public :: interp_spline_struct
  public :: meteo_struct
  public :: spval
  public :: variable_interface_cleanup_struct
  public :: variable_interface_setup_struct
  interface variable_interface_cleanup_struct
     module procedure finalize_bufr_struct
     module procedure finalize_error_struct
     module procedure finalize_hsa_struct
     module procedure finalize_interp_spline_struct
  end interface variable_interface_cleanup_struct
  interface variable_interface_setup_struct
     module procedure initialize_bufr_struct
     module procedure initialize_error_struct
     module procedure initialize_hsa_struct
     module procedure initialize_interp_spline_struct
  end interface variable_interface_setup_struct

  ! Define local variables

  type bufr_struct
     character(len=80)                                                  :: obstr
     character(len=80)                                                  :: hdstr
     character(len=80)                                                  :: qcstr
     character(len=80)                                                  :: oestr
     character(len=19)                                                  :: cdate
     character(len=8)                                                   :: subset
     real(r_double),            dimension(:,:),             allocatable :: obs
     real(r_double),            dimension(:,:),             allocatable :: qcf
     real(r_double),            dimension(:,:),             allocatable :: oer
     real(r_double),            dimension(:),               allocatable :: hdr
     integer                                                            :: mxmn
     integer                                                            :: mxlv
     integer                                                            :: idate
  end type bufr_struct            ! type bufr_struct
  type error_struct
     real(r_kind),              dimension(:),               allocatable :: plev
     real(r_kind),              dimension(:),               allocatable :: p
     real(r_kind),              dimension(:),               allocatable :: t
     real(r_kind),              dimension(:),               allocatable :: q
     real(r_kind),              dimension(:),               allocatable :: z
     real(r_kind),              dimension(:),               allocatable :: uv
     integer                                                            :: nz     
  end type error_struct           ! type error_struct
  type hsa_struct
     character(len=4),          dimension(:),               allocatable :: tail
     real(r_kind),              dimension(:),               allocatable :: lat
     real(r_kind),              dimension(:),               allocatable :: lon
     real(r_kind),              dimension(:),               allocatable :: p
     real(r_kind),              dimension(:),               allocatable :: rh
     real(r_kind),              dimension(:),               allocatable :: t
     real(r_kind),              dimension(:),               allocatable :: u
     real(r_kind),              dimension(:),               allocatable :: v
     real(r_kind),              dimension(:),               allocatable :: z
     real(r_kind),              dimension(:),               allocatable :: yymmdd
     integer,                   dimension(:),               allocatable :: gmt
     integer,                   dimension(:),               allocatable :: wx     
     integer                                                            :: nz
  end type hsa_struct             ! type hsa_struct
  type interp_spline_struct
     real(r_kind),              dimension(:),               allocatable :: xa
     real(r_kind),              dimension(:),               allocatable :: ya
     real(r_kind)                                                       :: x
     real(r_kind)                                                       :: y
     integer                                                            :: n    
  end type interp_spline_struct   ! type interp_spline_struct
  type meteo_struct
     real(r_double)                                                     :: p
     real(r_double)                                                     :: q
     real(r_double)                                                     :: rh
     real(r_double)                                                     :: t
     real(r_double)                                                     :: wvmxrt
  end type meteo_struct           ! type meteo_struct
  real(r_double), parameter                                             :: bufr_spval = 10.e10
  real(r_double), parameter                                             :: spval      = huge(1.0)
  real(r_kind),   parameter                                             :: hsa_spval  = -99.0
  integer,        parameter                                             :: bufr_mxlv  = 200
  integer,        parameter                                             :: bufr_mxmn  = 35

  !-----------------------------------------------------------------------

contains

  !=======================================================================

  ! SUBROUTINE:

  ! finalize_bufr_struct.f90

  ! DESCRIPTION:

  ! This subroutine deallocates memory for all arrays within the
  ! bufr_struct FORTRAN structure.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN bufr_struct variable.

  !-----------------------------------------------------------------------

  subroutine finalize_bufr_struct(grid)

    ! Define variables passed routine

    type(bufr_struct)                                                   :: grid

    !=====================================================================

    ! Deallocate memory for local variables

    if(allocated(grid%obs)) deallocate(grid%obs)
    if(allocated(grid%qcf)) deallocate(grid%qcf)
    if(allocated(grid%oer)) deallocate(grid%oer)
    if(allocated(grid%hdr)) deallocate(grid%hdr)

    !=====================================================================

  end subroutine finalize_bufr_struct

  !=======================================================================

  ! SUBROUTINE:

  ! finalize_error_struct.f90

  ! DESCRIPTION:

  ! This subroutine deallocates memory for all arrays within the
  ! error_struct FORTRAN structure.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN error_struct variable.

  !-----------------------------------------------------------------------

  subroutine finalize_error_struct(grid)

    ! Define variables passed routine

    type(error_struct)                                                  :: grid

    !=====================================================================

    ! Deallocate memory for local variables

    if(allocated(grid%plev)) deallocate(grid%plev)
    if(allocated(grid%p))    deallocate(grid%p)
    if(allocated(grid%t))    deallocate(grid%t)
    if(allocated(grid%q))    deallocate(grid%q)
    if(allocated(grid%z))    deallocate(grid%z)
    if(allocated(grid%uv))   deallocate(grid%uv)

    !=====================================================================

  end subroutine finalize_error_struct

  !=======================================================================

  ! SUBROUTINE:

  ! finalize_hsa_struct.f90

  ! DESCRIPTION:

  ! This subroutine deallocates memory for all arrays within the
  ! hsa_struct FORTRAN structure.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN hsa_struct variable.

  !-----------------------------------------------------------------------

  subroutine finalize_hsa_struct(grid)

    ! Define variables passed routine

    type(hsa_struct)                                                    :: grid

    !=====================================================================

    ! Deallocate memory for local variables

    if(allocated(grid%tail))   deallocate(grid%tail)
    if(allocated(grid%lat))    deallocate(grid%lat)
    if(allocated(grid%lon))    deallocate(grid%lon)
    if(allocated(grid%p))      deallocate(grid%p)
    if(allocated(grid%rh))     deallocate(grid%rh)
    if(allocated(grid%t))      deallocate(grid%t)
    if(allocated(grid%u))      deallocate(grid%u)
    if(allocated(grid%v))      deallocate(grid%v)
    if(allocated(grid%z))      deallocate(grid%z)
    if(allocated(grid%yymmdd)) deallocate(grid%yymmdd)
    if(allocated(grid%gmt))    deallocate(grid%gmt)
    if(allocated(grid%wx))     deallocate(grid%wx)

    !=====================================================================

  end subroutine finalize_hsa_struct

  !=======================================================================

  ! SUBROUTINE: 

  ! finalize_interp_spline_struct.f90

  ! DESCRIPTION:

  ! This subroutine deallocates memory for all arrays within the
  ! interp_spline_struct FORTRAN structure.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN interp_spline_struct variable.

  !-----------------------------------------------------------------------

  subroutine finalize_interp_spline_struct(grid)

    ! Define variables passed to routine

    type(interp_spline_struct)                                          :: grid

    !=====================================================================

    ! Deallocate memory for local variables

    if(allocated(grid%xa)) deallocate(grid%xa)
    if(allocated(grid%ya)) deallocate(grid%ya)

    !=====================================================================
    
  end subroutine finalize_interp_spline_struct

  !=======================================================================

  ! SUBROUTINE:

  ! initialize_bufr_struct.f90

  ! DESCRIPTION:

  ! This subroutine allocates memory for all arrays within the
  ! bufr_struct FORTRAN structure.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN bufr_struct variable containing the variables
  !   necessary to allocate and initialize the respective variable
  !   arrays.

  ! OUTPUT VARIABLES:

  ! * grid; a FORTRAN bufr_struct variable containing allocated and
  !   initialized variable arrays.

  !-----------------------------------------------------------------------

  subroutine initialize_bufr_struct(grid)

    ! Define variables passed routine

    type(bufr_struct)                                                   :: grid

    !=====================================================================

    ! Check local variable and proceed accordingly

    if(grid%mxmn .eq. 0) grid%mxmn = bufr_mxmn       
    if(grid%mxlv .eq. 0) grid%mxlv = bufr_mxlv

    ! Allocate memory for local variables

    if(.not. allocated(grid%obs)) allocate(grid%obs(grid%mxmn,grid%mxlv))
    if(.not. allocated(grid%qcf)) allocate(grid%qcf(grid%mxmn,grid%mxlv))
    if(.not. allocated(grid%oer)) allocate(grid%oer(grid%mxmn,grid%mxlv))
    if(.not. allocated(grid%hdr)) allocate(grid%hdr(grid%mxmn))

    ! Define local variables

    grid%obs = bufr_spval
    grid%qcf = bufr_spval
    grid%oer = bufr_spval
    grid%hdr = bufr_spval

    !=====================================================================

  end subroutine initialize_bufr_struct

  !=======================================================================

  ! SUBROUTINE:

  ! initialize_error_struct.f90

  ! DESCRIPTION:

  ! This subroutine allocates memory for all arrays within the
  ! error_struct FORTRAN structure.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN error_struct variable containing the variables
  !   necessary to allocate and initialize the respective variable
  !   arrays.

  ! OUTPUT VARIABLES:

  ! * grid; a FORTRAN error_struct variable containing allocated and
  !   initialized variable arrays.

  !-----------------------------------------------------------------------

  subroutine initialize_error_struct(grid)

    ! Define variables passed routine

    type(error_struct)                                                  :: grid

    !=====================================================================

    ! Allocate memory for local variables

    if(.not. allocated(grid%plev)) allocate(grid%plev(grid%nz))
    if(.not. allocated(grid%p))    allocate(grid%p(grid%nz))
    if(.not. allocated(grid%t))    allocate(grid%t(grid%nz))
    if(.not. allocated(grid%q))    allocate(grid%q(grid%nz))
    if(.not. allocated(grid%z))    allocate(grid%z(grid%nz))
    if(.not. allocated(grid%uv))   allocate(grid%uv(grid%nz))

    !=====================================================================

  end subroutine initialize_error_struct

  !=======================================================================

  ! SUBROUTINE:

  ! initialize_hsa_struct.f90

  ! DESCRIPTION:

  ! This subroutine allocates memory for all arrays within the
  ! hsa_struct FORTRAN structure.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN hsa_struct variable containing the variables
  !   necessary to allocate and initialize the respective variable
  !   arrays.

  ! OUTPUT VARIABLES:

  ! * grid; a FORTRAN hsa_struct variable containing allocated and
  !   initialized variable arrays.

  !-----------------------------------------------------------------------

  subroutine initialize_hsa_struct(grid)

    ! Define variables passed routine

    type(hsa_struct)                                                    :: grid

    !=====================================================================

    ! Allocate memory for local variables

    if(.not. allocated(grid%tail))   allocate(grid%tail(grid%nz))
    if(.not. allocated(grid%lat))    allocate(grid%lat(grid%nz))
    if(.not. allocated(grid%lon))    allocate(grid%lon(grid%nz))
    if(.not. allocated(grid%p))      allocate(grid%p(grid%nz))
    if(.not. allocated(grid%rh))     allocate(grid%rh(grid%nz))
    if(.not. allocated(grid%t))      allocate(grid%t(grid%nz))
    if(.not. allocated(grid%u))      allocate(grid%u(grid%nz))
    if(.not. allocated(grid%v))      allocate(grid%v(grid%nz))
    if(.not. allocated(grid%z))      allocate(grid%z(grid%nz))
    if(.not. allocated(grid%yymmdd)) allocate(grid%yymmdd(grid%nz))
    if(.not. allocated(grid%gmt))    allocate(grid%gmt(grid%nz))
    if(.not. allocated(grid%wx))     allocate(grid%wx(grid%nz))

    ! Define local variables

    grid%lat = hsa_spval
    grid%lon = hsa_spval
    grid%p   = hsa_spval
    grid%rh  = hsa_spval
    grid%t   = hsa_spval
    grid%u   = hsa_spval
    grid%v   = hsa_spval
    grid%z   = hsa_spval

    !=====================================================================

  end subroutine initialize_hsa_struct

  !=======================================================================

  ! SUBROUTINE: 
  
  ! initialize_interp_spline_struct.f90

  ! DESCRIPTION:

  ! This subroutine allocates memory for all arrays within the
  ! interp_spline_struct FORTRAN structure.

  ! INPUT VARIABLES:

  ! * grid; a FORTRAN interp_spline_struct variable.

  ! OUTPUT VARIABLES:

  ! * grid; a FORTRAN interp_spline_struct variable containing
  !   allocated variable arrays.

  !-----------------------------------------------------------------------

  subroutine initialize_interp_spline_struct(grid)

    ! Define variables passed to routine

    type(interp_spline_struct)                                          :: grid

    !=====================================================================

    ! Allocate memory for local variables

    if(.not. allocated(grid%xa)) allocate(grid%xa(grid%n))
    if(.not. allocated(grid%ya)) allocate(grid%ya(grid%n))

    !=====================================================================
    
  end subroutine initialize_interp_spline_struct

  !=======================================================================

end module variable_interface
