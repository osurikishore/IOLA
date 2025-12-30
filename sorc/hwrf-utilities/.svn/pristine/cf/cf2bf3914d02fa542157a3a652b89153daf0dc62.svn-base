module fileio_interface

  !=======================================================================

  !$$$ MODULE DOCUMENTATION BLOCK
  
  ! obs_to_bufr :: fileio_interface
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

  use bufrio_interface
  use kinds_interface
  use variable_interface

  ! Define interfaces and attributes for module routines

  implicit none
  private
  public :: fileio_interface_read
  public :: fileio_interface_write
  public :: nfiles
  public :: obs_filename
  interface fileio_interface_read
     module procedure read_convobserr
     module procedure read_hsa
     module procedure read_obs_filenames
  end interface fileio_interface_read
  interface fileio_interface_write
  end interface fileio_interface_write

  ! Define local variables

  character(len=500),           dimension(:),               allocatable :: obs_filename
  integer                                                               :: nfiles

  !-----------------------------------------------------------------------

contains

  !=======================================================================

  ! SUBROUTINE:

  ! read_convobserr.f90

  ! DESCRIPTION:

  ! This subroutine reads the user specified external file containing
  ! the observation-type, to be written to the BUFR file, observation
  ! error profiles (see the readme files accompanying this
  ! repository).
  
  ! INPUT VARIABLES:

  ! * filename; a FORTRAN character string specifying the full-path to
  !   the external file containing the observation-type observation
  !   error profiles.

  ! * grid; a FORTRAN error_struct variable.

  ! OUTPUT VARIABLES:

  ! * grid; a FORTRAN error_struct variable containing the
  !   observation-type observation error profiles.

  !-----------------------------------------------------------------------

  subroutine read_convobserr(filename,grid)

    ! Define variables passed to routine

    type(error_struct)                                                  :: grid
    character(len=500)                                                  :: filename

    ! Define variables computed within routine

    character(len=500)                                                  :: dummy

    ! Define counting variables

    integer                                                             :: i

    !=====================================================================

    ! Define local variables

    grid%nz = 0
    open(99,file=trim(adjustl(filename)),form='formatted')
1000 read(99,*,end=1001) dummy
    grid%nz = grid%nz + 1
    goto 1000
1001 continue
    close(99)
    call variable_interface_setup_struct(grid)
    open(99,file=trim(adjustl(filename)),form='formatted')

    ! Loop through local variable

    do i = 1, grid%nz

       ! Define local variables

       read(99,*) grid%plev(i), grid%p(i), grid%t(i), grid%q(i),           &
            & grid%z(i), grid%uv(i)

    end do ! do i = 1, grid%nz

    ! Define local variables

    close(99)

    !=====================================================================

  end subroutine read_convobserr

  !=======================================================================

  ! SUBROUTINE:

  ! read_hsa.f90

  ! DESCRIPTION:

  ! This subroutine reads a National Oceanic and Atmospheric
  ! Administration (NOAA) Atlantic Oceanographic and Meteorological
  ! Laboratory (AOML) Hurricane Research Division (HRD) HRD Spline
  ! Analysis (HSA) formatted observation files.

  ! INPUT VARIABLES:

  ! * filename; a FORTRAN character string specifying the full-path to
  !   the HSA formatted file.

  ! * grid; a FORTRAN hsa_struct variable.

  ! OUTPUT VARIABLES:

  ! * grid; a FORTRAN hsa_struct variable containing the HSA file
  !   contents.

  !-----------------------------------------------------------------------

  subroutine read_hsa(filename,grid)

    ! Define variables passed to routine

    type(hsa_struct)                                                    :: grid
    character(len=500)                                                  :: filename

    ! Define variables computed within routine

    real(r_kind)                                                        :: dummy

    ! Define counting variables

    integer                                                             :: i

    !=====================================================================

    ! Define local variables

    grid%nz = 0
    open(99,file=trim(adjustl(filename)),form='formatted')
1000 read(99,*,end=1001) dummy
    grid%nz = grid%nz + 1
    goto 1000
1001 continue
    close(99)
    call variable_interface_setup_struct(grid)
    open(99,file=trim(adjustl(filename)),form='formatted')

    ! Loop through local variable

    do i = 1, grid%nz

       ! Define local variables

       read(99,500,err=1002) grid%wx(i), grid%yymmdd(i), grid%gmt(i),     &
            & grid%lat(i), grid%lon(i), grid%p(i), grid%t(i), grid%rh(i), &
            & grid%z(i), grid%u(i), grid%v(i), grid%tail(i)

    end do ! do i = 1, grid%nz

    ! Define local variables

    close(99)
    return
1002 continue

    ! Deallocate memory for local variables

    call variable_interface_cleanup_struct(grid)

    ! Define local variables

    call variable_interface_setup_struct(grid)
500 format(i2,1x,f7.0,1x,i4,1x,2(f7.3,1x),3(f6.1,1x),f7.1,2(f6.1,1x),a4)

    !=====================================================================

  end subroutine read_hsa

  !=======================================================================

  ! SUBROUTINE:

  ! read_obs_filenames.f90

  ! DESCRIPTION:

  ! This subroutine defines the global variables 'nfiles' and
  ! 'obs_filename' array.

  ! INPUT VARIABLES:

  ! * filename; a FORTRAN character string containing the list of
  !   observation filename pathes.

  !-----------------------------------------------------------------------

  subroutine read_obs_filenames(filename)

    ! Define variables passed to routine

    character(len=500)                                                  :: filename

    ! Define variables computed within routine

    character(len=500)                                                  :: dummy
    
    ! Define counting variables

    integer                                                             :: i

    !=====================================================================

    ! Define local variables

    nfiles = 0
    open(99,file=trim(adjustl(filename)),form='formatted')
1000 read(99,*,end=1001) dummy
    nfiles = nfiles + 1
    goto 1000
1001 continue
    close(99)

    ! Allocate memory for local variables

    if(.not. allocated(obs_filename)) allocate(obs_filename(nfiles))

    ! Define local variables

    open(99,file=trim(adjustl(filename)),form='formatted')

    ! Loop through local variable

    do i = 1, nfiles

       ! Define local variables

       read(99,*) obs_filename(i)

    end do ! do i = 1, nfiles

    ! Define local variables

    close(99)

    !=====================================================================

  end subroutine read_obs_filenames

  !=======================================================================

end module fileio_interface
