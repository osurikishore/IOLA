module da_update_analysis_interface

  !=======================================================================

  !$$$ MODULE DOCUMENTATION BLOCK
  
  ! da-update-analysis :: da_update_analysis_interface
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

  use diagnostics_interface
  use hwrf_analysis_interface
  use kinds_interface
  use namelist_interface

  ! Define interfaces and attributes for module routines
  
  implicit none
  private
  public :: da_update_analysis

  !-----------------------------------------------------------------------

contains

  !=======================================================================

  ! SUBROUTINE:

  ! da_update_analysis.f90

  ! DESCRIPTION:

  ! This subroutine is the driver-level routine for all
  ! da_update_analysis routines.

  !-----------------------------------------------------------------------

  subroutine da_update_analysis()

    ! Define variables computed within routine

    real(r_kind)                                                        :: time_start

    !=====================================================================

    ! Define local variables

    call diagnostics_time_start(time_start)
    call namelist()

    ! Compute local variables

    if(is_hwrf) call hwrf_analysis()
    
    ! Define local variables

    call diagnostics_time_stop(time_start)

    !=====================================================================

  end subroutine da_update_analysis

  !=======================================================================

end module da_update_analysis_interface
