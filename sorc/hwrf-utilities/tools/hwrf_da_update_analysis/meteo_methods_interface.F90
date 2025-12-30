module meteo_methods_interface

  !=======================================================================

  !$$$ MODULE DOCUMENTATION BLOCK
  
  ! da-update-analysis :: meteo_methods_interface
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
  use variable_interface

  ! Define interfaces and attributes for module routines
  
  implicit none
  private
  public :: meteo_methods_moisture
  public :: meteo_methods_virttemp
  public :: meteo_methods_winds

  !-----------------------------------------------------------------------

contains

  !=======================================================================

  ! SUBROUTINE:

  ! ccesvp.f90

  ! DESCRIPTION:

  ! This subroutine computes the saturation vapor pressure and
  ! saturation water vapor mixing ratio in accordance with the
  ! Clausius-Clapeyron equation.

  ! INPUT VARIABLES:

  ! * moisture; a FORTRAN moisture_struct variable containing (at
  !   minimum) the temperature (t; Kelvin) profile.

  ! OUTPUT VARIABLES:

  ! * moisture; a FORTRAN moisture_struct variable containing the
  !   saturation vapor pressure (svp; Pascals) and the saturation
  !   water vapor mixing ratio (swvmxrt; kilograms per kilogram)
  !   computed from the Clausius-Clapeyron equation.

  !-----------------------------------------------------------------------

  subroutine ccesvp(moisture)

    ! Define variables passed to routine

    type(moisture_struct)                                               :: moisture

    ! Define variables computed within routine

    real(r_kind)                                                        :: e0
    real(r_kind)                                                        :: t0

    ! Define counting variables

    integer                                                             :: i

    !=====================================================================

    ! Define local variables

    e0 = 611.0
    t0 = 273.15

    ! Loop through local variable

    do i = 1, moisture%nz

       ! Compute local variables

       moisture%svp(:,i)     = e0*exp((ltnthtvpr/rv)*((1.0/t0) -           &
            & (1.0/moisture%t(:,i))))
       moisture%swvmxrt(:,i) = (moisture%svp(:,i)*rd)/(rv*(moisture%p(:,i) &
               & - moisture%svp(:,i)))

    end do ! do j = 1, moisture%nz

    !=====================================================================

  end subroutine ccesvp

  !=======================================================================

  ! SUBROUTINE:

  ! meteo_methods_moisture.f90 

  ! DESCRIPTION:

  ! This subroutine computes all moisture quantities in accordance
  ! with the information provided by the user.

  ! INPUT VARIABLES:

  ! * moisture; a FORTRAN moisture_struct variable containing (at
  !   minimum) the total pressure (p; Pascals) and temperature (t;
  !   Kelvin) profiles; the moisture quantities are calculated based
  !   on the moisture variable(s) provided by the user for either
  !   specific humidity (q; kilograms per kilogram), relative humidity
  !   (rh; ratio), and/or water vapor mixing ratio (wvmxrt; kilograms
  !   per kilogram).

  ! OUTPUT VARIABLES:

  ! * moisture; a FORTRAN moisture_struct variable containing all
  !   derived moisture quantities:

  !   + relative humidity (rh; ratio)
  !   + saturation vapor pressure (svp; Pascals)
  !   + specific humidity (q; kilograms per kilogram)
  !   + vapor pressure (vp; Pascals)
  !   + water vapor mixing ratio (wvmxrt; kilograms per kilogram)

  !-----------------------------------------------------------------------

  subroutine meteo_methods_moisture(moisture)

    ! Define variables passed to routine

    type(moisture_struct)                                               :: moisture

    !=====================================================================

    ! Compute local variables

    call ccesvp(moisture)

    ! Check local variable and proceed accordingly

    if(moisture%rh(1,1) .ne. spval) then

       ! Compute local variables

       call rh2wvmxrt(moisture)
       call wvmxrt2q(moisture)
       goto 1000

    end if ! if(moisture%rh(1,1) .ne. spval)

    ! Check local variable and proceed accordingly

    if(moisture%q(1,1) .ne. spval) then

       ! Compute local variables

       call q2wvmxrt(moisture)
       call wvmxrt2rh(moisture)
       goto 1000

    end if ! if(moisture%q(1,1) .ne. spval)

    ! Check local variable and proceed accordingly

    if(moisture%wvmxrt(1,1) .eq. spval) then

       ! Compute local variables

       call wvmxrt2q(moisture)
       call wvmxrt2rh(moisture)
       goto 1000

    end if ! if(moisture%wvmxrt(1,1) .eq. spval)

    ! Define local variables

1000 continue

    !=====================================================================

  end subroutine meteo_methods_moisture

  !=======================================================================

  ! SUBROUTINE:

  ! meteo_methods_virttemp.f90

  ! DESCRIPTION:

  ! This subroutine will compute the virtual temperature; units are
  ! Kelvin.

  ! INPUT VARIABLES:

  ! * meteo; a FORTRAN meteo_struct variable containing (at minimum)
  !   the temperature (t) and water-vapor mixing ratio (wvmxrt); if
  !   the water-vapor mixing ratio array (e.g., wvmxrt) is not
  !   composed of missing data values (i.e., spval), the water-vapor
  !   mixing ratio is not recomputed and the water-vapor mixing ratio
  !   passed within the meteo_struct variable is used.

  ! OUTPUT VARIABLES:
  
  ! * meteo; a FORTRAN meteo_struct variable containing the virtual
  !   temperature computed from the temperature and water-vapor mixing
  !   ratio values.

  !-----------------------------------------------------------------------

  subroutine meteo_methods_virttemp(meteo)

    ! Define variables passed to routine

    type(meteo_struct)                                                  :: meteo

    ! Define counting variables

    integer                                                             :: i

    !=====================================================================

    ! Check local variable and proceed accordingly

    if(meteo%wvmxrt(1,1) .eq. spval) then

       ! Compute local variables

       call meteo_methods_wvmxrt(meteo)

    end if ! if(meteo%wvmxrt(1,1) .eq. spval)

    ! Loop through local variable

    do i = 1, meteo%nz

       ! Compute local variables

       meteo%virtt(:,i) = meteo%t(:,i)*((meteo%wvmxrt(:,i) + (rd/rv))/     &
            & ((rd/rv)*(1.0 + meteo%wvmxrt(:,i))))

    end do ! do i = 1, meteo%nz

    !=====================================================================

  end subroutine meteo_methods_virttemp

  !=======================================================================

  ! SUBROUTINE:

  ! meteo_methods_winds.f90

  ! DESCRIPTION:

  ! This subroutine computes the attributes of the wind field in
  ! accordance with the user input variable (meteo_struct); if the
  ! zonal- (u-) and meridional (v-) wind components are equal to
  ! 'spval' upon entry, they are derived from the wind speed (spd) and
  ! direction (dir); if the wind speed and direction are equal to
  ! 'spval' upon entry, they are computed from the wind components (u-
  ! and v-).

  ! INPUT VARIABLES:

  ! * meteo; a FORTRAN meteo_struct variable containing the defined
  !   wind field attributes.

  ! OUTPUT VARIABLES:

  ! * meteo; a FORTRAN meteo_struct variable containing the computed
  !   wind field attributes

  !-----------------------------------------------------------------------

  subroutine meteo_methods_winds(meteo)

    ! Define variables passed to routine

    type(meteo_struct)                                                  :: meteo

    ! Define variables computed within routine

    type(winds_struct)                                                  :: winds

    ! Define counting variables

    integer                                                             :: i
  
    !=====================================================================

    ! Define local variables

    winds%nx = meteo%nx
    winds%ny = meteo%ny
    call variable_interface_setup_struct(winds)
    
    ! Loop through local variable

    do i = 1, meteo%nz
    
       ! Check local variable and proceed accordingly

       if((minval(meteo%u(:,i)) .eq. spval) .and. (maxval(meteo%u(:,i))    &
            & .eq. spval) .and. (minval(meteo%v(:,i)) .eq. spval) .and.    &
            & (maxval(meteo%v(:,i)) .eq. spval)) then

          ! Define local variables

          winds%dir = meteo%wdir(:,i)
          winds%spd = meteo%wspd(:,i)
          
          ! Compute local variables

          call wnd_comps(winds)

          ! Define local variables

          meteo%u(:,i) = winds%u
          meteo%v(:,i) = -1.0*winds%v

       end if ! if((minval(meteo%u(:,i)) .eq. spval)
              ! .and. (maxval(meteo%u(:,i)) .eq. spval)
              ! .and. (minval(meteo%v(:,i)) .eq. spval)
              ! .and. (maxval(meteo%v(:,i)) .eq. spval))

       ! Check local variable and proceed accordingly

       if((minval(meteo%wspd(:,i)) .eq. spval) .and.                       &
            & (maxval(meteo%wspd(:,i)) .eq. spval) .and.                   &
            & (minval(meteo%wdir(:,i)) .eq. spval) .and.                   &
            & (maxval(meteo%wdir(:,i)) .eq. spval)) then

          ! Define local variables

          winds%u = meteo%u(:,i)
          winds%v = meteo%v(:,i)
          
          ! Compute local variables

          call wnd_spddir(winds)

          ! Define local variables

          meteo%wdir(:,i) = winds%dir
          meteo%wspd(:,i) = winds%spd

       end if ! if((minval(meteo%wspd(:,i)) .eq. spval)
              ! .and. (maxval(meteo%wspd(:,i)) .eq. spval)
              ! .and. (minval(meteo%wdir(:,i)) .eq. spval)
              ! .and. (maxval(meteo%wdir(:,i)) .eq. spval))
          
    end do ! do i = 1, meteo%nz

    ! Deallocate memory for local variables

    call variable_interface_cleanup_struct(winds)

    !=====================================================================
    
  end subroutine meteo_methods_winds
    
  !=======================================================================

  ! SUBROUTINE:

  ! meteo_methods_wvmxrt.f90

  ! DESCRIPTION:

  ! This subroutine computes the water-vapor mixing ratio values
  ! corresponding to either the user provided specific (kilogram per
  ! kilogram) or relative humidity (ratio) values; units are kilograms
  ! per kilogram.

  ! NOTE: Currently only the specific humidity conversion is
  !       supported.

  ! INPUT VARIABLES:

  ! * meteo; a FORTRAN meteo_struct variable containing (at minimum)
  !   either the specific humidity (q) or relative humidity (ratio)
  !   values.

  ! OUTPUT VARIABLES:

  ! * meteo; a FORTRAN meteo_struct variable containing the
  !   water-vapor mixing ration values computed from the user
  !   specified moisture information.

  !-----------------------------------------------------------------------

  subroutine meteo_methods_wvmxrt(meteo)

    ! Define variables passed to routine

    type(meteo_struct)                                                  :: meteo

    ! Define counting variables

    integer                                                             :: i

    !=====================================================================

    ! Loop through local variable

    do i = 1, meteo%nz

       ! Check local variable and proceed accordingly

       if(meteo%q(1,i) .ne. spval) then
          
          ! Compute local variables
          
!          call q2wvmxrt(meteo%q(:,i),meteo%wvmxrt(:,i))

       end if ! if(meteo%q(1,j) .ne. spval)

    end do ! do i = 1, meteo%nz

    !=====================================================================

  end subroutine meteo_methods_wvmxrt

  !=======================================================================

  ! SUBROUTINE:

  ! q2wvmxrt.f90

  ! DESCRIPTION:

  ! This subroutine computes the specific humidity profile from the
  ! water-vapor mixing ratio values profile.

  ! INPUT VARIABLES:

  ! * moisture; a FORTRAN moisture_struct variable containing (at
  !   minimum) the specific humidity (q; kilograms per kilogram)
  !   profile.

  ! OUTPUT VARIABLES:

  ! * moisture; a FORTRAN moisture_struct variable containing the
  !   water-vapor mixing ratio (wvmxrt; kilograms per kilogram)
  !   profile.

  !-----------------------------------------------------------------------

  subroutine q2wvmxrt(moisture)

    ! Define variables passed to routine

    type(moisture_struct)                                               :: moisture

    ! Define variables computed within routine

    real(r_kind),               dimension(:,:),             allocatable :: mpi_wvmxrt
    integer                                                             :: mpi_id

    ! Define counting variables

    integer                                                             :: i, j

    !=====================================================================

    ! Loop through local variable
    
    do i = 1, moisture%nz

       ! Compute local variables

       moisture%wvmxrt(:,i) = moisture%q(:,i)/(1.0 - moisture%q(:,i))

       ! Check local variable and proceed accordingly

       where(moisture%wvmxrt(:,i) .lt. 0.0) moisture%wvmxrt(:,i) = clpval

    end do ! do j = 1, moisture%nz

    !=====================================================================

  end subroutine q2wvmxrt
  
  !=======================================================================

  ! SUBROUTINE:

  ! rh2wvmxrt.f90

  ! DESCRIPTION:

  ! This subroutine will compute the water-vapor mixing ratio profile
  ! from the relative humidity and saturation water-vapor mixing
  ! ratio profiles.

  ! INPUT VARIABLES:

  ! * moisture; a FORTRAN moisture_struct variable containing (at
  !   minimum) the relative humidity (rh; ratio) and saturation
  !   water-vapor mixing ratio (swvmxrt; kilograms per kilogram)
  !   profiles.

  ! OUTPUT VARIABLES:

  ! * moisture; a FORTRAN moisture_struct variable containing the
  !   water-vapor mixing ratio (wvmxrt; kilograms per kilogram)
  !   profile.

  !-----------------------------------------------------------------------

  subroutine rh2wvmxrt(moisture)

    ! Define variables passed to routine

    type(moisture_struct)                                               :: moisture

    ! Define counting variables

    integer                                                             :: i

    !=====================================================================

    ! Loop through local variable
    
    do i = 1, moisture%nz
       
       ! Compute local variables

       moisture%wvmxrt(:,i) = moisture%swvmxrt(:,i)*moisture%rh(:,i)

       ! Check local variable and proceed accordingly

       where(moisture%wvmxrt(:,i) .lt. 0.0) moisture%wvmxrt(:,i) = clpval

    end do ! do i = 1, moisture%nz

    !=====================================================================

  end subroutine rh2wvmxrt

  !=======================================================================

  ! SUBROUTINE:

  ! wnd_comps.f90

  ! DESCRIPTION:

  ! This subroutine computes the zonal- and meridional wind (u- and v,
  ! respectively) components from the wind speed (spd) and direction
  ! (dir).

  ! INPUT VARIABLES:

  ! * winds; a FORTRAN winds_struct variable containing (at minimum)
  !   the wind speed (spd) and direction (dir).

  ! OUTPUT VARIABLES:

  ! * winds; a FORTRAN winds_struct variable containing the computed
  !   zonal- (u-) and meridional (v-) wind components.

  !-----------------------------------------------------------------------

  subroutine wnd_comps(winds)

    ! Define variables passed to routine

    type(winds_struct)                                                  :: winds

    !=====================================================================

    ! Compute local variables

    winds%u = winds%spd*cos(winds%dir*deg2rad)
    winds%v = winds%spd*sin(winds%dir*deg2rad)

    !=====================================================================

  end subroutine wnd_comps

  !=======================================================================

  ! SUBROUTINE:

  ! wnd_spddir.f90

  ! DESCRIPTION:

  ! This subroutine computes the wind speed and direction (spd and
  ! dir, respectively) from the zonal- and meridional wind (u- and v,
  ! respectively) components.

  ! INPUT VARIABLES:

  ! * winds; a FORTRAN winds_struct variable containing (at minimum)
  !   the zonal (u) and meridional (v) components.

  ! OUTPUT VARIABLES:

  ! * winds; a FORTRAN winds_struct variable containing the computed
  !   wind speed (spd) and direction (dir).

  !-----------------------------------------------------------------------

  subroutine wnd_spddir(winds)

    ! Define variables passed to routine

    type(winds_struct)                                                  :: winds

    ! Define counting variables

    integer                                                             :: i

    !=====================================================================

    ! Loop through local variable

    do i = 1, winds%ncoords

       ! Compute local variables

       winds%spd(i) = sqrt(winds%u(i)*winds%u(i) + winds%v(i)*winds%v(i))

       ! Check local variable and proceed accordingly

       if(winds%spd(i) .lt. 1.e-10) then

          ! Define local variables

          winds%dir(i) = 0.0

       else   ! if(winds%spd(i) .lt. 1.e-10)

          ! Compute local variables

          winds%dir(i) = atan2(winds%u(i),winds%v(i))*rad2deg + 270.0

       end if ! if(winds%spd(i) .lt. 1.e-10)
                 
    end do ! do i = 1, winds%ncoords
       
    !=====================================================================

  end subroutine wnd_spddir

  !=======================================================================

  ! SUBROUTINE:

  ! wvmxrt2q.f90

  ! DESCRIPTION:

  ! This subroutine will compute the specific humidity from the
  ! water-vapor mixing ratio.

  ! INPUT VARIABLES:

  ! * moisture; a FORTRAN moisture_struct variable containing (at
  !   minimum) the water-vapor mixing ratio (wvmxrt; kilograms per
  !   kilogram) profile.

  ! OUTPUT VARIABLES:

  ! * moisture; a FORTRAN moisture_struct variable containing the
  !   specfic humidity (q; kilograms per kilogram) profile.

  !-----------------------------------------------------------------------

  subroutine wvmxrt2q(moisture)

    ! Define variables passed to routine

    type(moisture_struct)                                               :: moisture

    ! Define counting variables

    integer                                                             :: i

    !=====================================================================

    ! Loop through local variable

    do i = 1, moisture%nz

       ! Compute local variables

       moisture%q(:,i) = moisture%wvmxrt(:,i)/(1.0 + moisture%wvmxrt(:,i))

       ! Check local variable and proceed accordingly

       where(moisture%q(:,i) .lt. 0.0) moisture%q(:,i) = clpval

    end do ! do i = 1, moisture%nz

    !=====================================================================
    
  end subroutine wvmxrt2q

  !=======================================================================

  ! SUBROUTINE:

  ! wvmxrt2rh.f90

  ! DESCRIPTION:

  ! This subroutine will compute the relative humidity profile from
  ! the water-vapor mixing ratio profile and the vapor pressure
  ! profile from the relative humidity and saturation vapor pressure
  ! profiles.

  ! INPUT VARIABLES:

  ! * moisture; a FORTRAN moisture_struct variable containing (at
  !   minimum) the water-vapor mixing ratio (wvmxrt; kilograms per
  !   kilogram) profile.

  ! OUTPUT VARIABLES:

  ! * moisture; a FORTRAN moisture_struct variable containing the
  !   relative humidity (rh; ratio) and vapor pressure (vp; Pascals)
  !   profiles.

  !-----------------------------------------------------------------------

  subroutine wvmxrt2rh(moisture)

    ! Define variables passed to routine

    type(moisture_struct)                                               :: moisture

    ! Define counting variables

    integer                                                             :: i

    !=====================================================================

    ! Loop through local variable
    
    do i = 1, moisture%nz

       ! Compute local variables

       moisture%rh(:,i) = moisture%wvmxrt(:,i)/moisture%swvmxrt(:,i)

       ! Check local variable and proceed accordingly

       where(moisture%rh(:,i) .lt. 0.0) moisture%rh(:,i) = clpval
     
    end do ! do j = 1, moisture%nz

    !=====================================================================
    
  end subroutine wvmxrt2rh

  !=======================================================================

end module meteo_methods_interface
