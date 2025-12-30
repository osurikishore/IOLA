module obs_to_bufr_interface

  !=======================================================================

  !$$$ PROGRAM DOCUMENTATION BLOCK
  
  ! obs-to-bufr :: obs_to_bufr_interface
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

  ! Review the README, within the top-level directory, which provides
  ! relevant instructions and (any) references cited by algorithms
  ! within this software suite.

  !=======================================================================

  ! Define associated modules and subroutines

  use bufrio_interface
  use diagnostics_interface
  use fileio_interface
  use kinds_interface
  use math_methods_interface
  use meteo_methods_interface
  use namelist_interface
  use time_methods_interface
  use variable_interface

  ! Define interfaces and attributes for module routines
  
  implicit none
  private
  public :: obs_to_bufr

  ! Define local variables

  real(r_double)                                                        :: anljday
  real(r_double)                                                        :: obsjday
  integer                                                               :: yyyy
  integer                                                               :: mm
  integer                                                               :: dd
  integer                                                               :: hh
  integer                                                               :: nn
  integer                                                               :: ss

  !-----------------------------------------------------------------------

contains

  !=======================================================================

  ! SUBROUTINE:

  ! obs_to_bufr.f90

  ! DESCRIPTION:

  ! This subroutine is the driver routine to write observation-types
  ! to the BUFR format.

  !-----------------------------------------------------------------------

  subroutine obs_to_bufr()

    ! Define variables computed within routine

    type(error_struct)                                                  :: error
    real(r_kind)                                                        :: time_start

    !=====================================================================

    ! Define local variables

    call diagnostics_time_start(time_start)
    call namelist()
    call time_methods_date_attributes(analdate,yyyy,mm,dd,hh,nn,ss)

    ! Compute local variables

    call time_methods_julian_day(yyyy,mm,dd,hh,nn,ss,anljday)

    ! Define local variables

    call fileio_interface_read(obs_filepath)
    call fileio_interface_read(obserr_filepath,error)
    if(hsa) call hsa_to_bufr(error)

    ! Deallocate memory for local variables

    call variable_interface_cleanup_struct(error)

    ! Define local variables

    call diagnostics_time_stop(time_start)

    !=====================================================================

  end subroutine obs_to_bufr

  !=======================================================================

  ! SUBROUTINE:

  ! hsa_to_bufr.f90

  ! DESCRIPTION:

  ! This subroutine is the interface to the National Oceanic and
  ! Atmospheric Administration (NOAA) Atlantic Oceanographic and
  ! Meteorological Laboratory (AOML) Hurricane Research Division (HRD)
  ! HRD Spline Analysis (HSA) formatted files; the valid variable
  ! values, for 'SIGL' and 'MANL' type levels are written as
  ! individual records to the respective BUFR file such that the drift
  ! of the sonde (if applicable) can be appropriately accounted for.

  ! INPUT VARIABLES:

  ! * error; a FORTRAN error_struct variable containing the
  !  observation-type error profiles.

  !-----------------------------------------------------------------------

  subroutine hsa_to_bufr(error)

    ! Define variables passed to routine

    type(bufr_struct)                                                   :: bufr
    type(error_struct)                                                  :: error

    ! Define variables computed within routine

    type(hsa_struct)                                                    :: grid
    type(meteo_struct)                                                  :: meteo
    character(len=500)                                                  :: filename
    character(len=8)                                                    :: cacobid
    character(len=8)                                                    :: cdate
    character(len=6)                                                    :: ctime
    character(len=6)                                                    :: acid
    character(len=6)                                                    :: obsn
    character(len=2)                                                    :: msnid
    real(r_double)                                                      :: racobid
    real(r_kind)                                                        :: obs_err
    integer                                                             :: strstrt
    integer                                                             :: strstop

    ! Define counting variables

    integer                                                             :: i, j

    !=====================================================================

    ! Define local variables

    equivalence(racobid,cacobid)
    bufr%hdstr  = 'SID XOB YOB DHR TYP'
    bufr%obstr  = 'POB QOB TOB UOB VOB'
    bufr%qcstr  = 'PQM QQM TQM WQM'
    bufr%oestr  = 'POE QOE TOE WOE'
    bufr%subset = 'ADPUPA'
    bufr%mxmn   = 5
    bufr%mxlv   = 1
    call bufrio_interface_idate(analdate,bufr)
    call bufrio_interface_open(bufr_filepath,bufr_tblpath,bufr,.false.,    &
         & .true.)

    ! Loop through local variable

    do i = 1, nfiles

       ! Define local variables

       call fileio_interface_read(obs_filename(i),grid)
       strstrt  = index(trim(adjustl(obs_filename(i))),'/',back=.true.)    &
            & + 1
       strstop  = len(trim(adjustl(obs_filename(i))))
       filename = trim(adjustl(obs_filename(i)(strstrt:strstop)))
       msnid    = filename(1:2)
       acid     = filename(1:6) 
       obsn     = filename(7:8)

       ! Check local variable and proceed accordingly

       if(trim(adjustl(msnid)) .eq. 'AF') then

          ! Define local variables

          write(cacobid,500) trim(adjustl(acid(3:5))), trim(adjustl(obsn))

       end if ! if(trim(adjustl(msnid)) .eq. 'AF')

       ! Check local variable and proceed accordingly

       if(trim(adjustl(msnid)) .eq. 'NA') then

          ! Define local variables

          write(cacobid,501) trim(adjustl(obsn))

       end if ! if(trim(adjustl(msnid)) .eq. 'NA')

       ! Check local variable and proceed accordingly

       if(trim(adjustl(msnid)) .eq. 'NO') then

          ! Define local variables

          write(cacobid,502) trim(adjustl(acid(5:5))), trim(adjustl(obsn))

       end if ! if(trim(adjustl(msnid)) .eq. 'NO')

       ! Loop through local variable

       do j = 1, grid%nz

          ! Check local variable and proceed accordingly

          if(grid%yymmdd(j) .ne. hsa_spval .and. (grid%tail(j) .eq.        &
               & 'SIGL' .or. grid%tail(j) .eq. 'MANL')) then

             ! Define local variables

             call variable_interface_setup_struct(bufr)
             write(cdate,'(i8)')        (20000000 + (int(grid%yymmdd(j))))
             write(ctime,'(i4.4,i2.2)') grid%gmt(j), 0
             write(bufr%cdate,503) cdate(1:4), cdate(5:6), cdate(7:8),     &
                  & ctime(1:2), ctime(3:4), ctime(5:6)
             call time_methods_date_attributes(bufr%cdate,yyyy,mm,dd,hh,   &
                  & nn,ss)

             ! Compute local variables

             call time_methods_julian_day(yyyy,mm,dd,hh,nn,ss,obsjday)

             ! Define local variables

             bufr%hdr(1) = dble(racobid)
             bufr%hdr(2) = dble(-1.0*grid%lon(j) + 360.0)
             bufr%hdr(3) = dble(grid%lat(j))
             bufr%hdr(4) = (obsjday - anljday)*dble(24.0)
             bufr%hdr(5) = 137

             ! Check local variable and proceed accordingly

             if(grid%t(j) .ne. hsa_spval) then

                ! Define local variables

                bufr%obs(1,1) = grid%p(j)
                bufr%obs(3,1) = grid%t(j)
                call obs_error(error%nz,error%plev,error%t,grid%p(j),      &
                     & intrp_obserr,obs_err)
                bufr%oer(3,1) = obs_err
                bufr%qcf(1,1) = 2.0
                bufr%qcf(3,1) = 2.0

                ! Check local variable and proceed accordingly

                if(grid%rh(j) .ne. hsa_spval) then

                   ! Define local variables

                   meteo%p  = dble(grid%p(j)*100.0)
                   meteo%rh = dble(grid%rh(j))
                   meteo%t  = dble(grid%t(j) + 273.15)

                   ! Compute local variables

                   call meteo_methods_spechumd(meteo)

                   ! Define local variables
 
                   bufr%obs(2,1) = real(meteo%q)*1000.0*1000.0
                   call obs_error(error%nz,error%plev,error%q,grid%p(j),   &
                        & intrp_obserr,obs_err)
                   bufr%oer(2,1) = obs_err
                   bufr%qcf(2,1) = 2.0 

                end if ! if(grid%rh(j) .ne. hsa_spval)

                ! Define local variables

                call bufrio_interface_write(bufr)

             end if ! if(grid%t(j) .ne. hsa_spval)

	     ! Deallocate memory for local variables

	     call variable_interface_cleanup_struct(bufr)

	     ! Define local variables

	     call variable_interface_setup_struct(bufr)
             bufr%hdr(1) = dble(racobid)
             bufr%hdr(2) = dble(-1.0*grid%lon(j) + 360.0)
             bufr%hdr(3) = dble(grid%lat(j))
             bufr%hdr(4) = (obsjday - anljday)*dble(24.0)
             bufr%hdr(5) = 237

	     ! Check local variable and proceed accordingly

             if(grid%u(j) .ne. hsa_spval .and. grid%v(j) .ne. hsa_spval)   &
	          & then

                ! Define local variables

                bufr%obs(1,1) = grid%p(j)
                bufr%obs(4,1) = grid%u(j)
                bufr%obs(5,1) = grid%v(j)
                call obs_error(error%nz,error%plev,error%uv,grid%p(j),     &
                     & intrp_obserr,obs_err)
                bufr%oer(4,1) = obs_err
                bufr%qcf(1,1) = 2.0
                bufr%qcf(4,1) = 2.0
                call bufrio_interface_write(bufr)

             end if ! if(grid%u(j) .ne. hsa_spval .and. grid%v(j) .ne.     
                    ! hsa_spval)

             ! Deallocate memory for local variables
          
	     call variable_interface_cleanup_struct(bufr)

	  end if ! if(grid%yymmdd(j) .ne. hsa_spval .and. (grid%tail(j)    &
                 ! .eq. 'SIGL' .or. grid%tail(j) .eq. 'MANL')) 

       end do ! do j = 1, grid%nz

       ! Deallocate memory for local variables

       call variable_interface_cleanup_struct(grid)

    end do ! do i = 1, nfiles

    ! Define local variables

    call bufrio_interface_close(.false.,.true.)
500 format(a3,a2,'A')
501 format('872',a2,'A')
502 format('AA',a1,a2,'A')
503 format(a4,'-',a2,'-',a2,'_',a2,':',a2,':',a2)

    !=====================================================================

  end subroutine hsa_to_bufr

  !=======================================================================

  ! SUBROUTINE:

  ! obs_error.f90

  ! DESCRIPTION:

  ! This subroutine determines the observation error for a given
  ! observation provided the observation isobaric level.

  ! INPUT VARIABLES:

  ! * nz; a FORTRAN integer specifying the total number of levels in
  !   the observation error profile.

  ! * err_plev; a FORTRAN 4-byte real value array, of dimension nz,
  !   containing the isobaric levels for the observation error
  !   profile.

  ! * err_var; a FORTRAN 4-byte real value array, of dimension nz,
  !   containing the respective observation-type variable observation
  !   error profile.

  ! * obs_plev; a FORTRAN 4-byte real value specifying the observation
  !   isobaric level.

  ! * interp; a FORTRAN boolean value specifying whether to use spline
  !   interpolation or nearest-neighbor methods to define the
  !   observation error.

  ! * obs_err; a FORTRAN 4-byte real value to contain the observation
  !   error.

  ! OUTPUT VARIABLES:

  ! * obs_err; a FORTRAN 4-byte real value specifying the observation
  !   variable in accordance with the user instructions (e.g.,
  !   interp).

  !-----------------------------------------------------------------------

  subroutine obs_error(nz,err_plev,err_var,obs_plev,interp,obs_err)

    ! Define variables passed to subroutine

    logical                                                             :: interp
    integer                                                             :: nz
    real(r_kind),               dimension(nz)                           :: err_plev
    real(r_kind),               dimension(nz)                           :: err_var
    real(r_kind)                                                        :: obs_plev
    real(r_kind)                                                        :: obs_err

    ! Define variables computed within routine

    type(interp_spline_struct)                                          :: interp_spline

    ! Define counting variables

    integer                                                             :: i

    !=====================================================================

    ! Define local variables
    
    interp_spline%n  = nz
    call variable_interface_setup_struct(interp_spline)
    interp_spline%xa = err_plev
    interp_spline%ya = err_var

    ! Check local variable and proceed accordingly

    if(interp) then

       ! Define local variables

       interp_spline%x = obs_plev

       ! Compute local variables

       call math_methods_spline_interp(interp_spline)

       ! Define local variables

       obs_err = interp_spline%y

    end if ! if(interp)

    ! Check local variable and proceed accordingly

    if(.not. interp) then

       ! Define local variables
       
       call math_methods_sort_array(interp_spline,.true.,.false.)

       ! Loop through local variable

       do i = 1, (interp_spline%n - 1)

          ! Check local variable and proceed accordingly

          if(obs_plev .gt. interp_spline%xa(i) .and. obs_plev .gt.        &
               & obs_plev .le. interp_spline%xa(i+1)) then

             ! Define local variables

             obs_err = interp_spline%ya(i+1)

          end if ! if(obs_plev .gt. interp_spline%xa(i) .and. obs_plev
                 ! .gt. obs_plev .le. interp_spline%xa(i+1)

       end do ! do i = 1, (interp_spline%n - 1)

    end if ! if(.not. interp)
    
    ! Deallocate memory for local variables
    
    call variable_interface_cleanup_struct(interp_spline)

    !=====================================================================

  end subroutine obs_error

  !=======================================================================

end module obs_to_bufr_interface
