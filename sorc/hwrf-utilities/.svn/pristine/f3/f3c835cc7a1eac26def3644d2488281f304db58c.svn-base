module namelist_interface

  !=======================================================================

  !$$$ MODULE DOCUMENTATION BLOCK
  
  ! obs_to_bufr :: namelist_interface
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

  !-----------------------------------------------------------------------

  ! DESCRIPTION (alphabetized):

  ! * analdate; a FORTRAN character string specifying the analysis
  !   date about which to define the observation times; formatted as,
  !   assuming UNIX convention, ccyy-mm-dd_HH:MM:SS.

  ! * bufr_filepath; a FORTRAN character string specifying the
  !   full-path to the BUFR file to be written (or appended) to.

  ! * bufr_tblpath; a FORTRAN character string specifying the
  !   full-path to the external file containing the BUFR table to be
  !   written (or appended) to the output BUFR file.  error(s) for the
  !   respective observation type.

  ! * datapath; a FORTRAN character string specifying the full-path to
  !   the directory to contain output files written by the respective
  !   routines.

  ! * debug; a FORTRAN logical value specifying whether to include
  !   debug information during execution.

  ! * hsa; a FORTRAN logical value specifying whether the input data
  !   type is National Oceanic and Atmospheric Administration (NOAA)
  !   Atlantic Oceanographic and Meteorological Laboratory (AOML)
  !   Hurricane Research Division (HRD) spline analysis (HSA) format.

  ! * intrp_obserr; a FORTRAN logical value specifying whether to
  !   interpolate the user specified errors to the isobaric levels of
  !   the respective observations.

  ! * obs_filepath; a FORTRAN character string specifying the
  !   full-path to the external file containing a list of observation
  !   files to process; the respective files must comply with the
  !   logical variable indicating the input data type.

  ! * obserr_filepath; a FORTRAN character string specifying the
  !   full-path to the external file containing the respective
  !   observation-type errors; see the README files accompanying the
  !   software package.

  !-----------------------------------------------------------------------

  ! Define local variables

  character(len=500)                                                    :: &
       & bufr_filepath = 'NOT USED'
  character(len=500)                                                    :: &
       & bufr_tblpath = 'NOT USED'  
  character(len=500)                                                    :: &
       & datapath = './'
  character(len=500)                                                    :: &
       & obs_filepath = 'NOT USED'
  character(len=500)                                                    :: &
       & obserr_filepath = 'NOT USED'
  character(len=19)                                                     :: &
       & analdate = '2000-01-01_00:00:00'
  logical                                                               :: &
       & debug = .false.
  logical                                                               :: &
       & hsa = .false.
  logical                                                               :: &
       & intrp_obserr = .false.
  logical                                                               :: &
       & obs_flag = .false.  
  integer                                                               :: &
       & obs_type_flag(100) = -99
  namelist /share/ analdate, datapath, debug, obs_flag
  namelist /bufr/  bufr_filepath, bufr_tblpath
  namelist /obs/   hsa, intrp_obserr, obs_filepath, obserr_filepath
  !namelist /flag/  obs_type_flag
 
  !-----------------------------------------------------------------------

contains
  
  !=======================================================================

  ! SUBROUTINE: 

  ! namelist.f90

  ! DESCRIPTION: 

  ! This subroutine acts as the interface to the namelist file,
  ! provided as 'obs-to-bufr.input' by the user.

  !-----------------------------------------------------------------------

  subroutine namelist()

    ! Define variables computed within routine

    logical                                                             :: is_it_there
    integer                                                             :: unit_nml
   
    !=====================================================================

    ! Define local variables

    unit_nml     = 9
    is_it_there  = .false.
    inquire(file = 'obs-to-bufr.input',exist = is_it_there)

    ! Check local variable and proceed accordingly

    if(is_it_there) then

       ! Define local variables

       open(file   = 'obs-to-bufr.input',                                  &
            unit   = unit_nml        ,                                     &
            status = 'old'         ,                                       &
            form   = 'formatted'     ,                                     &
            action = 'read')
       read(unit_nml,NML = share)
       read(unit_nml,NML = bufr)
       read(unit_nml,NML = obs)
    !   read(unit_nml,NML = flag)
       close(unit_nml)

    else  ! if(is_it_there)

       ! Define local variables

       write(6,500)
       stop
       
    end if ! if(is_it_there)
       
    ! Define local variables
    
    write(6,*) '&SHARE'             
    write(6,*) 'DATAPATH                      = ',                         &
         & trim(adjustl(datapath))
    write(6,*) 'DEBUG                         = ', debug
    write(6,*) 'ANALDATE                      = ', analdate
    write(6,*) '/'
    write(6,*) '&BUFR'    
    write(6,*) 'BUFR_FILEPATH                 = ',                         &
         & trim(adjustl(bufr_filepath))
    write(6,*) 'BUFR_TBLPATH                  = ',                         &
         & trim(adjustl(bufr_tblpath))
    write(6,*) '/'
    write(6,*) '&OBS'  
    write(6,*) 'HSA                           = ', hsa
    write(6,*) 'INTRP_OBSERR                  = ', intrp_obserr
    write(6,*) 'OBS_FILEPATH                  = ',                         &
         & trim(adjustl(obs_filepath))
    write(6,*) 'OBSERR_FILEPATH               = ',                         &
         & trim(adjustl(obserr_filepath))
    write(6,*) '/'
500 format('NAMELISTPARAMS: obs-to-bufr.input not found in the ',          &
         & 'current working directory. ABORTING!!!!')

    !=====================================================================

  end subroutine namelist
  
  !=======================================================================

end module namelist_interface
