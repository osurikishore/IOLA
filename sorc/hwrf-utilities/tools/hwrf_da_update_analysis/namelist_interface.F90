module namelist_interface

  !=======================================================================

  !$$$ MODULE DOCUMENTATION BLOCK
  
  ! da-update-analysis :: namelist_interface
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

  ! * datapath; a FORTRAN character string specifying the full-path to
  !   the directory to contain output files written by the respective
  !   routines.

  ! * debug; a FORTRAN logical value specifying whether to include
  !   debug information during execution.

  ! * hwrf_an_filename; a FORTRAN character string specifying the
  !   full-path to the Hurricane WRF (HWRF) analysis variable file
  !   (e.g., the file updated by the respective data-assimilation
  !   methods).

  ! * hwrf_bdy_filename; a FORTRAN character string specifying the
  !   full-path to the Hurricane WRF (HWRF) boundary condition file.

  ! * hwrf_bdy_json_vtable; a FORTRAN character string specifying the
  !   full-path to the JSON-formatted for the Hurricane WRF (HWRF)
  !   prognostic state variables (and tendencies) to be updated using
  !   the boundary update methods within.

  ! * hwrf_fg_filename; a FORTRAN character string specifying the
  !   full-path to the Hurricane WRF (HWRF) first-guess variable file
  !   (e.g., the file containing the background-forecast for the
  !   respective data-assimilation methods).

  ! * hwrf_grid_json_vtable; a FORTRAN character string specifying the
  !   full-path to the JSON-formatted file for the HWRF variables
  !   describing the respective grid projections be operated upon.

  ! * hwrf_nest_filename; a FORTRAN character string specifying the
  !   full-path to the Hurricane WRF (HWRF) variable file for the
  !   respective nest region.

  ! * hwrf_parent_filename; a FORTRAN character string specifying the
  !   full-path to the Hurricane WRF (HWRF) variable file for the
  !   parent region within which a nest (hwrf_nest_filename) is to be
  !   embedded.

  ! * hwrf_var_json_vtable; a FORTRAN character string specifying the
  !   full-path to the JSON-formatted file for the HWRF variables to
  !   be operated upon.

  ! * is_bdy_update; a FORTRAN logical variable specifying whether to
  !   update the lateral boundary condition prognostic state variables
  !   and/or tendencies.

  ! * is_blend_env; a FORTRAN logical variable specifying whether to
  !   blend the data assimilation analysis increments with the
  !   environment (e.g., the first-guess) at a specified radius (see
  !   tc_large_scale_radius).

  ! * ic_blend_opt; a FORTRAN integer value variable specifying the
  !   blending option to be used for tropical cyclone (TC) inner-core
  !   increments; the following options are currently supported:

  !   0: use the full analysis increments.

  !   1: wave-number decomposition and reconstruction of analysis
  !   increments up to and including the maximum wave-number (see
  !   wndcmp_mxwvn) specified by the user.

  !   2: set all analysis increments within the user specified value
  !   for tc_ic_radius to 0.0.  

  ! * is_da; a FORTRAN logical variable specifying whether the
  !   algorithms are to be applied to the respective forecast model's
  !   data-assimilation defined variable states; currently only the
  !   Grid-point Statistical Interpolation (GSI) is supported.
  
  ! * is_ghst; a FORTRAN logical variable specifying whether the a
  !   'ghost domain' (e.g., for HWRF, this is the data-assimilation
  !   domain) is to be used to update/define a forecast model input
  !   file (e.g., for HWRF, wrfinput_d0*).

  ! * is_hwrf; a FORTRAN logical variable specifying whether the
  !   forecast model is the Hurricane WRF (HWRF).

  ! * is_llp; a FORTRAN logical variable specifying whether the
  !   vertical interpolation is performed using linear-log of pressure
  !   interpolation.

  ! * is_merge; a FORTRAN logical variable specifying whether the
  !   algorithms are to be applied to merge user specified variables
  !   within identical dimension forecast domains.

  ! * is_merge_tc; a FORTRAN logical variable specifying whether the
  !   algorithms are to be applied to merge user specified variables
  !   relative to records within the TC-vitals (see tcv_filename).

  ! * is_nest; a FORTRAN logical variable specifying whether the
  !   algorithms are to applied to blend a nested grid region within
  !   it's respective parent grid.

  ! * is_nest_cycle; a FORTRAN logical variable specifying whether the
  !   algorithms are to be applied to cycle respective forecast cycle
  !   nest's from (subsequent) forecast model cycles.

  ! * is_slint; a FORTRAN logical variable specifying whether the
  !   remapping coefficients are defined via the SLINT interpolation
  !   algorithm (Wang, 2006); currently this is the default and is not
  !   able to be modified by the user.

  ! * is_wnblnd; a FORTRAN logical variable specifying whether to
  !   apply a wave-number filtering for the analysis increments.

  ! * npts_blend; a FORTRAN integer value variable specifying the
  !   total number of point to be used for blending boundary values;
  !   this value is the total number which includes the relaxation
  !   values (see npts_relax).

  ! * npts_relax; a FORTRAN integer value variable specifying the
  !   total number of point to be used for relaxing values within the
  !   blending boundary values.

  ! * region_blend_ratio; a FORTRAN 4-byte real valued variable
  !   specifying the blending regions used to interpolate between a
  !   respective tropical cyclone (TC) inner-core region and
  !   environmental region (see tc_ic_radius and tc_env_radius below,
  !   respectively); for example, a value of 0.25 implies shrinking
  !   the inner-core radius by 25% and increasing the environment
  !   radius by 25% in order to provide enough values for a reasonable
  !   bi-linear interpolation of values; a smaller value means less
  !   points provided to the remapping algorithm (but potentially
  !   faster computations) while a larger value provides more points
  !   to the remapping algorithm at the cost of increased (e.g.,
  !   slower) computations.

  ! * tc_env_radius; a FORTRAN 4-byte real valued variable defining
  !   the radius, relative to a respective tropical cyclone (TC)
  !   position, of the environment surrounding the TC; this is used
  !   for the increment blending algorithms; units are meters.

  ! * tc_ic_radius; a FORTRAN 4-byte real valued variable defining the
  !   radius, relative to a respective tropical cyclone (TC) position,
  !   of the inner-core surrounding the TC; this is used for the
  !   increment blending algorithms; units are meters.

  ! * tc_large_scale_radius; a FORTRAN 4-byte real valued variable
  !   defining the radius, relative to a respective tropical cyclone
  !   (TC) position, at which to blend the analysis (determined from
  !   the updated data assimilation analysis increments) with the
  !   first-guess analysis variable; units are meters.

  ! * tc_nest_ratio_blend; a FORTRAN 4-byte real valued variable
  !   specifying the ratio of the blended nest and parent domains
  !   which contains the nest versus the parent domain; for example, a
  !   value of 0.95 implies that 95% of the total blended region is
  !   defined by the interpolated nested domain while the remaining 5%
  !   is a linear combination (as a function of radius) of the nest
  !   and parent grid solutions.
  
  ! * tc_region_area; a FORTRAN 4-byte real valued variable specifying
  !   the size of the region, relative to a TC-vitals record, which to
  !   merge into a different (i.e., parent) domain; units are degrees.

  ! * tcv_filename; a FORTRAN character string specifying the
  !   full-path to the TC-vitals formatted file.
  
  ! * use_smthr; a FORTRAN logical variable specifying whether to use
  !   a 9-point spatial smoother within blending regions.

  ! * wndcmp_area; a FORTRAN 4-byte real value variable specifying the
  !   area, relative to a respective tropical cyclone (TC), to be used
  !   for the wave-number decomposition; this value is only used when
  !   ic_blend_opt (see below) is equal to 1; units are degrees.

  ! * wndcmp_dangle; a FORTRAN 4-byte real value variable specifying
  !   azimuthal increment to be used for preparing the grid passed to
  !   the Fast-Fourier Transform (FFT) routines for the analysis
  !   increment wave-number decomposition; this value is only used
  !   when ic_blend_opt (see below) is equal to 1; units are degrees.

  ! * wndcmp_darea; a FORTRAN 4-byte real value variable specifying
  !   the areal spacing to be used for the wave-number decomposition;
  !   this value is only used when ic_blend_opt (see below) is equal
  !   to 1; units are degrees.

  ! * wndcmp_dradius; a FORTRAN 4-byte real value variable specifying
  !   radial increment to be used for preparing the grid passed to the
  !   Fast-Fourier Transform (FFT) routines for the analysis increment
  !   wave-number decomposition; this value is only used when
  !   ic_blend_opt (see below) is equal to 1; units are meters.

  ! * wndcmp_mxwvn; a FORTRAN integer value variable specifying the
  !   total number wave-numbers to retain (beginning with 0) for the
  !   tropical cyclone (TC) inner-core increment reconstruction; this
  !   value is only used when ic_blend_opt (see below) is equal to 1.

  !-----------------------------------------------------------------------

  ! Define local variables

  character(len=500)                                                    :: &
       & datapath = './'
  character(len=500)                                                    :: &
       & hwrf_an_filename = 'NOT USED'
  character(len=500)                                                    :: &
       & hwrf_bdy_filename = 'NOT USED'
  character(len=500)                                                    :: &
       & hwrf_bdy_json_vtable = 'NOT USED'
  character(len=500)                                                    :: &
       & hwrf_grid_json_vtable = 'NOT USED'
  character(len=500)                                                    :: &
       & hwrf_var_json_vtable = 'NOT USED'
  character(len=500)                                                    :: &
       & hwrf_fg_filename = 'NOT USED'
  character(len=500)                                                    :: &
       & hwrf_nest_filename = 'NOT USED'
  character(len=500)                                                    :: &
       & hwrf_parent_filename = 'NOT USED'
  character(len=500)                                                    :: &
       & tcv_filename = 'NOT USED'
  logical                                                               :: &
       & debug = .false.
  logical                                                               :: &
       & is_blend_env = .false.
  logical                                                               :: &
       & is_bdy_update = .false.
  logical                                                               :: &
       & is_da = .false.
  logical                                                               :: &
       & is_ghst = .false.
  logical                                                               :: &
       & is_hwrf = .false.
  logical                                                               :: &
       & is_llp = .false.
  logical                                                               :: &
       & is_merge = .false.
  logical                                                               :: &
       & is_merge_tc = .false.
  logical                                                               :: &
       & is_nest = .false.
  logical                                                               :: &
       & is_nest_cycle = .false.
  logical                                                               :: &
       & is_slint = .true. ! FOR NOW, THIS IS THE ONLY SUPPORTED
                           ! INTERPOLATION TYPE
  logical                                                               :: &
       & is_wnblnd = .false.
  logical                                                               :: &
       & use_smthr = .false.
  real(r_kind)                                                          :: &
       & interp_dradius = 25000.0
  real(r_kind)                                                          :: &
       & region_blend_ratio = 0.25
  real(r_kind)                                                          :: &
       & tc_env_radius = 600000.0
  real(r_kind)                                                          :: &
       & tc_ic_radius = 150000.0
  real(r_kind)                                                          :: &
       & tc_large_scale_radius = 1200000.0
  real(r_kind)                                                          :: &
       & tc_nest_ratio_blend = 1.0
  real(r_kind)                                                          :: &
       & tc_region_area = 11.0
  real(r_kind)                                                          :: &
       & wndcmp_area = 5.0
  real(r_kind)                                                          :: &
       & wndcmp_dangle = 15.0
  real(r_kind)                                                          :: &
       & wndcmp_darea = 0.1
  real(r_kind)                                                          :: &
       & wndcmp_dradius = 25000.0
  integer                                                               :: &
       & ic_blend_opt = 1
  integer                                                               :: &
       & npts_blend = 0
  integer                                                               :: &
       & npts_relax = 0 
  integer                                                               :: &
       & wndcmp_mxwvn = 1
  namelist /share/    debug, datapath, is_bdy_update, is_da, is_ghst,      &
       & is_hwrf, is_merge, is_nest, is_nest_cycle
  namelist /hwrf/     hwrf_an_filename, hwrf_bdy_filename,                 &
       & hwrf_bdy_json_vtable, hwrf_fg_filename, hwrf_grid_json_vtable,    &
       & hwrf_nest_filename, hwrf_parent_filename, hwrf_var_json_vtable
  namelist /incremnt/ ic_blend_opt, is_wnblnd
  namelist /interp/   is_llp, is_slint 
  namelist /remap/    npts_blend, npts_relax, region_blend_ratio,          &
       & use_smthr
  namelist /tc/       is_blend_env, tcv_filename, tc_ic_radius,            &
       & tc_env_radius, tc_large_scale_radius, tc_nest_ratio_blend,        &
       & tc_region_area
  namelist /wndcmp/   wndcmp_area, wndcmp_darea, wndcmp_dangle,            &
       & wndcmp_dradius, wndcmp_mxwvn

  !-----------------------------------------------------------------------

contains
  
  !=======================================================================

  ! SUBROUTINE: 

  ! namelist.f90

  ! DESCRIPTION: 

  ! This subroutine acts as the interface to the namelist file,
  ! provided as 'da-update-analysis.input' by the user.

  !-----------------------------------------------------------------------

  subroutine namelist()

    ! Define variables computed within routine

    character(len=500)                                                  :: nml_filename
    logical                                                             :: is_it_there
    integer                                                             :: unit_nml

    !=====================================================================
    
    ! Define local variables

    nml_filename = './da-update-analysis.input'
    unit_nml     = 9
    is_it_there  = .false.
    inquire(file = trim(adjustl(nml_filename)),exist = is_it_there)

    ! Check local variable and proceed accordingly

    if(is_it_there) then

       ! Define local variables

       open(file   = trim(adjustl(nml_filename)),                          &
            unit   = unit_nml        ,                                     &
            status = 'old'         ,                                       &
            form   = 'formatted'     ,                                     &
            action = 'read')
       read(unit_nml,NML = share)
       read(unit_nml,NML = hwrf)
       read(unit_nml,NML = incremnt)
       read(unit_nml,NML = interp)
       read(unit_nml,NML = remap)
       read(unit_nml,NML = tc)
       read(unit_nml,NML = wndcmp)
       close(unit_nml)

    else  ! if(is_it_there)

       ! Define local variables

       write(6,500) trim(adjustl(nml_filename))
       stop(99)
       
    end if ! if(is_it_there)
       
    ! Define local variables
    
    write(6,*) '&SHARE'                 
    write(6,*) 'DATAPATH                      = ',                         &
         & trim(adjustl(datapath))
    write(6,*) 'DEBUG                         = ', debug
    write(6,*) 'IS_BDY_UPDATE                 = ', is_bdy_update
    write(6,*) 'IS_DA                         = ', is_da
    write(6,*) 'IS_GHST                       = ', is_ghst
    write(6,*) 'IS_HWRF                       = ', is_hwrf
    write(6,*) 'IS_NEST                       = ', is_nest
    write(6,*) 'IS_NEST_CYCLE                 = ', is_nest_cycle
    write(6,*) '/'
    write(6,*) '&HWRF'
    write(6,*) 'HWRF_AN_FILENAME              = ',                         &
         & trim(adjustl(hwrf_an_filename))
    write(6,*) 'HWRF_BDY_FILENAME             = ',                         &
         & trim(adjustl(hwrf_bdy_filename))
    write(6,*) 'HWRF_BDY_JSON_VTABLE          = ',                         &
         & trim(adjustl(hwrf_bdy_json_vtable))
    write(6,*) 'HWRF_FG_FILENAME              = ',                         &
         & trim(adjustl(hwrf_fg_filename))
    write(6,*) 'HWRF_GRID_JSON_VTABLE         = ',                         &
         & trim(adjustl(hwrf_grid_json_vtable))
    write(6,*) 'HWRF_NEST_FILENAME            = ',                         &
         & trim(adjustl(hwrf_nest_filename))
    write(6,*) 'HWRF_PARENT_FILENAME          = ',                         &
         & trim(adjustl(hwrf_parent_filename))
    write(6,*) 'HWRF_VAR_JSON_VTABLE          = ',                         &
         & trim(adjustl(hwrf_var_json_vtable))
    write(6,*) '/'
    write(6,*) '&INCREMNT'
    write(6,*) 'IC_BLEND_OPT                  = ', ic_blend_opt
    write(6,*) 'IS_WNBLND                     = ', is_wnblnd
    write(6,*) '/' 
    write(6,*) '&INTERP'    
    write(6,*) 'IS_SLINT                      = ', is_slint
    write(6,*) 'IS_LLP                        = ', is_llp
    write(6,*) '/'
    write(6,*) '&REMAP'
    write(6,*) 'NPTS_BLEND                    = ', npts_blend
    write(6,*) 'NPTS_RELAX                    = ', npts_relax
    write(6,*) 'REGION_BLEND_RATIO            = ', region_blend_ratio
    write(6,*) '/'  
    write(6,*) '&TC'  
    write(6,*) 'IS_BLEND_ENV                  = ', is_blend_env
    write(6,*) 'TCV_FILENAME                  = ',                         &
         & trim(adjustl(tcv_filename))
    write(6,*) 'TC_ENV_RADIUS                 = ', tc_env_radius
    write(6,*) 'TC_IC_RADIUS                  = ', tc_ic_radius
    write(6,*) 'TC_LARGE_SCALE_RADIUS         = ', tc_large_scale_radius
    write(6,*) 'TC_NEST_RATIO_BLEND           = ', tc_nest_ratio_blend
    write(6,*) '/'
    write(6,*) '&WNDCMP'
    write(6,*) 'WNDCMP_AREA                   = ', wndcmp_area
    write(6,*) 'WNDCMP_DRADIUS                = ', wndcmp_dradius
    write(6,*) 'WNDCMP_DANGLE                 = ', wndcmp_dangle
    write(6,*) 'WNDCMP_MXWVN                  = ', wndcmp_mxwvn
    write(6,*) '/'
500 format('NAMELISTPARAMS: ', a, ' not found in the current working ',    &
         & 'directory. ABORTING!!!!')

    !=====================================================================

  end subroutine namelist

  !=======================================================================

end module namelist_interface
