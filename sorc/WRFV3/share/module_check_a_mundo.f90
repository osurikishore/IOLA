

   MODULE module_check_a_mundo


























      USE module_state_description
      USE module_model_constants
      USE module_wrf_error
      USE module_configure

      IMPLICIT NONE



   CONTAINS



   SUBROUTINE check_nml_consistency







      IMPLICIT NONE

      LOGICAL :: exists, vnest
      LOGICAL , EXTERNAL :: wrf_dm_on_monitor
      INTEGER :: i, j, oops, d1_value
      LOGICAL :: km_opt_already_done , diff_opt_already_done
      INTEGER :: count_opt
      LOGICAL :: lon_extent_is_global , lat_extent_is_global





      LOGICAL :: rinblw_already_done



      INTEGER :: count_fatal_error





   count_fatal_error = 0
   model_config_rec % wrf_hydro = 0




   do i=1,model_config_rec%max_dom
      if(model_config_rec%corral_x(i)<5) then
         call wrf_message("Corral X distance must be at least 5 due to intermediate domain halos.")
         model_config_rec%corral_x(i)=5
      endif
      if(model_config_rec%corral_y(i)<5) then
         call wrf_message("Corral Y distance must be at least 5 due to intermediate domain halos.")
         model_config_rec%corral_y(i)=5
      endif
    enddo






      DO i = 1, model_config_rec % max_dom
         IF ( model_config_rec % mp_physics(i) == nssl_2mom .OR. &
              model_config_rec % mp_physics(i) == nssl_2momccn .OR. &
              model_config_rec % mp_physics(i) == nssl_1mom .OR. &
              model_config_rec % mp_physics(i) == nssl_1momlfo .OR. &
              model_config_rec % mp_physics(i) == nssl_2momg ) THEN
            wrf_err_message = '--- ERROR: Chosen microphysics scheme cannot run with WRF-NMM '
            CALL wrf_message ( wrf_err_message )
            wrf_err_message = '--- Fix mp_physics in namelist.input '
         CALL wrf_debug ( 0, TRIM( wrf_err_message ) )
         count_fatal_error = count_fatal_error + 1
         END IF
      ENDDO





      DO i = 2, model_config_rec % max_dom
         IF ( model_config_rec % sf_surface_physics(i)     .NE. &
              model_config_rec % sf_surface_physics(1) ) THEN
            wrf_err_message = '--- ERROR: sf_surface_physics must be equal for all domains '
            CALL wrf_message ( wrf_err_message )
            wrf_err_message = '--- Fix sf_surface_physics in namelist.input '
         CALL wrf_debug ( 0, TRIM( wrf_err_message ) )
         count_fatal_error = count_fatal_error + 1
         END IF
      ENDDO






      DO i = 2, model_config_rec % max_dom
         IF ( model_config_rec % sf_sfclay_physics(i)     .NE. &
              model_config_rec % sf_sfclay_physics(1) ) THEN
            wrf_err_message = '--- ERROR: sf_sfclay_physics must be equal for all domains '
            CALL wrf_message ( wrf_err_message )
            wrf_err_message = '--- Fix sf_sfclay_physics in namelist.input '
         CALL wrf_debug ( 0, TRIM( wrf_err_message ) )
         count_fatal_error = count_fatal_error + 1
         END IF
      ENDDO






      DO i = 2, model_config_rec % max_dom
         IF ( model_config_rec % mp_physics(i)     .NE. &
              model_config_rec % mp_physics(1) ) THEN
            wrf_err_message = '--- NOTE: mp_physics must be equal for all domains '
            CALL wrf_debug ( 1, wrf_err_message )
            wrf_err_message = '--- NOTE:     ----> Setting all mp_physics entries to value defined in the inner most domain'
            CALL wrf_debug ( 1, wrf_err_message )
         END IF
      ENDDO
      d1_value = model_config_rec%mp_physics(model_config_rec % max_dom)
      DO i = 1, model_config_rec % max_dom-1
         model_config_rec%mp_physics(i) = d1_value
      END DO






      DO i = 2, model_config_rec % max_dom
         IF ( model_config_rec % ra_lw_physics(i)     .NE. &
              model_config_rec % ra_lw_physics(1) ) THEN
            wrf_err_message = '--- ERROR: ra_lw_physics must be equal for all domains '
            CALL wrf_message ( wrf_err_message )
            wrf_err_message = '--- Fix ra_lw_physics in namelist.input '
         CALL wrf_debug ( 0, TRIM( wrf_err_message ) )
         count_fatal_error = count_fatal_error + 1
         END IF
      ENDDO

      DO i = 2, model_config_rec % max_dom
         IF ( model_config_rec % ra_sw_physics(i)     .NE. &
              model_config_rec % ra_sw_physics(1) ) THEN
            wrf_err_message = '--- ERROR: ra_sw_physics must be equal for all domains '
            CALL wrf_message ( wrf_err_message )
            wrf_err_message = '--- Fix ra_sw_physics in namelist.input '
         CALL wrf_debug ( 0, TRIM( wrf_err_message ) )
         count_fatal_error = count_fatal_error + 1
         END IF
      ENDDO






         IF ( ( model_config_rec % use_wps_input == 0 ) .AND. &
              ( model_config_rec % time_step .EQ. -1 ) ) THEN

            wrf_err_message = '--- ERROR: Known problem.  time_step must be set to a positive integer'
         CALL wrf_debug ( 0, TRIM( wrf_err_message ) )
         count_fatal_error = count_fatal_error + 1

         END IF





      DO i = 2, model_config_rec % max_dom
         IF ( ( model_config_rec % bl_pbl_physics(i) .NE. model_config_rec % bl_pbl_physics(1) ) .AND. &
              ( model_config_rec % bl_pbl_physics(i) .NE. 0                                    ) ) THEN
            wrf_err_message = '--- ERROR: bl_pbl_physics must be equal for all domains (or = zero)'
            CALL wrf_message ( wrf_err_message )
            wrf_err_message = '--- Fix bl_pbl_physics in namelist.input '
         CALL wrf_debug ( 0, TRIM( wrf_err_message ) )
         count_fatal_error = count_fatal_error + 1
         END IF
      ENDDO







      DO i = 2, model_config_rec % max_dom
         IF ( ( model_config_rec % cu_physics(i) .NE. model_config_rec % cu_physics(1) ) .AND. &
              ( model_config_rec % cu_physics(i) .NE. 0                                ) ) THEN
            wrf_err_message = '--- ERROR: cu_physics must be equal for all domains (or = zero)'
            CALL wrf_message ( wrf_err_message )
            wrf_err_message = '--- Fix cu_physics in namelist.input '
         CALL wrf_debug ( 0, TRIM( wrf_err_message ) )
         count_fatal_error = count_fatal_error + 1
         END IF
      ENDDO







      IF ( ( model_config_rec%fractional_seaice .EQ. 0 ).AND. &
              ( model_config_rec%tice2tsk_if2cold ) ) THEN
            wrf_err_message = '--- WARNING: You set tice2tsk_if2cold = .true.,  but fractional_seaice = 0'
            CALL wrf_debug ( 1, wrf_err_message )
            wrf_err_message = '--- WARNING: tice2tsk_if2cold will have no effect on results.'
            CALL wrf_debug ( 1, wrf_err_message )
      END IF





      DO i = 1, model_config_rec % max_dom
         IF ( ( model_config_rec%fine_input_stream(i) .NE. 0 ).AND. &
              ( model_config_rec%io_form_auxinput2 .EQ. 0 ) ) THEN
            wrf_err_message = '--- ERROR: If fine_input_stream /= 0, io_form_auxinput2 must be /= 0'
            CALL wrf_message ( wrf_err_message )
            wrf_err_message = '--- Set io_form_auxinput2 in the time_control namelist (probably to 2).'
         CALL wrf_debug ( 0, TRIM( wrf_err_message ) )
         count_fatal_error = count_fatal_error + 1
         END IF
      ENDDO





      DO i = 1, model_config_rec % max_dom
         IF ( model_config_rec%sf_surface_physics(i) == NOAHMPSCHEME ) THEN
            WRITE(wrf_err_message, '(" --- ERROR:   Noah-MP LSM scheme (sf_surface_physics==", I2, ")")') NOAHMPSCHEME
            CALL wrf_message ( TRIM ( wrf_err_message ) )
            WRITE(wrf_err_message, '("              does not work with NMM ")')
            CALL wrf_message ( TRIM ( wrf_err_message ) )
            WRITE(wrf_err_message, '("Select a different LSM scheme ")')
         CALL wrf_debug ( 0, TRIM( wrf_err_message ) )
         count_fatal_error = count_fatal_error + 1
         END IF
      END DO




      DO i = 1, model_config_rec % max_dom
         IF ( model_config_rec%shcu_physics(i) == nscvshcuscheme .AND. model_config_rec%cu_physics(i) == nsasscheme) THEN
            WRITE(wrf_err_message, '(" --- ERROR: NSCV shallow convection scheme is already included in NSAS ")')
            CALL wrf_message ( TRIM ( wrf_err_message ) )
         CALL wrf_debug ( 0, TRIM( wrf_err_message ) )
         count_fatal_error = count_fatal_error + 1
         END IF
      END DO







   DO i = 1, model_config_rec % max_dom
         IF ( model_config_rec % sppt(i) .ne. 0)  then
           model_config_rec % sppt_on=1
           IF (( model_config_rec%KMINFORCT .ne. 1) .or. (model_config_rec%KMAXFORCT .ne. 1000000) .or.   &
               ( model_config_rec%LMINFORCT .ne. 1) .or. (model_config_rec%LMAXFORCT .ne. 1000000)) then
               wrf_err_message = '--- Warning: the namelist parameter "kminforct" etc. are for SKEBS only'
               CALL wrf_message ( wrf_err_message )
               wrf_err_message = '             and should not be changed from their default value for SPPT'
               CALL wrf_message ( wrf_err_message )
               wrf_err_message = '--- ERROR: If you really want to modify "kminforct" etc.,  edit module_check a_mundo.'
               CALL wrf_error_fatal3("<stdin>",306,&
TRIM( wrf_err_message ) )
           endif
         endif
   ENDDO
   DO i = 1, model_config_rec % max_dom
         IF ( model_config_rec % rand_perturb(i) .ne. 0)  then
           model_config_rec % rand_perturb_on=1
           IF (( model_config_rec%KMINFORCT .ne. 1) .or. (model_config_rec%KMAXFORCT .ne. 1000000) .or.   &
               ( model_config_rec%LMINFORCT .ne. 1) .or. (model_config_rec%LMAXFORCT .ne. 1000000)) then
               wrf_err_message = '--- Warning: the namelist parameter "kminforct" etc are for SKEBS only'
               CALL wrf_message ( wrf_err_message )
               wrf_err_message = '             and should not be changed from their default value for RAND_PERTURB'
               CALL wrf_message ( wrf_err_message )
               wrf_err_message = ' ABORT. If you really want to modify "kminforct" etc.,  edit module_check a_mundo.'
               CALL wrf_error_fatal3("<stdin>",321,&
TRIM( wrf_err_message ) )
           endif
         endif
   ENDDO
   DO i = 1, model_config_rec % max_dom
         IF ( model_config_rec % stoch_vertstruc_opt(i) ==1 )  then
           model_config_rec % skebs_vertstruc=1       
                                                      
           wrf_err_message = '--- WARNING: the namelist parameter "stoch_vertstruc_opt" is obsolete.'
           CALL wrf_message ( wrf_err_message )
           wrf_err_message = '             Please replace with namelist parameter "skebs_vertstruc" in V3.7 and later versions.'
           CALL wrf_message ( wrf_err_message )
         endif
  ENDDO

   DO i = 1, model_config_rec % max_dom
         IF ( model_config_rec % stoch_force_opt(i) ==1 )  then
           model_config_rec % skebs(i)=1    
                                            
           wrf_err_message = '--- WARNING: the namelist parameter "stoch_force_opt" is obsolete.'
           CALL wrf_message ( wrf_err_message )
           wrf_err_message = '             Please replace with namelist parameter "skebs" in V3.7 and later versions.'
           CALL wrf_message ( wrf_err_message )
         endif
   ENDDO
   DO i = 1, model_config_rec % max_dom
         IF ( model_config_rec % skebs(i) .ne. 0)  then
           model_config_rec % skebs_on=1
         endif
   ENDDO





   DO i = 1, model_config_rec % max_dom
         IF ( model_config_rec % spp_mp(i) .ne. 0)  then
           model_config_rec % spp_on=1
         endif
   ENDDO





   DO i = 1, model_config_rec % max_dom
         IF ((model_config_rec % isppt_mp(i) .ne. 0) .or.    &
             ( model_config_rec % isppt_pbl(i) .ne. 0) .or.   &
             ( model_config_rec % isppt_rad(i) .ne. 0) .or.   &
             ( model_config_rec % isppt_conv(i) .ne. 0 ))  then
               model_config_rec % isppt_on=1
               model_config_rec % sppt_on=0
         endif
   ENDDO





   IF ( model_config_rec % skebs_vertstruc     .ne. 99 )  then
      model_config_rec % num_stoch_levels = model_config_rec %e_vert(1)
   ENDIF
   IF ( model_config_rec % sppt_vertstruc      .ne. 99 )  then
      model_config_rec % num_stoch_levels = model_config_rec %e_vert(1)
   ENDIF
   IF ( model_config_rec % rand_pert_vertstruc .ne. 99 )  then
      model_config_rec % num_stoch_levels = model_config_rec %e_vert(1)
   ENDIF
   model_config_rec % num_stoch_levels = model_config_rec %e_vert(1)

















   IF ( model_config_rec%traj_opt /= 0 ) THEN
         wrf_err_message = 'Trajectories not supported in NMM core '
         CALL wrf_debug ( 0, TRIM( wrf_err_message ) )
         count_fatal_error = count_fatal_error + 1
   END IF








      IF ( model_config_rec%sst_update .EQ. 0 ) THEN
         model_config_rec%io_form_auxinput4 = 0
         DO i = 1, model_config_rec % max_dom
            wrf_err_message = '--- NOTE: sst_update is 0, ' // &
                  'setting io_form_auxinput4 = 0 and auxinput4_interval = 0 for all domains'
            CALL wrf_debug ( 1, TRIM( wrf_err_message ) )
            model_config_rec%auxinput4_interval(i)   = 0
            model_config_rec%auxinput4_interval_y(i) = 0
            model_config_rec%auxinput4_interval_d(i) = 0
            model_config_rec%auxinput4_interval_h(i) = 0
            model_config_rec%auxinput4_interval_m(i) = 0
            model_config_rec%auxinput4_interval_s(i) = 0
         ENDDO
      ELSE
         IF ( model_config_rec%io_form_auxinput4 .EQ. 0 ) THEN
            wrf_err_message = '--- ERROR: If sst_update /= 0, io_form_auxinput4 must be /= 0'
            CALL wrf_message ( wrf_err_message )
            wrf_err_message = '--- Set io_form_auxinput4 in the time_control namelist (probably to 2).'
         CALL wrf_debug ( 0, TRIM( wrf_err_message ) )
         count_fatal_error = count_fatal_error + 1
         END IF
      END IF







      IF ( model_config_rec%sst_update .EQ. 1 ) THEN
         IF ( model_config_rec%io_form_auxinput4 .EQ. 0 ) THEN
            wrf_err_message = '--- ERROR: If sst_update /= 0, io_form_auxinput4 must be /= 0'
            CALL wrf_debug ( 0, TRIM(wrf_err_message) )
            wrf_err_message = '--- Set io_form_auxinput4 in the time_control namelist (probably to 2).'
            CALL wrf_debug ( 0, TRIM( wrf_err_message ) )
            count_fatal_error = count_fatal_error + 1
         END IF

         IF ( ( model_config_rec%auxinput4_interval(1)   .EQ. 0 ) .AND. &
              ( model_config_rec%auxinput4_interval_y(1) .EQ. 0 ) .AND. &
              ( model_config_rec%auxinput4_interval_d(1) .EQ. 0 ) .AND. &
              ( model_config_rec%auxinput4_interval_h(1) .EQ. 0 ) .AND. &
              ( model_config_rec%auxinput4_interval_m(1) .EQ. 0 ) .AND. &
              ( model_config_rec%auxinput4_interval_s(1) .EQ. 0 ) ) THEN
            wrf_err_message = '--- ERROR: If sst_update /= 0, one of the auxinput4_interval settings must be /= 0'
            CALL wrf_debug ( 0, TRIM(wrf_err_message) )
            wrf_err_message = '--- Set auxinput4_interval_s to the same value as interval_seconds (usually a pretty good guess).'
            CALL wrf_debug ( 0, TRIM( wrf_err_message ) )
            count_fatal_error = count_fatal_error + 1
         END IF
      END IF





















      IF ( ( model_config_rec % ra_lw_physics(1) .EQ. RRTMG_LWSCHEME )  .OR. &
           ( model_config_rec % ra_sw_physics(1) .EQ. RRTMG_SWSCHEME )  .OR. &
           ( model_config_rec % ra_lw_physics(1) .EQ. RRTMK_LWSCHEME )  .OR. &
           ( model_config_rec % ra_lw_physics(1) .EQ. RRTMK_SWSCHEME )  .OR. &
           ( model_config_rec % ra_lw_physics(1) .EQ. RRTMG_LWSCHEME_FAST )  .OR. &
           ( model_config_rec % ra_sw_physics(1) .EQ. RRTMG_SWSCHEME_FAST )  ) THEN
         wrf_err_message = '--- NOTE: RRTMG radiation is used, namelist ' // &
                           'value for o3input (ozone input) is used '
         CALL wrf_debug ( 1, wrf_err_message )
      ELSE
         model_config_rec % o3input = 0
         wrf_err_message = '--- NOTE: RRTMG radiation is not used, setting:  ' // &
                           'o3input=0 to avoid data pre-processing'
         CALL wrf_debug ( 1, wrf_err_message )
      END IF











      IF ( model_config_rec%clean_atm_diag > 0 ) THEN

         wrf_err_message = '--- NOTE: "Clean" atmosphere diagnostics can only be used in WRF-Chem' 
         CALL wrf_debug ( 0, TRIM( wrf_err_message ) )
         model_config_rec%calc_clean_atm_diag = 0

      ENDIF






     IF ( model_config_rec%idcor .lt. 0 .OR. model_config_rec%idcor .gt. 1) THEN
         wrf_err_message = '--- ERROR: idcor value must be between 0 and 1'
         CALL wrf_message( wrf_err_message )
         count_fatal_error = count_fatal_error + 1
     END IF





      IF ( count_fatal_error .GT. 0 ) THEN
         WRITE (wrf_err_message, FMT='(A,I6, A)') 'NOTE:  ', count_fatal_error, &
                                            ' namelist settings are wrong. Please check and reset these options'
         CALL wrf_error_fatal3("<stdin>",545,&
wrf_err_message  )
      END IF

   END SUBROUTINE check_nml_consistency



   SUBROUTINE setup_physics_suite










      USE module_domain, ONLY : change_to_lower_case

      IMPLICIT NONE

   END SUBROUTINE setup_physics_suite



   SUBROUTINE set_physics_rconfigs










      IMPLICIT NONE

      INTEGER :: numsoiltemp , nummosaictemp
      INTEGER :: i






      IF ( model_config_rec % sf_surface_mosaic .EQ. 1 ) THEN
      
      numsoiltemp = model_config_rec % num_soil_layers
      nummosaictemp = model_config_rec % mosaic_cat
      
         model_config_rec % mosaic_cat_soil = numsoiltemp * nummosaictemp

         wrf_err_message = '--- NOTE: Noah-mosaic is in use, setting:  ' // &
                           'mosaic_cat_soil = mosaic_cat * num_soil_layers'
         CALL wrf_debug ( 1, wrf_err_message )

      END IF     
      

      





      IF (( model_config_rec % ra_lw_physics(1) .EQ. CAMLWSCHEME ) .OR. &
          ( model_config_rec % ra_sw_physics(1) .EQ. CAMSWSCHEME )) THEN
         model_config_rec % paerlev = 29
         model_config_rec % levsiz = 59
         model_config_rec % cam_abs_dim1 = 4
         model_config_rec % cam_abs_dim2 = model_config_rec % e_vert(1)

         wrf_err_message = '--- NOTE: CAM radiation is in use, setting:  ' // &
                           'paerlev=29, levsiz=59, cam_abs_dim1=4, cam_abs_dim2=e_vert'
         CALL wrf_debug ( 1, wrf_err_message )

      END IF
      






      DO i = 1, model_config_rec % max_dom
         IF ( ( model_config_rec % mp_physics(i) .EQ. MILBRANDT2MOM ) .OR. &
              ( model_config_rec % do_radar_ref  .EQ. 1             ) ) THEN
            model_config_rec % compute_radar_ref = 1
         END IF
      ENDDO












      IF (( model_config_rec % ra_lw_physics(1) .EQ. RRTMG_LWSCHEME      ) .OR. &
          ( model_config_rec % ra_sw_physics(1) .EQ. RRTMG_SWSCHEME      ) .OR. &
          ( model_config_rec % ra_lw_physics(1) .EQ. RRTMK_LWSCHEME      ) .OR. &
          ( model_config_rec % ra_sw_physics(1) .EQ. RRTMK_SWSCHEME      ) .OR. &
          ( model_config_rec % ra_lw_physics(1) .EQ. RRTMG_LWSCHEME_FAST ) .OR. &
          ( model_config_rec % ra_sw_physics(1) .EQ. RRTMG_SWSCHEME_FAST )) THEN
         model_config_rec % levsiz = 59
         model_config_rec % alevsiz = 12
         model_config_rec % no_src_types = 6

         wrf_err_message = '--- NOTE: One of the RRTMG radiation schemes is in use, setting:  ' // &
                           'levsiz=59, alevsiz=12, no_src_types=6'
         CALL wrf_debug ( 1, wrf_err_message )

      END IF






      IF      (   model_config_rec % sf_surface_physics(1) .EQ. NOLSMSCHEME  ) THEN
         model_config_rec % num_soil_layers = 5
      ELSE IF (   model_config_rec % sf_surface_physics(1) .EQ. SLABSCHEME   ) THEN
         model_config_rec % num_soil_layers = 5
      ELSE IF (   model_config_rec % sf_surface_physics(1) .EQ. LSMSCHEME    ) THEN
         model_config_rec % num_soil_layers = 4
      ELSE IF (   model_config_rec % sf_surface_physics(1) .EQ. NOAHMPSCHEME ) THEN
         model_config_rec % num_soil_layers = 4
      ELSE IF ( ( model_config_rec % sf_surface_physics(1) .EQ. RUCLSMSCHEME ) .AND. &
                ( model_config_rec % num_soil_layers .EQ. 6 ) ) THEN
         model_config_rec % num_soil_layers = 6
      ELSE IF ( ( model_config_rec % sf_surface_physics(1) .EQ. RUCLSMSCHEME ) .AND. &
                ( model_config_rec % num_soil_layers .EQ. 9 ) ) THEN
         model_config_rec % num_soil_layers = 9
      ELSE IF (   model_config_rec % sf_surface_physics(1) .EQ. RUCLSMSCHEME ) THEN
         model_config_rec % num_soil_layers = 6
      ELSE IF (   model_config_rec % sf_surface_physics(1) .EQ. PXLSMSCHEME  ) THEN
         model_config_rec % num_soil_layers = 2
      ELSE IF (   model_config_rec % sf_surface_physics(1) .EQ. CLMSCHEME    ) THEN
         model_config_rec % num_soil_layers = 10
      ELSE IF (   model_config_rec % sf_surface_physics(1) .EQ. SSIBSCHEME   ) THEN
         model_config_rec % num_soil_layers = 3
      ELSE IF (   model_config_rec % sf_surface_physics(1) .EQ. GFDLSLAB     ) THEN
         model_config_rec % num_soil_layers = 4
      ELSE
         CALL wrf_debug       ( 0 , '--- ERROR: Unknown sf_surface_physics has no associated number of soil levels' )
         WRITE (wrf_err_message, FMT='(A,I6)') '--- ERROR: sf_surface_physics = ' , model_config_rec % sf_surface_physics(1)
         CALL wrf_error_fatal3("<stdin>",698,&
TRIM(wrf_err_message) )
      END IF 

      WRITE (wrf_err_message, FMT='(A,I6)') '--- NOTE: num_soil_layers has been set to ', &
                                             model_config_rec % num_soil_layers
      CALL wrf_debug ( 1, wrf_err_message )

   END SUBROUTINE set_physics_rconfigs




   END MODULE module_check_a_mundo


