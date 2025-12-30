!STARTOFREGISTRYGENERATEDINCLUDE 'inc/in_use_for_config_d.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
IF(TRIM(vname).EQ.'dz_lake3d')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_lake_physics(id).EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_lake_physics(id).EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dz3d')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_lake_physics(id).EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_lake_physics(id).EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'distsq')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'distsq')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.5
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.5
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'distsq')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.6
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.6
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'distsq')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.7
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.7
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_re_cloud')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics_dfi(id).EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics_dfi(id).EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_re_ice')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics_dfi(id).EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics_dfi(id).EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_re_snow')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics_dfi(id).EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics_dfi(id).EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_re_cloud')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics_dfi(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics_dfi(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_re_ice')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics_dfi(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics_dfi(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_re_snow')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics_dfi(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics_dfi(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_re_cloud')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics_dfi(id).EQ.6
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics_dfi(id).EQ.6
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_re_ice')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics_dfi(id).EQ.6
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics_dfi(id).EQ.6
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_re_snow')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics_dfi(id).EQ.6
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics_dfi(id).EQ.6
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_re_cloud')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics_dfi(id).EQ.8
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics_dfi(id).EQ.8
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_re_ice')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics_dfi(id).EQ.8
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics_dfi(id).EQ.8
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_re_snow')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics_dfi(id).EQ.8
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics_dfi(id).EQ.8
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_re_cloud')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics_dfi(id).EQ.28
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics_dfi(id).EQ.28
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_re_ice')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics_dfi(id).EQ.28
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics_dfi(id).EQ.28
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_re_snow')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics_dfi(id).EQ.28
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics_dfi(id).EQ.28
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_re_cloud')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics_dfi(id).EQ.14
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics_dfi(id).EQ.14
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_re_ice')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics_dfi(id).EQ.14
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics_dfi(id).EQ.14
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_re_snow')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics_dfi(id).EQ.14
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics_dfi(id).EQ.14
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_re_cloud')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics_dfi(id).EQ.16
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics_dfi(id).EQ.16
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_re_ice')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics_dfi(id).EQ.16
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics_dfi(id).EQ.16
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_re_snow')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics_dfi(id).EQ.16
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics_dfi(id).EQ.16
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dqke')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%bl_mynn_tkebudget(id).EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%bl_mynn_tkebudget(id).EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'deeprechxy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'ddr_kf')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%kf_edrates(id).EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%kf_edrates(id).EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'der_kf')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%kf_edrates(id).EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%kf_edrates(id).EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfgdp')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%grav_settling(id).EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%grav_settling(id).EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfgdp')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%grav_settling(id).EQ.2
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%grav_settling(id).EQ.2
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_pd')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_pint')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_dwdt')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_t')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_q')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_u')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_v')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_q2')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_cwm')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_rrw')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_stc')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_smc')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_sh2o')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_snow')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_snowh')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_canwat')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_nmm_tsk')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_snowc')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_pd')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.2
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.2
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_pint')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.2
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.2
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_dwdt')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.2
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.2
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_t')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.2
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.2
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_q')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.2
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.2
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_u')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.2
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.2
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_v')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.2
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.2
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_q2')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.2
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.2
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_cwm')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.2
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.2
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_rrw')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.2
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.2
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_stc')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.2
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.2
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_smc')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.2
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.2
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_sh2o')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.2
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.2
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_snow')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.2
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.2
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_snowh')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.2
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.2
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_canwat')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.2
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.2
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_nmm_tsk')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.2
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.2
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_snowc')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.2
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.2
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_pd')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_pint')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_dwdt')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_t')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_q')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_u')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_v')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_q2')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_cwm')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_rrw')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_stc')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_smc')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_sh2o')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_snow')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_snowh')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_canwat')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_nmm_tsk')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'dfi_snowc')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%dfi_opt.EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%dfi_opt.EQ.3
  ENDIF
ENDIF
!ENDOFREGISTRYGENERATEDINCLUDE
