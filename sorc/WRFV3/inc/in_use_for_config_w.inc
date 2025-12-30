!STARTOFREGISTRYGENERATEDINCLUDE 'inc/in_use_for_config_w.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
IF(TRIM(vname).EQ.'watsat3d')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_lake_physics(id).EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_lake_physics(id).EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'windsq_swath')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%swath_mode.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%swath_mode.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'windsq_swath')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%swath_mode.EQ.2
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%swath_mode.EQ.2
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'wind10_ratio')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%swath_mode.EQ.2
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%swath_mode.EQ.2
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'weightout')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'wslakexy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'waxy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'wtxy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'woodxy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'wgapxy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
!ENDOFREGISTRYGENERATEDINCLUDE
