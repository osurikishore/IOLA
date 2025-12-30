!STARTOFREGISTRYGENERATEDINCLUDE 'inc/in_use_for_config_g.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
IF(TRIM(vname).EQ.'gppxy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'ghvxy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'ghbxy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'gddxy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'grainxy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'gecros_state')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'gd_cloud')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%cu_physics(id).EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%cu_physics(id).EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'gd_cloud2')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%cu_physics(id).EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%cu_physics(id).EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'gd_cloud')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%cu_diag(id).EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%cu_diag(id).EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'gd_cloud2')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%cu_diag(id).EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%cu_diag(id).EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'gd_cldfr')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%cu_diag(id).EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%cu_diag(id).EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'gd_cloud_a')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%cu_diag(id).EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%cu_diag(id).EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'gd_cloud2_a')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%cu_diag(id).EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%cu_diag(id).EQ.1
  ENDIF
ENDIF
!ENDOFREGISTRYGENERATEDINCLUDE
