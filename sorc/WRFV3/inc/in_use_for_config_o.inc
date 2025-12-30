!STARTOFREGISTRYGENERATEDINCLUDE 'inc/in_use_for_config_o.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
IF(TRIM(vname).EQ.'o3rad')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%ra_lw_physics(id).EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%ra_lw_physics(id).EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'o3rad')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%ra_lw_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%ra_lw_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'o3rad')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%ra_lw_physics(id).EQ.98
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%ra_lw_physics(id).EQ.98
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'o3rad')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%ra_sw_physics(id).EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%ra_sw_physics(id).EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'o3rad')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%ra_sw_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%ra_sw_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'o3rad')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%ra_sw_physics(id).EQ.98
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%ra_sw_physics(id).EQ.98
  ENDIF
ENDIF
!ENDOFREGISTRYGENERATEDINCLUDE
