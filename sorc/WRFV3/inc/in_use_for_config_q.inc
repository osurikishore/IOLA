!STARTOFREGISTRYGENERATEDINCLUDE 'inc/in_use_for_config_q.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
IF(TRIM(vname).EQ.'qshear')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%bl_mynn_tkebudget(id).EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%bl_mynn_tkebudget(id).EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'qbuoy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%bl_mynn_tkebudget(id).EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%bl_mynn_tkebudget(id).EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'qdiss')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%bl_mynn_tkebudget(id).EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%bl_mynn_tkebudget(id).EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'qwt')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%bl_mynn_tkebudget(id).EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%bl_mynn_tkebudget(id).EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'qc_bl')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%icloud_bl.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%icloud_bl.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'qsnowxy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'q2mvxy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'q2mbxy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'qinxy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'qke')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%bl_pbl_physics(id).EQ.5
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%bl_pbl_physics(id).EQ.5
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'qsq')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%bl_pbl_physics(id).EQ.5
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%bl_pbl_physics(id).EQ.5
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'qke')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%bl_pbl_physics(id).EQ.6
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%bl_pbl_physics(id).EQ.6
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'qsq')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%bl_pbl_physics(id).EQ.6
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%bl_pbl_physics(id).EQ.6
  ENDIF
ENDIF
!ENDOFREGISTRYGENERATEDINCLUDE
