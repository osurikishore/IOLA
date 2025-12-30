!STARTOFREGISTRYGENERATEDINCLUDE 'inc/in_use_for_config_k.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
IF(TRIM(vname).EQ.'k22_shallow')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%cu_physics(id).EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%cu_physics(id).EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'kbcon_shallow')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%cu_physics(id).EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%cu_physics(id).EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'ktop_shallow')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%cu_physics(id).EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%cu_physics(id).EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'ktop_deep')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%cu_physics(id).EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%cu_physics(id).EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'kbcon_deep')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%cu_diag(id).EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%cu_diag(id).EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'ktop_deep')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%cu_diag(id).EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%cu_diag(id).EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'k22_deep')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%cu_diag(id).EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%cu_diag(id).EQ.1
  ENDIF
ENDIF
!ENDOFREGISTRYGENERATEDINCLUDE
