!STARTOFREGISTRYGENERATEDINCLUDE 'inc/in_use_for_config_m.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
IF(TRIM(vname).EQ.'mslp_noisy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'m10wind')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.6
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.6
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'m10rv')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.6
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.6
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'m10wind')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.7
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.7
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'m10rv')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.7
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.7
  ENDIF
ENDIF
!ENDOFREGISTRYGENERATEDINCLUDE
