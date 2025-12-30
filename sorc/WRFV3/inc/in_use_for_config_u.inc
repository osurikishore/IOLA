!STARTOFREGISTRYGENERATEDINCLUDE 'inc/in_use_for_config_u.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
IF(TRIM(vname).EQ.'udr_kf')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%kf_edrates(id).EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%kf_edrates(id).EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'uer_kf')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%kf_edrates(id).EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%kf_edrates(id).EQ.1
  ENDIF
ENDIF
!ENDOFREGISTRYGENERATEDINCLUDE
