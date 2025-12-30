!STARTOFREGISTRYGENERATEDINCLUDE 'inc/in_use_for_config_v.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
IF(TRIM(vname).EQ.'vertstrucc')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sppt_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sppt_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'vertstrucs')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sppt_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sppt_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'vertampt')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sppt_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sppt_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'vertstrucc')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%rand_perturb_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%rand_perturb_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'vertstrucs')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%rand_perturb_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%rand_perturb_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'vertampt')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%rand_perturb_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%rand_perturb_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'vertstrucc')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%skebs_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%skebs_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'vertstrucs')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%skebs_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%skebs_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'vertampt')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%skebs_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%skebs_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'vertampuv')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%skebs_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%skebs_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'vertstrucc')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%spp_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%spp_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'vertstrucs')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%spp_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%spp_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'vertampt')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%spp_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%spp_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'vertstrucc')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%isppt_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%isppt_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'vertstrucs')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%isppt_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%isppt_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'vertampt')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%isppt_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%isppt_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'vdfg')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%grav_settling(id).EQ.0
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%grav_settling(id).EQ.0
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'vdfg')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%grav_settling(id).EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%grav_settling(id).EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'vdfg')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%grav_settling(id).EQ.2
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%grav_settling(id).EQ.2
  ENDIF
ENDIF
!ENDOFREGISTRYGENERATEDINCLUDE
