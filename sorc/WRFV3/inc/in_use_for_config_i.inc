!STARTOFREGISTRYGENERATEDINCLUDE 'inc/in_use_for_config_i.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
IF(TRIM(vname).EQ.'interesting')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%swath_mode.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%swath_mode.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'interesting')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%swath_mode.EQ.2
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%swath_mode.EQ.2
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'isnowxy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'irgxy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'ircxy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'irbxy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'icedepth')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%seaice_thickness_opt.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%seaice_thickness_opt.EQ.1
  ENDIF
ENDIF
!ENDOFREGISTRYGENERATEDINCLUDE
