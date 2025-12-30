!STARTOFREGISTRYGENERATEDINCLUDE 'inc/in_use_for_config_f.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
IF(TRIM(vname).EQ.'field_spp_mp')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%spp_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%spp_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'f_ice')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics(id).EQ.5
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics(id).EQ.5
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'f_rain')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics(id).EQ.5
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics(id).EQ.5
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'f_rimef')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics(id).EQ.5
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics(id).EQ.5
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'f_ice')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics(id).EQ.85
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics(id).EQ.85
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'f_rain')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics(id).EQ.85
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics(id).EQ.85
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'f_rimef')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics(id).EQ.85
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics(id).EQ.85
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'f_ice')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics(id).EQ.95
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics(id).EQ.95
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'f_rain')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics(id).EQ.95
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics(id).EQ.95
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'f_rimef')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics(id).EQ.95
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics(id).EQ.95
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'flx4')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.2
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.2
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'fvb')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.2
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.2
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'fbur')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.2
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.2
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'fgsn')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.2
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.2
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'fwetxy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'fastcpxy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'fvegxy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'fsaxy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'firaxy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'fgdp')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%grav_settling(id).EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%grav_settling(id).EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'fgdp')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%grav_settling(id).EQ.2
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%grav_settling(id).EQ.2
  ENDIF
ENDIF
!ENDOFREGISTRYGENERATEDINCLUDE
