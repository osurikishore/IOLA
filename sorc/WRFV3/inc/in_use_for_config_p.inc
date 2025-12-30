!STARTOFREGISTRYGENERATEDINCLUDE 'inc/in_use_for_config_p.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
IF(TRIM(vname).EQ.'pattern_spp_mp')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%spp_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%spp_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'pattern_isppt_mp')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%isppt_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%isppt_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'pattern_isppt_rad')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%isppt_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%isppt_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'pattern_isppt_conv')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%isppt_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%isppt_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'pattern_isppt_pbl')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%isppt_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%isppt_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'precip_swath')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%swath_mode.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%swath_mode.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'precip_swath')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%swath_mode.EQ.2
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%swath_mode.EQ.2
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'pdyn_parent')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.2
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.2
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'pdyn_smooth')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.2
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.2
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'pdyn_parent')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'pdyn_smooth')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'pdyn_parent')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'pdyn_smooth')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'pdyn_parent')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.5
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.5
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'pdyn_smooth')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.5
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.5
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'pdyn_parent')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.6
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.6
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'pdyn_smooth')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.6
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.6
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'p850rv')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.6
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.6
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'p700rv')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.6
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.6
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'p850wind')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.6
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.6
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'p700wind')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.6
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.6
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'p850z')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.6
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.6
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'p700z')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.6
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.6
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'pdyn_parent')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.7
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.7
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'pdyn_smooth')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.7
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.7
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'p850rv')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.7
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.7
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'p700rv')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.7
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.7
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'p850wind')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.7
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.7
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'p700wind')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.7
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.7
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'p850z')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.7
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.7
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'p700z')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.7
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.7
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'p500u')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.7
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.7
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'p500v')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.7
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.7
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'p700u')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.7
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.7
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'p700v')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.7
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.7
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'p850u')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.7
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.7
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'p850v')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.7
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.7
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'psnxy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'planting')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'pgsxy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
!ENDOFREGISTRYGENERATEDINCLUDE
