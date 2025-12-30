!STARTOFREGISTRYGENERATEDINCLUDE 'inc/in_use_for_config_r.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
IF(TRIM(vname).EQ.'rstoch')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sppt_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sppt_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'rand_real')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sppt_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sppt_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'rand_imag')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sppt_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sppt_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'rand_real_xxx')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sppt_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sppt_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'rand_real_yyy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sppt_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sppt_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'rand_imag_xxx')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sppt_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sppt_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'rand_imag_yyy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sppt_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sppt_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'rand_pert')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%rand_perturb_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%rand_perturb_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'rand_real')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%rand_perturb_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%rand_perturb_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'rand_imag')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%rand_perturb_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%rand_perturb_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'rand_real_xxx')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%rand_perturb_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%rand_perturb_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'rand_real_yyy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%rand_perturb_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%rand_perturb_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'rand_imag_xxx')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%rand_perturb_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%rand_perturb_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'rand_imag_yyy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%rand_perturb_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%rand_perturb_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'ru_tendf_stoch')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%skebs_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%skebs_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'rv_tendf_stoch')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%skebs_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%skebs_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'rt_tendf_stoch')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%skebs_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%skebs_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'rand_real')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%skebs_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%skebs_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'rand_imag')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%skebs_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%skebs_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'rand_real_xxx')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%skebs_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%skebs_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'rand_real_yyy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%skebs_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%skebs_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'rand_imag_xxx')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%skebs_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%skebs_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'rand_imag_yyy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%skebs_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%skebs_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'rand_real')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%spp_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%spp_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'rand_imag')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%spp_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%spp_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'rand_real_xxx')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%spp_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%spp_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'rand_real_yyy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%spp_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%spp_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'rand_imag_xxx')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%spp_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%spp_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'rand_imag_yyy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%spp_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%spp_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'rand_real')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%isppt_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%isppt_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'rand_imag')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%isppt_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%isppt_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'rand_real_xxx')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%isppt_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%isppt_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'rand_real_yyy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%isppt_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%isppt_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'rand_imag_xxx')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%isppt_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%isppt_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'rand_imag_yyy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%isppt_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%isppt_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'re_cloud')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics(id).EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics(id).EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'re_ice')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics(id).EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics(id).EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'re_snow')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics(id).EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics(id).EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'re_cloud')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'re_ice')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'re_snow')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'re_cloud')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics(id).EQ.6
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics(id).EQ.6
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'re_ice')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics(id).EQ.6
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics(id).EQ.6
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'re_snow')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics(id).EQ.6
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics(id).EQ.6
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'re_cloud')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics(id).EQ.8
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics(id).EQ.8
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'re_ice')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics(id).EQ.8
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics(id).EQ.8
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'re_snow')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics(id).EQ.8
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics(id).EQ.8
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'re_cloud')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics(id).EQ.28
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics(id).EQ.28
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'re_ice')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics(id).EQ.28
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics(id).EQ.28
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'re_snow')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics(id).EQ.28
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics(id).EQ.28
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'re_cloud')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics(id).EQ.14
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics(id).EQ.14
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'re_ice')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics(id).EQ.14
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics(id).EQ.14
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'re_snow')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics(id).EQ.14
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics(id).EQ.14
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'re_cloud')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics(id).EQ.16
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics(id).EQ.16
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'re_ice')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics(id).EQ.16
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics(id).EQ.16
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'re_snow')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics(id).EQ.16
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics(id).EQ.16
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'re_cloud')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics(id).EQ.17
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics(id).EQ.17
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'re_ice')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics(id).EQ.17
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics(id).EQ.17
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'re_snow')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics(id).EQ.17
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics(id).EQ.17
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'re_cloud')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics(id).EQ.18
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics(id).EQ.18
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'re_ice')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics(id).EQ.18
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics(id).EQ.18
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'re_snow')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics(id).EQ.18
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics(id).EQ.18
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'re_cloud')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics(id).EQ.22
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics(id).EQ.22
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'re_ice')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics(id).EQ.22
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics(id).EQ.22
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'re_snow')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%mp_physics(id).EQ.22
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%mp_physics(id).EQ.22
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'refl_10cm')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%compute_radar_ref.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%compute_radar_ref.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'refd_max')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%compute_radar_ref.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%compute_radar_ref.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'rtmassxy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'runsfxy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'runsbxy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'rssunxy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'rsshaxy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'rechxy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'rc2d')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%bl_pbl_physics(id).EQ.93
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%bl_pbl_physics(id).EQ.93
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'randstate1')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%cu_physics(id).EQ.95
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%cu_physics(id).EQ.95
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'randstate2')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%cu_physics(id).EQ.95
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%cu_physics(id).EQ.95
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'randstate3')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%cu_physics(id).EQ.95
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%cu_physics(id).EQ.95
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'randstate4')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%cu_physics(id).EQ.95
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%cu_physics(id).EQ.95
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'random')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%cu_physics(id).EQ.95
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%cu_physics(id).EQ.95
  ENDIF
ENDIF
!ENDOFREGISTRYGENERATEDINCLUDE
