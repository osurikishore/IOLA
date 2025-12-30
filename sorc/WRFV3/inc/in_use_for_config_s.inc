!STARTOFREGISTRYGENERATEDINCLUDE 'inc/in_use_for_config_s.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
IF(TRIM(vname).EQ.'savedtke12d')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_lake_physics(id).EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_lake_physics(id).EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'snowdp2d')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_lake_physics(id).EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_lake_physics(id).EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'snl2d')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_lake_physics(id).EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_lake_physics(id).EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'spptforcs')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sppt_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sppt_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'spptforcc')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sppt_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sppt_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'sppt_amp')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sppt_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sppt_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'spforcs')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%rand_perturb_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%rand_perturb_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'spforcc')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%rand_perturb_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%rand_perturb_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'sp_amp')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%rand_perturb_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%rand_perturb_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'spstreamforcc')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%skebs_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%skebs_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'spstreamforcs')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%skebs_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%skebs_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'sptforcc')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%skebs_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%skebs_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'sptforcs')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%skebs_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%skebs_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'spstream_amp')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%skebs_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%skebs_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'spt_amp')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%skebs_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%skebs_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'spforcs_spp_mp')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%spp_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%spp_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'spforcc_spp_mp')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%spp_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%spp_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'sp_amp_spp_mp')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%spp_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%spp_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'spforcs_isppt_mp')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%isppt_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%isppt_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'spforcc_isppt_mp')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%isppt_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%isppt_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'sp_amp_isppt_mp')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%isppt_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%isppt_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'spforcs_isppt_rad')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%isppt_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%isppt_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'spforcc_isppt_rad')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%isppt_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%isppt_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'sp_amp_isppt_rad')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%isppt_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%isppt_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'spforcs_isppt_conv')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%isppt_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%isppt_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'spforcc_isppt_conv')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%isppt_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%isppt_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'sp_amp_isppt_conv')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%isppt_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%isppt_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'spforcs_isppt_pbl')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%isppt_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%isppt_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'spforcc_isppt_pbl')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%isppt_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%isppt_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'sp_amp_isppt_pbl')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%isppt_on.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%isppt_on.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'suswind_swath')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%swath_mode.EQ.2
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%swath_mode.EQ.2
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'suswind')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%swath_mode.EQ.2
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%swath_mode.EQ.2
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'sp850rv')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.6
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.6
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'sp700rv')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.6
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.6
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'sp850wind')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.6
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.6
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'sp700wind')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.6
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.6
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'sp850z')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.6
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.6
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'sp700z')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.6
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.6
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'sm10wind')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.6
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.6
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'sm10rv')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.6
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.6
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'smslp')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.6
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.6
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'sp850rv')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.7
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.7
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'sp700rv')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.7
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.7
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'sp850z')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.7
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.7
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'sp700z')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.7
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.7
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'sm10rv')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.7
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.7
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'smslp')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%vortex_tracker(id).EQ.7
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%vortex_tracker(id).EQ.7
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'swupt')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%ra_lw_physics(id).EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%ra_lw_physics(id).EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'swuptc')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%ra_lw_physics(id).EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%ra_lw_physics(id).EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'swdnt')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%ra_lw_physics(id).EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%ra_lw_physics(id).EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'swdntc')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%ra_lw_physics(id).EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%ra_lw_physics(id).EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'swupb')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%ra_lw_physics(id).EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%ra_lw_physics(id).EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'swupbc')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%ra_lw_physics(id).EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%ra_lw_physics(id).EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'swdnb')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%ra_lw_physics(id).EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%ra_lw_physics(id).EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'swdnbc')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%ra_lw_physics(id).EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%ra_lw_physics(id).EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'swupt')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%ra_sw_physics(id).EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%ra_sw_physics(id).EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'swuptc')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%ra_sw_physics(id).EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%ra_sw_physics(id).EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'swdnt')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%ra_sw_physics(id).EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%ra_sw_physics(id).EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'swdntc')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%ra_sw_physics(id).EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%ra_sw_physics(id).EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'swupb')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%ra_sw_physics(id).EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%ra_sw_physics(id).EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'swupbc')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%ra_sw_physics(id).EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%ra_sw_physics(id).EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'swdnb')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%ra_sw_physics(id).EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%ra_sw_physics(id).EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'swdnbc')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%ra_sw_physics(id).EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%ra_sw_physics(id).EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'swupt')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%ra_sw_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%ra_sw_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'swuptc')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%ra_sw_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%ra_sw_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'swdnt')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%ra_sw_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%ra_sw_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'swdntc')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%ra_sw_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%ra_sw_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'swupb')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%ra_sw_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%ra_sw_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'swupbc')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%ra_sw_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%ra_sw_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'swdnb')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%ra_sw_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%ra_sw_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'swdnbc')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%ra_sw_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%ra_sw_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'sneqvoxy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'snicexy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'snliqxy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'stmassxy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'stblcpxy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'savxy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'sagxy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'shgxy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'shcxy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'shbxy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'smoiseq')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'smcwtdxy')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'season_gdd')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'soilcomp')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'soilcl1')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'soilcl2')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'soilcl3')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'soilcl4')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'sh3d')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%bl_pbl_physics(id).EQ.5
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%bl_pbl_physics(id).EQ.5
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'sh3d')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%bl_pbl_physics(id).EQ.6
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%bl_pbl_physics(id).EQ.6
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'snowsi')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%seaice_snowdepth_opt.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%seaice_snowdepth_opt.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'soldrain')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%wrf_hydro.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%wrf_hydro.EQ.1
  ENDIF
ENDIF
!ENDOFREGISTRYGENERATEDINCLUDE
