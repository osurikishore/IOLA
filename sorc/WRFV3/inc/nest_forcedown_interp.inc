!STARTOFREGISTRYGENERATEDINCLUDE 'inc/nest_forcedown_interp.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
IF ( SIZE( grid%rt_tendf_stoch, 1 ) * SIZE( grid%rt_tendf_stoch, 2 ) .GT. 1 ) THEN 
CALL force_sppt_nmm (  &         
                  grid%rt_tendf_stoch,   &       ! CD field
                 cids, cide, 1, config_flags%num_stoch_levels, cjds, cjde,   &         ! CD dims
                 cims, cime, 1, config_flags%num_stoch_levels, cjms, cjme,   &         ! CD dims
                 cips, cipe, 1, config_flags%num_stoch_levels, cjps, cjpe,   &         ! CD dims
                  ngrid%rt_tendf_stoch,  &   ! ND field
                 nids, nide, 1, config_flags%num_stoch_levels, njds, njde,   &         ! ND dims
                 nims, nime, 1, config_flags%num_stoch_levels, njms, njme,   &         ! ND dims
                 nips, nipe, 1, config_flags%num_stoch_levels, njps, njpe,   &         ! ND dims
                  config_flags%shw, ngrid%imask_nostag,         &         ! stencil half width
                  .FALSE., .FALSE.,                                                &         ! xstag, ystag
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%iih,ngrid%iih&
,grid%jjh,ngrid%jjh&
,grid%hbwgt1,ngrid%hbwgt1&
,grid%hbwgt2,ngrid%hbwgt2&
,grid%hbwgt3,ngrid%hbwgt3&
,grid%hbwgt4,ngrid%hbwgt4&
                  ) 
ENDIF
IF ( SIZE( grid%rand_pert, 1 ) * SIZE( grid%rand_pert, 2 ) .GT. 1 ) THEN 
CALL force_sppt_nmm (  &         
                  grid%rand_pert,   &       ! CD field
                 cids, cide, 1, config_flags%num_stoch_levels, cjds, cjde,   &         ! CD dims
                 cims, cime, 1, config_flags%num_stoch_levels, cjms, cjme,   &         ! CD dims
                 cips, cipe, 1, config_flags%num_stoch_levels, cjps, cjpe,   &         ! CD dims
                  ngrid%rand_pert,  &   ! ND field
                 nids, nide, 1, config_flags%num_stoch_levels, njds, njde,   &         ! ND dims
                 nims, nime, 1, config_flags%num_stoch_levels, njms, njme,   &         ! ND dims
                 nips, nipe, 1, config_flags%num_stoch_levels, njps, njpe,   &         ! ND dims
                  config_flags%shw, ngrid%imask_nostag,         &         ! stencil half width
                  .FALSE., .FALSE.,                                                &         ! xstag, ystag
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%iih,ngrid%iih&
,grid%jjh,ngrid%jjh&
,grid%hbwgt1,ngrid%hbwgt1&
,grid%hbwgt2,ngrid%hbwgt2&
,grid%hbwgt3,ngrid%hbwgt3&
,grid%hbwgt4,ngrid%hbwgt4&
                  ) 
ENDIF
IF ( SIZE( grid%rstoch, 1 ) * SIZE( grid%rstoch, 2 ) .GT. 1 ) THEN 
CALL force_sppt_nmm (  &         
                  grid%rstoch,   &       ! CD field
                 cids, cide, 1, config_flags%num_stoch_levels, cjds, cjde,   &         ! CD dims
                 cims, cime, 1, config_flags%num_stoch_levels, cjms, cjme,   &         ! CD dims
                 cips, cipe, 1, config_flags%num_stoch_levels, cjps, cjpe,   &         ! CD dims
                  ngrid%rstoch,  &   ! ND field
                 nids, nide, 1, config_flags%num_stoch_levels, njds, njde,   &         ! ND dims
                 nims, nime, 1, config_flags%num_stoch_levels, njms, njme,   &         ! ND dims
                 nips, nipe, 1, config_flags%num_stoch_levels, njps, njpe,   &         ! ND dims
                  config_flags%shw, ngrid%imask_nostag,         &         ! stencil half width
                  .FALSE., .FALSE.,                                                &         ! xstag, ystag
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%iih,ngrid%iih&
,grid%jjh,ngrid%jjh&
,grid%hbwgt1,ngrid%hbwgt1&
,grid%hbwgt2,ngrid%hbwgt2&
,grid%hbwgt3,ngrid%hbwgt3&
,grid%hbwgt4,ngrid%hbwgt4&
                  ) 
ENDIF
IF ( SIZE( grid%pattern_spp_mp, 1 ) * SIZE( grid%pattern_spp_mp, 2 ) .GT. 1 ) THEN 
CALL force_sppt_nmm (  &         
                  grid%pattern_spp_mp,   &       ! CD field
                 cids, cide, 1, config_flags%num_stoch_levels, cjds, cjde,   &         ! CD dims
                 cims, cime, 1, config_flags%num_stoch_levels, cjms, cjme,   &         ! CD dims
                 cips, cipe, 1, config_flags%num_stoch_levels, cjps, cjpe,   &         ! CD dims
                  ngrid%pattern_spp_mp,  &   ! ND field
                 nids, nide, 1, config_flags%num_stoch_levels, njds, njde,   &         ! ND dims
                 nims, nime, 1, config_flags%num_stoch_levels, njms, njme,   &         ! ND dims
                 nips, nipe, 1, config_flags%num_stoch_levels, njps, njpe,   &         ! ND dims
                  config_flags%shw, ngrid%imask_nostag,         &         ! stencil half width
                  .FALSE., .FALSE.,                                                &         ! xstag, ystag
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%iih,ngrid%iih&
,grid%jjh,ngrid%jjh&
,grid%hbwgt1,ngrid%hbwgt1&
,grid%hbwgt2,ngrid%hbwgt2&
,grid%hbwgt3,ngrid%hbwgt3&
,grid%hbwgt4,ngrid%hbwgt4&
                  ) 
ENDIF
IF ( SIZE( grid%pattern_isppt_mp, 1 ) * SIZE( grid%pattern_isppt_mp, 2 ) .GT. 1 ) THEN 
CALL force_sppt_nmm (  &         
                  grid%pattern_isppt_mp,   &       ! CD field
                 cids, cide, 1, config_flags%num_stoch_levels, cjds, cjde,   &         ! CD dims
                 cims, cime, 1, config_flags%num_stoch_levels, cjms, cjme,   &         ! CD dims
                 cips, cipe, 1, config_flags%num_stoch_levels, cjps, cjpe,   &         ! CD dims
                  ngrid%pattern_isppt_mp,  &   ! ND field
                 nids, nide, 1, config_flags%num_stoch_levels, njds, njde,   &         ! ND dims
                 nims, nime, 1, config_flags%num_stoch_levels, njms, njme,   &         ! ND dims
                 nips, nipe, 1, config_flags%num_stoch_levels, njps, njpe,   &         ! ND dims
                  config_flags%shw, ngrid%imask_nostag,         &         ! stencil half width
                  .FALSE., .FALSE.,                                                &         ! xstag, ystag
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%iih,ngrid%iih&
,grid%jjh,ngrid%jjh&
,grid%hbwgt1,ngrid%hbwgt1&
,grid%hbwgt2,ngrid%hbwgt2&
,grid%hbwgt3,ngrid%hbwgt3&
,grid%hbwgt4,ngrid%hbwgt4&
                  ) 
ENDIF
IF ( SIZE( grid%pattern_isppt_pbl, 1 ) * SIZE( grid%pattern_isppt_pbl, 2 ) .GT. 1 ) THEN 
CALL force_sppt_nmm (  &         
                  grid%pattern_isppt_pbl,   &       ! CD field
                 cids, cide, 1, config_flags%num_stoch_levels, cjds, cjde,   &         ! CD dims
                 cims, cime, 1, config_flags%num_stoch_levels, cjms, cjme,   &         ! CD dims
                 cips, cipe, 1, config_flags%num_stoch_levels, cjps, cjpe,   &         ! CD dims
                  ngrid%pattern_isppt_pbl,  &   ! ND field
                 nids, nide, 1, config_flags%num_stoch_levels, njds, njde,   &         ! ND dims
                 nims, nime, 1, config_flags%num_stoch_levels, njms, njme,   &         ! ND dims
                 nips, nipe, 1, config_flags%num_stoch_levels, njps, njpe,   &         ! ND dims
                  config_flags%shw, ngrid%imask_nostag,         &         ! stencil half width
                  .FALSE., .FALSE.,                                                &         ! xstag, ystag
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%iih,ngrid%iih&
,grid%jjh,ngrid%jjh&
,grid%hbwgt1,ngrid%hbwgt1&
,grid%hbwgt2,ngrid%hbwgt2&
,grid%hbwgt3,ngrid%hbwgt3&
,grid%hbwgt4,ngrid%hbwgt4&
                  ) 
ENDIF
IF ( SIZE( grid%pattern_isppt_rad, 1 ) * SIZE( grid%pattern_isppt_rad, 2 ) .GT. 1 ) THEN 
CALL force_sppt_nmm (  &         
                  grid%pattern_isppt_rad,   &       ! CD field
                 cids, cide, 1, config_flags%num_stoch_levels, cjds, cjde,   &         ! CD dims
                 cims, cime, 1, config_flags%num_stoch_levels, cjms, cjme,   &         ! CD dims
                 cips, cipe, 1, config_flags%num_stoch_levels, cjps, cjpe,   &         ! CD dims
                  ngrid%pattern_isppt_rad,  &   ! ND field
                 nids, nide, 1, config_flags%num_stoch_levels, njds, njde,   &         ! ND dims
                 nims, nime, 1, config_flags%num_stoch_levels, njms, njme,   &         ! ND dims
                 nips, nipe, 1, config_flags%num_stoch_levels, njps, njpe,   &         ! ND dims
                  config_flags%shw, ngrid%imask_nostag,         &         ! stencil half width
                  .FALSE., .FALSE.,                                                &         ! xstag, ystag
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%iih,ngrid%iih&
,grid%jjh,ngrid%jjh&
,grid%hbwgt1,ngrid%hbwgt1&
,grid%hbwgt2,ngrid%hbwgt2&
,grid%hbwgt3,ngrid%hbwgt3&
,grid%hbwgt4,ngrid%hbwgt4&
                  ) 
ENDIF
IF ( SIZE( grid%pattern_isppt_conv, 1 ) * SIZE( grid%pattern_isppt_conv, 2 ) .GT. 1 ) THEN 
CALL force_sppt_nmm (  &         
                  grid%pattern_isppt_conv,   &       ! CD field
                 cids, cide, 1, config_flags%num_stoch_levels, cjds, cjde,   &         ! CD dims
                 cims, cime, 1, config_flags%num_stoch_levels, cjms, cjme,   &         ! CD dims
                 cips, cipe, 1, config_flags%num_stoch_levels, cjps, cjpe,   &         ! CD dims
                  ngrid%pattern_isppt_conv,  &   ! ND field
                 nids, nide, 1, config_flags%num_stoch_levels, njds, njde,   &         ! ND dims
                 nims, nime, 1, config_flags%num_stoch_levels, njms, njme,   &         ! ND dims
                 nips, nipe, 1, config_flags%num_stoch_levels, njps, njpe,   &         ! ND dims
                  config_flags%shw, ngrid%imask_nostag,         &         ! stencil half width
                  .FALSE., .FALSE.,                                                &         ! xstag, ystag
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%iih,ngrid%iih&
,grid%jjh,ngrid%jjh&
,grid%hbwgt1,ngrid%hbwgt1&
,grid%hbwgt2,ngrid%hbwgt2&
,grid%hbwgt3,ngrid%hbwgt3&
,grid%hbwgt4,ngrid%hbwgt4&
                  ) 
ENDIF
IF ( SIZE( grid%lu_index, 1 ) * SIZE( grid%lu_index, 2 ) .GT. 1 ) THEN 
CONTINUE ! do not call nointerp
ENDIF
IF ( SIZE( grid%lu_mask, 1 ) * SIZE( grid%lu_mask, 2 ) .GT. 1 ) THEN 
CONTINUE ! do not call nointerp
ENDIF
IF ( SIZE( grid%pdyn_smooth, 1 ) * SIZE( grid%pdyn_smooth, 2 ) .GT. 1 ) THEN 
CONTINUE ! do not call nointerp
ENDIF
IF ( SIZE( grid%pdyn_parent, 1 ) * SIZE( grid%pdyn_parent, 2 ) .GT. 1 ) THEN 
CALL downaged2d (  &         
                  grid%pdyn_parent,   &       ! CD field
                 cids, cide, 1, 1, cjds, cjde,   &         ! CD dims
                 cims, cime, 1, 1, cjms, cjme,   &         ! CD dims
                 cips, cipe, 1, 1, cjps, cjpe,   &         ! CD dims
                  ngrid%pdyn_parent,  &   ! ND field
                 nids, nide, 1, 1, njds, njde,   &         ! ND dims
                 nims, nime, 1, 1, njms, njme,   &         ! ND dims
                 nips, nipe, 1, 1, njps, njpe,   &         ! ND dims
                  config_flags%shw, ngrid%imask_nostag,         &         ! stencil half width
                  .FALSE., .FALSE.,                                                &         ! xstag, ystag
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%pdyn_smooth_age&
,ngrid%pdyn_parent_age&
,grid%pdyn_smooth&
                  ) 
ENDIF
IF ( SIZE( grid%hres_fis, 1 ) * SIZE( grid%hres_fis, 2 ) .GT. 1 ) THEN 
CONTINUE ! do not call nointerp
ENDIF
IF ( SIZE( grid%hres_lnd, 1 ) * SIZE( grid%hres_lnd, 2 ) .GT. 1 ) THEN 
CONTINUE ! do not call nointerp
ENDIF
IF ( SIZE( grid%hres_lu_index, 1 ) * SIZE( grid%hres_lu_index, 2 ) .GT. 1 ) THEN 
CONTINUE ! do not call nointerp
ENDIF
IF ( SIZE( grid%sm, 1 ) * SIZE( grid%sm, 2 ) .GT. 1 ) THEN 
CONTINUE ! do not call nointerp
ENDIF
IF ( SIZE( grid%pd, 1 ) * SIZE( grid%pd, 2 ) .GT. 1 ) THEN 
CONTINUE ! do not call nointerp
ENDIF
IF ( SIZE( grid%fis, 1 ) * SIZE( grid%fis, 2 ) .GT. 1 ) THEN 
CONTINUE ! do not call nointerp
ENDIF
IF ( SIZE( grid%t, 1 ) * SIZE( grid%t, 2 ) .GT. 1 ) THEN 
CONTINUE ! do not call nointerp
ENDIF
IF ( SIZE( grid%q, 1 ) * SIZE( grid%q, 2 ) .GT. 1 ) THEN 
CONTINUE ! do not call nointerp
ENDIF
IF ( SIZE( grid%u, 1 ) * SIZE( grid%u, 2 ) .GT. 1 ) THEN 
CALL bdyvel (  &         
                  grid%u,   &       ! CD field
                 cids, cide, ckds, ckde, cjds, cjde,   &         ! CD dims
                 cims, cime, ckms, ckme, cjms, cjme,   &         ! CD dims
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         ! CD dims
                  ngrid%u,  &   ! ND field
                 nids, nide, nkds, nkde, njds, njde,   &         ! ND dims
                 nims, nime, nkms, nkme, njms, njme,   &         ! ND dims
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         ! ND dims
                  config_flags%shw, ngrid%imask_nostag,         &         ! stencil half width
                  .FALSE., .FALSE.,                                                &         ! xstag, ystag
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,dummy_xs ,ngrid%u_bxs &
,dummy_xe ,ngrid%u_bxe &
,dummy_ys ,ngrid%u_bys &
,dummy_ye ,ngrid%u_bye &
,dummy_xs ,ngrid%u_btxs &
,dummy_xe ,ngrid%u_btxe &
,dummy_ys ,ngrid%u_btys &
,dummy_ye ,ngrid%u_btye &
                  ) 
ENDIF
IF ( SIZE( grid%v, 1 ) * SIZE( grid%v, 2 ) .GT. 1 ) THEN 
CALL bdyvel (  &         
                  grid%v,   &       ! CD field
                 cids, cide, ckds, ckde, cjds, cjde,   &         ! CD dims
                 cims, cime, ckms, ckme, cjms, cjme,   &         ! CD dims
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         ! CD dims
                  ngrid%v,  &   ! ND field
                 nids, nide, nkds, nkde, njds, njde,   &         ! ND dims
                 nims, nime, nkms, nkme, njms, njme,   &         ! ND dims
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         ! ND dims
                  config_flags%shw, ngrid%imask_nostag,         &         ! stencil half width
                  .FALSE., .FALSE.,                                                &         ! xstag, ystag
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,dummy_xs ,ngrid%v_bxs &
,dummy_xe ,ngrid%v_bxe &
,dummy_ys ,ngrid%v_bys &
,dummy_ye ,ngrid%v_bye &
,dummy_xs ,ngrid%v_btxs &
,dummy_xe ,ngrid%v_btxe &
,dummy_ys ,ngrid%v_btys &
,dummy_ye ,ngrid%v_btye &
                  ) 
ENDIF
IF ( SIZE( grid%scurx, 1 ) * SIZE( grid%scurx, 2 ) .GT. 1 ) THEN 
CALL force_sst_nmm (  &         
                  grid%scurx,   &       ! CD field
                 cids, cide, 1, 1, cjds, cjde,   &         ! CD dims
                 cims, cime, 1, 1, cjms, cjme,   &         ! CD dims
                 cips, cipe, 1, 1, cjps, cjpe,   &         ! CD dims
                  ngrid%scurx,  &   ! ND field
                 nids, nide, 1, 1, njds, njde,   &         ! ND dims
                 nims, nime, 1, 1, njms, njme,   &         ! ND dims
                 nips, nipe, 1, 1, njps, njpe,   &         ! ND dims
                  config_flags%shw, ngrid%imask_nostag,         &         ! stencil half width
                  .FALSE., .FALSE.,                                                &         ! xstag, ystag
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%iih,ngrid%iih&
,grid%jjh,ngrid%jjh&
,grid%hbwgt1,ngrid%hbwgt1&
,grid%hbwgt2,ngrid%hbwgt2&
,grid%hbwgt3,ngrid%hbwgt3&
,grid%hbwgt4,ngrid%hbwgt4&
,grid%force_sst,ngrid%force_sst&
                  ) 
ENDIF
IF ( SIZE( grid%scury, 1 ) * SIZE( grid%scury, 2 ) .GT. 1 ) THEN 
CALL force_sst_nmm (  &         
                  grid%scury,   &       ! CD field
                 cids, cide, 1, 1, cjds, cjde,   &         ! CD dims
                 cims, cime, 1, 1, cjms, cjme,   &         ! CD dims
                 cips, cipe, 1, 1, cjps, cjpe,   &         ! CD dims
                  ngrid%scury,  &   ! ND field
                 nids, nide, 1, 1, njds, njde,   &         ! ND dims
                 nims, nime, 1, 1, njms, njme,   &         ! ND dims
                 nips, nipe, 1, 1, njps, njpe,   &         ! ND dims
                  config_flags%shw, ngrid%imask_nostag,         &         ! stencil half width
                  .FALSE., .FALSE.,                                                &         ! xstag, ystag
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%iih,ngrid%iih&
,grid%jjh,ngrid%jjh&
,grid%hbwgt1,ngrid%hbwgt1&
,grid%hbwgt2,ngrid%hbwgt2&
,grid%hbwgt3,ngrid%hbwgt3&
,grid%hbwgt4,ngrid%hbwgt4&
,grid%force_sst,ngrid%force_sst&
                  ) 
ENDIF
IF ( SIZE( grid%charn, 1 ) * SIZE( grid%charn, 2 ) .GT. 1 ) THEN 
CALL force_sst_nmm (  &         
                  grid%charn,   &       ! CD field
                 cids, cide, 1, 1, cjds, cjde,   &         ! CD dims
                 cims, cime, 1, 1, cjms, cjme,   &         ! CD dims
                 cips, cipe, 1, 1, cjps, cjpe,   &         ! CD dims
                  ngrid%charn,  &   ! ND field
                 nids, nide, 1, 1, njds, njde,   &         ! ND dims
                 nims, nime, 1, 1, njms, njme,   &         ! ND dims
                 nips, nipe, 1, 1, njps, njpe,   &         ! ND dims
                  config_flags%shw, ngrid%imask_nostag,         &         ! stencil half width
                  .FALSE., .FALSE.,                                                &         ! xstag, ystag
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%iih,ngrid%iih&
,grid%jjh,ngrid%jjh&
,grid%hbwgt1,ngrid%hbwgt1&
,grid%hbwgt2,ngrid%hbwgt2&
,grid%hbwgt3,ngrid%hbwgt3&
,grid%hbwgt4,ngrid%hbwgt4&
,grid%force_sst,ngrid%force_sst&
                  ) 
ENDIF
IF ( SIZE( grid%msang, 1 ) * SIZE( grid%msang, 2 ) .GT. 1 ) THEN 
CALL force_sst_nmm (  &         
                  grid%msang,   &       ! CD field
                 cids, cide, 1, 1, cjds, cjde,   &         ! CD dims
                 cims, cime, 1, 1, cjms, cjme,   &         ! CD dims
                 cips, cipe, 1, 1, cjps, cjpe,   &         ! CD dims
                  ngrid%msang,  &   ! ND field
                 nids, nide, 1, 1, njds, njde,   &         ! ND dims
                 nims, nime, 1, 1, njms, njme,   &         ! ND dims
                 nips, nipe, 1, 1, njps, njpe,   &         ! ND dims
                  config_flags%shw, ngrid%imask_nostag,         &         ! stencil half width
                  .FALSE., .FALSE.,                                                &         ! xstag, ystag
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%iih,ngrid%iih&
,grid%jjh,ngrid%jjh&
,grid%hbwgt1,ngrid%hbwgt1&
,grid%hbwgt2,ngrid%hbwgt2&
,grid%hbwgt3,ngrid%hbwgt3&
,grid%hbwgt4,ngrid%hbwgt4&
,grid%force_sst,ngrid%force_sst&
                  ) 
ENDIF
IF ( SIZE( grid%q2, 1 ) * SIZE( grid%q2, 2 ) .GT. 1 ) THEN 
CALL bdymass (  &         
                  grid%q2,   &       ! CD field
                 cids, cide, ckds, ckde, cjds, cjde,   &         ! CD dims
                 cims, cime, ckms, ckme, cjms, cjme,   &         ! CD dims
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         ! CD dims
                  ngrid%q2,  &   ! ND field
                 nids, nide, nkds, nkde, njds, njde,   &         ! ND dims
                 nims, nime, nkms, nkme, njms, njme,   &         ! ND dims
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         ! ND dims
                  config_flags%shw, ngrid%imask_nostag,         &         ! stencil half width
                  .FALSE., .FALSE.,                                                &         ! xstag, ystag
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,dummy_xs ,ngrid%q2_bxs &
,dummy_xe ,ngrid%q2_bxe &
,dummy_ys ,ngrid%q2_bys &
,dummy_ye ,ngrid%q2_bye &
,dummy_xs ,ngrid%q2_btxs &
,dummy_xe ,ngrid%q2_btxe &
,dummy_ys ,ngrid%q2_btys &
,dummy_ye ,ngrid%q2_btye &
,econst &
,0.000000000e+00 &
                  ) 
ENDIF
IF ( SIZE( grid%cwm, 1 ) * SIZE( grid%cwm, 2 ) .GT. 1 ) THEN 
CALL bdymass (  &         
                  grid%cwm,   &       ! CD field
                 cids, cide, ckds, ckde, cjds, cjde,   &         ! CD dims
                 cims, cime, ckms, ckme, cjms, cjme,   &         ! CD dims
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         ! CD dims
                  ngrid%cwm,  &   ! ND field
                 nids, nide, nkds, nkde, njds, njde,   &         ! ND dims
                 nims, nime, nkms, nkme, njms, njme,   &         ! ND dims
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         ! ND dims
                  config_flags%shw, ngrid%imask_nostag,         &         ! stencil half width
                  .FALSE., .FALSE.,                                                &         ! xstag, ystag
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,dummy_xs ,ngrid%cwm_bxs &
,dummy_xe ,ngrid%cwm_bxe &
,dummy_ys ,ngrid%cwm_bys &
,dummy_ye ,ngrid%cwm_bye &
,dummy_xs ,ngrid%cwm_btxs &
,dummy_xe ,ngrid%cwm_btxe &
,dummy_ys ,ngrid%cwm_btys &
,dummy_ye ,ngrid%cwm_btye &
,ecopy &
,0.000000000e+00 &
                  ) 
ENDIF
IF ( SIZE( grid%pint, 1 ) * SIZE( grid%pint, 2 ) .GT. 1 ) THEN 
CONTINUE ! do not call nointerp
ENDIF
DO itrace = PARAM_FIRST_SCALAR, num_moist
IF ( SIZE( moist, 1 ) * SIZE( moist, 2 ) .GT. 1 .and. (interp_mp .eqv. .true.) ) THEN 
CALL bdymass (  &         
                  moist(grid%sm31,grid%sm32,grid%sm33,itrace),   &       ! CD field
                 cids, cide, ckds, ckde, cjds, cjde,   &         ! CD dims
                 cims, cime, ckms, ckme, cjms, cjme,   &         ! CD dims
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         ! CD dims
                  ngrid%moist(ngrid%sm31,ngrid%sm32,ngrid%sm33,itrace),  &   ! ND field
                 nids, nide, nkds, nkde, njds, njde,   &         ! ND dims
                 nims, nime, nkms, nkme, njms, njme,   &         ! ND dims
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         ! ND dims
                  config_flags%shw, ngrid%imask_nostag,         &         ! stencil half width
                  .FALSE., .FALSE.,                                                &         ! xstag, ystag
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,moist_bxs(cjms,1,1,itrace) ,ngrid%moist_bxs(njms,1,1,itrace) &
,moist_bxe(cjms,1,1,itrace) ,ngrid%moist_bxe(njms,1,1,itrace) &
,moist_bys(cims,1,1,itrace) ,ngrid%moist_bys(nims,1,1,itrace) &
,moist_bye(cims,1,1,itrace) ,ngrid%moist_bye(nims,1,1,itrace) &
,moist_btxs(cjms,1,1,itrace) ,ngrid%moist_btxs(njms,1,1,itrace) &
,moist_btxe(cjms,1,1,itrace) ,ngrid%moist_btxe(njms,1,1,itrace) &
,moist_btys(cims,1,1,itrace) ,ngrid%moist_btys(nims,1,1,itrace) &
,moist_btye(cims,1,1,itrace) ,ngrid%moist_btye(nims,1,1,itrace) &
,ecopy &
,0.000000000e+00 &
                  ) 
ENDIF
ENDDO
DO itrace = PARAM_FIRST_SCALAR, num_scalar
IF ( SIZE( scalar, 1 ) * SIZE( scalar, 2 ) .GT. 1 .and. (interp_mp .eqv. .true.) ) THEN 
CALL bdymass (  &         
                  scalar(grid%sm31,grid%sm32,grid%sm33,itrace),   &       ! CD field
                 cids, cide, ckds, ckde, cjds, cjde,   &         ! CD dims
                 cims, cime, ckms, ckme, cjms, cjme,   &         ! CD dims
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         ! CD dims
                  ngrid%scalar(ngrid%sm31,ngrid%sm32,ngrid%sm33,itrace),  &   ! ND field
                 nids, nide, nkds, nkde, njds, njde,   &         ! ND dims
                 nims, nime, nkms, nkme, njms, njme,   &         ! ND dims
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         ! ND dims
                  config_flags%shw, ngrid%imask_nostag,         &         ! stencil half width
                  .FALSE., .FALSE.,                                                &         ! xstag, ystag
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,scalar_bxs(cjms,1,1,itrace) ,ngrid%scalar_bxs(njms,1,1,itrace) &
,scalar_bxe(cjms,1,1,itrace) ,ngrid%scalar_bxe(njms,1,1,itrace) &
,scalar_bys(cims,1,1,itrace) ,ngrid%scalar_bys(nims,1,1,itrace) &
,scalar_bye(cims,1,1,itrace) ,ngrid%scalar_bye(nims,1,1,itrace) &
,scalar_btxs(cjms,1,1,itrace) ,ngrid%scalar_btxs(njms,1,1,itrace) &
,scalar_btxe(cjms,1,1,itrace) ,ngrid%scalar_btxe(njms,1,1,itrace) &
,scalar_btys(cims,1,1,itrace) ,ngrid%scalar_btys(nims,1,1,itrace) &
,scalar_btye(cims,1,1,itrace) ,ngrid%scalar_btye(nims,1,1,itrace) &
,ecopy &
,0.000000000e+00 &
                  ) 
ENDIF
ENDDO
DO itrace = PARAM_FIRST_SCALAR, num_dfi_scalar
IF ( SIZE( dfi_scalar, 1 ) * SIZE( dfi_scalar, 2 ) .GT. 1 .and. (interp_mp .eqv. .true.) ) THEN 
CALL bdymass (  &         
                  dfi_scalar(grid%sm31,grid%sm32,grid%sm33,itrace),   &       ! CD field
                 cids, cide, ckds, ckde, cjds, cjde,   &         ! CD dims
                 cims, cime, ckms, ckme, cjms, cjme,   &         ! CD dims
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         ! CD dims
                  ngrid%dfi_scalar(ngrid%sm31,ngrid%sm32,ngrid%sm33,itrace),  &   ! ND field
                 nids, nide, nkds, nkde, njds, njde,   &         ! ND dims
                 nims, nime, nkms, nkme, njms, njme,   &         ! ND dims
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         ! ND dims
                  config_flags%shw, ngrid%imask_nostag,         &         ! stencil half width
                  .FALSE., .FALSE.,                                                &         ! xstag, ystag
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,dfi_scalar_bxs(cjms,1,1,itrace) ,ngrid%dfi_scalar_bxs(njms,1,1,itrace) &
,dfi_scalar_bxe(cjms,1,1,itrace) ,ngrid%dfi_scalar_bxe(njms,1,1,itrace) &
,dfi_scalar_bys(cims,1,1,itrace) ,ngrid%dfi_scalar_bys(nims,1,1,itrace) &
,dfi_scalar_bye(cims,1,1,itrace) ,ngrid%dfi_scalar_bye(nims,1,1,itrace) &
,dfi_scalar_btxs(cjms,1,1,itrace) ,ngrid%dfi_scalar_btxs(njms,1,1,itrace) &
,dfi_scalar_btxe(cjms,1,1,itrace) ,ngrid%dfi_scalar_btxe(njms,1,1,itrace) &
,dfi_scalar_btys(cims,1,1,itrace) ,ngrid%dfi_scalar_btys(nims,1,1,itrace) &
,dfi_scalar_btye(cims,1,1,itrace) ,ngrid%dfi_scalar_btye(nims,1,1,itrace) &
,ecopy &
,0.000000000e+00 &
                  ) 
ENDIF
ENDDO
IF ( SIZE( grid%sst, 1 ) * SIZE( grid%sst, 2 ) .GT. 1 ) THEN 
CALL force_sst_nmm (  &         
                  grid%sst,   &       ! CD field
                 cids, cide, 1, 1, cjds, cjde,   &         ! CD dims
                 cims, cime, 1, 1, cjms, cjme,   &         ! CD dims
                 cips, cipe, 1, 1, cjps, cjpe,   &         ! CD dims
                  ngrid%sst,  &   ! ND field
                 nids, nide, 1, 1, njds, njde,   &         ! ND dims
                 nims, nime, 1, 1, njms, njme,   &         ! ND dims
                 nips, nipe, 1, 1, njps, njpe,   &         ! ND dims
                  config_flags%shw, ngrid%imask_nostag,         &         ! stencil half width
                  .FALSE., .FALSE.,                                                &         ! xstag, ystag
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%iih,ngrid%iih&
,grid%jjh,ngrid%jjh&
,grid%hbwgt1,ngrid%hbwgt1&
,grid%hbwgt2,ngrid%hbwgt2&
,grid%hbwgt3,ngrid%hbwgt3&
,grid%hbwgt4,ngrid%hbwgt4&
,grid%force_sst,ngrid%force_sst&
                  ) 
ENDIF
!ENDOFREGISTRYGENERATEDINCLUDE
