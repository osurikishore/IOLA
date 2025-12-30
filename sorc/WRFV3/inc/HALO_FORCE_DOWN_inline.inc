!STARTOFREGISTRYGENERATEDINCLUDE 'inc/HALO_FORCE_DOWN_inline.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL wrf_debug(2,'calling inc/HALO_FORCE_DOWN_inline.inc')
CALL rsl_comm_iter_init(4,jps,jpe)
DO WHILE ( rsl_comm_iter( grid%id , grid%is_intermediate, 4 , &
                         0 , jds,jde,jps,jpe, grid%njds, grid%njde , & 
     rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p,   & 
     rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p    ))
 CALL RSL_LITE_INIT_EXCH ( local_communicator, 4, 0, &
     rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p,   & 
     rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p,   & 
     15  &
   + num_moist   &
   + num_scalar   &
   + num_dfi_scalar   &
     , 15, RWORDSIZE, &
     0, 0, IWORDSIZE, &
     0, 0, DWORDSIZE, &
      0,  0, LWORDSIZE, &
      mytask, ntasks, ntasks_x, ntasks_y,   &
      ips, ipe, jps, jpe, kps, MAX(1,1&
,config_flags%num_stoch_levels &
,kpe &
))
IF ( SIZE(grid%rt_tendf_stoch,1)*SIZE(grid%rt_tendf_stoch,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%rt_tendf_stoch, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1, config_flags%num_stoch_levels,             &
ims, ime, jms, jme, 1, config_flags%num_stoch_levels,             &
ips, ipe, jps, jpe, 1, config_flags%num_stoch_levels              )
ENDIF
IF ( SIZE(grid%rand_pert,1)*SIZE(grid%rand_pert,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%rand_pert, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1, config_flags%num_stoch_levels,             &
ims, ime, jms, jme, 1, config_flags%num_stoch_levels,             &
ips, ipe, jps, jpe, 1, config_flags%num_stoch_levels              )
ENDIF
IF ( SIZE(grid%rstoch,1)*SIZE(grid%rstoch,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%rstoch, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1, config_flags%num_stoch_levels,             &
ims, ime, jms, jme, 1, config_flags%num_stoch_levels,             &
ips, ipe, jps, jpe, 1, config_flags%num_stoch_levels              )
ENDIF
IF ( SIZE(grid%pattern_spp_mp,1)*SIZE(grid%pattern_spp_mp,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%pattern_spp_mp, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1, config_flags%num_stoch_levels,             &
ims, ime, jms, jme, 1, config_flags%num_stoch_levels,             &
ips, ipe, jps, jpe, 1, config_flags%num_stoch_levels              )
ENDIF
IF ( SIZE(grid%pattern_isppt_mp,1)*SIZE(grid%pattern_isppt_mp,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%pattern_isppt_mp, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1, config_flags%num_stoch_levels,             &
ims, ime, jms, jme, 1, config_flags%num_stoch_levels,             &
ips, ipe, jps, jpe, 1, config_flags%num_stoch_levels              )
ENDIF
IF ( SIZE(grid%pattern_isppt_pbl,1)*SIZE(grid%pattern_isppt_pbl,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%pattern_isppt_pbl, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1, config_flags%num_stoch_levels,             &
ims, ime, jms, jme, 1, config_flags%num_stoch_levels,             &
ips, ipe, jps, jpe, 1, config_flags%num_stoch_levels              )
ENDIF
IF ( SIZE(grid%pattern_isppt_rad,1)*SIZE(grid%pattern_isppt_rad,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%pattern_isppt_rad, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1, config_flags%num_stoch_levels,             &
ims, ime, jms, jme, 1, config_flags%num_stoch_levels,             &
ips, ipe, jps, jpe, 1, config_flags%num_stoch_levels              )
ENDIF
IF ( SIZE(grid%pattern_isppt_conv,1)*SIZE(grid%pattern_isppt_conv,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%pattern_isppt_conv, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1, config_flags%num_stoch_levels,             &
ims, ime, jms, jme, 1, config_flags%num_stoch_levels,             &
ips, ipe, jps, jpe, 1, config_flags%num_stoch_levels              )
ENDIF
IF ( SIZE(grid%lu_index,1)*SIZE(grid%lu_index,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%lu_index, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%lu_mask,1)*SIZE(grid%lu_mask,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%lu_mask, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%pdyn_smooth,1)*SIZE(grid%pdyn_smooth,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%pdyn_smooth, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%pdyn_parent,1)*SIZE(grid%pdyn_parent,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%pdyn_parent, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%hres_fis,1)*SIZE(grid%hres_fis,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%hres_fis, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%hres_lnd,1)*SIZE(grid%hres_lnd,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%hres_lnd, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%hres_lu_index,1)*SIZE(grid%hres_lu_index,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%hres_lu_index, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%sm,1)*SIZE(grid%sm,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%sm, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%pd,1)*SIZE(grid%pd,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%pd, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%fis,1)*SIZE(grid%fis,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%fis, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%t,1)*SIZE(grid%t,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%t, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%q,1)*SIZE(grid%q,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%q, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%u,1)*SIZE(grid%u,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%u, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%v,1)*SIZE(grid%v,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%v, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%scurx,1)*SIZE(grid%scurx,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%scurx, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%scury,1)*SIZE(grid%scury,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%scury, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%charn,1)*SIZE(grid%charn,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%charn, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%msang,1)*SIZE(grid%msang,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%msang, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%q2,1)*SIZE(grid%q2,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%q2, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%cwm,1)*SIZE(grid%cwm,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%cwm, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%pint,1)*SIZE(grid%pint,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%pint, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
DO itrace = PARAM_FIRST_SCALAR, num_moist
 IF ( SIZE(moist,1)*SIZE(moist,2) .GT. 1 ) THEN
  CALL RSL_LITE_PACK ( local_communicator,&
moist ( grid%sm31,grid%sm32,grid%sm33,itrace),4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
 ENDIF
ENDDO
DO itrace = PARAM_FIRST_SCALAR, num_scalar
 IF ( SIZE(scalar,1)*SIZE(scalar,2) .GT. 1 ) THEN
  CALL RSL_LITE_PACK ( local_communicator,&
scalar ( grid%sm31,grid%sm32,grid%sm33,itrace),4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
 ENDIF
ENDDO
DO itrace = PARAM_FIRST_SCALAR, num_dfi_scalar
 IF ( SIZE(dfi_scalar,1)*SIZE(dfi_scalar,2) .GT. 1 ) THEN
  CALL RSL_LITE_PACK ( local_communicator,&
dfi_scalar ( grid%sm31,grid%sm32,grid%sm33,itrace),4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
 ENDIF
ENDDO
IF ( SIZE(grid%sst,1)*SIZE(grid%sst,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%sst, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
   CALL RSL_LITE_EXCH_Y ( local_communicator , mytask, ntasks, ntasks_x, ntasks_y, &
                          rsl_sendw_m,  rsl_sendw_p, rsl_recvw_m,  rsl_recvw_p    )
IF ( SIZE(grid%rt_tendf_stoch,1)*SIZE(grid%rt_tendf_stoch,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%rt_tendf_stoch, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1, config_flags%num_stoch_levels,             &
ims, ime, jms, jme, 1, config_flags%num_stoch_levels,             &
ips, ipe, jps, jpe, 1, config_flags%num_stoch_levels              )
ENDIF
IF ( SIZE(grid%rand_pert,1)*SIZE(grid%rand_pert,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%rand_pert, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1, config_flags%num_stoch_levels,             &
ims, ime, jms, jme, 1, config_flags%num_stoch_levels,             &
ips, ipe, jps, jpe, 1, config_flags%num_stoch_levels              )
ENDIF
IF ( SIZE(grid%rstoch,1)*SIZE(grid%rstoch,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%rstoch, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1, config_flags%num_stoch_levels,             &
ims, ime, jms, jme, 1, config_flags%num_stoch_levels,             &
ips, ipe, jps, jpe, 1, config_flags%num_stoch_levels              )
ENDIF
IF ( SIZE(grid%pattern_spp_mp,1)*SIZE(grid%pattern_spp_mp,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%pattern_spp_mp, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1, config_flags%num_stoch_levels,             &
ims, ime, jms, jme, 1, config_flags%num_stoch_levels,             &
ips, ipe, jps, jpe, 1, config_flags%num_stoch_levels              )
ENDIF
IF ( SIZE(grid%pattern_isppt_mp,1)*SIZE(grid%pattern_isppt_mp,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%pattern_isppt_mp, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1, config_flags%num_stoch_levels,             &
ims, ime, jms, jme, 1, config_flags%num_stoch_levels,             &
ips, ipe, jps, jpe, 1, config_flags%num_stoch_levels              )
ENDIF
IF ( SIZE(grid%pattern_isppt_pbl,1)*SIZE(grid%pattern_isppt_pbl,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%pattern_isppt_pbl, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1, config_flags%num_stoch_levels,             &
ims, ime, jms, jme, 1, config_flags%num_stoch_levels,             &
ips, ipe, jps, jpe, 1, config_flags%num_stoch_levels              )
ENDIF
IF ( SIZE(grid%pattern_isppt_rad,1)*SIZE(grid%pattern_isppt_rad,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%pattern_isppt_rad, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1, config_flags%num_stoch_levels,             &
ims, ime, jms, jme, 1, config_flags%num_stoch_levels,             &
ips, ipe, jps, jpe, 1, config_flags%num_stoch_levels              )
ENDIF
IF ( SIZE(grid%pattern_isppt_conv,1)*SIZE(grid%pattern_isppt_conv,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%pattern_isppt_conv, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1, config_flags%num_stoch_levels,             &
ims, ime, jms, jme, 1, config_flags%num_stoch_levels,             &
ips, ipe, jps, jpe, 1, config_flags%num_stoch_levels              )
ENDIF
IF ( SIZE(grid%lu_index,1)*SIZE(grid%lu_index,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%lu_index, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%lu_mask,1)*SIZE(grid%lu_mask,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%lu_mask, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%pdyn_smooth,1)*SIZE(grid%pdyn_smooth,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%pdyn_smooth, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%pdyn_parent,1)*SIZE(grid%pdyn_parent,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%pdyn_parent, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%hres_fis,1)*SIZE(grid%hres_fis,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%hres_fis, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%hres_lnd,1)*SIZE(grid%hres_lnd,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%hres_lnd, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%hres_lu_index,1)*SIZE(grid%hres_lu_index,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%hres_lu_index, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%sm,1)*SIZE(grid%sm,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%sm, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%pd,1)*SIZE(grid%pd,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%pd, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%fis,1)*SIZE(grid%fis,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%fis, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%t,1)*SIZE(grid%t,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%t, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%q,1)*SIZE(grid%q,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%q, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%u,1)*SIZE(grid%u,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%u, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%v,1)*SIZE(grid%v,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%v, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%scurx,1)*SIZE(grid%scurx,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%scurx, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%scury,1)*SIZE(grid%scury,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%scury, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%charn,1)*SIZE(grid%charn,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%charn, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%msang,1)*SIZE(grid%msang,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%msang, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%q2,1)*SIZE(grid%q2,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%q2, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%cwm,1)*SIZE(grid%cwm,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%cwm, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%pint,1)*SIZE(grid%pint,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%pint, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
DO itrace = PARAM_FIRST_SCALAR, num_moist
 IF ( SIZE(moist,1)*SIZE(moist,2) .GT. 1 ) THEN
  CALL RSL_LITE_PACK ( local_communicator,&
moist ( grid%sm31,grid%sm32,grid%sm33,itrace),4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
 ENDIF
ENDDO
DO itrace = PARAM_FIRST_SCALAR, num_scalar
 IF ( SIZE(scalar,1)*SIZE(scalar,2) .GT. 1 ) THEN
  CALL RSL_LITE_PACK ( local_communicator,&
scalar ( grid%sm31,grid%sm32,grid%sm33,itrace),4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
 ENDIF
ENDDO
DO itrace = PARAM_FIRST_SCALAR, num_dfi_scalar
 IF ( SIZE(dfi_scalar,1)*SIZE(dfi_scalar,2) .GT. 1 ) THEN
  CALL RSL_LITE_PACK ( local_communicator,&
dfi_scalar ( grid%sm31,grid%sm32,grid%sm33,itrace),4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
 ENDIF
ENDDO
IF ( SIZE(grid%sst,1)*SIZE(grid%sst,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%sst, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 0, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
ENDDO
CALL rsl_comm_iter_init(4,ips,ipe)
DO WHILE ( rsl_comm_iter( grid%id , grid%is_intermediate, 4 , &
                         1 , ids,ide,ips,ipe, grid%nids, grid%nide , & 
     rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p,   & 
     rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p    ))
 CALL RSL_LITE_INIT_EXCH ( local_communicator, 4, 1, &
     rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p,   & 
     rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p,   & 
     15  &
   + num_moist   &
   + num_scalar   &
   + num_dfi_scalar   &
     , 15, RWORDSIZE, &
     0, 0, IWORDSIZE, &
     0, 0, DWORDSIZE, &
      0,  0, LWORDSIZE, &
      mytask, ntasks, ntasks_x, ntasks_y,   &
      ips, ipe, jps, jpe, kps, MAX(1,1&
,config_flags%num_stoch_levels &
,kpe &
))
IF ( SIZE(grid%rt_tendf_stoch,1)*SIZE(grid%rt_tendf_stoch,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%rt_tendf_stoch, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1, config_flags%num_stoch_levels,             &
ims, ime, jms, jme, 1, config_flags%num_stoch_levels,             &
ips, ipe, jps, jpe, 1, config_flags%num_stoch_levels              )
ENDIF
IF ( SIZE(grid%rand_pert,1)*SIZE(grid%rand_pert,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%rand_pert, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1, config_flags%num_stoch_levels,             &
ims, ime, jms, jme, 1, config_flags%num_stoch_levels,             &
ips, ipe, jps, jpe, 1, config_flags%num_stoch_levels              )
ENDIF
IF ( SIZE(grid%rstoch,1)*SIZE(grid%rstoch,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%rstoch, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1, config_flags%num_stoch_levels,             &
ims, ime, jms, jme, 1, config_flags%num_stoch_levels,             &
ips, ipe, jps, jpe, 1, config_flags%num_stoch_levels              )
ENDIF
IF ( SIZE(grid%pattern_spp_mp,1)*SIZE(grid%pattern_spp_mp,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%pattern_spp_mp, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1, config_flags%num_stoch_levels,             &
ims, ime, jms, jme, 1, config_flags%num_stoch_levels,             &
ips, ipe, jps, jpe, 1, config_flags%num_stoch_levels              )
ENDIF
IF ( SIZE(grid%pattern_isppt_mp,1)*SIZE(grid%pattern_isppt_mp,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%pattern_isppt_mp, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1, config_flags%num_stoch_levels,             &
ims, ime, jms, jme, 1, config_flags%num_stoch_levels,             &
ips, ipe, jps, jpe, 1, config_flags%num_stoch_levels              )
ENDIF
IF ( SIZE(grid%pattern_isppt_pbl,1)*SIZE(grid%pattern_isppt_pbl,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%pattern_isppt_pbl, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1, config_flags%num_stoch_levels,             &
ims, ime, jms, jme, 1, config_flags%num_stoch_levels,             &
ips, ipe, jps, jpe, 1, config_flags%num_stoch_levels              )
ENDIF
IF ( SIZE(grid%pattern_isppt_rad,1)*SIZE(grid%pattern_isppt_rad,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%pattern_isppt_rad, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1, config_flags%num_stoch_levels,             &
ims, ime, jms, jme, 1, config_flags%num_stoch_levels,             &
ips, ipe, jps, jpe, 1, config_flags%num_stoch_levels              )
ENDIF
IF ( SIZE(grid%pattern_isppt_conv,1)*SIZE(grid%pattern_isppt_conv,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%pattern_isppt_conv, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1, config_flags%num_stoch_levels,             &
ims, ime, jms, jme, 1, config_flags%num_stoch_levels,             &
ips, ipe, jps, jpe, 1, config_flags%num_stoch_levels              )
ENDIF
IF ( SIZE(grid%lu_index,1)*SIZE(grid%lu_index,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%lu_index, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%lu_mask,1)*SIZE(grid%lu_mask,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%lu_mask, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%pdyn_smooth,1)*SIZE(grid%pdyn_smooth,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%pdyn_smooth, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%pdyn_parent,1)*SIZE(grid%pdyn_parent,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%pdyn_parent, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%hres_fis,1)*SIZE(grid%hres_fis,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%hres_fis, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%hres_lnd,1)*SIZE(grid%hres_lnd,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%hres_lnd, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%hres_lu_index,1)*SIZE(grid%hres_lu_index,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%hres_lu_index, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%sm,1)*SIZE(grid%sm,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%sm, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%pd,1)*SIZE(grid%pd,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%pd, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%fis,1)*SIZE(grid%fis,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%fis, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%t,1)*SIZE(grid%t,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%t, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%q,1)*SIZE(grid%q,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%q, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%u,1)*SIZE(grid%u,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%u, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%v,1)*SIZE(grid%v,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%v, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%scurx,1)*SIZE(grid%scurx,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%scurx, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%scury,1)*SIZE(grid%scury,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%scury, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%charn,1)*SIZE(grid%charn,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%charn, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%msang,1)*SIZE(grid%msang,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%msang, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%q2,1)*SIZE(grid%q2,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%q2, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%cwm,1)*SIZE(grid%cwm,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%cwm, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%pint,1)*SIZE(grid%pint,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%pint, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
DO itrace = PARAM_FIRST_SCALAR, num_moist
 IF ( SIZE(moist,1)*SIZE(moist,2) .GT. 1 ) THEN
  CALL RSL_LITE_PACK ( local_communicator,&
moist ( grid%sm31,grid%sm32,grid%sm33,itrace),4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
 ENDIF
ENDDO
DO itrace = PARAM_FIRST_SCALAR, num_scalar
 IF ( SIZE(scalar,1)*SIZE(scalar,2) .GT. 1 ) THEN
  CALL RSL_LITE_PACK ( local_communicator,&
scalar ( grid%sm31,grid%sm32,grid%sm33,itrace),4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
 ENDIF
ENDDO
DO itrace = PARAM_FIRST_SCALAR, num_dfi_scalar
 IF ( SIZE(dfi_scalar,1)*SIZE(dfi_scalar,2) .GT. 1 ) THEN
  CALL RSL_LITE_PACK ( local_communicator,&
dfi_scalar ( grid%sm31,grid%sm32,grid%sm33,itrace),4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
 ENDIF
ENDDO
IF ( SIZE(grid%sst,1)*SIZE(grid%sst,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%sst, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
   CALL RSL_LITE_EXCH_X ( local_communicator , mytask, ntasks, ntasks_x, ntasks_y, &
                          rsl_sendw_m,  rsl_sendw_p, rsl_recvw_m,  rsl_recvw_p    )
IF ( SIZE(grid%rt_tendf_stoch,1)*SIZE(grid%rt_tendf_stoch,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%rt_tendf_stoch, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1, config_flags%num_stoch_levels,             &
ims, ime, jms, jme, 1, config_flags%num_stoch_levels,             &
ips, ipe, jps, jpe, 1, config_flags%num_stoch_levels              )
ENDIF
IF ( SIZE(grid%rand_pert,1)*SIZE(grid%rand_pert,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%rand_pert, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1, config_flags%num_stoch_levels,             &
ims, ime, jms, jme, 1, config_flags%num_stoch_levels,             &
ips, ipe, jps, jpe, 1, config_flags%num_stoch_levels              )
ENDIF
IF ( SIZE(grid%rstoch,1)*SIZE(grid%rstoch,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%rstoch, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1, config_flags%num_stoch_levels,             &
ims, ime, jms, jme, 1, config_flags%num_stoch_levels,             &
ips, ipe, jps, jpe, 1, config_flags%num_stoch_levels              )
ENDIF
IF ( SIZE(grid%pattern_spp_mp,1)*SIZE(grid%pattern_spp_mp,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%pattern_spp_mp, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1, config_flags%num_stoch_levels,             &
ims, ime, jms, jme, 1, config_flags%num_stoch_levels,             &
ips, ipe, jps, jpe, 1, config_flags%num_stoch_levels              )
ENDIF
IF ( SIZE(grid%pattern_isppt_mp,1)*SIZE(grid%pattern_isppt_mp,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%pattern_isppt_mp, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1, config_flags%num_stoch_levels,             &
ims, ime, jms, jme, 1, config_flags%num_stoch_levels,             &
ips, ipe, jps, jpe, 1, config_flags%num_stoch_levels              )
ENDIF
IF ( SIZE(grid%pattern_isppt_pbl,1)*SIZE(grid%pattern_isppt_pbl,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%pattern_isppt_pbl, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1, config_flags%num_stoch_levels,             &
ims, ime, jms, jme, 1, config_flags%num_stoch_levels,             &
ips, ipe, jps, jpe, 1, config_flags%num_stoch_levels              )
ENDIF
IF ( SIZE(grid%pattern_isppt_rad,1)*SIZE(grid%pattern_isppt_rad,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%pattern_isppt_rad, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1, config_flags%num_stoch_levels,             &
ims, ime, jms, jme, 1, config_flags%num_stoch_levels,             &
ips, ipe, jps, jpe, 1, config_flags%num_stoch_levels              )
ENDIF
IF ( SIZE(grid%pattern_isppt_conv,1)*SIZE(grid%pattern_isppt_conv,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%pattern_isppt_conv, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1, config_flags%num_stoch_levels,             &
ims, ime, jms, jme, 1, config_flags%num_stoch_levels,             &
ips, ipe, jps, jpe, 1, config_flags%num_stoch_levels              )
ENDIF
IF ( SIZE(grid%lu_index,1)*SIZE(grid%lu_index,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%lu_index, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%lu_mask,1)*SIZE(grid%lu_mask,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%lu_mask, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%pdyn_smooth,1)*SIZE(grid%pdyn_smooth,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%pdyn_smooth, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%pdyn_parent,1)*SIZE(grid%pdyn_parent,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%pdyn_parent, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%hres_fis,1)*SIZE(grid%hres_fis,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%hres_fis, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%hres_lnd,1)*SIZE(grid%hres_lnd,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%hres_lnd, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%hres_lu_index,1)*SIZE(grid%hres_lu_index,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%hres_lu_index, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%sm,1)*SIZE(grid%sm,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%sm, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%pd,1)*SIZE(grid%pd,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%pd, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%fis,1)*SIZE(grid%fis,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%fis, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%t,1)*SIZE(grid%t,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%t, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%q,1)*SIZE(grid%q,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%q, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%u,1)*SIZE(grid%u,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%u, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%v,1)*SIZE(grid%v,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%v, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%scurx,1)*SIZE(grid%scurx,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%scurx, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%scury,1)*SIZE(grid%scury,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%scury, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%charn,1)*SIZE(grid%charn,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%charn, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%msang,1)*SIZE(grid%msang,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%msang, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%q2,1)*SIZE(grid%q2,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%q2, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%cwm,1)*SIZE(grid%cwm,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%cwm, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%pint,1)*SIZE(grid%pint,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%pint, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
DO itrace = PARAM_FIRST_SCALAR, num_moist
 IF ( SIZE(moist,1)*SIZE(moist,2) .GT. 1 ) THEN
  CALL RSL_LITE_PACK ( local_communicator,&
moist ( grid%sm31,grid%sm32,grid%sm33,itrace),4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
 ENDIF
ENDDO
DO itrace = PARAM_FIRST_SCALAR, num_scalar
 IF ( SIZE(scalar,1)*SIZE(scalar,2) .GT. 1 ) THEN
  CALL RSL_LITE_PACK ( local_communicator,&
scalar ( grid%sm31,grid%sm32,grid%sm33,itrace),4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
 ENDIF
ENDDO
DO itrace = PARAM_FIRST_SCALAR, num_dfi_scalar
 IF ( SIZE(dfi_scalar,1)*SIZE(dfi_scalar,2) .GT. 1 ) THEN
  CALL RSL_LITE_PACK ( local_communicator,&
dfi_scalar ( grid%sm31,grid%sm32,grid%sm33,itrace),4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
 ENDIF
ENDDO
IF ( SIZE(grid%sst,1)*SIZE(grid%sst,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%sst, 4,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
RWORDSIZE, 1, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
    ENDDO
!ENDOFREGISTRYGENERATEDINCLUDE
