!STARTOFREGISTRYGENERATEDINCLUDE 'inc/nest_interpdown_unpack.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL rsl_lite_from_parent_info(pig,pjg,retval)
DO while ( retval .eq. 1 )
CALL rsl_lite_from_parent_msg(-1,xv) ! -1 = all data
ixv=0 ! Reset read/write counter
IF ( SIZE(grid%lakedepth2d) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%lakedepth2d(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%savedtke12d) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%savedtke12d(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%snowdp2d) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%snowdp2d(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%h2osno2d) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%h2osno2d(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%snl2d) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%snl2d(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%t_grnd2d) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%t_grnd2d(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%t_lake3d) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,10
grid%t_lake3d(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%lake_icefrac3d) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,10
grid%lake_icefrac3d(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%z_lake3d) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,10
grid%z_lake3d(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%dz_lake3d) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,10
grid%dz_lake3d(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%t_soisno3d) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,15
grid%t_soisno3d(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+15-1+1
ENDIF
IF ( SIZE(grid%h2osoi_ice3d) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,15
grid%h2osoi_ice3d(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+15-1+1
ENDIF
IF ( SIZE(grid%h2osoi_liq3d) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,15
grid%h2osoi_liq3d(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+15-1+1
ENDIF
IF ( SIZE(grid%h2osoi_vol3d) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,15
grid%h2osoi_vol3d(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+15-1+1
ENDIF
IF ( SIZE(grid%z3d) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,15
grid%z3d(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+15-1+1
ENDIF
IF ( SIZE(grid%dz3d) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,15
grid%dz3d(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+15-1+1
ENDIF
IF ( SIZE(grid%zi3d) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,16
grid%zi3d(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+16-1+1
ENDIF
IF ( SIZE(grid%watsat3d) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,10
grid%watsat3d(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%csol3d) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,10
grid%csol3d(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%tkmg3d) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,10
grid%tkmg3d(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%tkdry3d) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,10
grid%tkdry3d(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%tksatu3d) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,10
grid%tksatu3d(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%vertstrucc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,(ckde-1)
grid%vertstrucc(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%vertstrucs) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,(ckde-1)
grid%vertstrucs(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%ru_tendf_stoch) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,config_flags%num_stoch_levels
grid%ru_tendf_stoch(pig,pjg,k) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+config_flags%num_stoch_levels-1+1
ENDIF
IF ( SIZE(grid%rv_tendf_stoch) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,config_flags%num_stoch_levels
grid%rv_tendf_stoch(pig,pjg,k) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+config_flags%num_stoch_levels-1+1
ENDIF
IF ( SIZE(grid%lu_index) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%lu_index(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%lu_mask) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%lu_mask(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%cuprecip_swath) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%cuprecip_swath(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%precip_swath) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%precip_swath(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%windsq_swath) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%windsq_swath(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%suswind) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%suswind(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%wind10_ratio) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%wind10_ratio(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%suswind_swath) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%suswind_swath(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%membrane_mslp) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%membrane_mslp(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%pdyn_smooth) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%pdyn_smooth(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%pdyn_parent) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%pdyn_parent(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%best_mslp) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%best_mslp(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%hlon) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%hlon(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%hlat) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%hlat(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tg_max_m10wind) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%tg_max_m10wind(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tg_max_wwind) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%tg_max_wwind(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tg_min_wwind) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%tg_min_wwind(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tg_max_zhel_25) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%tg_max_zhel_25(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tg_min_zhel_25) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%tg_min_zhel_25(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tg_max_zhel_03) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%tg_max_zhel_03(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tg_min_zhel_03) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%tg_min_zhel_03(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tg_updhel25) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%tg_updhel25(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tg_max_updhel25) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%tg_max_updhel25(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tg_updhel03) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%tg_updhel03(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tg_max_updhel03) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%tg_max_updhel03(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tg_total_precip) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%tg_total_precip(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%hres_fis) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%hres_fis(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%hres_lnd) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%hres_lnd(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%hres_lu_index) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%hres_lu_index(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%sm) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%sm(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%sice) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%sice(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%pd) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%pd(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%fis) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%fis(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%t) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,(ckde-1)
grid%t(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%q) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,(ckde-1)
grid%q(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%u) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,(ckde-1)
grid%u(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%v) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,(ckde-1)
grid%v(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%ustar) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%ustar(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%z0) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%z0(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%scurx) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%scurx(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%scury) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%scury(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%charn) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%charn(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%msang) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%msang(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%ths) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%ths(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%qsh) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%qsh(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%taux) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%taux(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tauy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%tauy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%acprec) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%acprec(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%cuprec) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%cuprec(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%cldefi) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%cldefi(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%th10) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%th10(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%q10) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%q10(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%pshltr) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%pshltr(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tshltr) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%tshltr(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%qshltr) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%qshltr(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%q2) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,(ckde-1)
grid%q2(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%t_adj) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,(ckde-1)
grid%t_adj(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%cd_out) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%cd_out(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%ch_out) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%ch_out(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%albase) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%albase(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%nmm_tsk) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%nmm_tsk(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%mxsnal) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%mxsnal(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%sigt4) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%sigt4(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tg) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%tg(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%acswupt) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%acswupt(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%acswuptc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%acswuptc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%acswdnt) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%acswdnt(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%acswdntc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%acswdntc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%acswupb) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%acswupb(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%acswupbc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%acswupbc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%acswdnb) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%acswdnb(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%acswdnbc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%acswdnbc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%aclwupt) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%aclwupt(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%aclwuptc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%aclwuptc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%aclwdnt) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%aclwdnt(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%aclwdntc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%aclwdntc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%aclwupb) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%aclwupb(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%aclwupbc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%aclwupbc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%aclwdnb) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%aclwdnb(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%aclwdnbc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%aclwdnbc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%swupt) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%swupt(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%swuptc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%swuptc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%swdnt) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%swdnt(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%swdntc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%swdntc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%swupb) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%swupb(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%swupbc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%swupbc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%swdnb) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%swdnb(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%swdnbc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%swdnbc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%lwupt) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%lwupt(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%lwuptc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%lwuptc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%lwdnt) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%lwdnt(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%lwdntc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%lwdntc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%lwupb) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%lwupb(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%lwupbc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%lwupbc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%lwdnb) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%lwdnb(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%lwdnbc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%lwdnbc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%refl_10cm) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,(ckde-1)
grid%refl_10cm(pig,k,pjg) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%refd_max) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%refd_max(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%qnwfa2d) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%qnwfa2d(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%swddir) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%swddir(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%swddirc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%swddirc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%swddni) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%swddni(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%swddnic) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%swddnic(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%swddif) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%swddif(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%gx) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%gx(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%bx) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%bx(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%gg) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%gg(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%bb) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%bb(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%coszen_ref) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%coszen_ref(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%swdown_ref) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%swdown_ref(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%swddir_ref) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%swddir_ref(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%cwm) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,(ckde-1)
grid%cwm(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%f_ice) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,(ckde-1)
grid%f_ice(pig,k,pjg) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%f_rain) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,(ckde-1)
grid%f_rain(pig,k,pjg) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%f_rimef) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,(ckde-1)
grid%f_rimef(pig,k,pjg) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%cfrach) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%cfrach(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%cfracl) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%cfracl(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%cfracm) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%cfracm(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%islope) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%islope(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%cmc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%cmc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%soiltb) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%soiltb(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%vegfrc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%vegfrc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%shdmax) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%shdmax(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%shdmin) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%shdmin(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%sh2o) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,config_flags%num_soil_layers
grid%sh2o(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+config_flags%num_soil_layers-1+1
ENDIF
IF ( SIZE(grid%smc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,config_flags%num_soil_layers
grid%smc(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+config_flags%num_soil_layers-1+1
ENDIF
IF ( SIZE(grid%stc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,config_flags%num_soil_layers
grid%stc(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+config_flags%num_soil_layers-1+1
ENDIF
IF ( SIZE(grid%ctopo) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%ctopo(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%ctopo2) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%ctopo2(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%dwdt) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,(ckde-1)
grid%dwdt(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%pint) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,ckde
grid%pint(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+ckde-ckds+1
ENDIF
IF ( SIZE(grid%w) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,ckde
grid%w(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+ckde-ckds+1
ENDIF
IF ( SIZE(grid%w_tot) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,ckde
grid%w_tot(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+ckde-ckds+1
ENDIF
IF ( SIZE(grid%z) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,ckde
grid%z(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+ckde-ckds+1
ENDIF
IF ( SIZE(grid%rlwin) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%rlwin(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%rswin) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%rswin(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%rlwtt) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,(ckde-1)
grid%rlwtt(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%rswtt) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,(ckde-1)
grid%rswtt(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%winfo) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,ckde
grid%winfo(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+ckde-ckds+1
ENDIF
IF ( SIZE(grid%iinfo) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,ckde
grid%iinfo(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+ckde-ckds+1
ENDIF
IF ( SIZE(grid%landmask) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%landmask(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%toposoil) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%toposoil(pig,pjg) = xv(ixv)
ENDIF
if(interp_mp .eqv. .true.) then
DO itrace =  PARAM_FIRST_SCALAR, num_moist
DO k = ckds,(ckde-1)
moist(pig,pjg,k,itrace) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDDO
endif
if(interp_mp .eqv. .true.) then
DO itrace =  PARAM_FIRST_SCALAR, num_scalar
DO k = ckds,(ckde-1)
scalar(pig,pjg,k,itrace) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDDO
endif
if(interp_mp .eqv. .true.) then
DO itrace =  PARAM_FIRST_SCALAR, num_dfi_scalar
DO k = ckds,(ckde-1)
dfi_scalar(pig,pjg,k,itrace) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDDO
endif
IF ( SIZE(grid%lake_depth) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%lake_depth(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%u10) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%u10(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%v10) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%v10(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%xice) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%xice(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%icedepth) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%icedepth(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%albsi) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%albsi(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%snowsi) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%snowsi(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%ivgtyp) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%ivgtyp(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%isltyp) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%isltyp(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%vegfra) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%vegfra(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%sst) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%sst(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%weasd) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%weasd(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%thz0) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%thz0(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%qz0) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%qz0(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%uz0) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%uz0(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%vz0) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%vz0(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%htop) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%htop(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%hbot) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%hbot(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%htopr) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%htopr(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%hbotr) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%hbotr(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%cuppt) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%cuppt(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%snowh) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%snowh(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%rhosn) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%rhosn(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%isnowxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%isnowxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tvxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%tvxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tgxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%tgxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%canicexy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%canicexy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%canliqxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%canliqxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%eahxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%eahxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tahxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%tahxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%cmxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%cmxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%chxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%chxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%fwetxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%fwetxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%sneqvoxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%sneqvoxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%alboldxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%alboldxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%qsnowxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%qsnowxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%wslakexy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%wslakexy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%zwtxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%zwtxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%waxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%waxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%wtxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%wtxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tsnoxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,config_flags%num_snow_layers
grid%tsnoxy(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+config_flags%num_snow_layers-1+1
ENDIF
IF ( SIZE(grid%zsnsoxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,config_flags%num_snso_layers
grid%zsnsoxy(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+config_flags%num_snso_layers-1+1
ENDIF
IF ( SIZE(grid%snicexy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,config_flags%num_snow_layers
grid%snicexy(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+config_flags%num_snow_layers-1+1
ENDIF
IF ( SIZE(grid%snliqxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,config_flags%num_snow_layers
grid%snliqxy(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+config_flags%num_snow_layers-1+1
ENDIF
IF ( SIZE(grid%lfmassxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%lfmassxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%rtmassxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%rtmassxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%stmassxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%stmassxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%woodxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%woodxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%stblcpxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%stblcpxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%fastcpxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%fastcpxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%xsaixy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%xsaixy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%t2mvxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%t2mvxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%t2mbxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%t2mbxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%q2mvxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%q2mvxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%q2mbxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%q2mbxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tradxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%tradxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%neexy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%neexy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%gppxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%gppxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%nppxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%nppxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%fvegxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%fvegxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%qinxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%qinxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%runsfxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%runsfxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%runsbxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%runsbxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%ecanxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%ecanxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%edirxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%edirxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%etranxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%etranxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%fsaxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%fsaxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%firaxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%firaxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%aparxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%aparxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%psnxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%psnxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%savxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%savxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%sagxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%sagxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%rssunxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%rssunxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%rsshaxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%rsshaxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%bgapxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%bgapxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%wgapxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%wgapxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tgvxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%tgvxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tgbxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%tgbxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%chvxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%chvxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%chbxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%chbxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%shgxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%shgxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%shcxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%shcxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%shbxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%shbxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%evgxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%evgxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%evbxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%evbxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%ghvxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%ghvxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%ghbxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%ghbxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%irgxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%irgxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%ircxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%ircxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%irbxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%irbxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%trxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%trxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%evcxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%evcxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%chleafxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%chleafxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%chucxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%chucxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%chv2xy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%chv2xy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%chb2xy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%chb2xy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%chstarxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%chstarxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%soilcomp) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,8
grid%soilcomp(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+8-1+1
ENDIF
IF ( SIZE(grid%soilcl1) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%soilcl1(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%soilcl2) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%soilcl2(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%soilcl3) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%soilcl3(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%soilcl4) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%soilcl4(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%grainxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%grainxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%gddxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%gddxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%croptype) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,5
grid%croptype(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+5-1+1
ENDIF
IF ( SIZE(grid%planting) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%planting(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%harvest) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%harvest(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%season_gdd) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%season_gdd(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%gecros_state) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,60
grid%gecros_state(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+60-1+1
ENDIF
CALL rsl_lite_from_parent_info(pig,pjg,retval)
ENDDO
!ENDOFREGISTRYGENERATEDINCLUDE
