!STARTOFREGISTRYGENERATEDINCLUDE 'inc/nest_interpdown_pack.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
msize = (50 + ((num_szj - PARAM_FIRST_SCALAR + 1)) & 
 + ((num_s1z - PARAM_FIRST_SCALAR + 1)) & 
 + ((num_spz - PARAM_FIRST_SCALAR + 1)) & 
 + ((num_tcs - PARAM_FIRST_SCALAR + 1)) & 
 + ((num_chem - PARAM_FIRST_SCALAR + 1)) & 
 + ((num_ozmixm - PARAM_FIRST_SCALAR + 1)) )* nlev + 227
IF(interp_mp .eqv. .true.) then
    msize=msize + (0 + ((num_moist - PARAM_FIRST_SCALAR + 1)) & 
 + ((num_dfi_moist - PARAM_FIRST_SCALAR + 1)) & 
 + ((num_scalar - PARAM_FIRST_SCALAR + 1)) & 
 + ((num_dfi_scalar - PARAM_FIRST_SCALAR + 1)) )*nlev+0
ENDIF
CALL rsl_lite_to_child_info( msize*RWORDSIZE                               &
                        ,cips,cipe,cjps,cjpe                               &
                        ,iids,iide,ijds,ijde                               &
                        ,nids,nide,njds,njde                               &
                        ,pgr , sw                                          &
                        ,nest_task_offsets(ngrid%id)                      &
                        ,nest_pes_x(grid%id)                            &
                        ,nest_pes_y(grid%id)                            &
                        ,nest_pes_x(intermediate_grid%id)                 &
                        ,nest_pes_y(intermediate_grid%id)                 &
                        ,thisdomain_max_halo_width                         &
                        ,icoord,jcoord                                     &
                        ,idim_cd,jdim_cd                                   &
                        ,pig,pjg,retval )
DO while ( retval .eq. 1 )
ixv=0 ! Reset read/write counter
IF ( SIZE(grid%lakedepth2d) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%lakedepth2d(pig,pjg)
ENDIF
IF ( SIZE(grid%savedtke12d) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%savedtke12d(pig,pjg)
ENDIF
IF ( SIZE(grid%snowdp2d) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%snowdp2d(pig,pjg)
ENDIF
IF ( SIZE(grid%h2osno2d) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%h2osno2d(pig,pjg)
ENDIF
IF ( SIZE(grid%snl2d) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%snl2d(pig,pjg)
ENDIF
IF ( SIZE(grid%t_grnd2d) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%t_grnd2d(pig,pjg)
ENDIF
IF ( SIZE(grid%t_lake3d) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,10
xv(k+ixv+1-1)= grid%t_lake3d(pig,k,pjg)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%lake_icefrac3d) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,10
xv(k+ixv+1-1)= grid%lake_icefrac3d(pig,k,pjg)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%z_lake3d) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,10
xv(k+ixv+1-1)= grid%z_lake3d(pig,k,pjg)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%dz_lake3d) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,10
xv(k+ixv+1-1)= grid%dz_lake3d(pig,k,pjg)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%t_soisno3d) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,15
xv(k+ixv+1-1)= grid%t_soisno3d(pig,k,pjg)
ENDDO
ixv=ixv+15-1+1
ENDIF
IF ( SIZE(grid%h2osoi_ice3d) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,15
xv(k+ixv+1-1)= grid%h2osoi_ice3d(pig,k,pjg)
ENDDO
ixv=ixv+15-1+1
ENDIF
IF ( SIZE(grid%h2osoi_liq3d) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,15
xv(k+ixv+1-1)= grid%h2osoi_liq3d(pig,k,pjg)
ENDDO
ixv=ixv+15-1+1
ENDIF
IF ( SIZE(grid%h2osoi_vol3d) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,15
xv(k+ixv+1-1)= grid%h2osoi_vol3d(pig,k,pjg)
ENDDO
ixv=ixv+15-1+1
ENDIF
IF ( SIZE(grid%z3d) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,15
xv(k+ixv+1-1)= grid%z3d(pig,k,pjg)
ENDDO
ixv=ixv+15-1+1
ENDIF
IF ( SIZE(grid%dz3d) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,15
xv(k+ixv+1-1)= grid%dz3d(pig,k,pjg)
ENDDO
ixv=ixv+15-1+1
ENDIF
IF ( SIZE(grid%zi3d) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,16
xv(k+ixv+1-1)= grid%zi3d(pig,k,pjg)
ENDDO
ixv=ixv+16-1+1
ENDIF
IF ( SIZE(grid%watsat3d) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,10
xv(k+ixv+1-1)= grid%watsat3d(pig,k,pjg)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%csol3d) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,10
xv(k+ixv+1-1)= grid%csol3d(pig,k,pjg)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%tkmg3d) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,10
xv(k+ixv+1-1)= grid%tkmg3d(pig,k,pjg)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%tkdry3d) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,10
xv(k+ixv+1-1)= grid%tkdry3d(pig,k,pjg)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%tksatu3d) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,10
xv(k+ixv+1-1)= grid%tksatu3d(pig,k,pjg)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%vertstrucc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,(ckde-1)
xv(k+ixv+1-ckds)= grid%vertstrucc(pig,pjg,k)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%vertstrucs) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,(ckde-1)
xv(k+ixv+1-ckds)= grid%vertstrucs(pig,pjg,k)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%ru_tendf_stoch) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,config_flags%num_stoch_levels
xv(k+ixv+1-1)= grid%ru_tendf_stoch(pig,pjg,k)
ENDDO
ixv=ixv+config_flags%num_stoch_levels-1+1
ENDIF
IF ( SIZE(grid%rv_tendf_stoch) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,config_flags%num_stoch_levels
xv(k+ixv+1-1)= grid%rv_tendf_stoch(pig,pjg,k)
ENDDO
ixv=ixv+config_flags%num_stoch_levels-1+1
ENDIF
IF ( SIZE(grid%lu_index) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%lu_index(pig,pjg)
ENDIF
IF ( SIZE(grid%lu_mask) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%lu_mask(pig,pjg)
ENDIF
IF ( SIZE(grid%cuprecip_swath) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%cuprecip_swath(pig,pjg)
ENDIF
IF ( SIZE(grid%precip_swath) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%precip_swath(pig,pjg)
ENDIF
IF ( SIZE(grid%windsq_swath) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%windsq_swath(pig,pjg)
ENDIF
IF ( SIZE(grid%suswind) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%suswind(pig,pjg)
ENDIF
IF ( SIZE(grid%wind10_ratio) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%wind10_ratio(pig,pjg)
ENDIF
IF ( SIZE(grid%suswind_swath) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%suswind_swath(pig,pjg)
ENDIF
IF ( SIZE(grid%membrane_mslp) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%membrane_mslp(pig,pjg)
ENDIF
IF ( SIZE(grid%pdyn_smooth) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%pdyn_smooth(pig,pjg)
ENDIF
IF ( SIZE(grid%pdyn_parent) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%pdyn_parent(pig,pjg)
ENDIF
IF ( SIZE(grid%best_mslp) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%best_mslp(pig,pjg)
ENDIF
IF ( SIZE(grid%hlon) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%hlon(pig,pjg)
ENDIF
IF ( SIZE(grid%hlat) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%hlat(pig,pjg)
ENDIF
IF ( SIZE(grid%tg_max_m10wind) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%tg_max_m10wind(pig,pjg)
ENDIF
IF ( SIZE(grid%tg_max_wwind) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%tg_max_wwind(pig,pjg)
ENDIF
IF ( SIZE(grid%tg_min_wwind) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%tg_min_wwind(pig,pjg)
ENDIF
IF ( SIZE(grid%tg_max_zhel_25) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%tg_max_zhel_25(pig,pjg)
ENDIF
IF ( SIZE(grid%tg_min_zhel_25) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%tg_min_zhel_25(pig,pjg)
ENDIF
IF ( SIZE(grid%tg_max_zhel_03) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%tg_max_zhel_03(pig,pjg)
ENDIF
IF ( SIZE(grid%tg_min_zhel_03) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%tg_min_zhel_03(pig,pjg)
ENDIF
IF ( SIZE(grid%tg_updhel25) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%tg_updhel25(pig,pjg)
ENDIF
IF ( SIZE(grid%tg_max_updhel25) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%tg_max_updhel25(pig,pjg)
ENDIF
IF ( SIZE(grid%tg_updhel03) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%tg_updhel03(pig,pjg)
ENDIF
IF ( SIZE(grid%tg_max_updhel03) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%tg_max_updhel03(pig,pjg)
ENDIF
IF ( SIZE(grid%tg_total_precip) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%tg_total_precip(pig,pjg)
ENDIF
IF ( SIZE(grid%hres_fis) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%hres_fis(pig,pjg)
ENDIF
IF ( SIZE(grid%hres_lnd) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%hres_lnd(pig,pjg)
ENDIF
IF ( SIZE(grid%hres_lu_index) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%hres_lu_index(pig,pjg)
ENDIF
IF ( SIZE(grid%sm) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%sm(pig,pjg)
ENDIF
IF ( SIZE(grid%sice) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%sice(pig,pjg)
ENDIF
IF ( SIZE(grid%pd) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%pd(pig,pjg)
ENDIF
IF ( SIZE(grid%fis) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%fis(pig,pjg)
ENDIF
IF ( SIZE(grid%t) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,(ckde-1)
xv(k+ixv+1-ckds)= grid%t(pig,pjg,k)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%q) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,(ckde-1)
xv(k+ixv+1-ckds)= grid%q(pig,pjg,k)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%u) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,(ckde-1)
xv(k+ixv+1-ckds)= grid%u(pig,pjg,k)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%v) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,(ckde-1)
xv(k+ixv+1-ckds)= grid%v(pig,pjg,k)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%ustar) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%ustar(pig,pjg)
ENDIF
IF ( SIZE(grid%z0) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%z0(pig,pjg)
ENDIF
IF ( SIZE(grid%scurx) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%scurx(pig,pjg)
ENDIF
IF ( SIZE(grid%scury) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%scury(pig,pjg)
ENDIF
IF ( SIZE(grid%charn) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%charn(pig,pjg)
ENDIF
IF ( SIZE(grid%msang) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%msang(pig,pjg)
ENDIF
IF ( SIZE(grid%ths) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%ths(pig,pjg)
ENDIF
IF ( SIZE(grid%qsh) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%qsh(pig,pjg)
ENDIF
IF ( SIZE(grid%taux) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%taux(pig,pjg)
ENDIF
IF ( SIZE(grid%tauy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%tauy(pig,pjg)
ENDIF
IF ( SIZE(grid%acprec) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%acprec(pig,pjg)
ENDIF
IF ( SIZE(grid%cuprec) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%cuprec(pig,pjg)
ENDIF
IF ( SIZE(grid%cldefi) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%cldefi(pig,pjg)
ENDIF
IF ( SIZE(grid%th10) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%th10(pig,pjg)
ENDIF
IF ( SIZE(grid%q10) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%q10(pig,pjg)
ENDIF
IF ( SIZE(grid%pshltr) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%pshltr(pig,pjg)
ENDIF
IF ( SIZE(grid%tshltr) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%tshltr(pig,pjg)
ENDIF
IF ( SIZE(grid%qshltr) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%qshltr(pig,pjg)
ENDIF
IF ( SIZE(grid%q2) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,(ckde-1)
xv(k+ixv+1-ckds)= grid%q2(pig,pjg,k)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%t_adj) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,(ckde-1)
xv(k+ixv+1-ckds)= grid%t_adj(pig,pjg,k)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%cd_out) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%cd_out(pig,pjg)
ENDIF
IF ( SIZE(grid%ch_out) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%ch_out(pig,pjg)
ENDIF
IF ( SIZE(grid%albase) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%albase(pig,pjg)
ENDIF
IF ( SIZE(grid%nmm_tsk) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%nmm_tsk(pig,pjg)
ENDIF
IF ( SIZE(grid%mxsnal) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%mxsnal(pig,pjg)
ENDIF
IF ( SIZE(grid%sigt4) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%sigt4(pig,pjg)
ENDIF
IF ( SIZE(grid%tg) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%tg(pig,pjg)
ENDIF
IF ( SIZE(grid%acswupt) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%acswupt(pig,pjg)
ENDIF
IF ( SIZE(grid%acswuptc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%acswuptc(pig,pjg)
ENDIF
IF ( SIZE(grid%acswdnt) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%acswdnt(pig,pjg)
ENDIF
IF ( SIZE(grid%acswdntc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%acswdntc(pig,pjg)
ENDIF
IF ( SIZE(grid%acswupb) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%acswupb(pig,pjg)
ENDIF
IF ( SIZE(grid%acswupbc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%acswupbc(pig,pjg)
ENDIF
IF ( SIZE(grid%acswdnb) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%acswdnb(pig,pjg)
ENDIF
IF ( SIZE(grid%acswdnbc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%acswdnbc(pig,pjg)
ENDIF
IF ( SIZE(grid%aclwupt) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%aclwupt(pig,pjg)
ENDIF
IF ( SIZE(grid%aclwuptc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%aclwuptc(pig,pjg)
ENDIF
IF ( SIZE(grid%aclwdnt) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%aclwdnt(pig,pjg)
ENDIF
IF ( SIZE(grid%aclwdntc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%aclwdntc(pig,pjg)
ENDIF
IF ( SIZE(grid%aclwupb) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%aclwupb(pig,pjg)
ENDIF
IF ( SIZE(grid%aclwupbc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%aclwupbc(pig,pjg)
ENDIF
IF ( SIZE(grid%aclwdnb) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%aclwdnb(pig,pjg)
ENDIF
IF ( SIZE(grid%aclwdnbc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%aclwdnbc(pig,pjg)
ENDIF
IF ( SIZE(grid%swupt) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%swupt(pig,pjg)
ENDIF
IF ( SIZE(grid%swuptc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%swuptc(pig,pjg)
ENDIF
IF ( SIZE(grid%swdnt) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%swdnt(pig,pjg)
ENDIF
IF ( SIZE(grid%swdntc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%swdntc(pig,pjg)
ENDIF
IF ( SIZE(grid%swupb) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%swupb(pig,pjg)
ENDIF
IF ( SIZE(grid%swupbc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%swupbc(pig,pjg)
ENDIF
IF ( SIZE(grid%swdnb) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%swdnb(pig,pjg)
ENDIF
IF ( SIZE(grid%swdnbc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%swdnbc(pig,pjg)
ENDIF
IF ( SIZE(grid%lwupt) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%lwupt(pig,pjg)
ENDIF
IF ( SIZE(grid%lwuptc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%lwuptc(pig,pjg)
ENDIF
IF ( SIZE(grid%lwdnt) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%lwdnt(pig,pjg)
ENDIF
IF ( SIZE(grid%lwdntc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%lwdntc(pig,pjg)
ENDIF
IF ( SIZE(grid%lwupb) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%lwupb(pig,pjg)
ENDIF
IF ( SIZE(grid%lwupbc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%lwupbc(pig,pjg)
ENDIF
IF ( SIZE(grid%lwdnb) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%lwdnb(pig,pjg)
ENDIF
IF ( SIZE(grid%lwdnbc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%lwdnbc(pig,pjg)
ENDIF
IF ( SIZE(grid%refl_10cm) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,(ckde-1)
xv(k+ixv+1-ckds)= grid%refl_10cm(pig,k,pjg)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%refd_max) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%refd_max(pig,pjg)
ENDIF
IF ( SIZE(grid%qnwfa2d) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%qnwfa2d(pig,pjg)
ENDIF
IF ( SIZE(grid%swddir) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%swddir(pig,pjg)
ENDIF
IF ( SIZE(grid%swddirc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%swddirc(pig,pjg)
ENDIF
IF ( SIZE(grid%swddni) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%swddni(pig,pjg)
ENDIF
IF ( SIZE(grid%swddnic) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%swddnic(pig,pjg)
ENDIF
IF ( SIZE(grid%swddif) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%swddif(pig,pjg)
ENDIF
IF ( SIZE(grid%gx) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%gx(pig,pjg)
ENDIF
IF ( SIZE(grid%bx) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%bx(pig,pjg)
ENDIF
IF ( SIZE(grid%gg) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%gg(pig,pjg)
ENDIF
IF ( SIZE(grid%bb) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%bb(pig,pjg)
ENDIF
IF ( SIZE(grid%coszen_ref) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%coszen_ref(pig,pjg)
ENDIF
IF ( SIZE(grid%swdown_ref) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%swdown_ref(pig,pjg)
ENDIF
IF ( SIZE(grid%swddir_ref) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%swddir_ref(pig,pjg)
ENDIF
IF ( SIZE(grid%cwm) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,(ckde-1)
xv(k+ixv+1-ckds)= grid%cwm(pig,pjg,k)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%f_ice) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,(ckde-1)
xv(k+ixv+1-ckds)= grid%f_ice(pig,k,pjg)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%f_rain) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,(ckde-1)
xv(k+ixv+1-ckds)= grid%f_rain(pig,k,pjg)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%f_rimef) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,(ckde-1)
xv(k+ixv+1-ckds)= grid%f_rimef(pig,k,pjg)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%cfrach) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%cfrach(pig,pjg)
ENDIF
IF ( SIZE(grid%cfracl) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%cfracl(pig,pjg)
ENDIF
IF ( SIZE(grid%cfracm) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%cfracm(pig,pjg)
ENDIF
IF ( SIZE(grid%islope) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%islope(pig,pjg)
ENDIF
IF ( SIZE(grid%cmc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%cmc(pig,pjg)
ENDIF
IF ( SIZE(grid%soiltb) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%soiltb(pig,pjg)
ENDIF
IF ( SIZE(grid%vegfrc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%vegfrc(pig,pjg)
ENDIF
IF ( SIZE(grid%shdmax) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%shdmax(pig,pjg)
ENDIF
IF ( SIZE(grid%shdmin) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%shdmin(pig,pjg)
ENDIF
IF ( SIZE(grid%sh2o) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,config_flags%num_soil_layers
xv(k+ixv+1-1)= grid%sh2o(pig,k,pjg)
ENDDO
ixv=ixv+config_flags%num_soil_layers-1+1
ENDIF
IF ( SIZE(grid%smc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,config_flags%num_soil_layers
xv(k+ixv+1-1)= grid%smc(pig,k,pjg)
ENDDO
ixv=ixv+config_flags%num_soil_layers-1+1
ENDIF
IF ( SIZE(grid%stc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,config_flags%num_soil_layers
xv(k+ixv+1-1)= grid%stc(pig,k,pjg)
ENDDO
ixv=ixv+config_flags%num_soil_layers-1+1
ENDIF
IF ( SIZE(grid%ctopo) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%ctopo(pig,pjg)
ENDIF
IF ( SIZE(grid%ctopo2) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%ctopo2(pig,pjg)
ENDIF
IF ( SIZE(grid%dwdt) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,(ckde-1)
xv(k+ixv+1-ckds)= grid%dwdt(pig,pjg,k)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%pint) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,ckde
xv(k+ixv+1-ckds)= grid%pint(pig,pjg,k)
ENDDO
ixv=ixv+ckde-ckds+1
ENDIF
IF ( SIZE(grid%w) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,ckde
xv(k+ixv+1-ckds)= grid%w(pig,pjg,k)
ENDDO
ixv=ixv+ckde-ckds+1
ENDIF
IF ( SIZE(grid%w_tot) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,ckde
xv(k+ixv+1-ckds)= grid%w_tot(pig,pjg,k)
ENDDO
ixv=ixv+ckde-ckds+1
ENDIF
IF ( SIZE(grid%z) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,ckde
xv(k+ixv+1-ckds)= grid%z(pig,pjg,k)
ENDDO
ixv=ixv+ckde-ckds+1
ENDIF
IF ( SIZE(grid%rlwin) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%rlwin(pig,pjg)
ENDIF
IF ( SIZE(grid%rswin) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%rswin(pig,pjg)
ENDIF
IF ( SIZE(grid%rlwtt) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,(ckde-1)
xv(k+ixv+1-ckds)= grid%rlwtt(pig,pjg,k)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%rswtt) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,(ckde-1)
xv(k+ixv+1-ckds)= grid%rswtt(pig,pjg,k)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%winfo) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,ckde
xv(k+ixv+1-ckds)= grid%winfo(pig,pjg,k)
ENDDO
ixv=ixv+ckde-ckds+1
ENDIF
IF ( SIZE(grid%iinfo) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,ckde
xv(k+ixv+1-ckds)= grid%iinfo(pig,pjg,k)
ENDDO
ixv=ixv+ckde-ckds+1
ENDIF
IF ( SIZE(grid%landmask) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%landmask(pig,pjg)
ENDIF
IF ( SIZE(grid%toposoil) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%toposoil(pig,pjg)
ENDIF
if(interp_mp .eqv. .true.) then
DO itrace =  PARAM_FIRST_SCALAR, num_moist
DO k = ckds,(ckde-1)
xv(k+ixv+1-ckds)= moist(pig,pjg,k,itrace)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDDO
endif
if(interp_mp .eqv. .true.) then
DO itrace =  PARAM_FIRST_SCALAR, num_scalar
DO k = ckds,(ckde-1)
xv(k+ixv+1-ckds)= scalar(pig,pjg,k,itrace)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDDO
endif
if(interp_mp .eqv. .true.) then
DO itrace =  PARAM_FIRST_SCALAR, num_dfi_scalar
DO k = ckds,(ckde-1)
xv(k+ixv+1-ckds)= dfi_scalar(pig,pjg,k,itrace)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDDO
endif
IF ( SIZE(grid%lake_depth) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%lake_depth(pig,pjg)
ENDIF
IF ( SIZE(grid%u10) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%u10(pig,pjg)
ENDIF
IF ( SIZE(grid%v10) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%v10(pig,pjg)
ENDIF
IF ( SIZE(grid%xice) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%xice(pig,pjg)
ENDIF
IF ( SIZE(grid%icedepth) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%icedepth(pig,pjg)
ENDIF
IF ( SIZE(grid%albsi) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%albsi(pig,pjg)
ENDIF
IF ( SIZE(grid%snowsi) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%snowsi(pig,pjg)
ENDIF
IF ( SIZE(grid%ivgtyp) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%ivgtyp(pig,pjg)
ENDIF
IF ( SIZE(grid%isltyp) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%isltyp(pig,pjg)
ENDIF
IF ( SIZE(grid%vegfra) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%vegfra(pig,pjg)
ENDIF
IF ( SIZE(grid%sst) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%sst(pig,pjg)
ENDIF
IF ( SIZE(grid%weasd) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%weasd(pig,pjg)
ENDIF
IF ( SIZE(grid%thz0) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%thz0(pig,pjg)
ENDIF
IF ( SIZE(grid%qz0) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%qz0(pig,pjg)
ENDIF
IF ( SIZE(grid%uz0) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%uz0(pig,pjg)
ENDIF
IF ( SIZE(grid%vz0) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%vz0(pig,pjg)
ENDIF
IF ( SIZE(grid%htop) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%htop(pig,pjg)
ENDIF
IF ( SIZE(grid%hbot) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%hbot(pig,pjg)
ENDIF
IF ( SIZE(grid%htopr) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%htopr(pig,pjg)
ENDIF
IF ( SIZE(grid%hbotr) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%hbotr(pig,pjg)
ENDIF
IF ( SIZE(grid%cuppt) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%cuppt(pig,pjg)
ENDIF
IF ( SIZE(grid%snowh) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%snowh(pig,pjg)
ENDIF
IF ( SIZE(grid%rhosn) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%rhosn(pig,pjg)
ENDIF
IF ( SIZE(grid%isnowxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%isnowxy(pig,pjg)
ENDIF
IF ( SIZE(grid%tvxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%tvxy(pig,pjg)
ENDIF
IF ( SIZE(grid%tgxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%tgxy(pig,pjg)
ENDIF
IF ( SIZE(grid%canicexy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%canicexy(pig,pjg)
ENDIF
IF ( SIZE(grid%canliqxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%canliqxy(pig,pjg)
ENDIF
IF ( SIZE(grid%eahxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%eahxy(pig,pjg)
ENDIF
IF ( SIZE(grid%tahxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%tahxy(pig,pjg)
ENDIF
IF ( SIZE(grid%cmxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%cmxy(pig,pjg)
ENDIF
IF ( SIZE(grid%chxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%chxy(pig,pjg)
ENDIF
IF ( SIZE(grid%fwetxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%fwetxy(pig,pjg)
ENDIF
IF ( SIZE(grid%sneqvoxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%sneqvoxy(pig,pjg)
ENDIF
IF ( SIZE(grid%alboldxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%alboldxy(pig,pjg)
ENDIF
IF ( SIZE(grid%qsnowxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%qsnowxy(pig,pjg)
ENDIF
IF ( SIZE(grid%wslakexy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%wslakexy(pig,pjg)
ENDIF
IF ( SIZE(grid%zwtxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%zwtxy(pig,pjg)
ENDIF
IF ( SIZE(grid%waxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%waxy(pig,pjg)
ENDIF
IF ( SIZE(grid%wtxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%wtxy(pig,pjg)
ENDIF
IF ( SIZE(grid%tsnoxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,config_flags%num_snow_layers
xv(k+ixv+1-1)= grid%tsnoxy(pig,k,pjg)
ENDDO
ixv=ixv+config_flags%num_snow_layers-1+1
ENDIF
IF ( SIZE(grid%zsnsoxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,config_flags%num_snso_layers
xv(k+ixv+1-1)= grid%zsnsoxy(pig,k,pjg)
ENDDO
ixv=ixv+config_flags%num_snso_layers-1+1
ENDIF
IF ( SIZE(grid%snicexy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,config_flags%num_snow_layers
xv(k+ixv+1-1)= grid%snicexy(pig,k,pjg)
ENDDO
ixv=ixv+config_flags%num_snow_layers-1+1
ENDIF
IF ( SIZE(grid%snliqxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,config_flags%num_snow_layers
xv(k+ixv+1-1)= grid%snliqxy(pig,k,pjg)
ENDDO
ixv=ixv+config_flags%num_snow_layers-1+1
ENDIF
IF ( SIZE(grid%lfmassxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%lfmassxy(pig,pjg)
ENDIF
IF ( SIZE(grid%rtmassxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%rtmassxy(pig,pjg)
ENDIF
IF ( SIZE(grid%stmassxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%stmassxy(pig,pjg)
ENDIF
IF ( SIZE(grid%woodxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%woodxy(pig,pjg)
ENDIF
IF ( SIZE(grid%stblcpxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%stblcpxy(pig,pjg)
ENDIF
IF ( SIZE(grid%fastcpxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%fastcpxy(pig,pjg)
ENDIF
IF ( SIZE(grid%xsaixy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%xsaixy(pig,pjg)
ENDIF
IF ( SIZE(grid%t2mvxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%t2mvxy(pig,pjg)
ENDIF
IF ( SIZE(grid%t2mbxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%t2mbxy(pig,pjg)
ENDIF
IF ( SIZE(grid%q2mvxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%q2mvxy(pig,pjg)
ENDIF
IF ( SIZE(grid%q2mbxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%q2mbxy(pig,pjg)
ENDIF
IF ( SIZE(grid%tradxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%tradxy(pig,pjg)
ENDIF
IF ( SIZE(grid%neexy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%neexy(pig,pjg)
ENDIF
IF ( SIZE(grid%gppxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%gppxy(pig,pjg)
ENDIF
IF ( SIZE(grid%nppxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%nppxy(pig,pjg)
ENDIF
IF ( SIZE(grid%fvegxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%fvegxy(pig,pjg)
ENDIF
IF ( SIZE(grid%qinxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%qinxy(pig,pjg)
ENDIF
IF ( SIZE(grid%runsfxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%runsfxy(pig,pjg)
ENDIF
IF ( SIZE(grid%runsbxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%runsbxy(pig,pjg)
ENDIF
IF ( SIZE(grid%ecanxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%ecanxy(pig,pjg)
ENDIF
IF ( SIZE(grid%edirxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%edirxy(pig,pjg)
ENDIF
IF ( SIZE(grid%etranxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%etranxy(pig,pjg)
ENDIF
IF ( SIZE(grid%fsaxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%fsaxy(pig,pjg)
ENDIF
IF ( SIZE(grid%firaxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%firaxy(pig,pjg)
ENDIF
IF ( SIZE(grid%aparxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%aparxy(pig,pjg)
ENDIF
IF ( SIZE(grid%psnxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%psnxy(pig,pjg)
ENDIF
IF ( SIZE(grid%savxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%savxy(pig,pjg)
ENDIF
IF ( SIZE(grid%sagxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%sagxy(pig,pjg)
ENDIF
IF ( SIZE(grid%rssunxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%rssunxy(pig,pjg)
ENDIF
IF ( SIZE(grid%rsshaxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%rsshaxy(pig,pjg)
ENDIF
IF ( SIZE(grid%bgapxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%bgapxy(pig,pjg)
ENDIF
IF ( SIZE(grid%wgapxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%wgapxy(pig,pjg)
ENDIF
IF ( SIZE(grid%tgvxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%tgvxy(pig,pjg)
ENDIF
IF ( SIZE(grid%tgbxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%tgbxy(pig,pjg)
ENDIF
IF ( SIZE(grid%chvxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%chvxy(pig,pjg)
ENDIF
IF ( SIZE(grid%chbxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%chbxy(pig,pjg)
ENDIF
IF ( SIZE(grid%shgxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%shgxy(pig,pjg)
ENDIF
IF ( SIZE(grid%shcxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%shcxy(pig,pjg)
ENDIF
IF ( SIZE(grid%shbxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%shbxy(pig,pjg)
ENDIF
IF ( SIZE(grid%evgxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%evgxy(pig,pjg)
ENDIF
IF ( SIZE(grid%evbxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%evbxy(pig,pjg)
ENDIF
IF ( SIZE(grid%ghvxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%ghvxy(pig,pjg)
ENDIF
IF ( SIZE(grid%ghbxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%ghbxy(pig,pjg)
ENDIF
IF ( SIZE(grid%irgxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%irgxy(pig,pjg)
ENDIF
IF ( SIZE(grid%ircxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%ircxy(pig,pjg)
ENDIF
IF ( SIZE(grid%irbxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%irbxy(pig,pjg)
ENDIF
IF ( SIZE(grid%trxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%trxy(pig,pjg)
ENDIF
IF ( SIZE(grid%evcxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%evcxy(pig,pjg)
ENDIF
IF ( SIZE(grid%chleafxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%chleafxy(pig,pjg)
ENDIF
IF ( SIZE(grid%chucxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%chucxy(pig,pjg)
ENDIF
IF ( SIZE(grid%chv2xy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%chv2xy(pig,pjg)
ENDIF
IF ( SIZE(grid%chb2xy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%chb2xy(pig,pjg)
ENDIF
IF ( SIZE(grid%chstarxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%chstarxy(pig,pjg)
ENDIF
IF ( SIZE(grid%soilcomp) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,8
xv(k+ixv+1-1)= grid%soilcomp(pig,k,pjg)
ENDDO
ixv=ixv+8-1+1
ENDIF
IF ( SIZE(grid%soilcl1) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%soilcl1(pig,pjg)
ENDIF
IF ( SIZE(grid%soilcl2) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%soilcl2(pig,pjg)
ENDIF
IF ( SIZE(grid%soilcl3) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%soilcl3(pig,pjg)
ENDIF
IF ( SIZE(grid%soilcl4) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%soilcl4(pig,pjg)
ENDIF
IF ( SIZE(grid%grainxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%grainxy(pig,pjg)
ENDIF
IF ( SIZE(grid%gddxy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%gddxy(pig,pjg)
ENDIF
IF ( SIZE(grid%croptype) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,5
xv(k+ixv+1-1)= grid%croptype(pig,k,pjg)
ENDDO
ixv=ixv+5-1+1
ENDIF
IF ( SIZE(grid%planting) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%planting(pig,pjg)
ENDIF
IF ( SIZE(grid%harvest) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%harvest(pig,pjg)
ENDIF
IF ( SIZE(grid%season_gdd) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%season_gdd(pig,pjg)
ENDIF
IF ( SIZE(grid%gecros_state) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,60
xv(k+ixv+1-1)= grid%gecros_state(pig,k,pjg)
ENDDO
ixv=ixv+60-1+1
ENDIF
CALL rsl_lite_to_parent_msg(ixv*RWORDSIZE,xv)
CALL rsl_lite_to_child_info( msize*RWORDSIZE                               &
                        ,cips,cipe,cjps,cjpe                               &
                        ,iids,iide,ijds,ijde                               &
                        ,nids,nide,njds,njde                               &
                        ,pgr , sw                                          &
                        ,nest_task_offsets(ngrid%id)                      &
                        ,nest_pes_x(grid%id)                            &
                        ,nest_pes_y(grid%id)                            &
                        ,nest_pes_x(intermediate_grid%id)                 &
                        ,nest_pes_y(intermediate_grid%id)                 &
                        ,thisdomain_max_halo_width                         &
                        ,icoord,jcoord                                     &
                        ,idim_cd,jdim_cd                                   &
                        ,pig,pjg,retval )
ENDDO
!ENDOFREGISTRYGENERATEDINCLUDE
