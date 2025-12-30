!STARTOFREGISTRYGENERATEDINCLUDE 'inc/nest_forcedown_pack.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
msize = (15 + ((num_szj - PARAM_FIRST_SCALAR + 1)) & 
 + ((num_s1z - PARAM_FIRST_SCALAR + 1)) & 
 + ((num_spz - PARAM_FIRST_SCALAR + 1)) & 
 + ((num_tcs - PARAM_FIRST_SCALAR + 1)) & 
 + ((num_chem - PARAM_FIRST_SCALAR + 1)) & 
 + ((num_ozmixm - PARAM_FIRST_SCALAR + 1)) )* nlev + 15
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
IF ( SIZE(grid%rt_tendf_stoch) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,config_flags%num_stoch_levels
xv(k+ixv+1-1)= grid%rt_tendf_stoch(pig,pjg,k)
ENDDO
ixv=ixv+config_flags%num_stoch_levels-1+1
ENDIF
IF ( SIZE(grid%rand_pert) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,config_flags%num_stoch_levels
xv(k+ixv+1-1)= grid%rand_pert(pig,pjg,k)
ENDDO
ixv=ixv+config_flags%num_stoch_levels-1+1
ENDIF
IF ( SIZE(grid%rstoch) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,config_flags%num_stoch_levels
xv(k+ixv+1-1)= grid%rstoch(pig,pjg,k)
ENDDO
ixv=ixv+config_flags%num_stoch_levels-1+1
ENDIF
IF ( SIZE(grid%pattern_spp_mp) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,config_flags%num_stoch_levels
xv(k+ixv+1-1)= grid%pattern_spp_mp(pig,pjg,k)
ENDDO
ixv=ixv+config_flags%num_stoch_levels-1+1
ENDIF
IF ( SIZE(grid%pattern_isppt_mp) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,config_flags%num_stoch_levels
xv(k+ixv+1-1)= grid%pattern_isppt_mp(pig,pjg,k)
ENDDO
ixv=ixv+config_flags%num_stoch_levels-1+1
ENDIF
IF ( SIZE(grid%pattern_isppt_pbl) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,config_flags%num_stoch_levels
xv(k+ixv+1-1)= grid%pattern_isppt_pbl(pig,pjg,k)
ENDDO
ixv=ixv+config_flags%num_stoch_levels-1+1
ENDIF
IF ( SIZE(grid%pattern_isppt_rad) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,config_flags%num_stoch_levels
xv(k+ixv+1-1)= grid%pattern_isppt_rad(pig,pjg,k)
ENDDO
ixv=ixv+config_flags%num_stoch_levels-1+1
ENDIF
IF ( SIZE(grid%pattern_isppt_conv) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,config_flags%num_stoch_levels
xv(k+ixv+1-1)= grid%pattern_isppt_conv(pig,pjg,k)
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
IF ( SIZE(grid%pdyn_smooth) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%pdyn_smooth(pig,pjg)
ENDIF
IF ( SIZE(grid%pdyn_parent) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%pdyn_parent(pig,pjg)
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
IF ( SIZE(grid%q2) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,(ckde-1)
xv(k+ixv+1-ckds)= grid%q2(pig,pjg,k)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%cwm) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,(ckde-1)
xv(k+ixv+1-ckds)= grid%cwm(pig,pjg,k)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%pint) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,ckde
xv(k+ixv+1-ckds)= grid%pint(pig,pjg,k)
ENDDO
ixv=ixv+ckde-ckds+1
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
IF ( SIZE(grid%sst) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
xv(ixv)=grid%sst(pig,pjg)
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
