!STARTOFREGISTRYGENERATEDINCLUDE 'inc/nest_forcedown_unpack.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL rsl_lite_from_parent_info(pig,pjg,retval)
DO while ( retval .eq. 1 )
CALL rsl_lite_from_parent_msg(-1,xv) ! -1 = all data
ixv=0 ! Reset read/write counter
IF ( SIZE(grid%rt_tendf_stoch) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,config_flags%num_stoch_levels
grid%rt_tendf_stoch(pig,pjg,k) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+config_flags%num_stoch_levels-1+1
ENDIF
IF ( SIZE(grid%rand_pert) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,config_flags%num_stoch_levels
grid%rand_pert(pig,pjg,k) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+config_flags%num_stoch_levels-1+1
ENDIF
IF ( SIZE(grid%rstoch) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,config_flags%num_stoch_levels
grid%rstoch(pig,pjg,k) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+config_flags%num_stoch_levels-1+1
ENDIF
IF ( SIZE(grid%pattern_spp_mp) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,config_flags%num_stoch_levels
grid%pattern_spp_mp(pig,pjg,k) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+config_flags%num_stoch_levels-1+1
ENDIF
IF ( SIZE(grid%pattern_isppt_mp) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,config_flags%num_stoch_levels
grid%pattern_isppt_mp(pig,pjg,k) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+config_flags%num_stoch_levels-1+1
ENDIF
IF ( SIZE(grid%pattern_isppt_pbl) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,config_flags%num_stoch_levels
grid%pattern_isppt_pbl(pig,pjg,k) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+config_flags%num_stoch_levels-1+1
ENDIF
IF ( SIZE(grid%pattern_isppt_rad) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,config_flags%num_stoch_levels
grid%pattern_isppt_rad(pig,pjg,k) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+config_flags%num_stoch_levels-1+1
ENDIF
IF ( SIZE(grid%pattern_isppt_conv) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,config_flags%num_stoch_levels
grid%pattern_isppt_conv(pig,pjg,k) = xv(k-1+ixv+1)
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
IF ( SIZE(grid%pdyn_smooth) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%pdyn_smooth(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%pdyn_parent) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%pdyn_parent(pig,pjg) = xv(ixv)
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
IF ( SIZE(grid%q2) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,(ckde-1)
grid%q2(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%cwm) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,(ckde-1)
grid%cwm(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%pint) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,ckde
grid%pint(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+ckde-ckds+1
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
IF ( SIZE(grid%sst) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
ixv=ixv+1
grid%sst(pig,pjg) = xv(ixv)
ENDIF
CALL rsl_lite_from_parent_info(pig,pjg,retval)
ENDDO
!ENDOFREGISTRYGENERATEDINCLUDE
