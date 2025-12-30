














MODULE module_dm

   USE module_machine
   USE module_wrf_error
   USE module_driver_constants

   USE module_cpl, ONLY : coupler_on, cpl_init

   IMPLICIT NONE
   INCLUDE 'mpif.h'

   INTEGER, PARAMETER :: max_halo_width = 6

   INTEGER :: ips_save, ipe_save, jps_save, jpe_save, itrace
   INTEGER :: lats_to_mic, minx, miny

   INTEGER :: communicator_stack_cursor = 0
   INTEGER :: current_id  = 1
   INTEGER, DIMENSION(max_domains) ::  ntasks_stack, ntasks_y_stack          &
                                     , ntasks_x_stack, mytask_stack          &
                                     , mytask_x_stack, mytask_y_stack        &
                                     , id_stack                            
   INTEGER, DIMENSION(max_domains) ::  ntasks_store, ntasks_y_store          &
                                     , ntasks_x_store, mytask_store          &
                                     , mytask_x_store, mytask_y_store      
   INTEGER ntasks, ntasks_y, ntasks_x, mytask, mytask_x, mytask_y

   INTEGER, DIMENSION(max_domains) :: local_communicator_stack, local_communicator_periodic_stack &
                                     ,local_iocommunicator_stack                                  &
                                     ,local_communicator_x_stack, local_communicator_y_stack
   INTEGER, DIMENSION(max_domains) :: local_communicator_store, local_communicator_periodic_store &
                                     ,local_iocommunicator_store                                  &
                                     ,local_communicator_x_store, local_communicator_y_store

   INTEGER :: mpi_comm_allcompute         = MPI_UNDEFINED
   INTEGER :: local_communicator          = MPI_UNDEFINED
   INTEGER :: local_communicator_periodic = MPI_UNDEFINED
   INTEGER :: local_iocommunicator        = MPI_UNDEFINED
   INTEGER :: local_communicator_x        = MPI_UNDEFINED
   INTEGER :: local_communicator_y        = MPI_UNDEFINED 
   INTEGER :: local_quilt_comm            = MPI_UNDEFINED 
   LOGICAL :: dm_debug_flag = .FALSE.

   INTEGER intercomm_to_mom( max_domains ), intercomm_to_kid( max_nests, max_domains )
   INTEGER mpi_comm_to_mom( max_domains ), mpi_comm_to_kid( max_nests, max_domains )
   INTEGER which_kid(max_domains), nkids(max_domains)
   INTEGER nest_task_offsets(max_domains)
   LOGICAL intercomm_active( max_domains )
   LOGICAL domain_active_this_task( max_domains )

   INTEGER tasks_per_split
   INTEGER comm_start(max_domains)   




   INTEGER nest_pes_x(max_domains)   
   INTEGER nest_pes_y(max_domains)   
   INTEGER comms_i_am_in (max_domains)  
   INTEGER loc_comm(max_domains)
   LOGICAL poll_servers
   INTEGER nio_tasks_per_group(max_domains), nio_groups, num_io_tasks
   NAMELIST /dm_task_split/ tasks_per_split, comm_start, nest_pes_x, nest_pes_y
   NAMELIST /namelist_quilt/ nio_tasks_per_group, nio_groups, poll_servers



   INTERFACE wrf_dm_maxval
     MODULE PROCEDURE wrf_dm_maxval_real , wrf_dm_maxval_integer, wrf_dm_maxval_doubleprecision
   END INTERFACE

   INTERFACE wrf_dm_minval                       
     MODULE PROCEDURE wrf_dm_minval_real , wrf_dm_minval_integer, wrf_dm_minval_doubleprecision
   END INTERFACE

CONTAINS

   SUBROUTINE MPASPECT( P, MINM, MINN, PROCMIN_M, PROCMIN_N )
      IMPLICIT NONE
      INTEGER P, M, N, MINI, MINM, MINN, PROCMIN_M, PROCMIN_N
      MINI = 2*P
      MINM = 1
      MINN = P
      DO M = 1, P
        IF ( MOD( P, M ) .EQ. 0 ) THEN
          N = P / M
          IF ( ABS(M-N) .LT. MINI                &
               .AND. M .GE. PROCMIN_M            &
               .AND. N .GE. PROCMIN_N            &
             ) THEN
            MINI = ABS(M-N)
            MINM = M
            MINN = N
          END IF
        END IF
      END DO
      IF ( MINM .LT. PROCMIN_M .OR. MINN .LT. PROCMIN_N ) THEN
        WRITE( wrf_err_message , * )'MPASPECT: UNABLE TO GENERATE PROCESSOR MESH.  STOPPING.'
        CALL wrf_message ( TRIM ( wrf_err_message ) )
        WRITE( wrf_err_message , * )' PROCMIN_M ', PROCMIN_M
        CALL wrf_message ( TRIM ( wrf_err_message ) )
        WRITE( wrf_err_message , * )' PROCMIN_N ', PROCMIN_N
        CALL wrf_message ( TRIM ( wrf_err_message ) )
        WRITE( wrf_err_message , * )' P         ', P
        CALL wrf_message ( TRIM ( wrf_err_message ) )
        WRITE( wrf_err_message , * )' MINM      ', MINM
        CALL wrf_message ( TRIM ( wrf_err_message ) )
        WRITE( wrf_err_message , * )' MINN      ', MINN
        CALL wrf_message ( TRIM ( wrf_err_message ) )
        CALL wrf_error_fatal3("<stdin>",125,&
'module_dm: mpaspect' )
      END IF
   RETURN
   END SUBROUTINE MPASPECT

   SUBROUTINE compute_mesh( ntasks , ntasks_x, ntasks_y )
     IMPLICIT NONE
     INTEGER, INTENT(IN)  :: ntasks
     INTEGER, INTENT(OUT) :: ntasks_x, ntasks_y
     INTEGER lats_to_mic
     CALL nl_get_nproc_x ( 1, ntasks_x )
     CALL nl_get_nproc_y ( 1, ntasks_y )

     IF ( ntasks_x .GT. 0 .OR. ntasks_y .GT. 0 ) THEN
       
       IF      ( ntasks_x .GT. 0 .AND. ntasks_y .EQ. -1 ) THEN
         ntasks_y = ntasks / ntasks_x
       
       ELSE IF ( ntasks_x .EQ. -1 .AND. ntasks_y .GT. 0 ) THEN
         ntasks_x = ntasks / ntasks_y
       END IF

       IF ( ntasks_x * ntasks_y .NE. ntasks ) THEN
         WRITE( wrf_err_message , * )'WRF_DM_INITIALIZE (RSL_LITE): nproc_x * nproc_y in namelist ne ',ntasks
         CALL wrf_error_fatal3("<stdin>",150,&
wrf_err_message )
       END IF
     ELSE
       
       
       
       CALL mpaspect ( ntasks, ntasks_x, ntasks_y, 1, 1 )
     END IF
     ntasks_store(1) = ntasks
     ntasks_x_store(1) = ntasks_x
     ntasks_y_store(1) = ntasks_y
   END SUBROUTINE compute_mesh

   SUBROUTINE wrf_dm_initialize
      IMPLICIT NONE
      INTEGER :: local_comm_per, local_comm_x, local_comm_y, local_comm2, new_local_comm, group, newgroup, p, p1, ierr,itmp
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ranks
      INTEGER comdup
      INTEGER, DIMENSION(2) :: dims, coords
      LOGICAL, DIMENSION(2) :: isperiodic
      LOGICAL :: reorder_mesh

      CALL instate_communicators_for_domain(1)

      CALL wrf_get_dm_communicator ( new_local_comm )
      dims(1) = nest_pes_y(1)  
      dims(2) = nest_pes_x(1)  
      isperiodic(1) = .true.
      isperiodic(2) = .true.
      CALL mpi_cart_create( new_local_comm, 2, dims, isperiodic, .false., local_comm_per, ierr )
      local_communicator_periodic_store(1) = local_comm_per

      local_communicator_periodic_store = local_comm_per
      local_communicator_periodic = local_comm_per

      CALL nl_set_nproc_x ( 1, ntasks_x )
      CALL nl_set_nproc_y ( 1, ntasks_y )
      WRITE( wrf_err_message , * )'Ntasks in X ',ntasks_x,', ntasks in Y ',ntasks_y
      CALL wrf_message( wrf_err_message )
      RETURN
   END SUBROUTINE wrf_dm_initialize

   SUBROUTINE get_dm_max_halo_width( id, width )
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: id
     INTEGER, INTENT(OUT) :: width
     IF ( id .EQ. 1 ) THEN   
       width = max_halo_width
     ELSE
       width = max_halo_width + 3
     END IF
     RETURN
   END SUBROUTINE get_dm_max_halo_width

   SUBROUTINE patch_domain_rsl_lite( id  , parent, parent_id, &
                                sd1 , ed1 , sp1 , ep1 , sm1 , em1 ,        &
                                sd2 , ed2 , sp2 , ep2 , sm2 , em2 ,        &
                                sd3 , ed3 , sp3 , ep3 , sm3 , em3 ,        &
                                      sp1x , ep1x , sm1x , em1x , &
                                      sp2x , ep2x , sm2x , em2x , &
                                      sp3x , ep3x , sm3x , em3x , &
                                      sp1y , ep1y , sm1y , em1y , &
                                      sp2y , ep2y , sm2y , em2y , &
                                      sp3y , ep3y , sm3y , em3y , &
                                bdx , bdy )

      USE module_domain, ONLY : domain, head_grid, find_grid_by_id

      IMPLICIT NONE
      INTEGER, INTENT(IN)   :: sd1 , ed1 , sd2 , ed2 , sd3 , ed3 , bdx , bdy
      INTEGER, INTENT(OUT)  :: sp1 , ep1 , sp2 , ep2 , sp3 , ep3 , &
                               sm1 , em1 , sm2 , em2 , sm3 , em3
      INTEGER, INTENT(OUT)  :: sp1x , ep1x , sp2x , ep2x , sp3x , ep3x , &
                               sm1x , em1x , sm2x , em2x , sm3x , em3x
      INTEGER, INTENT(OUT)  :: sp1y , ep1y , sp2y , ep2y , sp3y , ep3y , &
                               sm1y , em1y , sm2y , em2y , sm3y , em3y
      INTEGER, INTENT(IN)   :: id, parent_id
      TYPE(domain),POINTER  :: parent


      INTEGER               :: ids, ide, jds, jde, kds, kde
      INTEGER               :: ims, ime, jms, jme, kms, kme
      INTEGER               :: ips, ipe, jps, jpe, kps, kpe
      INTEGER               :: imsx, imex, jmsx, jmex, kmsx, kmex
      INTEGER               :: ipsx, ipex, jpsx, jpex, kpsx, kpex
      INTEGER               :: imsy, imey, jmsy, jmey, kmsy, kmey
      INTEGER               :: ipsy, ipey, jpsy, jpey, kpsy, kpey

      INTEGER               :: c_sd1 , c_ed1 , c_sd2 , c_ed2 , c_sd3 , c_ed3
      INTEGER               :: c_sp1 , c_ep1 , c_sp2 , c_ep2 , c_sp3 , c_ep3 , &
                               c_sm1 , c_em1 , c_sm2 , c_em2 , c_sm3 , c_em3
      INTEGER               :: c_sp1x , c_ep1x , c_sp2x , c_ep2x , c_sp3x , c_ep3x , &
                               c_sm1x , c_em1x , c_sm2x , c_em2x , c_sm3x , c_em3x
      INTEGER               :: c_sp1y , c_ep1y , c_sp2y , c_ep2y , c_sp3y , c_ep3y , &
                               c_sm1y , c_em1y , c_sm2y , c_em2y , c_sm3y , c_em3y

      INTEGER               :: c_ids, c_ide, c_jds, c_jde, c_kds, c_kde
      INTEGER               :: c_ims, c_ime, c_jms, c_jme, c_kms, c_kme
      INTEGER               :: c_ips, c_ipe, c_jps, c_jpe, c_kps, c_kpe

      INTEGER               :: idim , jdim , kdim , rem , a, b
      INTEGER               :: i, j, ni, nj, Px, Py, P

      INTEGER               :: parent_grid_ratio, i_parent_start, j_parent_start
      INTEGER               :: shw
      INTEGER               :: idim_cd, jdim_cd, ierr
      INTEGER               :: max_dom


      TYPE(domain), POINTER :: intermediate_grid
      TYPE(domain), POINTER  :: nest_grid
      CHARACTER*256   :: mess

      INTEGER parent_max_halo_width
      INTEGER thisdomain_max_halo_width
      INTEGER lats_to_mic

     lats_to_mic=0
      IF ( lats_to_mic .GT. 0 ) THEN
        minx = -99  
        miny = lats_to_mic  
      ELSE
        minx = 1   
        miny = 1   
      END IF



      SELECT CASE ( model_data_order )
         
         CASE ( DATA_ORDER_ZXY )
            ids = sd2 ; ide = ed2
            jds = sd3 ; jde = ed3
            kds = sd1 ; kde = ed1
         CASE ( DATA_ORDER_XYZ )
            ids = sd1 ; ide = ed1
            jds = sd2 ; jde = ed2
            kds = sd3 ; kde = ed3
         CASE ( DATA_ORDER_XZY )
            ids = sd1 ; ide = ed1
            jds = sd3 ; jde = ed3
            kds = sd2 ; kde = ed2
         CASE ( DATA_ORDER_YXZ)
            ids = sd2 ; ide = ed2
            jds = sd1 ; jde = ed1
            kds = sd3 ; kde = ed3
      END SELECT

      CALL nl_get_max_dom( 1 , max_dom )

      CALL get_dm_max_halo_width( id , thisdomain_max_halo_width )
      IF ( id .GT. 1 ) THEN
        CALL get_dm_max_halo_width( parent%id , parent_max_halo_width )
      END IF

      CALL compute_memory_dims_rsl_lite ( id, thisdomain_max_halo_width, 0 , bdx, bdy,   &
                   ids,  ide,  jds,  jde,  kds,  kde, &
                   ims,  ime,  jms,  jme,  kms,  kme, &
                   imsx, imex, jmsx, jmex, kmsx, kmex, &
                   imsy, imey, jmsy, jmey, kmsy, kmey, &
                   ips,  ipe,  jps,  jpe,  kps,  kpe, &
                   ipsx, ipex, jpsx, jpex, kpsx, kpex, &
                   ipsy, ipey, jpsy, jpey, kpsy, kpey )

     
     
     
     

      IF ( id .GT. 1 ) THEN
         CALL nl_get_parent_grid_ratio( id, parent_grid_ratio )
         if ( mod(ime,parent_grid_ratio) .NE. 0 ) ime = ime + parent_grid_ratio - mod(ime,parent_grid_ratio)
         if ( mod(jme,parent_grid_ratio) .NE. 0 ) jme = jme + parent_grid_ratio - mod(jme,parent_grid_ratio)
      END IF

      SELECT CASE ( model_data_order )
         CASE ( DATA_ORDER_ZXY )
            sp2 = ips ; ep2 = ipe ; sm2 = ims ; em2 = ime
            sp3 = jps ; ep3 = jpe ; sm3 = jms ; em3 = jme
            sp1 = kps ; ep1 = kpe ; sm1 = kms ; em1 = kme
            sp2x = ipsx ; ep2x = ipex ; sm2x = imsx ; em2x = imex
            sp3x = jpsx ; ep3x = jpex ; sm3x = jmsx ; em3x = jmex
            sp1x = kpsx ; ep1x = kpex ; sm1x = kmsx ; em1x = kmex
            sp2y = ipsy ; ep2y = ipey ; sm2y = imsy ; em2y = imey
            sp3y = jpsy ; ep3y = jpey ; sm3y = jmsy ; em3y = jmey
            sp1y = kpsy ; ep1y = kpey ; sm1y = kmsy ; em1y = kmey
         CASE ( DATA_ORDER_ZYX )
            sp3 = ips ; ep3 = ipe ; sm3 = ims ; em3 = ime
            sp2 = jps ; ep2 = jpe ; sm2 = jms ; em2 = jme
            sp1 = kps ; ep1 = kpe ; sm1 = kms ; em1 = kme
            sp3x = ipsx ; ep3x = ipex ; sm3x = imsx ; em3x = imex
            sp2x = jpsx ; ep2x = jpex ; sm2x = jmsx ; em2x = jmex
            sp1x = kpsx ; ep1x = kpex ; sm1x = kmsx ; em1x = kmex
            sp3y = ipsy ; ep3y = ipey ; sm3y = imsy ; em3y = imey
            sp2y = jpsy ; ep2y = jpey ; sm2y = jmsy ; em2y = jmey
            sp1y = kpsy ; ep1y = kpey ; sm1y = kmsy ; em1y = kmey
         CASE ( DATA_ORDER_XYZ )
            sp1 = ips ; ep1 = ipe ; sm1 = ims ; em1 = ime
            sp2 = jps ; ep2 = jpe ; sm2 = jms ; em2 = jme
            sp3 = kps ; ep3 = kpe ; sm3 = kms ; em3 = kme
            sp1x = ipsx ; ep1x = ipex ; sm1x = imsx ; em1x = imex
            sp2x = jpsx ; ep2x = jpex ; sm2x = jmsx ; em2x = jmex
            sp3x = kpsx ; ep3x = kpex ; sm3x = kmsx ; em3x = kmex
            sp1y = ipsy ; ep1y = ipey ; sm1y = imsy ; em1y = imey
            sp2y = jpsy ; ep2y = jpey ; sm2y = jmsy ; em2y = jmey
            sp3y = kpsy ; ep3y = kpey ; sm3y = kmsy ; em3y = kmey
         CASE ( DATA_ORDER_YXZ)
            sp2 = ips ; ep2 = ipe ; sm2 = ims ; em2 = ime
            sp1 = jps ; ep1 = jpe ; sm1 = jms ; em1 = jme
            sp3 = kps ; ep3 = kpe ; sm3 = kms ; em3 = kme
            sp2x = ipsx ; ep2x = ipex ; sm2x = imsx ; em2x = imex
            sp1x = jpsx ; ep1x = jpex ; sm1x = jmsx ; em1x = jmex
            sp3x = kpsx ; ep3x = kpex ; sm3x = kmsx ; em3x = kmex
            sp2y = ipsy ; ep2y = ipey ; sm2y = imsy ; em2y = imey
            sp1y = jpsy ; ep1y = jpey ; sm1y = jmsy ; em1y = jmey
            sp3y = kpsy ; ep3y = kpey ; sm3y = kmsy ; em3y = kmey
         CASE ( DATA_ORDER_XZY )
            sp1 = ips ; ep1 = ipe ; sm1 = ims ; em1 = ime
            sp3 = jps ; ep3 = jpe ; sm3 = jms ; em3 = jme
            sp2 = kps ; ep2 = kpe ; sm2 = kms ; em2 = kme
            sp1x = ipsx ; ep1x = ipex ; sm1x = imsx ; em1x = imex
            sp3x = jpsx ; ep3x = jpex ; sm3x = jmsx ; em3x = jmex
            sp2x = kpsx ; ep2x = kpex ; sm2x = kmsx ; em2x = kmex
            sp1y = ipsy ; ep1y = ipey ; sm1y = imsy ; em1y = imey
            sp3y = jpsy ; ep3y = jpey ; sm3y = jmsy ; em3y = jmey
            sp2y = kpsy ; ep2y = kpey ; sm2y = kmsy ; em2y = kmey
         CASE ( DATA_ORDER_YZX )
            sp3 = ips ; ep3 = ipe ; sm3 = ims ; em3 = ime
            sp1 = jps ; ep1 = jpe ; sm1 = jms ; em1 = jme
            sp2 = kps ; ep2 = kpe ; sm2 = kms ; em2 = kme
            sp3x = ipsx ; ep3x = ipex ; sm3x = imsx ; em3x = imex
            sp1x = jpsx ; ep1x = jpex ; sm1x = jmsx ; em1x = jmex
            sp2x = kpsx ; ep2x = kpex ; sm2x = kmsx ; em2x = kmex
            sp3y = ipsy ; ep3y = ipey ; sm3y = imsy ; em3y = imey
            sp1y = jpsy ; ep1y = jpey ; sm1y = jmsy ; em1y = jmey
            sp2y = kpsy ; ep2y = kpey ; sm2y = kmsy ; em2y = kmey
      END SELECT

      IF ( id.EQ.1 ) THEN
         WRITE(wrf_err_message,*)'*************************************'
         CALL wrf_message( TRIM(wrf_err_message) )
         WRITE(wrf_err_message,*)'Parent domain'
         CALL wrf_message( TRIM(wrf_err_message) )
         WRITE(wrf_err_message,*)'ids,ide,jds,jde ',ids,ide,jds,jde
         CALL wrf_message( TRIM(wrf_err_message) )
         WRITE(wrf_err_message,*)'ims,ime,jms,jme ',ims,ime,jms,jme
         CALL wrf_message( TRIM(wrf_err_message) )
         WRITE(wrf_err_message,*)'ips,ipe,jps,jpe ',ips,ipe,jps,jpe
         CALL wrf_message( TRIM(wrf_err_message) )
         WRITE(wrf_err_message,*)'*************************************'
         CALL wrf_message( TRIM(wrf_err_message) )
      END IF

      IF ( id .GT. 1 ) THEN

         CALL nl_get_shw( id, shw )
         CALL nl_get_i_parent_start( id , i_parent_start )
         CALL nl_get_j_parent_start( id , j_parent_start )
         CALL nl_get_parent_grid_ratio( id, parent_grid_ratio )

         SELECT CASE ( model_data_order )
            CASE ( DATA_ORDER_ZXY )
               idim = ed2-sd2+1
               jdim = ed3-sd3+1
               kdim = ed1-sd1+1
               c_kds = sd1                ; c_kde = ed1
            CASE ( DATA_ORDER_ZYX )
               idim = ed3-sd3+1
               jdim = ed2-sd2+1
               kdim = ed1-sd1+1
               c_kds = sd1                ; c_kde = ed1
            CASE ( DATA_ORDER_XYZ )
               idim = ed1-sd1+1
               jdim = ed2-sd2+1
               kdim = ed3-sd3+1
               c_kds = sd3                ; c_kde = ed3
            CASE ( DATA_ORDER_YXZ)
               idim = ed2-sd2+1
               jdim = ed1-sd1+1
               kdim = ed3-sd3+1
               c_kds = sd3                ; c_kde = ed3
            CASE ( DATA_ORDER_XZY )
               idim = ed1-sd1+1
               jdim = ed3-sd3+1
               kdim = ed2-sd2+1
               c_kds = sd2                ; c_kde = ed2
            CASE ( DATA_ORDER_YZX )
               idim = ed3-sd3+1
               jdim = ed1-sd1+1
               kdim = ed2-sd2+1
               c_kds = sd2                ; c_kde = ed2
         END SELECT

         idim_cd = idim / parent_grid_ratio + 1 + 2*shw + 1
         jdim_cd = jdim / parent_grid_ratio + 1 + 2*shw + 1

         c_ids = i_parent_start-shw ; c_ide = c_ids + idim_cd - 1
         c_jds = j_parent_start-shw ; c_jde = c_jds + jdim_cd - 1

         
         

         c_ips = -1
         nj = ( c_jds - j_parent_start ) * parent_grid_ratio + 1 + 1 ;
         ierr = 0
         DO i = c_ids, c_ide
            ni = ( i - i_parent_start ) * parent_grid_ratio + 1 + 1 ;

            CALL task_for_point ( ni, nj, ids, ide, jds, jde, nest_pes_x(id), nest_pes_y(id),Px,Py, &
                                  minx, miny,  ierr )
            IF ( ierr .NE. 0 ) CALL wrf_error_fatal3("<stdin>",461,&
'error code returned by task_for_point in module_dm.F (a)')
            IF ( Px .EQ. mytask_x ) THEN
               c_ipe = i
               IF ( c_ips .EQ. -1 ) c_ips = i
            END IF
         END DO
         IF ( ierr .NE. 0 ) THEN
            CALL tfp_message("<stdin>",542)
         END IF
         IF (c_ips .EQ. -1 ) THEN
            c_ipe = -1
            c_ips = 0
         END IF

         c_jps = -1
         ni = ( c_ids - i_parent_start ) * parent_grid_ratio + 1 + 1 ;
         ierr = 0
         DO j = c_jds, c_jde
            nj = ( j - j_parent_start ) * parent_grid_ratio + 1 + 1 ;

            CALL task_for_point ( ni, nj, ids, ide, jds, jde, nest_pes_x(id), nest_pes_y(id), Px, Py, &
                                  minx, miny, ierr )
            IF ( ierr .NE. 0 ) CALL wrf_error_fatal3("<stdin>",484,&
'error code returned by task_for_point in module_dm.F (b)')


            IF ( Py .EQ. mytask_y ) THEN
               c_jpe = j
               IF ( c_jps .EQ. -1 ) c_jps = j
            END IF
         END DO
         IF ( ierr .NE. 0 ) THEN
            CALL tfp_message("<stdin>",566)
         END IF
         IF (c_jps .EQ. -1 ) THEN
            c_jpe = -1
            c_jps = 0
         END IF



         IF ( c_ips <= c_ipe ) THEN

           IF ( mytask_x .EQ. 0 ) THEN
             c_ips = c_ips - shw
           END IF

           IF ( mytask_x .EQ. nest_pes_x(id)-1 ) THEN
             c_ipe = c_ipe + shw
           END IF
           c_ims = max( c_ips - max(shw,thisdomain_max_halo_width), c_ids - bdx ) - 1
           c_ime = min( c_ipe + max(shw,thisdomain_max_halo_width), c_ide + bdx ) + 1
         ELSE
           c_ims = 0
           c_ime = 0
         END IF



         IF ( c_jps <= c_jpe ) THEN

           IF ( mytask_y .EQ. 0 ) THEN
              c_jps = c_jps - shw
           END IF

           IF ( mytask_y .EQ. nest_pes_y(id)-1 ) THEN
              c_jpe = c_jpe + shw
           END IF
           c_jms = max( c_jps - max(shw,thisdomain_max_halo_width), c_jds - bdx ) - 1
           c_jme = min( c_jpe + max(shw,thisdomain_max_halo_width), c_jde + bdx ) + 1

         ELSE
           c_jms = 0
           c_jme = 0
         END IF
         c_kps = 1
         c_kpe = c_kde
         c_kms = 1
         c_kme = c_kde


         c_sm1x = 1 ; c_em1x = 1 ; c_sm2x = 1 ; c_em2x = 1 ; c_sm3x = 1 ; c_em3x = 1
         c_sm1y = 1 ; c_em1y = 1 ; c_sm2y = 1 ; c_em2y = 1 ; c_sm3y = 1 ; c_em3y = 1


         WRITE(wrf_err_message,*)'*************************************'
         CALL wrf_message( TRIM(wrf_err_message) )
         WRITE(wrf_err_message,*)'Nesting domain'
         CALL wrf_message( TRIM(wrf_err_message) )
         WRITE(wrf_err_message,*)'ids,ide,jds,jde ',ids,ide,jds,jde
         CALL wrf_message( TRIM(wrf_err_message) )
         WRITE(wrf_err_message,*)'ims,ime,jms,jme ',ims,ime,jms,jme
         CALL wrf_message( TRIM(wrf_err_message) )
         WRITE(wrf_err_message,*)'ips,ipe,jps,jpe ',ips,ipe,jps,jpe
         CALL wrf_message( TRIM(wrf_err_message) )
         WRITE(wrf_err_message,*)'INTERMEDIATE domain'
         CALL wrf_message( TRIM(wrf_err_message) )
         WRITE(wrf_err_message,*)'ids,ide,jds,jde ',c_ids,c_ide,c_jds,c_jde
         CALL wrf_message( TRIM(wrf_err_message) )
         WRITE(wrf_err_message,*)'ims,ime,jms,jme ',c_ims,c_ime,c_jms,c_jme
         CALL wrf_message( TRIM(wrf_err_message) )
         WRITE(wrf_err_message,*)'ips,ipe,jps,jpe ',c_ips,c_ipe,c_jps,c_jpe
         CALL wrf_message( TRIM(wrf_err_message) )
         WRITE(wrf_err_message,*)'*************************************'
         CALL wrf_message( TRIM(wrf_err_message) )

         SELECT CASE ( model_data_order )
            CASE ( DATA_ORDER_ZXY )
               c_sd2 = c_ids ; c_ed2 = c_ide ; c_sp2 = c_ips ; c_ep2 = c_ipe ; c_sm2 = c_ims ; c_em2 = c_ime
               c_sd3 = c_jds ; c_ed3 = c_jde ; c_sp3 = c_jps ; c_ep3 = c_jpe ; c_sm3 = c_jms ; c_em3 = c_jme
               c_sd1 = c_kds ; c_ed1 = c_kde ; c_sp1 = c_kps ; c_ep1 = c_kpe ; c_sm1 = c_kms ; c_em1 = c_kme
            CASE ( DATA_ORDER_ZYX )
               c_sd3 = c_ids ; c_ed3 = c_ide ; c_sp3 = c_ips ; c_ep3 = c_ipe ; c_sm3 = c_ims ; c_em3 = c_ime
               c_sd2 = c_jds ; c_ed2 = c_jde ; c_sp2 = c_jps ; c_ep2 = c_jpe ; c_sm2 = c_jms ; c_em2 = c_jme
               c_sd1 = c_kds ; c_ed1 = c_kde ; c_sp1 = c_kps ; c_ep1 = c_kpe ; c_sm1 = c_kms ; c_em1 = c_kme
            CASE ( DATA_ORDER_XYZ )
               c_sd1 = c_ids ; c_ed1 = c_ide ; c_sp1 = c_ips ; c_ep1 = c_ipe ; c_sm1 = c_ims ; c_em1 = c_ime
               c_sd2 = c_jds ; c_ed2 = c_jde ; c_sp2 = c_jps ; c_ep2 = c_jpe ; c_sm2 = c_jms ; c_em2 = c_jme
               c_sd3 = c_kds ; c_ed3 = c_kde ; c_sp3 = c_kps ; c_ep3 = c_kpe ; c_sm3 = c_kms ; c_em3 = c_kme
            CASE ( DATA_ORDER_YXZ)
               c_sd2 = c_ids ; c_ed2 = c_ide ; c_sp2 = c_ips ; c_ep2 = c_ipe ; c_sm2 = c_ims ; c_em2 = c_ime
               c_sd1 = c_jds ; c_ed1 = c_jde ; c_sp1 = c_jps ; c_ep1 = c_jpe ; c_sm1 = c_jms ; c_em1 = c_jme
               c_sd3 = c_kds ; c_ed3 = c_kde ; c_sp3 = c_kps ; c_ep3 = c_kpe ; c_sm3 = c_kms ; c_em3 = c_kme
            CASE ( DATA_ORDER_XZY )
               c_sd1 = c_ids ; c_ed1 = c_ide ; c_sp1 = c_ips ; c_ep1 = c_ipe ; c_sm1 = c_ims ; c_em1 = c_ime
               c_sd3 = c_jds ; c_ed3 = c_jde ; c_sp3 = c_jps ; c_ep3 = c_jpe ; c_sm3 = c_jms ; c_em3 = c_jme
               c_sd2 = c_kds ; c_ed2 = c_kde ; c_sp2 = c_kps ; c_ep2 = c_kpe ; c_sm2 = c_kms ; c_em2 = c_kme
            CASE ( DATA_ORDER_YZX )
               c_sd3 = c_ids ; c_ed3 = c_ide ; c_sp3 = c_ips ; c_ep3 = c_ipe ; c_sm3 = c_ims ; c_em3 = c_ime
               c_sd1 = c_jds ; c_ed1 = c_jde ; c_sp1 = c_jps ; c_ep1 = c_jpe ; c_sm1 = c_jms ; c_em1 = c_jme
               c_sd2 = c_kds ; c_ed2 = c_kde ; c_sp2 = c_kps ; c_ep2 = c_kpe ; c_sm2 = c_kms ; c_em2 = c_kme
         END SELECT

         ALLOCATE ( intermediate_grid )
         ALLOCATE ( intermediate_grid%parents( max_parents ) )
         ALLOCATE ( intermediate_grid%nests( max_nests ) )
         intermediate_grid%allocated=.false.
         NULLIFY( intermediate_grid%sibling )
         DO i = 1, max_nests
            NULLIFY( intermediate_grid%nests(i)%ptr )
         END DO
         NULLIFY  (intermediate_grid%next)
         NULLIFY  (intermediate_grid%same_level)
         NULLIFY  (intermediate_grid%i_start)
         NULLIFY  (intermediate_grid%j_start)
         NULLIFY  (intermediate_grid%i_end)
         NULLIFY  (intermediate_grid%j_end)
         intermediate_grid%id = id   
         intermediate_grid%num_nests = 0
         intermediate_grid%num_siblings = 0
         intermediate_grid%num_parents = 1
         intermediate_grid%max_tiles   = 0
         intermediate_grid%num_tiles_spec   = 0
         CALL find_grid_by_id ( id, head_grid, nest_grid )

         nest_grid%intermediate_grid => intermediate_grid  
         intermediate_grid%parents(1)%ptr => nest_grid     
         intermediate_grid%num_parents = 1

         intermediate_grid%is_intermediate = .TRUE.
         SELECT CASE ( model_data_order )
            CASE ( DATA_ORDER_ZXY )
               intermediate_grid%nids = nest_grid%sd32 ; intermediate_grid%njds = nest_grid%sd33
               intermediate_grid%nide = nest_grid%ed32 ; intermediate_grid%njde = nest_grid%sd33
            CASE ( DATA_ORDER_ZYX )
               intermediate_grid%nids = nest_grid%sd33 ; intermediate_grid%njds = nest_grid%sd32
               intermediate_grid%nide = nest_grid%ed33 ; intermediate_grid%njde = nest_grid%sd32
            CASE ( DATA_ORDER_XYZ )
               intermediate_grid%nids = nest_grid%sd31 ; intermediate_grid%njds = nest_grid%sd32
               intermediate_grid%nide = nest_grid%ed31 ; intermediate_grid%njde = nest_grid%sd32
            CASE ( DATA_ORDER_YXZ)
               intermediate_grid%nids = nest_grid%sd32 ; intermediate_grid%njds = nest_grid%sd31
               intermediate_grid%nide = nest_grid%ed32 ; intermediate_grid%njde = nest_grid%sd31
            CASE ( DATA_ORDER_XZY )
               intermediate_grid%nids = nest_grid%sd31 ; intermediate_grid%njds = nest_grid%sd33
               intermediate_grid%nide = nest_grid%ed31 ; intermediate_grid%njde = nest_grid%sd33
            CASE ( DATA_ORDER_YZX )
               intermediate_grid%nids = nest_grid%sd33 ; intermediate_grid%njds = nest_grid%sd31
               intermediate_grid%nide = nest_grid%ed33 ; intermediate_grid%njde = nest_grid%sd31
         END SELECT
         intermediate_grid%nids = ids
         intermediate_grid%nide = ide
         intermediate_grid%njds = jds
         intermediate_grid%njde = jde

         intermediate_grid%sm31x                           = c_sm1x
         intermediate_grid%em31x                           = c_em1x
         intermediate_grid%sm32x                           = c_sm2x
         intermediate_grid%em32x                           = c_em2x
         intermediate_grid%sm33x                           = c_sm3x
         intermediate_grid%em33x                           = c_em3x
         intermediate_grid%sm31y                           = c_sm1y
         intermediate_grid%em31y                           = c_em1y
         intermediate_grid%sm32y                           = c_sm2y
         intermediate_grid%em32y                           = c_em2y
         intermediate_grid%sm33y                           = c_sm3y
         intermediate_grid%em33y                           = c_em3y


         intermediate_grid%sd31                            =   c_sd1
         intermediate_grid%ed31                            =   c_ed1
         intermediate_grid%sp31                            = c_sp1
         intermediate_grid%ep31                            = c_ep1
         intermediate_grid%sm31                            = c_sm1
         intermediate_grid%em31                            = c_em1
         intermediate_grid%sd32                            =   c_sd2
         intermediate_grid%ed32                            =   c_ed2
         intermediate_grid%sp32                            = c_sp2
         intermediate_grid%ep32                            = c_ep2
         intermediate_grid%sm32                            = c_sm2
         intermediate_grid%em32                            = c_em2
         intermediate_grid%sd33                            =   c_sd3
         intermediate_grid%ed33                            =   c_ed3
         intermediate_grid%sp33                            = c_sp3
         intermediate_grid%ep33                            = c_ep3
         intermediate_grid%sm33                            = c_sm3
         intermediate_grid%em33                            = c_em3

         CALL med_add_config_info_to_grid ( intermediate_grid )

         intermediate_grid%dx = parent%dx
         intermediate_grid%dy = parent%dy
         intermediate_grid%dt = parent%dt
      END IF

      RETURN
  END SUBROUTINE patch_domain_rsl_lite

  SUBROUTINE compute_memory_dims_rsl_lite  (      &
                   id , maxhalowidth ,            &
                   shw , bdx,  bdy ,              &
                   ids,  ide,  jds,  jde,  kds,  kde, &
                   ims,  ime,  jms,  jme,  kms,  kme, &
                   imsx, imex, jmsx, jmex, kmsx, kmex, &
                   imsy, imey, jmsy, jmey, kmsy, kmey, &
                   ips,  ipe,  jps,  jpe,  kps,  kpe, &
                   ipsx, ipex, jpsx, jpex, kpsx, kpex, &
                   ipsy, ipey, jpsy, jpey, kpsy, kpey )

    IMPLICIT NONE
    INTEGER, INTENT(IN)               ::  id , maxhalowidth
    INTEGER, INTENT(IN)               ::  shw, bdx, bdy
    INTEGER, INTENT(IN)     ::  ids, ide, jds, jde, kds, kde
    INTEGER, INTENT(OUT)    ::  ims, ime, jms, jme, kms, kme
    INTEGER, INTENT(OUT)    ::  imsx, imex, jmsx, jmex, kmsx, kmex
    INTEGER, INTENT(OUT)    ::  imsy, imey, jmsy, jmey, kmsy, kmey
    INTEGER, INTENT(OUT)    ::  ips, ipe, jps, jpe, kps, kpe
    INTEGER, INTENT(OUT)    ::  ipsx, ipex, jpsx, jpex, kpsx, kpex
    INTEGER, INTENT(OUT)    ::  ipsy, ipey, jpsy, jpey, kpsy, kpey

    INTEGER Px, Py, P, i, j, k, ierr




    ips = -1
    j = jds
    ierr = 0
    DO i = ids, ide-1

       CALL task_for_point ( i, j, ids, ide-1, jds, jde-1, nest_pes_x(id), nest_pes_y(id), Px, Py, &
                             minx, miny, ierr )
       IF ( Px .EQ. mytask_x ) THEN
          ipe = i

          IF ( Px .EQ. nest_pes_x(id)-1 ) ipe = ipe + 1
          IF ( ips .EQ. -1 ) ips = i
       END IF
    END DO
    IF ( ierr .NE. 0 ) THEN
       CALL tfp_message("<stdin>",1178)
    END IF
    jps = -1
    i = ids ;
    ierr = 0
    DO j = jds, jde-1

       CALL task_for_point ( i, j, ids, ide-1, jds, jde-1, nest_pes_x(id), nest_pes_y(id), Px, Py, &
                             minx, miny, ierr )
       IF ( Py .EQ. mytask_y ) THEN
          jpe = j

          IF ( Py .EQ. nest_pes_y(id)-1 ) jpe = jpe + 1
          IF ( jps .EQ. -1 ) jps = j
       END IF
    END DO
    IF ( ierr .NE. 0 ) THEN
       CALL tfp_message("<stdin>",1195)
    END IF




    kpsx = -1
    j = jds ;
    ierr = 0
    DO k = kds, kde

       CALL task_for_point ( k, j, kds, kde, jds, jde-1, nest_pes_x(id), nest_pes_y(id), Px, Py, &
                             minx, miny, ierr )
       IF ( ierr .NE. 0 ) CALL wrf_error_fatal3("<stdin>",762,&
'error code returned by task_for_point in module_dm.F (e)')
       IF ( Px .EQ. mytask_x ) THEN
          kpex = k
          IF ( kpsx .EQ. -1 ) kpsx = k
       END IF
    END DO
    IF ( ierr .NE. 0 ) THEN
       CALL tfp_message("<stdin>",1215)
    END IF



    IF (kpsx .EQ. -1 ) THEN
       kpex = -1
       kpsx = 0
    END IF
   jpsx = -1
    k = kds ;
    ierr = 0
    DO j = jds, jde-1

       CALL task_for_point ( k, j, kds, kde, jds, jde-1, nest_pes_x(id), nest_pes_y(id), Px, Py, &
                             minx, miny, ierr )
       IF ( ierr .NE. 0 ) CALL wrf_error_fatal3("<stdin>",786,&
'error code returned by task_for_point in module_dm.F (f)')
       IF ( Py .EQ. mytask_y ) THEN
          jpex = j
          IF ( Py .EQ. nest_pes_y(id)-1 ) jpex = jpex + 1 
          IF ( jpsx .EQ. -1 ) jpsx = j
       END IF
    END DO
    IF ( ierr .NE. 0 ) THEN
       CALL tfp_message("<stdin>",1239)
    END IF
    IF (jpsx .EQ. -1 ) THEN
       jpex = -1
       jpsx = 0
    END IF




    IF (jpex .EQ. -1) THEN
       ipex = -1
       ipsx = 0
       jpex = -1
       jpsx = 0
    END IF



    kpsy = kpsx   
    kpey = kpex   

    ipsy = -1
    k = kds ;
    ierr = 0
    DO i = ids, ide-1

       CALL task_for_point ( i, k, ids, ide-1, kds, kde, nest_pes_y(id), nest_pes_x(id), Py, Px, &
                             miny, minx, ierr )
       IF ( ierr .NE. 0 ) CALL wrf_error_fatal3("<stdin>",824,&
'error code returned by task_for_point in module_dm.F (g)')
       IF ( Py .EQ. mytask_y ) THEN
          ipey = i
          IF ( Py .EQ. nest_pes_y(id)-1 ) ipey = ipey + 1 
          IF ( ipsy .EQ. -1 ) ipsy = i
       END IF
    END DO
    IF ( ierr .NE. 0 ) THEN
       CALL tfp_message("<stdin>",1276)
    END IF
    IF (ipsy .EQ. -1 ) THEN
       ipey = -1
       ipsy = 0
    END IF




    IF ( ips < ipe .and. jps < jpe ) THEN           
       IF ( mytask_x .EQ. 0 ) THEN
          ips = ips - shw
          ipsy = ipsy - shw
       END IF

       IF ( mytask_x .EQ. nest_pes_x(id)-1 ) THEN
          ipe = ipe + shw
          ipey = ipey + shw
       END IF
       IF ( mytask_y .EQ. 0 ) THEN
          jps = jps - shw
          jpsx = jpsx - shw
       END IF

       IF ( mytask_y .EQ. nest_pes_y(id)-1 ) THEN
          jpe = jpe + shw
          jpex = jpex + shw
       END IF
    END IF                                           

    kps = 1
    kpe = kde-kds+1

    kms = 1
    kme = kpe
    kmsx = kpsx
    kmex = kpex
    kmsy = kpsy
    kmey = kpey

    
    IF ( kpsx .EQ. 0 .AND. kpex .EQ. -1 ) THEN
      kmsx = 0
      kmex = 0
    END IF
    IF ( kpsy .EQ. 0 .AND. kpey .EQ. -1 ) THEN
      kmsy = 0
      kmey = 0
    END IF

    IF ( (jps .EQ. 0 .AND. jpe .EQ. -1) .OR. (ips .EQ. 0 .AND. ipe .EQ. -1) ) THEN
      ims = 0
      ime = 0
    ELSE
      ims = max( ips - max(shw,maxhalowidth), ids - bdx ) - 1
      ime = min( ipe + max(shw,maxhalowidth), ide + bdx ) + 1
    END IF
    imsx = ids
    imex = ide
    ipsx = imsx
    ipex = imex
    
    IF ( ipsy .EQ. 0 .AND. ipey .EQ. -1 ) THEN
      imsy = 0
      imey = 0
    ELSE
      imsy = ipsy
      imey = ipey
    END IF

    IF ( (jps .EQ. 0 .AND. jpe .EQ. -1) .OR. (ips .EQ. 0 .AND. ipe .EQ. -1) ) THEN
      jms = 0
      jme = 0
    ELSE
      jms = max( jps - max(shw,maxhalowidth), jds - bdy ) - 1
      jme = min( jpe + max(shw,maxhalowidth), jde + bdy ) + 1
    END IF
    jmsx = jpsx
    jmex = jpex
    jmsy = jds
    jmey = jde
    
    IF ( jpsx .EQ. 0 .AND. jpex .EQ. -1 ) THEN
      jmsx = 0
      jmex = 0
      jpsy = 0
      jpey = -1
    ELSE
      jpsy = jmsy
      jpey = jmey
    END IF

  END SUBROUTINE compute_memory_dims_rsl_lite



   INTEGER function getrealmpitype()
      IMPLICIT NONE
      INTEGER rtypesize, dtypesize, ierr
      CALL mpi_type_size ( MPI_REAL, rtypesize, ierr )
      CALL mpi_type_size ( MPI_DOUBLE_PRECISION, dtypesize, ierr )
      IF ( 4 .EQ. rtypesize ) THEN
        getrealmpitype = MPI_REAL
      ELSE IF ( 4 .EQ. dtypesize ) THEN
        getrealmpitype = MPI_DOUBLE_PRECISION
      ELSE
        CALL wrf_error_fatal3("<stdin>",940,&
'RWORDSIZE or DWORDSIZE does not match any MPI type' )
      END IF
      RETURN
   END FUNCTION getrealmpitype

   REAL FUNCTION wrf_dm_max_int ( inval )
      IMPLICIT NONE
      INTEGER, intent(in) :: inval
      INTEGER :: ierr, retval
      CALL mpi_allreduce ( inval, retval , 1, MPI_INTEGER, MPI_MAX, local_communicator, ierr )
      wrf_dm_max_int = retval
   END FUNCTION wrf_dm_max_int

   REAL FUNCTION wrf_dm_max_real ( inval )
      IMPLICIT NONE
      REAL inval, retval
      INTEGER comm,ierr
      CALL wrf_get_dm_communicator(comm)
      CALL mpi_allreduce ( inval, retval , 1, getrealmpitype(), MPI_MAX, comm, ierr )
      wrf_dm_max_real = retval
   END FUNCTION wrf_dm_max_real

   REAL FUNCTION wrf_dm_min_real ( inval )
      IMPLICIT NONE
      REAL inval, retval
      INTEGER comm,ierr
      CALL wrf_get_dm_communicator(comm)
      CALL mpi_allreduce ( inval, retval , 1, getrealmpitype(), MPI_MIN, comm, ierr )
      wrf_dm_min_real = retval
   END FUNCTION wrf_dm_min_real

   SUBROUTINE wrf_dm_min_reals ( inval, retval, n )
      IMPLICIT NONE
      INTEGER n
      REAL inval(*)
      REAL retval(*)
      INTEGER comm,ierr
      CALL wrf_get_dm_communicator(comm)
      CALL mpi_allreduce ( inval, retval , n, getrealmpitype(), MPI_MIN, comm, ierr )
   END SUBROUTINE wrf_dm_min_reals

   FUNCTION wrf_dm_sum_real8 ( inval )
     
     
      IMPLICIT NONE
      REAL*8 inval, retval, wrf_dm_sum_real8
      INTEGER comm,ierr
      CALL wrf_get_dm_communicator(comm)
      CALL mpi_allreduce ( inval, retval , 1, MPI_REAL8, MPI_SUM, comm, ierr )
      wrf_dm_sum_real8 = retval
   END FUNCTION wrf_dm_sum_real8

   REAL FUNCTION wrf_dm_sum_real ( inval )
      IMPLICIT NONE
      REAL inval, retval
      INTEGER comm,ierr
      CALL wrf_get_dm_communicator(comm)
      CALL mpi_allreduce ( inval, retval , 1, getrealmpitype(), MPI_SUM, comm, ierr )
      wrf_dm_sum_real = retval
   END FUNCTION wrf_dm_sum_real

   SUBROUTINE wrf_dm_sum_reals (inval, retval)
      IMPLICIT NONE
      REAL, INTENT(IN)  :: inval(:)
      REAL, INTENT(OUT) :: retval(:)
      INTEGER comm,ierr
      CALL wrf_get_dm_communicator(comm)
      CALL mpi_allreduce ( inval, retval, SIZE(inval), getrealmpitype(), MPI_SUM, comm, ierr )
   END SUBROUTINE wrf_dm_sum_reals

   INTEGER FUNCTION wrf_dm_sum_integer ( inval )
      IMPLICIT NONE
      INTEGER inval, retval
      INTEGER comm,ierr
      CALL wrf_get_dm_communicator(comm)
      CALL mpi_allreduce ( inval, retval , 1, MPI_INTEGER, MPI_SUM, comm, ierr )
      wrf_dm_sum_integer = retval
   END FUNCTION wrf_dm_sum_integer

   SUBROUTINE wrf_dm_sum_integers (inval, retval)
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: inval(:)
      INTEGER, INTENT(OUT) :: retval(:)
      INTEGER comm,ierr
      CALL wrf_get_dm_communicator(comm)
      CALL mpi_allreduce ( inval, retval, SIZE(inval), MPI_INTEGER, MPI_SUM, comm, ierr )
   END SUBROUTINE wrf_dm_sum_integers

   SUBROUTINE wrf_dm_minloc_real ( val, lat, lon, z, idex, jdex )
      use mpi
      IMPLICIT NONE
      REAL val, lat, lon, z
      INTEGER idex, jdex, ierr, mrank, comm
      REAL inreduce(2), outreduce(2), bcast(5)

      inreduce=(/ val, real(mytask) /)
      CALL wrf_get_dm_communicator(comm)
      call MPI_Allreduce(inreduce,outreduce,1,MPI_2REAL,MPI_MINLOC,&
           comm,ierr)
      val=outreduce(1)
      mrank=outreduce(2)
      bcast=(/ lat,lon,z,real(idex),real(jdex) /)
      call MPI_Bcast(bcast,5,MPI_REAL,mrank,comm,ierr)
      lat=bcast(1)
      lon=bcast(2)
      z=bcast(3)
      idex=bcast(4)
      jdex=bcast(5)
   END SUBROUTINE wrf_dm_minloc_real
   SUBROUTINE wrf_dm_maxloc_real ( val, lat, lon, z, idex, jdex )
      use mpi
      IMPLICIT NONE
      REAL val, lat, lon, z
      INTEGER idex, jdex, ierr, mrank, comm
      REAL inreduce(2), outreduce(2), bcast(5)

      inreduce=(/ val, real(mytask) /)
      CALL wrf_get_dm_communicator(comm)
      call MPI_Allreduce(inreduce,outreduce,1,MPI_2REAL,MPI_MAXLOC,&
           comm,ierr)
      val=outreduce(1)
      mrank=outreduce(2)
      bcast=(/ lat,lon,z,real(idex),real(jdex) /)
      call MPI_Bcast(bcast,5,MPI_REAL,mrank,comm,ierr)
      lat=bcast(1)
      lon=bcast(2)
      z=bcast(3)
      idex=bcast(4)
      jdex=bcast(5)
   END SUBROUTINE wrf_dm_maxloc_real

   INTEGER FUNCTION wrf_dm_bxor_integer ( inval )
      IMPLICIT NONE
      INTEGER inval, retval
      INTEGER comm, ierr
      CALL wrf_get_dm_communicator(comm)
      CALL mpi_allreduce ( inval, retval , 1, MPI_INTEGER, MPI_BXOR, comm, ierr )
      wrf_dm_bxor_integer = retval
   END FUNCTION wrf_dm_bxor_integer


LOGICAL FUNCTION wrf_dm_lor_logical ( inval )
      IMPLICIT NONE
      LOGICAL inval, retval
      INTEGER comm, ierr
      CALL wrf_get_dm_communicator(comm)
      CALL mpi_allreduce ( inval, retval , 1, MPI_LOGICAL, MPI_LOR, comm, ierr )
      wrf_dm_lor_logical = retval
   END FUNCTION wrf_dm_lor_logical


LOGICAL FUNCTION wrf_dm_land_logical ( inval )
      IMPLICIT NONE
      LOGICAL inval, retval
      INTEGER comm, ierr
      CALL wrf_get_dm_communicator(comm)
      CALL mpi_allreduce ( inval, retval , 1, MPI_LOGICAL, MPI_LAND, comm, ierr )
      wrf_dm_land_logical = retval
   END FUNCTION wrf_dm_land_logical


   SUBROUTINE wrf_dm_maxval_real ( val, idex, jdex )
      use mpi
      IMPLICIT NONE
      REAL val
      INTEGER :: idex, jdex, i, comm
      INTEGER :: bcast(2),mrank
      REAL :: inreduce(2),outreduce(2)

      inreduce=(/ val, real(mytask) /)
      bcast=(/ idex,jdex /)
      CALL wrf_get_dm_communicator(comm)
      call MPI_Allreduce(inreduce,outreduce,1,MPI_2REAL,&
                         MPI_MAXLOC,comm,i)
      mrank=outreduce(2)
      val=outreduce(1)
      call MPI_Bcast(bcast,2,MPI_REAL,mrank,comm,i)
      idex=bcast(1)
      jdex=bcast(2)
    END SUBROUTINE wrf_dm_maxval_real

   SUBROUTINE wrf_dm_minval_real ( val, idex, jdex )
      use mpi
      IMPLICIT NONE
      REAL val
      INTEGER :: idex, jdex, i, comm
      INTEGER :: bcast(2),mrank
      REAL :: inreduce(2),outreduce(2)

      inreduce=(/ val, real(mytask) /)
      bcast=(/ idex,jdex /)
      CALL wrf_get_dm_communicator(comm)
      call MPI_Allreduce(inreduce,outreduce,1,MPI_2REAL,&
                         MPI_MINLOC,comm,i)
      mrank=outreduce(2)
      val=outreduce(1)
      call MPI_Bcast(bcast,2,MPI_REAL,mrank,comm,i)
      idex=bcast(1)
      jdex=bcast(2)
    END SUBROUTINE wrf_dm_minval_real

   SUBROUTINE wrf_dm_maxval_doubleprecision ( val, idex, jdex )
      use mpi
      IMPLICIT NONE
      DOUBLE PRECISION val
      INTEGER :: idex, jdex, i, comm
      INTEGER :: bcast(2),mrank
      DOUBLE PRECISION :: inreduce(2),outreduce(2)

      inreduce=(/ val, dble(mytask) /)
      bcast=(/ idex,jdex /)
      CALL wrf_get_dm_communicator(comm)
      call MPI_Allreduce(inreduce,outreduce,1,MPI_2DOUBLE_PRECISION,&
                         MPI_MAXLOC,comm,i)
      mrank=outreduce(2)
      val=outreduce(1)
      call MPI_Bcast(bcast,2,MPI_DOUBLE_PRECISION,mrank,comm,i)
      idex=bcast(1)
      jdex=bcast(2)
   END SUBROUTINE wrf_dm_maxval_doubleprecision

   SUBROUTINE wrf_dm_minval_doubleprecision ( val, idex, jdex )
      use mpi
      IMPLICIT NONE
      DOUBLE PRECISION val
      INTEGER :: idex, jdex, i, comm
      INTEGER :: bcast(2),mrank
      DOUBLE PRECISION :: inreduce(2),outreduce(2)

      inreduce=(/ val, dble(mytask) /)
      bcast=(/ idex,jdex /)
      CALL wrf_get_dm_communicator(comm)
      call MPI_Allreduce(inreduce,outreduce,1,MPI_2DOUBLE_PRECISION,&
                         MPI_MINLOC,comm,i)
      mrank=outreduce(2)
      val=outreduce(1)
      call MPI_Bcast(bcast,2,MPI_DOUBLE_PRECISION,mrank,comm,i)
      idex=bcast(1)
      jdex=bcast(2)
   END SUBROUTINE wrf_dm_minval_doubleprecision

   SUBROUTINE wrf_dm_maxval_integer ( val, idex, jdex )
      use mpi
      IMPLICIT NONE
      INTEGER val
      INTEGER :: idex, jdex, i, comm
      INTEGER :: bcast(2),mrank
      INTEGER :: inreduce(2),outreduce(2)

      inreduce=(/ val, mytask /)
      bcast=(/ idex,jdex /)
      CALL wrf_get_dm_communicator(comm)
      call MPI_Allreduce(inreduce,outreduce,1,MPI_2INTEGER,&
                         MPI_MAXLOC,comm,i)
      mrank=outreduce(2)
      val=outreduce(1)
      call MPI_Bcast(bcast,2,MPI_INTEGER,mrank,comm,i)
      idex=bcast(1)
      jdex=bcast(2)
    END SUBROUTINE wrf_dm_maxval_integer

   SUBROUTINE wrf_dm_minval_integer ( val, idex, jdex )
      use mpi
      IMPLICIT NONE
      INTEGER val
      INTEGER :: idex, jdex, i, comm
      INTEGER :: bcast(2),mrank
      INTEGER :: inreduce(2),outreduce(2)

      inreduce=(/ val, mytask /)
      bcast=(/ idex,jdex /)
      CALL wrf_get_dm_communicator(comm)
      call MPI_Allreduce(inreduce,outreduce,1,MPI_2INTEGER,&
                         MPI_MINLOC,comm,i)
      mrank=outreduce(2)
      val=outreduce(1)
      call MPI_Bcast(bcast,2,MPI_INTEGER,mrank,comm,i)
      idex=bcast(1)
      jdex=bcast(2)
    END SUBROUTINE wrf_dm_minval_integer

   SUBROUTINE hwrf_coupler_init
      IMPLICIT NONE
      LOGICAL mpi_inited
      INTEGER mpi_comm_here,ierr
      CALL MPI_INITIALIZED( mpi_inited, ierr )
      IF ( .NOT. mpi_inited ) THEN
        IF ( coupler_on ) THEN
           CALL cpl_init( mpi_comm_here )
        ELSE
           CALL mpi_init ( ierr )
           mpi_comm_here = MPI_COMM_WORLD
        END IF
        CALL atm_cmp_start( mpi_comm_here )
        CALL wrf_set_dm_communicator( mpi_comm_here )
        CALL wrf_termio_dup( mpi_comm_here )
      END IF
      RETURN
   END SUBROUTINE hwrf_coupler_init

   SUBROUTINE split_communicator
      IMPLICIT NONE
      LOGICAL mpi_inited

      INTEGER mpi_comm_here, mpi_comm_local, comdup, comdup2, origmytask,  ierr, io_status
      INTEGER mpi_comm_me_and_mom
      INTEGER coords(3)
      INTEGER mytask_local,ntasks_local,num_compute_tasks
      INTEGER i, j, k, x, y, n_x, n_y
      INTEGER iii
      INTEGER, ALLOCATABLE :: icolor(:),icolor2(:),idomain(:)
      INTEGER comm_id








































      INTEGER dims(3)

      INTEGER :: id
      INTEGER :: intercomm
      INTEGER :: domain_id,par_id,nest_id,kid_id
      INTEGER :: mytask_me_and_mom, ntasks_me_and_mom, remote_leader
      LOGICAL :: inthisone
      LOGICAL :: mytask_is_nest, mytask_is_par,isperiodic(3)

      LOGICAL :: quilting_is_turned_off









integer    :: first_item_in_struct
real , DIMENSION(max_domains) :: lakedepth_default
real , DIMENSION(max_domains) :: lake_min_elev
integer , DIMENSION(max_domains) :: use_lakedepth
integer :: nens
integer , DIMENSION(max_domains) :: skebs
integer , DIMENSION(max_domains) :: stoch_force_opt
integer :: skebs_vertstruc
integer , DIMENSION(max_domains) :: stoch_vertstruc_opt
real , DIMENSION(max_domains) :: tot_backscat_psi
real , DIMENSION(max_domains) :: tot_backscat_t
real :: ztau_psi
real :: ztau_t
real :: rexponent_psi
real :: rexponent_t
real :: zsigma2_eps
real :: zsigma2_eta
integer :: kminforc
integer :: lminforc
integer :: kminforct
integer :: lminforct
integer :: kmaxforc
integer :: lmaxforc
integer :: kmaxforct
integer :: lmaxforct
integer :: iseed_skebs
integer :: kmaxforch
integer :: lmaxforch
integer :: kmaxforcth
integer :: lmaxforcth
integer , DIMENSION(max_domains) :: sppt
real , DIMENSION(max_domains) :: gridpt_stddev_sppt
real , DIMENSION(max_domains) :: stddev_cutoff_sppt
real , DIMENSION(max_domains) :: lengthscale_sppt
real , DIMENSION(max_domains) :: timescale_sppt
integer :: sppt_vertstruc
integer :: iseed_sppt
integer , DIMENSION(max_domains) :: spp_mp
real , DIMENSION(max_domains) :: gridpt_stddev_spp_mp
real , DIMENSION(max_domains) :: stddev_cutoff_spp_mp
real , DIMENSION(max_domains) :: lengthscale_spp_mp
real , DIMENSION(max_domains) :: timescale_spp_mp
integer :: vertstruc_spp_mp
integer :: iseed_spp_mp
integer , DIMENSION(max_domains) :: isppt_mp
real , DIMENSION(max_domains) :: gridpt_stddev_isppt_mp
real , DIMENSION(max_domains) :: stddev_cutoff_isppt_mp
real , DIMENSION(max_domains) :: lengthscale_isppt_mp
real , DIMENSION(max_domains) :: timescale_isppt_mp
integer :: vertstruc_isppt_mp
integer :: iseed_isppt_mp
integer , DIMENSION(max_domains) :: isppt_pbl
real , DIMENSION(max_domains) :: gridpt_stddev_isppt_pbl
real , DIMENSION(max_domains) :: stddev_cutoff_isppt_pbl
real , DIMENSION(max_domains) :: lengthscale_isppt_pbl
real , DIMENSION(max_domains) :: timescale_isppt_pbl
integer :: vertstruc_isppt_pbl
integer :: iseed_isppt_pbl
integer , DIMENSION(max_domains) :: isppt_rad
real , DIMENSION(max_domains) :: gridpt_stddev_isppt_rad
real , DIMENSION(max_domains) :: stddev_cutoff_isppt_rad
real , DIMENSION(max_domains) :: lengthscale_isppt_rad
real , DIMENSION(max_domains) :: timescale_isppt_rad
integer :: vertstruc_isppt_rad
integer :: iseed_isppt_rad
integer , DIMENSION(max_domains) :: isppt_conv
real , DIMENSION(max_domains) :: gridpt_stddev_isppt_conv
real , DIMENSION(max_domains) :: stddev_cutoff_isppt_conv
real , DIMENSION(max_domains) :: lengthscale_isppt_conv
real , DIMENSION(max_domains) :: timescale_isppt_conv
integer :: vertstruc_isppt_conv
integer :: iseed_isppt_conv
integer , DIMENSION(max_domains) :: rand_perturb
real , DIMENSION(max_domains) :: gridpt_stddev_rand_pert
real , DIMENSION(max_domains) :: stddev_cutoff_rand_pert
real , DIMENSION(max_domains) :: lengthscale_rand_pert
real , DIMENSION(max_domains) :: timescale_rand_pert
integer :: rand_pert_vertstruc
integer :: iseed_rand_pert
integer :: skebs_on
integer :: sppt_on
integer :: isppt_on
integer :: spp_on
integer :: rand_perturb_on
integer :: num_stoch_levels
integer :: halo_debug
integer :: ntracers
integer , DIMENSION(max_domains) :: vortex_tracker
real , DIMENSION(max_domains) :: interest_rad_storm
real , DIMENSION(max_domains) :: interest_rad_parent
real , DIMENSION(max_domains) :: interest_rad_self
integer , DIMENSION(max_domains) :: interest_kids
integer , DIMENSION(max_domains) :: interest_self
integer , DIMENSION(max_domains) :: interest_storms
integer :: swath_mode
integer :: num_old_fixes
real , DIMENSION(max_domains) :: vt4_radius
real , DIMENSION(max_domains) :: vt4_weightexp
real , DIMENSION(max_domains) :: vt4_pmax
real , DIMENSION(max_domains) :: vt4_noise_pmax
real , DIMENSION(max_domains) :: vt4_noise_pmin
real , DIMENSION(max_domains) :: vt4_noise_dpdr
integer , DIMENSION(max_domains) :: vt4_noise_iter
real , DIMENSION(max_domains) :: nomove_freq
integer :: traj_opt
logical , DIMENSION(max_domains) :: dm_has_traj
integer :: num_traj
integer , DIMENSION(max_domains) :: corral_x
integer , DIMENSION(max_domains) :: corral_y
integer :: tg_reset_stream
integer :: tg_option
integer , DIMENSION(max_domains) :: ntornado
real , DIMENSION(max_domains) :: wbd0
real , DIMENSION(max_domains) :: sbd0
logical , DIMENSION(max_domains) :: analysis
logical , DIMENSION(max_domains) :: write_analysis
integer :: io_form_auxinput2
logical :: high_freq
integer :: high_dom
integer :: swint_opt
integer , DIMENSION(max_domains) :: aer_type
integer , DIMENSION(max_domains) :: aer_aod550_opt
integer , DIMENSION(max_domains) :: aer_angexp_opt
integer , DIMENSION(max_domains) :: aer_ssa_opt
integer , DIMENSION(max_domains) :: aer_asy_opt
real , DIMENSION(max_domains) :: aer_aod550_val
real , DIMENSION(max_domains) :: aer_angexp_val
real , DIMENSION(max_domains) :: aer_ssa_val
real , DIMENSION(max_domains) :: aer_asy_val
integer :: dveg
integer :: opt_crs
integer :: opt_btr
integer :: opt_run
integer :: opt_sfc
integer :: opt_frz
integer :: opt_inf
integer :: opt_rad
integer :: opt_alb
integer :: opt_snf
integer :: opt_tbot
integer :: opt_stc
integer :: opt_gla
integer :: opt_rsf
integer :: opt_soil
integer :: opt_pedo
integer :: opt_crop
integer :: wrf_hydro
integer :: run_days
integer :: run_hours
integer :: run_minutes
integer :: run_seconds
integer , DIMENSION(max_domains) :: start_year
integer , DIMENSION(max_domains) :: start_month
integer , DIMENSION(max_domains) :: start_day
integer , DIMENSION(max_domains) :: start_hour
integer , DIMENSION(max_domains) :: start_minute
integer , DIMENSION(max_domains) :: start_second
integer , DIMENSION(max_domains) :: end_year
integer , DIMENSION(max_domains) :: end_month
integer , DIMENSION(max_domains) :: end_day
integer , DIMENSION(max_domains) :: end_hour
integer , DIMENSION(max_domains) :: end_minute
integer , DIMENSION(max_domains) :: end_second
integer :: interval_seconds
logical , DIMENSION(max_domains) :: input_from_file
integer , DIMENSION(max_domains) :: fine_input_stream
character*256 :: auxinput1_inname
integer :: io_form_auxinput1
logical :: override_restart_timers
character*256 :: auxhist1_inname
character*256 :: auxhist1_outname
integer , DIMENSION(max_domains) :: auxhist1_interval_y
integer , DIMENSION(max_domains) :: auxhist1_interval_d
integer , DIMENSION(max_domains) :: auxhist1_interval_h
integer , DIMENSION(max_domains) :: auxhist1_interval_m
integer , DIMENSION(max_domains) :: auxhist1_interval_s
integer , DIMENSION(max_domains) :: auxhist1_interval
integer , DIMENSION(max_domains) :: auxhist1_begin_y
integer , DIMENSION(max_domains) :: auxhist1_begin_d
integer , DIMENSION(max_domains) :: auxhist1_begin_h
integer , DIMENSION(max_domains) :: auxhist1_begin_m
integer , DIMENSION(max_domains) :: auxhist1_begin_s
integer , DIMENSION(max_domains) :: auxhist1_begin
integer , DIMENSION(max_domains) :: auxhist1_end_y
integer , DIMENSION(max_domains) :: auxhist1_end_d
integer , DIMENSION(max_domains) :: auxhist1_end_h
integer , DIMENSION(max_domains) :: auxhist1_end_m
integer , DIMENSION(max_domains) :: auxhist1_end_s
integer , DIMENSION(max_domains) :: auxhist1_end
integer :: io_form_auxhist1
integer , DIMENSION(max_domains) :: frames_per_auxhist1
character*256 :: auxhist2_inname
character*256 :: auxhist2_outname
integer , DIMENSION(max_domains) :: auxhist2_interval_y
integer , DIMENSION(max_domains) :: auxhist2_interval_d
integer , DIMENSION(max_domains) :: auxhist2_interval_h
integer , DIMENSION(max_domains) :: auxhist2_interval_m
integer , DIMENSION(max_domains) :: auxhist2_interval_s
integer , DIMENSION(max_domains) :: auxhist2_interval
integer , DIMENSION(max_domains) :: auxhist2_begin_y
integer , DIMENSION(max_domains) :: auxhist2_begin_d
integer , DIMENSION(max_domains) :: auxhist2_begin_h
integer , DIMENSION(max_domains) :: auxhist2_begin_m
integer , DIMENSION(max_domains) :: auxhist2_begin_s
integer , DIMENSION(max_domains) :: auxhist2_begin
integer , DIMENSION(max_domains) :: auxhist2_end_y
integer , DIMENSION(max_domains) :: auxhist2_end_d
integer , DIMENSION(max_domains) :: auxhist2_end_h
integer , DIMENSION(max_domains) :: auxhist2_end_m
integer , DIMENSION(max_domains) :: auxhist2_end_s
integer , DIMENSION(max_domains) :: auxhist2_end
integer :: io_form_auxhist2
integer , DIMENSION(max_domains) :: frames_per_auxhist2
character*256 :: auxhist3_inname
character*256 :: auxhist3_outname
integer , DIMENSION(max_domains) :: auxhist3_interval_y
integer , DIMENSION(max_domains) :: auxhist3_interval_d
integer , DIMENSION(max_domains) :: auxhist3_interval_h
integer , DIMENSION(max_domains) :: auxhist3_interval_m
integer , DIMENSION(max_domains) :: auxhist3_interval_s
integer , DIMENSION(max_domains) :: auxhist3_interval
integer , DIMENSION(max_domains) :: auxhist3_begin_y
integer , DIMENSION(max_domains) :: auxhist3_begin_d
integer , DIMENSION(max_domains) :: auxhist3_begin_h
integer , DIMENSION(max_domains) :: auxhist3_begin_m
integer , DIMENSION(max_domains) :: auxhist3_begin_s
integer , DIMENSION(max_domains) :: auxhist3_begin
integer , DIMENSION(max_domains) :: auxhist3_end_y
integer , DIMENSION(max_domains) :: auxhist3_end_d
integer , DIMENSION(max_domains) :: auxhist3_end_h
integer , DIMENSION(max_domains) :: auxhist3_end_m
integer , DIMENSION(max_domains) :: auxhist3_end_s
integer , DIMENSION(max_domains) :: auxhist3_end
integer :: io_form_auxhist3
integer , DIMENSION(max_domains) :: frames_per_auxhist3
character*256 :: auxhist4_inname
character*256 :: auxhist4_outname
integer , DIMENSION(max_domains) :: auxhist4_interval_y
integer , DIMENSION(max_domains) :: auxhist4_interval_d
integer , DIMENSION(max_domains) :: auxhist4_interval_h
integer , DIMENSION(max_domains) :: auxhist4_interval_m
integer , DIMENSION(max_domains) :: auxhist4_interval_s
integer , DIMENSION(max_domains) :: auxhist4_interval
integer , DIMENSION(max_domains) :: auxhist4_begin_y
integer , DIMENSION(max_domains) :: auxhist4_begin_d
integer , DIMENSION(max_domains) :: auxhist4_begin_h
integer , DIMENSION(max_domains) :: auxhist4_begin_m
integer , DIMENSION(max_domains) :: auxhist4_begin_s
integer , DIMENSION(max_domains) :: auxhist4_begin
integer , DIMENSION(max_domains) :: auxhist4_end_y
integer , DIMENSION(max_domains) :: auxhist4_end_d
integer , DIMENSION(max_domains) :: auxhist4_end_h
integer , DIMENSION(max_domains) :: auxhist4_end_m
integer , DIMENSION(max_domains) :: auxhist4_end_s
integer , DIMENSION(max_domains) :: auxhist4_end
integer :: io_form_auxhist4
integer , DIMENSION(max_domains) :: frames_per_auxhist4
character*256 :: auxhist5_inname
character*256 :: auxhist5_outname
integer , DIMENSION(max_domains) :: auxhist5_interval_y
integer , DIMENSION(max_domains) :: auxhist5_interval_d
integer , DIMENSION(max_domains) :: auxhist5_interval_h
integer , DIMENSION(max_domains) :: auxhist5_interval_m
integer , DIMENSION(max_domains) :: auxhist5_interval_s
integer , DIMENSION(max_domains) :: auxhist5_interval
integer , DIMENSION(max_domains) :: auxhist5_begin_y
integer , DIMENSION(max_domains) :: auxhist5_begin_d
integer , DIMENSION(max_domains) :: auxhist5_begin_h
integer , DIMENSION(max_domains) :: auxhist5_begin_m
integer , DIMENSION(max_domains) :: auxhist5_begin_s
integer , DIMENSION(max_domains) :: auxhist5_begin
integer , DIMENSION(max_domains) :: auxhist5_end_y
integer , DIMENSION(max_domains) :: auxhist5_end_d
integer , DIMENSION(max_domains) :: auxhist5_end_h
integer , DIMENSION(max_domains) :: auxhist5_end_m
integer , DIMENSION(max_domains) :: auxhist5_end_s
integer , DIMENSION(max_domains) :: auxhist5_end
integer :: io_form_auxhist5
integer , DIMENSION(max_domains) :: frames_per_auxhist5
character*256 :: auxhist6_inname
character*256 :: auxhist6_outname
integer , DIMENSION(max_domains) :: auxhist6_interval_y
integer , DIMENSION(max_domains) :: auxhist6_interval_d
integer , DIMENSION(max_domains) :: auxhist6_interval_h
integer , DIMENSION(max_domains) :: auxhist6_interval_m
integer , DIMENSION(max_domains) :: auxhist6_interval_s
integer , DIMENSION(max_domains) :: auxhist6_interval
integer , DIMENSION(max_domains) :: auxhist6_begin_y
integer , DIMENSION(max_domains) :: auxhist6_begin_d
integer , DIMENSION(max_domains) :: auxhist6_begin_h
integer , DIMENSION(max_domains) :: auxhist6_begin_m
integer , DIMENSION(max_domains) :: auxhist6_begin_s
integer , DIMENSION(max_domains) :: auxhist6_begin
integer , DIMENSION(max_domains) :: auxhist6_end_y
integer , DIMENSION(max_domains) :: auxhist6_end_d
integer , DIMENSION(max_domains) :: auxhist6_end_h
integer , DIMENSION(max_domains) :: auxhist6_end_m
integer , DIMENSION(max_domains) :: auxhist6_end_s
integer , DIMENSION(max_domains) :: auxhist6_end
integer :: io_form_auxhist6
integer , DIMENSION(max_domains) :: frames_per_auxhist6
character*256 :: auxhist7_inname
character*256 :: auxhist7_outname
integer , DIMENSION(max_domains) :: auxhist7_interval_y
integer , DIMENSION(max_domains) :: auxhist7_interval_d
integer , DIMENSION(max_domains) :: auxhist7_interval_h
integer , DIMENSION(max_domains) :: auxhist7_interval_m
integer , DIMENSION(max_domains) :: auxhist7_interval_s
integer , DIMENSION(max_domains) :: auxhist7_interval
integer , DIMENSION(max_domains) :: auxhist7_begin_y
integer , DIMENSION(max_domains) :: auxhist7_begin_d
integer , DIMENSION(max_domains) :: auxhist7_begin_h
integer , DIMENSION(max_domains) :: auxhist7_begin_m
integer , DIMENSION(max_domains) :: auxhist7_begin_s
integer , DIMENSION(max_domains) :: auxhist7_begin
integer , DIMENSION(max_domains) :: auxhist7_end_y
integer , DIMENSION(max_domains) :: auxhist7_end_d
integer , DIMENSION(max_domains) :: auxhist7_end_h
integer , DIMENSION(max_domains) :: auxhist7_end_m
integer , DIMENSION(max_domains) :: auxhist7_end_s
integer , DIMENSION(max_domains) :: auxhist7_end
integer :: io_form_auxhist7
integer , DIMENSION(max_domains) :: frames_per_auxhist7
character*256 :: auxhist8_inname
character*256 :: auxhist8_outname
integer , DIMENSION(max_domains) :: auxhist8_interval_y
integer , DIMENSION(max_domains) :: auxhist8_interval_d
integer , DIMENSION(max_domains) :: auxhist8_interval_h
integer , DIMENSION(max_domains) :: auxhist8_interval_m
integer , DIMENSION(max_domains) :: auxhist8_interval_s
integer , DIMENSION(max_domains) :: auxhist8_interval
integer , DIMENSION(max_domains) :: auxhist8_begin_y
integer , DIMENSION(max_domains) :: auxhist8_begin_d
integer , DIMENSION(max_domains) :: auxhist8_begin_h
integer , DIMENSION(max_domains) :: auxhist8_begin_m
integer , DIMENSION(max_domains) :: auxhist8_begin_s
integer , DIMENSION(max_domains) :: auxhist8_begin
integer , DIMENSION(max_domains) :: auxhist8_end_y
integer , DIMENSION(max_domains) :: auxhist8_end_d
integer , DIMENSION(max_domains) :: auxhist8_end_h
integer , DIMENSION(max_domains) :: auxhist8_end_m
integer , DIMENSION(max_domains) :: auxhist8_end_s
integer , DIMENSION(max_domains) :: auxhist8_end
integer :: io_form_auxhist8
integer , DIMENSION(max_domains) :: frames_per_auxhist8
character*256 :: auxhist9_inname
character*256 :: auxhist9_outname
integer , DIMENSION(max_domains) :: auxhist9_interval_y
integer , DIMENSION(max_domains) :: auxhist9_interval_d
integer , DIMENSION(max_domains) :: auxhist9_interval_h
integer , DIMENSION(max_domains) :: auxhist9_interval_m
integer , DIMENSION(max_domains) :: auxhist9_interval_s
integer , DIMENSION(max_domains) :: auxhist9_interval
integer , DIMENSION(max_domains) :: auxhist9_begin_y
integer , DIMENSION(max_domains) :: auxhist9_begin_d
integer , DIMENSION(max_domains) :: auxhist9_begin_h
integer , DIMENSION(max_domains) :: auxhist9_begin_m
integer , DIMENSION(max_domains) :: auxhist9_begin_s
integer , DIMENSION(max_domains) :: auxhist9_begin
integer , DIMENSION(max_domains) :: auxhist9_end_y
integer , DIMENSION(max_domains) :: auxhist9_end_d
integer , DIMENSION(max_domains) :: auxhist9_end_h
integer , DIMENSION(max_domains) :: auxhist9_end_m
integer , DIMENSION(max_domains) :: auxhist9_end_s
integer , DIMENSION(max_domains) :: auxhist9_end
integer :: io_form_auxhist9
integer , DIMENSION(max_domains) :: frames_per_auxhist9
character*256 :: auxhist10_inname
character*256 :: auxhist10_outname
integer , DIMENSION(max_domains) :: auxhist10_interval_y
integer , DIMENSION(max_domains) :: auxhist10_interval_d
integer , DIMENSION(max_domains) :: auxhist10_interval_h
integer , DIMENSION(max_domains) :: auxhist10_interval_m
integer , DIMENSION(max_domains) :: auxhist10_interval_s
integer , DIMENSION(max_domains) :: auxhist10_interval
integer , DIMENSION(max_domains) :: auxhist10_begin_y
integer , DIMENSION(max_domains) :: auxhist10_begin_d
integer , DIMENSION(max_domains) :: auxhist10_begin_h
integer , DIMENSION(max_domains) :: auxhist10_begin_m
integer , DIMENSION(max_domains) :: auxhist10_begin_s
integer , DIMENSION(max_domains) :: auxhist10_begin
integer , DIMENSION(max_domains) :: auxhist10_end_y
integer , DIMENSION(max_domains) :: auxhist10_end_d
integer , DIMENSION(max_domains) :: auxhist10_end_h
integer , DIMENSION(max_domains) :: auxhist10_end_m
integer , DIMENSION(max_domains) :: auxhist10_end_s
integer , DIMENSION(max_domains) :: auxhist10_end
integer :: io_form_auxhist10
integer , DIMENSION(max_domains) :: frames_per_auxhist10
character*256 :: auxhist11_inname
character*256 :: auxhist11_outname
integer , DIMENSION(max_domains) :: auxhist11_interval_y
integer , DIMENSION(max_domains) :: auxhist11_interval_d
integer , DIMENSION(max_domains) :: auxhist11_interval_h
integer , DIMENSION(max_domains) :: auxhist11_interval_m
integer , DIMENSION(max_domains) :: auxhist11_interval_s
integer , DIMENSION(max_domains) :: auxhist11_interval
integer , DIMENSION(max_domains) :: auxhist11_begin_y
integer , DIMENSION(max_domains) :: auxhist11_begin_d
integer , DIMENSION(max_domains) :: auxhist11_begin_h
integer , DIMENSION(max_domains) :: auxhist11_begin_m
integer , DIMENSION(max_domains) :: auxhist11_begin_s
integer , DIMENSION(max_domains) :: auxhist11_begin
integer , DIMENSION(max_domains) :: auxhist11_end_y
integer , DIMENSION(max_domains) :: auxhist11_end_d
integer , DIMENSION(max_domains) :: auxhist11_end_h
integer , DIMENSION(max_domains) :: auxhist11_end_m
integer , DIMENSION(max_domains) :: auxhist11_end_s
integer , DIMENSION(max_domains) :: auxhist11_end
integer :: io_form_auxhist11
integer , DIMENSION(max_domains) :: frames_per_auxhist11
character*256 :: auxhist12_inname
character*256 :: auxhist12_outname
integer , DIMENSION(max_domains) :: auxhist12_interval_y
integer , DIMENSION(max_domains) :: auxhist12_interval_d
integer , DIMENSION(max_domains) :: auxhist12_interval_h
integer , DIMENSION(max_domains) :: auxhist12_interval_m
integer , DIMENSION(max_domains) :: auxhist12_interval_s
integer , DIMENSION(max_domains) :: auxhist12_interval
integer , DIMENSION(max_domains) :: auxhist12_begin_y
integer , DIMENSION(max_domains) :: auxhist12_begin_d
integer , DIMENSION(max_domains) :: auxhist12_begin_h
integer , DIMENSION(max_domains) :: auxhist12_begin_m
integer , DIMENSION(max_domains) :: auxhist12_begin_s
integer , DIMENSION(max_domains) :: auxhist12_begin
integer , DIMENSION(max_domains) :: auxhist12_end_y
integer , DIMENSION(max_domains) :: auxhist12_end_d
integer , DIMENSION(max_domains) :: auxhist12_end_h
integer , DIMENSION(max_domains) :: auxhist12_end_m
integer , DIMENSION(max_domains) :: auxhist12_end_s
integer , DIMENSION(max_domains) :: auxhist12_end
integer :: io_form_auxhist12
integer , DIMENSION(max_domains) :: frames_per_auxhist12
character*256 :: auxhist13_inname
character*256 :: auxhist13_outname
integer , DIMENSION(max_domains) :: auxhist13_interval_y
integer , DIMENSION(max_domains) :: auxhist13_interval_d
integer , DIMENSION(max_domains) :: auxhist13_interval_h
integer , DIMENSION(max_domains) :: auxhist13_interval_m
integer , DIMENSION(max_domains) :: auxhist13_interval_s
integer , DIMENSION(max_domains) :: auxhist13_interval
integer , DIMENSION(max_domains) :: auxhist13_begin_y
integer , DIMENSION(max_domains) :: auxhist13_begin_d
integer , DIMENSION(max_domains) :: auxhist13_begin_h
integer , DIMENSION(max_domains) :: auxhist13_begin_m
integer , DIMENSION(max_domains) :: auxhist13_begin_s
integer , DIMENSION(max_domains) :: auxhist13_begin
integer , DIMENSION(max_domains) :: auxhist13_end_y
integer , DIMENSION(max_domains) :: auxhist13_end_d
integer , DIMENSION(max_domains) :: auxhist13_end_h
integer , DIMENSION(max_domains) :: auxhist13_end_m
integer , DIMENSION(max_domains) :: auxhist13_end_s
integer , DIMENSION(max_domains) :: auxhist13_end
integer :: io_form_auxhist13
integer , DIMENSION(max_domains) :: frames_per_auxhist13
character*256 :: auxhist14_inname
character*256 :: auxhist14_outname
integer , DIMENSION(max_domains) :: auxhist14_interval_y
integer , DIMENSION(max_domains) :: auxhist14_interval_d
integer , DIMENSION(max_domains) :: auxhist14_interval_h
integer , DIMENSION(max_domains) :: auxhist14_interval_m
integer , DIMENSION(max_domains) :: auxhist14_interval_s
integer , DIMENSION(max_domains) :: auxhist14_interval
integer , DIMENSION(max_domains) :: auxhist14_begin_y
integer , DIMENSION(max_domains) :: auxhist14_begin_d
integer , DIMENSION(max_domains) :: auxhist14_begin_h
integer , DIMENSION(max_domains) :: auxhist14_begin_m
integer , DIMENSION(max_domains) :: auxhist14_begin_s
integer , DIMENSION(max_domains) :: auxhist14_begin
integer , DIMENSION(max_domains) :: auxhist14_end_y
integer , DIMENSION(max_domains) :: auxhist14_end_d
integer , DIMENSION(max_domains) :: auxhist14_end_h
integer , DIMENSION(max_domains) :: auxhist14_end_m
integer , DIMENSION(max_domains) :: auxhist14_end_s
integer , DIMENSION(max_domains) :: auxhist14_end
integer :: io_form_auxhist14
integer , DIMENSION(max_domains) :: frames_per_auxhist14
character*256 :: auxhist15_inname
character*256 :: auxhist15_outname
integer , DIMENSION(max_domains) :: auxhist15_interval_y
integer , DIMENSION(max_domains) :: auxhist15_interval_d
integer , DIMENSION(max_domains) :: auxhist15_interval_h
integer , DIMENSION(max_domains) :: auxhist15_interval_m
integer , DIMENSION(max_domains) :: auxhist15_interval_s
integer , DIMENSION(max_domains) :: auxhist15_interval
integer , DIMENSION(max_domains) :: auxhist15_begin_y
integer , DIMENSION(max_domains) :: auxhist15_begin_d
integer , DIMENSION(max_domains) :: auxhist15_begin_h
integer , DIMENSION(max_domains) :: auxhist15_begin_m
integer , DIMENSION(max_domains) :: auxhist15_begin_s
integer , DIMENSION(max_domains) :: auxhist15_begin
integer , DIMENSION(max_domains) :: auxhist15_end_y
integer , DIMENSION(max_domains) :: auxhist15_end_d
integer , DIMENSION(max_domains) :: auxhist15_end_h
integer , DIMENSION(max_domains) :: auxhist15_end_m
integer , DIMENSION(max_domains) :: auxhist15_end_s
integer , DIMENSION(max_domains) :: auxhist15_end
integer :: io_form_auxhist15
integer , DIMENSION(max_domains) :: frames_per_auxhist15
character*256 :: auxhist16_inname
character*256 :: auxhist16_outname
integer , DIMENSION(max_domains) :: auxhist16_interval_y
integer , DIMENSION(max_domains) :: auxhist16_interval_d
integer , DIMENSION(max_domains) :: auxhist16_interval_h
integer , DIMENSION(max_domains) :: auxhist16_interval_m
integer , DIMENSION(max_domains) :: auxhist16_interval_s
integer , DIMENSION(max_domains) :: auxhist16_interval
integer , DIMENSION(max_domains) :: auxhist16_begin_y
integer , DIMENSION(max_domains) :: auxhist16_begin_d
integer , DIMENSION(max_domains) :: auxhist16_begin_h
integer , DIMENSION(max_domains) :: auxhist16_begin_m
integer , DIMENSION(max_domains) :: auxhist16_begin_s
integer , DIMENSION(max_domains) :: auxhist16_begin
integer , DIMENSION(max_domains) :: auxhist16_end_y
integer , DIMENSION(max_domains) :: auxhist16_end_d
integer , DIMENSION(max_domains) :: auxhist16_end_h
integer , DIMENSION(max_domains) :: auxhist16_end_m
integer , DIMENSION(max_domains) :: auxhist16_end_s
integer , DIMENSION(max_domains) :: auxhist16_end
integer :: io_form_auxhist16
integer , DIMENSION(max_domains) :: frames_per_auxhist16
character*256 :: auxhist17_inname
character*256 :: auxhist17_outname
integer , DIMENSION(max_domains) :: auxhist17_interval_y
integer , DIMENSION(max_domains) :: auxhist17_interval_d
integer , DIMENSION(max_domains) :: auxhist17_interval_h
integer , DIMENSION(max_domains) :: auxhist17_interval_m
integer , DIMENSION(max_domains) :: auxhist17_interval_s
integer , DIMENSION(max_domains) :: auxhist17_interval
integer , DIMENSION(max_domains) :: auxhist17_begin_y
integer , DIMENSION(max_domains) :: auxhist17_begin_d
integer , DIMENSION(max_domains) :: auxhist17_begin_h
integer , DIMENSION(max_domains) :: auxhist17_begin_m
integer , DIMENSION(max_domains) :: auxhist17_begin_s
integer , DIMENSION(max_domains) :: auxhist17_begin
integer , DIMENSION(max_domains) :: auxhist17_end_y
integer , DIMENSION(max_domains) :: auxhist17_end_d
integer , DIMENSION(max_domains) :: auxhist17_end_h
integer , DIMENSION(max_domains) :: auxhist17_end_m
integer , DIMENSION(max_domains) :: auxhist17_end_s
integer , DIMENSION(max_domains) :: auxhist17_end
integer :: io_form_auxhist17
integer , DIMENSION(max_domains) :: frames_per_auxhist17
character*256 :: auxhist18_inname
character*256 :: auxhist18_outname
integer , DIMENSION(max_domains) :: auxhist18_interval_y
integer , DIMENSION(max_domains) :: auxhist18_interval_d
integer , DIMENSION(max_domains) :: auxhist18_interval_h
integer , DIMENSION(max_domains) :: auxhist18_interval_m
integer , DIMENSION(max_domains) :: auxhist18_interval_s
integer , DIMENSION(max_domains) :: auxhist18_interval
integer , DIMENSION(max_domains) :: auxhist18_begin_y
integer , DIMENSION(max_domains) :: auxhist18_begin_d
integer , DIMENSION(max_domains) :: auxhist18_begin_h
integer , DIMENSION(max_domains) :: auxhist18_begin_m
integer , DIMENSION(max_domains) :: auxhist18_begin_s
integer , DIMENSION(max_domains) :: auxhist18_begin
integer , DIMENSION(max_domains) :: auxhist18_end_y
integer , DIMENSION(max_domains) :: auxhist18_end_d
integer , DIMENSION(max_domains) :: auxhist18_end_h
integer , DIMENSION(max_domains) :: auxhist18_end_m
integer , DIMENSION(max_domains) :: auxhist18_end_s
integer , DIMENSION(max_domains) :: auxhist18_end
integer :: io_form_auxhist18
integer , DIMENSION(max_domains) :: frames_per_auxhist18
character*256 :: auxhist19_inname
character*256 :: auxhist19_outname
integer , DIMENSION(max_domains) :: auxhist19_interval_y
integer , DIMENSION(max_domains) :: auxhist19_interval_d
integer , DIMENSION(max_domains) :: auxhist19_interval_h
integer , DIMENSION(max_domains) :: auxhist19_interval_m
integer , DIMENSION(max_domains) :: auxhist19_interval_s
integer , DIMENSION(max_domains) :: auxhist19_interval
integer , DIMENSION(max_domains) :: auxhist19_begin_y
integer , DIMENSION(max_domains) :: auxhist19_begin_d
integer , DIMENSION(max_domains) :: auxhist19_begin_h
integer , DIMENSION(max_domains) :: auxhist19_begin_m
integer , DIMENSION(max_domains) :: auxhist19_begin_s
integer , DIMENSION(max_domains) :: auxhist19_begin
integer , DIMENSION(max_domains) :: auxhist19_end_y
integer , DIMENSION(max_domains) :: auxhist19_end_d
integer , DIMENSION(max_domains) :: auxhist19_end_h
integer , DIMENSION(max_domains) :: auxhist19_end_m
integer , DIMENSION(max_domains) :: auxhist19_end_s
integer , DIMENSION(max_domains) :: auxhist19_end
integer :: io_form_auxhist19
integer , DIMENSION(max_domains) :: frames_per_auxhist19
character*256 :: auxhist20_inname
character*256 :: auxhist20_outname
integer , DIMENSION(max_domains) :: auxhist20_interval_y
integer , DIMENSION(max_domains) :: auxhist20_interval_d
integer , DIMENSION(max_domains) :: auxhist20_interval_h
integer , DIMENSION(max_domains) :: auxhist20_interval_m
integer , DIMENSION(max_domains) :: auxhist20_interval_s
integer , DIMENSION(max_domains) :: auxhist20_interval
integer , DIMENSION(max_domains) :: auxhist20_begin_y
integer , DIMENSION(max_domains) :: auxhist20_begin_d
integer , DIMENSION(max_domains) :: auxhist20_begin_h
integer , DIMENSION(max_domains) :: auxhist20_begin_m
integer , DIMENSION(max_domains) :: auxhist20_begin_s
integer , DIMENSION(max_domains) :: auxhist20_begin
integer , DIMENSION(max_domains) :: auxhist20_end_y
integer , DIMENSION(max_domains) :: auxhist20_end_d
integer , DIMENSION(max_domains) :: auxhist20_end_h
integer , DIMENSION(max_domains) :: auxhist20_end_m
integer , DIMENSION(max_domains) :: auxhist20_end_s
integer , DIMENSION(max_domains) :: auxhist20_end
integer :: io_form_auxhist20
integer , DIMENSION(max_domains) :: frames_per_auxhist20
character*256 :: auxhist21_inname
character*256 :: auxhist21_outname
integer , DIMENSION(max_domains) :: auxhist21_interval_y
integer , DIMENSION(max_domains) :: auxhist21_interval_d
integer , DIMENSION(max_domains) :: auxhist21_interval_h
integer , DIMENSION(max_domains) :: auxhist21_interval_m
integer , DIMENSION(max_domains) :: auxhist21_interval_s
integer , DIMENSION(max_domains) :: auxhist21_interval
integer , DIMENSION(max_domains) :: auxhist21_begin_y
integer , DIMENSION(max_domains) :: auxhist21_begin_d
integer , DIMENSION(max_domains) :: auxhist21_begin_h
integer , DIMENSION(max_domains) :: auxhist21_begin_m
integer , DIMENSION(max_domains) :: auxhist21_begin_s
integer , DIMENSION(max_domains) :: auxhist21_begin
integer , DIMENSION(max_domains) :: auxhist21_end_y
integer , DIMENSION(max_domains) :: auxhist21_end_d
integer , DIMENSION(max_domains) :: auxhist21_end_h
integer , DIMENSION(max_domains) :: auxhist21_end_m
integer , DIMENSION(max_domains) :: auxhist21_end_s
integer , DIMENSION(max_domains) :: auxhist21_end
integer :: io_form_auxhist21
integer , DIMENSION(max_domains) :: frames_per_auxhist21
character*256 :: auxhist22_inname
character*256 :: auxhist22_outname
integer , DIMENSION(max_domains) :: auxhist22_interval_y
integer , DIMENSION(max_domains) :: auxhist22_interval_d
integer , DIMENSION(max_domains) :: auxhist22_interval_h
integer , DIMENSION(max_domains) :: auxhist22_interval_m
integer , DIMENSION(max_domains) :: auxhist22_interval_s
integer , DIMENSION(max_domains) :: auxhist22_interval
integer , DIMENSION(max_domains) :: auxhist22_begin_y
integer , DIMENSION(max_domains) :: auxhist22_begin_d
integer , DIMENSION(max_domains) :: auxhist22_begin_h
integer , DIMENSION(max_domains) :: auxhist22_begin_m
integer , DIMENSION(max_domains) :: auxhist22_begin_s
integer , DIMENSION(max_domains) :: auxhist22_begin
integer , DIMENSION(max_domains) :: auxhist22_end_y
integer , DIMENSION(max_domains) :: auxhist22_end_d
integer , DIMENSION(max_domains) :: auxhist22_end_h
integer , DIMENSION(max_domains) :: auxhist22_end_m
integer , DIMENSION(max_domains) :: auxhist22_end_s
integer , DIMENSION(max_domains) :: auxhist22_end
integer :: io_form_auxhist22
integer , DIMENSION(max_domains) :: frames_per_auxhist22
character*256 :: auxhist23_inname
character*256 :: auxhist23_outname
integer , DIMENSION(max_domains) :: auxhist23_interval_y
integer , DIMENSION(max_domains) :: auxhist23_interval_d
integer , DIMENSION(max_domains) :: auxhist23_interval_h
integer , DIMENSION(max_domains) :: auxhist23_interval_m
integer , DIMENSION(max_domains) :: auxhist23_interval_s
integer , DIMENSION(max_domains) :: auxhist23_interval
integer , DIMENSION(max_domains) :: auxhist23_begin_y
integer , DIMENSION(max_domains) :: auxhist23_begin_d
integer , DIMENSION(max_domains) :: auxhist23_begin_h
integer , DIMENSION(max_domains) :: auxhist23_begin_m
integer , DIMENSION(max_domains) :: auxhist23_begin_s
integer , DIMENSION(max_domains) :: auxhist23_begin
integer , DIMENSION(max_domains) :: auxhist23_end_y
integer , DIMENSION(max_domains) :: auxhist23_end_d
integer , DIMENSION(max_domains) :: auxhist23_end_h
integer , DIMENSION(max_domains) :: auxhist23_end_m
integer , DIMENSION(max_domains) :: auxhist23_end_s
integer , DIMENSION(max_domains) :: auxhist23_end
integer :: io_form_auxhist23
integer , DIMENSION(max_domains) :: frames_per_auxhist23
character*256 :: auxhist24_inname
character*256 :: auxhist24_outname
integer , DIMENSION(max_domains) :: auxhist24_interval_y
integer , DIMENSION(max_domains) :: auxhist24_interval_d
integer , DIMENSION(max_domains) :: auxhist24_interval_h
integer , DIMENSION(max_domains) :: auxhist24_interval_m
integer , DIMENSION(max_domains) :: auxhist24_interval_s
integer , DIMENSION(max_domains) :: auxhist24_interval
integer , DIMENSION(max_domains) :: auxhist24_begin_y
integer , DIMENSION(max_domains) :: auxhist24_begin_d
integer , DIMENSION(max_domains) :: auxhist24_begin_h
integer , DIMENSION(max_domains) :: auxhist24_begin_m
integer , DIMENSION(max_domains) :: auxhist24_begin_s
integer , DIMENSION(max_domains) :: auxhist24_begin
integer , DIMENSION(max_domains) :: auxhist24_end_y
integer , DIMENSION(max_domains) :: auxhist24_end_d
integer , DIMENSION(max_domains) :: auxhist24_end_h
integer , DIMENSION(max_domains) :: auxhist24_end_m
integer , DIMENSION(max_domains) :: auxhist24_end_s
integer , DIMENSION(max_domains) :: auxhist24_end
integer :: io_form_auxhist24
integer , DIMENSION(max_domains) :: frames_per_auxhist24
character*256 :: auxinput1_outname
integer , DIMENSION(max_domains) :: auxinput1_interval_y
integer , DIMENSION(max_domains) :: auxinput1_interval_d
integer , DIMENSION(max_domains) :: auxinput1_interval_h
integer , DIMENSION(max_domains) :: auxinput1_interval_m
integer , DIMENSION(max_domains) :: auxinput1_interval_s
integer , DIMENSION(max_domains) :: auxinput1_interval
integer , DIMENSION(max_domains) :: auxinput1_begin_y
integer , DIMENSION(max_domains) :: auxinput1_begin_d
integer , DIMENSION(max_domains) :: auxinput1_begin_h
integer , DIMENSION(max_domains) :: auxinput1_begin_m
integer , DIMENSION(max_domains) :: auxinput1_begin_s
integer , DIMENSION(max_domains) :: auxinput1_begin
integer , DIMENSION(max_domains) :: auxinput1_end_y
integer , DIMENSION(max_domains) :: auxinput1_end_d
integer , DIMENSION(max_domains) :: auxinput1_end_h
integer , DIMENSION(max_domains) :: auxinput1_end_m
integer , DIMENSION(max_domains) :: auxinput1_end_s
integer , DIMENSION(max_domains) :: auxinput1_end
integer , DIMENSION(max_domains) :: frames_per_auxinput1
character*256 :: auxinput2_inname
character*256 :: auxinput2_outname
integer , DIMENSION(max_domains) :: auxinput2_interval_y
integer , DIMENSION(max_domains) :: auxinput2_interval_d
integer , DIMENSION(max_domains) :: auxinput2_interval_h
integer , DIMENSION(max_domains) :: auxinput2_interval_m
integer , DIMENSION(max_domains) :: auxinput2_interval_s
integer , DIMENSION(max_domains) :: auxinput2_interval
integer , DIMENSION(max_domains) :: auxinput2_begin_y
integer , DIMENSION(max_domains) :: auxinput2_begin_d
integer , DIMENSION(max_domains) :: auxinput2_begin_h
integer , DIMENSION(max_domains) :: auxinput2_begin_m
integer , DIMENSION(max_domains) :: auxinput2_begin_s
integer , DIMENSION(max_domains) :: auxinput2_begin
integer , DIMENSION(max_domains) :: auxinput2_end_y
integer , DIMENSION(max_domains) :: auxinput2_end_d
integer , DIMENSION(max_domains) :: auxinput2_end_h
integer , DIMENSION(max_domains) :: auxinput2_end_m
integer , DIMENSION(max_domains) :: auxinput2_end_s
integer , DIMENSION(max_domains) :: auxinput2_end
integer , DIMENSION(max_domains) :: frames_per_auxinput2
character*256 :: auxinput3_inname
character*256 :: auxinput3_outname
integer , DIMENSION(max_domains) :: auxinput3_interval_y
integer , DIMENSION(max_domains) :: auxinput3_interval_d
integer , DIMENSION(max_domains) :: auxinput3_interval_h
integer , DIMENSION(max_domains) :: auxinput3_interval_m
integer , DIMENSION(max_domains) :: auxinput3_interval_s
integer , DIMENSION(max_domains) :: auxinput3_interval
integer , DIMENSION(max_domains) :: auxinput3_begin_y
integer , DIMENSION(max_domains) :: auxinput3_begin_d
integer , DIMENSION(max_domains) :: auxinput3_begin_h
integer , DIMENSION(max_domains) :: auxinput3_begin_m
integer , DIMENSION(max_domains) :: auxinput3_begin_s
integer , DIMENSION(max_domains) :: auxinput3_begin
integer , DIMENSION(max_domains) :: auxinput3_end_y
integer , DIMENSION(max_domains) :: auxinput3_end_d
integer , DIMENSION(max_domains) :: auxinput3_end_h
integer , DIMENSION(max_domains) :: auxinput3_end_m
integer , DIMENSION(max_domains) :: auxinput3_end_s
integer , DIMENSION(max_domains) :: auxinput3_end
integer :: io_form_auxinput3
integer , DIMENSION(max_domains) :: frames_per_auxinput3
character*256 :: auxinput4_inname
character*256 :: auxinput4_outname
integer , DIMENSION(max_domains) :: auxinput4_interval_y
integer , DIMENSION(max_domains) :: auxinput4_interval_d
integer , DIMENSION(max_domains) :: auxinput4_interval_h
integer , DIMENSION(max_domains) :: auxinput4_interval_m
integer , DIMENSION(max_domains) :: auxinput4_interval_s
integer , DIMENSION(max_domains) :: auxinput4_interval
integer , DIMENSION(max_domains) :: auxinput4_begin_y
integer , DIMENSION(max_domains) :: auxinput4_begin_d
integer , DIMENSION(max_domains) :: auxinput4_begin_h
integer , DIMENSION(max_domains) :: auxinput4_begin_m
integer , DIMENSION(max_domains) :: auxinput4_begin_s
integer , DIMENSION(max_domains) :: auxinput4_begin
integer , DIMENSION(max_domains) :: auxinput4_end_y
integer , DIMENSION(max_domains) :: auxinput4_end_d
integer , DIMENSION(max_domains) :: auxinput4_end_h
integer , DIMENSION(max_domains) :: auxinput4_end_m
integer , DIMENSION(max_domains) :: auxinput4_end_s
integer , DIMENSION(max_domains) :: auxinput4_end
integer :: io_form_auxinput4
integer , DIMENSION(max_domains) :: frames_per_auxinput4
character*256 :: auxinput5_inname
character*256 :: auxinput5_outname
integer , DIMENSION(max_domains) :: auxinput5_interval_y
integer , DIMENSION(max_domains) :: auxinput5_interval_d
integer , DIMENSION(max_domains) :: auxinput5_interval_h
integer , DIMENSION(max_domains) :: auxinput5_interval_m
integer , DIMENSION(max_domains) :: auxinput5_interval_s
integer , DIMENSION(max_domains) :: auxinput5_interval
integer , DIMENSION(max_domains) :: auxinput5_begin_y
integer , DIMENSION(max_domains) :: auxinput5_begin_d
integer , DIMENSION(max_domains) :: auxinput5_begin_h
integer , DIMENSION(max_domains) :: auxinput5_begin_m
integer , DIMENSION(max_domains) :: auxinput5_begin_s
integer , DIMENSION(max_domains) :: auxinput5_begin
integer , DIMENSION(max_domains) :: auxinput5_end_y
integer , DIMENSION(max_domains) :: auxinput5_end_d
integer , DIMENSION(max_domains) :: auxinput5_end_h
integer , DIMENSION(max_domains) :: auxinput5_end_m
integer , DIMENSION(max_domains) :: auxinput5_end_s
integer , DIMENSION(max_domains) :: auxinput5_end
integer :: io_form_auxinput5
integer , DIMENSION(max_domains) :: frames_per_auxinput5
character*256 :: auxinput6_inname
character*256 :: auxinput6_outname
integer , DIMENSION(max_domains) :: auxinput6_interval_y
integer , DIMENSION(max_domains) :: auxinput6_interval_d
integer , DIMENSION(max_domains) :: auxinput6_interval_h
integer , DIMENSION(max_domains) :: auxinput6_interval_m
integer , DIMENSION(max_domains) :: auxinput6_interval_s
integer , DIMENSION(max_domains) :: auxinput6_interval
integer , DIMENSION(max_domains) :: auxinput6_begin_y
integer , DIMENSION(max_domains) :: auxinput6_begin_d
integer , DIMENSION(max_domains) :: auxinput6_begin_h
integer , DIMENSION(max_domains) :: auxinput6_begin_m
integer , DIMENSION(max_domains) :: auxinput6_begin_s
integer , DIMENSION(max_domains) :: auxinput6_begin
integer , DIMENSION(max_domains) :: auxinput6_end_y
integer , DIMENSION(max_domains) :: auxinput6_end_d
integer , DIMENSION(max_domains) :: auxinput6_end_h
integer , DIMENSION(max_domains) :: auxinput6_end_m
integer , DIMENSION(max_domains) :: auxinput6_end_s
integer , DIMENSION(max_domains) :: auxinput6_end
integer :: io_form_auxinput6
integer , DIMENSION(max_domains) :: frames_per_auxinput6
character*256 :: auxinput7_inname
character*256 :: auxinput7_outname
integer , DIMENSION(max_domains) :: auxinput7_interval_y
integer , DIMENSION(max_domains) :: auxinput7_interval_d
integer , DIMENSION(max_domains) :: auxinput7_interval_h
integer , DIMENSION(max_domains) :: auxinput7_interval_m
integer , DIMENSION(max_domains) :: auxinput7_interval_s
integer , DIMENSION(max_domains) :: auxinput7_interval
integer , DIMENSION(max_domains) :: auxinput7_begin_y
integer , DIMENSION(max_domains) :: auxinput7_begin_d
integer , DIMENSION(max_domains) :: auxinput7_begin_h
integer , DIMENSION(max_domains) :: auxinput7_begin_m
integer , DIMENSION(max_domains) :: auxinput7_begin_s
integer , DIMENSION(max_domains) :: auxinput7_begin
integer , DIMENSION(max_domains) :: auxinput7_end_y
integer , DIMENSION(max_domains) :: auxinput7_end_d
integer , DIMENSION(max_domains) :: auxinput7_end_h
integer , DIMENSION(max_domains) :: auxinput7_end_m
integer , DIMENSION(max_domains) :: auxinput7_end_s
integer , DIMENSION(max_domains) :: auxinput7_end
integer :: io_form_auxinput7
integer , DIMENSION(max_domains) :: frames_per_auxinput7
character*256 :: auxinput8_inname
character*256 :: auxinput8_outname
integer , DIMENSION(max_domains) :: auxinput8_interval_y
integer , DIMENSION(max_domains) :: auxinput8_interval_d
integer , DIMENSION(max_domains) :: auxinput8_interval_h
integer , DIMENSION(max_domains) :: auxinput8_interval_m
integer , DIMENSION(max_domains) :: auxinput8_interval_s
integer , DIMENSION(max_domains) :: auxinput8_interval
integer , DIMENSION(max_domains) :: auxinput8_begin_y
integer , DIMENSION(max_domains) :: auxinput8_begin_d
integer , DIMENSION(max_domains) :: auxinput8_begin_h
integer , DIMENSION(max_domains) :: auxinput8_begin_m
integer , DIMENSION(max_domains) :: auxinput8_begin_s
integer , DIMENSION(max_domains) :: auxinput8_begin
integer , DIMENSION(max_domains) :: auxinput8_end_y
integer , DIMENSION(max_domains) :: auxinput8_end_d
integer , DIMENSION(max_domains) :: auxinput8_end_h
integer , DIMENSION(max_domains) :: auxinput8_end_m
integer , DIMENSION(max_domains) :: auxinput8_end_s
integer , DIMENSION(max_domains) :: auxinput8_end
integer :: io_form_auxinput8
integer , DIMENSION(max_domains) :: frames_per_auxinput8
character*256 :: auxinput9_inname
character*256 :: auxinput9_outname
integer , DIMENSION(max_domains) :: auxinput9_interval_y
integer , DIMENSION(max_domains) :: auxinput9_interval_d
integer , DIMENSION(max_domains) :: auxinput9_interval_h
integer , DIMENSION(max_domains) :: auxinput9_interval_m
integer , DIMENSION(max_domains) :: auxinput9_interval_s
integer , DIMENSION(max_domains) :: auxinput9_interval
integer , DIMENSION(max_domains) :: auxinput9_begin_y
integer , DIMENSION(max_domains) :: auxinput9_begin_d
integer , DIMENSION(max_domains) :: auxinput9_begin_h
integer , DIMENSION(max_domains) :: auxinput9_begin_m
integer , DIMENSION(max_domains) :: auxinput9_begin_s
integer , DIMENSION(max_domains) :: auxinput9_begin
integer , DIMENSION(max_domains) :: auxinput9_end_y
integer , DIMENSION(max_domains) :: auxinput9_end_d
integer , DIMENSION(max_domains) :: auxinput9_end_h
integer , DIMENSION(max_domains) :: auxinput9_end_m
integer , DIMENSION(max_domains) :: auxinput9_end_s
integer , DIMENSION(max_domains) :: auxinput9_end
integer :: io_form_auxinput9
integer , DIMENSION(max_domains) :: frames_per_auxinput9
character*256 :: auxinput10_inname
character*256 :: auxinput10_outname
integer , DIMENSION(max_domains) :: auxinput10_interval_y
integer , DIMENSION(max_domains) :: auxinput10_interval_d
integer , DIMENSION(max_domains) :: auxinput10_interval_h
integer , DIMENSION(max_domains) :: auxinput10_interval_m
integer , DIMENSION(max_domains) :: auxinput10_interval_s
integer , DIMENSION(max_domains) :: auxinput10_interval
integer , DIMENSION(max_domains) :: auxinput10_begin_y
integer , DIMENSION(max_domains) :: auxinput10_begin_d
integer , DIMENSION(max_domains) :: auxinput10_begin_h
integer , DIMENSION(max_domains) :: auxinput10_begin_m
integer , DIMENSION(max_domains) :: auxinput10_begin_s
integer , DIMENSION(max_domains) :: auxinput10_begin
integer , DIMENSION(max_domains) :: auxinput10_end_y
integer , DIMENSION(max_domains) :: auxinput10_end_d
integer , DIMENSION(max_domains) :: auxinput10_end_h
integer , DIMENSION(max_domains) :: auxinput10_end_m
integer , DIMENSION(max_domains) :: auxinput10_end_s
integer , DIMENSION(max_domains) :: auxinput10_end
integer :: io_form_auxinput10
integer , DIMENSION(max_domains) :: frames_per_auxinput10
character*256 :: auxinput11_inname
character*256 :: auxinput11_outname
integer , DIMENSION(max_domains) :: auxinput11_interval_y
integer , DIMENSION(max_domains) :: auxinput11_interval_d
integer , DIMENSION(max_domains) :: auxinput11_interval_h
integer , DIMENSION(max_domains) :: auxinput11_interval_m
integer , DIMENSION(max_domains) :: auxinput11_interval_s
integer , DIMENSION(max_domains) :: auxinput11_interval
integer , DIMENSION(max_domains) :: auxinput11_begin_y
integer , DIMENSION(max_domains) :: auxinput11_begin_d
integer , DIMENSION(max_domains) :: auxinput11_begin_h
integer , DIMENSION(max_domains) :: auxinput11_begin_m
integer , DIMENSION(max_domains) :: auxinput11_begin_s
integer , DIMENSION(max_domains) :: auxinput11_begin
integer , DIMENSION(max_domains) :: auxinput11_end_y
integer , DIMENSION(max_domains) :: auxinput11_end_d
integer , DIMENSION(max_domains) :: auxinput11_end_h
integer , DIMENSION(max_domains) :: auxinput11_end_m
integer , DIMENSION(max_domains) :: auxinput11_end_s
integer , DIMENSION(max_domains) :: auxinput11_end
integer :: io_form_auxinput11
integer , DIMENSION(max_domains) :: frames_per_auxinput11
character*256 :: auxinput12_inname
character*256 :: auxinput12_outname
integer , DIMENSION(max_domains) :: auxinput12_interval_y
integer , DIMENSION(max_domains) :: auxinput12_interval_d
integer , DIMENSION(max_domains) :: auxinput12_interval_h
integer , DIMENSION(max_domains) :: auxinput12_interval_m
integer , DIMENSION(max_domains) :: auxinput12_interval_s
integer , DIMENSION(max_domains) :: auxinput12_interval
integer , DIMENSION(max_domains) :: auxinput12_begin_y
integer , DIMENSION(max_domains) :: auxinput12_begin_d
integer , DIMENSION(max_domains) :: auxinput12_begin_h
integer , DIMENSION(max_domains) :: auxinput12_begin_m
integer , DIMENSION(max_domains) :: auxinput12_begin_s
integer , DIMENSION(max_domains) :: auxinput12_begin
integer , DIMENSION(max_domains) :: auxinput12_end_y
integer , DIMENSION(max_domains) :: auxinput12_end_d
integer , DIMENSION(max_domains) :: auxinput12_end_h
integer , DIMENSION(max_domains) :: auxinput12_end_m
integer , DIMENSION(max_domains) :: auxinput12_end_s
integer , DIMENSION(max_domains) :: auxinput12_end
integer :: io_form_auxinput12
integer , DIMENSION(max_domains) :: frames_per_auxinput12
character*256 :: auxinput13_inname
character*256 :: auxinput13_outname
integer , DIMENSION(max_domains) :: auxinput13_interval_y
integer , DIMENSION(max_domains) :: auxinput13_interval_d
integer , DIMENSION(max_domains) :: auxinput13_interval_h
integer , DIMENSION(max_domains) :: auxinput13_interval_m
integer , DIMENSION(max_domains) :: auxinput13_interval_s
integer , DIMENSION(max_domains) :: auxinput13_interval
integer , DIMENSION(max_domains) :: auxinput13_begin_y
integer , DIMENSION(max_domains) :: auxinput13_begin_d
integer , DIMENSION(max_domains) :: auxinput13_begin_h
integer , DIMENSION(max_domains) :: auxinput13_begin_m
integer , DIMENSION(max_domains) :: auxinput13_begin_s
integer , DIMENSION(max_domains) :: auxinput13_begin
integer , DIMENSION(max_domains) :: auxinput13_end_y
integer , DIMENSION(max_domains) :: auxinput13_end_d
integer , DIMENSION(max_domains) :: auxinput13_end_h
integer , DIMENSION(max_domains) :: auxinput13_end_m
integer , DIMENSION(max_domains) :: auxinput13_end_s
integer , DIMENSION(max_domains) :: auxinput13_end
integer :: io_form_auxinput13
integer , DIMENSION(max_domains) :: frames_per_auxinput13
character*256 :: auxinput14_inname
character*256 :: auxinput14_outname
integer , DIMENSION(max_domains) :: auxinput14_interval_y
integer , DIMENSION(max_domains) :: auxinput14_interval_d
integer , DIMENSION(max_domains) :: auxinput14_interval_h
integer , DIMENSION(max_domains) :: auxinput14_interval_m
integer , DIMENSION(max_domains) :: auxinput14_interval_s
integer , DIMENSION(max_domains) :: auxinput14_interval
integer , DIMENSION(max_domains) :: auxinput14_begin_y
integer , DIMENSION(max_domains) :: auxinput14_begin_d
integer , DIMENSION(max_domains) :: auxinput14_begin_h
integer , DIMENSION(max_domains) :: auxinput14_begin_m
integer , DIMENSION(max_domains) :: auxinput14_begin_s
integer , DIMENSION(max_domains) :: auxinput14_begin
integer , DIMENSION(max_domains) :: auxinput14_end_y
integer , DIMENSION(max_domains) :: auxinput14_end_d
integer , DIMENSION(max_domains) :: auxinput14_end_h
integer , DIMENSION(max_domains) :: auxinput14_end_m
integer , DIMENSION(max_domains) :: auxinput14_end_s
integer , DIMENSION(max_domains) :: auxinput14_end
integer :: io_form_auxinput14
integer , DIMENSION(max_domains) :: frames_per_auxinput14
character*256 :: auxinput15_inname
character*256 :: auxinput15_outname
integer , DIMENSION(max_domains) :: auxinput15_interval_y
integer , DIMENSION(max_domains) :: auxinput15_interval_d
integer , DIMENSION(max_domains) :: auxinput15_interval_h
integer , DIMENSION(max_domains) :: auxinput15_interval_m
integer , DIMENSION(max_domains) :: auxinput15_interval_s
integer , DIMENSION(max_domains) :: auxinput15_interval
integer , DIMENSION(max_domains) :: auxinput15_begin_y
integer , DIMENSION(max_domains) :: auxinput15_begin_d
integer , DIMENSION(max_domains) :: auxinput15_begin_h
integer , DIMENSION(max_domains) :: auxinput15_begin_m
integer , DIMENSION(max_domains) :: auxinput15_begin_s
integer , DIMENSION(max_domains) :: auxinput15_begin
integer , DIMENSION(max_domains) :: auxinput15_end_y
integer , DIMENSION(max_domains) :: auxinput15_end_d
integer , DIMENSION(max_domains) :: auxinput15_end_h
integer , DIMENSION(max_domains) :: auxinput15_end_m
integer , DIMENSION(max_domains) :: auxinput15_end_s
integer , DIMENSION(max_domains) :: auxinput15_end
integer :: io_form_auxinput15
integer , DIMENSION(max_domains) :: frames_per_auxinput15
character*256 :: auxinput16_inname
character*256 :: auxinput16_outname
integer , DIMENSION(max_domains) :: auxinput16_interval_y
integer , DIMENSION(max_domains) :: auxinput16_interval_d
integer , DIMENSION(max_domains) :: auxinput16_interval_h
integer , DIMENSION(max_domains) :: auxinput16_interval_m
integer , DIMENSION(max_domains) :: auxinput16_interval_s
integer , DIMENSION(max_domains) :: auxinput16_interval
integer , DIMENSION(max_domains) :: auxinput16_begin_y
integer , DIMENSION(max_domains) :: auxinput16_begin_d
integer , DIMENSION(max_domains) :: auxinput16_begin_h
integer , DIMENSION(max_domains) :: auxinput16_begin_m
integer , DIMENSION(max_domains) :: auxinput16_begin_s
integer , DIMENSION(max_domains) :: auxinput16_begin
integer , DIMENSION(max_domains) :: auxinput16_end_y
integer , DIMENSION(max_domains) :: auxinput16_end_d
integer , DIMENSION(max_domains) :: auxinput16_end_h
integer , DIMENSION(max_domains) :: auxinput16_end_m
integer , DIMENSION(max_domains) :: auxinput16_end_s
integer , DIMENSION(max_domains) :: auxinput16_end
integer :: io_form_auxinput16
integer , DIMENSION(max_domains) :: frames_per_auxinput16
character*256 :: auxinput17_inname
character*256 :: auxinput17_outname
integer , DIMENSION(max_domains) :: auxinput17_interval_y
integer , DIMENSION(max_domains) :: auxinput17_interval_d
integer , DIMENSION(max_domains) :: auxinput17_interval_h
integer , DIMENSION(max_domains) :: auxinput17_interval_m
integer , DIMENSION(max_domains) :: auxinput17_interval_s
integer , DIMENSION(max_domains) :: auxinput17_interval
integer , DIMENSION(max_domains) :: auxinput17_begin_y
integer , DIMENSION(max_domains) :: auxinput17_begin_d
integer , DIMENSION(max_domains) :: auxinput17_begin_h
integer , DIMENSION(max_domains) :: auxinput17_begin_m
integer , DIMENSION(max_domains) :: auxinput17_begin_s
integer , DIMENSION(max_domains) :: auxinput17_begin
integer , DIMENSION(max_domains) :: auxinput17_end_y
integer , DIMENSION(max_domains) :: auxinput17_end_d
integer , DIMENSION(max_domains) :: auxinput17_end_h
integer , DIMENSION(max_domains) :: auxinput17_end_m
integer , DIMENSION(max_domains) :: auxinput17_end_s
integer , DIMENSION(max_domains) :: auxinput17_end
integer :: io_form_auxinput17
integer , DIMENSION(max_domains) :: frames_per_auxinput17
character*256 :: auxinput18_inname
character*256 :: auxinput18_outname
integer , DIMENSION(max_domains) :: auxinput18_interval_y
integer , DIMENSION(max_domains) :: auxinput18_interval_d
integer , DIMENSION(max_domains) :: auxinput18_interval_h
integer , DIMENSION(max_domains) :: auxinput18_interval_m
integer , DIMENSION(max_domains) :: auxinput18_interval_s
integer , DIMENSION(max_domains) :: auxinput18_interval
integer , DIMENSION(max_domains) :: auxinput18_begin_y
integer , DIMENSION(max_domains) :: auxinput18_begin_d
integer , DIMENSION(max_domains) :: auxinput18_begin_h
integer , DIMENSION(max_domains) :: auxinput18_begin_m
integer , DIMENSION(max_domains) :: auxinput18_begin_s
integer , DIMENSION(max_domains) :: auxinput18_begin
integer , DIMENSION(max_domains) :: auxinput18_end_y
integer , DIMENSION(max_domains) :: auxinput18_end_d
integer , DIMENSION(max_domains) :: auxinput18_end_h
integer , DIMENSION(max_domains) :: auxinput18_end_m
integer , DIMENSION(max_domains) :: auxinput18_end_s
integer , DIMENSION(max_domains) :: auxinput18_end
integer :: io_form_auxinput18
integer , DIMENSION(max_domains) :: frames_per_auxinput18
character*256 :: auxinput19_inname
character*256 :: auxinput19_outname
integer , DIMENSION(max_domains) :: auxinput19_interval_y
integer , DIMENSION(max_domains) :: auxinput19_interval_d
integer , DIMENSION(max_domains) :: auxinput19_interval_h
integer , DIMENSION(max_domains) :: auxinput19_interval_m
integer , DIMENSION(max_domains) :: auxinput19_interval_s
integer , DIMENSION(max_domains) :: auxinput19_interval
integer , DIMENSION(max_domains) :: auxinput19_begin_y
integer , DIMENSION(max_domains) :: auxinput19_begin_d
integer , DIMENSION(max_domains) :: auxinput19_begin_h
integer , DIMENSION(max_domains) :: auxinput19_begin_m
integer , DIMENSION(max_domains) :: auxinput19_begin_s
integer , DIMENSION(max_domains) :: auxinput19_begin
integer , DIMENSION(max_domains) :: auxinput19_end_y
integer , DIMENSION(max_domains) :: auxinput19_end_d
integer , DIMENSION(max_domains) :: auxinput19_end_h
integer , DIMENSION(max_domains) :: auxinput19_end_m
integer , DIMENSION(max_domains) :: auxinput19_end_s
integer , DIMENSION(max_domains) :: auxinput19_end
integer :: io_form_auxinput19
integer , DIMENSION(max_domains) :: frames_per_auxinput19
character*256 :: auxinput20_inname
character*256 :: auxinput20_outname
integer , DIMENSION(max_domains) :: auxinput20_interval_y
integer , DIMENSION(max_domains) :: auxinput20_interval_d
integer , DIMENSION(max_domains) :: auxinput20_interval_h
integer , DIMENSION(max_domains) :: auxinput20_interval_m
integer , DIMENSION(max_domains) :: auxinput20_interval_s
integer , DIMENSION(max_domains) :: auxinput20_interval
integer , DIMENSION(max_domains) :: auxinput20_begin_y
integer , DIMENSION(max_domains) :: auxinput20_begin_d
integer , DIMENSION(max_domains) :: auxinput20_begin_h
integer , DIMENSION(max_domains) :: auxinput20_begin_m
integer , DIMENSION(max_domains) :: auxinput20_begin_s
integer , DIMENSION(max_domains) :: auxinput20_begin
integer , DIMENSION(max_domains) :: auxinput20_end_y
integer , DIMENSION(max_domains) :: auxinput20_end_d
integer , DIMENSION(max_domains) :: auxinput20_end_h
integer , DIMENSION(max_domains) :: auxinput20_end_m
integer , DIMENSION(max_domains) :: auxinput20_end_s
integer , DIMENSION(max_domains) :: auxinput20_end
integer :: io_form_auxinput20
integer , DIMENSION(max_domains) :: frames_per_auxinput20
character*256 :: auxinput21_inname
character*256 :: auxinput21_outname
integer , DIMENSION(max_domains) :: auxinput21_interval_y
integer , DIMENSION(max_domains) :: auxinput21_interval_d
integer , DIMENSION(max_domains) :: auxinput21_interval_h
integer , DIMENSION(max_domains) :: auxinput21_interval_m
integer , DIMENSION(max_domains) :: auxinput21_interval_s
integer , DIMENSION(max_domains) :: auxinput21_interval
integer , DIMENSION(max_domains) :: auxinput21_begin_y
integer , DIMENSION(max_domains) :: auxinput21_begin_d
integer , DIMENSION(max_domains) :: auxinput21_begin_h
integer , DIMENSION(max_domains) :: auxinput21_begin_m
integer , DIMENSION(max_domains) :: auxinput21_begin_s
integer , DIMENSION(max_domains) :: auxinput21_begin
integer , DIMENSION(max_domains) :: auxinput21_end_y
integer , DIMENSION(max_domains) :: auxinput21_end_d
integer , DIMENSION(max_domains) :: auxinput21_end_h
integer , DIMENSION(max_domains) :: auxinput21_end_m
integer , DIMENSION(max_domains) :: auxinput21_end_s
integer , DIMENSION(max_domains) :: auxinput21_end
integer :: io_form_auxinput21
integer , DIMENSION(max_domains) :: frames_per_auxinput21
character*256 :: auxinput22_inname
character*256 :: auxinput22_outname
integer , DIMENSION(max_domains) :: auxinput22_interval_y
integer , DIMENSION(max_domains) :: auxinput22_interval_d
integer , DIMENSION(max_domains) :: auxinput22_interval_h
integer , DIMENSION(max_domains) :: auxinput22_interval_m
integer , DIMENSION(max_domains) :: auxinput22_interval_s
integer , DIMENSION(max_domains) :: auxinput22_interval
integer , DIMENSION(max_domains) :: auxinput22_begin_y
integer , DIMENSION(max_domains) :: auxinput22_begin_d
integer , DIMENSION(max_domains) :: auxinput22_begin_h
integer , DIMENSION(max_domains) :: auxinput22_begin_m
integer , DIMENSION(max_domains) :: auxinput22_begin_s
integer , DIMENSION(max_domains) :: auxinput22_begin
integer , DIMENSION(max_domains) :: auxinput22_end_y
integer , DIMENSION(max_domains) :: auxinput22_end_d
integer , DIMENSION(max_domains) :: auxinput22_end_h
integer , DIMENSION(max_domains) :: auxinput22_end_m
integer , DIMENSION(max_domains) :: auxinput22_end_s
integer , DIMENSION(max_domains) :: auxinput22_end
integer :: io_form_auxinput22
integer , DIMENSION(max_domains) :: frames_per_auxinput22
character*256 :: auxinput23_inname
character*256 :: auxinput23_outname
integer , DIMENSION(max_domains) :: auxinput23_interval_y
integer , DIMENSION(max_domains) :: auxinput23_interval_d
integer , DIMENSION(max_domains) :: auxinput23_interval_h
integer , DIMENSION(max_domains) :: auxinput23_interval_m
integer , DIMENSION(max_domains) :: auxinput23_interval_s
integer , DIMENSION(max_domains) :: auxinput23_interval
integer , DIMENSION(max_domains) :: auxinput23_begin_y
integer , DIMENSION(max_domains) :: auxinput23_begin_d
integer , DIMENSION(max_domains) :: auxinput23_begin_h
integer , DIMENSION(max_domains) :: auxinput23_begin_m
integer , DIMENSION(max_domains) :: auxinput23_begin_s
integer , DIMENSION(max_domains) :: auxinput23_begin
integer , DIMENSION(max_domains) :: auxinput23_end_y
integer , DIMENSION(max_domains) :: auxinput23_end_d
integer , DIMENSION(max_domains) :: auxinput23_end_h
integer , DIMENSION(max_domains) :: auxinput23_end_m
integer , DIMENSION(max_domains) :: auxinput23_end_s
integer , DIMENSION(max_domains) :: auxinput23_end
integer :: io_form_auxinput23
integer , DIMENSION(max_domains) :: frames_per_auxinput23
character*256 :: auxinput24_inname
character*256 :: auxinput24_outname
integer , DIMENSION(max_domains) :: auxinput24_interval_y
integer , DIMENSION(max_domains) :: auxinput24_interval_d
integer , DIMENSION(max_domains) :: auxinput24_interval_h
integer , DIMENSION(max_domains) :: auxinput24_interval_m
integer , DIMENSION(max_domains) :: auxinput24_interval_s
integer , DIMENSION(max_domains) :: auxinput24_interval
integer , DIMENSION(max_domains) :: auxinput24_begin_y
integer , DIMENSION(max_domains) :: auxinput24_begin_d
integer , DIMENSION(max_domains) :: auxinput24_begin_h
integer , DIMENSION(max_domains) :: auxinput24_begin_m
integer , DIMENSION(max_domains) :: auxinput24_begin_s
integer , DIMENSION(max_domains) :: auxinput24_begin
integer , DIMENSION(max_domains) :: auxinput24_end_y
integer , DIMENSION(max_domains) :: auxinput24_end_d
integer , DIMENSION(max_domains) :: auxinput24_end_h
integer , DIMENSION(max_domains) :: auxinput24_end_m
integer , DIMENSION(max_domains) :: auxinput24_end_s
integer , DIMENSION(max_domains) :: auxinput24_end
integer :: io_form_auxinput24
integer , DIMENSION(max_domains) :: frames_per_auxinput24
integer , DIMENSION(max_domains) :: history_interval
integer , DIMENSION(max_domains) :: frames_per_outfile
logical :: restart
integer :: restart_interval
integer :: io_form_input
integer :: io_form_history
integer :: io_form_restart
integer :: io_form_boundary
integer :: debug_level
logical :: self_test_domain
character*256 :: history_outname
character*256 :: history_inname
logical :: use_netcdf_classic
integer , DIMENSION(max_domains) :: history_interval_d
integer , DIMENSION(max_domains) :: history_interval_h
integer , DIMENSION(max_domains) :: history_interval_m
integer , DIMENSION(max_domains) :: history_interval_s
integer , DIMENSION(max_domains) :: inputout_interval_d
integer , DIMENSION(max_domains) :: inputout_interval_h
integer , DIMENSION(max_domains) :: inputout_interval_m
integer , DIMENSION(max_domains) :: inputout_interval_s
integer , DIMENSION(max_domains) :: inputout_interval
integer :: restart_interval_d
integer :: restart_interval_h
integer :: restart_interval_m
integer :: restart_interval_s
integer , DIMENSION(max_domains) :: history_begin_y
integer , DIMENSION(max_domains) :: history_begin_d
integer , DIMENSION(max_domains) :: history_begin_h
integer , DIMENSION(max_domains) :: history_begin_m
integer , DIMENSION(max_domains) :: history_begin_s
integer , DIMENSION(max_domains) :: history_begin
integer , DIMENSION(max_domains) :: inputout_begin_y
integer , DIMENSION(max_domains) :: inputout_begin_d
integer , DIMENSION(max_domains) :: inputout_begin_h
integer , DIMENSION(max_domains) :: inputout_begin_m
integer , DIMENSION(max_domains) :: inputout_begin_s
integer :: restart_begin_y
integer :: restart_begin_d
integer :: restart_begin_h
integer :: restart_begin_m
integer :: restart_begin_s
integer :: restart_begin
integer , DIMENSION(max_domains) :: history_end_y
integer , DIMENSION(max_domains) :: history_end_d
integer , DIMENSION(max_domains) :: history_end_h
integer , DIMENSION(max_domains) :: history_end_m
integer , DIMENSION(max_domains) :: history_end_s
integer , DIMENSION(max_domains) :: history_end
integer , DIMENSION(max_domains) :: inputout_end_y
integer , DIMENSION(max_domains) :: inputout_end_d
integer , DIMENSION(max_domains) :: inputout_end_h
integer , DIMENSION(max_domains) :: inputout_end_m
integer , DIMENSION(max_domains) :: inputout_end_s
integer :: simulation_start_year
integer :: simulation_start_month
integer :: simulation_start_day
integer :: simulation_start_hour
integer :: simulation_start_minute
integer :: simulation_start_second
logical :: reset_simulation_start
integer , DIMENSION(max_domains) :: sr_x
integer , DIMENSION(max_domains) :: sr_y
character*256 , DIMENSION(max_domains) :: iofields_filename
logical :: ignore_iofields_warning
logical :: ncd_nofill
integer , DIMENSION(max_domains) :: julyr
integer , DIMENSION(max_domains) :: julday
real , DIMENSION(max_domains) :: gmt
character*256 :: high_freq_outname
character*256 :: partial_atcf_outname
character*256 :: input_inname
character*256 :: input_outname
character*256 :: bdy_inname
character*256 :: bdy_outname
character*256 :: rst_inname
character*256 :: rst_outname
character*256 , DIMENSION(max_domains) :: anl_outname
logical :: write_input
logical :: write_restart_at_0h
logical :: write_hist_at_0h_rst
logical :: adjust_output_times
logical :: adjust_input_times
real , DIMENSION(max_domains) :: tstart
logical :: nocolons
logical :: cycling
logical :: output_ready_flag
integer :: dfi_opt
integer :: dfi_savehydmeteors
integer :: dfi_nfilter
logical :: dfi_write_filtered_input
logical :: dfi_write_dfi_history
integer :: dfi_cutoff_seconds
integer :: dfi_time_dim
integer :: dfi_fwdstop_year
integer :: dfi_fwdstop_month
integer :: dfi_fwdstop_day
integer :: dfi_fwdstop_hour
integer :: dfi_fwdstop_minute
integer :: dfi_fwdstop_second
integer :: dfi_bckstop_year
integer :: dfi_bckstop_month
integer :: dfi_bckstop_day
integer :: dfi_bckstop_hour
integer :: dfi_bckstop_minute
integer :: dfi_bckstop_second
integer :: time_step
integer :: time_step_fract_num
integer :: time_step_fract_den
integer :: time_step_dfi
integer :: max_dom
integer , DIMENSION(max_domains) :: s_we
integer , DIMENSION(max_domains) :: e_we
integer , DIMENSION(max_domains) :: s_sn
integer , DIMENSION(max_domains) :: e_sn
integer , DIMENSION(max_domains) :: s_vert
integer , DIMENSION(max_domains) :: e_vert
integer :: num_metgrid_soil_levels
real , DIMENSION(max_domains) :: dx
real , DIMENSION(max_domains) :: dy
integer , DIMENSION(max_domains) :: grid_id
logical , DIMENSION(max_domains) :: grid_allowed
integer , DIMENSION(max_domains) :: parent_id
integer , DIMENSION(max_domains) :: i_parent_start
integer , DIMENSION(max_domains) :: j_parent_start
integer , DIMENSION(max_domains) :: parent_grid_ratio
integer , DIMENSION(max_domains) :: parent_time_step_ratio
integer :: feedback
integer :: smooth_option
real , DIMENSION(max_domains) :: ztop
integer , DIMENSION(max_domains) :: moad_grid_ratio
integer , DIMENSION(max_domains) :: moad_time_step_ratio
integer , DIMENSION(max_domains) :: shw
integer :: tile_sz_x
integer :: tile_sz_y
integer :: numtiles
integer :: numtiles_inc
integer :: numtiles_x
integer :: numtiles_y
integer :: tile_strategy
integer :: nproc_x
integer :: nproc_y
integer :: irand
real , DIMENSION(max_domains) :: dt
integer :: ts_buf_size
integer :: max_ts_locs
logical :: ocean_cpl
integer , DIMENSION(max_domains) :: num_moves
integer , DIMENSION(max_domains) :: vortex_interval
integer , DIMENSION(max_domains) :: corral_dist
logical :: multi_storm
logical :: no_ocean
integer , DIMENSION(max_moves) :: move_id
integer , DIMENSION(max_moves) :: move_interval
integer , DIMENSION(max_moves) :: move_cd_x
integer , DIMENSION(max_moves) :: move_cd_y
logical , DIMENSION(max_domains) :: swap_x
logical , DIMENSION(max_domains) :: swap_y
logical , DIMENSION(max_domains) :: cycle_x
logical , DIMENSION(max_domains) :: cycle_y
logical :: reorder_mesh
logical :: perturb_input
real , DIMENSION(max_eta) :: eta_levels
real :: ptsgm
integer :: num_metgrid_levels
real :: p_top_requested
logical :: use_prep_hybrid
logical :: force_read_thompson
logical :: write_thompson_tables
real , DIMENSION(max_domains) :: nssl_cccn
real , DIMENSION(max_domains) :: nssl_alphah
real , DIMENSION(max_domains) :: nssl_alphahl
real , DIMENSION(max_domains) :: nssl_cnoh
real , DIMENSION(max_domains) :: nssl_cnohl
real , DIMENSION(max_domains) :: nssl_cnor
real , DIMENSION(max_domains) :: nssl_cnos
real , DIMENSION(max_domains) :: nssl_rho_qh
real , DIMENSION(max_domains) :: nssl_rho_qhl
real , DIMENSION(max_domains) :: nssl_rho_qs
integer :: elec_physics
integer , DIMENSION(max_domains) :: nssl_ipelec
integer :: nssl_isaund
integer , DIMENSION(max_domains) :: mp_physics
real , DIMENSION(max_domains) :: mommix
logical , DIMENSION(max_domains) :: disheat
integer :: do_radar_ref
integer :: compute_radar_ref
integer , DIMENSION(max_domains) :: ra_lw_physics
integer , DIMENSION(max_domains) :: ra_sw_physics
real , DIMENSION(max_domains) :: radt
integer , DIMENSION(max_domains) :: sf_sfclay_physics
integer , DIMENSION(max_domains) :: sf_surface_physics
integer , DIMENSION(max_domains) :: bl_pbl_physics
integer , DIMENSION(max_domains) :: bl_mynn_tkebudget
integer , DIMENSION(max_domains) :: ysu_topdown_pblmix
integer , DIMENSION(max_domains) :: shinhong_tke_diag
integer :: ysu_firstorder
real :: ysu_brcr_ub
integer :: ysu_cgrad_scalar
integer :: ysu_cgrad_mom
integer , DIMENSION(max_domains) :: windfarm_opt
integer :: windfarm_ij
integer , DIMENSION(max_domains) :: mfshconv
real , DIMENSION(max_domains) :: bldt
integer , DIMENSION(max_domains) :: cu_physics
integer , DIMENSION(max_domains) :: shcu_physics
integer , DIMENSION(max_domains) :: cu_diag
integer , DIMENSION(max_domains) :: kf_edrates
logical , DIMENSION(max_domains) :: bl_mynn_tkeadvect
integer :: bl_mynn_cloudpdf
integer :: bl_mynn_mixlength
integer , DIMENSION(max_domains) :: bl_mynn_edmf
integer , DIMENSION(max_domains) :: bl_mynn_edmf_mom
integer , DIMENSION(max_domains) :: bl_mynn_edmf_tke
integer , DIMENSION(max_domains) :: bl_mynn_edmf_part
integer , DIMENSION(max_domains) :: bl_mynn_cloudmix
integer , DIMENSION(max_domains) :: bl_mynn_mixqt
integer :: icloud_bl
real , DIMENSION(max_domains) :: gfs_alpha
real , DIMENSION(max_domains) :: cudt
real , DIMENSION(max_domains) :: gsmdt
integer :: isfflx
integer :: ideal_xland
integer :: ifsnow
integer :: icloud
integer :: cldovrlp
integer :: idcor
real :: swrad_scat
integer :: surface_input_source
integer :: num_soil_layers
integer :: num_urban_layers
integer :: sf_surface_mosaic
integer :: mosaic_cat
integer :: mosaic_cat_soil
integer :: num_urban_hi
integer :: mosaic_lu
integer :: mosaic_soil
integer :: maxiens
integer :: maxens
integer :: maxens2
integer :: maxens3
integer :: ensdim
integer , DIMENSION(max_domains) :: ishallow
integer :: clos_choice
integer , DIMENSION(max_domains) :: chem_opt
integer :: num_land_cat
integer :: num_soil_cat
integer , DIMENSION(max_domains) :: topo_wind
integer :: mp_zero_out
real :: mp_zero_out_thresh
real :: seaice_threshold
integer :: fractional_seaice
integer :: seaice_albedo_opt
real :: seaice_albedo_default
integer :: seaice_snowdepth_opt
real :: seaice_snowdepth_max
real :: seaice_snowdepth_min
integer :: seaice_thickness_opt
real :: seaice_thickness_default
logical :: tice2tsk_if2cold
integer :: sst_update
integer , DIMENSION(max_domains) :: sf_urban_physics
logical :: usemonalb
logical :: rdmaxalb
logical :: rdlai2d
logical :: ua_phys
integer :: opt_thcnd
integer , DIMENSION(max_domains) :: gwd_opt
integer :: isftcflx
integer :: iz0tlnd
real , DIMENSION(max_domains) :: sas_pgcon
real , DIMENSION(max_domains) :: sas_shal_pgcon
integer , DIMENSION(max_domains) :: sas_shal_conv
real , DIMENSION(max_domains) :: sas_mass_flux
integer , DIMENSION(max_domains) :: grav_settling
real :: var_ric
real :: coef_ric_l
real :: coef_ric_s
integer , DIMENSION(max_domains) :: random_seed
integer , DIMENSION(max_domains) :: icoef_sf
integer , DIMENSION(max_domains) :: iwavecpl
logical , DIMENSION(max_domains) :: lcurr_sf
integer , DIMENSION(max_domains) :: ens_random_seed
integer , DIMENSION(max_domains) :: ensda_physics_pert
logical :: pert_sas
logical :: pert_pbl
logical :: pert_cd
real , DIMENSION(max_domains) :: ens_sasamp
real , DIMENSION(max_domains) :: ens_pblamp
real , DIMENSION(max_domains) :: ens_cdamp
integer , DIMENSION(max_domains) :: idtad
integer , DIMENSION(max_domains) :: nsoil
integer , DIMENSION(max_domains) :: nphs
integer , DIMENSION(max_domains) :: ncnvc
integer , DIMENSION(max_domains) :: nrand
integer , DIMENSION(max_domains) :: nrads
integer , DIMENSION(max_domains) :: nradl
real , DIMENSION(max_domains) :: tprec
real , DIMENSION(max_domains) :: theat
real , DIMENSION(max_domains) :: tclod
real , DIMENSION(max_domains) :: trdsw
real , DIMENSION(max_domains) :: trdlw
real , DIMENSION(max_domains) :: tsrfc
logical , DIMENSION(max_domains) :: pcpflg
integer , DIMENSION(max_domains) :: sigma
real , DIMENSION(max_domains) :: sfenth
integer :: co2tf
integer :: ra_call_offset
real :: cam_abs_freq_s
integer :: levsiz
integer :: paerlev
integer :: cam_abs_dim1
integer :: cam_abs_dim2
integer :: no_src_types
integer :: alevsiz
integer :: o3input
integer :: aer_opt
logical , DIMENSION(max_domains) :: cu_rad_feedback
integer :: aercu_opt
real :: aercu_fct
integer :: icloud_cu
real , DIMENSION(max_domains) :: h_diff
integer , DIMENSION(max_domains) :: ntrack
integer :: num_snso_layers
integer :: num_snow_layers
logical :: use_aero_icbc
integer :: use_mp_re
real :: ccn_conc
integer :: hail_opt
integer :: morr_rimed_ice
integer :: clean_atm_diag
integer :: calc_clean_atm_diag
integer , DIMENSION(max_domains) :: sf_lake_physics
integer :: dyn_opt
integer :: rk_ord
integer :: w_damping
integer , DIMENSION(max_domains) :: diff_opt
integer , DIMENSION(max_domains) :: km_opt
integer :: damp_opt
real , DIMENSION(max_domains) :: zdamp
real :: base_pres
real :: base_temp
real :: base_lapse
real :: iso_temp
real , DIMENSION(max_domains) :: dampcoef
real , DIMENSION(max_domains) :: khdif
real , DIMENSION(max_domains) :: kvdif
real , DIMENSION(max_domains) :: c_s
real , DIMENSION(max_domains) :: c_k
real , DIMENSION(max_domains) :: smdiv
real , DIMENSION(max_domains) :: emdiv
real , DIMENSION(max_domains) :: epssm
integer , DIMENSION(max_domains) :: nhmove
logical , DIMENSION(max_domains) :: non_hydrostatic
integer , DIMENSION(max_domains) :: time_step_sound
integer , DIMENSION(max_domains) :: h_mom_adv_order
integer , DIMENSION(max_domains) :: v_mom_adv_order
integer , DIMENSION(max_domains) :: h_sca_adv_order
integer , DIMENSION(max_domains) :: v_sca_adv_order
logical , DIMENSION(max_domains) :: top_radiation
real , DIMENSION(max_domains) :: tke_upper_bound
real , DIMENSION(max_domains) :: tke_drag_coefficient
real , DIMENSION(max_domains) :: tke_heat_flux
logical , DIMENSION(max_domains) :: pert_coriolis
logical :: euler_adv
integer :: idtadt
integer :: idtadc
real , DIMENSION(max_domains) :: codamp
real , DIMENSION(max_domains) :: coac
real , DIMENSION(max_domains) :: slophc
real , DIMENSION(max_domains) :: wp
real , DIMENSION(max_domains) :: dwdt_damping_lev
integer :: terrain_smoothing
integer :: spec_bdy_width
integer :: spec_zone
integer :: relax_zone
logical , DIMENSION(max_domains) :: specified
logical , DIMENSION(max_domains) :: periodic_x
logical , DIMENSION(max_domains) :: symmetric_xs
logical , DIMENSION(max_domains) :: symmetric_xe
logical , DIMENSION(max_domains) :: open_xs
logical , DIMENSION(max_domains) :: open_xe
logical , DIMENSION(max_domains) :: periodic_y
logical , DIMENSION(max_domains) :: symmetric_ys
logical , DIMENSION(max_domains) :: symmetric_ye
logical , DIMENSION(max_domains) :: open_ys
logical , DIMENSION(max_domains) :: open_ye
logical , DIMENSION(max_domains) :: polar
logical , DIMENSION(max_domains) :: nested
integer :: real_data_init_type
integer :: background_proc_id
integer :: forecast_proc_id
integer :: production_status
integer :: compression
real , DIMENSION(max_domains) :: cen_lat
real , DIMENSION(max_domains) :: cen_lon
real , DIMENSION(max_domains) :: truelat1
real , DIMENSION(max_domains) :: truelat2
real , DIMENSION(max_domains) :: moad_cen_lat
real , DIMENSION(max_domains) :: stand_lon
integer :: flag_metgrid
integer :: flag_snow
integer :: flag_psfc
integer :: flag_sm000010
integer :: flag_sm010040
integer :: flag_sm040100
integer :: flag_sm100200
integer :: flag_st000010
integer :: flag_st010040
integer :: flag_st040100
integer :: flag_st100200
integer :: flag_slp
integer :: flag_soilhgt
integer :: flag_mf_xy
real , DIMENSION(max_domains) :: bdyfrq
character*256 , DIMENSION(max_domains) :: mminlu
integer , DIMENSION(max_domains) :: iswater
integer , DIMENSION(max_domains) :: islake
integer , DIMENSION(max_domains) :: isice
integer , DIMENSION(max_domains) :: isurban
integer , DIMENSION(max_domains) :: isoilwater
integer , DIMENSION(max_domains) :: map_proj
integer :: dfi_stage
integer , DIMENSION(max_domains) :: mp_physics_dfi
integer :: use_wps_input
integer :: maxpatch
integer    :: last_item_in_struct









NAMELIST /physics/ lakedepth_default
NAMELIST /physics/ lake_min_elev
NAMELIST /physics/ use_lakedepth
NAMELIST /stoch/ nens
NAMELIST /stoch/ skebs
NAMELIST /stoch/ stoch_force_opt
NAMELIST /stoch/ skebs_vertstruc
NAMELIST /stoch/ stoch_vertstruc_opt
NAMELIST /stoch/ tot_backscat_psi
NAMELIST /stoch/ tot_backscat_t
NAMELIST /stoch/ ztau_psi
NAMELIST /stoch/ ztau_t
NAMELIST /stoch/ rexponent_psi
NAMELIST /stoch/ rexponent_t
NAMELIST /stoch/ zsigma2_eps
NAMELIST /stoch/ zsigma2_eta
NAMELIST /stoch/ kminforc
NAMELIST /stoch/ lminforc
NAMELIST /stoch/ kminforct
NAMELIST /stoch/ lminforct
NAMELIST /stoch/ kmaxforc
NAMELIST /stoch/ lmaxforc
NAMELIST /stoch/ kmaxforct
NAMELIST /stoch/ lmaxforct
NAMELIST /stoch/ iseed_skebs
NAMELIST /stoch/ sppt
NAMELIST /stoch/ gridpt_stddev_sppt
NAMELIST /stoch/ stddev_cutoff_sppt
NAMELIST /stoch/ lengthscale_sppt
NAMELIST /stoch/ timescale_sppt
NAMELIST /stoch/ sppt_vertstruc
NAMELIST /stoch/ iseed_sppt
NAMELIST /stoch/ spp_mp
NAMELIST /stoch/ gridpt_stddev_spp_mp
NAMELIST /stoch/ stddev_cutoff_spp_mp
NAMELIST /stoch/ lengthscale_spp_mp
NAMELIST /stoch/ timescale_spp_mp
NAMELIST /stoch/ vertstruc_spp_mp
NAMELIST /stoch/ iseed_spp_mp
NAMELIST /stoch/ isppt_mp
NAMELIST /stoch/ gridpt_stddev_isppt_mp
NAMELIST /stoch/ stddev_cutoff_isppt_mp
NAMELIST /stoch/ lengthscale_isppt_mp
NAMELIST /stoch/ timescale_isppt_mp
NAMELIST /stoch/ vertstruc_isppt_mp
NAMELIST /stoch/ iseed_isppt_mp
NAMELIST /stoch/ isppt_pbl
NAMELIST /stoch/ gridpt_stddev_isppt_pbl
NAMELIST /stoch/ stddev_cutoff_isppt_pbl
NAMELIST /stoch/ lengthscale_isppt_pbl
NAMELIST /stoch/ timescale_isppt_pbl
NAMELIST /stoch/ vertstruc_isppt_pbl
NAMELIST /stoch/ iseed_isppt_pbl
NAMELIST /stoch/ isppt_rad
NAMELIST /stoch/ gridpt_stddev_isppt_rad
NAMELIST /stoch/ stddev_cutoff_isppt_rad
NAMELIST /stoch/ lengthscale_isppt_rad
NAMELIST /stoch/ timescale_isppt_rad
NAMELIST /stoch/ vertstruc_isppt_rad
NAMELIST /stoch/ iseed_isppt_rad
NAMELIST /stoch/ isppt_conv
NAMELIST /stoch/ gridpt_stddev_isppt_conv
NAMELIST /stoch/ stddev_cutoff_isppt_conv
NAMELIST /stoch/ lengthscale_isppt_conv
NAMELIST /stoch/ timescale_isppt_conv
NAMELIST /stoch/ vertstruc_isppt_conv
NAMELIST /stoch/ iseed_isppt_conv
NAMELIST /stoch/ rand_perturb
NAMELIST /stoch/ gridpt_stddev_rand_pert
NAMELIST /stoch/ stddev_cutoff_rand_pert
NAMELIST /stoch/ lengthscale_rand_pert
NAMELIST /stoch/ timescale_rand_pert
NAMELIST /stoch/ rand_pert_vertstruc
NAMELIST /stoch/ iseed_rand_pert
NAMELIST /domains/ halo_debug
NAMELIST /physics/ ntracers
NAMELIST /physics/ vortex_tracker
NAMELIST /physics/ interest_rad_storm
NAMELIST /physics/ interest_rad_parent
NAMELIST /physics/ interest_rad_self
NAMELIST /physics/ interest_kids
NAMELIST /physics/ interest_self
NAMELIST /physics/ interest_storms
NAMELIST /physics/ swath_mode
NAMELIST /physics/ num_old_fixes
NAMELIST /physics/ vt4_radius
NAMELIST /physics/ vt4_weightexp
NAMELIST /physics/ vt4_pmax
NAMELIST /physics/ vt4_noise_pmax
NAMELIST /physics/ vt4_noise_pmin
NAMELIST /physics/ vt4_noise_dpdr
NAMELIST /physics/ vt4_noise_iter
NAMELIST /physics/ nomove_freq
NAMELIST /physics/ traj_opt
NAMELIST /physics/ dm_has_traj
NAMELIST /domains/ num_traj
NAMELIST /domains/ corral_x
NAMELIST /domains/ corral_y
NAMELIST /time_control/ tg_reset_stream
NAMELIST /physics/ tg_option
NAMELIST /physics/ ntornado
NAMELIST /time_control/ analysis
NAMELIST /time_control/ write_analysis
NAMELIST /time_control/ io_form_auxinput2
NAMELIST /time_control/ high_freq
NAMELIST /time_control/ high_dom
NAMELIST /physics/ swint_opt
NAMELIST /physics/ aer_type
NAMELIST /physics/ aer_aod550_opt
NAMELIST /physics/ aer_angexp_opt
NAMELIST /physics/ aer_ssa_opt
NAMELIST /physics/ aer_asy_opt
NAMELIST /physics/ aer_aod550_val
NAMELIST /physics/ aer_angexp_val
NAMELIST /physics/ aer_ssa_val
NAMELIST /physics/ aer_asy_val
NAMELIST /noah_mp/ dveg
NAMELIST /noah_mp/ opt_crs
NAMELIST /noah_mp/ opt_btr
NAMELIST /noah_mp/ opt_run
NAMELIST /noah_mp/ opt_sfc
NAMELIST /noah_mp/ opt_frz
NAMELIST /noah_mp/ opt_inf
NAMELIST /noah_mp/ opt_rad
NAMELIST /noah_mp/ opt_alb
NAMELIST /noah_mp/ opt_snf
NAMELIST /noah_mp/ opt_tbot
NAMELIST /noah_mp/ opt_stc
NAMELIST /noah_mp/ opt_gla
NAMELIST /noah_mp/ opt_rsf
NAMELIST /noah_mp/ opt_soil
NAMELIST /noah_mp/ opt_pedo
NAMELIST /noah_mp/ opt_crop
NAMELIST /time_control/ run_days
NAMELIST /time_control/ run_hours
NAMELIST /time_control/ run_minutes
NAMELIST /time_control/ run_seconds
NAMELIST /time_control/ start_year
NAMELIST /time_control/ start_month
NAMELIST /time_control/ start_day
NAMELIST /time_control/ start_hour
NAMELIST /time_control/ start_minute
NAMELIST /time_control/ start_second
NAMELIST /time_control/ end_year
NAMELIST /time_control/ end_month
NAMELIST /time_control/ end_day
NAMELIST /time_control/ end_hour
NAMELIST /time_control/ end_minute
NAMELIST /time_control/ end_second
NAMELIST /time_control/ interval_seconds
NAMELIST /time_control/ input_from_file
NAMELIST /time_control/ fine_input_stream
NAMELIST /time_control/ auxinput1_inname
NAMELIST /time_control/ io_form_auxinput1
NAMELIST /time_control/ override_restart_timers
NAMELIST /time_control/ auxhist1_inname
NAMELIST /time_control/ auxhist1_outname
NAMELIST /time_control/ auxhist1_interval_y
NAMELIST /time_control/ auxhist1_interval_d
NAMELIST /time_control/ auxhist1_interval_h
NAMELIST /time_control/ auxhist1_interval_m
NAMELIST /time_control/ auxhist1_interval_s
NAMELIST /time_control/ auxhist1_interval
NAMELIST /time_control/ auxhist1_begin_y
NAMELIST /time_control/ auxhist1_begin_d
NAMELIST /time_control/ auxhist1_begin_h
NAMELIST /time_control/ auxhist1_begin_m
NAMELIST /time_control/ auxhist1_begin_s
NAMELIST /time_control/ auxhist1_begin
NAMELIST /time_control/ auxhist1_end_y
NAMELIST /time_control/ auxhist1_end_d
NAMELIST /time_control/ auxhist1_end_h
NAMELIST /time_control/ auxhist1_end_m
NAMELIST /time_control/ auxhist1_end_s
NAMELIST /time_control/ auxhist1_end
NAMELIST /time_control/ io_form_auxhist1
NAMELIST /time_control/ frames_per_auxhist1
NAMELIST /time_control/ auxhist2_inname
NAMELIST /time_control/ auxhist2_outname
NAMELIST /time_control/ auxhist2_interval_y
NAMELIST /time_control/ auxhist2_interval_d
NAMELIST /time_control/ auxhist2_interval_h
NAMELIST /time_control/ auxhist2_interval_m
NAMELIST /time_control/ auxhist2_interval_s
NAMELIST /time_control/ auxhist2_interval
NAMELIST /time_control/ auxhist2_begin_y
NAMELIST /time_control/ auxhist2_begin_d
NAMELIST /time_control/ auxhist2_begin_h
NAMELIST /time_control/ auxhist2_begin_m
NAMELIST /time_control/ auxhist2_begin_s
NAMELIST /time_control/ auxhist2_begin
NAMELIST /time_control/ auxhist2_end_y
NAMELIST /time_control/ auxhist2_end_d
NAMELIST /time_control/ auxhist2_end_h
NAMELIST /time_control/ auxhist2_end_m
NAMELIST /time_control/ auxhist2_end_s
NAMELIST /time_control/ auxhist2_end
NAMELIST /time_control/ io_form_auxhist2
NAMELIST /time_control/ frames_per_auxhist2
NAMELIST /time_control/ auxhist3_inname
NAMELIST /time_control/ auxhist3_outname
NAMELIST /time_control/ auxhist3_interval_y
NAMELIST /time_control/ auxhist3_interval_d
NAMELIST /time_control/ auxhist3_interval_h
NAMELIST /time_control/ auxhist3_interval_m
NAMELIST /time_control/ auxhist3_interval_s
NAMELIST /time_control/ auxhist3_interval
NAMELIST /time_control/ auxhist3_begin_y
NAMELIST /time_control/ auxhist3_begin_d
NAMELIST /time_control/ auxhist3_begin_h
NAMELIST /time_control/ auxhist3_begin_m
NAMELIST /time_control/ auxhist3_begin_s
NAMELIST /time_control/ auxhist3_begin
NAMELIST /time_control/ auxhist3_end_y
NAMELIST /time_control/ auxhist3_end_d
NAMELIST /time_control/ auxhist3_end_h
NAMELIST /time_control/ auxhist3_end_m
NAMELIST /time_control/ auxhist3_end_s
NAMELIST /time_control/ auxhist3_end
NAMELIST /time_control/ io_form_auxhist3
NAMELIST /time_control/ frames_per_auxhist3
NAMELIST /time_control/ auxhist4_inname
NAMELIST /time_control/ auxhist4_outname
NAMELIST /time_control/ auxhist4_interval_y
NAMELIST /time_control/ auxhist4_interval_d
NAMELIST /time_control/ auxhist4_interval_h
NAMELIST /time_control/ auxhist4_interval_m
NAMELIST /time_control/ auxhist4_interval_s
NAMELIST /time_control/ auxhist4_interval
NAMELIST /time_control/ auxhist4_begin_y
NAMELIST /time_control/ auxhist4_begin_d
NAMELIST /time_control/ auxhist4_begin_h
NAMELIST /time_control/ auxhist4_begin_m
NAMELIST /time_control/ auxhist4_begin_s
NAMELIST /time_control/ auxhist4_begin
NAMELIST /time_control/ auxhist4_end_y
NAMELIST /time_control/ auxhist4_end_d
NAMELIST /time_control/ auxhist4_end_h
NAMELIST /time_control/ auxhist4_end_m
NAMELIST /time_control/ auxhist4_end_s
NAMELIST /time_control/ auxhist4_end
NAMELIST /time_control/ io_form_auxhist4
NAMELIST /time_control/ frames_per_auxhist4
NAMELIST /time_control/ auxhist5_inname
NAMELIST /time_control/ auxhist5_outname
NAMELIST /time_control/ auxhist5_interval_y
NAMELIST /time_control/ auxhist5_interval_d
NAMELIST /time_control/ auxhist5_interval_h
NAMELIST /time_control/ auxhist5_interval_m
NAMELIST /time_control/ auxhist5_interval_s
NAMELIST /time_control/ auxhist5_interval
NAMELIST /time_control/ auxhist5_begin_y
NAMELIST /time_control/ auxhist5_begin_d
NAMELIST /time_control/ auxhist5_begin_h
NAMELIST /time_control/ auxhist5_begin_m
NAMELIST /time_control/ auxhist5_begin_s
NAMELIST /time_control/ auxhist5_begin
NAMELIST /time_control/ auxhist5_end_y
NAMELIST /time_control/ auxhist5_end_d
NAMELIST /time_control/ auxhist5_end_h
NAMELIST /time_control/ auxhist5_end_m
NAMELIST /time_control/ auxhist5_end_s
NAMELIST /time_control/ auxhist5_end
NAMELIST /time_control/ io_form_auxhist5
NAMELIST /time_control/ frames_per_auxhist5
NAMELIST /time_control/ auxhist6_inname
NAMELIST /time_control/ auxhist6_outname
NAMELIST /time_control/ auxhist6_interval_y
NAMELIST /time_control/ auxhist6_interval_d
NAMELIST /time_control/ auxhist6_interval_h
NAMELIST /time_control/ auxhist6_interval_m
NAMELIST /time_control/ auxhist6_interval_s
NAMELIST /time_control/ auxhist6_interval
NAMELIST /time_control/ auxhist6_begin_y
NAMELIST /time_control/ auxhist6_begin_d
NAMELIST /time_control/ auxhist6_begin_h
NAMELIST /time_control/ auxhist6_begin_m
NAMELIST /time_control/ auxhist6_begin_s
NAMELIST /time_control/ auxhist6_begin
NAMELIST /time_control/ auxhist6_end_y
NAMELIST /time_control/ auxhist6_end_d
NAMELIST /time_control/ auxhist6_end_h
NAMELIST /time_control/ auxhist6_end_m
NAMELIST /time_control/ auxhist6_end_s
NAMELIST /time_control/ auxhist6_end
NAMELIST /time_control/ io_form_auxhist6
NAMELIST /time_control/ frames_per_auxhist6
NAMELIST /time_control/ auxhist7_inname
NAMELIST /time_control/ auxhist7_outname
NAMELIST /time_control/ auxhist7_interval_y
NAMELIST /time_control/ auxhist7_interval_d
NAMELIST /time_control/ auxhist7_interval_h
NAMELIST /time_control/ auxhist7_interval_m
NAMELIST /time_control/ auxhist7_interval_s
NAMELIST /time_control/ auxhist7_interval
NAMELIST /time_control/ auxhist7_begin_y
NAMELIST /time_control/ auxhist7_begin_d
NAMELIST /time_control/ auxhist7_begin_h
NAMELIST /time_control/ auxhist7_begin_m
NAMELIST /time_control/ auxhist7_begin_s
NAMELIST /time_control/ auxhist7_begin
NAMELIST /time_control/ auxhist7_end_y
NAMELIST /time_control/ auxhist7_end_d
NAMELIST /time_control/ auxhist7_end_h
NAMELIST /time_control/ auxhist7_end_m
NAMELIST /time_control/ auxhist7_end_s
NAMELIST /time_control/ auxhist7_end
NAMELIST /time_control/ io_form_auxhist7
NAMELIST /time_control/ frames_per_auxhist7
NAMELIST /time_control/ auxhist8_inname
NAMELIST /time_control/ auxhist8_outname
NAMELIST /time_control/ auxhist8_interval_y
NAMELIST /time_control/ auxhist8_interval_d
NAMELIST /time_control/ auxhist8_interval_h
NAMELIST /time_control/ auxhist8_interval_m
NAMELIST /time_control/ auxhist8_interval_s
NAMELIST /time_control/ auxhist8_interval
NAMELIST /time_control/ auxhist8_begin_y
NAMELIST /time_control/ auxhist8_begin_d
NAMELIST /time_control/ auxhist8_begin_h
NAMELIST /time_control/ auxhist8_begin_m
NAMELIST /time_control/ auxhist8_begin_s
NAMELIST /time_control/ auxhist8_begin
NAMELIST /time_control/ auxhist8_end_y
NAMELIST /time_control/ auxhist8_end_d
NAMELIST /time_control/ auxhist8_end_h
NAMELIST /time_control/ auxhist8_end_m
NAMELIST /time_control/ auxhist8_end_s
NAMELIST /time_control/ auxhist8_end
NAMELIST /time_control/ io_form_auxhist8
NAMELIST /time_control/ frames_per_auxhist8
NAMELIST /time_control/ auxhist9_inname
NAMELIST /time_control/ auxhist9_outname
NAMELIST /time_control/ auxhist9_interval_y
NAMELIST /time_control/ auxhist9_interval_d
NAMELIST /time_control/ auxhist9_interval_h
NAMELIST /time_control/ auxhist9_interval_m
NAMELIST /time_control/ auxhist9_interval_s
NAMELIST /time_control/ auxhist9_interval
NAMELIST /time_control/ auxhist9_begin_y
NAMELIST /time_control/ auxhist9_begin_d
NAMELIST /time_control/ auxhist9_begin_h
NAMELIST /time_control/ auxhist9_begin_m
NAMELIST /time_control/ auxhist9_begin_s
NAMELIST /time_control/ auxhist9_begin
NAMELIST /time_control/ auxhist9_end_y
NAMELIST /time_control/ auxhist9_end_d
NAMELIST /time_control/ auxhist9_end_h
NAMELIST /time_control/ auxhist9_end_m
NAMELIST /time_control/ auxhist9_end_s
NAMELIST /time_control/ auxhist9_end
NAMELIST /time_control/ io_form_auxhist9
NAMELIST /time_control/ frames_per_auxhist9
NAMELIST /time_control/ auxhist10_inname
NAMELIST /time_control/ auxhist10_outname
NAMELIST /time_control/ auxhist10_interval_y
NAMELIST /time_control/ auxhist10_interval_d
NAMELIST /time_control/ auxhist10_interval_h
NAMELIST /time_control/ auxhist10_interval_m
NAMELIST /time_control/ auxhist10_interval_s
NAMELIST /time_control/ auxhist10_interval
NAMELIST /time_control/ auxhist10_begin_y
NAMELIST /time_control/ auxhist10_begin_d
NAMELIST /time_control/ auxhist10_begin_h
NAMELIST /time_control/ auxhist10_begin_m
NAMELIST /time_control/ auxhist10_begin_s
NAMELIST /time_control/ auxhist10_begin
NAMELIST /time_control/ auxhist10_end_y
NAMELIST /time_control/ auxhist10_end_d
NAMELIST /time_control/ auxhist10_end_h
NAMELIST /time_control/ auxhist10_end_m
NAMELIST /time_control/ auxhist10_end_s
NAMELIST /time_control/ auxhist10_end
NAMELIST /time_control/ io_form_auxhist10
NAMELIST /time_control/ frames_per_auxhist10
NAMELIST /time_control/ auxhist11_inname
NAMELIST /time_control/ auxhist11_outname
NAMELIST /time_control/ auxhist11_interval_y
NAMELIST /time_control/ auxhist11_interval_d
NAMELIST /time_control/ auxhist11_interval_h
NAMELIST /time_control/ auxhist11_interval_m
NAMELIST /time_control/ auxhist11_interval_s
NAMELIST /time_control/ auxhist11_interval
NAMELIST /time_control/ auxhist11_begin_y
NAMELIST /time_control/ auxhist11_begin_d
NAMELIST /time_control/ auxhist11_begin_h
NAMELIST /time_control/ auxhist11_begin_m
NAMELIST /time_control/ auxhist11_begin_s
NAMELIST /time_control/ auxhist11_begin
NAMELIST /time_control/ auxhist11_end_y
NAMELIST /time_control/ auxhist11_end_d
NAMELIST /time_control/ auxhist11_end_h
NAMELIST /time_control/ auxhist11_end_m
NAMELIST /time_control/ auxhist11_end_s
NAMELIST /time_control/ auxhist11_end
NAMELIST /time_control/ io_form_auxhist11
NAMELIST /time_control/ frames_per_auxhist11
NAMELIST /time_control/ auxhist12_inname
NAMELIST /time_control/ auxhist12_outname
NAMELIST /time_control/ auxhist12_interval_y
NAMELIST /time_control/ auxhist12_interval_d
NAMELIST /time_control/ auxhist12_interval_h
NAMELIST /time_control/ auxhist12_interval_m
NAMELIST /time_control/ auxhist12_interval_s
NAMELIST /time_control/ auxhist12_interval
NAMELIST /time_control/ auxhist12_begin_y
NAMELIST /time_control/ auxhist12_begin_d
NAMELIST /time_control/ auxhist12_begin_h
NAMELIST /time_control/ auxhist12_begin_m
NAMELIST /time_control/ auxhist12_begin_s
NAMELIST /time_control/ auxhist12_begin
NAMELIST /time_control/ auxhist12_end_y
NAMELIST /time_control/ auxhist12_end_d
NAMELIST /time_control/ auxhist12_end_h
NAMELIST /time_control/ auxhist12_end_m
NAMELIST /time_control/ auxhist12_end_s
NAMELIST /time_control/ auxhist12_end
NAMELIST /time_control/ io_form_auxhist12
NAMELIST /time_control/ frames_per_auxhist12
NAMELIST /time_control/ auxhist13_inname
NAMELIST /time_control/ auxhist13_outname
NAMELIST /time_control/ auxhist13_interval_y
NAMELIST /time_control/ auxhist13_interval_d
NAMELIST /time_control/ auxhist13_interval_h
NAMELIST /time_control/ auxhist13_interval_m
NAMELIST /time_control/ auxhist13_interval_s
NAMELIST /time_control/ auxhist13_interval
NAMELIST /time_control/ auxhist13_begin_y
NAMELIST /time_control/ auxhist13_begin_d
NAMELIST /time_control/ auxhist13_begin_h
NAMELIST /time_control/ auxhist13_begin_m
NAMELIST /time_control/ auxhist13_begin_s
NAMELIST /time_control/ auxhist13_begin
NAMELIST /time_control/ auxhist13_end_y
NAMELIST /time_control/ auxhist13_end_d
NAMELIST /time_control/ auxhist13_end_h
NAMELIST /time_control/ auxhist13_end_m
NAMELIST /time_control/ auxhist13_end_s
NAMELIST /time_control/ auxhist13_end
NAMELIST /time_control/ io_form_auxhist13
NAMELIST /time_control/ frames_per_auxhist13
NAMELIST /time_control/ auxhist14_inname
NAMELIST /time_control/ auxhist14_outname
NAMELIST /time_control/ auxhist14_interval_y
NAMELIST /time_control/ auxhist14_interval_d
NAMELIST /time_control/ auxhist14_interval_h
NAMELIST /time_control/ auxhist14_interval_m
NAMELIST /time_control/ auxhist14_interval_s
NAMELIST /time_control/ auxhist14_interval
NAMELIST /time_control/ auxhist14_begin_y
NAMELIST /time_control/ auxhist14_begin_d
NAMELIST /time_control/ auxhist14_begin_h
NAMELIST /time_control/ auxhist14_begin_m
NAMELIST /time_control/ auxhist14_begin_s
NAMELIST /time_control/ auxhist14_begin
NAMELIST /time_control/ auxhist14_end_y
NAMELIST /time_control/ auxhist14_end_d
NAMELIST /time_control/ auxhist14_end_h
NAMELIST /time_control/ auxhist14_end_m
NAMELIST /time_control/ auxhist14_end_s
NAMELIST /time_control/ auxhist14_end
NAMELIST /time_control/ io_form_auxhist14
NAMELIST /time_control/ frames_per_auxhist14
NAMELIST /time_control/ auxhist15_inname
NAMELIST /time_control/ auxhist15_outname
NAMELIST /time_control/ auxhist15_interval_y
NAMELIST /time_control/ auxhist15_interval_d
NAMELIST /time_control/ auxhist15_interval_h
NAMELIST /time_control/ auxhist15_interval_m
NAMELIST /time_control/ auxhist15_interval_s
NAMELIST /time_control/ auxhist15_interval
NAMELIST /time_control/ auxhist15_begin_y
NAMELIST /time_control/ auxhist15_begin_d
NAMELIST /time_control/ auxhist15_begin_h
NAMELIST /time_control/ auxhist15_begin_m
NAMELIST /time_control/ auxhist15_begin_s
NAMELIST /time_control/ auxhist15_begin
NAMELIST /time_control/ auxhist15_end_y
NAMELIST /time_control/ auxhist15_end_d
NAMELIST /time_control/ auxhist15_end_h
NAMELIST /time_control/ auxhist15_end_m
NAMELIST /time_control/ auxhist15_end_s
NAMELIST /time_control/ auxhist15_end
NAMELIST /time_control/ io_form_auxhist15
NAMELIST /time_control/ frames_per_auxhist15
NAMELIST /time_control/ auxhist16_inname
NAMELIST /time_control/ auxhist16_outname
NAMELIST /time_control/ auxhist16_interval_y
NAMELIST /time_control/ auxhist16_interval_d
NAMELIST /time_control/ auxhist16_interval_h
NAMELIST /time_control/ auxhist16_interval_m
NAMELIST /time_control/ auxhist16_interval_s
NAMELIST /time_control/ auxhist16_interval
NAMELIST /time_control/ auxhist16_begin_y
NAMELIST /time_control/ auxhist16_begin_d
NAMELIST /time_control/ auxhist16_begin_h
NAMELIST /time_control/ auxhist16_begin_m
NAMELIST /time_control/ auxhist16_begin_s
NAMELIST /time_control/ auxhist16_begin
NAMELIST /time_control/ auxhist16_end_y
NAMELIST /time_control/ auxhist16_end_d
NAMELIST /time_control/ auxhist16_end_h
NAMELIST /time_control/ auxhist16_end_m
NAMELIST /time_control/ auxhist16_end_s
NAMELIST /time_control/ auxhist16_end
NAMELIST /time_control/ io_form_auxhist16
NAMELIST /time_control/ frames_per_auxhist16
NAMELIST /time_control/ auxhist17_inname
NAMELIST /time_control/ auxhist17_outname
NAMELIST /time_control/ auxhist17_interval_y
NAMELIST /time_control/ auxhist17_interval_d
NAMELIST /time_control/ auxhist17_interval_h
NAMELIST /time_control/ auxhist17_interval_m
NAMELIST /time_control/ auxhist17_interval_s
NAMELIST /time_control/ auxhist17_interval
NAMELIST /time_control/ auxhist17_begin_y
NAMELIST /time_control/ auxhist17_begin_d
NAMELIST /time_control/ auxhist17_begin_h
NAMELIST /time_control/ auxhist17_begin_m
NAMELIST /time_control/ auxhist17_begin_s
NAMELIST /time_control/ auxhist17_begin
NAMELIST /time_control/ auxhist17_end_y
NAMELIST /time_control/ auxhist17_end_d
NAMELIST /time_control/ auxhist17_end_h
NAMELIST /time_control/ auxhist17_end_m
NAMELIST /time_control/ auxhist17_end_s
NAMELIST /time_control/ auxhist17_end
NAMELIST /time_control/ io_form_auxhist17
NAMELIST /time_control/ frames_per_auxhist17
NAMELIST /time_control/ auxhist18_inname
NAMELIST /time_control/ auxhist18_outname
NAMELIST /time_control/ auxhist18_interval_y
NAMELIST /time_control/ auxhist18_interval_d
NAMELIST /time_control/ auxhist18_interval_h
NAMELIST /time_control/ auxhist18_interval_m
NAMELIST /time_control/ auxhist18_interval_s
NAMELIST /time_control/ auxhist18_interval
NAMELIST /time_control/ auxhist18_begin_y
NAMELIST /time_control/ auxhist18_begin_d
NAMELIST /time_control/ auxhist18_begin_h
NAMELIST /time_control/ auxhist18_begin_m
NAMELIST /time_control/ auxhist18_begin_s
NAMELIST /time_control/ auxhist18_begin
NAMELIST /time_control/ auxhist18_end_y
NAMELIST /time_control/ auxhist18_end_d
NAMELIST /time_control/ auxhist18_end_h
NAMELIST /time_control/ auxhist18_end_m
NAMELIST /time_control/ auxhist18_end_s
NAMELIST /time_control/ auxhist18_end
NAMELIST /time_control/ io_form_auxhist18
NAMELIST /time_control/ frames_per_auxhist18
NAMELIST /time_control/ auxhist19_inname
NAMELIST /time_control/ auxhist19_outname
NAMELIST /time_control/ auxhist19_interval_y
NAMELIST /time_control/ auxhist19_interval_d
NAMELIST /time_control/ auxhist19_interval_h
NAMELIST /time_control/ auxhist19_interval_m
NAMELIST /time_control/ auxhist19_interval_s
NAMELIST /time_control/ auxhist19_interval
NAMELIST /time_control/ auxhist19_begin_y
NAMELIST /time_control/ auxhist19_begin_d
NAMELIST /time_control/ auxhist19_begin_h
NAMELIST /time_control/ auxhist19_begin_m
NAMELIST /time_control/ auxhist19_begin_s
NAMELIST /time_control/ auxhist19_begin
NAMELIST /time_control/ auxhist19_end_y
NAMELIST /time_control/ auxhist19_end_d
NAMELIST /time_control/ auxhist19_end_h
NAMELIST /time_control/ auxhist19_end_m
NAMELIST /time_control/ auxhist19_end_s
NAMELIST /time_control/ auxhist19_end
NAMELIST /time_control/ io_form_auxhist19
NAMELIST /time_control/ frames_per_auxhist19
NAMELIST /time_control/ auxhist20_inname
NAMELIST /time_control/ auxhist20_outname
NAMELIST /time_control/ auxhist20_interval_y
NAMELIST /time_control/ auxhist20_interval_d
NAMELIST /time_control/ auxhist20_interval_h
NAMELIST /time_control/ auxhist20_interval_m
NAMELIST /time_control/ auxhist20_interval_s
NAMELIST /time_control/ auxhist20_interval
NAMELIST /time_control/ auxhist20_begin_y
NAMELIST /time_control/ auxhist20_begin_d
NAMELIST /time_control/ auxhist20_begin_h
NAMELIST /time_control/ auxhist20_begin_m
NAMELIST /time_control/ auxhist20_begin_s
NAMELIST /time_control/ auxhist20_begin
NAMELIST /time_control/ auxhist20_end_y
NAMELIST /time_control/ auxhist20_end_d
NAMELIST /time_control/ auxhist20_end_h
NAMELIST /time_control/ auxhist20_end_m
NAMELIST /time_control/ auxhist20_end_s
NAMELIST /time_control/ auxhist20_end
NAMELIST /time_control/ io_form_auxhist20
NAMELIST /time_control/ frames_per_auxhist20
NAMELIST /time_control/ auxhist21_inname
NAMELIST /time_control/ auxhist21_outname
NAMELIST /time_control/ auxhist21_interval_y
NAMELIST /time_control/ auxhist21_interval_d
NAMELIST /time_control/ auxhist21_interval_h
NAMELIST /time_control/ auxhist21_interval_m
NAMELIST /time_control/ auxhist21_interval_s
NAMELIST /time_control/ auxhist21_interval
NAMELIST /time_control/ auxhist21_begin_y
NAMELIST /time_control/ auxhist21_begin_d
NAMELIST /time_control/ auxhist21_begin_h
NAMELIST /time_control/ auxhist21_begin_m
NAMELIST /time_control/ auxhist21_begin_s
NAMELIST /time_control/ auxhist21_begin
NAMELIST /time_control/ auxhist21_end_y
NAMELIST /time_control/ auxhist21_end_d
NAMELIST /time_control/ auxhist21_end_h
NAMELIST /time_control/ auxhist21_end_m
NAMELIST /time_control/ auxhist21_end_s
NAMELIST /time_control/ auxhist21_end
NAMELIST /time_control/ io_form_auxhist21
NAMELIST /time_control/ frames_per_auxhist21
NAMELIST /time_control/ auxhist22_inname
NAMELIST /time_control/ auxhist22_outname
NAMELIST /time_control/ auxhist22_interval_y
NAMELIST /time_control/ auxhist22_interval_d
NAMELIST /time_control/ auxhist22_interval_h
NAMELIST /time_control/ auxhist22_interval_m
NAMELIST /time_control/ auxhist22_interval_s
NAMELIST /time_control/ auxhist22_interval
NAMELIST /time_control/ auxhist22_begin_y
NAMELIST /time_control/ auxhist22_begin_d
NAMELIST /time_control/ auxhist22_begin_h
NAMELIST /time_control/ auxhist22_begin_m
NAMELIST /time_control/ auxhist22_begin_s
NAMELIST /time_control/ auxhist22_begin
NAMELIST /time_control/ auxhist22_end_y
NAMELIST /time_control/ auxhist22_end_d
NAMELIST /time_control/ auxhist22_end_h
NAMELIST /time_control/ auxhist22_end_m
NAMELIST /time_control/ auxhist22_end_s
NAMELIST /time_control/ auxhist22_end
NAMELIST /time_control/ io_form_auxhist22
NAMELIST /time_control/ frames_per_auxhist22
NAMELIST /time_control/ auxhist23_inname
NAMELIST /time_control/ auxhist23_outname
NAMELIST /time_control/ auxhist23_interval_y
NAMELIST /time_control/ auxhist23_interval_d
NAMELIST /time_control/ auxhist23_interval_h
NAMELIST /time_control/ auxhist23_interval_m
NAMELIST /time_control/ auxhist23_interval_s
NAMELIST /time_control/ auxhist23_interval
NAMELIST /time_control/ auxhist23_begin_y
NAMELIST /time_control/ auxhist23_begin_d
NAMELIST /time_control/ auxhist23_begin_h
NAMELIST /time_control/ auxhist23_begin_m
NAMELIST /time_control/ auxhist23_begin_s
NAMELIST /time_control/ auxhist23_begin
NAMELIST /time_control/ auxhist23_end_y
NAMELIST /time_control/ auxhist23_end_d
NAMELIST /time_control/ auxhist23_end_h
NAMELIST /time_control/ auxhist23_end_m
NAMELIST /time_control/ auxhist23_end_s
NAMELIST /time_control/ auxhist23_end
NAMELIST /time_control/ io_form_auxhist23
NAMELIST /time_control/ frames_per_auxhist23
NAMELIST /time_control/ auxhist24_inname
NAMELIST /time_control/ auxhist24_outname
NAMELIST /time_control/ auxhist24_interval_y
NAMELIST /time_control/ auxhist24_interval_d
NAMELIST /time_control/ auxhist24_interval_h
NAMELIST /time_control/ auxhist24_interval_m
NAMELIST /time_control/ auxhist24_interval_s
NAMELIST /time_control/ auxhist24_interval
NAMELIST /time_control/ auxhist24_begin_y
NAMELIST /time_control/ auxhist24_begin_d
NAMELIST /time_control/ auxhist24_begin_h
NAMELIST /time_control/ auxhist24_begin_m
NAMELIST /time_control/ auxhist24_begin_s
NAMELIST /time_control/ auxhist24_begin
NAMELIST /time_control/ auxhist24_end_y
NAMELIST /time_control/ auxhist24_end_d
NAMELIST /time_control/ auxhist24_end_h
NAMELIST /time_control/ auxhist24_end_m
NAMELIST /time_control/ auxhist24_end_s
NAMELIST /time_control/ auxhist24_end
NAMELIST /time_control/ io_form_auxhist24
NAMELIST /time_control/ frames_per_auxhist24
NAMELIST /time_control/ auxinput1_outname
NAMELIST /time_control/ auxinput1_interval_y
NAMELIST /time_control/ auxinput1_interval_d
NAMELIST /time_control/ auxinput1_interval_h
NAMELIST /time_control/ auxinput1_interval_m
NAMELIST /time_control/ auxinput1_interval_s
NAMELIST /time_control/ auxinput1_interval
NAMELIST /time_control/ auxinput1_begin_y
NAMELIST /time_control/ auxinput1_begin_d
NAMELIST /time_control/ auxinput1_begin_h
NAMELIST /time_control/ auxinput1_begin_m
NAMELIST /time_control/ auxinput1_begin_s
NAMELIST /time_control/ auxinput1_begin
NAMELIST /time_control/ auxinput1_end_y
NAMELIST /time_control/ auxinput1_end_d
NAMELIST /time_control/ auxinput1_end_h
NAMELIST /time_control/ auxinput1_end_m
NAMELIST /time_control/ auxinput1_end_s
NAMELIST /time_control/ auxinput1_end
NAMELIST /time_control/ frames_per_auxinput1
NAMELIST /time_control/ auxinput2_inname
NAMELIST /time_control/ auxinput2_outname
NAMELIST /time_control/ auxinput2_interval_y
NAMELIST /time_control/ auxinput2_interval_d
NAMELIST /time_control/ auxinput2_interval_h
NAMELIST /time_control/ auxinput2_interval_m
NAMELIST /time_control/ auxinput2_interval_s
NAMELIST /time_control/ auxinput2_interval
NAMELIST /time_control/ auxinput2_begin_y
NAMELIST /time_control/ auxinput2_begin_d
NAMELIST /time_control/ auxinput2_begin_h
NAMELIST /time_control/ auxinput2_begin_m
NAMELIST /time_control/ auxinput2_begin_s
NAMELIST /time_control/ auxinput2_begin
NAMELIST /time_control/ auxinput2_end_y
NAMELIST /time_control/ auxinput2_end_d
NAMELIST /time_control/ auxinput2_end_h
NAMELIST /time_control/ auxinput2_end_m
NAMELIST /time_control/ auxinput2_end_s
NAMELIST /time_control/ auxinput2_end
NAMELIST /time_control/ frames_per_auxinput2
NAMELIST /time_control/ auxinput3_inname
NAMELIST /time_control/ auxinput3_outname
NAMELIST /time_control/ auxinput3_interval_y
NAMELIST /time_control/ auxinput3_interval_d
NAMELIST /time_control/ auxinput3_interval_h
NAMELIST /time_control/ auxinput3_interval_m
NAMELIST /time_control/ auxinput3_interval_s
NAMELIST /time_control/ auxinput3_interval
NAMELIST /time_control/ auxinput3_begin_y
NAMELIST /time_control/ auxinput3_begin_d
NAMELIST /time_control/ auxinput3_begin_h
NAMELIST /time_control/ auxinput3_begin_m
NAMELIST /time_control/ auxinput3_begin_s
NAMELIST /time_control/ auxinput3_begin
NAMELIST /time_control/ auxinput3_end_y
NAMELIST /time_control/ auxinput3_end_d
NAMELIST /time_control/ auxinput3_end_h
NAMELIST /time_control/ auxinput3_end_m
NAMELIST /time_control/ auxinput3_end_s
NAMELIST /time_control/ auxinput3_end
NAMELIST /time_control/ io_form_auxinput3
NAMELIST /time_control/ frames_per_auxinput3
NAMELIST /time_control/ auxinput4_inname
NAMELIST /time_control/ auxinput4_outname
NAMELIST /time_control/ auxinput4_interval_y
NAMELIST /time_control/ auxinput4_interval_d
NAMELIST /time_control/ auxinput4_interval_h
NAMELIST /time_control/ auxinput4_interval_m
NAMELIST /time_control/ auxinput4_interval_s
NAMELIST /time_control/ auxinput4_interval
NAMELIST /time_control/ auxinput4_begin_y
NAMELIST /time_control/ auxinput4_begin_d
NAMELIST /time_control/ auxinput4_begin_h
NAMELIST /time_control/ auxinput4_begin_m
NAMELIST /time_control/ auxinput4_begin_s
NAMELIST /time_control/ auxinput4_begin
NAMELIST /time_control/ auxinput4_end_y
NAMELIST /time_control/ auxinput4_end_d
NAMELIST /time_control/ auxinput4_end_h
NAMELIST /time_control/ auxinput4_end_m
NAMELIST /time_control/ auxinput4_end_s
NAMELIST /time_control/ auxinput4_end
NAMELIST /time_control/ io_form_auxinput4
NAMELIST /time_control/ frames_per_auxinput4
NAMELIST /time_control/ auxinput5_inname
NAMELIST /time_control/ auxinput5_outname
NAMELIST /time_control/ auxinput5_interval_y
NAMELIST /time_control/ auxinput5_interval_d
NAMELIST /time_control/ auxinput5_interval_h
NAMELIST /time_control/ auxinput5_interval_m
NAMELIST /time_control/ auxinput5_interval_s
NAMELIST /time_control/ auxinput5_interval
NAMELIST /time_control/ auxinput5_begin_y
NAMELIST /time_control/ auxinput5_begin_d
NAMELIST /time_control/ auxinput5_begin_h
NAMELIST /time_control/ auxinput5_begin_m
NAMELIST /time_control/ auxinput5_begin_s
NAMELIST /time_control/ auxinput5_begin
NAMELIST /time_control/ auxinput5_end_y
NAMELIST /time_control/ auxinput5_end_d
NAMELIST /time_control/ auxinput5_end_h
NAMELIST /time_control/ auxinput5_end_m
NAMELIST /time_control/ auxinput5_end_s
NAMELIST /time_control/ auxinput5_end
NAMELIST /time_control/ io_form_auxinput5
NAMELIST /time_control/ frames_per_auxinput5
NAMELIST /time_control/ auxinput6_inname
NAMELIST /time_control/ auxinput6_outname
NAMELIST /time_control/ auxinput6_interval_y
NAMELIST /time_control/ auxinput6_interval_d
NAMELIST /time_control/ auxinput6_interval_h
NAMELIST /time_control/ auxinput6_interval_m
NAMELIST /time_control/ auxinput6_interval_s
NAMELIST /time_control/ auxinput6_interval
NAMELIST /time_control/ auxinput6_begin_y
NAMELIST /time_control/ auxinput6_begin_d
NAMELIST /time_control/ auxinput6_begin_h
NAMELIST /time_control/ auxinput6_begin_m
NAMELIST /time_control/ auxinput6_begin_s
NAMELIST /time_control/ auxinput6_begin
NAMELIST /time_control/ auxinput6_end_y
NAMELIST /time_control/ auxinput6_end_d
NAMELIST /time_control/ auxinput6_end_h
NAMELIST /time_control/ auxinput6_end_m
NAMELIST /time_control/ auxinput6_end_s
NAMELIST /time_control/ auxinput6_end
NAMELIST /time_control/ io_form_auxinput6
NAMELIST /time_control/ frames_per_auxinput6
NAMELIST /time_control/ auxinput7_inname
NAMELIST /time_control/ auxinput7_outname
NAMELIST /time_control/ auxinput7_interval_y
NAMELIST /time_control/ auxinput7_interval_d
NAMELIST /time_control/ auxinput7_interval_h
NAMELIST /time_control/ auxinput7_interval_m
NAMELIST /time_control/ auxinput7_interval_s
NAMELIST /time_control/ auxinput7_interval
NAMELIST /time_control/ auxinput7_begin_y
NAMELIST /time_control/ auxinput7_begin_d
NAMELIST /time_control/ auxinput7_begin_h
NAMELIST /time_control/ auxinput7_begin_m
NAMELIST /time_control/ auxinput7_begin_s
NAMELIST /time_control/ auxinput7_begin
NAMELIST /time_control/ auxinput7_end_y
NAMELIST /time_control/ auxinput7_end_d
NAMELIST /time_control/ auxinput7_end_h
NAMELIST /time_control/ auxinput7_end_m
NAMELIST /time_control/ auxinput7_end_s
NAMELIST /time_control/ auxinput7_end
NAMELIST /time_control/ io_form_auxinput7
NAMELIST /time_control/ frames_per_auxinput7
NAMELIST /time_control/ auxinput8_inname
NAMELIST /time_control/ auxinput8_outname
NAMELIST /time_control/ auxinput8_interval_y
NAMELIST /time_control/ auxinput8_interval_d
NAMELIST /time_control/ auxinput8_interval_h
NAMELIST /time_control/ auxinput8_interval_m
NAMELIST /time_control/ auxinput8_interval_s
NAMELIST /time_control/ auxinput8_interval
NAMELIST /time_control/ auxinput8_begin_y
NAMELIST /time_control/ auxinput8_begin_d
NAMELIST /time_control/ auxinput8_begin_h
NAMELIST /time_control/ auxinput8_begin_m
NAMELIST /time_control/ auxinput8_begin_s
NAMELIST /time_control/ auxinput8_begin
NAMELIST /time_control/ auxinput8_end_y
NAMELIST /time_control/ auxinput8_end_d
NAMELIST /time_control/ auxinput8_end_h
NAMELIST /time_control/ auxinput8_end_m
NAMELIST /time_control/ auxinput8_end_s
NAMELIST /time_control/ auxinput8_end
NAMELIST /time_control/ io_form_auxinput8
NAMELIST /time_control/ frames_per_auxinput8
NAMELIST /time_control/ auxinput9_inname
NAMELIST /time_control/ auxinput9_outname
NAMELIST /time_control/ auxinput9_interval_y
NAMELIST /time_control/ auxinput9_interval_d
NAMELIST /time_control/ auxinput9_interval_h
NAMELIST /time_control/ auxinput9_interval_m
NAMELIST /time_control/ auxinput9_interval_s
NAMELIST /time_control/ auxinput9_interval
NAMELIST /time_control/ auxinput9_begin_y
NAMELIST /time_control/ auxinput9_begin_d
NAMELIST /time_control/ auxinput9_begin_h
NAMELIST /time_control/ auxinput9_begin_m
NAMELIST /time_control/ auxinput9_begin_s
NAMELIST /time_control/ auxinput9_begin
NAMELIST /time_control/ auxinput9_end_y
NAMELIST /time_control/ auxinput9_end_d
NAMELIST /time_control/ auxinput9_end_h
NAMELIST /time_control/ auxinput9_end_m
NAMELIST /time_control/ auxinput9_end_s
NAMELIST /time_control/ auxinput9_end
NAMELIST /time_control/ io_form_auxinput9
NAMELIST /time_control/ frames_per_auxinput9
NAMELIST /time_control/ auxinput10_inname
NAMELIST /time_control/ auxinput10_outname
NAMELIST /time_control/ auxinput10_interval_y
NAMELIST /time_control/ auxinput10_interval_d
NAMELIST /time_control/ auxinput10_interval_h
NAMELIST /time_control/ auxinput10_interval_m
NAMELIST /time_control/ auxinput10_interval_s
NAMELIST /time_control/ auxinput10_interval
NAMELIST /time_control/ auxinput10_begin_y
NAMELIST /time_control/ auxinput10_begin_d
NAMELIST /time_control/ auxinput10_begin_h
NAMELIST /time_control/ auxinput10_begin_m
NAMELIST /time_control/ auxinput10_begin_s
NAMELIST /time_control/ auxinput10_begin
NAMELIST /time_control/ auxinput10_end_y
NAMELIST /time_control/ auxinput10_end_d
NAMELIST /time_control/ auxinput10_end_h
NAMELIST /time_control/ auxinput10_end_m
NAMELIST /time_control/ auxinput10_end_s
NAMELIST /time_control/ auxinput10_end
NAMELIST /time_control/ io_form_auxinput10
NAMELIST /time_control/ frames_per_auxinput10
NAMELIST /time_control/ auxinput11_inname
NAMELIST /time_control/ auxinput11_outname
NAMELIST /time_control/ auxinput11_interval_y
NAMELIST /time_control/ auxinput11_interval_d
NAMELIST /time_control/ auxinput11_interval_h
NAMELIST /time_control/ auxinput11_interval_m
NAMELIST /time_control/ auxinput11_interval_s
NAMELIST /time_control/ auxinput11_interval
NAMELIST /time_control/ auxinput11_begin_y
NAMELIST /time_control/ auxinput11_begin_d
NAMELIST /time_control/ auxinput11_begin_h
NAMELIST /time_control/ auxinput11_begin_m
NAMELIST /time_control/ auxinput11_begin_s
NAMELIST /time_control/ auxinput11_begin
NAMELIST /time_control/ auxinput11_end_y
NAMELIST /time_control/ auxinput11_end_d
NAMELIST /time_control/ auxinput11_end_h
NAMELIST /time_control/ auxinput11_end_m
NAMELIST /time_control/ auxinput11_end_s
NAMELIST /time_control/ auxinput11_end
NAMELIST /time_control/ io_form_auxinput11
NAMELIST /time_control/ frames_per_auxinput11
NAMELIST /time_control/ auxinput12_inname
NAMELIST /time_control/ auxinput12_outname
NAMELIST /time_control/ auxinput12_interval_y
NAMELIST /time_control/ auxinput12_interval_d
NAMELIST /time_control/ auxinput12_interval_h
NAMELIST /time_control/ auxinput12_interval_m
NAMELIST /time_control/ auxinput12_interval_s
NAMELIST /time_control/ auxinput12_interval
NAMELIST /time_control/ auxinput12_begin_y
NAMELIST /time_control/ auxinput12_begin_d
NAMELIST /time_control/ auxinput12_begin_h
NAMELIST /time_control/ auxinput12_begin_m
NAMELIST /time_control/ auxinput12_begin_s
NAMELIST /time_control/ auxinput12_begin
NAMELIST /time_control/ auxinput12_end_y
NAMELIST /time_control/ auxinput12_end_d
NAMELIST /time_control/ auxinput12_end_h
NAMELIST /time_control/ auxinput12_end_m
NAMELIST /time_control/ auxinput12_end_s
NAMELIST /time_control/ auxinput12_end
NAMELIST /time_control/ io_form_auxinput12
NAMELIST /time_control/ frames_per_auxinput12
NAMELIST /time_control/ auxinput13_inname
NAMELIST /time_control/ auxinput13_outname
NAMELIST /time_control/ auxinput13_interval_y
NAMELIST /time_control/ auxinput13_interval_d
NAMELIST /time_control/ auxinput13_interval_h
NAMELIST /time_control/ auxinput13_interval_m
NAMELIST /time_control/ auxinput13_interval_s
NAMELIST /time_control/ auxinput13_interval
NAMELIST /time_control/ auxinput13_begin_y
NAMELIST /time_control/ auxinput13_begin_d
NAMELIST /time_control/ auxinput13_begin_h
NAMELIST /time_control/ auxinput13_begin_m
NAMELIST /time_control/ auxinput13_begin_s
NAMELIST /time_control/ auxinput13_begin
NAMELIST /time_control/ auxinput13_end_y
NAMELIST /time_control/ auxinput13_end_d
NAMELIST /time_control/ auxinput13_end_h
NAMELIST /time_control/ auxinput13_end_m
NAMELIST /time_control/ auxinput13_end_s
NAMELIST /time_control/ auxinput13_end
NAMELIST /time_control/ io_form_auxinput13
NAMELIST /time_control/ frames_per_auxinput13
NAMELIST /time_control/ auxinput14_inname
NAMELIST /time_control/ auxinput14_outname
NAMELIST /time_control/ auxinput14_interval_y
NAMELIST /time_control/ auxinput14_interval_d
NAMELIST /time_control/ auxinput14_interval_h
NAMELIST /time_control/ auxinput14_interval_m
NAMELIST /time_control/ auxinput14_interval_s
NAMELIST /time_control/ auxinput14_interval
NAMELIST /time_control/ auxinput14_begin_y
NAMELIST /time_control/ auxinput14_begin_d
NAMELIST /time_control/ auxinput14_begin_h
NAMELIST /time_control/ auxinput14_begin_m
NAMELIST /time_control/ auxinput14_begin_s
NAMELIST /time_control/ auxinput14_begin
NAMELIST /time_control/ auxinput14_end_y
NAMELIST /time_control/ auxinput14_end_d
NAMELIST /time_control/ auxinput14_end_h
NAMELIST /time_control/ auxinput14_end_m
NAMELIST /time_control/ auxinput14_end_s
NAMELIST /time_control/ auxinput14_end
NAMELIST /time_control/ io_form_auxinput14
NAMELIST /time_control/ frames_per_auxinput14
NAMELIST /time_control/ auxinput15_inname
NAMELIST /time_control/ auxinput15_outname
NAMELIST /time_control/ auxinput15_interval_y
NAMELIST /time_control/ auxinput15_interval_d
NAMELIST /time_control/ auxinput15_interval_h
NAMELIST /time_control/ auxinput15_interval_m
NAMELIST /time_control/ auxinput15_interval_s
NAMELIST /time_control/ auxinput15_interval
NAMELIST /time_control/ auxinput15_begin_y
NAMELIST /time_control/ auxinput15_begin_d
NAMELIST /time_control/ auxinput15_begin_h
NAMELIST /time_control/ auxinput15_begin_m
NAMELIST /time_control/ auxinput15_begin_s
NAMELIST /time_control/ auxinput15_begin
NAMELIST /time_control/ auxinput15_end_y
NAMELIST /time_control/ auxinput15_end_d
NAMELIST /time_control/ auxinput15_end_h
NAMELIST /time_control/ auxinput15_end_m
NAMELIST /time_control/ auxinput15_end_s
NAMELIST /time_control/ auxinput15_end
NAMELIST /time_control/ io_form_auxinput15
NAMELIST /time_control/ frames_per_auxinput15
NAMELIST /time_control/ auxinput16_inname
NAMELIST /time_control/ auxinput16_outname
NAMELIST /time_control/ auxinput16_interval_y
NAMELIST /time_control/ auxinput16_interval_d
NAMELIST /time_control/ auxinput16_interval_h
NAMELIST /time_control/ auxinput16_interval_m
NAMELIST /time_control/ auxinput16_interval_s
NAMELIST /time_control/ auxinput16_interval
NAMELIST /time_control/ auxinput16_begin_y
NAMELIST /time_control/ auxinput16_begin_d
NAMELIST /time_control/ auxinput16_begin_h
NAMELIST /time_control/ auxinput16_begin_m
NAMELIST /time_control/ auxinput16_begin_s
NAMELIST /time_control/ auxinput16_begin
NAMELIST /time_control/ auxinput16_end_y
NAMELIST /time_control/ auxinput16_end_d
NAMELIST /time_control/ auxinput16_end_h
NAMELIST /time_control/ auxinput16_end_m
NAMELIST /time_control/ auxinput16_end_s
NAMELIST /time_control/ auxinput16_end
NAMELIST /time_control/ io_form_auxinput16
NAMELIST /time_control/ frames_per_auxinput16
NAMELIST /time_control/ auxinput17_inname
NAMELIST /time_control/ auxinput17_outname
NAMELIST /time_control/ auxinput17_interval_y
NAMELIST /time_control/ auxinput17_interval_d
NAMELIST /time_control/ auxinput17_interval_h
NAMELIST /time_control/ auxinput17_interval_m
NAMELIST /time_control/ auxinput17_interval_s
NAMELIST /time_control/ auxinput17_interval
NAMELIST /time_control/ auxinput17_begin_y
NAMELIST /time_control/ auxinput17_begin_d
NAMELIST /time_control/ auxinput17_begin_h
NAMELIST /time_control/ auxinput17_begin_m
NAMELIST /time_control/ auxinput17_begin_s
NAMELIST /time_control/ auxinput17_begin
NAMELIST /time_control/ auxinput17_end_y
NAMELIST /time_control/ auxinput17_end_d
NAMELIST /time_control/ auxinput17_end_h
NAMELIST /time_control/ auxinput17_end_m
NAMELIST /time_control/ auxinput17_end_s
NAMELIST /time_control/ auxinput17_end
NAMELIST /time_control/ io_form_auxinput17
NAMELIST /time_control/ frames_per_auxinput17
NAMELIST /time_control/ auxinput18_inname
NAMELIST /time_control/ auxinput18_outname
NAMELIST /time_control/ auxinput18_interval_y
NAMELIST /time_control/ auxinput18_interval_d
NAMELIST /time_control/ auxinput18_interval_h
NAMELIST /time_control/ auxinput18_interval_m
NAMELIST /time_control/ auxinput18_interval_s
NAMELIST /time_control/ auxinput18_interval
NAMELIST /time_control/ auxinput18_begin_y
NAMELIST /time_control/ auxinput18_begin_d
NAMELIST /time_control/ auxinput18_begin_h
NAMELIST /time_control/ auxinput18_begin_m
NAMELIST /time_control/ auxinput18_begin_s
NAMELIST /time_control/ auxinput18_begin
NAMELIST /time_control/ auxinput18_end_y
NAMELIST /time_control/ auxinput18_end_d
NAMELIST /time_control/ auxinput18_end_h
NAMELIST /time_control/ auxinput18_end_m
NAMELIST /time_control/ auxinput18_end_s
NAMELIST /time_control/ auxinput18_end
NAMELIST /time_control/ io_form_auxinput18
NAMELIST /time_control/ frames_per_auxinput18
NAMELIST /time_control/ auxinput19_inname
NAMELIST /time_control/ auxinput19_outname
NAMELIST /time_control/ auxinput19_interval_y
NAMELIST /time_control/ auxinput19_interval_d
NAMELIST /time_control/ auxinput19_interval_h
NAMELIST /time_control/ auxinput19_interval_m
NAMELIST /time_control/ auxinput19_interval_s
NAMELIST /time_control/ auxinput19_interval
NAMELIST /time_control/ auxinput19_begin_y
NAMELIST /time_control/ auxinput19_begin_d
NAMELIST /time_control/ auxinput19_begin_h
NAMELIST /time_control/ auxinput19_begin_m
NAMELIST /time_control/ auxinput19_begin_s
NAMELIST /time_control/ auxinput19_begin
NAMELIST /time_control/ auxinput19_end_y
NAMELIST /time_control/ auxinput19_end_d
NAMELIST /time_control/ auxinput19_end_h
NAMELIST /time_control/ auxinput19_end_m
NAMELIST /time_control/ auxinput19_end_s
NAMELIST /time_control/ auxinput19_end
NAMELIST /time_control/ io_form_auxinput19
NAMELIST /time_control/ frames_per_auxinput19
NAMELIST /time_control/ auxinput20_inname
NAMELIST /time_control/ auxinput20_outname
NAMELIST /time_control/ auxinput20_interval_y
NAMELIST /time_control/ auxinput20_interval_d
NAMELIST /time_control/ auxinput20_interval_h
NAMELIST /time_control/ auxinput20_interval_m
NAMELIST /time_control/ auxinput20_interval_s
NAMELIST /time_control/ auxinput20_interval
NAMELIST /time_control/ auxinput20_begin_y
NAMELIST /time_control/ auxinput20_begin_d
NAMELIST /time_control/ auxinput20_begin_h
NAMELIST /time_control/ auxinput20_begin_m
NAMELIST /time_control/ auxinput20_begin_s
NAMELIST /time_control/ auxinput20_begin
NAMELIST /time_control/ auxinput20_end_y
NAMELIST /time_control/ auxinput20_end_d
NAMELIST /time_control/ auxinput20_end_h
NAMELIST /time_control/ auxinput20_end_m
NAMELIST /time_control/ auxinput20_end_s
NAMELIST /time_control/ auxinput20_end
NAMELIST /time_control/ io_form_auxinput20
NAMELIST /time_control/ frames_per_auxinput20
NAMELIST /time_control/ auxinput21_inname
NAMELIST /time_control/ auxinput21_outname
NAMELIST /time_control/ auxinput21_interval_y
NAMELIST /time_control/ auxinput21_interval_d
NAMELIST /time_control/ auxinput21_interval_h
NAMELIST /time_control/ auxinput21_interval_m
NAMELIST /time_control/ auxinput21_interval_s
NAMELIST /time_control/ auxinput21_interval
NAMELIST /time_control/ auxinput21_begin_y
NAMELIST /time_control/ auxinput21_begin_d
NAMELIST /time_control/ auxinput21_begin_h
NAMELIST /time_control/ auxinput21_begin_m
NAMELIST /time_control/ auxinput21_begin_s
NAMELIST /time_control/ auxinput21_begin
NAMELIST /time_control/ auxinput21_end_y
NAMELIST /time_control/ auxinput21_end_d
NAMELIST /time_control/ auxinput21_end_h
NAMELIST /time_control/ auxinput21_end_m
NAMELIST /time_control/ auxinput21_end_s
NAMELIST /time_control/ auxinput21_end
NAMELIST /time_control/ io_form_auxinput21
NAMELIST /time_control/ frames_per_auxinput21
NAMELIST /time_control/ auxinput22_inname
NAMELIST /time_control/ auxinput22_outname
NAMELIST /time_control/ auxinput22_interval_y
NAMELIST /time_control/ auxinput22_interval_d
NAMELIST /time_control/ auxinput22_interval_h
NAMELIST /time_control/ auxinput22_interval_m
NAMELIST /time_control/ auxinput22_interval_s
NAMELIST /time_control/ auxinput22_interval
NAMELIST /time_control/ auxinput22_begin_y
NAMELIST /time_control/ auxinput22_begin_d
NAMELIST /time_control/ auxinput22_begin_h
NAMELIST /time_control/ auxinput22_begin_m
NAMELIST /time_control/ auxinput22_begin_s
NAMELIST /time_control/ auxinput22_begin
NAMELIST /time_control/ auxinput22_end_y
NAMELIST /time_control/ auxinput22_end_d
NAMELIST /time_control/ auxinput22_end_h
NAMELIST /time_control/ auxinput22_end_m
NAMELIST /time_control/ auxinput22_end_s
NAMELIST /time_control/ auxinput22_end
NAMELIST /time_control/ io_form_auxinput22
NAMELIST /time_control/ frames_per_auxinput22
NAMELIST /time_control/ auxinput23_inname
NAMELIST /time_control/ auxinput23_outname
NAMELIST /time_control/ auxinput23_interval_y
NAMELIST /time_control/ auxinput23_interval_d
NAMELIST /time_control/ auxinput23_interval_h
NAMELIST /time_control/ auxinput23_interval_m
NAMELIST /time_control/ auxinput23_interval_s
NAMELIST /time_control/ auxinput23_interval
NAMELIST /time_control/ auxinput23_begin_y
NAMELIST /time_control/ auxinput23_begin_d
NAMELIST /time_control/ auxinput23_begin_h
NAMELIST /time_control/ auxinput23_begin_m
NAMELIST /time_control/ auxinput23_begin_s
NAMELIST /time_control/ auxinput23_begin
NAMELIST /time_control/ auxinput23_end_y
NAMELIST /time_control/ auxinput23_end_d
NAMELIST /time_control/ auxinput23_end_h
NAMELIST /time_control/ auxinput23_end_m
NAMELIST /time_control/ auxinput23_end_s
NAMELIST /time_control/ auxinput23_end
NAMELIST /time_control/ io_form_auxinput23
NAMELIST /time_control/ frames_per_auxinput23
NAMELIST /time_control/ auxinput24_inname
NAMELIST /time_control/ auxinput24_outname
NAMELIST /time_control/ auxinput24_interval_y
NAMELIST /time_control/ auxinput24_interval_d
NAMELIST /time_control/ auxinput24_interval_h
NAMELIST /time_control/ auxinput24_interval_m
NAMELIST /time_control/ auxinput24_interval_s
NAMELIST /time_control/ auxinput24_interval
NAMELIST /time_control/ auxinput24_begin_y
NAMELIST /time_control/ auxinput24_begin_d
NAMELIST /time_control/ auxinput24_begin_h
NAMELIST /time_control/ auxinput24_begin_m
NAMELIST /time_control/ auxinput24_begin_s
NAMELIST /time_control/ auxinput24_begin
NAMELIST /time_control/ auxinput24_end_y
NAMELIST /time_control/ auxinput24_end_d
NAMELIST /time_control/ auxinput24_end_h
NAMELIST /time_control/ auxinput24_end_m
NAMELIST /time_control/ auxinput24_end_s
NAMELIST /time_control/ auxinput24_end
NAMELIST /time_control/ io_form_auxinput24
NAMELIST /time_control/ frames_per_auxinput24
NAMELIST /time_control/ history_interval
NAMELIST /time_control/ frames_per_outfile
NAMELIST /time_control/ restart
NAMELIST /time_control/ restart_interval
NAMELIST /time_control/ io_form_input
NAMELIST /time_control/ io_form_history
NAMELIST /time_control/ io_form_restart
NAMELIST /time_control/ io_form_boundary
NAMELIST /time_control/ debug_level
NAMELIST /time_control/ self_test_domain
NAMELIST /time_control/ history_outname
NAMELIST /time_control/ history_inname
NAMELIST /time_control/ use_netcdf_classic
NAMELIST /time_control/ history_interval_d
NAMELIST /time_control/ history_interval_h
NAMELIST /time_control/ history_interval_m
NAMELIST /time_control/ history_interval_s
NAMELIST /time_control/ inputout_interval_d
NAMELIST /time_control/ inputout_interval_h
NAMELIST /time_control/ inputout_interval_m
NAMELIST /time_control/ inputout_interval_s
NAMELIST /time_control/ inputout_interval
NAMELIST /time_control/ restart_interval_d
NAMELIST /time_control/ restart_interval_h
NAMELIST /time_control/ restart_interval_m
NAMELIST /time_control/ restart_interval_s
NAMELIST /time_control/ history_begin_y
NAMELIST /time_control/ history_begin_d
NAMELIST /time_control/ history_begin_h
NAMELIST /time_control/ history_begin_m
NAMELIST /time_control/ history_begin_s
NAMELIST /time_control/ history_begin
NAMELIST /time_control/ inputout_begin_y
NAMELIST /time_control/ inputout_begin_d
NAMELIST /time_control/ inputout_begin_h
NAMELIST /time_control/ inputout_begin_m
NAMELIST /time_control/ inputout_begin_s
NAMELIST /time_control/ restart_begin_y
NAMELIST /time_control/ restart_begin_d
NAMELIST /time_control/ restart_begin_h
NAMELIST /time_control/ restart_begin_m
NAMELIST /time_control/ restart_begin_s
NAMELIST /time_control/ restart_begin
NAMELIST /time_control/ history_end_y
NAMELIST /time_control/ history_end_d
NAMELIST /time_control/ history_end_h
NAMELIST /time_control/ history_end_m
NAMELIST /time_control/ history_end_s
NAMELIST /time_control/ history_end
NAMELIST /time_control/ inputout_end_y
NAMELIST /time_control/ inputout_end_d
NAMELIST /time_control/ inputout_end_h
NAMELIST /time_control/ inputout_end_m
NAMELIST /time_control/ inputout_end_s
NAMELIST /time_control/ reset_simulation_start
NAMELIST /domains/ sr_x
NAMELIST /domains/ sr_y
NAMELIST /time_control/ iofields_filename
NAMELIST /time_control/ ignore_iofields_warning
NAMELIST /time_control/ ncd_nofill
NAMELIST /time_control/ julyr
NAMELIST /time_control/ julday
NAMELIST /time_control/ gmt
NAMELIST /time_control/ high_freq_outname
NAMELIST /time_control/ partial_atcf_outname
NAMELIST /time_control/ input_inname
NAMELIST /time_control/ input_outname
NAMELIST /time_control/ bdy_inname
NAMELIST /time_control/ bdy_outname
NAMELIST /time_control/ rst_inname
NAMELIST /time_control/ rst_outname
NAMELIST /time_control/ anl_outname
NAMELIST /time_control/ write_input
NAMELIST /time_control/ write_restart_at_0h
NAMELIST /time_control/ write_hist_at_0h_rst
NAMELIST /time_control/ adjust_output_times
NAMELIST /time_control/ adjust_input_times
NAMELIST /time_control/ tstart
NAMELIST /time_control/ nocolons
NAMELIST /time_control/ cycling
NAMELIST /time_control/ output_ready_flag
NAMELIST /dfi_control/ dfi_opt
NAMELIST /dfi_control/ dfi_savehydmeteors
NAMELIST /dfi_control/ dfi_nfilter
NAMELIST /dfi_control/ dfi_write_filtered_input
NAMELIST /dfi_control/ dfi_write_dfi_history
NAMELIST /dfi_control/ dfi_cutoff_seconds
NAMELIST /dfi_control/ dfi_time_dim
NAMELIST /dfi_control/ dfi_fwdstop_year
NAMELIST /dfi_control/ dfi_fwdstop_month
NAMELIST /dfi_control/ dfi_fwdstop_day
NAMELIST /dfi_control/ dfi_fwdstop_hour
NAMELIST /dfi_control/ dfi_fwdstop_minute
NAMELIST /dfi_control/ dfi_fwdstop_second
NAMELIST /dfi_control/ dfi_bckstop_year
NAMELIST /dfi_control/ dfi_bckstop_month
NAMELIST /dfi_control/ dfi_bckstop_day
NAMELIST /dfi_control/ dfi_bckstop_hour
NAMELIST /dfi_control/ dfi_bckstop_minute
NAMELIST /dfi_control/ dfi_bckstop_second
NAMELIST /domains/ time_step
NAMELIST /domains/ time_step_fract_num
NAMELIST /domains/ time_step_fract_den
NAMELIST /domains/ time_step_dfi
NAMELIST /domains/ max_dom
NAMELIST /domains/ s_we
NAMELIST /domains/ e_we
NAMELIST /domains/ s_sn
NAMELIST /domains/ e_sn
NAMELIST /domains/ s_vert
NAMELIST /domains/ e_vert
NAMELIST /domains/ num_metgrid_soil_levels
NAMELIST /domains/ dx
NAMELIST /domains/ dy
NAMELIST /domains/ grid_id
NAMELIST /domains/ grid_allowed
NAMELIST /domains/ parent_id
NAMELIST /domains/ i_parent_start
NAMELIST /domains/ j_parent_start
NAMELIST /domains/ parent_grid_ratio
NAMELIST /domains/ parent_time_step_ratio
NAMELIST /domains/ feedback
NAMELIST /domains/ smooth_option
NAMELIST /domains/ ztop
NAMELIST /domains/ moad_grid_ratio
NAMELIST /domains/ moad_time_step_ratio
NAMELIST /domains/ shw
NAMELIST /domains/ tile_sz_x
NAMELIST /domains/ tile_sz_y
NAMELIST /domains/ numtiles
NAMELIST /domains/ numtiles_inc
NAMELIST /domains/ numtiles_x
NAMELIST /domains/ numtiles_y
NAMELIST /domains/ tile_strategy
NAMELIST /domains/ nproc_x
NAMELIST /domains/ nproc_y
NAMELIST /domains/ irand
NAMELIST /domains/ ts_buf_size
NAMELIST /domains/ max_ts_locs
NAMELIST /domains/ ocean_cpl
NAMELIST /domains/ num_moves
NAMELIST /domains/ vortex_interval
NAMELIST /domains/ corral_dist
NAMELIST /domains/ no_ocean
NAMELIST /domains/ move_id
NAMELIST /domains/ move_interval
NAMELIST /domains/ move_cd_x
NAMELIST /domains/ move_cd_y
NAMELIST /domains/ swap_x
NAMELIST /domains/ swap_y
NAMELIST /domains/ cycle_x
NAMELIST /domains/ cycle_y
NAMELIST /domains/ reorder_mesh
NAMELIST /domains/ perturb_input
NAMELIST /domains/ eta_levels
NAMELIST /domains/ ptsgm
NAMELIST /domains/ num_metgrid_levels
NAMELIST /domains/ p_top_requested
NAMELIST /domains/ use_prep_hybrid
NAMELIST /physics/ force_read_thompson
NAMELIST /physics/ write_thompson_tables
NAMELIST /physics/ nssl_cccn
NAMELIST /physics/ nssl_alphah
NAMELIST /physics/ nssl_alphahl
NAMELIST /physics/ nssl_cnoh
NAMELIST /physics/ nssl_cnohl
NAMELIST /physics/ nssl_cnor
NAMELIST /physics/ nssl_cnos
NAMELIST /physics/ nssl_rho_qh
NAMELIST /physics/ nssl_rho_qhl
NAMELIST /physics/ nssl_rho_qs
NAMELIST /physics/ elec_physics
NAMELIST /physics/ nssl_ipelec
NAMELIST /physics/ nssl_isaund
NAMELIST /physics/ mp_physics
NAMELIST /physics/ mommix
NAMELIST /physics/ disheat
NAMELIST /physics/ do_radar_ref
NAMELIST /physics/ ra_lw_physics
NAMELIST /physics/ ra_sw_physics
NAMELIST /physics/ radt
NAMELIST /physics/ sf_sfclay_physics
NAMELIST /physics/ sf_surface_physics
NAMELIST /physics/ bl_pbl_physics
NAMELIST /physics/ bl_mynn_tkebudget
NAMELIST /physics/ ysu_topdown_pblmix
NAMELIST /physics/ shinhong_tke_diag
NAMELIST /physics/ ysu_firstorder
NAMELIST /physics/ ysu_brcr_ub
NAMELIST /physics/ ysu_cgrad_scalar
NAMELIST /physics/ ysu_cgrad_mom
NAMELIST /physics/ windfarm_opt
NAMELIST /physics/ windfarm_ij
NAMELIST /physics/ mfshconv
NAMELIST /physics/ bldt
NAMELIST /physics/ cu_physics
NAMELIST /physics/ shcu_physics
NAMELIST /physics/ cu_diag
NAMELIST /physics/ kf_edrates
NAMELIST /physics/ bl_mynn_tkeadvect
NAMELIST /physics/ bl_mynn_cloudpdf
NAMELIST /physics/ bl_mynn_mixlength
NAMELIST /physics/ bl_mynn_edmf
NAMELIST /physics/ bl_mynn_edmf_mom
NAMELIST /physics/ bl_mynn_edmf_tke
NAMELIST /physics/ bl_mynn_edmf_part
NAMELIST /physics/ bl_mynn_cloudmix
NAMELIST /physics/ bl_mynn_mixqt
NAMELIST /physics/ icloud_bl
NAMELIST /physics/ gfs_alpha
NAMELIST /physics/ cudt
NAMELIST /physics/ gsmdt
NAMELIST /physics/ isfflx
NAMELIST /physics/ ideal_xland
NAMELIST /physics/ ifsnow
NAMELIST /physics/ icloud
NAMELIST /physics/ cldovrlp
NAMELIST /physics/ idcor
NAMELIST /physics/ swrad_scat
NAMELIST /physics/ surface_input_source
NAMELIST /physics/ num_soil_layers
NAMELIST /physics/ num_urban_layers
NAMELIST /physics/ sf_surface_mosaic
NAMELIST /physics/ mosaic_cat
NAMELIST /physics/ num_urban_hi
NAMELIST /physics/ mosaic_lu
NAMELIST /physics/ mosaic_soil
NAMELIST /physics/ maxiens
NAMELIST /physics/ maxens
NAMELIST /physics/ maxens2
NAMELIST /physics/ maxens3
NAMELIST /physics/ ensdim
NAMELIST /physics/ ishallow
NAMELIST /physics/ clos_choice
NAMELIST /physics/ chem_opt
NAMELIST /physics/ num_land_cat
NAMELIST /physics/ num_soil_cat
NAMELIST /physics/ topo_wind
NAMELIST /physics/ mp_zero_out
NAMELIST /physics/ mp_zero_out_thresh
NAMELIST /physics/ seaice_threshold
NAMELIST /physics/ fractional_seaice
NAMELIST /physics/ seaice_albedo_opt
NAMELIST /physics/ seaice_albedo_default
NAMELIST /physics/ seaice_snowdepth_opt
NAMELIST /physics/ seaice_snowdepth_max
NAMELIST /physics/ seaice_snowdepth_min
NAMELIST /physics/ seaice_thickness_opt
NAMELIST /physics/ seaice_thickness_default
NAMELIST /physics/ tice2tsk_if2cold
NAMELIST /physics/ sst_update
NAMELIST /physics/ sf_urban_physics
NAMELIST /physics/ usemonalb
NAMELIST /physics/ rdmaxalb
NAMELIST /physics/ rdlai2d
NAMELIST /physics/ ua_phys
NAMELIST /physics/ opt_thcnd
NAMELIST /physics/ gwd_opt
NAMELIST /physics/ isftcflx
NAMELIST /physics/ iz0tlnd
NAMELIST /physics/ sas_pgcon
NAMELIST /physics/ sas_shal_pgcon
NAMELIST /physics/ sas_shal_conv
NAMELIST /physics/ sas_mass_flux
NAMELIST /physics/ grav_settling
NAMELIST /physics/ var_ric
NAMELIST /physics/ coef_ric_l
NAMELIST /physics/ coef_ric_s
NAMELIST /physics/ random_seed
NAMELIST /physics/ icoef_sf
NAMELIST /physics/ iwavecpl
NAMELIST /physics/ lcurr_sf
NAMELIST /physics/ ens_random_seed
NAMELIST /physics/ ensda_physics_pert
NAMELIST /physics/ pert_sas
NAMELIST /physics/ pert_pbl
NAMELIST /physics/ pert_cd
NAMELIST /physics/ ens_sasamp
NAMELIST /physics/ ens_pblamp
NAMELIST /physics/ ens_cdamp
NAMELIST /physics/ idtad
NAMELIST /physics/ nsoil
NAMELIST /physics/ nphs
NAMELIST /physics/ ncnvc
NAMELIST /physics/ nrand
NAMELIST /physics/ nrads
NAMELIST /physics/ nradl
NAMELIST /physics/ tprec
NAMELIST /physics/ theat
NAMELIST /physics/ tclod
NAMELIST /physics/ trdsw
NAMELIST /physics/ trdlw
NAMELIST /physics/ tsrfc
NAMELIST /physics/ pcpflg
NAMELIST /physics/ sigma
NAMELIST /physics/ sfenth
NAMELIST /physics/ co2tf
NAMELIST /physics/ ra_call_offset
NAMELIST /physics/ cam_abs_freq_s
NAMELIST /physics/ levsiz
NAMELIST /physics/ paerlev
NAMELIST /physics/ cam_abs_dim1
NAMELIST /physics/ cam_abs_dim2
NAMELIST /physics/ no_src_types
NAMELIST /physics/ alevsiz
NAMELIST /physics/ o3input
NAMELIST /physics/ aer_opt
NAMELIST /physics/ cu_rad_feedback
NAMELIST /physics/ aercu_opt
NAMELIST /physics/ aercu_fct
NAMELIST /physics/ h_diff
NAMELIST /physics/ ntrack
NAMELIST /physics/ num_snso_layers
NAMELIST /physics/ num_snow_layers
NAMELIST /physics/ use_aero_icbc
NAMELIST /physics/ use_mp_re
NAMELIST /physics/ ccn_conc
NAMELIST /physics/ hail_opt
NAMELIST /physics/ morr_rimed_ice
NAMELIST /physics/ clean_atm_diag
NAMELIST /physics/ sf_lake_physics
NAMELIST /dynamics/ dyn_opt
NAMELIST /dynamics/ rk_ord
NAMELIST /dynamics/ w_damping
NAMELIST /dynamics/ diff_opt
NAMELIST /dynamics/ km_opt
NAMELIST /dynamics/ damp_opt
NAMELIST /dynamics/ zdamp
NAMELIST /dynamics/ base_pres
NAMELIST /dynamics/ base_temp
NAMELIST /dynamics/ base_lapse
NAMELIST /dynamics/ iso_temp
NAMELIST /dynamics/ dampcoef
NAMELIST /dynamics/ khdif
NAMELIST /dynamics/ kvdif
NAMELIST /dynamics/ c_s
NAMELIST /dynamics/ c_k
NAMELIST /dynamics/ smdiv
NAMELIST /dynamics/ emdiv
NAMELIST /dynamics/ epssm
NAMELIST /dynamics/ nhmove
NAMELIST /dynamics/ non_hydrostatic
NAMELIST /dynamics/ time_step_sound
NAMELIST /dynamics/ h_mom_adv_order
NAMELIST /dynamics/ v_mom_adv_order
NAMELIST /dynamics/ h_sca_adv_order
NAMELIST /dynamics/ v_sca_adv_order
NAMELIST /dynamics/ top_radiation
NAMELIST /dynamics/ tke_upper_bound
NAMELIST /dynamics/ tke_drag_coefficient
NAMELIST /dynamics/ tke_heat_flux
NAMELIST /dynamics/ pert_coriolis
NAMELIST /dynamics/ euler_adv
NAMELIST /dynamics/ idtadt
NAMELIST /dynamics/ idtadc
NAMELIST /dynamics/ codamp
NAMELIST /dynamics/ coac
NAMELIST /dynamics/ slophc
NAMELIST /dynamics/ wp
NAMELIST /dynamics/ dwdt_damping_lev
NAMELIST /dynamics/ terrain_smoothing
NAMELIST /bdy_control/ spec_bdy_width
NAMELIST /bdy_control/ spec_zone
NAMELIST /bdy_control/ relax_zone
NAMELIST /bdy_control/ specified
NAMELIST /bdy_control/ periodic_x
NAMELIST /bdy_control/ symmetric_xs
NAMELIST /bdy_control/ symmetric_xe
NAMELIST /bdy_control/ open_xs
NAMELIST /bdy_control/ open_xe
NAMELIST /bdy_control/ periodic_y
NAMELIST /bdy_control/ symmetric_ys
NAMELIST /bdy_control/ symmetric_ye
NAMELIST /bdy_control/ open_ys
NAMELIST /bdy_control/ open_ye
NAMELIST /bdy_control/ polar
NAMELIST /bdy_control/ nested
NAMELIST /bdy_control/ real_data_init_type
NAMELIST /grib2/ background_proc_id
NAMELIST /grib2/ forecast_proc_id
NAMELIST /grib2/ production_status
NAMELIST /grib2/ compression
NAMELIST /physics/ maxpatch


      CALL MPI_INITIALIZED( mpi_inited, ierr )
      IF ( .NOT. mpi_inited ) THEN
        IF ( coupler_on ) THEN
           CALL cpl_init( mpi_comm_here )
        ELSE
           CALL mpi_init ( ierr )
           mpi_comm_here = MPI_COMM_WORLD
        END IF


        CALL wrf_set_dm_communicator( mpi_comm_here )
        CALL wrf_termio_dup( mpi_comm_here )
      END IF

      CALL wrf_get_dm_communicator( mpi_comm_here )

      CALL MPI_Comm_rank ( mpi_comm_here, mytask_local, ierr ) ;
      CALL MPI_Comm_size ( mpi_comm_here, ntasks_local, ierr ) ;
      mpi_comm_allcompute = mpi_comm_here

      IF ( mytask_local .EQ. 0 ) THEN
        max_dom = 1
        OPEN ( unit=27, file="namelist.input", form="formatted", status="old" )
        READ ( UNIT = 27 , NML = domains , IOSTAT=io_status )
        REWIND(27)
        nio_groups = 1
        nio_tasks_per_group  = 0
        poll_servers = .false.
        READ ( 27 , NML = namelist_quilt, IOSTAT=io_status )
        CLOSE(27)
      END IF
      CALL mpi_bcast( nio_tasks_per_group  , max_domains , MPI_INTEGER , 0 , mpi_comm_here, ierr )
      CALL mpi_bcast( nio_groups , 1 , MPI_INTEGER , 0 , mpi_comm_here, ierr )
      CALL mpi_bcast( max_dom, 1 , MPI_INTEGER , 0 , mpi_comm_here, ierr )
      CALL mpi_bcast( parent_id, max_domains , MPI_INTEGER , 0 , mpi_comm_here, ierr )

      i = MAXVAL(nio_tasks_per_group(1:max_dom))
      IF ( i .GT. 0 .AND. nio_groups .GT. 0 ) THEN
        DO id = 1, max_dom
          IF ( nio_tasks_per_group(id) .LE. 0 ) THEN
             CALL wrf_error_fatal3("<stdin>",4603,&
'If nio_tasks_per_group in namelist.input is non-zero for any domain, every active domain must have a non-zero value in nio_tasks_per_group')
          END IF
        END DO
      END IF

      num_io_tasks = 0
      DO id = 1, max_dom
        num_io_tasks = num_io_tasks + nio_tasks_per_group(id)*nio_groups
      END DO
      CALL nl_set_max_dom(1,max_dom)  

      IF ( mytask_local .EQ. 0 ) THEN
        OPEN ( unit=27, file="namelist.input", form="formatted", status="old" )

        nproc_x = -1
        nproc_y = -1
        READ ( 27 , NML = domains, IOSTAT=io_status )
        CLOSE ( 27 )
        OPEN ( unit=27, file="namelist.input", form="formatted", status="old" )
        tasks_per_split = ntasks_local



        nest_pes_x = 0    
        nest_pes_y = 0
        IF ( nproc_x .EQ. -1 .OR. nproc_y .EQ. -1 ) THEN
          CALL compute_mesh( ntasks_local, n_x, n_y )
        ELSE
          n_x = nproc_x
          n_y = nproc_y
        END IF
        comm_start = 0   
        nest_pes_x(1:max_dom) = n_x
        nest_pes_y(1:max_dom) = n_y
        READ ( 27 , NML = dm_task_split, IOSTAT=io_status )
        CLOSE ( 27 )
      END IF
      CALL mpi_bcast( io_status, 1 , MPI_INTEGER , 0 , mpi_comm_here, ierr )
      IF ( io_status .NE. 0 ) THEN

      END IF
      CALL mpi_bcast( tasks_per_split, 1 , MPI_INTEGER , 0 , mpi_comm_here, ierr )
      CALL mpi_bcast( nproc_x, 1 , MPI_INTEGER , 0 , mpi_comm_here, ierr )
      CALL mpi_bcast( nproc_y, 1 , MPI_INTEGER , 0 , mpi_comm_here, ierr )
      CALL mpi_bcast( comm_start, max_domains , MPI_INTEGER , 0 , mpi_comm_here, ierr )
      CALL mpi_bcast( nest_pes_x, max_domains , MPI_INTEGER , 0 , mpi_comm_here, ierr )
      CALL mpi_bcast( nest_pes_y, max_domains , MPI_INTEGER , 0 , mpi_comm_here, ierr )

      nkids = 1
      which_kid = 0
      DO i = 2, max_dom
        IF ( 1 .le. parent_id(i) .AND. parent_id(i) .LE. max_domains ) THEN
          which_kid(i) = nkids(parent_id(i))
          nkids(parent_id(i)) = nkids(parent_id(i)) + 1
        ELSE
          WRITE(wrf_err_message,*)'invalid parent id for domain ',i
          CALL wrf_error_fatal3("<stdin>",4660,&
TRIM(wrf_err_message))
        END IF
      END DO

      num_compute_tasks = -99
      DO nest_id = 1,max_dom
        IF ( nest_id .EQ. 1 ) THEN
          nest_task_offsets(nest_id) = comm_start(nest_id)
        ELSE
          IF ( comm_start(nest_id) .LT. comm_start(parent_id(nest_id)) ) THEN
            WRITE(wrf_err_message,&
        "('nest domain ',i3,'comm_start (',i3,') lt parent ',i3,' comm_start (',i3,')')") &
                   nest_id,comm_start,parent_id(nest_id),comm_start(parent_id(nest_id))
            CALL wrf_error_fatal3("<stdin>",4674,&
TRIM(wrf_err_message))
          ELSE IF ( comm_start(nest_id) .LT. &
                    comm_start(parent_id(nest_id)) &
                   +nest_pes_x(parent_id(nest_id))*nest_pes_y(parent_id(nest_id))) THEN
            nest_task_offsets(nest_id) = comm_start(nest_id)-comm_start(parent_id(nest_id))
          ELSE
            nest_task_offsets(nest_id) = nest_pes_x(parent_id(nest_id))*nest_pes_y(parent_id(nest_id))
          END IF
        END IF
        IF ((comm_start(nest_id)+nest_pes_x(nest_id)*nest_pes_y(nest_id)) .GT. num_compute_tasks ) THEN
          num_compute_tasks = (comm_start(nest_id)+nest_pes_x(nest_id)*nest_pes_y(nest_id))
        END IF
      END DO

      IF ( .TRUE. ) THEN










        CALL MPI_Comm_rank ( mpi_comm_here, mytask_local, ierr ) ;
        CALL MPI_Comm_rank ( mpi_comm_here, origmytask, ierr ) ;
        CALL mpi_comm_size ( mpi_comm_here, ntasks_local, ierr ) ;
        ALLOCATE( icolor(ntasks_local) )
        ALLOCATE( icolor2(ntasks_local) )
        ALLOCATE( idomain(ntasks_local) )
        k = 0



        comms_i_am_in = MPI_UNDEFINED
        DO i = 1, max_dom
          inthisone = .FALSE.
          icolor = 0
          DO j = comm_start(i), comm_start(i)+nest_pes_x(i)*nest_pes_y(i)-1
            IF ( j+1 .GT. ntasks_local ) THEN
              WRITE(wrf_err_message,*)"check comm_start, nest_pes_x, nest_pes_y settings in namelist for comm ",i
              CALL wrf_error_fatal3("<stdin>",4717,&
wrf_err_message)
            END IF
            icolor(j+1) = 1
          END DO
          IF ( icolor(mytask_local+1) .EQ. 1 ) inthisone = .TRUE.
          CALL MPI_Comm_dup(mpi_comm_here,comdup,ierr)
          CALL MPI_Comm_split(comdup,icolor(mytask_local+1),mytask_local,mpi_comm_local,ierr)
          IF ( inthisone ) THEN
            dims(1) = nest_pes_y(i) 
            dims(2) = nest_pes_x(i)  
            isperiodic(1) = .false.
            isperiodic(2) = .false.
            CALL mpi_cart_create( mpi_comm_local, 2, dims, isperiodic, .false., comms_i_am_in(i), ierr )
          END IF
        END DO


        local_communicator = MPI_UNDEFINED
        DO i = 1, max_dom
          local_communicator_store(i) = comms_i_am_in(i)
          domain_active_this_task(i) = ( local_communicator_store(i) .NE. MPI_UNDEFINED )
          IF ( local_communicator_store(i) .NE. MPI_UNDEFINED ) THEN
             CALL MPI_Comm_size( local_communicator_store(i), ntasks_store(i), ierr )
             CALL MPI_Comm_rank( local_communicator_store(i), mytask_store(i), ierr )
             CALL mpi_cart_coords( local_communicator_store(i), mytask_store(i), 2, coords, ierr )
             IF ( ierr .NE. 0 ) CALL wrf_error_fatal3("<stdin>",4743,&
'MPI_cart_coords fails ')
             mytask_y_store(i) = coords(1)   
             mytask_x_store(i) = coords(2)   
             CALL MPI_Comm_dup( local_communicator_store(i), comdup2, ierr )
             IF ( ierr .NE. 0 ) CALL wrf_error_fatal3("<stdin>",4748,&
'MPI_Comm_dup fails ')

             CALL MPI_Comm_split(comdup2,mytask_y_store(i),mytask_store(i),local_communicator_x_store(i),ierr)
             IF ( ierr .NE. 0 ) CALL wrf_error_fatal3("<stdin>",4752,&
'MPI_Comm_split fails for y ')

             CALL MPI_Comm_split(comdup2,mytask_x_store(i),mytask_store(i),local_communicator_y_store(i),ierr)
             IF ( ierr .NE. 0 ) CALL wrf_error_fatal3("<stdin>",4756,&
'MPI_Comm_split fails for x ')

             CALL MPI_Comm_size( local_communicator_x_store(i), ntasks_x_store(i), ierr )
             CALL MPI_Comm_rank( local_communicator_x_store(i), mytask_x_store(i), ierr )
             CALL MPI_Comm_size( local_communicator_y_store(i), ntasks_y_store(i), ierr )
             CALL MPI_Comm_rank( local_communicator_y_store(i), mytask_y_store(i), ierr )
          END IF
        END DO

        intercomm_active  = .FALSE.
        
        
        
        
        
        


        ntasks_local = num_compute_tasks
        DO nest_id = 2, max_dom
           par_id  = parent_id(nest_id)
           icolor2 = 0
           DO j = 1,ntasks_local 
             IF ( local_communicator_store( par_id ) .NE. MPI_UNDEFINED .OR. local_communicator_store( nest_id ) .NE. MPI_UNDEFINED ) icolor2(j)=1
           END DO
        
           icolor2 = 0
           mytask_is_nest = .FALSE.
           mytask_is_par = .FALSE.
           DO j = 1,ntasks_local

             IF ( comm_start(nest_id) .LE. j-1 .AND. j-1 .LT. comm_start(nest_id) + nest_pes_x(nest_id)*nest_pes_y(nest_id) )  THEN
               icolor2(j)=1
               if ( j-1 .EQ. mytask_local ) mytask_is_nest=.TRUE.
             END IF
             IF ( comm_start(par_id ) .LE. j-1 .AND. j-1 .LT. comm_start(par_id ) + nest_pes_x(par_id )*nest_pes_y(par_id ) )  THEN
               icolor2(j)=1
               if ( j-1 .EQ. mytask_local ) mytask_is_par=.TRUE.
             END IF
           END DO

           i = icolor2(mytask_local+1)
           CALL MPI_Comm_dup(mpi_comm_here,comdup,ierr)
           CALL MPI_Comm_split(comdup,i,origmytask,mpi_comm_me_and_mom,ierr)

           IF ( mytask_is_nest  ) THEN
              intercomm_active(nest_id)  = .TRUE.
              mpi_comm_to_mom(nest_id)   =  mpi_comm_me_and_mom
           END IF
           IF ( mytask_is_par ) THEN
              intercomm_active(par_id)              = .TRUE.
              mpi_comm_to_kid(which_kid(nest_id),par_id)  =  mpi_comm_me_and_mom
           END IF
        END DO
        DEALLOCATE( icolor )
        DEALLOCATE( icolor2 )
        DEALLOCATE( idomain )

      ELSE IF ( ( tasks_per_split .LE. ntasks_local .AND. tasks_per_split .LE. 0 ) ) THEN
        domain_active_this_task(1) = .TRUE.
        IF ( mod( ntasks_local, tasks_per_split ) .NE. 0 ) THEN
          CALL wrf_message( 'WARNING: tasks_per_split does not evenly divide ntasks. Some tasks will be wasted.' )
        END IF

        ALLOCATE( icolor(ntasks_local) )
        j = 0
        DO WHILE ( j .LT. ntasks_local / tasks_per_split )
          DO i = 1, tasks_per_split
            icolor( i + j * tasks_per_split ) = j
          END DO
          j = j + 1
        END DO

        CALL MPI_Comm_dup(mpi_comm_here,comdup,ierr)
        CALL MPI_Comm_split(comdup,icolor(mytask_local+1),mytask_local,mpi_comm_local,ierr)
        CALL wrf_set_dm_communicator( mpi_comm_local )
        CALL store_communicators_for_domain(1)
        DEALLOCATE( icolor )
      ELSE
        domain_active_this_task(1) = .TRUE.
        mpi_comm_local = mpi_comm_here
        CALL wrf_set_dm_communicator( mpi_comm_local )
        CALL store_communicators_for_domain(1)
      END IF

      CALL instate_communicators_for_domain(1)

   END SUBROUTINE split_communicator

   SUBROUTINE init_module_dm
      IMPLICIT NONE
      INTEGER mpi_comm_local, mpi_comm_here, ierr, mytask, nproc
      LOGICAL mpi_inited
      CALL mpi_initialized( mpi_inited, ierr )
      IF ( .NOT. mpi_inited ) THEN
        
        
        
        
        
        CALL mpi_init ( ierr )
        mpi_comm_here = MPI_COMM_WORLD
        CALL wrf_set_dm_communicator ( mpi_comm_here )
      END IF
      CALL wrf_get_dm_communicator( mpi_comm_local )
   END SUBROUTINE init_module_dm


   SUBROUTINE wrf_dm_move_nest ( parent, nest, dx, dy )
      USE module_domain, ONLY : domain
      IMPLICIT NONE
      TYPE (domain), INTENT(INOUT) :: parent, nest
      INTEGER, INTENT(IN)          :: dx,dy
      RETURN
   END SUBROUTINE wrf_dm_move_nest


   SUBROUTINE get_full_obs_vector( nsta, nerrf, niobf,          &
                                   mp_local_uobmask,            &
                                   mp_local_vobmask,            &
                                   mp_local_cobmask, errf )
      





        
    INTEGER, INTENT(IN)   :: nsta                
    INTEGER, INTENT(IN)   :: nerrf               
    INTEGER, INTENT(IN)   :: niobf               
    LOGICAL, INTENT(IN)   :: MP_LOCAL_UOBMASK(NIOBF)
    LOGICAL, INTENT(IN)   :: MP_LOCAL_VOBMASK(NIOBF)
    LOGICAL, INTENT(IN)   :: MP_LOCAL_COBMASK(NIOBF)
    REAL, INTENT(INOUT)   :: errf(nerrf, niobf)

        

    integer i, n, nlocal_dot, nlocal_crs
    REAL UVT_BUFFER(NIOBF)    
    REAL QRK_BUFFER(NIOBF)    
    REAL SFP_BUFFER(NIOBF)    
    REAL PBL_BUFFER(NIOBF)    
    REAL QATOB_BUFFER(NIOBF)  
    INTEGER N_BUFFER(NIOBF)
    REAL FULL_BUFFER(NIOBF)
    INTEGER IFULL_BUFFER(NIOBF)
    INTEGER IDISPLACEMENT(1024)   
    INTEGER ICOUNT(1024)          

    INTEGER :: MPI_COMM_COMP      
    INTEGER :: NPROCS             
    INTEGER :: IERR               


    CALL WRF_GET_DM_COMMUNICATOR(MPI_COMM_COMP)


    CALL MPI_COMM_SIZE( MPI_COMM_COMP, NPROCS, IERR )


   NLOCAL_DOT = 0
   DO N = 1, NSTA
     IF ( MP_LOCAL_UOBMASK(N) ) THEN      
       NLOCAL_DOT = NLOCAL_DOT + 1
       UVT_BUFFER(NLOCAL_DOT) = ERRF(1,N)        
       SFP_BUFFER(NLOCAL_DOT) = ERRF(7,N)        
       QRK_BUFFER(NLOCAL_DOT) = ERRF(9,N)        
       N_BUFFER(NLOCAL_DOT) = N
     END IF
   END DO
   CALL MPI_ALLGATHER(NLOCAL_DOT,1,MPI_INTEGER, &
                      ICOUNT,1,MPI_INTEGER,     &
                      MPI_COMM_COMP,IERR)
   I = 1

   IDISPLACEMENT(1) = 0
   DO I = 2, NPROCS
     IDISPLACEMENT(I) = IDISPLACEMENT(I-1) + ICOUNT(I-1)
   END DO
   CALL MPI_ALLGATHERV( N_BUFFER, NLOCAL_DOT, MPI_INTEGER,    &
                        IFULL_BUFFER, ICOUNT, IDISPLACEMENT,  &
                        MPI_INTEGER, MPI_COMM_COMP, IERR)

   CALL MPI_ALLGATHERV( UVT_BUFFER, NLOCAL_DOT, MPI_REAL,     &
                        FULL_BUFFER, ICOUNT, IDISPLACEMENT,   &
                        MPI_REAL, MPI_COMM_COMP, IERR)
   DO N = 1, NSTA
     ERRF(1,IFULL_BUFFER(N)) = FULL_BUFFER(N)
   END DO

   CALL MPI_ALLGATHERV( SFP_BUFFER, NLOCAL_DOT, MPI_REAL,     &
                        FULL_BUFFER, ICOUNT, IDISPLACEMENT,   &
                        MPI_REAL, MPI_COMM_COMP, IERR)
   DO N = 1, NSTA
     ERRF(7,IFULL_BUFFER(N)) = FULL_BUFFER(N)
   END DO

   CALL MPI_ALLGATHERV( QRK_BUFFER, NLOCAL_DOT, MPI_REAL,     &
                        FULL_BUFFER, ICOUNT, IDISPLACEMENT,   &
                        MPI_REAL, MPI_COMM_COMP, IERR)
   DO N = 1, NSTA
     ERRF(9,IFULL_BUFFER(N)) = FULL_BUFFER(N)
   END DO


   NLOCAL_DOT = 0
   DO N = 1, NSTA
     IF ( MP_LOCAL_VOBMASK(N) ) THEN         
       NLOCAL_DOT = NLOCAL_DOT + 1
       UVT_BUFFER(NLOCAL_DOT) = ERRF(2,N)    
       SFP_BUFFER(NLOCAL_DOT) = ERRF(8,N)    
       N_BUFFER(NLOCAL_DOT) = N
     END IF
   END DO
   CALL MPI_ALLGATHER(NLOCAL_DOT,1,MPI_INTEGER, &
                      ICOUNT,1,MPI_INTEGER,     &
                      MPI_COMM_COMP,IERR)
   I = 1

   IDISPLACEMENT(1) = 0
   DO I = 2, NPROCS
     IDISPLACEMENT(I) = IDISPLACEMENT(I-1) + ICOUNT(I-1)
   END DO
   CALL MPI_ALLGATHERV( N_BUFFER, NLOCAL_DOT, MPI_INTEGER,    &
                        IFULL_BUFFER, ICOUNT, IDISPLACEMENT,  &
                        MPI_INTEGER, MPI_COMM_COMP, IERR)

   CALL MPI_ALLGATHERV( UVT_BUFFER, NLOCAL_DOT, MPI_REAL,     &
                        FULL_BUFFER, ICOUNT, IDISPLACEMENT,   &
                        MPI_REAL, MPI_COMM_COMP, IERR)
   DO N = 1, NSTA
     ERRF(2,IFULL_BUFFER(N)) = FULL_BUFFER(N)
   END DO

   CALL MPI_ALLGATHERV( SFP_BUFFER, NLOCAL_DOT, MPI_REAL,     &
                        FULL_BUFFER, ICOUNT, IDISPLACEMENT,   &
                        MPI_REAL, MPI_COMM_COMP, IERR)
   DO N = 1, NSTA
     ERRF(8,IFULL_BUFFER(N)) = FULL_BUFFER(N)
   END DO


   NLOCAL_CRS = 0
   DO N = 1, NSTA
     IF ( MP_LOCAL_COBMASK(N) ) THEN       
       NLOCAL_CRS = NLOCAL_CRS + 1
       UVT_BUFFER(NLOCAL_CRS) = ERRF(3,N)     
       QRK_BUFFER(NLOCAL_CRS) = ERRF(4,N)     
       PBL_BUFFER(NLOCAL_CRS) = ERRF(5,N)     
       SFP_BUFFER(NLOCAL_CRS) = ERRF(6,N)     
       QATOB_BUFFER(NLOCAL_CRS) = ERRF(10,N)     
       N_BUFFER(NLOCAL_CRS) = N
     END IF
   END DO
   CALL MPI_ALLGATHER(NLOCAL_CRS,1,MPI_INTEGER, &
                      ICOUNT,1,MPI_INTEGER,     &
                      MPI_COMM_COMP,IERR)
   IDISPLACEMENT(1) = 0
   DO I = 2, NPROCS
     IDISPLACEMENT(I) = IDISPLACEMENT(I-1) + ICOUNT(I-1)
   END DO
   CALL MPI_ALLGATHERV( N_BUFFER, NLOCAL_CRS, MPI_INTEGER,    &
                        IFULL_BUFFER, ICOUNT, IDISPLACEMENT,  &
                        MPI_INTEGER, MPI_COMM_COMP, IERR)

   CALL MPI_ALLGATHERV( UVT_BUFFER, NLOCAL_CRS, MPI_REAL,     &
                        FULL_BUFFER, ICOUNT, IDISPLACEMENT,   &
                        MPI_REAL, MPI_COMM_COMP, IERR)

   DO N = 1, NSTA
     ERRF(3,IFULL_BUFFER(N)) = FULL_BUFFER(N)
   END DO

   CALL MPI_ALLGATHERV( QRK_BUFFER, NLOCAL_CRS, MPI_REAL,     &
                        FULL_BUFFER, ICOUNT, IDISPLACEMENT,   &
                        MPI_REAL, MPI_COMM_COMP, IERR)
   DO N = 1, NSTA
     ERRF(4,IFULL_BUFFER(N)) = FULL_BUFFER(N)
   END DO

   CALL MPI_ALLGATHERV( PBL_BUFFER, NLOCAL_CRS, MPI_REAL,     &
                        FULL_BUFFER, ICOUNT, IDISPLACEMENT,   &
                        MPI_REAL, MPI_COMM_COMP, IERR)
   DO N = 1, NSTA
     ERRF(5,IFULL_BUFFER(N)) = FULL_BUFFER(N)
   END DO

   CALL MPI_ALLGATHERV( SFP_BUFFER, NLOCAL_CRS, MPI_REAL,     &
                        FULL_BUFFER, ICOUNT, IDISPLACEMENT,   &
                        MPI_REAL, MPI_COMM_COMP, IERR)
   DO N = 1, NSTA
     ERRF(6,IFULL_BUFFER(N)) = FULL_BUFFER(N)
   END DO


   CALL MPI_ALLGATHERV( QATOB_BUFFER, NLOCAL_CRS, MPI_REAL,     &
                        FULL_BUFFER, ICOUNT, IDISPLACEMENT,   &
                        MPI_REAL, MPI_COMM_COMP, IERR)
   DO N = 1, NSTA
     ERRF(10,IFULL_BUFFER(N)) = FULL_BUFFER(N)
   END DO

   END SUBROUTINE get_full_obs_vector



   SUBROUTINE wrf_dm_maxtile_real ( val , tile)
      IMPLICIT NONE
      REAL val, val_all( ntasks )
      INTEGER tile
      INTEGER ierr






      INTEGER i, comm

      CALL wrf_get_dm_communicator ( comm )
      CALL mpi_allgather ( val, 1, getrealmpitype(), val_all , 1, getrealmpitype(), comm, ierr )
      val = val_all(1)
      tile = 1
      DO i = 2, ntasks
        IF ( val_all(i) .GT. val ) THEN
           tile = i
           val = val_all(i)
        END IF
      END DO
   END SUBROUTINE wrf_dm_maxtile_real


   SUBROUTINE wrf_dm_mintile_real ( val , tile)
      IMPLICIT NONE
      REAL val, val_all( ntasks )
      INTEGER tile
      INTEGER ierr






      INTEGER i, comm

      CALL wrf_get_dm_communicator ( comm )
      CALL mpi_allgather ( val, 1, getrealmpitype(), val_all , 1, getrealmpitype(), comm, ierr )
      val = val_all(1)
      tile = 1
      DO i = 2, ntasks
        IF ( val_all(i) .LT. val ) THEN
           tile = i
           val = val_all(i)
        END IF
      END DO
   END SUBROUTINE wrf_dm_mintile_real


   SUBROUTINE wrf_dm_mintile_double ( val , tile)
      IMPLICIT NONE
      DOUBLE PRECISION val, val_all( ntasks )
      INTEGER tile
      INTEGER ierr






      INTEGER i, comm

      CALL wrf_get_dm_communicator ( comm )
      CALL mpi_allgather ( val, 1, MPI_DOUBLE_PRECISION, val_all , 1, MPI_DOUBLE_PRECISION, comm, ierr )
      val = val_all(1)
      tile = 1
      DO i = 2, ntasks
        IF ( val_all(i) .LT. val ) THEN
           tile = i
           val = val_all(i)
        END IF
      END DO
   END SUBROUTINE wrf_dm_mintile_double


   SUBROUTINE wrf_dm_tile_val_int ( val , tile)
      IMPLICIT NONE
      INTEGER val, val_all( ntasks )
      INTEGER tile
      INTEGER ierr





      INTEGER i, comm

      CALL wrf_get_dm_communicator ( comm )
      CALL mpi_allgather ( val, 1, MPI_INTEGER, val_all , 1, MPI_INTEGER, comm, ierr )
      val = val_all(tile)
   END SUBROUTINE wrf_dm_tile_val_int

   SUBROUTINE wrf_get_hostname  ( str )
      CHARACTER*(*) str
      CHARACTER tmp(512)
      INTEGER i , n, cs
      CALL rsl_lite_get_hostname( tmp, 512, n, cs )
      DO i = 1, n
        str(i:i) = tmp(i)
      END DO
      RETURN
   END SUBROUTINE wrf_get_hostname

   SUBROUTINE wrf_get_hostid  ( hostid )
      INTEGER hostid
      CHARACTER tmp(512)
      INTEGER i, sz, n, cs
      CALL rsl_lite_get_hostname( tmp, 512, n, cs )
      hostid = cs
      RETURN
   END SUBROUTINE wrf_get_hostid

END MODULE module_dm


   SUBROUTINE push_communicators_for_domain( id )
      USE module_dm
      INTEGER, INTENT(IN) :: id   
      IF ( communicator_stack_cursor .GE. max_domains ) CALL wrf_error_fatal3("<stdin>",5185,&
"push_communicators_for_domain would excede stacksize")
      communicator_stack_cursor = communicator_stack_cursor + 1

      id_stack(communicator_stack_cursor) = current_id
      local_communicator_stack( communicator_stack_cursor )    =    local_communicator
      local_communicator_periodic_stack( communicator_stack_cursor )  =    local_communicator_periodic
      local_iocommunicator_stack( communicator_stack_cursor )  =    local_iocommunicator
      local_communicator_x_stack( communicator_stack_cursor )  =    local_communicator_x
      local_communicator_y_stack( communicator_stack_cursor )  =    local_communicator_y
      ntasks_stack( communicator_stack_cursor )        =    ntasks
      ntasks_y_stack( communicator_stack_cursor )      =    ntasks_y
      ntasks_x_stack( communicator_stack_cursor )      =    ntasks_x
      mytask_stack( communicator_stack_cursor )        =    mytask
      mytask_x_stack( communicator_stack_cursor )       =    mytask_x
      mytask_y_stack( communicator_stack_cursor )       =    mytask_y

      CALL instate_communicators_for_domain( id )
   END SUBROUTINE push_communicators_for_domain
   SUBROUTINE pop_communicators_for_domain
      USE module_dm
      IMPLICIT NONE
      IF ( communicator_stack_cursor .LT. 1 ) CALL wrf_error_fatal3("<stdin>",5207,&
"pop_communicators_for_domain on empty stack")
      current_id = id_stack(communicator_stack_cursor)
      local_communicator = local_communicator_stack( communicator_stack_cursor )
      local_communicator_periodic = local_communicator_periodic_stack( communicator_stack_cursor )
      local_iocommunicator = local_iocommunicator_stack( communicator_stack_cursor )
      local_communicator_x = local_communicator_x_stack( communicator_stack_cursor )
      local_communicator_y = local_communicator_y_stack( communicator_stack_cursor )
      ntasks = ntasks_stack( communicator_stack_cursor )
      ntasks_y = ntasks_y_stack( communicator_stack_cursor )
      ntasks_x = ntasks_x_stack( communicator_stack_cursor )
      mytask = mytask_stack( communicator_stack_cursor )
      mytask_x = mytask_x_stack( communicator_stack_cursor )
      mytask_y = mytask_y_stack( communicator_stack_cursor )
      communicator_stack_cursor = communicator_stack_cursor - 1
   END SUBROUTINE pop_communicators_for_domain
   SUBROUTINE instate_communicators_for_domain( id )
      USE module_dm
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: id
      INTEGER ierr
      current_id = id
      local_communicator          = local_communicator_store( id )
      local_communicator_periodic = local_communicator_periodic_store( id )
      local_iocommunicator        = local_iocommunicator_store( id )
      local_communicator_x        = local_communicator_x_store( id )
      local_communicator_y        = local_communicator_y_store( id )
      ntasks         = ntasks_store( id )
      mytask         = mytask_store( id )
      ntasks_x       = ntasks_x_store( id )
      ntasks_y       = ntasks_y_store( id )
      mytask_x       = mytask_x_store( id )
      mytask_y       = mytask_y_store( id )
   END SUBROUTINE instate_communicators_for_domain
   SUBROUTINE store_communicators_for_domain( id )
      USE module_dm
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: id
      local_communicator_store( id )    =    local_communicator
      local_communicator_periodic_store( id )  =    local_communicator_periodic
      local_iocommunicator_store( id )  =    local_iocommunicator
      local_communicator_x_store( id )  =    local_communicator_x
      local_communicator_y_store( id )  =    local_communicator_y
      ntasks_store( id )        =    ntasks
      ntasks_x_store( id )      =    ntasks_x
      ntasks_y_store( id )      =    ntasks_y
      mytask_store( id )        =    mytask
      mytask_x_store( id )      =    mytask_x
      mytask_y_store( id )      =    mytask_y
   END SUBROUTINE store_communicators_for_domain





SUBROUTINE wrf_dm_patch_domain ( id  , domdesc , parent_id , parent_domdesc , &
                          sd1 , ed1 , sp1 , ep1 , sm1 , em1 , &
                          sd2 , ed2 , sp2 , ep2 , sm2 , em2 , &
                          sd3 , ed3 , sp3 , ep3 , sm3 , em3 , &
                                      sp1x , ep1x , sm1x , em1x , &
                                      sp2x , ep2x , sm2x , em2x , &
                                      sp3x , ep3x , sm3x , em3x , &
                                      sp1y , ep1y , sm1y , em1y , &
                                      sp2y , ep2y , sm2y , em2y , &
                                      sp3y , ep3y , sm3y , em3y , &
                          bdx , bdy )
   USE module_domain, ONLY : domain, head_grid, find_grid_by_id
   USE module_dm, ONLY : patch_domain_rsl_lite  
   IMPLICIT NONE

   INTEGER, INTENT(IN)   :: sd1 , ed1 , sd2 , ed2 , sd3 , ed3 , bdx , bdy
   INTEGER, INTENT(OUT)  :: sp1 , ep1 , sp2 , ep2 , sp3 , ep3 , &
                            sm1 , em1 , sm2 , em2 , sm3 , em3
   INTEGER               :: sp1x , ep1x , sp2x , ep2x , sp3x , ep3x , &
                            sm1x , em1x , sm2x , em2x , sm3x , em3x
   INTEGER               :: sp1y , ep1y , sp2y , ep2y , sp3y , ep3y , &
                            sm1y , em1y , sm2y , em2y , sm3y , em3y
   INTEGER, INTENT(INOUT):: id  , domdesc , parent_id , parent_domdesc

   TYPE(domain), POINTER :: parent
   TYPE(domain), POINTER :: grid_ptr

   
   
   
   

   NULLIFY( parent )
   grid_ptr => head_grid
   CALL find_grid_by_id( parent_id , grid_ptr , parent )

   CALL push_communicators_for_domain(id)

   CALL patch_domain_rsl_lite ( id  , parent, parent_id , &
                           sd1 , ed1 , sp1 , ep1 , sm1 , em1 , &
                           sd2 , ed2 , sp2 , ep2 , sm2 , em2 , &
                           sd3 , ed3 , sp3 , ep3 , sm3 , em3 , &
                                      sp1x , ep1x , sm1x , em1x , &
                                      sp2x , ep2x , sm2x , em2x , &
                                      sp3x , ep3x , sm3x , em3x , &
                                      sp1y , ep1y , sm1y , em1y , &
                                      sp2y , ep2y , sm2y , em2y , &
                                      sp3y , ep3y , sm3y , em3y , &
                           bdx , bdy )

   CALL pop_communicators_for_domain

   RETURN
END SUBROUTINE wrf_dm_patch_domain

SUBROUTINE wrf_termio_dup( comm )
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: comm
  INTEGER mytask, ntasks
  INTEGER ierr
  INCLUDE 'mpif.h'
  CALL mpi_comm_size(comm, ntasks, ierr )
  CALL mpi_comm_rank(comm, mytask, ierr )
  write(0,*)'starting wrf task ',mytask,' of ',ntasks
  CALL rsl_error_dup1( mytask )
END SUBROUTINE wrf_termio_dup

SUBROUTINE wrf_get_myproc( myproc )
  USE module_dm , ONLY : mytask
  IMPLICIT NONE
  INTEGER myproc
  myproc = mytask
  RETURN
END SUBROUTINE wrf_get_myproc

SUBROUTINE wrf_get_nproc( nproc )
  USE module_dm , ONLY : ntasks
  IMPLICIT NONE
  INTEGER nproc
  nproc = ntasks
  RETURN
END SUBROUTINE wrf_get_nproc

SUBROUTINE wrf_get_nprocx( nprocx )
  USE module_dm , ONLY : ntasks_x
  IMPLICIT NONE
  INTEGER nprocx
  nprocx = ntasks_x
  RETURN
END SUBROUTINE wrf_get_nprocx

SUBROUTINE wrf_get_nprocy( nprocy )
  USE module_dm , ONLY : ntasks_y
  IMPLICIT NONE
  INTEGER nprocy
  nprocy = ntasks_y
  RETURN
END SUBROUTINE wrf_get_nprocy

SUBROUTINE wrf_dm_bcast_bytes ( buf , size )
   USE module_dm , ONLY : local_communicator
   IMPLICIT NONE
   INCLUDE 'mpif.h'
   INTEGER size
   INTEGER*1 BUF(size)
   CALL BYTE_BCAST ( buf , size, local_communicator )
   RETURN
END SUBROUTINE wrf_dm_bcast_bytes

SUBROUTINE wrf_dm_bcast_string( BUF, N1 )
   IMPLICIT NONE
   INTEGER n1




   CHARACTER*(*) buf
   INTEGER ibuf(256),i,n
   CHARACTER*256 tstr
   n = n1
   
   
   CALL wrf_dm_bcast_integer( n , 1 )
   IF (n .GT. 256) n = 256
   IF (n .GT. 0 ) then
     DO i = 1, n
       ibuf(I) = ichar(buf(I:I))
     END DO
     CALL wrf_dm_bcast_integer( ibuf, n )
     buf = ''
     DO i = 1, n
       buf(i:i) = char(ibuf(i))
     END DO
   END IF
   RETURN
END SUBROUTINE wrf_dm_bcast_string

SUBROUTINE wrf_dm_bcast_string_comm( BUF, N1, COMM )
   IMPLICIT NONE
   INTEGER n1
   INTEGER COMM




   CHARACTER*(*) buf
   INTEGER ibuf(256),i,n
   CHARACTER*256 tstr
   n = n1
   
   
   CALL BYTE_BCAST( n, 4, COMM )
   IF (n .GT. 256) n = 256
   IF (n .GT. 0 ) then
     DO i = 1, n
       ibuf(I) = ichar(buf(I:I))
     END DO
     CALL BYTE_BCAST( ibuf, N*4, COMM )
     buf = ''
     DO i = 1, n
       buf(i:i) = char(ibuf(i))
     END DO
   END IF
   RETURN
END SUBROUTINE wrf_dm_bcast_string_comm

SUBROUTINE wrf_dm_bcast_integer( BUF, N1 )
   IMPLICIT NONE
   INTEGER n1
   INTEGER  buf(*)
   CALL wrf_dm_bcast_bytes ( BUF , N1 * 4 )
   RETURN
END SUBROUTINE wrf_dm_bcast_integer

SUBROUTINE wrf_dm_bcast_double( BUF, N1 )
   IMPLICIT NONE
   INTEGER n1




   REAL  buf(*)
   CALL wrf_dm_bcast_bytes ( BUF , N1 * 8 )
   RETURN
END SUBROUTINE wrf_dm_bcast_double

SUBROUTINE wrf_dm_bcast_real( BUF, N1 )
   IMPLICIT NONE
   INTEGER n1
   REAL  buf(*)
   CALL wrf_dm_bcast_bytes ( BUF , N1 * 4 )
   RETURN
END SUBROUTINE wrf_dm_bcast_real

SUBROUTINE wrf_dm_bcast_logical( BUF, N1 )
   IMPLICIT NONE
   INTEGER n1
   LOGICAL  buf(*)
   CALL wrf_dm_bcast_bytes ( BUF , N1 * 4 )
   RETURN
END SUBROUTINE wrf_dm_bcast_logical

SUBROUTINE write_68( grid, v , s , &
                   ids, ide, jds, jde, kds, kde, &
                   ims, ime, jms, jme, kms, kme, &
                   its, ite, jts, jte, kts, kte )
  USE module_domain, ONLY : domain
  IMPLICIT NONE
  TYPE(domain) , INTENT (INOUT) :: grid
  CHARACTER *(*) s
  INTEGER ids, ide, jds, jde, kds, kde, &
          ims, ime, jms, jme, kms, kme, &
          its, ite, jts, jte, kts, kte
  REAL, DIMENSION( ims:ime , kms:kme, jms:jme ) :: v

  INTEGER i,j,k,ierr

  logical, external :: wrf_dm_on_monitor
  real globbuf( ids:ide, kds:kde, jds:jde )
  character*3 ord, stag

  if ( kds == kde ) then
    ord = 'xy'
    stag = 'xy'
  CALL wrf_patch_to_global_real ( v, globbuf, grid%domdesc, stag, ord, &
                     ids, ide, jds, jde, kds, kde, &
                     ims, ime, jms, jme, kms, kme, &
                     its, ite, jts, jte, kts, kte )
  else

    stag = 'xyz'
    ord = 'xzy'
  CALL wrf_patch_to_global_real ( v, globbuf, grid%domdesc, stag, ord, &
                     ids, ide, kds, kde, jds, jde, &
                     ims, ime, kms, kme, jms, jme, &
                     its, ite, kts, kte, jts, jte )
  endif


  if ( wrf_dm_on_monitor() ) THEN
    WRITE(68,*) ide-ids+1, jde-jds+1 , s
    DO j = jds, jde
    DO i = ids, ide
       WRITE(68,*) globbuf(i,1,j)
    END DO
    END DO
  endif

  RETURN
END

   SUBROUTINE wrf_abort

      USE module_cpl, ONLY : coupler_on, cpl_abort

      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER ierr
      IF ( coupler_on ) THEN
         CALL cpl_abort( 'wrf_abort', 'look for abort message in rsl* files' )
      ELSE
         CALL mpi_abort(MPI_COMM_WORLD,1,ierr)
      END IF
   END SUBROUTINE wrf_abort

   SUBROUTINE wrf_dm_shutdown
      IMPLICIT NONE
      INTEGER ierr
      CALL MPI_FINALIZE( ierr )
      RETURN
   END SUBROUTINE wrf_dm_shutdown

   LOGICAL FUNCTION wrf_dm_on_monitor()
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER tsk, ierr, mpi_comm_local
      CALL wrf_get_dm_communicator( mpi_comm_local )
      IF ( mpi_comm_local .NE. MPI_UNDEFINED ) THEN
        CALL mpi_comm_rank ( mpi_comm_local, tsk , ierr )
        wrf_dm_on_monitor = tsk .EQ. 0
      ELSE
        wrf_dm_on_monitor = .FALSE.
      END IF
      RETURN
   END FUNCTION wrf_dm_on_monitor

   SUBROUTINE rsl_comm_iter_init(shw,ps,pe)
      INTEGER shw, ps, pe
      INTEGER iter, plus_send_start, plus_recv_start, &
                    minus_send_start, minus_recv_start
      COMMON /rcii/ iter, plus_send_start, plus_recv_start, &
                          minus_send_start, minus_recv_start
      iter = 0
      minus_send_start = ps
      minus_recv_start = ps-1
      plus_send_start = pe
      plus_recv_start = pe+1
   END SUBROUTINE rsl_comm_iter_init

   LOGICAL FUNCTION rsl_comm_iter ( id , is_intermediate,                     &
                                    shw ,  xy , ds, de_in, ps, pe, nds,nde, &
                                    sendbeg_m, sendw_m, sendbeg_p, sendw_p,   &
                                    recvbeg_m, recvw_m, recvbeg_p, recvw_p    )
      USE module_dm, ONLY : ntasks_x, ntasks_y, mytask_x, mytask_y, minx, miny, &
                            nest_pes_x, nest_pes_y
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: id,shw,xy,ds,de_in,ps,pe,nds,nde
      LOGICAL, INTENT(IN)  :: is_intermediate  
      INTEGER, INTENT(OUT) :: sendbeg_m, sendw_m, sendbeg_p, sendw_p
      INTEGER, INTENT(OUT) :: recvbeg_m, recvw_m, recvbeg_p, recvw_p
      INTEGER k, kn, ni, nj, de, Px, Py, nt, ntx, nty, me, lb, ub, ierr
      INTEGER dum
      LOGICAL went
      INTEGER iter, plus_send_start, plus_recv_start, &
                    minus_send_start, minus_recv_start
      INTEGER parent_grid_ratio, parent_start
      COMMON /rcii/ iter, plus_send_start, plus_recv_start, &
                          minus_send_start, minus_recv_start



      de = de_in - 1
      ntx = nest_pes_x(id)
      nty = nest_pes_y(id)
      IF ( xy .EQ. 1 ) THEN  
        nt = ntasks_x
        me = mytask_x
        dum = 2 * nty  
        IF ( is_intermediate ) THEN
           CALL nl_get_i_parent_start(id,parent_start)
           CALL nl_get_parent_grid_ratio(id,parent_grid_ratio)
        END IF
      ELSE
        nt = ntasks_y
        me = mytask_y
        dum = 2 * ntx  
        IF ( is_intermediate ) THEN
           CALL nl_get_j_parent_start(id,parent_start)
           CALL nl_get_parent_grid_ratio(id,parent_grid_ratio)
        END IF
      END IF
      iter = iter + 1

      went = .FALSE.
      
      sendw_m = 0
      sendbeg_m = 1
      IF ( me .GT. 0 ) THEN
        lb = minus_send_start
        sendbeg_m = lb-ps+1
        DO k = lb,ps+shw-1
          went = .TRUE.
          IF ( xy .eq. 1 ) THEN
            IF ( is_intermediate ) THEN
              kn =  ( k - parent_start ) * parent_grid_ratio + 1 + 1 ;
              CALL task_for_point (kn,1,nds,nde,1,dum,ntx,nty,Px,Py,minx,miny,ierr) 
              IF ( ierr .NE. 0 ) CALL wrf_error_fatal3("<stdin>",5618,&
'error code returned by task_for_point in module_dm.F (h)')
            ELSE
              CALL task_for_point (k,1,ds,de,1,dum,ntx,nty,Px,Py,minx,miny,ierr) 
              IF ( ierr .NE. 0 ) CALL wrf_error_fatal3("<stdin>",5622,&
'error code returned by task_for_point in module_dm.F (i)')
            END IF
            IF ( Px .NE. me+(iter-1) ) THEN
              exit
            END IF
          ELSE
            IF ( is_intermediate ) THEN
              kn =  ( k - parent_start ) * parent_grid_ratio + 1 + 1 ;
              CALL task_for_point (1,kn,1,dum,nds,nde,ntx,nty,Px,Py,minx,miny,ierr) 
              IF ( ierr .NE. 0 ) CALL wrf_error_fatal3("<stdin>",5632,&
'error code returned by task_for_point in module_dm.F (h)')
            ELSE
              CALL task_for_point (1,k,1,dum,ds,de,ntx,nty,Px,Py,minx,miny,ierr) 
              IF ( ierr .NE. 0 ) CALL wrf_error_fatal3("<stdin>",5636,&
'error code returned by task_for_point in module_dm.F (i)')
            END IF
            IF ( Py .NE. me+(iter-1) ) THEN
              exit
            END IF
          END IF
          minus_send_start = minus_send_start+1
          sendw_m = sendw_m + 1
        END DO
      END IF
      
      recvw_m = 0
      recvbeg_m = 1
      IF ( me .GT. 0 ) THEN
        ub = minus_recv_start
        recvbeg_m = ps - ub
        DO k = minus_recv_start,ps-shw,-1
          went = .TRUE.
          IF ( xy .eq. 1 ) THEN
          IF ( is_intermediate ) THEN
            kn =  ( k - parent_start ) * parent_grid_ratio + 1 + 1 ;
              CALL task_for_point (kn,1,nds,nde,1,dum,ntx,nty,Px,Py,minx,miny,ierr) 
              IF ( ierr .NE. 0 ) CALL wrf_error_fatal3("<stdin>",5659,&
'error code returned by task_for_point in module_dm.F (j)')
          ELSE
              CALL task_for_point (k,1,ds,de,1,dum,ntx,nty,Px,Py,minx,miny,ierr) 
              IF ( ierr .NE. 0 ) CALL wrf_error_fatal3("<stdin>",5663,&
'error code returned by task_for_point in module_dm.F (k)')
          END IF
          IF ( Px .NE. me-iter ) THEN
            exit
          END IF
          ELSE
            IF ( is_intermediate ) THEN
              kn =  ( k - parent_start ) * parent_grid_ratio + 1 + 1 ;
              CALL task_for_point (1,kn,1,dum,nds,nde,ntx,nty,Px,Py,minx,miny,ierr) 
              IF ( ierr .NE. 0 ) CALL wrf_error_fatal3("<stdin>",5673,&
'error code returned by task_for_point in module_dm.F (j)')
            ELSE
              CALL task_for_point (1,k,1,dum,ds,de,ntx,nty,Px,Py,minx,miny,ierr) 
              IF ( ierr .NE. 0 ) CALL wrf_error_fatal3("<stdin>",5677,&
'error code returned by task_for_point in module_dm.F (k)')
            END IF
            IF ( Py .NE. me-iter ) THEN
              exit
            END IF
          END IF
          minus_recv_start = minus_recv_start-1
          recvw_m = recvw_m + 1
        END DO
      END IF

      
      sendw_p = 0
      sendbeg_p = 1
      IF ( ( xy .eq. 1 .and. me .LT. ntx-1 ) .OR. ( xy .eq. 0 .and. me .LT. nty-1 ) ) THEN
        ub = plus_send_start
        sendbeg_p = pe - ub + 1
        DO k = ub,pe-shw+1,-1
          went = .TRUE.
          IF ( xy .eq. 1 ) THEN
          IF ( is_intermediate ) THEN
            kn =  ( k - parent_start ) * parent_grid_ratio + 1 + 1 ;
              CALL task_for_point (kn,1,nds,nde,1,dum,ntx,nty,Px,Py,minx,miny,ierr) 
              IF ( ierr .NE. 0 ) CALL wrf_error_fatal3("<stdin>",5701,&
'error code returned by task_for_point in module_dm.F (l)')
          ELSE
              CALL task_for_point (k,1,ds,de,1,dum,ntx,nty,Px,Py,minx,miny,ierr) 
              IF ( ierr .NE. 0 ) CALL wrf_error_fatal3("<stdin>",5705,&
'error code returned by task_for_point in module_dm.F (m)')
          END IF
          IF ( Px .NE. me-(iter-1) ) THEN
            exit
          END IF
          ELSE
            IF ( is_intermediate ) THEN
              kn =  ( k - parent_start ) * parent_grid_ratio + 1 + 1 ;
              CALL task_for_point (1,kn,1,dum,nds,nde,ntx,nty,Px,Py,minx,miny,ierr) 
              IF ( ierr .NE. 0 ) CALL wrf_error_fatal3("<stdin>",5715,&
'error code returned by task_for_point in module_dm.F (l)')
            ELSE
              CALL task_for_point (1,k,1,dum,ds,de,ntx,nty,Px,Py,minx,miny,ierr) 
              IF ( ierr .NE. 0 ) CALL wrf_error_fatal3("<stdin>",5719,&
'error code returned by task_for_point in module_dm.F (m)')
            END IF
            IF ( Py .NE. me-(iter-1) ) THEN
              exit
            END IF
          END IF
          plus_send_start = plus_send_start - 1
          sendw_p = sendw_p + 1
        END DO
      END IF
      
      recvw_p = 0
      recvbeg_p = 1
      IF ( ( xy .eq. 1 .and. me .LT. ntx-1 ) .OR. ( xy .eq. 0 .and. me .LT. nty-1 ) ) THEN
        lb = plus_recv_start
        recvbeg_p = lb - pe
        DO k = lb,pe+shw
          went = .TRUE.
          IF ( xy .eq. 1 ) THEN
          IF ( is_intermediate ) THEN
            kn =  ( k - parent_start ) * parent_grid_ratio + 1 + 1 ;
              CALL task_for_point (kn,1,nds,nde,1,dum,ntx,nty,Px,Py,minx,miny,ierr) 
              IF ( ierr .NE. 0 ) CALL wrf_error_fatal3("<stdin>",5742,&
'error code returned by task_for_point in module_dm.F (n)')
          ELSE

              CALL task_for_point (k,1,ds,de,1,dum,ntx,nty,Px,Py,minx,miny,ierr) 

              IF ( ierr .NE. 0 ) CALL wrf_error_fatal3("<stdin>",5748,&
'error code returned by task_for_point in module_dm.F (o)')
          END IF
          IF ( Px .NE. me+iter ) THEN
            exit
          END IF
          ELSE
            IF ( is_intermediate ) THEN
              kn =  ( k - parent_start ) * parent_grid_ratio + 1 + 1 ;
              CALL task_for_point (1,kn,1,dum,nds,nde,ntx,nty,Px,Py,minx,miny,ierr) 
              IF ( ierr .NE. 0 ) CALL wrf_error_fatal3("<stdin>",5758,&
'error code returned by task_for_point in module_dm.F (n)')
            ELSE
              CALL task_for_point (1,k,1,dum,ds,de,ntx,nty,Px,Py,minx,miny,ierr) 
              IF ( ierr .NE. 0 ) CALL wrf_error_fatal3("<stdin>",5762,&
'error code returned by task_for_point in module_dm.F (o)')
            END IF
            IF ( Py .NE. me+iter ) THEN
              exit
            END IF
          END IF
          plus_recv_start = plus_recv_start + 1
          recvw_p = recvw_p + 1
        END DO
      END IF
      



      
      
      rsl_comm_iter = went
   END FUNCTION rsl_comm_iter

   INTEGER FUNCTION wrf_dm_monitor_rank()
      IMPLICIT NONE
      wrf_dm_monitor_rank = 0
      RETURN
   END FUNCTION wrf_dm_monitor_rank


   SUBROUTINE wrf_get_dm_communicator_for_id ( id, communicator )
      USE module_dm , ONLY : local_communicator_store, mpi_comm_allcompute
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: id
      INTEGER , INTENT(OUT) :: communicator
      IF ( id .le. 0 ) THEN
        communicator = mpi_comm_allcompute
      ELSE
        communicator = local_communicator_store(id)
      END IF
      RETURN
   END SUBROUTINE wrf_get_dm_communicator_for_id

   SUBROUTINE wrf_get_dm_communicator ( communicator )
      USE module_dm , ONLY : local_communicator
      IMPLICIT NONE
      INTEGER , INTENT(OUT) :: communicator
      communicator = local_communicator
      RETURN
   END SUBROUTINE wrf_get_dm_communicator

   SUBROUTINE wrf_get_dm_communicator_x ( communicator )
      USE module_dm , ONLY : local_communicator_x
      IMPLICIT NONE
      INTEGER , INTENT(OUT) :: communicator
      communicator = local_communicator_x
      RETURN
   END SUBROUTINE wrf_get_dm_communicator_x

   SUBROUTINE wrf_get_dm_communicator_y ( communicator )
      USE module_dm , ONLY : local_communicator_y
      IMPLICIT NONE
      INTEGER , INTENT(OUT) :: communicator
      communicator = local_communicator_y
      RETURN
   END SUBROUTINE wrf_get_dm_communicator_y

   SUBROUTINE wrf_get_dm_iocommunicator ( iocommunicator )
      USE module_dm , ONLY : local_iocommunicator
      IMPLICIT NONE
      INTEGER , INTENT(OUT) :: iocommunicator
      iocommunicator = local_iocommunicator
      RETURN
   END SUBROUTINE wrf_get_dm_iocommunicator

   SUBROUTINE wrf_set_dm_communicator ( communicator )
      USE module_dm , ONLY : local_communicator
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: communicator
      local_communicator = communicator
      RETURN
   END SUBROUTINE wrf_set_dm_communicator

   SUBROUTINE wrf_set_dm_iocommunicator ( iocommunicator )
      USE module_dm , ONLY : local_iocommunicator
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: iocommunicator
      local_iocommunicator = iocommunicator
      RETURN
   END SUBROUTINE wrf_set_dm_iocommunicator

   SUBROUTINE wrf_get_dm_ntasks_x ( retval )
      USE module_dm , ONLY : ntasks_x
      IMPLICIT NONE
      INTEGER , INTENT(OUT) :: retval
      retval = ntasks_x
      RETURN
   END SUBROUTINE wrf_get_dm_ntasks_x

   SUBROUTINE wrf_get_dm_ntasks_y ( retval )
      USE module_dm , ONLY : ntasks_y
      IMPLICIT NONE
      INTEGER , INTENT(OUT) :: retval
      retval = ntasks_y
      RETURN
   END SUBROUTINE wrf_get_dm_ntasks_y


   SUBROUTINE wrf_set_dm_quilt_comm ( communicator )
      USE module_dm , ONLY : local_quilt_comm
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: communicator
      local_quilt_comm = communicator
      RETURN
   END SUBROUTINE wrf_set_dm_quilt_comm

   SUBROUTINE wrf_get_dm_quilt_comm ( communicator )
      USE module_dm , ONLY : local_quilt_comm
      IMPLICIT NONE
      INTEGER , INTENT(OUT) :: communicator
      communicator = local_quilt_comm
      RETURN
   END SUBROUTINE wrf_get_dm_quilt_comm




   SUBROUTINE wrf_patch_to_global_real (buf,globbuf,domdesc,stagger,ordering,&
                                       DS1,DE1,DS2,DE2,DS3,DE3,&
                                       MS1,ME1,MS2,ME2,MS3,ME3,&
                                       PS1,PE1,PS2,PE2,PS3,PE3 )
       IMPLICIT NONE
       INTEGER                         DS1,DE1,DS2,DE2,DS3,DE3,&
                                       MS1,ME1,MS2,ME2,MS3,ME3,&
                                       PS1,PE1,PS2,PE2,PS3,PE3
       CHARACTER *(*) stagger,ordering
       INTEGER fid,domdesc
       REAL globbuf(*)
       REAL buf(*)

       CALL wrf_patch_to_global_generic (buf,globbuf,domdesc,stagger,ordering,4,&
                                         DS1,DE1,DS2,DE2,DS3,DE3,&
                                         MS1,ME1,MS2,ME2,MS3,ME3,&
                                         PS1,PE1,PS2,PE2,PS3,PE3 )

       RETURN
   END SUBROUTINE wrf_patch_to_global_real

   SUBROUTINE wrf_patch_to_global_double (buf,globbuf,domdesc,stagger,ordering,&
                                       DS1,DE1,DS2,DE2,DS3,DE3,&
                                       MS1,ME1,MS2,ME2,MS3,ME3,&
                                       PS1,PE1,PS2,PE2,PS3,PE3 )
       IMPLICIT NONE
       INTEGER                         DS1,DE1,DS2,DE2,DS3,DE3,&
                                       MS1,ME1,MS2,ME2,MS3,ME3,&
                                       PS1,PE1,PS2,PE2,PS3,PE3
       CHARACTER *(*) stagger,ordering
       INTEGER fid,domdesc




       REAL globbuf(*)
       REAL buf(*)

       CALL wrf_patch_to_global_generic (buf,globbuf,domdesc,stagger,ordering,8,&
                                         DS1,DE1,DS2,DE2,DS3,DE3,&
                                         MS1,ME1,MS2,ME2,MS3,ME3,&
                                         PS1,PE1,PS2,PE2,PS3,PE3 )

       RETURN
   END SUBROUTINE wrf_patch_to_global_double


   SUBROUTINE wrf_patch_to_global_integer (buf,globbuf,domdesc,stagger,ordering,&
                                       DS1,DE1,DS2,DE2,DS3,DE3,&
                                       MS1,ME1,MS2,ME2,MS3,ME3,&
                                       PS1,PE1,PS2,PE2,PS3,PE3 )
       IMPLICIT NONE
       INTEGER                         DS1,DE1,DS2,DE2,DS3,DE3,&
                                       MS1,ME1,MS2,ME2,MS3,ME3,&
                                       PS1,PE1,PS2,PE2,PS3,PE3
       CHARACTER *(*) stagger,ordering
       INTEGER fid,domdesc
       INTEGER globbuf(*)
       INTEGER buf(*)

       CALL wrf_patch_to_global_generic (buf,globbuf,domdesc,stagger,ordering,4,&
                                         DS1,DE1,DS2,DE2,DS3,DE3,&
                                         MS1,ME1,MS2,ME2,MS3,ME3,&
                                         PS1,PE1,PS2,PE2,PS3,PE3 )

       RETURN
   END SUBROUTINE wrf_patch_to_global_integer


   SUBROUTINE wrf_patch_to_global_logical (buf,globbuf,domdesc,stagger,ordering,&
                                       DS1,DE1,DS2,DE2,DS3,DE3,&
                                       MS1,ME1,MS2,ME2,MS3,ME3,&
                                       PS1,PE1,PS2,PE2,PS3,PE3 )
       IMPLICIT NONE
       INTEGER                         DS1,DE1,DS2,DE2,DS3,DE3,&
                                       MS1,ME1,MS2,ME2,MS3,ME3,&
                                       PS1,PE1,PS2,PE2,PS3,PE3
       CHARACTER *(*) stagger,ordering
       INTEGER fid,domdesc
       LOGICAL globbuf(*)
       LOGICAL buf(*)

       CALL wrf_patch_to_global_generic (buf,globbuf,domdesc,stagger,ordering,4,&
                                         DS1,DE1,DS2,DE2,DS3,DE3,&
                                         MS1,ME1,MS2,ME2,MS3,ME3,&
                                         PS1,PE1,PS2,PE2,PS3,PE3 )

       RETURN
   END SUBROUTINE wrf_patch_to_global_logical


   SUBROUTINE wrf_patch_to_global_generic (buf,globbuf,domdesc,stagger,ordering,typesize,&
                                       DS1a,DE1a,DS2a,DE2a,DS3a,DE3a,&
                                       MS1a,ME1a,MS2a,ME2a,MS3a,ME3a,&
                                       PS1a,PE1a,PS2a,PE2a,PS3a,PE3a )
       USE module_driver_constants
       USE module_timing
       USE module_wrf_error, ONLY : wrf_at_debug_level
       USE module_dm, ONLY : local_communicator, ntasks

       IMPLICIT NONE
       INTEGER                         DS1a,DE1a,DS2a,DE2a,DS3a,DE3a,&
                                       MS1a,ME1a,MS2a,ME2a,MS3a,ME3a,&
                                       PS1a,PE1a,PS2a,PE2a,PS3a,PE3A
       CHARACTER *(*) stagger,ordering
       INTEGER domdesc,typesize,ierr
       REAL globbuf(*)
       REAL buf(*)
       INTEGER                         DS1,DE1,DS2,DE2,DS3,DE3,&
                                       MS1,ME1,MS2,ME2,MS3,ME3,&
                                       PS1,PE1,PS2,PE2,PS3,PE3
       INTEGER                         ids,ide,jds,jde,kds,kde,&
                                       ims,ime,jms,jme,kms,kme,&
                                       ips,ipe,jps,jpe,kps,kpe
       LOGICAL, EXTERNAL :: wrf_dm_on_monitor, has_char

       INTEGER i, j, k,  ndim
       INTEGER  Patch(3,2), Gpatch(3,2,ntasks)
    
       REAL, ALLOCATABLE :: tmpbuf( : )
       REAL locbuf( (PE1a-PS1a+1)*(PE2a-PS2a+1)*(PE3a-PS3a+1)/4*typesize+32 )

       DS1 = DS1a ; DE1 = DE1a ; DS2=DS2a ; DE2 = DE2a ; DS3 = DS3a ; DE3 = DE3a
       MS1 = MS1a ; ME1 = ME1a ; MS2=MS2a ; ME2 = ME2a ; MS3 = MS3a ; ME3 = ME3a
       PS1 = PS1a ; PE1 = PE1a ; PS2=PS2a ; PE2 = PE2a ; PS3 = PS3a ; PE3 = PE3a

       SELECT CASE ( TRIM(ordering) )
         CASE ( 'xy', 'yx' )
           ndim = 2
         CASE DEFAULT
           ndim = 3   
       END SELECT

       SELECT CASE ( TRIM(ordering) )
         CASE ( 'xyz','xy' )
            
            
            
           IF ( .NOT. has_char( stagger, 'x' ) ) DE1 = DE1+1
           IF ( .NOT. has_char( stagger, 'y' ) ) DE2 = DE2+1
           IF ( ndim .EQ. 3 .AND. .NOT. has_char( stagger, 'z' ) ) DE3 = DE3+1
         CASE ( 'yxz','yx' )
           IF ( .NOT. has_char( stagger, 'x' ) ) DE2 = DE2+1
           IF ( .NOT. has_char( stagger, 'y' ) ) DE1 = DE1+1
           IF ( ndim .EQ. 3 .AND. .NOT. has_char( stagger, 'z' ) ) DE3 = DE3+1
         CASE ( 'zxy' )
           IF ( .NOT. has_char( stagger, 'x' ) ) DE2 = DE2+1
           IF ( .NOT. has_char( stagger, 'y' ) ) DE3 = DE3+1
           IF ( ndim .EQ. 3 .AND. .NOT. has_char( stagger, 'z' ) ) DE1 = DE1+1
         CASE ( 'xzy' )
           IF ( .NOT. has_char( stagger, 'x' ) ) DE1 = DE1+1
           IF ( .NOT. has_char( stagger, 'y' ) ) DE3 = DE3+1
           IF ( ndim .EQ. 3 .AND. .NOT. has_char( stagger, 'z' ) ) DE2 = DE2+1
         CASE DEFAULT
       END SELECT

     
       IF ( wrf_dm_on_monitor() ) THEN
         ALLOCATE ( tmpbuf ( (DE1-DS1+1)*(DE2-DS2+1)*(DE3-DS3+1)/4*typesize+32 ), STAT=ierr )
       ELSE
         ALLOCATE ( tmpbuf ( 1 ), STAT=ierr )
       END IF
       IF ( ierr .ne. 0 ) CALL wrf_error_fatal3("<stdin>",6048,&
'allocating tmpbuf in wrf_patch_to_global_generic')

       Patch(1,1) = ps1 ; Patch(1,2) = pe1    
       Patch(2,1) = ps2 ; Patch(2,2) = pe2
       Patch(3,1) = ps3 ; Patch(3,2) = pe3

       IF      ( typesize .EQ. 4 ) THEN
         CALL just_patch_r ( buf , locbuf , size(locbuf)*4/typesize, &
                                   PS1, PE1, PS2, PE2, PS3, PE3 , &
                                   MS1, ME1, MS2, ME2, MS3, ME3   )
       ELSE IF ( typesize .EQ. 4 ) THEN
         CALL just_patch_i ( buf , locbuf , size(locbuf)*4/typesize, &
                                   PS1, PE1, PS2, PE2, PS3, PE3 , &
                                   MS1, ME1, MS2, ME2, MS3, ME3   )
       ELSE IF ( typesize .EQ. 8 ) THEN
         CALL just_patch_d ( buf , locbuf , size(locbuf)*4/typesize, &
                                   PS1, PE1, PS2, PE2, PS3, PE3 , &
                                   MS1, ME1, MS2, ME2, MS3, ME3   )
       ELSE IF ( typesize .EQ. 4 ) THEN
         CALL just_patch_l ( buf , locbuf , size(locbuf)*4/typesize, &
                                   PS1, PE1, PS2, PE2, PS3, PE3 , &
                                   MS1, ME1, MS2, ME2, MS3, ME3   )
       END IF


       CALL collect_on_comm0 (  local_communicator , 4 ,  &
                                Patch , 6 ,                       &
                                GPatch , 6*ntasks                 )

       CALL collect_on_comm0 (  local_communicator , typesize ,  &
                                locbuf , (pe1-ps1+1)*(pe2-ps2+1)*(pe3-ps3+1),   &
                                tmpbuf  , (de1-ds1+1)*(de2-ds2+1)*(de3-ds3+1) )

       ndim = len(TRIM(ordering))

       IF ( wrf_at_debug_level(500) ) THEN
         CALL start_timing
       END IF

       IF ( ndim .GE. 2 .AND. wrf_dm_on_monitor() ) THEN

         IF      ( typesize .EQ. 4 ) THEN
           CALL patch_2_outbuf_r ( tmpbuf  , globbuf ,             &
                                   DS1, DE1, DS2, DE2, DS3, DE3 , &
                                   GPATCH                         )
         ELSE IF ( typesize .EQ. 4 ) THEN
           CALL patch_2_outbuf_i ( tmpbuf  , globbuf ,             &
                                   DS1, DE1, DS2, DE2, DS3, DE3 , &
                                   GPATCH                         )
         ELSE IF ( typesize .EQ. 8 ) THEN
           CALL patch_2_outbuf_d ( tmpbuf  , globbuf ,             &
                                   DS1, DE1, DS2, DE2, DS3, DE3 , &
                                   GPATCH                         )
         ELSE IF ( typesize .EQ. 4 ) THEN
           CALL patch_2_outbuf_l ( tmpbuf  , globbuf ,             &
                                   DS1, DE1, DS2, DE2, DS3, DE3 , &
                                   GPATCH                         )
         END IF

       END IF

       IF ( wrf_at_debug_level(500) ) THEN
         CALL end_timing('wrf_patch_to_global_generic')
       END IF
       DEALLOCATE( tmpbuf )
       RETURN
    END SUBROUTINE wrf_patch_to_global_generic

  SUBROUTINE just_patch_i ( inbuf , outbuf, noutbuf,     &
                               PS1,PE1,PS2,PE2,PS3,PE3,  &
                               MS1,ME1,MS2,ME2,MS3,ME3   )
    IMPLICIT NONE
    INTEGER                         , INTENT(IN)  :: noutbuf
    INTEGER    , DIMENSION(noutbuf) , INTENT(OUT) :: outbuf
    INTEGER   MS1,ME1,MS2,ME2,MS3,ME3
    INTEGER   PS1,PE1,PS2,PE2,PS3,PE3
    INTEGER    , DIMENSION( MS1:ME1,MS2:ME2,MS3:ME3 ) , INTENT(IN) :: inbuf

    INTEGER               :: i,j,k,n   ,  icurs
    icurs = 1
      DO k = PS3, PE3
        DO j = PS2, PE2
          DO i = PS1, PE1
            outbuf( icurs )  = inbuf( i, j, k )
            icurs = icurs + 1
          END DO
        END DO
      END DO
    RETURN
  END SUBROUTINE just_patch_i

  SUBROUTINE just_patch_r ( inbuf , outbuf, noutbuf,     &
                               PS1,PE1,PS2,PE2,PS3,PE3,  &
                               MS1,ME1,MS2,ME2,MS3,ME3   )
    IMPLICIT NONE
    INTEGER                      , INTENT(IN)  :: noutbuf
    REAL    , DIMENSION(noutbuf) , INTENT(OUT) :: outbuf
    INTEGER   MS1,ME1,MS2,ME2,MS3,ME3
    INTEGER   PS1,PE1,PS2,PE2,PS3,PE3
    REAL    , DIMENSION( MS1:ME1,MS2:ME2,MS3:ME3 ) , INTENT(in) :: inbuf

    INTEGER               :: i,j,k   ,  icurs
    icurs = 1
      DO k = PS3, PE3
        DO j = PS2, PE2
          DO i = PS1, PE1
            outbuf( icurs )  = inbuf( i, j, k )
            icurs = icurs + 1
          END DO
        END DO
      END DO
    RETURN
  END SUBROUTINE just_patch_r

  SUBROUTINE just_patch_d ( inbuf , outbuf, noutbuf,     &
                               PS1,PE1,PS2,PE2,PS3,PE3,  &
                               MS1,ME1,MS2,ME2,MS3,ME3   )
    IMPLICIT NONE
    INTEGER                                  , INTENT(IN)  :: noutbuf
    DOUBLE PRECISION    , DIMENSION(noutbuf) , INTENT(OUT) :: outbuf
    INTEGER   MS1,ME1,MS2,ME2,MS3,ME3
    INTEGER   PS1,PE1,PS2,PE2,PS3,PE3
    DOUBLE PRECISION    , DIMENSION( MS1:ME1,MS2:ME2,MS3:ME3 ) , INTENT(in) :: inbuf

    INTEGER               :: i,j,k,n   ,  icurs
    icurs = 1
      DO k = PS3, PE3
        DO j = PS2, PE2
          DO i = PS1, PE1
            outbuf( icurs )  = inbuf( i, j, k )
            icurs = icurs + 1
          END DO
        END DO
      END DO
    RETURN
  END SUBROUTINE just_patch_d

  SUBROUTINE just_patch_l ( inbuf , outbuf, noutbuf,     &
                               PS1,PE1,PS2,PE2,PS3,PE3,  &
                               MS1,ME1,MS2,ME2,MS3,ME3   )
    IMPLICIT NONE
    INTEGER                         , INTENT(IN)  :: noutbuf
    LOGICAL    , DIMENSION(noutbuf) , INTENT(OUT) :: outbuf
    INTEGER   MS1,ME1,MS2,ME2,MS3,ME3
    INTEGER   PS1,PE1,PS2,PE2,PS3,PE3
    LOGICAL    , DIMENSION( MS1:ME1,MS2:ME2,MS3:ME3 ) , INTENT(in) :: inbuf

    INTEGER               :: i,j,k,n   ,  icurs
    icurs = 1
      DO k = PS3, PE3
        DO j = PS2, PE2
          DO i = PS1, PE1
            outbuf( icurs )  = inbuf( i, j, k )
            icurs = icurs + 1
          END DO
        END DO
      END DO
    RETURN
  END SUBROUTINE just_patch_l


  SUBROUTINE patch_2_outbuf_r( inbuf, outbuf,            &
                               DS1,DE1,DS2,DE2,DS3,DE3,  &
                               GPATCH )
    USE module_dm, ONLY : ntasks
    IMPLICIT NONE
    REAL    , DIMENSION(*) , INTENT(IN) :: inbuf
    INTEGER   DS1,DE1,DS2,DE2,DS3,DE3,GPATCH(3,2,ntasks)
    REAL    , DIMENSION( DS1:DE1,DS2:DE2,DS3:DE3 ) , INTENT(out) :: outbuf

    INTEGER               :: i,j,k,n   ,  icurs
    icurs = 1
    DO n = 1, ntasks
      DO k = GPATCH( 3,1,n ), GPATCH( 3,2,n )
        DO j = GPATCH( 2,1,n ), GPATCH( 2,2,n )
          DO i = GPATCH( 1,1,n ), GPATCH( 1,2,n )
            outbuf( i, j, k ) = inbuf( icurs )
            icurs = icurs + 1
          END DO
        END DO
      END DO
    END DO

    RETURN
  END SUBROUTINE patch_2_outbuf_r

  SUBROUTINE patch_2_outbuf_i( inbuf, outbuf,         &
                               DS1,DE1,DS2,DE2,DS3,DE3,&
                               GPATCH )
    USE module_dm, ONLY : ntasks
    IMPLICIT NONE
    INTEGER    , DIMENSION(*) , INTENT(IN) :: inbuf
    INTEGER   DS1,DE1,DS2,DE2,DS3,DE3,GPATCH(3,2,ntasks)
    INTEGER    , DIMENSION( DS1:DE1,DS2:DE2,DS3:DE3 ) , INTENT(out) :: outbuf

    INTEGER               :: i,j,k,n   ,  icurs
    icurs = 1
    DO n = 1, ntasks
      DO k = GPATCH( 3,1,n ), GPATCH( 3,2,n )
        DO j = GPATCH( 2,1,n ), GPATCH( 2,2,n )
          DO i = GPATCH( 1,1,n ), GPATCH( 1,2,n )
            outbuf( i, j, k ) = inbuf( icurs )
            icurs = icurs + 1
          END DO
        END DO
      END DO
    END DO
    RETURN
  END SUBROUTINE patch_2_outbuf_i

  SUBROUTINE patch_2_outbuf_d( inbuf, outbuf,         &
                               DS1,DE1,DS2,DE2,DS3,DE3,&
                               GPATCH )
    USE module_dm, ONLY : ntasks
    IMPLICIT NONE
    DOUBLE PRECISION    , DIMENSION(*) , INTENT(IN) :: inbuf
    INTEGER   DS1,DE1,DS2,DE2,DS3,DE3,GPATCH(3,2,ntasks)
    DOUBLE PRECISION    , DIMENSION( DS1:DE1,DS2:DE2,DS3:DE3 ) , INTENT(out) :: outbuf

    INTEGER               :: i,j,k,n   ,  icurs
    icurs = 1
    DO n = 1, ntasks
      DO k = GPATCH( 3,1,n ), GPATCH( 3,2,n )
        DO j = GPATCH( 2,1,n ), GPATCH( 2,2,n )
          DO i = GPATCH( 1,1,n ), GPATCH( 1,2,n )
            outbuf( i, j, k ) = inbuf( icurs )
            icurs = icurs + 1
          END DO
        END DO
      END DO
    END DO
    RETURN
  END SUBROUTINE patch_2_outbuf_d

  SUBROUTINE patch_2_outbuf_l( inbuf, outbuf,         &
                               DS1,DE1,DS2,DE2,DS3,DE3,&
                               GPATCH )
    USE module_dm, ONLY : ntasks
    IMPLICIT NONE
    LOGICAL    , DIMENSION(*) , INTENT(IN) :: inbuf
    INTEGER   DS1,DE1,DS2,DE2,DS3,DE3,GPATCH(3,2,ntasks)
    LOGICAL    , DIMENSION( DS1:DE1,DS2:DE2,DS3:DE3 ) , INTENT(out) :: outbuf

    INTEGER               :: i,j,k,n   ,  icurs
    icurs = 1
    DO n = 1, ntasks
      DO k = GPATCH( 3,1,n ), GPATCH( 3,2,n )
        DO j = GPATCH( 2,1,n ), GPATCH( 2,2,n )
          DO i = GPATCH( 1,1,n ), GPATCH( 1,2,n )
            outbuf( i, j, k ) = inbuf( icurs )
            icurs = icurs + 1
          END DO
        END DO
      END DO
    END DO
    RETURN
  END SUBROUTINE patch_2_outbuf_l



    SUBROUTINE wrf_global_to_patch_real (globbuf,buf,domdesc,stagger,ordering,&
                                       DS1,DE1,DS2,DE2,DS3,DE3,&
                                       MS1,ME1,MS2,ME2,MS3,ME3,&
                                       PS1,PE1,PS2,PE2,PS3,PE3 )
       IMPLICIT NONE
       INTEGER                         DS1,DE1,DS2,DE2,DS3,DE3,&
                                       MS1,ME1,MS2,ME2,MS3,ME3,&
                                       PS1,PE1,PS2,PE2,PS3,PE3
       CHARACTER *(*) stagger,ordering
       INTEGER fid,domdesc
       REAL globbuf(*)
       REAL buf(*)

       CALL wrf_global_to_patch_generic (globbuf,buf,domdesc,stagger,ordering,4,&
                                       DS1,DE1,DS2,DE2,DS3,DE3,&
                                       MS1,ME1,MS2,ME2,MS3,ME3,&
                                       PS1,PE1,PS2,PE2,PS3,PE3 )
       RETURN
    END SUBROUTINE wrf_global_to_patch_real

    SUBROUTINE wrf_global_to_patch_double (globbuf,buf,domdesc,stagger,ordering,&
                                       DS1,DE1,DS2,DE2,DS3,DE3,&
                                       MS1,ME1,MS2,ME2,MS3,ME3,&
                                       PS1,PE1,PS2,PE2,PS3,PE3 )
       IMPLICIT NONE
       INTEGER                         DS1,DE1,DS2,DE2,DS3,DE3,&
                                       MS1,ME1,MS2,ME2,MS3,ME3,&
                                       PS1,PE1,PS2,PE2,PS3,PE3
       CHARACTER *(*) stagger,ordering
       INTEGER fid,domdesc




       REAL globbuf(*)
       REAL buf(*)

       CALL wrf_global_to_patch_generic (globbuf,buf,domdesc,stagger,ordering,8,&
                                       DS1,DE1,DS2,DE2,DS3,DE3,&
                                       MS1,ME1,MS2,ME2,MS3,ME3,&
                                       PS1,PE1,PS2,PE2,PS3,PE3 )
       RETURN
    END SUBROUTINE wrf_global_to_patch_double


    SUBROUTINE wrf_global_to_patch_integer (globbuf,buf,domdesc,stagger,ordering,&
                                       DS1,DE1,DS2,DE2,DS3,DE3,&
                                       MS1,ME1,MS2,ME2,MS3,ME3,&
                                       PS1,PE1,PS2,PE2,PS3,PE3 )
       IMPLICIT NONE
       INTEGER                         DS1,DE1,DS2,DE2,DS3,DE3,&
                                       MS1,ME1,MS2,ME2,MS3,ME3,&
                                       PS1,PE1,PS2,PE2,PS3,PE3
       CHARACTER *(*) stagger,ordering
       INTEGER fid,domdesc
       INTEGER globbuf(*)
       INTEGER buf(*)

       CALL wrf_global_to_patch_generic (globbuf,buf,domdesc,stagger,ordering,4,&
                                       DS1,DE1,DS2,DE2,DS3,DE3,&
                                       MS1,ME1,MS2,ME2,MS3,ME3,&
                                       PS1,PE1,PS2,PE2,PS3,PE3 )
       RETURN
    END SUBROUTINE wrf_global_to_patch_integer

    SUBROUTINE wrf_global_to_patch_logical (globbuf,buf,domdesc,stagger,ordering,&
                                       DS1,DE1,DS2,DE2,DS3,DE3,&
                                       MS1,ME1,MS2,ME2,MS3,ME3,&
                                       PS1,PE1,PS2,PE2,PS3,PE3 )
       IMPLICIT NONE
       INTEGER                         DS1,DE1,DS2,DE2,DS3,DE3,&
                                       MS1,ME1,MS2,ME2,MS3,ME3,&
                                       PS1,PE1,PS2,PE2,PS3,PE3
       CHARACTER *(*) stagger,ordering
       INTEGER fid,domdesc
       LOGICAL globbuf(*)
       LOGICAL buf(*)

       CALL wrf_global_to_patch_generic (globbuf,buf,domdesc,stagger,ordering,4,&
                                       DS1,DE1,DS2,DE2,DS3,DE3,&
                                       MS1,ME1,MS2,ME2,MS3,ME3,&
                                       PS1,PE1,PS2,PE2,PS3,PE3 )
       RETURN
    END SUBROUTINE wrf_global_to_patch_logical

    SUBROUTINE wrf_global_to_patch_generic (globbuf,buf,domdesc,stagger,ordering,typesize,&
                                       DS1a,DE1a,DS2a,DE2a,DS3a,DE3a,&
                                       MS1a,ME1a,MS2a,ME2a,MS3a,ME3a,&
                                       PS1a,PE1a,PS2a,PE2a,PS3a,PE3a )
       USE module_dm, ONLY : local_communicator, ntasks
       USE module_driver_constants
       IMPLICIT NONE
       INTEGER                         DS1a,DE1a,DS2a,DE2a,DS3a,DE3a,&
                                       MS1a,ME1a,MS2a,ME2a,MS3a,ME3a,&
                                       PS1a,PE1a,PS2a,PE2a,PS3a,PE3A
       CHARACTER *(*) stagger,ordering
       INTEGER domdesc,typesize,ierr
       REAL globbuf(*)
       REAL buf(*)
       INTEGER                         DS1,DE1,DS2,DE2,DS3,DE3,&
                                       MS1,ME1,MS2,ME2,MS3,ME3,&
                                       PS1,PE1,PS2,PE2,PS3,PE3
       LOGICAL, EXTERNAL :: wrf_dm_on_monitor, has_char

       INTEGER i,j,k,ord,ord2d,ndim
       INTEGER  Patch(3,2), Gpatch(3,2,ntasks)
       REAL, ALLOCATABLE :: tmpbuf( : )
       REAL locbuf( (PE1a-PS1a+1)*(PE2a-PS2a+1)*(PE3a-PS3a+1)/4*typesize+32 )

       DS1 = DS1a ; DE1 = DE1a ; DS2=DS2a ; DE2 = DE2a ; DS3 = DS3a ; DE3 = DE3a
       MS1 = MS1a ; ME1 = ME1a ; MS2=MS2a ; ME2 = ME2a ; MS3 = MS3a ; ME3 = ME3a
       PS1 = PS1a ; PE1 = PE1a ; PS2=PS2a ; PE2 = PE2a ; PS3 = PS3a ; PE3 = PE3a

       SELECT CASE ( TRIM(ordering) )
         CASE ( 'xy', 'yx' )
           ndim = 2
         CASE DEFAULT
           ndim = 3   
       END SELECT

       SELECT CASE ( TRIM(ordering) )
         CASE ( 'xyz','xy' )
            
            
            
           IF ( .NOT. has_char( stagger, 'x' ) ) DE1 = DE1+1
           IF ( .NOT. has_char( stagger, 'y' ) ) DE2 = DE2+1
           IF ( ndim .EQ. 3 .AND. .NOT. has_char( stagger, 'z' ) ) DE3 = DE3+1
         CASE ( 'yxz','yx' )
           IF ( .NOT. has_char( stagger, 'x' ) ) DE2 = DE2+1
           IF ( .NOT. has_char( stagger, 'y' ) ) DE1 = DE1+1
           IF ( ndim .EQ. 3 .AND. .NOT. has_char( stagger, 'z' ) ) DE3 = DE3+1
         CASE ( 'zxy' )
           IF ( .NOT. has_char( stagger, 'x' ) ) DE2 = DE2+1
           IF ( .NOT. has_char( stagger, 'y' ) ) DE3 = DE3+1
           IF ( ndim .EQ. 3 .AND. .NOT. has_char( stagger, 'z' ) ) DE1 = DE1+1
         CASE ( 'xzy' )
           IF ( .NOT. has_char( stagger, 'x' ) ) DE1 = DE1+1
           IF ( .NOT. has_char( stagger, 'y' ) ) DE3 = DE3+1
           IF ( ndim .EQ. 3 .AND. .NOT. has_char( stagger, 'z' ) ) DE2 = DE2+1
         CASE DEFAULT
       END SELECT

     
       IF ( wrf_dm_on_monitor() ) THEN
         ALLOCATE ( tmpbuf ( (DE1-DS1+1)*(DE2-DS2+1)*(DE3-DS3+1)/4*typesize+32 ), STAT=ierr )
       ELSE
         ALLOCATE ( tmpbuf ( 1 ), STAT=ierr )
       END IF
       IF ( ierr .ne. 0 ) CALL wrf_error_fatal3("<stdin>",6458,&
'allocating tmpbuf in wrf_global_to_patch_generic')

       Patch(1,1) = ps1 ; Patch(1,2) = pe1    
       Patch(2,1) = ps2 ; Patch(2,2) = pe2
       Patch(3,1) = ps3 ; Patch(3,2) = pe3


       CALL collect_on_comm0 (  local_communicator , 4 ,  &
                                Patch , 6 ,                       &
                                GPatch , 6*ntasks                 )
       ndim = len(TRIM(ordering))

       IF ( wrf_dm_on_monitor() .AND. ndim .GE. 2 ) THEN
         IF      ( typesize .EQ. 4 ) THEN
           CALL outbuf_2_patch_r ( globbuf , tmpbuf  ,    &
                                   DS1, DE1, DS2, DE2, DS3, DE3 , &
                                   MS1, ME1, MS2, ME2, MS3, ME3 , &
                                   GPATCH                         )
         ELSE IF ( typesize .EQ. 4 ) THEN
           CALL outbuf_2_patch_i ( globbuf , tmpbuf  ,    &
                                   DS1, DE1, DS2, DE2, DS3, DE3 , &
                                   GPATCH                         )
         ELSE IF ( typesize .EQ. 8 ) THEN
           CALL outbuf_2_patch_d ( globbuf , tmpbuf  ,    &
                                   DS1, DE1, DS2, DE2, DS3, DE3 , &
                                   GPATCH                         )
         ELSE IF ( typesize .EQ. 4 ) THEN
           CALL outbuf_2_patch_l ( globbuf , tmpbuf  ,    &
                                   DS1, DE1, DS2, DE2, DS3, DE3 , &
                                   GPATCH                         )
         END IF
       END IF

       CALL dist_on_comm0 (  local_communicator , typesize ,  &
                             tmpbuf  , (de1-ds1+1)*(de2-ds2+1)*(de3-ds3+1) , &
                             locbuf    , (pe1-ps1+1)*(pe2-ps2+1)*(pe3-ps3+1) )

       IF      ( typesize .EQ. 4 ) THEN
         CALL all_sub_r ( locbuf , buf ,             &
                                   PS1, PE1, PS2, PE2, PS3, PE3 , &
                                   MS1, ME1, MS2, ME2, MS3, ME3   )

       ELSE IF ( typesize .EQ. 4 ) THEN
         CALL all_sub_i ( locbuf , buf ,             &
                                   PS1, PE1, PS2, PE2, PS3, PE3 , &
                                   MS1, ME1, MS2, ME2, MS3, ME3   )
       ELSE IF ( typesize .EQ. 8 ) THEN
         CALL all_sub_d ( locbuf , buf ,             &
                                   PS1, PE1, PS2, PE2, PS3, PE3 , &
                                   MS1, ME1, MS2, ME2, MS3, ME3   )
       ELSE IF ( typesize .EQ. 4 ) THEN
         CALL all_sub_l ( locbuf , buf ,             &
                                   PS1, PE1, PS2, PE2, PS3, PE3 , &
                                   MS1, ME1, MS2, ME2, MS3, ME3   )
       END IF


       DEALLOCATE ( tmpbuf )
       RETURN
    END SUBROUTINE wrf_global_to_patch_generic

  SUBROUTINE all_sub_i ( inbuf , outbuf,              &
                               PS1,PE1,PS2,PE2,PS3,PE3,  &
                               MS1,ME1,MS2,ME2,MS3,ME3   )
    IMPLICIT NONE
    INTEGER    , DIMENSION(*) , INTENT(IN) :: inbuf
    INTEGER   MS1,ME1,MS2,ME2,MS3,ME3
    INTEGER   PS1,PE1,PS2,PE2,PS3,PE3
    INTEGER    , DIMENSION( MS1:ME1,MS2:ME2,MS3:ME3 ) , INTENT(OUT) :: outbuf

    INTEGER               :: i,j,k,n   ,  icurs
    icurs = 1
      DO k = PS3, PE3
        DO j = PS2, PE2
          DO i = PS1, PE1
            outbuf( i, j, k )  = inbuf ( icurs )
            icurs = icurs + 1
          END DO
        END DO
      END DO
    RETURN
  END SUBROUTINE all_sub_i

  SUBROUTINE all_sub_r ( inbuf , outbuf,              &
                               PS1,PE1,PS2,PE2,PS3,PE3,  &
                               MS1,ME1,MS2,ME2,MS3,ME3   )
    IMPLICIT NONE
    REAL       , DIMENSION(*) , INTENT(IN) :: inbuf
    INTEGER   MS1,ME1,MS2,ME2,MS3,ME3
    INTEGER   PS1,PE1,PS2,PE2,PS3,PE3
    REAL       , DIMENSION( MS1:ME1,MS2:ME2,MS3:ME3 ) , INTENT(OUT) :: outbuf

    INTEGER               :: i,j,k,n   ,  icurs
    icurs = 1
      DO k = PS3, PE3
        DO j = PS2, PE2
          DO i = PS1, PE1
            outbuf( i, j, k )  = inbuf ( icurs )
            icurs = icurs + 1
          END DO
        END DO
      END DO

    RETURN
  END SUBROUTINE all_sub_r

  SUBROUTINE all_sub_d ( inbuf , outbuf,              &
                               PS1,PE1,PS2,PE2,PS3,PE3,  &
                               MS1,ME1,MS2,ME2,MS3,ME3   )
    IMPLICIT NONE
    DOUBLE PRECISION    , DIMENSION(*) , INTENT(IN) :: inbuf
    INTEGER   MS1,ME1,MS2,ME2,MS3,ME3
    INTEGER   PS1,PE1,PS2,PE2,PS3,PE3
    DOUBLE PRECISION    , DIMENSION( MS1:ME1,MS2:ME2,MS3:ME3 ) , INTENT(OUT) :: outbuf

    INTEGER               :: i,j,k,n   ,  icurs
    icurs = 1
      DO k = PS3, PE3
        DO j = PS2, PE2
          DO i = PS1, PE1
            outbuf( i, j, k )  = inbuf ( icurs )
            icurs = icurs + 1
          END DO
        END DO
      END DO
    RETURN
  END SUBROUTINE all_sub_d

  SUBROUTINE all_sub_l ( inbuf , outbuf,              &
                               PS1,PE1,PS2,PE2,PS3,PE3,  &
                               MS1,ME1,MS2,ME2,MS3,ME3   )
    IMPLICIT NONE
    LOGICAL    , DIMENSION(*) , INTENT(IN) :: inbuf
    INTEGER   MS1,ME1,MS2,ME2,MS3,ME3
    INTEGER   PS1,PE1,PS2,PE2,PS3,PE3
    LOGICAL    , DIMENSION( MS1:ME1,MS2:ME2,MS3:ME3 ) , INTENT(OUT) :: outbuf

    INTEGER               :: i,j,k,n   ,  icurs
    icurs = 1
      DO k = PS3, PE3
        DO j = PS2, PE2
          DO i = PS1, PE1
            outbuf( i, j, k )  = inbuf ( icurs )
            icurs = icurs + 1
          END DO
        END DO
      END DO
    RETURN
  END SUBROUTINE all_sub_l

  SUBROUTINE outbuf_2_patch_r( inbuf, outbuf,         &
                               DS1,DE1,DS2,DE2,DS3,DE3, &
                               MS1, ME1, MS2, ME2, MS3, ME3 , &
                               GPATCH )
    USE module_dm, ONLY : ntasks
    IMPLICIT NONE
    REAL    , DIMENSION(*) , INTENT(OUT) :: outbuf
    INTEGER   DS1,DE1,DS2,DE2,DS3,DE3,GPATCH(3,2,ntasks)
    INTEGER   MS1,ME1,MS2,ME2,MS3,ME3
    REAL    , DIMENSION( DS1:DE1,DS2:DE2,DS3:DE3 ) , INTENT(IN) :: inbuf

    INTEGER               :: i,j,k,n   ,  icurs

    icurs = 1
    DO n = 1, ntasks
      DO k = GPATCH( 3,1,n ), GPATCH( 3,2,n )
        DO j = GPATCH( 2,1,n ), GPATCH( 2,2,n )
          DO i = GPATCH( 1,1,n ), GPATCH( 1,2,n )
            outbuf( icurs ) = inbuf( i,j,k )
            icurs = icurs + 1
          END DO
        END DO
      END DO
    END DO
    RETURN
  END SUBROUTINE outbuf_2_patch_r

  SUBROUTINE outbuf_2_patch_i( inbuf, outbuf,         &
                               DS1,DE1,DS2,DE2,DS3,DE3,&
                               GPATCH )
    USE module_dm, ONLY : ntasks
    IMPLICIT NONE
    INTEGER    , DIMENSION(*) , INTENT(OUT) :: outbuf
    INTEGER   DS1,DE1,DS2,DE2,DS3,DE3,GPATCH(3,2,ntasks)
    INTEGER    , DIMENSION( DS1:DE1,DS2:DE2,DS3:DE3 ) , INTENT(IN) :: inbuf

    INTEGER               :: i,j,k,n   ,  icurs
    icurs = 1
    DO n = 1, ntasks
      DO k = GPATCH( 3,1,n ), GPATCH( 3,2,n )
        DO j = GPATCH( 2,1,n ), GPATCH( 2,2,n )
          DO i = GPATCH( 1,1,n ), GPATCH( 1,2,n )
            outbuf( icurs ) = inbuf( i,j,k )
            icurs = icurs + 1
          END DO
        END DO
      END DO
    END DO
    RETURN
  END SUBROUTINE outbuf_2_patch_i

  SUBROUTINE outbuf_2_patch_d( inbuf, outbuf,         &
                               DS1,DE1,DS2,DE2,DS3,DE3,&
                               GPATCH )
    USE module_dm, ONLY : ntasks
    IMPLICIT NONE
    DOUBLE PRECISION    , DIMENSION(*) , INTENT(OUT) :: outbuf
    INTEGER   DS1,DE1,DS2,DE2,DS3,DE3,GPATCH(3,2,ntasks)
    DOUBLE PRECISION    , DIMENSION( DS1:DE1,DS2:DE2,DS3:DE3 ) , INTENT(IN) :: inbuf

    INTEGER               :: i,j,k,n   ,  icurs
    icurs = 1
    DO n = 1, ntasks
      DO k = GPATCH( 3,1,n ), GPATCH( 3,2,n )
        DO j = GPATCH( 2,1,n ), GPATCH( 2,2,n )
          DO i = GPATCH( 1,1,n ), GPATCH( 1,2,n )
            outbuf( icurs ) = inbuf( i,j,k )
            icurs = icurs + 1
          END DO
        END DO
      END DO
    END DO
    RETURN
  END SUBROUTINE outbuf_2_patch_d

  SUBROUTINE outbuf_2_patch_l( inbuf, outbuf,         &
                               DS1,DE1,DS2,DE2,DS3,DE3,&
                               GPATCH )
    USE module_dm, ONLY : ntasks
    IMPLICIT NONE
    LOGICAL    , DIMENSION(*) , INTENT(OUT) :: outbuf
    INTEGER   DS1,DE1,DS2,DE2,DS3,DE3,GPATCH(3,2,ntasks)
    LOGICAL    , DIMENSION( DS1:DE1,DS2:DE2,DS3:DE3 ) , INTENT(IN) :: inbuf

    INTEGER               :: i,j,k,n   ,  icurs
    icurs = 1
    DO n = 1, ntasks
      DO k = GPATCH( 3,1,n ), GPATCH( 3,2,n )
        DO j = GPATCH( 2,1,n ), GPATCH( 2,2,n )
          DO i = GPATCH( 1,1,n ), GPATCH( 1,2,n )
            outbuf( icurs ) = inbuf( i,j,k )
            icurs = icurs + 1
          END DO
        END DO
      END DO
    END DO
    RETURN
  END SUBROUTINE outbuf_2_patch_l


  SUBROUTINE wrf_dm_nestexchange_init
      CALL rsl_lite_nesting_reset
  END SUBROUTINE wrf_dm_nestexchange_init









   SUBROUTINE before_interp_halos_nmm(grid,config_flags &







,szj,s1z,spz,tcs,moist,moist_bxs,moist_bxe,moist_bys,moist_bye,moist_btxs,moist_btxe,moist_btys,moist_btye,dfi_moist, &
dfi_moist_bxs,dfi_moist_bxe,dfi_moist_bys,dfi_moist_bye,dfi_moist_btxs,dfi_moist_btxe,dfi_moist_btys,dfi_moist_btye,scalar, &
scalar_bxs,scalar_bxe,scalar_bys,scalar_bye,scalar_btxs,scalar_btxe,scalar_btys,scalar_btye,dfi_scalar,dfi_scalar_bxs, &
dfi_scalar_bxe,dfi_scalar_bys,dfi_scalar_bye,dfi_scalar_btxs,dfi_scalar_btxe,dfi_scalar_btys,dfi_scalar_btye,chem,ozmixm &


      )
     
     
     
      USE module_state_description
      USE module_domain, ONLY : domain, get_ijk_from_grid
      USE module_configure, ONLY : grid_config_rec_type
      USE module_dm, ONLY : ntasks, ntasks_x, ntasks_y, itrace, local_communicator, mytask, &
                            ipe_save, jpe_save, ips_save, jps_save, get_dm_max_halo_width,  &
                            nest_pes_x, nest_pes_y 
                             
      USE module_comm_dm, ONLY : HALO_NMM_WEIGHTS_sub

      IMPLICIT NONE

      TYPE(domain), POINTER :: grid          
      TYPE (grid_config_rec_type)            :: config_flags






real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_szj)           :: szj
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_s1z)           :: s1z
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_spz)           :: spz
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_tcs)           :: tcs
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_moist)           :: moist
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_bxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_bxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_bys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_bye
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_btxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_btxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_btys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_btye
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_dfi_moist)           :: dfi_moist
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_bxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_bxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_bys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_bye
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_btxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_btxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_btys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_btye
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_scalar)           :: scalar
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_bxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_bxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_bys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_bye
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_btxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_btxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_btys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_btye
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_dfi_scalar)           :: dfi_scalar
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_bxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_bxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_bys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_bye
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_btxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_btxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_btys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_btye
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%sm32:grid%em32,num_chem)           :: chem
real      ,DIMENSION(grid%sm31:grid%em31,1:grid%levsiz,grid%sm32:grid%em32,num_ozmixm)           :: ozmixm

     INTEGER :: IDS,IDE,JDS,JDE,KDS,KDE, &
                IMS,IME,JMS,JME,KMS,KME, &
                IPS,IPE,JPS,JPE,KPS,KPE












     
      IDS=-1; IDE=-1; JDS=-1; JDE=-1; KDS=-1; KDE=-1
      IMS=-1; IME=-1; JMS=-1; JME=-1; KMS=-1; KME=-1
      IPS=-1; IPE=-1; JPS=-1; JPE=-1; KPS=-1; KPE=-1
      CALL GET_IJK_FROM_GRID(GRID                                       &
     &                      ,IDS,IDE,JDS,JDE,KDS,KDE                    &
     &                      ,IMS,IME,JMS,JME,KMS,KME                    &
     &                      ,IPS,IPE,JPS,JPE,KPS,KPE )

     CALL push_communicators_for_domain(grid%id)






CALL HALO_NMM_WEIGHTS_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )

     CALL pop_communicators_for_domain

   END SUBROUTINE before_interp_halos_nmm

   SUBROUTINE interp_domain_nmm_part1 ( grid, intermediate_grid, ngrid, config_flags    &







,szj,s1z,spz,tcs,moist,moist_bxs,moist_bxe,moist_bys,moist_bye,moist_btxs,moist_btxe,moist_btys,moist_btye,dfi_moist, &
dfi_moist_bxs,dfi_moist_bxe,dfi_moist_bys,dfi_moist_bye,dfi_moist_btxs,dfi_moist_btxe,dfi_moist_btys,dfi_moist_btye,scalar, &
scalar_bxs,scalar_bxe,scalar_bys,scalar_bye,scalar_btxs,scalar_btxe,scalar_btys,scalar_btye,dfi_scalar,dfi_scalar_bxs, &
dfi_scalar_bxe,dfi_scalar_bys,dfi_scalar_bye,dfi_scalar_btxs,dfi_scalar_btxe,dfi_scalar_btys,dfi_scalar_btye,chem,ozmixm &


                 )
      USE module_state_description
      USE module_domain, ONLY : domain, get_ijk_from_grid
      USE module_configure, ONLY : grid_config_rec_type
      USE module_dm, ONLY : ntasks, ntasks_x, ntasks_y, itrace, local_communicator, mytask, &
                            ipe_save, jpe_save, ips_save, jps_save, get_dm_max_halo_width,  &
                            nest_pes_x, nest_pes_y
      USE module_dm, ONLY : ntasks, ntasks_x, ntasks_y, local_communicator, mytask, &
                            nest_pes_x, nest_pes_y,                                 &
                            intercomm_active, nest_task_offsets,                    &
                            mpi_comm_to_mom, mpi_comm_to_kid, which_kid 
                             

      USE module_timing
      IMPLICIT NONE

      TYPE(domain), POINTER :: grid          
      TYPE(domain), POINTER :: intermediate_grid
      TYPE(domain), POINTER :: ngrid






real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_szj)           :: szj
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_s1z)           :: s1z
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_spz)           :: spz
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_tcs)           :: tcs
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_moist)           :: moist
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_bxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_bxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_bys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_bye
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_btxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_btxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_btys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_btye
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_dfi_moist)           :: dfi_moist
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_bxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_bxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_bys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_bye
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_btxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_btxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_btys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_btye
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_scalar)           :: scalar
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_bxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_bxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_bys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_bye
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_btxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_btxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_btys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_btye
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_dfi_scalar)           :: dfi_scalar
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_bxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_bxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_bys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_bye
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_btxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_btxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_btys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_btye
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%sm32:grid%em32,num_chem)           :: chem
real      ,DIMENSION(grid%sm31:grid%em31,1:grid%levsiz,grid%sm32:grid%em32,num_ozmixm)           :: ozmixm

      INTEGER nlev, msize
      INTEGER i,j,pig,pjg,cm,cn,nig,njg,retval,k
      INTEGER iparstrt,jparstrt,sw
      TYPE (grid_config_rec_type)            :: config_flags
      REAL xv(2000000)
      INTEGER ixv
      INTEGER       ::          cids, cide, cjds, cjde, ckds, ckde,    &
                                cims, cime, cjms, cjme, ckms, ckme,    &
                                cips, cipe, cjps, cjpe, ckps, ckpe
      INTEGER       ::          iids, iide, ijds, ijde, ikds, ikde,    &
                                iims, iime, ijms, ijme, ikms, ikme,    &
                                iips, iipe, ijps, ijpe, ikps, ikpe
      INTEGER       ::          nids, nide, njds, njde, nkds, nkde,    &
                                nims, nime, njms, njme, nkms, nkme,    &
                                nips, nipe, njps, njpe, nkps, nkpe

      INTEGER idim1,idim2,idim3,idim4,idim5,idim6,idim7
      LOGICAL feedback_flag, feedback_flag_v
      INTEGER icoord, jcoord, idim_cd, jdim_cd, pgr
      INTEGER local_comm, ioffset, myproc, nproc, ierr
      INTEGER thisdomain_max_halo_width

      LOGICAL interp_mp
      interp_mp=grid%interp_mp .or. ngrid%interp_mp

      CALL wrf_get_dm_communicator ( local_comm )
      CALL wrf_get_myproc( myproc )
      CALL wrf_get_nproc( nproc )

      CALL get_ijk_from_grid (  grid ,                   &
                                cids, cide, cjds, cjde, ckds, ckde,    &
                                cims, cime, cjms, cjme, ckms, ckme,    &
                                cips, cipe, cjps, cjpe, ckps, ckpe    )
      CALL get_ijk_from_grid (  intermediate_grid ,              &
                                iids, iide, ijds, ijde, ikds, ikde,    &
                                iims, iime, ijms, ijme, ikms, ikme,    &
                                iips, iipe, ijps, ijpe, ikps, ikpe    )
      CALL get_ijk_from_grid (  ngrid ,              &
                                nids, nide, njds, njde, nkds, nkde,    &
                                nims, nime, njms, njme, nkms, nkme,    &
                                nips, nipe, njps, njpe, nkps, nkpe    )

      CALL nl_get_parent_grid_ratio ( ngrid%id, pgr )
      CALL nl_get_i_parent_start ( intermediate_grid%id, iparstrt )
      CALL nl_get_j_parent_start ( intermediate_grid%id, jparstrt )
      CALL nl_get_shw            ( intermediate_grid%id, sw )
      icoord =    iparstrt - sw
      jcoord =    jparstrt - sw
      idim_cd = iide - iids + 1
      jdim_cd = ijde - ijds + 1

      nlev  = ckde - ckds + 1

      
      CALL get_dm_max_halo_width ( ngrid%id , thisdomain_max_halo_width )

      IF ( grid%active_this_task ) THEN






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
CALL rsl_lite_to_child_info( msize*4                               &
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
ixv=0 
IF ( SIZE(grid%lakedepth2d) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%lakedepth2d(pig,pjg)
ENDIF
IF ( SIZE(grid%savedtke12d) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%savedtke12d(pig,pjg)
ENDIF
IF ( SIZE(grid%snowdp2d) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%snowdp2d(pig,pjg)
ENDIF
IF ( SIZE(grid%h2osno2d) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%h2osno2d(pig,pjg)
ENDIF
IF ( SIZE(grid%snl2d) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%snl2d(pig,pjg)
ENDIF
IF ( SIZE(grid%t_grnd2d) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%t_grnd2d(pig,pjg)
ENDIF
IF ( SIZE(grid%t_lake3d) .GT. 1 ) THEN 
DO k = 1,10
xv(k+ixv+1-1)= grid%t_lake3d(pig,k,pjg)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%lake_icefrac3d) .GT. 1 ) THEN 
DO k = 1,10
xv(k+ixv+1-1)= grid%lake_icefrac3d(pig,k,pjg)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%z_lake3d) .GT. 1 ) THEN 
DO k = 1,10
xv(k+ixv+1-1)= grid%z_lake3d(pig,k,pjg)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%dz_lake3d) .GT. 1 ) THEN 
DO k = 1,10
xv(k+ixv+1-1)= grid%dz_lake3d(pig,k,pjg)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%t_soisno3d) .GT. 1 ) THEN 
DO k = 1,15
xv(k+ixv+1-1)= grid%t_soisno3d(pig,k,pjg)
ENDDO
ixv=ixv+15-1+1
ENDIF
IF ( SIZE(grid%h2osoi_ice3d) .GT. 1 ) THEN 
DO k = 1,15
xv(k+ixv+1-1)= grid%h2osoi_ice3d(pig,k,pjg)
ENDDO
ixv=ixv+15-1+1
ENDIF
IF ( SIZE(grid%h2osoi_liq3d) .GT. 1 ) THEN 
DO k = 1,15
xv(k+ixv+1-1)= grid%h2osoi_liq3d(pig,k,pjg)
ENDDO
ixv=ixv+15-1+1
ENDIF
IF ( SIZE(grid%h2osoi_vol3d) .GT. 1 ) THEN 
DO k = 1,15
xv(k+ixv+1-1)= grid%h2osoi_vol3d(pig,k,pjg)
ENDDO
ixv=ixv+15-1+1
ENDIF
IF ( SIZE(grid%z3d) .GT. 1 ) THEN 
DO k = 1,15
xv(k+ixv+1-1)= grid%z3d(pig,k,pjg)
ENDDO
ixv=ixv+15-1+1
ENDIF
IF ( SIZE(grid%dz3d) .GT. 1 ) THEN 
DO k = 1,15
xv(k+ixv+1-1)= grid%dz3d(pig,k,pjg)
ENDDO
ixv=ixv+15-1+1
ENDIF
IF ( SIZE(grid%zi3d) .GT. 1 ) THEN 
DO k = 1,16
xv(k+ixv+1-1)= grid%zi3d(pig,k,pjg)
ENDDO
ixv=ixv+16-1+1
ENDIF
IF ( SIZE(grid%watsat3d) .GT. 1 ) THEN 
DO k = 1,10
xv(k+ixv+1-1)= grid%watsat3d(pig,k,pjg)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%csol3d) .GT. 1 ) THEN 
DO k = 1,10
xv(k+ixv+1-1)= grid%csol3d(pig,k,pjg)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%tkmg3d) .GT. 1 ) THEN 
DO k = 1,10
xv(k+ixv+1-1)= grid%tkmg3d(pig,k,pjg)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%tkdry3d) .GT. 1 ) THEN 
DO k = 1,10
xv(k+ixv+1-1)= grid%tkdry3d(pig,k,pjg)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%tksatu3d) .GT. 1 ) THEN 
DO k = 1,10
xv(k+ixv+1-1)= grid%tksatu3d(pig,k,pjg)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%vertstrucc) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
xv(k+ixv+1-ckds)= grid%vertstrucc(pig,pjg,k)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%vertstrucs) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
xv(k+ixv+1-ckds)= grid%vertstrucs(pig,pjg,k)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%ru_tendf_stoch) .GT. 1 ) THEN 
DO k = 1,config_flags%num_stoch_levels
xv(k+ixv+1-1)= grid%ru_tendf_stoch(pig,pjg,k)
ENDDO
ixv=ixv+config_flags%num_stoch_levels-1+1
ENDIF
IF ( SIZE(grid%rv_tendf_stoch) .GT. 1 ) THEN 
DO k = 1,config_flags%num_stoch_levels
xv(k+ixv+1-1)= grid%rv_tendf_stoch(pig,pjg,k)
ENDDO
ixv=ixv+config_flags%num_stoch_levels-1+1
ENDIF
IF ( SIZE(grid%lu_index) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%lu_index(pig,pjg)
ENDIF
IF ( SIZE(grid%lu_mask) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%lu_mask(pig,pjg)
ENDIF
IF ( SIZE(grid%cuprecip_swath) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%cuprecip_swath(pig,pjg)
ENDIF
IF ( SIZE(grid%precip_swath) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%precip_swath(pig,pjg)
ENDIF
IF ( SIZE(grid%windsq_swath) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%windsq_swath(pig,pjg)
ENDIF
IF ( SIZE(grid%suswind) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%suswind(pig,pjg)
ENDIF
IF ( SIZE(grid%wind10_ratio) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%wind10_ratio(pig,pjg)
ENDIF
IF ( SIZE(grid%suswind_swath) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%suswind_swath(pig,pjg)
ENDIF
IF ( SIZE(grid%membrane_mslp) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%membrane_mslp(pig,pjg)
ENDIF
IF ( SIZE(grid%pdyn_smooth) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%pdyn_smooth(pig,pjg)
ENDIF
IF ( SIZE(grid%pdyn_parent) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%pdyn_parent(pig,pjg)
ENDIF
IF ( SIZE(grid%best_mslp) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%best_mslp(pig,pjg)
ENDIF
IF ( SIZE(grid%hlon) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%hlon(pig,pjg)
ENDIF
IF ( SIZE(grid%hlat) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%hlat(pig,pjg)
ENDIF
IF ( SIZE(grid%tg_max_m10wind) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%tg_max_m10wind(pig,pjg)
ENDIF
IF ( SIZE(grid%tg_max_wwind) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%tg_max_wwind(pig,pjg)
ENDIF
IF ( SIZE(grid%tg_min_wwind) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%tg_min_wwind(pig,pjg)
ENDIF
IF ( SIZE(grid%tg_max_zhel_25) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%tg_max_zhel_25(pig,pjg)
ENDIF
IF ( SIZE(grid%tg_min_zhel_25) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%tg_min_zhel_25(pig,pjg)
ENDIF
IF ( SIZE(grid%tg_max_zhel_03) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%tg_max_zhel_03(pig,pjg)
ENDIF
IF ( SIZE(grid%tg_min_zhel_03) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%tg_min_zhel_03(pig,pjg)
ENDIF
IF ( SIZE(grid%tg_updhel25) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%tg_updhel25(pig,pjg)
ENDIF
IF ( SIZE(grid%tg_max_updhel25) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%tg_max_updhel25(pig,pjg)
ENDIF
IF ( SIZE(grid%tg_updhel03) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%tg_updhel03(pig,pjg)
ENDIF
IF ( SIZE(grid%tg_max_updhel03) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%tg_max_updhel03(pig,pjg)
ENDIF
IF ( SIZE(grid%tg_total_precip) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%tg_total_precip(pig,pjg)
ENDIF
IF ( SIZE(grid%hres_fis) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%hres_fis(pig,pjg)
ENDIF
IF ( SIZE(grid%hres_lnd) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%hres_lnd(pig,pjg)
ENDIF
IF ( SIZE(grid%hres_lu_index) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%hres_lu_index(pig,pjg)
ENDIF
IF ( SIZE(grid%sm) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%sm(pig,pjg)
ENDIF
IF ( SIZE(grid%sice) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%sice(pig,pjg)
ENDIF
IF ( SIZE(grid%pd) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%pd(pig,pjg)
ENDIF
IF ( SIZE(grid%fis) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%fis(pig,pjg)
ENDIF
IF ( SIZE(grid%t) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
xv(k+ixv+1-ckds)= grid%t(pig,pjg,k)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%q) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
xv(k+ixv+1-ckds)= grid%q(pig,pjg,k)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%u) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
xv(k+ixv+1-ckds)= grid%u(pig,pjg,k)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%v) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
xv(k+ixv+1-ckds)= grid%v(pig,pjg,k)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%ustar) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%ustar(pig,pjg)
ENDIF
IF ( SIZE(grid%z0) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%z0(pig,pjg)
ENDIF
IF ( SIZE(grid%scurx) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%scurx(pig,pjg)
ENDIF
IF ( SIZE(grid%scury) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%scury(pig,pjg)
ENDIF
IF ( SIZE(grid%charn) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%charn(pig,pjg)
ENDIF
IF ( SIZE(grid%msang) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%msang(pig,pjg)
ENDIF
IF ( SIZE(grid%ths) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%ths(pig,pjg)
ENDIF
IF ( SIZE(grid%qsh) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%qsh(pig,pjg)
ENDIF
IF ( SIZE(grid%taux) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%taux(pig,pjg)
ENDIF
IF ( SIZE(grid%tauy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%tauy(pig,pjg)
ENDIF
IF ( SIZE(grid%acprec) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%acprec(pig,pjg)
ENDIF
IF ( SIZE(grid%cuprec) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%cuprec(pig,pjg)
ENDIF
IF ( SIZE(grid%cldefi) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%cldefi(pig,pjg)
ENDIF
IF ( SIZE(grid%th10) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%th10(pig,pjg)
ENDIF
IF ( SIZE(grid%q10) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%q10(pig,pjg)
ENDIF
IF ( SIZE(grid%pshltr) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%pshltr(pig,pjg)
ENDIF
IF ( SIZE(grid%tshltr) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%tshltr(pig,pjg)
ENDIF
IF ( SIZE(grid%qshltr) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%qshltr(pig,pjg)
ENDIF
IF ( SIZE(grid%q2) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
xv(k+ixv+1-ckds)= grid%q2(pig,pjg,k)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%t_adj) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
xv(k+ixv+1-ckds)= grid%t_adj(pig,pjg,k)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%cd_out) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%cd_out(pig,pjg)
ENDIF
IF ( SIZE(grid%ch_out) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%ch_out(pig,pjg)
ENDIF
IF ( SIZE(grid%albase) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%albase(pig,pjg)
ENDIF
IF ( SIZE(grid%nmm_tsk) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%nmm_tsk(pig,pjg)
ENDIF
IF ( SIZE(grid%mxsnal) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%mxsnal(pig,pjg)
ENDIF
IF ( SIZE(grid%sigt4) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%sigt4(pig,pjg)
ENDIF
IF ( SIZE(grid%tg) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%tg(pig,pjg)
ENDIF
IF ( SIZE(grid%acswupt) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%acswupt(pig,pjg)
ENDIF
IF ( SIZE(grid%acswuptc) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%acswuptc(pig,pjg)
ENDIF
IF ( SIZE(grid%acswdnt) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%acswdnt(pig,pjg)
ENDIF
IF ( SIZE(grid%acswdntc) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%acswdntc(pig,pjg)
ENDIF
IF ( SIZE(grid%acswupb) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%acswupb(pig,pjg)
ENDIF
IF ( SIZE(grid%acswupbc) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%acswupbc(pig,pjg)
ENDIF
IF ( SIZE(grid%acswdnb) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%acswdnb(pig,pjg)
ENDIF
IF ( SIZE(grid%acswdnbc) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%acswdnbc(pig,pjg)
ENDIF
IF ( SIZE(grid%aclwupt) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%aclwupt(pig,pjg)
ENDIF
IF ( SIZE(grid%aclwuptc) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%aclwuptc(pig,pjg)
ENDIF
IF ( SIZE(grid%aclwdnt) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%aclwdnt(pig,pjg)
ENDIF
IF ( SIZE(grid%aclwdntc) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%aclwdntc(pig,pjg)
ENDIF
IF ( SIZE(grid%aclwupb) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%aclwupb(pig,pjg)
ENDIF
IF ( SIZE(grid%aclwupbc) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%aclwupbc(pig,pjg)
ENDIF
IF ( SIZE(grid%aclwdnb) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%aclwdnb(pig,pjg)
ENDIF
IF ( SIZE(grid%aclwdnbc) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%aclwdnbc(pig,pjg)
ENDIF
IF ( SIZE(grid%swupt) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%swupt(pig,pjg)
ENDIF
IF ( SIZE(grid%swuptc) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%swuptc(pig,pjg)
ENDIF
IF ( SIZE(grid%swdnt) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%swdnt(pig,pjg)
ENDIF
IF ( SIZE(grid%swdntc) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%swdntc(pig,pjg)
ENDIF
IF ( SIZE(grid%swupb) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%swupb(pig,pjg)
ENDIF
IF ( SIZE(grid%swupbc) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%swupbc(pig,pjg)
ENDIF
IF ( SIZE(grid%swdnb) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%swdnb(pig,pjg)
ENDIF
IF ( SIZE(grid%swdnbc) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%swdnbc(pig,pjg)
ENDIF
IF ( SIZE(grid%lwupt) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%lwupt(pig,pjg)
ENDIF
IF ( SIZE(grid%lwuptc) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%lwuptc(pig,pjg)
ENDIF
IF ( SIZE(grid%lwdnt) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%lwdnt(pig,pjg)
ENDIF
IF ( SIZE(grid%lwdntc) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%lwdntc(pig,pjg)
ENDIF
IF ( SIZE(grid%lwupb) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%lwupb(pig,pjg)
ENDIF
IF ( SIZE(grid%lwupbc) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%lwupbc(pig,pjg)
ENDIF
IF ( SIZE(grid%lwdnb) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%lwdnb(pig,pjg)
ENDIF
IF ( SIZE(grid%lwdnbc) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%lwdnbc(pig,pjg)
ENDIF
IF ( SIZE(grid%refl_10cm) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
xv(k+ixv+1-ckds)= grid%refl_10cm(pig,k,pjg)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%refd_max) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%refd_max(pig,pjg)
ENDIF
IF ( SIZE(grid%qnwfa2d) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%qnwfa2d(pig,pjg)
ENDIF
IF ( SIZE(grid%swddir) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%swddir(pig,pjg)
ENDIF
IF ( SIZE(grid%swddirc) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%swddirc(pig,pjg)
ENDIF
IF ( SIZE(grid%swddni) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%swddni(pig,pjg)
ENDIF
IF ( SIZE(grid%swddnic) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%swddnic(pig,pjg)
ENDIF
IF ( SIZE(grid%swddif) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%swddif(pig,pjg)
ENDIF
IF ( SIZE(grid%gx) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%gx(pig,pjg)
ENDIF
IF ( SIZE(grid%bx) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%bx(pig,pjg)
ENDIF
IF ( SIZE(grid%gg) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%gg(pig,pjg)
ENDIF
IF ( SIZE(grid%bb) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%bb(pig,pjg)
ENDIF
IF ( SIZE(grid%coszen_ref) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%coszen_ref(pig,pjg)
ENDIF
IF ( SIZE(grid%swdown_ref) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%swdown_ref(pig,pjg)
ENDIF
IF ( SIZE(grid%swddir_ref) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%swddir_ref(pig,pjg)
ENDIF
IF ( SIZE(grid%cwm) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
xv(k+ixv+1-ckds)= grid%cwm(pig,pjg,k)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%f_ice) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
xv(k+ixv+1-ckds)= grid%f_ice(pig,k,pjg)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%f_rain) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
xv(k+ixv+1-ckds)= grid%f_rain(pig,k,pjg)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%f_rimef) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
xv(k+ixv+1-ckds)= grid%f_rimef(pig,k,pjg)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%cfrach) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%cfrach(pig,pjg)
ENDIF
IF ( SIZE(grid%cfracl) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%cfracl(pig,pjg)
ENDIF
IF ( SIZE(grid%cfracm) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%cfracm(pig,pjg)
ENDIF
IF ( SIZE(grid%islope) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%islope(pig,pjg)
ENDIF
IF ( SIZE(grid%cmc) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%cmc(pig,pjg)
ENDIF
IF ( SIZE(grid%soiltb) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%soiltb(pig,pjg)
ENDIF
IF ( SIZE(grid%vegfrc) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%vegfrc(pig,pjg)
ENDIF
IF ( SIZE(grid%shdmax) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%shdmax(pig,pjg)
ENDIF
IF ( SIZE(grid%shdmin) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%shdmin(pig,pjg)
ENDIF
IF ( SIZE(grid%sh2o) .GT. 1 ) THEN 
DO k = 1,config_flags%num_soil_layers
xv(k+ixv+1-1)= grid%sh2o(pig,k,pjg)
ENDDO
ixv=ixv+config_flags%num_soil_layers-1+1
ENDIF
IF ( SIZE(grid%smc) .GT. 1 ) THEN 
DO k = 1,config_flags%num_soil_layers
xv(k+ixv+1-1)= grid%smc(pig,k,pjg)
ENDDO
ixv=ixv+config_flags%num_soil_layers-1+1
ENDIF
IF ( SIZE(grid%stc) .GT. 1 ) THEN 
DO k = 1,config_flags%num_soil_layers
xv(k+ixv+1-1)= grid%stc(pig,k,pjg)
ENDDO
ixv=ixv+config_flags%num_soil_layers-1+1
ENDIF
IF ( SIZE(grid%ctopo) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%ctopo(pig,pjg)
ENDIF
IF ( SIZE(grid%ctopo2) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%ctopo2(pig,pjg)
ENDIF
IF ( SIZE(grid%dwdt) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
xv(k+ixv+1-ckds)= grid%dwdt(pig,pjg,k)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%pint) .GT. 1 ) THEN 
DO k = ckds,ckde
xv(k+ixv+1-ckds)= grid%pint(pig,pjg,k)
ENDDO
ixv=ixv+ckde-ckds+1
ENDIF
IF ( SIZE(grid%w) .GT. 1 ) THEN 
DO k = ckds,ckde
xv(k+ixv+1-ckds)= grid%w(pig,pjg,k)
ENDDO
ixv=ixv+ckde-ckds+1
ENDIF
IF ( SIZE(grid%w_tot) .GT. 1 ) THEN 
DO k = ckds,ckde
xv(k+ixv+1-ckds)= grid%w_tot(pig,pjg,k)
ENDDO
ixv=ixv+ckde-ckds+1
ENDIF
IF ( SIZE(grid%z) .GT. 1 ) THEN 
DO k = ckds,ckde
xv(k+ixv+1-ckds)= grid%z(pig,pjg,k)
ENDDO
ixv=ixv+ckde-ckds+1
ENDIF
IF ( SIZE(grid%rlwin) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%rlwin(pig,pjg)
ENDIF
IF ( SIZE(grid%rswin) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%rswin(pig,pjg)
ENDIF
IF ( SIZE(grid%rlwtt) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
xv(k+ixv+1-ckds)= grid%rlwtt(pig,pjg,k)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%rswtt) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
xv(k+ixv+1-ckds)= grid%rswtt(pig,pjg,k)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%winfo) .GT. 1 ) THEN 
DO k = ckds,ckde
xv(k+ixv+1-ckds)= grid%winfo(pig,pjg,k)
ENDDO
ixv=ixv+ckde-ckds+1
ENDIF
IF ( SIZE(grid%iinfo) .GT. 1 ) THEN 
DO k = ckds,ckde
xv(k+ixv+1-ckds)= grid%iinfo(pig,pjg,k)
ENDDO
ixv=ixv+ckde-ckds+1
ENDIF
IF ( SIZE(grid%landmask) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%landmask(pig,pjg)
ENDIF
IF ( SIZE(grid%toposoil) .GT. 1 ) THEN 
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
IF ( SIZE(grid%lake_depth) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%lake_depth(pig,pjg)
ENDIF
IF ( SIZE(grid%u10) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%u10(pig,pjg)
ENDIF
IF ( SIZE(grid%v10) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%v10(pig,pjg)
ENDIF
IF ( SIZE(grid%xice) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%xice(pig,pjg)
ENDIF
IF ( SIZE(grid%icedepth) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%icedepth(pig,pjg)
ENDIF
IF ( SIZE(grid%albsi) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%albsi(pig,pjg)
ENDIF
IF ( SIZE(grid%snowsi) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%snowsi(pig,pjg)
ENDIF
IF ( SIZE(grid%ivgtyp) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%ivgtyp(pig,pjg)
ENDIF
IF ( SIZE(grid%isltyp) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%isltyp(pig,pjg)
ENDIF
IF ( SIZE(grid%vegfra) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%vegfra(pig,pjg)
ENDIF
IF ( SIZE(grid%sst) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%sst(pig,pjg)
ENDIF
IF ( SIZE(grid%weasd) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%weasd(pig,pjg)
ENDIF
IF ( SIZE(grid%thz0) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%thz0(pig,pjg)
ENDIF
IF ( SIZE(grid%qz0) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%qz0(pig,pjg)
ENDIF
IF ( SIZE(grid%uz0) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%uz0(pig,pjg)
ENDIF
IF ( SIZE(grid%vz0) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%vz0(pig,pjg)
ENDIF
IF ( SIZE(grid%htop) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%htop(pig,pjg)
ENDIF
IF ( SIZE(grid%hbot) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%hbot(pig,pjg)
ENDIF
IF ( SIZE(grid%htopr) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%htopr(pig,pjg)
ENDIF
IF ( SIZE(grid%hbotr) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%hbotr(pig,pjg)
ENDIF
IF ( SIZE(grid%cuppt) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%cuppt(pig,pjg)
ENDIF
IF ( SIZE(grid%snowh) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%snowh(pig,pjg)
ENDIF
IF ( SIZE(grid%rhosn) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%rhosn(pig,pjg)
ENDIF
IF ( SIZE(grid%isnowxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%isnowxy(pig,pjg)
ENDIF
IF ( SIZE(grid%tvxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%tvxy(pig,pjg)
ENDIF
IF ( SIZE(grid%tgxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%tgxy(pig,pjg)
ENDIF
IF ( SIZE(grid%canicexy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%canicexy(pig,pjg)
ENDIF
IF ( SIZE(grid%canliqxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%canliqxy(pig,pjg)
ENDIF
IF ( SIZE(grid%eahxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%eahxy(pig,pjg)
ENDIF
IF ( SIZE(grid%tahxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%tahxy(pig,pjg)
ENDIF
IF ( SIZE(grid%cmxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%cmxy(pig,pjg)
ENDIF
IF ( SIZE(grid%chxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%chxy(pig,pjg)
ENDIF
IF ( SIZE(grid%fwetxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%fwetxy(pig,pjg)
ENDIF
IF ( SIZE(grid%sneqvoxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%sneqvoxy(pig,pjg)
ENDIF
IF ( SIZE(grid%alboldxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%alboldxy(pig,pjg)
ENDIF
IF ( SIZE(grid%qsnowxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%qsnowxy(pig,pjg)
ENDIF
IF ( SIZE(grid%wslakexy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%wslakexy(pig,pjg)
ENDIF
IF ( SIZE(grid%zwtxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%zwtxy(pig,pjg)
ENDIF
IF ( SIZE(grid%waxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%waxy(pig,pjg)
ENDIF
IF ( SIZE(grid%wtxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%wtxy(pig,pjg)
ENDIF
IF ( SIZE(grid%tsnoxy) .GT. 1 ) THEN 
DO k = 1,config_flags%num_snow_layers
xv(k+ixv+1-1)= grid%tsnoxy(pig,k,pjg)
ENDDO
ixv=ixv+config_flags%num_snow_layers-1+1
ENDIF
IF ( SIZE(grid%zsnsoxy) .GT. 1 ) THEN 
DO k = 1,config_flags%num_snso_layers
xv(k+ixv+1-1)= grid%zsnsoxy(pig,k,pjg)
ENDDO
ixv=ixv+config_flags%num_snso_layers-1+1
ENDIF
IF ( SIZE(grid%snicexy) .GT. 1 ) THEN 
DO k = 1,config_flags%num_snow_layers
xv(k+ixv+1-1)= grid%snicexy(pig,k,pjg)
ENDDO
ixv=ixv+config_flags%num_snow_layers-1+1
ENDIF
IF ( SIZE(grid%snliqxy) .GT. 1 ) THEN 
DO k = 1,config_flags%num_snow_layers
xv(k+ixv+1-1)= grid%snliqxy(pig,k,pjg)
ENDDO
ixv=ixv+config_flags%num_snow_layers-1+1
ENDIF
IF ( SIZE(grid%lfmassxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%lfmassxy(pig,pjg)
ENDIF
IF ( SIZE(grid%rtmassxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%rtmassxy(pig,pjg)
ENDIF
IF ( SIZE(grid%stmassxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%stmassxy(pig,pjg)
ENDIF
IF ( SIZE(grid%woodxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%woodxy(pig,pjg)
ENDIF
IF ( SIZE(grid%stblcpxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%stblcpxy(pig,pjg)
ENDIF
IF ( SIZE(grid%fastcpxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%fastcpxy(pig,pjg)
ENDIF
IF ( SIZE(grid%xsaixy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%xsaixy(pig,pjg)
ENDIF
IF ( SIZE(grid%t2mvxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%t2mvxy(pig,pjg)
ENDIF
IF ( SIZE(grid%t2mbxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%t2mbxy(pig,pjg)
ENDIF
IF ( SIZE(grid%q2mvxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%q2mvxy(pig,pjg)
ENDIF
IF ( SIZE(grid%q2mbxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%q2mbxy(pig,pjg)
ENDIF
IF ( SIZE(grid%tradxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%tradxy(pig,pjg)
ENDIF
IF ( SIZE(grid%neexy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%neexy(pig,pjg)
ENDIF
IF ( SIZE(grid%gppxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%gppxy(pig,pjg)
ENDIF
IF ( SIZE(grid%nppxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%nppxy(pig,pjg)
ENDIF
IF ( SIZE(grid%fvegxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%fvegxy(pig,pjg)
ENDIF
IF ( SIZE(grid%qinxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%qinxy(pig,pjg)
ENDIF
IF ( SIZE(grid%runsfxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%runsfxy(pig,pjg)
ENDIF
IF ( SIZE(grid%runsbxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%runsbxy(pig,pjg)
ENDIF
IF ( SIZE(grid%ecanxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%ecanxy(pig,pjg)
ENDIF
IF ( SIZE(grid%edirxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%edirxy(pig,pjg)
ENDIF
IF ( SIZE(grid%etranxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%etranxy(pig,pjg)
ENDIF
IF ( SIZE(grid%fsaxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%fsaxy(pig,pjg)
ENDIF
IF ( SIZE(grid%firaxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%firaxy(pig,pjg)
ENDIF
IF ( SIZE(grid%aparxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%aparxy(pig,pjg)
ENDIF
IF ( SIZE(grid%psnxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%psnxy(pig,pjg)
ENDIF
IF ( SIZE(grid%savxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%savxy(pig,pjg)
ENDIF
IF ( SIZE(grid%sagxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%sagxy(pig,pjg)
ENDIF
IF ( SIZE(grid%rssunxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%rssunxy(pig,pjg)
ENDIF
IF ( SIZE(grid%rsshaxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%rsshaxy(pig,pjg)
ENDIF
IF ( SIZE(grid%bgapxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%bgapxy(pig,pjg)
ENDIF
IF ( SIZE(grid%wgapxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%wgapxy(pig,pjg)
ENDIF
IF ( SIZE(grid%tgvxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%tgvxy(pig,pjg)
ENDIF
IF ( SIZE(grid%tgbxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%tgbxy(pig,pjg)
ENDIF
IF ( SIZE(grid%chvxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%chvxy(pig,pjg)
ENDIF
IF ( SIZE(grid%chbxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%chbxy(pig,pjg)
ENDIF
IF ( SIZE(grid%shgxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%shgxy(pig,pjg)
ENDIF
IF ( SIZE(grid%shcxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%shcxy(pig,pjg)
ENDIF
IF ( SIZE(grid%shbxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%shbxy(pig,pjg)
ENDIF
IF ( SIZE(grid%evgxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%evgxy(pig,pjg)
ENDIF
IF ( SIZE(grid%evbxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%evbxy(pig,pjg)
ENDIF
IF ( SIZE(grid%ghvxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%ghvxy(pig,pjg)
ENDIF
IF ( SIZE(grid%ghbxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%ghbxy(pig,pjg)
ENDIF
IF ( SIZE(grid%irgxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%irgxy(pig,pjg)
ENDIF
IF ( SIZE(grid%ircxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%ircxy(pig,pjg)
ENDIF
IF ( SIZE(grid%irbxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%irbxy(pig,pjg)
ENDIF
IF ( SIZE(grid%trxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%trxy(pig,pjg)
ENDIF
IF ( SIZE(grid%evcxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%evcxy(pig,pjg)
ENDIF
IF ( SIZE(grid%chleafxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%chleafxy(pig,pjg)
ENDIF
IF ( SIZE(grid%chucxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%chucxy(pig,pjg)
ENDIF
IF ( SIZE(grid%chv2xy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%chv2xy(pig,pjg)
ENDIF
IF ( SIZE(grid%chb2xy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%chb2xy(pig,pjg)
ENDIF
IF ( SIZE(grid%chstarxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%chstarxy(pig,pjg)
ENDIF
IF ( SIZE(grid%soilcomp) .GT. 1 ) THEN 
DO k = 1,8
xv(k+ixv+1-1)= grid%soilcomp(pig,k,pjg)
ENDDO
ixv=ixv+8-1+1
ENDIF
IF ( SIZE(grid%soilcl1) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%soilcl1(pig,pjg)
ENDIF
IF ( SIZE(grid%soilcl2) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%soilcl2(pig,pjg)
ENDIF
IF ( SIZE(grid%soilcl3) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%soilcl3(pig,pjg)
ENDIF
IF ( SIZE(grid%soilcl4) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%soilcl4(pig,pjg)
ENDIF
IF ( SIZE(grid%grainxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%grainxy(pig,pjg)
ENDIF
IF ( SIZE(grid%gddxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%gddxy(pig,pjg)
ENDIF
IF ( SIZE(grid%croptype) .GT. 1 ) THEN 
DO k = 1,5
xv(k+ixv+1-1)= grid%croptype(pig,k,pjg)
ENDDO
ixv=ixv+5-1+1
ENDIF
IF ( SIZE(grid%planting) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%planting(pig,pjg)
ENDIF
IF ( SIZE(grid%harvest) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%harvest(pig,pjg)
ENDIF
IF ( SIZE(grid%season_gdd) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%season_gdd(pig,pjg)
ENDIF
IF ( SIZE(grid%gecros_state) .GT. 1 ) THEN 
DO k = 1,60
xv(k+ixv+1-1)= grid%gecros_state(pig,k,pjg)
ENDDO
ixv=ixv+60-1+1
ENDIF
CALL rsl_lite_to_parent_msg(ixv*4,xv)
CALL rsl_lite_to_child_info( msize*4                               &
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

      END IF

      
      IF ( intercomm_active( grid%id ) ) THEN        
        local_comm = mpi_comm_to_kid( which_kid(ngrid%id), grid%id )
        ioffset = nest_task_offsets(ngrid%id)
      ELSE IF ( intercomm_active( ngrid%id ) ) THEN  
        local_comm = mpi_comm_to_mom( ngrid%id )
        ioffset = nest_task_offsets(ngrid%id)
      END IF

      IF ( grid%active_this_task .OR. ngrid%active_this_task ) THEN
        CALL mpi_comm_rank(local_comm,myproc,ierr)
        CALL mpi_comm_size(local_comm,nproc,ierr)

        CALL rsl_lite_bcast_msgs( myproc, nest_pes_x(grid%id)*nest_pes_y(grid%id),         &
                                        nest_pes_x(ngrid%id)*nest_pes_y(ngrid%id),       &
                                        ioffset, local_comm )
      END IF

      RETURN
   END SUBROUTINE interp_domain_nmm_part1



   SUBROUTINE interp_domain_nmm_part2 ( grid, ngrid, config_flags    &







,szj,s1z,spz,tcs,moist,moist_bxs,moist_bxe,moist_bys,moist_bye,moist_btxs,moist_btxe,moist_btys,moist_btye,dfi_moist, &
dfi_moist_bxs,dfi_moist_bxe,dfi_moist_bys,dfi_moist_bye,dfi_moist_btxs,dfi_moist_btxe,dfi_moist_btys,dfi_moist_btye,scalar, &
scalar_bxs,scalar_bxe,scalar_bys,scalar_bye,scalar_btxs,scalar_btxe,scalar_btys,scalar_btye,dfi_scalar,dfi_scalar_bxs, &
dfi_scalar_bxe,dfi_scalar_bys,dfi_scalar_bye,dfi_scalar_btxs,dfi_scalar_btxe,dfi_scalar_btys,dfi_scalar_btye,chem,ozmixm &


                 )
      USE module_state_description
      USE module_domain, ONLY : domain, get_ijk_from_grid
      USE module_configure, ONLY : grid_config_rec_type
      USE module_dm, ONLY : ntasks, ntasks_x, ntasks_y, itrace, local_communicator, mytask, &
                            ipe_save, jpe_save, ips_save, jps_save, get_dm_max_halo_width ,  &
                            nest_task_offsets
                             
      USE module_comm_nesting_dm, ONLY : halo_interp_down_sub
      IMPLICIT NONE

      TYPE(domain), POINTER :: grid          
      TYPE(domain), POINTER :: ngrid






real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_szj)           :: szj
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_s1z)           :: s1z
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_spz)           :: spz
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_tcs)           :: tcs
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_moist)           :: moist
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_bxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_bxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_bys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_bye
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_btxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_btxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_btys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_btye
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_dfi_moist)           :: dfi_moist
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_bxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_bxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_bys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_bye
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_btxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_btxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_btys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_btye
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_scalar)           :: scalar
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_bxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_bxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_bys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_bye
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_btxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_btxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_btys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_btye
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_dfi_scalar)           :: dfi_scalar
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_bxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_bxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_bys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_bye
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_btxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_btxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_btys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_btye
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%sm32:grid%em32,num_chem)           :: chem
real      ,DIMENSION(grid%sm31:grid%em31,1:grid%levsiz,grid%sm32:grid%em32,num_ozmixm)           :: ozmixm

      INTEGER nlev, msize
      INTEGER i,j,pig,pjg,cm,cn,nig,njg,retval,k
      TYPE (grid_config_rec_type)            :: config_flags
      REAL xv(2000000)
      INTEGER ixv
      INTEGER       ::          cids, cide, cjds, cjde, ckds, ckde,    &
                                cims, cime, cjms, cjme, ckms, ckme,    &
                                cips, cipe, cjps, cjpe, ckps, ckpe
      INTEGER       ::          nids, nide, njds, njde, nkds, nkde,    &
                                nims, nime, njms, njme, nkms, nkme,    &
                                nips, nipe, njps, njpe, nkps, nkpe
      INTEGER       ::          ids, ide, jds, jde, kds, kde,    &
                                ims, ime, jms, jme, kms, kme,    &
                                ips, ipe, jps, jpe, kps, kpe

      INTEGER idim1,idim2,idim3,idim4,idim5,idim6,idim7
      LOGICAL feedback_flag, feedback_flag_v
      INTEGER myproc
      INTEGER ierr

      integer, parameter :: EConst=0, ECopy=1, EExtrap=2 
      LOGICAL interp_mp




      interp_mp=grid%interp_mp .or. ngrid%interp_mp

      IF ( ngrid%active_this_task ) THEN

        CALL get_ijk_from_grid (  grid ,                   &
                                  cids, cide, cjds, cjde, ckds, ckde,    &
                                  cims, cime, cjms, cjme, ckms, ckme,    &
                                  cips, cipe, cjps, cjpe, ckps, ckpe    )
        CALL get_ijk_from_grid (  ngrid ,              &
                                  nids, nide, njds, njde, nkds, nkde,    &
                                  nims, nime, njms, njme, nkms, nkme,    &
                                  nips, nipe, njps, njpe, nkps, nkpe    )


        nlev  = ckde - ckds + 1







CALL rsl_lite_from_parent_info(pig,pjg,retval)
DO while ( retval .eq. 1 )
CALL rsl_lite_from_parent_msg(-1,xv) 
ixv=0 
IF ( SIZE(grid%lakedepth2d) .GT. 1 ) THEN 
ixv=ixv+1
grid%lakedepth2d(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%savedtke12d) .GT. 1 ) THEN 
ixv=ixv+1
grid%savedtke12d(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%snowdp2d) .GT. 1 ) THEN 
ixv=ixv+1
grid%snowdp2d(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%h2osno2d) .GT. 1 ) THEN 
ixv=ixv+1
grid%h2osno2d(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%snl2d) .GT. 1 ) THEN 
ixv=ixv+1
grid%snl2d(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%t_grnd2d) .GT. 1 ) THEN 
ixv=ixv+1
grid%t_grnd2d(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%t_lake3d) .GT. 1 ) THEN 
DO k = 1,10
grid%t_lake3d(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%lake_icefrac3d) .GT. 1 ) THEN 
DO k = 1,10
grid%lake_icefrac3d(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%z_lake3d) .GT. 1 ) THEN 
DO k = 1,10
grid%z_lake3d(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%dz_lake3d) .GT. 1 ) THEN 
DO k = 1,10
grid%dz_lake3d(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%t_soisno3d) .GT. 1 ) THEN 
DO k = 1,15
grid%t_soisno3d(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+15-1+1
ENDIF
IF ( SIZE(grid%h2osoi_ice3d) .GT. 1 ) THEN 
DO k = 1,15
grid%h2osoi_ice3d(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+15-1+1
ENDIF
IF ( SIZE(grid%h2osoi_liq3d) .GT. 1 ) THEN 
DO k = 1,15
grid%h2osoi_liq3d(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+15-1+1
ENDIF
IF ( SIZE(grid%h2osoi_vol3d) .GT. 1 ) THEN 
DO k = 1,15
grid%h2osoi_vol3d(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+15-1+1
ENDIF
IF ( SIZE(grid%z3d) .GT. 1 ) THEN 
DO k = 1,15
grid%z3d(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+15-1+1
ENDIF
IF ( SIZE(grid%dz3d) .GT. 1 ) THEN 
DO k = 1,15
grid%dz3d(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+15-1+1
ENDIF
IF ( SIZE(grid%zi3d) .GT. 1 ) THEN 
DO k = 1,16
grid%zi3d(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+16-1+1
ENDIF
IF ( SIZE(grid%watsat3d) .GT. 1 ) THEN 
DO k = 1,10
grid%watsat3d(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%csol3d) .GT. 1 ) THEN 
DO k = 1,10
grid%csol3d(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%tkmg3d) .GT. 1 ) THEN 
DO k = 1,10
grid%tkmg3d(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%tkdry3d) .GT. 1 ) THEN 
DO k = 1,10
grid%tkdry3d(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%tksatu3d) .GT. 1 ) THEN 
DO k = 1,10
grid%tksatu3d(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%vertstrucc) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
grid%vertstrucc(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%vertstrucs) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
grid%vertstrucs(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%ru_tendf_stoch) .GT. 1 ) THEN 
DO k = 1,config_flags%num_stoch_levels
grid%ru_tendf_stoch(pig,pjg,k) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+config_flags%num_stoch_levels-1+1
ENDIF
IF ( SIZE(grid%rv_tendf_stoch) .GT. 1 ) THEN 
DO k = 1,config_flags%num_stoch_levels
grid%rv_tendf_stoch(pig,pjg,k) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+config_flags%num_stoch_levels-1+1
ENDIF
IF ( SIZE(grid%lu_index) .GT. 1 ) THEN 
ixv=ixv+1
grid%lu_index(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%lu_mask) .GT. 1 ) THEN 
ixv=ixv+1
grid%lu_mask(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%cuprecip_swath) .GT. 1 ) THEN 
ixv=ixv+1
grid%cuprecip_swath(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%precip_swath) .GT. 1 ) THEN 
ixv=ixv+1
grid%precip_swath(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%windsq_swath) .GT. 1 ) THEN 
ixv=ixv+1
grid%windsq_swath(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%suswind) .GT. 1 ) THEN 
ixv=ixv+1
grid%suswind(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%wind10_ratio) .GT. 1 ) THEN 
ixv=ixv+1
grid%wind10_ratio(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%suswind_swath) .GT. 1 ) THEN 
ixv=ixv+1
grid%suswind_swath(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%membrane_mslp) .GT. 1 ) THEN 
ixv=ixv+1
grid%membrane_mslp(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%pdyn_smooth) .GT. 1 ) THEN 
ixv=ixv+1
grid%pdyn_smooth(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%pdyn_parent) .GT. 1 ) THEN 
ixv=ixv+1
grid%pdyn_parent(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%best_mslp) .GT. 1 ) THEN 
ixv=ixv+1
grid%best_mslp(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%hlon) .GT. 1 ) THEN 
ixv=ixv+1
grid%hlon(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%hlat) .GT. 1 ) THEN 
ixv=ixv+1
grid%hlat(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tg_max_m10wind) .GT. 1 ) THEN 
ixv=ixv+1
grid%tg_max_m10wind(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tg_max_wwind) .GT. 1 ) THEN 
ixv=ixv+1
grid%tg_max_wwind(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tg_min_wwind) .GT. 1 ) THEN 
ixv=ixv+1
grid%tg_min_wwind(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tg_max_zhel_25) .GT. 1 ) THEN 
ixv=ixv+1
grid%tg_max_zhel_25(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tg_min_zhel_25) .GT. 1 ) THEN 
ixv=ixv+1
grid%tg_min_zhel_25(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tg_max_zhel_03) .GT. 1 ) THEN 
ixv=ixv+1
grid%tg_max_zhel_03(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tg_min_zhel_03) .GT. 1 ) THEN 
ixv=ixv+1
grid%tg_min_zhel_03(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tg_updhel25) .GT. 1 ) THEN 
ixv=ixv+1
grid%tg_updhel25(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tg_max_updhel25) .GT. 1 ) THEN 
ixv=ixv+1
grid%tg_max_updhel25(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tg_updhel03) .GT. 1 ) THEN 
ixv=ixv+1
grid%tg_updhel03(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tg_max_updhel03) .GT. 1 ) THEN 
ixv=ixv+1
grid%tg_max_updhel03(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tg_total_precip) .GT. 1 ) THEN 
ixv=ixv+1
grid%tg_total_precip(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%hres_fis) .GT. 1 ) THEN 
ixv=ixv+1
grid%hres_fis(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%hres_lnd) .GT. 1 ) THEN 
ixv=ixv+1
grid%hres_lnd(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%hres_lu_index) .GT. 1 ) THEN 
ixv=ixv+1
grid%hres_lu_index(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%sm) .GT. 1 ) THEN 
ixv=ixv+1
grid%sm(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%sice) .GT. 1 ) THEN 
ixv=ixv+1
grid%sice(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%pd) .GT. 1 ) THEN 
ixv=ixv+1
grid%pd(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%fis) .GT. 1 ) THEN 
ixv=ixv+1
grid%fis(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%t) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
grid%t(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%q) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
grid%q(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%u) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
grid%u(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%v) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
grid%v(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%ustar) .GT. 1 ) THEN 
ixv=ixv+1
grid%ustar(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%z0) .GT. 1 ) THEN 
ixv=ixv+1
grid%z0(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%scurx) .GT. 1 ) THEN 
ixv=ixv+1
grid%scurx(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%scury) .GT. 1 ) THEN 
ixv=ixv+1
grid%scury(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%charn) .GT. 1 ) THEN 
ixv=ixv+1
grid%charn(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%msang) .GT. 1 ) THEN 
ixv=ixv+1
grid%msang(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%ths) .GT. 1 ) THEN 
ixv=ixv+1
grid%ths(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%qsh) .GT. 1 ) THEN 
ixv=ixv+1
grid%qsh(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%taux) .GT. 1 ) THEN 
ixv=ixv+1
grid%taux(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tauy) .GT. 1 ) THEN 
ixv=ixv+1
grid%tauy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%acprec) .GT. 1 ) THEN 
ixv=ixv+1
grid%acprec(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%cuprec) .GT. 1 ) THEN 
ixv=ixv+1
grid%cuprec(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%cldefi) .GT. 1 ) THEN 
ixv=ixv+1
grid%cldefi(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%th10) .GT. 1 ) THEN 
ixv=ixv+1
grid%th10(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%q10) .GT. 1 ) THEN 
ixv=ixv+1
grid%q10(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%pshltr) .GT. 1 ) THEN 
ixv=ixv+1
grid%pshltr(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tshltr) .GT. 1 ) THEN 
ixv=ixv+1
grid%tshltr(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%qshltr) .GT. 1 ) THEN 
ixv=ixv+1
grid%qshltr(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%q2) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
grid%q2(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%t_adj) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
grid%t_adj(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%cd_out) .GT. 1 ) THEN 
ixv=ixv+1
grid%cd_out(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%ch_out) .GT. 1 ) THEN 
ixv=ixv+1
grid%ch_out(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%albase) .GT. 1 ) THEN 
ixv=ixv+1
grid%albase(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%nmm_tsk) .GT. 1 ) THEN 
ixv=ixv+1
grid%nmm_tsk(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%mxsnal) .GT. 1 ) THEN 
ixv=ixv+1
grid%mxsnal(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%sigt4) .GT. 1 ) THEN 
ixv=ixv+1
grid%sigt4(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tg) .GT. 1 ) THEN 
ixv=ixv+1
grid%tg(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%acswupt) .GT. 1 ) THEN 
ixv=ixv+1
grid%acswupt(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%acswuptc) .GT. 1 ) THEN 
ixv=ixv+1
grid%acswuptc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%acswdnt) .GT. 1 ) THEN 
ixv=ixv+1
grid%acswdnt(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%acswdntc) .GT. 1 ) THEN 
ixv=ixv+1
grid%acswdntc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%acswupb) .GT. 1 ) THEN 
ixv=ixv+1
grid%acswupb(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%acswupbc) .GT. 1 ) THEN 
ixv=ixv+1
grid%acswupbc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%acswdnb) .GT. 1 ) THEN 
ixv=ixv+1
grid%acswdnb(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%acswdnbc) .GT. 1 ) THEN 
ixv=ixv+1
grid%acswdnbc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%aclwupt) .GT. 1 ) THEN 
ixv=ixv+1
grid%aclwupt(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%aclwuptc) .GT. 1 ) THEN 
ixv=ixv+1
grid%aclwuptc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%aclwdnt) .GT. 1 ) THEN 
ixv=ixv+1
grid%aclwdnt(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%aclwdntc) .GT. 1 ) THEN 
ixv=ixv+1
grid%aclwdntc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%aclwupb) .GT. 1 ) THEN 
ixv=ixv+1
grid%aclwupb(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%aclwupbc) .GT. 1 ) THEN 
ixv=ixv+1
grid%aclwupbc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%aclwdnb) .GT. 1 ) THEN 
ixv=ixv+1
grid%aclwdnb(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%aclwdnbc) .GT. 1 ) THEN 
ixv=ixv+1
grid%aclwdnbc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%swupt) .GT. 1 ) THEN 
ixv=ixv+1
grid%swupt(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%swuptc) .GT. 1 ) THEN 
ixv=ixv+1
grid%swuptc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%swdnt) .GT. 1 ) THEN 
ixv=ixv+1
grid%swdnt(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%swdntc) .GT. 1 ) THEN 
ixv=ixv+1
grid%swdntc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%swupb) .GT. 1 ) THEN 
ixv=ixv+1
grid%swupb(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%swupbc) .GT. 1 ) THEN 
ixv=ixv+1
grid%swupbc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%swdnb) .GT. 1 ) THEN 
ixv=ixv+1
grid%swdnb(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%swdnbc) .GT. 1 ) THEN 
ixv=ixv+1
grid%swdnbc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%lwupt) .GT. 1 ) THEN 
ixv=ixv+1
grid%lwupt(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%lwuptc) .GT. 1 ) THEN 
ixv=ixv+1
grid%lwuptc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%lwdnt) .GT. 1 ) THEN 
ixv=ixv+1
grid%lwdnt(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%lwdntc) .GT. 1 ) THEN 
ixv=ixv+1
grid%lwdntc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%lwupb) .GT. 1 ) THEN 
ixv=ixv+1
grid%lwupb(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%lwupbc) .GT. 1 ) THEN 
ixv=ixv+1
grid%lwupbc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%lwdnb) .GT. 1 ) THEN 
ixv=ixv+1
grid%lwdnb(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%lwdnbc) .GT. 1 ) THEN 
ixv=ixv+1
grid%lwdnbc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%refl_10cm) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
grid%refl_10cm(pig,k,pjg) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%refd_max) .GT. 1 ) THEN 
ixv=ixv+1
grid%refd_max(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%qnwfa2d) .GT. 1 ) THEN 
ixv=ixv+1
grid%qnwfa2d(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%swddir) .GT. 1 ) THEN 
ixv=ixv+1
grid%swddir(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%swddirc) .GT. 1 ) THEN 
ixv=ixv+1
grid%swddirc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%swddni) .GT. 1 ) THEN 
ixv=ixv+1
grid%swddni(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%swddnic) .GT. 1 ) THEN 
ixv=ixv+1
grid%swddnic(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%swddif) .GT. 1 ) THEN 
ixv=ixv+1
grid%swddif(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%gx) .GT. 1 ) THEN 
ixv=ixv+1
grid%gx(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%bx) .GT. 1 ) THEN 
ixv=ixv+1
grid%bx(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%gg) .GT. 1 ) THEN 
ixv=ixv+1
grid%gg(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%bb) .GT. 1 ) THEN 
ixv=ixv+1
grid%bb(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%coszen_ref) .GT. 1 ) THEN 
ixv=ixv+1
grid%coszen_ref(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%swdown_ref) .GT. 1 ) THEN 
ixv=ixv+1
grid%swdown_ref(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%swddir_ref) .GT. 1 ) THEN 
ixv=ixv+1
grid%swddir_ref(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%cwm) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
grid%cwm(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%f_ice) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
grid%f_ice(pig,k,pjg) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%f_rain) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
grid%f_rain(pig,k,pjg) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%f_rimef) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
grid%f_rimef(pig,k,pjg) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%cfrach) .GT. 1 ) THEN 
ixv=ixv+1
grid%cfrach(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%cfracl) .GT. 1 ) THEN 
ixv=ixv+1
grid%cfracl(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%cfracm) .GT. 1 ) THEN 
ixv=ixv+1
grid%cfracm(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%islope) .GT. 1 ) THEN 
ixv=ixv+1
grid%islope(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%cmc) .GT. 1 ) THEN 
ixv=ixv+1
grid%cmc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%soiltb) .GT. 1 ) THEN 
ixv=ixv+1
grid%soiltb(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%vegfrc) .GT. 1 ) THEN 
ixv=ixv+1
grid%vegfrc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%shdmax) .GT. 1 ) THEN 
ixv=ixv+1
grid%shdmax(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%shdmin) .GT. 1 ) THEN 
ixv=ixv+1
grid%shdmin(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%sh2o) .GT. 1 ) THEN 
DO k = 1,config_flags%num_soil_layers
grid%sh2o(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+config_flags%num_soil_layers-1+1
ENDIF
IF ( SIZE(grid%smc) .GT. 1 ) THEN 
DO k = 1,config_flags%num_soil_layers
grid%smc(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+config_flags%num_soil_layers-1+1
ENDIF
IF ( SIZE(grid%stc) .GT. 1 ) THEN 
DO k = 1,config_flags%num_soil_layers
grid%stc(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+config_flags%num_soil_layers-1+1
ENDIF
IF ( SIZE(grid%ctopo) .GT. 1 ) THEN 
ixv=ixv+1
grid%ctopo(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%ctopo2) .GT. 1 ) THEN 
ixv=ixv+1
grid%ctopo2(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%dwdt) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
grid%dwdt(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%pint) .GT. 1 ) THEN 
DO k = ckds,ckde
grid%pint(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+ckde-ckds+1
ENDIF
IF ( SIZE(grid%w) .GT. 1 ) THEN 
DO k = ckds,ckde
grid%w(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+ckde-ckds+1
ENDIF
IF ( SIZE(grid%w_tot) .GT. 1 ) THEN 
DO k = ckds,ckde
grid%w_tot(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+ckde-ckds+1
ENDIF
IF ( SIZE(grid%z) .GT. 1 ) THEN 
DO k = ckds,ckde
grid%z(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+ckde-ckds+1
ENDIF
IF ( SIZE(grid%rlwin) .GT. 1 ) THEN 
ixv=ixv+1
grid%rlwin(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%rswin) .GT. 1 ) THEN 
ixv=ixv+1
grid%rswin(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%rlwtt) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
grid%rlwtt(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%rswtt) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
grid%rswtt(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%winfo) .GT. 1 ) THEN 
DO k = ckds,ckde
grid%winfo(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+ckde-ckds+1
ENDIF
IF ( SIZE(grid%iinfo) .GT. 1 ) THEN 
DO k = ckds,ckde
grid%iinfo(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+ckde-ckds+1
ENDIF
IF ( SIZE(grid%landmask) .GT. 1 ) THEN 
ixv=ixv+1
grid%landmask(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%toposoil) .GT. 1 ) THEN 
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
IF ( SIZE(grid%lake_depth) .GT. 1 ) THEN 
ixv=ixv+1
grid%lake_depth(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%u10) .GT. 1 ) THEN 
ixv=ixv+1
grid%u10(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%v10) .GT. 1 ) THEN 
ixv=ixv+1
grid%v10(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%xice) .GT. 1 ) THEN 
ixv=ixv+1
grid%xice(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%icedepth) .GT. 1 ) THEN 
ixv=ixv+1
grid%icedepth(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%albsi) .GT. 1 ) THEN 
ixv=ixv+1
grid%albsi(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%snowsi) .GT. 1 ) THEN 
ixv=ixv+1
grid%snowsi(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%ivgtyp) .GT. 1 ) THEN 
ixv=ixv+1
grid%ivgtyp(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%isltyp) .GT. 1 ) THEN 
ixv=ixv+1
grid%isltyp(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%vegfra) .GT. 1 ) THEN 
ixv=ixv+1
grid%vegfra(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%sst) .GT. 1 ) THEN 
ixv=ixv+1
grid%sst(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%weasd) .GT. 1 ) THEN 
ixv=ixv+1
grid%weasd(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%thz0) .GT. 1 ) THEN 
ixv=ixv+1
grid%thz0(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%qz0) .GT. 1 ) THEN 
ixv=ixv+1
grid%qz0(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%uz0) .GT. 1 ) THEN 
ixv=ixv+1
grid%uz0(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%vz0) .GT. 1 ) THEN 
ixv=ixv+1
grid%vz0(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%htop) .GT. 1 ) THEN 
ixv=ixv+1
grid%htop(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%hbot) .GT. 1 ) THEN 
ixv=ixv+1
grid%hbot(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%htopr) .GT. 1 ) THEN 
ixv=ixv+1
grid%htopr(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%hbotr) .GT. 1 ) THEN 
ixv=ixv+1
grid%hbotr(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%cuppt) .GT. 1 ) THEN 
ixv=ixv+1
grid%cuppt(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%snowh) .GT. 1 ) THEN 
ixv=ixv+1
grid%snowh(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%rhosn) .GT. 1 ) THEN 
ixv=ixv+1
grid%rhosn(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%isnowxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%isnowxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tvxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%tvxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tgxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%tgxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%canicexy) .GT. 1 ) THEN 
ixv=ixv+1
grid%canicexy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%canliqxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%canliqxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%eahxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%eahxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tahxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%tahxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%cmxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%cmxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%chxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%chxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%fwetxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%fwetxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%sneqvoxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%sneqvoxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%alboldxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%alboldxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%qsnowxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%qsnowxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%wslakexy) .GT. 1 ) THEN 
ixv=ixv+1
grid%wslakexy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%zwtxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%zwtxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%waxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%waxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%wtxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%wtxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tsnoxy) .GT. 1 ) THEN 
DO k = 1,config_flags%num_snow_layers
grid%tsnoxy(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+config_flags%num_snow_layers-1+1
ENDIF
IF ( SIZE(grid%zsnsoxy) .GT. 1 ) THEN 
DO k = 1,config_flags%num_snso_layers
grid%zsnsoxy(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+config_flags%num_snso_layers-1+1
ENDIF
IF ( SIZE(grid%snicexy) .GT. 1 ) THEN 
DO k = 1,config_flags%num_snow_layers
grid%snicexy(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+config_flags%num_snow_layers-1+1
ENDIF
IF ( SIZE(grid%snliqxy) .GT. 1 ) THEN 
DO k = 1,config_flags%num_snow_layers
grid%snliqxy(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+config_flags%num_snow_layers-1+1
ENDIF
IF ( SIZE(grid%lfmassxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%lfmassxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%rtmassxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%rtmassxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%stmassxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%stmassxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%woodxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%woodxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%stblcpxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%stblcpxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%fastcpxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%fastcpxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%xsaixy) .GT. 1 ) THEN 
ixv=ixv+1
grid%xsaixy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%t2mvxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%t2mvxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%t2mbxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%t2mbxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%q2mvxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%q2mvxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%q2mbxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%q2mbxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tradxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%tradxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%neexy) .GT. 1 ) THEN 
ixv=ixv+1
grid%neexy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%gppxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%gppxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%nppxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%nppxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%fvegxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%fvegxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%qinxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%qinxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%runsfxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%runsfxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%runsbxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%runsbxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%ecanxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%ecanxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%edirxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%edirxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%etranxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%etranxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%fsaxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%fsaxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%firaxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%firaxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%aparxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%aparxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%psnxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%psnxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%savxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%savxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%sagxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%sagxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%rssunxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%rssunxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%rsshaxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%rsshaxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%bgapxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%bgapxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%wgapxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%wgapxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tgvxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%tgvxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tgbxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%tgbxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%chvxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%chvxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%chbxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%chbxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%shgxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%shgxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%shcxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%shcxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%shbxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%shbxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%evgxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%evgxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%evbxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%evbxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%ghvxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%ghvxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%ghbxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%ghbxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%irgxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%irgxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%ircxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%ircxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%irbxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%irbxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%trxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%trxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%evcxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%evcxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%chleafxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%chleafxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%chucxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%chucxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%chv2xy) .GT. 1 ) THEN 
ixv=ixv+1
grid%chv2xy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%chb2xy) .GT. 1 ) THEN 
ixv=ixv+1
grid%chb2xy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%chstarxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%chstarxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%soilcomp) .GT. 1 ) THEN 
DO k = 1,8
grid%soilcomp(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+8-1+1
ENDIF
IF ( SIZE(grid%soilcl1) .GT. 1 ) THEN 
ixv=ixv+1
grid%soilcl1(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%soilcl2) .GT. 1 ) THEN 
ixv=ixv+1
grid%soilcl2(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%soilcl3) .GT. 1 ) THEN 
ixv=ixv+1
grid%soilcl3(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%soilcl4) .GT. 1 ) THEN 
ixv=ixv+1
grid%soilcl4(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%grainxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%grainxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%gddxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%gddxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%croptype) .GT. 1 ) THEN 
DO k = 1,5
grid%croptype(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+5-1+1
ENDIF
IF ( SIZE(grid%planting) .GT. 1 ) THEN 
ixv=ixv+1
grid%planting(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%harvest) .GT. 1 ) THEN 
ixv=ixv+1
grid%harvest(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%season_gdd) .GT. 1 ) THEN 
ixv=ixv+1
grid%season_gdd(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%gecros_state) .GT. 1 ) THEN 
DO k = 1,60
grid%gecros_state(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+60-1+1
ENDIF
CALL rsl_lite_from_parent_info(pig,pjg,retval)
ENDDO


        CALL get_ijk_from_grid (  grid ,              &
                                  ids, ide, jds, jde, kds, kde,    &
                                  ims, ime, jms, jme, kms, kme,    &
                                  ips, ipe, jps, jpe, kps, kpe    )

        CALL push_communicators_for_domain( grid%id )






CALL HALO_INTERP_DOWN_sub ( grid, &
  config_flags, &
  num_moist, &
  moist, &
  num_scalar, &
  scalar, &
  num_dfi_scalar, &
  dfi_scalar, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )


      




        call store_interp_info(ngrid,grid)
        call ext_c2n_fulldom(ngrid%IIH,ngrid%JJH,ngrid%HBWGT1, &
             ngrid%HBWGT2,ngrid%HBWGT3,ngrid%HBWGT4,         &
             ngrid%deta1,ngrid%deta2,ngrid%eta1,             &
             ngrid%eta2,ngrid%pt,ngrid%pdtop,                &
             grid%pint,grid%t,grid%pd,grid%q,       &
             cims, cime, cjms, cjme, ckms, ckme,             &
             ngrid%pint,ngrid%t,ngrid%pd,ngrid%q,&
             ngrid%iinfo,ngrid%winfo,ngrid%imask_nostag, &
             nids, nide, njds, njde, nkds, nkde,             &
             nims, nime, njms, njme, nkms, nkme,             &
             nips, nipe, njps, njpe, nkps, nkpe)







IF ( .TRUE. ) THEN 
CALL interp_mask_water_field (  &         
  ( SIZE( grid%lakedepth2d , 1 )*SIZE( grid%lakedepth2d , 2 ) .GT. 1 ), & 
                  grid%lakedepth2d,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%lakedepth2d,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%islake,ngrid%islake&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_water_field (  &         
  ( SIZE( grid%savedtke12d , 1 )*SIZE( grid%savedtke12d , 2 ) .GT. 1 ), & 
                  grid%savedtke12d,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%savedtke12d,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%islake,ngrid%islake&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_water_field (  &         
  ( SIZE( grid%snowdp2d , 1 )*SIZE( grid%snowdp2d , 2 ) .GT. 1 ), & 
                  grid%snowdp2d,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%snowdp2d,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%islake,ngrid%islake&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_water_field (  &         
  ( SIZE( grid%h2osno2d , 1 )*SIZE( grid%h2osno2d , 2 ) .GT. 1 ), & 
                  grid%h2osno2d,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%h2osno2d,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%islake,ngrid%islake&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_water_field (  &         
  ( SIZE( grid%snl2d , 1 )*SIZE( grid%snl2d , 2 ) .GT. 1 ), & 
                  grid%snl2d,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%snl2d,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%islake,ngrid%islake&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_water_field (  &         
  ( SIZE( grid%t_grnd2d , 1 )*SIZE( grid%t_grnd2d , 2 ) .GT. 1 ), & 
                  grid%t_grnd2d,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%t_grnd2d,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%islake,ngrid%islake&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_water_field (  &         
  ( SIZE( grid%t_lake3d , 1 )*SIZE( grid%t_lake3d , 3 ) .GT. 1 ), & 
                  grid%t_lake3d,   &       
                 cids, cide, 1, 10, cjds, cjde,   &         
                 cims, cime, 1, 10, cjms, cjme,   &         
                 cips, cipe, 1, 10, cjps, cjpe,   &         
                  ngrid%t_lake3d,  &   
                 nids, nide, 1, 10, njds, njde,   &         
                 nims, nime, 1, 10, njms, njme,   &         
                 nips, nipe, 1, 10, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%islake,ngrid%islake&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_water_field (  &         
  ( SIZE( grid%lake_icefrac3d , 1 )*SIZE( grid%lake_icefrac3d , 3 ) .GT. 1 ), & 
                  grid%lake_icefrac3d,   &       
                 cids, cide, 1, 10, cjds, cjde,   &         
                 cims, cime, 1, 10, cjms, cjme,   &         
                 cips, cipe, 1, 10, cjps, cjpe,   &         
                  ngrid%lake_icefrac3d,  &   
                 nids, nide, 1, 10, njds, njde,   &         
                 nims, nime, 1, 10, njms, njme,   &         
                 nips, nipe, 1, 10, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%islake,ngrid%islake&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_water_field (  &         
  ( SIZE( grid%z_lake3d , 1 )*SIZE( grid%z_lake3d , 3 ) .GT. 1 ), & 
                  grid%z_lake3d,   &       
                 cids, cide, 1, 10, cjds, cjde,   &         
                 cims, cime, 1, 10, cjms, cjme,   &         
                 cips, cipe, 1, 10, cjps, cjpe,   &         
                  ngrid%z_lake3d,  &   
                 nids, nide, 1, 10, njds, njde,   &         
                 nims, nime, 1, 10, njms, njme,   &         
                 nips, nipe, 1, 10, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%islake,ngrid%islake&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_water_field (  &         
  ( SIZE( grid%dz_lake3d , 1 )*SIZE( grid%dz_lake3d , 3 ) .GT. 1 ), & 
                  grid%dz_lake3d,   &       
                 cids, cide, 1, 10, cjds, cjde,   &         
                 cims, cime, 1, 10, cjms, cjme,   &         
                 cips, cipe, 1, 10, cjps, cjpe,   &         
                  ngrid%dz_lake3d,  &   
                 nids, nide, 1, 10, njds, njde,   &         
                 nims, nime, 1, 10, njms, njme,   &         
                 nips, nipe, 1, 10, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%islake,ngrid%islake&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_water_field (  &         
  ( SIZE( grid%t_soisno3d , 1 )*SIZE( grid%t_soisno3d , 3 ) .GT. 1 ), & 
                  grid%t_soisno3d,   &       
                 cids, cide, 1, 15, cjds, cjde,   &         
                 cims, cime, 1, 15, cjms, cjme,   &         
                 cips, cipe, 1, 15, cjps, cjpe,   &         
                  ngrid%t_soisno3d,  &   
                 nids, nide, 1, 15, njds, njde,   &         
                 nims, nime, 1, 15, njms, njme,   &         
                 nips, nipe, 1, 15, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%islake,ngrid%islake&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_water_field (  &         
  ( SIZE( grid%h2osoi_ice3d , 1 )*SIZE( grid%h2osoi_ice3d , 3 ) .GT. 1 ), & 
                  grid%h2osoi_ice3d,   &       
                 cids, cide, 1, 15, cjds, cjde,   &         
                 cims, cime, 1, 15, cjms, cjme,   &         
                 cips, cipe, 1, 15, cjps, cjpe,   &         
                  ngrid%h2osoi_ice3d,  &   
                 nids, nide, 1, 15, njds, njde,   &         
                 nims, nime, 1, 15, njms, njme,   &         
                 nips, nipe, 1, 15, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%islake,ngrid%islake&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_water_field (  &         
  ( SIZE( grid%h2osoi_liq3d , 1 )*SIZE( grid%h2osoi_liq3d , 3 ) .GT. 1 ), & 
                  grid%h2osoi_liq3d,   &       
                 cids, cide, 1, 15, cjds, cjde,   &         
                 cims, cime, 1, 15, cjms, cjme,   &         
                 cips, cipe, 1, 15, cjps, cjpe,   &         
                  ngrid%h2osoi_liq3d,  &   
                 nids, nide, 1, 15, njds, njde,   &         
                 nims, nime, 1, 15, njms, njme,   &         
                 nips, nipe, 1, 15, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%islake,ngrid%islake&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_water_field (  &         
  ( SIZE( grid%h2osoi_vol3d , 1 )*SIZE( grid%h2osoi_vol3d , 3 ) .GT. 1 ), & 
                  grid%h2osoi_vol3d,   &       
                 cids, cide, 1, 15, cjds, cjde,   &         
                 cims, cime, 1, 15, cjms, cjme,   &         
                 cips, cipe, 1, 15, cjps, cjpe,   &         
                  ngrid%h2osoi_vol3d,  &   
                 nids, nide, 1, 15, njds, njde,   &         
                 nims, nime, 1, 15, njms, njme,   &         
                 nips, nipe, 1, 15, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%islake,ngrid%islake&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_water_field (  &         
  ( SIZE( grid%z3d , 1 )*SIZE( grid%z3d , 3 ) .GT. 1 ), & 
                  grid%z3d,   &       
                 cids, cide, 1, 15, cjds, cjde,   &         
                 cims, cime, 1, 15, cjms, cjme,   &         
                 cips, cipe, 1, 15, cjps, cjpe,   &         
                  ngrid%z3d,  &   
                 nids, nide, 1, 15, njds, njde,   &         
                 nims, nime, 1, 15, njms, njme,   &         
                 nips, nipe, 1, 15, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%islake,ngrid%islake&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_water_field (  &         
  ( SIZE( grid%dz3d , 1 )*SIZE( grid%dz3d , 3 ) .GT. 1 ), & 
                  grid%dz3d,   &       
                 cids, cide, 1, 15, cjds, cjde,   &         
                 cims, cime, 1, 15, cjms, cjme,   &         
                 cips, cipe, 1, 15, cjps, cjpe,   &         
                  ngrid%dz3d,  &   
                 nids, nide, 1, 15, njds, njde,   &         
                 nims, nime, 1, 15, njms, njme,   &         
                 nips, nipe, 1, 15, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%islake,ngrid%islake&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_water_field (  &         
  ( SIZE( grid%zi3d , 1 )*SIZE( grid%zi3d , 3 ) .GT. 1 ), & 
                  grid%zi3d,   &       
                 cids, cide, 1, 16, cjds, cjde,   &         
                 cims, cime, 1, 16, cjms, cjme,   &         
                 cips, cipe, 1, 16, cjps, cjpe,   &         
                  ngrid%zi3d,  &   
                 nids, nide, 1, 16, njds, njde,   &         
                 nims, nime, 1, 16, njms, njme,   &         
                 nips, nipe, 1, 16, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%islake,ngrid%islake&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_water_field (  &         
  ( SIZE( grid%watsat3d , 1 )*SIZE( grid%watsat3d , 3 ) .GT. 1 ), & 
                  grid%watsat3d,   &       
                 cids, cide, 1, 10, cjds, cjde,   &         
                 cims, cime, 1, 10, cjms, cjme,   &         
                 cips, cipe, 1, 10, cjps, cjpe,   &         
                  ngrid%watsat3d,  &   
                 nids, nide, 1, 10, njds, njde,   &         
                 nims, nime, 1, 10, njms, njme,   &         
                 nips, nipe, 1, 10, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%islake,ngrid%islake&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_water_field (  &         
  ( SIZE( grid%csol3d , 1 )*SIZE( grid%csol3d , 3 ) .GT. 1 ), & 
                  grid%csol3d,   &       
                 cids, cide, 1, 10, cjds, cjde,   &         
                 cims, cime, 1, 10, cjms, cjme,   &         
                 cips, cipe, 1, 10, cjps, cjpe,   &         
                  ngrid%csol3d,  &   
                 nids, nide, 1, 10, njds, njde,   &         
                 nims, nime, 1, 10, njms, njme,   &         
                 nips, nipe, 1, 10, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%islake,ngrid%islake&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_water_field (  &         
  ( SIZE( grid%tkmg3d , 1 )*SIZE( grid%tkmg3d , 3 ) .GT. 1 ), & 
                  grid%tkmg3d,   &       
                 cids, cide, 1, 10, cjds, cjde,   &         
                 cims, cime, 1, 10, cjms, cjme,   &         
                 cips, cipe, 1, 10, cjps, cjpe,   &         
                  ngrid%tkmg3d,  &   
                 nids, nide, 1, 10, njds, njde,   &         
                 nims, nime, 1, 10, njms, njme,   &         
                 nips, nipe, 1, 10, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%islake,ngrid%islake&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_water_field (  &         
  ( SIZE( grid%tkdry3d , 1 )*SIZE( grid%tkdry3d , 3 ) .GT. 1 ), & 
                  grid%tkdry3d,   &       
                 cids, cide, 1, 10, cjds, cjde,   &         
                 cims, cime, 1, 10, cjms, cjme,   &         
                 cips, cipe, 1, 10, cjps, cjpe,   &         
                  ngrid%tkdry3d,  &   
                 nids, nide, 1, 10, njds, njde,   &         
                 nims, nime, 1, 10, njms, njme,   &         
                 nips, nipe, 1, 10, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%islake,ngrid%islake&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_water_field (  &         
  ( SIZE( grid%tksatu3d , 1 )*SIZE( grid%tksatu3d , 3 ) .GT. 1 ), & 
                  grid%tksatu3d,   &       
                 cids, cide, 1, 10, cjds, cjde,   &         
                 cims, cime, 1, 10, cjms, cjme,   &         
                 cips, cipe, 1, 10, cjps, cjpe,   &         
                  ngrid%tksatu3d,  &   
                 nids, nide, 1, 10, njds, njde,   &         
                 nims, nime, 1, 10, njms, njme,   &         
                 nips, nipe, 1, 10, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%islake,ngrid%islake&
                  ) 
ENDIF
IF ( SIZE( grid%vertstrucc, 1 ) * SIZE( grid%vertstrucc, 2 ) .GT. 1 ) THEN 
CALL DownCopy (  &         
                  grid%vertstrucc,   &       
                 cids, cide, ckds, ckde, cjds, cjde,   &         
                 cims, cime, ckms, ckme, cjms, cjme,   &         
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         
                  ngrid%vertstrucc,  &   
                 nids, nide, nkds, nkde, njds, njde,   &         
                 nims, nime, nkms, nkme, njms, njme,   &         
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%vertstrucs, 1 ) * SIZE( grid%vertstrucs, 2 ) .GT. 1 ) THEN 
CALL DownCopy (  &         
                  grid%vertstrucs,   &       
                 cids, cide, ckds, ckde, cjds, cjde,   &         
                 cims, cime, ckms, ckme, cjms, cjme,   &         
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         
                  ngrid%vertstrucs,  &   
                 nids, nide, nkds, nkde, njds, njde,   &         
                 nims, nime, nkms, nkme, njms, njme,   &         
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%ru_tendf_stoch, 1 ) * SIZE( grid%ru_tendf_stoch, 2 ) .GT. 1 ) THEN 
CALL downvel (  &         
                  grid%ru_tendf_stoch,   &       
                 cids, cide, 1, config_flags%num_stoch_levels, cjds, cjde,   &         
                 cims, cime, 1, config_flags%num_stoch_levels, cjms, cjme,   &         
                 cips, cipe, 1, config_flags%num_stoch_levels, cjps, cjpe,   &         
                  ngrid%ru_tendf_stoch,  &   
                 nids, nide, 1, config_flags%num_stoch_levels, njds, njde,   &         
                 nims, nime, 1, config_flags%num_stoch_levels, njms, njme,   &         
                 nips, nipe, 1, config_flags%num_stoch_levels, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_xstag,         &         
                  .TRUE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%rv_tendf_stoch, 1 ) * SIZE( grid%rv_tendf_stoch, 2 ) .GT. 1 ) THEN 
CALL downvel (  &         
                  grid%rv_tendf_stoch,   &       
                 cids, cide, 1, config_flags%num_stoch_levels, cjds, cjde,   &         
                 cims, cime, 1, config_flags%num_stoch_levels, cjms, cjme,   &         
                 cips, cipe, 1, config_flags%num_stoch_levels, cjps, cjpe,   &         
                  ngrid%rv_tendf_stoch,  &   
                 nids, nide, 1, config_flags%num_stoch_levels, njds, njde,   &         
                 nims, nime, 1, config_flags%num_stoch_levels, njms, njme,   &         
                 nips, nipe, 1, config_flags%num_stoch_levels, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_ystag,         &         
                  .FALSE., .TRUE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%lu_index, 1 ) * SIZE( grid%lu_index, 2 ) .GT. 1 ) THEN 
CONTINUE 
ENDIF
IF ( SIZE( grid%lu_mask, 1 ) * SIZE( grid%lu_mask, 2 ) .GT. 1 ) THEN 
CONTINUE 
ENDIF
IF ( SIZE( grid%cuprecip_swath, 1 ) * SIZE( grid%cuprecip_swath, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%cuprecip_swath,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%cuprecip_swath,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%precip_swath, 1 ) * SIZE( grid%precip_swath, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%precip_swath,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%precip_swath,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%windsq_swath, 1 ) * SIZE( grid%windsq_swath, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%windsq_swath,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%windsq_swath,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%suswind, 1 ) * SIZE( grid%suswind, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%suswind,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%suswind,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%wind10_ratio, 1 ) * SIZE( grid%wind10_ratio, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%wind10_ratio,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%wind10_ratio,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%suswind_swath, 1 ) * SIZE( grid%suswind_swath, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%suswind_swath,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%suswind_swath,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%membrane_mslp, 1 ) * SIZE( grid%membrane_mslp, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%membrane_mslp,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%membrane_mslp,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%pdyn_smooth, 1 ) * SIZE( grid%pdyn_smooth, 2 ) .GT. 1 ) THEN 
CONTINUE 
ENDIF
IF ( SIZE( grid%pdyn_parent, 1 ) * SIZE( grid%pdyn_parent, 2 ) .GT. 1 ) THEN 
CALL downaged2d (  &         
                  grid%pdyn_parent,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%pdyn_parent,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,0 &
,ngrid%pdyn_parent_age&
,grid%pdyn_smooth&
                  ) 
ENDIF
IF ( SIZE( grid%best_mslp, 1 ) * SIZE( grid%best_mslp, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%best_mslp,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%best_mslp,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%hlon, 1 ) * SIZE( grid%hlon, 2 ) .GT. 1 ) THEN 
CONTINUE 
ENDIF
IF ( SIZE( grid%hlat, 1 ) * SIZE( grid%hlat, 2 ) .GT. 1 ) THEN 
CONTINUE 
ENDIF
IF ( SIZE( grid%tg_max_m10wind, 1 ) * SIZE( grid%tg_max_m10wind, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%tg_max_m10wind,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%tg_max_m10wind,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%tg_max_wwind, 1 ) * SIZE( grid%tg_max_wwind, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%tg_max_wwind,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%tg_max_wwind,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%tg_min_wwind, 1 ) * SIZE( grid%tg_min_wwind, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%tg_min_wwind,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%tg_min_wwind,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%tg_max_zhel_25, 1 ) * SIZE( grid%tg_max_zhel_25, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%tg_max_zhel_25,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%tg_max_zhel_25,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%tg_min_zhel_25, 1 ) * SIZE( grid%tg_min_zhel_25, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%tg_min_zhel_25,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%tg_min_zhel_25,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%tg_max_zhel_03, 1 ) * SIZE( grid%tg_max_zhel_03, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%tg_max_zhel_03,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%tg_max_zhel_03,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%tg_min_zhel_03, 1 ) * SIZE( grid%tg_min_zhel_03, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%tg_min_zhel_03,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%tg_min_zhel_03,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%tg_updhel25, 1 ) * SIZE( grid%tg_updhel25, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%tg_updhel25,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%tg_updhel25,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%tg_max_updhel25, 1 ) * SIZE( grid%tg_max_updhel25, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%tg_max_updhel25,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%tg_max_updhel25,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%tg_updhel03, 1 ) * SIZE( grid%tg_updhel03, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%tg_updhel03,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%tg_updhel03,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%tg_max_updhel03, 1 ) * SIZE( grid%tg_max_updhel03, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%tg_max_updhel03,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%tg_max_updhel03,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%tg_total_precip, 1 ) * SIZE( grid%tg_total_precip, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%tg_total_precip,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%tg_total_precip,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%hres_fis, 1 ) * SIZE( grid%hres_fis, 2 ) .GT. 1 ) THEN 
CONTINUE 
ENDIF
IF ( SIZE( grid%hres_lnd, 1 ) * SIZE( grid%hres_lnd, 2 ) .GT. 1 ) THEN 
CONTINUE 
ENDIF
IF ( SIZE( grid%hres_lu_index, 1 ) * SIZE( grid%hres_lu_index, 2 ) .GT. 1 ) THEN 
CONTINUE 
ENDIF
IF ( SIZE( grid%sm, 1 ) * SIZE( grid%sm, 2 ) .GT. 1 ) THEN 
CONTINUE 
ENDIF
IF ( SIZE( grid%sice, 1 ) * SIZE( grid%sice, 2 ) .GT. 1 ) THEN 
CALL downnearmaskdef (  &         
                  grid%sice,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%sice,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%iswater,ngrid%iswater&
,0.000000000e+00 &
,0.000000000e+00 &
                  ) 
ENDIF
IF ( SIZE( grid%pd, 1 ) * SIZE( grid%pd, 2 ) .GT. 1 ) THEN 
CONTINUE 
ENDIF
IF ( SIZE( grid%fis, 1 ) * SIZE( grid%fis, 2 ) .GT. 1 ) THEN 
CONTINUE 
ENDIF
IF ( SIZE( grid%t, 1 ) * SIZE( grid%t, 2 ) .GT. 1 ) THEN 
CONTINUE 
ENDIF
IF ( SIZE( grid%q, 1 ) * SIZE( grid%q, 2 ) .GT. 1 ) THEN 
CONTINUE 
ENDIF
IF ( SIZE( grid%u, 1 ) * SIZE( grid%u, 2 ) .GT. 1 ) THEN 
CALL downvel (  &         
                  grid%u,   &       
                 cids, cide, ckds, ckde, cjds, cjde,   &         
                 cims, cime, ckms, ckme, cjms, cjme,   &         
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         
                  ngrid%u,  &   
                 nids, nide, nkds, nkde, njds, njde,   &         
                 nims, nime, nkms, nkme, njms, njme,   &         
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%v, 1 ) * SIZE( grid%v, 2 ) .GT. 1 ) THEN 
CALL downvel (  &         
                  grid%v,   &       
                 cids, cide, ckds, ckde, cjds, cjde,   &         
                 cims, cime, ckms, ckme, cjms, cjme,   &         
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         
                  ngrid%v,  &   
                 nids, nide, nkds, nkde, njds, njde,   &         
                 nims, nime, nkms, nkme, njms, njme,   &         
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%ustar, 1 ) * SIZE( grid%ustar, 2 ) .GT. 1 ) THEN 
CALL downnearmask (  &         
                  grid%ustar,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%ustar,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%iswater,ngrid%iswater&
                  ) 
ENDIF
IF ( SIZE( grid%z0, 1 ) * SIZE( grid%z0, 2 ) .GT. 1 ) THEN 
CALL downnearmask (  &         
                  grid%z0,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%z0,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%iswater,ngrid%iswater&
                  ) 
ENDIF
IF ( SIZE( grid%scurx, 1 ) * SIZE( grid%scurx, 2 ) .GT. 1 ) THEN 
CALL downnearmask (  &         
                  grid%scurx,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%scurx,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%iswater,ngrid%iswater&
                  ) 
ENDIF
IF ( SIZE( grid%scury, 1 ) * SIZE( grid%scury, 2 ) .GT. 1 ) THEN 
CALL downnearmask (  &         
                  grid%scury,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%scury,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%iswater,ngrid%iswater&
                  ) 
ENDIF
IF ( SIZE( grid%charn, 1 ) * SIZE( grid%charn, 2 ) .GT. 1 ) THEN 
CALL downnearmask (  &         
                  grid%charn,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%charn,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%iswater,ngrid%iswater&
                  ) 
ENDIF
IF ( SIZE( grid%msang, 1 ) * SIZE( grid%msang, 2 ) .GT. 1 ) THEN 
CALL downnearmask (  &         
                  grid%msang,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%msang,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%iswater,ngrid%iswater&
                  ) 
ENDIF
IF ( SIZE( grid%ths, 1 ) * SIZE( grid%ths, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%ths,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%ths,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%qsh, 1 ) * SIZE( grid%qsh, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%qsh,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%qsh,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%taux, 1 ) * SIZE( grid%taux, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%taux,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%taux,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%tauy, 1 ) * SIZE( grid%tauy, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%tauy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%tauy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%acprec, 1 ) * SIZE( grid%acprec, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%acprec,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%acprec,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%cuprec, 1 ) * SIZE( grid%cuprec, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%cuprec,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%cuprec,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%cldefi, 1 ) * SIZE( grid%cldefi, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%cldefi,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%cldefi,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%th10, 1 ) * SIZE( grid%th10, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%th10,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%th10,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%q10, 1 ) * SIZE( grid%q10, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%q10,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%q10,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%pshltr, 1 ) * SIZE( grid%pshltr, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%pshltr,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%pshltr,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%tshltr, 1 ) * SIZE( grid%tshltr, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%tshltr,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%tshltr,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%qshltr, 1 ) * SIZE( grid%qshltr, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%qshltr,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%qshltr,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%q2, 1 ) * SIZE( grid%q2, 2 ) .GT. 1 ) THEN 
CALL downmass (  &         
                  grid%q2,   &       
                 cids, cide, ckds, ckde, cjds, cjde,   &         
                 cims, cime, ckms, ckme, cjms, cjme,   &         
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         
                  ngrid%q2,  &   
                 nids, nide, nkds, nkde, njds, njde,   &         
                 nims, nime, nkms, nkme, njms, njme,   &         
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,econst &
,0.000000000e+00 &
                  ) 
ENDIF
IF ( SIZE( grid%t_adj, 1 ) * SIZE( grid%t_adj, 2 ) .GT. 1 ) THEN 
CALL downnear (  &         
                  grid%t_adj,   &       
                 cids, cide, ckds, ckde, cjds, cjde,   &         
                 cims, cime, ckms, ckme, cjms, cjme,   &         
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         
                  ngrid%t_adj,  &   
                 nids, nide, nkds, nkde, njds, njde,   &         
                 nims, nime, nkms, nkme, njms, njme,   &         
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%cd_out, 1 ) * SIZE( grid%cd_out, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%cd_out,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%cd_out,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%ch_out, 1 ) * SIZE( grid%ch_out, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%ch_out,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%ch_out,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%albase, 1 ) * SIZE( grid%albase, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%albase,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%albase,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%nmm_tsk, 1 ) * SIZE( grid%nmm_tsk, 2 ) .GT. 1 ) THEN 
CALL downnearmask (  &         
                  grid%nmm_tsk,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%nmm_tsk,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%iswater,ngrid%iswater&
                  ) 
ENDIF
IF ( SIZE( grid%mxsnal, 1 ) * SIZE( grid%mxsnal, 2 ) .GT. 1 ) THEN 
CALL downnearmask (  &         
                  grid%mxsnal,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%mxsnal,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%iswater,ngrid%iswater&
                  ) 
ENDIF
IF ( SIZE( grid%sigt4, 1 ) * SIZE( grid%sigt4, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%sigt4,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%sigt4,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%tg, 1 ) * SIZE( grid%tg, 2 ) .GT. 1 ) THEN 
CALL downnearmask (  &         
                  grid%tg,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%tg,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%iswater,ngrid%iswater&
                  ) 
ENDIF
IF ( SIZE( grid%acswupt, 1 ) * SIZE( grid%acswupt, 2 ) .GT. 1 ) THEN 
CALL DownCopy (  &         
                  grid%acswupt,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%acswupt,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%acswuptc, 1 ) * SIZE( grid%acswuptc, 2 ) .GT. 1 ) THEN 
CALL DownCopy (  &         
                  grid%acswuptc,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%acswuptc,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%acswdnt, 1 ) * SIZE( grid%acswdnt, 2 ) .GT. 1 ) THEN 
CALL DownCopy (  &         
                  grid%acswdnt,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%acswdnt,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%acswdntc, 1 ) * SIZE( grid%acswdntc, 2 ) .GT. 1 ) THEN 
CALL DownCopy (  &         
                  grid%acswdntc,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%acswdntc,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%acswupb, 1 ) * SIZE( grid%acswupb, 2 ) .GT. 1 ) THEN 
CALL DownCopy (  &         
                  grid%acswupb,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%acswupb,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%acswupbc, 1 ) * SIZE( grid%acswupbc, 2 ) .GT. 1 ) THEN 
CALL DownCopy (  &         
                  grid%acswupbc,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%acswupbc,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%acswdnb, 1 ) * SIZE( grid%acswdnb, 2 ) .GT. 1 ) THEN 
CALL DownCopy (  &         
                  grid%acswdnb,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%acswdnb,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%acswdnbc, 1 ) * SIZE( grid%acswdnbc, 2 ) .GT. 1 ) THEN 
CALL DownCopy (  &         
                  grid%acswdnbc,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%acswdnbc,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%aclwupt, 1 ) * SIZE( grid%aclwupt, 2 ) .GT. 1 ) THEN 
CALL DownCopy (  &         
                  grid%aclwupt,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%aclwupt,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%aclwuptc, 1 ) * SIZE( grid%aclwuptc, 2 ) .GT. 1 ) THEN 
CALL DownCopy (  &         
                  grid%aclwuptc,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%aclwuptc,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%aclwdnt, 1 ) * SIZE( grid%aclwdnt, 2 ) .GT. 1 ) THEN 
CALL DownCopy (  &         
                  grid%aclwdnt,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%aclwdnt,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%aclwdntc, 1 ) * SIZE( grid%aclwdntc, 2 ) .GT. 1 ) THEN 
CALL DownCopy (  &         
                  grid%aclwdntc,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%aclwdntc,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%aclwupb, 1 ) * SIZE( grid%aclwupb, 2 ) .GT. 1 ) THEN 
CALL DownCopy (  &         
                  grid%aclwupb,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%aclwupb,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%aclwupbc, 1 ) * SIZE( grid%aclwupbc, 2 ) .GT. 1 ) THEN 
CALL DownCopy (  &         
                  grid%aclwupbc,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%aclwupbc,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%aclwdnb, 1 ) * SIZE( grid%aclwdnb, 2 ) .GT. 1 ) THEN 
CALL DownCopy (  &         
                  grid%aclwdnb,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%aclwdnb,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%aclwdnbc, 1 ) * SIZE( grid%aclwdnbc, 2 ) .GT. 1 ) THEN 
CALL DownCopy (  &         
                  grid%aclwdnbc,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%aclwdnbc,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%swupt, 1 ) * SIZE( grid%swupt, 2 ) .GT. 1 ) THEN 
CALL DownCopy (  &         
                  grid%swupt,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%swupt,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%swuptc, 1 ) * SIZE( grid%swuptc, 2 ) .GT. 1 ) THEN 
CALL DownCopy (  &         
                  grid%swuptc,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%swuptc,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%swdnt, 1 ) * SIZE( grid%swdnt, 2 ) .GT. 1 ) THEN 
CALL DownCopy (  &         
                  grid%swdnt,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%swdnt,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%swdntc, 1 ) * SIZE( grid%swdntc, 2 ) .GT. 1 ) THEN 
CALL DownCopy (  &         
                  grid%swdntc,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%swdntc,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%swupb, 1 ) * SIZE( grid%swupb, 2 ) .GT. 1 ) THEN 
CALL DownCopy (  &         
                  grid%swupb,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%swupb,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%swupbc, 1 ) * SIZE( grid%swupbc, 2 ) .GT. 1 ) THEN 
CALL DownCopy (  &         
                  grid%swupbc,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%swupbc,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%swdnb, 1 ) * SIZE( grid%swdnb, 2 ) .GT. 1 ) THEN 
CALL DownCopy (  &         
                  grid%swdnb,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%swdnb,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%swdnbc, 1 ) * SIZE( grid%swdnbc, 2 ) .GT. 1 ) THEN 
CALL DownCopy (  &         
                  grid%swdnbc,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%swdnbc,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%lwupt, 1 ) * SIZE( grid%lwupt, 2 ) .GT. 1 ) THEN 
CALL DownCopy (  &         
                  grid%lwupt,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%lwupt,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%lwuptc, 1 ) * SIZE( grid%lwuptc, 2 ) .GT. 1 ) THEN 
CALL DownCopy (  &         
                  grid%lwuptc,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%lwuptc,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%lwdnt, 1 ) * SIZE( grid%lwdnt, 2 ) .GT. 1 ) THEN 
CALL DownCopy (  &         
                  grid%lwdnt,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%lwdnt,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%lwdntc, 1 ) * SIZE( grid%lwdntc, 2 ) .GT. 1 ) THEN 
CALL DownCopy (  &         
                  grid%lwdntc,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%lwdntc,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%lwupb, 1 ) * SIZE( grid%lwupb, 2 ) .GT. 1 ) THEN 
CALL DownCopy (  &         
                  grid%lwupb,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%lwupb,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%lwupbc, 1 ) * SIZE( grid%lwupbc, 2 ) .GT. 1 ) THEN 
CALL DownCopy (  &         
                  grid%lwupbc,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%lwupbc,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%lwdnb, 1 ) * SIZE( grid%lwdnb, 2 ) .GT. 1 ) THEN 
CALL DownCopy (  &         
                  grid%lwdnb,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%lwdnb,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%lwdnbc, 1 ) * SIZE( grid%lwdnbc, 2 ) .GT. 1 ) THEN 
CALL DownCopy (  &         
                  grid%lwdnbc,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%lwdnbc,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%refl_10cm, 1 ) * SIZE( grid%refl_10cm, 3 ) .GT. 1 ) THEN 
CALL downmassikj (  &         
                  grid%refl_10cm,   &       
                 cids, cide, ckds, ckde, cjds, cjde,   &         
                 cims, cime, ckms, ckme, cjms, cjme,   &         
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         
                  ngrid%refl_10cm,  &   
                 nids, nide, nkds, nkde, njds, njde,   &         
                 nims, nime, nkms, nkme, njms, njme,   &         
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,ecopy &
,-3.500000000e+01 &
                  ) 
ENDIF
IF ( SIZE( grid%refd_max, 1 ) * SIZE( grid%refd_max, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%refd_max,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%refd_max,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%qnwfa2d, 1 ) * SIZE( grid%qnwfa2d, 2 ) .GT. 1 ) THEN 
CALL DownCopy (  &         
                  grid%qnwfa2d,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%qnwfa2d,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%swddir, 1 ) * SIZE( grid%swddir, 2 ) .GT. 1 ) THEN 
CALL DownCopy (  &         
                  grid%swddir,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%swddir,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%swddirc, 1 ) * SIZE( grid%swddirc, 2 ) .GT. 1 ) THEN 
CALL DownCopy (  &         
                  grid%swddirc,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%swddirc,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%swddni, 1 ) * SIZE( grid%swddni, 2 ) .GT. 1 ) THEN 
CALL DownCopy (  &         
                  grid%swddni,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%swddni,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%swddnic, 1 ) * SIZE( grid%swddnic, 2 ) .GT. 1 ) THEN 
CALL DownCopy (  &         
                  grid%swddnic,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%swddnic,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%swddif, 1 ) * SIZE( grid%swddif, 2 ) .GT. 1 ) THEN 
CALL DownCopy (  &         
                  grid%swddif,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%swddif,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%gx, 1 ) * SIZE( grid%gx, 2 ) .GT. 1 ) THEN 
CALL DownCopy (  &         
                  grid%gx,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%gx,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%bx, 1 ) * SIZE( grid%bx, 2 ) .GT. 1 ) THEN 
CALL DownCopy (  &         
                  grid%bx,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%bx,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%gg, 1 ) * SIZE( grid%gg, 2 ) .GT. 1 ) THEN 
CALL DownCopy (  &         
                  grid%gg,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%gg,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%bb, 1 ) * SIZE( grid%bb, 2 ) .GT. 1 ) THEN 
CALL DownCopy (  &         
                  grid%bb,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%bb,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%coszen_ref, 1 ) * SIZE( grid%coszen_ref, 2 ) .GT. 1 ) THEN 
CALL DownCopy (  &         
                  grid%coszen_ref,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%coszen_ref,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%swdown_ref, 1 ) * SIZE( grid%swdown_ref, 2 ) .GT. 1 ) THEN 
CALL DownCopy (  &         
                  grid%swdown_ref,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%swdown_ref,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%swddir_ref, 1 ) * SIZE( grid%swddir_ref, 2 ) .GT. 1 ) THEN 
CALL DownCopy (  &         
                  grid%swddir_ref,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%swddir_ref,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%cwm, 1 ) * SIZE( grid%cwm, 2 ) .GT. 1 ) THEN 
CALL downmass (  &         
                  grid%cwm,   &       
                 cids, cide, ckds, ckde, cjds, cjde,   &         
                 cims, cime, ckms, ckme, cjms, cjme,   &         
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         
                  ngrid%cwm,  &   
                 nids, nide, nkds, nkde, njds, njde,   &         
                 nims, nime, nkms, nkme, njms, njme,   &         
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,ecopy &
,0.000000000e+00 &
                  ) 
ENDIF
IF ( SIZE( grid%f_ice, 1 ) * SIZE( grid%f_ice, 3 ) .GT. 1 ) THEN 
CALL downmassikj (  &         
                  grid%f_ice,   &       
                 cids, cide, ckds, ckde, cjds, cjde,   &         
                 cims, cime, ckms, ckme, cjms, cjme,   &         
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         
                  ngrid%f_ice,  &   
                 nids, nide, nkds, nkde, njds, njde,   &         
                 nims, nime, nkms, nkme, njms, njme,   &         
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,eextrap &
,0.000000000e+00 &
                  ) 
ENDIF
IF ( SIZE( grid%f_rain, 1 ) * SIZE( grid%f_rain, 3 ) .GT. 1 ) THEN 
CALL downmassikj (  &         
                  grid%f_rain,   &       
                 cids, cide, ckds, ckde, cjds, cjde,   &         
                 cims, cime, ckms, ckme, cjms, cjme,   &         
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         
                  ngrid%f_rain,  &   
                 nids, nide, nkds, nkde, njds, njde,   &         
                 nims, nime, nkms, nkme, njms, njme,   &         
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,eextrap &
,0.000000000e+00 &
                  ) 
ENDIF
IF ( SIZE( grid%f_rimef, 1 ) * SIZE( grid%f_rimef, 3 ) .GT. 1 ) THEN 
CALL downmassikj (  &         
                  grid%f_rimef,   &       
                 cids, cide, ckds, ckde, cjds, cjde,   &         
                 cims, cime, ckms, ckme, cjms, cjme,   &         
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         
                  ngrid%f_rimef,  &   
                 nids, nide, nkds, nkde, njds, njde,   &         
                 nims, nime, nkms, nkme, njms, njme,   &         
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,eextrap &
,1.000000000e+00 &
                  ) 
ENDIF
IF ( SIZE( grid%cfrach, 1 ) * SIZE( grid%cfrach, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%cfrach,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%cfrach,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%cfracl, 1 ) * SIZE( grid%cfracl, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%cfracl,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%cfracl,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%cfracm, 1 ) * SIZE( grid%cfracm, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%cfracm,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%cfracm,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%islope, 1 ) * SIZE( grid%islope, 2 ) .GT. 1 ) THEN 
CALL downinearmask (  &         
                  grid%islope,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%islope,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%iswater,ngrid%iswater&
                  ) 
ENDIF
IF ( SIZE( grid%cmc, 1 ) * SIZE( grid%cmc, 2 ) .GT. 1 ) THEN 
CALL downnearmask (  &         
                  grid%cmc,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%cmc,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%iswater,ngrid%iswater&
                  ) 
ENDIF
IF ( SIZE( grid%soiltb, 1 ) * SIZE( grid%soiltb, 2 ) .GT. 1 ) THEN 
CALL downnearmask (  &         
                  grid%soiltb,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%soiltb,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%iswater,ngrid%iswater&
                  ) 
ENDIF
IF ( SIZE( grid%vegfrc, 1 ) * SIZE( grid%vegfrc, 2 ) .GT. 1 ) THEN 
CALL downnearmask (  &         
                  grid%vegfrc,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%vegfrc,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%iswater,ngrid%iswater&
                  ) 
ENDIF
IF ( SIZE( grid%shdmax, 1 ) * SIZE( grid%shdmax, 2 ) .GT. 1 ) THEN 
CALL downnearmask (  &         
                  grid%shdmax,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%shdmax,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%iswater,ngrid%iswater&
                  ) 
ENDIF
IF ( SIZE( grid%shdmin, 1 ) * SIZE( grid%shdmin, 2 ) .GT. 1 ) THEN 
CALL downnearmask (  &         
                  grid%shdmin,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%shdmin,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%iswater,ngrid%iswater&
                  ) 
ENDIF
IF ( SIZE( grid%sh2o, 1 ) * SIZE( grid%sh2o, 3 ) .GT. 1 ) THEN 
CALL downnearmaskikj (  &         
                  grid%sh2o,   &       
                 cids, cide, 1, config_flags%num_soil_layers, cjds, cjde,   &         
                 cims, cime, 1, config_flags%num_soil_layers, cjms, cjme,   &         
                 cips, cipe, 1, config_flags%num_soil_layers, cjps, cjpe,   &         
                  ngrid%sh2o,  &   
                 nids, nide, 1, config_flags%num_soil_layers, njds, njde,   &         
                 nims, nime, 1, config_flags%num_soil_layers, njms, njme,   &         
                 nips, nipe, 1, config_flags%num_soil_layers, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%iswater,ngrid%iswater&
                  ) 
ENDIF
IF ( SIZE( grid%smc, 1 ) * SIZE( grid%smc, 3 ) .GT. 1 ) THEN 
CALL downnearmaskikj (  &         
                  grid%smc,   &       
                 cids, cide, 1, config_flags%num_soil_layers, cjds, cjde,   &         
                 cims, cime, 1, config_flags%num_soil_layers, cjms, cjme,   &         
                 cips, cipe, 1, config_flags%num_soil_layers, cjps, cjpe,   &         
                  ngrid%smc,  &   
                 nids, nide, 1, config_flags%num_soil_layers, njds, njde,   &         
                 nims, nime, 1, config_flags%num_soil_layers, njms, njme,   &         
                 nips, nipe, 1, config_flags%num_soil_layers, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%iswater,ngrid%iswater&
                  ) 
ENDIF
IF ( SIZE( grid%stc, 1 ) * SIZE( grid%stc, 3 ) .GT. 1 ) THEN 
CALL downnearmaskikj (  &         
                  grid%stc,   &       
                 cids, cide, 1, config_flags%num_soil_layers, cjds, cjde,   &         
                 cims, cime, 1, config_flags%num_soil_layers, cjms, cjme,   &         
                 cips, cipe, 1, config_flags%num_soil_layers, cjps, cjpe,   &         
                  ngrid%stc,  &   
                 nids, nide, 1, config_flags%num_soil_layers, njds, njde,   &         
                 nims, nime, 1, config_flags%num_soil_layers, njms, njme,   &         
                 nips, nipe, 1, config_flags%num_soil_layers, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%iswater,ngrid%iswater&
                  ) 
ENDIF
IF ( SIZE( grid%ctopo, 1 ) * SIZE( grid%ctopo, 2 ) .GT. 1 ) THEN 
CALL DownCopy (  &         
                  grid%ctopo,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%ctopo,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%ctopo2, 1 ) * SIZE( grid%ctopo2, 2 ) .GT. 1 ) THEN 
CALL DownCopy (  &         
                  grid%ctopo2,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%ctopo2,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%dwdt, 1 ) * SIZE( grid%dwdt, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%dwdt,   &       
                 cids, cide, ckds, ckde, cjds, cjde,   &         
                 cims, cime, ckms, ckme, cjms, cjme,   &         
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         
                  ngrid%dwdt,  &   
                 nids, nide, nkds, nkde, njds, njde,   &         
                 nims, nime, nkms, nkme, njms, njme,   &         
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%pint, 1 ) * SIZE( grid%pint, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%pint,   &       
                 cids, cide, ckds, ckde, cjds, cjde,   &         
                 cims, cime, ckms, ckme, cjms, cjme,   &         
                 cips, cipe, ckps, MIN( ckde, ckpe ), cjps, cjpe,   &         
                  ngrid%pint,  &   
                 nids, nide, nkds, nkde, njds, njde,   &         
                 nims, nime, nkms, nkme, njms, njme,   &         
                 nips, nipe, nkps, MIN( nkde, nkpe ), njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%w, 1 ) * SIZE( grid%w, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%w,   &       
                 cids, cide, ckds, ckde, cjds, cjde,   &         
                 cims, cime, ckms, ckme, cjms, cjme,   &         
                 cips, cipe, ckps, MIN( ckde, ckpe ), cjps, cjpe,   &         
                  ngrid%w,  &   
                 nids, nide, nkds, nkde, njds, njde,   &         
                 nims, nime, nkms, nkme, njms, njme,   &         
                 nips, nipe, nkps, MIN( nkde, nkpe ), njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%w_tot, 1 ) * SIZE( grid%w_tot, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%w_tot,   &       
                 cids, cide, ckds, ckde, cjds, cjde,   &         
                 cims, cime, ckms, ckme, cjms, cjme,   &         
                 cips, cipe, ckps, MIN( ckde, ckpe ), cjps, cjpe,   &         
                  ngrid%w_tot,  &   
                 nids, nide, nkds, nkde, njds, njde,   &         
                 nims, nime, nkms, nkme, njms, njme,   &         
                 nips, nipe, nkps, MIN( nkde, nkpe ), njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%z, 1 ) * SIZE( grid%z, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%z,   &       
                 cids, cide, ckds, ckde, cjds, cjde,   &         
                 cims, cime, ckms, ckme, cjms, cjme,   &         
                 cips, cipe, ckps, MIN( ckde, ckpe ), cjps, cjpe,   &         
                  ngrid%z,  &   
                 nids, nide, nkds, nkde, njds, njde,   &         
                 nims, nime, nkms, nkme, njms, njme,   &         
                 nips, nipe, nkps, MIN( nkde, nkpe ), njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%rlwin, 1 ) * SIZE( grid%rlwin, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%rlwin,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%rlwin,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%rswin, 1 ) * SIZE( grid%rswin, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%rswin,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%rswin,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%rlwtt, 1 ) * SIZE( grid%rlwtt, 2 ) .GT. 1 ) THEN 
CALL downnear (  &         
                  grid%rlwtt,   &       
                 cids, cide, ckds, ckde, cjds, cjde,   &         
                 cims, cime, ckms, ckme, cjms, cjme,   &         
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         
                  ngrid%rlwtt,  &   
                 nids, nide, nkds, nkde, njds, njde,   &         
                 nims, nime, nkms, nkme, njms, njme,   &         
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%rswtt, 1 ) * SIZE( grid%rswtt, 2 ) .GT. 1 ) THEN 
CALL downnear (  &         
                  grid%rswtt,   &       
                 cids, cide, ckds, ckde, cjds, cjde,   &         
                 cims, cime, ckms, ckme, cjms, cjme,   &         
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         
                  ngrid%rswtt,  &   
                 nids, nide, nkds, nkde, njds, njde,   &         
                 nims, nime, nkms, nkme, njms, njme,   &         
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%winfo, 1 ) * SIZE( grid%winfo, 2 ) .GT. 1 ) THEN 
CONTINUE 
ENDIF
IF ( SIZE( grid%iinfo, 1 ) * SIZE( grid%iinfo, 2 ) .GT. 1 ) THEN 
CONTINUE 
ENDIF
IF ( SIZE( grid%landmask, 1 ) * SIZE( grid%landmask, 2 ) .GT. 1 ) THEN 
CONTINUE 
ENDIF
IF ( SIZE( grid%toposoil, 1 ) * SIZE( grid%toposoil, 2 ) .GT. 1 ) THEN 
CALL downnearmask (  &         
                  grid%toposoil,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%toposoil,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%iswater,ngrid%iswater&
                  ) 
ENDIF
DO itrace = PARAM_FIRST_SCALAR, num_moist
IF ( SIZE( moist, 1 ) * SIZE( moist, 2 ) .GT. 1 .and. (interp_mp .eqv. .true.) ) THEN 
CALL downmass (  &         
                  moist(grid%sm31,grid%sm32,grid%sm33,itrace),   &       
                 cids, cide, ckds, ckde, cjds, cjde,   &         
                 cims, cime, ckms, ckme, cjms, cjme,   &         
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         
                  ngrid%moist(ngrid%sm31,ngrid%sm32,ngrid%sm33,itrace),  &   
                 nids, nide, nkds, nkde, njds, njde,   &         
                 nims, nime, nkms, nkme, njms, njme,   &         
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,ecopy &
,0.000000000e+00 &
                  ) 
ENDIF
ENDDO
DO itrace = PARAM_FIRST_SCALAR, num_scalar
IF ( SIZE( scalar, 1 ) * SIZE( scalar, 2 ) .GT. 1 .and. (interp_mp .eqv. .true.) ) THEN 
CALL downmass (  &         
                  scalar(grid%sm31,grid%sm32,grid%sm33,itrace),   &       
                 cids, cide, ckds, ckde, cjds, cjde,   &         
                 cims, cime, ckms, ckme, cjms, cjme,   &         
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         
                  ngrid%scalar(ngrid%sm31,ngrid%sm32,ngrid%sm33,itrace),  &   
                 nids, nide, nkds, nkde, njds, njde,   &         
                 nims, nime, nkms, nkme, njms, njme,   &         
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,ecopy &
,0.000000000e+00 &
                  ) 
ENDIF
ENDDO
DO itrace = PARAM_FIRST_SCALAR, num_dfi_scalar
IF ( SIZE( dfi_scalar, 1 ) * SIZE( dfi_scalar, 2 ) .GT. 1 .and. (interp_mp .eqv. .true.) ) THEN 
CALL downmass (  &         
                  dfi_scalar(grid%sm31,grid%sm32,grid%sm33,itrace),   &       
                 cids, cide, ckds, ckde, cjds, cjde,   &         
                 cims, cime, ckms, ckme, cjms, cjme,   &         
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         
                  ngrid%dfi_scalar(ngrid%sm31,ngrid%sm32,ngrid%sm33,itrace),  &   
                 nids, nide, nkds, nkde, njds, njde,   &         
                 nims, nime, nkms, nkme, njms, njme,   &         
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,ecopy &
,0.000000000e+00 &
                  ) 
ENDIF
ENDDO
IF ( .TRUE. ) THEN 
CALL interp_mask_water_field (  &         
  ( SIZE( grid%lake_depth , 1 )*SIZE( grid%lake_depth , 2 ) .GT. 1 ), & 
                  grid%lake_depth,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%lake_depth,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%iswater,ngrid%iswater&
                  ) 
ENDIF
IF ( SIZE( grid%u10, 1 ) * SIZE( grid%u10, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%u10,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%u10,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%v10, 1 ) * SIZE( grid%v10, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%v10,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%v10,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%xice, 1 ) * SIZE( grid%xice, 2 ) .GT. 1 ) THEN 
CALL downnearmask (  &         
                  grid%xice,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%xice,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%iswater,ngrid%iswater&
                  ) 
ENDIF
IF ( SIZE( grid%icedepth, 1 ) * SIZE( grid%icedepth, 2 ) .GT. 1 ) THEN 
CALL downnearmask (  &         
                  grid%icedepth,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%icedepth,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%iswater,ngrid%iswater&
                  ) 
ENDIF
IF ( SIZE( grid%albsi, 1 ) * SIZE( grid%albsi, 2 ) .GT. 1 ) THEN 
CALL downnearmask (  &         
                  grid%albsi,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%albsi,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%iswater,ngrid%iswater&
                  ) 
ENDIF
IF ( SIZE( grid%snowsi, 1 ) * SIZE( grid%snowsi, 2 ) .GT. 1 ) THEN 
CALL downnearmask (  &         
                  grid%snowsi,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%snowsi,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%iswater,ngrid%iswater&
                  ) 
ENDIF
IF ( SIZE( grid%ivgtyp, 1 ) * SIZE( grid%ivgtyp, 2 ) .GT. 1 ) THEN 
CALL downinearmask (  &         
                  grid%ivgtyp,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%ivgtyp,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%iswater,ngrid%iswater&
                  ) 
ENDIF
IF ( SIZE( grid%isltyp, 1 ) * SIZE( grid%isltyp, 2 ) .GT. 1 ) THEN 
CALL downinearmaskdef (  &         
                  grid%isltyp,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%isltyp,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%iswater,ngrid%iswater&
,14 &
,8 &
                  ) 
ENDIF
IF ( SIZE( grid%vegfra, 1 ) * SIZE( grid%vegfra, 2 ) .GT. 1 ) THEN 
CALL downnearmask (  &         
                  grid%vegfra,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%vegfra,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%iswater,ngrid%iswater&
                  ) 
ENDIF
IF ( SIZE( grid%sst, 1 ) * SIZE( grid%sst, 2 ) .GT. 1 ) THEN 
CALL downnearmask (  &         
                  grid%sst,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%sst,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%iswater,ngrid%iswater&
                  ) 
ENDIF
IF ( SIZE( grid%weasd, 1 ) * SIZE( grid%weasd, 2 ) .GT. 1 ) THEN 
CALL downnearmask (  &         
                  grid%weasd,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%weasd,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
,grid%iswater,ngrid%iswater&
                  ) 
ENDIF
IF ( SIZE( grid%thz0, 1 ) * SIZE( grid%thz0, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%thz0,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%thz0,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%qz0, 1 ) * SIZE( grid%qz0, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%qz0,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%qz0,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%uz0, 1 ) * SIZE( grid%uz0, 2 ) .GT. 1 ) THEN 
CALL downvel (  &         
                  grid%uz0,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%uz0,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%vz0, 1 ) * SIZE( grid%vz0, 2 ) .GT. 1 ) THEN 
CALL downvel (  &         
                  grid%vz0,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%vz0,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%htop, 1 ) * SIZE( grid%htop, 2 ) .GT. 1 ) THEN 
CALL downnear (  &         
                  grid%htop,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%htop,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%hbot, 1 ) * SIZE( grid%hbot, 2 ) .GT. 1 ) THEN 
CALL downnear (  &         
                  grid%hbot,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%hbot,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%htopr, 1 ) * SIZE( grid%htopr, 2 ) .GT. 1 ) THEN 
CALL downnear (  &         
                  grid%htopr,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%htopr,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%hbotr, 1 ) * SIZE( grid%hbotr, 2 ) .GT. 1 ) THEN 
CALL downnear (  &         
                  grid%hbotr,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%hbotr,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%cuppt, 1 ) * SIZE( grid%cuppt, 2 ) .GT. 1 ) THEN 
CALL downnear (  &         
                  grid%cuppt,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%cuppt,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%snowh, 1 ) * SIZE( grid%snowh, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%snowh,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%snowh,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%rhosn, 1 ) * SIZE( grid%rhosn, 2 ) .GT. 1 ) THEN 
CALL downcopy (  &         
                  grid%rhosn,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%rhosn,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%isnowxy , 1 )*SIZE( grid%isnowxy , 2 ) .GT. 1 ), & 
                  grid%isnowxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%isnowxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%tvxy , 1 )*SIZE( grid%tvxy , 2 ) .GT. 1 ), & 
                  grid%tvxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%tvxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%tgxy , 1 )*SIZE( grid%tgxy , 2 ) .GT. 1 ), & 
                  grid%tgxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%tgxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%canicexy , 1 )*SIZE( grid%canicexy , 2 ) .GT. 1 ), & 
                  grid%canicexy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%canicexy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%canliqxy , 1 )*SIZE( grid%canliqxy , 2 ) .GT. 1 ), & 
                  grid%canliqxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%canliqxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%eahxy , 1 )*SIZE( grid%eahxy , 2 ) .GT. 1 ), & 
                  grid%eahxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%eahxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%tahxy , 1 )*SIZE( grid%tahxy , 2 ) .GT. 1 ), & 
                  grid%tahxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%tahxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%cmxy , 1 )*SIZE( grid%cmxy , 2 ) .GT. 1 ), & 
                  grid%cmxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%cmxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%chxy , 1 )*SIZE( grid%chxy , 2 ) .GT. 1 ), & 
                  grid%chxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%chxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%fwetxy , 1 )*SIZE( grid%fwetxy , 2 ) .GT. 1 ), & 
                  grid%fwetxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%fwetxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%sneqvoxy , 1 )*SIZE( grid%sneqvoxy , 2 ) .GT. 1 ), & 
                  grid%sneqvoxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%sneqvoxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%alboldxy , 1 )*SIZE( grid%alboldxy , 2 ) .GT. 1 ), & 
                  grid%alboldxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%alboldxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%qsnowxy , 1 )*SIZE( grid%qsnowxy , 2 ) .GT. 1 ), & 
                  grid%qsnowxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%qsnowxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%wslakexy , 1 )*SIZE( grid%wslakexy , 2 ) .GT. 1 ), & 
                  grid%wslakexy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%wslakexy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%zwtxy , 1 )*SIZE( grid%zwtxy , 2 ) .GT. 1 ), & 
                  grid%zwtxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%zwtxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%waxy , 1 )*SIZE( grid%waxy , 2 ) .GT. 1 ), & 
                  grid%waxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%waxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%wtxy , 1 )*SIZE( grid%wtxy , 2 ) .GT. 1 ), & 
                  grid%wtxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%wtxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%tsnoxy , 1 )*SIZE( grid%tsnoxy , 3 ) .GT. 1 ), & 
                  grid%tsnoxy,   &       
                 cids, cide, 1, config_flags%num_snow_layers, cjds, cjde,   &         
                 cims, cime, 1, config_flags%num_snow_layers, cjms, cjme,   &         
                 cips, cipe, 1, config_flags%num_snow_layers, cjps, cjpe,   &         
                  ngrid%tsnoxy,  &   
                 nids, nide, 1, config_flags%num_snow_layers, njds, njde,   &         
                 nims, nime, 1, config_flags%num_snow_layers, njms, njme,   &         
                 nips, nipe, 1, config_flags%num_snow_layers, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%zsnsoxy , 1 )*SIZE( grid%zsnsoxy , 3 ) .GT. 1 ), & 
                  grid%zsnsoxy,   &       
                 cids, cide, 1, config_flags%num_snso_layers, cjds, cjde,   &         
                 cims, cime, 1, config_flags%num_snso_layers, cjms, cjme,   &         
                 cips, cipe, 1, config_flags%num_snso_layers, cjps, cjpe,   &         
                  ngrid%zsnsoxy,  &   
                 nids, nide, 1, config_flags%num_snso_layers, njds, njde,   &         
                 nims, nime, 1, config_flags%num_snso_layers, njms, njme,   &         
                 nips, nipe, 1, config_flags%num_snso_layers, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%snicexy , 1 )*SIZE( grid%snicexy , 3 ) .GT. 1 ), & 
                  grid%snicexy,   &       
                 cids, cide, 1, config_flags%num_snow_layers, cjds, cjde,   &         
                 cims, cime, 1, config_flags%num_snow_layers, cjms, cjme,   &         
                 cips, cipe, 1, config_flags%num_snow_layers, cjps, cjpe,   &         
                  ngrid%snicexy,  &   
                 nids, nide, 1, config_flags%num_snow_layers, njds, njde,   &         
                 nims, nime, 1, config_flags%num_snow_layers, njms, njme,   &         
                 nips, nipe, 1, config_flags%num_snow_layers, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%snliqxy , 1 )*SIZE( grid%snliqxy , 3 ) .GT. 1 ), & 
                  grid%snliqxy,   &       
                 cids, cide, 1, config_flags%num_snow_layers, cjds, cjde,   &         
                 cims, cime, 1, config_flags%num_snow_layers, cjms, cjme,   &         
                 cips, cipe, 1, config_flags%num_snow_layers, cjps, cjpe,   &         
                  ngrid%snliqxy,  &   
                 nids, nide, 1, config_flags%num_snow_layers, njds, njde,   &         
                 nims, nime, 1, config_flags%num_snow_layers, njms, njme,   &         
                 nips, nipe, 1, config_flags%num_snow_layers, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%lfmassxy , 1 )*SIZE( grid%lfmassxy , 2 ) .GT. 1 ), & 
                  grid%lfmassxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%lfmassxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%rtmassxy , 1 )*SIZE( grid%rtmassxy , 2 ) .GT. 1 ), & 
                  grid%rtmassxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%rtmassxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%stmassxy , 1 )*SIZE( grid%stmassxy , 2 ) .GT. 1 ), & 
                  grid%stmassxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%stmassxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%woodxy , 1 )*SIZE( grid%woodxy , 2 ) .GT. 1 ), & 
                  grid%woodxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%woodxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%stblcpxy , 1 )*SIZE( grid%stblcpxy , 2 ) .GT. 1 ), & 
                  grid%stblcpxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%stblcpxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%fastcpxy , 1 )*SIZE( grid%fastcpxy , 2 ) .GT. 1 ), & 
                  grid%fastcpxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%fastcpxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%xsaixy , 1 )*SIZE( grid%xsaixy , 2 ) .GT. 1 ), & 
                  grid%xsaixy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%xsaixy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%t2mvxy , 1 )*SIZE( grid%t2mvxy , 2 ) .GT. 1 ), & 
                  grid%t2mvxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%t2mvxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%t2mbxy , 1 )*SIZE( grid%t2mbxy , 2 ) .GT. 1 ), & 
                  grid%t2mbxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%t2mbxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%q2mvxy , 1 )*SIZE( grid%q2mvxy , 2 ) .GT. 1 ), & 
                  grid%q2mvxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%q2mvxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%q2mbxy , 1 )*SIZE( grid%q2mbxy , 2 ) .GT. 1 ), & 
                  grid%q2mbxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%q2mbxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%tradxy , 1 )*SIZE( grid%tradxy , 2 ) .GT. 1 ), & 
                  grid%tradxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%tradxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%neexy , 1 )*SIZE( grid%neexy , 2 ) .GT. 1 ), & 
                  grid%neexy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%neexy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%gppxy , 1 )*SIZE( grid%gppxy , 2 ) .GT. 1 ), & 
                  grid%gppxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%gppxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%nppxy , 1 )*SIZE( grid%nppxy , 2 ) .GT. 1 ), & 
                  grid%nppxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%nppxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%fvegxy , 1 )*SIZE( grid%fvegxy , 2 ) .GT. 1 ), & 
                  grid%fvegxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%fvegxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%qinxy , 1 )*SIZE( grid%qinxy , 2 ) .GT. 1 ), & 
                  grid%qinxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%qinxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%runsfxy , 1 )*SIZE( grid%runsfxy , 2 ) .GT. 1 ), & 
                  grid%runsfxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%runsfxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%runsbxy , 1 )*SIZE( grid%runsbxy , 2 ) .GT. 1 ), & 
                  grid%runsbxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%runsbxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%ecanxy , 1 )*SIZE( grid%ecanxy , 2 ) .GT. 1 ), & 
                  grid%ecanxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%ecanxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%edirxy , 1 )*SIZE( grid%edirxy , 2 ) .GT. 1 ), & 
                  grid%edirxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%edirxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%etranxy , 1 )*SIZE( grid%etranxy , 2 ) .GT. 1 ), & 
                  grid%etranxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%etranxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%fsaxy , 1 )*SIZE( grid%fsaxy , 2 ) .GT. 1 ), & 
                  grid%fsaxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%fsaxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%firaxy , 1 )*SIZE( grid%firaxy , 2 ) .GT. 1 ), & 
                  grid%firaxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%firaxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%aparxy , 1 )*SIZE( grid%aparxy , 2 ) .GT. 1 ), & 
                  grid%aparxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%aparxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%psnxy , 1 )*SIZE( grid%psnxy , 2 ) .GT. 1 ), & 
                  grid%psnxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%psnxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%savxy , 1 )*SIZE( grid%savxy , 2 ) .GT. 1 ), & 
                  grid%savxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%savxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%sagxy , 1 )*SIZE( grid%sagxy , 2 ) .GT. 1 ), & 
                  grid%sagxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%sagxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%rssunxy , 1 )*SIZE( grid%rssunxy , 2 ) .GT. 1 ), & 
                  grid%rssunxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%rssunxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%rsshaxy , 1 )*SIZE( grid%rsshaxy , 2 ) .GT. 1 ), & 
                  grid%rsshaxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%rsshaxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%bgapxy , 1 )*SIZE( grid%bgapxy , 2 ) .GT. 1 ), & 
                  grid%bgapxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%bgapxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%wgapxy , 1 )*SIZE( grid%wgapxy , 2 ) .GT. 1 ), & 
                  grid%wgapxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%wgapxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%tgvxy , 1 )*SIZE( grid%tgvxy , 2 ) .GT. 1 ), & 
                  grid%tgvxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%tgvxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%tgbxy , 1 )*SIZE( grid%tgbxy , 2 ) .GT. 1 ), & 
                  grid%tgbxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%tgbxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%chvxy , 1 )*SIZE( grid%chvxy , 2 ) .GT. 1 ), & 
                  grid%chvxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%chvxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%chbxy , 1 )*SIZE( grid%chbxy , 2 ) .GT. 1 ), & 
                  grid%chbxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%chbxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%shgxy , 1 )*SIZE( grid%shgxy , 2 ) .GT. 1 ), & 
                  grid%shgxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%shgxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%shcxy , 1 )*SIZE( grid%shcxy , 2 ) .GT. 1 ), & 
                  grid%shcxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%shcxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%shbxy , 1 )*SIZE( grid%shbxy , 2 ) .GT. 1 ), & 
                  grid%shbxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%shbxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%evgxy , 1 )*SIZE( grid%evgxy , 2 ) .GT. 1 ), & 
                  grid%evgxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%evgxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%evbxy , 1 )*SIZE( grid%evbxy , 2 ) .GT. 1 ), & 
                  grid%evbxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%evbxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%ghvxy , 1 )*SIZE( grid%ghvxy , 2 ) .GT. 1 ), & 
                  grid%ghvxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%ghvxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%ghbxy , 1 )*SIZE( grid%ghbxy , 2 ) .GT. 1 ), & 
                  grid%ghbxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%ghbxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%irgxy , 1 )*SIZE( grid%irgxy , 2 ) .GT. 1 ), & 
                  grid%irgxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%irgxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%ircxy , 1 )*SIZE( grid%ircxy , 2 ) .GT. 1 ), & 
                  grid%ircxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%ircxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%irbxy , 1 )*SIZE( grid%irbxy , 2 ) .GT. 1 ), & 
                  grid%irbxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%irbxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%trxy , 1 )*SIZE( grid%trxy , 2 ) .GT. 1 ), & 
                  grid%trxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%trxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%evcxy , 1 )*SIZE( grid%evcxy , 2 ) .GT. 1 ), & 
                  grid%evcxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%evcxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%chleafxy , 1 )*SIZE( grid%chleafxy , 2 ) .GT. 1 ), & 
                  grid%chleafxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%chleafxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%chucxy , 1 )*SIZE( grid%chucxy , 2 ) .GT. 1 ), & 
                  grid%chucxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%chucxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%chv2xy , 1 )*SIZE( grid%chv2xy , 2 ) .GT. 1 ), & 
                  grid%chv2xy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%chv2xy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%chb2xy , 1 )*SIZE( grid%chb2xy , 2 ) .GT. 1 ), & 
                  grid%chb2xy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%chb2xy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%chstarxy , 1 )*SIZE( grid%chstarxy , 2 ) .GT. 1 ), & 
                  grid%chstarxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%chstarxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%soilcomp , 1 )*SIZE( grid%soilcomp , 3 ) .GT. 1 ), & 
                  grid%soilcomp,   &       
                 cids, cide, 1, 8, cjds, cjde,   &         
                 cims, cime, 1, 8, cjms, cjme,   &         
                 cips, cipe, 1, 8, cjps, cjpe,   &         
                  ngrid%soilcomp,  &   
                 nids, nide, 1, 8, njds, njde,   &         
                 nims, nime, 1, 8, njms, njme,   &         
                 nips, nipe, 1, 8, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%soilcl1 , 1 )*SIZE( grid%soilcl1 , 2 ) .GT. 1 ), & 
                  grid%soilcl1,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%soilcl1,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%soilcl2 , 1 )*SIZE( grid%soilcl2 , 2 ) .GT. 1 ), & 
                  grid%soilcl2,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%soilcl2,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%soilcl3 , 1 )*SIZE( grid%soilcl3 , 2 ) .GT. 1 ), & 
                  grid%soilcl3,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%soilcl3,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%soilcl4 , 1 )*SIZE( grid%soilcl4 , 2 ) .GT. 1 ), & 
                  grid%soilcl4,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%soilcl4,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%grainxy , 1 )*SIZE( grid%grainxy , 2 ) .GT. 1 ), & 
                  grid%grainxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%grainxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%gddxy , 1 )*SIZE( grid%gddxy , 2 ) .GT. 1 ), & 
                  grid%gddxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%gddxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%croptype , 1 )*SIZE( grid%croptype , 3 ) .GT. 1 ), & 
                  grid%croptype,   &       
                 cids, cide, 1, 5, cjds, cjde,   &         
                 cims, cime, 1, 5, cjms, cjme,   &         
                 cips, cipe, 1, 5, cjps, cjpe,   &         
                  ngrid%croptype,  &   
                 nids, nide, 1, 5, njds, njde,   &         
                 nims, nime, 1, 5, njms, njme,   &         
                 nips, nipe, 1, 5, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%planting , 1 )*SIZE( grid%planting , 2 ) .GT. 1 ), & 
                  grid%planting,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%planting,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%harvest , 1 )*SIZE( grid%harvest , 2 ) .GT. 1 ), & 
                  grid%harvest,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%harvest,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%season_gdd , 1 )*SIZE( grid%season_gdd , 2 ) .GT. 1 ), & 
                  grid%season_gdd,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%season_gdd,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF
IF ( .TRUE. ) THEN 
CALL interp_mask_land_field (  &         
  ( SIZE( grid%gecros_state , 1 )*SIZE( grid%gecros_state , 3 ) .GT. 1 ), & 
                  grid%gecros_state,   &       
                 cids, cide, 1, 60, cjds, cjde,   &         
                 cims, cime, 1, 60, cjms, cjme,   &         
                 cips, cipe, 1, 60, cjps, cjpe,   &         
                  ngrid%gecros_state,  &   
                 nids, nide, 1, 60, njds, njde,   &         
                 nims, nime, 1, 60, njms, njme,   &         
                 nips, nipe, 1, 60, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%lu_index,ngrid%lu_index&
                  ) 
ENDIF



        CALL pop_communicators_for_domain

      END IF

      RETURN
   END SUBROUTINE interp_domain_nmm_part2



   SUBROUTINE force_domain_nmm_part1 ( grid, intermediate_grid, ngrid, config_flags    &







,szj,s1z,spz,tcs,moist,moist_bxs,moist_bxe,moist_bys,moist_bye,moist_btxs,moist_btxe,moist_btys,moist_btye,dfi_moist, &
dfi_moist_bxs,dfi_moist_bxe,dfi_moist_bys,dfi_moist_bye,dfi_moist_btxs,dfi_moist_btxe,dfi_moist_btys,dfi_moist_btye,scalar, &
scalar_bxs,scalar_bxe,scalar_bys,scalar_bye,scalar_btxs,scalar_btxe,scalar_btys,scalar_btye,dfi_scalar,dfi_scalar_bxs, &
dfi_scalar_bxe,dfi_scalar_bys,dfi_scalar_bye,dfi_scalar_btxs,dfi_scalar_btxe,dfi_scalar_btys,dfi_scalar_btye,chem,ozmixm &


                 )
      USE module_state_description
      USE module_domain, ONLY : domain, get_ijk_from_grid
      USE module_configure, ONLY : grid_config_rec_type
      USE module_dm, ONLY : ntasks, ntasks_x, ntasks_y, itrace, local_communicator, mytask, &
                            ipe_save, jpe_save, ips_save, jps_save, get_dm_max_halo_width,  &
                            nest_pes_x, nest_pes_y,                                 &
                            intercomm_active, nest_task_offsets,                    &
                            mpi_comm_to_mom, mpi_comm_to_kid, which_kid 
                            

      USE module_timing
      IMPLICIT NONE

      TYPE(domain), POINTER :: grid          
      TYPE(domain), POINTER :: intermediate_grid
      TYPE(domain), POINTER :: ngrid






real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_szj)           :: szj
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_s1z)           :: s1z
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_spz)           :: spz
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_tcs)           :: tcs
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_moist)           :: moist
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_bxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_bxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_bys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_bye
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_btxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_btxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_btys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_btye
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_dfi_moist)           :: dfi_moist
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_bxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_bxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_bys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_bye
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_btxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_btxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_btys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_btye
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_scalar)           :: scalar
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_bxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_bxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_bys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_bye
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_btxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_btxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_btys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_btye
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_dfi_scalar)           :: dfi_scalar
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_bxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_bxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_bys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_bye
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_btxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_btxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_btys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_btye
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%sm32:grid%em32,num_chem)           :: chem
real      ,DIMENSION(grid%sm31:grid%em31,1:grid%levsiz,grid%sm32:grid%em32,num_ozmixm)           :: ozmixm

      INTEGER nlev, msize
      INTEGER i,j,pig,pjg,cm,cn,nig,njg,retval,k
      INTEGER iparstrt,jparstrt,sw
      TYPE (grid_config_rec_type)            :: config_flags
      REAL xv(2000000)
      INTEGER ixv
      INTEGER       ::          cids, cide, cjds, cjde, ckds, ckde,    &
                                cims, cime, cjms, cjme, ckms, ckme,    &
                                cips, cipe, cjps, cjpe, ckps, ckpe
      INTEGER       ::          iids, iide, ijds, ijde, ikds, ikde,    &
                                iims, iime, ijms, ijme, ikms, ikme,    &
                                iips, iipe, ijps, ijpe, ikps, ikpe
      INTEGER       ::          nids, nide, njds, njde, nkds, nkde,    &
                                nims, nime, njms, njme, nkms, nkme,    &
                                nips, nipe, njps, njpe, nkps, nkpe

      INTEGER idim1,idim2,idim3,idim4,idim5,idim6,idim7
      LOGICAL feedback_flag, feedback_flag_v
      INTEGER icoord, jcoord, idim_cd, jdim_cd, pgr
      INTEGER local_comm, ioffset, myproc, nproc, ierr
      INTEGER thisdomain_max_halo_width
      LOGICAL interp_mp

      interp_mp=grid%interp_mp .or. ngrid%interp_mp

      CALL wrf_get_dm_communicator ( local_comm )
      CALL wrf_get_myproc( myproc )
      CALL wrf_get_nproc( nproc )

      CALL get_ijk_from_grid (  grid ,                   &
                                cids, cide, cjds, cjde, ckds, ckde,    &
                                cims, cime, cjms, cjme, ckms, ckme,    &
                                cips, cipe, cjps, cjpe, ckps, ckpe    )
      CALL get_ijk_from_grid (  intermediate_grid ,              &
                                iids, iide, ijds, ijde, ikds, ikde,    &
                                iims, iime, ijms, ijme, ikms, ikme,    &
                                iips, iipe, ijps, ijpe, ikps, ikpe    )
      CALL get_ijk_from_grid (  ngrid ,              &
                                nids, nide, njds, njde, nkds, nkde,    &
                                nims, nime, njms, njme, nkms, nkme,    &
                                nips, nipe, njps, njpe, nkps, nkpe    )

      CALL nl_get_parent_grid_ratio ( ngrid%id, pgr )
      CALL nl_get_i_parent_start ( intermediate_grid%id, iparstrt )
      CALL nl_get_j_parent_start ( intermediate_grid%id, jparstrt )
      CALL nl_get_shw            ( intermediate_grid%id, sw )
      icoord =    iparstrt - sw
      jcoord =    jparstrt - sw
      idim_cd = iide - iids + 1
      jdim_cd = ijde - ijds + 1

      nlev  = ckde - ckds + 1

      CALL get_dm_max_halo_width ( ngrid%id , thisdomain_max_halo_width )

      IF ( grid%active_this_task ) THEN






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
CALL rsl_lite_to_child_info( msize*4                               &
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
ixv=0 
IF ( SIZE(grid%rt_tendf_stoch) .GT. 1 ) THEN 
DO k = 1,config_flags%num_stoch_levels
xv(k+ixv+1-1)= grid%rt_tendf_stoch(pig,pjg,k)
ENDDO
ixv=ixv+config_flags%num_stoch_levels-1+1
ENDIF
IF ( SIZE(grid%rand_pert) .GT. 1 ) THEN 
DO k = 1,config_flags%num_stoch_levels
xv(k+ixv+1-1)= grid%rand_pert(pig,pjg,k)
ENDDO
ixv=ixv+config_flags%num_stoch_levels-1+1
ENDIF
IF ( SIZE(grid%rstoch) .GT. 1 ) THEN 
DO k = 1,config_flags%num_stoch_levels
xv(k+ixv+1-1)= grid%rstoch(pig,pjg,k)
ENDDO
ixv=ixv+config_flags%num_stoch_levels-1+1
ENDIF
IF ( SIZE(grid%pattern_spp_mp) .GT. 1 ) THEN 
DO k = 1,config_flags%num_stoch_levels
xv(k+ixv+1-1)= grid%pattern_spp_mp(pig,pjg,k)
ENDDO
ixv=ixv+config_flags%num_stoch_levels-1+1
ENDIF
IF ( SIZE(grid%pattern_isppt_mp) .GT. 1 ) THEN 
DO k = 1,config_flags%num_stoch_levels
xv(k+ixv+1-1)= grid%pattern_isppt_mp(pig,pjg,k)
ENDDO
ixv=ixv+config_flags%num_stoch_levels-1+1
ENDIF
IF ( SIZE(grid%pattern_isppt_pbl) .GT. 1 ) THEN 
DO k = 1,config_flags%num_stoch_levels
xv(k+ixv+1-1)= grid%pattern_isppt_pbl(pig,pjg,k)
ENDDO
ixv=ixv+config_flags%num_stoch_levels-1+1
ENDIF
IF ( SIZE(grid%pattern_isppt_rad) .GT. 1 ) THEN 
DO k = 1,config_flags%num_stoch_levels
xv(k+ixv+1-1)= grid%pattern_isppt_rad(pig,pjg,k)
ENDDO
ixv=ixv+config_flags%num_stoch_levels-1+1
ENDIF
IF ( SIZE(grid%pattern_isppt_conv) .GT. 1 ) THEN 
DO k = 1,config_flags%num_stoch_levels
xv(k+ixv+1-1)= grid%pattern_isppt_conv(pig,pjg,k)
ENDDO
ixv=ixv+config_flags%num_stoch_levels-1+1
ENDIF
IF ( SIZE(grid%lu_index) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%lu_index(pig,pjg)
ENDIF
IF ( SIZE(grid%lu_mask) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%lu_mask(pig,pjg)
ENDIF
IF ( SIZE(grid%pdyn_smooth) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%pdyn_smooth(pig,pjg)
ENDIF
IF ( SIZE(grid%pdyn_parent) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%pdyn_parent(pig,pjg)
ENDIF
IF ( SIZE(grid%hres_fis) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%hres_fis(pig,pjg)
ENDIF
IF ( SIZE(grid%hres_lnd) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%hres_lnd(pig,pjg)
ENDIF
IF ( SIZE(grid%hres_lu_index) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%hres_lu_index(pig,pjg)
ENDIF
IF ( SIZE(grid%sm) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%sm(pig,pjg)
ENDIF
IF ( SIZE(grid%pd) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%pd(pig,pjg)
ENDIF
IF ( SIZE(grid%fis) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%fis(pig,pjg)
ENDIF
IF ( SIZE(grid%t) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
xv(k+ixv+1-ckds)= grid%t(pig,pjg,k)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%q) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
xv(k+ixv+1-ckds)= grid%q(pig,pjg,k)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%u) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
xv(k+ixv+1-ckds)= grid%u(pig,pjg,k)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%v) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
xv(k+ixv+1-ckds)= grid%v(pig,pjg,k)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%scurx) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%scurx(pig,pjg)
ENDIF
IF ( SIZE(grid%scury) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%scury(pig,pjg)
ENDIF
IF ( SIZE(grid%charn) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%charn(pig,pjg)
ENDIF
IF ( SIZE(grid%msang) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%msang(pig,pjg)
ENDIF
IF ( SIZE(grid%q2) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
xv(k+ixv+1-ckds)= grid%q2(pig,pjg,k)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%cwm) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
xv(k+ixv+1-ckds)= grid%cwm(pig,pjg,k)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%pint) .GT. 1 ) THEN 
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
IF ( SIZE(grid%sst) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)=grid%sst(pig,pjg)
ENDIF
CALL rsl_lite_to_parent_msg(ixv*4,xv)
CALL rsl_lite_to_child_info( msize*4                               &
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

      END IF

      
      IF ( intercomm_active( grid%id ) ) THEN        
        local_comm = mpi_comm_to_kid( which_kid(ngrid%id), grid%id )
        ioffset = nest_task_offsets(ngrid%id)
      ELSE IF ( intercomm_active( ngrid%id ) ) THEN  
        local_comm = mpi_comm_to_mom( ngrid%id )
        ioffset = nest_task_offsets(ngrid%id)
      END IF

      IF ( grid%active_this_task .OR. ngrid%active_this_task ) THEN
        CALL mpi_comm_rank(local_comm,myproc,ierr)
        CALL mpi_comm_size(local_comm,nproc,ierr)
        CALL rsl_lite_bcast_msgs( myproc, nest_pes_x(grid%id)*nest_pes_y(grid%id),         &
                                        nest_pes_x(ngrid%id)*nest_pes_y(ngrid%id),       &
                                        ioffset, local_comm )
      END IF

      RETURN
      END SUBROUTINE force_domain_nmm_part1



   SUBROUTINE force_domain_nmm_part2 ( grid, ngrid, config_flags    &







,szj,s1z,spz,tcs,moist,moist_bxs,moist_bxe,moist_bys,moist_bye,moist_btxs,moist_btxe,moist_btys,moist_btye,dfi_moist, &
dfi_moist_bxs,dfi_moist_bxe,dfi_moist_bys,dfi_moist_bye,dfi_moist_btxs,dfi_moist_btxe,dfi_moist_btys,dfi_moist_btye,scalar, &
scalar_bxs,scalar_bxe,scalar_bys,scalar_bye,scalar_btxs,scalar_btxe,scalar_btys,scalar_btye,dfi_scalar,dfi_scalar_bxs, &
dfi_scalar_bxe,dfi_scalar_bys,dfi_scalar_bye,dfi_scalar_btxs,dfi_scalar_btxe,dfi_scalar_btys,dfi_scalar_btye,chem,ozmixm &


                 )
      USE module_state_description
      USE module_domain, ONLY : domain, get_ijk_from_grid
      USE module_configure, ONLY : grid_config_rec_type
      USE module_dm, ONLY : ntasks, ntasks_x, ntasks_y, itrace, local_communicator, mytask, &
                            ipe_save, jpe_save, ips_save, jps_save, get_dm_max_halo_width,  &
                            nest_pes_x, nest_pes_y 
                             

      USE module_comm_nesting_dm, ONLY : halo_force_down_sub
      use module_comm_dm, only: HALO_NMM_INTERP_INFO_sub
      use module_comm_dm, only: HALO_NMM_FORCE_DOWN_SST_sub
      IMPLICIT NONE

      TYPE(domain), POINTER :: grid          
      TYPE(domain), POINTER :: ngrid,cgrid






real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_szj)           :: szj
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_s1z)           :: s1z
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_spz)           :: spz
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_tcs)           :: tcs
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_moist)           :: moist
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_bxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_bxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_bys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_bye
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_btxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_btxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_btys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_btye
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_dfi_moist)           :: dfi_moist
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_bxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_bxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_bys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_bye
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_btxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_btxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_btys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_btye
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_scalar)           :: scalar
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_bxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_bxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_bys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_bye
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_btxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_btxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_btys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_btye
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_dfi_scalar)           :: dfi_scalar
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_bxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_bxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_bys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_bye
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_btxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_btxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_btys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_btye
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%sm32:grid%em32,num_chem)           :: chem
real      ,DIMENSION(grid%sm31:grid%em31,1:grid%levsiz,grid%sm32:grid%em32,num_ozmixm)           :: ozmixm

      INTEGER nlev, msize
      INTEGER i,j,pig,pjg,cm,cn,nig,njg,retval,k
      TYPE (grid_config_rec_type)            :: config_flags
      REAL xv(2000000)
      INTEGER ixv
      INTEGER       ::          cids, cide, cjds, cjde, ckds, ckde,    &
                                cims, cime, cjms, cjme, ckms, ckme,    &
                                cips, cipe, cjps, cjpe, ckps, ckpe
      INTEGER       ::          nids, nide, njds, njde, nkds, nkde,    &
                                nims, nime, njms, njme, nkms, nkme,    &
                                nips, nipe, njps, njpe, nkps, nkpe
      INTEGER       ::          ids, ide, jds, jde, kds, kde,    &
                                ims, ime, jms, jme, kms, kme,    &
                                ips, ipe, jps, jpe, kps, kpe

      INTEGER idim1,idim2,idim3,idim4,idim5,idim6,idim7
      REAL  dummy_xs, dummy_xe, dummy_ys, dummy_ye
      LOGICAL feedback_flag, feedback_flag_v

      LOGICAL interp_mp
      integer, parameter :: EConst=0, ECopy=1, EExtrap=2 


      interp_mp=grid%interp_mp .or. ngrid%interp_mp

IF ( ngrid%active_this_task ) THEN
      CALL get_ijk_from_grid (  grid ,                   &
                                cids, cide, cjds, cjde, ckds, ckde,    &
                                cims, cime, cjms, cjme, ckms, ckme,    &
                                cips, cipe, cjps, cjpe, ckps, ckpe    )
      CALL get_ijk_from_grid (  ngrid ,              &
                                nids, nide, njds, njde, nkds, nkde,    &
                                nims, nime, njms, njme, nkms, nkme,    &
                                nips, nipe, njps, njpe, nkps, nkpe    )






      cgrid=>grid
      nlev  = ckde - ckds + 1







CALL rsl_lite_from_parent_info(pig,pjg,retval)
DO while ( retval .eq. 1 )
CALL rsl_lite_from_parent_msg(-1,xv) 
ixv=0 
IF ( SIZE(grid%rt_tendf_stoch) .GT. 1 ) THEN 
DO k = 1,config_flags%num_stoch_levels
grid%rt_tendf_stoch(pig,pjg,k) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+config_flags%num_stoch_levels-1+1
ENDIF
IF ( SIZE(grid%rand_pert) .GT. 1 ) THEN 
DO k = 1,config_flags%num_stoch_levels
grid%rand_pert(pig,pjg,k) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+config_flags%num_stoch_levels-1+1
ENDIF
IF ( SIZE(grid%rstoch) .GT. 1 ) THEN 
DO k = 1,config_flags%num_stoch_levels
grid%rstoch(pig,pjg,k) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+config_flags%num_stoch_levels-1+1
ENDIF
IF ( SIZE(grid%pattern_spp_mp) .GT. 1 ) THEN 
DO k = 1,config_flags%num_stoch_levels
grid%pattern_spp_mp(pig,pjg,k) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+config_flags%num_stoch_levels-1+1
ENDIF
IF ( SIZE(grid%pattern_isppt_mp) .GT. 1 ) THEN 
DO k = 1,config_flags%num_stoch_levels
grid%pattern_isppt_mp(pig,pjg,k) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+config_flags%num_stoch_levels-1+1
ENDIF
IF ( SIZE(grid%pattern_isppt_pbl) .GT. 1 ) THEN 
DO k = 1,config_flags%num_stoch_levels
grid%pattern_isppt_pbl(pig,pjg,k) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+config_flags%num_stoch_levels-1+1
ENDIF
IF ( SIZE(grid%pattern_isppt_rad) .GT. 1 ) THEN 
DO k = 1,config_flags%num_stoch_levels
grid%pattern_isppt_rad(pig,pjg,k) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+config_flags%num_stoch_levels-1+1
ENDIF
IF ( SIZE(grid%pattern_isppt_conv) .GT. 1 ) THEN 
DO k = 1,config_flags%num_stoch_levels
grid%pattern_isppt_conv(pig,pjg,k) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+config_flags%num_stoch_levels-1+1
ENDIF
IF ( SIZE(grid%lu_index) .GT. 1 ) THEN 
ixv=ixv+1
grid%lu_index(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%lu_mask) .GT. 1 ) THEN 
ixv=ixv+1
grid%lu_mask(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%pdyn_smooth) .GT. 1 ) THEN 
ixv=ixv+1
grid%pdyn_smooth(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%pdyn_parent) .GT. 1 ) THEN 
ixv=ixv+1
grid%pdyn_parent(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%hres_fis) .GT. 1 ) THEN 
ixv=ixv+1
grid%hres_fis(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%hres_lnd) .GT. 1 ) THEN 
ixv=ixv+1
grid%hres_lnd(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%hres_lu_index) .GT. 1 ) THEN 
ixv=ixv+1
grid%hres_lu_index(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%sm) .GT. 1 ) THEN 
ixv=ixv+1
grid%sm(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%pd) .GT. 1 ) THEN 
ixv=ixv+1
grid%pd(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%fis) .GT. 1 ) THEN 
ixv=ixv+1
grid%fis(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%t) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
grid%t(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%q) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
grid%q(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%u) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
grid%u(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%v) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
grid%v(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%scurx) .GT. 1 ) THEN 
ixv=ixv+1
grid%scurx(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%scury) .GT. 1 ) THEN 
ixv=ixv+1
grid%scury(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%charn) .GT. 1 ) THEN 
ixv=ixv+1
grid%charn(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%msang) .GT. 1 ) THEN 
ixv=ixv+1
grid%msang(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%q2) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
grid%q2(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%cwm) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
grid%cwm(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%pint) .GT. 1 ) THEN 
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
IF ( SIZE(grid%sst) .GT. 1 ) THEN 
ixv=ixv+1
grid%sst(pig,pjg) = xv(ixv)
ENDIF
CALL rsl_lite_from_parent_info(pig,pjg,retval)
ENDDO


      CALL get_ijk_from_grid (  grid ,              &
                                ids, ide, jds, jde, kds, kde,    &
                                ims, ime, jms, jme, kms, kme,    &
                                ips, ipe, jps, jpe, kps, kpe    )

      CALL push_communicators_for_domain( grid%id )

        IF(ngrid%force_sst(1) == 1) then






CALL HALO_NMM_FORCE_DOWN_SST_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )

        END IF







CALL HALO_FORCE_DOWN_sub ( grid, &
  config_flags, &
  num_moist, &
  moist, &
  num_scalar, &
  scalar, &
  num_dfi_scalar, &
  dfi_scalar, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )

      CALL pop_communicators_for_domain

      call store_interp_info(ngrid,grid)


      call ext_c2b_fulldom(ngrid%IIH,ngrid%JJH,ngrid%HBWGT1, &
           ngrid%HBWGT2,ngrid%HBWGT3,ngrid%HBWGT4,         &
           ngrid%deta1,ngrid%deta2,ngrid%eta1,             &
           ngrid%eta2,ngrid%pt,ngrid%pdtop,                &
           grid%pint,grid%t,grid%pd,grid%q,                &
           cims, cime, cjms, cjme, ckms, ckme,             &
           nids, nide, njds, njde, nkds, nkde,             &
           nims, nime, njms, njme, nkms, nkme,             &
           nips, nipe, njps, njpe, nkps, nkpe,             &
           ngrid%iinfo_bxs, ngrid%iinfo_bxe,               &
           ngrid%iinfo_bys, ngrid%iinfo_bye,               &
           ngrid%winfo_bxs, ngrid%winfo_bxe,               &
           ngrid%winfo_bys, ngrid%winfo_bye,               &
           ngrid%pd_bxs, ngrid%pd_bxe,             &
           ngrid%pd_bys, ngrid%pd_bye,             &
           ngrid%t_bxs, ngrid%t_bxe,               &
           ngrid%t_bys, ngrid%t_bye,               &
           ngrid%q_bxs, ngrid%q_bxe,               &
           ngrid%q_bys, ngrid%q_bye)

      
      
      grid=>ngrid
      CALL get_ijk_from_grid (  grid ,              &
                                ids, ide, jds, jde, kds, kde,    &
                                ims, ime, jms, jme, kms, kme,    &
                                ips, ipe, jps, jpe, kps, kpe    )
      CALL push_communicators_for_domain( grid%id )






CALL HALO_NMM_INTERP_INFO_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )

      CALL pop_communicators_for_domain

      grid=>cgrid
      CALL get_ijk_from_grid (  grid ,              &
                                ids, ide, jds, jde, kds, kde,    &
                                ims, ime, jms, jme, kms, kme,    &
                                ips, ipe, jps, jpe, kps, kpe    )


      
      CALL push_communicators_for_domain( grid%id )






IF ( SIZE( grid%rt_tendf_stoch, 1 ) * SIZE( grid%rt_tendf_stoch, 2 ) .GT. 1 ) THEN 
CALL force_sppt_nmm (  &         
                  grid%rt_tendf_stoch,   &       
                 cids, cide, 1, config_flags%num_stoch_levels, cjds, cjde,   &         
                 cims, cime, 1, config_flags%num_stoch_levels, cjms, cjme,   &         
                 cips, cipe, 1, config_flags%num_stoch_levels, cjps, cjpe,   &         
                  ngrid%rt_tendf_stoch,  &   
                 nids, nide, 1, config_flags%num_stoch_levels, njds, njde,   &         
                 nims, nime, 1, config_flags%num_stoch_levels, njms, njme,   &         
                 nips, nipe, 1, config_flags%num_stoch_levels, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
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
                  grid%rand_pert,   &       
                 cids, cide, 1, config_flags%num_stoch_levels, cjds, cjde,   &         
                 cims, cime, 1, config_flags%num_stoch_levels, cjms, cjme,   &         
                 cips, cipe, 1, config_flags%num_stoch_levels, cjps, cjpe,   &         
                  ngrid%rand_pert,  &   
                 nids, nide, 1, config_flags%num_stoch_levels, njds, njde,   &         
                 nims, nime, 1, config_flags%num_stoch_levels, njms, njme,   &         
                 nips, nipe, 1, config_flags%num_stoch_levels, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
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
                  grid%rstoch,   &       
                 cids, cide, 1, config_flags%num_stoch_levels, cjds, cjde,   &         
                 cims, cime, 1, config_flags%num_stoch_levels, cjms, cjme,   &         
                 cips, cipe, 1, config_flags%num_stoch_levels, cjps, cjpe,   &         
                  ngrid%rstoch,  &   
                 nids, nide, 1, config_flags%num_stoch_levels, njds, njde,   &         
                 nims, nime, 1, config_flags%num_stoch_levels, njms, njme,   &         
                 nips, nipe, 1, config_flags%num_stoch_levels, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
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
                  grid%pattern_spp_mp,   &       
                 cids, cide, 1, config_flags%num_stoch_levels, cjds, cjde,   &         
                 cims, cime, 1, config_flags%num_stoch_levels, cjms, cjme,   &         
                 cips, cipe, 1, config_flags%num_stoch_levels, cjps, cjpe,   &         
                  ngrid%pattern_spp_mp,  &   
                 nids, nide, 1, config_flags%num_stoch_levels, njds, njde,   &         
                 nims, nime, 1, config_flags%num_stoch_levels, njms, njme,   &         
                 nips, nipe, 1, config_flags%num_stoch_levels, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
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
                  grid%pattern_isppt_mp,   &       
                 cids, cide, 1, config_flags%num_stoch_levels, cjds, cjde,   &         
                 cims, cime, 1, config_flags%num_stoch_levels, cjms, cjme,   &         
                 cips, cipe, 1, config_flags%num_stoch_levels, cjps, cjpe,   &         
                  ngrid%pattern_isppt_mp,  &   
                 nids, nide, 1, config_flags%num_stoch_levels, njds, njde,   &         
                 nims, nime, 1, config_flags%num_stoch_levels, njms, njme,   &         
                 nips, nipe, 1, config_flags%num_stoch_levels, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
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
                  grid%pattern_isppt_pbl,   &       
                 cids, cide, 1, config_flags%num_stoch_levels, cjds, cjde,   &         
                 cims, cime, 1, config_flags%num_stoch_levels, cjms, cjme,   &         
                 cips, cipe, 1, config_flags%num_stoch_levels, cjps, cjpe,   &         
                  ngrid%pattern_isppt_pbl,  &   
                 nids, nide, 1, config_flags%num_stoch_levels, njds, njde,   &         
                 nims, nime, 1, config_flags%num_stoch_levels, njms, njme,   &         
                 nips, nipe, 1, config_flags%num_stoch_levels, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
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
                  grid%pattern_isppt_rad,   &       
                 cids, cide, 1, config_flags%num_stoch_levels, cjds, cjde,   &         
                 cims, cime, 1, config_flags%num_stoch_levels, cjms, cjme,   &         
                 cips, cipe, 1, config_flags%num_stoch_levels, cjps, cjpe,   &         
                  ngrid%pattern_isppt_rad,  &   
                 nids, nide, 1, config_flags%num_stoch_levels, njds, njde,   &         
                 nims, nime, 1, config_flags%num_stoch_levels, njms, njme,   &         
                 nips, nipe, 1, config_flags%num_stoch_levels, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
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
                  grid%pattern_isppt_conv,   &       
                 cids, cide, 1, config_flags%num_stoch_levels, cjds, cjde,   &         
                 cims, cime, 1, config_flags%num_stoch_levels, cjms, cjme,   &         
                 cips, cipe, 1, config_flags%num_stoch_levels, cjps, cjpe,   &         
                  ngrid%pattern_isppt_conv,  &   
                 nids, nide, 1, config_flags%num_stoch_levels, njds, njde,   &         
                 nims, nime, 1, config_flags%num_stoch_levels, njms, njme,   &         
                 nips, nipe, 1, config_flags%num_stoch_levels, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
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
CONTINUE 
ENDIF
IF ( SIZE( grid%lu_mask, 1 ) * SIZE( grid%lu_mask, 2 ) .GT. 1 ) THEN 
CONTINUE 
ENDIF
IF ( SIZE( grid%pdyn_smooth, 1 ) * SIZE( grid%pdyn_smooth, 2 ) .GT. 1 ) THEN 
CONTINUE 
ENDIF
IF ( SIZE( grid%pdyn_parent, 1 ) * SIZE( grid%pdyn_parent, 2 ) .GT. 1 ) THEN 
CALL downaged2d (  &         
                  grid%pdyn_parent,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%pdyn_parent,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,grid%pdyn_smooth_age&
,ngrid%pdyn_parent_age&
,grid%pdyn_smooth&
                  ) 
ENDIF
IF ( SIZE( grid%hres_fis, 1 ) * SIZE( grid%hres_fis, 2 ) .GT. 1 ) THEN 
CONTINUE 
ENDIF
IF ( SIZE( grid%hres_lnd, 1 ) * SIZE( grid%hres_lnd, 2 ) .GT. 1 ) THEN 
CONTINUE 
ENDIF
IF ( SIZE( grid%hres_lu_index, 1 ) * SIZE( grid%hres_lu_index, 2 ) .GT. 1 ) THEN 
CONTINUE 
ENDIF
IF ( SIZE( grid%sm, 1 ) * SIZE( grid%sm, 2 ) .GT. 1 ) THEN 
CONTINUE 
ENDIF
IF ( SIZE( grid%pd, 1 ) * SIZE( grid%pd, 2 ) .GT. 1 ) THEN 
CONTINUE 
ENDIF
IF ( SIZE( grid%fis, 1 ) * SIZE( grid%fis, 2 ) .GT. 1 ) THEN 
CONTINUE 
ENDIF
IF ( SIZE( grid%t, 1 ) * SIZE( grid%t, 2 ) .GT. 1 ) THEN 
CONTINUE 
ENDIF
IF ( SIZE( grid%q, 1 ) * SIZE( grid%q, 2 ) .GT. 1 ) THEN 
CONTINUE 
ENDIF
IF ( SIZE( grid%u, 1 ) * SIZE( grid%u, 2 ) .GT. 1 ) THEN 
CALL bdyvel (  &         
                  grid%u,   &       
                 cids, cide, ckds, ckde, cjds, cjde,   &         
                 cims, cime, ckms, ckme, cjms, cjme,   &         
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         
                  ngrid%u,  &   
                 nids, nide, nkds, nkde, njds, njde,   &         
                 nims, nime, nkms, nkme, njms, njme,   &         
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
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
                  grid%v,   &       
                 cids, cide, ckds, ckde, cjds, cjde,   &         
                 cims, cime, ckms, ckme, cjms, cjme,   &         
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         
                  ngrid%v,  &   
                 nids, nide, nkds, nkde, njds, njde,   &         
                 nims, nime, nkms, nkme, njms, njme,   &         
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
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
                  grid%scurx,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%scurx,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
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
                  grid%scury,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%scury,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
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
                  grid%charn,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%charn,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
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
                  grid%msang,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%msang,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
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
                  grid%q2,   &       
                 cids, cide, ckds, ckde, cjds, cjde,   &         
                 cims, cime, ckms, ckme, cjms, cjme,   &         
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         
                  ngrid%q2,  &   
                 nids, nide, nkds, nkde, njds, njde,   &         
                 nims, nime, nkms, nkme, njms, njme,   &         
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
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
                  grid%cwm,   &       
                 cids, cide, ckds, ckde, cjds, cjde,   &         
                 cims, cime, ckms, ckme, cjms, cjme,   &         
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         
                  ngrid%cwm,  &   
                 nids, nide, nkds, nkde, njds, njde,   &         
                 nims, nime, nkms, nkme, njms, njme,   &         
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
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
CONTINUE 
ENDIF
DO itrace = PARAM_FIRST_SCALAR, num_moist
IF ( SIZE( moist, 1 ) * SIZE( moist, 2 ) .GT. 1 .and. (interp_mp .eqv. .true.) ) THEN 
CALL bdymass (  &         
                  moist(grid%sm31,grid%sm32,grid%sm33,itrace),   &       
                 cids, cide, ckds, ckde, cjds, cjde,   &         
                 cims, cime, ckms, ckme, cjms, cjme,   &         
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         
                  ngrid%moist(ngrid%sm31,ngrid%sm32,ngrid%sm33,itrace),  &   
                 nids, nide, nkds, nkde, njds, njde,   &         
                 nims, nime, nkms, nkme, njms, njme,   &         
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
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
                  scalar(grid%sm31,grid%sm32,grid%sm33,itrace),   &       
                 cids, cide, ckds, ckde, cjds, cjde,   &         
                 cims, cime, ckms, ckme, cjms, cjme,   &         
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         
                  ngrid%scalar(ngrid%sm31,ngrid%sm32,ngrid%sm33,itrace),  &   
                 nids, nide, nkds, nkde, njds, njde,   &         
                 nims, nime, nkms, nkme, njms, njme,   &         
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
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
                  dfi_scalar(grid%sm31,grid%sm32,grid%sm33,itrace),   &       
                 cids, cide, ckds, ckde, cjds, cjde,   &         
                 cims, cime, ckms, ckme, cjms, cjme,   &         
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         
                  ngrid%dfi_scalar(ngrid%sm31,ngrid%sm32,ngrid%sm33,itrace),  &   
                 nids, nide, nkds, nkde, njds, njde,   &         
                 nims, nime, nkms, nkme, njms, njme,   &         
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
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
                  grid%sst,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%sst,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
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

      CALL pop_communicators_for_domain

END IF

      RETURN
   END SUBROUTINE force_domain_nmm_part2












   SUBROUTINE feedback_nest_prep_nmm ( grid, config_flags    &







,szj,s1z,spz,tcs,moist,moist_bxs,moist_bxe,moist_bys,moist_bye,moist_btxs,moist_btxe,moist_btys,moist_btye,dfi_moist, &
dfi_moist_bxs,dfi_moist_bxe,dfi_moist_bys,dfi_moist_bye,dfi_moist_btxs,dfi_moist_btxe,dfi_moist_btys,dfi_moist_btye,scalar, &
scalar_bxs,scalar_bxe,scalar_bys,scalar_bye,scalar_btxs,scalar_btxe,scalar_btys,scalar_btye,dfi_scalar,dfi_scalar_bxs, &
dfi_scalar_bxe,dfi_scalar_bys,dfi_scalar_bye,dfi_scalar_btxs,dfi_scalar_btxe,dfi_scalar_btys,dfi_scalar_btye,chem,ozmixm &


)
      USE module_state_description
      USE module_domain, ONLY : domain, get_ijk_from_grid
      USE module_configure, ONLY : grid_config_rec_type
      USE module_dm, ONLY : ntasks, ntasks_x, ntasks_y, itrace, local_communicator, mytask, &
                            ipe_save, jpe_save, ips_save, jps_save, get_dm_max_halo_width,  &
                            nest_pes_x, nest_pes_y
                            
      USE module_comm_dm, ONLY : HALO_NMM_WEIGHTS_sub
      USE module_comm_nesting_dm, ONLY : HALO_INTERP_UP_sub

      IMPLICIT NONE

      TYPE(domain), TARGET :: grid          
      TYPE (grid_config_rec_type) :: config_flags 
                                                  
                                                  






real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_szj)           :: szj
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_s1z)           :: s1z
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_spz)           :: spz
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_tcs)           :: tcs
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_moist)           :: moist
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_bxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_bxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_bys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_bye
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_btxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_btxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_btys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_btye
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_dfi_moist)           :: dfi_moist
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_bxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_bxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_bys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_bye
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_btxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_btxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_btys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_btye
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_scalar)           :: scalar
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_bxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_bxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_bys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_bye
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_btxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_btxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_btys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_btye
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_dfi_scalar)           :: dfi_scalar
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_bxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_bxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_bys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_bye
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_btxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_btxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_btys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_btye
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%sm32:grid%em32,num_chem)           :: chem
real      ,DIMENSION(grid%sm31:grid%em31,1:grid%levsiz,grid%sm32:grid%em32,num_ozmixm)           :: ozmixm


      INTEGER       ::          ids, ide, jds, jde, kds, kde,    &
                                ims, ime, jms, jme, kms, kme,    &
                                ips, ipe, jps, jpe, kps, kpe

      INTEGER idim1,idim2,idim3,idim4,idim5,idim6,idim7

      INTEGER       :: idum1, idum2
      LOGICAL :: interp_mp
      interp_mp=.true.



      CALL get_ijk_from_grid (  grid ,              &
                                ids, ide, jds, jde, kds, kde,    &
                                ims, ime, jms, jme, kms, kme,    &
                                ips, ipe, jps, jpe, kps, kpe    )
    IF ( grid%active_this_task ) THEN
      CALL push_communicators_for_domain( grid%id )






CALL HALO_INTERP_UP_sub ( grid, &
  config_flags, &
  num_moist, &
  moist, &
  num_scalar, &
  scalar, &
  num_dfi_scalar, &
  dfi_scalar, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )







CALL HALO_NMM_WEIGHTS_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )

      CALL pop_communicators_for_domain
    END IF

   END SUBROUTINE feedback_nest_prep_nmm






   SUBROUTINE force_intermediate_nmm ( grid, ngrid, config_flags    &







,szj,s1z,spz,tcs,moist,moist_bxs,moist_bxe,moist_bys,moist_bye,moist_btxs,moist_btxe,moist_btys,moist_btye,dfi_moist, &
dfi_moist_bxs,dfi_moist_bxe,dfi_moist_bys,dfi_moist_bye,dfi_moist_btxs,dfi_moist_btxe,dfi_moist_btys,dfi_moist_btye,scalar, &
scalar_bxs,scalar_bxe,scalar_bys,scalar_bye,scalar_btxs,scalar_btxe,scalar_btys,scalar_btye,dfi_scalar,dfi_scalar_bxs, &
dfi_scalar_bxe,dfi_scalar_bys,dfi_scalar_bye,dfi_scalar_btxs,dfi_scalar_btxe,dfi_scalar_btys,dfi_scalar_btye,chem,ozmixm &


                 )
      USE module_state_description
      USE module_domain, ONLY : domain, get_ijk_from_grid
      USE module_configure, ONLY : grid_config_rec_type
      USE module_dm, ONLY : ntasks, ntasks_x, ntasks_y, itrace, local_communicator, mytask, &
                            ipe_save, jpe_save, ips_save, jps_save, get_dm_max_halo_width
      USE module_comm_nesting_dm, ONLY : halo_force_down_sub
      IMPLICIT NONE

      TYPE(domain), POINTER :: grid          
      TYPE(domain), POINTER :: cgrid
      TYPE(domain), POINTER :: ngrid






real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_szj)           :: szj
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_s1z)           :: s1z
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_spz)           :: spz
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_tcs)           :: tcs
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_moist)           :: moist
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_bxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_bxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_bys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_bye
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_btxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_btxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_btys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_btye
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_dfi_moist)           :: dfi_moist
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_bxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_bxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_bys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_bye
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_btxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_btxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_btys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_btye
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_scalar)           :: scalar
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_bxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_bxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_bys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_bye
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_btxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_btxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_btys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_btye
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_dfi_scalar)           :: dfi_scalar
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_bxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_bxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_bys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_bye
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_btxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_btxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_btys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_btye
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%sm32:grid%em32,num_chem)           :: chem
real      ,DIMENSION(grid%sm31:grid%em31,1:grid%levsiz,grid%sm32:grid%em32,num_ozmixm)           :: ozmixm

      INTEGER nlev, msize
      INTEGER i,j,pig,pjg,cm,cn,nig,njg,retval,k
      TYPE (grid_config_rec_type)            :: config_flags
      REAL xv(2000000)
      INTEGER ixv
      INTEGER       ::          cids, cide, cjds, cjde, ckds, ckde,    &
                                cims, cime, cjms, cjme, ckms, ckme,    &
                                cips, cipe, cjps, cjpe, ckps, ckpe
      INTEGER       ::          nids, nide, njds, njde, nkds, nkde,    &
                                nims, nime, njms, njme, nkms, nkme,    &
                                nips, nipe, njps, njpe, nkps, nkpe
      INTEGER       ::          ids, ide, jds, jde, kds, kde,    &
                                ims, ime, jms, jme, kms, kme,    &
                                ips, ipe, jps, jpe, kps, kpe

      INTEGER idim1,idim2,idim3,idim4,idim5,idim6,idim7
      REAL  dummy_xs, dummy_xe, dummy_ys, dummy_ye
      LOGICAL feedback_flag, feedback_flag_v

integer myproc
      LOGICAL interp_mp









      interp_mp=grid%interp_mp .or. ngrid%interp_mp


      CALL get_ijk_from_grid (  grid ,                   &
                                cids, cide, cjds, cjde, ckds, ckde,    &
                                cims, cime, cjms, cjme, ckms, ckme,    &
                                cips, cipe, cjps, cjpe, ckps, ckpe    )
      CALL get_ijk_from_grid (  ngrid ,              &
                                nids, nide, njds, njde, nkds, nkde,    &
                                nims, nime, njms, njme, nkms, nkme,    &
                                nips, nipe, njps, njpe, nkps, nkpe    )

      cgrid=>grid
      nlev  = ckde - ckds + 1







CALL rsl_lite_from_parent_info(pig,pjg,retval)
DO while ( retval .eq. 1 )
CALL rsl_lite_from_parent_msg(-1,xv) 
ixv=0 
IF ( SIZE(grid%lakedepth2d) .GT. 1 ) THEN 
ixv=ixv+1
grid%lakedepth2d(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%savedtke12d) .GT. 1 ) THEN 
ixv=ixv+1
grid%savedtke12d(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%snowdp2d) .GT. 1 ) THEN 
ixv=ixv+1
grid%snowdp2d(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%h2osno2d) .GT. 1 ) THEN 
ixv=ixv+1
grid%h2osno2d(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%snl2d) .GT. 1 ) THEN 
ixv=ixv+1
grid%snl2d(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%t_grnd2d) .GT. 1 ) THEN 
ixv=ixv+1
grid%t_grnd2d(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%t_lake3d) .GT. 1 ) THEN 
DO k = 1,10
grid%t_lake3d(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%lake_icefrac3d) .GT. 1 ) THEN 
DO k = 1,10
grid%lake_icefrac3d(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%z_lake3d) .GT. 1 ) THEN 
DO k = 1,10
grid%z_lake3d(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%dz_lake3d) .GT. 1 ) THEN 
DO k = 1,10
grid%dz_lake3d(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%t_soisno3d) .GT. 1 ) THEN 
DO k = 1,15
grid%t_soisno3d(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+15-1+1
ENDIF
IF ( SIZE(grid%h2osoi_ice3d) .GT. 1 ) THEN 
DO k = 1,15
grid%h2osoi_ice3d(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+15-1+1
ENDIF
IF ( SIZE(grid%h2osoi_liq3d) .GT. 1 ) THEN 
DO k = 1,15
grid%h2osoi_liq3d(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+15-1+1
ENDIF
IF ( SIZE(grid%h2osoi_vol3d) .GT. 1 ) THEN 
DO k = 1,15
grid%h2osoi_vol3d(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+15-1+1
ENDIF
IF ( SIZE(grid%z3d) .GT. 1 ) THEN 
DO k = 1,15
grid%z3d(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+15-1+1
ENDIF
IF ( SIZE(grid%dz3d) .GT. 1 ) THEN 
DO k = 1,15
grid%dz3d(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+15-1+1
ENDIF
IF ( SIZE(grid%zi3d) .GT. 1 ) THEN 
DO k = 1,16
grid%zi3d(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+16-1+1
ENDIF
IF ( SIZE(grid%watsat3d) .GT. 1 ) THEN 
DO k = 1,10
grid%watsat3d(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%csol3d) .GT. 1 ) THEN 
DO k = 1,10
grid%csol3d(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%tkmg3d) .GT. 1 ) THEN 
DO k = 1,10
grid%tkmg3d(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%tkdry3d) .GT. 1 ) THEN 
DO k = 1,10
grid%tkdry3d(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%tksatu3d) .GT. 1 ) THEN 
DO k = 1,10
grid%tksatu3d(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%vertstrucc) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
grid%vertstrucc(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%vertstrucs) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
grid%vertstrucs(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%ru_tendf_stoch) .GT. 1 ) THEN 
DO k = 1,config_flags%num_stoch_levels
grid%ru_tendf_stoch(pig,pjg,k) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+config_flags%num_stoch_levels-1+1
ENDIF
IF ( SIZE(grid%rv_tendf_stoch) .GT. 1 ) THEN 
DO k = 1,config_flags%num_stoch_levels
grid%rv_tendf_stoch(pig,pjg,k) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+config_flags%num_stoch_levels-1+1
ENDIF
IF ( SIZE(grid%lu_index) .GT. 1 ) THEN 
ixv=ixv+1
grid%lu_index(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%lu_mask) .GT. 1 ) THEN 
ixv=ixv+1
grid%lu_mask(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%cuprecip_swath) .GT. 1 ) THEN 
ixv=ixv+1
grid%cuprecip_swath(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%precip_swath) .GT. 1 ) THEN 
ixv=ixv+1
grid%precip_swath(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%windsq_swath) .GT. 1 ) THEN 
ixv=ixv+1
grid%windsq_swath(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%suswind) .GT. 1 ) THEN 
ixv=ixv+1
grid%suswind(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%wind10_ratio) .GT. 1 ) THEN 
ixv=ixv+1
grid%wind10_ratio(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%suswind_swath) .GT. 1 ) THEN 
ixv=ixv+1
grid%suswind_swath(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%membrane_mslp) .GT. 1 ) THEN 
ixv=ixv+1
grid%membrane_mslp(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%pdyn_smooth) .GT. 1 ) THEN 
ixv=ixv+1
grid%pdyn_smooth(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%pdyn_parent) .GT. 1 ) THEN 
ixv=ixv+1
grid%pdyn_parent(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%best_mslp) .GT. 1 ) THEN 
ixv=ixv+1
grid%best_mslp(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%hlon) .GT. 1 ) THEN 
ixv=ixv+1
grid%hlon(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%hlat) .GT. 1 ) THEN 
ixv=ixv+1
grid%hlat(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tg_max_m10wind) .GT. 1 ) THEN 
ixv=ixv+1
grid%tg_max_m10wind(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tg_max_wwind) .GT. 1 ) THEN 
ixv=ixv+1
grid%tg_max_wwind(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tg_min_wwind) .GT. 1 ) THEN 
ixv=ixv+1
grid%tg_min_wwind(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tg_max_zhel_25) .GT. 1 ) THEN 
ixv=ixv+1
grid%tg_max_zhel_25(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tg_min_zhel_25) .GT. 1 ) THEN 
ixv=ixv+1
grid%tg_min_zhel_25(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tg_max_zhel_03) .GT. 1 ) THEN 
ixv=ixv+1
grid%tg_max_zhel_03(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tg_min_zhel_03) .GT. 1 ) THEN 
ixv=ixv+1
grid%tg_min_zhel_03(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tg_updhel25) .GT. 1 ) THEN 
ixv=ixv+1
grid%tg_updhel25(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tg_max_updhel25) .GT. 1 ) THEN 
ixv=ixv+1
grid%tg_max_updhel25(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tg_updhel03) .GT. 1 ) THEN 
ixv=ixv+1
grid%tg_updhel03(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tg_max_updhel03) .GT. 1 ) THEN 
ixv=ixv+1
grid%tg_max_updhel03(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tg_total_precip) .GT. 1 ) THEN 
ixv=ixv+1
grid%tg_total_precip(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%hres_fis) .GT. 1 ) THEN 
ixv=ixv+1
grid%hres_fis(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%hres_lnd) .GT. 1 ) THEN 
ixv=ixv+1
grid%hres_lnd(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%hres_lu_index) .GT. 1 ) THEN 
ixv=ixv+1
grid%hres_lu_index(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%sm) .GT. 1 ) THEN 
ixv=ixv+1
grid%sm(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%sice) .GT. 1 ) THEN 
ixv=ixv+1
grid%sice(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%pd) .GT. 1 ) THEN 
ixv=ixv+1
grid%pd(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%fis) .GT. 1 ) THEN 
ixv=ixv+1
grid%fis(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%t) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
grid%t(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%q) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
grid%q(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%u) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
grid%u(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%v) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
grid%v(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%ustar) .GT. 1 ) THEN 
ixv=ixv+1
grid%ustar(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%z0) .GT. 1 ) THEN 
ixv=ixv+1
grid%z0(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%scurx) .GT. 1 ) THEN 
ixv=ixv+1
grid%scurx(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%scury) .GT. 1 ) THEN 
ixv=ixv+1
grid%scury(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%charn) .GT. 1 ) THEN 
ixv=ixv+1
grid%charn(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%msang) .GT. 1 ) THEN 
ixv=ixv+1
grid%msang(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%ths) .GT. 1 ) THEN 
ixv=ixv+1
grid%ths(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%qsh) .GT. 1 ) THEN 
ixv=ixv+1
grid%qsh(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%taux) .GT. 1 ) THEN 
ixv=ixv+1
grid%taux(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tauy) .GT. 1 ) THEN 
ixv=ixv+1
grid%tauy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%acprec) .GT. 1 ) THEN 
ixv=ixv+1
grid%acprec(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%cuprec) .GT. 1 ) THEN 
ixv=ixv+1
grid%cuprec(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%cldefi) .GT. 1 ) THEN 
ixv=ixv+1
grid%cldefi(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%th10) .GT. 1 ) THEN 
ixv=ixv+1
grid%th10(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%q10) .GT. 1 ) THEN 
ixv=ixv+1
grid%q10(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%pshltr) .GT. 1 ) THEN 
ixv=ixv+1
grid%pshltr(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tshltr) .GT. 1 ) THEN 
ixv=ixv+1
grid%tshltr(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%qshltr) .GT. 1 ) THEN 
ixv=ixv+1
grid%qshltr(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%q2) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
grid%q2(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%t_adj) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
grid%t_adj(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%cd_out) .GT. 1 ) THEN 
ixv=ixv+1
grid%cd_out(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%ch_out) .GT. 1 ) THEN 
ixv=ixv+1
grid%ch_out(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%albase) .GT. 1 ) THEN 
ixv=ixv+1
grid%albase(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%nmm_tsk) .GT. 1 ) THEN 
ixv=ixv+1
grid%nmm_tsk(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%mxsnal) .GT. 1 ) THEN 
ixv=ixv+1
grid%mxsnal(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%sigt4) .GT. 1 ) THEN 
ixv=ixv+1
grid%sigt4(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tg) .GT. 1 ) THEN 
ixv=ixv+1
grid%tg(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%acswupt) .GT. 1 ) THEN 
ixv=ixv+1
grid%acswupt(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%acswuptc) .GT. 1 ) THEN 
ixv=ixv+1
grid%acswuptc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%acswdnt) .GT. 1 ) THEN 
ixv=ixv+1
grid%acswdnt(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%acswdntc) .GT. 1 ) THEN 
ixv=ixv+1
grid%acswdntc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%acswupb) .GT. 1 ) THEN 
ixv=ixv+1
grid%acswupb(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%acswupbc) .GT. 1 ) THEN 
ixv=ixv+1
grid%acswupbc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%acswdnb) .GT. 1 ) THEN 
ixv=ixv+1
grid%acswdnb(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%acswdnbc) .GT. 1 ) THEN 
ixv=ixv+1
grid%acswdnbc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%aclwupt) .GT. 1 ) THEN 
ixv=ixv+1
grid%aclwupt(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%aclwuptc) .GT. 1 ) THEN 
ixv=ixv+1
grid%aclwuptc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%aclwdnt) .GT. 1 ) THEN 
ixv=ixv+1
grid%aclwdnt(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%aclwdntc) .GT. 1 ) THEN 
ixv=ixv+1
grid%aclwdntc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%aclwupb) .GT. 1 ) THEN 
ixv=ixv+1
grid%aclwupb(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%aclwupbc) .GT. 1 ) THEN 
ixv=ixv+1
grid%aclwupbc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%aclwdnb) .GT. 1 ) THEN 
ixv=ixv+1
grid%aclwdnb(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%aclwdnbc) .GT. 1 ) THEN 
ixv=ixv+1
grid%aclwdnbc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%swupt) .GT. 1 ) THEN 
ixv=ixv+1
grid%swupt(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%swuptc) .GT. 1 ) THEN 
ixv=ixv+1
grid%swuptc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%swdnt) .GT. 1 ) THEN 
ixv=ixv+1
grid%swdnt(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%swdntc) .GT. 1 ) THEN 
ixv=ixv+1
grid%swdntc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%swupb) .GT. 1 ) THEN 
ixv=ixv+1
grid%swupb(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%swupbc) .GT. 1 ) THEN 
ixv=ixv+1
grid%swupbc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%swdnb) .GT. 1 ) THEN 
ixv=ixv+1
grid%swdnb(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%swdnbc) .GT. 1 ) THEN 
ixv=ixv+1
grid%swdnbc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%lwupt) .GT. 1 ) THEN 
ixv=ixv+1
grid%lwupt(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%lwuptc) .GT. 1 ) THEN 
ixv=ixv+1
grid%lwuptc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%lwdnt) .GT. 1 ) THEN 
ixv=ixv+1
grid%lwdnt(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%lwdntc) .GT. 1 ) THEN 
ixv=ixv+1
grid%lwdntc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%lwupb) .GT. 1 ) THEN 
ixv=ixv+1
grid%lwupb(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%lwupbc) .GT. 1 ) THEN 
ixv=ixv+1
grid%lwupbc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%lwdnb) .GT. 1 ) THEN 
ixv=ixv+1
grid%lwdnb(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%lwdnbc) .GT. 1 ) THEN 
ixv=ixv+1
grid%lwdnbc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%refl_10cm) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
grid%refl_10cm(pig,k,pjg) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%refd_max) .GT. 1 ) THEN 
ixv=ixv+1
grid%refd_max(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%qnwfa2d) .GT. 1 ) THEN 
ixv=ixv+1
grid%qnwfa2d(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%swddir) .GT. 1 ) THEN 
ixv=ixv+1
grid%swddir(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%swddirc) .GT. 1 ) THEN 
ixv=ixv+1
grid%swddirc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%swddni) .GT. 1 ) THEN 
ixv=ixv+1
grid%swddni(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%swddnic) .GT. 1 ) THEN 
ixv=ixv+1
grid%swddnic(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%swddif) .GT. 1 ) THEN 
ixv=ixv+1
grid%swddif(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%gx) .GT. 1 ) THEN 
ixv=ixv+1
grid%gx(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%bx) .GT. 1 ) THEN 
ixv=ixv+1
grid%bx(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%gg) .GT. 1 ) THEN 
ixv=ixv+1
grid%gg(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%bb) .GT. 1 ) THEN 
ixv=ixv+1
grid%bb(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%coszen_ref) .GT. 1 ) THEN 
ixv=ixv+1
grid%coszen_ref(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%swdown_ref) .GT. 1 ) THEN 
ixv=ixv+1
grid%swdown_ref(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%swddir_ref) .GT. 1 ) THEN 
ixv=ixv+1
grid%swddir_ref(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%cwm) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
grid%cwm(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%f_ice) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
grid%f_ice(pig,k,pjg) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%f_rain) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
grid%f_rain(pig,k,pjg) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%f_rimef) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
grid%f_rimef(pig,k,pjg) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%cfrach) .GT. 1 ) THEN 
ixv=ixv+1
grid%cfrach(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%cfracl) .GT. 1 ) THEN 
ixv=ixv+1
grid%cfracl(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%cfracm) .GT. 1 ) THEN 
ixv=ixv+1
grid%cfracm(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%islope) .GT. 1 ) THEN 
ixv=ixv+1
grid%islope(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%cmc) .GT. 1 ) THEN 
ixv=ixv+1
grid%cmc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%soiltb) .GT. 1 ) THEN 
ixv=ixv+1
grid%soiltb(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%vegfrc) .GT. 1 ) THEN 
ixv=ixv+1
grid%vegfrc(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%shdmax) .GT. 1 ) THEN 
ixv=ixv+1
grid%shdmax(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%shdmin) .GT. 1 ) THEN 
ixv=ixv+1
grid%shdmin(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%sh2o) .GT. 1 ) THEN 
DO k = 1,config_flags%num_soil_layers
grid%sh2o(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+config_flags%num_soil_layers-1+1
ENDIF
IF ( SIZE(grid%smc) .GT. 1 ) THEN 
DO k = 1,config_flags%num_soil_layers
grid%smc(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+config_flags%num_soil_layers-1+1
ENDIF
IF ( SIZE(grid%stc) .GT. 1 ) THEN 
DO k = 1,config_flags%num_soil_layers
grid%stc(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+config_flags%num_soil_layers-1+1
ENDIF
IF ( SIZE(grid%ctopo) .GT. 1 ) THEN 
ixv=ixv+1
grid%ctopo(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%ctopo2) .GT. 1 ) THEN 
ixv=ixv+1
grid%ctopo2(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%dwdt) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
grid%dwdt(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%pint) .GT. 1 ) THEN 
DO k = ckds,ckde
grid%pint(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+ckde-ckds+1
ENDIF
IF ( SIZE(grid%w) .GT. 1 ) THEN 
DO k = ckds,ckde
grid%w(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+ckde-ckds+1
ENDIF
IF ( SIZE(grid%w_tot) .GT. 1 ) THEN 
DO k = ckds,ckde
grid%w_tot(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+ckde-ckds+1
ENDIF
IF ( SIZE(grid%z) .GT. 1 ) THEN 
DO k = ckds,ckde
grid%z(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+ckde-ckds+1
ENDIF
IF ( SIZE(grid%rlwin) .GT. 1 ) THEN 
ixv=ixv+1
grid%rlwin(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%rswin) .GT. 1 ) THEN 
ixv=ixv+1
grid%rswin(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%rlwtt) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
grid%rlwtt(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%rswtt) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
grid%rswtt(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%winfo) .GT. 1 ) THEN 
DO k = ckds,ckde
grid%winfo(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+ckde-ckds+1
ENDIF
IF ( SIZE(grid%iinfo) .GT. 1 ) THEN 
DO k = ckds,ckde
grid%iinfo(pig,pjg,k) = xv(k-ckds+ixv+1)
ENDDO
ixv=ixv+ckde-ckds+1
ENDIF
IF ( SIZE(grid%landmask) .GT. 1 ) THEN 
ixv=ixv+1
grid%landmask(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%toposoil) .GT. 1 ) THEN 
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
IF ( SIZE(grid%lake_depth) .GT. 1 ) THEN 
ixv=ixv+1
grid%lake_depth(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%u10) .GT. 1 ) THEN 
ixv=ixv+1
grid%u10(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%v10) .GT. 1 ) THEN 
ixv=ixv+1
grid%v10(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%xice) .GT. 1 ) THEN 
ixv=ixv+1
grid%xice(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%icedepth) .GT. 1 ) THEN 
ixv=ixv+1
grid%icedepth(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%albsi) .GT. 1 ) THEN 
ixv=ixv+1
grid%albsi(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%snowsi) .GT. 1 ) THEN 
ixv=ixv+1
grid%snowsi(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%ivgtyp) .GT. 1 ) THEN 
ixv=ixv+1
grid%ivgtyp(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%isltyp) .GT. 1 ) THEN 
ixv=ixv+1
grid%isltyp(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%vegfra) .GT. 1 ) THEN 
ixv=ixv+1
grid%vegfra(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%sst) .GT. 1 ) THEN 
ixv=ixv+1
grid%sst(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%weasd) .GT. 1 ) THEN 
ixv=ixv+1
grid%weasd(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%thz0) .GT. 1 ) THEN 
ixv=ixv+1
grid%thz0(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%qz0) .GT. 1 ) THEN 
ixv=ixv+1
grid%qz0(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%uz0) .GT. 1 ) THEN 
ixv=ixv+1
grid%uz0(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%vz0) .GT. 1 ) THEN 
ixv=ixv+1
grid%vz0(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%htop) .GT. 1 ) THEN 
ixv=ixv+1
grid%htop(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%hbot) .GT. 1 ) THEN 
ixv=ixv+1
grid%hbot(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%htopr) .GT. 1 ) THEN 
ixv=ixv+1
grid%htopr(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%hbotr) .GT. 1 ) THEN 
ixv=ixv+1
grid%hbotr(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%cuppt) .GT. 1 ) THEN 
ixv=ixv+1
grid%cuppt(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%snowh) .GT. 1 ) THEN 
ixv=ixv+1
grid%snowh(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%rhosn) .GT. 1 ) THEN 
ixv=ixv+1
grid%rhosn(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%isnowxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%isnowxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tvxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%tvxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tgxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%tgxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%canicexy) .GT. 1 ) THEN 
ixv=ixv+1
grid%canicexy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%canliqxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%canliqxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%eahxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%eahxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tahxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%tahxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%cmxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%cmxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%chxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%chxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%fwetxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%fwetxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%sneqvoxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%sneqvoxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%alboldxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%alboldxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%qsnowxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%qsnowxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%wslakexy) .GT. 1 ) THEN 
ixv=ixv+1
grid%wslakexy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%zwtxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%zwtxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%waxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%waxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%wtxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%wtxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tsnoxy) .GT. 1 ) THEN 
DO k = 1,config_flags%num_snow_layers
grid%tsnoxy(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+config_flags%num_snow_layers-1+1
ENDIF
IF ( SIZE(grid%zsnsoxy) .GT. 1 ) THEN 
DO k = 1,config_flags%num_snso_layers
grid%zsnsoxy(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+config_flags%num_snso_layers-1+1
ENDIF
IF ( SIZE(grid%snicexy) .GT. 1 ) THEN 
DO k = 1,config_flags%num_snow_layers
grid%snicexy(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+config_flags%num_snow_layers-1+1
ENDIF
IF ( SIZE(grid%snliqxy) .GT. 1 ) THEN 
DO k = 1,config_flags%num_snow_layers
grid%snliqxy(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+config_flags%num_snow_layers-1+1
ENDIF
IF ( SIZE(grid%lfmassxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%lfmassxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%rtmassxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%rtmassxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%stmassxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%stmassxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%woodxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%woodxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%stblcpxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%stblcpxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%fastcpxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%fastcpxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%xsaixy) .GT. 1 ) THEN 
ixv=ixv+1
grid%xsaixy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%t2mvxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%t2mvxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%t2mbxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%t2mbxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%q2mvxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%q2mvxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%q2mbxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%q2mbxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tradxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%tradxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%neexy) .GT. 1 ) THEN 
ixv=ixv+1
grid%neexy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%gppxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%gppxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%nppxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%nppxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%fvegxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%fvegxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%qinxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%qinxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%runsfxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%runsfxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%runsbxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%runsbxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%ecanxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%ecanxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%edirxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%edirxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%etranxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%etranxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%fsaxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%fsaxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%firaxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%firaxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%aparxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%aparxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%psnxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%psnxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%savxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%savxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%sagxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%sagxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%rssunxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%rssunxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%rsshaxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%rsshaxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%bgapxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%bgapxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%wgapxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%wgapxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tgvxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%tgvxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%tgbxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%tgbxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%chvxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%chvxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%chbxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%chbxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%shgxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%shgxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%shcxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%shcxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%shbxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%shbxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%evgxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%evgxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%evbxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%evbxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%ghvxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%ghvxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%ghbxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%ghbxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%irgxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%irgxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%ircxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%ircxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%irbxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%irbxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%trxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%trxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%evcxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%evcxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%chleafxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%chleafxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%chucxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%chucxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%chv2xy) .GT. 1 ) THEN 
ixv=ixv+1
grid%chv2xy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%chb2xy) .GT. 1 ) THEN 
ixv=ixv+1
grid%chb2xy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%chstarxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%chstarxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%soilcomp) .GT. 1 ) THEN 
DO k = 1,8
grid%soilcomp(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+8-1+1
ENDIF
IF ( SIZE(grid%soilcl1) .GT. 1 ) THEN 
ixv=ixv+1
grid%soilcl1(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%soilcl2) .GT. 1 ) THEN 
ixv=ixv+1
grid%soilcl2(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%soilcl3) .GT. 1 ) THEN 
ixv=ixv+1
grid%soilcl3(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%soilcl4) .GT. 1 ) THEN 
ixv=ixv+1
grid%soilcl4(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%grainxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%grainxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%gddxy) .GT. 1 ) THEN 
ixv=ixv+1
grid%gddxy(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%croptype) .GT. 1 ) THEN 
DO k = 1,5
grid%croptype(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+5-1+1
ENDIF
IF ( SIZE(grid%planting) .GT. 1 ) THEN 
ixv=ixv+1
grid%planting(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%harvest) .GT. 1 ) THEN 
ixv=ixv+1
grid%harvest(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%season_gdd) .GT. 1 ) THEN 
ixv=ixv+1
grid%season_gdd(pig,pjg) = xv(ixv)
ENDIF
IF ( SIZE(grid%gecros_state) .GT. 1 ) THEN 
DO k = 1,60
grid%gecros_state(pig,k,pjg) = xv(k-1+ixv+1)
ENDDO
ixv=ixv+60-1+1
ENDIF
CALL rsl_lite_from_parent_info(pig,pjg,retval)
ENDDO


      CALL get_ijk_from_grid (  grid ,              &
                                ids, ide, jds, jde, kds, kde,    &
                                ims, ime, jms, jme, kms, kme,    &
                                ips, ipe, jps, jpe, kps, kpe    )







CALL HALO_FORCE_DOWN_sub ( grid, &
  config_flags, &
  num_moist, &
  moist, &
  num_scalar, &
  scalar, &
  num_dfi_scalar, &
  dfi_scalar, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )


      RETURN
    END SUBROUTINE force_intermediate_nmm



   SUBROUTINE feedback_domain_nmm_part1 ( grid, ngrid, config_flags    &







,szj,s1z,spz,tcs,moist,moist_bxs,moist_bxe,moist_bys,moist_bye,moist_btxs,moist_btxe,moist_btys,moist_btye,dfi_moist, &
dfi_moist_bxs,dfi_moist_bxe,dfi_moist_bys,dfi_moist_bye,dfi_moist_btxs,dfi_moist_btxe,dfi_moist_btys,dfi_moist_btye,scalar, &
scalar_bxs,scalar_bxe,scalar_bys,scalar_bye,scalar_btxs,scalar_btxe,scalar_btys,scalar_btye,dfi_scalar,dfi_scalar_bxs, &
dfi_scalar_bxe,dfi_scalar_bys,dfi_scalar_bye,dfi_scalar_btxs,dfi_scalar_btxe,dfi_scalar_btys,dfi_scalar_btye,chem,ozmixm &


                 )
      USE module_state_description
      USE module_domain, ONLY : domain, get_ijk_from_grid
      USE module_configure, ONLY : grid_config_rec_type, model_config_rec, model_to_grid_config_rec
      USE module_dm, ONLY : ntasks, ntasks_x, ntasks_y, itrace, local_communicator, mytask, &
                            ipe_save, jpe_save, ips_save, jps_save, get_dm_max_halo_width,  &
                            nest_pes_x, nest_pes_y
      IMPLICIT NONE

      TYPE(domain), POINTER :: grid          
      TYPE(domain), POINTER :: ngrid






real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_szj)           :: szj
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_s1z)           :: s1z
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_spz)           :: spz
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_tcs)           :: tcs
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_moist)           :: moist
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_bxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_bxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_bys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_bye
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_btxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_btxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_btys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_btye
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_dfi_moist)           :: dfi_moist
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_bxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_bxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_bys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_bye
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_btxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_btxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_btys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_btye
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_scalar)           :: scalar
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_bxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_bxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_bys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_bye
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_btxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_btxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_btys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_btye
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_dfi_scalar)           :: dfi_scalar
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_bxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_bxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_bys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_bye
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_btxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_btxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_btys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_btye
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%sm32:grid%em32,num_chem)           :: chem
real      ,DIMENSION(grid%sm31:grid%em31,1:grid%levsiz,grid%sm32:grid%em32,num_ozmixm)           :: ozmixm

      INTEGER nlev, msize, i_parent_start, j_parent_start
      INTEGER i,j,pig,pjg,cm,cn,nig,njg,retval,k
      TYPE(domain), POINTER :: xgrid
      TYPE (grid_config_rec_type)            :: config_flags, nconfig_flags
      REAL xv(2000000)
      INTEGER ixv
      INTEGER       ::          cids, cide, cjds, cjde, ckds, ckde,    &
                                cims, cime, cjms, cjme, ckms, ckme,    &
                                cips, cipe, cjps, cjpe, ckps, ckpe
      INTEGER       ::          nids, nide, njds, njde, nkds, nkde,    &
                                nims, nime, njms, njme, nkms, nkme,    &
                                nips, nipe, njps, njpe, nkps, nkpe

      INTEGER idim1,idim2,idim3,idim4,idim5,idim6,idim7

      INTEGER local_comm, myproc, nproc, idum1, idum2

      integer, parameter :: EConst=0, ECopy=1, EExtrap=2 

      LOGICAL interp_mp

      INTERFACE
          SUBROUTINE feedback_nest_prep_nmm ( grid, config_flags    &







,szj,s1z,spz,tcs,moist,moist_bxs,moist_bxe,moist_bys,moist_bye,moist_btxs,moist_btxe,moist_btys,moist_btye,dfi_moist, &
dfi_moist_bxs,dfi_moist_bxe,dfi_moist_bys,dfi_moist_bye,dfi_moist_btxs,dfi_moist_btxe,dfi_moist_btys,dfi_moist_btye,scalar, &
scalar_bxs,scalar_bxe,scalar_bys,scalar_bye,scalar_btxs,scalar_btxe,scalar_btys,scalar_btye,dfi_scalar,dfi_scalar_bxs, &
dfi_scalar_bxe,dfi_scalar_bys,dfi_scalar_bye,dfi_scalar_btxs,dfi_scalar_btxe,dfi_scalar_btys,dfi_scalar_btye,chem,ozmixm &


)
             USE module_state_description
             USE module_domain, ONLY : domain
             USE module_configure, ONLY : grid_config_rec_type

             TYPE (grid_config_rec_type)            :: config_flags
             TYPE(domain), TARGET                   :: grid






real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_szj)           :: szj
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_s1z)           :: s1z
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_spz)           :: spz
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_tcs)           :: tcs
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_moist)           :: moist
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_bxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_bxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_bys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_bye
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_btxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_btxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_btys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_btye
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_dfi_moist)           :: dfi_moist
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_bxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_bxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_bys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_bye
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_btxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_btxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_btys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_btye
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_scalar)           :: scalar
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_bxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_bxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_bys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_bye
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_btxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_btxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_btys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_btye
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_dfi_scalar)           :: dfi_scalar
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_bxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_bxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_bys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_bye
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_btxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_btxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_btys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_btye
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%sm32:grid%em32,num_chem)           :: chem
real      ,DIMENSION(grid%sm31:grid%em31,1:grid%levsiz,grid%sm32:grid%em32,num_ozmixm)           :: ozmixm

          END SUBROUTINE feedback_nest_prep_nmm
      END INTERFACE


      interp_mp=grid%interp_mp .or. ngrid%interp_mp
      CALL wrf_get_dm_communicator ( local_comm )
      CALL wrf_get_myproc( myproc )
      CALL wrf_get_nproc( nproc )




      CALL get_ijk_from_grid (  grid ,                   &
                                cids, cide, cjds, cjde, ckds, ckde,    &
                                cims, cime, cjms, cjme, ckms, ckme,    &
                                cips, cipe, cjps, cjpe, ckps, ckpe    )

      CALL get_ijk_from_grid (  ngrid ,                  &
                                nids, nide, njds, njde, nkds, nkde,    &
                                nims, nime, njms, njme, nkms, nkme,    &
                                nips, nipe, njps, njpe, nkps, nkpe    )

      nlev  = ckde - ckds + 1

      ips_save = ngrid%i_parent_start  
      jps_save = ngrid%j_parent_start  
      ipe_save = ngrid%i_parent_start + (nide-nids) / ngrid%parent_grid_ratio - 1
      jpe_save = ngrid%j_parent_start + (njde-njds) / ngrid%parent_grid_ratio - 1








      CALL model_to_grid_config_rec ( ngrid%id , model_config_rec , nconfig_flags )
      CALL set_scalar_indices_from_config ( ngrid%id , idum1 , idum2 )
      xgrid => grid
      grid => ngrid

      CALL feedback_nest_prep_nmm ( grid, config_flags    &







,grid%szj,grid%s1z,grid%spz,grid%tcs,grid%moist,grid%moist_bxs,grid%moist_bxe,grid%moist_bys,grid%moist_bye,grid%moist_btxs, &
grid%moist_btxe,grid%moist_btys,grid%moist_btye,grid%dfi_moist,grid%dfi_moist_bxs,grid%dfi_moist_bxe,grid%dfi_moist_bys, &
grid%dfi_moist_bye,grid%dfi_moist_btxs,grid%dfi_moist_btxe,grid%dfi_moist_btys,grid%dfi_moist_btye,grid%scalar,grid%scalar_bxs, &
grid%scalar_bxe,grid%scalar_bys,grid%scalar_bye,grid%scalar_btxs,grid%scalar_btxe,grid%scalar_btys,grid%scalar_btye, &
grid%dfi_scalar,grid%dfi_scalar_bxs,grid%dfi_scalar_bxe,grid%dfi_scalar_bys,grid%dfi_scalar_bye,grid%dfi_scalar_btxs, &
grid%dfi_scalar_btxe,grid%dfi_scalar_btys,grid%dfi_scalar_btye,grid%chem,grid%ozmixm &


)



      grid => xgrid
      CALL set_scalar_indices_from_config ( grid%id , idum1 , idum2 )



      


      call store_interp_info(ngrid,grid)
      call ext_n2c_fulldom(&
           ngrid%deta1,ngrid%deta2,ngrid%eta1,             &
           ngrid%eta2,ngrid%pt,ngrid%pdtop,                &
           grid%pint,grid%t,grid%pd,grid%q,                &
           cids, cide, cjds, cjde, ckds, ckde,             &
           cims, cime, cjms, cjme, ckms, ckme,             &
           cips, cipe, cjps, cjpe, ckps, ckpe,             &
           ngrid%pint,ngrid%t,                             &
           ngrid%pd,ngrid%q,                               &
           ngrid%i_parent_start, ngrid%j_parent_start,     &
           grid%iinfo,grid%winfo,                          &
           nids, nide, njds, njde, nkds, nkde,             &
           nims, nime, njms, njme, nkms, nkme,             &
           nips, nipe, njps, njpe, nkps, nkpe)


      






IF ( SIZE( grid%lakedepth2d, 1 ) * SIZE( grid%lakedepth2d, 2 ) .GT. 1 ) THEN 
CALL UpNear (  &         
                  grid%lakedepth2d,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%lakedepth2d,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%savedtke12d, 1 ) * SIZE( grid%savedtke12d, 2 ) .GT. 1 ) THEN 
CALL UpNear (  &         
                  grid%savedtke12d,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%savedtke12d,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%snowdp2d, 1 ) * SIZE( grid%snowdp2d, 2 ) .GT. 1 ) THEN 
CALL UpNear (  &         
                  grid%snowdp2d,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%snowdp2d,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%h2osno2d, 1 ) * SIZE( grid%h2osno2d, 2 ) .GT. 1 ) THEN 
CALL UpNear (  &         
                  grid%h2osno2d,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%h2osno2d,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%snl2d, 1 ) * SIZE( grid%snl2d, 2 ) .GT. 1 ) THEN 
CALL UpNear (  &         
                  grid%snl2d,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%snl2d,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%t_grnd2d, 1 ) * SIZE( grid%t_grnd2d, 2 ) .GT. 1 ) THEN 
CALL UpNear (  &         
                  grid%t_grnd2d,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%t_grnd2d,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%t_lake3d, 1 ) * SIZE( grid%t_lake3d, 3 ) .GT. 1 ) THEN 
CALL UpNear (  &         
                  grid%t_lake3d,   &       
                 cids, cide, 1, 10, cjds, cjde,   &         
                 cims, cime, 1, 10, cjms, cjme,   &         
                 cips, cipe, 1, 10, cjps, cjpe,   &         
                  ngrid%t_lake3d,  &   
                 nids, nide, 1, 10, njds, njde,   &         
                 nims, nime, 1, 10, njms, njme,   &         
                 nips, nipe, 1, 10, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%lake_icefrac3d, 1 ) * SIZE( grid%lake_icefrac3d, 3 ) .GT. 1 ) THEN 
CALL UpNear (  &         
                  grid%lake_icefrac3d,   &       
                 cids, cide, 1, 10, cjds, cjde,   &         
                 cims, cime, 1, 10, cjms, cjme,   &         
                 cips, cipe, 1, 10, cjps, cjpe,   &         
                  ngrid%lake_icefrac3d,  &   
                 nids, nide, 1, 10, njds, njde,   &         
                 nims, nime, 1, 10, njms, njme,   &         
                 nips, nipe, 1, 10, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%z_lake3d, 1 ) * SIZE( grid%z_lake3d, 3 ) .GT. 1 ) THEN 
CALL UpNear (  &         
                  grid%z_lake3d,   &       
                 cids, cide, 1, 10, cjds, cjde,   &         
                 cims, cime, 1, 10, cjms, cjme,   &         
                 cips, cipe, 1, 10, cjps, cjpe,   &         
                  ngrid%z_lake3d,  &   
                 nids, nide, 1, 10, njds, njde,   &         
                 nims, nime, 1, 10, njms, njme,   &         
                 nips, nipe, 1, 10, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%dz_lake3d, 1 ) * SIZE( grid%dz_lake3d, 3 ) .GT. 1 ) THEN 
CALL UpNear (  &         
                  grid%dz_lake3d,   &       
                 cids, cide, 1, 10, cjds, cjde,   &         
                 cims, cime, 1, 10, cjms, cjme,   &         
                 cips, cipe, 1, 10, cjps, cjpe,   &         
                  ngrid%dz_lake3d,  &   
                 nids, nide, 1, 10, njds, njde,   &         
                 nims, nime, 1, 10, njms, njme,   &         
                 nips, nipe, 1, 10, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%t_soisno3d, 1 ) * SIZE( grid%t_soisno3d, 3 ) .GT. 1 ) THEN 
CALL UpNear (  &         
                  grid%t_soisno3d,   &       
                 cids, cide, 1, 15, cjds, cjde,   &         
                 cims, cime, 1, 15, cjms, cjme,   &         
                 cips, cipe, 1, 15, cjps, cjpe,   &         
                  ngrid%t_soisno3d,  &   
                 nids, nide, 1, 15, njds, njde,   &         
                 nims, nime, 1, 15, njms, njme,   &         
                 nips, nipe, 1, 15, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%h2osoi_ice3d, 1 ) * SIZE( grid%h2osoi_ice3d, 3 ) .GT. 1 ) THEN 
CALL UpNear (  &         
                  grid%h2osoi_ice3d,   &       
                 cids, cide, 1, 15, cjds, cjde,   &         
                 cims, cime, 1, 15, cjms, cjme,   &         
                 cips, cipe, 1, 15, cjps, cjpe,   &         
                  ngrid%h2osoi_ice3d,  &   
                 nids, nide, 1, 15, njds, njde,   &         
                 nims, nime, 1, 15, njms, njme,   &         
                 nips, nipe, 1, 15, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%h2osoi_liq3d, 1 ) * SIZE( grid%h2osoi_liq3d, 3 ) .GT. 1 ) THEN 
CALL UpNear (  &         
                  grid%h2osoi_liq3d,   &       
                 cids, cide, 1, 15, cjds, cjde,   &         
                 cims, cime, 1, 15, cjms, cjme,   &         
                 cips, cipe, 1, 15, cjps, cjpe,   &         
                  ngrid%h2osoi_liq3d,  &   
                 nids, nide, 1, 15, njds, njde,   &         
                 nims, nime, 1, 15, njms, njme,   &         
                 nips, nipe, 1, 15, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%h2osoi_vol3d, 1 ) * SIZE( grid%h2osoi_vol3d, 3 ) .GT. 1 ) THEN 
CALL UpNear (  &         
                  grid%h2osoi_vol3d,   &       
                 cids, cide, 1, 15, cjds, cjde,   &         
                 cims, cime, 1, 15, cjms, cjme,   &         
                 cips, cipe, 1, 15, cjps, cjpe,   &         
                  ngrid%h2osoi_vol3d,  &   
                 nids, nide, 1, 15, njds, njde,   &         
                 nims, nime, 1, 15, njms, njme,   &         
                 nips, nipe, 1, 15, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%z3d, 1 ) * SIZE( grid%z3d, 3 ) .GT. 1 ) THEN 
CALL UpNear (  &         
                  grid%z3d,   &       
                 cids, cide, 1, 15, cjds, cjde,   &         
                 cims, cime, 1, 15, cjms, cjme,   &         
                 cips, cipe, 1, 15, cjps, cjpe,   &         
                  ngrid%z3d,  &   
                 nids, nide, 1, 15, njds, njde,   &         
                 nims, nime, 1, 15, njms, njme,   &         
                 nips, nipe, 1, 15, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%dz3d, 1 ) * SIZE( grid%dz3d, 3 ) .GT. 1 ) THEN 
CALL UpNear (  &         
                  grid%dz3d,   &       
                 cids, cide, 1, 15, cjds, cjde,   &         
                 cims, cime, 1, 15, cjms, cjme,   &         
                 cips, cipe, 1, 15, cjps, cjpe,   &         
                  ngrid%dz3d,  &   
                 nids, nide, 1, 15, njds, njde,   &         
                 nims, nime, 1, 15, njms, njme,   &         
                 nips, nipe, 1, 15, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%zi3d, 1 ) * SIZE( grid%zi3d, 3 ) .GT. 1 ) THEN 
CALL UpNear (  &         
                  grid%zi3d,   &       
                 cids, cide, 1, 16, cjds, cjde,   &         
                 cims, cime, 1, 16, cjms, cjme,   &         
                 cips, cipe, 1, 16, cjps, cjpe,   &         
                  ngrid%zi3d,  &   
                 nids, nide, 1, 16, njds, njde,   &         
                 nims, nime, 1, 16, njms, njme,   &         
                 nips, nipe, 1, 16, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%watsat3d, 1 ) * SIZE( grid%watsat3d, 3 ) .GT. 1 ) THEN 
CALL UpNear (  &         
                  grid%watsat3d,   &       
                 cids, cide, 1, 10, cjds, cjde,   &         
                 cims, cime, 1, 10, cjms, cjme,   &         
                 cips, cipe, 1, 10, cjps, cjpe,   &         
                  ngrid%watsat3d,  &   
                 nids, nide, 1, 10, njds, njde,   &         
                 nims, nime, 1, 10, njms, njme,   &         
                 nips, nipe, 1, 10, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%csol3d, 1 ) * SIZE( grid%csol3d, 3 ) .GT. 1 ) THEN 
CALL UpNear (  &         
                  grid%csol3d,   &       
                 cids, cide, 1, 10, cjds, cjde,   &         
                 cims, cime, 1, 10, cjms, cjme,   &         
                 cips, cipe, 1, 10, cjps, cjpe,   &         
                  ngrid%csol3d,  &   
                 nids, nide, 1, 10, njds, njde,   &         
                 nims, nime, 1, 10, njms, njme,   &         
                 nips, nipe, 1, 10, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%tkmg3d, 1 ) * SIZE( grid%tkmg3d, 3 ) .GT. 1 ) THEN 
CALL UpNear (  &         
                  grid%tkmg3d,   &       
                 cids, cide, 1, 10, cjds, cjde,   &         
                 cims, cime, 1, 10, cjms, cjme,   &         
                 cips, cipe, 1, 10, cjps, cjpe,   &         
                  ngrid%tkmg3d,  &   
                 nids, nide, 1, 10, njds, njde,   &         
                 nims, nime, 1, 10, njms, njme,   &         
                 nips, nipe, 1, 10, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%tkdry3d, 1 ) * SIZE( grid%tkdry3d, 3 ) .GT. 1 ) THEN 
CALL UpNear (  &         
                  grid%tkdry3d,   &       
                 cids, cide, 1, 10, cjds, cjde,   &         
                 cims, cime, 1, 10, cjms, cjme,   &         
                 cips, cipe, 1, 10, cjps, cjpe,   &         
                  ngrid%tkdry3d,  &   
                 nids, nide, 1, 10, njds, njde,   &         
                 nims, nime, 1, 10, njms, njme,   &         
                 nips, nipe, 1, 10, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%tksatu3d, 1 ) * SIZE( grid%tksatu3d, 3 ) .GT. 1 ) THEN 
CALL UpNear (  &         
                  grid%tksatu3d,   &       
                 cids, cide, 1, 10, cjds, cjde,   &         
                 cims, cime, 1, 10, cjms, cjme,   &         
                 cips, cipe, 1, 10, cjps, cjpe,   &         
                  ngrid%tksatu3d,  &   
                 nids, nide, 1, 10, njds, njde,   &         
                 nims, nime, 1, 10, njms, njme,   &         
                 nips, nipe, 1, 10, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%lu_index, 1 ) * SIZE( grid%lu_index, 2 ) .GT. 1 ) THEN 
CONTINUE 
ENDIF
IF ( SIZE( grid%lu_mask, 1 ) * SIZE( grid%lu_mask, 2 ) .GT. 1 ) THEN 
CONTINUE 
ENDIF
IF ( SIZE( grid%cuprecip_swath, 1 ) * SIZE( grid%cuprecip_swath, 2 ) .GT. 1 ) THEN 
CALL upmax (  &         
                  grid%cuprecip_swath,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%cuprecip_swath,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%precip_swath, 1 ) * SIZE( grid%precip_swath, 2 ) .GT. 1 ) THEN 
CALL upmax (  &         
                  grid%precip_swath,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%precip_swath,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%windsq_swath, 1 ) * SIZE( grid%windsq_swath, 2 ) .GT. 1 ) THEN 
CALL upmax (  &         
                  grid%windsq_swath,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%windsq_swath,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%suswind, 1 ) * SIZE( grid%suswind, 2 ) .GT. 1 ) THEN 
CALL upmax (  &         
                  grid%suswind,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%suswind,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%suswind_swath, 1 ) * SIZE( grid%suswind_swath, 2 ) .GT. 1 ) THEN 
CALL upmax (  &         
                  grid%suswind_swath,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%suswind_swath,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%pdyn_smooth, 1 ) * SIZE( grid%pdyn_smooth, 2 ) .GT. 1 ) THEN 
CONTINUE 
ENDIF
IF ( SIZE( grid%pdyn_parent, 1 ) * SIZE( grid%pdyn_parent, 2 ) .GT. 1 ) THEN 
CONTINUE 
ENDIF
IF ( SIZE( grid%tg_max_m10wind, 1 ) * SIZE( grid%tg_max_m10wind, 2 ) .GT. 1 ) THEN 
CALL upmax (  &         
                  grid%tg_max_m10wind,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%tg_max_m10wind,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%tg_max_wwind, 1 ) * SIZE( grid%tg_max_wwind, 2 ) .GT. 1 ) THEN 
CALL upmax (  &         
                  grid%tg_max_wwind,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%tg_max_wwind,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%tg_min_wwind, 1 ) * SIZE( grid%tg_min_wwind, 2 ) .GT. 1 ) THEN 
CALL upmax (  &         
                  grid%tg_min_wwind,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%tg_min_wwind,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%tg_max_zhel_25, 1 ) * SIZE( grid%tg_max_zhel_25, 2 ) .GT. 1 ) THEN 
CALL upmax (  &         
                  grid%tg_max_zhel_25,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%tg_max_zhel_25,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%tg_min_zhel_25, 1 ) * SIZE( grid%tg_min_zhel_25, 2 ) .GT. 1 ) THEN 
CALL upmax (  &         
                  grid%tg_min_zhel_25,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%tg_min_zhel_25,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%tg_max_zhel_03, 1 ) * SIZE( grid%tg_max_zhel_03, 2 ) .GT. 1 ) THEN 
CALL upmax (  &         
                  grid%tg_max_zhel_03,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%tg_max_zhel_03,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%tg_min_zhel_03, 1 ) * SIZE( grid%tg_min_zhel_03, 2 ) .GT. 1 ) THEN 
CALL upmax (  &         
                  grid%tg_min_zhel_03,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%tg_min_zhel_03,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%tg_updhel25, 1 ) * SIZE( grid%tg_updhel25, 2 ) .GT. 1 ) THEN 
CALL upcopy (  &         
                  grid%tg_updhel25,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%tg_updhel25,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%tg_max_updhel25, 1 ) * SIZE( grid%tg_max_updhel25, 2 ) .GT. 1 ) THEN 
CALL upmax (  &         
                  grid%tg_max_updhel25,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%tg_max_updhel25,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%tg_updhel03, 1 ) * SIZE( grid%tg_updhel03, 2 ) .GT. 1 ) THEN 
CALL upcopy (  &         
                  grid%tg_updhel03,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%tg_updhel03,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%tg_max_updhel03, 1 ) * SIZE( grid%tg_max_updhel03, 2 ) .GT. 1 ) THEN 
CALL upmax (  &         
                  grid%tg_max_updhel03,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%tg_max_updhel03,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%tg_total_precip, 1 ) * SIZE( grid%tg_total_precip, 2 ) .GT. 1 ) THEN 
CALL upcopy (  &         
                  grid%tg_total_precip,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%tg_total_precip,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%hres_fis, 1 ) * SIZE( grid%hres_fis, 2 ) .GT. 1 ) THEN 
CONTINUE 
ENDIF
IF ( SIZE( grid%hres_lnd, 1 ) * SIZE( grid%hres_lnd, 2 ) .GT. 1 ) THEN 
CONTINUE 
ENDIF
IF ( SIZE( grid%hres_lu_index, 1 ) * SIZE( grid%hres_lu_index, 2 ) .GT. 1 ) THEN 
CONTINUE 
ENDIF
IF ( SIZE( grid%sm, 1 ) * SIZE( grid%sm, 2 ) .GT. 1 ) THEN 
CONTINUE 
ENDIF
IF ( SIZE( grid%pd, 1 ) * SIZE( grid%pd, 2 ) .GT. 1 ) THEN 
CONTINUE 
ENDIF
IF ( SIZE( grid%fis, 1 ) * SIZE( grid%fis, 2 ) .GT. 1 ) THEN 
CONTINUE 
ENDIF
IF ( SIZE( grid%t, 1 ) * SIZE( grid%t, 2 ) .GT. 1 ) THEN 
CONTINUE 
ENDIF
IF ( SIZE( grid%q, 1 ) * SIZE( grid%q, 2 ) .GT. 1 ) THEN 
CONTINUE 
ENDIF
IF ( SIZE( grid%u, 1 ) * SIZE( grid%u, 2 ) .GT. 1 ) THEN 
CALL upvel (  &         
                  grid%u,   &       
                 cids, cide, ckds, ckde, cjds, cjde,   &         
                 cims, cime, ckms, ckme, cjms, cjme,   &         
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         
                  ngrid%u,  &   
                 nids, nide, nkds, nkde, njds, njde,   &         
                 nims, nime, nkms, nkme, njms, njme,   &         
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%v, 1 ) * SIZE( grid%v, 2 ) .GT. 1 ) THEN 
CALL upvel (  &         
                  grid%v,   &       
                 cids, cide, ckds, ckde, cjds, cjde,   &         
                 cims, cime, ckms, ckme, cjms, cjme,   &         
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         
                  ngrid%v,  &   
                 nids, nide, nkds, nkde, njds, njde,   &         
                 nims, nime, nkms, nkme, njms, njme,   &         
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%q2, 1 ) * SIZE( grid%q2, 2 ) .GT. 1 ) THEN 
CALL upmass (  &         
                  grid%q2,   &       
                 cids, cide, ckds, ckde, cjds, cjde,   &         
                 cims, cime, ckms, ckme, cjms, cjme,   &         
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         
                  ngrid%q2,  &   
                 nids, nide, nkds, nkde, njds, njde,   &         
                 nims, nime, nkms, nkme, njms, njme,   &         
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,econst &
,0.000000000e+00 &
                  ) 
ENDIF
IF ( SIZE( grid%acswupt, 1 ) * SIZE( grid%acswupt, 2 ) .GT. 1 ) THEN 
CALL UpCopy (  &         
                  grid%acswupt,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%acswupt,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%acswuptc, 1 ) * SIZE( grid%acswuptc, 2 ) .GT. 1 ) THEN 
CALL UpCopy (  &         
                  grid%acswuptc,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%acswuptc,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%acswdnt, 1 ) * SIZE( grid%acswdnt, 2 ) .GT. 1 ) THEN 
CALL UpCopy (  &         
                  grid%acswdnt,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%acswdnt,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%acswdntc, 1 ) * SIZE( grid%acswdntc, 2 ) .GT. 1 ) THEN 
CALL UpCopy (  &         
                  grid%acswdntc,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%acswdntc,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%acswupb, 1 ) * SIZE( grid%acswupb, 2 ) .GT. 1 ) THEN 
CALL UpCopy (  &         
                  grid%acswupb,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%acswupb,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%acswupbc, 1 ) * SIZE( grid%acswupbc, 2 ) .GT. 1 ) THEN 
CALL UpCopy (  &         
                  grid%acswupbc,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%acswupbc,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%acswdnb, 1 ) * SIZE( grid%acswdnb, 2 ) .GT. 1 ) THEN 
CALL UpCopy (  &         
                  grid%acswdnb,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%acswdnb,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%acswdnbc, 1 ) * SIZE( grid%acswdnbc, 2 ) .GT. 1 ) THEN 
CALL UpCopy (  &         
                  grid%acswdnbc,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%acswdnbc,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%aclwupt, 1 ) * SIZE( grid%aclwupt, 2 ) .GT. 1 ) THEN 
CALL UpCopy (  &         
                  grid%aclwupt,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%aclwupt,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%aclwuptc, 1 ) * SIZE( grid%aclwuptc, 2 ) .GT. 1 ) THEN 
CALL UpCopy (  &         
                  grid%aclwuptc,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%aclwuptc,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%aclwdnt, 1 ) * SIZE( grid%aclwdnt, 2 ) .GT. 1 ) THEN 
CALL UpCopy (  &         
                  grid%aclwdnt,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%aclwdnt,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%aclwdntc, 1 ) * SIZE( grid%aclwdntc, 2 ) .GT. 1 ) THEN 
CALL UpCopy (  &         
                  grid%aclwdntc,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%aclwdntc,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%aclwupb, 1 ) * SIZE( grid%aclwupb, 2 ) .GT. 1 ) THEN 
CALL UpCopy (  &         
                  grid%aclwupb,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%aclwupb,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%aclwupbc, 1 ) * SIZE( grid%aclwupbc, 2 ) .GT. 1 ) THEN 
CALL UpCopy (  &         
                  grid%aclwupbc,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%aclwupbc,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%aclwdnb, 1 ) * SIZE( grid%aclwdnb, 2 ) .GT. 1 ) THEN 
CALL UpCopy (  &         
                  grid%aclwdnb,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%aclwdnb,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%aclwdnbc, 1 ) * SIZE( grid%aclwdnbc, 2 ) .GT. 1 ) THEN 
CALL UpCopy (  &         
                  grid%aclwdnbc,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%aclwdnbc,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%swupt, 1 ) * SIZE( grid%swupt, 2 ) .GT. 1 ) THEN 
CALL UpCopy (  &         
                  grid%swupt,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%swupt,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%swuptc, 1 ) * SIZE( grid%swuptc, 2 ) .GT. 1 ) THEN 
CALL UpCopy (  &         
                  grid%swuptc,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%swuptc,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%swdnt, 1 ) * SIZE( grid%swdnt, 2 ) .GT. 1 ) THEN 
CALL UpCopy (  &         
                  grid%swdnt,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%swdnt,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%swdntc, 1 ) * SIZE( grid%swdntc, 2 ) .GT. 1 ) THEN 
CALL UpCopy (  &         
                  grid%swdntc,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%swdntc,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%swupb, 1 ) * SIZE( grid%swupb, 2 ) .GT. 1 ) THEN 
CALL UpCopy (  &         
                  grid%swupb,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%swupb,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%swupbc, 1 ) * SIZE( grid%swupbc, 2 ) .GT. 1 ) THEN 
CALL UpCopy (  &         
                  grid%swupbc,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%swupbc,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%swdnb, 1 ) * SIZE( grid%swdnb, 2 ) .GT. 1 ) THEN 
CALL UpCopy (  &         
                  grid%swdnb,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%swdnb,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%swdnbc, 1 ) * SIZE( grid%swdnbc, 2 ) .GT. 1 ) THEN 
CALL UpCopy (  &         
                  grid%swdnbc,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%swdnbc,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%lwupt, 1 ) * SIZE( grid%lwupt, 2 ) .GT. 1 ) THEN 
CALL UpCopy (  &         
                  grid%lwupt,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%lwupt,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%lwuptc, 1 ) * SIZE( grid%lwuptc, 2 ) .GT. 1 ) THEN 
CALL UpCopy (  &         
                  grid%lwuptc,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%lwuptc,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%lwdnt, 1 ) * SIZE( grid%lwdnt, 2 ) .GT. 1 ) THEN 
CALL UpCopy (  &         
                  grid%lwdnt,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%lwdnt,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%lwdntc, 1 ) * SIZE( grid%lwdntc, 2 ) .GT. 1 ) THEN 
CALL UpCopy (  &         
                  grid%lwdntc,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%lwdntc,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%lwupb, 1 ) * SIZE( grid%lwupb, 2 ) .GT. 1 ) THEN 
CALL UpCopy (  &         
                  grid%lwupb,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%lwupb,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%lwupbc, 1 ) * SIZE( grid%lwupbc, 2 ) .GT. 1 ) THEN 
CALL UpCopy (  &         
                  grid%lwupbc,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%lwupbc,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%lwdnb, 1 ) * SIZE( grid%lwdnb, 2 ) .GT. 1 ) THEN 
CALL UpCopy (  &         
                  grid%lwdnb,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%lwdnb,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%lwdnbc, 1 ) * SIZE( grid%lwdnbc, 2 ) .GT. 1 ) THEN 
CALL UpCopy (  &         
                  grid%lwdnbc,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%lwdnbc,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%qnwfa2d, 1 ) * SIZE( grid%qnwfa2d, 2 ) .GT. 1 ) THEN 
CALL UpCopy (  &         
                  grid%qnwfa2d,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%qnwfa2d,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%cwm, 1 ) * SIZE( grid%cwm, 2 ) .GT. 1 ) THEN 
CALL upmass (  &         
                  grid%cwm,   &       
                 cids, cide, ckds, ckde, cjds, cjde,   &         
                 cims, cime, ckms, ckme, cjms, cjme,   &         
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         
                  ngrid%cwm,  &   
                 nids, nide, nkds, nkde, njds, njde,   &         
                 nims, nime, nkms, nkme, njms, njme,   &         
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,ecopy &
,0.000000000e+00 &
                  ) 
ENDIF
IF ( SIZE( grid%f_ice, 1 ) * SIZE( grid%f_ice, 3 ) .GT. 1 ) THEN 
CALL upmassikj (  &         
                  grid%f_ice,   &       
                 cids, cide, ckds, ckde, cjds, cjde,   &         
                 cims, cime, ckms, ckme, cjms, cjme,   &         
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         
                  ngrid%f_ice,  &   
                 nids, nide, nkds, nkde, njds, njde,   &         
                 nims, nime, nkms, nkme, njms, njme,   &         
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,eextrap &
,0.000000000e+00 &
                  ) 
ENDIF
IF ( SIZE( grid%f_rain, 1 ) * SIZE( grid%f_rain, 3 ) .GT. 1 ) THEN 
CALL upmassikj (  &         
                  grid%f_rain,   &       
                 cids, cide, ckds, ckde, cjds, cjde,   &         
                 cims, cime, ckms, ckme, cjms, cjme,   &         
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         
                  ngrid%f_rain,  &   
                 nids, nide, nkds, nkde, njds, njde,   &         
                 nims, nime, nkms, nkme, njms, njme,   &         
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,eextrap &
,0.000000000e+00 &
                  ) 
ENDIF
IF ( SIZE( grid%f_rimef, 1 ) * SIZE( grid%f_rimef, 3 ) .GT. 1 ) THEN 
CALL upmassikj (  &         
                  grid%f_rimef,   &       
                 cids, cide, ckds, ckde, cjds, cjde,   &         
                 cims, cime, ckms, ckme, cjms, cjme,   &         
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         
                  ngrid%f_rimef,  &   
                 nids, nide, nkds, nkde, njds, njde,   &         
                 nims, nime, nkms, nkme, njms, njme,   &         
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,eextrap &
,1.000000000e+00 &
                  ) 
ENDIF
IF ( SIZE( grid%ctopo, 1 ) * SIZE( grid%ctopo, 2 ) .GT. 1 ) THEN 
CALL UpCopy (  &         
                  grid%ctopo,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%ctopo,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%ctopo2, 1 ) * SIZE( grid%ctopo2, 2 ) .GT. 1 ) THEN 
CALL UpCopy (  &         
                  grid%ctopo2,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%ctopo2,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%pint, 1 ) * SIZE( grid%pint, 2 ) .GT. 1 ) THEN 
CONTINUE 
ENDIF
IF ( SIZE( grid%winfo, 1 ) * SIZE( grid%winfo, 2 ) .GT. 1 ) THEN 
CONTINUE 
ENDIF
IF ( SIZE( grid%iinfo, 1 ) * SIZE( grid%iinfo, 2 ) .GT. 1 ) THEN 
CONTINUE 
ENDIF
DO itrace = PARAM_FIRST_SCALAR, num_moist
IF ( SIZE( moist, 1 ) * SIZE( moist, 2 ) .GT. 1 .and. (interp_mp .eqv. .true.) ) THEN 
CALL upmass (  &         
                  moist(grid%sm31,grid%sm32,grid%sm33,itrace),   &       
                 cids, cide, ckds, ckde, cjds, cjde,   &         
                 cims, cime, ckms, ckme, cjms, cjme,   &         
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         
                  ngrid%moist(ngrid%sm31,ngrid%sm32,ngrid%sm33,itrace),  &   
                 nids, nide, nkds, nkde, njds, njde,   &         
                 nims, nime, nkms, nkme, njms, njme,   &         
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,ecopy &
,0.000000000e+00 &
                  ) 
ENDIF
ENDDO
DO itrace = PARAM_FIRST_SCALAR, num_scalar
IF ( SIZE( scalar, 1 ) * SIZE( scalar, 2 ) .GT. 1 .and. (interp_mp .eqv. .true.) ) THEN 
CALL upmass (  &         
                  scalar(grid%sm31,grid%sm32,grid%sm33,itrace),   &       
                 cids, cide, ckds, ckde, cjds, cjde,   &         
                 cims, cime, ckms, ckme, cjms, cjme,   &         
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         
                  ngrid%scalar(ngrid%sm31,ngrid%sm32,ngrid%sm33,itrace),  &   
                 nids, nide, nkds, nkde, njds, njde,   &         
                 nims, nime, nkms, nkme, njms, njme,   &         
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,ecopy &
,0.000000000e+00 &
                  ) 
ENDIF
ENDDO
DO itrace = PARAM_FIRST_SCALAR, num_dfi_scalar
IF ( SIZE( dfi_scalar, 1 ) * SIZE( dfi_scalar, 2 ) .GT. 1 .and. (interp_mp .eqv. .true.) ) THEN 
CALL upmass (  &         
                  dfi_scalar(grid%sm31,grid%sm32,grid%sm33,itrace),   &       
                 cids, cide, ckds, ckde, cjds, cjde,   &         
                 cims, cime, ckms, ckme, cjms, cjme,   &         
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         
                  ngrid%dfi_scalar(ngrid%sm31,ngrid%sm32,ngrid%sm33,itrace),  &   
                 nids, nide, nkds, nkde, njds, njde,   &         
                 nims, nime, nkms, nkme, njms, njme,   &         
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
,ecopy &
,0.000000000e+00 &
                  ) 
ENDIF
ENDDO
IF ( SIZE( grid%isnowxy, 1 ) * SIZE( grid%isnowxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%isnowxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%isnowxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%tvxy, 1 ) * SIZE( grid%tvxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%tvxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%tvxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%tgxy, 1 ) * SIZE( grid%tgxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%tgxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%tgxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%canicexy, 1 ) * SIZE( grid%canicexy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%canicexy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%canicexy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%canliqxy, 1 ) * SIZE( grid%canliqxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%canliqxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%canliqxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%eahxy, 1 ) * SIZE( grid%eahxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%eahxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%eahxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%tahxy, 1 ) * SIZE( grid%tahxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%tahxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%tahxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%cmxy, 1 ) * SIZE( grid%cmxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%cmxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%cmxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%chxy, 1 ) * SIZE( grid%chxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%chxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%chxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%fwetxy, 1 ) * SIZE( grid%fwetxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%fwetxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%fwetxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%sneqvoxy, 1 ) * SIZE( grid%sneqvoxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%sneqvoxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%sneqvoxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%alboldxy, 1 ) * SIZE( grid%alboldxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%alboldxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%alboldxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%qsnowxy, 1 ) * SIZE( grid%qsnowxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%qsnowxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%qsnowxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%wslakexy, 1 ) * SIZE( grid%wslakexy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%wslakexy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%wslakexy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%zwtxy, 1 ) * SIZE( grid%zwtxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%zwtxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%zwtxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%waxy, 1 ) * SIZE( grid%waxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%waxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%waxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%wtxy, 1 ) * SIZE( grid%wtxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%wtxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%wtxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%tsnoxy, 1 ) * SIZE( grid%tsnoxy, 3 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%tsnoxy,   &       
                 cids, cide, 1, config_flags%num_snow_layers, cjds, cjde,   &         
                 cims, cime, 1, config_flags%num_snow_layers, cjms, cjme,   &         
                 cips, cipe, 1, config_flags%num_snow_layers, cjps, cjpe,   &         
                  ngrid%tsnoxy,  &   
                 nids, nide, 1, config_flags%num_snow_layers, njds, njde,   &         
                 nims, nime, 1, config_flags%num_snow_layers, njms, njme,   &         
                 nips, nipe, 1, config_flags%num_snow_layers, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%zsnsoxy, 1 ) * SIZE( grid%zsnsoxy, 3 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%zsnsoxy,   &       
                 cids, cide, 1, config_flags%num_snso_layers, cjds, cjde,   &         
                 cims, cime, 1, config_flags%num_snso_layers, cjms, cjme,   &         
                 cips, cipe, 1, config_flags%num_snso_layers, cjps, cjpe,   &         
                  ngrid%zsnsoxy,  &   
                 nids, nide, 1, config_flags%num_snso_layers, njds, njde,   &         
                 nims, nime, 1, config_flags%num_snso_layers, njms, njme,   &         
                 nips, nipe, 1, config_flags%num_snso_layers, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%snicexy, 1 ) * SIZE( grid%snicexy, 3 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%snicexy,   &       
                 cids, cide, 1, config_flags%num_snow_layers, cjds, cjde,   &         
                 cims, cime, 1, config_flags%num_snow_layers, cjms, cjme,   &         
                 cips, cipe, 1, config_flags%num_snow_layers, cjps, cjpe,   &         
                  ngrid%snicexy,  &   
                 nids, nide, 1, config_flags%num_snow_layers, njds, njde,   &         
                 nims, nime, 1, config_flags%num_snow_layers, njms, njme,   &         
                 nips, nipe, 1, config_flags%num_snow_layers, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%snliqxy, 1 ) * SIZE( grid%snliqxy, 3 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%snliqxy,   &       
                 cids, cide, 1, config_flags%num_snow_layers, cjds, cjde,   &         
                 cims, cime, 1, config_flags%num_snow_layers, cjms, cjme,   &         
                 cips, cipe, 1, config_flags%num_snow_layers, cjps, cjpe,   &         
                  ngrid%snliqxy,  &   
                 nids, nide, 1, config_flags%num_snow_layers, njds, njde,   &         
                 nims, nime, 1, config_flags%num_snow_layers, njms, njme,   &         
                 nips, nipe, 1, config_flags%num_snow_layers, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%lfmassxy, 1 ) * SIZE( grid%lfmassxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%lfmassxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%lfmassxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%rtmassxy, 1 ) * SIZE( grid%rtmassxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%rtmassxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%rtmassxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%stmassxy, 1 ) * SIZE( grid%stmassxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%stmassxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%stmassxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%woodxy, 1 ) * SIZE( grid%woodxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%woodxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%woodxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%stblcpxy, 1 ) * SIZE( grid%stblcpxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%stblcpxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%stblcpxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%fastcpxy, 1 ) * SIZE( grid%fastcpxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%fastcpxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%fastcpxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%xsaixy, 1 ) * SIZE( grid%xsaixy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%xsaixy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%xsaixy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%t2mvxy, 1 ) * SIZE( grid%t2mvxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%t2mvxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%t2mvxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%t2mbxy, 1 ) * SIZE( grid%t2mbxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%t2mbxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%t2mbxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%q2mvxy, 1 ) * SIZE( grid%q2mvxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%q2mvxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%q2mvxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%q2mbxy, 1 ) * SIZE( grid%q2mbxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%q2mbxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%q2mbxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%tradxy, 1 ) * SIZE( grid%tradxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%tradxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%tradxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%neexy, 1 ) * SIZE( grid%neexy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%neexy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%neexy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%gppxy, 1 ) * SIZE( grid%gppxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%gppxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%gppxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%nppxy, 1 ) * SIZE( grid%nppxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%nppxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%nppxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%fvegxy, 1 ) * SIZE( grid%fvegxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%fvegxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%fvegxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%qinxy, 1 ) * SIZE( grid%qinxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%qinxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%qinxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%runsfxy, 1 ) * SIZE( grid%runsfxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%runsfxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%runsfxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%runsbxy, 1 ) * SIZE( grid%runsbxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%runsbxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%runsbxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%ecanxy, 1 ) * SIZE( grid%ecanxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%ecanxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%ecanxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%edirxy, 1 ) * SIZE( grid%edirxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%edirxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%edirxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%etranxy, 1 ) * SIZE( grid%etranxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%etranxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%etranxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%fsaxy, 1 ) * SIZE( grid%fsaxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%fsaxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%fsaxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%firaxy, 1 ) * SIZE( grid%firaxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%firaxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%firaxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%aparxy, 1 ) * SIZE( grid%aparxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%aparxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%aparxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%psnxy, 1 ) * SIZE( grid%psnxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%psnxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%psnxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%savxy, 1 ) * SIZE( grid%savxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%savxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%savxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%sagxy, 1 ) * SIZE( grid%sagxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%sagxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%sagxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%rssunxy, 1 ) * SIZE( grid%rssunxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%rssunxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%rssunxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%rsshaxy, 1 ) * SIZE( grid%rsshaxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%rsshaxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%rsshaxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%bgapxy, 1 ) * SIZE( grid%bgapxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%bgapxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%bgapxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%wgapxy, 1 ) * SIZE( grid%wgapxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%wgapxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%wgapxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%tgvxy, 1 ) * SIZE( grid%tgvxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%tgvxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%tgvxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%tgbxy, 1 ) * SIZE( grid%tgbxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%tgbxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%tgbxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%chvxy, 1 ) * SIZE( grid%chvxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%chvxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%chvxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%chbxy, 1 ) * SIZE( grid%chbxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%chbxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%chbxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%shgxy, 1 ) * SIZE( grid%shgxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%shgxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%shgxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%shcxy, 1 ) * SIZE( grid%shcxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%shcxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%shcxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%shbxy, 1 ) * SIZE( grid%shbxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%shbxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%shbxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%evgxy, 1 ) * SIZE( grid%evgxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%evgxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%evgxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%evbxy, 1 ) * SIZE( grid%evbxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%evbxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%evbxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%ghvxy, 1 ) * SIZE( grid%ghvxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%ghvxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%ghvxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%ghbxy, 1 ) * SIZE( grid%ghbxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%ghbxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%ghbxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%irgxy, 1 ) * SIZE( grid%irgxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%irgxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%irgxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%ircxy, 1 ) * SIZE( grid%ircxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%ircxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%ircxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%irbxy, 1 ) * SIZE( grid%irbxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%irbxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%irbxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%trxy, 1 ) * SIZE( grid%trxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%trxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%trxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%evcxy, 1 ) * SIZE( grid%evcxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%evcxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%evcxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%chleafxy, 1 ) * SIZE( grid%chleafxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%chleafxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%chleafxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%chucxy, 1 ) * SIZE( grid%chucxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%chucxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%chucxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%chv2xy, 1 ) * SIZE( grid%chv2xy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%chv2xy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%chv2xy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%chb2xy, 1 ) * SIZE( grid%chb2xy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%chb2xy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%chb2xy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%chstarxy, 1 ) * SIZE( grid%chstarxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%chstarxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%chstarxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%soilcomp, 1 ) * SIZE( grid%soilcomp, 3 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%soilcomp,   &       
                 cids, cide, 1, 8, cjds, cjde,   &         
                 cims, cime, 1, 8, cjms, cjme,   &         
                 cips, cipe, 1, 8, cjps, cjpe,   &         
                  ngrid%soilcomp,  &   
                 nids, nide, 1, 8, njds, njde,   &         
                 nims, nime, 1, 8, njms, njme,   &         
                 nips, nipe, 1, 8, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%soilcl1, 1 ) * SIZE( grid%soilcl1, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%soilcl1,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%soilcl1,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%soilcl2, 1 ) * SIZE( grid%soilcl2, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%soilcl2,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%soilcl2,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%soilcl3, 1 ) * SIZE( grid%soilcl3, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%soilcl3,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%soilcl3,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%soilcl4, 1 ) * SIZE( grid%soilcl4, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%soilcl4,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%soilcl4,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%grainxy, 1 ) * SIZE( grid%grainxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%grainxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%grainxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%gddxy, 1 ) * SIZE( grid%gddxy, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%gddxy,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%gddxy,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%croptype, 1 ) * SIZE( grid%croptype, 3 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%croptype,   &       
                 cids, cide, 1, 5, cjds, cjde,   &         
                 cims, cime, 1, 5, cjms, cjme,   &         
                 cips, cipe, 1, 5, cjps, cjpe,   &         
                  ngrid%croptype,  &   
                 nids, nide, 1, 5, njds, njde,   &         
                 nims, nime, 1, 5, njms, njme,   &         
                 nips, nipe, 1, 5, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%planting, 1 ) * SIZE( grid%planting, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%planting,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%planting,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%harvest, 1 ) * SIZE( grid%harvest, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%harvest,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%harvest,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%season_gdd, 1 ) * SIZE( grid%season_gdd, 2 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%season_gdd,   &       
                 cids, cide, 1, 1, cjds, cjde,   &         
                 cims, cime, 1, 1, cjms, cjme,   &         
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  ngrid%season_gdd,  &   
                 nids, nide, 1, 1, njds, njde,   &         
                 nims, nime, 1, 1, njms, njme,   &         
                 nips, nipe, 1, 1, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
IF ( SIZE( grid%gecros_state, 1 ) * SIZE( grid%gecros_state, 3 ) .GT. 1 ) THEN 
CALL upnear (  &         
                  grid%gecros_state,   &       
                 cids, cide, 1, 60, cjds, cjde,   &         
                 cims, cime, 1, 60, cjms, cjme,   &         
                 cips, cipe, 1, 60, cjps, cjpe,   &         
                  ngrid%gecros_state,  &   
                 nids, nide, 1, 60, njds, njde,   &         
                 nims, nime, 1, 60, njms, njme,   &         
                 nips, nipe, 1, 60, njps, njpe,   &         
                  config_flags%shw, ngrid%imask_nostag,         &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF


      RETURN
   END SUBROUTINE feedback_domain_nmm_part1



   SUBROUTINE feedback_domain_nmm_part2 ( grid, intermediate_grid, ngrid , config_flags    &







,szj,s1z,spz,tcs,moist,moist_bxs,moist_bxe,moist_bys,moist_bye,moist_btxs,moist_btxe,moist_btys,moist_btye,dfi_moist, &
dfi_moist_bxs,dfi_moist_bxe,dfi_moist_bys,dfi_moist_bye,dfi_moist_btxs,dfi_moist_btxe,dfi_moist_btys,dfi_moist_btye,scalar, &
scalar_bxs,scalar_bxe,scalar_bys,scalar_bye,scalar_btxs,scalar_btxe,scalar_btys,scalar_btye,dfi_scalar,dfi_scalar_bxs, &
dfi_scalar_bxe,dfi_scalar_bys,dfi_scalar_bye,dfi_scalar_btxs,dfi_scalar_btxe,dfi_scalar_btys,dfi_scalar_btye,chem,ozmixm &


                 )
      USE module_state_description
      USE module_domain, ONLY : domain, domain_clock_get, get_ijk_from_grid
      USE module_configure, ONLY : grid_config_rec_type
      USE module_dm, ONLY : get_dm_max_halo_width, ips_save, ipe_save,              &
                            jps_save, jpe_save, ntasks, mytask, ntasks_x, ntasks_y, &
                            local_communicator, itrace,                             &
                            nest_pes_x, nest_pes_y,                                 &
                            intercomm_active, nest_task_offsets,                    &
                            mpi_comm_to_mom, mpi_comm_to_kid, which_kid 
                            

      USE module_comm_nesting_dm, ONLY : halo_interp_up_sub
      USE module_utility
      IMPLICIT NONE


      TYPE(domain), POINTER :: grid          
      TYPE(domain), POINTER :: intermediate_grid
      TYPE(domain), POINTER :: ngrid
      TYPE(domain), POINTER :: parent_grid







real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_szj)           :: szj
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_s1z)           :: s1z
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_spz)           :: spz
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_tcs)           :: tcs
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_moist)           :: moist
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_bxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_bxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_bys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_bye
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_btxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_btxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_btys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_btye
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_dfi_moist)           :: dfi_moist
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_bxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_bxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_bys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_bye
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_btxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_btxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_btys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_btye
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_scalar)           :: scalar
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_bxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_bxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_bys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_bye
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_btxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_btxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_btys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_btye
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_dfi_scalar)           :: dfi_scalar
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_bxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_bxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_bys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_bye
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_btxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_btxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_btys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_btye
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%sm32:grid%em32,num_chem)           :: chem
real      ,DIMENSION(grid%sm31:grid%em31,1:grid%levsiz,grid%sm32:grid%em32,num_ozmixm)           :: ozmixm

      INTEGER nlev, msize
      INTEGER i,j,pig,pjg,cm,cn,nig,njg,retval,k
      TYPE (grid_config_rec_type)            :: config_flags
      REAL xv(2000000)
      INTEGER ixv
      INTEGER       ::          cids, cide, cjds, cjde, ckds, ckde,    &
                                cims, cime, cjms, cjme, ckms, ckme,    &
                                cips, cipe, cjps, cjpe, ckps, ckpe
      INTEGER       ::          nids, nide, njds, njde, nkds, nkde,    &
                                nims, nime, njms, njme, nkms, nkme,    &
                                nips, nipe, njps, njpe, nkps, nkpe
      INTEGER       ::          xids, xide, xjds, xjde, xkds, xkde,    &
                                xims, xime, xjms, xjme, xkms, xkme,    &
                                xips, xipe, xjps, xjpe, xkps, xkpe
      INTEGER       ::          ids, ide, jds, jde, kds, kde,    &
                                ims, ime, jms, jme, kms, kme,    &
                                ips, ipe, jps, jpe, kps, kpe

      INTEGER idim1,idim2,idim3,idim4,idim5,idim6,idim7

      INTEGER icoord, jcoord, idim_cd, jdim_cd
      INTEGER local_comm, myproc, nproc
      INTEGER iparstrt, jparstrt, sw
      INTEGER thisdomain_max_halo_width

      character*256 :: timestr
      integer ioffset, ierr

      REAL    nest_influence
      LOGICAL feedback_flag, feedback_flag_v
      LOGICAL, EXTERNAL  :: cd_feedback_mask
      LOGICAL, EXTERNAL  :: cd_feedback_mask_v

      LOGICAL interp_mp















      interp_mp=grid%interp_mp .or. ngrid%interp_mp
      nest_influence = 0.5


      CALL domain_clock_get( grid, current_timestr=timestr )

      CALL get_ijk_from_grid (  intermediate_grid ,                   &
                                cids, cide, cjds, cjde, ckds, ckde,    &
                                cims, cime, cjms, cjme, ckms, ckme,    &
                                cips, cipe, cjps, cjpe, ckps, ckpe    )
      CALL get_ijk_from_grid (  grid ,              &
                                nids, nide, njds, njde, nkds, nkde,    &
                                nims, nime, njms, njme, nkms, nkme,    &
                                nips, nipe, njps, njpe, nkps, nkpe    )
      CALL get_ijk_from_grid (  ngrid ,              &
                                xids, xide, xjds, xjde, xkds, xkde,    &
                                xims, xime, xjms, xjme, xkms, xkme,    &
                                xips, xipe, xjps, xjpe, xkps, xkpe    )

      ips_save = ngrid%i_parent_start
      jps_save = ngrid%j_parent_start
      ipe_save = ngrid%i_parent_start + (xide-xids) / ngrid%parent_grid_ratio - 1
      jpe_save = ngrid%j_parent_start + (xjde-xjds) / ngrid%parent_grid_ratio - 1


      nide = nide - 1   
      njde = njde - 1   

IF ( ngrid%active_this_task ) THEN
    CALL push_communicators_for_domain( ngrid%id )
      CALL nl_get_i_parent_start ( intermediate_grid%id, iparstrt )
      CALL nl_get_j_parent_start ( intermediate_grid%id, jparstrt )
      CALL nl_get_shw            ( intermediate_grid%id, sw )
      icoord =    iparstrt  - sw
      jcoord =    jparstrt  - sw
      idim_cd = cide - cids + 1
      jdim_cd = cjde - cjds + 1

      nlev  = ckde - ckds + 1

      CALL get_dm_max_halo_width ( ngrid%id , thisdomain_max_halo_width )
      parent_grid => grid
      grid => ngrid






msize = (35 + ((num_szj - PARAM_FIRST_SCALAR + 1)) & 
 + ((num_s1z - PARAM_FIRST_SCALAR + 1)) & 
 + ((num_spz - PARAM_FIRST_SCALAR + 1)) & 
 + ((num_tcs - PARAM_FIRST_SCALAR + 1)) & 
 + ((num_chem - PARAM_FIRST_SCALAR + 1)) & 
 + ((num_ozmixm - PARAM_FIRST_SCALAR + 1)) )* nlev + 147
IF(interp_mp .eqv. .true.) then
    msize=msize + (0 + ((num_moist - PARAM_FIRST_SCALAR + 1)) & 
 + ((num_dfi_moist - PARAM_FIRST_SCALAR + 1)) & 
 + ((num_scalar - PARAM_FIRST_SCALAR + 1)) & 
 + ((num_dfi_scalar - PARAM_FIRST_SCALAR + 1)) )*nlev+0
ENDIF
CALL rsl_lite_to_parent_info( msize*4                               &
                        ,cips,cipe,cjps,cjpe                               &
                        ,nids,nide,njds,njde                               &
                        ,nest_task_offsets(ngrid%id)                      &
                        ,nest_pes_x(parent_grid%id)                            &
                        ,nest_pes_y(parent_grid%id)                            &
                        ,nest_pes_x(intermediate_grid%id)                 &
                        ,nest_pes_y(intermediate_grid%id)                 &
                        ,thisdomain_max_halo_width                         &
                        ,icoord,jcoord                                     &
                        ,idim_cd,jdim_cd                                   &
                        ,pig,pjg,retval )
DO while ( retval .eq. 1 )
ixv=0 
IF ( SIZE(grid%lakedepth2d) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%lakedepth2d(pig,pjg)
ENDIF
IF ( SIZE(grid%savedtke12d) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%savedtke12d(pig,pjg)
ENDIF
IF ( SIZE(grid%snowdp2d) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%snowdp2d(pig,pjg)
ENDIF
IF ( SIZE(grid%h2osno2d) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%h2osno2d(pig,pjg)
ENDIF
IF ( SIZE(grid%snl2d) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%snl2d(pig,pjg)
ENDIF
IF ( SIZE(grid%t_grnd2d) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%t_grnd2d(pig,pjg)
ENDIF
IF ( SIZE(grid%t_lake3d) .GT. 1 ) THEN 
DO k = 1,10
xv(ixv+1+k-1)= intermediate_grid%t_lake3d(pig,k,pjg)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%lake_icefrac3d) .GT. 1 ) THEN 
DO k = 1,10
xv(ixv+1+k-1)= intermediate_grid%lake_icefrac3d(pig,k,pjg)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%z_lake3d) .GT. 1 ) THEN 
DO k = 1,10
xv(ixv+1+k-1)= intermediate_grid%z_lake3d(pig,k,pjg)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%dz_lake3d) .GT. 1 ) THEN 
DO k = 1,10
xv(ixv+1+k-1)= intermediate_grid%dz_lake3d(pig,k,pjg)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%t_soisno3d) .GT. 1 ) THEN 
DO k = 1,15
xv(ixv+1+k-1)= intermediate_grid%t_soisno3d(pig,k,pjg)
ENDDO
ixv=ixv+15-1+1
ENDIF
IF ( SIZE(grid%h2osoi_ice3d) .GT. 1 ) THEN 
DO k = 1,15
xv(ixv+1+k-1)= intermediate_grid%h2osoi_ice3d(pig,k,pjg)
ENDDO
ixv=ixv+15-1+1
ENDIF
IF ( SIZE(grid%h2osoi_liq3d) .GT. 1 ) THEN 
DO k = 1,15
xv(ixv+1+k-1)= intermediate_grid%h2osoi_liq3d(pig,k,pjg)
ENDDO
ixv=ixv+15-1+1
ENDIF
IF ( SIZE(grid%h2osoi_vol3d) .GT. 1 ) THEN 
DO k = 1,15
xv(ixv+1+k-1)= intermediate_grid%h2osoi_vol3d(pig,k,pjg)
ENDDO
ixv=ixv+15-1+1
ENDIF
IF ( SIZE(grid%z3d) .GT. 1 ) THEN 
DO k = 1,15
xv(ixv+1+k-1)= intermediate_grid%z3d(pig,k,pjg)
ENDDO
ixv=ixv+15-1+1
ENDIF
IF ( SIZE(grid%dz3d) .GT. 1 ) THEN 
DO k = 1,15
xv(ixv+1+k-1)= intermediate_grid%dz3d(pig,k,pjg)
ENDDO
ixv=ixv+15-1+1
ENDIF
IF ( SIZE(grid%zi3d) .GT. 1 ) THEN 
DO k = 1,16
xv(ixv+1+k-1)= intermediate_grid%zi3d(pig,k,pjg)
ENDDO
ixv=ixv+16-1+1
ENDIF
IF ( SIZE(grid%watsat3d) .GT. 1 ) THEN 
DO k = 1,10
xv(ixv+1+k-1)= intermediate_grid%watsat3d(pig,k,pjg)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%csol3d) .GT. 1 ) THEN 
DO k = 1,10
xv(ixv+1+k-1)= intermediate_grid%csol3d(pig,k,pjg)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%tkmg3d) .GT. 1 ) THEN 
DO k = 1,10
xv(ixv+1+k-1)= intermediate_grid%tkmg3d(pig,k,pjg)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%tkdry3d) .GT. 1 ) THEN 
DO k = 1,10
xv(ixv+1+k-1)= intermediate_grid%tkdry3d(pig,k,pjg)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%tksatu3d) .GT. 1 ) THEN 
DO k = 1,10
xv(ixv+1+k-1)= intermediate_grid%tksatu3d(pig,k,pjg)
ENDDO
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%lu_index) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%lu_index(pig,pjg)
ENDIF
IF ( SIZE(grid%lu_mask) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%lu_mask(pig,pjg)
ENDIF
IF ( SIZE(grid%cuprecip_swath) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%cuprecip_swath(pig,pjg)
ENDIF
IF ( SIZE(grid%precip_swath) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%precip_swath(pig,pjg)
ENDIF
IF ( SIZE(grid%windsq_swath) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%windsq_swath(pig,pjg)
ENDIF
IF ( SIZE(grid%suswind) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%suswind(pig,pjg)
ENDIF
IF ( SIZE(grid%suswind_swath) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%suswind_swath(pig,pjg)
ENDIF
IF ( SIZE(grid%pdyn_smooth) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%pdyn_smooth(pig,pjg)
ENDIF
IF ( SIZE(grid%pdyn_parent) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%pdyn_parent(pig,pjg)
ENDIF
IF ( SIZE(grid%tg_max_m10wind) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%tg_max_m10wind(pig,pjg)
ENDIF
IF ( SIZE(grid%tg_max_wwind) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%tg_max_wwind(pig,pjg)
ENDIF
IF ( SIZE(grid%tg_min_wwind) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%tg_min_wwind(pig,pjg)
ENDIF
IF ( SIZE(grid%tg_max_zhel_25) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%tg_max_zhel_25(pig,pjg)
ENDIF
IF ( SIZE(grid%tg_min_zhel_25) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%tg_min_zhel_25(pig,pjg)
ENDIF
IF ( SIZE(grid%tg_max_zhel_03) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%tg_max_zhel_03(pig,pjg)
ENDIF
IF ( SIZE(grid%tg_min_zhel_03) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%tg_min_zhel_03(pig,pjg)
ENDIF
IF ( SIZE(grid%tg_updhel25) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%tg_updhel25(pig,pjg)
ENDIF
IF ( SIZE(grid%tg_max_updhel25) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%tg_max_updhel25(pig,pjg)
ENDIF
IF ( SIZE(grid%tg_updhel03) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%tg_updhel03(pig,pjg)
ENDIF
IF ( SIZE(grid%tg_max_updhel03) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%tg_max_updhel03(pig,pjg)
ENDIF
IF ( SIZE(grid%tg_total_precip) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%tg_total_precip(pig,pjg)
ENDIF
IF ( SIZE(grid%hres_fis) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%hres_fis(pig,pjg)
ENDIF
IF ( SIZE(grid%hres_lnd) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%hres_lnd(pig,pjg)
ENDIF
IF ( SIZE(grid%hres_lu_index) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%hres_lu_index(pig,pjg)
ENDIF
IF ( SIZE(grid%sm) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%sm(pig,pjg)
ENDIF
IF ( SIZE(grid%pd) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%pd(pig,pjg)
ENDIF
IF ( SIZE(grid%fis) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%fis(pig,pjg)
ENDIF
IF ( SIZE(grid%t) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
xv(ixv+1+k-ckds)= intermediate_grid%t(pig,pjg,k)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%q) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
xv(ixv+1+k-ckds)= intermediate_grid%q(pig,pjg,k)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%u) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
xv(ixv+1+k-ckds)= intermediate_grid%u(pig,pjg,k)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%v) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
xv(ixv+1+k-ckds)= intermediate_grid%v(pig,pjg,k)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%q2) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
xv(ixv+1+k-ckds)= intermediate_grid%q2(pig,pjg,k)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%acswupt) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%acswupt(pig,pjg)
ENDIF
IF ( SIZE(grid%acswuptc) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%acswuptc(pig,pjg)
ENDIF
IF ( SIZE(grid%acswdnt) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%acswdnt(pig,pjg)
ENDIF
IF ( SIZE(grid%acswdntc) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%acswdntc(pig,pjg)
ENDIF
IF ( SIZE(grid%acswupb) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%acswupb(pig,pjg)
ENDIF
IF ( SIZE(grid%acswupbc) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%acswupbc(pig,pjg)
ENDIF
IF ( SIZE(grid%acswdnb) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%acswdnb(pig,pjg)
ENDIF
IF ( SIZE(grid%acswdnbc) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%acswdnbc(pig,pjg)
ENDIF
IF ( SIZE(grid%aclwupt) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%aclwupt(pig,pjg)
ENDIF
IF ( SIZE(grid%aclwuptc) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%aclwuptc(pig,pjg)
ENDIF
IF ( SIZE(grid%aclwdnt) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%aclwdnt(pig,pjg)
ENDIF
IF ( SIZE(grid%aclwdntc) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%aclwdntc(pig,pjg)
ENDIF
IF ( SIZE(grid%aclwupb) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%aclwupb(pig,pjg)
ENDIF
IF ( SIZE(grid%aclwupbc) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%aclwupbc(pig,pjg)
ENDIF
IF ( SIZE(grid%aclwdnb) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%aclwdnb(pig,pjg)
ENDIF
IF ( SIZE(grid%aclwdnbc) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%aclwdnbc(pig,pjg)
ENDIF
IF ( SIZE(grid%swupt) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%swupt(pig,pjg)
ENDIF
IF ( SIZE(grid%swuptc) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%swuptc(pig,pjg)
ENDIF
IF ( SIZE(grid%swdnt) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%swdnt(pig,pjg)
ENDIF
IF ( SIZE(grid%swdntc) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%swdntc(pig,pjg)
ENDIF
IF ( SIZE(grid%swupb) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%swupb(pig,pjg)
ENDIF
IF ( SIZE(grid%swupbc) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%swupbc(pig,pjg)
ENDIF
IF ( SIZE(grid%swdnb) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%swdnb(pig,pjg)
ENDIF
IF ( SIZE(grid%swdnbc) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%swdnbc(pig,pjg)
ENDIF
IF ( SIZE(grid%lwupt) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%lwupt(pig,pjg)
ENDIF
IF ( SIZE(grid%lwuptc) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%lwuptc(pig,pjg)
ENDIF
IF ( SIZE(grid%lwdnt) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%lwdnt(pig,pjg)
ENDIF
IF ( SIZE(grid%lwdntc) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%lwdntc(pig,pjg)
ENDIF
IF ( SIZE(grid%lwupb) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%lwupb(pig,pjg)
ENDIF
IF ( SIZE(grid%lwupbc) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%lwupbc(pig,pjg)
ENDIF
IF ( SIZE(grid%lwdnb) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%lwdnb(pig,pjg)
ENDIF
IF ( SIZE(grid%lwdnbc) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%lwdnbc(pig,pjg)
ENDIF
IF ( SIZE(grid%qnwfa2d) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%qnwfa2d(pig,pjg)
ENDIF
IF ( SIZE(grid%cwm) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
xv(ixv+1+k-ckds)= intermediate_grid%cwm(pig,pjg,k)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%f_ice) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
xv(ixv+1+k-ckds)= intermediate_grid%f_ice(pig,k,pjg)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%f_rain) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
xv(ixv+1+k-ckds)= intermediate_grid%f_rain(pig,k,pjg)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%f_rimef) .GT. 1 ) THEN 
DO k = ckds,(ckde-1)
xv(ixv+1+k-ckds)= intermediate_grid%f_rimef(pig,k,pjg)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%ctopo) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%ctopo(pig,pjg)
ENDIF
IF ( SIZE(grid%ctopo2) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%ctopo2(pig,pjg)
ENDIF
IF ( SIZE(grid%winfo) .GT. 1 ) THEN 
DO k = ckds,ckde
xv(ixv+1+k-ckds)= intermediate_grid%winfo(pig,pjg,k)
ENDDO
ixv=ixv+ckde-ckds+1
ENDIF
IF ( SIZE(grid%iinfo) .GT. 1 ) THEN 
DO k = ckds,ckde
xv(ixv+1+k-ckds)= intermediate_grid%iinfo(pig,pjg,k)
ENDDO
ixv=ixv+ckde-ckds+1
ENDIF
if(interp_mp .eqv. .true.) then
DO itrace =  PARAM_FIRST_SCALAR, num_moist
DO k = ckds,(ckde-1)
xv(ixv+1+k-ckds)= intermediate_grid%moist(pig,pjg,k,itrace)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDDO
endif
if(interp_mp .eqv. .true.) then
DO itrace =  PARAM_FIRST_SCALAR, num_scalar
DO k = ckds,(ckde-1)
xv(ixv+1+k-ckds)= intermediate_grid%scalar(pig,pjg,k,itrace)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDDO
endif
if(interp_mp .eqv. .true.) then
DO itrace =  PARAM_FIRST_SCALAR, num_dfi_scalar
DO k = ckds,(ckde-1)
xv(ixv+1+k-ckds)= intermediate_grid%dfi_scalar(pig,pjg,k,itrace)
ENDDO
ixv=ixv+(ckde-1)-ckds+1
ENDDO
endif
IF ( SIZE(grid%isnowxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%isnowxy(pig,pjg)
ENDIF
IF ( SIZE(grid%tvxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%tvxy(pig,pjg)
ENDIF
IF ( SIZE(grid%tgxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%tgxy(pig,pjg)
ENDIF
IF ( SIZE(grid%canicexy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%canicexy(pig,pjg)
ENDIF
IF ( SIZE(grid%canliqxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%canliqxy(pig,pjg)
ENDIF
IF ( SIZE(grid%eahxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%eahxy(pig,pjg)
ENDIF
IF ( SIZE(grid%tahxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%tahxy(pig,pjg)
ENDIF
IF ( SIZE(grid%cmxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%cmxy(pig,pjg)
ENDIF
IF ( SIZE(grid%chxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%chxy(pig,pjg)
ENDIF
IF ( SIZE(grid%fwetxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%fwetxy(pig,pjg)
ENDIF
IF ( SIZE(grid%sneqvoxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%sneqvoxy(pig,pjg)
ENDIF
IF ( SIZE(grid%alboldxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%alboldxy(pig,pjg)
ENDIF
IF ( SIZE(grid%qsnowxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%qsnowxy(pig,pjg)
ENDIF
IF ( SIZE(grid%wslakexy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%wslakexy(pig,pjg)
ENDIF
IF ( SIZE(grid%zwtxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%zwtxy(pig,pjg)
ENDIF
IF ( SIZE(grid%waxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%waxy(pig,pjg)
ENDIF
IF ( SIZE(grid%wtxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%wtxy(pig,pjg)
ENDIF
IF ( SIZE(grid%tsnoxy) .GT. 1 ) THEN 
DO k = 1,config_flags%num_snow_layers
xv(ixv+1+k-1)= intermediate_grid%tsnoxy(pig,k,pjg)
ENDDO
ixv=ixv+config_flags%num_snow_layers-1+1
ENDIF
IF ( SIZE(grid%zsnsoxy) .GT. 1 ) THEN 
DO k = 1,config_flags%num_snso_layers
xv(ixv+1+k-1)= intermediate_grid%zsnsoxy(pig,k,pjg)
ENDDO
ixv=ixv+config_flags%num_snso_layers-1+1
ENDIF
IF ( SIZE(grid%snicexy) .GT. 1 ) THEN 
DO k = 1,config_flags%num_snow_layers
xv(ixv+1+k-1)= intermediate_grid%snicexy(pig,k,pjg)
ENDDO
ixv=ixv+config_flags%num_snow_layers-1+1
ENDIF
IF ( SIZE(grid%snliqxy) .GT. 1 ) THEN 
DO k = 1,config_flags%num_snow_layers
xv(ixv+1+k-1)= intermediate_grid%snliqxy(pig,k,pjg)
ENDDO
ixv=ixv+config_flags%num_snow_layers-1+1
ENDIF
IF ( SIZE(grid%lfmassxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%lfmassxy(pig,pjg)
ENDIF
IF ( SIZE(grid%rtmassxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%rtmassxy(pig,pjg)
ENDIF
IF ( SIZE(grid%stmassxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%stmassxy(pig,pjg)
ENDIF
IF ( SIZE(grid%woodxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%woodxy(pig,pjg)
ENDIF
IF ( SIZE(grid%stblcpxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%stblcpxy(pig,pjg)
ENDIF
IF ( SIZE(grid%fastcpxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%fastcpxy(pig,pjg)
ENDIF
IF ( SIZE(grid%xsaixy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%xsaixy(pig,pjg)
ENDIF
IF ( SIZE(grid%t2mvxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%t2mvxy(pig,pjg)
ENDIF
IF ( SIZE(grid%t2mbxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%t2mbxy(pig,pjg)
ENDIF
IF ( SIZE(grid%q2mvxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%q2mvxy(pig,pjg)
ENDIF
IF ( SIZE(grid%q2mbxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%q2mbxy(pig,pjg)
ENDIF
IF ( SIZE(grid%tradxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%tradxy(pig,pjg)
ENDIF
IF ( SIZE(grid%neexy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%neexy(pig,pjg)
ENDIF
IF ( SIZE(grid%gppxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%gppxy(pig,pjg)
ENDIF
IF ( SIZE(grid%nppxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%nppxy(pig,pjg)
ENDIF
IF ( SIZE(grid%fvegxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%fvegxy(pig,pjg)
ENDIF
IF ( SIZE(grid%qinxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%qinxy(pig,pjg)
ENDIF
IF ( SIZE(grid%runsfxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%runsfxy(pig,pjg)
ENDIF
IF ( SIZE(grid%runsbxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%runsbxy(pig,pjg)
ENDIF
IF ( SIZE(grid%ecanxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%ecanxy(pig,pjg)
ENDIF
IF ( SIZE(grid%edirxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%edirxy(pig,pjg)
ENDIF
IF ( SIZE(grid%etranxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%etranxy(pig,pjg)
ENDIF
IF ( SIZE(grid%fsaxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%fsaxy(pig,pjg)
ENDIF
IF ( SIZE(grid%firaxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%firaxy(pig,pjg)
ENDIF
IF ( SIZE(grid%aparxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%aparxy(pig,pjg)
ENDIF
IF ( SIZE(grid%psnxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%psnxy(pig,pjg)
ENDIF
IF ( SIZE(grid%savxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%savxy(pig,pjg)
ENDIF
IF ( SIZE(grid%sagxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%sagxy(pig,pjg)
ENDIF
IF ( SIZE(grid%rssunxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%rssunxy(pig,pjg)
ENDIF
IF ( SIZE(grid%rsshaxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%rsshaxy(pig,pjg)
ENDIF
IF ( SIZE(grid%bgapxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%bgapxy(pig,pjg)
ENDIF
IF ( SIZE(grid%wgapxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%wgapxy(pig,pjg)
ENDIF
IF ( SIZE(grid%tgvxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%tgvxy(pig,pjg)
ENDIF
IF ( SIZE(grid%tgbxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%tgbxy(pig,pjg)
ENDIF
IF ( SIZE(grid%chvxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%chvxy(pig,pjg)
ENDIF
IF ( SIZE(grid%chbxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%chbxy(pig,pjg)
ENDIF
IF ( SIZE(grid%shgxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%shgxy(pig,pjg)
ENDIF
IF ( SIZE(grid%shcxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%shcxy(pig,pjg)
ENDIF
IF ( SIZE(grid%shbxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%shbxy(pig,pjg)
ENDIF
IF ( SIZE(grid%evgxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%evgxy(pig,pjg)
ENDIF
IF ( SIZE(grid%evbxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%evbxy(pig,pjg)
ENDIF
IF ( SIZE(grid%ghvxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%ghvxy(pig,pjg)
ENDIF
IF ( SIZE(grid%ghbxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%ghbxy(pig,pjg)
ENDIF
IF ( SIZE(grid%irgxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%irgxy(pig,pjg)
ENDIF
IF ( SIZE(grid%ircxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%ircxy(pig,pjg)
ENDIF
IF ( SIZE(grid%irbxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%irbxy(pig,pjg)
ENDIF
IF ( SIZE(grid%trxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%trxy(pig,pjg)
ENDIF
IF ( SIZE(grid%evcxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%evcxy(pig,pjg)
ENDIF
IF ( SIZE(grid%chleafxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%chleafxy(pig,pjg)
ENDIF
IF ( SIZE(grid%chucxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%chucxy(pig,pjg)
ENDIF
IF ( SIZE(grid%chv2xy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%chv2xy(pig,pjg)
ENDIF
IF ( SIZE(grid%chb2xy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%chb2xy(pig,pjg)
ENDIF
IF ( SIZE(grid%chstarxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%chstarxy(pig,pjg)
ENDIF
IF ( SIZE(grid%soilcomp) .GT. 1 ) THEN 
DO k = 1,8
xv(ixv+1+k-1)= intermediate_grid%soilcomp(pig,k,pjg)
ENDDO
ixv=ixv+8-1+1
ENDIF
IF ( SIZE(grid%soilcl1) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%soilcl1(pig,pjg)
ENDIF
IF ( SIZE(grid%soilcl2) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%soilcl2(pig,pjg)
ENDIF
IF ( SIZE(grid%soilcl3) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%soilcl3(pig,pjg)
ENDIF
IF ( SIZE(grid%soilcl4) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%soilcl4(pig,pjg)
ENDIF
IF ( SIZE(grid%grainxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%grainxy(pig,pjg)
ENDIF
IF ( SIZE(grid%gddxy) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%gddxy(pig,pjg)
ENDIF
IF ( SIZE(grid%croptype) .GT. 1 ) THEN 
DO k = 1,5
xv(ixv+1+k-1)= intermediate_grid%croptype(pig,k,pjg)
ENDDO
ixv=ixv+5-1+1
ENDIF
IF ( SIZE(grid%planting) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%planting(pig,pjg)
ENDIF
IF ( SIZE(grid%harvest) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%harvest(pig,pjg)
ENDIF
IF ( SIZE(grid%season_gdd) .GT. 1 ) THEN 
ixv=ixv+1
xv(ixv)= intermediate_grid%season_gdd(pig,pjg)
ENDIF
IF ( SIZE(grid%gecros_state) .GT. 1 ) THEN 
DO k = 1,60
xv(ixv+1+k-1)= intermediate_grid%gecros_state(pig,k,pjg)
ENDDO
ixv=ixv+60-1+1
ENDIF
CALL rsl_lite_to_child_msg(ixv*4,xv)
CALL rsl_lite_to_parent_info( msize*4                               &
                        ,cips,cipe,cjps,cjpe                               &
                        ,nids,nide,njds,njde                               &
                        ,nest_task_offsets(ngrid%id)                      &
                        ,nest_pes_x(parent_grid%id)                            &
                        ,nest_pes_y(parent_grid%id)                            &
                        ,nest_pes_x(intermediate_grid%id)                 &
                        ,nest_pes_y(intermediate_grid%id)                 &
                        ,thisdomain_max_halo_width                         &
                        ,icoord,jcoord                                     &
                        ,idim_cd,jdim_cd                                   &
                        ,pig,pjg,retval )
ENDDO

      grid => parent_grid
    CALL pop_communicators_for_domain
END IF





      
      IF ( intercomm_active( grid%id ) ) THEN        
        local_comm = mpi_comm_to_kid( which_kid(ngrid%id), grid%id )
        ioffset = nest_task_offsets(ngrid%id)
      ELSE IF ( intercomm_active( ngrid%id ) ) THEN  
        local_comm = mpi_comm_to_mom( ngrid%id )
        ioffset = nest_task_offsets(ngrid%id)
      END IF

      IF ( grid%active_this_task .OR. ngrid%active_this_task ) THEN
        CALL mpi_comm_rank(local_comm,myproc,ierr)
        CALL mpi_comm_size(local_comm,nproc,ierr)
        CALL rsl_lite_merge_msgs( myproc, nest_pes_x(grid%id)*nest_pes_y(grid%id),         &
                                          nest_pes_x(ngrid%id)*nest_pes_y(ngrid%id),       &
                                          ioffset, local_comm )
      END IF

IF ( grid%active_this_task ) THEN
    CALL push_communicators_for_domain( grid%id )







CALL rsl_lite_from_child_info(pig,pjg,retval)
DO while ( retval .eq. 1 )
feedback_flag=cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. )
feedback_flag_v=cd_feedback_mask_v( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. )
CALL rsl_lite_from_child_msg(-1,xv) 
ixv=0 
IF ( SIZE(grid%lakedepth2d) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%lakedepth2d(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%lakedepth2d(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%savedtke12d) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%savedtke12d(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%savedtke12d(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%snowdp2d) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%snowdp2d(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%snowdp2d(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%h2osno2d) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%h2osno2d(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%h2osno2d(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%snl2d) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%snl2d(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%snl2d(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%t_grnd2d) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%t_grnd2d(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%t_grnd2d(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%t_lake3d) .GT. 1 ) THEN 
IF(feedback_flag) THEN
DO k = 1,10
grid%t_lake3d(pig,k,pjg) = nest_influence*(xv(ixv+1+k-1)) + (1.0-nest_influence)*(grid%t_lake3d(pig,k,pjg))
ENDDO
ENDIF
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%lake_icefrac3d) .GT. 1 ) THEN 
IF(feedback_flag) THEN
DO k = 1,10
grid%lake_icefrac3d(pig,k,pjg) = nest_influence*(xv(ixv+1+k-1)) + (1.0-nest_influence)*(grid%lake_icefrac3d(pig,k,pjg))
ENDDO
ENDIF
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%z_lake3d) .GT. 1 ) THEN 
IF(feedback_flag) THEN
DO k = 1,10
grid%z_lake3d(pig,k,pjg) = nest_influence*(xv(ixv+1+k-1)) + (1.0-nest_influence)*(grid%z_lake3d(pig,k,pjg))
ENDDO
ENDIF
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%dz_lake3d) .GT. 1 ) THEN 
IF(feedback_flag) THEN
DO k = 1,10
grid%dz_lake3d(pig,k,pjg) = nest_influence*(xv(ixv+1+k-1)) + (1.0-nest_influence)*(grid%dz_lake3d(pig,k,pjg))
ENDDO
ENDIF
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%t_soisno3d) .GT. 1 ) THEN 
IF(feedback_flag) THEN
DO k = 1,15
grid%t_soisno3d(pig,k,pjg) = nest_influence*(xv(ixv+1+k-1)) + (1.0-nest_influence)*(grid%t_soisno3d(pig,k,pjg))
ENDDO
ENDIF
ixv=ixv+15-1+1
ENDIF
IF ( SIZE(grid%h2osoi_ice3d) .GT. 1 ) THEN 
IF(feedback_flag) THEN
DO k = 1,15
grid%h2osoi_ice3d(pig,k,pjg) = nest_influence*(xv(ixv+1+k-1)) + (1.0-nest_influence)*(grid%h2osoi_ice3d(pig,k,pjg))
ENDDO
ENDIF
ixv=ixv+15-1+1
ENDIF
IF ( SIZE(grid%h2osoi_liq3d) .GT. 1 ) THEN 
IF(feedback_flag) THEN
DO k = 1,15
grid%h2osoi_liq3d(pig,k,pjg) = nest_influence*(xv(ixv+1+k-1)) + (1.0-nest_influence)*(grid%h2osoi_liq3d(pig,k,pjg))
ENDDO
ENDIF
ixv=ixv+15-1+1
ENDIF
IF ( SIZE(grid%h2osoi_vol3d) .GT. 1 ) THEN 
IF(feedback_flag) THEN
DO k = 1,15
grid%h2osoi_vol3d(pig,k,pjg) = nest_influence*(xv(ixv+1+k-1)) + (1.0-nest_influence)*(grid%h2osoi_vol3d(pig,k,pjg))
ENDDO
ENDIF
ixv=ixv+15-1+1
ENDIF
IF ( SIZE(grid%z3d) .GT. 1 ) THEN 
IF(feedback_flag) THEN
DO k = 1,15
grid%z3d(pig,k,pjg) = nest_influence*(xv(ixv+1+k-1)) + (1.0-nest_influence)*(grid%z3d(pig,k,pjg))
ENDDO
ENDIF
ixv=ixv+15-1+1
ENDIF
IF ( SIZE(grid%dz3d) .GT. 1 ) THEN 
IF(feedback_flag) THEN
DO k = 1,15
grid%dz3d(pig,k,pjg) = nest_influence*(xv(ixv+1+k-1)) + (1.0-nest_influence)*(grid%dz3d(pig,k,pjg))
ENDDO
ENDIF
ixv=ixv+15-1+1
ENDIF
IF ( SIZE(grid%zi3d) .GT. 1 ) THEN 
IF(feedback_flag) THEN
DO k = 1,16
grid%zi3d(pig,k,pjg) = nest_influence*(xv(ixv+1+k-1)) + (1.0-nest_influence)*(grid%zi3d(pig,k,pjg))
ENDDO
ENDIF
ixv=ixv+16-1+1
ENDIF
IF ( SIZE(grid%watsat3d) .GT. 1 ) THEN 
IF(feedback_flag) THEN
DO k = 1,10
grid%watsat3d(pig,k,pjg) = nest_influence*(xv(ixv+1+k-1)) + (1.0-nest_influence)*(grid%watsat3d(pig,k,pjg))
ENDDO
ENDIF
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%csol3d) .GT. 1 ) THEN 
IF(feedback_flag) THEN
DO k = 1,10
grid%csol3d(pig,k,pjg) = nest_influence*(xv(ixv+1+k-1)) + (1.0-nest_influence)*(grid%csol3d(pig,k,pjg))
ENDDO
ENDIF
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%tkmg3d) .GT. 1 ) THEN 
IF(feedback_flag) THEN
DO k = 1,10
grid%tkmg3d(pig,k,pjg) = nest_influence*(xv(ixv+1+k-1)) + (1.0-nest_influence)*(grid%tkmg3d(pig,k,pjg))
ENDDO
ENDIF
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%tkdry3d) .GT. 1 ) THEN 
IF(feedback_flag) THEN
DO k = 1,10
grid%tkdry3d(pig,k,pjg) = nest_influence*(xv(ixv+1+k-1)) + (1.0-nest_influence)*(grid%tkdry3d(pig,k,pjg))
ENDDO
ENDIF
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%tksatu3d) .GT. 1 ) THEN 
IF(feedback_flag) THEN
DO k = 1,10
grid%tksatu3d(pig,k,pjg) = nest_influence*(xv(ixv+1+k-1)) + (1.0-nest_influence)*(grid%tksatu3d(pig,k,pjg))
ENDDO
ENDIF
ixv=ixv+10-1+1
ENDIF
IF ( SIZE(grid%lu_index) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%lu_index(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%lu_index(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%lu_mask) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%lu_mask(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%lu_mask(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%cuprecip_swath) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%cuprecip_swath(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%cuprecip_swath(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%precip_swath) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%precip_swath(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%precip_swath(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%windsq_swath) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%windsq_swath(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%windsq_swath(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%suswind) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%suswind(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%suswind(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%suswind_swath) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%suswind_swath(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%suswind_swath(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%pdyn_smooth) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%pdyn_smooth(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%pdyn_smooth(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%pdyn_parent) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%pdyn_parent(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%pdyn_parent(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%tg_max_m10wind) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%tg_max_m10wind(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%tg_max_m10wind(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%tg_max_wwind) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%tg_max_wwind(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%tg_max_wwind(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%tg_min_wwind) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%tg_min_wwind(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%tg_min_wwind(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%tg_max_zhel_25) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%tg_max_zhel_25(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%tg_max_zhel_25(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%tg_min_zhel_25) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%tg_min_zhel_25(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%tg_min_zhel_25(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%tg_max_zhel_03) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%tg_max_zhel_03(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%tg_max_zhel_03(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%tg_min_zhel_03) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%tg_min_zhel_03(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%tg_min_zhel_03(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%tg_updhel25) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%tg_updhel25(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%tg_updhel25(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%tg_max_updhel25) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%tg_max_updhel25(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%tg_max_updhel25(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%tg_updhel03) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%tg_updhel03(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%tg_updhel03(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%tg_max_updhel03) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%tg_max_updhel03(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%tg_max_updhel03(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%tg_total_precip) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%tg_total_precip(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%tg_total_precip(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%hres_fis) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%hres_fis(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%hres_fis(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%hres_lnd) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%hres_lnd(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%hres_lnd(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%hres_lu_index) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%hres_lu_index(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%hres_lu_index(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%sm) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%sm(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%sm(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%pd) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%pd(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%pd(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%fis) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%fis(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%fis(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%t) .GT. 1 ) THEN 
IF(feedback_flag) THEN
DO k = ckds,(ckde-1)
grid%t(pig,pjg,k) = nest_influence*(xv(ixv+1+k-ckds)) + (1.0-nest_influence)*(grid%t(pig,pjg,k))
ENDDO
ENDIF
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%q) .GT. 1 ) THEN 
IF(feedback_flag) THEN
DO k = ckds,(ckde-1)
grid%q(pig,pjg,k) = nest_influence*(xv(ixv+1+k-ckds)) + (1.0-nest_influence)*(grid%q(pig,pjg,k))
ENDDO
ENDIF
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%u) .GT. 1 ) THEN 
IF(feedback_flag_v) THEN
DO k = ckds,(ckde-1)
grid%u(pig,pjg,k) = nest_influence*(xv(ixv+1+k-ckds)) + (1.0-nest_influence)*(grid%u(pig,pjg,k))
ENDDO
ENDIF
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%v) .GT. 1 ) THEN 
IF(feedback_flag_v) THEN
DO k = ckds,(ckde-1)
grid%v(pig,pjg,k) = nest_influence*(xv(ixv+1+k-ckds)) + (1.0-nest_influence)*(grid%v(pig,pjg,k))
ENDDO
ENDIF
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%q2) .GT. 1 ) THEN 
IF(feedback_flag) THEN
DO k = ckds,(ckde-1)
grid%q2(pig,pjg,k) = nest_influence*(xv(ixv+1+k-ckds)) + (1.0-nest_influence)*(grid%q2(pig,pjg,k))
ENDDO
ENDIF
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%acswupt) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%acswupt(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%acswupt(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%acswuptc) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%acswuptc(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%acswuptc(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%acswdnt) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%acswdnt(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%acswdnt(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%acswdntc) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%acswdntc(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%acswdntc(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%acswupb) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%acswupb(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%acswupb(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%acswupbc) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%acswupbc(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%acswupbc(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%acswdnb) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%acswdnb(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%acswdnb(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%acswdnbc) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%acswdnbc(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%acswdnbc(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%aclwupt) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%aclwupt(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%aclwupt(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%aclwuptc) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%aclwuptc(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%aclwuptc(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%aclwdnt) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%aclwdnt(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%aclwdnt(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%aclwdntc) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%aclwdntc(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%aclwdntc(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%aclwupb) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%aclwupb(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%aclwupb(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%aclwupbc) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%aclwupbc(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%aclwupbc(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%aclwdnb) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%aclwdnb(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%aclwdnb(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%aclwdnbc) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%aclwdnbc(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%aclwdnbc(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%swupt) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%swupt(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%swupt(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%swuptc) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%swuptc(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%swuptc(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%swdnt) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%swdnt(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%swdnt(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%swdntc) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%swdntc(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%swdntc(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%swupb) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%swupb(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%swupb(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%swupbc) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%swupbc(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%swupbc(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%swdnb) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%swdnb(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%swdnb(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%swdnbc) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%swdnbc(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%swdnbc(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%lwupt) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%lwupt(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%lwupt(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%lwuptc) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%lwuptc(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%lwuptc(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%lwdnt) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%lwdnt(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%lwdnt(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%lwdntc) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%lwdntc(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%lwdntc(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%lwupb) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%lwupb(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%lwupb(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%lwupbc) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%lwupbc(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%lwupbc(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%lwdnb) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%lwdnb(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%lwdnb(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%lwdnbc) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%lwdnbc(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%lwdnbc(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%qnwfa2d) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%qnwfa2d(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%qnwfa2d(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%cwm) .GT. 1 ) THEN 
IF(feedback_flag) THEN
DO k = ckds,(ckde-1)
grid%cwm(pig,pjg,k) = nest_influence*(xv(ixv+1+k-ckds)) + (1.0-nest_influence)*(grid%cwm(pig,pjg,k))
ENDDO
ENDIF
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%f_ice) .GT. 1 ) THEN 
IF(feedback_flag) THEN
DO k = ckds,(ckde-1)
grid%f_ice(pig,k,pjg) = nest_influence*(xv(ixv+1+k-ckds)) + (1.0-nest_influence)*(grid%f_ice(pig,k,pjg))
ENDDO
ENDIF
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%f_rain) .GT. 1 ) THEN 
IF(feedback_flag) THEN
DO k = ckds,(ckde-1)
grid%f_rain(pig,k,pjg) = nest_influence*(xv(ixv+1+k-ckds)) + (1.0-nest_influence)*(grid%f_rain(pig,k,pjg))
ENDDO
ENDIF
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%f_rimef) .GT. 1 ) THEN 
IF(feedback_flag) THEN
DO k = ckds,(ckde-1)
grid%f_rimef(pig,k,pjg) = nest_influence*(xv(ixv+1+k-ckds)) + (1.0-nest_influence)*(grid%f_rimef(pig,k,pjg))
ENDDO
ENDIF
ixv=ixv+(ckde-1)-ckds+1
ENDIF
IF ( SIZE(grid%ctopo) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%ctopo(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%ctopo(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%ctopo2) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%ctopo2(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%ctopo2(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%winfo) .GT. 1 ) THEN 
IF(feedback_flag) THEN
DO k = ckds,ckde
grid%winfo(pig,pjg,k) = nest_influence*(xv(ixv+1+k-ckds)) + (1.0-nest_influence)*(grid%winfo(pig,pjg,k))
ENDDO
ENDIF
ixv=ixv+ckde-ckds+1
ENDIF
IF ( SIZE(grid%iinfo) .GT. 1 ) THEN 
IF(feedback_flag) THEN
DO k = ckds,ckde
grid%iinfo(pig,pjg,k) = nest_influence*(xv(ixv+1+k-ckds)) + (1.0-nest_influence)*(grid%iinfo(pig,pjg,k))
ENDDO
ENDIF
ixv=ixv+ckde-ckds+1
ENDIF
if(interp_mp .eqv. .true.) then
DO itrace =  PARAM_FIRST_SCALAR, num_moist
IF(feedback_flag) THEN
DO k = ckds,(ckde-1)
moist(pig,pjg,k,itrace) = nest_influence*(xv(ixv+1+k-ckds)) + (1.0-nest_influence)*(moist(pig,pjg,k,itrace))
ENDDO
ENDIF
ixv=ixv+(ckde-1)-ckds+1
ENDDO
endif
if(interp_mp .eqv. .true.) then
DO itrace =  PARAM_FIRST_SCALAR, num_scalar
IF(feedback_flag) THEN
DO k = ckds,(ckde-1)
scalar(pig,pjg,k,itrace) = nest_influence*(xv(ixv+1+k-ckds)) + (1.0-nest_influence)*(scalar(pig,pjg,k,itrace))
ENDDO
ENDIF
ixv=ixv+(ckde-1)-ckds+1
ENDDO
endif
if(interp_mp .eqv. .true.) then
DO itrace =  PARAM_FIRST_SCALAR, num_dfi_scalar
IF(feedback_flag) THEN
DO k = ckds,(ckde-1)
dfi_scalar(pig,pjg,k,itrace) = nest_influence*(xv(ixv+1+k-ckds)) + (1.0-nest_influence)*(dfi_scalar(pig,pjg,k,itrace))
ENDDO
ENDIF
ixv=ixv+(ckde-1)-ckds+1
ENDDO
endif
IF ( SIZE(grid%isnowxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%isnowxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%isnowxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%tvxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%tvxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%tvxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%tgxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%tgxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%tgxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%canicexy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%canicexy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%canicexy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%canliqxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%canliqxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%canliqxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%eahxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%eahxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%eahxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%tahxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%tahxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%tahxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%cmxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%cmxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%cmxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%chxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%chxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%chxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%fwetxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%fwetxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%fwetxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%sneqvoxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%sneqvoxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%sneqvoxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%alboldxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%alboldxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%alboldxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%qsnowxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%qsnowxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%qsnowxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%wslakexy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%wslakexy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%wslakexy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%zwtxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%zwtxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%zwtxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%waxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%waxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%waxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%wtxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%wtxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%wtxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%tsnoxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
DO k = 1,config_flags%num_snow_layers
grid%tsnoxy(pig,k,pjg) = nest_influence*(xv(ixv+1+k-1)) + (1.0-nest_influence)*(grid%tsnoxy(pig,k,pjg))
ENDDO
ENDIF
ixv=ixv+config_flags%num_snow_layers-1+1
ENDIF
IF ( SIZE(grid%zsnsoxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
DO k = 1,config_flags%num_snso_layers
grid%zsnsoxy(pig,k,pjg) = nest_influence*(xv(ixv+1+k-1)) + (1.0-nest_influence)*(grid%zsnsoxy(pig,k,pjg))
ENDDO
ENDIF
ixv=ixv+config_flags%num_snso_layers-1+1
ENDIF
IF ( SIZE(grid%snicexy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
DO k = 1,config_flags%num_snow_layers
grid%snicexy(pig,k,pjg) = nest_influence*(xv(ixv+1+k-1)) + (1.0-nest_influence)*(grid%snicexy(pig,k,pjg))
ENDDO
ENDIF
ixv=ixv+config_flags%num_snow_layers-1+1
ENDIF
IF ( SIZE(grid%snliqxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
DO k = 1,config_flags%num_snow_layers
grid%snliqxy(pig,k,pjg) = nest_influence*(xv(ixv+1+k-1)) + (1.0-nest_influence)*(grid%snliqxy(pig,k,pjg))
ENDDO
ENDIF
ixv=ixv+config_flags%num_snow_layers-1+1
ENDIF
IF ( SIZE(grid%lfmassxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%lfmassxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%lfmassxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%rtmassxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%rtmassxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%rtmassxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%stmassxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%stmassxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%stmassxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%woodxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%woodxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%woodxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%stblcpxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%stblcpxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%stblcpxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%fastcpxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%fastcpxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%fastcpxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%xsaixy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%xsaixy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%xsaixy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%t2mvxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%t2mvxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%t2mvxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%t2mbxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%t2mbxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%t2mbxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%q2mvxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%q2mvxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%q2mvxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%q2mbxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%q2mbxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%q2mbxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%tradxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%tradxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%tradxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%neexy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%neexy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%neexy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%gppxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%gppxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%gppxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%nppxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%nppxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%nppxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%fvegxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%fvegxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%fvegxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%qinxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%qinxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%qinxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%runsfxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%runsfxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%runsfxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%runsbxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%runsbxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%runsbxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%ecanxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%ecanxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%ecanxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%edirxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%edirxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%edirxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%etranxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%etranxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%etranxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%fsaxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%fsaxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%fsaxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%firaxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%firaxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%firaxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%aparxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%aparxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%aparxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%psnxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%psnxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%psnxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%savxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%savxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%savxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%sagxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%sagxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%sagxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%rssunxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%rssunxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%rssunxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%rsshaxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%rsshaxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%rsshaxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%bgapxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%bgapxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%bgapxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%wgapxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%wgapxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%wgapxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%tgvxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%tgvxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%tgvxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%tgbxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%tgbxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%tgbxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%chvxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%chvxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%chvxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%chbxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%chbxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%chbxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%shgxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%shgxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%shgxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%shcxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%shcxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%shcxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%shbxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%shbxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%shbxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%evgxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%evgxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%evgxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%evbxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%evbxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%evbxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%ghvxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%ghvxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%ghvxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%ghbxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%ghbxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%ghbxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%irgxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%irgxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%irgxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%ircxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%ircxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%ircxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%irbxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%irbxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%irbxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%trxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%trxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%trxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%evcxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%evcxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%evcxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%chleafxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%chleafxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%chleafxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%chucxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%chucxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%chucxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%chv2xy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%chv2xy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%chv2xy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%chb2xy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%chb2xy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%chb2xy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%chstarxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%chstarxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%chstarxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%soilcomp) .GT. 1 ) THEN 
IF(feedback_flag) THEN
DO k = 1,8
grid%soilcomp(pig,k,pjg) = nest_influence*(xv(ixv+1+k-1)) + (1.0-nest_influence)*(grid%soilcomp(pig,k,pjg))
ENDDO
ENDIF
ixv=ixv+8-1+1
ENDIF
IF ( SIZE(grid%soilcl1) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%soilcl1(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%soilcl1(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%soilcl2) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%soilcl2(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%soilcl2(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%soilcl3) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%soilcl3(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%soilcl3(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%soilcl4) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%soilcl4(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%soilcl4(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%grainxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%grainxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%grainxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%gddxy) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%gddxy(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%gddxy(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%croptype) .GT. 1 ) THEN 
IF(feedback_flag) THEN
DO k = 1,5
grid%croptype(pig,k,pjg) = nest_influence*(xv(ixv+1+k-1)) + (1.0-nest_influence)*(grid%croptype(pig,k,pjg))
ENDDO
ENDIF
ixv=ixv+5-1+1
ENDIF
IF ( SIZE(grid%planting) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%planting(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%planting(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%harvest) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%harvest(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%harvest(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%season_gdd) .GT. 1 ) THEN 
IF(feedback_flag) THEN
ixv=ixv+1
grid%season_gdd(pig,pjg) = nest_influence*(xv(ixv)) + (1.0-nest_influence)*(grid%season_gdd(pig,pjg))
ELSE
ixv=ixv+1
ENDIF
ENDIF
IF ( SIZE(grid%gecros_state) .GT. 1 ) THEN 
IF(feedback_flag) THEN
DO k = 1,60
grid%gecros_state(pig,k,pjg) = nest_influence*(xv(ixv+1+k-1)) + (1.0-nest_influence)*(grid%gecros_state(pig,k,pjg))
ENDDO
ENDIF
ixv=ixv+60-1+1
ENDIF
CALL rsl_lite_from_child_info(pig,pjg,retval)
ENDDO



      

      CALL get_ijk_from_grid (  ngrid,                                 &
                                nids, nide, njds, njde, nkds, nkde,    &
                                nims, nime, njms, njme, nkms, nkme,    &
                                nips, nipe, njps, njpe, nkps, nkpe     )
      CALL get_ijk_from_grid (  grid ,              &
                                ids, ide, jds, jde, kds, kde,    &
                                ims, ime, jms, jme, kms, kme,    &
                                ips, ipe, jps, jpe, kps, kpe    )

     before_smooth_halo: if(config_flags%smooth_option/=0) then






CALL HALO_INTERP_UP_sub ( grid, &
  config_flags, &
  num_moist, &
  moist, &
  num_scalar, &
  scalar, &
  num_dfi_scalar, &
  dfi_scalar, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )

     endif before_smooth_halo

      CALL get_ijk_from_grid (  grid ,                   &
                                cids, cide, cjds, cjde, ckds, ckde,    &
                                cims, cime, cjms, cjme, ckms, ckme,    &
                                cips, cipe, cjps, cjpe, ckps, ckpe    )

     smoother: if(config_flags%smooth_option/=0) then






DO itrace = PARAM_FIRST_SCALAR, num_dfi_scalar
IF ( SIZE( dfi_scalar, 1 ) * SIZE( dfi_scalar, 2 ) .GT. 1 .and. (interp_mp .eqv. .true.) ) THEN 
CALL nmm_smoother_ijk (  &         
                  dfi_scalar(grid%sm31,grid%sm32,grid%sm33,itrace),   &       
                 cids, cide, ckds, ckde, cjds, cjde,   &         
                 cims, cime, ckms, ckme, cjms, cjme,   &         
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         
                 nids, nide, nkds, nkde, njds, njde,   &         
                 nims, nime, nkms, nkme, njms, njme,   &         
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         
                  .FALSE., .FALSE.,                                                &         
                  ngrid%i_parent_start, ngrid%j_parent_start,                     &
                  ngrid%parent_grid_ratio, ngrid%parent_grid_ratio                &
                  ) 
ENDIF
ENDDO

     endif smoother

    CALL pop_communicators_for_domain
END IF


      RETURN
   END SUBROUTINE feedback_domain_nmm_part2







   SUBROUTINE wrf_gatherv_real (Field, field_ofst,            &
                                my_count ,                    &    
                                globbuf, glob_ofst ,          &    
                                counts                      , &    
                                displs                      , &    
                                root                        , &    
                                communicator                , &    
                                ierr )
   USE module_dm, ONLY : getrealmpitype
   IMPLICIT NONE
   INTEGER field_ofst, glob_ofst
   INTEGER my_count, communicator, root, ierr
   INTEGER , DIMENSION(*) :: counts, displs
   REAL, DIMENSION(*) :: Field, globbuf
   INCLUDE 'mpif.h'

           CALL mpi_gatherv( Field( field_ofst ),      &    
                            my_count ,                       &    
                            getrealmpitype() ,               &    
                            globbuf( glob_ofst ) ,                 &    
                            counts                         , &    
                            displs                         , &    
                            getrealmpitype()               , &    
                            root                           , &    
                            communicator                   , &    
                            ierr )

   END SUBROUTINE wrf_gatherv_real

   SUBROUTINE wrf_gatherv_double (Field, field_ofst,            &
                                my_count ,                    &    
                                globbuf, glob_ofst ,          &    
                                counts                      , &    
                                displs                      , &    
                                root                        , &    
                                communicator                , &    
                                ierr )

   IMPLICIT NONE
   INTEGER field_ofst, glob_ofst
   INTEGER my_count, communicator, root, ierr
   INTEGER , DIMENSION(*) :: counts, displs




   REAL, DIMENSION(*) :: Field, globbuf
   INCLUDE 'mpif.h'

           CALL mpi_gatherv( Field( field_ofst ),      &    
                            my_count ,                       &    
                            MPI_DOUBLE_PRECISION         ,               &    
                            globbuf( glob_ofst ) ,                 &    
                            counts                         , &    
                            displs                         , &    
                            MPI_DOUBLE_PRECISION                       , &    
                            root                           , &    
                            communicator                   , &    
                            ierr )

   END SUBROUTINE wrf_gatherv_double

   SUBROUTINE wrf_gatherv_integer (Field, field_ofst,            &
                                my_count ,                    &    
                                globbuf, glob_ofst ,          &    
                                counts                      , &    
                                displs                      , &    
                                root                        , &    
                                communicator                , &    
                                ierr )
   IMPLICIT NONE
   INTEGER field_ofst, glob_ofst
   INTEGER my_count, communicator, root, ierr
   INTEGER , DIMENSION(*) :: counts, displs
   INTEGER, DIMENSION(*) :: Field, globbuf
   INCLUDE 'mpif.h'

           CALL mpi_gatherv( Field( field_ofst ),      &    
                            my_count ,                       &    
                            MPI_INTEGER         ,               &    
                            globbuf( glob_ofst ) ,                 &    
                            counts                         , &    
                            displs                         , &    
                            MPI_INTEGER                       , &    
                            root                           , &    
                            communicator                   , &    
                            ierr )

   END SUBROUTINE wrf_gatherv_integer


   SUBROUTINE wrf_scatterv_real (                             &
                                globbuf, glob_ofst ,          &    
                                counts                      , &    
                                Field, field_ofst,            &
                                my_count ,                    &    
                                displs                      , &    
                                root                        , &    
                                communicator                , &    
                                ierr )
   USE module_dm, ONLY : getrealmpitype
   IMPLICIT NONE
   INTEGER field_ofst, glob_ofst
   INTEGER my_count, communicator, root, ierr
   INTEGER , DIMENSION(*) :: counts, displs
   REAL, DIMENSION(*) :: Field, globbuf
   INCLUDE 'mpif.h'

           CALL mpi_scatterv(                                &
                            globbuf( glob_ofst ) ,           &    
                            counts                         , &    
                            displs                         , &    
                            getrealmpitype()               , &    
                            Field( field_ofst ),             &    
                            my_count ,                       &    
                            getrealmpitype() ,               &    
                            root                           , &    
                            communicator                   , &    
                            ierr )

   END SUBROUTINE wrf_scatterv_real

   SUBROUTINE wrf_scatterv_double (                           &
                                globbuf, glob_ofst ,          &    
                                counts                      , &    
                                Field, field_ofst,            &
                                my_count ,                    &    
                                displs                      , &    
                                root                        , &    
                                communicator                , &    
                                ierr )
   IMPLICIT NONE
   INTEGER field_ofst, glob_ofst
   INTEGER my_count, communicator, root, ierr
   INTEGER , DIMENSION(*) :: counts, displs
   REAL, DIMENSION(*) :: Field, globbuf
   INCLUDE 'mpif.h'





           CALL mpi_scatterv(                                &
                            globbuf( glob_ofst ) ,           &    
                            counts                         , &    
                            displs                         , &    
                            MPI_DOUBLE_PRECISION           , &    
                            Field( field_ofst ),             &    
                            my_count ,                       &    
                            MPI_DOUBLE_PRECISION         ,   &    
                            root                           , &    
                            communicator                   , &    
                            ierr )

   END SUBROUTINE wrf_scatterv_double

   SUBROUTINE wrf_scatterv_integer (                          &
                                globbuf, glob_ofst ,          &    
                                counts                      , &    
                                Field, field_ofst,            &
                                my_count ,                    &    
                                displs                      , &    
                                root                        , &    
                                communicator                , &    
                                ierr )
   IMPLICIT NONE
   INTEGER field_ofst, glob_ofst
   INTEGER my_count, communicator, root, ierr
   INTEGER , DIMENSION(*) :: counts, displs
   INTEGER, DIMENSION(*) :: Field, globbuf
   INCLUDE 'mpif.h'

           CALL mpi_scatterv(                                &
                            globbuf( glob_ofst ) ,           &    
                            counts                         , &    
                            displs                         , &    
                            MPI_INTEGER                    , &    
                            Field( field_ofst ),             &    
                            my_count ,                       &    
                            MPI_INTEGER         ,            &    
                            root                           , &    
                            communicator                   , &    
                            ierr )

   END SUBROUTINE wrf_scatterv_integer


     SUBROUTINE wrf_dm_gatherv ( v, elemsize , km_s, km_e, wordsz )
      IMPLICIT NONE
      INTEGER  elemsize, km_s, km_e, wordsz
      REAL v(*)
      IF ( wordsz .EQ. 8 ) THEN
         CALL wrf_dm_gatherv_double(v, elemsize , km_s, km_e)
      ELSE
         CALL wrf_dm_gatherv_single(v, elemsize , km_s, km_e)
      END IF
     END SUBROUTINE wrf_dm_gatherv

     SUBROUTINE wrf_dm_gatherv_double ( v, elemsize , km_s, km_e )
      IMPLICIT NONE
      INTEGER  elemsize, km_s, km_e
      REAL*8 v(0:*)
      REAL*8 v_local((km_e-km_s+1)*elemsize)
      INTEGER, DIMENSION(:), ALLOCATABLE :: recvcounts, displs
      INTEGER send_type, myproc, nproc, local_comm, ierr, i
   INCLUDE 'mpif.h'
      send_type = MPI_DOUBLE_PRECISION
      CALL wrf_get_dm_communicator ( local_comm )
      CALL wrf_get_nproc( nproc )
      CALL wrf_get_myproc( myproc )
      ALLOCATE( recvcounts(nproc), displs(nproc) )
      i = (km_e-km_s+1)*elemsize
      CALL mpi_allgather( i,1,MPI_INTEGER,recvcounts,1,MPI_INTEGER,local_comm,ierr) ;
      i = (km_s)*elemsize
      CALL mpi_allgather( i,1,MPI_INTEGER,displs,1,MPI_INTEGER,local_comm,ierr) ;
      DO i = 1,elemsize*(km_e-km_s+1)
        v_local(i) = v(i+elemsize*km_s-1)
      END DO
      CALL mpi_allgatherv( v_local,                                       &
                           (km_e-km_s+1)*elemsize,                        &
                           send_type,                                     &
                           v,                                             &
                           recvcounts,                                    &
                           displs,                                        &
                           send_type,                                     &
                           local_comm,                                    &
                           ierr )
      DEALLOCATE(recvcounts)
      DEALLOCATE(displs)
      return
     END SUBROUTINE wrf_dm_gatherv_double

     SUBROUTINE wrf_dm_gatherv_single ( v, elemsize , km_s, km_e )
      IMPLICIT NONE
      INTEGER  elemsize, km_s, km_e
      REAL*4 v(0:*)
      REAL*4 v_local((km_e-km_s+1)*elemsize)
      INTEGER, DIMENSION(:), ALLOCATABLE :: recvcounts, displs
      INTEGER send_type, myproc, nproc, local_comm, ierr, i
   INCLUDE 'mpif.h'
      send_type = MPI_REAL
      CALL wrf_get_dm_communicator ( local_comm )
      CALL wrf_get_nproc( nproc )
      CALL wrf_get_myproc( myproc )
      ALLOCATE( recvcounts(nproc), displs(nproc) )
      i = (km_e-km_s+1)*elemsize
      CALL mpi_allgather( i,1,MPI_INTEGER,recvcounts,1,MPI_INTEGER,local_comm,ierr) ;
      i = (km_s)*elemsize
      CALL mpi_allgather( i,1,MPI_INTEGER,displs,1,MPI_INTEGER,local_comm,ierr) ;
      DO i = 1,elemsize*(km_e-km_s+1)
        v_local(i) = v(i+elemsize*km_s-1)
      END DO
      CALL mpi_allgatherv( v_local,                                       &
                           (km_e-km_s+1)*elemsize,                        &
                           send_type,                                     &
                           v,                                             &
                           recvcounts,                                    &
                           displs,                                        &
                           send_type,                                     &
                           local_comm,                                    &
                           ierr )
      DEALLOCATE(recvcounts)
      DEALLOCATE(displs)
      return
     END SUBROUTINE wrf_dm_gatherv_single

      SUBROUTINE wrf_dm_decomp1d( nt, km_s, km_e )
       IMPLICIT NONE
       INTEGER, INTENT(IN)  :: nt
       INTEGER, INTENT(OUT) :: km_s, km_e
     
       INTEGER nn, nnp,  na, nb
       INTEGER myproc, nproc

       CALL wrf_get_myproc(myproc)
       CALL wrf_get_nproc(nproc)
       nn = nt / nproc           
       nnp = nn
       if ( myproc .lt. mod( nt, nproc ) )   nnp = nnp + 1 

       na = min( myproc, mod(nt,nproc) ) 
       nb = max( 0, myproc - na )        
       km_s = na * ( nn+1) + nb * nn     
       km_e = km_s + nnp - 1             
      END SUBROUTINE wrf_dm_decomp1d


SUBROUTINE wrf_dm_define_comms ( grid )
   USE module_domain, ONLY : domain
   IMPLICIT NONE
   TYPE(domain) , INTENT (INOUT) :: grid
   RETURN
END SUBROUTINE wrf_dm_define_comms

SUBROUTINE tfp_message( fname, lno )
   CHARACTER*(*) fname
   INTEGER lno
   CHARACTER*1024 mess
   WRITE(mess,*)'tfp_message: ',trim(fname),lno
   CALL wrf_message(mess)
     CALL wrf_error_fatal3("<stdin>",23202,&
mess)
END SUBROUTINE tfp_message

   SUBROUTINE set_dm_debug
    USE module_dm, ONLY : dm_debug_flag
    IMPLICIT NONE
    dm_debug_flag = .TRUE.
   END SUBROUTINE set_dm_debug
   SUBROUTINE reset_dm_debug
    USE module_dm, ONLY : dm_debug_flag
    IMPLICIT NONE
    dm_debug_flag = .FALSE.
   END SUBROUTINE reset_dm_debug
   SUBROUTINE get_dm_debug ( arg )
    USE module_dm, ONLY : dm_debug_flag
    IMPLICIT NONE
    LOGICAL arg
    arg = dm_debug_flag
   END SUBROUTINE get_dm_debug

