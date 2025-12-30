






MODULE module_wrf_quilt



















































  USE module_internal_header_util
  USE module_timing
  USE module_driver_constants 
  USE module_dm, ONLY : current_id
  USE module_cpl, ONLY : coupler_on, cpl_set_dm_communicator, cpl_finalize

  INTEGER, PARAMETER :: int_num_handles = 99
  INTEGER, PARAMETER :: max_servers = int_num_handles+1  
  LOGICAL, DIMENSION(0:int_num_handles) :: okay_to_write, int_handle_in_use, okay_to_commit
  INTEGER, DIMENSION(0:int_num_handles) :: int_num_bytes_to_write, io_form
  INTEGER, DIMENSION(0:int_num_handles) :: which_grid_is_handle
  INTEGER, DIMENSION(0:int_num_handles) :: prev_server_for_handle
  REAL, POINTER,SAVE :: int_local_output_buffer(:)
  INTEGER,      SAVE :: int_local_output_cursor
  LOGICAL          :: quilting_enabled
  LOGICAL          :: disable_quilt = .FALSE.
  INTEGER          :: server_for_handle(int_num_handles,max_domains)
  INTEGER          :: reduced(2), reduced_dummy(2)
  LOGICAL, EXTERNAL :: wrf_dm_on_monitor

  INTEGER :: mpi_comm_avail(max_domains),availrank(max_domains)
  LOGICAL :: in_avail=.false., poll_servers=.false.

  INTEGER, ALLOCATABLE :: role_for_task(:)
  INTEGER nio_groups
  INTEGER :: mpi_comm_local, mpi_comm_local_io_server_tmp
  LOGICAL :: compute_group_master(max_servers,max_domains)
  INTEGER :: mpi_comm_io_groups(max_servers,max_domains)
  INTEGER :: nio_tasks_per_group(max_domains)
  INTEGER :: ntasks
  INTEGER :: mytask

  INTEGER, PARAMETER           :: onebyte = 1
  INTEGER comm_io_servers, iserver, hdrbufsize, obufsize
  INTEGER, DIMENSION(4096)     :: hdrbuf
  INTEGER, DIMENSION(int_num_handles)     :: handle


  CONTAINS

    SUBROUTINE mpi_x_comm_size(i,j,k)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: i
        INTEGER, INTENT(OUT) :: j, k
        call Mpi_Comm_Size ( i,j,k )
    END SUBROUTINE mpi_x_comm_size

    INTEGER FUNCTION get_server_id ( dhandle )









      IMPLICIT NONE
      INTEGER, INTENT(IN) :: dhandle

      INTEGER :: id
      
      IF ( dhandle .GE. 1 .AND. dhandle .LE. int_num_handles ) THEN
        id = which_grid_is_handle(dhandle)
        IF ( id .LT. 1 .OR. id .GT. max_domains ) THEN
          CALL wrf_error_fatal3("<stdin>",125,&
"module_io_quilt: get_server_id, bad grid id stored with handle")
        ENDIF
        IF ( server_for_handle ( dhandle, id ) .GE. 1 ) THEN
          get_server_id = server_for_handle ( dhandle, id )
        ELSE
           IF(poll_servers) THEN
              
              call wrf_quilt_find_server(server_for_handle(dhandle, id))
           ELSE
              
              prev_server_for_handle(id) = mod ( prev_server_for_handle(id) + 1 , nio_groups )
              server_for_handle( dhandle, id ) = prev_server_for_handle(id)+1
           ENDIF
           get_server_id=server_for_handle(dhandle, id)
        ENDIF
      ELSE
         CALL wrf_message('module_io_quilt: get_server_id bad dhandle' )
      ENDIF
    END FUNCTION get_server_id

    SUBROUTINE set_server_id ( dhandle, value )
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: dhandle, value
       INTEGER id
       IF ( dhandle .GE. 1 .AND. dhandle .LE. int_num_handles ) THEN
         id = which_grid_is_handle(dhandle)
         IF ( id .GE. 1 .AND. ID .LE. max_domains ) THEN
           server_for_handle(dhandle,id) = value
         ELSE
           CALL wrf_message('module_io_quilt: set_server_id bad grid id stored with handle' )
         ENDIF
       ELSE
         CALL wrf_message('module_io_quilt: set_server_id bad dhandle' )
       ENDIF
    END SUBROUTINE set_server_id

    LOGICAL FUNCTION get_poll_servers() 
      implicit none
      get_poll_servers=poll_servers
    end FUNCTION get_poll_servers

    SUBROUTINE int_get_fresh_handle( retval )










      INTEGER, INTENT(OUT) :: retval
      INTEGER i
      retval = -1
      DO i = 1, int_num_handles
        IF ( .NOT. int_handle_in_use(i) )  THEN
          retval = i
          GOTO 33
        ENDIF
      ENDDO
33    CONTINUE
      IF ( retval < 0 )  THEN
        CALL wrf_error_fatal3("<stdin>",189,&
"frame/module_io_quilt.F: int_get_fresh_handle() can not")
      ENDIF
      int_handle_in_use(i) = .TRUE.
      NULLIFY ( int_local_output_buffer )
    END SUBROUTINE int_get_fresh_handle

    SUBROUTINE setup_quilt_servers ( id, nio_tasks_per_group, &
                                     role_for_task,           &
                                     num_io_tasks,            &
                                     ncompute_tasks,          &
                                     mytask,                  &
                                     ntasks,                  &
                                     n_groups_arg,            &
                                     mpi_comm_wrld,           &
                                     mpi_comm_local,          &
                                     mpi_comm_io_groups,      &
                                     compute_node   )




























































      USE module_configure
      USE module_dm, ONLY : compute_mesh,nest_pes_x,nest_pes_y,domain_active_this_task,&
                            tasks_per_split,comm_start,dm_task_split
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER,                      INTENT(IN)  :: id, nio_tasks_per_group(:), mytask, ntasks, &
                                                   n_groups_arg, mpi_comm_wrld
      INTEGER,  INTENT(IN)                      :: role_for_task(:)
      INTEGER,  INTENT(IN)                      :: ncompute_tasks, num_io_tasks
      INTEGER,  INTENT(OUT)                     :: mpi_comm_local
      INTEGER, DIMENSION(100,max_domains),      INTENT(OUT) :: mpi_comm_io_groups
      LOGICAL, INTENT(OUT)   :: compute_node

      INTEGER                     :: i, j, ii, comdup, ierr, niotasks, n_groups, iisize, itask, nio
      INTEGER, DIMENSION(ntasks)  :: icolor
      CHARACTER*128 mess
      INTEGER :: io_form_setting
      INTEGER :: me
      INTEGER :: k, m, nprocx, nprocy, found, found_io
      LOGICAL :: reorder_mesh
      CHARACTER*256 message


      n_groups = n_groups_arg
      IF ( n_groups .LT. 1 ) n_groups = 1









      nio = nio_tasks_per_group(id)
      IF ( nio .LT. 0 ) THEN
        nio = 0
      ENDIF
      IF ( nio .EQ. 0 ) THEN
        quilting_enabled = .FALSE.
        compute_node = .TRUE.
        mpi_comm_local = mpi_comm_wrld
        mpi_comm_io_groups(:,id) = mpi_comm_wrld
        RETURN
      ENDIF
      found = comm_start(id)
      found_io = -99
      DO i=1,ntasks
         IF ( role_for_task(i) .EQ. id ) THEN
           found_io = i-1         
           exit
         ENDIF
      ENDDO
      IF ( found_io .eq. -99 ) THEN
        CALL wrf_error_fatal3("<stdin>",321,&
'setup_quilt_servers: found_io not found, should not happen (but there it is)')
      ENDIF

      quilting_enabled = .TRUE.

      icolor = -99  
      ii = 0
      DO i = 1, ntasks
        IF ( role_for_task( i ) .EQ. id ) THEN
          icolor(i) = ii / nio
          ii = ii + 1
        ENDIF
      ENDDO
      CALL MPI_Comm_dup(mpi_comm_wrld,comdup,ierr)
      CALL MPI_Comm_split(comdup,icolor(mytask+1),mytask,mpi_comm_local,ierr)
      IF ( icolor(mytask+1) .NE. -99 ) THEN
          mpi_comm_local_io_server_tmp = mpi_comm_local
      ENDIF





      nprocx = nest_pes_x(id)
      nprocy = nest_pes_y(id)

      IF ( nio .GT. nprocy*nprocx ) THEN
        CALL wrf_error_fatal3("<stdin>",349,&
'more io tasks than compute tasks specified' ) 
      ENDIF
  
      m = mod(nprocy*nprocx,nio)  
      ii = 0
      j = 1
      do while ( j .le. nprocx*nprocy )
        do i = 0, (nprocx*nprocy)/nio+min(m,1)-1
           icolor(j+comm_start(id)) = ii
           j = j + 1
        enddo
        ii = ii + 1
        m = max(m-1,0)
      enddo


      DO j = 1, n_groups
        DO i = found_io+1,ntasks
          icolor(i) = -99
        ENDDO
        ii = 0
        DO i = found_io+(j-1)*nio+1,found_io+j*nio
          icolor(i) = ii
          ii = ii+1
        ENDDO
        CALL MPI_Comm_dup(mpi_comm_wrld,comdup,ierr)
        CALL MPI_Comm_split(comdup,icolor(mytask+1),mytask, &
                            mpi_comm_io_groups(j,id),ierr)
      ENDDO

      if(poll_servers) then
         poll_servers=.false.
         call wrf_message('Warning: server polling does not work with pnetcdf_quilt.  Disabled poll_servers.')
      else
         if(n_groups==1) then
            poll_servers=.false.
            call wrf_message('Server polling is useless with one io group.  Disabled poll_servers.')
         endif
      endif


      if(poll_servers) then
         
         
         
         
         

         

         call mpi_comm_rank(mpi_comm_wrld,me,ierr)

         icolor=-99
         in_avail=.false.

         IF ( found .eq. me ) THEN  
           in_avail=.true.
         ENDIF
         icolor(1)=1

         nio = nio_tasks_per_group(id)
         do j=1,n_groups
            i=ncompute_tasks+j*nio-1
            if(me+1==i) then
               in_avail=.true. 
            endif
            icolor(i)=1
         enddo

         CALL MPI_Comm_dup(mpi_comm_wrld,comdup,ierr)
         CALL MPI_Comm_split(comdup,icolor(me+1),me, &
                             mpi_comm_avail(id),ierr)

         availrank=MPI_UNDEFINED
         IF(in_avail) THEN
            call mpi_comm_rank(mpi_comm_avail(id),availrank(id),ierr)
         ENDIF
         IF ( role_for_task(me+1) .GT. 1000 ) THEN   
            mpi_comm_avail(1) = mpi_comm_avail(id)
            availrank(1) = availrank(id)
         ENDIF

      endif




      DO j = 1, n_groups

        IF ((found .LE. mytask .AND. mytask .LT.  found +ncompute_tasks ).OR. &    
             (found_io+(j-1)*nio .LE. mytask .AND. mytask .LT. found_io+j*nio) &    
            ) THEN

         CALL MPI_Comm_Size( mpi_comm_io_groups(j,id) , iisize, ierr )
         
         
         CALL MPI_Comm_Rank( mpi_comm_io_groups(j,id) , me , ierr )


         
         

         IF ( (found .LE. mytask .AND. mytask .LT.  ncompute_tasks ) ) THEN
            compute_node = .TRUE.
            
            
            
            compute_group_master(j,id) = (me .EQ. 0)
         ELSE
           IF (found_io+(j-1)*nio .LE. mytask .AND. mytask .LT. found_io+j*nio) THEN
             mpi_comm_io_groups(1,1) = mpi_comm_io_groups(j,id)
           ENDIF
         ENDIF
        ENDIF

      ENDDO

    END SUBROUTINE setup_quilt_servers

    SUBROUTINE sokay ( stream, io_form )
    USE module_state_description
    CHARACTER*(*) stream
    CHARACTER*256 mess
    INTEGER io_form

    SELECT CASE (io_form)
      CASE ( IO_NETCDF   )
         RETURN
      CASE ( IO_INTIO   )
         RETURN
      CASE ( IO_GRIB1 )
         RETURN
      CASE ( IO_GRIB2 )
         RETURN
      CASE (0)
         RETURN
      CASE DEFAULT
         WRITE(mess,*)' An output format has been specified that is incompatible with quilting: io_form: ',io_form,' ',TRIM(stream)
         CALL wrf_error_fatal3("<stdin>",488,&
mess)
    END SELECT
    END SUBROUTINE sokay

    SUBROUTINE quilt













      USE module_state_description
      USE module_quilt_outbuf_ops
      USE module_configure, only : grid_config_rec_type, model_config_rec, model_to_grid_config_rec
      IMPLICIT NONE
      INCLUDE 'mpif.h'
  INTEGER, PARAMETER ::  int_ioexit			=  	     10
  INTEGER, PARAMETER ::  int_open_for_write_begin	=  	     20
  INTEGER, PARAMETER ::  int_open_for_write_commit	=  	     30
  INTEGER, PARAMETER ::  int_open_for_read 		=  	     40
  INTEGER, PARAMETER ::  int_inquire_opened 		=  	     60
  INTEGER, PARAMETER ::  int_inquire_filename 		=  	     70
  INTEGER, PARAMETER ::  int_iosync 			=  	     80
  INTEGER, PARAMETER ::  int_ioclose 			=  	     90
  INTEGER, PARAMETER ::  int_next_time 			=  	    100
  INTEGER, PARAMETER ::  int_set_time 			=  	    110
  INTEGER, PARAMETER ::  int_next_var 			=  	    120
  INTEGER, PARAMETER ::  int_dom_ti_real 		=  	    140
  INTEGER, PARAMETER ::  int_dom_ti_double 		=  	    160
  INTEGER, PARAMETER ::  int_dom_ti_integer 		=  	    180
  INTEGER, PARAMETER ::  int_dom_ti_logical 		=  	    200
  INTEGER, PARAMETER ::  int_dom_ti_char 		=  	    220
  INTEGER, PARAMETER ::  int_dom_td_real 		=  	    240
  INTEGER, PARAMETER ::  int_dom_td_double 		=  	    260
  INTEGER, PARAMETER ::  int_dom_td_integer 		=  	    280
  INTEGER, PARAMETER ::  int_dom_td_logical 		=  	    300
  INTEGER, PARAMETER ::  int_dom_td_char 		=  	    320
  INTEGER, PARAMETER ::  int_var_ti_real 		=  	    340
  INTEGER, PARAMETER ::  int_var_ti_double 		=  	    360
  INTEGER, PARAMETER ::  int_var_ti_integer 		=  	    380
  INTEGER, PARAMETER ::  int_var_ti_logical 		=  	    400
  INTEGER, PARAMETER ::  int_var_ti_char 		=  	    420
  INTEGER, PARAMETER ::  int_var_td_real 		=  	    440
  INTEGER, PARAMETER ::  int_var_td_double 		=  	    460
  INTEGER, PARAMETER ::  int_var_td_integer 		=  	    480
  INTEGER, PARAMETER ::  int_var_td_logical 		=  	    500
  INTEGER, PARAMETER ::  int_var_td_char 		=  	    520
  INTEGER, PARAMETER ::  int_field 			=  	    530
  INTEGER, PARAMETER ::  int_var_info 			=  	    540
  INTEGER, PARAMETER ::  int_noop 			=  	    550
      integer, parameter  :: WRF_FILE_NOT_OPENED                  = 100
      integer, parameter  :: WRF_FILE_OPENED_NOT_COMMITTED        = 101
      integer, parameter  :: WRF_FILE_OPENED_FOR_WRITE            = 102
      integer, parameter  :: WRF_FILE_OPENED_FOR_READ             = 103
      integer, parameter  :: WRF_REAL                             = 104
      integer, parameter  :: WRF_DOUBLE                           = 105
      integer, parameter  :: WRF_FLOAT=WRF_REAL
      integer, parameter  :: WRF_INTEGER                          = 106
      integer, parameter  :: WRF_LOGICAL                          = 107
      integer, parameter  :: WRF_COMPLEX                          = 108
      integer, parameter  :: WRF_DOUBLE_COMPLEX                   = 109
      integer, parameter  :: WRF_FILE_OPENED_FOR_UPDATE           = 110


      integer, parameter  :: WRF_FILE_OPENED_AND_COMMITTED        = 102
      TYPE (grid_config_rec_type)  :: config_flags
      INTEGER itag, ninbuf, ntasks_io_group, ntasks_local_group, mytask_local, ierr
      INTEGER istat
      INTEGER mytask_io_group
      INTEGER   :: nout_set = 0
      INTEGER   :: obufsize, bigbufsize, chunksize, sz
      REAL, DIMENSION(1)      :: dummy
      INTEGER, ALLOCATABLE, DIMENSION(:) :: obuf, bigbuf
      REAL,    ALLOCATABLE, DIMENSION(:) :: RDATA
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IDATA
      CHARACTER (LEN=512) :: CDATA
      CHARACTER (LEN=80) :: fname
      INTEGER icurs, hdrbufsize, itypesize, ftypesize, rtypesize, Status, fstat, io_form_arg
      INTEGER :: DataHandle, FieldType, Comm, IOComm, DomainDesc, code, Count
      INTEGER, DIMENSION(3) :: DomainStart , DomainEnd , MemoryStart , MemoryEnd , PatchStart , PatchEnd
      INTEGER :: dummybuf(1)
      INTEGER :: num_noops, num_commit_messages, num_field_training_msgs, hdr_tag
      CHARACTER (len=256) :: DateStr , Element, VarName, MemoryOrder , Stagger , DimNames(3), FileName, SysDepInfo, mess
      INTEGER, EXTERNAL :: use_package
      LOGICAL           :: stored_write_record, retval
      INTEGER iii, jjj, vid, dom_id
      LOGICAL           :: call_server_ready

logical okay_to_w
character*120 sysline

      dom_id = 1 
      CALL model_to_grid_config_rec ( dom_id , model_config_rec , config_flags )







      CALL quilt_pnc()
      RETURN


      SysDepInfo = " "
      if ( config_flags%use_netcdf_classic ) SysDepInfo="use_netcdf_classic"
      CALL ext_ncd_ioinit( SysDepInfo, ierr )
      SysDepInfo = " "
      CALL ext_int_ioinit( SysDepInfo, ierr )
      CALL ext_gr1_ioinit( SysDepInfo, ierr)
      CALL ext_gr2_ioinit( SysDepInfo, ierr)

      call_server_ready = .true. 

      okay_to_commit = .false.
      stored_write_record = .false.
      ninbuf = 0
      
      
      
      
      
      CALL mpi_x_comm_size( mpi_comm_io_groups(1,1), ntasks_io_group,    ierr )
      CALL MPI_COMM_RANK( mpi_comm_io_groups(1,1), mytask_io_group,    ierr )
      CALL mpi_x_comm_size( mpi_comm_local,        ntasks_local_group, ierr )
      CALL MPI_COMM_RANK( mpi_comm_local,        mytask_local,       ierr )

      CALL MPI_TYPE_SIZE( MPI_INTEGER, itypesize, ierr )
      IF ( itypesize <= 0 ) THEN
        CALL wrf_error_fatal3("<stdin>",625,&
"external/RSL/module_dm.F: quilt: type size <= 0 invalid")
      ENDIF










okay_to_w = .false.
      DO WHILE (.TRUE.)  














         if(poll_servers .and. call_server_ready) then
            call_server_ready=.false.

            
            call wrf_quilt_server_ready()
         endif


        

        
        ! if needed (currently needed only for ioclose).
        reduced_dummy = 0


        CALL MPI_Reduce( reduced_dummy, reduced, 2, MPI_INTEGER, MPI_SUM, mytask_io_group, mpi_comm_io_groups(1,1), ierr )
        obufsize = reduced(1)




        IF ( obufsize .LT. 0 ) THEN
          IF ( obufsize .EQ. -100 ) THEN         
            CALL ext_ncd_ioexit( Status )
            CALL ext_int_ioexit( Status )
            CALL ext_gr1_ioexit( Status )
            CALL ext_gr2_ioexit( Status )
            CALL wrf_message ( 'I/O QUILT SERVERS DONE' )
            IF (coupler_on) THEN 
               CALL cpl_finalize() 
            ELSE
               CALL mpi_finalize(ierr)
            END IF
            STOP
          ELSE
            WRITE(mess,*)'Possible 32-bit overflow on output server. Try larger nio_tasks_per_group in namelist.'
            CALL wrf_error_fatal3("<stdin>",690,&
mess)
          ENDIF
        ENDIF







        IF ( obufsize .GT. 0 ) THEN
          ALLOCATE( obuf( (obufsize+1)/itypesize ) )


          CALL collect_on_comm_debug("module_io_quilt.G",753, mpi_comm_io_groups(1,1),        &
                                onebyte,                      &
                                dummy, 0,                     &
                                obuf, obufsize )

        ELSE
          ! Necessarily, the compute processes send the ioclose signal,
          
          ! will stall on the ioclose message waiting for the quilt 
          
          
          
          
          
          ! Then a header representing the ioclose message is constructed
          
          
          
          
          
          
          ALLOCATE( obuf( 4096 ) )
          
          CALL int_gen_handle_header( obuf, obufsize, itypesize, &
                                      reduced(2) , int_ioclose )

          if(poll_servers) then 


             call_server_ready=.true.
          endif
        ENDIF











        CALL init_store_piece_of_field
        CALL mpi_type_size ( MPI_INTEGER , itypesize , ierr )



        vid = 0
        icurs = itypesize
        num_noops = 0 
        num_commit_messages = 0 
        num_field_training_msgs = 0 
        DO WHILE ( icurs .lt. obufsize ) 
          hdr_tag = get_hdr_tag( obuf ( icurs / itypesize ) )
          SELECT CASE ( hdr_tag )
            CASE ( int_field )
              CALL int_get_write_field_header ( obuf(icurs/itypesize), hdrbufsize, itypesize, ftypesize,  &
                                                DataHandle , DateStr , VarName , Dummy , FieldType , Comm , IOComm, &
                                                DomainDesc , MemoryOrder , Stagger , DimNames ,              &
                                                DomainStart , DomainEnd ,                                    &
                                                MemoryStart , MemoryEnd ,                                    &
                                                PatchStart , PatchEnd )
              chunksize = (PatchEnd(1)-PatchStart(1)+1)*(PatchEnd(2)-PatchStart(2)+1)* &
                          (PatchEnd(3)-PatchStart(3)+1)*ftypesize

              IF ( DomainDesc .EQ. 333933 ) THEN  
                 IF ( num_field_training_msgs .EQ. 0 ) THEN
                   call add_to_bufsize_for_field( VarName, hdrbufsize )

                 ENDIF
                 num_field_training_msgs = num_field_training_msgs + 1
              ELSE
                 call add_to_bufsize_for_field( VarName, hdrbufsize )

              ENDIF
              icurs = icurs + hdrbufsize



              
              
              IF ( DomainDesc .NE. 333933 ) THEN   

                call add_to_bufsize_for_field( VarName, chunksize )
                icurs = icurs + chunksize
              ENDIF
            CASE ( int_open_for_write_commit )  
              hdrbufsize = obuf(icurs/itypesize)
              IF (num_commit_messages.EQ.0) THEN
                call add_to_bufsize_for_field( 'COMMIT', hdrbufsize )
              ENDIF
              num_commit_messages = num_commit_messages + 1
              icurs = icurs + hdrbufsize
            CASE DEFAULT
              hdrbufsize = obuf(icurs/itypesize)




































              IF ((hdr_tag.EQ.int_noop.AND.mytask_local.NE.0.AND.num_noops.LE.0)  &
                  .OR.hdr_tag.NE.int_noop) THEN
                write(VarName,'(I5.5)')vid 

                call add_to_bufsize_for_field( VarName, hdrbufsize )
                vid = vid+1
              ENDIF
              IF ( hdr_tag .EQ. int_noop ) num_noops = num_noops + 1
              icurs = icurs + hdrbufsize
          END SELECT
        ENDDO 



        vid = 0
        icurs = itypesize
        num_noops = 0 
        num_commit_messages = 0 
        num_field_training_msgs = 0 
        DO WHILE ( icurs .lt. obufsize ) 

          hdr_tag = get_hdr_tag( obuf ( icurs / itypesize ) )
          SELECT CASE ( hdr_tag )
            CASE ( int_field )
              CALL int_get_write_field_header ( obuf(icurs/itypesize), hdrbufsize, itypesize, ftypesize,  &
                                                DataHandle , DateStr , VarName , Dummy , FieldType , Comm , IOComm, &
                                                DomainDesc , MemoryOrder , Stagger , DimNames ,              &
                                                DomainStart , DomainEnd ,                                    &
                                                MemoryStart , MemoryEnd ,                                    &
                                                PatchStart , PatchEnd )
              chunksize = (PatchEnd(1)-PatchStart(1)+1)*(PatchEnd(2)-PatchStart(2)+1)* &
                          (PatchEnd(3)-PatchStart(3)+1)*ftypesize

              IF ( DomainDesc .EQ. 333933 ) THEN  
                 IF ( num_field_training_msgs .EQ. 0 ) THEN
                   call store_piece_of_field( obuf(icurs/itypesize), VarName, hdrbufsize )

                 ENDIF
                 num_field_training_msgs = num_field_training_msgs + 1
              ELSE
                 call store_piece_of_field( obuf(icurs/itypesize), VarName, hdrbufsize )

              ENDIF
              icurs = icurs + hdrbufsize
              
              
              IF ( DomainDesc .NE. 333933 ) THEN   

                call store_piece_of_field( obuf(icurs/itypesize), VarName, chunksize )
                icurs = icurs + chunksize
              ENDIF
            CASE ( int_open_for_write_commit )  
              hdrbufsize = obuf(icurs/itypesize)
              IF (num_commit_messages.EQ.0) THEN
                call store_piece_of_field( obuf(icurs/itypesize), 'COMMIT', hdrbufsize )
              ENDIF
              num_commit_messages = num_commit_messages + 1
              icurs = icurs + hdrbufsize
            CASE DEFAULT
              hdrbufsize = obuf(icurs/itypesize)
              IF ((hdr_tag.EQ.int_noop.AND.mytask_local.NE.0.AND.num_noops.LE.0)  &
                  .OR.hdr_tag.NE.int_noop) THEN
                write(VarName,'(I5.5)')vid 

                call store_piece_of_field( obuf(icurs/itypesize), VarName, hdrbufsize )
                vid = vid+1
              ENDIF
              IF ( hdr_tag .EQ. int_noop ) num_noops = num_noops + 1
              icurs = icurs + hdrbufsize
          END SELECT
        ENDDO 



        CALL init_retrieve_pieces_of_field


        CALL retrieve_pieces_of_field ( obuf , VarName, obufsize, sz, retval )


        CALL MPI_Reduce( sz, bigbufsize, 1, MPI_INTEGER, MPI_SUM, ntasks_local_group-1, mpi_comm_local, ierr )



        DO WHILE ( retval ) 



          IF ( mytask_local .EQ. ntasks_local_group-1 ) THEN
            ALLOCATE( bigbuf( (bigbufsize+1)/itypesize ) )
         else
            ALLOCATE( bigbuf(1) )
          ENDIF



          CALL collect_on_comm_debug2("module_io_quilt.G",983,Trim(VarName),        &
                                get_hdr_tag(obuf),sz,get_hdr_rec_size(obuf),  &
                                mpi_comm_local,                               &
                                onebyte,                                      &
                                obuf, sz,                                     &
                                bigbuf, bigbufsize )


          IF ( mytask_local .EQ. ntasks_local_group-1 ) THEN






            icurs = itypesize  

            stored_write_record = .false.


            DO WHILE ( icurs .lt. bigbufsize ) 
              CALL mpi_type_size ( MPI_INTEGER , itypesize , ierr )





              SELECT CASE ( get_hdr_tag( bigbuf(icurs/itypesize) ) )




                CASE ( int_noop )
                  CALL int_get_noop_header( bigbuf(icurs/itypesize), hdrbufsize, itypesize )
                  icurs = icurs + hdrbufsize


                CASE ( int_dom_td_real )
                  CALL mpi_type_size( MPI_REAL, ftypesize, ierr )
                  ALLOCATE( RData( bigbuf(icurs/itypesize + 4 ) ) )      
                  CALL int_get_td_header( bigbuf(icurs/itypesize:), hdrbufsize, itypesize, ftypesize, &
                                          DataHandle, DateStr, Element, RData, Count, code )
                  icurs = icurs + hdrbufsize

                  SELECT CASE (use_package(io_form(DataHandle)))
                    CASE ( IO_NETCDF   )
                      CALL ext_ncd_put_dom_td_real( handle(DataHandle),TRIM(Element),TRIM(DateStr),RData, Count, Status )
                    CASE ( IO_INTIO   )
                      CALL ext_int_put_dom_td_real( handle(DataHandle),TRIM(Element),TRIM(DateStr),RData, Count, Status )
                 CASE ( IO_GRIB1 )
                    CALL ext_gr1_put_dom_td_real( handle(DataHandle),TRIM(Element),TRIM(DateStr),RData, Count, Status )
                 CASE ( IO_GRIB2 )
                    CALL ext_gr2_put_dom_td_real( handle(DataHandle),TRIM(Element),TRIM(DateStr),RData, Count, Status )
                     CASE DEFAULT
                      Status = 0
                  END SELECT

                  DEALLOCATE( RData )

                CASE ( int_dom_ti_real )
                  CALL mpi_type_size( MPI_REAL, ftypesize, ierr )
                  ALLOCATE( RData( bigbuf(icurs/itypesize + 4 ) ) )      
                  CALL int_get_ti_header( bigbuf(icurs/itypesize:), hdrbufsize, itypesize, ftypesize, &
                                          DataHandle, Element, RData, Count, code )

                  icurs = icurs + hdrbufsize

                  SELECT CASE (use_package(io_form(DataHandle)))
                    CASE ( IO_NETCDF   )
                      CALL ext_ncd_put_dom_ti_real( handle(DataHandle),TRIM(Element), RData, Count, Status )

                    CASE ( IO_INTIO   )
                      CALL ext_int_put_dom_ti_real( handle(DataHandle),TRIM(Element), RData, Count, Status )
                 CASE ( IO_GRIB1 )
                    CALL ext_gr1_put_dom_ti_real( handle(DataHandle),TRIM(Element), RData, Count, Status )
                 CASE ( IO_GRIB2 )
                    CALL ext_gr2_put_dom_ti_real( handle(DataHandle),TRIM(Element), RData, Count, Status )
                    CASE DEFAULT
                      Status = 0
                  END SELECT

                  DEALLOCATE( RData )


                CASE ( int_dom_td_integer )
                  CALL mpi_type_size( MPI_INTEGER, ftypesize, ierr )
                  ALLOCATE( IData( bigbuf(icurs/itypesize + 4 ) ) )      
                  CALL int_get_td_header( bigbuf(icurs/itypesize:), hdrbufsize, itypesize, ftypesize, &
                                          DataHandle, DateStr, Element, IData, Count, code )

                  icurs = icurs + hdrbufsize

                  SELECT CASE (use_package(io_form(DataHandle)))
                    CASE ( IO_NETCDF   )
                      CALL ext_ncd_put_dom_td_integer( handle(DataHandle),TRIM(Element), Trim(DateStr), IData, Count, Status )
                    CASE ( IO_INTIO   )
                      CALL ext_int_put_dom_td_integer( handle(DataHandle),TRIM(Element), Trim(DateStr), IData, Count, Status )
                 CASE ( IO_GRIB1 )
                    CALL ext_gr1_put_dom_td_integer( handle(DataHandle),TRIM(Element), Trim(DateStr), IData, Count, Status )
                 CASE ( IO_GRIB2 )
                    CALL ext_gr2_put_dom_td_integer( handle(DataHandle),TRIM(Element), Trim(DateStr), IData, Count, Status )
                    CASE DEFAULT
                      Status = 0
                  END SELECT

                  DEALLOCATE( IData )


                CASE ( int_dom_ti_integer )

                  CALL mpi_type_size( MPI_INTEGER, ftypesize, ierr )
                  ALLOCATE( IData( bigbuf(icurs/itypesize + 4 ) ) )      
                  CALL int_get_ti_header( bigbuf(icurs/itypesize:), hdrbufsize, itypesize, ftypesize, &
                                          DataHandle, Element, IData, Count, code )

                  icurs = icurs + hdrbufsize
                  SELECT CASE (use_package(io_form(DataHandle)))
                    CASE ( IO_NETCDF   )
                      CALL ext_ncd_put_dom_ti_integer( handle(DataHandle),TRIM(Element), IData, Count, Status )

                    CASE ( IO_INTIO   )
                      CALL ext_int_put_dom_ti_integer( handle(DataHandle),TRIM(Element), IData, Count, Status )
                 CASE ( IO_GRIB1 )
                    CALL ext_gr1_put_dom_ti_integer( handle(DataHandle),TRIM(Element), IData, Count, Status )
                 CASE ( IO_GRIB2 )
                    CALL ext_gr2_put_dom_ti_integer( handle(DataHandle),TRIM(Element), IData, Count, Status )

                    CASE DEFAULT
                      Status = 0
                  END SELECT

                  DEALLOCATE( IData)
 

                CASE ( int_set_time )
                  CALL int_get_ti_header_char( bigbuf(icurs/itypesize), hdrbufsize, itypesize, &
                                               DataHandle, Element, VarName, CData, code )

                  SELECT CASE (use_package(io_form(DataHandle)))
                    CASE ( IO_INTIO   )
                      CALL ext_int_set_time ( handle(DataHandle), TRIM(CData), Status)
                    CASE DEFAULT
                      Status = 0
                  END SELECT

                  icurs = icurs + hdrbufsize


                CASE ( int_dom_ti_char )
                  CALL int_get_ti_header_char( bigbuf(icurs/itypesize), hdrbufsize, itypesize, &
                                               DataHandle, Element, VarName, CData, code )


                  SELECT CASE (use_package(io_form(DataHandle)))
                    CASE ( IO_NETCDF   )
                      CALL ext_ncd_put_dom_ti_char ( handle(DataHandle), TRIM(Element), TRIM(CData), Status)
                    CASE ( IO_INTIO   )
                      CALL ext_int_put_dom_ti_char ( handle(DataHandle), TRIM(Element), TRIM(CData), Status)
                 CASE ( IO_GRIB1 )
                    CALL ext_gr1_put_dom_ti_char ( handle(DataHandle), TRIM(Element), TRIM(CData), Status)
                 CASE ( IO_GRIB2 )
                    CALL ext_gr2_put_dom_ti_char ( handle(DataHandle), TRIM(Element), TRIM(CData), Status)
                    CASE DEFAULT
                      Status = 0
                  END SELECT

                  icurs = icurs + hdrbufsize


                CASE ( int_var_ti_char )
                  CALL int_get_ti_header_char( bigbuf(icurs/itypesize), hdrbufsize, itypesize, &
                                               DataHandle, Element, VarName, CData, code )


                  SELECT CASE (use_package(io_form(DataHandle)))
                    CASE ( IO_NETCDF   )
                      CALL ext_ncd_put_var_ti_char ( handle(DataHandle), TRIM(Element), TRIM(VarName), TRIM(CData), Status)
                    CASE ( IO_INTIO   )
                      CALL ext_int_put_var_ti_char ( handle(DataHandle), TRIM(Element), TRIM(VarName), TRIM(CData), Status)
                 CASE ( IO_GRIB1 )
                    CALL ext_gr1_put_var_ti_char ( handle(DataHandle), TRIM(Element), TRIM(VarName), TRIM(CData), Status)
                 CASE ( IO_GRIB2 )
                    CALL ext_gr2_put_var_ti_char ( handle(DataHandle), TRIM(Element), TRIM(VarName), TRIM(CData), Status)
                    CASE DEFAULT
                      Status = 0
                  END SELECT

                  icurs = icurs + hdrbufsize

                CASE ( int_ioexit )

                  CALL wrf_error_fatal3("<stdin>",1124,&
                         "quilt: should have handled int_ioexit already")
! The I/O server "root" handles the "ioclose" request.
                CASE ( int_ioclose )
                  CALL int_get_handle_header( bigbuf(icurs/itypesize), hdrbufsize, itypesize, &
                                              DataHandle , code )
                  icurs = icurs + hdrbufsize

                  IF ( DataHandle .GE. 1 ) THEN


                  SELECT CASE (use_package(io_form(DataHandle)))
                    CASE ( IO_NETCDF   )
                      CALL ext_ncd_inquire_filename( handle(DataHandle), fname, fstat, Status )
                      IF ( fstat .EQ. WRF_FILE_OPENED_FOR_WRITE .OR. fstat .EQ. WRF_FILE_OPENED_NOT_COMMITTED ) THEN
                        CALL ext_ncd_ioclose(handle(DataHandle),Status)
                      ENDIF
                    CASE ( IO_PNETCDF   )
                      CALL ext_pnc_inquire_filename( handle(DataHandle), fname, fstat, Status )
                      IF ( fstat .EQ. WRF_FILE_OPENED_FOR_WRITE .OR. fstat .EQ. WRF_FILE_OPENED_NOT_COMMITTED ) THEN
                        CALL ext_pnc_ioclose(handle(DataHandle),Status)
                      ENDIF
                    CASE ( IO_INTIO   )
                      CALL ext_int_inquire_filename( handle(DataHandle), fname, fstat, Status )
                      IF ( fstat .EQ. WRF_FILE_OPENED_FOR_WRITE .OR. fstat .EQ. WRF_FILE_OPENED_NOT_COMMITTED ) THEN
                        CALL ext_int_ioclose(handle(DataHandle),Status)
                      ENDIF
                 CASE ( IO_GRIB1 )
                    CALL ext_gr1_inquire_filename( handle(DataHandle), fname, fstat, Status )
                    IF ( fstat .EQ. WRF_FILE_OPENED_FOR_WRITE .OR. fstat .EQ. WRF_FILE_OPENED_NOT_COMMITTED ) THEN
                      CALL ext_gr1_ioclose(handle(DataHandle),Status)
                    ENDIF
                 CASE ( IO_GRIB2 )
                    CALL ext_gr2_inquire_filename( handle(DataHandle), fname, fstat, Status )
                    IF ( fstat .EQ. WRF_FILE_OPENED_FOR_WRITE .OR. fstat .EQ. WRF_FILE_OPENED_NOT_COMMITTED ) THEN
                      CALL ext_gr2_ioclose(handle(DataHandle),Status)
                    ENDIF
                    CASE DEFAULT
                      Status = 0
                  END SELECT
                  ENDIF



                  IF (fname(1:6) .EQ. 'wrfout' .AND. config_flags%output_ready_flag ) THEN
                    OPEN (unit=99,file='wrfoutReady' // fname(7:30), status='unknown', access='sequential')
                    CLOSE (99)
                  ENDIF


                CASE ( int_open_for_write_begin )

                  CALL int_get_ofwb_header( bigbuf(icurs/itypesize), hdrbufsize, itypesize, &
                                            FileName,SysDepInfo,io_form_arg,DataHandle )





                  icurs = icurs + hdrbufsize

                
                  io_form(DataHandle) = io_form_arg


                  SELECT CASE (use_package(io_form(DataHandle)))
                    CASE ( IO_NETCDF   )
                      CALL ext_ncd_open_for_write_begin(FileName,Comm,IOComm,SysDepInfo,handle(DataHandle),Status)

                    CASE ( IO_INTIO   )
                      CALL ext_int_open_for_write_begin(FileName,Comm,IOComm,SysDepInfo,handle(DataHandle),Status)
                    CASE ( IO_GRIB1 )
                       CALL ext_gr1_open_for_write_begin(FileName,Comm,IOComm,SysDepInfo,handle(DataHandle),Status)
                    CASE ( IO_GRIB2 )
                       CALL ext_gr2_open_for_write_begin(FileName,Comm,IOComm,SysDepInfo,handle(DataHandle),Status)
                    CASE DEFAULT
                      Status = 0
                  END SELECT
                
                  okay_to_write(DataHandle) = .false.





                CASE ( int_open_for_write_commit )

                  CALL int_get_handle_header( bigbuf(icurs/itypesize), hdrbufsize, itypesize, &
                                              DataHandle , code )
                  icurs = icurs + hdrbufsize
                  okay_to_commit(DataHandle) = .true.











                CASE ( int_field )
                  CALL mpi_type_size( MPI_INTEGER, ftypesize, ierr )
                  CALL int_get_write_field_header ( bigbuf(icurs/itypesize), hdrbufsize, itypesize, ftypesize,  &
                                                    DataHandle , DateStr , VarName , Dummy , FieldType , Comm , IOComm, &
                                                    DomainDesc , MemoryOrder , Stagger , DimNames ,              &
                                                    DomainStart , DomainEnd ,                                    &
                                                    MemoryStart , MemoryEnd ,                                    &
                                                    PatchStart , PatchEnd )

                  icurs = icurs + hdrbufsize

                  IF ( okay_to_write(DataHandle) ) THEN




                    IF ( FieldType .EQ. WRF_FLOAT .OR. FieldType .EQ. WRF_DOUBLE)  THEN
                      
                      
                      IF ( FieldType .EQ. WRF_DOUBLE)  THEN

                        CALL mpi_type_size( MPI_DOUBLE_PRECISION, ftypesize, ierr )
                      ELSE
                        CALL mpi_type_size( MPI_REAL, ftypesize, ierr )
                      ENDIF
                      stored_write_record = .true.
                      CALL store_patch_in_outbuf ( bigbuf(icurs/itypesize), dummybuf, TRIM(DateStr), TRIM(VarName) , &
                                                   FieldType, TRIM(MemoryOrder), TRIM(Stagger), DimNames, &
                                                   DomainStart , DomainEnd , &
                                                   MemoryStart , MemoryEnd , &
                                                   PatchStart , PatchEnd )

                    ELSE IF ( FieldType .EQ. WRF_INTEGER ) THEN
                      CALL mpi_type_size( MPI_INTEGER, ftypesize, ierr )
                      stored_write_record = .true.
                      CALL store_patch_in_outbuf ( dummybuf, bigbuf(icurs/itypesize), TRIM(DateStr), TRIM(VarName) , &
                                                   FieldType, TRIM(MemoryOrder), TRIM(Stagger), DimNames, &
                                                   DomainStart , DomainEnd , &
                                                   MemoryStart , MemoryEnd , &
                                                   PatchStart , PatchEnd )
                    ELSE IF ( FieldType .EQ. WRF_LOGICAL ) THEN
                      ftypesize = 4
                    ENDIF
                    icurs = icurs + (PatchEnd(1)-PatchStart(1)+1)*(PatchEnd(2)-PatchStart(2)+1)* &
                                    (PatchEnd(3)-PatchStart(3)+1)*ftypesize
                  ELSE
                    SELECT CASE (use_package(io_form(DataHandle)))
                      CASE ( IO_NETCDF   )
                        CALL ext_ncd_write_field ( handle(DataHandle) , TRIM(DateStr) ,         &
                                   TRIM(VarName) , dummy , FieldType , Comm , IOComm,           &
                                   DomainDesc , TRIM(MemoryOrder) , TRIM(Stagger) , DimNames ,  &
                                   DomainStart , DomainEnd ,                                    &
                                   DomainStart , DomainEnd ,                                    &
                                   DomainStart , DomainEnd ,                                    &
                                   Status )
                      CASE DEFAULT
                        Status = 0
                    END SELECT
                  ENDIF
                CASE ( int_iosync )
                  CALL int_get_handle_header( bigbuf(icurs/itypesize), hdrbufsize, itypesize, &
                                            DataHandle , code )
                  icurs = icurs + hdrbufsize
                CASE DEFAULT
                  WRITE(mess,*)'quilt: bad tag: ',get_hdr_tag( bigbuf(icurs/itypesize) ),' icurs ',icurs/itypesize
                  CALL wrf_error_fatal3("<stdin>",1291,&
mess )
              END SELECT

            ENDDO 



            IF (stored_write_record) THEN







              CALL write_outbuf ( handle(DataHandle), use_package(io_form(DataHandle))) 

            ENDIF




            IF (okay_to_commit(DataHandle)) THEN

              SELECT CASE (use_package(io_form(DataHandle)))
                CASE ( IO_NETCDF   )
                  CALL ext_ncd_inquire_filename( handle(DataHandle), fname, fstat, Status )
                  IF ( fstat .EQ. WRF_FILE_OPENED_NOT_COMMITTED ) THEN
                    CALL ext_ncd_open_for_write_commit(handle(DataHandle),Status)
                    okay_to_write(DataHandle) = .true.
                  ENDIF
                CASE ( IO_INTIO   )
                  CALL ext_int_inquire_filename( handle(DataHandle), fname, fstat, Status )
                  IF ( fstat .EQ. WRF_FILE_OPENED_NOT_COMMITTED ) THEN
                    CALL ext_int_open_for_write_commit(handle(DataHandle),Status)
                    okay_to_write(DataHandle) = .true.
                  ENDIF
                 CASE ( IO_GRIB1 )
                    CALL ext_gr1_inquire_filename( handle(DataHandle), fname, fstat, Status )
                    IF ( fstat .EQ. WRF_FILE_OPENED_NOT_COMMITTED ) THEN
                       CALL ext_gr1_open_for_write_commit(handle(DataHandle),Status)
                       okay_to_write(DataHandle) = .true.
                    ENDIF
                 CASE ( IO_GRIB2 )
                    CALL ext_gr2_inquire_filename( handle(DataHandle), fname, fstat, Status )
                    IF ( fstat .EQ. WRF_FILE_OPENED_NOT_COMMITTED ) THEN
                       CALL ext_gr2_open_for_write_commit(handle(DataHandle),Status)
                       okay_to_write(DataHandle) = .true.
                    ENDIF

                CASE DEFAULT
                  Status = 0
              END SELECT

            okay_to_commit(DataHandle) = .false.
          ENDIF
          DEALLOCATE( bigbuf )
        ENDIF
        if(allocated(bigbuf)) deallocate(bigbuf)


        CALL retrieve_pieces_of_field ( obuf , VarName, obufsize, sz, retval )


        CALL MPI_Reduce( sz, bigbufsize, 1, MPI_INTEGER,MPI_SUM, ntasks_local_group-1,mpi_comm_local, ierr )



      END DO 

      DEALLOCATE( obuf )

      
      IF (stored_write_record) THEN

        SELECT CASE ( use_package(io_form) )
          CASE ( IO_NETCDF   )
            CALL ext_ncd_iosync( handle(DataHandle), Status )
          CASE ( IO_GRIB1   )
            CALL ext_gr1_iosync( handle(DataHandle), Status )
          CASE ( IO_GRIB2   )
            CALL ext_gr2_iosync( handle(DataHandle), Status )
          CASE ( IO_INTIO   )
            CALL ext_int_iosync( handle(DataHandle), Status )
          CASE DEFAULT
            Status = 0
        END SELECT

      ENDIF

      END DO 

    END SUBROUTINE quilt

    SUBROUTINE quilt_pnc





      USE module_state_description
      USE module_quilt_outbuf_ops
      IMPLICIT NONE
      INCLUDE 'mpif.h'
  INTEGER, PARAMETER ::  int_ioexit			=  	     10
  INTEGER, PARAMETER ::  int_open_for_write_begin	=  	     20
  INTEGER, PARAMETER ::  int_open_for_write_commit	=  	     30
  INTEGER, PARAMETER ::  int_open_for_read 		=  	     40
  INTEGER, PARAMETER ::  int_inquire_opened 		=  	     60
  INTEGER, PARAMETER ::  int_inquire_filename 		=  	     70
  INTEGER, PARAMETER ::  int_iosync 			=  	     80
  INTEGER, PARAMETER ::  int_ioclose 			=  	     90
  INTEGER, PARAMETER ::  int_next_time 			=  	    100
  INTEGER, PARAMETER ::  int_set_time 			=  	    110
  INTEGER, PARAMETER ::  int_next_var 			=  	    120
  INTEGER, PARAMETER ::  int_dom_ti_real 		=  	    140
  INTEGER, PARAMETER ::  int_dom_ti_double 		=  	    160
  INTEGER, PARAMETER ::  int_dom_ti_integer 		=  	    180
  INTEGER, PARAMETER ::  int_dom_ti_logical 		=  	    200
  INTEGER, PARAMETER ::  int_dom_ti_char 		=  	    220
  INTEGER, PARAMETER ::  int_dom_td_real 		=  	    240
  INTEGER, PARAMETER ::  int_dom_td_double 		=  	    260
  INTEGER, PARAMETER ::  int_dom_td_integer 		=  	    280
  INTEGER, PARAMETER ::  int_dom_td_logical 		=  	    300
  INTEGER, PARAMETER ::  int_dom_td_char 		=  	    320
  INTEGER, PARAMETER ::  int_var_ti_real 		=  	    340
  INTEGER, PARAMETER ::  int_var_ti_double 		=  	    360
  INTEGER, PARAMETER ::  int_var_ti_integer 		=  	    380
  INTEGER, PARAMETER ::  int_var_ti_logical 		=  	    400
  INTEGER, PARAMETER ::  int_var_ti_char 		=  	    420
  INTEGER, PARAMETER ::  int_var_td_real 		=  	    440
  INTEGER, PARAMETER ::  int_var_td_double 		=  	    460
  INTEGER, PARAMETER ::  int_var_td_integer 		=  	    480
  INTEGER, PARAMETER ::  int_var_td_logical 		=  	    500
  INTEGER, PARAMETER ::  int_var_td_char 		=  	    520
  INTEGER, PARAMETER ::  int_field 			=  	    530
  INTEGER, PARAMETER ::  int_var_info 			=  	    540
  INTEGER, PARAMETER ::  int_noop 			=  	    550
      integer, parameter  :: WRF_FILE_NOT_OPENED                  = 100
      integer, parameter  :: WRF_FILE_OPENED_NOT_COMMITTED        = 101
      integer, parameter  :: WRF_FILE_OPENED_FOR_WRITE            = 102
      integer, parameter  :: WRF_FILE_OPENED_FOR_READ             = 103
      integer, parameter  :: WRF_REAL                             = 104
      integer, parameter  :: WRF_DOUBLE                           = 105
      integer, parameter  :: WRF_FLOAT=WRF_REAL
      integer, parameter  :: WRF_INTEGER                          = 106
      integer, parameter  :: WRF_LOGICAL                          = 107
      integer, parameter  :: WRF_COMPLEX                          = 108
      integer, parameter  :: WRF_DOUBLE_COMPLEX                   = 109
      integer, parameter  :: WRF_FILE_OPENED_FOR_UPDATE           = 110


      integer, parameter  :: WRF_FILE_OPENED_AND_COMMITTED        = 102
      INTEGER itag, ninbuf, ntasks_io_group, ntasks_local_group, mytask_local, ierr
      INTEGER istat
      INTEGER mytask_io_group
      INTEGER   :: nout_set = 0
      INTEGER   :: obufsize, bigbufsize, chunksize, sz
      REAL,                 DIMENSION(1) :: dummy
      INTEGER, ALLOCATABLE, DIMENSION(:) :: obuf, bigbuf
      REAL,    ALLOCATABLE, DIMENSION(:) :: RDATA
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IDATA
      CHARACTER (LEN=512) :: CDATA
      CHARACTER (LEN=80) :: fname
      INTEGER icurs, hdrbufsize, itypesize, ftypesize, rtypesize, Status, fstat, io_form_arg
      INTEGER :: DataHandle, FieldType, Comm, IOComm, DomainDesc, code, Count
      INTEGER, DIMENSION(3) :: DomainStart , DomainEnd , MemoryStart , MemoryEnd , PatchStart , PatchEnd
      INTEGER :: dummybuf(1)
      INTEGER :: num_noops, num_commit_messages, num_field_training_msgs, hdr_tag
      CHARACTER (len=256) :: DateStr , Element, VarName, MemoryOrder , Stagger , DimNames(3), FileName, SysDepInfo, mess
      INTEGER, EXTERNAL :: use_package
      LOGICAL           :: stored_write_record, retval, written_record
      INTEGER iii, jjj, vid





      SysDepInfo = " "
      CALL ext_ncd_ioinit( SysDepInfo, ierr)
      CALL ext_pnc_ioinit( SysDepInfo, ierr)
      CALL ext_int_ioinit( SysDepInfo, ierr )
      CALL ext_gr1_ioinit( SysDepInfo, ierr)
      CALL ext_gr2_ioinit( SysDepInfo, ierr)

      okay_to_commit = .false.
      stored_write_record = .false.
      ninbuf = 0
      
      
      CALL mpi_x_comm_size( mpi_comm_io_groups(1,1), ntasks_io_group,    ierr )
      CALL MPI_COMM_RANK( mpi_comm_io_groups(1,1), mytask_io_group,    ierr )
      CALL mpi_x_comm_size( mpi_comm_local,        ntasks_local_group, ierr )
      CALL MPI_COMM_RANK( mpi_comm_local,        mytask_local,       ierr )

      CALL MPI_TYPE_SIZE( MPI_INTEGER, itypesize, ierr )
      IF ( itypesize <= 0 ) THEN
        CALL wrf_error_fatal3("<stdin>",1489,&
"external/RSL/module_dm.F: quilt: type size <= 0 invalid")
      ENDIF











      DO WHILE (.TRUE.)  














        

        
        ! if needed (currently needed only for ioclose).
        reduced_dummy = 0
        CALL MPI_Reduce( reduced_dummy, reduced, 2, MPI_INTEGER, MPI_SUM, mytask_io_group, mpi_comm_io_groups(1,1), ierr )
        obufsize = reduced(1)



        IF ( obufsize .LT. 0 ) THEN
          IF ( obufsize .EQ. -100 ) THEN         
            CALL ext_ncd_ioexit( Status )
            CALL ext_pnc_ioexit( Status )
            CALL ext_int_ioexit( Status )
            CALL ext_gr1_ioexit( Status )
            CALL ext_gr2_ioexit( Status )
            CALL wrf_message ( 'I/O QUILT SERVERS DONE' )
            IF (coupler_on) THEN
               CALL cpl_finalize()
            ELSE
               CALL mpi_finalize(ierr)
            END IF
            STOP
          ELSE
            WRITE(mess,*)'Possible 32-bit overflow on output server. Try larger nio_tasks_per_group in namelist.'
            CALL wrf_error_fatal3("<stdin>",1544,&
mess)
          ENDIF
        ENDIF







        IF ( obufsize .GT. 0 ) THEN
          ALLOCATE( obuf( (obufsize+1)/itypesize ) )


          CALL collect_on_comm_debug("module_io_quilt.G",1748, mpi_comm_io_groups(1,1),        &
                                onebyte,                      &
                                dummy, 0,                     &
                                obuf, obufsize )

        ELSE
          ! Necessarily, the compute processes send the ioclose signal,
          
          ! will stall on the ioclose message waiting for the quilt 
          
          
          
          
          
          ! Then a header representing the ioclose message is constructed
          
          
          
          
          
          
          ALLOCATE( obuf( 4096 ) )
          
          CALL int_gen_handle_header( obuf, obufsize, itypesize, &
                                      reduced(2) , int_ioclose )
        ENDIF












        CALL init_store_piece_of_field
        CALL mpi_type_size ( MPI_INTEGER , itypesize , ierr )



        vid = 0
        icurs = itypesize
        num_noops = 0 
        num_commit_messages = 0 
        num_field_training_msgs = 0 
        DO WHILE ( icurs .lt. obufsize ) 
          hdr_tag = get_hdr_tag( obuf ( icurs / itypesize ) )
          SELECT CASE ( hdr_tag )
            CASE ( int_field )
              CALL int_get_write_field_header ( obuf(icurs/itypesize), hdrbufsize, itypesize, ftypesize,  &
                                                DataHandle , DateStr , VarName , Dummy , FieldType , Comm , IOComm, &
                                                DomainDesc , MemoryOrder , Stagger , DimNames ,              &
                                                DomainStart , DomainEnd ,                                    &
                                                MemoryStart , MemoryEnd ,                                    &
                                                PatchStart , PatchEnd )
              chunksize = (PatchEnd(1)-PatchStart(1)+1)*(PatchEnd(2)-PatchStart(2)+1)* &
                          (PatchEnd(3)-PatchStart(3)+1)*ftypesize

              IF ( DomainDesc .EQ. 333933 ) THEN  
                 IF ( num_field_training_msgs .EQ. 0 ) THEN
                   call add_to_bufsize_for_field( VarName, hdrbufsize )

                 ENDIF
                 num_field_training_msgs = num_field_training_msgs + 1
              ELSE
                 call add_to_bufsize_for_field( VarName, hdrbufsize )

              ENDIF
              icurs = icurs + hdrbufsize



              
              
              IF ( DomainDesc .NE. 333933 ) THEN   

                call add_to_bufsize_for_field( VarName, chunksize )
                icurs = icurs + chunksize
              ENDIF
            CASE ( int_open_for_write_commit )  
              hdrbufsize = obuf(icurs/itypesize)
              IF (num_commit_messages.EQ.0) THEN
                call add_to_bufsize_for_field( 'COMMIT', hdrbufsize )
              ENDIF
              num_commit_messages = num_commit_messages + 1
              icurs = icurs + hdrbufsize
            CASE DEFAULT
              hdrbufsize = obuf(icurs/itypesize)





























              IF (hdr_tag.NE.int_noop) THEN
                write(VarName,'(I5.5)')vid 

                call add_to_bufsize_for_field( VarName, hdrbufsize )
                vid = vid+1
              ENDIF
              IF ( hdr_tag .EQ. int_noop ) num_noops = num_noops + 1
              icurs = icurs + hdrbufsize

          END SELECT
        ENDDO 



        vid = 0
        icurs = itypesize
        num_noops = 0 
        num_commit_messages = 0 
        num_field_training_msgs = 0 
        DO WHILE ( icurs .lt. obufsize ) 

          hdr_tag = get_hdr_tag( obuf ( icurs / itypesize ) )
          SELECT CASE ( hdr_tag )
            CASE ( int_field )
              CALL int_get_write_field_header ( obuf(icurs/itypesize), hdrbufsize, itypesize, ftypesize,  &
                                                DataHandle , DateStr , VarName , Dummy , FieldType , Comm , IOComm, &
                                                DomainDesc , MemoryOrder , Stagger , DimNames ,              &
                                                DomainStart , DomainEnd ,                                    &
                                                MemoryStart , MemoryEnd ,                                    &
                                                PatchStart , PatchEnd )
              chunksize = (PatchEnd(1)-PatchStart(1)+1)*(PatchEnd(2)-PatchStart(2)+1)* &
                          (PatchEnd(3)-PatchStart(3)+1)*ftypesize

              IF ( DomainDesc .EQ. 333933 ) THEN  
                 IF ( num_field_training_msgs .EQ. 0 ) THEN
                   call store_piece_of_field( obuf(icurs/itypesize), VarName, hdrbufsize )

                 ENDIF
                 num_field_training_msgs = num_field_training_msgs + 1
              ELSE
                 call store_piece_of_field( obuf(icurs/itypesize), VarName, hdrbufsize )

              ENDIF
              icurs = icurs + hdrbufsize
              
              
              IF ( DomainDesc .NE. 333933 ) THEN   
                call store_piece_of_field( obuf(icurs/itypesize), VarName, chunksize )
                icurs = icurs + chunksize

              ENDIF
            CASE ( int_open_for_write_commit )  
              hdrbufsize = obuf(icurs/itypesize)
              IF (num_commit_messages.EQ.0) THEN
                call store_piece_of_field( obuf(icurs/itypesize), 'COMMIT', hdrbufsize )
              ENDIF
              num_commit_messages = num_commit_messages + 1
              icurs = icurs + hdrbufsize
            CASE DEFAULT
              hdrbufsize = obuf(icurs/itypesize)
              IF (hdr_tag.NE.int_noop) THEN

                write(VarName,'(I5.5)')vid 

                call store_piece_of_field( obuf(icurs/itypesize), VarName, hdrbufsize )
                vid = vid+1
              ENDIF
              IF ( hdr_tag .EQ. int_noop ) num_noops = num_noops + 1
              icurs = icurs + hdrbufsize
          END SELECT
       ENDDO 



       CALL init_retrieve_pieces_of_field


       CALL retrieve_pieces_of_field ( obuf , VarName, obufsize, sz, retval )
       written_record = .false.


       DO WHILE ( retval ) 




            icurs = itypesize  

            stored_write_record = .false.



            DO WHILE ( icurs .lt. sz)




              SELECT CASE ( get_hdr_tag( obuf(icurs/itypesize) ) )


                CASE ( int_noop )
                  CALL int_get_noop_header( obuf(icurs/itypesize), &
                                            hdrbufsize, itypesize )
                  icurs = icurs + hdrbufsize


                CASE ( int_dom_td_real )
                  CALL mpi_type_size( MPI_REAL, ftypesize, ierr )
                  ALLOCATE( RData( obuf(icurs/itypesize + 4 ) ) )      
                  CALL int_get_td_header( obuf(icurs/itypesize:), hdrbufsize, itypesize, ftypesize, &
                                          DataHandle, DateStr, Element, RData, Count, code )
                  icurs = icurs + hdrbufsize

                  SELECT CASE (use_package(io_form(DataHandle)))
                    CASE (IO_PNETCDF  )
                      CALL ext_pnc_put_dom_td_real( handle(DataHandle),TRIM(Element),TRIM(DateStr),RData, Count, Status )
                    CASE ( IO_NETCDF   )
                      CALL ext_ncd_put_dom_td_real( handle(DataHandle),TRIM(Element),TRIM(DateStr),RData, Count, Status )
                    CASE ( IO_INTIO   )
                      CALL ext_int_put_dom_td_real( handle(DataHandle),TRIM(Element),TRIM(DateStr),RData, Count, Status )
                 CASE ( IO_GRIB1 )
                    CALL ext_gr1_put_dom_td_real( handle(DataHandle),TRIM(Element),TRIM(DateStr),RData, Count, Status )
                 CASE ( IO_GRIB2 )
                    CALL ext_gr2_put_dom_td_real( handle(DataHandle),TRIM(Element),TRIM(DateStr),RData, Count, Status )
                     CASE DEFAULT
                      Status = 0
                  END SELECT

                  DEALLOCATE( RData )

                CASE ( int_dom_ti_real )

                  CALL mpi_type_size( MPI_REAL, ftypesize, ierr )
                  ALLOCATE( RData( obuf(icurs/itypesize + 4 ) ) )      
                  CALL int_get_ti_header( obuf(icurs/itypesize:), hdrbufsize, itypesize, ftypesize, &
                                          DataHandle, Element, RData, Count, code )
                  icurs = icurs + hdrbufsize

                  SELECT CASE (use_package(io_form(DataHandle)))
                    CASE (IO_PNETCDF  )
                      CALL ext_pnc_put_dom_ti_real( handle(DataHandle),TRIM(Element), RData, Count, Status )
                    CASE ( IO_NETCDF   )
                      CALL ext_ncd_put_dom_ti_real( handle(DataHandle),TRIM(Element), RData, Count, Status )
                    CASE ( IO_INTIO   )
                      CALL ext_int_put_dom_ti_real( handle(DataHandle),TRIM(Element), RData, Count, Status )
                 CASE ( IO_GRIB1 )
                    CALL ext_gr1_put_dom_ti_real( handle(DataHandle),TRIM(Element), RData, Count, Status )
                 CASE ( IO_GRIB2 )
                    CALL ext_gr2_put_dom_ti_real( handle(DataHandle),TRIM(Element), RData, Count, Status )
                    CASE DEFAULT
                      Status = 0
                  END SELECT

                  DEALLOCATE( RData )


                CASE ( int_dom_td_integer )

                  CALL mpi_type_size( MPI_INTEGER, ftypesize, ierr )
                  ALLOCATE( IData( obuf(icurs/itypesize + 4 ) ) )      
                  CALL int_get_td_header( obuf(icurs/itypesize:), hdrbufsize, itypesize, ftypesize, &
                                          DataHandle, DateStr, Element, IData, Count, code )
                  icurs = icurs + hdrbufsize

                  SELECT CASE (use_package(io_form(DataHandle)))
                  CASE (IO_PNETCDF  )
                      CALL ext_pnc_put_dom_td_integer( handle(DataHandle),TRIM(Element), Trim(DateStr), IData, Count, Status )
                   CASE ( IO_NETCDF   )
                      CALL ext_ncd_put_dom_td_integer( handle(DataHandle),TRIM(Element), Trim(DateStr), IData, Count, Status )
                   CASE ( IO_INTIO   )
                      CALL ext_int_put_dom_td_integer( handle(DataHandle),TRIM(Element), Trim(DateStr), IData, Count, Status )
                   CASE ( IO_GRIB1 )
                      CALL ext_gr1_put_dom_td_integer( handle(DataHandle),TRIM(Element), Trim(DateStr), IData, Count, Status )
                   CASE ( IO_GRIB2 )
                      CALL ext_gr2_put_dom_td_integer( handle(DataHandle),TRIM(Element), Trim(DateStr), IData, Count, Status )
                   CASE DEFAULT
                      Status = 0
                   END SELECT

                   DEALLOCATE( IData )


                CASE ( int_dom_ti_integer )

                  CALL mpi_type_size( MPI_INTEGER, ftypesize, ierr )
                  ALLOCATE( IData( obuf(icurs/itypesize + 4 ) ) )      
                  CALL int_get_ti_header( obuf(icurs/itypesize:), hdrbufsize, itypesize, ftypesize, &
                                          DataHandle, Element, IData, Count, code )
                  icurs = icurs + hdrbufsize
                  SELECT CASE (use_package(io_form(DataHandle)))
                    CASE (IO_PNETCDF  )
                      CALL ext_pnc_put_dom_ti_integer( handle(DataHandle),TRIM(Element), IData, Count, Status )
                    CASE ( IO_NETCDF   )
                      CALL ext_ncd_put_dom_ti_integer( handle(DataHandle),TRIM(Element), IData, Count, Status )
                    CASE ( IO_INTIO   )
                      CALL ext_int_put_dom_ti_integer( handle(DataHandle),TRIM(Element), IData, Count, Status )
                 CASE ( IO_GRIB1 )
                    CALL ext_gr1_put_dom_ti_integer( handle(DataHandle),TRIM(Element), IData, Count, Status )
                 CASE ( IO_GRIB2 )
                    CALL ext_gr2_put_dom_ti_integer( handle(DataHandle),TRIM(Element), IData, Count, Status )

                    CASE DEFAULT
                      Status = 0
                  END SELECT

                  DEALLOCATE( IData)
 

                CASE ( int_set_time )

                  CALL int_get_ti_header_char( obuf(icurs/itypesize), hdrbufsize, itypesize, &
                                               DataHandle, Element, VarName, CData, code )
                  SELECT CASE (use_package(io_form(DataHandle)))
                    CASE ( IO_INTIO   )
                      CALL ext_int_set_time ( handle(DataHandle), TRIM(CData), Status)
                    CASE DEFAULT
                      Status = 0
                  END SELECT

                  icurs = icurs + hdrbufsize


                CASE ( int_dom_ti_char )

                  CALL int_get_ti_header_char( obuf(icurs/itypesize), hdrbufsize, itypesize, &
                                               DataHandle, Element, VarName, CData, code )

                  SELECT CASE (use_package(io_form(DataHandle)))
                    CASE (IO_PNETCDF  )
                      CALL ext_pnc_put_dom_ti_char ( handle(DataHandle), TRIM(Element), Trim(CData), Status)
                    CASE ( IO_NETCDF   )
                      CALL ext_ncd_put_dom_ti_char ( handle(DataHandle), TRIM(Element), TRIM(CData), Status)
                    CASE ( IO_INTIO   )
                      CALL ext_int_put_dom_ti_char ( handle(DataHandle), TRIM(Element), TRIM(CData), Status)
                   CASE ( IO_GRIB1 )
                      CALL ext_gr1_put_dom_ti_char ( handle(DataHandle), TRIM(Element), TRIM(CData), Status)
                   CASE ( IO_GRIB2 )
                      CALL ext_gr2_put_dom_ti_char ( handle(DataHandle), TRIM(Element), TRIM(CData), Status)
                   CASE DEFAULT
                      Status = 0
                   END SELECT

                  icurs = icurs + hdrbufsize


                CASE ( int_var_ti_char )

                  CALL int_get_ti_header_char( obuf(icurs/itypesize), hdrbufsize, itypesize, &
                                               DataHandle, Element, VarName, CData, code )

                  SELECT CASE (use_package(io_form(DataHandle)))
                    CASE (IO_PNETCDF  )
                      CALL ext_pnc_put_var_ti_char ( handle(DataHandle), TRIM(Element), TRIM(VarName), TRIM(CData), Status )
                    CASE ( IO_NETCDF   )
                      CALL ext_ncd_put_var_ti_char ( handle(DataHandle), TRIM(Element), TRIM(VarName), TRIM(CData), Status)
                    CASE ( IO_INTIO   )
                      CALL ext_int_put_var_ti_char ( handle(DataHandle), TRIM(Element), TRIM(VarName), TRIM(CData), Status)
                   CASE ( IO_GRIB1 )
                      CALL ext_gr1_put_var_ti_char ( handle(DataHandle), TRIM(Element), TRIM(VarName), TRIM(CData), Status)
                   CASE ( IO_GRIB2 )
                      CALL ext_gr2_put_var_ti_char ( handle(DataHandle), TRIM(Element), TRIM(VarName), TRIM(CData), Status)
                   CASE DEFAULT
                      Status = 0
                   END SELECT

                  icurs = icurs + hdrbufsize

                CASE ( int_ioexit )

                  CALL wrf_error_fatal3("<stdin>",1948,&
                         "quilt: should have handled int_ioexit already")
! Every I/O server handles the "ioclose" request.
                CASE ( int_ioclose )
                  CALL int_get_handle_header( obuf(icurs/itypesize), hdrbufsize, itypesize, &
                                              DataHandle , code )
                  icurs = icurs + hdrbufsize

                  IF ( DataHandle .GE. 1 ) THEN

                     SELECT CASE (use_package(io_form(DataHandle)))
                    CASE ( IO_PNETCDF   )
                      CALL ext_pnc_inquire_filename( handle(DataHandle), fname, fstat, Status )
                      IF ( fstat .EQ. WRF_FILE_OPENED_FOR_WRITE .OR. fstat .EQ. WRF_FILE_OPENED_NOT_COMMITTED ) THEN
                        CALL ext_pnc_ioclose(handle(DataHandle),Status)
                      ENDIF
                     CASE ( IO_NETCDF   )
                        CALL ext_ncd_inquire_filename( handle(DataHandle), fname, fstat, Status )
                        IF ( fstat .EQ. WRF_FILE_OPENED_FOR_WRITE .OR. fstat .EQ. WRF_FILE_OPENED_NOT_COMMITTED ) THEN
                           CALL ext_ncd_ioclose(handle(DataHandle),Status)
                        ENDIF
                     CASE ( IO_INTIO   )
                        CALL ext_int_inquire_filename( handle(DataHandle), fname, fstat, Status )
                        IF ( fstat .EQ. WRF_FILE_OPENED_FOR_WRITE .OR. fstat .EQ. WRF_FILE_OPENED_NOT_COMMITTED ) THEN
                           CALL ext_int_ioclose(handle(DataHandle),Status)
                        ENDIF
                     CASE ( IO_GRIB1 )
                        CALL ext_gr1_inquire_filename( handle(DataHandle), fname, fstat, Status )
                        IF ( fstat .EQ. WRF_FILE_OPENED_FOR_WRITE .OR. fstat .EQ. WRF_FILE_OPENED_NOT_COMMITTED ) THEN
                           CALL ext_gr1_ioclose(handle(DataHandle),Status)
                        ENDIF
                     CASE ( IO_GRIB2 )
                        CALL ext_gr2_inquire_filename( handle(DataHandle), fname, fstat, Status )
                        IF ( fstat .EQ. WRF_FILE_OPENED_FOR_WRITE .OR. fstat .EQ. WRF_FILE_OPENED_NOT_COMMITTED ) THEN
                           CALL ext_gr2_ioclose(handle(DataHandle),Status)
                        ENDIF
                     CASE DEFAULT
                        Status = 0
                     END SELECT
                  ENDIF


                CASE ( int_open_for_write_begin )

                  CALL int_get_ofwb_header( obuf(icurs/itypesize), hdrbufsize, itypesize, &
                                            FileName,SysDepInfo,io_form_arg,DataHandle )





                  icurs = icurs + hdrbufsize

                
                  io_form(DataHandle) = io_form_arg

                  SELECT CASE (use_package(io_form(DataHandle)))
                    CASE (IO_PNETCDF  )
                      CALL ext_pnc_open_for_write_begin(FileName,mpi_comm_local,mpi_comm_local,SysDepInfo,handle(DataHandle),Status )
                    CASE ( IO_NETCDF   )
                      CALL ext_ncd_open_for_write_begin(FileName,Comm,IOComm,SysDepInfo,handle(DataHandle),Status)

                    CASE ( IO_INTIO   )
                      CALL ext_int_open_for_write_begin(FileName,Comm,IOComm,SysDepInfo,handle(DataHandle),Status)
                    CASE ( IO_GRIB1 )
                       CALL ext_gr1_open_for_write_begin(FileName,Comm,IOComm,SysDepInfo,handle(DataHandle),Status)
                    CASE ( IO_GRIB2 )
                       CALL ext_gr2_open_for_write_begin(FileName,Comm,IOComm,SysDepInfo,handle(DataHandle),Status)
                    CASE DEFAULT
                      Status = 0
                  END SELECT
                
                  okay_to_write(DataHandle) = .false.





                CASE ( int_open_for_write_commit )

                  CALL int_get_handle_header( obuf(icurs/itypesize), hdrbufsize, itypesize, &
                                              DataHandle , code )
                  icurs = icurs + hdrbufsize
                  okay_to_commit(DataHandle) = .true.










                CASE ( int_field )
                  CALL mpi_type_size( MPI_INTEGER, ftypesize, ierr )
                  CALL int_get_write_field_header ( obuf(icurs/itypesize), hdrbufsize, itypesize, ftypesize,  &
                                                    DataHandle , DateStr , VarName , Dummy , FieldType , Comm , IOComm, &
                                                    DomainDesc , MemoryOrder , Stagger , DimNames ,              &
                                                    DomainStart , DomainEnd ,                                    &
                                                    MemoryStart , MemoryEnd ,                                    &
                                                    PatchStart , PatchEnd )

                  icurs = icurs + hdrbufsize

                  IF ( okay_to_write(DataHandle) ) THEN









                    IF ( FieldType .EQ. WRF_FLOAT .OR. FieldType .EQ. WRF_DOUBLE)  THEN
                      
                      
                      IF ( FieldType .EQ. WRF_DOUBLE)  THEN

                        CALL mpi_type_size( MPI_DOUBLE_PRECISION, ftypesize, ierr )
                      ELSE
                        CALL mpi_type_size( MPI_REAL, ftypesize, ierr )
                      ENDIF







                      CALL store_patch_in_outbuf_pnc(obuf(icurs/itypesize), &
                                                     dummybuf, TRIM(DateStr), &
                                                     TRIM(VarName) , &
                                                     FieldType,      &
                                                     TRIM(MemoryOrder), &
                                                     TRIM(Stagger), &
                                                     DimNames, &
                                                     DomainStart , DomainEnd ,&
                                                     MemoryStart , MemoryEnd ,&
                                                     PatchStart , PatchEnd, &
                                                     ntasks_io_group-1 )
                      stored_write_record = .true.







                    ELSE IF ( FieldType .EQ. WRF_INTEGER ) THEN
                      CALL mpi_type_size( MPI_INTEGER, ftypesize, ierr )
                      CALL store_patch_in_outbuf_pnc ( dummybuf,             & 
                                                   obuf(icurs/itypesize) ,   &
                                                   TRIM(DateStr) ,           &
                                                   TRIM(VarName) ,           &
                                                   FieldType,                &
                                                   TRIM(MemoryOrder) ,       &
                                                   TRIM(Stagger), DimNames,  &
                                                   DomainStart , DomainEnd , &
                                                   MemoryStart , MemoryEnd , &
                                                   PatchStart , PatchEnd   , &
                                                   ntasks_io_group-1 )
                      stored_write_record = .true.
                    ELSE IF ( FieldType .EQ. WRF_LOGICAL ) THEN
                      ftypesize = 4
                    ENDIF

                    icurs = icurs + (PatchEnd(1)-PatchStart(1)+1)* &
                                    (PatchEnd(2)-PatchStart(2)+1)* &
                                    (PatchEnd(3)-PatchStart(3)+1)*ftypesize

                  ELSE 

                    SELECT CASE (use_package(io_form(DataHandle)))

                      CASE ( IO_PNETCDF )
                        CALL ext_pnc_write_field ( handle(DataHandle) , TRIM(DateStr),        &
                                   TRIM(VarName) , dummy , FieldType , mpi_comm_local , mpi_comm_local,         &
                                   DomainDesc , TRIM(MemoryOrder) , TRIM(Stagger), DimNames , &
                                   DomainStart , DomainEnd ,                                  &
                                   MemoryStart , MemoryEnd ,                                  &
                                   PatchStart ,  PatchEnd,                                  &
                                   Status )
                      CASE ( IO_NETCDF   )
                        CALL ext_ncd_write_field ( handle(DataHandle) , TRIM(DateStr) ,         &
                                   TRIM(VarName) , dummy , FieldType , Comm , IOComm,           &
                                   DomainDesc , TRIM(MemoryOrder) , TRIM(Stagger) , DimNames ,  &
                                   DomainStart , DomainEnd ,                                    &
                                   DomainStart , DomainEnd ,                                    &
                                   DomainStart , DomainEnd ,                                    &
                                   Status )
                      CASE DEFAULT
                        Status = 0
                    END SELECT
                  ENDIF
                CASE ( int_iosync )
                  CALL int_get_handle_header( obuf(icurs/itypesize), hdrbufsize, itypesize, &
                                            DataHandle , code )
                  icurs = icurs + hdrbufsize
                CASE DEFAULT
                  WRITE(mess,*)'quilt: bad tag: ',                            &
                               get_hdr_tag( obuf(icurs/itypesize) ),' icurs ',&
                               icurs/itypesize
                  CALL wrf_error_fatal3("<stdin>",2152,&
mess )
              END SELECT

            ENDDO 



            IF (stored_write_record) THEN







              CALL write_outbuf_pnc( handle(DataHandle), &
                                     use_package(io_form(DataHandle)), &
                                     mpi_comm_local, mytask_local,     &
                                     ntasks_local_group) 

              stored_write_record = .false.
              written_record = .true.
            ENDIF




            IF (okay_to_commit(DataHandle)) THEN

              SELECT CASE (use_package(io_form(DataHandle)))
                CASE ( IO_PNETCDF   )
                  CALL ext_pnc_inquire_filename( handle(DataHandle), fname, fstat, Status )
                  IF ( fstat .EQ. WRF_FILE_OPENED_NOT_COMMITTED ) THEN
                    CALL ext_pnc_open_for_write_commit(handle(DataHandle),Status)
                    okay_to_write(DataHandle) = .true.
                  ENDIF
                CASE ( IO_NETCDF   )
                  CALL ext_ncd_inquire_filename( handle(DataHandle), fname, fstat, Status )
                  IF ( fstat .EQ. WRF_FILE_OPENED_NOT_COMMITTED ) THEN
                    CALL ext_ncd_open_for_write_commit(handle(DataHandle),Status)
                    okay_to_write(DataHandle) = .true.
                  ENDIF
                CASE ( IO_INTIO   )
                  CALL ext_int_inquire_filename( handle(DataHandle), fname, fstat, Status )
                  IF ( fstat .EQ. WRF_FILE_OPENED_NOT_COMMITTED ) THEN
                    CALL ext_int_open_for_write_commit(handle(DataHandle),Status)
                    okay_to_write(DataHandle) = .true.
                  ENDIF
                 CASE ( IO_GRIB1 )
                    CALL ext_gr1_inquire_filename( handle(DataHandle), fname, fstat, Status )
                    IF ( fstat .EQ. WRF_FILE_OPENED_NOT_COMMITTED ) THEN
                       CALL ext_gr1_open_for_write_commit(handle(DataHandle),Status)
                       okay_to_write(DataHandle) = .true.
                    ENDIF
                 CASE ( IO_GRIB2 )
                    CALL ext_gr2_inquire_filename( handle(DataHandle), fname, fstat, Status )
                    IF ( fstat .EQ. WRF_FILE_OPENED_NOT_COMMITTED ) THEN
                       CALL ext_gr2_open_for_write_commit(handle(DataHandle),Status)
                       okay_to_write(DataHandle) = .true.
                    ENDIF

                CASE DEFAULT
                  Status = 0
              END SELECT

            okay_to_commit(DataHandle) = .false.
          ENDIF




        CALL retrieve_pieces_of_field ( obuf , VarName, obufsize, sz, retval )
      END DO 

      DEALLOCATE( obuf )

      
      IF (written_record) THEN

        SELECT CASE ( use_package(io_form) )
          CASE ( IO_PNETCDF   )
            CALL ext_pnc_iosync( handle(DataHandle), Status )
          CASE DEFAULT
            Status = 0
        END SELECT
        written_record = .false.

      ENDIF

      END DO 

    END SUBROUTINE quilt_pnc



    SUBROUTINE init_module_wrf_quilt
      USE module_wrf_error, only: init_module_wrf_error
      USE module_dm, ONLY : compute_mesh,nest_pes_x, nest_pes_y, domain_active_this_task








      USE module_dm, ONLY : compute_mesh,nest_pes_x,nest_pes_y,domain_active_this_task,&
                            tasks_per_split,comm_start,dm_task_split
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER i,j
      NAMELIST /namelist_quilt/ nio_tasks_per_group, nio_groups, poll_servers
      INTEGER ntasks, mytask, ierr, io_status, id, itask
      INTEGER num_io_tasks, ncompute_tasks
      INTEGER n_x, n_y
      INTEGER mpi_comm_here, temp_poll, comdup
      INTEGER, ALLOCATABLE ::  icolor(:)
      LOGICAL mpi_inited
      LOGICAL compute_node
      LOGICAL esmf_coupling
      CHARACTER*256 message









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



      esmf_coupling = .FALSE.
      quilting_enabled = .FALSE.
      IF ( disable_quilt ) RETURN

      DO i = 1,int_num_handles
        okay_to_write(i) = .FALSE.
        int_handle_in_use(i) = .FALSE.
        which_grid_is_handle(i) = -1
        prev_server_for_handle(i) = -1
        int_num_bytes_to_write(i) = 0
      ENDDO
      DO j = 1, max_domains
        DO i = 1,int_num_handles
          server_for_handle(i,j) = 0 
        ENDDO
      ENDDO

      CALL MPI_INITIALIZED( mpi_inited, ierr )
      IF ( .NOT. mpi_inited ) THEN
        CALL mpi_init ( ierr )
        CALL wrf_set_dm_communicator( MPI_COMM_WORLD )
        CALL wrf_termio_dup(MPI_COMM_WORLD)
      ENDIF
      CALL wrf_get_dm_communicator( mpi_comm_here )   

      CALL MPI_Comm_rank ( mpi_comm_here, mytask, ierr ) ;
      CALL MPI_Comm_size ( mpi_comm_here, ntasks, ierr ) ;

      IF ( mytask .EQ. 0 ) THEN
        OPEN ( unit=27, file="namelist.input", form="formatted", status="old" )
        tasks_per_split = ntasks


        READ ( UNIT = 27 , NML = domains , IOSTAT=io_status )
        REWIND(27)
        nio_groups = 1
        nio_tasks_per_group  = 0
        poll_servers = .false.
        READ ( 27 , NML = namelist_quilt, IOSTAT=io_status )
        IF (io_status .NE. 0) THEN
          CALL wrf_error_fatal3("<stdin>",5575,&
"ERROR reading namelist namelist_quilt" )
        ENDIF
        CLOSE ( 27 )
        IF ( esmf_coupling ) THEN
          IF ( any ( nio_tasks_per_group > 0 ) ) THEN
            CALL wrf_error_fatal3("<stdin>",5581,&
"frame/module_io_quilt.F: cannot use "// &
                                 "ESMF coupling with quilt tasks") ;
          ENDIF
        ENDIF
        if(poll_servers) then
           temp_poll=1
        else
           temp_poll=0
        endif
      ENDIF
      CALL mpi_bcast( parent_id, max_domains , MPI_INTEGER , 0 , mpi_comm_here, ierr )
      CALL mpi_bcast( max_dom, 1 , MPI_INTEGER , 0 , mpi_comm_here, ierr )
      CALL nl_set_max_dom(1,max_dom)
      CALL mpi_bcast( nio_tasks_per_group  , max_domains , MPI_INTEGER , 0 , mpi_comm_here, ierr )
      CALL mpi_bcast( nio_groups , 1 , MPI_INTEGER , 0 , mpi_comm_here, ierr )
      CALL mpi_bcast( temp_poll , 1 , MPI_INTEGER , 0 , mpi_comm_here, ierr )

      i = MAXVAL(nio_tasks_per_group(1:max_dom))
      IF ( i .GT. 0 .AND. nio_groups .GT. 0 ) THEN
        DO id = 1, max_dom
          IF ( nio_tasks_per_group(id) .LE. 0 ) THEN
             CALL wrf_error_fatal3("<stdin>",5603,&
'If nio_tasks_per_group in namelist.input is non-zero for any domain, every active domain must have a non-zero value in nio_tasks_per_group')
          ENDIF
        ENDDO
      ENDIF

      num_io_tasks = 0
      DO id = 1, max_dom
        num_io_tasks = num_io_tasks + nio_tasks_per_group(id)*nio_groups
      ENDDO
      IF ( mytask .EQ. 0 ) THEN
        OPEN ( unit=27, file="namelist.input", form="formatted", status="old" )
        comm_start = -1   
        nest_pes_x = 0    
        nest_pes_y = 0
        CALL compute_mesh( ntasks-num_io_tasks, n_x, n_y )
        comm_start = 0   
        nest_pes_x(1:max_dom) = n_x
        nest_pes_y(1:max_dom) = n_y
        READ ( 27 , NML = dm_task_split, IOSTAT=io_status )



        CLOSE ( 27 )
      ENDIF

      CALL mpi_bcast( io_status, 1 , MPI_INTEGER , 0 , mpi_comm_here, ierr )

      IF ( io_status .NE. 0 ) THEN
      ENDIF
      CALL mpi_bcast( tasks_per_split, 1 , MPI_INTEGER , 0 , mpi_comm_here, ierr )
      CALL mpi_bcast( comm_start, max_domains , MPI_INTEGER , 0 , mpi_comm_here, ierr )
      CALL mpi_bcast( nest_pes_x, max_domains , MPI_INTEGER , 0 , mpi_comm_here, ierr )
      CALL mpi_bcast( nest_pes_y, max_domains , MPI_INTEGER , 0 , mpi_comm_here, ierr )

      CALL MPI_Comm_rank ( mpi_comm_here, mytask, ierr ) ;
      CALL mpi_x_comm_size( mpi_comm_here, ntasks, ierr ) ;






      ALLOCATE(role_for_task(ntasks))
      
      ncompute_tasks = 0
      DO id = 1, max_dom
         IF ( ncompute_tasks < comm_start(id)+nest_pes_x(id)*nest_pes_y(id) ) THEN
              ncompute_tasks = comm_start(id)+nest_pes_x(id)*nest_pes_y(id)
         ENDIF
      ENDDO
      num_io_tasks = 0
      DO id = 1, max_dom
        num_io_tasks = num_io_tasks + nio_tasks_per_group(id)*nio_groups
      ENDDO


      IF ( ncompute_tasks + num_io_tasks .GT. ntasks ) THEN
        WRITE(message,"('ncompute_tasks(',i9,')+num_io_tasks(',i9,') .GT. ntasks(',i9,')')")ncompute_tasks,num_io_tasks,ntasks
        CALL wrf_error_fatal3("<stdin>",5662,&
TRIM(message))
      ENDIF
      DO itask = 1, ncompute_tasks
        role_for_task(itask) = 0
      ENDDO
      itask = ncompute_tasks + 1
      DO id = 1, max_dom
        DO i = 1, nio_tasks_per_group(id)*nio_groups
          role_for_task(itask) = id  
          itask = itask + 1
        ENDDO
      ENDDO
      ntasks = itask - 1


      poll_servers = (temp_poll == 1)

      compute_group_master = .FALSE.
      compute_node         = .FALSE.

DO id = 1, max_dom







      CALL setup_quilt_servers( id, nio_tasks_per_group,  &
                                role_for_task,            &  
                                num_io_tasks,             &
                                ncompute_tasks,           &
                                mytask,                   &
                                ntasks,                   &
                                nio_groups,               &
                                mpi_comm_here,            &
                                mpi_comm_local,           &  
                                mpi_comm_io_groups,       &
                                compute_node      )

ENDDO

      call init_module_wrf_error(on_io_server=.true.)


      CALL MPI_Comm_dup( mpi_comm_here, comdup, ierr )
      
      CALL MPI_Comm_split(comdup,role_for_task(mytask+1),mytask,mpi_comm_local, ierr )


      CALL wrf_set_dm_communicator( mpi_comm_local ) 
      IF ( compute_node ) THEN
          IF (coupler_on) CALL cpl_set_dm_communicator( mpi_comm_local )
          call ATM_SET_COMM(mpi_comm_local)
       ELSE
          call ATM_LEAVE_COUPLING()
          IF (coupler_on) CALL cpl_set_dm_communicator( MPI_COMM_NULL )
          mpi_comm_local = mpi_comm_local_io_server_tmp  
          CALL quilt    
       ENDIF

      RETURN
    END SUBROUTINE init_module_wrf_quilt




END MODULE module_wrf_quilt







SUBROUTINE disable_quilting




  USE module_wrf_quilt
  disable_quilt = .TRUE.
  RETURN
END SUBROUTINE disable_quilting

LOGICAL FUNCTION  use_output_servers_for(ioform)







  USE module_wrf_quilt
  integer, intent(in) :: ioform
  use_output_servers_for = quilting_enabled
  use_output_servers_for = ( use_output_servers_for .and. ioform<100 )
  RETURN
END FUNCTION use_output_servers_for

LOGICAL FUNCTION  use_output_servers()




  USE module_wrf_quilt
  use_output_servers = quilting_enabled
  RETURN
END FUNCTION use_output_servers

LOGICAL FUNCTION  use_input_servers()




  USE module_wrf_quilt
  use_input_servers = .FALSE.
  RETURN
END FUNCTION use_input_servers

SUBROUTINE wrf_quilt_open_for_write_begin( FileName , gridid, Comm_compute, Comm_io, SysDepInfo, &
                                     DataHandle , io_form_arg, Status )





  USE module_wrf_quilt
  USE module_state_description, ONLY: IO_PNETCDF
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  INTEGER, PARAMETER ::  int_ioexit			=  	     10
  INTEGER, PARAMETER ::  int_open_for_write_begin	=  	     20
  INTEGER, PARAMETER ::  int_open_for_write_commit	=  	     30
  INTEGER, PARAMETER ::  int_open_for_read 		=  	     40
  INTEGER, PARAMETER ::  int_inquire_opened 		=  	     60
  INTEGER, PARAMETER ::  int_inquire_filename 		=  	     70
  INTEGER, PARAMETER ::  int_iosync 			=  	     80
  INTEGER, PARAMETER ::  int_ioclose 			=  	     90
  INTEGER, PARAMETER ::  int_next_time 			=  	    100
  INTEGER, PARAMETER ::  int_set_time 			=  	    110
  INTEGER, PARAMETER ::  int_next_var 			=  	    120
  INTEGER, PARAMETER ::  int_dom_ti_real 		=  	    140
  INTEGER, PARAMETER ::  int_dom_ti_double 		=  	    160
  INTEGER, PARAMETER ::  int_dom_ti_integer 		=  	    180
  INTEGER, PARAMETER ::  int_dom_ti_logical 		=  	    200
  INTEGER, PARAMETER ::  int_dom_ti_char 		=  	    220
  INTEGER, PARAMETER ::  int_dom_td_real 		=  	    240
  INTEGER, PARAMETER ::  int_dom_td_double 		=  	    260
  INTEGER, PARAMETER ::  int_dom_td_integer 		=  	    280
  INTEGER, PARAMETER ::  int_dom_td_logical 		=  	    300
  INTEGER, PARAMETER ::  int_dom_td_char 		=  	    320
  INTEGER, PARAMETER ::  int_var_ti_real 		=  	    340
  INTEGER, PARAMETER ::  int_var_ti_double 		=  	    360
  INTEGER, PARAMETER ::  int_var_ti_integer 		=  	    380
  INTEGER, PARAMETER ::  int_var_ti_logical 		=  	    400
  INTEGER, PARAMETER ::  int_var_ti_char 		=  	    420
  INTEGER, PARAMETER ::  int_var_td_real 		=  	    440
  INTEGER, PARAMETER ::  int_var_td_double 		=  	    460
  INTEGER, PARAMETER ::  int_var_td_integer 		=  	    480
  INTEGER, PARAMETER ::  int_var_td_logical 		=  	    500
  INTEGER, PARAMETER ::  int_var_td_char 		=  	    520
  INTEGER, PARAMETER ::  int_field 			=  	    530
  INTEGER, PARAMETER ::  int_var_info 			=  	    540
  INTEGER, PARAMETER ::  int_noop 			=  	    550
  CHARACTER *(*), INTENT(IN)  :: FileName
  INTEGER ,       INTENT(IN)  :: gridid
  INTEGER ,       INTENT(IN)  :: Comm_compute , Comm_io
  CHARACTER *(*), INTENT(IN)  :: SysDepInfo
  INTEGER ,       INTENT(OUT) :: DataHandle
  INTEGER ,       INTENT(IN)  :: io_form_arg
  INTEGER ,       INTENT(OUT) :: Status

  CHARACTER*132   :: locFileName, locSysDepInfo
  INTEGER i, itypesize, tasks_in_group, ierr, comm_io_group
  REAL dummy
  INTEGER, EXTERNAL :: use_package

  CALL wrf_debug ( 50, 'in wrf_quilt_open_for_write_begin' ) 
  CALL int_get_fresh_handle(i)
  okay_to_write(i) = .false.
  DataHandle = i

  locFileName = FileName
  locSysDepInfo = SysDepInfo

  CALL MPI_TYPE_SIZE( MPI_INTEGER, itypesize, ierr )

  SELECT CASE(use_package(io_form_arg))

  CASE(IO_PNETCDF)
     IF(compute_group_master(1,current_id)) THEN
        CALL int_gen_ofwb_header( hdrbuf, hdrbufsize, itypesize, &
                                  locFileName,locSysDepInfo,io_form_arg,&
                                  DataHandle )
     ELSE
        CALL int_gen_noop_header( hdrbuf, hdrbufsize, itypesize )
     END IF
  CASE DEFAULT

     IF ( wrf_dm_on_monitor() ) THEN
        CALL int_gen_ofwb_header( hdrbuf, hdrbufsize, itypesize, &
                                  locFileName,locSysDepInfo,io_form_arg,DataHandle )
     ELSE
        CALL int_gen_noop_header( hdrbuf, hdrbufsize, itypesize )
     ENDIF

  END SELECT

  which_grid_is_handle(DataHandle) = gridid
  iserver = get_server_id ( DataHandle )




  CALL get_mpi_comm_io_groups( comm_io_group , iserver )


  CALL mpi_x_comm_size( comm_io_group , tasks_in_group , ierr )



  
  reduced = 0
  reduced(1) = hdrbufsize 
  IF ( compute_group_master(1,current_id) ) reduced(2) = i

  CALL MPI_Reduce( reduced, reduced_dummy, 2, MPI_INTEGER, MPI_SUM, tasks_in_group-1, comm_io_group, ierr )


  
  CALL collect_on_comm_debug("module_io_quilt.G",3178, comm_io_group,            &
                        onebyte,                       &
                        hdrbuf, hdrbufsize , &
                        dummy, 0 )

  Status = 0


  RETURN  
END SUBROUTINE wrf_quilt_open_for_write_begin

SUBROUTINE wrf_quilt_open_for_write_commit( DataHandle , Status )







  USE module_wrf_quilt
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  INTEGER, PARAMETER ::  int_ioexit			=  	     10
  INTEGER, PARAMETER ::  int_open_for_write_begin	=  	     20
  INTEGER, PARAMETER ::  int_open_for_write_commit	=  	     30
  INTEGER, PARAMETER ::  int_open_for_read 		=  	     40
  INTEGER, PARAMETER ::  int_inquire_opened 		=  	     60
  INTEGER, PARAMETER ::  int_inquire_filename 		=  	     70
  INTEGER, PARAMETER ::  int_iosync 			=  	     80
  INTEGER, PARAMETER ::  int_ioclose 			=  	     90
  INTEGER, PARAMETER ::  int_next_time 			=  	    100
  INTEGER, PARAMETER ::  int_set_time 			=  	    110
  INTEGER, PARAMETER ::  int_next_var 			=  	    120
  INTEGER, PARAMETER ::  int_dom_ti_real 		=  	    140
  INTEGER, PARAMETER ::  int_dom_ti_double 		=  	    160
  INTEGER, PARAMETER ::  int_dom_ti_integer 		=  	    180
  INTEGER, PARAMETER ::  int_dom_ti_logical 		=  	    200
  INTEGER, PARAMETER ::  int_dom_ti_char 		=  	    220
  INTEGER, PARAMETER ::  int_dom_td_real 		=  	    240
  INTEGER, PARAMETER ::  int_dom_td_double 		=  	    260
  INTEGER, PARAMETER ::  int_dom_td_integer 		=  	    280
  INTEGER, PARAMETER ::  int_dom_td_logical 		=  	    300
  INTEGER, PARAMETER ::  int_dom_td_char 		=  	    320
  INTEGER, PARAMETER ::  int_var_ti_real 		=  	    340
  INTEGER, PARAMETER ::  int_var_ti_double 		=  	    360
  INTEGER, PARAMETER ::  int_var_ti_integer 		=  	    380
  INTEGER, PARAMETER ::  int_var_ti_logical 		=  	    400
  INTEGER, PARAMETER ::  int_var_ti_char 		=  	    420
  INTEGER, PARAMETER ::  int_var_td_real 		=  	    440
  INTEGER, PARAMETER ::  int_var_td_double 		=  	    460
  INTEGER, PARAMETER ::  int_var_td_integer 		=  	    480
  INTEGER, PARAMETER ::  int_var_td_logical 		=  	    500
  INTEGER, PARAMETER ::  int_var_td_char 		=  	    520
  INTEGER, PARAMETER ::  int_field 			=  	    530
  INTEGER, PARAMETER ::  int_var_info 			=  	    540
  INTEGER, PARAMETER ::  int_noop 			=  	    550
  INTEGER ,       INTENT(IN ) :: DataHandle
  INTEGER ,       INTENT(OUT) :: Status
  INTEGER i, itypesize, tasks_in_group, ierr, comm_io_group
  REAL dummy

  CALL wrf_debug ( 50, 'in wrf_quilt_open_for_write_commit' ) 
  IF ( DataHandle .GE. 1 .AND. DataHandle .LE. int_num_handles ) THEN
    IF ( int_handle_in_use( DataHandle ) ) THEN
      okay_to_write( DataHandle ) = .true.
    ENDIF
  ENDIF

  CALL MPI_TYPE_SIZE( MPI_INTEGER, itypesize, ierr )



  IF(compute_group_master(1,current_id)) THEN
     CALL int_gen_handle_header( hdrbuf, hdrbufsize, itypesize, &
                                 DataHandle, int_open_for_write_commit )
  ELSE
     CALL int_gen_noop_header( hdrbuf, hdrbufsize, itypesize )
  END IF

  iserver = get_server_id ( DataHandle )
  CALL get_mpi_comm_io_groups( comm_io_group , iserver )

  CALL mpi_x_comm_size( comm_io_group , tasks_in_group , ierr )


  
  reduced = 0
  reduced(1) = hdrbufsize 
  IF ( compute_group_master(1,current_id) ) reduced(2) = DataHandle
  CALL MPI_Reduce( reduced, reduced_dummy, 2, MPI_INTEGER, MPI_SUM, tasks_in_group-1, comm_io_group, ierr )


  
  CALL collect_on_comm_debug("module_io_quilt.G",3254, comm_io_group,            &
                        onebyte,                       &
                        hdrbuf, hdrbufsize , &
                        dummy, 0 )

  Status = 0

  RETURN  
END SUBROUTINE wrf_quilt_open_for_write_commit

SUBROUTINE wrf_quilt_open_for_read ( FileName , Comm_compute, Comm_io, SysDepInfo, &
                               DataHandle , Status )





  IMPLICIT NONE
  CHARACTER *(*), INTENT(IN)  :: FileName
  INTEGER ,       INTENT(IN)  :: Comm_compute , Comm_io
  CHARACTER *(*), INTENT(IN)  :: SysDepInfo
  INTEGER ,       INTENT(OUT) :: DataHandle
  INTEGER ,       INTENT(OUT) :: Status

  CALL wrf_debug ( 50, 'in wrf_quilt_open_for_read' ) 
  DataHandle = -1
  Status = -1
  CALL wrf_error_fatal3("<stdin>",6014,&
"frame/module_io_quilt.F: wrf_quilt_open_for_read not yet supported" )
  RETURN  
END SUBROUTINE wrf_quilt_open_for_read

SUBROUTINE wrf_quilt_inquire_opened ( DataHandle, FileName , FileStatus, Status )





  USE module_wrf_quilt
  IMPLICIT NONE
      integer, parameter  :: WRF_FILE_NOT_OPENED                  = 100
      integer, parameter  :: WRF_FILE_OPENED_NOT_COMMITTED        = 101
      integer, parameter  :: WRF_FILE_OPENED_FOR_WRITE            = 102
      integer, parameter  :: WRF_FILE_OPENED_FOR_READ             = 103
      integer, parameter  :: WRF_REAL                             = 104
      integer, parameter  :: WRF_DOUBLE                           = 105
      integer, parameter  :: WRF_FLOAT=WRF_REAL
      integer, parameter  :: WRF_INTEGER                          = 106
      integer, parameter  :: WRF_LOGICAL                          = 107
      integer, parameter  :: WRF_COMPLEX                          = 108
      integer, parameter  :: WRF_DOUBLE_COMPLEX                   = 109
      integer, parameter  :: WRF_FILE_OPENED_FOR_UPDATE           = 110


      integer, parameter  :: WRF_FILE_OPENED_AND_COMMITTED        = 102
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER *(*), INTENT(IN)  :: FileName
  INTEGER ,       INTENT(OUT) :: FileStatus
  INTEGER ,       INTENT(OUT) :: Status

  Status = 0

  CALL wrf_debug ( 50, 'in wrf_quilt_inquire_opened' ) 
  IF ( DataHandle .GE. 1 .AND. DataHandle .LE. int_num_handles ) THEN
    IF ( int_handle_in_use( DataHandle ) ) THEN
      IF ( okay_to_write( DataHandle ) ) THEN
        FileStatus = WRF_FILE_OPENED_FOR_WRITE
      ENDIF
    ENDIF
  ENDIF
  Status = 0
  
  RETURN
END SUBROUTINE wrf_quilt_inquire_opened

SUBROUTINE wrf_quilt_inquire_filename ( DataHandle, FileName , FileStatus, Status )










  USE module_wrf_quilt
  IMPLICIT NONE
      integer, parameter  :: WRF_FILE_NOT_OPENED                  = 100
      integer, parameter  :: WRF_FILE_OPENED_NOT_COMMITTED        = 101
      integer, parameter  :: WRF_FILE_OPENED_FOR_WRITE            = 102
      integer, parameter  :: WRF_FILE_OPENED_FOR_READ             = 103
      integer, parameter  :: WRF_REAL                             = 104
      integer, parameter  :: WRF_DOUBLE                           = 105
      integer, parameter  :: WRF_FLOAT=WRF_REAL
      integer, parameter  :: WRF_INTEGER                          = 106
      integer, parameter  :: WRF_LOGICAL                          = 107
      integer, parameter  :: WRF_COMPLEX                          = 108
      integer, parameter  :: WRF_DOUBLE_COMPLEX                   = 109
      integer, parameter  :: WRF_FILE_OPENED_FOR_UPDATE           = 110


      integer, parameter  :: WRF_FILE_OPENED_AND_COMMITTED        = 102
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER *(*), INTENT(OUT) :: FileName
  INTEGER ,       INTENT(OUT) :: FileStatus
  INTEGER ,       INTENT(OUT) :: Status
  CALL wrf_debug ( 50, 'in wrf_quilt_inquire_filename' ) 
  Status = 0
  IF ( DataHandle .GE. 1 .AND. DataHandle .LE. int_num_handles ) THEN
    IF ( int_handle_in_use( DataHandle ) ) THEN
      IF ( okay_to_write( DataHandle ) ) THEN
        FileStatus = WRF_FILE_OPENED_FOR_WRITE
      ELSE
        FileStatus = WRF_FILE_OPENED_NOT_COMMITTED
      ENDIF
    ELSE
        FileStatus = WRF_FILE_NOT_OPENED
    ENDIF
    Status = 0
    FileName = "bogusfornow"
  ELSE
    Status = -1
  ENDIF
  RETURN
END SUBROUTINE wrf_quilt_inquire_filename

SUBROUTINE wrf_quilt_iosync ( DataHandle, Status )



















  USE module_wrf_quilt
  IMPLICIT NONE
  include "mpif.h"
  INTEGER ,       INTENT(IN)  :: DataHandle
  INTEGER ,       INTENT(OUT) :: Status

  INTEGER locsize , itypesize
  INTEGER ierr, tasks_in_group, comm_io_group, dummy, i, rank

  CALL wrf_debug ( 50, 'in wrf_quilt_iosync' ) 


  IF ( associated ( int_local_output_buffer ) ) THEN

    iserver = get_server_id ( DataHandle )
    CALL get_mpi_comm_io_groups( comm_io_group , iserver )

    CALL mpi_x_comm_size( comm_io_group , tasks_in_group , ierr )

    locsize = int_num_bytes_to_write(DataHandle)


    
    reduced = 0
    reduced(1) = locsize 

    IF ( compute_group_master(1,current_id) ) reduced(2) = DataHandle
    CALL MPI_Reduce( reduced, reduced_dummy, 2, MPI_INTEGER, MPI_SUM, tasks_in_group-1, comm_io_group, ierr )


    
    CALL collect_on_comm_debug("module_io_quilt.G",3421, comm_io_group,            &
                          onebyte,                       &
                          int_local_output_buffer, locsize , &
                          dummy, 0 )

    int_local_output_cursor = 1

    DEALLOCATE ( int_local_output_buffer )
    NULLIFY ( int_local_output_buffer )
  ELSE
    CALL wrf_message ("frame/module_io_quilt.F: wrf_quilt_iosync: no buffer allocated")
  ENDIF

  Status = 0
  RETURN
END SUBROUTINE wrf_quilt_iosync

SUBROUTINE wrf_quilt_ioclose ( DataHandle, Status )







  USE module_wrf_quilt
  USE module_timing
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  INTEGER, PARAMETER ::  int_ioexit			=  	     10
  INTEGER, PARAMETER ::  int_open_for_write_begin	=  	     20
  INTEGER, PARAMETER ::  int_open_for_write_commit	=  	     30
  INTEGER, PARAMETER ::  int_open_for_read 		=  	     40
  INTEGER, PARAMETER ::  int_inquire_opened 		=  	     60
  INTEGER, PARAMETER ::  int_inquire_filename 		=  	     70
  INTEGER, PARAMETER ::  int_iosync 			=  	     80
  INTEGER, PARAMETER ::  int_ioclose 			=  	     90
  INTEGER, PARAMETER ::  int_next_time 			=  	    100
  INTEGER, PARAMETER ::  int_set_time 			=  	    110
  INTEGER, PARAMETER ::  int_next_var 			=  	    120
  INTEGER, PARAMETER ::  int_dom_ti_real 		=  	    140
  INTEGER, PARAMETER ::  int_dom_ti_double 		=  	    160
  INTEGER, PARAMETER ::  int_dom_ti_integer 		=  	    180
  INTEGER, PARAMETER ::  int_dom_ti_logical 		=  	    200
  INTEGER, PARAMETER ::  int_dom_ti_char 		=  	    220
  INTEGER, PARAMETER ::  int_dom_td_real 		=  	    240
  INTEGER, PARAMETER ::  int_dom_td_double 		=  	    260
  INTEGER, PARAMETER ::  int_dom_td_integer 		=  	    280
  INTEGER, PARAMETER ::  int_dom_td_logical 		=  	    300
  INTEGER, PARAMETER ::  int_dom_td_char 		=  	    320
  INTEGER, PARAMETER ::  int_var_ti_real 		=  	    340
  INTEGER, PARAMETER ::  int_var_ti_double 		=  	    360
  INTEGER, PARAMETER ::  int_var_ti_integer 		=  	    380
  INTEGER, PARAMETER ::  int_var_ti_logical 		=  	    400
  INTEGER, PARAMETER ::  int_var_ti_char 		=  	    420
  INTEGER, PARAMETER ::  int_var_td_real 		=  	    440
  INTEGER, PARAMETER ::  int_var_td_double 		=  	    460
  INTEGER, PARAMETER ::  int_var_td_integer 		=  	    480
  INTEGER, PARAMETER ::  int_var_td_logical 		=  	    500
  INTEGER, PARAMETER ::  int_var_td_char 		=  	    520
  INTEGER, PARAMETER ::  int_field 			=  	    530
  INTEGER, PARAMETER ::  int_var_info 			=  	    540
  INTEGER, PARAMETER ::  int_noop 			=  	    550
  INTEGER ,       INTENT(IN)  :: DataHandle
  INTEGER ,       INTENT(OUT) :: Status
  INTEGER i, itypesize, tasks_in_group, comm_io_group, ierr
  REAL dummy


  CALL wrf_debug ( 50, 'in wrf_quilt_ioclose' ) 
  CALL MPI_TYPE_SIZE( MPI_INTEGER, itypesize, ierr )




  IF ( compute_group_master(1,current_id) )THEN
     CALL int_gen_handle_header( hdrbuf, hdrbufsize, itypesize, &
                                 DataHandle, int_ioclose )
  ELSE
     CALL int_gen_noop_header( hdrbuf, hdrbufsize, itypesize )
  ENDIF

  iserver = get_server_id ( DataHandle )
  CALL get_mpi_comm_io_groups( comm_io_group , iserver )

  CALL mpi_x_comm_size( comm_io_group , tasks_in_group , ierr )


  
  reduced = 0


  IF ( compute_group_master(1,current_id) ) reduced(2) = DataHandle
  CALL MPI_Reduce( reduced, reduced_dummy, 2, MPI_INTEGER, MPI_SUM, tasks_in_group-1, comm_io_group, ierr )
!!JMTIMING   CALL end_timing("MPI_Reduce in ioclose")


  int_handle_in_use(DataHandle) = .false.
  CALL set_server_id( DataHandle, 0 ) 
  okay_to_write(DataHandle) = .false.
  okay_to_commit(DataHandle) = .false.
  int_local_output_cursor = 1
  int_num_bytes_to_write(DataHandle) = 0
  IF ( associated ( int_local_output_buffer ) ) THEN
    DEALLOCATE ( int_local_output_buffer )
    NULLIFY ( int_local_output_buffer )
  ENDIF

  Status = 0
!!JMTIMING   CALL end_timing( "wrf_quilt_ioclose" )

  RETURN
END SUBROUTINE wrf_quilt_ioclose

SUBROUTINE wrf_quilt_ioexit( Status )





  USE module_wrf_quilt
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  INTEGER, PARAMETER ::  int_ioexit			=  	     10
  INTEGER, PARAMETER ::  int_open_for_write_begin	=  	     20
  INTEGER, PARAMETER ::  int_open_for_write_commit	=  	     30
  INTEGER, PARAMETER ::  int_open_for_read 		=  	     40
  INTEGER, PARAMETER ::  int_inquire_opened 		=  	     60
  INTEGER, PARAMETER ::  int_inquire_filename 		=  	     70
  INTEGER, PARAMETER ::  int_iosync 			=  	     80
  INTEGER, PARAMETER ::  int_ioclose 			=  	     90
  INTEGER, PARAMETER ::  int_next_time 			=  	    100
  INTEGER, PARAMETER ::  int_set_time 			=  	    110
  INTEGER, PARAMETER ::  int_next_var 			=  	    120
  INTEGER, PARAMETER ::  int_dom_ti_real 		=  	    140
  INTEGER, PARAMETER ::  int_dom_ti_double 		=  	    160
  INTEGER, PARAMETER ::  int_dom_ti_integer 		=  	    180
  INTEGER, PARAMETER ::  int_dom_ti_logical 		=  	    200
  INTEGER, PARAMETER ::  int_dom_ti_char 		=  	    220
  INTEGER, PARAMETER ::  int_dom_td_real 		=  	    240
  INTEGER, PARAMETER ::  int_dom_td_double 		=  	    260
  INTEGER, PARAMETER ::  int_dom_td_integer 		=  	    280
  INTEGER, PARAMETER ::  int_dom_td_logical 		=  	    300
  INTEGER, PARAMETER ::  int_dom_td_char 		=  	    320
  INTEGER, PARAMETER ::  int_var_ti_real 		=  	    340
  INTEGER, PARAMETER ::  int_var_ti_double 		=  	    360
  INTEGER, PARAMETER ::  int_var_ti_integer 		=  	    380
  INTEGER, PARAMETER ::  int_var_ti_logical 		=  	    400
  INTEGER, PARAMETER ::  int_var_ti_char 		=  	    420
  INTEGER, PARAMETER ::  int_var_td_real 		=  	    440
  INTEGER, PARAMETER ::  int_var_td_double 		=  	    460
  INTEGER, PARAMETER ::  int_var_td_integer 		=  	    480
  INTEGER, PARAMETER ::  int_var_td_logical 		=  	    500
  INTEGER, PARAMETER ::  int_var_td_char 		=  	    520
  INTEGER, PARAMETER ::  int_field 			=  	    530
  INTEGER, PARAMETER ::  int_var_info 			=  	    540
  INTEGER, PARAMETER ::  int_noop 			=  	    550
  INTEGER ,       INTENT(OUT) :: Status
  INTEGER                     :: DataHandle, actual_iserver
  INTEGER i, itypesize, tasks_in_group, comm_io_group, me, ierr 
  REAL dummy

  CALL wrf_debug ( 50, 'in wrf_quilt_ioexit' ) 
  CALL MPI_TYPE_SIZE( MPI_INTEGER, itypesize, ierr )




  IF( compute_group_master(1,current_id) ) THEN
     CALL int_gen_handle_header( hdrbuf, hdrbufsize, itypesize, &
                                 DataHandle, int_ioexit )
  ELSE
     CALL int_gen_noop_header( hdrbuf, hdrbufsize, itypesize )
  END IF

  DO iserver = 1, nio_groups
    if(poll_servers) then

       
       

       call wrf_quilt_find_server(actual_iserver)

       

       
    else
       
       actual_iserver=iserver
    endif

    CALL get_mpi_comm_io_groups( comm_io_group , actual_iserver )

    CALL mpi_x_comm_size( comm_io_group , tasks_in_group , ierr )
    CALL mpi_comm_rank( comm_io_group , me , ierr )


    hdrbufsize = -100 
    reduced = 0
    IF ( me .eq. 0 ) reduced(1) = hdrbufsize 
    CALL MPI_Reduce( reduced, reduced_dummy, 2, MPI_INTEGER, MPI_SUM, tasks_in_group-1, comm_io_group, ierr )

  ENDDO
  Status = 0

  RETURN  
END SUBROUTINE wrf_quilt_ioexit

SUBROUTINE wrf_quilt_get_next_time ( DataHandle, DateStr, Status )





  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*)               :: DateStr
  INTEGER                     :: Status
  RETURN
END SUBROUTINE wrf_quilt_get_next_time

SUBROUTINE wrf_quilt_get_previous_time ( DataHandle, DateStr, Status )





  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*)               :: DateStr
  INTEGER                     :: Status
  RETURN
END SUBROUTINE wrf_quilt_get_previous_time

SUBROUTINE wrf_quilt_set_time ( DataHandle, Data,  Status )





  USE module_wrf_quilt
  USE module_state_description, ONLY: IO_PNETCDF
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  INTEGER, PARAMETER ::  int_ioexit			=  	     10
  INTEGER, PARAMETER ::  int_open_for_write_begin	=  	     20
  INTEGER, PARAMETER ::  int_open_for_write_commit	=  	     30
  INTEGER, PARAMETER ::  int_open_for_read 		=  	     40
  INTEGER, PARAMETER ::  int_inquire_opened 		=  	     60
  INTEGER, PARAMETER ::  int_inquire_filename 		=  	     70
  INTEGER, PARAMETER ::  int_iosync 			=  	     80
  INTEGER, PARAMETER ::  int_ioclose 			=  	     90
  INTEGER, PARAMETER ::  int_next_time 			=  	    100
  INTEGER, PARAMETER ::  int_set_time 			=  	    110
  INTEGER, PARAMETER ::  int_next_var 			=  	    120
  INTEGER, PARAMETER ::  int_dom_ti_real 		=  	    140
  INTEGER, PARAMETER ::  int_dom_ti_double 		=  	    160
  INTEGER, PARAMETER ::  int_dom_ti_integer 		=  	    180
  INTEGER, PARAMETER ::  int_dom_ti_logical 		=  	    200
  INTEGER, PARAMETER ::  int_dom_ti_char 		=  	    220
  INTEGER, PARAMETER ::  int_dom_td_real 		=  	    240
  INTEGER, PARAMETER ::  int_dom_td_double 		=  	    260
  INTEGER, PARAMETER ::  int_dom_td_integer 		=  	    280
  INTEGER, PARAMETER ::  int_dom_td_logical 		=  	    300
  INTEGER, PARAMETER ::  int_dom_td_char 		=  	    320
  INTEGER, PARAMETER ::  int_var_ti_real 		=  	    340
  INTEGER, PARAMETER ::  int_var_ti_double 		=  	    360
  INTEGER, PARAMETER ::  int_var_ti_integer 		=  	    380
  INTEGER, PARAMETER ::  int_var_ti_logical 		=  	    400
  INTEGER, PARAMETER ::  int_var_ti_char 		=  	    420
  INTEGER, PARAMETER ::  int_var_td_real 		=  	    440
  INTEGER, PARAMETER ::  int_var_td_double 		=  	    460
  INTEGER, PARAMETER ::  int_var_td_integer 		=  	    480
  INTEGER, PARAMETER ::  int_var_td_logical 		=  	    500
  INTEGER, PARAMETER ::  int_var_td_char 		=  	    520
  INTEGER, PARAMETER ::  int_field 			=  	    530
  INTEGER, PARAMETER ::  int_var_info 			=  	    540
  INTEGER, PARAMETER ::  int_noop 			=  	    550
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Data
  INTEGER                     :: Status
  INTEGER i, itypesize, tasks_in_group, ierr, comm_io_group
  REAL dummy
  INTEGER                 :: Count
  INTEGER, EXTERNAL       :: use_package

  CALL wrf_debug ( 50, 'in wrf_quilt_set_time' )

  IF ( DataHandle .GE. 1 .AND. DataHandle .LE. int_num_handles ) THEN
    IF ( int_handle_in_use( DataHandle ) ) THEN
      CALL MPI_TYPE_SIZE( MPI_INTEGER, itypesize, ierr )
      Count = 0   



      IF(compute_group_master(1,current_id) )THEN



         CALL int_gen_ti_header_char( hdrbuf, hdrbufsize, itypesize, &
                                      DataHandle, "TIMESTAMP", "", Data, int_set_time )
      ELSE
         CALL int_gen_noop_header( hdrbuf, hdrbufsize, itypesize )
      END IF

      iserver = get_server_id ( DataHandle )
      CALL get_mpi_comm_io_groups( comm_io_group , iserver )
      CALL mpi_x_comm_size( comm_io_group , tasks_in_group , ierr )

      
      reduced = 0
      reduced(1) = hdrbufsize 
      IF ( compute_group_master(1,current_id) ) reduced(2) = DataHandle
      CALL MPI_Reduce( reduced, reduced_dummy, 2, MPI_INTEGER, MPI_SUM, tasks_in_group-1, comm_io_group, ierr )
      
      CALL collect_on_comm_debug("module_io_quilt.G",3693, comm_io_group,            &
                            onebyte,                       &
                            hdrbuf, hdrbufsize , &
                            dummy, 0 )
    ENDIF
  ENDIF

RETURN
END SUBROUTINE wrf_quilt_set_time

SUBROUTINE wrf_quilt_get_next_var ( DataHandle, VarName, Status )






  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*)               :: VarName
  INTEGER                     :: Status
  RETURN
END SUBROUTINE wrf_quilt_get_next_var

SUBROUTINE wrf_quilt_get_dom_ti_real ( DataHandle,Element,   Data, Count, Outcount, Status )











  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  REAL,            INTENT(IN) :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: Outcount
  INTEGER                     :: Status
  CALL wrf_message('wrf_quilt_get_dom_ti_real not supported yet')
RETURN
END SUBROUTINE wrf_quilt_get_dom_ti_real 

SUBROUTINE wrf_quilt_put_dom_ti_real ( DataHandle,Element,   Data, Count,  Status )








  USE module_wrf_quilt
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  INTEGER, PARAMETER ::  int_ioexit			=  	     10
  INTEGER, PARAMETER ::  int_open_for_write_begin	=  	     20
  INTEGER, PARAMETER ::  int_open_for_write_commit	=  	     30
  INTEGER, PARAMETER ::  int_open_for_read 		=  	     40
  INTEGER, PARAMETER ::  int_inquire_opened 		=  	     60
  INTEGER, PARAMETER ::  int_inquire_filename 		=  	     70
  INTEGER, PARAMETER ::  int_iosync 			=  	     80
  INTEGER, PARAMETER ::  int_ioclose 			=  	     90
  INTEGER, PARAMETER ::  int_next_time 			=  	    100
  INTEGER, PARAMETER ::  int_set_time 			=  	    110
  INTEGER, PARAMETER ::  int_next_var 			=  	    120
  INTEGER, PARAMETER ::  int_dom_ti_real 		=  	    140
  INTEGER, PARAMETER ::  int_dom_ti_double 		=  	    160
  INTEGER, PARAMETER ::  int_dom_ti_integer 		=  	    180
  INTEGER, PARAMETER ::  int_dom_ti_logical 		=  	    200
  INTEGER, PARAMETER ::  int_dom_ti_char 		=  	    220
  INTEGER, PARAMETER ::  int_dom_td_real 		=  	    240
  INTEGER, PARAMETER ::  int_dom_td_double 		=  	    260
  INTEGER, PARAMETER ::  int_dom_td_integer 		=  	    280
  INTEGER, PARAMETER ::  int_dom_td_logical 		=  	    300
  INTEGER, PARAMETER ::  int_dom_td_char 		=  	    320
  INTEGER, PARAMETER ::  int_var_ti_real 		=  	    340
  INTEGER, PARAMETER ::  int_var_ti_double 		=  	    360
  INTEGER, PARAMETER ::  int_var_ti_integer 		=  	    380
  INTEGER, PARAMETER ::  int_var_ti_logical 		=  	    400
  INTEGER, PARAMETER ::  int_var_ti_char 		=  	    420
  INTEGER, PARAMETER ::  int_var_td_real 		=  	    440
  INTEGER, PARAMETER ::  int_var_td_double 		=  	    460
  INTEGER, PARAMETER ::  int_var_td_integer 		=  	    480
  INTEGER, PARAMETER ::  int_var_td_logical 		=  	    500
  INTEGER, PARAMETER ::  int_var_td_char 		=  	    520
  INTEGER, PARAMETER ::  int_field 			=  	    530
  INTEGER, PARAMETER ::  int_var_info 			=  	    540
  INTEGER, PARAMETER ::  int_noop 			=  	    550
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  REAL ,          INTENT(IN)  :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: Status

  CHARACTER*132   :: locElement
  INTEGER i, typesize, itypesize, tasks_in_group, ierr, comm_io_group
  REAL dummy



  CALL wrf_debug ( 50, 'in wrf_quilt_put_dom_ti_real' ) 
  CALL MPI_TYPE_SIZE( MPI_INTEGER, itypesize, ierr )
  locElement = Element

  IF ( DataHandle .GE. 1 .AND. DataHandle .LE. int_num_handles ) THEN
    IF ( int_handle_in_use( DataHandle ) ) THEN
      CALL MPI_TYPE_SIZE( MPI_INTEGER, itypesize, ierr )
      CALL MPI_TYPE_SIZE( MPI_REAL, typesize, ierr )

      IF ( compute_group_master(1,current_id) ) THEN
         CALL int_gen_ti_header( hdrbuf, hdrbufsize, itypesize, typesize, &
                                 DataHandle, locElement, Data, Count, int_dom_ti_real )
      ELSE
         CALL int_gen_noop_header( hdrbuf, hdrbufsize, itypesize )
      ENDIF

      iserver = get_server_id ( DataHandle )
      CALL get_mpi_comm_io_groups( comm_io_group , iserver )
      CALL mpi_x_comm_size( comm_io_group , tasks_in_group , ierr )


      
      reduced = 0
      reduced(1) = hdrbufsize 
      IF( compute_group_master(1,current_id) )  reduced(2) = DataHandle
      CALL MPI_Reduce( reduced, reduced_dummy, 2, MPI_INTEGER, MPI_SUM, tasks_in_group-1, comm_io_group, ierr )

      
      CALL collect_on_comm_debug("module_io_quilt.G",3812, comm_io_group,            &
                            onebyte,                       &
                            hdrbuf, hdrbufsize , &
                            dummy, 0 )
    ENDIF
  ENDIF

  Status = 0

RETURN
END SUBROUTINE wrf_quilt_put_dom_ti_real 

SUBROUTINE wrf_quilt_get_dom_ti_double ( DataHandle,Element,   Data, Count, Outcount, Status )











  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  real*8                      :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: OutCount
  INTEGER                     :: Status
  CALL wrf_error_fatal3("<stdin>",6644,&
'wrf_quilt_get_dom_ti_double not supported yet')
RETURN
END SUBROUTINE wrf_quilt_get_dom_ti_double 

SUBROUTINE wrf_quilt_put_dom_ti_double ( DataHandle,Element,   Data, Count,  Status )










  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  REAL*8 ,        INTENT(IN)  :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: Status
  CALL wrf_error_fatal3("<stdin>",6666,&
'wrf_quilt_put_dom_ti_double not supported yet')
RETURN
END SUBROUTINE wrf_quilt_put_dom_ti_double 

SUBROUTINE wrf_quilt_get_dom_ti_integer ( DataHandle,Element,   Data, Count, Outcount, Status )











  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  integer                     :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                      :: OutCount
  INTEGER                     :: Status
  CALL wrf_message('wrf_quilt_get_dom_ti_integer not supported yet')
RETURN
END SUBROUTINE wrf_quilt_get_dom_ti_integer 

SUBROUTINE wrf_quilt_put_dom_ti_integer ( DataHandle,Element,   Data, Count,  Status )








  USE module_wrf_quilt
  USE module_state_description, ONLY: IO_PNETCDF
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  INTEGER, PARAMETER ::  int_ioexit			=  	     10
  INTEGER, PARAMETER ::  int_open_for_write_begin	=  	     20
  INTEGER, PARAMETER ::  int_open_for_write_commit	=  	     30
  INTEGER, PARAMETER ::  int_open_for_read 		=  	     40
  INTEGER, PARAMETER ::  int_inquire_opened 		=  	     60
  INTEGER, PARAMETER ::  int_inquire_filename 		=  	     70
  INTEGER, PARAMETER ::  int_iosync 			=  	     80
  INTEGER, PARAMETER ::  int_ioclose 			=  	     90
  INTEGER, PARAMETER ::  int_next_time 			=  	    100
  INTEGER, PARAMETER ::  int_set_time 			=  	    110
  INTEGER, PARAMETER ::  int_next_var 			=  	    120
  INTEGER, PARAMETER ::  int_dom_ti_real 		=  	    140
  INTEGER, PARAMETER ::  int_dom_ti_double 		=  	    160
  INTEGER, PARAMETER ::  int_dom_ti_integer 		=  	    180
  INTEGER, PARAMETER ::  int_dom_ti_logical 		=  	    200
  INTEGER, PARAMETER ::  int_dom_ti_char 		=  	    220
  INTEGER, PARAMETER ::  int_dom_td_real 		=  	    240
  INTEGER, PARAMETER ::  int_dom_td_double 		=  	    260
  INTEGER, PARAMETER ::  int_dom_td_integer 		=  	    280
  INTEGER, PARAMETER ::  int_dom_td_logical 		=  	    300
  INTEGER, PARAMETER ::  int_dom_td_char 		=  	    320
  INTEGER, PARAMETER ::  int_var_ti_real 		=  	    340
  INTEGER, PARAMETER ::  int_var_ti_double 		=  	    360
  INTEGER, PARAMETER ::  int_var_ti_integer 		=  	    380
  INTEGER, PARAMETER ::  int_var_ti_logical 		=  	    400
  INTEGER, PARAMETER ::  int_var_ti_char 		=  	    420
  INTEGER, PARAMETER ::  int_var_td_real 		=  	    440
  INTEGER, PARAMETER ::  int_var_td_double 		=  	    460
  INTEGER, PARAMETER ::  int_var_td_integer 		=  	    480
  INTEGER, PARAMETER ::  int_var_td_logical 		=  	    500
  INTEGER, PARAMETER ::  int_var_td_char 		=  	    520
  INTEGER, PARAMETER ::  int_field 			=  	    530
  INTEGER, PARAMETER ::  int_var_info 			=  	    540
  INTEGER, PARAMETER ::  int_noop 			=  	    550
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  INTEGER ,       INTENT(IN)  :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: Status

  CHARACTER*132   :: locElement
  INTEGER i, typesize, itypesize, tasks_in_group, ierr, comm_io_group
  REAL dummy
  INTEGER, EXTERNAL :: use_package



  locElement = Element

  CALL wrf_debug ( 50, 'in wrf_quilt_put_dom_ti_integer' ) 


  IF ( DataHandle .GE. 1 .AND. DataHandle .LE. int_num_handles ) THEN
    IF ( int_handle_in_use( DataHandle ) ) THEN
      CALL MPI_TYPE_SIZE( MPI_INTEGER, itypesize, ierr )
      CALL MPI_TYPE_SIZE( MPI_INTEGER, typesize, ierr )



      IF ( compute_group_master(1,current_id) )THEN
         CALL int_gen_ti_header( hdrbuf, hdrbufsize, itypesize, typesize, &
                                 DataHandle, locElement, Data, Count,     &
                                 int_dom_ti_integer )
      ELSE
         CALL int_gen_noop_header( hdrbuf, hdrbufsize, itypesize )
      ENDIF

      iserver = get_server_id ( DataHandle )
      CALL get_mpi_comm_io_groups( comm_io_group , iserver )
      CALL mpi_x_comm_size( comm_io_group , tasks_in_group , ierr )


      
      reduced = 0
      reduced(1) = hdrbufsize 
      IF ( compute_group_master(1,current_id) ) reduced(2) = DataHandle
      CALL MPI_Reduce( reduced, reduced_dummy, 2, MPI_INTEGER, MPI_SUM, tasks_in_group-1, comm_io_group, ierr )


      
      CALL collect_on_comm_debug("module_io_quilt.G",3973, comm_io_group,            &
                            onebyte,                       &
                            hdrbuf, hdrbufsize , &
                            dummy, 0 )
    ENDIF
  ENDIF
  CALL wrf_debug ( 50, 'returning from wrf_quilt_put_dom_ti_integer' ) 


RETURN
END SUBROUTINE wrf_quilt_put_dom_ti_integer 

SUBROUTINE wrf_quilt_get_dom_ti_logical ( DataHandle,Element,   Data, Count, Outcount, Status )











  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  logical                     :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                      :: OutCount
  INTEGER                     :: Status

RETURN
END SUBROUTINE wrf_quilt_get_dom_ti_logical 

SUBROUTINE wrf_quilt_put_dom_ti_logical ( DataHandle,Element,   Data, Count,  Status )










  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  logical ,            INTENT(IN) :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: Status

  INTEGER i
  INTEGER one_or_zero(Count)

  DO i = 1, Count
    IF ( Data(i) ) THEN
      one_or_zero(i) = 1
    ELSE
      one_or_zero(i) = 0
    ENDIF
  ENDDO

  CALL wrf_quilt_put_dom_ti_integer ( DataHandle,Element,   one_or_zero, Count,  Status )
RETURN
END SUBROUTINE wrf_quilt_put_dom_ti_logical 

SUBROUTINE wrf_quilt_get_dom_ti_char ( DataHandle,Element,   Data,  Status )










  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*)               :: Data
  INTEGER                     :: Status
  CALL wrf_message('wrf_quilt_get_dom_ti_char not supported yet')
RETURN
END SUBROUTINE wrf_quilt_get_dom_ti_char 

SUBROUTINE wrf_quilt_put_dom_ti_char ( DataHandle, Element,  Data,  Status )








  USE module_wrf_quilt
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  INTEGER, PARAMETER ::  int_ioexit			=  	     10
  INTEGER, PARAMETER ::  int_open_for_write_begin	=  	     20
  INTEGER, PARAMETER ::  int_open_for_write_commit	=  	     30
  INTEGER, PARAMETER ::  int_open_for_read 		=  	     40
  INTEGER, PARAMETER ::  int_inquire_opened 		=  	     60
  INTEGER, PARAMETER ::  int_inquire_filename 		=  	     70
  INTEGER, PARAMETER ::  int_iosync 			=  	     80
  INTEGER, PARAMETER ::  int_ioclose 			=  	     90
  INTEGER, PARAMETER ::  int_next_time 			=  	    100
  INTEGER, PARAMETER ::  int_set_time 			=  	    110
  INTEGER, PARAMETER ::  int_next_var 			=  	    120
  INTEGER, PARAMETER ::  int_dom_ti_real 		=  	    140
  INTEGER, PARAMETER ::  int_dom_ti_double 		=  	    160
  INTEGER, PARAMETER ::  int_dom_ti_integer 		=  	    180
  INTEGER, PARAMETER ::  int_dom_ti_logical 		=  	    200
  INTEGER, PARAMETER ::  int_dom_ti_char 		=  	    220
  INTEGER, PARAMETER ::  int_dom_td_real 		=  	    240
  INTEGER, PARAMETER ::  int_dom_td_double 		=  	    260
  INTEGER, PARAMETER ::  int_dom_td_integer 		=  	    280
  INTEGER, PARAMETER ::  int_dom_td_logical 		=  	    300
  INTEGER, PARAMETER ::  int_dom_td_char 		=  	    320
  INTEGER, PARAMETER ::  int_var_ti_real 		=  	    340
  INTEGER, PARAMETER ::  int_var_ti_double 		=  	    360
  INTEGER, PARAMETER ::  int_var_ti_integer 		=  	    380
  INTEGER, PARAMETER ::  int_var_ti_logical 		=  	    400
  INTEGER, PARAMETER ::  int_var_ti_char 		=  	    420
  INTEGER, PARAMETER ::  int_var_td_real 		=  	    440
  INTEGER, PARAMETER ::  int_var_td_double 		=  	    460
  INTEGER, PARAMETER ::  int_var_td_integer 		=  	    480
  INTEGER, PARAMETER ::  int_var_td_logical 		=  	    500
  INTEGER, PARAMETER ::  int_var_td_char 		=  	    520
  INTEGER, PARAMETER ::  int_field 			=  	    530
  INTEGER, PARAMETER ::  int_var_info 			=  	    540
  INTEGER, PARAMETER ::  int_noop 			=  	    550
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: Data
  INTEGER                     :: Status
  INTEGER i, itypesize, tasks_in_group, ierr, comm_io_group, me
  REAL dummy


  CALL wrf_debug ( 50, 'in wrf_quilt_put_dom_ti_char' ) 


  IF ( DataHandle .GE. 1 .AND. DataHandle .LE. int_num_handles ) THEN
    IF ( int_handle_in_use( DataHandle ) ) THEN
      CALL MPI_TYPE_SIZE( MPI_INTEGER, itypesize, ierr )



      IF(compute_group_master(1,current_id))THEN
         CALL int_gen_ti_header_char( hdrbuf, hdrbufsize, itypesize, &
                                      DataHandle, Element, "", Data, &
                                      int_dom_ti_char )
      ELSE
         CALL int_gen_noop_header( hdrbuf, hdrbufsize, itypesize )
      END IF

      iserver = get_server_id ( DataHandle )

      CALL get_mpi_comm_io_groups( comm_io_group , iserver )
      CALL mpi_x_comm_size( comm_io_group , tasks_in_group , ierr )
      







      
      reduced_dummy = 0 
      reduced = 0
      reduced(1) = hdrbufsize 
      IF(compute_group_master(1,current_id))    reduced(2) = DataHandle


      CALL MPI_Reduce( reduced, reduced_dummy, 2, MPI_INTEGER, MPI_SUM, tasks_in_group-1, comm_io_group, ierr )


      


      CALL collect_on_comm_debug("module_io_quilt.G",4145, comm_io_group,            &
                            onebyte,                       &
                            hdrbuf, hdrbufsize , &
                            dummy, 0 )

    ENDIF
  ENDIF


RETURN
END SUBROUTINE wrf_quilt_put_dom_ti_char 

SUBROUTINE wrf_quilt_get_dom_td_real ( DataHandle,Element, DateStr,  Data, Count, Outcount, Status )











  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: DateStr
  real                        :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: OutCount
  INTEGER                     :: Status
RETURN
END SUBROUTINE wrf_quilt_get_dom_td_real 

SUBROUTINE wrf_quilt_put_dom_td_real ( DataHandle,Element, DateStr,  Data, Count,  Status )










  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: DateStr
  real ,            INTENT(IN) :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: Status
RETURN
END SUBROUTINE wrf_quilt_put_dom_td_real 

SUBROUTINE wrf_quilt_get_dom_td_double ( DataHandle,Element, DateStr,  Data, Count, Outcount, Status )











  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: DateStr
  real*8                          :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                      :: OutCount
  INTEGER                     :: Status
  CALL wrf_error_fatal3("<stdin>",7047,&
'wrf_quilt_get_dom_td_double not supported yet')
RETURN
END SUBROUTINE wrf_quilt_get_dom_td_double 

SUBROUTINE wrf_quilt_put_dom_td_double ( DataHandle,Element, DateStr,  Data, Count,  Status )










  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: DateStr
  real*8 ,            INTENT(IN) :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: Status
  CALL wrf_error_fatal3("<stdin>",7070,&
'wrf_quilt_put_dom_td_double not supported yet')
RETURN
END SUBROUTINE wrf_quilt_put_dom_td_double 

SUBROUTINE wrf_quilt_get_dom_td_integer ( DataHandle,Element, DateStr,  Data, Count, Outcount, Status )











  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: DateStr
  integer                          :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                      :: OutCount
  INTEGER                     :: Status
RETURN
END SUBROUTINE wrf_quilt_get_dom_td_integer 

SUBROUTINE wrf_quilt_put_dom_td_integer ( DataHandle,Element, DateStr,  Data, Count,  Status )










  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: DateStr
  integer ,            INTENT(IN) :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: Status
RETURN
END SUBROUTINE wrf_quilt_put_dom_td_integer 

SUBROUTINE wrf_quilt_get_dom_td_logical ( DataHandle,Element, DateStr,  Data, Count, Outcount, Status )











  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: DateStr
  logical                          :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                      :: OutCount
  INTEGER                     :: Status
RETURN
END SUBROUTINE wrf_quilt_get_dom_td_logical 

SUBROUTINE wrf_quilt_put_dom_td_logical ( DataHandle,Element, DateStr,  Data, Count,  Status )










  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: DateStr
  logical ,            INTENT(IN) :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: Status
RETURN
END SUBROUTINE wrf_quilt_put_dom_td_logical 

SUBROUTINE wrf_quilt_get_dom_td_char ( DataHandle,Element, DateStr,  Data,  Status )










  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: DateStr
  CHARACTER*(*)               :: Data
  INTEGER                     :: Status
RETURN
END SUBROUTINE wrf_quilt_get_dom_td_char 

SUBROUTINE wrf_quilt_put_dom_td_char ( DataHandle,Element, DateStr,  Data,  Status )










  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: DateStr
  CHARACTER*(*) , INTENT(IN) :: Data
  INTEGER                          :: Status
RETURN
END SUBROUTINE wrf_quilt_put_dom_td_char 

SUBROUTINE wrf_quilt_get_var_ti_real ( DataHandle,Element,  Varname, Data, Count, Outcount, Status )











  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: VarName 
  real                          :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: OutCount
  INTEGER                     :: Status
RETURN
END SUBROUTINE wrf_quilt_get_var_ti_real 

SUBROUTINE wrf_quilt_put_var_ti_real ( DataHandle,Element,  Varname, Data, Count,  Status )










  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: VarName 
  real ,            INTENT(IN) :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: Status
RETURN
END SUBROUTINE wrf_quilt_put_var_ti_real 

SUBROUTINE wrf_quilt_get_var_ti_double ( DataHandle,Element,  Varname, Data, Count, Outcount, Status )











  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: VarName 
  real*8                      :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: OutCount
  INTEGER                     :: Status
  CALL wrf_error_fatal3("<stdin>",7267,&
'wrf_quilt_get_var_ti_double not supported yet')
RETURN
END SUBROUTINE wrf_quilt_get_var_ti_double 

SUBROUTINE wrf_quilt_put_var_ti_double ( DataHandle,Element,  Varname, Data, Count,  Status )










  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: VarName 
  real*8 ,        INTENT(IN) :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: Status
  CALL wrf_error_fatal3("<stdin>",7290,&
'wrf_quilt_put_var_ti_double not supported yet')
RETURN
END SUBROUTINE wrf_quilt_put_var_ti_double 

SUBROUTINE wrf_quilt_get_var_ti_integer ( DataHandle,Element,  Varname, Data, Count, Outcount, Status )











  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: VarName 
  integer                     :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: OutCount
  INTEGER                     :: Status
RETURN
END SUBROUTINE wrf_quilt_get_var_ti_integer 

SUBROUTINE wrf_quilt_put_var_ti_integer ( DataHandle,Element,  Varname, Data, Count,  Status )










  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: VarName 
  integer ,            INTENT(IN) :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: Status
RETURN
END SUBROUTINE wrf_quilt_put_var_ti_integer 

SUBROUTINE wrf_quilt_get_var_ti_logical ( DataHandle,Element,  Varname, Data, Count, Outcount, Status )











  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: VarName 
  logical                     :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: OutCount
  INTEGER                     :: Status
RETURN
END SUBROUTINE wrf_quilt_get_var_ti_logical 

SUBROUTINE wrf_quilt_put_var_ti_logical ( DataHandle,Element,  Varname, Data, Count,  Status )










  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: VarName 
  logical ,            INTENT(IN) :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: Status
RETURN
END SUBROUTINE wrf_quilt_put_var_ti_logical 

SUBROUTINE wrf_quilt_get_var_ti_char ( DataHandle,Element,  Varname, Data,  Status )










  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: VarName 
  CHARACTER*(*)               :: Data
  INTEGER                     :: Status
RETURN
END SUBROUTINE wrf_quilt_get_var_ti_char 

SUBROUTINE wrf_quilt_put_var_ti_char ( DataHandle,Element,  Varname, Data,  Status )









  USE module_wrf_quilt
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  INTEGER, PARAMETER ::  int_ioexit			=  	     10
  INTEGER, PARAMETER ::  int_open_for_write_begin	=  	     20
  INTEGER, PARAMETER ::  int_open_for_write_commit	=  	     30
  INTEGER, PARAMETER ::  int_open_for_read 		=  	     40
  INTEGER, PARAMETER ::  int_inquire_opened 		=  	     60
  INTEGER, PARAMETER ::  int_inquire_filename 		=  	     70
  INTEGER, PARAMETER ::  int_iosync 			=  	     80
  INTEGER, PARAMETER ::  int_ioclose 			=  	     90
  INTEGER, PARAMETER ::  int_next_time 			=  	    100
  INTEGER, PARAMETER ::  int_set_time 			=  	    110
  INTEGER, PARAMETER ::  int_next_var 			=  	    120
  INTEGER, PARAMETER ::  int_dom_ti_real 		=  	    140
  INTEGER, PARAMETER ::  int_dom_ti_double 		=  	    160
  INTEGER, PARAMETER ::  int_dom_ti_integer 		=  	    180
  INTEGER, PARAMETER ::  int_dom_ti_logical 		=  	    200
  INTEGER, PARAMETER ::  int_dom_ti_char 		=  	    220
  INTEGER, PARAMETER ::  int_dom_td_real 		=  	    240
  INTEGER, PARAMETER ::  int_dom_td_double 		=  	    260
  INTEGER, PARAMETER ::  int_dom_td_integer 		=  	    280
  INTEGER, PARAMETER ::  int_dom_td_logical 		=  	    300
  INTEGER, PARAMETER ::  int_dom_td_char 		=  	    320
  INTEGER, PARAMETER ::  int_var_ti_real 		=  	    340
  INTEGER, PARAMETER ::  int_var_ti_double 		=  	    360
  INTEGER, PARAMETER ::  int_var_ti_integer 		=  	    380
  INTEGER, PARAMETER ::  int_var_ti_logical 		=  	    400
  INTEGER, PARAMETER ::  int_var_ti_char 		=  	    420
  INTEGER, PARAMETER ::  int_var_td_real 		=  	    440
  INTEGER, PARAMETER ::  int_var_td_double 		=  	    460
  INTEGER, PARAMETER ::  int_var_td_integer 		=  	    480
  INTEGER, PARAMETER ::  int_var_td_logical 		=  	    500
  INTEGER, PARAMETER ::  int_var_td_char 		=  	    520
  INTEGER, PARAMETER ::  int_field 			=  	    530
  INTEGER, PARAMETER ::  int_var_info 			=  	    540
  INTEGER, PARAMETER ::  int_noop 			=  	    550
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: VarName 
  CHARACTER*(*) , INTENT(IN)  :: Data
  INTEGER                     :: Status
  INTEGER i, itypesize, tasks_in_group, ierr, comm_io_group
  REAL dummy



  CALL wrf_debug ( 50, 'in wrf_quilt_put_var_ti_char' ) 


  IF ( DataHandle .GE. 1 .AND. DataHandle .LE. int_num_handles ) THEN
    IF ( int_handle_in_use( DataHandle ) ) THEN
      CALL MPI_TYPE_SIZE( MPI_INTEGER, itypesize, ierr )

      IF ( compute_group_master(1,current_id) ) THEN
         CALL int_gen_ti_header_char( hdrbuf, hdrbufsize, itypesize, &
                                      DataHandle, TRIM(Element),     &
                                      TRIM(VarName), TRIM(Data), int_var_ti_char )
      ELSE
         CALL int_gen_noop_header( hdrbuf, hdrbufsize, itypesize )
      ENDIF

      iserver = get_server_id ( DataHandle )
      CALL get_mpi_comm_io_groups( comm_io_group , iserver )
      CALL mpi_x_comm_size( comm_io_group , tasks_in_group , ierr )


      
      reduced = 0
      reduced(1) = hdrbufsize 
      IF ( compute_group_master(1,current_id) ) reduced(2) = DataHandle
      CALL MPI_Reduce( reduced, reduced_dummy, 2, MPI_INTEGER, MPI_SUM, tasks_in_group-1, comm_io_group, ierr )

      
      CALL collect_on_comm_debug("module_io_quilt.G",4678, comm_io_group,            &
                            onebyte,                       &
                            hdrbuf, hdrbufsize , &
                            dummy, 0 )
    ENDIF
  ENDIF


RETURN
END SUBROUTINE wrf_quilt_put_var_ti_char 

SUBROUTINE wrf_quilt_get_var_td_real ( DataHandle,Element,  DateStr,Varname, Data, Count, Outcount, Status )











  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: DateStr
  CHARACTER*(*) , INTENT(IN)  :: VarName 
  real                        :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: OutCount
  INTEGER                     :: Status
RETURN
END SUBROUTINE wrf_quilt_get_var_td_real 

SUBROUTINE wrf_quilt_put_var_td_real ( DataHandle,Element,  DateStr,Varname, Data, Count,  Status )










  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: DateStr
  CHARACTER*(*) , INTENT(IN)  :: VarName 
  real ,            INTENT(IN) :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: Status
RETURN
END SUBROUTINE wrf_quilt_put_var_td_real 

SUBROUTINE wrf_quilt_get_var_td_double ( DataHandle,Element,  DateStr,Varname, Data, Count, Outcount, Status )











  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: DateStr
  CHARACTER*(*) , INTENT(IN)  :: VarName 
  real*8                      :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: OutCount
  INTEGER                     :: Status
  CALL wrf_error_fatal3("<stdin>",7565,&
'wrf_quilt_get_var_td_double not supported yet')
RETURN
END SUBROUTINE wrf_quilt_get_var_td_double 

SUBROUTINE wrf_quilt_put_var_td_double ( DataHandle,Element,  DateStr,Varname, Data, Count,  Status )










  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: DateStr
  CHARACTER*(*) , INTENT(IN)  :: VarName 
  real*8 ,            INTENT(IN) :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: Status
  CALL wrf_error_fatal3("<stdin>",7589,&
'wrf_quilt_put_var_td_double not supported yet')
RETURN
END SUBROUTINE wrf_quilt_put_var_td_double 

SUBROUTINE wrf_quilt_get_var_td_integer ( DataHandle,Element,  DateStr,Varname, Data, Count, Outcount,Status)











  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: DateStr
  CHARACTER*(*) , INTENT(IN)  :: VarName 
  integer                     :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: OutCount
  INTEGER                     :: Status
RETURN
END SUBROUTINE wrf_quilt_get_var_td_integer 

SUBROUTINE wrf_quilt_put_var_td_integer ( DataHandle,Element,  DateStr,Varname, Data, Count,  Status )










  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: DateStr
  CHARACTER*(*) , INTENT(IN)  :: VarName 
  integer ,       INTENT(IN)  :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: Status
RETURN
END SUBROUTINE wrf_quilt_put_var_td_integer 

SUBROUTINE wrf_quilt_get_var_td_logical ( DataHandle,Element,  DateStr,Varname, Data, Count, Outcount, Status )











  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: DateStr
  CHARACTER*(*) , INTENT(IN)  :: VarName 
  logical                          :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                      :: OutCount
  INTEGER                     :: Status
RETURN
END SUBROUTINE wrf_quilt_get_var_td_logical 

SUBROUTINE wrf_quilt_put_var_td_logical ( DataHandle,Element,  DateStr,Varname, Data, Count,  Status )










  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: DateStr
  CHARACTER*(*) , INTENT(IN)  :: VarName 
  logical ,            INTENT(IN) :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: Status
RETURN
END SUBROUTINE wrf_quilt_put_var_td_logical 

SUBROUTINE wrf_quilt_get_var_td_char ( DataHandle,Element,  DateStr,Varname, Data,  Status )










  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: DateStr
  CHARACTER*(*) , INTENT(IN)  :: VarName 
  CHARACTER*(*)               :: Data
  INTEGER                     :: Status
RETURN
END SUBROUTINE wrf_quilt_get_var_td_char 

SUBROUTINE wrf_quilt_put_var_td_char ( DataHandle,Element,  DateStr,Varname, Data,  Status )










  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: DateStr
  CHARACTER*(*) , INTENT(IN)  :: VarName 
  CHARACTER*(*) , INTENT(IN) :: Data
  INTEGER                    :: Status
RETURN
END SUBROUTINE wrf_quilt_put_var_td_char 

SUBROUTINE wrf_quilt_read_field ( DataHandle , DateStr , VarName , Field , FieldType , Comm , IOComm, &
                            DomainDesc , MemoryOrder , Stagger , DimNames ,              &
                            DomainStart , DomainEnd ,                                    &
                            MemoryStart , MemoryEnd ,                                    &
                            PatchStart , PatchEnd ,                                      &
                            Status )







  IMPLICIT NONE
  INTEGER ,       INTENT(IN)    :: DataHandle 
  CHARACTER*(*) , INTENT(INOUT) :: DateStr
  CHARACTER*(*) , INTENT(INOUT) :: VarName
  INTEGER ,       INTENT(INOUT) :: Field(*)
  integer                       ,intent(in)    :: FieldType
  integer                       ,intent(inout) :: Comm
  integer                       ,intent(inout) :: IOComm
  integer                       ,intent(in)    :: DomainDesc
  character*(*)                 ,intent(in)    :: MemoryOrder
  character*(*)                 ,intent(in)    :: Stagger
  character*(*) , dimension (*) ,intent(in)    :: DimNames
  integer ,dimension(*)         ,intent(in)    :: DomainStart, DomainEnd
  integer ,dimension(*)         ,intent(in)    :: MemoryStart, MemoryEnd
  integer ,dimension(*)         ,intent(in)    :: PatchStart,  PatchEnd
  integer                       ,intent(out)   :: Status
  Status = 0
RETURN
END SUBROUTINE wrf_quilt_read_field

SUBROUTINE wrf_quilt_write_field ( DataHandle , DateStr , VarName , Field , FieldType , Comm , IOComm,  &
                             DomainDesc , MemoryOrder , Stagger , DimNames ,              &
                             DomainStart , DomainEnd ,                                    &
                             MemoryStart , MemoryEnd ,                                    &
                             PatchStart , PatchEnd ,                                      &
                             Status )


















  USE module_state_description
  USE module_wrf_quilt
  IMPLICIT NONE
  INCLUDE 'mpif.h'
      integer, parameter  :: WRF_FILE_NOT_OPENED                  = 100
      integer, parameter  :: WRF_FILE_OPENED_NOT_COMMITTED        = 101
      integer, parameter  :: WRF_FILE_OPENED_FOR_WRITE            = 102
      integer, parameter  :: WRF_FILE_OPENED_FOR_READ             = 103
      integer, parameter  :: WRF_REAL                             = 104
      integer, parameter  :: WRF_DOUBLE                           = 105
      integer, parameter  :: WRF_FLOAT=WRF_REAL
      integer, parameter  :: WRF_INTEGER                          = 106
      integer, parameter  :: WRF_LOGICAL                          = 107
      integer, parameter  :: WRF_COMPLEX                          = 108
      integer, parameter  :: WRF_DOUBLE_COMPLEX                   = 109
      integer, parameter  :: WRF_FILE_OPENED_FOR_UPDATE           = 110


      integer, parameter  :: WRF_FILE_OPENED_AND_COMMITTED        = 102
  INTEGER ,       INTENT(IN)    :: DataHandle 
  CHARACTER*(*) , INTENT(IN)    :: DateStr
  CHARACTER*(*) , INTENT(IN)    :: VarName

  integer                       ,intent(in)    :: FieldType
  integer                       ,intent(inout) :: Comm
  integer                       ,intent(inout) :: IOComm
  integer                       ,intent(in)    :: DomainDesc
  character*(*)                 ,intent(in)    :: MemoryOrder
  character*(*)                 ,intent(in)    :: Stagger
  character*(*) , dimension (*) ,intent(in)    :: DimNames
  integer ,dimension(*)         ,intent(in)    :: DomainStart, DomainEnd
  integer ,dimension(*)         ,intent(in)    :: MemoryStart, MemoryEnd
  integer ,dimension(*)         ,intent(in)    :: PatchStart,  PatchEnd
  integer                       ,intent(out)   :: Status

  integer ii,jj,kk,myrank

  REAL, DIMENSION( MemoryStart(1):MemoryEnd(1), &
                   MemoryStart(2):MemoryEnd(2), &
                   MemoryStart(3):MemoryEnd(3) ) :: Field
  INTEGER locsize , typesize, itypesize
  INTEGER ierr, tasks_in_group, comm_io_group, dummy, i
  INTEGER, EXTERNAL :: use_package


  CALL wrf_debug ( 50, 'in wrf_quilt_write_field' ) 

  IF ( .NOT. (DataHandle .GE. 1 .AND. DataHandle .LE. int_num_handles) ) THEN
    CALL wrf_error_fatal3("<stdin>",7833,&
"frame/module_io_quilt.F: wrf_quilt_write_field: invalid data handle" )
  ENDIF
  IF ( .NOT. int_handle_in_use( DataHandle ) ) THEN
    CALL wrf_error_fatal3("<stdin>",7837,&
"frame/module_io_quilt.F: wrf_quilt_write_field: DataHandle not opened" )
  ENDIF

  locsize = (PatchEnd(1)-PatchStart(1)+1)* &
            (PatchEnd(2)-PatchStart(2)+1)* &
            (PatchEnd(3)-PatchStart(3)+1)

  CALL mpi_type_size( MPI_INTEGER, itypesize, ierr )
  
  
  IF ( FieldType .EQ. WRF_DOUBLE ) THEN
    CALL mpi_type_size( MPI_DOUBLE_PRECISION, typesize, ierr )
  ELSE IF ( FieldType .EQ. WRF_FLOAT ) THEN
    CALL mpi_type_size( MPI_REAL, typesize, ierr )
  ELSE IF ( FieldType .EQ. WRF_INTEGER ) THEN
    CALL mpi_type_size( MPI_INTEGER, typesize, ierr )
  ELSE IF ( FieldType .EQ. WRF_LOGICAL ) THEN
    CALL mpi_type_size( MPI_LOGICAL, typesize, ierr )
  ENDIF

  IF ( .NOT. okay_to_write( DataHandle ) ) THEN

      
      
      
      

      CALL int_gen_write_field_header ( hdrbuf, hdrbufsize, itypesize, typesize,           &
                               DataHandle , DateStr , VarName , Field , FieldType , Comm , IOComm,  &
                               333933         , MemoryOrder , Stagger , DimNames ,              &   
                               DomainStart , DomainEnd ,                                    &
                               MemoryStart , MemoryEnd ,                                    &
                               PatchStart , PatchEnd )

      int_num_bytes_to_write(DataHandle) = int_num_bytes_to_write(DataHandle) + locsize * typesize + hdrbufsize

      

      iserver = get_server_id ( DataHandle )
      CALL get_mpi_comm_io_groups( comm_io_group , iserver )
      

      CALL mpi_x_comm_size( comm_io_group , tasks_in_group , ierr )




      
      reduced = 0
      reduced(1) = hdrbufsize 
      IF ( compute_group_master(1,current_id) ) reduced(2) = DataHandle
      CALL MPI_Reduce( reduced, reduced_dummy, 2, MPI_INTEGER, MPI_SUM, tasks_in_group-1, comm_io_group, ierr )

      

      CALL collect_on_comm_debug("module_io_quilt.G",5099, comm_io_group,                   &
                            onebyte,                          &
                            hdrbuf, hdrbufsize ,                 &
                            dummy, 0 )

  ELSE

    IF ( .NOT. associated( int_local_output_buffer ) ) THEN
      ALLOCATE ( int_local_output_buffer( (int_num_bytes_to_write( DataHandle )+1)/itypesize ), Stat=ierr )
      IF(ierr /= 0)THEN
         CALL wrf_error_fatal3("<stdin>",7903,&
"frame/module_io_quilt.F: wrf_quilt_write_field: allocate of int_local_output_buffer failed" )
      END IF
      int_local_output_cursor = 1
    ENDIF
      iserver = get_server_id ( DataHandle )

    
    CALL int_gen_write_field_header ( hdrbuf, hdrbufsize, itypesize, typesize,           &
                             DataHandle , DateStr , VarName , Field , FieldType , Comm , IOComm,  &
                             0          , MemoryOrder , Stagger , DimNames ,              &   
                             DomainStart , DomainEnd ,                                    &
                             MemoryStart , MemoryEnd ,                                    &
                             PatchStart , PatchEnd )

    
    
    CALL int_pack_data ( hdrbuf , hdrbufsize , int_local_output_buffer, int_local_output_cursor )

    
    
    CALL int_pack_data ( Field(PatchStart(1):PatchEnd(1),PatchStart(2):PatchEnd(2),PatchStart(3):PatchEnd(3) ), &
                                  locsize * typesize , int_local_output_buffer, int_local_output_cursor )

  ENDIF
  Status = 0


  RETURN
END SUBROUTINE wrf_quilt_write_field

SUBROUTINE wrf_quilt_get_var_info ( DataHandle , VarName , NDim , MemoryOrder , Stagger , &
                              DomainStart , DomainEnd , Status )







  IMPLICIT NONE
  integer               ,intent(in)     :: DataHandle
  character*(*)         ,intent(in)     :: VarName
  integer                               :: NDim
  character*(*)                         :: MemoryOrder
  character*(*)                         :: Stagger
  integer ,dimension(*)                 :: DomainStart, DomainEnd
  integer                               :: Status
RETURN
END SUBROUTINE wrf_quilt_get_var_info

subroutine wrf_quilt_find_server(iserver)

  
  
  

  
  
  

  use module_wrf_quilt, only : in_avail, mpi_comm_avail, mpi_comm_local

  implicit none
  INCLUDE 'mpif.h'
  integer, intent(inout) :: iserver
  integer :: ierr
  character(255) :: message

  call wrf_message('Polling I/O servers...')

  if(in_avail) then
     call mpi_recv(iserver,1,MPI_INTEGER,MPI_ANY_SOURCE,0,mpi_comm_avail,MPI_STATUS_IGNORE,ierr)
     if(ierr/=0) then
        call wrf_error_fatal3("<stdin>",7977,&
'mpi_recv failed in wrf_quilt_find_server')
     endif
  endif

  call mpi_bcast(iserver,1,MPI_INTEGER,0,mpi_comm_local,ierr)
  if(ierr/=0) then
     call wrf_error_fatal3("<stdin>",7984,&
'mpi_bcast failed in wrf_quilt_find_server')
  endif

  write(message,'("I/O server ",I0," is ready for operations.")') iserver
  call wrf_message(message)


end subroutine wrf_quilt_find_server
subroutine wrf_quilt_server_ready()


  
  
  

  
  
  
  

  use module_wrf_quilt, only : mpi_comm_local, in_avail, availrank, mpi_comm_avail

  implicit none
  INCLUDE 'mpif.h'
  integer :: ierr
  character*255 :: message

  write(message,*) 'Entering wrf_quilt_server_ready.'
  call wrf_debug(1,message)

  call mpi_barrier(mpi_comm_local,ierr)
  if(ierr/=0) then
     call wrf_error_fatal3("<stdin>",8017,&
'mpi_barrier failed in wrf_quilt_server_ready')
  endif

  if(in_avail) then
     write(message,'("mpi_ssend ioserver=",I0," in wrf_quilt_server_ready")') availrank
     call wrf_debug(1,message)
     call mpi_ssend(availrank,1,MPI_INTEGER,0,0,mpi_comm_avail,ierr)
     if(ierr/=0) then
        call wrf_error_fatal3("<stdin>",8026,&
'mpi_ssend failed in wrf_quilt_server_ready')
     endif
  endif

  call mpi_barrier(mpi_comm_local,ierr)
  if(ierr/=0) then
     call wrf_error_fatal3("<stdin>",8033,&
'mpi_barrier failed in wrf_quilt_server_ready')
  endif

  write(message,*) 'Leaving wrf_quilt_server_ready.'
  call wrf_debug(1,message)

end subroutine wrf_quilt_server_ready

SUBROUTINE get_mpi_comm_io_groups( retval, isrvr )





      USE module_wrf_quilt
      USE module_dm, ONLY : current_id  
      IMPLICIT NONE
      INTEGER, INTENT(IN ) :: isrvr
      INTEGER, INTENT(OUT) :: retval
      retval = mpi_comm_io_groups(isrvr,current_id)
      RETURN
END SUBROUTINE get_mpi_comm_io_groups

SUBROUTINE get_nio_tasks_in_group( id, retval )





      USE module_wrf_quilt
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: id
      INTEGER, INTENT(OUT) :: retval
      retval = nio_tasks_per_group(id)
      RETURN
END SUBROUTINE get_nio_tasks_in_group

SUBROUTINE collect_on_comm_debug(file,line, comm_io_group,   &
                        sze,                                 &
                        hdrbuf, hdrbufsize ,                 &
                        outbuf, outbufsize                   )
  IMPLICIT NONE
  CHARACTER*(*) file
  INTEGER line
  INTEGER comm_io_group
  INTEGER sze
  INTEGER hdrbuf(*), outbuf(*)
  INTEGER hdrbufsize, outbufsize 


  CALL collect_on_comm( comm_io_group,                       &
                        sze,                                 &
                        hdrbuf, hdrbufsize ,                 &
                        outbuf, outbufsize                   )

  RETURN
END


SUBROUTINE collect_on_comm_debug2(file,line,var,tag,sz,hdr_rec_size, &
                        comm_io_group,                       &
                        sze,                                 &
                        hdrbuf, hdrbufsize ,                 &
                        outbuf, outbufsize                   )
  IMPLICIT NONE
  CHARACTER*(*) file,var
  INTEGER line,tag,sz,hdr_rec_size
  INTEGER comm_io_group
  INTEGER sze
  INTEGER hdrbuf(*), outbuf(*)
  INTEGER hdrbufsize, outbufsize


  CALL collect_on_comm( comm_io_group,                       &
                        sze,                                 &
                        hdrbuf, hdrbufsize ,                 &
                        outbuf, outbufsize                   )

  RETURN
END
