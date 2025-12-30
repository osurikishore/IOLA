module mpi_interface

  !===================================================================

  !$$$ Module documentation block
  !

  !$$$

  !===================================================================

  use kinds

  !-------------------------------------------------------------------

  implicit none

  !-------------------------------------------------------------------

  ! Define necessary include files

  include "mpif.h"

  !-------------------------------------------------------------------

  ! Define all global MPI variables

  character                                              :: mpi_nodename(mpi_max_processor_name)
  character                                              :: mpi_noderequest
  integer(kind=4),           dimension(:),   allocatable :: mpi_ranks
  integer(kind=4)                                        :: mpi_errorstatus(mpi_status_size)
  integer(kind=4)                                        :: mpi_masternode
  integer(kind=4)                                        :: mpi_slavenode
  integer(kind=4)                                        :: mpi_ierror
  integer(kind=4)                                        :: mpi_ierrorcode
  integer(kind=4)                                        :: mpi_procid
  integer(kind=4)                                        :: mpi_nprocs
  integer(kind=4)                                        :: mpi_node_source
  integer(kind=4)                                        :: mpi_node_destination
  integer(kind=4)                                        :: mpi_loopcount
  integer(kind=4)                                        :: mpi_request
  integer(kind=4)                                        :: mpi_group_user
  integer(kind=4)                                        :: mpi_group_nprocs
  integer(kind=4)                                        :: mpi_group_procid
  integer(kind=4)                                        :: mpi_group_begin
  integer(kind=4)                                        :: mpi_group_end

  !-------------------------------------------------------------------

contains

  !===================================================================
  
  ! mpi_interface_initialize.f90:

  !-------------------------------------------------------------------

  subroutine mpi_interface_initialize()

    ! Initialize MPI session

    call mpi_init(mpi_ierror)

    ! Define rank for all nodes requested by user

    call mpi_comm_rank(mpi_comm_world,mpi_procid,mpi_ierror)

    ! Define the total number of nodes requested by user

    call mpi_comm_size(mpi_comm_world,mpi_nprocs,mpi_ierror)

    ! Define global variables

    mpi_masternode = 0

    ! Initialize global variables


  end subroutine mpi_interface_initialize

  !===================================================================

  ! mpi_interface_terminate.f90:

  !-------------------------------------------------------------------

  subroutine mpi_interface_terminate()

    ! Terminate MPI session

    call mpi_finalize(mpi_ierror)
    
  end subroutine mpi_interface_terminate

  !===================================================================

  ! mpi_interface_define_comm.f90:

  !-------------------------------------------------------------------

  subroutine mpi_interface_define_comm()

    ! Define variables computed within routine

    integer(kind=4),     dimension(:),   allocatable :: mpi_processes
    integer(kind=4)                                  :: mpi_worldgroup
    integer(kind=4)                                  :: mpi_newgroup

    ! Define counting variables

    integer                                          :: i, j, k
    integer                                          :: count

    !=================================================================

    ! Compute local variable

    mpi_group_nprocs = (mpi_group_end - mpi_group_begin) + 1

    ! Allocate memory for local variable

    if(.not. allocated(mpi_processes))                                 &
         & allocate(mpi_processes(mpi_group_nprocs))

    !-----------------------------------------------------------------

    ! Define local variable

    mpi_worldgroup = mpi_comm_world

    ! Initialize local variable

    count = 0

    ! Loop through each processor and define local variable

    do k = 1, mpi_nprocs

       ! Define local variable accordingly

       if(k .le. mpi_group_end .and. k .ge. mpi_group_begin) then

          ! Update local variable

          count = count + 1

          ! Define local variable

          mpi_processes(count) = k - 1

       end if ! if(k .le. mpi_group_end .and. k .ge. mpi_group_begin)

    end do ! do k = 1, mpi_nprocs

    ! Define local variables

    call mpi_comm_group(mpi_comm_world,mpi_worldgroup,mpi_ierror)
    call mpi_group_incl(mpi_worldgroup,mpi_group_nprocs,               &
         & mpi_processes(1:mpi_group_nprocs),mpi_newgroup,mpi_ierror)
    call mpi_comm_create(mpi_comm_world,mpi_newgroup,mpi_group_user,   &
         & mpi_ierror)

    !-----------------------------------------------------------------

    ! Deallocate memory for local variable

    if(allocated(mpi_processes)) deallocate(mpi_processes)

    !=================================================================

    ! Return calculated values

    return

    !=================================================================

  end subroutine mpi_interface_define_comm

  !===================================================================

  ! mpi_interface_deadlock.f90:

  !------------------------------------------------------------------- 

  subroutine mpi_interface_deadlock()

    real(r_kind)                                     :: mpi_sx
    real(r_kind)                                     :: mpi_rx
    integer(kind=4)                                  :: mpi_sproc
    integer(kind=4)                                  :: mpi_rproc

    !=================================================================

    ! Define local variables

    mpi_sproc = mod(mpi_procid + 1,mpi_nprocs)
    mpi_rproc = mod(mpi_procid + mpi_nprocs - 1,mpi_nprocs)
    mpi_sx    = real(mpi_procid)
    mpi_rx    = -1.0

    ! Check local variable and proceed accordingly

    if(mod(mpi_procid,2) .eq. 0) then

       ! Send message to nearest neighbor

       call mpi_send(mpi_sx,1,mpi_real,mpi_sproc,0,mpi_comm_world,     &
            & mpi_ierror)

       ! Receive message from nearest neightbor

       call mpi_recv(mpi_rx,1,mpi_real,mpi_rproc,0,mpi_comm_world,     &
            & mpi_errorstatus,mpi_ierror)

    else   ! if(mod(mpi_procid,2) .eq. 0)

       ! Receive message on task

       call mpi_recv(mpi_rx,1,mpi_real,mpi_rproc,0,mpi_comm_world,     &
            & mpi_errorstatus,mpi_ierror)

       ! Send message to nearest neighbor

       call mpi_send(mpi_sx,1,mpi_real,mpi_sproc,0,mpi_comm_world,     &
            & mpi_ierror)

    end if ! if(mod(mpi_procid,2) .eq. 0)

    ! Print message to user

    write(6,500) mpi_procid, mpi_sx, mpi_rx

    call mpi_finalize(mpi_ierror)

    !================================================================= 

    ! Define format statements

500 format('MPI_INTERFACE_DEADLOCK: Processor ID: ', i6, ' Sent: ', f, &
         & ' Received: ', f)

    !=================================================================

  end subroutine mpi_interface_deadlock

  !===================================================================

end module mpi_interface
