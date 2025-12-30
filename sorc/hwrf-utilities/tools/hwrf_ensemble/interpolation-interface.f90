module interpolation_interface

  !=======================================================================

  ! Define associated modules and subroutines

  !-----------------------------------------------------------------------

  use constants
  use kinds

  !-----------------------------------------------------------------------

  use kdtree2_module
  use mpi_interface
  use namelist
  use netcdf

  !-----------------------------------------------------------------------

  implicit none

  !-----------------------------------------------------------------------

  ! Define all data and structure types for routine; these variables
  ! are variables required by the subroutines within this module

  type interpgrid
     type(kdtree2),        pointer                                     :: kdtree_grid
     real(r_kind),                      dimension(:,:,:),  allocatable :: weights
     real(r_kind),                      dimension(:),      allocatable :: filter
     real(r_kind),                      dimension(:),      allocatable :: xlong
     real(r_kind),                      dimension(:),      allocatable :: xlat
     real(r_kind),                      dimension(:),      allocatable :: scutoff
     real(r_kind)                                                      :: distthresh
     real(r_kind),                      dimension(:,:),    allocatable :: grdloc
     integer,                           dimension(:,:),    allocatable :: grdnbors
     integer,                           dimension(:),      allocatable :: mpi_count_begin
     integer,                           dimension(:),      allocatable :: mpi_count_end
     integer                                                           :: mpi_maxprocid
     integer                                                           :: ncoords
     integer                                                           :: nsig
     integer                                                           :: npasses
     integer                                                           :: neighbors
     integer                                                           :: xdim
     integer                                                           :: ydim
  end type interpgrid

  type relaxmask
     real(r_kind),                      dimension(:,:),    allocatable :: mask
     real(r_kind),                      dimension(:,:),    allocatable :: radius
     real(r_kind),                      dimension(:,:),    allocatable :: relax
     real(r_kind)                                                      :: scale
     integer                                                           :: xdim
     integer                                                           :: ydim
  end type relaxmask

  ! Define global variables

  real(r_kind),                         dimension(:),      allocatable :: srcgrid_var
  real(r_kind),                         dimension(:),      allocatable :: dstgrid_var
  real(r_kind),                         dimension(:),      allocatable :: blndgrid_b
  real(r_kind),                         dimension(:),      allocatable :: blndgrid_a
  real(r_kind),                         dimension(:),      allocatable :: blndgrid_r

  !-----------------------------------------------------------------------

contains

  !=======================================================================

  ! interpolation_initialize_relaxation_mask.f90:

  !-----------------------------------------------------------------------

  subroutine interpolation_initialize_relaxation_mask(grid)

    ! Define variables passed to routine

    type(relaxmask)                                          :: grid

    !=====================================================================

    ! Deallocate memory for local variables

    call interpolation_cleanup_relaxation_mask(grid)

    ! Allocate memory for local variables

    if(.not. allocated(grid%mask))                                         &
         & allocate(grid%mask(grid%xdim,grid%ydim))
    if(.not. allocated(grid%radius))                                       &
         & allocate(grid%radius(grid%xdim,grid%ydim))
    if(.not. allocated(grid%relax))                                        &
         & allocate(grid%relax(grid%xdim,grid%ydim))

    !=====================================================================

  end subroutine interpolation_initialize_relaxation_mask

  !=======================================================================

  ! interpolation_define_relaxation_mask.f90:

  !-----------------------------------------------------------------------

  subroutine interpolation_define_relaxation_mask(grid)

    ! Define variables passed to routine

    type(relaxmask)                                                :: grid

    ! Define variables computed within routine

    real(r_kind)                                                   :: maxradius
    real(r_kind)                                                   :: minradius
    real(r_kind)                                                   :: scale_radius
    real(r_kind)                                                   :: scale_distance
    integer                                                        :: cenxcoord
    integer                                                        :: cenycoord

    ! Define counting variables

    integer                                                        :: i, j, k

    !=====================================================================

    ! Initialize local variables

    cenxcoord = grid%xdim
    cenycoord = grid%ydim
    minradius = 1.e30

    ! Loop through local variable

    do j = 1, grid%ydim

       ! Loop through local variable

       do i = 1, grid%xdim

          ! Define local variable

          minradius = min(grid%radius(i,j),minradius)

       end do ! do i = 1, grid%xdim

    end do ! do j = 1, grid%ydim

    ! Loop through local variable

    do j = 1, grid%ydim

       ! Loop through local variable

       do i = 1, grid%xdim

          ! Check local variable and proceed accordingly

          if(grid%radius(i,j) .eq. minradius) then

             ! Define local variables

             cenxcoord = i
             cenycoord = j

          end if ! if(grid%radius(i,j) .eq. minradius)

       end do ! do i = 1, grid%xdim

    end do ! do j = 1, grid%ydim

    !---------------------------------------------------------------------

    ! Initialize local variables

    grid%relax = 0.0

    ! Loop through local variable

    do j = cenycoord, 1, -1

       ! Initialize local variable

       maxradius = -1.e30

       ! Loop through local variable

       do i = cenxcoord, 1, -1

          ! Check local variable and proceed accordingly

          if(grid%mask(i,j) .ne. 0.0) then

             ! Define local variables

             maxradius = max(maxradius,grid%radius(i,j))
 
          end if ! if(grid%mask(i,j) .ne. 0.0)
       
       end do ! do i = cenxcoord, 1, -1

      ! Loop through local variable

       do i = cenxcoord, 1, -1

          ! Define local variables

          scale_radius = grid%scale*maxradius

          ! Compute local variables

          grid%relax(i,j) = 1.0 - (grid%radius(i,j)/maxradius)*            &
               & exp(-1.0*(((grid%radius(i,j) - maxradius)*                &
               & (grid%radius(i,j) - maxradius))/(maxradius*maxradius)))

          ! Define local variables

          if(grid%radius(i,j) .lt. scale_radius) grid%relax(i,j) = 1.0

          ! Check local variable and proceed accordingly

          if(grid%radius(i,j) .ge. scale_radius .and. grid%radius(i,j)     &
               & .le. maxradius) then
             
             ! Define local variables

             scale_distance  = maxradius - scale_radius

             ! Compute local variables

             grid%relax(i,j) = 1.0 - ((grid%radius(i,j) - scale_radius)/   &
                  & (scale_distance))*exp(-1.0*(((grid%radius(i,j) -       &
                  & scale_radius) - scale_distance)*((grid%radius(i,j) -   &
                  & scale_radius) - scale_distance))/(scale_distance*      &
                  & scale_distance))

          end if ! if(grid%radius(i,j) .ge. scale_radius .and.             &
                 ! grid%radius(i,j) .le. maxradius)

          ! Define local variables

          if(grid%radius(i,j) .gt. maxradius) grid%relax(i,j) = 0.0

       end do ! do i = cenxcoord, 1, -1

    end do ! do j = cenycoord, 1, -1

    ! Loop through local variable

    do j = cenycoord, 1, -1

       ! Initialize local variable

       maxradius = -1.e30

       ! Loop through local variable

       do i = cenxcoord, grid%xdim

          ! Check local variable and proceed accordingly

          if(grid%mask(i,j) .ne. 0.0) then

             ! Define local variables

             maxradius = max(maxradius,grid%radius(i,j))
 
          end if ! if(grid%mask(i,j) .ne. 0.0)
       
       end do ! do i = cenxcoord, grid%xdim

       ! Loop through local variable

       do i = cenxcoord, grid%xdim

          ! Define local variables

          scale_radius = grid%scale*maxradius

          ! Compute local variables

          grid%relax(i,j) = 1.0 - (grid%radius(i,j)/maxradius)*            &
               & exp(-1.0*(((grid%radius(i,j) - maxradius)*                &
               & (grid%radius(i,j) - maxradius))/(maxradius*maxradius)))

          ! Define local variables

          if(grid%radius(i,j) .lt. scale_radius) grid%relax(i,j) = 1.0

          ! Check local variable and proceed accordingly

          if(grid%radius(i,j) .ge. scale_radius .and. grid%radius(i,j)     &
               & .le. maxradius) then
             
             ! Define local variables

             scale_distance  = maxradius - scale_radius

             ! Compute local variables

             grid%relax(i,j) = 1.0 - ((grid%radius(i,j) - scale_radius)/   &
                  & (scale_distance))*exp(-1.0*(((grid%radius(i,j) -       &
                  & scale_radius) - scale_distance)*((grid%radius(i,j) -   &
                  & scale_radius) - scale_distance))/(scale_distance*      &
                  & scale_distance))

          end if ! if(grid%radius(i,j) .ge. scale_radius .and.             &
                 ! grid%radius(i,j) .le. maxradius)

          ! Define local variables

          if(grid%radius(i,j) .gt. maxradius) grid%relax(i,j) = 0.0

       end do ! do i = cenxcoord, grid%xdim

    end do ! do j = cenycoord, 1, -1

    ! Loop through local variable

    do j = cenycoord, grid%ydim

       ! Initialize local variable

       maxradius = -1.e30

       ! Loop through local variable

       do i = cenxcoord, 1, -1

          ! Check local variable and proceed accordingly

          if(grid%mask(i,j) .ne. 0.0) then

             ! Define local variables

             maxradius = max(maxradius,grid%radius(i,j))
 
          end if ! if(grid%mask(i,j) .ne. 0.0)
       
       end do ! do i = cenxcoord, 1, -1

       ! Loop through local variable

       do i = cenxcoord, 1, -1

          ! Define local variables

          scale_radius = grid%scale*maxradius

          ! Compute local variables

          grid%relax(i,j) = 1.0 - (grid%radius(i,j)/maxradius)*            &
               & exp(-1.0*(((grid%radius(i,j) - maxradius)*                &
               & (grid%radius(i,j) - maxradius))/(maxradius*maxradius)))

          ! Define local variables

          if(grid%radius(i,j) .lt. scale_radius) grid%relax(i,j) = 1.0

          ! Check local variable and proceed accordingly

          if(grid%radius(i,j) .ge. scale_radius .and. grid%radius(i,j)     &
               & .le. maxradius) then
             
             ! Define local variables

             scale_distance  = maxradius - scale_radius

             ! Compute local variables

             grid%relax(i,j) = 1.0 - ((grid%radius(i,j) - scale_radius)/   &
                  & (scale_distance))*exp(-1.0*(((grid%radius(i,j) -       &
                  & scale_radius) - scale_distance)*((grid%radius(i,j) -   &
                  & scale_radius) - scale_distance))/(scale_distance*      &
                  & scale_distance))

          end if ! if(grid%radius(i,j) .ge. scale_radius .and.             &
                 ! grid%radius(i,j) .le. maxradius)

          ! Define local variables

          if(grid%radius(i,j) .gt. maxradius) grid%relax(i,j) = 0.0

       end do ! do i = cenxcoord, 1, -1

    end do ! do j = cenycoord, grid%ydim

    ! Loop through local variable

    do j = cenycoord, grid%ydim

       ! Initialize local variable

       maxradius = -1.e30

       ! Loop through local variable

       do i = cenxcoord, grid%xdim

          ! Check local variable and proceed accordingly

          if(grid%mask(i,j) .ne. 0.0) then

             ! Define local variables

             maxradius = max(maxradius,grid%radius(i,j))
 
          end if ! if(grid%mask(i,j) .ne. 0.0)
       
       end do ! do i = cenxcoord, grid%xdim

       ! Loop through local variable

       do i = cenxcoord, grid%xdim

          ! Define local variables

          scale_radius = grid%scale*maxradius

          ! Compute local variables

          grid%relax(i,j) = 1.0 - (grid%radius(i,j)/maxradius)*            &
               & exp(-1.0*(((grid%radius(i,j) - maxradius)*                &
               & (grid%radius(i,j) - maxradius))/(maxradius*maxradius)))

          ! Define local variables

          if(grid%radius(i,j) .lt. scale_radius) grid%relax(i,j) = 1.0

          ! Check local variable and proceed accordingly

          if(grid%radius(i,j) .ge. scale_radius .and. grid%radius(i,j)     &
               & .le. maxradius) then
             
             ! Define local variables

             scale_distance  = maxradius - scale_radius

             ! Compute local variables

             grid%relax(i,j) = 1.0 - ((grid%radius(i,j) - scale_radius)/   &
                  & (scale_distance))*exp(-1.0*(((grid%radius(i,j) -       &
                  & scale_radius) - scale_distance)*((grid%radius(i,j) -   &
                  & scale_radius) - scale_distance))/(scale_distance*      &
                  & scale_distance))

          end if ! if(grid%radius(i,j) .ge. scale_radius .and.             &
                 ! grid%radius(i,j) .le. maxradius)

          ! Define local variables

          if(grid%radius(i,j) .gt. maxradius) grid%relax(i,j) = 0.0

       end do ! do i = cenxcoord, grid%xdim

    end do ! do j = cenycoord, grid%ydim

    ! Define local variables

    where(grid%relax .gt. 1.0) grid%relax = 1.0
    where(grid%relax .lt. 0.0) grid%relax = 0.0

    !=====================================================================

  end subroutine interpolation_define_relaxation_mask

  !=======================================================================

  ! interpolation_cleanup_relaxation_mask.f90:

  !-----------------------------------------------------------------------

  subroutine interpolation_cleanup_relaxation_mask(grid)

    ! Define variables passed to routine

    type(relaxmask)                                          :: grid

    !=====================================================================

    ! Deallocate memory for local variables

    if(allocated(grid%mask))   deallocate(grid%mask)
    if(allocated(grid%radius)) deallocate(grid%radius)
    if(allocated(grid%relax))  deallocate(grid%relax)

    !=====================================================================

  end subroutine interpolation_cleanup_relaxation_mask

  !=======================================================================

  ! interpolation_initialize_task_balance.f90:

  !-----------------------------------------------------------------------

  subroutine interpolation_initialize_task_balance(grid)

    ! Define variables passed to routine

    type(interpgrid)                                                     :: grid

    ! Define variables computed within routine

    integer                                                              :: mpi_count_interval

    ! Define counting variables

    integer                                                              :: i, j, k, l
    integer                                                              :: count

    !=====================================================================

    ! Allocate memory for local variables

    if(.not. allocated(grid%mpi_count_begin))                              &
         & allocate(grid%mpi_count_begin((mpi_nprocs)))
    if(.not. allocated(grid%mpi_count_end))                                &
         & allocate(grid%mpi_count_end((mpi_nprocs)))

    !---------------------------------------------------------------------

    ! If on master (root) node (task), define problem and broadcast
    ! variables to each slave (compute) node (task)

    if(mpi_procid .eq. mpi_masternode) then

       ! Initialize counting variable

       count = 1

       ! Compute local variable

       mpi_count_interval = grid%ncoords/(mpi_nprocs - 1)

       ! Initialize local variables

       grid%mpi_count_begin    = 0
       grid%mpi_count_end      = 0
       grid%mpi_count_begin(1) = 1
       grid%mpi_count_end(1)   = grid%mpi_count_begin(1) +                 &
            & mpi_count_interval

       ! Loop through total number of processors

       do l = 2, mpi_nprocs - 1

          ! Define local variables

          grid%mpi_count_begin(l) = grid%mpi_count_end(l-1) + 1
          grid%mpi_count_end(l)   = grid%mpi_count_begin(l) +              &
               & mpi_count_interval

          ! Check local variables and proceed accordingly

          if(grid%mpi_count_begin(l) .gt. grid%ncoords) then

             ! Define local variables

             grid%mpi_count_begin(l) = grid%ncoords
             grid%mpi_count_end(l)   = grid%ncoords
             grid%mpi_maxprocid      = l

             ! Define exit from loop

             goto 1000

          end if ! if(grid%mpi_count_begin(l) .gt. grid%ncoords)

          ! Check local variables and proceed accordingly

          if(grid%mpi_count_end(l) .gt. grid%ncoords) then

             ! Define local variables

             grid%mpi_count_end(l) = grid%ncoords
             grid%mpi_maxprocid    = l

             ! Define exit from loop

             goto 1000

          end if ! if(grid%mpi_count_end(l) .gt. grid%ncoords)

       end do ! do l = 2, (mpi_nprocs - 1)

       ! Define exit from loop

1000   continue

       ! Loop through local variable and proceed accordingly

       do l = 1, grid%mpi_maxprocid

          ! Print message to user

          if(debug) write(6,500) grid%ncoords, l, grid%mpi_count_begin(l), &
               & grid%mpi_count_end(l)

       end do ! do l = 1, grid%mpi_maxprocid

    end if ! if(mpi_procid .eq. mpi_masternode)

    ! Enable the root task to catch up from I/O and calculations

    call mpi_barrier(mpi_comm_world,mpi_ierror)

    !---------------------------------------------------------------------

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)

    call mpi_bcast(grid%mpi_count_begin,mpi_nprocs,mpi_integer,            &
         & mpi_masternode,mpi_comm_world,mpi_ierror)
    call mpi_bcast(grid%mpi_count_end,mpi_nprocs,mpi_integer,              &
         & mpi_masternode,mpi_comm_world,mpi_ierror)
    call mpi_bcast(grid%mpi_maxprocid,1,mpi_integer,mpi_masternode,        &
         & mpi_comm_world,mpi_ierror)

    !=====================================================================

    ! Define format statements

500 format('INTERPOLATION_INITIALIZE_TASK_BALANCE: (grid size/',           &
         & 'task ID/tile min/tile max) : ', i9, i6, i9, i9)

    !=====================================================================

  end subroutine interpolation_initialize_task_balance

  !=======================================================================

  ! interpolation_cleanup_task_balance.f90:

  !-----------------------------------------------------------------------

  subroutine interpolation_cleanup_task_balance(grid)

    ! Define variables passed to routine

    type(interpgrid)                                                     :: grid

    !=====================================================================

    ! Deallocate memory for local variable

    if(allocated(grid%mpi_count_begin))                                    &
         & deallocate(grid%mpi_count_begin)
    if(allocated(grid%mpi_count_end))                                      &
         & deallocate(grid%mpi_count_end)

    !=====================================================================

  end subroutine interpolation_cleanup_task_balance

  !=======================================================================

  ! interpolation_initialize_grid.f90:

  !-----------------------------------------------------------------------

  subroutine interpolation_initialize_grid(grid)

    ! Define variables passed to routine

    type(interpgrid)                                                     :: grid

    !=====================================================================

    ! Deallocate memory for local variables

    call interpolation_cleanup_grid(grid)

    !---------------------------------------------------------------------

    ! Allocate memory for local variables

    if(.not. allocated(grid%grdloc))                                       &
         & allocate(grid%grdloc(3,grid%ncoords))
    if(.not. allocated(grid%filter))                                       &
         & allocate(grid%filter(grid%ncoords))
    if(.not. allocated(grid%xlong))                                        &
         & allocate(grid%xlong(grid%ncoords))
    if(.not. allocated(grid%xlat))                                         &
         & allocate(grid%xlat(grid%ncoords))
    if(.not. allocated(grid%grdnbors))                                     &
         & allocate(grid%grdnbors(grid%ncoords,grid%neighbors))
    if(.not. allocated(grid%scutoff))                                      &
         & allocate(grid%scutoff(grid%npasses))

    !=====================================================================

  end subroutine interpolation_initialize_grid

  !=======================================================================

  ! interpolation_cleanup_grid.f90:

  !-----------------------------------------------------------------------

  subroutine interpolation_cleanup_grid(grid)

    ! Define variables passed to routine

    type(interpgrid)                                                     :: grid

    !=====================================================================

    ! Deallocate memory for local variables

    if(allocated(grid%grdloc))           deallocate(grid%grdloc)
    if(allocated(grid%weights))          deallocate(grid%weights)
    if(allocated(grid%filter))           deallocate(grid%filter)
    if(allocated(grid%xlong))            deallocate(grid%xlong)
    if(allocated(grid%xlat))             deallocate(grid%xlat)
    if(allocated(grid%grdnbors))         deallocate(grid%grdnbors)
    if(allocated(grid%scutoff))          deallocate(grid%scutoff)
    if(mpi_procid .eq. mpi_masternode .and. associated(grid%kdtree_grid))  &
         & call kdtree2_destroy(grid%kdtree_grid)

    !=====================================================================

  end subroutine interpolation_cleanup_grid

  !=======================================================================

  ! interpolation_barnes_analysis_mpi.f90:

  !-----------------------------------------------------------------------

  subroutine interpolation_barnes_analysis_mpi(srcgrid,dstgrid)

    ! Define variable passed to routine

    type(interpgrid)                                                     :: srcgrid

    ! Define variable returned by routine

    type(interpgrid)                                                     :: dstgrid

    ! Define variables computed within routine

    real(r_kind),               dimension(:),                allocatable :: mpi_dstgrid_var
    real(r_kind),               dimension(:),                allocatable :: workgrid
    real(r_kind)                                                         :: weights_sum

    ! Define counting variables

    integer                                                              :: i, j, k, l

    !=====================================================================

    ! Allocate memory for local variables

    if(.not. allocated(mpi_dstgrid_var))                                   &
         & allocate(mpi_dstgrid_var(dstgrid%ncoords))
    if(.not. allocated(workgrid))                                          &
         & allocate(workgrid(dstgrid%ncoords))

    !---------------------------------------------------------------------

    ! Initialize local variables

    mpi_dstgrid_var = 0.0

    ! Check local variable and proceed accordingly

    if(mpi_procid .ne. mpi_masternode .and. mpi_procid .le.                &
         & dstgrid%mpi_maxprocid) then

       ! Loop through local variable and proceed accordingly
       
       do k = dstgrid%mpi_count_begin(mpi_procid),                         &
            & dstgrid%mpi_count_end(mpi_procid)

          ! Loop through the total number of analysis passes
          
          do j = 1, dstgrid%npasses
             
             ! Initialize local variable
             
             weights_sum = 0.0
             
             ! Loop through total number of neighboring points on
             ! source grid

             do i = 1, dstgrid%neighbors
             
                ! Define local variable
                
                weights_sum = weights_sum + dstgrid%weights(k,i,j)
             
                ! Define local variable
             
                mpi_dstgrid_var(k) = mpi_dstgrid_var(k) +                  &
                     & (srcgrid_var(dstgrid%grdnbors(k,i)) -               &
                     & workgrid(k))*dstgrid%weights(k,i,j)
                
             end do ! do i = 1, dstgrid%neighbors

             ! Define local variable accordingly

             if(weights_sum .gt. barnes_weights_threshold .and. j .eq.     &
                  & 1) then

                ! Define local variable

                mpi_dstgrid_var(k) = mpi_dstgrid_var(k)/weights_sum

             else  ! if(weights_sum .gt. barnes_weights_threshold          &
                   ! .and. k .eq. 1)
                
                ! Define local variable

                mpi_dstgrid_var(k) = srcgrid_var(dstgrid%grdnbors(k,1))

             end if ! if(weights_sum .gt. barnes_weights_threshold         &
                    ! .and. k .eq. 1)
             
             ! Define local variable
             
             workgrid(k) = mpi_dstgrid_var(k)
                
          end do ! do j = 1, dstgrid%npasses

       end do ! do k = dstgrid%mpi_count_begin(mpi_procid),                &
              ! dstgrid%mpi_count_end(mpi_procid)

    endif ! if(mpi_procid .ne. mpi_masternode .and. mpi_procid             &
          ! .le. dstgrid%mpi_maxprocid)

    !---------------------------------------------------------------------

    ! Define local variable

    call mpi_reduce(mpi_dstgrid_var(1:dstgrid%ncoords),                    &
         & dstgrid_var(1:dstgrid%ncoords),dstgrid%ncoords,mpi_real,        &
         & mpi_sum,mpi_masternode,mpi_comm_world,mpi_ierror)

    !---------------------------------------------------------------------

    ! Deallocate memory for local variable

    if(allocated(mpi_dstgrid_var)) deallocate(mpi_dstgrid_var)
    if(allocated(workgrid))        deallocate(workgrid)

    !=====================================================================

    ! Return calculated values

    return

    !=====================================================================

  end subroutine interpolation_barnes_analysis_mpi

  !=======================================================================

  ! interpolation_nearest_neighbor_mpi.f90:

  !-----------------------------------------------------------------------

  subroutine interpolation_nearest_neighbor_mpi(srcgrid,dstgrid)

    ! Define variable passed to routine

    type(interpgrid)                                                     :: srcgrid

    ! Define variable returned by routine

    type(interpgrid)                                                     :: dstgrid

    ! Define variables computed within routine

    real(r_kind),               dimension(:),                allocatable :: mpi_dstgrid_var

    ! Define counting variables

    integer                                                              :: i, j, k, l

    !=====================================================================

    ! Allocate memory for local variable

    if(.not. allocated(mpi_dstgrid_var))                                   &
         & allocate(mpi_dstgrid_var(dstgrid%ncoords))

    !---------------------------------------------------------------------

    ! Initialize local variables

    mpi_dstgrid_var = 0.0

    ! Check local variable and proceed accordingly
       
    if(mpi_procid .ne. mpi_masternode .and. mpi_procid .le.                &
         & dstgrid%mpi_maxprocid) then

       ! Loop through local variable and proceed accordingly
       
       do k = dstgrid%mpi_count_begin(mpi_procid),                         &
            & dstgrid%mpi_count_end(mpi_procid)

          ! Define local variable

          mpi_dstgrid_var(k) = srcgrid_var(dstgrid%grdnbors(k,1))

       end do ! do k = dstgrid%mpi_count_begin(mpi_procid),                &
              ! dstgrid%mpi_count_end(mpi_procid)

    endif ! if(mpi_procid .ne. mpi_masternode .and. mpi_procid .le.        &
          ! dstgrid%mpi_maxprocid)

    !---------------------------------------------------------------------

    ! Define local variable

    call mpi_reduce(mpi_dstgrid_var(1:dstgrid%ncoords),                    &
         & dstgrid_var(1:dstgrid%ncoords),dstgrid%ncoords,mpi_real,        &
         & mpi_sum,mpi_masternode,mpi_comm_world,mpi_ierror)

    !---------------------------------------------------------------------

    ! Deallocate memory for local variable

    if(allocated(mpi_dstgrid_var)) deallocate(mpi_dstgrid_var)

    !=====================================================================

    ! Return calculated values

    return

    !=====================================================================

  end subroutine interpolation_nearest_neighbor_mpi

  !=======================================================================

  ! interpolation_define_kdtree_mpi.f90:

  !-----------------------------------------------------------------------

  subroutine interpolation_define_kdtree_mpi(grid)

    ! Define variables passed to routine

    type(interpgrid)                                         :: grid

    ! Define variables computed within routine

    real(r_kind),              dimension(:,:),   allocatable :: mpi_grdloc

    ! Define counting variables

    integer                                                  :: i, j, k

    !=====================================================================

    ! Initialize local variables

    call init_constants_derived()

    !---------------------------------------------------------------------

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)

    call mpi_bcast(grid%xlong,grid%ncoords,mpi_real,mpi_masternode,        &
         & mpi_comm_world,mpi_ierror)
    call mpi_bcast(grid%xlat,grid%ncoords,mpi_real,mpi_masternode,         &
         & mpi_comm_world,mpi_ierror)

    !---------------------------------------------------------------------

    ! Allocate memory for local variable

    if(.not. allocated(mpi_grdloc)) allocate(mpi_grdloc(3,grid%ncoords))

    ! Initialize local variables

    mpi_grdloc  = 0.0
    grid%grdloc = 0.0

    ! If on slave (compute) node (task), receive variables, compute
    ! variables, and send variables to master (root) node (task)
    
    if(mpi_procid .ne. mpi_masternode) then

       ! Check local variable and proceed accordingly

       if(grid%mpi_count_begin(mpi_procid) .ne. 0 .and.                    &
            & grid%mpi_count_end(mpi_procid) .ne. 0) then

          ! Loop through local variable and proceed accordingly
       
          do k = grid%mpi_count_begin(mpi_procid),                         &
               & grid%mpi_count_end(mpi_procid)
          
             ! Compute local variables

             mpi_grdloc(1,k) = rearth_equator*cos(grid%xlat(k))*           &
                  & cos(grid%xlong(k))
             mpi_grdloc(2,k) = rearth_equator*cos(grid%xlat(k))*           &
                  & sin(grid%xlong(k))
             mpi_grdloc(3,k) = rearth_equator*sin(grid%xlat(k))

          end do ! do k = grid%mpi_count_begin(mpi_procid),                &
                 ! grid%mpi_count_end(mpi_procid)

       end if ! if(grid%mpi_count_begin(mpi_procid) .ne. 0 .and.           &
              ! grid%mpi_count_end(mpi_procid) .ne. 0)

    end if ! if(mpi_procid .ne. mpi_masternode)

    ! Define local variable

    call mpi_reduce(mpi_grdloc(1:3,1:grid%ncoords),grid%grdloc(1:3,1:      &
         & grid%ncoords),(3*grid%ncoords),mpi_real,mpi_sum,                &
         & mpi_masternode,mpi_comm_world,mpi_ierror)  

    ! Deallocate memory for local variable

    if(allocated(mpi_grdloc)) deallocate(mpi_grdloc)

    !---------------------------------------------------------------------

    ! If on master (root) node (task), define problem and broadcast
    ! variables to each slave (compute) node (task)

    if(mpi_procid .eq. mpi_masternode) then

       ! Initialize local variable

       grid%kdtree_grid => kdtree2_create(grid%grdloc,sort=.true.,         &
            & rearrange=.true.)

    end if ! if(mpi_procid .eq. mpi_masternode)

    ! Enable the root task to catch up from I/O and calculations

    call mpi_barrier(mpi_comm_world,mpi_ierror)

    !=====================================================================

  end subroutine interpolation_define_kdtree_mpi

  !=======================================================================

  ! interpolation_define_weights_mpi.f90:

  !-----------------------------------------------------------------------

  subroutine interpolation_define_weights_mpi(srcgrid,dstgrid)

    ! Define variables passed to routine

    type(interpgrid)                                                     :: srcgrid
    type(interpgrid)                                                     :: dstgrid

    ! Define variables computed within routine

    type(kdtree2_result),       dimension(srcgrid%neighbors)             :: sresults
    real(r_kind),               dimension(:,:,:),            allocatable :: mpi_weights
    real(r_kind),               dimension(:,:),              allocatable :: mpi_sresults_dis
    real(r_kind),               dimension(:,:),              allocatable :: mpi_scutoff
    real(r_kind)                                                         :: mpi_distthresh

    ! Define counting variables

    integer                                                              :: i, j, k, l 

    !=====================================================================

    ! Allocate memory for local variables

    if(.not. allocated(srcgrid%weights))                                   &
         & allocate(srcgrid%weights(srcgrid%ncoords,srcgrid%neighbors,     &
         & srcgrid%npasses))
    if(.not. allocated(dstgrid%weights))                                   & 
         & allocate(dstgrid%weights(dstgrid%ncoords,dstgrid%neighbors,     &
         & dstgrid%npasses))
    if(.not. allocated(mpi_weights))                                       &     
         & allocate(mpi_weights(dstgrid%ncoords,dstgrid%neighbors,         &
         & dstgrid%npasses))
    if(.not. allocated(mpi_scutoff))                                       &
         & allocate(mpi_scutoff(dstgrid%ncoords,dstgrid%npasses))
    if(.not. allocated(mpi_sresults_dis))                                  &
         & allocate(mpi_sresults_dis(dstgrid%ncoords,dstgrid%neighbors))

    !---------------------------------------------------------------------

    ! Initialize local variables

    dstgrid%weights = 0.0
    mpi_weights     = 0.0

    !---------------------------------------------------------------------

    ! If on master (root) node (task), define problem and broadcast
    ! variables to each slave (compute) node (task)
       
    if(mpi_procid .eq. mpi_masternode) then
          
       ! Loop through local variable and proceed accordingly

       do i = 1, dstgrid%ncoords

          ! Define local variable
             
          call kdtree2_n_nearest(tp=srcgrid%kdtree_grid,                   &
               & qv=dstgrid%grdloc(:,i),nn=dstgrid%neighbors,              &
               & results=sresults)
             
          ! Define local variables
             
          mpi_distthresh                          =                        &
               & dstgrid%distthresh
          mpi_scutoff(i,1:dstgrid%npasses)        =                        &
               & dstgrid%scutoff(1:dstgrid%npasses)
          mpi_sresults_dis(i,1:dstgrid%neighbors) =                        &
               & sresults(1:dstgrid%neighbors)%dis
          dstgrid%grdnbors(i,1:dstgrid%neighbors) =                        &
               & sresults(1:dstgrid%neighbors)%idx

       end do ! do i = 1, dstgrid%ncoords
       
    end if ! if(mpi_procid .eq. mpi_masternode)
       
    ! Enable the root task to catch up from I/O and calculations

    call mpi_barrier(mpi_comm_world,mpi_ierror)

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)

    call mpi_bcast(mpi_distthresh,1,mpi_real,mpi_masternode,               &
         & mpi_comm_world,mpi_ierror)
    call mpi_bcast(mpi_scutoff,(dstgrid%ncoords*dstgrid%npasses),          &
         & mpi_real,mpi_masternode,mpi_comm_world,mpi_ierror)
    call mpi_bcast(mpi_sresults_dis,(dstgrid%ncoords*dstgrid%neighbors),   &
         & mpi_real,mpi_masternode,mpi_comm_world,mpi_ierror)
    call mpi_bcast(dstgrid%ncoords,1,mpi_integer,mpi_masternode,           &
         & mpi_comm_world,mpi_ierror)
    call mpi_bcast(dstgrid%neighbors,1,mpi_integer,mpi_masternode,         &
         & mpi_comm_world,mpi_ierror)
    call mpi_bcast(dstgrid%grdnbors,(dstgrid%ncoords*dstgrid%neighbors),   &
         & mpi_integer,mpi_masternode,mpi_comm_world,mpi_ierror)

    ! If on slave (compute) node (task), receive variables, compute
    ! variables, and send variables to master (root) node (task)
    
    if(mpi_procid .ne. mpi_masternode) then

       ! Check local variable and proceed accordingly

       if(dstgrid%mpi_count_begin(mpi_procid) .ne. 0 .and.                 &
            & dstgrid%mpi_count_end(mpi_procid) .ne. 0) then

          ! Loop through local variable and proceed accordingly
       
          do k = dstgrid%mpi_count_begin(mpi_procid),                      &
               & dstgrid%mpi_count_end(mpi_procid)
          
             ! Loop through total number of neighboring points on
             ! destination grid

             do i = 1, dstgrid%neighbors

                ! Loop through the total number of analysis passes to
                ! perform

                do j = 1, dstgrid%npasses

                   ! Compute local variable

                   mpi_weights(k,i,j) = exp((-1.0)*(                       &
                        & sqrt(mpi_sresults_dis(k,i))*                     &
                        & sqrt(mpi_sresults_dis(k,i)))/(4.0*               &
                        & mpi_scutoff(k,j)*mpi_distthresh*                 &
                        & mpi_distthresh))

                end do ! do j = 1, dstgrid%npasses

             end do ! do i = 1, dstgrid%neighbors

          end do ! do k = dstgrid%mpi_count_begin(mpi_procid),             &
                 ! dstgrid%mpi_count_end(mpi_procid)

       endif ! if(dstgrid%mpi_count_begin(mpi_procid) .ne. 0 .and.         &
             ! dstgrid%mpi_count_end(mpi_procid) .ne. 0)

    end if ! if(mpi_procid .ne. mpi_masternode)

    ! Define local variable

    call mpi_reduce(mpi_weights(1:dstgrid%ncoords,1:dstgrid%neighbors,1:   &
         & dstgrid%npasses),dstgrid%weights(1:dstgrid%ncoords,1:           &
         & dstgrid%neighbors,1:dstgrid%npasses),(dstgrid%ncoords*          &
         & dstgrid%neighbors*dstgrid%npasses),mpi_real,mpi_sum,            &
         & mpi_masternode,mpi_comm_world,mpi_ierror)

    !---------------------------------------------------------------------

    ! Deallocate memory for local variable on all compute tasks

    if(allocated(mpi_weights))      deallocate(mpi_weights)
    if(allocated(mpi_sresults_dis)) deallocate(mpi_sresults_dis)
    if(allocated(mpi_scutoff))      deallocate(mpi_scutoff)

    !=====================================================================

    ! Return calculated values

    return

    !=====================================================================

  end subroutine interpolation_define_weights_mpi

  !=======================================================================

  ! interpolation_blend_analysis.f90:

  !-----------------------------------------------------------------------

  subroutine interpolation_blend_analysis(dstgrid)

    ! Define variables passed to routine

    type(interpgrid)                                                     :: dstgrid

    ! Define variables computed within routine
    
    real(r_kind),                 dimension(dstgrid%xdim,dstgrid%ydim)   :: relax_t
    real(r_kind),                 dimension(dstgrid%xdim*dstgrid%ydim)   :: workgrid
    real(r_kind)                                                         :: relax_width
    real(r_kind)                                                         :: relax_center
    integer                                                              :: mpi_xlong_min_idx
    integer                                                              :: mpi_xlat_min_idx
    integer                                                              :: mpi_xlong_max_idx
    integer                                                              :: mpi_xlat_max_idx
    integer                                                              :: tmp_mpi_xlong_min_idx
    integer                                                              :: tmp_mpi_xlat_min_idx
    integer                                                              :: tmp_mpi_xlong_max_idx
    integer                                                              :: tmp_mpi_xlat_max_idx
    integer                                                              :: idx_cnt
    integer                                                              :: idx_ii
    integer                                                              :: idx_jj

    ! Define counting variables

    integer                                                              :: i, j, k

    !=====================================================================

    ! Initialize local variables

    call init_constants(.true.)
    call init_constants_derived()
    mpi_xlong_min_idx = dstgrid%xdim
    mpi_xlong_max_idx = 1
    mpi_xlat_min_idx  = dstgrid%ydim
    mpi_xlat_max_idx  = 1
    idx_cnt           = 0

    ! Loop through local variable

    do j = 1, dstgrid%ydim

       ! Loop through local variable

       do i = 1, dstgrid%xdim

          ! Update local variable

          idx_cnt = idx_cnt + 1

          ! Check local variable and proceed accordingly

          if(dstgrid_var(idx_cnt) .gt. -1.e20) then

             ! Define local variables

             mpi_xlong_min_idx = min(mpi_xlong_min_idx,i)
             mpi_xlong_max_idx = max(mpi_xlong_max_idx,i)
             mpi_xlat_min_idx  = min(mpi_xlat_min_idx,j)
             mpi_xlat_max_idx  = max(mpi_xlat_max_idx,j)

          end if ! if(dstgrid_var .gt. -1.e20)

       end do ! do i = 1, dstgrid%xdim

    end do ! do j = 1, dstgrid%ydim

    ! Define local variables

    dstgrid_var = blndgrid_a

    ! If on master (root) node (task), define problem and broadcast
    ! variables to each slave (compute) node (task)
       
    if(mpi_procid .eq. mpi_masternode) then

       ! Define local variables

       workgrid = blndgrid_a

       ! Loop through local variable

       do k = 1, blndgrid_nsmooth
    
          ! Compute local variables
       
          call interpolation_ninepointsmoother(dstgrid%xdim,dstgrid%ydim,  &
               & workgrid,2,dstgrid%xdim-1,2,dstgrid%ydim-1)

       end do ! do k = 1, blndgrid_nsmooth

       ! Initialize local variables

       relax_t = 0.0

       ! Loop through local variable

       do j = mpi_xlat_min_idx + 1, mpi_xlat_max_idx - 1

          ! Compute local variables

          relax_center = 0.5*(max(mpi_xlong_min_idx - blndgrid_nhalo,1) +  &
               & min(mpi_xlong_min_idx + blndgrid_nhalo,dstgrid%xdim))
          relax_width  = sqrt(real(blndgrid_nhalo))

          ! Loop through local variable
          
          do i = max(mpi_xlong_min_idx - blndgrid_nhalo,1),                &
               & min(mpi_xlong_min_idx + blndgrid_nhalo,dstgrid%xdim)

             ! Compute local variables

             relax_t(i,j) = exp(-((i - relax_center)**2.0)/                &
                  & (2.0*((relax_width))**2.0))

          end do ! do i = max(mpi_xlong_min_idx - blndgrid_nhalo,1),       &
                 ! min(mpi_xlong_min_idx + blndgrid_nhalo,dstgrid%xdim)

       end do ! do j = mpi_xlat_min_idx, mpi_xlat_max_idx

       ! Compute local variables

       dstgrid_var = dstgrid_var*(1.0 - reshape(relax_t,                   &
            & (/dstgrid%ncoords/))) + reshape(relax_t,                     &
            & (/dstgrid%ncoords/))*workgrid

       ! Initialize local variables

       relax_t = 0.0

       ! Loop through local variable

       do j = mpi_xlat_min_idx + 1, mpi_xlat_max_idx - 1

          ! Compute local variables

          relax_center = 0.5*(max(mpi_xlong_max_idx - blndgrid_nhalo,1) +  &
               & min(mpi_xlong_max_idx + blndgrid_nhalo,dstgrid%xdim))
          relax_width  = sqrt(real(blndgrid_nhalo))

          ! Loop through local variable
          
          do i = max(mpi_xlong_max_idx - blndgrid_nhalo,1),                &
               & min(mpi_xlong_max_idx + blndgrid_nhalo,dstgrid%xdim)

             ! Compute local variables

             relax_t(i,j) = exp(-((i - relax_center)**2.0)/                &
                  & (2.0*((relax_width))**2.0))

          end do ! do i = max(mpi_xlong_max_idx - blndgrid_nhalo,1),       &
                 ! min(mpi_xlong_max_idx + blndgrid_nhalo,dstgrid%xdim)

       end do ! do j = mpi_xlat_min_idx, mpi_xlat_max_idx

       ! Compute local variables

       dstgrid_var = dstgrid_var*(1.0 - reshape(relax_t,                   &
            & (/dstgrid%ncoords/))) + reshape(relax_t,                     &
            & (/dstgrid%ncoords/))*workgrid

       ! Initialize local variables

       relax_t = 0.0

       ! Loop through local variable

       do i = mpi_xlong_min_idx + 1, mpi_xlong_max_idx - 1

          ! Compute local variables

          relax_center = 0.5*(max(mpi_xlat_min_idx - blndgrid_nhalo,1) +   &
               & min(mpi_xlat_min_idx + blndgrid_nhalo,dstgrid%ydim))
          relax_width  = sqrt(real(blndgrid_nhalo))

          ! Loop through local variable

          do j = max(mpi_xlat_min_idx - blndgrid_nhalo,1),                 &
               & min(mpi_xlat_min_idx + blndgrid_nhalo,dstgrid%ydim)

             ! Compute local variables

             relax_t(i,j) = exp(-((j - relax_center)**2.0)/                &
                  & (2.0*((relax_width))**2.0))

          end do ! do j = max(mpi_xlat_min_idx - blndgrid_nhalo,1),        &
                 ! min(mpi_xlat_min_idx + blndgrid_nhalo,dstgrid%ydim)

       end do ! do i = mpi_xlong_min_idx, mpi_xlong_max_idx 

       ! Compute local variables

       dstgrid_var = dstgrid_var*(1.0 - reshape(relax_t,                   &
            & (/dstgrid%ncoords/))) + reshape(relax_t,                     &
            & (/dstgrid%ncoords/))*workgrid

       ! Initialize local variables

       relax_t = 0.0

       ! Loop through local variable

       do i = mpi_xlong_min_idx + 1, mpi_xlong_max_idx - 1

          ! Compute local variables

          relax_center = 0.5*(max(mpi_xlat_max_idx - blndgrid_nhalo,1) +   &
               & min(mpi_xlat_max_idx + blndgrid_nhalo,dstgrid%ydim))
          relax_width  = sqrt(real(blndgrid_nhalo))

          ! Loop through local variable

          do j = max(mpi_xlat_max_idx - blndgrid_nhalo,1),                 &
               & min(mpi_xlat_max_idx + blndgrid_nhalo,dstgrid%ydim)

             ! Compute local variables

             relax_t(i,j) = exp(-((j - relax_center)**2.0)/                &
                  & (2.0*((relax_width))**2.0))  

          end do ! do j = max(mpi_xlat_max_idx - blndgrid_nhalo,1),        &
                 ! min(mpi_xlat_max_idx + blndgrid_nhalo,dstgrid%ydim)

       end do ! do i = mpi_xlong_min_idx, mpi_xlong_max_idx

       ! Compute local variables

       dstgrid_var = dstgrid_var*(1.0 - reshape(relax_t,                   &
            & (/dstgrid%ncoords/))) + reshape(relax_t,                     &
            & (/dstgrid%ncoords/))*workgrid

    end if ! if(mpi_procid .eq. mpi_masternode)

    ! Enable the root task to catch up from I/O and calculations

    call mpi_barrier(mpi_comm_world,mpi_ierror)

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)

    call mpi_bcast(dstgrid_var,(dstgrid%ncoords),mpi_real,                 &
         & mpi_masternode,mpi_comm_world,mpi_ierror)

    !=====================================================================

  end subroutine interpolation_blend_analysis

  !=======================================================================

  ! interpolation_ninepointsmoother.f90:

  !-----------------------------------------------------------------------

  subroutine interpolation_ninepointsmoother(xdim,ydim,array,xmin,xmax,    &
       & ymin,ymax)

    ! Define array dimension variables

    integer                                                              :: xdim
    integer                                                              :: ydim

    ! Define variables passed to routine

    real(r_kind),                 dimension(xdim*ydim)                   :: array
    integer                                                              :: xmin
    integer                                                              :: xmax
    integer                                                              :: ymin
    integer                                                              :: ymax

    ! Define variables computed within routine

    real(r_kind),                 dimension(xdim,ydim)                   :: workgrid
    real(r_kind),                 dimension(xdim,ydim)                   :: vargrid

    ! Define counting variables

    integer                                                              :: i, j, k

    !=====================================================================

    ! Define local variables

    vargrid = reshape(array,(/xdim,ydim/))

    ! Initialize local variables

    workgrid = vargrid

    ! Loop through local variable

    do j = ymin, ymax

       ! Loop through local variable

       do i = xmin, xmax

          ! Compute local variables

          workgrid(i,j) = ((vargrid(i-1,j-1) + vargrid(i,j-1) +             &
               & vargrid(i+1,j-1) + vargrid(i-1,j) + vargrid(i,j) +         &
               & vargrid(i+1,j) + vargrid(i-1,j+1) + vargrid(i,j+1) +       &
               & vargrid(i+1,j+1))/9.0)

       end do ! do i = xmin, xmax

    end do ! do j = ymin, ymax

    ! Define local variables

    vargrid = workgrid
    array   = reshape(vargrid,(/xdim*ydim/))

    !=====================================================================

  end subroutine interpolation_ninepointsmoother

  !=======================================================================

  ! define_scaling_coefficients.f90:

  !-----------------------------------------------------------------------

  subroutine define_scaling_coefficients(grid)

    ! Define variable passed to routine

    type(interpgrid)                                         :: grid

    ! Define counting variables

    integer                                                  :: i, j, k

    !=====================================================================

    ! Initialize local variable

    grid%scutoff(1) = 1.0

    ! Loop through total number of threshold cutoff values

    do k = 1, grid%npasses

       ! Compute local variable

       grid%scutoff(k) = grid%scutoff(1)/(10**(real(k-1)))

    end do ! do k = 1, grid%npasses

    !=====================================================================

  end subroutine define_scaling_coefficients

  !=======================================================================

  ! latlondist.f90:

  !-----------------------------------------------------------------------

  subroutine latlondist(xlat1,xlong1,xlat2,xlong2,distance)

    ! Define variables passed to routine

    real(r_kind)                                                    :: xlat1
    real(r_kind)                                                    :: xlong1
    real(r_kind)                                                    :: xlat2
    real(r_kind)                                                    :: xlong2

    ! Define variables returned by routine

    real(r_kind)                                                    :: distance

    ! Define variables computed within routine

    real(r_double)                                                  :: xlat1_deg
    real(r_double)                                                  :: xlong1_deg
    real(r_double)                                                  :: xlat2_deg
    real(r_double)                                                  :: xlong2_deg
    real(r_double)                                                  :: xlat1_rad
    real(r_double)                                                  :: xlong1_rad
    real(r_double)                                                  :: xlat2_rad
    real(r_double)                                                  :: xlong2_rad
    real(r_double)                                                  :: x1
    real(r_double)                                                  :: y1
    real(r_double)                                                  :: z1
    real(r_double)                                                  :: x2
    real(r_double)                                                  :: y2
    real(r_double)                                                  :: z2
    real(r_kind)                                                    :: dr

    !=====================================================================

    ! Initialize local variables

    call init_constants(.true.)
    call init_constants_derived()

    !---------------------------------------------------------------------

    ! Compute local variables 

    xlong1_deg = xlong1
    xlong1_deg = mod(xlong1_deg,360.d0)
    if(xlong1_deg .lt. 0.d0) then
       xlong1_deg = xlong1_deg + 360.d0
    end if ! if(xlong1_deg .lt. 0.d0)
    xlong2_deg = xlong2
    xlong2_deg = mod(xlong2_deg,360.d0)
    if(xlong2_deg .lt. 0.d0) then
       xlong2_deg = xlong2_deg + 360.d0
    end if ! if(xlong2_deg .lt. 0.d0)
    if(xlat1 .lt. xlat2) then
       xlong1_rad = xlong1_deg*deg2rad
       xlat1_rad = (90.d0 - xlat1)*deg2rad
       xlong2_rad = xlong2_deg*deg2rad
       xlat2_rad = (90.d0 - xlat2)*deg2rad
    else if(xlat1 .eq. xlat2 .and. xlong1_deg .le. xlong2_deg) then
       xlong1_rad = xlong1_deg*deg2rad
       xlat1_rad = (90.d0 - xlat1)*deg2rad
       xlong2_rad = xlong2_deg*deg2rad
       xlat2_rad = (90.d0 - xlat2)*deg2rad
    else
       xlong1_rad = xlong2_deg*deg2rad
       xlat1_rad = (90.d0 - xlat2)*deg2rad
       xlong2_rad = xlong1_deg*deg2rad
       xlat2_rad = (90.d0 - xlat1)*deg2rad
    end if ! if(xlat1 .lt. xlat2)

    ! Define local variables

    x1 = sin(xlat1_rad)*cos(xlong1_rad)
    y1 = sin(xlat1_rad)*sin(xlong1_rad)
    z1 = cos(xlat1_rad)
    x2 = sin(xlat2_rad)*cos(xlong2_rad)
    y2 = sin(xlat2_rad)*sin(xlong2_rad)
    z2 = cos(xlat2_rad)

    ! Compute local variables

    dr       = acos(min(1.d0,x1*x2 + y1*y2 + z1*z2))
    distance = dr*rearth_equator

    !=====================================================================

  end subroutine latlondist

  !=======================================================================

end module interpolation_interface
