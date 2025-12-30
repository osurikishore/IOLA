program wrf_ensemble_main

  !=====================================================================

  !$$$ PROGRAM DOCUMENTATION BLOCK
  !
  ! ABSTRACT:
  !
  !  This routine will ingest WRF NMM input, output, and/or restart
  !  files, typically produced for the NCEP GSI and ENKF applications,
  !  and can complete the following tasks:
  !
  !  (1) Compute the ensemble mean for the NMM state variables
  !      (if_ensemble_mean = .true.)
  !
  !  (2) Recenter the ensemble about (typically, a hybrid-ensemble) an
  !      analysis; in this instance, if the analysis about which the
  !      ensemble is recentered is of a different resolution, the
  !      respective analysis is interpolated first to ensemble grid
  !      and then the recentering is performed (if_ensemble_recenter =
  !      .true.)
  !
  ! PRGMMR: Winterbottom        
  ! ORG:    ESRL/PSD1       
  ! DATE:   2015-01-02
  !
  ! PROGRAM HISTORY LOG:
  !
  !  2015-01-02 Initial version. Henry R. Winterbottom
  !
  ! LANGUAGE ATTRIBUTES: FORTRAN 90, MPICH
  ! COMPILER ATTRIBUTES: INTEL, MVAPICH2
  !
  ! EXTERNAL I/O ROUTINES:
  !
  !  mpi-interface.f90, namelist.f90
  !
  ! EXTERNAL MODULES:
  !  
  !  constants.f90, kinds.f90, kdtree2.f90, interpolation-interface,
  !  ensemble-mean-interface.f90, ensemble-recenter-interface.f90,
  !  wrf-ensemble-interface.f90
  !
  !$$$

  !=====================================================================

  ! Define associated modules and subroutines

  !---------------------------------------------------------------------

  use kinds

  !---------------------------------------------------------------------

  use wrf_ensemble_interface
  use mpi_interface

  !---------------------------------------------------------------------

  implicit none

  !=====================================================================

  ! Define variables computed within routine

  real(r_kind)                                             :: exectime_start
  real(r_kind)                                             :: exectime_finish

  !=====================================================================

  ! Initialize MPI session

  call mpi_interface_initialize()

  !---------------------------------------------------------------------

  ! If on root (master) task, perform all necessary tasks

  if(mpi_procid .eq. mpi_masternode) then

     ! Define local variable

     call cpu_time(exectime_start)

  end if ! if(mpi_procid .eq. mpi_masternode)

  ! Enable the root task to catch up from I/O and calculations

  call mpi_barrier(mpi_comm_world,mpi_ierror)

  !---------------------------------------------------------------------

  ! Compute local variables

  call wrf_ensemble()

  !---------------------------------------------------------------------

  ! If on root (master) task, perform all necessary tasks

  if(mpi_procid .eq. mpi_masternode) then

     ! Define local variable

     call cpu_time(exectime_finish)

     ! Print message to user

     write(6,500) exectime_finish - exectime_start

  end if ! if(mpi_procid .eq. mpi_masternode)

  ! Enable the root task to catch up from I/O and calculations

  call mpi_barrier(mpi_comm_world,mpi_ierror)

  ! Finalize MPI session

  call mpi_interface_terminate()

  !=====================================================================

  ! Define format statements

500 format('MAIN: Execution time: ', f13.5, ' seconds.')

  !=====================================================================

end program wrf_ensemble_main
