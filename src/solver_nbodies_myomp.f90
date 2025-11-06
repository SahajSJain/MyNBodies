program MyOMPNbodiesSolver
  use nbodies
  use myomp_nbodies
  use serial_nbodies
  use omp_lib
  implicit none

  type(NBodies_type) :: NB
  integer :: num_threads
  real(real64) :: cpu_start, cpu_end, cpu_total_time
  omp_tforces = 0.0_wp
  omp_tverlet = 0.0_wp
  omp_tmetrics = 0.0_wp
  ! Get number of threads from environment or set default
  !$omp parallel
  !$omp master
  num_threads = omp_get_num_threads()
  !$omp end master
  !$omp end parallel
  ! Initialize the simulation (this reads input and sets up everything)
  print *, "=== My OpenMP N-Body Simulation ==="
  print '(A,I4)', "Number of OpenMP threads: ", num_threads
  print *, "Maximum threads available:", omp_get_max_threads()
  print *, "Number of processors:", omp_get_num_procs() 
  print *, "================================"
  call initialize_simulation(NB)

  idump = -1 ; ! so that first dump is at 0 step
  call compute_forces_serial(NB)
  if(Force_Diagnose == 1) then
    call compute_forces_serial_diagnosis(NB)
  end if
  call get_global_metrics_serial(NB)
  NB%total_energy_initial = NB%total_energy
  call output_diagnostics(NB)
  call dump_data(NB) ! dump initial state
  call softening_length_setup_myomp(NB)
  ! Time integration loop
  cpu_start = omp_get_wtime() 
  do istep = 1, nsteps
    ! Compute gravitational forces
    call compute_forces_myomp(NB)
    ! Update positions and velocities
    call update_positions_myomp(NB)

    ! Update simulation time
    time = time + dt
    if(mod(istep, ndump) == 0 .or. mod(istep, nprint) == 0) then
      call get_global_metrics_myomp(NB)
    end if

    if(mod(istep, ndump) == 0) then
      call dump_data(NB)
    end if
    ! Output diagnostics at specified intervals
    if (mod(istep, nprint) == 0) then
      ! call adjust_timestep_acceleration(NB)
      call output_diagnostics(NB)
    end if
  end do
  cpu_end = omp_get_wtime()
  print *, ""
  print *, "=== Simulation Complete ==="
  print '(A,F15.8,A)', "Total wall time: ", cpu_end - cpu_start, " seconds"
  print '(A,F15.8,A)', "Time per step: ", (cpu_end - cpu_start)/nsteps, " seconds"
  print '(A,F15.8,A)', "Time in compute_forces_myomp: ", omp_tforces, " seconds"
  print '(A,F15.8,A)', "Time in update_positions_myomp: ", omp_tverlet, " seconds"
  print '(A,F15.8,A)', "Time in get_global_metrics_myomp: ", omp_tmetrics, " seconds"
  print *, "==========================="

  ! Finalize the simulation
  call finalize_simulation(NB)

end program MyOMPNbodiesSolver
