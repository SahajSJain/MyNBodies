program OMPNbodiesSolver
  use nbodies
  use omp_nbodies
  use omp_lib
  implicit none

  type(NBodies_type) :: NB
  integer :: num_threads
  real(real64) :: start_time, end_time

  ! Get number of threads from environment or set default
  !$omp parallel
  !$omp master
  num_threads = omp_get_num_threads()
  !$omp end master
  !$omp end parallel

  print *, "=== OpenMP N-Body Simulation ==="
  print '(A,I4)', "Number of OpenMP threads: ", num_threads
  print *, "================================"

  ! Initialize the simulation (this reads input and sets up everything)
  call initialize_simulation(NB)

  idump = -1 ! so that first dump is at 0 step
  call compute_forces_omp(NB)
  if (Force_Diagnose == 1) then
    call compute_forces_serial_diagnosis(NB)
  end if
  call get_global_metrics_omp(NB)
  NB%total_energy_initial = NB%total_energy
  call output_diagnostics(NB)
  call dump_data(NB) ! dump initial state
  ! Time integration loop
  start_time = omp_get_wtime()
  do istep = 1, nsteps
    ! Compute gravitational forces
    call compute_forces_omp(NB)
    ! Update positions and velocities
    ! call update_positions_omp(NB)
    call update_positions_omp(NB)
    ! Update simulation time
    time = time + dt
    if (mod(istep, ndump) == 0 .or. mod(istep, nprint) == 0) then
      call get_global_metrics_omp(NB)
    end if

    if (mod(istep, ndump) == 0) then
      call dump_data(NB)
    end if
    ! Output diagnostics at specified intervals
    if (mod(istep, nprint) == 0) then
      ! call adjust_timestep_acceleration(NB)
      call output_diagnostics(NB)
    end if
  end do

  ! End timing
  end_time = omp_get_wtime()

  print *, ""
  print *, "=== Simulation Complete ==="
  print '(A,F10.3,A)', "Total wall time: ", end_time - start_time, " seconds"
  print '(A,F10.3,A)', "Time per step: ", (end_time - start_time)/nsteps, " seconds"
  print *, "==========================="

  ! Finalize the simulation
  call finalize_simulation(NB)

end program OMPNbodiesSolver
