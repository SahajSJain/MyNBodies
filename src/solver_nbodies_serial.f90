program SerialNbodiesSolver
  use nbodies
  use serial_nbodies
  use omp_lib
  implicit none

  type(NBodies_type) :: NB
  real(8) :: start_time, end_time, elapsed_time  ! Changed to real(8) for double precision
  s_tforces = 0.0_wp
  s_tverlet = 0.0_wp
  s_tmetrics = 0.0_wp
  ! Initialize the simulation (this reads input and sets up everything)
  call initialize_simulation(NB)

  idump = -1; ! so that first dump is at 0 step
  call compute_forces_serial(NB)
  if (Force_Diagnose == 1) then
    call compute_forces_serial_diagnosis(NB)
  end if
  call get_global_metrics_serial(NB)
  NB%total_energy_initial = NB%total_energy
  call output_diagnostics(NB)
  call dump_data(NB) ! dump initial state

  ! Time integration loop
  start_time = omp_get_wtime()  ! Start timing here

  do istep = 1, nsteps
    ! Compute gravitational forces
    call compute_forces_serial(NB)
    ! Update positions and velocities
    call update_positions_serial(NB)

    ! Update simulation time
    time = time + dt
    if (mod(istep, ndump) == 0 .or. mod(istep, nprint) == 0) then
      call get_global_metrics_serial(NB)
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

  end_time = omp_get_wtime()  ! End timing here
  elapsed_time = end_time - start_time

  print *, ""
  print *, "=== Simulation Complete ==="
  print '(A,F15.8,A)', "Total wall time: ", elapsed_time, " seconds"
  print '(A,F15.8,A)', "Time per step: ", elapsed_time/real(nsteps, 8), " seconds"
  print '(A,F15.8,A)', "Steps per second: ", real(nsteps, 8)/elapsed_time, " steps/s"
  print '(A,F15.8,A)', "Time in compute_forces_serial: ", s_tforces, " seconds"
  print '(A,F15.8,A)', "Time in update_positions_serial: ", s_tverlet, " seconds"
  print '(A,F15.8,A)', "Time in get_global_metrics_serial: ", s_tmetrics, " seconds"
  print *, "==========================="

  ! Finalize the simulation
  call finalize_simulation(NB)

end program SerialNbodiesSolver
