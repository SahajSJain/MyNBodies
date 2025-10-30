program SerialNbodiesSolver
  use nbodies
  use serial_nbodies
  implicit none

  type(NBodies_type) :: NB

  ! Initialize the simulation (this reads input and sets up everything)
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
  call softening_length_setup_serial(NB)
  ! Time integration loop
  do istep = 1, nsteps
    ! Compute gravitational forces
    call compute_forces_serial(NB)

    ! Update positions and velocities
    call update_positions_serial(NB)

    ! Compute global metrics for diagnostics
    call get_global_metrics_serial(NB)

    ! Update simulation time
    time = time + dt
    ! adjust time step based on aarseth acceleration criterion
    call adjust_timestep_acceleration(NB)
    if(mod(istep, ndump) == 0) then
      call dump_data(NB)
    end if
    ! Output diagnostics at specified intervals
    if (mod(istep, nprint) == 0) then
      ! call adjust_timestep_acceleration(NB)
      call output_diagnostics(NB)
    end if
  end do

  ! Finalize the simulation
  call finalize_simulation(NB)

end program SerialNbodiesSolver
