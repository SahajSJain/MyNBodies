program LLMACCNbodiesSolver
  use nbodies
  use llmacc_nbodies
  use serial_nbodies
  use acc_nbodies
  implicit none

  type(NBodies_type) :: NB
  integer :: num_devices, device_id
  real(8) :: start_time, end_time
  
  ! For timing without OpenACC runtime functions
  integer :: count_rate, count_max
  integer(8) :: count_start, count_end
  
  ! Initialize the simulation (this reads input and sets up everything)
  print *, "=== LLMACC (OpenACC) N-Body Simulation ==="
  print *, "=========================================="
  
  call initialize_simulation(NB)

  idump = -1 ! so that first dump is at 0 step
  
  ! Initial force computation on CPU for setup
  call compute_forces_serial(NB)
  if(Force_Diagnose == 1) then
    call compute_forces_serial_diagnosis(NB)
  end if
  call get_global_metrics_serial(NB)
  NB%total_energy_initial = NB%total_energy
  call output_diagnostics(NB)
  call dump_data(NB) ! dump initial state
  call softening_length_setup_llmacc(NB)
  call move_data_to_device_acc(NB)
  ! Time integration loop
  call system_clock(count_start, count_rate, count_max)
  
  ! Create persistent data region for frequently used arrays
  !$acc data copy(NB%pos, NB%vel, NB%pos_0, NB%vel_0, NB%force, NB%accel) &
  !$acc      copyin(NB%mass, NB%mass_inv) &
  !$acc      create(NB%vel_h, NB%potential_energy, NB%kinetic_energy, &
  !$acc             NB%sum_energy, NB%bc_factor)
  
  do istep = 1, nsteps
    ! Compute gravitational forces on GPU
    call compute_forces_llmacc(NB)
    
    ! Update positions and velocities on GPU
    call update_positions_llmacc(NB)

    ! Update simulation time
    time = time + dt
    
    ! Periodically update metrics and output
    if(mod(istep, ndump) == 0 .or. mod(istep, nprint) == 0) then
      !$acc update host(NB%pos, NB%vel, NB%potential_energy, NB%kinetic_energy)
      call get_global_metrics_llmacc(NB)
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
  
  !$acc end data
  
  call system_clock(count_end, count_rate, count_max)
  
  ! Calculate elapsed time
  if (count_end >= count_start) then
    start_time = 0.0d0
    end_time = real(count_end - count_start, kind=8) / real(count_rate, kind=8)
  else
    ! Handle clock wraparound
    start_time = 0.0d0
    end_time = real(count_max - count_start + count_end + 1, kind=8) / real(count_rate, kind=8)
  end if
  
  print *, ""
  print *, "=== Simulation Complete ==="
  print '(A,F10.3,A)', "Total wall time: ", end_time - start_time, " seconds"
  print '(A,F10.3,A)', "Time per step: ", (end_time - start_time)/nsteps, " seconds"
  print '(A,F10.6,A)', "Time per step per body: ", (end_time - start_time)/(nsteps*NB%n_bodies), " seconds"
  print *, "======================================="

  ! Finalize the simulation
  call finalize_simulation(NB)

end program LLMACCNbodiesSolver