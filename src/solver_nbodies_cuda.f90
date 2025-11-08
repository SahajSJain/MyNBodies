program CudaNbodiesSolver
  use nbodies
  use cuda_nbodies
  use omp_lib
  use gpu_utils
  implicit none

  type(NBodies_type) :: NB
  type(NBodies_type_device) :: d_NB
  real(8) :: start_time, end_time, elapsed_time
  logical :: gpu_available
  integer :: num_devices, device_id

  ! Check GPU availability
  call print_gpu_summary()
  call check_gpu_environment(gpu_available, num_devices, device_id)

  if (.not. gpu_available) then
    stop "ERROR: This is a CUDA program - GPU is required!"
  end if

  ! Initialize the simulation
  call initialize_simulation(NB)
  call initialize_cuda_timers()

  ! Move data to GPU
  print *, "=== Initializing GPU Computation ==="
  call allocate_and_copy_data_to_device_cuda(NB, d_NB)

  idump = -1

  ! Initial computations
  call compute_forces_cuda(NB, d_NB)
  call get_global_metrics_cuda(NB, d_NB)
  call update_host_for_diagnostics_cuda(NB, d_NB)
  
  call update_host_for_dump_cuda(NB, d_NB)
  call output_diagnostics(NB)
  call dump_data(NB)
  NB%total_energy_initial = NB%total_energy
  print *, ""
  print*, " INITIAL TOTAL ENERGY (GPU): ", NB%total_energy_initial 
  print *, "" 
  ! Main simulation loop
  print *, ""
  print *, "=== Starting Main Simulation Loop ==="
  print '(A,I8,A)', "  Total steps: ", nsteps, " (GPU)"
  print *, ""

  start_time = omp_get_wtime()

  do istep = 1, nsteps
    ! Core computation
    call compute_forces_cuda(NB, d_NB)
    call update_positions_cuda(NB, d_NB)
    time = time + dt

    ! Periodic metrics calculation
    if (mod(istep, ndump) == 0 .or. mod(istep, nprint) == 0) then
      call get_global_metrics_cuda(NB, d_NB)
    end if

    ! Periodic data dump
    if (mod(istep, ndump) == 0) then
      call update_host_for_dump_cuda(NB, d_NB)
      call dump_data(NB)
    end if

    ! Periodic diagnostics
    if (mod(istep, nprint) == 0) then
      call update_host_for_diagnostics_cuda(NB, d_NB)
      call output_diagnostics(NB)
    end if
  end do

  end_time = omp_get_wtime()
  elapsed_time = end_time - start_time

  ! Final summary
  print *, ""
  print *, "=== Simulation Complete ==="
  print '(A,F15.8,A)', "Total wall time: ", elapsed_time, " seconds"
  print '(A,F15.8,A)', "Time per step: ", elapsed_time/nsteps, " seconds"
  print *, ""
  
  ! Print CUDA timing breakdown
  call print_cuda_timers()

  ! Cleanup
  call cleanup_device_cuda(d_NB)
  call finalize_simulation(NB)

end program CudaNbodiesSolver