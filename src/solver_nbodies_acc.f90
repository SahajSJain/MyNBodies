program AccNbodiesSolver
  use nbodies
  use acc_nbodies
  use omp_lib
  use gpu_utils
  implicit none

  type(NBodies_type) :: NB
  real(8) :: start_time, end_time, elapsed_time, intermediate_time
  real(8) :: gpu_start, gpu_end, gpu_total_time
  integer :: gpu_call_count
  logical :: gpu_available
  integer :: num_devices, device_id

  ! Print comprehensive GPU information
  call print_gpu_summary()

  ! Quick GPU check for use in the program
  call check_gpu_environment(gpu_available, num_devices, device_id)

  if (.not. gpu_available) then
    print *, "WARNING: Proceeding without GPU acceleration"
    print *, ""
  end if

  ! Initialize timing variables
  gpu_total_time = 0.0d0
  gpu_call_count = 0
  ! Initialize the simulation
  call initialize_simulation(NB)
  total_time_forces = 0.0d0
  total_time_verlet = 0.0d0
  total_time_update_host = 0.0d0
  total_time_upload = 0.0d0
  total_time_metrics = 0.0d0
  ! Move data to GPU with timing
  print *, "=== Initializing GPU Computation ==="
  gpu_start = omp_get_wtime()
  call move_data_to_device_acc(NB)
  !$acc wait
  gpu_end = omp_get_wtime()
  print '(A,F15.8,A)', "  Data transfer to GPU: ", gpu_end - gpu_start, " seconds"

  idump = -1

  ! Initial computations
  gpu_start = omp_get_wtime()
  call compute_forces_acc(NB)
  call get_global_metrics_acc(NB)
  !$acc wait
  gpu_end = omp_get_wtime()
  gpu_total_time = gpu_total_time + (gpu_end - gpu_start)
  gpu_call_count = gpu_call_count + 2

  !$acc kernels present(NB)
  NB%total_energy_initial = NB%total_energy
  !$acc end kernels

  call update_host_for_diagnostics_acc(NB)
  call update_host_for_dump_acc(NB)
  call output_diagnostics(NB)
  call dump_data(NB)

  ! Main simulation loop
  print *, ""
  print *, "=== Starting Main Simulation Loop ==="
  print '(A,I8,A)', "  Total steps: ", nsteps, merge(" (GPU)", " (CPU)", gpu_available)
  print *, ""

  start_time = omp_get_wtime()

  do istep = 1, nsteps
    ! GPU computations with timing for performance tracking
    gpu_start = omp_get_wtime()

    call compute_forces_acc(NB)
    call update_positions_acc(NB)

    !$acc wait
    gpu_end = omp_get_wtime()
    gpu_total_time = gpu_total_time + (gpu_end - gpu_start)
    gpu_call_count = gpu_call_count + 2

    time = time + dt

    if (mod(istep, ndump) == 0 .or. mod(istep, nprint) == 0) then
      gpu_start = omp_get_wtime()
      call get_global_metrics_acc(NB)
      !$acc wait
      gpu_end = omp_get_wtime()
      gpu_total_time = gpu_total_time + (gpu_end - gpu_start)
      gpu_call_count = gpu_call_count + 1
    end if
    if (mod(istep, ndump) == 0) then
      call update_host_for_dump_acc(NB)
      call dump_data(NB)
    end if

    if (mod(istep, nprint) == 0) then
      call update_host_for_diagnostics_acc(NB)
      call output_diagnostics(NB)
      intermediate_time = omp_get_wtime()
      elapsed_time = intermediate_time - start_time
      print *, "   Time Info: "
      print '(A,F15.8,A)', "    Total time till now: ", elapsed_time, " seconds"
      print '(A,F15.8,A)', "    Total Forces time: ", total_time_forces, " seconds"
      print '(A,F15.8,A)', "    Total Verlet Integration time: ", total_time_verlet, " seconds"
      print '(A,F15.8,A)', "    Total Metrics time: ", total_time_metrics, " seconds"
      print '(A,F15.8,A)', "    Total Host Update time: ", total_time_update_host, " seconds"
      print '(A,F15.8,A)', "    Total Data Upload time: ", total_time_upload, " seconds"
      print '(A)', "---------------------------------------------------------------------------"

    end if
  end do

  end_time = omp_get_wtime()
  elapsed_time = end_time - start_time

  ! Final summary
  print *, ""
  print *, "=== Simulation Complete ==="
  print '(A,F15.8,A)', "Total wall time: ", end_time - start_time, " seconds"
  print '(A,F15.8,A)', "Time per step: ", (end_time - start_time)/nsteps, " seconds"
  print '(A,F15.8,A)', "Total Forces time: ", total_time_forces, " seconds"
  print '(A,F15.8,A)', "Total Verlet Integration time: ", total_time_verlet, " seconds"
  print '(A,F15.8,A)', "Total Metrics time: ", total_time_metrics, " seconds"
  print '(A,F15.8,A)', "Total Host Update time: ", total_time_update_host, " seconds"
  print '(A,F15.8,A)', "Total Data Upload time: ", total_time_upload, " seconds"
  print *, "==========================="

  if (gpu_available .and. gpu_call_count > 0) then
    print *, ""
    print *, "┌─── GPU Performance Metrics ──────────────────────────────────┐"
    print '(A,F15.8,A)', " │ Total GPU compute time: ", gpu_total_time, " seconds"
    print '(A,F15.8,A)', " │ GPU utilization:        ", (gpu_total_time/elapsed_time)*100.0, "%"
    print '(A,F15.8,A)', " │ Avg GPU call time:      ", gpu_total_time/real(gpu_call_count, 8), " seconds"
    print *, "└──────────────────────────────────────────────────────────────┘"
  end if
  print *, ""

  ! Cleanup
  call cleanup_device_acc(NB)
  call finalize_simulation(NB)

end program AccNbodiesSolver
