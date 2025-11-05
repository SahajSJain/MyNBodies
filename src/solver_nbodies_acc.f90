module gpu_utils
  !!!! gpu-utils was entirely LLM generated for the purpose of debugging GPU environment detection and timing.
  !!!! It provides routines to check for GPU availability, print system and GPU info, and time GPU operations.
  use omp_lib
  implicit none
  
  private
  public :: check_gpu_environment, print_gpu_summary, time_gpu_operation
  
  ! If the above doesn't work, we need to explicitly declare the interfaces
  interface
    integer function acc_get_num_devices(devicetype)
      integer :: devicetype
    end function acc_get_num_devices
    
    integer function acc_get_device_num(devicetype)
      integer :: devicetype
    end function acc_get_device_num
    
    subroutine acc_set_device_num(devicenum, devicetype)
      integer :: devicenum, devicetype
    end subroutine acc_set_device_num
    
    subroutine acc_init(devicetype)
      integer :: devicetype
    end subroutine acc_init
  end interface
  
  ! Define the device type constant
  integer, parameter :: acc_device_nvidia = 1
  
contains
  
  subroutine check_gpu_environment(gpu_available, num_devices, device_id)
    implicit none
    logical, intent(out) :: gpu_available
    integer, intent(out) :: num_devices, device_id
    
    ! Check for NVIDIA GPUs
    num_devices = acc_get_num_devices(acc_device_nvidia)
    
    if (num_devices > 0) then
      gpu_available = .true.
      device_id = acc_get_device_num(acc_device_nvidia)
      
      ! Set and initialize the device explicitly
      call acc_set_device_num(0, acc_device_nvidia)
      call acc_init(acc_device_nvidia)
    else
      gpu_available = .false.
      device_id = -1
    end if
  end subroutine check_gpu_environment
  
  subroutine print_gpu_summary()
    implicit none
    logical :: gpu_available
    integer :: num_devices, device_id, max_threads
    
    call check_gpu_environment(gpu_available, num_devices, device_id)
    max_threads = omp_get_max_threads()
    
    print *, ""
    print *, "╔══════════════════════════════════════════════════════════════╗"
    print *, "║           SYSTEM AND DEVICE INFORMATION                      ║"
    print *, "╚══════════════════════════════════════════════════════════════╝"
    print *, ""
    print *, "┌─── HOST INFORMATION ─────────────────────────────────────────┐"
    print '(A,I4)', " │ CPU Threads Available: ", max_threads
    call system("echo -n ' │ Hostname: ' && hostname -f")
    call system("echo -n ' │ CPU Model: ' && lscpu | grep 'Model name' | head -1 | cut -d':' -f2 | xargs")
    print *, "└──────────────────────────────────────────────────────────────┘"
    print *, ""
    
    if (gpu_available) then
      print *, "┌─── GPU DEVICE INFORMATION ───────────────────────────────────┐"
      print '(A,I2)', " │ Number of GPUs Found: ", num_devices
      print '(A,I2)', " │ Active GPU Device ID: ", device_id
      print *, " │"
      call system("echo -n ' │ GPU Name: ' && nvidia-smi --query-gpu=name --format=csv,noheader -i 0")
      call system("echo -n ' │ Compute Capability: ' && nvidia-smi --query-gpu=compute_cap --format=csv,noheader -i 0")
      call system("echo -n ' │ GPU Memory: ' && nvidia-smi --query-gpu=memory.total --format=csv,noheader -i 0")
      call system("echo -n ' │ Driver Version: ' && nvidia-smi --query-gpu=driver_version --format=csv,noheader -i 0")
      print *, " │"
      print *, " │ ✓ GPU Acceleration ENABLED"
      print *, "└──────────────────────────────────────────────────────────────┘"
    else
      print *, "┌─── GPU DEVICE INFORMATION ───────────────────────────────────┐"
      print *, " │ ✗ No GPU devices found!"
      print *, " │ Running in CPU-only mode"
      print *, "└──────────────────────────────────────────────────────────────┘"
    end if
    print *, ""
  end subroutine print_gpu_summary
  
  function time_gpu_operation(operation_name) result(elapsed_time)
    implicit none
    character(len=*), intent(in) :: operation_name
    real(8) :: elapsed_time
    real(8) :: start_time
    
    start_time = omp_get_wtime()
    
    ! Your GPU operation would go here
    !$acc wait
    
    elapsed_time = omp_get_wtime() - start_time
    
    if (len_trim(operation_name) > 0) then
      print '(A,A,F8.4,A)', "  GPU ", trim(operation_name), elapsed_time, " seconds"
    end if
  end function time_gpu_operation
  
end module gpu_utils
program AccNbodiesSolver
  use nbodies
  use acc_nbodies
  use omp_lib
  use gpu_utils
  implicit none

  type(NBodies_type) :: NB
  real(8) :: start_time, end_time, elapsed_time
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
  
  ! Move data to GPU with timing
  print *, "=== Initializing GPU Computation ==="
  gpu_start = omp_get_wtime()
  call move_data_to_device_acc(NB)
  !$acc wait
  gpu_end = omp_get_wtime()
  print '(A,F8.4,A)', "  Data transfer to GPU: ", gpu_end - gpu_start, " seconds"
  
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
  
  call softening_length_setup_acc(NB)

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
    call get_global_metrics_acc(NB)
    
    !$acc wait
    gpu_end = omp_get_wtime()
    gpu_total_time = gpu_total_time + (gpu_end - gpu_start)
    gpu_call_count = gpu_call_count + 3

    time = time + dt

    if(mod(istep, ndump) == 0) then
      call update_host_for_dump_acc(NB)
      call dump_data(NB)
    end if

    if (mod(istep, nprint) == 0) then
      call update_host_for_diagnostics_acc(NB)
      call output_diagnostics(NB)
    end if
  end do

  end_time = omp_get_wtime()
  elapsed_time = end_time - start_time

  ! Final summary
  print *, ""
  print *, "=== Simulation Complete ==="
  print '(A,F10.3,A)', "Total wall time: ", end_time - start_time, " seconds"
  print '(A,F10.3,A)', "Time per step: ", (end_time - start_time)/nsteps, " seconds"
  print *, "==========================="
  
  if (gpu_available .and. gpu_call_count > 0) then
    print *, ""
    print *, "┌─── GPU Performance Metrics ──────────────────────────────────┐"
    print '(A,F12.3,A)', " │ Total GPU compute time: ", gpu_total_time, " seconds"
    print '(A,F12.2,A)', " │ GPU utilization:        ", (gpu_total_time/elapsed_time)*100.0, "%"
    print '(A,F12.6,A)', " │ Avg GPU call time:      ", gpu_total_time/real(gpu_call_count,8), " seconds"
    print *, "└──────────────────────────────────────────────────────────────┘"
  end if
  print *, ""
  
  ! Cleanup
  call cleanup_device_acc(NB)
  call finalize_simulation(NB)

end program AccNbodiesSolver