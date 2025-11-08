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
      print '(A,A,F15.8,A)', "  GPU ", trim(operation_name), elapsed_time, " seconds"
    end if
  end function time_gpu_operation

end module gpu_utils