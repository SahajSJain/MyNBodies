module cuda_nbodies
  use nbodies
  use omp_lib
  use cudafor
  use openacc
  implicit none
  ! Number of threads per block for CUDA kernels
  integer, parameter :: num_threads = 256
  ! Timing variables
  real(8) :: cuda_time_forces, cuda_time_verlet, cuda_time_update_host, cuda_time_upload, cuda_time_metrics
  real(8) :: t_cuda_start, t_cuda_end
  ! cuda specific variables
  real(kind=wp), device :: d_total_KE, d_total_PE
  real(kind=wp), device :: d_total_COM(3), d_total_COM_vel(3)
  real(kind=wp), device :: d_total_mass, d_total_mass_inv

  integer, device :: d_ndim
  real(kind=wp), device :: d_dt, d_dt_half, d_Gravitational_Constant
  logical, device :: d_REFLECTIVE_BC
  integer, device :: d_iReflective_BC(3)
  real(kind=wp), device :: d_L_bound(3)

  ! Define a device-compatible version of the NBodies_type structure
  type :: NBodies_type_device
    ! Number of bodies
    integer, device, allocatable :: n_bodies(:) ! size 1 integer

    ! DYNAMICAL ARRAYS - BODY PROPERTIES
    ! Position and velocity arrays (n_bodies x ndim)
    real(kind=wp), device, allocatable :: pos(:, :)       ! Current positions
    real(kind=wp), device, allocatable :: pos_0(:, :)     ! Previous positions
    real(kind=wp), device, allocatable :: pos_d(:, :)     ! donor positions to make life easier in compute forces
    real(kind=wp), device, allocatable :: vel(:, :)       ! Current velocities
    real(kind=wp), device, allocatable :: vel_0(:, :)     ! Previous velocities
    real(kind=wp), device, allocatable :: vel_h(:, :)     ! Half-step velocities

    ! Force and acceleration arrays (n_bodies x ndim)
    real(kind=wp), device, allocatable :: accel(:, :)     ! Accelerations
    real(kind=wp), device, allocatable :: force(:, :)     ! Force vectors
    real(kind=wp), device, allocatable :: bc_factor(:, :) ! Boundary condition factors

    ! Mass and geometry arrays (n_bodies)
    real(kind=wp), device, allocatable :: mass(:)        ! Body masses
    real(kind=wp), device, allocatable :: mass_inv(:)    ! Inverse masses
    real(kind=wp), device, allocatable :: radius(:)      ! Body radii
    real(kind=wp), device, allocatable :: radius2(:)      ! Body radii squared
    ! ENERGY AND GLOBAL PROPERTIES
    ! Per-body energy arrays (n_bodies)
    real(kind=wp), device, allocatable :: kinetic_energy(:)    ! KE per body
    real(kind=wp), device, allocatable :: potential_energy(:)  ! PE per body
    real(kind=wp), device, allocatable :: sum_energy(:)        ! Total E per body

    ! (to,from,ndim)
  end type NBodies_type_device

contains

  subroutine allocate_and_copy_data_to_device_cuda(h_NB, d_NB)
    implicit none
    type(NBodies_type), intent(inout) :: h_NB
    type(NBodies_type_device), intent(out) :: d_NB
    integer :: N, nd
    real(8) :: t_start, t_end
    N = h_NB%N_bodies
    nd = ndim
    t_start = omp_get_wtime()

    ! Allocate and set n_bodies
    allocate (d_NB%n_bodies(1))
    d_NB%n_bodies(1) = N
    d_ndim = ndim

    ! ~~~~~~~~~~~~~~~~~~~~~ MOVE ENTIRE STRUCTURE TO DEVICE ~~~~~~~~~~~~~~~~~~~~~~
    ! Now copy all the allocatable arrays within the structure
    ! Arrays that are initialized on host
    allocate (d_NB%pos(N, nd), d_NB%pos_0(N, nd), d_NB%pos_d(N, nd), d_NB%vel(N, nd), d_NB%vel_0(N, nd), &
              d_NB%mass(N), d_NB%mass_inv(N), d_NB%radius(N), d_NB%radius2(N))

    ! Copy data from host to device
    d_NB%pos(1:N, 1:nd) = h_NB%pos(1:N, 1:nd)
    d_NB%pos_0(1:N, 1:nd) = h_NB%pos_0(1:N, 1:nd)
    d_NB%pos_d(1:N, 1:nd) = h_NB%pos(1:N, 1:nd) ! initialize donor positions to current positions
    d_NB%vel(1:N, 1:nd) = h_NB%vel(1:N, 1:nd)
    d_NB%vel_0(1:N, 1:nd) = h_NB%vel_0(1:N, 1:nd)
    d_NB%mass(1:N) = h_NB%mass(1:N)
    d_NB%mass_inv(1:N) = h_NB%mass_inv(1:N)
    d_NB%radius(1:N) = h_NB%radius(1:N)
    d_NB%radius2(1:N) = h_NB%radius2(1:N) ! new array for softening

    ! Arrays that will be computed on device - create them
    allocate (d_NB%vel_h(N, nd), d_NB%accel(N, nd), d_NB%force(N, 3), d_NB%bc_factor(N, nd), &
              d_NB%kinetic_energy(N), d_NB%potential_energy(N), d_NB%sum_energy(N))

    ! ~~~~~~~~~~~~~~~~~~~~~ OTHER GLOBAL VARIABLES ~~~~~~~~~~~~~~~~~~~~~~
    ! Set module-level device variables for scalars
    d_total_mass = h_NB%total_mass
    d_total_mass_inv = h_NB%total_mass_inv
    
    ! Set other module-level device variables
    d_ndim = ndim
    d_Gravitational_Constant = Gravitational_Constant
    d_dt = dt
    d_dt_half = dt_half
    d_REFLECTIVE_BC = REFLECTIVE_BC
    d_iReflective_BC(1:3) = iReflective_BC(1:3)
    d_L_bound(1:3) = L_bound(1:3)

    ! Initialize global accumulators
    d_total_KE = 0.0_wp
    d_total_PE = 0.0_wp
    d_total_COM = 0.0_wp
    d_total_COM_vel = 0.0_wp

    t_end = omp_get_wtime()
    cuda_time_upload = cuda_time_upload + (t_end - t_start)
  end subroutine allocate_and_copy_data_to_device_cuda

  subroutine update_host_for_dump_cuda(h_NB, d_NB)
    implicit none
    type(NBodies_type), intent(inout) :: h_NB
    type(NBodies_type_device), intent(inout) :: d_NB
    integer :: N, nd
    real(8) :: t_start, t_end
    N = h_NB%N_bodies
    nd = ndim
    t_start = omp_get_wtime()
    ! Update arrays needed for dump_data
    h_NB%pos(1:N, 1:nd) = d_NB%pos(1:N, 1:nd)
    h_NB%vel(1:N, 1:nd) = d_NB%vel(1:N, 1:nd)
    h_NB%accel(1:N, 1:nd) = d_NB%accel(1:N, 1:nd)
    h_NB%force(1:N, 1:3) = d_NB%force(1:N, 1:3)
    h_NB%kinetic_energy(1:N) = d_NB%kinetic_energy(1:N)
    h_NB%potential_energy(1:N) = d_NB%potential_energy(1:N)
    h_NB%sum_energy(1:N) = d_NB%sum_energy(1:N)
    t_end = omp_get_wtime()
    cuda_time_update_host = cuda_time_update_host + (t_end - t_start)
  end subroutine update_host_for_dump_cuda

  subroutine update_host_for_diagnostics_cuda(h_NB, d_NB)
    implicit none
    type(NBodies_type), intent(inout) :: h_NB
    type(NBodies_type_device), intent(inout) :: d_NB
    integer :: nd
    real(8) :: t_start, t_end

    nd = ndim
    t_start = omp_get_wtime()
    ! Update scalar values from module-level device variables
    h_NB%total_kinetic_energy = d_total_KE
    h_NB%total_potential_energy = d_total_PE
    h_NB%total_energy = h_NB%total_kinetic_energy + h_NB%total_potential_energy
    h_NB%center_of_mass(1:nd) = d_total_COM(1:nd)
    h_NB%center_of_mass_velocity(1:nd) = d_total_COM_vel(1:nd)
    h_NB%total_mass = d_total_mass
    h_NB%total_mass_inv = d_total_mass_inv
    t_end = omp_get_wtime()
    cuda_time_update_host = cuda_time_update_host + (t_end - t_start)
  end subroutine update_host_for_diagnostics_cuda

  subroutine cleanup_device_cuda(d_NB)
    ! delete all cuda arrays
    implicit none
    type(NBodies_type_device), intent(inout) :: d_NB
    integer :: istat
    ! ~~~~~~~~~~~~~~~~~~~~~ CLEANUP ARRAYS ~~~~~~~~~~~~~~~~~~~~~~
    if (allocated(d_NB%n_bodies)) deallocate (d_NB%n_bodies, stat=istat)
    if (allocated(d_NB%pos)) deallocate (d_NB%pos, stat=istat)
    if (allocated(d_NB%pos_0)) deallocate (d_NB%pos_0, stat=istat)
    if (allocated(d_NB%pos_d)) deallocate (d_NB%pos_d, stat=istat)
    if (allocated(d_NB%vel)) deallocate (d_NB%vel, stat=istat)
    if (allocated(d_NB%vel_0)) deallocate (d_NB%vel_0, stat=istat)
    if (allocated(d_NB%vel_h)) deallocate (d_NB%vel_h, stat=istat)
    if (allocated(d_NB%accel)) deallocate (d_NB%accel, stat=istat)
    if (allocated(d_NB%force)) deallocate (d_NB%force, stat=istat)
    if (allocated(d_NB%bc_factor)) deallocate (d_NB%bc_factor, stat=istat)
    if (allocated(d_NB%mass)) deallocate (d_NB%mass, stat=istat)
    if (allocated(d_NB%mass_inv)) deallocate (d_NB%mass_inv, stat=istat)
    if (allocated(d_NB%radius)) deallocate (d_NB%radius, stat=istat)
    if (allocated(d_NB%radius2)) deallocate (d_NB%radius2, stat=istat)
    if (allocated(d_NB%kinetic_energy)) deallocate (d_NB%kinetic_energy, stat=istat)
    if (allocated(d_NB%potential_energy)) deallocate (d_NB%potential_energy, stat=istat)
    if (allocated(d_NB%sum_energy)) deallocate (d_NB%sum_energy, stat=istat)
  end subroutine cleanup_device_cuda

  attributes(global) subroutine zero_forces_kernel(forces, potential_energy, Nbodies, ndim)
    implicit none
    integer, value :: Nbodies, ndim
    real(kind=wp), device :: forces(Nbodies, ndim)
    real(kind=wp), device :: potential_energy(Nbodies)
    integer :: i, nd
    i = blockDim%x*(blockIdx%x - 1) + threadIdx%x
    if (i <= Nbodies) then
      do nd = 1, ndim
        forces(i, nd) = 0.0_wp
      end do
      potential_energy(i) = 0.0_wp
    end if
  end subroutine zero_forces_kernel

  attributes(global) subroutine compute_forces_kernel(pos, mass, force, potential_energy, Nbodies, ndim, Gravitational_Constant, radius2)
    implicit none
    integer, value :: Nbodies, ndim
    real(kind=wp), value :: Gravitational_Constant
    real(kind=wp), device :: pos(Nbodies, ndim)
    real(kind=wp), device :: mass(Nbodies)
    real(kind=wp), device :: force(Nbodies, ndim)
    real(kind=wp), device :: potential_energy(Nbodies)
    real(kind=wp), device :: radius2(Nbodies)
    integer :: i_r, i_d, nd
    real(kind=wp) :: pos_diff_1, pos_diff_2, pos_diff_3
    real(kind=wp) :: dist, dist_squared, dist_squared_soft, dist_cubed, dist_cubed_inv
    real(kind=wp) :: force_local_1, force_local_2, force_local_3, pe_local
    real(kind=wp) :: mass_i, mass_product, f_mag

    ! Initialize local variables
    pos_diff_1 = 0.0_wp
    pos_diff_2 = 0.0_wp
    pos_diff_3 = 0.0_wp
    force_local_1 = 0.0_wp
    force_local_2 = 0.0_wp
    force_local_3 = 0.0_wp
    pe_local = 0.0_wp

i_r = blockDim%x*(blockIdx%x - 1) + threadIdx%x
    if (i_r <= Nbodies) then
      ! Cache frequently accessed values
      mass_i = mass(i_r)

      do i_d = 1, Nbodies
        if (i_r /= i_d) then
          ! Compute displacement vector
          pos_diff_1 = pos(i_d, 1) - pos(i_r, 1)
          if (ndim >= 2) then
            pos_diff_2 = pos(i_d, 2) - pos(i_r, 2)
          end if
          if (ndim == 3) then
            pos_diff_3 = pos(i_d, 3) - pos(i_r, 3)
          end if

          ! Compute squared distance with softening
          dist_squared = pos_diff_1**2 + pos_diff_2**2 + pos_diff_3**2
          dist_squared_soft = dist_squared + radius2(i_r) + radius2(i_d)
          dist = sqrt(dist_squared_soft)
          dist_cubed = dist_squared_soft*dist
          dist_cubed_inv = 1.0_wp/dist_cubed

          ! Compute force magnitude
          mass_product = mass_i*mass(i_d)
          f_mag = Gravitational_Constant*mass_product*dist_cubed_inv

          ! Update forces
          force_local_1 = force_local_1 + f_mag*pos_diff_1
          if (ndim >= 2) then
            force_local_2 = force_local_2 + f_mag*pos_diff_2
          end if
          if (ndim == 3) then
            force_local_3 = force_local_3 + f_mag*pos_diff_3
          end if

          ! Update potential energy
          pe_local = pe_local - 0.5_wp*Gravitational_Constant*mass_product/dist
        end if
      end do

      ! Write back computed forces and potential energy
      force(i_r, 1) = force_local_1
      if (ndim >= 2) then
        force(i_r, 2) = force_local_2
      end if
      if (ndim == 3) then
        force(i_r, 3) = force_local_3
      end if
      potential_energy(i_r) = pe_local
    end if

  end subroutine compute_forces_kernel

  attributes(global) subroutine calculate_particle_energy_kernel(kinetic_energy, potential_energy, sum_energy, vel, mass, Nbodies, ndim)
    implicit none
    integer, value :: Nbodies, ndim
    real(kind=wp), device :: kinetic_energy(Nbodies)
    real(kind=wp), device :: potential_energy(Nbodies)
    real(kind=wp), device :: sum_energy(Nbodies)
    real(kind=wp), device :: vel(Nbodies, ndim)
    real(kind=wp), device :: mass(Nbodies)
    integer :: i, nd
    real(kind=wp) :: v2
    i = blockDim%x*(blockIdx%x - 1) + threadIdx%x
    if (i <= Nbodies) then
      ! Calculate kinetic energy
      v2 = 0.0_wp
      do nd = 1, ndim
        v2 = v2 + vel(i, nd)**2
      end do
      kinetic_energy(i) = 0.5_wp*mass(i)*v2
      ! Total energy per particle
      sum_energy(i) = kinetic_energy(i) + potential_energy(i)
    end if
  end subroutine calculate_particle_energy_kernel

  attributes(global) subroutine sum_energy_kernel(COM, COM_vel, total_KE, total_PE, kinetic_energy, potential_energy, pos, vel, mass, total_mass_inv, Nbodies, ndim)
    ! sum up energies and compute COM
    implicit none
    integer, value :: Nbodies, ndim
    real(kind=wp), value :: total_mass_inv
    real(kind=wp), device :: COM(3)
    real(kind=wp), device :: COM_vel(3)
    real(kind=wp), device :: total_KE, total_PE
    real(kind=wp), device :: kinetic_energy(Nbodies)
    real(kind=wp), device :: potential_energy(Nbodies)
    real(kind=wp), device :: pos(Nbodies, ndim)
    real(kind=wp), device :: vel(Nbodies, ndim)
    real(kind=wp), device :: mass(Nbodies)
    integer :: i, nd
    real(kind=wp) :: mass_i, weighted_pos, weighted_vel
    real(kind=wp) :: atomic_tmp

    i = blockDim%x*(blockIdx%x - 1) + threadIdx%x
    if (i <= Nbodies) then
      ! Cache mass value
      mass_i = mass(i)

      ! Use atomic operations to safely accumulate totals
      atomic_tmp = atomicAdd(total_KE, kinetic_energy(i))
      atomic_tmp = atomicAdd(total_PE, potential_energy(i))

      do nd = 1, ndim
        weighted_pos = mass_i*pos(i, nd)*total_mass_inv
        weighted_vel = mass_i*vel(i, nd)*total_mass_inv
        atomic_tmp = atomicAdd(COM(nd), weighted_pos)
        atomic_tmp = atomicAdd(COM_vel(nd), weighted_vel)
      end do
    end if
  end subroutine sum_energy_kernel

  attributes(global) subroutine update_positions_kernel(Nbodies, ndim, force, accel, mass_inv, pos, vel, pos_0, vel_0, dt, dt_half)
    implicit none
    integer, value :: Nbodies, ndim
    real(kind=wp), device :: force(Nbodies, ndim)
    real(kind=wp), device :: accel(Nbodies, ndim)
    real(kind=wp), device :: mass_inv(Nbodies)
    real(kind=wp), device :: pos(Nbodies, ndim)
    real(kind=wp), device :: vel(Nbodies, ndim)
    real(kind=wp), device :: pos_0(Nbodies, ndim)
    real(kind=wp), device :: vel_0(Nbodies, ndim)
    real(kind=wp), value :: dt, dt_half
    integer :: i, nd
    real(kind=wp) :: vel_half(3)
    ! update position using velocity Verlet
    i = blockDim%x*(blockIdx%x - 1) + threadIdx%x
    if (i <= Nbodies) then
      do nd = 1, ndim
        ! Compute acceleration
        accel(i, nd) = force(i, nd)*mass_inv(i)
        ! Update velocity at half step
        vel_half(nd) = vel_0(i, nd) + accel(i, nd)*dt_half
        vel(i, nd) = vel_0(i, nd) + accel(i, nd)*dt
        ! Update position
        pos(i, nd) = pos_0(i, nd) + vel_half(nd)*dt
        ! store old position and velocity for next step
        pos_0(i, nd) = pos(i, nd)
        vel_0(i, nd) = vel(i, nd)
      end do
    end if
  end subroutine update_positions_kernel

  subroutine compute_forces_cuda(h_NB, d_NB)
    implicit none
    type(NBodies_type), intent(inout) :: h_NB
    type(NBodies_type_device), intent(inout) :: d_NB
    integer :: N, nd, num_blocks
    integer :: ierr
    N = h_NB%N_bodies
    nd = ndim
    ! Determine number of blocks
    num_blocks = (N + num_threads - 1)/num_threads
    t_cuda_start = omp_get_wtime()
    ! Zero forces and potential energy
    call zero_forces_kernel<<<num_blocks, num_threads>>>(d_NB%force, d_NB%potential_energy, N, nd)
    ierr = cudaDeviceSynchronize()
    ! Compute forces
    call compute_forces_kernel<<<num_blocks, num_threads>>>(d_NB%pos, d_NB%mass, d_NB%force, d_NB%potential_energy, N, nd, Gravitational_Constant, d_NB%radius2)
    ! Synchronize to ensure completion
    ierr = cudaDeviceSynchronize()
    t_cuda_end = omp_get_wtime()
    cuda_time_forces = cuda_time_forces + (t_cuda_end - t_cuda_start)
  end subroutine compute_forces_cuda

subroutine get_global_metrics_cuda(h_NB, d_NB)
    implicit none
    type(NBodies_type), intent(inout) :: h_NB
    type(NBodies_type_device), intent(inout) :: d_NB
    integer :: N, nd, num_blocks
    integer :: ierr
    real(kind=wp) :: total_mass_inv_value  ! Local host variable

    N = h_NB%N_bodies
    nd = ndim
    
    ! Get the value from host
    total_mass_inv_value = h_NB%total_mass_inv
    
    ! Initialize accumulators to zero
    d_total_KE = 0.0_wp
    d_total_PE = 0.0_wp
    d_total_COM = 0.0_wp
    d_total_COM_vel = 0.0_wp
    
    ! Determine number of blocks
    num_blocks = (N + num_threads - 1)/num_threads
    t_cuda_start = omp_get_wtime()
    ! First compute per-particle energies
    call calculate_particle_energy_kernel<<<num_blocks, num_threads>>>(d_NB%kinetic_energy, d_NB%potential_energy, d_NB%sum_energy, &
                                                                       d_NB%vel, d_NB%mass, N, nd)
    ! Synchronize to ensure completion
    ierr = cudaDeviceSynchronize()
    ! sum up energies and compute COM - pass total_mass_inv as a value
    call sum_energy_kernel<<<num_blocks, num_threads>>>(d_total_COM, d_total_COM_vel, &
                                                        d_total_KE, d_total_PE, &
                                                        d_NB%kinetic_energy, d_NB%potential_energy, d_NB%pos, d_NB%vel, d_NB%mass, &
                                                        total_mass_inv_value, N, nd)
    ! Synchronize to ensure completion
    ierr = cudaDeviceSynchronize()
    t_cuda_end = omp_get_wtime()
    cuda_time_metrics = cuda_time_metrics + (t_cuda_end - t_cuda_start)
  end subroutine get_global_metrics_cuda

subroutine update_positions_cuda(h_NB, d_NB)
    implicit none
    type(NBodies_type), intent(inout) :: h_NB
    type(NBodies_type_device), intent(inout) :: d_NB
    integer :: N, nd, num_blocks
    integer :: ierr
    real(8) :: t_start, t_end
    real(kind=wp) :: dt_value, dt_half_value  ! Local host variables
    
    N = h_NB%N_bodies
    nd = ndim
    
    ! Get values from host module variables
    dt_value = dt
    dt_half_value = dt_half
    
    ! Determine number of blocks
    num_blocks = (N + num_threads - 1)/num_threads
    t_start = omp_get_wtime()
    ! Update positions using velocity Verlet
    call update_positions_kernel<<<num_blocks, num_threads>>>(N, nd, &
                                                              d_NB%force, d_NB%accel, d_NB%mass_inv, &
                                                              d_NB%pos, d_NB%vel, d_NB%pos_0, d_NB%vel_0, &
                                                              dt_value, dt_half_value)
    ! Synchronize to ensure completion
    ierr = cudaDeviceSynchronize()
    t_end = omp_get_wtime()
    cuda_time_verlet = cuda_time_verlet + (t_end - t_start)
  end subroutine update_positions_cuda

  subroutine initialize_cuda_timers()
    implicit none
    cuda_time_forces = 0.0d0
    cuda_time_verlet = 0.0d0
    cuda_time_update_host = 0.0d0
    cuda_time_upload = 0.0d0
    cuda_time_metrics = 0.0d0
  end subroutine initialize_cuda_timers

  subroutine print_cuda_timers()
    implicit none
    real(8) :: total_cuda_time
    total_cuda_time = cuda_time_forces + cuda_time_verlet + cuda_time_update_host + cuda_time_upload + cuda_time_metrics
    
    write(*,'(A)') "CUDA Timing Summary:"
    write(*,'(A,F12.6,A)') "  Data Upload:        ", cuda_time_upload, " seconds"
    write(*,'(A,F12.6,A)') "  Force Calculation:  ", cuda_time_forces, " seconds"
    write(*,'(A,F12.6,A)') "  Position Update:    ", cuda_time_verlet, " seconds"
    write(*,'(A,F12.6,A)') "  Metrics Calc:       ", cuda_time_metrics, " seconds"
    write(*,'(A,F12.6,A)') "  Host Update:        ", cuda_time_update_host, " seconds"
    write(*,'(A,F12.6,A)') "  Total CUDA Time:    ", total_cuda_time, " seconds"
  end subroutine print_cuda_timers

end module cuda_nbodies