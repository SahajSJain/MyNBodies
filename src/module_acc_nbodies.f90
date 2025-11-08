module acc_nbodies
  use nbodies
  use omp_lib
  implicit none

  real(8) :: total_time_forces, total_time_verlet, t_start, t_end, total_time_update_host, total_time_upload, total_time_metrics

contains
  subroutine move_data_to_device_acc(NB)
    implicit none
    type(NBodies_type), intent(inout) :: NB
    integer :: N, nd
    N = NB%N_bodies
    nd = ndim
    t_start = omp_get_wtime()
    ! ~~~~~~~~~~~~~~~~~~~~~ MOVE ENTIRE STRUCTURE TO DEVICE ~~~~~~~~~~~~~~~~~~~~~~
    ! First, copy the entire NB structure to the device
    !$acc enter data copyin(NB) async(1)

    ! Now copy all the allocatable arrays within the structure
    ! Arrays that are initialized on host
    !$acc enter data copyin(NB%pos(1:N,1:nd)) async(1)
    !$acc enter data copyin(NB%pos_0(1:N,1:nd)) async(1)
    !$acc enter data copyin(NB%vel(1:N,1:nd)) async(1)
    !$acc enter data copyin(NB%vel_0(1:N,1:nd)) async(1)
    !$acc enter data copyin(NB%mass(1:N)) async(1)
    !$acc enter data copyin(NB%mass_inv(1:N)) async(1)
    !$acc enter data copyin(NB%radius(1:N)) async(1)
    !$acc enter data copyin(NB%radius2(1:N)) async(1) ! new array for softening
    ! Arrays that will be computed on device - create them
    !$acc enter data create(NB%vel_h(1:N,1:nd)) async(1)
    !$acc enter data create(NB%accel(1:N,1:nd)) async(1)
    !$acc enter data create(NB%force(1:N,1:3)) async(1)
    !$acc enter data create(NB%bc_factor(1:N,1:nd)) async(1)
    !$acc enter data create(NB%kinetic_energy(1:N)) async(1)
    !$acc enter data create(NB%potential_energy(1:N)) async(1)
    !$acc enter data create(NB%sum_energy(1:N)) async(1)
    !$acc enter data create(NB%center_of_mass(1:3)) async(1)
    !$acc enter data create(NB%center_of_mass_velocity(1:3)) async(1)

    ! ~~~~~~~~~~~~~~~~~~~~~ OTHER GLOBAL VARIABLES ~~~~~~~~~~~~~~~~~~~~~~
    ! Copy other module variables needed on device
    !$acc enter data copyin(NB%N_bodies) async(1)
    !$acc enter data copyin(NB%total_mass) async(1)
    !$acc enter data copyin(NB%total_mass_inv) async(1)
    !$acc enter data copyin(ndim) async(1)
    !$acc enter data copyin(Gravitational_Constant) async(1)
    ! small_number is a parameter - don't copy it
    !$acc enter data copyin(dt, dt_half) async(1)
    !$acc enter data copyin(REFLECTIVE_BC) async(1)
    !$acc enter data copyin(iReflective_BC(1:3)) async(1)
    !$acc enter data copyin(L_bound(1:3)) async(1)
    ! Create temporary variables for reductions
    !$acc enter data create(total_KE, total_PE)   async(1)
    !$acc enter data create(total_COM(1:3), total_COM_vel(1:3))   async(1)
    !$acc wait(1)
    t_end = omp_get_wtime()
    total_time_upload = total_time_upload + (t_end - t_start)
  end subroutine move_data_to_device_acc

  subroutine update_host_for_dump_acc(NB)
    implicit none
    type(NBodies_type), intent(inout) :: NB
    integer :: N, nd
    N = NB%N_bodies
    nd = ndim
    t_start = omp_get_wtime()
    ! Update arrays needed for dump_data
    !$acc update host(NB%pos(1:N,1:nd)) async(1)
    !$acc update host(NB%vel(1:N,1:nd)) async(1)
    !$acc update host(NB%accel(1:N,1:nd)) async(1)
    !$acc update host(NB%force(1:N,1:3)) async(1)
    !$acc update host(NB%kinetic_energy(1:N)) async(1)
    !$acc update host(NB%potential_energy(1:N)) async(1)
    !$acc update host(NB%sum_energy(1:N)) async(1)
    !$acc wait(1)
    t_end = omp_get_wtime()
    total_time_update_host = total_time_update_host + (t_end - t_start)
  end subroutine update_host_for_dump_acc

  subroutine update_host_for_diagnostics_acc(NB)
    implicit none
    type(NBodies_type), intent(inout) :: NB
    integer :: nd
    nd = ndim
    t_start = omp_get_wtime()
    ! Update scalar values in the structure
    !$acc update host(NB%total_kinetic_energy) async(1)
    !$acc update host(NB%total_potential_energy) async(1)
    !$acc update host(NB%total_energy) async(1)
    !$acc update host(NB%center_of_mass(1:3)) async(1)
    !$acc update host(NB%center_of_mass_velocity(1:3)) async(1)
    !$acc update host(NB%total_energy_initial) async(1)
    !$acc wait(1)
    t_end = omp_get_wtime()
    total_time_update_host = total_time_update_host + (t_end - t_start)
  end subroutine update_host_for_diagnostics_acc

  subroutine cleanup_device_acc(NB)
    implicit none
    type(NBodies_type), intent(inout) :: NB
    integer :: N, nd
    N = NB%N_bodies
    nd = ndim; 
    ! ~~~~~~~~~~~~~~~~~~~~~ CLEANUP ARRAYS ~~~~~~~~~~~~~~~~~~~~~~
    ! Delete allocatable arrays first
    !$acc exit data delete(NB%pos(1:N,1:nd))
    !$acc exit data delete(NB%pos_0(1:N,1:nd))
    !$acc exit data delete(NB%vel(1:N,1:nd))
    !$acc exit data delete(NB%vel_0(1:N,1:nd))
    !$acc exit data delete(NB%vel_h(1:N,1:nd))
    !$acc exit data delete(NB%accel(1:N,1:nd))
    !$acc exit data delete(NB%force(1:N,1:3))
    !$acc exit data delete(NB%bc_factor(1:N,1:nd))
    !$acc exit data delete(NB%mass(1:N))
    !$acc exit data delete(NB%mass_inv(1:N))
    !$acc exit data delete(NB%radius(1:N))
    !$acc exit data delete(NB%radius2(1:N)) ! new array for softening
    !$acc exit data delete(NB%kinetic_energy(1:N))
    !$acc exit data delete(NB%potential_energy(1:N))
    !$acc exit data delete(NB%sum_energy(1:N))
    !$acc exit data delete(NB%center_of_mass(1:nd))
    !$acc exit data delete(NB%center_of_mass_velocity(1:nd))

    ! ~~~~~~~~~~~~~~~~~~~~~ CLEANUP STRUCTURE ~~~~~~~~~~~~~~~~~~~~~~
    ! Delete the structure itself
    !$acc exit data delete(NB)

    ! ~~~~~~~~~~~~~~~~~~~~~ CLEANUP OTHER VARIABLES ~~~~~~~~~~~~~~~~~~~~~~
    !$acc exit data delete(ndim)
    !$acc exit data delete(Gravitational_Constant)
    !$acc exit data delete(dt, dt_half)
    !$acc exit data delete(REFLECTIVE_BC)
    !$acc exit data delete(iReflective_BC(1:3))
    !$acc exit data delete(L_bound(1:3))
    !$acc exit data delete(total_KE, total_PE)
    !$acc exit data delete(total_COM(1:3), total_COM_vel(1:3))
    !$acc exit data delete(NB%N_bodies)
    !$acc exit data delete(NB%total_mass)
    !$acc exit data delete(NB%total_mass_inv)
    !$acc exit data delete(NB%total_energy_initial)
  end subroutine cleanup_device_acc

  subroutine compute_forces_acc(NB)
    implicit none
    type(NBodies_type), intent(inout) :: NB
    integer :: i_r, i_d, i, nd, k
    real(kind=wp) :: pos_diff_1, pos_diff_2, pos_diff_3
    real(kind=wp) :: force_local_1, force_local_2, force_local_3
    real(kind=wp) :: dist, dist_cubed
    real(kind=wp) :: dist_squared, dist_squared_soft
    real(kind=wp) :: prefactor
    real(kind=wp) :: pe_local

    t_start = omp_get_wtime()

    ! Initialize forces to zero
    !$acc kernels present(NB, NB%total_potential_energy, NB%total_kinetic_energy, NB%total_energy)
    NB%total_potential_energy = 0.0_wp
    NB%total_kinetic_energy = 0.0_wp
    NB%total_energy = 0.0_wp
    !$acc end kernels

    !$acc parallel loop gang vector private(i,nd) &
    !$acc& present(NB, NB%force, NB%potential_energy, NB%N_bodies, ndim, NB%kinetic_energy) default(none)
    do i = 1, NB%N_bodies
      !$acc loop seq
      do nd = 1, ndim
        NB%force(i, nd) = 0.0_wp
      end do
      NB%potential_energy(i) = 0.0_wp
      NB%kinetic_energy(i) = 0.0_wp
    end do
    !$acc end parallel loop

!$acc parallel loop gang vector present(NB, NB%pos, NB%mass, NB%force, NB%potential_energy, &
!$acc& NB%N_bodies, ndim, Gravitational_Constant, NB%radius2) &
!$acc& private(i_r,i_d,k,dist_squared,dist_squared_soft,dist,dist_cubed,pe_local, pos_diff_1, pos_diff_2, pos_diff_3, force_local_1, force_local_2, force_local_3) default(none)
    do i_r = 1, NB%n_bodies
      ! Initialize with scalar operations
      force_local_1 = 0.0_wp
      force_local_2 = 0.0_wp
      force_local_3 = 0.0_wp
      pos_diff_1 = 0.0_wp
      pos_diff_2 = 0.0_wp
      pos_diff_3 = 0.0_wp
      pe_local = 0.0_wp

      !$acc loop seq
      do i_d = 1, NB%n_bodies
        if (i_r /= i_d) then
          ! Calculate everything with explicit scalar operations
          pos_diff_1 = NB%pos(i_d, 1) - NB%pos(i_r, 1)
          if (ndim >= 2) then
            pos_diff_2 = NB%pos(i_d, 2) - NB%pos(i_r, 2)
          end if
          if (ndim == 3) then
            pos_diff_3 = NB%pos(i_d, 3) - NB%pos(i_r, 3)
          end if
          dist_squared = pos_diff_1**2 + pos_diff_2**2 + pos_diff_3**2

          dist_squared_soft = dist_squared + NB%radius2(i_r) + NB%radius2(i_d)
          dist = sqrt(dist_squared_soft)
          dist_cubed = dist_squared_soft*dist

          force_local_1 = force_local_1 + &
                          Gravitational_Constant*NB%mass(i_r)*NB%mass(i_d)*pos_diff_1/dist_cubed
          if (ndim >= 2) then
            force_local_2 = force_local_2 + &
                            Gravitational_Constant*NB%mass(i_r)*NB%mass(i_d)*pos_diff_2/dist_cubed
          end if
          if (ndim == 3) then
            force_local_3 = force_local_3 + &
                            Gravitational_Constant*NB%mass(i_r)*NB%mass(i_d)*pos_diff_3/dist_cubed
          end if
          ! Update potential energy with scalar operations
          pe_local = pe_local - &
                     0.50_wp*Gravitational_Constant*NB%mass(i_r)*NB%mass(i_d)/dist
        end if
      end do
      ! Store results with scalar operations
      NB%force(i_r, 1) = force_local_1
      if (ndim >= 2) then
        NB%force(i_r, 2) = force_local_2
      end if
      if (ndim == 3) then
        NB%force(i_r, 3) = force_local_3
      end if
      NB%potential_energy(i_r) = pe_local
    end do
!$acc end parallel loop
    t_end = omp_get_wtime()
    total_time_forces = total_time_forces + (t_end - t_start)
  end subroutine compute_forces_acc

  subroutine get_global_metrics_acc(NB)
    implicit none
    type(NBodies_type), intent(inout) :: NB
    integer :: i, nd
    t_start = omp_get_wtime()
    ! Initialize totals
    !$acc kernels &
    !$acc& present(NB, NB%total_potential_energy, NB%total_kinetic_energy, NB%total_energy, NB%center_of_mass, NB%center_of_mass_velocity, ndim,total_KE, total_PE, total_COM, total_COM_vel)
    NB%total_kinetic_energy = 0.0_wp
    NB%total_potential_energy = 0.0_wp
    NB%total_energy = 0.0_wp
    NB%center_of_mass(1:ndim) = 0.0_wp
    NB%center_of_mass_velocity(1:ndim) = 0.0_wp
    total_KE = 0.0_wp
    total_PE = 0.0_wp
    total_COM(1:3) = 0.0_wp
    total_COM_vel(1:3) = 0.0_wp
    !$acc end kernels

    ! First calculate kinetic energy and sum_energy per body
!$acc parallel present(NB, NB%kinetic_energy, NB%potential_energy, NB%sum_energy, NB%vel, NB%mass, &
    !$acc&                 NB%N_bodies, ndim, NB%total_kinetic_energy, NB%total_potential_energy, &
    !$acc&                 NB%center_of_mass, NB%pos, NB%total_mass_inv, NB%center_of_mass_velocity,total_com_vel, total_com, total_PE, total_KE) &
    !$acc& default(none)
    !$acc loop gang vector private(i)
    do i = 1, NB%N_bodies
      NB%kinetic_energy(i) = 0.5_wp*NB%mass(i)*sum(NB%vel(i, 1:ndim)**2)
      NB%sum_energy(i) = NB%kinetic_energy(i) + NB%potential_energy(i)
    end do
    !$acc end loop

    ! Now sum total energies and compute center of mass and center of mass velocity
    !$acc loop gang vector private(i) reduction(+:total_KE)
    do i = 1, NB%N_bodies
      total_KE = total_KE + NB%kinetic_energy(i)
    end do
    !$acc end loop

    !$acc loop gang vector private(i) reduction(+:total_PE)
    do i = 1, NB%N_bodies
      total_PE = total_PE + NB%potential_energy(i)
    end do
    !$acc end loop

    !$acc loop gang vector private(i, nd) reduction(+:total_COM(1:ndim))
    do i = 1, NB%N_bodies
      !$acc loop seq
      do nd = 1, ndim
        total_COM(nd) = total_COM(nd) + &
                        NB%mass(i)*NB%pos(i, nd)*NB%total_mass_inv
      end do
    end do
    !$acc end loop

    !$acc loop gang vector private(i, nd) reduction(+:total_COM_vel(1:ndim))
    do i = 1, NB%N_bodies
      !$acc loop seq
      do nd = 1, ndim
        total_COM_vel(nd) = total_COM_vel(nd) + &
                            NB%mass(i)*NB%vel(i, nd)*NB%total_mass_inv
      end do
    end do
    !$acc end loop
    !$acc end parallel

    !$acc kernels &
    !$acc& present(NB, NB%total_energy, NB%total_kinetic_energy, NB%total_potential_energy, total_KE, total_PE, total_COM, total_COM_vel, ndim, NB%center_of_mass, NB%center_of_mass_velocity)
    NB%total_kinetic_energy = total_KE
    NB%total_potential_energy = total_PE
    NB%total_energy = NB%total_kinetic_energy + NB%total_potential_energy
    NB%center_of_mass(1:ndim) = total_COM(1:ndim)
    NB%center_of_mass_velocity(1:ndim) = total_COM_vel(1:ndim)
    !$acc end kernels

    t_end = omp_get_wtime()
    total_time_metrics = total_time_metrics + (t_end - t_start)
  end subroutine get_global_metrics_acc

  subroutine update_positions_acc(NB)
    implicit none
    type(NBodies_type), intent(inout) :: NB
    integer :: i, nd
    ! Use Verlet integration to update positions and velocities
    ! First update accelerations based on forces
    ! Update velocities and positions based on computed forces
    ! in kick-drift-kick scheme form
    ! First update accelerations:
    ! Note: We have dimensional arrays:
    !     force(n_bodies, ndim)
    !     accel(n_bodies, ndim)
    !     vel(n_bodies, ndim)
    !     vel_0(n_bodies, ndim)
    !     vel_h(n_bodies, ndim)
    !     pos(n_bodies, ndim)
    !     pos_0(n_bodies, ndim)
    !   and non-dimensional arrays:
    !     mass(n_bodies)
    !     mass_inv(n_bodies)
    ! Create OMP PARALLEL region first to avoid overhead in loops
    ! Calculate acceleration from forces
    t_start = omp_get_wtime()
!$acc parallel present(NB, NB%force, NB%accel, NB%mass_inv, NB%N_bodies, ndim, &
!$acc&                 NB%vel_h, NB%vel, NB%vel_0, dt_half, dt, &
!$acc&                 NB%bc_factor, REFLECTIVE_BC, iReflective_BC, &
!$acc&                 NB%pos, L_bound, NB%pos_0) default(none)

!$acc loop gang vector collapse(2) private(i,nd)
    do i = 1, NB%N_bodies
      do nd = 1, ndim
        NB%accel(i, nd) = NB%force(i, nd)*NB%mass_inv(i); 
        ! Update velocities (full step for diagnostics, half step for position update)
        ! Half-step velocity for position update
        NB%vel_h(i, nd) = NB%vel_0(i, nd) + NB%accel(i, nd)*dt_half; 
        ! Full-step velocity for next iteration and diagnostics
        NB%vel(i, nd) = NB%vel_0(i, nd) + NB%accel(i, nd)*dt
        NB%pos(i, nd) = NB%pos_0(i, nd) + NB%vel_h(i, nd)*dt
        ! Store current values as old values for next iteration
        NB%pos_0(i, nd) = NB%pos(i, nd)
        NB%vel_0(i, nd) = NB%vel(i, nd)
      end do
    end do
!$acc end loop
!$acc end parallel
    t_end = omp_get_wtime()
    total_time_verlet = total_time_verlet + (t_end - t_start)
  end subroutine update_positions_acc

end module acc_nbodies

