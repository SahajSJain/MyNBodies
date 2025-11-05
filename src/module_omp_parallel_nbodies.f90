module omp_nbodies
  !!! This specific OMP module was entirely LLM generated to provide OpenMP parallel implementations 
  !!! For my personal implementation, look at module_myomp_nbodies.f90 
  !!! This was generated for comparison and educational purposes.
  use omp_lib
  use nbodies
  implicit none
contains
  subroutine compute_forces_omp(NB)
    implicit none
    type(NBodies_type), intent(inout) :: NB
    integer :: i, j
    real(kind=wp) :: pos_diff(ndim)
    real(kind=wp) :: dist, dist_cubed ! distance and distance cubed between bodies
    real(kind=wp) :: dist_squared, dist_squared_soft ! distance squared and softened distance squared
    real(kind=wp) :: force_local(ndim)
    real(kind=wp) :: pe_local

    ! Initialize forces to zero
    !$omp parallel do default(shared) private(i)
    do i = 1, NB%n_bodies
      NB%force(i, :) = 0.0_wp
      NB%potential_energy(i) = 0.0_wp
      NB%kinetic_energy(i) = 0.0_wp
    end do
    !$omp end parallel do

    NB%total_potential_energy = 0.0_wp
    NB%total_kinetic_energy = 0.0_wp
    NB%total_energy = 0.0_wp

    ! Compute pairwise gravitational forces
    !$omp parallel default(shared) &
    !$omp private(i, j, pos_diff, dist_squared, dist_squared_soft, dist, dist_cubed, force_local, pe_local)
    !$omp do
    do i = 1, NB%n_bodies
      force_local = 0.0_wp
      pe_local = 0.0_wp

      do j = 1, NB%n_bodies
        if (i /= j) then
          ! Calculate position difference
          pos_diff = NB%pos(j, 1:ndim) - NB%pos(i, 1:ndim)

          ! Calculate distance squared first (more efficient)
          dist_squared = sum(pos_diff(1:ndim)**2)

          dist_squared_soft = dist_squared + softening_length_squared
          dist = sqrt(dist_squared_soft)
          dist_cubed = dist_squared_soft * dist  ! More numerically stable

          ! Calculate gravitational force
          ! F = G * m1 * m2 * r_vec / r^3
          force_local = force_local + &
            Gravitational_Constant * NB%mass(i) * NB%mass(j) * pos_diff(1:ndim) / dist_cubed

          ! Calculate potential energy
          ! Half to avoid double counting
          pe_local = pe_local - 0.5_wp * Gravitational_Constant * NB%mass(i) * NB%mass(j) / dist
        end if
      end do

      !$omp critical
      NB%force(i, 1:ndim) = NB%force(i, 1:ndim) + force_local
      NB%potential_energy(i) = NB%potential_energy(i) + pe_local
      !$omp end critical
    end do
    !$omp end do
    !$omp end parallel

    ! Calculate total potential energy
    NB%total_potential_energy = sum(NB%potential_energy)

  end subroutine compute_forces_omp

  subroutine get_global_metrics_omp(NB)
    implicit none
    type(NBodies_type), intent(inout) :: NB
    integer :: i
    real(kind=wp) :: ke_temp
    real(kind=wp) :: com(ndim), comv(ndim)

    NB%total_kinetic_energy = 0.0_wp
    NB%total_potential_energy = 0.0_wp
    NB%total_energy = 0.0_wp
    NB%center_of_mass = 0.0_wp
    NB%center_of_mass_velocity = 0.0_wp

    !$omp parallel do default(shared) private(i, ke_temp)
    do i = 1, NB%n_bodies
      ke_temp = 0.5_wp * NB%mass(i) * sum(NB%vel(i, 1:ndim)**2)
      NB%kinetic_energy(i) = ke_temp
      NB%sum_energy(i) = ke_temp + NB%potential_energy(i)
    end do
    !$omp end parallel do

    total_ke = 0.0_wp
    total_pe = 0.0_wp
    com = 0.0_wp
    comv = 0.0_wp

    !$omp parallel default(shared) private(i) &
    !$omp reduction(+:total_ke, total_pe, com, comv)
    !$omp do
    do i = 1, NB%n_bodies
      total_ke = total_ke + NB%kinetic_energy(i)
      ! no need to multiply by 0.5 here since already done in compute_forces_omp
      total_pe = total_pe + NB%potential_energy(i)

      com = com + NB%mass(i) * NB%pos(i, 1:ndim) * NB%total_mass_inv
      comv = comv + NB%mass(i) * NB%vel(i, 1:ndim) * NB%total_mass_inv
    end do
    !$omp end do
    !$omp end parallel

    NB%total_kinetic_energy = total_ke
    NB%total_potential_energy = total_pe
    NB%center_of_mass = com
    NB%center_of_mass_velocity = comv
    NB%total_energy = NB%total_kinetic_energy + NB%total_potential_energy

  end subroutine get_global_metrics_omp

  subroutine update_positions_omp(NB)
    implicit none
    type(NBodies_type), intent(inout) :: NB
    integer :: i, nd
    real(kind=wp) :: ke_temp

    ! Reset energies
    !$omp parallel do default(shared) private(i)
    do i = 1, NB%n_bodies
      NB%kinetic_energy(i) = 0.0_wp
      NB%sum_energy(i) = 0.0_wp
    end do
    !$omp end parallel do

    ! Calculate acceleration from forces
    !$omp parallel do default(shared) private(i)
    do i = 1, NB%n_bodies
      NB%accel(i, 1:ndim) = NB%force(i, 1:ndim) * NB%mass_inv(i)
    end do
    !$omp end parallel do

    ! Update velocities (full step for diagnostics, half step for position update)
    !$omp parallel do default(shared) private(i)
    do i = 1, NB%n_bodies
      ! Half-step velocity for position update
      NB%vel_h(i, 1:ndim) = NB%vel_0(i, 1:ndim) + NB%accel(i, 1:ndim) * dt_half
      ! Full-step velocity for energy calculation
      NB%vel(i, 1:ndim) = NB%vel_0(i, 1:ndim) + NB%accel(i, 1:ndim) * dt
    end do
    !$omp end parallel do

    if (REFLECTIVE_BC) then
      !$omp parallel do default(shared) private(i, nd) collapse(2)
      do i = 1, NB%n_bodies
        do nd = 1, ndim
          NB%bc_factor(i, nd) = (1.0_wp &
            - iReflective_BC(nd) * 2.0_wp * merge(1.0_wp, 0.0_wp, &
            ((NB%pos(i, nd) > +L_bound(nd) * 0.50_wp .and. NB%vel(i, nd) > 0.0_wp) &
            .or. &
            (NB%pos(i, nd) < -L_bound(nd) * 0.50_wp .and. NB%vel(i, nd) < 0.0_wp))))
          NB%vel_h(i, nd) = NB%vel_h(i, nd) * NB%bc_factor(i, nd)
          NB%vel(i, nd) = NB%vel(i, nd) * NB%bc_factor(i, nd)
        end do
      end do
      !$omp end parallel do
    end if

    ! Update positions using half-step velocities
    !$omp parallel do default(shared) private(i)
    do i = 1, NB%n_bodies
      NB%pos(i, 1:ndim) = NB%pos_0(i, 1:ndim) + NB%vel_h(i, 1:ndim) * dt
    end do
    !$omp end parallel do

    ! Calculate kinetic energy using full-step velocities
    !$omp parallel do default(shared) private(i, ke_temp)
    do i = 1, NB%n_bodies
      ke_temp = 0.5_wp * NB%mass(i) * sum(NB%vel(i, 1:ndim)**2)
      NB%kinetic_energy(i) = ke_temp
      NB%sum_energy(i) = ke_temp + NB%potential_energy(i)
    end do
    !$omp end parallel do

    ! Store current values as old values for next iteration
    !$omp parallel do default(shared) private(i)
    do i = 1, NB%n_bodies
      NB%pos_0(i, 1:ndim) = NB%pos(i, 1:ndim)
      NB%vel_0(i, 1:ndim) = NB%vel(i, 1:ndim)
    end do
    !$omp end parallel do
  end subroutine update_positions_omp
  subroutine softening_length_setup_omp(NB)
    implicit none
    type(NBodies_type), intent(inout) :: NB
    real(kind=wp) :: E_total, KE_max, KE_initial, PE_initial
    real(kind=wp) :: v_max, r_min
! BASIC IDEA FOR SOFTENING LENGTH ESTIMATION:
!     # Estimate based on energy
!     Estimate initial total kinetic energy (KE) and potential energy (PE)
!     and E_total = KE_initial + PE_initial
!       Then say KE_Max :
!    if E_total < 0:  # Bound system
!        KE_max = abs(E_total) + KE_initial
!    else:  # Unbound system
!        KE_max = KE_initial + abs(PE_initial)
!     if M=total mass
!     v²_max ≈ 2 * KE_max / M       # Maximum velocity
! # Minimum approach distance
! r_min ≈ G * m / v²_max = G * m² * R / (2 * M_total²)
! # Softening
! ε = 0.1 * r_min
    KE_initial = NB%total_kinetic_energy
    PE_initial = NB%total_potential_energy
    E_total = KE_initial + PE_initial
    if (E_total < 0.0_wp) then
      KE_max = abs(E_total) + KE_initial
    else
      KE_max = KE_initial + abs(PE_initial)
    end if
    ! Estimate maximum velocity
    v_max = sqrt(2.0_wp * KE_max * NB%total_mass_inv) ! avoid sqrt(0)
    r_min = Gravitational_Constant * maxval(NB%mass(1:NB%n_bodies)) / (v_max**2 + small_number)
    softening_length = 100.0_wp * r_min
    ! Set softening length squared based on radius_0
    softening_length_squared = softening_length**2
    print *, "Estimated softening length: ", softening_length
  end subroutine softening_length_setup_omp
end module omp_nbodies
