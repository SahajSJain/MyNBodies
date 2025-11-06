module llmacc_nbodies
  use nbodies
  implicit none
contains
  subroutine compute_forces_llmacc(NB)
    implicit none
    type(NBodies_type), intent(inout) :: NB
    integer :: i_r, i_d
    real(kind=wp) :: pos_diff(ndim)
    real(kind=wp) :: dist, dist_cubed ! distance and distance cubed between bodies
    real(kind=wp) :: dist_squared, dist_squared_soft ! distance squared and softened distance squared
    real(kind=wp) :: force_local(ndim), pe_local
    
    ! Initialize forces to zero
    NB%force(:, :) = 0.0_wp
    NB%potential_energy = 0.0_wp ! Reset potential energy array
    NB%total_potential_energy = 0.0_wp ! Reset total potential energy
    NB%kinetic_energy = 0.0_wp ! Reset kinetic energy array
    NB%total_kinetic_energy = 0.0_wp ! Reset total kinetic energy
    NB%total_energy = 0.0_wp ! Reset total energy
    
    ! Compute pairwise gravitational forces
    ! Copy data to device and create data environment
    !$acc data copyin(NB%pos, NB%mass, NB%n_bodies) &
    !$acc      copyout(NB%force, NB%potential_energy) &
    !$acc      copyin(Gravitational_Constant, softening_length_squared)
    
    !$acc parallel loop gang vector private(i_r,i_d,pos_diff,dist_squared,dist_squared_soft,dist,dist_cubed,force_local,pe_local)
    do i_r = 1, NB%n_bodies ! i_r is the receiving body
      force_local = 0.0_wp
      pe_local = 0.0_wp
      
      !$acc loop seq
      do i_d = 1, NB%n_bodies ! i_d is the exerting body
        if (i_r /= i_d) then
          ! Calculate position difference
          pos_diff = NB%pos(i_d, 1:ndim) - NB%pos(i_r, 1:ndim)

          ! Calculate distance squared first (more efficient)
          dist_squared = sum(pos_diff(1:ndim)**2)

          dist_squared_soft = dist_squared + softening_length_squared
          dist = sqrt(dist_squared_soft)
          dist_cubed = dist_squared_soft * dist  ! More numerically stable

          ! Calculate gravitational force
          ! F = G * m1 * m2 * r_vec / r^3
          force_local(1:ndim) = force_local(1:ndim) + &
            Gravitational_Constant * NB%mass(i_r) * NB%mass(i_d) * pos_diff(1:ndim) / dist_cubed

          ! Calculate potential energy: half to avoid double counting
          pe_local = pe_local - &
            0.50_wp*Gravitational_Constant * NB%mass(i_r) * NB%mass(i_d) / dist ! Potential energy
          ! 0.50_wp to avoid double counting
        end if
      end do
      
      NB%force(i_r, 1:ndim) = force_local(1:ndim)
      NB%potential_energy(i_r) = pe_local
    end do
    !$acc end parallel loop
    
    !$acc end data

    ! Calculate total potential energy on host
    NB%total_potential_energy = sum(NB%potential_energy)

  end subroutine compute_forces_llmacc

  subroutine get_global_metrics_llmacc(NB)
    implicit none
    type(NBodies_type), intent(inout) :: NB
    integer :: i
    real(kind=wp) :: total_kinetic_energy, total_potential_energy, com_local(ndim), com_vel_local(ndim)
    
    NB%total_kinetic_energy = 0.0_wp
    NB%total_potential_energy = 0.0_wp
    NB%total_energy = 0.0_wp
    NB%center_of_mass = 0.0_wp
    NB%center_of_mass_velocity = 0.0_wp
    
    total_kinetic_energy = 0.0_wp
    total_potential_energy = 0.0_wp
    com_local = 0.0_wp
    com_vel_local = 0.0_wp
    
    ! Create ACC data environment and compute kinetic energies
    !$acc data copyin(NB%mass, NB%vel, NB%pos, NB%potential_energy, NB%n_bodies) &
    !$acc      copyout(NB%kinetic_energy, NB%sum_energy) &
    !$acc      copy(total_kinetic_energy, total_potential_energy, com_local, com_vel_local)
    
    !$acc parallel loop gang vector reduction(+:total_kinetic_energy, total_potential_energy, com_local, com_vel_local)
    do i = 1, NB%n_bodies
      NB%kinetic_energy(i) = 0.5_wp * NB%mass(i) * sum(NB%vel(i, 1:ndim)**2)
      NB%sum_energy(i) = NB%kinetic_energy(i) + NB%potential_energy(i)
      total_kinetic_energy = total_kinetic_energy + NB%kinetic_energy(i)
      total_potential_energy = total_potential_energy + NB%potential_energy(i)
      com_local(1:ndim) = com_local(1:ndim) + &
        NB%mass(i) * NB%pos(i, 1:ndim)
      com_vel_local(1:ndim) = com_vel_local(1:ndim) + &
        NB%mass(i) * NB%vel(i, 1:ndim)
    end do
    !$acc end parallel loop
    
    !$acc end data
    
    NB%total_kinetic_energy = total_kinetic_energy
    NB%total_potential_energy = total_potential_energy
    NB%center_of_mass(1:ndim) = com_local(1:ndim)/NB%total_mass
    NB%center_of_mass_velocity(1:ndim) = com_vel_local(1:ndim)/NB%total_mass
    NB%total_energy = NB%total_kinetic_energy + NB%total_potential_energy
  end subroutine get_global_metrics_llmacc

  subroutine update_positions_llmacc(NB)
    implicit none
    type(NBodies_type), intent(inout) :: NB
    integer :: i, nd
    
    ! Reset energies
    NB%kinetic_energy = 0.0_wp
    NB%sum_energy = 0.0_wp
    
    ! Create ACC data environment
    !$acc data copyin(NB%force, NB%mass_inv, NB%vel_0, NB%pos_0, NB%potential_energy, &
    !$acc            NB%mass, NB%n_bodies, dt, dt_half, REFLECTIVE_BC, iReflective_BC, L_bound) &
    !$acc      copyout(NB%accel, NB%vel, NB%vel_h, NB%pos, NB%bc_factor, &
    !$acc              NB%kinetic_energy, NB%sum_energy) &
    !$acc      copy(NB%pos_0, NB%vel_0)
    
    !$acc parallel loop gang vector private(i,nd)
    do i = 1, NB%n_bodies
      NB%accel(i, 1:ndim) = NB%force(i, 1:ndim) * NB%mass_inv(i)
      
      ! Update velocities (full step for diagnostics, half step for position update)
      ! Half-step velocity for position update
      NB%vel_h(i, 1:ndim) = NB%vel_0(i, 1:ndim) + NB%accel(i, 1:ndim) * dt_half
      ! Full-step velocity for energy calculation
      NB%vel(i, 1:ndim) = NB%vel_0(i, 1:ndim) + NB%accel(i, 1:ndim) * dt
      
      if (REFLECTIVE_BC) then
        !$acc loop seq
        do nd = 1, ndim
          NB%bc_factor(i, nd) = (1.0_wp &
            - iReflective_BC(nd) * 2.0_wp * merge(1.0_wp, 0.0_wp, &
            ((NB%pos(i, nd) > +L_bound(nd) * 0.50_wp .and. NB%vel(i, nd) > &
            0.0_wp) &
            .or. &
            (NB%pos(i, nd) < -L_bound(nd) * 0.50_wp .and. NB%vel(i, nd) < &
            0.0_wp))))
          NB%vel_h(i, nd) = NB%vel_h(i, nd) * NB%bc_factor(i, nd)
          NB%vel(i, nd) = NB%vel(i, nd) * NB%bc_factor(i, nd)
        end do
      end if
      
      ! Update positions using half-step velocities
      NB%pos(i, 1:ndim) = NB%pos_0(i, 1:ndim) + NB%vel_h(i, 1:ndim) * dt
      
      ! Calculate kinetic energy using full-step velocities
      NB%kinetic_energy(i) = 0.5_wp * NB%mass(i) * sum(NB%vel(i, 1:ndim)**2)
      NB%sum_energy(i) = NB%kinetic_energy(i) + NB%potential_energy(i)
      
      ! Store current values as old values for next iteration
      NB%pos_0(i, 1:ndim) = NB%pos(i, 1:ndim)
      NB%vel_0(i, 1:ndim) = NB%vel(i, 1:ndim)
    end do
    !$acc end parallel loop
    
    !$acc end data
  end subroutine update_positions_llmacc

  subroutine softening_length_setup_llmacc(NB)
    implicit none
    type(NBodies_type), intent(inout) :: NB
    real(kind=wp) :: E_total, KE_max, KE_initial, PE_initial
    real(kind=wp) :: v_max, r_min
    
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
    softening_length = 0.1_wp * r_min
    
    ! Set softening length squared based on radius_0
    softening_length_squared = softening_length**2
    print *, "Estimated softening length: ", softening_length
  end subroutine softening_length_setup_llmacc

end module llmacc_nbodies