module serial_nbodies
  use nbodies
  implicit none
contains
  subroutine compute_forces_serial(NB)
    implicit none
    type(NBodies_type), intent(inout) :: NB
    integer :: i, j
    real(kind=wp) :: pos_diff(ndim)
    real(kind=wp) :: dist, dist_cubed ! distance and distance cubed between bodies
    real(kind=wp) :: dist_squared, dist_squared_soft ! distance squared and softened distance squared
    ! Initialize forces to zero
    NB%force(:, :) = 0.0_wp
    NB%potential_energy = 0.0_wp ! Reset potential energy array
    NB%total_potential_energy = 0.0_wp ! Reset total potential energy
    NB%kinetic_energy = 0.0_wp ! Reset kinetic energy array
    NB%total_kinetic_energy = 0.0_wp ! Reset total kinetic energy
    NB%total_energy = 0.0_wp ! Reset total energy
    ! Compute pairwise gravitational forces
    do i = 1, NB%n_bodies
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
          NB%force(i, 1:ndim) = NB%force(i, 1:ndim) + &
            Gravitational_Constant * NB%mass(i) * NB%mass(j) * pos_diff(1:ndim) / dist_cubed

          ! Calculate potential energy: half to avoid double counting
          NB%potential_energy(i) = NB%potential_energy(i) - &
            0.50_wp*Gravitational_Constant * NB%mass(i) * NB%mass(j) / dist
        end if
      end do
    end do

    ! Calculate total potential energy
    NB%total_potential_energy = sum(NB%potential_energy)

  end subroutine compute_forces_serial

  subroutine get_global_metrics_serial(NB)
    implicit none
    type(NBodies_type), intent(inout) :: NB
    integer :: i

    NB%total_kinetic_energy = 0.0_wp
    NB%total_potential_energy = 0.0_wp
    NB%total_energy = 0.0_wp
    NB%center_of_mass = 0.0_wp
    NB%center_of_mass_velocity = 0.0_wp

    do i = 1, NB%n_bodies
      NB%kinetic_energy(i) = 0.5_wp * NB%mass(i) * sum(NB%vel(i, 1:ndim)**2)
      NB%sum_energy(i) = NB%kinetic_energy(i) + NB%potential_energy(i)
    end do

    do i = 1, NB%n_bodies
      NB%total_kinetic_energy = NB%total_kinetic_energy + NB%kinetic_energy(i)
      ! multiply potential energy by 0.5 to account for double counting
      NB%total_potential_energy = NB%total_potential_energy &
        + NB%potential_energy(i) ! already accounted for in compute_forces_serial
      NB%center_of_mass(1:ndim) = NB%center_of_mass(1:ndim) + &
        NB%mass(i) * NB%pos(i, 1:ndim) * NB%total_mass_inv
      NB%center_of_mass_velocity(1:ndim) = NB%center_of_mass_velocity(1:ndim) + &
        NB%mass(i) * NB%vel(i, 1:ndim) * NB%total_mass_inv
    end do
    NB%total_energy = NB%total_kinetic_energy + NB%total_potential_energy
  end subroutine get_global_metrics_serial

  subroutine update_positions_serial(NB)
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
    ! Reset energies
    NB%kinetic_energy = 0.0_wp
    NB%sum_energy = 0.0_wp
    ! Calculate acceleration from forces
    do i = 1, NB%n_bodies
      NB%accel(i, 1:ndim) = NB%force(i, 1:ndim) * NB%mass_inv(i)
    end do
    ! Update velocities (full step for diagnostics, half step for position update)
    do i = 1, NB%n_bodies
      ! Half-step velocity for position update
      NB%vel_h(i, 1:ndim) = NB%vel_0(i, 1:ndim) + NB%accel(i, 1:ndim) * dt_half
      ! Full-step velocity for energy calculation
      NB%vel(i, 1:ndim) = NB%vel_0(i, 1:ndim) + NB%accel(i, 1:ndim) * dt
    end do
    if (REFLECTIVE_BC) then
      do i = 1, NB%n_bodies
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
      end do
    end if
    ! Update positions using half-step velocities
    do i = 1, NB%n_bodies
      NB%pos(i, 1:ndim) = NB%pos_0(i, 1:ndim) + NB%vel_h(i, 1:ndim) * dt
    end do

    ! Calculate kinetic energy using full-step velocities
    do i = 1, NB%n_bodies
      NB%kinetic_energy(i) = 0.5_wp * NB%mass(i) * sum(NB%vel(i, 1:ndim)**2)
      NB%sum_energy(i) = NB%kinetic_energy(i) + NB%potential_energy(i)
    end do

! Store current values as old values for next iteration
    do i = 1, NB%n_bodies
      NB%pos_0(i, 1:ndim) = NB%pos(i, 1:ndim)
      NB%vel_0(i, 1:ndim) = NB%vel(i, 1:ndim)
    end do
  end subroutine update_positions_serial
  subroutine softening_length_setup_serial(NB)
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
    softening_length = 0.1_wp * r_min
    ! Set softening length squared based on radius_0
    softening_length_squared = softening_length**2
    print *, "Estimated softening length: ", softening_length
  end subroutine softening_length_setup_serial

end module serial_nbodies

 ! Notes for reflective BCs:
 ! Now based on old positions, if out of domain, apply reflective BCs
 ! Note we only have reflective BCs here for simplicity since I dont know
 ! how to model
 ! gravitational
 ! forces with periodic BCs lol
 ! if x position is greater than Lx, and velocity is positive, reverse
 ! velocity
 ! but not using if statement to avoid branching
 ! Merge:     RESULT = MERGE(TSOURCE, FSOURCE, MASK)
 !            Select values from two arrays according to a logical mask.
 ! The result is equal to
 ! TSOURCE if MASK is .TRUE., or equal to FSOURCE if it is .FALSE..
 ! bc_factor = 1.0  means no change, bc_factor = -1.0 means reverse
 ! velocity
 ! first bracket : x(i) > Lx/2 .and. vx_temp > 0.0_wp (heading out of box on
 ! right towards +ve side)
 ! second bracket: x(i) < -Lx/2 .and. vx_temp < 0.0_wp (or heading out of
 ! box on left towards -ve side)
 ! if either condition is true, we want to reverse velocity:bc_factor =
 ! 1-2=-1.0_wp
 ! else bc_factor = 1-2*0=1.0_wp
 !     if (REFLECTIVE_BC) then
 !   do i = 1, NB%n_bodies
 !     do nd = 1, ndim
 !       NB%bc_factor(i, nd) = (1.0_wp &
 !         - iReflective_BC(nd) * 2.0_wp * merge(1.0_wp, 0.0_wp, &
 !         ((NB%pos(i, nd) > +L_bound(nd) * 0.50_wp .and. NB%vel(i, nd) > &
 !         0.0_wp) &
 !         .or. &
 !         (NB%pos(i, nd) < -L_bound(nd) * 0.50_wp .and. NB%vel(i, nd) < &
 !         0.0_wp))))
 !       NB%vel_h(i, nd) = NB%vel_h(i, nd) * NB%bc_factor(i, nd)
 !       NB%vel(i, nd) = NB%vel(i, nd) * NB%bc_factor(i, nd)
 !     end do
 !   end do
 ! end if
