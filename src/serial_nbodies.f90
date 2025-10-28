module serial_nbodies
  use nbodies
  implicit none
contains
  subroutine compute_forces_serial()
    integer :: i, j
    real(kind=wp) :: pos_diff(ndim)
    real(kind=wp) :: dist, dist_cubed ! distance and distance cubed between bodies
    ! Initialize forces to zero
    force(:, :) = 0.0_wp
    ! Compute pairwise gravitational forces
    do i = 1, n_bodies
      do j = 1, n_bodies
        if (i /= j) then
          ! force on particle i due to particle j
          ! i will be pulled towards j
          ! if x(j) > x(i) then fx(i) += positive value i.e. pull right towards
          ! j
          ! therefore x_diff = x(j) - x(i) is correct.
          ! alternatively we can do dx = x(i) - x(j) and then subtract the
          ! force.
          ! in vectorized form:
          pos_diff = pos(j, 1:ndim) - pos(i, 1:ndim)
          dist = sqrt(sum(pos_diff(1:ndim)**2) + epsilon) ! Softening to avoid singularity
          dist_cubed = dist**3

          force(i, 1:ndim) = force(i, 1:ndim) &
            + Gravitational_Constant * mass(i) * mass(j) * (pos_diff(1:ndim) / dist_cubed)
        end if
      end do
    end do 
  end subroutine compute_forces_serial

  subroutine update_positions_serial()
    integer :: i 
    total_kinetic_energy = 0.0_wp 
    ! Use Verlet integration to update positions and velocities
    ! First update accelerations based on forces
    ! Update velocities and positions based on computed forces
    ! in kick-drift-kick scheme form
    ! First update accelerations:
    ! Note: We have dimensional arrays: 
    !     force(ndim,n_bodies) 
    !     accel(ndim,n_bodies)
    !     vel(ndim,n_bodies)
    !     vel_0(ndim,n_bodies)
    !     vel_h(ndim,n_bodies) 
    !     pos(ndim,n_bodies) 
    !     pos_0(ndim,n_bodies) 
    !   and non-dimensional arrays:
    !     mass(n_bodies) 
    !     mass_inv(n_bodies) 
    do i = 1, n_bodies
      ! acceleration = force / mass 
      accel(i, 1:ndim) = force(i, 1:ndim) * mass_inv(i)
    end do
    ! Now update half step velocities and full step velocities:
    do i=1, n_bodies
      vel_h(i, 1:ndim) = vel_0(i, 1:ndim) + accel(i, 1:ndim) * dt_half 
      vel(i, 1:ndim) = vel_0(i, 1:ndim) + accel(i, 1:ndim) * dt 
    end do 
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
    ! if either condition is true, we want to reverse velocity:bc_factor = 1-2=-1.0_wp
    ! else bc_factor = 1-2*0=1.0_wp
    IF(REFLECTIVE_BC) THEN
      do i = 1, n_bodies
        do nd = 1, ndim 
          bc_factor(i, nd) = (1.0_wp - 2.0_wp * merge(1.0_wp, 0.0_wp, &
            ((pos(i, nd) > + L_bound(nd)*0.50_wp .and. vel(i, nd) > 0.0_wp) .or. &
             (pos(i, nd) < - L_bound(nd)*0.50_wp .and. vel(i, nd) < 0.0_wp))))   
          vel_h(i, nd) = vel_h(i, nd) * bc_factor(i, nd) 
          vel(i, nd) = vel(i, nd) * bc_factor(i, nd) 
        enddo 
      enddo 
    end if 
    ! Now update positions based on half-step velocities
    do i = 1, n_bodies
      ! Using temporary variables for clarity
      pos(i, 1:ndim) = pos_0(i, 1:ndim) + vel_h(i, 1:ndim) * dt 
    end do 
    ! Now compute kinetic energy for diagnostics 
    do i = 1, n_bodies
      kinetic_energy(i) = 0.5_wp * mass(i) * sum(vel(i, 1:ndim)**2)
      total_kinetic_energy = total_kinetic_energy + kinetic_energy(i) 
    end do  
    ! Finally, update old positions and velocities for next iteration
    do i = 1, n_bodies
      pos_0(i, 1:ndim) = pos(i, 1:ndim)
      vel_0(i, 1:ndim) = vel(i, 1:ndim)
    enddo 
  end subroutine update_positions_serial
end module serial_nbodies
