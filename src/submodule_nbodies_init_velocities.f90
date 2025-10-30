submodule (nbodies) init_vel
  implicit none
contains

  module subroutine initialize_velocities(NB)
    implicit none
    type(NBodies_type), intent(inout) :: NB

    ! Initialize velocities based on VELOCITY_INITIALIZATION_TYPE
    select case (VELOCITY_INITIALIZATION_TYPE)
     case (0)
      call read_velocities_from_file(NB)
     case (1)
      call random_velocities_random_speed(NB)
     case (2)
      call solid_body_rotation_velocities(NB)
     case (3)
      call random_velocities_fixed_speed(NB)
     case default
      print *, "Error: Unknown VELOCITY_INITIALIZATION_TYPE ", &
        VELOCITY_INITIALIZATION_TYPE
      stop
    end select
    ! Initialize old velocities
    NB%vel_0 = NB%vel
  end subroutine initialize_velocities

  subroutine solid_body_rotation_velocities(NB)
    implicit none
    type(NBodies_type), intent(inout) :: NB

    ! This is only for 2D and 3D cases
    integer :: i
    real(kind=wp) :: vel_magnitude ! magnitude of velocity at position

    do i = 1, NB%n_bodies
      ! u = -omega * y, v = omega * x for 2D rotation about z-axis
      NB%vel(i, 1) = -omega_init * NB%pos(i, 2)
      NB%vel(i, 2) = omega_init * NB%pos(i, 1)
      if (ndim == 3) then
        ! w = 0 for rotation about z-axis
        NB%vel(i, 3) = 0.0_wp
      end if
      vel_magnitude = sqrt(sum(NB%vel(i, 1:ndim)**2))
      ! now add size_var based random perturbation to velocities
      call random_number(NB%rand_vec)
      ! add random perturbation scaled by vel_var (between 0 and 1) of the local
      ! velocity magnitude
      NB%vel(i, 1:ndim) = NB%vel(i, 1:ndim) + &
        vel_var * vel_magnitude * (NB%rand_vec(1:ndim) - 0.5_wp) * 2.0_wp
    end do
  end subroutine solid_body_rotation_velocities

  subroutine random_velocities_fixed_speed(NB)
    implicit none
    type(NBodies_type), intent(inout) :: NB

    ! This is only for 2D and 3D cases
    integer :: i
    real(kind=wp) :: vel_magnitude ! magnitude of velocity at position

    do i = 1, NB%n_bodies
      ! add random perturbation vel_var
      ! if vel_var = 2.0, then velocity magnitude range from -2.0 to +2.0
      ! Procedure:
      ! First generate random vector in [-1,1] in each dimension
      call random_number(NB%rand_vec)
      NB%vel(i, 1:ndim) = vel_var * (NB%rand_vec(1:ndim) - 0.5_wp) * 2.0_wp
      ! Then normalize to 1
      vel_magnitude = sqrt(sum(NB%vel(i, 1:ndim)**2)) + small_number ! avoid division by zero
      NB%vel(i, 1:ndim) = NB%vel(i, 1:ndim) / vel_magnitude
      ! Finally scale to desired vel_var magnitude
      NB%vel(i, 1:ndim) = NB%vel(i, 1:ndim) * vel_var
    end do
  end subroutine random_velocities_fixed_speed

  subroutine random_velocities_random_speed(NB)
    implicit none
    type(NBodies_type), intent(inout) :: NB

    ! This is only for 2D and 3D cases
    integer :: i
    real(kind=wp) :: vel_magnitude, rand_val ! magnitude of velocity at position

    do i = 1, NB%n_bodies
      ! add random perturbation vel_var
      ! if vel_var = 2.0, then velocity magnitude range from -2.0 to +2.0
      ! Procedure:
      ! First generate random vector in [-1,1] in each dimension
      call random_number(NB%rand_vec)
      NB%vel(i, 1:ndim) = vel_var * (NB%rand_vec(1:ndim) - 0.5_wp) * 2.0_wp
      ! Then normalize to 1
      vel_magnitude = sqrt(sum(NB%vel(i, 1:ndim)**2)) + small_number ! avoid division by zero
      NB%vel(i, 1:ndim) = NB%vel(i, 1:ndim) / vel_magnitude
      ! Then select a random speed in [0, vel_var]
      call random_number(rand_val) ! rand_val in [0,1]
      vel_magnitude = rand_val * vel_var
      NB%vel(i, 1:ndim) = NB%vel(i, 1:ndim) * vel_magnitude
    end do
  end subroutine random_velocities_random_speed
  subroutine read_velocities_from_file(NB)
    implicit none
    type(NBodies_type), intent(inout) :: NB

    ! Local variables
    integer :: i, ios
    integer :: file_unit
    character(len=256) :: filename

    ! Set filename
    filename = 'data/velocities_ascii.dat'

    ! Open velocity file
    open(newunit=file_unit, file=trim(filename), status='old', action='read', iostat=ios)

    ! Read velocities
    do i = 1, NB%n_bodies
      read(file_unit, *) NB%vel(i, 1), NB%vel(i, 2), NB%vel(i, 3)
    end do

    ! Close file
    close(file_unit)

    ! Initialize old velocities
    NB%vel_0 = NB%vel

  end subroutine read_velocities_from_file
end submodule init_vel
