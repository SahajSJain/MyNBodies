module nbodies
  use, intrinsic :: iso_fortran_env
  implicit none
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Precision stuff:
  ! Define working precision based on preprocessor flags
#ifdef P16
  ! Half precision (if supported)
  integer, parameter :: wp = real16
  character(len=*), parameter :: precision_name = "half (16-bit)"
#elif defined(P32)
  ! Single precision
  integer, parameter :: wp = real32
  character(len=*), parameter :: precision_name = "single (32-bit)"
#elif defined(P64)
  ! Double precision (default)
  integer, parameter :: wp = real64
  character(len=*), parameter :: precision_name = "double (64-bit)"
#elif defined(P128)
  ! Quadruple precision
  integer, parameter :: wp = real128
  character(len=*), parameter :: precision_name = "quadruple (128-bit)"
#else
  ! Default to double precision
  integer, parameter :: wp = real64
  character(len=*), parameter :: precision_name = "double (64-bit) [default]"
#endif

  ! Useful constants in working precision
  real(kind=wp), parameter :: pi = 4.0_wp * atan(1.0_wp)
  real(kind=wp), parameter :: e = exp(1.0_wp)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! DIMENSION DEPENDENT VARIABLES: :
  ! Dimension of the simulation (2D or 3D)
  integer :: ndim
  integer :: n_bodies ! Number of bodies in the simulation
  ! Positions x, y, z
  ! first dimension is body index, second dimension is x,y,(z)
  ! These are 2d arrays for easier extension to 3D
  ! First index: body number, second index: dimension (1=x,2=y,3=z)
  real(kind=wp), allocatable :: pos(:, :), pos_0(:, :) ! current and old positions
  real(kind=wp), allocatable :: vel(:, :), vel_0(:, :), vel_h(:, :) ! velocities: current, old, half-step
  real(kind=wp), allocatable :: accel(:, :) ! accelerations
  real(kind=wp), allocatable :: force(:, :) ! force vectors
  real(kind=wp), allocatable :: bc_factor(:, :) ! boundary condition factors
  ! random number temporary variables
  real(kind=wp), allocatable :: rand_vec(:) ! temporary random vector for initialization
  integer :: INITIALIZATION_TYPE ! 1=random, 2=ordered placement
  integer :: nbx_init(1:3) ! number of bodies in each dimension for ordered placement
  ! With these, the domain ranges from (-L/2) to (L/2) in each dimension
  ! instead of (0 to L) for easier centering
  real(kind=wp) :: L_init(1:3) ! box dimensions for placement
  real(kind=wp) :: L_bound(1:3) ! box dimensions for boundary conditions
  ! NON-DIMENSION DEPENDENT VARIABLES:
  real(kind=wp), allocatable :: mass(:), mass_inv(:), radius(:) !m: mass, radius: size of each body
  real(KIND=wp), allocatable :: kinetic_energy(:) ! Kinetic energy at each time step 
  real(kind=wp) :: total_kinetic_energy ! total kinetic energy at each time step 
  real(kind=wp) :: vel_init ! initial velocity scale for random initialization
  real(kind=wp) :: mass_0 ! mass input from the user
  real(kind=wp) :: radius_0 ! initial mean radius input from the user
  real(kind=wp) :: radius_var ! size variation scale for random initialization:
  ! size = size_0 * rand number ranging from (1 - radius_var) to (1 + radius_var)
  real(kind=wp) :: Gravitational_Constant ! Gravitational constant in working precision
  real(kind=wp) :: dt, dt_half ! time step for integration
  real(kind=wp) :: nsteps ! number of time steps to simulate
  real(kind=wp), parameter :: epsilon = 1.0e-10_wp ! small number to avoid singularities
  integer :: BOUNDARY_CONDITION_TYPE ! 1=open, 2=reflective, 3=periodic (not implemented yet)
  logical :: REFLECTIVE_BC ! flag for reflective boundary conditions

contains
  ! Additional module procedures can be added here
  ! Allocation= subroutine to allocate arrays based on n_bodies
  ! this will also initialize masses and radii based on radius_0 and radius_var
  subroutine allocate_arrays(num_bodies)
    integer, intent(in) :: num_bodies
    integer :: i
    n_bodies = num_bodies
    if (ndim /= 2 .and. ndim /= 3) then
      ndim = 2 ! default to 2D
    end if
    ! Dimension dependent allocations
    allocate(pos(n_bodies, ndim), pos_0(n_bodies, ndim), &
      vel(n_bodies, ndim), vel_0(n_bodies, ndim), vel_h(n_bodies, ndim), &
      accel(n_bodies, ndim), force(n_bodies, ndim), &
      bc_factor(n_bodies, ndim))
    ! Non-dimension dependent allocations
    allocate(mass(n_bodies), mass_inv(n_bodies), radius(n_bodies))
    allocate(rand_vec(ndim)) ! random vector for initialization
    ! Initialize radii with variation
    radius_var = max(0.0_wp, min(radius_var, 0.9_wp)) ! clamp radius_var between 0 and 0.9
    do i = 1, n_bodies
      call random_number(rand_val) ! rand_val in [0,1]
      ! Scale to range [(1-radius_var), (1+radius_var)]
      radius(i) = radius_0 &
        * ((1.0_wp - radius_var) + rand_val * 2.0_wp * radius_var)
      ! initialize velocities too between -v0 to v0
      call random_number(rand_vec)
      vel(i, 1:ndim) = (rand_vec - 0.5_wp) * 2.0_wp * vel_init
    end do
    ! for mass, body of radius r has mass proportional to r^3 (volume)
    ! and body of radius r has mass proportional to r^2 (area) in 2D
    ! and body of radius radius_0 has mass mass_0
    do i = 1, n_bodies
      mass(i) = mass_0 * (radius(i) / radius_0)**ndim
    end do
    mass_inv = 1.0_wp / (mass + epsilon) ! to avoid division by zero
  end subroutine allocate_arrays

  ! Deallocation subroutine
  subroutine deallocate_arrays()
    ! Deallocate dimension dependent arrays
    deallocate(pos, pos_0, vel, vel_0, vel_h, accel, force, bc_factor)
    ! Deallocate non-dimension dependent arrays
    deallocate(mass, mass_inv, radius)
    ! Deallocate random vector
    deallocate(rand_vec)
    ! Reset n_bodies to 0
    n_bodies = 0
  end subroutine deallocate_arrays

  subroutine ordered_placement()
    ! Example subroutine for ordered placement of bodies
    integer :: i, j, k, idx
    integer :: count
    ! Ensure ndim is valid and change n_bodies for lattice placement 
    if (ndim < 1) ndim = 1 ! minimum dimension is 1D
    if (ndim > 3) ndim = 3 ! maximum dimension is 3D 
    ! Change number of bodies based on user input for lattice size
    if (ndim < 3) then
      nbx_init(3) = 1 ! for 2D case
    end if
    if (ndim < 2) then
      nbx_init(2) = 1 ! for 1D case (not really used here)
    end if
    n_bodies = nbx_init(1) * nbx_init(2) * nbx_init(3) ! total number of bodies

    call allocate_arrays(n_bodies)
    count = 0

    select case (ndim)
     case (1)
      do i = 0, nbx_init(1) - 1
        count = count + 1
        pos(count, 1) = (i + 0.5_wp) * (L_init(1) / nbx_init(1))
      end do
     case (2)
      do i = 0, nbx_init(1) - 1
        do j = 0, nbx_init(2) - 1
          count = count + 1
          pos(count, 1) = (i + 0.5_wp) * (L_init(1) / nbx_init(1))
          pos(count, 2) = (j + 0.5_wp) * (L_init(2) / nbx_init(2))
        end do
      end do
     case (3)
      do i = 0, nbx_init(1) - 1
        do j = 0, nbx_init(2) - 1
          do k = 0, nbx_init(3) - 1
            count = count + 1
            pos(count, 1) = (i + 0.5_wp) * (L_init(1) / nbx_init(1))
            pos(count, 2) = (j + 0.5_wp) * (L_init(2) / nbx_init(2))
            pos(count, 3) = (k + 0.5_wp) * (L_init(3) / nbx_init(3))
          end do
        end do
      end do
     case default
      print *, "Error: ndim must be 1, 2, or 3"
      stop
    end select
    ! Shift positions from [0, L] to [-L/2, L/2]
    do i = 1, ndim
      pos(:, i) = pos(:, i) - L_init(i) / 2.0_wp
    end do
    ! initialize old positions
    pos_0 = pos
  end subroutine ordered_placement

  subroutine random_placement()
    ! Example subroutine for random placement of bodies
    integer :: i
    if (ndim < 1) ndim = 1 ! minimum dimension is 1D
    if (ndim > 3) ndim = 3 ! maximum dimension is 3D 
    call allocate_arrays(n_bodies)
    do i = 1, n_bodies
      call random_number(rand_vec(1:ndim))
      pos(i, 1:ndim) = rand_vec(1:ndim) * L_init(1:ndim)
    end do
    ! Shift positions from [0, L] to [-L/2, L/2]
    do i = 1, ndim
      pos(:, i) = pos(:, i) - L_init(i) / 2.0_wp
    end do
    ! lets not worry about overlaps for now
    ! Smoothing and random velocity initialization might ensure things dont get
    ! too close.

    ! initialize old positions
    pos_0 = pos
  end subroutine random_placement
end module nbodies
