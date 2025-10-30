module nbodies
  use, intrinsic :: iso_fortran_env
  implicit none

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! PRECISION CONFIGURATION
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

  ! Mathematical constants in working precision
  real(kind=wp), parameter :: pi = 4.0_wp * atan(1.0_wp)
  real(kind=wp), parameter :: e = exp(1.0_wp)
  real(kind=wp), parameter :: small_number = 2.0_wp * epsilon(1.0_wp) ! 2x machine small_number for working precision
  real(kind=wp) :: softening_length = 1.e-5_wp ! default softening length
  real(kind=wp) :: softening_length_squared = 1.e-10_wp ! default softening length squared
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! SIMULATION PARAMETERS (kept as module variables)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Dimension
  integer :: ndim                   ! Dimension of simulation (2D or 3D)
  ! input-output file units
  ! we will use newunit feature of Fortran 2008 to avoid unit conflicts
  integer :: ifinput_unit, iftimehistory_unit, ifdump_unit, ifforcematrix_unit
  ! Physical constants
  real(kind=wp) :: Gravitational_Constant  ! Gravitational constant
  ! Time integration parameters
  real(kind=wp) :: dt               ! Time step
  real(kind=wp) :: dt_half          ! Half time step (for leapfrog)
  real(KIND=wp) :: dt_squared        !
  real(kind=wp) :: dt_min, dt_max   ! Minimum and maximum time step
  integer :: nsteps                 ! Number of time steps
  integer :: ndump                  ! Dump frequency for output
  integer :: nprint                 ! Print frequency for diagnostics
  real(kind=wp) :: time             ! Current simulation time
  integer :: istep, idump           ! current time step and dump count

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! INITIALIZATION PARAMETERS: all in input.dat
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Position initialization
  integer :: INITIALIZATION_TYPE      ! 1=random, 2=ordered placement
  integer :: nbx_init(1:3)            ! Number of bodies in each dimension for ordered
  real(kind=wp) :: L_init(1:3)        ! Box dimensions for placement

  ! Velocity initialization
  integer :: VELOCITY_INITIALIZATION_TYPE   ! 1=random, 2=solid body rotation
  real(kind=wp) :: vel_var                  ! Initial velocity scale for random
  real(kind=wp) :: omega_init               ! Angular velocity for rotation

  ! Mass and size parameters
  real(kind=wp) :: mass_0           ! Base mass value
  real(kind=wp) :: radius_0         ! Base radius value
  real(kind=wp) :: radius_var       ! Radius variation scale

  ! Boundary conditions
  integer :: BOUNDARY_CONDITION_TYPE  ! 1=open, 2=reflective, 3=periodic
  logical :: REFLECTIVE_BC            ! Flag for reflective BC
  real(kind=wp) :: L_bound(1:3)       ! Box dimensions for boundaries
  real(kind=wp) :: iReflective_BC(1:3)      ! Reflective BC flags per dimension
  ! Force Diagnose
  integer :: Force_Diagnose           ! 0=off, 1=on (only first timestep)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NBODIES TYPE DEFINITION
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  type :: NBodies_type
    ! Number of bodies
    integer :: n_bodies

    ! DYNAMICAL ARRAYS - BODY PROPERTIES
    ! Position and velocity arrays (n_bodies x ndim)
    real(kind=wp), allocatable :: pos(:,:)       ! Current positions
    real(kind=wp), allocatable :: pos_0(:,:)     ! Previous positions
    real(kind=wp), allocatable :: pos_00(:,:)    ! Previous positions
    real(kind=wp), allocatable :: vel(:,:)       ! Current velocities
    real(kind=wp), allocatable :: vel_0(:,:)     ! Previous velocities
    real(kind=wp), allocatable :: vel_00(:,:)    ! Previous velocities
    real(kind=wp), allocatable :: vel_h(:,:)     ! Half-step velocities

    ! Force and acceleration arrays (n_bodies x ndim)
    real(kind=wp), allocatable :: accel(:,:)     ! Accelerations
    real(kind=wp), allocatable :: accel_0(:,:)    ! Previous accelerations
    real(kind=wp), allocatable :: accel_00(:,:)   ! Previous accelerations
    real(kind=wp), allocatable :: force(:,:)     ! Force vectors
    real(kind=wp), allocatable :: bc_factor(:,:) ! Boundary condition factors

    ! Mass and geometry arrays (n_bodies)
    real(kind=wp), allocatable :: mass(:)        ! Body masses
    real(kind=wp), allocatable :: mass_inv(:)    ! Inverse masses
    real(kind=wp), allocatable :: radius(:)      ! Body radii

    ! ENERGY AND GLOBAL PROPERTIES
    ! Per-body energy arrays (n_bodies)
    real(kind=wp), allocatable :: kinetic_energy(:)    ! KE per body
    real(kind=wp), allocatable :: potential_energy(:)  ! PE per body
    real(kind=wp), allocatable :: sum_energy(:)        ! Total E per body

    ! Global system properties
    real(kind=wp) :: total_kinetic_energy      ! Total system KE
    real(kind=wp) :: total_potential_energy    ! Total system PE
    real(kind=wp) :: total_energy              ! Total system energy
    real(kind=wp) :: total_mass                ! Total system mass
    real(kind=wp) :: total_mass_inv            ! Inverse total mass
    real(kind=wp) :: total_energy_initial      ! Maximum velocity magnitude
    ! Center of mass arrays (ndim)
    real(kind=wp), allocatable :: center_of_mass(:)          ! COM position
    real(kind=wp), allocatable :: center_of_mass_velocity(:) ! COM velocity

    ! Temporary arrays for initialization
    real(kind=wp), allocatable :: rand_vec(:)   ! Random vector (ndim)
    real(kind=wp) :: min_vel, max_vel, min_accel, max_accel ! min-max magnitudes
    real(kind=wp) :: char_dist                   ! min-max distances
    ! Force Diagnoses
    real(kind=wp), allocatable :: Force_matrix(:,:,:)  ! matrix to store forces for diagnosis
    ! (to,from,ndim)
  contains
    procedure :: allocate => allocate_nbodies
    procedure :: deallocate => deallocate_nbodies
    final :: finalize_nbodies
  end type NBodies_type
  interface
    module subroutine initialize_velocities(NB)
      type(NBodies_type), intent(inout) :: NB
    end subroutine initialize_velocities
    module subroutine initialize_positions(NB)
      type(NBodies_type), intent(inout) :: NB
    end subroutine initialize_positions
  end interface

contains

  subroutine allocate_nbodies(this, num_bodies)
    implicit none
    class(NBodies_type), intent(inout) :: this
    integer, intent(in) :: num_bodies
    integer :: i
    real(kind=wp) :: rand_val

    this%n_bodies = num_bodies

    ! Validate ndim: not necessary but keep it.
    if (ndim /= 2 .and. ndim /= 3) then
      ndim = 2 ! default to 2D
    end if

    ! Allocate position and velocity arrays
    allocate(this%pos(this%n_bodies, ndim), source=0.0_wp)
    allocate(this%pos_0(this%n_bodies, ndim), source=0.0_wp)
    allocate(this%pos_00(this%n_bodies, ndim), source=0.0_wp)
    allocate(this%vel(this%n_bodies, ndim), source=0.0_wp)
    allocate(this%vel_0(this%n_bodies, ndim), source=0.0_wp)
    allocate(this%vel_00(this%n_bodies, ndim), source=0.0_wp)
    allocate(this%vel_h(this%n_bodies, ndim), source=0.0_wp)

    ! Allocate force and acceleration arrays
    allocate(this%accel(this%n_bodies, ndim), source=0.0_wp)
    allocate(this%accel_0(this%n_bodies, ndim), source=0.0_wp)
    allocate(this%accel_00(this%n_bodies, ndim), source=0.0_wp)
    allocate(this%force(this%n_bodies, ndim), source=0.0_wp)
    allocate(this%bc_factor(this%n_bodies, ndim), source=0.0_wp)

    ! Allocate mass and geometry arrays
    allocate(this%mass(this%n_bodies), source=0.0_wp)
    allocate(this%mass_inv(this%n_bodies), source=0.0_wp)
    allocate(this%radius(this%n_bodies), source=0.0_wp)

    ! Allocate energy arrays
    allocate(this%kinetic_energy(this%n_bodies), source=0.0_wp)
    allocate(this%potential_energy(this%n_bodies), source=0.0_wp)
    allocate(this%sum_energy(this%n_bodies), source=0.0_wp)

    ! Allocate center of mass arrays
    allocate(this%center_of_mass(ndim), source=0.0_wp)
    allocate(this%center_of_mass_velocity(ndim), source=0.0_wp)

    ! Allocate temporary arrays
    allocate(this%rand_vec(ndim), source=0.0_wp)
    ! Allocate force diagnosis matrix if needed
    if(Force_Diagnose == 1) then
      allocate(this%Force_matrix(this%n_bodies, this%n_bodies, ndim), source=0.0_wp)
    end if
    ! Initialize radii with variation
    radius_var = max(0.0_wp, min(radius_var, 0.9_wp)) ! clamp between 0 and 0.9
    do i = 1, this%n_bodies
      call random_number(rand_val) ! rand_val in [0,1]
      ! Scale to range [(1-radius_var), (1+radius_var)]
      this%radius(i) = radius_0 * ((1.0_wp - radius_var) + rand_val * 2.0_wp * radius_var)
    end do

    ! Initialize masses based on radii (assuming constant density)
    do i = 1, this%n_bodies
      this%mass(i) = mass_0 * (this%radius(i) / radius_0)**ndim
    end do

    ! Calculate inverse masses and total mass
    this%mass_inv = 1.0_wp / (this%mass + small_number)
    this%total_mass = sum(this%mass)
    this%total_mass_inv = 1.0_wp / (this%total_mass + small_number)

    ! Initialize global properties
    this%total_kinetic_energy = 0.0_wp
    this%total_potential_energy = 0.0_wp
    this%total_energy = 0.0_wp

  end subroutine allocate_nbodies

  subroutine deallocate_nbodies(this)
    implicit none
    class(NBodies_type), intent(inout) :: this

    ! Deallocate position and velocity arrays
    if (allocated(this%pos)) deallocate(this%pos)
    if (allocated(this%pos_0)) deallocate(this%pos_0)
    if (allocated(this%vel)) deallocate(this%vel)
    if (allocated(this%vel_0)) deallocate(this%vel_0)
    if (allocated(this%vel_h)) deallocate(this%vel_h)

    ! Deallocate force and acceleration arrays
    if (allocated(this%accel)) deallocate(this%accel)
    if (allocated(this%force)) deallocate(this%force)
    if (allocated(this%bc_factor)) deallocate(this%bc_factor)

    ! Deallocate mass and geometry arrays
    if (allocated(this%mass)) deallocate(this%mass)
    if (allocated(this%mass_inv)) deallocate(this%mass_inv)
    if (allocated(this%radius)) deallocate(this%radius)

    ! Deallocate energy arrays
    if (allocated(this%kinetic_energy)) deallocate(this%kinetic_energy)
    if (allocated(this%potential_energy)) deallocate(this%potential_energy)
    if (allocated(this%sum_energy)) deallocate(this%sum_energy)

    ! Deallocate center of mass and temporary arrays
    if (allocated(this%center_of_mass)) deallocate(this%center_of_mass)
    if (allocated(this%center_of_mass_velocity)) deallocate(this%center_of_mass_velocity)
    if (allocated(this%rand_vec)) deallocate(this%rand_vec)

    ! Reset n_bodies
    this%n_bodies = 0
  end subroutine deallocate_nbodies

  subroutine finalize_nbodies(this)
    implicit none
    type(NBodies_type), intent(inout) :: this
    call this%deallocate()
  end subroutine finalize_nbodies

  subroutine initialize_simulation(NB)
    implicit none
    logical :: dir_exists
    integer :: ios
    type(NBodies_type), intent(inout) :: NB

    ! Read parameters from input.dat
    call read_input_parameters(NB)

    ! Initialize positions - this will also allocate arrays based on initialization type
    call initialize_positions(NB)

    ! Initialize velocities
    call initialize_velocities(NB)
    ! initialize masses
    call normalize_mass(NB)
    ! Create output directory if it doesn't exist
    INQUIRE (FILE='./NBDUMP/', EXIST=dir_exists)
    IF (.NOT. dir_exists) THEN
      CALL EXECUTE_COMMAND_LINE('mkdir -p ./NBDUMP/')
    END IF
    INQUIRE (FILE='./NBDUMPCSV/', EXIST=dir_exists)
    IF (.NOT. dir_exists) THEN
      CALL EXECUTE_COMMAND_LINE('mkdir -p ./NBDUMPCSV/')
    END IF
    open(newunit=iftimehistory_unit, file='Time_History.dat', status='replace', action='write', iostat=ios)
    if (ios /= 0) then
      print *, "Error: Cannot open input file ", trim('Time_History.dat')
      stop
    end if
    ! Write header for time history file (do this once before the time loop)
    if(ndim == 2) then
      write(iftimehistory_unit, '(A)') 'step, time, COM_x, COM_y, ' // &
        'COM_vx, COM_vy, KE_total, PE_total, E_total'
    else if(ndim == 3) then
      write(iftimehistory_unit, '(A)') 'step, time, COM_x, COM_y, COM_z, ' // &
        'COM_vx, COM_vy, COM_vz, KE_total, PE_total, E_total'
    end if

    ! Initialize half-step velocities for leapfrog integration
    NB%vel_h = NB%vel
    NB%pos_0 = NB%pos
    NB%pos_00 = NB%pos_0
    NB%vel_0 = NB%vel
    NB%vel_00 = NB%vel_0
    NB%accel_0 = NB%accel
    NB%accel_00 = NB%accel_0
    ! Initialize time variables
    time = 0.0_wp
    istep = 0
    dt_half = 0.5_wp * dt
    dt_squared = dt * dt
    dt_min = dt * 0.01_wp
    dt_max = dt * 100.0_wp
    ! Print summary of parameters read
    print *, "=== Input Parameters Read Successfully ==="
    print '(A,I3,A,I6)', "Dimensions: ", ndim, ", Bodies: ", NB%n_bodies
    print '(A,F8.3)', "Gravitational constant: ", Gravitational_Constant
    print '(A,F15.8,A,F8.3,A,F8.3)', "Mass: ", mass_0, ", Radius: ", radius_0, &
      ", Radius variation: ", radius_var
    print '(A,F8.6,A,I8)', "Time step: ", dt, ", Number of steps: ", nsteps
    print '(A,I2)', "Initialization type: ", INITIALIZATION_TYPE
    print '(A,I2)', "Velocity initialization type: ", VELOCITY_INITIALIZATION_TYPE
    print '(A,I2)', "Boundary condition type: ", BOUNDARY_CONDITION_TYPE
    print *, "========================================="

    print *, "=== Simulation Initialized Successfully ==="
    print '(A)', " Precision: " // trim(precision_name)
    print '(A,I8)', " Number of bodies: ", NB%n_bodies
    print '(A,I2)', " Dimensions: ", ndim
    print '(A,F12.6)', " Time step (dt): ", dt
    print '(A,F12.6)', " Half time step (dt/2): ", dt_half
    print '(A,F12.6)', " Minimum time step (dt_min): ", dt_min
    print '(A,F12.6)', " Maximum time step (dt_max): ", dt_max
    print '(A,I8)', " Total steps: ", nsteps
    print '(A,I8)', " Dump frequency: ", ndump
    print *, "==========================================="
  end subroutine initialize_simulation

  subroutine read_input_parameters(NB)
    implicit none
    type(NBodies_type), intent(inout) :: NB
    integer :: ios
    character(len=200) :: line

    open(newunit=ifinput_unit, file='input.dat', status='old', action='read', iostat=ios)
    if (ios /= 0) then
      print *, "Error: Cannot open input file ", trim('input.dat')
      stop
    end if

    ! Read dimensions and number of bodies
    read(ifinput_unit, *) ! Skip header line
    read(ifinput_unit, *) ndim, NB%n_bodies, Force_Diagnose

    ! Read gravitational constant
    read(ifinput_unit, *) ! Skip header line
    read(ifinput_unit, *) Gravitational_Constant

    ! Read mass and radius parameters
    read(ifinput_unit, *) ! Skip header line
    read(ifinput_unit, *) mass_0, radius_0, radius_var

    ! Read time parameters
    read(ifinput_unit, *) ! Skip header line
    read(ifinput_unit, *) dt, nsteps, nprint, ndump

    ! Read initialization type
    read(ifinput_unit, *) ! Skip header line
    read(ifinput_unit, *) INITIALIZATION_TYPE

    ! Read initial grid dimensions
    read(ifinput_unit, *) ! Skip header line
    read(ifinput_unit, *) nbx_init(1), nbx_init(2), nbx_init(3)

    ! Read initial box dimensions
    read(ifinput_unit, *) ! Skip header line
    read(ifinput_unit, *) L_init(1), L_init(2), L_init(3)

    ! Read velocity initialization type
    read(ifinput_unit, *) ! Skip header line
    read(ifinput_unit, *) VELOCITY_INITIALIZATION_TYPE

    ! Read velocity parameters
    read(ifinput_unit, *) ! Skip header line
    read(ifinput_unit, *) vel_var, omega_init

    ! Read boundary condition type
    read(ifinput_unit, *) ! Skip header line
    read(ifinput_unit, *) BOUNDARY_CONDITION_TYPE

    ! Read boundary dimensions
    read(ifinput_unit, *) ! Skip header line
    read(ifinput_unit, *) iReflective_BC(1), iReflective_BC(2), iReflective_BC(3)
    read(ifinput_unit, *) ! Skip header line
    read(ifinput_unit, *) L_bound(1), L_bound(2), L_bound(3)

    ! Close file
    close(ifinput_unit)

    ! Set derived parameters
    dt_half = 0.5_wp * dt
    dt_squared = dt * dt
    REFLECTIVE_BC = (BOUNDARY_CONDITION_TYPE == 2)
    if(REFLECTIVE_BC .eq. .false.) then
      iReflective_BC = 0.0_wp
    end if

  end subroutine read_input_parameters

  subroutine dump_data(NB)
    implicit none
    type(NBodies_type), intent(inout) :: NB

    character(len=256) :: filename
    character(len=7) :: step_str
    integer :: i, ios
    idump = idump + 1
    ! Create filename with 7-digit timestep number
    write(step_str, '(I7.7)') idump
    filename = './NBDUMP/nb.' // step_str // '.dump'

    ! Open file for writing
    open(newunit=ifdump_unit, file=trim(filename), status='replace', action='write', iostat=ios)
    if (ios /= 0) then
      print *, "Error: Cannot open dump file ", trim(filename)
      return
    end if

    ! Write LAMMPS dump format header
    ! write(ifdump_unit, '(A)') 'ITEM: TIMESTEP'
    ! write(ifdump_unit, '(I0)') istep

    write(ifdump_unit, '(A)') 'ITEM: NUMBER OF ATOMS'
    write(ifdump_unit, '(I0)') NB%n_bodies

    write(ifdump_unit, '(A)') 'ITEM: BOX BOUNDS pp pp pp'
    ! X bounds
    if (abs(L_bound(1)) < 1e-10_wp) then
      write(ifdump_unit, '(2F12.6)') minval(NB%pos(:,1)) - 1.0_wp, maxval(NB%pos(:,1)) + 1.0_wp
    else
      write(ifdump_unit, '(2F12.6)') -L_bound(1)*0.5_wp, L_bound(1)*0.5_wp
    end if
    ! Y bounds
    if (abs(L_bound(2)) < 1e-10_wp) then
      write(ifdump_unit, '(2F12.6)') minval(NB%pos(:,2)) - 1.0_wp, maxval(NB%pos(:,2)) + 1.0_wp
    else
      write(ifdump_unit, '(2F12.6)') -L_bound(2)*0.5_wp, L_bound(2)*0.5_wp
    end if
    ! Z bounds
    if (ndim == 3) then
      if (abs(L_bound(3)) < 1e-10_wp) then
        write(ifdump_unit, '(2F12.6)') minval(NB%pos(:,3)) - 1.0_wp, maxval(NB%pos(:,3)) + 1.0_wp
      else
        write(ifdump_unit, '(2F12.6)') -L_bound(3)*0.5_wp, L_bound(3)*0.5_wp
      end if
    else
      write(ifdump_unit, '(2F12.6)') 0.0_wp, 0.0_wp
    end if

    ! Write atom data header based on dimension
    if (ndim == 3) then
      write(ifdump_unit, '(A)') 'ITEM: ATOMS id type x y vx vy ax ay fx fy mass radius ke pe te'
    else if (ndim == 2) then
      write(ifdump_unit, '(A)') 'ITEM: ATOMS id type x y z vx vy vz ax ay az fx fy fz mass radius ke pe te'
    end if

    ! Write atom data
    do i = 1, NB%n_bodies
      if (ndim == 3) then
        write(ifdump_unit, '(I0,1X,I0,1X,17ES14.6)') i, 1, &
          NB%pos(i,1), NB%pos(i,2), NB%pos(i,3), &
          NB%vel(i,1), NB%vel(i,2), NB%vel(i,3), &
          NB%accel(i,1), NB%accel(i,2), NB%accel(i,3), &
          NB%force(i,1), NB%force(i,2), NB%force(i,3), &
          NB%mass(i), NB%radius(i), &
          NB%kinetic_energy(i), NB%potential_energy(i), NB%sum_energy(i)
      else if (ndim == 2) then
        write(ifdump_unit, '(I0,1X,I0,1X,14ES14.6)') i, 1, &
          NB%pos(i,1), NB%pos(i,2), 0.0_wp, &
          NB%vel(i,1), NB%vel(i,2), 0.0_wp, &
          NB%accel(i,1), NB%accel(i,2), 0.0_wp, &
          NB%force(i,1), NB%force(i,2), 0.0_wp, &
          NB%mass(i), NB%radius(i), &
          NB%kinetic_energy(i), NB%potential_energy(i), NB%sum_energy(i)
      end if
    end do
    ! Close file
    close(ifdump_unit)
    ! Create filename with 7-digit timestep number
    write(step_str, '(I7.7)') idump
    filename = './NBDUMPCSV/nb.' // step_str // '.csv'

    ! Open file for writing
    open(newunit=ifdump_unit, file=trim(filename), status='replace', action='write', iostat=ios)
    if (ios /= 0) then
      print *, "Error: Cannot open CSV dump file ", trim(filename)
      return
    end if

    ! Write CSV header based on dimension
    if (ndim == 3) then
      write(ifdump_unit, '(A)') 'id,type,x,y,z,vx,vy,vz,ax,ay,az,fx,fy,fz,mass,radius,ke,pe,te'
    else if (ndim == 2) then
      write(ifdump_unit, '(A)') 'id,type,x,y,vx,vy,ax,ay,fx,fy,mass,radius,ke,pe,te'
    else ! ndim == 1
      write(ifdump_unit, '(A)') 'id,type,x,vx,ax,fx,mass,radius,ke,pe,te'
    end if

    ! Write particle data
    do i = 1, NB%n_bodies
      if (ndim == 3) then
        write(ifdump_unit, '(I0,A,I0,17(A,ES14.6))') i, ',', 1, &
          ',', NB%pos(i,1), ',', NB%pos(i,2), ',', NB%pos(i,3), &
          ',', NB%vel(i,1), ',', NB%vel(i,2), ',', NB%vel(i,3), &
          ',', NB%accel(i,1), ',', NB%accel(i,2), ',', NB%accel(i,3), &
          ',', NB%force(i,1), ',', NB%force(i,2), ',', NB%force(i,3), &
          ',', NB%mass(i), ',', NB%radius(i), &
          ',', NB%kinetic_energy(i), ',', NB%potential_energy(i), ',', NB%sum_energy(i)
      else if (ndim == 2) then
        write(ifdump_unit, '(I0,A,I0,14(A,ES14.6))') i, ',', 1, &
          ',', NB%pos(i,1), ',', NB%pos(i,2), &
          ',', NB%vel(i,1), ',', NB%vel(i,2), &
          ',', NB%accel(i,1), ',', NB%accel(i,2), &
          ',', NB%force(i,1), ',', NB%force(i,2), &
          ',', NB%mass(i), ',', NB%radius(i), &
          ',', NB%kinetic_energy(i), ',', NB%potential_energy(i), ',', NB%sum_energy(i)
      else ! ndim == 1
        write(ifdump_unit, '(I0,A,I0,10(A,ES14.6))') i, ',', 1, &
          ',', NB%pos(i,1), &
          ',', NB%vel(i,1), &
          ',', NB%accel(i,1), &
          ',', NB%force(i,1), &
          ',', NB%mass(i), ',', NB%radius(i), &
          ',', NB%kinetic_energy(i), ',', NB%potential_energy(i), ',', NB%sum_energy(i)
      end if
    end do

    ! Close file
    close(ifdump_unit)
  end subroutine dump_data

  subroutine finalize_simulation(NB)
    implicit none
    type(NBodies_type), intent(inout) :: NB

    ! Output final state if needed
    print *, "Simulation completed successfully!"
    print *, "Final time: ", time
    print *, "Total energy: ", NB%total_energy
    print *, "Energy conservation: ", abs(NB%total_energy - NB%total_energy) ! You'd need initial energy stored

    ! Deallocate arrays
    call NB%deallocate()

  end subroutine finalize_simulation

  subroutine output_diagnostics(NB)
    implicit none
    type(NBodies_type), intent(in) :: NB

    print '(A,I8,A,F12.6)', "Step: ", istep, ", Time: ", time
    print '(A,F12.6)',  "  dt:                ", dt
    print '(A,ES15.8)', "  Total KE:          ", NB%total_kinetic_energy
    print '(A,ES15.8)', "  Total PE:          ", NB%total_potential_energy
    print '(A,ES15.8)', "  Total E:           ", NB%total_energy
    print '(A,ES15.8,A)', "  Delta E:           ", &
      ((NB%total_energy - NB%total_energy_initial) / ABS(NB%total_energy_initial)) * 100.0_wp," %"
    print '(A,ES15.8)', "  Virial F (2KE+PE): ", NB%total_kinetic_energy*2.0_wp + NB%total_potential_energy
    ! In the time loop:
    if(ndim == 2) then
      print '(A,2ES12.4)', "  COM position: ", NB%center_of_mass(1:ndim)
      print '(A,2ES12.4)', "  COM velocity: ", NB%center_of_mass_velocity(1:ndim)
      write(iftimehistory_unit, '(I8,1X,F12.6,1X,8ES16.8)') &
        istep, time, &
        NB%center_of_mass(1), NB%center_of_mass(2), &
        NB%center_of_mass_velocity(1), NB%center_of_mass_velocity(2), &
        NB%total_kinetic_energy, NB%total_potential_energy, NB%total_energy
    else if(ndim == 3) then
      print '(A,3ES12.4)', "  COM position: ", NB%center_of_mass(1:ndim)
      print '(A,3ES12.4)', "  COM velocity: ", NB%center_of_mass_velocity(1:ndim)
      write(iftimehistory_unit, '(I8,1X,F12.6,1X,9ES16.8)') &
        istep, time, &
        NB%center_of_mass(1), NB%center_of_mass(2), NB%center_of_mass(3), &
        NB%center_of_mass_velocity(1), NB%center_of_mass_velocity(2), NB%center_of_mass_velocity(3), &
        NB%total_kinetic_energy, NB%total_potential_energy, NB%total_energy
    end if

  end subroutine output_diagnostics
  subroutine compute_forces_serial_diagnosis(NB)
    ! Force_diagnose = 1: basic force computation with diagnosis
    ! and Force_matrix should be allocated
    implicit none
    type(NBodies_type), intent(inout) :: NB

    integer :: i, j
    real(kind=wp) :: pos_diff(ndim)
    real(kind=wp) :: dist, dist_cubed ! distance and distance cubed between bodies
    character(len=256) :: filename
    integer :: ios
    character(len=7) :: step_str
    logical :: dir_exists
    ! Initialize forces to zero
    NB%force(:, :) = 0.0_wp
    NB%potential_energy = 0.0_wp ! Reset potential energy array
    NB%total_potential_energy = 0.0_wp ! Reset total potential energy
    NB%kinetic_energy = 0.0_wp ! Reset kinetic energy array
    NB%total_kinetic_energy = 0.0_wp ! Reset total kinetic energy
    NB%total_energy = 0.0_wp ! Reset total energy
    NB%Force_matrix(:, :, :) = 0.0_wp ! Reset force matrix for diagnosis
    ! Compute pairwise gravitational forces
    do i = 1, NB%n_bodies
      do j = 1, NB%n_bodies
        if (i /= j) then
          ! force on particle i due to particle j
          ! i will be pulled towards j
          ! if x(j) > x(i) then fx(i) += positive value i.e. pull right towards
          ! j
          ! therefore x_diff = x(j) - x(i) is correct.
          ! alternatively we can do dx = x(i) - x(j) and then subtract the
          ! force.
          ! in vectorized form:
          pos_diff = NB%pos(j, 1:ndim) - NB%pos(i, 1:ndim)
          dist = sqrt(sum(pos_diff(1:ndim)**2)) + small_number ! Softening to avoid singularity
          dist_cubed = dist**3

          NB%force(i, 1:ndim) = NB%force(i, 1:ndim) &
            + Gravitational_Constant * NB%mass(i) * NB%mass(j) * (pos_diff(1:ndim) &
            / dist_cubed)

          NB%Force_matrix(i, j, 1:ndim) = Gravitational_Constant * NB%mass(i) &
            * NB%mass(j) * (pos_diff(1:ndim) / dist_cubed)
          ! Calculate potential energy for diagnostics
          ! Note: Using 0.5 factor to avoid double counting since we loop over
          ! all pairs
          NB%potential_energy(i) = NB%potential_energy(i) - &
            (Gravitational_Constant * NB%mass(i) * NB%mass(j) / dist) * 0.5_wp
        end if
      end do
    end do

    ! Calculate total potential energy
    NB%total_potential_energy = sum(NB%potential_energy)
    INQUIRE (FILE='./ForceMatrixCSV/', EXIST=dir_exists)
    IF (.NOT. dir_exists) THEN
      CALL EXECUTE_COMMAND_LINE('mkdir -p ./ForceMatrixCSV/')
    END IF

    do i=1, NB%n_bodies
      write(step_str, '(I7.7)') i
      filename = './ForceMatrixCSV/nb.' // step_str // '.csv'
      ! Open file for writing
      open(newunit=ifforcematrix_unit, file=trim(filename), status='replace', action='write', iostat=ios)
      if (ios /= 0) then
        print *, "Error: Cannot open Force Matrix CSV file ", trim(filename)
        return
      end if

      ! Write CSV header and data based on dimension
      if(ndim == 2) then
        write(ifforcematrix_unit, '(A)') 'from_id,dx,dy,mt,mf,fx,fy'
        do j=1, NB%n_bodies
          if (i /= j) then
            ! vector is: to -> from .
            ! ie. goes from point "to" (i) to point "from" (j)
            pos_diff = NB%pos(j, 1:ndim) - NB%pos(i, 1:ndim) ! vector going as: from - to
            write(ifforcematrix_unit, '(I0,6(",",ES14.6))') &
              j, pos_diff(1), pos_diff(2), NB%mass(j), NB%mass(i), &
              NB%Force_matrix(i,j,1), NB%Force_matrix(i,j,2)
          end if
        end do
      else if(ndim == 3) then
        write(ifforcematrix_unit, '(A)') 'from_id,dx,dy,dz,mt,mf,fx,fy,fz'
        do j=1, NB%n_bodies
          if (i /= j) then
            pos_diff = NB%pos(j, 1:ndim) - NB%pos(i, 1:ndim) ! vector going as: from - to
            write(ifforcematrix_unit, '(I0,8(",",ES14.6))') &
              j, pos_diff(1), pos_diff(2), pos_diff(3), NB%mass(j), NB%mass(i), &
              NB%Force_matrix(i,j,1), NB%Force_matrix(i,j,2), NB%Force_matrix(i,j,3)
          end if
        end do
      end if

      ! Close file
      close(ifforcematrix_unit)
    enddo

  end subroutine compute_forces_serial_diagnosis
  subroutine adjust_timestep_acceleration(NB)
    ! Based on :arXiv:2401.02849v1 [astro-ph.EP] 05 Jan 2024
    ! A new timestep criterion for N-body simulations - P, R, S
    implicit none
    type(NBodies_type), intent(inout) :: NB

    real(wp) :: tau_prs, tau_min, tau_local
    real(wp) :: accel_mag, snap_mag
    real(wp) :: snap(ndim), accel(ndim), accel_0(ndim), accel_00(ndim)
    real(wp) :: jerk(ndim), jerk_mag
    real(wp) :: dt_calculated, dt_close, min_distance, dist, relative_vel
    real(wp), parameter :: eta_aarseth = 0.02_wp  ! Aarseth accuracy parameter
    real(wp), parameter :: eta_close = 0.1_wp
    real(wp), parameter :: safety_factor = 0.5_wp
    real(wp), parameter :: eta_prs = (5.04*(10**(-6)))**(1.0_wp/7.0_wp) ! From the paper
    real(wp), parameter :: epsilon_prs = (eta_prs**7.0_wp)/5040.0_wp
    real(wp) :: numerator, denominator
    integer :: i, j

    ! Initialize
    tau_min = huge(1.0_wp)
    min_distance = huge(1.0_wp)
    dt_close = dt_max

    ! Calculate derivatives for Aarseth criterion
    !$omp parallel do default(shared) &
    !$omp private(i, accel_mag, snap_mag, jerk_mag, accel, accel_0, accel_00, &
    !$omp         jerk, snap, numerator, denominator, tau_prs) &
    !$omp reduction(min:tau_min)
    do i = 1, NB%n_bodies
      ! Compute acceleration magnitude
      accel_mag = sqrt(sum(NB%accel(i, 1:ndim)**2))

      ! Compute jerk and its magnitude
      accel = NB%accel(i, 1:ndim)
      accel_0 = NB%accel_0(i, 1:ndim)
      jerk = (accel - accel_0) / dt
      jerk_mag = sqrt(sum(jerk**2))

      ! Compute snap and its magnitude
      accel_00 = NB%accel_00(i, 1:ndim)
      snap = (accel - 2.0_wp * accel_0 + accel_00) / dt**2
      snap_mag = sqrt(sum(snap**2))

      ! Use the formula from the image to compute tau_prs
      if (snap_mag > epsilon(1.0_wp) .and. jerk_mag > epsilon(1.0_wp) .and. accel_mag > epsilon(1.0_wp)) then
        numerator = 2.0_wp * accel_mag**2
        denominator = jerk_mag**2 + accel_mag * snap_mag + small_number
        if (denominator > epsilon(1.0_wp)) then
          tau_prs = sqrt(2.0_wp * (numerator / denominator))  ! From formula (16) scaled with sqrt(2)

          ! Find minimum tau_prs across all bodies through OpenMP reduction
          tau_min = min(tau_min, tau_prs)
        end if
      end if
    end do
    !$omp end parallel do

    dt_calculated = (5040.0_wp*epsilon_prs)**(1.0_wp/7.0_wp) * tau_min

    ! Apply safety factor and limits
    ! dt_calculated = safety_factor * dt_calculated

    ! Apply absolute limits
    dt = min(max(dt_calculated, dt_min), dt_max)
    dt = dt_min*100.0_wp ! For testing, fix dt to a constant value
    dt_half = 0.5_wp * dt
    dt_squared = dt * dt
    
    ! Diagnostic output
    if (mod(istep, nprint) == 0) then
      print '(A,4ES12.5)', "  Timestep info: dt, dt_arseth, dt_min, dt_max = ", &
        dt, dt_calculated, dt_min, dt_max
    end if
  end subroutine adjust_timestep_acceleration
  subroutine normalize_mass(NB)
    implicit none
    type(NBodies_type), intent(inout) :: NB
    integer :: i, j
    real(kind=wp) :: pos_diff(ndim)
    real(kind=wp) :: dist
    real(kind=wp) :: v_squared_sum, pe_coefficient
    real(kind=wp) :: KE_total, PE_total

    v_squared_sum = 0.0_wp
    pe_coefficient = 0.0_wp

    ! Calculate sum of velocity squared (for all particles)
    do i = 1, NB%n_bodies
      v_squared_sum = v_squared_sum + sum(NB%vel(i, 1:ndim)**2)
    end do

    ! Calculate PE coefficient (sum of 1/r_ij)
    do i = 1, NB%n_bodies
      do j = i+1, NB%n_bodies
        pos_diff = NB%pos(j, 1:ndim) - NB%pos(i, 1:ndim)
        dist = sqrt(sum(pos_diff(1:ndim)**2) + softening_length_squared)
        pe_coefficient = pe_coefficient - 1.0_wp / dist  ! Negative because PE is negative
      end do
    end do

    ! For virial equilibrium: 2*KE + PE = 0, or KE = -PE/2
    ! Total KE = 0.5 * sum(m_i * v_i^2) = 0.5 * m * sum(v_i^2)  (for equal masses)
    ! Total PE = G * sum(m_i * m_j / r_ij) = G * m^2 * sum(1/r_ij)  (for equal masses)
    !
    ! So: 0.5 * m * v_squared_sum = -0.5 * G * m^2 * pe_coefficient
    ! Therefore: m = -v_squared_sum / (G * pe_coefficient)

    mass_0 = -v_squared_sum / (Gravitational_Constant * pe_coefficient)

    ! Set all masses to mass_0
    do i = 1, NB%n_bodies
      NB%mass(i) = mass_0
      NB%mass_inv(i) = 1.0_wp / NB%mass(i)
    end do

    NB%total_mass = NB%n_bodies * mass_0
    NB%total_mass_inv = 1.0_wp / NB%total_mass

    ! Verify virial equilibrium
    KE_total = 0.5_wp * mass_0 * v_squared_sum
    PE_total = Gravitational_Constant * mass_0**2 * pe_coefficient

    print *, "Mass per particle:", mass_0
    print *, "Total KE:", KE_total
    print *, "Total PE:", PE_total
    print *, "Virial ratio KE/|PE|:", KE_total/abs(PE_total)
    print *, "2*KE + PE:", 2.0_wp*KE_total + PE_total

  end subroutine normalize_mass
end module nbodies
