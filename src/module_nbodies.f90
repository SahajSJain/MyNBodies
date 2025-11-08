module nbodies
  use, intrinsic :: iso_fortran_env
  implicit none

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! PRECISION CONFIGURATION
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! #ifdef P16
!   ! Half precision - using single precision as fallback
!   ! Note: Fortran doesn't have native 16-bit float
!   integer, parameter :: wp = real32
!   character(len=*), parameter :: precision_name = "single (32-bit) - P16 fallback"
! #elif defined(P32)
  ! Single precision
  integer, parameter :: wp = real32
  character(len=*), parameter :: precision_name = "single (32-bit)"
! #elif defined(P64)
!   ! Double precision
!   integer, parameter :: wp = real64
!   character(len=*), parameter :: precision_name = "double (64-bit)"
! #elif defined(P128)
!   ! Quadruple precision
!   integer, parameter :: wp = real128
!   character(len=*), parameter :: precision_name = "quadruple (128-bit)"
! #else
!   ! Default to double precision
!   integer, parameter :: wp = real64
!   character(len=*), parameter :: precision_name = "double (64-bit) [default]"
! #endif

  ! Mathematical constants in working precision
  real(kind=wp), parameter :: pi = 4.0_wp*atan(1.0_wp)
  real(kind=wp), parameter :: e = exp(1.0_wp)
  real(kind=wp), parameter :: small_number = 2.0_wp*epsilon(1.0_wp) ! 2x machine small_number for working precision
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
  ! Mass initialization
  integer :: MASS_INITIALIZATION_TYPE ! 0=from file, 1=normalized
  ! Velocity initialization
  integer :: VELOCITY_INITIALIZATION_TYPE    ! 1=random, 2=solid body rotation
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
  real(kind=wp) :: total_ke, total_pe, total_COM(1:3), total_COM_vel(1:3) ! Tolerance for force diagnosis
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NBODIES TYPE DEFINITION
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  type :: NBodies_type
    ! Number of bodies
    integer :: n_bodies

    ! DYNAMICAL ARRAYS - BODY PROPERTIES
    ! Position and velocity arrays (n_bodies x ndim)
    real(kind=wp), allocatable :: pos(:, :)       ! Current positions
    real(kind=wp), allocatable :: pos_0(:, :)     ! Previous positions
    real(kind=wp), allocatable :: vel(:, :)       ! Current velocities
    real(kind=wp), allocatable :: vel_0(:, :)     ! Previous velocities
    real(kind=wp), allocatable :: vel_h(:, :)     ! Half-step velocities

    ! Force and acceleration arrays (n_bodies x ndim)
    real(kind=wp), allocatable :: accel(:, :)     ! Accelerations
    real(kind=wp), allocatable :: force(:, :)     ! Force vectors
    real(kind=wp), allocatable :: bc_factor(:, :) ! Boundary condition factors

    ! Mass and geometry arrays (n_bodies)
    real(kind=wp), allocatable :: mass(:)        ! Body masses
    real(kind=wp), allocatable :: mass_inv(:)    ! Inverse masses
    real(kind=wp), allocatable :: radius(:)      ! Body radii
    real(kind=wp), allocatable :: radius2(:)      ! Body radii squared
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
    real(kind=wp), allocatable :: Force_matrix(:, :, :)  ! matrix to store forces for diagnosis
    ! (to,from,ndim)
  contains
    procedure :: allocate => allocate_nbodies
    procedure :: deallocate => deallocate_nbodies
    procedure :: estimate_memory => estimate_memory_nbodies
    final :: finalize_nbodies
  end type NBodies_type
  interface
    module subroutine initialize_velocities(NB)
      type(NBodies_type), intent(inout) :: NB
    end subroutine initialize_velocities
    module subroutine initialize_positions(NB)
      type(NBodies_type), intent(inout) :: NB
    end subroutine initialize_positions
    module subroutine read_input_parameters(NB)
      type(NBodies_type), intent(inout) :: NB
    end subroutine read_input_parameters
    module subroutine dump_data(NB)
      type(NBodies_type), intent(inout) :: NB
    end subroutine dump_data
    module subroutine output_diagnostics(NB)
      type(NBodies_type), intent(in) :: NB
    end subroutine output_diagnostics
  end interface

contains
  subroutine estimate_memory_nbodies(this)
    class(NBodies_type), intent(in) :: this
    integer(kind=8) :: total_bytes
    integer :: n

    n = this%n_bodies
    total_bytes = 0

    ! Position and velocity arrays (n_bodies x ndim)
    total_bytes = total_bytes + 5*n*ndim*storage_size(1.0_wp)/8  ! pos, pos_0, vel, vel_0, vel_h

    ! Force and acceleration arrays (n_bodies x ndim)
    total_bytes = total_bytes + 3*n*ndim*storage_size(1.0_wp)/8  ! accel, force, bc_factor

    ! Mass and geometry arrays (n_bodies)
    total_bytes = total_bytes + 4*n*storage_size(1.0_wp)/8  ! mass, mass_inv, radius, radius2

    ! Energy arrays (n_bodies)
    total_bytes = total_bytes + 3*n*storage_size(1.0_wp)/8  ! kinetic_energy, potential_energy, sum_energy

    ! Center of mass arrays (ndim)
    total_bytes = total_bytes + 2*ndim*storage_size(1.0_wp)/8  ! center_of_mass, center_of_mass_velocity

    ! Random vector (ndim)
    total_bytes = total_bytes + ndim*storage_size(1.0_wp)/8  ! rand_vec

    ! Force matrix (n_bodies x n_bodies x ndim)
    total_bytes = total_bytes + n*n*ndim*storage_size(1.0_wp)/8  ! Force_matrix

    ! Print results
    print '(A)', '--- Memory Estimate for NBodies_type ---'
    print '(A,I0)', 'Number of bodies: ', n
    print '(A,I0)', 'Number of dimensions: ', ndim
    print '(A,I0,A)', 'Total memory required: ', total_bytes, ' bytes'
    print '(A,F0.2,A)', 'Total memory required: ', real(total_bytes)/1024.0, ' KB'
    print '(A,F0.2,A)', 'Total memory required: ', real(total_bytes)/1048576.0, ' MB'

  end subroutine estimate_memory_nbodies
  subroutine allocate_nbodies(this, num_bodies)
    implicit none
    class(NBodies_type), intent(inout) :: this
    integer, intent(in) :: num_bodies
    integer :: i
    real(kind=wp) :: rand_val

    this%n_bodies = num_bodies

    ! Validate ndim
    if (ndim /= 2 .and. ndim /= 3) then
      ndim = 2 ! default to 2D
    end if

    ! Allocate position and velocity arrays
    allocate (this%pos(this%n_bodies, ndim), source=0.0_wp)
    allocate (this%pos_0(this%n_bodies, ndim), source=0.0_wp)
    allocate (this%vel(this%n_bodies, ndim), source=0.0_wp)
    allocate (this%vel_0(this%n_bodies, ndim), source=0.0_wp)
    allocate (this%vel_h(this%n_bodies, ndim), source=0.0_wp)

    ! Allocate force and acceleration arrays
    allocate (this%accel(this%n_bodies, ndim), source=0.0_wp)
    allocate (this%force(this%n_bodies, 3), source=0.0_wp)
    allocate (this%bc_factor(this%n_bodies, ndim), source=0.0_wp)

    ! Allocate mass and geometry arrays
    allocate (this%mass(this%n_bodies), source=0.0_wp)
    allocate (this%mass_inv(this%n_bodies), source=0.0_wp)
    allocate (this%radius(this%n_bodies), source=0.0_wp)
    allocate (this%radius2(this%n_bodies), source=0.0_wp)

    ! Allocate energy arrays
    allocate (this%kinetic_energy(this%n_bodies), source=0.0_wp)
    allocate (this%potential_energy(this%n_bodies), source=0.0_wp)
    allocate (this%sum_energy(this%n_bodies), source=0.0_wp)

    ! Allocate center of mass arrays
    allocate (this%center_of_mass(3), source=0.0_wp)
    allocate (this%center_of_mass_velocity(3), source=0.0_wp)

    ! Allocate temporary arrays
    allocate (this%rand_vec(ndim), source=0.0_wp)

    ! Allocate force diagnosis matrix if needed
    if (Force_Diagnose == 1) then
      allocate (this%Force_matrix(this%n_bodies, this%n_bodies, ndim), source=0.0_wp)
    end if

    ! Initialize radii with variation
    radius_var = max(0.0_wp, min(radius_var, 0.9_wp)) ! clamp between 0 and 0.9
    do i = 1, this%n_bodies
      call random_number(rand_val)
      this%radius(i) = radius_0*((1.0_wp - radius_var) + rand_val*2.0_wp*radius_var)
      this%radius2(i) = this%radius(i)*this%radius(i) ! used in softening for force calculation
    end do

    ! Initialize masses based on radii (assuming constant density)
    do i = 1, this%n_bodies
      this%mass(i) = mass_0*(this%radius(i)/radius_0)**ndim
    end do

    ! Calculate inverse masses and total mass
    this%mass_inv = 1.0_wp/(this%mass + small_number)
    this%total_mass = sum(this%mass)
    this%total_mass_inv = 1.0_wp/(this%total_mass + small_number)

    ! Initialize global properties
    this%total_kinetic_energy = 0.0_wp
    this%total_potential_energy = 0.0_wp
    this%total_energy = 0.0_wp
    this%total_energy_initial = 0.0_wp
    this%min_vel = 0.0_wp
    this%max_vel = 0.0_wp
    this%min_accel = 0.0_wp
    this%max_accel = 0.0_wp
    this%char_dist = 0.0_wp

  end subroutine allocate_nbodies

  subroutine deallocate_nbodies(this)
    implicit none
    class(NBodies_type), intent(inout) :: this
    integer :: stat
    logical :: debug = .true.  ! Set to .true. for debugging

    ! Check if already deallocated by checking n_bodies
    if (this%n_bodies == 0) then
      if (debug) print *, "NBodies already deallocated, skipping"
      return
    end if

    if (debug) print *, "Starting deallocation for n_bodies =", this%n_bodies

    ! Deallocate position arrays
    if (allocated(this%pos)) then
      deallocate (this%pos, stat=stat)
      if (stat /= 0 .and. debug) print *, "Warning: Error deallocating pos, stat =", stat
    end if
    if (allocated(this%pos_0)) then
      deallocate (this%pos_0, stat=stat)
      if (stat /= 0 .and. debug) print *, "Warning: Error deallocating pos_0, stat =", stat
    end if
    ! Deallocate velocity arrays
    if (allocated(this%vel)) then
      deallocate (this%vel, stat=stat)
      if (stat /= 0 .and. debug) print *, "Warning: Error deallocating vel, stat =", stat
    end if
    if (allocated(this%vel_0)) then
      deallocate (this%vel_0, stat=stat)
      if (stat /= 0 .and. debug) print *, "Warning: Error deallocating vel_0, stat =", stat
    end if
    if (allocated(this%vel_h)) then
      deallocate (this%vel_h, stat=stat)
      if (stat /= 0 .and. debug) print *, "Warning: Error deallocating vel_h, stat =", stat
    end if

    ! Deallocate acceleration arrays
    if (allocated(this%accel)) then
      deallocate (this%accel, stat=stat)
      if (stat /= 0 .and. debug) print *, "Warning: Error deallocating accel, stat =", stat
    end if

    ! Deallocate force arrays
    if (allocated(this%force)) then
      deallocate (this%force, stat=stat)
      if (stat /= 0 .and. debug) print *, "Warning: Error deallocating force, stat =", stat
    end if
    if (allocated(this%bc_factor)) then
      deallocate (this%bc_factor, stat=stat)
      if (stat /= 0 .and. debug) print *, "Warning: Error deallocating bc_factor, stat =", stat
    end if

    ! Deallocate mass and geometry arrays
    if (allocated(this%mass)) then
      deallocate (this%mass, stat=stat)
      if (stat /= 0 .and. debug) print *, "Warning: Error deallocating mass, stat =", stat
    end if
    if (allocated(this%mass_inv)) then
      deallocate (this%mass_inv, stat=stat)
      if (stat /= 0 .and. debug) print *, "Warning: Error deallocating mass_inv, stat =", stat
    end if
    if (allocated(this%radius)) then
      deallocate (this%radius, stat=stat)
      if (stat /= 0 .and. debug) print *, "Warning: Error deallocating radius, stat =", stat
    end if
    if (allocated(this%radius2)) then
      deallocate (this%radius2, stat=stat)
      if (stat /= 0 .and. debug) print *, "Warning: Error deallocating radius2, stat =", stat
    end if
    ! Deallocate energy arrays
    if (allocated(this%kinetic_energy)) then
      deallocate (this%kinetic_energy, stat=stat)
      if (stat /= 0 .and. debug) print *, "Warning: Error deallocating kinetic_energy, stat =", stat
    end if
    if (allocated(this%potential_energy)) then
      deallocate (this%potential_energy, stat=stat)
      if (stat /= 0 .and. debug) print *, "Warning: Error deallocating potential_energy, stat =", stat
    end if
    if (allocated(this%sum_energy)) then
      deallocate (this%sum_energy, stat=stat)
      if (stat /= 0 .and. debug) print *, "Warning: Error deallocating sum_energy, stat =", stat
    end if

    ! Deallocate center of mass arrays
    if (allocated(this%center_of_mass)) then
      deallocate (this%center_of_mass, stat=stat)
      if (stat /= 0 .and. debug) print *, "Warning: Error deallocating center_of_mass, stat =", stat
    end if
    if (allocated(this%center_of_mass_velocity)) then
      deallocate (this%center_of_mass_velocity, stat=stat)
      if (stat /= 0 .and. debug) print *, "Warning: Error deallocating center_of_mass_velocity, stat =", stat
    end if

    ! Deallocate temporary arrays
    if (allocated(this%rand_vec)) then
      deallocate (this%rand_vec, stat=stat)
      if (stat /= 0 .and. debug) print *, "Warning: Error deallocating rand_vec, stat =", stat
    end if
    if (allocated(this%Force_matrix)) then
      deallocate (this%Force_matrix, stat=stat)
      if (stat /= 0 .and. debug) print *, "Warning: Error deallocating Force_matrix, stat =", stat
    end if

    ! Reset all scalar values to indicate deallocation
    this%n_bodies = 0
    this%total_kinetic_energy = 0.0_wp
    this%total_potential_energy = 0.0_wp
    this%total_energy = 0.0_wp
    this%total_mass = 0.0_wp
    this%total_mass_inv = 0.0_wp
    this%total_energy_initial = 0.0_wp
    this%min_vel = 0.0_wp
    this%max_vel = 0.0_wp
    this%min_accel = 0.0_wp
    this%max_accel = 0.0_wp
    this%char_dist = 0.0_wp

    if (debug) print *, "Deallocation complete"

  end subroutine deallocate_nbodies

  subroutine finalize_nbodies(this)
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
    call set_mass(NB) ! normalize masses only if not read in from file
    ! Create output directory if it doesn't exist
    INQUIRE (FILE='./NBDUMP/', EXIST=dir_exists)
    IF (.NOT. dir_exists) THEN
      CALL EXECUTE_COMMAND_LINE('mkdir -p ./NBDUMP/')
    END IF
    INQUIRE (FILE='./NBDUMPCSV/', EXIST=dir_exists)
    IF (.NOT. dir_exists) THEN
      CALL EXECUTE_COMMAND_LINE('mkdir -p ./NBDUMPCSV/')
    END IF
    open (newunit=iftimehistory_unit, file='Time_History.dat', status='replace', action='write', iostat=ios)
    if (ios /= 0) then
      print *, "Error: Cannot open input file ", trim('Time_History.dat')
      stop
    end if
    ! Write header for time history file (do this once before the time loop)
    if (ndim == 2) then
      write (iftimehistory_unit, '(A)') 'step, time, COM_x, COM_y, '// &
        'COM_vx, COM_vy, KE_total, PE_total, E_total'
    else if (ndim == 3) then
      write (iftimehistory_unit, '(A)') 'step, time, COM_x, COM_y, COM_z, '// &
        'COM_vx, COM_vy, COM_vz, KE_total, PE_total, E_total'
    end if

    ! Initialize half-step velocities for leapfrog integration
    NB%vel_h = NB%vel
    NB%pos_0 = NB%pos
    NB%vel_0 = NB%vel
    ! Initialize time variables
    time = 0.0_wp
    istep = 0
    dt_half = 0.5_wp*dt
    dt_squared = dt*dt
    dt_min = dt*0.01_wp
    dt_max = dt*100.0_wp
    ! Print summary of parameters read
    print *, "=== Input Parameters Read Successfully ==="
    print '(A,I3,A,I6)', "Dimensions: ", ndim, ", Bodies: ", NB%n_bodies
    print '(A,F8.3)', "Gravitational constant: ", Gravitational_Constant
    print '(A,F15.8,A,F8.3,A,F8.3,A,F8.3)', "Mass: ", mass_0, ", Radius: ", radius_0, &
      ", Radius variation: ", radius_var, "radius^2_0 used in softening: ", radius_0*radius_0
    print '(A,F8.6,A,I8)', "Time step: ", dt, ", Number of steps: ", nsteps
    print '(A,I2)', "Initialization type: ", INITIALIZATION_TYPE
    print '(A,I2)', "Velocity initialization type: ", VELOCITY_INITIALIZATION_TYPE
    print '(A,I2)', "Boundary condition type: ", BOUNDARY_CONDITION_TYPE
    print *, "========================================="

    print *, "=== Simulation Initialized Successfully ==="
    print '(A)', " Precision: "//trim(precision_name)
    print '(A,I8)', " Number of bodies: ", NB%n_bodies
    print '(A,I2)', " Dimensions: ", ndim
    print '(A,F12.6)', " Time step (dt): ", dt
    print '(A,F12.6)', " Half time step (dt/2): ", dt_half
    print '(A,F12.6)', " Minimum time step (dt_min): ", dt_min
    print '(A,F12.6)', " Maximum time step (dt_max): ", dt_max
    print '(A,I8)', " Total steps: ", nsteps
    print '(A,I8)', " Dump frequency: ", ndump
    print *, "==========================================="
    call estimate_memory_nbodies(NB)
  end subroutine initialize_simulation

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
          dist = sqrt(sum(pos_diff(1:ndim)**2)) + NB%radius2(i) + NB%radius2(j) ! Softening to avoid singularity
          dist_cubed = dist**3

          NB%force(i, 1:ndim) = NB%force(i, 1:ndim) &
                                + Gravitational_Constant*NB%mass(i)*NB%mass(j)*(pos_diff(1:ndim) &
                                                                                /dist_cubed)

          NB%Force_matrix(i, j, 1:ndim) = Gravitational_Constant*NB%mass(i) &
                                          *NB%mass(j)*(pos_diff(1:ndim)/dist_cubed)
          ! Calculate potential energy for diagnostics
          ! Note: Using 0.5 factor to avoid double counting since we loop over
          ! all pairs
          NB%potential_energy(i) = NB%potential_energy(i) - &
                                   (Gravitational_Constant*NB%mass(i)*NB%mass(j)/dist)*0.5_wp
        end if
      end do
    end do

    ! Calculate total potential energy
    NB%total_potential_energy = sum(NB%potential_energy)
    INQUIRE (FILE='./ForceMatrixCSV/', EXIST=dir_exists)
    IF (.NOT. dir_exists) THEN
      CALL EXECUTE_COMMAND_LINE('mkdir -p ./ForceMatrixCSV/')
    END IF

    do i = 1, NB%n_bodies
      write (step_str, '(I7.7)') i
      filename = './ForceMatrixCSV/nb.'//step_str//'.csv'
      ! Open file for writing
      open (newunit=ifforcematrix_unit, file=trim(filename), status='replace', action='write', iostat=ios)
      if (ios /= 0) then
        print *, "Error: Cannot open Force Matrix CSV file ", trim(filename)
        return
      end if

      ! Write CSV header and data based on dimension
      if (ndim == 2) then
        write (ifforcematrix_unit, '(A)') 'from_id,dx,dy,mt,mf,fx,fy'
        do j = 1, NB%n_bodies
          if (i /= j) then
            ! vector is: to -> from .
            ! ie. goes from point "to" (i) to point "from" (j)
            pos_diff = NB%pos(j, 1:ndim) - NB%pos(i, 1:ndim) ! vector going as: from - to
            write (ifforcematrix_unit, '(I0,6(",",ES14.6))') &
              j, pos_diff(1), pos_diff(2), NB%mass(j), NB%mass(i), &
              NB%Force_matrix(i, j, 1), NB%Force_matrix(i, j, 2)
          end if
        end do
      else if (ndim == 3) then
        write (ifforcematrix_unit, '(A)') 'from_id,dx,dy,dz,mt,mf,fx,fy,fz'
        do j = 1, NB%n_bodies
          if (i /= j) then
            pos_diff = NB%pos(j, 1:ndim) - NB%pos(i, 1:ndim) ! vector going as: from - to
            write (ifforcematrix_unit, '(I0,8(",",ES14.6))') &
              j, pos_diff(1), pos_diff(2), pos_diff(3), NB%mass(j), NB%mass(i), &
              NB%Force_matrix(i, j, 1), NB%Force_matrix(i, j, 2), NB%Force_matrix(i, j, 3)
          end if
        end do
      end if

      ! Close file
      close (ifforcematrix_unit)
    end do

  end subroutine compute_forces_serial_diagnosis

  subroutine set_mass(NB)
    implicit none
    type(NBodies_type), intent(inout) :: NB
    integer :: i, j
    real(kind=wp) :: pos_diff(ndim)
    real(kind=wp) :: dist
    real(kind=wp) :: v_squared_sum, pe_coefficient
    real(kind=wp) :: KE_total, PE_total
    ! Local variables
    integer :: ios
    integer :: n_bodies_from_file
    integer :: file_unit
    character(len=256) :: filename
    integer :: step_size, sample_index

    if (MASS_INITIALIZATION_TYPE == 0) then
      ! Masses are read from file, do not modify
      ! Also read radius from file
      ! read masses from file
      filename = 'data/masses_ascii.dat'
      open (newunit=file_unit, file=trim(filename), status='old', action='read', iostat=ios)
      do i = 1, NB%n_bodies
        read (file_unit, *) NB%mass(i)
        NB%mass_inv(i) = 1.0_wp/NB%mass(i)
      end do
      NB%total_mass = sum(NB%mass(1:NB%n_bodies))
      NB%total_mass_inv = 1.0_wp/NB%total_mass
      close (file_unit)

      filename = 'data/radii_ascii.dat'
      open (newunit=file_unit, file=trim(filename), status='old', action='read', iostat=ios)
      do i = 1, NB%n_bodies
        read (file_unit, *) NB%radius(i)
        NB%radius2(i) = NB%radius(i)*NB%radius(i) ! used in softening for force calculation
      end do
      close (file_unit)

      ! Print sample masses and radii
      write (*, *) ''
      write (*, *) 'Sample masses and radii (10 equally spaced particles):'
      write (*, *) '--------------------------------------------------------'
      write (*, '(A10,A15,A15)') 'Index', 'Mass', 'Radius'
      write (*, *) '--------------------------------------------------------'

      ! Calculate step size to get 10 samples
      step_size = max(1, NB%n_bodies/10)

      do i = 1, 10
        sample_index = min(i*step_size, NB%n_bodies)
        write (*, '(I10,2E15.6)') sample_index, NB%mass(sample_index), NB%radius(sample_index)
      end do

      write (*, *) '--------------------------------------------------------'
      write (*, '(A,E15.6)') 'Total mass: ', NB%total_mass
      write (*, '(A,I10)') 'Total particles: ', NB%n_bodies
      write (*, *) ''

      return
    end if
    v_squared_sum = 0.0_wp
    pe_coefficient = 0.0_wp

    ! Calculate sum of velocity squared (for all particles)
    do i = 1, NB%n_bodies
      v_squared_sum = v_squared_sum + sum(NB%vel(i, 1:ndim)**2)
    end do

    ! Calculate PE coefficient (sum of 1/r_ij)
    do i = 1, NB%n_bodies
      do j = i + 1, NB%n_bodies ! avoid double counting
        pos_diff = NB%pos(j, 1:ndim) - NB%pos(i, 1:ndim)
        dist = sqrt(sum(pos_diff(1:ndim)**2) + NB%radius(i)*NB%radius(j))
        pe_coefficient = pe_coefficient - 1.0_wp/dist  ! Negative because PE is negative
      end do
    end do

    ! For virial equilibrium: 2*KE + PE = 0, or KE = -PE/2
    ! Total KE = 0.5 * sum(m_i * v_i^2) = 0.5 * m * sum(v_i^2)  (for equal masses)
    ! Total PE = G * sum(m_i * m_j / r_ij) = G * m^2 * sum(1/r_ij)  (for equal masses)
    !
    ! So: 0.5 * m * v_squared_sum = -0.5 * G * m^2 * pe_coefficient
    ! Therefore: m = -v_squared_sum / (G * pe_coefficient)

    mass_0 = -v_squared_sum/(Gravitational_Constant*pe_coefficient)

    ! Set all masses to mass_0
    do i = 1, NB%n_bodies
      NB%mass(i) = mass_0
      NB%mass_inv(i) = 1.0_wp/NB%mass(i)
    end do

    NB%total_mass = NB%n_bodies*mass_0
    NB%total_mass_inv = 1.0_wp/NB%total_mass

    ! Verify virial equilibrium
    KE_total = 0.5_wp*mass_0*v_squared_sum
    PE_total = Gravitational_Constant*mass_0**2*pe_coefficient

    print *, "Mass per particle:", mass_0
    print *, "Total KE:", KE_total
    print *, "Total PE:", PE_total
    print *, "Virial ratio KE/|PE|:", KE_total/abs(PE_total)
    print *, "2*KE + PE:", 2.0_wp*KE_total + PE_total

  end subroutine set_mass
end module nbodies

! subroutine adjust_timestep_acceleration(NB)
!   ! Based on :arXiv:2401.02849v1 [astro-ph.EP] 05 Jan 2024
!   ! A new timestep criterion for N-body simulations - P, R, S
!   implicit none
!   type(NBodies_type), intent(inout) :: NB

!   real(wp) :: tau_prs, tau_min, tau_local
!   real(wp) :: accel_mag, snap_mag
!   real(wp) :: snap(ndim), accel(ndim), accel_0(ndim), accel_00(ndim)
!   real(wp) :: jerk(ndim), jerk_mag
!   real(wp) :: dt_calculated, dt_close, min_distance, dist, relative_vel
!   real(wp), parameter :: eta_aarseth = 0.02_wp  ! Aarseth accuracy parameter
!   real(wp), parameter :: eta_close = 0.1_wp
!   real(wp), parameter :: safety_factor = 0.5_wp
!   real(wp), parameter :: eta_prs = (5.04*(10**(-6)))**(1.0_wp/7.0_wp) ! From the paper
!   real(wp), parameter :: epsilon_prs = (eta_prs**7.0_wp)/5040.0_wp
!   real(wp) :: numerator, denominator
!   integer :: i, j

!   ! Initialize
!   tau_min = huge(1.0_wp)
!   min_distance = huge(1.0_wp)
!   dt_close = dt_max

!   ! Calculate derivatives for Aarseth criterion
!   !$omp parallel do default(shared) &
!   !$omp private(i, accel_mag, snap_mag, jerk_mag, accel, accel_0, accel_00, &
!   !$omp         jerk, snap, numerator, denominator, tau_prs) &
!   !$omp reduction(min:tau_min)
!   do i = 1, NB%n_bodies
!     ! Compute acceleration magnitude
!     accel_mag = sqrt(sum(NB%accel(i, 1:ndim)**2))

!     ! Compute jerk and its magnitude
!     accel = NB%accel(i, 1:ndim)
!     accel_0 = NB%accel_0(i, 1:ndim)
!     jerk = (accel - accel_0) / dt
!     jerk_mag = sqrt(sum(jerk**2))

!     ! Compute snap and its magnitude
!     accel_00 = NB%accel_00(i, 1:ndim)
!     snap = (accel - 2.0_wp * accel_0 + accel_00) / dt**2
!     snap_mag = sqrt(sum(snap**2))

!     ! Use the formula from the image to compute tau_prs
!     if (snap_mag > epsilon(1.0_wp) .and. jerk_mag > epsilon(1.0_wp) .and. accel_mag > epsilon(1.0_wp)) then
!       numerator = 2.0_wp * accel_mag**2
!       denominator = jerk_mag**2 + accel_mag * snap_mag + small_number
!       if (denominator > epsilon(1.0_wp)) then
!         tau_prs = sqrt(2.0_wp * (numerator / denominator))  ! From formula (16) scaled with sqrt(2)

!         ! Find minimum tau_prs across all bodies through OpenMP reduction
!         tau_min = min(tau_min, tau_prs)
!       end if
!     end if
!   end do
!   !$omp end parallel do

!   dt_calculated = (5040.0_wp*epsilon_prs)**(1.0_wp/7.0_wp) * tau_min

!   ! Apply safety factor and limits
!   ! dt_calculated = safety_factor * dt_calculated

!   ! Apply absolute limits
!   dt = min(max(dt_calculated, dt_min), dt_max)
!   dt = dt_min*100.0_wp ! For testing, fix dt to a constant value
!   dt_half = 0.5_wp * dt
!   dt_squared = dt * dt

!   ! Diagnostic output
!   if (mod(istep, nprint) == 0) then
!     print '(A,4ES12.5)', "  Timestep info: dt, dt_arseth, dt_min, dt_max = ", &
!       dt, dt_calculated, dt_min, dt_max
!   end if
! end subroutine adjust_timestep_acceleration
