submodule(nbodies) init_pos
  implicit none
contains
  module subroutine initialize_positions(NB)
    implicit none
    type(NBodies_type), intent(inout) :: NB

    ! Initialize positions based on INITIALIZATION_TYPE
    select case (INITIALIZATION_TYPE)
    case (0)
      call read_positions_from_file(NB)
    case (1)
      call random_placement(NB)
    case (2)
      call random_spherical_placement(NB)
    case (3)
      call random_cylindrical_placement(NB)
    case (4)
      call ordered_placement(NB)
    case (5)
      call spherical_placement(NB)
    case (6)
      call cylindrical_placement(NB)
    case default
      print *, "Error: Unknown INITIALIZATION_TYPE ", INITIALIZATION_TYPE
      stop
    end select
  end subroutine initialize_positions

  subroutine read_positions_from_file(NB)
    implicit none
    type(NBodies_type), intent(inout) :: NB

    ! Local variables
    integer :: i, ios
    integer :: n_bodies_from_file
    integer :: file_unit
    character(len=256) :: filename
    filename = 'data/nstars.dat'

    open (newunit=file_unit, file=trim(filename), status='old', action='read', iostat=ios)
    read (file_unit, *) n_bodies_from_file
    read (file_unit, *) mass_0
    close (file_unit)

    ! Check if NB is already allocated with different size
    if (allocated(NB%pos)) then
      if (NB%n_bodies /= n_bodies_from_file) then
        call NB%deallocate()
        call NB%allocate(n_bodies_from_file)
      end if
    else
      call NB%allocate(n_bodies_from_file)
    end if
    ! Set filename
    filename = 'data/positions_ascii.dat'

    ! Open position file
    open (newunit=file_unit, file=trim(filename), status='old', action='read', iostat=ios)

    ! Read positions
    do i = 1, NB%n_bodies
      read (file_unit, *) NB%pos(i, 1), NB%pos(i, 2), NB%pos(i, 3)
    end do

    ! Close file
    close (file_unit)
    ! read masses from file
    ! this has been moved to set_mass subroutine
    ! Initialize old positions
    NB%pos_0 = NB%pos

  end subroutine read_positions_from_file

  subroutine ordered_placement(NB)
    implicit none
    type(NBodies_type), intent(inout) :: NB

    ! Local variables
    integer :: i, j, k, idx
    integer :: n_bodies_calc
    real(wp) :: dx, dy, dz

    ! Ensure ndim is valid
    if (ndim < 1) ndim = 1 ! minimum dimension is 1D
    if (ndim > 3) ndim = 3 ! maximum dimension is 3D

    ! Set unused dimensions to 1 for proper calculation
    if (ndim < 3) nbx_init(3) = 1
    if (ndim < 2) nbx_init(2) = 1

    ! Calculate total number of bodies
    n_bodies_calc = nbx_init(1)*nbx_init(2)*nbx_init(3)

    ! Check if NB is already allocated with different size
    if (allocated(NB%pos)) then
      if (NB%n_bodies /= n_bodies_calc) then
        call NB%deallocate()
        call NB%allocate(n_bodies_calc)
      end if
    else
      call NB%allocate(n_bodies_calc)
    end if

    ! Calculate grid spacings
    dx = L_init(1)/real(nbx_init(1), wp)
    dy = L_init(2)/real(nbx_init(2), wp)
    dz = L_init(3)/real(nbx_init(3), wp)

    ! Fill positions array
    idx = 0
    do k = 0, nbx_init(3) - 1
      do j = 0, nbx_init(2) - 1
        do i = 0, nbx_init(1) - 1
          idx = idx + 1

          ! Place particles at cell centers, then shift to [-L/2, L/2]
          NB%pos(idx, 1) = (i + 0.5_wp)*dx - L_init(1)/2.0_wp

          if (ndim >= 2) then
            NB%pos(idx, 2) = (j + 0.5_wp)*dy - L_init(2)/2.0_wp
          end if
          if (ndim >= 3) then
            NB%pos(idx, 3) = (k + 0.5_wp)*dz - L_init(3)/2.0_wp
          end if

        end do
      end do
    end do

    ! Initialize old positions
    NB%pos_0 = NB%pos

  end subroutine ordered_placement

  subroutine random_placement(NB)
    implicit none
    type(NBodies_type), intent(inout) :: NB

    integer :: i, j, iter
    real(wp) :: min_dist, dist
    logical :: too_close

    if (ndim < 1) ndim = 1
    if (ndim > 3) ndim = 3

    ! Check if NB needs allocation
    if (.not. allocated(NB%pos)) then
      call NB%allocate(NB%n_bodies)
    end if

    ! Minimum distance between particles
    min_dist = L_init(1)*0.0001

    ! Place particles with overlap checking
    do i = 1, NB%n_bodies
      iter = 0
      too_close = .true.

      do while (too_close .and. iter < 10)
        iter = iter + 1

        ! Generate random position
        call random_number(NB%rand_vec(1:ndim))
        NB%pos(i, 1:ndim) = (NB%rand_vec(1:ndim) - 0.5_wp)*L_init(1:ndim)

        ! Set unused dimensions to 0
        if (ndim < 3) NB%pos(i, 3) = 0.0_wp
        if (ndim < 2) NB%pos(i, 2) = 0.0_wp

        ! Check distance to all previously placed particles
        too_close = .false.
        do j = 1, i - 1
          dist = sqrt(sum((NB%pos(i, 1:3) - NB%pos(j, 1:3))**2))
          if (dist < min_dist) then
            too_close = .true.
            exit
          end if
        end do
      end do

      ! Warning if couldn't find non-overlapping position
      if (too_close) then
        print '(A,I0,A)', "Warning: Particle ", i, " may be too close to others after 10 attempts"
      end if
    end do

    ! Initialize old positions
    NB%pos_0 = NB%pos
  end subroutine random_placement

  subroutine random_spherical_placement(NB)
    implicit none
    type(NBodies_type), intent(inout) :: NB

    integer :: i, j, iter
    real(wp) :: min_dist, dist
    real(wp) :: r, theta, phi
    real(wp) :: radius_max
    logical :: too_close

    ! Check if NB needs allocation
    if (.not. allocated(NB%pos)) then
      call NB%allocate(NB%n_bodies)
    end if

    ! Maximum radius (half of smallest box dimension)
    radius_max = 0.5_wp*minval(L_init(1:3))

    ! Minimum distance between particles
    min_dist = L_init(1)*0.0001

    ! Place particles in sphere with overlap checking
    do i = 1, NB%n_bodies
      iter = 0
      too_close = .true.

      do while (too_close .and. iter < 10)
        iter = iter + 1

        ! Generate random spherical coordinates
        call random_number(NB%rand_vec(1:3))

        ! Random radius with uniform volume distribution
        r = radius_max*(NB%rand_vec(1))**(1.0_wp/3.0_wp)

        ! Random angles
        theta = acos(2.0_wp*NB%rand_vec(2) - 1.0_wp)  ! polar angle [0, pi]
        phi = 2.0_wp*pi*NB%rand_vec(3)               ! azimuthal angle [0, 2pi]

        ! Convert to Cartesian coordinates
        NB%pos(i, 1) = r*sin(theta)*cos(phi)
        NB%pos(i, 2) = r*sin(theta)*sin(phi)
        NB%pos(i, 3) = r*cos(theta)

        ! Check distance to all previously placed particles
        too_close = .false.
        do j = 1, i - 1
          dist = sqrt(sum((NB%pos(i, 1:3) - NB%pos(j, 1:3))**2))
          if (dist < min_dist) then
            too_close = .true.
            exit
          end if
        end do
      end do

      ! Warning if couldn't find non-overlapping position
      if (too_close) then
        print '(A,I0,A)', "Warning: Particle ", i, " may be too close to others after 10 attempts"
      end if
    end do

    ! Initialize old positions
    NB%pos_0 = NB%pos
  end subroutine random_spherical_placement

  subroutine random_cylindrical_placement(NB)
    implicit none
    type(NBodies_type), intent(inout) :: NB

    integer :: i, j, iter
    real(wp) :: min_dist, dist
    real(wp) :: r, theta, z
    real(wp) :: radius_max, height
    logical :: too_close

    ! Check if NB needs allocation
    if (.not. allocated(NB%pos)) then
      call NB%allocate(NB%n_bodies)
    end if

    ! Cylinder parameters (aligned along z-axis)
    radius_max = 0.5_wp*min(L_init(1), L_init(2))
    height = L_init(3)

    ! Minimum distance between particles
    min_dist = L_init(1)*0.0001

    ! Place particles in cylinder with overlap checking
    do i = 1, NB%n_bodies
      iter = 0
      too_close = .true.

      do while (too_close .and. iter < 10)
        iter = iter + 1

        ! Generate random cylindrical coordinates
        call random_number(NB%rand_vec(1:3))

        ! Random radius with uniform area distribution
        r = radius_max*sqrt(NB%rand_vec(1))

        ! Random angle
        theta = 2.0_wp*pi*NB%rand_vec(2)

        ! Random height
        z = (NB%rand_vec(3) - 0.5_wp)*height

        ! Convert to Cartesian coordinates
        NB%pos(i, 1) = r*cos(theta)
        NB%pos(i, 2) = r*sin(theta)
        NB%pos(i, 3) = z

        ! Check distance to all previously placed particles
        too_close = .false.
        do j = 1, i - 1
          dist = sqrt(sum((NB%pos(i, 1:3) - NB%pos(j, 1:3))**2))
          if (dist < min_dist) then
            too_close = .true.
            exit
          end if
        end do
      end do

      ! Warning if couldn't find non-overlapping position
      if (too_close) then
        print '(A,I0,A)', "Warning: Particle ", i, " may be too close to others after 10 attempts"
      end if
    end do

    ! Initialize old positions
    NB%pos_0 = NB%pos
  end subroutine random_cylindrical_placement

  subroutine spherical_placement(NB)
    implicit none
    type(NBodies_type), intent(inout) :: NB

    ! generated by Claude AI
    ! Spherical/ellipsoidal coordinate placement of bodies (3D only)
    integer :: i, j, k, count
    real(kind=wp) :: r, theta, phi
    real(kind=wp) :: dr, dtheta, dphi
    real(kind=wp) :: a, b, c ! semi-axes of ellipsoid
    integer :: n_bodies_calc
    ! Ensure ndim is 3
    if (ndim /= 3) then
      print *, "Error: Spherical placement requires ndim = 3"
      stop
    end if
    ! Set up grid parameters
    n_bodies_calc = nbx_init(1)*nbx_init(2)*nbx_init(3) ! r, theta, phi
    ! Check if NB is already allocated with different size
    if (allocated(NB%pos)) then
      if (NB%n_bodies /= n_bodies_calc) then
        call NB%deallocate()
        call NB%allocate(n_bodies_calc)
      end if
    else
      call NB%allocate(n_bodies_calc)
    end if
    ! Semi-axes of ellipsoid
    a = L_init(1)/2.0_wp
    b = L_init(2)/2.0_wp
    c = L_init(3)/2.0_wp
    count = 0
    ! 3D ellipsoidal coordinates (r, theta, phi)
    dr = 1.0_wp/real(nbx_init(1), wp) ! r goes from 0 to 1
    dtheta = pi/real(nbx_init(2), wp) ! polar angle [0, pi]
    dphi = 2.0_wp*pi/real(nbx_init(3), wp) ! azimuthal angle [0, 2pi]

    do i = 0, nbx_init(1) - 1
      do j = 0, nbx_init(2) - 1
        do k = 0, nbx_init(3) - 1
          count = count + 1
          r = (i + 0.5_wp)*dr ! This ensures no particle at r=0
          theta = (j + 0.5_wp)*dtheta ! This ensures no particle at theta=0 or pi
          phi = (k + 0.5_wp)*dphi ! This ensures no particle at phi=0 or 2pi
          ! Ellipsoid parametric equations
          NB%pos(count, 1) = r*a*sin(theta)*cos(phi)
          NB%pos(count, 2) = r*b*sin(theta)*sin(phi)
          NB%pos(count, 3) = r*c*cos(theta)
        end do
      end do
    end do
    ! Initialize old positions
    NB%pos_0 = NB%pos
  end subroutine spherical_placement

  subroutine cylindrical_placement(NB)
    implicit none
    type(NBodies_type), intent(inout) :: NB
    ! generated by Claude AI
    ! Cylindrical coordinate placement of bodies (2D/3D)
    integer :: i, j, k, count
    real(kind=wp) :: r, theta, z
    real(kind=wp) :: dr, dtheta, dz
    real(kind=wp) :: a, b ! semi-axes of ellipse
    integer :: n_bodies_calc
    ! Ensure ndim is valid
    if (ndim < 2) then
      print *, "Error: Cylindrical placement requires ndim >= 2"
      stop
    end if
    if (ndim > 3) ndim = 3
    ! Set up grid parameters
    if (ndim == 2) then
      n_bodies_calc = nbx_init(1)*nbx_init(2) ! r and theta
    else ! ndim == 3
      n_bodies_calc = nbx_init(1)*nbx_init(2)*nbx_init(3) ! r, theta, z
    end if
    ! Check if NB is already allocated with different size
    if (allocated(NB%pos)) then
      if (NB%n_bodies /= n_bodies_calc) then
        call NB%deallocate()
        call NB%allocate(n_bodies_calc)
      end if
    else
      call NB%allocate(n_bodies_calc)
    end if
    ! Semi-axes of ellipse
    a = L_init(1)/2.0_wp
    b = L_init(2)/2.0_wp
    count = 0
    if (ndim == 2) then
      ! 2D elliptical coordinates (r, theta)
      dr = 1.0_wp/real(nbx_init(1), wp) ! r goes from 0 to 1
      dtheta = 2.0_wp*pi/real(nbx_init(2), wp)
      do i = 0, nbx_init(1) - 1
        do j = 0, nbx_init(2) - 1
          count = count + 1
          r = (i + 0.5_wp)*dr ! This ensures no particle at r=0
          theta = (j + 0.5_wp)*dtheta ! This ensures no particle at theta=0 or 2pi
          ! Ellipse parametric equations
          NB%pos(count, 1) = r*a*cos(theta)
          NB%pos(count, 2) = r*b*sin(theta)
        end do
      end do
    else ! ndim == 3
      ! 3D cylindrical coordinates (r, theta, z)
      dr = 1.0_wp/real(nbx_init(1), wp) ! r goes from 0 to 1
      dtheta = 2.0_wp*pi/real(nbx_init(2), wp)
      dz = L_init(3)/real(nbx_init(3), wp)
      do i = 0, nbx_init(1) - 1
        do j = 0, nbx_init(2) - 1
          do k = 0, nbx_init(3) - 1
            count = count + 1
            r = (i + 0.5_wp)*dr ! This ensures no particle at r=0
            theta = (j + 0.5_wp)*dtheta ! This ensures no particle at theta=0 or 2pi
            z = (k + 0.5_wp)*dz - L_init(3)/2.0_wp
            ! Elliptical cylinder parametric equations
            NB%pos(count, 1) = r*a*cos(theta)
            NB%pos(count, 2) = r*b*sin(theta)
            NB%pos(count, 3) = z
          end do
        end do
      end do
    end if
    ! Initialize old positions
    NB%pos_0 = NB%pos
  end subroutine cylindrical_placement

end submodule init_pos
