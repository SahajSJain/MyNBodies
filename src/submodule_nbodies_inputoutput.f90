submodule (nbodies) input_output
  implicit none
contains
  module   subroutine read_input_parameters(NB)
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

  module subroutine dump_data(NB)
    ! Now outputs only to CSV format as for some reason paraview cant parse the LAMMPS dump format
    implicit none
    type(NBodies_type), intent(inout) :: NB
    character(len=256) :: filename
    character(len=7) :: step_str
    integer :: i, ios
    idump = idump + 1
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

  ! module subroutine dump_data_old(NB)
  !   implicit none
  !   type(NBodies_type), intent(inout) :: NB

  !   character(len=256) :: filename
  !   character(len=7) :: step_str
  !   integer :: i, ios
  !   idump = idump + 1
  !   ! Create filename with 7-digit timestep number
  !   write(step_str, '(I7.7)') idump
  !   filename = './NBDUMP/nb.' // step_str // '.dump'

  !   ! Open file for writing
  !   open(newunit=ifdump_unit, file=trim(filename), status='replace', action='write', iostat=ios)
  !   if (ios /= 0) then
  !     print *, "Error: Cannot open dump file ", trim(filename)
  !     return
  !   end if

  !   ! Write LAMMPS dump format header
  !   ! write(ifdump_unit, '(A)') 'ITEM: TIMESTEP'
  !   ! write(ifdump_unit, '(I0)') istep

  !   write(ifdump_unit, '(A)') 'ITEM: NUMBER OF ATOMS'
  !   write(ifdump_unit, '(I0)') NB%n_bodies

  !   write(ifdump_unit, '(A)') 'ITEM: BOX BOUNDS pp pp pp'
  !   ! X bounds
  !   if (abs(L_bound(1)) < 1e-10_wp) then
  !     write(ifdump_unit, '(2F12.6)') minval(NB%pos(:,1)) - 1.0_wp, maxval(NB%pos(:,1)) + 1.0_wp
  !   else
  !     write(ifdump_unit, '(2F12.6)') -L_bound(1)*0.5_wp, L_bound(1)*0.5_wp
  !   end if
  !   ! Y bounds
  !   if (abs(L_bound(2)) < 1e-10_wp) then
  !     write(ifdump_unit, '(2F12.6)') minval(NB%pos(:,2)) - 1.0_wp, maxval(NB%pos(:,2)) + 1.0_wp
  !   else
  !     write(ifdump_unit, '(2F12.6)') -L_bound(2)*0.5_wp, L_bound(2)*0.5_wp
  !   end if
  !   ! Z bounds
  !   if (ndim == 3) then
  !     if (abs(L_bound(3)) < 1e-10_wp) then
  !       write(ifdump_unit, '(2F12.6)') minval(NB%pos(:,3)) - 1.0_wp, maxval(NB%pos(:,3)) + 1.0_wp
  !     else
  !       write(ifdump_unit, '(2F12.6)') -L_bound(3)*0.5_wp, L_bound(3)*0.5_wp
  !     end if
  !   else
  !     write(ifdump_unit, '(2F12.6)') 0.0_wp, 0.0_wp
  !   end if

  !   ! Write atom data header based on dimension
  !   if (ndim == 3) then
  !     write(ifdump_unit, '(A)') 'ITEM: ATOMS id type x y vx vy ax ay fx fy mass radius ke pe te'
  !   else if (ndim == 2) then
  !     write(ifdump_unit, '(A)') 'ITEM: ATOMS id type x y z vx vy vz ax ay az fx fy fz mass radius ke pe te'
  !   end if

  !   ! Write atom data
  !   do i = 1, NB%n_bodies
  !     if (ndim == 3) then
  !       write(ifdump_unit, '(I0,1X,I0,1X,17ES14.6)') i, 1, &
  !         NB%pos(i,1), NB%pos(i,2), NB%pos(i,3), &
  !         NB%vel(i,1), NB%vel(i,2), NB%vel(i,3), &
  !         NB%accel(i,1), NB%accel(i,2), NB%accel(i,3), &
  !         NB%force(i,1), NB%force(i,2), NB%force(i,3), &
  !         NB%mass(i), NB%radius(i), &
  !         NB%kinetic_energy(i), NB%potential_energy(i), NB%sum_energy(i)
  !     else if (ndim == 2) then
  !       write(ifdump_unit, '(I0,1X,I0,1X,17ES14.6)') i, 1, &
  !         NB%pos(i,1), NB%pos(i,2), 0.0_wp, &
  !         NB%vel(i,1), NB%vel(i,2), 0.0_wp, &
  !         NB%accel(i,1), NB%accel(i,2), 0.0_wp, &
  !         NB%force(i,1), NB%force(i,2), 0.0_wp, &
  !         NB%mass(i), NB%radius(i), &
  !         NB%kinetic_energy(i), NB%potential_energy(i), NB%sum_energy(i)
  !     end if
  !   end do
  !   ! Close file
  !   close(ifdump_unit)
  !   ! Create filename with 7-digit timestep number
  !   write(step_str, '(I7.7)') idump
  !   filename = './NBDUMPCSV/nb.' // step_str // '.csv'

  !   ! Open file for writing
  !   open(newunit=ifdump_unit, file=trim(filename), status='replace', action='write', iostat=ios)
  !   if (ios /= 0) then
  !     print *, "Error: Cannot open CSV dump file ", trim(filename)
  !     return
  !   end if

  !   ! Write CSV header based on dimension
  !   if (ndim == 3) then
  !     write(ifdump_unit, '(A)') 'id,type,x,y,z,vx,vy,vz,ax,ay,az,fx,fy,fz,mass,radius,ke,pe,te'
  !   else if (ndim == 2) then
  !     write(ifdump_unit, '(A)') 'id,type,x,y,vx,vy,ax,ay,fx,fy,mass,radius,ke,pe,te'
  !   else ! ndim == 1
  !     write(ifdump_unit, '(A)') 'id,type,x,vx,ax,fx,mass,radius,ke,pe,te'
  !   end if

  !   ! Write particle data
  !   do i = 1, NB%n_bodies
  !     if (ndim == 3) then
  !       write(ifdump_unit, '(I0,A,I0,17(A,ES14.6))') i, ',', 1, &
  !         ',', NB%pos(i,1), ',', NB%pos(i,2), ',', NB%pos(i,3), &
  !         ',', NB%vel(i,1), ',', NB%vel(i,2), ',', NB%vel(i,3), &
  !         ',', NB%accel(i,1), ',', NB%accel(i,2), ',', NB%accel(i,3), &
  !         ',', NB%force(i,1), ',', NB%force(i,2), ',', NB%force(i,3), &
  !         ',', NB%mass(i), ',', NB%radius(i), &
  !         ',', NB%kinetic_energy(i), ',', NB%potential_energy(i), ',', NB%sum_energy(i)
  !     else if (ndim == 2) then
  !       write(ifdump_unit, '(I0,A,I0,14(A,ES14.6))') i, ',', 1, &
  !         ',', NB%pos(i,1), ',', NB%pos(i,2), &
  !         ',', NB%vel(i,1), ',', NB%vel(i,2), &
  !         ',', NB%accel(i,1), ',', NB%accel(i,2), &
  !         ',', NB%force(i,1), ',', NB%force(i,2), &
  !         ',', NB%mass(i), ',', NB%radius(i), &
  !         ',', NB%kinetic_energy(i), ',', NB%potential_energy(i), ',', NB%sum_energy(i)
  !     else ! ndim == 1
  !       write(ifdump_unit, '(I0,A,I0,10(A,ES14.6))') i, ',', 1, &
  !         ',', NB%pos(i,1), &
  !         ',', NB%vel(i,1), &
  !         ',', NB%accel(i,1), &
  !         ',', NB%force(i,1), &
  !         ',', NB%mass(i), ',', NB%radius(i), &
  !         ',', NB%kinetic_energy(i), ',', NB%potential_energy(i), ',', NB%sum_energy(i)
  !     end if
  !   end do

  !   ! Close file
  !   close(ifdump_unit)
  ! end subroutine dump_data_old


  module  subroutine output_diagnostics(NB)
    implicit none
    type(NBodies_type), intent(in) :: NB

    print '(A,I8,A,F12.6)', "Step: ", istep, ", Time: ", time
    print '(A,F12.6)',  "  dt:                ", dt
    print '(A,ES15.8)', "  Total KE:          ", NB%total_kinetic_energy
    print '(A,ES15.8)', "  Total PE:          ", NB%total_potential_energy
    print '(A,ES15.8)', "  Total E:           ", NB%total_energy
    print '(A,F15.8,A)', "  Delta E:           ", &
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


end submodule input_output

