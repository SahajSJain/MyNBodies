
using LinearAlgebra
using Random
using Printf
using Base.Threads
using DelimitedFiles
using Statistics

# === Configuration
const total_particles = 21000  # Divisible by 3
const velocity_perturbation = 0.01  # 1% perturbation for stability
const plummer_scale_radius = 0.5   # Scale radius for each body
const separation_radius = 10.0     # Distance from COM to each body center
const softening = 0.01            # Softening parameter

function generate_plummer_sphere(n_points::Int, a::Float64=0.5, center::Vector{Float64}=[0.0, 0.0, 0.0])
    """
    Generate positions following Plummer model distribution centered at given position
    """
    positions = zeros(n_points, 3)

    @threads for i in 1:n_points
        # Generate radius using inverse transform sampling for Plummer model
        u = rand()
        r = a / sqrt(u^(-2/3) - 1)

        # Generate random direction
        cos_theta = 2 * rand() - 1
        sin_theta = sqrt(1 - cos_theta^2)
        phi = 2π * rand()

        # Convert to Cartesian coordinates and add center offset
        positions[i, 1] = r * sin_theta * cos(phi) + center[1]
        positions[i, 2] = r * sin_theta * sin(phi) + center[2]
        positions[i, 3] = r * cos_theta + center[3]
    end

    return positions
end

function generate_plummer_velocities(positions::Matrix{Float64}, center::Vector{Float64},
                                   M_total::Float64, a::Float64=0.5, G::Float64=1.0)
    """
    Generate velocities for Plummer sphere with internal virial equilibrium
    """
    n_particles = size(positions, 1)
    velocities = zeros(n_particles, 3)

    # Velocity dispersion for Plummer model
    sigma = sqrt(G * M_total / (6 * a))

    @threads for i in 1:n_particles
        # Distance from sphere center
        r_vec = positions[i, :] - center
        r = norm(r_vec)

        # Local escape velocity for Plummer model
        v_esc = sqrt(2 * G * M_total / sqrt(r^2 + a^2))

        # Use rejection sampling to get velocity magnitude
        v_mag = 0.0
        accepted = false

        while !accepted
            # Sample from Maxwell-Boltzmann distribution
            v_trial = sigma * sqrt(-2 * log(1 - rand()))

            # Accept if below escape velocity
            if v_trial < v_esc
                v_mag = v_trial
                accepted = true
            end
        end

        # Random direction for velocity (in sphere's rest frame)
        cos_theta = 2 * rand() - 1
        sin_theta = sqrt(1 - cos_theta^2)
        phi = 2π * rand()

        # Velocity components in sphere's rest frame
        velocities[i, 1] = v_mag * sin_theta * cos(phi)
        velocities[i, 2] = v_mag * sin_theta * sin(phi)
        velocities[i, 3] = v_mag * cos_theta
    end

    return velocities
end

function calculate_lagrange_configuration(R::Float64, M_body::Float64, G::Float64=1.0)
    """
    Calculate positions and velocities for Lagrange equilateral triangle configuration
    """
    # Positions of three bodies at vertices of equilateral triangle
    centers = [
        [R, 0.0, 0.0],                           # Body 1
        [-R/2, R*sqrt(3)/2, 0.0],               # Body 2
        [-R/2, -R*sqrt(3)/2, 0.0]               # Body 3
    ]

    # Angular velocity for circular orbits
    # For three equal masses in equilateral triangle: ω = sqrt(3GM/(R³))
    omega = sqrt(3 * G * M_body / R^3)

    # Velocity of each body (perpendicular to position vector)
    center_velocities = [
        [0.0, omega * R, 0.0],                   # Body 1
        [-omega * R * sqrt(3)/2, -omega * R/2, 0.0],    # Body 2
        [omega * R * sqrt(3)/2, -omega * R/2, 0.0]      # Body 3
    ]

    # Verify COM is at origin
    com_pos = (centers[1] + centers[2] + centers[3]) / 3
    com_vel = (center_velocities[1] + center_velocities[2] + center_velocities[3]) / 3

    @printf("COM position: [%.2e, %.2e, %.2e]\n", com_pos[1], com_pos[2], com_pos[3])
    @printf("COM velocity: [%.2e, %.2e, %.2e]\n", com_vel[1], com_vel[2], com_vel[3])
    @printf("Angular velocity: %.6f\n", omega)
    @printf("Orbital period: %.6f\n", 2π/omega)

    return centers, center_velocities, omega
end

function generate_three_body_system(perturbation::Float64=0.01)
    # Create data folder if it doesn't exist
    mkpath("data")

    n_threads = Threads.nthreads()
    println("Generating three-body Lagrange system using $n_threads threads...")
    @printf("Number of particles: %d\n", total_particles)
    @printf("Particles per body: %d\n", total_particles ÷ 3)
    @printf("Plummer scale radius: %.3f\n", plummer_scale_radius)
    @printf("Separation radius: %.3f\n", separation_radius)
    @printf("Velocity perturbation: %.1f%%\n", perturbation*100)

    # Ensure total particles is divisible by 3
    if total_particles % 3 != 0
        error("Total particles must be divisible by 3")
    end

    n_per_body = total_particles ÷ 3

    # First, calculate mass per body for desired configuration
    # We'll use unit mass per particle initially
    M_body_initial = Float64(n_per_body)

    # Get Lagrange configuration
    println("\nCalculating Lagrange configuration...")
    centers, center_velocities, omega = calculate_lagrange_configuration(separation_radius, M_body_initial, 1.0)

    # Generate three Plummer spheres
    println("\nGenerating three Plummer spheres...")
    all_positions = zeros(total_particles, 3)
    all_velocities = zeros(total_particles, 3)

    for i in 1:3
        println("  Generating body $i...")
        start_idx = (i-1) * n_per_body + 1
        end_idx = i * n_per_body

        # Generate positions for this body
        body_positions = generate_plummer_sphere(n_per_body, plummer_scale_radius, centers[i])
        all_positions[start_idx:end_idx, :] = body_positions

        # Generate internal velocities (in body's rest frame)
        internal_velocities = generate_plummer_velocities(body_positions, centers[i],
                                                        M_body_initial, plummer_scale_radius, 1.0)

        # Add bulk motion of the body
        for j in 1:n_per_body
            all_velocities[start_idx+j-1, :] = internal_velocities[j, :] + center_velocities[i]
        end
    end

    # Add velocity perturbation if requested
    if perturbation > 0
        println("\nAdding velocity perturbations...")
        for i in 1:total_particles
            v_mag = norm(all_velocities[i, :])
            if v_mag > 0
                random_perturbation = randn(3)
                perturbation_norm = norm(random_perturbation)
                if perturbation_norm > 0
                    random_perturbation /= perturbation_norm
                    random_perturbation *= v_mag * perturbation
                    all_velocities[i, :] += random_perturbation
                end
            end
        end
    end

    # Calculate mass for overall stability
    println("\nCalculating mass for system stability...")

    # For the three-body system, we need to balance:
    # 1. Internal binding energy of each Plummer sphere
    # 2. Orbital kinetic energy
    # 3. Gravitational interaction between bodies

    # Calculate kinetic energy coefficient
    v_squared_sum = sum(all_velocities.^2)
    ke_coefficient = 0.5 * v_squared_sum

    # Calculate potential energy coefficient
    pe_coefficient = 0.0
    pe_lock = ReentrantLock()

    # Internal PE of each body (Plummer sphere)
    for body in 1:3
        start_idx = (body-1) * n_per_body + 1
        end_idx = body * n_per_body

        @threads for i in start_idx:end_idx
            local_pe = 0.0
            for j in (i+1):end_idx
                r_vec = all_positions[i, :] - all_positions[j, :]
                r_squared = sum(r_vec.^2) + softening^2
                r = sqrt(r_squared)
                local_pe += 1.0 / r
            end

            lock(pe_lock) do
                pe_coefficient += local_pe
            end
        end
    end

    # PE between different bodies
    for body1 in 1:2
        for body2 in (body1+1):3
            start1 = (body1-1) * n_per_body + 1
            end1 = body1 * n_per_body
            start2 = (body2-1) * n_per_body + 1
            end2 = body2 * n_per_body

            @threads for i in start1:end1
                local_pe = 0.0
                for j in start2:end2
                    r_vec = all_positions[i, :] - all_positions[j, :]
                    r_squared = sum(r_vec.^2) + softening^2
                    r = sqrt(r_squared)
                    local_pe += 1.0 / r
                end

                lock(pe_lock) do
                    pe_coefficient += local_pe
                end
            end
        end
    end

    pe_coefficient *= -1.0  # G = 1

    # For stable system, we want the orbital motion to be correct
    # The mass should be set so that ω² = 3GM/(R³)
    # We already have the velocities set up for this
    mass_per_particle = (omega^2 * separation_radius^3) / (3 * n_per_body)

    @printf("\nMass calculation:\n")
    @printf("Mass per particle: %.6e\n", mass_per_particle)
    @printf("Mass per body: %.6f\n", mass_per_particle * n_per_body)
    @printf("Total system mass: %.6f\n", mass_per_particle * total_particles)

    # Verify energy balance
    total_ke = 0.5 * mass_per_particle * v_squared_sum
    total_pe = mass_per_particle^2 * pe_coefficient
    total_energy = total_ke + total_pe

    @printf("\nEnergy analysis:\n")
    @printf("Total KE: %.6f\n", total_ke)
    @printf("Total PE: %.6f\n", total_pe)
    @printf("Total Energy: %.6f\n", total_energy)
    @printf("Virial ratio (2KE/|PE|): %.6f\n", 2*total_ke/abs(total_pe))

    # Create masses array
    masses = fill(mass_per_particle, total_particles)

    # Final COM check
    com_pos = sum(masses .* all_positions, dims=1) / sum(masses)
    com_vel = sum(masses .* all_velocities, dims=1) / sum(masses)

    # Save all data
    println("\nSaving data...")

    # Positions
    open("data/positions.dat", "w") do io
        write(io, Matrix(all_positions'))
    end
    writedlm("data/positions_ascii.dat", all_positions, ' ')

    # Velocities
    open("data/velocities.dat", "w") do io
        write(io, Matrix(all_velocities'))
    end
    writedlm("data/velocities_ascii.dat", all_velocities, ' ')

    # Masses
    open("data/masses.dat", "w") do io
        write(io, masses)
    end
    writedlm("data/masses_ascii.dat", masses, ' ')

    # Particle types (1, 2, or 3 for each body)
    particle_types = vcat([fill(i, n_per_body) for i in 1:3]...)
    writedlm("data/particle_types.dat", particle_types, ' ')

    # Metadata
    open("data/nstars.dat", "w") do io
        println(io, total_particles)
        println(io, mass_per_particle)
        println(io, n_per_body)
        println(io, plummer_scale_radius)
        println(io, separation_radius)
        println(io, omega)
    end

    # Summary
    println("\n=== Three-Body System Summary ===")
    println("Total particles: ", total_particles)
    println("Particles per body: ", n_per_body)
    println("Mass per particle: ", mass_per_particle)
    println("Plummer scale radius: ", plummer_scale_radius)
    println("Separation radius: ", separation_radius)
    println("Angular velocity: ", omega)
    println("Orbital period: ", 2π/omega)
    println("System COM: [", @sprintf("%.2e", com_pos[1]), ", ",
            @sprintf("%.2e", com_pos[2]), ", ", @sprintf("%.2e", com_pos[3]), "]")
    println("System COM velocity: [", @sprintf("%.2e", com_vel[1]), ", ",
            @sprintf("%.2e", com_vel[2]), ", ", @sprintf("%.2e", com_vel[3]), "]")

    println("\nData generation complete!")
end

# Main execution
if abspath(PROGRAM_FILE) == @__FILE__
    generate_three_body_system(velocity_perturbation)
end
