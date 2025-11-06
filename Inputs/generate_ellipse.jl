using LinearAlgebra
using Random
using Printf
using Base.Threads
using DelimitedFiles
using Statistics

# === Configuration
const hidden_nstars = 20000
const velocity_perturbation = 0.01  # 1% perturbation by default
const plummer_scale_radius = 1.0   # Plummer scale radius
const axis_ratios = (1.0, 1.0, 1.5)  # (x, y, z) axis ratios

function generate_ellipsoidal_plummer_sphere(n_points::Int, a::Float64=1.0, 
                                            axis_ratios::Tuple{Float64,Float64,Float64}=(1.0, 1.0, 1.5),
                                            M_total::Float64=1.0)
    """
    Generate positions following ellipsoidal Plummer model distribution
    
    Parameters:
    -----------
    n_points : Int
        Number of particles
    a : Float64
        Scale radius (default=1.0)
    axis_ratios : Tuple{Float64,Float64,Float64}
        Axis ratios (ax, ay, az) for ellipsoid
    M_total : Float64
        Total mass (used for velocity scaling, default=1.0)
    
    Returns:
    --------
    positions : Matrix{Float64}, shape (n_points, 3)
        Positions of particles
    """
    
    positions = zeros(n_points, 3)
    ax, ay, az = axis_ratios
    
    @threads for i in 1:n_points
        # Generate radius using inverse transform sampling for Plummer model
        u = rand()
        r = a / sqrt(u^(-2/3) - 1)
        
        # Generate random direction on unit sphere
        cos_theta = 2 * rand() - 1
        sin_theta = sqrt(1 - cos_theta^2)
        phi = 2π * rand()
        
        # Convert to Cartesian coordinates on unit sphere
        x_sphere = sin_theta * cos(phi)
        y_sphere = sin_theta * sin(phi)
        z_sphere = cos_theta
        
        # Scale by ellipsoid axes and radius
        positions[i, 1] = r * ax * x_sphere
        positions[i, 2] = r * ay * y_sphere
        positions[i, 3] = r * az * z_sphere
    end
    
    # Verify COM is at origin
    com = mean(positions, dims=1)
    com_error = norm(vec(com))
    @printf("Ellipsoidal Plummer sphere COM error: %.2e\n", com_error)
    
    # Print some statistics
    radii = [sqrt(positions[i,1]^2/ax^2 + positions[i,2]^2/ay^2 + positions[i,3]^2/az^2) for i in 1:n_points]
    @printf("Mean ellipsoidal radius: %.3f, Median: %.3f, Max: %.3f\n", 
            mean(radii), median(radii), maximum(radii))
    
    return positions
end

function calculate_rotation_parameters(positions::Matrix{Float64}, 
                                     axis_ratios::Tuple{Float64,Float64,Float64},
                                     M_total::Float64=1.0, G::Float64=1.0)
    """
    Calculate appropriate rotation parameters for ellipsoidal system
    
    For a rotating ellipsoid, we'll rotate around the z-axis (longest axis)
    The angular velocity is chosen to provide partial support against gravity
    """
    
    ax, ay, az = axis_ratios
    n_particles = size(positions, 1)
    
    # Calculate moment of inertia tensor components
    Ixx = sum(positions[:, 2].^2 + positions[:, 3].^2) / n_particles
    Iyy = sum(positions[:, 1].^2 + positions[:, 3].^2) / n_particles
    Izz = sum(positions[:, 1].^2 + positions[:, 2].^2) / n_particles
    
    # For stability, choose angular velocity that provides partial centrifugal support
    # This is a fraction of the critical angular velocity
    # We use a conservative fraction to ensure stability
    stability_factor = 0.3  # 30% of critical rotation
    
    # Estimate characteristic angular velocity from virial considerations
    # For a Plummer model, use the characteristic frequency at the scale radius
    omega_char = sqrt(G * M_total / plummer_scale_radius^3)
    
    # Scale by ellipticity and stability factor
    # Rotate around z-axis (longest axis for stability)
    omega_z = stability_factor * omega_char * sqrt(ax/az)
    
    @printf("\nRotation parameters:\n")
    @printf("Characteristic frequency: %.3f\n", omega_char)
    @printf("Angular velocity (z-axis): %.3f\n", omega_z)
    @printf("Rotation period: %.3f\n", 2π/omega_z)
    
    return omega_z
end

function generate_rotating_plummer_velocities(positions::Matrix{Float64}, 
                                            omega_z::Float64,
                                            a::Float64=1.0, 
                                            axis_ratios::Tuple{Float64,Float64,Float64}=(1.0, 1.0, 1.5),
                                            M_total::Float64=1.0, 
                                            G::Float64=1.0, 
                                            perturbation::Float64=0.01)
    """
    Generate velocities for rotating ellipsoidal Plummer model
    
    Parameters:
    -----------
    positions : Matrix{Float64}, shape (n, 3)
        Positions of particles
    omega_z : Float64
        Angular velocity around z-axis
    a : Float64
        Plummer scale radius
    axis_ratios : Tuple{Float64,Float64,Float64}
        Axis ratios for ellipsoid
    M_total : Float64
        Total mass
    G : Float64
        Gravitational constant
    perturbation : Float64
        Fractional velocity perturbation
    
    Returns:
    --------
    velocities : Matrix{Float64}, shape (n, 3)
        Velocities of particles
    """
    
    n_particles = size(positions, 1)
    velocities = zeros(n_particles, 3)
    ax, ay, az = axis_ratios
    
    # Velocity dispersion for Plummer model (adjusted for ellipsoid)
    # Use geometric mean of axes for effective scale
    a_eff = a * (ax * ay * az)^(1/3)
    sigma = sqrt(G * M_total / (6 * a_eff))
    
    # Reduce velocity dispersion to account for rotational support
    sigma *= sqrt(1 - omega_z^2 * a_eff^2 / (G * M_total / a_eff))
    sigma = max(sigma, 0.1 * sqrt(G * M_total / (6 * a_eff)))  # Ensure positive
    
    @threads for i in 1:n_particles
        # Calculate ellipsoidal radius
        r_ell = sqrt(positions[i,1]^2/ax^2 + positions[i,2]^2/ay^2 + positions[i,3]^2/az^2)
        
        # Local escape velocity for ellipsoidal Plummer model
        v_esc = sqrt(2 * G * M_total / sqrt(r_ell^2 + a^2))
        
        # Add rotational velocity (solid body rotation around z-axis)
        v_rot_x = -omega_z * positions[i, 2]
        v_rot_y = omega_z * positions[i, 1]
        v_rot_z = 0.0
        
        # Use rejection sampling for random velocity component
        v_rand_mag = 0.0
        accepted = false
        
        while !accepted
            # Sample from Maxwell-Boltzmann distribution
            v_trial = sigma * sqrt(-2 * log(1 - rand()))
            
            # Accept if total velocity (including rotation) is below escape velocity
            if v_trial < v_esc * 0.8  # Extra safety margin
                v_rand_mag = v_trial
                accepted = true
            end
        end
        
        # Random direction for velocity dispersion
        cos_theta = 2 * rand() - 1
        sin_theta = sqrt(1 - cos_theta^2)
        phi = 2π * rand()
        
        # Random velocity components
        v_rand_x = v_rand_mag * sin_theta * cos(phi)
        v_rand_y = v_rand_mag * sin_theta * sin(phi)
        v_rand_z = v_rand_mag * cos_theta
        
        # Total velocity = rotation + random
        velocities[i, 1] = v_rot_x + v_rand_x
        velocities[i, 2] = v_rot_y + v_rand_y
        velocities[i, 3] = v_rot_z + v_rand_z
        
        # Add perturbation if requested
        if perturbation > 0
            v_total = norm(velocities[i, :])
            random_perturbation = randn(3)
            perturbation_norm = norm(random_perturbation)
            if perturbation_norm > 0
                random_perturbation /= perturbation_norm
                random_perturbation *= v_total * perturbation
                velocities[i, :] += random_perturbation
            end
        end
    end
    
    # Remove net linear momentum (but keep angular momentum)
    mean_velocity = mean(velocities, dims=1)
    velocities .-= mean_velocity
    
    # Calculate angular momentum
    L_total = zeros(3)
    for i in 1:n_particles
        L_total += cross(positions[i, :], velocities[i, :])
    end
    L_total /= n_particles
    
    @printf("\nVelocity statistics:\n")
    @printf("Mean velocity removed: [%.2e, %.2e, %.2e]\n", 
            mean_velocity[1], mean_velocity[2], mean_velocity[3])
    @printf("Angular momentum per particle: [%.3f, %.3f, %.3f]\n",
            L_total[1], L_total[2], L_total[3])
    @printf("L_z / (M * omega_z * <R²>): %.3f (should be ~1 for solid body)\n",
            L_total[3] / (omega_z * mean(positions[:, 1].^2 + positions[:, 2].^2)))
    
    return velocities
end

function calculate_mass_for_virial_equilibrium(points::Matrix{Float64}, velocities::Matrix{Float64}, 
                                              G::Float64=1.0, softening::Float64=0.01)
    """
    Calculate the mass per particle needed to achieve virial equilibrium: 2*KE + PE = 0
    
    Parameters:
    -----------
    points : Matrix{Float64}, shape (n, 3)
        Positions of particles
    velocities : Matrix{Float64}, shape (n, 3)
        Velocities of particles
    G : Float64
        Gravitational constant (default=1.0)
    softening : Float64
        Softening length to avoid singularities (default=0.01)
    
    Returns:
    --------
    mass : Float64
        Mass per particle to achieve virial equilibrium
    """
    n_particles = size(points, 1)

    if n_particles < 2
        error("Need at least 2 particles")
    end

    # Calculate total kinetic energy coefficient (without mass)
    v_squared_sum = sum(velocities.^2)
    ke_coefficient = 0.5 * v_squared_sum

    # Calculate potential energy coefficient in parallel
    pe_coefficient = 0.0
    pe_lock = ReentrantLock()
    
    # Split work among threads more efficiently
    n_threads = Threads.nthreads()
    chunk_size = ceil(Int, n_particles / n_threads)
    
    @threads for thread_id in 1:n_threads
        start_idx = (thread_id - 1) * chunk_size + 1
        end_idx = min(thread_id * chunk_size, n_particles)
        
        local_pe = 0.0
        for i in start_idx:end_idx
            for j in (i+1):n_particles
                r_vec = points[i, :] - points[j, :]
                r_squared = sum(r_vec.^2) + softening^2
                r = sqrt(r_squared)
                local_pe += 1.0 / r
            end
        end
        
        # Thread-safe accumulation
        lock(pe_lock) do
            pe_coefficient += local_pe
        end
    end
    
    pe_coefficient *= -G

    if pe_coefficient >= 0
        error("Potential energy coefficient should be negative")
    end

    # Calculate mass for virial equilibrium
    mass = -v_squared_sum / (G * pe_coefficient)

    # Verify the calculation
    total_ke = 0.5 * mass * v_squared_sum
    total_pe = G * mass^2 * pe_coefficient
    virial = 2 * total_ke + total_pe

    @printf("\nVirial Equilibrium Analysis:\n")
    @printf("Calculated mass per particle: %.6f\n", mass)
    @printf("Total mass: %.6f\n", mass * n_particles)
    @printf("Total KE: %.6f\n", total_ke)
    @printf("Total PE: %.6f\n", total_pe)
    @printf("2*KE + PE (virial): %.6e (should be ≈ 0)\n", virial)
    @printf("KE/|PE| ratio: %.6f (should be ≈ 0.5)\n", total_ke/abs(total_pe))

    return mass
end

function generate_data(perturbation::Float64=0.01)
    # Create data folder if it doesn't exist
    mkpath("data")
    
    n_threads = Threads.nthreads()
    println("Generating rotating ellipsoidal Plummer model data using $n_threads threads...")
    @printf("Number of particles: %d\n", hidden_nstars)
    @printf("Plummer scale radius: %.3f\n", plummer_scale_radius)
    @printf("Axis ratios (x:y:z): %.1f:%.1f:%.1f\n", axis_ratios...)
    @printf("Velocity perturbation: %.1f%%\n", perturbation*100)
    
    # Generate ellipsoidal Plummer sphere positions
    println("\nGenerating ellipsoidal Plummer sphere positions...")
    pos_data = generate_ellipsoidal_plummer_sphere(hidden_nstars, plummer_scale_radius, axis_ratios)
    
    # Save positions as binary
    open("data/positions.dat", "w") do io
        write(io, Matrix(pos_data'))  # Transpose for backward compatibility
    end
    
    # Save positions as ASCII
    writedlm("data/positions_ascii.dat", pos_data, ' ')
    
    # Calculate appropriate rotation rate
    println("\nCalculating rotation parameters...")
    omega_z = calculate_rotation_parameters(pos_data, axis_ratios, 1.0, 1.0)
    
    # Generate velocities with rotation
    println("\nGenerating rotating Plummer model velocities...")
    vel_data = generate_rotating_plummer_velocities(pos_data, omega_z, plummer_scale_radius, 
                                                   axis_ratios, 1.0, 1.0, perturbation)
   
# Save velocities as binary
    open("data/velocities.dat", "w") do io
        write(io, Matrix(vel_data'))  # Transpose for backward compatibility
    end
    
    # Save velocities as ASCII
    writedlm("data/velocities_ascii.dat", vel_data, ' ')
    
    # Calculate mass for virial equilibrium (accounting for rotation)
    println("\nCalculating mass for virial equilibrium...")
    softening = 0.01 * plummer_scale_radius  # Softening as fraction of scale radius
    mass = calculate_mass_for_virial_equilibrium(pos_data, vel_data, 1.0, softening)
    
    # Create masses array (all particles have the same mass)
    masses = fill(mass, hidden_nstars)
    
    # Save masses as binary
    open("data/masses.dat", "w") do io
        write(io, masses)
    end
    
    # Save masses as ASCII
    writedlm("data/masses_ascii.dat", masses, ' ')
    
    # Save metadata
    open("data/nstars.dat", "w") do io
        println(io, hidden_nstars)
        println(io, mass)
        println(io, plummer_scale_radius)
        println(io, softening)
        println(io, "# Axis ratios: ", join(axis_ratios, ", "))
        println(io, "# Angular velocity (z-axis): ", omega_z)
    end
    
    # Calculate and display additional system properties
    println("\n=== System Properties ===")
    
    # Calculate total angular momentum
    L_total = zeros(3)
    for i in 1:hidden_nstars
        L_total += mass * cross(pos_data[i, :], vel_data[i, :])
    end
    L_mag = norm(L_total)
    
    # Calculate rotational and random kinetic energies
    KE_rot = 0.0
    KE_rand = 0.0
    for i in 1:hidden_nstars
        v_rot = [-omega_z * pos_data[i, 2], omega_z * pos_data[i, 1], 0.0]
        v_rand = vel_data[i, :] - v_rot
        KE_rot += 0.5 * mass * dot(v_rot, v_rot)
        KE_rand += 0.5 * mass * dot(v_rand, v_rand)
    end
    
    @printf("Total angular momentum: %.3e\n", L_mag)
    @printf("Angular momentum direction: [%.3f, %.3f, %.3f]\n", 
            L_total[1]/L_mag, L_total[2]/L_mag, L_total[3]/L_mag)
    @printf("Rotational KE: %.3f\n", KE_rot)
    @printf("Random KE: %.3f\n", KE_rand)
    @printf("KE_rot/KE_total ratio: %.3f\n", KE_rot/(KE_rot + KE_rand))
    
    # Print summary statistics
    println("\n=== Summary ===")
    println("No. of stars = ", hidden_nstars)
    println("Mass per star = ", mass)
    println("Total mass = ", mass * hidden_nstars)
    println("Plummer scale radius = ", plummer_scale_radius)
    println("Axis ratios (x:y:z) = ", axis_ratios)
    println("Angular velocity (z-axis) = ", omega_z)
    println("Rotation period = ", 2π/omega_z)
    println("Softening length = ", softening)
    println("Data generation complete!")
end

# Main execution
if abspath(PROGRAM_FILE) == @__FILE__
    # You can change the perturbation here (0.01 = 1%)
    generate_data(velocity_perturbation)
end
