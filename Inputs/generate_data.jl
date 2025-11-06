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

function generate_plummer_sphere(n_points::Int, a::Float64=1.0, M_total::Float64=1.0)
    """
    Generate positions following Plummer model distribution
    
    Parameters:
    -----------
    n_points : Int
        Number of particles
    a : Float64
        Scale radius (default=1.0)
    M_total : Float64
        Total mass (used for velocity scaling, default=1.0)
    
    Returns:
    --------
    positions : Matrix{Float64}, shape (n_points, 3)
        Positions of particles
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
        
        # Convert to Cartesian coordinates
        positions[i, 1] = r * sin_theta * cos(phi)
        positions[i, 2] = r * sin_theta * sin(phi)
        positions[i, 3] = r * cos_theta
    end
    
    # Verify COM is at origin
    com = mean(positions, dims=1)
    com_error = norm(vec(com))
    @printf("Plummer sphere COM error: %.2e\n", com_error)
    
    # Print some statistics
    radii = [norm(positions[i, :]) for i in 1:n_points]
    @printf("Mean radius: %.3f, Median radius: %.3f, Max radius: %.3f\n", 
            mean(radii), median(radii), maximum(radii))
    
    return positions
end

function generate_plummer_velocities(positions::Matrix{Float64}, a::Float64=1.0, 
                                   M_total::Float64=1.0, G::Float64=1.0, 
                                   perturbation::Float64=0.01)
    """
    Generate velocities for Plummer model using rejection sampling
    
    Parameters:
    -----------
    positions : Matrix{Float64}, shape (n, 3)
        Positions of particles
    a : Float64
        Plummer scale radius
    M_total : Float64
        Total mass (for velocity scaling)
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
    
    # Velocity dispersion for Plummer model
    sigma = sqrt(G * M_total / (6 * a))
    
    @threads for i in 1:n_particles
        r = norm(positions[i, :])
        
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
        
        # Random direction for velocity
        cos_theta = 2 * rand() - 1
        sin_theta = sqrt(1 - cos_theta^2)
        phi = 2π * rand()
        
        # Base velocity components
        vx = v_mag * sin_theta * cos(phi)
        vy = v_mag * sin_theta * sin(phi)
        vz = v_mag * cos_theta
        
        velocities[i, :] = [vx, vy, vz]
        
        # Add perturbation if requested
        if perturbation > 0
            random_perturbation = randn(3)
            perturbation_norm = norm(random_perturbation)
            if perturbation_norm > 0
                random_perturbation /= perturbation_norm
                random_perturbation *= v_mag * perturbation
                velocities[i, :] += random_perturbation
            end
        end
    end
    
    # Remove net momentum to ensure COM velocity is zero
    mean_velocity = mean(velocities, dims=1)
    velocities .-= mean_velocity
    
    @printf("Mean velocity removed: [%.2e, %.2e, %.2e]\n", 
            mean_velocity[1], mean_velocity[2], mean_velocity[3])
    
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
    println("Generating Plummer model data using $n_threads threads...")
    @printf("Number of particles: %d\n", hidden_nstars)
    @printf("Plummer scale radius: %.3f\n", plummer_scale_radius)
    @printf("Velocity perturbation: %.1f%%\n", perturbation*100)
    
    # Generate Plummer sphere positions
    println("\nGenerating Plummer sphere positions...")
    pos_data = generate_plummer_sphere(hidden_nstars, plummer_scale_radius)
    
    # Save positions as binary
    open("data/positions.dat", "w") do io
        write(io, Matrix(pos_data'))  # Transpose for backward compatibility
    end
    
    # Save positions as ASCII
    writedlm("data/positions_ascii.dat", pos_data, ' ')
    
    # Generate velocities consistent with Plummer model
    println("\nGenerating Plummer model velocities...")
    # Use unit mass initially for velocity generation
    vel_data = generate_plummer_velocities(pos_data, plummer_scale_radius, 1.0, 1.0, perturbation)
    
    # Save velocities as binary
    open("data/velocities.dat", "w") do io
        write(io, Matrix(vel_data'))  # Transpose for backward compatibility
    end
    
    # Save velocities as ASCII
    writedlm("data/velocities_ascii.dat", vel_data, ' ')
    
    # Calculate mass for virial equilibrium
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
    end
    
    # Print summary statistics
    println("\n=== Summary ===")
    println("No. of stars = ", hidden_nstars)
    println("Mass per star = ", mass)
    println("Total mass = ", mass * hidden_nstars)
    println("Plummer scale radius = ", plummer_scale_radius)
    println("Softening length = ", softening)
    println("Data generation complete!")
end

# Main execution
if abspath(PROGRAM_FILE) == @__FILE__
    # You can change the perturbation here (0.01 = 1%)
    generate_data(velocity_perturbation)
end
