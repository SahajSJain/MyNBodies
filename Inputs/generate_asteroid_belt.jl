using LinearAlgebra
using Random
using Printf
using Base.Threads
using DelimitedFiles
using Statistics

# === Configuration
const total_particles = 20000
const velocity_perturbation = 0.01  # 1% perturbation for stability
const eccentricity_max = 0.2  # Maximum eccentricity for asteroids
const inclination_max = 0.1  # Maximum inclination (radians, ~5.7 degrees)

function generate_central_mass_plummer(n_points::Int, a::Float64=0.5)
    """
    Generate central mass (sun) using Plummer model with small scale radius
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
    @printf("Central mass COM error: %.2e\n", com_error)
    
    return positions
end

function generate_asteroid_belt_distribution(n_points::Int, r_inner::Float64=8.0, r_outer::Float64=12.0,
                                           e_max::Float64=0.2, i_max::Float64=0.1)
    """
    Generate asteroid belt with realistic orbital distribution
    Including eccentricity and inclination variations
    """
    positions = zeros(n_points, 3)
    orbital_elements = zeros(n_points, 6)  # Store a, e, i, Ω, ω, ν for each particle
    
    @threads for i in 1:n_points
        # Semi-major axis - uniform in area
        a = sqrt(rand() * (r_outer^2 - r_inner^2) + r_inner^2)
        
        # Eccentricity - beta distribution favoring circular orbits
        e = e_max * rand()^2  # More circular orbits
        
        # Inclination - exponential distribution favoring low inclinations
        inc = i_max * sqrt(-log(1 - rand() * (1 - exp(-1))))
        if inc > i_max
            inc = i_max * rand()  # Cap at maximum
        end
        
        # Orbital angles
        Ω = 2π * rand()  # Longitude of ascending node
        ω = 2π * rand()  # Argument of periapsis
        ν = 2π * rand()  # True anomaly
        
        # Store orbital elements
        orbital_elements[i, :] = [a, e, inc, Ω, ω, ν]
        
        # Calculate position in orbital plane
        r = a * (1 - e^2) / (1 + e * cos(ν))
        
        # Position in orbital plane
        x_orb = r * cos(ν)
        y_orb = r * sin(ν)
        
        # Rotation matrices for 3D position
        # First rotate by argument of periapsis
        x_temp = x_orb * cos(ω) - y_orb * sin(ω)
        y_temp = x_orb * sin(ω) + y_orb * cos(ω)
        
        # Then by inclination
        x_temp2 = x_temp
        y_temp2 = y_temp * cos(inc)
        z_temp = y_temp * sin(inc)
        
        # Finally by longitude of ascending node
        positions[i, 1] = x_temp2 * cos(Ω) - y_temp2 * sin(Ω)
        positions[i, 2] = x_temp2 * sin(Ω) + y_temp2 * cos(Ω)
        positions[i, 3] = z_temp
    end
    
    # Print statistics
    mean_a = mean(orbital_elements[:, 1])
    mean_e = mean(orbital_elements[:, 2])
    mean_i = mean(orbital_elements[:, 3])
    @printf("Belt orbital elements - Mean a: %.2f, Mean e: %.3f, Mean i: %.3f rad (%.1f°)\n", 
            mean_a, mean_e, mean_i, rad2deg(mean_i))
    
    return positions, orbital_elements
end

function generate_keplerian_velocities(pos_data::Matrix{Float64}, orbital_elements::Matrix{Float64},
                                     M_central::Float64, G::Float64=1.0, perturbation::Float64=0.01)
    """
    Generate Keplerian velocities for particles orbiting central mass
    Using proper orbital mechanics
    """
    n_particles = size(pos_data, 1)
    velocities = zeros(n_particles, 3)
    
    @threads for i in 1:n_particles
        a, e, inc, Ω, ω, ν = orbital_elements[i, :]
        
        # Calculate velocity components in orbital plane
        r = a * (1 - e^2) / (1 + e * cos(ν))
        h = sqrt(G * M_central * a * (1 - e^2))  # Specific angular momentum
        
        # Velocity in orbital plane
        vr = h * e * sin(ν) / (r * (1 - e^2))  # Radial velocity
        vt = h / r  # Tangential velocity
        
        # Convert to Cartesian in orbital plane
        vx_orb = vr * cos(ν) - vt * sin(ν)
        vy_orb = vr * sin(ν) + vt * cos(ν)
        
        # Apply rotation transformations
        # Rotate by argument of periapsis
        vx_temp = vx_orb * cos(ω) - vy_orb * sin(ω)
        vy_temp = vx_orb * sin(ω) + vy_orb * cos(ω)
        
        # Rotate by inclination
        vx_temp2 = vx_temp
        vy_temp2 = vy_temp * cos(inc)
        vz_temp = vy_temp * sin(inc)
        
        # Rotate by longitude of ascending node
        velocities[i, 1] = vx_temp2 * cos(Ω) - vy_temp2 * sin(Ω)
        velocities[i, 2] = vx_temp2 * sin(Ω) + vy_temp2 * cos(Ω)
        velocities[i, 3] = vz_temp
        
        # Add perturbation
        if perturbation > 0
            v_mag = norm(velocities[i, :])
            if v_mag > 0
                random_perturbation = randn(3)
                perturbation_norm = norm(random_perturbation)
                if perturbation_norm > 0
                    random_perturbation /= perturbation_norm
                    random_perturbation *= v_mag * perturbation
                    velocities[i, :] += random_perturbation
                end
            end
        end
    end
    
    return velocities
end

function generate_plummer_velocities(positions::Matrix{Float64}, M_total::Float64, 
                                   a::Float64=0.5, G::Float64=1.0, perturbation::Float64=0.01)
    """
    Generate velocities for Plummer sphere (central mass)
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
        
        # Velocity components
        velocities[i, 1] = v_mag * sin_theta * cos(phi)
        velocities[i, 2] = v_mag * sin_theta * sin(phi)
        velocities[i, 3] = v_mag * cos_theta
        
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
    
    return velocities
end

function calculate_masses_for_stable_system(sun_points::Matrix{Float64}, belt_points::Matrix{Float64},
                                          sun_velocities::Matrix{Float64}, belt_velocities::Matrix{Float64},
                                          G::Float64=1.0, softening::Float64=0.01)
    """
    Calculate masses for a stable sun-asteroid belt system
    """
    n_sun = size(sun_points, 1)
    n_belt = size(belt_points, 1)
    
    # Calculate mean orbital properties of belt
    belt_radii = [norm(belt_points[i, :]) for i in 1:n_belt]
    mean_belt_radius = mean(belt_radii)
    
    # Calculate mean orbital velocity of belt
    belt_v_squared = [sum(belt_velocities[i, :].^2) for i in 1:n_belt]
    mean_v_squared = mean(belt_v_squared)
    
    # For Keplerian orbits: v² = GM/r
    M_sun_total = mean_v_squared * mean_belt_radius / G
    
    # Mass per sun particle
    m_sun = M_sun_total / n_sun
    
    # Belt particles should have much smaller mass
    # Total belt mass ~ 0.001 * M_sun (realistic ratio)
    M_belt_total = 0.001 * M_sun_total
    m_belt = M_belt_total / n_belt
    
    # Calculate energy balance
    # Kinetic energy
    sun_ke = 0.5 * m_sun * sum(sun_velocities.^2)
    belt_ke = 0.5 * m_belt * sum(belt_velocities.^2)
    total_ke = sun_ke + belt_ke
    
    # Potential energy (approximate - sun internal + sun-belt interaction)
    # This is a rough estimate for verification
    sun_internal_pe = -0.5 * G * (m_sun * n_sun)^2 / (0.5 + softening)  # Plummer scale
    sun_belt_pe = -G * M_sun_total * M_belt_total / mean_belt_radius
    total_pe = sun_internal_pe + sun_belt_pe
    
    @printf("\nMass and Energy Analysis:\n")
    @printf("Mean belt radius: %.3f\n", mean_belt_radius)
    @printf("Mean belt velocity²: %.3f\n", mean_v_squared)
    @printf("Total sun mass: %.6f\n", M_sun_total)
    @printf("Total belt mass: %.6f (%.3f%% of sun)\n", M_belt_total, 100*M_belt_total/M_sun_total)
    @printf("Mass per sun particle: %.6e\n", m_sun)
    @printf("Mass per belt particle: %.6e\n", m_belt)
    @printf("\nEnergy components:\n")
    @printf("Sun KE: %.6f\n", sun_ke)
    @printf("Belt KE: %.6f\n", belt_ke)
    @printf("Total KE: %.6f\n", total_ke)
    @printf("Approx Total PE: %.6f\n", total_pe)
    @printf("Total Energy: %.6f\n", total_ke + total_pe)
    
    return m_sun, m_belt
end

function generate_asteroid_belt_data(perturbation::Float64=0.01)
    # Create data folder if it doesn't exist
    mkpath("data")
    
    n_threads = Threads.nthreads()
    println("Generating asteroid belt simulation data using $n_threads threads...")
    @printf("Velocity perturbation: %.1f%%\n", perturbation*100)
    
    # Distribution of particles: 10% for sun, 90% for belt
    n_sun = total_particles ÷ 10
    n_belt = total_particles - n_sun
    
    println("\nGenerating central mass (sun) with $n_sun particles...")
    # Generate sun particles using Plummer model
    sun_points = generate_central_mass_plummer(n_sun, 0.5)
    
    println("\nGenerating asteroid belt with $n_belt particles...")
    # Generate belt particles with realistic distribution
    belt_points, orbital_elements = generate_asteroid_belt_distribution(n_belt, 8.0, 12.0, 
                                                                      eccentricity_max, inclination_max)
    
    # Combine all positions
    pos_data = vcat(sun_points, belt_points)
    
    # Generate velocities
    println("\nGenerating velocities...")
    
    # First calculate sun mass based on Kepler's laws
    # Assume unit system where G = 1
    M_sun_assumed = 1.0
    
    # Sun particles: Plummer sphere velocities
    sun_velocities = generate_plummer_velocities(sun_points, M_sun_assumed, 0.5, 1.0, perturbation)
    
    # Belt particles: Keplerian velocities
    belt_velocities = generate_keplerian_velocities(belt_points, orbital_elements, 
                                                  M_sun_assumed, 1.0, perturbation)
    
    # Combine all velocities
    vel_data = vcat(sun_velocities, belt_velocities)
    
    # Calculate appropriate masses
    println("\nCalculating masses for stable system...")
    m_sun, m_belt = calculate_masses_for_stable_system(sun_points, belt_points,
                                                       sun_velocities, belt_velocities, 1.0, 0.01)
    
    # Create masses array
    masses = vcat(fill(m_sun, n_sun), fill(m_belt, n_belt))
    
    # Remove center of mass motion
    total_mass = sum(masses)
    com_velocity = sum(masses .* vel_data, dims=1) / total_mass
    vel_data .-= com_velocity
   
# Save positions
    open("data/positions.dat", "w") do io
        write(io, Matrix(pos_data'))  # Transpose for backward compatibility
    end
    writedlm("data/positions_ascii.dat", pos_data, ' ')

    # Save velocities
    open("data/velocities.dat", "w") do io
        write(io, Matrix(vel_data'))  # Transpose for backward compatibility
    end
    writedlm("data/velocities_ascii.dat", vel_data, ' ')

    # Save masses
    open("data/masses.dat", "w") do io
        write(io, masses)
    end
    writedlm("data/masses_ascii.dat", masses, ' ')

    # Save metadata
    open("data/nstars.dat", "w") do io
        println(io, total_particles)
        println(io, m_sun)
        println(io, m_belt)
        println(io, n_sun)
        println(io, n_belt)
    end

    # Save orbital elements for belt particles (useful for analysis)
    if n_belt > 0
        writedlm("data/orbital_elements.dat", orbital_elements, ' ')
    end

    # Print summary statistics
    println("\n=== Simulation Summary ===")
    println("Total particles: ", total_particles)
    println("Sun particles: ", n_sun, " (", round(100*n_sun/total_particles, digits=1), "%)")
    println("Belt particles: ", n_belt, " (", round(100*n_belt/total_particles, digits=1), "%)")
    println("Sun particle mass: ", m_sun)
    println("Belt particle mass: ", m_belt)
    println("Total sun mass: ", m_sun * n_sun)
    println("Total belt mass: ", m_belt * n_belt)
    println("Mass ratio (belt/sun): ", (m_belt * n_belt)/(m_sun * n_sun))
    println("\nSun properties:")
    println("  Plummer scale radius: 0.5")
    println("  Mean radius: ", mean([norm(sun_points[i, :]) for i in 1:n_sun]))
    println("\nBelt properties:")
    println("  Inner radius: 8.0")
    println("  Outer radius: 12.0")
    println("  Max eccentricity: ", eccentricity_max)
    println("  Max inclination: ", inclination_max, " rad (", round(rad2deg(inclination_max), digits=1), "°)")

    # Verify center of mass
    com = sum(masses .* pos_data, dims=1) / total_mass
    println("\nSystem center of mass: [", @sprintf("%.2e", com[1]), ", ",
            @sprintf("%.2e", com[2]), ", ", @sprintf("%.2e", com[3]), "]")
    com_vel = sum(masses .* vel_data, dims=1) / total_mass
    println("System COM velocity: [", @sprintf("%.2e", com_vel[1]), ", ",
            @sprintf("%.2e", com_vel[2]), ", ", @sprintf("%.2e", com_vel[3]), "]")

    # Calculate and display orbital periods
    if n_belt > 0
        mean_a = mean(orbital_elements[:, 1])
        min_a = minimum(orbital_elements[:, 1])
        max_a = maximum(orbital_elements[:, 1])

        # Kepler's third law: T = 2π√(a³/GM)
        mean_period = 2π * sqrt(mean_a^3 / (m_sun * n_sun))
        min_period = 2π * sqrt(min_a^3 / (m_sun * n_sun))
        max_period = 2π * sqrt(max_a^3 / (m_sun * n_sun))

        println("\nOrbital periods:")
        println("  Mean: ", round(mean_period, digits=2), " time units")
        println("  Range: [", round(min_period, digits=2), ", ", round(max_period, digits=2), "]")

        # Velocity statistics
        belt_speeds = [norm(belt_velocities[i, :]) for i in 1:n_belt]
        println("\nBelt velocity statistics:")
        println("  Mean speed: ", round(mean(belt_speeds), digits=3))
        println("  Min speed: ", round(minimum(belt_speeds), digits=3))
        println("  Max speed: ", round(maximum(belt_speeds), digits=3))
    end

    # Display particle type indicators (for visualization)
    particle_types = vcat(fill(1, n_sun), fill(2, n_belt))  # 1 = sun, 2 = belt
    writedlm("data/particle_types.dat", particle_types, ' ')

    println("\nData generation complete!")
    println("Files created in 'data/' directory:")
    println("  - positions.dat (binary)")
    println("  - positions_ascii.dat")
    println("  - velocities.dat (binary)")
    println("  - velocities_ascii.dat")
    println("  - masses.dat (binary)")
    println("  - masses_ascii.dat")
    println("  - nstars.dat (metadata)")
    println("  - orbital_elements.dat (belt orbital elements)")
    println("  - particle_types.dat (1=sun, 2=belt)")
end

# Main execution
if abspath(PROGRAM_FILE) == @__FILE__
    # You can change the perturbation here (0.01 = 1%)
    generate_asteroid_belt_data(velocity_perturbation)
end
