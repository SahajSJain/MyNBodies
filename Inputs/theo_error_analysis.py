
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import glob
from scipy.integrate import solve_ivp, cumulative_trapezoid as cumtrapz
import warnings
warnings.filterwarnings('ignore')

# Create output directory
os.makedirs('./Plots', exist_ok=True)

# Set matplotlib parameters for better plots
plt.rcParams['font.size'] = 14
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['legend.fontsize'] = 12
plt.rcParams['figure.figsize'] = (6, 6)

def read_time_history():
    """Read the Time_History.dat file"""
    data = []
    with open('Time_History.dat', 'r') as f:
        # Skip header
        f.readline()
        for line in f:
            parts = line.strip().split()
            if len(parts) > 0:
                data.append([int(parts[0]), float(parts[1])] + [float(p) for p in parts[2:]])

    df = pd.DataFrame(data, columns=['step', 'time', 'COM_x', 'COM_y', 'COM_z',
                                     'COM_vx', 'COM_vy', 'COM_vz', 'KE_total',
                                     'PE_total', 'E_total'])
    return df

def read_csv_files():
    """Read all nb.*.csv files and organize by timestep"""
    csv_files = sorted(glob.glob('NBDUMPCSV/nb.*.csv'))

    all_data = {}
    for file in csv_files:
        # Extract file number
        file_num = int(file.split('.')[-2])

        # Read CSV
        df = pd.read_csv(file)
        all_data[file_num] = df

    return all_data

def get_initial_conditions(csv_data):
    """Extract initial conditions from the first CSV file"""
    df = csv_data[0]

    # Get particle 1 data
    p1 = df[df['id'] == 1].iloc[0]
    m1 = p1['mass']
    r1_0 = np.array([p1['x'], p1['y'], p1['z']])
    v1_0 = np.array([p1['vx'], p1['vy'], p1['vz']])

    # Get particle 2 data
    p2 = df[df['id'] == 2].iloc[0]
    m2 = p2['mass']
    r2_0 = np.array([p2['x'], p2['y'], p2['z']])
    v2_0 = np.array([p2['vx'], p2['vy'], p2['vz']])

    return {
        'm1': m1, 'm2': m2,
        'r1_0': r1_0, 'r2_0': r2_0,
        'v1_0': v1_0, 'v2_0': v2_0
    }

def two_body_derivatives(t, state, m1, m2, G=1):
    """
    Calculate derivatives for the two-body problem
    state = [x1, y1, z1, x2, y2, z2, vx1, vy1, vz1, vx2, vy2, vz2]
    """
    # Unpack positions and velocities
    r1 = state[0:3]
    r2 = state[3:6]
    v1 = state[6:9]
    v2 = state[9:12]

    # Calculate relative position and distance
    r = r2 - r1
    r_mag = np.linalg.norm(r)

    # Calculate accelerations
    a1 = G * m2 * r / r_mag**3
    a2 = -G * m1 * r / r_mag**3

    # Return derivatives [velocities, accelerations]
    return np.concatenate([v1, v2, a1, a2])

def calculate_theoretical_trajectory(initial_conditions, times, G=1):
    """
    Calculate theoretical trajectory using numerical integration of the two-body problem
    """
    # Extract initial conditions
    m1 = initial_conditions['m1']
    m2 = initial_conditions['m2']
    r1_0 = initial_conditions['r1_0']
    r2_0 = initial_conditions['r2_0']
    v1_0 = initial_conditions['v1_0']
    v2_0 = initial_conditions['v2_0']

    # Initial state vector
    y0 = np.concatenate([r1_0, r2_0, v1_0, v2_0])

    # Solve ODE
    sol = solve_ivp(two_body_derivatives, [times[0], times[-1]], y0,
                    t_eval=times, args=(m1, m2, G),
                    method='DOP853', rtol=1e-12, atol=1e-14)

    # Extract results
    results = {
        'particle1': {
            'x': sol.y[0], 'y': sol.y[1], 'z': sol.y[2],
            'vx': sol.y[6], 'vy': sol.y[7], 'vz': sol.y[8]
        },
        'particle2': {
            'x': sol.y[3], 'y': sol.y[4], 'z': sol.y[5],
            'vx': sol.y[9], 'vy': sol.y[10], 'vz': sol.y[11]
        }
    }

    # Calculate accelerations, forces, and energies
    for i, t in enumerate(times):
        state = sol.y[:, i]
        derivs = two_body_derivatives(t, state, m1, m2, G)

        # Accelerations
        results['particle1']['ax'] = results['particle1'].get('ax', [])
        results['particle1']['ay'] = results['particle1'].get('ay', [])
        results['particle1']['az'] = results['particle1'].get('az', [])
        results['particle1']['ax'] = np.append(results['particle1']['ax'], derivs[6])
        results['particle1']['ay'] = np.append(results['particle1']['ay'], derivs[7])
        results['particle1']['az'] = np.append(results['particle1']['az'], derivs[8])

        results['particle2']['ax'] = results['particle2'].get('ax', [])
        results['particle2']['ay'] = results['particle2'].get('ay', [])
        results['particle2']['az'] = results['particle2'].get('az', [])
        results['particle2']['ax'] = np.append(results['particle2']['ax'], derivs[9])
        results['particle2']['ay'] = np.append(results['particle2']['ay'], derivs[10])
        results['particle2']['az'] = np.append(results['particle2']['az'], derivs[11])

        # Forces
        results['particle1']['fx'] = results['particle1'].get('fx', [])
        results['particle1']['fy'] = results['particle1'].get('fy', [])
        results['particle1']['fz'] = results['particle1'].get('fz', [])
        results['particle1']['fx'] = np.append(results['particle1']['fx'], m1 * derivs[6])
        results['particle1']['fy'] = np.append(results['particle1']['fy'], m1 * derivs[7])
        results['particle1']['fz'] = np.append(results['particle1']['fz'], m1 * derivs[8])

        results['particle2']['fx'] = results['particle2'].get('fx', [])
        results['particle2']['fy'] = results['particle2'].get('fy', [])
        results['particle2']['fz'] = results['particle2'].get('fz', [])
        results['particle2']['fx'] = np.append(results['particle2']['fx'], m2 * derivs[9])
        results['particle2']['fy'] = np.append(results['particle2']['fy'], m2 * derivs[10])
        results['particle2']['fz'] = np.append(results['particle2']['fz'], m2 * derivs[11])

        # Kinetic energies
        v1_mag = np.sqrt(state[6]**2 + state[7]**2 + state[8]**2)
        v2_mag = np.sqrt(state[9]**2 + state[10]**2 + state[11]**2)

        results['particle1']['ke'] = results['particle1'].get('ke', [])
        results['particle2']['ke'] = results['particle2'].get('ke', [])
        results['particle1']['ke'] = np.append(results['particle1']['ke'], 0.5 * m1 * v1_mag**2)
        results['particle2']['ke'] = np.append(results['particle2']['ke'], 0.5 * m2 * v2_mag**2)

        # Potential energy (shared equally between particles)
        r = np.sqrt((state[3]-state[0])**2 + (state[4]-state[1])**2 + (state[5]-state[2])**2)
        pe_total = -G * m1 * m2 / r

        results['particle1']['pe'] = results['particle1'].get('pe', [])
        results['particle2']['pe'] = results['particle2'].get('pe', [])
        results['particle1']['pe'] = np.append(results['particle1']['pe'], pe_total / 2)
        results['particle2']['pe'] = np.append(results['particle2']['pe'], pe_total / 2)

        # Total energy
        results['particle1']['te'] = results['particle1'].get('te', [])
        results['particle2']['te'] = results['particle2'].get('te', [])
        results['particle1']['te'] = np.append(results['particle1']['te'],
                                               results['particle1']['ke'][i] + results['particle1']['pe'][i])
        results['particle2']['te'] = np.append(results['particle2']['te'],
                                               results['particle2']['ke'][i] + results['particle2']['pe'][i])

    return results

def calculate_errors(sim_data, theory_data, times):
    """Calculate mean and RMS integrated errors"""
    diff = sim_data - theory_data

    # Integrate absolute error over time
    integrated_error = cumtrapz(np.abs(diff), times, initial=0)
    mean_integrated_error = integrated_error[-1] / (times[-1] - times[0])

    # RMS integrated error
    integrated_squared_error = cumtrapz(diff**2, times, initial=0)
    rms_integrated_error = np.sqrt(integrated_squared_error[-1] / (times[-1] - times[0]))

    return mean_integrated_error, rms_integrated_error, diff

def plot_comparison(times, sim_data, theory_data, diff, ylabel, filename, title):
    """Create comparison plot with dual y-axes"""
    fig, ax1 = plt.subplots(figsize=(6, 6))

    # Left y-axis: theoretical vs simulation
    color1 = 'tab:blue'
    color2 = 'tab:red'
    ax1.set_xlabel('Time')
    ax1.set_ylabel(ylabel, color='black')
    line1 = ax1.plot(times, theory_data, color=color1, linewidth=2, label='Theoretical')
    line2 = ax1.plot(times, sim_data, color=color2, linewidth=2, linestyle='--', label='Simulation')
    ax1.tick_params(axis='y')
    ax1.grid(True, alpha=0.3)

    # Right y-axis: difference
    ax2 = ax1.twinx()
    color3 = 'tab:green'
    ax2.set_ylabel('Difference', color=color3)
    line3 = ax2.plot(times, diff, color=color3, linewidth=1, alpha=0.7, label='Difference')
    ax2.tick_params(axis='y', labelcolor=color3)

    # Combine legends
    lines = line1 + line2 + line3
    labels = [l.get_label() for l in lines]
    ax1.legend(lines, labels, loc='best')

    plt.title(title)
    plt.tight_layout()
    plt.savefig(f'./Plots/{filename}', dpi=300, bbox_inches='tight')
    plt.close()

def main():
    print("Reading Time_History.dat...")
    time_history = read_time_history()

    print("Reading CSV files...")
    csv_data = read_csv_files()

    print("Extracting initial conditions...")
    initial_conditions = get_initial_conditions(csv_data)

    print(f"Initial conditions:")
    print(f"  Particle 1: mass = {initial_conditions['m1']}")
    print(f"              position = {initial_conditions['r1_0']}")
    print(f"              velocity = {initial_conditions['v1_0']}")
    print(f"  Particle 2: mass = {initial_conditions['m2']}")
    print(f"              position = {initial_conditions['r2_0']}")
    print(f"              velocity = {initial_conditions['v2_0']}")

    # Initialize storage for time series data
    times = []
    particle1_data = {key: [] for key in ['x', 'y', 'z', 'vx', 'vy', 'vz',
                                          'ax', 'ay', 'az', 'fx', 'fy', 'fz',
                                          'ke', 'pe', 'te']}
    particle2_data = {key: [] for key in particle1_data.keys()}

    # Extract time series data from simulation
    for file_num in sorted(csv_data.keys()):
        # Get corresponding time from time_history
        if file_num < len(time_history):
            time = time_history.iloc[file_num]['time']
            times.append(time)

            # Get particle data
            df = csv_data[file_num]

            # Particle 1 (id=1)
            p1 = df[df['id'] == 1].iloc[0]
            for key in particle1_data.keys():
                particle1_data[key].append(p1[key])

            # Particle 2 (id=2)
            p2 = df[df['id'] == 2].iloc[0]
            for key in particle2_data.keys():
                particle2_data[key].append(p2[key])

    times = np.array(times)

    # Convert to numpy arrays
    for key in particle1_data.keys():
        particle1_data[key] = np.array(particle1_data[key])
        particle2_data[key] = np.array(particle2_data[key])

    print("Calculating theoretical trajectories...")
    theory_results = calculate_theoretical_trajectory(initial_conditions, times)
    with open('./Plots/error_analysis.dat', 'w') as f:
        f.write("Error Analysis for N-body Simulation vs Theoretical Predictions\n")
        f.write("="*70 + "\n\n")
        f.write(f"Initial Conditions:\n")
        f.write(f"  Particle 1: mass = {initial_conditions['m1']:.6f}\n")
        f.write(f"              position = ({initial_conditions['r1_0'][0]:.6e}, {initial_conditions['r1_0'][1]:.6e}, {initial_conditions['r1_0'][2]:.6e})\n")
        f.write(f"              velocity = ({initial_conditions['v1_0'][0]:.6e}, {initial_conditions['v1_0'][1]:.6e}, {initial_conditions['v1_0'][2]:.6e})\n")
        f.write(f"  Particle 2: mass = {initial_conditions['m2']:.6f}\n")
        f.write(f"              position = ({initial_conditions['r2_0'][0]:.6e}, {initial_conditions['r2_0'][1]:.6e}, {initial_conditions['r2_0'][2]:.6e})\n")
        f.write(f"              velocity = ({initial_conditions['v2_0'][0]:.6e}, {initial_conditions['v2_0'][1]:.6e}, {initial_conditions['v2_0'][2]:.6e})\n")
        f.write("\n" + "="*70 + "\n\n")
        
        # Analyze and plot for each quantity
        quantities = [
            ('x', 'Position X (m)', 'position_x'),
            ('y', 'Position Y (m)', 'position_y'),
            ('z', 'Position Z (m)', 'position_z'),
            ('vx', 'Velocity X (m/s)', 'velocity_x'),
            ('vy', 'Velocity Y (m/s)', 'velocity_y'),
            ('vz', 'Velocity Z (m/s)', 'velocity_z'),
            ('ke', 'Kinetic Energy', 'kinetic_energy'),
            ('pe', 'Potential Energy', 'potential_energy'),
            ('te', 'Total Energy', 'total_energy'),
            ('fx', 'Force X', 'force_x'),
            ('fy', 'Force Y', 'force_y'),
            ('fz', 'Force Z', 'force_z')
        ]
        
        for key, ylabel, filename in quantities:
            # Particle 1
            mean_err1, rms_err1, diff1 = calculate_errors(
                particle1_data[key], theory_results['particle1'][key], times
            )
            plot_comparison(times, particle1_data[key], theory_results['particle1'][key], 
                          diff1, ylabel, f'{filename}_particle1.png', 
                          f'Particle 1: {ylabel}')
            
            # Particle 2
            mean_err2, rms_err2, diff2 = calculate_errors(
                particle2_data[key], theory_results['particle2'][key], times
            )
            plot_comparison(times, particle2_data[key], theory_results['particle2'][key], 
                          diff2, ylabel, f'{filename}_particle2.png', 
                          f'Particle 2: {ylabel}')
            
            # Write errors to file
            f.write(f"{ylabel}:\n")
            f.write(f"  Particle 1 - Mean Integrated Error: {mean_err1:.6e}, RMS Integrated Error: {rms_err1:.6e}\n")
            f.write(f"  Particle 2 - Mean Integrated Error: {mean_err2:.6e}, RMS Integrated Error: {rms_err2:.6e}\n")
            f.write(f"  Average   - Mean Integrated Error: {(mean_err1+mean_err2)/2:.6e}, RMS Integrated Error: {(rms_err1+rms_err2)/2:.6e}\n")
            f.write("\n")
            
            print(f"Processed {ylabel}")
        
        # Plot trajectories in xy plane
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
        
        # Particle 1 trajectory
        ax1.plot(theory_results['particle1']['x'], theory_results['particle1']['y'], 
                'b-', linewidth=2, label='Theoretical')
        ax1.plot(particle1_data['x'], particle1_data['y'], 'r--', linewidth=2, label='Simulation')
        ax1.scatter(initial_conditions['r1_0'][0], initial_conditions['r1_0'][1], 
                   color='green', s=100, marker='o', label='Start', zorder=5)
        ax1.set_xlabel('X Position')
        ax1.set_ylabel('Y Position')
        ax1.set_title('Particle 1 Trajectory')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        ax1.set_aspect('equal')
        
        # Particle 2 trajectory
        ax2.plot(theory_results['particle2']['x'], theory_results['particle2']['y'], 
                'b-', linewidth=2, label='Theoretical')
        ax2.plot(particle2_data['x'], particle2_data['y'], 'r--', linewidth=2, label='Simulation')
        ax2.scatter(initial_conditions['r2_0'][0], initial_conditions['r2_0'][1], 
                   color='green', s=100, marker='o', label='Start', zorder=5)
        ax2.set_xlabel('X Position')
        ax2.set_ylabel('Y Position')
        ax2.set_title('Particle 2 Trajectory')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        ax2.set_aspect('equal')
        
        plt.tight_layout()
        plt.savefig('./Plots/trajectories_xy.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # Plot combined trajectory
        fig, ax = plt.subplots(figsize=(8, 8))
        
        # Theoretical
        ax.plot(theory_results['particle1']['x'], theory_results['particle1']['y'], 
               'b-', linewidth=2, label='Particle 1 (Theory)')
        ax.plot(theory_results['particle2']['x'], theory_results['particle2']['y'], 
               'g-', linewidth=2, label='Particle 2 (Theory)')
        
        # Simulation
        ax.plot(particle1_data['x'], particle1_data['y'], 
               'b--', linewidth=2, alpha=0.7, label='Particle 1 (Sim)')
        ax.plot(particle2_data['x'], particle2_data['y'], 
               'g--', linewidth=2, alpha=0.7, label='Particle 2 (Sim)')
        
        # Mark initial positions
        ax.scatter(initial_conditions['r1_0'][0], initial_conditions['r1_0'][1], 
                  color='blue', s=100, marker='o', edgecolor='black', linewidth=2, zorder=5)
        ax.scatter(initial_conditions['r2_0'][0], initial_conditions['r2_0'][1], 
                  color='green', s=100, marker='o', edgecolor='black', linewidth=2, zorder=5)
        
        ax.set_xlabel('X Position')
        ax.set_ylabel('Y Position')
        ax.set_title('Two-Body System Trajectories')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_aspect('equal')
        
        plt.tight_layout()
        plt.savefig('./Plots/trajectories_combined.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # Calculate trajectory errors
        traj_err1 = np.sqrt((particle1_data['x'] - theory_results['particle1']['x'])**2 + 
                           (particle1_data['y'] - theory_results['particle1']['y'])**2 +
                           (particle1_data['z'] - theory_results['particle1']['z'])**2)
        traj_err2 = np.sqrt((particle2_data['x'] - theory_results['particle2']['x'])**2 + 
                           (particle2_data['y'] - theory_results['particle2']['y'])**2 +
                           (particle2_data['z'] - theory_results['particle2']['z'])**2)
        
        mean_traj_err1 = np.mean(traj_err1)
        mean_traj_err2 = np.mean(traj_err2)
        max_traj_err1 = np.max(traj_err1)
        max_traj_err2 = np.max(traj_err2)
        
        f.write("\nTrajectory Error (Euclidean distance in 3D):\n")
        f.write(f"  Particle 1 - Mean Error: {mean_traj_err1:.6e}, Max Error: {max_traj_err1:.6e}\n")
        f.write(f"  Particle 2 - Mean Error: {mean_traj_err2:.6e}, Max Error: {max_traj_err2:.6e}\n")
        f.write(f"  Average   - Mean Error: {(mean_traj_err1+mean_traj_err2)/2:.6e}, Max Error: {(max_traj_err1+max_traj_err2)/2:.6e}\n")
        
        # Plot trajectory error over time
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.plot(times, traj_err1, 'b-', linewidth=2, label='Particle 1')
        ax.plot(times, traj_err2, 'g-', linewidth=2, label='Particle 2')
        ax.set_xlabel('Time')
        ax.set_ylabel('Position Error (Euclidean distance)')
        ax.set_title('Trajectory Error vs Time')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_yscale('log')
        
        plt.tight_layout()
        plt.savefig('./Plots/trajectory_error_vs_time.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # Energy conservation analysis
        f.write("\n" + "="*70 + "\n")
        f.write("Energy Conservation Analysis:\n\n")
        
        # Total energy from simulation
        total_energy_sim = particle1_data['ke'] + particle1_data['pe'] + \
                          particle2_data['ke'] + particle2_data['pe']
        
        # Total energy from theory
        total_energy_theory = theory_results['particle1']['ke'] + theory_results['particle1']['pe'] + \
                             theory_results['particle2']['ke'] + theory_results['particle2']['pe']
        
        # Energy drift
        energy_drift_sim = (total_energy_sim - total_energy_sim[0]) / np.abs(total_energy_sim[0])
        energy_drift_theory = (total_energy_theory - total_energy_theory[0]) / np.abs(total_energy_theory[0])
        
        f.write(f"Initial Total Energy:\n")
        f.write(f"  Simulation: {total_energy_sim[0]:.10e}\n")
        f.write(f"  Theory:     {total_energy_theory[0]:.10e}\n")
        f.write(f"\nFinal Total Energy:\n")
        f.write(f"  Simulation: {total_energy_sim[-1]:.10e}\n")
        f.write(f"  Theory:     {total_energy_theory[-1]:.10e}\n")
        f.write(f"\nRelative Energy Drift:\n")
        f.write(f"  Simulation: {energy_drift_sim[-1]:.6e}\n")
        f.write(f"  Theory:     {energy_drift_theory[-1]:.6e}\n")
        
        # Plot energy conservation
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 10))
        
        # Total energy
        ax1.plot(times, total_energy_sim, 'r-', linewidth=2, label='Simulation')
        ax1.plot(times, total_energy_theory, 'b--', linewidth=2, label='Theory')
        ax1.set_xlabel('Time')
        ax1.set_ylabel('Total Energy')
        ax1.set_title('Total Energy vs Time')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Energy drift
        ax2.plot(times, energy_drift_sim * 100, 'r-', linewidth=2, label='Simulation')
        ax2.plot(times, energy_drift_theory * 100, 'b--', linewidth=2, label='Theory')
        ax2.set_xlabel('Time')
        ax2.set_ylabel('Relative Energy Drift (%)')
        ax2.set_title('Energy Conservation')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('./Plots/energy_conservation.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    print("\n Analysis complete! Results saved in ./Plots/")
    print("- error_analysis.dat contains numerical error analysis")
    print("- PNG files show visual comparisons")

if __name__ == "__main__":
    main()
